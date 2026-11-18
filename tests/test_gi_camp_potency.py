import csv
import io
import json
import time
import zlib

import numpy as np
import pytest

from drugdevgym import ENVS, knowledge, make
from drugdevgym.analysis import (ec_fraction, estimate_pec50, estimate_pec50_anchored, standard_curve,
                                 titration_curve, titration_level, zprime)
from drugdevgym.chem import Kit, accumulated, hill, uninhibited_fraction
from drugdevgym.envs import GiCampPotency
from drugdevgym.envs import gi_camp_potency as gi
from drugdevgym.envs.dose_response import plate_cost
from drugdevgym.policies import (GI_CELL_RANGE, GI_DEFAULT_FORSKOLIN_UM, CAMP_EXPERT, GI_HANDBOOK_CELLS_K, POLICIES,
                                 GiCampPolicy, gi_conditions, random_gi_plan, run_episode)
from drugdevgym.rig import EDGE, ROW_NAMES
from drugdevgym.sim2real import worklist_csv

NAME = "gi_camp_potency"
HIDDEN_WORDS = ("truth", "pec50", "emax", "ec50", "fmol", "pde", "basal", "gi_max", "sticky", "kit_k")
TITRATION = {"cells_k": 10.0, "stim_min": 30.0, "ibmx": True, "replicates": 2}
PLATE = {**CAMP_EXPERT, "forskolin_uM": 10.0, "cells_k": 20.0}
ASSAY = gi.GiAssay(Kit(k=5.0), fsk_ec50=3.0, rate=2.0, basal_frac=0.02, pde=0.3, gi_max=0.9, cell_cv=0.05)


def play(policy, seed, rig="nominal"):
    return run_episode(make(NAME, seed=seed, rig=rig), POLICIES[NAME][policy](seed))


def drive(policy, seed, rig="nominal"):
    """One episode, keeping every observation, the final action and the final (reward, info)."""
    env, pol = make(NAME, seed=seed, rig=rig), POLICIES[NAME][policy](seed)
    obs, seen = env.reset(), []
    while True:
        seen.append(obs)
        action = pol.act(obs)
        obs, reward, done, info = env.step(action)
        if done:
            return env, seen + [obs], action, reward, info
        assert info == {}


# ---------------------------------------------------------------- registration

def test_the_environment_is_registered_with_its_policies_and_rules():
    assert ENVS[NAME] is GiCampPotency and make(NAME).name == NAME
    assert {"expert", "null", "random"} <= set(POLICIES[NAME])
    traps = [p for p in POLICIES[NAME] if p not in knowledge.NOT_KNOWLEDGE]
    assert 4 <= len(traps) <= 8
    rules = [r for r in knowledge.load_rules() if r["env"] == NAME]
    assert sorted(t for r in rules for t in r["traps"]) == sorted(traps)
    assert all(r["status"] == "draft" for r in rules)
    assert knowledge.orphan_traps(knowledge.load_rules()) == []


def test_spec_documents_actions_reward_and_is_json():
    spec = make(NAME).spec()
    assert spec["env"] == NAME and set(spec["actions"]) == {"run_titration", "run_plate", "report"}
    assert set(spec["actions"]["run_plate"]["design"]) == set(gi.DESIGN_SPEC)
    assert set(spec["actions"]["run_titration"]["design"]) == set(gi.TITRATION_SPEC)
    json.dumps(spec)


def test_spec_states_the_whole_reward_rule():
    text = make(NAME).spec()["reward"]
    for part in ("accuracy * gate - 0.15 * cost / 300",  # how the three terms combine
                 "max(0, 1 - |pEC50 error| / 0.5)", "tolerance is 0.5 log units", "left out of the report scores zero",
                 "Analysis: none is run for the agent",  # the reported numbers are scored as given
                 "Gate: 0.5, so accuracy is halved, unless at least one potency plate was run",
                 "and every potency plate run has at least 3 controls a side (pos and neg wells)",
                 "Z' >= 0.5 on the raw ratio; then it is 1",
                 "Z' = 1 - 3 * (sd of pos + sd of neg) / |mean of pos - mean of neg|", "sample standard deviations",
                 "A titration plate is not a potency plate and has no gate",
                 "10 per plate + 1 per well + 0.1 per tip",
                 "5 invalid actions end the episode with the cost term alone"):
        assert part in text, part
    assert gi.TOLERANCE == 0.5 and gi.ZPRIME_GATE == 0.5 and gi.MIN_CONTROLS == 3
    assert plate_cost(60, 76) == 10 + 60 + 7.6 and GiCampPotency.cost_reference == 300.0
    assert "between 0 and 14" in make(NAME).spec()["actions"]["report"]["pEC50"]


def test_five_invalid_actions_after_a_plate_leave_only_the_cost_term():
    env = make(NAME, seed=0)
    env.step({"type": "run_plate", "design": PLATE})
    for _ in range(4):
        env.step({"type": "dance"})
    _, reward, done, info = env.step({"type": "dance"})
    assert done and info == {"truncated": True} and env.cost > 0
    assert np.isclose(reward, -0.15 * env.cost / 300.0)


# ---------------------------------------------------------------- seeds, truth, observations

@pytest.mark.parametrize("policy", list(POLICIES[NAME]))
def test_same_seed_gives_the_same_episode(policy):
    a, b = drive(policy, 3), drive(policy, 3)
    assert json.dumps(a[1:]) == json.dumps(b[1:])  # observations, final action, reward and info
    assert a[0].step_log == b[0].step_log


def test_the_hidden_truth_is_salted_with_the_environment_name():
    assert GiCampPotency.salted
    for seed in range(5):
        first_draw = make(NAME, seed=seed).assay.kit.k
        salted = np.random.default_rng([seed, 1, zlib.crc32(NAME.encode())])
        unsalted = np.random.default_rng([seed, 1])  # the stream the four original environments draw from
        assert first_draw == float(salted.uniform(*gi.KIT_K_NM)) != float(unsalted.uniform(*gi.KIT_K_NM))


def test_editing_an_observation_cannot_reach_the_plates_the_gate_is_scored_on():
    reference = make(NAME, seed=0)
    reference.step({"type": "run_plate", "design": PLATE})
    honest = reference.step({"type": "report", "pEC50": {}})[3]
    assert honest["qc_pass"]
    env = make(NAME, seed=0)
    obs = env.step({"type": "run_plate", "design": PLATE})[0]
    assert obs["wells"] is not env.potency_plates[0] and obs["wells"][0] is not env.potency_plates[0][0]
    for w in obs["wells"]:  # wreck the controls in the agent's copy
        w["ratio"], w["role"] = 0.0, "cpd"
    obs["wells"].clear()
    info = env.step({"type": "report", "pEC50": {}})[3]
    assert info["qc_pass"] and info["zprime"] == honest["zprime"]
    env = make(NAME, seed=0)  # and the other way: a plate without controls cannot be edited into a pass
    obs = env.step({"type": "run_plate", "design": {**PLATE, "controls_per_side": 0}})[0]
    for i, w in enumerate(obs["wells"][:12]):
        w["role"], w["ratio"] = ("pos", 5000.0 + i) if i % 2 else ("neg", 1000.0 + i)
    assert not env.step({"type": "report", "pEC50": {}})[3]["qc_pass"]
    start = make(NAME, seed=0)
    first = start.reset()
    first["compounds"].clear()
    first["spec"]["labware"]["standard_series_nM"].clear()
    assert len(start.reset()["compounds"]) == 4 and len(start.spec()["labware"]["standard_series_nM"]) == 8
    assert len(gi.STD_SERIES_NM) == 8 and len(start.compounds) == 4


def test_the_truth_depends_on_the_seed_and_not_on_the_policy_or_the_rig():
    reference = play("expert", 9)[1]["truth"]
    assert set(reference) >= {"pEC50", "emax", "forskolin_ec50_uM", "camp_fmol_per_k_cells_min", "pde_per_min"}
    for policy in POLICIES[NAME]:
        assert play(policy, 9)[1]["truth"] == reference, policy
    assert play("expert", 9, rig="harsh")[1]["truth"] == reference
    assert play("expert", 10)[1]["truth"] != reference


@pytest.mark.parametrize("policy", list(POLICIES[NAME]))
@pytest.mark.parametrize("rig", ["nominal", "randomized", "harsh"])
def test_no_truth_field_names_in_any_observation(policy, rig):
    for seed in (5, 6):
        _env, seen, _action, reward, info = drive(policy, seed, rig)
        for obs in seen:
            text = json.dumps({k: v for k, v in obs.items() if k != "spec"}).lower()  # plain JSON, spec is public
            assert not any(word in text for word in HIDDEN_WORDS), (policy, rig, seed)
        assert "truth" in info and reward <= 1.0 and seen[-1] == {"phase": "done"}


def test_plates_carry_only_what_a_reader_gives():
    env = make(NAME, seed=0)
    env.reset()
    titration = env.step({"type": "run_titration", "design": TITRATION})[0]
    assert titration["kind"] == "titration" and titration["plates_left"] == 2
    assert {tuple(sorted(w)) for w in titration["wells"]} == {("forskolin_uM", "ratio", "role", "std_nM", "well")}
    roles = [w["role"] for w in titration["wells"]]
    assert roles.count("std") == 8 and roles.count("fsk") == 16 and roles.count("unstim") == 2
    assert sorted({w["forskolin_uM"] for w in titration["wells"] if w["role"] == "fsk"}) == sorted(gi.FSK_SERIES_UM)
    potency = env.step({"type": "run_plate", "design": PLATE})[0]
    assert potency["kind"] == "potency" and len(potency["wells"]) == 60
    assert {tuple(sorted(w)) for w in potency["wells"]} == {("conc_nM", "id", "ratio", "role", "std_nM", "well")}


# ---------------------------------------------------------------- malformed input

HUGE = 10 ** 400
BAD_NUMBERS = (float("nan"), float("inf"), -float("inf"), HUGE, -HUGE, [[1.0]], [1.0, 2.0], "10", None, {}, True)
BAD_BOOLS = ("yes", 1, 0.0, None, [True], float("nan"))


def rejected(action, seed=1):
    env = make(NAME, seed=seed)
    env.reset()
    obs, reward, done, info = env.step(action)
    return "error" in obs and reward == 0.0 and not done and info == {} and env.cost == 0.0 and env.plates == 0


@pytest.mark.parametrize("kind,good,spec", [("run_titration", TITRATION, gi.TITRATION_SPEC),
                                            ("run_plate", PLATE, gi.DESIGN_SPEC)])
def test_malformed_designs_are_rejected_without_exceptions(kind, good, spec):
    assert not rejected({"type": kind, "design": good})  # the baseline is valid
    for field, (ftype, lo, hi, _help) in spec.items():
        for bad in (BAD_BOOLS if ftype == "bool" else BAD_NUMBERS):
            assert rejected({"type": kind, "design": {**good, field: bad}}), (field, bad)
        if ftype != "bool":
            for bad in (lo - 1.0, hi + 1.0):
                assert rejected({"type": kind, "design": {**good, field: bad}}), (field, bad)
        if ftype == "int":
            assert rejected({"type": kind, "design": {**good, field: lo + 0.5}}), field
        assert rejected({"type": kind, "design": {k: v for k, v in good.items() if k != field}}), field
    assert rejected({"type": kind, "design": {**good, "lid": False}})  # unknown field
    for bad in (None, [], [good], "design", 3.0, float("nan")):
        assert rejected({"type": kind, "design": bad}), bad
    assert rejected({"type": kind})
    assert rejected({"type": kind, "design": good, "note": "x"})  # unknown field next to the design


def test_malformed_actions_are_rejected_without_exceptions():
    for bad in (None, [], "report", 3, float("nan"), {}, {"design": PLATE}, {"type": "dance"}, {"type": None},
                {"type": ["run_plate"], "design": PLATE}, {"type": {}}, {"type": float("nan")}, {"type": HUGE}):
        assert rejected(bad), bad


def test_malformed_reports_are_rejected_without_exceptions():
    for bad in (None, [], "7.5", 7.5, [["cpd1", 7.5]]):
        assert rejected({"type": "report", "pEC50": bad}), bad
    assert rejected({"type": "report"})
    assert rejected({"type": "report", "pEC50": {}, "plates_used": [1]})
    assert rejected({"type": "report", "pEC50": {"cpd9": 7.0}})
    assert rejected({"type": "report", "pEC50": {1: 7.0}})
    for bad in BAD_NUMBERS + (-1.0, 14.5):
        assert rejected({"type": "report", "pEC50": {"cpd1": bad}}), bad
    env = make(NAME, seed=1)
    _, reward, done, info = env.step({"type": "report", "pEC50": {"cpd1": np.float64(7.0), "cpd2": 7}})
    assert done and reward <= 0.5 and set(info["errors"]) == {"cpd1", "cpd2", "cpd3", "cpd4"}
    assert info["errors"]["cpd3"] is None  # a compound left out scores zero, it is not an error


def test_five_invalid_actions_end_the_episode_at_no_cost():
    env = make(NAME, seed=0)
    env.reset()
    for _ in range(4):
        obs, reward, done, info = env.step({"type": "run_plate", "design": {**PLATE, "n_points": 99}})
        assert "error" in obs and reward == 0.0 and not done
    obs, reward, done, info = env.step({"type": "dance"})
    assert done and reward == 0.0 and info == {"truncated": True}
    with pytest.raises(RuntimeError):
        env.step({"type": "report", "pEC50": {}})


def test_designs_that_do_not_fit_the_labware_are_rejected():
    assert rejected({"type": "run_plate", "design": {**PLATE, "n_points": 12, "replicates": 2}})  # 116 wells
    assert rejected({"type": "run_plate", "design": {**PLATE, "transfer_ul": 100.0, "dilution_factor": 10.0}})
    env = make(NAME, seed=0)
    assert "wells" in env.step({"type": "run_plate", "design": {**PLATE, "replicates": 4}})[0]["error"]


def test_three_plates_then_only_a_report():
    env = make(NAME, seed=0)
    for kind, design in (("run_titration", TITRATION), ("run_titration", TITRATION), ("run_plate", PLATE)):
        obs = env.step({"type": kind, "design": design})[0]
    assert obs["plates_left"] == 0
    for kind, design in (("run_titration", TITRATION), ("run_plate", PLATE)):
        assert "no plates left" in env.step({"type": kind, "design": design})[0]["error"]
    _, _, done, info = env.step({"type": "report", "pEC50": {}})
    assert done and len(info["zprime"]) == 1


def test_constructor_limits():
    with pytest.raises(ValueError):
        GiCampPotency(seed=0, n_compounds=9)
    with pytest.raises(ValueError):
        GiCampPotency(seed=-1)
    assert len(GiCampPotency(seed=0, n_compounds=2).reset()["compounds"]) == 2


# ---------------------------------------------------------------- reward and gate

def expected_reward(potency_plates, answers, info, cost):
    """The arithmetic of spec()["reward"], recomputed from what the agent saw and the final truth."""
    truth = info["truth"]["pEC50"]
    accuracy = np.mean([max(0.0, 1.0 - abs(answers[c] - truth[c]) / 0.5) if c in answers else 0.0 for c in truth])
    gate = len(potency_plates) >= 1
    for wells in potency_plates:
        pos = [w["ratio"] for w in wells if w["role"] == "pos"]
        neg = [w["ratio"] for w in wells if w["role"] == "neg"]
        gate = gate and len(pos) >= 3 and len(neg) >= 3 and zprime(pos, neg) >= 0.5
    return accuracy, gate, accuracy * (1.0 if gate else 0.5) - 0.15 * cost / 300.0


@pytest.mark.parametrize("policy", list(POLICIES[NAME]))
def test_reward_and_gate_arithmetic_match_the_spec_text(policy):
    text = make(NAME).spec()["reward"]
    assert "max(0, 1 - |pEC50 error| / 0.5)" in text and "halved, unless" in text
    assert "at least 3 controls a side" in text and "Z' >= 0.5" in text and "0.15 * cost / 300" in text
    for seed in range(6):
        env, seen, action, reward, info = drive(policy, seed)
        plates = [obs["wells"] for obs in seen if obs.get("kind") == "potency"]
        cost = sum(plate_cost(len(obs["wells"]), 0) for obs in seen if "wells" in obs)
        tips = sum(1 for s in env.step_log if s["op"] == "transfer")  # every transfer takes a fresh tip
        cost += 0.1 * tips
        accuracy, gate, expected = expected_reward(plates, action["pEC50"], info, cost)
        assert np.isclose(info["cost"], cost) and info["qc_pass"] == gate
        assert abs(info["accuracy"] - accuracy) < 2e-3  # the truth in info is rounded to 3 decimals
        assert abs(reward - expected) < 2e-3 and reward <= 1.0
        assert np.isclose(reward, info["accuracy"] * (1.0 if gate else 0.5) - 0.15 * info["cost"] / 300.0)


def test_the_gate_needs_a_potency_plate_with_enough_controls_that_pass():
    assert not play("null", 0)[1]["qc_pass"]  # no plate at all
    env = make(NAME, seed=0)
    env.step({"type": "run_titration", "design": TITRATION})
    assert not env.step({"type": "report", "pEC50": {}})[3]["qc_pass"]  # a titration is not a potency plate
    env = make(NAME, seed=0)
    env.step({"type": "run_plate", "design": {**PLATE, "controls_per_side": 2}})
    _, reward, _, info = env.step({"type": "report", "pEC50": {}})
    assert not info["qc_pass"] and np.isnan(info["zprime"][0]) and np.isclose(reward, -0.15 * info["cost"] / 300.0)
    env = make(NAME, seed=0)  # a later good plate does not excuse an earlier one without controls
    env.step({"type": "run_plate", "design": {**PLATE, "controls_per_side": 0}})
    env.step({"type": "run_plate", "design": PLATE})
    info = env.step({"type": "report", "pEC50": {}})[3]
    assert not info["qc_pass"] and len(info["zprime"]) == 2 and info["zprime"][1] >= 0.5
    env = make(NAME, seed=0)
    env.step({"type": "run_plate", "design": PLATE})
    assert env.step({"type": "report", "pEC50": {}})[3]["qc_pass"]


# ---------------------------------------------------------------- policies

def test_each_trap_is_the_expert_with_exactly_one_change():
    expert = POLICIES[NAME]["expert"](0)
    flags = ("use_standards", "anchor")
    assert expert.plan == CAMP_EXPERT and set(vars(expert)) == {"plan", *flags}
    assert expert.use_standards and expert.anchor
    # one named mistake that moves two plan fields: with no titration plate, both settings the
    # titration would have given (forskolin and cells) are fixed numbers
    no_titration = {"skip_titration": ["forskolin_uM", "cells_k"]}
    for name in POLICIES[NAME]:
        if name in knowledge.NOT_KNOWLEDGE:
            continue
        trap = POLICIES[NAME][name](0)
        assert type(trap) is GiCampPolicy and set(vars(trap)) == set(vars(expert)) and set(trap.plan) == set(CAMP_EXPERT)
        changed = [k for k in CAMP_EXPERT if trap.plan[k] != CAMP_EXPERT[k]]
        changed += [f for f in flags if getattr(trap, f) != getattr(expert, f)]
        if name in no_titration:
            assert changed == no_titration[name], (name, changed)
            assert all(isinstance(trap.plan[k], float) for k in changed)  # numbers, not rules read off a titration
        else:
            assert len(changed) == 1, (name, changed)
    assert set(no_titration) <= set(POLICIES[NAME])


def test_skip_titration_runs_one_potency_plate_on_handbook_conditions_and_nothing_else():
    env, seen, _action, _reward, _info = drive("skip_titration", 0)
    assert [obs.get("kind") for obs in seen] == [None, "potency", None] and {s["run"] for s in env.step_log} == {1}
    labels = {s["what"] for s in env.step_log if s["op"] == "fill"}
    assert f"forskolin {GI_DEFAULT_FORSKOLIN_UM * 5:g} uM" in labels
    assert f"cells {GI_HANDBOOK_CELLS_K:g}k per well with IBMX" in labels
    plan = POLICIES[NAME]["skip_titration"](0).plan
    assert plan == {**CAMP_EXPERT, "forskolin_uM": GI_DEFAULT_FORSKOLIN_UM, "cells_k": GI_HANDBOOK_CELLS_K}
    for name in POLICIES[NAME]:  # every other trap keeps the titration plate
        if name not in knowledge.NOT_KNOWLEDGE and name != "skip_titration":
            assert [obs.get("kind") for obs in drive(name, 0)[1]] == [None, "titration", "potency", None], name


def test_every_gi_rule_names_exactly_one_trap_and_every_trap_one_rule():
    rules = [r for r in knowledge.load_rules() if r["env"] == NAME]
    assert all(len(r["traps"]) == 1 for r in rules)
    named = [r["traps"][0] for r in rules]
    assert len(named) == len(set(named))
    by_trap = {r["traps"][0]: r for r in rules}
    assert "titrat" in by_trap["skip_titration"]["said"].lower()
    assert "titrat" not in by_trap["forskolin_too_high"]["said"].lower()
    assert "standard curve" not in by_trap["forskolin_too_high"]["why"]  # that clause was not the mechanism


def test_expert_beats_null_on_held_out_seeds_by_a_clear_margin():
    seeds = range(30000, 30060)
    expert = np.array([play("expert", s)[0] for s in seeds])
    null = np.array([play("null", s)[0] for s in seeds])
    assert expert.mean() > null.mean() + 0.3
    assert np.mean(expert > null) > 0.9


def test_null_runs_nothing_and_random_plans_always_fit_the_plate():
    reward, info = play("null", 4)
    assert info["cost"] == 0.0 and not info["qc_pass"] and reward <= 0.5
    for seed in range(200):
        plan = random_gi_plan(np.random.default_rng([seed, 99]))
        fixed = {**plan, "forskolin_uM": 10.0 if plan["forskolin_uM"] == "EC80" else plan["forskolin_uM"],
                 "cells_k": 10.0 if plan["cells_k"] == "titrate" else plan["cells_k"]}
        assert gi.check_design(fixed, 4)[1] is None, plan
    for seed in range(20):
        _env, seen, _action, _reward, info = drive("random", seed)
        assert not any("error" in obs for obs in seen) and "truth" in info


def test_the_expert_sets_forskolin_and_cells_from_its_titration():
    env, seen, _action, _reward, info = drive("expert", 0)
    assert [obs.get("kind") for obs in seen] == [None, "titration", "potency", None]
    forskolin, cells = gi_conditions(CAMP_EXPERT, seen[1]["wells"])
    labels = {s["what"] for s in env.step_log if s["op"] == "fill"}
    assert f"forskolin {forskolin * 5:g} uM" in labels and f"cells {cells:g}k per well with IBMX" in labels
    assert 0.0 < forskolin <= gi.FSK_MAX_UM and GI_CELL_RANGE[0] <= cells <= GI_CELL_RANGE[1]
    true_ec80 = 4.0 * info["truth"]["forskolin_ec50_uM"]
    assert true_ec80 / 4.0 < forskolin < true_ec80 * 4.0
    assert gi_conditions({**CAMP_EXPERT, "forskolin_uM": 100.0, "cells_k": 80.0}, seen[1]["wells"]) == (100.0, 80.0)
    no_standards = [w for w in seen[1]["wells"] if w["role"] != "std"]
    assert gi_conditions(CAMP_EXPERT, no_standards) == (10.0, 10.0)  # nothing to read: defaults


# ---------------------------------------------------------------- chemistry

def test_accumulation_is_linear_without_loss_and_plateaus_with_it():
    assert accumulated(0.0, 30.0) == 30.0
    assert np.isclose(accumulated(1e-9, 30.0), 30.0)
    assert np.isclose(accumulated(0.3, 1e4), 1.0 / 0.3)
    assert np.isclose(accumulated(0.3, 5.0), (1 - np.exp(-1.5)) / 0.3)
    assert accumulated(0.3, 10.0) < accumulated(0.3, 20.0) < accumulated(0.01, 20.0) < 20.0


def test_competitive_inhibitor_leaves_one_over_one_plus_i_over_ki():
    assert uninhibited_fraction(0.0, 20.0) == 1.0 and uninhibited_fraction(-5.0, 20.0) == 1.0
    assert np.isclose(uninhibited_fraction(20.0, 20.0), 0.5)
    assert np.isclose(uninhibited_fraction(gi.IBMX_UM, gi.IBMX_KI_UM), 1.0 / 26.0)


def test_forskolin_raises_camp_and_the_agonist_removes_a_fixed_fraction_of_that_rise():
    def made(forskolin, inhibition, ibmx=gi.IBMX_UM, minutes=30.0, cells=10.0):
        return gi.camp_fmol(ASSAY, cells, forskolin, inhibition, ibmx, minutes)

    unstimulated = made(0.0, 0.0)
    assert unstimulated < made(1.0, 0.0) < made(10.0, 0.0) < made(100.0, 0.0)
    assert made(0.0, 0.9) == unstimulated  # nothing for the Gi pathway to pull down
    for forskolin in (1.0, 3.0, 100.0):  # noncompetitive: the same fraction at any forskolin
        rise = made(forskolin, 0.0) - unstimulated
        assert np.isclose(made(forskolin, 0.6) - unstimulated, 0.4 * rise)
    assert np.isclose(made(3.0, 0.0, cells=20.0), 2.0 * made(3.0, 0.0))
    assert np.isclose(gi.cyclase_drive(ASSAY, 3.0, 0.0), 0.02 + 0.98 * 0.5)
    assert np.isclose(float(hill(3.0, ASSAY.fsk_ec50)), 0.5)


def test_without_ibmx_the_signal_collapses_to_a_steady_state():
    def made(ibmx, minutes):
        return gi.camp_fmol(ASSAY, 10.0, 10.0, 0.0, ibmx, minutes)

    assert made(0.0, 30.0) < 0.2 * made(gi.IBMX_UM, 30.0)
    assert np.isclose(made(0.0, 60.0), made(0.0, 120.0), rtol=1e-6)  # flat: synthesis equals breakdown
    assert made(gi.IBMX_UM, 120.0) > 2.5 * made(gi.IBMX_UM, 30.0)  # with IBMX it keeps piling up


def test_hidden_truth_is_drawn_inside_the_documented_ranges():
    for seed in range(50):
        env = make(NAME, seed=seed)
        a = env.assay
        for value, (lo, hi) in ((a.kit.k, gi.KIT_K_NM), (a.fsk_ec50, gi.FSK_EC50_UM), (a.rate, gi.CAMP_RATE),
                                (a.basal_frac, gi.BASAL_FRAC), (a.pde, gi.PDE_PER_MIN), (a.gi_max, gi.GI_MAX),
                                (a.cell_cv, gi.CELL_CV)):
            assert lo <= value <= hi
        for c in env.compounds:
            assert gi.PEC50[0] <= c.pec50 <= gi.PEC50[1] and gi.HILL[0] <= c.hill <= gi.HILL[1]
            assert a.gi_max * gi.EFFICACY[0] <= c.emax <= a.gi_max * gi.EFFICACY[1] and c.sticky > 0


# ---------------------------------------------------------------- the plates, through the rig

def ratios(wells, role, **match):
    return [w["ratio"] for w in wells if w["role"] == role and all(w[k] == v for k, v in match.items())]


def test_the_ratio_falls_with_forskolin_and_rises_as_the_agonist_acts():
    env = make(NAME, seed=0)
    titration = env.step({"type": "run_titration", "design": TITRATION})[0]["wells"]
    assert np.mean(ratios(titration, "fsk", forskolin_uM=100.0)) < np.mean(ratios(titration, "unstim")) - 500
    potency = env.step({"type": "run_plate", "design": PLATE})[0]["wells"]
    assert np.mean(ratios(potency, "pos")) > np.mean(ratios(potency, "neg")) + 500
    for cid in ("cpd1", "cpd2", "cpd3", "cpd4"):
        curve = sorted((w["conc_nM"], w["ratio"]) for w in potency if w["id"] == cid)
        assert np.mean([r for _, r in curve[-2:]]) > np.mean([r for _, r in curve[:2]])


def test_no_forskolin_means_no_window_and_no_ibmx_means_less_camp():
    def controls(**changes):
        env = make(NAME, seed=0)
        wells = env.step({"type": "run_plate", "design": {**PLATE, **changes}})[0]["wells"]
        neg, pos = ratios(wells, "neg"), ratios(wells, "pos")
        return np.mean(neg), np.mean(pos), zprime(pos, neg)

    neg, pos, z = controls()
    assert pos - neg > 500 and z >= 0.5
    neg0, pos0, z0 = controls(forskolin_uM=0.0)
    assert abs(pos0 - neg0) < 300 and z0 < 0.0 and neg0 > neg + 1000  # both sit at the unstimulated level
    neg_open, _pos, _z = controls(ibmx=False)
    assert neg_open > neg + 500  # less cAMP in the forskolin control: a higher ratio


def test_a_plate_fills_the_inner_wells_first_and_spills_onto_the_outer_ring():
    def edge(well):
        return bool(EDGE[ROW_NAMES.index(well[0]), int(well[1:]) - 1])

    env = make(NAME, seed=0)
    small = env.step({"type": "run_plate", "design": PLATE})[0]["wells"]
    assert len(small) == 60 and not any(edge(w["well"]) for w in small)
    big = env.step({"type": "run_plate", "design": {**PLATE, "n_points": 12}})[0]["wells"]
    assert len(big) == 68 and not any(edge(w["well"]) for w in big[:60]) and all(edge(w["well"]) for w in big[60:])
    assert len({w["well"] for w in big}) == 68


def test_the_protocol_runs_on_the_rig_and_the_worklist_names_every_reagent():
    env = make(NAME, seed=0)
    run_episode(env, POLICIES[NAME]["expert"](0))
    rows = list(csv.DictReader(io.StringIO(worklist_csv(env.step_log))))
    assert len(rows) == len(env.step_log) > 300 and {r["run"] for r in rows} == {"1", "2"}
    assert {r["op"] for r in rows} == {"fill", "transfer", "incubate", "read"}
    labels = {s["what"] for s in env.step_log if s["op"] == "fill"}
    assert {"diluent", "vehicle", "reference agonist", "detection mix", "stimulation buffer"} <= labels
    for start in ("cells ", "forskolin ", "cAMP standard ", "cpd1 at "):
        assert any(x.startswith(start) for x in labels), start
    assert all(s.get("what") for s in env.step_log if s["op"] == "fill")
    assert all(s["lid"] for s in env.step_log if s["op"] == "incubate")
    assert all(s["new_tip"] for s in env.step_log if s["op"] == "transfer")
    ladder = [s for s in env.step_log if s["op"] == "transfer" and s["dst_plate"] == "dilution" and s["run"] == 2]
    assert len(ladder) == 4 * 9 and all(s["mix_after"] and s["ul"] == CAMP_EXPERT["transfer_ul"] for s in ladder)


def test_dead_volume_and_mixing_act_through_the_rig():
    def mean_accuracy(**changes):
        out = []
        for seed in range(8):
            env = make(NAME, seed=seed)
            wells = env.step({"type": "run_plate", "design": {**PLATE, **changes}})[0]["wells"]
            out.append(env.step({"type": "report", "pEC50": estimate_pec50_anchored(wells)})[3]["accuracy"])
        return np.mean(out)

    sound = mean_accuracy()
    assert mean_accuracy(transfer_ul=5.0) < 0.1 < sound  # a 20 uL ladder well is all dead volume after the transfer
    assert mean_accuracy(mix=False) < sound - 0.2


def test_one_episode_is_fast():
    start = time.perf_counter()
    for seed in range(10):
        play("expert", seed)
    assert (time.perf_counter() - start) / 10 < 0.25  # measured near 0.02 s; the bound leaves room for a slow runner


# ---------------------------------------------------------------- analysis

def ideal_titration(ec50=3.0, top=20.0, bottom=0.5, kit=Kit(k=5.0)):
    wells = [{"well": "", "role": "std", "forskolin_uM": None, "std_nM": s, "ratio": float(kit.ratio(s))}
             for s in gi.STD_SERIES_NM]
    for f in gi.FSK_SERIES_UM:
        level = bottom + (top - bottom) * float(hill(f, ec50))
        wells.append({"well": "", "role": "fsk", "forskolin_uM": f, "std_nM": None, "ratio": float(kit.ratio(level))})
    wells.append({"well": "", "role": "unstim", "forskolin_uM": 0.0, "std_nM": None, "ratio": float(kit.ratio(bottom))})
    return wells


def test_ec_fraction_is_the_hill_inverse():
    assert np.isclose(ec_fraction(10.0, 1.0, 0.8), 40.0) and np.isclose(ec_fraction(10.0, 2.0, 0.8), 20.0)
    assert np.isclose(ec_fraction(10.0, 1.3, 0.5), 10.0)
    assert np.isclose(float(hill(ec_fraction(7.0, 1.4, 0.8), 7.0, 1.4)), 0.8)
    for ec50, slope, fraction in ((1.0, 1.0, 1.0), (1.0, 1.0, 0.0), (1.0, 1.0, 1.2), (1.0, 1.0, -0.1),
                                  (1.0, 0.0, 0.8), (1.0, -1.0, 0.8), (1.0, float("nan"), 0.8),
                                  (1.0, 1.0, float("nan"))):
        with pytest.raises(ValueError):
            ec_fraction(ec50, slope, fraction)


def test_standard_curve_needs_four_standards():
    wells = ideal_titration()
    fit = standard_curve(wells)
    assert abs(10.0 ** fit["logec50"] - 5.0) < 0.5
    assert standard_curve([w for w in wells if w["role"] != "std"][:3]) is None
    assert standard_curve([w for w in wells if w["role"] == "std"][:3]) is None


def test_titration_curve_recovers_the_forskolin_response_in_analyte_units():
    curve = titration_curve(ideal_titration())
    assert curve["fitted"] and abs(np.log10(curve["ec50_uM"] / 3.0)) < 0.1 and abs(curve["hill"] - 1.0) <= 0.2
    assert abs(curve["kit_mid_nM"] - 5.0) < 0.5 and abs(curve["top_nM"] - 20.0) < 2.5
    assert abs(curve["unstim_nM"] - 0.5) < 0.1 and abs(curve["top_wells_nM"] - 19.4) < 2.0
    assert abs(titration_level(curve, 3.0) - 10.25) < 1.5 and titration_level(curve, 0.0) == curve["unstim_nM"]
    assert titration_level(curve, 1.0) < titration_level(curve, 10.0) < titration_level(curve, 100.0)


def test_titration_curve_says_when_there_is_nothing_to_fit():
    wells = ideal_titration()
    assert titration_curve([w for w in wells if w["role"] != "std"]) is None
    flat = [{**w, "ratio": 5000.0} if w["role"] == "fsk" else w for w in wells]
    curve = titration_curve(flat)
    assert not curve["fitted"] and curve["ec50_uM"] is None and curve["top_wells_nM"] > 0
    assert titration_level(curve, 10.0) == curve["top_wells_nM"]  # the only estimate there is
    only_standards = titration_curve([w for w in wells if w["role"] == "std"])
    assert only_standards["top_wells_nM"] is None and titration_level(only_standards, 0.0) is None
    falling = [{**w, "ratio": 6600.0 - w["ratio"]} if w["role"] == "fsk" else w for w in wells]
    assert not titration_curve(falling)["fitted"]


def inhibition_plate(pec50=9.2, emax=0.6, noise=0.0, seed=0, controls=6, kit=Kit(k=5.0), level=15.0):
    """A potency plate whose ladder stops short of the uninhibited plateau of a potent agonist."""
    rng = np.random.default_rng(seed)
    wells = [{"role": "std", "id": None, "conc_nM": None, "std_nM": s, "ratio": float(kit.ratio(s))}
             for s in gi.STD_SERIES_NM]

    def read(camp):
        return float(kit.ratio(camp * (1 + rng.normal(0.0, noise))))

    for _ in range(controls):
        wells.append({"role": "neg", "id": None, "conc_nM": None, "std_nM": None, "ratio": read(level)})
        wells.append({"role": "pos", "id": None, "conc_nM": None, "std_nM": None, "ratio": read(level * 0.1)})
    for j in range(10):
        conc = 30000.0 / 4.0 ** j
        camp = level * (1 - emax * float(hill(conc, 10.0 ** (9.0 - pec50))))
        wells.append({"role": "cpd", "id": "c", "conc_nM": conc, "std_nM": None, "ratio": read(camp)})
    return wells


def test_anchoring_on_the_vehicle_controls_steadies_a_potent_inhibition_curve():
    clean = estimate_pec50_anchored(inhibition_plate())["c"]
    assert abs(clean - 9.2) < 0.1
    free, anchored = [], []
    for seed in range(40):
        wells = inhibition_plate(noise=0.08, seed=seed)
        free.append(abs(estimate_pec50(wells).get("c", 12.0) - 9.2))
        anchored.append(abs(estimate_pec50_anchored(wells).get("c", 12.0) - 9.2))
    assert np.mean(anchored) < np.mean(free)


def test_anchored_fit_without_vehicle_controls_is_the_plain_fit():
    wells = inhibition_plate(pec50=7.0, noise=0.05, controls=0)
    assert estimate_pec50_anchored(wells) == estimate_pec50(wells)
    assert estimate_pec50_anchored(wells, use_standards=False) == estimate_pec50(wells, use_standards=False)
    flat = [{**w, "ratio": 3000.0} if w["role"] == "cpd" else w for w in inhibition_plate()]
    assert estimate_pec50_anchored(flat) == {}  # the plateau is pinned and nothing moves: no curve


def test_raw_ratio_puts_the_midpoint_in_the_wrong_place():
    wells = inhibition_plate(pec50=7.0, emax=0.9)
    through_curve = estimate_pec50_anchored(wells, use_standards=True)["c"]
    raw = estimate_pec50_anchored(wells, use_standards=False)["c"]
    assert abs(through_curve - 7.0) < 0.1 and abs(raw - 7.0) > 0.2
