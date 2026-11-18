import csv
import io
import json
import time

import numpy as np
import pytest

from drugdevgym import envqc, knowledge, make
from drugdevgym.analysis import enzyme_activity, estimate_fluid_half_life, response
from drugdevgym.chem import acid_catalysis_factor, enzyme_exposure, q10_factor, remaining
from drugdevgym.envs import GiFluidStability
from drugdevgym.envs import gi_fluid_stability as gi
from drugdevgym.policies import (GI_EXPERT, GI_INTERACTIONS, GI_WARM_OVERNIGHT, POLICIES, GiFluidPolicy,
                                 run_episode)
from drugdevgym.rig import RigProfile
from drugdevgym.sim2real import worklist_csv

NAME = "gi_fluid_stability"
TRAPS = [p for p in POLICIES[NAME] if p not in knowledge.NOT_KNOWLEDGE]
SINGLE_TRAPS = [p for p in TRAPS if p not in GI_INTERACTIONS]  # the expert with exactly one change each
# names of everything hidden: the truth block and the attributes it is built from
HIDDEN_WORDS = ("truth", "enzymatic", "acid_min", "acid_half_life", "lot_activity")
BOTH = {"cpd1": ["sgf", "sif"], "cpd2": ["sgf", "sif"], "cpd3": ["sgf", "sif"]}
DESIGN = {**GI_EXPERT, "fluids": BOTH}


def play(policy, seed, rig="nominal"):
    return run_episode(make(NAME, seed=seed, rig=rig), POLICIES[NAME][policy](seed))


def walk(policy_name, seed, rig="nominal"):
    """One episode step by step: (observations, actions, final reward, final info, env)."""
    env, policy = make(NAME, seed=seed, rig=rig), POLICIES[NAME][policy_name](seed)
    obs, done, seen, actions = env.reset(), False, [], []
    while not done:
        seen.append(obs)
        actions.append(policy.act(obs))
        obs, reward, done, info = env.step(actions[-1])
    seen.append(obs)
    return seen, actions, reward, info, env


def simulate(peptides, seed=0, lot=None, **changes):
    design, err = gi.check_design({**GI_EXPERT, "fluids": {p.cid: list(gi.FLUIDS) for p in peptides}, **changes},
                                  [p.cid for p in peptides])
    assert err is None, err
    return gi.simulate_gi(RigProfile(), np.random.default_rng(seed), peptides, lot or {"sgf": 1.0, "sif": 1.0}, design)


def t0_response(samples, cid, fluid):
    return float(np.mean([response(s) for s in samples if (s["id"], s["fluid"], s["time_min"]) == (cid, fluid, 0.0)]))


# ---------------------------------------------------------------- the contract

@pytest.mark.parametrize("policy", ["expert", "random", "acid_quench_in_sgf"])
def test_same_seed_gives_the_same_episode(policy):
    first, second = walk(policy, 3), walk(policy, 3)
    assert first[0] == second[0] and first[1] == second[1]  # every observation and every action
    assert first[2] == second[2] and first[3] == second[3]  # reward and final info
    assert walk(policy, 3)[3]["truth"] != walk(policy, 4)[3]["truth"]


def test_the_truth_depends_on_the_seed_and_not_on_the_policy():
    truths = [play(policy, 9)[1]["truth"] for policy in POLICIES[NAME]]
    assert all(t == truths[0] for t in truths)
    assert set(truths[0]) == {"half_life_min", "acid_half_life_min", "lot_activity"}
    assert set(truths[0]["half_life_min"]["cpd1"]) == set(gi.FLUIDS)
    assert play("expert", 10)[1]["truth"] != truths[0]
    for rig in ("randomized", "harsh"):
        assert play("expert", 9, rig)[1]["truth"] == truths[0]  # the rig is not part of the truth


@pytest.mark.parametrize("policy", list(POLICIES[NAME]))
@pytest.mark.parametrize("rig", ["nominal", "randomized", "harsh"])
def test_no_observation_carries_a_truth_field(policy, rig):
    for seed in (5, 6):
        seen, _, reward, info, _ = walk(policy, seed, rig)
        for obs in seen:
            text = json.dumps(obs).lower()  # plain JSON, the task text included
            assert not any(word in text for word in HIDDEN_WORDS), (policy, rig)
        assert "truth" in info and reward <= 1.0
        assert all(word in json.dumps(info) for word in ("truth", "acid_half_life", "lot_activity"))


def test_editing_an_observation_cannot_flip_the_qc_gate_or_touch_the_data():
    # the cheap run: one replicate and no activity control fails the gate
    cheap = {**DESIGN, "replicates": 1, "activity_control": False, "timepoints_min": [0, 5, 30, 240]}
    honest = make(NAME, seed=0)
    honest.step({"type": "run", "design": cheap})
    assert not honest.step({"type": "report", "half_life_min": {}})[3]["qc_pass"]
    # the same run, with the observation doctored into one that would pass
    env = make(NAME, seed=0)
    obs = env.step({"type": "run", "design": cheap})[0]
    kept = json.dumps(env.samples)
    assert obs["samples"] == env.samples and obs["samples"] is not env.samples
    assert all(a is not b for a, b in zip(obs["samples"], env.samples))
    doctored = obs["samples"]
    doctored += [{**s, "replicate": 2} for s in doctored if s["time_min"] == 0.0]
    for fluid in gi.FLUIDS:
        for t in (0.0, 15.0):
            doctored += [{"id": "control", "fluid": fluid, "replicate": r, "time_min": t,
                          "area": 20000.0 * 0.5 ** (t / 15.0), "is_area": 100.0} for r in (1, 2)]
    assert gi.run_qc(doctored)["qc_pass"]  # the forged table would pass
    info = env.step({"type": "report", "half_life_min": {}})[3]
    assert not info["qc_pass"] and info["enzyme_activity"] == {} and json.dumps(env.samples) == kept
    wiped = make(NAME, seed=0)
    wiped.step({"type": "run", "design": DESIGN})[0]["samples"].clear()
    untouched = make(NAME, seed=0)
    untouched.step({"type": "run", "design": DESIGN})
    assert wiped.step({"type": "report", "half_life_min": {}})[3] == untouched.step({"type": "report", "half_life_min": {}})[3]
    # reset and spec hand out fresh objects too
    env = make(NAME, seed=0)
    first = env.reset()
    first["compounds"].clear(), first["fluids"].clear(), first["spec"]["actions"].clear()
    first["spec"]["fluids"]["sif"]["enzyme_activity_half_life_37C_min"] = 1.0
    again = env.reset()
    assert again["compounds"] == ["cpd1", "cpd2", "cpd3"] and again["fluids"] == ["sgf", "sif"]
    assert again["spec"] == make(NAME, seed=0).spec() and gi.ACTIVITY_HL_37C_MIN["sif"] == 240.0


def test_the_hidden_truth_is_salted_so_another_environment_does_not_give_it_away():
    assert GiFluidStability.salted
    lo, hi = np.log10(gi.ENZYMATIC_HL_RANGE_MIN["sgf"])
    for seed in range(40000, 40020):
        # stability reveals its first hidden draw at the end of a free episode; unsalted, that draw is this one too
        draw = (np.log10(make("stability", seed=seed).half_lives["cpd1"]) - np.log10(5.0)) / np.log10(2000.0 / 5.0)
        would_leak = 10 ** (lo + draw * (hi - lo))
        assert not np.isclose(make(NAME, seed=seed).peptides[0].sgf_enzymatic_min, would_leak, rtol=1e-3), seed


BAD_DESIGNS = [
    None, [], "design", 3, {},
    {**DESIGN, "extra": 1}, {**DESIGN, 1: 2},
    {k: v for k, v in DESIGN.items() if k != "replicates"}, {k: v for k, v in DESIGN.items() if k != "sgf_quench"},
    {k: v for k, v in DESIGN.items() if k != "fluids"},
    {**DESIGN, "replicates": float("nan")}, {**DESIGN, "replicates": float("inf")}, {**DESIGN, "replicates": 10 ** 400},
    {**DESIGN, "replicates": 2.5}, {**DESIGN, "replicates": "2"}, {**DESIGN, "replicates": True},
    {**DESIGN, "replicates": [2]}, {**DESIGN, "replicates": None}, {**DESIGN, "replicates": 0},
    {**DESIGN, "timepoints_min": [[0], [5]]}, {**DESIGN, "timepoints_min": [0, float("nan")]},
    {**DESIGN, "timepoints_min": [0, 10 ** 400]}, {**DESIGN, "timepoints_min": [0]},
    {**DESIGN, "timepoints_min": [0, 0.5]},
    {**DESIGN, "timepoints_min": "0,5"}, {**DESIGN, "timepoints_min": [0, "5"]},
    {**DESIGN, "timepoints_min": list(range(9))},
    {**DESIGN, "timepoints_min": [0, -1]}, {**DESIGN, "timepoints_min": [0, 500]}, {**DESIGN, "timepoints_min": {}},
    {**DESIGN, "timepoints_min": [0, True]},
    {**DESIGN, "internal_standard": 1}, {**DESIGN, "fresh_enzyme": "yes"}, {**DESIGN, "activity_control": None},
    {**DESIGN, "prequenched_t0": [True]},
    {**DESIGN, "hold_min": -1}, {**DESIGN, "hold_min": float("nan")}, {**DESIGN, "hold_min": float("inf")},
    {**DESIGN, "hold_min": "30"}, {**DESIGN, "hold_min": 10 ** 400}, {**DESIGN, "hold_min": [30]},
    {**DESIGN, "hold_temp_c": 0}, {**DESIGN, "hold_temp_c": 99}, {**DESIGN, "hold_temp_c": float("nan")},
    {**DESIGN, "sgf_quench": "bleach"}, {**DESIGN, "sgf_quench": 3}, {**DESIGN, "sgf_quench": None},
    {**DESIGN, "sif_quench": ["acid"]}, {**DESIGN, "sif_quench": {"acid": 1}}, {**DESIGN, "sif_quench": float("nan")},
    {**DESIGN, "fluids": None}, {**DESIGN, "fluids": []}, {**DESIGN, "fluids": "both"}, {**DESIGN, "fluids": {}},
    {**DESIGN, "fluids": {"cpd9": ["sgf"]}}, {**DESIGN, "fluids": {"cpd1": "sgf"}}, {**DESIGN, "fluids": {"cpd1": ["bile"]}},
    {**DESIGN, "fluids": {"cpd1": ["sgf", "sgf"]}}, {**DESIGN, "fluids": {"cpd1": [["sgf"]]}},
    {**DESIGN, "fluids": {"cpd1": [1]}}, {**DESIGN, "fluids": {"cpd1": [float("nan")]}}, {**DESIGN, "fluids": {"cpd1": []}},
    {**DESIGN, "fluids": {1: ["sgf"]}}, {**DESIGN, "fluids": {"cpd1": {"sgf": True}}}, {**DESIGN, "fluids": 10 ** 400},
]


@pytest.mark.parametrize("bad", BAD_DESIGNS, ids=range(len(BAD_DESIGNS)))
def test_malformed_designs_are_rejected_without_exceptions(bad):
    env = make(NAME, seed=1)
    obs, reward, done, info = env.step({"type": "run", "design": bad})
    assert "error" in obs and reward == 0.0 and not done and info == {}
    assert env.cost == 0.0 and env.samples is None and env.step_log == []
    json.dumps(obs)


def test_malformed_actions_are_rejected_and_five_end_the_episode():
    env = make(NAME, seed=1)
    for bad in (None, [], "run", {}):
        obs, reward, done, info = env.step(bad)
        assert "error" in obs and reward == 0.0 and not done and info == {}
    obs, reward, done, info = env.step({"type": ["run"]})
    assert done and reward == 0.0 and info == {"truncated": True} and "truth" not in json.dumps(obs)
    with pytest.raises(RuntimeError):
        env.step({"type": "report"})
    for bad in ({"type": None}, {"type": "dance"}, {"type": "run"}, {"type": 10 ** 400}, {"type": 10 ** 5000},
                {"type": float("nan")}, {"type": {"run": 1}}):
        obs, reward, done, info = make(NAME, seed=1).step(bad)
        assert "error" in obs and reward == 0.0 and not done and info == {}
        json.dumps(obs)
    env = make(NAME, seed=1)
    assert "error" not in env.step({"type": "run", "design": DESIGN})[0]
    assert "has been run" in env.step({"type": "run", "design": DESIGN})[0]["error"]  # one time course per episode


@pytest.mark.parametrize("answers", [None, [], "fast", 3, {"cpd1": 5.0}, {"cpd1": [5.0, 5.0]},
                                     {"cpd1": {"sgf": float("nan")}}, {"cpd1": {"sgf": float("inf")}},
                                     {"cpd1": {"sgf": 10 ** 400}}, {"cpd1": {"sgf": 10 ** 5000}}, {"cpd1": {"sgf": -5}},
                                     {"cpd1": {"sgf": "12"}}, {"cpd1": {"sgf": True}}, {"cpd1": {"sgf": [12]}},
                                     {"cpd1": {"sgf": {}}}, {"cpd1": {"bile": 3.0}}, {"cpd9": {"sgf": 3.0}}])
def test_malformed_answers_score_zero_without_exceptions(answers):
    obs, reward, done, info = make(NAME, seed=1).step({"type": "report", "half_life_min": answers})
    assert done and obs == {"phase": "done"} and info["accuracy"] == 0.0 and reward == 0.0 and not info["qc_pass"]


def spec_gate(samples):
    """The QC gate exactly as spec()["reward"] words it, with the numbers written out."""
    if not samples:
        return False
    activity = enzyme_activity(samples, "control", 15.0, {"sgf": 1440.0, "sif": 240.0}, use_lloq=True, lloq_frac=0.02)
    for fluid in {s["fluid"] for s in samples}:
        rows = [s for s in samples if s["fluid"] == fluid]
        for cid in {s["id"] for s in rows}:
            y0 = [response(s) for s in rows if s["id"] == cid and s["time_min"] == 0]
            if len(y0) < 2 or np.mean(y0) <= 0 or np.std(y0, ddof=1) / np.mean(y0) > 0.15:
                return False
        if fluid not in activity or not 0.8 <= activity[fluid] <= 1.25:
            return False
    return True


def spec_reward(env, samples, answers, cost):
    """accuracy, gate and reward recomputed from the words of spec()["reward"] alone."""
    gate = spec_gate(samples)
    ran = {(s["id"], s["fluid"]) for s in samples or []}
    scores, credited = [], []
    for p in env.peptides:
        for fluid in ("sgf", "sif"):
            est = answers.get(p.cid, {}).get(fluid)
            score = 0.0
            if est is not None:
                error = abs(np.log10(np.clip(est, 1.0, 960.0)) - np.log10(np.clip(gi.true_half_life(p, fluid), 1.0, 960.0)))
                score = max(0.0, 1.0 - error / 0.5)
            scores.append(score)
            credited.append(score if gate and (p.cid, fluid) in ran else score / 2)
    return float(np.mean(scores)), gate, float(np.mean(credited)) - 0.15 * cost / 500.0


def test_reward_and_gate_arithmetic_match_the_spec_text():
    text = make(NAME).spec()["reward"]
    for phrase in ("1 - |log10 half-life error| / 0.5", "clipped to 1 to 960 min", "no answer scoring 0",
                   "at specification enzyme activity", "halved unless that compound was run in that fluid and the run "
                   "passes its QC gate", "2+ replicates sampled at time 0", "is <= 15%", "0.8 to 1.25",
                   "analysis.enzyme_activity", "log-linear regression", "chem.enzyme_exposure", "at or under 2%",
                   "minus 0.15 * cost / 500", "10 + 2 per LC-MS sample + 0.1 per tip + 20 per fresh enzyme preparation"):
        assert phrase in text, phrase
    gates = set()
    for policy in POLICIES[NAME]:
        for seed in range(6):
            seen, actions, reward, info, env = walk(policy, seed)
            samples = seen[1].get("samples") if len(seen) > 2 else None
            accuracy, gate, expected = spec_reward(env, samples, actions[-1]["half_life_min"], info["cost"])
            assert info["qc_pass"] == gate, (policy, seed)
            assert np.isclose(info["accuracy"], accuracy), (policy, seed)
            assert np.isclose(reward, expected), (policy, seed)
            assert reward <= 1.0
            gates.add(gate)
    assert gates == {True, False}  # both branches of the gate were exercised


def test_the_cost_is_the_one_the_spec_states():
    for policy in ("expert", "stale_enzyme", "random"):
        seen, _, _, info, env = walk(policy, 2)
        samples = seen[1]["samples"]
        tips = sum(s["op"] == "transfer" for s in env.step_log)  # every transfer takes a fresh tip
        fresh = len({s["fluid"] for s in samples}) if any("made fresh" in str(s.get("what")) for s in env.step_log) else 0
        assert np.isclose(info["cost"], 10 + 2 * len(samples) + 0.1 * tips + 20 * fresh), policy


def test_a_guessed_pair_never_counts_in_full():
    # one arm run well enough to pass the gate, the true value reported for every pair
    env = make(NAME, seed=1)
    token = {**GI_EXPERT, "timepoints_min": [0, 5, 30], "fluids": {"cpd1": ["sgf"]}}
    env.step({"type": "run", "design": token})
    truth = {p.cid: {f: gi.true_half_life(p, f) for f in gi.FLUIDS} for p in env.peptides}
    _, reward, _, info = env.step({"type": "report", "half_life_min": truth})
    assert info["qc_pass"] and info["accuracy"] == 1.0 and info["pairs_run"] == ["cpd1/sgf"]
    assert np.isclose(reward, (1.0 + 5 * 0.5) / 6 - 0.15 * info["cost"] / 500.0)
    # and a token run does not lift constant guesses for the pairs it did not measure
    def mean_reward(run_first):
        total = []
        for seed in range(30):
            env = make(NAME, seed=seed)
            guess = {p.cid: {"sgf": 60.0, "sif": 10.0} for p in env.peptides}
            if run_first:
                obs = env.step({"type": "run", "design": {**token, "fresh_enzyme": False}})[0]
                guess["cpd1"]["sgf"] = GiFluidPolicy(token).act(obs)["half_life_min"]["cpd1"]["sgf"]
            total.append(env.step({"type": "report", "half_life_min": guess})[1])
        return np.mean(total)

    assert mean_reward(True) < 0.3 and mean_reward(True) - mean_reward(False) < 1.0 / 6


def changes(trap):
    """Design fields and analysis flags in which a policy differs from the expert."""
    expert, policy = POLICIES[NAME]["expert"](0), POLICIES[NAME][trap](0)
    assert set(policy.design) == set(GI_EXPERT) and policy.fluids == expert.fluids and policy.pick is None, trap
    return ([k for k in GI_EXPERT if policy.design[k] != GI_EXPERT[k]]
            + [k for k in expert.flags if policy.flags[k] != expert.flags[k]])


def test_each_trap_is_the_expert_with_exactly_one_change():
    assert 4 <= len(SINGLE_TRAPS) <= 7 and {"expert", "null", "random"} <= set(POLICIES[NAME])
    expert = POLICIES[NAME]["expert"](0)
    assert expert.design == GI_EXPERT and all(expert.flags.values())
    assert gi.check_design({**GI_EXPERT, "fluids": BOTH}, list(BOTH))[1] is None
    changed = set()
    for trap in SINGLE_TRAPS:
        assert len(changes(trap)) == 1, (trap, changes(trap))
        changed.add(changes(trap)[0])
    assert len(changed) == len(SINGLE_TRAPS)  # and no two traps make the same mistake
    assert POLICIES[NAME]["null"](0).design is None


def test_an_interaction_is_a_single_mistake_plus_the_warm_overnight_wait_and_is_labelled_as_one():
    assert set(GI_INTERACTIONS) == set(TRAPS) - set(SINGLE_TRAPS) and all(GI_EXPERT[k] != v for k, v in GI_WARM_OVERNIGHT.items())
    for both, single in GI_INTERACTIONS.items():
        assert single in SINGLE_TRAPS
        assert sorted(changes(both)) == sorted(changes(single) + list(GI_WARM_OVERNIGHT))
        assert POLICIES[NAME][both](0).design == {**POLICIES[NAME][single](0).design, **GI_WARM_OVERNIGHT}
        rule = [r for r in knowledge.load_rules() if r["env"] == NAME and both in r["traps"]]
        assert len(rule) == 1 and single in rule[0]["traps"]  # filed under the rule its single mistake breaks


def test_every_trap_has_its_own_rule():
    rules = [r for r in knowledge.load_rules() if r["env"] == NAME]
    assert knowledge.orphan_traps(knowledge.load_rules()) == []
    assert sorted(t for r in rules for t in r["traps"]) == sorted(TRAPS)
    assert all(len([t for t in r["traps"] if t in SINGLE_TRAPS]) == 1 and r["status"] == "draft" for r in rules)
    card = envqc.qc_card(NAME, n_seeds=4)
    rows = knowledge.coverage(rules, {NAME: card})
    assert [r["trap"] for r in rows] == [t for r in rules for t in r["traps"]]


def test_expert_beats_null_by_a_clear_margin_on_held_out_seeds():
    seeds = range(30000, 30060)
    expert = np.mean([play("expert", s)[0] for s in seeds])
    null = np.mean([play("null", s)[0] for s in seeds])
    assert expert >= 0.6 and null <= 0.3 and expert > null + 0.5


def test_one_episode_takes_well_under_50_ms():
    best = np.inf
    for seed in range(5):
        start = time.perf_counter()
        play("expert", seed)
        best = min(best, time.perf_counter() - start)
    assert best < 0.05


# ---------------------------------------------------------------- the environment's own rules

def test_spec_documents_every_design_field_and_is_json():
    env = make(NAME, seed=0)
    spec = env.reset()["spec"]
    json.dumps(spec)
    fields = spec["actions"]["run"]["design"]
    assert set(fields) == set(DESIGN) and set(fields["sgf_quench"]["options"]) == set(gi.QUENCHES)
    assert spec["activity_control"] == {"id": "control", "reference_half_life_min": 15.0, "acceptance_window": [0.8, 1.25]}
    assert spec["fluids"]["sif"]["enzyme_activity_half_life_37C_min"] == 240.0
    assert env.reset()["compounds"] == ["cpd1", "cpd2", "cpd3"]


def test_constructor_and_labware_limits():
    for n in (0, 8):
        with pytest.raises(ValueError):
            GiFluidStability(seed=0, n_compounds=n)
    with pytest.raises(ValueError):
        make(NAME, seed=-1)
    env = GiFluidStability(seed=0, n_compounds=7)
    everyone = {f"cpd{i + 1}": ["sgf", "sif"] for i in range(7)}
    error = env.step({"type": "run", "design": {**GI_EXPERT, "fluids": everyone}})[0]["error"]
    assert "quench plate is full" in error
    fewer = {**GI_EXPERT, "fluids": everyone, "replicates": 1}
    assert "error" not in env.step({"type": "run", "design": fewer})[0]


def test_only_the_chosen_arms_are_run_and_charged():
    env = make(NAME, seed=2)
    chosen = {"cpd1": ["sgf"], "cpd2": ["sif", "sgf"], "cpd3": []}
    obs = env.step({"type": "run", "design": {**GI_EXPERT, "fluids": chosen}})[0]
    arms = {(s["id"], s["fluid"]) for s in obs["samples"]}
    assert arms == {("cpd1", "sgf"), ("cpd2", "sgf"), ("cpd2", "sif"), ("control", "sgf"), ("control", "sif")}
    n, tips = len(obs["samples"]), 5 * 2 + len(obs["samples"])  # one tip per spike and per aliquot
    assert n == 5 * 2 * 8 and np.isclose(env.cost, 10.0 + 2.0 * n + 0.1 * tips + 2 * 20.0)
    # a pair that was not run can still be answered, but an answer with no measurement behind it counts at half
    answer = {p.cid: {f: gi.true_half_life(p, f) for f in gi.FLUIDS} for p in env.peptides}
    _, reward, _, info = env.step({"type": "report", "half_life_min": answer})
    assert info["accuracy"] == 1.0 and info["qc_pass"] and info["pairs_run"] == ["cpd1/sgf", "cpd2/sgf", "cpd2/sif"]
    assert np.isclose(reward, (3 * 1.0 + 3 * 0.5) / 6 - 0.15 * env.cost / 500.0)


def test_fresh_enzyme_is_what_costs_more():
    costs = {}
    for fresh in (True, False):
        env = make(NAME, seed=0)
        env.step({"type": "run", "design": {**DESIGN, "fresh_enzyme": fresh}})
        costs[fresh] = env.cost
    assert np.isclose(costs[True], 10.0 + 2.0 * 128 + 0.1 * 144 + 40.0) and np.isclose(costs[True] - costs[False], 40.0)


def test_worklist_names_every_reagent_and_the_hold():
    env = make(NAME, seed=0)
    run_episode(env, POLICIES[NAME]["expert"](0))
    fills = [s for s in env.step_log if s["op"] == "fill"]
    assert all(s.get("what") for s in fills)
    assert {"SGF with pepsin, made fresh for this run", "SIF with pancreatin, made fresh for this run", "cpd1 stock",
            "control substrate stock", "base quench with internal standard",
            "organic quench with internal standard"} <= {s["what"] for s in fills}
    holds = [s for s in env.step_log if s["op"] == "hold"]
    assert [h["dst_plate"] for h in holds] == ["sgf quench", "sif quench"] and holds[0]["minutes"] == 30.0
    rows = list(csv.DictReader(io.StringIO(worklist_csv(env.step_log))))
    assert len(rows) == len(env.step_log) and {r["run"] for r in rows} == {"1"}
    assert sum(r["op"] == "transfer" and r["dst_plate"].endswith("quench") for r in rows) == 128


def test_the_exported_worklist_tells_the_runs_apart():
    def exported(**changes):
        env = make(NAME, seed=0)
        env.step({"type": "run", "design": {**DESIGN, **changes}})
        return list(csv.DictReader(io.StringIO(worklist_csv(env.step_log))))

    def whats(rows, op):
        return {r["what"] for r in rows if r["op"] == op}

    expert = exported()
    assert whats(expert, "hold") == {"quenched samples wait for injection at 4 C"}
    assert whats(exported(hold_temp_c=25.0, hold_min=960.0), "hold") == {"quenched samples wait for injection at 25 C"}
    stale = whats(exported(fresh_enzyme=False), "fill")
    assert {"SGF with pepsin, made 24 h ago and kept at 4 C", "SIF with pancreatin, made 24 h ago and kept at 4 C"} <= stale
    assert not any("made fresh" in w for w in stale)
    # a pre-quenched time zero is named on its quench fills: one per incubation well, 4 arms x 2 replicates a fluid
    zero = "base quench with internal standard, time zero quenched before the enzyme acts"
    assert sum(r["what"] == zero for r in expert) == 8
    late_zero = exported(prequenched_t0=False)
    assert not any("time zero quenched" in r["what"] for r in late_zero)
    assert [r["minutes"] for r in late_zero if r["op"] == "incubate"][0] == "1.0"  # it digests for the handling time first


def test_qc_gate_needs_replicates_a_time_zero_and_a_control_in_window():
    def gate(**changes):
        env = make(NAME, seed=1)  # a seed whose fresh lots sit inside the window
        env.step({"type": "run", "design": {**DESIGN, **changes}})
        return env.step({"type": "report", "half_life_min": {}})[3]

    assert gate()["qc_pass"]
    assert not gate(replicates=1)["qc_pass"]
    assert not gate(timepoints_min=[2, 5, 10, 20, 45])["qc_pass"]
    assert not gate(activity_control=False)["qc_pass"]
    stale = gate(fresh_enzyme=False)
    assert not stale["qc_pass"] and stale["enzyme_activity"]["sif"] < 0.8 <= stale["enzyme_activity"]["sgf"]
    assert not make(NAME, seed=1).step({"type": "report", "half_life_min": {}})[3]["qc_pass"]  # nothing was run


# ---------------------------------------------------------------- chemistry

def test_rate_laws():
    assert q10_factor(37.0) == 1.0 and np.isclose(q10_factor(27.0), 0.5) and np.isclose(q10_factor(4.0, q10=2.0), 2 ** -3.3)
    assert np.isclose(acid_catalysis_factor(1.2, 1.2), 1.0) and np.isclose(acid_catalysis_factor(2.2, 1.2), 0.1)
    assert enzyme_exposure(10.0, 25.0) == 15.0 and enzyme_exposure(0.0, 30.0, float("inf")) == 30.0
    half_life = 240.0
    assert np.isclose(enzyme_exposure(0.0, 1e6, half_life), half_life / np.log(2.0))  # all the enzyme will ever do
    assert np.isclose(enzyme_exposure(0.0, 60.0, half_life) + enzyme_exposure(60.0, 240.0, half_life),
                      enzyme_exposure(0.0, 240.0, half_life))
    assert enzyme_exposure(0.0, 240.0, half_life) < 240.0


def test_truth_is_the_half_life_at_specification_activity_with_acid_loss_counted():
    stable = gi.Peptide("x", 100.0, 40.0, None)
    labile = gi.Peptide("y", 100.0, 40.0, 100.0)
    assert np.isclose(gi.true_half_life(stable, "sgf"), 100.0) and np.isclose(gi.true_half_life(stable, "sif"), 40.0)
    assert np.isclose(gi.true_half_life(labile, "sgf"), 50.0)  # two equal first-order routes
    assert np.isclose(gi.true_half_life(labile, "sif"), 40.0, rtol=1e-4)  # no acid at pH 6.8
    assert np.isclose(gi.stored_activity("sif", 24.0), remaining(240.0 / 2 ** -3.3, 1440.0))
    assert 0.6 < gi.stored_activity("sif", 24.0) < 0.7 < 0.9 < gi.stored_activity("sgf", 24.0) < 1.0
    rng = np.random.default_rng(0)
    drawn = [gi.draw_peptide(rng, f"c{i}") for i in range(2000)]
    assert 0.25 < np.mean([p.acid_min is not None for p in drawn]) < 0.35
    assert np.mean([gi.true_half_life(p, "sif") < 5.0 for p in drawn]) > 0.25  # many are gone within minutes


# ---------------------------------------------------------------- analysis

def clean_samples(cid, fluid, half_life, times, activity=1.0, enzyme_half_life=None, reps=1):
    return [{"id": cid, "fluid": fluid, "replicate": r + 1, "time_min": float(t),
             "area": 1000.0 * remaining(half_life / activity, enzyme_exposure(0.0, t, enzyme_half_life))}
            for t in times for r in range(reps)]


def test_fluid_half_life_is_fitted_per_fluid():
    times = [0, 2, 5, 10, 20, 45, 120, 240]
    samples = clean_samples("c", "sgf", 80.0, times) + clean_samples("c", "sif", 3.0, times)
    est = estimate_fluid_half_life(samples)
    assert np.isclose(est["c"]["sgf"], 80.0, rtol=0.01) and np.isclose(est["c"]["sif"], 3.0, rtol=0.01)
    assert estimate_fluid_half_life(clean_samples("c", "sgf", 1e6, times))["c"]["sgf"] == 960.0  # stable: upper limit
    assert estimate_fluid_half_life(clean_samples("c", "sif", 0.2, times))["c"]["sif"] == 1.0  # gone: lower limit
    assert estimate_fluid_half_life(clean_samples("c", "sif", 30.0, [0])) == {}  # no time course, no answer


def test_exposure_time_corrects_for_an_enzyme_that_weakens_during_the_run():
    times = [0, 2, 5, 10, 20, 45, 120, 240]
    samples = clean_samples("c", "sif", 300.0, times, enzyme_half_life=240.0)
    naive = estimate_fluid_half_life(samples)["c"]["sif"]
    corrected = estimate_fluid_half_life(samples, {"sif": 240.0})["c"]["sif"]
    assert np.isclose(corrected, 300.0, rtol=0.01) and naive > 1.25 * corrected
    # a fluid the mapping does not name is fitted against the clock
    assert estimate_fluid_half_life(samples, {"sgf": 1440.0})["c"]["sif"] == naive


def test_the_control_substrate_reads_activity_and_scales_the_answer():
    times = [0, 2, 5, 10, 20, 45, 120, 240]
    samples = (clean_samples("control", "sif", 15.0, times, activity=0.6)
               + clean_samples("c", "sif", 50.0, times, activity=0.6) + clean_samples("c", "sgf", 50.0, times))
    activity = enzyme_activity(samples, "control", 15.0)
    assert set(activity) == {"sif"} and np.isclose(activity["sif"], 0.6, rtol=0.01)
    raw = estimate_fluid_half_life(samples)["c"]
    fixed = estimate_fluid_half_life(samples, relative_activity=activity)["c"]
    assert np.isclose(raw["sif"], 50.0 / 0.6, rtol=0.01) and np.isclose(fixed["sif"], 50.0, rtol=0.01)
    assert fixed["sgf"] == raw["sgf"]  # no control in that fluid, so no correction
    for useless in (0.0, float("nan"), -1.0):
        assert estimate_fluid_half_life(samples, relative_activity={"sif": useless})["c"]["sif"] == raw["sif"]
    dead = clean_samples("control", "sgf", 1e9, times)
    assert enzyme_activity(dead, "control", 15.0)["sgf"] < 1e-3  # a control that does not move: no activity
    assert enzyme_activity(samples, "nobody", 15.0) == {}


def test_lloq_is_applied_per_time_point_and_everything_after_is_dropped():
    times = [0, 2, 5, 10, 20, 45]
    samples = clean_samples("c", "sif", 2.0, times, reps=2)
    lucky = [dict(s) for s in samples]
    for s in lucky:  # one replicate at 20 min and both at 45 min read far above the noise by accident
        if (s["time_min"], s["replicate"]) in ((20.0, 1), (45.0, 1), (45.0, 2)):
            s["area"] = 30.0
    assert np.isclose(estimate_fluid_half_life(samples)["c"]["sif"], 2.0, rtol=0.01)
    assert np.isclose(estimate_fluid_half_life(lucky)["c"]["sif"], 2.0, rtol=0.01)
    assert estimate_fluid_half_life(lucky, use_lloq=False)["c"]["sif"] > 4.0  # the noise points flatten the slope
    # only time zero is quantifiable: the first dropped time point gives an upper bound
    bound = estimate_fluid_half_life(clean_samples("c", "sif", 1.0, [0, 15, 30, 60]))["c"]["sif"]
    assert np.isclose(bound, 15.0 / np.log2(50.0))
    with_is = [{**s, "is_area": 0.0} for s in samples]
    assert estimate_fluid_half_life(with_is) == {}  # no usable response at all


# ---------------------------------------------------------------- mechanisms in the simulator

FAST = gi.Peptide("cpd1", 5.0, 5.0, None)


def test_an_acid_quench_leaves_pepsin_running_and_stops_pancreatin():
    warm = {"hold_temp_c": 25.0, "hold_min": 240.0}
    base, _ = simulate([FAST], sgf_quench="base", sif_quench="acid", **warm)
    acid, _ = simulate([FAST], sgf_quench="acid", sif_quench="acid", **warm)
    wrong_way, _ = simulate([FAST], sgf_quench="base", sif_quench="base", **warm)
    assert t0_response(acid, "cpd1", "sgf") < 0.05 * t0_response(base, "cpd1", "sgf")  # digested while it waited
    assert np.isclose(t0_response(acid, "cpd1", "sif"), t0_response(base, "cpd1", "sif"), rtol=0.15)
    assert t0_response(wrong_way, "cpd1", "sif") < 0.05 * t0_response(base, "cpd1", "sif")  # base does not stop pancreatin
    # with a quench that works, the same warm wait costs nothing
    cold, _ = simulate([FAST], sgf_quench="base", sif_quench="acid")
    assert np.isclose(t0_response(base, "cpd1", "sgf"), t0_response(cold, "cpd1", "sgf"), rtol=0.15)


def test_the_wrong_quench_is_survivable_cold_and_ruinous_warm_overnight():
    def mean_reward(policy):
        return np.mean([run_episode(make(NAME, seed=s), policy)[0] for s in range(40)])

    expert, cold = mean_reward(POLICIES[NAME]["expert"](0)), mean_reward(POLICIES[NAME]["acid_quench_in_sgf"](0))
    warm_overnight = mean_reward(POLICIES[NAME]["acid_quench_in_sgf_warm_overnight"](0))
    right_quench_warm_overnight = mean_reward(GiFluidPolicy({**GI_EXPERT, **GI_WARM_OVERNIGHT}))
    assert expert - cold < 0.15 and expert - right_quench_warm_overnight < 0.15  # neither half alone
    assert expert - warm_overnight > 0.4  # the two together


def test_a_cold_plate_slows_what_the_quench_did_not_stop():
    cold, _ = simulate([FAST], sgf_quench="acid", hold_temp_c=4.0, hold_min=240.0)
    warm, _ = simulate([FAST], sgf_quench="acid", hold_temp_c=25.0, hold_min=240.0)
    assert t0_response(warm, "cpd1", "sgf") < 0.2 * t0_response(cold, "cpd1", "sgf")


def test_an_acid_labile_compound_keeps_hydrolysing_in_an_acidic_sample():
    labile = gi.Peptide("cpd1", 2000.0, 2000.0, 30.0)
    wait = {"hold_temp_c": 25.0, "hold_min": 960.0}
    organic, _ = simulate([labile], sgf_quench="organic", sif_quench="organic", **wait)
    base, _ = simulate([labile], sgf_quench="base", sif_quench="organic", **wait)
    acid, _ = simulate([labile], sgf_quench="base", sif_quench="acid", **wait)
    assert t0_response(organic, "cpd1", "sgf") < 0.5 * t0_response(base, "cpd1", "sgf")  # still at pH 1.9
    assert t0_response(acid, "cpd1", "sif") < 0.8 * t0_response(organic, "cpd1", "sif")  # acidified to pH 2.5
    est = estimate_fluid_half_life(base, gi.ACTIVITY_HL_37C_MIN)["cpd1"]
    assert est["sgf"] < 40.0 < 500.0 < est["sif"]  # acid loss is part of the SGF half-life and absent in SIF


def test_stale_enzyme_reads_long_and_only_the_control_shows_it():
    medium = gi.Peptide("cpd1", 20.0, 20.0, None)
    activity, raw, fixed = [], [], []
    for seed in range(6):
        samples, _ = simulate([medium], seed=seed, fresh_enzyme=False)
        a = enzyme_activity(samples, "control", gi.CONTROL_HALF_LIFE_MIN, gi.ACTIVITY_HL_37C_MIN)
        compounds = [s for s in samples if s["id"] != "control"]
        activity.append(a["sif"])
        raw.append(estimate_fluid_half_life(compounds, gi.ACTIVITY_HL_37C_MIN)["cpd1"]["sif"])
        fixed.append(estimate_fluid_half_life(compounds, gi.ACTIVITY_HL_37C_MIN, a)["cpd1"]["sif"])
    assert np.isclose(np.mean(activity), gi.stored_activity("sif", gi.STALE_AGE_H), rtol=0.05)
    assert np.isclose(np.mean(raw), 20.0 / gi.stored_activity("sif", gi.STALE_AGE_H), rtol=0.06)  # 1.5 times too long
    assert np.isclose(np.mean(fixed), 20.0, rtol=0.06)
    without, _ = simulate([medium], fresh_enzyme=False, activity_control=False)
    assert "control" not in {s["id"] for s in without}
    assert enzyme_activity(without, "control", gi.CONTROL_HALF_LIFE_MIN) == {}  # nothing in the data shows it


def test_a_fresh_lot_off_specification_is_read_by_the_control():
    samples, _ = simulate([FAST], lot={"sgf": 0.7, "sif": 1.3})
    activity = enzyme_activity(samples, "control", gi.CONTROL_HALF_LIFE_MIN, gi.ACTIVITY_HL_37C_MIN)
    assert np.isclose(activity["sgf"], 0.7, rtol=0.08) and np.isclose(activity["sif"], 1.3, rtol=0.08)
    assert not gi.run_qc(samples)["qc_pass"]


def test_a_time_zero_taken_after_the_spike_has_already_digested():
    one_minute = gi.Peptide("cpd1", 1.0, 1.0, None)
    ratio = []
    for seed in range(6):
        true_zero, _ = simulate([one_minute], seed=seed)
        late_zero, rig = simulate([one_minute], seed=seed, prequenched_t0=False)
        ratio.append(t0_response(late_zero, "cpd1", "sif") / t0_response(true_zero, "cpd1", "sif"))
        assert rig.log[[s["op"] for s in rig.log].index("incubate")]["minutes"] == gi.HANDLING_MIN
    assert np.isclose(np.mean(ratio), 0.5, atol=0.08)  # one half-life gone before "time zero"


def test_a_late_first_point_only_bounds_a_fast_peptide():
    one_minute = gi.Peptide("cpd1", 1.0, 1.0, None)
    early, _ = simulate([one_minute])
    late, _ = simulate([one_minute], timepoints_min=[0, 15, 30, 45, 60, 90, 120, 240])
    fit = dict(activity_half_life_min=gi.ACTIVITY_HL_37C_MIN)
    assert np.isclose(estimate_fluid_half_life(early, **fit)["cpd1"]["sif"], 1.0, rtol=0.2)
    assert estimate_fluid_half_life(late, **fit)["cpd1"]["sif"] > 2.0


def test_the_enzyme_weakening_during_the_run_makes_slow_peptides_read_long_against_the_clock():
    slow = gi.Peptide("cpd1", 3000.0, 300.0, None)
    naive, corrected = [], []
    for seed in range(8):
        samples, _ = simulate([slow], seed=seed, activity_control=False)
        naive.append(estimate_fluid_half_life(samples)["cpd1"]["sif"])
        corrected.append(estimate_fluid_half_life(samples, gi.ACTIVITY_HL_37C_MIN)["cpd1"]["sif"])
    assert np.isclose(np.mean(corrected), 300.0, rtol=0.12) and np.mean(naive) > 1.2 * np.mean(corrected)
    # half the pancreatin and most of the pepsin are still there when a four-hour run ends: nothing ran out
    assert np.isclose(remaining(gi.ACTIVITY_HL_37C_MIN["sif"], gi.HORIZON_MIN), 0.5)
    assert remaining(gi.ACTIVITY_HL_37C_MIN["sgf"], gi.HORIZON_MIN) > 0.85


def test_no_internal_standard_fails_the_time_zero_gate_most_of_the_time():
    passed = [play("no_internal_standard", s)[1]["qc_pass"] for s in range(20)]
    assert np.mean(passed) < 0.3
