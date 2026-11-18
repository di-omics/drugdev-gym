import copy
import json
import time
from pathlib import Path

import numpy as np
import pytest

from drugdevgym import analysis, knowledge, make
from drugdevgym.analysis import estimate_cc50, fit_cc50, normalize_viability, viability_zprime
from drugdevgym.chem import dissolved, grown, haber_shift, survival
from drugdevgym.envs import cytotox_window as cw
from drugdevgym.policies import CX_CHOICES, CX_EXPERT, POLICIES, CytotoxPolicy, random_cx_design, run_episode
from drugdevgym.rig import EDGE, RigProfile
from drugdevgym.sim2real import worklist_csv

NAME = "cytotox_window"
ROOT = Path(__file__).resolve().parents[1]
TRAPS = [p for p in POLICIES[NAME] if p not in ("expert", "random", "null")]
TRUTH_WORDS = ("truth", "cc50", "kinetics", "solubility", "luciferase", "hill", "ceiling", "window_log10", "slow",
               "sticky", "cell_cv")


def play(policy, seed, rig="nominal"):
    return run_episode(make(NAME, seed=seed, rig=rig), POLICIES[NAME][policy](seed))


def compound(cid, cc50, hill=1.5, slow=False, solubility=1e4, luc=None, sticky=1.0):
    return cw.Compound(cid, 0.01, cc50, hill, slow, solubility, luc, sticky)


def run(compounds, design=None, seed=0, profile=None):
    wells, rig, plates = cw.simulate_run(profile or RigProfile(), np.random.default_rng(seed), compounds, 0.05,
                                         {**CX_EXPERT, **(design or {})})
    return wells, rig, plates


def spec_rule(said_um, at_least, true_um, ceiling_um, tested_um):
    """The per-compound rule, written again from the words in spec()["reward"]. tested_um None = no run.

    T, S, the "no run" line and the cap of 0.5 were added after the first card, following review.
    """
    x, big_l, c = np.log10(said_um), np.log10(true_um), np.log10(ceiling_um)
    if not at_least:
        return max(-0.5, 1 - abs(x - big_l) / 0.5)
    if tested_um is None:
        return -0.5
    t = np.log10(tested_um)
    b = min(big_l, c)
    s = min(b, t)
    score = max(0.0, 1 - (b - x) / 0.5) if x <= s else max(-0.5, max(0.0, 1 - (b - s) / 0.5) - (x - s) / 0.5)
    return min(score, 0.5) if big_l < min(c, t) else score


# ---------------------------------------------------------------- the contract

@pytest.mark.parametrize("policy", ["expert", "random", "edge_wells", "unmatched_dmso"])
def test_same_seed_gives_the_same_episode(policy):
    assert play(policy, 7) == play(policy, 7)
    a, b = make(NAME, seed=7), make(NAME, seed=7)
    assert a.reset() == b.reset()
    action = {"type": "run", "design": CX_EXPERT}
    assert a.step(action) == b.step(action)


def test_truth_depends_on_the_seed_and_not_on_the_policy():
    truths = [play(p, 9)[1]["truth"] for p in POLICIES[NAME]]
    assert all(t == truths[0] for t in truths)
    assert play("expert", 9)[1]["truth"] != play("expert", 10)[1]["truth"]
    assert cw.CytotoxWindow.salted
    env = make(NAME, seed=9)  # salted: not the first draw an unsalted environment makes on this seed
    assert env.cell_cv != float(np.random.default_rng([9, 1]).uniform(*cw.CELL_CV))


@pytest.mark.parametrize("policy", list(POLICIES[NAME]))
@pytest.mark.parametrize("rig", ["nominal", "harsh"])
def test_no_truth_field_names_in_any_observation(policy, rig):
    for seed in (5, 6):
        env, pol = make(NAME, seed=seed, rig=rig), POLICIES[NAME][policy](seed)
        obs, done = env.reset(), False
        while not done:
            text = json.dumps({k: v for k, v in obs.items() if k != "spec"}).lower()
            assert not any(word in text for word in TRUTH_WORDS), (policy, seed)
            obs, reward, done, info = env.step(pol.act(obs))
            assert done or info == {}
        assert set(info["truth"]["cpd1"]) == {"cc50_uM", "window_log10", "ceiling_uM", "solubility_uM", "kinetics",
                                              "viability_hill", "luciferase_ic50_uM"}
        assert obs == {"phase": "done"} and reward <= 1.0


def test_the_spec_gives_no_hidden_value_away():
    env = make(NAME, seed=3)
    text = json.dumps(env.reset())
    for c in env.compounds:
        for hidden in (c.cc50_um, c.solubility_um, c.hill, c.sticky, c.luc_ic50_um):
            assert hidden is None or repr(hidden) not in text
    assert json.dumps(make(NAME, seed=3).spec()) == json.dumps(make(NAME, seed=4).spec())


HUGE = 10 ** 400
BAD_NUMBERS = (float("nan"), float("inf"), -float("inf"), HUGE, -HUGE, [8], [[8]], None, "8", {}, True)
BAD_FLAGS = ("yes", 1, 0, None, [True], float("nan"))


def bad_designs():
    for field, (kind, *_rest) in cw.DESIGN_SPEC.items():
        for bad in (BAD_FLAGS if kind == "bool" else BAD_NUMBERS):
            yield {**CX_EXPERT, field: bad}
    yield {**CX_EXPERT, "extra": 1}
    yield {k: v for k, v in CX_EXPERT.items() if k != "lid"}
    yield {**CX_EXPERT, "n_points": 8.5}
    yield {**CX_EXPERT, "exposure_h": 36}
    yield {**CX_EXPERT, "top_conc_uM": 100.5}
    yield {**CX_EXPERT, "top_conc_uM": 300.0}  # allowed between the second and the third pass only
    yield {**CX_EXPERT, "untreated_controls": 5}  # a field of the second pass, removed in the third
    yield {**CX_EXPERT, "top_conc_uM": 0.5}
    yield {**CX_EXPERT, 5: 1}
    yield {**CX_EXPERT, HUGE ** 20: 1}  # a key too long to print
    yield from (None, [], "design", 5, [CX_EXPERT], {})


def test_malformed_designs_are_rejected_without_exceptions():
    for design in bad_designs():
        env = make(NAME, seed=1)
        obs, reward, done, info = env.step({"type": "run", "design": design})
        assert "error" in obs and reward == 0.0 and not done and info == {} and env.cost == 0.0, design
        assert env.wells is None and env.step_log == []
    env = make(NAME, seed=1)
    assert "error" in env.step({"type": "run", "design": CX_EXPERT, "also": 1})[0]
    assert "error" in env.step({"type": "run"})[0]


def bad_reports():
    ok = {"cc50_uM": 10.0, "at_least": True}
    yield from (None, [], "cpd1", 3, [["cpd1", ok]])
    yield {"cpd9": ok}
    yield {1: ok}
    yield {HUGE ** 20: ok}
    for answer in (None, [], 10.0, "10", [ok], {}, {"cc50_uM": 10.0}, {"at_least": True}, {**ok, "window": 2.0},
                   {**ok, 7: 1}):
        yield {"cpd1": answer}
    for bad in (*BAD_NUMBERS, 0.0, -3.0, 1e7):
        yield {"cpd1": {**ok, "cc50_uM": bad}}
    for bad in BAD_FLAGS:
        yield {"cpd1": {**ok, "at_least": bad}}


def test_malformed_reports_are_rejected_without_exceptions():
    for results in bad_reports():
        env = make(NAME, seed=1)
        obs, reward, done, info = env.step({"type": "report", "results": results})
        assert "error" in obs and reward == 0.0 and not done and info == {}, results
    env = make(NAME, seed=1)
    assert "error" in env.step({"type": "report"})[0]
    assert "error" in env.step({"type": "report", "results": {}, "plates_used": [1]})[0]
    assert "error" in env.step({"type": "titrate"})[0]


def test_five_invalid_actions_end_the_episode_and_one_run_is_all_there_is():
    env = make(NAME, seed=2)
    for _ in range(4):
        assert not env.step({"type": "run", "design": {}})[2]
    obs, reward, done, info = env.step({"type": "run", "design": {}})
    assert done and info == {"truncated": True} and reward == 0.0
    env = make(NAME, seed=2)
    assert env.step({"type": "run", "design": CX_EXPERT})[0]["phase"] == "data"
    cost = env.cost
    assert "error" in env.step({"type": "run", "design": CX_EXPERT})[0] and env.cost == cost
    _, _, done, info = env.step({"type": "report", "results": {}})
    assert done and info["quality"] == 0.0
    with pytest.raises(RuntimeError):
        env.step({"type": "report", "results": {}})


def test_designs_that_do_not_fit_a_plate_are_rejected():
    env = make(NAME, seed=0)
    crowded = {**CX_EXPERT, "n_points": 12, "replicates": 3, "vehicle_controls": 16, "kill_controls": 16}
    assert "a plate has 60" in env.step({"type": "run", "design": crowded})[0]["error"]
    assert cw.plan_plates(CX_EXPERT, 6)[0] == [[0, 1, 2, 3], [4, 5]]
    assert cw.plan_plates({**CX_EXPERT, "avoid_edges": False}, 6)[0] == [[0, 1, 2, 3, 4, 5]]
    assert cw.plan_plates({**CX_EXPERT, "replicates": 2, "n_points": 8}, 6)[0] == [[0, 1], [2, 3], [4, 5]]
    assert cw.check_design({k: v[0] for k, v in CX_CHOICES.items()}, 6)[1] is None
    # the largest choice of every field still fits the inner wells, so every random design is valid as drawn
    worst = {k: v[-1] for k, v in CX_CHOICES.items()}
    assert worst["avoid_edges"] and cw.check_design(worst, 6) == (worst, None)
    assert max(CX_CHOICES["top_conc_uM"]) == 100.0 and "untreated_controls" not in CX_CHOICES
    for seed in range(300):
        design = random_cx_design(np.random.default_rng([seed, 99]))
        assert cw.check_design(design, 6) == (design, None)
        assert POLICIES[NAME]["random"](seed).design == design
    with pytest.raises(ValueError):
        cw.CytotoxWindow(seed=0, n_compounds=9)


def test_editing_an_observation_cannot_change_the_score():
    design = {**CX_EXPERT, "vehicle_controls": 3, "kill_controls": 3}
    plain, edited = make(NAME, seed=4), make(NAME, seed=4)
    first = edited.reset()
    first["compounds"].clear()
    first["spec"]["reward"] = "always 1"
    first["spec"]["actions"]["run"]["design"].clear()
    obs_plain, *_ = plain.step({"type": "run", "design": design})
    obs, *_ = edited.step({"type": "run", "design": design})
    answer = {"type": "report", "results": estimate_cc50(obs["wells"])}
    for w in obs["wells"]:  # turn every control into a perfect one, then throw the plate away
        w["lum"] = 1e6 if w["role"] == "vehicle" else 0.0
        w["role"] = "vehicle"
    obs["wells"].clear()
    assert edited.reset()["compounds"] == plain.reset()["compounds"] and len(edited.reset()["compounds"]) == 6
    assert edited.spec() == plain.spec() and len(edited.wells) == len(obs_plain["wells"])
    assert edited.step(copy.deepcopy(answer)) == plain.step(copy.deepcopy(answer))
    # the gate is computed from the environment's own copy: a failing plate stays failing
    env = make(NAME, seed=4)
    obs, *_ = env.step({"type": "run", "design": {**CX_EXPERT, "kill_controls": 0}})
    for w in obs["wells"][:6]:
        w["role"], w["lum"] = "kill", 0.0
    assert not env.step({"type": "report", "results": {}})[3]["qc_pass"]
    env = make(NAME, seed=4)  # nor can an unmatched run be relabelled as matched
    obs, *_ = env.step({"type": "run", "design": {**CX_EXPERT, "match_dmso": False}})
    for w in obs["wells"]:
        w["dmso_pct"] = 1.0
    info = env.step({"type": "report", "results": {}})[3]
    assert not info["qc_pass"] and info["vehicle_matched"] == [False, False]


def test_score_claim_is_the_rule_in_the_spec():
    def said(x, at_least, true, ceiling, tested=100.0):
        return cw.score_claim(x, at_least, true, ceiling, tested)

    assert said(20.0, False, 20.0, 100.0) == 1.0
    assert np.isclose(said(20.0 * 10 ** 0.25, False, 20.0, 100.0), 0.5)
    assert np.isclose(said(20.0 * 10 ** 0.5, False, 20.0, 100.0), 0.0, atol=1e-9)
    assert said(2.0, False, 500.0, 100.0) == -0.5  # a false alarm
    assert said(90.0, False, 5.0, 100.0) == -0.5  # an understated liability
    assert said(100.0, True, 500.0, 100.0) == 1.0  # clean: the bound at the ceiling, tested
    assert said(30.0, True, 500.0, 30.0) == 1.0  # clean, but only 30 uM could be dissolved
    assert np.isclose(said(100.0, True, 500.0, 30.0), max(-0.5, 1 - np.log10(100 / 30) / 0.5))  # above solubility
    assert said(100.0, True, 5.0, 100.0) == -0.5  # a missed liability
    assert said(0.001, True, 500.0, 100.0) == 0.0  # a very safe bound earns nothing, and costs nothing
    assert said(100.0, False, 3000.0, 100.0) == -0.5  # "clean" given as a number
    assert said(20.0, False, 20.0, 100.0, None) == 1.0  # a measured value is scored against the truth alone
    rng = np.random.default_rng(0)
    for _ in range(2000):
        x, true, ceiling, tested = (float(10 ** rng.uniform(-1, 3.5)) for _ in range(4))
        for at_least in (False, True):
            for t in (tested, None):
                got = cw.score_claim(x, at_least, true, ceiling, t)
                assert np.isclose(got, spec_rule(x, at_least, true, ceiling, t)) and -0.5 <= got <= 1.0


def test_a_bound_is_scored_against_what_the_run_tested():
    """Added after the first card, following review: a bound above the highest dose tested was paid in full."""
    clean = dict(true_cc50_um=2000.0, ceiling_um=100.0)
    honest = cw.score_claim(30.0, True, tested_um=30.0, **clean)
    assert np.isclose(honest, 0.0, atol=1e-9)  # under-dosing still earns little: 30 uM is 0.52 log under the ceiling
    assert cw.score_claim(100.0, True, tested_um=100.0, **clean) == 1.0
    for claim in (31.0, 40.0, 55.0, 100.0, 300.0):  # claiming more than was tested never beats the honest bound
        assert cw.score_claim(claim, True, tested_um=30.0, **clean) < honest
    assert cw.score_claim(100.0, True, tested_um=30.0, **clean) == -0.5
    assert np.isclose(cw.score_claim(70.0, True, tested_um=50.0, **clean),
                      (1 - np.log10(100 / 50) / 0.5) - np.log10(70 / 50) / 0.5)  # the credit at 50 uM, less the excess
    # the function is defined for any tested dose, though no design can dose above 100 uM: a tested, true bound
    # above the ceiling is still above the ceiling
    assert np.isclose(cw.score_claim(300.0, True, tested_um=300.0, **clean), 1 - np.log10(3.0) / 0.5)
    assert cw.score_claim(100.0, True, tested_um=300.0, **clean) == 1.0
    for claim in (1e-6, 1.0, 100.0, 1e6):  # with no run nothing was tested
        assert cw.score_claim(claim, True, tested_um=None, **clean) == -0.5
    xs = 10 ** np.linspace(0, 3, 601)  # no jump where the claim crosses what was tested
    for tested in (10.0, 30.0, 100.0, 300.0):
        scores = np.array([cw.score_claim(float(x), True, tested_um=tested, **clean) for x in xs])
        best = cw.score_claim(min(tested, 100.0), True, tested_um=tested, **clean)  # the honest bound is the best one
        assert np.abs(np.diff(scores)).max() < 0.011 and scores.max() <= best + 1e-9
    # in an episode: a run to 30 uM, then every bound raised to "at least 100 uM"
    low = {**CX_EXPERT, "top_conc_uM": 30.0}
    rewards = {"honest": [], "raised": []}
    for seed in range(30100, 30130):
        for name in rewards:
            env = make(NAME, seed=seed)
            obs, *_ = env.step({"type": "run", "design": low})
            results = estimate_cc50(obs["wells"])
            if name == "raised":
                results = {k: {"cc50_uM": 100.0, "at_least": True} if v["at_least"] else v for k, v in results.items()}
            reward, info = env.step({"type": "report", "results": results})[1::2]
            assert info["tested_top_uM"] == 30.0
            rewards[name].append(reward)
    assert np.mean(rewards["honest"]) > np.mean(rewards["raised"]) + 0.3
    assert all(h >= r for h, r in zip(rewards["honest"], rewards["raised"]))


def test_a_bound_label_is_not_free_on_a_compound_that_lost_viability():
    """Added after the first card, following review: relabelling a measured liability as a bound cost nothing."""
    assert cw.HEDGE_CAP == 0.5
    liable = dict(true_cc50_um=20.0, ceiling_um=100.0)
    assert np.isclose(cw.score_claim(18.0, False, tested_um=100.0, **liable), 1 - np.log10(20 / 18) / 0.5)
    assert cw.score_claim(18.0, True, tested_um=100.0, **liable) == 0.5  # true, close, and capped
    assert cw.score_claim(20.0, True, tested_um=100.0, **liable) == 0.5
    assert np.isclose(cw.score_claim(8.0, True, tested_um=100.0, **liable), 1 - np.log10(20 / 8) / 0.5)  # under the cap
    assert np.isclose(cw.score_claim(100.0, True, tested_um=100.0, **liable), 1 - np.log10(5.0) / 0.5)  # still a miss
    # no cap when the run stopped under the CC50, or the CC50 is above the ceiling: nothing was seen to die
    assert np.isclose(cw.score_claim(10.0, True, tested_um=10.0, **liable), 1 - np.log10(2.0) / 0.5)
    assert cw.score_claim(18.0, True, tested_um=18.0, **liable) > 0.9
    assert cw.score_claim(100.0, True, 150.0, 100.0, 300.0) == 1.0
    assert cw.score_claim(25.0, True, 60.0, 30.0, 100.0) > 0.8  # CC50 above the solubility limit
    lost = []  # in an episode: the expert's own answers, every one relabelled as a bound
    for seed in range(30100, 30140):
        rewards = []
        for relabel in (False, True):
            env = make(NAME, seed=seed)
            obs, *_ = env.step({"type": "run", "design": CX_EXPERT})
            results = estimate_cc50(obs["wells"])
            if relabel:
                results = {k: {**v, "at_least": True} for k, v in results.items()}
            rewards.append(env.step({"type": "report", "results": results})[1])
        lost.append(rewards[0] - rewards[1])
    assert np.mean(lost) > 0.05 and np.mean(np.array(lost) >= 0) > 0.9


def test_the_spec_states_the_numbers_the_code_uses():
    text = make(NAME).spec()["reward"]
    for phrase in ("max(-0.5, 1 - |X - L| / 0.5)", "B = min(L, C)", "S = min(B, T)", "max(0, 1 - (B - X) / 0.5)",
                   "max(-0.5, max(0, 1 - (B - S) / 0.5) - (X - S) / 0.5)",
                   "the lower of 100 uM (the dose limit, the highest dose the stock reaches at 1% DMSO) and the "
                   "compound's solubility limit",
                   "T is the highest nominal dose the run put on cells", "every bound scores -0.5",
                   "(L < min(C, T))", "scores at most 0.5",
                   "at least 3 vehicle and 3 full-kill wells and a Z' of at least 0.5 between their luminescence",
                   "By convention: a solvent-matched vehicle",
                   "the vehicle wells carry the DMSO percentage of the top dose on that plate",
                   "their nominal dmso_pct equals that of the plate's highest-dose wells",
                   "what a cytotoxicity guideline asks of a solvent control",
                   "It is a rule, not a simulated effect",
                   "a run with match_dmso false leaves its vehicle wells without DMSO and cannot pass the gate",
                   "a quality above 0 is halved", "minus 0.15 * cost / 600",
                   "10 per assay plate, 1 per assay well (controls and cell-free wells included) and 0.1 per tip",
                   "the DMSO ladder plates and the dosing plates are not charged",
                   "2 * n_points - 1 tips per compound", "an unmatched run n_points",
                   "one tip per dosed assay well, cell-free wells included",
                   "left out scores 0", "true CC50 after 72 h"):
        assert phrase in text, phrase
    assert (cw.TOLERANCE, cw.PENALTY_FLOOR, cw.ZPRIME_GATE, cw.MIN_CONTROLS) == (0.5, -0.5, 0.5, 3)
    assert (cw.PLATE_COST, cw.WELL_COST, cw.TIP_COST, cw.CytotoxWindow.cost_reference) == (10.0, 1.0, 0.1, 600.0)
    assert cw.MATCHED_DMSO_PCT == 1.0 and cw.TOP_LIMIT_UM == 100.0 and cw.REFERENCE_EXPOSURE_H == 72.0
    assert cw.HEDGE_CAP == 0.5  # added after the first card, following review
    # added in the second pass and removed in the third: no constant, design field or info key is left behind
    for gone in ("MAX_DMSO_PCT", "VEHICLE_HEALTH_GATE", "STANDARD_TOP_UM", "STANDARD_DMSO_PCT", "CEILING_UM",
                 "top_dmso_pct"):
        assert not hasattr(cw, gone), gone
    assert "untreated_controls" not in cw.DESIGN_SPEC and cw.DESIGN_SPEC["top_conc_uM"][2] == 100.0


def test_what_changed_after_the_first_card_is_said_where_the_agent_and_the_reader_look():
    text = make(NAME).spec()["reward"]
    assert "Changed after the first card, following two adversarial reviews" in text
    for phrase in ("a bound is scored against what the run tested", "the cap of 0.5 on a bound",
                   "the solvent-matched vehicle term of the gate are new",
                   "untreated wells, a vehicle-health term of the gate and doses up to 300 uM",
                   "they were measured to do no job and were removed",
                   "the dose limit is 100 uM as on the first card",
                   "one well per dose and 6 dose points instead of 2 and 8",
                   "the 5 untreated wells per plate it carried between the second and the third pass are gone"):
        assert phrase in text, phrase
    doc = (ROOT / "docs" / "ASSUMPTIONS.md").read_text()
    section = doc[doc.index("## cytotox_window: what was fixed first"):]
    for phrase in ("After an adversarial review", "After the re-review", "HEDGE_CAP", "seeds 20000 to 20199",
                   "solvent-matched vehicle", "punished by the gate alone", "as no_kill_control is",
                   "With the term switched off", "VEHICLE_HEALTH_GATE", "MAX_DMSO_PCT", "untreated_controls",
                   "were removed", "unmatched_dmso", "edge_wells"):
        assert phrase in section, phrase
    # the reference design itself: changed after the first card (2 wells per dose and 8 points before); the
    # 5 untreated wells of the second pass are gone again
    assert (CX_EXPERT["replicates"], CX_EXPERT["n_points"]) == (1, 6) and "untreated_controls" not in CX_EXPERT
    assert (CX_EXPERT["vehicle_controls"], CX_EXPERT["kill_controls"]) == (5, 5)
    assert (CX_EXPERT["top_conc_uM"], CX_EXPERT["match_dmso"], CX_EXPERT["exposure_h"]) == (100.0, True, 72)


@pytest.mark.parametrize("policy", ["expert", "no_kill_control", "edge_wells", "exact_not_bound", "unmatched_dmso",
                                    "random", "null"])
def test_reward_and_gate_arithmetic_match_the_spec_text(policy):
    for seed in range(6):
        env = make(NAME, seed=seed)
        pol = POLICIES[NAME][policy](seed)
        reward, info = run_episode(env, pol)
        tested = None if pol.design is None else pol.design["top_conc_uM"]
        assert info["tested_top_uM"] == tested
        scores = []
        for cid, truth in info["truth"].items():
            said = info["reported"][cid]
            ceiling = min(100.0, truth["solubility_uM"])
            assert truth["ceiling_uM"] == ceiling
            scores.append(0.0 if said is None else
                          spec_rule(said["cc50_uM"], said["at_least"], truth["cc50_uM"], ceiling, tested))
            assert np.isclose(info["scores"][cid], scores[-1], atol=1e-3)
            if said is not None:
                assert np.isclose(said["window_log10"] - truth["window_log10"],
                                  np.log10(said["cc50_uM"] / truth["cc50_uM"]), atol=1e-3)
        quality = float(np.mean(scores))
        wells = env.wells or []
        plates = sorted({w["plate"] for w in wells})
        passed = bool(plates)
        for i, p in enumerate(plates):
            on = [w for w in wells if w["plate"] == p]
            vehicle, kill = (np.array([w["lum"] for w in on if w["role"] == role]) for role in ("vehicle", "kill"))
            # the convention: the vehicle wells carry the DMSO percentage of the top dose on that plate
            top = max(w["conc_uM"] for w in on if w["role"] == "cpd")
            top_dmso = {w["dmso_pct"] for w in on if w["role"] == "cpd" and w["conc_uM"] == top}
            matched = len(vehicle) > 0 and {w["dmso_pct"] for w in on if w["role"] == "vehicle"} == top_dmso
            assert info["vehicle_matched"][i] == matched == pol.design["match_dmso"] or len(vehicle) == 0
            if len(vehicle) < 3 or len(kill) < 3:
                passed = False
                continue
            z = 1 - 3 * (vehicle.std(ddof=1) + kill.std(ddof=1)) / abs(vehicle.mean() - kill.mean())
            passed = passed and z >= 0.5 and matched
        # the cost line: assay plates, assay wells, and the tips the spec says each protocol takes
        tips = sum(1 for s in env.step_log if s["op"] == "transfer")
        if pol.design is not None:
            d = pol.design
            before_stamps = 2 * d["n_points"] - 1 if d["match_dmso"] else d["n_points"]
            dosed = sum(w["role"] in ("cpd", "cell_free") for w in wells)
            assert tips == 6 * before_stamps + dosed == 6 * cw.tips_per_compound(d)
        cost = 10.0 * len(plates) + 1.0 * len(wells) + 0.1 * tips
        assert np.isclose(info["cost"], cost) and info["qc_pass"] == passed
        assert np.isclose(info["quality"], quality)
        assert np.isclose(reward, (quality if passed or quality <= 0 else quality / 2) - 0.15 * cost / 600.0)


def test_the_gate_halves_a_good_answer_and_never_lifts_a_bad_one():
    reward, info = play("no_kill_control", 3)
    assert not info["qc_pass"] and info["zprime"] == [None, None] and info["quality"] > 0.5
    assert info["vehicle_matched"] == [True, True]  # it fails on the missing full-kill wells alone
    assert np.isclose(reward, 0.5 * info["quality"] - 0.15 * info["cost"] / 600.0)
    env = make(NAME, seed=3)  # no run at all: the gate fails, and a quality below zero is not halved
    wrong = {c["id"]: {"cc50_uM": 1e-3, "at_least": False} for c in env.reset()["compounds"]}
    _, reward, _, info = env.step({"type": "report", "results": wrong})
    assert not info["qc_pass"] and info["zprime"] == [] and info["quality"] == -0.5 and reward == -0.5
    assert info["vehicle_matched"] == [] and info["tested_top_uM"] is None
    reward, info = play("expert", 3)
    assert info["qc_pass"] and all(z >= 0.5 for z in info["zprime"]) and info["vehicle_matched"] == [True, True]
    assert np.isclose(reward, info["quality"] - 0.15 * info["cost"] / 600.0)


def plate_map(vehicle_dmso, dosed_dmso, plate=1, vehicle_lum=20000.0, kill_lum=100.0):
    """Hand-made control and dosed wells of one plate; dosed_dmso lists (conc_uM, dmso_pct), any order."""
    rows = [{"plate": plate, "role": "vehicle", "conc_uM": None, "dmso_pct": d, "lum": vehicle_lum + 300.0 * i}
            for i, d in enumerate(vehicle_dmso)]
    rows += [{"plate": plate, "role": "kill", "conc_uM": None, "dmso_pct": d, "lum": kill_lum + 5.0 * i}
             for i, d in enumerate(vehicle_dmso)]
    rows += [{"plate": plate, "role": "cpd", "conc_uM": c, "dmso_pct": d, "lum": 15000.0} for c, d in dosed_dmso]
    return rows


def test_the_gate_asks_for_a_solvent_matched_vehicle():
    """Added after the first card, in the third pass (the re-review's fix): a plate passes its gate only if its
    vehicle wells carry the DMSO percentage of the top dose on that plate. A convention, not a simulated effect."""
    matched_ladder = [(100.0, 1.0), (50.0, 1.0), (25.0, 1.0)]
    diluted_ladder = [(25.0, 0.25), (100.0, 1.0), (50.0, 0.5)]
    ok, z, matched = cw.gate(plate_map([1.0, 1.0, 1.0], matched_ladder))
    assert ok and z[0] > 0.8 and matched == [True]
    ok, z, matched = cw.gate(plate_map([0.0, 0.0, 0.0], diluted_ladder))  # a tight Z' does not buy the plate back
    assert not ok and z[0] > 0.8 and matched == [False]
    assert cw.gate(plate_map([1.0, 1.0, 1.0], diluted_ladder))[0]  # what counts is the vehicle against the top dose
    assert not cw.gate(plate_map([0.5, 0.5, 0.5], diluted_ladder))[0]  # ... not against a lower one
    assert not cw.gate(plate_map([1.0, 1.0, 0.0], matched_ladder))[0]  # every vehicle well, not most of them
    assert not cw.gate(plate_map([1.0, 1.0, 1.0], [(100.0, 0.5), (50.0, 0.5)]))[0]  # more DMSO than the top dose
    two = plate_map([1.0, 1.0, 1.0], matched_ladder) + plate_map([0.0, 0.0, 0.0], diluted_ladder, plate=2)
    assert cw.gate(two) == (False, cw.gate(two)[1], [True, False])  # one unmatched plate fails the run
    assert not cw.vehicle_is_matched([]) and not cw.vehicle_is_matched(plate_map([1.0] * 3, []))
    assert not cw.vehicle_is_matched(plate_map([], matched_ladder))
    assert cw.gate([]) == (False, [], []) and cw.gate(None) == (False, [], [])
    # it reads the plate map, never the luminescence: dead vehicle wells still count as matched (and fail Z')
    ok, z, matched = cw.gate(plate_map([1.0, 1.0, 1.0], matched_ladder, vehicle_lum=100.0))
    assert not ok and matched == [True]
    # in an episode, on every rig: a matched run passes this term and an unmatched run cannot
    for rig in ("nominal", "randomized", "harsh"):
        for seed in range(30300, 30306):
            for match in (True, False):
                env = make(NAME, seed=seed, rig=rig)
                obs, *_ = env.step({"type": "run", "design": {**CX_EXPERT, "match_dmso": match}})
                assert {w["dmso_pct"] for w in obs["wells"] if w["role"] == "vehicle"} == {1.0 if match else 0.0}
                assert max(w["dmso_pct"] for w in obs["wells"] if w["role"] == "cpd") == 1.0
                _, reward, _, info = env.step({"type": "report", "results": estimate_cc50(obs["wells"])})
                assert info["vehicle_matched"] == [match, match]
                if not match:
                    assert not info["qc_pass"]
                    gated = info["quality"] if info["quality"] <= 0 else 0.5 * info["quality"]
                    assert np.isclose(reward, gated - 0.15 * info["cost"] / 600.0)


def test_unmatched_dmso_is_punished_by_the_gate_alone(monkeypatch):
    """The honest reading of the unmatched_dmso row: with the convention switched off in process the trap is
    level with the expert or ahead of it, because at 1% DMSO the mismatch costs little and the route saves tips."""
    seeds = range(30300, 30360)
    expert = np.array([play("expert", s)[0] for s in seeds])
    as_coded = [play("unmatched_dmso", s) for s in seeds]
    trap = np.array([r for r, _ in as_coded])
    assert all(not info["qc_pass"] and min(info["zprime"]) >= 0.5 for _, info in as_coded)  # Z' is fine: the rule fails
    assert expert.mean() > trap.mean() + 0.3 and np.mean(expert > trap) > 0.95
    quality_expert = np.mean([play("expert", s)[1]["quality"] for s in seeds])
    quality_trap = np.mean([info["quality"] for _, info in as_coded])
    assert abs(quality_expert - quality_trap) < 0.02  # the answers are as good either way
    monkeypatch.setattr(cw, "vehicle_is_matched", lambda on_plate: True)
    physics = [play("unmatched_dmso", s) for s in seeds]
    assert all(info["qc_pass"] for _, info in physics)
    free = np.array([r for r, _ in physics])
    assert abs(expert.mean() - free.mean()) < 0.02 and np.mean(expert > free) < 0.6  # no simulated effect to see


def test_each_trap_differs_from_the_expert_in_exactly_one_choice():
    expert = POLICIES[NAME]["expert"](0)
    assert expert.design == CX_EXPERT and expert.use_turbidity and expert.bound_when_clean
    assert 4 <= len(TRAPS) <= 7
    changed = []
    for trap in TRAPS:
        pol = POLICIES[NAME][trap](0)
        assert isinstance(pol, CytotoxPolicy) and set(pol.design) == set(CX_EXPERT)
        diffs = [k for k in CX_EXPERT if pol.design[k] != CX_EXPERT[k]]
        diffs += [flag for flag in ("use_turbidity", "bound_when_clean") if getattr(pol, flag) != getattr(expert, flag)]
        assert len(diffs) == 1, (trap, diffs)
        changed.append(diffs[0])
    assert len(set(changed)) == len(TRAPS)  # and no two traps make the same mistake
    assert cw.check_design(CX_EXPERT, 6) == (CX_EXPERT, None)


def test_every_trap_has_a_rule_and_the_null_runs_nothing():
    rules = [r for r in knowledge.load_rules() if r["env"] == NAME]
    assert sorted(t for r in rules for t in r["traps"]) == sorted(TRAPS)
    assert all(r["status"] == "draft" and len(r["traps"]) == 1 for r in rules)
    assert knowledge.orphan_traps(knowledge.load_rules()) == []
    env = make(NAME, seed=0)
    reward, info = run_episode(env, POLICIES[NAME]["null"](0))
    assert env.step_log == [] and info["cost"] == 0.0 and not info["qc_pass"]


def test_the_expert_beats_the_null_on_fresh_seeds_by_a_clear_margin():
    seeds = range(30000, 30060)
    expert = np.array([play("expert", s)[0] for s in seeds])
    null = np.array([play("null", s)[0] for s in seeds])
    assert expert.mean() > 0.6 and null.mean() < 0.3
    assert expert.mean() > null.mean() + 0.4 and np.mean(expert > null) > 0.95


def test_no_constant_answer_and_no_safe_bound_scores_near_the_expert():
    seeds = range(30000, 30040)
    expert = np.mean([play("expert", s)[0] for s in seeds])
    token = {**CX_EXPERT, "n_points": 4, "replicates": 1, "interference_control": False, "vehicle_controls": 3,
             "kill_controls": 3, "avoid_edges": False, "exposure_h": 24}

    def constant(value, at_least, design=None):
        rewards = []
        for s in seeds:
            env = make(NAME, seed=s)
            ids = [c["id"] for c in env.reset()["compounds"]]
            if design is not None:
                env.step({"type": "run", "design": design})
            rewards.append(env.step({"type": "report", "results": {i: {"cc50_uM": value, "at_least": at_least}
                                                                   for i in ids}})[1])
        return float(np.mean(rewards))

    for value in (1e-6, 1.0, 10.0, 30.0, 50.0, 100.0, 300.0, 1e6):
        # with no run every bound is unsupported (changed after the first card: "at least 100 uM" used to earn 0.11)
        assert constant(value, True) == -0.5 and constant(value, False) < 0.1, value
    env = make(NAME, seed=30000)  # with no run, the best report is the empty one
    assert env.step({"type": "report", "results": {}})[1] == 0.0
    # a token run that only buys the QC gate, then the best constant: still far from the expert
    assert constant(1e-6, True, token) < 0.0 < constant(100.0, True, token)  # the always-safe bound earns nothing
    assert max(constant(v, True, token) for v in (50.0, 100.0)) < 0.35 < expert - 0.35


def test_every_fill_is_labelled_and_the_run_is_recorded():
    for design in (CX_EXPERT, {**CX_EXPERT, "match_dmso": False}):
        env = make(NAME, seed=0)
        env.step({"type": "run", "design": design})
        fills = [s for s in env.step_log if s["op"] == "fill"]
        assert fills and all(s["run"] == 1 for s in env.step_log)
        labels = {s["what"] for s in fills}
        assert not labels & {"buffer", "DMSO+cpd1", "ATP"} and all(isinstance(x, str) and len(x) > 3 for x in labels)
        assert {"medium", "ATP detection reagent", "ATP standard in medium, no cells", "cells, 4k per well"} <= labels
        assert any(x.startswith("cpd1 at 10000 uM in DMSO") for x in labels)
        assert any(x.startswith("vehicle") for x in labels) and any(x.startswith("full-kill") for x in labels)
        assert {s["op"] for s in env.step_log} == {"fill", "transfer", "incubate", "read"}
        assert len(worklist_csv(env.step_log).splitlines()) == len(env.step_log) + 1
        incubations = [s for s in env.step_log if s["op"] == "incubate"]
        assert len(incubations) == len({w["plate"] for w in env.wells}) <= 2
        assert all(s["minutes"] == 72 * 60 and s["lid"] for s in incubations)
        vehicle = "vehicle: medium with 10% DMSO" if design["match_dmso"] else "vehicle: medium"
        assert vehicle in labels and not any("untreated" in x for x in labels)


def test_one_episode_takes_well_under_fifty_milliseconds():
    play("expert", 0)
    times = []
    for seed in range(1, 11):
        start = time.perf_counter()
        play("expert", seed)
        times.append(time.perf_counter() - start)
    assert np.median(times) < 0.05


def test_every_constant_has_a_row_in_the_assumptions():
    text = (ROOT / "docs" / "ASSUMPTIONS.md").read_text()
    names = [n for n, v in vars(cw).items() if n.isupper() and isinstance(v, (int, float, tuple))]
    names += ["TURBIDITY_SD_MULT", "INTERFERENCE_FLAG", "MIN_READOUT_FRACTION", "cost_reference"]
    assert len(names) > 40
    borrowed = {"EDGE", "ROWS", "COST_WEIGHT"}  # from rig.py and base.py, not new here
    missing = [n for n in names if n not in borrowed and n not in text]
    assert missing == []


# ---------------------------------------------------------------- the physics, on hand-made compounds

def viability(wells, cid):
    """Mean normalised viability per dose, top dose first."""
    rows = [w for w in wells if w["id"] == cid and w["role"] == "cpd"]
    plate = rows[0]["plate"]
    vehicle = [w["lum"] for w in wells if w["plate"] == plate and w["role"] == "vehicle"]
    kill = [w["lum"] for w in wells if w["plate"] == plate and w["role"] == "kill"]
    doses = sorted({w["conc_uM"] for w in rows}, reverse=True)
    return [float(np.mean(normalize_viability([w["lum"] for w in rows if w["conc_uM"] == d], vehicle, kill)))
            for d in doses]


def test_the_layout_and_nominal_doses_follow_the_design():
    wells, rig, plates = run([compound(f"c{i}", 1e4) for i in range(6)])
    assert plates == 2 and len(wells) == (10 + 2 + 4 * 12) + (10 + 2 + 2 * 12) and rig.short_aspirates == 0
    first = [w for w in wells if w["plate"] == 1]
    assert [w["role"] for w in first[:12]] == ["vehicle", "kill"] * 5 + ["cell_free_ref"] * 2
    assert {w["role"] for w in wells} == {"vehicle", "kill", "cell_free_ref", "cpd", "cell_free"}
    assert sorted({w["conc_uM"] for w in first if w["id"] == "c0"}, reverse=True) == [100 / 2 ** j for j in range(6)]
    assert sum(w["role"] == "cpd" for w in first) == 24 and sum(w["role"] == "cell_free" for w in first) == 24
    assert {w["dmso_pct"] for w in wells} == {1.0}
    assert set(wells[0]) == {"plate", "well", "role", "id", "conc_uM", "dmso_pct", "lum", "turbidity"}
    wells, _, plates = run([compound(f"c{i}", 1e4) for i in range(6)], {"match_dmso": False, "avoid_edges": False})
    assert plates == 1 and {w["dmso_pct"] for w in wells if w["role"] in ("vehicle", "kill", "cell_free_ref")} == {0.0}
    assert sorted({w["dmso_pct"] for w in wells if w["id"] == "c0"}, reverse=True)[:3] == [1.0, 0.5, 0.25]
    wells, _, plates = run([compound(f"c{i}", 1e4) for i in range(6)], {"replicates": 2, "n_points": 8})
    assert plates == 3 and len(wells) == 3 * (10 + 2 + 2 * 24)  # the first-card design: two compounds a plate


def test_a_measurable_liability_is_measured_and_a_clean_compound_gets_a_bound():
    wells, *_ = run([compound("toxic", 20.0), compound("clean", 2000.0)])
    est = estimate_cc50(wells)
    assert not est["toxic"]["at_least"] and abs(np.log10(est["toxic"]["cc50_uM"] / 20.0)) < 0.15
    assert est["clean"] == {"cc50_uM": 100.0, "at_least": True}
    assert estimate_cc50(wells, bound_when_clean=False)["clean"] == {"cc50_uM": 100.0, "at_least": False}
    assert all(z > 0.5 for z in viability_zprime(wells).values())


def test_a_slow_compound_is_missed_after_one_day():
    pair = [compound("slow", 20.0, slow=True), compound("fast", 20.0)]
    day, _, _ = run(pair, {"exposure_h": 24})
    three, _, _ = run(pair, {"exposure_h": 72})
    assert estimate_cc50(day)["slow"] == {"cc50_uM": 100.0, "at_least": True}
    assert not estimate_cc50(three)["slow"]["at_least"]
    assert abs(np.log10(estimate_cc50(day)["fast"]["cc50_uM"] / estimate_cc50(three)["fast"]["cc50_uM"])) < 0.15
    missed = cw.score_claim(100.0, True, 20.0, 100.0, 100.0)  # a missed liability: -0.4
    assert np.isclose(missed, 1 - np.log10(5.0) / 0.5)


def test_unmatched_dmso_lowers_only_the_top_wells_and_both_protocols_use_the_same_transfers_per_dose():
    six = [compound(f"c{i}", 1e5) for i in range(6)]
    matched = np.mean([viability(run(six, seed=s)[0], "c0")[:2] for s in range(8)], axis=0)
    unmatched = np.mean([viability(run(six, {"match_dmso": False}, seed=s)[0], "c0")[:2] for s in range(8)], axis=0)
    assert abs(matched[0] - 1.0) < 0.05 and 0.82 < unmatched[0] < 0.95 and unmatched[1] > 0.93
    assert np.isclose(cw.dmso_survival(1.0), 1 - 1 / 9) and cw.dmso_survival(2.0) == 0.5 and cw.dmso_survival(0.0) == 1.0
    tips = {m: run(six, {"match_dmso": m})[1].tips for m in (True, False)}
    stamps = 6 * 6 * 2
    assert tips[True] == stamps + 6 * (5 + 6) and tips[False] == stamps + 6 * (5 + 1)
    assert cw.tips_per_compound(CX_EXPERT) == 23 and cw.tips_per_compound({**CX_EXPERT, "match_dmso": False}) == 18


def test_the_dose_limit_is_100_uM_and_the_top_well_holds_1_percent_dmso():
    """The second pass let a design dose to 300 uM with up to 3% DMSO and gated on vehicle health against
    untreated wells. The third pass measured that neither did a job and removed both: the first-card limit is back."""
    assert cw.TOP_LIMIT_UM == 100.0 == cw.STOCK_MM * 1000.0 * cw.MATCHED_DMSO_PCT / 100.0
    for top in (100.5, 150.0, 300.0):
        assert "top_conc_uM must be between" in cw.check_design({**CX_EXPERT, "top_conc_uM": top}, 6)[1]
    six = [compound(f"c{i}", 1e5) for i in range(6)]  # nothing here is toxic or insoluble
    for top in (10.0, 30.0, 100.0):  # the stock is pre-diluted in DMSO, so the top well holds 1% at any top dose
        wells, rig, _ = run(six, {"top_conc_uM": top})
        assert rig.short_aspirates == 0 and {w["dmso_pct"] for w in wells} == {1.0}
        assert sorted({w["conc_uM"] for w in wells if w["id"] == "c0"})[-1] == top
        assert cw.gate(wells)[0] and abs(np.mean(viability(wells, "c0")) - 1.0) < 0.1
        assert estimate_cc50(wells)["c0"] == {"cc50_uM": top, "at_least": True}
        assert any(s["op"] == "fill" and s["what"] == "vehicle: medium with 10% DMSO" for s in rig.log)
        wells, rig, _ = run(six, {"top_conc_uM": top, "match_dmso": False})
        assert rig.short_aspirates == 0 and not cw.gate(wells)[0] and cw.gate(wells)[2] == [False, False]
        assert sorted({w["dmso_pct"] for w in wells if w["id"] == "c0"}, reverse=True)[:2] == [1.0, 0.5]
    env = make(NAME, seed=0)
    assert set(env.step({"type": "report", "results": {}})[3]) == {
        "quality", "qc_pass", "zprime", "vehicle_matched", "tested_top_uM", "cost", "scores", "reported", "truth"}


def test_the_reference_design_is_the_one_the_reward_supports():
    """Changed after the first card, following review: the first design (2 wells per dose, 8 points) was dominated."""
    seeds = range(30200, 30240)
    now = np.array([play("expert", s)[0] for s in seeds])
    first_card = {**CX_EXPERT, "replicates": 2, "n_points": 8}
    before = np.array([run_episode(make(NAME, seed=s), CytotoxPolicy(first_card))[0] for s in seeds])
    assert now.mean() > before.mean() + 0.02 and np.mean(now > before) > 0.7
    for change in ({"replicates": 2}, {"n_points": 8}):  # one step back toward it does not pay either
        other = np.array([run_episode(make(NAME, seed=s), CytotoxPolicy({**CX_EXPERT, **change}))[0] for s in seeds])
        assert now.mean() > other.mean()


def test_edge_wells_lose_volume_and_cells_over_three_days():
    six = [compound(f"c{i}", 1e5) for i in range(6)]
    wells, *_ = run(six, {"avoid_edges": False, "vehicle_controls": 8, "kill_controls": 8})
    name = lambda r, c: "ABCDEFGH"[r] + str(c + 1)  # noqa: E731
    edge = {name(r, c) for r in range(8) for c in range(12) if EDGE[r, c]}
    vehicle = [w for w in wells if w["role"] == "vehicle"]
    on_edge = np.mean([w["lum"] for w in vehicle if w["well"] in edge])
    inner = np.mean([w["lum"] for w in vehicle if w["well"] not in edge])
    assert on_edge < 0.85 * inner
    assert np.isclose(cw.osmotic_survival(1.4), 0.5) and cw.osmotic_survival(1.0) == 1.0 == cw.osmotic_survival(0.9)
    open_plate, *_ = run(six, {"lid": False})
    closed, *_ = run(six)
    lum = lambda ws: np.mean([w["lum"] for w in ws if w["role"] == "vehicle"])  # noqa: E731
    assert lum(open_plate) < 0.9 * lum(closed)


def test_confluent_wells_hide_toxicity():
    pair = [compound("toxic", 20.0), compound("clean", 2000.0)]
    sparse = estimate_cc50(run(pair, {"cells_k": 4.0})[0])["toxic"]["cc50_uM"]
    crowded = estimate_cc50(run(pair, {"cells_k": 20.0})[0])["toxic"]["cc50_uM"]
    assert crowded > 2 * sparse
    assert grown(4000.0, 72.0, 24.0, 50000.0) == 32000.0 and grown(20000.0, 72.0, 24.0, 50000.0) == 50000.0


def test_turbidity_shows_precipitate_and_the_bound_stops_at_the_last_clear_dose():
    wells, *_ = run([compound("crashes", 2000.0, solubility=30.0), compound("clear", 2000.0)])
    top = [w["turbidity"] for w in wells if w["id"] == "crashes" and w["role"] == "cpd" and w["conc_uM"] == 100.0]
    low = [w["turbidity"] for w in wells if w["id"] == "crashes" and w["role"] == "cpd" and w["conc_uM"] == 25.0]
    assert min(top) > 1300 and max(low) < 650
    assert estimate_cc50(wells)["crashes"] == {"cc50_uM": 25.0, "at_least": True}
    assert estimate_cc50(wells, use_turbidity=False)["crashes"] == {"cc50_uM": 100.0, "at_least": True}
    assert cw.score_claim(25.0, True, 2000.0, 30.0, 100.0) > 0.8
    assert cw.score_claim(100.0, True, 2000.0, 30.0, 100.0) < 0.0
    assert float(dissolved(100.0, 30.0)) == 30.0 and float(dissolved(10.0, 30.0)) == 10.0


def test_readout_interference_looks_like_toxicity_until_the_cell_free_control_is_used():
    pair = [compound("inhibitor", 2000.0, luc=10.0), compound("clean", 2000.0)]
    wells, *_ = run(pair)
    cell_free = sorted((w["conc_uM"], w["lum"]) for w in wells if w["id"] == "inhibitor" and w["role"] == "cell_free")
    ref = np.mean([w["lum"] for w in wells if w["role"] == "cell_free_ref"])
    assert cell_free[-1][1] < 0.15 * ref and cell_free[0][1] > 0.7 * ref  # 9% left at 100 uM, 76% at 3.1 uM
    assert estimate_cc50(wells)["inhibitor"] == {"cc50_uM": 100.0, "at_least": True}
    without, *_ = run(pair, {"interference_control": False})
    fooled = estimate_cc50(without)["inhibitor"]
    assert not fooled["at_least"] and 5.0 < fooled["cc50_uM"] < 20.0
    assert np.isclose(cw.luciferase_activity(10.0, 10.0), 0.5)


# ---------------------------------------------------------------- analysis and chemistry functions

def test_chemistry_functions():
    assert survival(0.0, 5.0, 2.0) == 1.0 and np.isclose(survival(5.0, 5.0, 2.0), 0.5)
    assert np.isclose(survival(10.0, 5.0, 2.0), 0.2)
    assert haber_shift(72.0, 72.0, 2.0) == 1.0 and np.isclose(haber_shift(48.0, 72.0, 2.0), 2.25)
    assert np.isclose(haber_shift(24.0, 72.0, 2.0), 9.0) and haber_shift(24.0, 72.0, 0.0) == 1.0
    assert np.allclose(dissolved([5.0, 50.0, -1.0], 20.0), [5.0, 20.0, 0.0])
    assert np.allclose(grown([1000.0, 10000.0], 48.0, 24.0, 30000.0), [4000.0, 30000.0])


def test_normalize_viability_uses_both_controls():
    assert np.allclose(normalize_viability([100.0, 550.0, 1000.0], [1000.0, 1000.0], [100.0, 100.0]), [0.0, 0.5, 1.0])
    assert np.allclose(normalize_viability([500.0], [1000.0]), [0.5])  # no full-kill wells: zero signal is the zero
    assert np.isnan(normalize_viability([500.0], [100.0], [200.0])).all()
    assert np.isnan(normalize_viability([500.0], [])).all()


def test_fit_cc50_recovers_a_curve_and_extrapolates_flat_data_past_the_top_dose():
    x = 100.0 / 2.0 ** np.arange(8)
    fit = fit_cc50(np.repeat(x, 2), np.repeat(survival(x, 12.0, 1.8), 2))
    assert abs(fit["logcc50"] - np.log10(12.0)) <= 0.02 and abs(fit["hill"] - 1.8) <= 0.1 and not fit["at_edge"]
    rng = np.random.default_rng(0)
    flat = fit_cc50(x, 1.0 + rng.normal(0, 0.05, 8))
    assert flat["logcc50"] > np.log10(100.0)
    partial = fit_cc50(x, survival(x, 150.0, 1.5))  # 65% left at the top dose: no measured CC50
    assert partial["logcc50"] > np.log10(100.0)
    assert fit_cc50([1.0, 10.0], [1.0, 0.5]) is None and fit_cc50([5.0, 5.0, 5.0], [1.0, 0.9, 0.8]) is None
    assert fit_cc50([1.0, 10.0, float("nan"), 0.0], [1.0, 0.5, 0.2, 0.3]) is None


def plate_rows(curves, vehicle=20000.0, kill=100.0, cloudy=(), cell_free=None):
    rows = [{"plate": 1, "role": "vehicle", "id": None, "conc_uM": None, "lum": vehicle + d, "turbidity": 500.0 + t}
            for d, t in ((-400.0, -20.0), (0.0, 0.0), (400.0, 20.0))]
    rows += [{"plate": 1, "role": "kill", "id": None, "conc_uM": None, "lum": kill, "turbidity": 500.0 + t}
             for t in (-10.0, 0.0, 10.0)]
    for cid, (doses, viab) in curves.items():
        for d, v in zip(doses, viab):
            rows.append({"plate": 1, "role": "cpd", "id": cid, "conc_uM": d, "lum": kill + v * (vehicle - kill),
                         "turbidity": 2000.0 if (cid, d) in cloudy else 500.0})
    if cell_free:
        rows += [{"plate": 1, "role": "cell_free_ref", "id": None, "conc_uM": None, "lum": 10100.0, "turbidity": 500.0}]
        for cid, fractions in cell_free.items():
            for d, f in zip(curves[cid][0], fractions):
                rows.append({"plate": 1, "role": "cell_free", "id": cid, "conc_uM": d, "lum": kill + f * 10000.0,
                             "turbidity": 500.0})
    return rows


def test_estimate_cc50_on_hand_made_plates():
    doses = [100.0, 50.0, 25.0, 12.5, 6.25, 3.125]
    flat = [1.0] * 6
    toxic = list(survival(np.array(doses), 10.0, 2.0))
    rows = plate_rows({"toxic": (doses, toxic), "clean": (doses, flat)})
    est = estimate_cc50(rows)
    assert not est["toxic"]["at_least"] and np.isclose(est["toxic"]["cc50_uM"], 10.0, rtol=0.05)
    assert est["clean"] == {"cc50_uM": 100.0, "at_least": True}
    # precipitate at the two top doses: they are dropped, the bound moves down, and a later cloudy dose is ignored
    rows = plate_rows({"clean": (doses, flat)}, cloudy={("clean", 100.0), ("clean", 50.0), ("clean", 6.25)})
    assert estimate_cc50(rows)["clean"] == {"cc50_uM": 25.0, "at_least": True}
    assert estimate_cc50(rows, use_turbidity=False)["clean"] == {"cc50_uM": 100.0, "at_least": True}
    # a compound that only hits the readout: corrected dose by dose; unreadable doses are dropped
    inhibition = [1 / (1 + d / 10.0) for d in doses]
    rows = plate_rows({"clean": (doses, inhibition)}, cell_free={"clean": inhibition})
    assert estimate_cc50(rows)["clean"] == {"cc50_uM": 100.0, "at_least": True}
    strong = [1 / (1 + d / 4.0) for d in doses]
    assert strong[0] < analysis.MIN_READOUT_FRACTION < strong[1] and strong[-1] < analysis.INTERFERENCE_FLAG
    rows = plate_rows({"clean": (doses, strong)}, cell_free={"clean": strong})
    assert estimate_cc50(rows)["clean"] == {"cc50_uM": 50.0, "at_least": True}
    # a control within the flag is left alone, so its noise is not added
    rows = plate_rows({"clean": (doses, flat)}, cell_free={"clean": [0.9, 1.1, 1.0, 1.0, 1.0, 1.0]})
    assert estimate_cc50(rows)["clean"] == {"cc50_uM": 100.0, "at_least": True}
    # no vehicle wells: nothing to normalise to, so no answer; no full-kill wells: zero signal is the zero
    assert estimate_cc50([r for r in plate_rows({"toxic": (doses, toxic)}) if r["role"] != "vehicle"]) == {}
    no_kill = estimate_cc50([r for r in plate_rows({"toxic": (doses, toxic)}) if r["role"] != "kill"])
    assert np.isclose(no_kill["toxic"]["cc50_uM"], 10.0, rtol=0.05)
    # everything cloudy, or toxic at the only dose left: no answer
    all_cloudy = plate_rows({"clean": (doses, flat)}, cloudy={("clean", d) for d in doses})
    assert estimate_cc50(all_cloudy) == {}
    one_left = plate_rows({"toxic": (doses, [0.1] * 6)}, cloudy={("toxic", d) for d in doses[:5]})
    assert estimate_cc50(one_left) == {}
    assert analysis.TURBIDITY_SD_MULT == 4.0
    assert viability_zprime(plate_rows({"clean": (doses, flat)}))[1] > 0.9


def test_spec_documents_every_design_field_and_is_plain_json():
    spec = make(NAME, seed=0).spec()
    assert set(spec["actions"]["run"]["design"]) == set(cw.DESIGN_SPEC) == set(CX_EXPERT) == set(CX_CHOICES)
    assert set(spec["actions"]["report"]["results"]["<compound id>"]) == {"cc50_uM", "at_least"}
    assert spec["labware"]["top_limit_uM"] == 100.0 and spec["labware"]["doubling_time_h"] == 24.0
    assert spec["labware"]["matched_dmso_pct"] == 1.0 and spec["actions"]["run"]["design"]["top_conc_uM"]["max"] == 100.0
    assert "untreated" not in spec["observation"]["wells"] and "untreated" not in spec["labware"]["layout"]
    assert "fails the QC gate" in spec["actions"]["run"]["design"]["match_dmso"]["help"]
    assert "a solvent control has to carry the solvent of the dosed wells" in spec["task"]
    assert json.loads(json.dumps(spec)) == spec
    obs = make(NAME, seed=0).reset()
    assert [c["id"] for c in obs["compounds"]] == [f"cpd{i}" for i in range(1, 7)]
    assert all(0.001 <= c["ec50_uM"] <= 1.0 for c in obs["compounds"])
