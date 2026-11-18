import ast
import collections
import copy
import json
from pathlib import Path

import numpy as np
import pytest

from drugdevgym import envqc, knowledge, make
from drugdevgym.analysis import (ACTIVE_FRACTION, activity_fraction, estimate_potency, fit_4pl, fit_4pl_capped,
                                 selectivity_window)
from drugdevgym.chem import Kit, hill
from drugdevgym.envs import selectivity_panel as sp
from drugdevgym.envs.base import COST_WEIGHT
from drugdevgym.envs.dose_response import STD_SERIES_NM, ZPRIME_GATE, plate_cost
from drugdevgym.policies import (POLICIES, SP_INTERACTIONS, SP_OFF, SP_ON, SP_OVERCLAIM_LOG, SP_RANDOM_WINDOW,
                                 SelectivityPolicy, SelectivityRandom, run_episode)

NAME = "selectivity_panel"
REFERENCES = ("expert", "random", "null")
TRAPS = [p for p in POLICIES[NAME] if p not in REFERENCES]
SINGLE_TRAPS = [p for p in TRAPS if p not in SP_INTERACTIONS]  # the expert with exactly one change each
HIDDEN_WORDS = ("truth", "pec50", "emax", "hill", "sticky", "margin", "liab", "pharmacology", "cell_line",
                "achievable", "supported", "basal")
SCREEN = {**SP_OFF, "n_points": 1, "standards": False, "controls_per_side": 3}  # the cheapest plate that passes QC


def play(policy, seed, rig="nominal"):
    return run_episode(make(NAME, seed=seed, rig=rig), POLICIES[NAME][policy](seed))


def trajectory(policy_name, seed, rig="nominal"):
    """Every action, observation, reward and info of one episode, plus the environment."""
    env, policy = make(NAME, seed=seed, rig=rig), POLICIES[NAME][policy_name](seed)
    obs, steps, done = env.reset(), [], False
    while not done:
        action = policy.act(obs)
        sent = copy.deepcopy(action)
        obs, reward, done, info = env.step(action)
        steps.append((sent, copy.deepcopy(obs), reward, done, copy.deepcopy(info)))
    return env, steps


def truth_of(env, cid):
    """(true window W, achievable window A, on-target pEC50) from the environment's hidden state."""
    at = env.pharmacology[cid]
    on = at["on"].pec50
    return on - max(at[t].pec50 for t in sp.OFF_TARGETS), on - sp.P_FLOOR, on


def run_all(env, off_design=SP_OFF, on_design=SP_ON):
    """One plate per target. Returns the last observation."""
    obs = None
    for target in sp.TARGETS:
        obs, reward, done, info = env.step({"type": "run_plate", "target": target,
                                            "design": on_design if target == "on" else off_design})
        assert "error" not in obs and reward == 0.0 and not done and info == {}
    return obs


def perfect_report(env):
    """The best possible report, written from the hidden truth: exact where measurable, the bound otherwise."""
    windows = {}
    for cid in env.pharmacology:
        w, a, _ = truth_of(env, cid)
        windows[cid] = {"window": min(w, a), "at_least": bool(w > a)}
    return {"type": "report", "windows": windows}


# ---------------------------------------------------------------- determinism and hidden truth

def test_same_seed_gives_the_same_episode():
    for policy in ("expert", "random", "single_point"):
        a, b = trajectory(policy, 3)[1], trajectory(policy, 3)[1]
        assert json.dumps(a, sort_keys=True) == json.dumps(b, sort_keys=True)
    assert play("expert", 3)[1]["truth"] != play("expert", 4)[1]["truth"]


def test_the_truth_depends_on_the_seed_and_not_on_the_policy_or_the_rig():
    for seed in (7, 8):
        truths = [play(policy, seed)[1]["truth"] for policy in POLICIES[NAME]]
        assert all(t == truths[0] for t in truths)
        assert play("expert", seed, "harsh")[1]["truth"] == truths[0]
    assert play("null", 7)[1]["truth"] != play("null", 8)[1]["truth"]


def test_hidden_draws_are_not_shared_with_dose_response_on_the_same_seed():
    for seed in range(5):
        assert make(NAME, seed=seed).cell_lines["on"].kit.k != make("dose_response", seed=seed).assay.kit.k


def test_the_hidden_truth_has_the_documented_shape():
    on, margins, liable, inactive = [], [], 0, []
    for seed in range(300):
        for at in make(NAME, seed=seed).pharmacology.values():
            p = at["on"].pec50
            m = sorted(p - at[t].pec50 for t in sp.OFF_TARGETS)
            on.append(p)
            margins.append(m)
            liable += m[0] < sp.LIABILITY_MARGIN_RANGE[1]
            inactive += [at[t].pec50 < sp.P_FLOOR for t in sp.OFF_TARGETS]
            assert len({at[t].sticky for t in sp.TARGETS}) == 1  # stickiness belongs to the compound
            assert sp.LIABILITY_MARGIN_RANGE[0] <= m[0] <= sp.MARGIN_RANGE[1]
            assert all(sp.MARGIN_RANGE[0] <= x <= sp.MARGIN_RANGE[1] for x in m[1:])  # at most one close off-target
    assert sp.ON_PEC50_RANGE[0] <= min(on) and max(on) <= sp.ON_PEC50_RANGE[1]
    assert abs(liable / len(on) - sp.LIABILITY_RATE) < 0.05  # a minority carries one close liability
    assert 0.65 < np.mean(inactive) < 0.9  # most off-target pairs show nothing at the labware limit
    assert np.isclose(sp.P_FLOOR, 4.0) and sp.MAX_TOP_NM == 100000.0


@pytest.mark.parametrize("policy_name", list(POLICIES[NAME]))
@pytest.mark.parametrize("rig", ["nominal", "harsh"])
def test_no_observation_names_a_truth_field(policy_name, rig):
    _, steps = trajectory(policy_name, 5, rig)
    start = make(NAME, seed=5, rig=rig).reset()
    seen = [start] + [obs for _, obs, _, _, _ in steps]
    for obs in seen:
        text = json.dumps({k: v for k, v in obs.items() if k != "spec"}).lower()  # the spec is the public task text
        assert not any(word in text for word in HIDDEN_WORDS), (policy_name, text[:200])
    assert all(info == {} for _, _, _, done, info in steps if not done)
    final = steps[-1]
    assert final[3] and "truth" in final[4] and final[2] <= 1.0 and final[1] == {"phase": "done"}
    assert set(final[4]) == {"quality", "qc_pass", "zprime", "plates_used", "cost", "detail", "truth"}
    assert all(set(d) == {"claim", "score", "on_target_tested", "tested_limit"} for d in final[4]["detail"].values())


def test_the_spec_is_json_and_states_the_actions_and_the_reward():
    spec = make(NAME, seed=0).spec()
    json.dumps(spec)
    assert spec["env"] == NAME and set(spec["actions"]) == {"run_plate", "report"}
    assert set(spec["actions"]["run_plate"]["design"]) == set(SP_OFF)
    assert spec["actions"]["run_plate"]["design"]["n_points"]["min"] == 1  # a single-point screen is allowed
    assert spec["labware"]["plate_budget"] == sp.PLATE_BUDGET == 6
    text = json.dumps(spec)
    assert not any(trap in text for trap in TRAPS)  # the task text never names a reference mistake


# ---------------------------------------------------------------- malformed input

NAN, INF, HUGE = float("nan"), float("inf"), 10 ** 400
BAD_NUMBERS = (NAN, INF, -INF, HUGE, -HUGE, None, "3", [3], [[3]], {"a": 1}, True)
BAD_PLATES = [
    {"type": "run_plate"},
    {"type": "run_plate", "target": "on"},
    {"type": "run_plate", "design": SP_ON},
    {"type": "run_plate", "target": "off9", "design": SP_OFF},
    {"type": "run_plate", "target": None, "design": SP_OFF},
    {"type": "run_plate", "target": ["on"], "design": SP_OFF},
    {"type": "run_plate", "target": 1, "design": SP_OFF},
    {"type": "run_plate", "target": HUGE, "design": SP_OFF},
    {"type": "run_plate", "target": "on", "design": None},
    {"type": "run_plate", "target": "on", "design": [SP_ON]},
    {"type": "run_plate", "target": "on", "design": "default"},
    {"type": "run_plate", "target": "on", "design": {}},
    {"type": "run_plate", "target": "on", "design": {**SP_ON, "extra": 1}},
    {"type": "run_plate", "target": "on", "design": {1: 2, "b": 3}},
    {"type": "run_plate", "target": "on", "design": {k: v for k, v in SP_ON.items() if k != "lid"}},
    {"type": "run_plate", "target": "on", "design": {**SP_ON, "mix": 1}},
    {"type": "run_plate", "target": "on", "design": {**SP_ON, "lid": None}},
    {"type": "run_plate", "target": "on", "design": {**SP_ON, "n_points": 2.5}},
    {"type": "run_plate", "target": "on", "design": {**SP_ON, "n_points": 0}},
    {"type": "run_plate", "target": "on", "design": {**SP_ON, "n_points": 13}},
    {"type": "run_plate", "target": "on", "design": {**SP_ON, "top_conc_nM": 200000.0}},
    {"type": "run_plate", "target": "on", "design": {**SP_ON, "replicates": 4}},  # does not fit the plate
    {"type": "run_plate", "target": "on", "design": {**SP_ON, "transfer_ul": 100.0, "dilution_factor": 10.0}},
    {"type": "run_plate", "target": "on", "design": SP_ON, "plate": 2},
    {"type": "run_plate", "target": "on", "design": SP_ON, 7: 7},
    *[{"type": "run_plate", "target": "on", "design": {**SP_ON, field: bad}}
      for field in ("top_conc_nM", "dilution_factor", "n_points", "replicates", "transfer_ul", "controls_per_side")
      for bad in BAD_NUMBERS],
    *[{"type": "run_plate", "target": "on", "design": SP_ON, "compounds": bad}
      for bad in (None, "cpd1", [], ["cpd9"], ["cpd1", "cpd1"], [["cpd1"]], [1], [None], [NAN], {"cpd1": 1}, 4, HUGE)],
]
GOOD = {"window": 3.0, "at_least": True}
BAD_REPORTS = [
    {"type": "report"},
    {"type": "report", "windows": None},
    {"type": "report", "windows": [GOOD]},
    {"type": "report", "windows": "cpd1"},
    {"type": "report", "windows": 3.0},
    {"type": "report", "windows": {}, "plates_used": [1]},
    {"type": "report", "windows": {"cpd9": GOOD}},
    {"type": "report", "windows": {1: GOOD}},
    {"type": "report", "windows": {"cpd1": 3.0}},
    {"type": "report", "windows": {"cpd1": [3.0, True]}},
    {"type": "report", "windows": {"cpd1": "at least 3"}},
    {"type": "report", "windows": {"cpd1": {}}},
    {"type": "report", "windows": {"cpd1": {"window": 3.0}}},
    {"type": "report", "windows": {"cpd1": {"at_least": True}}},
    {"type": "report", "windows": {"cpd1": {**GOOD, "pEC50": 8.0}}},
    {"type": "report", "windows": {"cpd1": {"window": 3.0, "at_least": "yes"}}},
    {"type": "report", "windows": {"cpd1": {"window": 3.0, "at_least": 1}}},
    {"type": "report", "windows": {"cpd1": {"window": 3.0, "at_least": None}}},
    {"type": "report", "windows": {"cpd1": {"window": 12.5, "at_least": True}}},
    *[{"type": "report", "windows": {"cpd1": {"window": bad, "at_least": True}}} for bad in BAD_NUMBERS],
]


@pytest.mark.parametrize("action", BAD_PLATES, ids=range(len(BAD_PLATES)))
def test_a_malformed_plate_is_rejected_without_an_exception_or_a_cost(action):
    env = make(NAME, seed=1)
    obs, reward, done, info = env.step(action)
    assert "error" in obs and reward == 0.0 and not done and info == {}
    assert env.cost == 0.0 and env.plates == [] and env.step_log == []


@pytest.mark.parametrize("action", BAD_REPORTS, ids=range(len(BAD_REPORTS)))
def test_a_malformed_report_is_rejected_and_the_episode_stays_open(action):
    env = make(NAME, seed=1)
    obs, reward, done, info = env.step(action)
    assert "error" in obs and reward == 0.0 and not done and info == {}
    _, _, done, info = env.step({"type": "report", "windows": {"cpd1": None, "cpd2": {"window": 3, "at_least": False}}})
    assert done and "truth" in info  # null is "no answer", and a whole number is a number


def test_a_malformed_claim_is_named_as_a_claim_not_as_a_plate_design():
    for entry in (3.0, [1, 2], "at least 3", True):
        obs, *_ = make(NAME, seed=1).step({"type": "report", "windows": {"cpd1": entry}})
        assert obs["error"] == "cpd1: must be {window, at_least} or null"
    obs, *_ = make(NAME, seed=1).step({"type": "report", "windows": {"cpd1": {"window": 3.0}}})
    assert obs["error"] == "cpd1: missing field: at_least"


def test_five_invalid_actions_end_the_episode():
    env = make(NAME, seed=1)
    run_all(env)
    for _ in range(4):
        obs, reward, done, _ = env.step({"type": "report", "windows": None})
        assert "error" in obs and reward == 0.0 and not done
    obs, reward, done, info = env.step({"type": "dance"})
    assert done and info == {"truncated": True} and np.isclose(reward, -COST_WEIGHT * env.cost / sp.COST_REFERENCE)
    with pytest.raises(RuntimeError):
        env.step({"type": "report", "windows": {}})


def test_the_plate_budget_is_enforced():
    env = make(NAME, seed=2)
    for i in range(sp.PLATE_BUDGET):
        obs, *_ = env.step({"type": "run_plate", "target": "off1", "design": SCREEN})
        assert obs["plate"] == i + 1 and obs["plates_left"] == sp.PLATE_BUDGET - i - 1
    obs, reward, done, _ = env.step({"type": "run_plate", "target": "off1", "design": SCREEN})
    assert "no plates left" in obs["error"] and not done and len(env.plates) == sp.PLATE_BUDGET


# ---------------------------------------------------------------- what a plate is

def test_a_plate_carries_the_chosen_compounds_on_one_target():
    env = make(NAME, seed=4)
    obs, *_ = env.step({"type": "run_plate", "target": "off2", "compounds": ["cpd3", "cpd1"],
                        "design": {**SP_OFF, "replicates": 2}})
    assert obs["target"] == "off2" and obs["compounds"] == ["cpd1", "cpd3"]
    cpd = [w for w in obs["wells"] if w["role"] == "cpd"]
    assert {w["id"] for w in cpd} == {"cpd1", "cpd3"} and len(cpd) == 2 * 10 * 2
    assert set(obs["wells"][0]) == {"well", "role", "id", "conc_nM", "std_nM", "ratio"}
    assert np.isclose(obs["cost"], round(plate_cost(len(obs["wells"]), 2 * (9 + 20)), 2))
    obs, *_ = env.step({"type": "run_plate", "target": "on", "design": SCREEN})
    cpd = [w for w in obs["wells"] if w["role"] == "cpd"]
    assert len(cpd) == 4 and {w["conc_nM"] for w in cpd} == {SCREEN["top_conc_nM"]}  # a single-point screen


def test_every_plate_goes_through_the_rig_and_every_fill_is_labelled():
    env, steps = trajectory("expert", 0)
    plates = len(env.plates)
    assert plates >= 5 and {s["run"] for s in env.step_log} == set(range(1, plates + 1))
    assert {s["op"] for s in env.step_log} == {"fill", "transfer", "incubate", "read"}
    labels = {s.get("what") for s in env.step_log if s["op"] == "fill"}
    assert None not in labels and "buffer" not in labels
    assert {"cells", "diluent", "vehicle", "reference agonist", "detection mix"} <= labels
    assert any(str(x).startswith("analyte standard") for x in labels) and any(str(x).startswith("cpd1 at") for x in labels)


def test_a_starved_ladder_reads_as_inactive_through_the_rig():
    # dead volume acts on this environment too: an off-target curve the compound never reached is flat
    seed = next(s for s in range(200) if truth_of(make(NAME, seed=s), "cpd1")[0] < 2.0)
    env = make(NAME, seed=seed)
    at = env.pharmacology["cpd1"]
    target = max(sp.OFF_TARGETS, key=lambda t: at[t].pec50)  # the close liability
    good, *_ = env.step({"type": "run_plate", "target": target, "design": SP_OFF})
    dry, *_ = env.step({"type": "run_plate", "target": target, "design": {**SP_OFF, "transfer_ul": 5.0}})
    assert "pEC50" in estimate_potency(good["wells"])["cpd1"]
    assert "below" in estimate_potency(dry["wells"])["cpd1"]


# ---------------------------------------------------------------- observations are copies

def test_editing_an_observation_cannot_change_the_score():
    def run(tamper: bool):
        env, policy = make(NAME, seed=6), POLICIES[NAME]["expert"](6)
        obs, done = env.reset(), False
        while not done:
            action = policy.act(copy.deepcopy(obs))  # the policy keeps its own copy of what it saw
            sent = copy.deepcopy(action)
            if tamper:  # wreck every object the environment handed out
                obs.get("spec", {}).clear()
                for key in ("compounds", "targets", "wells"):
                    for item in obs.get(key, []):
                        if isinstance(item, dict):
                            item["role"], item["ratio"], item["conc_nM"] = "pos", 0.0, 1.0
                    obs.get(key, []).clear()
                obs.clear()
            obs, reward, done, info = env.step(sent)
            if tamper:  # and every object it was given
                for part in (sent.get("design"), sent.get("compounds"), sent.get("windows")):
                    if part is not None:
                        part.clear()
                sent.clear()
        return reward, info

    assert run(False) == run(True)


def test_faking_controls_in_an_observation_cannot_open_the_gate():
    env = make(NAME, seed=6)
    obs = run_all(env, off_design={**SP_OFF, "controls_per_side": 0})
    for w in obs["wells"][:12]:
        w["role"] = "neg" if w["role"] != "pos" else "pos"
    reward, info = env.step(perfect_report(env))[1::2]
    assert not info["qc_pass"] and np.isclose(info["quality"], 1.0)
    assert np.isclose(reward, 0.5 - COST_WEIGHT * info["cost"] / sp.COST_REFERENCE)
    assert make(NAME, seed=6).reset()["spec"] == make(NAME, seed=6).spec()


# ---------------------------------------------------------------- the reward is what spec() says

def test_the_reward_text_carries_the_numbers_the_code_uses():
    text = make(NAME, seed=0).spec()["reward"]
    assert (sp.TOLERANCE, sp.SUPPORT_MARGIN, sp.FALSE_CLAIM_MARGIN, sp.CLAIM_PENALTY, sp.P_FLOOR) == (0.5, 0.5, 1.0, 0.5, 4.0)
    assert (COST_WEIGHT, sp.COST_REFERENCE, ZPRIME_GATE, sp.MIN_CONTROLS) == (0.15, 600.0, 0.5, 3)
    for phrase in ("A = true on-target pEC50 minus 4,", "nothing seen up to 100000 nM",
                   "a plate passes when it has at least 3 controls a side and control Z' >= 0.5 on the raw ratio",
                   "rests on one plate: among the plates that tested the pair, those that pass plate QC if any do "
                   "(otherwise all of them), the one with the highest top_conc_nM (the earliest if tied)",
                   "the highest, over the 4 off-targets, of 9 - log10(top_conc_nM of the plate the pair rests on)",
                   "if the compound was never tested on the on-target or on any one of the off-targets, no window is "
                   "supported",
                   "max(0, 1 - |x - W| / 0.5)", "max(0, 1 - max(|x - A|, x - W) / 0.5)",
                   "if x > S + 0.5 (unsupported by the tested range) or x > W + 1 (more selectivity than the compound "
                   "has, by over a log unit), the score is -0.5 instead",
                   "No answer scores 0", "Quality = mean over compounds; it can be negative",
                   "halved when it is positive and the QC gate fails", "minus 0.15 * cost / 600",
                   "at least one plate was run and, for every compound and target pair that was tested, the plate the "
                   "pair rests on passes plate QC"):
        assert phrase in text, phrase
    assert "max(0, mean" not in text and "latest plate" not in text  # the two rules the review replaced


def recompute(env, steps):
    """The reward rebuilt from spec()['reward'] alone: hidden truth, the actions sent and the wells observed."""
    plates = [(action, obs) for action, obs, *_ in steps if action["type"] == "run_plate"]
    passes = []
    for _, obs in plates:
        pos = np.array([w["ratio"] for w in obs["wells"] if w["role"] == "pos"])
        neg = np.array([w["ratio"] for w in obs["wells"] if w["role"] == "neg"])
        ok = len(pos) >= 3 and len(neg) >= 3
        passes.append(bool(ok and 1 - 3 * (pos.std(ddof=1) + neg.std(ddof=1)) / abs(pos.mean() - neg.mean()) >= 0.5))
    rests = {}  # (compound, target) -> index of the plate the pair rests on
    for i, (action, _) in enumerate(plates):
        for cid in action.get("compounds", list(env.pharmacology)):
            tested = rests.setdefault((cid, action["target"]), [])
            tested.append(i)
    for pair, tested in rests.items():
        pool = [i for i in tested if passes[i]] or tested
        top = max(plates[i][0]["design"]["top_conc_nM"] for i in pool)
        rests[pair] = min(i for i in pool if plates[i][0]["design"]["top_conc_nM"] == top)
    gate = bool(rests) and all(passes[i] for i in rests.values())
    claims = steps[-1][0].get("windows", {})
    scores = []
    for cid in env.pharmacology:
        w, a, on = truth_of(env, cid)
        claim = claims.get(cid)
        if claim is None:
            scores.append(0.0)
            continue
        x = claim["window"]
        if all((cid, t) in rests for t in sp.TARGETS):
            s = on - max(9 - np.log10(plates[rests[(cid, t)]][0]["design"]["top_conc_nM"]) for t in sp.OFF_TARGETS)
        else:
            s = -np.inf
        if x > s + 0.5 or x > w + 1:
            scores.append(-0.5)
        elif claim["at_least"]:
            scores.append(max(0.0, 1 - max(abs(x - a), x - w) / 0.5))
        else:
            scores.append(max(0.0, 1 - abs(x - w) / 0.5))
    quality = float(np.mean(scores))
    cost = plates[-1][1]["cost"] if plates else 0.0
    halved = quality > 0 and not gate
    return quality * (0.5 if halved else 1.0) - 0.15 * cost / 600.0, quality, gate, scores


@pytest.mark.parametrize("policy_name", list(POLICIES[NAME]))
def test_the_reward_and_gate_arithmetic_match_the_spec_text(policy_name):
    for seed, rig in ((0, "nominal"), (1, "nominal"), (2, "nominal"), (3, "harsh")):
        env, steps = trajectory(policy_name, seed, rig)
        reward, quality, gate, scores = recompute(env, steps)
        info = steps[-1][4]
        assert sum(r for _, _, r, _, _ in steps[:-1]) == 0.0
        assert np.isclose(steps[-1][2], reward, atol=1e-4), (policy_name, seed)  # the observed cost is rounded
        assert np.isclose(info["quality"], quality) and info["qc_pass"] == gate
        assert np.allclose([d["score"] for d in info["detail"].values()], scores)


def test_the_spec_text_reproduces_random_episodes_with_repeats_subsets_and_mixed_claims():
    # the registered policies never test a pair twice at different doses; random scripts do, with and without controls
    from drugdevgym.policies import random_sp_design
    arms = collections.Counter()
    for seed in range(120):
        rng = np.random.default_rng([seed, 7])
        env = make(NAME, seed=seed, rig="harsh" if seed % 3 == 0 else "nominal")
        steps, plan = [], []
        if seed % 4:  # most scripts cover the panel, in a random order, on good or random plates; then one more plate
            for i in rng.permutation(len(sp.TARGETS)):
                design = random_sp_design(rng) if rng.random() < 0.3 else (SP_ON if sp.TARGETS[i] == "on" else SP_OFF)
                plan.append((sp.TARGETS[i], list(env.pharmacology), design))
        for _ in range(int(rng.integers(0, sp.PLATE_BUDGET + 1 - len(plan)))):  # anything goes
            chosen = [c for c in env.pharmacology if rng.random() < 0.7] or ["cpd1"]
            plan.append((sp.TARGETS[int(rng.integers(len(sp.TARGETS)))], chosen, random_sp_design(rng, len(chosen))))
        for target, chosen, design in plan:
            action = {"type": "run_plate", "target": target, "compounds": chosen, "design": dict(design)}
            obs, reward, done, info = env.step(copy.deepcopy(action))
            assert "error" not in obs
            steps.append((action, obs, reward, done, info))
        windows = {}
        for cid in env.pharmacology:
            w, a, _ = truth_of(env, cid)
            kind = int(rng.integers(4))
            if kind == 0:
                continue  # no answer
            x = (w, a, a - 1.5, float(rng.uniform(-1.0, 6.0)))[int(rng.integers(4))] + float(rng.normal(0.0, 0.2))
            windows[cid] = {"window": float(np.clip(x, *sp.WINDOW_RANGE)), "at_least": bool(kind % 2)}
        action = {"type": "report", "windows": windows}
        obs, reward, done, info = env.step(copy.deepcopy(action))
        steps.append((action, obs, reward, done, info))
        again, quality, gate, scores = recompute(env, steps)
        assert np.isclose(reward, again, atol=1e-4), seed
        assert np.isclose(info["quality"], quality) and info["qc_pass"] == gate, seed
        assert np.allclose([d["score"] for d in info["detail"].values()], scores), seed
        arms.update("penalised" if x < 0 else "credit" if x > 0 else "zero" for x in scores)
        arms.update(["gate passes" if gate else "gate fails", "negative quality" if quality < 0 else "quality >= 0"])
    assert all(arms[k] >= 10 for k in ("penalised", "credit", "zero", "gate passes", "gate fails", "negative quality",
                                       "quality >= 0")), arms


def test_every_arm_of_the_claim_score():
    w, a, s = 3.0, 4.0, 4.0  # a measurable off-target, everything tested to the labware limit
    assert sp.score_claim(3.0, False, w, a, s) == 1.0
    assert np.isclose(sp.score_claim(3.25, False, w, a, s), 0.5) and sp.score_claim(3.6, False, w, a, s) == 0.0
    assert sp.score_claim(3.0, True, w, a, s) == 0.0  # hedging a measured window as a bound, a log from the limit
    w = 5.2  # nothing measurable: the right answer is the bound at the labware limit
    assert sp.score_claim(4.0, True, w, a, s) == 1.0
    assert np.isclose(sp.score_claim(3.75, True, w, a, s), 0.5)  # true but looser than the labware allows
    assert sp.score_claim(3.0, True, w, a, s) == 0.0 and sp.score_claim(-5.0, True, w, a, s) == 0.0
    assert np.isclose(sp.score_claim(4.25, True, w, a, s), 0.5)  # past the limit, inside the margin
    assert sp.score_claim(4.51, True, w, a, s) == -0.5 and sp.score_claim(5.2, False, w, a, s) == -0.5  # unsupported
    assert sp.score_claim(4.0, False, w, a, s) == 0.0  # the bound written as a number
    assert np.isclose(sp.score_claim(3.9, True, 3.8, a, s), 0.8)  # a bound that is false by 0.1
    s = 2.5  # tested only to a low top dose: the tested range supports a narrow window
    assert sp.score_claim(2.5, True, w, a, s) == 0.0 and sp.score_claim(4.0, True, w, a, s) == -0.5
    assert sp.score_claim(-5.0, True, w, a, float("-inf")) == -0.5  # an untested off-target supports nothing
    w, s = 2.0, 4.0  # a close liability, tested to the labware limit: the range supports the claim, the truth does not
    assert sp.score_claim(4.0, True, w, a, s) == -0.5 and sp.score_claim(3.01, False, w, a, s) == -0.5  # false by over a log
    assert sp.score_claim(3.0, False, w, a, s) == 0.0 and sp.score_claim(2.9, True, w, a, s) == 0.0  # wrong, not penalised
    assert sp.score_claim(0.5, False, w, a, s) == 0.0  # too little selectivity claimed is wrong, never penalised
    assert sp.supported_window(8.0, [100000.0, 100000.0, None, 100000.0]) == float("-inf")
    assert sp.supported_window(8.0, [100000.0] * 4, on_tested=False) == float("-inf")  # nor does an untested on-target
    assert np.isclose(sp.supported_window(8.0, [100000.0, 3000.0, 100000.0, 100000.0]), 8.0 - (9 - np.log10(3000.0)))
    assert np.isclose(sp.tested_limit(100000.0), 4.0) and np.isclose(sp.true_window(8.0, [3.0, 6.5, 1.0, 4.0]), 1.5)


def test_a_perfect_report_scores_one_minus_cost_and_claims_are_scored_per_compound():
    env = make(NAME, seed=11)
    run_all(env)
    report = perfect_report(env)
    kinds = {cid: c["at_least"] for cid, c in report["windows"].items()}
    _, reward, done, info = env.step(report)
    assert done and info["qc_pass"] and np.isclose(info["quality"], 1.0)
    assert np.isclose(reward, 1.0 - COST_WEIGHT * info["cost"] / sp.COST_REFERENCE) and reward < 1.0
    assert all(d["score"] == 1.0 for d in info["detail"].values())
    assert {cid: d["claim"]["at_least"] for cid, d in info["detail"].items()} == kinds


def test_an_unsupported_claim_is_penalised_and_a_wrong_one_is_not():
    def scores(seed, edit, plates=run_all):
        env = make(NAME, seed=seed)
        plates(env)
        report = perfect_report(env)
        edit(env, report["windows"])
        return env.step(report)[3]

    def far_past_the_limit(seed):
        w, a, _ = truth_of(make(NAME, seed=seed), "cpd1")
        return w > a + 1.2

    clean = next(s for s in range(100) if far_past_the_limit(s))

    def past_the_range(env, windows):
        windows["cpd1"]["window"] += 1.0  # true, but a log beyond the top dose tested

    info = scores(clean, past_the_range)
    assert info["detail"]["cpd1"]["score"] == -0.5
    assert np.isclose(info["quality"], (3 * 1.0 - 0.5) / 4)

    def as_a_number(env, windows):
        windows["cpd1"]["at_least"] = False  # the bound written as an exact window

    assert scores(clean, as_a_number)["detail"]["cpd1"]["score"] == 0.0

    def off1_at_low_dose(env):
        for target in sp.TARGETS:
            env.step({"type": "run_plate", "target": target,
                      "design": {**SP_OFF, "top_conc_nM": 1000.0} if target == "off1" else SP_OFF})

    info = scores(clean, lambda env, windows: None, off1_at_low_dose)  # one off-target tested only to 1 uM
    assert info["detail"]["cpd1"]["score"] == -0.5 and np.isclose(info["detail"]["cpd1"]["tested_limit"], 6.0)


def test_a_false_claim_of_selectivity_inside_the_tested_range_is_penalised():
    # a close liability sits on the agent's own off-target curve; calling the compound selective anyway costs -0.5
    seed, cid = next((s, c) for s in range(100) for c in ("cpd1", "cpd2", "cpd3", "cpd4")
                     if truth_of(make(NAME, seed=s), c)[0] < 2.0)
    for at_least in (True, False):
        env = make(NAME, seed=seed)
        run_all(env)
        report = perfect_report(env)
        w, a, _ = truth_of(env, cid)
        report["windows"][cid] = {"window": a, "at_least": at_least}  # the labware limit, claimed without reading the curve
        info = env.step(report)[3]
        assert a > w + sp.FALSE_CLAIM_MARGIN and np.isclose(info["detail"][cid]["tested_limit"], sp.P_FLOOR)
        assert info["detail"][cid]["score"] == -0.5 and np.isclose(info["quality"], (3 * 1.0 - 0.5) / 4)
    env = make(NAME, seed=seed)
    run_all(env)
    report = perfect_report(env)
    w, a, _ = truth_of(env, cid)
    report["windows"][cid] = {"window": w - 1.5, "at_least": False}  # far too little selectivity claimed: wrong, not penalised
    assert env.step(report)[3]["detail"][cid]["score"] == 0.0


def test_a_report_of_penalised_claims_scores_below_no_report_and_a_failed_gate_never_halves_a_penalty():
    def run(windows, plates=()):
        env = make(NAME, seed=14)
        for target, design in plates:
            env.step({"type": "run_plate", "target": target, "design": design})
        _, reward, _, info = env.step({"type": "report", "windows": windows(env)})
        return reward, info

    def claim_all(x):
        return lambda env: {cid: {"window": x, "at_least": True} for cid in env.pharmacology}

    reward, info = run(lambda env: {})  # nothing run, nothing claimed
    assert reward == 0.0 and info["quality"] == 0.0 and not info["qc_pass"]
    for x in (2.0, 12.0):  # nothing run, selectivity claimed for every compound
        reward, info = run(claim_all(x))
        assert info["quality"] == -0.5 and not info["qc_pass"] and reward == -0.5  # the failed gate does not halve it
    five_on = [("on", SP_ON)] * 5  # the whole budget on the on-target, then an untested claim
    reward, info = run(claim_all(12.0), five_on)
    silent = run(lambda env: {}, five_on)[0]
    assert info["qc_pass"] and info["quality"] == -0.5 and np.isclose(reward, silent - 0.5)
    assert np.isclose(silent, -COST_WEIGHT * info["cost"] / sp.COST_REFERENCE)
    # a mixed report: the penalties are not floored away by the credit either
    env = make(NAME, seed=14)
    run_all(env)
    report = perfect_report(env)
    for cid in ("cpd1", "cpd2", "cpd3"):
        report["windows"][cid] = {"window": 12.0, "at_least": True}
    info = env.step(report)[3]
    assert np.isclose(info["quality"], (1.0 - 3 * 0.5) / 4) and info["quality"] < 0


def test_nothing_is_supported_until_every_off_target_is_tested():
    env = make(NAME, seed=12)
    for target in ("on", "off1", "off2", "off3"):
        env.step({"type": "run_plate", "target": target, "design": SP_OFF})
    windows = {cid: {"window": -5.0, "at_least": True} for cid in env.pharmacology}
    _, reward, _, info = env.step({"type": "report", "windows": windows})
    assert all(d["score"] == -0.5 and d["tested_limit"] is None and d["on_target_tested"] for d in info["detail"].values())
    assert info["quality"] == -0.5 and info["qc_pass"] and np.isclose(reward, -0.5 - COST_WEIGHT * info["cost"] / 600.0)


def test_nothing_is_supported_until_the_on_target_is_tested():
    # four token off-target screens that pass plate QC, no on-target plate, one constant answer
    for claim in ({"window": 3.5, "at_least": False}, {"window": 3.5, "at_least": True}, {"window": -5.0, "at_least": True}):
        env = make(NAME, seed=12)
        for target in sp.OFF_TARGETS:
            env.step({"type": "run_plate", "target": target, "design": SCREEN})
        _, reward, _, info = env.step({"type": "report", "windows": {cid: dict(claim) for cid in env.pharmacology}})
        assert info["qc_pass"] and info["quality"] == -0.5
        assert all(d["score"] == -0.5 and not d["on_target_tested"] and np.isclose(d["tested_limit"], sp.P_FLOOR)
                   for d in info["detail"].values())
        assert np.isclose(reward, -0.5 - COST_WEIGHT * info["cost"] / 600.0)
    # the on-target tested for one compound only supports that compound only
    env = make(NAME, seed=12)
    run_all(env, on_design=SCREEN)  # a screen is not a curve, but the pair was tested
    env2 = make(NAME, seed=12)
    env2.step({"type": "run_plate", "target": "on", "compounds": ["cpd2"], "design": SP_ON})
    for target in sp.OFF_TARGETS:
        env2.step({"type": "run_plate", "target": target, "design": SP_OFF})
    info = env2.step(perfect_report(env2))[3]
    assert [d["on_target_tested"] for d in info["detail"].values()] == [False, True, False, False]
    assert [d["score"] for d in info["detail"].values()] == [-0.5, 1.0, -0.5, -0.5]


def test_the_gate_and_the_support_rest_on_the_best_passing_plate_of_each_pair():
    def episode(extra, first=sp.TARGETS):
        env = make(NAME, seed=13)
        for target in first:
            weak = {**SP_OFF, "controls_per_side": 2} if target == "off3" else SP_OFF
            env.step({"type": "run_plate", "target": target, "design": weak})
        for action in extra:
            assert "error" not in env.step(action)[0]
        return env.step(perfect_report(env))[1::2]

    reward, info = episode([])
    assert not info["qc_pass"] and np.isclose(reward, 0.5 * 1.0 - COST_WEIGHT * info["cost"] / 600.0)
    assert info["plates_used"] == [1, 2, 3, 4, 5] and set(info["zprime"]) == {"1", "2", "3", "4", "5"}
    reward, info = episode([{"type": "run_plate", "target": "off3", "design": SP_OFF}])  # the repeat replaces plate 4
    assert info["qc_pass"] and info["plates_used"] == [1, 2, 3, 5, 6]
    assert np.isclose(reward, 1.0 - COST_WEIGHT * info["cost"] / 600.0)
    # a repeat that covers one compound leaves the failed plate in force for the other three
    reward, info = episode([{"type": "run_plate", "target": "off3", "compounds": ["cpd1"], "design": SP_OFF}])
    assert not info["qc_pass"] and info["plates_used"] == [1, 2, 3, 4, 5, 6]

    def after_the_panel(extra):
        env = make(NAME, seed=13)
        run_all(env)
        assert "error" not in env.step(extra)[0]
        return env.step(perfect_report(env))[1::2]

    alone = make(NAME, seed=13)
    run_all(alone)
    reward0, info0 = alone.step(perfect_report(alone))[1::2]
    assert info0["qc_pass"] and np.isclose(info0["quality"], 1.0)
    # a later plate at a lower top dose takes nothing away: the pair still rests on its 100 uM plate
    low = {"type": "run_plate", "target": "off1", "compounds": ["cpd2"], "design": {**SP_OFF, "top_conc_nM": 100.0}}
    reward, info = after_the_panel(low)
    assert all(np.isclose(d["tested_limit"], sp.P_FLOOR) and d["score"] == 1.0 for d in info["detail"].values())
    assert info["qc_pass"] and info["plates_used"] == [1, 2, 3, 4, 5]
    assert np.isclose(reward0 - reward, COST_WEIGHT * (info["cost"] - info0["cost"]) / 600.0)  # it costs its price only
    # nor does a later plate without controls: one passing plate is enough for the pair
    bare = {"type": "run_plate", "target": "off1", "design": {**SP_OFF, "controls_per_side": 0}}
    reward, info = after_the_panel(bare)
    assert info["qc_pass"] and np.isclose(info["quality"], 1.0) and info["plates_used"] == [1, 2, 3, 4, 5]
    assert np.isnan(info["zprime"]["6"])  # every plate's Z' is reported, used or not
    # a plain repeat changes nothing either: the earliest of two equal plates is kept
    reward, info = after_the_panel({"type": "run_plate", "target": "off1", "design": SP_OFF})
    assert info["qc_pass"] and np.isclose(info["quality"], 1.0) and info["plates_used"] == [1, 2, 3, 4, 5]
    assert all(np.isclose(t["achievable"], t["on"] - sp.P_FLOOR, atol=2e-3) for t in info["truth"].values())
    assert all(np.isclose(t["window"], t["on"] - max(t["off"].values()), atol=2e-3) for t in info["truth"].values())


def test_a_plate_that_fails_its_qc_supports_a_claim_only_while_no_passing_plate_tested_the_pair():
    def episode(off1_plates):
        env = make(NAME, seed=13)
        for target in ("on", "off2", "off3", "off4"):
            env.step({"type": "run_plate", "target": target, "design": SP_OFF})
        for design in off1_plates:
            assert "error" not in env.step({"type": "run_plate", "target": "off1", "design": design})[0]
        return env.step(perfect_report(env))[3]

    bare_high = {**SP_OFF, "controls_per_side": 0}  # 100 uM, no controls: fails plate QC
    low = {**SP_OFF, "top_conc_nM": 1000.0}  # 1 uM with controls: passes
    info = episode([bare_high])  # only a failing plate: it carries the support and the gate fails
    assert not info["qc_pass"] and all(np.isclose(d["tested_limit"], sp.P_FLOOR) for d in info["detail"].values())
    for order in ([bare_high, low], [low, bare_high]):  # a passing plate exists: the pair rests on it, whatever the order
        info = episode(order)
        assert info["qc_pass"] and all(np.isclose(d["tested_limit"], 6.0) for d in info["detail"].values())
        assert sorted(info["plates_used"]) == [1, 2, 3, 4, 5 + order.index(low)]
    info = episode([low, SP_OFF])  # two passing plates: the higher top dose counts
    assert info["qc_pass"] and all(np.isclose(d["tested_limit"], sp.P_FLOOR) for d in info["detail"].values())


# ---------------------------------------------------------------- policies

EXPECTED_CHANGE = {"on_target_range": "off.top_conc_nM", "single_point": "off.n_points",
                   "no_controls": "off.controls_per_side", "no_standards": "off.standards", "on_target_only": "plan",
                   "raw_ratio": "use_standards", "claim_past_range": "overclaim_log", "inactive_as_number": "bound_as_number"}


def changed_settings(trap):
    expert, other = POLICIES[NAME]["expert"](0).settings(), POLICIES[NAME][trap](0)
    assert isinstance(other, SelectivityPolicy) and set(other.settings()) == set(expert)
    return [k for k in expert if other.settings()[k] != expert[k]]


def test_each_trap_differs_from_the_expert_in_exactly_one_choice():
    expert = POLICIES[NAME]["expert"](0).settings()
    assert set(EXPECTED_CHANGE) == set(SINGLE_TRAPS) and len(SINGLE_TRAPS) >= 4
    for trap in SINGLE_TRAPS:
        assert changed_settings(trap) == [EXPECTED_CHANGE[trap]], trap
    assert expert["on.top_conc_nM"] < expert["off.top_conc_nM"] == sp.MAX_TOP_NM
    assert POLICIES[NAME]["on_target_range"](0).settings()["off.top_conc_nM"] == expert["on.top_conc_nM"]
    assert POLICIES[NAME]["claim_past_range"](0).settings()["overclaim_log"] == SP_OVERCLAIM_LOG > sp.SUPPORT_MARGIN


def test_the_labelled_interaction_is_one_single_mistake_plus_one_reporting_change():
    assert set(SP_INTERACTIONS) == set(TRAPS) - set(SINGLE_TRAPS) == {"single_point_as_curve"}
    for both, single in SP_INTERACTIONS.items():
        assert single in SINGLE_TRAPS
        assert sorted(changed_settings(both)) == sorted([EXPECTED_CHANGE[single], "hit_as_potency"])
    # the reporting change on its own: an active single point becomes a potency at the dose tested, nothing else moves
    screen = synthetic_plate({"hit": (6.0, 0.7, 1.0), "miss": (3.0, 1.0, 1.0)}, n=1, reps=3)
    honest, literal = SelectivityPolicy(), SelectivityPolicy(hit_as_potency=True)
    honest.data = literal.data = {"off1": screen}
    assert honest._calls()["off1"]["hit"] == {"above": 4.0, "top_nM": 100000.0}
    assert literal._calls()["off1"]["hit"] == {"pEC50": 4.0, "top_nM": 100000.0}
    assert literal._calls()["off1"]["miss"] == honest._calls()["off1"]["miss"] == {"below": 4.0, "top_nM": 100000.0}


def test_every_trap_traces_back_to_one_draft_rule():
    rules = [r for r in knowledge.load_rules() if r["env"] == NAME]
    assert sorted(t for r in rules for t in r["traps"]) == sorted(TRAPS)
    assert all(r["status"] == "draft" for r in rules)
    shared = {r["id"]: r["traps"] for r in rules if len(r["traps"]) != 1}
    assert shared == {"a-curve-where-you-need-a-potency": ["single_point", "single_point_as_curve"],
                      "compare-potencies-in-analyte-units": ["raw_ratio", "no_standards"]}
    assert knowledge.orphan_traps(knowledge.load_rules()) == []


def test_the_rules_join_to_a_card_of_this_environment():
    # results/envqc.json has no card for this environment until one is committed; with a card every rule is measured
    rules = [r for r in knowledge.load_rules() if r["env"] == NAME]
    assert all(r["teaches"] == "not measured" for r in knowledge.coverage(rules, {}))
    card = envqc.qc_card(NAME, n_seeds=3)
    assert [r["policy"] for r in card["rows"]] == list(POLICIES[NAME]) and card["rows"][0]["qc_pass_rate"] is not None
    rows = knowledge.coverage(rules, {NAME: card})
    assert len(rows) == len(TRAPS) and all(r["verdict"] is not None for r in rows)


def test_the_null_claims_untested_and_the_random_policy_runs_valid_plates_and_always_claims():
    env, steps = trajectory("null", 9)
    assert len(steps) == 1 and env.plates == [] and not steps[0][4]["qc_pass"]
    assert steps[0][2] == -0.5 and all(d["score"] == -0.5 for d in steps[0][4]["detail"].values())  # worse than silence
    kinds = set()
    for seed in range(40):
        policy = POLICIES[NAME]["random"](seed)
        assert isinstance(policy, SelectivityRandom) and 1 <= len(policy.plan) <= sp.PLATE_BUDGET
        assert sp.check_plate(policy.on_design, 4)[1] is None and sp.check_plate(policy.off_design, 4)[1] is None
        env, steps = trajectory("random", seed)
        claims = steps[-1][0]["windows"]
        assert set(claims) == set(env.pharmacology) and len(env.plates) == len(policy.plan)
        assert all(SP_RANDOM_WINDOW[0] <= c["window"] <= SP_RANDOM_WINDOW[1] for c in claims.values())
        kinds |= {c["at_least"] for c in claims.values()}
        assert all(d["claim"] is not None for d in steps[-1][4]["detail"].values())  # the row exercises the scoring rule
    assert kinds == {True, False}


def test_no_top_level_name_is_defined_twice_in_the_shared_modules():
    # two environments once defined the same names in policies.py and the later one silently replaced the earlier
    src = Path(sp.__file__).resolve().parents[1]
    for module in ("policies.py", "analysis.py", "chem.py"):
        names = []
        for node in ast.parse((src / module).read_text()).body:
            if isinstance(node, (ast.FunctionDef, ast.ClassDef)):
                names.append(node.name)
            elif isinstance(node, ast.Assign):
                names += [t.id for t in node.targets if isinstance(t, ast.Name)]
            elif isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name):
                names.append(node.target.id)
            elif isinstance(node, (ast.Import, ast.ImportFrom)):
                names += [(a.asname or a.name).split(".")[0] for a in node.names]
        twice = sorted(n for n, count in collections.Counter(names).items() if count > 1)
        assert twice == [], (module, twice)


FRESH = range(30000, 30060)


@pytest.fixture(scope="module")
def fresh():
    return {name: np.array([play(name, s)[0] for s in FRESH]) for name in POLICIES[NAME]}


def test_the_expert_beats_the_null_by_a_clear_margin_on_fresh_seeds(fresh):
    assert fresh["null"].mean() <= 0.3 and fresh["expert"].mean() >= 0.5
    assert fresh["expert"].mean() > fresh["null"].mean() + 0.4
    assert np.mean(fresh["expert"] > fresh["null"]) > 0.95
    assert np.all(fresh["null"] == -0.5)  # the untested claim; an empty report scores 0.0 and is the real floor
    assert np.mean(fresh["expert"] > 0.0) > 0.95  # ... which the expert also clears


def test_the_literal_single_point_trap_scores_below_the_honest_one_on_fresh_seeds(fresh):
    # reading an active single point as a curve must not pay better than leaving it unanswered
    honest, literal = fresh["single_point"], fresh["single_point_as_curve"]
    assert honest.mean() > literal.mean() + 0.05
    assert np.mean(honest > literal) + 0.5 * np.mean(honest == literal) > 0.6


def test_no_trap_and_no_reference_policy_outscores_the_expert_on_fresh_seeds(fresh):
    for name, rewards in fresh.items():
        if name != "expert":
            assert fresh["expert"].mean() > rewards.mean() + 0.1, name


class Probe:
    """Cheap token testing, then a fixed claim for every compound: the answers a reward must not pay."""

    def __init__(self, claim, plan=sp.OFF_TARGETS, on_curve=False):
        self.claim, self.queue, self.on_curve, self.on = claim, list(plan), on_curve, {}

    def act(self, obs):
        if obs.get("phase") == "data" and obs["target"] == "on":
            self.on = {c: v["pEC50"] for c, v in estimate_potency(obs["wells"]).items() if "pEC50" in v}
        if obs.get("phase") == "start":
            self.compounds = obs["compounds"]
            if self.on_curve:
                return {"type": "run_plate", "target": "on", "design": SP_ON}
        if self.queue:
            return {"type": "run_plate", "target": self.queue.pop(0), "design": SCREEN}
        return {"type": "report", "windows": {c: w for c in self.compounds if (w := self.claim(self.on.get(c))) is not None}}


PROBES = {
    **{f"untested_bound_{x}": (lambda x=x: Probe(lambda on: {"window": x, "at_least": True}, plan=()))
       for x in (-5.0, 0.0, 2.0, 4.0)},
    **{f"screened_bound_{x}": (lambda x=x: Probe(lambda on: {"window": x, "at_least": True}))
       for x in (-5.0, 0.0, 2.0, 3.0, 3.5, 4.0, 4.5)},
    **{f"screened_number_{x}": (lambda x=x: Probe(lambda on: {"window": x, "at_least": False}))
       for x in (1.0, 2.0, 3.0, 3.5, 4.0, 4.5)},
    # the strongest of them: a real on-target curve, token screens that are never looked at, the limit as the bound
    "always_at_the_limit": lambda: Probe(lambda on: None if on is None else {"window": on - sp.P_FLOOR, "at_least": True},
                                         on_curve=True),
}


def test_no_constant_answer_and_no_always_safe_bound_scores_near_the_expert(fresh):
    seeds = FRESH[:30]
    expert = fresh["expert"][:30].mean()
    means = {name: np.mean([run_episode(make(NAME, seed=s), build())[0] for s in seeds]) for name, build in PROBES.items()}
    assert all(m < expert - 0.25 for m in means.values()), means
    assert all(means[k] == -0.5 for k in means if k.startswith("untested"))  # worse than an empty report
    assert all(v < 0.0 for k, v in means.items() if k.startswith("screened")), means  # no on-target plate: unsupported
    assert means["always_at_the_limit"] < 0.25, means


# ---------------------------------------------------------------- analysis functions

def synthetic_plate(curves, top=100000.0, n=10, fold=4.0, reps=1, controls=6, standards=True, k=30.0, low=6.0, high=150.0):
    """A noise-free plate: curves maps compound id to (pEC50 or None for inactive, emax, hill)."""
    kit = Kit(k=k)
    wells = []
    for _ in range(controls):
        wells.append({"role": "neg", "id": None, "conc_nM": None, "std_nM": None, "ratio": float(kit.ratio(low))})
        wells.append({"role": "pos", "id": None, "conc_nM": None, "std_nM": None, "ratio": float(kit.ratio(high))})
    if standards:
        wells += [{"role": "std", "id": None, "conc_nM": None, "std_nM": s, "ratio": float(kit.ratio(s))}
                  for s in STD_SERIES_NM]
    for cid, (p, emax, slope) in curves.items():
        for j in range(n):
            conc = top / fold ** j
            frac = 0.0 if p is None else emax * float(hill(conc, 10.0 ** (9.0 - p), slope))
            wells += [{"role": "cpd", "id": cid, "conc_nM": conc, "std_nM": None,
                       "ratio": float(kit.ratio(low + (high - low) * frac))}] * reps
    return wells


def test_capped_fit_matches_the_free_fit_on_a_complete_curve_and_holds_the_ceiling():
    logx = np.linspace(-1, 4, 10)
    y = 0.05 + 0.8 * hill(10.0 ** logx, 10.0 ** 1.7, 1.2)
    free, capped = fit_4pl(logx, y), fit_4pl_capped(logx, y, 1.0)
    assert set(capped) == set(free) and abs(capped["logec50"] - 1.7) <= 0.02 and np.isclose(capped["top"], free["top"])
    held = fit_4pl_capped(logx, 2.0 * y, 1.0)  # the data rise to 1.7 but the plateau may not pass 1.0
    assert held["top"] <= 1.0 + 1e-9 and held["sse"] > fit_4pl(logx, 2.0 * y)["sse"]
    assert fit_4pl_capped([0, 1, 2], [1, 2, 3], 1.0) is None
    # an incomplete curve with noise: the free plateau runs away, the capped midpoint stays near the truth
    rng = np.random.default_rng(1)
    logx = np.log10(100000.0 / 4.0 ** np.arange(10))
    misses = {"free": [], "capped": []}
    for _ in range(40):
        y = 0.9 * hill(10.0 ** logx, 10.0 ** 4.8, 1.0) + rng.normal(0.0, 0.03, 10)
        misses["free"].append(abs(fit_4pl(logx, y, step=0.05)["logec50"] - 4.8))
        misses["capped"].append(abs(fit_4pl_capped(logx, y, 1.0, step=0.05)["logec50"] - 4.8))
    assert np.mean(misses["capped"]) < np.mean(misses["free"]) and np.mean(misses["capped"]) < 0.2


def test_activity_fraction_is_zero_at_the_negative_and_one_at_the_positive_control():
    wells = synthetic_plate({"full": (7.0, 1.0, 1.0), "flat": (None, 1.0, 1.0)}, reps=2)
    frac = activity_fraction(wells)
    assert list(frac["full"]) == sorted(frac["full"], reverse=True) and len(frac["full"]) == 10
    assert np.isclose(frac["full"][100000.0], 1.0, atol=0.02) and abs(frac["flat"][100000.0]) < 0.02
    assert np.isclose(frac["full"][100000.0 / 4 ** 5], 0.5, atol=0.05)  # 98 nM against an EC50 of 100 nM
    assert activity_fraction(synthetic_plate({"full": (7.0, 1.0, 1.0)}, controls=0)) == {}
    raw = activity_fraction(wells, use_standards=False)["full"][100000.0 / 4 ** 5]
    assert raw > 0.7  # the raw ratio saturates early: the same well reads far more than half way


def test_estimate_potency_gives_a_number_a_bound_or_neither():
    wells = synthetic_plate({"potent": (7.0, 1.0, 1.0), "partial": (6.5, 0.6, 1.2), "weak": (4.5, 0.9, 1.0),
                             "foot": (3.6, 1.0, 1.0)})
    calls = estimate_potency(wells)
    assert abs(calls["potent"]["pEC50"] - 7.0) < 0.1 and abs(calls["partial"]["pEC50"] - 6.5) < 0.1
    assert abs(calls["weak"]["pEC50"] - 4.5) < 0.2  # an incomplete curve, plateau held at the control
    assert calls["foot"] == {"below": 4.0, "top_nM": 100000.0}  # 28% at the top dose: a bound, never a number
    assert all(c["top_nM"] == 100000.0 for c in calls.values())
    flat = estimate_potency(synthetic_plate({"flat": (None, 1.0, 1.0)}, top=3000.0))
    assert np.isclose(flat["flat"]["below"], 9 - np.log10(3000.0)) and "pEC50" not in flat["flat"]
    assert ACTIVE_FRACTION == 0.5


def test_a_single_point_can_bound_but_cannot_give_a_potency():
    screen = synthetic_plate({"hit": (6.0, 0.7, 1.0), "weak_hit": (4.1, 1.0, 1.0), "miss": (3.0, 1.0, 1.0)}, n=1, reps=3)
    calls = estimate_potency(screen)
    assert calls["hit"] == {"above": 4.0, "top_nM": 100000.0} and calls["weak_hit"] == {"above": 4.0, "top_nM": 100000.0}
    assert calls["miss"] == {"below": 4.0, "top_nM": 100000.0}
    saturated = estimate_potency(synthetic_plate({"c": (9.5, 1.0, 1.0)}, top=100000.0, n=6))  # flat at the top plateau
    assert "above" in saturated["c"]


def test_without_controls_nothing_can_be_called_inactive():
    wells = synthetic_plate({"potent": (7.0, 1.0, 1.0), "flat": (None, 1.0, 1.0)}, controls=0)
    calls = estimate_potency(wells)
    assert set(calls) == {"potent"} and abs(calls["potent"]["pEC50"] - 7.0) < 0.1


def test_raw_ratio_fits_read_too_potent():
    wells = synthetic_plate({"c": (7.0, 1.0, 1.0)})
    assert estimate_potency(wells, use_standards=False)["c"]["pEC50"] > estimate_potency(wells)["c"]["pEC50"] + 0.3
    no_std = synthetic_plate({"c": (7.0, 1.0, 1.0)}, standards=False)  # no standards on the plate: raw is all there is
    assert estimate_potency(no_std) == estimate_potency(no_std, use_standards=False)


def test_selectivity_window_is_a_number_a_bound_or_nothing():
    below = {"below": 4.0, "top_nM": 100000.0}
    assert selectivity_window(8.0, [below] * 4) == {"window": 4.0, "at_least": True}
    measured = [{"pEC50": 6.5, "top_nM": 100000.0}, below, {"pEC50": 4.6, "top_nM": 100000.0}, below]
    assert selectivity_window(8.0, measured) == {"window": 1.5, "at_least": False}
    narrow = {"below": 9 - np.log10(3000.0), "top_nM": 3000.0}  # an off-target tested only to 3 uM
    out = selectivity_window(8.0, [{"pEC50": 4.6, "top_nM": 100000.0}, narrow, below, below])
    assert out["at_least"] and np.isclose(out["window"], 8.0 - (9 - np.log10(3000.0)))
    assert selectivity_window(8.0, [below, below, below, None]) is None  # an untested off-target
    assert selectivity_window(8.0, [below, below, below, {"above": 4.0, "top_nM": 100000.0}]) is None
    assert selectivity_window(8.0, [below, below, below, {}]) is None and selectivity_window(8.0, []) is None
    for bad in (None, float("nan"), float("inf"), True):
        assert selectivity_window(bad, [below] * 4) is None


def test_task_id_and_constructor_limits():
    env = make(f"{NAME}/randomized/17")
    assert env.seed == 17 and env.profile.name == "randomized" and len(env.pharmacology) == 4
    for bad in (0, 9):
        with pytest.raises(ValueError):
            sp.SelectivityPanel(seed=0, n_compounds=bad)
    assert len(sp.SelectivityPanel(seed=0, n_compounds=2).reset()["compounds"]) == 2
