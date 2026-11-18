import json

import numpy as np
import pytest

from drugdevgym import ENVS, make
from drugdevgym.policies import DR_EXPERT, POLICIES, ST_EXPERT, run_episode

HIDDEN_WORDS = ("truth", "pec50", "emax", "half_lives", "per_k", "sticky", "basal")


def play(env_name, policy, seed, rig="nominal"):
    return run_episode(make(env_name, seed=seed, rig=rig), POLICIES[env_name][policy](seed))


@pytest.mark.parametrize("env_name", list(ENVS))
def test_same_seed_same_episode(env_name):
    assert play(env_name, "expert", 3)[0] == play(env_name, "expert", 3)[0]
    assert play(env_name, "expert", 3)[0] != play(env_name, "expert", 4)[0]


EVERY_POLICY = [(name, policy) for name in ENVS for policy in POLICIES[name]]


@pytest.mark.parametrize("env_name,policy_name", EVERY_POLICY)
@pytest.mark.parametrize("rig", ["nominal", "harsh"])
def test_observations_are_json_and_never_carry_truth_fields(env_name, policy_name, rig):
    env, policy = make(env_name, seed=5, rig=rig), POLICIES[env_name][policy_name](5)
    obs, done = env.reset(), False
    while not done:
        json.dumps(obs)  # the whole observation is plain JSON
        text = json.dumps({k: v for k, v in obs.items() if k != "spec"}).lower()  # spec is public task text
        assert not any(word in text for word in HIDDEN_WORDS), (env_name, policy_name)
        obs, reward, done, info = env.step(policy.act(obs))
        assert done or info == {}
    assert "truth" in info and reward <= 1.0
    assert "truth" not in json.dumps(obs)


@pytest.mark.parametrize("env_name", list(ENVS))
def test_step_after_the_end_raises(env_name):
    env = make(env_name, seed=0)
    run_episode(env, POLICIES[env_name]["expert"](0))
    with pytest.raises(RuntimeError):
        env.step({"type": "report"})


def test_hidden_truth_depends_on_the_seed_not_the_policy():
    assert play("dose_response", "expert", 9)[1]["truth"] == play("dose_response", "null", 9)[1]["truth"]
    assert play("stability", "expert", 9)[1]["truth"] == play("stability", "two_point", 9)[1]["truth"]


def test_invalid_designs_cost_nothing_and_five_end_the_episode():
    env = make("dose_response", seed=0)
    env.reset()
    for bad in ({}, {**DR_EXPERT, "n_points": 99}, {**DR_EXPERT, "mix": "yes"}, {**DR_EXPERT, "extra": 1}):
        obs, reward, done, _ = env.step({"type": "run_plate", "design": bad})
        assert "error" in obs and reward == 0.0 and not done and env.cost == 0.0
    obs, reward, done, info = env.step({"type": "dance"})
    assert done and info == {"truncated": True}


def test_designs_that_do_not_fit_the_labware_are_rejected():
    env = make("dose_response", seed=0)
    too_many = {**DR_EXPERT, "replicates": 4}
    assert "wells" in env.step({"type": "run_plate", "design": too_many})[0]["error"]
    overflow = {**DR_EXPERT, "transfer_ul": 100.0, "dilution_factor": 10.0}
    assert "ladder well" in env.step({"type": "run_plate", "design": overflow})[0]["error"]
    st = make("stability", seed=0)
    assert "error" in st.step({"type": "run", "design": {**ST_EXPERT, "timepoints_min": [0]}})[0]


def test_dose_response_qc_gate_halves_an_ungated_answer():
    reward, info = play("dose_response", "few_controls", 2)
    assert not info["qc_pass"]
    assert np.isclose(reward, 0.5 * info["accuracy"] - 0.15 * info["cost"] / 200.0)


def test_a_starved_ladder_ruins_the_curve():
    rewards = [play("dose_response", "thin_ladder", s)[0] for s in range(10)]
    assert np.mean(rewards) < 0.2


@pytest.mark.parametrize("env_name,weaker", [("dose_response", "null"), ("stability", "null"),
                                          ("assay_window", "random"), ("lead_opt", "null")])
def test_expert_beats_the_floor_policy(env_name, weaker):
    seeds = range(12)
    expert = np.mean([play(env_name, "expert", s)[0] for s in seeds])
    floor = np.mean([play(env_name, weaker, s)[0] for s in seeds])
    assert expert > floor + 0.1


def test_lead_opt_quality_is_relative_to_the_library_best():
    reward, info = play("lead_opt", "expert", 1)
    assert 0.0 <= info["quality"] <= 1.0 and info["best"] > 0
    assert np.isclose(reward, info["quality"] - 0.15 * info["cost"] / 2000.0)


def test_task_ids_round_trip():
    from drugdevgym import task_id
    assert task_id("stability", "randomized", 17) == "stability/randomized/17"
    env = make(task_id("stability", "randomized", 17))
    assert env.seed == 17 and env.rig_mode == "randomized" and env.profile.name == "randomized"
    with pytest.raises(KeyError):
        make("alchemy")


def test_every_env_documents_its_actions_and_reward():
    for name in ENVS:
        spec = make(name).spec()
        assert spec["env"] == name and spec["actions"] and "reward" in spec
        json.dumps(spec)


def test_malformed_numbers_are_rejected_not_crashes():
    huge = 10 ** 400
    for bad in ([[0], [5]], [[0, 5, 15]], [{}, {}], [0, float("nan")], [0, huge]):
        obs, _, done, _ = make("stability", seed=1).step({"type": "run", "design": {**ST_EXPERT, "timepoints_min": bad}})
        assert "error" in obs and not done, bad
    for field, bad in (("n_points", float("nan")), ("n_points", float("inf")), ("n_points", huge),
                       ("controls_per_side", float("nan")), ("top_conc_nM", float("nan"))):
        obs, _, done, _ = make("dose_response", seed=1).step({"type": "run_plate", "design": {**DR_EXPERT, field: bad}})
        assert "error" in obs and not done, (field, bad)


def test_an_absurd_answer_scores_zero_instead_of_crashing():
    env = make("dose_response", seed=1)
    _, reward, done, info = env.step({"type": "report", "pEC50": {"cpd1": 10 ** 400, "cpd2": float("nan")}})
    assert done and info["accuracy"] == 0.0


def test_plates_used_must_name_real_plates():
    env = make("dose_response", seed=1)
    env.step({"type": "run_plate", "design": DR_EXPERT})
    for bad in ([1.5], ["1"], [3], [float("nan")]):
        obs, _, done, _ = env.step({"type": "report", "pEC50": {}, "plates_used": bad})
        assert "error" in obs and not done, bad
    _, _, done, info = env.step({"type": "report", "pEC50": {}, "plates_used": [1.0]})
    assert done and info["qc_pass"]


def test_constructor_limits():
    from drugdevgym.envs import DoseResponse, Stability
    for cls in (DoseResponse, Stability):
        with pytest.raises(ValueError):
            cls(seed=0, n_compounds=9)
    with pytest.raises(ValueError):
        make("dose_response", seed=-1)


def test_a_saturated_assay_earns_nothing():
    # seed 329 has a very productive cell line: 80k cells for 180 min is far off the standard curve
    env = make("assay_window", seed=329)
    _, reward, done, info = env.step({"type": "commit", "cells_k": 80.0, "stim_min": 180.0, "lid": True,
                                      "avoid_edges": True})
    assert done and reward <= 0.0 and not info["qc_pass"]
    assert info["zprime_raw_ratio"] > 0.5  # the raw ratio would have passed it


def test_the_gate_bypass_exploit_is_still_open_and_on_the_card():
    assert "gate_bypass" in POLICIES["dose_response"]
    _, info = play("dose_response", "gate_bypass", 3)
    assert info["qc_pass"] and len(info["zprime"]) == 1


def test_worklist_names_what_goes_into_each_well():
    env = make("dose_response", seed=0)
    run_episode(env, POLICIES["dose_response"]["expert"](0))
    labels = {s.get("what") for s in env.step_log if s["op"] == "fill"}
    assert {"cells", "diluent", "vehicle", "reference agonist", "detection mix"} <= labels
    assert any(str(x).startswith("analyte standard") for x in labels)
