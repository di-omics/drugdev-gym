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


@pytest.mark.parametrize("env_name", list(ENVS))
@pytest.mark.parametrize("rig", ["nominal", "randomized", "harsh"])
def test_observations_are_json_and_never_leak_the_truth(env_name, rig):
    env, policy = make(env_name, seed=5, rig=rig), POLICIES[env_name]["expert"](5)
    obs, done = env.reset(), False
    while not done:
        json.dumps(obs)  # the whole observation is plain JSON
        text = json.dumps({k: v for k, v in obs.items() if k != "spec"}).lower()  # spec is public task text
        assert not any(word in text for word in HIDDEN_WORDS), env_name
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
    env = make("stability/randomized/17")
    assert env.seed == 17 and env.rig_mode == "randomized" and env.profile.name == "randomized"
    with pytest.raises(KeyError):
        make("alchemy")


def test_every_env_documents_its_actions_and_reward():
    for name in ENVS:
        spec = make(name).spec()
        assert spec["env"] == name and spec["actions"] and "reward" in spec
        json.dumps(spec)
