import json

import pytest

from drugdevgym import ENVS, rollouts
from drugdevgym.policies import POLICIES

HIDDEN = ("truth", "pec50\": ", "should_repeat", "fault")


@pytest.mark.parametrize("env_name", list(ENVS))
def test_an_episode_is_a_complete_json_record(env_name):
    record = rollouts.episode(env_name, "expert", 100000)
    assert json.loads(json.dumps(record)) == record
    assert record["task"] == f"{env_name}/nominal/100000" and record["steps"]
    assert "truth" in record["info"] and record["reward"] <= 1.0
    assert record == rollouts.episode(env_name, "expert", 100000)


@pytest.mark.parametrize("env_name", list(ENVS))
def test_sft_prompts_never_contain_the_truth_or_the_outcome(env_name):
    for policy in POLICIES[env_name]:
        record = rollouts.episode(env_name, policy, 100001)
        pairs = rollouts.sft_pairs(record)
        assert len(pairs) == len(record["steps"])
        for i, pair in enumerate(pairs):
            prompt = json.loads(pair["prompt"])
            assert set(prompt) == {"task", "step", "history", "observation"} and json.loads(pair["response"])["type"]
            text = json.dumps([prompt["observation"], prompt["history"]])
            assert not any(word in text for word in HIDDEN), (env_name, policy)
            assert policy not in pair["prompt"] or policy in ("random", "null", "expert")
            # the prompt carries the episode so far: every earlier observation with the action taken on it
            assert prompt["history"] == [{"observation": s["observation"], "action": s["action"]}
                                         for s in record["steps"][:i]]


def test_a_report_after_several_runs_is_asked_from_all_the_data_it_used():
    record = rollouts.episode("selectivity_panel", "expert", 100000)
    last = json.loads(rollouts.sft_pairs(record)[-1]["prompt"])
    plates = [h["observation"] for h in last["history"]] + [last["observation"]]
    targets = [o["target"] for o in plates if "wells" in o]
    assert len(set(targets)) == 5 and json.loads(rollouts.sft_pairs(record)[-1]["response"])["type"] == "report"


def test_write_filters_sft_by_reward_and_counts(tmp_path):
    path = tmp_path / "all.jsonl"
    stats = rollouts.write(path, ["stability"], ["expert", "null"], n_seeds=3)
    assert stats["lines"] == stats["episodes_run"] == 6 and len(path.read_text().splitlines()) == 6
    sft = tmp_path / "sft.jsonl"
    kept = rollouts.write(sft, ["stability"], ["expert", "null"], n_seeds=3, fmt="sft", min_reward=0.6)
    assert kept["episodes_kept"] == 3 and kept["mean_reward"]["stability/expert"] > kept["mean_reward"]["stability/null"]
    rows = [json.loads(line) for line in sft.read_text().splitlines()]
    assert rows and all(r["episode_reward"] >= 0.6 for r in rows)
    with pytest.raises(ValueError):
        rollouts.write(sft, ["stability"], fmt="parquet")
