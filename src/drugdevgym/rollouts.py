"""Rollouts: turn environments into training data.

An environment is only useful for post-training if it produces data cheaply. This module
runs policies through environments and writes what happened, in two shapes:

- episodes: every step of an episode (observation, action, reward) plus the final info with
  the hidden truth. This is the raw record, for RL or for analysis.
- sft: prompt and response pairs for supervised fine-tuning, one per decision, taken only
  from episodes that scored at least min_reward. The prompt is the task and everything the
  agent had seen and done up to that decision; the response is the action it took. The
  hidden truth is never in a prompt.

The policies in this repo are scripted, so these rollouts are synthetic demonstrations: a
model fine-tuned on them learns the scripted expert and nothing more. Their use is as a
starting point and as negatives (the traps come with their rewards); real post-training
would let the model act in the environment and learn from the reward itself.
"""
from __future__ import annotations

import copy
import json
import math

from .policies import POLICIES
from .registry import make, task_id


def _plain(value):
    """JSON-safe copy: NaN and infinity become null."""
    if isinstance(value, dict):
        return {str(k): _plain(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_plain(v) for v in value]
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if hasattr(value, "item") and not isinstance(value, (str, bytes)):
        return _plain(value.item())
    return value


def episode(env_name: str, policy_name: str, seed: int, rig: str = "nominal") -> dict:
    """One full episode as a JSON-ready record."""
    env = make(env_name, seed=seed, rig=rig)
    policy = POLICIES[env_name][policy_name](seed)
    obs = env.reset()
    spec = obs.pop("spec", None)
    steps, total = [], 0.0
    while True:
        action = policy.act(obs)
        seen = copy.deepcopy(obs)
        obs, reward, done, info = env.step(action)
        obs.pop("spec", None)
        total += reward
        steps.append({"observation": _plain(seen), "action": _plain(action), "reward": reward})
        if done:
            break
    return {"task": task_id(env_name, rig, seed), "env": env_name, "policy": policy_name, "spec": _plain(spec),
            "steps": steps, "reward": total, "info": _plain(info)}


def sft_pairs(record: dict) -> list[dict]:
    """Prompt and response pairs from one episode. No truth, no reward, no policy name in the prompt.

    Each prompt carries the episode so far: history holds every earlier observation with the action
    taken on it. A report at the end of several runs is computed from all of them, so a prompt with
    only the last observation would ask for numbers from data it does not show.
    """
    pairs, history = [], []
    for i, step in enumerate(record["steps"]):
        prompt = {"task": record["spec"], "step": i + 1, "history": list(history), "observation": step["observation"]}
        pairs.append({"task": record["task"], "prompt": json.dumps(prompt, sort_keys=True),
                      "response": json.dumps(step["action"], sort_keys=True), "episode_reward": record["reward"]})
        history.append({"observation": step["observation"], "action": step["action"]})
    return pairs


def generate(env_names, policy_names=None, n_seeds: int = 100, first_seed: int = 100000, rig: str = "nominal"):
    """Yield episode records. Training seeds start far from the seeds any card is reported on."""
    for env_name in env_names:
        for policy_name in (policy_names or list(POLICIES[env_name])):
            if policy_name not in POLICIES[env_name]:
                continue
            for seed in range(first_seed, first_seed + n_seeds):
                yield episode(env_name, policy_name, seed, rig)


def write(path, env_names, policy_names=None, n_seeds: int = 100, first_seed: int = 100000, rig: str = "nominal",
          fmt: str = "episodes", min_reward: float = 0.0) -> dict:
    """Write JSON lines to path. Returns counts and the mean reward per (env, policy)."""
    if fmt not in ("episodes", "sft"):
        raise ValueError("fmt must be episodes or sft")
    lines, kept, rewards = 0, 0, {}
    with open(path, "w") as out:
        for record in generate(env_names, policy_names, n_seeds, first_seed, rig):
            rewards.setdefault((record["env"], record["policy"]), []).append(record["reward"])
            if fmt == "episodes":
                out.write(json.dumps(record) + "\n")
                lines += 1
                kept += 1
            elif record["reward"] >= min_reward:
                kept += 1
                for pair in sft_pairs(record):
                    out.write(json.dumps(pair) + "\n")
                    lines += 1
    return {"lines": lines, "episodes_kept": kept, "episodes_run": sum(len(v) for v in rewards.values()),
            "mean_reward": {f"{e}/{p}": sum(v) / len(v) for (e, p), v in rewards.items()}}
