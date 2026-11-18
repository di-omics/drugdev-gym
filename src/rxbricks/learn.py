"""A small learner: the cross-entropy method over dose-response protocol choices.

It starts from uniform random protocols and learns from reward alone. Training on the
nominal rig versus on randomized rigs, then testing on a held-out harsh rig, is a
sim-to-sim stand-in for the sim-to-real question. No hardware is involved.
"""
from __future__ import annotations

import numpy as np

from .policies import DR_CHOICES, DosePolicy, run_episode
from .registry import make


def design_reward(design: dict, seeds, rig) -> float:
    policy_std = bool(design["standards"])
    return float(np.mean([run_episode(make("dose_response", seed=s, rig=rig), DosePolicy(design, policy_std))[0]
                          for s in seeds]))


def greedy(probs: dict) -> dict:
    return {k: DR_CHOICES[k][int(np.argmax(p))] for k, p in probs.items()}


def cem(train_rig="nominal", iterations: int = 12, population: int = 40, seeds_per_design: int = 4,
        elite_frac: float = 0.2, smoothing: float = 0.6, seed: int = 0) -> dict:
    """Returns the learning curve and the final greedy design."""
    rng = np.random.default_rng([seed, 7])
    probs = {k: np.full(len(v), 1.0 / len(v)) for k, v in DR_CHOICES.items()}
    next_seed = 1_000_000 * (seed + 1)  # training seeds never overlap the QC or test seeds
    curve = []
    for it in range(iterations):
        designs = [{k: v[int(rng.choice(len(v), p=probs[k]))] for k, v in DR_CHOICES.items()}
                   for _ in range(population)]
        scores = []
        for d in designs:
            seeds = range(next_seed, next_seed + seeds_per_design)
            next_seed += seeds_per_design
            scores.append(design_reward(d, seeds, train_rig))
        order = np.argsort(scores)[::-1][: max(2, int(elite_frac * population))]
        for k, values in DR_CHOICES.items():
            freq = np.array([sum(designs[i][k] == v for i in order) for v in values], float) + 0.25
            probs[k] = (1 - smoothing) * probs[k] + smoothing * freq / freq.sum()
        curve.append({"iteration": it + 1, "population_mean": float(np.mean(scores)),
                      "elite_mean": float(np.mean([scores[i] for i in order]))})
    return {"train_rig": train_rig, "curve": curve, "design": greedy(probs)}


def transfer_table(designs: dict[str, dict], rigs=("nominal", "randomized", "harsh"), n_seeds: int = 100,
                   first_seed: int = 5000) -> list[dict]:
    """Mean reward of each named design on each rig, on seeds no learner has seen."""
    seeds = range(first_seed, first_seed + n_seeds)
    return [{"design": name, **{rig: design_reward(d, seeds, rig) for rig in rigs}} for name, d in designs.items()]
