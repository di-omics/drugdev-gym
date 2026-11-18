"""Environment QC: validate a reward the way an assay is validated, before scaling it.

For each brick, run the expert and every trap over the same seeds and ask whether the
reward separates them. Z' between two reward distributions is the strict reading (3 SD of
clearance); the paired win rate (same seed, same hidden truth) is the lenient one. A trap
the reward cannot see is reported as "blind", not hidden.
"""
from __future__ import annotations

import numpy as np

from .policies import POLICIES, run_episode
from .registry import make


def evaluate(brick: str, policy: str, seeds, rig="nominal") -> tuple[np.ndarray, np.ndarray]:
    """Rewards and QC-gate outcomes (NaN where a brick has no gate) for one policy over seeds."""
    rewards, gates = [], []
    for seed in seeds:
        reward, info = run_episode(make(brick, seed=seed, rig=rig), POLICIES[brick][policy](seed))
        rewards.append(reward)
        gates.append(float(info["qc_pass"]) if "qc_pass" in info else np.nan)
    return np.array(rewards), np.array(gates)


def verdict(gap: float, z: float, wins: float) -> str:
    """clean > visible > weak > blind. See the legend in card_markdown."""
    if gap > 0 and z >= 0.5:
        return "clean"
    if wins >= 0.9:
        return "visible"
    if wins >= 0.6 or gap >= 0.1:
        return "weak"
    return "blind"


def separation(expert: np.ndarray, other: np.ndarray) -> dict:
    gap = float(expert.mean() - other.mean())
    spread = expert.std(ddof=1) + other.std(ddof=1)
    z = float(1 - 3 * spread / abs(gap)) if gap != 0 else float("-inf")
    wins = float(np.mean(expert > other) + 0.5 * np.mean(expert == other))
    return {"gap": gap, "zprime": z, "win_rate": wins, "verdict": verdict(gap, z, wins)}


def qc_card(brick: str, n_seeds: int = 200, rig="nominal", first_seed: int = 0) -> dict:
    seeds = range(first_seed, first_seed + n_seeds)
    runs = {name: evaluate(brick, name, seeds, rig) for name in POLICIES[brick]}
    expert = runs["expert"][0]
    rows = []
    for name, (rewards, gates) in runs.items():
        row = {"policy": name, "mean": float(rewards.mean()), "sd": float(rewards.std(ddof=1)),
               "p10": float(np.percentile(rewards, 10)), "p90": float(np.percentile(rewards, 90)),
               "qc_pass_rate": None if np.isnan(gates).all() else float(np.nanmean(gates))}
        if name != "expert":
            row.update(separation(expert, rewards))
        rows.append(row)
    return {"brick": brick, "rig": rig if isinstance(rig, str) else rig.name, "n_seeds": n_seeds,
            "first_seed": first_seed, "rows": rows}


def card_markdown(card: dict) -> str:
    lines = [f"# Environment QC: {card['brick']}", "",
             f"Rig: {card['rig']}. Seeds: {card['first_seed']} to {card['first_seed'] + card['n_seeds'] - 1} "
             f"({card['n_seeds']} episodes per policy, same hidden truth across policies).", "",
             "| policy | mean reward | sd | p10 | p90 | QC gate pass | gap to expert | Z' vs expert | paired win rate | verdict |",
             "|---|---|---|---|---|---|---|---|---|---|"]
    for r in card["rows"]:
        gate = "n/a" if r["qc_pass_rate"] is None else f"{100 * r['qc_pass_rate']:.0f}%"
        if "verdict" in r:
            z = "n/a" if not np.isfinite(r["zprime"]) else f"{r['zprime']:.2f}"
            tail = f"{r['gap']:.3f} | {z} | {100 * r['win_rate']:.0f}% | {r['verdict']}"
        else:
            tail = "reference | | | |"
        lines.append(f"| {r['policy']} | {r['mean']:.3f} | {r['sd']:.3f} | {r['p10']:.3f} | {r['p90']:.3f} | {gate} | {tail}")
    lines += ["", "Verdicts: clean = Z' >= 0.5 between the two reward distributions; visible = the expert wins on "
              ">= 90% of seeds; weak = the expert wins on 60 to 90% of seeds, or the mean gap is >= 0.1 (punished on "
              "average, not reliably per episode); blind = neither (the reward cannot see this mistake).", ""]
    return "\n".join(lines)
