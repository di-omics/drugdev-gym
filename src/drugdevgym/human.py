"""Capture a human expert's judgment on review cases, and score it against planted truth.

export_cases writes what a reviewer may see (never the truth). The expert judges each case:
accept or repeat, the most likely cause, how sure they are, whether the data looks like
something a real bench would produce, and a note in their own words. score_judgments replays
those calls through the same environment a model would face, next to the scripted QC rules.
The notes are the raw material for new rules; "not realistic" flags are simulator to-dos.
"""
from __future__ import annotations

from collections import defaultdict

import numpy as np

from .policies import POLICIES, run_episode
from .registry import make, task_id

REVIEW_ENVS = ("plate_review", "timecourse_review")


def balanced_seeds(env_name: str, first_seed: int, clean: int, per_fault: int, rig: str = "nominal") -> list[int]:
    """Seeds giving `clean` fault-free cases and `per_fault` of every fault, in a shuffled order.

    A human set is small, so it is stratified to cover every fault. Its base rates therefore
    differ from the environment's own (four in ten clean)."""
    want: dict[str, int] = {}
    picked, seed = [], first_seed
    while True:
        env = make(env_name, seed=seed, rig=rig)
        if not want:
            want = {f: (clean if f == "none" else per_fault) for f in env.faults}
        if want[env.fault] > 0:
            want[env.fault] -= 1
            picked.append(seed)
        if not any(want.values()):
            break
        seed += 1
    order = np.random.default_rng([first_seed, 5]).permutation(len(picked))
    return [picked[i] for i in order]


def export_cases(env_name: str, n: int = 12, first_seed: int = 20000, rig: str = "nominal",
                 seeds: list[int] | None = None) -> dict:
    if env_name not in REVIEW_ENVS:
        raise ValueError(f"cases come from a review environment: {', '.join(REVIEW_ENVS)}")
    cases, spec = [], None
    for seed in (seeds if seeds is not None else range(first_seed, first_seed + n)):
        obs = make(env_name, seed=seed, rig=rig).reset()
        spec = obs.pop("spec")
        obs.pop("phase")
        cases.append({"id": task_id(env_name, rig, seed), **obs})
    return {"env": env_name, "task": spec["task"], "intended_protocol": spec["intended_protocol"],
            "causes": spec["actions"]["judge"]["cause"], "cases": cases}


def score_judgments(judgments: list[dict]) -> dict:
    """Each judgment: {"id": task id, "verdict": ..., "cause": ..., optional confidence, realistic, note}."""
    rows = []
    for j in judgments:
        env = make(j["id"])
        _, reward, done, info = env.step({"type": "judge", "verdict": j.get("verdict"), "cause": j.get("cause")})
        if not done:
            raise ValueError(f"{j['id']}: not a valid judgment: {j.get('verdict')!r}, {j.get('cause')!r}")
        name, _rig, seed = j["id"].split("/")
        rules_reward, rules_info = run_episode(make(j["id"]), POLICIES[name]["expert"](int(seed)))
        truth = info["truth"]
        rows.append({"id": j["id"], "env": name, "fault": truth["fault"], "should_repeat": truth["should_repeat"],
                     "verdict": j["verdict"], "cause": j["cause"], "reward": reward,
                     "verdict_ok": info["verdict_ok"], "cause_ok": info["cause_ok"],
                     "rules_reward": rules_reward, "rules_verdict_ok": rules_info["verdict_ok"],
                     "rules_cause_ok": rules_info["cause_ok"], "confidence": j.get("confidence"),
                     "realistic": j.get("realistic"), "note": (j.get("note") or "").strip()})
    return {"n": len(rows), "rows": rows, "by_env": _by(rows, "env"), "by_fault": _by(rows, "fault")}


def _by(rows: list[dict], key: str) -> dict:
    groups = defaultdict(list)
    for r in rows:
        groups[r[key]].append(r)

    def mean(items, field):
        return sum(float(i[field]) for i in items) / len(items)

    return {k: {"n": len(v), "reward": mean(v, "reward"), "verdict_ok": mean(v, "verdict_ok"),
                "cause_ok": mean(v, "cause_ok"), "rules_reward": mean(v, "rules_reward"),
                "rules_verdict_ok": mean(v, "rules_verdict_ok"), "rules_cause_ok": mean(v, "rules_cause_ok")}
            for k, v in sorted(groups.items())}


def markdown(score: dict, name: str = "expert") -> str:
    lines = [f"# Human review: {name}", "",
             f"{score['n']} cases judged blind, scored against planted truth. 'Rules' is the scripted QC "
             "reviewer on the same cases.", ""]
    for title, table in (("By environment", score["by_env"]), ("By planted fault", score["by_fault"])):
        lines += [f"## {title}", "", f"| | n | {name} reward | verdict right | cause right | rules reward | rules verdict | rules cause |",
                  "|---|---|---|---|---|---|---|---|"]
        for k, s in table.items():
            lines.append(f"| {k} | {s['n']} | {s['reward']:.2f} | {100 * s['verdict_ok']:.0f}% | {100 * s['cause_ok']:.0f}% | "
                         f"{s['rules_reward']:.2f} | {100 * s['rules_verdict_ok']:.0f}% | {100 * s['rules_cause_ok']:.0f}% |")
        lines.append("")
    unreal = [r for r in score["rows"] if r["realistic"] is False]
    noted = [r for r in score["rows"] if r["note"]]
    lines += ["## Flagged as not realistic (simulator to-dos)", ""]
    lines += [f"- {r['id']} (planted: {r['fault']}): {r['note'] or 'no note'}" for r in unreal] or ["None."]
    lines += ["", "## Notes in the expert's words (candidate rules)", ""]
    lines += [f"- {r['id']} (planted: {r['fault']}; called {r['verdict']}, {r['cause']}): {r['note']}" for r in noted] or ["None."]
    lines.append("")
    return "\n".join(lines)
