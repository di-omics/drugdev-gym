"""Capture a human expert's judgment on review cases, and score it against planted truth.

export_cases writes what a reviewer may see (never the truth). The expert judges each case:
accept or repeat, the most likely cause, how sure they are, whether the data looks like
something a real bench would produce, and a note in their own words. score_judgments scores
those calls next to the scripted QC rules. The notes are the raw material for new rules;
"not realistic" flags are simulator to-dos.

A case id is the task id, so anyone holding the package can recompute a case's truth. For a
blind study give the expert opaque ids and keep the mapping private, or freeze the truth at
export (frozen_truth) and score against that file; frozen scoring also survives later changes
to the simulator.
"""
from __future__ import annotations

from collections import defaultdict

import numpy as np

from .envs.review import CAUSE_WEIGHT, VERDICT_WEIGHT
from .policies import POLICIES, run_episode
from .registry import make, task_id

REVIEW_ENVS = ("plate_review", "timecourse_review")
DATA_KEY = {"plate_review": "wells", "timecourse_review": "samples"}


def balanced_seeds(env_name: str, first_seed: int, clean: int, per_fault: int, rig: str = "nominal") -> list[int]:
    """Seeds giving `clean` fault-free cases and `per_fault` of every fault, in a shuffled order.

    A human set is small, so it is stratified to cover every fault. Its base rates therefore
    differ from the environment's own (four in ten clean)."""
    if clean < 0 or per_fault < 0 or clean + per_fault == 0:
        raise ValueError("clean and per_fault must be non-negative and not both zero")
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
    seeds = list(seeds) if seeds is not None else list(range(first_seed, first_seed + n))
    if not seeds:
        raise ValueError("no cases to export")
    cases, spec = [], None
    for seed in seeds:
        obs = make(env_name, seed=seed, rig=rig).reset()
        spec = obs.pop("spec")
        obs.pop("phase")
        cases.append({"id": task_id(env_name, rig, seed), **obs})
    return {"env": env_name, "task": spec["task"], "intended_protocol": spec["intended_protocol"],
            "causes": spec["actions"]["judge"]["cause"], "label": spec["label"], "scoring": spec["reward"],
            "cases": cases}


def frozen_truth(pack: dict) -> dict:
    """The hidden answers of a pack, to be stored away from the person judging it."""
    return {case["id"]: make(case["id"]).truth() for case in pack["cases"]}


def _canonical(case_id) -> tuple[str, str, int]:
    try:
        name, rig, seed = str(case_id).split("/")
        return name, rig, int(seed)
    except ValueError:
        raise ValueError(f"not a task id: {case_id!r}") from None


def _check(j: dict) -> None:
    for field in ("verdict", "cause"):
        if not isinstance(j.get(field), str):
            raise ValueError(f"{j.get('id')}: {field} must be text")
    if j.get("note") is not None and not isinstance(j["note"], str):
        raise ValueError(f"{j.get('id')}: note must be text")
    if j.get("realistic") is not None and not isinstance(j["realistic"], bool):
        raise ValueError(f"{j.get('id')}: realistic must be true, false or absent")
    if j.get("confidence") is not None and j["confidence"] not in (1, 2, 3):
        raise ValueError(f"{j.get('id')}: confidence must be 1, 2 or 3")


def score_judgments(judgments: list[dict], packs: list[dict] | None = None, truth: dict | None = None) -> dict:
    """Score calls of the form {"id", "verdict", "cause", optional confidence, realistic, note}.

    With packs and truth (as frozen at export) the case data and answers come from those files.
    Without them each case is regenerated from its id with the current simulator.
    """
    if (packs is None) != (truth is None):
        raise ValueError("give both packs and truth, or neither")
    case_data, causes = {}, {}
    for pack in packs or []:
        for case in pack["cases"]:
            case_data[case["id"]] = (pack["env"], case)
        causes[pack["env"]] = set(pack["causes"])
    rows, seen = [], set()
    for j in judgments:
        _check(j)
        name, rig, seed = _canonical(j["id"])
        cid = task_id(name, rig, seed)
        if cid in seen:
            raise ValueError(f"{cid}: judged more than once")
        seen.add(cid)
        if name not in REVIEW_ENVS:
            raise ValueError(f"{cid}: not a review case")
        if j["verdict"] not in ("accept", "repeat"):
            raise ValueError(f"{cid}: verdict must be accept or repeat")
        expert = POLICIES[name]["expert"](seed)
        if packs is not None:
            if cid not in case_data or cid not in truth:
                raise ValueError(f"{cid}: not in the case pack")
            if j["cause"] not in causes[name]:
                raise ValueError(f"{cid}: unknown cause {j['cause']!r}")
            _env, case = case_data[cid]
            answer = truth[cid]
            rules = expert.act({"phase": "review", DATA_KEY[name]: case[DATA_KEY[name]], "reference": case["reference"]})
        else:
            env = make(cid)
            if j["cause"] not in env.faults:
                raise ValueError(f"{cid}: unknown cause {j['cause']!r}")
            answer = env.truth()
            rules = expert.act(env.reset())
        verdict_ok = (j["verdict"] == "repeat") == answer["should_repeat"]
        cause_ok = j["cause"] == answer["fault"]
        rules_verdict_ok = (rules["verdict"] == "repeat") == answer["should_repeat"]
        rules_cause_ok = rules["cause"] == answer["fault"]
        rows.append({"id": cid, "env": name, "fault": answer["fault"], "should_repeat": answer["should_repeat"],
                     "verdict": j["verdict"], "cause": j["cause"],
                     "reward": VERDICT_WEIGHT * verdict_ok + CAUSE_WEIGHT * cause_ok,
                     "verdict_ok": verdict_ok, "cause_ok": cause_ok,
                     "rules_reward": VERDICT_WEIGHT * rules_verdict_ok + CAUSE_WEIGHT * rules_cause_ok,
                     "rules_verdict_ok": rules_verdict_ok, "rules_cause_ok": rules_cause_ok,
                     "confidence": j.get("confidence"), "realistic": j.get("realistic"),
                     "note": " ".join((j.get("note") or "").split())})
    for r in rows:
        r["group"] = f"{r['env']}: {r['fault']}"
    return {"n": len(rows), "of": len(case_data) if packs is not None else None, "rows": rows,
            "by_env": _by(rows, "env"), "by_fault": _by(rows, "group")}


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
    of = f" of {score['of']}" if score.get("of") else ""
    lines = [f"# Human review: {name}", "",
             f"{score['n']}{of} cases judged blind, scored against planted truth. 'Rules' is the scripted QC "
             "reviewer on the same cases.", ""]
    for title, table in (("By environment", score["by_env"]), ("By planted deviation", score["by_fault"])):
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
