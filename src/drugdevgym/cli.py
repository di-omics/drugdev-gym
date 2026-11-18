"""Command line: drugdev-gym list | play | qc | learn | knowledge | cases | score."""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

from . import envqc, human, knowledge, learn
from .policies import DR_EXPERT, POLICIES, run_episode
from .registry import ENVS, RIGS, make


def _list(_args) -> None:
    for name, cls in ENVS.items():
        print(f"{name}: {sys.modules[cls.__module__].__doc__.strip().splitlines()[0]}")
        print(f"  policies: {', '.join(POLICIES[name])}")
    print(f"rigs: {', '.join(RIGS)}; task id = env/rig/seed, any non-negative integer seed")


def _positive(text: str) -> int:
    value = int(text)
    if value < 1:
        raise argparse.ArgumentTypeError("must be 1 or more")
    return value


def _plain(value):
    """JSON-safe copy: NaN and infinity become null."""
    if isinstance(value, dict):
        return {k: _plain(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_plain(v) for v in value]
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def _play(args) -> None:
    if args.policy not in POLICIES[args.env]:
        raise SystemExit(f"unknown policy for {args.env}: {args.policy}; have {', '.join(POLICIES[args.env])}")
    env = make(args.env, seed=args.seed, rig=args.rig)
    reward, info = run_episode(env, POLICIES[args.env][args.policy](args.seed))
    print(json.dumps({"task": f"{args.env}/{args.rig}/{args.seed}", "policy": args.policy,
                      "reward": round(reward, 4), "info": _plain(info)}, indent=2, default=float))


def _qc(args) -> None:
    out = Path(args.out)
    (out / "cards").mkdir(parents=True, exist_ok=True)
    summary = {}
    for name in (args.env or list(ENVS)):
        card = envqc.qc_card(name, n_seeds=args.seeds, rig=args.rig, first_seed=args.first_seed)
        suffix = "" if args.rig == "nominal" else f"_{args.rig}"
        (out / "cards" / f"{name}{suffix}.md").write_text(envqc.card_markdown(card))
        summary[name] = card
        print(envqc.card_markdown(card))
    suffix = "" if args.rig == "nominal" else f"_{args.rig}"
    path = out / f"envqc{suffix}.json"
    if args.env and path.exists():  # a partial run updates its own cards and keeps the others
        summary = {**json.loads(path.read_text()), **summary}
    path.write_text(json.dumps(summary, indent=2) + "\n")


def _learn(args) -> None:
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    runs = {rig: learn.cem(train_rig=rig, iterations=args.iterations, seed=args.seed) for rig in ("nominal", "randomized")}
    designs = {"learned_on_nominal": runs["nominal"]["design"], "learned_on_randomized": runs["randomized"]["design"],
               "expert": DR_EXPERT}
    table = learn.transfer_table(designs, n_seeds=args.test_seeds)
    (out / "learning.json").write_text(json.dumps({"runs": runs, "transfer": table}, indent=2) + "\n")
    lines = ["# Learning and transfer: dose_response", "",
             "Cross-entropy method over protocol choices, reward only. Test seeds are unseen by the learner.", "",
             "## Learning curve (population mean reward)", "", "| iteration | trained on nominal | trained on randomized |",
             "|---|---|---|"]
    for a, b in zip(runs["nominal"]["curve"], runs["randomized"]["curve"]):
        lines.append(f"| {a['iteration']} | {a['population_mean']:.3f} | {b['population_mean']:.3f} |")
    lines += ["", f"## Transfer (mean reward over {args.test_seeds} unseen seeds per rig)", "",
              "| design | nominal rig | randomized rigs | harsh rig (held out) |", "|---|---|---|---|"]
    for row in table:
        lines.append(f"| {row['design']} | {row['nominal']:.3f} | {row['randomized']:.3f} | {row['harsh']:.3f} |")
    lines += ["", "## Learned designs", ""]
    for name in ("learned_on_nominal", "learned_on_randomized"):
        lines.append(f"- {name}: `{json.dumps(designs[name])}`")
    lines.append("")
    (out / "learning.md").write_text("\n".join(lines))
    print("\n".join(lines))


def _knowledge(args) -> None:
    out = Path(args.out)
    cards = json.loads((out / "envqc.json").read_text())
    learned = None
    if (out / "learning.json").exists():
        runs = json.loads((out / "learning.json").read_text())["runs"]
        learned = {f"learned on {rig}": run["design"] for rig, run in runs.items()}
    rules = knowledge.load_rules()
    rows = knowledge.coverage(rules, cards, learned)
    (out / "knowledge.md").write_text(knowledge.markdown(rows))
    (out / "knowledge.json").write_text(json.dumps({"summary": knowledge.summary(rows), "rows": rows}, indent=2) + "\n")
    print(knowledge.markdown(rows))


def _cases(args) -> None:
    seeds = human.balanced_seeds(args.env, args.first_seed, args.clean, args.per_fault) if args.balanced else None
    data = human.export_cases(args.env, args.n, args.first_seed, seeds=seeds)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(data) + "\n")
    print(f"wrote {len(data['cases'])} {args.env} cases to {out}")
    if args.truth:
        Path(args.truth).write_text(json.dumps(human.frozen_truth(data), indent=1) + "\n")
        print(f"wrote the hidden answers to {args.truth}: keep this file away from the person judging")


def _score(args) -> None:
    judgments = json.loads(Path(args.judgments).read_text())
    if bool(args.cases) != bool(args.truth):
        raise SystemExit("give --cases and --truth together, or neither")
    if "/" in args.name or "\\" in args.name:
        raise SystemExit("--name must be a plain name")
    packs = [json.loads(Path(p).read_text()) for p in args.cases] if args.cases else None
    truth = None
    if args.truth:
        truth = {}
        for p in args.truth:
            loaded = json.loads(Path(p).read_text())
            truth.update(loaded.get("truth", loaded))
    score = human.score_judgments(judgments, packs=packs, truth=truth)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    (out / f"{args.name}.md").write_text(human.markdown(score, args.name))
    (out / f"{args.name}.json").write_text(json.dumps(_plain(score), indent=2) + "\n")
    print(human.markdown(score, args.name))


def main(argv=None) -> None:
    parser = argparse.ArgumentParser(prog="drugdev-gym", description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("list").set_defaults(func=_list)
    p = sub.add_parser("play")
    p.add_argument("env", choices=list(ENVS))
    p.add_argument("--policy", default="expert")
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--rig", default="nominal", choices=RIGS)
    p.set_defaults(func=_play)
    p = sub.add_parser("qc")
    p.add_argument("--env", action="append", choices=list(ENVS))
    p.add_argument("--seeds", type=int, default=200)
    p.add_argument("--first-seed", type=int, default=0)
    p.add_argument("--rig", default="nominal", choices=RIGS)
    p.add_argument("--out", default="results")
    p.set_defaults(func=_qc)
    p = sub.add_parser("cases", help="export review cases for a human expert (no truth inside)")
    p.add_argument("--env", required=True, choices=human.REVIEW_ENVS)
    p.add_argument("--n", type=_positive, default=12, help="number of consecutive seeds (ignored with --balanced)")
    p.add_argument("--first-seed", type=int, default=20000)
    p.add_argument("--balanced", action="store_true", help="cover every fault: --clean clean cases, --per-fault of each fault")
    p.add_argument("--clean", type=_positive, default=4)
    p.add_argument("--per-fault", type=_positive, default=2)
    p.add_argument("--out", required=True)
    p.add_argument("--truth", help="also write the hidden answers to this file, for frozen scoring later")
    p.set_defaults(func=_cases)
    p = sub.add_parser("score", help="score a human expert's judgments against planted truth")
    p.add_argument("judgments")
    p.add_argument("--cases", action="append", help="case pack the calls were made on (repeatable)")
    p.add_argument("--truth", action="append", help="hidden answers frozen at export (repeatable)")
    p.add_argument("--name", default="expert")
    p.add_argument("--out", default="results/human")
    p.set_defaults(func=_score)
    p = sub.add_parser("learn")
    p.add_argument("--iterations", type=int, default=12)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--test-seeds", type=int, default=100)
    p.add_argument("--out", default="results")
    p.set_defaults(func=_learn)
    p = sub.add_parser("knowledge")
    p.add_argument("--out", default="results")
    p.set_defaults(func=_knowledge)
    args = parser.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
