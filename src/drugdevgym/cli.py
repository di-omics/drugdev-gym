"""Command line: drugdev-gym list | play | qc | learn."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from . import envqc, learn
from .policies import DR_EXPERT, POLICIES, run_episode
from .registry import ENVS, RIGS, make


def _list(_args) -> None:
    for name, cls in ENVS.items():
        print(f"{name}: {sys.modules[cls.__module__].__doc__.strip().splitlines()[0]}")
        print(f"  policies: {', '.join(POLICIES[name])}")
    print(f"rigs: {', '.join(RIGS)}; task id = env/rig/seed, any integer seed")


def _play(args) -> None:
    env = make(args.env, seed=args.seed, rig=args.rig)
    reward, info = run_episode(env, POLICIES[args.env][args.policy](args.seed))
    print(json.dumps({"task": f"{args.env}/{args.rig}/{args.seed}", "policy": args.policy,
                      "reward": round(reward, 4), "info": info}, indent=2, default=float))


def _qc(args) -> None:
    out = Path(args.out)
    (out / "cards").mkdir(parents=True, exist_ok=True)
    summary = {}
    for name in (args.env or list(ENVS)):
        card = envqc.qc_card(name, n_seeds=args.seeds, rig=args.rig)
        suffix = "" if args.rig == "nominal" else f"_{args.rig}"
        (out / "cards" / f"{name}{suffix}.md").write_text(envqc.card_markdown(card))
        summary[name] = card
        print(envqc.card_markdown(card))
    suffix = "" if args.rig == "nominal" else f"_{args.rig}"
    (out / f"envqc{suffix}.json").write_text(json.dumps(summary, indent=2) + "\n")


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
    p.add_argument("--rig", default="nominal", choices=RIGS)
    p.add_argument("--out", default="results")
    p.set_defaults(func=_qc)
    p = sub.add_parser("learn")
    p.add_argument("--iterations", type=int, default=12)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--test-seeds", type=int, default=100)
    p.add_argument("--out", default="results")
    p.set_defaults(func=_learn)
    args = parser.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
