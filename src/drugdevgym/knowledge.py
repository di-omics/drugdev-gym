"""Expert knowledge as something to hill-climb on.

The loop: an expert states a rule in their own words (rules.toml); the rule is written as a
trap, the expert protocol with that one rule broken (policies.py); environment QC measures
whether the reward punishes the trap (envqc.py); and a learner trained on reward alone shows
which rules it rediscovered (learn.py). The coverage table puts those side by side, so a
rule the reward cannot teach is visible as the next thing to build.
"""
from __future__ import annotations

import math
import tomllib
from pathlib import Path

from .policies import POLICIES

RULES_PATH = Path(__file__).with_name("rules.toml")
NOT_KNOWLEDGE = {"expert", "random", "null", "gate_bypass", "constant_repeat"}  # reference points and exploits, not broken rules
TEACHES = {"clean": "yes", "visible": "yes", "weak": "on average", "blind": "no", "inverted": "no, it rewards the mistake"}

# Does a dose-response design follow the rule? Used to ask what a learner rediscovered.
CHECKS = {
    "mix": lambda d: bool(d["mix"]),
    "new_tips": lambda d: bool(d["new_tips"]),
    "standards": lambda d: bool(d["standards"]),
    "avoid_edges": lambda d: bool(d["avoid_edges"]),
    "lid": lambda d: bool(d["lid"]),
    "controls": lambda d: d["controls_per_side"] >= 3,
    "span": lambda d: d["top_conc_nM"] >= 10000 and (d["n_points"] - 1) * math.log10(d["dilution_factor"]) >= 4.0,
    "ladder_volume": lambda d: d["transfer_ul"] * (d["dilution_factor"] - 1) >= 10.0 * d["replicates"] + 25.0,
}


def load_rules(path: Path = RULES_PATH) -> list[dict]:
    """Read and validate the rules: every env, trap and check they name must exist."""
    rules = tomllib.loads(Path(path).read_text())["rule"]
    seen = set()
    for r in rules:
        missing = {"id", "env", "said", "why", "traps", "status"} - set(r)
        if missing:
            raise ValueError(f"rule {r.get('id')}: missing {sorted(missing)}")
        if r["id"] in seen:
            raise ValueError(f"duplicate rule id: {r['id']}")
        seen.add(r["id"])
        if r["env"] not in POLICIES:
            raise ValueError(f"rule {r['id']}: unknown env {r['env']}")
        for trap in r["traps"]:
            if trap not in POLICIES[r["env"]] or trap in NOT_KNOWLEDGE:
                raise ValueError(f"rule {r['id']}: {trap} is not a trap of {r['env']}")
        if "check" in r and r["check"] not in CHECKS:
            raise ValueError(f"rule {r['id']}: unknown check {r['check']}")
        if r["status"] not in ("draft", "expert"):
            raise ValueError(f"rule {r['id']}: status must be draft or expert")
    return rules


def orphan_traps(rules: list[dict]) -> list[tuple[str, str]]:
    """Traps no rule explains. Every broken-rule policy should trace back to stated knowledge."""
    used = {(r["env"], t) for r in rules for t in r["traps"]}
    return [(env, t) for env, pols in POLICIES.items() for t in pols if t not in NOT_KNOWLEDGE and (env, t) not in used]


def coverage(rules: list[dict], cards: dict, learned: dict[str, dict] | None = None) -> list[dict]:
    """One row per (rule, trap): the QC verdict, and whether each learned design follows the rule.

    cards is the content of results/envqc.json. learned maps a label to a dose-response design.
    """
    rows = []
    for r in rules:
        by_policy = {row["policy"]: row for row in cards.get(r["env"], {}).get("rows", [])}
        for trap in r["traps"]:
            if trap not in by_policy:  # no card yet for this environment or trap: say so, do not guess
                rows.append({"env": r["env"], "id": r["id"], "said": r["said"], "status": r["status"], "trap": trap,
                             "verdict": None, "teaches": "not measured", "win_rate": None, "gap": None, "learner": None})
                continue
            q = by_policy[trap]
            row = {"env": r["env"], "id": r["id"], "said": r["said"], "status": r["status"], "trap": trap,
                   "verdict": q["verdict"], "teaches": TEACHES[q["verdict"]], "win_rate": q["win_rate"],
                   "gap": q["gap"], "learner": None}
            if learned and "check" in r:
                row["learner"] = {label: bool(CHECKS[r["check"]](design)) for label, design in learned.items()}
            rows.append(row)
    return rows


def summary(rows: list[dict]) -> dict:
    counts = {k: sum(r["teaches"] == v for r in rows) for k, v in
              (("yes", "yes"), ("on average", "on average"), ("no", "no"), ("inverted", TEACHES["inverted"]),
               ("not_measured", "not measured"))}
    checked = [r for r in rows if r["learner"]]
    labels = sorted({k for r in checked for k in r["learner"]})
    found = {label: sum(r["learner"][label] for r in checked) for label in labels}
    missed = sorted({r["id"] for r in checked if not all(r["learner"].values())})
    return {"rule_traps": len(rows), "rules": len({r["id"] for r in rows}), **counts,
            "checked": len(checked), "rediscovered": found, "missed": missed,
            "draft": len({r["id"] for r in rows if r["status"] == "draft"})}


def markdown(rows: list[dict]) -> str:
    s = summary(rows)
    lines = ["# Expert rules: what the rewards can teach", "",
             f"{s['rules']} rules, {s['rule_traps']} rule-trap pairs. The reward teaches {s['yes']} of them "
             f"(the expert wins on at least 90% of seeds), {s['on average']} only on average, is blind to {s['no']}, "
             f"and rewards the mistake in {s['inverted']}.",
             f"{s['draft']} of the {s['rules']} rules are drafts seeded from general bench practice and have not "
             "been reviewed by the expert yet.", "",
             "| env | rule | status | trap | reward teaches it | expert wins | reward gap |",
             "|---|---|---|---|---|---|---|"]
    for r in rows:
        if r["win_rate"] is None:
            lines.append(f"| {r['env']} | {r['said']} | {r['status']} | {r['trap']} | not measured | n/a | n/a |")
            continue
        lines.append(f"| {r['env']} | {r['said']} | {r['status']} | {r['trap']} | {r['teaches']} | "
                     f"{100 * r['win_rate']:.1f}% | {r['gap']:.3f} |")
    checked = [r for r in rows if r["learner"]]
    if checked:
        labels = sorted(checked[0]["learner"])
        lines += ["", "## What a learner rediscovered from reward alone", "",
                  "Cross-entropy learner on dose_response, one run per training rig (see learning.md). "
                  "A rule counts as rediscovered if the learner's final design follows it.", "",
                  "| rule | reward teaches it | " + " | ".join(labels) + " |", "|---|---|" + "---|" * len(labels)]
        for r in checked:
            marks = " | ".join("yes" if r["learner"][label] else "no" for label in labels)
            lines.append(f"| {r['said']} | {r['teaches']} | {marks} |")
        found = ", ".join(f"{label} {n} of {s['checked']}" for label, n in s["rediscovered"].items())
        lines += ["", f"Rediscovered: {found}."]
    lines.append("")
    return "\n".join(lines)
