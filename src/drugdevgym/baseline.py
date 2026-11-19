"""Baseline: what a model scores on the environments, and which named mistakes it makes.

Two questions decide which mistakes are worth pricing on a real rig. How well does a current
model do on a task, and which of the rules does it actually break? A mistake no model makes
is not worth a plate. This module runs a model over seeded tasks (agent.run_episode) and
reports, per environment, the reward next to the scripted expert and null on the same
seeds, how the episodes ended, and for every rule-trap pair how often the model made that
mistake.

How a mistake is detected, and its limits:
  check       dose_response rules have a check on the design (knowledge.CHECKS).
  field       the model set the design field the trap sets to the trap's value. For a
              number or a list of time points this is an exact match, so it undercounts:
              a different wrong value is not seen.
  tool flag   the model's last call of the analysis tool chose the reading the trap uses.
              A model that does its own arithmetic is not seen.
  not detected  reporting conventions, campaign strategy and judgment calls.

No model has been run with this yet. The scripted policies are used to check the detector:
each scripted trap, run as if it were a model, has to be caught making its own mistake.
"""
from __future__ import annotations

import json

import numpy as np

from . import agent, campaign, knowledge
from .policies import GI_INTERACTIONS, POLICIES, SP_INTERACTIONS, run_episode
from .registry import make

FIRST_SEED = 300000  # away from the cards (0+), the review cases (10000+), rollouts (100000+) and the learner
PASS_AT = 0.5  # a stated line, not a validated pass rule: the share of episodes scoring at least this is reported
# A reading the trap uses, as the flags of the analysis tool that makes it.
READING = {
    ("stability", "no_lloq"): ("estimate_half_life", {"use_lloq": False}),
    ("gi_fluid_stability", "ignore_enzyme_decay"): ("estimate_half_life", {"correct_decay": False}),
    ("gi_camp_potency", "raw_ratio"): ("estimate_pec50", {"use_standards": False}),
    ("gi_camp_potency", "free_plateau"): ("estimate_pec50", {"anchor": False}),
    ("selectivity_panel", "raw_ratio"): ("selectivity_windows", {"use_standards": False}),
    ("cytotox_window", "nominal_above_solubility"): ("estimate_cc50", {"use_turbidity": False}),
    ("cytotox_window", "exact_not_bound"): ("estimate_cc50", {"bound_when_clean": False}),
}
FIELD_KINDS = ("data", "qc", "dose", "condition", "handling", "position", "upstream")


REFERENCE_SEEDS = (0, 1, 2)  # where each trap is compared with the expert to find what it sets
# A trap that is one mistake plus a second change is counted on the single mistake's own row.
INTERACTIONS = {("gi_fluid_stability", t): single for t, single in GI_INTERACTIONS.items()}
INTERACTIONS.update({("selectivity_panel", t): single for t, single in SP_INTERACTIONS.items()})


def _run_key(action: dict) -> str:
    """Which run an action is: its type, and the target it serves where an environment has several."""
    return str(action.get("type")) + (f":{action['target']}" if "target" in action else "")


def _values(actions: list[dict]) -> dict[tuple[str, str], set[str]]:
    """Every value each design field took, per run: (run, field) -> set of JSON values."""
    out: dict[tuple[str, str], set[str]] = {}
    for action in actions:
        for name, value in campaign._fields(action).items():
            if name not in ("target", "compounds"):
                out.setdefault((_run_key(action), name), set()).add(value)
    return out


def _is_number(value: str) -> bool:
    """A number or a list of numbers: matched exactly, so a different wrong value is not seen."""
    return not (value in ("true", "false", "null") or value.startswith('"'))


def _ran(env_name: str, policy: str, seed: int, traces: dict) -> list[dict]:
    key = (env_name, policy, seed)
    if key not in traces:
        traces[key] = [json.loads(a) for a in campaign.trace(env_name, policy, seed)["ran"]]
    return traces[key]


def trap_sets(env_name: str, trap: str, traces: dict) -> dict:
    """What a trap sets, found by comparing it with the expert on the reference seeds.

    fields: (run, field) -> the values the trap uses there, kept only where they are the same on
    every reference seed and differ from the expert's on at least one. A field that moves with the
    data (a cell number read off a titration) is a consequence, not what the trap sets, and is
    left out. dropped: runs the expert makes and the trap never does.
    """
    key = ("sets", env_name, trap)
    if key in traces:
        return traces[key]
    kind = campaign.TRAP_KIND[env_name][trap]
    per_seed, differs, dropped = [], set(), None
    for seed in REFERENCE_SEEDS:
        expert, theirs = _ran(env_name, "expert", seed, traces), _ran(env_name, trap, seed, traces)
        expert_values, trap_values = _values(expert), _values(theirs)
        per_seed.append(trap_values)
        differs |= {k for k, v in trap_values.items() if v - expert_values.get(k, set())}
        gone = {_run_key(a) for a in expert} - {_run_key(a) for a in theirs}
        dropped = gone if dropped is None else dropped & gone
    fields = {}
    for k in sorted(differs):
        field_kind = campaign.field_kind(k[1])
        mine = field_kind != "strategy" if kind == "upstream" else field_kind == kind
        if mine and all(other.get(k) == per_seed[0].get(k) for other in per_seed[1:]):
            fields[k] = per_seed[0][k]
    dropped = sorted(k for k in dropped if campaign.RUN_KIND.get(k.split(":")[0]) == kind)
    traces[key] = {"fields": fields, "dropped": dropped}
    return traces[key]


def mistakes(episode: dict, rules: list[dict] | None = None, traces: dict | None = None) -> dict[str, dict]:
    """For every rule-trap pair of the episode's environment: did the model make that mistake?

    Returns trap -> {"made": True, False or None, "how": ...}. None means it cannot be told from this
    episode: the model ran nothing, never set that field, never called that tool, the expert makes
    the same choice on this seed, or the pair is of a kind this harness does not detect.
    """
    env_name, seed = episode["env"], episode["seed"]
    rules = knowledge.load_rules() if rules is None else rules
    traces = {} if traces is None else traces
    runs = [campaign._canon(a["action"]) for a in episode["actions"]
            if a["accepted"] and isinstance(a["action"], dict) and a["action"].get("type") not in campaign.TERMINAL]
    mine, my_runs = _values(runs), {_run_key(a) for a in runs}
    out = {}
    for rule in rules:
        if rule["env"] != env_name:
            continue
        moot = knowledge.moot_traps(rule)
        for trap in rule["traps"]:
            if trap in moot:
                continue
            kind = campaign.TRAP_KIND[env_name][trap]
            if (env_name, trap) in READING:
                tool, flags = READING[(env_name, trap)]
                calls = [c for c in episode["tool_calls"] if c["tool"] == tool and c["ok"]]
                made = all(calls[-1]["args"].get(k) == v for k, v in flags.items()) if calls else None
                out[trap] = {"made": made, "how": "tool flag"}
            elif (env_name, trap) in INTERACTIONS:
                out[trap] = {"made": None, "how": f"counted as {INTERACTIONS[(env_name, trap)]}"}
            elif kind not in FIELD_KINDS:
                out[trap] = {"made": None, "how": "not detected"}
            elif not runs:
                out[trap] = {"made": None, "how": "ran nothing"}
            elif env_name == "dose_response" and "check" in rule:
                design = next((a.get("design") for a in runs if a.get("type") == "run_plate"), None)
                try:
                    out[trap] = {"made": not knowledge.CHECKS[rule["check"]](design), "how": "check"}
                except (KeyError, TypeError):
                    out[trap] = {"made": None, "how": "check"}
            else:
                sets = trap_sets(env_name, trap, traces)
                how = "field, exact value" if any(_is_number(v) for only in sets["fields"].values() for v in only) else "field"
                here = _values(_ran(env_name, "expert", seed, traces))
                # the fields the model set, leaving out any where the expert makes the trap's choice on this seed
                seen = [k for k in sets["fields"] if k in mine and not sets["fields"][k] <= here.get(k, set())]
                if not sets["fields"] and not sets["dropped"]:
                    out[trap] = {"made": None, "how": "follows the data"}
                elif not seen and not sets["dropped"]:
                    out[trap] = {"made": None, "how": how}  # nothing here tells the model's choice from the expert's
                else:
                    made = all(mine[k] & sets["fields"][k] for k in seen) and not any(k in my_runs for k in sets["dropped"])
                    out[trap] = {"made": bool(made), "how": how}
    return out


def run(make_model, envs: list[str], tasks: int = 20, attempts: int = 1, first_seed: int = FIRST_SEED,
        rig: str = "nominal", tools: bool = True, keep_transcripts: bool = True, log=None) -> list[dict]:
    """Run a model over tasks x attempts per environment. make_model(env_name, seed) returns the model."""
    episodes = []
    for env_name in envs:
        for seed in range(first_seed, first_seed + tasks):
            for attempt in range(attempts):
                episode = agent.run_episode(env_name, seed, make_model(env_name, seed), rig=rig, attempt=attempt, tools=tools)
                if not keep_transcripts:
                    episode.pop("transcript", None)
                episodes.append(episode)
                if log:
                    log(f"{episode['task']} attempt {attempt}: {episode['stop']}, reward {episode['reward']:.3f}")
    return episodes


def _reference(env_name: str, policy: str, seeds, rig: str) -> float | None:
    if policy not in POLICIES[env_name]:
        return None
    return float(np.mean([run_episode(make(env_name, seed=s, rig=rig), POLICIES[env_name][policy](s))[0] for s in seeds]))


def profile(episodes: list[dict], pass_at: float = PASS_AT) -> dict:
    """Per environment: reward, how episodes ended, cost in tokens, and the rate of every named mistake."""
    rules, traces, out = knowledge.load_rules(), {}, {}
    said = {(r["env"], t): r["said"] for r in rules for t in r["traps"]}
    for env_name in dict.fromkeys(e["env"] for e in episodes):
        mine = [e for e in episodes if e["env"] == env_name]
        rewards = np.array([e["reward"] for e in mine])
        seeds, rig = sorted({e["seed"] for e in mine}), mine[0]["rig"]
        by_task = [np.mean([e["reward"] for e in mine if e["seed"] == s]) for s in seeds]
        spread = [np.ptp([e["reward"] for e in mine if e["seed"] == s]) for s in seeds]
        gates = [e["info"]["qc_pass"] for e in mine if e["info"] and "qc_pass" in e["info"]]
        stops: dict[str, int] = {}
        for e in mine:
            stops[e["stop"]] = stops.get(e["stop"], 0) + 1
        found: dict[str, dict] = {}
        for e in mine:
            for trap, verdict in mistakes(e, rules, traces).items():
                cell = found.setdefault(trap, {"said": said[(env_name, trap)], "kind": campaign.TRAP_KIND[env_name][trap],
                                               "how": verdict["how"], "seen": 0, "made": 0})
                if verdict["made"] is not None:
                    cell["seen"] += 1
                    cell["made"] += bool(verdict["made"])
                    cell["how"] = verdict["how"]
        for cell in found.values():
            cell["rate"] = cell["made"] / cell["seen"] if cell["seen"] else None
        out[env_name] = {
            "episodes": len(mine), "tasks": len(seeds), "attempts": len(mine) // max(1, len(seeds)),
            "mean_reward": float(rewards.mean()), "task_sd": float(np.std(by_task, ddof=1)) if len(by_task) > 1 else None,
            "spread_within_task": float(np.mean(spread)), "pass_at": pass_at, "pass_rate": float(np.mean(rewards >= pass_at)),
            "expert_reward": _reference(env_name, "expert", seeds, rig), "null_reward": _reference(env_name, "null", seeds, rig),
            "qc_pass_rate": float(np.mean(gates)) if gates else None, "stops": stops,
            "finished_rate": float(np.mean([e["finished"] for e in mine])),
            "turns": float(np.mean([e["turns"] for e in mine])),
            "tool_calls": float(np.mean([len(e["tool_calls"]) for e in mine])),
            "invalid_replies": float(np.mean([e["invalid_replies"] for e in mine])),
            "rejected_actions": float(np.mean([e["rejected_actions"] for e in mine])),
            "input_tokens": float(np.mean([e["input_tokens"] for e in mine])),
            "output_tokens": float(np.mean([e["output_tokens"] for e in mine])),
            "mistakes": found,
        }
    return {"model": episodes[0]["model"] if episodes else None, "environments": out}


def _pct(x) -> str:
    return "n/a" if x is None else f"{100 * x:.0f}%"


def _f(x) -> str:
    return "n/a" if x is None else f"{x:.3f}"


def markdown(p: dict) -> str:
    lines = [f"# Baseline: {p['model']}", "",
             "Reward next to the scripted expert and the null policy on the same seeds. An episode that did not end "
             "by its final action scores 0. The pass line is a stated number, not a validated pass rule.", "",
             "| environment | tasks | attempts | mean reward | expert | null | at or above the pass line | finished | "
             "turns | tool calls | tokens in | tokens out |", "|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for env_name, e in p["environments"].items():
        lines.append(f"| {env_name} | {e['tasks']} | {e['attempts']} | {_f(e['mean_reward'])} | {_f(e['expert_reward'])} | "
                     f"{_f(e['null_reward'])} | {_pct(e['pass_rate'])} at {e['pass_at']:g} | {_pct(e['finished_rate'])} | "
                     f"{e['turns']:.1f} | {e['tool_calls']:.1f} | {e['input_tokens']:.0f} | {e['output_tokens']:.0f} |")
    lines += ["", "## How episodes ended", "", "| environment | ended as | episodes |", "|---|---|---|"]
    for env_name, e in p["environments"].items():
        lines += [f"| {env_name} | {stop} | {n} |" for stop, n in sorted(e["stops"].items())]
    lines += ["", "## Which named mistakes it made", "",
              "The rate is over the episodes where the mistake could be told. A field marked exact value counts only "
              "the trap's own number, so it undercounts.", "",
              "| environment | rule | mistake | kind | made it | could be told in | detected by |", "|---|---|---|---|---|---|---|"]
    for env_name, e in p["environments"].items():
        for trap, cell in e["mistakes"].items():
            rate = _pct(cell["rate"]) if cell["seen"] else "n/a"
            lines.append(f"| {env_name} | {cell['said']} | {trap} | {cell['kind']} | {rate} | {cell['seen']} of {e['episodes']} | "
                         f"{cell['how']} |")
    lines.append("")
    return "\n".join(lines)
