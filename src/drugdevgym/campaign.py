"""Campaign: the real measurements that would price each bench mistake.

Environment QC says which mistakes the simulated reward can see. A lab in the loop has to
answer the next question: what real data would measure what each one costs? That depends
on the kind of choice the mistake changes.

From complete records of correct runs
  reading     same liquid, read differently. A complete raw record prices it, read again.
  data        data left off the run: standards, time points, a plate, the internal
              standard, the time zero. A record designed for it prices it, with the same
              data left out.
From grids
  dose        which concentrations were tested. A dense grid where every well is dosed on
              its own holds every sparser design. On a serial ladder a design is its own
              ladder.
  condition   a setting of the assay: cells, time, stimulant, exposure. A grid of
              conditions prices it, repeated over days and cell batches, because the right
              setting moves.
On the rig
  handling    how the liquid was moved or the sample treated. Measure the physical
              quantity it changes directly on the rig, then confirm with a few same-day
              paired runs.
  position    where on the plate the wells sit. Uniformity plates and a split layout.
Only when something is wrong
  qc          run QC left off: controls, replicates, an activity control. A correct run
              has nothing for them to catch, so their value shows on runs that carry a
              fault or a liability, and depends on how often faults happen.
Not a plate question
  convention  what may be claimed where there is no data, such as a bound past the top
              dose. A rule, with no price on a plate.
  upstream    an assay mistake made inside a campaign. Priced by the same mistake in its
              own assay environment.
  strategy    which experiments a campaign does.
  judgment    a call on finished data. Priced by review cases with a planted answer.

The kind of every design field is written down once (FIELD_KIND) and the kind of every trap
once (TRAP_KIND). The simulator checks the second against the first (check), by four rules.
A reading or a convention asks the rig for exactly what the expert asks for. A judgment is
made where nothing is run at all. An upstream trap changes the same assay field, to the
same value, as the assay trap it names. Any other trap changes a design field of its own
kind, or drops a run of that kind, and every other field it changes is named and explained
(ALSO). plan() refuses to build while any check fails.

That check is about what a trap changes in the simulator. It cannot tell a reading from a
convention, since both run what the expert runs; that split is the author's call. And
whether the named data would really measure the cost on a bench is a judgment no check
here can make. The kinds are drafts until the expert has read them. NEEDS holds what an
automated review, with one reviewer prompted to read as a bench scientist, said each pair
depends on. Nothing here is measured on hardware.
"""
from __future__ import annotations

import itertools
import json
import math
from statistics import NormalDist

from . import knowledge
from .policies import POLICIES
from .registry import ENVS, make

GROUPS = {
    "records": ("reading", "data"),
    "grids": ("dose", "condition"),
    "rig": ("handling", "position"),
    "faults": ("qc",),
    "not_a_plate": ("convention", "upstream", "strategy", "judgment"),
}
GROUP_TITLE = {"records": "from complete records of correct runs", "grids": "from grids", "rig": "on the rig",
               "faults": "only when something is wrong", "not_a_plate": "not a plate question"}
KINDS = tuple(kind for kinds in GROUPS.values() for kind in kinds)
GROUP_OF = {kind: group for group, kinds in GROUPS.items() for kind in kinds}
PRICED_BY = {
    "reading": "a complete raw record of a correct run, read again",
    "data": "a record of a correct run designed for it, with the same data left out",
    "dose": "a dense grid of concentrations where every well is dosed on its own; on a serial ladder, a second ladder",
    "condition": "a grid of assay conditions, repeated over days and cell batches",
    "handling": "the physical quantity it changes, measured directly on the rig; then a few same-day paired runs",
    "position": "uniformity plates and a split layout on one plate",
    "qc": "the same item left out of runs that carry a fault or a liability, and how often faults happen",
    "convention": "nothing on a plate: a rule about what may be claimed",
    "upstream": "the size of the error from its own assay environment; what it costs a campaign, only from a campaign",
    "strategy": "a whole campaign of experiments, not a plate",
    "judgment": "review cases with a planted answer, not a new run",
}
SAME_RUNS = ("reading", "convention")  # the trap asks the rig for exactly what the expert asks for

TERMINAL = {"report", "nominate", "judge"}  # actions that ask the rig for nothing
MAX_PAIRED = 7  # most unmatched runs of one type that _pair_up will pair
CONTAINERS = ("design", "dose_response", "stability")  # where an action nests a design

# What kind of choice each design field is. Every field an agent can set in an action that
# runs something must be listed (a test checks it against every environment's spec), so a
# new environment has to say what its choices are.
FIELD_KIND = {
    # which concentrations
    "top_conc_nM": "dose", "top_conc_uM": "dose", "dilution_factor": "dose", "n_points": "dose",
    # data on the run
    "standards": "data", "timepoints_min": "data", "internal_standard": "data", "prequenched_t0": "data",
    "fluids": "data", "target": "data", "compounds": "data",
    # run QC
    "replicates": "qc", "controls_per_side": "qc", "activity_control": "qc", "vehicle_controls": "qc",
    "kill_controls": "qc", "interference_control": "qc",
    # settings of the assay
    "cells_k": "condition", "stim_min": "condition", "conditions": "condition", "exposure_h": "condition",
    "forskolin_uM": "condition", "ibmx": "condition",
    # how the liquid is moved and the sample treated
    "transfer_ul": "handling", "mix": "handling", "new_tips": "handling", "lid": "handling",
    "mix_after_spike": "handling", "match_dmso": "handling", "fresh_enzyme": "handling",
    "hold_min": "handling", "hold_temp_c": "handling", "sgf_quench": "handling", "sif_quench": "handling",
    # where on the plate
    "avoid_edges": "position",
    # which experiments to do at all
    "analogs": "strategy",
}
# A whole run that only one side made: what kind of choice dropping or adding it is.
RUN_KIND = {"run_plate": "data", "run": "data", "run_pilot": "condition", "run_titration": "condition",
            "commit": "condition", "make_test": "strategy"}

# What kind of mistake each trap is. Drafts, checked against the simulator by check().
TRAP_KIND = {
    "dose_response": {
        "no_mix": "handling", "reuse_tips": "handling", "thin_ladder": "handling", "no_lid": "handling",
        "edge_wells": "position",
        "raw_ratio": "data",  # the standard curve is left off the plate and out of the fit
        "few_controls": "qc",
        "narrow_range": "dose",
    },
    "stability": {
        "no_mix": "handling", "no_lid": "handling", "no_lid_on_edge_wells": "handling",
        "edge_wells": "position",
        "no_internal_standard": "data", "two_point": "data", "early_only": "data",
        "single_replicate": "qc",
        "no_lloq": "reading",
    },
    "assay_window": {
        "max_cells": "condition", "min_cells": "condition", "kit_default": "condition",
        "raw_ratio_zprime": "condition",  # the wrong reading of the pilots shows up as a different condition chosen
        "no_lid": "handling",
        "edge_wells": "position",
    },
    "lead_opt": {
        "upstream_no_mix": "upstream", "upstream_raw_ratio": "upstream", "upstream_two_point": "upstream",
        "upstream_no_internal_standard": "upstream",
        "no_model": "strategy", "one_round": "strategy",
    },
    "plate_review": {"always_accept": "judgment", "always_repeat": "judgment"},
    "timecourse_review": {"always_accept": "judgment", "always_repeat": "judgment"},
    "gi_fluid_stability": {
        "acid_quench_in_sgf": "handling", "acid_quench_in_sgf_warm_overnight": "handling", "stale_enzyme": "handling",
        "no_true_time_zero": "data", "late_first_point": "data", "no_internal_standard": "data",
        "no_activity_control": "qc",
        "ignore_enzyme_decay": "reading",
    },
    "gi_camp_potency": {
        "forskolin_too_high": "condition", "skip_titration": "condition", "no_forskolin": "condition",
        "no_ibmx": "condition", "too_many_cells": "condition",
        "no_mix": "handling",
        "raw_ratio": "reading", "free_plateau": "reading",
    },
    "selectivity_panel": {
        "on_target_range": "dose", "single_point": "dose", "single_point_as_curve": "dose",
        "no_standards": "data",
        "on_target_only": "data",  # the off-target plates are left out
        "no_controls": "qc",
        "raw_ratio": "reading",
        "claim_past_range": "convention", "inactive_as_number": "convention",
    },
    "cytotox_window": {
        "unmatched_dmso": "handling",
        "short_exposure": "condition",
        "edge_wells": "position",
        "no_kill_control": "qc", "no_interference_control": "qc",
        "nominal_above_solubility": "reading",
        "exact_not_bound": "convention",
    },
}
# A field a trap changes that is not of the trap's own kind, and why it moves. Where a trap changes
# fields of two kinds, the entry here is what records which one the author took as the mistake.
ALSO = {
    ("dose_response", "narrow_range"): {"transfer_ul": "keeps the wells of the 2-fold ladder above the dead volume"},
    ("stability", "no_lid_on_edge_wells"): {"avoid_edges": "the mistake is the 2 together"},
}
# An upstream trap is this assay mistake, made inside a campaign.
UPSTREAM = {
    ("lead_opt", "upstream_no_mix"): ("dose_response", "no_mix"),
    ("lead_opt", "upstream_raw_ratio"): ("dose_response", "raw_ratio"),
    ("lead_opt", "upstream_two_point"): ("stability", "two_point"),
    ("lead_opt", "upstream_no_internal_standard"): ("stability", "no_internal_standard"),
}
# What each pair's stated data depends on, or the direct measurement that prices it. Drafts, from an
# automated review with one reviewer prompted to read as a bench scientist. A pair priced from
# records with no entry here is taken to hold on a plain complete record.
_TIME_COURSE = "one well per time point or a large incubation, and the real quench time of every sample"
_PLANTED_FAULT = ("a run with a planted cell, reader or reagent fault; a ladder mistake does not touch the control "
                  "wells, so the handling arms cannot show it")
_DYE_MIX = "a dye ladder with and without mixing"
_UNIFORMITY = "uniformity plates and a split inner and outer layout"
_SECOND_LADDER = "on a serial ladder this is a second ladder, run beside the first"
_TOP_WELL = "it is the top well of the expert's own ladder, so a plain record of the expert's plate holds it today"
NEEDS = {
    # from records
    ("stability", "two_point"): _TIME_COURSE,
    ("stability", "early_only"): _TIME_COURSE,
    ("gi_fluid_stability", "late_first_point"): _TIME_COURSE,
    ("stability", "no_internal_standard"): "a fixed injection order: the price holds for that order only",
    ("gi_fluid_stability", "no_internal_standard"): "a fixed injection order: the price holds for that order only",
    ("gi_fluid_stability", "no_true_time_zero"): "both zeros on the run: the pre-quenched one and the earliest sample "
                                                 "after the spike, with its real quench time",
    ("stability", "no_lloq"): "blanks and low calibrators, so the limit of quantification is measured",
    ("cytotox_window", "nominal_above_solubility"): "a turbidity read or a solubility measurement in the same medium",
    ("gi_fluid_stability", "ignore_enzyme_decay"): "the activity of each fluid lot followed over the incubation: a "
                                                   "series of its own, since the control substrate is gone too soon",
    # run QC
    ("dose_response", "few_controls"): _PLANTED_FAULT,
    ("selectivity_panel", "no_controls"): _PLANTED_FAULT,
    ("cytotox_window", "no_kill_control"): _PLANTED_FAULT,
    ("stability", "single_replicate"): "runs with an unmixed spike, to see whether one replicate still flags them",
    ("gi_fluid_stability", "no_activity_control"): "an aged or off-specification enzyme lot among the runs",
    ("cytotox_window", "no_interference_control"): "a reference compound that interferes with the readout",
    # grids
    ("dose_response", "narrow_range"): _SECOND_LADDER,
    ("selectivity_panel", "on_target_range"): _SECOND_LADDER,
    ("selectivity_panel", "single_point"): _TOP_WELL,
    ("selectivity_panel", "single_point_as_curve"): _TOP_WELL + "; it also misreports the hit, which no rig field shows",
    ("assay_window", "raw_ratio_zprime"): "what the trap changes is how the pilots are read; the condition it commits "
                                          "follows from that, and on some seeds it is the expert's",
    ("cytotox_window", "short_exposure"): "staggered dosing or a readout that does not end the plate, and a slow-acting "
                                          "reference compound",
    # on the rig: the direct measurement
    ("dose_response", "no_mix"): _DYE_MIX,
    ("gi_camp_potency", "no_mix"): _DYE_MIX,
    ("dose_response", "thin_ladder"): "a dye ladder at the small and the large transfer volume",
    ("dose_response", "reuse_tips"): "the ladder wells of a sticky reference compound read by LC-MS, with fresh and "
                                     "reused tips; a dye does not stick to a tip",
    ("stability", "no_mix"): "time-zero samples alone, mixed and unmixed",
    ("stability", "no_lid_on_edge_wells"): "plate weights or a stable marker in the matrix; and the expert's answer "
                                           "first, since the lid alone is ruled moot",
    ("gi_fluid_stability", "acid_quench_in_sgf"): "each gastric sample split into the acid and the base quench",
    ("gi_fluid_stability", "acid_quench_in_sgf_warm_overnight"): "the same split, injected again after the hold, with "
                                                                 "the temperature of the held plate",
    ("gi_fluid_stability", "stale_enzyme"): "a control substrate in fresh and aged fluid, the age and storage stated",
    ("cytotox_window", "unmatched_dmso"): "a DMSO tolerance curve on the cells; here the mistake is punished by a "
                                          "plate-map rule, not by a simulated effect",
    ("dose_response", "edge_wells"): _UNIFORMITY,
    ("assay_window", "edge_wells"): _UNIFORMITY,
    ("cytotox_window", "edge_wells"): _UNIFORMITY,
    ("stability", "edge_wells"): "a stable marker compound in inner and outer wells of the incubation plate; there is "
                                 "no plate read here, so the uniformity study does not cover it",
}
# Environments whose reported result is not a log10 value; the paired-run arithmetic is borrowed for them.
NOT_A_LOG_VALUE = {"assay_window": "a Z-prime"}


def _canon(value):
    """A plain copy that compares by value: numbers as floats, tuples and arrays as lists."""
    if isinstance(value, dict):
        return {str(k): _canon(v) for k, v in value.items()}
    if hasattr(value, "tolist") and not isinstance(value, (str, bytes)):
        value = value.tolist()  # numpy scalar or array
    if isinstance(value, (list, tuple)):
        return [_canon(v) for v in value]
    if isinstance(value, bool) or value is None or isinstance(value, str):
        return value
    if isinstance(value, (int, float)):
        return float(value)
    return str(value)


def measured(step_log: list[dict]) -> dict:
    """What a step log measured: plate wells read, and samples quenched for LC-MS."""
    read = {(s["run"], s["dst_plate"]) for s in step_log if s.get("op") == "read"}
    wells = {(s["run"], s["dst_plate"], s["dst_well"]) for s in step_log
             if "dst_well" in s and (s["run"], s["dst_plate"]) in read}
    samples = {(s["run"], s["dst_plate"], s["dst_well"]) for s in step_log
               if "dst_well" in s and str(s.get("dst_plate", "")).endswith("quench")}
    return {"wells": len(wells), "samples": len(samples)}


def trace(env_name: str, policy_name: str, seed: int, rig="nominal") -> dict:
    """One episode: the distinct runs the policy got the rig to do, and what one pass of them measures.

    ran is a sorted list of canonical JSON strings, one per distinct action that made the rig do
    something. An action the environment rejects is not a run. A repeated plate of the same design
    appears once and is measured once: repeating is not a different choice. rig_runs is the number
    of runs as executed, repeats included; two policies that differ only in how often they repeat
    one design compare as the same.
    """
    env = make(env_name, seed=seed, rig=rig)
    policy = POLICIES[env_name][policy_name](seed)
    obs, ran, kept, rig_runs = env.reset(), [], [], 0
    while True:
        action = policy.act(obs)
        before = len(env.step_log)
        obs, _reward, done, _info = env.step(action)
        if action.get("type") not in TERMINAL and len(env.step_log) > before:
            rig_runs += 1
            key = json.dumps(_canon(action), sort_keys=True)
            if key not in ran:
                ran.append(key)
                kept += env.step_log[before:]
        if done:
            break
    return {"ran": sorted(ran), "rig_runs": rig_runs, **measured(kept)}


def _fields(action: dict) -> dict[str, str]:
    """An action as field -> value. A nested design is opened up; "design." is dropped from the name."""
    out = {}
    for key, value in action.items():
        if key == "type":
            continue
        if key in CONTAINERS and isinstance(value, dict):
            prefix = "" if key == "design" else key + "."
            out.update({prefix + k: json.dumps(v, sort_keys=True) for k, v in value.items()})
        else:
            out[key] = json.dumps(value, sort_keys=True)
    return out


def _pair_up(mine: list[dict], theirs: list[dict]) -> tuple[set[str], int]:
    """Pair runs of one type so that the fewest fields differ. Returns (fields that differ, runs left unpaired)."""
    short, long = (mine, theirs) if len(mine) <= len(theirs) else (theirs, mine)
    if not short:
        return set(), len(long)
    flat_short, flat_long = [_fields(a) for a in short], [_fields(a) for a in long]

    def diff(a: dict, b: dict) -> set[str]:
        return {k for k in set(a) | set(b) if a.get(k) != b.get(k)}

    if len(long) > MAX_PAIRED:  # every assignment is tried, so keep the search small; no real episode comes close
        raise ValueError(f"{len(long)} unmatched {long[0]['type']} runs on one side: too many to pair")
    best = None
    for picks in itertools.permutations(range(len(long)), len(short)):
        diffs = [diff(flat_short[i], flat_long[j]) for i, j in enumerate(picks)]
        cost = (sum(len(d) for d in diffs), sorted(set().union(*diffs)))
        if best is None or cost < best[0]:
            best = (cost, set().union(*diffs))
    return best[1], len(long) - len(short)


def _differences(expert_ran: list[str], trap_ran: list[str]) -> set[str]:
    """Fields in which two sets of runs differ. A run only one side made is "+type" or "-type"."""
    only_expert = [json.loads(s) for s in expert_ran if s not in set(trap_ran)]
    only_trap = [json.loads(s) for s in trap_ran if s not in set(expert_ran)]
    fields: set[str] = set()
    for kind in sorted({a["type"] for a in only_expert + only_trap}):
        mine = [a for a in only_trap if a["type"] == kind]
        theirs = [a for a in only_expert if a["type"] == kind]
        differ, left = _pair_up(mine, theirs)
        fields |= differ
        if left:
            fields.add(("+" if len(mine) > len(theirs) else "-") + kind)
    return fields


def compare(env_name: str, trap: str, seeds=range(3), traces: dict | None = None) -> dict:
    """What the trap's runs change against the expert's on the nominal rig, over the same seeds.

    fields is the union over seeds. differs_on counts the seeds on which anything differed: a trap
    whose choice depends on the data can be the expert on some seeds.
    """
    seeds = list(seeds)
    if not seeds:
        raise ValueError("give at least one seed")
    traces = {} if traces is None else traces

    def get(policy: str, seed: int) -> dict:
        key = (env_name, policy, seed)
        if key not in traces:
            traces[key] = trace(env_name, policy, seed)
        return traces[key]

    fields: set[str] = set()
    differs_on = 0
    for seed in seeds:
        found = _differences(get("expert", seed)["ran"], get(trap, seed)["ran"])
        differs_on += bool(found)
        fields |= found
    return {"fields": sorted(fields), "differs_on": differs_on, "seeds": len(seeds),
            "expert_runs": any(get("expert", seed)["ran"] for seed in seeds)}


def changed(env_name: str, trap: str, seeds=range(3), traces: dict | None = None) -> list[str]:
    """Fields in which the trap's runs differ from the expert's. Empty: it runs what the expert runs."""
    return compare(env_name, trap, seeds, traces)["fields"]


def field_kind(field: str) -> str:
    """Kind of one changed field, or of a whole run that was dropped or added."""
    if field[0] in "+-":
        if field[1:] not in RUN_KIND:
            raise KeyError(f"no RUN_KIND for a dropped or added {field[1:]}")
        return RUN_KIND[field[1:]]
    leaf = field.rsplit(".", 1)[-1]
    if leaf not in FIELD_KIND:
        raise KeyError(f"no FIELD_KIND for the design field {leaf}")
    return FIELD_KIND[leaf]


def _values(env_name: str, policy: str, seeds, traces: dict) -> dict[str, set[str]]:
    """Every value each design field took in a policy's runs, over seeds (field -> set of JSON values)."""
    out: dict[str, set[str]] = {}
    for seed in seeds:
        for action in traces[(env_name, policy, seed)]["ran"]:
            for field, value in _fields(json.loads(action)).items():
                out.setdefault(field, set()).add(value)
    return out


def check(env_name: str, trap: str, seeds=range(3), traces: dict | None = None) -> str | None:
    """None when the stated kind of a trap agrees with what the simulator shows it changes, else what is wrong.

    reading, convention: the trap runs exactly what the expert runs (the check cannot tell the two
    apart). judgment: nothing is run at all. upstream: the assay fields it changes sit in the design
    of the environment its named assay trap belongs to, are the fields that trap changes, and take
    the same values; anything else it changes is strategy. Any other kind: at least one changed
    field, or dropped run, is of that kind, and every other changed field is named in ALSO with
    its reason. A note in ALSO that explains nothing is itself an error.
    """
    if trap not in TRAP_KIND.get(env_name, {}):
        return "no kind stated"
    kind = TRAP_KIND[env_name][trap]
    seeds = list(seeds)
    traces = {} if traces is None else traces
    seen = compare(env_name, trap, seeds, traces)
    fields, also = seen["fields"], ALSO.get((env_name, trap), {})
    kinds = {f: field_kind(f) for f in fields}
    if any(not str(reason).strip() for reason in also.values()):
        return "ALSO gives no reason for " + ", ".join(sorted(f for f, r in also.items() if not str(r).strip()))
    used = kind not in SAME_RUNS + ("judgment", "upstream")  # the kinds whose check reads ALSO
    stale = sorted(f for f in also if not used or kinds.get(f, kind) == kind)
    if stale:
        return f"ALSO explains {', '.join(stale)}, which this trap does not change as a field of another kind"
    if kind == "judgment":
        return None if not seen["expert_runs"] and not fields else "called judgment, but the environment runs something"
    if not seen["expert_runs"]:
        return f"called {kind}, but nothing is run in this environment: that is judgment"
    if kind in SAME_RUNS:
        return None if not fields else f"called {kind}, but its runs differ in {', '.join(fields)}"
    if not fields:
        return (f"called {kind}, but on seeds {', '.join(str(s) for s in seeds)} it runs exactly what the expert runs; "
                "if its runs follow from the data, check more seeds")
    if kind == "upstream":
        source = UPSTREAM.get((env_name, trap))
        if source is None:
            return "called upstream, but no assay trap is named for it in UPSTREAM"
        source_kind = TRAP_KIND.get(source[0], {}).get(source[1])
        if source_kind in SAME_RUNS + ("judgment", "upstream", None):
            return f"called upstream of {'/'.join(source)}, which is {source_kind}: name an assay trap that changes a run"
        assay = sorted(f for f, k in kinds.items() if k != "strategy")
        if not assay:
            return "called upstream, but it changes no assay field: that is strategy"
        if any(not f.startswith(source[0] + ".") for f in assay):
            return f"called upstream of {'/'.join(source)}, but it changes {', '.join(assay)}, outside that assay's design"
        theirs = changed(*source, seeds, traces)
        if [f.split(".", 1)[1] for f in assay] != theirs:
            return f"called upstream of {'/'.join(source)}, but it changes {', '.join(assay)} and that trap changes {', '.join(theirs)}"
        mine, source_values = _values(env_name, trap, seeds, traces), _values(source[0], source[1], seeds, traces)
        differ = sorted(f for f in assay if mine[f] != source_values[f.split(".", 1)[1]])
        if differ:
            return f"called upstream of {'/'.join(source)}, but it sets {', '.join(differ)} to another value than that trap"
        return None
    if kind not in kinds.values():
        return f"called {kind}, but it changes {', '.join(fields)} ({', '.join(sorted(set(kinds.values())))})"
    unexplained = sorted(f for f, k in kinds.items() if k != kind and f not in also)
    if unexplained:
        return f"called {kind}, but it also changes {', '.join(unexplained)}, which ALSO does not explain"
    return None


def pairs_per_arm(sd_within: float = 0.15, delta: float = 0.25, alpha: float = 0.05, power: float = 0.8) -> int:
    """Same-day pairs (the run done right and the run done wrong, side by side) to see a mean shift of delta.

    Paired t-test, by the normal approximation with Guenther's correction:
    n = ((z(1 - alpha / 2) + z(power)) * sd_d / delta) ** 2 + z(1 - alpha / 2) ** 2 / 2, rounded up, where
    sd_d = sqrt(2) * sd_within is the scatter of the difference between the two runs of a pair.
    Inputs are in log10 units of the reported value, and both defaults are assumptions: 0.15 for
    the scatter of one run within a day (in dose_response as simulated the median for one expert
    potency is about 0.13), 0.25 for half the 0.5 reward tolerance. It sizes a shift in the mean
    only: a mistake that mostly widens the scatter, or fails a run now and then, needs more. The
    approximation can fall 1 pair short of the exact test, mostly below about 7 pairs.
    """
    values = (sd_within, delta, alpha, power)
    if not all(isinstance(v, (int, float)) and not isinstance(v, bool) and math.isfinite(v) for v in values):
        raise ValueError("sd, delta, alpha and power must be finite numbers")
    if not (sd_within > 0 and delta > 0 and 0 < alpha < 1 and alpha / 2 < power < 1 and 1e-6 <= sd_within / delta <= 1e3):
        raise ValueError("sd and delta must be positive and within a factor of a million of each other, alpha between "
                         "0 and 1, power between alpha / 2 and 1")
    z_alpha = NormalDist().inv_cdf(1 - alpha / 2)
    z = z_alpha + NormalDist().inv_cdf(power)
    return max(3, math.ceil((z * math.sqrt(2) * sd_within / delta) ** 2 + z_alpha ** 2 / 2 - 1e-9))


def design_fields() -> dict[str, set[str]]:
    """Every design field an agent can set in an action that runs something, per environment."""
    out = {}
    for name in ENVS:
        fields: set[str] = set()
        for action_type, spec in make(name, seed=0).spec()["actions"].items():
            if action_type in TERMINAL:
                continue
            for field, value in spec.items():
                fields |= set(value) if field in CONTAINERS and isinstance(value, dict) else {field}
        out[name] = fields
    return out


def plan(cards: dict, seeds=range(3), sd_within: float = 0.15, delta: float = 0.25, alpha: float = 0.05,
         power: float = 0.8, rules: list[dict] | None = None) -> dict:
    """The campaign: every rule-trap pair with its kind, what a record must keep, and the arms on the rig.

    cards is the content of results/envqc.json (for what the simulated reward teaches today). Raises
    ValueError when a stated kind disagrees with what the simulator shows the trap changes.
    """
    n = pairs_per_arm(sd_within, delta, alpha, power)  # first: a bad number should fail before any episode runs
    rules = knowledge.load_rules() if rules is None else rules
    seeds = list(seeds)
    traces: dict = {}

    def size(env_name: str, policy: str) -> dict:
        runs = [traces[(env_name, policy, s)] for s in seeds]
        return {key: sum(r[key] for r in runs) / len(runs) for key in ("wells", "samples", "rig_runs")}

    rows, wrong = [], []
    for c in knowledge.coverage(rules, cards):
        env_name, trap = c["env"], c["trap"]
        if trap not in TRAP_KIND.get(env_name, {}):
            raise KeyError(f"no TRAP_KIND for {env_name}/{trap}: say what kind of choice this mistake changes")
        kind = TRAP_KIND[env_name][trap]
        message = check(env_name, trap, seeds, traces)
        if message:
            wrong.append(f"{env_name}/{trap}: {message}")
        seen = compare(env_name, trap, seeds, traces)
        for f in seen["fields"]:
            field_kind(f)  # every changed field has to have a kind
        source = UPSTREAM.get((env_name, trap))
        rows.append({"env": env_name, "id": c["id"], "said": c["said"], "trap": trap, "moot": c["moot"], "kind": kind,
                     "group": GROUP_OF[kind], "priced_by": PRICED_BY[kind], "needs": NEEDS.get((env_name, trap)),
                     "also": ALSO.get((env_name, trap)), "upstream_of": "/".join(source) if source else None,
                     "teaches": c["teaches"], "win_rate": c["win_rate"], "changes": seen["fields"],
                     "differs_on": seen["differs_on"], **size(env_name, trap)})
    if wrong:
        raise ValueError("stated kinds disagree with the simulator: " + "; ".join(wrong))
    live = [r for r in rows if not r["moot"]]

    verdicts = ("yes", "on average", "no", knowledge.TEACHES["inverted"], "not measured")
    by_kind = {k: {"pairs": 0, **{v: 0 for v in verdicts}} for k in KINDS}
    for r in live:
        by_kind[r["kind"]]["pairs"] += 1
        by_kind[r["kind"]][r["teaches"]] += 1
    by_group = {g: sum(by_kind[k]["pairs"] for k in kinds) for g, kinds in GROUPS.items()}
    env_order = list(dict.fromkeys(r["env"] for r in rows))
    by_env = {e: {g: sum(r["env"] == e and r["group"] == g for r in live) for g in GROUPS} for e in env_order}

    records = [r for r in live if r["group"] == "records"]
    times = {}
    for env_name in env_order:
        traps = [r["trap"] for r in records if r["env"] == env_name and "timepoints_min" in r["changes"]]
        if not traps:
            continue
        seen_times: set[float] = set()
        for policy in ["expert"] + traps:
            for s in seeds:
                for action in traces[(env_name, policy, s)]["ran"]:
                    seen_times |= set(json.loads(action).get("design", {}).get("timepoints_min", []))
        times[env_name] = sorted(seen_times)

    grids = {}
    for r in live:
        if r["group"] == "grids":
            for f in r["changes"]:
                if f[0] not in "+-" and field_kind(f) == r["kind"]:
                    grids.setdefault(r["env"], {}).setdefault(f, []).append(r["trap"])

    arms = []
    for env_name in env_order:
        mine = [r for r in live if r["env"] == env_name and r["kind"] == "handling"]
        if not mine:
            continue
        expert = size(env_name, "expert")
        arms.append({"env": env_name, "arm": "expert", "kind": "reference", "teaches": None,
                     "wells": expert["wells"], "samples": expert["samples"]})
        arms += [{"env": r["env"], "arm": r["trap"], "kind": r["kind"], "teaches": r["teaches"],
                  "wells": r["wells"], "samples": r["samples"]} for r in mine]
    totals = {"arms": sum(a["kind"] != "reference" for a in arms), "reference_arms": sum(a["kind"] == "reference" for a in arms),
              "pairs_per_arm": n, "runs": len(arms) * n,
              "wells": round(sum(a["wells"] for a in arms) * n), "samples": round(sum(a["samples"] for a in arms) * n)}

    # what no trap exercises: the kind of these is stated but nothing checks it
    exercised = {(r["env"], f) for r in rows for f in r["changes"] if f[0] not in "+-"}
    unexercised = {e: sorted(f for f in fields if (e, f) not in exercised) for e, fields in design_fields().items()}
    unexercised = {e: fields for e, fields in unexercised.items() if fields and e != "lead_opt"}  # lead_opt nests the assay designs
    dropped = {f[1:] for r in rows for f in r["changes"] if f[0] in "+-"}
    unexercised_runs = sorted(set(RUN_KIND) - dropped)

    plain = sum(1 for r in records if not r["needs"])
    top_well = sum(1 for r in live if r["kind"] == "dose" and (r["needs"] or "").startswith(_TOP_WELL))
    return {"rows": rows, "by_kind": by_kind, "by_group": by_group, "by_env": by_env, "times": times, "grids": grids,
            "arms": arms, "totals": totals, "unexercised": unexercised, "unexercised_runs": unexercised_runs,
            "seeds": seeds, "assumptions": {"sd_within": sd_within, "delta": delta, "alpha": alpha, "power": power},
            "counts": {"pairs": len(live), "moot": len(rows) - len(live), **by_group,
                       "records_plain": plain, "records_conditional": len(records) - plain, "dose_top_well": top_well}}


def _num(x: float) -> str:
    return f"{x:.0f}" if abs(x - round(x)) < 0.05 else f"{x:.1f}"


def markdown(p: dict) -> str:
    c, t, a, k = p["counts"], p["totals"], p["assumptions"], p["by_kind"]
    live = [r for r in p["rows"] if not r["moot"]]
    inverted = knowledge.TEACHES["inverted"]
    lines = [
        "# Campaign: what would price each mistake on a real rig", "",
        f"{c['pairs']} rule-trap pairs count. Sorted by the kind of choice the mistake changes:", "",
        f"- {c['records']} can be priced from complete records of correct runs. {k['reading']['pairs']} read the same "
        f"raw data differently and {k['data']['pairs']} leave data off the run. Of the {c['records']}, "
        f"{c['records_plain']} hold on a plain complete record and {c['records_conditional']} need the record to "
        "keep something specific, listed below.",
        f"- {c['grids']} are a point on a grid. {k['condition']['pairs']} are a setting of the assay, which moves "
        f"between days and cell batches, so the grid has to be repeated. {k['dose']['pairs']} are a choice of "
        f"concentrations: {c['dose_top_well']} are the top well of the expert's own ladder and hold on a plain "
        f"record today, and {k['dose']['pairs'] - c['dose_top_well']} need a second ladder until every well is "
        "dosed on its own.",
        f"- {c['rig']} change something physical on the rig: {k['handling']['pairs']} how the liquid or the sample is "
        f"handled, {k['position']['pairs']} where the wells sit. Direct measurement comes first, paired runs second.",
        f"- {c['faults']} leave out run QC. A correct run has nothing for a control to catch, so these show their "
        "value only on a run that carries the fault the control exists for.",
        f"- {c['not_a_plate']} are not plate questions: {k['convention']['pairs']} reporting conventions, "
        f"{k['upstream']['pairs']} assay mistakes made inside a campaign, {k['strategy']['pairs']} campaign strategy, "
        f"{k['judgment']['pairs']} judgment on finished data.", "",
        "The kinds are drafts the expert has not read. The simulator checks what each trap changes, by 4 rules. A "
        "reading or a convention asks the rig for exactly what the expert asks for; the check cannot tell those 2 "
        "apart, so that split is the author's call. A judgment is made where nothing is run. An upstream trap "
        "changes the same assay field, to the same value, as the assay trap it names. Any other trap changes a "
        "design field of its own kind, or drops a run of that kind, and every other field it changes is explained "
        "at the end. That is a check on the simulator's actions. Whether the named data would measure the cost on a "
        "bench is a judgment. The note on each pair comes from an automated review, with one reviewer prompted to "
        "read as a bench scientist. Nothing here is measured on hardware.", "",
        "## By kind", "",
        "| group | kind | what prices it | pairs | reward teaches it today | only on average | blind |",
        "|---|---|---|---|---|---|---|"]
    for kind in KINDS:
        cell = k[kind]
        lines.append(f"| {GROUP_TITLE[GROUP_OF[kind]]} | {kind} | {PRICED_BY[kind]} | {cell['pairs']} | {cell['yes']} | "
                     f"{cell['on average']} | {cell['no']} |")
    odd = {v: sum(k[kind][v] for kind in KINDS) for v in (inverted, "not measured")}
    if any(odd.values()):
        lines += ["", f"Not in the table: {odd[inverted]} pairs where the reward prefers the mistake and "
                  f"{odd['not measured']} with no QC card yet."]
    lines += ["", "## By environment", "",
              "| environment | " + " | ".join(GROUP_TITLE[g] for g in GROUPS) + " |", "|---|" + "---|" * len(GROUPS)]
    for env_name, counts in p["by_env"].items():
        lines.append(f"| {env_name} | " + " | ".join(str(counts[g]) for g in GROUPS) + " |")

    def table(group: str, last_column: str) -> list[str]:
        out = ["| environment | trap | kind | what it changes | " + last_column + " |", "|---|---|---|---|---|"]
        for r in live:
            if r["group"] != group:
                continue
            what = ", ".join(r["changes"]) or "nothing the rig is asked for"
            if r["upstream_of"]:
                last = f"the same mistake as {r['upstream_of']}: {r['priced_by']}"
            elif group == "not_a_plate":
                last = r["priced_by"]
            else:
                last = r["needs"] or "nothing more"
            out.append(f"| {r['env']} | {r['trap']} | {r['kind']} | {what} | {last} |")
        return out

    lines += ["", "## From complete records of correct runs", "",
              "A record prices one of these only if it kept, as raw values, what the mistake leaves out or reads "
              "differently.", ""] + table("records", "the record also needs")
    for env_name, points in p["times"].items():
        lines += ["", f"A {env_name} record that holds every time-point design here samples at "
                  f"{', '.join(_num(x) for x in points)} min. That is more samples than one incubation well gives."]
    lines += ["", "## From grids", "",
              "| environment | grid over | prices |", "|---|---|---|"]
    for env_name, fields in p["grids"].items():
        traps = sorted({trap for names in fields.values() for trap in names})
        lines.append(f"| {env_name} | {', '.join(fields)} | {', '.join(traps)} |")
    lines += ["", "A setting chosen for one day's cells is the answer for that day. The price of not titrating is "
              "how far the best setting moves between days and cell batches, so a condition grid is repeated, and "
              "a setting chosen from a pilot can move a second setting with it. A dose grid is a lookup only where "
              "every well is dosed on its own. Every environment here doses by serial ladder."]
    noted = [r for r in live if r["group"] == "grids" and r["needs"]]
    if noted:
        lines += [""] + [f"- {r['env']}/{r['trap']}: {r['needs']}." for r in noted]
    lines += ["", "## On the rig", "",
              "Handling: measure the physical quantity the mistake changes, then confirm with paired runs, the run "
              "done right and the run done wrong on the same day. Position is an effect inside a plate.", "",
              "| environment | trap | kind | reward teaches it today | direct measurement first | plate wells per run | "
              "LC-MS samples per run |", "|---|---|---|---|---|---|---|"]
    for r in live:
        if r["group"] == "rig":
            lines.append(f"| {r['env']} | {r['trap']} | {r['kind']} | {r['teaches']} | {r['needs'] or 'not stated'} | "
                         f"{_num(r['wells'])} | {_num(r['samples'])} |")
    seeds = ", ".join(str(s) for s in p["seeds"])
    lines += ["", f"If each of the {t['arms']} handling mistakes were also confirmed by paired runs, with the expert "
              f"run beside it in each of the {t['reference_arms']} environments, the assumed numbers give about "
              f"{t['pairs_per_arm']} same-day pairs per arm: {t['runs']} runs, {t['wells']} plate wells and "
              f"{t['samples']} LC-MS samples. The {t['pairs_per_arm']} is for a shift of {a['delta']:g} log units in "
              f"the mean when one run scatters by {a['sd_within']:g} within a day (two-sided {a['alpha']:g}, power "
              f"{a['power']:g}, a paired t-test by a normal approximation that can fall 1 pair short). Both numbers "
              "are assumptions, so this is a starting size and not a bound. A smaller shift, a mistake that mostly "
              "widens the scatter or fails a run now and then, or an arm that often returns no fitted value needs "
              f"more, and {t['arms']} comparisons read at {a['alpha']:g} each will show a false difference somewhere "
              f"more often than that. Run sizes are one pass of each simulated design in the 96-well format (seeds "
              f"{seeds})."]
    odd_units = sorted({arm["env"] for arm in p["arms"]} & set(NOT_A_LOG_VALUE))
    lines += [f"{e} reports {NOT_A_LOG_VALUE[e]}, not a log value: the pair count is borrowed for it." for e in odd_units]
    lines += ["", "## Only when something is wrong", ""] + table("faults", "needs")
    lines += ["", "## Not a plate question", ""] + table("not_a_plate", "priced by")
    lines += ["", "## Every pair", "",
              "| environment | trap | kind | what it changes | differs from the expert on | reward teaches it today |",
              "|---|---|---|---|---|---|"]
    for r in live:
        what = ", ".join(r["changes"]) or "nothing the rig is asked for"
        often = f"{r['differs_on']} of {len(p['seeds'])} seeds" if r["changes"] else "no seed"
        lines.append(f"| {r['env']} | {r['trap']} | {r['kind']} | {what} | {often} | {r['teaches']} |")
    explained = [r for r in p["rows"] if r["also"]]
    if explained:
        lines += ["", "## Fields a trap changes outside its own kind", "",
                  "Where a trap changes fields of 2 kinds, this is the author's choice of which one is the mistake.", ""]
        lines += [f"- {r['env']}/{r['trap']} ({r['kind']}) also changes {field}: {reason}."
                  for r in explained for field, reason in r["also"].items()]
    if p["unexercised"] or p["unexercised_runs"]:
        lines += ["", "## Kinds nothing checks yet", "",
                  "No trap changes these design fields, so their kind is stated and unchecked. lead_opt is left out: "
                  "it nests the 2 assay designs.", ""]
        lines += [f"- {env_name}: {', '.join(fields)}" for env_name, fields in p["unexercised"].items()]
        if p["unexercised_runs"]:
            lines.append(f"- whole runs no trap drops or adds: {', '.join(p['unexercised_runs'])}")
    moot = [r for r in p["rows"] if r["moot"]]
    if moot:
        lines += ["", "## Ruled moot by the expert", "", "Left out of every count above.", "",
                  "| environment | trap | kind |", "|---|---|---|"]
        lines += [f"| {r['env']} | {r['trap']} | {r['kind']} |" for r in moot]
    lines.append("")
    return "\n".join(lines)
