"""Tools: the analysis a model may call on data it has already been shown.

An environment returns raw wells and raw peak areas and scores a reported number. The fit in
between is real work, and a model with no way to compute would be scored on arithmetic. These
tools are that fit: the same functions the scripted experts use, called by name on the
observations of the episode. They read data only. Nothing here touches the simulator, the
hidden truth or the reward.

Every choice of reading is a required argument with both options described in the same
plain way, so that calling a tool does not make the choice for the model.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable

from . import analysis, judges
from .envs import gi_fluid_stability as gi
from .envs import selectivity_panel as sp
from .rollouts import _plain

USE_STANDARDS = ("bool", "true: convert each ratio to analyte through the plate's standard curve before fitting; "
                         "false: fit the raw ratio")
USE_LLOQ = ("bool", "true: drop samples below the limit of quantification before fitting; false: fit every sample")
RUN = ("int", "which data observation to use, 1 for the first; leave out for the latest")


@dataclass(frozen=True)
class Tool:
    name: str
    help: str
    fn: Callable
    args: dict = field(default_factory=dict)  # name -> (kind, help); every one is required except "run"
    needs: str = "wells"  # the observation key the tool reads

    def describe(self) -> dict:
        return {"name": self.name, "help": self.help,
                "args": {k: {"type": kind, "help": text, "required": k != "run"} for k, (kind, text) in {**self.args, "run": RUN}.items()}}


def _flat(obs: dict) -> dict:
    """lead_opt nests its two assays; the tools read them as wells and samples like any other."""
    if "dose_response" in obs and "stability" in obs:
        return {**obs, "wells": obs["dose_response"]["wells"], "samples": obs["stability"]["samples"]}
    return obs


def _pick(seen: list[dict], key: str, run) -> dict:
    with_data = [o for o in map(_flat, seen) if key in o]
    if not with_data:
        raise ValueError(f"no observation with {key} yet")
    if run is None:
        return with_data[-1]
    if isinstance(run, bool) or not isinstance(run, int) or not 1 <= run <= len(with_data):
        raise ValueError(f"run must be a whole number from 1 to {len(with_data)}")
    return with_data[run - 1]


def call_tool(toolset: list[Tool], name, args, seen: list[dict]) -> tuple[object, bool]:
    """Run one tool call. Returns (result, True) or (what was wrong, False). Never raises."""
    tool = next((t for t in toolset if t.name == name), None)
    if tool is None:
        return f"no such tool: {name}; have {', '.join(t.name for t in toolset) or 'none'}", False
    args = {} if args is None else args
    if not isinstance(args, dict):
        return "args must be an object", False
    extra = sorted(set(args) - set(tool.args) - {"run"})
    missing = sorted(set(tool.args) - set(args))
    if extra or missing:
        return (f"unknown argument: {extra[0]}" if extra else f"missing argument: {missing[0]}"), False
    for key, (kind, _text) in tool.args.items():
        value = args[key]
        if kind == "bool" and not isinstance(value, bool):
            return f"{key} must be true or false", False
        if kind.startswith("one of") and value not in kind[len("one of "):].split(", "):
            return f"{key} must be {kind}", False
    try:
        obs = _pick(seen, tool.needs, args.get("run"))
        return _plain(tool.fn(obs, seen, **{k: args[k] for k in tool.args})), True
    except Exception as err:  # bad data for this analysis: say so, do not end the episode
        return f"{name} failed: {err}", False


# ---------------------------------------------------------------- plate assays

def _pec50(obs, seen, use_standards):
    return analysis.estimate_pec50(obs["wells"], use_standards)


def _pec50_gi(obs, seen, use_standards, anchor):
    fit = analysis.estimate_pec50_anchored if anchor else analysis.estimate_pec50
    return fit(obs["wells"], use_standards)


def _titration(obs, seen):
    from .policies import gi_conditions  # the same arithmetic the scripted expert does on its titration plate
    forskolin, cells = gi_conditions({"forskolin_uM": "EC80", "cells_k": "titrate"}, obs["wells"])
    return {"curve": analysis.titration_curve(obs["wells"]), "forskolin_ec80_uM": forskolin, "cells_k_for_mid_curve": cells}


def _windows(obs, seen, use_standards):
    latest = {o["target"]: o["wells"] for o in seen if "target" in o and "wells" in o}  # the latest plate per target
    calls = {target: analysis.estimate_potency(wells, use_standards) for target, wells in latest.items()}
    out = {}
    for cid in sorted({c for by_cid in calls.values() for c in by_cid}):
        on = calls.get(sp.ON_TARGET, {}).get(cid, {}).get("pEC50")
        window = analysis.selectivity_window(on, [calls.get(t, {}).get(cid) for t in sp.OFF_TARGETS])
        if window is not None:
            out[cid] = window
    return {"windows": out, "targets_with_data": sorted(latest)}


CONTROL_ZPRIME = Tool("control_zprime", "Z-prime of the negative and positive control wells of one plate.",
                      lambda obs, seen: analysis.control_zprime(obs["wells"]))
ESTIMATE_PEC50 = Tool("estimate_pec50", "Fit one curve per compound on one plate and return each pEC50.", _pec50,
                      {"use_standards": USE_STANDARDS})
ESTIMATE_PEC50_GI = Tool("estimate_pec50", "Fit one inhibition curve per compound on one potency plate and return each pEC50.",
                         _pec50_gi, {"use_standards": USE_STANDARDS,
                                     "anchor": ("bool", "true: pin the uninhibited top of every curve on the forskolin-only "
                                                        "controls; false: leave the top free")})
TITRATION = Tool("titration_summary", "Fit the forskolin titration plate: the curve, the forskolin EC80, and the cell number "
                                      "that would put the forskolin control mid-way up the standard curve.", _titration)
ESTIMATE_POTENCY = Tool("estimate_potency", "Per compound on one plate: a fitted pEC50, or what the plate does and does not show.",
                        lambda obs, seen, use_standards: analysis.estimate_potency(obs["wells"], use_standards),
                        {"use_standards": USE_STANDARDS})
SELECTIVITY_WINDOWS = Tool("selectivity_windows", "The selectivity window of every compound from the latest plate of each target "
                                                  "run so far: exact where an off-target curve was measured, a lower bound "
                                                  "where nothing was seen up to the top dose.", _windows,
                           {"use_standards": USE_STANDARDS})
VIABILITY_ZPRIME = Tool("viability_zprime", "Z-prime of the vehicle and full-kill controls, per plate.",
                        lambda obs, seen: analysis.viability_zprime(obs["wells"]))
ESTIMATE_CC50 = Tool("estimate_cc50", "Fit a viability curve per compound and return each CC50 or bound.",
                     lambda obs, seen, use_turbidity, bound_when_clean: analysis.estimate_cc50(obs["wells"], use_turbidity,
                                                                                               bound_when_clean),
                     {"use_turbidity": ("bool", "true: treat doses that read turbid as not dissolved; false: take every "
                                                "nominal dose as exposure"),
                      "bound_when_clean": ("bool", "true: give a compound that shows no loss of viability as a lower bound at "
                                                   "its top dose; false: give that top dose as its CC50")})


# ---------------------------------------------------------------- time courses

def _half_life_gi(obs, seen, correct_decay, correct_activity, use_lloq):
    samples = obs["samples"]
    decay = gi.ACTIVITY_HL_37C_MIN if correct_decay else None
    activity = None
    if correct_activity:
        activity = analysis.enzyme_activity(samples, gi.CONTROL_ID, gi.CONTROL_HALF_LIFE_MIN, decay, use_lloq)
    compounds = [s for s in samples if s["id"] != gi.CONTROL_ID]
    return analysis.estimate_fluid_half_life(compounds, decay, activity, use_lloq, lo=gi.HL_LO, hi=gi.HL_HI)


ESTIMATE_HALF_LIFE = Tool("estimate_half_life", "Fit first-order decay per compound and return each half-life in minutes.",
                          lambda obs, seen, use_lloq: analysis.estimate_half_life(obs["samples"], use_lloq),
                          {"use_lloq": USE_LLOQ}, needs="samples")
ESTIMATE_HALF_LIFE_GI = Tool("estimate_half_life", "Fit decay per peptide and fluid and return each half-life in minutes.",
                             _half_life_gi,
                             {"correct_decay": ("bool", "true: fit against minutes of digestion at the starting enzyme activity; "
                                                        "false: fit against the clock"),
                              "correct_activity": ("bool", "true: scale to specification activity using the control substrate; "
                                                           "false: leave as measured"),
                              "use_lloq": USE_LLOQ}, needs="samples")


# ---------------------------------------------------------------- assay development and review

def _condition_zprime(obs, seen, metric):
    out = []
    for c in obs["conditions"]:
        z = (analysis.analyte_zprime(obs["standards"], c["pos"], c["neg"]) if metric == "analyte"
             else analysis.zprime(c["pos"], c["neg"]))
        out.append({"cells_k": c["cells_k"], "stim_min": c["stim_min"], "zprime": z})
    return out


def _target_cells(obs, seen):
    from .policies import WindowSearch
    return WindowSearch._target_cells(obs)


CONDITION_ZPRIME = Tool("condition_zprime", "Z-prime of every condition on one pilot plate.", _condition_zprime,
                        {"metric": ("one of analyte, raw", "analyte: after back-calculation through the standard curve; "
                                                           "raw: on the raw ratio")}, needs="conditions")
TARGET_CELLS = Tool("target_cells", "The cell number that would centre the two controls on the middle of the standard curve, "
                                    "read from one pilot plate, or null when the plate cannot say.", _target_cells,
                    needs="conditions")
PLATE_FEATURES = Tool("plate_features", "QC features of the plate under review: Z-prime, control window, drift across rows, "
                                        "disagreement between duplicates, outer-ring excess, and how far the reference "
                                        "compound reads from its known value.",
                      lambda obs, seen: judges.plate_features(obs["wells"], obs.get("reference")))
TIMECOURSE_FEATURES = Tool("timecourse_features", "QC features of the time course under review, including how far the "
                                                  "reference compound reads from its known value.",
                           lambda obs, seen: judges.timecourse_features(obs["samples"], obs.get("reference")), needs="samples")

TOOLS = {
    "dose_response": [CONTROL_ZPRIME, ESTIMATE_PEC50],
    "stability": [ESTIMATE_HALF_LIFE],
    "assay_window": [CONDITION_ZPRIME, TARGET_CELLS],
    "lead_opt": [CONTROL_ZPRIME, ESTIMATE_PEC50, ESTIMATE_HALF_LIFE],
    "plate_review": [PLATE_FEATURES, CONTROL_ZPRIME, ESTIMATE_PEC50],
    "timecourse_review": [TIMECOURSE_FEATURES, ESTIMATE_HALF_LIFE],
    "gi_fluid_stability": [ESTIMATE_HALF_LIFE_GI],
    "gi_camp_potency": [CONTROL_ZPRIME, TITRATION, ESTIMATE_PEC50_GI],
    "selectivity_panel": [CONTROL_ZPRIME, ESTIMATE_POTENCY, SELECTIVITY_WINDOWS],
    "cytotox_window": [VIABILITY_ZPRIME, ESTIMATE_CC50],
}
