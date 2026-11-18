"""Selectivity panel: spend a plate budget on five cell assays, then report each selectivity window or its honest bound.

There is one on-target and four related off-targets. Every plate carries one cell line (one
target) and is run with the dose_response plate physics: a serial dilution on the rig, a
Gq-style cell response, a competitive TR-FRET readout. The agent decides which target each
plate serves, which compounds go on it, single-point screen or full curve, the dose range, and
the usual protocol choices. It then reports, per compound, the selectivity window in log
units: on-target pEC50 minus the pEC50 at its most potent off-target.

Off-target activity is usually far weaker than on-target activity, so often nothing is seen up
to the highest concentration tested. The honest report is then a lower bound ("at least X").
The reward pays a bound that is true and sits at the labware limit, pays an exact window for
its accuracy, and penalises a window claimed beyond what the tested range supports or more than
a log unit beyond what the compound has. A report made of such claims scores below no report.

What the plate budget does and does not decide: with four compounds it never binds. All four fit
on one 10-point curve plate per target, so one curve plate per target (five plates) is the only
sensible plan and the sixth plate is a spare for a repeat. The budget only forbids an uncovered
panel (two wasted plates leave an off-target untested) and a full screen-then-curve plan (nine
plates); it never forces a choice between compounds or targets. The decisions that carry the
reward are the off-target dose range, curve or single point, the controls, and how each window
is reported.
"""
from __future__ import annotations

import copy
from dataclasses import replace

import numpy as np

from ..analysis import zprime
from ..rig import ROWS, RigProfile, usable_wells
from .base import COST_WEIGHT, Env, spec_help, validate
from .dose_response import (CELLS_UL, DESIGN_SPEC, LADDER_MAX_UL, STAMP_UL, STD_SERIES_NM, STIM_MIN, ZPRIME_GATE,
                            CellAssay, Compound, draw_assay, draw_compound, plate_cost, simulate_plate)

ON_TARGET = "on"
OFF_TARGETS = ("off1", "off2", "off3", "off4")
TARGETS = (ON_TARGET,) + OFF_TARGETS

PLATE_BUDGET = 6  # plates per episode: one per target plus one spare for a repeat or a follow-up
COST_REFERENCE = 100.0 * PLATE_BUDGET  # 100 per budgeted plate, the same proportion as dose_response (200 for 2 plates)
ON_PEC50_RANGE = (7.0, 9.0)  # hidden on-target pEC50: leads between 100 nM and 1 nM
MARGIN_RANGE = (3.0, 8.0)  # log units an ordinary off-target is weaker than the on-target (weak or inactive)
LIABILITY_RATE = 0.3  # chance that a compound carries one close off-target liability
LIABILITY_MARGIN_RANGE = (0.0, 2.0)  # log units the liability off-target is weaker than the on-target
MAX_TOP_NM = DESIGN_SPEC["top_conc_nM"][2]  # highest final concentration the labware allows (100 uM)
P_FLOOR = 9.0 - float(np.log10(MAX_TOP_NM))  # weakest pEC50 the labware can show: half response at the top dose
TOLERANCE = 0.5  # log units of window error at which credit reaches zero (the repository's reward tolerance)
SUPPORT_MARGIN = 0.5  # log units a claimed window may exceed the supported window before it counts as unsupported
FALSE_CLAIM_MARGIN = 1.0  # log units a claimed window may exceed the true window before it counts as a false claim
CLAIM_PENALTY = 0.5  # score of an unsupported or false claim is minus this: worse than no answer
MIN_CONTROLS = 3  # controls a side a plate needs before its Z' exists (analysis.zprime)
WINDOW_RANGE = (-12.0, 12.0)  # accepted range of a reported window, wider than any possible truth

PLATE_SPEC = {**DESIGN_SPEC,
              "n_points": ("int", 1, 12, "concentrations per compound; 1 is a single-point screen at top_conc_nM")}
REPORT_SPEC = {
    "window": ("float", WINDOW_RANGE[0], WINDOW_RANGE[1], "selectivity window in log units"),
    "at_least": ("bool", None, None, "true: a lower bound (the window is at least this); false: an exact window"),
}


def draw_cell_lines(rng: np.random.Generator) -> dict[str, CellAssay]:
    """One cell line per target, each with its own window and noise, all read with one kit lot."""
    first = draw_assay(rng)
    lines = {ON_TARGET: first}
    for target in OFF_TARGETS:
        other = draw_assay(rng)
        scale = first.kit.k / other.kit.k
        lines[target] = CellAssay(first.kit, basal=other.basal * scale, top=other.top * scale, cell_cv=other.cell_cv)
    return lines


def draw_pharmacology(rng: np.random.Generator, cid: str) -> dict[str, Compound]:
    """One compound at every target: its own potency, efficacy and slope per target, one stickiness.

    Ordinary off-targets sit MARGIN_RANGE log units below the on-target potency. With probability
    LIABILITY_RATE one off-target, chosen at random, sits only LIABILITY_MARGIN_RANGE below it.
    """
    on = replace(draw_compound(rng, cid), pec50=float(rng.uniform(*ON_PEC50_RANGE)))
    margins = rng.uniform(*MARGIN_RANGE, size=len(OFF_TARGETS))
    liable = rng.random() < LIABILITY_RATE
    which = int(rng.integers(len(OFF_TARGETS)))
    close = float(rng.uniform(*LIABILITY_MARGIN_RANGE))
    if liable:
        margins[which] = close
    out = {ON_TARGET: on}
    for target, margin in zip(OFF_TARGETS, margins):
        out[target] = replace(draw_compound(rng, cid), pec50=on.pec50 - float(margin), sticky=on.sticky)
    return out


def tested_limit(top_conc_nM: float) -> float:
    """Weakest pEC50 a test up to top_conc_nM can show: 9 - log10(top concentration in nM)."""
    return 9.0 - float(np.log10(top_conc_nM))


def true_window(on_pec50: float, off_pec50s) -> float:
    """On-target pEC50 minus the pEC50 at the most potent off-target."""
    return float(on_pec50 - max(off_pec50s))


def supported_window(on_pec50: float, tested_tops, on_tested: bool = True) -> float:
    """Widest window the tested range can support for one compound.

    tested_tops holds, per off-target, the top concentration (nM) of the plate that pair rests on
    (see reward_text), or None if the compound was never tested there. One untested off-target
    supports nothing, and neither does an untested on-target: a window needs both of its ends.
    """
    if not on_tested or any(t is None for t in tested_tops):
        return float("-inf")
    return float(on_pec50 - max(tested_limit(t) for t in tested_tops))


def score_claim(window: float, at_least: bool, true: float, achievable: float, supported: float) -> float:
    """Score one reported window. See reward_text() for the rule in words."""
    if window > supported + SUPPORT_MARGIN or window > true + FALSE_CLAIM_MARGIN:
        return -CLAIM_PENALTY
    miss = max(abs(window - achievable), window - true) if at_least else abs(window - true)
    return float(np.clip(1.0 - miss / TOLERANCE, 0.0, 1.0))


def reward_text() -> str:
    return (
        "Per compound: W = true on-target pEC50 minus the highest true off-target pEC50 (the true window). "
        f"A = true on-target pEC50 minus {P_FLOOR:g}, the widest window the labware can show (nothing seen up to "
        f"{MAX_TOP_NM:g} nM). "
        f"Plate QC: a plate passes when it has at least {MIN_CONTROLS} controls a side and control Z' >= "
        f"{ZPRIME_GATE:g} on the raw ratio. "
        "Every compound and target pair that was tested rests on one plate: among the plates that tested the pair, "
        "those that pass plate QC if any do (otherwise all of them), the one with the highest top_conc_nM (the "
        "earliest if tied). "
        f"S = true on-target pEC50 minus the highest, over the {len(OFF_TARGETS)} off-targets, of "
        "9 - log10(top_conc_nM of the plate the pair rests on); if the compound was never tested on the on-target "
        "or on any one of the off-targets, no window is supported (S is minus infinity). "
        f"Exact claim {{window: x, at_least: false}}: max(0, 1 - |x - W| / {TOLERANCE:g}). "
        f"Lower bound {{window: x, at_least: true}}: max(0, 1 - max(|x - A|, x - W) / {TOLERANCE:g}), so a bound earns "
        "full credit when it is true and sits at the labware limit, and loses credit when it is looser than that or "
        "false. "
        f"Penalised claim, either kind: if x > S + {SUPPORT_MARGIN:g} (unsupported by the tested range) or "
        f"x > W + {FALSE_CLAIM_MARGIN:g} (more selectivity than the compound has, by over a log unit), the score is "
        f"-{CLAIM_PENALTY:g} instead. "
        "No answer scores 0. Quality = mean over compounds; it can be negative. "
        f"Reward = quality, halved when it is positive and the QC gate fails, minus {COST_WEIGHT:g} * cost / "
        f"{COST_REFERENCE:g}. "
        "QC gate: at least one plate was run and, for every compound and target pair that was tested, the plate "
        "the pair rests on passes plate QC."
    )


def check_plate(design, n_compounds: int) -> tuple[dict | None, str | None]:
    """Validate one plate design for n_compounds compounds. Returns (clean, error)."""
    clean, err = validate(PLATE_SPEC, design)
    if err:
        return None, err
    if clean["transfer_ul"] * clean["dilution_factor"] > LADDER_MAX_UL:
        return None, f"ladder well would hold more than {LADDER_MAX_UL:.0f} uL; lower transfer_ul or dilution_factor"
    need = n_compounds * clean["n_points"] * clean["replicates"] + 2 * clean["controls_per_side"]
    need += len(STD_SERIES_NM) if clean["standards"] else 0
    have = len(usable_wells(clean["avoid_edges"]))
    if need > have:
        return None, f"design needs {need} wells but the plate has {have}"
    return clean, None


class SelectivityPanel(Env):
    name = "selectivity_panel"
    salted = True
    cost_reference = COST_REFERENCE
    MAX_PLATES = PLATE_BUDGET

    def __init__(self, seed: int = 0, rig="nominal", n_compounds: int = 4):
        super().__init__(seed, rig)
        if not 1 <= n_compounds <= ROWS:
            raise ValueError(f"n_compounds must be between 1 and {ROWS}")
        self.cell_lines = draw_cell_lines(self._truth_rng)
        self.pharmacology = {f"cpd{i + 1}": draw_pharmacology(self._truth_rng, f"cpd{i + 1}")
                             for i in range(n_compounds)}
        self.plates: list[dict] = []  # per plate: target, compounds, top_conc_nM, wells

    def spec(self) -> dict:
        return {
            "env": self.name,
            "task": ("Each compound is a potent on-target lead. Measure how selective it is over four related "
                     "off-targets. Every target has its own cell line, read by the same Gq-style cell assay and "
                     "competitive TR-FRET immunoassay (the ratio falls as the analyte rises), and a plate carries "
                     f"one cell line. Run up to {self.MAX_PLATES} plates, then report one selectivity window per "
                     "compound in log units: on-target pEC50 minus the pEC50 at its most potent off-target. Where "
                     "no off-target activity is seen up to the highest concentration tested, report a lower bound "
                     "instead of a number."),
            "targets": {"on_target": ON_TARGET, "off_targets": list(OFF_TARGETS)},
            "actions": {
                "run_plate": {"target": "one of " + ", ".join(TARGETS),
                              "compounds": "optional list of compound ids to put on the plate (default: all)",
                              "design": spec_help(PLATE_SPEC)},
                "report": {"windows": "object mapping compound id to {window: x, at_least: true or false}, or to "
                                      "null for no answer; compounds left out are not answered",
                           "window": spec_help(REPORT_SPEC)},
            },
            "labware": {"stamp_ul": STAMP_UL, "cells_ul": CELLS_UL, "ladder_well_max_ul": LADDER_MAX_UL,
                        "nominal_dead_volume_ul": RigProfile().dead_volume_ul, "stimulation_min": STIM_MIN,
                        "standard_series_nM": list(STD_SERIES_NM), "max_top_conc_nM": MAX_TOP_NM,
                        "plate_budget": self.MAX_PLATES},
            "analysis": ("drugdevgym.analysis.estimate_potency(wells) gives one call per compound on a plate (a "
                         "potency, a bound, or active without a potency); selectivity_window(on_potency, off_calls) "
                         "turns the calls of the four off-targets into a window or a lower bound."),
            "reward": reward_text(),
        }

    def reset(self) -> dict:
        return {"phase": "start", "compounds": list(self.pharmacology), "targets": list(TARGETS),
                "plates_left": self.MAX_PLATES - len(self.plates), "spec": self.spec()}

    def step(self, action: dict):
        err = self._check_open(action)
        if err:
            return self._reject(err)
        if action["type"] == "run_plate":
            return self._run_plate(action)
        if action["type"] == "report":
            return self._score(action)
        return self._reject(f"unknown action type: {action['type']}")

    def _compounds(self, action: dict) -> tuple[list[str] | None, str | None]:
        if "compounds" not in action:
            return list(self.pharmacology), None
        chosen = action["compounds"]
        if not isinstance(chosen, list) or not chosen or any(not isinstance(c, str) for c in chosen):
            return None, "compounds must be a non-empty list of compound ids"
        if any(c not in self.pharmacology for c in chosen):
            return None, "compounds names an unknown compound id"
        if len(set(chosen)) != len(chosen):
            return None, "compounds lists a compound twice"
        return [c for c in self.pharmacology if c in chosen], None

    def _run_plate(self, action: dict):
        extra = sorted(str(k) for k in set(action) - {"type", "target", "compounds", "design"})
        if extra:
            return self._reject(f"unknown field: {extra[0]}")
        if len(self.plates) >= self.MAX_PLATES:
            return self._reject("no plates left; report")
        target = action.get("target")
        if not isinstance(target, str) or target not in TARGETS:
            return self._reject("target must be one of: " + ", ".join(TARGETS))
        chosen, err = self._compounds(action)
        if err:
            return self._reject(err)
        design, err = check_plate(action.get("design"), len(chosen))
        if err:
            return self._reject(err)
        on_plate = [self.pharmacology[c][target] for c in chosen]
        wells, rig = simulate_plate(self.profile, self._noise, self.cell_lines[target], on_plate, design)
        self.plates.append({"target": target, "compounds": chosen, "top_conc_nM": design["top_conc_nM"],
                            "wells": wells})
        self._record(rig.log)
        self.cost += plate_cost(len(wells), rig.tips)
        return ({"phase": "data", "plate": len(self.plates), "target": target, "compounds": list(chosen),
                 "wells": copy.deepcopy(wells), "cost": round(self.cost, 2),
                 "plates_left": self.MAX_PLATES - len(self.plates)}, 0.0, False, {})

    def _basis(self, passes: list[bool]) -> dict[tuple[str, str], int]:
        """Index of the one plate each tested (compound, target) pair rests on.

        Among the plates that tested the pair: a plate that passes plate QC before one that fails,
        then the highest top concentration, then the earliest.
        """
        basis: dict[tuple[str, str], int] = {}
        for i, plate in enumerate(self.plates):
            rank = (passes[i], plate["top_conc_nM"])
            for cid in plate["compounds"]:
                held = basis.get((cid, plate["target"]))
                if held is None or rank > (passes[held], self.plates[held]["top_conc_nM"]):
                    basis[(cid, plate["target"])] = i
        return basis

    def _score(self, action: dict):
        extra = sorted(str(k) for k in set(action) - {"type", "windows"})
        if extra:
            return self._reject(f"unknown field: {extra[0]}")
        windows = action.get("windows")
        if not isinstance(windows, dict):
            return self._reject("windows must be an object mapping compound id to {window, at_least}")
        claims: dict[str, dict] = {}
        for cid, entry in windows.items():
            if not isinstance(cid, str) or cid not in self.pharmacology:
                return self._reject("windows names an unknown compound id")
            if entry is None:
                continue
            if not isinstance(entry, dict):
                return self._reject(f"{cid}: must be {{window, at_least}} or null")
            clean, err = validate(REPORT_SPEC, entry)
            if err:
                return self._reject(f"{cid}: {err}")
            claims[cid] = clean

        z = [zprime([w["ratio"] for w in plate["wells"] if w["role"] == "pos"],
                    [w["ratio"] for w in plate["wells"] if w["role"] == "neg"]) for plate in self.plates]
        passes = [bool(np.isfinite(v) and v >= ZPRIME_GATE) for v in z]
        basis = self._basis(passes)
        used = sorted(set(basis.values()))
        qc_pass = bool(used) and all(passes[i] for i in used)

        scores, detail, truth = [], {}, {}
        for cid, at in self.pharmacology.items():
            on = at[ON_TARGET].pec50
            window = true_window(on, [at[t].pec50 for t in OFF_TARGETS])
            achievable = on - P_FLOOR
            tops = [self.plates[basis[(cid, t)]]["top_conc_nM"] if (cid, t) in basis else None for t in OFF_TARGETS]
            supported = supported_window(on, tops, (cid, ON_TARGET) in basis)
            claim = claims.get(cid)
            score = 0.0 if claim is None else score_claim(claim["window"], claim["at_least"], window, achievable,
                                                           supported)
            scores.append(score)
            # tested_limit comes from the agent's own plates: the supported window is the on-target pEC50 minus it
            detail[cid] = {"claim": None if claim is None else dict(claim), "score": score,
                           "on_target_tested": (cid, ON_TARGET) in basis,
                           "tested_limit": None if None in tops else round(max(tested_limit(t) for t in tops), 3)}
            truth[cid] = {"on": round(on, 3), "off": {t: round(at[t].pec50, 3) for t in OFF_TARGETS},
                          "window": round(window, 3), "achievable": round(achievable, 3)}
        quality = float(np.mean(scores))  # no floor: a report of penalised claims scores below no report
        halved = quality > 0 and not qc_pass  # a failed gate halves credit, never a penalty
        reward = quality * (0.5 if halved else 1.0) - self._cost_penalty()
        self.done = True
        info = {"quality": quality, "qc_pass": qc_pass, "zprime": {str(i + 1): float(v) for i, v in enumerate(z)},
                "plates_used": [i + 1 for i in used], "cost": self.cost, "detail": detail, "truth": truth}
        return {"phase": "done"}, reward, True, info
