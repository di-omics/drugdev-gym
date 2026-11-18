"""Review environments: expert judgment as a verifiable task.

A run was executed from a fixed, sound protocol, possibly with one hidden deviation: a step
done wrong (no mixing, a reused tip, a lid left off, an under-filled well) or a condition that
went wrong (weak cells, an edge effect, reader drift, an uneven internal standard). The agent
sees what a reviewer sees (the intended protocol, the raw data, and the known value of one
reference compound) and must decide accept or repeat and name the most likely cause.

Truth is known by construction:

- "repeat" is right when the standard analysis of this data (named in spec()) is off by more
  than the tolerance for any compound, leaves a compound unanswered, or the run fails its
  QC gate. A deviation that leaves every answer within tolerance should be accepted.
- the cause is the deviation that was planted ("none" when there was none).

The intended protocol is sound on the nominal and randomized rigs, where a clean run is
labelled repeat a few percent of the time. On the harsh rig it is not: a large share of clean
runs fail on their own, so results there say little about judgment.

The same cases can be shown to a human expert, so their calls become the reference level.
"""
from __future__ import annotations

import copy
from dataclasses import replace

import numpy as np

from ..analysis import MIN_R2, MIN_SPAN, control_zprime, estimate_half_life, estimate_pec50
from .base import Env
from .dose_response import ZPRIME_GATE, Compound, draw_assay, simulate_plate
from .dose_response import check_design as check_dr
from .stability import HL_HI, HL_LO, T0_CV_GATE, simulate_stability, t0_cv
from .stability import check_design as check_st

VERDICT_WEIGHT, CAUSE_WEIGHT = 0.6, 0.4

PLATE_PROTOCOL = {"top_conc_nM": 30000.0, "dilution_factor": 4.0, "n_points": 9, "replicates": 2,
                  "transfer_ul": 25.0, "mix": True, "new_tips": True, "avoid_edges": False, "lid": True,
                  "standards": True, "controls_per_side": 6}
PLATE_FAULTS = ("none", "unmixed_ladder", "starved_ladder", "tip_carryover", "edge_effect", "reader_drift",
                "weak_cells")
PLATE_FAULT_HELP = {
    "none": "the run went as written",
    "unmixed_ladder": "the serial dilution was not mixed between steps",
    "starved_ladder": "the ladder wells held too little volume: one replicate was short-dosed, the other got almost nothing",
    "tip_carryover": "one tip was reused down each ladder instead of a fresh tip per step",
    "edge_effect": "cells on the outer ring responded much less than the rest",
    "reader_drift": "the signal drifts strongly from row A to row H",
    "weak_cells": "the cells barely respond, so the assay window has collapsed",
}
STARVED_TRANSFER_UL = 7.0  # ladder transfer volume in a starved run (intended: 25)
EDGE_FAULT_LOSS = 0.5  # edge-ring cell response lost in an edge-effect run (rig default: 0.15)
DRIFT_FAULT_GRADIENT = 0.35  # row A to H signal drift in a reader-drift run (rig default: 0.04)

TIMECOURSE_PROTOCOL = {"timepoints_min": [0.0, 15.0, 30.0, 60.0, 120.0, 240.0], "replicates": 2,
                       "internal_standard": True, "mix_after_spike": True, "lid": True, "avoid_edges": True}
TIMECOURSE_FAULTS = ("none", "unmixed_spike", "evaporation", "internal_standard_failure", "short_aliquots")
TIMECOURSE_FAULT_HELP = {
    "none": "the run went as written",
    "unmixed_spike": "the matrix was not mixed after spiking, before time zero",
    "evaporation": "the plate sat unlidded on the outer ring",
    "internal_standard_failure": "the internal standard was dispensed unevenly",
    "short_aliquots": "the incubation wells were under-filled, so the last aliquots came up short or empty",
}
UNDERFILLED_MATRIX_UL = 105.0  # matrix volume in a short-aliquot run (intended: 200)
FAILED_IS_SD = 0.3  # internal-standard CV in a failed-dispense run (normal: 0.03)
RESPONSE_SD = 0.5  # log-scale spread of per-compound MS response, so absolute area carries no information
IS_LEVEL_SD = 0.3  # log-scale spread of the internal-standard level from run to run


class _Review(Env):
    """Shared one-step judging logic. Subclasses build the case and its truth."""

    salted = True  # the planted deviation must not be predictable from another environment's seed
    faults: tuple[str, ...] = ()
    fault_help: dict[str, str] = {}
    protocol: dict = {}
    tolerance = 0.5
    task = ""
    data_key = ""
    label_rule = ""

    def __init__(self, seed: int = 0, rig="nominal"):
        super().__init__(seed, rig)
        r = self._truth_rng
        # four in ten runs are clean; the rest carry one deviation each
        self.fault = "none" if r.random() < 0.4 else str(r.choice(self.faults[1:]))
        self.data, self.reference, self.should_repeat, self.detail, self._executed = self._build()

    def _build(self):
        raise NotImplementedError

    def spec(self) -> dict:
        return {
            "env": self.name,
            "task": self.task,
            "actions": {"judge": {"verdict": "accept or repeat", "cause": dict(self.fault_help)}},
            "intended_protocol": copy.deepcopy(self.protocol),
            "label": self.label_rule,
            "reward": f"{VERDICT_WEIGHT} for the right verdict plus {CAUSE_WEIGHT} for the right cause",
        }

    def reset(self) -> dict:
        # the reference compound's true value is known to the reviewer, as on a real plate
        return {"phase": "review", self.data_key: copy.deepcopy(self.data), "reference": dict(self.reference),
                "spec": self.spec()}

    def truth(self) -> dict:
        """The planted deviation and whether the data misleads. Never part of an observation."""
        return {"fault": self.fault, "should_repeat": self.should_repeat, **self.detail}

    def step(self, action: dict):
        err = self._check_open(action)
        if err:
            return self._reject(err)
        if action["type"] != "judge":
            return self._reject(f"unknown action type: {action['type']}")
        verdict, cause = action.get("verdict"), action.get("cause")
        if not isinstance(verdict, str) or verdict not in ("accept", "repeat"):
            return self._reject("verdict must be accept or repeat")
        if not isinstance(cause, str) or cause not in self.faults:
            return self._reject(f"cause must be one of: {', '.join(self.faults)}")
        verdict_ok = (verdict == "repeat") == self.should_repeat
        cause_ok = cause == self.fault
        self._record(self._executed)  # the executed run is only disclosed once the call is in
        self.done = True
        info = {"verdict_ok": verdict_ok, "cause_ok": cause_ok, "truth": self.truth()}
        return {"phase": "done"}, VERDICT_WEIGHT * verdict_ok + CAUSE_WEIGHT * cause_ok, True, info

    def _worst(self, errors: dict) -> float | None:
        known = [abs(e) for e in errors.values() if e is not None]
        return round(max(known), 3) if known else None


class PlateReview(_Review):
    name = "plate_review"
    faults, fault_help, protocol = PLATE_FAULTS, PLATE_FAULT_HELP, PLATE_PROTOCOL
    data_key = "wells"
    task = ("Review one dose-response plate (4 compounds, duplicates, competitive TR-FRET readout). It was "
            "meant to be run with the intended protocol; one thing may have gone wrong. Accept or repeat, "
            "and name the most likely cause. cpd1 is a reference compound whose true pEC50 is given.")
    label_rule = (f"Repeat is right when any of these holds. (1) The control Z' on the raw ratio is below "
                  f"{ZPRIME_GATE}. (2) The standard analysis leaves a compound unanswered: each compound's ratios are "
                  f"converted to analyte through the on-plate standard curve and fitted with a 4-parameter curve, "
                  f"and a fit is dropped when its R2 is below {MIN_R2}, its midpoint sits on the edge of the search "
                  f"range, or its span is under {MIN_SPAN} of the control window. (3) Any fitted pEC50 is more than "
                  f"0.5 log units from the truth. Otherwise accept is right, even if something went wrong.")

    def _build(self):
        r = self._truth_rng
        assay = draw_assay(r)
        # potencies the intended range covers well, so a clean run is almost always right
        compounds = [Compound(f"cpd{i + 1}", float(r.uniform(6.0, 8.8)), float(r.uniform(0.6, 1.0)),
                              float(r.uniform(0.8, 1.3)), float(np.exp(r.normal(np.log(2.0), 0.8))))
                     for i in range(4)]
        # cpd1 is the reference: mid-range potency, full agonist, not sticky, value known to the reviewer
        compounds[0] = replace(compounds[0], pec50=float(r.uniform(7.0, 8.0)), emax=1.0, hill=1.0, sticky=1.0)
        reference = {"id": "cpd1", "pEC50": round(compounds[0].pec50, 2)}
        design, profile = dict(PLATE_PROTOCOL), self.profile
        if self.fault == "unmixed_ladder":
            design["mix"] = False
        elif self.fault == "starved_ladder":
            design["transfer_ul"] = STARVED_TRANSFER_UL
        elif self.fault == "tip_carryover":
            design["new_tips"] = False
        elif self.fault == "edge_effect":
            profile = replace(profile, edge_cell_loss=EDGE_FAULT_LOSS)
        elif self.fault == "reader_drift":
            profile = replace(profile, row_gradient=DRIFT_FAULT_GRADIENT)
        elif self.fault == "weak_cells":
            assay = replace(assay, top=assay.basal * float(r.uniform(1.3, 2.0)))
        design, err = check_dr(design, len(compounds))
        assert err is None, err
        wells, rig = simulate_plate(profile, self._noise, assay, compounds, design)
        estimates = estimate_pec50(wells, use_standards=True)
        errors = {c.cid: (round(estimates[c.cid] - c.pec50, 3) if c.cid in estimates else None) for c in compounds}
        z = control_zprime(wells)
        misleading = any(e is None or abs(e) > self.tolerance for e in errors.values())
        qc_fail = not np.isfinite(z) or z < ZPRIME_GATE
        detail = {"errors": errors, "worst_error": self._worst(errors), "zprime": round(float(z), 3),
                  "qc_fail": bool(qc_fail)}
        return wells, reference, bool(misleading or qc_fail), detail, rig.log


class TimecourseReview(_Review):
    name = "timecourse_review"
    faults, fault_help, protocol = TIMECOURSE_FAULTS, TIMECOURSE_FAULT_HELP, TIMECOURSE_PROTOCOL
    data_key = "samples"
    tolerance = 0.3
    task = ("Review one matrix stability time course (3 compounds, duplicate incubations, LC-MS with an "
            "internal standard). It was meant to be run with the intended protocol; one thing may have gone "
            "wrong. Accept or repeat, and name the most likely cause. cpd1 is a control compound whose true "
            "half-life is given.")
    label_rule = (f"Repeat is right when any of these holds. (1) The time-zero replicate CV of analyte over internal "
                  f"standard is above {T0_CV_GATE:.0%} for any compound. (2) The standard analysis leaves a compound "
                  f"unanswered or is off: half-life by log-linear regression of analyte over internal standard on "
                  f"time, dropping points under 2% of the time-zero mean, clipped to {HL_LO:g} to {HL_HI:g} min, more "
                  f"than 0.3 log units (2-fold) from the truth. Otherwise accept is right, even if something went wrong.")

    def _build(self):
        r = self._truth_rng
        half_lives = {f"cpd{i + 1}": float(10 ** r.uniform(1.0, np.log10(600.0))) for i in range(3)}
        # cpd1 is the control compound: its half-life is known to the reviewer
        half_lives["cpd1"] = float(10 ** r.uniform(np.log10(40.0), np.log10(200.0)))
        reference = {"id": "cpd1", "half_life_min": round(half_lives["cpd1"], 1)}
        response = {cid: float(np.exp(r.normal(0.0, RESPONSE_SD))) for cid in half_lives}
        is_level = float(np.exp(r.normal(0.0, IS_LEVEL_SD)))
        design, is_sd, matrix_ul = dict(TIMECOURSE_PROTOCOL), 0.03, 200.0
        if self.fault == "unmixed_spike":
            design["mix_after_spike"] = False
        elif self.fault == "evaporation":
            design["lid"], design["avoid_edges"] = False, False
        elif self.fault == "internal_standard_failure":
            is_sd = FAILED_IS_SD
        elif self.fault == "short_aliquots":
            matrix_ul = UNDERFILLED_MATRIX_UL
        design, err = check_st(design, len(half_lives))
        assert err is None, err
        samples, rig = simulate_stability(self.profile, self._noise, half_lives, design, is_sd=is_sd,
                                          matrix_ul=matrix_ul)
        for s in samples:  # compound-specific MS response and a run-specific internal-standard level
            s["area"] = round(s["area"] * response[s["id"]], 1)
            s["is_area"] = round(s["is_area"] * is_level, 1)
        estimates = estimate_half_life(samples)
        errors = {}
        for cid, hl in half_lives.items():
            est = estimates.get(cid)
            errors[cid] = None if est is None else round(float(np.log10(est) - np.log10(np.clip(hl, HL_LO, HL_HI))), 3)
        cvs = t0_cv(samples)
        misleading = any(e is None or abs(e) > self.tolerance for e in errors.values())
        qc_fail = any(not np.isfinite(v) or v > T0_CV_GATE for v in cvs.values())
        detail = {"errors": errors, "worst_error": self._worst(errors),
                  "t0_cv": {k: round(v, 3) for k, v in cvs.items()}, "qc_fail": bool(qc_fail)}
        return samples, reference, bool(misleading or qc_fail), detail, rig.log
