"""Review environments: expert judgment as a verifiable task.

A run was executed from a fixed, sound protocol, possibly with one hidden execution fault.
The agent sees what a reviewer sees (the intended protocol, the raw data, and the known value
of one reference compound) and must decide accept or repeat and name the most likely cause. Truth is known by construction:

- "repeat" is right when the standard analysis of this data misses the true answer by more
  than the tolerance, or the run fails its own QC. A fault that does no harm should be accepted.
- the cause is the fault that was planted ("none" when there was none).

The same cases can be shown to a human expert, so their calls become the reference level.
"""
from __future__ import annotations

from dataclasses import replace

import numpy as np

from ..analysis import control_zprime, estimate_half_life, estimate_pec50
from .base import Env
from .dose_response import Compound, draw_assay, simulate_plate
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
    "starved_ladder": "the ladder wells held too little volume, so assay wells were short-dosed",
    "tip_carryover": "one tip was reused down the ladder on a rig with heavy carryover",
    "edge_effect": "cells on the outer ring responded much less than the rest",
    "reader_drift": "the signal drifts strongly from row A to row H",
    "weak_cells": "the cells barely respond, so the assay window has collapsed",
}

TIMECOURSE_PROTOCOL = {"timepoints_min": [0.0, 15.0, 30.0, 60.0, 120.0, 240.0], "replicates": 2,
                       "internal_standard": True, "mix_after_spike": True, "lid": True, "avoid_edges": True}
TIMECOURSE_FAULTS = ("none", "unmixed_spike", "evaporation", "internal_standard_failure", "short_aliquots")
TIMECOURSE_FAULT_HELP = {
    "none": "the run went as written",
    "unmixed_spike": "the matrix was not mixed after spiking, before time zero",
    "evaporation": "the plate sat unlidded on the outer ring in a dry incubator",
    "internal_standard_failure": "the internal standard was dispensed unevenly",
    "short_aliquots": "the incubation wells ran low, so the later aliquots were short",
}


class _Review(Env):
    """Shared one-step judging logic. Subclasses build the case and its truth."""

    faults: tuple[str, ...] = ()
    fault_help: dict[str, str] = {}
    protocol: dict = {}
    reference: dict = {}
    tolerance = 0.5
    task = ""
    data_key = ""

    def __init__(self, seed: int = 0, rig="nominal"):
        super().__init__(seed, rig)
        r = self._truth_rng
        # four in ten runs are clean; the rest carry one fault each
        self.fault = "none" if r.random() < 0.4 else str(r.choice(self.faults[1:]))
        self.data, self.reference, self.should_repeat, self.detail = self._build()

    def _build(self):
        raise NotImplementedError

    def spec(self) -> dict:
        return {
            "env": self.name,
            "task": self.task,
            "actions": {"judge": {"verdict": "accept or repeat", "cause": dict(self.fault_help)}},
            "intended_protocol": self.protocol,
            "reward": (f"{VERDICT_WEIGHT} for the right verdict plus {CAUSE_WEIGHT} for the right cause. Repeat is "
                       f"right when the standard analysis of this data is off by more than {self.tolerance} log "
                       "units for any compound, leaves a compound unanswered, or the run fails its own QC."),
        }

    def reset(self) -> dict:
        # the reference compound's true value is known to the reviewer, as on a real plate
        return {"phase": "review", self.data_key: self.data, "reference": self.reference, "spec": self.spec()}

    def truth(self) -> dict:
        """The planted fault and whether the data misleads. Never part of an observation."""
        return {"fault": self.fault, "should_repeat": self.should_repeat, **self.detail}

    def step(self, action: dict):
        err = self._check_open(action)
        if err:
            return self._reject(err)
        if action["type"] != "judge":
            return self._reject(f"unknown action type: {action['type']}")
        verdict, cause = action.get("verdict"), action.get("cause")
        if verdict not in ("accept", "repeat"):
            return self._reject("verdict must be accept or repeat")
        if cause not in self.faults:
            return self._reject(f"cause must be one of: {', '.join(self.faults)}")
        verdict_ok = (verdict == "repeat") == self.should_repeat
        cause_ok = cause == self.fault
        self.done = True
        info = {"verdict_ok": verdict_ok, "cause_ok": cause_ok, "truth": self.truth()}
        return {"phase": "done"}, VERDICT_WEIGHT * verdict_ok + CAUSE_WEIGHT * cause_ok, True, info


class PlateReview(_Review):
    name = "plate_review"
    faults, fault_help, protocol = PLATE_FAULTS, PLATE_FAULT_HELP, PLATE_PROTOCOL
    data_key = "wells"
    task = ("Review one dose-response plate (4 compounds, duplicates, competitive TR-FRET readout). It was "
            "meant to be run with the intended protocol; one thing may have gone wrong. Accept or repeat, "
            "and name the most likely cause.")

    def _build(self):
        r = self._truth_rng
        assay = draw_assay(r)
        # potencies the intended range covers well, so a clean run is almost always right
        compounds = [Compound(f"cpd{i + 1}", float(r.uniform(6.0, 8.8)), float(r.uniform(0.6, 1.0)),
                              float(r.uniform(0.8, 1.3)), float(np.exp(r.normal(np.log(2.0), 0.8))))
                     for i in range(4)]
        # cpd1 is the reference compound: mid-range potency, full agonist, value known to the reviewer
        compounds[0] = replace(compounds[0], pec50=float(r.uniform(7.0, 8.0)), emax=1.0, hill=1.0)
        reference = {"id": "cpd1", "pEC50": round(compounds[0].pec50, 2)}
        design, profile = dict(PLATE_PROTOCOL), self.profile
        if self.fault == "unmixed_ladder":
            design["mix"] = False
        elif self.fault == "starved_ladder":
            design["transfer_ul"] = 7.0
        elif self.fault == "tip_carryover":
            design["new_tips"] = False
            profile = replace(profile, carryover=0.05)
        elif self.fault == "edge_effect":
            profile = replace(profile, edge_cell_loss=0.5)
        elif self.fault == "reader_drift":
            profile = replace(profile, row_gradient=0.35)
        elif self.fault == "weak_cells":
            assay = replace(assay, top=assay.basal * float(r.uniform(1.3, 2.0)))
        design, err = check_dr(design, len(compounds))
        assert err is None, err
        wells, rig = simulate_plate(profile, self._noise, assay, compounds, design)
        self._record(rig.log)
        estimates = estimate_pec50(wells, use_standards=True)
        errors = {c.cid: (round(estimates[c.cid] - c.pec50, 3) if c.cid in estimates else None) for c in compounds}
        z = control_zprime(wells)
        misleading = any(e is None or abs(e) > self.tolerance for e in errors.values())
        qc_fail = not np.isfinite(z) or z < 0.5
        detail = {"errors": errors, "zprime": round(float(z), 3), "qc_fail": bool(qc_fail)}
        return wells, reference, bool(misleading or qc_fail), detail


class TimecourseReview(_Review):
    name = "timecourse_review"
    faults, fault_help, protocol = TIMECOURSE_FAULTS, TIMECOURSE_FAULT_HELP, TIMECOURSE_PROTOCOL
    data_key = "samples"
    tolerance = 0.3
    task = ("Review one matrix stability time course (3 compounds, duplicate incubations, LC-MS with an "
            "internal standard). It was meant to be run with the intended protocol; one thing may have gone "
            "wrong. Accept or repeat, and name the most likely cause.")

    def _build(self):
        r = self._truth_rng
        half_lives = {f"cpd{i + 1}": float(10 ** r.uniform(1.0, np.log10(600.0))) for i in range(3)}
        # cpd1 is the control compound: its half-life is known to the reviewer
        half_lives["cpd1"] = float(10 ** r.uniform(np.log10(40.0), np.log10(200.0)))
        reference = {"id": "cpd1", "half_life_min": round(half_lives["cpd1"], 1)}
        design, profile, is_sd = dict(TIMECOURSE_PROTOCOL), self.profile, 0.03
        if self.fault == "unmixed_spike":
            design["mix_after_spike"] = False
        elif self.fault == "evaporation":
            design["lid"], design["avoid_edges"] = False, False
            profile = replace(profile, evap_edge_ul_per_h=10.0)
        elif self.fault == "internal_standard_failure":
            is_sd = 0.3
        elif self.fault == "short_aliquots":
            profile = replace(profile, dead_volume_ul=110.0)
        design, err = check_st(design, len(half_lives))
        assert err is None, err
        samples, rig = simulate_stability(profile, self._noise, half_lives, design, is_sd=is_sd)
        self._record(rig.log)
        estimates = estimate_half_life(samples)
        errors = {}
        for cid, hl in half_lives.items():
            est = estimates.get(cid)
            errors[cid] = None if est is None else round(float(np.log10(est) - np.log10(np.clip(hl, HL_LO, HL_HI))), 3)
        cvs = t0_cv(samples)
        misleading = any(e is None or abs(e) > self.tolerance for e in errors.values())
        qc_fail = any(not np.isfinite(v) or v > T0_CV_GATE for v in cvs.values())
        detail = {"errors": errors, "t0_cv": {k: round(v, 3) for k, v in cvs.items()}, "qc_fail": bool(qc_fail)}
        return samples, reference, bool(misleading or qc_fail), detail
