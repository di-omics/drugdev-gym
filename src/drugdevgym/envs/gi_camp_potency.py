"""Gi cAMP potency environment: set up a forskolin-stimulated cAMP assay and report agonist pEC50s.

The receptor is Gi-coupled: forskolin raises cAMP and the agonist lowers it. cAMP is read by a
competitive TR-FRET immunoassay, so the ratio RISES as the agonist acts, the opposite direction
to dose_response. The agent chooses the forskolin concentration (it may spend a plate on a
forskolin titration first), whether the PDE inhibitor IBMX is in the buffer, cells per well,
stimulation time, the agonist ladder, standards and controls. The reward is the accuracy of the
reported pEC50s against the hidden truth, halved when a potency plate fails its own QC
(control Z' below 0.5), minus a small cost term.

The cells are plated the day before at cells_k per well and adhere overnight, which is why the
outer ring responds less (the rig's edge_cell_loss). On the day the medium is replaced by 30 uL
of stimulation buffer, with the IBMX in it when IBMX is used. The simulator runs the plating
and the buffer exchange as one 30 uL fill that carries the cells and the IBMX.
"""
from __future__ import annotations

import copy
from dataclasses import dataclass

import numpy as np

from ..analysis import zprime
from ..chem import Kit, accumulated, hill, uninhibited_fraction
from ..rig import COLS, EDGE, ROWS, Plate, Rig, RigProfile, usable_wells, well_name
from .base import COST_WEIGHT, Env, log_error_score, spec_help, validate
from .dose_response import Compound, plate_cost

# protocol volumes for one well of a 96-well plate
CELLS_UL = 30.0  # stimulation buffer over the cells plated the day before; it carries the IBMX when IBMX is used
STAMP_UL = 10.0  # agonist from the ladder, vehicle, or reference agonist
FSK_UL = 10.0  # forskolin, added last: it starts the stimulation
STIM_UL = CELLS_UL + STAMP_UL + FSK_UL
DETECT_UL = 50.0  # lysis buffer with the two TR-FRET detection reagents, one stimulation volume
FINAL_UL = STIM_UL + DETECT_UL
INTO_STIM = STIM_UL / STAMP_UL  # agonist ladder wells are 5x the stimulation concentration
FSK_INTO_STIM = STIM_UL / FSK_UL  # so is the forskolin stock
LADDER_MAX_UL = 300.0  # working volume of a dilution-plate well

STD_SERIES_NM = (0.04, 0.17, 0.7, 2.8, 11.0, 44.0, 178.0, 712.0)  # cAMP standards, nM in the read well, 4-fold
FSK_MAX_UM = 100.0  # highest forskolin an agent can ask for, and the top of the titration
FSK_SERIES_UM = tuple(round(FSK_MAX_UM / 10.0 ** (j / 2.0), 4) for j in range(8))  # half-log titration series
FSK_FOLD = 10.0 ** 0.5  # fold change between titration points
FSK_LADDER_UL = 40.0  # transfer volume of the forskolin titration ladder
IBMX_UM = 500.0  # IBMX in the stimulation when it is included (0.5 mM, the usual kit recommendation)
IBMX_KI_UM = 20.0  # lumped inhibition constant of IBMX over the phosphodiesterases of the cell
FSK_HILL = 1.0  # Hill slope of the forskolin response
ZPRIME_GATE = 0.5  # control Z' a potency plate must reach
MIN_CONTROLS = 3  # controls a side below which a plate has no Z'
ANSWER_RANGE = (0.0, 14.0)  # a reported pEC50 outside this is malformed
TOLERANCE = 0.5  # pEC50 error (log units) at which the credit for a compound reaches zero

# hidden truth, drawn once per seed
KIT_K_NM = (2.5, 7.5)  # cAMP at half displacement of the kit, nM in the read well
FSK_EC50_UM = (0.5, 10.0)  # forskolin EC50 of the cell line, log-uniform
CAMP_RATE = (0.5, 8.0)  # fmol cAMP per 1000 cells per minute at saturating forskolin, log-uniform
BASAL_FRAC = (0.01, 0.05)  # cyclase rate without forskolin, as a fraction of the forskolin maximum
PDE_PER_MIN = (0.1, 0.5)  # first-order cAMP breakdown with no inhibitor (half-life 1.4 to 7 min)
GI_MAX = (0.7, 0.95)  # fraction of the forskolin response the reference full agonist inhibits
EFFICACY = (0.6, 1.0)  # efficacy of a test agonist relative to the reference agonist
CELL_CV = (0.04, 0.08)  # well-to-well cell variation
PEC50 = (5.5, 9.5)  # agonist potency
HILL = (0.8, 1.3)  # agonist Hill slope
STICKY_MEDIAN, STICKY_LOG_SD = 2.0, 0.8  # tip-carryover multiplier, log-normal as in dose_response

TITRATION_SPEC = {
    "cells_k": ("float", 1.0, 80.0, "thousand cells per well"),
    "stim_min": ("float", 5.0, 120.0, "stimulation time in minutes"),
    "ibmx": ("bool", None, None, "IBMX (PDE inhibitor) in the stimulation buffer"),
    "replicates": ("int", 1, 4, "wells per forskolin concentration"),
}

DESIGN_SPEC = {
    "forskolin_uM": ("float", 0.0, FSK_MAX_UM, "forskolin in the stimulation, in every cell well; 0 leaves it out"),
    "ibmx": ("bool", None, None, "IBMX (PDE inhibitor) in the stimulation buffer"),
    "cells_k": ("float", 1.0, 80.0, "thousand cells per well"),
    "stim_min": ("float", 5.0, 120.0, "stimulation time in minutes"),
    "top_conc_nM": ("float", 1.0, 100000.0, "highest agonist concentration in the stimulation"),
    "dilution_factor": ("float", 1.5, 10.0, "fold change between ladder points"),
    "n_points": ("int", 4, 12, "ladder points per compound"),
    "replicates": ("int", 1, 4, "assay wells per ladder point"),
    "transfer_ul": ("float", 2.0, 100.0, "ladder transfer volume; each ladder well holds transfer_ul * dilution_factor"),
    "mix": ("bool", None, None, "mix each ladder well after the transfer into it"),
    "standards": ("bool", None, None, "include the 8-point cAMP standard curve (8 wells)"),
    "controls_per_side": ("int", 0, 16, "forskolin-only (neg) and forskolin plus reference agonist (pos) wells, each"),
}


@dataclass(frozen=True)
class GiAssay:
    """One cell line on one kit. rate is fmol cAMP per 1000 cells per minute at saturating forskolin."""

    kit: Kit
    fsk_ec50: float
    rate: float
    basal_frac: float
    pde: float
    gi_max: float
    cell_cv: float


def _log_uniform(rng: np.random.Generator, lo: float, hi: float) -> float:
    return float(np.exp(rng.uniform(np.log(lo), np.log(hi))))


def draw_gi_assay(rng: np.random.Generator) -> GiAssay:
    return GiAssay(Kit(k=float(rng.uniform(*KIT_K_NM))), fsk_ec50=_log_uniform(rng, *FSK_EC50_UM),
                   rate=_log_uniform(rng, *CAMP_RATE), basal_frac=float(rng.uniform(*BASAL_FRAC)),
                   pde=float(rng.uniform(*PDE_PER_MIN)), gi_max=float(rng.uniform(*GI_MAX)),
                   cell_cv=float(rng.uniform(*CELL_CV)))


def draw_agonist(rng: np.random.Generator, cid: str, gi_max: float) -> Compound:
    """emax is the fraction of the forskolin response this agonist inhibits at saturation."""
    return Compound(cid, pec50=float(rng.uniform(*PEC50)), emax=gi_max * float(rng.uniform(*EFFICACY)),
                    hill=float(rng.uniform(*HILL)), sticky=float(np.exp(rng.normal(np.log(STICKY_MEDIAN), STICKY_LOG_SD))))


def cyclase_drive(assay: GiAssay, forskolin_um: float, inhibition: float) -> float:
    """Synthesis rate as a fraction of the forskolin maximum.

    Basal plus the forskolin response, of which the Gi pathway removes the fraction
    `inhibition`. The inhibition is noncompetitive with forskolin: the same fraction at any
    forskolin concentration.
    """
    stimulated = float(hill(forskolin_um, assay.fsk_ec50, FSK_HILL))
    return assay.basal_frac + (1.0 - assay.basal_frac) * stimulated * (1.0 - inhibition)


def camp_fmol(assay: GiAssay, cells_k: float, forskolin_um: float, inhibition: float, ibmx_um: float,
              minutes: float) -> float:
    """cAMP in one well at lysis (fmol): constant synthesis against first-order PDE breakdown.

    IBMX leaves 1 / (1 + [IBMX] / Ki) of the PDE activity. With IBMX cAMP piles up almost
    linearly in time; without it the well sits at the steady state rate / k within minutes.
    """
    loss = assay.pde * uninhibited_fraction(ibmx_um, IBMX_KI_UM)
    return cells_k * assay.rate * cyclase_drive(assay, forskolin_um, inhibition) * accumulated(loss, minutes)


def plate_order() -> list[tuple[int, int]]:
    """Wells in the order a plate is filled: the inner 60 first, the outer ring only when they run out."""
    return usable_wells(True) + [w for w in usable_wells(False) if EDGE[w]]


def check_titration(design) -> tuple[dict | None, str | None]:
    return validate(TITRATION_SPEC, design)


def check_design(design, n_compounds: int) -> tuple[dict | None, str | None]:
    clean, err = validate(DESIGN_SPEC, design)
    if err:
        return None, err
    if clean["transfer_ul"] * clean["dilution_factor"] > LADDER_MAX_UL:
        return None, f"ladder well would hold more than {LADDER_MAX_UL:.0f} uL; lower transfer_ul or dilution_factor"
    need = n_compounds * clean["n_points"] * clean["replicates"] + 2 * clean["controls_per_side"]
    need += len(STD_SERIES_NM) if clean["standards"] else 0
    if need > ROWS * COLS:
        return None, f"design needs {need} wells but the plate has {ROWS * COLS}"
    return clean, None


def _add_cells(rig: Rig, plate: Plate, w, cells_k: float, ibmx: bool) -> None:
    solutes = {"cells": cells_k / CELLS_UL}
    if ibmx:
        solutes["ibmx"] = IBMX_UM * STIM_UL / CELLS_UL
    rig.fill(plate, w, CELLS_UL, solutes, what=f"cells {cells_k:g}k per well {'with' if ibmx else 'without'} IBMX")


def _add_standard(rig: Rig, plate: Plate, w, std_nm: float) -> None:
    rig.fill(plate, w, STIM_UL, {"cAMP": std_nm * FINAL_UL / STIM_UL}, what=f"cAMP standard {std_nm:g} nM")


def _stimulate_and_read(rig: Rig, rng: np.random.Generator, assay: GiAssay, plate: Plate, filled: list,
                        cell_wells: list, minutes: float, inhibition_of) -> np.ndarray:
    """Lidded stimulation, cAMP made in every cell well, detection mix into every well, one read.

    The outer ring responds less by the rig's edge_cell_loss: the cells were plated the day before.
    """
    before = plate.vol.copy()
    rig.incubate(plate, minutes, lid=True)
    volume = (before + plate.vol) / 2.0  # what the cells saw, averaged over the stimulation
    camp, cells = plate.a("cAMP"), plate.a("cells")
    fsk, ibmx = plate.a("fsk"), plate.a("ibmx")
    for w in cell_wells:
        edge = 1.0 - rig.p.edge_cell_loss if EDGE[w] else 1.0
        responding = cells[w] * max(0.0, 1.0 + rng.normal(0.0, assay.cell_cv)) * edge
        camp[w] = camp_fmol(assay, responding, fsk[w] / volume[w], inhibition_of(w, volume[w]), ibmx[w] / volume[w],
                            minutes)
    for w in filled:
        rig.fill(plate, w, DETECT_UL, what="detection mix")
    conc = np.where(plate.vol > 0, camp / np.maximum(plate.vol, 1e-9), 0.0)
    return rig.read(assay.kit.ratio(conc))


def simulate_titration(profile: RigProfile, rng: np.random.Generator, assay: GiAssay,
                       design: dict) -> tuple[list[dict], Rig]:
    """Forskolin titration with no agonist: the half-log series, unstimulated wells and the standard curve."""
    rig = Rig(profile, rng)
    dil, plate = Plate("dilution"), Plate("assay")
    free = iter(plate_order())
    layout: list[tuple] = [(next(free), "std", None, s) for s in STD_SERIES_NM]  # (well, role, forskolin uM, standard nM)
    fsk_wells: dict[int, list] = {}
    for j, f in enumerate(FSK_SERIES_UM):
        for _ in range(design["replicates"]):
            w = next(free)
            layout.append((w, "fsk", f, None))
            fsk_wells.setdefault(j, []).append(w)
    for _ in range(design["replicates"]):
        layout.append((next(free), "unstim", 0.0, None))

    # forskolin ladder, high to low, mixed, a fresh tip per step
    diluent = FSK_LADDER_UL * (FSK_FOLD - 1.0)
    rig.fill(dil, (0, 0), diluent + FSK_LADDER_UL, {"fsk": FSK_MAX_UM * FSK_INTO_STIM},
             what=f"forskolin {FSK_MAX_UM * FSK_INTO_STIM:g} uM")
    for j in range(1, len(FSK_SERIES_UM)):
        rig.fill(dil, (0, j), diluent, what="diluent")
    for j in range(1, len(FSK_SERIES_UM)):
        rig.transfer(dil, (0, j - 1), dil, (0, j), FSK_LADDER_UL, new_tip=True, mix_after=True)

    cell_wells = []
    for w, role, _f, std in layout:
        if role == "std":
            _add_standard(rig, plate, w, std)
            continue
        cell_wells.append(w)
        _add_cells(rig, plate, w, design["cells_k"], design["ibmx"])
        rig.fill(plate, w, STAMP_UL, what="vehicle")
        if role == "unstim":
            rig.fill(plate, w, FSK_UL, what="stimulation buffer")
    for j in reversed(range(len(FSK_SERIES_UM))):
        for w in fsk_wells[j]:
            rig.transfer(dil, (0, j), plate, w, FSK_UL, new_tip=True, mix_after=True)

    ratio = _stimulate_and_read(rig, rng, assay, plate, [w for w, *_ in layout], cell_wells, design["stim_min"],
                                lambda w, volume: 0.0)
    wells = [{"well": well_name(w), "role": role, "forskolin_uM": f, "std_nM": std,
              "ratio": round(float(ratio[w]), 1)} for w, role, f, std in layout]
    return wells, rig


def simulate_potency(profile: RigProfile, rng: np.random.Generator, assay: GiAssay, compounds: list[Compound],
                     design: dict) -> tuple[list[dict], Rig]:
    """Run one validated potency design on the rig. Returns the observed wells and the rig (log, tips)."""
    rig = Rig(profile, rng, sticky={c.cid: c.sticky for c in compounds})
    n, reps, t, df = design["n_points"], design["replicates"], design["transfer_ul"], design["dilution_factor"]
    diluent = t * (df - 1.0)
    dil, plate = Plate("dilution"), Plate("assay")
    free = iter(plate_order())
    layout: list[tuple] = []  # (well, role, cid, nominal nM, standard nM)
    for _ in range(design["controls_per_side"]):
        layout.append((next(free), "neg", None, None, None))
        layout.append((next(free), "pos", None, None, None))
    if design["standards"]:
        for s in STD_SERIES_NM:
            layout.append((next(free), "std", None, None, s))
    nominal = [design["top_conc_nM"] / df ** j for j in range(n)]
    cpd_wells: dict[tuple[str, int], list] = {}
    for c in compounds:
        for j in range(n):
            for _ in range(reps):
                w = next(free)
                layout.append((w, "cpd", c.cid, nominal[j], None))
                cpd_wells.setdefault((c.cid, j), []).append(w)

    # agonist ladders, one row of the dilution plate per compound, high to low, a fresh tip per step
    for k, c in enumerate(compounds):
        rig.fill(dil, (k, 0), diluent + t, {c.cid: design["top_conc_nM"] * INTO_STIM},
                 what=f"{c.cid} at {design['top_conc_nM'] * INTO_STIM:g} nM")
        for j in range(1, n):
            rig.fill(dil, (k, j), diluent, what="diluent")
        for j in range(1, n):
            rig.transfer(dil, (k, j - 1), dil, (k, j), t, new_tip=True, mix_after=design["mix"])

    # assay plate: cells, controls, standards, agonist stamps from low to high, then forskolin
    cell_wells = []
    for w, role, _cid, _conc, std in layout:
        if role == "std":
            _add_standard(rig, plate, w, std)
            continue
        cell_wells.append(w)
        _add_cells(rig, plate, w, design["cells_k"], design["ibmx"])
        if role == "neg":
            rig.fill(plate, w, STAMP_UL, what="vehicle")
        elif role == "pos":
            rig.fill(plate, w, STAMP_UL, {"ref": 1.0}, what="reference agonist")
    for k, c in enumerate(compounds):
        for j in reversed(range(n)):
            for w in cpd_wells[(c.cid, j)]:
                rig.transfer(dil, (k, j), plate, w, STAMP_UL, new_tip=True)
    for w in cell_wells:
        if design["forskolin_uM"] > 0:
            rig.fill(plate, w, FSK_UL, {"fsk": design["forskolin_uM"] * FSK_INTO_STIM},
                     what=f"forskolin {design['forskolin_uM'] * FSK_INTO_STIM:g} uM")
        else:
            rig.fill(plate, w, FSK_UL, what="stimulation buffer")

    truth = {c.cid: c for c in compounds}
    roles = {w: (role, cid) for w, role, cid, _conc, _std in layout}

    def inhibition_of(w, volume: float) -> float:
        role, cid = roles[w]
        if role == "neg":
            return 0.0
        if role == "pos":
            return assay.gi_max
        c = truth[cid]
        return c.emax * float(hill(plate.a(cid)[w] / volume, 10.0 ** (9.0 - c.pec50), c.hill))

    ratio = _stimulate_and_read(rig, rng, assay, plate, [w for w, *_ in layout], cell_wells, design["stim_min"],
                                inhibition_of)
    wells = [{"well": well_name(w), "role": role, "id": cid, "conc_nM": conc_nm, "std_nM": std,
              "ratio": round(float(ratio[w]), 1)} for w, role, cid, conc_nm, std in layout]
    return wells, rig


class GiCampPotency(Env):
    name = "gi_camp_potency"
    salted = True  # the hidden truth is not shared with any other environment on the same seed
    cost_reference = 300.0
    MAX_PLATES = 3

    def __init__(self, seed: int = 0, rig="nominal", n_compounds: int = 4):
        super().__init__(seed, rig)
        if not 1 <= n_compounds <= ROWS:
            raise ValueError(f"n_compounds must be between 1 and {ROWS}")
        self.assay = draw_gi_assay(self._truth_rng)
        self.compounds = [draw_agonist(self._truth_rng, f"cpd{i + 1}", self.assay.gi_max) for i in range(n_compounds)]
        self.plates = 0
        self.potency_plates: list[list[dict]] = []

    def spec(self) -> dict:
        return {
            "env": self.name,
            "task": ("Measure the pEC50 of each agonist at a Gi-coupled receptor. Forskolin raises cAMP and the "
                     "agonist lowers it. cAMP is read by a competitive TR-FRET immunoassay (the ratio falls as cAMP "
                     "rises, so it RISES as the agonist acts). Run up to "
                     f"{self.MAX_PLATES} plates, forskolin titrations or potency plates in any order, then report."),
            "actions": {
                "run_titration": {"design": spec_help(TITRATION_SPEC)},
                "run_plate": {"design": spec_help(DESIGN_SPEC)},
                "report": {"pEC50": "object mapping compound id to pEC50 (molar), each a number between "
                                    f"{ANSWER_RANGE[0]:g} and {ANSWER_RANGE[1]:g}; a compound left out scores zero"},
            },
            "labware": {"cells_ul": CELLS_UL, "stamp_ul": STAMP_UL, "forskolin_ul": FSK_UL, "detection_ul": DETECT_UL,
                        "ladder_well_max_ul": LADDER_MAX_UL, "nominal_dead_volume_ul": RigProfile().dead_volume_ul,
                        "ibmx_uM": IBMX_UM, "forskolin_series_uM": list(FSK_SERIES_UM),
                        "standard_series_nM": list(STD_SERIES_NM),
                        "fixed": "cells are plated the day before at cells_k per well, so the outer ring responds "
                                 "less; fresh tip for every transfer; lid on during the stimulation; wells fill the "
                                 "inner 60 first and spill onto the outer ring after that",
                        "roles": {"std": "cAMP standard, no cells", "fsk": "titration well at forskolin_uM",
                                  "unstim": "titration well without forskolin", "neg": "forskolin plus vehicle",
                                  "pos": "forskolin plus a saturating reference agonist",
                                  "cpd": "forskolin plus the agonist at conc_nM"}},
            "reward": f"accuracy * gate - {COST_WEIGHT} * cost / {self.cost_reference:g}. "
                      f"Accuracy: mean over compounds of max(0, 1 - |pEC50 error| / {TOLERANCE}), where the error is "
                      f"the reported pEC50 minus the true one; the tolerance is {TOLERANCE} log units and a compound "
                      "left out of the report scores zero. "
                      "Analysis: none is run for the agent; the reported numbers are scored as given, so turning "
                      "ratios into a pEC50 is part of the task. "
                      "Gate: 0.5, so accuracy is halved, unless at least one potency plate was run and every "
                      f"potency plate run has at least {MIN_CONTROLS} controls a side (pos and neg wells) and control "
                      f"Z' >= {ZPRIME_GATE} on the raw ratio; then it is 1. Z' = 1 - 3 * (sd of pos + sd of neg) / "
                      "|mean of pos - mean of neg| with sample standard deviations. A titration plate is not a "
                      "potency plate and has no gate. "
                      "Cost: 10 per plate + 1 per well + 0.1 per tip, summed over every plate run, titrations "
                      f"included. {self.MAX_INVALID} invalid actions end the episode with the cost term alone.",
        }

    def reset(self) -> dict:
        return {"phase": "start", "compounds": [c.cid for c in self.compounds], "plates_left": self.MAX_PLATES,
                "spec": self.spec()}

    def step(self, action: dict):
        err = self._check_open(action)
        if err:
            return self._reject(err)
        kind = action["type"]
        if kind in ("run_titration", "run_plate"):
            extra = sorted(str(k) for k in set(action) - {"type", "design"})
            if extra:
                return self._reject(f"unknown field: {extra[0]}")
            if self.plates >= self.MAX_PLATES:
                return self._reject("no plates left; report")
            if kind == "run_titration":
                design, err = check_titration(action.get("design"))
                if err:
                    return self._reject(err)
                wells, rig = simulate_titration(self.profile, self._noise, self.assay, design)
            else:
                design, err = check_design(action.get("design"), len(self.compounds))
                if err:
                    return self._reject(err)
                wells, rig = simulate_potency(self.profile, self._noise, self.assay, self.compounds, design)
                self.potency_plates.append(wells)
            self.plates += 1
            self._record(rig.log)
            self.cost += plate_cost(len(wells), rig.tips)
            # the observation gets its own copy: editing it must not reach the plates the gate is scored on
            return ({"phase": "data", "plate": self.plates, "kind": "titration" if kind == "run_titration" else "potency",
                     "wells": copy.deepcopy(wells), "cost": round(self.cost, 2),
                     "plates_left": self.MAX_PLATES - self.plates}, 0.0, False, {})
        if kind == "report":
            return self._score(action)
        return self._reject("unknown action type")

    def _score(self, action: dict):
        extra = sorted(str(k) for k in set(action) - {"type", "pEC50"})
        if extra:
            return self._reject(f"unknown field: {extra[0]}")
        answers = action.get("pEC50")
        if not isinstance(answers, dict):
            return self._reject("the report needs an object mapping compound id to a number")
        ids = {c.cid for c in self.compounds}
        lo, hi = ANSWER_RANGE
        for key, value in answers.items():
            if key not in ids:
                return self._reject("the report names a compound that is not on the plate")
            if isinstance(value, bool) or not isinstance(value, (int, float, np.integer, np.floating)) \
                    or not lo <= value <= hi:  # also rejects NaN, infinity and absurd integers before any conversion
                return self._reject(f"the answer for {key} must be a number between {lo:g} and {hi:g}")
        z = [zprime([w["ratio"] for w in wells if w["role"] == "pos"], [w["ratio"] for w in wells if w["role"] == "neg"])
             for wells in self.potency_plates]
        enough = all(sum(w["role"] == side for w in wells) >= MIN_CONTROLS
                     for wells in self.potency_plates for side in ("pos", "neg"))
        qc_pass = bool(z) and enough and all(np.isfinite(v) and v >= ZPRIME_GATE for v in z)
        errors, scores = {}, []
        for c in self.compounds:
            est = float(answers[c.cid]) if c.cid in answers else None
            scores.append(log_error_score(est, c.pec50, TOLERANCE))
            errors[c.cid] = None if est is None else round(est - c.pec50, 3)
        accuracy = float(np.mean(scores))
        reward = accuracy * (1.0 if qc_pass else 0.5) - self._cost_penalty()
        self.done = True
        a = self.assay
        info = {"accuracy": accuracy, "qc_pass": qc_pass, "zprime": z, "cost": self.cost, "errors": errors,
                "truth": {"pEC50": {c.cid: round(c.pec50, 3) for c in self.compounds},
                          "emax": {c.cid: round(c.emax, 3) for c in self.compounds},
                          "forskolin_ec50_uM": round(a.fsk_ec50, 3), "camp_fmol_per_k_cells_min": round(a.rate, 3),
                          "pde_per_min": round(a.pde, 3), "basal_frac": round(a.basal_frac, 4),
                          "reference_emax": round(a.gi_max, 3), "kit_k_nM": round(a.kit.k, 3)}}
        return {"phase": "done"}, reward, True, info
