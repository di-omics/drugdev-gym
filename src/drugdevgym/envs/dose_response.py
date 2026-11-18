"""Dose-response environment: design a serial dilution on a cell assay plate and report pEC50s.

The agent chooses the protocol; the environment runs it on the simulated rig. The reward is the
accuracy of the reported pEC50s against the hidden truth, halved when the plate fails its
own QC (control Z' below 0.5), minus a small cost term.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..analysis import zprime
from ..chem import Kit, hill
from ..rig import EDGE, ROWS, Plate, Rig, RigProfile, usable_wells, well_name
from .base import Env, log_error_score, spec_help, validate

STAMP_UL = 10.0
CELLS_UL = 40.0
DETECT_UL = 20.0
FINAL_UL = CELLS_UL + STAMP_UL + DETECT_UL
INTO_ASSAY = (CELLS_UL + STAMP_UL) / STAMP_UL  # the stamp is diluted 5x by the cells
LADDER_MAX_UL = 300.0
STD_SERIES_NM = (1.0, 3.0, 10.0, 30.0, 100.0, 300.0, 1000.0, 3000.0)
STIM_MIN = 60.0
ZPRIME_GATE = 0.5

DESIGN_SPEC = {
    "top_conc_nM": ("float", 1.0, 100000.0, "highest final concentration in the assay well"),
    "dilution_factor": ("float", 1.5, 10.0, "fold change between ladder points"),
    "n_points": ("int", 4, 12, "ladder points per compound"),
    "replicates": ("int", 1, 4, "assay wells per ladder point"),
    "transfer_ul": ("float", 2.0, 100.0, "ladder transfer volume; each ladder well holds transfer_ul * dilution_factor"),
    "mix": ("bool", None, None, "mix each ladder well after the transfer into it"),
    "new_tips": ("bool", None, None, "fresh tip for every transfer instead of two tips per compound (one for the ladder, one for stamping)"),
    "avoid_edges": ("bool", None, None, "leave the outer ring of the assay plate empty (60 wells instead of 96)"),
    "lid": ("bool", None, None, "lid on the assay plate during the 37 C stimulation"),
    "standards": ("bool", None, None, "include the 8-point analyte standard curve (8 wells)"),
    "controls_per_side": ("int", 0, 16, "negative and positive control wells, each"),
}


@dataclass(frozen=True)
class Compound:
    cid: str
    pec50: float
    emax: float
    hill: float
    sticky: float


@dataclass(frozen=True)
class CellAssay:
    """basal and top are analyte nM at the nominal final volume, unstimulated and fully stimulated."""

    kit: Kit
    basal: float
    top: float
    cell_cv: float


def draw_assay(rng: np.random.Generator) -> CellAssay:
    k = float(rng.uniform(20.0, 60.0))
    return CellAssay(Kit(k=k), basal=k * float(rng.uniform(0.1, 0.3)), top=k * float(rng.uniform(3.0, 8.0)),
                     cell_cv=float(rng.uniform(0.04, 0.08)))


def draw_compound(rng: np.random.Generator, cid: str) -> Compound:
    return Compound(cid, pec50=float(rng.uniform(5.5, 9.5)), emax=float(rng.uniform(0.6, 1.0)),
                    hill=float(rng.uniform(0.8, 1.3)), sticky=float(np.exp(rng.normal(np.log(2.0), 0.8))))


def check_design(design, n_compounds: int) -> tuple[dict | None, str | None]:
    clean, err = validate(DESIGN_SPEC, design)
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


def simulate_plate(profile: RigProfile, rng: np.random.Generator, assay: CellAssay,
                   compounds: list[Compound], design: dict) -> tuple[list[dict], Rig]:
    """Run one validated design on the rig. Returns the observed wells and the rig (log, tips)."""
    rig = Rig(profile, rng, sticky={c.cid: c.sticky for c in compounds})
    n, reps, t, df = design["n_points"], design["replicates"], design["transfer_ul"], design["dilution_factor"]
    diluent = t * (df - 1.0)
    dil, plate = Plate("dilution"), Plate("assay")
    free = iter(usable_wells(design["avoid_edges"]))
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

    # serial dilution ladders, one row of the dilution plate per compound, high to low
    for k, c in enumerate(compounds):
        rig.fill(dil, (k, 0), diluent + t, {c.cid: design["top_conc_nM"] * INTO_ASSAY},
                 what=f"{c.cid} at {design['top_conc_nM'] * INTO_ASSAY:g} nM")
        for j in range(1, n):
            rig.fill(dil, (k, j), diluent, what="diluent")
        for j in range(1, n):
            rig.transfer(dil, (k, j - 1), dil, (k, j), t, new_tip=design["new_tips"] or j == 1,
                         mix_after=design["mix"])

    # assay plate: cells, controls, standards, then compound stamps from low to high
    cellf: dict[tuple[int, int], float] = {}
    for w, role, *_ in layout:
        if role == "std":
            continue
        rig.fill(plate, w, CELLS_UL, what="cells")
        edge = 1.0 - profile.edge_cell_loss if EDGE[w] else 1.0
        cellf[w] = max(0.0, 1.0 + rng.normal(0.0, assay.cell_cv)) * edge
    for w, role, _cid, _conc, std in layout:
        if role == "neg":
            rig.fill(plate, w, STAMP_UL, what="vehicle")
        elif role == "pos":
            rig.fill(plate, w, STAMP_UL, {"ref": 1.0}, what="reference agonist")
        elif role == "std":
            rig.fill(plate, w, CELLS_UL + STAMP_UL, {"M": std * FINAL_UL / (CELLS_UL + STAMP_UL)},
                     what=f"analyte standard {std:g} nM")
    for k, c in enumerate(compounds):
        first = True
        for j in reversed(range(n)):
            for w in cpd_wells[(c.cid, j)]:
                rig.transfer(dil, (k, j), plate, w, STAMP_UL, new_tip=design["new_tips"] or first)
                first = False

    before = plate.vol.copy()
    rig.incubate(plate, STIM_MIN, design["lid"])
    effective_vol = (before + plate.vol) / 2.0
    truth = {c.cid: c for c in compounds}
    analyte = plate.a("M")
    for w, role, cid, _conc, _std in layout:
        if role == "std":
            continue
        if role == "neg":
            frac = 0.0
        elif role == "pos":
            frac = 1.0
        else:
            c = truth[cid]
            frac = c.emax * float(hill(plate.a(cid)[w] / effective_vol[w], 10.0 ** (9.0 - c.pec50), c.hill))
        analyte[w] = cellf[w] * (assay.basal + (assay.top - assay.basal) * frac) * FINAL_UL
    for w, *_ in layout:
        rig.fill(plate, w, DETECT_UL, what="detection mix")
    conc = np.where(plate.vol > 0, analyte / np.maximum(plate.vol, 1e-9), 0.0)
    ratio = rig.read(assay.kit.ratio(conc))
    wells = [{"well": well_name(w), "role": role, "id": cid, "conc_nM": conc_nm, "std_nM": std,
              "ratio": round(float(ratio[w]), 1)} for w, role, cid, conc_nm, std in layout]
    return wells, rig


def plate_cost(n_wells: int, tips: int) -> float:
    return 10.0 + 1.0 * n_wells + 0.1 * tips


class DoseResponse(Env):
    name = "dose_response"
    cost_reference = 200.0
    MAX_PLATES = 2

    def __init__(self, seed: int = 0, rig="nominal", n_compounds: int = 4):
        super().__init__(seed, rig)
        if not 1 <= n_compounds <= ROWS:
            raise ValueError(f"n_compounds must be between 1 and {ROWS}")
        self.assay = draw_assay(self._truth_rng)
        self.compounds = [draw_compound(self._truth_rng, f"cpd{i + 1}") for i in range(n_compounds)]
        self.plates: list[list[dict]] = []

    def spec(self) -> dict:
        return {
            "env": self.name,
            "task": ("Measure the pEC50 of each compound on a Gq-style cell assay read by a competitive TR-FRET "
                     "immunoassay (the ratio falls as the analyte rises). Run up to "
                     f"{self.MAX_PLATES} plates, then report."),
            "actions": {
                "run_plate": {"design": spec_help(DESIGN_SPEC)},
                "report": {"pEC50": "object mapping compound id to pEC50 (molar)",
                           "plates_used": "optional list of plate numbers the answer rests on (default: all run). "
                                          "Self-declared: the QC gate is evaluated on these plates and nothing "
                                          "checks that the answer was computed from them."},
            },
            "labware": {"stamp_ul": STAMP_UL, "cells_ul": CELLS_UL, "ladder_well_max_ul": LADDER_MAX_UL,
                        "nominal_dead_volume_ul": RigProfile().dead_volume_ul, "stimulation_min": STIM_MIN,
                        "standard_series_nM": list(STD_SERIES_NM)},
            "reward": "mean over compounds of max(0, 1 - |pEC50 error| / 0.5); halved unless every plate used "
                      "has at least 3 controls a side and control Z' >= 0.5; minus 0.15 * cost / 200",
        }

    def reset(self) -> dict:
        return {"phase": "start", "compounds": [c.cid for c in self.compounds], "plates_left": self.MAX_PLATES,
                "spec": self.spec()}

    def step(self, action: dict):
        err = self._check_open(action)
        if err:
            return self._reject(err)
        if action["type"] == "run_plate":
            if len(self.plates) >= self.MAX_PLATES:
                return self._reject("no plates left; report")
            design, err = check_design(action.get("design"), len(self.compounds))
            if err:
                return self._reject(err)
            wells, rig = simulate_plate(self.profile, self._noise, self.assay, self.compounds, design)
            self.plates.append(wells)
            self._record(rig.log)
            self.cost += plate_cost(len(wells), rig.tips)
            return ({"phase": "data", "plate": len(self.plates), "wells": wells, "cost": round(self.cost, 2),
                     "plates_left": self.MAX_PLATES - len(self.plates)}, 0.0, False, {})
        if action["type"] == "report":
            return self._score(action)
        return self._reject(f"unknown action type: {action['type']}")

    def _score(self, action: dict):
        answers = action.get("pEC50")
        answers = answers if isinstance(answers, dict) else {}
        used = action.get("plates_used")
        if used is None:
            used = list(range(1, len(self.plates) + 1))
        elif not isinstance(used, list) or any(
                isinstance(p, bool) or not isinstance(p, (int, float)) or not 1 <= p <= len(self.plates) or p != int(p)
                for p in used):
            return self._reject(f"plates_used must list plate numbers between 1 and {len(self.plates)}")
        used = sorted({int(p) for p in used})
        z = [zprime([w["ratio"] for w in self.plates[p - 1] if w["role"] == "pos"],
                    [w["ratio"] for w in self.plates[p - 1] if w["role"] == "neg"]) for p in used]
        qc_pass = bool(used) and all(np.isfinite(v) and v >= ZPRIME_GATE for v in z)
        errors = {}
        scores = []
        for c in self.compounds:
            est = answers.get(c.cid)
            try:
                est = float(est) if isinstance(est, (int, float, np.floating)) and not isinstance(est, bool) else None
            except (OverflowError, ValueError):
                est = None
            scores.append(log_error_score(est, c.pec50))
            errors[c.cid] = None if est is None else round(est - c.pec50, 3)
        accuracy = float(np.mean(scores))
        reward = accuracy * (1.0 if qc_pass else 0.5) - self._cost_penalty()
        self.done = True
        info = {"accuracy": accuracy, "qc_pass": qc_pass, "zprime": z, "cost": self.cost, "errors": errors,
                "truth": {c.cid: round(c.pec50, 3) for c in self.compounds}}
        return {"phase": "done"}, reward, True, info
