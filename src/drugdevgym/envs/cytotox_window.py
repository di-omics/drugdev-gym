"""Cytotoxicity-window environment: design a cell viability dose-response and report each compound's CC50.

Safety, part two: the in-vitro safety window. The functional potency (EC50) of every compound is
given. The agent designs a viability dose-response read by ATP luminescence (the signal falls as
cells die), the environment runs it on the simulated rig, and the agent reports a CC50 per compound,
either as a measured value or as an honest lower bound ("at least X") when no loss of viability
was seen. The window is log10(CC50 / EC50).

Hidden per compound: the CC50 after 72 h (most are above the 100 uM dose limit), whether it acts
fast or slowly, a solubility limit, and for a minority direct inhibition of the luciferase readout.
The reward gives credit for a correct number and for a bound that is true, was tested, and is close
to the highest value the assay could support. It takes credit away for a false alarm (a liability
that is not there), for a missed liability (calling a compound clean above its real CC50, or above
what could be dissolved) and for a bound above what the run tested. It is halved when a plate fails
its QC gate, minus a small cost term.

The QC gate has two terms. One is measured: the control Z' of every plate. The other is a
convention, the one a cytotoxicity guideline asks of a solvent control: the vehicle wells of a
plate must carry the DMSO percentage of the top dose on that plate. A run whose vehicle wells hold
no DMSO fails the gate by rule, whatever its data look like. That is how the unmatched-DMSO mistake
is punished here. The simulated effect of the mismatch (about 11% of the signal in the top well,
at 1% DMSO) is too small for a CC50 to show, and the unmatched route saves tips; see
docs/ASSUMPTIONS.md for the numbers with the convention switched off.
"""
from __future__ import annotations

import copy
from dataclasses import dataclass

import numpy as np

from ..analysis import zprime
from ..chem import dissolved, grown, haber_shift, survival
from ..rig import EDGE, ROWS, Plate, Rig, RigProfile, usable_wells, well_name
from .base import COST_WEIGHT, Env, spec_help, validate

# ---- labware and protocol (fixed; printed in spec()["labware"])
STOCK_MM = 10.0  # compound stock in neat DMSO
CELLS_UL = 90.0  # cells in medium per assay well, plated before dosing
DOSE_UL = 10.0  # 10x dosing solution stamped into each assay well
FINAL_UL = CELLS_UL + DOSE_UL
INTO_ASSAY = FINAL_UL / DOSE_UL  # the stamp is diluted 10x by the cells
INTER_DMSO_UL = 10.0  # aliquot of a DMSO solution taken into medium to make a 10x dosing solution
INTER_MEDIUM_UL = 90.0  # medium it is diluted into
INTO_DOSING = (INTER_DMSO_UL + INTER_MEDIUM_UL) / INTER_DMSO_UL
MATCHED_DMSO_PCT = 100.0 / (INTO_DOSING * INTO_ASSAY)  # 1%: every well of a matched run, the top dose of any run
TOP_LIMIT_UM = STOCK_MM * 1000.0 / (INTO_DOSING * INTO_ASSAY)  # highest dose the stock allows: 100 uM at 1% DMSO
LADDER_KEEP_UL = 45.0  # left in a DMSO ladder well after the transfer on: one 10 uL aliquot, dead volume, spare
MEDIUM_LADDER_KEEP_UL = 70.0  # same for a ladder made in medium, which also feeds up to four stamps
REAGENT_UL = 100.0  # ATP detection reagent, one well volume
CELL_FREE_REF_WELLS = 2  # cell-free ATP wells with vehicle on each plate that carries the interference control
EXPOSURES_H = (24, 48, 72)
REFERENCE_EXPOSURE_H = 72.0  # the exposure the true CC50 is defined at

# ---- cells and readout (assumed; see docs/ASSUMPTIONS.md)
DOUBLING_H = 24.0  # cell doubling time
CONFLUENT_CELLS = 50000.0  # cells a 96-well holds at confluence; growth stops there
CELL_CV = (0.04, 0.08)  # well-to-well cell-number CV, drawn per seed
LUM_PER_CELL = 1.0  # luminescence counts per live cell at the read
LUM_BACKGROUND = 100.0  # counts from reagent and medium with no ATP
ATP_STANDARD_LUM = 20000.0  # counts from the cell-free ATP standard, about 20,000 cells' worth
DMSO_IC50_PCT = 2.0  # DMSO (% v/v) that halves the final cell number
DMSO_HILL = 3.0  # steepness of the DMSO effect: about 11% lost at 1%, 1.5% at 0.5%
OSMOTIC_HALF = 0.4  # fractional rise in medium concentration by evaporation that halves the final cell number
OSMOTIC_HILL = 3.0  # steepness of that effect
SLOW_EXPONENT = 2.0  # slow compounds: CC50(t) = CC50(72 h) * (72 / t) ** 2, so 2.25x at 48 h and 9x at 24 h
TURBIDITY_BLANK = 500.0  # turbidity pre-read of a clear well, in nephelometry units
TURBIDITY_PER_UM = 20.0  # extra turbidity per uM of compound out of solution

# ---- hidden truth per compound (assumed distributions)
LOG_EC50_UM = (-3.0, 0.0)  # given potency: 1 nM to 1 uM
P_LIABILITY = 0.35  # chance that the true CC50 lies inside the 1 to 100 uM dose range
LOG_CC50_LIABLE_UM = (0.0, 2.0)  # a liability: 1 to 100 uM
LOG_CC50_CLEAN_UM = (2.0, 3.5)  # no liability in range: 100 to about 3000 uM
VIABILITY_HILL = (1.0, 2.5)  # slope of the viability curve
P_SLOW = 0.5  # chance that a compound acts slowly
LOG_SOLUBILITY_UM = (1.0, 3.0)  # solubility limit in the assay well: 10 to 1000 uM
P_INTERFERENCE = 0.15  # chance that a compound inhibits the luciferase readout
LOG_LUC_IC50_UM = (0.0, 1.7)  # its IC50 against the readout: 1 to 50 uM
STICKY_MEDIAN, STICKY_SD = 2.0, 0.8  # multiplier on tip carryover, log-normal, as in dose_response

# ---- reward
TOLERANCE = 0.5  # log10 units at which credit reaches zero
PENALTY_FLOOR = -0.5  # lowest score for one compound: a false alarm, a missed liability, an unsupported bound
HEDGE_CAP = 0.5  # the most a bound can score for a compound whose CC50 the run dosed past
ZPRIME_GATE = 0.5  # control Z' every plate must reach
MIN_CONTROLS = 3  # vehicle and full-kill wells a plate needs for a Z'
CC50_REPORT_RANGE_UM = (1e-6, 1e6)  # accepted range of a reported CC50
PLATE_COST, WELL_COST, TIP_COST = 10.0, 1.0, 0.1  # same round numbers as dose_response

DESIGN_SPEC = {
    "top_conc_uM": ("float", 1.0, TOP_LIMIT_UM, "highest final concentration in the assay well"),
    "dilution_factor": ("float", 1.5, 10.0, "fold change between dose points"),
    "n_points": ("int", 4, 12, "dose points per compound"),
    "replicates": ("int", 1, 3, "cell wells per dose point"),
    "exposure_h": ("int", 24, 72, "hours of compound exposure at 37 C: 24, 48 or 72"),
    "cells_k": ("float", 0.5, 40.0, "thousand cells seeded per well"),
    "match_dmso": ("bool", None, None, "true: serial dilution in DMSO, then every point 1:10 into medium, so every "
                                       "well, vehicle included, holds 1% DMSO; false: the top point goes 1:10 into "
                                       "medium and the serial dilution is made in medium, so DMSO falls with the "
                                       "dose and the vehicle wells get none. By convention a plate whose vehicle "
                                       "wells do not carry the DMSO of its top dose fails the QC gate (see reward)"),
    "vehicle_controls": ("int", 0, 16, "vehicle control wells per plate (full viability)"),
    "kill_controls": ("int", 0, 16, "full-kill control wells per plate (no viability)"),
    "interference_control": ("bool", None, None, "also stamp every dose point into one cell-free well holding an "
                                                 "ATP standard, plus two cell-free reference wells per plate"),
    "avoid_edges": ("bool", None, None, "leave the outer ring of every assay plate empty (60 wells instead of 96)"),
    "lid": ("bool", None, None, "lid on the assay plates during the exposure"),
}
ANSWER_SPEC = {
    "cc50_uM": ("float", CC50_REPORT_RANGE_UM[0], CC50_REPORT_RANGE_UM[1], "the CC50, or the bound, in uM"),
    "at_least": ("bool", None, None, "true: a lower bound, the CC50 is at least this; false: a measured CC50"),
}


@dataclass(frozen=True)
class Compound:
    """cc50_um is the true CC50 after 72 h, as dissolved concentration; luc_ic50_um is None without interference."""

    cid: str
    ec50_um: float
    cc50_um: float
    hill: float
    slow: bool
    solubility_um: float
    luc_ic50_um: float | None
    sticky: float

    @property
    def ceiling_um(self) -> float:
        """Highest concentration the assay can expose cells to: the dose limit or the solubility limit."""
        return min(TOP_LIMIT_UM, self.solubility_um)


def draw_compound(rng: np.random.Generator, cid: str) -> Compound:
    ec50 = float(10.0 ** rng.uniform(*LOG_EC50_UM))
    liable, slow, interferes = (bool(rng.random() < p) for p in (P_LIABILITY, P_SLOW, P_INTERFERENCE))
    cc50 = float(10.0 ** rng.uniform(*(LOG_CC50_LIABLE_UM if liable else LOG_CC50_CLEAN_UM)))
    luc = float(10.0 ** rng.uniform(*LOG_LUC_IC50_UM))
    return Compound(cid, ec50_um=ec50, cc50_um=cc50, hill=float(rng.uniform(*VIABILITY_HILL)), slow=slow,
                    solubility_um=float(10.0 ** rng.uniform(*LOG_SOLUBILITY_UM)),
                    luc_ic50_um=luc if interferes else None,
                    sticky=float(np.exp(rng.normal(np.log(STICKY_MEDIAN), STICKY_SD))))


def dmso_survival(pct):
    """Fraction of the final cell number left by DMSO at `pct` percent (v/v)."""
    return survival(pct, DMSO_IC50_PCT, DMSO_HILL)


def osmotic_survival(concentration_ratio):
    """Fraction of the final cell number left when evaporation has concentrated the medium by this ratio."""
    return survival(np.asarray(concentration_ratio, float) - 1.0, OSMOTIC_HALF, OSMOTIC_HILL)


def luciferase_activity(conc_um, ic50_um: float):
    """Fraction of the luminescence left by a compound that inhibits the luciferase reagent."""
    return survival(conc_um, ic50_um, 1.0)


def plan_plates(design: dict, n_compounds: int) -> tuple[list[list[int]] | None, str | None]:
    """Which compounds go on which assay plate. Every plate carries its own controls."""
    capacity = len(usable_wells(design["avoid_edges"]))
    extra = 1 if design["interference_control"] else 0
    fixed = design["vehicle_controls"] + design["kill_controls"] + CELL_FREE_REF_WELLS * extra
    per_compound = design["n_points"] * (design["replicates"] + extra)
    per_plate = min((capacity - fixed) // per_compound, ROWS)
    if per_plate < 1:
        return None, (f"one compound needs {per_compound} wells next to {fixed} control wells "
                      f"but a plate has {capacity}")
    return [list(range(i, min(i + per_plate, n_compounds))) for i in range(0, n_compounds, per_plate)], None


def check_fields(spec: dict, value, what: str) -> tuple[dict | None, str | None]:
    """base.validate, after making sure every key is text (a non-text key cannot come from JSON)."""
    if not isinstance(value, dict) or any(not isinstance(k, str) for k in value):
        return None, f"{what} must be an object with text field names"
    return validate(spec, value)


def check_design(design, n_compounds: int) -> tuple[dict | None, str | None]:
    clean, err = check_fields(DESIGN_SPEC, design, "design")
    if err:
        return None, err
    if clean["exposure_h"] not in EXPOSURES_H:
        return None, "exposure_h must be 24, 48 or 72"
    _, err = plan_plates(clean, n_compounds)
    return (None, err) if err else (clean, None)


def simulate_run(profile: RigProfile, rng: np.random.Generator, compounds: list[Compound], cell_cv: float,
                 design: dict) -> tuple[list[dict], Rig, int]:
    """Run one validated design on the rig. Returns the observed wells, the rig (log, tips) and the plate count."""
    rig = Rig(profile, rng, sticky={c.cid: c.sticky for c in compounds})
    n, reps, df, top = design["n_points"], design["replicates"], design["dilution_factor"], design["top_conc_uM"]
    hours, cells_k, matched = float(design["exposure_h"]), design["cells_k"], design["match_dmso"]
    cell_free = design["interference_control"]
    doses = [top / df ** j for j in range(n)]
    stock = top * INTO_ASSAY * INTO_DOSING  # the DMSO solution both protocols start from: the neat stock for 100 uM
    plates, err = plan_plates(design, len(compounds))
    assert err is None, err
    wells: list[dict] = []
    for p, members in enumerate(plates):
        assay, dosing = Plate(f"assay{p + 1}"), Plate(f"dosing{p + 1}")
        free = iter(usable_wells(design["avoid_edges"]))
        layout: list[tuple] = []  # (well, role, compound index, dose index)
        for i in range(max(design["vehicle_controls"], design["kill_controls"])):
            if i < design["vehicle_controls"]:
                layout.append((next(free), "vehicle", None, None))
            if i < design["kill_controls"]:
                layout.append((next(free), "kill", None, None))
        if cell_free:
            layout += [(next(free), "cell_free_ref", None, None) for _ in range(CELL_FREE_REF_WELLS)]
        for k in members:
            layout += [(next(free), "cpd", k, j) for j in range(n) for _ in range(reps)]
            if cell_free:
                layout += [(next(free), "cell_free", k, j) for j in range(n)]

        # 10x dosing solutions, one row of the dosing plate per compound, high to low
        ladder = Plate(f"dmso{p + 1}")
        for row, k in enumerate(members):
            c = compounds[k]
            if matched:  # serial dilution in neat DMSO, then every point 1:10 into medium
                t = LADDER_KEEP_UL / (df - 1.0)
                rig.fill(ladder, (row, 0), t * df, {c.cid: stock, "DMSO": 100.0},
                         what=f"{c.cid} at {stock:g} uM in DMSO")
                for j in range(1, n):
                    rig.fill(ladder, (row, j), t * (df - 1.0), {"DMSO": 100.0}, what="DMSO")
                for j in range(1, n):
                    rig.transfer(ladder, (row, j - 1), ladder, (row, j), t, new_tip=True, mix_after=True)
                for j in range(n):
                    rig.fill(dosing, (row, j), INTER_MEDIUM_UL, what="medium")
                    rig.transfer(ladder, (row, j), dosing, (row, j), INTER_DMSO_UL, new_tip=True, mix_after=True)
            else:  # the top point 1:10 into medium, then serial dilution in medium: DMSO is diluted with the compound
                t = MEDIUM_LADDER_KEEP_UL / (df - 1.0)
                x = t * df / INTO_DOSING
                rig.fill(ladder, (row, 0), x + LADDER_KEEP_UL, {c.cid: stock, "DMSO": 100.0},
                         what=f"{c.cid} at {stock:g} uM in DMSO")
                rig.fill(dosing, (row, 0), t * df - x, what="medium")
                rig.transfer(ladder, (row, 0), dosing, (row, 0), x, new_tip=True, mix_after=True)
                for j in range(1, n):
                    rig.fill(dosing, (row, j), t * (df - 1.0), what="medium")
                for j in range(1, n):
                    rig.transfer(dosing, (row, j - 1), dosing, (row, j), t, new_tip=True, mix_after=True)

        # assay plate: cells (or the cell-free ATP standard), controls, then the compound stamps
        cellf: dict[tuple[int, int], float] = {}
        for w, role, _k, _j in layout:
            if role in ("cell_free", "cell_free_ref"):
                rig.fill(assay, w, CELLS_UL, {"ATP": 1.0}, what="ATP standard in medium, no cells")
            else:
                rig.fill(assay, w, CELLS_UL, what=f"cells, {cells_k:g}k per well")
                edge = 1.0 - profile.edge_cell_loss if EDGE[w] else 1.0
                cellf[w] = max(0.0, 1.0 + rng.normal(0.0, cell_cv)) * edge
        vehicle = {"DMSO": MATCHED_DMSO_PCT * INTO_ASSAY} if matched else None
        vehicle_name = f"medium with {MATCHED_DMSO_PCT * INTO_ASSAY:g}% DMSO" if matched else "medium"
        for w, role, _k, _j in layout:
            if role in ("vehicle", "cell_free_ref"):
                rig.fill(assay, w, DOSE_UL, vehicle, what=f"vehicle: {vehicle_name}")
            elif role == "kill":
                rig.fill(assay, w, DOSE_UL, vehicle, what=f"full-kill control: lysis agent in {vehicle_name}")
        for row, k in enumerate(members):
            for w, _role, kk, j in layout:
                if kk == k:
                    rig.transfer(dosing, (row, j), assay, w, DOSE_UL, new_tip=True, mix_after=True)

        # turbidity pre-read right after dosing: compound above its solubility limit scatters light
        vol0 = assay.vol.copy()
        safe0 = np.maximum(vol0, 1e-9)
        cloud = np.full(vol0.shape, TURBIDITY_BLANK)
        for k in members:
            c = compounds[k]
            total = assay.a(c.cid) / safe0
            cloud += TURBIDITY_PER_UM * (total - dissolved(total, c.solubility_um))
        turbidity = rig.read(cloud, plate=assay.name)

        rig.incubate(assay, hours * 60.0, design["lid"])
        mean_vol = np.maximum((vol0 + assay.vol) / 2.0, 1e-9)  # what the cells saw on average
        signal = np.zeros(vol0.shape)
        for w, role, k, _j in layout:
            alive, readout = 1.0, 1.0
            if k is not None:
                c = compounds[k]
                in_solution = float(dissolved(assay.a(c.cid)[w] / mean_vol[w], c.solubility_um))
                shift = haber_shift(hours, REFERENCE_EXPOSURE_H, SLOW_EXPONENT) if c.slow else 1.0
                alive = float(survival(in_solution, c.cc50_um * shift, c.hill))
                if c.luc_ic50_um is not None:
                    readout = float(luciferase_activity(in_solution, c.luc_ic50_um))
            if role in ("cell_free", "cell_free_ref"):
                atp = ATP_STANDARD_LUM * assay.a("ATP")[w] / CELLS_UL
            else:
                if role == "kill":
                    alive = 0.0
                alive *= float(dmso_survival(assay.a("DMSO")[w] / mean_vol[w]))
                alive *= float(osmotic_survival(vol0[w] / mean_vol[w]))
                seeded = cells_k * 1000.0 * cellf[w]
                atp = LUM_PER_CELL * float(grown(seeded * alive, hours, DOUBLING_H, CONFLUENT_CELLS))
            signal[w] = atp * readout + LUM_BACKGROUND
        for w, *_ in layout:
            rig.fill(assay, w, REAGENT_UL, what="ATP detection reagent")
        lum = rig.read(signal, plate=assay.name)
        for w, role, k, j in layout:
            dosed = k is not None
            # nominal DMSO: 1% everywhere when matched; otherwise diluted with the dose, and none in the controls
            dmso = MATCHED_DMSO_PCT if matched else (MATCHED_DMSO_PCT / df ** j if dosed else 0.0)
            wells.append({"plate": p + 1, "well": well_name(w), "role": role,
                          "id": compounds[k].cid if dosed else None, "conc_uM": doses[j] if dosed else None,
                          "dmso_pct": round(float(dmso), 4), "lum": round(float(lum[w]), 1),
                          "turbidity": round(float(turbidity[w]), 1)})
    return wells, rig, len(plates)


def run_cost(n_plates: int, n_wells: int, tips: int) -> float:
    """Assay plates, assay wells and tips. The DMSO ladder plates and the dosing plates are not charged."""
    return PLATE_COST * n_plates + WELL_COST * n_wells + TIP_COST * tips


def tips_per_compound(design: dict) -> int:
    """Tips one compound takes: one per transfer. Stated in spec()["reward"] so cost can be planned."""
    n = design["n_points"]
    ladder = 2 * n - 1 if design["match_dmso"] else n  # DMSO ladder plus every point into medium, or one medium ladder
    return ladder + n * (design["replicates"] + (1 if design["interference_control"] else 0))


def score_claim(cc50_um: float, at_least: bool, true_cc50_um: float, ceiling_um: float,
                tested_um: float | None) -> float:
    """Score one reported CC50 against the truth. See spec()["reward"] for the rule in words.

    tested_um is the highest nominal dose the run put on cells; None when nothing was run.
    """
    said, truth = float(np.log10(cc50_um)), float(np.log10(true_cc50_um))
    if not at_least:  # a measured value is a claim about the true CC50
        return float(max(PENALTY_FLOOR, 1.0 - abs(said - truth) / TOLERANCE))
    if tested_um is None:  # nothing was tested, so no bound is supported
        return float(PENALTY_FLOOR)
    ceiling, tested = float(np.log10(ceiling_um)), float(np.log10(tested_um))
    best = min(truth, ceiling)  # the highest bound that is both true and within the ceiling
    supported = min(best, tested)  # the highest bound this run supports
    at_supported = max(0.0, 1.0 - (best - supported) / TOLERANCE)
    if said <= supported:
        score = max(0.0, 1.0 - (best - said) / TOLERANCE)
    else:  # false, above the ceiling or above what was tested: the credit at `supported`, less the excess
        score = max(PENALTY_FLOOR, at_supported - (said - supported) / TOLERANCE)
    if truth < min(ceiling, tested):  # the run dosed past the true CC50: the compound did lose viability
        score = min(score, HEDGE_CAP)
    return float(score)


def vehicle_is_matched(on_plate: list[dict]) -> bool:
    """The solvent-control convention for one plate, read from its plate map (the nominal dmso_pct of each well):
    every vehicle well carries the DMSO percentage of the plate's top dose. False without vehicle or dosed wells.

    A convention, not a measurement: nothing here looks at the luminescence.
    """
    vehicle = [w["dmso_pct"] for w in on_plate if w["role"] == "vehicle"]
    dosed = [w for w in on_plate if w["role"] == "cpd"]
    if not vehicle or not dosed:
        return False
    top = max(w["conc_uM"] for w in dosed)
    top_dmso = max(w["dmso_pct"] for w in dosed if w["conc_uM"] == top)
    return all(v == top_dmso for v in vehicle)


def gate(wells: list[dict] | None) -> tuple[bool, list[float | None], list[bool]]:
    """QC gate from the run's own data, per plate: the control Z', whether the vehicle is solvent-matched,
    and whether every plate passes both."""
    if not wells:
        return False, [], []
    values: list[float | None] = []
    matched: list[bool] = []
    ok = True
    for plate in sorted({w["plate"] for w in wells}):
        on_plate = [w for w in wells if w["plate"] == plate]
        vehicle = [w["lum"] for w in on_plate if w["role"] == "vehicle"]
        kill = [w["lum"] for w in on_plate if w["role"] == "kill"]
        z = zprime(vehicle, kill) if min(len(vehicle), len(kill)) >= MIN_CONTROLS else float("nan")
        values.append(round(float(z), 3) if np.isfinite(z) else None)
        matched.append(bool(vehicle_is_matched(on_plate)))
        ok = ok and bool(np.isfinite(z) and z >= ZPRIME_GATE) and matched[-1]
    return ok, values, matched


class CytotoxWindow(Env):
    name = "cytotox_window"
    salted = True
    cost_reference = 600.0

    def __init__(self, seed: int = 0, rig="nominal", n_compounds: int = 6):
        super().__init__(seed, rig)
        if not 1 <= n_compounds <= ROWS:
            raise ValueError(f"n_compounds must be between 1 and {ROWS}")
        self.cell_cv = float(self._truth_rng.uniform(*CELL_CV))
        self.compounds = [draw_compound(self._truth_rng, f"cpd{i + 1}") for i in range(n_compounds)]
        self.wells: list[dict] | None = None
        self.tested_um: float | None = None  # highest nominal dose the run put on cells

    def spec(self) -> dict:
        return {
            "env": self.name,
            "task": ("Measure general cytotoxicity for each compound in a cell viability assay read by ATP "
                     "luminescence (the signal falls as cells die) and report its CC50 in uM. The functional EC50 "
                     "of each compound is given, and the safety window is log10(CC50 / EC50). The CC50 that counts "
                     "is the one after 72 h of exposure; some compounds act slowly and show little at shorter "
                     "exposures. Most compounds are not toxic up to the 100 uM dose limit: for those the right "
                     "answer is a lower bound, not a number. Every compound has an unknown solubility limit, above "
                     "which the nominal concentration is not reached, and a few inhibit the luciferase readout "
                     "itself, which looks like lost viability. DMSO lowers viability above about 1%, and a solvent "
                     "control has to carry the solvent of the dosed wells. Design one run (as many plates as it "
                     "needs), then report."),
            "actions": {
                "run": {"design": spec_help(DESIGN_SPEC)},
                "report": {"results": {"<compound id>": spec_help(ANSWER_SPEC)},
                           "unanswered": "leave a compound out of results to give no answer for it"},
            },
            "observation": {
                "compounds": "id and given ec50_uM of each compound",
                "wells": ("one row per assay well: plate, well, role (vehicle, kill, cpd, cell_free, cell_free_ref), "
                          "compound id, nominal conc_uM, nominal dmso_pct, lum (luminescence at the end of the "
                          "exposure) and turbidity (a pre-read taken right after dosing; it rises when a well "
                          "holds compound out of solution)"),
            },
            "labware": {"stock_mM": STOCK_MM, "top_limit_uM": TOP_LIMIT_UM, "cells_ul": CELLS_UL, "dose_ul": DOSE_UL,
                        "matched_dmso_pct": MATCHED_DMSO_PCT, "exposures_h": list(EXPOSURES_H),
                        "doubling_time_h": DOUBLING_H, "confluent_cells_per_well": CONFLUENT_CELLS,
                        "cell_free_reference_wells": CELL_FREE_REF_WELLS,
                        "nominal_dead_volume_ul": RigProfile().dead_volume_ul,
                        "layout": "wells fill column by column: controls first (vehicle and full-kill alternating), "
                                  "then the cell-free reference wells, then each compound from its top dose down, "
                                  "followed by its cell-free wells"},
            "reward": (
                f"Per compound, with L, C, T and X each the log10 of a concentration in uM. L is the true CC50 after "
                f"72 h. C is the assay ceiling: the lower of {TOP_LIMIT_UM:g} uM (the dose limit, the highest dose "
                f"the stock reaches at {MATCHED_DMSO_PCT:g}% DMSO) and the compound's solubility limit. T is the "
                f"highest nominal dose the run put on cells (the design's top_conc_uM). X is the reported cc50_uM. "
                f"A measured value (at_least false) scores max({PENALTY_FLOOR:g}, 1 - |X - L| / {TOLERANCE:g}): 1 when "
                f"exact, 0 when {TOLERANCE:g} log units off, and below 0 for a false alarm or an understated "
                f"liability. A lower bound (at_least true) is scored against B = min(L, C), the highest bound that "
                f"is both true and within the ceiling, and S = min(B, T), the highest bound this run supports. At or "
                f"below S it scores max(0, 1 - (B - X) / {TOLERANCE:g}), so a supported bound earns credit only when "
                f"it is close to B and a very safe bound earns nothing. Above S it scores max({PENALTY_FLOOR:g}, "
                f"max(0, 1 - (B - S) / {TOLERANCE:g}) - (X - S) / {TOLERANCE:g}): the credit of the bound at S, less "
                f"the excess, so calling a compound clean above its real CC50, above what could be dissolved or "
                f"above the highest dose tested is penalised. With no run nothing was tested and "
                f"every bound scores {PENALTY_FLOOR:g}. When the run dosed past the true CC50 (L < min(C, T)) the "
                f"compound did lose viability, and a bound for it scores at most {HEDGE_CAP:g}. A compound left out "
                f"scores 0. The EC50 is given, so the error in the window log10(CC50 / EC50) is the same number. "
                f"Quality is the mean over compounds. QC gate, two terms per plate. Measured: at least "
                f"{MIN_CONTROLS} vehicle and {MIN_CONTROLS} full-kill wells and a Z' of at least {ZPRIME_GATE:g} "
                f"between their luminescence. By convention: a solvent-matched vehicle, meaning the vehicle wells "
                f"carry the DMSO percentage of the top dose on that plate (their nominal dmso_pct equals that of "
                f"the plate's highest-dose wells). This is what a cytotoxicity guideline asks of a solvent control. "
                f"It is a rule, not a simulated effect: a run with match_dmso false leaves its vehicle wells without "
                f"DMSO and cannot pass the gate, whatever its data look like, although at {MATCHED_DMSO_PCT:g}% DMSO "
                f"the mismatch itself changes the answers very little. "
                f"When the gate fails, or nothing was run, a quality above 0 is halved (a quality below 0 is kept). "
                f"Reward = that quality minus {COST_WEIGHT:g} * cost / {self.cost_reference:g}, where cost is "
                f"{PLATE_COST:g} per assay plate, {WELL_COST:g} per assay well (controls and cell-free wells "
                f"included) and {TIP_COST:g} per tip; the DMSO ladder plates and the dosing plates are not charged. "
                f"Tips: one per transfer. Before the stamps a matched run takes 2 * n_points - 1 tips per compound "
                f"(the DMSO ladder, then every point into medium) and an unmatched run n_points; both then take one "
                f"tip per dosed assay well, cell-free wells included. "
                f"Changed after the first card, following two adversarial reviews: T and S (a bound is scored "
                f"against what the run tested), the cap of {HEDGE_CAP:g} on a bound for a compound that lost "
                f"viability, and the solvent-matched vehicle term of the gate are new. A second pass had added "
                f"untreated wells, a vehicle-health term of the gate and doses up to 300 uM; once the "
                f"solvent-matched vehicle term was in they were measured to do no job and were removed, so the "
                f"dose limit is {TOP_LIMIT_UM:g} uM as on the first card. The scripted reference design changed: "
                f"one well per dose and 6 dose points instead of 2 and 8; the 5 untreated wells per plate it "
                f"carried between the second and the third pass are gone. See docs/ASSUMPTIONS.md."),
        }

    def reset(self) -> dict:
        return {"phase": "start", "compounds": [{"id": c.cid, "ec50_uM": c.ec50_um} for c in self.compounds],
                "spec": self.spec()}

    def step(self, action: dict):
        err = self._check_open(action)
        if err:
            return self._reject(err)
        if action["type"] == "run":
            if set(action) - {"type", "design"}:
                return self._reject("a run takes one field: design")
            if self.wells is not None:
                return self._reject("the plates have been run; report")
            design, err = check_design(action.get("design"), len(self.compounds))
            if err:
                return self._reject(err)
            self.wells, rig, n_plates = simulate_run(self.profile, self._noise, self.compounds, self.cell_cv, design)
            self.tested_um = design["top_conc_uM"]
            self._record(rig.log)
            self.cost += run_cost(n_plates, len(self.wells), rig.tips)
            return ({"phase": "data", "plates": n_plates, "wells": copy.deepcopy(self.wells),
                     "cost": round(self.cost, 2)}, 0.0, False, {})
        if action["type"] == "report":
            if set(action) - {"type", "results"}:
                return self._reject("a report takes one field: results")
            return self._score(action.get("results"))
        return self._reject(f"unknown action type: {action['type']}")

    def _score(self, results):
        if not isinstance(results, dict):
            return self._reject("results must be an object mapping compound id to an answer")
        known = {c.cid for c in self.compounds}
        answers: dict[str, dict] = {}
        for cid, answer in results.items():
            if not isinstance(cid, str) or cid not in known:
                return self._reject("results names a compound that is not in this task")
            clean, err = check_fields(ANSWER_SPEC, answer, "an answer")
            if err:
                return self._reject(f"{cid}: {err}")
            answers[cid] = clean
        scores, reported, truth = {}, {}, {}
        for c in self.compounds:
            a = answers.get(c.cid)
            scores[c.cid] = 0.0 if a is None else score_claim(a["cc50_uM"], a["at_least"], c.cc50_um, c.ceiling_um,
                                                              self.tested_um)
            reported[c.cid] = None if a is None else {"cc50_uM": a["cc50_uM"], "at_least": a["at_least"],
                                                      "window_log10": round(float(np.log10(a["cc50_uM"] / c.ec50_um)), 3)}
            truth[c.cid] = {"cc50_uM": c.cc50_um, "window_log10": float(np.log10(c.cc50_um / c.ec50_um)),
                            "ceiling_uM": c.ceiling_um, "solubility_uM": c.solubility_um,
                            "kinetics": "slow" if c.slow else "fast", "viability_hill": c.hill,
                            "luciferase_ic50_uM": c.luc_ic50_um}
        quality = float(np.mean(list(scores.values())))
        qc_pass, z, matched = gate(self.wells)
        gated = quality if qc_pass or quality <= 0.0 else 0.5 * quality
        reward = gated - self._cost_penalty()
        self.done = True
        info = {"quality": quality, "qc_pass": qc_pass, "zprime": z, "vehicle_matched": matched,
                "tested_top_uM": self.tested_um, "cost": self.cost,
                "scores": {k: round(v, 3) for k, v in scores.items()}, "reported": reported, "truth": truth}
        return {"phase": "done"}, reward, True, info
