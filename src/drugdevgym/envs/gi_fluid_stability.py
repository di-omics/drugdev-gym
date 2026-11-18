"""GI fluid stability environment: digest peptides in simulated gastric and intestinal fluid and report half-lives.

SGF is pepsin at pH 1.2, SIF is pancreatin at pH 6.8, both at 37 C. Each replicate is its own
incubation well. At every time point an aliquot goes into a quench plate, and the quenched
samples wait on that plate and in the injection queue until their LC-MS injection. Three things
make this harder than a plasma time course:

- the enzyme loses activity in storage and during the run. Yesterday's pancreatin has about two
  thirds of its activity, so intestinal half-lives read about 1.5 times too long (pepsin loses
  under 10%), and only a control substrate shows it;
- digestion goes on in a quenched sample unless the quench stops that enzyme: an acid stop
  works for pancreatin and does nothing to pepsin;
- half-lives run from under a minute to days, and some compounds also fall apart in acid
  with no enzyme at all.

The reward is the accuracy of the reported half-lives at specification enzyme activity, on a
log scale. A pair counts at half when the run fails its own QC (time-zero replicate CV and the
activity control) or when that compound was not run in that fluid.
"""
from __future__ import annotations

import copy
from dataclasses import dataclass

import numpy as np

from ..analysis import enzyme_activity
from ..chem import acid_catalysis_factor, enzyme_exposure, q10_factor, remaining
from ..rig import ROWS, Plate, Rig, RigProfile, usable_wells
from .base import COST_WEIGHT, Env, log_error_score, spec_help, validate
from .stability import INJECTION_SD, IS_SD, NOISE_FLOOR, T0_CV_GATE, t0_cv

FLUID_UL = 200.0
SPIKE_UL = 2.0
ALIQUOT_UL = 20.0
QUENCH_UL = 80.0
HORIZON_MIN = 240.0
HL_LO, HL_HI = 1.0, 960.0  # quantifiable range; answers and truth are clipped to it
TOLERANCE = 0.5  # log10 half-life error at which a pair's credit reaches zero
LLOQ_FRAC = 0.02  # limit of quantification of the gate's control fit, as a fraction of the time-zero mean

FLUIDS = ("sgf", "sif")
ENZYME = {"sgf": "pepsin", "sif": "pancreatin"}
FLUID_NAME = {"sgf": "SGF with pepsin", "sif": "SIF with pancreatin"}
FLUID_PH = {"sgf": 1.2, "sif": 6.8}  # USP simulated gastric and intestinal fluid
ACTIVITY_HL_37C_MIN = {"sgf": 1440.0, "sif": 240.0}  # half-life of proteolytic activity at 37 C (self-digestion)
Q10 = 2.0  # every rate here (digestion, enzyme inactivation, acid hydrolysis) halves per 10 C of cooling
STORAGE_TEMP_C = 4.0  # where a preparation made ahead is kept
STALE_AGE_H = 24.0  # age of a preparation that was not made for this run: yesterday's, from the fridge
LOT_ACTIVITY_SD = 0.10  # log-scale spread of a fresh preparation's activity around specification (lot potency, weighing)
HANDLING_MIN = 1.0  # spike, mix and return: earliest later time point, and the real age of a time zero not pre-quenched
INJECTION_MIN = 3.0  # LC-MS cycle time; samples are injected in plate order, the SGF plate first

QUENCHES = ("acid", "base", "organic", "none")
# Fraction of the enzyme's activity left in the quenched sample, before dilution and cooling are counted.
# Pepsin works at low pH and is destroyed above pH 8; pancreatic proteases are the reverse; four volumes of
# acetonitrile precipitate most of either; "none" is diluent only.
RESIDUAL_ACTIVITY = {"sgf": {"acid": 1.0, "base": 0.0, "organic": 0.02, "none": 1.0},
                     "sif": {"acid": 0.01, "base": 1.0, "organic": 0.05, "none": 1.0}}
# pH of the quenched sample, which sets how fast an acid-labile compound keeps hydrolysing while it waits.
# Unbuffered gastric acid diluted five-fold reads about 1.9 whatever non-basic quench is used.
QUENCHED_PH = {"sgf": {"acid": 1.9, "base": 9.0, "organic": 1.9, "none": 1.9},
               "sif": {"acid": 2.5, "base": 9.0, "organic": 6.8, "none": 6.8}}

CONTROL_ID = "control"
CONTROL_HALF_LIFE_MIN = 15.0  # reference half-life of the activity control substrate at specification activity, either fluid
ACTIVITY_GATE = (0.8, 1.25)  # acceptance window for the activity read from the control substrate

ENZYMATIC_HL_RANGE_MIN = {"sgf": (1.0, 3000.0), "sif": (0.5, 1000.0)}  # log-uniform hidden truth at specification activity
ACID_LABILE_FRACTION = 0.3  # share of compounds that also hydrolyse in acid with no enzyme
ACID_HL_RANGE_MIN = (30.0, 600.0)  # their non-enzymatic half-life at pH 1.2 and 37 C, log-uniform
ENZYME_PREP_COST = 20.0  # one fresh enzyme preparation, per fluid
COST_REFERENCE = 500.0

LN2 = float(np.log(2.0))

REWARD_RULE = (
    f"Each compound-fluid pair scores max(0, 1 - |log10 half-life error| / {TOLERANCE:g}): the reported half-life "
    f"against the true one, both clipped to {HL_LO:g} to {HL_HI:g} min, a pair with no answer scoring 0. The true "
    "half-life is the one in that fluid at specification enzyme activity, with any loss that needs no enzyme "
    "counted in. A pair's score is halved unless that compound was run in that fluid and the run passes its QC "
    f"gate. The reward is the mean over every pair, minus {COST_WEIGHT:g} * cost / {COST_REFERENCE:g}, where cost = "
    f"10 + 2 per LC-MS sample + 0.1 per tip + {ENZYME_PREP_COST:g} per fresh enzyme preparation. "
    "The gate passes when both of these hold. (1) Every arm that was run, the control substrate included, has 2+ "
    "replicates sampled at time 0, and their CV (sample SD over mean of area over internal standard, or of area "
    f"alone without one) is <= {T0_CV_GATE:.0%}. (2) In every fluid that was run the control substrate was run too "
    f"and reads an enzyme activity of {ACTIVITY_GATE[0]:g} to {ACTIVITY_GATE[1]:g}; a control that cannot be "
    "fitted fails. Activity is the control's reference half-life over its fitted half-life "
    "(analysis.enzyme_activity): log-linear regression of its response on minutes of digestion at the starting "
    "activity (chem.enzyme_exposure with the listed enzyme activity half-life), after dropping the first time point "
    f"whose replicate mean is at or under {LLOQ_FRAC:.0%} of the time-zero mean and every later one, and any "
    "response that is not above zero. If only time zero is left, the fitted half-life is that first dropped time "
    f"point's digestion time over log2(1 / {LLOQ_FRAC:g})."
)

DESIGN_SPEC = {
    "timepoints_min": ("times", 0.0, HORIZON_MIN, "2 to 8 sampling times in minutes, shared by both fluids; "
                                                   f"a time other than 0 cannot be earlier than {HANDLING_MIN:g}"),
    "replicates": ("int", 1, 3, "incubation wells per compound and fluid"),
    "internal_standard": ("bool", None, None, "internal standard in the quench solution"),
    "fresh_enzyme": ("bool", None, None, "dissolve the enzymes for this run (one preparation charged per fluid) "
                                         f"instead of using the batch made {STALE_AGE_H:g} h ago and kept at "
                                         f"{STORAGE_TEMP_C:g} C (free)"),
    "activity_control": ("bool", None, None, "also run the control substrate in every fluid used"),
    "prequenched_t0": ("bool", None, None, "make time zero by quenching the fluid before the compound meets live "
                                           "enzyme, instead of sampling after the spike and mix"),
    "hold_min": ("float", 0.0, 1440.0, "minutes between the last time point and the first injection"),
    "hold_temp_c": ("float", 4.0, 37.0, "temperature of the quench plate from quench to injection"),
}
CHOICE_SPEC = {
    "sgf_quench": "quench solution for gastric samples",
    "sif_quench": "quench solution for intestinal samples",
}
QUENCH_HELP = {"acid": "formic acid, to about pH 2", "base": "bicarbonate buffer, to about pH 9",
               "organic": "four volumes of acetonitrile", "none": "diluent only"}
FLUIDS_HELP = ("object mapping compound id to the list of fluids to run it in, from: " + ", ".join(FLUIDS)
               + "; a compound left out or given an empty list is not run")


@dataclass(frozen=True)
class Peptide:
    """One compound's hidden kinetics. The enzymatic half-lives are at specification activity."""

    cid: str
    sgf_enzymatic_min: float
    sif_enzymatic_min: float
    acid_min: float | None  # non-enzymatic half-life at pH 1.2, or None when the compound is acid-stable


CONTROL = Peptide(CONTROL_ID, CONTROL_HALF_LIFE_MIN, CONTROL_HALF_LIFE_MIN, None)


def draw_peptide(rng: np.random.Generator, cid: str) -> Peptide:
    def log_uniform(lo: float, hi: float) -> float:
        return float(10 ** rng.uniform(np.log10(lo), np.log10(hi)))

    sgf = log_uniform(*ENZYMATIC_HL_RANGE_MIN["sgf"])
    sif = log_uniform(*ENZYMATIC_HL_RANGE_MIN["sif"])
    labile, acid = rng.random() < ACID_LABILE_FRACTION, log_uniform(*ACID_HL_RANGE_MIN)
    return Peptide(cid, sgf, sif, acid if labile else None)


def enzymatic_rate(p: Peptide, fluid: str) -> float:
    """First-order rate constant (1/min) of digestion at specification activity."""
    return LN2 / (p.sgf_enzymatic_min if fluid == "sgf" else p.sif_enzymatic_min)


def acid_rate(p: Peptide, ph: float) -> float:
    """First-order rate constant (1/min) of non-enzymatic acid hydrolysis at this pH and 37 C."""
    return 0.0 if p.acid_min is None else LN2 / p.acid_min * acid_catalysis_factor(ph, FLUID_PH["sgf"])


def true_half_life(p: Peptide, fluid: str) -> float:
    """Half-life in the fluid at specification enzyme activity: enzymatic and acid loss together."""
    return LN2 / (enzymatic_rate(p, fluid) + acid_rate(p, FLUID_PH[fluid]))


def stored_activity(fluid: str, age_h: float, temp_c: float = STORAGE_TEMP_C) -> float:
    """Fraction of its activity an enzyme preparation keeps after age_h hours at temp_c."""
    return remaining(ACTIVITY_HL_37C_MIN[fluid] / q10_factor(temp_c, 37.0, Q10), 60.0 * age_h)


def _check_fluids(value, compounds: list[str]) -> tuple[dict | None, str | None]:
    if not isinstance(value, dict):
        return None, "fluids must be an object mapping compound id to a list of fluids"
    out: dict[str, list[str]] = {}
    for cid, chosen in value.items():
        if not isinstance(cid, str) or cid not in compounds:
            return None, "fluids names a compound that is not in this task"
        if not isinstance(chosen, (list, tuple)) or any(not isinstance(f, str) or f not in FLUIDS for f in chosen):
            return None, f"fluids for {cid} must be a list drawn from: {', '.join(FLUIDS)}"
        if len(set(chosen)) != len(chosen):
            return None, f"fluids for {cid} lists a fluid twice"
        out[cid] = [f for f in FLUIDS if f in chosen]
    if not any(out.values()):
        return None, "fluids selects nothing to run"
    return out, None


def check_design(design, compounds: list[str]) -> tuple[dict | None, str | None]:
    if not isinstance(design, dict):
        return None, "design must be an object"
    if any(not isinstance(k, str) for k in design):
        return None, "design field names must be text"
    if "fluids" not in design:
        return None, "missing field: fluids"
    rest = dict(design)
    fluids, err = _check_fluids(rest.pop("fluids"), compounds)
    if err:
        return None, err
    choices = {}
    for key in CHOICE_SPEC:
        if key not in rest:
            return None, f"missing field: {key}"
        value = rest.pop(key)
        if not isinstance(value, str) or value not in QUENCHES:
            return None, f"{key} must be one of: {', '.join(QUENCHES)}"
        choices[key] = value
    clean, err = validate(DESIGN_SPEC, rest)
    if err:
        return None, err
    if any(0.0 < t < HANDLING_MIN for t in clean["timepoints_min"]):
        return None, f"a time point other than 0 cannot be earlier than {HANDLING_MIN:g} min"
    for fluid in FLUIDS:
        arms = sum(fluid in chosen for chosen in fluids.values())
        arms += 1 if arms and clean["activity_control"] else 0
        if arms * clean["replicates"] * len(clean["timepoints_min"]) > 96:
            return None, f"more than 96 {fluid} samples; that quench plate is full"
    return {**clean, **choices, "fluids": fluids}, None


def simulate_gi(profile: RigProfile, rng: np.random.Generator, peptides: list[Peptide], lot_activity: dict[str, float],
                design: dict, is_sd: float = IS_SD) -> tuple[list[dict], Rig]:
    """Run one validated design. lot_activity maps fluid to the activity of a fresh preparation (1 = specification)."""
    rig = Rig(profile, rng)
    by_id = {p.cid: p for p in peptides}
    arms = {f: [p.cid for p in peptides if f in design["fluids"].get(p.cid, ())] for f in FLUIDS}
    fluids = [f for f in FLUIDS if arms[f]]
    if design["activity_control"]:
        by_id[CONTROL_ID] = CONTROL
        for f in fluids:
            arms[f].append(CONTROL_ID)
    potency = {f: lot_activity[f] * (1.0 if design["fresh_enzyme"] else stored_activity(f, STALE_AGE_H)) for f in fluids}
    made = ("made fresh for this run" if design["fresh_enzyme"]
            else f"made {STALE_AGE_H:g} h ago and kept at {STORAGE_TEMP_C:g} C")

    stock = Plate("stock")
    stock_well: dict[str, tuple] = {}
    for cid in by_id:
        if any(cid in arms[f] for f in fluids):
            stock_well[cid] = (len(stock_well), 0)
            rig.fill(stock, stock_well[cid], 300.0, {cid: 1000.0 * (FLUID_UL + SPIKE_UL) / SPIKE_UL},
                     what="control substrate stock" if cid == CONTROL_ID else f"{cid} stock")
    inc = {f: Plate(f"{f} incubation") for f in fluids}
    quench = {f: Plate(f"{f} quench") for f in fluids}
    wells: dict[str, list[tuple[str, int, tuple]]] = {f: [] for f in fluids}
    for f in fluids:
        free = iter(usable_wells(True))
        for cid in arms[f]:
            for rep in range(design["replicates"]):
                w = next(free)
                rig.fill(inc[f], w, FLUID_UL, {ENZYME[f]: 1.0}, what=f"{FLUID_NAME[f]}, {made}")
                rig.transfer(stock, stock_well[cid], inc[f], w, SPIKE_UL, new_tip=True, mix_after=True)
                wells[f].append((cid, rep + 1, w))

    def digest(f: str, t0: float, t1: float) -> None:
        """Loss in the incubation wells between clock times t0 and t1: enzyme at its decaying activity, plus acid."""
        plate = inc[f]
        enzyme = np.where(plate.vol > 0, plate.a(ENZYME[f]) / np.maximum(plate.vol, 1e-9), 0.0)
        exposure = potency[f] * enzyme_exposure(t0, t1, ACTIVITY_HL_37C_MIN[f])
        for cid in arms[f]:
            p = by_id[cid]
            plate.a(cid)[:] *= np.exp(-(enzymatic_rate(p, f) * exposure * enzyme + acid_rate(p, FLUID_PH[f]) * (t1 - t0)))

    quench_free = {f: iter(usable_wells(False)) for f in fluids}
    taken: dict[str, list[tuple]] = {f: [] for f in fluids}  # (id, replicate, labelled time, quench well, clock at quench)
    clock = 0.0
    for t in design["timepoints_min"]:
        due = HANDLING_MIN if t == 0.0 and not design["prequenched_t0"] else t
        if due > clock:
            for f in fluids:
                rig.incubate(inc[f], due - clock, True)
                digest(f, clock, due)
            clock = due
        for f in fluids:
            method = design[f"{f}_quench"]
            label = f"{method} quench" + (" with internal standard" if design["internal_standard"] else "")
            if t == 0.0 and design["prequenched_t0"]:
                label += ", time zero quenched before the enzyme acts"
            for cid, rep, w in wells[f]:
                q = next(quench_free[f])
                rig.fill(quench[f], q, QUENCH_UL, {"IS": 1.0} if design["internal_standard"] else None, what=label)
                rig.transfer(inc[f], w, quench[f], q, ALIQUOT_UL, new_tip=True, mix_after=True)
                taken[f].append((cid, rep, t, q, clock))

    # the quenched samples wait: on the plate until the last time point, through the hold, then in the queue
    cold = q10_factor(design["hold_temp_c"], 37.0, Q10)
    nominal_t0 = 1000.0 * ALIQUOT_UL  # 1 uM in the fluid, as amount in one aliquot
    samples: list[dict] = []
    position = 0
    for f in fluids:
        method = design[f"{f}_quench"]
        rig.log.append({"op": "hold", "dst_plate": quench[f].name, "minutes": design["hold_min"],
                        "temp_c": design["hold_temp_c"],
                        "what": f"quenched samples wait for injection at {design['hold_temp_c']:g} C"})
        for cid, rep, t, q, quenched_at in taken[f]:
            wait = (clock - quenched_at) + design["hold_min"] + position * INJECTION_MIN
            position += 1
            p = by_id[cid]
            live = (potency[f] * remaining(ACTIVITY_HL_37C_MIN[f], quenched_at) * RESIDUAL_ACTIVITY[f][method]
                    * quench[f].conc(ENZYME[f], q))
            loss = (enzymatic_rate(p, f) * live * cold * enzyme_exposure(0.0, wait, ACTIVITY_HL_37C_MIN[f] / cold)
                    + acid_rate(p, QUENCHED_PH[f][method]) * cold * wait)
            amount = quench[f].a(cid)[q] * float(np.exp(-loss))
            injection = float(np.exp(rng.normal(0.0, INJECTION_SD)))
            area = amount * injection + rng.normal(0.0, NOISE_FLOOR * nominal_t0)
            sample = {"id": cid, "fluid": f, "replicate": rep, "time_min": t, "area": round(float(area), 1)}
            if design["internal_standard"]:
                is_area = quench[f].a("IS")[q] * 100.0 * injection * max(1 + rng.normal(0.0, is_sd), 0.05)
                sample["is_area"] = round(float(is_area), 1)
            samples.append(sample)
    return samples, rig


def run_cost(n_samples: int, tips: int, fresh_preparations: int) -> float:
    return 10.0 + 2.0 * n_samples + 0.1 * tips + ENZYME_PREP_COST * fresh_preparations


def run_qc(samples: list[dict] | None) -> dict:
    """The run's own QC, from its data alone: time-zero replicate CV per arm and the activity control per fluid."""
    cvs = {f: t0_cv([s for s in samples or [] if s["fluid"] == f]) for f in FLUIDS}
    cvs = {f: v for f, v in cvs.items() if v}
    activity = enzyme_activity(samples or [], CONTROL_ID, CONTROL_HALF_LIFE_MIN, ACTIVITY_HL_37C_MIN,
                               use_lloq=True, lloq_frac=LLOQ_FRAC)
    cv_ok = all(np.isfinite(v) and v <= T0_CV_GATE for per in cvs.values() for v in per.values())
    activity_ok = all(f in activity and ACTIVITY_GATE[0] <= activity[f] <= ACTIVITY_GATE[1] for f in cvs)
    return {"qc_pass": bool(cvs) and cv_ok and activity_ok,
            "t0_cv": {f: {k: (round(v, 4) if np.isfinite(v) else None) for k, v in per.items()} for f, per in cvs.items()},
            "enzyme_activity": {f: round(v, 4) for f, v in activity.items()}}


def score_half_life(estimate, truth: float) -> float:
    if isinstance(estimate, bool) or not isinstance(estimate, (int, float, np.integer, np.floating)):
        return 0.0
    try:
        estimate = float(estimate)
    except (OverflowError, ValueError):
        return 0.0
    if not np.isfinite(estimate) or estimate <= 0:
        return 0.0
    return log_error_score(float(np.log10(np.clip(estimate, HL_LO, HL_HI))), float(np.log10(np.clip(truth, HL_LO, HL_HI))),
                           TOLERANCE)


class GiFluidStability(Env):
    name = "gi_fluid_stability"
    cost_reference = COST_REFERENCE
    salted = True  # hidden draws depend on the environment's name: no other environment's truth gives this one away

    def __init__(self, seed: int = 0, rig="nominal", n_compounds: int = 3):
        super().__init__(seed, rig)
        if not 1 <= n_compounds <= ROWS - 1:
            raise ValueError(f"n_compounds must be between 1 and {ROWS - 1}")
        self.peptides = [draw_peptide(self._truth_rng, f"cpd{i + 1}") for i in range(n_compounds)]
        self.lot_activity = {f: float(np.exp(self._truth_rng.normal(0.0, LOT_ACTIVITY_SD))) for f in FLUIDS}
        self.samples: list[dict] | None = None

    def spec(self) -> dict:
        return {
            "env": self.name,
            "task": ("Measure the half-life of each peptide in simulated gastric fluid (SGF: pepsin, pH 1.2) and "
                     "simulated intestinal fluid (SIF: pancreatin, pH 6.8) at 37 C, read by LC-MS. Run one time "
                     "course, then report half-lives in minutes. The half-life asked for is the one at "
                     "specification enzyme activity. Enzymes lose activity in storage and during the run, and a "
                     "quenched sample keeps digesting while it waits unless the quench stops that enzyme."),
            "actions": {
                "run": {"design": {**spec_help(DESIGN_SPEC),
                                   **{k: {"type": "choice", "options": dict(QUENCH_HELP), "help": text}
                                      for k, text in CHOICE_SPEC.items()},
                                   "fluids": {"type": "object", "help": FLUIDS_HELP}}},
                "report": {"half_life_min": "object mapping compound id to an object of fluid -> half-life in "
                                            "minutes; a pair left out scores 0, and an answer for a pair that "
                                            "was not run counts at half"},
            },
            "fluids": {f: {"enzyme": ENZYME[f], "ph": FLUID_PH[f],
                           "enzyme_activity_half_life_37C_min": ACTIVITY_HL_37C_MIN[f]} for f in FLUIDS},
            "activity_control": {"id": CONTROL_ID, "reference_half_life_min": CONTROL_HALF_LIFE_MIN,
                                 "acceptance_window": list(ACTIVITY_GATE)},
            "labware": {"fluid_ul": FLUID_UL, "aliquot_ul": ALIQUOT_UL, "quench_ul": QUENCH_UL,
                        "horizon_min": HORIZON_MIN, "quantifiable_range_min": [HL_LO, HL_HI],
                        "injection_cycle_min": INJECTION_MIN, "injection_order": "plate order, the SGF plate first",
                        "nominal_dead_volume_ul": RigProfile().dead_volume_ul},
            "reward": REWARD_RULE,
        }

    def reset(self) -> dict:
        return {"phase": "start", "compounds": [p.cid for p in self.peptides], "fluids": list(FLUIDS),
                "spec": self.spec()}

    def truth(self) -> dict:
        return {"half_life_min": {p.cid: {f: round(true_half_life(p, f), 2) for f in FLUIDS} for p in self.peptides},
                "acid_half_life_min": {p.cid: None if p.acid_min is None else round(p.acid_min, 1) for p in self.peptides},
                "lot_activity": {f: round(v, 3) for f, v in self.lot_activity.items()}}

    def step(self, action: dict):
        err = self._check_open(action)
        if err:
            return self._reject(err)
        if action["type"] == "run":
            if self.samples is not None:
                return self._reject("the time course has been run; report")
            design, err = check_design(action.get("design"), [p.cid for p in self.peptides])
            if err:
                return self._reject(err)
            self.samples, rig = simulate_gi(self.profile, self._noise, self.peptides, self.lot_activity, design)
            self._record(rig.log)
            fresh = len({s["fluid"] for s in self.samples}) if design["fresh_enzyme"] else 0
            self.cost += run_cost(len(self.samples), rig.tips, fresh)
            return {"phase": "data", "samples": copy.deepcopy(self.samples), "cost": round(self.cost, 2)}, 0.0, False, {}
        if action["type"] == "report":
            answers = action.get("half_life_min")
            answers = answers if isinstance(answers, dict) else {}
            qc = run_qc(self.samples)
            ran = {(s["id"], s["fluid"]) for s in self.samples or []}
            scores, credited = [], []
            for p in self.peptides:
                row = answers.get(p.cid)
                row = row if isinstance(row, dict) else {}
                for f in FLUIDS:
                    scores.append(score_half_life(row.get(f), true_half_life(p, f)))
                    # an answer for a pair that was not measured is a guess: it never counts in full
                    credited.append(scores[-1] * (1.0 if qc["qc_pass"] and (p.cid, f) in ran else 0.5))
            accuracy = float(np.mean(scores))
            reward = float(np.mean(credited)) - self._cost_penalty()
            self.done = True
            info = {"accuracy": accuracy, **qc, "pairs_run": sorted(f"{cid}/{f}" for cid, f in ran if cid != CONTROL_ID),
                    "cost": self.cost, "truth": self.truth()}
            return {"phase": "done"}, reward, True, info
        return self._reject("unknown action type; use run or report")  # the type is not echoed: it may be anything
