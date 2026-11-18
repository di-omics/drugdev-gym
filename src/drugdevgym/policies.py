"""Reference policies: one expert per environment, named bench mistakes ("traps"), and a null.

A trap is the expert protocol with one named mistake a real bench makes. Environment QC
(envqc.py) then measures whether the reward can tell the expert from each trap. A policy is
any object with act(observation) -> action.
"""
from __future__ import annotations

import itertools

import numpy as np

from .analysis import (analyte_zprime, control_zprime, enzyme_activity, estimate_fluid_half_life, estimate_half_life,
                       estimate_pec50, fit_4pl, invert_4pl, zprime)
from .envs import gi_fluid_stability as gi
from .envs.lead_opt import OPTIONS, POSITIONS, analog_id, desirability, parse_id
from .envs.review import PLATE_FAULTS, TIMECOURSE_FAULTS
from .judges import FixedJudge, ForcedVerdict, PlateJudge, RandomJudge, TimecourseJudge

# ---------------------------------------------------------------- dose response

DR_EXPERT = {"top_conc_nM": 30000.0, "dilution_factor": 4.0, "n_points": 10, "replicates": 1,
             "transfer_ul": 15.0, "mix": True, "new_tips": True, "avoid_edges": True, "lid": True,
             "standards": True, "controls_per_side": 6}

DR_CHOICES = {
    "top_conc_nM": [100.0, 1000.0, 10000.0, 30000.0, 100000.0],
    "dilution_factor": [2.0, 3.16, 4.0, 10.0],
    "n_points": [6, 8, 10, 12],
    "replicates": [1, 2],
    "transfer_ul": [5.0, 10.0, 15.0, 25.0, 50.0],
    "mix": [False, True],
    "new_tips": [False, True],
    "avoid_edges": [False, True],
    "lid": [False, True],
    "standards": [False, True],
    "controls_per_side": [0, 2, 4, 6, 8],
}


class DosePolicy:
    """Run one design, repeat once if the controls fail Z' 0.5, then fit and report."""

    def __init__(self, design: dict | None, use_standards: bool = True):
        self.design = design
        self.use_standards = use_standards

    def act(self, obs: dict) -> dict:
        if obs.get("phase") == "start":
            if self.design is None:
                return {"type": "report", "pEC50": {c: 7.5 for c in obs["compounds"]}}
            return {"type": "run_plate", "design": self.design}
        if obs.get("phase") == "data":
            z = control_zprime(obs["wells"])
            if obs["plates_left"] > 0 and np.isfinite(z) and z < 0.5:
                return {"type": "run_plate", "design": self.design}
            return {"type": "report", "pEC50": estimate_pec50(obs["wells"], self.use_standards),
                    "plates_used": [obs["plate"]]}
        return {"type": "report", "pEC50": {}}


class GateBypass:
    """A known exploit of the dose_response reward, kept on the QC card so it stays visible.

    Answer from a plate with no controls, then run a cheap token plate that has controls and
    declare it as the plate the answer rests on. The QC gate only sees the declared plate.
    """

    def __init__(self):
        self.answer: dict = {}

    def act(self, obs: dict) -> dict:
        if obs.get("phase") == "start":
            return {"type": "run_plate", "design": {**DR_EXPERT, "controls_per_side": 0}}
        if obs.get("phase") == "data" and obs["plate"] == 1:
            self.answer = estimate_pec50(obs["wells"], True)
            token = {**DR_EXPERT, "n_points": 4, "standards": False, "new_tips": False, "controls_per_side": 8}
            return {"type": "run_plate", "design": token}
        return {"type": "report", "pEC50": self.answer, "plates_used": [2]}


def random_dr_design(rng: np.random.Generator) -> dict:
    return {k: v[int(rng.integers(len(v)))] for k, v in DR_CHOICES.items()}


def _dr(**changes):
    return lambda seed: DosePolicy({**DR_EXPERT, **changes})


DOSE_RESPONSE = {
    "expert": _dr(),
    "no_mix": _dr(mix=False),
    "reuse_tips": _dr(new_tips=False),
    "edge_wells": _dr(avoid_edges=False),
    "no_lid": _dr(lid=False),
    "raw_ratio": lambda seed: DosePolicy({**DR_EXPERT, "standards": False}, use_standards=False),
    "narrow_range": _dr(top_conc_nM=1000.0, dilution_factor=2.0, n_points=8, transfer_ul=40.0),
    "thin_ladder": _dr(transfer_ul=5.0),
    "few_controls": _dr(controls_per_side=2),
    "random": lambda seed: DosePolicy(random_dr_design(np.random.default_rng([seed, 99]))),
    "null": lambda seed: DosePolicy(None),
    "gate_bypass": lambda seed: GateBypass(),
}

# ---------------------------------------------------------------- stability

ST_EXPERT = {"timepoints_min": [0, 5, 15, 30, 60, 120, 240, 360], "replicates": 2, "internal_standard": True,
             "mix_after_spike": True, "lid": True, "avoid_edges": True}


class StabilityPolicy:
    def __init__(self, design: dict | None, use_lloq: bool = True):
        self.design = design
        self.use_lloq = use_lloq

    def act(self, obs: dict) -> dict:
        if obs.get("phase") == "start":
            if self.design is None:
                return {"type": "report", "half_life_min": {c: 60.0 for c in obs["compounds"]}}
            return {"type": "run", "design": self.design}
        if obs.get("phase") == "data":
            return {"type": "report", "half_life_min": estimate_half_life(obs["samples"], self.use_lloq)}
        return {"type": "report", "half_life_min": {}}


def _st(**changes):
    return lambda seed: StabilityPolicy({**ST_EXPERT, **changes})


def _random_st(seed: int) -> StabilityPolicy:
    rng = np.random.default_rng([seed, 99])
    n = int(rng.integers(2, 9))
    times = sorted({float(t) for t in rng.choice([0, 5, 10, 15, 30, 45, 60, 90, 120, 180, 240, 300, 360], n, replace=False)})
    return StabilityPolicy({"timepoints_min": times, "replicates": int(rng.integers(1, 4)),
                            "internal_standard": bool(rng.integers(2)), "mix_after_spike": bool(rng.integers(2)),
                            "lid": bool(rng.integers(2)), "avoid_edges": bool(rng.integers(2))},
                           use_lloq=bool(rng.integers(2)))


STABILITY = {
    "expert": _st(),
    "two_point": _st(timepoints_min=[0, 360]),
    "early_only": _st(timepoints_min=[0, 5, 10, 15, 20, 30]),
    "no_internal_standard": _st(internal_standard=False),
    "no_mix": _st(mix_after_spike=False),
    "edge_wells": _st(avoid_edges=False),
    "no_lid": _st(lid=False),
    "no_lid_on_edge_wells": _st(lid=False, avoid_edges=False),  # the two together: an interaction
    "single_replicate": _st(replicates=1),
    "no_lloq": lambda seed: StabilityPolicy(dict(ST_EXPERT), use_lloq=False),
    "random": _random_st,
    "null": lambda seed: StabilityPolicy(None),
}

# ---------------------------------------------------------------- assay window


class WindowSearch:
    """Coarse cell-density pilot, a refining pilot around a chemistry-guided centre, then commit the
    best condition of the refining pilot.

    metric="analyte" scores each condition by Z' after back-calculation through the standard
    curve; metric="raw" scores the raw ratio, which is the wrong window for this assay format.
    """

    def __init__(self, commit_lid: bool = True, commit_avoid_edges: bool = True, metric: str = "analyte"):
        self.commit = {"lid": commit_lid, "avoid_edges": commit_avoid_edges}
        self.metric = metric

    def _score(self, obs: dict, c: dict) -> float:
        z = analyte_zprime(obs["standards"], c["pos"], c["neg"]) if self.metric == "analyte" else zprime(c["pos"], c["neg"])
        return z if np.isfinite(z) else -np.inf

    @staticmethod
    def _target_cells(obs: dict) -> float | None:
        """Chemistry-guided: the density that centres the two controls on the standard curve midpoint."""
        curve = fit_4pl(np.log10([s["std_nM"] for s in obs["standards"]]), [s["ratio"] for s in obs["standards"]])
        if curve is None:
            return None
        lo, hi = sorted((curve["bottom"], curve["top"]))
        rates: dict[str, list[float]] = {"pos": [], "neg": []}
        for c in obs["conditions"]:
            for side in rates:
                r = float(np.median(c[side]))
                if lo + 0.1 * (hi - lo) <= r <= hi - 0.1 * (hi - lo):  # only where the curve can be read
                    rates[side].append(float(invert_4pl(curve, r)) / c["cells_k"])
        if not rates["pos"] or not rates["neg"]:
            return None
        return float(10.0 ** curve["logec50"] / np.sqrt(np.median(rates["pos"]) * np.median(rates["neg"])))

    def act(self, obs: dict) -> dict:
        if obs.get("phase") == "start":
            return {"type": "run_pilot", "lid": True, "avoid_edges": False,
                    "conditions": [{"cells_k": c, "stim_min": 60.0} for c in (2.5, 7.0, 20.0, 50.0, 80.0)]}
        if obs.get("phase") == "data":
            best = max((self._score(obs, c), c["cells_k"], c["stim_min"]) for c in obs["conditions"])
            if obs["pilot"] == 1:
                b = self._target_cells(obs)
                b = best[1] if b is None else b
                grid = [(b, 60.0), (b * 0.6, 60.0), (b * 1.6, 60.0), (b, 30.0), (b, 120.0)]
                return {"type": "run_pilot", "lid": True, "avoid_edges": False,
                        "conditions": [{"cells_k": float(np.clip(c, 1.0, 80.0)), "stim_min": s} for c, s in grid]}
            return {"type": "commit", "cells_k": best[1], "stim_min": best[2], **self.commit}
        return {"type": "commit", "cells_k": 20.0, "stim_min": 60.0, "lid": True, "avoid_edges": True}


class WindowFixed:
    def __init__(self, cells_k: float, stim_min: float, lid: bool = True, avoid_edges: bool = True):
        self.action = {"type": "commit", "cells_k": cells_k, "stim_min": stim_min, "lid": lid, "avoid_edges": avoid_edges}

    def act(self, obs: dict) -> dict:
        return self.action


def _random_window(seed: int) -> WindowFixed:
    rng = np.random.default_rng([seed, 99])
    return WindowFixed(float(np.exp(rng.uniform(0.0, np.log(80.0)))), float(rng.uniform(15.0, 180.0)),
                       bool(rng.integers(2)), bool(rng.integers(2)))


ASSAY_WINDOW = {
    "expert": lambda seed: WindowSearch(),
    "edge_wells": lambda seed: WindowSearch(commit_avoid_edges=False),
    "no_lid": lambda seed: WindowSearch(commit_lid=False),
    "raw_ratio_zprime": lambda seed: WindowSearch(metric="raw"),
    "max_cells": lambda seed: WindowFixed(80.0, 60.0),
    "min_cells": lambda seed: WindowFixed(2.5, 60.0),
    "kit_default": lambda seed: WindowFixed(20.0, 60.0),
    "random": _random_window,
}

# ---------------------------------------------------------------- lead optimization

ALL_ANALOGS = list(itertools.product(range(OPTIONS), repeat=POSITIONS))
LATIN = [tuple((m * j) % OPTIONS for m in range(1, POSITIONS + 1)) for j in range(OPTIONS)]


def _one_hot(analogs) -> np.ndarray:
    x = np.zeros((len(analogs), POSITIONS * OPTIONS))
    for row, a in enumerate(analogs):
        for i, opt in enumerate(a):
            x[row, i * OPTIONS + opt] = 1.0
    return x


def additive_model(analogs, y, ridge: float = 1.0) -> np.ndarray:
    """Ridge fit of an additive (Free-Wilson style) model; returns predictions for the whole library."""
    x, y = _one_hot(analogs), np.asarray(y, float)
    mean = y.mean()
    w = np.linalg.solve(x.T @ x + ridge * np.eye(x.shape[1]), x.T @ (y - mean))
    return _one_hot(ALL_ANALOGS) @ w + mean


class CampaignPolicy:
    """Latin first batch, then half model-guided and half exploratory batches; nominate the best measured."""

    def __init__(self, seed: int, dr_design: dict, use_standards: bool, st_design: dict, use_lloq: bool,
                 model: bool = True, rounds: int = 4):
        self.rng = np.random.default_rng([seed, 99])
        self.dr, self.use_standards, self.st, self.use_lloq = dr_design, use_standards, st_design, use_lloq
        self.model, self.rounds = model, rounds
        self.measured: dict[tuple, tuple[float, float]] = {}
        self.pending: list[tuple] = []

    def _untested(self, n: int, exclude=()) -> list[tuple]:
        pool = [a for a in ALL_ANALOGS if a not in self.measured and a not in exclude]
        return [pool[i] for i in self.rng.choice(len(pool), n, replace=False)]

    def _batch(self) -> list[tuple]:
        if not self.measured:
            return (LATIN + self._untested(1, LATIN)) if self.model else self._untested(6)
        if not self.model:
            return self._untested(6)
        analogs = list(self.measured)
        pred = desirability(additive_model(analogs, [self.measured[a][0] for a in analogs]),
                            additive_model(analogs, [self.measured[a][1] for a in analogs]))
        ranked = [ALL_ANALOGS[i] for i in np.argsort(-pred) if ALL_ANALOGS[i] not in self.measured]
        return ranked[:3] + self._untested(3, ranked[:3])

    def _make(self) -> dict:
        self.pending = self._batch()
        return {"type": "make_test", "analogs": [list(a) for a in self.pending],
                "dose_response": self.dr, "stability": self.st}

    def act(self, obs: dict) -> dict:
        if obs.get("phase") == "start":
            if self.rounds == 0:
                return {"type": "nominate", "analog": list(self._untested(1)[0])}
            return self._make()
        if obs.get("phase") == "data":
            potency = estimate_pec50(obs["dose_response"]["wells"], self.use_standards)
            half_life = estimate_half_life(obs["stability"]["samples"], self.use_lloq)
            for a in self.pending:
                cid = analog_id(a)
                if cid in potency and cid in half_life:
                    self.measured[parse_id(cid)] = (potency[cid], float(np.log10(half_life[cid])))
            if obs["rounds_left"] > 0 and obs["round"] < self.rounds:
                return self._make()
        if not self.measured:
            return {"type": "nominate", "analog": list(ALL_ANALOGS[0])}
        best = max(self.measured, key=lambda a: float(desirability(*self.measured[a])))
        return {"type": "nominate", "analog": list(best)}


LO_DR = {**DR_EXPERT, "avoid_edges": False}  # six analogs need the outer ring
LEAD_OPT = {
    "expert": lambda seed: CampaignPolicy(seed, LO_DR, True, ST_EXPERT, True),
    # one upstream measurement mistake each, with the campaign strategy unchanged
    "upstream_no_mix": lambda seed: CampaignPolicy(seed, {**LO_DR, "mix": False}, True, ST_EXPERT, True),
    "upstream_raw_ratio": lambda seed: CampaignPolicy(seed, {**LO_DR, "standards": False}, False, ST_EXPERT, True),
    "upstream_two_point": lambda seed: CampaignPolicy(seed, LO_DR, True, {**ST_EXPERT, "timepoints_min": [0, 360]}, True),
    "upstream_no_internal_standard": lambda seed: CampaignPolicy(seed, LO_DR, True,
                                                                 {**ST_EXPERT, "internal_standard": False}, True),
    "no_model": lambda seed: CampaignPolicy(seed, LO_DR, True, ST_EXPERT, True, model=False),
    "one_round": lambda seed: CampaignPolicy(seed, LO_DR, True, ST_EXPERT, True, rounds=1),
    "null": lambda seed: CampaignPolicy(seed, LO_DR, True, ST_EXPERT, True, rounds=0),
}

# ---------------------------------------------------------------- review (judgment)

PLATE_REVIEW = {
    "expert": lambda seed: PlateJudge(),  # scripted QC rules; a human expert's calls are scored separately
    "always_accept": lambda seed: ForcedVerdict(PlateJudge(), "accept"),
    "always_repeat": lambda seed: ForcedVerdict(PlateJudge(), "repeat"),
    "random": lambda seed: RandomJudge(seed, PLATE_FAULTS),
    "null": lambda seed: FixedJudge("accept", "none"),
    "constant_repeat": lambda seed: FixedJudge("repeat", "none"),
}
TIMECOURSE_REVIEW = {
    "expert": lambda seed: TimecourseJudge(),
    "always_accept": lambda seed: ForcedVerdict(TimecourseJudge(), "accept"),
    "always_repeat": lambda seed: ForcedVerdict(TimecourseJudge(), "repeat"),
    "random": lambda seed: RandomJudge(seed, TIMECOURSE_FAULTS),
    "null": lambda seed: FixedJudge("accept", "none"),
    "constant_repeat": lambda seed: FixedJudge("repeat", "none"),
}

# ---------------------------------------------------------------- GI fluid stability

GI_EXPERT = {"timepoints_min": [0, 2, 5, 10, 20, 45, 120, 240], "replicates": 2, "internal_standard": True,
             "fresh_enzyme": True, "activity_control": True, "prequenched_t0": True, "sgf_quench": "base",
             "sif_quench": "organic", "hold_min": 30.0, "hold_temp_c": 4.0}
GI_NULL_GUESS = {"sgf": 60.0, "sif": 10.0}
GI_WARM_OVERNIGHT = {"hold_min": 960.0, "hold_temp_c": 25.0}  # the quench plate left on the bench until the morning
GI_TIME_CHOICES = [0, 1, 2, 5, 10, 15, 20, 30, 45, 60, 90, 120, 180, 240]


class GiFluidPolicy:
    """Run one design in the chosen fluids, read the enzyme activity off the control, fit and report.

    correct_decay fits against minutes of digestion at the starting activity instead of the clock,
    for the compounds and for the control substrate alike; correct_activity scales the result to
    specification activity with the control substrate.
    """

    def __init__(self, design: dict | None, fluids=gi.FLUIDS, correct_decay: bool = True, correct_activity: bool = True,
                 use_lloq: bool = True, pick: np.random.Generator | None = None):
        self.design = design
        self.fluids = tuple(fluids)
        self.flags = {"correct_decay": correct_decay, "correct_activity": correct_activity, "use_lloq": use_lloq}
        self.pick = pick  # when set, each compound gets its own random set of fluids

    def _fluids(self, compounds: list[str]) -> dict:
        if self.pick is None:
            return {c: list(self.fluids) for c in compounds}
        options = [["sgf"], ["sif"], ["sgf", "sif"]]
        return {c: options[int(self.pick.integers(len(options)))] for c in compounds}

    def act(self, obs: dict) -> dict:
        if obs.get("phase") == "start":
            if self.design is None:
                return {"type": "report", "half_life_min": {c: dict(GI_NULL_GUESS) for c in obs["compounds"]}}
            return {"type": "run", "design": {**self.design, "fluids": self._fluids(obs["compounds"])}}
        if obs.get("phase") == "data":
            samples, lloq = obs["samples"], self.flags["use_lloq"]
            decay = gi.ACTIVITY_HL_37C_MIN if self.flags["correct_decay"] else None
            activity = None
            if self.flags["correct_activity"]:
                activity = enzyme_activity(samples, gi.CONTROL_ID, gi.CONTROL_HALF_LIFE_MIN, decay, lloq)
            compounds = [s for s in samples if s["id"] != gi.CONTROL_ID]
            return {"type": "report", "half_life_min": estimate_fluid_half_life(compounds, decay, activity, lloq,
                                                                              lo=gi.HL_LO, hi=gi.HL_HI)}
        return {"type": "report", "half_life_min": {}}


def _gi(**changes):
    return lambda seed: GiFluidPolicy({**GI_EXPERT, **changes})


def _random_gi(seed: int) -> GiFluidPolicy:
    rng = np.random.default_rng([seed, 99])
    n = int(rng.integers(2, 9))
    design = {"timepoints_min": sorted(float(t) for t in rng.choice(GI_TIME_CHOICES, n, replace=False)),
              "replicates": int(rng.integers(1, 4)), "internal_standard": bool(rng.integers(2)),
              "fresh_enzyme": bool(rng.integers(2)), "activity_control": bool(rng.integers(2)),
              "prequenched_t0": bool(rng.integers(2)), "sgf_quench": str(rng.choice(gi.QUENCHES)),
              "sif_quench": str(rng.choice(gi.QUENCHES)), "hold_min": float(rng.choice([0.0, 30.0, 120.0, 480.0, 960.0])),
              "hold_temp_c": float(rng.choice([4.0, 10.0, 25.0]))}
    return GiFluidPolicy(design, correct_decay=bool(rng.integers(2)), correct_activity=bool(rng.integers(2)),
                         use_lloq=bool(rng.integers(2)), pick=rng)


GI_FLUID_STABILITY = {
    "expert": _gi(),
    "acid_quench_in_sgf": _gi(sgf_quench="acid"),
    # the wrong quench plus the wait that lets it act: an interaction, not a single mistake (see GI_INTERACTIONS)
    "acid_quench_in_sgf_warm_overnight": _gi(sgf_quench="acid", **GI_WARM_OVERNIGHT),
    "stale_enzyme": _gi(fresh_enzyme=False),
    "no_activity_control": _gi(activity_control=False),
    "late_first_point": _gi(timepoints_min=[0, 15, 30, 45, 60, 90, 120, 240]),
    "no_true_time_zero": _gi(prequenched_t0=False),
    "no_internal_standard": _gi(internal_standard=False),
    "ignore_enzyme_decay": lambda seed: GiFluidPolicy(dict(GI_EXPERT), correct_decay=False),  # an analysis choice
    "random": _random_gi,
    "null": lambda seed: GiFluidPolicy(None),
}
# Labelled interactions: the policy on the left is the single mistake on the right plus GI_WARM_OVERNIGHT.
# Every other GI trap is the expert with exactly one change.
GI_INTERACTIONS = {"acid_quench_in_sgf_warm_overnight": "acid_quench_in_sgf"}

POLICIES = {"dose_response": DOSE_RESPONSE, "stability": STABILITY, "assay_window": ASSAY_WINDOW,
            "lead_opt": LEAD_OPT, "plate_review": PLATE_REVIEW, "timecourse_review": TIMECOURSE_REVIEW,
            "gi_fluid_stability": GI_FLUID_STABILITY}

# ---------------------------------------------------------------- Gi cAMP potency

from .analysis import ec_fraction, estimate_pec50_anchored, titration_curve, titration_level  # noqa: E402
from .envs.gi_camp_potency import DESIGN_SPEC as GI_SPEC  # noqa: E402
from .envs.gi_camp_potency import FSK_MAX_UM, LADDER_MAX_UL, STD_SERIES_NM  # noqa: E402

# The expert plan. forskolin_uM and cells_k are either a number or a rule applied to the titration:
# "EC80" = the fitted forskolin EC80; "titrate" = the cell number that puts the forskolin control
# at GI_TARGET_FOLD times the standard-curve midpoint.
CAMP_EXPERT = {"forskolin_uM": "EC80", "cells_k": "titrate", "ibmx": True, "stim_min": 30.0,
             "top_conc_nM": 30000.0, "dilution_factor": 4.0, "n_points": 10, "replicates": 1,
             "transfer_ul": 15.0, "mix": True, "standards": True, "controls_per_side": 6}
GI_PILOT_CELLS_K = 10.0  # cell number of the forskolin titration
GI_PILOT_REPLICATES = 2
GI_TARGET_FOLD = 3.0  # forskolin control at 3x the kit midpoint, so the inhibited wells fall through the middle
GI_DEFAULT_FORSKOLIN_UM = 10.0  # used when the titration gives no curve to read an EC80 from
GI_HANDBOOK_CELLS_K = 30.0  # cell number of the skip_titration trap: a handbook figure, not read from these cells
GI_CELL_RANGE = GI_SPEC["cells_k"][1:3]  # the cell numbers a design may ask for

GI_CHOICES = {
    "forskolin_uM": ["EC80", 0.0, 1.0, 10.0, 100.0],
    "cells_k": ["titrate", 2.0, 10.0, 40.0, 80.0],
    "ibmx": [False, True],
    "stim_min": [10.0, 30.0, 60.0, 120.0],
    "top_conc_nM": [100.0, 1000.0, 10000.0, 30000.0, 100000.0],
    "dilution_factor": [2.0, 3.16, 4.0, 10.0],
    "n_points": [6, 8, 10, 12],
    "replicates": [1, 2],
    "transfer_ul": [5.0, 10.0, 15.0, 25.0, 50.0],
    "mix": [False, True],
    "standards": [False, True],
    "controls_per_side": [0, 2, 4, 6, 8],
}


def gi_conditions(plan: dict, wells: list[dict], pilot_cells_k: float = GI_PILOT_CELLS_K) -> tuple[float, float]:
    """Forskolin (uM) and cells (thousand per well) for the potency plate, from a plan and a titration plate."""
    curve = titration_curve(wells)
    forskolin = plan["forskolin_uM"]
    if forskolin == "EC80":
        forskolin = GI_DEFAULT_FORSKOLIN_UM
        if curve is not None and curve["fitted"]:
            forskolin = float(np.clip(ec_fraction(curve["ec50_uM"], curve["hill"], 0.8), 0.0, FSK_MAX_UM))
    cells = plan["cells_k"]
    if cells == "titrate":
        cells = pilot_cells_k
        level = titration_level(curve, forskolin) if curve is not None else None
        if level is not None and level > 0:  # cAMP scales with the cell number
            cells = float(np.clip(pilot_cells_k * GI_TARGET_FOLD * curve["kit_mid_nM"] / level, *GI_CELL_RANGE))
    return float(forskolin), float(cells)


class GiCampPolicy:
    """Titrate forskolin, set forskolin and cells from it, run one potency plate, fit and report.

    plan is CAMP_EXPERT with any changes. Two analysis choices: use_standards (cAMP through the
    standard curve, or the raw ratio) and anchor (pin the uninhibited plateau of every curve on
    the forskolin-only controls, or leave it free). A plan with fixed forskolin and fixed cells
    skips the titration. The anchor was added to the expert after its first card; the free fit it
    replaced is the free_plateau policy (see docs/ASSUMPTIONS.md).
    """

    def __init__(self, plan: dict | None, use_standards: bool = True, anchor: bool = True):
        self.plan = plan
        self.use_standards = use_standards
        self.anchor = anchor

    def _plate(self, forskolin: float, cells: float) -> dict:
        return {"type": "run_plate", "design": {**self.plan, "forskolin_uM": forskolin, "cells_k": cells}}

    def act(self, obs: dict) -> dict:
        if obs.get("phase") == "start":
            if self.plan is None:
                return {"type": "report", "pEC50": {c: 7.5 for c in obs["compounds"]}}
            if self.plan["forskolin_uM"] == "EC80" or self.plan["cells_k"] == "titrate":
                return {"type": "run_titration", "design": {
                    "cells_k": GI_PILOT_CELLS_K, "stim_min": self.plan["stim_min"], "ibmx": self.plan["ibmx"],
                    "replicates": GI_PILOT_REPLICATES}}
            return self._plate(self.plan["forskolin_uM"], self.plan["cells_k"])
        if obs.get("phase") == "data" and obs["kind"] == "titration":
            return self._plate(*gi_conditions(self.plan, obs["wells"]))
        if obs.get("phase") == "data":
            fit = estimate_pec50_anchored if self.anchor else estimate_pec50
            return {"type": "report", "pEC50": fit(obs["wells"], self.use_standards)}
        return {"type": "report", "pEC50": {}}


def random_gi_plan(rng: np.random.Generator, n_compounds: int = 4) -> dict:
    """A random plan that always fits the labware: transfer volume, replicates and controls give way."""
    plan = {k: v[int(rng.integers(len(v)))] for k, v in GI_CHOICES.items()}
    plan["transfer_ul"] = max(t for t in GI_CHOICES["transfer_ul"]
                              if t <= plan["transfer_ul"] and t * plan["dilution_factor"] <= LADDER_MAX_UL)
    room = 96 - (len(STD_SERIES_NM) if plan["standards"] else 0)  # wells left on the plate
    if n_compounds * plan["n_points"] * plan["replicates"] > room:
        plan["replicates"] = 1
    left = room - n_compounds * plan["n_points"] * plan["replicates"]
    plan["controls_per_side"] = int(min(plan["controls_per_side"], left // 2))
    return plan


def _camp(**changes):
    return lambda seed: GiCampPolicy({**CAMP_EXPERT, **changes})


def _random_camp(seed: int) -> GiCampPolicy:
    rng = np.random.default_rng([seed, 99])
    return GiCampPolicy(random_gi_plan(rng), use_standards=bool(rng.integers(2)), anchor=bool(rng.integers(2)))


GI_CAMP_POTENCY = {
    "expert": _camp(),
    "forskolin_too_high": _camp(forskolin_uM=FSK_MAX_UM),
    # one mistake, no titration plate: both settings the titration would give are handbook numbers instead
    "skip_titration": _camp(forskolin_uM=GI_DEFAULT_FORSKOLIN_UM, cells_k=GI_HANDBOOK_CELLS_K),
    "no_forskolin": _camp(forskolin_uM=0.0),
    "no_ibmx": _camp(ibmx=False),
    "raw_ratio": lambda seed: GiCampPolicy(dict(CAMP_EXPERT), use_standards=False),
    "free_plateau": lambda seed: GiCampPolicy(dict(CAMP_EXPERT), anchor=False),
    "too_many_cells": _camp(cells_k=GI_CELL_RANGE[1]),
    "no_mix": _camp(mix=False),
    "random": _random_camp,
    "null": lambda seed: GiCampPolicy(None),
}
POLICIES["gi_camp_potency"] = GI_CAMP_POTENCY


def run_episode(env, policy) -> tuple[float, dict]:
    """Drive one environment to the end with one policy. Returns (total reward, final info)."""
    obs = env.reset()
    total = 0.0
    while True:
        obs, reward, done, info = env.step(policy.act(obs))
        total += reward
        if done:
            return total, info
