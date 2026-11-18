"""Stability brick: design a matrix stability time course (plasma-style) and report half-lives.

Each replicate is its own incubation well. At every time point an aliquot goes into a quench
plate and is read by a simulated LC-MS injection. The reward is the accuracy of the reported
half-lives on a log scale, halved when the run fails its own QC (time-zero replicate CV).
"""
from __future__ import annotations

import numpy as np

from ..analysis import response
from ..chem import remaining
from ..rig import Plate, Rig, RigProfile, usable_wells
from .base import Brick, log_error_score, spec_help, validate

MATRIX_UL = 200.0
SPIKE_UL = 2.0
ALIQUOT_UL = 20.0
QUENCH_UL = 80.0
HORIZON_MIN = 360.0
HL_LO, HL_HI = 5.0, 1440.0  # quantifiable range; answers and truth are clipped to it
INJECTION_SD = 0.18  # log-scale injection-to-injection variability without an internal standard
IS_SD = 0.03
NOISE_FLOOR = 0.004  # detector noise as a fraction of the nominal time-zero area
T0_CV_GATE = 0.15

DESIGN_SPEC = {
    "timepoints_min": ("times", 0.0, HORIZON_MIN, "2 to 8 sampling times in minutes"),
    "replicates": ("int", 1, 3, "incubation wells per compound"),
    "internal_standard": ("bool", None, None, "internal standard in the quench solution"),
    "mix_after_spike": ("bool", None, None, "mix the matrix after spiking the compound"),
    "lid": ("bool", None, None, "lid on the incubation plate at 37 C"),
    "avoid_edges": ("bool", None, None, "keep incubation wells off the outer ring"),
}


def check_design(design, n_compounds: int) -> tuple[dict | None, str | None]:
    clean, err = validate(DESIGN_SPEC, design)
    if err:
        return None, err
    if n_compounds * clean["replicates"] * len(clean["timepoints_min"]) > 96:
        return None, "more than 96 samples; the quench plate is full"
    return clean, None


def simulate_stability(profile: RigProfile, rng: np.random.Generator, half_lives: dict[str, float],
                       design: dict) -> tuple[list[dict], Rig]:
    """Run one validated design. half_lives maps compound id to true half-life in minutes."""
    rig = Rig(profile, rng)
    stock, inc, quench = Plate("stock"), Plate("incubation"), Plate("quench")
    inc_free = iter(usable_wells(design["avoid_edges"]))
    quench_free = iter(usable_wells(False))
    nominal_t0 = 1000.0 * ALIQUOT_UL  # 1 uM in the matrix, as amount in one aliquot
    wells: list[tuple[str, int, tuple]] = []
    for k, cid in enumerate(half_lives):
        rig.fill(stock, (k, 0), 100.0, {cid: 1000.0 * (MATRIX_UL + SPIKE_UL) / SPIKE_UL})
        for rep in range(design["replicates"]):
            w = next(inc_free)
            rig.fill(inc, w, MATRIX_UL)
            rig.transfer(stock, (k, 0), inc, w, SPIKE_UL, new_tip=True, mix_after=design["mix_after_spike"])
            wells.append((cid, rep + 1, w))
    samples: list[dict] = []
    clock = 0.0
    for t in design["timepoints_min"]:
        if t > clock:
            rig.incubate(inc, t - clock, design["lid"])
            for cid, hl in half_lives.items():
                inc.a(cid)[:] *= remaining(hl, t - clock)
            clock = t
        for cid, rep, w in wells:
            q = next(quench_free)
            rig.fill(quench, q, QUENCH_UL, {"IS": 1.0} if design["internal_standard"] else None)
            rig.transfer(inc, w, quench, q, ALIQUOT_UL, new_tip=True, mix_after=True)
            injection = float(np.exp(rng.normal(0.0, INJECTION_SD)))
            area = quench.a(cid)[q] * injection + rng.normal(0.0, NOISE_FLOOR * nominal_t0)
            sample = {"id": cid, "replicate": rep, "time_min": t, "area": round(float(area), 1)}
            if design["internal_standard"]:
                is_area = quench.a("IS")[q] * 100.0 * injection * (1 + rng.normal(0.0, IS_SD))
                sample["is_area"] = round(float(is_area), 1)
            samples.append(sample)
    return samples, rig


def run_cost(n_samples: int, tips: int) -> float:
    return 10.0 + 2.0 * n_samples + 0.1 * tips


def t0_cv(samples: list[dict]) -> dict[str, float]:
    """Time-zero replicate CV per compound; NaN without a time zero or with one replicate."""
    out = {}
    for cid in sorted({s["id"] for s in samples}):
        y = np.array([response(s) for s in samples if s["id"] == cid and s["time_min"] == 0.0])
        out[cid] = float(y.std(ddof=1) / y.mean()) if len(y) >= 2 and y.mean() > 0 else float("nan")
    return out


def score_half_life(estimate, truth: float) -> float:
    if isinstance(estimate, bool) or not isinstance(estimate, (int, float, np.floating)) or not estimate > 0:
        return 0.0
    return log_error_score(float(np.log10(np.clip(estimate, HL_LO, HL_HI))), float(np.log10(np.clip(truth, HL_LO, HL_HI))))


class Stability(Brick):
    name = "stability"
    cost_reference = 300.0

    def __init__(self, seed: int = 0, rig="nominal", n_compounds: int = 3):
        super().__init__(seed, rig)
        self.half_lives = {f"cpd{i + 1}": float(10 ** self._truth_rng.uniform(np.log10(5.0), np.log10(2000.0)))
                           for i in range(n_compounds)}
        self.samples: list[dict] | None = None

    def spec(self) -> dict:
        return {
            "brick": self.name,
            "task": ("Measure the half-life of each compound in a 37 C matrix incubation read by LC-MS. "
                     "Run one time course, then report half-lives in minutes."),
            "actions": {"run": {"design": spec_help(DESIGN_SPEC)},
                        "report": {"half_life_min": "object mapping compound id to half-life in minutes"}},
            "labware": {"matrix_ul": MATRIX_UL, "aliquot_ul": ALIQUOT_UL, "horizon_min": HORIZON_MIN,
                        "quantifiable_range_min": [HL_LO, HL_HI],
                        "nominal_dead_volume_ul": RigProfile().dead_volume_ul},
            "reward": "mean over compounds of max(0, 1 - |log10 half-life error| / 0.5), both clipped to the "
                      "quantifiable range; halved unless there are 2+ replicates, a time zero, and time-zero "
                      "CV <= 15% for every compound; minus 0.15 * cost / 300",
        }

    def reset(self) -> dict:
        return {"phase": "start", "compounds": list(self.half_lives), "spec": self.spec()}

    def step(self, action: dict):
        err = self._check_open(action)
        if err:
            return self._reject(err)
        if action["type"] == "run":
            if self.samples is not None:
                return self._reject("the time course has been run; report")
            design, err = check_design(action.get("design"), len(self.half_lives))
            if err:
                return self._reject(err)
            self.samples, rig = simulate_stability(self.profile, self._noise, self.half_lives, design)
            self.last_log = rig.log
            self.cost += run_cost(len(self.samples), rig.tips)
            return {"phase": "data", "samples": self.samples, "cost": round(self.cost, 2)}, 0.0, False, {}
        if action["type"] == "report":
            answers = action.get("half_life_min")
            answers = answers if isinstance(answers, dict) else {}
            scores = [score_half_life(answers.get(cid), hl) for cid, hl in self.half_lives.items()]
            cvs = t0_cv(self.samples) if self.samples else {}
            qc_pass = bool(cvs) and all(np.isfinite(v) and v <= T0_CV_GATE for v in cvs.values())
            accuracy = float(np.mean(scores))
            reward = accuracy * (1.0 if qc_pass else 0.5) - self._cost_penalty()
            self.done = True
            info = {"accuracy": accuracy, "qc_pass": qc_pass, "t0_cv": cvs, "cost": self.cost,
                    "truth": {k: round(v, 1) for k, v in self.half_lives.items()}}
            return {"phase": "done"}, reward, True, info
        return self._reject(f"unknown action type: {action['type']}")
