"""Scripted reviewers for the review environments: plain QC rules, written down.

These are the reference level a human expert or a model is compared against. Each feature is
something a reviewer can compute from the raw data and the known reference value alone. A
judge makes two separate calls: the cause (which signature does the data carry?) and the
verdict (does it matter: is the cause one that usually does harm, or does the run fail its QC,
lose a compound, or misread the reference by more than a limit?). The thresholds and the list
of serious causes were set on seeds 0 to 399; report on other seeds.
"""
from __future__ import annotations

import numpy as np

from .analysis import control_zprime, estimate_half_life, estimate_pec50, fit_4pl, invert_4pl, response
from .envs.dose_response import ZPRIME_GATE
from .envs.stability import T0_CV_GATE, t0_cv
from .rig import EDGE, ROW_NAMES


def _pos(w: dict) -> tuple[int, int]:
    return ROW_NAMES.index(w["well"][0]), int(w["well"][1:]) - 1


def plate_features(wells: list[dict], reference: dict | None = None) -> dict:
    """window: negative over positive control signal. row_slope: fractional drift per row in the
    negative controls. replicate_gap: median relative disagreement between duplicate wells.
    edge_gap: mean signal excess of an outer-ring well over its inner duplicate. reference_shift:
    measured minus known pEC50 of the reference. reference_top: how far the reference's top dose
    moves the signal across the control window (1 = all the way)."""
    neg = [w for w in wells if w["role"] == "neg"]
    pos = [w for w in wells if w["role"] == "pos"]
    std = [w for w in wells if w["role"] == "std"]
    cpd = [w for w in wells if w["role"] == "cpd"]
    neg_mean, pos_mean = float(np.mean([w["ratio"] for w in neg])), float(np.mean([w["ratio"] for w in pos]))
    out = {"zprime": control_zprime(wells), "window": neg_mean / pos_mean}
    rows = np.array([_pos(w)[0] for w in neg], float)
    signal = np.array([w["ratio"] for w in neg], float)
    out["row_slope"] = float(np.polyfit(rows, signal / signal.mean(), 1)[0])
    pairs: dict[tuple, list[dict]] = {}
    for w in cpd:
        pairs.setdefault((w["id"], w["conc_nM"]), []).append(w)
    gaps, edge_gaps = [], []
    for members in pairs.values():
        if len(members) != 2:
            continue
        a, b = members
        gaps.append(abs(a["ratio"] - b["ratio"]) / ((a["ratio"] + b["ratio"]) / 2))
        on_edge = [bool(EDGE[_pos(m)]) for m in members]
        if on_edge[0] != on_edge[1]:  # exactly one of the two duplicates sits on the outer ring
            edge, inner = (a, b) if on_edge[0] else (b, a)
            edge_gaps.append((edge["ratio"] - inner["ratio"]) / inner["ratio"])
    out["replicate_gap"] = float(np.median(gaps)) if gaps else float("nan")
    out["edge_gap"] = float(np.mean(edge_gaps)) if edge_gaps else 0.0
    out["reference_shift"], out["reference_top"] = float("nan"), float("nan")
    if reference:
        curve = fit_4pl(np.log10([w["std_nM"] for w in std]), [w["ratio"] for w in std]) if len(std) >= 4 else None
        pts = [(w["conc_nM"], w["ratio"]) for w in cpd if w["id"] == reference["id"]]
        top = max(p[0] for p in pts)
        top_mean = float(np.mean([p[1] for p in pts if p[0] == top]))
        if neg_mean != pos_mean:
            out["reference_top"] = (neg_mean - top_mean) / (neg_mean - pos_mean)
        y = np.array([p[1] for p in pts], float)
        fit = fit_4pl(np.log10([p[0] for p in pts]), invert_4pl(curve, y) if curve is not None else y)
        if fit is not None and fit["r2"] >= 0.8 and not fit["at_edge"]:
            out["reference_shift"] = 9.0 - fit["logec50"] - reference["pEC50"]
    return out


class PlateJudge:
    """QC rules for a dose-response plate."""

    SHIFT_LIMIT = 0.4  # a reference misread by more than this is taken to mean the plate misleads
    SERIOUS = {"weak_cells", "reader_drift", "starved_ladder", "unmixed_ladder", "edge_effect"}  # usually do harm

    def cause(self, f: dict) -> str:
        if f["window"] < 1.6:
            return "weak_cells"
        if abs(f["row_slope"]) > 0.03:
            return "reader_drift"
        if f["replicate_gap"] > 0.2 or f["reference_top"] < 0.75:
            return "starved_ladder"
        if f["edge_gap"] > 0.18:
            return "edge_effect"
        if f["reference_shift"] < -0.3:
            return "unmixed_ladder"  # the reference reads less potent than it is
        if f["reference_shift"] > 0.08:
            return "tip_carryover"  # the reference reads more potent than it is
        return "none"

    def harm(self, obs: dict, f: dict) -> bool:
        cpds = {w["id"] for w in obs["wells"] if w["role"] == "cpd"}
        fitted = estimate_pec50(obs["wells"], use_standards=True)
        return (not np.isfinite(f["zprime"]) or f["zprime"] < ZPRIME_GATE or set(fitted) != cpds
                or not np.isfinite(f["reference_shift"]) or abs(f["reference_shift"]) > self.SHIFT_LIMIT)

    def act(self, obs: dict) -> dict:
        f = plate_features(obs["wells"], obs.get("reference"))
        cause = self.cause(f)
        repeat = cause in self.SERIOUS or self.harm(obs, f)
        return {"type": "judge", "verdict": "repeat" if repeat else "accept", "cause": cause}


def timecourse_features(samples: list[dict], reference: dict | None = None) -> dict:
    """t0_deficit: how far time zero reads below the line through the later points (log units).
    replicate_gap: median log disagreement between duplicate incubations.
    reference_lift: how far the control compound's later points sit above (+) or below (-) the
    decay its known half-life predicts. reference_error: log10 of measured over known half-life."""
    deficits, gaps, lift, ref_error = [], [], float("nan"), float("nan")
    estimates = estimate_half_life(samples)
    for cid in sorted({s["id"] for s in samples}):
        rows = [s for s in samples if s["id"] == cid]
        t = np.array([s["time_min"] for s in rows], float)
        y = np.array([response(s) for s in rows], float)
        times = np.unique(t)
        mean = {x: float(np.mean(y[t == x])) for x in times}
        y0 = mean[times[0]]
        if not y0 > 0:
            continue
        for x in times:
            pair = y[t == x]
            if len(pair) == 2 and (pair > 0.02 * y0).all():
                gaps.append(abs(float(np.log(pair[0] / pair[1]))))
        later = (t >= 30.0) & (y > 0.02 * y0)
        if len(np.unique(t[later])) >= 2:
            _slope, intercept = np.polyfit(t[later], np.log(y[later]), 1)
            deficits.append(float(intercept - np.log(y0)))
        if reference and cid == reference["id"]:
            expected = {x: 0.5 ** (x / reference["half_life_min"]) for x in times}
            # choose points by where the control SHOULD still be readable, not by what was observed
            late = [x for x in times[len(times) // 2:] if expected[x] > 0.02]
            if late:
                lift = float(np.mean([np.log(max(mean[x], 1e-9 * y0) / y0 / expected[x]) for x in late]))
            if cid in estimates:
                ref_error = float(np.log10(estimates[cid] / reference["half_life_min"]))
    return {"t0_deficit": max(deficits, default=float("nan")),
            "replicate_gap": float(np.median(gaps)) if gaps else float("nan"),
            "reference_lift": lift, "reference_error": ref_error, "answered": len(estimates),
            "compounds": len({s["id"] for s in samples})}


class TimecourseJudge:
    """QC rules for a stability time course."""

    ERROR_LIMIT = 0.3  # the control compound's half-life read more than 2-fold off
    SERIOUS = {"internal_standard_failure", "short_aliquots", "unmixed_spike"}  # usually do harm

    def cause(self, f: dict) -> str:
        if f["replicate_gap"] > 0.15:
            return "internal_standard_failure"
        if f["reference_lift"] < -0.1:
            return "short_aliquots"  # the control compound disappears faster than it can decay
        if f["t0_deficit"] > 0.2:
            return "unmixed_spike"
        if f["reference_lift"] > 0.06:
            return "evaporation"  # the control compound decays slower than it can
        return "none"

    def harm(self, obs: dict, f: dict) -> bool:
        cvs = t0_cv(obs["samples"])
        gate_fail = any(not np.isfinite(v) or v > T0_CV_GATE for v in cvs.values())
        return (gate_fail or f["answered"] != f["compounds"] or not np.isfinite(f["reference_error"])
                or abs(f["reference_error"]) > self.ERROR_LIMIT)

    def act(self, obs: dict) -> dict:
        f = timecourse_features(obs["samples"], obs.get("reference"))
        cause = self.cause(f)
        repeat = cause in self.SERIOUS or self.harm(obs, f)
        return {"type": "judge", "verdict": "repeat" if repeat else "accept", "cause": cause}


class ForcedVerdict:
    """A judge that looks at the data for the cause but always gives the same verdict.

    This is the single mistake of never questioning data (accept) or trusting nothing (repeat)."""

    def __init__(self, judge, verdict: str):
        self.judge, self.verdict = judge, verdict

    def act(self, obs: dict) -> dict:
        return {**self.judge.act(obs), "verdict": self.verdict}


class FixedJudge:
    """A constant call made without looking at the data."""

    def __init__(self, verdict: str, cause: str):
        self.action = {"type": "judge", "verdict": verdict, "cause": cause}

    def act(self, obs: dict) -> dict:
        return self.action


class RandomJudge:
    def __init__(self, seed: int, causes):
        rng = np.random.default_rng([seed, 99])
        self.action = {"type": "judge", "verdict": str(rng.choice(["accept", "repeat"])), "cause": str(rng.choice(causes))}

    def act(self, obs: dict) -> dict:
        return self.action
