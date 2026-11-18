"""Scripted reviewers for the review environments: plain QC rules, written down.

These are the reference level a human expert or a model is compared against. Each feature is
something a reviewer can compute from the raw data alone.
"""
from __future__ import annotations

import numpy as np

from .analysis import control_zprime, fit_4pl, invert_4pl, response
from .rig import ROW_NAMES


def _row(w: dict) -> int:
    return ROW_NAMES.index(w["well"][0])


def plate_features(wells: list[dict], reference: dict | None = None) -> dict:
    neg = [w for w in wells if w["role"] == "neg"]
    pos = [w for w in wells if w["role"] == "pos"]
    std = [w for w in wells if w["role"] == "std"]
    cpd = [w for w in wells if w["role"] == "cpd"]
    out = {"zprime": control_zprime(wells),
           "window": float(np.mean([w["ratio"] for w in neg]) / np.mean([w["ratio"] for w in pos]))}
    rows = np.array([_row(w) for w in neg], float)
    signal = np.array([w["ratio"] for w in neg], float)
    out["row_slope"] = float(np.polyfit(rows, signal / signal.mean(), 1)[0])  # fractional change per row
    pairs: dict[tuple, list[dict]] = {}
    for w in cpd:
        pairs.setdefault((w["id"], w["conc_nM"]), []).append(w)
    gaps, edge_gaps = [], []
    for members in pairs.values():
        if len(members) != 2:
            continue
        a, b = members
        gaps.append(abs(a["ratio"] - b["ratio"]) / ((a["ratio"] + b["ratio"]) / 2))
        on_edge = [_row(m) in (0, 7) for m in members]
        if on_edge[0] != on_edge[1]:  # one replicate on row A or H, its partner one row in
            edge, inner = (a, b) if on_edge[0] else (b, a)
            edge_gaps.append((edge["ratio"] - inner["ratio"]) / inner["ratio"])
    out["replicate_gap"] = float(np.median(gaps)) if gaps else float("nan")
    out["edge_gap"] = float(np.mean(edge_gaps)) if edge_gaps else 0.0
    out["reference_shift"] = float("nan")  # measured minus known pEC50 of the reference compound
    if reference:
        curve = fit_4pl(np.log10([w["std_nM"] for w in std]), [w["ratio"] for w in std]) if len(std) >= 4 else None
        pts = [(w["conc_nM"], w["ratio"]) for w in cpd if w["id"] == reference["id"]]
        y = np.array([p[1] for p in pts], float)
        fit = fit_4pl(np.log10([p[0] for p in pts]), invert_4pl(curve, y) if curve is not None else y)
        if fit is not None and fit["r2"] >= 0.8 and not fit["at_edge"]:
            out["reference_shift"] = 9.0 - fit["logec50"] - reference["pEC50"]
    return out


class PlateJudge:
    """QC rules for a dose-response plate, checked in order of how unambiguous they are."""

    def act(self, obs: dict) -> dict:
        f = plate_features(obs["wells"], obs.get("reference"))
        cause = "none"
        if f["window"] < 1.6:
            cause = "weak_cells"
        elif abs(f["row_slope"]) > 0.03:
            cause = "reader_drift"
        elif f["replicate_gap"] > 0.2:
            cause = "starved_ladder"
        elif f["edge_gap"] > 0.18:
            cause = "edge_effect"
        elif f["reference_shift"] < -0.38:
            cause = "unmixed_ladder"  # the reference reads less potent than it is
        elif f["reference_shift"] > 0.1:
            cause = "tip_carryover"  # the reference reads more potent than it is
        failed_qc = not np.isfinite(f["zprime"]) or f["zprime"] < 0.5 or not np.isfinite(f["reference_shift"])
        return {"type": "judge", "verdict": "repeat" if cause != "none" or failed_qc else "accept", "cause": cause}


def timecourse_features(samples: list[dict], reference: dict | None = None) -> dict:
    """t0_deficit: how far time zero reads below the line through the later points (log units).
    replicate_gap: median log disagreement between duplicate incubations.
    reference_lift: how far the control compound's later points sit above (+) or below (-) the
    decay its known half-life predicts."""
    deficits, gaps, lift = [], [], float("nan")
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
            late = [x for x in times[len(times) // 2:] if mean[x] > 0.02 * y0]
            if late:
                lift = float(np.mean([np.log(mean[x] / y0 / 0.5 ** (x / reference["half_life_min"])) for x in late]))
    return {"t0_deficit": max(deficits, default=float("nan")),
            "replicate_gap": float(np.median(gaps)) if gaps else float("nan"), "reference_lift": lift}


class TimecourseJudge:
    """QC rules for a stability time course."""

    def act(self, obs: dict) -> dict:
        f = timecourse_features(obs["samples"], obs.get("reference"))
        cause = "none"
        if f["replicate_gap"] > 0.15:
            cause = "internal_standard_failure"
        elif f["reference_lift"] < -0.1:
            cause = "short_aliquots"  # the control compound disappears faster than it can decay
        elif f["t0_deficit"] > 0.2:
            cause = "unmixed_spike"
        elif f["reference_lift"] > 0.12:
            cause = "evaporation"  # the control compound decays slower than it can
        unreadable = not np.isfinite(f["reference_lift"])
        return {"type": "judge", "verdict": "repeat" if cause != "none" or unreadable else "accept", "cause": cause}


class FixedJudge:
    """Always the same call. always_accept never questions data; always_repeat trusts nothing."""

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
