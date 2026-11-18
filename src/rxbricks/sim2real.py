"""Sim-to-real hooks. Nothing here has touched hardware; these are the joints built for it.

1. worklist_csv: the steps an agent's protocol expands to, as a robot-neutral worklist.
2. RigProfile.from_qc (in rig.py): fit the simulated pipetting error to real QC data.
3. fingerprint + gap: summarise a real control plate, and measure where it falls inside the
   distribution the simulator predicts for the same protocol. That number is the gap.
"""
from __future__ import annotations

import csv
import io

import numpy as np

from .analysis import zprime
from .rig import EDGE, ROW_NAMES

WORKLIST_FIELDS = ["step", "op", "src_plate", "src_well", "dst_plate", "dst_well", "ul", "new_tip", "mix_after",
                   "what", "minutes", "lid"]


def worklist_csv(log: list[dict]) -> str:
    """Turn a brick's step log (brick.last_log) into a CSV worklist."""
    buf = io.StringIO()
    writer = csv.DictWriter(buf, fieldnames=WORKLIST_FIELDS, extrasaction="ignore", lineterminator="\n")
    writer.writeheader()
    for i, step in enumerate(log, 1):
        writer.writerow({"step": i, **step})
    return buf.getvalue()


def _cv(x) -> float:
    x = np.asarray(x, float)
    return float(x.std(ddof=1) / x.mean()) if len(x) >= 2 and x.mean() != 0 else float("nan")


def fingerprint(wells: list[dict]) -> dict:
    """QC statistics of the control wells on one plate observation (simulated or real)."""
    neg = [w for w in wells if w["role"] == "neg"]
    pos = [w for w in wells if w["role"] == "pos"]
    n, p = [w["ratio"] for w in neg], [w["ratio"] for w in pos]
    out = {"neg_cv": _cv(n), "pos_cv": _cv(p), "zprime": zprime(p, n),
           "window": float(np.mean(n) / np.mean(p)) if n and p else float("nan")}

    def is_edge(w) -> bool:
        return bool(EDGE[ROW_NAMES.index(w["well"][0]), int(w["well"][1:]) - 1])

    edge, inner = [w["ratio"] for w in neg if is_edge(w)], [w["ratio"] for w in neg if not is_edge(w)]
    out["edge_over_inner"] = float(np.mean(edge) / np.mean(inner)) if edge and inner else float("nan")
    rows = np.array([ROW_NAMES.index(w["well"][0]) for w in neg], float)
    if len(set(rows)) >= 3:
        out["row_slope"] = float(np.polyfit(rows, np.array(n) / np.mean(n), 1)[0])
    else:
        out["row_slope"] = float("nan")
    return out


def gap(real: dict, simulated: list[dict], z_limit: float = 2.0) -> dict:
    """Where each real QC statistic sits in the simulated distribution, as a z-score.

    covered = |z| <= z_limit for every statistic both sides could compute. An uncovered
    statistic names the physics the simulator is missing or has mis-sized.
    """
    report, covered = {}, True
    for key, value in real.items():
        sims = np.array([s[key] for s in simulated if np.isfinite(s.get(key, np.nan))], float)
        if not np.isfinite(value) or len(sims) < 5:
            report[key] = {"real": value, "z": None, "covered": None}
            continue
        sd = sims.std(ddof=1)
        z = float((value - sims.mean()) / sd) if sd > 0 else float("inf")
        ok = abs(z) <= z_limit
        covered = covered and ok
        report[key] = {"real": float(value), "sim_mean": float(sims.mean()), "sim_sd": float(sd), "z": z, "covered": ok}
    return {"covered": covered, "stats": report}
