"""Analysis toolkit an agent may call: Z-prime, 4PL fits, standard curves, half-lives.

Numpy only. The 4PL fit is a variable-projection grid: for every (midpoint, slope) pair the
two plateaus are solved exactly, so it is deterministic and has no starting-value failures.
"""
from __future__ import annotations

import numpy as np


def zprime(pos, neg) -> float:
    """Z' = 1 - 3 (sd_pos + sd_neg) / |mean_pos - mean_neg|. NaN with fewer than 3 wells a side."""
    pos, neg = np.asarray(pos, float), np.asarray(neg, float)
    if len(pos) < 3 or len(neg) < 3:
        return float("nan")
    gap = abs(pos.mean() - neg.mean())
    if gap == 0:
        return float("-inf")
    return float(1 - 3 * (pos.std(ddof=1) + neg.std(ddof=1)) / gap)


def fit_4pl(logx, y, step: float = 0.02, pad: float = 1.0):
    """Fit y = bottom + (top - bottom) / (1 + 10 ** (hill * (logec50 - logx))).

    `bottom` is the plateau at low x and `top` at high x, so falling curves work too.
    Returns None when there are fewer than 4 usable points.
    """
    logx, y = np.asarray(logx, float), np.asarray(y, float)
    ok = np.isfinite(logx) & np.isfinite(y)
    logx, y = logx[ok], y[ok]
    if len(y) < 4 or np.ptp(logx) == 0:
        return None
    mids = np.arange(logx.min() - pad, logx.max() + pad + step, step)
    hills = np.arange(0.5, 2.51, 0.1)
    m, h = np.meshgrid(mids, hills, indexing="ij")
    g = 1.0 / (1.0 + 10.0 ** (h[..., None] * (m[..., None] - logx)))
    a11 = ((1 - g) ** 2).sum(-1) + 1e-9
    a22 = (g ** 2).sum(-1) + 1e-9
    a12 = ((1 - g) * g).sum(-1)
    b1 = ((1 - g) * y).sum(-1)
    b2 = (g * y).sum(-1)
    det = a11 * a22 - a12 ** 2
    det = np.where(np.abs(det) < 1e-12, 1e-12, det)
    bottom = (b1 * a22 - b2 * a12) / det
    top = (a11 * b2 - a12 * b1) / det
    sse = ((y - (bottom[..., None] * (1 - g) + top[..., None] * g)) ** 2).sum(-1)
    i, j = np.unravel_index(np.argmin(sse), sse.shape)
    return {"bottom": float(bottom[i, j]), "top": float(top[i, j]), "logec50": float(m[i, j]),
            "hill": float(h[i, j]), "sse": float(sse[i, j])}


def invert_4pl(fit: dict, y):
    """Back-calculate x from y on a fitted curve (standard-curve interpolation)."""
    g = (np.asarray(y, float) - fit["bottom"]) / (fit["top"] - fit["bottom"])
    g = np.clip(g, 0.005, 0.995)
    return 10.0 ** (fit["logec50"] - np.log10(1.0 / g - 1.0) / fit["hill"])


def estimate_pec50(wells: list[dict], use_standards: bool = True) -> dict[str, float]:
    """pEC50 per compound from one dose-response plate observation.

    With use_standards the ratios are first converted to analyte concentration through the
    on-plate standard curve; without it the raw ratio is fitted directly.
    """
    std = [(w["std_nM"], w["ratio"]) for w in wells if w["role"] == "std"]
    curve = None
    if use_standards and len(std) >= 4:
        curve = fit_4pl(np.log10([s for s, _ in std]), [r for _, r in std])
    out: dict[str, float] = {}
    for cid in sorted({w["id"] for w in wells if w["role"] == "cpd"}):
        pts = [(w["conc_nM"], w["ratio"]) for w in wells if w["role"] == "cpd" and w["id"] == cid]
        y = np.array([p[1] for p in pts], float)
        if curve is not None:
            y = invert_4pl(curve, y)
        fit = fit_4pl(np.log10([p[0] for p in pts]), y)
        if fit is not None:
            out[cid] = 9.0 - fit["logec50"]
    return out


def analyte_zprime(standards: list[dict], pos, neg) -> float:
    """Z' in analyte units: controls are back-calculated through the standard curve first.

    This is the window potency is actually read in. A saturated competitive assay can keep a
    good raw-ratio Z' while its back-calculated window has collapsed.
    """
    curve = fit_4pl(np.log10([s["std_nM"] for s in standards]), [s["ratio"] for s in standards])
    if curve is None:
        return float("nan")
    return zprime(invert_4pl(curve, pos), invert_4pl(curve, neg))


def control_zprime(wells: list[dict]) -> float:
    return zprime([w["ratio"] for w in wells if w["role"] == "pos"], [w["ratio"] for w in wells if w["role"] == "neg"])


def response(sample: dict) -> float:
    """MS response: analyte area over internal standard area when one was used."""
    return sample["area"] / sample["is_area"] if sample.get("is_area") else sample["area"]


def estimate_half_life(samples: list[dict], use_lloq: bool = True, lloq_frac: float = 0.02,
                       lo: float = 5.0, hi: float = 1440.0) -> dict[str, float]:
    """Half-life (min) per compound from a stability time course, by log-linear regression.

    Points below lloq_frac of the time-zero mean are dropped when use_lloq is set. If only
    one time point survives, the limit of quantification gives an upper bound, reported as is.
    """
    out: dict[str, float] = {}
    for cid in sorted({s["id"] for s in samples}):
        rows = [s for s in samples if s["id"] == cid]
        t = np.array([s["time_min"] for s in rows], float)
        y = np.array([response(s) for s in rows], float)
        y0 = y[t == t.min()].mean()
        keep = y > (lloq_frac * y0 if use_lloq else 0.0)
        if len(np.unique(t[keep])) < 2:
            later = t[(~keep) & (t > t.min())]
            hl = later.min() / np.log2(1.0 / lloq_frac) if len(later) else hi
        else:
            k = -np.polyfit(t[keep], np.log(y[keep]), 1)[0]
            hl = np.log(2.0) / k if k > 1e-9 else hi
        out[cid] = float(np.clip(hl, lo, hi))
    return out
