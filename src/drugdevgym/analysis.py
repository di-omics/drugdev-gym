"""Analysis toolkit an agent may call: Z-prime, 4PL fits, standard curves, half-lives.

Numpy only. The 4PL fit is a variable-projection grid: for every (midpoint, slope) pair the
two plateaus are solved exactly, so it is deterministic and has no starting-value failures.
"""
from __future__ import annotations

import numpy as np


MIN_R2 = 0.8  # below this a dose-response fit is treated as "no curve"
MIN_SPAN = 0.2  # ... and so is a fitted span under this fraction of the control window


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
    sst = float(((y - y.mean()) ** 2).sum())
    return {"bottom": float(bottom[i, j]), "top": float(top[i, j]), "logec50": float(m[i, j]),
            "hill": float(h[i, j]), "sse": float(sse[i, j]),
            "r2": 1.0 - float(sse[i, j]) / sst if sst > 0 else 0.0,
            "at_edge": bool(i == 0 or i == len(mids) - 1)}


def invert_4pl(fit: dict, y):
    """Back-calculate x from y on a fitted curve (standard-curve interpolation)."""
    span = fit["top"] - fit["bottom"]
    if span == 0:
        return np.full(np.shape(y), np.nan)
    g = np.clip((np.asarray(y, float) - fit["bottom"]) / span, 0.005, 0.995)
    return 10.0 ** (fit["logec50"] - np.log10(1.0 / g - 1.0) / fit["hill"])


def estimate_pec50(wells: list[dict], use_standards: bool = True) -> dict[str, float]:
    """pEC50 per compound from one dose-response plate observation.

    With use_standards the ratios are first converted to analyte concentration through the
    on-plate standard curve; without it the raw ratio is fitted directly. A compound is left
    out when there is no curve to fit: the fit explains too little of the variance (flat
    data), its midpoint sits on the edge of the search range, or its span is a small fraction
    of the control window.
    """
    std = [(w["std_nM"], w["ratio"]) for w in wells if w["role"] == "std"]
    curve = None
    if use_standards and len(std) >= 4:
        curve = fit_4pl(np.log10([s for s, _ in std]), [r for _, r in std])
    def convert(ratios):
        ratios = np.array(ratios, float)
        return invert_4pl(curve, ratios) if curve is not None else ratios

    neg = [w["ratio"] for w in wells if w["role"] == "neg"]
    pos = [w["ratio"] for w in wells if w["role"] == "pos"]
    window = abs(convert(pos).mean() - convert(neg).mean()) if neg and pos else None
    out: dict[str, float] = {}
    for cid in sorted({w["id"] for w in wells if w["role"] == "cpd"}):
        pts = [(w["conc_nM"], w["ratio"]) for w in wells if w["role"] == "cpd" and w["id"] == cid]
        fit = fit_4pl(np.log10([p[0] for p in pts]), convert([p[1] for p in pts]))
        if fit is None or fit["r2"] < MIN_R2 or fit["at_edge"]:
            continue
        if window is not None and abs(fit["top"] - fit["bottom"]) < MIN_SPAN * window:
            continue
        out[cid] = 9.0 - fit["logec50"]
    return out


def analyte_zprime(standards: list[dict], pos, neg) -> float:
    """Z' in analyte units: controls are back-calculated through the standard curve first.

    This is the window potency is actually read in. A saturated competitive assay can keep a
    good raw-ratio Z' while its back-calculated window has collapsed.
    Returns -inf when a control group sits off the top of the standard curve (more than 10% of
    its wells at the saturated end of the fit, or its median above the highest standard): a saturated
    group back-calculates to one clipped value, which would otherwise look like zero spread.
    """
    top_std = max(s["std_nM"] for s in standards)
    curve = fit_4pl(np.log10([s["std_nM"] for s in standards]), [s["ratio"] for s in standards])
    if curve is None:
        return float("nan")
    for group in (pos, neg):
        g = (np.asarray(group, float) - curve["bottom"]) / (curve["top"] - curve["bottom"])
        if np.mean(g >= 0.995) > 0.1 or float(np.median(invert_4pl(curve, group))) > top_std:
            return float("-inf")
    return zprime(invert_4pl(curve, pos), invert_4pl(curve, neg))


def control_zprime(wells: list[dict]) -> float:
    return zprime([w["ratio"] for w in wells if w["role"] == "pos"], [w["ratio"] for w in wells if w["role"] == "neg"])


def response(sample: dict) -> float:
    """MS response: analyte area over internal standard area when one was used."""
    if sample.get("is_area") is None:
        return sample["area"]
    return sample["area"] / sample["is_area"] if sample["is_area"] else float("nan")


def estimate_half_life(samples: list[dict], use_lloq: bool = True, lloq_frac: float = 0.02,
                       lo: float = 5.0, hi: float = 1440.0) -> dict[str, float]:
    """Half-life (min) per compound from a stability time course, by log-linear regression.

    Points below lloq_frac of the time-zero mean are dropped when use_lloq is set. If only
    one time point survives, the limit of quantification gives an upper bound, reported as is.
    Non-finite responses are ignored, and a compound with fewer than two time points is left out.
    """
    out: dict[str, float] = {}
    for cid in sorted({s["id"] for s in samples}):
        rows = [s for s in samples if s["id"] == cid]
        t = np.array([s["time_min"] for s in rows], float)
        y = np.array([response(s) for s in rows], float)
        t, y = t[np.isfinite(y)], y[np.isfinite(y)]
        if len(np.unique(t)) < 2:
            continue  # no time course: leave the compound unanswered rather than guess
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
