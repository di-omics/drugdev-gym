"""Analysis toolkit an agent may call: Z-prime, 4PL fits, standard curves, half-lives.

Numpy only. The 4PL fit is a variable-projection grid: for every (midpoint, slope) pair the
two plateaus are solved exactly, so it is deterministic and has no starting-value failures.
"""
from __future__ import annotations

import numpy as np

from .chem import enzyme_exposure


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
        if y0 <= 0:
            continue  # no usable time zero
        keep = y > (lloq_frac * y0 if use_lloq else 0.0)
        if len(np.unique(t[keep])) < 2:
            later = t[(~keep) & (t > t.min())]
            hl = later.min() / np.log2(1.0 / lloq_frac) if len(later) else hi
        else:
            k = -np.polyfit(t[keep], np.log(y[keep]), 1)[0]
            hl = np.log(2.0) / k if k > 1e-9 else hi
        out[cid] = float(np.clip(hl, lo, hi))
    return out


def estimate_fluid_half_life(samples: list[dict], activity_half_life_min: dict[str, float] | None = None,
                             relative_activity: dict[str, float] | None = None, use_lloq: bool = True,
                             lloq_frac: float = 0.02, lo: float = 1.0, hi: float = 960.0) -> dict[str, dict[str, float]]:
    """Half-life (min) per compound and per fluid from a time course in enzymatic fluids.

    Every sample names its "fluid"; each compound is fitted in each fluid by log-linear regression.
    activity_half_life_min (fluid -> minutes) corrects for an enzyme that loses activity during
    the run: the clock is converted to minutes of digestion at the starting activity
    (chem.enzyme_exposure) before the fit, so a slow compound is not read too long because the
    enzyme was weaker at the late time points.
    relative_activity (fluid -> fraction, from enzyme_activity) then scales the result to
    specification activity: a digest at 70% activity takes 1/0.7 times as long. A fluid with no
    usable entry in either mapping is left uncorrected. The clip to lo..hi is applied last.

    The limit of quantification is applied per time point, not per injection: the first time
    point whose replicate mean is under lloq_frac of the time-zero mean is dropped with every
    later one. Keeping single replicates that happen to read above the limit would flatten the
    slope. If only time zero survives, the dropped time point gives an upper bound, reported as is.
    A compound with no time course or no usable time zero is left out.
    """
    out: dict[str, dict[str, float]] = {}
    for fluid in sorted({s["fluid"] for s in samples}):
        decay = (activity_half_life_min or {}).get(fluid)
        scale = (relative_activity or {}).get(fluid)
        scale = float(scale) if scale is not None and np.isfinite(scale) and scale > 0 else 1.0
        in_fluid = [s for s in samples if s["fluid"] == fluid]
        for cid in sorted({s["id"] for s in in_fluid}):
            rows = [s for s in in_fluid if s["id"] == cid]
            t = np.array([enzyme_exposure(0.0, s["time_min"], decay) for s in rows], float)
            y = np.array([response(s) for s in rows], float)
            t, y = t[np.isfinite(y)], y[np.isfinite(y)]
            if len(np.unique(t)) < 2:
                continue
            y0 = y[t == t.min()].mean()
            if y0 <= 0:
                continue
            keep = y > 0
            if use_lloq:
                low = [u for u in np.unique(t) if u > t.min() and y[t == u].mean() <= lloq_frac * y0]
                if low:
                    keep &= t < low[0]
            if len(np.unique(t[keep])) < 2:
                later = t[(~keep) & (t > t.min())]
                hl = later.min() / np.log2(1.0 / lloq_frac) if len(later) else float("inf")
            else:
                k = -np.polyfit(t[keep], np.log(y[keep]), 1)[0]
                hl = np.log(2.0) / k if k > 1e-9 else float("inf")
            out.setdefault(cid, {})[fluid] = float(np.clip(hl * scale, lo, hi))
    return out


def enzyme_activity(samples: list[dict], control_id: str, reference_half_life_min: float,
                    activity_half_life_min: dict[str, float] | None = None, use_lloq: bool = True,
                    lloq_frac: float = 0.02) -> dict[str, float]:
    """Relative enzyme activity per fluid, read from a control substrate run alongside the compounds.

    Activity = the control's reference half-life (at specification activity) over the half-life it
    showed in this run. 1.0 is a fluid at specification, 0.5 an enzyme that has lost half its
    activity, 0.0 a control that did not move. A fluid with no usable control time course is left out.
    """
    control = [s for s in samples if s["id"] == control_id]
    fit = estimate_fluid_half_life(control, activity_half_life_min, None, use_lloq, lloq_frac, lo=0.0, hi=float("inf"))
    return {fluid: float(reference_half_life_min / hl) for fluid, hl in fit.get(control_id, {}).items() if hl > 0}


def standard_curve(wells: list[dict]) -> dict | None:
    """4PL through the on-plate standards (role "std") of one plate observation; None with fewer than 4."""
    std = [(w["std_nM"], w["ratio"]) for w in wells if w["role"] == "std"]
    if len(std) < 4:
        return None
    return fit_4pl(np.log10([s for s, _ in std]), [r for _, r in std])


def ec_fraction(ec50: float, hill: float = 1.0, fraction: float = 0.8) -> float:
    """Concentration that gives `fraction` of the maximal effect on a Hill curve (EC80 = EC50 * 4 ** (1 / hill)).

    Raises ValueError unless 0 < fraction < 1 and hill > 0: the curve never reaches 0 or 1, and a
    flat curve has no such concentration.
    """
    if not 0.0 < fraction < 1.0:
        raise ValueError("fraction must lie strictly between 0 and 1")
    if not hill > 0.0:
        raise ValueError("hill must be positive")
    return float(ec50 * (fraction / (1.0 - fraction)) ** (1.0 / hill))


def titration_curve(wells: list[dict]) -> dict | None:
    """Stimulant dose-response in analyte units from one titration plate observation.

    Wells with role "fsk" carry forskolin_uM and a ratio. Ratios are converted to analyte
    through the on-plate standard curve, then analyte against log10 concentration is fitted
    with a rising 4PL. Returns None without a standard curve. Otherwise returns kit_mid_nM
    (the midpoint of the standard curve), top_wells_nM (mean analyte at the highest
    concentration tested), unstim_nM (mean analyte of the role "unstim" wells, or None) and
    fitted (whether a rising curve could be placed). When fitted, ec50_uM, hill, bottom_nM
    and top_nM describe it; they are None when the data are flat, fall instead of rise, or
    put the midpoint on the edge of the search range.
    """
    curve = standard_curve(wells)
    if curve is None:
        return None
    pts = [(w["forskolin_uM"], w["ratio"]) for w in wells if w["role"] == "fsk" and w["forskolin_uM"] > 0]
    unstim = [w["ratio"] for w in wells if w["role"] == "unstim"]
    out = {"kit_mid_nM": float(10.0 ** curve["logec50"]), "top_wells_nM": None,
           "unstim_nM": float(np.mean(invert_4pl(curve, unstim))) if unstim else None,
           "fitted": False, "ec50_uM": None, "hill": None, "bottom_nM": None, "top_nM": None}
    if not pts:
        return out
    conc = np.array([p[0] for p in pts], float)
    analyte = invert_4pl(curve, [p[1] for p in pts])
    out["top_wells_nM"] = float(analyte[conc == conc.max()].mean())
    fit = fit_4pl(np.log10(conc), analyte)
    if fit is None or fit["r2"] < MIN_R2 or fit["at_edge"] or fit["top"] <= fit["bottom"]:
        return out
    out.update(fitted=True, ec50_uM=float(10.0 ** fit["logec50"]), hill=fit["hill"], bottom_nM=fit["bottom"],
               top_nM=fit["top"])
    return out


def titration_level(curve: dict, forskolin_uM: float) -> float | None:
    """Analyte (nM) a titration_curve result predicts at one concentration, at the titration's cell number.

    Zero concentration returns the unstimulated wells. Without a fitted curve the mean of the
    highest concentration tested is the only estimate there is. None when neither exists.
    """
    if forskolin_uM <= 0:
        return curve["unstim_nM"] if curve["unstim_nM"] is not None else curve["bottom_nM"]
    if not curve["fitted"]:
        return curve["top_wells_nM"]
    x = (forskolin_uM / curve["ec50_uM"]) ** curve["hill"]
    return float(curve["bottom_nM"] + (curve["top_nM"] - curve["bottom_nM"]) * x / (1.0 + x))


def estimate_pec50_anchored(wells: list[dict], use_standards: bool = True, anchor_logs: float = 2.0) -> dict[str, float]:
    """pEC50 per compound with the vehicle controls pinning the zero-dose plateau.

    The conversion and the no-curve rules are those of estimate_pec50. The difference is that
    the vehicle control wells (role "neg") join every compound's curve as zero-dose points,
    placed anchor_logs below its lowest concentration. In an inhibition assay the zero-dose
    plateau is the high, noisy end of the curve and a potent compound may not reach it on the
    ladder, so a free plateau lets the midpoint drift. A midpoint that lands below the lowest
    concentration tested sits in the gap between the anchor and the ladder, where no well
    constrains it, so that compound is left out. Without vehicle controls this is
    estimate_pec50.
    """
    curve = standard_curve(wells) if use_standards else None

    def convert(ratios):
        ratios = np.array(ratios, float)
        return invert_4pl(curve, ratios) if curve is not None else ratios

    neg = convert([w["ratio"] for w in wells if w["role"] == "neg"])
    pos = convert([w["ratio"] for w in wells if w["role"] == "pos"])
    window = abs(pos.mean() - neg.mean()) if len(neg) and len(pos) else None
    out: dict[str, float] = {}
    for cid in sorted({w["id"] for w in wells if w["role"] == "cpd"}):
        pts = [(w["conc_nM"], w["ratio"]) for w in wells if w["role"] == "cpd" and w["id"] == cid]
        logx, y = np.log10([p[0] for p in pts]), convert([p[1] for p in pts])
        lowest = logx.min()
        if len(neg):
            logx = np.concatenate([np.full(len(neg), lowest - anchor_logs), logx])
            y = np.concatenate([neg, y])
        fit = fit_4pl(logx, y)
        if fit is None or fit["r2"] < MIN_R2 or fit["at_edge"]:
            continue
        if len(neg) and fit["logec50"] < lowest:
            continue
        if window is not None and abs(fit["top"] - fit["bottom"]) < MIN_SPAN * window:
            continue
        out[cid] = 9.0 - fit["logec50"]
    return out
