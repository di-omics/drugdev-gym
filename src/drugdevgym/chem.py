"""Chemistry layer: the textbook models every environment is built on.

Receptor response is a Hill curve. Detection is a competitive TR-FRET immunoassay, where
the signal ratio FALLS as the analyte rises (for a Gq assay the analyte is IP1, the
accumulated metabolite of the second messenger IP3), so raw ratios are not
proportional to response. Degradation is first order.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class Kit:
    """Competitive TR-FRET immunoassay kit. k is the analyte concentration at half displacement."""

    k: float = 30.0
    rmin: float = 600.0
    rmax: float = 6000.0
    slope: float = 1.0

    def ratio(self, analyte):
        m = np.maximum(np.asarray(analyte, float), 0.0)
        return self.rmin + (self.rmax - self.rmin) / (1.0 + (m / self.k) ** self.slope)


def hill(conc, ec50: float, n: float = 1.0):
    """Fractional occupancy-style response in [0, 1]."""
    x = (np.maximum(np.asarray(conc, float), 0.0) / ec50) ** n
    return x / (1.0 + x)


def remaining(half_life_min: float, minutes: float) -> float:
    """First-order decay: fraction left after `minutes`."""
    return float(0.5 ** (minutes / half_life_min))


def q10_factor(temp_c: float, ref_c: float = 37.0, q10: float = 2.0) -> float:
    """Q10 rule: rate at temp_c relative to the rate at ref_c, q10 ** ((temp_c - ref_c) / 10)."""
    return float(q10 ** ((temp_c - ref_c) / 10.0))


def acid_catalysis_factor(ph: float, ref_ph: float) -> float:
    """Specific acid catalysis, first order in [H+]: rate at ph relative to the rate at ref_ph."""
    return float(10.0 ** (ref_ph - ph))


def enzyme_exposure(t0: float, t1: float, activity_half_life_min: float | None = None) -> float:
    """Minutes of digestion at the time-zero activity that the clock interval t0..t1 is worth.

    The enzyme loses activity first order, a(t) = 0.5 ** (t / activity_half_life_min), and a substrate
    far below Km is digested at a rate proportional to a(t), so its loss depends on the integral of
    a(t). None means a stable enzyme: the exposure is the clock time.
    """
    if activity_half_life_min is None or not np.isfinite(activity_half_life_min):
        return float(t1 - t0)
    lam = np.log(2.0) / activity_half_life_min
    return float((np.exp(-lam * t0) - np.exp(-lam * t1)) / lam)


def accumulated(loss_per_min: float, minutes: float) -> float:
    """Zero-order production against first-order loss: the minutes of production still present.

    Solves dC/dt = 1 - k C from C(0) = 0. With no loss this is the elapsed time; with fast
    loss it plateaus at the steady state 1 / k.
    """
    if loss_per_min <= 0:
        return float(minutes)
    return float(-np.expm1(-loss_per_min * minutes) / loss_per_min)


def uninhibited_fraction(inhibitor: float, ki: float) -> float:
    """Enzyme activity left under a competitive inhibitor with the substrate far below its Km: 1 / (1 + [I] / Ki)."""
    return float(1.0 / (1.0 + max(float(inhibitor), 0.0) / ki))


def survival(conc, half_conc: float, n: float = 1.0):
    """Fraction of cells surviving an exposure: 1 minus a Hill response, 0.5 at `half_conc`."""
    return 1.0 - hill(conc, half_conc, n)


def haber_shift(hours: float, reference_hours: float, exponent: float) -> float:
    """Haber-type exposure rule (C ** n * t = constant): fold rise of the half-effect concentration
    when the exposure lasts `hours` instead of `reference_hours`. exponent = 1 / n; 0 = no time dependence."""
    return float((reference_hours / hours) ** exponent)


def dissolved(conc, solubility: float):
    """Concentration actually in solution: the total, capped at the solubility limit."""
    return np.minimum(np.maximum(np.asarray(conc, float), 0.0), solubility)


def grown(cells, hours: float, doubling_h: float, capacity: float):
    """Exponential growth up to a carrying capacity: a confluent well cannot hold more cells."""
    return np.minimum(np.maximum(np.asarray(cells, float), 0.0) * 2.0 ** (hours / doubling_h), capacity)
