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
