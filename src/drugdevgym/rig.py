"""Physical layer: plates, a liquid handler and a plate reader with explicit error sources.

Every environment runs its protocol through this layer, so an agent's choices (volumes, mixing,
tips, lids, layout) change the data the same way they would at the bench.

Every number in RigProfile is an assumed order-of-magnitude default, NOT a measurement.
Replace the pipetting terms with your own QC through RigProfile.from_qc before trusting a
policy on hardware. See docs/ASSUMPTIONS.md.
"""
from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np

ROWS, COLS = 8, 12
ROW_NAMES = "ABCDEFGH"
EDGE = np.zeros((ROWS, COLS), dtype=bool)
EDGE[0, :] = EDGE[-1, :] = True
EDGE[:, 0] = EDGE[:, -1] = True
MIX_BY_DIFFUSION_MIN = 20.0


def well_name(w) -> str:
    r, c = w
    return f"{ROW_NAMES[r]}{c + 1}"


def usable_wells(avoid_edges: bool) -> list[tuple[int, int]]:
    """Wells in column-major order, optionally skipping the outer ring."""
    return [(r, c) for c in range(COLS) for r in range(ROWS) if not (avoid_edges and EDGE[r, c])]


@dataclass(frozen=True)
class RigProfile:
    """One liquid handler plus incubator plus reader. These are the sim-to-real parameters."""

    name: str = "nominal"
    cv_floor: float = 0.01  # random pipetting CV at large volume
    cv_k: float = 0.04  # extra CV times uL: cv(v) = cv_floor + cv_k / v
    bias_k: float = -0.03  # bias times uL: bias(v) = bias_k / v (small volumes under-deliver)
    dead_volume_ul: float = 15.0  # volume a tip cannot reach in a source well
    carryover: float = 0.01  # fraction of an aspirate left in a reused tip, times stickiness
    unmixed_recovery: float = 0.7  # solute fraction recovered when aspirating an unmixed well
    unmixed_sd: float = 0.08
    evap_edge_ul_per_h: float = 4.0  # open plate at 37 C, outer ring
    evap_inner_ul_per_h: float = 0.5
    lid_factor: float = 0.1  # a lid keeps this fraction of the open-plate loss
    edge_cell_loss: float = 0.15  # weaker cell response in the outer ring after overnight plating
    reader_cv: float = 0.02
    reader_floor: float = 30.0  # additive reader noise in ratio units (background)
    row_gradient: float = 0.04  # fractional signal drift from row A to row H

    def pipette_cv(self, v: float) -> float:
        return self.cv_floor + self.cv_k / max(v, 1e-9)

    def pipette_bias(self, v: float) -> float:
        return self.bias_k / max(v, 1e-9)

    @classmethod
    def randomized(cls, rng: np.random.Generator, spread: float = 1.0) -> "RigProfile":
        """Domain randomization: one plausible rig drawn around the nominal one."""
        b = cls()

        def j(x: float) -> float:
            return float(x * np.exp(rng.normal(0.0, 0.35 * spread)))

        return replace(
            b,
            name="randomized",
            cv_floor=j(b.cv_floor),
            cv_k=j(b.cv_k),
            bias_k=j(b.bias_k),
            dead_volume_ul=j(b.dead_volume_ul),
            carryover=j(b.carryover),
            unmixed_recovery=float(np.clip(b.unmixed_recovery + rng.normal(0.0, 0.1 * spread), 0.4, 0.95)),
            evap_edge_ul_per_h=j(b.evap_edge_ul_per_h),
            evap_inner_ul_per_h=j(b.evap_inner_ul_per_h),
            edge_cell_loss=min(0.5, j(b.edge_cell_loss)),
            reader_cv=j(b.reader_cv),
            reader_floor=j(b.reader_floor),
            row_gradient=j(b.row_gradient),
        )

    @classmethod
    def harsh(cls) -> "RigProfile":
        """A deliberately worse rig, held out from training, for transfer tests."""
        return replace(
            cls(),
            name="harsh",
            cv_floor=0.02,
            cv_k=0.08,
            bias_k=-0.06,
            dead_volume_ul=25.0,
            carryover=0.02,
            unmixed_recovery=0.55,
            evap_edge_ul_per_h=8.0,
            evap_inner_ul_per_h=1.0,
            edge_cell_loss=0.25,
            reader_cv=0.04,
            reader_floor=60.0,
            row_gradient=0.08,
        )

    @classmethod
    def from_qc(cls, volumes_ul, measured_cv, measured_bias=None, base: "RigProfile | None" = None, name="calibrated"):
        """Real-to-sim: fit the pipetting terms to a gravimetric or dye QC table.

        volumes_ul: tested volumes. measured_cv: CV at each volume (fraction, not percent).
        measured_bias: mean fractional error at each volume (optional).
        """
        v = np.asarray(volumes_ul, float)
        cv = np.asarray(measured_cv, float)
        design = np.column_stack([np.ones_like(v), 1.0 / v])
        floor, k = np.linalg.lstsq(design, cv, rcond=None)[0]
        out = replace(base or cls(), name=name, cv_floor=float(max(floor, 0.0)), cv_k=float(max(k, 0.0)))
        if measured_bias is not None:
            inv = 1.0 / v
            bias_k = float((inv * np.asarray(measured_bias, float)).sum() / (inv * inv).sum())
            out = replace(out, bias_k=bias_k)
        return out


class Plate:
    """A 96-well plate: volume, solute amounts (concentration times uL) and mixing state per well."""

    def __init__(self, name: str):
        self.name = name
        self.vol = np.zeros((ROWS, COLS))
        self.amt: dict[str, np.ndarray] = {}
        self.mixed = np.ones((ROWS, COLS), dtype=bool)
        self.settle = np.zeros((ROWS, COLS))

    def a(self, solute: str) -> np.ndarray:
        if solute not in self.amt:
            self.amt[solute] = np.zeros((ROWS, COLS))
        return self.amt[solute]

    def conc(self, solute: str, w) -> float:
        v = self.vol[w]
        return float(self.a(solute)[w] / v) if v > 0 else 0.0


class Rig:
    """Executes protocol steps on plates and keeps a step log that exports as a worklist."""

    def __init__(self, profile: RigProfile, rng: np.random.Generator, sticky: dict[str, float] | None = None):
        self.p = profile
        self.rng = rng
        self.sticky = sticky or {}
        self.tips = 0
        self.hours = 0.0
        self.short_aspirates = 0
        self.log: list[dict] = []
        self._residue: dict[str, float] = {}

    def _actual(self, v: float) -> float:
        if v <= 0:
            return 0.0
        return max(0.0, v * (1 + self.p.pipette_bias(v)) * (1 + self.rng.normal(0.0, self.p.pipette_cv(v))))

    def fill(self, plate: Plate, w, v: float, solutes: dict[str, float] | None = None, cv: float = 0.02) -> None:
        """Bulk reagent dispenser (cells, buffer, detection mix): its own CV, no tips, jets mix."""
        got = max(0.0, v * (1 + self.rng.normal(0.0, cv)))
        plate.vol[w] += got
        for s, c in (solutes or {}).items():
            plate.a(s)[w] += c * got
        plate.mixed[w] = True
        self.log.append({"op": "fill", "dst_plate": plate.name, "dst_well": well_name(w), "ul": v,
                         "what": "+".join(solutes) if solutes else "buffer"})

    def transfer(self, src: Plate, sw, dst: Plate, dw, v: float, new_tip: bool = True, mix_after: bool = False) -> None:
        """Tip-based transfer with volume error, dead volume, unmixed sampling and tip carryover."""
        if new_tip:
            self.tips += 1
            self._residue = {}
        want = self._actual(v)
        got = min(want, max(0.0, src.vol[sw] - self.p.dead_volume_ul))
        if got < want - 1e-9:
            self.short_aspirates += 1  # air in the tip: silent on a rig without pressure monitoring
        frac = got / src.vol[sw] if src.vol[sw] > 0 else 0.0
        rec = 1.0
        if not src.mixed[sw]:
            rec = float(np.clip(self.rng.normal(self.p.unmixed_recovery, self.p.unmixed_sd), 0.2, 1.2))
        residue: dict[str, float] = {}
        for s in list(src.amt):
            a = src.amt[s]
            moved = min(a[sw], a[sw] * frac * rec)
            a[sw] -= moved
            keep = moved * min(0.5, self.p.carryover * self.sticky.get(s, 1.0))
            residue[s] = keep
            dst.a(s)[dw] += moved - keep
        for s, r in self._residue.items():
            dst.a(s)[dw] += r
        self._residue = residue
        src.vol[sw] -= got
        dst.vol[dw] += got
        dst.mixed[dw] = bool(mix_after)
        dst.settle[dw] = 0.0
        self.log.append({"op": "transfer", "src_plate": src.name, "src_well": well_name(sw), "dst_plate": dst.name,
                         "dst_well": well_name(dw), "ul": v, "new_tip": bool(new_tip), "mix_after": bool(mix_after)})

    def incubate(self, plate: Plate, minutes: float, lid: bool = True) -> None:
        """37 C incubation: evaporation (worse on the outer ring, mostly stopped by a lid) and diffusion."""
        h = minutes / 60.0
        rate = np.where(EDGE, self.p.evap_edge_ul_per_h, self.p.evap_inner_ul_per_h)
        jitter = np.clip(1 + self.rng.normal(0.0, 0.15, size=rate.shape), 0.0, None)
        loss = rate * h * (self.p.lid_factor if lid else 1.0) * jitter
        plate.vol = np.where(plate.vol > 0, np.maximum(plate.vol - loss, 0.3 * plate.vol), plate.vol)
        plate.settle = plate.settle + minutes
        plate.mixed = plate.mixed | (plate.settle >= MIX_BY_DIFFUSION_MIN)
        self.hours += h
        self.log.append({"op": "incubate", "dst_plate": plate.name, "minutes": minutes, "lid": bool(lid)})

    def read(self, true_signal: np.ndarray) -> np.ndarray:
        """Plate reader: proportional noise, a row gradient, and an additive background floor."""
        noise = 1 + self.rng.normal(0.0, self.p.reader_cv, size=true_signal.shape)
        grad = 1 - self.p.row_gradient * (np.arange(ROWS)[:, None] / (ROWS - 1))
        floor = self.rng.normal(0.0, self.p.reader_floor, size=true_signal.shape)
        self.log.append({"op": "read"})
        return true_signal * noise * grad + floor
