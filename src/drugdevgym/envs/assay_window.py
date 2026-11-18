"""Assay-window environment: find cell density and stimulation time that give a usable window.

Sequential: up to three pilot plates, then commit one condition. The environment runs a fresh
validation plate at the committed condition and the reward is that plate's Z' in ANALYTE
units (controls back-calculated through the on-plate standard curve, clipped to 0..1)
minus cost. Detection is a competitive immunoassay, so more cells is not better: too much
analyte pushes the signal onto the flat end of the standard curve and the window potency is
read in collapses, even while the raw-ratio Z' still looks fine.
"""
from __future__ import annotations

import numpy as np

from ..analysis import analyte_zprime, zprime
from ..chem import Kit
from ..rig import EDGE, Plate, Rig, usable_wells
from .base import Env, spec_help, validate
from .dose_response import CELLS_UL, DETECT_UL, FINAL_UL, STAMP_UL, STD_SERIES_NM

PILOT_PER_SIDE = 8
VALIDATION_PER_SIDE = 16

CONDITION_SPEC = {
    "cells_k": ("float", 1.0, 80.0, "thousand cells per well"),
    "stim_min": ("float", 15.0, 180.0, "stimulation time at 37 C in minutes"),
}
PLATE_SPEC = {
    "lid": ("bool", None, None, "lid on during stimulation"),
    "avoid_edges": ("bool", None, None, "leave the outer ring empty"),
}


class AssayWindow(Env):
    name = "assay_window"
    cost_reference = 800.0
    MAX_PILOTS = 3

    def __init__(self, seed: int = 0, rig="nominal"):
        super().__init__(seed, rig)
        r = self._truth_rng
        self.kit = Kit(k=float(r.uniform(20.0, 60.0)))
        self.per_k = float(np.exp(r.uniform(np.log(0.02), np.log(5.0))))  # analyte per 1000 cells, in kit-k units
        self.basal_frac = float(r.uniform(0.04, 0.15))
        self.tau = float(r.uniform(30.0, 90.0))
        self.cell_cv = float(r.uniform(0.04, 0.10))
        self.pilots = 0

    def spec(self) -> dict:
        return {
            "env": self.name,
            "task": ("Develop a cell assay read by a competitive TR-FRET immunoassay. Run up to "
                     f"{self.MAX_PILOTS} pilot plates of control wells, then commit one condition. "
                     "A fresh validation plate is run at the committed condition."),
            "actions": {
                "run_pilot": {"conditions": f"list of conditions, each {spec_help(CONDITION_SPEC)}; "
                                            f"{PILOT_PER_SIDE} negative and {PILOT_PER_SIDE} positive wells each; "
                                            "every plate also carries the 8-well standard curve",
                              **spec_help(PLATE_SPEC)},
                "commit": {**spec_help(CONDITION_SPEC), **spec_help(PLATE_SPEC)},
            },
            "labware": {"standard_series_nM": list(STD_SERIES_NM)},
            "reward": "Z' of the validation plate in analyte units (controls back-calculated through the "
                      "standard curve), clipped to 0..1, minus 0.15 * cost / 800",
        }

    def reset(self) -> dict:
        return {"phase": "start", "pilots_left": self.MAX_PILOTS, "spec": self.spec()}

    def _run(self, conditions: list[dict], lid: bool, avoid_edges: bool, per_side: int) -> tuple[dict, Rig]:
        rig = Rig(self.profile, self._noise)
        free = iter(usable_wells(avoid_edges))
        std_plate = Plate("assay")
        std_wells = []
        for s in STD_SERIES_NM:  # standards are added at detection, so they do not incubate
            w = next(free)
            std_wells.append((w, s))
            rig.fill(std_plate, w, CELLS_UL + STAMP_UL, {"M": s * FINAL_UL / (CELLS_UL + STAMP_UL)},
                     what=f"analyte standard {s:g} nM")
            rig.fill(std_plate, w, DETECT_UL, what="detection mix")
        conc = np.where(std_plate.vol > 0, std_plate.a("M") / np.maximum(std_plate.vol, 1e-9), 0.0)
        ratio = rig.read(self.kit.ratio(conc))
        standards = [{"std_nM": s, "ratio": round(float(ratio[w]), 1)} for w, s in std_wells]
        out = []
        for cond in conditions:
            plate = Plate("assay")
            wells = []
            for _ in range(per_side):
                for role in ("neg", "pos"):
                    w = next(free)
                    wells.append((w, role))
                    rig.fill(plate, w, CELLS_UL, what=f"cells {cond['cells_k']:g}k per well")
                    rig.fill(plate, w, STAMP_UL, {"ref": 1.0} if role == "pos" else None,
                             what="reference agonist" if role == "pos" else "vehicle")
            rig.incubate(plate, cond["stim_min"], lid)
            accumulated = 1.0 - np.exp(-cond["stim_min"] / self.tau)
            analyte = plate.a("M")
            for w, role in wells:
                edge = 1.0 - self.profile.edge_cell_loss if EDGE[w] else 1.0
                cellf = max(0.0, 1.0 + self._noise.normal(0.0, self.cell_cv)) * edge
                level = self.kit.k * self.per_k * cond["cells_k"] * accumulated
                analyte[w] = cellf * level * (1.0 if role == "pos" else self.basal_frac) * FINAL_UL
                rig.fill(plate, w, DETECT_UL, what="detection mix")
            conc = np.where(plate.vol > 0, analyte / np.maximum(plate.vol, 1e-9), 0.0)
            ratio = rig.read(self.kit.ratio(conc))
            out.append({"cells_k": cond["cells_k"], "stim_min": cond["stim_min"],
                        "neg": [round(float(ratio[w]), 1) for w, role in wells if role == "neg"],
                        "pos": [round(float(ratio[w]), 1) for w, role in wells if role == "pos"]})
            self.cost += 0.01 * cond["cells_k"] * len(wells) + 1.0 * len(wells)
        self.cost += 10.0 + 1.0 * len(std_wells)
        return {"standards": standards, "conditions": out}, rig

    @staticmethod
    def _wells_needed(n_conditions: int, per_side: int) -> int:
        return n_conditions * 2 * per_side + len(STD_SERIES_NM)

    def step(self, action: dict):
        err = self._check_open(action)
        if err:
            return self._reject(err)
        if action["type"] == "run_pilot":
            if self.pilots >= self.MAX_PILOTS:
                return self._reject("no pilot plates left; commit")
            plate, err = validate(PLATE_SPEC, {k: action.get(k) for k in PLATE_SPEC})
            if err:
                return self._reject(err)
            raw = action.get("conditions")
            if not isinstance(raw, list) or not raw:
                return self._reject("conditions must be a non-empty list")
            conditions = []
            for c in raw:
                clean, err = validate(CONDITION_SPEC, c)
                if err:
                    return self._reject(err)
                conditions.append(clean)
            need, have = self._wells_needed(len(conditions), PILOT_PER_SIDE), len(usable_wells(plate["avoid_edges"]))
            if need > have:
                return self._reject(f"{len(conditions)} conditions need {need} wells but the plate has {have}")
            data, rig = self._run(conditions, plate["lid"], plate["avoid_edges"], PILOT_PER_SIDE)
            self.pilots += 1
            self._record(rig.log)
            return ({"phase": "data", "pilot": self.pilots, **data, "cost": round(self.cost, 2),
                     "pilots_left": self.MAX_PILOTS - self.pilots}, 0.0, False, {})
        if action["type"] == "commit":
            clean, err = validate({**CONDITION_SPEC, **PLATE_SPEC}, {k: v for k, v in action.items() if k != "type"})
            if err:
                return self._reject(err)
            data, rig = self._run([clean], clean["lid"], clean["avoid_edges"], VALIDATION_PER_SIDE)
            self._record(rig.log)
            cond = data["conditions"][0]
            z = analyte_zprime(data["standards"], cond["pos"], cond["neg"])
            z = z if np.isfinite(z) else -1.0
            reward = float(np.clip(z, 0.0, 1.0)) - self._cost_penalty()
            self.done = True
            info = {"zprime_analyte": z, "zprime_raw_ratio": zprime(cond["pos"], cond["neg"]),
                    "qc_pass": bool(z >= 0.5), "cost": self.cost,
                    "truth": {"per_k": self.per_k, "basal_frac": self.basal_frac, "tau": self.tau}}
            return {"phase": "done", "validation_zprime_analyte": round(z, 3)}, reward, True, info
        return self._reject(f"unknown action type: {action['type']}")
