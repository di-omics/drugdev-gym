"""Lead-optimization campaign: the dose-response and stability bricks snapped into a DMTA loop.

A hidden structure-activity landscape covers 5^4 = 625 analogs. Each round the agent picks
up to six analogs to make and chooses the protocols used to test them; the raw plate and
LC-MS data come back exactly as the single bricks return them. After at most four rounds the
agent nominates one analog. The reward is the TRUE desirability of the nominee over the best
in the library, so sloppy measurement upstream costs reward downstream.
"""
from __future__ import annotations

import itertools

import numpy as np

from .base import Brick, spec_help
from .dose_response import DESIGN_SPEC as DR_SPEC
from .dose_response import Compound, draw_assay, plate_cost, simulate_plate
from .dose_response import check_design as check_dr
from .stability import DESIGN_SPEC as ST_SPEC
from .stability import check_design as check_st
from .stability import run_cost, simulate_stability

POSITIONS, OPTIONS = 4, 5


def desirability(pec50, log_hl):
    """Geometric mean of a potency score (1 uM to 1 nM) and a stability score (10 to 316 min)."""
    potency = np.clip((np.asarray(pec50, float) - 6.0) / 3.0, 0.0, 1.0)
    stable = np.clip((np.asarray(log_hl, float) - 1.0) / 1.5, 0.0, 1.0)
    return np.sqrt(potency * stable)


def analog_id(analog) -> str:
    return "a" + "".join(str(int(i)) for i in analog)


def parse_id(cid: str) -> tuple[int, ...]:
    return tuple(int(ch) for ch in cid[1:])


class LeadOpt(Brick):
    name = "lead_opt"
    cost_reference = 2000.0
    MAX_ROUNDS = 4
    BATCH = 6

    def __init__(self, seed: int = 0, rig="nominal"):
        super().__init__(seed, rig)
        r = self._truth_rng
        self.assay = draw_assay(r)
        zp = r.normal(0.0, 1.0, (POSITIONS, OPTIONS))
        zs = r.normal(0.0, 1.0, (POSITIONS, OPTIONS))
        self._add_p = 0.55 * zp
        self._add_s = 0.30 * (-0.4 * zp + np.sqrt(1 - 0.16) * zs)  # potency and stability trade off
        self._epi = r.normal(0.0, 0.3, (OPTIONS, OPTIONS))  # pairwise term between positions 0 and 1
        self._sticky = np.exp(r.normal(np.log(2.0), 0.8, size=(OPTIONS,) * POSITIONS))
        self.rounds = 0

    def properties(self, analog) -> tuple[float, float]:
        """True (pEC50, log10 half-life in minutes) of one analog."""
        a = tuple(int(i) for i in analog)
        pec50 = 6.3 + sum(self._add_p[i, a[i]] for i in range(POSITIONS)) + self._epi[a[0], a[1]]
        log_hl = 1.6 + sum(self._add_s[i, a[i]] for i in range(POSITIONS))
        return float(np.clip(pec50, 4.5, 10.5)), float(np.clip(log_hl, 0.7, 3.3))

    def _best(self) -> float:
        return max(float(desirability(*self.properties(a)))
                   for a in itertools.product(range(OPTIONS), repeat=POSITIONS))

    def spec(self) -> dict:
        return {
            "brick": self.name,
            "task": (f"Find the best analog in a library of {OPTIONS ** POSITIONS}: an analog is {POSITIONS} "
                     f"positions with {OPTIONS} options each, written as a list such as [0, 3, 1, 4]. Up to "
                     f"{self.MAX_ROUNDS} rounds of up to {self.BATCH} analogs, then nominate one."),
            "objective": "desirability = sqrt(clip((pEC50 - 6) / 3, 0, 1) * clip((log10(half_life_min) - 1) / 1.5, 0, 1))",
            "actions": {
                "make_test": {"analogs": "list of analogs", "dose_response": spec_help(DR_SPEC),
                              "stability": spec_help(ST_SPEC)},
                "nominate": {"analog": "the analog to advance"},
            },
            "reward": "true desirability of the nominee / best desirability in the library, minus 0.15 * cost / 2000",
        }

    def reset(self) -> dict:
        return {"phase": "start", "rounds_left": self.MAX_ROUNDS, "spec": self.spec()}

    @staticmethod
    def _clean_analog(a):
        ok = (isinstance(a, (list, tuple)) and len(a) == POSITIONS
              and all(isinstance(i, (int, np.integer)) and not isinstance(i, bool) and 0 <= i < OPTIONS for i in a))
        return tuple(int(i) for i in a) if ok else None

    def step(self, action: dict):
        err = self._check_open(action)
        if err:
            return self._reject(err)
        if action["type"] == "make_test":
            if self.rounds >= self.MAX_ROUNDS:
                return self._reject("no rounds left; nominate")
            raw = action.get("analogs")
            if not isinstance(raw, list) or not 1 <= len(raw) <= self.BATCH:
                return self._reject(f"analogs must list 1 to {self.BATCH} analogs")
            analogs = [self._clean_analog(a) for a in raw]
            if any(a is None for a in analogs) or len(set(analogs)) != len(analogs):
                return self._reject(f"each analog must be {POSITIONS} distinct-per-batch integers in 0..{OPTIONS - 1}")
            dr, err = check_dr(action.get("dose_response"), len(analogs))
            if err:
                return self._reject("dose_response: " + err)
            st, err = check_st(action.get("stability"), len(analogs))
            if err:
                return self._reject("stability: " + err)
            compounds, half_lives = [], {}
            for a in analogs:
                pec50, log_hl = self.properties(a)
                compounds.append(Compound(analog_id(a), pec50, 0.9, 1.0, float(self._sticky[a])))
                half_lives[analog_id(a)] = 10.0 ** log_hl
            wells, rig_dr = simulate_plate(self.profile, self._noise, self.assay, compounds, dr)
            samples, rig_st = simulate_stability(self.profile, self._noise, half_lives, st)
            self.last_log = rig_dr.log + rig_st.log
            self.cost += 20.0 * len(analogs)  # synthesis
            self.cost += plate_cost(len(wells), rig_dr.tips) + run_cost(len(samples), rig_st.tips)
            self.rounds += 1
            return ({"phase": "data", "round": self.rounds, "dose_response": {"wells": wells},
                     "stability": {"samples": samples}, "cost": round(self.cost, 2),
                     "rounds_left": self.MAX_ROUNDS - self.rounds}, 0.0, False, {})
        if action["type"] == "nominate":
            analog = self._clean_analog(action.get("analog"))
            if analog is None:
                return self._reject("analog must be a list of position choices")
            pec50, log_hl = self.properties(analog)
            score = float(desirability(pec50, log_hl))
            best = self._best()
            quality = score / best if best > 0 else 0.0
            self.done = True
            info = {"quality": quality, "desirability": score, "best": best, "cost": self.cost,
                    "truth": {"pEC50": round(pec50, 3), "half_life_min": round(10.0 ** log_hl, 1)}}
            return {"phase": "done"}, quality - self._cost_penalty(), True, info
        return self._reject(f"unknown action type: {action['type']}")
