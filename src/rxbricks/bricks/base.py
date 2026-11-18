"""Brick: one seeded, simulated stage of drug development with a verifiable reward.

One brick instance is one episode. Observations and actions are plain JSON-style dicts so
the same brick can be driven by a script, a learner, or an LLM agent with a tool call.
The hidden truth never appears in an observation; it is returned in `info` only at the end.
"""
from __future__ import annotations

import numpy as np

from ..rig import RigProfile

COST_WEIGHT = 0.15  # reward = quality - COST_WEIGHT * cost / cost_reference


def resolve_rig(rig, rng: np.random.Generator) -> RigProfile:
    if isinstance(rig, RigProfile):
        return rig
    if rig == "nominal":
        return RigProfile()
    if rig == "randomized":
        return RigProfile.randomized(rng)
    if rig == "harsh":
        return RigProfile.harsh()
    raise ValueError(f"unknown rig: {rig!r}")


def validate(spec: dict, design) -> tuple[dict | None, str | None]:
    """Check a design against a spec of field -> (kind, lo, hi, help). Returns (clean, error)."""
    if not isinstance(design, dict):
        return None, "design must be an object"
    extra = sorted(set(design) - set(spec))
    if extra:
        return None, f"unknown field: {extra[0]}"
    clean: dict = {}
    for key, (kind, lo, hi, _help) in spec.items():
        if key not in design:
            return None, f"missing field: {key}"
        v = design[key]
        if kind == "bool":
            if not isinstance(v, (bool, np.bool_)):
                return None, f"{key} must be true or false"
            clean[key] = bool(v)
        elif kind == "times":
            if not isinstance(v, (list, tuple)) or not 2 <= len(set(v)) <= 8:
                return None, f"{key} must list 2 to 8 distinct time points"
            if any(isinstance(t, bool) or not isinstance(t, (int, float, np.integer, np.floating)) for t in v):
                return None, f"{key} must be numbers"
            if any(not lo <= t <= hi for t in v):
                return None, f"{key} must lie between {lo} and {hi}"
            clean[key] = sorted({float(t) for t in v})
        else:
            if isinstance(v, bool) or not isinstance(v, (int, float, np.integer, np.floating)):
                return None, f"{key} must be a number"
            if kind == "int" and float(v) != int(v):
                return None, f"{key} must be a whole number"
            if not lo <= v <= hi:
                return None, f"{key} must be between {lo} and {hi}"
            clean[key] = int(v) if kind == "int" else float(v)
    return clean, None


def spec_help(spec: dict) -> dict:
    return {k: {"type": kind, "min": lo, "max": hi, "help": text} for k, (kind, lo, hi, text) in spec.items()}


def log_error_score(estimate: float, truth: float, tol: float = 0.5) -> float:
    """1 at a perfect answer, falling linearly to 0 at `tol` log10 units away."""
    if estimate is None or not np.isfinite(estimate):
        return 0.0
    return float(np.clip(1.0 - abs(estimate - truth) / tol, 0.0, 1.0))


class Brick:
    name = "brick"
    cost_reference = 100.0
    MAX_INVALID = 5

    def __init__(self, seed: int = 0, rig="nominal"):
        self.seed = int(seed)
        self.rig_mode = rig if isinstance(rig, str) else rig.name
        self._truth_rng = np.random.default_rng([self.seed, 1])
        self.profile = resolve_rig(rig, np.random.default_rng([self.seed, 2]))
        self._noise = np.random.default_rng([self.seed, 3])
        self.cost = 0.0
        self.done = False
        self.last_log: list[dict] = []
        self._invalid = 0

    def spec(self) -> dict:
        raise NotImplementedError

    def reset(self) -> dict:
        raise NotImplementedError

    def step(self, action: dict) -> tuple[dict, float, bool, dict]:
        raise NotImplementedError

    def _check_open(self, action) -> str | None:
        if self.done:
            raise RuntimeError("episode is over; make a new brick")
        if not isinstance(action, dict) or "type" not in action:
            return "action must be an object with a 'type'"
        return None

    def _cost_penalty(self) -> float:
        return COST_WEIGHT * self.cost / self.cost_reference

    def _reject(self, message: str) -> tuple[dict, float, bool, dict]:
        """Invalid actions cost nothing but are counted; too many end the episode."""
        self._invalid += 1
        if self._invalid >= self.MAX_INVALID:
            self.done = True
            return {"error": message, "note": "too many invalid actions"}, -self._cost_penalty(), True, {"truncated": True}
        return {"error": message}, 0.0, False, {}
