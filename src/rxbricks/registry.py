"""Registry: every task is (brick, rig, seed). Tasks are generated, not hand-written."""
from __future__ import annotations

from .bricks import AssayWindow, DoseResponse, LeadOpt, Stability

BRICKS = {"assay_window": AssayWindow, "dose_response": DoseResponse, "stability": Stability, "lead_opt": LeadOpt}
RIGS = ("nominal", "randomized", "harsh")


def make(task: str, seed: int = 0, rig="nominal"):
    """make("dose_response", seed=3) or make("dose_response/randomized/3")."""
    if "/" in task:
        name, rig, seed = task.split("/")
        seed = int(seed)
    else:
        name = task
    if name not in BRICKS:
        raise KeyError(f"unknown brick: {name}; have {sorted(BRICKS)}")
    return BRICKS[name](seed=seed, rig=rig)


def task_id(brick: str, rig: str, seed: int) -> str:
    return f"{brick}/{rig}/{seed}"
