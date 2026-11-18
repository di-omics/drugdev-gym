"""Registry: every task is (environment, rig, seed). Tasks are generated, not hand-written."""
from __future__ import annotations

from .envs import AssayWindow, DoseResponse, LeadOpt, PlateReview, Stability, TimecourseReview

ENVS = {"assay_window": AssayWindow, "dose_response": DoseResponse, "stability": Stability, "lead_opt": LeadOpt,
        "plate_review": PlateReview, "timecourse_review": TimecourseReview}
RIGS = ("nominal", "randomized", "harsh")


def make(task: str, seed: int = 0, rig="nominal"):
    """make("dose_response", seed=3) or make("dose_response/randomized/3")."""
    if "/" in task:
        name, rig, seed = task.split("/")
        seed = int(seed)
    else:
        name = task
    if name not in ENVS:
        raise KeyError(f"unknown environment: {name}; have {sorted(ENVS)}")
    return ENVS[name](seed=seed, rig=rig)


def task_id(env_name: str, rig: str, seed: int) -> str:
    return f"{env_name}/{rig}/{seed}"
