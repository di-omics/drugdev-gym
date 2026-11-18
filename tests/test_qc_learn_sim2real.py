import csv
import io
import subprocess
from pathlib import Path

import numpy as np

from drugdevgym import envqc, learn, make
from drugdevgym.policies import DR_CHOICES, DR_EXPERT, POLICIES, run_episode
from drugdevgym.sim2real import fingerprint, gap, worklist_csv

ROOT = Path(__file__).resolve().parents[1]


def test_separation_verdicts():
    rng = np.random.default_rng(0)
    high, low = 0.9 + rng.normal(0, 0.01, 200), 0.2 + rng.normal(0, 0.01, 200)
    assert envqc.separation(high, low)["verdict"] == "clean"
    assert envqc.separation(high, high.copy())["verdict"] == "blind"
    noisy = high - 0.05 + rng.normal(0, 0.2, 200)
    assert envqc.separation(high, noisy)["zprime"] < 0.5


def test_a_large_mean_gap_is_never_called_blind():
    assert envqc.verdict(0.27, -4.0, 0.57) == "weak"  # bimodal trap: lucky half the time, ruined the rest
    assert envqc.verdict(0.03, -20.0, 0.59) == "blind"
    assert envqc.verdict(0.03, -20.0, 0.95) == "visible"
    assert envqc.verdict(-0.05, -20.0, 0.30) == "inverted"  # the reward prefers the mistake


def test_qc_card_covers_every_policy_on_the_same_seeds():
    card = envqc.qc_card("stability", n_seeds=4)
    assert [r["policy"] for r in card["rows"]] == list(POLICIES["stability"])
    assert "verdict" not in card["rows"][0] and all("verdict" in r for r in card["rows"][1:])
    assert card["first_seed"] == 0 and card["n_seeds"] == 4
    text = envqc.card_markdown(card)
    null = next(r for r in card["rows"] if r["policy"] == "null")
    assert "| expert |" in text and f"| {null['verdict']}" in text and null["verdict"] in ("clean", "visible")


def test_worklist_export_has_one_row_per_step():
    env = make("dose_response", seed=0)
    run_episode(env, POLICIES["dose_response"]["expert"](0))
    rows = list(csv.DictReader(io.StringIO(worklist_csv(env.step_log))))
    assert len(rows) == len(env.step_log) > 100
    transfers = [r for r in rows if r["op"] == "transfer"]
    assert transfers[0]["src_plate"] == "dilution" and float(transfers[0]["ul"]) == DR_EXPERT["transfer_ul"]
    assert {r["op"] for r in rows} == {"fill", "transfer", "incubate", "read"}


def _control_plate(seed, rig):
    env = make("dose_response", seed=seed, rig=rig)
    design = {**DR_EXPERT, "avoid_edges": False, "controls_per_side": 16, "n_points": 4}
    return env.step({"type": "run_plate", "design": design})[0]["wells"]


def test_fingerprint_and_gap():
    fp = fingerprint(_control_plate(0, "nominal"))
    assert set(fp) == {"neg_cv", "pos_cv", "zprime", "window", "edge_over_inner", "row_slope"}
    assert fp["zprime"] > 0.5 and fp["window"] > 1.5
    sims = [{"zprime": 0.8 + 0.01 * i, "neg_cv": 0.03} for i in range(10)]
    assert gap({"zprime": 0.84, "neg_cv": 0.03}, sims)["stats"]["zprime"]["covered"]
    off = gap({"zprime": 0.2, "neg_cv": 0.03}, sims)
    assert not off["covered"] and off["stats"]["zprime"]["z"] < -2
    assert off["stats"]["neg_cv"]["z"] == 0.0  # exact match with no simulated spread
    assert gap({"zprime": float("nan")}, sims)["covered"] is None  # nothing could be compared


def test_gap_flags_a_rig_the_simulator_was_not_calibrated_to():
    sims = [fingerprint(_control_plate(s, "nominal")) for s in range(40)]
    harsh = [gap(fingerprint(_control_plate(s, "harsh")), sims)["covered"] for s in range(100, 110)]
    nominal = [gap(fingerprint(_control_plate(s, "nominal")), sims)["covered"] for s in range(100, 110)]
    assert sum(harsh) < sum(nominal)


def test_cem_returns_a_valid_design_and_transfer_table():
    run = learn.cem(iterations=1, population=4, seeds_per_design=1)
    assert set(run["design"]) == set(DR_CHOICES) and len(run["curve"]) == 1
    table = learn.transfer_table({"expert": DR_EXPERT}, n_seeds=2)
    assert set(table[0]) == {"design", "nominal", "randomized", "harsh"}


def test_tracked_text_is_ascii_with_no_long_dashes():
    files = subprocess.run(["git", "ls-files"], cwd=ROOT, capture_output=True, text=True).stdout.split()
    assert files, "run inside the git repository"
    for name in files:
        data = (ROOT / name).read_bytes()
        assert all(b < 128 for b in data), f"non-ASCII byte in {name}"
