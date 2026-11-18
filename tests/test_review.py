import json

import numpy as np
import pytest

from drugdevgym import human, make
from drugdevgym.envs.review import PLATE_FAULTS, TIMECOURSE_FAULTS
from drugdevgym.policies import POLICIES, run_episode

REVIEWS = {"plate_review": ("wells", PLATE_FAULTS), "timecourse_review": ("samples", TIMECOURSE_FAULTS)}


@pytest.mark.parametrize("name", list(REVIEWS))
def test_a_right_call_earns_full_reward_and_a_wrong_one_none(name):
    truth = make(name, seed=3).truth()
    right = {"type": "judge", "verdict": "repeat" if truth["should_repeat"] else "accept", "cause": truth["fault"]}
    wrong_cause = next(f for f in REVIEWS[name][1] if f != truth["fault"])
    wrong = {"type": "judge", "verdict": "accept" if truth["should_repeat"] else "repeat", "cause": wrong_cause}
    assert make(name, seed=3).step(right)[1] == 1.0
    assert make(name, seed=3).step(wrong)[1] == 0.0
    half = {**right, "cause": wrong_cause}
    assert np.isclose(make(name, seed=3).step(half)[1], 0.6)


@pytest.mark.parametrize("name", list(REVIEWS))
def test_the_observation_is_the_same_shape_whatever_the_fault(name):
    key, faults = REVIEWS[name]
    shapes, seen = set(), set()
    for seed in range(60):
        env = make(name, seed=seed)
        obs = env.reset()
        seen.add(env.fault)
        assert set(obs) == {"phase", key, "reference", "spec"}
        rows = obs[key]
        shapes.add((len(rows), tuple(sorted(rows[0]))))
        assert obs["spec"]["intended_protocol"] == make(name, seed=0).spec()["intended_protocol"]
        text = json.dumps({k: v for k, v in obs.items() if k != "spec"})
        assert not any(f in text for f in faults[1:]) and "should_repeat" not in text and "fault" not in text
    assert len(shapes) == 1 and seen == set(faults)


@pytest.mark.parametrize("name", list(REVIEWS))
def test_plate_layout_and_nominal_doses_do_not_depend_on_the_fault(name):
    if name != "plate_review":
        pytest.skip("layout check is for plates")
    layouts = set()
    for seed in range(40):
        wells = make(name, seed=seed).reset()["wells"]
        layouts.add(tuple((w["well"], w["role"], w["id"], w["conc_nM"], w["std_nM"]) for w in wells))
    assert len(layouts) == 1


@pytest.mark.parametrize("name", list(REVIEWS))
def test_clean_runs_are_rarely_misleading_and_bad_calls_are_rejected(name):
    clean = [make(name, seed=s) for s in range(120)]
    clean = [e for e in clean if e.fault == "none"]
    assert np.mean([e.should_repeat for e in clean]) < 0.1
    env = make(name, seed=0)
    assert "error" in env.step({"type": "judge", "verdict": "maybe", "cause": "none"})[0]
    assert "error" in env.step({"type": "judge", "verdict": "accept", "cause": "gremlins"})[0]


@pytest.mark.parametrize("name", list(REVIEWS))
def test_scripted_rules_beat_fixed_calls_on_fresh_seeds(name):
    seeds = range(30000, 30060)
    mean = {p: np.mean([run_episode(make(name, seed=s), POLICIES[name][p](s))[0] for s in seeds])
            for p in ("expert", "always_accept", "always_repeat")}
    assert mean["expert"] > mean["always_accept"] + 0.2 and mean["expert"] > mean["always_repeat"] + 0.2


def test_exported_cases_carry_no_truth_and_scoring_matches_the_environment():
    pack = human.export_cases("plate_review", n=6, first_seed=20000)
    text = json.dumps(pack["cases"])
    assert "should_repeat" not in text and "fault" not in text and len(pack["cases"]) == 6
    assert set(pack["causes"]) == set(PLATE_FAULTS)
    calls = []
    for case in pack["cases"]:
        truth = make(case["id"]).truth()
        calls.append({"id": case["id"], "verdict": "repeat" if truth["should_repeat"] else "accept",
                      "cause": truth["fault"], "note": "looks fine" if truth["fault"] == "none" else "",
                      "realistic": truth["fault"] != "reader_drift"})
    score = human.score_judgments(calls)
    assert score["n"] == 6 and all(r["reward"] == 1.0 for r in score["rows"])
    text = human.markdown(score, "tester")
    assert "# Human review: tester" in text and "| plate_review | 6 | 1.00 | 100% | 100% |" in text
    with pytest.raises(ValueError):
        human.score_judgments([{"id": pack["cases"][0]["id"], "verdict": "shrug", "cause": "none"}])
    with pytest.raises(ValueError):
        human.export_cases("dose_response", n=1)


def test_a_balanced_human_set_covers_every_fault():
    seeds = human.balanced_seeds("timecourse_review", 20000, clean=2, per_fault=1)
    faults = [make("timecourse_review", seed=s).fault for s in seeds]
    assert sorted(faults) == sorted(["none", "none", *TIMECOURSE_FAULTS[1:]])
    assert seeds == human.balanced_seeds("timecourse_review", 20000, clean=2, per_fault=1)
    pack = human.export_cases("timecourse_review", seeds=seeds)
    assert [c["id"] for c in pack["cases"]] == [f"timecourse_review/nominal/{s}" for s in seeds]
