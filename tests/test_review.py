import copy
import json

import numpy as np
import pytest

from drugdevgym import human, make
from drugdevgym.envs import review
from drugdevgym.envs.review import PLATE_FAULTS, TIMECOURSE_FAULTS
from drugdevgym.judges import PlateJudge, TimecourseJudge
from drugdevgym.policies import POLICIES, run_episode

REVIEWS = {"plate_review": ("wells", PLATE_FAULTS), "timecourse_review": ("samples", TIMECOURSE_FAULTS)}
PROTOCOLS = {"plate_review": copy.deepcopy(review.PLATE_PROTOCOL),
             "timecourse_review": copy.deepcopy(review.TIMECOURSE_PROTOCOL)}
KINDS = {"clean": lambda t: t["fault"] == "none" and not t["should_repeat"],
         "harmful": lambda t: t["fault"] != "none" and t["should_repeat"],
         "within_tolerance": lambda t: t["fault"] != "none" and not t["should_repeat"]}


def find_seed(name, want):
    for seed in range(3000):
        if want(make(name, seed=seed).truth()):
            return seed
    raise AssertionError("no such case in 3000 seeds")


def judge(name, seed, verdict, cause):
    return make(name, seed=seed).step({"type": "judge", "verdict": verdict, "cause": cause})


@pytest.mark.parametrize("name", list(REVIEWS))
@pytest.mark.parametrize("kind", list(KINDS))
def test_all_four_reward_arms(name, kind):
    seed = find_seed(name, KINDS[kind])
    truth = make(name, seed=seed).truth()
    right, wrong = ("repeat", "accept") if truth["should_repeat"] else ("accept", "repeat")
    other = next(f for f in REVIEWS[name][1] if f != truth["fault"])
    assert judge(name, seed, right, truth["fault"])[1] == 1.0
    assert np.isclose(judge(name, seed, right, other)[1], 0.6)
    assert np.isclose(judge(name, seed, wrong, truth["fault"])[1], 0.4)
    assert judge(name, seed, wrong, other)[1] == 0.0


@pytest.mark.parametrize("name", list(REVIEWS))
def test_the_observation_is_the_same_shape_whatever_the_fault(name):
    key, faults = REVIEWS[name]
    shapes, seen = set(), set()
    for seed in range(60):
        env = make(name, seed=seed)
        obs = env.reset()
        seen.add(env.fault)
        assert set(obs) == {"phase", key, "reference", "spec"}
        shapes.add((len(obs[key]), frozenset(tuple(sorted(row)) for row in obs[key])))
        assert obs["spec"]["intended_protocol"] == PROTOCOLS[name]
        text = json.dumps({k: v for k, v in obs.items() if k != "spec"})
        assert not any(f in text for f in faults[1:]) and "should_repeat" not in text and "fault" not in text
    assert len(shapes) == 1 and seen == set(faults)


def test_plate_layout_and_nominal_doses_do_not_depend_on_the_fault():
    layouts = set()
    for seed in range(40):
        wells = make("plate_review", seed=seed).reset()["wells"]
        layouts.add(tuple((w["well"], w["role"], w["id"], w["conc_nM"], w["std_nM"]) for w in wells))
    assert len(layouts) == 1


@pytest.mark.parametrize("name", list(REVIEWS))
@pytest.mark.parametrize("rig", ["nominal", "randomized"])
def test_clean_runs_are_rarely_misleading(name, rig):
    clean = [e for e in (make(name, seed=s, rig=rig) for s in range(150)) if e.fault == "none"]
    assert len(clean) > 30 and np.mean([e.should_repeat for e in clean]) < 0.1


@pytest.mark.parametrize("name", list(REVIEWS))
def test_bad_calls_are_rejected(name):
    class Anything:
        def __eq__(self, other):
            return True

        __hash__ = None

    env = make(name, seed=0)
    for verdict, cause in (("maybe", "none"), ("accept", "gremlins"), (Anything(), "none"), ("accept", Anything()),
                           (None, "none"), ("accept", ["none"])):
        obs, reward, done, _ = make(name, seed=0).step({"type": "judge", "verdict": verdict, "cause": cause})
        assert "error" in obs and reward == 0.0 and not done
    assert "error" in env.step({"type": "guess"})[0]


@pytest.mark.parametrize("name", list(REVIEWS))
def test_editing_an_observation_cannot_change_later_cases(name):
    key, _ = REVIEWS[name]
    before = make(name, seed=11).truth()
    obs = make(name, seed=5).reset()
    obs["spec"]["intended_protocol"]["replicates"] = 1
    obs["spec"]["intended_protocol"].pop("lid")
    obs[key][0].clear()
    obs["reference"].clear()
    assert make(name, seed=11).truth() == before
    assert make(name, seed=5).reset()["spec"]["intended_protocol"] == PROTOCOLS[name]


@pytest.mark.parametrize("name", list(REVIEWS))
def test_the_executed_run_is_hidden_until_the_call_is_in(name):
    env = make(name, seed=find_seed(name, KINDS["harmful"]))
    env.reset()
    assert env.step_log == []
    env.step({"type": "judge", "verdict": "accept", "cause": "none"})
    assert len(env.step_log) > 20


def test_the_two_review_environments_do_not_share_their_fault_on_a_seed():
    both = [(make("plate_review", seed=s).fault == "none", make("timecourse_review", seed=s).fault == "none")
            for s in range(300)]
    agree = np.mean([a == b for a, b in both])
    assert 0.4 < agree < 0.65  # independent draws at a 40% clean rate agree 52% of the time


@pytest.mark.parametrize("name", list(REVIEWS))
def test_scripted_rules_beat_every_reference_policy_on_fresh_seeds(name):
    seeds = range(30000, 30080)
    mean = {p: np.mean([run_episode(make(name, seed=s), POLICIES[name][p](s))[0] for s in seeds])
            for p in POLICIES[name]}
    assert all(mean["expert"] > mean[p] + 0.15 for p in mean if p != "expert"), mean
    assert mean["null"] < 0.6 and mean["constant_repeat"] < 0.6


@pytest.mark.parametrize("name,judge_cls", [("plate_review", PlateJudge), ("timecourse_review", TimecourseJudge)])
def test_forced_verdict_traps_change_the_verdict_and_nothing_else(name, judge_cls):
    for seed in range(30):
        obs = make(name, seed=seed).reset()
        base = judge_cls().act(obs)
        for trap, verdict in (("always_accept", "accept"), ("always_repeat", "repeat")):
            call = POLICIES[name][trap](seed).act(obs)
            assert call["verdict"] == verdict and call["cause"] == base["cause"]


def test_exported_cases_carry_no_truth_but_do_carry_the_rule():
    pack = human.export_cases("plate_review", n=6, first_seed=20000)
    text = json.dumps(pack["cases"])
    assert "should_repeat" not in text and "fault" not in text and len(pack["cases"]) == 6
    assert set(pack["causes"]) == set(PLATE_FAULTS)
    assert "Z'" in pack["label"] and "0.5 log units" in pack["label"] and "0.6" in pack["scoring"]
    with pytest.raises(ValueError):
        human.export_cases("dose_response", n=1)
    with pytest.raises(ValueError):
        human.export_cases("plate_review", seeds=[])


def test_scoring_matches_the_environment_call_by_call():
    pack = human.export_cases("timecourse_review", n=8, first_seed=20000)
    truth = human.frozen_truth(pack)
    ids = [c["id"] for c in pack["cases"]]
    faulty = next(i for i in ids if truth[i]["fault"] != "none")
    calls = []
    for i in ids:
        t = truth[i]
        right, wrong = ("repeat", "accept") if t["should_repeat"] else ("accept", "repeat")
        other = next(f for f in TIMECOURSE_FAULTS if f != t["fault"])
        if i == ids[0]:
            calls.append({"id": i, "verdict": wrong, "cause": other, "realistic": False, "note": "looks\nfake"})
        elif i == ids[1]:
            calls.append({"id": i, "verdict": right, "cause": other, "note": "close call", "confidence": 2})
        else:
            calls.append({"id": i, "verdict": right, "cause": t["fault"]})
    for score in (human.score_judgments(calls), human.score_judgments(calls, packs=[pack], truth=truth)):
        rows = {r["id"]: r for r in score["rows"]}
        assert rows[ids[0]]["reward"] == 0.0 and np.isclose(rows[ids[1]]["reward"], 0.6)
        assert all(rows[i]["reward"] == 1.0 for i in ids[2:])
        for i in ids:
            seed = int(i.split("/")[2])
            expected = run_episode(make(i), POLICIES["timecourse_review"]["expert"](seed))[0]
            assert np.isclose(rows[i]["rules_reward"], expected)
        assert f"timecourse_review: {truth[faulty]['fault']}" in score["by_fault"]
        text = human.markdown(score, "tester")
        flagged = text.split("## Flagged as not realistic")[1].split("## Notes")[0]
        assert ids[0] in flagged and "looks fake" in flagged and ids[1] not in flagged
        assert "close call" in text.split("## Notes")[1]
    assert human.score_judgments(calls[:3], packs=[pack], truth=truth)["of"] == 8


def test_scoring_refuses_bad_input():
    pack = human.export_cases("plate_review", n=2, first_seed=20000)
    truth = human.frozen_truth(pack)
    ok = {"id": pack["cases"][0]["id"], "verdict": "accept", "cause": "none"}
    for bad in ([ok, dict(ok)], [{**ok, "id": "plate_review/nominal/020000"}, ok], [{**ok, "verdict": "shrug"}],
                [{**ok, "cause": "gremlins"}], [{**ok, "note": 5}], [{**ok, "realistic": "no"}],
                [{**ok, "confidence": 9}], [{**ok, "id": "dose_response/nominal/1"}], [{**ok, "id": "nonsense"}]):
        with pytest.raises(ValueError):
            human.score_judgments(bad)
    with pytest.raises(ValueError):
        human.score_judgments([{**ok, "id": "plate_review/nominal/99"}], packs=[pack], truth=truth)
    with pytest.raises(ValueError):
        human.score_judgments([ok], packs=[pack])


def test_a_balanced_human_set_covers_every_fault():
    seeds = human.balanced_seeds("timecourse_review", 20000, clean=2, per_fault=1)
    faults = [make("timecourse_review", seed=s).fault for s in seeds]
    assert sorted(faults) == sorted(["none", "none", *TIMECOURSE_FAULTS[1:]])
    assert seeds == human.balanced_seeds("timecourse_review", 20000, clean=2, per_fault=1)
    pack = human.export_cases("timecourse_review", seeds=seeds)
    assert [c["id"] for c in pack["cases"]] == [f"timecourse_review/nominal/{s}" for s in seeds]
    with pytest.raises(ValueError):
        human.balanced_seeds("timecourse_review", 20000, clean=-1, per_fault=1)
