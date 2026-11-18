import json
from pathlib import Path

import pytest

from drugdevgym import knowledge
from drugdevgym.policies import DR_EXPERT, POLICIES

ROOT = Path(__file__).resolve().parents[1]


def test_rules_load_and_name_only_real_traps():
    rules = knowledge.load_rules()
    assert len(rules) >= 15
    assert all(t in POLICIES[r["env"]] for r in rules for t in r["traps"])


def test_every_trap_traces_back_to_a_stated_rule():
    assert knowledge.orphan_traps(knowledge.load_rules()) == []


def test_the_expert_protocol_follows_every_checkable_rule():
    for name, check in knowledge.CHECKS.items():
        assert check(DR_EXPERT), name


def test_each_dose_response_trap_breaks_its_own_rule_and_no_other():
    rules = [r for r in knowledge.load_rules() if r["env"] == "dose_response"]
    for r in rules:
        for trap in r["traps"]:
            design = POLICIES["dose_response"][trap](0).design
            assert not knowledge.CHECKS[r["check"]](design), (r["id"], trap)
            others = [o["id"] for o in rules if o is not r and not knowledge.CHECKS[o["check"]](design)]
            assert others == [], (trap, others)


def test_bad_rules_are_rejected(tmp_path):
    good = '[[rule]]\nid = "a"\nenv = "stability"\nsaid = "s"\nwhy = "w"\ntraps = ["two_point"]\nstatus = "draft"\n'
    (tmp_path / "ok.toml").write_text(good)
    assert len(knowledge.load_rules(tmp_path / "ok.toml")) == 1
    for bad in (good.replace("two_point", "expert"), good.replace("stability", "alchemy"),
                good.replace('"draft"', '"gospel"'), good + good, good + 'moot = ["no_lid"]\n'):
        (tmp_path / "bad.toml").write_text(bad)
        with pytest.raises(ValueError):
            knowledge.load_rules(tmp_path / "bad.toml")


def test_coverage_joins_rules_to_the_committed_qc_results():
    cards = json.loads((ROOT / "results" / "envqc.json").read_text())
    rules = knowledge.load_rules()
    learned = {"a": DR_EXPERT, "b": {**DR_EXPERT, "mix": False}}
    rows = knowledge.coverage(rules, cards, learned)
    assert len(rows) == sum(len(r["traps"]) for r in rules)
    mix = next(r for r in rows if r["id"] == "mix-every-ladder-step")
    assert mix["teaches"] == "yes" and mix["learner"] == {"a": True, "b": False}
    s = knowledge.summary(rows)
    assert s["yes"] + s["on average"] + s["no"] + s["inverted"] + s["moot"] == len(rows)
    assert s["rediscovered"] == {"a": 7, "b": 6} and s["missed"] == ["mix-every-ladder-step"]
    text = knowledge.markdown(rows)
    assert "| dose_response | Mix after every transfer" in text and "Rediscovered: a 7 of 7, b 6 of 7." in text


def test_a_pair_the_expert_ruled_moot_is_kept_but_not_counted():
    cards = json.loads((ROOT / "results" / "envqc.json").read_text())
    rules = knowledge.load_rules()
    lid = next(r for r in rules if r["id"] == "lid-the-plate")
    part = next(r for r in rules if r["id"] == "cover-long-incubations")
    assert lid["status"] == "moot" and knowledge.moot_traps(lid) == {"no_lid"}
    assert part["status"] == "draft" and knowledge.moot_traps(part) == {"no_lid"}  # only that one trap of the rule
    rows = knowledge.coverage(rules, cards, {"a": DR_EXPERT})
    moot = [r for r in rows if r["moot"]]
    assert sorted((r["env"], r["trap"]) for r in moot) == [("assay_window", "no_lid"), ("dose_response", "no_lid"),
                                                            ("stability", "no_lid")]
    assert all(r["teaches"] != "not measured" and r["learner"] is None and r["note"] for r in moot)  # still measured
    s = knowledge.summary(rows)
    assert s["moot"] == 3 and s["checked"] == 7 and "lid-the-plate" not in s["missed"]
    with_them = knowledge.summary([{**r, "moot": False} for r in rows])
    assert with_them["no"] == s["no"] + sum(r["teaches"] == "no" for r in moot)
    text = knowledge.markdown(rows)
    main, tail = text.split("## Ruled moot by the expert")
    assert "no_lid |" not in main.split("## What a learner")[0] and "Lid the plate" not in main
    assert tail.count("| no_lid |") == 3 and "left out of the counts" in main


def test_a_rule_without_a_card_is_reported_as_not_measured():
    rules = knowledge.load_rules()
    rows = knowledge.coverage(rules, {})
    assert rows and all(r["teaches"] == "not measured" for r in rows)
    s = knowledge.summary(rows)
    assert s["not_measured"] + s["moot"] == len(rows) and s["moot"] == 3
    assert "not measured" in knowledge.markdown(rows)
