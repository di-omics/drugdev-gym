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
                good.replace('"draft"', '"gospel"'), good + good):
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
    assert s["yes"] + s["on average"] + s["no"] + s["inverted"] == len(rows)
    assert s["rediscovered"] == {"a": 8, "b": 7} and s["missed"] == ["mix-every-ladder-step"]
    text = knowledge.markdown(rows)
    assert "| dose_response | Mix after every transfer" in text and "Rediscovered: a 8 of 8, b 7 of 8." in text
