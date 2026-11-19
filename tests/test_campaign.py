import json
from pathlib import Path

import numpy as np
import pytest

from drugdevgym import campaign, cli, knowledge
from drugdevgym.policies import DR_EXPERT, POLICIES, DosePolicy
from drugdevgym.registry import ENVS, make

ROOT = Path(__file__).resolve().parents[1]
CARDS = ROOT / "results" / "envqc.json"
CUT = "## Every pair"  # below this line the report says what each seed's data did; above it, what the designs are


def _pairs():
    return [(r["env"], t) for r in knowledge.load_rules() for t in r["traps"]]


def test_every_rule_trap_pair_has_a_kind_and_nothing_else_does():
    stated = {(env, trap) for env, traps in campaign.TRAP_KIND.items() for trap in traps}
    assert stated == set(_pairs())
    assert all(kind in campaign.KINDS for traps in campaign.TRAP_KIND.values() for kind in traps.values())
    assert set(campaign.PRICED_BY) == set(campaign.KINDS) == set(campaign.GROUP_OF)
    assert sorted(campaign.KINDS) == sorted(k for kinds in campaign.GROUPS.values() for k in kinds)
    assert set(campaign.GROUP_TITLE) == set(campaign.GROUPS)
    # the notes name real pairs only, and say something
    assert set(campaign.NEEDS) <= stated and set(campaign.ALSO) <= stated
    assert all(note.strip() for note in campaign.NEEDS.values())
    upstream = {pair for pair in stated if campaign.TRAP_KIND[pair[0]][pair[1]] == "upstream"}
    assert set(campaign.UPSTREAM) == upstream and set(campaign.UPSTREAM.values()) <= stated
    # every pair priced on the rig names its direct measurement (the pairs the expert ruled moot need none)
    moot = {(r["env"], t) for r in knowledge.load_rules() for t in knowledge.moot_traps(r)}
    assert all(pair in campaign.NEEDS for pair in stated - moot
               if campaign.TRAP_KIND[pair[0]][pair[1]] in campaign.GROUPS["rig"])


def test_every_design_field_an_agent_can_set_has_a_kind():
    field_kinds = {"dose", "data", "qc", "condition", "handling", "position", "strategy"}
    assert set(campaign.FIELD_KIND.values()) == field_kinds and set(campaign.RUN_KIND.values()) <= field_kinds
    fields_by_env = campaign.design_fields()
    assert fields_by_env["plate_review"] == set()  # a reviewer sets nothing on the rig
    assert {"mix", "standards", "controls_per_side"} <= fields_by_env["dose_response"] <= set(campaign.FIELD_KIND)
    for name, fields in fields_by_env.items():
        for leaf in fields:
            assert campaign.FIELD_KIND.get(leaf) in field_kinds, (name, leaf)
    for name in ENVS:
        for action_type in make(name, seed=0).spec()["actions"]:
            assert action_type in campaign.TERMINAL or action_type in campaign.RUN_KIND, (name, action_type)


def test_each_stated_kind_agrees_with_what_the_simulator_shows_the_trap_changes():
    traces: dict = {}
    wrong = {}
    for env, trap in _pairs():
        message = campaign.check(env, trap, traces=traces)
        if message:
            wrong[(env, trap)] = message
    assert wrong == {}


@pytest.mark.parametrize("env, trap, kind", [
    ("dose_response", "no_mix", "reading"),  # its runs differ
    ("dose_response", "no_mix", "data"),  # it changes a handling field, not a data field
    ("stability", "no_lloq", "handling"),  # it runs exactly what the expert runs
    ("stability", "no_lloq", "judgment"),  # and the environment does run something
    ("plate_review", "always_accept", "reading"),  # nothing is run there at all
    ("dose_response", "narrow_range", "handling"),  # the dose fields it changes would be unexplained
    ("stability", "no_lid_on_edge_wells", "position"),  # the lid would be unexplained
    ("lead_opt", "upstream_no_mix", "strategy"),  # the assay field it changes would be unexplained
    ("lead_opt", "upstream_no_mix", "handling"),  # the analogs would be unexplained
    ("lead_opt", "no_model", "upstream"),  # it changes no assay field: that is strategy
])
def test_a_wrong_kind_is_caught(monkeypatch, env, trap, kind):
    if kind == "upstream":
        monkeypatch.setitem(campaign.UPSTREAM, (env, trap), ("dose_response", "no_mix"))
    monkeypatch.setitem(campaign.TRAP_KIND[env], trap, kind)
    assert campaign.check(env, trap), (env, trap, kind)


@pytest.mark.parametrize("source, message", [
    (("stability", "two_point"), "outside that assay's design"),  # another assay's field of another name
    (("gi_camp_potency", "no_mix"), "outside that assay's design"),  # the same field name in an assay it does not nest
    (("dose_response", "reuse_tips"), "and that trap changes new_tips"),  # the right assay, another field
    (("stability", "no_lloq"), "which is reading"),  # a trap that changes no run
    (("lead_opt", "upstream_raw_ratio"), "which is upstream"),
])
def test_an_upstream_trap_has_to_name_the_assay_trap_it_repeats(monkeypatch, source, message):
    monkeypatch.setitem(campaign.UPSTREAM, ("lead_opt", "upstream_no_mix"), source)
    assert message in campaign.check("lead_opt", "upstream_no_mix")


def test_an_upstream_trap_has_to_set_the_same_value(monkeypatch):
    # early_only changes the same field as two_point, to other time points
    monkeypatch.setitem(campaign.UPSTREAM, ("lead_opt", "upstream_two_point"), ("stability", "early_only"))
    assert "to another value" in campaign.check("lead_opt", "upstream_two_point")
    monkeypatch.delitem(campaign.UPSTREAM, ("lead_opt", "upstream_two_point"))
    assert "no assay trap is named" in campaign.check("lead_opt", "upstream_two_point")


def test_a_wrong_note_is_caught(monkeypatch):
    assert campaign.check("stability", "not_a_trap") == "no kind stated"
    for pair in (("dose_response", "no_mix"), ("stability", "no_lloq"), ("plate_review", "always_accept"),
                 ("lead_opt", "upstream_no_mix")):  # a note that explains nothing, under every kind of rule
        monkeypatch.setitem(campaign.ALSO, pair, {"lid": "a stale note"})
        assert "ALSO explains lid" in campaign.check(*pair), pair
        monkeypatch.delitem(campaign.ALSO, pair)
    monkeypatch.setitem(campaign.ALSO, ("dose_response", "narrow_range"), {"transfer_ul": " "})
    assert "ALSO gives no reason for transfer_ul" in campaign.check("dose_response", "narrow_range")
    monkeypatch.delitem(campaign.ALSO, ("dose_response", "narrow_range"))
    assert "also changes transfer_ul" in campaign.check("dose_response", "narrow_range")
    with pytest.raises(KeyError):
        campaign.field_kind("a_field_nobody_classified")
    with pytest.raises(KeyError):
        campaign.field_kind("-a_run_nobody_classified")
    with pytest.raises(ValueError):
        campaign.compare("dose_response", "no_mix", seeds=[])


def test_a_trap_that_follows_the_data_can_equal_the_expert_on_some_seeds():
    seen = campaign.compare("gi_camp_potency", "too_many_cells", seeds=range(6))
    assert seen["seeds"] == 6 and 0 < seen["differs_on"] < 6 and seen["fields"] == ["cells_k"]
    assert campaign.compare("stability", "no_lloq")["differs_on"] == 0
    # on seed 1 alone the expert itself picks the trap's cell number, and the message says so
    assert "check more seeds" in campaign.check("gi_camp_potency", "too_many_cells", seeds=[1])
    assert campaign.check("gi_camp_potency", "too_many_cells", seeds=iter(range(3))) is None  # any iterable of seeds


def _plate(**design):
    return json.dumps({"type": "run_plate", "design": {k: campaign._canon(v) for k, v in design.items()}}, sort_keys=True)


def test_runs_are_compared_as_sets_of_distinct_designs():
    plate, unmixed = _plate(mix=True, n_points=10), _plate(mix=False, n_points=10)
    pilot = json.dumps({"type": "run_titration", "design": {"cells_k": 10.0}}, sort_keys=True)
    assert campaign._differences([plate], [plate]) == set()
    assert campaign._differences([plate], [unmixed]) == {"mix"}
    assert campaign._differences(sorted([plate, pilot]), [plate]) == {"-run_titration"}
    assert campaign._differences([plate], sorted([plate, pilot])) == {"+run_titration"}
    nested = json.dumps({"type": "make_test", "analogs": [[0.0]], "stability": {"lid": True}}, sort_keys=True)
    other = json.dumps({"type": "make_test", "analogs": [[0.0]], "stability": {"lid": False}}, sort_keys=True)
    assert campaign._differences([nested], [other]) == {"stability.lid"}
    assert campaign.field_kind("stability.lid") == "handling" and campaign.field_kind("-run_titration") == "condition"


def test_runs_are_paired_so_that_the_fewest_fields_differ():
    # pairing in sort order would also report n_points; the best pairing does not
    expert = sorted([_plate(mix=False, n_points=10, replicates=1), _plate(mix=False, n_points=8, replicates=2)])
    trap = sorted([_plate(mix=False, n_points=8, replicates=1), _plate(mix=True, n_points=10, replicates=1)])
    assert campaign._differences(expert, trap) == {"mix", "replicates"}
    # every assignment is tried, so more unmatched runs than MAX_PAIRED on one side is refused, not guessed
    n = campaign.MAX_PAIRED
    assert campaign._differences([_plate(n_points=i) for i in range(n)], [_plate(n_points=i, mix=True) for i in range(n)]) == {"mix"}
    with pytest.raises(ValueError, match="too many to pair"):
        campaign._differences([_plate(n_points=i) for i in range(n + 1)], [_plate(n_points=i, mix=True) for i in range(n + 1)])


class _Twice:
    """Runs one plate, runs the same plate again, then reports."""

    def __init__(self):
        self.plates = 0

    def act(self, obs):
        if obs.get("phase") == "data" and self.plates == 2:
            return {"type": "report", "pEC50": {}}
        self.plates += 1
        return {"type": "run_plate", "design": DR_EXPERT}


def test_trace_counts_what_an_episode_ran_and_measured(monkeypatch):
    dose = campaign.trace("dose_response", "expert", 0)
    assert len(dose["ran"]) == 1 and dose["wells"] == 60 and dose["samples"] == 0 and dose["rig_runs"] == 1
    time_course = campaign.trace("stability", "expert", 0)
    assert time_course["wells"] == 0 and time_course["samples"] == 48
    review = campaign.trace("plate_review", "expert", 0)  # a reviewer runs nothing; the plate shown is not theirs
    assert review == {"ran": [], "rig_runs": 0, "wells": 0, "samples": 0}
    assert campaign.changed("stability", "no_lloq") == [] and campaign.changed("dose_response", "no_mix") == ["mix"]
    # a repeated plate is one design, measured once, and counted as run twice
    monkeypatch.setitem(POLICIES["dose_response"], "twice", lambda seed: _Twice())
    assert campaign.trace("dose_response", "twice", 0) == {**dose, "rig_runs": 2}
    # an action the environment rejects is not a run
    monkeypatch.setitem(POLICIES["dose_response"], "rejected", lambda seed: DosePolicy({**DR_EXPERT, "n_points": 99}))
    assert campaign.trace("dose_response", "rejected", 0) == {"ran": [], "rig_runs": 0, "wells": 0, "samples": 0}
    # a real trap that repeats one design: every plate of on_target_only is the same on-target plate
    repeat = campaign.trace("selectivity_panel", "on_target_only", 0)
    assert len(repeat["ran"]) == 1 and repeat["wells"] == 60 and repeat["rig_runs"] == 5


def test_canon_compares_numbers_by_value():
    assert campaign._canon({"a": np.int64(6), "b": (1, 2.0), "c": np.bool_(True), "d": None, "e": np.array([1, 2])}) == \
        {"a": 6.0, "b": [1.0, 2.0], "c": True, "d": None, "e": [1.0, 2.0]}


def test_pairs_per_arm_is_the_paired_arithmetic():
    assert campaign.pairs_per_arm() == campaign.pairs_per_arm(0.15, 0.25) == 8
    assert campaign.pairs_per_arm(0.1, 0.25) == 5
    assert campaign.pairs_per_arm(0.05, 0.25) == 3
    assert campaign.pairs_per_arm(0.2, 0.25) == 12  # the approximation's value; the exact paired test asks for 13
    assert campaign.pairs_per_arm(0.01, 1.0) == 3  # never fewer than three pairs
    assert campaign.pairs_per_arm(0.1, 0.25, power=0.9) > campaign.pairs_per_arm(0.1, 0.25)
    for bad in ((0.0, 0.25), (0.1, 0.0), (float("inf"), 0.25), (0.1, float("inf")), (float("nan"), 0.25), (True, 0.25),
                (0.1, 1e-200), (1e200, 0.25), (0.1, 0.25, 1.5), (0.1, 0.25, 0.0), (0.1, 0.25, 0.05, 0.0),
                (0.1, 0.25, 0.05, 1.0), (0.1, 0.25, 0.05, 0.01)):
        with pytest.raises(ValueError):
            campaign.pairs_per_arm(*bad)


@pytest.mark.parametrize("sd, pairs", [(0.15, 8), (0.1, 5)])
def test_the_paired_arithmetic_reaches_its_power_in_simulation(sd, pairs):
    assert campaign.pairs_per_arm(sd, 0.25) == pairs
    rng = np.random.default_rng(0)

    def power(n):
        null = rng.standard_normal((200000, n))
        critical = np.quantile(np.abs(null.mean(1) / (null.std(1, ddof=1) / np.sqrt(n))), 0.95)
        d = rng.normal(0.25, sd * np.sqrt(2), (200000, n))  # n pairs, one run scatters by sd, shift 0.25
        return np.mean(np.abs(d.mean(1) / (d.std(1, ddof=1) / np.sqrt(n))) > critical)

    assert power(pairs) > 0.8 > power(pairs - 1)


def test_the_plan_adds_up_and_matches_the_committed_report():
    p = campaign.plan(json.loads(CARDS.read_text()))
    c, t = p["counts"], p["totals"]
    assert sum(c[g] for g in campaign.GROUPS) == c["pairs"] == sum(cell["pairs"] for cell in p["by_kind"].values())
    assert c["pairs"] + c["moot"] == len(_pairs()) == len(p["rows"])
    assert c["records_plain"] + c["records_conditional"] == c["records"]
    assert 0 < c["dose_top_well"] < p["by_kind"]["dose"]["pairs"]
    assert all(sum(cell.values()) == 2 * cell["pairs"] for cell in p["by_kind"].values())  # pairs plus its verdict cells
    arms = [a for a in p["arms"] if a["kind"] != "reference"]
    assert len(arms) == t["arms"] == p["by_kind"]["handling"]["pairs"] and all(a["kind"] == "handling" for a in arms)
    assert t["reference_arms"] == len({a["env"] for a in p["arms"]})  # one reference arm per environment with an arm
    assert t["pairs_per_arm"] == campaign.pairs_per_arm() and t["runs"] == (t["arms"] + t["reference_arms"]) * t["pairs_per_arm"]
    assert t["wells"] == round(sum(a["wells"] for a in p["arms"]) * t["pairs_per_arm"]) > 0 and t["samples"] > 0
    assert all(r["wells"] == r["samples"] == r["rig_runs"] == 0 for r in p["rows"] if r["kind"] == "judgment")
    assert all(r["rig_runs"] >= 1 for r in p["rows"] if r["kind"] != "judgment")
    assert "replicates" in p["unexercised"]["dose_response"] and "lead_opt" not in p["unexercised"]
    assert set(p["unexercised_runs"]) < set(campaign.RUN_KIND)
    # the committed report: the text above the cut exactly; below it, every row
    text, committed = campaign.markdown(p), (ROOT / "results" / "campaign.md").read_text()
    assert CUT in text and text.split(CUT)[0] == committed.split(CUT)[0]
    stored = json.loads((ROOT / "results" / "campaign.json").read_text())
    fresh = json.loads(json.dumps(p))
    for key in ("counts", "by_kind", "by_group", "by_env", "totals", "arms", "times", "grids", "assumptions", "seeds"):
        assert stored[key] == fresh[key], key
    follows_data = ("differs_on", "changes")  # what a trap that follows the data did on these seeds
    for old, new in zip(stored["rows"], fresh["rows"], strict=True):
        skip = follows_data if new["kind"] == "condition" else ()
        assert {k: v for k, v in old.items() if k not in skip} == {k: v for k, v in new.items() if k not in skip}
    for heading in ("## Fields a trap changes outside its own kind", "## Kinds nothing checks yet", "## Ruled moot by the expert"):
        assert heading in text and heading in committed
    assert text.split(CUT)[1].count("\n| ") == committed.split(CUT)[1].count("\n| ")


def test_the_plan_refuses_a_missing_or_a_wrong_kind(monkeypatch):
    cards = json.loads(CARDS.read_text())
    monkeypatch.setitem(campaign.TRAP_KIND["dose_response"], "no_mix", "reading")
    with pytest.raises(ValueError, match="dose_response/no_mix: called reading"):
        campaign.plan(cards)
    monkeypatch.delitem(campaign.TRAP_KIND["dose_response"], "no_mix")
    with pytest.raises(KeyError, match="no TRAP_KIND for dose_response/no_mix"):
        campaign.plan(cards)


def test_a_pair_with_no_card_is_not_counted_as_blind():
    p = campaign.plan({})
    assert all(cell["not measured"] == cell["pairs"] and cell["no"] == 0 for cell in p["by_kind"].values())
    assert "with no QC card yet" in campaign.markdown(p)


def test_an_arm_that_does_not_report_a_log_value_says_so(monkeypatch):
    monkeypatch.setitem(campaign.NOT_A_LOG_VALUE, "dose_response", "a Z-prime")
    text = campaign.markdown(campaign.plan(json.loads(CARDS.read_text())))
    assert "dose_response reports a Z-prime, not a log value" in text


def test_the_command_writes_the_committed_report(tmp_path, capsys):
    cli.main(["campaign", "--cards", str(CARDS), "--out", str(tmp_path / "out")])
    capsys.readouterr()
    assert (tmp_path / "out" / "campaign.md").read_text().split(CUT)[0] == \
        (ROOT / "results" / "campaign.md").read_text().split(CUT)[0]
    assert json.loads((tmp_path / "out" / "campaign.json").read_text())["totals"] == \
        json.loads((ROOT / "results" / "campaign.json").read_text())["totals"]


def test_the_command_passes_its_numbers_on(tmp_path, capsys):
    cli.main(["campaign", "--cards", str(CARDS), "--out", str(tmp_path / "other"), "--seeds", "2", "--sd", "0.1",
              "--delta", "0.3", "--alpha", "0.1", "--power", "0.9"])
    capsys.readouterr()
    stored = json.loads((tmp_path / "other" / "campaign.json").read_text())
    assert stored["assumptions"] == {"sd_within": 0.1, "delta": 0.3, "alpha": 0.1, "power": 0.9}
    assert stored["seeds"] == [0, 1] and stored["totals"]["pairs_per_arm"] == campaign.pairs_per_arm(0.1, 0.3, 0.1, 0.9)
    assert "seeds 0, 1" in (tmp_path / "other" / "campaign.md").read_text()


def test_the_command_rejects_bad_input_and_keeps_the_committed_report(tmp_path, monkeypatch):
    (tmp_path / "list.json").write_text("[1, 2]")
    (tmp_path / "broken.json").write_text("{")
    (tmp_path / "a_file").write_text("x")
    for bad in (["--sd", "0"], ["--delta", "inf"], ["--delta", "1e-200"], ["--power", "1"], ["--alpha", "0"],
                ["--cards", str(tmp_path / "nothing.json")], ["--cards", str(tmp_path)],
                ["--cards", str(tmp_path / "list.json")], ["--cards", str(tmp_path / "broken.json")]):
        args = ["campaign", "--out", str(tmp_path / "bad")] + (bad if "--cards" in bad else ["--cards", str(CARDS)] + bad)
        with pytest.raises(SystemExit, match="campaign: "):
            cli.main(args)
    assert not (tmp_path / "bad").exists()
    with pytest.raises(SystemExit, match="is not a directory"):
        cli.main(["campaign", "--cards", str(CARDS), "--out", str(tmp_path / "a_file")])
    # with a number changed, the default output directory is refused: results/ holds the committed report
    monkeypatch.chdir(tmp_path)
    for other in (["--sd", "0.2"], ["--seeds", "5"], ["--power", "0.9"]):
        with pytest.raises(SystemExit, match="give --out"):
            cli.main(["campaign", "--cards", str(CARDS)] + other)
    assert not (tmp_path / "results").exists()
