import json
from types import SimpleNamespace

import pytest

from drugdevgym import agent, baseline, campaign, cli, knowledge
from drugdevgym.agent import Reply
from drugdevgym.policies import POLICIES, run_episode
from drugdevgym.registry import ENVS, make
from drugdevgym.rollouts import _plain
from drugdevgym.tools import TOOLS, call_tool

SEED = baseline.FIRST_SEED


def _expert_data(env_name: str, seed: int = SEED):
    """Run the scripted expert up to its final action. Returns (observations seen, final action)."""
    env, policy = make(env_name, seed=seed), POLICIES[env_name]["expert"](seed)
    obs = env.reset()
    obs.pop("spec")
    seen = [obs]
    while True:
        action = policy.act(obs)
        obs, _reward, done, _info = env.step(action)
        if done:
            return seen, action
        seen.append(obs)


class ToolUser:
    """A fake model: runs a scripted policy's experiments, asks one tool for the numbers, reports them."""
    name = "fake:tool-user"

    def __init__(self, env_name, seed, tool, args, key, policy="expert"):
        self.policy, self.tool, self.args, self.key = POLICIES[env_name][policy](seed), tool, args, key
        self.obs, self.asked = {}, False

    def observe(self, obs):
        self.obs = obs

    def __call__(self, system, messages, attempt):
        last = json.loads(messages[-1]["content"])
        if "tool_result" in last:
            return Reply(json.dumps({"action": {"type": "report", self.key: last["tool_result"]}}))
        action = self.policy.act(self.obs)
        if action["type"] == "report" and not self.asked:
            self.asked = True
            return Reply("I will fit the data first.\n" + json.dumps({"tool": self.tool, "args": self.args}))
        return Reply(json.dumps({"action": action}, default=agent._default))


@pytest.mark.parametrize("env_name", list(ENVS))
def test_a_scripted_policy_scores_the_same_through_the_protocol(env_name):
    for policy in POLICIES[env_name]:
        want, info = run_episode(make(env_name, seed=SEED), POLICIES[env_name][policy](SEED))
        got = agent.run_episode(env_name, SEED, agent.ScriptedModel(env_name, policy, SEED))
        assert got["reward"] == want and got["finished"] and got["stop"] == "done", (env_name, policy)
        assert got["model"] == f"scripted:{policy}" and got["task"] == f"{env_name}/nominal/{SEED}"
        assert json.loads(json.dumps(got)) == got  # the whole record is plain JSON


@pytest.mark.parametrize("env_name", list(ENVS))
def test_the_model_is_shown_the_task_and_the_data_and_never_the_truth(env_name):
    record = agent.run_episode(env_name, SEED, agent.ScriptedModel(env_name, "expert", SEED))
    assert "truth" in record["info"]  # kept in the record, for the profile
    for message in record["transcript"]:
        if message["role"] != "user":
            continue
        said = json.loads(message["content"])
        assert set(said) <= {"task", "tools", "observation", "tool_result", "error"}
        assert "truth" not in json.dumps(said.get("observation")) and "reward" not in said
    first = json.loads(record["transcript"][0]["content"])
    assert first["task"]["env"] == env_name and [t["name"] for t in first["tools"]] == [t.name for t in TOOLS[env_name]]


@pytest.mark.parametrize("env_name, tool, args, key", [
    ("dose_response", "estimate_pec50", {"use_standards": True}, "pEC50"),
    ("stability", "estimate_half_life", {"use_lloq": True}, "half_life_min"),
    ("gi_fluid_stability", "estimate_half_life", {"correct_decay": True, "correct_activity": True, "use_lloq": True}, "half_life_min"),
    ("gi_camp_potency", "estimate_pec50", {"use_standards": True, "anchor": True}, "pEC50"),
    ("cytotox_window", "estimate_cc50", {"use_turbidity": True, "bound_when_clean": True}, "results"),
    ("selectivity_panel", "selectivity_windows", {"use_standards": True}, "windows"),
])
def test_a_tool_gives_the_numbers_the_scripted_expert_reports(env_name, tool, args, key):
    seen, final = _expert_data(env_name)
    result, ok = call_tool(TOOLS[env_name], tool, args, seen)
    assert ok, result
    assert (result["windows"] if key == "windows" else result) == _plain(final[key])


def test_every_tool_runs_on_the_data_of_its_environment():
    for env_name, toolset in TOOLS.items():
        seen, _final = _expert_data(env_name)
        for tool in toolset:
            args = {k: (True if kind == "bool" else kind[len("one of "):].split(", ")[0]) for k, (kind, _text) in tool.args.items()}
            result, ok = call_tool(toolset, tool.name, args, seen)
            assert ok, (env_name, tool.name, result)
            json.dumps(result)
            described = tool.describe()
            assert described["args"]["run"]["required"] is False and all(a["help"] for a in described["args"].values())
    assert set(TOOLS) == set(ENVS)


def test_a_tool_call_is_checked_and_never_raises():
    seen, _final = _expert_data("dose_response")
    tools = TOOLS["dose_response"]
    assert call_tool(tools, "estimate_pec50", {}, seen) == ("missing argument: use_standards", False)
    assert call_tool(tools, "estimate_pec50", {"use_standards": "yes"}, seen) == ("use_standards must be true or false", False)
    assert call_tool(tools, "estimate_pec50", {"use_standards": True, "plate": 1}, seen) == ("unknown argument: plate", False)
    assert call_tool(tools, "estimate_pec50", "x", seen) == ("args must be an object", False)
    assert call_tool(tools, "fit_it", {}, seen)[0].startswith("no such tool: fit_it")
    assert call_tool(tools, "control_zprime", {"run": 7}, seen) == ("control_zprime failed: run must be a whole number from 1 to 1", False)
    assert call_tool(tools, "control_zprime", {}, seen[:1]) == ("control_zprime failed: no observation with wells yet", False)
    assert call_tool(TOOLS["assay_window"], "condition_zprime", {"metric": "best"}, []) == ("metric must be one of analyte, raw", False)
    first, ok = call_tool(tools, "control_zprime", {"run": 1}, seen)
    assert ok and first == call_tool(tools, "control_zprime", None, seen)[0]


def test_the_last_json_object_with_an_action_or_a_tool_is_the_reply():
    assert agent.extract_json('I think {"x": 1} then {"action": {"type": "report", "pEC50": {"a": 7.1}}}') == \
        {"action": {"type": "report", "pEC50": {"a": 7.1}}}
    assert agent.extract_json('{"tool": "a", "args": {}} no, rather {"tool": "b", "args": {"run": 1}}\n```') == \
        {"tool": "b", "args": {"run": 1}}
    assert agent.extract_json('{"action": {"type": "run", "design": {"lid": true}}} trailing words') is not None
    for nothing in ("", "no json here", '{"type": "report"}', "{broken", '["action"]'):
        assert agent.extract_json(nothing) is None


def test_a_model_that_uses_a_tool_is_scored_and_its_reading_is_seen():
    lloq = dict(env_name="stability", seed=SEED, tool="estimate_half_life", key="half_life_min")
    want, _info = run_episode(make("stability", seed=SEED), POLICIES["stability"]["expert"](SEED))
    good = agent.run_episode("stability", SEED, ToolUser(args={"use_lloq": True}, **lloq))
    assert good["reward"] == want and good["turns"] == 3 and good["tool_calls"] == [
        {"tool": "estimate_half_life", "args": {"use_lloq": True}, "ok": True}]
    assert baseline.mistakes(good)["no_lloq"] == {"made": False, "how": "tool flag"}
    bad = agent.run_episode("stability", SEED, ToolUser(args={"use_lloq": False}, **lloq))
    assert baseline.mistakes(bad)["no_lloq"] == {"made": True, "how": "tool flag"}
    free = agent.run_episode("gi_camp_potency", SEED, ToolUser("gi_camp_potency", SEED, "estimate_pec50",
                                                             {"use_standards": True, "anchor": False}, "pEC50"))
    found = baseline.mistakes(free)
    assert found["free_plateau"]["made"] is True and found["raw_ratio"]["made"] is False
    # a model that never calls the tool cannot be told
    silent = agent.run_episode("stability", SEED, agent.ScriptedModel("stability", "expert", SEED))
    assert baseline.mistakes(silent)["no_lloq"] == {"made": None, "how": "tool flag"}


def test_an_episode_that_does_not_finish_scores_nothing():
    def run(reply, **kwargs):
        return agent.run_episode("dose_response", SEED, lambda system, messages, attempt: reply, **kwargs)
    lost = run(Reply("I am not sure what to do."))
    assert (lost["stop"], lost["finished"], lost["reward"], lost["invalid_replies"]) == ("invalid_replies", False, 0.0, agent.MAX_INVALID_REPLIES)
    assert json.loads(lost["transcript"][2]["content"])["error"].startswith("end your reply with one JSON object")
    refused = run(Reply("", 10, 0, "refusal", "bio"))
    assert (refused["stop"], refused["refusal"], refused["reward"], refused["input_tokens"]) == ("refused", "bio", 0.0, 10)
    cut = run(Reply("thinking about the plate", stop="max_tokens"))
    assert cut["stop"] == "max_tokens" and cut["turns"] == 1
    looping = run(Reply('{"tool": "control_zprime", "args": {}}', 5, 7), max_turns=3)
    assert (looping["stop"], looping["turns"], looping["output_tokens"]) == ("turn_limit", 3, 21)
    assert all(not call["ok"] for call in looping["tool_calls"])  # no plate has been run yet
    wrong = run(Reply('{"action": {"type": "dance"}}'))
    assert wrong["stop"] == "invalid_actions" and wrong["rejected_actions"] == make("dose_response", seed=SEED).MAX_INVALID
    odd = run(Reply('{"action": "run the plate"}'), max_turns=2)
    assert odd["rejected_actions"] == 2 and odd["stop"] == "turn_limit"
    bare = agent.run_episode("dose_response", SEED, agent.ScriptedModel("dose_response", "expert", SEED), tools=False)
    assert json.loads(bare["transcript"][0]["content"])["tools"] == [] and bare["finished"]


def test_each_scripted_trap_is_caught_making_its_own_mistake_and_the_expert_none():
    rules, traces, caught = knowledge.load_rules(), {}, 0
    for rule in rules:
        for trap in set(rule["traps"]) - knowledge.moot_traps(rule):
            env_name = rule["env"]
            episode = agent.run_episode(env_name, SEED, agent.ScriptedModel(env_name, trap, SEED))
            verdict = baseline.mistakes(episode, rules, traces)[trap]
            if campaign.TRAP_KIND[env_name][trap] not in baseline.FIELD_KINDS or (env_name, trap) in baseline.INTERACTIONS:
                assert verdict["made"] is None, (env_name, trap)
            elif verdict["how"] != "follows the data":
                assert verdict["made"] is True, (env_name, trap, verdict)
                caught += 1
    assert caught >= 40
    for env_name in ENVS:
        episode = agent.run_episode(env_name, SEED, agent.ScriptedModel(env_name, "expert", SEED))
        assert not [t for t, v in baseline.mistakes(episode, rules, traces).items() if v["made"]], env_name
    # the two labelled interactions are counted on the single mistake they contain
    warm = agent.run_episode("gi_fluid_stability", SEED, agent.ScriptedModel("gi_fluid_stability", "acid_quench_in_sgf_warm_overnight", SEED))
    found = baseline.mistakes(warm, rules, traces)
    assert found["acid_quench_in_sgf"]["made"] is True and found["acid_quench_in_sgf_warm_overnight"]["how"] == "counted as acid_quench_in_sgf"
    # a trap is told apart by what it sets, per run: only the off-target plates are dosed too low here
    low = agent.run_episode("selectivity_panel", SEED, agent.ScriptedModel("selectivity_panel", "on_target_range", SEED))
    sets = baseline.trap_sets("selectivity_panel", "on_target_range", traces)
    assert {run for run, _field in sets["fields"]} == {"run_plate:off1", "run_plate:off2", "run_plate:off3", "run_plate:off4"}
    assert baseline.mistakes(low, rules, traces)["on_target_range"] == {"made": True, "how": "field, exact value"}
    # what moves with the data is not what a trap sets
    assert set(baseline.trap_sets("gi_camp_potency", "no_ibmx", traces)["fields"]) == {("run_titration", "ibmx"), ("run_plate", "ibmx")}
    assert baseline.trap_sets("gi_camp_potency", "skip_titration", traces)["dropped"] == ["run_titration"]
    nothing = agent.run_episode("stability", SEED, agent.ScriptedModel("stability", "null", SEED))
    assert baseline.mistakes(nothing, rules, traces)["no_mix"] == {"made": None, "how": "ran nothing"}


def test_the_profile_puts_a_model_next_to_the_scripted_expert():
    make_model = agent.model_factory("scripted:no_mix")
    episodes = baseline.run(make_model, ["dose_response", "stability"], tasks=3, attempts=2, keep_transcripts=False)
    assert len(episodes) == 12 and "transcript" not in episodes[0] and {e["attempt"] for e in episodes} == {0, 1}
    report = baseline.profile(episodes)
    dose = report["environments"]["dose_response"]
    assert report["model"] == "scripted:no_mix" and (dose["tasks"], dose["attempts"], dose["episodes"]) == (3, 2, 6)
    assert dose["mean_reward"] < dose["expert_reward"] and dose["null_reward"] < dose["expert_reward"]
    assert dose["spread_within_task"] == 0.0 and dose["stops"] == {"done": 6} and dose["finished_rate"] == 1.0
    assert dose["mistakes"]["no_mix"]["rate"] == 1.0 and dose["mistakes"]["reuse_tips"]["rate"] == 0.0
    assert report["environments"]["stability"]["mistakes"]["no_mix"] == {
        "said": "Mix the matrix after spiking, before the time-zero sample.", "kind": "handling", "how": "field",
        "seen": 6, "made": 6, "rate": 1.0}
    assert report["environments"]["stability"]["mistakes"]["no_lloq"]["rate"] is None
    expert = baseline.profile(baseline.run(agent.model_factory("scripted:expert"), ["dose_response"], tasks=3))
    env = expert["environments"]["dose_response"]
    assert env["mean_reward"] == pytest.approx(env["expert_reward"]) and 0.0 <= env["pass_rate"] <= 1.0
    text = baseline.markdown(report)
    assert "# Baseline: scripted:no_mix" in text and "| dose_response | 3 | 2 |" in text
    assert "| stability | Mix the matrix after spiking, before the time-zero sample. | no_mix | handling | 100% | 6 of 6 | field |" in text
    assert baseline.profile([]) == {"model": None, "environments": {}}


class _FakeMessages:
    def __init__(self, responses):
        self.responses, self.calls = list(responses), []

    def create(self, **kwargs):
        self.calls.append(kwargs)
        return self.responses.pop(0)


def _response(stop_reason="end_turn", text='{"action": {"type": "report", "pEC50": {}}}', category=None):
    blocks = [SimpleNamespace(type="thinking", thinking=""), SimpleNamespace(type="text", text=text)]
    usage = SimpleNamespace(input_tokens=100, output_tokens=40, cache_read_input_tokens=900, cache_creation_input_tokens=None)
    details = SimpleNamespace(category=category) if stop_reason == "refusal" else None
    return SimpleNamespace(content=blocks, usage=usage, stop_reason=stop_reason, stop_details=details)


def test_the_anthropic_model_sends_a_plain_request_and_reads_the_reply():
    messages = _FakeMessages([_response(), _response("refusal", category="bio"), _response("max_tokens", text="cut")])
    model = agent.AnthropicModel("claude-opus-5-5", max_tokens=2000, effort="low", client=SimpleNamespace(messages=messages))
    history = [{"role": "user", "content": "task"}]
    reply = model("system text", history, 0)
    assert messages.calls[0] == {"model": "claude-opus-5-5", "max_tokens": 2000, "system": "system text", "messages": history,
                                 "cache_control": {"type": "ephemeral"}, "output_config": {"effort": "low"}}
    assert reply.text == '{"action": {"type": "report", "pEC50": {}}}' and reply.stop == "end"
    assert (reply.input_tokens, reply.output_tokens) == (1000, 40) and len(reply.raw) == 2  # the thinking block goes back too
    refused = model("system text", history, 0)
    assert (refused.stop, refused.detail, refused.text, refused.raw) == ("refusal", "bio", "", None)
    assert model("system text", history, 0).stop == "max_tokens" and model.name == "anthropic:claude-opus-5-5"
    plain = agent.AnthropicModel("claude-opus-5-5", client=SimpleNamespace(messages=_FakeMessages([_response()])))
    plain("s", history, 0)
    assert "output_config" not in plain.client.messages.calls[0] and "temperature" not in plain.client.messages.calls[0]
    # the loop replays the whole content as the assistant turn and ends the episode on the report
    looped = _FakeMessages([_response(text='{"action": {"type": "report", "pEC50": {}}}')])
    record = agent.run_episode("dose_response", SEED, agent.AnthropicModel("m", client=SimpleNamespace(messages=looped)))
    assert record["stop"] == "done" and record["input_tokens"] == 1000 and record["transcript"][1]["role"] == "assistant"


def test_a_model_is_named_by_its_provider_and_id():
    assert agent.model_factory("scripted:expert")("stability", 3).name == "scripted:expert"
    for bad in ("claude", "scripted:", "openai:gpt", ":x"):
        with pytest.raises(ValueError):
            agent.model_factory(bad)
    with pytest.raises(ValueError, match="assay_window has no scripted policy no_mix"):
        agent.model_factory("scripted:no_mix")("assay_window", 0)


def test_the_command_writes_a_report_and_the_transcripts(tmp_path, capsys):
    out = tmp_path / "baseline"
    cli.main(["baseline", "--model", "scripted:expert", "--env", "stability", "--tasks", "2", "--out", str(out)])
    printed = capsys.readouterr()
    assert "# Baseline: scripted:expert" in printed.out and "2 episodes of scripted:expert" in printed.err
    report = json.loads((out / "scripted_expert.json").read_text())
    assert report["environments"]["stability"]["episodes"] == 2
    lines = (out / "scripted_expert.transcripts.jsonl").read_text().splitlines()
    assert len(lines) == 2 and json.loads(lines[0])["transcript"][0]["role"] == "user"
    assert (out / "scripted_expert.md").read_text() == printed.out.rstrip("\n") + "\n"
    for bad in (["--model", "scripted:no_mix"], ["--model", "gpt"], ["--model", "scripted:expert", "--name", "a/b"]):
        with pytest.raises(SystemExit, match="baseline: "):
            cli.main(["baseline", "--tasks", "1", "--out", str(tmp_path / "bad")] + bad)
    assert not (tmp_path / "bad").exists()
