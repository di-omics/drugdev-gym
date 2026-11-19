"""Agent: run a language model through an environment, one JSON reply at a time.

The protocol is plain text, so any chat model can be measured the same way. The model is
given the task (the environment's own spec), the analysis tools it may call, and the first
observation. Every reply has to end in one JSON object:

    {"action": {...}}                    act in the environment
    {"tool": "name", "args": {...}}      analyse data it has already been shown

An action is answered with the next observation, a tool call with its result. The hidden
truth and the reward never go to the model.

A model is any callable `model(system, messages, attempt) -> Reply`. ScriptedModel wraps a
reference policy, for tests and dry runs. AnthropicModel calls the Anthropic API and needs
the optional SDK (`pip install -e ".[llm]"`) and credentials the SDK can find. Nothing in
this module runs a model by itself, and no model has been run on these environments yet.
"""
from __future__ import annotations

import json
from dataclasses import dataclass

from .policies import POLICIES
from .registry import make, task_id
from .rollouts import _plain
from .tools import TOOLS, call_tool

SYSTEM = (
    "You are running one laboratory assay task in a simulator. Work only from the task description and the "
    "data you are shown. You may reason in plain text first, but every reply must end with exactly one JSON "
    "object and nothing after it. To act in the environment, send {\"action\": {...}} with an action the task "
    "describes. To analyse data you have already been shown, send {\"tool\": \"<name>\", \"args\": {...}} with one "
    "of the listed tools. The next message is the next observation or the tool result. The episode ends when "
    "you send the final action the task describes."
)
MAX_TURNS = 24  # model replies per episode
MAX_INVALID_REPLIES = 4  # replies with no usable JSON object before the episode is stopped


@dataclass
class Reply:
    """One model reply. stop is "end", "max_tokens" or "refusal"; raw is what to send back as the
    assistant turn when the provider wants more than the text (for example its thinking blocks)."""
    text: str
    input_tokens: int = 0
    output_tokens: int = 0
    stop: str = "end"
    detail: str | None = None
    raw: object = None


def _default(value):
    if hasattr(value, "tolist"):
        return value.tolist()  # numpy scalar or array
    raise TypeError(f"not JSON: {type(value).__name__}")


def extract_json(text: str) -> dict | None:
    """The last JSON object in a reply that has an "action" or a "tool" key, or None."""
    decoder = json.JSONDecoder()
    start = len(text)
    while True:
        start = text.rfind("{", 0, start)
        if start < 0:
            return None
        try:
            value, _end = decoder.raw_decode(text, start)
        except ValueError:
            continue
        if isinstance(value, dict) and ("action" in value or "tool" in value):
            return value


class ScriptedModel:
    """A reference policy answering through the protocol. It is handed each observation as an object
    (observe), so what it sends is the policy's own action, not a reading of the text."""

    def __init__(self, env_name: str, policy_name: str, seed: int):
        self.name = f"scripted:{policy_name}"
        self.policy = POLICIES[env_name][policy_name](seed)
        self.obs: dict = {}

    def observe(self, obs: dict) -> None:
        self.obs = obs

    def __call__(self, system: str, messages: list[dict], attempt: int) -> Reply:
        return Reply(json.dumps({"action": self.policy.act(self.obs)}, default=_default))


class AnthropicModel:
    """A Claude model through the Anthropic Messages API.

    Thinking is left to the model's default and no sampling parameters are sent, which is what the
    current models accept. A refusal is returned as a refusal and nothing is retried on another
    model: a baseline measures one model. The whole response content is sent back as the assistant
    turn, so thinking blocks are replayed unchanged.
    """

    def __init__(self, model: str, max_tokens: int = 16000, effort: str | None = None, client=None):
        if client is None:
            try:
                import anthropic
            except ImportError as err:
                raise RuntimeError('the Anthropic SDK is not installed: pip install -e ".[llm]"') from err
            client = anthropic.Anthropic()  # the SDK finds the credentials: an API key, a token or a login profile
        self.name, self.model, self.max_tokens, self.effort, self.client = f"anthropic:{model}", model, max_tokens, effort, client

    def __call__(self, system: str, messages: list[dict], attempt: int) -> Reply:
        kwargs = {"model": self.model, "max_tokens": self.max_tokens, "system": system, "messages": messages,
                  "cache_control": {"type": "ephemeral"}}  # the growing history is the same prefix every turn
        if self.effort:
            kwargs["output_config"] = {"effort": self.effort}
        response = self.client.messages.create(**kwargs)
        usage = response.usage
        read = (usage.input_tokens + (getattr(usage, "cache_read_input_tokens", 0) or 0)
                + (getattr(usage, "cache_creation_input_tokens", 0) or 0))
        if response.stop_reason == "refusal":  # check before reading content
            details = getattr(response, "stop_details", None)
            return Reply("", read, usage.output_tokens, "refusal", getattr(details, "category", None))
        text = "".join(block.text for block in response.content if block.type == "text")
        stop = "max_tokens" if response.stop_reason == "max_tokens" else "end"
        return Reply(text, read, usage.output_tokens, stop, raw=response.content)


def model_factory(spec: str, max_tokens: int = 16000, effort: str | None = None):
    """A function (env_name, seed) -> model, from "scripted:<policy>" or "anthropic:<model id>"."""
    kind, _, name = spec.partition(":")
    if kind == "scripted" and name:
        def scripted(env_name: str, seed: int):
            if name not in POLICIES[env_name]:
                raise ValueError(f"{env_name} has no scripted policy {name}")
            return ScriptedModel(env_name, name, seed)
        return scripted
    if kind == "anthropic" and name:
        shared = AnthropicModel(name, max_tokens, effort)
        return lambda env_name, seed: shared
    raise ValueError('model must be "scripted:<policy>" or "anthropic:<model id>"')


def run_episode(env_name: str, seed: int, model, rig: str = "nominal", attempt: int = 0, tools: bool = True,
                max_turns: int = MAX_TURNS) -> dict:
    """One episode of one model on one task. Returns the record: what it did, what it scored, how it ended.

    stop is "done" (the environment ended the episode), "invalid_actions" (the environment gave up
    on it), "invalid_replies", "turn_limit", "refused" or "max_tokens". An episode that does not
    end by "done" scores 0.
    """
    env = make(env_name, seed=seed, rig=rig)
    obs = env.reset()
    spec = obs.pop("spec", None)
    toolset = TOOLS.get(env_name, []) if tools else []
    first = {"task": spec, "tools": [t.describe() for t in toolset], "observation": _plain(obs)}
    messages = [{"role": "user", "content": json.dumps(first)}]
    transcript = [dict(messages[0])]
    seen = [obs]  # every observation so far, for the tools
    record = {"task": task_id(env_name, rig, seed), "env": env_name, "rig": rig, "seed": seed, "attempt": attempt,
              "model": getattr(model, "name", type(model).__name__), "stop": None, "finished": False, "reward": 0.0,
              "turns": 0, "invalid_replies": 0, "rejected_actions": 0, "input_tokens": 0, "output_tokens": 0,
              "actions": [], "tool_calls": [], "info": None, "refusal": None}
    total, fresh = 0.0, True

    def say(payload: dict) -> None:
        text = json.dumps(_plain(payload))
        messages.append({"role": "user", "content": text})
        transcript.append({"role": "user", "content": text})

    while record["stop"] is None:
        if record["turns"] >= max_turns:
            record["stop"] = "turn_limit"
            break
        if fresh and hasattr(model, "observe"):
            model.observe(seen[-1])
        fresh = False
        reply = model(SYSTEM, messages, attempt)
        record["turns"] += 1
        record["input_tokens"] += reply.input_tokens
        record["output_tokens"] += reply.output_tokens
        messages.append({"role": "assistant", "content": reply.raw if reply.raw is not None else reply.text})
        transcript.append({"role": "assistant", "content": reply.text})
        if reply.stop == "refusal":
            record["stop"], record["refusal"] = "refused", reply.detail
            break
        parsed = extract_json(reply.text)
        if parsed is None:
            if reply.stop == "max_tokens":
                record["stop"] = "max_tokens"
                break
            record["invalid_replies"] += 1
            if record["invalid_replies"] >= MAX_INVALID_REPLIES:
                record["stop"] = "invalid_replies"
                break
            say({"error": 'end your reply with one JSON object: {"action": {...}} or {"tool": "<name>", "args": {...}}'})
            continue
        if "tool" in parsed:
            result, ok = call_tool(toolset, parsed.get("tool"), parsed.get("args"), seen)
            record["tool_calls"].append({"tool": parsed.get("tool"), "args": _plain(parsed.get("args")), "ok": ok})
            say({"tool_result": result} if ok else {"error": result})
            continue
        action = parsed["action"]
        try:
            obs, reward, done, info = env.step(action)
        except Exception as err:  # an action no environment check foresaw: tell the model, do not crash the run
            obs, reward, done, info = {"error": f"the environment could not read that action ({type(err).__name__})"}, 0.0, False, {}
        accepted = "error" not in obs
        record["actions"].append({"action": _plain(action), "accepted": accepted})
        record["rejected_actions"] += not accepted
        total += reward
        if done:
            truncated = bool(info.get("truncated"))
            record["stop"] = "invalid_actions" if truncated else "done"
            record["finished"], record["info"] = not truncated, _plain(info)
            break
        if accepted:
            seen.append(obs)
            fresh = True
        say({"observation": obs})
    record["reward"] = float(total) if record["finished"] else 0.0
    record["transcript"] = transcript
    return record
