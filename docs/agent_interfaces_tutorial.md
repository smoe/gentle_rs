# GENtle Agent Interfaces Tutorial

Last updated: 2026-03-04

This is a self-standing tutorial for using GENtle with AI/agent systems.

Goal: make it explicit who runs what, where, for which purpose, and how the
interfaces differ:

- deterministic CLI/shared shell
- MCP (`gentle_mcp`)
- in-app Agent Assistant
- external coding agents (for example OpenAI Codex asked to "use GENtle")

This document is conceptual + operational. Protocol-level details remain in
`docs/protocol.md`.

Plain-language note used in this file:

- "fixed command format" means a command/tool with defined inputs and
  predictable outputs (usually JSON).

## 1) One mental model first

GENtle has one deterministic execution core (engine + shared shell commands).

Agent interfaces are wrappers around that core.

- deterministic layer:
  - `gentle_cli ...`
  - `gentle_cli shell '...'`
  - GUI `Shell` panel
  - MCP tools calling shared engine commands
- conversational layer:
  - Agent Assistant (GUI or `agents ask`)
  - external coding agents (Codex, etc.) that decide what commands to run

Important: reproducibility comes from fixed commands with explicit arguments, not
from free-form prompting.

## 2) Who runs what where (exact role map)

### A) Bench scientist in GENtle GUI

Runs in: GENtle desktop app window.

Uses:

- GUI controls/windows
- GUI `Shell` for deterministic commands
- Agent Assistant window for planning + command suggestions

Purpose:

- interactive exploration and visualization
- human-in-the-loop command execution

### B) Analyst/engineer in terminal

Runs in: terminal (bash/zsh/etc.).

Uses:

- `gentle_cli` direct commands
- `gentle_cli shell '...'` shared parser

Purpose:

- scripts, CI, reproducible pipelines, explicit state files

### C) External AI orchestrator with MCP support

Runs in: MCP-capable client/runner.

Uses:

- `gentle_mcp` JSON-RPC `tools/list` and `tools/call`

Purpose:

- controlled tool-calling loops with structured envelopes
- capability discovery/negotiation before execution (`tools/list`,
  `capabilities`, `help`)
- integration with tool-based agent frameworks

### D) External coding agent (for example Codex) asked to "use GENtle"

Runs in: coding-agent environment (editor/terminal automation context), not
inside GENtle itself.

Uses:

- shell commands (`cargo`, `gentle_cli`, file edits)
- optionally MCP if explicitly wired by user

Purpose:

- automate repo tasks, run commands, patch docs/code, orchestrate workflows

Note: this is not the same as GENtle Agent Assistant. It is a separate agent
runtime controlling your environment.

## 3) How AI environments discover GENtle functionality

This question has two parts: static discovery and runtime discovery.

Static discovery (what exists and how it is intended to be used):

- repository instructions and architecture docs:
  - `AGENTS.md`
  - `docs/architecture.md`
  - `docs/protocol.md`
  - `docs/cli.md`
  - `docs/gui.md`
  - `docs/glossary.json` (canonical command catalog source)

Runtime discovery (what this build can do right now):

- CLI/shared shell discovery:
  - `gentle_cli capabilities`
  - `gentle_cli help --format markdown --interface cli-direct`
  - `gentle_cli help --format json --interface all`
  - `gentle_cli agents list`
- MCP discovery:
  - `tools/list` against `gentle_mcp`
  - `tools/call(name="capabilities", ...)`
  - `tools/call(name="help", arguments={...})`

For Codex-style coding agents specifically:

- in-repo: they usually learn first from `AGENTS.md` + docs, then verify using
  CLI/MCP discovery calls above.
- outside repo context: they only know what you provide at runtime, so include
  explicit commands and state paths in prompts.

## 4) Local LLMs (Msty, Jan, Ollama/OpenAI-compatible)

Local models are deployment choices, not a separate GENtle operation model. They
use the same interfaces as hosted models.

Primary route (recommended for most users):

- Agent Assistant (GUI) or `gentle_cli agents ask ...`
- transport: `native_openai_compat`
- endpoint: local OpenAI-compatible HTTP base URL (for example Jan/Msty/Ollama
  gateway)
- model: concrete local model id (for example `deepseek-r1:8b`)

GUI flow:

1. Open `Tools -> Agent Assistant...`.
2. Select local system template (for example `Msty Local (template)`).
3. Set/confirm base URL.
4. Click `Discover models`.
5. Pick one discovered model (model `unspecified` is intentionally blocked for
   execution).
6. Ask the assistant as usual.

CLI example:

```bash
gentle_cli agents ask local_llama_compat --prompt "summarize current project and suggest next command" --base-url http://localhost:11964 --model deepseek-r1:8b
```

Alternative route for tool-calling stacks:

- If your local-agent runtime supports MCP, connect it to `gentle_mcp` and use
  `tools/list`/`tools/call`.
- If it does not support MCP tools, use chat/planning mode and have it emit
  explicit `gentle_cli` commands for you to run.

## 5) Fast comparison matrix

| Topic | CLI/shared shell | MCP | Agent Assistant (GUI/`agents ask`) | External coding agent (Codex-style) |
|---|---|---|---|---|
| Where it runs | Terminal or GUI Shell | MCP client + `gentle_mcp` server | Inside GENtle (or CLI `agents ask`) | External agent runtime |
| Primary unit | command | tool call (`tools/call`) | prompt + suggestions | free-form tasks, then commands |
| Output form | direct JSON/text/markdown | JSON-RPC envelope + structured content | assistant text/questions/suggested commands | whatever the agent emits |
| Determinism | high | high | medium for planning, high once commands are executed | depends on how strictly it uses commands |
| Safety gate | command intent | `confirm=true` on mutating tools | per-suggestion execution (`ask`/`auto`) | governed by that runtime/user policy |
| Best use | scripts/reproducibility | tool ecosystems | user guidance + quick planning | code/docs automation + orchestration |

## 6) Tutorial flow: same task through each interface

Task example: prepare a helper genome and run BLAST for one query sequence.

### Path A: deterministic CLI/shared shell

Run directly:

```bash
gentle_cli helpers prepare "Plasmid pUC19 (local)" --cache-dir data/helper_genomes --timeout-secs 600
gentle_cli helpers blast "Plasmid pUC19 (local)" ACGTACGT --task blastn-short --max-hits 20 --cache-dir data/helper_genomes
```

Or via shared shell route:

```bash
gentle_cli shell 'helpers prepare "Plasmid pUC19 (local)" --cache-dir data/helper_genomes --timeout-secs 600'
gentle_cli shell 'helpers blast "Plasmid pUC19 (local)" ACGTACGT --task blastn-short --max-hits 20 --cache-dir data/helper_genomes'
```

Who executes: you (or script/CI).

### Path B: MCP tool-calling

Client calls `gentle_mcp` tools (abbreviated concept flow):

1. `tools/list`
2. `tools/call(name="capabilities", arguments={...})`
3. `tools/call(name="op", arguments={..., "confirm": true})` for mutating ops
4. `tools/call(name="blast_async_start", ...)`
5. poll `blast_async_status`

Who executes: MCP client orchestrator.

Where state lives: resolved MCP `state_path` and engine persistence rules.

### Path C: GENtle Agent Assistant

1. Open `Tools -> Agent Assistant...`
2. Select system from catalog.
3. Optionally set API key/base URL/model/timeouts.
4. Paste a structured prompt (or use `Insert` template).
5. Review suggestions.
6. Execute selected suggestions (`Run` per row or explicit execute options).

Who executes:

- model generates suggestions
- GENtle executes only selected/allowed suggestions through shared shell

### Path D: external coding agent (Codex-style) asked to use GENtle

Typical user request:

- "Run helper BLAST via gentle_cli and summarize top hits."

Agent behavior:

- chooses terminal commands (for example `gentle_cli ...`)
- runs them in your environment
- summarizes outputs / edits files

Key distinction:

- this is outside GENtle’s own assistant transport and governance UI
- reliability depends on how strictly the coding agent is instructed to use
  deterministic commands and explicit state paths

## 7) MCP vs command line (practical difference)

Both can be deterministic and adapter-equivalent. The difference is integration
style.

- command line:
  - best for humans/scripts/CI using plain process invocations
  - minimal envelope overhead
- MCP:
  - best for tool-calling systems that already use JSON-RPC (the standard MCP
    request/response format)
  - explicit mutating gate (`confirm=true`) and typed tool names

If you already have shell automation, CLI is simpler.
If you are integrating into an MCP-native agent stack, MCP is simpler.

## 8) MCP vs "ask Codex to use GENtle"

MCP:

- GENtle is a declared tool endpoint
- narrow, explicit input/output format per tool
- easier to bound and audit tool invocations

Codex-style external agent:

- general-purpose environment automation
- can do much more than GENtle (file edits, build/test/docs, etc.)
- not inherently constrained to GENtle tool schemas unless you enforce it

Use MCP when you want strict tool interface discipline.
Use external coding agent when you want broader repo/system automation.

## 9) Agent Assistant vs "GENtle Shell" vs "prompting in Agent Assistant"

These are often conflated:

- GUI `Shell`:
  - deterministic command executor
  - you type exact command syntax
- Agent Assistant prompt box:
  - conversational planning input to a model
  - model returns text/questions/suggested commands
- Agent Assistant execution:
  - optional step that runs selected suggestions via shared shell commands

So:

- prompting in Agent Assistant is not the same as executing shell commands
- execution determinism starts when concrete commands are run

## 10) Minimal-success guidance for colleagues

For internal preview, default to this policy:

1. Ask Assistant for explicit commands only.
2. Keep execution mode at `ask-before-run`.
3. Run one suggestion at a time.
4. Keep a state file/project path explicit for reproducibility.
5. If a run matters, re-run same command through `gentle_cli` and record output.

## 11) Common confusion checklist

- "I used Agent Assistant, why is this less reproducible?"
  - Because prompt text is not a stable execution format; commands are.
- "MCP and CLI gave different results."
  - They should not for equivalent routed commands; this indicates a parity bug.
- "Codex succeeded but I cannot replay."
  - Ask for explicit final `gentle_cli` commands and state path, then rerun.
- "My local LLM cannot execute tool calls."
  - Use Agent Assistant/`agents ask` with `native_openai_compat`, or have the
    model output explicit `gentle_cli` commands.

## 12) When to choose which interface

Choose CLI/shared shell when:

- you need exact replayability and scripts.

Choose MCP when:

- your orchestrator is tool-calling-first and expects JSON-RPC envelopes.

Choose Agent Assistant when:

- you need help translating goals into concrete GENtle commands.

Choose external coding agent (Codex-style) when:

- you want broader engineering automation around GENtle, not only GENtle tool
  calls.

Choose local OpenAI-compatible Agent Assistant when:

- you want private/local inference while keeping the same GENtle execution
  routes.

## Related docs

- `docs/agent_interface.md` (concise route overview)
- `docs/gui.md` (Agent Assistant UI operations)
- `docs/cli.md` (command syntax and examples)
- `docs/protocol.md` (schemas and field definitions)
