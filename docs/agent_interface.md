# GENtle Agent Interface

Last updated: 2026-02-27

This guide explains how agents can control GENtle and how the available
interfaces differ.

## Why this exists

GENtle exposes multiple agent-facing routes:

- deterministic command execution through CLI/shared shell
- deterministic tool execution through MCP
- guided assistant interaction through the Agent Assistant bridge

These routes are related, but they are not identical in transport, ergonomics,
and safety flow.

## Agent-facing routes

### 1) CLI and shared shell

Primary binaries and routes:

- `gentle_cli` direct commands (`op`, `workflow`, `genomes`, `tracks`, etc.)
- `gentle_cli shell '...'` shared shell command parser/executor
- GUI Shell panel (same shared shell parser/executor as `gentle_cli shell`)

What it is good for:

- deterministic automation scripts
- CI jobs and reproducible local workflows
- explicit operation sequencing with low overhead

Key properties:

- same engine contract as GUI/JS/Lua
- machine-readable JSON outputs
- strong reproducibility and auditability

### 2) MCP (`gentle_mcp`)

MCP is the tool-based route for external AI clients that speak JSON-RPC over
stdio.

Current MCP tool families include:

- engine/state/help tools:
  - `capabilities`
  - `state_summary`
  - `op` (requires `confirm=true`)
  - `workflow` (requires `confirm=true`)
  - `help`
- UI-intent tools:
  - `ui_intents`
  - `ui_intent`
  - `ui_prepared_genomes`
  - `ui_latest_prepared`

What it is good for:

- integrating GENtle into MCP-compatible agent runners
- explicit tool-call loops (`tools/list`, `tools/call`)
- deterministic tool result handling with standard MCP envelopes

Key properties:

- MCP tools are thin wrappers over existing shared engine/shell paths
- no MCP-only biology logic branch
- mutating tools require explicit confirmation (`confirm=true`)
- UI-intent tools are currently non-mutating query/intent routes

### 3) Agent Assistant bridge (`agents ...` and GUI Agent Assistant)

Agent Assistant runs configured external/internal AI systems and can return:

- assistant text
- follow-up questions
- suggested shell commands (with execution intent: `chat|ask|auto`)

Entry points:

- CLI/shared shell: `agents list`, `agents ask`
- GUI: `Tools -> Agent Assistant...`

What it is good for:

- interactive planning and translation from human request to deterministic
  commands
- human-in-the-loop execution control
- quick exploration with optional state summary context

Key properties:

- not a direct replacement for deterministic interfaces
- produces suggestions that can be executed through shared shell contracts
- recursion guardrail blocks nested `agents ask` execution from suggested
  commands

## CLI vs MCP vs Agent Assistant

| Topic | CLI/shared shell | MCP | Agent Assistant |
|---|---|---|---|
| Transport | process args/stdin/stdout | JSON-RPC over stdio | catalog-driven agent transport + shared shell execution |
| Best use | scripts and deterministic automation | tool-based external agent integration | interactive assistant with suggestion flow |
| Mutating safety gate | command-level intent | explicit `confirm=true` on mutating tools | per-suggestion execution policy (`ask`/`auto`) |
| Output model | direct JSON/text/markdown | MCP envelope + structuredContent | assistant payload + optional execution reports |
| Deterministic parity target | canonical | canonical via wrappers | uses canonical shell routes for execution |

## The "agent prompt" offered by GENtle

GENtle provides prompt templates in the GUI Agent Assistant. These templates
help users write better requests, but they are not an execution interface by
themselves.

Current template set includes:

- Structured (recommended)
- Candidate between anchors
- BLAST specificity check
- Track import + prioritization
- Macro/template authoring

Important distinction:

- Interface (CLI/MCP/shared shell): executable contract and deterministic result
- Prompt template: guidance text for the model request

In short:

- prompt quality affects assistant output quality
- execution determinism comes from the underlying shell/engine contracts

## Typical usage patterns

### Pattern A: deterministic script-first

1. Use `gentle_cli` or `gentle_cli shell` directly.
2. Store commands/workflows in version control.
3. Re-run unchanged for reproducible outputs.

### Pattern B: external MCP agent orchestrator

1. Connect to `gentle_mcp`.
2. Discover tools with `tools/list`.
3. Call tools deterministically (`ui_*`, `state_summary`, `op`, `workflow`).
4. Require explicit `confirm=true` for mutating calls.

### Pattern C: interactive assistant then execute

1. Use GUI Agent Assistant or `agents ask`.
2. Ask for explicit command suggestions.
3. Execute selected suggestions through shared shell.
4. Keep `ask-before-run` by default for safety.

## Safety and governance notes

- Keep business/biology logic in shared engine contracts only.
- Keep adapters thin (CLI/MCP/GUI assistant should not fork biology behavior).
- Use explicit confirmation for mutating routes.
- Prefer deterministic machine-readable outputs for automation.

## Related manuals

- `docs/gui.md` (GUI usage and Agent Assistant UI details)
- `docs/cli.md` (CLI and MCP operational commands)
- `docs/protocol.md` (protocol contracts and schema-level details)
- `docs/architecture.md` (architecture invariants and parity rules)
