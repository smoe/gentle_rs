# Quick-Start: Driving GENtle From Claude

Claude can help with GENtle in two useful ways. Inside GENtle, Claude acts as
an Agent Assistant: the user types an intent, Claude returns
`gentle.agent_response.v1`, and GENtle executes selected shared-shell command
suggestions internally. Outside GENtle, Claude Code or another Claude surface
can drive the working-tree CLI with `scripts/dev-gentle-cli`, inspect JSON
results, and help extend GENtle when a route is missing. Both paths rely on the
same engine and shared-shell contracts described by the
[adapter parity policy](architecture.md#1-product-intent).

## Which Path To Use

Internal Agent Assistant / `agents ask` keeps GENtle in control: Claude returns
shared-shell command suggestions and GENtle reviews and executes them. External
Claude Code uses `scripts/dev-gentle-cli` to discover routes, run commands,
inspect JSON, and edit source. MCP is only for agents that already want a typed
tool palette through `gentle_mcp`.

## Setup In Ten Minutes

- Clone the GENtle repository and enter the repo root.
- Confirm Rust is available:

```sh
rustc --version
```

- Warm the checkout once; the first build can take minutes:

```sh
cargo build
```

- Have Claude ready as either a configured GENtle agent system, Claude Code in
  the checkout, or a browser/API surface for copy-paste one-offs.

Important: `assets/agent_systems.json` now ships a native Claude entry,
`anthropic_claude_sonnet_native`. Real internal Claude use needs
`ANTHROPIC_API_KEY` or an Anthropic Console API key pasted in the GUI. Claude
Code or Claude.ai subscription/login tokens are not Anthropic API keys. The
offline `builtin_echo` system proves the GENtle loop, but it is not Claude.

Do not start by wiring `gentle_mcp` into `~/.claude.json`. That can be useful
for a dedicated MCP agent, but it pins Claude to an installed binary. The
checkout loop below uses the current source tree, so improvements made during a
session are immediately available.

## Pre-Flight

For the checkout-driven path, run:

```sh
scripts/dev-gentle-cli doctor --agent
```

Recorded in a Codex sandbox with a writable temp target; normal users do not need that environment variable.

```json
{
  "schema": "gentle.agent_dev_doctor.v1",
  "ok": true,
  "cargo": {"ok": true, "program": "cargo", "version": "cargo 1.96.0-beta.6 (9fb171546 2026-04-14)", "error": null},
  "binary": {"ok": true, "binary_exists": true, "fresh": true, "error": null},
  "state": {"ok": true, "path": ".gentle_state.json", "parent_writable": true, "file_writable_if_exists": true, "error": null},
  "examples_docs_check": {"ok": true, "status_code": 0, "stdout_tail": "{\n  \"example_count\": 36,\n  \"source_dir\": \"docs/examples/workflows\",\n  \"status\": \"ok\"\n}", "error": null}
}
```

A healthy reply has `ok: true`, a fresh binary, writable state path, and
`examples_docs_check.ok: true`.

For the internal Agent Assistant path, prove the local agent catalog first:

```sh
scripts/dev-gentle-cli agents preflight builtin_echo
```

```json
{"schema":"gentle.agent_preflight.v1","system_id":"builtin_echo","system_label":"Built-in Echo (demo)","transport":"builtin_echo","available":true,"catalog_path":"assets/agent_systems.json","warnings":[]}
```

For real Claude, use `anthropic_claude_sonnet_native` and provide
`ANTHROPIC_API_KEY`. A configuration-only preflight should succeed before trying
live requests.

## Scenario 1: Claude Inside GENtle

The internal route is not "Claude runs a terminal command." Instead, Claude
returns a structured assistant response. GENtle then displays any suggested
commands and executes selected rows through the shared shell.

The response shape is:

```json
{
  "schema": "gentle.agent_response.v1",
  "assistant_message": "Short explanation for the user.",
  "questions": [],
  "suggested_commands": [{
    "title": "Scan restriction sites",
    "rationale": "EcoRI and SmaI are requested on the inline sequence.",
    "command": "features restriction-scan --sequence-text ATGCGTGAATTCTTAGGCCCGGGAAATTTCCCGGGATCGATCGAATTCTTACGATCGATCGATCG --topology linear --id-hint toy_restriction_70bp --enzyme EcoRI --enzyme SmaI --max-sites-per-enzyme 5",
    "execution": "ask"
  }]
}
```

The `command` value is a GENtle shared-shell command. It is not an OS shell
command and it should not include `gentle_cli`. Prefer `execution: "ask"` so
the user reviews before running. Use `execution: "auto"` only for safe,
read-only commands where the system explicitly allows auto execution. Use
`execution: "chat"` when the command is explanatory and should not run.

Agent-suggested commands must not recursively invoke `agents ask`,
`agents plan`, or `agents execute-plan`; GENtle blocks those nested routes.

The offline demo proves the suggestion shape:

```sh
scripts/dev-gentle-cli agents ask builtin_echo \
  --prompt "ask: capabilities" \
  --no-state-summary
```

It returns one suggested command:

```json
{
  "schema": "gentle.agent_ask_result.v1",
  "summary": {"suggested_command_count": 1, "executed_count": 0},
  "invocation": {
    "response": {
      "schema": "gentle.agent_response.v1",
      "suggested_commands": [{"title": "Confirm suggestion (demo)", "command": "capabilities", "execution": "ask"}]
    }
  }
}
```

Because the suggestion is `ask`, GENtle records it but does not execute it until
the user selects it.

## Worked Example A: Inline Restriction Scan

User intent:

```text
Find EcoRI and SmaI sites in this synthetic 65 bp sequence:
ATGCGTGAATTCTTAGGCCCGGGAAATTTCCCGGGATCGATCGAATTCTTACGATCGATCGATCG
```

Claude should first learn whether the operation exists:

```sh
scripts/dev-gentle-cli capabilities |
  jq '.capability_registry[]
      | select(.name=="FindRestrictionSites")
      | {name, source, mutating, cli, mcp, js, lua, inline_operand_ok, engine_operations}'
```

Real output:

```json
{
  "name": "FindRestrictionSites",
  "source": "engine_operation",
  "mutating": "false",
  "cli": "shell_passthrough",
  "mcp": "shell_passthrough",
  "js": "shell_passthrough",
  "lua": "shell_passthrough",
  "inline_operand_ok": true,
  "engine_operations": ["FindRestrictionSites"]
}
```

Then Claude suggests this shared-shell command:

```sh
features restriction-scan --sequence-text ATGCGTGAATTCTTAGGCCCGGGAAATTTCCCGGGATCGATCGAATTCTTACGATCGATCGATCG --topology linear --id-hint toy_restriction_70bp --enzyme EcoRI --enzyme SmaI --max-sites-per-enzyme 5
```

For checkout verification, run the same command through the wrapper:

```sh
scripts/dev-gentle-cli shell "features restriction-scan --sequence-text ATGCGTGAATTCTTAGGCCCGGGAAATTTCCCGGGATCGATCGAATTCTTACGATCGATCGATCG --topology linear --id-hint toy_restriction_70bp --enzyme EcoRI --enzyme SmaI --max-sites-per-enzyme 5"
```

Real result, trimmed to the matched rows:

```json
{
  "result": {
    "messages": ["Restriction-site scan on 'toy_restriction_70bp' matched 4 site(s) across 2 enzyme(s) over 0..65"],
    "restriction_site_scan": {
      "schema": "gentle.restriction_site_scan.v1",
      "matched_site_count": 4,
      "rows": [
        {"enzyme_name": "EcoRI", "recognition_sequence": "GAATTC", "recognition_start_0based": 6, "recognition_end_0based_exclusive": 12},
        {"enzyme_name": "SmaI", "recognition_sequence": "CCCGGG", "recognition_start_0based": 17, "recognition_end_0based_exclusive": 23},
        {"enzyme_name": "SmaI", "recognition_sequence": "CCCGGG", "recognition_start_0based": 29, "recognition_end_0based_exclusive": 35},
        {"enzyme_name": "EcoRI", "recognition_sequence": "GAATTC", "recognition_start_0based": 42, "recognition_end_0based_exclusive": 48}
      ]
    }
  }
}
```

Claude's user-facing summary can be short: EcoRI appears at `6..12` and
`42..48`; SmaI appears at `17..23` and `29..35`, using 0-based half-open
coordinates.

## Scenario 2: Claude Outside GENtle

Claude Code can drive the checkout directly. Open it in the repo root and give
it this instruction:

```text
For GENtle actions in this conversation, use scripts/dev-gentle-cli.
Use scripts/dev-gentle-cli capabilities to discover routes.
Use scripts/dev-gentle-cli shell "..." for shared-shell commands.
Use --state PATH when persistence is needed.
Do not use gentle_mcp unless I explicitly ask for the MCP route.
```

A short turn is: user asks for an EcoRI/SmaI scan, Claude inspects
`capabilities`, runs `scripts/dev-gentle-cli shell "features restriction-scan
..."`, and reports the four sites. This path is best for agents that can read
source, run tests, and extend GENtle.

## Worked Example B: When A Command Is Missing

For a tiny GC-content command, Claude should first verify the route is missing:

```sh
rg -n "gc-content|gc_content|GcContents" src/bin/gentle_cli.rs src/gc_contents.rs
```

Then it should read the existing helper:

```sh
sed -n '1,80p' src/gc_contents.rs
```

This sketch was compile-checked locally with `cargo check -q` and then
restored; it is not committed by this quick-start:

```diff
+use gentle::gc_contents::GcContents;
+
+        "gc-content" => {
+            let sequence = args.get(cmd_idx + 1)
+                .ok_or_else(|| "gc-content requires SEQUENCE".to_string())?;
+            let gc = GcContents::new_from_sequence(sequence.as_bytes());
+            print_json(&json!({
+                "schema": "gentle.gc_content.v1",
+                "length_bp": sequence.len(),
+                "regions": gc.regions(),
+            }))
+        }
```

For a real change, Claude would keep the route in the shared shell if it should
be available through GUI Shell, CLI shell, MCP shell-like routes, JS, Lua, and
Agent Assistant suggestions. Direct CLI-only commands are appropriate only for
local utility behavior.

## Differences From Driving Via Codex

Codex and Claude use the same GENtle command language: discover routes through
`capabilities`, execute deterministic work through shared shell or typed
operation/workflow calls, and follow `AGENTS.md`. Claude Code often benefits
from explicit repo instructions; Codex tends to take longer self-directed runs.

## Other Claude Surfaces

Claude in a browser can produce shared-shell commands that you paste into
GENtle's Shell or Agent Assistant review flow. Anthropic API integrations are
viable for custom bridges; a GENtle Agent Assistant bridge must return
`gentle.agent_response.v1`.

## Why MCP Is Not The Invariant Route Here

`gentle_mcp` presents GENtle as a typed tool palette. That is valuable when an
external agent already speaks MCP and should call stable tools with explicit
JSON schemas.

The Agent Assistant loop is different: Claude returns shared-shell suggestions
for GENtle to review and execute internally. The development loop is also
different: Claude drives the current checkout through `scripts/dev-gentle-cli`.
Those two loops stay close to the source tree and to the GUI-as-inspection
direction.

## Where To Go Next

Use [Agent development loop](agent_dev_loop.md) for the external checkout
loop, [Agent interface](agent_interface.md) for schema details,
[Architecture](architecture.md) for parity policy, [Protocol](protocol.md) for
operation contracts, [ClawBio integration](../integrations/clawbio/skills/gentle-cloning/SKILL.md)
for curated agent-to-agent invocation, and [AGENTS.md](../AGENTS.md) for
repository-local agent rules.
