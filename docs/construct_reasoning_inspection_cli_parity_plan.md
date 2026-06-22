# Construct-Reasoning Inspection-Action CLI / ClawBio Parity — Status (as-built)

> **Status: DONE — already implemented and shipped.**
> This document was originally written as a forward-looking implementation plan,
> but a re-check against the tree shows the capability already exists across every
> surface (and goes beyond the original plan). The plan is therefore **superseded**
> and rewritten below as an as-built record. Honest note: the first draft of this
> plan was authored without noticing the feature was already present — the code was
> at `HEAD` the whole time.

Goal (unchanged): let headless callers (a) list the recommended inspection actions
for a construct-reasoning graph, optionally narrowed to one fact / annotation /
candidate / summary, and (b) compute / render the recommended dotplot from an
`action_id`, at parity with the GUI.

## 1. Where it landed (commits)

All present in the working tree at `HEAD = 3f717f62` (2026-06-21); tree clean for
these files.

| Commit | What it added |
| --- | --- |
| `795ff3c5` "feat (reasoning): Implemented the narrow Phase 1 bridge." | Shared-shell `ShellCommand` variants + parser + execute arms; **MCP tools**; glossary rows; ClawBio modes |
| `f0a18416` "GUI↔CLI dotplot window-math parity for construct-reasoning inspection actions" | Shared engine helpers (`bounded_center_window`, `construct_reasoning_action_dotplot_request`); GUI refactored to use them |
| `3f717f62` "feat: reasoning for lua+js" | JS + Lua bindings |

## 2. As-built map (file:line)

**Engine helpers** — `src/engine.rs`
- `pub fn bounded_center_window(..)` — `src/engine.rs:242` (now a shared free fn).
- `pub fn construct_reasoning_action_dotplot_request(action, fallback_seq_id, sequence_len)`
  — `src/engine.rs:301`; guards `action.action_kind == Dotplot` and reproduces the
  focus→viewport window math. This is the single source of truth shared with the GUI.

**Shared shell** — `src/engine_shell.rs`
- Enum variants: `ConstructReasoningListInspectionActions` (`:2483`),
  `ConstructReasoningRunInspectionAction` (`:2493`).
- `describe()` arms: `:11034`, `:11054`.
- Read-only vs mutating classification: list is read-only; run is mutating
  (`:11850`, `:40488`–`:40489`, `:40771`).
- Execute arms + output schemas:
  - list → `:38060`, schema `gentle.construct_reasoning_inspection_action_list.v1` (`:38139`).
  - run → `:38155`, schema `gentle.construct_reasoning_inspection_action_dotplot_run.v1` (`:38291`).

**Parser** — `src/engine_shell/command_parsers.rs`
- `parse_construct_reasoning_command`: `"list-inspection-actions" | "list-actions"`
  at `:8066`, run arm following; subcommand-list error updated at `:7829`.
- List filters (richer than the original plan): `--fact-id`, `--annotation-id`,
  `--candidate-id`, `--evidence-id`, `--seq-id`, `--action-kind`, `--summary-id`.
- Run flags: `--word-size` (default 12), `--step`/`--step-bp` (default 2),
  `--max-mismatches`, `--tile-bp`, `--id`/`--dotplot-id`,
  `--render-svg`/`--output-svg`.

**GUI** — refactored to share the math (true single-source parity)
- `open_construct_reasoning_dotplot_action` (`src/main_area_dna.rs:13014`) now calls
  `construct_reasoning_action_dotplot_request` (`:13023`).
- `bounded_center_window` in `src/main_area_dna/auxiliary_workspaces.rs:1380`
  delegates to `crate::engine::bounded_center_window`.

**MCP** — `src/mcp_server.rs`
- Tools `construct_reasoning_inspection_actions` (`:822`) and
  `construct_reasoning_run_inspection_action` (`:854`); shell-contract mapping at
  `:1206`–`:1210`.

**JS / Lua**
- JS: `src/js_interface.rs:248` (list), `:277` (run); test `:2072`.
- Lua: `src/lua_interface.rs:375` (list), `:418` (run); test `:2748`.

**ClawBio** — `integrations/clawbio/skills/gentle-cloning/`
- First-class request modes (not raw `mode=shell`): `construct-reasoning-list-inspections`
  and `construct-reasoning-run-inspection` — `gentle_cloning.py:42`–`43`, routing
  at `:888` / `:1684` / `:1705`.
- `INTENTS.json:1356` (list intent), `:1392` (run intent).
- Example requests: `examples/request_construct_reasoning_list_inspections.json`,
  `examples/request_construct_reasoning_run_inspection_dotplot.json`.
- `README.md`, `SKILL.md`, and `tests/test_gentle_cloning.py` updated.

**Docs / generated**
- `docs/glossary.json:3877` (list), `:3898` (run).
- `docs/gui_cli_mcp_parity.md:316`–`317` (matrix rows present, regenerated).
- `docs/cli.md:627`, `:1112`–`:1116` (full flag list for both commands).

**Tests** — `src/engine_shell/tests.rs`
- Parser round-trips: `:1895` (list, exercises all 7 filters), `:1920` (run).
- Execution: `execute_construct_reasoning_inspection_action_commands_list_and_run_dotplot`
  at `:2453` — builds a fixture graph, lists + filters actions, runs the dotplot,
  asserts the focus range is windowed and the run schema.
- Glossary smoke tests auto-cover the new rows; `tests/parity_matrix_freshness.rs`
  covers the matrix.

## 3. How the result compares to the original plan

The implementation **meets or exceeds** every point of the original plan:

| Original plan said | Reality |
| --- | --- |
| Add 2 shell subcommands (single source of truth) | Done — `list-inspection-actions` (+ `list-actions` alias) and `run-inspection-action`. |
| Filters `--fact/--annotation/--summary/--kind` | Done and richer: also `--candidate-id`, `--evidence-id`, `--seq-id`, `--action-kind`. |
| Engine helpers for lookup + window math; share `bounded_center_window` | Done — `construct_reasoning_action_dotplot_request` + shared `bounded_center_window` in the engine. |
| GUI refactor to share the math — "recommended but optional" | Done (commit `f0a18416`): GUI now calls the shared helper. |
| MCP / JS / Lua — marked **out of scope / optional** | Done — full MCP tools + JS + Lua bindings. |
| ClawBio via `mode=shell`, plus examples/intents | Exceeded — first-class ClawBio modes, intents, and example requests. |
| CLI stays shell-only (not forwarded) | As decided — reachable via `gentle_cli shell 'construct-reasoning ...'`, consistent with the family. |
| Run output schema `…inspection_action_run.v1` (guessed) | Actual: `gentle.construct_reasoning_inspection_action_dotplot_run.v1`. |

## 4. Verification (confirm green)

Nothing functional remains; this is verification only.

```bash
cargo test -p gentle --lib engine_shell::tests::execute_construct_reasoning_inspection_action 2>&1 | tail -20
cargo test -p gentle --lib engine_shell::tests::glossary 2>&1 | tail -20
cargo test --test parity_matrix_freshness
# ClawBio skill tests
python -m pytest integrations/clawbio/skills/gentle-cloning/tests/test_gentle_cloning.py -k inspection

# Manual smoke (needs a state with a built graph + its sequence loaded)
gentle_cli --state /tmp/demo.json shell 'construct-reasoning list-inspection-actions <GRAPH_ID> --fact-id <FACT_ID>'
gentle_cli --state /tmp/demo.json shell 'construct-reasoning run-inspection-action <GRAPH_ID> <ACTION_ID> --render-svg /tmp/r.svg'
```

## 5. Optional polish (not required)

- `docs/glossary.json` documents the common filter subset
  (`--fact-id/--annotation-id/--summary-id`) for the list command; `docs/cli.md`
  carries the full set. The glossary usage could be expanded to all seven filters
  for completeness — but the smoke tests already pass on the documented subset, so
  this is cosmetic.
- No new `action_kind` beyond `Dotplot` exists; the run helper already guards on it,
  so future kinds slot in cleanly.
