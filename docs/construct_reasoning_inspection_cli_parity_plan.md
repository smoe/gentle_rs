# Construct-Reasoning Inspection-Action CLI / ClawBio Parity Plan

> Implementation plan for Codex. Goal: bring the **CLI and ClawBio** surfaces to
> parity with the **GUI** for construct-reasoning *inspection actions* — i.e. let
> a headless caller (a) list the recommended inspection actions for a reasoning
> graph (optionally narrowed to one fact / annotation / summary), and (b)
> compute / render the recommended dotplot directly from an `action_id`.

## 1. Goal & scope

Today the GUI already knows how to surface and act on the dotplot inspection
actions that the engine attaches to every construct-reasoning graph. The CLI and
the ClawBio skill have **no** way to do this — there is no shell command that
exposes `ConstructReasoningGraph.inspection_actions`, and no way to turn an
`action_id` into a computed/rendered dotplot.

Add two new `construct-reasoning` shared-shell subcommands (the single source of
truth that the CLI and ClawBio both route through):

1. `construct-reasoning list-inspection-actions GRAPH_ID [filters]`
   — list the recommended inspection actions for a graph, optionally filtered to
   the actions that back a given fact / annotation / summary.
2. `construct-reasoning run-inspection-action GRAPH_ID ACTION_ID [--render-svg PATH] [tuning]`
   — resolve the action, compute the recommended dotplot (mirroring the GUI's
   focus→viewport window math), store it, and optionally render an SVG.

Then wire the doc/test/ClawBio surfaces so the new commands are discoverable and
covered.

### Out of scope
- No new GUI behaviour (the GUI already does this; we only optionally refactor it
  to share the window math — see §5.1).
- No new `action_kind` values. `ConstructReasoningInspectionActionKind::Dotplot`
  is the only kind today; keep the design open to future kinds but do not add any.
- No changes to how graphs are *built* — we only read existing
  `inspection_actions`.

## 2. Key facts about the current code (verified)

Data model — `crates/gentle-protocol/src/construct_reasoning.rs`:
- `ConstructReasoningGraph` carries `inspection_actions: Vec<ConstructReasoningInspectionAction>`
  plus `facts`, `annotation_candidates`, `annotation_candidate_summaries`, `seq_id`, `graph_id`.
- `ConstructReasoningInspectionAction` fields used here:
  `action_id`, `action_kind` (`Dotplot`), `button_label`, `hover_text`, `seq_id`,
  `mode: DotplotMode`, `focus_start_0based`, `focus_end_0based_exclusive`,
  `driving_evidence_ids`, `source_fact_ids`, `source_annotation_ids`,
  `source_summary_ids`, `context_tags`, `repeat_family_provenance`.
- Schema const: `CONSTRUCT_REASONING_INSPECTION_ACTION_SCHEMA = "gentle.construct_reasoning_inspection_action.v1"`.

Engine — `src/engine.rs`:
- Actions are generated during `normalize_construct_reasoning_graph`
  (`graph.inspection_actions = Self::construct_reasoning_build_inspection_actions(&graph)`,
  ~line 17398), so any stored/returned graph already has them populated.
- Read accessors already exist:
  - `pub fn construct_reasoning_graph(&self, graph_id: &str) -> Result<ConstructReasoningGraph, EngineError>` (~17558)
  - `pub fn list_construct_reasoning_graph_summaries(&self, seq_id_filter: Option<&str>)` (~17486)
  - `pub fn list_dotplot_views(&self, seq_id_filter: Option<&str>) -> Vec<DotplotViewSummary>` (~23938)
  - `pub fn get_dotplot_view(&self, dotplot_id: &str) -> Result<DotplotView, EngineError>` (~23986)
- `ProjectState.sequences: HashMap<SeqId, DNAsequence>` (~1212) → sequence length
  via `self.state().sequences.get(seq_id).map(|d| d.len())`.
- `Operation::ComputeDotplot { seq_id, reference_seq_id, span_start_0based,
  span_end_0based, reference_span_start_0based, reference_span_end_0based, mode,
  word_size, step_bp, max_mismatches, tile_bp, store_as }` (~3711).
- `Operation::RenderDotplotSvg { .. }` exists (~2574) and is what the
  `render-dotplot-svg` direct CLI command + `RenderDotplotSvg` shell variant use.

GUI reference behaviour (what we are matching) — `src/main_area_dna.rs` &
`src/main_area_dna/auxiliary_workspaces.rs`:
- Per-source filtering helpers (match these semantics exactly):
  - `construct_reasoning_dotplot_actions_for_fact` → keep actions whose
    `source_fact_ids` contains the fact id.
  - `construct_reasoning_dotplot_actions_for_annotation` → keep actions whose
    `source_annotation_ids` contains the annotation id **and** whose
    `source_summary_ids` is empty.
  - `construct_reasoning_dotplot_actions_for_summary` → keep actions whose
    `source_summary_ids` contains the summary id.
- `open_construct_reasoning_dotplot_action` (~src/main_area_dna.rs:12611) is the
  authoritative focus→dotplot window derivation. Reproduce it (see §5.1):
  ```text
  len            = sequence length (bp); error if 0
  focus_start    = action.focus_start_0based.min(len-1)
  focus_end      = action.focus_end_0based_exclusive.max(focus_start+1).min(len)
  focus_span     = (focus_end - focus_start).max(1)
  target_span_bp = if len <= 200 { len } else { (focus_span*3).clamp(200, len) }
  half_window    = (target_span_bp - 1) / 2
  center         = (focus_start + focus_span/2).min(len-1)
  (vp_start, vp_end) = bounded_center_window(len, center, half_window)  // see below
  span_start_0based = vp_start ; span_end_0based = vp_end
  mode           = action.mode
  store_as       = "{normalize(seq_id)}_reasoning_{focus_start+1}_{focus_end}_{mode_tag}"
                   where mode_tag ∈ {self, revcomp, pair_forward, pair_revcomp}
  ```
- `bounded_center_window(len, center, half_window) -> Option<(start, end_excl)>`
  lives at `src/main_area_dna/auxiliary_workspaces.rs:1375` (currently
  `pub(super)`); it returns a window of width `2*half_window+1` clamped to
  `[0, len)` and shifted to stay in-bounds near the edges.
- Default `word_size = 12`, `step_bp = 2` (GUI `dotplot_ui` defaults; also the
  `dotplot compute` parser defaults at `src/engine_shell/command_parsers.rs:6384`).

Shell plumbing — the `construct-reasoning` family already exists but is **not**
top-level forwarded (only reachable via `gentle_cli shell '...'`), exactly like
its sibling subcommands. Touch points:
- `ShellCommand` enum variants: `src/engine_shell.rs:2441`–`2456`
  (`ConstructReasoningListGraphs`, `ConstructReasoningShowGraph`, …).
- `describe()` arm (human text): `src/engine_shell.rs:10877`+.
- Parser: `parse_construct_reasoning_command` at
  `src/engine_shell/command_parsers.rs:7826` (dispatched from the
  `"construct-reasoning" | "construct_reasoning" | "constructreasoning"` arm
  ~line 23878).
- Execute arm: `src/engine_shell.rs:36630`+ (`ConstructReasoningListGraphs`, …);
  note the read-only vs mutating classification lists at ~38835 and ~39116.

Doc/test machinery (must stay green):
- `docs/glossary.json` is the hand-maintained source of truth (included via
  `include_str!`). Shell help + capability registry + parity matrix derive from it.
- `src/engine_shell/tests.rs::glossary_cli_usage_smoke_commands_parse` and
  `glossary_cli_usage_flags_parse_one_by_one` parse **every** glossary usage
  string (and each flag) through `parse_shell_line`. New glossary rows must parse.
- Parity matrix: `docs/gui_cli_mcp_parity.md` is generated by
  `gentle_protocol::render_gui_cli_mcp_parity_matrix_markdown()` from
  `docs/glossary.json` + `docs/parity_matrix_overrides.json`;
  `tests/parity_matrix_freshness.rs` asserts freshness. Regenerate with
  `scripts/regenerate_parity_matrix.sh` after editing the glossary.
- `src/bin/gentle_cli.rs::test_shell_forwarded_commands_are_documented_in_glossary`
  only constrains the `SHELL_FORWARDED_COMMANDS` list — irrelevant unless we add
  forwarding (we are not, see §6).

ClawBio — `integrations/clawbio/skills/gentle-cloning/`:
- The Python skill (`gentle_cloning.py`) supports `mode=shell`, which forwards an
  arbitrary shared-shell line to `gentle_cli shell '...'`. So **both new commands
  become reachable from ClawBio for free** once they exist in the shared shell.
  Parity work is discoverability: example requests, `INTENTS.json`, README.
- The in-app bridge `src/app/clawbio_bridge.rs` is a generic transport and needs
  **no** change.

## 3. New command grammar

### 3.1 `list-inspection-actions`
```
construct-reasoning list-inspection-actions GRAPH_ID
    [--fact FACT_ID] [--annotation ANNOTATION_ID] [--summary SUMMARY_ID]
    [--kind dotplot]
```
- `GRAPH_ID` required, non-empty.
- The three source filters are optional and may be combined; an action matches if
  it satisfies **all** supplied filters, using the GUI per-source semantics above
  (annotation filter additionally requires `source_summary_ids` empty).
- `--kind` optional; defaults to all kinds. Only `dotplot` is valid today; parse
  via `ConstructReasoningInspectionActionKind` (reject unknown values).
- Output (read-only, `state_changed = false`):
  ```json
  {
    "schema": "gentle.construct_reasoning_inspection_action_list.v1",
    "graph_id": "...",
    "seq_id": "...",
    "filters": { "fact_id": null, "annotation_id": null, "summary_id": null, "kind": null },
    "action_count": 2,
    "actions": [ <ConstructReasoningInspectionAction>, ... ]
  }
  ```
- Errors: empty `GRAPH_ID` → InvalidInput; unknown graph → propagate the
  `engine.construct_reasoning_graph` NotFound error.

### 3.2 `run-inspection-action`
```
construct-reasoning run-inspection-action GRAPH_ID ACTION_ID
    [--word-size N] [--step N] [--max-mismatches N] [--tile-bp N]
    [--id DOTPLOT_ID]
    [--render-svg PATH]
    [--flex-track ID] [--display-threshold N] [--intensity-gain N]
```
- `GRAPH_ID`, `ACTION_ID` required, non-empty.
- Defaults: `word_size = 12`, `step_bp = 2`, `max_mismatches = 0`, `tile_bp = None`
  (match `dotplot compute`).
- `--id` overrides the derived `store_as` dotplot id.
- `--render-svg PATH` (optional): after compute, apply `Operation::RenderDotplotSvg`
  for the resolved seq/dotplot id to `PATH`; `--flex-track/--display-threshold/
  --intensity-gain` only meaningful with `--render-svg` (mirror
  `render-dotplot-svg` flag semantics in
  `src/bin/gentle_cli/rendering.rs` / `RenderDotplotSvg` parser).
- Behaviour:
  1. `graph = engine.construct_reasoning_graph(GRAPH_ID)?`
  2. find `action` in `graph.inspection_actions` by `action_id` (NotFound error,
     message lists a few valid ids if missing).
  3. resolve `seq_id` (the action's `seq_id`; fall back to `graph.seq_id` if the
     action's is empty), get `len` from `state().sequences`; error if absent/empty.
  4. derive `ComputeDotplot` params via the shared helper (§5.1).
  5. apply `ComputeDotplot`, then select the stored view (same selection logic as
     `DotplotCompute` handler: by `--id` if given else newest for seq).
  6. if `--render-svg`, apply `RenderDotplotSvg` and capture the output path.
- Output (`state_changed = true` when the dotplot store changed; compare the
  `DOTPLOT_ANALYSIS_METADATA_KEY` metadata before/after exactly like the existing
  `DotplotCompute` arm):
  ```json
  {
    "schema": "gentle.construct_reasoning_inspection_action_run.v1",
    "graph_id": "...",
    "action": <ConstructReasoningInspectionAction>,
    "resolved": {
      "seq_id": "...", "mode": "self_forward",
      "span_start_0based": 0, "span_end_0based": 123,
      "word_size": 12, "step_bp": 2, "dotplot_id": "..."
    },
    "result": <ComputeDotplot op result>,
    "dotplot": <DotplotViewSummary | null>,
    "svg_path": "artifacts/....svg"   // present only with --render-svg
  }
  ```
- Note on "open": opening a window is GUI-only and already implemented. For the
  headless surfaces "open" maps to "compute + return the stored view"; `--render-svg`
  is the materialised-artifact form ClawBio consumes.

## 4. Engine API additions (`src/engine.rs`)

Add small, testable, `pub` helpers near the existing construct-reasoning
accessors (~17558). These centralise lookup + window math so the GUI and the
shell command cannot drift.

1. Action lookup:
   ```rust
   pub fn construct_reasoning_inspection_action(
       &self, graph_id: &str, action_id: &str,
   ) -> Result<(ConstructReasoningGraph, ConstructReasoningInspectionAction), EngineError>
   ```
   Loads the graph, finds the action by `action_id` (trimmed), NotFound otherwise.

2. Filtered listing (pure, mirrors GUI filters) — either a `pub` method or a free
   function reused by both list command and tests:
   ```rust
   pub fn construct_reasoning_inspection_actions_filtered(
       graph: &ConstructReasoningGraph,
       fact_id: Option<&str>,
       annotation_id: Option<&str>,
       summary_id: Option<&str>,
       kind: Option<ConstructReasoningInspectionActionKind>,
   ) -> Vec<ConstructReasoningInspectionAction>
   ```

3. Window derivation (the §2 math), pure and unit-testable:
   ```rust
   pub struct ConstructReasoningDotplotRequest {
       pub seq_id: SeqId,
       pub mode: DotplotMode,
       pub span_start_0based: usize,
       pub span_end_0based: usize,
       pub store_as: String,
   }

   pub fn construct_reasoning_action_dotplot_request(
       action: &ConstructReasoningInspectionAction,
       fallback_seq_id: &str,
       sequence_len: usize,        // 0 => Err
   ) -> Result<ConstructReasoningDotplotRequest, EngineError>
   ```
   Move `bounded_center_window` into a shared spot (engine or
   `gentle-protocol`) so this function and the GUI use one implementation. Keep
   the exact arithmetic from `open_construct_reasoning_dotplot_action` and
   `bounded_center_window` so outputs are identical.

> If a smaller blast radius is preferred, the window math may live entirely inside
> the shell handler — but then add a focused unit test asserting it matches the
> GUI helper for a few cases. Sharing via the engine is the recommended path.

## 5. Implementation steps

### 5.1 Share the window math (recommended)
- Promote `bounded_center_window` to a shared, testable function (engine assoc.
  fn or `gentle-protocol` free fn). Update `auxiliary_workspaces.rs:1375` to
  delegate (keep `pub(super)` wrapper if convenient).
- Implement `construct_reasoning_action_dotplot_request` (§4.3).
- Optional but ideal: refactor `open_construct_reasoning_dotplot_action`
  (`src/main_area_dna.rs:12611`) to call the new helper for span/mode/store_as,
  keeping the viewport-setting + `compute_primary_dotplot`/`open_dotplot_window`
  GUI steps. This guarantees true GUI↔CLI parity and is the cleanest outcome.

### 5.2 Shell command variants — `src/engine_shell.rs`
- Add two `ShellCommand` variants next to the other `ConstructReasoning*`
  variants (~2441):
  ```rust
  ConstructReasoningListInspectionActions {
      graph_id: String,
      fact_id: Option<String>,
      annotation_id: Option<String>,
      summary_id: Option<String>,
      kind: Option<ConstructReasoningInspectionActionKind>,
  },
  ConstructReasoningRunInspectionAction {
      graph_id: String,
      action_id: String,
      word_size: usize,
      step_bp: usize,
      max_mismatches: usize,
      tile_bp: Option<usize>,
      dotplot_id: Option<String>,
      render_svg_path: Option<String>,
      flex_track_id: Option<String>,
      display_threshold: Option<f64>,
      intensity_gain: Option<f64>,
  },
  ```
- Add `describe()` arms (~10877) with human-readable strings (follow the existing
  `ConstructReasoningShowGraph` phrasing).
- Add the execute arms (~36630, in the same `match` as the other
  `ConstructReasoning*` handlers):
  - List: load graph, filter via §4.2 helper, emit the §3.1 JSON,
    `state_changed = false`.
  - Run: follow §3.2; emulate the `DotplotCompute` arm's before/after
    `DOTPLOT_ANALYSIS_METADATA_KEY` comparison for `state_changed`, and its
    "select stored view" logic; render SVG when requested.
- Add both variants to the read-only / mutating classification helper lists
  (~38835 and ~39116): list is **read-only**, run is **mutating** (it writes a
  stored dotplot view), matching how `DotplotCompute` is treated.

### 5.3 Parser — `src/engine_shell/command_parsers.rs`
- Extend `parse_construct_reasoning_command` (7826): add `"list-inspection-actions"`
  and `"run-inspection-action"` match arms; update the subcommand-list error
  strings (7829, 8129) to mention them.
- Use the existing `parse_option_path` helper for value flags and the existing
  numeric-parse idioms (see the `build-protein-dna-handoff` arm and the
  `dotplot compute` parser for `--word-size/--step/--max-mismatches/--tile-bp`).
- Add a kind parser (small local fn) mapping `"dotplot"` →
  `ConstructReasoningInspectionActionKind::Dotplot`, rejecting others.
- `parse_dotplot_mode` already exists if needed.

### 5.4 Glossary — `docs/glossary.json`
Add two command entries adjacent to the existing `construct-reasoning *` rows
(after `export-graph`). Mirror the exact flag spellings the parser accepts so the
glossary smoke tests pass. Suggested entries:
```json
{
  "aliases": ["construct_reasoning list-inspection-actions"],
  "path": "construct-reasoning list-inspection-actions",
  "summary": "List the recommended dotplot inspection actions attached to a construct-reasoning graph, optionally narrowed to one fact, annotation, or summary.",
  "usage": "construct-reasoning list-inspection-actions GRAPH_ID [--fact FACT_ID] [--annotation ANNOTATION_ID] [--summary SUMMARY_ID] [--kind dotplot]"
},
{
  "aliases": ["construct_reasoning run-inspection-action"],
  "path": "construct-reasoning run-inspection-action",
  "summary": "Compute (and optionally render) the recommended dotplot for one construct-reasoning inspection action id.",
  "usage": "construct-reasoning run-inspection-action GRAPH_ID ACTION_ID [--word-size N] [--step N] [--max-mismatches N] [--tile-bp N] [--id DOTPLOT_ID] [--render-svg OUTPUT.svg] [--flex-track ID] [--display-threshold F] [--intensity-gain F]"
}
```
- Match the surrounding entries' field set/order exactly (check whether sibling
  rows carry `interfaces`/`category` keys and replicate).
- If the smoke value generator can't satisfy a flag (e.g. it feeds a non-numeric
  token to a numeric flag), add a targeted skip in
  `skip_glossary_flag_parse` (`src/engine_shell/tests.rs:456`) — but first confirm
  it's actually needed; the existing `dotplot compute` numeric flags pass, so the
  generator already handles `N`-style placeholders.

### 5.5 Regenerate parity matrix + (optional) overrides
- If sibling `construct-reasoning` rows have entries in
  `docs/parity_matrix_overrides.json`, add matching override rows for the two new
  paths (same ClawBio "curated skill intents" note) for a consistent matrix;
  otherwise the defaults render fine.
- Run `scripts/regenerate_parity_matrix.sh` and commit the updated
  `docs/gui_cli_mcp_parity.md` so `tests/parity_matrix_freshness.rs` passes.

### 5.6 ClawBio parity — `integrations/clawbio/skills/gentle-cloning/`
- Add example requests under `examples/` using `mode=shell`, e.g.:
  - `request_construct_reasoning_list_inspection_actions.json`
    ```json
    { "mode": "shell", "shell_line": "construct-reasoning list-inspection-actions GRAPH_ID" }
    ```
  - `request_construct_reasoning_run_inspection_action_render_svg.json`
    ```json
    { "mode": "shell", "shell_line": "construct-reasoning run-inspection-action GRAPH_ID ACTION_ID --render-svg artifacts/reasoning_dotplot.svg" }
    ```
    (use a real graph/action id from a fixture if the example suite is executed in
    CI; otherwise keep placeholders consistent with other `mode=shell` examples).
- Register them where the example catalog is enumerated:
  `INTENTS.json`, `catalog_entry.json`, and/or the experimental-followup catalog
  (`experimental_followup_request_catalog.json`) following the pattern of an
  existing `mode=shell` example. Update `SKILL.md` / `README.md` prose to mention
  "inspect reasoning graph dotplot recommendations".
- No change to `gentle_cloning.py` routing is required (the `shell` mode already
  forwards arbitrary shared-shell lines). If the skill maintains a curated
  allow/suggestion list of shell lines, add the two new lines there.
- ClawBio already post-processes `*.svg` artifacts into PNGs; `--render-svg`
  output will flow through that path automatically.

### 5.7 Docs
- `docs/cli.md`: document both subcommands under the construct-reasoning section.
- `docs/construct_reasoning_plan.md`: tick the inspection-action CLI/ClawBio
  parity item (it already references inspection actions and the self-dotplot path).

## 6. CLI exposure decision

Keep the new subcommands **shell-only** (reachable via
`gentle_cli shell 'construct-reasoning ...'`), consistent with the entire existing
`construct-reasoning` family — `construct-reasoning` is deliberately **not** in
`SHELL_FORWARDED_COMMANDS` (`src/bin/gentle_cli.rs:507`). This is the lower-risk
choice and matches sibling commands, so "CLI parity" is satisfied without new
forwarding.

> Optional follow-up (separate change, not required here): add `construct-reasoning`
> to `SHELL_FORWARDED_COMMANDS` to allow `gentle_cli construct-reasoning ...`
> directly. That would forward the *whole* family and require every
> `construct-reasoning *` glossary path to satisfy
> `test_shell_forwarded_commands_are_documented_in_glossary` — broader than this
> task. Recommend deferring.

## 7. Tests

- `src/engine/tests.rs`: 
  - `construct_reasoning_action_dotplot_request` window math (edge cases: short
    sequence ≤200 bp, focus near start/end, revcomp/pair modes, empty sequence →
    Err). If the GUI is refactored, assert the helper matches a known
    `open_construct_reasoning_dotplot_action` expectation.
  - filtered listing semantics (fact / annotation-with-summary-exclusion / summary).
- `src/engine_shell/tests.rs`:
  - parser round-trips for both new subcommands incl. each flag.
  - execution: build a graph with `build-protein-dna-handoff` (or a fixture graph
    with `inspection_actions`), then `list-inspection-actions` returns the actions
    and filters narrow them; `run-inspection-action` produces a stored dotplot view
    and (with `--render-svg` to a tempdir) writes an SVG file.
  - the glossary smoke tests (`glossary_cli_usage_smoke_commands_parse`,
    `glossary_cli_usage_flags_parse_one_by_one`) automatically cover the new
    glossary rows once added — make sure they stay green.
- `tests/parity_matrix_freshness.rs`: passes after regeneration.
- ClawBio: if `integrations/clawbio/skills/gentle-cloning/tests/` validates example
  files against the schema, ensure the new examples conform.

## 8. Verification checklist (commands to run)

```bash
# Build + targeted tests
cargo test -p gentle --lib engine_shell:: 2>&1 | tail -40
cargo test -p gentle --lib engine::tests 2>&1 | tail -40
cargo test --test parity_matrix_freshness

# Regenerate generated docs
scripts/regenerate_parity_matrix.sh
git diff --stat docs/gui_cli_mcp_parity.md

# Manual smoke (needs a project state with a built graph + its sequence loaded)
gentle_cli --state /tmp/demo.json shell 'construct-reasoning list-graphs'
gentle_cli --state /tmp/demo.json shell 'construct-reasoning list-inspection-actions <GRAPH_ID>'
gentle_cli --state /tmp/demo.json shell 'construct-reasoning list-inspection-actions <GRAPH_ID> --fact <FACT_ID>'
gentle_cli --state /tmp/demo.json shell 'construct-reasoning run-inspection-action <GRAPH_ID> <ACTION_ID> --render-svg /tmp/r.svg'

# ClawBio (mode=shell) — same lines flow through the skill
```

Acceptance: the two `shell` lines above return the §3 JSON shapes; the `run` line
writes `/tmp/r.svg`; `cargo test` and the freshness test pass; ClawBio examples
validate.

## 9. Risks & notes
- **Window-math drift** is the main correctness risk — mitigated by extracting one
  shared helper (§5.1) and unit-testing it.
- **Glossary smoke parsing**: keep usage flag spellings identical to the parser;
  reuse existing numeric-flag idioms so the smoke value generator is satisfied.
- **`run-inspection-action` preconditions**: the action's `seq_id` must be loaded
  in the project state for `ComputeDotplot` to succeed — return a clear error if
  not (same failure mode as `dotplot compute`).
- **`state_changed`**: compute writes a stored dotplot view, so the run command is
  mutating; the list command is read-only. Classify both correctly so
  `--state` persistence and any read-only guards behave like `dotplot compute` /
  `dotplot list`.
