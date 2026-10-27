# Macros and Routines

<!-- manual-counterparts:start -->
Related GUI reference: [Routine Assistant](../gui/routine-assistant.md).
These topic links do not imply identical feature coverage across interfaces.
<!-- manual-counterparts:end -->

- Shared shell:
  - Commands:
    - `macros run [--transactional] [--file PATH | SCRIPT_OR_@FILE]`
    - `macros instance-list`
    - `macros instance-show MACRO_INSTANCE_ID`
    - `macros template-list`
    - `macros template-show TEMPLATE_NAME`
    - `macros template-put TEMPLATE_NAME (--script SCRIPT_OR_@FILE|--file PATH) [--description TEXT] [--details-url URL] [--param NAME|NAME=DEFAULT ...] [--input-port PORT_ID:KIND[:one|many][:required|optional][:description]] [--output-port PORT_ID:KIND[:one|many][:required|optional][:description]]`
    - `macros template-delete TEMPLATE_NAME`
    - `macros template-import PATH`
    - `macros template-run TEMPLATE_NAME [--bind KEY=VALUE ...] [--transactional] [--validate-only]`
    - `routines list [--catalog PATH] [--family NAME] [--status NAME] [--tag TAG] [--query TEXT] [--seq-id SEQ_ID]`
    - `routines explain ROUTINE_ID [--catalog PATH] [--seq-id SEQ_ID]`
    - `routines compare ROUTINE_A ROUTINE_B [--catalog PATH] [--seq-id SEQ_ID]`

Workflow macro commands (`gentle_cli shell 'macros ...'`):

- `macros run [--transactional] [--file PATH | SCRIPT_OR_@FILE]`
  - Executes semicolon/newline-separated shell statements.
  - Existing file paths are auto-loaded even without `@` (shebang-friendly).
  - Supports transactional rollback (`--transactional`) when any statement fails.
    This includes parse/nested-macro errors after earlier steps: pre-run project
    state, journal and undo/redo history are restored, while revisions advance
    to invalidate stale work. External files, resources and running jobs are
    not rolled back. The same guarantee applies to `candidates macro` and
    transactional template execution.
  - Designed for full cloning workflows through `op ...` and `workflow ...`
    statements (Digest/Ligation/PCR/ExtractRegion/container ops, etc.).
  - All runs persist a lineage macro-instance record:
    - success: status `ok`
    - failure: status `failed` (or `cancelled` when cancellation-like error text is detected)
  - Successful runs return `macro_instance_id`; failed runs include
    `macro_instance_id=...` in the error message.
- `macros instance-list`
  - Lists recorded macro-instance lineage rows.
  - Response schema: `gentle.lineage_macro_instances.v1`.
- `macros instance-show MACRO_INSTANCE_ID`
  - Shows one recorded macro-instance lineage row.
  - Response schema: `gentle.lineage_macro_instance.v1`.
- `macros template-list`
  - Lists persisted workflow macro templates.
- `macros template-show TEMPLATE_NAME`
  - Shows one persisted workflow template definition.
- `macros template-put TEMPLATE_NAME (--script SCRIPT_OR_@FILE|--file PATH) [--description TEXT] [--details-url URL] [--param NAME|NAME=DEFAULT ...] [--input-port PORT_ID:KIND[:one|many][:required|optional][:description]] [--output-port PORT_ID:KIND[:one|many][:required|optional][:description]]`
  - Creates/updates a named workflow macro template in project metadata.
  - Placeholders in script use `${param_name}` and must be declared via `--param`.
  - Optional `--details-url URL` records external protocol/reference details for
    catalog display.
  - Optional typed port contracts (`--input-port` / `--output-port`) are stored
    in template metadata and used as first-class preflight contract source.
- `macros template-delete TEMPLATE_NAME`
  - Deletes one persisted workflow template.
- `macros template-import PATH`
  - Imports workflow macro templates from:
    - one pack JSON file (`gentle.cloning_patterns.v1`)
    - one single-template JSON file (`gentle.cloning_pattern_template.v1`)
    - one directory tree (recursive `*.json` import)
  - If one template fails validation, no imported template changes are kept.
- `macros template-run TEMPLATE_NAME [--bind KEY=VALUE ...] [--transactional] [--validate-only]`
  - Expands a named template with provided bindings/defaults, then executes it as
    a workflow macro script.
  - Executes typed preflight checks before mutation using template port
    contracts when present, otherwise routine catalog mapping.
  - Preflight now includes richer semantics:
    - cross-port alias/collision checks,
    - input sequence vs input container consistency checks,
    - sequence-anchor semantic checks against the bound input sequence when unambiguous,
    - Gibson-family overlap checks (adjacent suffix/prefix validation against
      configured overlap length on Gibson routines),
    - restriction-family digest checks (enzyme-name resolution, duplicate-enzyme
      misuse, recognition-site presence across bound input sequences, and
      digest parameter sanity for fragment indices/extract range).
  - `--validate-only` runs preflight only and never mutates state.
  - Mutating runs now always record a macro-instance lineage row:
    - success: `ok`
    - preflight/execute failure: `failed`/`cancelled`
  - Successful runs return `macro_instance_id`; failed runs include
    `macro_instance_id=...` in the error message.

Typed routine catalog command (`gentle_cli routines ...` or `gentle_cli shell 'routines ...'`):

- `routines list [--catalog PATH] [--family NAME] [--status NAME] [--tag TAG] [--query TEXT] [--seq-id SEQ_ID]`
  - Lists typed cloning routines from catalog JSON (`gentle.cloning_routines.v1`).
  - `--family`, `--status`, `--tag`: exact case-insensitive filters.
  - `--query`: case-insensitive substring match across id/title/family/status/template/tags/summary plus explainability metadata fields.
  - `--seq-id`: refreshes construct reasoning for that sequence and feeds any
    variant-derived assay preferences into the planning estimate / ranking.
  - Default catalog path: `assets/cloning_routines.json`.
- `routines explain ROUTINE_ID [--catalog PATH] [--seq-id SEQ_ID]`
  - Returns structured explainability payload for one routine.
  - Response schema: `gentle.cloning_routine_explain.v1`.
  - Includes resolved alternatives plus purpose/mechanism/requirements, contraindications, disambiguation questions, and failure modes.
  - `--seq-id` adds the same construct-reasoning planning context used by
    `routines list` / `routines compare`, including any variant-derived assay
    preferences and the resulting planning estimate for the explained routine.
  - `planning.protocol_overview` contains source-bound JSON and `svg`, schema
    `gentle.planning_protocol_overview.v1`. It shows recorded synthetic bench
    prerequisites, parallel branches and the typical critical path, not
    executed lineage or verified molecule transfers. Columns are dependency
    depth, not time. Unmodelled routines return `unavailable`; legacy reports
    without this field do not imply an empty valid procedure. The same shared
    report is used by Routine Assistant and external/inner agents. See
    [tutorial 09.04](../tutorial/09-04_cloning_time_cost_and_protocol_suppliers.md)
    for extracting its SVG without executing a routine.
- `routines compare ROUTINE_A ROUTINE_B [--catalog PATH] [--seq-id SEQ_ID]`
  - Returns deterministic side-by-side comparison payload for two routines.
  - Response schema: `gentle.cloning_routine_compare.v1`.
  - Includes shared/unique vocabulary tags, difference-matrix rows, and merged disambiguation questions.
  - Also includes planning-aware estimate rows:
    `estimated_time_hours`, `estimated_cost`, `local_fit_score`,
    `composite_meta_score`.
    Additional axes expose hands-on/unattended/elapsed ranges, time basis and
    stock-draw, missing-material purchasing proxy, machine and base costs.
    Unknown ranges print `unknown`; a legacy scalar and step-model path have
    different labelled bases, not like-for-like duration evidence.
  - `--seq-id` applies the same sequence-aware construct-reasoning preference
    context used by `routines list`.

<!-- manual-page-nav:start -->

---

[CLI Manual](../cli.md) | [Previous: Batch Workflows](batch.md) | [Next: Feature Annotation](features.md)
<!-- manual-page-nav:end -->
