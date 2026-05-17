---
chapter_id: "pcr_selection_batch_primer_pairs_offline"
title: "Selection-first PCR batch primer design (offline)"
tier: "core"
example_id: "pcr_selection_batch_primer_pairs_offline"
source_example: "docs/examples/workflows/pcr_selection_batch_primer_pairs_offline.json"
example_test_mode: "always"
executed_during_generation: true
---

# Selection-first PCR batch primer design (offline)

Use paint-first PCR Designer controls on local TP73 context to queue ROIs and run deterministic batch primer-pair reports, with optional extracted-region copies.

This chapter is PCR-focused and biology-anchored without requiring genome mapping mechanics. You work directly on local TP73 locus sequence (`test_files/tp73.ncbi.gb`), open the dedicated PCR Designer, paint semantic regions on the linear map (`ROI` green, upstream primer window red, downstream primer window blue), queue promoter-proximal targets around TP73-AS2, and run one shared primer-constraint set across queued regions. Results remain deterministic (`..._r01`, `..._r02`) and can be shown/exported for handoff before entering the separate luciferase cloning tutorial.

## When This Routine Is Useful

- You want a mouse-first way to assign PCR ROI and primer-window coordinates with immediate visual feedback.
- You want to combine one ad hoc painted selection with one or more annotation-driven queue rows in one PCR batch run.
- You want to amplify TP73-AS2 promoter-proximal windows from local TP73 sequence context.
- You want to design primer pairs for several selected regions without retyping ROI coordinates.
- You want PCR-first preparation that cleanly feeds later promoter-luciferase planning.
- You need copied sequence artifacts for each PCR target before forwarding to other workflows.
- You want deterministic report IDs for review and export.

## What You Learn

- Capture PCR queue regions from paint-first map interaction and selected-feature fallback actions.
- Use semantic painted intervals to preview pair-PCR geometry before execution.
- Run deterministic multi-region primer design in one batch action with stable `_rNN` report suffixes.
- Use batch-results actions (`Show`, `Export`, `Open`) for validation and downstream handoff.

## Concepts

- **Shared Engine Contract** (`shared_engine_contract`): GUI, CLI, shell, and scripting interfaces execute the same operation semantics.
- **Deterministic Workflows** (`deterministic_workflows`): Operation chains should produce stable IDs and comparable outputs across repeated runs.
- **Sequence Lineage** (`sequence_lineage`): Derived sequences are explicit products linked to upstream inputs and operations.

## At a Glance

1. Load test_files/tp73.ncbi.gb, then open Patterns -> PCR Designer... (or comma...
2. Keep the map in linear mode. Choose paint role ROI (green) and drag over one ...
3. Choose paint role Upstream (red) and drag one upstream primer-window interval...
4. Use post-drag actions (Set PCR ROI, Add ROI to Queue, Open PCR Designer) or S...
5. Add one annotation-driven queue row via DNA-window fallback (PCR ROI -> Add s...
6. In Design primer pairs, set shared constraints and a report_id base (for exam...
7. Optionally enable Also create extracted region copies, then run Design Primer...
8. In PCR batch results, use Show for report summary, Export for JSON output, an...

## GUI First

1. Load `test_files/tp73.ncbi.gb`, then open `Patterns -> PCR Designer...` (or command palette `PCR Designer`).
2. Keep the map in linear mode. Choose paint role `ROI` (green) and drag over one TP73-AS2 promoter-proximal window (for example around `61720..62120`).
3. Choose paint role `Upstream` (red) and drag one upstream primer-window interval. Choose `Downstream` (blue) and drag one downstream primer-window interval.
4. Use post-drag actions (`Set PCR ROI`, `Add ROI to Queue`, `Open PCR Designer`) or `Shift+drag` on ROI for immediate queue append.
5. Add one annotation-driven queue row via DNA-window fallback (`PCR ROI -> Add selected feature(s) to queue`) so queue sources include both painted and feature-driven entries.
6. In `Design primer pairs`, set shared constraints and a `report_id` base (for example `tp73_as2_promoter_batch`). The right pane live cartoon should reflect painted ROI/window geometry.
7. Optionally enable `Also create extracted region copies`, then run `Design Primer Pairs for queued regions`.
8. In `PCR batch results`, use `Show` for report summary, `Export` for JSON output, and `Open` to jump to the extracted copy (or template when copy mode is off).

## Command Equivalent (After GUI)

Run the same routine non-interactively once the GUI flow is clear:

```bash
cargo run --bin gentle_cli -- workflow @docs/examples/workflows/pcr_selection_batch_primer_pairs_offline.json
cargo run --bin gentle_cli -- shell 'workflow @docs/examples/workflows/pcr_selection_batch_primer_pairs_offline.json'
```

## Parameters That Matter

- `TP73-AS2 promoter-proximal ROI window (for example `61720..62120`)` (where used: painted ROI role and Design primer pairs ROI fields around TP73-AS2)
  - Why it matters: Keeps PCR targets tied to a biologically meaningful promoter-proximal context that can later feed luciferase planning.
  - How to derive it: Anchor on TP73-AS2 annotation in `test_files/tp73.ncbi.gb` (gene/ncRNA at lines ~1257-1270), then choose a local promoter-proximal window around the antisense 5' boundary.
- `Paint roles (`ROI` green, `Upstream` red, `Downstream` blue)` (where used: PCR Designer map interaction + live pair-PCR cartoon geometry)
  - Why it matters: Separates biological target interval from allowed primer-placement windows and makes constraint intent visible before design execution.
  - How to derive it: Paint ROI over the target biology first, then paint upstream/downstream windows that bound forward/reverse primer search space.
- `DesignPrimerPairs.report_id base` (where used: batch execution in Design primer pairs form)
  - Why it matters: Batch runs derive deterministic report IDs from this base with `_rNN` suffixes.
  - How to derive it: Pick one short base tied to experiment context (for example `tp73_as2_promoter_batch`).
- `PCR queue copy toggle (`Also create extracted region copies`)` (where used: batch execution in Design primer pairs form)
  - Why it matters: Enables per-region `ExtractRegion` artifacts before/alongside primer report generation.
  - How to derive it: Enable when downstream steps need sequence artifacts; leave off for report-only planning (default off).
- `Queue source mode (`painted ROI` vs `selected feature` fallback)` (where used: queue rows before batch execution)
  - Why it matters: Preserves coordinate provenance for each batch row and keeps mixed ad hoc + annotation flows inspectable.
  - How to derive it: Use paint-first for custom windows; use selected-feature fallback for annotation-defined intervals.

## Follow-up Commands

```bash
cargo run --bin gentle_cli -- shell 'ui open pcr-design'
cargo run --bin gentle_cli -- shell 'ui focus pcr-design'
cargo run --bin gentle_cli -- workflow @docs/examples/workflows/pcr_selection_batch_primer_pairs_offline.json
cargo run --bin gentle_cli -- shell 'primers list-reports'
cargo run --bin gentle_cli -- shell 'primers show-report tp73_as2_promoter_batch_r01'
cargo run --bin gentle_cli -- shell 'primers export-report tp73_as2_promoter_batch_r01 tp73_as2_promoter_batch_r01.json'
```

## Checkpoints

- PCR queue rows include both `painted ROI` and selected-feature source labels.
- Live pair-PCR cartoon preview updates when ROI/upstream/downstream paint intervals change.
- Two primer reports are created with deterministic `_r01`/`_r02` suffixes and appear in batch results.
- If your runtime falls back to the internal primer backend, reports may validly contain zero accepted pairs under strict constraints; report creation/export behavior remains deterministic.
- If copy toggle is enabled, each batch-results row includes a copy ID and `Open` jumps to that extracted sequence.
- At least one report can be shown and exported from the batch-results row actions without manual JSON editing.

## Canonical Source

- Chapter id: `pcr_selection_batch_primer_pairs_offline`
- Tier: `core`
- Example id: `pcr_selection_batch_primer_pairs_offline`
- Workflow file: `docs/examples/workflows/pcr_selection_batch_primer_pairs_offline.json`
- Example test_mode: `always`
- Executed during generation: `yes`
- Inspect this JSON file directly when you need full option-level detail.
