---
chapter_id: "guides_export_csv_and_protocol"
title: "Guide oligo export (CSV + protocol)"
tier: "advanced"
example_id: "guides_export_csv_and_protocol"
source_example: "docs/examples/workflows/guides_export_csv_and_protocol.json"
example_test_mode: "skip"
executed_during_generation: true
---

# Guide oligo export (CSV + protocol)

Export representative machine-readable and human-readable artifacts.

Operational work is only useful if outputs can be shared with collaborators and ordering pipelines. This routine keeps one representative CSV and one protocol text file so the tutorial remains readable while still proving export behavior.

**Prerequisites:** Read [Chapter 5: Guide practical filtering and oligo generation](./05_guides_filter_and_generate_oligos.md) first.

## Parameters That Matter

- `ExportGuideOligos.format` (where used: operation 3)
  - Why it matters: Format must match the receiving workflow (spreadsheet import, plate workflow, or FASTA).
  - How to derive it: Pick `csv_table` for review/order sheets, `plate_csv` for plate automation, `fasta` for sequence-oriented tools.
- `ExportGuideProtocolText.include_qc_checklist` (where used: operation 4)
  - Why it matters: Controls whether QC reminders are embedded in the generated protocol text.
  - How to derive it: Enable for handoff to wet-lab execution; disable only for compact machine-only summaries.

## When This Routine Is Useful

- You need oligo tables for ordering and a protocol summary for bench execution.
- You want reproducible artifacts tied to explicit operation history.
- You need a concise output bundle for review without committing redundant files.

## What You Learn

- Export guide outputs in machine-readable and human-readable forms.
- Understand selective artifact retention for tutorial readability.
- Map retained artifacts back to the operation chain that produced them.

## Applied Concepts

- **Guide Design Pipeline** (`guide_design_pipeline`): Guide sets can be created, filtered, expanded to oligos, and exported with protocol context.
- **Artifact Exports** (`artifact_exports`): Representative outputs (CSV/protocol/SVG/text) are retained for auditability and sharing.

## GUI First

CLI snippets use GENtle's default `.gentle_state.json` state unless they say otherwise. Add `--state PATH` or `--project PATH` when you want an explicit sandboxed state file for copied commands.

### Step 1: Open the guide/oligo export controls after guide generation

GUI: Open the guide/oligo export controls after guide generation.

CLI:

```bash
cargo run --bin gentle_cli -- guides put demo_guides --json '[{"guide_id":"demo_1","seq_id":"target_demo","start_0based":10,"end_0based_exclusive":30,"strand":"+","protospacer":"GACCTGTTGACGATGTTCCA","pam":"AGG","nuclease":"SpCas9","cut_offset_from_protospacer_start":17,"rank":1}]'
cargo run --bin gentle_cli -- guides oligos-generate demo_guides lenti_bsmbi_u6_default --apply-5prime-g-extension --output-oligo-set demo_lenti
```

> Expected: The guide set `demo_guides` and oligo set `demo_lenti` are stored with deterministic ids.

### Step 2: Export one CSV table and one protocol text file to verify both machine and hu...

GUI: Export one CSV table and one protocol text file to verify both machine and human output forms.

CLI:

```bash
cargo run --bin gentle_cli -- guides oligos-export demo_guides exports/demo_guides.csv --format csv_table --oligo-set demo_lenti
cargo run --bin gentle_cli -- guides protocol-export demo_guides exports/demo_guides.protocol.txt --oligo-set demo_lenti
```

> Expected: The CSV and protocol text files are written under `exports/` for machine and bench-facing review.

### Step 3: Inspect the exported files and confirm they match current guide/oligo set IDs

GUI: Inspect the exported files and confirm they match current guide/oligo set IDs.

CLI:

```bash
cargo run --bin gentle_cli -- guides oligos-show demo_lenti
```

> Expected: The oligo-set inspection output matches the exported set id and guide provenance.


## Checkpoints

- CSV export file exists and contains guide rows.
- Protocol text export exists and contains checklist content.

## What This Chapter Produces

- [`artifacts/guides_export_csv_and_protocol/exports/demo_guides.csv`](../artifacts/guides_export_csv_and_protocol/exports/demo_guides.csv) - `guide_id,rank,forward_oligo,reverse_oligo,notes`
- [`artifacts/guides_export_csv_and_protocol/exports/demo_guides.protocol.txt`](../artifacts/guides_export_csv_and_protocol/exports/demo_guides.protocol.txt) - `GENtle Guide Oligo Protocol`

## Canonical Source

- Chapter id: `guides_export_csv_and_protocol`
- Tier: `advanced`
- Example id: `guides_export_csv_and_protocol`
- Workflow file: `docs/examples/workflows/guides_export_csv_and_protocol.json`
- Example test_mode: `skip`
- Executed during generation: `yes`
- Inspect this JSON file directly when you need full option-level detail.
