# Guide oligo export (CSV + protocol)

- Chapter id: `guides_export_csv_and_protocol`
- Tier: `advanced`
- Example id: `guides_export_csv_and_protocol`
- Source example: `docs/examples/workflows/guides_export_csv_and_protocol.json`
- Example test_mode: `skip`
- Executed during generation: `yes`

Export representative machine-readable and human-readable artifacts.

Operational work is only useful if outputs can be shared with collaborators and ordering pipelines. This routine keeps one representative CSV and one protocol text file so the tutorial remains readable while still proving export behavior.

## When This Routine Is Useful

- You need oligo tables for ordering and a protocol summary for bench execution.
- You want reproducible artifacts tied to explicit operation history.
- You need a concise output bundle for review without committing redundant files.

## What You Learn

- Export guide outputs in machine-readable and human-readable forms.
- Understand selective artifact retention for tutorial readability.
- Map retained artifacts back to the operation chain that produced them.

## Concepts and Recurrence

- **Guide Design Pipeline** (`guide_design_pipeline`): Guide sets can be created, filtered, expanded to oligos, and exported with protocol context.
  - Status: reinforced from [Chapter 3: Guide practical filtering and oligo generation](./03_guides_filter_and_generate_oligos.md).
  - Reoccurs in: no later chapter.
- **Artifact Exports** (`artifact_exports`): Representative outputs (CSV/protocol/SVG/text) are retained for auditability and sharing.
  - Status: introduced in this chapter.
  - Reoccurs in: no later chapter.

## GUI First

1. Open the guide/oligo export controls after guide generation.
2. Export one CSV table and one protocol text file to verify both machine and human output forms.
3. Inspect the exported files and confirm they match current guide/oligo set IDs.

## Command Equivalent (After GUI)

Run the same routine non-interactively once the GUI flow is clear:

```bash
cargo run --bin gentle_cli -- workflow @docs/examples/workflows/guides_export_csv_and_protocol.json
cargo run --bin gentle_cli -- shell 'workflow @docs/examples/workflows/guides_export_csv_and_protocol.json'
```

## Parameters That Matter

- `ExportGuideOligos.format` (where used: operation 3)
  - Why it matters: Format must match the receiving workflow (spreadsheet import, plate workflow, or FASTA).
  - How to derive it: Pick `csv_table` for review/order sheets, `plate_csv` for plate automation, `fasta` for sequence-oriented tools.
- `ExportGuideProtocolText.include_qc_checklist` (where used: operation 4)
  - Why it matters: Controls whether QC reminders are embedded in the generated protocol text.
  - How to derive it: Enable for handoff to wet-lab execution; disable only for compact machine-only summaries.

## Checkpoints

- CSV export file exists and contains guide rows.
- Protocol text export exists and contains checklist content.

## Retained Outputs

- [`artifacts/guides_export_csv_and_protocol/exports/demo_guides.csv`](../artifacts/guides_export_csv_and_protocol/exports/demo_guides.csv)
- [`artifacts/guides_export_csv_and_protocol/exports/demo_guides.protocol.txt`](../artifacts/guides_export_csv_and_protocol/exports/demo_guides.protocol.txt)

## Canonical Source

- Workflow file: `docs/examples/workflows/guides_export_csv_and_protocol.json`
- Inspect this JSON file directly when you need full option-level detail.
