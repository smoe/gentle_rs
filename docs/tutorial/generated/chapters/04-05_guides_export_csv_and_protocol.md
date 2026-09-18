---
chapter_id: "guides_export_csv_and_protocol"
title: "Guide oligo export (CSV + protocol)"
tier: "advanced"
example_id: "guides_export_csv_and_protocol"
source_example: "docs/examples/workflows/guides_export_csv_and_protocol.json"
example_test_mode: "skip"
executed_during_generation: true
automated_status: "passing"
review_status: "unreviewed"
review_stale: false
codex_reviewed_at: null
human_reviewed_at: null
human_reviewer: null
review_stale_reason: null
review_issue_template: "Tutorial confusion"
review_issue_template_path: ".github/ISSUE_TEMPLATE/tutorial-confusion.md"
generated_artifact_dir: "docs/tutorial/generated/artifacts/guides_export_csv_and_protocol"
---

# Guide oligo export (CSV + protocol)

Export representative machine-readable and human-readable artifacts.

Operational work is only useful if outputs can be shared with collaborators and ordering pipelines. This routine keeps one representative CSV and one protocol text file so the tutorial remains readable while still proving export behavior.

## What You Will Accomplish

- Export guide outputs in machine-readable and human-readable forms.
- Understand selective artifact retention for tutorial readability.
- Map retained artifacts back to the operation chain that produced them.

## Before You Start

**Prerequisites:** Read [Chapter 5: Guide practical filtering and oligo generation](./04-04_guides_filter_and_generate_oligos.md) first.

**Useful when:**

- You need oligo tables for ordering and a protocol summary for bench execution.
- You want reproducible artifacts tied to explicit operation history.
- You need a concise output bundle for review without committing redundant files.

## Walkthrough: GUI, CLI and Inner Agent

The three routes below describe the same operation. CLI snippets assume an installed `gentle_cli` and use GENtle's default `.gentle_state.json` unless stated otherwise. From a source checkout, replace `gentle_cli` with `cargo run --bin gentle_cli --`. Add `--state PATH` or `--project PATH` for an explicit sandbox. Inner-agent examples request a proposal for review; they are not executed during tutorial generation.

### Step 1: Open the guide/oligo export controls after guide generation

**GUI**

Open the guide/oligo export controls after guide generation.

**CLI / GUI Shell**

```bash
gentle_cli guides put demo_guides --json '[{"guide_id":"demo_1","seq_id":"target_demo","start_0based":10,"end_0based_exclusive":30,"strand":"+","protospacer":"GACCTGTTGACGATGTTCCA","pam":"AGG","nuclease":"SpCas9","cut_offset_from_protospacer_start":17,"rank":1}]'
gentle_cli guides oligos-generate demo_guides lenti_bsmbi_u6_default --apply-5prime-g-extension --output-oligo-set demo_lenti
```

**Ask the inner agent**

> In the current GENtle project, help me perform this tutorial step: Open the guide/oligo export controls after guide generation. Show the exact GENtle operation or command and its expected result for my review. State any missing input. Do not execute it until I approve.

**Expected**

> The guide set `demo_guides` and oligo set `demo_lenti` are stored with deterministic ids.

### Step 2: Export one CSV table and one protocol text file to

**GUI**

Export one CSV table and one protocol text file to verify both machine and human output forms.

**CLI / GUI Shell**

```bash
gentle_cli guides oligos-export demo_guides exports/demo_guides.csv --format csv_table --oligo-set demo_lenti
gentle_cli guides protocol-export demo_guides exports/demo_guides.protocol.txt --oligo-set demo_lenti
```

**Ask the inner agent**

> In the current GENtle project, help me perform this tutorial step: Export one CSV table and one protocol text file to verify both machine and human output forms. Show the exact GENtle operation or command and its expected result for my review. State any missing input. Do not execute it until I approve.

**Expected**

> The CSV and protocol text files are written under `exports/` for machine and bench-facing review.

### Step 3: Inspect the exported files and confirm they match current guide/oligo

**GUI**

Inspect the exported files and confirm they match current guide/oligo set IDs.

**CLI / GUI Shell**

```bash
gentle_cli guides oligos-show demo_lenti
```

**Ask the inner agent**

> In the current GENtle project, help me perform this tutorial step: Inspect the exported files and confirm they match current guide/oligo set IDs. Show the exact GENtle operation or command and its expected result for my review. State any missing input. Do not execute it until I approve.

**Expected**

> The oligo-set inspection output matches the exported set id and guide provenance.


## Interpretation and Reference

## Parameters That Matter

- `ExportGuideOligos.format` (where used: operation 3)
  - Why it matters: Format must match the receiving workflow (spreadsheet import, plate workflow, or FASTA).
  - How to derive it: Pick `csv_table` for review/order sheets, `plate_csv` for plate automation, `fasta` for sequence-oriented tools.
- `ExportGuideProtocolText.include_qc_checklist` (where used: operation 4)
  - Why it matters: Controls whether QC reminders are embedded in the generated protocol text.
  - How to derive it: Enable for handoff to wet-lab execution; disable only for compact machine-only summaries.

## Applied Concepts

- **Guide Design Pipeline** (`guide_design_pipeline`): Guide sets can be created, filtered, expanded to oligos, and exported with protocol context.
- **Artifact Exports** (`artifact_exports`): Representative outputs (CSV/protocol/SVG/text) are retained for auditability and sharing.

## Checkpoints

- CSV export file exists and contains guide rows.
- Protocol text export exists and contains checklist content.

## What This Chapter Produces

- [`artifacts/guides_export_csv_and_protocol/exports/demo_guides.csv`](../artifacts/guides_export_csv_and_protocol/exports/demo_guides.csv) - `guide_id,rank,forward_oligo,reverse_oligo,notes`
- [`artifacts/guides_export_csv_and_protocol/exports/demo_guides.protocol.txt`](../artifacts/guides_export_csv_and_protocol/exports/demo_guides.protocol.txt) - `GENtle Guide Oligo Protocol`

## Tutorial Provenance

- Chapter id: `guides_export_csv_and_protocol`
- Tier: `advanced`
- Example id: `guides_export_csv_and_protocol`
- Tutorial source JSON: `docs/tutorial/sources/04-05_guides_export_csv_and_protocol.json`
- Workflow file: `docs/examples/workflows/guides_export_csv_and_protocol.json`
- Generated artifact dir: `docs/tutorial/generated/artifacts/guides_export_csv_and_protocol`
- Example test_mode: `skip`
- Executed during generation: `yes`
- Automated status: `passing`
- Review status: `unreviewed`
- Codex reviewed at: `not recorded`
- Human reviewed at: `not recorded`
- Inspect the source JSON when you need full option-level detail.

## Feedback

If this tutorial is confusing, execution-stale, biologically suspect, or missing a useful figure, please open the matching tutorial issue template and include the context below.

- Tutorial title: `Guide oligo export (CSV + protocol)`
- Tutorial/chapter id: `guides_export_csv_and_protocol`
- Step reached:
- Expected vs. actual:
- Interface used: GUI / CLI / Agent Assistant / ClawBio

Paste the Tutorial feedback context here:

```text

```
