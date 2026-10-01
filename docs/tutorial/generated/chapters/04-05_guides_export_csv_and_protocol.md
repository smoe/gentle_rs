---
chapter_id: "guides_export_csv_and_protocol"
title: "Guide oligo export (CSV + protocol)"
tier: "advanced"
example_id: "guides_export_csv_and_protocol"
source_example: "docs/examples/workflows/guides_export_csv_and_protocol.json"
example_test_mode: "skip"
executed_during_generation: true
automated_status: "passing"
review_status: "codex_reviewed"
review_stale: false
codex_reviewed_at: "2026-10-01"
human_reviewed_at: null
human_reviewer: null
review_stale_reason: null
review_issue_template: null
review_issue_template_path: null
generated_artifact_dir: "docs/tutorial/generated/artifacts/guides_export_csv_and_protocol"
---

# Guide oligo export (CSV + protocol)

Export one declared guide's template-formatted oligos as a retained CSV and a human-readable protocol draft, then compare both artifacts with their source set.

Operational work becomes reviewable when exact output bytes can be handed to collaborators. This routine declares one synthetic guide record, formats it with `lenti_bsmbi_u6_default`, and retains one CSV plus one protocol text file. The CSV proves deterministic tabular export and the text proves checklist-bearing protocol rendering. Neither file validates the declared target coordinates, authorizes ordering, or turns the generated suggestions into a laboratory-approved SOP.

## What You Will Accomplish

- Export guide outputs in machine-readable and human-readable forms without confusing export with ordering or bench authorization.
- Understand selective artifact retention for tutorial readability.
- Map retained artifacts back to the operation chain that produced them.

## Before You Start

**Prerequisites:** Read [Chapter 5: Filter declared guides and format cloning oligos](./04-04_guides_filter_and_generate_oligos.md) first.

**Useful when:**

- You need an oligo table and protocol draft for review before any ordering or bench authorization.
- You want reproducible artifacts tied to explicit operation history.
- You need a concise output bundle for review without committing redundant files.

## Walkthrough: GUI, CLI and Inner Agent

Each step pairs GUI instructions with related terminal commands or guidance and a review-only inner-agent prompt. Some steps require GUI interaction; a listing command only inspects results, it does not perform the design. CLI snippets assume an installed `gentle_cli` and use GENtle's default `.gentle_state.json` unless stated otherwise. From a source checkout, replace `gentle_cli` with `cargo run --bin gentle_cli --`. Add `--state PATH` or `--project PATH` for an explicit sandbox; a separate CLI process does not inherit the open GUI project's unsaved state.

In the **GUI Shell**, enter only the shared command inside `gentle_cli shell '...'`, without the executable prefix or outer quotes. Run UI-opening commands there to open windows: a headless CLI returns the UI intent but does not open a GUI. Other terminal commands are not automatically GUI Shell commands. Inner-agent examples request a proposal for review; they are not executed during tutorial generation.

### Step 1: Open the guide/oligo export controls after guide generation

**GUI**

Open the guide/oligo export controls after guide generation.

**CLI (terminal)**

```bash
gentle_cli guides put demo_guides --json '[{"guide_id":"demo_1","seq_id":"target_demo","start_0based":10,"end_0based_exclusive":30,"strand":"+","protospacer":"GACCTGTTGACGATGTTCCA","pam":"AGG","nuclease":"SpCas9","cut_offset_from_protospacer_start":17,"rank":1}]'
gentle_cli guides oligos-generate demo_guides lenti_bsmbi_u6_default --apply-5prime-g-extension --output-oligo-set demo_lenti
```

**Ask the inner agent**

> Propose `UpsertGuideSet` for the single declared record `demo_1` under `demo_guides`, followed by `GenerateGuideOligos` with template `lenti_bsmbi_u6_default`, 5-prime-G extension enabled and output id `demo_lenti`. Show the exact record and predicted forward/reverse oligos; state that no target reference is loaded, and wait for approval before mutating the registry.

**Expected**

> The guide set `demo_guides` and oligo set `demo_lenti` are stored with deterministic ids.

### Step 2: Export one CSV table and one protocol text file to verify both machine and human output forms

**GUI**

Export one CSV table and one protocol text file to verify both machine and human output forms.

**CLI (terminal)**

```bash
gentle_cli guides oligos-export demo_guides exports/demo_guides.csv --format csv_table --oligo-set demo_lenti
gentle_cli guides protocol-export demo_guides exports/demo_guides.protocol.txt --oligo-set demo_lenti
```

**Ask the inner agent**

> Propose two separate file writes: `ExportGuideOligos` in `csv_table` format to `exports/demo_guides.csv`, and `ExportGuideProtocolText` with the QC checklist to `exports/demo_guides.protocol.txt`. Show both resolved paths and source set ids before execution; do not treat approval to export as approval to order or run the protocol.

**Expected**

> The CSV and protocol text files are written under `exports/` for machine-readable and human-readable review.

### Step 3: Inspect the exported files and confirm they match current guide/oligo set IDs

**GUI**

Inspect the exported files and confirm they match current guide/oligo set IDs.

**CLI (terminal)**

```bash
gentle_cli guides oligos-show demo_lenti
```

**Ask the inner agent**

> Read the retained CSV and protocol text without modifying state. Confirm exactly one `demo_1` row, oligo set `demo_lenti`, template `lenti_bsmbi_u6_default`, matching forward/reverse sequences and the checklist. Return checksums or exact bytes for handoff, and label the text a protocol draft rather than a validated SOP.

**Expected**

> The oligo-set inspection output and both retained files match set id `demo_lenti`, guide `demo_1`, template `lenti_bsmbi_u6_default` and the exact oligo sequences.


## Ask an Outer Agent (MCP or ClawBio/OpenClaw)

An outer agent does not inherit the unsaved GUI project. Give it this chapter's canonical workflow and an explicit disposable state path; ask it to retain the structured result, artifacts and reproducibility receipt instead of replacing them with prose.

> Use GENtle's `gentle-cloning` skill to replay `docs/examples/workflows/guides_export_csv_and_protocol.json` against a new disposable state. First report the exact workflow, inputs, state path, outputs and whether the selected route needs confirmation. Do not infer state from an open GUI. Return the structured result, produced artifacts and reproducibility receipt, and state any unmet prerequisite.

Equivalent direct structured request for the generic wrapper:

```json
{
  "schema": "gentle.clawbio_skill_request.v1",
  "mode": "workflow",
  "state_path": "/tmp/gentle-guides-export-csv-and-protocol.state.json",
  "workflow_path": "docs/examples/workflows/guides_export_csv_and_protocol.json",
  "timeout_secs": 300
}
```

Submitting a direct structured request is an explicit wrapper invocation. If natural language selects a narrower delegated skill and that route mutates state, selects biological material or writes artifacts, the caller must preserve that skill's proposal/approval boundary and approve only the exact bound digest. This tutorial generation step does not invoke an agent or grant approval.

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

- CSV export contains exactly one `demo_1` row with the generated forward and reverse oligos.
- Protocol text names `demo_guides`, `demo_lenti` and `lenti_bsmbi_u6_default`, repeats the same oligos and contains checklist content.
- The chapter labels both files as review artifacts rather than target validation, ordering approval or a laboratory-approved SOP.

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
- Review status: `codex_reviewed`
- Codex reviewed at: `2026-10-01`
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
