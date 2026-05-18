---
chapter_id: "tp53_isoform_architecture_online"
title: "TP53 isoform architecture expert panel (online)"
tier: "online"
example_id: "tp53_isoform_architecture_online"
source_example: "docs/examples/workflows/tp53_isoform_architecture_online.json"
example_test_mode: "online"
executed_during_generation: false
---

# TP53 isoform architecture expert panel (online)

Build a TP53 locus project and export Figure-1-style transcript/protein isoform architecture from one deterministic engine route.

This chapter demonstrates a publication-oriented use case: derive TP53 from a prepared GRCh38 reference, import a curated isoform panel, and render transcript/protein architecture from the same expert-view payload used by GUI and shell interfaces. The focus is parity and provenance, not manual figure drawing.

> **How to Run This Locally**
> Set `GENTLE_TEST_ONLINE=1` and run from the repository root. The workflow prepares/extracts `Human GRCh38 Ensembl 116` from Ensembl FTP, then imports the local curated panel `assets/panels/tp53_isoforms_v1.json` and writes `exports/tp53_isoform_architecture.svg`.

## Parameters That Matter

- `ImportIsoformPanel.panel_path / panel_id / strict` (where used: operation 3)
  - Why it matters: Defines which curated panel is loaded and whether transcript mapping mismatches should fail hard.
  - How to derive it: Use `assets/panels/tp53_isoforms_v1.json` and keep `strict=false` for exploratory mapping; switch to `strict=true` when curation is locked.
- `RenderIsoformArchitectureSvg.path` (where used: operation 4)
  - Why it matters: Controls deterministic export location used for tutorial artifact retention and figure review.
  - How to derive it: Use a stable project-relative path such as `exports/tp53_isoform_architecture.svg`.

## When This Routine Is Useful

- You want transcript and protein isoform architecture from one sequence context with explicit panel provenance.
- You need a deterministic SVG baseline for Figure-1-style TP53 isoform presentation.
- You want to preserve the exact panel import and rendering operations in sequence history.

## What You Learn

- Execute an end-to-end isoform-panel workflow using shared engine operations.
- Understand how curated panel JSON joins genome-derived transcripts without GUI-only logic.
- Export deterministic architecture SVG suitable for downstream figure styling.

## Applied Concepts

- **Shared Engine Contract** (`shared_engine_contract`): GUI, CLI, shell, and scripting interfaces execute the same operation semantics.
- **Genome Catalog Targeting** (`genome_catalog_targeting`): Prepared genome catalogs, annotation-based gene filters, and anchor extension connect imported entries to genomic context.
- **Isoform Architecture Panels** (`isoform_architecture_panels`): Curated transcript/protein architecture overlays can be imported and rendered as deterministic expert-view SVG outputs.
- **Online Opt-in Execution** (`online_opt_in`): Network-dependent chapters remain explicit opt-in and do not break offline default CI.
- **Artifact Exports** (`artifact_exports`): Representative outputs (CSV/protocol/SVG/text) are retained for auditability and sharing.

## GUI First

1. Prepare `Human GRCh38 Ensembl 116` and extract gene `TP53` into `grch38_tp53`.
2. Open DNA window `Engine Ops -> Isoform architecture panels`, import `assets/panels/tp53_isoforms_v1.json`.
3. Open `Isoform Expert` and export SVG from the same panel context.

## Command Equivalent (After GUI)

Run the same routine non-interactively once the GUI flow is clear:

```bash
cargo run --bin gentle_cli -- workflow @docs/examples/workflows/tp53_isoform_architecture_online.json
cargo run --bin gentle_cli -- shell 'workflow @docs/examples/workflows/tp53_isoform_architecture_online.json'
```

## Follow-up Commands

```bash
cargo run --bin gentle_cli -- shell 'panels inspect-isoform grch38_tp53 tp53_isoforms_v1'
cargo run --bin gentle_cli -- save-project tp53_isoform_architecture.project.gentle.json
```

## Checkpoints

- Panel import reports mapped transcript lanes and any unresolved transcript warnings.
- Isoform architecture SVG export succeeds with deterministic lane ordering.
- Saved project contains TP53 sequence plus imported panel metadata for replay.

## Canonical Source

- Chapter id: `tp53_isoform_architecture_online`
- Tier: `online`
- Example id: `tp53_isoform_architecture_online`
- Workflow file: `docs/examples/workflows/tp53_isoform_architecture_online.json`
- Example test_mode: `online`
- Executed during generation: `no`
- Execution note: set `GENTLE_TEST_ONLINE=1` before `tutorial-generate` to execute this chapter.
- Inspect this JSON file directly when you need full option-level detail.
