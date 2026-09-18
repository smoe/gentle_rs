---
chapter_id: "tp53_isoform_architecture_online"
title: "TP53 isoform architecture expert panel (online)"
tier: "online"
example_id: "tp53_isoform_architecture_online"
source_example: "docs/examples/workflows/tp53_isoform_architecture_online.json"
example_test_mode: "online"
executed_during_generation: false
automated_status: "skipped_online"
review_status: "unreviewed"
review_stale: false
codex_reviewed_at: null
human_reviewed_at: null
human_reviewer: null
review_stale_reason: null
review_issue_template: "Tutorial confusion"
review_issue_template_path: ".github/ISSUE_TEMPLATE/tutorial-confusion.md"
generated_artifact_dir: "docs/tutorial/generated/artifacts/tp53_isoform_architecture_online"
---

# TP53 isoform architecture expert panel (online)

Build a TP53 locus project and export Figure-1-style transcript/protein isoform architecture from one deterministic engine route.

This chapter demonstrates a publication-oriented use case: derive TP53 from a prepared GRCh38 reference, import a curated isoform panel, and render transcript/protein architecture from the same expert-view payload used by GUI and shell interfaces. The focus is parity and provenance, not manual figure drawing.

## What You Will Accomplish

- Execute an end-to-end isoform-panel workflow using shared engine operations.
- Understand how curated panel JSON joins genome-derived transcripts without GUI-only logic.
- Export deterministic architecture SVG suitable for downstream figure styling.

## Before You Start

**Prerequisites:** Read [Chapter 9: Prepare a reference genome cache (online)](./05-02_prepare_reference_genome_online.md) first.

> **How to Run This Locally**
> Set `GENTLE_TEST_ONLINE=1` and run from the repository root. The workflow prepares/extracts `Human GRCh38 Ensembl 116` from Ensembl FTP, then imports the local curated panel `assets/panels/tp53_isoforms_v1.json` and writes `exports/tp53_isoform_architecture.svg`.

**Useful when:**

- You want transcript and protein isoform architecture from one sequence context with explicit panel provenance.
- You need a deterministic SVG baseline for Figure-1-style TP53 isoform presentation.
- You want to preserve the exact panel import and rendering operations in sequence history.

## Walkthrough: GUI, CLI and Inner Agent

Each step pairs GUI instructions with related terminal commands or guidance and a review-only inner-agent prompt. Some steps require GUI interaction; a listing command only inspects results, it does not perform the design. CLI snippets assume an installed `gentle_cli` and use GENtle's default `.gentle_state.json` unless stated otherwise. From a source checkout, replace `gentle_cli` with `cargo run --bin gentle_cli --`. Add `--state PATH` or `--project PATH` for an explicit sandbox; a separate CLI process does not inherit the open GUI project's unsaved state.

In the **GUI Shell**, enter only the shared command inside `gentle_cli shell '...'`, without the executable prefix or outer quotes. Run UI-opening commands there to open windows: a headless CLI returns the UI intent but does not open a GUI. Other terminal commands are not automatically GUI Shell commands. Inner-agent examples request a proposal for review; they are not executed during tutorial generation.

### Step 1: Prepare Human GRCh38 Ensembl 116 and extract gene TP53 into grch38_tp53

**GUI**

Prepare `Human GRCh38 Ensembl 116` and extract gene `TP53` into `grch38_tp53`.

**CLI (terminal)**

```bash
GENTLE_TEST_ONLINE=1 gentle_cli genomes prepare "Human GRCh38 Ensembl 116" --catalog assets/genomes.json --cache-dir data/genomes --timeout-secs 3600
gentle_cli genomes extract-gene "Human GRCh38 Ensembl 116" TP53 --occurrence 1 --output-id grch38_tp53 --catalog assets/genomes.json --cache-dir data/genomes
```

**Ask the inner agent**

> In the current GENtle project, help me perform this tutorial step: Prepare `Human GRCh38 Ensembl 116` and extract gene `TP53` into `grch38_tp53`. Show the exact GENtle operation or command and its expected result for my review. State any missing input. Do not execute it until I approve.

**Expected**

> The reference is prepared if needed and TP53 is extracted into the stable anchored sequence id `grch38_tp53`.

### Step 2: Open DNA window Engine Ops -> Isoform architecture panels, import assets/panels/tp53_isoforms_v1.json

**GUI**

Open DNA window `Engine Ops -> Isoform architecture panels`, import `assets/panels/tp53_isoforms_v1.json`.

**CLI (terminal)**

```bash
gentle_cli shell 'panels import-isoform grch38_tp53 assets/panels/tp53_isoforms_v1.json --panel-id tp53_isoforms_v1'
gentle_cli shell 'panels inspect-isoform grch38_tp53 tp53_isoforms_v1'
```

**Ask the inner agent**

> In the current GENtle project, help me perform this tutorial step: Open DNA window `Engine Ops -> Isoform architecture panels`, import `assets/panels/tp53_isoforms_v1.json`. Show the exact GENtle operation or command and its expected result for my review. State any missing input. Do not execute it until I approve.

**Expected**

> The panel import stores `tp53_isoforms_v1`; the inspect route returns the same isoform-architecture payload that the GUI expert view renders.

### Step 3: Open Isoform Expert and export SVG from the same panel context

**GUI**

Open `Isoform Expert` and export SVG from the same panel context.

**CLI (terminal)**

```bash
gentle_cli shell 'panels render-isoform-svg grch38_tp53 tp53_isoforms_v1 exports/tp53_isoform_architecture.svg'
```

**Ask the inner agent**

> In the current GENtle project, help me perform this tutorial step: Open `Isoform Expert` and export SVG from the same panel context. Show the exact GENtle operation or command and its expected result for my review. State any missing input. Do not execute it until I approve.

**Expected**

> The renderer writes `exports/tp53_isoform_architecture.svg` with deterministic isoform lane ordering.


## Interpretation and Reference

## Parameters That Matter

- `ImportIsoformPanel.panel_path / panel_id / strict` (where used: operation 3)
  - Why it matters: Defines which curated panel is loaded and whether transcript mapping mismatches should fail hard.
  - How to derive it: Use `assets/panels/tp53_isoforms_v1.json` and keep `strict=false` for exploratory mapping; switch to `strict=true` when curation is locked.
- `RenderIsoformArchitectureSvg.path` (where used: operation 4)
  - Why it matters: Controls deterministic export location used for tutorial artifact retention and figure review.
  - How to derive it: Use a stable project-relative path such as `exports/tp53_isoform_architecture.svg`.

## Applied Concepts

- **Shared Engine Contract** (`shared_engine_contract`): GUI, CLI, shell, and scripting interfaces execute the same operation semantics.
- **Genome Catalog Targeting** (`genome_catalog_targeting`): Prepared genome catalogs, annotation-based gene filters, and anchor extension connect imported entries to genomic context.
- **Isoform Architecture Panels** (`isoform_architecture_panels`): Curated transcript/protein architecture overlays can be imported and rendered as deterministic expert-view SVG outputs.
- **Online Opt-in Execution** (`online_opt_in`): Network-dependent chapters remain explicit opt-in and do not break offline default CI.
- **Artifact Exports** (`artifact_exports`): Representative outputs (CSV/protocol/SVG/text) are retained for auditability and sharing.

## Follow-up Commands

```bash
gentle_cli shell 'panels inspect-isoform grch38_tp53 tp53_isoforms_v1'
gentle_cli save-project tp53_isoform_architecture.project.gentle.json
```

## Checkpoints

- Panel import reports mapped transcript lanes and any unresolved transcript warnings.
- Isoform architecture SVG export succeeds with deterministic lane ordering.
- Saved project contains TP53 sequence plus imported panel metadata for replay.

## Tutorial Provenance

- Chapter id: `tp53_isoform_architecture_online`
- Tier: `online`
- Example id: `tp53_isoform_architecture_online`
- Tutorial source JSON: `docs/tutorial/sources/06-03_tp53_isoform_architecture_online.json`
- Workflow file: `docs/examples/workflows/tp53_isoform_architecture_online.json`
- Generated artifact dir: `docs/tutorial/generated/artifacts/tp53_isoform_architecture_online`
- Example test_mode: `online`
- Executed during generation: `no`
- Automated status: `skipped_online`
- Review status: `unreviewed`
- Codex reviewed at: `not recorded`
- Human reviewed at: `not recorded`
- Execution note: set `GENTLE_TEST_ONLINE=1` before `tutorial-generate` to execute this chapter.
- Inspect the source JSON when you need full option-level detail.

## Feedback

If this tutorial is confusing, execution-stale, biologically suspect, or missing a useful figure, please open the matching tutorial issue template and include the context below.

- Tutorial title: `TP53 isoform architecture expert panel (online)`
- Tutorial/chapter id: `tp53_isoform_architecture_online`
- Step reached:
- Expected vs. actual:
- Interface used: GUI / CLI / Agent Assistant / ClawBio

Paste the Tutorial feedback context here:

```text

```
