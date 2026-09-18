---
chapter_id: "tp63_anchor_extension_online"
title: "Retrieve TP63 and extend the displayed region by +/-2 kb (online)"
tier: "online"
example_id: "tp63_extend_anchor_online"
source_example: "docs/examples/workflows/tp63_extend_anchor_online.json"
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
generated_artifact_dir: "docs/tutorial/generated/artifacts/tp63_anchor_extension_online"
---

# Retrieve TP63 and extend the displayed region by +/-2 kb (online)

Prepare Human GRCh38, inspect TP63 genomic coordinates, extract TP63, and extend the anchored sequence by 2000 bp on each side.

This chapter focuses on day-to-day genome-anchored sequence inspection: identify TP63 in GRCh38, verify coordinates before extraction, and then widen the visible anchored region directly in the DNA sequence window. The key point is that extension remains deterministic and provenance-preserving, whether triggered in the GUI or via shell/CLI operations.

## What You Will Accomplish

- Use annotation-backed gene retrieval to inspect coordinates before extracting a sequence.
- Apply anchored-region extension from the DNA sequence viewer without losing genome-anchor provenance.
- Map GUI actions to equivalent deterministic `ExtendGenomeAnchor` operations.

## Before You Start

**Prerequisites:** Read [Chapter 9: Prepare a reference genome cache (online)](./05-02_prepare_reference_genome_online.md) first.

> **How to Run This Locally**
> Set `GENTLE_TEST_ONLINE=1` and run from the repository root. The workflow prepares `Human GRCh38 Ensembl 116` using Ensembl FTP FASTA/GTF endpoints from `assets/genomes.json`, then extracts TP63 and extends the anchored locus locally.

**Useful when:**

- You want to inspect promoter-proximal and downstream context around TP63 without manually typing coordinates.
- You need to confirm annotated TP63 coordinates before creating an anchored sequence.
- You want a reproducible +/-2 kb extension workflow that can be replayed by GUI and CLI users.

## Walkthrough: GUI, CLI and Inner Agent

Each step pairs GUI instructions with related terminal commands or guidance and a review-only inner-agent prompt. Some steps require GUI interaction; a listing command only inspects results, it does not perform the design. CLI snippets assume an installed `gentle_cli` and use GENtle's default `.gentle_state.json` unless stated otherwise. From a source checkout, replace `gentle_cli` with `cargo run --bin gentle_cli --`. Add `--state PATH` or `--project PATH` for an explicit sandbox; a separate CLI process does not inherit the open GUI project's unsaved state.

In the **GUI Shell**, enter only the shared command inside `gentle_cli shell '...'`, without the executable prefix or outer quotes. Run UI-opening commands there to open windows: a headless CLI returns the UI intent but does not open a GUI. Other terminal commands are not automatically GUI Shell commands. Inner-agent examples request a proposal for review; they are not executed during tutorial generation.

### Step 1: Open File -> Prepare Reference Genome... and prepare Human GRCh38 Ensembl 116 from assets/genomes.json

**GUI**

Open `File -> Prepare Reference Genome...` and prepare `Human GRCh38 Ensembl 116` from `assets/genomes.json`.

**CLI (terminal)**

```bash
GENTLE_TEST_ONLINE=1 gentle_cli genomes prepare "Human GRCh38 Ensembl 116" --catalog assets/genomes.json --cache-dir data/genomes --timeout-secs 3600
```

**Ask the inner agent**

> In the current GENtle project, help me perform this tutorial step: Open `File -> Prepare Reference Genome...` and prepare `Human GRCh38 Ensembl 116` from `assets/genomes.json`. Show the exact GENtle operation or command and its expected result for my review. State any missing input. Do not execute it until I approve.

**Expected**

> The reference-preparation status becomes ready or reports a reusable prepared cache for GRCh38.

### Step 2: Open File -> Retrieve Genome Sequence..., set gene query/filter to TP63, review the displayed TP63 coordinate hit, and extract the first TP63 entry

**GUI**

Open `File -> Retrieve Genome Sequence...`, set gene query/filter to `TP63`, review the displayed TP63 coordinate hit, and extract the first TP63 entry.

**CLI (terminal)**

```bash
gentle_cli genomes genes "Human GRCh38 Ensembl 116" --catalog assets/genomes.json --cache-dir data/genomes --filter "^TP63$" --limit 20
gentle_cli genomes extract-gene "Human GRCh38 Ensembl 116" TP63 --occurrence 1 --output-id grch38_tp63 --catalog assets/genomes.json --cache-dir data/genomes
```

**Ask the inner agent**

> In the current GENtle project, help me perform this tutorial step: Open `File -> Retrieve Genome Sequence...`, set gene query/filter to `TP63`, review the displayed TP63 coordinate hit, and extract the first TP63 entry. Show the exact GENtle operation or command and its expected result for my review. State any missing input. Do not execute it until I approve.

**Expected**

> The gene listing shows the TP63 coordinate candidate(s), and extraction creates the anchored sequence id `grch38_tp63`.

### Step 3: In the resulting DNA sequence window (grch38_tp63), use Extend 5' with 2000 bp, then Extend 3' with 2000 bp to produce the +/-2 kb context sequence

**GUI**

In the resulting DNA sequence window (`grch38_tp63`), use `Extend 5'` with `2000 bp`, then `Extend 3'` with `2000 bp` to produce the +/-2 kb context sequence.

**CLI (terminal)**

```bash
gentle_cli genomes extend-anchor grch38_tp63 5p 2000 --output-id grch38_tp63_ext5_2kb --catalog assets/genomes.json --cache-dir data/genomes
gentle_cli genomes extend-anchor grch38_tp63_ext5_2kb 3p 2000 --output-id grch38_tp63_ext5_ext3_2kb --catalog assets/genomes.json --cache-dir data/genomes
```

**Ask the inner agent**

> In the current GENtle project, help me perform this tutorial step: In the resulting DNA sequence window (`grch38_tp63`), use `Extend 5'` with `2000 bp`, then `Extend 3'` with `2000 bp` to produce the +/-2 kb context sequence. Show the exact GENtle operation or command and its expected result for my review. State any missing input. Do not execute it until I approve.

**Expected**

> The final extension step creates or refreshes `grch38_tp63_ext5_ext3_2kb` while preserving genome-anchor provenance.


## Interpretation and Reference

## Parameters That Matter

- `ExtractGenomeGene.gene_query / occurrence` (where used: operation 2 and Retrieve Genome Sequence dialog)
  - Why it matters: Controls which TP63 transcript/gene match is selected when multiple annotation entries exist.
  - How to derive it: Start with `gene_query=TP63` and `occurrence=1`; inspect coordinate listing first and adjust occurrence only if you need a non-primary hit.
- `ExtendGenomeAnchor.side / length_bp` (where used: operations 3 and 4 + DNA window Extend controls)
  - Why it matters: Defines biological flank direction and exact context length added around the anchored region.
  - How to derive it: Use `five_prime,2000` then `three_prime,2000` for symmetric +/-2 kb context expansion.

## Applied Concepts

- **Shared Engine Contract** (`shared_engine_contract`): GUI, CLI, shell, and scripting interfaces execute the same operation semantics.
- **Genome Catalog Targeting** (`genome_catalog_targeting`): Prepared genome catalogs, annotation-based gene filters, and anchor extension connect imported entries to genomic context.
- **Sequence Lineage** (`sequence_lineage`): Derived sequences are explicit products linked to upstream inputs and operations.
- **Online Opt-in Execution** (`online_opt_in`): Network-dependent chapters remain explicit opt-in and do not break offline default CI.

## Follow-up Commands

```bash
gentle_cli genomes genes "Human GRCh38 Ensembl 116" --catalog assets/genomes.json --cache-dir data/genomes --filter "^TP63$" --limit 20
gentle_cli genomes extract-gene "Human GRCh38 Ensembl 116" TP63 --occurrence 1 --output-id grch38_tp63 --catalog assets/genomes.json --cache-dir data/genomes
gentle_cli genomes extend-anchor grch38_tp63 5p 2000 --output-id grch38_tp63_ext5_2kb --catalog assets/genomes.json --cache-dir data/genomes
gentle_cli genomes extend-anchor grch38_tp63_ext5_2kb 3p 2000 --output-id grch38_tp63_ext5_ext3_2kb --catalog assets/genomes.json --cache-dir data/genomes
```

## Checkpoints

- TP63 coordinate candidates are visible before extraction in Retrieve Genome Sequence workflows.
- ExtractGenomeGene produces `grch38_tp63` as anchored sequence context.
- Sequential 5'/3' extension by 2000 bp yields `grch38_tp63_ext5_ext3_2kb` with preserved anchor provenance.

## Tutorial Provenance

- Chapter id: `tp63_anchor_extension_online`
- Tier: `online`
- Example id: `tp63_extend_anchor_online`
- Tutorial source JSON: `docs/tutorial/sources/05-03_tp63_anchor_extension_online.json`
- Workflow file: `docs/examples/workflows/tp63_extend_anchor_online.json`
- Generated artifact dir: `docs/tutorial/generated/artifacts/tp63_anchor_extension_online`
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

- Tutorial title: `Retrieve TP63 and extend the displayed region by +/-2 kb (online)`
- Tutorial/chapter id: `tp63_anchor_extension_online`
- Step reached:
- Expected vs. actual:
- Interface used: GUI / CLI / Agent Assistant / ClawBio

Paste the Tutorial feedback context here:

```text

```
