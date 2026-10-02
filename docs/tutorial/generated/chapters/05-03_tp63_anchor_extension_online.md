---
chapter_id: "tp63_anchor_extension_online"
title: "Retrieve TP63 and extend the displayed region by +/-2 kb (online)"
tier: "online"
example_id: "tp63_extend_anchor_online"
source_example: "docs/examples/workflows/tp63_extend_anchor_online.json"
example_test_mode: "online"
executed_during_generation: false
automated_status: "skipped_online"
review_status: "codex_reviewed"
review_stale: false
codex_reviewed_at: "2026-10-01"
human_reviewed_at: null
human_reviewer: null
review_stale_reason: null
review_issue_template: null
review_issue_template_path: null
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
> Complete chapter 05.02 first, then run from the repository root with the same catalog and effective cache. `GENTLE_TEST_ONLINE=1` permits the canonical all-in-one workflow to prepare a missing reference, so prefer the explicit status/list/extract/extend commands when reviewing approval boundaries. Current acceptance is offline/`skipped_online`: public TP63 metadata and local extension semantics were checked, but no large reference download or typed GUI capture was performed.

**Useful when:**

- You want to inspect promoter-proximal and downstream context around TP63 without manually typing coordinates.
- You need to confirm annotated TP63 coordinates before creating an anchored sequence.
- You want a reproducible +/-2 kb extension workflow that can be replayed by GUI and CLI users.

## Walkthrough: GUI, CLI and Inner Agent

Each step pairs GUI instructions with related terminal commands or guidance and a review-only inner-agent prompt. Some steps require GUI interaction; a listing command only inspects results, it does not perform the design. CLI snippets assume an installed `gentle_cli` and use GENtle's default `.gentle_state.json` unless stated otherwise. From a source checkout, replace `gentle_cli` with `cargo run --bin gentle_cli --`. Add `--state PATH` or `--project PATH` for an explicit sandbox; a separate CLI process does not inherit the open GUI project's unsaved state.

In the **GUI Shell**, enter only the shared command inside `gentle_cli shell '...'`, without the executable prefix or outer quotes. Run UI-opening commands there to open windows: a headless CLI returns the UI intent but does not open a GUI. Other terminal commands are not automatically GUI Shell commands. Inner-agent examples request a proposal for review; they are not executed during tutorial generation.

### Step 1: Confirm that chapter 05.02 completed for exactly Human GRCh38 Ensembl 116; in File -> Prepare Reference Genome..., inspect the prepared status and effective cache without starting another download

**GUI**

Confirm that chapter 05.02 completed for exactly `Human GRCh38 Ensembl 116`; in `File -> Prepare Reference Genome...`, inspect the prepared status and effective cache without starting another download.

**CLI (terminal)**

```bash
gentle_cli genomes status "Human GRCh38 Ensembl 116" --catalog assets/genomes.json --cache-dir data/genomes
```

**Ask the inner agent**

> Run only `genomes status` for exact catalog key `Human GRCh38 Ensembl 116` with catalog `assets/genomes.json` and cache `data/genomes`. Report `prepared`, lifecycle, effective cache and source identities. If it is not prepared, stop and point to chapter 05.02; do not start a download or silently use Ensembl 113/RefSeq.

**Expected**

> The exact Ensembl-116 reference is already prepared in the reported effective cache; this chapter does not hide preparation inside retrieval.

### Step 2: Open File -> Retrieve Genome Sequence..., filter to exact symbol TP63, and bind the candidate before extraction: ENSG00000073282, GRCh38 chromosome 3, plus strand, 189631252..189897276 (1-based inclusive), protein coding

**GUI**

Open `File -> Retrieve Genome Sequence...`, filter to exact symbol `TP63`, and bind the candidate before extraction: `ENSG00000073282`, GRCh38 chromosome 3, plus strand, 189631252..189897276 (1-based inclusive), protein coding.

**CLI (terminal)**

```bash
gentle_cli genomes genes "Human GRCh38 Ensembl 116" --catalog assets/genomes.json --cache-dir data/genomes --filter "^TP63$" --biotype protein_coding --limit 20
gentle_cli genomes extract-gene "Human GRCh38 Ensembl 116" TP63 --occurrence 1 --output-id grch38_tp63 --catalog assets/genomes.json --cache-dir data/genomes
```

**Ask the inner agent**

> On the same prepared reference, list exact-symbol, protein-coding TP63 candidates before proposing extraction. Bind the proposal to `ENSG00000073282`, chromosome 3, plus strand and 189631252..189897276 (1-based inclusive), occurrence 1 and output `grch38_tp63`; if the installed annotation disagrees, stop for review rather than selecting by row order. Execute only after approval.

**Expected**

> The listing binds `ENSG00000073282` on chromosome 3, plus strand, 189631252..189897276 before extraction creates anchored sequence `grch38_tp63`.

### Step 3: In the resulting DNA sequence window (grch38_tp63), use Extend 5' with 2000 bp, verify start 189629252, then use Extend 3' with 2000 bp and verify final anchor 3:189629252..189899276 with length 270025 bp

**GUI**

In the resulting DNA sequence window (`grch38_tp63`), use `Extend 5'` with `2000 bp`, verify start 189629252, then use `Extend 3'` with `2000 bp` and verify final anchor 3:189629252..189899276 with length 270025 bp.

**CLI (terminal)**

```bash
gentle_cli genomes extend-anchor grch38_tp63 5p 2000 --output-id grch38_tp63_ext5_2kb --catalog assets/genomes.json --cache-dir data/genomes
gentle_cli genomes extend-anchor grch38_tp63_ext5_2kb 3p 2000 --output-id grch38_tp63_ext5_ext3_2kb --catalog assets/genomes.json --cache-dir data/genomes
```

**Ask the inner agent**

> For approved `grch38_tp63`, propose two ordered `ExtendGenomeAnchor` operations: five-prime 2000 to `grch38_tp63_ext5_2kb`, then three-prime 2000 to `grch38_tp63_ext5_ext3_2kb`. Explain that plus-strand 5' lowers the genomic start and 3' raises the end. After approval, verify final anchor 3:189629252..189899276, length 270025 bp, parent chain and verified-anchor provenance; do not describe genomic left/right as 5'/3' without the strand.

**Expected**

> The five-prime result starts at 189629252; the final three-prime result ends at 189899276, has 270025 bp and preserves verified genome-anchor provenance.


## Ask an Outer Agent (MCP or ClawBio/OpenClaw)

An outer agent does not inherit the unsaved GUI project. Give it this chapter's canonical workflow and an explicit disposable state path; ask it to retain the structured result, artifacts and reproducibility receipt instead of replacing them with prose.

> Use GENtle's `gentle-cloning` skill to replay `docs/examples/workflows/tp63_extend_anchor_online.json` against a new disposable state. First report the exact workflow, inputs, state path, outputs and whether the selected route needs confirmation. Do not infer state from an open GUI. Return the structured result, produced artifacts and reproducibility receipt, and state any unmet prerequisite.

Equivalent direct structured request for the generic wrapper:

```json
{
  "schema": "gentle.clawbio_skill_request.v1",
  "mode": "workflow",
  "state_path": "/tmp/gentle-tp63-anchor-extension-online.state.json",
  "workflow_path": "docs/examples/workflows/tp63_extend_anchor_online.json",
  "timeout_secs": 7200
}
```

Submitting a direct structured request is an explicit wrapper invocation. If natural language selects a narrower delegated skill and that route mutates state, selects biological material or writes artifacts, the caller must preserve that skill's proposal/approval boundary and approve only the exact bound digest. This tutorial generation step does not invoke an agent or grant approval.

## Interpretation and Reference

## Parameters That Matter

- `ExtractGenomeGene.gene_query / occurrence` (where used: operation 2 and Retrieve Genome Sequence dialog)
  - Why it matters: Controls which TP63 transcript/gene match is selected when multiple annotation entries exist.
  - How to derive it: Require the reviewed Ensembl-116 gene row `ENSG00000073282`, chromosome 3, plus strand, 189631252..189897276. Treat `occurrence=1` as a replay selector only after those fields match; never use row order as biological identity.
- `ExtendGenomeAnchor.side / length_bp` (where used: operations 3 and 4 + DNA window Extend controls)
  - Why it matters: Defines biological flank direction and exact context length added around the anchored region.
  - How to derive it: Use `five_prime,2000` then `three_prime,2000`. For this plus-strand anchor, five-prime lowers start from 189631252 to 189629252 and three-prime raises end from 189897276 to 189899276; reverse-strand targets invert that genomic direction.

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

- The exact prepared Ensembl-116 cache and TP63 gene identity are visible before extraction; a missing or mismatched prerequisite stops the route.
- ExtractGenomeGene produces `grch38_tp63` anchored to plus-strand `ENSG00000073282` at 3:189631252..189897276.
- Sequential 5'/3' extension yields `grch38_tp63_ext5_ext3_2kb` at 3:189629252..189899276, length 270025 bp, with preserved verified-anchor provenance.

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
- Review status: `codex_reviewed`
- Codex reviewed at: `2026-10-01`
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
