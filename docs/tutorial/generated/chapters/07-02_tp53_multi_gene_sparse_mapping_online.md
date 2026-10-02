---
chapter_id: "tp53_multi_gene_sparse_mapping_online"
title: "Probe TP53-family sparse indexing from a TP53 locus (online)"
tier: "online"
example_id: "tp53_multi_gene_sparse_mapping_online"
source_example: "docs/examples/workflows/tp53_multi_gene_sparse_mapping_online.json"
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
generated_artifact_dir: "docs/tutorial/generated/artifacts/tp53_multi_gene_sparse_mapping_online"
---

# Probe TP53-family sparse indexing from a TP53 locus (online)

Prepare GRCh38, extract TP53, and inspect how a multi_gene_sparse request reports TP63/TP73 as unavailable in the local TP53-only annotation.

This chapter extends the TP53 genome-targeting path toward read-origin mapping in one deterministic route. It also documents an important current boundary: `multi_gene_sparse` expands transcript templates only from annotations attached to the selected DNA sequence; target-gene ids are not a request to download or join other loci. The canonical workflow extracts TP53 alone, so the reviewed Ensembl-116 run matches TP53, reports TP63 and TP73 as absent from local annotation, and adds no extra transcript lanes. This is therefore a reproducible warning/provenance exercise, not a completed TP53/TP63/TP73 family comparison.

The 2026-10-01 acceptance run reused a prepared Ensembl-116 cache, processed the committed 10,000-read fixture, seed-passed 11 reads, and aligned none because the workflow stops after phase 1. It also exposed and repaired a real parity defect: seed feature 0 was the gene row and produced a report that the documented mRNA workspace could not open. The workflow now binds feature 1, the first TP53 mRNA, so the native report view and CLI address the same saved report.

## What You Will Accomplish

- Map one TP53 locus from reference-genome preparation through phase-1 read interpretation in a single operation chain.
- Understand that `multi_gene_sparse` expands only from transcript annotations already attached to the selected sequence.
- Interpret deterministic report provenance for the seed feature, requested/matched/missing genes, phase-1 counts, and planned ROI capture flag.

## Before You Start

**Prerequisites:** Read [Chapter 9: Prepare a reference genome cache (online)](./05-02_prepare_reference_genome_online.md) first.

> **How to Run This Locally**
> Set `GENTLE_TEST_ONLINE=1` and run from the repository root. The default workflow may download and prepare roughly 1 GB of compressed GRCh38 sequence plus annotation unless the exact cache directory already contains a valid Ensembl-116 manifest. A relocated prepared cache must be supplied through an explicitly reviewed workflow/cache path; do not assume that matching files at a different path will be reused. The run indexes TP53-local transcripts only and reports TP63/TP73 missing.

**Useful when:**

- You want to verify that target-gene ids do not silently fetch or merge annotations from other loci.
- You need a reproducible TP53 baseline for comparing single-gene and sparse-origin request provenance.
- You want GUI and CLI routes to expose the same missing-local-annotation boundary before constructing a true multi-locus fixture.

## Walkthrough: GUI, CLI and Inner Agent

Each step pairs GUI instructions with related terminal commands or guidance and a review-only inner-agent prompt. Some steps require GUI interaction; a listing command only inspects results, it does not perform the design. CLI snippets assume an installed `gentle_cli` and use GENtle's default `.gentle_state.json` unless stated otherwise. From a source checkout, replace `gentle_cli` with `cargo run --bin gentle_cli --`. Add `--state PATH` or `--project PATH` for an explicit sandbox; a separate CLI process does not inherit the open GUI project's unsaved state.

In the **GUI Shell**, enter only the shared command inside `gentle_cli shell '...'`, without the executable prefix or outer quotes. Run UI-opening commands there to open windows: a headless CLI returns the UI intent but does not open a GUI. Other terminal commands are not automatically GUI Shell commands. Inner-agent examples request a proposal for review; they are not executed during tutorial generation.

### Step 1: Prepare or reuse the Ensembl-116 TP53 locus

**GUI**

Prepare or explicitly reuse `Human GRCh38 Ensembl 116`, verify the cache manifest, and extract `TP53` into `grch38_tp53`.

**CLI (terminal)**

```bash
GENTLE_TEST_ONLINE=1 gentle_cli workflow @docs/examples/workflows/tp53_multi_gene_sparse_mapping_online.json
```

**Ask the inner agent**

> Inspect the exact `Human GRCh38 Ensembl 116` catalog entry and requested cache path before proposing preparation. Report whether GENtle will reuse a checksum-verified manifest or perform network downloads; do not treat a similarly named cache at another path as reusable without verification.

**Expected**

> The workflow reuses or prepares the exact Ensembl-116 cache, extracts 41 TP53 features, and binds interpretation to mRNA feature `1` rather than the gene row.

![The native project graph retains the Ensembl-116 TP53 locus, its exon-concatenated derivative, and the 10,000-read sparse-origin report as separate provenance nodes. The report records 11 seed-passed reads and zero aligned reads because this workflow runs phase 1 only.](../../../screenshots/tp53_multi_gene_sparse_mapping_online/01-project-graph.png)

*Figure: The native project graph retains the Ensembl-116 TP53 locus, its exon-concatenated derivative, and the 10,000-read sparse-origin report as separate provenance nodes. The report records 11 seed-passed reads and zero aligned reads because this workflow runs phase 1 only. Screenshot captured 2026-10-01.*

### Step 2: Bind the report to TP53 mRNA feature 1

**GUI**

Open Splicing Expert for TP53 mRNA feature `1` (`ENST00000269305` in the reviewed Ensembl-116 extraction), set `Origin mode` to `multi_gene_sparse`, and set `Target genes` to `TP53, TP63, TP73`.

**CLI (terminal)**

```bash
gentle_cli shell 'rna-reads list-reports grch38_tp53'
```

**Ask the inner agent**

> After extraction, query the TP53 features and verify that feature `1` is mRNA `ENST00000269305` in this pinned release. Propose the interpretation request with `origin_mode=multi_gene_sparse` and requested targets TP53, TP63 and TP73; never substitute the gene feature at index `0` for the documented mRNA workspace.

**Expected**

> The report list includes `tp53_family_sparse_template` with `origin_mode=multi_gene_sparse`, three requested target ids, 10,000 total reads, 11 seed-passed reads, and zero aligned reads.

### Step 3: Inspect the sparse-origin boundary and phase-1 counts

**GUI**

Run Nanopore phase-1 interpretation, load report `tp53_family_sparse_template`, and inspect the counts and sparse-origin warnings. Confirm that TP53 matched locally, TP63/TP73 did not, no extra lanes were added, 11/10,000 reads seed-passed, and alignment remains deferred.

**CLI (terminal)**

```bash
gentle_cli shell 'rna-reads show-report tp53_family_sparse_template'
```

**Ask the inner agent**

> Inspect the saved structured report before making a biological statement. Keep requested genes separate from matched genes and missing local annotations; report 11 seed-passed of 10,000, zero aligned, and phase-1-only status. Do not call this a TP53-family comparison or infer TP63/TP73 absence from the reads.

**Expected**

> The report warns that no additional lanes were added, TP53 matched, and TP63/TP73 were not found in the selected sequence's local annotation.

![The native RNA-read Mapping workspace opens the report through TP53 mRNA feature 1 and shows that multi_gene_sparse added no local lanes beyond TP53. The structured report additionally records TP63 and TP73 as absent from the extracted sequence's local annotation.](../../../screenshots/tp53_multi_gene_sparse_mapping_online/02-sparse-report-warnings.png)

*Figure: The native RNA-read Mapping workspace opens the report through TP53 mRNA feature 1 and shows that multi_gene_sparse added no local lanes beyond TP53. The structured report additionally records TP63 and TP73 as absent from the extracted sequence's local annotation. Screenshot captured 2026-10-01.*


## Ask an Outer Agent (MCP or ClawBio/OpenClaw)

An outer agent does not inherit the unsaved GUI project. Give it this chapter's canonical workflow and an explicit disposable state path; ask it to retain the structured result, artifacts and reproducibility receipt instead of replacing them with prose.

> Use GENtle's `gentle-cloning` skill to replay `docs/examples/workflows/tp53_multi_gene_sparse_mapping_online.json` against a new disposable state. First report the exact workflow, inputs, state path, outputs and whether the selected route needs confirmation. Do not infer state from an open GUI. Return the structured result, produced artifacts and reproducibility receipt, and state any unmet prerequisite.

Equivalent direct structured request for the generic wrapper:

```json
{
  "schema": "gentle.clawbio_skill_request.v1",
  "mode": "workflow",
  "state_path": "/tmp/gentle-tp53-multi-gene-sparse-mapping-online.state.json",
  "workflow_path": "docs/examples/workflows/tp53_multi_gene_sparse_mapping_online.json",
  "timeout_secs": 7200
}
```

Submitting a direct structured request is an explicit wrapper invocation. If natural language selects a narrower delegated skill and that route mutates state, selects biological material or writes artifacts, the caller must preserve that skill's proposal/approval boundary and approve only the exact bound digest. This tutorial generation step does not invoke an agent or grant approval.

## Interpretation and Reference

## Parameters That Matter

- `InterpretRnaReads.seed_feature_id` (where used: operation 3)
  - Why it matters: Defines which feature seeds the splicing view and baseline transcript scope for sparse expansion.
  - How to derive it: Query the pinned extracted sequence. In the reviewed Ensembl-116 result feature 0 is the gene and feature 1 is mRNA ENST00000269305; the canonical workflow must use 1 so the native mRNA workspace can reopen the report.
- `InterpretRnaReads.origin_mode / target_gene_ids` (where used: operation 3)
  - Why it matters: Controls whether indexing stays baseline (`single_gene`) or is expanded with matched target-gene transcripts (`multi_gene_sparse`).
  - How to derive it: Use TP53, TP63 and TP73 here to exercise provenance and missing-local-annotation reporting. A true contrast requires those transcript annotations on the selected sequence; ids alone do not fetch or combine loci.
- `InterpretRnaReads.roi_seed_capture_enabled` (where used: operation 3)
  - Why it matters: Tracks request intent for future annotation-independent ROI capture layer.
  - How to derive it: Keep `false` for current runtime behavior; set `true` only when you want the deterministic pending-feature warning in provenance.

## Applied Concepts

- **Shared Engine Contract** (`shared_engine_contract`): GUI, CLI, shell, and scripting interfaces execute the same operation semantics.
- **Deterministic Workflows** (`deterministic_workflows`): Operation chains should produce stable IDs and comparable outputs across repeated runs.
- **Genome Catalog Targeting** (`genome_catalog_targeting`): Prepared genome catalogs, annotation-based gene filters, and anchor extension connect imported entries to genomic context.
- **Online Opt-in Execution** (`online_opt_in`): Network-dependent chapters remain explicit opt-in and do not break offline default CI.

## Follow-up Commands

```bash
gentle_cli workflow @docs/examples/workflows/tp53_multi_gene_sparse_mapping_online.json
gentle_cli shell 'rna-reads list-reports grch38_tp53'
gentle_cli shell 'rna-reads show-report tp53_family_sparse_template'
```

## Checkpoints

- TP53 locus extraction and RNA-read interpretation operations are both present in one workflow run history, and the report is bound to mRNA feature 1.
- Report summary shows `origin_mode=multi_gene_sparse`, requested targets TP53/TP63/TP73, 10,000 total reads, 11 seed-passed reads and zero aligned reads.
- Warnings explicitly state that TP53 matched, TP63/TP73 were absent from local annotation, and no additional transcript lanes were added.

## Tutorial Provenance

- Chapter id: `tp53_multi_gene_sparse_mapping_online`
- Tier: `online`
- Example id: `tp53_multi_gene_sparse_mapping_online`
- Tutorial source JSON: `docs/tutorial/sources/07-02_tp53_multi_gene_sparse_mapping_online.json`
- Workflow file: `docs/examples/workflows/tp53_multi_gene_sparse_mapping_online.json`
- Generated artifact dir: `docs/tutorial/generated/artifacts/tp53_multi_gene_sparse_mapping_online`
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

- Tutorial title: `Probe TP53-family sparse indexing from a TP53 locus (online)`
- Tutorial/chapter id: `tp53_multi_gene_sparse_mapping_online`
- Step reached:
- Expected vs. actual:
- Interface used: GUI / CLI / Agent Assistant / ClawBio

Paste the Tutorial feedback context here:

```text

```
