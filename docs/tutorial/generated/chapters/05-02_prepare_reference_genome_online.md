---
chapter_id: "prepare_reference_genome_online"
title: "Prepare a reference genome cache (online)"
tier: "online"
example_id: "prepare_reference_genome_online"
source_example: "docs/examples/workflows/prepare_reference_genome_online.json"
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
generated_artifact_dir: "docs/tutorial/generated/artifacts/prepare_reference_genome_online"
---

# Prepare a reference genome cache (online)

Document network-dependent preparation without destabilizing offline default checks.

Reference-genome preparation is crucial for genome-anchored cloning interpretation, but it depends on online resources and local tool setup. This chapter keeps that path explicit and opt-in so default tutorial checks remain robust.

## What You Will Accomplish

- Recognize which tutorial flows require network access and external tools.
- Use `GENTLE_TEST_ONLINE` to opt into online chapter execution.
- Preserve offline CI reliability while still documenting online capabilities.

## Before You Start

> **How to Run This Locally**
> Run `genomes status` first and read `effective_cache_dir`; with catalog `assets/genomes.json`, relative cache `data/genomes` currently resolves to `assets/data/genomes`. Only then set `GENTLE_TEST_ONLINE=1`. The workflow downloads the GRCh38 Ensembl 116 soft-masked FASTA and GTF from the standard `https://ftp.ensembl.org/pub/release-116/` FASTA/GTF trees; ensure the effective cache filesystem has enough space and interrupted downloads can be retried.

**Useful when:**

- You need genome-anchored extraction around gene/promoter context.
- You want to prepare local cache/index assets for repeated anchor operations.
- You need to understand which routines are intentionally online-only in CI defaults.

## Walkthrough: GUI, CLI and Inner Agent

Each step pairs GUI instructions with related terminal commands or guidance and a review-only inner-agent prompt. Some steps require GUI interaction; a listing command only inspects results, it does not perform the design. CLI snippets assume an installed `gentle_cli` and use GENtle's default `.gentle_state.json` unless stated otherwise. From a source checkout, replace `gentle_cli` with `cargo run --bin gentle_cli --`. Add `--state PATH` or `--project PATH` for an explicit sandbox; a separate CLI process does not inherit the open GUI project's unsaved state.

In the **GUI Shell**, enter only the shared command inside `gentle_cli shell '...'`, without the executable prefix or outer quotes. Run UI-opening commands there to open windows: a headless CLI returns the UI intent but does not open a GUI. Other terminal commands are not automatically GUI Shell commands. Inner-agent examples request a proposal for review; they are not executed during tutorial generation.

### Step 1: Open File -> Prepare Reference Genome..., select catalog assets/genomes.json and inspect the status for Human GRCh38 Ensembl 116, including both remote sources and the effective cache directory

**GUI**

Open `File -> Prepare Reference Genome...`, select catalog `assets/genomes.json` and inspect the status for `Human GRCh38 Ensembl 116`, including both remote sources and the effective cache directory.

**CLI (terminal)**

```bash
gentle_cli genomes status "Human GRCh38 Ensembl 116" --catalog assets/genomes.json --cache-dir data/genomes
```

**Ask the inner agent**

> Read `assets/genomes.json` and run only `genomes status` for `Human GRCh38 Ensembl 116` with cache `data/genomes`. Return the exact Ensembl-116 soft-masked FASTA and GTF URLs, lifecycle state and `effective_cache_dir`; explain that the relative cache is resolved against the catalog location. Do not download anything.

**Expected**

> The preflight status names `Human GRCh38 Ensembl 116`, the exact Ensembl-116 FASTA/GTF sources and the effective cache directory without downloading.

### Step 2: After explicit approval for network transfer, disk use and cache writes, prepare exactly Human GRCh38 Ensembl 116 with the reviewed catalog/cache settings

**GUI**

After explicit approval for network transfer, disk use and cache writes, prepare exactly `Human GRCh38 Ensembl 116` with the reviewed catalog/cache settings.

**CLI (terminal)**

```bash
GENTLE_TEST_ONLINE=1 gentle_cli workflow @docs/examples/workflows/prepare_reference_genome_online.json
```

**Ask the inner agent**

> If the status is missing, propose exactly one `PrepareGenome` operation for `Human GRCh38 Ensembl 116`, catalog `assets/genomes.json`, cache `data/genomes` and timeout 3600 seconds. Estimate the network/disk boundary from the catalog/status, expose all cache writes, and wait for explicit approval before network access or mutation; do not substitute Ensembl 113, NCBI RefSeq or another GRCh38 family member.

**Expected**

> Only the explicitly opted-in workflow prepares the selected cache target; offline tutorial generation leaves it skipped.

### Step 3: Query status again and confirm prepared component metadata before attempting extraction workflows

**GUI**

Query status again and confirm prepared component metadata before attempting extraction workflows.

**CLI (terminal)**

```bash
gentle_cli genomes status "Human GRCh38 Ensembl 116" --catalog assets/genomes.json --cache-dir data/genomes
```

**Ask the inner agent**

> After preparation completes, run the same read-only status request. Confirm `prepared=true`, the exact requested catalog key, effective cache, source URLs and component metadata before proposing any gene extraction; a successful transfer alone is not evidence that a later biological target is correct.

**Expected**

> Postflight status reports `prepared=true` and component metadata before extraction or promoter chapters depend on this reference.


## Ask an Outer Agent (MCP or ClawBio/OpenClaw)

An outer agent does not inherit the unsaved GUI project. Give it this chapter's canonical workflow and an explicit disposable state path; ask it to retain the structured result, artifacts and reproducibility receipt instead of replacing them with prose.

> Use GENtle's `gentle-cloning` skill to replay `docs/examples/workflows/prepare_reference_genome_online.json` against a new disposable state. First report the exact workflow, inputs, state path, outputs and whether the selected route needs confirmation. Do not infer state from an open GUI. Return the structured result, produced artifacts and reproducibility receipt, and state any unmet prerequisite.

Equivalent direct structured request for the generic wrapper:

```json
{
  "schema": "gentle.clawbio_skill_request.v1",
  "mode": "workflow",
  "state_path": "/tmp/gentle-prepare-reference-genome-online.state.json",
  "workflow_path": "docs/examples/workflows/prepare_reference_genome_online.json",
  "timeout_secs": 7200
}
```

Submitting a direct structured request is an explicit wrapper invocation. If natural language selects a narrower delegated skill and that route mutates state, selects biological material or writes artifacts, the caller must preserve that skill's proposal/approval boundary and approve only the exact bound digest. This tutorial generation step does not invoke an agent or grant approval.

## Interpretation and Reference

## Parameters That Matter

- `PrepareGenome.genome_id` (where used: operation 1)
  - Why it matters: Selects the exact reference build and annotation set.
  - How to derive it: Choose the genome build matching your experimental system and downstream coordinate system.
- `PrepareGenome.catalog_path / cache_dir` (where used: operation 1)
  - Why it matters: Controls source catalog and local cache destination.
  - How to derive it: Use repository defaults unless your environment requires custom catalogs or cache locations.

## Applied Concepts

- **Deterministic Workflows** (`deterministic_workflows`): Operation chains should produce stable IDs and comparable outputs across repeated runs.
- **Online Opt-in Execution** (`online_opt_in`): Network-dependent chapters remain explicit opt-in and do not break offline default CI.
- **Genome Catalog Targeting** (`genome_catalog_targeting`): Prepared genome catalogs, annotation-based gene filters, and anchor extension connect imported entries to genomic context.

## Follow-up Commands

```bash
GENTLE_TEST_ONLINE=1 cargo run --bin gentle_examples_docs -- tutorial-generate
```

## Checkpoints

- Preflight status resolves the requested Ensembl-116 sources and effective cache path without network mutation.
- Genome preparation runs only when `GENTLE_TEST_ONLINE` is enabled.
- Offline generation still emits the chapter as `skipped_online`; a later online run must prove `prepared=true` and component metadata separately.

## Tutorial Provenance

- Chapter id: `prepare_reference_genome_online`
- Tier: `online`
- Example id: `prepare_reference_genome_online`
- Tutorial source JSON: `docs/tutorial/sources/05-02_prepare_reference_genome_online.json`
- Workflow file: `docs/examples/workflows/prepare_reference_genome_online.json`
- Generated artifact dir: `docs/tutorial/generated/artifacts/prepare_reference_genome_online`
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

- Tutorial title: `Prepare a reference genome cache (online)`
- Tutorial/chapter id: `prepare_reference_genome_online`
- Step reached:
- Expected vs. actual:
- Interface used: GUI / CLI / Agent Assistant / ClawBio

Paste the Tutorial feedback context here:

```text

```
