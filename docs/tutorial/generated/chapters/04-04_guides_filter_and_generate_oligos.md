---
chapter_id: "guides_filter_and_generate_oligos"
title: "Filter declared guides and format cloning oligos"
tier: "core"
example_id: "guides_filter_and_generate_oligos"
source_example: "docs/examples/workflows/guides_filter_and_generate_oligos.json"
example_test_mode: "always"
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
generated_artifact_dir: "docs/tutorial/generated/artifacts/guides_filter_and_generate_oligos"
---

# Filter declared guides and format cloning oligos

Apply explicit practical filters to three declared guide candidates and format only the passing guides with one named cloning template.

For CRISPR-style cloning, practical sequence filtering can prevent avoidable downstream failures, but it is only one review layer. This offline registry demonstration declares three guide records without loading or validating a TP73 reference sequence. It then records why `g2` fails the U6 `TTTT` rule, retains `g1` and `g3`, and formats those two candidates with `lenti_bsmbi_u6_default`. The resulting oligos are deterministic template-formatted candidates, not proof of target identity, genome-wide specificity, nuclease activity, cloning success or order readiness.

## What You Will Accomplish

- Create and filter guide sets through explicit engine operations.
- Generate template-formatted oligo records from filtered guide candidates without confusing formatting with target validation.
- Connect guide workflows to the same deterministic operation model used for cloning steps.

## Before You Start

**Prerequisites:** Read [Chapter 1: Load FASTA, branch, and reverse-complement](./02-01_load_branch_reverse_complement_pgex_fasta.md) first.

**Useful when:**

- You need a transparent guide filtering step before oligo ordering.
- You want to document why specific guides were excluded.
- You need deterministic oligo generation from a reusable guide set.

## Walkthrough: GUI, CLI and Inner Agent

Each step pairs GUI instructions with related terminal commands or guidance and a review-only inner-agent prompt. Some steps require GUI interaction; a listing command only inspects results, it does not perform the design. CLI snippets assume an installed `gentle_cli` and use GENtle's default `.gentle_state.json` unless stated otherwise. From a source checkout, replace `gentle_cli` with `cargo run --bin gentle_cli --`. Add `--state PATH` or `--project PATH` for an explicit sandbox; a separate CLI process does not inherit the open GUI project's unsaved state.

In the **GUI Shell**, enter only the shared command inside `gentle_cli shell '...'`, without the executable prefix or outer quotes. Run UI-opening commands there to open windows: a headless CLI returns the UI intent but does not open a GUI. Other terminal commands are not automatically GUI Shell commands. Inner-agent examples request a proposal for review; they are not executed during tutorial generation.

### Step 1: Open the guides workflow controls in GENtle and create/import a guide set for a target region

**GUI**

Open the guides workflow controls in GENtle and create/import a guide set for a target region.

**CLI (terminal)**

```bash
gentle_cli guides put tp73_guides --json '[{"guide_id":"g1","seq_id":"tp73","start_0based":100,"end_0based_exclusive":120,"strand":"+","protospacer":"GACCTGTTGACGATGTTCCA","pam":"AGG","nuclease":"SpCas9","cut_offset_from_protospacer_start":17,"rank":1},{"guide_id":"g2","seq_id":"tp73","start_0based":220,"end_0based_exclusive":240,"strand":"+","protospacer":"TTTTGCCATGTTGACCTGAA","pam":"TGG","nuclease":"SpCas9","cut_offset_from_protospacer_start":17,"rank":2},{"guide_id":"g3","seq_id":"tp73","start_0based":340,"end_0based_exclusive":360,"strand":"-","protospacer":"GGTACCGATGTTGCCAGTAA","pam":"CGG","nuclease":"SpCas9","cut_offset_from_protospacer_start":17,"rank":3}]'
```

**Ask the inner agent**

> Propose `UpsertGuideSet` for exactly `g1`, `g2` and `g3` as declared demo records under `tp73_guides`. State that no TP73 reference is loaded and that the supplied `seq_id` and coordinates are not validated in this chapter; show the full records and wait for approval before registry mutation.

**Expected**

> The guide registry contains `tp73_guides` with three ranked guide candidates.

### Step 2: Apply practical filters (GC range, homopolymer limits, U6 terminator avoidance)

**GUI**

Apply practical filters (GC range, homopolymer limits, U6 terminator avoidance).

**CLI (terminal)**

```bash
gentle_cli guides filter tp73_guides --config '{"gc_min":0.3,"gc_max":0.7,"max_homopolymer_run":4,"reject_ambiguous_bases":true,"avoid_u6_terminator_tttt":true,"u6_terminator_window":"spacer_plus_tail","required_5prime_base":"G","allow_5prime_g_extension":true}' --output-set tp73_guides_pass
```

**Ask the inner agent**

> Propose `FilterGuidesPractical` on `tp73_guides` with GC 0.30..0.70, maximum homopolymer 4, ambiguous-base rejection, U6 `TTTT` avoidance over `spacer_plus_tail`, required 5-prime G and allowed G extension, writing subset `tp73_guides_pass`. Show the configuration first; after approval, report that `g1` and `g3` pass while `g2` fails `u6_terminator_t4`.

**Expected**

> The filter report records that `g1` and `g3` pass and `g2` fails `u6_terminator_t4`, and writes the passing subset as `tp73_guides_pass`.

### Step 3: Generate oligos from passed guides and inspect the resulting oligo set IDs

**GUI**

Generate oligos from passed guides and inspect the resulting oligo set IDs.

**CLI (terminal)**

```bash
gentle_cli guides oligos-generate tp73_guides lenti_bsmbi_u6_default --apply-5prime-g-extension --output-oligo-set tp73_lenti --passed-only
```

**Ask the inner agent**

> Propose `GenerateGuideOligos` from `tp73_guides` with template `lenti_bsmbi_u6_default`, `passed_only=true`, 5-prime-G extension enabled and output id `tp73_lenti`. Explain that `passed_only` reads the stored filter report on the source set, predict records only for `g1` and `g3`, and wait for separate approval. Return the exact oligo records without calling them target-validated or order-ready.

**Expected**

> The oligo registry contains two `tp73_lenti` records, generated only for `g1` and `g3` by consulting the filter report stored on `tp73_guides`.


## Ask an Outer Agent (MCP or ClawBio/OpenClaw)

An outer agent does not inherit the unsaved GUI project. Give it this chapter's canonical workflow and an explicit disposable state path; ask it to retain the structured result, artifacts and reproducibility receipt instead of replacing them with prose.

> Use GENtle's `gentle-cloning` skill to replay `docs/examples/workflows/guides_filter_and_generate_oligos.json` against a new disposable state. First report the exact workflow, inputs, state path, outputs and whether the selected route needs confirmation. Do not infer state from an open GUI. Return the structured result, produced artifacts and reproducibility receipt, and state any unmet prerequisite.

Equivalent direct structured request for the generic wrapper:

```json
{
  "schema": "gentle.clawbio_skill_request.v1",
  "mode": "workflow",
  "state_path": "/tmp/gentle-guides-filter-and-generate-oligos.state.json",
  "workflow_path": "docs/examples/workflows/guides_filter_and_generate_oligos.json",
  "timeout_secs": 300
}
```

Submitting a direct structured request is an explicit wrapper invocation. If natural language selects a narrower delegated skill and that route mutates state, selects biological material or writes artifacts, the caller must preserve that skill's proposal/approval boundary and approve only the exact bound digest. This tutorial generation step does not invoke an agent or grant approval.

## Interpretation and Reference

## Parameters That Matter

- `FilterGuidesPractical.config.gc_min / gc_max` (where used: operation 2)
  - Why it matters: GC bounds trade off stability and synthesis/efficiency behavior.
  - How to derive it: Start with literature-typical bounds (e.g., 0.3-0.7) and tighten per assay constraints.
- `GenerateGuideOligos.template_id` (where used: operation 3)
  - Why it matters: Template controls overhang/adaptor context for your cloning backbone.
  - How to derive it: Select the template matching your vector and cloning strategy.

## Applied Concepts

- **Deterministic Workflows** (`deterministic_workflows`): Operation chains should produce stable IDs and comparable outputs across repeated runs.
- **Guide Design Pipeline** (`guide_design_pipeline`): Guide sets can be created, filtered, expanded to oligos, and exported with protocol context.

## Checkpoints

- Guide set and passed guide set are present in metadata.
- Oligo set generation yields exactly two records for `g1` and `g3`; `g2` remains excluded by its U6 `TTTT` failure.
- The chapter explicitly states that its declared coordinates are not checked against a loaded TP73 reference and that local filtering does not establish off-target safety, activity, cloning success or order readiness.

## Tutorial Provenance

- Chapter id: `guides_filter_and_generate_oligos`
- Tier: `core`
- Example id: `guides_filter_and_generate_oligos`
- Tutorial source JSON: `docs/tutorial/sources/04-04_guides_filter_and_generate_oligos.json`
- Workflow file: `docs/examples/workflows/guides_filter_and_generate_oligos.json`
- Generated artifact dir: `docs/tutorial/generated/artifacts/guides_filter_and_generate_oligos`
- Example test_mode: `always`
- Executed during generation: `yes`
- Automated status: `passing`
- Review status: `codex_reviewed`
- Codex reviewed at: `2026-10-01`
- Human reviewed at: `not recorded`
- Inspect the source JSON when you need full option-level detail.

## Feedback

If this tutorial is confusing, execution-stale, biologically suspect, or missing a useful figure, please open the matching tutorial issue template and include the context below.

- Tutorial title: `Filter declared guides and format cloning oligos`
- Tutorial/chapter id: `guides_filter_and_generate_oligos`
- Step reached:
- Expected vs. actual:
- Interface used: GUI / CLI / Agent Assistant / ClawBio

Paste the Tutorial feedback context here:

```text

```
