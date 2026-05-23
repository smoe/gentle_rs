---
chapter_id: "guides_filter_and_generate_oligos"
title: "Guide practical filtering and oligo generation"
tier: "core"
example_id: "guides_filter_and_generate_oligos"
source_example: "docs/examples/workflows/guides_filter_and_generate_oligos.json"
example_test_mode: "always"
executed_during_generation: true
automated_status: "passing"
review_status: "unreviewed"
review_stale: false
codex_reviewed_at: null
human_reviewed_at: null
human_reviewer: null
generated_artifact_dir: "docs/tutorial/generated/artifacts/guides_filter_and_generate_oligos"
---

# Guide practical filtering and oligo generation

Apply practical guide constraints and produce cloning-ready oligo candidates.

For CRISPR-style cloning, guide quality control is where many downstream failures are prevented. This routine demonstrates how to encode practical constraints directly in operations, then generate oligos from the passed candidates.

**Prerequisites:** Read [Chapter 1: Load FASTA, branch, and reverse-complement](./02-01_load_branch_reverse_complement_pgex_fasta.md) first.

## Parameters That Matter

- `FilterGuidesPractical.config.gc_min / gc_max` (where used: operation 2)
  - Why it matters: GC bounds trade off stability and synthesis/efficiency behavior.
  - How to derive it: Start with literature-typical bounds (e.g., 0.3-0.7) and tighten per assay constraints.
- `GenerateGuideOligos.template_id` (where used: operation 3)
  - Why it matters: Template controls overhang/adaptor context for your cloning backbone.
  - How to derive it: Select the template matching your vector and cloning strategy.

## When This Routine Is Useful

- You need a transparent guide filtering step before oligo ordering.
- You want to document why specific guides were excluded.
- You need deterministic oligo generation from a reusable guide set.

## What You Learn

- Create and filter guide sets through explicit engine operations.
- Generate oligo records from filtered guide candidates.
- Connect guide workflows to the same deterministic operation model used for cloning steps.

## Applied Concepts

- **Deterministic Workflows** (`deterministic_workflows`): Operation chains should produce stable IDs and comparable outputs across repeated runs.
- **Guide Design Pipeline** (`guide_design_pipeline`): Guide sets can be created, filtered, expanded to oligos, and exported with protocol context.

## GUI First

CLI snippets use GENtle's default `.gentle_state.json` state unless they say otherwise. Add `--state PATH` or `--project PATH` when you want an explicit sandboxed state file for copied commands.

### Step 1: Open the guides workflow controls in GENtle and create/import a guide set for...

GUI: Open the guides workflow controls in GENtle and create/import a guide set for a target region.

CLI:

```bash
cargo run --bin gentle_cli -- guides put tp73_guides --json '[{"guide_id":"g1","seq_id":"tp73","start_0based":100,"end_0based_exclusive":120,"strand":"+","protospacer":"GACCTGTTGACGATGTTCCA","pam":"AGG","nuclease":"SpCas9","cut_offset_from_protospacer_start":17,"rank":1},{"guide_id":"g2","seq_id":"tp73","start_0based":220,"end_0based_exclusive":240,"strand":"+","protospacer":"TTTTGCCATGTTGACCTGAA","pam":"TGG","nuclease":"SpCas9","cut_offset_from_protospacer_start":17,"rank":2},{"guide_id":"g3","seq_id":"tp73","start_0based":340,"end_0based_exclusive":360,"strand":"-","protospacer":"GGTACCGATGTTGCCAGTAA","pam":"CGG","nuclease":"SpCas9","cut_offset_from_protospacer_start":17,"rank":3}]'
```

> Expected: The guide registry contains `tp73_guides` with three ranked guide candidates.

### Step 2: Apply practical filters (GC range, homopolymer limits, U6 terminator avoidance)

GUI: Apply practical filters (GC range, homopolymer limits, U6 terminator avoidance).

CLI:

```bash
cargo run --bin gentle_cli -- guides filter tp73_guides --config '{"gc_min":0.3,"gc_max":0.7,"max_homopolymer_run":4,"reject_ambiguous_bases":true,"avoid_u6_terminator_tttt":true,"u6_terminator_window":"spacer_plus_tail","required_5prime_base":"G","allow_5prime_g_extension":true}' --output-set tp73_guides_pass
```

> Expected: The filter report records pass/fail decisions and writes the passing subset as `tp73_guides_pass`.

### Step 3: Generate oligos from passed guides and inspect the resulting oligo set IDs

GUI: Generate oligos from passed guides and inspect the resulting oligo set IDs.

CLI:

```bash
cargo run --bin gentle_cli -- guides oligos-generate tp73_guides lenti_bsmbi_u6_default --apply-5prime-g-extension --output-oligo-set tp73_lenti --passed-only
```

> Expected: The oligo registry contains `tp73_lenti`, generated only from guides that passed the practical filter.


## Checkpoints

- Guide set and passed guide set are present in metadata.
- Oligo set generation succeeds for passed guides.

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
- Review status: `unreviewed`
- Codex reviewed at: `not recorded`
- Human reviewed at: `not recorded`
- Inspect the source JSON when you need full option-level detail.

## Feedback

If this tutorial is confusing, execution-stale, biologically suspect, or missing a useful figure, please open the matching tutorial issue template and include the context below.

- Tutorial title: `Guide practical filtering and oligo generation`
- Tutorial/chapter id: `guides_filter_and_generate_oligos`
- Step reached:
- Expected vs. actual:
- Interface used: GUI / CLI / Agent Assistant / ClawBio

Paste the Tutorial feedback context here:

```text

```
