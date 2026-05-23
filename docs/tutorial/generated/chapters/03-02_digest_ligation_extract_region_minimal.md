---
chapter_id: "digest_ligation_extract_region_minimal"
title: "Digest -> Ligation -> ExtractRegion minimal slice"
tier: "core"
example_id: "digest_ligation_extract_region_minimal"
source_example: "docs/examples/workflows/digest_ligation_extract_region_minimal.json"
example_test_mode: "always"
executed_during_generation: true
automated_status: "passing"
review_status: "unreviewed"
codex_reviewed_at: null
human_reviewed_at: null
generated_artifact_dir: "docs/tutorial/generated/artifacts/digest_ligation_extract_region_minimal"
---

# Digest -> Ligation -> ExtractRegion minimal slice

Run a full mini-loop from fragment production to assembled product extraction.

This chapter models a compact molecular cloning routine in one chain: digest source material, produce a ligation product, and extract a target segment for subsequent validation or design. It is the smallest end-to-end routine that still reflects real bench-side reasoning.

**Prerequisites:** Read [Chapter 3: Load pGEX and digest with BamHI/EcoRI](./03-01_load_and_digest_pgex.md) first.

## Parameters That Matter

- `Ligation.protocol` (where used: operation 3)
  - Why it matters: Protocol controls compatibility logic (e.g., blunt vs sticky behavior).
  - How to derive it: Choose based on end chemistry produced by upstream digest products.
- `ExtractRegion.from / to` (where used: operation 4)
  - Why it matters: Defines the exact segment handed to downstream interpretation.
  - How to derive it: Derive boundaries from feature coordinates or expected amplicon/design window.

## When This Routine Is Useful

- You want to test whether your planned digest/ligation sequence is internally consistent.
- You need a deterministic extracted segment for primer design or annotation checks.
- You want a regression slice that exercises core cloning operations together.

## What You Learn

- Execute a minimal end-to-end cloning chain.
- Track how intermediate IDs are consumed by downstream operations.
- Use ExtractRegion output as a stable hand-off point for later analyses.

## Applied Concepts

- **Deterministic Workflows** (`deterministic_workflows`): Operation chains should produce stable IDs and comparable outputs across repeated runs.
- **Sequence Lineage** (`sequence_lineage`): Derived sequences are explicit products linked to upstream inputs and operations.

## GUI First

CLI snippets use GENtle's default `.gentle_state.json` state unless they say otherwise. Add `--state PATH` or `--project PATH` when you want an explicit sandboxed state file for copied commands.

### Step 1: Start from the loaded FASTA/plasmid sequence in the GUI

GUI: Start from the loaded FASTA/plasmid sequence in the GUI.

CLI:

```bash
cargo run --bin gentle_cli -- op '{"LoadFile":{"path":"test_files/pGEX_3X.fa","as_id":"pgex_fasta"}}'
```

> Expected: The source sequence is loaded as `pgex_fasta`, ready to feed cloning operations.

### Step 2: Run digest with selected enzymes and inspect available products

GUI: Run digest with selected enzymes and inspect available products.

CLI:

```bash
cargo run --bin gentle_cli -- op '{"Digest":{"input":"pgex_fasta","enzymes":["BamHI","EcoRI"],"output_prefix":"d"}}'
```

> Expected: The digest step creates deterministic fragment IDs with the `d` prefix.

### Step 3: Run ligation with the intended inputs, then extract a focused region from the...

GUI: Run ligation with the intended inputs, then extract a focused region from the ligation result.

CLI:

```bash
cargo run --bin gentle_cli -- workflow @docs/examples/workflows/digest_ligation_extract_region_minimal.json
```

> Expected: The full workflow creates ligation product `lig_1` and extracted handoff sequence `lig_extract`.


## Checkpoints

- Digest creates at least two fragments for ligation input.
- Ligation creates deterministic output IDs.
- ExtractRegion creates `lig_extract`.

## Tutorial Provenance

- Chapter id: `digest_ligation_extract_region_minimal`
- Tier: `core`
- Example id: `digest_ligation_extract_region_minimal`
- Tutorial source JSON: `docs/tutorial/sources/03-02_digest_ligation_extract_region_minimal.json`
- Workflow file: `docs/examples/workflows/digest_ligation_extract_region_minimal.json`
- Generated artifact dir: `docs/tutorial/generated/artifacts/digest_ligation_extract_region_minimal`
- Example test_mode: `always`
- Executed during generation: `yes`
- Automated status: `passing`
- Review status: `unreviewed`
- Codex reviewed at: `not recorded`
- Human reviewed at: `not recorded`
- Inspect the source JSON when you need full option-level detail.

## Feedback

If this tutorial is confusing, execution-stale, biologically suspect, or missing a useful figure, please open the matching tutorial issue template and include the context below.

- Tutorial title: `Digest -> Ligation -> ExtractRegion minimal slice`
- Tutorial/chapter id: `digest_ligation_extract_region_minimal`
- Step reached:
- Expected vs. actual:
- Interface used: GUI / CLI / Agent Assistant / ClawBio

Paste the Tutorial feedback context here:

```text

```
