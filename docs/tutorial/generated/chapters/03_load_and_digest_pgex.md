---
chapter_id: "load_and_digest_pgex"
title: "Load pGEX and digest with BamHI/EcoRI"
tier: "core"
example_id: "load_and_digest_pgex"
source_example: "docs/examples/workflows/load_and_digest_pgex.json"
example_test_mode: "always"
executed_during_generation: true
---

# Load pGEX and digest with BamHI/EcoRI

Introduce restriction digest planning and deterministic fragment products.

Restriction digest is a core molecular cloning routine used for vector linearization, insert release, and diagnostic fragment checks. This chapter focuses on how digest parameters map to reproducible fragment sets that can later feed ligation or analysis steps.

## Parameters That Matter

- `Digest.enzymes` (where used: operation 2)
  - Why it matters: The enzyme set determines cut positions and resulting fragment repertoire.
  - How to derive it: Choose enzymes based on cloning strategy (diagnostic digest vs insert release vs vector opening).
- `Digest.output_prefix` (where used: operation 2)
  - Why it matters: Controls deterministic fragment ID namespace.
  - How to derive it: Use a short routine-specific prefix (e.g., `frag`, `d`, `eco_bam`).
  - Omit when: Omit only if auto-generated IDs are acceptable for ad hoc inspection.

## When This Routine Is Useful

- You want to verify expected restriction fragments before ordering primers or designing ligations.
- You need a reproducible digest baseline to compare against wet-lab gel expectations.
- You plan to reuse fragment IDs in later operations.

## What You Learn

- Execute digest operations and inspect created fragment IDs.
- Reason about multi-product lineage from one parent sequence.
- Identify why stable IDs matter for follow-up ligation/extraction steps.

## Applied Concepts

- **Deterministic Workflows** (`deterministic_workflows`): Operation chains should produce stable IDs and comparable outputs across repeated runs.
- **Sequence Lineage** (`sequence_lineage`): Derived sequences are explicit products linked to upstream inputs and operations.

## GUI First

1. Load `test_files/pGEX-3X.gb` in the GUI and inspect annotated features.
2. Run Digest from the DNA window using enzymes `BamHI` and `EcoRI`.
3. Review created fragment entries and confirm they are stored as independent sequence products.

## Command Equivalent (After GUI)

Run the same routine non-interactively once the GUI flow is clear:

```bash
cargo run --bin gentle_cli -- workflow @docs/examples/workflows/load_and_digest_pgex.json
cargo run --bin gentle_cli -- shell 'workflow @docs/examples/workflows/load_and_digest_pgex.json'
```

## Checkpoints

- Digest operation completes and creates fragment sequence IDs.
- Fragment IDs are deterministic across repeated runs.

## Canonical Source

- Chapter id: `load_and_digest_pgex`
- Tier: `core`
- Example id: `load_and_digest_pgex`
- Workflow file: `docs/examples/workflows/load_and_digest_pgex.json`
- Example test_mode: `always`
- Executed during generation: `yes`
- Inspect this JSON file directly when you need full option-level detail.
