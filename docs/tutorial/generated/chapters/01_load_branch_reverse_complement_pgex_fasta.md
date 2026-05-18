---
chapter_id: "load_branch_reverse_complement_pgex_fasta"
title: "Load FASTA, branch, and reverse-complement"
tier: "core"
example_id: "load_branch_reverse_complement_pgex_fasta"
source_example: "docs/examples/workflows/load_branch_reverse_complement_pgex_fasta.json"
example_test_mode: "always"
executed_during_generation: true
---

# Load FASTA, branch, and reverse-complement

Start with a simple cloning-prep routine to learn product IDs and sequence lineage.

In real cloning projects, people often say "the sequence" even when they mean different biological records: the sequence as downloaded from a public database or read from a local sample, versus a derivative changed by an operation such as introducing a nonsense mutation. GENtle keeps those records apart. Each conceptual or hands-on sequence is represented as a vial node in a lineage graph, and each operation is represented as a transaction node connecting biological inputs to outputs. Even this small routine can branch: the original FASTA remains the provenance-preserving input, while the branch and reverse-complement become separate products. Larger projects can branch one product into several downstream operations or merge independent products into one process.

## Parameters That Matter

- `LoadFile.path` (where used: operation 1)
  - Why it matters: Defines the exact source sequence. Wrong file path means a different biological starting point.
  - How to derive it: Use the file you just loaded in the GUI for this routine (`test_files/pGEX_3X.fa`).
- `Branch.output_id / ReverseComplement.output_id` (where used: operations 2 and 3)
  - Why it matters: Stable IDs make downstream commands unambiguous.
  - How to derive it: Pick short descriptive IDs tied to intent (e.g., branch copy vs reverse-complement product).
  - Omit when: You can omit IDs only when you accept auto-generated names.

## When This Routine Is Useful

- You received a FASTA plasmid or amplicon sequence and want a reproducible project starting point.
- You need a reverse-complemented working copy but want to preserve the original sequence identity.
- You are onboarding to GENtle and need to understand how derived IDs are created and reused.

## What You Learn

- Run one canonical workflow from source JSON and interpret the result IDs.
- Understand that branch and reverse-complement operations create explicit derivative sequences.
- Recognize the shared-engine contract across GUI, CLI, and shell entry points.

## Applied Concepts

- **Shared Engine Contract** (`shared_engine_contract`): GUI, CLI, shell, and scripting interfaces execute the same operation semantics.
- **Deterministic Workflows** (`deterministic_workflows`): Operation chains should produce stable IDs and comparable outputs across repeated runs.
- **Sequence Lineage** (`sequence_lineage`): Derived sequences are explicit products linked to upstream inputs and operations.

## GUI First

1. Open GENtle and load `test_files/pGEX_3X.fa` via `File -> Open`.
2. In the DNA window, create a branch copy from the loaded sequence (Branch action).
3. Apply reverse-complement to the branch and confirm a new sequence entry appears in lineage/table views.

## Command Equivalent (After GUI)

Run the same routine non-interactively once the GUI flow is clear:

```bash
cargo run --bin gentle_cli -- workflow @docs/examples/workflows/load_branch_reverse_complement_pgex_fasta.json
cargo run --bin gentle_cli -- shell 'workflow @docs/examples/workflows/load_branch_reverse_complement_pgex_fasta.json'
```

## Checkpoints

- Workflow executes without warnings or errors.
- Derived sequence IDs include a branch and reverse-complement product.

## Canonical Source

- Chapter id: `load_branch_reverse_complement_pgex_fasta`
- Tier: `core`
- Example id: `load_branch_reverse_complement_pgex_fasta`
- Workflow file: `docs/examples/workflows/load_branch_reverse_complement_pgex_fasta.json`
- Example test_mode: `always`
- Executed during generation: `yes`
- Inspect this JSON file directly when you need full option-level detail.
