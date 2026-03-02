# Digest -> Ligation -> ExtractRegion minimal slice

- Chapter id: `digest_ligation_extract_region_minimal`
- Tier: `core`
- Example id: `digest_ligation_extract_region_minimal`
- Source example: `docs/examples/workflows/digest_ligation_extract_region_minimal.json`
- Example test_mode: `always`
- Executed during generation: `yes`

Run a full mini-loop from fragment production to assembled product extraction.

This chapter models a compact molecular cloning routine in one chain: digest source material, produce a ligation product, and extract a target segment for subsequent validation or design. It is the smallest end-to-end routine that still reflects real bench-side reasoning.

## When This Routine Is Useful

- You want to test whether your planned digest/ligation sequence is internally consistent.
- You need a deterministic extracted segment for primer design or annotation checks.
- You want a regression slice that exercises core cloning operations together.

## What You Learn

- Execute a minimal end-to-end cloning chain.
- Track how intermediate IDs are consumed by downstream operations.
- Use ExtractRegion output as a stable hand-off point for later analyses.

## Concepts and Recurrence

- **Deterministic Workflows** (`deterministic_workflows`): Operation chains should produce stable IDs and comparable outputs across repeated runs.
  - Status: reinforced from [Chapter 1: Load FASTA, branch, and reverse-complement](./01_load_branch_reverse_complement_pgex_fasta.md), [Chapter 2: Load pGEX and digest with BamHI/EcoRI](./02_load_and_digest_pgex.md), [Chapter 3: Guide practical filtering and oligo generation](./03_guides_filter_and_generate_oligos.md).
  - Reoccurs in: [Chapter 7: Prepare a reference genome cache (online)](./07_prepare_reference_genome_online.md).
- **Sequence Lineage** (`sequence_lineage`): Derived sequences are explicit products linked to upstream inputs and operations.
  - Status: reinforced from [Chapter 1: Load FASTA, branch, and reverse-complement](./01_load_branch_reverse_complement_pgex_fasta.md), [Chapter 2: Load pGEX and digest with BamHI/EcoRI](./02_load_and_digest_pgex.md).
  - Reoccurs in: no later chapter.

## GUI First

1. Start from the loaded FASTA/plasmid sequence in the GUI.
2. Run digest with selected enzymes and inspect available products.
3. Run ligation with the intended inputs, then extract a focused region from the ligation result.

## Command Equivalent (After GUI)

Run the same routine non-interactively once the GUI flow is clear:

```bash
cargo run --bin gentle_cli -- workflow @docs/examples/workflows/digest_ligation_extract_region_minimal.json
cargo run --bin gentle_cli -- shell 'workflow @docs/examples/workflows/digest_ligation_extract_region_minimal.json'
```

## Parameters That Matter

- `Ligation.protocol` (where used: operation 3)
  - Why it matters: Protocol controls compatibility logic (e.g., blunt vs sticky behavior).
  - How to derive it: Choose based on end chemistry produced by upstream digest products.
- `ExtractRegion.from / to` (where used: operation 4)
  - Why it matters: Defines the exact segment handed to downstream interpretation.
  - How to derive it: Derive boundaries from feature coordinates or expected amplicon/design window.

## Checkpoints

- Digest creates at least two fragments for ligation input.
- Ligation creates deterministic output IDs.
- ExtractRegion creates `lig_extract`.

## Retained Outputs

- None for this chapter.

## Canonical Source

- Workflow file: `docs/examples/workflows/digest_ligation_extract_region_minimal.json`
- Inspect this JSON file directly when you need full option-level detail.
