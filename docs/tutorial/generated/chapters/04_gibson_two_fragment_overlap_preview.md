# Gibson two-fragment overlap planning baseline

- Chapter id: `gibson_two_fragment_overlap_preview`
- Tier: `core`
- Example id: `gibson_two_fragment_overlap_preview`
- Source example: `docs/examples/workflows/gibson_two_fragment_overlap_preview.json`
- Example test_mode: `always`
- Executed during generation: `yes`

Use the built-in Gibson routine baseline to validate overlap assumptions and generate deterministic preview assemblies.

Gibson assembly in practice depends on overlap design quality. This chapter anchors that design logic in GENtle through one repeatable workflow: build two overlapping fragments, run Gibson-specific preflight checks on overlap compatibility, and produce forward-order preview sequences for review and communication.

## When This Routine Is Useful

- You want reproducible documentation of intended Gibson fragment order and overlap assumptions.
- You need fast preflight feedback before primer ordering or wet-lab execution.
- You want one routine that can be explained to collaborators through GUI + shell parity.

## What You Learn

- Understand where Gibson overlap assumptions are checked in shared-shell preflight.
- Interpret overlap mismatch diagnostics before state mutation.
- Use deterministic output IDs for clear cloning communication.

## Concepts and Recurrence

- **Shared Engine Contract** (`shared_engine_contract`): GUI, CLI, shell, and scripting interfaces execute the same operation semantics.
  - Status: reinforced from [Chapter 1: Load FASTA, branch, and reverse-complement](./01_load_branch_reverse_complement_pgex_fasta.md), [Chapter 2: Find and extend the right genomic target (local catalog)](./02_find_and_extend_genomic_target_local_catalog.md).
  - Reoccurs in: [Chapter 8: Contribute to GENtle development](./08_contribute_to_gentle_development.md), [Chapter 10: TP53 isoform architecture expert panel (online)](./10_tp53_isoform_architecture_online.md), [Chapter 11: Retrieve TP63 and extend the displayed region by +/-2 kb (online)](./11_tp63_anchor_extension_online.md).
- **Deterministic Workflows** (`deterministic_workflows`): Operation chains should produce stable IDs and comparable outputs across repeated runs.
  - Status: reinforced from [Chapter 1: Load FASTA, branch, and reverse-complement](./01_load_branch_reverse_complement_pgex_fasta.md), [Chapter 2: Find and extend the right genomic target (local catalog)](./02_find_and_extend_genomic_target_local_catalog.md), [Chapter 3: Load pGEX and digest with BamHI/EcoRI](./03_load_and_digest_pgex.md).
  - Reoccurs in: [Chapter 5: Guide practical filtering and oligo generation](./05_guides_filter_and_generate_oligos.md), [Chapter 6: Digest -> Ligation -> ExtractRegion minimal slice](./06_digest_ligation_extract_region_minimal.md), [Chapter 9: Prepare a reference genome cache (online)](./09_prepare_reference_genome_online.md).
- **Sequence Lineage** (`sequence_lineage`): Derived sequences are explicit products linked to upstream inputs and operations.
  - Status: reinforced from [Chapter 1: Load FASTA, branch, and reverse-complement](./01_load_branch_reverse_complement_pgex_fasta.md), [Chapter 2: Find and extend the right genomic target (local catalog)](./02_find_and_extend_genomic_target_local_catalog.md), [Chapter 3: Load pGEX and digest with BamHI/EcoRI](./03_load_and_digest_pgex.md).
  - Reoccurs in: [Chapter 6: Digest -> Ligation -> ExtractRegion minimal slice](./06_digest_ligation_extract_region_minimal.md), [Chapter 11: Retrieve TP63 and extend the displayed region by +/-2 kb (online)](./11_tp63_anchor_extension_online.md).

## GUI First

1. Load `test_files/pGEX_3X.fa` and create two overlapping fragments (`gibson_left`, `gibson_right`) using region extraction.
2. Open `Patterns -> Routine catalog`, locate `Gibson Two-Fragment Overlap Preview`, and import it.
3. Run from `Shell`: `macros template-run gibson_two_fragment_overlap_preview --bind left_seq_id=gibson_left --bind right_seq_id=gibson_right --bind overlap_bp=20 --bind assembly_prefix=gibson_demo --bind output_id=gibson_demo_forward --validate-only`, then execute without `--validate-only`.

## Command Equivalent (After GUI)

Run the same routine non-interactively once the GUI flow is clear:

```bash
cargo run --bin gentle_cli -- workflow @docs/examples/workflows/gibson_two_fragment_overlap_preview.json
cargo run --bin gentle_cli -- shell 'workflow @docs/examples/workflows/gibson_two_fragment_overlap_preview.json'
```

## Parameters That Matter

- `left_seq_id / right_seq_id` (where used: template preflight + run)
  - Why it matters: Order defines the expected assembled product direction and adjacent overlap checks.
  - How to derive it: Choose fragment IDs in intended 5'->3' assembly order.
- `overlap_bp` (where used: Gibson family preflight check)
  - Why it matters: Controls suffix/prefix overlap length validation between adjacent fragments.
  - How to derive it: Use planned primer-homology length (commonly 15-40 bp for many Gibson workflows).
- `assembly_prefix / output_id` (where used: template run output naming)
  - Why it matters: Stable IDs make review/export/communication unambiguous.
  - How to derive it: Use project-specific names (for example `tp73_gibson_round1`).
  - Omit when: Omit only when default IDs are acceptable for exploratory runs.

## Follow-up Commands

```bash
cargo run --bin gentle_cli -- shell 'macros template-import assets/cloning_patterns_catalog/gibson/overlap_assembly/gibson_two_fragment_overlap_preview.json'
cargo run --bin gentle_cli -- shell 'macros template-run gibson_two_fragment_overlap_preview --bind left_seq_id=gibson_left --bind right_seq_id=gibson_right --bind overlap_bp=20 --bind assembly_prefix=gibson_demo --bind output_id=gibson_demo_forward --validate-only'
cargo run --bin gentle_cli -- shell 'macros template-run gibson_two_fragment_overlap_preview --bind left_seq_id=gibson_left --bind right_seq_id=gibson_right --bind overlap_bp=20 --bind assembly_prefix=gibson_demo --bind output_id=gibson_demo_forward --transactional'
```

## Checkpoints

- Validate-only run returns `can_execute=true` for matching overlap fragments.
- Mismatch overlap inputs return explicit Gibson overlap diagnostics.
- Mutating run creates deterministic preview outputs (`${assembly_prefix}_*` and `${output_id}`).

## Retained Outputs

- None for this chapter.

## Canonical Source

- Workflow file: `docs/examples/workflows/gibson_two_fragment_overlap_preview.json`
- Inspect this JSON file directly when you need full option-level detail.
