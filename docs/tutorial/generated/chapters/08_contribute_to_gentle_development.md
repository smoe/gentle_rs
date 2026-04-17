# Contribute to GENtle development

- Chapter id: `contribute_to_gentle_development`
- Tier: `advanced`
- Example id: `contribute_gentle_development_baseline`
- Source example: `docs/examples/workflows/contribute_gentle_development_baseline.json`
- Example test_mode: `always`
- Executed during generation: `yes`

Executable contributor onboarding from baseline run to test/documentation checks.

Contributing effectively requires keeping biological behavior, command contracts, and docs in sync. This chapter provides a concrete contributor routine: run a baseline operation chain, then validate examples, tutorial generation, and tests before opening a change request.

## When This Routine Is Useful

- You want to submit a feature or fix without breaking cross-interface behavior.
- You need to verify that tutorial/docs stay synchronized with executable workflows.
- You want a deterministic pre-PR checklist for local validation.

## What You Learn

- Understand the expected contributor loop from change to verification.
- Use tutorial and example checks to prevent documentation drift.
- Identify the minimum validation commands expected before handoff or PR.

## Concepts and Recurrence

- **Shared Engine Contract** (`shared_engine_contract`): GUI, CLI, shell, and scripting interfaces execute the same operation semantics.
  - Status: reinforced from [Chapter 1: Load FASTA, branch, and reverse-complement](./01_load_branch_reverse_complement_pgex_fasta.md), [Chapter 2: Find and extend the right genomic target (local catalog)](./02_find_and_extend_genomic_target_local_catalog.md), [Chapter 4: Gibson two-fragment overlap planning baseline](./04_gibson_two_fragment_overlap_preview.md).
  - Reoccurs in: [Chapter 10: TP53 isoform architecture expert panel (online)](./10_tp53_isoform_architecture_online.md), [Chapter 11: Retrieve TP63 and extend the displayed region by +/-2 kb (online)](./11_tp63_anchor_extension_online.md), [Chapter 12: Map TP53 locus reads with multi-gene sparse indexing (online)](./12_tp53_multi_gene_sparse_mapping_online.md), [Chapter 13: Selection-first PCR batch primer design (offline)](./13_pcr_selection_batch_primer_pairs_offline.md), [Chapter 14: Compare TP73 cDNA against TP73 genomic context via dotplot (online)](./14_tp73_cdna_genomic_dotplot_online.md), [Chapter 15: Gibson Specialist Starter Project (offline)](./15_gibson_specialist_testing_baseline.md), [Chapter 16: Gibson Arrangements Starter Project (offline)](./16_gibson_arrangements_baseline.md), [Chapter 17: TP53 UniProt domain mapping and feature-coding DNA query (online)](./17_tp53_uniprot_projection_online.md), [Chapter 18: Simple PCR From a Selected Core Region](./18_simple_pcr_selection_gui.md), [Chapter 19: Audit a TP73 UniProt Projection Against Ensembl and Derived Coding Sequence (CLI Tutorial)](./19_tp73_uniprot_projection_audit_cli.md).
- **Tutorial Drift Checks** (`tutorial_drift_checks`): Tutorial content is generated from executable examples and verified in automated checks.
  - Status: introduced in this chapter.
  - Reoccurs in: no later chapter.
- **Contribution Loop** (`contribution_loop`): Contributions should couple code edits with docs updates and deterministic tests.
  - Status: introduced in this chapter.
  - Reoccurs in: no later chapter.

## GUI First

1. Run the baseline sequence routine in the GUI and inspect resulting lineage entries.
2. Locate the same routine as a canonical workflow JSON example in `docs/examples/workflows`.
3. Use this mapping to validate that your planned code change affects shared engine behavior, not only one interface.

## Command Equivalent (After GUI)

Run the same routine non-interactively once the GUI flow is clear:

```bash
cargo run --bin gentle_cli -- workflow @docs/examples/workflows/contribute_gentle_development_baseline.json
cargo run --bin gentle_cli -- shell 'workflow @docs/examples/workflows/contribute_gentle_development_baseline.json'
```

## Parameters That Matter

- `Workflow file path @docs/examples/workflows/...` (where used: command equivalent)
  - Why it matters: Ensures you are validating the same canonical routine used by docs/tests.
  - How to derive it: Select the chapter's `example_id` and open the matching JSON file under `docs/examples/workflows`.

## Follow-up Commands

```bash
cargo check -q
cargo test -q workflow_examples -- --test-threads=1
cargo run --bin gentle_examples_docs -- --check
cargo run --bin gentle_examples_docs -- tutorial-check
cargo run --bin gentle_cli -- capabilities
```

## Checkpoints

- Baseline contribution workflow executes without failures.
- Contributor validation commands pass in a clean environment.
- Tutorial-generated content remains synchronized with source manifests/examples.

## Retained Outputs

- None for this chapter.

## Canonical Source

- Workflow file: `docs/examples/workflows/contribute_gentle_development_baseline.json`
- Inspect this JSON file directly when you need full option-level detail.
