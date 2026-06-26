---
chapter_id: "promoter_gene_set_ortholog_cohort_offline"
title: "Promoter Cohorts With Gene-Set and Ortholog Relationship Flags"
tier: "core"
example_id: "promoter_gene_set_ortholog_cohort_offline"
source_example: "docs/examples/workflows/promoter_gene_set_ortholog_cohort_offline.json"
example_test_mode: "always"
executed_during_generation: true
automated_status: "passing"
review_status: "unreviewed"
review_stale: false
codex_reviewed_at: null
human_reviewed_at: null
human_reviewer: null
review_stale_reason: null
review_issue_template: "Tutorial confusion"
review_issue_template_path: ".github/ISSUE_TEMPLATE/tutorial-confusion.md"
generated_artifact_dir: "docs/tutorial/generated/artifacts/promoter_gene_set_ortholog_cohort_offline"
---

# Promoter Cohorts With Gene-Set and Ortholog Relationship Flags

Run an offline TP73-centered promoter-cohort demo that records a declared co-regulation expectation, surfaces a deliberately synthetic unexpected-divergence motif flag, and repeats the idea across a local ortholog resource.

Promoter cohorts are useful only when GENtle keeps two ideas separate: how the member list was obtained, and what relationship the user expects among the promoters. A gene-set source can be explicit, catalog-backed, ontology-derived, neighboring, or random; the relationship expectation can independently be manual, co-regulated, or anti-co-regulated. This chapter uses a small artificial human-labeled genome where TP73 and E2F1 have similar SP1-rich promoter windows while PATZ1 is intentionally motif-poor. The declared expectation is co-regulation, so GENtle records the low TFBS-score-track similarity involving PATZ1 as an unexpected-divergence evidence flag, not as proof that the genes are or are not co-regulated.

The second half resolves TP73 ortholog promoter windows from a local `gentle.ortholog_resource.v1` file for synthetic human, mouse, and rat genomes. The mouse-labeled Trp73 promoter is deliberately different from the human- and rat-labeled TP73/Tp73 promoters, so the ortholog comparison emits relationship flags under the same conservative evidence-triage framing. All data are synthetic and local; the purpose is to teach the contract, not biology.

**Prerequisites:** Read [Chapter 24: Promoter Design Artifact Slice (Offline Synthetic TP73 Locus)](./08-03_promoter_design_artifact_slice_offline.md) first.

## Parameters That Matter

- `relationship=co_regulated` (where used: BuildGeneSetPromoterCohort, SummarizePromoterCohortComparison, ResolveOrthologPromoterCohort, SummarizeOrthologPromoterComparison)
  - Why it matters: The relationship is a declared expectation used to label evidence consistency or divergence; GENtle does not claim the promoters are biologically co-regulated.
  - How to derive it: Use `co_regulated` only when the cohort was selected for expected similarity, then inspect any unexpected-divergence flags as prompts for review.
- `upstream_bp=40 / downstream_bp=10` (where used: All promoter-window operations in this tutorial)
  - Why it matters: The local chromosomes are artificial and short; these values create compact 51 bp promoter windows around the synthetic TSS positions.
  - How to derive it: Use the fixed tutorial values so the inserted motif-rich and motif-poor segments dominate each promoter window.
- `motif SP1` (where used: Promoter cohort and ortholog promoter comparisons)
  - Why it matters: The motif-rich synthetic windows contain GC-rich SP1-like sequence, while the motif-poor windows do not, producing deterministic similarity differences.
  - How to derive it: Use the exact `SP1` token so the chapter stays compact and exercises one motif-track relationship flag.
- `local ortholog resource` (where used: ResolveOrthologPromoterCohort)
  - Why it matters: The tutorial demonstrates offline-first ortholog promoter resolution without online orthology services.
  - How to derive it: Use `docs/examples/assets/promoter_cohort_ortholog_demo/ortholog_resource.json`; real analyses should replace it with reviewed ortholog mappings.

## When This Routine Is Useful

- You want an offline example of gene-set promoter cohort comparison before prepared genome resources are available.
- You want to see how `--relationship co-regulated` changes report interpretation without becoming a biological verdict.
- You want to inspect an unexpected-divergence flag caused by one cohort member with a different promoter motif profile.
- You want a small local ortholog resource that resolves promoter windows without online genome or orthology downloads.
- You want generated JSON artifacts that can be compared by GUI, CLI, and automation tests.

## What You Learn

- Distinguish a gene-set source kind from a promoter relationship expectation.
- Explain why `relationship_flags` are evidence-triage annotations and not causal claims.
- Run a co-regulated promoter cohort comparison with a synthetic unexpected-divergence member.
- Resolve an ortholog promoter cohort from a local `gentle.ortholog_resource.v1` fixture.
- Inspect cross-species promoter relationship flags without requiring online ortholog or genome downloads.

## Applied Concepts

- **Shared Engine Contract** (`shared_engine_contract`): GUI, CLI, shell, and scripting interfaces execute the same operation semantics.
- **Genome Catalog Targeting** (`genome_catalog_targeting`): Prepared genome catalogs, annotation-based gene filters, and anchor extension connect imported entries to genomic context.
- **Promoter Motif Controls** (`promoter_motif_controls`): Foreground promoter motif signals should be compared with matched controls before being treated as candidate enrichment, depletion, or co-occurrence evidence.
- **Artifact Exports** (`artifact_exports`): Representative outputs (CSV/protocol/SVG/text) are retained for auditability and sharing.
- **Deterministic Workflows** (`deterministic_workflows`): Operation chains should produce stable IDs and comparable outputs across repeated runs.

## At a Glance

1. Open the tutorial chapter and inspect the retained JSON artifacts under docs/...
2. Use the Promoter design and gene-set inspection concepts from the previous pr...
3. Compare the relationship field in gene_set_promoter_cohort.json with the rela...
4. Open ortholog_promoter_cohort.json and confirm that the local ortholog resour...
5. Open ortholog_promoter_comparison.json and inspect the relationship_flags tha...

## GUI First

CLI snippets use GENtle's default `.gentle_state.json` state unless they say otherwise. Add `--state PATH` or `--project PATH` when you want an explicit sandboxed state file for copied commands.

### Step 1: Open the tutorial chapter and inspect the retained JSON artifacts under docs/...

GUI: Open the tutorial chapter and inspect the retained JSON artifacts under `docs/tutorial/generated/artifacts/promoter_gene_set_ortholog_cohort_offline/`.

CLI:

```bash
cargo run --bin gentle_cli -- workflow @docs/examples/workflows/promoter_gene_set_ortholog_cohort_offline.json
```

> Expected: The canonical workflow prepares the three synthetic local genomes into a temporary tutorial cache and writes four retained JSON artifacts.

### Step 2: Use the Promoter design and gene-set inspection concepts from the previous pr...

GUI: Use the Promoter design and gene-set inspection concepts from the previous promoter chapters; this slice is headless because it demonstrates cohort contracts rather than a new GUI-only workflow.

CLI:

```bash
cargo run --bin gentle_cli -- shell 'gene-sets promoter-cohort HumanPromoterToy --members TP73,E2F1,PATZ1 --relationship co-regulated --upstream-bp 40 --downstream-bp 10 --genome-catalog docs/examples/assets/promoter_cohort_ortholog_demo/genomes.json --cache-dir /tmp/gentle-promoter-cohort-demo-cache --output /tmp/gene_set_promoter_cohort.json'
```

> Expected: `gene_set_promoter_cohort.json` records an explicit-member source and `relationship: co_regulated` for TP73, E2F1, and PATZ1 promoter windows.

### Step 3: Compare the relationship field in gene_set_promoter_cohort.json with the rela...

GUI: Compare the `relationship` field in `gene_set_promoter_cohort.json` with the `relationship_flags` in `gene_set_promoter_cohort_comparison.json`.

CLI:

```bash
cargo run --bin gentle_cli -- shell 'genomes promoter-cohort-comparison HumanPromoterToy --cohort-label synthetic_tp73_e2f1_patz1_co_regulated_review --cohort-kind co-regulated --gene TP73 --gene E2F1 --gene PATZ1 --motif SP1 --upstream-bp 40 --downstream-bp 10 --catalog docs/examples/assets/promoter_cohort_ortholog_demo/genomes.json --cache-dir /tmp/gentle-promoter-cohort-demo-cache --path /tmp/gene_set_promoter_cohort_comparison.json'
```

> Expected: `gene_set_promoter_cohort_comparison.json` contains at least one `unexpected_divergence` flag involving the deliberately motif-poor PATZ1 promoter.

### Step 4: Open ortholog_promoter_cohort.json and confirm that the local ortholog resour...

GUI: Open `ortholog_promoter_cohort.json` and confirm that the local ortholog resource resolved human, mouse, and rat promoter rows.

CLI:

```bash
cargo run --bin gentle_cli -- shell 'orthologs resolve-promoter-cohort --anchor-species "Homo sapiens" --anchor-genome HumanPromoterToy --anchor-gene TP73 --target-species "Mus musculus" --target-species "Rattus norvegicus" --target-genome "Mus musculus=MousePromoterToy" --target-genome "Rattus norvegicus=RatPromoterToy" --transcript "Homo sapiens=TX_TP73_HUMAN" --transcript "Mus musculus=TX_TRP73_MOUSE" --transcript "Rattus norvegicus=TX_TP73_RAT" --orthologs docs/examples/assets/promoter_cohort_ortholog_demo/ortholog_resource.json --relationship co-regulated --upstream-bp 40 --downstream-bp 10 --catalog docs/examples/assets/promoter_cohort_ortholog_demo/genomes.json --cache-dir /tmp/gentle-promoter-cohort-demo-cache --path /tmp/ortholog_promoter_cohort.json'
```

> Expected: `ortholog_promoter_cohort.json` uses the local `gentle.ortholog_resource.v1` fixture to resolve Homo sapiens TP73, Mus musculus Trp73, and Rattus norvegicus Tp73 promoter windows.

### Step 5: Open ortholog_promoter_comparison.json and inspect the relationship_flags tha...

GUI: Open `ortholog_promoter_comparison.json` and inspect the `relationship_flags` that mark the synthetic mouse-labeled Trp73 promoter as divergent from the co-regulated expectation.

CLI:

```bash
cargo run --bin gentle_cli -- shell 'orthologs promoter-comparison --cohort /tmp/ortholog_promoter_cohort.json --motif SP1 --relationship co-regulated --path /tmp/ortholog_promoter_comparison.json'
```

> Expected: `ortholog_promoter_comparison.json` contains co-regulated relationship flags where the synthetic mouse promoter diverges from the human/rat motif profile.


## Follow-up Commands

```bash
cargo run --bin gentle_cli -- workflow @docs/examples/workflows/promoter_gene_set_ortholog_cohort_offline.json
cargo run --bin gentle_cli -- shell 'gene-sets promoter-cohort HumanPromoterToy --members TP73,E2F1,PATZ1 --relationship co-regulated --upstream-bp 40 --downstream-bp 10 --genome-catalog docs/examples/assets/promoter_cohort_ortholog_demo/genomes.json --cache-dir /tmp/gentle-promoter-cohort-demo-cache --output /tmp/gene_set_promoter_cohort.json'
cargo run --bin gentle_cli -- shell 'genomes promoter-cohort-comparison HumanPromoterToy --cohort-kind co-regulated --gene TP73 --gene E2F1 --gene PATZ1 --motif SP1 --upstream-bp 40 --downstream-bp 10 --catalog docs/examples/assets/promoter_cohort_ortholog_demo/genomes.json --cache-dir /tmp/gentle-promoter-cohort-demo-cache --path /tmp/gene_set_promoter_cohort_comparison.json'
cargo run --bin gentle_cli -- shell 'orthologs resolve-promoter-cohort --anchor-species "Homo sapiens" --anchor-genome HumanPromoterToy --anchor-gene TP73 --target-species "Mus musculus" --target-species "Rattus norvegicus" --target-genome "Mus musculus=MousePromoterToy" --target-genome "Rattus norvegicus=RatPromoterToy" --orthologs docs/examples/assets/promoter_cohort_ortholog_demo/ortholog_resource.json --relationship co-regulated --upstream-bp 40 --downstream-bp 10 --catalog docs/examples/assets/promoter_cohort_ortholog_demo/genomes.json --cache-dir /tmp/gentle-promoter-cohort-demo-cache --path /tmp/ortholog_promoter_cohort.json'
cargo run --bin gentle_cli -- shell 'orthologs promoter-comparison --cohort /tmp/ortholog_promoter_cohort.json --motif SP1 --relationship co-regulated --path /tmp/ortholog_promoter_comparison.json'
```

## Checkpoints

- The workflow executes offline and does not require prepared Ensembl/UCSC resources.
- `gene_set_promoter_cohort.json` returns three synthetic promoter windows and carries `relationship: co_regulated`.
- `gene_set_promoter_cohort_comparison.json` reports `unexpected_divergence` rather than silently treating all members as equivalent.
- `ortholog_promoter_cohort.json` reports three resolved promoter rows and zero unresolved rows.
- `ortholog_promoter_comparison.json` reports `relationship_flags` under the same evidence-not-verdict framing.

## What This Chapter Produces

- [`artifacts/promoter_gene_set_ortholog_cohort_offline/artifacts/promoter_cohort_ortholog_demo/gene_set_promoter_cohort.json`](../artifacts/promoter_gene_set_ortholog_cohort_offline/artifacts/promoter_cohort_ortholog_demo/gene_set_promoter_cohort.json) - schema: `gentle.gene_set_promoter_cohort.v1`
- [`artifacts/promoter_gene_set_ortholog_cohort_offline/artifacts/promoter_cohort_ortholog_demo/gene_set_promoter_cohort_comparison.json`](../artifacts/promoter_gene_set_ortholog_cohort_offline/artifacts/promoter_cohort_ortholog_demo/gene_set_promoter_cohort_comparison.json) - schema: `gentle.promoter_cohort_comparison.v1`
- [`artifacts/promoter_gene_set_ortholog_cohort_offline/artifacts/promoter_cohort_ortholog_demo/ortholog_promoter_cohort.json`](../artifacts/promoter_gene_set_ortholog_cohort_offline/artifacts/promoter_cohort_ortholog_demo/ortholog_promoter_cohort.json) - schema: `gentle.ortholog_promoter_cohort.v1`
- [`artifacts/promoter_gene_set_ortholog_cohort_offline/artifacts/promoter_cohort_ortholog_demo/ortholog_promoter_comparison.json`](../artifacts/promoter_gene_set_ortholog_cohort_offline/artifacts/promoter_cohort_ortholog_demo/ortholog_promoter_comparison.json) - schema: `gentle.ortholog_promoter_comparison.v1`

## Tutorial Provenance

- Chapter id: `promoter_gene_set_ortholog_cohort_offline`
- Tier: `core`
- Example id: `promoter_gene_set_ortholog_cohort_offline`
- Tutorial source JSON: `docs/tutorial/sources/08-07_promoter_gene_set_ortholog_cohort_offline.json`
- Workflow file: `docs/examples/workflows/promoter_gene_set_ortholog_cohort_offline.json`
- Generated artifact dir: `docs/tutorial/generated/artifacts/promoter_gene_set_ortholog_cohort_offline`
- Example test_mode: `always`
- Executed during generation: `yes`
- Automated status: `passing`
- Review status: `unreviewed`
- Codex reviewed at: `not recorded`
- Human reviewed at: `not recorded`
- Inspect the source JSON when you need full option-level detail.

## Feedback

If this tutorial is confusing, execution-stale, biologically suspect, or missing a useful figure, please open the matching tutorial issue template and include the context below.

- Tutorial title: `Promoter Cohorts With Gene-Set and Ortholog Relationship Flags`
- Tutorial/chapter id: `promoter_gene_set_ortholog_cohort_offline`
- Step reached:
- Expected vs. actual:
- Interface used: GUI / CLI / Agent Assistant / ClawBio

Paste the Tutorial feedback context here:

```text

```
