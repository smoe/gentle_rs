---
chapter_id: "gene_set_ortholog_promoter_cohorts_offline"
title: "Gene-Set and Ortholog Promoter Cohorts (Offline Synthetic TP73 Demo)"
tier: "core"
example_id: "gene_set_ortholog_promoter_cohorts_offline"
source_example: "docs/examples/workflows/gene_set_ortholog_promoter_cohorts_offline.json"
example_test_mode: "always"
executed_during_generation: true
automated_status: "passing"
review_status: "codex_reviewed"
review_stale: false
codex_reviewed_at: "2026-06-25"
human_reviewed_at: null
human_reviewer: null
review_stale_reason: null
review_issue_template: null
review_issue_template_path: null
generated_artifact_dir: "docs/tutorial/generated/artifacts/gene_set_ortholog_promoter_cohorts_offline"
---

# Gene-Set and Ortholog Promoter Cohorts (Offline Synthetic TP73 Demo)

Run a fully offline promoter-cohort demo that resolves a reviewed gene set, records a co-regulated relationship expectation, shows an unexpected TFBS-divergence flag, then repeats the same expectation framing for a local TP73/Trp73 ortholog promoter comparison.

Relationship flags are triage aids, not biological verdicts. In this chapter, GENtle prepares two tiny synthetic genome fixtures, resolves a reviewed human gene set (`TP73` plus a deliberately divergent synthetic partner), and records that the set is expected to be co-regulated. The follow-up promoter comparison finds that the available SP1 score-track evidence is not concordant, so it emits an `unexpected_divergence` flag.

The second half uses a local `gentle.ortholog_resource.v1` fixture mapping synthetic human `TP73` to mouse `Trp73`. It resolves strand-aware promoter windows from local prepared genome indexes and then compares promoter TFBS evidence with the same conservative co-regulated expectation. This is cross-species association evidence only; the fixture is artificial and does not prove regulation, conservation, or orthology in real organisms.

**Prerequisites:** Read [Chapter 24: Promoter Design Artifact Slice (Offline Synthetic TP73 Locus)](./08-03_promoter_design_artifact_slice_offline.md) first.

## Parameters That Matter

- `upstream_bp=40 / downstream_bp=10` (where used: BuildGeneSetPromoterCohort, SummarizePromoterCohortComparison, ResolveOrthologPromoterCohort)
  - Why it matters: The synthetic chromosomes are only a few hundred bases long; compact promoter windows make the TSS geometry and motif contrast inspectable.
  - How to derive it: Use the fixed tutorial TSS at 151 bp, which yields promoter windows spanning 111..161 on the plus strand.
- `relationship=co_regulated` (where used: Gene-set promoter cohort, same-genome promoter comparison, ortholog promoter cohort, ortholog promoter comparison)
  - Why it matters: The relationship is a declared expectation carried into reports; it does not itself assert regulation.
  - How to derive it: Use `co_regulated` here because the fixture is deliberately arranged to show what an unexpected divergence flag looks like.
- `motif=SP1` (where used: Same-genome and ortholog promoter comparisons)
  - Why it matters: The synthetic TP73 promoter contains one GC-rich SP1-like patch while the comparison promoters are T-rich.
  - How to derive it: Use the exact local fixture sequences in `docs/examples/assets/promoter_cohort_tutorial_*.fa`.

## When This Routine Is Useful

- You want a deterministic offline smoke test for gene-set promoter cohort artifacts.
- You want to see how a declared co-regulated expectation is carried separately from evidence-derived relationship flags.
- You want a tiny local ortholog-resource example before using a larger reviewed ortholog table.
- You need generated JSON artifacts that can be inspected by GUI, CLI, MCP, or agent consumers without network access.

## What You Learn

- Distinguish a resolved gene-set promoter cohort from a promoter TFBS comparison over that cohort.
- Read `relationship` as a declared expectation and `relationship_flags[]` as non-blocking evidence triage.
- Use a local `gentle.ortholog_resource.v1` fixture for offline ortholog promoter resolution.
- Keep cross-species promoter evidence conservative: association evidence, not proof of regulation.
- Verify that generated tutorial artifacts are deterministic and reviewable.

## Applied Concepts

- **Shared Engine Contract** (`shared_engine_contract`): GUI, CLI, shell, and scripting interfaces execute the same operation semantics.
- **Deterministic Workflows** (`deterministic_workflows`): Operation chains should produce stable IDs and comparable outputs across repeated runs.
- **Genome Catalog Targeting** (`genome_catalog_targeting`): Prepared genome catalogs, annotation-based gene filters, and anchor extension connect imported entries to genomic context.
- **Promoter Motif Controls** (`promoter_motif_controls`): Foreground promoter motif signals should be compared with matched controls before being treated as candidate enrichment, depletion, or co-occurrence evidence.
- **Tutorial Drift Checks** (`tutorial_drift_checks`): Tutorial content is generated from executable examples and verified in automated checks.

## At a Glance

1. Open the generated artifact directory for this chapter after running tutorial...
2. Inspect promoter_cohort_tutorial.gene_set_promoter_cohort.json to see the res...
3. Inspect promoter_cohort_tutorial.gene_set_promoter_comparison.json to see the...
4. Inspect promoter_cohort_tutorial.ortholog_promoter_cohort.json to compare hum...
5. Inspect promoter_cohort_tutorial.ortholog_promoter_comparison.json to see cro...

## GUI First

CLI snippets use GENtle's default `.gentle_state.json` state unless they say otherwise. Add `--state PATH` or `--project PATH` when you want an explicit sandboxed state file for copied commands.

### Step 1: Open the generated artifact directory for this chapter after running tutorial...

GUI: Open the generated artifact directory for this chapter after running tutorial generation.

CLI:

```bash
cargo run --bin gentle_cli -- workflow @docs/examples/workflows/gene_set_ortholog_promoter_cohorts_offline.json
```

> Expected: The canonical workflow prepares `TutorialHumanPromoterToy` and `TutorialMousePromoterToy` from local FASTA/GTF fixtures.

### Step 2: Inspect promoter_cohort_tutorial.gene_set_promoter_cohort.json to see the res...

GUI: Inspect `promoter_cohort_tutorial.gene_set_promoter_cohort.json` to see the resolved gene-set members and promoter windows.

CLI:

```bash
cargo run --bin gentle_cli -- shell 'gene-sets promoter-cohort TutorialHumanPromoterToy --group tutorial_p73_promoter_cohort --relationship co-regulated --upstream-bp 40 --downstream-bp 10 --gene-group-catalog docs/examples/assets/promoter_cohort_tutorial_gene_groups.json --genome-catalog docs/examples/assets/promoter_cohort_tutorial_genomes.json --cache-dir /tmp/gentle-promoter-cohort-cache --path /tmp/gene_set_promoter_cohort.json'
```

> Expected: `gene_set_promoter_cohort.json` resolves two reviewed gene-set members and records `relationship: co_regulated`.

### Step 3: Inspect promoter_cohort_tutorial.gene_set_promoter_comparison.json to see the...

GUI: Inspect `promoter_cohort_tutorial.gene_set_promoter_comparison.json` to see the co-regulated expectation and the `unexpected_divergence` relationship flag.

CLI:

```bash
cargo run --bin gentle_cli -- shell 'genomes promoter-cohort-comparison TutorialHumanPromoterToy --cohort-label tutorial_p73_gene_set_co_regulated --cohort-kind co_regulated --gene TP73 --gene TP73D --motif SP1 --upstream-bp 40 --downstream-bp 10 --score-kind llr_background_tail_log10 --catalog docs/examples/assets/promoter_cohort_tutorial_genomes.json --cache-dir /tmp/gentle-promoter-cohort-cache --path /tmp/gene_set_promoter_comparison.json'
```

> Expected: `gene_set_promoter_comparison.json` contains a non-empty `relationship_flags[]` row with `flag_kind: unexpected_divergence`.

### Step 4: Inspect promoter_cohort_tutorial.ortholog_promoter_cohort.json to compare hum...

GUI: Inspect `promoter_cohort_tutorial.ortholog_promoter_cohort.json` to compare human TP73 and mouse Trp73 promoter window geometry.

CLI:

```bash
cargo run --bin gentle_cli -- shell 'orthologs resolve-promoter-cohort --anchor-species human --anchor-genome TutorialHumanPromoterToy --anchor-gene TP73 --target-species mouse --target-genome mouse=TutorialMousePromoterToy --transcript human=TX_HUMAN_TP73 --transcript mouse=TX_MOUSE_TRP73 --orthologs docs/examples/assets/promoter_cohort_tutorial_orthologs.json --relationship co-regulated --upstream-bp 40 --downstream-bp 10 --catalog docs/examples/assets/promoter_cohort_tutorial_genomes.json --cache-dir /tmp/gentle-promoter-cohort-cache --path /tmp/ortholog_promoter_cohort.json'
```

> Expected: `ortholog_promoter_cohort.json` resolves human TP73 and mouse Trp73 promoter windows from the local ortholog resource.

### Step 5: Inspect promoter_cohort_tutorial.ortholog_promoter_comparison.json to see cro...

GUI: Inspect `promoter_cohort_tutorial.ortholog_promoter_comparison.json` to see cross-species TFBS evidence kept separate from the relationship flag.

CLI:

```bash
cargo run --bin gentle_cli -- shell 'orthologs promoter-comparison --cohort /tmp/ortholog_promoter_cohort.json --motif SP1 --score-kind llr_background_tail_log10 --relationship co-regulated --path /tmp/ortholog_promoter_comparison.json'
```

> Expected: `ortholog_promoter_comparison.json` contains the ortholog `relationship: co_regulated` expectation and an `unexpected_divergence` relationship flag.


## Follow-up Commands

```bash
cargo run --bin gentle_examples_docs -- tutorial-generate
cargo run --bin gentle_examples_docs -- tutorial-check
cargo run --bin gentle_cli -- workflow @docs/examples/workflows/gene_set_ortholog_promoter_cohorts_offline.json
```

## Checkpoints

- The workflow executes offline with no genome download or ortholog API call.
- The gene-set cohort artifact reports exactly two returned promoter windows.
- The same-genome promoter comparison reports at least one `unexpected_divergence` flag under the co-regulated expectation.
- The ortholog cohort artifact reports two resolved promoter rows and zero unresolved targets.
- The ortholog comparison artifact reports at least one `unexpected_divergence` relationship flag and keeps the wording evidence-based.

## What This Chapter Produces

- [`artifacts/gene_set_ortholog_promoter_cohorts_offline/artifacts/promoter_cohort_tutorial.gene_set_promoter_cohort.json`](../artifacts/gene_set_ortholog_promoter_cohorts_offline/artifacts/promoter_cohort_tutorial.gene_set_promoter_cohort.json) - schema: `gentle.gene_set_promoter_cohort.v1`
- [`artifacts/gene_set_ortholog_promoter_cohorts_offline/artifacts/promoter_cohort_tutorial.gene_set_promoter_comparison.json`](../artifacts/gene_set_ortholog_promoter_cohorts_offline/artifacts/promoter_cohort_tutorial.gene_set_promoter_comparison.json) - schema: `gentle.promoter_cohort_comparison.v1`
- [`artifacts/gene_set_ortholog_promoter_cohorts_offline/artifacts/promoter_cohort_tutorial.ortholog_promoter_cohort.json`](../artifacts/gene_set_ortholog_promoter_cohorts_offline/artifacts/promoter_cohort_tutorial.ortholog_promoter_cohort.json) - schema: `gentle.ortholog_promoter_cohort.v1`
- [`artifacts/gene_set_ortholog_promoter_cohorts_offline/artifacts/promoter_cohort_tutorial.ortholog_promoter_comparison.json`](../artifacts/gene_set_ortholog_promoter_cohorts_offline/artifacts/promoter_cohort_tutorial.ortholog_promoter_comparison.json) - schema: `gentle.ortholog_promoter_comparison.v1`

## Tutorial Provenance

- Chapter id: `gene_set_ortholog_promoter_cohorts_offline`
- Tier: `core`
- Example id: `gene_set_ortholog_promoter_cohorts_offline`
- Tutorial source JSON: `docs/tutorial/sources/08-08_gene_set_ortholog_promoter_cohorts_offline.json`
- Workflow file: `docs/examples/workflows/gene_set_ortholog_promoter_cohorts_offline.json`
- Generated artifact dir: `docs/tutorial/generated/artifacts/gene_set_ortholog_promoter_cohorts_offline`
- Example test_mode: `always`
- Executed during generation: `yes`
- Automated status: `passing`
- Review status: `codex_reviewed`
- Codex reviewed at: `2026-06-25`
- Human reviewed at: `not recorded`
- Inspect the source JSON when you need full option-level detail.

## Feedback

If this tutorial is confusing, execution-stale, biologically suspect, or missing a useful figure, please open the matching tutorial issue template and include the context below.

- Tutorial title: `Gene-Set and Ortholog Promoter Cohorts (Offline Synthetic TP73 Demo)`
- Tutorial/chapter id: `gene_set_ortholog_promoter_cohorts_offline`
- Step reached:
- Expected vs. actual:
- Interface used: GUI / CLI / Agent Assistant / ClawBio

Paste the Tutorial feedback context here:

```text

```
