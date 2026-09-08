---
chapter_id: "region_homology_promoter_modules_offline"
title: "Conserved Blocks and Testable Promoter-Module Hypotheses"
tier: "advanced"
example_id: "region_homology_promoter_modules_offline"
source_example: "docs/examples/workflows/region_homology_promoter_modules_offline.json"
example_test_mode: "optional_blast"
executed_during_generation: false
automated_status: "skipped"
review_status: "codex_reviewed"
review_stale: false
codex_reviewed_at: "2026-09-07"
human_reviewed_at: null
human_reviewer: null
review_stale_reason: null
review_issue_template: null
review_issue_template_path: null
generated_artifact_dir: "docs/tutorial/generated/artifacts/region_homology_promoter_modules_offline"
---

# Conserved Blocks and Testable Promoter-Module Hypotheses

Compare one saved genomic region with validated local genomic indexes, inspect a query-referenced alignment without insertion columns, and turn explicitly selected evidence spans into conservative reporter-module hypotheses.

A conserved promoter segment can be a useful fragment candidate, but sequence similarity does not prove that the segment works alone. This chapter keeps three questions separate: whether an explicitly expected ortholog locus supports each query base, whether other species contain merely unassigned similarity, and whether the same genome contains competing copies. The synthetic target includes an inserted base; GENtle records it in provenance but does not add a display column, so every row remains aligned to the query. The Conservation workspace can then evaluate saved CUT&RUN, motif, provider-annotation, or reporter-candidate spans against the exact-support blocks. Its standalone, paired-context, repetitive, and insufficient outcomes are traceable design hypotheses, not regulatory verdicts.

**Prerequisites:** Read [Chapter 30: Save and Share Genomic Regions (Offline SERPINE1 Example)](./08-10_portable_genomic_regions_offline.md) first.

> **How to Run This Locally**
> This offline example still needs already-installed BLAST+ (makeblastdb, blastdbcmd and blastn). Nothing installs tools automatically. Documentation generation does not execute it; tutorial-check runs its starter workflow only when all three tools are available, otherwise reporting an explicit skip rather than BLAST acceptance. For the full walkthrough, run the companion command below after making BLAST+ available.

## Parameters That Matter

- `targets[].role and expected_loci[]` (where used: ScreenGenomicRegionHomology)
  - Why it matters: An expected ortholog label is accepted only when an accepted locus overlaps a caller-supplied, evidence-identified expected locus.
  - How to derive it: Use reviewed orthology resources or an explicit expected locus; leave unrelated cross-species searches unassigned.
- `min_identity_percent=70, min_alignment_length_bp=24` (where used: Synthetic homology screen)
  - Why it matters: These default thresholds retain the intended ortholog and 24 bp same-genome copy while excluding incidental 12 bp matches in the hand-crafted teaching sequences.
  - How to derive it: Choose thresholds before search and preserve them in the effective request; do not tune them after seeing a preferred locus.
- `max_hsps_per_target=1000` (where used: Synthetic homology screen)
  - Why it matters: The cap is a post-BLAST processing budget and never a `-max_target_seqs` completeness shortcut.
  - How to derive it: Use a bounded value suitable for the query and fail explicitly when highly repetitive output exceeds it.
- `selected_evidence_spans[]` (where used: AssessPromoterConservedModules)
  - Why it matters: These are caller-selected CUT&RUN, motif, annotation, or reporter spans that the proposed fragment should retain; they are not inferred from conservation.
  - How to derive it: Save reviewed evidence as portable regions or supply content-bound query-local spans with independent provenance.

## When This Routine Is Useful

- You want to inspect whether a saved promoter candidate contains cross-species exact-support blocks.
- You want same-genome repetition reported separately from cross-species conservation.
- You want a multiple-alignment-like display whose columns always remain query coordinates.
- You want an inspectable reason for testing one reporter block alone or together with a partner block.
- You need ordinary GENtle use to remain unaffected when optional genomic BLAST indexes are absent.

## What You Learn

- Distinguish explicit orthology evidence from BLAST similarity.
- Read a query-referenced projection that suppresses target insertion columns without losing insertion provenance.
- Interpret exact-support blocks separately by evidence class and available-genome denominator.
- Use same-genome repetition as an ambiguity signal rather than a functional off-target claim.
- Explain why a reporter-module assessment proposes experiments but cannot establish autonomous regulatory function.

## Applied Concepts

- **Portable Genomic Regions** (`portable_genomic_regions`): Assembly-bound regions retain one canonical coordinate and evidence identity across project, CLI, GUI, and exported representations.
- **Local BLAST Indexes** (`local_blast_indexes`): Validated local indexes are content-bound optional resources; their absence remains explicit without disabling unrelated work.
- **Query-Referenced Alignment** (`query_referenced_alignment`): The query defines every display column while target insertions remain separate structured provenance.
- **Explicit Orthology Evidence** (`explicit_orthology`): An ortholog label requires a declared, provenance-bearing expected locus and is never inferred from similarity rank alone.
- **Reporter-Module Hypotheses** (`reporter_module_hypotheses`): Conservation and independently selected evidence produce traceable fragment-testing hypotheses, never proof of autonomous regulatory function.
- **Deterministic Workflows** (`deterministic_workflows`): Operation chains should produce stable IDs and comparable outputs across repeated runs.

## At a Glance

1. Run the complete companion command shown below into an empty directory. It pr...
2. Open the DNA viewer, choose Regions..., and use Conservation... on synthetic_...
3. Choose Open report... and load homology_report.json from the companion output...
4. Expand Search request to inspect target IDs, required flags, orthology eviden...
5. Inspect target readiness and confirm that expected-ortholog, unassigned cross...
6. Inspect the query-referenced alignment: dots are exact bases, letters are sub...
7. Select only synthetic_occupancy_anchor, choose Assess selected evidence, and ...
8. Optionally save a selected conserved block as a new portable region. This exp...

## GUI First

CLI snippets use GENtle's default `.gentle_state.json` state unless they say otherwise. Add `--state PATH` or `--project PATH` when you want an explicit sandboxed state file for copied commands.

### Step 1: Run the complete companion command shown below into an empty directory. It pr...

GUI: Run the complete companion command shown below into an empty directory. It prepares tiny local indexes, saves the query and two independent evidence spans, then generates the homology SVG and four assessed module reports. Open its `tutorial.project.json` in GENtle.

CLI:

```bash
cargo build --locked --bin gentle_cli
```

> Expected: The workflow uses only hand-crafted local FASTA files and requires no network access; it never downloads or indexes an undeclared genome.

### Step 2: Open the DNA viewer, choose Regions..., and use Conservation... on synthetic_...

GUI: Open the DNA viewer, choose `Regions...`, and use `Conservation...` on `synthetic_promoter_query`.

CLI:

```bash
python3 docs/examples/run_region_homology_tutorial.py --gentle target/debug/gentle_cli --output /tmp/gentle-conservation-tutorial
```

> Expected: Every projected alignment row has exactly the query length even though one target contains an insertion.

### Step 3: Choose Open report... and load homology_report.json from the companion output...

GUI: Choose `Open report...` and load `homology_report.json` from the companion output directory; GENtle validates its digest and exact saved-region binding. `Import request...` loads that directory's `homology_request.json` without executing it.

CLI:

```bash
target/debug/gentle_cli --state /tmp/gentle-conservation-tutorial/report_only.project.json shell 'regions render-homology-svg @/tmp/gentle-conservation-tutorial/homology_report.json /tmp/gentle-conservation-tutorial/replayed.svg'
```

> Expected: The inserted target base is retained under `omitted_insertions` with query anchor, target coordinates, strand, and HSP identity.

### Step 4: Expand Search request to inspect target IDs, required flags, orthology eviden...

GUI: Expand `Search request` to inspect target IDs, required flags, orthology evidence and search thresholds. `Run local screen` submits these exact settings; empty targets instead select all validated local indexes. During BLAST, elapsed time and heartbeats remain visible, and `Cancel` stops the active search without publishing an incomplete report.

CLI:

```bash
target/debug/gentle_cli --state /tmp/gentle-conservation-tutorial/report_only.project.json shell 'promoters assess-conserved-modules @/tmp/gentle-conservation-tutorial/paired.request.json'
```

> Expected: Only the locus backed by `synthetic_declared_orthology` is labelled `expected_ortholog`; BLAST rank alone never establishes orthology.

### Step 5: Inspect target readiness and confirm that expected-ortholog, unassigned cross...

GUI: Inspect target readiness and confirm that expected-ortholog, unassigned cross-species, and same-genome evidence remain separate.

> Expected: Same-genome non-self similarity has its own loci, support blocks, coverage percentage, and ambiguity rule.

### Step 6: Inspect the query-referenced alignment: dots are exact bases, letters are sub...

GUI: Inspect the query-referenced alignment: dots are exact bases, letters are substitutions, dashes are target deletions, and omitted insertions remain in JSON rather than adding columns.

> Expected: Module outcomes retain all thresholds, evidence IDs, passed and failed rules, alternatives, and explicit non-claims.

### Step 7: Select only synthetic_occupancy_anchor, choose Assess selected evidence, and ...

GUI: Select only `synthetic_occupancy_anchor`, choose `Assess selected evidence`, and compare the standalone decision with `standalone.json`. Add `synthetic_left_motif` and reassess: the paired decision must cite a shared ortholog locus and compatible target-coordinate gaps, not merely proximity in the query. Select a block to jump to its virtualized alignment tile.

> Expected: The companion asserts standalone and paired results, an insufficient result for evidence in the deliberately substituted gap, and a repetitive result under an explicitly strict 1% policy. Its receipt records the exact commands, binary and artifact hashes. GENtle alone computes these outcomes. Report-only steps use an empty scratch project because their inputs are self-contained; the original tutorial project is not reopened or rewritten.

### Step 8: Optionally save a selected conserved block as a new portable region. This exp...

GUI: Optionally save a selected conserved block as a new portable region. This explicit mutation preserves the report digest and non-claims.

> Expected: Missing or unrequested same-genome evidence cannot pass uniqueness. Query and target gaps must both satisfy the declared bounds; exact spacing is the default unless a gap-difference tolerance is supplied.


## Command Equivalent (After GUI)

Run the same routine non-interactively once the GUI flow is clear:

```bash
cargo run --bin gentle_cli -- workflow @docs/examples/workflows/region_homology_promoter_modules_offline.json
cargo run --bin gentle_cli -- shell 'workflow @docs/examples/workflows/region_homology_promoter_modules_offline.json'
```

## Follow-up Commands

```bash
python3 docs/examples/run_region_homology_tutorial.py --gentle target/debug/gentle_cli --output /tmp/gentle-conservation-tutorial-repeat
```

## Checkpoints

- The homology report binds query region identity, query sequence SHA-256, effective request, tool version, and every genomic database content fingerprint.
- All target rows are query-length projections; insertions are omitted from display columns and retained structurally.
- Unavailable optional indexes remain typed unavailable and do not impair unrelated GENtle use.
- The module decision trace distinguishes standalone, paired-context, repetitive-or-ambiguous, and insufficient-evidence hypotheses.
- Every presentation repeats that conservation does not prove autonomous promoter function.

## Tutorial Provenance

- Chapter id: `region_homology_promoter_modules_offline`
- Tier: `advanced`
- Example id: `region_homology_promoter_modules_offline`
- Tutorial source JSON: `docs/tutorial/sources/08-12_region_homology_promoter_modules_offline.json`
- Workflow file: `docs/examples/workflows/region_homology_promoter_modules_offline.json`
- Generated artifact dir: `docs/tutorial/generated/artifacts/region_homology_promoter_modules_offline`
- Example test_mode: `optional_blast`
- Executed during generation: `no`
- Automated status: `skipped`
- Review status: `codex_reviewed`
- Codex reviewed at: `2026-09-07`
- Human reviewed at: `not recorded`
- Inspect the source JSON when you need full option-level detail.

## Feedback

If this tutorial is confusing, execution-stale, biologically suspect, or missing a useful figure, please open the matching tutorial issue template and include the context below.

- Tutorial title: `Conserved Blocks and Testable Promoter-Module Hypotheses`
- Tutorial/chapter id: `region_homology_promoter_modules_offline`
- Step reached:
- Expected vs. actual:
- Interface used: GUI / CLI / Agent Assistant / ClawBio

Paste the Tutorial feedback context here:

```text

```
