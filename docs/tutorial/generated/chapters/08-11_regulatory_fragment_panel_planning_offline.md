---
chapter_id: "regulatory_fragment_panel_planning_offline"
title: "Plan Exact A/B Regulatory-Fragment Contrasts (Offline)"
tier: "core"
example_id: "regulatory_fragment_panel_planning_offline"
source_example: "docs/examples/workflows/regulatory_fragment_panel_planning_offline.json"
example_test_mode: "always"
executed_during_generation: true
automated_status: "passing"
review_status: "codex_reviewed"
review_stale: false
codex_reviewed_at: "2026-09-08"
human_reviewed_at: null
human_reviewer: null
review_stale_reason: null
review_issue_template: null
review_issue_template_path: null
generated_artifact_dir: "docs/tutorial/generated/artifacts/regulatory_fragment_panel_planning_offline"
---

# Plan Exact A/B Regulatory-Fragment Contrasts (Offline)

Bind three current genomic regions to candidate, partner, and promoter-context roles; request only named A/B order, orientation, and spacing comparisons; inspect the deterministic minimal panel and its independent evidence lanes; and export a digest-valid review figure without materializing constructs.

A fragment cannot be called sufficient or partner-dependent from its sequence, annotation, or occupancy signal. Those claims require contrasts. GENtle therefore starts with exact, persisted genomic regions and asks which reporter members are needed to make each caller-declared question identifiable. Candidate A, partner B, and a promoter-context control are roles in that design. They are not inferred relationships.

This offline example uses a pinned public Ensembl-116 SERPINE1 locus only as a reproducible coordinate and sequence substrate. The three tutorial spans are hand-selected, and the reporter backbone is the repository-owned synthetic MCS fixture rather than a commercial vector. The request asks for A alone, B alone, the reference A+B geometry, one reversed-order construct, one reversed-orientation construct, one controlled-spacing construct, and promoterless/minimal-promoter controls. The resulting eight-member panel is the deterministic bounded cover of those questions, not a Cartesian expansion. Five sequence-computable evidence lanes are evaluated independently; Ensembl Regulation, TFBS, and CUT&RUN/chromatin lanes remain `not_evaluated` because a report citation alone is not verifiable report content.

**Prerequisites:** Read [Chapter 29: Promoter-Reporter Panel Planning (Offline Approval-Gated Demo)](./08-09_promoter_reporter_panel_planning_offline.md), [Chapter 30: Save and Share Genomic Regions (Offline SERPINE1 Example)](./08-10_portable_genomic_regions_offline.md) first.

## Parameters That Matter

- `candidate_alone / partner_alone and minimal_promoter` (where used: single-fragment test members and promoter-context control)
  - Why it matters: In this example, A alone means A followed by the declared promoter-context span, without B; B alone uses the same context without A. Neither label means promoter-independent activity. The promoterless control has an empty insert, while the minimal-promoter control contains only the declared context span.
  - How to derive it: Inspect each member's instances array. Treat the hand-selected context span as an instructional control, not a validated minimal promoter suitable for an experiment.
- `helper_catalog_path and output artifacts` (where used: workflow execution and proposal provenance)
  - Why it matters: The generated chapter retains the portable locus SVG. Each execution also writes plan JSON, and the final command renders its review SVG. Those review artifacts bind the resolved local helper-catalog path, so their digest can differ across checkouts; they are not committed as path-independent tutorial snapshots.
  - How to derive it: Run the workflow and renderer together in your checkout. Retain the resulting JSON and matching SVG for that execution; identical inputs and paths produce identical plan bytes.
- `region_set_content_sha256 and ROI identity/content digests` (where used: each fragments[] binding)
  - Why it matters: Identity pins the genomic region, while content also pins its current projection, evidence, notes, and presentation-bearing record. A changed set or ROI must not inherit an old review.
  - How to derive it: Select persisted regions from one exported `gentle.genomic_region_set.v1`; copy the engine-generated set and ROI digests rather than recomputing them by hand.
- `reference_combination.instances[]` (where used: A+B reference member)
  - Why it matters: The array is the exact 5'-to-3' insert order. Each `spacer_before` belongs to that instance, and the first instance must have an empty spacer.
  - How to derive it: Declare the intended experimental geometry explicitly; GENtle does not infer a partner, order, orientation, or gap from genomic proximity.
- `requested_variants` (where used: order, orientation, and spacing questions)
  - Why it matters: Only listed geometries are eligible. This prevents a silent Cartesian expansion and keeps the resulting panel reviewable.
  - How to derive it: Add one exact geometry for each scientific contrast you intend to test and omit combinations that do not answer a declared question.
- `max_panel_members=8` (where used: deterministic minimal-cover selection)
  - Why it matters: The bound is part of the request and digest. If it prevents complete question coverage, the plan names the uncovered questions instead of exceeding it.
  - How to derive it: Count mandatory controls and the smallest endpoint set for the requested contrasts, then choose an explicit practical upper bound.
- `proposal_digest` (where used: JSON/SVG review identity)
  - Why it matters: It binds normalized policy, source regions and sequences, vector identity, exact ordered members, evidence, and planned operations. It records what was reviewed, not whether the biology is valid.
  - How to derive it: Use the digest emitted by a fresh plan; never edit a plan document and retain its old digest.

## When This Routine Is Useful

- You have persisted genomic regions and want to turn a declared A/B hypothesis into explicit reporter contrasts.
- You need order, orientation, and spacer DNA recorded per construct rather than implied in prose.
- You want a bounded panel that names uncovered questions instead of silently dropping them or generating a Cartesian product.
- You want sequence similarity, repeats, junction uniqueness, cloning risk, and unavailable biological evidence kept in separate lanes.
- You need the same content-addressed plan in GUI and CLI, with exact-product creation kept behind a separate proposal and approval.

## What You Learn

- Distinguish a caller-declared experimental role from an inferred biological relationship.
- Explain why A alone, B alone, and A+B are needed before a partner-dependence question becomes experimentally identifiable.
- Read ordered fragment instances, reverse-complement orientation, and exact spacer DNA from the portable plan.
- Recognize deterministic minimal coverage as different from a full Cartesian construct library.
- Keep sequence, cloning, provider annotation, motif-model, and occupancy evidence in separate states and never turn `not_evaluated` into a pass.
- Use the same read-only operation through the Promoter design GUI and `promoters regulatory-panel-*` shell routes.
- Understand why panel review does not authorize biology or mutation: creating ordered multi-fragment products requires a separate exact-product proposal and approval, while cloning and functional validation remain separate.

## Applied Concepts

- **Shared Engine Contract** (`shared_engine_contract`): GUI, CLI, shell, and scripting interfaces execute the same operation semantics.
- **Deterministic Workflows** (`deterministic_workflows`): Operation chains should produce stable IDs and comparable outputs across repeated runs.
- **Genome Catalog Targeting** (`genome_catalog_targeting`): Prepared genome catalogs, annotation-based gene filters, and anchor extension connect imported entries to genomic context.
- **Artifact Exports** (`artifact_exports`): Representative outputs (CSV/protocol/SVG/text) are retained for auditability and sharing.
- **Tutorial Drift Checks** (`tutorial_drift_checks`): Tutorial content is generated from executable examples and verified in automated checks.

## At a Glance

1. Run the offline workflow once using a new state-file path and open that saved...
2. Open regulatory_panel_locus, choose Promoter design, and expand Regulatory-fr...
3. Enter loaded vector id regulatory_panel_vector, vector catalog id Synthetic p...
4. Inspect A, B, A+B, the three explicitly requested geometry variants, and both...
5. Expand the five sequence-computable evidence lanes. Treat evaluated as an ass...
6. Inspect Ensembl Regulation, TFBS/model-score, and CUT&RUN/chromatin lanes. Al...
7. Review the panel digest, blockers, and non-claims. This walkthrough stops at ...
8. Export JSON and SVG. The SVG is a view of the same digest-valid plan; changin...

## GUI First

CLI snippets use GENtle's default `.gentle_state.json` state unless they say otherwise. Add `--state PATH` or `--project PATH` when you want an explicit sandboxed state file for copied commands.

### Step 1: Run the offline workflow once using a new state-file path and open that saved...

GUI: Run the offline workflow once using a new state-file path and open that saved project. It contains both the pinned locus and the synthetic vector. Reusing an already imported region-set id is deliberately rejected; use a different state file for a fresh run.

CLI:

```bash
cargo run --bin gentle_cli -- --state /tmp/gentle-regulatory-fragment-tutorial.json workflow @docs/examples/workflows/regulatory_fragment_panel_planning_offline.json
```

> Expected: The workflow resolves the pinned locus, imports the content-addressed region set, validates the exact synthetic vector, and plans entirely offline.

### Step 2: Open regulatory_panel_locus, choose Promoter design, and expand Regulatory-fr...

GUI: Open `regulatory_panel_locus`, choose `Promoter design`, and expand `Regulatory-fragment panel`. Select `Candidate A tutorial span`, `Partner B tutorial span`, and `Minimal promoter-context tutorial span` for their respective roles; leave reference control empty. The selectors show coordinates and `Current` projection status; exact region-set and ROI digests are available in exported JSON.

CLI:

```bash
jq '.request.fragments[] | {fragment_id,role,reference_release,region_id:.region.region_id,interval:.region.interval,projection:.region.local_projection,region_set_content_sha256}' artifacts/regulatory_fragment_panel.plan.json
```

> Expected: All three bound regions retain GRCh38 coordinates, Ensembl-116 release context, source-sequence digests, and `current` local projections.

![The pinned Ensembl-116 SERPINE1 locus used as a coordinate substrate, with candidate A, caller-declared partner B, and the promoter-context control retained as exact genomic regions. The tutorial roles define comparisons; they do not establish regulatory activity or partnership.](../artifacts/regulatory_fragment_panel_planning_offline/artifacts/regulatory_fragment_panel.locus.svg)

*Figure: The pinned Ensembl-116 SERPINE1 locus used as a coordinate substrate, with candidate A, caller-declared partner B, and the promoter-context control retained as exact genomic regions. The tutorial roles define comparisons; they do not establish regulatory activity or partnership. Regenerate with `cargo run --bin gentle_examples_docs -- tutorial-generate`.*

> SVG text labels: `SERPINE1 locus evidence | 5' -> 3' display, genomic 101121158 -> 101139263 | strand + | assembly GRCh38 | annotation Ensembl 116 pinned public fixture | upstream 1000 bp | downs...`. If the embedded preview omits text in the GUI, open the linked SVG or use these labels as the figure legend.

### Step 3: Enter loaded vector id regulatory_panel_vector, vector catalog id Synthetic p...

GUI: Enter loaded vector id `regulatory_panel_vector`, vector catalog id `Synthetic panel vector`, helper catalog path `docs/examples/assets/promoter_reporter_panel_demo_helper_vectors.json`, and reference release `Ensembl 116 pinned public fixture`. Enable partner, order, orientation, and spacing comparisons, enter spacer `GCGCGC`, set the member bound to eight and insert-length bound to 1,000 bp, then click `Plan exact panel`. Confirm eight members and no uncovered questions.

CLI:

```bash
jq '{plan_id,proposal_digest,planning_label,member_count:(.members|length),contrast_count:(.contrasts|length),uncovered_questions,materialization_supported}' artifacts/regulatory_fragment_panel.plan.json
```

> Expected: The planner returns exactly eight members and covers all five requested questions without generating unrequested construct combinations.

### Step 4: Inspect A, B, A+B, the three explicitly requested geometry variants, and both...

GUI: Inspect A, B, A+B, the three explicitly requested geometry variants, and both controls. Only the controlled-spacing member should have `GCGCGC` before B; the other geometries have no spacer. GUI-generated member names and request digests may differ from the workflow's names, but these comparison geometries should agree.

CLI:

```bash
jq '.members[] | {member_id,construct_kind,insert_length_bp,instances:[.instances[]|{fragment_id,orientation,spacer_before,assembled_start_0based,assembled_end_0based_exclusive}],contrast_ids}' artifacts/regulatory_fragment_panel.plan.json
```

> Expected: Every ordered instance records genomic and assembled coordinates, orientation, and exact spacer-before DNA.

### Step 5: Expand the five sequence-computable evidence lanes. Treat evaluated as an ass...

GUI: Expand the five sequence-computable evidence lanes. Treat `evaluated` as an assessment state, not a pass. Reference uniqueness here concerns the loaded, ROI-bound locus only, not the complete genome.

CLI:

```bash
jq '.evidence_dimensions[] | {kind,state,observations:(.observations|length),blockers,warnings}' artifacts/regulatory_fragment_panel.plan.json
```

> Expected: Reference uniqueness within the loaded locus, panel/vector similarity, repeats/low complexity, junction uniqueness, and cloning risk are evaluated as independent lanes; no genome-wide uniqueness claim is made.

### Step 6: Inspect Ensembl Regulation, TFBS/model-score, and CUT&RUN/chromatin lanes. Al...

GUI: Inspect Ensembl Regulation, TFBS/model-score, and CUT&RUN/chromatin lanes. All three remain `not_evaluated`; the imported coordinates are not a substitute for typed biological evidence reports.

CLI:

```bash
jq '.evidence_dimensions[] | select(.kind=="ensembl_regulatory_overlap" or .kind=="tfbs_model_score_context" or .kind=="cutrun_and_chromatin_context") | {kind,state,detail}' artifacts/regulatory_fragment_panel.plan.json
```

> Expected: Ensembl Regulation, TFBS/model-score, and CUT&RUN/chromatin lanes are explicitly `not_evaluated`; this state is not a pass or evidence of absence.

### Step 7: Review the panel digest, blockers, and non-claims. This walkthrough stops at ...

GUI: Review the panel digest, blockers, and non-claims. This walkthrough stops at review. The separate `Exact design products` section can prepare full product sequences for another review, but creation requires that product proposal's own digest, not the panel digest. A designed molecule is not a validated cloning reaction.

CLI:

```bash
jq '{proposal_digest,approval_required,materialization_supported,blockers,warnings,nonclaims}' artifacts/regulatory_fragment_panel.plan.json
```

> Expected: The plan requires review under its exact digest and advertises a separate exact-product proposal via `materialization_supported=true`; planning alone adds no construct or primer to project state.

### Step 8: Export JSON and SVG. The SVG is a view of the same digest-valid plan; changin...

GUI: Export JSON and SVG. The SVG is a view of the same digest-valid plan; changing a bound region, sequence, vector, geometry, or member order requires replanning. Keep each plan with its own figure rather than substituting the CLI workflow's digest for a GUI-authored request.

CLI:

```bash
cargo run --bin gentle_cli -- --state /tmp/gentle-regulatory-fragment-tutorial.json shell 'promoters regulatory-panel-render @artifacts/regulatory_fragment_panel.plan.json --path artifacts/regulatory_fragment_panel.plan.svg'
```

> Expected: The renderer accepts the unchanged plan and refuses edited plan content that no longer matches its proposal digest.


## Follow-up Commands

```bash
cargo run --bin gentle_cli -- --state /tmp/gentle-regulatory-fragment-tutorial.json shell 'promoters regulatory-panel-render @artifacts/regulatory_fragment_panel.plan.json --path artifacts/regulatory_fragment_panel.plan.svg'
cargo run --bin gentle_examples_docs -- tutorial-generate
cargo run --bin gentle_examples_docs -- tutorial-check
```

## Checkpoints

- The retained locus SVG and imported region-set asset are deterministic projections of the pinned offline sequence fixture.
- The generated plan contains exactly eight members, no uncovered question, and no Cartesian expansion.
- The reference A+B member records candidate A, partner B, and the promoter-context control with no implicit spacer DNA.
- The controlled-spacing member alone records `GCGCGC` before partner B; that spacer is not propagated into the other geometries.
- Only the requested reverse-order, reverse-partner-orientation, and six-base-spacing variants appear.
- Every evidence dimension remains separately inspectable; the three opaque external-report dimensions stay `not_evaluated`.
- Planning writes a review artifact but performs no materialization; no result text claims sufficiency, enhancer/silencer activity, or partner dependence.

## What This Chapter Produces

- [`artifacts/regulatory_fragment_panel_planning_offline/artifacts/regulatory_fragment_panel.locus.svg`](../artifacts/regulatory_fragment_panel_planning_offline/artifacts/regulatory_fragment_panel.locus.svg)

  - Embedded above near Step 2; kept here as an audit link.

> SVG text labels: `SERPINE1 locus evidence | 5' -> 3' display, genomic 101121158 -> 101139263 | strand + | assembly GRCh38 | annotation Ensembl 116 pinned public fixture | upstream 1000 bp | downs...`. If this embedded preview omits text in the GUI, open the linked SVG or use these labels as the figure legend.


## Tutorial Provenance

- Chapter id: `regulatory_fragment_panel_planning_offline`
- Tier: `core`
- Example id: `regulatory_fragment_panel_planning_offline`
- Tutorial source JSON: `docs/tutorial/sources/08-11_regulatory_fragment_panel_planning_offline.json`
- Workflow file: `docs/examples/workflows/regulatory_fragment_panel_planning_offline.json`
- Generated artifact dir: `docs/tutorial/generated/artifacts/regulatory_fragment_panel_planning_offline`
- Example test_mode: `always`
- Executed during generation: `yes`
- Automated status: `passing`
- Review status: `codex_reviewed`
- Codex reviewed at: `2026-09-08`
- Human reviewed at: `not recorded`
- Inspect the source JSON when you need full option-level detail.

## Feedback

If this tutorial is confusing, execution-stale, biologically suspect, or missing a useful figure, please open the matching tutorial issue template and include the context below.

- Tutorial title: `Plan Exact A/B Regulatory-Fragment Contrasts (Offline)`
- Tutorial/chapter id: `regulatory_fragment_panel_planning_offline`
- Step reached:
- Expected vs. actual:
- Interface used: GUI / CLI / Agent Assistant / ClawBio

Paste the Tutorial feedback context here:

```text

```
