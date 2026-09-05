# Regulatory-fragment panel planner: Slice 1 design

Status: implemented for review; later evidence, presentation, GUI, and tutorial
slices are intentionally not part of this change.

## Purpose

The Slice 1 planner turns exact, persisted genomic regions into the smallest
bounded set of reporter constructs and pairwise contrasts needed to address a
declared experimental question. It supports one candidate region, an optional
partner, one required minimal-promoter region under the default policy, and an
optional reference control.

This is experimental-design support. A planned construct is a testable
candidate, not evidence that a region is sufficient, an enhancer, a silencer,
or partner-dependent.

## Baseline verification

The prompt was inspected at `6cb8f829`; implementation started after rebasing
onto `00cb5af2`. The relevant contracts had not drifted:

- `PromoterReporterPanelMemberProposal` still represents exactly one source
  interval per member.
- `PromoterReporterPanelFragmentRole` remains the closed `core`/`extended`
  vocabulary and was not widened.
- `PromoterReporterPanelMutationPolicy` remains the closed existing mutation
  vocabulary and is reused for the optional motif control.
- persisted regions still provide identity/content digests and the five
  `GenomicRegionLocalProjectionStatus` values used by the fail-closed gate.

The only baseline issue encountered was operational: a clean target directory
made the initial filtered test command spend more than thirteen minutes
rebuilding native macOS/GUI dependencies. It was stopped without a diagnostic;
subsequent incremental checks are the verification basis for this slice.

## Contracts

`gentle.regulatory_fragment_panel_request.v1` binds:

- the exact region-set id and content digest;
- the complete ROI snapshot, including identity/content digests, assembly,
  contig, 0-based half-open interval, strand, and local projection;
- an explicit annotation/reference release for cross-region compatibility;
- one exact catalog-owned vector and its validated multiple-cloning region;
- requested questions and exact opt-in order/orientation/spacing/boundary
  geometries;
- optional report/row evidence citations by exact digest;
- panel and construct bounds; and
- the complete effective policy, including defaults.

`gentle.regulatory_fragment_panel_plan.v1` records the normalized request,
source-state digest, vector context, selected constructs, omitted variants,
induced pairwise contrasts, exact insert sequences/digests, cloning
feasibility, eight independent evidence lanes, blockers, warnings, non-claims,
and an approval digest.

ROI identity and content remain distinct. `identity_sha256` identifies the
stable region geometry/provenance identity defined by the existing ROI
contract. `content_sha256` binds the complete current ROI record. The planner
requires both to match the persisted region and also binds the containing
region-set `content_sha256`.

## Validation order

Normalization and validation complete before candidate construction or
cloning simulation:

1. Validate schema, ids, bounds, effective policy, role cardinalities,
   questions, exact geometries, and evidence digests.
2. Resolve each bound set and ROI and compare exact identities/content.
3. Reject absent projections and every non-`current` stored or recomputed
   projection status by its exact status name.
4. Verify source-sequence existence and digest.
5. Reject mixed assemblies or releases; no liftover or substitution occurs.
6. Resolve the exact vector through the existing reporter catalog validation
   and bind its observed multiple-cloning-region context.
7. Only then generate candidates, contrasts, bounded selection, and detached
   cloning feasibility.

Role/purpose disagreement is retained as a named warning because purpose and
experimental role are independent evidence. It is never silently corrected.

## Deterministic selection

The planner generates only the closed base family and explicitly requested
variants:

- promoterless and minimal-promoter controls under policy;
- candidate alone and, when bound, partner alone;
- one exact declared A+B reference geometry;
- only the exact requested order, orientation, spacing, boundary, or tiling
  geometries;
- an optional exact reference control; and
- an optional motif edit produced solely by
  `p53_family_core_disruption_v1` at an exact fragment-local interval.

It does not form a Cartesian product. Duplicate exact insert geometries retain
the first declared candidate and receive one typed omission record.

Each requested question maps to exact member pairs. Selection exhaustively
searches the bounded candidate family and maximizes covered questions, then
minimizes member count. Ties follow normalized declared order, construct-role
rank, and member id. Required controls are retained when the bound permits.
When no complete cover fits, the plan stays within `max_panel_members`, uses
the same deterministic partial-order rule, and emits
`unresolved_experimental_comparison_required` with every uncovered question.

Every selected member has typed inclusion reasons. Every generated but omitted
member has a typed omission reason.

## Evidence and scientific language

Slice 1 preserves these independent lanes:

1. reference-genomic uniqueness;
2. panel/vector sequence similarity;
3. repeats and low complexity;
4. pair-specific junction uniqueness;
5. restriction and cloning risk;
6. Ensembl regulatory overlap;
7. TFBS/model-score context; and
8. CUT&RUN/chromatin context.

All are explicitly `not_evaluated` in this slice, including lanes with cited
upstream report identities. `not_evaluated` is never a pass. Candidate/partner
global similarity is used only to emit the conservative
`context_or_geometry_confounded` planning label.

## Approval and mutation boundary

Planning runs against live state but performs construct feasibility in a
detached engine. It does not mutate project state or write files. Its digest
binds the normalized request, effective policy, persisted source state,
evidence citations, vector context, selected member order, exact insert
sequences, contrasts, and planned operation payloads.

`validate_regulatory_fragment_panel_approval` checks the embedded plan digest,
the supplied approval digest, and a freshly recomputed plan against current
state. Approval is therefore exact and becomes stale when any bound input or
ordered output changes. It is authorization to proceed with a reviewed plan,
not validation of biological claims.

Materialization is deliberately unavailable in Slice 1. The existing
`PromoterReporterPanelMemberProposal` cannot encode ordered multi-fragment
instances, per-instance orientation, or spacers. The plan therefore reports
`materialization_supported=false` and carries exact insert-creation payloads
for audit, but does not route them through the legacy materializer or create a
parallel mutation path. A later reviewed transition must either add a
compatible materialization contract or explicitly restrict materialization to
the subset representable by the legacy model.

## Relationship to the study composer

The roadmap's regulatory-reporter study composer remains upstream. It may
translate a perturbation-response cohort and explicit TSS/evidence policies
into candidate ROI sets and request inputs. This planner starts only after the
user has persisted and exactly selected those regions; it does not discover
partners, infer TSSs, or rank biological activity. The two components therefore
compose without duplicating either candidate discovery or panel selection.

## Follow-on slices

- Slice 2: populate evidence lanes by composing existing engine reports, then
  add a shared machine-readable presentation/figure.
- Slice 3: add thin GUI and adapter entry points over the same operation.
- Slice 4: add the offline tutorial and broader user documentation.

