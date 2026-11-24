# GENtle Roadmap

Last session: 2026-11-24. Obtain required updated-lock CI/platform checks;
manual native Windows GUI acceptance is non-blocking. Do not repeat the scoped
macOS, termination VM or process checks.
See the [Glen handoff](glen_next_verification_handoff.md).

Purpose: fast session orientation, followed by optional phase and parking-lot
navigation. The opening release/priorities overview should take about two
minutes; the whole file has an advisory 2,000-word budget. Completed work belongs
in [CHANGELOG.md](CHANGELOG.md), durable rules in [decisions.md](decisions.md),
and contracts in [protocol.md](protocol.md).

## Release Gate

Published baseline: `v0.1.0-internal.12`; retain its
[publication ledger](release_notes/release_notes_v0.1.0-internal.12.md#publication-evidence).
Current candidate: `v0.1.0-internal.13` (development version only), with **no
selected release SHA or publication**. Development continues without interim freezes.

At actual release preparation, complete the
[nine-gate acceptance ledger](release_notes/release_notes_v0.1.0-internal.13.md#internal13-acceptance-ledger):

1. Select one clean, pushed candidate SHA, lockfile and profile; never combine
   different revisions' passes.
2. Obtain green Windows, macOS and full Linux CI at that SHA, including locked
   checks, tests and canonical tutorial/parity freshness.
3. Build all three native packages with receipts and extracted smokes; retain
   the headless container identity/receipt and no-network CLI/MCP smokes.
4. Complete rollover/version/tutorial/parity and the A/B/C functional/native
   checks below, then obtain Glen's scoped readiness verdict.
5. Obtain explicit owner tag/publication approval. After confirmed release,
   perform the next development-version rollover.

Use the [owner procedure](internal_13_preparation_handoff.md#owner-procedure)
and [release process](release.md): `package-opt2` first; select the separate
opt1 fallback only after a retained hosted failure, on a new candidate.
These instructions authorize no push, dispatch, tag, build or publication.

Only manual Windows GUI acceptance is non-blocking. Windows/macOS/Linux CI and
all three native packages remain required. See the [owner clarification](decisions.md#windows-gui-acceptance-is-not-a-release-gate).

Historical `.10`/`.11` scientific ledgers, native latency, real-reference import
performance, profile comparisons and live-agent research carry forward without
becoming new `.13` blockers. Older passes never certify newer source; headless
checks do not replace graphical acceptance. Keep detailed scope distinctions in
the [verification navigation](glen_next_verification_handoff.md#roadmap-detail-navigation-2026-11-20).
The [historical CI correction](https://github.com/smoe/gentle_rs/blob/7c52082c24fca85817d70956146623eef6f24607/docs/internal_12_candidate_preparation.md#corrected-ci-evidence)
remains source-bound, not a current-source pass.

## Next Session Priorities

1. Obtain the updated lock's required cross-platform acceptance. Keep manual
   Windows GUI/focus checks non-blocking; do not replay completed focused tests.
   [Remaining checks](glen_next_verification_handoff.md#native-pcr-safeguards-and-viewports-2026-11-22).
2. Prepare A/B/C release evidence only when a candidate is selected; retain
   historical package and scientific receipts without reissuing completed work.
   [Acceptance ledger](release_notes/release_notes_v0.1.0-internal.13.md#internal13-acceptance-ledger).
3. Carry only remaining read-end/exon-locus, C1-C4, RNA and native-platform gaps;
   do not repeat the completed reader, barcode, VM or bound-artifact campaigns.
   [Detailed boundaries](glen_next_verification_handoff.md#roadmap-detail-navigation-2026-11-20).
4. Complete TSS/tutorial native inspection, source-coherent exports and
   live-host acceptance without treating synthetic geometry as scientific proof.
   [Native checklist](glen_final_candidate_execution_handoff.md#remaining-native-and-human-checks).
5. Verify responsive jobs, cancellation and lifecycle save/preview guards;
   general queues and durable cross-process admission remain separate work.
   [Migration plan](asynchronous_command_execution_plan.md#reconciled-codex-plan-2026-09-14).
6. Review authentic-reference primer-PAIR studies and specificity imports,
   retaining ambiguity, source bindings and separate order readiness.
   [Specialist scope](primer_design_specialist_plan.md).

## .13 Priorities

Owner-selected aim, 2026-10-15: minimal working paths, not feature breadth.
Use the ledger above for candidate-bound acceptance and responsibility.

- A: review early containers, pools/aliquots, arrangements/gels, racks/plates
  and 3D exports. Check fresh GENtle-exported SCAD, nine legacy baselines and
  native controls; template renders do not establish physical fit.
  [Labware scope](../data/labware/README.md).
- B: verify provider-neutral protein-expression drafts and tutorial 09.03:
  explicit provider choice, unanswered inputs and no submission/order/contact.
  Inspector refinements are stretch-only; ranking is later work.
  [Expression workflow](tutorial/09-03_protein_expression_handoff_cli.md).
- C: verify explicit complete-set/source-tube binding, one serial/default
  storage rack, provenance, save/reopen, undo/redo and fresh exports; no implicit
  placement or new freezer mapper.
  [Storage workflow](tutorial/03-10_gene_set_storage_projection.md).

The first official `0.1` direction remains the separately owner-earmarked
[Glen release aims](official_0_1_release_aims.md), not an expansion of `.13`.

## Phase A: AI Communication And Safety Plane

Keep GUI, CLI, MCP, JS/Lua and ClawBio/OpenClaw thin over shared engine/shell
contracts. Prioritize reachability, fact-annotated discovery, explicit approvals
and typed UI intents, not tool counts or autonomous biology.
See [agent contracts and exclusions](agent_interface.md).

Validate descriptor/runtime/tutorial parity before live model-quality review.
New skill activation, alias routing, commercial APIs and credentials require
separate contracts; provider handoffs are not submission or procurement.
See [capability maturity](capability_maturity.md).

Verify the shared expression inspector and explicit provider gates before later
chassis/vector/provider ranking; execute actual JS/Lua VMs rather than infer
acceptance from compilation. Check catalog imports and tutorial 09.04 without
inventing commercial offers or calibrated timing.
See [external-service contracts](catalog_updates.md).

## Phase B: Cloning Routine Standardization

Extend routine-family preflight, replayable macros and protocol packs only
within the selected story. Preserve portable records and candidate auditability;
do not add biology in adapters or silently optimize constructs.
See [routine catalog plan](cloning_routine_catalog_plan.md).

## Phase C: Engine And Protocol Parity Hardening

Harden shared file, sequencing, primer, projection and RNA contracts with
targeted tests; defer broad crate surgery. Separate strict engine/reachability
parity from optional first-class surfacing.
See [parity policy](architecture.md#1-product-intent).

Follow the [inline-operand audit](inline_operand_audit.md) order rather than
adding one-off adapters. Keep Clariom/probe/transcript evidence ambiguity,
provenance and unknown states; overlap/intensity are not isoform validation or
antibody compatibility.

### Reserved For Anze: Repeat/Similarity Inspection Follow-Up

Extend collection inspection only through portable, objective-aware evidence
records; curated repeat families and broader sequence-set mappings stay scoped.

## Phase D: Visualization And Workflow UX

Verify offline manual navigation/contextual Help, dense DNA/splicing inspection,
scroll/zoom, deterministic exports and context-click actions on real hosts.
No unrelated redesign; source tests and paint tests are not native acceptance.
See [manual navigation](manual_navigation_plan.md).

Retain microarray captures as attributed raw evidence; high-contrast teaching
derivatives and parent-overview polish remain separate. Keep cryptic-splicing
v2, persisted report selectors and richer provider GUI refinements deferred.
See [contextual availability](contextual_action_availability_plan.md).

## Phase E: Integration Polish And Deferred Policy Items

Reduce distribution risk through reproducible native/container receipts,
headless MCP smokes and truthful installation docs. Capture/viewport consent,
translations and small-window fonts require native review; no new screenshots
are authorized by this roadmap.
See [release procedures](release.md).

## Phase F: Interpretation Later

Keep interpretation inspectable, bounded and separately approved before
mutation. For read-end/exon-locus and C1-C4, verify only unproved native/data/
release boundaries; preserve original hypotheses, counts, bindings and caps.
See [remaining evidence navigation](glen_next_verification_handoff.md#roadmap-detail-navigation-2026-11-20).

Image-to-state proposals, richer findings, cohort reasoning and wet-lab digital
twins are later work. Rendered geometry is not physical fit or biological
validity; no private archive regeneration or new data run is authorized.
See [robotic simulation plan](robotic_wetlab_simulation_plan.md).

## Parking Lot

These are unselected follow-ups, not release commitments or execution authority.
Existing implementations retain their own verification ledgers.

- Sequence-design GUI, codon objectives and a verified standalone compiler
  contract: [scope](dna_sequence_optimization_plan.md).
- TSS original-input rescoring, shared-scale/PDF checks and source acceptance:
  [profiles](tss_tfbs_profiles.md#rescoring-existing-reports).
- A bounded representative Simple-PCR native smoke preset.
- Stream/index UniProt transcript inventory and bound finalization progress.
- Extend semantic GUI IDs and candidate tutorial gates; review Gibson/dotplot.
- Provenance-complete array evidence and approved transcript-study execution.
- Readiness-bound agent tutorial actions and explicit receipt replay.
- Optional OS credential storage and a bounded quota-generation probe.
- Separate opt-in coding-agent mode with workspace scope, diff review and
  rollback, never inherited Assistant execution authority.
- Opt-in Mistral conformance without implicit model calls.
- Restriction cut-offset/usage audit, detail popovers and linear-SVG clipping.
- Portable report/finding inspection and exhaustive operation-effect metadata.
- Larger cohort/locus navigation without collapsing reads or hard-coded genes.
- Explicit observed-junction `relative_to`; separate spliced-view conflicts.
- Vendor protocols, scheduling, quantities and dated offers beyond synthetic
  time/cost breakdowns.
- Cross-gene homology anchors, extensible gene groups and ontology bridges.
- Provider/CRO integrations and redistributable Luebeck vector resources.
- VKORC1/Factor X companion after constructs, controls and readouts are agreed.
- Bind SERPINE1 comparisons to actual insert/vector and CUT&RUN manifests.
- Primer walking, contig management and helper-construct terminology cleanup.
- CUT&RUN controls, replicates, calibrated differentials and whole-genome input.
- Quantitative Ensembl/ENCODE lanes with assay, units and context retained.
- Deeper Primer3/preflight, virtual PCR, multiplex, LAMP and allele-specific
  assays: [audit](primer3_parity_audit_20261103.md).
- Global RT-PCR redesign over genome/transcriptome evidence and explicit QC:
  [staged plan](rt_capture_primer_redesign_plan.md); supplied-pair checks follow
  the [compatibility plan](primer_blast_local_compatibility_plan_20261116.md#revised-minimal-commit-order).
- Advanced/Mutagenesis handoffs, persisted views and derived single-primer
  inspection: [residual scope](primer_design_specialist_plan.md).
- GuideRNA ranking, macro packaging and separately approved base/prime S2:
  [editing preparation](base_prime_editing_preparation.md).
- Lifecycle extensions beyond remaining Slice 3-6 verification: opaque wrappers,
  queues, durability and broader exporter safety remain outside those slices.
- Cross-tool comparison and separately assessed JS plain-JSON ingress.
- Measured gels: remaining native checks, coefficients/masks, scaled references,
  vendor traces, D5000, ratios and replicate inference; no enrichment claim.
  [Cross-ladder permission](decisions.md#cross-ladder-residual-policy) does not
  approve real coefficients or masks; see [gel boundaries](gel_workbench.md).
- Weekly/monthly maintenance automation, without creating a release treadmill.
- Local-first server/browser ownership; WebAssembly portability separately.
- More collection BLAST/arrangement/storage adapters over shared policies:
  [collection plan](gui_gene_set_collection_operations_plan.md).
- Gene-set enrichment with explicit universe, mappings, context, provenance,
  multiple-testing correction and unresolved counts; no causal claim.
- Reporter boundary/native-cloning QA with real peaks, controls and cohorts:
  [selector](reporter_fragment_selection.md).
- Regulatory-partner calibration, background cohorts, ontology and distal
  evidence, preserving the exact-hit audit layer.
- Paralog resolution only with a concrete consumer and relationship provenance.
- RNA panel throughput, queues/attachment and primer-payload calibration:
  [processing plan](rna_panel_gui_processing_plan.md).
- Observed-structure/J1 affected-source verification; J2/J3 stay owner-gated:
  [junction plan](observed_junction_consequences_plan.md).
- Synonymous restriction-site introduction with verification-digest proof.
- Report-first feature-library annotation; termination catalogs/calibration,
  circular-complete scans and dedicated GUI remain separate:
  [termination plan](transcription_termination_annotation_plan.md).
- Glossary/help inversion from typed descriptors:
  [remaining discovery work](introspection.md#12-deferred-follow-up).

## Scope And Session Close

- Do not add private proposal, grant, sequencing or unrelated paper/figure work.
- No broad engine extraction, GUI redesign or speculative protocol expansion
  without a selected task; keep business logic engine-owned.
- Do not add a `Done` entry. Move completed work directly to
  [CHANGELOG.md](CHANGELOG.md); durable rules belong in [decisions.md](decisions.md).
- Migrate detailed evidence before expanding this file; retain open requirements
  through [handoff navigation](glen_next_verification_handoff.md#roadmap-detail-navigation-2026-11-20).
- Run `scripts/maintenance_chore.py session-close --plan PATH_OR_LABEL`;
  report warnings and actual tests, without changing source-bound verdicts.
