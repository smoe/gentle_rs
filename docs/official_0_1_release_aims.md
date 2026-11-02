# First Official 0.1 Release Aims

Status: owner-earmarked direction, 2026-11-02. Detailed implementation slices,
acceptance thresholds, release date and candidate are not selected here.
This records Steffen's proposed use of Glen's prioritized requests for the
first official `0.1` release; it is not a claim that these capabilities ship or
are verified. The current version remains `0.1.0-internal.13`.

Glen's source review was bound to `7f73a42a`, not later local main. Preserve his
[positive discovery/decision verdict](GENtle_discovery_verification_7f73a42a.md)
separately from the failing glossary gate and newer egui changes. The typed SAM
fixture repair is already integrated in `1e46525c`; do not merge the audit branch
again. Follow the [current repair handoff](glen_next_verification_handoff.md#egui-037-and-ported-glossary-repair-2026-11-02)
for pending compilation, focused tests and native verification.

## Immediate Prerequisites

1. Restore green main by verifying the existing glossary/egui repairs and fixing
   any independently reproduced remaining failures. Do not weaken typed parsing
   or treat a scoped discovery pass as an overall CI pass.
2. Prepare a disciplined `.13` candidate under its existing
   [scope and acceptance ledger](release_notes/release_notes_v0.1.0-internal.13.md).
   At actual release preparation, require matching Windows/macOS/full-Linux CI,
   native packages, container, extracted-package smokes and the documented
   GUI/tutorial checks. Candidate selection, dispatch, tags and publication
   remain owner-authorized actions. The aims below do not expand `.13`'s scope.

## Earmarked 0.1 Outcomes

All six outcomes remain open for implementation/verification planning. Existing
code and historical evidence are starting points, not blanket acceptance.

### 1. Installable OpenClaw Skill

Extend the existing [ClawBio integration](../integrations/clawbio/README.md),
not a second command language. Provide a documented install/update route,
pinned GENtle binary or container identity, a readiness check, and portable
discovery for OpenClaw, Codex and Claude. Keep wrappers thin; do not add more
biology-specific modes. Results distinguish read-only inspection, proposals,
confirmed mutations and refusal, with artifact checksums and receipts.

Evidence needed: clean installation and upgrade on explicitly supported hosts,
matching descriptor/runtime identities, discovery/execution through shared
contracts, approval refusals and reopenable artifact-first results. Host versions
and the supported installation mechanism still need a bounded implementation plan.

### 2. Evidence-Backed Golden Path

Populate [capability maturity](capability_maturity.md) for approximately 8-12
selected capabilities, rather than promoting every route to a green badge.
The initial selection should cover runtime/readiness, sequence-context
inspection, gene/isoform and Splicing Expert rendering, PCR/qPCR design,
complete specificity inspection, cloning-method consultation, sequencing
confirmation, one provider-neutral handoff and one RNA-read evidence workflow.

Evidence needed per capability: actual execution, source/binary/platform/profile
and input/resource hashes, limitations, and native GUI or biological review where
the claim requires it. Keep availability, adapter coverage, model behavior and
scientific maturity separate. Unknown or untested scopes remain unassessed.

### 3. Flagship Global RT-PCR Specificity And Redesign

Build on the [primer specialist plan](primer_design_specialist_plan.md) and
[reference-bound specificity tutorial](tutorial/04-09_reference_bound_primer_specificity.md).
Combine prepared genome and transcriptome searches, explicit intended-product,
isoform and paralogue policies, complete/incomplete search provenance, Primer3,
exact-run QC and ranked redesign candidates in one review bundle. Link primers,
products, off-targets, exclusions and hashes; do not hard-wire the current panel.

Evidence needed: bounded real-reference import performance, independently checked
synthetic cases, authentic-reference execution, ranking/exclusion explanations
and reproducible review bundles. An incomplete search must not become a clean
specificity claim. This supports an informed primer-order review; it neither
authorizes ordering nor proves experimental success.

### 4. Long-Read Cohorts And Evidence-To-Curation

Extend the existing [RNA-panel contracts](rna_panel_gui_processing_plan.md) and
[tutorial 10.04](tutorial/10-04_targeted_nanopore_fastq_panel.md) with sample-sheet
condition, time point, replicate and extraction/library notes; calibrated
background/control handling; scalable multi-sample/gene presentation; explicit
missing-target and origin-confidence accounting. Produce content-bound curation
proposals with separately confirmed application.

Evidence needed: synthetic cohort/control oracles, explicit calibration limits,
missing-versus-zero and ambiguity controls, streaming/cancellation behavior,
tamper/stale refusal, and usable cohort views. Real private inputs require their
own authorization. Observed structures, annotation compatibility and biological
consequence remain separate; incomplete sequence/model/phase evidence yields
`not_assessable`, not invented isoforms, abundance or consequence.

### 5. Agent-Safe Operation Lifecycle

Continue the [engine-owned lifecycle work](operation_lifecycle_policy_proposal.md)
with exhaustive operation-effect metadata, portable findings/artifact inspection,
preview digests, confirmation, atomic project apply and explicit rollback limits.
Give journals, lineage, checkpoints and files clear effect owners. Preserve the
implemented safeguards in [Slices 3](operation_lifecycle_slice3_plan.md),
[4](operation_lifecycle_slice4_plan.md), [5](operation_lifecycle_slice5_plan.md)
and [6](operation_lifecycle_slice6_plan.md), rather than replacing them in wrappers.

Evidence needed: effect coverage, confirmation-before-I/O, cancellation/staleness
and recovery tests across GUI, CLI, MCP, JS/Lua and ClawBio. Distinguish project
rollback from separately published files and retain orphan/failure receipts.
Safety applies throughout the other aims; this is not permission to defer it
until the end or enable autonomous editing now.

### 6. Targeted GUI Performance And Inspection

Measure feature density as well as sequence length using the existing
[latency plan](dna_feature_rendering_latency_plan.md). Prioritize responsive long
commands/cancellation, large-locus Splicing Expert readability, cohort
virtualization, deterministic exports and save/reopen. Avoid a broad redesign.

Evidence needed: representative retained inputs, same-host/profile measurements,
native readability and responsiveness checks, and reproducible exports/reopening.
Numeric budgets and representative datasets must be chosen before performance
acceptance; synthetic frame tests are not native timing evidence.

## Delivery And Acceptance Boundary

Follow Glen's sequence: main stabilization, `.13` release preparation, installable
skill plus golden-path evidence, flagship specificity/redesign, long-read cohorts,
then deeper lifecycle/GUI refinement. Safety remains a prerequisite at every step.
Select narrow implementation slices; do not start all workstreams at once.

At official-release preparation, select one candidate and consolidate affected
cross-platform, installer/container, extracted-package, tutorial, native and
readiness evidence. Every verdict retains its actual source/binary/input scope;
an earlier pass does not silently accept a changed candidate. Detailed acceptance
criteria and any justified scope reduction require an explicit owner decision.
There is no selected release SHA, schedule, version rollover or publication here.

ClawBio/Glen owns conversational strategy, literature/predictor synthesis and
cross-tool orchestration. GENtle owns deterministic sequence reasoning, project
effects, evidence artifacts, receipts, refusal and enforcement of approval gates.
External orchestration cannot grant itself GENtle mutation authority. Direct
vendor ordering, clinical interpretation and autonomous biological conclusions
remain out of scope; no software evidence claims human scientific approval.

The existing [verification cadence](decisions.md#development-and-release-verification-cadence)
still applies: no interim development freezes, recurring SHA-sync requests or
bundle treadmill. This earmarking authorizes no push, CI dispatch, provider/model
call, private-data analysis, screenshot capture, package publication or release.
