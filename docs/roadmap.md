# GENtle Roadmap

Last session: 2026-05-10 - sequence loader engine-GUI seam

Purpose: fast session orientation. This file answers "what next?" and should be
readable in under two minutes. Completed work belongs in
[`CHANGELOG.md`](CHANGELOG.md). Durable implementation constraints belong in
[`decisions.md`](decisions.md). Protocol and schema contracts remain in
[`protocol.md`](protocol.md).

Maintenance invariant:

- Do not add a `Done` entry to this roadmap. Move completed work directly to
  [`CHANGELOG.md`](CHANGELOG.md) in the same session.
- The roadmap may temporarily contain at most one active session's in-flight
  status notes; migrate them before handoff.
- If a note is a durable rule, move it to [`decisions.md`](decisions.md).
- If a note is a speculative idea, keep it as a one-line parking-lot entry.

## Release Gate

Current rule: pre-release work is limited to release confidence,
reproducibility, correctness regressions, and one persuasive end-to-end proof
path. New biology features are frozen by default until the release proof is
settled.

The release proof path is TP73 pancreatic cancer Nanopore cDNA benchmarking via
[`tp73_pancreas_benchmark_runbook.md`](tp73_pancreas_benchmark_runbook.md).
Preserve the generated preflight summary, final report summary, TSV exports,
SVG target-quality export, logs, resolved gene-group catalog snapshots, and
evidence-bundle note.

If cohort-level abundance variation is needed, use the fixed-parameter route in
[`tp73_pancreas_cohort_batch_runbook.md`](tp73_pancreas_cohort_batch_runbook.md).
Use the preflight-derived seed filter consistently, keep per-run state copies
isolated, and merge sample sheets for comparison.

Keep CUT&RUN work to the release smoke/proof slice in
[`cutrun_release_smoke.md`](cutrun_release_smoke.md): prepared-dataset
projection, ROI read interpretation, regulatory-support inspection, existing
TFBS scan/score-track surfaces, and the shared-report GUI inspector. Defer
de-novo motif discovery, differential-expression integration, and improved
scoring models.

Before tagging, run the release-signoff slice: targeted RNA-read/preflight
tests for TP73 positives and TP53/TP63 controls, shell/CLI command-glossary
parity tests, workflow/example tests touched by the proof path,
`cargo check -q`, and one broader `cargo test -q` when local time permits.

Before tagging, run a small manual GUI smoke pass: launch the app, open a known
TP73 sequence/project, open or focus RNA-read Mapping, inspect a saved RNA-read
report, export at least one evidence artifact, and confirm scrolling/redraw is
demonstration-safe.

Keep release notes and install/run documentation aligned with the actual
shipped artifacts and the headless/agent-first positioning.

Continue only when:

- The TP73 proof artifacts can be regenerated or traced from documented commands.
- The GUI smoke path can be demonstrated without relying on hidden local state.
- CLI, shell, MCP, and agent wording describe the same operation contracts.
- Release-note claims name artifacts that exist in the repository or run output.
- Any failing targeted test has a clear blocker note and owner-facing next step.

Stop and fix before tag when:

- A release proof command becomes nondeterministic across repeated runs.
- A GUI action mutates state without the same shell/engine confirmation boundary.
- A release note describes a feature that is not backed by a deterministic path.
- A proof artifact depends on an untracked local-only file.
- A regression affects sequence correctness, report provenance, or safety guards.

Artifacts to name in the release handoff:

- TP73 pancreatic benchmark state path and report summary.
- Preflight summary and seed-filter parameters.
- Final report TSV/SVG exports and evidence-bundle note.
- Any cohort batch sample sheet or merged comparison sheet used for claims.
- Manual GUI smoke notes with observed app version and input path.
- Release-signoff command list with failures, skips, and local constraints.
- Updated install/run/release-note files that describe the demonstrated path.

## Next Session Priorities

1. Complete or verify the TP73 pancreatic benchmark proof bundle and record the
   exact artifact paths needed for release review.
2. Run the release-signoff test slice and note any failing or skipped coverage
   as release blockers.
3. Perform the manual GUI smoke pass against the release proof path and record
   only actionable regressions.
4. Align release notes, install/run docs, and ClawBio/headless handoff wording
   with the artifacts that actually exist.
5. Fix only release-risk regressions or reproducibility gaps before the tag;
   defer broader feature work to post-release phases.

Current non-goals:

- Do not add new wet-lab biology modeling before the release gate is settled.
- Do not widen ClawBio beyond catalog exposure and handoff-oriented actions.
- Do not refactor engine modules solely to make the split cleaner.
- Do not add GUI-only business logic for proof-path convenience.
- Do not promote speculative assistant ideas into protocols without confirmation.

Useful session close:

- Roadmap still names only next work and open phase direction.
- Completed outcomes from the session are in [`CHANGELOG.md`](CHANGELOG.md).
- Any new durable rule is in [`decisions.md`](decisions.md).
- Release blockers are concrete, reproducible, and tied to file paths or commands.
- Unrelated local paper/figure work remains untouched unless explicitly requested.

## Phase A: AI Communication And Safety Plane

Keep MCP, agent, shell, ClawBio/OpenClaw, and UI-intent surfaces thin and
adapter-equivalent. The next work is breadth and safety, not new biology:
expand deterministic tool coverage over shared shell/engine routes, harden
mutating-intent confirmation consistently across agent, voice, and MCP paths,
and preserve the shared `ui ...` intent catalog across menu, command palette,
shell, agent, MCP, and ClawBio discovery. Useful work here is parity,
discoverability, confirmation wording, and `ui_intent`-shaped operator handoff;
defer autonomous planning policy or biology-specific reasoning that is not
backed by current engine records.

## Phase B: Cloning Routine Standardization

Continue routine catalog and macro-box work after the release gate. Priorities
are richer routine-family preflight, replay-friendly macro instances,
dense-lineage controls, and protocol packs that make Golden Gate, Gateway,
TOPO, TA/GC, In-Fusion, NEBuilder HiFi, and deeper restriction variants feel
first-class without duplicating biology in adapters. Useful work here extends
engine-owned descriptors, portable macro records, candidate auditability, and
protocol-pack wording; defer new cloning families unless they reuse the current
routine framework.

## Phase C: Engine And Protocol Parity Hardening

Promote adapter-level helpers into engine-owned operations where they matter,
keep process/run-bundle artifacts portable, and converge persisted
computational conclusions on one provenance vocabulary. Sequence file loading
now routes through the shared DNA sequence module rather than the GUI top-level
type; continue XML/SnapGene, sequencing-confirmation, primer/qPCR, projection,
construct-reasoning, RNA-read, and ClawBio/MCP parity work only through shared
contracts. Useful work here is contract hardening, adapter-helper promotion,
stack-safe helper splits, protocol snapshots, and targeted parity tests; defer
broad crate surgery that is not needed for the current release proof. Phase 2
has not started; the `gentle-engine` crate scaffold has been removed and will
be re-introduced when the first execution-side module is ready to move.

## Phase D: Visualization And Workflow UX

Use the GUI as the human inspection surface for engine-owned evidence. Continue
dense DNA-map readability, alternative-splicing view polish, gel arrangement
editing, feature editing, contextual interpretation links, visual regression
fixtures, and scroll/zoom hardening after the release proof is stable. Useful
work here improves inspection clarity, deterministic exports, contextual links
to evidence records, and manual-smoke reliability; defer large visual redesigns
unrelated to release confidence.

## Phase E: Integration Polish And Deferred Policy Items

Keep integration work focused on reproducible distribution and policy clarity:
OCI/Apptainer validation, release attachment decisions, cross-application
clipboard contracts, documentation automation, GUI screenshot atlas work after
explicit screenshot approval, and broader adapter/documentation polish. Useful
work here validates headless MCP deployment, aligns install docs with actual
outputs, and decides release attachments versus git-tracked assets; defer
infrastructure expansion that does not reduce release risk.

## Phase F: Interpretation Later

Later interpretation work includes image/sketch-to-state patch proposals,
richer assistant-facing findings records, deeper wet-lab process modeling,
repeat/mobile-element curation, and cohort-scale comparative reasoning. These
must remain explicit, inspectable proposals with confirmation before mutation.
Useful work here improves proposal formats, evidence records, reproducible
cohort comparisons, and rollback paths for state patches; defer autonomous
wet-lab conclusions or unconfirmed mutations.

## Parking Lot

- Optional OS credential-store persistence for Agent Assistant API keys.
- Optional tiny generation probe for quota verification.
- Supplemental restriction-enzyme usage annotations beyond REBASE.
- Floating restriction-site detail popover/window if the Description panel is
  too easy to miss.
- SnapGene-style plasmid-map presentation parity and dense selected-site polish.
- Engine-owned portable findings/artifact inspection for agent-driven work.
- Vendor-protocol and deeper wet-lab process modeling primitives.
- Automatic cross-gene homology anchors on top of explicit `query_anchor_bp`.
- Catalog-extensible gene-group and ontology bridge after release proof work.
- External-service provider/CRO integration once deterministic local contracts
  are stable enough to wrap.
- Primer-walking support and iterative read/contig data management.
- Helper-construct terminology migration away from legacy "helper genome"
  wording.
- Broader CUT&RUN development beyond the release smoke slice.
- Primer3 parity, virtual PCR/off-target filtering, multiplex tiling, LAMP, and
  allele-specific assay families.
- GuideRNA off-target ranking and macro-template packaging.
- Cross-tool parity synthesis for Serial Cloner, MacVector, and SnapGene.
- Post-delivery hardening backlog that is not part of the release gate.
- Weekly/monthly maintenance chore automation rollout from
  [`maintenance_chore_plan.md`](maintenance_chore_plan.md).
- Browser/WebAssembly frontend portability after core/headless contracts settle.
