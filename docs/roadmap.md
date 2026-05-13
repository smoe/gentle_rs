# GENtle Roadmap

Last session: 2026-05-12 - interim release shipped

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

## Current Status

An interim release is out. The previous pre-tag checklist has been retired from
this roadmap.

The next release scope is intentionally open until it is discussed and narrowed
to one primary story, one proof path, and a compact test/documentation slice.
Once selected, that scope should be recorded here as the next active release
plan instead of reviving the old TP73 pancreatic benchmark checklist.

## Next Session Priorities

1. Choose the next release aim: one primary user-facing story, one deterministic
   proof path, and a small set of test/documentation acceptance criteria.
2. Convert the chosen aim into a compact release plan with explicit artifacts,
   commands, GUI smoke expectations, and adapter-parity checks.
3. Select at most one or two implementation threads from the active phases below
   so the next release remains coherent rather than feature-scattered.
4. Keep private grant material out of the GENtle source tree; only general,
   reusable GENtle improvements should return here.
5. Continue small correctness, reproducibility, and parity fixes when they
   support the chosen release aim or prevent obvious regressions.

Current non-goals:

- Do not add private proposal or unpublished grant content to this repository.
- Do not start broad engine extraction, GUI redesign, or infrastructure work
  before it is tied to the selected next-release story.
- Do not add GUI-only business logic for convenience; keep shared
  engine/shell/protocol contracts as the source of truth.
- Do not promote speculative assistant ideas into protocols without
  confirmation.

Useful session close:

- Roadmap names the selected next-release aim or explicitly says it remains
  undecided.
- Completed outcomes from the session are in [`CHANGELOG.md`](CHANGELOG.md).
- Any new durable rule is in [`decisions.md`](decisions.md).
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

For ClawBio/OpenClaw specifically, keep new integration work centered on
descriptor/runtime parity (`INTENTS.json`, `mode=intents`, examples, and
trigger-keyword drift checks) plus explicit scope/presentation contracts rather
than adding more biology-specific wrapper modes.

## Phase B: Cloning Routine Standardization

Continue routine catalog and macro-box work when it fits the selected release
scope. Priorities are richer routine-family preflight, replay-friendly macro instances,
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
broad crate surgery that is not tied to the selected release story. Phase 2
has not started; the `gentle-engine` crate scaffold has been removed and will
be re-introduced when the first execution-side module is ready to move.

## Phase D: Visualization And Workflow UX

Use the GUI as the human inspection surface for engine-owned evidence. Continue
dense DNA-map readability, alternative-splicing view polish, gel arrangement
editing, feature editing, contextual interpretation links, visual regression
fixtures, and scroll/zoom hardening when they fit the selected release story. Useful
work here improves inspection clarity, deterministic exports, contextual links
to evidence records, and manual-smoke reliability; defer large visual redesigns
unrelated to the next release aim.

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
- GUI batch display for RNA-read gene-cohort alignments with virtualization for
  thousands of per-read aligned-column blocks.
- Vendor-protocol and deeper wet-lab process modeling primitives.
- Automatic cross-gene homology anchors on top of explicit `query_anchor_bp`.
- Catalog-extensible gene-group and ontology bridge after the next release scope
  is selected.
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
- Post-release hardening backlog that is not tied to the next release scope.
- Weekly/monthly maintenance chore automation rollout from
  [`maintenance_chore_plan.md`](maintenance_chore_plan.md).
- Browser/WebAssembly frontend portability after core/headless contracts settle.
