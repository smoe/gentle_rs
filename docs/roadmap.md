# GENtle Roadmap

Last session: 2026-05-16 - gentle-engine Phase 2 extraction began with
`iupac_code`; TP73 evidence-viewer release aim remains active

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

Active next-release aim: TP73 genome-anchored evidence viewer.

Release story: GENtle can open the public GRCh38.p14 TP73 locus and let a user
inspect exons, introns, repeats, array/evidence tracks, CUT&RUN-style BED
intervals, TFBS/regulatory annotations, and coordinate-build provenance in the
DNA viewer.

Smoke status: manual GUI smoke for the TP73 evidence-viewer runbook is still
pending in this worktree. Keep the release-facing notes short until that smoke
path has passed locally.

Proof path:

- Headless regeneration:
  `docs/examples/workflows/tp73_genome_evidence_viewer_release_proof.json`
  loads `test_files/tp73.ncbi.gb`, overlays tiny local repeat, Clariom-style
  array, CUT&RUN-style BED, and TFBS fixtures, and emits SVG/report artifacts.
- Public runbook: `docs/tp73_genome_evidence_viewer_runbook.md`.
- Fixture bundle: `test_files/fixtures/evidence_viewer/` with provenance and
  deterministic regeneration notes.

Release acceptance:

- The proof workflow remains offline-safe and writes non-empty sequence,
  splicing-expert, TFBS SVG, and repeat-materialization JSON artifacts.
- The GUI smoke path opens the proof state, shows anchor/build status, toggles
  repeat and array layers, and exposes informative feature details for one
  exon, one transcript/intron context, one repeat, one array row, and one
  CUT&RUN-style interval.
- Full UCSC `rmsk`, raw CEL, full SRA, and genome downloads remain optional
  external resources; CI uses only tiny local fixtures.

Release cut line:

Post-release scope:

## Next Session Priorities

1. Keep the TP73 evidence-viewer proof workflow green and offline-safe while
   finishing release-facing documentation.
2. Run the manual GUI smoke from the runbook and fix only evidence-viewer
   inspection problems that block the release story.
3. Preserve headless/GUI parity for repeat, array, CUT&RUN-style BED, TFBS, and
   feature-detail views; avoid adding broader biology conclusions.
4. Keep private grant material out of the GENtle source tree; only general,
   reusable GENtle improvements should return here.
5. After the smoke path passes locally, prepare concise release notes around
   the TP73 viewer story as the next concrete release-facing artifact.

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
than adding more biology-specific wrapper modes. After the first GENtle-side
ClawBio panel, the next bridge work should stay contract-led: reserve any
action-level `skill_alias` field on the ClawBio side before GENtle honors it,
and avoid adding routing or planner semantics to GENtle.

External-service bridge status: the provider-neutral schema/shell layer is now
concrete through `services providers list`, `services providers doctor`,
`services delivery-route`, `services project-preflight`, and
`services project-quote`. Built-in provider configuration covers GeneArt quote
handoff for DNA fragments/cloned genes/reorders/mutagenesis/protein expression
and Metabion handoff routes for single-tube DNA oligos and m-block fragments.
Generic "deliver this sequence" wording should go through delivery-route
classification first so sequence kind, length, and construct/protein context
choose the provider/service route before any quote packet is prepared. Agent
work should stay focused on ClawBio intent/example parity and clear "handoff,
not submission" wording; direct GeneArt API use, WOP/OCI automation,
credential handling, and commercial order state remain later explicitly
confirmed integration phases.

Near-term consult rule: keep `planning consult cloning` narrow and
deterministic until richer construct reasoning lands. The preferred v1 surface
is one best candidate per the 11 catalogued routine families plus structured
helper/vector ranking; `seq_id` stays traceability-only and unresolved
marker/promoter/MCS constraints stay explicit questions rather than hidden
heuristics.
High-yield protein-production requests should go through the separate
read-only `planning protein-expression-handoff` route, which surfaces product
context, chassis/route candidates, GeneArt service preflight scaffolding, and
the required yield/folding/chassis review questions without creating
constructs.

MCP follow-up: the typed `tools/list` surface is curated rather than a
one-command/one-tool mirror. Keep `docs/agent_interface.md` as the exclusion
ledger for shell commands without dedicated MCP tools, and promote a command to
a route-specific MCP tool only when its JSON input/output schema and
mutating/external safety semantics are stable.

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
broad crate surgery that is not tied to the selected release story.
- Phase 2 has begun with `crates/gentle-engine` re-introduced and
  `iupac_code` moved as the first root-independent execution-side module behind
  the existing root compatibility shim.
- New operations adopt the three-level parity model: engine and reachability
  parity from the parity matrix in [`gui_cli_mcp_parity.md`](gui_cli_mcp_parity.md);
  first-class surfacing decisions land per adapter and may be staged.
- The parity matrix tests reachability, not surfacing; surfacing decisions are
  recorded per row with a one-line justification.
- Inline/stateless operand follow-up order from
  [`inline_operand_audit.md`](inline_operand_audit.md): `RenderSequenceSvg`,
  `RenderRnaStructureSvg`, `RenderTfbsScoreTrackCorrelationSvg`,
  `ComputeDotplot`, then `ComputeFlexibilityTrack`.
- Clariom D and splice-isoform evidence need an engine-owned interpretation
  layer beyond the current prepared-track projection. Current explicit APT
  imports can preserve supplied PM probe-level intensity matrices in
  `probe_intensity_chrom_order.csv` and project rows marked
  `probe_level_input` into genome-anchored array features; future work should
  use those mapped probe oligos/probesets with explicit build, liftover,
  multi-hit, and ambiguity provenance to identify compatible splice isoforms
  and rule out incompatible ones without hiding uncertain evidence.

## Phase D: Visualization And Workflow UX

Use the GUI as the human inspection surface for engine-owned evidence. Continue
dense DNA-map readability, alternative-splicing view polish, gel arrangement
editing, feature editing, contextual interpretation links, visual regression
fixtures, and scroll/zoom hardening when they fit the selected release story. Useful
work here improves inspection clarity, deterministic exports, contextual links
to evidence records, and manual-smoke reliability; defer large visual redesigns
unrelated to the next release aim.
- Evidence-viewer follow-up: make Clariom D probe/probeset evidence and
  splice-isoform constraints inspectable together, so a user can see which
  oligos, exons, junctions, and mapped intervals support, constrain, or exclude
  each transcript isoform.
- External-services GUI follow-up: keep the window a thin inspector over the
  shared provider/preflight/quote contracts, then improve product-specific
  starter templates, validation previews, and exported-bundle review affordances
  before considering provider-specific portal/API actions.

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
- Luebeck S1 vector-resource follow-up: resolve sequence/procurement URLs and
  redistribution clearance before bundling actual sequences for the metadata
  candidates.
- VKORC1 / warfarin Factor X follow-up tutorial idea: after the maintained
  VKORC1 rs9923231 promoter-reporter tutorial, design a separate companion
  story for downstream warfarin-relevant Factor X target-site constructs once
  the construct family, controls, readout, and bench-facing success criteria
  are stable enough for a validated walkthrough.
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
