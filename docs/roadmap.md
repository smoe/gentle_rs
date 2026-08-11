# GENtle Roadmap

Last session: 2026-08-03 - made isoform-study normalization directly reusable
and added one content-bound, prevalidated multi-study second-stage approval;
manual GUI smoke for the unreleased `v0.1.0-internal.10` release remains before
the tag

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

Active next-release aim: TP73 genome-anchored evidence viewer.

Release story: GENtle can open the public GRCh38.p14 TP73 locus and let a user
inspect exons, introns, repeats, array/evidence tracks, prepared CUT&RUN
intervals, paired-read ROI support, regulatory/motif-context reasoning,
TFBS annotations, and coordinate-build provenance in the DNA viewer.

Smoke status: deterministic headless proofs and automated release checks are
green-path requirements for `v0.1.0-internal.10`. Before tagging, run the
manual DNA-viewer and Gene Locus Evidence composer smoke so the graphical
release story is inspected as well as regenerated.

Proof path:

- [Headless regeneration workflow](examples/workflows/tp73_genome_evidence_viewer_release_proof.json)
  loads `test_files/tp73.ncbi.gb`, overlays tiny local repeat and Clariom-style
  array fixtures, exercises CUT&RUN V1-V3 with synthetic paired reads, adds
  TFBS fixtures, and emits SVG/JSON/TSV proof artifacts.
- Public [TP73 evidence-viewer runbook](tp73_genome_evidence_viewer_runbook.md).
- Offline [PATZ1 Gene Locus Evidence workflow](examples/workflows/patz1_gene_locus_evidence_offline.json)
  and generated tutorial align transcript models, PSR/JUC effects, grouped
  occupancy, TP73 motif scoring, provenance, and assay continuation on one
  negative-strand axis.
- [Fixture provenance and regeneration notes](../test_files/fixtures/evidence_viewer/README.md).
- Draft [`v0.1.0-internal.10` release notes](release_notes/release_notes_v0.1.0-internal.10.md).

Release acceptance:

- The proof workflow remains offline-safe and writes non-empty sequence,
  splicing-expert, TFBS SVG, repeat-materialization JSON, CUT&RUN regulatory
  JSON, and CUT&RUN coverage/cut-site/fragment TSV artifacts.
- Version metadata, generated documentation, capability/parity checks, and the
  deterministic proof workflow pass the pre-tag validation recorded in the
  versioned release notes.
- Full UCSC `rmsk`, raw CEL, full SRA, and genome downloads remain optional
  external resources; CI uses tiny local fixtures without restricting the
  engine or GUI contracts to those fixture sizes.

Release cut line:

- `v0.1.0-internal.10` may be tagged after the automated pre-tag matrix and the
  manual DNA-viewer/Splicing Expert locus-composer smoke are recorded.

Pre-release finishing scope:

- Run the manual GUI smoke from the TP73 and PATZ1 runbooks and fix only
  evidence-inspection or composition problems supported by that review.
- Discuss the deprecated ClawBio shell-normalizer modes before their earliest
  possible removal in `v0.1.0-internal.11`.

## Next Session Priorities

1. Keep the TP73 evidence-viewer and PATZ1 locus-composer proof workflows green
   and offline-safe through the `v0.1.0-internal.10` tag.
2. Run the manual GUI smoke from both runbooks, including file relocation,
   graphical preview, JSON/SVG/PDF export, and evidence-to-assay continuation.
3. Review the ClawBio shell-normalizer deprecations before deciding whether any
   compatibility modes should be removed in `v0.1.0-internal.11`.
4. Preserve headless/GUI parity for repeat, array, CUT&RUN V1-V3, TFBS, and
   feature-detail views without promoting evidence overlap into conclusions.

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
read-only `planning protein-expression-handoff` route, which now surfaces
product-readiness context for stored DNA/protein sequences, chassis/route
candidates, GeneArt service preflight scaffolding, and the required
yield/folding/chassis review questions without creating constructs.
`product_definition.readiness` now drives `suggested_next_actions[]`; text and
JSON outputs expose the same provider actions. A provider-ready selected
product also produces a product-specific GeneArt review draft referencing its
project `seq_id` or `protein_seq_id`, while undefined products retain only the
clearly labelled bundled example and no provider action. The route now accepts
the versioned `gentle.protein_expression_requirements.v1` record; reviewed
topic records resolve matching questions, and an explicit in-house-only choice
withholds provider actions.

Keep the next synthetic-biology bridge slices in this order:

1. Derive a fully populated, review-only provider-neutral external-service
   request draft from
   the selected product plus those reviewed requirements, then project it to
   GeneArt only after provider selection. Keep submission and ordering out of
   scope.
2. Add a thin Synthetic Biology / Protein Expression Handoff inspector over
   the shared report. It should display readiness, unanswered requirements,
   sequence/CDS/tag evidence, candidates, and shared suggested actions; the
   same window can later expose existing reporter recommendation and reporter
   construct-handoff reports as sibling views.
3. Rank chassis, vector/tag
   routes, and provider candidates from structured biology, provider
   capabilities, and local inventory. Emit score components and rejection
   reasons; do not silently codon-optimize or create constructs.

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
- Clariom D and splice-isoform evidence now have engine-owned interpretation
  layers beyond prepared-track projection. Current explicit APT
  imports can preserve supplied PM probe-level intensity matrices in
  `probe_intensity_chrom_order.csv` and project rows marked
  `probe_level_input` into genome-anchored array features. The first read-only
  `InterpretProbeRegionEvidence` report now compares projected array features
  with transcript/exon geometry and records shared transcript, parent probeset,
  multi-hit-not-assessed, and coordinate-projection ambiguity; the
  sequence-window GUI panel can run and preview that report after projection.
  The report now includes explicit per-evidence transcript mappings with exon
  ordinals, exon ranges, junction spans, overlap base counts, and conservative
  geometry score/basis fields plus review-only transcript labels for unique,
  shared, constraining, or absent geometry. Those audited records now feed the
  read-only `gentle.gene_isoform_evidence.v2` ledger, which composes curated
  transcript families, RNA-read/cDNA/EST support, expression, probe
  constraints, selected projected occupancy tracks, existing qPCR candidates,
  per-source/contrast measurements, annotation-derived protein identities, and
  rule-based assay-triage tiers while retaining unknown/not-evaluated states.
  `gentle.gene_transcript_assay_routine.v1` can now join that exported ledger
  with persisted common-control, junction, and endpoint assay panels by
  digest, without rerunning them. Future work may add
  persisted probe/expression report stores and separately sourced
  antibody-epitope evidence; it must not turn overlap, protein mass, or array
  intensity into an isoform-validation or antibody-compatibility claim.

### Reserved For Anze: Repeat/Similarity Inspection Follow-Up

A later phase stays intentionally separate: broader mapping of single-sequence
operations across sequence sets. Portable inspection actions, structured GUI
evidence panes, cross-adapter parity, curated repeat-family provenance, and
objective-aware severity presentation are recorded in `docs/CHANGELOG.md`.

## Phase D: Visualization And Workflow UX

Use the GUI as the human inspection surface for engine-owned evidence. Continue
dense DNA-map readability, alternative-splicing view polish, gel arrangement
editing, feature editing, contextual interpretation links, visual regression
fixtures, and scroll/zoom hardening when they fit the selected release story. Useful
work here improves inspection clarity, deterministic exports, contextual links
to evidence records, and manual-smoke reliability; defer large visual redesigns
unrelated to the next release aim.
- Evidence-viewer follow-up: keep the Splicing Expert evidence ledger readable
  on large real loci and add direct report-store selectors if probe/expression
  reports become persisted project objects.
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
- Implement the opt-in native-Mistral inner-agent conformance routine described
  in `docs/testing.md`; self-skip unless `MISTRAL_API_KEY` is supplied.
- Supplemental restriction-enzyme usage annotations beyond REBASE.
- Floating restriction-site detail popover/window if the Description panel is
  too easy to miss.
- SnapGene-style plasmid-map presentation parity and dense selected-site polish.
- Engine-owned portable findings/artifact inspection for agent-driven work.
- Engine-owned exhaustive operation-effect metadata so new filesystem-writing
  operations cannot bypass rollback-safety classification or provenance path
  collection when the `Operation` enum grows.
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
- CUT&RUN controls/replicate comparison, calibrated differential support,
  formal peak calling, and whole-genome read processing beyond the ROI-first
  release scope.
- Primer3 parity, virtual PCR/off-target filtering, multiplex tiling, LAMP, and
  allele-specific assay families.
- First-class global RT-PCR specificity/redesign workflow building on the
  landed local assembly-aware primer variant screen and combining prepared
  genome and transcriptome searches, intended-product/paralogue policy,
  junction evidence, Primer3 plus exact-run QC, ranked redesign candidates,
  and one provenance-bound review bundle.
- GuideRNA off-target ranking and macro-template packaging.
- Cross-tool parity synthesis for Serial Cloner, MacVector, and SnapGene.
- Post-release hardening backlog that is not tied to the next release scope.
- Weekly/monthly maintenance chore automation rollout from
  [`maintenance_chore_plan.md`](maintenance_chore_plan.md).
- Browser/WebAssembly frontend portability after core/headless contracts settle.
- Expand the landed policy-driven collection surface beyond promoter
  derivation, primer specificity, restriction-site and TFBS hit scanning,
  fingerprint-locked restriction digestion, and atomic physical-container
  pool export plus explicit gene-set-to-pool aliquoting and registry-declared
  serial-arrangement gel rendering to later shared BLAST, generic arrangement,
  and storage-projection adapters described in
  [`gui_gene_set_collection_operations_plan.md`](gui_gene_set_collection_operations_plan.md).
- Gene-set enrichment analysis over resolved logical sets, with an explicit
  background universe, identifier namespace/mapping audit, ontology and cache
  provenance, biological-context compatibility, multiple-testing correction,
  unresolved-member accounting, and no causal-regulation claim inferred from
  enrichment alone.
- Add paralog mapping/resolution only with a concrete engine consumer,
  preserving relationship provenance and one-to-many ambiguity without
  inferring functional equivalence from matching symbols.
- Full glossary/help-generation inversion so `docs/glossary.json` becomes a
  generated or validated projection of engine/protocol-side descriptors.
