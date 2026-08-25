# GENtle Changelog

- Grounded every Agent Assistant request in a bounded, path-free inventory of
  manifest-backed local references, and supplied Pi Local with the same shared
  parser contract as native providers. Agents are now instructed to prefer a
  prepared-genome `extract-gene` plus strand-aware `extend-anchor` chain before
  suggesting network retrieval. Direct Ensembl gene fetch remains available as
  a confirmation-gated fallback and now accepts an assembly assertion plus
  symmetric or directional 5-prime/3-prime flanks; imported annotations retain
  the actual expanded genomic span on both strands.

- Fixed the DNA-viewer opening status so unrelated pending hosted-window probes
  cannot keep `Opening DNA Sequence Viewer…` visible after the sequence viewer
  has completed its first render.

- Kept committed tutorial assay artifacts deterministic across Git checkouts by
  normalizing their selection-audit generator revision to the package version;
  ordinary engine reports still retain the complete source revision.

- Routed the separate Sequence tools surface through GENtle's hosted-window
  compatibility wrapper instead of constructing a direct `egui::Window`.

- Fixed Array setup navigation after sequence-specific controls moved into the
  separate Sequence tools window: the action now opens and scrolls directly to
  the Clariom D panel. Added a dedicated read-only `Microarrays` Configuration
  tab with file-role/source guidance and bounded checkout-local status for
  Thermo Fisher support archives, E-MTAB-14704 inputs, analyzed TP73 tables, and
  the synthetic demonstration output. The sequence panel reuses those checks
  for explicit local suggestions without downloading files, overwriting user
  paths, or presenting synthetic fixtures as biological evidence.

- Prevented automatic first-frame lineage panel fitting from marking a freshly
  opened, saved project as modified. Sequence-node opening now shows the
  translated `Opening DNA Sequence Viewer…` status while the viewer is queued
  or completing its first render, and removes it once the viewer is open.

- Moved sequence-specific configuration and execution controls out of the DNA
  map toolbar into a separate resizable Sequence tools window. Engine Ops,
  Shell, genome-anchor extension, repeat-index configuration, and the Clariom D
  probe-region workflow no longer consume map height; the new workspace has
  two-dimensional scrolling and German/English labels for its main controls
  and Clariom workflow stages.

- Added explicit, bounded local-document grounding for the Agent Assistant.
  Exact supported text paths in GUI/CLI prompts (or the GUI `Attach text
  document...` picker) now produce validated `gentle.agent_local_documents.v1`
  request context with UTF-8 checks, byte/truncation limits, SHA-256 provenance,
  typed read warnings, and one prioritized local Markdown guide-link level.
  Pi remains tool-free: it can use supplied roadmap/runbook text to propose
  reviewed GUI intents or coach manual observations without browsing the
  project tree. The UI warns that copied contents are sent to the selected
  provider.

- Added live Pi command-shape probing to Agent Assistant setup. `agents preflight
  pi_local_stdio --live` and GUI `Test Setup` now invoke the selected Pi binary
  with GENtle's exact ephemeral no-tools/no-session flags plus `--help`, catching
  unsupported local Pi CLI flags without sending a prompt to a model.

- Added explicit, window-bound screenshot help for Agent Assistant. `Agent help`
  captures the GENtle viewport where it was clicked before opening the
  assistant, presents a local preview, and sends that image only after the
  user submits the prompt to an image-capable provider. Optional macOS native
  full-window capture uses ScreenCaptureKit, is restricted to GENtle-owned
  windows, and fails clearly when Screen Recording permission is unavailable.
  Codex/Pi bridges and native multimodal providers share the attachment
  contract; stored conversation turns retain path-free provenance only.

- Fixed invisible Agent Assistant results for GENtle-local read-only commands.
  `/help` now opens the built-in Help window at the Shell Commands tab, and
  `/list` renders a compact current-project overview instead of recording only
  `success - executed`. The latest local command result remains available in a
  dedicated panel with complete JSON copy support. `/list` object rows now expose
  contextual shared-shell actions and can focus the complete main-window project
  overview rather than duplicating its full lineage/analysis model; shared shell
  payloads and headless behavior are unchanged.

- Added Pi as a local, multi-provider Agent Assistant harness. The co-shipped
  bridge runs an ephemeral tool-free Pi request in an isolated directory,
  keeps conversation state in GENtle, discovers provider-qualified models via
  `pi --list-models`, forwards explicit model choices, and never reads Pi's
  credential store. The generic stdio response boundary now also accepts one
  complete top-level Markdown JSON fence while retaining strict schema,
  unknown-field, and surrounding-prose rejection.

- Clarified the DNA viewer's microarray state and made nearest-array navigation
  reusable. An empty layer now opens `Array setup...` directly at the built-in
  Clariom D Human profile, populated layers distinguish in-view from
  full-sequence counts, and the context menu can focus the nearest projected
  interval. The additive shared feature-query contract and shell option
  `--nearest-to` expose the same deterministic distance ordering to CLI and
  embedding clients; transcript association remains an explicit probe-region
  interpretation step.

- Hardened the approval-bound transcript-assay fallback adapters after a real
  IRF9 replay. The shared shell worker now reserves enough stack for the nested
  strict/fallback reports, selection audit, and combined gel payload. Ordinary
  Cargo builds also bind the repository HEAD into one package-plus-commit
  source revision used by both fallback execution and selection-audit
  provenance; source archives without Git remain explicitly package-bound.

- Added an approval-bound informative fallback for infeasible strict
  transcript-primer panels. The engine preserves the typed strict coverage
  failure, derives a separate best-effort operation only for that failure
  class, and binds the parent operation, policy, search context, exact diff,
  engine revision, and outputs in a persisted execution audit. The new
  budgeted selector counts both presence/absence and threshold-resolvable
  product-size distinctions, reports a bounded-greedy rather than optimality
  claim, and exposes marginal/redundancy/unresolved-target evidence. Every
  fallback stays visibly partial and includes a combined virtual gel or typed
  no-products reason. Planner/workflow digests, shared shell/MCP/direct JSON,
  PCR Designer controls, publication projection, and synthetic IRF9 acceptance
  evidence use the same engine contract; legacy requests remain fail-closed.

- Completed replayable transcript-assay selection provenance across
  interpretation, planning, saved reports, and publication. Probe-region
  interpretations can bind an effect threshold to its policy source/digest
  and normalized request digest; JUC evidence receives typed preferred,
  required, below-threshold, or incomplete dispositions. Selection audits now
  carry computed/legacy status, method, generator revision, stable target
  identity, and separate observation/projection/matched-target counts. Legacy
  artifacts never borrow later thresholds or present absent audits as
  computed zeroes. Gene-set publication can digest-bind a primer row to one
  transcript-panel assay and render only its sanitized engine decision.

- Clarified transcript-assay selection reports for differential Clariom JUC
  validation. Pair summaries now distinguish geometry-only junction evidence
  from rows that pass an explicit contrast/effect threshold, retain the
  qualifying contrast and threshold, and expose incremental, leave-one-out,
  and redundant binary transcript-detection value. Required response-
  associated junction assays with no exclusive discrimination are now stated
  as separately retained validation obligations rather than appearing to add
  isoform discrimination; descriptive effects remain explicitly distinct
  from significance or isoform-abundance claims.

- Made exhaustive cDNA primer-specificity handoff import scale with complete
  BLAST HSP count. New handoffs retain validated `qlen`/`qseq`/`sseq` fields
  and convert complete-query HSPs directly into persisted full-primer
  evidence; legacy and partial HSPs keep the established subject-window
  semiglobal fallback. Existing subject-grouped bounded product pairing and
  all native completeness, intended-target, allowlist, and provenance gates
  remain authoritative. Independent linear-fallback reports are not promoted
  into native GENtle passes.

- Extended terminal-exon sequence-specific RT-primer pools with optional,
  exhaustive local genomic-DNA BLAST evidence. Reports now retain exact
  full-length hits, raw HSP counts, genomic-uniqueness readiness, and prepared
  database/tool/content-fingerprint provenance; the PCR ClawBio skill exposes
  the same reviewed request as an explicit delegated intent.

- Improved terminal-exon sequence-specific RT-primer pools from greedy caller-
  order selection to deterministic bounded pool-wide selection. Reports retain
  the exact policy and evaluated-state count, and the existing all-against-all
  interaction matrix now describes the jointly selected complete oligos and
  their variable segments.

- Made biological workflow extension an explicit GENtle development path.
  A new guide classifies related requests as parameterization, engine-service
  composition, or a missing primitive; documents reusable private
  transcript/window/coordinate/oligo seams; and provides parity, provenance,
  non-claim, and edge-test checklists. Agent Assistant now offers an `Adapt a
  biological workflow` prompt and returns an implementation brief instead of
  inventing commands for unsupported biology. The terminal-exon RT-primer
  handler was decomposed onto named private transcript-resolution,
  canonical-window, and deterministic-ranking helpers without changing its
  public request/report.

- Added a first-class terminal-exon sequence-specific RT-primer-pool workflow.
  The shared engine now composes Splicing Expert transcript resolution,
  strand-aware mature-transcript/exon mapping, fixed-adapter and multiplex
  complementarity ranking, persisted provenance, undo, and JSON export. PCR
  Designer, shared shell/direct CLI, MCP introspection, JavaScript, Lua, and
  Python operation adapters consume the same
  `gentle.terminal_exon_rt_primer_pool.v1` report; Tm remains descriptive and
  global specificity is explicitly not inferred.

- Added non-mutating construct-reasoning inspection across project-sequence
  selections and homogeneous gene-set resolutions. The shared
  `collections run construct-reasoning-inspection` route returns per-member
  actions, fingerprints, successes, and failures without persisting temporary
  graph snapshots or publishing dangling graph report IDs.

- Added a content-bound `gentle.transcript_assay_cdna_similarity_map.v1`
  planning seam for transcript-assay design. Whole-cDNA similarity intervals
  can now rank bounded Primer3 records toward lower-risk regions while exact
  operation digests, template coordinates, BLAST/database provenance, emitted
  Primer3 constraints, and final complete-cDNA/genomic specificity gates stay
  explicit and unchanged.
- Added collapsed structured evidence panes beside construct-reasoning
  inspection buttons. The GUI now presents action geometry, readable dotplot
  mode, rationale, discrete evidence/context rows, repeat-family provenance,
  and task-severity fields without flattening optional values into prose;
  legacy provenance fallback and stale-action refusal remain unchanged.

- Added deterministic related-sequence group-target primer design and explicit
  primer-pair Pareto alternatives. `DesignPrimerGroupTarget` aligns loaded
  members, derives exact conserved binding intervals, bounds work before
  Primer3, and reports a complete pair-by-member product matrix under strict or
  explicit best-effort coverage. Existing primer and transcript-panel reports
  now expose bounded non-dominated tradeoffs without changing their selection
  or specificity gates.

- Extended Primer-BLAST-style primer planning without surrendering GENtle's
  transcript-aware contracts. GENtle can now build content-bound whole-cDNA
  similarity maps, combine variant/repeat/similarity intervals as explicit
  hard or soft pre-Primer3 evidence, emit reviewed chemistry settings, cap
  junction-side matches, realign complete primers after BLAST, retain narrowly
  reviewed off-target products, render stored alignments as HTML, and propose
  content-bound replacements after panel specificity failure. Advisory design
  evidence never substitutes for the final complete-cDNA/genomic gates.

- Grounded Agent Assistant and `agents ask` requests in engine-owned
  introspection whenever project context is enabled. The additive
  `gentle.agent_introspection_context.v1` extension carries a bounded current
  fact projection, complete fact-type counts, explicit truncation metadata,
  and read-only retrieval routes without granting providers filesystem access
  or changing execution authority. The agent-interface tutorial now includes a
  Mistral-compatible facts-before-action exercise that tests open-world absence
  reasoning on a deterministic local sequence.

- Expanded the inner Agent Assistant tutorial with two introspection-complete,
  offline exercises and documented the precise determinism boundary between
  variable model proposals and GENtle's fixed parser, readiness, execution,
  and effect-verification contracts. The testing guide now specifies an
  opt-in, self-skipping native-Mistral conformance routine.

- Grouped all project lifecycle commands at the top of the File menu, with
  `Save Project...` immediately beside `Close Project` before sequence and
  resource commands.

- Made the JASPAR Expert catalog sortable by every displayed column, show
  cached species names directly, and fetch missing metadata for only the
  on-screen rows through its existing background operation. Always-visible
  horizontal and vertical scroll bars keep the full controls, motif details,
  sequence logo, matrices, and score distributions reachable at smaller
  window sizes.

- Replaced the sparse standard macOS About panel with one localized,
  GENtle-owned About experience shared by the native application menu and the
  in-app Help menu. The compact Finder-inspired card now presents version,
  build, authorship, preview status, and licence; `More Info …` expands to the
  mascot, reproducibility context, project/documentation/credit links, and a
  copyable build identity. The German native menu item is now `Über GENtle`.

- Added an orthogonal transcript-assay coverage-universe contract. Existing
  all-annotation plans remain byte-compatible, while explicit transcript sets
  and content-bound UniProt projection audits can define mandatory targets for
  any existing panel objective. UniProt universes now bind the complete named
  protein-isoform inventory rather than multiplying mapped transcript rows;
  reports group every mapping beneath its protein target, retain annotation/CDS
  discordance, and record coverage separately from distinction. Unmapped named
  isoforms stop before Primer3 without weakening `require_all` or `best_effort`.

- Replaced the PCR skill's hidden 30-minute isoform-workflow ceiling with an
  approval-bound `execution_timeout_secs` slot. Single-study design defaults to
  a two-hour allowance and ordered multi-study execution to eight hours, while
  callers may review and override either ceiling before approval; documentation
  now distinguishes this wrapper policy from Primer3 runtime and records the
  workload dimensions that affect transcript-panel cost.

- Added manifest-driven multi-gene publication bundles. A shared resolved
  report now drives a responsive PARK7-style HTML overview and a printable
  companion containing the complete primer-pair list and all selected locus,
  primer/exon-map, and virtual-gel figures. The new
  `gentle_publication_report` binary makes this reusable for arbitrary gene
  sets and emits an auditable generation receipt. Bundle filenames are confined
  to the output root, colliding asset basenames are rejected, and optional PDF
  tool failures identify their executable override without leaving a false PDF
  link in the resolved report.

This file is the auditable record of completed work that used to accumulate in
[`roadmap.md`](roadmap.md). Keep entries brief and outcome-focused. Durable
implementation rules belong in [`decisions.md`](decisions.md), not here.

Maintenance rule:

- Do not add completed-work bullets to [`roadmap.md`](roadmap.md).
- At the end of each session, move completed roadmap items here.
- Prefer a date and one short outcome sentence. Include command names,
  document names, schemas, or feature names only when they help a reader
understand what changed.

## 2026-08-09

- Hardened ClawBio CLI supervision. On POSIX, `SIGUSR1` now reaches only the
  exact active GENtle child, while `SIGINT`/`SIGTERM` stop the child-owned
  execution group with bounded escalation. Timeout and interruption remain
  failed wrapper outcomes even if a child handles the signal and exits zero;
  receipts retain the actual child status, signal events, complete output, and
  hashes as orchestration provenance rather than scientific claims.
- Hardened bounded Primer3 integration around platform and signal behavior.
  Engine-owned executable discovery now honors Windows `PATHEXT`, and a
  progress-capable binary is not sent `SIGUSR1` until it has emitted its first
  valid progress row. Repository-owned fake-binary protocol tests always run;
  only separately marked live-binary integration tests remain opt-in.
- Per-gene continuation now refuses approved operations with declared external
  output paths before executing anything, because filesystem writes cannot be
  revoked by its in-memory rollback. Rolling back a failed gene also restores
  `last_checkpoint_manifest` to the last valid whole-gene boundary. Continue
  mode now writes checkpoints only at complete-gene boundaries, so a
  rolled-back partial gene emits no new discoverable checkpoint; abort mode
  retains its existing per-operation recovery checkpoints.

## 2026-08-06

- A primer-pair QC status of `warning` is no longer reported as a failed
  `critical_qc` gate. A failed gate blocks regardless of whether the policy
  requires it, so treating a caution as a failure made it unwaivable and
  indistinguishable from real QC failure. `warning` is now `incomplete`: it
  still blocks under the default policy, and a recorded `fail` still blocks
  unconditionally.
- Regenerated the PATZ1 transcript-assay tutorial artifacts, which had not been
  refreshed when bounded Primer3 search landed and therefore no longer matched
  the generator.
- GENtle now chooses between multiple installed Primer3 builds instead of
  taking the first `PATH` match: `--progress` support wins first, then the
  highest version, then the newest file. An explicitly configured path is still
  used unchanged, and `primers preflight` reports every candidate and why one
  was selected.
- Primer3 protocol, progress, runtime-reduction, preflight, and candidate
  ranking tests use repository-owned deterministic fake executables and run
  without a system Primer3 installation. The separate tests that intentionally
  exercise a genuine installed binary remain explicit opt-in integration
  checks.
- Added opt-in per-gene continuation to approved multi-gene study batches.
  `primers execute-gene-isoform-study-workflow-batch --on-gene-failure continue`
  rolls a failing gene back to its own boundary and carries on with the other
  precomputed genes, moving atomicity from the batch down to the gene rather
  than removing it. Approval, verification, and coverage semantics are
  unchanged, the reuse checkpoint freezes at the last contiguous prefix, and
  the report names every skipped gene.
- Isoform-assay dossiers can now be published while a study is incomplete.
  Genes carry a `resolved`/`pending`/`unresolved` status with a reason, derived
  from the gene's own contents when the request does not declare one. The index
  lists what is outstanding and every incomplete page opens with a notice, so
  biologists can start on the finished genes and the dossier is regenerated
  once the remaining results arrive.

## 2026-08-05

- Added sequence-aware bounded Primer3 search plans for transcript-assay
  operations. Endpoint and junction targets now expose annotation-anchored
  allowed/rejected intervals, deterministic work budgets, exact Primer3 fields,
  and fail-before-launch statuses; suspicious progress can stop one bounded
  record without weakening strict coverage. ClawBio process supervision now
  forwards status/termination signals to the active process group and records
  exact child output hashes.
- Made approved multi-gene isoform-assay execution truly atomic in memory and
  on disk, added gene/workflow/operation ordinals to runtime/SIGUSR1 status,
  and added opt-in content-bound operation checkpoints with a separately
  approved exact-prefix reuse proposal. Compatible work can now survive a
  later batch failure or an appended gene without automatic reuse or weakened
  approval/provenance checks.
- Added deterministic pre-Primer3 feasibility for endpoint transcript end
  matrices. Exact approved operations now expose stable reaction/class ids and
  structural blockers through `gentle.transcript_assay_panel_feasibility.v1`;
  strict single- or multi-study execution rejects impossible geometry before
  any operation runs, while explicit best-effort reports retain skipped
  reactions. Endpoint searches now request four candidate pairs per reaction
  by default instead of the previous tenfold candidate surplus.

## 2026-08-04

- Integrated Primer3 bounded-progress reporting across engine, GUI, and CLI.
  GENtle now detects the `--progress` capability from cached help probes,
  streams candidate-work counters when available, safely requests refreshed
  Unix status with `SIGUSR1`, preserves cancellation without an accidental
  `auto` fallback, and runs older Primer3 binaries directly in legacy mode.

## 2026-08-03

- Made isoform-study normalization directly reusable by emitting the normalized
  request as stdout, and added content-bound multi-study workflow composition
  and execution. The `gentle-pcr-primer-design` 0.6.0 descriptor now offers one
  second-stage approval that prevalidates the complete ordered set before any
  shared project state is changed.
- Closed the two-stage isoform-study approval chain: plans now bind the exact
  emitted workflow bytes, `primers execute-gene-isoform-study-workflow`
  verifies both workflow and ordered-operation digests before execution, and
  the PCR ClawBio descriptor exposes distinct read-only normalization, first
  approval, and second approval routes.
- Added non-destructive `genomes adopt-blast-resource` metadata repair for an
  existing transcriptome-cDNA BLAST database. Adoption validates the expected
  index fingerprint, records or verifies the prepared FASTA identity, can
  compare all accession/length rows, regenerates only subject annotations,
  preserves legacy manifests, and never invokes `makeblastdb`; ordinary
  preparation now protects content-bound indices from implicit replacement.
- Expanded gene isoform-assay dossier publication with GENtle-owned quality-
  assurance and provenance projections in every standard profile. Gene pages
  now expose exact readiness gates, evidence ids, blockers/warnings, and a
  role-described plan-to-handoff-to-order chain with schema/report identities,
  policy/panel bindings, source paths, and SHA-256 digests.
- Dossier PDF rendering now discovers Chromium, Google Chrome, and Microsoft
  Edge through `PATH` and conventional Linux, macOS, and Windows locations,
  while preserving `GENTLE_BROWSER_BIN` as the authoritative override.
- Content-bound legacy gene-set publication bundles with a deterministic
  `bundle-manifest.json`, while preserving the v1 generation schema and
  `copied_files[]`. Added regression coverage for browser-PDF failure,
  presentation profile/block filtering, iterative study overrides, and v1
  projection-option rejection.
- Added deterministic gene isoform-assay study planning, two digest-bound
  approval subjects (normalized evidence/policy inputs and exact ordered design
  operations), readiness-bound oligo-order creation, and canonical multi-gene
  dossier publication. GENtle now owns the meta-page, one HTML page per gene,
  declared presentation blocks/profiles, print HTML/PDF projection, and
  provenance-rich order sheets; OpenClaw may select and host those projections
  but cannot rewrite their scientific content.
- Added a fail-closed two-phase approval contract for confirmation-gated
  ClawBio delegation. Scientific PCR routes now emit a content-addressed
  `gentle.clawbio_execution_proposal.v1`, execute only an exact caller-approved
  proposal, and reject route/runtime/state/input/backend drift; list, show, and
  preflight diagnostics remain automatic, while direct structured GENtle use
  is unchanged. Local Cargo launchers are pinned to their built executable and
  OCI execution requires an immutable image digest before approval.

## 2026-08-02

- Added a bounded external-primer handoff mode to the co-shipped ClawBio skill.
  Explicit qPCR/endpoint-PCR records now retain GENtle's native provenance and
  return an input-to-result ID/hash join with verified report, state, runtime,
  and artifact bindings; cloning/sequencing oligos remain typed
  `not_submitted`, with no discussion-state or interpretation logic added to
  the engine.
- Declared the existing serial-arrangement gel renderer as a
  context-agnostic collection `arrange` while retaining container `combine`
  semantics and explicit physical-pool rejections for project selections and
  logical gene sets. Arrangement gel preview/export now shares the engine's
  serial-mode readiness check, preserves tested lane order, and stays disabled
  for plate layouts without affecting their rack/plate actions.
- Added explicit physical pooling for complete resolved gene sets.
  `CreateGeneSetPool` and `gene-sets create-pool` preview a deterministic,
  fingerprint-locked plan from one exclusive singleton source container per
  member; apply derives lineage-linked aliquot sequences into one physical
  pool while retaining every source tube. The Gene Set Inspector uses the same
  typed command and exposes per-member source/output provenance.
- Hardened physical-container pool export in the Lineage GUI: export work now
  runs in a detached engine snapshot rather than serializing sequence JSON
  under the egui-thread write lock, and completion requires a typed,
  schema-checked `gentle.collection_pool_export.v1` report whose container and
  artifact path match the request. Missing container members are also covered
  by an atomic no-output regression test.
- Prevented in-silico `Selection` containers from being exported as physical
  `gentle.pool.v1` artifacts; engine and Lineage GUI readiness now accept only
  exhaustive physical singleton or pool containers.
- Fixed Gene Set Inspector collection readiness so promoter forms follow actual
  gene-set selection changes, unsupported GUI context modes are identified as
  adapter gaps rather than biological policy rejections, and the selected
  operation reuses its already computed readiness.

## 2026-08-01

- Added restriction digestion as the first mutating collection `map`
  operation. Shared shell/CLI/MCP and the Gene Set Inspector now preview exact
  source- and namespace-locked fragment plans before apply; successful apply
  preserves circular/linear direct-digest semantics, creates individually
  openable sequences with per-source lineage, and never implies a pooled
  physical container.
- Made direct and collection digest convergence independent of machine load by
  replacing the per-enzyme wall-clock deadline with deterministic fragment,
  round, repeated-state, and progress guards.
- Standardized A5-A7 collection output wrappers as additive, version-tolerant
  report contracts while keeping their member-binding request records strict.
- Added the direct shared-shell `digest` parity route, made unknown restriction
  enzymes fail closed before mutation, and removed digest progress text from
  stdout so MCP stdio framing remains valid.
- Added a local-only E-MTAB-14704 Clariom D probe-set activity renderer and
  runbook. It emits provenance-bound descriptive raw-PM tables and figures with
  explicit E1/E2/E3 pairing, bounded streamed inputs, deterministic SVG/PDF
  metadata, and clear non-significance/non-specificity/non-isoform caveats.
- Added non-mutating TFBS hit scanning as a collection `map` operation.
  Project sequences and explicitly sequence-bound gene-set members now share
  direct-scan motif/threshold semantics across engine, shell/CLI/MCP, and the
  Gene Set Inspector; the portable wrapper distinguishes retained aggregate
  hits from complete counts when members fail or a cap stops motif scanning.
- Added physical-container pool export as an atomic collection `combine`.
  `ExportPoolCollection` and `collections run export-pool` reuse the existing
  `gentle.pool.v1` exporter, return `gentle.collection_pool_export.v1`
  provenance and membership-lock metadata, and reject non-exclusive containers
  before writing because the pool artifact cannot preserve incomplete contents.
  The Lineage container row calls the same shell/engine route, and the lift
  registry now also declares the already-shipped container pooled-gel action.
- Added restriction-site scanning as the second collection `map` operation.
  Project sequences and explicitly sequence-bound gene-set members now use the
  same topology-aware scanner through engine, shell/CLI/MCP, and the Gene Set
  Inspector; a portable wrapper carries child reports and aggregate site
  counts without dangling report ids.
- Replaced the Gene Set Inspector's one-off primer action with a
  registry-projected collection-operation launcher. Typed promoter-cohort and
  primer-specificity adapters now share detached execution and subject checks,
  while physical pool, pooled-gel, and serial-arrangement policies remain
  visible with their canonical typed rejection reasons.

## 2026-07-31

- Added all five gene-set source kinds to `Genome > Gene Set Inspector...`:
  the GUI now constructs the same `ResolveGeneSet` payloads as the shell,
  commits resolutions through detached engine execution, auto-selects current
  results, and exposes resolved/unresolved membership, provenance, warnings,
  and portable JSON before existing collection operations.
- Added depth-aware, auditable sequence-representation verdicts to
  `gentle.rna_allele_hash_screen.v3`: phased variant/transcript/block summaries
  and the gene roll-up now distinguish preferred, balanced, low-depth, and
  unphased evidence while retaining RNA-report expression depth as labeled
  context rather than a causation or significance claim.
- Added an in-figure transcript-feature legend to the shared gene-locus SVG
  renderer, explicitly identifying exon, intron, CDS, annotated start-codon,
  and annotated stop-codon glyphs.
- Split Clariom evidence in the shared gene-locus report into separately scaled
  abundance/splice-geometry and differential-activity lanes. Tables may now
  supply both `log2_mean_*` and `log2_*_minus_*` columns; the figure retains
  explicit cautions that array signal prioritizes assay regions but does not
  prove PCR-primer binding, significance, or isoform identity.
- Added project-backed and Salmon-aware read sourcing to
  `rna-reads allele-hash-screen`: target-gene accepted reads can flow directly
  from a persisted RNA-read report, Salmon outputs select IDs from explicit
  sequence files, and `gentle.rna_allele_hash_screen.v2` records per-read and
  aggregate source provenance while retaining v1 deserialization.
- Added offline, assembly-aware primer/probe variation screening through
  `ScreenPrimerVariants` and `primers screen-variants`: one local VCF pass now
  emits source-fingerprinted, strand-aware, handoff-ready
  `gentle.primer_variant_evidence.v1` reports without treating unknown allele
  frequency or incompatible reference evidence as clear. An optional
  frequency-gated follow-up now proposes newly identified mixed-IUPAC primer
  pairs for frequent simple SNPs, preserves selected VCF annotations for
  inspection, and reuses IUPAC-aware cDNA matching without treating a synthesis
  mixture as a genotype or splice-effect correction.
- Linked reasoning-guided dotplots back to their recommending
  construct-reasoning graphs in interactive and exported lineage views;
  lineage details now distinguish verified, failed, unknown, and absent/manual
  citations and expose the stored rationale and driving evidence ids.
- Kept tutorial evidence digests byte-exact while tolerating Windows CRLF
  checkouts: the tutorial runner stages an LF-normalized UTF-8 input only when
  that normalized content exactly matches the pinned SHA-256 digest. Windows
  builds also no longer import the Unix-only SIGUSR1 `AtomicBool`.
- Persisted construct-reasoning inspection provenance on stored dotplots:
  payloads and list summaries now retain operation/run identity plus a verified
  graph/action/fact/evidence/request citation, while manual absence remains
  distinct from failed or unknown verification.
- Preserved an existing gene-set resolution's operation/run identity when
  deriving a promoter cohort, preventing one logical source set from appearing
  as duplicate lineage nodes. Anonymous inline resolutions still receive a
  stable identity when first persisted.
- Fixed generic `r_oligo` platform annotation discovery: the helper now
  resolves the selected pdInfo package's bundled SQLite database instead of
  hard-coding Clariom D Human, allows only an unambiguous package-local
  fallback, and refuses probeset runs that would otherwise emit coordinate
  columns without annotations.
- Completed Promoter design parity for ortholog ambiguity review and export:
  the GUI now exposes `reject|first|preserve`, renders preserved mapping
  candidates, and exports cohort/comparison JSON through the shared engine
  operations.
- Added an additive, fail-closed normalized CUT&RUN contract to
  `gentle.ortholog_promoter_comparison.v1`. Quantitative cross-species rows
  require explicit method, unit, shared reference, provenance, one
  source-bound value per resolved promoter, and matching selected evidence;
  otherwise GENtle retains qualitative states and reports `not_comparable`
  without comparing raw intensity.
- Added the descriptor-only `gentle-pcr-primer-design` OpenClaw/ClawBio skill.
  Explicit PCR, RT-PCR, SYBR, TaqMan, imported-primer, and specificity intents
  now delegate to the existing `gentle-cloning` runtime while preserving
  GENtle-owned evidence, provenance, and review states. Delegation now fails
  closed on a stale or partial two-skill deployment and records only the
  selected route, without claiming reproducible natural-language routing.
- Added the topic-neutral `gentle.clawbio_execution_manifest.v1` provider
  receipt. It content-binds requests and immutable inputs before and after a
  run, mutable state snapshots, wrapper/runtime identity, command outputs,
  native result identifiers/statuses, and artifacts while keeping execution
  outcome separate from biological verdicts. PCR routes also use exact
  pointer/hash/projection-backed `gentle.clawbio_claim_ledger.v1` attribution;
  caller-supplied source prefixes cannot claim GENtle authority.
- Added an optional, post-run discussion-moderation handoff adapter. It verifies
  exact caller input and provider-result hashes, emits a typed computed-evidence
  object plus a canonical analysis-run intake, and leaves permission checks,
  recording, replay, freshness, and interpretation to the caller-owned ledger.
- Extended the experimental assay handoff with a GENtle-attributed virtual gel
  containing one comparable lane per primer pair and explicit predicted-empty
  lanes.
  Transcript-panel requests now preserve oligo-dT priming, cap-dependent
  5-prime capture, and supplied completeness evidence as separate `[input]`
  claims, without promoting those preparation assumptions to GENtle findings.

## 2026-07-30

- Added `rna-reads verify-dexseq` and
  `gentle.rna_read_dexseq_verification.v1`: GENtle now exports a matched
  flattened GFF/count pair, performs bounded non-downloading R/DEXSeq
  preflight, and can prove the pair loads through the real
  `DEXSeqDataSetFromHTSeq()` contract.
- Hardened preserved ortholog ambiguity reports: species-only mappings now
  retain organism context without fabricating a genome id, candidate labels
  expose provider source when available, and the closed-policy v1
  compatibility boundary is explicit.
- Added `ambiguity_policy=preserve` for offline ortholog promoter resolution.
  It leaves ambiguous targets unresolved while retaining ordered, structured
  mapping candidates with portable context and provider evidence; existing
  `reject` and `first` behavior is unchanged.
- Bound the existing offline ortholog resource and promoter-cohort reports to
  portable biological-context registries. Orthology type/confidence now use
  open typed vocabularies that preserve provider-specific legacy strings;
  context organism/genome conflicts fail before catalog access, while
  symbol-only mappings and the existing `reject`/`first` ambiguity behavior
  remain compatible.
- Unified protein-gel, peptide-gel, 2D-gel, and isoform-evidence molecular
  weights on the amino-acid residue-mass model. All routes now add one terminal
  water per chain and reject ambiguous residues instead of silently producing
  a zero-mass estimate; gel errors consistently identify empty sequences or
  list the ambiguous/unsupported residues.
- Made probe-region R readiness use the same explicit, repeatable
  `--r-library-path` contract as the generated oligo/affy helpers. Preflight now
  batches and bounds direct package/version checks, records effective
  `.libPaths()`, and names those paths in diagnostics when sandbox, user, and
  system R installations disagree.
- Added report-owned biological-context registries for resolved gene sets and
  portable collection-operation reports. Promoter derivation and gene-set
  primer-specificity mapping now require a homogeneous context matching the
  requested genome and reject missing or mixed contexts before coordinate or
  BLAST work, while legacy report-level context fields remain readable.
- Added `Genome > Gene Set Inspector...` as the first prominent collection GUI
  consumer. It requires explicit persisted gene-to-primer-report bindings,
  executes the shared `collections run primer-specificity` route in a detached
  background snapshot, and keeps member execution outcomes distinct from child
  reports' biological specificity verdicts.
- Added primer-selection provenance to `gentle.primer_design_report.v1`:
  exact additive score terms, bounded deterministic evaluated near misses with
  explicit capture completeness, and report-fingerprinted construct-reasoning
  decisions whose rejected intervals reuse the existing linear-map overlay.
- Extended the same descriptive provenance to
  `gentle.qpcr_design_report.v1`: exact primer/probe score terms, bounded
  evaluated assay near misses, report-bound reasoning graphs, and cached
  score/capture inspection in the existing PCR Designer. Region-level repeat,
  variant, and paralogue exclusions remain explicitly unconsulted and are
  recorded as `not_run` rather than implied clear.
- Added engine-owned collection `map` execution for primer specificity through
  `AssessPrimerPairSpecificityCollection` and
  `collections run primer-specificity`. Persisted gene sets use explicit
  member-to-primer-report bindings, project sequences may resolve one unique
  report, and the aggregate `gentle.collection_operation.v1` keeps execution
  failures separate from each child report's biological verdict.
- Upgraded the gene isoform evidence ledger to
  `gentle.gene_isoform_evidence.v2`, retaining every contrast/source
  measurement without cross-unit magnitude selection, adding typed PSR/JUC
  provenance, ambiguity-safe protein identity/mass, and rule-based assay
  recommendation tiers. Added pure-read
  `gentle.gene_transcript_assay_routine.v1` composition for existing evidence,
  panel, specificity, and order-table artifacts, plus R/oligo package-version
  and input-fingerprint inspection without package installation.
- Added the first engine-owned collection-lifting slice:
  `gentle.collection_lift_policy_registry.v1` capability policies,
  subject-aware canonical membership fingerprints, and portable
  `gentle.collection_operation.v1` reports. Gene-set promoter cohorts now link
  source genes to derived windows through typed per-member rows and persisted
  report ids, while logical gene sets are explicitly rejected as implicit
  physical pools or gel lanes.

## 2026-07-29

- Preserved gene-level identity for prepared transcriptome BLAST resources in
  checksum-bound `gentle.blast_subject_annotation_index.v1` sidecars. Primer
  specificity hits and products now retain optional Ensembl/GenBank subject
  annotation in compact and full reports without allowing gene identity to
  affect transcript-level intended-product matching; readiness ceilings also
  record their human-readable derivation.
- Completed the primer-specificity follow-up with narrow Ensembl transcript
  version matching, persisted panel genome anchors, distinct external
  `pass|specificity_fail|not_assessed|incomplete` outcomes, panel-inherited
  readiness ceilings plus long-product warnings, indexed/compact v4 reports,
  complete shell help, and catalog-driven preparation/inspection of
  `transcriptome_cdna` BLAST resources.
- Hardened persisted primer-specificity artifacts as
  `gentle.primer_specificity_report.v3`: each HSP now exposes unaligned query
  bases and effective mismatches, expected products are explicit per genomic or
  transcriptome target space, intended isoform coverage and
  junction/genomic-carryover have independent four-state results, and cited
  design pairs are content-fingerprinted against the oligos actually assessed.
- Completed engine-owned feature-record curation with strict Split/Merge
  preview/apply variants in `gentle.feature_record_curation.v2`, shared
  `features split` / `features merge` routes, and matching Feature Editor tabs.
  Split is limited to an internal boundary in one exact simple feature; Merge
  requires exactly touching simple records with identical kind, strand, and
  ordered qualifiers. Both use per-feature plus annotation-state locks and
  normal undo/redo, while ambiguous topology or metadata reconciliation is
  rejected.
- Extended primer-specificity v2 with transcript-set intended products,
  GENtle-derived target geometry for imported pairs, effective mismatches for
  partial HSPs, and independent genomic-carryover and transcriptome-specificity
  readiness gates that retain failed or incomplete evidence.

## 2026-07-28

- Added engine-owned `gentle.feature_record_curation.v1` Create/Delete
  preview/apply operations, portable annotation-state locks, ordered
  duplicate/valueless qualifier preservation, full undo/redo, shared
  `features create` / `features delete` routes, and Create/Delete tabs in the
  DNA Feature Editor. Overlap and shared INSDC identifiers remain
  informational review evidence; automatic propagation remains deferred.
- Extended the engine-owned Feature Location Editor and shared
  `features edit-location` route to edit one exact child in flat
  `Join`/`Order` annotations, including outer-complement reverse-strand
  compounds. Segment edits preserve the complete location topology and
  qualifiers, report stored versus biological segment numbering, expose typed
  overlap/direction/CDS-length warnings and cross-role related boundaries, and
  keep nested, fuzzy, topology-changing, and circular cross-origin locations
  read-only.
- Corrected primer-specificity interpretation for minus-strand targets and
  per-HSP query coverage, and made exhaustive BLAST coverage enforceable by
  setting a subject cap above the validated database sequence count; direct and
  handoff reports now expose `search_completeness`, and unproven searches yield
  `incomplete` rather than a specificity pass.
- Promoted the existing primer-specificity report family to persisted
  `gentle.primer_specificity_report.v2` computational artifacts with stable
  ids, operation/run and database provenance, explicit independent assessment
  states, design-report citations, CLI inspection/export, and GUI lineage
  reopening; raw primer pairs remain explicit about absent selection
  provenance.

## 2026-07-27

- Added engine-owned exact feature-location preview/apply operations,
  `features edit-location`, optimistic-concurrency fingerprints, undo/redo,
  adapter introspection, and a GUI Feature Location Editor reached from the
  command palette and feature context menus; complex/fuzzy locations remain
  explicitly read-only and related annotations are review-only.

## 2026-07-23

- Added provenance-aware JSON/TSV external primer-pair import with stable
  sequence-derived identities, duplicate-source retention, shared Tm/GC/oligo-QC
  metrics, cDNA product and genomic-carryover evaluation, optional prepared-genome
  specificity and product gels, and an explicit rule that provider claims never
  count as transcript-coverage or specificity evidence.

## 2026-07-22

- Added a policy-aware PrimerBank adapter with live or saved-HTML search,
  typed `gentle.primerbank_search.v1` records, exact pair lookup, JSON export,
  and an explicit continuation into GENtle's transcript-aware cDNA PCR test
  without inferring genomic specificity or experimental validation; typed
  requested-versus-observed species checks are available to every adapter and
  cDNA continuation now refuses mismatched or unresolved catalog records and
  independently rejects a conflicting project-sequence organism annotation.
- Added deterministic per-panel experimental assay handoffs with canonical
  oligo/pair identities, automatic pair-linked cDNA tests, exact product-
  sequence classes, versioned readiness gates, optional provenance-bound
  variant evidence, procurement formulation projection, and JSON/TSV exports
  across the shared engine, PCR Designer, shell/CLI, and generic adapter
  surfaces.
- Added experimental-practicality provenance to multi-transcript primer panels:
  annotation-only routine common-region screening, separate PSR/JUC support,
  preferred-versus-allowed product ranges, explicit long-range fallbacks,
  bounded alternatives, and conservative endpoint-gel interpretation now flow
  through the shared engine, shell/CLI, GUI, and JSON adapters.
- Rejected opposing PCR primer pairs whose binding footprints overlap, aligned
  transcript-panel selection with its product matrix, and regenerated the PATZ1
  SYBR tutorial so every selected pair detects its own design transcript;
  legacy summary enrichment now preserves unknown origin and restores only
  explicitly recorded requested-junction overlaps.
- Completed the offline CUT&RUN V1-V3 release path with a built-in synthetic
  TP73 catalog dataset, prepared BED and paired-read evidence, zero-flank
  mapping against an imported genome anchor, custom catalog/cache replay in V3
  and its GUI inspector, and deterministic JSON/TSV workflow artifacts.
- Completed the `gentle.primer_pair_summary.v2` communication contract with
  canonical design-amplicon coordinates/length, explicit separation from
  cross-transcript predicted products, blank-specificity normalization, and
  exact regression checks across primer, QC, backend, junction, and export
  fields.
- Extended `gentle.primer_pair_summary.v2` with immutable sequence-derived
  primer identities, transcript/exon display labels, aliases, typed origin and
  selection provenance, three distinct junction semantics, and structured
  Clariom region evidence that never implies exact probe-sequence reuse.
- Updated the Rust dependency graph to the latest Rust 1.98-beta-compatible
  releases, including Deno 0.408, V8 149, Temporal 0.2.4, ICU 2.2, Serde
  1.0.229, and libc 0.2.189; direct requirements now use normal semver ranges
  instead of exact or wildcard pins, while the embedded JavaScript shell keeps
  its synchronous interface through an owned Tokio runtime and a V8 execution
  smoke test.
- Added `blastdbcmd`-backed BLAST database component inspection and
  content-bound primer-specificity handoffs; specificity now searches all
  database subjects, normalizes prepared-FASTA identifiers, distinguishes
  genomic carryover from transcriptome cross-amplification, and models
  junction-spanning assays without inventing contiguous genomic targets.
- Added `gentle.primer_pair_summary.v1` to selected transcript-panel assays so
  CLI/MCP/API and JSON exports carry one joined, assay-neutral communication
  view of the exact primer sequences, lengths, melting temperatures, GC values,
  pair delta-Tm, transcript products, QC reasons, specificity state, and tool
  provenance without guessing a PCR annealing temperature.

## 2026-07-21

- Made serial/pool-gel SVG labels collision-aware: dense gels now keep short
  lane names horizontal, wrap isolated long names, angle difficult dense names,
  and move band text that cannot fit between lanes into the fragment table;
  callers can override both lane- and band-label placement explicitly.
- Added redundant transcript-isoform identity cues to serial/pool gels:
  recognized Ensembl or RefSeq transcript accessions now retain a deterministic
  color and relative marker position across lanes and merged bands, while a
  side legend and fragment rows provide an `O`/`I` binary text fallback.
- Added aggregate transcript-assay panel specificity handoff/finalization:
  mechanical wrappers now return per-command exit and output identities, while
  GENtle alone classifies complete all-assay evidence as `pass`,
  `specificity_fail`, or `incomplete` and atomically persists only `pass`.
- Added deterministic `primers specificity-plan` / `specificity-import`
  handoffs so external schedulers can own BLAST process completion while GENtle
  retains query, database, policy, and result-interpretation provenance; the
  shared `PreparePrimerPairSpecificityHandoff` /
  `ImportPrimerPairSpecificityHandoff` operations expose the same lifecycle to
  CLI workflows, MCP, JavaScript, and Lua, and the transcript-assay panel can
  also request inline report-only or required-pass local BLAST confirmation.
- Added a deterministic offline CLI tutorial for the synthetic PATZ1 endpoint
  RT-PCR and primer-only SYBR transcript panels, including full-operation
  `@FILE` execution, persisted report inspection/export, oligo-dT interpretation
  safeguards, and required Clariom-style junction review.
- Added a full-operation JSON form to
  `primers design-transcript-assay-panel`, so the complete
  `DesignTranscriptAssayPanel` request is reusable unchanged across direct CLI,
  shared shell, workflows, MCP, JavaScript, and Lua while retaining the concise
  flag form.

## 2026-07-20

- Added GUI parity for `gentle.transcript_assay_panel.v2` in PCR Designer:
  endpoint/SYBR/TaqMan setup from Splicing Expert, strict or best-effort
  coverage, explicit/Clariom junction inputs, background execution, persisted
  report load/export, and graphical product/band matrices with per-cell
  oligo-dT reach interpretation.
- Extended `gentle.transcript_assay_panel.v2` with backward-compatible endpoint
  RT-PCR, primer-only SYBR, and TaqMan modes; annotated first x terminal reaction
  and band matrices; required/preferred explicit or Clariom JUC targets;
  Primer3 junction tags; order/provenance/specificity rows; oligo-dT cautions;
  per-cell annotation-derived RT-reach assessments with optional experimental
  risk thresholds; and a deterministic synthetic PATZ1 endpoint-plus-SYBR
  workflow.
- Added persisted `gentle.transcript_assay_panel.v2` reports with exact mature-cDNA
  equivalence classes, transcript-by-assay product matrices, strict-by-default
  pan-transcript coverage, explicit best-effort partials, and shared shell
  list/show/export lifecycle commands; the legacy fixed-component panel now
  uses canonical primer candidate generation and labels byte-identical models.
- Kept portable ClawBio parity coverage active on Windows while excluding the
  two POSIX-only Bash launcher tests that cannot run from a normal CRLF
  checkout; macOS and Linux retain launcher coverage.

## 2026-07-19

- Made `gentle_examples_docs tutorial-check` treat CRLF and LF as equivalent
  for known text artifacts while retaining byte-exact drift checks for binary
  files, so normal Windows checkouts pass tutorial validation.
- Added a graphical `Locus figure` tab to Splicing Expert that composes the
  shared `gentle.gene_locus_evidence_display.v1` report, previews it in the
  GUI, reports local-resource readiness, supports relocation and JSON/SVG/PDF
  export, and continues selected junction/qPCR evidence into existing assay
  design and reviewed oligo-order paths.
- Added a deterministic offline PATZ1 locus-composer workflow, synthetic
  PSR/JUC and grouped-occupancy fixtures, generated tutorial/SVG, structural
  visual-regression checks, and focused Puffin scopes for composition costs.
- Standardized GUI evidence wording around observed evidence, candidate
  association, design constraint, not evaluated, and unresolved evidence
  without changing the stable wire spellings.

## 2026-07-18

- Upgraded the optional embedded Lua adapter to `mlua 0.12`/`mlua-sys 0.11`
  and the current vendored Lua and LuaJIT sources while retaining Lua 5.4,
  Serde conversion, and the shared engine/shell behavior.
- Distinguished generated probe-region evidence from prepared array datasets
  in live feature-tree grouping, and pinned backward-compatible inspection of
  coordinate-less v1 output while continuing to refuse its genome projection.
- Added the engine-owned `gentle.protein_expression_requirements.v1` contract
  and `planning protein-expression-handoff --requirements`, with strict schema
  validation, topic-specific missing-question resolution, report replay, and an
  explicit guard that withholds provider actions for in-house-only work.
- Made protein-expression provider handoffs readiness-consistent across text
  and JSON, and replaced tutorial-protein commands with product-specific,
  project-referenced GeneArt review drafts when a selected DNA/protein product
  is ready for provider preflight.

## 2026-07-14

- Extended `gentle.gene_locus_evidence_display.v1` with a generic,
  anchor-verified Clariom-style probe-effect overlay that preserves PSR
  intervals, JUC spans, contrast values, source provenance, and one
  zero-centered SVG scale; the committed 24-row PATZ1 fixture now substitutes
  for the local CEL-analysis tree in CI.

## 2026-07-13 - `v0.1.0-internal.10` release candidate

- Hardened construct reasoning with typed objective task intent, separate
  intrinsic/applicability/effective severity fields, typed multi-family repeat
  provenance, and SHA-256 input/rule fingerprints; GUI and shell readers now
  report stale snapshots and refuse outdated inspection actions until rebuild.
- Added a read-only `gentle.gene_isoform_evidence.v1` feature-expert target,
  deterministic SVG, PATZ1 minus-strand fixture, and Splicing Expert Evidence
  tab that keep transcript geometry, RNA/cDNA support, probe constraints,
  expression, and existing qPCR candidates as separately auditable layers.
- Extended that inspector with selected projected BED/BigWig occupancy lanes,
  shared-scale transcript-aligned SVG rendering, GUI/CLI selection, and a local
  TA/DN CUT&RUN figure runbook without treating locus occupancy as isoform
  regulation.
- Added the pure-read `gentle.gene_locus_evidence_display.v1` composition and
  deterministic SVG target, aligning strand-aware transcript/CDS metrics,
  annotated start/stop glyphs, explicitly grouped occupancy lanes, continuous
  JASPAR motif scores, and deduplicated junction-qPCR markers; added a PATZ1
  layout separating two cell lines and reproducible SVG/PNG/PDF commands.
- Prepared the `v0.1.0-internal.10` release boundary for changes since
  `v0.1.0-internal.9` (2026-06-05), with the genome-anchored TP73 evidence
  viewer as the primary release story and a self-contained summary in
  [`release_notes_v0.1.0-internal.10.md`](release_notes/release_notes_v0.1.0-internal.10.md).
- Consolidated this cut's other major outcomes around Clariom/probe-region and
  RNA evidence, fact-aware agent introspection, gene-set/promoter reasoning,
  review-gated material/service handoffs, and extensive GUI responsiveness and
  background-execution hardening.
- Aligned the interim release gate and release-process documentation, retained
  deprecated ClawBio normalizer modes through `.10`, added indexed release-note
  navigation and a version-consistency guard, and improved package metadata.
- Repaired Codex Local structured responses by making the bridge schema strict-
  compatible, carrying human-readable command preconditions/outcomes, and
  distinguishing schema rejection from genuine network failure.
- Separated Agent Assistant conversation from provider setup in
  `Configuration -> Agent Systems`, and added a bounded, project-stored
  `gentle.agent_conversation.v1` transcript so ephemeral Codex/local-model
  requests retain explicit follow-up context without persisting credentials.
- Added protocol-owned `/history`, `/undo`, and `/redo` aliases. GUI history
  aliases use guarded transitions with full view/cache refresh, while
  agent-proposed undo/redo always requires explicit user confirmation.
- Reordered Agent Systems configuration around the persistent provider
  dropdown, collapsed secondary quick-start/catalog controls, and added
  session-only native-provider token-file fallback with explicit source and
  file-permission diagnostics.
- Replaced the Agent Assistant's horizontally scrolling suggestion table with
  wrapping cards that keep each `Run` control and validation reason visible;
  invalid, empty, and chat-only suggestions remain non-executable.

## 2026-07-18

- Made RNA target-rescue junction and exon-intron-boundary catalogs
  boundary-specific, separated RNA-anchored from unanchored retained-intron
  evidence, and added lockstep paired-end fragment screening with explicit
  physical-read and evidence-unit counts.

## 2026-07-12

- Added the same content inset to the DNA viewer's feature-description pane as
  the feature tree above it, keeping long sequence descriptions clear of the
  hosted-window edge.
- Moved projected microarray intervals from the detached regulatory ceiling
  into the coordinate-aligned feature stack in both the live DNA viewer and
  linear SVG export.
- Extended Codex Local discovery to the current macOS ChatGPT bundled CLI and
  app executable while retaining explicit `CODEX_BIN`, `PATH`, and legacy
  Codex app support.
- Added immediate Codex Local model selection backed by the CLI's visible local
  model metadata, with shared GUI/CLI discovery and `codex --model` forwarding.

## 2026-07-13

- Hardened the offline RNA allele hash screen with strict transcript-reference
  validation, phase-block and unphased allele-level reporting, bounded inline
  read calls with complete streamed TSV output, paired-fragment accounting,
  and local PASS biallelic VCF projection through an explicit transcript map.
- Extended the standalone target rescue screen with optional local exon,
  junction, exon-intron-boundary, intron, and genomic-region k-mer catalogs;
  matches are reported as conservative structural candidates rather than
  transcript, fusion, or novelty calls.

## 2026-07-11

- Reworked the repository landing documentation into a concise, visually
  structured README; moved the detailed cloning, genome-context, primer/qPCR,
  physical-carrier, and agent narratives into `docs/showcase.md`; and
  centralized optional-tool, RNA rescue-screen, figure-provenance, licensing,
  contributor-agent guidance, and the inner/outer agent distinction in their
  dedicated manuals and landing-page overview.
- Closed GUI performance and state-tracking gaps: read-only RNA/Agent workers
  exclude inherited engine history, metadata-only track/planning edits mark the
  project dirty and invalidate stale redo, genome-track auto-sync avoids
  unchanged-frame serialization, idle RNA mapping uses the slower repaint
  cadence, BigWig preflight runs in the background, and cached splicing
  transitions now scale with supported transitions rather than exon pairs.

## 2026-07-10

- Hardened GUI responsiveness: long engine jobs now execute on detached,
  revision-checked snapshots; dirty checks use a constant-time state revision;
  layer/GC/lineage/title/JASPAR derivations are cached; and dense JASPAR and
  feature-tree lists use virtualized rows.
- Made unchanged sequence-window display synchronization a constant-time cache
  hit, while display, structure, feature, topology, sequence-replacement, and
  presentation changes still force a complete settings refresh.
- Moved GenBank, UniProt, Ensembl-protein, UniProt-linked GenBank, remote
  JASPAR, all-anchored track import, tracked-file reapply, and track auto-sync
  work off the egui thread with typed completion, stale-result protection,
  duplicate suppression, and explicit cancel/stop-waiting states.

## 2026-07-09

- Persisted the RNA-read disjoint exonic-part partition and added
  `rna-reads export-dexseq-annotation-gff` plus
  `rna-reads export-dexseq-counts-tsv`, producing a mutually joinable DEXSeq
  flattened annotation and strict two-column HTSeq-style count table without
  changing existing exon-path or abundance exports.

## 2026-07-08

- Made undo/redo checkpoints for display-only DNA-view operations cheap.
  `SetDisplayVisibility` and `SetLinearViewport` now snapshot only
  `DisplaySettings` plus history bookkeeping instead of cloning full sequence
  state, while full checkpoints remain the default for data-mutating
  operations.

## 2026-07-05

- Added a deterministic allele-aware RNA-read hash screen
  (`gentle.rna_allele_hash_screen.v1`) with a standalone
  `gentle_cli allele-hash-screen` adapter, shared
  `rna-reads allele-hash-screen` shell route, and synthetic FUS fixture.

## 2026-07-04

- Removed the in-process Rust `sha1` dependency. Local regenerated
  fingerprints now use the shared SHA-256 helper, while legacy SHA-1 download
  verification remains delegated to external platform tools.

## 2026-07-03

- Updated the direct macOS Objective-C bridge dependencies to the current
  `objc2` 0.6 / `objc2-app-kit` 0.3 / `objc2-foundation` 0.3 stack, migrated
  the native menu bridge to the new class-definition API, and kept the lockfile
  refresh for `num-bigint` and `rustc-hash`.
- Delegated legacy SHA-1 download verification to external platform tools,
  exposed `legacy_sha1` in resource status, and documented
  `GENTLE_SHA1_TOOL` / `GENTLE_DISABLE_LEGACY_SHA1` for optional verification
  control.

## 2026-07-02

- Added an in-memory runtime activity stack exposed through
  `introspect runtime`, MCP `runtime_status`, and Unix SIGUSR1 STDERR dumps.
  Shared-shell dispatch, Agent Assistant requests, genome preparation, and GUI
  background jobs now mirror current phase/progress without writing a status
  file.
- Extended `gentle.runtime_status.v1` with observed `activities[]` sourced
  from existing genome-prepare, CUT&RUN shared-asset, and BLAST async ledgers,
  with explicit process-local/cross-process/stale labelling.
- Added a compact GUI status-bar runtime indicator while live process-local
  runtime frames are active.
- Split genome-prepare transcript indexing progress from gene indexing. Tabular
  annotation prepares now report `index_transcripts` /
  `reuse_transcript_index` with the `Transcript Index` step, while gene-index
  progress remains `index_genes` / `reuse_gene_index`. Existing-manifest
  prepares also report a phase-only `verify_checksums` heartbeat before source
  checksum verification.

## 2026-06-28

- Projected built-in shell help from the shared protocol capability registry.
  Glossary command capability descriptors now carry usage, interfaces, and
  aliases, so `help --format json` and topic help expose one descriptor-backed
  command surface while `docs/glossary.json` remains the transitional
  compile-time seed.
- Extended `introspect capabilities` registry metadata so glossary-backed rows
  expose the same descriptor-backed `usage`, `interfaces`, and `aliases` fields
  used by built-in help.
- Promoted `docs/introspection.md` from proposal wording to an implemented
  contract summary, including the `introspect verify-effects` route and the
  actual headless `ui.host_available` fact shape.
- Moved the full glossary/help-generation inversion out of the implemented
  introspection contract and into the roadmap parking lot as deferred work.
- Added a first glossary-to-protocol projection bridge for JSON help output:
  `help --format json` and topic help rows now include the matching shared
  protocol `capability` descriptor for glossary commands when available.
- Completed fact-aware introspection coverage for the current shared
  capability registry. The last pool/container, primer-specificity, and
  reserved screenshot rows now have descriptors. List-valued pool/container
  rows now use `foreach_arg` readiness atoms so every supplied sequence or
  container id is checked before execution when arguments are bound.
- Added fact-aware introspection for CUT&RUN and RNA-read interpretation
  routes, including dataset preparation/projection, read interpretation,
  regulatory-support reports, gene-set regulatory support, and RNA-read
  batch-map outputs.
- Added fact-aware introspection for external resource/report routes, including
  array probe-region conversion/planning, ortholog promoter reports, miRNA
  target scans, RNA transcript-index export, service project-source routing,
  and batch plan/run commands.
- Added fact-aware introspection for explicit cloning/handoff rows, including
  Gibson preview/apply payload routes, protein-to-DNA construct-reasoning
  handoffs, reporter construct handoffs, and the raw restriction-cloning PCR
  handoff operation.
- Added fact-aware introspection for gene-group/gene-set resource rows,
  including draft gene-group reports, gene-set resolution/producer routes, and
  promoter-cohort construction.
- Added fact-aware introspection for adapter/file utility rows, including
  run-bundle shell export, pool import helpers, RNA-read sample-sheet shell
  export, GenBank writer helpers, and VCF display-filter configuration helpers.
- Added fact-aware introspection for repeat/similarity evidence routes:
  repeat annotation queries, repeat overlap/materialization reports,
  repeat-window TFBS summaries, and raw TFBS hit scanning.
- Added fact-aware introspection for report/export/resource/promoter summary
  operations, including process/lab exports, JASPAR entry inspection, RNA-read
  sample-sheet export, multi-gene promoter TFBS summary/render, and scalar
  promoter evidence summaries.
- Added fact-aware introspection for shell-level prepared-genome UI query
  routes (`ui prepared-genomes`, `ui latest-prepared`) and the lower-case
  `set_parameter` adapter alias.
- Added fact-aware introspection for reference/helper genome extraction and
  anchor extension/verification routes. Extraction effects verify loaded
  `sequence.exists(OUTPUT_ID)` results while prepared-cache validation remains
  execution-time behavior.
- Added fact-aware introspection for raw single-container transform operation
  rows: `DigestContainer`, `LigationContainer`, and
  `FilterContainerByMolecularWeight`. These rows require
  `container.exists(CONTAINER_ID)` and model product creation as
  `may_on_success` because product ids are execution-derived.
- Added fact-aware introspection for RNA-read report inspection shell routes:
  `rna-reads show-alignment`, `rna-reads show-alignments`,
  `rna-reads summarize-gene-support`, `rna-reads inspect-gene-support`,
  `rna-reads inspect-alignments`, and `rna-reads inspect-concatemers`.
- Added closed-world `isoform_panel.exists` and `isoform_panel.seq_id`
  introspection facts for imported curated isoform panels. Panel import,
  inspect, render, validation, and raw render operation rows now expose
  fact-aware readiness/effect descriptors that preserve the sequence binding.
- Added closed-world `exon_skip_plan.exists` introspection for persisted
  exon-skip selection plans. Transcript derivation aliases, splicing-reference
  derivation, residue-coordinate lookup, exon-skip planning/materialization,
  and raw exon-skip operation rows now expose fact-aware readiness/effect
  descriptors.
- Added fact-aware introspection for the glossary `digest` alias and raw
  `SelectCandidate` operation row. Both require
  `sequence.exists(INPUT_SEQ_ID)`; `SelectCandidate` can verify a deterministic
  `sequence.exists(OUTPUT_ID)` result.
- Added fact-aware introspection for raw `ExtractAnchoredRegion` and the
  `render-dotplot-svg` shell alias. Anchored extraction is readiness-only
  because candidate ids are prefix/rank-derived; the dotplot alias requires the
  loaded sequence and persisted dotplot facts and models SVG output as an
  external handoff.
- Added fact-aware introspection for raw `ProjectGenomeInterval`,
  `TestCdnaPcr`, `TestCdnaQpcr`, and `TestCdnaQpcrFasta`. Genome interval
  projection and FASTA screening remain execution-time external-file reads,
  while transcript-derived cDNA assay tests require `sequence.exists(SEQ_ID)`
  and declare conservative optional-output effects.
- Added matching fact-aware introspection for the shell cDNA assay-test routes:
  `primers test-cdna-pcr`, `primers test-cdna-qpcr`, and
  `primers test-cdna-qpcr-fasta`.

## 2026-06-27

- Added fact-aware introspection for sequence-gated variant promoter/reporter
  routes and raw operation rows: `variant annotate-promoters`,
  `AnnotatePromoterWindows`, `variant promoter-context`,
  `SummarizeVariantPromoterContext`, `variant reporter-fragments`,
  `SuggestPromoterReporterFragments`, and `AnnotateTfbs`. Optional JSON output
  paths are modeled as `artifact.written` external handoffs rather than
  invented persistent project reports.
- Added closed-world `dotplot.exists` and `flexibility_track.exists` project
  facts plus fact-aware introspection for dotplot/flex compute, show, overlay,
  and SVG render routes. Deterministic `DOTPLOT_ID`/`TRACK_ID` bindings can now
  be verified after compute operations, while dotplot SVG export remains an
  external handoff.
- Added closed-world `candidate_set.exists` introspection for persisted
  candidate-window sets and fact-aware readiness/effect descriptors for
  candidate generation, show/metrics, score, filter, top-k, Pareto, set algebra,
  and delete routes. Delete rows are readiness-gated and verify closed-world
  absence after successful deletion.
- Added closed-world `guide_set.exists`, `guide_filter_report.exists`, and
  `guide_oligo_set.exists` introspection for guide-design workflows, including
  fact-aware descriptors for guide-set upsert/show/filter/delete, oligo
  generation/list/show, and guide oligo/protocol exports.
- Added closed-world `container.exists`, `arrangement.exists`, and
  `rack.exists` introspection for persisted wet-lab container/rack authoring,
  including readiness/effect descriptors for arrangement creation, rack
  placement/mutation, and rack SVG/OpenSCAD/simulation exports.
- Added closed-world `workflow_macro_template.exists`,
  `candidate_macro_template.exists`, and `macro_instance.exists`
  introspection for persisted macro templates and recorded macro lineage rows,
  including readiness/effect descriptors for template show/upsert/delete/run
  and macro-instance inspection.
- Added fact-aware external sequence creation introspection for `LoadFile`,
  `genbank fetch`, `FetchGenBankAccession`, `ensembl-region fetch`,
  `FetchEnsemblRegion`, `dbsnp fetch`, `FetchDbSnpRegion`,
  `FetchUniprotLinkedGenBank`, `ImportUniprotEntrySequence`, and Ensembl
  gene/protein entry sequence import routes. These routes have no project-state
  preconditions, and post-run verification checks `sequence.exists(OUTPUT_ID)`
  when a deterministic id was supplied.
- Added fact-aware introspection for tracked genome signal-file subscription
  routes: `tracks tracked list`, `tracks tracked add`,
  `tracks tracked remove`, `tracks tracked clear`, and
  `tracks tracked apply`.
- Added fact-aware introspection for planning consultation, protein-expression
  handoff, planning profile/objective mutation, suggestion resolution, and
  sync status/pull/push routes.
- Added fact-aware introspection for resource sync/import/install/benchmark
  routes plus local cache cleanup and publication-dataset preparation routes.
- Added fact-aware introspection for raw core sequence operation rows:
  `SaveFile`, `Digest`, `Pcr`, `PcrAdvanced`, `PcrMutagenesis`, and
  `PcrOverlapExtensionMutagenesis`. PCR rows require a loaded template and can
  verify deterministic product ids; digest and overlap-extension rows are
  readiness-only until prefix/rank-derived products are projected as facts.
- Added closed-world `uniprot_entry.exists`,
  `ensembl_gene_entry.exists`, and `ensembl_protein_entry.exists`
  introspection for stored protein/gene metadata entries. UniProt and Ensembl
  fetch/import routes can verify explicit entry ids, show/import-sequence
  routes require stored metadata facts, and metadata-backed sequence imports can
  still verify deterministic `sequence.exists(OUTPUT_ID)` products.
- Added closed-world `uniprot_projection.exists` introspection for persisted
  UniProt genome projections. Projection-specific show, feature-coding,
  Ensembl-link resolution, transcript-accounting, Ensembl-comparison, and
  audit-generation routes now use that fact for readiness; projection
  generation can verify a deterministic `PROJECTION_ID`.
- Tightened construct-reasoning repeat-family wording so rmsk-backed Alu-like
  evidence uses curated-family support language while rmsk-absent heuristic
  calls retain the soft-catalog caveat.
- Extended fact/introspection projection with report, view, and configuration
  facts plus glossary/parity metadata for `facts` and `introspect` shell routes.
- Added Agent Assistant hover help for prompt templates, project-state context,
  and auto-run suggestions, and fixed External Services window routing so
  native child viewports render directly instead of appearing as nested hosted
  windows.

## 2026-06-26

- Added the first shared `introspect facts|capabilities|readiness|all`
  implementation, including additive fact domains, explicit headless
  `ui.host_available=false`, deterministic agent-system
  `host.tool_available` facts, registry-backed capability discovery,
  fact-annotated readiness descriptors, and `--seq-id` / `--readiness`
  introspection scoping.
- Added `introspect verify-effects` to verify fact-annotated
  `must_on_success` effects against the current fact graph plus supplied
  evidence, starting with restriction-scan report verification.
- Added fact-aware introspection for `sequence create`, including readiness
  with no fact preconditions and post-run verification that the requested
  output sequence id exists.
- Added fact-aware introspection for `variant materialize-allele`, including
  input sequence readiness and post-run verification of the materialized output
  sequence id.
- Added fact-aware readiness introspection for `align compute`, checking that
  both loaded query and target sequence ids exist before alignment.
- Added matching fact-aware introspection for the raw `AlignSequences`
  operation row so registry-driven adapters can make the same query/target
  readiness check as the shell command.
- Added fact-aware readiness introspection for `render-svg`, `render-rna-svg`,
  and `render-lineage-svg`, checking sequence inputs where applicable and
  modeling external SVG files as `artifact.written` `external_handoff` effects.
- Added matching fact-aware introspection for the raw render operation rows
  `RenderSequenceSvg`, `RenderRnaStructureSvg`, `RenderFeatureExpertSvg`,
  `RenderTfbsScoreTrackCorrelationSvg`, and `RenderLineageSvg`, keeping
  registry-driven adapters aligned with the shell render commands.
- Added fact-aware introspection for unambiguous raw persisted-report operation
  rows: `ListSequencingConfirmationReports`, `ListCutRunReadReports`,
  `ListRnaReadReports`, `ExportSequencingConfirmationReport`,
  `ExportSequencingConfirmationSupportTsv`, `ExportCutRunReadCoverage`, and
  `ExportRnaReadReport`.
- Added fact-aware introspection for raw persisted-report show operations:
  `ShowSequencingConfirmationReport`, `ShowCutRunReadReport`, and
  `ShowRnaReadReport`.
- Updated `introspect readiness` to evaluate full descriptor
  `precondition_expr` trees, including `any` branches, and added fact-aware
  introspection for the shared raw `ExportPrimerDesignReport` operation across
  primer and qPCR design reports.
- Added fact-aware introspection for raw primer/qPCR design operations:
  `DesignPrimerPairs`, `DesignInsertionPrimerPairs`, and `DesignQpcrAssays`.
- Added fact-aware introspection for raw `ReverseTranslateProteinSequence`,
  including protein-kind input readiness and deterministic output sequence
  effect verification.
- Added fact-aware introspection for raw `MaterializeVariantAllele`, including
  input-sequence readiness and deterministic output sequence effect
  verification.
- Added fact-aware readiness introspection for the read-only `rna-info`
  sequence inspector.
- Added fact-aware introspection for `reverse-translate run`, including
  protein-kind input readiness and post-run verification of the requested
  coding-DNA output sequence id.
- Added `report.exists` fact projection for persisted reverse-translation
  reports and fact-aware readiness introspection for
  `reverse-translate show-report`.
- Added `report.exists` fact projection for persisted primer/qPCR design
  reports and fact-aware introspection for `primers design` and
  `primers design-qpcr`, including template-sequence readiness and post-run
  verification of the requested report id.
- Added `report.exists` fact projection and fact-aware introspection for
  `primers prepare-restriction-cloning`, checking template, primer-report, and
  destination-vector readiness and verifying the generated handoff report.
- Added fact-aware readiness introspection for persisted primer, qPCR, and
  restriction-cloning report inspection commands:
  `primers show-report`, `primers show-qpcr-report`, and
  `primers show-restriction-cloning-handoff`.
- Added fact-aware, no-precondition catalog readiness for report-list commands:
  `primers list-reports`, `primers list-qpcr-reports`,
  `primers list-restriction-cloning-handoffs`, and
  `reverse-translate list-reports`.
- Added fact-aware readiness introspection for `features tfbs-summary`, treating
  it as a read-only sequence inspector whose readiness depends on the inspected
  `sequence.exists` fact.
- Added matching fact-aware introspection for the raw `SummarizeTfbsRegion`
  operation row with the same sequence-readiness model.
- Added fact-aware introspection for the raw
  `QueryProteinResidueGenomicCoordinates` operation row as a read-only
  sequence inspector.
- Added fact-aware catalog readiness for self-description commands
  `capabilities` and `state-summary`, both with empty project preconditions.
- Added fact-aware readiness introspection for `features query` and
  `features export-bed`; feature queries are read-only sequence inspectors,
  while BED exports model the output path as an `artifact.written` external
  handoff.
- Added matching fact-aware introspection for the raw `ExportFeaturesBed`
  operation row, keeping registry-driven adapters aligned with the shell BED
  export command.
- Added fact-aware introspection for raw sequence-context operations:
  `InspectSequenceContextView` checks sequence readiness, and
  `ExportSequenceContextBundle` additionally models its output directory as an
  external artifact handoff.
- Added fact-aware readiness introspection for `inspect-feature-expert` and
  `render-feature-expert-svg`; both require the inspected sequence to exist,
  and the SVG renderer models its output path as an `artifact.written`
  external handoff.
- Added matching fact-aware introspection for the MCP `restriction_site_detail`
  tool row, preserving the same sequence-readiness semantics as its shared
  `inspect-feature-expert ... restriction` shell route.
- Added fact-aware introspection for genome preparation/cache maintenance and
  BLAST job-control routes (`genomes|helpers prepare`, install/remove,
  synchronous BLAST, async BLAST start/cancel, `PrepareGenome`,
  `prepare_genome`, and JS/Lua BLAST helpers). These rows intentionally model
  catalog/cache/job effects as external local state rather than project facts.
- Added fact-aware introspection for generic execution routes
  (`apply_operation`, `apply_workflow`, MCP `op`/`workflow`, `macros run`, and
  `candidates macro`) using command-dependent `may_on_success` effects instead
  of pretending their concrete project effects are knowable before parsing the
  supplied payload.
- Registered `view.selection` and `view.visible_tracks` fact vocabulary entries
  and added fact-aware catalog readiness for the `display` visibility command
  with a `view.visible_tracks` view-session effect.
- Added fact-aware catalog readiness for `ui intents`, exposing GUI intent
  discovery as a no-precondition view-intent route.
- Added fact-aware catalog readiness for `help`, exposing shell help/topic
  lookup as a no-precondition command-catalog route.
- Added fact-aware catalog readiness for `history status`, exposing
  undo/redo availability as a no-precondition status route.
- Registered `facts graph` and `facts eval` in the shared glossary and added
  fact-aware catalog readiness for both fact-layer routes.
- Added fact-aware readiness introspection for `reverse-translate export-report`,
  checking the source report and modeling the output JSON path as an
  `artifact.written` external handoff.
- Added fact-aware readiness introspection for `primers export-report` and
  `primers export-qpcr-report`, using persisted report facts and external JSON
  handoff effects.
- Added fact-aware readiness introspection for
  `primers export-restriction-cloning-handoff`, using the persisted handoff
  report fact and an external JSON handoff effect.
- Added `report.exists` fact projection for persisted sequencing-confirmation
  reports and fact-aware readiness introspection for `seq-confirm list-reports`,
  `seq-confirm show-report`, `seq-confirm export-report`, and
  `seq-confirm export-support-tsv`.
- Added fact-aware, no-precondition catalog readiness for
  `cutrun list-read-reports` and `rna-reads list-reports`.
- Extended `config.param` fact projection to cover persisted display settings,
  metadata-backed BLAST option parameters, and accepted `set-param` aliases, so
  Agent/MCP readiness and effect checks can reason about display configuration
  with the same closed-world facts as engine parameters.
- Added closed-world absence effect verification for candidate-set, guide-set,
  workflow macro-template, and candidate macro-template delete routes using
  `not(...exists)` hard effects.
- Added `report.exists` fact projection for persisted CUT&RUN read and RNA-read
  interpretation reports, plus fact-aware readiness introspection for
  `cutrun show-read-report`, `cutrun export-coverage`,
  `rna-reads show-report`, and `rna-reads export-report`.
- Added fact-aware introspection for raw RNA-read gene-support operations:
  `SummarizeRnaReadGeneSupport` and `InspectRnaReadGeneSupport` require an
  existing `rna_read` report and treat optional JSON output paths as external
  artifact handoffs.
- Added fact-aware introspection for no-project local catalog/report operations
  `SummarizeJasparEntries`, `BenchmarkJasparRegistry`, `ListJasparCatalog`,
  `ResolveTfQueries`, `ListReporterCatalog`, and `RecommendReporters`, treating
  optional JSON output paths as external artifact handoffs.
- Added `config.param` fact projection for engine-owned configuration
  parameters and fact-aware `set-param` introspection, including JSON-value
  effect verification through `introspect verify-effects`.
- Added fact-aware introspection for direct sequence-derivation operation rows
  `Reverse`, `Complement`, `ReverseComplement`, `Branch`, and `ExtractRegion`,
  checking the input sequence and verifying deterministic output sequence ids.
- Added fact-aware introspection for the raw `SetTopology` operation row,
  checking the target sequence and verifying the projected `sequence.circular`
  boolean fact.
- Added readiness-only fact-aware introspection for the raw `RecomputeFeatures`
  operation row, checking that the target sequence exists without inventing a
  computed-feature freshness fact.
- Added the closed-world `view.viewport` fact projection and fact-aware
  introspection for the raw `SetLinearViewport` operation row, including nested
  argument binding for exact viewport verification.
- Added closed-world `view.visible_tracks` fact projection from persisted
  display-layer booleans so words-only clients can read back visible track
  state, not only see display intent descriptors.
- Added fact-aware introspection for raw `SetDisplayVisibility` and
  `SetParameter` operation rows, mirroring the existing `display` and
  `set-param` shell descriptors without introducing duplicate fact names.
- Added fact-aware introspection for specialized RNA-read artifact exports
  such as hit FASTA, target-quality, exon-path, exon-abundance, score-density,
  alignment TSV, isoform-triage, and alignment-dotplot exports, including raw
  engine operation rows.
- Added fact-aware introspection for RNA-read alignment, isoform preflight, and
  hit-sequence materialization routes, with conservative effect declarations
  for selection-dependent materialized sequence ids.
- Added fact-aware introspection for built-in DNA/RNA ladder catalog list and
  export routes, including JS/Lua helper names and raw ladder export operation
  rows.
- Added fact-aware introspection for agent-system catalog list routes
  `agents list`, `agent_systems`, and `list_agent_systems`, keeping configured
  system enumeration distinct from live adapter/model-discovery readiness.
- Added fact-aware introspection for the raw `agent_preflight` MCP/tool
  operation row, mirroring `agents preflight` through the shared
  `host.tool_available(SYSTEM_ID)` readiness fact.
- Added fact-aware introspection for protocol-cartoon catalog, render, template
  validation, and template export routes, including the raw engine operation
  rows and external SVG/JSON artifact handoff effects.
- Added fact-aware introspection for no-project external inspection routes:
  prepared-cache inspection, CUT&RUN dataset catalog/status inspection, and
  array helper inspection/rendering, with SVG artifact handoff modeling where
  applicable.
- Added fact-aware readiness descriptors for genome-track imports, BLAST-track
  imports, and array projection routes. These rows require the loaded target
  sequence and intentionally remain readiness-only until feature
  freshness/track-update facts are projected.
- Added fact-aware introspection for no-project catalog/list routes covering
  candidate sets, candidate macro templates, guide sets, workflow macro
  instances/templates, and routine catalog list/explain/compare operations.
- Added fact-aware introspection for construct-reasoning graph list routes,
  treating optional sequence ids as filters rather than readiness
  preconditions.
- Added closed-world `construct_reasoning_graph.exists` introspection for
  persisted construct-reasoning graphs, including readiness/effect descriptors
  for named graph inspection, inspection-action listing/running, annotation
  status/writeback routes, and graph JSON export.
- Added fact-aware introspection for persisted dotplot and flexibility-track
  list routes, treating optional sequence ids as filters rather than readiness
  preconditions.
- Added fact-aware introspection for local metadata/catalog routes covering
  reference/helper catalog listing, stored Ensembl gene/protein metadata lists,
  and gene-group list/show/resolve/doctor routes with optional JSON artifact
  handoffs.
- Added fact-aware introspection for shell-level resource/catalog inspection
  routes covering JASPAR summary/list/inspect/TF-query resolution, resource
  status, UCSC rmsk indexing suggestions, publication dataset list/status,
  reference/helper catalog validation, and helper vocabulary list/doctor
  commands.
- Added fact-aware introspection for reporter catalog shell routes
  `reporters list`, `reporters recommend`, `reporters export-corpus`, and raw
  `ExportReporterCorpus`, keeping shell and operation-level catalog readiness
  aligned.
- Added fact-aware introspection for service readiness/provider catalog routes
  `services status`, `services providers list`, and
  `services providers doctor`, with optional doctor JSON output modeled as an
  external artifact handoff.
- Added fact-aware introspection for planning read-back routes
  `planning profile show`, `planning objective show`, and
  `planning suggestions list`, keeping mutating planning updates separate.
- Added fact-aware introspection for host/helper/protease/microRNA catalog
  helper routes, including JS/Lua/MCP helper catalog row names and optional
  protease JSON artifact handoffs.
- Added fact-aware introspection for reference/helper genome cache inspection
  and adapter readbacks (`genomes status|genes`, `helpers status|genes`,
  `list_reference_genomes`, `list_reference_catalog_entries`,
  `is_reference_genome_prepared`, `list_reference_genome_genes`), keeping
  concrete catalog/cache validation as execution-time behavior.
- Added fact-aware descriptors for the top-level read-only introspection shell
  routes: `introspect facts`, `introspect capabilities`,
  `introspect readiness`, `introspect verify-effects`, and `introspect all`.
- Added fact-aware introspection for async BLAST status/list readbacks:
  `genomes blast-status`, `helpers blast-status`, `genomes blast-list`,
  `helpers blast-list`, `blast_async_status`, and `blast_async_list`, while
  leaving start/cancel/execution rows registry-only until their external
  runtime effects are modeled.
- Added fact-aware introspection for agent model-discovery routes
  `agents discover-models` and `agent_models`, using the same
  `host.tool_available(SYSTEM_ID)` readiness model as agent preflight while
  declaring no project effects for model enumeration.
- Added fact-aware introspection for adapter parity aliases:
  `state_summary`, `reference_catalog_entries`, `ui_intents`,
  `ui_prepared_genomes`, and `ui_latest_prepared`.
- Added fact-aware introspection for generic GUI intent requests:
  `ui open`, `ui focus`, `ui close`, and `ui_intent`, with readiness gated by
  `ui.host_available`.
- Added fact-aware introspection for sequence-scan report/render routes:
  `FindRestrictionSites`, `features tfbs-score-tracks-svg`,
  `RenderTfbsScoreTracksSvg`, `SummarizeTfbsScoreTracks`,
  `features tfbs-track-similarity`, and `SummarizeTfbsTrackSimilarity`.
- Added fact-aware introspection for local external-service handoff routes:
  `services delivery-route`, `services project-preflight`,
  `services project-quote`, `services handoff`, and `services guide`, while
  keeping `services route-project-source` separate until conditional
  project-object preconditions are modeled.
- Added fact-aware introspection for agent dispatch rows: `agents ask`,
  `ask_agent_system`, `agents plan`, and `agent_plan` now use
  `host.tool_available(SYSTEM_ID)` readiness, while `agents execute-plan` and
  `agent_execute_plan` model supplied-plan execution as payload-ready with
  command-dependent `may_on_success` effects.
- Added fact-aware introspection for SRA/read-acquisition job-control rows:
  `reads acquire status`, `reads acquire prepare`, `reads acquire inspect`,
  `reads acquire cancel`, and raw `ReadAcquire*` rows are payload/path-ready,
  with prepare/cancel external state updates modeled as `artifact.written`
  handoffs.
- Added fact-aware introspection for Ensembl discovery/catalog-maintenance
  routes: availability, installable-list, preview, and update-spec rows are
  ready without project state, while update rows model external catalog writes
  as `artifact.written` handoffs.
- Added fact-aware introspection for project/session utility rows:
  `history undo`, `history redo`, `load-project`, `load_project`,
  `save-project`, `save_project`, and `load_dna` now expose conservative
  readiness/effect descriptors.
- Added closed-world `sequencing_trace.exists` introspection for imported
  sequencing-trace evidence records, including readiness/effect descriptors for
  trace import/list/show shell routes and raw operation rows.
- Added fact-aware introspection for sequencing-confirmation execution and
  primer overlay suggestion rows: `seq-confirm run`, `ConfirmConstructReads`,
  `seq-primer suggest`, and `SuggestSequencingPrimers` now model expected
  construct sequence readiness plus deterministic confirmation-report effect
  verification when explicit report ids are supplied.
- Added fact-aware introspection for protease catalog/digest and protein-gel
  rendering routes. Protease digest readiness now requires an existing
  protein-kind sequence, persisted protein-derivation reports project as
  `report.exists == protein_derivation`, and SVG render routes model their
  output paths as `artifact.written` handoffs.
- Added fact-aware introspection for raw transcript/protein/splicing
  derivation rows: `DeriveTranscriptSequences`, `DeriveProteinSequences`, and
  `DeriveSplicingReferences` now require `sequence.exists(SEQ_ID)`, declare
  conservative derived-sequence effects, and verify protein-derivation reports
  when deterministic report ids are supplied.
- Added fact-aware introspection for non-mutating primer helper readbacks:
  Primer3/backend preflight, feature/splicing ROI seed helpers, restriction
  cloning vector suggestions, and restriction-cloning handoff request seeding.

## 2026-06-24

- Added shared Agent Assistant/shell controls for DNA sequence-window selection
  (`ui selection sequence-window ...`) and project display-layer visibility
  (`display show|hide|visibility`), keeping window/selection actions separate
  from sequence deletion.
- Added Promoter design GUI parity for offline ortholog promoter cohorts and
  comparisons, plus ortholog relationship expectations that emit non-blocking
  TFBS/CUT&RUN divergence or concordance flags.

## 2026-06-22

- Added retrieval-producer metadata to `gentle.gene_set_resolution.v1` and the
  offline `gene-sets produce direct-list` route for local JSON/TSV gene-list
  caches, resolving candidates through the existing explicit-member resolver.
- Added the offline `gene-sets produce ontology-assignment` route for local
  ontology assignment JSON/TSV caches, keeping provider term membership distinct
  from local gene-group `external_mapping` resolution.
- Added the offline `gene-sets produce co-regulated` route for local
  evidence-derived cohort caches with explicit dataset/contrast, score
  threshold, sign-direction, and relationship-expectation metadata.
- Made produced gene-set resolutions, promoter cohorts, and gene-set CUT&RUN
  support reports persist as logical lineage artifacts rendered as `GeneSet`
  nodes, linked from producer operation to downstream analysis.

## 2026-06-21

- Added shared operand metavariable conventions for glossary `usage` rows and
  updated Agent Assistant prompt guidance so inner helpers consult GENtle docs
  before proposing commands and ask instead of guessing ambiguous operands.
- Added JS and Lua adapter bindings for construct-reasoning inspection actions,
  so scripts can list and run the same portable dotplot recommendations exposed
  through CLI, ClawBio, and MCP shell routes.
- Made construct-reasoning repeat/similarity task severity quantitative and
  objective-specific with protocol `score` fields, explicit score-to-bucket
  thresholds, and visible objective boost/down-weight rationales.

## 2026-06-19

- Deepened construct-reasoning repeat-family interpretation with a shared
  Alu/SINE/LINE/LTR/satellite taxonomy, stricter internal-vs-curated
  corroboration, and confidence tiers for repeat-family provenance.
- Added rule-based `task_severities[]` to construct-reasoning facts so
  repeat/similarity warnings can report PCR, nanopore, read-mapping, cloning
  stability, and construct-maintenance severity without creating extra map
  overlays.
- Completed the Clariom/probe-region real-data-to-figure bridge slice with
  coordinate-consistent `gentle.probe_region_evidence_interpretation.v2`,
  gated `arrays run-probe-region-backend`, GUI shared-capability surfacing, and
  runbook/docs coverage for the explicit E-MTAB-14704 TP73 loop.

## 2026-06-18

- Integrated materialized RepeatMasker/UCSC `rmsk`-style repeat annotations
  into construct reasoning so overlapping curated Alu/SINE repeat-family rows
  back soft internal repeat/mobile-element calls without duplicating fact rows.
- Made native HTTP agent transports tolerant of local models that omit or
  mis-shape the `schema` field when the returned JSON otherwise matches
  `gentle.agent_response.v1`; external stdio adapters remain strict.
- Added an explicit Msty MLX OpenAI-compatible agent template for
  `http://localhost:11973/v1`, with GUI/CLI docs that distinguish it from the
  `11964` Msty gateway when that gateway reports no model ids.
- Added construct-reasoning inspection-action rationale to the portable action
  payload and surfaced action mode, focus, evidence ids, and rationale in the
  existing GUI inspector rows.
- Added offline-first ortholog promoter reasoning with local
  `gentle.ortholog_resource.v1` mapping resources,
  `ResolveOrthologPromoterCohort`, `SummarizeOrthologPromoterComparison`, and
  `orthologs resolve-promoter-cohort` / `orthologs promoter-comparison` shell
  routes, including species/genome-matched CUT&RUN motif-support states.
- Added additive `relationship_flags[]` to
  `gentle.promoter_cohort_comparison.v1`, surfacing unexpected TFBS-track
  divergence for declared co-regulated promoter cohorts and unexpected
  concordance for declared anti-co-regulated cohorts as non-blocking review
  cues.

## 2026-06-17

- Documented sequence collection subjects as the shared model for gene sets,
  pools, arrangements, alignments, derived collections, and storage projections,
  with a GUI implementation plan for adapter-parity operation lifting.
- Completed gene-set cohort relationship support by deriving evaluated-only
  CUT&RUN occupancy consistency flags for declared co-regulated and
  anti-co-regulated promoter cohorts.
- Extended `planning protein-expression-handoff --seq-id` with read-only
  sequence/readiness/CDS/tag context so high-yield protein-expression handoffs
  report product uncertainty before review-gated provider or construct work.
- Added optional gene-set promoter-cohort relationship expectations
  (`manual`, `co_regulated`, `anti_co_regulated`) to promoter cohort and
  CUT&RUN aggregate reports, with non-blocking expectation flags and
  `--relationship` shell parsing.
- Routed protein-expression `suggested_next_actions[]` from
  `product_definition.readiness.status`, distinguishing CDS candidates,
  protein-only targets, and sequences that need CDS/ORF boundary review.
- Added Promoter design GUI parity for promoter expression evidence and
  cached CUT&RUN regulatory-support reports, additive four-state TFBS support
  status, and promoter cohort comparison; `genomes promoter-cohort-comparison`
  emits `gentle.promoter_cohort_comparison.v1`.
- Added portable construct-reasoning inspection actions for repeat/similarity
  dotplot handoffs, carrying deterministic action ids and driving evidence ids
  in the graph payload instead of deriving GUI buttons from labels.
- Added `construct-reasoning list-inspection-actions` and
  `run-inspection-action` so CLI/ClawBio users can list the same portable
  repeat/similarity dotplot recommendations and compute/render the selected
  action through shared dotplot operations; listing now supports fact,
  annotation/candidate, evidence, sequence, action-kind, and summary filters,
  and ClawBio exposes typed request modes
  `construct-reasoning-list-inspections` and
  `construct-reasoning-run-inspection` over the same shell contract.
- Added optional true PM probe-intensity matrix input to
  `arrays import-apt-probe-region-output` and the Clariom D GUI panel, with
  probe-level sample, condition-summary, and logFC columns preserved in
  `probe_intensity_chrom_order.csv`.
- Added `--level pm_probe` projection for probe-region helper outputs, so true
  PM probe rows marked `probe_level_input` can become genome-anchored array
  features without promoting parent-summary fallback rows.
- Added `arrays interpret-probe-region-evidence` /
  `InterpretProbeRegionEvidence` to compare projected probe/probeset-region
  array features with transcript/exon geometry while preserving shared
  transcript, multi-hit, and coordinate-projection ambiguity.
- Added first-class sequence-window controls to run and preview
  `InterpretProbeRegionEvidence` reports from the Clariom D / probe-region
  evidence panel, including bounded evidence and transcript geometry tables.
- Extended `InterpretProbeRegionEvidence` with explicit per-evidence
  transcript mappings that record exon ordinals, exon ranges, junction spans,
  and overlap base counts without turning array evidence into isoform calls.
- Added conservative geometry scores and score-basis guardrails to
  `InterpretProbeRegionEvidence`, with transcript-level unique/shared/
  constraining score sums for review-only probeset evidence triage.
- Added transcript-level `review_status` labels to probe-region interpretation
  reports so GUI/CLI users can triage unique, shared, constraining, and absent
  geometry without treating the result as an isoform call.

## 2026-06-16

- Added `services route-project-source` plus External Services GUI helpers so
  selected sequences/spans, persisted oligo order forms, and primer report pair
  ranks can produce delivery-route candidates while preserving duplicate-review
  gates before quote handoff.
- Added first-class sequence-window GUI controls for importing explicit APT
  probe-region summaries plus annotation tables into inspectable GENtle helper
  output, reusing the shared shell import/inspect/export/project contracts.
- Added Phase B oligo order handoff routes: `primers oligo-order route` and
  `quote` now reuse external-service delivery/quote contracts, block unreviewed
  duplicate forms before quote handoff, and preserve oligo line-item
  provenance/modification fields in normalized JSON/CSV bundles.

## 2026-06-15

- Added Phase A first-class oligo order forms under `primers oligo-order`,
  with deterministic line ids, primer/qPCR report provenance, duplicate/reuse
  grouping, review marking, list/show/export, and persistence inside the
  existing primer-design report store.
- Added the command-palette `Evidence Preparation` assistant for the TP73
  evidence-viewer proof path, reusing shared operations for local repeat,
  array, BED, TFBS, and proof-export preparation while keeping CEL/R/vendor
  steps as explicit copy-command handoffs.
- Added `services delivery-route` and the
  `gentle.external_service_delivery_route.v1` contract so generic "deliver this
  sequence" wording is classified by sequence kind before selecting Metabion or
  GeneArt quote-handoff routes.
- Added ClawBio external-service intents and request examples for provider
  catalog/doctor checks plus review-only Metabion oligo/m-block and GeneArt
  cloned-gene/protein-expression preflight and quote handoff routes.
- Routed high-yield protein-expression ClawBio intents such as "maximal amount
  of protein" to the read-only `planning protein-expression-handoff` request
  example and added the same scenario to the experimental follow-up catalog,
  graph, and human guide, with a review-gated GeneArt protein-expression quote
  packet as the downstream provider handoff.
- Added an explicit `services project-quote @...` suggested next action to the
  protein-expression handoff report so GeneArt preflight and quote-packet
  preparation remain separate review stages.
- Added a native GUI inspector for completed probe-region helper output, plus
  coordinate/build provenance, bounded preview rows, and projection-readiness
  blockers in `gentle.probe_region_output_inspection.v1`.
- Added `arrays render-probe-region-output-svg` and a matching GUI export
  control for deterministic native SVG plots of inspected `mean_log2_*` and
  `log2FC_*` probe-region helper tracks.
- Added `arrays project-probe-region-output` / `ProjectProbeRegionOutput` for
  direct-coordinate-compatible projection of inspected helper `log2FC_*` rows
  into genome-anchored array features.
- Added first-class sequence-window controls for projecting inspected
  probe-region helper output into genome-anchored array features without
  leaving the Clariom D evidence panel.
- Added explicit `coordinate_projections[]` / `projection_maps[]` support for
  probe-region helper-output projection, reusing the existing interval-map
  contract to preserve native and displayed array coordinates.
- Added explicit Affymetrix Power Tools command planning to
  `arrays probe-regions` when user-supplied PGF/CLF and optional MPS library
  files are present.
- Added `arrays import-apt-probe-region-output` to convert explicit APT summary
  output plus an explicit annotation/NetAffx coordinate table into GENtle's
  probe-region helper-output directory contract.
- Extended explicit APT probe-region imports with optional sample metadata
  matching that writes `mean_log2_*`, `sd_log2_*`, and default `log2FC_*`
  tracks for native inspection, plotting, and projection.
- Added optional `probe_intensity_chrom_order.csv` output for explicit APT
  imports when annotation rows provide PM probe coordinates, marking values as
  parent probeset-summary intensities when no explicit PM probe matrix is
  supplied.

## 2026-06-14

- Added the read-only `planning protein-expression-handoff` route emitting
  `gentle.protein_expression_handoff.v1` with product context, chassis/route
  candidates, high-yield missing questions, and a GeneArt protein-expression
  preflight scaffold.
- Added tutorial `09.03` for using the protein-expression handoff route from
  CLI/agent workflows, plus a ClawBio example request.

## 2026-06-13

- Added `arrays inspect-probe-region-output` to validate completed
  `probe_regions_oligo.R` output directories and summarize region-table,
  sample, condition, logFC, chromosome/gene, manifest, and provenance status as
  `gentle.probe_region_output_inspection.v1`.

## 2026-06-12

- Added `arrays probe-regions` as a shared-shell, GUI-shell-visible
  Affymetrix CEL probe/probeset-region preflight, emitting
  `gentle.probe_region_plan.v1` with CEL, metadata, annotation/library,
  platform, dependency, backend-candidate, metadata/contrast, output, and
  cache-readiness checks.
- Added `scripts/probe_regions_oligo.R` as a generic explicit R/oligo helper
  for the `arrays probe-regions` plan, producing chromosome-ordered intensity
  CSVs, expression/feature TSVs, limma contrast tables, provenance, and an
  RMA-normalized matrix manifest from user-supplied CEL files.

## 2026-06-06

- Added `protein_expression_max_yield` as a planning consultation
  `biological_intent`, with explicit high-yield protein-expression questions
  before any construct/vector route is accepted.

## 2026-06-05

- Switched push/PR CI to one commit-sampled platform per run, with manual
  platform override, a PowerShell-native Windows validation pass matching the
  Unix test shape, clearer `macOS` naming, a single Linux runner when Linux is
  selected, and no multi-GB `target` cache restores.
- Limited the container-image workflow to release tags and manual dispatch so
  ordinary pushes and PRs no longer build Docker CLI/GUI images.
- Added explicit `biological_intent` fields to reporter recommendation and
  reporter construct handoff reports, and staged generated synthetic-biology
  chef background assets for a future synthetic-biology window.

## 2026-06-04

- Wired the ClawBio `gentle-cloning` skill toward the existing reporter
  catalog, recommender, corpus export, and reporter construct handoff routes
  so synthetic-biology follow-up quotes GENtle reports instead of improvising.

## 2026-05-29

- Added opt-in GUI frame profiling behind `--features gui-profiler` and
  `GENTLE_GUI_PROFILE=1`, with coarse Puffin spans for app, DNA-window,
  DNA-map, sequence-text, and feature-tree latency investigations.

## 2026-05-28

- Switched the startup splash screen to the app-ready transparent mascot asset
  at `assets/mascots/Mascot_transparent.png` and enlarged its splash
  presentation.
- Defaulted macOS back to native OS child viewports on `egui/eframe` `0.34.3`,
  with `GENTLE_MACOS_HOSTED_CHILD_VIEWPORTS=1` left as the explicit hosted
  fallback and root-window fullscreen restoration disabled on macOS.
- Delayed new macOS native sequence windows until the root window leaves
  fullscreen/maximized state, avoiding the unstable Split View-style flicker
  path when opening the first child window.
- Tightened hosted egui window frame-drag ownership so resize-edge drags keep
  priority over lower hosted-window body hits while ordinary hosted-window body
  interactions such as DNA selection are no longer treated as frame drags.
- Added a six-well cell-culture plate rack profile/template, culture-well SVG
  rendering, arrangement-labelled top-down `racks hero-svg`, `svg-pdf`
  conversion, and README figure assets for plate-layout presentation.

## 2026-05-27

- Added the `codex_local_stdio` Agent Assistant provider and
  `scripts/codex-agent-bridge`, allowing GENtle to route inner-agent requests
  through a logged-in local Codex CLI/App account without requiring
  `OPENAI_API_KEY`.
- Hardened hosted egui window drag ownership during embedded-window resize/move
  interactions so lower hosted windows do not react to a drag that began on the
  active window frame.
- Switched the egui/eframe family from the temporary upstream git override to
  the published `0.34.3` crates while preserving GENtle's hosted-window
  title-bar movement and lower-window drag lockout in the compatibility wrapper.
- Added the shared `mirna` target-site scan service with built-in
  `hsa-miR-96-5p` seed catalog support, JSON schema
  `gentle.mirna_target_scan.v1`, and CLI/shared-shell commands for seed
  explanation, catalog inspection, and annotated target scanning.
- Added `Patterns -> microRNA Target Scan...` as a graphical GUI wrapper over
  the shared scan command, including seed-pairing drawings, region-specific
  splicing interpretation, and side-by-side ortholog/candidate snippet scans.
- Hardened the microRNA scanner with typed evidence tags, a reverse-strand
  coordinate regression, and warnings when a known catalog name is paired with a
  non-canonical mature-sequence override.

## 2026-05-26

- Refreshed `Cargo.lock` after a pre-release `cargo update` and verified the
  release-facing examples/tutorial checks against the updated patch versions.
- Fixed linear DNA-map drag selection so drags that originate on splitters or
  hosted-window decorations cannot become sequence selections after entering
  the map canvas.
- Added `integrations/clawbio/local_agent_handoff.md` as a shared routing note
  for Codex, Claude, OpenClaw, and other local agents: use GENtle through the
  known ClawBio `gentle-cloning` runner, inspect `result.json`/`report.md` and
  artifacts, and avoid inventing a second GENtle command surface.
- Tightened DNA-map drag ownership checks so hosted auxiliary window drags
  cannot pass through into DNA selection.

## 2026-05-24

- Added `File -> New Sequence...` and `File -> New Sequence from Clipboard...`
  so typed or clipboard IUPAC DNA can become an ordinary project sequence
  through `CreateSequenceFromText`.
- Added a phase-1 GENtle-side ClawBio bridge: compact GUI context export,
  `CancelToken`-based subprocess transport, `result.json` parsing, and an
  updated federation/boundary contract in `docs/clawbio_gentle_integration_onepager.md`.
- Added the phase-2 ClawBio GUI panel under `Services -> ClawBio...`, with
  worker-thread dispatch, cancellation, stdout/stderr output-bundle logs,
  artifact/report links, and verbatim suggested-action request dispatch.
- Fixed the External Services provider/service picker so changing from
  metabion to GeneArt refreshes the editable request JSON and clears stale
  preflight/quote previews.

## 2026-05-23

- Started Phase 4 of the tutorial presentation overhaul by surfacing
  `review_manifest.json` status in `catalog.json`, generated tutorial
  front matter, generated and human tutorial hubs, and GUI tutorial hover text.
- Completed the second Phase 4 tutorial review slice: catalog/generated
  tutorial outputs now carry dependency-aware stale reasons, feedback issue
  template links, and the two May 18 human reviews are recorded as stale after
  the subsequent tutorial source changes.
- Began Phase 5 of the tutorial presentation overhaul by carrying declared
  `graphics[]` metadata through catalog/manifest/report outputs, embedding the
  promoter-design generated TFBS SVG beside its tutorial step, and tracking one
  existing screenshot-backed dotplot tutorial figure as review-staleness input.

## 2026-05-22

- Added a Phase 0 tutorial coverage audit for the decimal-numbering overhaul,
  identifying uncataloged tutorial-like Markdown pages and the TP73/VKORC1
  promoter-luciferase decision gate before schema or rename work begins.
- Began Phase 1 of the tutorial presentation overhaul with v2 catalog/manifest
  schemas, v4 source-unit grouping/graphics fields, a 10-group tutorial
  taxonomy, and two newly cataloged guided walkthroughs for stateless sequence
  inspection and TFBS similarity ranking.
- Completed the Phase 1 tutorial group-assignment review gate: all 40 source
  units now use `gentle.tutorial_source.v4`, 38 tutorial units carry derived
  decimal ids, and the generated hub plus landscape overview remain grouped but
  intentionally unnumbered.
- Began Phase 2 tutorial presentation regrouping: the human tutorial hub,
  generated reference hub, Help tutorial picker, and `File -> Open Tutorial
  Project...` menu now present tutorials by content group and decimal id rather
  than flat lists or tier-only buckets.
- Began Phase 3 tutorial filename convergence: numbered tutorial sources,
  generated chapters, and hand-written walkthroughs now use decimal
  `<GG>-<PP>_<unit_id>` filenames, while unnumbered reference units keep
  explicit `ref_*.json` source names.
- Archived superseded TP73 promoter/luciferase tutorial pointers and moved the
  planned VKORC1 / Factor X companion tutorial note into the roadmap parking
  lot instead of active tutorial numbering.
- Added deterministic helper vector-card and catalog-doctor inspection routes
  with GUI Planning-window parity for metadata-only vector candidates.

## 2026-05-21

- Added metadata-only Luebeck S1 vector/resource candidates to the helper
  catalog, with explicit sequence availability, redistribution, biosafety,
  procurement, host-compatibility, marker/origin, and empty-backbone fields.

## 2026-05-19

- Added GUI export for lab-assistant cloning handoff reports and upgraded
  `ExportLabAssistantInstructions` to `gentle.lab_assistant_instructions.v2`
  with ODT/DOCX editable reports plus embedded lineage graphics where
  rasterization is available.
- Made `gentle_examples_docs tutorial-check` failures include paste-ready
  tutorial feedback context with chapter/source/workflow/artifact paths and a
  suggested issue-template category.
- Added `gentle.tutorial_review_manifest.v1` review metadata with nonfatal
  tutorial-check warnings for missing, unknown, or stale tutorial review
  entries.
- Regenerated executable tutorial chapters with review/automation provenance
  front matter and a standardized tutorial feedback section.
- Added a GUI Tutorial help button that copies issue-ready feedback context
  with tutorial id, source JSON, workflow, artifact, version, platform, and
  current search context.
- Grouped the in-app Tutorial topic picker by catalog audience buckets and
  added standardized feedback sections to hand-written tutorial pages.
- Added the read-only `planning consult cloning` route and
  `gentle.planning_cloning_consultation.v1` report so Agent Assistant and
  ClawBio can quote deterministic cloning-strategy/vector advice instead of
  improvising prose.
- Pinned `planning consult cloning` v1 to the 11 catalogued routine families,
  kept `--seq-id` traceability-only, and left marker/promoter/MCS gaps as
  explicit questions instead of narrative ranking heuristics.
- Clarified the Agent Assistant tutorial's same-live-state relationship with
  the GUI and GUI Shell, and added `CreateSequenceFromText` /
  `sequence create --sequence-text ...` so agent-proposed sequence text can be
  materialized as an ordinary project sequence.
- Added direct OpenAI usage and billing hyperlinks to GUI quota-failure
  messages from the Agent Assistant, including `Test Setup` live-probe results
  whose static configuration is otherwise available.
- Tightened the native Agent Assistant system prompt so providers are told to
  suggest only GENtle shared-shell commands, avoid invented `fs.*`/gateway-style
  verbs, ask for local file paths when discovery is outside the GENtle contract,
  and use ASCII punctuation for egui-safe rendering.
- Documented the current non-OpenClaw command boundary for the in-app Agent
  Assistant, including Finder/Spotlight as the macOS path-discovery fallback
  before handing a selected file path to GENtle.
- Added a protocol-owned GENtle-local slash alias registry covering `/help`,
  `/list`, GUI file opening, exact file import, pasted sequence creation, and
  explicit external fetch aliases while keeping vague filesystem/OpenClaw-style
  commands rejected with typed alternatives.
- Added `Ctrl+Return` submission for the Agent Assistant prompt editor.
- Added a `Copy Response JSON` button for the latest Agent Assistant response.
- The Agent Assistant suggestion table now validates suggested commands before
  enabling `Run`, marking hallucinated or unsupported commands as invalid.

## 2026-05-18

- Split the tutorial hub into hand-written guided walkthroughs versus
  machine-generated executable reference chapters, with reciprocal links for
  generated chapters that have guided counterparts.
- Added native Mistral Agent Assistant integration with a quick-start GUI route,
  `MISTRAL_API_KEY` session/env handling, model discovery, and live setup
  preflight parity.
- Added tutorial stewardship rules so recurring maintenance checks tutorial
  implementation drift, wet-lab readability, biological teaching depth, and
  Codex-versus-human review sign-off boundaries.
- Kept the Agent Assistant and Agent Interfaces tutorial pinned in the
  in-app `Help -> Tutorials` picker even if catalog discovery falls back to
  markdown scanning.
- Enriched the first generated tutorial batch with v3 per-step CLI snippets,
  expected-result callouts, and prerequisite links.
- Enriched the Gibson generated tutorial batch with v3 per-step CLI snippets,
  expected-result callouts, and prerequisite links.
- Enriched the online generated tutorial batch with v3 per-step CLI snippets,
  expected-result callouts, and prerequisite links to reference preparation.
- Completed generated tutorial-source enrichment so every executable chapter
  now has v3 per-step CLI/expected-result teaching fields.
- Polished generated tutorial rendering with quieter provenance in Help,
  applied-concept wording, earlier parameter explanation, CLI state guidance,
  SVG label fallbacks, and capped tutorial-project progress before final open.
- Surfaced documentation-only guided walkthroughs, including the Agent
  Assistant and Agent Interfaces tutorial, from `File -> Open Tutorial
  Project...` so readers can find them from the tutorial-opening menu as well
  as `Help -> Tutorials`.
- Fixed `File -> Open Tutorial Project... -> Guided walkthroughs` to read all
  `manual/reference` tutorial entries from `docs/tutorial/catalog.json`
  instead of showing only the hard-coded agent tutorial entry.

## 2026-05-17

- Added v3 tutorial-source teaching fields for generated chapters, including
  per-step CLI/expected-outcome rendering, prerequisite links, and online
  local-execution callouts while keeping v2 sources loadable during migration.
- Added a reporter construct handoff tutorial that shows how to generate and
  inspect `gentle.reporter_construct_handoff.v1` from a saved
  promoter-reporter candidate set before running any macro commands.

## 2026-05-16

- Restored shared-shell reachability for the reporter catalog/recommender
  family so `reporters list`, `reporters recommend`, and
  `reporters export-corpus` now work through GUI Shell and
  `gentle_cli shell ...` as documented.

## 2026-05-15

- Regenerated the GUI/CLI/MCP parity matrix from `CapabilityDescriptor`
  surfacing metadata and added a freshness guard plus regeneration script.
- Replaced binary capability adapter exposure with per-adapter surfacing states
  and added `EngineError.cause_chain` preservation across adapter boundaries.
- Made engine-backed missing adapter routes render as parity gaps unless a
  curated not-applicable override supplies a real justification.
- Added `gui-menu` glossary surfacing and GUI menu/palette oracle checks so the
  parity matrix can distinguish first-class GUI affordances from GUI shell reachability.
- Added conservative RNA-read isoform triage TSV export for known-isoform,
  ambiguous, gene-supported/no-call, and off-target/bad-seed read bins without
  calling novel isoforms.
- Added Splicing Expert isoform read-support inspection using saved RNA-read
  mapped-isoform triage counts, with audit/export links back to contributing
  aligned reads.
- Fixed Agent/Routine Assistant hosted-window focus lookup so embedded macOS
  windows are raised by their stable hosted ids instead of stale title layers.
- Added an offline reporter-recommender V1 with a provenance-gated local
  reporter catalog, deterministic constraint ranking, rejected-candidate
  reasons, and JSON/JSONL corpus export for local AI retrieval or training prep.
- Added a read-only reporter construct handoff plan that joins saved
  promoter-reporter candidate reports, the offline luciferase recommender, and
  the existing `allele_paired_promoter_luciferase_reporter` macro without
  creating constructs automatically.

## 2026-05-14

- Made `AlignSequences` accept inline ASCII `SequenceScanTarget` operands for
  stateless pairwise alignment while preserving legacy `query_seq_id` /
  `target_seq_id` workflows.
- Added a protocol-owned capability registry for CLI/MCP/JavaScript/Lua
  discovery and normalized CLI/MCP adapter failures into structured
  `EngineError` payloads.
- Added a short in-tree agent development loop with `scripts/dev-gentle-cli`,
  `gentle_cli doctor --agent`, and a dedicated `docs/agent_dev_loop.md`.
- Hardened MCP capability discovery so `tools/list` carries glossary-backed
  descriptions, input/output schema metadata, mutating/external safety flags,
  and a documented exclusion ledger for shell commands without dedicated tools.
- Added generated tutorial chapters that still need explicit human functional
  confirmation to the tutorial catalog and Help -> Tutorials picker with the
  `generated+checked/human-pending` status.

## 2026-05-13

- Selected the TP73 genome-anchored evidence viewer as the next release aim and
  added the offline proof workflow, public runbook, fixture provenance, and
  generic genome-track detail polish needed for the first smoke path.
- Added a GUI/CLI/MCP parity audit and brought the glossary-backed help catalog
  in sync for implemented agent preflight/model/planning routes plus
  construct-reasoning graph routes.
- Fixed hosted Promoter design foregrounding so DNA-owned auxiliary windows now
  open above their parent sequence window instead of behind it.
- Added a dedicated macOS `egui` viewport investigation pack plus a separate
  source folder for the minimal repro harness.
- Added an optional `--expression-tsv` heatmap overlay to isoform architecture
  SVG exports, keyed by `isoform_id`/`sample_label` values.

## 2026-05-12

- Cleared the completed interim-release block from `docs/roadmap.md` and reset
  the roadmap to choosing the next release aim.

## 2026-05-11

- Revised the GUI egui stack from `0.34.1` to `0.34.2` across `eframe`,
  `egui`, `egui_extras`, and the `gentle-gui` helper crate.
- Added the ClawBio `mode=intents` runtime descriptor surface, trigger-keyword
  generation/drift checks, and a provider-neutral default scope for
  `services handoff`.
- Added `rna-reads show-alignments` for headless batch export of the same
  per-read alignment display payload used by `rna-reads show-alignment` and the
  GUI splicing-expert detail pane.

## 2026-05-10

- Moved sequence file-loading dispatch into `dna_sequence::load_from_file` so
  engine, JS, Lua, and tests no longer depend on `GENtleApp` for parser
  selection.
- Added a TP73 pancreas proof-bundle `audit` handoff to
  `scripts/tp73_pancreas_rna_mapping.sh`, including gene-group snapshot capture
  and expanded evidence-bundle paths for release review.

## 2026-05-09

- Aligned the `docs/release.md` Local Pre-Tag Smoke Checklist with the roadmap
  release gate by adding `cargo check -q` before the release build matrix.
- Split roadmap maintenance into `roadmap.md` for next work, this changelog for
  completed work, and `decisions.md` for durable implementation decisions.
- Added [`maintenance_chore_plan.md`](maintenance_chore_plan.md) to define
  five recurring repository hygiene chores: session close, daily bug scan,
  weekly drift/test/provenance scan, release-gate readiness, and monthly
  decisions/parity review.
- Added `scripts/maintenance_chore.py session-close` as the first runnable
  maintenance chore implementation, covering roadmap size, completed-history
  bullets, invariant links, conflict markers, worktree/artifact warnings, and
  plan-fidelity handoff evidence.
- Polished the session-close runner after external review: conflict-marker
  detection now uses the exact seven-equals Git separator, split-doc existence
  is checked explicitly, and plan fidelity is reported as a manual reminder.
- Added `scripts/maintenance_chore.py drift-scan --base-ref REF` as the first
  weekly drift/test/provenance automation slice, warning when shared
  engine/protocol/shell-contract files change without obvious test evidence in
  the same diff.
- Polished `drift-scan` after external review by narrowing the warning to
  shared-contract code without separate test-file evidence, recognizing plural
  Rust `_tests.rs` files, normalizing rename status display, and documenting
  inline-test and pure protocol-doc limits.
- Added `scripts/maintenance_chore.py release-gate` as the first
  release-readiness automation slice, checking compact roadmap release-gate
  structure, local release-gate links, release-facing document presence, and
  `cargo check` alignment with `docs/release.md`.
- Hardened the release-gate scan after external review so a missing roadmap
  returns a structured failure and titled Markdown links do not produce false
  broken-link findings.
- ClawBio/OpenClaw `gentle-cloning` capabilities replies gained shared
  UI-intent catalog handoff support through `ui_intent_catalog`,
  `ui_intent_catalog_error`, and `kind = ui_intent` suggested actions.
- Stash-pop conflicts around UI-intent roadmap/docs wording were resolved by
  preserving both shared GUI catalog consumption and ClawBio/OpenClaw handoff
  documentation.

## 2026-05-08

- Release planning narrowed pre-release work to TP73 pancreatic Nanopore cDNA
  proof artifacts, release-signoff tests, GUI smoke, and documentation
  alignment.
- Gene-agnostic pancreas RNA screening was routed through
  `scripts/pancreas_gene_rna_screen.sh` for seed-only cohort triage and
  optional retained-read alignment.

## 2026-05-02

- Restriction-site inspection gained pinned Description-panel details plus
  `Copy Summary` and `Copy Detail JSON` actions for protocol notes and agent
  handoff.
- Restriction-site context menus gained pin/copy actions for hovered or
  selected labels without changing the rendering layer.

## 2026-05-01

- Restriction-site hover tooltips began reusing the shared
  `RestrictionSiteExpertView`, including enzyme grouping, cut geometry,
  tooltip lines, and REBASE links.
- `restriction_site_detail` exposed the same restriction-site detail record
  through MCP for non-rendered agent inspection.
- GUI Agent Assistant suggested-command execution and shared plan execution
  reject nested `agents ...` invocations across the full ask/plan/execute
  family.

## 2026-04-28

- RepeatMasker/UCSC `rmsk`-style repeat annotations gained deterministic
  labels and subtype colors across GUI maps, feature-tree grouping/filtering,
  and SVG export.

## 2026-04-21

- Repeat/mobile-element reasoning gained direct Dotplot and RevComp Dotplot
  actions from the DNA-window inspector.

## 2026-04-18

- Construct-reasoning annotation candidates gained accepted/rejected/locked
  review state, shared shell/CLI mutation routes, GUI curation, and JS/Lua/MCP
  parity.
- Accepted or locked construct-reasoning annotation candidates can be written
  back as ordinary sequence features through a shared engine report.
- Construct-reasoning graphs gained portable annotation-candidate summaries so
  long genomic views can show collapsed context instead of raw overlap lists.
- Similarity/repeat/mobile-element reasoning began emitting generated
  repeat-region annotations and operational-risk facts for PCR, mapping,
  nanopore, and cloning review.

## 2026-04-17

- Construct-reasoning overlays were narrowed to annotation-grade and
  decision-linked spans so raw restriction/TFBS evidence no longer floods long
  genomic windows by default.
- Construct-reasoning graphs gained a portable annotation-candidates layer
  shared by the sequence overlay and inspector.

## 2026-04-16

- Adapter/linker restriction-capture reasoning was added to construct
  objectives, GUI summaries, run bundles, shared shell/CLI, JS, Lua, and MCP.

## 2026-04-13

- Variant markers became first-class construct-reasoning evidence with
  promoter/TFBS/CDS/UTR/splice effect hypotheses and first assay-family
  suggestions.
- Routine planning began consuming transcript-aware variant summaries and
  construct-reasoning context for `routines list`, `routines explain`, and
  `routines compare`.

## 2026-04-12

- Construct-reasoning contracts reserved host/helper/growth fields and began
  emitting non-sequence context facts and decisions.
- Host-route restriction/methylation, selection/complementation, helper-backed
  selection, and growth-condition interpretation were added to shared
  construct-reasoning outputs.
- Process run-bundles gained a portable `construct_reasoning` section for
  ClawBio/OpenClaw offline inspection.
- DNA-window inspection gained non-sequence construct-reasoning summaries while
  keeping overlays sequence-only.
- A starter host-profile catalog and GUI browser were added.

## 2026-04-10

- Construct-reasoning overlays gained hover/click inspection, side-panel
  evidence details, and GUI filters without mutating the engine-owned graph.

## 2026-04-09

- Construct-reasoning graph records and metadata storage were introduced for
  deterministic read-only reasoning over existing sequence evidence.
- The DNA window began auto-refreshing a read-only reasoning graph and painting
  selected reasoning spans through a `Reasoning` layer.

## 2026-03-19

- Async BLAST job durability was hardened with persisted job snapshots,
  deterministic restart/reload recovery, cancellation semantics, and
  conformance coverage.
- Prepare-job cancellation and timeout classification became engine-owned and
  deterministic.
- BLAST and related genome routes gained external-binary preflight diagnostics
  for `blastn` and `makeblastdb`.

## 2026-03-04

- Malformed GTF/GFF annotation reporting began summarizing non-fatal parse
  issues with file/line context in prepare/report payloads.

## 2026-03-03

- JS adapter parity tests began comparing shared-shell-backed import-pool,
  REBASE sync, and JASPAR sync outcomes.
- Lua adapter parity tests began comparing import-pool, REBASE sync, and
  JASPAR sync outcomes.
- CLI forwarded-dispatch parity tests began verifying import-pool and resource
  sync routes against shared-shell execution.

## 2026-02-27

- MCP UI-intent parity tests began comparing discovery, prepared-query,
  latest-selection, and open-intent outputs against shared-shell `ui ...`
  execution.

## 2026-02-24

- Shared-shell execution tests began covering local-fixture REBASE and JASPAR
  resource sync without network dependency.

## Undated Historical Baseline Imported From Roadmap

- Shared engine operation/workflow execution exists in `src/engine.rs` with
  adapters for GUI, CLI, shared shell, JS, Lua, Python wrapper, and MCP.
- Build tooling defaults to per-worktree Cargo target directories, while
  release-installer CI avoids stale `target/` caches for tag/manual builds.
- Multi-crate extraction has started around `gentle-protocol`,
  `gentle-engine`, `gentle-render`, `gentle-shell`, and `gentle-gui`, with
  portable contracts and selected render/shell/gui helpers already moved.
- Engine, engine-shell, and GUI module decomposition has begun, including
  extracted tests, helper modules, and navigation-oriented module docs.
- Shared shell and CLI cover capabilities, state summary, genomes, helpers,
  resources, tracks, ladders, candidates, services, protocol cartoons,
  RNA-read routes, CUT&RUN routes, sequencing confirmation, UI intents, and
  many direct rendering/export paths.
- Python integration exists as `integrations/python/gentle_py`, wrapping
  deterministic CLI contracts.
- ClawBio/OpenClaw integration exists as a copy-ready `gentle-cloning` skill
  scaffold with local checkout and Apptainer launchers, example requests,
  reproducibility bundles, PNG-first graphics bundling, lifecycle-aware
  suggested actions, and first-run bootstrap documentation.
- MCP server baseline exposes deterministic tool routes including
  `capabilities`, `state_summary`, guarded `op`/`workflow`, `help`,
  UI-intent tools, prepared-genome helpers, BLAST async tools, and
  restriction-site detail.
- GUI Agent Assistant and typed planner/execution boundaries exist; machine
  callers are expected to prefer typed `agents plan` / `agents execute-plan`
  style contracts over chat-only flows.
- GUI command palette, menus, shell, agent, MCP, and ClawBio discovery now
  consume the shared UI-intent target catalog instead of maintaining separate
  target lists.
- Feature expert views for TFBS, restriction, splicing, isoform panels, and
  related SVG exports are engine-owned and exposed through multiple adapters.
- Candidate-set scoring/filter/set algebra, guide-set operations, practical
  guide filtering, oligo/protocol export, and design-constraint filtering are
  first-class engine operations.
- Reference/helper genome catalog, preparation, extraction, BLAST, anchored
  track import, subscription, and cache cleanup surfaces are implemented across
  engine/shell/CLI/GUI paths.
- Genome-anchor resolution supports deterministic fallback policy controls.
- BigWig/VCF/BED/BLAST genome-track imports and display filters are shared
  engine/display contracts.
- SnapGene `.dna` import is supported through a reusable headless parser crate
  with synthetic fixtures alongside GenBank/EMBL/XML parity fixtures.
- XML import is additive and GenBank-like in semantics, with XML dialect
  rejection/parity tests.
- Process run-bundle export, protocol-cartoon rendering, lab-assistant
  instruction export, routine decision traces, and protocol template export
  have shared record/export baselines.
- Cloning routine catalog, macro templates, macro validation, routine explain
  and compare, Routine Assistant GUI baseline, macro lineage visualization, and
  baseline restriction/Gibson/Golden Gate/Gateway/TOPO/TA/GC/In-Fusion/
  NEBuilder HiFi packs are present.
- Gibson specialist, Gibson preview/apply, arrangement/rack/ladder state,
  rack-profile authoring, label export, fabrication/export, OpenSCAD export,
  and rack simulation JSON paths are present.
- Virtual pool gels, protein gels, protease digest gels, 2D protein gels,
  arrangement-aware gel exports, and text band rows are present.
- Primer design, qPCR assay design, primer specificity assessment, Primer3
  backend baseline, backend preflight, and progress parity exist, with richer
  constraint parity still tracked as future work.
- PCR, advanced PCR, mutagenesis, overlap-extension mutagenesis,
  insertion-primer design, restriction-cloning PCR handoff, transcript-derived
  cDNA assays, and qPCR setup/promotion helpers are present.
- UniProt, Ensembl gene/protein, GenBank, dbSNP, and related projection/import
  paths exist with persisted reports and provenance for several artifact
  families.
- Alternative-splicing interpretation baseline exists with transcript lanes,
  boundaries, event/evidence summaries, primary map read-only mode, and expert
  SVG parity.
- RNA-read interpretation includes seed/preflight workflows, retained-read
  alignment, saved report stores, GUI mapping workspace, CLI/ClawBio examples,
  target-gene cohort summaries/audits, score-density export, and batch mapping
  helpers.
- CUT&RUN support includes dataset status/prepare/project/read interpretation,
  ROI read reports, coverage/export, regulatory-support inspection, and GUI
  inspector baselines.
- Isoform architecture panels, curated TP53/TP73 panel resources, online/offline
  workflows, expert SVG exports, protein/domain lanes, and validation reports
  are present.
- TFBS/JASPAR resources include catalog listing, remote metadata sync,
  deterministic registry benchmarks, single-entry inspection, promoter-window
  summaries, multi-gene promoter views, score tracks, score-track similarity,
  and query resolution.
- Promoter and variant context work includes variant promoter summaries,
  reporter-fragment suggestions, allele materialization, luciferase planning
  examples, and regulatory evidence ledgers.
- Dense rendering controls, adaptive linear DNA letter routing, scroll/zoom
  policy, display visibility, topology, viewport controls, and GUI declutter
  baselines are present.
- Agent recursion guardrails block nested `agents ...` execution directly and
  through macro expansion paths.
- Test-data provenance policy, remote-resource skip policy, source-documentation
  pilot docs, and cargo-doc module documentation goals are established.
- Screenshot artifact support exists historically but remains disabled by
  security policy unless explicit opt-in/approval is restored.
