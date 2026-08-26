# Release Notes / Changelog: `v0.1.0-internal.10` (draft)

| Release field | Value |
| --- | --- |
| Status | Draft through 2026-08-14, pending automated pre-tag validation |
| Target date | August 2026 |
| Previous tag | `v0.1.0-internal.9` (2026-06-05) |
| Primary story | Genome-anchored TP73 evidence viewer and transcript-aware primer workflows |
| Manual GUI smoke | Deliberately deferred until after this interim release |

This internal release covers the work after `v0.1.0-internal.9`, tagged on
2026-06-05. Once the new tag is cut, the auditable Git comparison is:

```text
v0.1.0-internal.9..v0.1.0-internal.10
```

[Open the GitHub comparison after the tag is published.](https://github.com/smoe/gentle_rs/compare/v0.1.0-internal.9...v0.1.0-internal.10)

The main release story is a genome-anchored TP73 evidence viewer. GENtle can
open the GRCh38.p14 TP73 locus and inspect exon and transcript structure,
repeats, projected microarray intervals, prepared CUT&RUN intervals,
paired-read ROI support, regulatory/motif-context reasoning, TFBS annotations,
and coordinate/build provenance in one project. The same proof is reproducible
through a committed headless workflow and inspectable through GUI, CLI/shared
shell, MCP, and agent-facing paths.

This is also a substantial infrastructure cut. It expands transcript-aware
primer and isoform-assay planning, fact-aware capability introspection, local
and external agent integration, probe-region and RNA evidence handling,
gene-set/promoter reasoning, review-gated material handoffs, and the amount of
expensive GUI work that runs safely away from the egui thread.

## Highlights

- **Genome-anchored TP73 evidence viewer:** a public, offline-safe proof joins
  GRCh38.p14 sequence context with repeat, Clariom-style array, BED, TFBS, and
  splicing evidence, with deterministic SVG/JSON artifacts and a manual GUI
  runbook.
- **Microarray and probe-region bridge:** explicit Affymetrix CEL/APT planning,
  helper-output inspection, plotting, coordinate projection, PM-probe input,
  transcript/exon geometry interpretation, and Clariom D GUI support now form
  a shared engine/shell workflow rather than a GUI-only path.
- **Splicing and RNA evidence:** junction probes are visible in the DNA and
  Splicing Expert views; transcript/exon and exon/exon summaries gained clearer
  support/provenance presentation; RNA workflows gained target-region rescue,
  allele-aware hash screening, persisted disjoint exonic parts, and DEXSeq
  annotation/count exports.
- **Transcript-aware primer workflows:** bounded and progress-aware Primer3
  searches, explicit transcript-coverage universes, whole-cDNA similarity
  guidance, terminal-exon RT-primer pools, related-sequence common primers,
  exhaustive cDNA/genomic specificity checks, and digest-bound multi-gene
  isoform-study execution now share auditable engine contracts.
- **Fact-aware automation:** project facts, capability requirements, expected
  effects, readiness evaluation, and post-run effect verification now cover
  the shared capability registry broadly. JSON help and capability
  introspection derive from the same descriptors.
- **Inner and outer agents:** Agent Assistant project guidance, shared UI
  intents, local OpenAI-compatible providers, Ollama/Msty guidance, Codex Local
  discovery, visible local `/help` and `/list` results, and ClawBio/OpenClaw
  request paths were hardened without creating a second biology command
  surface. Fresh sessions offer direct open-document and Configuration actions;
  typed suggestion preconditions now dim and block dependent rows until the
  live project fact exists, including a fail-closed check for automatic runs.
  DNA-window intents also expose their queued/opening state in Agent Assistant;
  annotation-rich lazy windows postpone their feature browser until requested.
  Explicit `Agent help` controls can attach the exact GENtle viewport
  being diagnosed with a local preview-before-send; optional macOS native
  full-window capture is permission-gated and restricted to GENtle-owned
  windows. Image-capable agents may also explain why they need one screenshot;
  GENtle then requires one-shot approval for a user-selected registered window
  and still previews the result before the separate send action.
- **Codex Local usability:** GENtle can discover the current macOS ChatGPT
  bundled Codex executable, presents a model selector immediately, reads the
  visible local Codex model metadata, and forwards an explicit choice through
  `codex --model`.
- **Promoter, gene-set, and construct reasoning:** offline gene-set producers,
  promoter cohorts, ortholog promoter comparison, relationship expectations,
  CUT&RUN support summaries, repeat-family interpretation, and portable
  construct-reasoning inspection actions are available through shared
  contracts.
- **Review-gated experimental handoffs:** first-class oligo order forms,
  duplicate/reuse review, delivery routing, protein-expression planning, and
  Metabion/GeneArt-oriented preflight/quote packets keep planning outputs
  separate from confirmed purchasing actions.
- **GUI responsiveness:** detached revision-checked engine work, cheap display
  checkpoints, history-excluding read-only worker clones, presentation caches,
  virtualized dense lists, slower idle repainting, and background network,
  track, BigWig, and model-discovery work substantially reduce idle and active
  GUI overhead.
- **Runtime and dependency hardening:** process-local runtime activity is
  visible through `introspect runtime`, MCP, SIGUSR1, and the GUI status bar;
  local fingerprints use SHA-256; legacy SHA-1 download checks are delegated
  to optional platform tools; and the macOS Objective-C bridge uses the current
  `objc2` API family.

## Notable Changes by Area

### 1. TP73 Genome Evidence and Presentation

- Added the release-proof workflow
  `docs/examples/workflows/tp73_genome_evidence_viewer_release_proof.json` and
  the public `docs/tp73_genome_evidence_viewer_runbook.md`.
- Added tiny, provenance-documented fixtures that deterministically produce
  three repeat features, four projected array intervals, two prepared/projected
  CUT&RUN BED features, four concordant paired-read fragments in two support
  windows, V3 regulatory-support JSON, TFBS annotations, and TFBS score-track
  output without a prepared whole genome or SRA download.
- Added headless queries for repeat, array, and BED qualifiers so an inner
  agent, MCP client, or CLI user can verify the same evidence shown in the GUI.
- Improved linear and Splicing Expert SVG legends, provenance text, vertical
  layout, support-frequency color coding, and long-text containment.
- Added inspectable/copyable feature details for repeats, array intervals, and
  BED tracks, including coordinates, score/statistical fields, source files,
  projection state, and assembly context.
- Added a sequence-visibility control, clearer array-effect color semantics,
  and coordinate-aligned array placement in both the live DNA view and linear
  SVG export.
- Added junction-probe evidence to Splicing Expert and the DNA sequence viewer,
  using split probe geometry and connection arcs where appropriate.
- Improved exon/intron continuity and phase cues, including CDS-phase flank
  coloring and tooltips for frame interpretation.

### 2. Array and Probe-Region Evidence

- Added `arrays probe-regions` as a shared preflight for explicit CEL,
  metadata, annotation/library, platform, dependency, contrast, output, and
  cache checks.
- Added explicit R/oligo and Affymetrix Power Tools handoff paths without
  embedding vendor analysis logic in the GUI.
- Added inspection, native SVG rendering, direct and mapped coordinate
  projection, and GUI controls for completed helper output.
- Added optional sample metadata and true PM-probe matrices, preserving sample,
  condition summary, logFC, probeset/transcript-cluster, exon, and coordinate
  provenance.
- Added transcript/exon geometry interpretation with explicit overlap,
  junction, ambiguity, score-basis, and review-status records. These are
  evidence-triage aids, not automatic isoform calls.
- Added a registry-first Clariom D platform preparation path and a small
  E-MTAB-14704 TP73 validation fixture with deterministic provenance.

### 3. RNA Read Interpretation

- Added a standalone target-region rescue screen for finding reads overlooked
  by an initial mapping while preserving explicit target and read provenance.
- Added the deterministic `gentle.rna_allele_hash_screen.v1` route and
  standalone `gentle_cli allele-hash-screen` adapter for small allele-aware
  RNA-read screens.
- Persisted disjoint exonic-part partitions and added mutually joinable DEXSeq
  flattened-annotation GFF and strict two-column count exports.
- Improved transcript exon-path handling and reduced RNA mapping/inspection UI
  repaint and presentation work.

### 4. August Primer and Isoform-Assay Work

- Added deterministic, digest-bound planning and execution for single- and
  multi-gene isoform-assay studies. Two explicit approval stages bind the
  normalized biological request and the exact ordered operations before shared
  project state changes.
- Added atomic multi-gene execution, content-bound checkpoints, separately
  approved exact-prefix reuse, and opt-in per-gene continuation. A failed gene
  rolls back to its own boundary; external output paths remain fail-closed
  because filesystem writes cannot be revoked by an in-memory rollback.
- Added sequence-aware bounded Primer3 plans with annotation-anchored allowed
  and rejected intervals, deterministic work budgets, pre-launch feasibility,
  candidate counters, cancellation, and progress streaming. GENtle detects
  `--progress`, can request an updated Unix status with `SIGUSR1`, and retains a
  legacy path for unextended Primer3 binaries.
- Added capability-aware Primer3 discovery. A progress-capable build is
  preferred, then the highest version and newest executable; an explicit user
  override remains authoritative. Windows `PATHEXT` discovery and deterministic
  fake-binary tests cover the cross-platform behavior.
- Added explicit transcript-coverage universes, including content-bound UniProt
  isoform inventories. Coverage and isoform distinction remain separate, and
  unresolved mandatory targets stop strict searches before Primer3.
- Added content-bound whole-cDNA resources and similarity maps that guide
  Primer3 toward lower-risk intervals. Variant, repeat, and cDNA-similarity
  evidence can be reviewed as hard or soft design constraints without replacing
  final complete-cDNA or genomic specificity gates.
- Added Primer-BLAST-style review material: bounded junction-side matching,
  complete-primer realignment, narrowly reviewed off-target products, stored
  alignment rendering, chemistry settings, and content-bound redesign proposals
  after specificity failure.
- Added deterministic common-primer design for related loaded sequences,
  complete pair-by-member product matrices, explicit strict versus best-effort
  coverage, and bounded Pareto alternatives. Existing selection and specificity
  gates are unchanged.
- Added `gentle.terminal_exon_rt_primer_pool.v1` for sequence-specific RT oligos
  composed from a fixed adapter and terminal-exon-derived variable segment.
  Strand-aware transcript mapping, bounded pool-wide interaction selection,
  persisted provenance, JSON export, and optional exhaustive local genomic
  BLAST evidence are shared across PCR Designer, shell/CLI, MCP, Python,
  JavaScript, and Lua.
- Improved exhaustive cDNA specificity import by consuming validated complete
  BLAST HSP alignments directly while retaining the established semiglobal
  fallback for partial or legacy rows. Subject-grouped product pairing,
  completeness, intended-target, allowlist, and provenance gates remain
  authoritative.
- Retained the linear cDNA evaluator as independently labelled fallback
  evidence only. It can help inspect very large handoffs, but it is never
  promoted into a native GENtle-finalized specificity pass.
- Added replayable selection provenance and clearer differential-junction
  explanations, including qualifying contrasts, thresholds, incremental and
  redundant detection value, and explicit validation obligations that do not
  claim isoform abundance.
- Added an approval-bound informative fallback for strict transcript-panel
  designs that are structurally infeasible. The original typed failure remains
  intact; any derived best-effort panel is visibly partial, budgeted, audited,
  and never presented as an optimal or complete strict solution.
- Expanded isoform-assay dossiers and publication bundles with readiness gates,
  blockers, evidence and report identities, order provenance, complete primer
  tables, locus/exon/gel figures, and explicit `resolved`, `pending`, or
  `unresolved` status for partially completed studies.

### 5. Introspection, Agents, and Interface Parity

- Added the first shared project fact graph and fact domains for sequences,
  reports, displays, persisted workflow objects, host tools, and GUI
  availability.
- Added capability preconditions, expected-effect modalities, argument
  binding, catalog versus bound readiness, and post-run verification for the
  shared operation/shell registry.
- Expanded fact-aware descriptors across sequence operations, genomes, BLAST,
  tracks, reports, primers, candidate/guide sets, containers/racks, macros,
  gene sets, external resources, UI intents, agents, and service handoffs.
- Added `introspect facts`, `introspect capabilities`,
  `introspect readiness`, `introspect verify-effects`, `introspect all`, and
  direct CLI introspection access.
- Projected built-in shell help and JSON help rows from shared capability
  descriptors while retaining `docs/glossary.json` as the transitional source
  for glossary wording.
- Added Agent Assistant controls for sequence-window selection and display
  visibility, while preserving the rule that selection/window intents do not
  imply sequence deletion or arbitrary execution.
- Hardened local model response parsing and Msty/Ollama configuration guidance,
  including Msty's gateway/model-server distinction.
- Added JS and Lua bindings for construct-reasoning inspection actions and
  maintained CLI/MCP/ClawBio parity checks over shared routes.

### 6. Gene Sets, Promoters, and Construct Reasoning

- Added first-class gene sets across protocol, engine, shell, GUI, CUT&RUN,
  promoter TFBS, lineage, and documentation paths.
- Added offline direct-list, ontology-assignment, and co-regulated gene-set
  producers with explicit source and relationship metadata.
- Added offline ortholog promoter cohorts and comparisons using local mapping
  resources, including species/genome-matched evidence states.
- Added co-regulated and anti-co-regulated relationship expectations with
  non-blocking divergence/concordance review flags.
- Added portable construct-reasoning inspection actions and shared dotplot
  execution routes, including objective-specific task severity and curated
  repeat-family corroboration.

### 7. Materials and External-Service Handoffs

- Added persisted oligo order forms with deterministic line ids, source primer
  or qPCR report provenance, duplicate/reuse groups, explicit review, and
  JSON/CSV export.
- Added delivery routing from selected sequences/spans, oligo forms, and primer
  report pairs, with review gates before quote preparation.
- Added read-only protein-expression handoff reports with sequence/CDS/tag
  context, high-yield questions, candidate routes, and separate GeneArt
  preflight versus quote-packet actions.
- Added ClawBio/OpenClaw request examples for provider diagnostics, Metabion
  oligo/m-block handoffs, GeneArt cloned-gene/protein-expression handoffs, and
  generic sequence-delivery routing.

### 8. GUI Responsiveness and Correctness

- Moved long mutating GUI jobs onto engine-owned detached forks and retained
  stale-result, display-merge, auxiliary-metadata conflict, revision, and undo
  semantics when committing results.
- Prevented detached and read-only workers from cloning inherited undo/redo
  history while preserving undo across live and detached operations.
- Made metadata-only track/planning edits dirty the project and invalidate
  stale redo checkpoints.
- Cached root-frame presentation state, linear rendering derivations,
  sequence-window display synchronization, and Splicing Expert models.
- Virtualized dense JASPAR, feature-tree, and sequencing-confirmation lists and
  made Splicing Expert transition caching sparse rather than exon-count
  squared.
- Moved GenBank, UniProt, Ensembl-protein, JASPAR, tracked-file, track auto-sync,
  BigWig preflight, and related I/O away from the egui thread with typed task
  completion and cancellation/stop-waiting behavior.
- Reduced forced repaint cadence for idle RNA workspaces and prevented
  unchanged per-frame genome-track serialization/synchronization.
- Added constant-time dirty-state and unchanged-display checks while retaining
  complete refreshes for structural, feature, topology, sequence, and visible
  presentation changes.

### 9. Runtime, Security, Dependencies, and Documentation

- Added `gentle.runtime_status.v1` and a process-local activity stack exposed
  through shared introspection, MCP, Unix SIGUSR1 dumps, and a compact GUI
  status indicator.
- Improved genome preparation progress with separate transcript/gene indexing
  and checksum-verification phases.
- Removed the in-process Rust SHA-1 implementation. Local generated identities
  use SHA-256; legacy external-download SHA-1 verification uses optional
  platform tools and reports unavailable verification explicitly.
- Updated the egui/resvg and macOS Objective-C integration stack and migrated
  native menu code to current `objc2` class-definition APIs.
- Reworked the repository README into a concise overview, moved detailed
  demonstrations to `docs/showcase.md`, and documented the complementary inner
  Agent Assistant and outer-agent/MCP/ClawBio paths with generated TP73 figures.
- Hardened Windows examples/tutorial execution against small-stack overflow and
  kept capability, glossary, MCP, fixture-provenance, and generated-tutorial
  checks aligned with the expanded command surface.

## Compatibility and Contract Notes

- The shared engine remains the source of biology and workflow behavior across
  GUI, CLI/shared shell, MCP, Python, JavaScript, Lua, and ClawBio/OpenClaw
  adapters.
- Most new schemas and fields in this interval are additive. Callers should
  continue to validate the advertised schema id/version and ignore documented
  optional fields they do not consume.
- ClawBio shell-normalizer compatibility modes deprecated in `.8` remain
  accepted in `.10` with warnings. Their earliest possible removal is `.11`,
  after a separate compatibility discussion.
- Agent-suggested commands remain subject to confirmation/execution policy;
  chat replies and clarification questions do not imply an executable action.
- Vendor/service output remains a reviewable handoff packet. GENtle does not
  submit a purchase merely because an order form, preflight, or quote request
  was generated.
- Microarray, repeat, TFBS, and CUT&RUN tracks are evidence layers. Their
  presence or geometric overlap is not itself a biological conclusion.
- Primer search intervals, similarity maps, Pareto alternatives, and fallback
  panels are design evidence. They do not replace native specificity
  finalization or experimental primer validation.

## Release-Facing Known Limitations

- This remains an internal release, not a stable file-format or public API
  commitment.
- The macOS bundle identifier remains the interim `com.example.gentle` value
  pending coordination with Magnus Manske on a stable project identity and
  URL. Changing it is intentionally deferred from this tag.
- The committed TP73 proof uses tiny local fixtures. Full UCSC `rmsk`, raw CEL,
  full SRA, BigWig, prepared genomes, and vendor resources remain optional
  external inputs.
- The vendor-derived Clariom identifiers/coordinates in the TP73 proof are
  paired with synthetic expression/statistical values; the fixture must not be
  interpreted as a biological differential-expression result.
- Probe-region transcript geometry is deliberately conservative and does not
  infer an expressed isoform from overlap alone.
- Codex Local model enumeration depends on a readable local Codex metadata
  cache. `Codex default` remains available and lets the installed CLI choose a
  model when no explicit model is selected.
- BLAST+, Primer3, ViennaRNA/RNAPKIN, `bigWigToBedGraph`, R/oligo/APT helpers,
  and legacy SHA-1 tools remain optional and are reported through preflight or
  runtime status when absent.
- Terminal-exon RT-primer melting temperatures remain descriptive. Global
  genomic uniqueness is reported only when the optional exhaustive genomic
  BLAST step was requested and completed against an identified prepared
  database.
- Independent linear cDNA fallback reports remain separately labelled evidence
  and cannot establish a native GENtle specificity pass.

## Suggested Pre-Tag Validation

Run the complete workspace suite when time permits, then retain the focused
release-story checks in the tag record:

```bash
cargo check -q
cargo test -q --test release_version_consistency
cargo test -q terminal_exon_rt_primer_pool --lib
cargo test -q transcript_assay --lib
cargo test -q primer3 --lib
cargo test --workspace --no-fail-fast
cargo run --quiet --bin gentle_examples_docs -- --check
cargo run --quiet --bin gentle_examples_docs -- tutorial-check
cargo test -q workflow_examples_tp73_cutrun_release_proof_writes_artifacts_and_features
python3 -m pytest -q tests/test_codex_agent_bridge.py
cargo run --release --bin gentle -- --version
cargo run --release --bin gentle_cli -- capabilities
git diff --check
```

Automated pre-tag validation results: _pending_.

## Post-Release GUI Smoke

Manual GUI smoke is deliberately deferred until after this interim tag. It
should then follow `docs/tp73_genome_evidence_viewer_runbook.md`, including
anchor/build status, repeat/array/BED/TFBS visibility and details,
feature-detail copy controls, Splicing Expert presentation, and the first 1200
bp linear viewport. Also open the Agent Assistant, select Codex Local, and
confirm that the model selector is visible before model discovery completes.
From Splicing Expert, click `Agent help` and confirm that Agent Assistant opens
with a preview labelled for that source window, that removing the image works,
and that no request is sent before `Ask Agent`. On macOS, right-clicking the
same control should either capture the complete GENtle window after Screen
Recording permission is granted or explain the missing permission without
affecting normal viewport capture.
Also ask an image-capable agent to request visual context. Confirm that its
reason/provider and registered-window selector appear before any capture, that
`Decline` produces no attachment, and that `Allow one screenshot` produces one
local preview but sends nothing until `Ask Agent`. Change the provider or close
the selected window before approval and confirm that the request expires.
For the August primer tranche, open PCR Designer and confirm that Primer3
preflight explains the selected executable, a terminal-exon RT-primer-pool
request exposes its fixed adapter and complete selected oligos, and a bounded
Primer3 run reports progress and remains cancellable when a progress-capable
binary is available. A machine without Primer3 or BLAST+ should instead show a
clear preflight limitation without disabling unrelated project work.

Manual GUI smoke: _deferred until after `v0.1.0-internal.10`_.
