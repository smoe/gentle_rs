# GENtle Changelog

This file is the auditable record of completed work that used to accumulate in
[`roadmap.md`](roadmap.md). Keep entries brief and outcome-focused. Durable
implementation rules belong in [`decisions.md`](decisions.md), not here.

Maintenance rule:

- Do not add completed-work bullets to [`roadmap.md`](roadmap.md).
- At the end of each session, move completed roadmap items here.
- Prefer a date and one short outcome sentence. Include command names,
  document names, schemas, or feature names only when they help a reader
  understand what changed.

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
