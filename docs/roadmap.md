# GENtle Roadmap

Last session: 2026-10-18 - Integrate gentle_rs_2_main's bounded rack-safety and source-backed tube references into the [next-check handoff and return path for Glen](glen_next_verification_handoff.md), alongside the pending sequence-design follow-up; carrier support remains undecided. Continue development on main under the owner's [release-only coordination cadence](decisions.md#development-and-release-verification-cadence): no interim SHA freezes, recurring bundles or synchronized acceptance waits. Normal checks record their tested source without blocking every subsequent change. Keep the [six-item ledger](glen_audit_followups_20261008.md), remaining full-suite, native/live-agent/scroll/package and scientific gates for actual release preparation; earlier evidence is not a newer-source pass. No local Rust, Claude, push, new dispatch, publication or paper/output/history edits; `.11` DMG remains deferred.

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

Published baseline: `v0.1.0-internal.12` at `5e6201a69a599cba87491c86cfd7311cb6fee173`, published 2026-10-13 with all three native `package-opt1` packages, four JSON receipts and GHCR publication. Preserve the [release ledger](release_notes/release_notes_v0.1.0-internal.12.md#publication-evidence) and separately verified [non-publishing headless receipt](audits/glen_followups_20261008/container_validate_5e6201a6/README.md), including the initial Ubuntu setup failure and open independent native acceptance; the `.10`/`.11` scientific ledgers stay historical.
Current candidate: `v0.1.0-internal.13` is the development version, with no selected release SHA or publication. Validate version/tutorial metadata through the [Glen handoff](internal_13_rollover_glen_handoff.md); consolidate fresh `package-opt2` native/container/extracted-package and full acceptance at actual release preparation, not every intermediate commit. Root `Cargo.toml` owns compiler compatibility; the sequence-design crate has no independent Rust 1.85 floor. Retain the owner-completed dependency update's recorded lock hash and older verdicts at their actual source/profile boundaries. External Rust verification remains required without becoming an interim development freeze.
Active aim (owner, 2026-10-15): review early containers/pools/arrangements/racks/plates and 3D exports; ship minimal provider-neutral protein-expression and explicit collection-to-storage paths with tutorials. Require one candidate's green Windows/macOS/full-Linux CI, build-only installers/container with extracted smokes, Glen's rollover/tutorial/readiness verdict and owner approval. Native latency, real-reference primer-import performance, profile comparisons, live-agent and historical scientific ledgers carry forward without blocking `.13`; retain the focused Linux worker result without inferring Windows/full-suite acceptance.

Retained scientific story, not a selected `.13` aim: source-bound TP73/PATZ1 evidence, conservation navigation and approval-bound reporter plans. Designed products are not validated assays; similarity does not prove sufficiency.

Smoke status: headless proofs do not replace graphical acceptance. All three
Linux chapters (branch/reverse complement, digest, compact Simple PCR) must pass,
without waiving the ten-minute budget. Glen reports green at `3de4803a`; retain the newer Simple-PCR evidence at `9070cb1d` separately. Recheck at the merged candidate and improve the result capture's clipped lower rows without treating five saved pairs as five visibly complete rows.
Glen validates Linux and BLAST (missing-tool skips are not acceptance); CI owns Windows/macOS and containers.

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
- Published [`.10` release notes and pending exact-candidate gate ledger](release_notes/release_notes_v0.1.0-internal.10.md).

Release acceptance:

- The proof workflow remains offline-safe and writes non-empty sequence,
  splicing-expert, TFBS SVG, repeat-materialization JSON, CUT&RUN regulatory
  JSON, and CUT&RUN coverage/cut-site/fragment TSV artifacts.
- Version metadata, LF/CRLF tutorial replays, capability/parity checks, and the
  deterministic proof workflow pass the pre-tag validation recorded in the
  versioned release notes.
- Require Windows CI on every push plus one selected Linux/macOS job; before tagging, retain both Unix verdicts at the same SHA. Confirm the `.11` retained tutorial-version replay on native Windows, not only locally.
  Track job `107075366344` at `e4c94bf3`'s legacy SHA-1 CertUtil probe sharing violation separately; it was identical in both outputs, not the parity failure. The historical CLI-lock recheck is recorded in the [archived handoff](https://github.com/smoe/gentle_rs/blob/7c52082c24fca85817d70956146623eef6f24607/docs/internal_12_candidate_preparation.md#corrected-ci-evidence), not a `.12` pass.
- Full UCSC `rmsk`, raw CEL, full SRA, and genome downloads remain optional
  external resources; CI uses tiny local fixtures without restricting the
  engine or GUI contracts to those fixture sizes.
- Bind all verdicts to one SHA, lockfile, toolchain/profile and input hashes; rerun gates on candidate changes, never combine nearby revisions' passes.
- For each new candidate verify the [package-opt2 recipe](release.md): opt-level 2, all LTO off, 256 codegen units, assertions/overflow checks on, unwinding, no stripping, incremental/debug off, one build job, five native/three container binaries and unoptimized helpers. Keep `package-opt1` unchanged. The owner authorizes this one recipe change, not builds/publication; fresh exact-source receipts and extracted smokes remain required. Glen compares same-SHA opt1/opt2 resources, startup/PATZ1 and output parity before any further optimization; `.12` package success is not a speedup verdict.
- Require fresh headless `runtime-cli` receipts when source/profile changes, preserving RNAPKIN's corrected helper lock, SVG/PNG and unprivileged no-network smokes. Released `.12` container publication passed at `5e6201a6`; that evidence does not accept newer source, independent native behavior or scientific results.
- Owner candidate order ([commands](internal_13_preparation_handoff.md#owner-procedure)): complete generated follow-up; owner-approved push of main including the integrated Ubuntu 24.04 pin before October 19; select the full SHA; `gh workflow run ci.yml -R smoe/gentle_rs --ref main -f platform=linux` (also macOS if unsampled); dispatch `release.yml` and `container.yml` with `tag=v0.1.0-internal.13`, full `candidate_sha`, `publish=false`; if hosted opt2 fails, select the separate opt1 fallback and repeat on a new SHA; Glen's rollover/tutorial/A-B-C/readiness; owner tag/publication approval; .14 rollover only after confirmed release. Hosted execution is pending; opt1 fallback remains unmerged, no Ubuntu 26 qualification.

Release cut line:

- Carry unrecorded `.10`/`.11` scientific and GUI acceptance forward without erasing historical ledgers; bind each new verdict to its explicitly selected source SHA.
- Preserve the successful `.12` CI-native package and container cycle at `5e6201a6`, including DMG/extracted-package smokes. Independent downloaded-package/native checks remain pending; later candidates require fresh receipts, Glen's readiness verdict and owner publication approval.

Outstanding acceptance scope (including [native Windows repairs and the expanded input audit](testing.md#september-2026-failure-analysis), and [independent anchor-to-assembly authority](architecture.md#source-coherent-locus-transcript-presentation)):
- Glen reports five-gene export QC at `d978b0c4` with LOCUS-column fix `d37bd427`; recheck the [selected annotated downloads](integrated_locus_tss_profiles.md), receipts and raw signal/motif provenance at the final candidate. Before offering EMBL, re-export into fresh directories and independently compare GenBank/EMBL/FASTA sequences and annotations. [Synthetic format-parity coverage](sequence_format_parity.md) does not replace his September 13 real-data acceptance. Run `dbsnp_spdi_` GenBank re-import regressions via Glen/CI; live dbSNP/native C/T acceptance remains separate.
- Run TP73/PATZ1 and promoter-matrix GUI smoke, including cancellation and [DNA-map/Structure source comparison](transcript_source_presentation.md#dna-map-and-structure-tab) with unchanged primer targets; obtain canonical
  CD44/TGFB1/SERPINE1 region sets before real-candidate analysis; verify overlapping TSS bands, zero-position peaks and per-TSS CUT&RUN/TFBS arrays against continuation-page IDs in the exact supplied bundle, never infer missing ROIs.
- Repeat the copied-state IRF9/Q00978 acceptance and require the one-assay panel
  to report one of three linked cDNAs and two of four linked records covered,
  with only the true patch record genomically unassessed.
- Decide `.13` retention/removal of deprecated ClawBio shell-normalizer modes
  using the [compatibility inventory](clawbio_gentle_integration_onepager.md#compatibility-review-inventory); `.11` eligibility did not authorize removal.

## Next Session Priorities

1. Have Glen recheck three introductory GUI chapters plus promoter-similarity navigation and TATA evidence review/materialization at the merged SHA, retaining outcomes/screenshots. Run the 10-step learner path in [synthetic TSS acceptance](tutorial/08-15_tss_collection_gui.md) first; retain Glen's `1fcee6c0` evidence separately from merged-candidate acceptance. Repeated-open reuse, Forget/Undo, restart and stale-member editing remain separate contributor checks, not part of that learner verdict. Verify compact validation, cascade and first-member native focus (including mixed open/new members) without harness refocusing; multi-window clutter and individual-member navigation remain follow-ups. Exercise the [TSS collection workflow](tss_workspace.md#acceptance) on TP73: membership, strands/flanks, separate-process approval, cache refresh, DeltaNp73-only cropped locus, validated scans, repeated opening and inner-agent handoff. Compare symbol/ID memberships with mixed or case-variant gene labels; verify historical collections stay intact while new grouping requires fresh approval. Include stored-collection inspection and confirmed stale-registry removal without deleting sequences. Check unsupported-index warnings and locus-level unassigned transcripts remain distinct from loaded-feature coverage and gene-specific exclusions.
2. Carry historical `.10`/`.11` acceptance into one named `.13` candidate: locked tests, five native package entrypoints plus separate optional scripting checks, tutorial/parity, IRF9/Q00978, benchmarks, circular ORF/feature regressions and TP73/PATZ1 TFBS/reporter/CUT&RUN/locus proofs with exact JASPAR 2026 PFMs. Keep both tagged-source ledgers unchanged.
   Keep tutorial 08.13's live GUI capture and real Parquet acceptance separate from the [synthetic TSS tutorial replay](tss_tfbs_profiles.md#rescoring-existing-reports) and mocked tests. For private five-gene rescoring, first obtain the old scored report/receipt and exact matrix/context bindings; the [historical input inventory](tss_tfbs_profiles.md#check-input-completeness-before-regeneration) alone is not acceptance.
3. Preserve the `.12` [published package ledger](release_notes/release_notes_v0.1.0-internal.12.md#publication-evidence); do not redispatch that completed CI artifact cycle as unfinished work. After explicit push/dispatch approval, require fresh macOS/Linux/Windows, headless and build-only package evidence for any new candidate. Keep independent native GUI/tutorial checks and source/profile identity separate; preserve all existing tags/assets.
4. Have Glen refresh the `a5b01843` receipt with native interval/score and catalog-bound assembly checks; preserve verified figures unless mismatched. Validate the [candidate-first view and peak-centred follow-up](cofactor_candidate_interface_plan.md); keep offline regeneration, live Conservation/partner-evidence proof and ClawBio compatibility pending removal review.
5. Have Glen run the [promoter-cofactor tutorial](tutorial/08-14_promoter_cofactor_browser.md) against the exact candidate: live stale-form/copy binding, saved evidence, verified plus/minus feature preview/confirmation/undo, then the optional IRF9 package. Human scientific approval remains pending; unstranded evidence remains regions and raw coverage reconstruction is not supported by the reduced package.
   Validate DuckDB boundaries/limits, scan filters and selected-package GUI staleness; have Glen inspect [imported-hit triangles](tss_tfbs_profiles.md#imported-duckdb-motif-hits) on both strands and run the [regulatory/TSS-subset reader](regulatory_motif_subset.md) against the actual genome-wide delivery (including MT and source floors); synthetic fixture acceptance is not production acceptance. Native TSS follow-up (not an added release gate): live acceptance of [report attachment](gui.md#attach-the-quantitative-report) and [view SVG export](gui.md#export-the-native-tss-view), including CD44/TGFB1 gaps, local-strand legends, curve colours, span-only navigation, missing units, axes, selected lanes and sparse-hit coverage; reporter rows remain open. Validate approved [local scoring](gui.md#compute-local-tss-scores): InlineSequence/matrix provenance, separate score kinds/scales/cache, cancellation, stale-form handling and genuine zero versus unavailable masks. Its agent/headless routes now exist (`ui open tss-view --local-score`, `promoters tss-view-svg`) with a real-locus walkthrough in [tutorial 08.17](tutorial/08-17_tp73_dnp73_factor_curves.md); retain Glen's Linux/manual 08.16/08.17 screenshots at `8652ecf2` separately from merged-candidate equivalence, live hosted-intent adoption and Windows/macOS acceptance per [handoff](tp73_dnp73_factor_curves_glen_handoff.md). Span-scoped headless local scoring is implemented pending tests/acceptance: complete footprints inside the requested span, full-window admission, explicit span counts/hashes and no hidden GUI/report rescoring. Native TSS export follow-ups are implemented pending tests/acceptance: disclosed default cap of 128 hover titles per lane/strand without dropping curve points or validity bands, top-three within-matrix raw positions, documented prepared-reference P1 access and validated collection SVG/index/hash receipts with per-window/per-lane scales (no new calibration). Retain current paint/file-size refusals; Glen measures the revised title counts and checks receipts, cancellation and GUI/headless equivalence before any limit increase. New live database queries still require exact reference-bound interval admission; retain [flat-file limits and readability plan](tss_annotated_sequence_readability_plan.md). Glen: regenerate [source-coherent Ensembl/RefSeq pages](transcript_source_presentation.md) from one enriched locus and verify mappings/PDF/SVG/HTML receipts; synthetic coverage is not real-data acceptance. Keep the [standalone comparison](transcript_source_presentation.md#standalone-comparison-review) deferred rather than adding a second source-join implementation.
6. Continue the [responsive-command migration](asynchronous_command_execution_plan.md#reconciled-codex-plan-2026-09-14): validate BLAST start-v2 deferred diagnostics, bounded FIFO/cancellation and observation-only status at the merged revision. Preserve DEC-026 merging and [DEC-048](decisions.md#dec-048-command-submission-and-observation-never-wait-for-execution); audit general managed-command queuing/dependencies and result/file/job semantics before expanding beyond the BLAST queue.
   Have Glen validate genome preparation, Run/slash, live navigation/conversation and cancellation at one revision; extend the synthetic growing-journal cost probe to a representative project before claiming bounded UI lock latency. Keep cross-process control, durable admission and other adapters pending. No private T2T/model run in development; general asynchronous migration is not complete.
7. Carry into `.13` the [gene-informed primer-PAIR study GUI](gene_assay_study_gui_plan.md), not single-primer capture. Have Glen annotate the authentic-reference [PATZ1 tutorial 04.08](tutorial/04-08_gene_assay_study_gui.md) at G1-G6, including the 13-Ensembl/four-RefSeq source comparison and unresolved discrimination matrix; authentic whole-reference specificity and order approval remain pending. Verify reference-bound specificity, companion selection receipts, standalone multi-reference import/manifest/parity and legacy fail-closed regressions at the merged SHA using [tutorial 04.09's verification commands](tutorial/04-09_reference_bound_primer_specificity.md#verification-scope); Rust execution and real-reference acceptance remain with Glen/CI. Validate read-only `primers transcript-assay-specificity-matrix` receipt history, binding failures and unlinked oligo citations at that SHA; fan-out remains a [conditional Slice B](primer_design_specialist_plan.md#conditional-slice-b-bound-panel-fan-out), with import profiling first and separate approval; readiness attachment remains a later policy decision. The [reconciled specialist plan](primer_design_specialist_plan.md) does not authorize new GUI/filtered-view scope. Supply an original typed paper-facing dossier/source bundle. Include coverage-scope interpretation and last-page navigation for large transcript/assay panels in that live check; synthetic pagination tests do not establish real-locus GUI performance. Retain the historical [reconciled `.11` scope](internal_11_plan.md#review-reconciliation) at the new candidate: package parity, exact-SHA TP73/DeltaNp73/tutorial acceptance and independent tutorial/vector-PDF verification. Keep TSS managed-worker migration and capture-retention requirements optional; defer capture specificity/thermodynamics/order readiness until lab inputs are retained and calibrated cofactor ranking until anchor/expression/background contracts exist. Retain separate RefSeq accession matching and real-image gel acceptance; automatic peak detection remains deferred.

Current non-goals:

- Do not add private proposal or unpublished grant content to this repository.
- Do not start broad engine extraction, GUI redesign, or infrastructure work
  before it is tied to the selected next-release story.
- Do not add GUI-only business logic for convenience; keep shared
  engine/shell/protocol contracts as the source of truth.
- Do not promote speculative assistant ideas into protocols without
  confirmation.

Useful session close:
- Name the next-release aim or state that it remains undecided.
- Keep completed outcomes in [`CHANGELOG.md`](CHANGELOG.md) and durable rules in [`decisions.md`](decisions.md).
- Leave unrelated local paper/figure work untouched unless explicitly requested.

## .13 Priorities
Owner-selected aim, 2026-10-15; use the [release notes](release_notes/release_notes_v0.1.0-internal.13.md) for the acceptance ledger. Minimal working paths, not breadth or a new performance claim:
- A: review early containers, pools/aliquots, arrangements/gels, racks/plates/templates and 3D exports plus early cloning/analysis features; verify the authored 1,536-slot guards and [source-backed Safe-Lock reference models](../data/labware/README.md) with Glen/CI, including legacy carrier byte baselines and native review. Carrier support/pitch, 96-position interpretation and storage boxes remain separate owner decisions, not new `.13` acceptance.
- B: finish the shared provider-neutral protein-expression draft and tutorial 09.03; preserve explicit provider selection, no submission/order/contact and conservative product/readiness boundaries. Inspector is stretch-only; ranking stays later.
- C: verify the minimal `CreateGeneSetArrangement`/generic-op path and [registered tutorial 03.10](tutorial/03-10_gene_set_storage_projection.md): explicit complete-set/source-tube binding, one serial/default storage rack, provenance and existing exports. Sorted fallback IDs repair Glen's reported `b2b6b30a` parity failure; all four storage and both reconstruction regressions pass in Linux default-feature CI at `1c5be317`. Headless filters, CLI/parity and `scripts.test_release_workflow_tutorials` with a built CLI, generator freshness, fresh exports and native acceptance remain Glen/CI work; no implicit placement, new freezer template or GUI mapper.
- Release blockers: green Windows/macOS/full-Linux CI, build-only native/container/extracted smokes, Glen's rollover/version/tutorial/readiness checks and owner approval at one selected candidate; package-opt2 first, prepared opt1 fallback if hosted builds fail. Ordinary [CI at `1c5be317`](glen_audit_followups_20261008.md#ordinary-ci-observation-at-1c5be317-2026-10-17) is terminal/failure: Linux/Windows hit the same README metadata guard after 4,182/4,121 root-library passes; `57a2403a` restores the label and adds an offline guard, not yet externally accepted. Optional audits skip; no release freeze is selected.
- Non-blocking carry-forward: [S6/native latency](dna_feature_rendering_latency_plan.md), [real-reference primer import](glen_post_blast_import_phase_trace.md), same-host profile comparisons, live agent and the [remaining native/human checks](glen_final_candidate_execution_handoff.md#remaining-native-and-human-checks). No earlier evidence becomes a current-source pass.
- Historical scientific/FLNA gates remain open at their true sources, without blocking `.13`; preserve human teaching-exon/version-label decisions and the inactive draft, plus the [`.10`](release_notes/release_notes_v0.1.0-internal.10.md)/[`.11`](release_notes/release_notes_v0.1.0-internal.11.md) ledgers. VKORC1: focused Linux repair passed at `6b7109db`; all five [stripped-child families now pass Linux/Windows at `1c5be317`](glen_audit_followups_20261008.md#completed-windows-job-and-terminal-run-2026-10-18). New-source public/native and physical GUI acceptance remain with the owner/Glen. Package/extracted-installer acceptance stays separate; no earlier verdict accepts newer source.

## Phase A: AI Communication And Safety Plane

Keep MCP, agent, shell, ClawBio/OpenClaw, and UI-intent surfaces thin and
adapter-equivalent. The next work is breadth and safety, not new biology:
expand deterministic tool coverage over shared shell/engine routes, harden
mutating-intent confirmation consistently across agent, voice, and MCP paths,
and preserve the shared `ui ...` intent catalog across menu, command palette,
shell, agent, MCP, and ClawBio discovery. Useful work here is parity,
discoverability, confirmation wording, and `ui_intent`-shaped operator handoff;
defer autonomous planning policy or biology-specific reasoning that is not
backed by current engine records. Carry forward to `.13` Glen's validation of explicit Splicing Expert open/focus/close and the [transcript-assay playbook](ai_task_playbooks.md#playbook-12-gene-informed-primer-pairs-and-isoform-discrimination) with a live inner agent: real feature IDs, assay objectives, coverage and separate QA gates. This is not a retroactive `.10` gate or order approval. Run `tutorial-check` for declared command/confirmation parity before model review; live review remains for model quality, voice and broader MCP parity.

For ClawBio/OpenClaw specifically, keep new integration work centered on
descriptor/runtime parity (`INTENTS.json`, `mode=intents`, examples, and
trigger-keyword drift checks) plus explicit scope/presentation contracts rather
than adding more biology-specific wrapper modes. Following 08.15's six TSS descriptor routes, extend the existing 08.16/08.17 cases plus the new `promoters tss-view-svg --collection ID OUTPUT_DIR` route through the same delegate with artifact-write/selection approvals, using the settled per-window scale policy and ordered SVG/index/receipt contract; retain model-quality and GUI-host acceptance separately. After the first GENtle-side
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
candidates and yield/folding/chassis questions without creating constructs.
The v2 shared report separates the canonical neutral delivery request from
optional provider projection. Supplied reviewed requirements and annotation
evidence populate the draft; missing facts remain unanswered, never vendor
example defaults. Provider binding requires explicit selection, permitted
outsourcing, an unambiguous product and catalog capability/allowlist checks.
No usable product means no draft; quote suggestions also require local preflight.
V1 stays readable; external Rust checks and static checks do not prove wet-lab yield.
Keep text, JSON, CLI and agent consumers on the same report and action rows.
Do not add submission, ordering or automatic biological optimization.

Keep the next synthetic-biology bridge slices in this order:

1. Verify the read-only Protein Expression Handoff inspector and tutorial
   09.03 on the integrated source: shared-report parity, explicit selection,
   typed withholding, stale-result rejection and unchanged project/history.
   Glen/CI must execute authored Rust/CLI checks and native window review;
   offline source checks do not establish runtime or native acceptance.
   Include all seven `protein_expression_inspector_` adapter regressions.
   Preserve absent product/requirements as unanswered, not inferred defaults.
   Keep actions copy-only; edits and close must discard stale worker results.
   Check shared UI-intent registration and canonical tutorial/parity freshness.
   Reporter sibling views and ranking remain separate follow-ups.
2. Rank chassis, vector/tag
   routes, and provider candidates from structured biology, provider
   capabilities, and local inventory. Emit score components and rejection
   reasons; do not silently codon-optimize or create constructs.

MCP follow-up: `tools/list` is curated, not a one-command/one-tool mirror.
Keep `docs/agent_interface.md` as the exclusion ledger; add tools only with
stable typed schemas and mutating/external safety semantics. Prioritize the
remaining parity-matrix gaps by user workflow, not tool counts or exclusions.

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

The first collection mapping now exposes construct-reasoning inspection actions
across project sequences and homogeneous gene-set resolutions without persisting
batch-built graphs. Structured GUI evidence panes have also landed. A later
phase stays intentionally separate: broader mapping of additional
single-sequence operations across sequence sets. Portable inspection actions,
cross-adapter parity, curated repeat-family provenance, and objective-aware
severity presentation are recorded in `docs/CHANGELOG.md`.

## Phase D: Visualization And Workflow UX

Use the GUI as the human inspection surface for engine-owned evidence. Continue dense DNA-map readability, alternative-splicing polish, gel/feature editing, contextual links, visual regression fixtures and scroll/zoom hardening when they fit the selected release story.
Improve inspection clarity, deterministic exports and manual-smoke reliability; defer unrelated visual redesigns. Rendering-latency work belongs to [.13 Priorities](#13-priorities), not here. Next [contextual availability](contextual_action_availability_plan.md): Glen's live native/hosted check of item-first lineage context menus (especially clicked-DNA PCR binding), compact palette guidance and Edit-menu/palette Undo/Redo; prioritize map/result/empty-space context-click guidance next. Native-menu and agent/host-probe integration remain deferred.

- Evidence-viewer follow-up: keep the Splicing Expert evidence ledger readable
  on large real loci and add direct report-store selectors if probe/expression
  reports become persisted project objects.
- Cryptic-splicing v2 follow-up: consider explicitly versioned segmented spans
  for origin-crossing circular sequences, non-GT/AG boundary policies, and an
  opt-in typed ORF-inference contract. Do not fold those inferences into the v1
  annotation-bound CDS or evidence-layer semantics.
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
- Agent GUI acceptance: test screenshot consent on Linux/Xvfb and native macOS;
  review Assistant/setup translations and small-window layout, especially
  Chinese/Japanese fonts and Latin. Native child capture must match its viewport;
  hosted/embedded fallback captures the GENtle root. Consider a content-layer
  crop only if review finds the selected hosted window insufficiently clear.

## Phase F: Interpretation Later

Later interpretation work includes image/sketch-to-state patch proposals,
richer assistant-facing findings records, deeper wet-lab process modeling,
repeat/mobile-element curation, and cohort-scale comparative reasoning. These
must remain explicit, inspectable proposals with confirmation before mutation.
Useful work here improves proposal formats, evidence records, reproducible
cohort comparisons, and rollback paths for state patches; defer autonomous wet-lab conclusions or unconfirmed mutations.
## Parking Lot
- Experimental [synonymous motif/GC design](dna_sequence_optimization_plan.md) stays outside `.12`; after the pinned [repair/GC acceptance](sequence_design_acceptance.md#executed-audit-at-de9df8af), retain the [first-slice audit](glen_sequence_design_handoff_20261005.md) unchanged and require current-source checks, native Windows/macOS/Pi, usefulness/live-model and biological verdicts separately. Verify the private core's standalone compiler contract before extraction/publication; do not restore an untested floor or infer a workspace minimum from core-only tests. Dedicated GUI, codon objectives and richer specifications need their own contracts; bounded Linux/X11/process evidence implies no general performance or suitability claim. Run the new pruning/tie-break, exhaustive one-/two-symbol IUPAC, 256 fixed-seed differential, 4-6 six-synonym and match-limit retention regressions through the focused core tests; these authored Rust checks and Cargo check remain unrun locally, without an interim acceptance freeze.
- Accession-pinned [TSS TFBS profiles](tss_tfbs_profiles.md#rescoring-existing-reports): Glen to rescore original five-gene inputs after the inclusive-tail correction, then repeat PDF pixel/hash, shared-scale and sparse-peak acceptance; the `44b73e4a` audit does not cover this scoring method. Keep the `.10` gate; native wizard and independent reference verification remain follow-ups.
- Give TP73 Simple-PCR live acceptance a representative small fixture or bounded fast preset; exact interactive defaults exceed the ten-minute smoke budget.
- Replace whole-transcriptome FASTA materialization in UniProt linked-transcript
  inventory construction with indexed/streaming retrieval, and expose bounded
  progress during large specificity finalization joins.
- Extend the feature-gated semantic GUI test identifier registry beyond the
  bounded TP73/PATZ1 DNA-viewer, Splicing Expert, and Agent Assistant release
  controls as additional native acceptance stories require it.
- Have Glen run the candidate-bound tutorial gate for `smoke` plus conservation's
  view-only `offline-core`; extend typed coverage next to Gibson and dotplot.
  MCP/Feature Editor, transcript-assay, vendor and trace guides need Glen's review; keep biological verdicts separate.
- Extend the transcript follow-up with provenance-complete array evidence for positive differential selection and approved-study execution.
- Extend inner-agent context with readiness-bound tutorial actions and explicit headless receipt replay; deferred-job observation is part of active DEC-048 work. Keep approvals and undo history unchanged.
  Defer skill activation and inverse discovery until readiness is safeguarded.
- Optional OS credential-store persistence for Agent Assistant API keys.
- Separate opt-in inner coding-agent mode for Codex, Pi, or comparable local
  harnesses, with explicit workspace scope, tool permissions, diff review,
  bounded tests, and rollback. Do not inherit Agent Assistant command-execution
  permissions or present model prose as an applied source change.
- Optional tiny generation probe for quota verification.
- Implement opt-in Mistral conformance (`docs/testing.md`); skip without `MISTRAL_API_KEY`.
- Audit restriction-catalog cut offsets (especially Type IIS) and extend usage annotations beyond REBASE; recheck older cloning absence/uniqueness claims before laboratory use.
- Floating restriction-site detail popover/window if the Description panel is
  too easy to miss.
- Audit general linear-SVG feature clipping (inclusive-end helper versus half-open feature bounds); retain broader plasmid-map presentation parity and dense selected-site polish.
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
- External-service provider/CRO integration once deterministic local contracts are stable enough to wrap.
- Luebeck S1 vector-resource follow-up: resolve sequence/procurement URLs and
  redistribution clearance before bundling actual sequences for the metadata
  candidates.
- VKORC1 / warfarin Factor X follow-up tutorial idea: after the maintained
  VKORC1 rs9923231 promoter-reporter tutorial, design a separate companion
  story for downstream warfarin-relevant Factor X target-site constructs once
  the construct family, controls, readout, and bench-facing success criteria
  are stable enough for a validated walkthrough.
- Bind the transcript-aware SERPINE1 architecture comparison to the actual
  laboratory insert/vector sequence and prepared Rostock CUT&RUN manifests
  before any construct assignment, occupancy comparison, or materialization.
  The public pinned SERPINE1 plus clearly synthetic evidence test proves the
  integrated presentation only; it does not substitute for those materials.
- Primer-walking support and iterative read/contig data management.
- Helper-construct terminology migration away from legacy "helper genome" wording.
- CUT&RUN controls/replicate comparison, calibrated differential support, formal
  peak calling, and whole-genome read processing beyond the ROI-first release scope.
- Add viewport-bounded quantitative Ensembl/ENCODE regulatory-signal lanes, retaining assay/biosample/condition/replicate/units and separate activity rules.
- Primer3 parity, virtual PCR/off-target filtering, multiplex tiling, LAMP, and
  allele-specific assay families.
- First-class global RT-PCR specificity/redesign workflow building on the
  landed local assembly-aware primer variant screen and combining prepared
  genome and transcriptome searches, intended-product/paralogue policy,
  junction evidence, Primer3 plus exact-run QC, ranked redesign candidates,
  and one provenance-bound review bundle.
- GuideRNA off-target ranking and macro-template packaging.
- Cross-tool parity synthesis for Serial Cloner, MacVector, and SnapGene.
- Measured gels: [verified ladders and real-image/native-editor acceptance](gel_image_analysis.md#remaining-work) for the manual editor included in `.10`; automatic peak detection remains deferred.
- Weekly/monthly maintenance chore automation rollout from [`maintenance_chore_plan.md`](maintenance_chore_plan.md).
- First regular release: local-first GENtle server with optional browser interaction and shared project/job ownership; deferred beyond internal releases, with WebAssembly portability considered separately.
- Expand the landed policy-driven collection surface beyond promoter
  derivation, primer specificity, construct-reasoning inspection,
  restriction-site and TFBS hit scanning,
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
- Extend the [evidence-guided boundary selector](reporter_fragment_selection.md) with real called-peak/cohort/GUI acceptance, independently verified control/replicate metadata, calibrated local support and native-cloning simulation;
  keep broad-window enrichment separate from insert-local proof and approved downstream plans.
- Extend regulatory-fragment exact design products with explicit cloning-reaction
  simulation and final-product QA; retain live GUI acceptance of exact-product
  review, digest approval, receipt inspection and one-step undo before release.
- Extend the promoter regulatory-partner screen beyond its exact motif-tuple
  ledger with explicit matched/background cohorts, calibrated enrichment and
  multiple-testing correction, local ontology/hallmark interpretation, and
  distal-element evidence. Preserve the landed per-hit coordinates and
  replayable decision traces as the audit layer for those later statistics.
- Add paralog mapping/resolution only with a concrete engine consumer,
  preserving relationship provenance and one-to-many ambiguity without
  inferring functional equivalence from matching symbols.
- Annotation-independent RNA ROI capture and any reusable sparse-posting resource remain proposed; see [origin status](rna_read_origin_sparse_index_plan.md).
- Synonymous restriction-site introduction with verification-digest proof, distinct from motif-avoidance sequence design.
- Report-first curated feature-library auto-annotation with a separately confirmed apply step.
- Full glossary/help-generation inversion so `docs/glossary.json` becomes a
  generated or validated projection of engine/protocol-side descriptors. Resolve duplicate registry `(source, name)` identities and derive a typed parity backlog before regenerating the matrix; see [verified gap](introspection.md#12-deferred-follow-up).
