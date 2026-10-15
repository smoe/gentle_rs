# GENtle `v0.1.0-internal.12` Release Notes

Status: [published prerelease](https://github.com/smoe/gentle_rs/releases/tag/v0.1.0-internal.12),
not a draft, at `5e6201a69a599cba87491c86cfd7311cb6fee173`.
The owner published it on **2026-10-13 at 01:02:04 UTC**; native asset upload
completed at 10:05 UTC after the installer rerun. Read-only GitHub metadata,
job results and the small release receipt were rechecked on 2026-10-13.
Publication closes the artifact cycle, not scientific, native-GUI, live-agent
or performance acceptance. The [`.11` notes](release_notes_v0.1.0-internal.11.md)
retain that release's separate history. No tag or asset was changed here.

## Publication Evidence

[Installer run 38330176948](https://github.com/smoe/gentle_rs/actions/runs/38330176948)
used the released SHA. Attempt 1's Ubuntu job
[115077994549](https://github.com/smoe/gentle_rs/actions/runs/38330176948/job/115077994549)
failed in `Set up Rust`, 45 seconds into the job; compilation and package
creation were skipped. This was not a compiler or package-smoke failure.
The successful macOS/Windows results retain their initial-attempt times in
attempt 2; only Ubuntu needed another build.

| Attempt 2 Job | Result And UTC Interval, 2026-10-13 |
| --- | --- |
| [Ubuntu 115222015196](https://github.com/smoe/gentle_rs/actions/runs/38330176948/job/115222015196) | Success, 09:17:54-10:04:13; entrypoint and extracted Linux package smokes pass |
| [macOS 115222017090](https://github.com/smoe/gentle_rs/actions/runs/38330176948/job/115222017090) | Success, 01:02:34-02:39:55; entrypoint and installer smokes pass |
| [Windows 115222017419](https://github.com/smoe/gentle_rs/actions/runs/38330176948/job/115222017419) | Success, 01:02:28-02:13:17; entrypoint and extracted Windows package smokes pass |
| [Validate packages 115240491353](https://github.com/smoe/gentle_rs/actions/runs/38330176948/job/115240491353) | Success, 10:04:16-10:04:51; receipts validated against the selected candidate |
| [Publish release 115240737681](https://github.com/smoe/gentle_rs/actions/runs/38330176948/job/115240737681) | Success, 10:04:54-10:05:44 |

The release contains exactly these seven assets:

- `gentle-v0.1.0-internal.12-linux-x64-package-opt1.tar.gz`
- `gentle-v0.1.0-internal.12-macos-arm64-package-opt1.dmg`
- `gentle-v0.1.0-internal.12-windows-x64-package-opt1.zip`
- `gentle-v0.1.0-internal.12-linux-x64.build.json`
- `gentle-v0.1.0-internal.12-macos-arm64.build.json`
- `gentle-v0.1.0-internal.12-windows-x64.build.json`
- `gentle-v0.1.0-internal.12-release-attributes.json`

The published [release-attributes receipt](https://github.com/smoe/gentle_rs/releases/download/v0.1.0-internal.12/gentle-v0.1.0-internal.12-release-attributes.json)
records all three archive hashes/sizes and all fifteen native binary identities,
source/workflow SHA `5e6201a6`, `mode=publish`, `profile=package-opt1`, and lock
SHA-256 `ead4421545a125267e33af84aeeea39f0e7ea2f60717a59ef0271d0e71855f1b`.
The receipt and GitHub asset hashes agree; no installer was downloaded or
independently executed for this documentation update.

[Release-event container run 38330177090](https://github.com/smoe/gentle_rs/actions/runs/38330177090)
also passed at `5e6201a6`, including the network-disabled headless smoke and
[Publish GHCR Image job 115078248256](https://github.com/smoe/gentle_rs/actions/runs/38330177090/job/115078248256).
It publishes `ghcr.io/smoe/gentle_rs:v0.1.0-internal.12`; the earlier
[tag-push run 38322028567](https://github.com/smoe/gentle_rs/actions/runs/38322028567)
remains a separate validate-only record, not publication evidence.

**Review boundary:** Glen's [final review](../audits/glen_followups_20261008/candidate_f35a58e6/README.md)
is at `f35a58e6d2bdd7e0a5551e4531c336d7949bfaf0`. Git comparison with the tag
shows identical `src/`, `crates/`, `packages/`, `Cargo.toml`, `Cargo.lock` and
`.github/`. Documentation, bundled tutorial resources and a Python checkout
test differ. His recorded results cover that matching Rust/workflow source
only within their original scope, not the exact downloaded packages or their
different help resources. Independent extracted-package, native-GUI, scientific,
live-agent and same-SHA `dev`/`package-opt1` performance checks remain pending.

## First Optimization Step

Native installers and the headless container select **`package-opt1`**, a
dedicated Cargo profile inheriting `dev`. It sets `opt-level=1` and `lto="off"`;
`lto=false` would allow local ThinLTO and is not equivalent. This is a controlled
first step, not a return to the former fully optimized release recipe.

- Keep 256 codegen units, no incremental compilation, no debug information,
  debug assertions and overflow checks enabled, and panic unwinding. Local
  `dev` builds remain unoptimized. Normal GUI/export content remains identical
  across profiles by contract; CI package smokes passed for the released SHA,
  while independent native/profile comparison remains pending.
- Keep one Cargo build job, the existing job allowance and resource diagnostics.
  Native packages contain the same five binaries; the container contains the
  same three headless binaries. JS/Lua remain optional source interfaces.
- Keep `cargo-bundle` and `rnapkin` on their pinned, unoptimized helper recipes.
  Preserve RNAPKIN's corrected lockfile and SVG/PNG smokes.
- Name archives with `-package-opt1`, stage from `target/package-opt1`, and bind
  the explicit optimization and safety settings into build receipts. Reject
  missing/mismatched settings rather than combining different recipes.
- Keep `bench-audit` and `release-fast` unchanged for historical comparisons;
  neither selects the new package build. Higher optimization, thin LTO, reduced
  codegen-unit counts, stripping and panic-abort are not enabled by this change.

## Learning And Agent Workflows

- Tutorial 04.08 now teaches **Find Primer-Pair Markers for PATZ1 Isoforms**:
  compare annotation sources, review the request, and explain why seven pairs
  cover 13 cDNA classes yet leave nine class pairs unresolved in the retained
  example. Planning/readiness/dossier steps are optional continuations. The
  endpoint-only feasibility command is no longer prescribed for this SYBR case.
- Approval-bound outer routing for TSS collections remains on the shared
  command/proposal path. Retained branch-local evidence does not imply live
  inner-agent or scientific acceptance of this new candidate.
- Tutorial review projections now use Git history and repository-relative
  paths, never checkout timestamps. Unavailable history is an explicit check
  error; stored catalog age follows the revision date, while live age warnings
  remain current-date-based. No human review or historical panel is relabelled.

## Historical Pre-Publication Evidence

**Validation pending:** historical synthetic transcript-panel reports retain
their true `.11` provenance. That is no longer a rollover error: the revised
contract checks historical integrity separately from a fresh `.12` result
replay, allowing only the two scoped generator fields and their verified ledger
hashes to differ. No stored report bytes or hashes were changed. The
[historical baseline result](../internal_12_tutorial_baseline_acceptance.md)
records the scoped comparison and remaining current-candidate/path checks. The replay evidence below covers only the
reported checks, not package or scientific acceptance.

**Reported replay evidence, 2026-10-02:** Glen reproduced the checkout failure
at `f1c904a4` and repaired the generated 08.03 chapter, tutorial hub and ledger
in `e462de1f`, now integrated. With Rust 1.96.1 he reports successful exact-commit
LF and CRLF replays, each covering parity, 51 examples and 29 chapters, plus
17/17 Python checkout regressions. The historical `.11` panel bytes and hashes
remain unchanged. The two stale human-review warnings remain expected. This
is externally reported evidence, not a local rerun or native macOS/Windows
Actions, `package-opt1` build, or scientific acceptance. The remaining handoff
checks and native-platform results still need to be recorded.

**CI follow-up, 2026-10-05:** the October 4 run `37235516663` at `f523a757`
finished with Linux and Windows catalog-test failures; their preceding drift
checks passed, Windows's full suite was skipped and macOS was not selected.
`19889dae` already corrects that catalog projection. The separate Git-history
guard, checkout mock repair and twentieth TSS case counts are locally validated;
they do not establish current native Actions or package acceptance. The later
run `37292637621` at `462b5ffe` also exposes a Windows rooted-path admission
failure in tutorial validation; at that inspection, its bounded repair still
required new native Windows CI, not a macOS-only verdict. The
[archived candidate-preparation handoff](https://github.com/smoe/gentle_rs/blob/7c52082c24fca85817d70956146623eef6f24607/docs/internal_12_candidate_preparation.md)
records that historical evidence and scope; its local checks are not current
`.12` acceptance. No workflow was dispatched for that handoff.

## Remaining Independent Acceptance

Current [build-only instructions](../release.md#build-only-candidate-verification),
[runner inventory and integration](../release.md#ubuntu-runner-decision),
and [compatibility review scope](../clawbio_gentle_integration_onepager.md#compatibility-review-inventory)
live in their topic owners. These pointers grant no acceptance or dispatch authority.

1. Build the `.12` generator and verify replay against the unchanged historical
   baselines. Regenerate only for a reviewed result change, never merely to
   relabel `.11` provenance. Run version, tutorial, LF/CRLF and packaging-policy
   gates on the committed candidate.
2. Retain the successful released-SHA native/container builds and CI package
   smokes above. Glen's independent extracted-package and native checks remain
   open; any newer candidate requires its own receipts. No workflow dispatch,
   push or tag mutation is implied here.
3. Have Glen compare `dev` and `package-opt1` at the same SHA/toolchain/host,
   recording build time/peak memory/package size plus repeated startup and
   PATZ1 open/resize timings. Keep `bench-audit` evidence separate. Decide on
   `opt-level=2` or thin LTO only after this measured result.
4. Carry forward the [roadmap](../roadmap.md)'s exact-candidate GUI, biological
   and external-tool checks. The existing linear-locus envelope (250 kbp /
   5,000 features) and circular correctness checks are not newly passed gates.

## CI Follow-Up And Owner Decisions

**Historical snapshot, 2026-10-07:** the sections below preserve the cleanup's
pre-publication scope, not current release status. The publication/package
results above supersede its pending build gates. Current `.13` priorities,
Ubuntu 24.04 integration and open acceptance remain with the [roadmap](../roadmap.md)
and [release process](../release.md); no gate is reopened or backfilled here.

The original candidate-preparation request was retired in `afbef444`; Git
history retains its exact command and evidence inventory. The newer
[historical pointer](../internal_12_candidate_preparation.md) remains intact. The Windows workspace pass at
`a51bbc06` (job `109794599955`) closes the historical JASPAR-lock recheck,
not acceptance of a later candidate. Preserve the separate CertUtil SHA-1 sharing-
violation note for job `107075366344` at `e4c94bf3`; no new verdict follows
from moving these links.

Retention/removal of deprecated ClawBio modes remains an owner decision.
The maintained [deprecation table](../clawbio_gentle_integration_onepager.md#deprecation-table)
and `DEPRECATED_REQUEST_MODE_EQUIVALENTS` in
`integrations/clawbio/skills/gentle-cloning/gentle_cloning.py` identify the
modes. Any future removal must audit that skill's `INTENTS.json`, examples,
routing, warnings and tests; an earliest-removal marker is not authorization.

The retired handoff recorded the October 19 rolling Ubuntu migration warning.
Runner selection was open at that historical review; the owner has since
integrated the Ubuntu 24.04 pin. Ubuntu 26 qualification remains separate.
This historical documentation section changes no workflow. [Build-only commands](../release.md#build-only-candidate-verification)
remain instructions requiring a frozen pushed SHA and explicit authorization.

## Carried Acceptance Scope (2026-10-07)

Moved from the roadmap to keep its front page concise. These are retained
requirements and reported historical evidence, not new gates, reruns or accepted
results. Optional follow-ups remain optional where stated.

Current Rust verification also remains pending for `b1c47ed6`'s transcript-CDS/
exon-skip repair: 1,353/16/1,337/1,571-bp, ambiguous identity, both-strand phase/
persistence and removed-CDS regressions. Retain tutorial 06.07 replay/capture,
exact UniProt/Splicing Expert and live FLNA/readability checks, sequence-design
provenance/tampering and GC reconciliation separately; moving this scope does
not turn experimental sequence design into a new `.12` gate.

Release story: inspect TP73 GRCh38.p14 exon, repeat, array, CUT&RUN, paired-read,
motif and coordinate-provenance tracks, plus PATZ1 negative-strand transcript/assay
evidence. Conservation adds saved-region request exchange, alignment navigation
and cancellation. Reporter plans bind evidence before approved atomic materialization.
Designed products are not validated assays; similarity does not prove sufficiency.

Smoke status: headless proofs do not replace graphical acceptance. All three
Linux chapters (branch/reverse complement, digest, compact Simple PCR) must pass,
without waiving the ten-minute budget. Glen reports green at `3de4803a`; retain the newer Simple-PCR evidence at `9070cb1d` separately. Recheck at the merged candidate and improve the result capture's clipped lower rows without treating five saved pairs as five visibly complete rows.
Glen validates Linux and BLAST (missing-tool skips are not acceptance); CI owns Windows/macOS and containers.

Proof path:

- [Headless regeneration workflow](../examples/workflows/tp73_genome_evidence_viewer_release_proof.json)
  loads `test_files/tp73.ncbi.gb`, overlays tiny local repeat and Clariom-style
  array fixtures, exercises CUT&RUN V1-V3 with synthetic paired reads, adds
  TFBS fixtures, and emits SVG/JSON/TSV proof artifacts.
- Public [TP73 evidence-viewer runbook](../tp73_genome_evidence_viewer_runbook.md).
- Offline [PATZ1 Gene Locus Evidence workflow](../examples/workflows/patz1_gene_locus_evidence_offline.json)
  and generated tutorial align transcript models, PSR/JUC effects, grouped
  occupancy, TP73 motif scoring, provenance, and assay continuation on one
  negative-strand axis.
- [Fixture provenance and regeneration notes](../../test_files/fixtures/evidence_viewer/README.md).
- Published [`.10` release notes and pending exact-candidate gate ledger](../release_notes/release_notes_v0.1.0-internal.10.md).

Release acceptance:

- The proof workflow remains offline-safe and writes non-empty sequence,
  splicing-expert, TFBS SVG, repeat-materialization JSON, CUT&RUN regulatory
  JSON, and CUT&RUN coverage/cut-site/fragment TSV artifacts.
- Version metadata, LF/CRLF tutorial replays, capability/parity checks, and the
  deterministic proof workflow pass the pre-tag validation recorded in the
  versioned release notes.
- Require Windows CI on every push plus one selected Linux/macOS job; before tagging, retain both Unix verdicts at the same SHA. Confirm the `.11` retained tutorial-version replay on native Windows, not only locally.
  Track job `107075366344` at `e4c94bf3`'s legacy SHA-1 CertUtil probe sharing violation separately; it was identical in both outputs, not the parity failure. The historical CLI-lock recheck is recorded in the [handoff](../release_notes/release_notes_v0.1.0-internal.12.md#ci-follow-up-and-owner-decisions), not a `.12` pass.
- Full UCSC `rmsk`, raw CEL, full SRA, and genome downloads remain optional
  external resources; CI uses tiny local fixtures without restricting the
  engine or GUI contracts to those fixture sizes.
- Bind all verdicts to one SHA, lockfile, toolchain/profile and input hashes; rerun gates on candidate changes, never combine nearby revisions' passes.
- CI must prove the new [package-opt1 recipe](../release.md) at one SHA: opt-level 1, all LTO off, 256 codegen units, assertions/overflow checks on, unwinding, no stripping, `incremental=false`, `debug=0`, one build job, five native/three container binaries (no JS/Lua). Keep Rust-built helpers unoptimized. Require fresh receipts and extracted-package smokes. Glen should compare same-SHA `dev`/`package-opt1` build resources and startup/PATZ1 latency before proposing opt-level 2 or thin LTO. Profile-independent UI/export content remains required; `.11` success does not certify `.12` or establish causes for historical compiler SIGKILL/runner shutdown.
- Re-run headless `runtime-cli` validation with `package-opt1`, preserving RNAPKIN's separately pinned corrected helper lock, early SVG/PNG smoke and unprivileged no-network smoke. `.11`'s successful container publication `36714939309` at `a51bbc06` is the baseline, not acceptance of this changed GENtle recipe. Retain exact source/profile/lock receipts and compiler RSS, pressure/swap and live logs; native packages and GUI acceptance remain separate.
- Owner decision: pin the [rolling Ubuntu hosts](../release_notes/release_notes_v0.1.0-internal.12.md#ci-follow-up-and-owner-decisions) before the October 19 Ubuntu 26 migration; Linux installer builds already pin `ubuntu-24.04`. No workflow change or Ubuntu 26 qualification is implied.

Release cut line:

- Carry unrecorded `.10`/`.11` scientific and GUI acceptance forward without erasing historical ledgers; bind each new verdict to the selected `.12` SHA.
- Require fresh `.12` native proof, including the [pinned macOS bundler preflight](../release.md), DMG and extracted-package smokes. The successful `.11` installer run `36714939304` at `a51bbc06` does not validate custom-profile bundling. Retain its tag/assets unchanged; require Glen's readiness verdict plus owner approval for later publication.

Outstanding acceptance scope (including [native Windows repairs and the expanded input audit](../testing.md#september-2026-failure-analysis), and [independent anchor-to-assembly authority](../architecture.md#source-coherent-locus-transcript-presentation)):
- Glen reports five-gene export QC at `d978b0c4` with LOCUS-column fix `d37bd427`; recheck the [selected annotated downloads](../integrated_locus_tss_profiles.md), receipts and raw signal/motif provenance at the final candidate. Before offering EMBL, re-export into fresh directories and independently compare GenBank/EMBL/FASTA sequences and annotations. [Synthetic format-parity coverage](../sequence_format_parity.md) does not replace his September 13 real-data acceptance.
- Run TP73/PATZ1 and promoter-matrix GUI smoke, including cancellation and [DNA-map/Structure source comparison](../transcript_source_presentation.md#dna-map-and-structure-tab) with unchanged primer targets; obtain canonical
  CD44/TGFB1/SERPINE1 region sets before real-candidate analysis; verify overlapping TSS bands, zero-position peaks and per-TSS CUT&RUN/TFBS arrays against continuation-page IDs in the exact supplied bundle, never infer missing ROIs.
- Repeat the copied-state IRF9/Q00978 acceptance and require the one-assay panel
  to report one of three linked cDNAs and two of four linked records covered,
  with only the true patch record genomically unassessed.
- Decide `.12` retention/removal of deprecated ClawBio shell-normalizer modes
  using the [compatibility inventory](../clawbio_gentle_integration_onepager.md#deprecation-table); `.11` eligibility did not authorize removal.

### Workflow Acceptance Details

1. Have Glen recheck three introductory GUI chapters plus promoter-similarity navigation and TATA evidence review/materialization at the merged SHA, retaining outcomes/screenshots. Run the 10-step learner path in [synthetic TSS acceptance](../tutorial/08-15_tss_collection_gui.md) first; retain Glen's `1fcee6c0` evidence separately from merged-candidate acceptance. Repeated-open reuse, Forget/Undo, restart and stale-member editing remain separate contributor checks, not part of that learner verdict. Verify compact validation, cascade and first-member native focus (including mixed open/new members) without harness refocusing; multi-window clutter and individual-member navigation remain follow-ups. Exercise the [TSS collection workflow](../tss_workspace.md#acceptance) on TP73: membership, strands/flanks, separate-process approval, cache refresh, DeltaNp73-only cropped locus, validated scans, repeated opening and inner-agent handoff. Compare symbol/ID memberships with mixed or case-variant gene labels; verify historical collections stay intact while new grouping requires fresh approval. Include stored-collection inspection and confirmed stale-registry removal without deleting sequences. Check unsupported-index warnings and locus-level unassigned transcripts remain distinct from loaded-feature coverage and gene-specific exclusions.
2. Carry historical `.10`/`.11` acceptance into one named `.12` candidate: locked tests, five native package entrypoints plus separate optional scripting checks, tutorial/parity, IRF9/Q00978, benchmarks, circular ORF/feature regressions and TP73/PATZ1 TFBS/reporter/CUT&RUN/locus proofs with exact JASPAR 2026 PFMs. Keep both tagged-source ledgers unchanged.
   Keep tutorial 08.13's live GUI capture and real Parquet acceptance separate from the [synthetic TSS tutorial replay](../tss_tfbs_profiles.md#rescoring-existing-reports) and mocked tests. For private five-gene rescoring, first obtain the old scored report/receipt and exact matrix/context bindings; the [historical input inventory](../tss_tfbs_profiles.md#check-input-completeness-before-regeneration) alone is not acceptance.
3. After explicit push/dispatch approval, run `.12` macOS/Linux/Windows CI (including all-feature checks), headless container validation and [build-only packaging](../release.md#build-only-candidate-verification) at one frozen SHA on a compatible pushed ref in `smoe/gentle_rs` (or an explicitly approved fork ref, not its stale `main`). Check the logged candidate SHA, not only the workflow revision; retain receipts and repository/run URLs before approval. Native GUI launch/tutorial navigation remain separate. Preserve the historical GUI-only `.10` Windows ZIP; `.12` must prove its five-binary packages separately. No further feature scope until these gates; tagging and upload remain owner-managed.
4. Have Glen refresh the `a5b01843` receipt with native interval/score and catalog-bound assembly checks; preserve verified figures unless mismatched. Validate the [candidate-first view and peak-centred follow-up](../cofactor_candidate_interface_plan.md); keep offline regeneration, live Conservation/partner-evidence proof and ClawBio compatibility pending removal review.
5. Have Glen run the [promoter-cofactor tutorial](../tutorial/08-14_promoter_cofactor_browser.md) against the exact candidate: live stale-form/copy binding, saved evidence, verified plus/minus feature preview/confirmation/undo, then the optional IRF9 package. Human scientific approval remains pending; unstranded evidence remains regions and raw coverage reconstruction is not supported by the reduced package.
   Validate DuckDB boundaries/limits, scan filters and selected-package GUI staleness; have Glen inspect [imported-hit triangles](../tss_tfbs_profiles.md#imported-duckdb-motif-hits) on both strands and run the [regulatory/TSS-subset reader](../regulatory_motif_subset.md) against the actual genome-wide delivery (including MT and source floors); synthetic fixture acceptance is not production acceptance. Native TSS follow-up (not an added release gate): live acceptance of [report attachment](../gui.md#attach-the-quantitative-report) and [view SVG export](../gui.md#export-the-native-tss-view), including CD44/TGFB1 gaps, local-strand legends, curve colours, span-only navigation, missing units, axes, selected lanes and sparse-hit coverage; reporter rows remain open. Validate approved [local scoring](../gui.md#compute-local-tss-scores): InlineSequence/matrix provenance, separate score kinds/scales/cache, cancellation, stale-form handling and genuine zero versus unavailable masks. Its agent/headless routes now exist (`ui open tss-view --local-score`, `promoters tss-view-svg`) with a real-locus walkthrough in [tutorial 08.17](../tutorial/08-17_tp73_dnp73_factor_curves.md); retain Glen's Linux/manual 08.16/08.17 screenshots at `8652ecf2` separately from merged-candidate equivalence, live hosted-intent adoption and Windows/macOS acceptance per [handoff](../tp73_dnp73_factor_curves_glen_handoff.md). Span-scoped headless local scoring is implemented pending tests/acceptance: complete footprints inside the requested span, full-window admission, explicit span counts/hashes and no hidden GUI/report rescoring. Native TSS export follow-ups are implemented pending tests/acceptance: disclosed default cap of 128 hover titles per lane/strand without dropping curve points or validity bands, top-three within-matrix raw positions, documented prepared-reference P1 access and validated collection SVG/index/hash receipts with per-window/per-lane scales (no new calibration). Retain current paint/file-size refusals; Glen measures the revised title counts and checks receipts, cancellation and GUI/headless equivalence before any limit increase. New live database queries still require exact reference-bound interval admission; retain [flat-file limits and readability plan](../tss_annotated_sequence_readability_plan.md). Glen: regenerate [source-coherent Ensembl/RefSeq pages](../transcript_source_presentation.md) from one enriched locus and verify mappings/PDF/SVG/HTML receipts; synthetic coverage is not real-data acceptance. Keep the [standalone comparison](../transcript_source_presentation.md#standalone-comparison-review) deferred rather than adding a second source-join implementation.
6. Continue the [responsive-command migration](../asynchronous_command_execution_plan.md#reconciled-codex-plan-2026-09-14): validate BLAST start-v2 deferred diagnostics, bounded FIFO/cancellation and observation-only status at the merged revision. Preserve DEC-026 merging and [DEC-048](../decisions.md#dec-048-command-submission-and-observation-never-wait-for-execution); audit general managed-command queuing/dependencies and result/file/job semantics before expanding beyond the BLAST queue.
   Have Glen validate genome preparation, Run/slash, live navigation/conversation and cancellation at one revision; extend the synthetic growing-journal cost probe to a representative project before claiming bounded UI lock latency. Keep cross-process control, durable admission and other adapters pending. No private T2T/model run in development; general asynchronous migration is not complete.
7. Carry into `.12` the [gene-informed primer-PAIR study GUI](../gene_assay_study_gui_plan.md), not single-primer capture. Have Glen annotate the authentic-reference [PATZ1 tutorial 04.08](../tutorial/04-08_gene_assay_study_gui.md) at G1-G6, including the 13-Ensembl/four-RefSeq source comparison and unresolved discrimination matrix; authentic whole-reference specificity and order approval remain pending. Supply an original typed paper-facing dossier/source bundle. Include coverage-scope interpretation and last-page navigation for large transcript/assay panels in that live check; synthetic pagination tests do not establish real-locus GUI performance. Retain the historical [reconciled `.11` scope](../internal_11_plan.md#review-reconciliation) at the new candidate: package parity, exact-SHA TP73/DeltaNp73/tutorial acceptance and independent tutorial/vector-PDF verification. Keep TSS managed-worker migration and capture-retention requirements optional; defer capture specificity/thermodynamics/order readiness until lab inputs are retained and calibrated cofactor ranking until anchor/expression/background contracts exist. Retain separate RefSeq accession matching and real-image gel acceptance; automatic peak detection remains deferred.
