# Glen: Consolidated Final-Candidate Checks

## Current Status: Do Not Reissue The Full Request

Reconciled with Glen's 2026-10-13 clarification and the retained
[original execution report](audits/glen_followups_20261008/candidate_f35a58e6/REPORT.md).
The software/learner part of the historical request below was executed at
`f35a58e6d2bdd7e0a5551e4531c336d7949bfaf0`, not the earlier `5e6201a6`.
Preserve that exact-source verdict; do not relay the old worktree's full list
as a new assignment.

Already evidenced at `f35a58e6`: one actual 08.16 DuckDB regression and 3/3/2
learner hits with preserved geometry/scores; the 08.18 LocalProject pipeline
and distinct opposite-strand `AGCCAG` loci; 07.04 microarray, 26 gel tests and
separate-process CLI save/reopen; FLNA source/hash/UniProt joins and the
2400x14808 rendering; focused 04.09/parity and complete real, hash-bound BLAST
commands. The 60-second real imports remain performance-blocked, with no
completed result/state and no CPU-stack attribution. See the
[separate bounded profiling receipt](audits/glen_followups_20261008/post_blast_profile_f35a58e6/README.md).

Do not rerun 08.18 CLI acceptance or invent a cohort GUI; none exists. RNA
ranking stays unimplemented. Tutorial-CI/review-registration/LF-CRLF and link
repairs remain closed by retained `5383f7ec` evidence. A newly observed defect
is a new regression, not erasure of those passes.

## Remaining Native And Human Checks

These four residual tutorial checks are the active consolidated request, not
the historical command list below. Use only an explicitly named source and
its same-SHA GUI binary on an approved native host, retaining input/binary
hashes and platform/profile. Closing a residual at `f35a58e6` does not accept a
newer main; testing a newer source does not rewrite its historical verdict.
Missing host, capture consent or human decision remains blocked/not-run.

| Area | Remaining Check Only |
| --- | --- |
| 08.16 | Open the retained complete/partial imported lanes in the native view, including plus/minus, clipped footprints, negative/real-zero values, units/scales and coverage warnings. Check normal-size label/legend readability. Retain approved original captures; do not rerun DuckDB queries or alter original score arrays to get a screenshot. |
| 07.04 / 10.03 | Native navigation, readable configuration/result views and approved screenshots; save/reopen a disposable copy through the GUI. The existing CLI reopen evidence is already passed. Keep assembly/mapping and synthetic-statistics warnings visible; outside-calibration gel values remain null with `outside_calibrated_range`. |
| FLNA | Inspect the normal-size native transcript/exon/UniProt presentation and both genomic-minus/local-plus orientation claims. Obtain a human teaching-exon choice with its intended lesson, consequence/uncertainty and reviewer; Codex/Glen must not manufacture approval. Seek historical annotation-release evidence only if a genuine contemporaneous record is available. The 2400-px render/source joins are already evidenced; keep the draft inactive. |
| 04.09 | Inspect the existing native Shell/report presentation, without claiming a new specialist GUI or a completed real-import verdict. Complete real-reference import/save/reopen acceptance only after an evidenced performance repair, at that repair's exact SHA, with unchanged raw HSPs, policies and reference bindings. A timed-out attempt is not a failed biological verdict. |

Do not enable the disabled screenshot bridge, bypass isolation/host policy,
reuse render-only images as native captures or substitute today's Ensembl
release for missing historical provenance. Preserve original inputs/evidence;
keep annotated teaching crops separate from captures.

## Separate New-Code Diagnostics

Local runtime commit `52908b2fd0a24489156413c48e5edd9eb5df8a56` contains the
merged panel matrix, its MCP path-test repair and opt-in import phase tracing.
It is **not** the executed `f35a58e6` candidate or a measured performance fix.
Only its focused matrix/trace/import regressions and bounded diagnostic pairs
are requested by the [new phase-trace procedure](glen_post_blast_import_phase_trace.md).
Those timed diagnostic attempts do not request a full real-import acceptance
before repair, or repeat the completed tutorial suites. No Rust execution of
that new source is claimed by Codex, and no push is authorized.

Return each residual's pass/fail/not-run/blocked with exact source, binary/input
hashes, host/profile, observations and retained capture/state/log locations.
Human decisions must name the actual reviewer and remaining uncertainty.
The broader six-item ledger remains the status owner: 08.04 isolation, live
agent, native latency and `package-opt1`/extracted-package gates stay separate
and open; this clarification neither dispatches them nor grants publication.

## Historical Full Request: Reference Only

The original coordination text and command list follow unchanged for audit
context. Completed `f35a58e6` checks above must not be redispatched from it.

Coordination update, 2026-10-13. The relayed `gentle_rs_2_main` handoff is
integrated into this request, alongside the merged visual/provenance follow-ups.
The pre-integration baseline was
`5e6201a69a599cba87491c86cfd7311cb6fee173`; do not substitute it for the full
merged SHA in Steffen's forwarded message. No candidate acceptance is implied.

## Freeze Before Execution

Steffen will name the **full final pushed SHA** in the forwarded message, after
all intended integrations. If it is missing or unavailable, ask for it; do not
substitute moving `main`, a tag, an unmerged branch or nearby binaries. Read this
document and its linked procedures at that same SHA, then use a clean, isolated
detached checkout. Leave Glen's private/dirty primary checkout untouched.

Record HEAD, tree, clean status, raw `Cargo.lock` SHA-256, toolchain, platform,
features/profile, executable self-identities and binary hashes before execution.
Use only binaries built from that SHA and keep caches/new evidence outside the
checkout. Rust execution belongs on Glen's approved external hosts; Codex ran
no local Rust build/check/test for this handoff. No Claude consultation is needed.

## Current Tutorial Checks

1. **08.16: imported TSS evidence.** Follow
   [the companion procedure and geometry oracle](glen_tutorial_gap_preparation_20261012.md#a-nonempty-imported-hits-in-tutorial-0816)
   and the candidate's [learner commands](tutorial/08-16_tss_regulatory_view_gui.md).

   ```sh
   python3 -m unittest scripts.test_regulatory_motif_fixture
   cargo test --locked --lib genomic_motif_evidence::regulatory::tests::regulatory_tss_view_companion_query_export_preserves_geometry -- --ignored --exact
   ```

   Require the real DuckDB executable (`GENTLE_DUCKDB_BIN` if necessary) and
   exactly one executed Rust regression, not a zero-test pass. Execute the
   full/partial minus and plus queries into fresh directories (3/3/2 hits).
   Inspect native complete/partial lanes, strand inversion, clipped footprints,
   signed and genuine-zero scores, per-source units/scales and coverage warnings.
   Attaching imported evidence must leave original score arrays, context,
   annotated inputs and historical reports/receipts/captures unchanged. The older
   `a403d2eb` A1 package is not this newer companion's acceptance.

2. **08.18: repeat cohorts.** Follow
   [the LocalProject procedure](glen_tutorial_gap_preparation_20261012.md#b-repeat-environment-cohorts-in-tutorial-0818)
   and [all learner commands](tutorial/08-18_repeat_environment_cohorts.md).

   ```sh
   cargo test --locked --lib repeat_cohort_local_project_shell_pipeline_keeps_distinct_loci
   ```

   Require nonzero test execution, two actual midpoint windows and retained
   identical DNA at distinct loci (`AGCCAG`, opposite repeat strands), with
   nonempty scores and genomic peak positions. Exercise all four geometry modes;
   missing transcript projection stays explicitly unavailable, not a fabricated
   anchor or evidence of absent CDS. The bacterial sequence is genuine but the
   repeat labels are synthetic. RNA ranking is **not implemented**; retain unset
   ranks/empty support, not measured zeros. Do not invent a cohort GUI or RNA route.

3. **07.04 microarrays and 10.03 gels.** Follow
   [the native/learner acceptance procedure](glen_microarray_gel_tutorial_review_20261012.md)
   using [07.04](tutorial/07-04_microarray_evidence_gui_cli.md) and
   [10.03](tutorial/10-03_measured_gel_sizing_gui_cli.md) at the candidate SHA;
   earlier preparation baselines in the procedure are historical, not this target.

   ```sh
   cargo test --locked --lib project_microarray_track_uses_vendor_subset_on_tp73_genbank_anchor
   python3 docs/tutorial/inputs/make_gel_sizing_fixture.py docs/tutorial/inputs/gel_sizing --check
   cargo test --locked --lib gel_image
   ```

   Execute learner imports, inspection/projection, analyses and exports. Capture
   readable native configuration/result views, navigation and save/reopen.
   Preserve assembly/mapping warnings and synthetic-statistics labels. Both
   out-of-calibration gel outcomes remain null with `outside_calibrated_range`;
   no extrapolation, zero substitution or experimental-calibration claim.

4. **FLNA: remaining evidence only.** Genuine archived sources, hashes, stable
   selector and transcript/UniProt joins are already
   [independently checked](glen_tutorial_inputs_review_20261012.md#a2-genuine-retained-flna-snapshot).
   Seek the missing historical annotation-release record in retained acquisition
   logs; do not label the old snapshot with today's release. If unavailable,
   report that gap and propose a separately reviewed dated-snapshot or newly
   pinned-input route. Obtain human teaching-exon selection and record the
   observed consequence/uncertainty before asserting an expected outcome.
   Perform the [remaining normal-size native and 2400-px rendered-view checks](glen_tutorial_gap_preparation_20261012.md#e-supply-retained-genuine-flna-inputs-before-activation),
   keeping genomic minus/local plus and unknown joins explicit. Keep the
   [draft](tutorial/drafts/flna_splicing_interpretation.md) inactive; source
   integrity or a rendering pass does not supply scientific approval.

5. **04.09: reference-bound primer specificity.** Execute every focused Rust,
   parity and Python command in the candidate's
   [Verification Scope](tutorial/04-09_reference_bound_primer_specificity.md#verification-scope),
   with nonzero test counts and explicit platform/tool skips. Inspect the
   supported native Shell/report presentation; no new specialist GUI is claimed.
   Separately execute the learner handoff/external-runner/import path against
   actual prepared references, retaining commands, raw TSVs, execution manifest,
   reference/annotation bindings and independent genomic/transcriptome verdicts.
   Test required failure, missing/stale/incomplete evidence, retained history and
   save/reopen. A positive Claude source review or Unix fake-tool fixture is not
   real BLAST, native Windows, scientific specificity or order acceptance.

## Integrated Visual And Provenance Checks

Read the [nineteen visual-placement decisions](tutorial/visual_forward_port_review.md)
and run `python3 -m unittest scripts.test_tutorial_checkouts` at the merged SHA.
Verify catalog/manifest and generated-chapter freshness with the exact-source
docs binary, then inspect the added graphics in their learner/native context.
Twelve placements belong in the proposed chapters; the five dotplot captures
stay in 02.03, and the two unsuitable agent/variant illustrations stay omitted.
Historical screenshots, teaching cartoons and saved-state tester examples must
not become fresh GUI, fabrication, assay-specificity or scientific evidence.

Check the minimal [E-MTAB-15709 README clarification](../data/resources/cutandrun_20250602_noDuplicates/README.md)
against its pinned public SDRF: all three named `pos_*_SKMel29` libraries are
H3K4me3, not HDAC2 or an unspecified target. Keep this source binding separate
from broader DN terminology or claims about unbound local filenames. Do not
import scientific report branches, RefSeq alias hardcodes or dirty-worktree data.

## Six Gates And Historical Separation

Keep the [six-item audit ledger](glen_audit_followups_20261008.md) as the status
owner: audit preservation, 08.04 native/agent contracts, genuine live agent,
honest visual evidence, 24 scroll observations and `.12` package acceptance.
Use the [six-gate execution procedures](glen_e981469e_execution_handoff.md#six-follow-up-gates-e981469e)
for fresh candidate checks, **not their historical SHA/bundle defaults**.
For any separately authorized non-publishing package run, `candidate_sha` must
be Steffen's final SHA, `publish=false`, profile `package-opt1`; retain actual
workflow revision, extracted-package and container execution receipts.

The explicit frozen `e981469e` and `dfadf2a5` requests in that older handoff
remain separate audits at their original SHAs. Attribute `a403d2eb` retained
inputs/execution separately as well. None automatically accepts the final
candidate. Keep experimental sequence design outside the `.12` release gate.

Tutorial-CI/review-registration/LF-CRLF and link repairs are **closed** by
retained evidence at
`5383f7eceeb761d778b7d26efc7d0f5c4486d941`; do not reopen them merely because
older handoffs list them. Run ordinary candidate freshness/regression checks,
including exact-source `tutorial-check` and
`python3 scripts/check_tutorial_checkouts.py --binary /EXACT/gentle_examples_docs`.
Report a new concrete failure as a new candidate regression, not erasure of
that historical pass. Skipped native/agent/scroll/package gates stay open.

## Return Through Steffen

Return per-item **pass/fail/not-run/blocked**, with source/tree/lockfile,
binary/input hashes, exact commands/exit codes/test counts, platform/profile,
tool versions and retained artifact locations. Include complete logs, raw native
captures, states/exports/receipts and a per-file SHA-256 manifest plus archive
hash/size. Distinguish product defects, environment/data limitations and missing
human decisions. State untested platforms and real-reference/native skips.

Capture only under existing approved host/project policy; do not enable the
disabled screenshot bridge or bypass consent. Keep annotated teaching crops
separate from original images. Preserve `paper`, `output/`, private studies and
historical evidence. No push, workflow dispatch, tag or release/image publication
is authorized by this request; those require separate owner authorization.
Return failures before source repairs, which need their own commit and reruns.
No scientific approval is implied, and no credentials belong in the evidence.
