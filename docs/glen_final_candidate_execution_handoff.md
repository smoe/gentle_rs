# Glen: Consolidated Final-Candidate Checks

Coordination update, 2026-10-13. The relayed `gentle_rs_2_main` handoff is
integrated into this request. Its integration baseline is
`5e6201a69a599cba87491c86cfd7311cb6fee173`, shared by local `main`,
`gentle_rs_2_main` and last-fetched `origin/main`. This is not selection or
acceptance of the final candidate.

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
