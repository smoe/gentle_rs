# Claude Review of the TP73 CUT&RUN Promoter Comparison

Review the current, corrected TP73/CUT&RUN promoter-comparison implementation
and retained evidence in `/Users/u005069/GitHub/gentle_rs`. Work read-only with
Read, Grep and Glob only. Do not edit, build, run tests, use shell/network tools,
download private inputs, regenerate reports, change hashes, push or publish.
Give concise findings with file/line evidence and a revised minimal follow-up
plan. This is a review request, not authorization to implement that plan.

TP73 is the perturbation/occupancy context here; the compared TSS windows belong
to CD44, TGFB1 and SERPINE1, not to the TP73 locus itself.

## Revision and Integration Context

The owner asked to merge and commit Glen's older work before this review:

- Branch: `smoe-bot/glen/tp73-cutrun-promoter-comparison-20260909`.
- Branch tip: `b8b731586a31498a89ec14457d0cf20165de40ed`.
- Pre-merge main: `a02ef7779a226e3366ba3ac7d53a84e4d366a591`.
- Committed merge and review baseline:
  `4a8e81d39fa1f61ca32309969ffa50b2df7b55fc`.

This is an ancestry-only `ours` merge, not a restoration of the old analysis.
Its tree equals pre-merge main's tree. Glen's first commit, `54c8924a`, has
21 non-changelog outputs identical to the already-integrated `c21298fb`.
The different changelog wording in `c21298fb` deliberately qualifies the first
run as exploratory/uncorrected rather than complete. `b8b73158` is patch-equivalent
to the already-integrated ordering fix `dcc1a2b9`.

Do not review only the empty merge diff or propose reapplying these commits.
Review the actual files below, including their later corrections. Subsequent
work includes exact reference/report/SVG bindings, gene-ID exclusions, raw
target/HSP cap auditing, strand-aware layouts, chromosome-general BigWig
conversion, native interval/score verification and shared engine TSS geometry.
The retained real-data producer remains `ff44ebde`; the ancestry merge is not
a new scientific replay or acceptance verdict.

## Original Codex Plan and Outcomes

1. Inspect the branch and main before merging; preserve newer implementation,
   evidence and scientific qualifications instead of overwriting them.
2. Commit the integration without changing unrelated `paper` or `output/`,
   release gates, source code or retained report bytes. This is completed.
3. Check Python-only regressions and retained bundle integrity. Defer Rust/native
   execution and source-data replay under the owner's no-local-build restriction.
4. Prepare this scoped review prompt. Claude has not been invoked.

## Files to Read

Paths below are relative to the repository above.

- `scripts/prepare_tp73_cutrun_promoter_candidates.py` and
  `scripts/prepare_transcript_promoterome.py`: supported TSS selection and
  prepared promoterome/reference identity.
- `scripts/prepare_tss_regulatory_similarity_candidates.py`: selected transcript
  membership, assembly-bound feature intersections, assembly-forward DNA and
  delegation to `ComputeTssWindowGeometry` without a Python fallback.
- `crates/gentle-engine/src/tss_window_geometry.rs` and
  `crates/gentle-protocol/src/tss_window_geometry.rs`: shared geometry contract;
  root dispatch is in `src/engine/ops/operation_handlers.rs`.
- `scripts/compare_candidates_to_promoterome.py`: other-gene/self-locus exclusion,
  distinct gene/window accounting, HSP/query coverage and search-cap semantics.
- `scripts/tss_regulatory_report_binding.py`,
  `scripts/render_integrated_tss_regulatory_report.py`,
  `scripts/append_tss_similarity_to_locus_report.py` and
  `scripts/render_tp73_cutrun_promoter_comparison.py`: bound input resolution,
  recurrence rendering, block orientation, missing-data states and provenance.
- `scripts/bigwig_to_bedgraph.py` and
  `scripts/validate_tp73_locus_cutrun_lanes.py`: chromosome-general conversion,
  exact source/assembly bindings, clipping/filtering/score precision and
  sequence-anchor orientation rather than assumed gene orientation.
- `scripts/test_tp73_cutrun_promoter_comparison.py`,
  `scripts/test_tss_regulatory_integrated_report.py` and
  `scripts/test_tp73_locus_cutrun_lane_validation.py`: deterministic coverage.
- `docs/examples/regulatory_region_comparison/tp73_tss_local_integrated/`:
  README, selected/candidate JSON/FASTA, comparison, summaries, interpretation,
  validation/renderer receipts, SHA256SUMS and tall SVG/PDF/PNG derivatives.
- `docs/examples/regulatory_region_comparison/tp73_cutrun_supported_first_comparison/`:
  historical evidence, not a current completeness or functional verdict.
- `.gitattributes`, `scripts/test_tutorial_checkouts.py`, `docs/testing.md`,
  `docs/architecture.md`, `docs/decisions.md` DEC-044 through DEC-047,
  `docs/roadmap.md` and relevant September 9-13 entries in `docs/CHANGELOG.md`.

## Evidence and Known Limits

At the exact merged baseline, macOS arm64 / Python 3.14.7 passed 29 checks:

```sh
PYTHONDONTWRITEBYTECODE=1 python3 -B -m unittest \
  scripts.test_tp73_cutrun_promoter_comparison \
  scripts.test_tss_regulatory_integrated_report.FrequencyTests \
  scripts.test_tss_regulatory_integrated_report.RetainedBundleTests
```

These cover synthetic selection/comparison, both-axis SVG rendering, block
ordering and current retained manifest/derivative integrity. They do not replay
the private locus reports, BigWigs or the roughly 109 MiB raw comparison tables.
The seven BigWig checks could not import `pyBigWig` in this environment; this is
an unavailable dependency, not a biological failure. The SourceBindingTests and
TssWindowTests groups (14 checks) were deliberately not run; executing both
groups in full requires a current real CLI.
No Cargo build/check/test, native Windows/Linux or live GUI acceptance is claimed.

Codex reproduced one portability issue without changing the checkout:
`git check-attr text eol` reports no policy for the retained candidate JSON.
Its LF SHA-256 is
`dd467f4d33d463b11e95f697b69361169a74e14cfce6cb2db1b06deee7ff129a`.
`git -c core.autocrlf=true cat-file --filters` on the same tracked blob yields
`f7a9be2cba6915af0eb6efc466cd815ec855d7eec72765752897f9f4d7f0905d`.
That disagrees with the retained manifest and raw-byte receipt. The existing
fast TSS checkout test covers GUI screenshot evidence, not this analysis bundle.
Assess a scoped LF policy and LF/CRLF regression with a missing-policy negative
control; do not suggest normalizing inputs or refreshing scientific hashes.

The bundle README explicitly says that the retained `ff44ebde` lane receipt
predates stronger catalog-assembly and complete native interval/score checks.
Refreshing it requires the original request/report JSON and source BigWigs.
The existing 36-lane / 11,797-interval receipt is historical: a valid hash or
nonempty plotted lane cannot certify these newer checks. Do not relabel it or
request redrawing figures unless exact replay identifies a mismatch.

## Review Questions and Requested Response

Check for remaining concrete bugs in source identity, assembly/strand/coordinate
conversion, recurrence/counting/cap semantics, bound rendering and evidence
claims. Distinguish annotations, predicted TFBS, occupancy/enrichment and sequence
recurrence from causal regulation or reporter sufficiency. In particular, a
missing hit, exhausted cap, motif-scale non-search or unavailable signal must
not become evidence of uniqueness, completeness, activity or a measured zero.

Separate already-repaired historical defects from current code defects and
acceptance gaps. Does the evidence justify only structural reporter hypotheses,
and are the data/model/plot boundaries compatible with the shared-engine rules?
Keep any architectural follow-up minimal; do not propose a new parallel GUI,
publication model or wholesale Rust port merely because historical scripts exist.

Return risk-ordered findings with exact current paths/lines, a brief no-finding
note for checked boundaries, and a revised minimal plan with deterministic
acceptance checks. Route real source-data refresh to Glen at one frozen SHA;
keep native/platform, scientific and `.12` acceptance separate. Say which
questions cannot be settled from retained files instead of inventing a pass.
