# Glen: Scrutinize the Remaining Tutorial Slices

Prepared 2026-10-12 from local `150a339074e7a2d282567d791a8f50834b70ff9c`,
after Steffen supplied Claude's read-only review and authorized implementation.
Navigation, microarrays and measured-gel teaching are already separated in
`d95286e31caa96a59e711c7504347c0a21009b58`; their
[specific review request](glen_microarray_gel_tutorial_review_20261012.md) remains
valid. This second request covers A/B and the missing real inputs for E.
Audit the exact commit containing this note, not the old preparation baseline.
Steffen can forward it; no agent-to-agent message has been sent.

## Common Instructions

Use a clean isolated worktree and one exact-source CLI/GUI/docs-helper set.
Leave private studies and the primary working copy untouched. Report full
source SHA, binary self-identities/hashes, Cargo.lock hash, toolchain, effective
profile/features, OS/display environment and drift before evaluating anything.
Do not combine old binaries with a new source verdict. No push, tag, publication
or replacement of historical evidence is requested.

Read the learner prose as a wet-lab biologist, then execute the exact shared
commands. Separate product defects, missing/incompatible data, environment
problems and human-scientific-review-needed. A missing tool is not a biological
zero, and source checks are not a live GUI pass. Preserve raw downloaded and
normalized sequence digests as distinct provenance roles. Never assign an old
checksum to a newly retrieved substitute.

Local preparation is deliberately limited: Python source/fixture checks and
real DuckDB 1.5.4 package generation pass on macOS (84 scoped tests, including
two real I/O tests). Formatting/whitespace and catalog/manifest freshness pass;
the latter use the existing unchanged-source `150a3390` docs helper, not a new
candidate binary. Comparing the unchanged
default generator with `150a3390` gives 15 equal Parquet tables, identical
Parquet bytes and five equal catalog tables; this does not claim identical
DuckDB catalog bytes across runs/versions. New Rust regressions, GENtle query/
export commands and native screenshots have **not** run locally, per the owner.
No scientific, release, timing or native-Windows acceptance follows.

## A: Nonempty Imported Hits in Tutorial 08.16

Read the optional companion in
`docs/tutorial/08-16_tss_regulatory_view_gui.md` and
`test_files/fixtures/regulatory_motif_subset/README.md`.
The existing generator's `--example tss-view` preset matches the retained
08.16 report exactly: reference/assembly/contig, MA0004.1 six-base footprint,
plus `[96,107)` and minus `[292,303)` BED windows. Its five scores are explicitly
hand-crafted `log2_relative_risk`, not calculated PFM scores or experimental
occupancy. Audit/atlas placeholder fields attest no executed provider audit.

First run the fixture check and actual shared-query/export regression:

```sh
python3 -m unittest scripts.test_regulatory_motif_fixture
cargo test --locked --lib genomic_motif_evidence::regulatory::tests::regulatory_tss_view_companion_query_export_preserves_geometry -- --ignored --exact
```

Verify **one test actually runs**, not a zero-test pass. Set
`GENTLE_DUCKDB_BIN` to your real executable if needed. This test is ignored by
default because it needs Python and actual DuckDB, not mock Parquets.

Execute every companion command into a fresh `/tmp/tss-import-lesson` parent;
the package/export subdirectories must be new. Retain actual full, partial and
plus query JSON, fresh export receipt, enriched `report.json` and imported-hit
table. Never manufacture these outputs by editing the old report. Both minus
queries should return three hits; the plus query two. The narrower query
`[294,301)` does not cover the full displayed minus window, even when the hit
counts match. Verify actual coverage and raw/canonical source hashes.

Expected minus-axis projection:

| Genomic footprint / strand | Local full footprint / strand | Clipped display |
| --- | --- | --- |
| `[294,300)` / `-` | `[3,9)` / `+` | unchanged |
| `[296,302)` / `+` | `[1,7)` / `-` | unchanged |
| `[290,296)` / `+` | `[7,13)` / `-` | `[7,11)` |

Check genuine zero, negative raw values and separate per-source/matrix scales.
Imported `log2_relative_risk` must never be relabelled as the bound report's
`llr_bits`. Confirm attachment leaves original score arrays,
context, annotated input and old receipts/screenshots unchanged.

For native acceptance, open the original annotated minus record, attach the
fresh enriched report through the GUI-hosted Agent Assistant route (or loader),
and capture complete/partial imported lanes with readable hover/side summary,
strand inversion, clipped footprints and interval-only navigation. Check that
a report-only hit does not require an annotation feature ID. Preserve the
unavailable control; no new control biology is invented. 08.17 and 08.05/08.09
teach local scoring and reporter continuation separately.

The retained old screenshots do not prove this companion. Return raw native
images, semantic snapshots where available and a hash manifest; any annotated
teaching crop is a separate derivative. Do not replace historical captures.

## B: Repeat Environment Cohorts in Tutorial 08.18

Read `docs/tutorial/08-18_repeat_environment_cohorts.md` and
`test_files/fixtures/repeat_cohort/README.md`.
The local reference is genuine retained bacterial AB011549.2. The three rmsk
rows and LINE/L1/SINE/Alu labels are **invented software inputs**, not mammalian
repeat evidence. Do not name the bacterial output a human LINE cohort.

```sh
cargo test --locked --lib repeat_cohort_local_project_shell_pipeline_keeps_distinct_loci
```

Run the learner commands in one fresh state and isolated cache. The alias should
select two of three rows; forward-genomic substrings at `[1843,1849)` and
`[3877,3883)` are both `AGCCAG` with opposite repeat strands. Keep both locations
as separate rows. Require two **actual** midpoint windows, nonempty score arrays
and genomic peak positions; zero windows with warnings is not success.
The test uses paths containing spaces/apostrophes to exercise shared quoting.

Rebuild all four geometry modes into separate files. Selection must remain two
loci; this GenBank reference lacks the prepared transcript-index projection,
so UTR, promoter-upstream and stop modes should have explicit reasons and no
invented anchors. Unavailable projection is not absence of source CDS.

If compatible genuine public genome/GTF/rmsk resources are already prepared,
optionally replay the human example and retain their source/assembly/digests.
Show available plus/minus transcript contexts and a missing-CDS case if the
inputs contain them. Do not download a full reference just for a screenshot.
Keep source labels and truncation visible; nearest-gene assignment is not proof
of regulatory targeting. RNA ranking is **unimplemented**: unset rank fields and
empty support lists are not measured zeros. Do not invent `--rna-report` or
an external ranking algorithm. No dedicated cohort table/row action is claimed.

Capture only actually supported GUI follow-ups on a genuine anchored locus.
Retain CLI reports when there is no dedicated cohort GUI. Separate software
contract acceptance from real repeat biology and RNA interpretation.

## E: Supply Retained Genuine FLNA Inputs Before Activation

`docs/tutorial/drafts/flna_splicing_interpretation.md` is not in the catalog.
The historical temporary FLNA state/view are gone; numeric feature `300` is
not a stable tutorial selector. Please locate original retained public evidence
and verify hashes/redistribution permission, or retrieve fresh public inputs
with fresh provenance. Do not fabricate or rename a synthetic substitute FLNA.

Prefer release-pinned archive/GTF data. Historical gene ENSG00000196924,
canonical ENST00000369850.10, reviewed UniProt P21333 and 36 transcripts are
context, not a timeless expected count/version. If using live REST, retain raw
provider responses and `/info/data` before/after; that alone does not pin each
response to a historical release.

The shared acquisition routes in a fresh state are:

```sh
gentle_cli --state /fresh/flna/state.json ensembl-gene fetch FLNA --species homo_sapiens --assembly GRCh38 --flank-bp 2000 --entry-id flna
gentle_cli --state /fresh/flna/state.json ensembl-gene import-sequence flna --output-id flna_locus
gentle_cli --state /fresh/flna/state.json shell 'uniprot fetch FLNA_HUMAN --entry-id flna_uniprot'
```

Create the parent first. Discover the desired transcript's current feature ID
from its stable annotation identity in the retained bundle. Resolve that ID
before using the existing shared commands (placeholders below are not IDs):

```text
inspect-feature-expert flna_locus splicing FEATURE_ID
render-feature-expert-svg flna_locus splicing FEATURE_ID /fresh/flna/splicing.svg
ui open splicing-expert flna_locus FEATURE_ID
```

The first two run through `gentle_cli --state ... shell '...'`; the last requires
GUI-hosted adoption. Retain the offline-reopenable state, source inventory,
release/accession versions, raw and sequence hashes, exact loaded-xref evidence,
view fingerprint and actual command outputs. Show any reduction preserves exon,
CDS and transcript geometry plus exact UniProt joins.

Choose one scientifically informative exon selection and record the observed
frame/protein consequence and uncertainty, not an invented expected product.
Capture readable native overview/detail plus a separately rendered 2400-px SVG
image. Explain boundary motifs, both canonical flags, exceptional markers,
loaded-evidence `?` diagnostics, extra-entry counts and unknown joins. Local
plus and genomic minus must remain distinct. Unlinked rows are not global
UniProt absence and a projected exon skip is not an observed splice event.
No active FLNA learner guide or human sign-off is requested before these inputs.

## Verification and Return

```sh
python3 -m unittest scripts.test_tutorial_gap_lessons scripts.test_regulatory_motif_fixture scripts.test_tutorial_checkouts scripts.test_tutorial_acceptance
gentle_examples_docs tutorial-catalog-check
gentle_examples_docs tutorial-manifest-check
gentle_examples_docs tutorial-check
cargo check -q --locked
cargo fmt --check
git diff --check
```

Confirm the catalog has 65 entries and the executable manifest still has 30
chapters/six GUI contracts. 08.18 is manual/hybrid, not a fabricated executed
contract. Historical human review dates and all old scientific bundles remain
unchanged. Missing real DuckDB tests must remain skipped/not-run explicitly.

Reply **through Steffen** with findings by severity, exact SHA, per-A/B/E
pass/fail/not-run/blocked, commands, hashes and retained artifact paths. Separate
input preparation, Rust/shared-command execution, native interaction,
screenshot suitability and scientific interpretation. Say which platforms were
not tested. Return a branch/full SHA only if separately authorized; otherwise
return the checksum-bound evidence and minimal reproductions without pushing.
