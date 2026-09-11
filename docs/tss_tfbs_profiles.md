# Accession-Pinned TSS Profiles

This development workflow plots predicted TF-binding scores over already
transcript-oriented TSS windows. It does not infer TSSs or retrieve a genome.
It leaves the older locus-evidence figures and their defaults unchanged.
Native wizard support and exact-candidate release acceptance remain separate
work; this is not a new `.10` scientific release requirement.

## Offline Example

From the checkout root, using a built `gentle_cli`:

```sh
gentle_cli features tss-tfbs-profiles \
  --manifest test_files/fixtures/tss_profiles/manifest.json \
  --panel test_files/fixtures/tss_profiles/panel.json \
  --selection selection.json \
  --expected-genome-id synthetic-genome-v1 \
  --expected-assembly synthetic-assembly-v1 \
  --expected-annotation-release synthetic-annotation-v1 \
  --output-dir tss-example \
  --formats svg,png,pdf
```

The output parent must already exist; choose a fresh destination. This tiny
plus/minus example tests software behavior, not promoter biology. Its 11-base
windows deliberately include an ambiguous base and are not adequate biological
promoter windows. The pinned Arnt matrix comes from the bundled JASPAR registry;
the synthetic DNA does not come from an organism. See the fixture
[provenance and header specification](../test_files/fixtures/tss_profiles/README.md).

To change presentation without rescoring:

```sh
gentle_cli features tss-tfbs-profiles-export \
  --report tss-example/report.json \
  --output-dir tss-reexport \
  --formats svg --panels-per-page 1
```

The same `features ...` commands work in the GUI Shell. JSON operation/workflow
adapters use `ComputeTssTfbsProfiles` and `ExportTssTfbsProfiles`; MCP exposes
these through `op` with explicit `confirm: true`. Neither operation changes
source sequences. Both are conservatively classified as external effects,
because they may publish files. Input availability and scientific compatibility
are validated at execution, not inferred from an empty project's readiness.

## CUT&RUN, Gene Structure And TATA Context

Add `--context-manifest FILE` to either command above to put context on each
matching gene's **detailed TSS pages**, immediately above its TF scores. This
is optional: old reports and commands retain their score-only layout. A context
entry applies to all the gene's TSS windows, including unselected ones; selected
pages can still be assembled into the integrated per-gene PDF.

The shared engine reads a hash-bound locus report, its full source-locus FASTA,
and an optional TATA screen. It does not rerun TF scoring, call peaks, retrieve
data or change the project. The details show:

- CUT&RUN/chromatin lanes and their sample/control labels, retaining the source
  locus's scale and raw interval scores. Cropping does not renormalize a signal.
  Values beyond a fixed source scale are clipped visually but retained in JSON.
  Gaps and missing controls are labelled, never filled with measured zeros.
- Exons for the exact TSS transcript members, with biological exon ordinals,
  thicker CDS intervals, and green translation-start / red translation-stop
  markers when the source annotation supplies them. Other transcript models
  remain in the full locus overview. Markers outside this window are not moved
  into it; unknown CDS geometry is not reconstructed by finding an arbitrary ATG.
  Source ranges classified as inferred ORFs or lacking a coding classification
  are explicitly warned about, not displayed as annotated CDS.
- TATA **source annotations**, **TBP predictions**, and **EPD classifications**
  as separate evidence rows. An EPD classification marks its TSS, not an inferred
  exact TATA interval. No supplied screen, or no reported row in this window,
  does not establish TATA absence.

All tracks use the same transcript-oriented 5'-to-3' axis as the TF scores,
including decreasing genomic coordinates for negative-strand genes. The context
and its provenance repeat on continuation pages. Large contexts fail the existing
readable-page limits explicitly instead of silently dropping lanes.

For a fully offline demonstration, after the first example above:

```sh
gentle_cli features tss-tfbs-profiles-export \
  --report tss-example/report.json \
  --context-manifest test_files/fixtures/tss_profiles/context/manifest.json \
  --output-dir tss-context-example --formats svg,png,pdf
```

Only the synthetic minus-strand gene gets context in this example. Its 11-base
window, exon/CDS markers, signal and three TATA evidence types are software
fixtures, not biological findings. See the [context fixture provenance](../test_files/fixtures/tss_profiles/context/README.md).

For real data, prepare a manifest with this structure; replace every placeholder
with the exact values from the existing TSS/locus reports and input-file hashes:

```json
{
  "schema": "gentle.tss_detail_context_inputs.v1",
  "reference": {
    "genome_id": "EXACT_GENOME_ID",
    "assembly": "EXACT_ASSEMBLY",
    "annotation_release": null
  },
  "genes": [{
    "gene_id": "EXACT_GENE_ID_FROM_TSS_REPORT",
    "locus_report": {"path": "gene.locus.json", "sha256": "FILE_SHA256"},
    "locus_fasta": {"path": "gene.locus.fa", "sha256": "FILE_SHA256"},
    "tata_report": {"path": "gene.tata.json", "sha256": "FILE_SHA256"}
  }]
}
```

Paths resolve relative to the manifest (explicit absolute paths also work).
Hash exact file bytes with `shasum -a 256` or Linux `sha256sum`. The full locus
FASTA must contain one uppercase IUPAC DNA record in the **loaded sequence's**
orientation, matching `sequence_binding` in the locus report. This is not the
short TSS FASTA. Export it from the same loaded DNA sequence that produced the
locus report; do not independently retrieve a potentially different version.
The engine slices that sequence and compares the oriented bases with every
matching TSS record's SHA-256 before attaching any context. It also checks gene
symbol, explicit gene ID selection, chromosome, anchor strand, interval and
reference labels. Exact versioned transcript IDs are used, without alias guessing.

The manifest `reference` must equal the TSS report's reference, including a
declared annotation release. If a release is declared, the locus must agree.
If it is null, no release is invented from a genome ID. For compatibility with
existing locus composers, `isoform_evidence.assembly` may equal either the exact
TSS assembly or its exact genome ID. Leading `chr` is the only chromosome-label
normalization. Hash/geometry agreement is consistency evidence, not independent
authentication of a reference or an experimental sample.

Generate a TATA report using the existing [TATA screen](tata_box_evidence.md),
on that same loaded sequence, for example in the shared Shell:

```text
op '{"SaveFile":{"seq_id":"EXACT_LOCUS_SEQ_ID","path":"gene.locus.fa","format":"Fasta"}}'
promoters tata-screen '{"seq_id":"EXACT_LOCUS_SEQ_ID"}' --output gene.tata.json
```

Review the screen's TSS restrictions and thresholds; supply explicit TSSs or
opt into an unrestricted scan only when intended. Omit `tata_report` when no
screen is available. Do not replace that state with an empty synthetic report.
The context resolver checks the TATA report's own content hash, sequence hash
and complete genome anchor. Its source warnings and classification remain visible.

The enriched `report.json` stores `windows[].detail_context`
(`gentle.tss_detail_context.v1`), original genomic intervals, clipped local
intervals, source metadata and file digests. The export receipt includes these
context inputs. Subsequent export of that report needs **no context manifest or
source files**; all original TF scores and comparisons are retained. To change
context, enrich the original score-only report into a fresh directory. An
already attached context is never silently replaced. No approval or checkpoint
binding is weakened.

For Glen's existing reports: use the same locus JSON that supplies the overview,
including its CUT&RUN groups; attach it with the matching source FASTA and optional
TATA report, then give the new report/index/receipt paths to the
[integrated PDF composer](integrated_locus_tss_profiles.md). The composer rejects
an enriched detail page whose locus-report hash differs from the overview's.
This does not rescore the expensive TF tracks or rewrite the existing export.

## CLI Output

After a successful export, CLI stdout contains
`tss_tfbs_profile_summary` with schema `gentle.tss_tfbs_profile_cli_summary.v1`,
instead of `tss_tfbs_profiles` and the complete `tss_tfbs_profile_receipt`.
The normal outer `result` wrapper (for shell routes), operation ID, messages and
warnings are retained. The summary links `report.json`, `index.json` and
`receipt.json`, with the report/index hashes, TSS/page/file counts, source and
producer/exporter revisions, report warnings and verification limits. It does
not contain score vectors, matrix PFMs or the per-page font audit.

The **complete** scored report and receipt remain in the export directory.
Existing scripts reading full JSON from stdout must opt in with the global flag:

```sh
gentle_cli --full-report features tss-tfbs-profiles-export \
  --report tss-example/report.json --output-dir tss-full-output --formats svg
```

Place `--full-report` before `features`, `shell`, `op` or `workflow`, not inside
a shared-shell command or JSON request. Default summaries apply to direct TSS
commands, CLI `shell` TSS/`op` commands, and exported TSS operations in the CLI's
direct JSON `op`/`workflow` entrypoints. Unexported compute results remain full:
there is no saved copy to point to. GUI Shell, MCP and other adapters retain
complete typed results; unrelated CLI results are unchanged.

Computation summaries copy the report's `verification` verbatim, including
`prepared_reference_not_assessed`. Report-only export summaries say
`not_reassessed_report_only_export`; inspect the saved report for its original
verification/reference/warnings. Neither mode upgrades source authenticity.

## Exact Inputs

The reader supports the existing `gentle.target_tss_fasta_export.v1` input
format identified in Glen's pinned input revision
`a106cbbd5223f8c4be55a7d846b8416bc4b3ae33`. That revision is not merged by this
feature. Its manifest binds the exact `SHA256SUMS` bytes; both the manifest and
checksum inventory must agree with each FASTA. The profile report/receipt bind
the manifest's digest and its declared source revision separately from GENtle's
profile-producer revision. Original FASTA files are not rewritten.

For that format, callers must supply exact `--expected-genome-id`,
`--expected-assembly` and `--expected-dataset-id` strings. For Glen's input these
are `Human GRCh38 Ensembl 116`, `GRCh38`, and
`human_grch38_ensembl116_promoterome_2000_200_v1`. Substrings such as `116` do not
match. The manifest has no separately declared annotation-release field; that
field remains null, rather than being guessed from the genome/dataset names.
The geometry uses the manifest's extents, not gene-specific runtime defaults.

With `BUNDLE`, `PANEL` and `SELECTION` set to absolute paths to that pinned input
bundle directory, panel JSON and selection JSON respectively:

```sh
gentle_cli features tss-tfbs-profiles \
  --manifest "$BUNDLE/manifest.json" --panel "$PANEL" \
  --selection "$SELECTION" \
  --expected-genome-id 'Human GRCh38 Ensembl 116' \
  --expected-assembly GRCh38 \
  --expected-dataset-id human_grch38_ensembl116_promoterome_2000_200_v1 \
  --output-dir five-gene-tss-profiles --formats svg,png,pdf
```

Use a committed producer and a fresh output directory for an audit. This is an
explicit input example, not a built-in gene or factor preset. Inspect the new
receipt rather than reusing a successful receipt from another revision.

The small synthetic example uses a second explicit development format:

`gentle.tss_fasta_bundle.v1` contains `reference` (separate `genome_id`,
`assembly`, optional `annotation_release`), `fasta_files` (relative name to SHA-256), and
`records`. Each record has `promoter_id`, `gene_id`, `gene_symbol`, `geometry`,
`transcripts`, and a normalized `sequence_sha256`. Geometry specifies chromosome,
strand, inclusive genomic bounds/TSS, and upstream/downstream extents. The
sequence length is upstream + 1 + downstream. No clipping at contig ends is
silently performed. The fixture manifest is a complete concrete example.

For this synthetic format, the manifest directory is the bundle root and
`SHA256SUMS` binds exact manifest, FASTA and requested selection bytes.
Repeated `--fasta` arguments, if supplied,
must name exactly the manifest's FASTA set; selection and FASTA paths remain
inside the bundle. Absolute contained paths work, but traversal and escaping
symlinks fail. Headers must agree with the manifest. Only ASCII-whitespace
removal and uppercasing A/C/G/T/N are allowed before sequence hashing. Missing
fields, unknown settings, duplicate identifiers and contradictory checksums fail.

The existing `gentle.regulatory_region_comparison_sequences.v1` selection
format joins `regions[].promoterome_id` to the TSS manifest's `promoter_id`.
The gene, chromosome, TSS and strand must match; selection transcript IDs must
be a subset of the manifest's memberships. Evidence-window lengths, genomic
extents and sequence hashes are deliberately not compared: a -2000/+200
selection can refer to the same TSS as a -500/+200 display window. The explicitly
supplied selection file is read-only and its exact bytes enter the receipt.
Selected panels carry the recorded factor/criterion and a legend explaining
that selection does not establish TSS usage, direct binding or promoter activity.

For simple synthetic selection, `gentle.tss_profile_selection.v1` has the same `reference` and a `selected`
array of exact `{promoter_id, gene_id}` records. Selection changes ordering,
not which windows are scored. Without it, every record is unselected. Windows
are grouped by gene ID, selected first, then chromosome, ascending interval
start and promoter ID. Equal DNA at different positions remains distinct.
Repeated physical TSSs currently fail rather than silently losing memberships;
multi-gene ownership of one physical TSS needs an additive contract extension.

The existing `gentle.jaspar_target_panel.v1` accepts `tracks[]`, each carrying
its exact `source_ids`, factor identity and display/scoring policies. Its generic
typed parser can represent mixed score kinds, but the TSS-profile gate requires
one score kind and one clipping policy. Unsupported settings are rejected, not
silently discarded. No factor or accession list is a production default.

The synthetic shorthand of the same schema has `panel_id`, `label`, `score_kind`,
`clip_negative`, `scale_mode`, `strand_policy`, `calibration_state`,
`calibration_statement`, `top_hit_count` and `factors`. Each `factors` entry is
one matrix, with exact versioned `source_id`, case-sensitive expected registry
`factor_id`, presentation `label`, ascending unique `display_order`, and optional
`color_hint` (`#RRGGBB`). Multiple accessions for one factor remain separate.
Unknown accessions, aliases, duplicate accessions, incomplete PFMs, name/order
conflicts and mixed score kinds fail. Names alone do not verify species.

## Scores And Geometry

Both motif orientations use the same transcript-oriented **window-start** axis.
For index `i`, relative coordinate is `i - upstream`; genomic coordinate is
`TSS + strand_sign * (i - upstream)`. A minus-strand input is already oriented:
it is not reverse-complemented again. Its genomic labels decrease left to right.
A reverse-motif hit's 5-prime end differs from that shared plotting coordinate;
TSV includes ascending motif interval bounds, local/genomic strand and endpoints.

Choose one existing score kind for the whole panel. `llr_bits` is the sum of
log2(model base probability / background base probability). GENtle's separately
named `true_log_odds_bits` uses per-base odds ratios. Quantile variants and
`*_background_tail_log10` variants are not raw log-odds scores: tail scores use
a negative log10 background-tail probability and the existing 0.95 quantile
display threshold. The report retains background, pseudocount, quantization,
random-seed and threshold metadata. Raw PWM scores and the pseudocount policy
are unchanged by the September 11 tail-calculation correction below. Neither a
tail probability nor a quantile is a biochemical binding probability.
The existing quantile display gate still suppresses lower-ranked tail-track
values to zero: such zeros do not mean `Ptail=1` or absence of binding. Exported
`raw_score` means the selected `score_kind` before negative clipping, **not**
necessarily LLR bits. Tail-only arrays are insufficient to recover original
bit scores; recompute those from the bound sequence and PFM when needed.

### Inclusive Background Tails

`uniform_iid_quantized_conservative_survival_v2` models independent uniform
A/C/G/T windows. Each column score is rounded to 0.001 bits for the dynamic
program. Let `E` be the sum of the largest absolute rounding error in each
column, plus a floating-point summation guard. For an observed **raw** sum `s`,
the reported inclusive tail sums bins `q >= s - E`. Every sequence with raw
score at least `s` is included. Near the threshold, additional sequences can be
included: this is a conservative **upper bound on Ptail**, hence a lower bound
on `-log10(Ptail)`, not an exact continuous-score distribution. The discretization
part of `E` is at most `motif_length * 0.0005` bits. Midrank background quantiles
use the same uncertainty interval; they are approximate ranks, not tail scores.

Survival masses are accumulated from the high end in log space, never as
`1 - CDF`. Motifs longer than 500 bases also use log-domain dynamic-program
masses. The significance score has no artificial 300 cap. A separately exported
linear probability can still underflow below the representable f64 range;
its log score remains finite. Ordinary floating-point rounding remains, and
the model is not a GC-matched genome, multiple-testing correction or affinity
model. For any achievable length-L word, significance cannot exceed
`L * log10(4)` (up to floating-point tolerance).

Regression example, exact bundled TP73 **MA0861.2**:
`ACATGTCTGGACATGT` retains raw LLR **19.543680326691 bits**. Its unique maximum
has inclusive `Ptail = 4^-16 = 2.328306436539e-10`, scoring **9.632959861247**,
instead of the old false zero tail and display score 300.

Raw arrays are stored before optional negative clipping. A positive-only figure
therefore does not redefine the inputs to correlations. Ambiguous-base windows
are `null`, not zero. The last motif-length-minus-one bases have no complete
window and are visibly unscored, never stretched across the axis.

Rows default to independent numeric scales with labeled units and an upfront
`SCALE` notice on every TSS/continuation. Equal heights across matrices or TSS
pages do not establish equal scores. `--scale-mode shared_across_tss` gives each
exact accession one numeric range across **all windows in the supplied report**,
all genes, selected and unselected, and both strands. A peak of 2 then occupies
1% of a 0..200 axis when that same accession reaches 200 elsewhere. Other
matrices retain separate ranges; no cross-matrix calibration is implied.
The effective choice is bound by `export-request.json`, `receipt.json` rendering
options and SVG `data-scale-mode`; the panel default remains in `report.json`.
For the same comparison on re-export, use the complete report, not a per-gene
subset: changing the supplied window set can change the common range.

The older `--scale-mode shared` shares a range across matrices **within each
TSS only**, not across TSSs. Cross-matrix shared scales require
`calibration_state: cross_source_calibrated`, a `calibration_id` and a lowercase
`calibration_sha256`; prose alone is insufficient. This records the caller's
explicit calibration declaration, not independently proven comparability.
PFM-derived logos show information content, not a fabricated consensus.

Compressed gene-locus overviews use one unconnected min/max whisker and mean
mark per pixel bucket, separately for forward/reverse scores. For repeating
`[0,0,0,5]` buckets the range is 0..5 and mean 1.25, not a connected plateau at
5. Native-resolution samples remain connected, with gaps at missing values.
Labels identify this aggregation. Means summarize displayed scores, not binding
occupancy. Whiskers describe finite samples in a bucket,
not confidence intervals or exact subpixel peak locations; missing samples are
not measured zero. Original arrays remain available in JSON for zoom/reanalysis.

Within-factor comparisons use common valid start coordinates, independently for
forward/forward and reverse/reverse. Pearson and Spearman consume unclipped,
unsmoothed scores; Spearman uses average ranks for ties. Reports include paired
and excluded counts. Insufficient or constant signals have undefined values
and reasons, not a misleading zero correlation. Maximum scores and separated
ranked peaks are distinct; `top_hit_count` controls the latter.

## Rescoring Existing Reports

Re-export alone cannot correct stored background-derived values. Rescore original
sequence/PFM inputs for `gentle.tss_tfbs_profiles.v1` with background-tail or
background-quantile score kinds; regenerate maxima, ranked peaks and comparisons
as well. New reports record `score_policy.modeled_tail_method`. Missing metadata
is visibly warned about on re-export; existing scores are not repaired in place.

The same shared calculation supplies `TfbsScoreTrackReport`
(`gentle.tfbs_score_tracks.v1`), its normalization references,
`gentle.window_cohort_tfbs.v1`, `gentle.promoter_cohort_comparison.v1`, and
JASPAR-derived tracks embedded in `gentle.gene_locus_evidence_display.v1` and
reporter-comparison reports. Regenerate
these background-derived tracks and dependent decisions from their original
requests and exact sequence/matrix bindings. Even raw-bit track reports need
refreshed **normalization-reference statistics** if those statistics are used;
their raw LLR/log-odds arrays do not change. Externally supplied scores are not
rescored by GENtle and require their producer's policy audit. Reports with only
raw-bit/empirical-quantile arrays and no affected background statistics need only
a presentation re-export for the scale/bucket changes. CUT&RUN signal, geometry,
FASTA, and PFM logos are not numerically changed by this fix.

Committed tutorial outputs also retain older normalization references in
`promoter_design_artifact_slice_offline`, `promoter_gene_set_ortholog_cohort_offline`
and `gene_set_ortholog_promoter_cohorts_offline`. They were deliberately not
overwritten in this correction; an explicit reviewed tutorial regeneration is
needed before treating generated-artifact freshness as a release pass.

The private original five-gene scored inputs are not available in this checkout.
The published PDF/selected-FASTA bundle is not enough to reconstruct all TSSs,
matrices and source context. Glen should use the original hash-bound manifest,
full FASTA inputs, exact panel and selection, with the corrected binary:

```sh
gentle_cli features tss-tfbs-profiles \
  --manifest /absolute/path/to/original/manifest.json \
  --panel /absolute/path/to/original/panel.json \
  --selection /absolute/path/to/original/selection.json \
  --expected-genome-id 'EXACT_ORIGINAL_GENOME_ID' \
  --expected-assembly 'EXACT_ORIGINAL_ASSEMBLY' \
  --scale-mode shared_across_tss --formats svg,png,pdf \
  --output-dir /absolute/path/to/new-rescored-bundle
```

Add the original annotation-release/dataset checks and FASTA overrides if the
original request supplied them. Keep accession versions and PFM hashes fixed;
do not substitute newly downloaded matrices. Regenerate affected locus TFBS
tracks first, keeping their original source evidence; then supply the newly
hash-bound `--context-manifest` and follow the
[integrated report regeneration checklist](integrated_locus_tss_profiles.md#refreshing-the-september-10-bundle).
Use a fresh destination; retain the old report/receipts for comparison. Verify
the scoring-method marker, matrix/sequence hashes, TP73 maximum regression,
both strands and common numeric ranges before publishing replacement figures.

## Outputs And Limits

`gentle.tss_tfbs_profiles.v1` retains matrix counts and digests, actual registry
source bindings, input digests, reference geometry, raw arrays, maxima, peaks,
comparisons, scoring policy, input source revision, profile-producer executable
digest and lockfile digest. JSON preserves floating-point score bits on replay.
Bundle checks
establish internal consistency only: independent prepared-reference sequence
verification is explicitly not assessed in this version.

Exports include a complete `report.json`, per-gene data and pages, long-form
TSV, comparisons, index, methods README and
`gentle.tss_tfbs_profile_receipt.v1`. SVG is rendered from the report; optional
PNG and single-page raster-backed PDF use the existing in-process renderer.
PDF RGB streams use lossless FlateDecode/zlib compression: every pixel, page
dimension, link annotation and raster resolution is retained, not JPEG-compressed
or downsampled. The receipt records `pdf_image_encoding`. Older uncompressed PDFs
remain valid historical artifacts; re-export changes PDF and receipt hashes,
not scientific scores. Long-form TSV serialization is unchanged. For transfer,
archive the **whole** directory (including the receipt) rather than replacing
individual files inside a hash-bound export with compressed variants.
Pages default to one TSS with full-size rows; `--panels-per-page` is bounded to
1..32. Very tall exports may exceed raster limits rather than silently shrink.

Publication uses a fresh staging directory and verifies output hashes before
publishing. Failure or cancellation does not leave a success receipt. Receipts
bind inputs, outputs, producer/exporter revision, executable, lockfile and
rendering metadata, including fonts actually used for PNG/PDF glyphs. SVG text
also depends on the viewer's fonts; cross-host byte identity is not promised.
Receipts are integrity records, not digital signatures or
proof of biological correctness. The receipt does not hash itself.

Input and output sizes are bounded; computation allows at most 4096 TSS windows,
256 matrices of up to 64 columns, and ten million strand/position scores per
request. Larger jobs must be split explicitly; oversized matrices are rejected.
Predictions are not measured binding, occupancy, affinity or promoter activity.
Matrix correlation is not evidence of co-regulation or a preferred model.
Glen reports a successful five-gene replay at
`44b73e4ac3296a9ba10caa4136a65ae49d89f289`: 58 TSSs, 13 selected windows,
30 exact matrices, 194/194 output hashes and 58 readable, visually checked PDFs.
That audit exposed uncompressed PDF size and excessive CLI stdout, addressed
by the presentation changes above. It does not certify this newer revision:
Glen must repeat exact-candidate size, hash, pixel and stdout checks. Independent
reference extraction, optional smoothing and a native wizard remain follow-ups.

## Developer Verification

Focused tests keep scientific values, malformed inputs, rendering and adapter
consent separate:

```sh
cargo test --locked -p gentle-protocol --lib tss_profiles
cargo test --locked -p gentle-engine --lib tss_profiles
cargo test --locked -p gentle-render --lib tss_profiles
cargo test --locked --no-default-features --lib tss_ -- --test-threads=1
cargo test --locked --no-default-features --lib tfbs_track_panel::tests
cargo test --locked --no-default-features --lib svg_png::tests
cargo test --locked --no-default-features --lib svg_pdf::tests
cargo test --locked --no-default-features --bin gentle_cli tss_
cargo test --locked --no-default-features --bin gentle_cli test_parse_global_args
```

The real-bundle reader check is deliberately ignored by ordinary tests: it needs
the original files rather than a download hidden inside a test. Materialize the
three documented input paths from `a106cbbd5223f8c4be55a7d846b8416bc4b3ae33`
under a separate directory, preserving `docs/examples/regulatory_region_comparison/`,
then run:

```sh
GENTLE_TSS_REAL_INPUTS=/absolute/path/to/materialized-inputs \
  cargo test --locked --no-default-features --lib \
  target_real_58_tsss_and_13_selected_read_only_acceptance -- --ignored
```

This checks all 58 records, transcript memberships and 13 selection joins. It is
not a substitute for generating and visually reviewing the five-gene outputs on
the final producer revision. The full release and native GUI gates remain with
the release auditor.

For Glen's next acceptance, use a fresh directory and retain the original
`44b73e4a` bundle. First use the report-only export command above with
`--formats svg,png,pdf`: `report.json` and all score/comparison TSV bytes must be
unchanged. Decode old/new PDF image streams and compare RGB pixels/dimensions
under the same recorded fonts; review representative plus/minus pages and the
tallest page. Verify every new receipt hash and record PDF/TSV/bundle byte totals.
Then rerun the original compute command at the exact new candidate, measure
stdout bytes, confirm the summary's verification limits and artifact hashes,
and compare saved score bits/comparisons against the retained report. Producer
revision/executable changes are expected; changed scientific values are not.
