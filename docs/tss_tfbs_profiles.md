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
random-seed and threshold metadata. No existing scoring formula is changed.

Raw arrays are stored before optional negative clipping. A positive-only figure
therefore does not redefine the inputs to correlations. Ambiguous-base windows
are `null`, not zero. The last motif-length-minus-one bases have no complete
window and are visibly unscored, never stretched across the axis.

Rows default to independent numeric scales with labeled units. Equal heights
across matrices or TSS pages do not establish equal scores. Shared scales require
`calibration_state: cross_source_calibrated`, a `calibration_id` and a lowercase
`calibration_sha256`; prose alone is insufficient. This records the caller's
explicit calibration declaration, not independently proven comparability.
PFM-derived logos show information content, not a fabricated consensus.

Within-factor comparisons use common valid start coordinates, independently for
forward/forward and reverse/reverse. Pearson and Spearman consume unclipped,
unsmoothed scores; Spearman uses average ranks for ties. Reports include paired
and excluded counts. Insufficient or constant signals have undefined values
and reasons, not a misleading zero correlation. Maximum scores and separated
ranked peaks are distinct; `top_hit_count` controls the latter.

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
