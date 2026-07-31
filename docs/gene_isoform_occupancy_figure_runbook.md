# Gene Locus Evidence Figure Runbook

GENtle can compose transcript architecture, annotation-derived transcript/CDS
metrics, selected projected BED or BigWig occupancy tracks, continuous motif
score tracks, coordinate-aligned PSR/JUC probe-effect rows, and existing
junction-qPCR candidates into one
`gentle.gene_locus_evidence_display.v1` SVG. This is intended for selected-gene
CUT&RUN figures such as TAp73alpha versus DNp73beta at PATZ1. It is reusable for
other genes and occupancy assays.

## Interpretation boundary

The figure places transcript geometry, locus-level occupancy, motif scores, and
assay candidates on one strand-aware 5'-to-3' coordinate axis. It does not
assign a peak to one transcript, infer isoform-specific regulation, or convert
spatial overlap into causality. Expression, RNA-read, probe, cDNA/EST,
occupancy, motif, and qPCR evidence remain separate report layers.

Before projection, verify and record the track assembly independently. BigWig
chromosome names can detect obvious chromosome mismatches but cannot establish
GRCh37 versus GRCh38 identity. The open sequence should carry a verified genome
anchor for the same assembly.

## Prepare the project

1. Load or extract one genome-anchored gene sequence containing its transcript
   annotations.
2. Import an isoform panel whose transcript ids map to those annotations:

```text
panels import-isoform GENE_SEQ_ID PANEL.json --panel-id GENE_PANEL
```

3. Register descriptive, condition-specific track names. The example paths are
   local external resources and are not committed:

```text
tracks tracked add data/resources/cutandrun_20250602_noDuplicates/tp73_saos2_TA_R1.bigWig --source bigwig --name "SAOS-2 TP73 TAp73alpha"
tracks tracked add data/resources/cutandrun_20250602_noDuplicates/tp73_saos2_DN_R1.bigWig --source bigwig --name "SAOS-2 TP73 DNp73beta"
tracks tracked add data/resources/cutandrun_20250602_noDuplicates/tp73_skmel29_2_TA_R1.bigWig --source bigwig --name "SK-MEL-29 TP73 TAp73alpha"
tracks tracked add data/resources/cutandrun_20250602_noDuplicates/tp73_skmel29_2_DN_R1.bigWig --source bigwig --name "SK-MEL-29 TP73 DNp73beta"
tracks tracked apply
```

BigWig projection uses `bigWigToBedGraph`; set
`GENTLE_BIGWIG_TO_BEDGRAPH_BIN` when the executable is not on `PATH`. Peak-only
figures may register the matching `.clipped.clean.bed` files with
`--source bed` instead.

### Prepare probe-effect evidence

The committed PATZ1 development table substitutes for the otherwise local
`analysis/e_mtab_14704_tp73_microarray/` tree:

```text
test_files/fixtures/gene_locus_evidence/patz1_probe_effects/patz1_clariom_probe_effects.tsv
```

It contains 24 Clariom D PSR/JUC rows, 130 underlying PM probes, explicit hg38
coordinates, three condition-wise log2 abundance means, and three pairwise raw
activity differences. The fixture README records its public dataset and
derived-table provenance. GENtle uses this compact table directly; the
uncommitted analysis directory is only an optional source-level cross-check and
is not required by the command or tests.

## Offline proof before real data

Run the committed synthetic negative-strand acceptance case first:

```text
cargo run --bin gentle_cli -- workflow @docs/examples/workflows/patz1_gene_locus_evidence_offline.json
```

It writes `patz1_gene_locus_evidence.svg` without network access. The workflow
uses one synthetic PSR row, one synthetic JUC row, two raw-effect contrasts,
four small BED tracks in two declared groups, and the local TP73 JASPAR matrix.
Fixture origin and deterministic recreation are documented in
`test_files/fixtures/gene_locus_evidence/patz1_offline_composer/README.md`.
The generated tutorial chapter embeds the same SVG. This proof checks
composition and provenance only; it is not experimental PATZ1 evidence.

## Inspect and render

Use the exact registered track names in an explicit occupancy layout. The
checked-in example keeps Saos-2 and SK-MEL-29 in separate groups and uses one
shared scale within each group:

```text
inspect-feature-expert GENE_SEQ_ID gene-locus-evidence GENE_PANEL --annotation-release Ensembl116 --probe-effect-table test_files/fixtures/gene_locus_evidence/patz1_probe_effects/patz1_clariom_probe_effects.tsv --probe-effect-contrast TAp73alpha-GFP --probe-effect-contrast DNp73beta-GFP --probe-effect-coordinate-system GRCh38.p14 --upstream-bp 5000 --downstream-bp 1000 --occupancy-layout @docs/examples/gene_locus_evidence/patz1_cutrun_layout.json --motif TP73 --score-kind llr_background_tail_log10 --motif-threshold 0 --motif-top-hits 5 --qpcr-report-id QPCR_REPORT_ID

render-feature-expert-svg GENE_SEQ_ID gene-locus-evidence GENE_PANEL --annotation-release Ensembl116 --probe-effect-table test_files/fixtures/gene_locus_evidence/patz1_probe_effects/patz1_clariom_probe_effects.tsv --probe-effect-contrast TAp73alpha-GFP --probe-effect-contrast DNp73beta-GFP --probe-effect-coordinate-system GRCh38.p14 --upstream-bp 5000 --downstream-bp 1000 --occupancy-layout @docs/examples/gene_locus_evidence/patz1_cutrun_layout.json --motif TP73 --score-kind llr_background_tail_log10 --motif-threshold 0 --motif-top-hits 5 --qpcr-report-id QPCR_REPORT_ID patz1_locus_evidence.svg

svg-png patz1_locus_evidence.svg patz1_locus_evidence.png --scale 2
svg-pdf patz1_locus_evidence.svg patz1_locus_evidence.pdf
```

Omit `--qpcr-report-id` when no persisted assay report should be shown. Repeat
`--motif` for additional matrices. `--motif-threshold 0` labels only scores
strictly above zero; the continuous score trace is retained regardless. Motif
provenance records the local JASPAR matrix id, score kind, and clipping policy.

For a quick ungrouped inspection, repeat `--occupancy-track NAME` instead of
using `--occupancy-layout`; `--occupancy-track '*'` includes every projected
BED/BigWig track overlapping the locus. Explicit layout files are preferable
for publication figures because they prevent unrelated tracks from entering,
declare lane roles, and avoid mixing distinct biological contexts.

Each occupancy group can use `shared_group`, `independent`, `fixed`, or
`shared_all` scaling. `shared_all` should only be used when preprocessing and
normalization make groups quantitatively comparable; include
`cross_group_scale_justification` in that case. Unscored BED peaks are shown as
interval marks. The JSON report retains local and genomic interval coordinates,
source paths, score ranges, and explicit non-causal evidence notes.

Probe-effect tables require an explicit coordinate system matching the open
sequence anchor. Each PSR is drawn once as a locus interval and each JUC once as
a junction arc. `log2_mean_*` columns form an abundance lane with sequential
scaling; `log2_*_minus_*` columns form a separately scaled, zero-centered
differential-activity lane. This separation prevents the larger abundance
values from flattening the contrast colors. Missing values stay `NA`; genomic
coordinates and table-row provenance remain in JSON and SVG attributes. Array
signal and PSR/JUC geometry can prioritize assay regions and flag alternative
structure, but do not prove PCR-primer binding, formal differential-expression
significance, or isoform identity.

Transcript metrics report spliced exon length, CDS length, expected peptide
length, coding status, and retained-intron/incomplete flags. Green translation
start and red stop glyphs are drawn only when supported by an annotated CDS and
its derived translation. Junction-qPCR candidates are deduplicated into one
junction-spanning marker carrying all supported transcript families.

The shared Shell panel can run the same commands in the GUI. The Splicing
Expert now also provides a graphical `Locus figure` tab:

1. Enter the imported panel and optional RNA/cDNA, probe, expression, and qPCR
   sources under `Evidence`.
2. In `Locus figure`, add probe-effect table paths and contrasts, the explicit
   coordinate system, occupancy-layout path, motifs, flanks, and score options.
3. Review the resource table. Missing local files can be relocated with the
   adjacent browse action; panel and projected-track identities remain
   engine-validated project objects.
4. Use `Compose / refresh` to run the shared pure-read operation and inspect
   the in-GUI SVG preview, metrics, effects, warnings, and provenance.
5. Export report JSON, SVG, or PDF. Existing qPCR rows can open their stored
   report or create a deterministic oligo-order form; junction rows can seed
   the established transcript-aware qPCR designer.

The detailed `Evidence` tab remains the audit surface and `Locus figure` the
compact composition surface. Relocation changes the requested source path; it
does not bypass schema, anchor, or coordinate checks. Assay continuations are
planning actions and do not turn a displayed association into validation.

## Reproducibility checklist

- genome assembly and sequence anchor verified;
- annotation release and isoform-panel source recorded;
- track processing, normalization, clipping, and control-channel meanings
  recorded outside the figure;
- track names and the layout identify cell line, antibody/target channel,
  condition, replicate, lane role, and scale policy;
- input/negative/GFP controls inspected alongside experimental tracks where
  scientifically appropriate;
- JSON report, occupancy-layout JSON, and SVG retained together.
