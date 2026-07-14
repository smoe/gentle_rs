# Gene Locus Evidence Figure Runbook

GENtle can compose transcript architecture, annotation-derived transcript/CDS
metrics, selected projected BED or BigWig occupancy tracks, continuous motif
score tracks, and existing junction-qPCR candidates into one
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

## Inspect and render

Use the exact registered track names in an explicit occupancy layout. The
checked-in example keeps Saos-2 and SK-MEL-29 in separate groups and uses one
shared scale within each group:

```text
inspect-feature-expert GENE_SEQ_ID gene-locus-evidence GENE_PANEL --annotation-release Ensembl116 --upstream-bp 5000 --downstream-bp 1000 --occupancy-layout @docs/examples/gene_locus_evidence/patz1_cutrun_layout.json --motif TP73 --score-kind llr_background_tail_log10 --motif-threshold 0 --motif-top-hits 5 --qpcr-report-id QPCR_REPORT_ID

render-feature-expert-svg GENE_SEQ_ID gene-locus-evidence GENE_PANEL --annotation-release Ensembl116 --upstream-bp 5000 --downstream-bp 1000 --occupancy-layout @docs/examples/gene_locus_evidence/patz1_cutrun_layout.json --motif TP73 --score-kind llr_background_tail_log10 --motif-threshold 0 --motif-top-hits 5 --qpcr-report-id QPCR_REPORT_ID patz1_locus_evidence.svg

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

Transcript metrics report spliced exon length, CDS length, expected peptide
length, coding status, and retained-intron/incomplete flags. Green translation
start and red stop glyphs are drawn only when supported by an annotated CDS and
its derived translation. Junction-qPCR candidates are deduplicated into one
junction-spanning marker carrying all supported transcript families.

The shared Shell panel can run the same commands in the GUI. The Splicing
Expert Evidence tab remains the detailed evidence-ledger view; the composed
publication SVG is the compact locus-level view.

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
