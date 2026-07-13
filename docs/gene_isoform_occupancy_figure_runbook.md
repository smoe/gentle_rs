# Gene Isoform And Occupancy Figure Runbook

GENtle can align selected projected BED or BigWig occupancy tracks beneath all
transcript models in a `gentle.gene_isoform_evidence.v1` SVG. This is intended
for selected-gene CUT&RUN figures such as TAp73alpha versus DNp73beta at PATZ1.
It is reusable for other genes and occupancy assays.

## Interpretation boundary

The figure places transcript geometry and locus-level occupancy on one
strand-aware coordinate axis. It does not assign a peak to one transcript,
infer isoform-specific regulation, or convert spatial overlap into causality.
Expression, RNA-read, probe, cDNA/EST, occupancy, and qPCR evidence remain
separate report layers.

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
tracks tracked apply
```

BigWig projection uses `bigWigToBedGraph`; set
`GENTLE_BIGWIG_TO_BEDGRAPH_BIN` when the executable is not on `PATH`. Peak-only
figures may register the matching `.clipped.clean.bed` files with
`--source bed` instead.

## Inspect and render

Use the exact registered track names to control lane inclusion and order:

```text
inspect-feature-expert GENE_SEQ_ID isoform-evidence GENE_PANEL --annotation-release Ensembl116 --occupancy-track "SAOS-2 TP73 TAp73alpha" --occupancy-track "SAOS-2 TP73 DNp73beta"

render-feature-expert-svg GENE_SEQ_ID isoform-evidence GENE_PANEL --annotation-release Ensembl116 --occupancy-track "SAOS-2 TP73 TAp73alpha" --occupancy-track "SAOS-2 TP73 DNp73beta" gene_isoforms_cutrun.svg
```

`--occupancy-track '*'` includes every projected BED/BigWig track overlapping
the inspected gene span. Explicit names are preferable for publication figures
because they prevent unrelated project tracks from entering the result.

The SVG uses one shared absolute score scale across lanes. Unscored BED peaks
are shown as interval marks. The JSON report retains local and genomic interval
coordinates, source paths, score ranges, and an explicit non-causal evidence
note. In the GUI, the same request is available in **Splicing Expert ->
Evidence** through the **Occupancy track names** field; use `*` there to include
all projected tracks.

## Reproducibility checklist

- genome assembly and sequence anchor verified;
- annotation release and isoform-panel source recorded;
- track processing, normalization, clipping, and control-channel meanings
  recorded outside the figure;
- track names encode cell line, antibody/target channel, condition, and
  replicate;
- input/negative/GFP controls inspected alongside experimental tracks where
  scientifically appropriate;
- JSON report and SVG retained together.
