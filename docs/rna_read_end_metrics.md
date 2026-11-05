# Transcript-oriented RNA read-end metrics

GENtle's read-end utilities separate three questions that must not be collapsed:

1. **whole-read length** — a library and sequencing control;
2. **5′-UTR reach** — a useful truncation/control comparison; and
3. **3′-UTR reach** — the annotation-relative measurement of interest.

The utilities are deliberately downstream of the existing observed-SAM engine.
They do not align reads, choose a transcript silently, or infer a poly(A) site.

## 1. Freeze the transcript model

`scripts/project_rna_read_end_metrics.py` requires a CDS-aware model document
with schema `gentle.rna_read_end_transcript_models.v1`. Coordinates are genomic,
zero-based and half-open. Each transcript has stable transcript and gene IDs,
contig, strand, exons and CDS intervals. The document also binds the assembly,
reference FASTA SHA-256 and annotation SHA-256 already present in the completed
`gentle.observed_sam_structures_data.v1` artifact.

This extra binding is necessary because the observed-SAM report retains exon
comparisons but does not itself define CDS boundaries. A mismatched model is
refused.

Prepare the model from an explicit GTF and the exact reference FASTA. Supply
the annotation digest recorded in the completed observed-SAM analysis; it is
kept separate from the raw GTF digest recorded as model provenance.

```sh
python3 scripts/prepare_rna_read_end_transcript_models.py \
  --gtf annotation.gtf.gz \
  --reference-fasta reference.fa \
  --assembly GRCh38 \
  --source-id ensembl-116 \
  --annotation-sha256 sha256:... \
  --output frozen-transcript-models.json
```

## 2. Project one sample

![Transcript-oriented UTR reach and right-censoring](figures/rna_read_end_projection_geometry.svg)

The solid read segment is projected through annotated exons. The 5′ and 3′
measurements begin at the corresponding CDS boundary. A clipped end is a lower
bound; it is not an observed short UTR.

```sh
python3 scripts/project_rna_read_end_metrics.py \
  --observed-sam-data sample2.observed-sam.data.json \
  --transcript-models frozen-transcript-models.json \
  --sample-id sample2 \
  --output-dir sample2-read-end-metrics
```

The fresh output directory contains:

- `metrics.tsv`: the four-column input expected by
  `compare_rna_read_end_metrics.py`;
- `projections.tsv`: per-read metrics plus transcript and gene identities; and
- `report.json`: exact input hashes, policy, counts and nonclaims.

Only saved optimal hypotheses and resolved `exact_full_chain`,
`exact_partial_chain` or `within_exon` transcript relations are eligible. A UTR
value is retained across transcript ambiguity only when all eligible paths and
models give the same value and status.

An aligned end that is clipped is **right-censored**: its value is a lower
bound, not an exact length. An alignment that does not cross the relevant CDS
boundary is **missing**, not zero.

## 3. Compare all six samples globally

Concatenate the six `metrics.tsv` files (retaining one header) and run:

```sh
python3 scripts/compare_rna_read_end_metrics.py \
  --input six-sample-metrics.tsv \
  --output-dir six-sample-qq
```

This produces the 15 sample-pair comparisons for each metric. Read the three
Q–Q matrices together: a 3′-UTR displacement is more interesting when the
whole-read and 5′-UTR controls do not show the same displacement.

## 4. Rank genes for the TP73-low sample

If `sample2` is the TP73-low sample, use the detailed projection files:

```sh
python3 scripts/rank_rna_read_end_genes.py \
  --projection sample1-read-end-metrics/projections.tsv \
  --projection sample2-read-end-metrics/projections.tsv \
  --projection sample3-read-end-metrics/projections.tsv \
  --projection sample4-read-end-metrics/projections.tsv \
  --projection sample5-read-end-metrics/projections.tsv \
  --projection sample6-read-end-metrics/projections.tsv \
  --focus-sample sample2 \
  --minimum-observed 10 \
  --maximum-censored-fraction 0.25 \
  --output-dir sample2-gene-ranking
```

The score is the median across five **sample-level** differences of read-level
medians. Individual reads are not treated as biological replicates. The report
keeps censored and missing counts and excludes a comparison when its evidence
is too sparse or too censored. A gene enters the ranking only when at least
half of the available comparison samples are eligible (three of five in the
six-sample design).

The resulting candidates can motivate follow-up around the TP73–MDM2/MDM4
hypothesis. They do not demonstrate that TP73, MDM2 or MDM4 caused a change.
Library preparation, read truncation, mapping, transcript selection and the
frozen annotation remain alternative explanations.

## Interpretation picture

![Synthetic gene-wise read-end comparison focused on sample2](figures/rna_read_end_gene_ranking.synthetic.svg)

The committed figure is synthetic. Regenerate its six-sample fixture and SVG:

```sh
python3 scripts/build_rna_read_end_tutorial.py
```

The three panels deliberately use separate scales. Compare direction and
candidate ordering, not the apparent horizontal distance between panels.
