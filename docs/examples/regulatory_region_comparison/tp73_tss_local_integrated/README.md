# TSS-local integrated regulatory comparison

This retained analysis revises the first TP73 CUT&RUN-supported promoterome
comparison around transcript-oriented **−500/+200 bp** windows. For each of the
seven selected `CD44`, `TGFB1`, and `SERPINE1` TSS windows it intersects every
overlapping Ensembl Regulation feature with the displayed genomic stretch.
Feature similarity is therefore computed only for the sequence visible in this
TSS-local report, not for a distant portion of a larger eMAR annotation.

Each connected stretch has one shared 15-factor JASPAR TFBS layer and the same
12 TP73 CUT&RUN/H3K4me3 BigWig lanes used by the prior luciferase-planning
reports. Beneath them, one recurrence panel per Ensembl feature shows matches
to the 389,722-window GRCh38/Ensembl-116 promoterome. Self-locus and same-gene
targets are excluded. Rows represent distinct genomic promoter windows; genes
and transcripts sharing one TSS stay attached to that single occurrence.

## Reproduce

Prepare the exact intersections from the selected TSS bundle, hash-verified
GENtle locus reports, and GRCh38 reference:

```bash
python3 scripts/prepare_tss_regulatory_similarity_candidates.py \
  --selected-tss /path/to/candidate_regions.json \
  --locus-report /path/to/CD44.report.json \
  --locus-report /path/to/TGFB1.report.json \
  --locus-report /path/to/SERPINE1.report.json \
  --reference-fasta /path/to/GRCh38.fa \
  --output /path/to/tss-local/candidates \
  --upstream-bp 500 --downstream-bp 200
```

Compare all 15 feature intersections with the transcript-linked promoterome:

```bash
python3 scripts/compare_candidates_to_promoterome.py \
  --candidates-json /path/to/tss-local/candidates/candidate_regions.json \
  --query-fasta /path/to/tss-local/candidates/candidate_regions.fa \
  --promoterome /path/to/human-grch38-ensembl116-promoterome \
  --output /path/to/tss-local/blastn-40bp-80pct \
  --task blastn --min-alignment-bp 40 --min-identity-pct 80 \
  --max-evalue 1e-5 --max-target-seqs 1000000 --max-hsps 10
```

Render the integrated report:

```bash
python3 scripts/render_integrated_tss_regulatory_report.py \
  --candidates-json /path/to/tss-local/candidates/candidate_regions.json \
  --comparison /path/to/tss-local/blastn-40bp-80pct/comparison.json \
  --matches /path/to/tss-local/blastn-40bp-80pct/matches.blastn.tsv \
  --hits /path/to/tss-local/blastn-40bp-80pct/hits.blastn.tsv \
  --locus-report /path/to/CD44.report.json \
  --locus-report /path/to/TGFB1.report.json \
  --locus-report /path/to/SERPINE1.report.json \
  --output /path/to/tss-local/reports
```

The 108 MiB full BLAST tables are regenerable and intentionally not versioned.
`receipt.json` binds them by SHA-256. The repository retains the exact feature
candidates, compact comparison metadata/tables, integrated PDFs/PNGs, and
interpretation.

## First result

- `CD44`: three of four feature intersections have no qualifying other-gene
  hit. `ENSR11_9HXKD` has two short matches, but neither covers 25% of the
  501-bp promoter feature.
- `TGFB1`: neither the 740-bp clipped eMAR nor the 271-bp promoter has a
  qualifying other-gene hit. The 22-bp CTCF feature is displayed but is below
  the declared 40-bp fragment-search minimum.
- Distal `SERPINE1`: a repeat-rich approximately 100-bp tract at
  **−334..−235 bp relative to TSS 101122158** recurs widely. It lies in promoter
  `ENSR7_93H5NS` and its enclosing eMAR. The leading `ZNF586` promoter matches
  cover 28.9% of the 346-bp promoter; 12 other genes reach at least 25%, and
  none reaches 50%.
- Proximal `SERPINE1`: the eMAR, two promoters, and enhancer have no qualifying
  other-gene recurrence. Its 22-bp CTCF feature is motif-scale and untested by
  this fragment gate.

This result motivates a `SERPINE1` reporter deletion/split contrast around the
recurrent distal tract while preserving the local CUT&RUN and cofactor-TFBS
context. It does not establish that the tract is functional, dispensable, or
sufficient. Ensembl feature classes are annotations, JASPAR sites are
predictions, and CUT&RUN is occupancy/enrichment evidence rather than proof of
direct binding.
