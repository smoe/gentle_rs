# TP73-supported TSS-local regulatory comparison

## Scope

This report integrates the prior `CD44`, `TGFB1`, and `SERPINE1` evidence on
the same genomic axis. It uses every Ensembl Regulation feature overlapping a
selected transcript-oriented **−500/+200 bp** TSS window. Adjacent selected
windows are combined only for display; each feature is intersected with that
displayed stretch before sequence comparison.

For each connected stretch, the report shows:

1. selected TSS windows and all overlapping Ensembl regulatory features;
2. the 15 candidate/cofactor JASPAR TFBS tracks as one shared layer;
3. the 12 TP73 CUT&RUN and H3K4me3 BigWig lanes from the prior report; and
4. one promoterome-recurrence strip per regulatory-feature intersection.

The similarity search compares those exact feature intersections with 389,722
distinct GRCh38/Ensembl-116 transcript-linked promoter windows. It requires at
least 40 aligned bp, at least 80% identity, and E ≤ 1e−5. Self-locus and
same-gene targets are excluded. Counts are complete: the target cap was not
reached.

## Main observations

### CD44

Four Ensembl feature intersections span the two selected TSS windows. The
eMAR and two promoter intersections have no qualifying other-gene recurrence.
Promoter `ENSR11_9HXKD` has two short matches—to `MIR3201` and `RSL24D1P9`—but
neither covers 25% of its 501-bp sequence.

### TGFB1

The 740-bp clipped eMAR and 271-bp promoter have no qualifying other-gene
recurrence. The 22-bp CTCF feature remains visible as biological context but is
shorter than the declared fragment-similarity gate.

### SERPINE1

The distal TSS stretch contains the only broad recurrence signal. A
repeat-rich approximately 100-bp tract at **−334..−235 bp relative to TSS
101122158** lies within promoter `ENSR7_93H5NS` and its enclosing eMAR. The top
`ZNF586` promoter matches cover 28.9% of the 346-bp promoter intersection.
Twelve other genes reach at least 25% coverage; none reaches 50%.

The adjacent enhancer `ENSR7_BQ3T2` has only two short other-gene matches and
no match covering 25%. In the proximal SERPINE1 stretch, the eMAR, two
promoters, and enhancer have no qualifying recurrence. Its 22-bp CTCF feature
is again motif-scale rather than a ≥40-bp fragment query.

## Reporter-design consequence

The distal `SERPINE1` recurrent tract is the clearest boundary for an explicit
deletion/split contrast. A reporter panel could preserve the local CUT&RUN and
cofactor-TFBS context while comparing constructs that retain or omit this
tract. The relatively sequence-restricted proximal features are useful
contrasts, not automatically better promoters.

These results do **not** establish promoter activity, direct TP73 binding, or
fragment sufficiency. Ensembl feature types are annotations, JASPAR sites are
predictions, CUT&RUN is occupancy/enrichment evidence, and sequence recurrence
is structural evidence. Reporter activity still requires experimental testing.
