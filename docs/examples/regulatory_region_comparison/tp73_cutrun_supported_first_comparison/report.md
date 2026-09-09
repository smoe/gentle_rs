# First TP73 CUT&RUN-supported promoterome comparison

This run compares every distinct selected-gene TSS window that passes the declared matched-control TP73 BigWig rule against 389,722 GRCh38/Ensembl-116 transcript-linked promoter windows.

## Parameters

- Window: transcript-oriented −2,000/+200 bp around each annotated TSS.
- Support: at least one TAp73α or DNp73β mean exceeds the matched GFP-control mean in the same cell line across that exact window.
- Similarity: BLASTN, ≥40 aligned bp, ≥80% identity, E ≤ 1e−5; dust and soft masking enabled.
- Counting: self-locus overlaps excluded; promoter windows, genes, and transcripts counted separately.

## First observations

- The SERPINE1 ENST00000950058 window contains a highly recurrent component: its strongest MAN2A2 promoter matches cover about 47% of the proposed window.
- The two promoter-proximal SERPINE1 windows have many short/repeat-rich hits, but only tens of other genes reach 25% aggregate coverage and none reaches 50%.
- CD44 and TGFB1 have only sparse short matches under this threshold; no other gene reaches 25% aggregate window coverage.
- The CD44 ENST00000428726 window has no qualifying other-gene promoter hit in this BLASTN pass.

These are sequence-recurrence observations, not reporter recommendations. The next decision layer should intersect the recurrent blocks with localized CUT&RUN peaks, Ensembl Regulation features, motif/module evidence, and planned deletion/split constructs.
