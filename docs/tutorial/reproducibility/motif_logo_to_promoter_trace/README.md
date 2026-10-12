# Motif-score tutorial reproducibility inputs

All DNA and the four-column count matrix used in the teaching calculation are
synthetic. The counts, in A/C/G/T order by column, are
`[8,1,1,0] [1,7,1,1] [0,1,8,1] [1,1,1,7]`. They were written for this
tutorial on 2026-09-11 and carry no biological claim or third-party sequence
license. `worked_scores.tsv` is recreated by the formulae printed in the
tutorial and checked by `scripts/test_motif_score_tutorial.py`.
The table separates `score` from `units`: log2 weights use bits, while the
background-tail transformation uses `-log10(probability)`. The diagram's ten
sites are seven copies of `ACGT`, then `AAGA`, `CGCC` and `GTTG`; recounting them
recreates the PFM. Its information stacks use unsmoothed `p=count/10`,
`IC=2+sum(p*log2(p))` and base heights `50*p*IC` pixels, with zero bases omitted.
These are teaching illustrations, not a reconstruction of biological training
sequences.

`synthetic_4bp.pfm`, `tiny.fa` and `tiny.fa.fai` are the minimal independent
scanner inputs. The FASTA contains `ACGTACGCTCGAACGT` on a single line and the
FAI byte offset/line geometry is therefore deterministic. Run the scanner from
its own checkout as shown in the tutorial, substituting these absolute paths.

`shared_across_tss_panel.json` contains no copied matrix counts. It pins the
real JASPAR accession `MA0004.1` and resolves its versioned PFM from GENtle's
bundled `assets/jaspar.motifs.json` at runtime. The two plus/minus DNA windows
come from `test_files/fixtures/tss_profiles/`; their origin, hashes and
deterministic recreation are documented in that fixture's README.

The exact TP73 row uses the bundled JASPAR 2026 `MA0861.2` PFM. The matrix
identity and counts are bound by GENtle's runtime matrix digest; the tutorial
does not duplicate them. The independent scanner comparison was audited at
`IEGT/jaspar-mapping` revision
`0883ec719abee70bafbb2e8abfcae45e5bc9bcd9`.

`provenance.json` records the actual source revisions, executable SHA-256
digests, matrix/panel hashes and effective parameters used for the two verified
replays. `generated/receipt.json` independently binds the shared-engine TSS
output inventory.

Generated SVG teaching diagrams are explanatory vector illustrations, not GUI
screenshots. No private promoter package, genome download or genome-wide hit
table is part of this directory.
