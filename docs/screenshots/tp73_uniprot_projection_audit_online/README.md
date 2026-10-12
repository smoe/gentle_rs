# TP73 UniProt projection audit screenshots

Native Linux/Xvfb captures from GENtle built from this review branch on
2026-10-01. The loaded project used Human GRCh38 Ensembl 116, full TP73
annotation, reviewed UniProt O15350 and the corrected genomic-CDS accounting
path.

| File | What it proves | SHA-256 |
| --- | --- | --- |
| `01-audit-mismatch-vs-missing-evidence.png` | The native audit distinguishes six protein-length mismatches from four rows limited by missing Ensembl evidence, exposes corrected nucleotide/AA accounting, and labels the maintainer text as an unsent draft. | `f7fb5b7528496fcd5d639d0646b8e9c5293597ad1ecd252589f2bf8531783cff` |
| `02-direct-vs-composed-parity.png` | The stored parity report contains ten rows, zero divergent rows, and a matching draft-transcript set. | `6ea28522c36d6996faae08272eb79ee15183b8121c44c0bf019e7182af9b59ff` |

Both images are 1920x1080 PNG captures. They are display acceptance evidence;
the structured report remains the semantic source of truth.
