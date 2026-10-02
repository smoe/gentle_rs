# TP73 UniProt projection audit review evidence

This bounded bundle records the 2026-10-01 review of tutorial 06.05. It does
not contain the 9.15 GB prepared reference cache or the 31 MB transient project
state.

The reviewed run used:

- `Human GRCh38 Ensembl 116`
- sequence SHA-1 `294c432b69a88521437ffe2a9080e08d3f84902c`
- annotation SHA-1 `1cc995089649e7d6bc291b12947f2fe98989975f`
- reviewed UniProt `O15350` (`P73_HUMAN`, 636 aa)
- full annotation: 20 transcripts, 222 exons and 187 CDS features
- ten projected UniProt-linked transcripts

`review-summary.json` retains the exact compact acceptance counts and
per-transcript accounting. `tp73_uniprot_projection.svg` is the deterministic
shared expert export (SHA-256
`7cb5c32e7e0c9f605e55fca781cc3cce14028bbc0fd75ba19240ef53a95d8267`).

The SVG is included in the LF/CRLF checkout regression
`scripts.test_tutorial_checkouts.TutorialCheckoutTests.test_uniprot_review_evidence_hashes_survive_both_checkout_modes`.
The test copies the retained bytes into a disposable repository and checks this
hash, including a negative control with the scoped LF rule removed.

Ensembl REST lookup and sequence endpoints returned HTTP 500 during review.
Consequently, external Ensembl protein evidence was not fabricated: six rows
are genuine length mismatches and four are `missing_evidence`. The generated
maintainer text remained local and unsent.
