# Synthetic Allele-Pair Guard

`multiallelic.gb` is hand-crafted test data, not a human genomic sequence,
clinical allele report or promoter-activity evidence. Its labels deliberately
exercise the same default GUI IDs as Tutorial 08.04; they do not confer
authentic VKORC1 provenance.

Recreate the payload as six `A` bases, one `C`, then thirteen `A` bases (20 bp).
Add one forward `variation` at GenBank base 7, with the literal qualifiers in
the file: `label=rs9923231`, `gene=VKORC1`, `vcf_ref=C`, `vcf_alt=A,G,T`.
The single feature is index 0. No imported genome, transcript, TFBS, population
or experimental annotation is implied.

Used by `vkorc1_pair_gui_starter` and `vkorc1_pair_gui_oracle` through the normal
GenBank loader and shared `MaterializeVariantAllele` engine operations. The
08.04 offline-core GUI contract starts without outputs and must create a
reviewed reference/T pair matching the independent oracle. The expected
genomic-forward difference is `C -> T` at zero-based position 6 only. The
adapter's synthetic refusal tests additionally require no project/history
mutation when the alternate is omitted or invalid. Native acceptance remains
separate from those deterministic tests and from the full online tutorial.
