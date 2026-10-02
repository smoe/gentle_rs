# TP53 UniProt projection GUI evidence

These are direct Linux/X11 captures from GENtle opened on the disposable state
described in
`docs/tutorial/reproducibility/tp53_uniprot_projection_online/README.md`.

- `01-protein-expert.png`: 1100x800 native Protein Expert showing the TP53
  transcript/protein rows, mapped status and explicit missing-transcript
  warnings; SHA-256
  `cc2e9f2e0bd7bdaaee7f58c09b0c29411f18afe80a6b23f577f829487fd0ed3e`
- `02-feature-coding-dna.png`: 1100x800 native Protein Evidence view with the
  `DNA-binding` result expanded to show amino acids 368-387, exon 11, the exact
  60-base genomic coding DNA and its separately labelled optimized
  alternative; SHA-256
  `4e53c82c4d4c85d3aed67c56d642383dd4eb74624dcf40d545e5b1ccea6a8744`
- source base: `80260957a6912dd2aedbca724354c8c8727e6822`
- GUI binary revision shown at capture: `v0.1.1790819522996`
- environment: isolated 1920x1080 Xvfb display; child windows captured only

The native query required selecting both pieces of existing context: **Use**
on the recent projection row and **Select** on the imported `P04637` row. The
screenshots are evidence of GUI read-back, not evidence that Ensembl reference
preparation succeeded during this review.
