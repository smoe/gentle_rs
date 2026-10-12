# Local processed CUT&RUN tracks

This directory is the expected local staging location for the processed
BigWig, bedGraph, and BED files supplied for the Rostock p73 CUT&RUN analysis.
The multi-gigabyte track files are intentionally not committed. Only this
location contract and usage guidance live in the repository.

Expected names include condition-specific tracks such as:

- `tp73_saos2_TA_R1.bigWig`
- `tp73_saos2_DN_R1.bigWig`
- `tp73_skmel29_2_TA_R1.bigWig`
- `tp73_skmel29_2_DN_R1.bigWig`
- matching `input_*`, `pos_*`, and `neg_*` controls where available

For public study **E-MTAB-15709**, `pos_DNp73b_SKMel29`, `pos_GFP_SKMel29`
and `pos_TAp73a_SKMel29` are **H3K4me3** libraries: every corresponding row in
the [public SDRF](https://www.ebi.ac.uk/biostudies/files/E-MTAB-15709/E-MTAB-15709.sdrf.txt)
records `Factor Value[immunoprecipitate] = H3K4me3`. These are neither HDAC2
libraries nor controls with an unspecified target. This annotation concerns
those public libraries; bind local processed tracks to their actual source
library rather than inferring control identity from a filename alone.
The SDRF checked on 2026-10-13 is 11,940 bytes, SHA-256
`60a3cdd3b5ba8c3e2b7f3fa823f3911dbf5091c61b90a31dd453bc2efad1377d`.

Before interpreting or publishing a figure, record the alignment assembly,
normalization and clipping procedure, and the meaning of the `tp73`, `pos`,
`neg`, and `input` channels. A matching chromosome name is not by itself proof
that a track and an open sequence use the same assembly.

See [`docs/gene_isoform_occupancy_figure_runbook.md`](../../../docs/gene_isoform_occupancy_figure_runbook.md)
for projection and combined transcript/occupancy rendering commands. Public raw
data provenance remains recorded under
`data/publication_resources/rostock_p73_cutrun_e_mtab_15709/`.
