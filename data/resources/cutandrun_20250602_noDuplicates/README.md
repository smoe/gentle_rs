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

Before interpreting or publishing a figure, record the alignment assembly,
normalization and clipping procedure, and the meaning of the `tp73`, `pos`,
`neg`, and `input` channels. A matching chromosome name is not by itself proof
that a track and an open sequence use the same assembly.

See [`docs/gene_isoform_occupancy_figure_runbook.md`](../../../docs/gene_isoform_occupancy_figure_runbook.md)
for projection and combined transcript/occupancy rendering commands. Public raw
data provenance remains recorded under
`data/publication_resources/rostock_p73_cutrun_e_mtab_15709/`.
