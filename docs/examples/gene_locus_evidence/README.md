# Gene-locus occupancy layout example

`patz1_cutrun_layout.json` is a deterministic display-layout example for the
local PATZ1 CUT&RUN workflow. It contains no experimental measurements. Track
names must match projected GENtle track names exactly; the BigWig/BED files and
their biological provenance remain external resources.

The two cell lines intentionally use separate `shared_group` scales. Do not
switch to `shared_all` unless preprocessing and normalization make the groups
quantitatively comparable and the layout records that rationale in
`cross_group_scale_justification`.
