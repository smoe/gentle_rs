# Tutorial Input Provenance

These files are pinned local tutorial inputs for offline, reproducible GENtle
walkthroughs.

## `gentle_mammalian_luciferase_backbone_v1.gb`

- Status: synthetic, repository-owned tutorial asset
- Intended use:
  - offline VKORC1 / `rs9923231` promoter-reporter planning
  - deterministic circular-map rendering and construct preview without live
    GenBank fetches
- Biological role:
  - promoterless mammalian luciferase reporter backbone for transient
    transfection planning in human cells
- Important limitation:
  - this is a teaching/demo backbone architecture, not a claim that the file is
    an exact deposited plasmid sequence

When the tutorial wants a reproducible local baseline, it should import this
record through `LoadFile` instead of relying on online `AY738222` retrieval.
