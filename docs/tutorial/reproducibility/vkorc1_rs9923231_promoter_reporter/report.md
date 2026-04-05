# VKORC1 / rs9923231 promoter-reporter handoff

## Starting point from ClawBio

ClawBio interprets a pharmacogenomic alert around warfarin sensitivity and
`VKORC1 rs9923231`. The follow-up question taken into GENtle is regulatory:
does a `VKORC1` promoter fragment containing `rs9923231` justify allele-matched
reporter constructs for a human-cell assay?

## What GENtle is responsible for

GENtle is used here for deterministic build/render work:

- retrieve the `VKORC1` locus from a prepared GRCh38 reference,
- make the reverse-strand promoter geometry explicit,
- extract one reproducible promoter fragment with the SNP inside the insert,
- import one mammalian luciferase reporter backbone,
- preview one construct,
- export one readable construct map.

## Chosen design

- assembly: `GRCh38`
- prepared genome: `Human GRCh38 Ensembl 116`
- variant: `rs9923231`
- fetched context: `vkorc1_rs9923231_context` (`+/- 3000 bp`)
- default promoter fragment: `chr16:31095780..31096868`
- extracted insert id: `vkorc1_rs9923231_promoter_ref`
- reporter backbone id: `promega_luciferase_ay738222`
- construct preview id: `vkorc1_rs9923231_luc_construct_preview`

## Biological assumptions

- this is a human-cell regulatory assay story,
- the baseline reporter is a promoterless mammalian luciferase backbone,
- the promoter fragment keeps `rs9923231` inside the insert rather than at the
  edge,
- retaining about `200 bp` past the TSS is intentional and is meant to avoid a
  cut exactly at the annotated TSS,
- adenoviral delivery is not the baseline and would only be considered later if
  transfection efficiency becomes the bottleneck.

## Non-goals

- no claim of wet-lab validation,
- no direct claim that the construct alone proves warfarin response,
- no bacterial-expression framing.

## Deterministic replay material

- tutorial page:
  `docs/tutorial/vkorc1_warfarin_promoter_luciferase_gui.md`
- workflow skeleton:
  `docs/examples/workflows/vkorc1_rs9923231_promoter_luciferase_assay_planning.json`
- commands:
  `docs/tutorial/reproducibility/vkorc1_rs9923231_promoter_reporter/commands.sh`
- structured summary:
  `docs/tutorial/reproducibility/vkorc1_rs9923231_promoter_reporter/result.json`

## Bench-facing handoff

The immediate next practical step is to build the reference construct exactly as
documented here, then branch an alternate-allele construct with the same
fragment geometry so promoter output can be compared without changing the
insert boundaries.

