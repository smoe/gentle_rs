# VKORC1 / rs9923231 promoter-reporter handoff

## Interpretation coming from ClawBio

ClawBio interprets a pharmacogenomic alert around warfarin and `VKORC1
rs9923231` and proposes a regulatory follow-up.

The narrow question handed to GENtle is:

> Is `rs9923231` promoter-proximal enough that it is worth building matched
> promoter-reporter constructs for follow-up in human cells?

ClawBio is therefore responsible for:

- the pharmacogenomic interpretation
- the motivation for a regulatory assay
- the choice to follow up in a human-cell luciferase reporter format

## Deterministic build/render work done in GENtle

GENtle is responsible for:

- retrieving the prepared GRCh38 locus
- deriving transcript-TSS-centered promoter windows
- classifying the SNP as promoter-overlapping context
- suggesting the reporter fragment
- materializing matched reference and alternate inserts
- loading a pinned local mammalian luciferase backbone
- previewing reference and alternate reporter constructs
- exporting reviewable artifacts

## Chosen baseline design

- prepared genome: `Human GRCh38 Ensembl 116`
- context sequence id: `vkorc1_rs9923231_context`
- variant: `rs9923231`
- chosen gene: `VKORC1`
- promoter window defaults:
  - upstream = `1000 bp`
  - downstream = `200 bp`
- recommended promoter fragment interval:
  - local `2412..3501`
  - extracted id = `vkorc1_rs9923231_promoter_fragment`
- matched inserts:
  - `vkorc1_rs9923231_promoter_reference`
  - `vkorc1_rs9923231_promoter_alternate`
- pinned local backbone:
  - `gentle_mammalian_luciferase_backbone_v1`
- construct previews:
  - `vkorc1_rs9923231_reporter_reference`
  - `vkorc1_rs9923231_reporter_alternate`
  - SVG export targets are defined by the workflow, but a clean live replay is
    still required before those SVGs can be committed into this bundle

## Why this backbone is appropriate for human-cell work

- it is a promoterless mammalian luciferase reporter architecture
- it keeps the readout focused on promoter output rather than bacterial
  protein-expression logic
- it matches transient transfection planning in human cells
- it gives the tutorial one pinned local backbone so the canonical handoff does
  not depend on live GenBank retrieval

## Important assumptions

- this is a human-cell regulatory assay story
- this bundle stops at a reproducible design/handoff point
- it does **not** claim wet-lab validation
- it does **not** claim that the construct alone proves warfarin response
- adenoviral delivery is deferred as a later escalation only if transfection
  efficiency becomes the bottleneck

## Artifacts

- tutorial:
  [docs/tutorial/vkorc1_warfarin_promoter_luciferase_gui.md](/Users/u005069/.codex/worktrees/47dd/gentle_rs/docs/tutorial/vkorc1_warfarin_promoter_luciferase_gui.md)
- workflow:
  [docs/examples/workflows/vkorc1_rs9923231_promoter_luciferase_assay_planning.json](/Users/u005069/.codex/worktrees/47dd/gentle_rs/docs/examples/workflows/vkorc1_rs9923231_promoter_luciferase_assay_planning.json)
- promoter-context JSON:
  [variant_promoter_context.json](/Users/u005069/.codex/worktrees/47dd/gentle_rs/docs/tutorial/reproducibility/vkorc1_rs9923231_promoter_reporter/variant_promoter_context.json)
- promoter-candidate JSON:
  [promoter_reporter_candidates.json](/Users/u005069/.codex/worktrees/47dd/gentle_rs/docs/tutorial/reproducibility/vkorc1_rs9923231_promoter_reporter/promoter_reporter_candidates.json)
- promoter-context SVG:
  [vkorc1_rs9923231_promoter_context.svg](/Users/u005069/.codex/worktrees/47dd/gentle_rs/docs/tutorial/reproducibility/vkorc1_rs9923231_promoter_reporter/vkorc1_rs9923231_promoter_context.svg)
- expected reference construct SVG path after full replay:
  `docs/tutorial/reproducibility/vkorc1_rs9923231_promoter_reporter/vkorc1_rs9923231_reporter_reference.svg`
- expected alternate construct SVG path after full replay:
  `docs/tutorial/reproducibility/vkorc1_rs9923231_promoter_reporter/vkorc1_rs9923231_reporter_alternate.svg`
- commands:
  [commands.sh](/Users/u005069/.codex/worktrees/47dd/gentle_rs/docs/tutorial/reproducibility/vkorc1_rs9923231_promoter_reporter/commands.sh)
- structured summary:
  [result.json](/Users/u005069/.codex/worktrees/47dd/gentle_rs/docs/tutorial/reproducibility/vkorc1_rs9923231_promoter_reporter/result.json)

## Bench-facing next actions

1. Build the reference and alternate inserts with identical boundaries.
2. Keep the mammalian reporter backbone constant between alleles.
3. Verify insert orientation and junction integrity.
4. Choose one human cell model and a normalization control for the later assay.
5. Treat warfarin exposure as a later experimental condition layered on top of
   the finished promoter-reporter pair.
