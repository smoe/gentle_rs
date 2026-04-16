# VKORC1 / rs9923231 PGx Alert -> Mammalian Luciferase Reporter Handoff

> Type: `GUI walkthrough + CLI parity`
> Status: `shared-engine baseline`
> Goal: start from one pharmacogenomic alert around warfarin and end with one
> biologically sensible, reproducible promoter-reporter design handoff in
> GENtle.

This tutorial is intentionally modest.

It does **not** claim wet-lab validation, direct drug-response proof, or a
finished screening campaign. It ends when we have one reproducible promoter-SNP
reporter design that a bench-facing colleague could build and compare in a
human-cell luciferase assay.

## What ClawBio Does vs What GENtle Does

ClawBio and GENtle play different roles here.

- ClawBio interprets the pharmacogenomic alert:
  - warfarin sensitivity points us to `VKORC1`
  - `rs9923231` is upstream of `VKORC1`
  - that makes a regulatory follow-up more plausible than a coding follow-up
- GENtle does deterministic build/render work:
  - fetch the locus
  - derive strand-aware promoter windows from transcript TSS geometry
  - summarize whether the SNP falls into promoter context
  - suggest a reporter fragment
  - materialize matched reference/alternate inserts
  - place them into a mammalian luciferase backbone
  - export reviewable artifacts

That distinction matters because the interpretation and the sequence-building
should stay reproducible but conceptually separate.

## Assay Logic

The assay question is:

> Does a `VKORC1` promoter fragment carrying `rs9923231` change reporter output
> in transfected human cells?

The design logic is:

1. retrieve the local `VKORC1` locus around `rs9923231`
2. confirm that the SNP is promoter-proximal on the reverse strand
3. keep the SNP *inside* the promoter fragment, not at the border
4. make two matched inserts:
   - reference allele
   - alternate allele
5. place each insert upstream of luciferase in the same backbone
6. compare reporter output later without changing fragment boundaries

This is a human-cell regulatory assay story. It is **not** a bacterial
expression story. Adenoviral delivery is not the baseline; it only becomes
relevant later if transfection efficiency turns into the bottleneck.

## Backbone Choice

Primary backbone used here:

- `gentle_mammalian_luciferase_backbone_v1`
- file:
  [data/tutorial_inputs/gentle_mammalian_luciferase_backbone_v1.gb](/Users/u005069/.codex/worktrees/47dd/gentle_rs/data/tutorial_inputs/gentle_mammalian_luciferase_backbone_v1.gb)

Why this is the right baseline for human-cell work:

- it is a promoterless mammalian luciferase reporter architecture
- the assay readout is promoter activity, not protein production in bacteria
- it fits transient transfection planning in human cells
- it gives us one pinned local input, so the canonical tutorial no longer
  depends on live GenBank retrieval

Reasonable later alternatives:

- `pGL4.23[luc2/minP]`-class enhancer follow-up if the question shifts from
  promoter sufficiency to enhancer contribution
- NanoLuc-class mammalian reporters if sensitivity becomes the main limitation

## Default Biological Assumptions

This tutorial uses one fixed default interpretation:

- assembly: `GRCh38`
- prepared genome: `Human GRCh38 Ensembl 116`
- SNP: `rs9923231`
- dbSNP fetch span: `+/- 3000 bp`
- selected gene: `VKORC1`
- promoter window default:
  - `1000 bp` upstream of TSS
  - `200 bp` downstream of TSS
- reporter-fragment heuristic:
  - keep `200 bp` downstream of the selected TSS
  - keep `500 bp` beyond the SNP on the biologically upstream side

For the current prepared reference, the recommended fragment is:

- parent context sequence: `vkorc1_rs9923231_context`
- extracted fragment interval: local `2412..3501` (`0-based [from,to)`)
- extracted fragment id: `vkorc1_rs9923231_promoter_fragment`

Why retain `~200 bp` past the TSS:

- it avoids cutting exactly at the annotated TSS
- it preserves a little immediate promoter-proximal / `5' UTR` context
- it still keeps the insert interpretable as a promoter fragment rather than a
  whole-gene fragment

## Prerequisites

1. GENtle desktop application running.
2. Genome catalog available at
   [assets/genomes.json](/Users/u005069/.codex/worktrees/47dd/gentle_rs/assets/genomes.json).
3. The active instance can prepare or already has `Human GRCh38 Ensembl 116`.
4. dbSNP resolution is reachable for `FetchDbSnpRegion`.
5. The local tutorial backbone file exists:
   [data/tutorial_inputs/gentle_mammalian_luciferase_backbone_v1.gb](/Users/u005069/.codex/worktrees/47dd/gentle_rs/data/tutorial_inputs/gentle_mammalian_luciferase_backbone_v1.gb)

## Step 1: Prepare the Reference and Fetch the SNP Locus

GUI:

1. `File -> Prepare Reference Genome...`
2. choose `Human GRCh38 Ensembl 116`
3. prepare it if needed
4. `File -> Fetch GenBank / dbSNP...`
5. rsID = `rs9923231`
6. genome = `Human GRCh38 Ensembl 116`
7. `+/- flank bp` = `3000`
8. output id = `vkorc1_rs9923231_context`
9. click `Fetch Region`

The status footer should now communicate staged progress rather than just
showing a stale warning:

- contacting NCBI Variation
- waiting for response
- parsing placement
- resolving assembly-compatible chromosome
- extracting annotated slice from the prepared genome

CLI parity:

```bash
cargo run --quiet --bin gentle_cli -- \
  op '{"FetchDbSnpRegion":{"rs_id":"rs9923231","genome_id":"Human GRCh38 Ensembl 116","flank_bp":3000,"output_id":"vkorc1_rs9923231_context","annotation_scope":"full","catalog_path":"assets/genomes.json","cache_dir":"data/genomes"}}' \
  --confirm
```

## Step 2: Let GENtle Classify the Variant as Promoter-Proximal

This is the important new part. We no longer start by hand-picking
coordinates. We first let GENtle derive promoter windows from transcript TSS
geometry and summarize the variant context.

GUI:

1. open `vkorc1_rs9923231_context`
2. confirm `Variation`, `Gene`, and `mRNA` are visible
3. open `Engine Ops`
4. run `AnnotatePromoterWindows` with:
   - `input = vkorc1_rs9923231_context`
   - `gene_label = VKORC1`
   - `upstream_bp = 1000`
   - `downstream_bp = 200`
   - `collapse_mode = transcript`
5. run `SummarizeVariantPromoterContext` with:
   - `input = vkorc1_rs9923231_context`
   - `variant_label_or_id = rs9923231`
   - `gene_label = VKORC1`

What to look for:

- `VKORC1` is reverse-strand
- the summary should classify the SNP as promoter-overlapping / promoter
  candidate context
- the signed TSS distance should make sense for the reverse-strand geometry
- suggested assay ids should include
  `allele_paired_promoter_luciferase_reporter`

Shell parity:

```bash
cargo run --quiet --bin gentle_cli -- \
  variant annotate-promoters vkorc1_rs9923231_context \
  --gene-label VKORC1 \
  --upstream-bp 1000 \
  --downstream-bp 200

cargo run --quiet --bin gentle_cli -- \
  variant promoter-context vkorc1_rs9923231_context \
  --variant rs9923231 \
  --gene-label VKORC1 \
  --path docs/tutorial/reproducibility/vkorc1_rs9923231_promoter_reporter/variant_promoter_context.json
```

## Step 3: Let GENtle Suggest the Reporter Fragment

GUI:

1. stay on `vkorc1_rs9923231_context`
2. in `Engine Ops`, run `SuggestPromoterReporterFragments` with:
   - `input = vkorc1_rs9923231_context`
   - `variant_label_or_id = rs9923231`
   - `gene_label = VKORC1`
   - `retain_downstream_from_tss_bp = 200`
   - `retain_upstream_beyond_variant_bp = 500`
   - `max_candidates = 5`
3. inspect the recommended top candidate

Current default recommended interval:

- start = `2412`
- end = `3501`

For this baseline tutorial, we still use `ExtractRegion` to materialize that
candidate, but the geometry is now justified by the engine rather than chosen
manually first.

CLI parity:

```bash
cargo run --quiet --bin gentle_cli -- \
  variant reporter-fragments vkorc1_rs9923231_context \
  --variant rs9923231 \
  --gene-label VKORC1 \
  --retain-downstream-from-tss-bp 200 \
  --retain-upstream-beyond-variant-bp 500 \
  --path docs/tutorial/reproducibility/vkorc1_rs9923231_promoter_reporter/promoter_reporter_candidates.json

cargo run --quiet --bin gentle_cli -- \
  op '{"ExtractRegion":{"input":"vkorc1_rs9923231_context","from":2412,"to":3501,"output_id":"vkorc1_rs9923231_promoter_fragment"}}' \
  --confirm
```

## Step 4: Materialize Matched Reference and Alternate Inserts

Now we turn one promoter fragment into an allele-matched pair without changing
the boundaries.

GUI:

1. in `Engine Ops`, run `MaterializeVariantAllele` twice:
   - reference -> `vkorc1_rs9923231_promoter_reference`
   - alternate -> `vkorc1_rs9923231_promoter_alternate`

CLI parity:

```bash
cargo run --quiet --bin gentle_cli -- \
  variant materialize-allele vkorc1_rs9923231_promoter_fragment \
  --variant rs9923231 \
  --allele reference \
  --output-id vkorc1_rs9923231_promoter_reference

cargo run --quiet --bin gentle_cli -- \
  variant materialize-allele vkorc1_rs9923231_promoter_fragment \
  --variant rs9923231 \
  --allele alternate \
  --output-id vkorc1_rs9923231_promoter_alternate
```

This is a key reproducibility point:

- same fragment geometry
- same backbone later
- only the allele changes

## Step 5: Load the Local Mammalian Reporter Backbone

GUI:

1. `File -> Open Sequence...`
2. open
   [data/tutorial_inputs/gentle_mammalian_luciferase_backbone_v1.gb](/Users/u005069/.codex/worktrees/47dd/gentle_rs/data/tutorial_inputs/gentle_mammalian_luciferase_backbone_v1.gb)
3. use sequence id `gentle_mammalian_luciferase_backbone_v1`

CLI parity:

```bash
cargo run --quiet --bin gentle_cli -- \
  op '{"LoadFile":{"path":"data/tutorial_inputs/gentle_mammalian_luciferase_backbone_v1.gb","as_id":"gentle_mammalian_luciferase_backbone_v1"}}' \
  --confirm
```

## Step 6: Preview the Reporter Pair

There are two equivalent ways to do this now.

### Direct operation path

Run one ligation preview per allele and branch the first preview into stable ids:

```bash
cargo run --quiet --bin gentle_cli -- \
  op '{"Ligation":{"inputs":["vkorc1_rs9923231_promoter_reference","gentle_mammalian_luciferase_backbone_v1"],"circularize_if_possible":false,"protocol":"Blunt","output_prefix":"vkorc1_rs9923231_reporter_reference_assembly","unique":false}}' \
  --confirm

cargo run --quiet --bin gentle_cli -- \
  op '{"Branch":{"input":"vkorc1_rs9923231_reporter_reference_assembly_1","output_id":"vkorc1_rs9923231_reporter_reference"}}' \
  --confirm

cargo run --quiet --bin gentle_cli -- \
  op '{"Ligation":{"inputs":["vkorc1_rs9923231_promoter_alternate","gentle_mammalian_luciferase_backbone_v1"],"circularize_if_possible":false,"protocol":"Blunt","output_prefix":"vkorc1_rs9923231_reporter_alternate_assembly","unique":false}}' \
  --confirm

cargo run --quiet --bin gentle_cli -- \
  op '{"Branch":{"input":"vkorc1_rs9923231_reporter_alternate_assembly_1","output_id":"vkorc1_rs9923231_reporter_alternate"}}' \
  --confirm
```

### Shared macro-template path

The repository now also ships a workflow macro template:

- [assets/cloning_patterns_catalog/reporter/promoter_luciferase/allele_paired_promoter_luciferase_reporter.json](/Users/u005069/.codex/worktrees/47dd/gentle_rs/assets/cloning_patterns_catalog/reporter/promoter_luciferase/allele_paired_promoter_luciferase_reporter.json)

Import and run it through the shared shell:

```bash
cargo run --quiet --bin gentle_cli -- \
  shell 'macros template-import assets/cloning_patterns_catalog'

cargo run --quiet --bin gentle_cli -- \
  shell "macros template-run allele_paired_promoter_luciferase_reporter --bind reference_fragment_seq_id=vkorc1_rs9923231_promoter_reference --bind alternate_fragment_seq_id=vkorc1_rs9923231_promoter_alternate --bind reporter_backbone_seq_id=gentle_mammalian_luciferase_backbone_v1 --bind overlap_bp=20 --bind output_prefix=vkorc1_rs9923231_reporter_pair --transactional"
```

## Step 7: Export Reviewable Artifacts

The current baseline bundle contains:

- promoter-context JSON
- promoter-reporter candidate JSON
- one promoter-context SVG
- planned reference construct SVG target
- planned alternate construct SVG target
- report + result + commands

Workflow replay:

```bash
cargo run --quiet --bin gentle_cli -- \
  workflow @docs/examples/workflows/vkorc1_rs9923231_promoter_luciferase_assay_planning.json
```

Key output files:

- promoter-context SVG:
  [vkorc1_rs9923231_promoter_context.svg](/Users/u005069/.codex/worktrees/47dd/gentle_rs/docs/tutorial/reproducibility/vkorc1_rs9923231_promoter_reporter/vkorc1_rs9923231_promoter_context.svg)
- expected reference construct SVG path after full replay:
  `docs/tutorial/reproducibility/vkorc1_rs9923231_promoter_reporter/vkorc1_rs9923231_reporter_reference.svg`
- expected alternate construct SVG path after full replay:
  `docs/tutorial/reproducibility/vkorc1_rs9923231_promoter_reporter/vkorc1_rs9923231_reporter_alternate.svg`

At the moment, the promoter-context half of the workflow is committed and
replay-verified in-tree. The paired construct SVG targets are still owned by
the same shared workflow, but they need one clean end-to-end smoke run before
they can be committed alongside the rest of the bundle.

## Reproducibility Bundle

The handoff bundle for this tutorial lives in:

- [report.md](/Users/u005069/.codex/worktrees/47dd/gentle_rs/docs/tutorial/reproducibility/vkorc1_rs9923231_promoter_reporter/report.md)
- [result.json](/Users/u005069/.codex/worktrees/47dd/gentle_rs/docs/tutorial/reproducibility/vkorc1_rs9923231_promoter_reporter/result.json)
- [commands.sh](/Users/u005069/.codex/worktrees/47dd/gentle_rs/docs/tutorial/reproducibility/vkorc1_rs9923231_promoter_reporter/commands.sh)

## Bench-Facing Next Actions

1. Build the reference and alternate promoter fragments with identical
   boundaries.
2. Keep the mammalian reporter backbone constant between the two constructs.
3. Confirm junctions and insert orientation before comparing reporter output.
4. Choose one human cell model and one normalization strategy for the later
   assay run.
5. Treat warfarin exposure as a later experimental condition, not as the first
   claim of this design handoff.
