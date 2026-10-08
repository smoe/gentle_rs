# VKORC1 / rs9923231 PGx Alert -> Mammalian Luciferase Reporter Handoff

> Type: `GUI walkthrough + CLI parity`
> Status: `shared-engine baseline`
> Goal: start from one pharmacogenomic alert around warfarin and end with one
> biologically sensible, reproducible promoter-reporter design handoff in
> GENtle.

If you want the shortest GUI path first, use:

- [docs/tutorial/08-02_vkorc1_variant_followup_expert_gui.md](./08-02_vkorc1_variant_followup_expert_gui.md)

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
  [`data/tutorial_inputs/gentle_mammalian_luciferase_backbone_v1.gb`](../../data/tutorial_inputs/gentle_mammalian_luciferase_backbone_v1.gb)

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

### Offline Allele-Choice Guard

The generated 08.04 companion is a **synthetic 20-base GUI exercise**, not the
human locus below. Its incomplete starter and independent C/T oracle test the
explicit-allele part of step 4 without network access. Open the fragment,
expand `variation`, select its exact feature, and use `Open Promoter Design`.
Scroll to `Alternate base`, enter `T`, then scroll to
`Make reference/alternate inserts`. Persisted source and output bases and
annotations must match the separate oracle; the sole base difference is
zero-based position 6. A bare alternate must be refused, with no partial pair.

The [fixture provenance](../examples/assets/allele_pair_guard/README.md) states
exactly how that non-biological input was constructed. This guard does not
certify online genome retrieval, the full reporter handoff, a real model's
planning, native package behavior, or human scientific approval. Authored
contracts and parser checks are not an executed native pass.

### Online Tutorial

1. GENtle desktop application running.
2. Genome catalog available at
   [`assets/genomes.json`](../../assets/genomes.json).
3. The active instance can prepare or already has `Human GRCh38 Ensembl 116`.
4. dbSNP resolution is reachable for `FetchDbSnpRegion`.
5. The local tutorial backbone file exists:
   [`data/tutorial_inputs/gentle_mammalian_luciferase_backbone_v1.gb`](../../data/tutorial_inputs/gentle_mammalian_luciferase_backbone_v1.gb)

### Review-First Shared-Shell Commands

These are the exact parser-backed commands an Agent Assistant may propose for
human review on the **online** project. They are not executed by tutorial
validation. Never infer `T` from a bare `alternate`, and never recode the
natural assay target as part of this workflow.

```text
variant annotate-promoters vkorc1_rs9923231_context --gene-label VKORC1 --upstream-bp 1000 --downstream-bp 200
variant promoter-context vkorc1_rs9923231_context --variant rs9923231 --gene-label VKORC1
variant reporter-fragments vkorc1_rs9923231_context --variant rs9923231 --gene-label VKORC1 --retain-downstream-from-tss-bp 200 --retain-upstream-beyond-variant-bp 500
variant materialize-allele vkorc1_rs9923231_promoter_fragment --variant rs9923231 --allele reference --output-id vkorc1_rs9923231_promoter_reference
variant materialize-allele vkorc1_rs9923231_promoter_fragment --variant rs9923231 --allele alternate --alternate-base T --output-id vkorc1_rs9923231_promoter_alternate
```

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

Before proceeding, inspect the fetched marker's `dbsnp_reference_check` and
`dbsnp_assembly_check`: both must be `match`. A fallback warning means the
available placement was not bound to the selected assembly family. Matching a
single base does not fix that coordinate mismatch; retain it for inspection,
but fetch compatible placement evidence before making allele inserts.

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
3. select the `variation` marker for `rs9923231`
4. open `Promoter design`
   - from the description pane, map context menu, or feature-tree context menu
5. in the dedicated window:
   - leave `variant_label_or_id = rs9923231`
   - keep `gene_label = VKORC1`
   - keep promoter window defaults `1000 / 200`
6. click `Annotate promoter windows`
7. click `Summarize promoter context`

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

1. stay in `Promoter design`
2. keep the default fragment heuristic:
   - `retain_downstream_from_tss_bp = 200`
   - `retain_upstream_beyond_variant_bp = 500`
   - `max_candidates = 5`
3. click `Propose reporter fragment`
4. inspect the recommended top candidate
5. click `Extract recommended fragment`

Current default recommended interval:

- start = `2412`
- end = `3501`

For this baseline tutorial, the GUI expert still calls the shared
`ExtractRegion` operation under the hood, but the geometry is now justified by
the engine rather than chosen manually first.

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

1. stay in `Promoter design`
2. confirm the default fragment id:
   - `vkorc1_rs9923231_promoter_fragment`
3. set `Alternate base` to `T`
4. click `Make reference/alternate inserts`
5. expect:
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
  --alternate-base T \
  --output-id vkorc1_rs9923231_promoter_alternate
```

Glen's retained WIP report records genomic-forward `C` as the GRCh38 reference
and `A,G,T` as alternate candidates for this refSNP. This tutorial explicitly
selects `T` from that reported set, not from clinical allele nomenclature or
the gene strand. GENtle therefore refuses an
ambiguous bare `alternate` request and requires the explicit base when more
than one candidate is present.

The GUI checks both inserts before creating either. If `Alternate base` is
empty or invalid, correct the choice and retry: the failed validation creates
no reference copy and does not change the requested output names. The two
CLI commands above are separate operations, so review both before running them.

These are alleles on the ascending genomic-forward slice, not on the
reverse-strand transcript. Recheck the fetched marker before proceeding;
the tutorial is not permission to substitute `T` into a different assembly
or a reverse-complemented sequence with stale allele qualifiers.

This is a key reproducibility point:

- same fragment geometry
- same backbone later
- only the allele changes

## Step 5: Use the Local Mammalian Reporter Backbone

GUI:

1. stay in `Promoter design`
2. confirm the pinned local backbone fields:
   - sequence id `gentle_mammalian_luciferase_backbone_v1`
   - path
     [`data/tutorial_inputs/gentle_mammalian_luciferase_backbone_v1.gb`](../../data/tutorial_inputs/gentle_mammalian_luciferase_backbone_v1.gb)
3. the GUI expert loads this backbone automatically on first preview if it is
   not already present in the current project state

CLI parity:

```bash
cargo run --quiet --bin gentle_cli -- \
  op '{"LoadFile":{"path":"data/tutorial_inputs/gentle_mammalian_luciferase_backbone_v1.gb","as_id":"gentle_mammalian_luciferase_backbone_v1"}}' \
  --confirm
```

## Step 6: Preview the Reporter Pair

There are two equivalent ways to do this now.

GUI:

1. stay in `Promoter design`
2. click `Preview luciferase pair`
3. expect two derived preview ids:
   - `vkorc1_rs9923231_reporter_reference`
   - `vkorc1_rs9923231_reporter_alternate`
4. these are still preview/build artifacts, not a wet-lab validation claim
5. the current GUI expert uses the same shared `LoadFile`, `Ligation`, and
   `Branch` operations described below

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

- [`assets/cloning_patterns_catalog/reporter/promoter_luciferase/allele_paired_promoter_luciferase_reporter.json`](../../assets/cloning_patterns_catalog/reporter/promoter_luciferase/allele_paired_promoter_luciferase_reporter.json)

Before running the macro, generate the read-only reporter construct handoff.
This consumes the saved promoter-fragment candidate report, keeps the reporter
recommendation offline, names the exact macro template, and reports which
fragment/backbone inputs are ready, derivable, or still need to be loaded. If
the source refSNP is multiallelic, review the generated materialization command
and add the same explicit `alternate_allele` choice used above:

For a single alternate, a handoff can include the exact base only when its
source marker is loaded and passes the shared validator. A fresh stateless
run or ambiguous/unverified source produces `Review required:` with no chosen
alternate. Candidate geometry and ready macro ports do not choose an allele.

```bash
cargo run --quiet --bin gentle_cli -- \
  reporters plan-handoff docs/tutorial/reproducibility/vkorc1_rs9923231_promoter_reporter/promoter_reporter_candidates.json \
  --backbone-seq-id gentle_mammalian_luciferase_backbone_v1 \
  --backbone-path data/tutorial_inputs/gentle_mammalian_luciferase_backbone_v1.gb \
  --reference-fragment-seq-id vkorc1_rs9923231_promoter_reference \
  --alternate-fragment-seq-id vkorc1_rs9923231_promoter_alternate \
  --output docs/tutorial/reproducibility/vkorc1_rs9923231_promoter_reporter/reporter_construct_handoff.json
```

For a slower walkthrough of how to inspect that JSON plan and use its
validate-only command safely, see
[`docs/tutorial/08-05_reporter_construct_handoff_cli.md`](./08-05_reporter_construct_handoff_cli.md).

Import and run it through the shared shell:

```bash
cargo run --quiet --bin gentle_cli -- \
  shell 'macros template-import assets/cloning_patterns_catalog'

cargo run --quiet --bin gentle_cli -- \
  shell "macros template-run allele_paired_promoter_luciferase_reporter --bind reference_fragment_seq_id=vkorc1_rs9923231_promoter_reference --bind alternate_fragment_seq_id=vkorc1_rs9923231_promoter_alternate --bind reporter_backbone_seq_id=gentle_mammalian_luciferase_backbone_v1 --bind overlap_bp=20 --bind output_prefix=vkorc1_rs9923231_reporter_pair --transactional"
```

## Agent Assistant Parity

The inner Agent Assistant should inspect live project state before proposing
mutations. A useful prompt after Step 1 is:

> For `vkorc1_rs9923231_context`, verify the `rs9923231` marker alleles, derive
> the VKORC1 promoter reporter fragment, and propose matched reference and
> genomic-forward T constructs. Show every mutating command and wait for my
> approval.

The review must expose `vcf_ref=C`, `vcf_alt=A,G,T`, explain why this tutorial
chooses `T`, and propose the same shared-shell operations used above. In
particular, the alternate command must contain `--alternate-base T`; an agent
must not silently choose one allele from a multiallelic marker.

## Step 7: Export Reviewable Artifacts

GUI:

1. stay in `Promoter design`
2. click `Export handoff bundle`
3. choose a parent folder
4. GENtle creates one bundle directory containing:
   - promoter-context JSON
   - promoter-candidate JSON
   - promoter-context SVG
   - reference construct SVG
   - alternate construct SVG
   - `report.md`
   - `result.json`
   - `commands.sh`

The current baseline bundle format contains:

- promoter-context JSON
- promoter-reporter candidate JSON
- one promoter-context SVG
- one reference construct SVG
- one alternate construct SVG
- report + result + commands

Workflow replay:

```bash
cargo run --quiet --bin gentle_cli -- \
  workflow @docs/examples/workflows/vkorc1_rs9923231_promoter_luciferase_assay_planning.json
```

Key output files:

- promoter-context SVG:
  [`vkorc1_rs9923231_promoter_context.svg`](./reproducibility/vkorc1_rs9923231_promoter_reporter/vkorc1_rs9923231_promoter_context.svg)
- reference construct SVG path:
  `vkorc1_rs9923231_reporter_reference.svg`
- alternate construct SVG path:
  `vkorc1_rs9923231_reporter_alternate.svg`

The GUI expert now writes this bundle directly. The workflow is the parity and
replay path for the same story; its dbSNP lookup and first genome preparation
are online, while the pinned backbone and later build steps are local.

## Reproducibility Bundle

The handoff bundle for this tutorial lives in:

- [`report.md`](./reproducibility/vkorc1_rs9923231_promoter_reporter/report.md)
- [`result.json`](./reproducibility/vkorc1_rs9923231_promoter_reporter/result.json)
- [`commands.sh`](./reproducibility/vkorc1_rs9923231_promoter_reporter/commands.sh)

The retained SVGs and JSON are **historical WIP previews**, not an accepted
replay at the current GENtle revision. See the
[artifact provenance and pending checks](./reproducibility/vkorc1_rs9923231_promoter_reporter/README.md).
The synthetic backbone is a map/planning placeholder, not an exact functional
luciferase plasmid. In particular, selecting blunt-ligation output `_1` and
drawing it circularly does not verify promoter-to-luciferase orientation,
insertion-site geometry or molecular topology. Those checks remain required
before treating the preview as an assay construct.

## Bench-Facing Next Actions

1. Build the reference and alternate promoter fragments with identical
   boundaries.
2. Keep the mammalian reporter backbone constant between the two constructs.
3. Confirm junctions and insert orientation before comparing reporter output.
4. Choose one human cell model and one normalization strategy for the later
   assay run.
5. Treat warfarin exposure as a later experimental condition, not as the first
   claim of this design handoff.

## Feedback

If this tutorial is confusing, execution-stale, biologically suspect, or missing a useful figure, please open the matching tutorial issue template and include the context copied from GENtle Help -> Tutorial -> Copy Feedback Context.

- Tutorial title:
- Tutorial/chapter id:
- Step reached:
- Expected vs. actual:
- Interface used: GUI / CLI / Agent Assistant / ClawBio

Paste the Tutorial feedback context here:

```text

```
