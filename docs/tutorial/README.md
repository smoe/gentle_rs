# GENtle Tutorial Guide

This is the canonical entry page for GENtle tutorials.

Machine-readable catalog source:

- [`docs/tutorial/catalog.json`](./catalog.json)
- generated from:
  - [`docs/tutorial/sources/catalog_meta.json`](./sources/catalog_meta.json)
  - [`docs/tutorial/sources/*.json`](./sources/)

GENtle currently has two tutorial families:

- executable tutorials generated from versioned workflow examples
- hand-written walkthroughs and reference guides for GUI-first and agent-facing use

This page keeps those together so users do not need to know the repository
layout first.

## Start Here

If you are new to GENtle and want the most reproducible path first:

- Start with the executable tutorial hub:
  [`docs/tutorial/generated/README.md`](./generated/README.md)

If you want the big-picture map first:

- Use the tutorial landscape overview:
  [`docs/tutorial/landscape_overview.md`](./landscape_overview.md)

If you want the simplest possible PCR walkthrough first:

- Use the simple-PCR selection tutorial:
  [`docs/tutorial/simple_pcr_selection_gui.md`](./simple_pcr_selection_gui.md)

If you want to manually sanity-check the newer protein workflows:

- Start with the transcript-native Protein Expert tutorial:
  [`docs/tutorial/protein_transcript_native_expert_gui.md`](./protein_transcript_native_expert_gui.md)
- Then continue with the reverse-translation and lineage-audit tutorial:
  [`docs/tutorial/protein_reverse_translation_gui.md`](./protein_reverse_translation_gui.md)

If you are GUI-first and want one concrete pharmacogenomic handoff into cloning:

- Use the VKORC1 / rs9923231 PGx-alert-to-mammalian-reporter walkthrough:
  [`docs/tutorial/vkorc1_warfarin_promoter_luciferase_gui.md`](./vkorc1_warfarin_promoter_luciferase_gui.md)

If you want to test the destination-first Gibson specialist end-to-end:

- Use the Gibson specialist testing tutorial:
  [`docs/tutorial/gibson_specialist_testing_gui.md`](./gibson_specialist_testing_gui.md)
- fastest setup path:
  [`docs/tutorial/generated/chapters/15_gibson_specialist_testing_baseline.md`](./generated/chapters/15_gibson_specialist_testing_baseline.md)
  - in the app, `File -> Open Tutorial Project...` -> `Gibson Specialist Starter Project (offline)`
    now opens the starter project and the matching Help/Tutorial guide together
- current guardrail:
  multi-insert Gibson execution currently requires a defined destination
  opening; `existing_termini` remains the single-fragment handoff path

If you want to inspect the arrangement that Gibson apply creates and export one
reusable gel lane set from it:

- Use the Gibson arrangements tutorial:
  [`docs/tutorial/gibson_arrangements_gui.md`](./gibson_arrangements_gui.md)
- fastest setup path:
  [`docs/tutorial/generated/chapters/16_gibson_arrangements_baseline.md`](./generated/chapters/16_gibson_arrangements_baseline.md)
  - in the app, `File -> Open Tutorial Project...` -> `Gibson Arrangements Starter Project (offline)`
    now opens a distinct arrangement-ready Gibson result and lands directly on
    the arrangement guide in Help/Tutorial
- this starts from the deterministic cloned result plus stored arrangement, then
  continues into `Containers`, `Arrangements`, graph inspection, and
  arrangement-level gel export

If you want to continue from that arrangement into a physical carrier and one
README-ready rack figure:

- Use the Gibson physical rack tutorial:
  [`docs/tutorial/gibson_physical_rack_gui.md`](./gibson_physical_rack_gui.md)
- fastest setup path:
  start from the same arrangement-ready tutorial project as above:
  `File -> Open Tutorial Project...` -> `Gibson Arrangements Starter Project (offline)`
- this continues into:
  - linked physical rack inspection
  - pseudo-3D/isometric SVG export
  - carrier labels
  - fabrication SVG
  - OpenSCAD

If you want a sequence-analysis example with screenshots:

- Use the TP73 cDNA-vs-genomic dotplot walkthrough:
  [`docs/tutorial/two_sequence_dotplot_gui.md`](./two_sequence_dotplot_gui.md)

If you want to use GENtle with agents, MCP, or the command line:

- Use the agent interfaces tutorial:
  [`docs/agent_interfaces_tutorial.md`](../agent_interfaces_tutorial.md)

If you want to validate a construct from imported ABI/AB1/SCF evidence in the
GUI today:

- Use the sequencing-confirmation GUI tutorial:
  [`docs/tutorial/sequencing_confirmation_gui.md`](./sequencing_confirmation_gui.md)

If you want the matching shell/CLI parity route:

- Use the sequencing-trace confirmation tutorial:
  [`docs/tutorial/sequencing_confirmation_trace_cli.md`](./sequencing_confirmation_trace_cli.md)

If you want to batch-compare many cDNA `fasta.gz` samples for one target gene:

- Use the RNA-read batch gene-support tutorial:
  [`docs/tutorial/rna_read_batch_gene_support_cli.md`](./rna_read_batch_gene_support_cli.md)

## Confidence Map

Treat the tutorial catalog as three confidence bands:

### Recommended now

- The generated executable tutorial hub:
  [`docs/tutorial/generated/README.md`](./generated/README.md)
- The single-insert Gibson specialist testing path, especially when paired
  with the offline starter project and `gibson preview` parity check.
- Tutorials that point back to shared exported artifacts or deterministic
  workflows rather than screenshots alone.

### Works with caveats

- Manual/hybrid tutorials that are intentionally mapped to engine or CLI
  routes, but still depend on live GUI wording and layout.
- Multi-insert Gibson walkthrough work for preview/review, with the current
  apply guardrail on destination opening still in force.
- Trace-aware sequencing confirmation via both the GUI specialist and
  CLI/shared shell is ready for deterministic local use.
- Current sequencing-confirmation caveat:
  chromatogram review is still variant-focused rather than a full whole-trace
  browser.

### Exploratory / drift-prone

- Purely manual GUI walkthroughs where screenshots and control labels can drift
  faster than generated runtime chapters.
- Topics whose underlying product areas are still explicitly marked as open
  roadmap gaps rather than stable daily-use paths.

## Tutorial Catalog

| Tutorial | Type | Status | Best for | Notes |
| --- | --- | --- | --- | --- |
| [`docs/tutorial/generated/README.md`](./generated/README.md) | Executable tutorial collection | `generated+checked` | Reproducible learning paths, CLI parity, CI-backed examples | Generated from `docs/tutorial/sources/` and executable workflows through the runtime manifest; validated by `tutorial-check`. |
| [`docs/tutorial/simple_pcr_selection_gui.md`](./simple_pcr_selection_gui.md) | GUI walkthrough + shared-engine parity | `manual/hybrid` | Beginner PCR setup from one selected core region | Hand-written short path for core ROI + maximum primer distance + maximum amplicon using the selection context menu and PCR Designer starter block. |
| [`docs/tutorial/protein_transcript_native_expert_gui.md`](./protein_transcript_native_expert_gui.md) | GUI walkthrough + shared-engine parity | `manual/hybrid` | Local transcript-native protein expert validation | Hand-written manual check for `Protein Evidence... -> Open Derived Protein Expert` and derived-only SVG export without external protein evidence. |
| [`docs/tutorial/protein_reverse_translation_gui.md`](./protein_reverse_translation_gui.md) | GUI walkthrough + shared-engine parity | `manual/hybrid` | Protein import, reverse translation, provenance inspection, lineage audit | Hand-written manual check for Ensembl protein import, reverse translation result inspection, and reverse-translation lineage reopen. |
| [`docs/tutorial/vkorc1_warfarin_promoter_luciferase_gui.md`](./vkorc1_warfarin_promoter_luciferase_gui.md) | GUI walkthrough + CLI mapping | `manual/hybrid` | GUI-first pharmacogenomic reporter planning | Hand-written narrative, but intentionally mapped to engine/CLI operations and linked to executable workflow material. |
| [`docs/tutorial/gibson_specialist_testing_gui.md`](./gibson_specialist_testing_gui.md) | GUI walkthrough + CLI parity | `manual/hybrid` | Gibson specialist testing, preview/export parity, contributor verification | Hand-written end-to-end test script for `Patterns -> Gibson...` using local inputs plus `gibson preview`; documents the current multi-insert `defined opening` guardrail. |
| [`docs/tutorial/gibson_arrangements_gui.md`](./gibson_arrangements_gui.md) | GUI walkthrough + CLI parity | `manual/hybrid` | Arrangement reuse, Gibson output inspection, gel-lane planning | Hand-written walkthrough for the arrangement that Gibson apply creates automatically, including singleton output containers and arrangement-level gel export. |
| [`docs/tutorial/gibson_physical_rack_gui.md`](./gibson_physical_rack_gui.md) | GUI walkthrough + CLI parity | `manual/hybrid` | Physical rack export, README-grade isometric figure generation, carrier-label/OpenSCAD handoff | Hand-written walkthrough for taking the arrangement-ready Gibson starter into the linked rack layer and exporting one pseudo-3D/isometric hero SVG plus the other physical carrier projections. |
| [`docs/tutorial/two_sequence_dotplot_gui.md`](./two_sequence_dotplot_gui.md) | GUI walkthrough | `manual` | Screenshot-backed visual analysis tutorial | Good for interface learning; more exposed to UI drift than generated chapters. |
| [`docs/tutorial/sequencing_confirmation_gui.md`](./sequencing_confirmation_gui.md) | GUI walkthrough + shared-engine parity | `manual/hybrid` | GUI-first trace import, intended-edit review, chromatogram inspection | Hand-written walkthrough for `Patterns -> Sequencing Confirmation...`, including raw trace import, baseline-aware variant classification, and lineage reopen. |
| [`docs/tutorial/sequencing_confirmation_trace_cli.md`](./sequencing_confirmation_trace_cli.md) | CLI walkthrough + shell parity | `manual/hybrid` | Imported trace inspection, trace-aware sequencing confirmation, report export | Hand-written local walkthrough for `seq-trace import|list|show` plus `seq-confirm run --trace-id ...` on one deterministic bundled ABI fixture. |
| [`docs/tutorial/rna_read_batch_gene_support_cli.md`](./rna_read_batch_gene_support_cli.md) | CLI walkthrough + shared-shell parity | `manual/hybrid` | Batch cDNA cohort comparison for one target gene | Hand-written batch workflow for many `fa.gz` / `fasta.gz` inputs, report alignment, target-gene abundance summaries, exon-pair co-presence, and assigned-read mean-length export. |
| [`docs/tutorial/landscape_overview.md`](./landscape_overview.md) | Operational reference tutorial | `manual/reference` | Tutorial planning, onboarding sequencing, contributor navigation | Graphical map of current tutorial dependencies, roadmap-grounded additions, and heuristic human-feedback intensity. |
| [`docs/agent_interfaces_tutorial.md`](../agent_interfaces_tutorial.md) | Operational reference tutorial | `manual/reference` | CLI, MCP, in-app agent assistant, external coding agents | Conceptual and operational guide rather than an executable biology walkthrough. |

## Type and Status Labels

Use the labels above as trust/maintenance signals:

- `generated+checked`
  - produced from executable examples and validated by
    `cargo run --bin gentle_examples_docs -- tutorial-check`
  - best choice when reproducibility matters most
- `manual`
  - hand-written walkthrough
  - useful for GUI and screenshot-driven learning
  - higher drift risk when controls or layout change
- `manual/hybrid`
  - hand-written tutorial with explicit links to executable workflows or shared
    engine command paths
  - good bridge between GUI learning and reproducible automation
- `manual/reference`
  - explanatory or operational guide
  - intended for concepts, interfaces, and mental models rather than
    end-to-end biological execution

## Suggested Learning Paths

### Path A: New to GENtle

1. Read [`docs/tutorial/generated/README.md`](./generated/README.md)
2. Run one `core` chapter
3. Open the matching state/artifact outputs

### Path B: GUI-first cloning planning

1. Read [`docs/tutorial/vkorc1_warfarin_promoter_luciferase_gui.md`](./vkorc1_warfarin_promoter_luciferase_gui.md)
2. Use the executable PCR chapter it references
3. Save the resulting project state for later CLI/agent replay

### Path B1: GUI-first simple PCR

1. Read [`docs/tutorial/simple_pcr_selection_gui.md`](./simple_pcr_selection_gui.md)
2. Start from one map selection and the `Simple PCR from selection` context-menu action
3. Use the richer executable PCR chapter afterward if you need batch queueing or painted primer windows

### Path B2: Protein workflow sanity checks

1. Read [`docs/tutorial/protein_transcript_native_expert_gui.md`](./protein_transcript_native_expert_gui.md)
2. Verify transcript-native Protein Expert open/export from one local project
3. Continue with [`docs/tutorial/protein_reverse_translation_gui.md`](./protein_reverse_translation_gui.md)
4. Confirm reverse-translation provenance and lineage reopen behavior

### Path C: Gibson specialist testing

1. Read [`docs/tutorial/gibson_specialist_testing_gui.md`](./gibson_specialist_testing_gui.md)
2. Run the GUI preview/export steps with local inputs
3. Replay the exported plan through `gibson preview` for parity checking

### Path D: Gibson arrangements and gel export

1. Read [`docs/tutorial/gibson_arrangements_gui.md`](./gibson_arrangements_gui.md)
2. Optionally begin from [`docs/tutorial/generated/chapters/16_gibson_arrangements_baseline.md`](./generated/chapters/16_gibson_arrangements_baseline.md)
3. Inspect the prebuilt singleton output containers and reusable arrangement
4. Export one arrangement-level gel with ladders flanking the sample lanes
5. Use the Gibson specialist tutorial separately if you also want to replay the
   cloning step itself

### Path J: Gibson physical rack and README figure export

1. Read [`docs/tutorial/gibson_physical_rack_gui.md`](./gibson_physical_rack_gui.md)
2. Start from the arrangement-ready Gibson tutorial project
3. Open the linked rack draft from the arrangement row
4. Export the pseudo-3D/isometric rack SVG intended for README reuse
5. Optionally export fabrication SVG, carrier labels, and OpenSCAD from the
   same rack state

### Path E: Sequence-analysis and visualization

1. Read [`docs/tutorial/two_sequence_dotplot_gui.md`](./two_sequence_dotplot_gui.md)
2. Compare its screenshots to the live GUI state
3. Re-run the corresponding workflow example if you want a deterministic audit
   path

### Path F: Agents and automation

1. Read [`docs/agent_interfaces_tutorial.md`](../agent_interfaces_tutorial.md)
2. Verify capabilities with `gentle_cli help` or MCP `tools/list`
3. Use executable tutorials as the reproducible backing layer

### Path G: GUI-first sequencing confirmation

1. Read [`docs/tutorial/sequencing_confirmation_gui.md`](./sequencing_confirmation_gui.md)
2. Load the tiny expected construct and baseline sequence
3. Import the bundled `3100.ab1` trace inside `Patterns -> Sequencing Confirmation...`
4. Run confirmation and inspect the intended-edit chromatogram
5. Reopen the stored report from lineage

### Path H: Trace-aware sequencing confirmation via shell/CLI

1. Read [`docs/tutorial/sequencing_confirmation_trace_cli.md`](./sequencing_confirmation_trace_cli.md)
2. Import the bundled `3100.ab1` trace into a dedicated temporary CLI state
3. Run `seq-confirm ... --trace-id ...` and inspect the stored report
4. Export JSON + TSV evidence artifacts for handoff or regression checks

### Path I: Batch RNA-read target-gene cohort comparison

1. Read [`docs/tutorial/rna_read_batch_gene_support_cli.md`](./rna_read_batch_gene_support_cli.md)
2. Start from one prepared target locus and identify the seed `mRNA` feature id
3. Interpret one saved RNA-read report per gzipped cDNA FASTA input
4. Align each saved report with `rna-reads align-report --selection all`
5. Export one target-gene sample sheet for abundance, exon-pair co-presence,
   and mean-length comparison

### Path J: Tutorial planning and dependency triage

1. Read [`docs/tutorial/landscape_overview.md`](./landscape_overview.md)
2. Choose whether you want generated executable depth, GUI walkthrough depth,
   or contributor/reference material first
3. Use the dependency map to decide which starter chapter or companion tutorial
   should come before the page you actually want to follow

## Contributor Notes

When adding a tutorial, decide and document all three of these explicitly:

1. Tutorial type: walkthrough, executable chapter, or reference guide
2. Status: `generated+checked`, `manual`, `manual/hybrid`, or
   `manual/reference`
3. Primary audience: GUI-first, CLI-first, agent-facing, contributor, or
   biology-specific

Contributors should link new tutorials from this page in the same change so
discovery remains centralized.

Machine-readable discovery should also be updated in the same change:

- [`docs/tutorial/catalog.json`](./catalog.json)
- [`docs/tutorial/sources/catalog_meta.json`](./sources/catalog_meta.json)
- relevant [`docs/tutorial/sources/*.json`](./sources/)

## Maintenance Direction

Current chosen direction:

- keep tutorial authoring grounded in per-tutorial source units under
  `docs/tutorial/sources/`
- generate both:
  - discovery catalog: `docs/tutorial/catalog.json`
  - executable runtime manifest: `docs/tutorial/manifest.json`
- keep workflow execution grounded in `docs/examples/workflows/*.json`
- keep manual/tutorial-reference pages, but surface them through this one
  canonical landing page
- keep the distinction between tutorial, recipe, and reference material
  explicit in page metadata and catalog text

Planned next documentation evolution:

- keep reducing drift between source units and hand-written narrative pages
- explore stronger validation for manual tutorial bodies and page-local metadata
