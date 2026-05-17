# GENtle Tutorial Guide

This is the canonical entry page for GENtle tutorials.

Machine-readable catalog source:

- [`docs/tutorial/catalog.json`](./catalog.json)
- generated from:
  - [`docs/tutorial/sources/catalog_meta.json`](./sources/catalog_meta.json)
  - [`docs/tutorial/sources/*.json`](./sources/)

GENtle has two tutorial tracks:

- guided walkthroughs: hand-written pages for first-time learning, GUI context,
  and topic-level orientation
- executable reference chapters: machine-generated, numbered chapters for
  reproducible workflow replay and artifact auditing

Use the guided walkthroughs first when you want a human teaching path. Use the
executable reference chapters when you want to run or inspect the exact
workflow operations behind a topic.

## Start Here

If you are new to GENtle:

- Start with the guided walkthroughs below, especially
  [`docs/tutorial/simple_pcr_selection_gui.md`](./simple_pcr_selection_gui.md)
  for a compact GUI-first path.

If you want the big-picture map first:

- Use the tutorial landscape overview:
  [`docs/tutorial/landscape_overview.md`](./landscape_overview.md)

If you want machine-checked reproducibility first:

- Use the executable reference hub:
  [`docs/tutorial/generated/README.md`](./generated/README.md)

## Guided Walkthroughs (hand-written)

These pages are written as tutorials rather than indexes. They are the best
entry point when you want context, screenshots, GUI orientation, or a
bench-user narrative before running the exact workflow.

### Orientation And Interfaces

- [`docs/tutorial/landscape_overview.md`](./landscape_overview.md) - Use this when: you want the tutorial map, dependencies, and suggested learning order before choosing a path.
- [`docs/agent_interfaces_tutorial.md`](../agent_interfaces_tutorial.md) - Use this when: you want to operate GENtle through CLI, MCP, the in-app assistant, or external coding agents.

### PCR, qPCR, And Direct Sequence Inspection

- [`docs/tutorial/simple_pcr_selection_gui.md`](./simple_pcr_selection_gui.md) - Use this when: you want the shortest GUI path from one selected core region to PCR primer candidates.
- [`docs/tutorial/qpcr_exon_junctions_gui.md`](./qpcr_exon_junctions_gui.md) - Use this when: you want qPCR assays from transcript or exon-junction context.
- [`docs/tutorial/stateless_sequence_inspection_gui_cli.md`](./stateless_sequence_inspection_gui_cli.md) - Use this when: you want to inspect one pasted or local sequence across GUI, CLI, and ClawBio without building a project.
- [`docs/tutorial/tfbs_similarity_ranking_gui.md`](./tfbs_similarity_ranking_gui.md) - Use this when: you want to sign off anchor-vs-candidate TFBS similarity ranking in the GUI.
- [`docs/tutorial/two_sequence_dotplot_gui.md`](./two_sequence_dotplot_gui.md) - Use this when: you want a screenshot-backed TP73 cDNA-vs-genomic dotplot walkthrough.

### Cloning, Reporter Design, And External Handoffs

- [`docs/tutorial/vkorc1_variant_followup_expert_gui.md`](./vkorc1_variant_followup_expert_gui.md) - Use this when: you want the fastest GUI path from one SNP to one portable promoter-design bundle.
- [`docs/tutorial/vkorc1_warfarin_promoter_luciferase_gui.md`](./vkorc1_warfarin_promoter_luciferase_gui.md) - Use this when: you want the longer VKORC1/rs9923231 PGx-alert-to-mammalian-reporter story with CLI mapping.
- [`docs/tutorial/reporter_construct_handoff_cli.md`](./reporter_construct_handoff_cli.md) - Use this when: you want a saved promoter-candidate report turned into a read-only macro-readiness plan.
- [`docs/tutorial/metabion_external_service_handoff_gui_cli.md`](./metabion_external_service_handoff_gui_cli.md) - Use this when: you want to rehearse a safe Metabion oligo or m-block quote handoff.
- [`docs/tutorial/geneart_external_service_handoff_gui_cli.md`](./geneart_external_service_handoff_gui_cli.md) - Use this when: you want to rehearse a safe GeneArt cloned-gene or protein-expression quote handoff.

### Gibson Assembly And Physical Layout

- [`docs/tutorial/gibson_specialist_testing_gui.md`](./gibson_specialist_testing_gui.md) - Use this when: you want to test the destination-first Gibson specialist end to end.
- [`docs/tutorial/gibson_arrangements_gui.md`](./gibson_arrangements_gui.md) - Use this when: you want to inspect the arrangement created by Gibson apply and export reusable gel lanes.
- [`docs/tutorial/gibson_physical_rack_gui.md`](./gibson_physical_rack_gui.md) - Use this when: you want to continue from a Gibson arrangement into rack, label, fabrication SVG, and OpenSCAD exports.

### Protein, RNA, And Projection Audits

- [`docs/tutorial/protein_transcript_native_expert_gui.md`](./protein_transcript_native_expert_gui.md) - Use this when: you want a local transcript-native Protein Expert sanity check.
- [`docs/tutorial/protein_reverse_translation_gui.md`](./protein_reverse_translation_gui.md) - Use this when: you want protein import, reverse translation, provenance inspection, and lineage audit.
- [`docs/tutorial/rna_read_batch_gene_support_cli.md`](./rna_read_batch_gene_support_cli.md) - Use this when: you want to batch-compare many cDNA `fasta.gz` samples for one target gene.
- [`docs/tutorial/tp73_uniprot_projection_audit_cli.md`](./tp73_uniprot_projection_audit_cli.md) - Use this when: you want a CLI-first TP73 protein/transcript audit against UniProt and Ensembl.

### Sequencing Confirmation

- [`docs/tutorial/sequencing_confirmation_gui.md`](./sequencing_confirmation_gui.md) - Use this when: you want GUI-first trace import, intended-edit review, and chromatogram inspection.
- [`docs/tutorial/sequencing_confirmation_trace_cli.md`](./sequencing_confirmation_trace_cli.md) - Use this when: you want the matching shell/CLI route for imported trace inspection and report export.

## Executable Reference Chapters (machine-generated)

[`docs/tutorial/generated/README.md`](./generated/README.md) lists the numbered
chapters generated from `docs/examples/workflows/`. Use these when you want to
run a workflow non-interactively or audit what each operation produces.

## Confidence Map

Treat the tutorial catalog as three confidence bands:

### Recommended now

- Guided walkthroughs that point back to shared exported artifacts or
  deterministic workflows.
- The generated executable reference hub when reproducibility or artifact audit
  matters most:
  [`docs/tutorial/generated/README.md`](./generated/README.md)
- The single-insert Gibson specialist testing path, especially when paired
  with the offline starter project and `gibson preview` parity check.

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

## Machine-Readable Catalog

The discoverable tutorial catalog remains machine-readable and is generated
from `docs/tutorial/sources/`:

- [`docs/tutorial/catalog.json`](./catalog.json) lists both hand-written
  walkthroughs and the generated hub for application discovery.
- [`docs/tutorial/manifest.json`](./manifest.json) lists executable reference
  chapters and their workflow examples for tutorial project generation.

## Type and Status Labels

Use the labels above as trust/maintenance signals:

- `generated+checked`
  - produced from executable examples and validated by
    `cargo run --bin gentle_examples_docs -- tutorial-check`
  - best choice when reproducibility matters most
- `generated+checked/human-pending`
  - generated and machine-checked, but listed separately for explicit human
    functional confirmation
  - use this queue when you want to pick one tutorial and validate that the
    prose, GUI affordances, and expected artifacts still line up
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

1. Read [`docs/tutorial/simple_pcr_selection_gui.md`](./simple_pcr_selection_gui.md)
2. Use [`docs/tutorial/landscape_overview.md`](./landscape_overview.md) to choose your next topic.
3. Open [`docs/tutorial/generated/README.md`](./generated/README.md) when you want the matching executable reference chapter.

### Path B: GUI-first cloning planning

1. Read [`docs/tutorial/vkorc1_variant_followup_expert_gui.md`](./vkorc1_variant_followup_expert_gui.md)
2. Continue with [`docs/tutorial/vkorc1_warfarin_promoter_luciferase_gui.md`](./vkorc1_warfarin_promoter_luciferase_gui.md) if you want the fuller background and CLI mapping
3. Use [`docs/tutorial/reporter_construct_handoff_cli.md`](./reporter_construct_handoff_cli.md) to turn the saved promoter-candidate report into a reviewed macro-readiness plan
4. Save the resulting project state and exported bundle for later CLI/agent replay

### Path B1: GUI-first simple PCR

1. Open `File -> Open Tutorial Project...` -> `Core (9)` -> `18. Simple PCR From a Selected Core Region`
2. Follow [`docs/tutorial/simple_pcr_selection_gui.md`](./simple_pcr_selection_gui.md) from the loaded TP73 starter project
3. Use the richer executable PCR chapter afterward if you need batch queueing or painted primer windows

### Path B3: TFBS similarity sign-off

1. Read [`docs/tutorial/tfbs_similarity_ranking_gui.md`](./tfbs_similarity_ranking_gui.md)
2. Load the tiny synthetic FASTA and verify the GUI `TFBS similarity` path end-to-end
3. Replay the matching offline workflow and ClawBio request if you want parity before release sign-off

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

### Path F2: External-service handoff rehearsal

1. Read [`docs/tutorial/metabion_external_service_handoff_gui_cli.md`](./metabion_external_service_handoff_gui_cli.md)
2. Read [`docs/tutorial/geneart_external_service_handoff_gui_cli.md`](./geneart_external_service_handoff_gui_cli.md)
3. Run `services providers doctor` and `services providers list`
4. Preflight the bundled oligo, m-block, cloned-gene, and protein-expression
   request JSON examples
5. Generate quote-handoff bundles and inspect normalized JSON/CSV, email draft,
   WOP checklist, warnings, and required follow-up
6. Repeat the same review in `Services -> External Services...`

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

### Path J: TP73 UniProt/Ensembl projection audit

1. Optionally begin from the executable companion:
   [`docs/tutorial/generated/chapters/19_tp73_uniprot_projection_audit_cli.md`](./generated/chapters/19_tp73_uniprot_projection_audit_cli.md)
2. Read [`docs/tutorial/tp73_uniprot_projection_audit_cli.md`](./tp73_uniprot_projection_audit_cli.md)
3. Run the integrated `uniprot audit-projection ...` path on the stored TP73 projection
4. Inspect the saved audit report and local unsent maintainer-email draft
5. Rebuild the same reasoning from `resolve-ensembl-links`, `transcript-accounting`,
   `compare-ensembl-exons`, and `compare-ensembl-peptide`
6. Run `uniprot audit-parity ...` to compare the integrated Rust audit against
   the public primitive composition

### Path K: Tutorial planning and dependency triage

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
