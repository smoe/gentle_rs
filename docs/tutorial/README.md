# GENtle Tutorial Guide

This is the canonical entry page for GENtle tutorials.

Machine-readable catalog source:

- [`docs/tutorial/catalog.json`](./catalog.json)
- generated from:
  - [`docs/tutorial/sources/catalog_meta.json`](./sources/catalog_meta.json)
  - [`docs/tutorial/sources/*.json`](./sources/)
- review metadata:
  - [`docs/tutorial/review_manifest.json`](./review_manifest.json)
- presentation-overhaul audit:
  - [`docs/tutorial/coverage_audit.md`](./coverage_audit.md)

`review_manifest.json` is hand-maintained tutorial review metadata, not test
fixture data from an external source. It is deterministically recreated from
the current tutorial ids in `docs/tutorial/catalog.json` and
`docs/tutorial/manifest.json`, then filled with review dates as tutorials are
read by Codex or a human reviewer. The tutorial checker treats missing, stale,
or unknown review entries as warnings so review provenance cannot mask
execution drift.

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

## Tutorials By Content Group

Tutorial numbers now describe the topic group first and the learning position second. Guided walkthroughs and executable reference chapters are shown together inside each content group, while the track marker tells you whether a link opens a hand-written guide, a generated executable chapter, or a reference hub.

### Getting Started & Interfaces

- `01.01` [GENtle Agent Assistant and Agent Interfaces Tutorial](../agent_interfaces_tutorial.md) - reference; status `manual/reference`. Use this when: Practical guide for the in-app Agent Assistant, provider quick starts, reviewed shared-shell suggestions, CLI/shared shell, MCP, and external coding agents.
- `01.02` [Contribute to GENtle development](./generated/chapters/08_contribute_to_gentle_development.md) - executable reference; status `generated+checked/human-pending`. Use this when: Machine-checked generated chapter; listed in Tutorials for explicit human functional confirmation.
- `reference` [GENtle Tutorial (Generated)](./generated/README.md) - reference hub; status `generated+checked`. Use this when: Generated from docs/tutorial/sources plus executable workflows; validated by tutorial-check.
- `reference` [GENtle Tutorial Landscape Overview](./landscape_overview.md) - reference; status `manual/reference`. Use this when: Graphical overview of current tutorial dependencies, roadmap-grounded additions, and heuristic human-feedback intensity.

### Sequence Basics & Lineage

- `02.01` [Load FASTA, branch, and reverse-complement](./generated/chapters/01_load_branch_reverse_complement_pgex_fasta.md) - executable reference; status `generated+checked/human-pending`. Use this when: Machine-checked generated chapter; listed in Tutorials for explicit human functional confirmation.
- `02.02` [Stateless Sequence Inspection Tutorial](./stateless_sequence_inspection_gui_cli.md) - guided GUI/CLI; status `manual/hybrid`. Use this when: Parity walkthrough for direct restriction-site, TFBS hit, and TFBS score-track inspection from pasted or local DNA without first creating project state.
- `02.03` [Retrieve a cDNA and Compare It Against the Genomic Sequence (Dotplot Tutorial)](./two_sequence_dotplot_gui.md) - guided GUI; status `manual`. Use this when: Screenshot-backed manual tutorial for cDNA-vs-genomic comparison.
- `02.04` [Compare TP73 cDNA against TP73 genomic context via dotplot (online)](./generated/chapters/14_tp73_cdna_genomic_dotplot_online.md) - executable reference; status `generated+checked/human-pending`. Use this when: Machine-checked generated chapter; listed in Tutorials for explicit human functional confirmation. Requires explicit online setup for full execution.

### Cloning & Assembly

- `03.01` [Load pGEX and digest with BamHI/EcoRI](./generated/chapters/03_load_and_digest_pgex.md) - executable reference; status `generated+checked/human-pending`. Use this when: Machine-checked generated chapter; listed in Tutorials for explicit human functional confirmation.
- `03.02` [Digest -> Ligation -> ExtractRegion minimal slice](./generated/chapters/06_digest_ligation_extract_region_minimal.md) - executable reference; status `generated+checked/human-pending`. Use this when: Machine-checked generated chapter; listed in Tutorials for explicit human functional confirmation.
- `03.03` [Gibson two-fragment overlap planning baseline](./generated/chapters/04_gibson_two_fragment_overlap_preview.md) - executable reference; status `generated+checked/human-pending`. Use this when: Machine-checked generated chapter; listed in Tutorials for explicit human functional confirmation.
- `03.04` [Gibson Specialist Starter Project (offline)](./generated/chapters/15_gibson_specialist_testing_baseline.md) - executable reference; status `generated+checked/human-pending`. Use this when: Machine-checked generated chapter; listed in Tutorials for explicit human functional confirmation.
- `03.05` [Gibson Specialist Testing Tutorial](./gibson_specialist_testing_gui.md) - guided GUI; status `manual/hybrid`. Use this when: Hands-on test script for the destination-first Gibson specialist with shared CLI parity.
- `03.06` [Gibson Arrangements Starter Project (offline)](./generated/chapters/16_gibson_arrangements_baseline.md) - executable reference; status `generated+checked/human-pending`. Use this when: Machine-checked generated chapter; listed in Tutorials for explicit human functional confirmation.
- `03.07` [Gibson Arrangements Tutorial](./gibson_arrangements_gui.md) - guided GUI; status `manual/hybrid`. Use this when: Hands-on walkthrough for Gibson-created arrangements, lane inspection, and arrangement-level gel export.
- `03.08` [Gibson Physical Rack Tutorial](./gibson_physical_rack_gui.md) - guided GUI; status `manual/hybrid`. Use this when: Hands-on walkthrough for turning one arrangement-ready Gibson result into a linked physical rack, README-grade isometric SVG, carrier labels, and OpenSCAD export.

### Primers, PCR & qPCR

- `04.01` [Simple PCR From a Selected Core Region](./simple_pcr_selection_gui.md) - guided GUI; status `manual/hybrid`. Use this when: Beginner-friendly selection-first PCR walkthrough centered on core ROI, maximum primer distance from the core, and maximum amplicon length.
- `04.02` [Selection-first PCR batch primer design (offline)](./generated/chapters/13_pcr_selection_batch_primer_pairs_offline.md) - executable reference; status `generated+checked/human-pending`. Use this when: Machine-checked generated chapter; listed in Tutorials for explicit human functional confirmation.
- `04.03` [qPCR Across Exon Junctions](./qpcr_exon_junctions_gui.md) - guided GUI; status `manual/hybrid`. Use this when: GUI-first walkthrough for seeding shared-transcript or transcript-specific qPCR from Splicing Expert, including junction-only transcript evidence and shared-shell parity commands.
- `04.04` [Guide practical filtering and oligo generation](./generated/chapters/05_guides_filter_and_generate_oligos.md) - executable reference; status `generated+checked/human-pending`. Use this when: Machine-checked generated chapter; listed in Tutorials for explicit human functional confirmation.
- `04.05` [Guide oligo export (CSV + protocol)](./generated/chapters/07_guides_export_csv_and_protocol.md) - executable reference; status `generated+checked/human-pending`. Use this when: Machine-checked generated chapter; listed in Tutorials for explicit human functional confirmation.

### Genome Context & Coordinates

- `05.01` [Find and extend the right genomic target (local catalog)](./generated/chapters/02_find_and_extend_genomic_target_local_catalog.md) - executable reference; status `generated+checked/human-pending`. Use this when: Machine-checked generated chapter; listed in Tutorials for explicit human functional confirmation.
- `05.02` [Prepare a reference genome cache (online)](./generated/chapters/09_prepare_reference_genome_online.md) - executable reference; status `generated+checked/human-pending`. Use this when: Machine-checked generated chapter; listed in Tutorials for explicit human functional confirmation. Requires explicit online setup for full execution.
- `05.03` [Retrieve TP63 and extend the displayed region by +/-2 kb (online)](./generated/chapters/11_tp63_anchor_extension_online.md) - executable reference; status `generated+checked/human-pending`. Use this when: Machine-checked generated chapter; listed in Tutorials for explicit human functional confirmation. Requires explicit online setup for full execution.

### Transcript, Protein & Projection

- `06.01` [Reverse Translate an Imported Protein and Audit the Result](./protein_reverse_translation_gui.md) - guided GUI; status `manual/hybrid`. Use this when: GUI-first manual check for Ensembl protein import, reverse translation, provenance inspection, and lineage reopen.
- `06.02` [Transcript-Native Protein Expert Sanity Check](./protein_transcript_native_expert_gui.md) - guided GUI; status `manual/hybrid`. Use this when: Local-first manual check for transcript-native Protein Expert opening, provenance visibility, and derived-only SVG export.
- `06.03` [TP53 isoform architecture expert panel (online)](./generated/chapters/10_tp53_isoform_architecture_online.md) - executable reference; status `generated+checked/human-pending`. Use this when: Machine-checked generated chapter; listed in Tutorials for explicit human functional confirmation. Requires explicit online setup for full execution.
- `06.04` [TP53 UniProt domain mapping and feature-coding DNA query (online)](./generated/chapters/17_tp53_uniprot_projection_online.md) - executable reference; status `generated+checked/human-pending`. Use this when: Machine-checked generated chapter; listed in Tutorials for explicit human functional confirmation. Requires explicit online setup for full execution.
- `06.05` [Audit a TP73 UniProt Projection Against Ensembl and Derived Coding Sequence (CLI Tutorial)](./tp73_uniprot_projection_audit_cli.md) - guided CLI; status `manual/hybrid`. Use this when: CLI-first walkthrough for the integrated TP73 UniProt/Ensembl projection audit, the reusable low-level primitives, parity comparison, and the unsent maintainer-email draft.

### RNA Reads, Splicing & Expression

- `07.01` [Batch-Compare cDNA FASTA.GZ Samples for One Target Gene (CLI Tutorial)](./rna_read_batch_gene_support_cli.md) - guided CLI; status `manual/hybrid`. Use this when: Hands-on batch workflow for interpreting many gzipped cDNA FASTA files, aligning saved reports, and exporting one target-gene cohort sample sheet.
- `07.02` [Map TP53 locus reads with multi-gene sparse indexing (online)](./generated/chapters/12_tp53_multi_gene_sparse_mapping_online.md) - executable reference; status `generated+checked/human-pending`. Use this when: Machine-checked generated chapter; listed in Tutorials for explicit human functional confirmation. Requires explicit online setup for full execution.

### Regulatory, TFBS & Reporter Design

- `08.01` [TFBS Similarity Ranking Tutorial](./tfbs_similarity_ranking_gui.md) - guided GUI; status `manual/hybrid`. Use this when: GUI sign-off tutorial for ranking candidate transcription-factor motifs against one anchor motif over the same DNA span, with JSON export and workflow/ClawBio replay links.
- `08.02` [Promoter Design for VKORC1 / rs9923231](./vkorc1_variant_followup_expert_gui.md) - guided GUI; status `manual/hybrid`. Use this when: Shortest GUI path for the current ClawBio-facing story: start from the dbSNP marker, drive the dedicated Promoter design window, and export one portable promoter-reporter handoff bundle.
- `08.03` [Promoter Design Artifact Slice (Offline Synthetic TP73 Locus)](./generated/chapters/24_promoter_design_artifact_slice_offline.md) - executable reference; status `generated+checked/human-pending`. Use this when: Machine-checked generated chapter; listed in Tutorials for explicit human functional confirmation.
- `08.04` [VKORC1 / rs9923231 PGx Alert -> Mammalian Luciferase Reporter (GUI Tutorial with Matching CLI Commands)](./vkorc1_warfarin_promoter_luciferase_gui.md) - guided GUI; status `manual/hybrid`. Use this when: Hand-written GUI-first tutorial for turning the VKORC1/rs9923231 pharmacogenomic alert into one mammalian promoter-reporter planning workflow with explicit engine and CLI mapping.
- `08.05` [Plan a Reporter Construct Handoff from a Saved Candidate Set](./reporter_construct_handoff_cli.md) - guided CLI; status `manual/hybrid`. Use this when: CLI/agent walkthrough for turning a saved promoter-reporter candidate set into a read-only luciferase macro handoff plan with typed readiness and explicit follow-up commands.

### External Services & Handoffs

- `09.01` [Prepare a Metabion Handoff from Shared External-Service Contracts](./metabion_external_service_handoff_gui_cli.md) - guided GUI/CLI; status `manual/hybrid`. Use this when: Safe handoff rehearsal for Metabion oligo and m-block requests through the shared provider catalog, preflight, and quote bundle contracts; no vendor submission.
- `09.02` [Prepare a GeneArt Handoff from Shared External-Service Contracts](./geneart_external_service_handoff_gui_cli.md) - guided GUI/CLI; status `manual/hybrid`. Use this when: Safe GeneArt cloned-gene and protein-expression quote-handoff rehearsal through the shared provider catalog, preflight, quote bundle, and return_spec contracts; no vendor submission.

### Sequencing Confirmation & QC

- `10.01` [Confirm a Construct from an Imported Sequencing Trace (CLI Tutorial)](./sequencing_confirmation_trace_cli.md) - guided CLI; status `manual/hybrid`. Use this when: Hands-on local walkthrough for seq-trace import, trace-aware seq-confirm run, and shared report export.
- `10.02` [Inspect an Imported Sequencing Trace and Confirm a Construct (GUI Tutorial)](./sequencing_confirmation_gui.md) - guided GUI; status `manual/hybrid`. Use this when: GUI-first walkthrough for raw trace import, baseline-aware confirmation, chromatogram inspection, and lineage reopen.

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

1. Open `File -> Open Tutorial Project...` -> `Primers, PCR & qPCR` -> `04.01 Simple PCR From a Selected Core Region`
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
