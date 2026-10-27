# Isoform Architecture and Genomic Regions

<!-- manual-counterparts:start -->
Related GUI reference: [Saved Genomic Regions](../gui/saved-genomic-regions.md) | [Isoform Architecture Panels (Engine Ops)](../gui/isoform-architecture-panels-engine-ops.md).
These topic links do not imply identical feature coverage across interfaces.
<!-- manual-counterparts:end -->

Isoform architecture panel workflow:

- canonical panel resource:
  - `assets/panels/tp53_isoforms_v1.json`
- operation route (JSON):
  - `ImportIsoformPanel` + `RenderIsoformArchitectureSvg`
- shared shell route:
  - `panels import-isoform ...`
  - `panels inspect-isoform ...`
  - `panels render-isoform-svg ... [--expression-tsv PATH]`
  - `panels validate-isoform ...`
- optional expression overlay:
  - TSV header: `isoform_id	sample_label	value`
  - rows with unknown `isoform_id` are ignored with SVG warnings
  - duplicate cells or non-finite/negative values fail deterministically
- direct expert route:
  - `inspect-feature-expert SEQ_ID isoform PANEL_ID`
  - `render-feature-expert-svg SEQ_ID isoform PANEL_ID OUTPUT.svg`
  - gene-level isoform evidence ledger:
    - `inspect-feature-expert SEQ_ID isoform-evidence PANEL_ID [--annotation-release LABEL] [--rna-read-report-id ID]... [--probe-evidence PATH]... [--cdna-est-resource PATH]... [--expression-tsv PATH] [--occupancy-track NAME]... [--qpcr-report-id ID]...`
    - `render-feature-expert-svg SEQ_ID isoform-evidence PANEL_ID [same evidence options] OUTPUT.svg`
    - the panel and sequence annotation define transcript geometry; optional
      report ids and files add independent RNA-read, cDNA/EST, array,
      expression, projected occupancy, and existing qPCR-candidate layers
    - the JSON/SVG distinguish transcript 5'->3' ordering from genomic-ascending
      ordering and retain transcript-oriented splice donor/acceptor positions
    - probe overlap is reported as an assay constraint, not validation; missing
      evidence remains unknown/not evaluated
    - repeat `--occupancy-track NAME` to align selected projected BED/BigWig
      lanes below the transcript models; `--occupancy-track '*'` includes all
      overlapping projected tracks. Lanes share one absolute score scale and
      remain locus-level evidence rather than isoform-specific regulation
    - see `docs/gene_isoform_occupancy_figure_runbook.md` for a local
      TAp73alpha/DNp73beta CUT&RUN example
  - publication-oriented gene-locus evidence composition:
    - `gene-locus prepare REQUEST_JSON_OR_@FILE`
    - `inspect-feature-expert SEQ_ID gene-locus-evidence PANEL_ID [isoform-evidence options] [--probe-effect-table PATH]... [--probe-effect-contrast TOKEN]... [--probe-effect-coordinate-system ID] [--upstream-bp N] [--downstream-bp N] [--occupancy-layout JSON_OR_@FILE | --occupancy-track NAME ...] [--motif TOKEN]... [--score-kind KIND] [--motif-threshold N] [--motif-top-hits N] [--regulatory-score-tracks JSON_OR_@FILE] [--scale-bar hidden|auto|fixed] [--scale-bar-bp N] [--include-local-source-paths] [--allow-negative]`
    - `render-feature-expert-svg SEQ_ID gene-locus-evidence PANEL_ID [same options] OUTPUT.svg`
    - this pure-read route aligns transcript/CDS metrics, annotation-backed
      start/stop glyphs, grouped projected occupancy, continuous motif scores,
      optional coordinate-aligned PSR/JUC effect rows, and one marker per
      qPCR-targeted junction on a strand-aware 5'->3' axis
    - `gene-locus prepare` is the mutating convenience route. Its typed request
      composes an explicitly network-approved or offline Ensembl gene entry,
      automatic versioned isoform panel, declared BED/BigWig tracks, the same
      report/renderer, and optional PNG/PDF conversion. It returns a
      hash-bound `gentle.gene_locus_evidence_preparation_receipt.v1`; the
      committed offline example is
      `test_files/fixtures/gene_locus_evidence/general_locus_demo/preparation_request.json`
    - the preparation request may include `reporter_architecture_request` and
      `reporter_report_path`. GENtle binds the nested comparison to the newly
      prepared sequence and normalized locus request, then renders its
      canonical architecture rows in the same SVG/PNG/PDF as occupancy,
      chromatin, TF-score, transcript, and scale-bar evidence. Nested callers
      must not supply a competing `locus_evidence_request`
    - an optional `ensembl_regulation` member performs a bounded, source-bound
      overlap query against an explicitly installed Ensembl Regulation index.
      It adds only Ensembl-annotated regulatory-region rows to the canonical
      report and combined figure; generic materialized `regulatory_region`
      features, including SCREEN rows, are not enumerated. The default
      `retain_unavailable` policy keeps missing evidence visible without
      affecting unrelated locus preparation, while `required` fails closed
    - combined reporter rows include a material legend for genomic source,
      spliced exon/cDNA, and schematic `LUC` segments. The offline LOCUSDEMO
      request demonstrates independent TP73, PATZ1, E2F1, and SP1 JASPAR rows
    - repeat `--probe-effect-table` for compact TSV sources whose abundance
      columns are named `log2_mean_*` and whose differential columns are named
      `log2_*_minus_*`; repeat `--probe-effect-contrast` to select a subset.
      `--probe-effect-coordinate-system` is required whenever a table is
      supplied and must match the sequence genome anchor. PSR rows are rendered
      as intervals, JUC rows as junction arcs, and abundance and differential
      values use separate sequential and zero-centered scales. GENtle labels
      these as prioritization evidence, not primer-binding proof,
      differential-expression significance, or isoform validation
    - use a `gentle.gene_locus_occupancy_layout.v1` JSON file to keep cell
      lines/conditions in separate declared groups and to choose
      `shared_group`, `independent`, `fixed`, or explicitly justified
      `shared_all` scaling; GENtle does not infer these meanings from filenames.
      Each requested lane remains visible in JSON/SVG as `available`,
      `not_prepared`, `no_compatible_interval`, or `assembly_mismatch`, and may
      carry stable source id/digest/assembly plus explicit assay, mark, and
      factor metadata
    - `--regulatory-score-tracks` accepts independent JASPAR or normalized
      offline external-model tracks with their own factor/source identity,
      score semantics, strand policy, thresholds, top sites, and scales.
      Multi-source requests can bind factors per exact matrix/model through
      `source_factor_bindings`. Different matrices, models, or providers remain
      non-comparable unless `cross_source_calibrated` carries the same explicit
      calibration id and SHA-256 digest on every member plus a shared-scale
      justification; calibration prose alone never authorizes comparison
    - `--scale-bar` defaults to backward-compatible `hidden`; `auto` selects a
      deterministic 1/2/5 x 10^n length and `fixed` requires a positive
      `--scale-bar-bp`. Local source paths are omitted from portable reports by
      default; `--include-local-source-paths` is an explicit opt-in
    - repeat `--homology-report PATH` to project digest-valid exact-support
      blocks from `gentle.genomic_region_homology_screen.v1` reports onto the
      same locus coordinate axis
    - see `docs/gene_isoform_occupancy_figure_runbook.md` and
      `docs/examples/gene_locus_evidence/patz1_cutrun_layout.json`
    - the fully offline structural/visual proof is:
      `cargo run --bin gentle_cli -- workflow @docs/examples/workflows/patz1_gene_locus_evidence_offline.json`;
      it writes `patz1_gene_locus_evidence.svg` from committed synthetic
      negative-strand PATZ1, PSR/JUC, and occupancy fixtures
    - the general preparation fixture writes and verifies SVG, PNG, and PDF;
      the pinned public Ensembl-116 SERPINE1 acceptance test additionally
      verifies that its reporter architectures and clearly synthetic occupancy
      context are rendered in one canonical figure
  - portable, assembly-bound genomic region sets use one JSON-request command
    family across direct CLI, shared shell, MCP, and adapters:
    - `regions create REQUEST_JSON_OR_@FILE`
    - `regions capture REQUEST_JSON_OR_@FILE`
    - `regions preview-feature REQUEST_JSON_OR_@FILE`
    - `regions materialize-feature REQUEST_JSON_OR_@FILE`
    - `regions list [REQUEST_JSON_OR_@FILE]`
    - `regions inspect REQUEST_JSON_OR_@FILE`
    - `regions update REQUEST_JSON_OR_@FILE`
    - `regions derive REQUEST_JSON_OR_@FILE`
    - `regions import REQUEST_JSON_OR_@FILE`
    - `regions export REQUEST_JSON_OR_@FILE`
    - `regions homology-screen REQUEST_JSON_OR_@FILE`
    - `regions render-homology-svg REPORT_JSON_OR_@FILE OUTPUT.svg`
    - `promoters assess-conserved-modules REQUEST_JSON_OR_@FILE`
    - `region-sets` is accepted as an alias. Create/capture/update/derive/import
      mutate project metadata; list/inspect/export are read-only apart from
      explicit export files
    - `preview-feature` verifies a saved, stranded cofactor hit against exact
      prepared local reference DNA without mutation. Request fields: `set_id`,
      `region_id`, `seq_id`, `expected_region_content_sha256`; optional
      `catalog_path`/`cache_dir` override anchor provenance. Copy the returned
      `feature_materialization.approval_sha256` into `expected_approval_sha256`
      before explicitly calling `materialize-feature`. Apply is undoable and
      rejects changed evidence, reference, DNA or annotation state and duplicate
      attachment. No reference fallback/download or occupancy inference occurs.
      Unstranded evidence stays in the region ledger. See the
      [offline feature-attachment replay](../tutorial/08-14_promoter_cofactor_browser.md).
    - example manual capture:
      `regions capture '{"set_id":"candidates","source":{"source_kind":"sequence_selection","seq_id":"locus","local_start_0based":20,"local_end_0based_exclusive":80,"reference_override":{"species_scientific_name":"Homo sapiens","taxon_id":9606,"assembly_name":"GRCh38","assembly_accession":"GCA_000001405.15","contig_name":"7"}},"purpose":"reporter_candidate"}'`
    - lossless export:
      `regions export '{"set_id":"candidates","json_path":"candidates.json","bed_path":"candidates.bed"}'`.
      The BED writes a sibling manifest by default. Manifest-bound import checks
      both the BED SHA-256 and regenerated canonical rows. Bare BED instead
      requires `bare_bed_reference` plus `set_id_override`; no assembly is
      guessed and no liftover is attempted
    - run the complete offline example with
      `cargo run --bin gentle_cli -- workflow @docs/examples/workflows/portable_genomic_regions_offline.json`
    - homology screening is read-only and uses only already validated local
      genomic-DNA BLAST indexes. Empty `targets[]` resolves all currently
      available indexes into the effective request; explicit optional targets
      remain `unavailable`, while an unavailable required target fails before
      BLAST. Transcriptome/cDNA indexes are rejected, and GENtle never downloads
      or prepares a target implicitly
    - `regions homology-screen` applies the HSP processing budget after the
      complete BLAST search rather than using `-max_target_seqs` as a
      completeness filter. Expected-ortholog labels require explicit expected
      loci with evidence IDs; similarity rank alone remains
      `cross_species_unassigned`
    - query-referenced alignment rows always have the exact query length.
      Target insertions are omitted from display columns but retained under
      `omitted_insertions[]`; substitutions and target deletions remain visible
    - `promoters assess-conserved-modules` consumes a content-valid homology
      report plus independently sourced query-local evidence spans. Its result
      is a reporter-fragment hypothesis and does not claim autonomous promoter
      or enhancer function
    - the synthetic local-BLAST walkthrough is
      `docs/examples/workflows/region_homology_promoter_modules_offline.json`;
      it requires installed BLAST+ tools but no network access
      - complete it with
        `python3 docs/examples/run_region_homology_tutorial.py --gentle target/debug/gentle_cli --output /tmp/gentle-conservation-tutorial`
      - the empty output directory receives the project, exact request/workflow
        JSONs, SVG, four GENtle module assessments, command outputs, and a
        hash-bound acceptance receipt. The companion only forwards requests
        and checks facts; it never computes or rewrites scientific verdicts
    - `scripts/prepare_regulatory_region_indexes.py` prepares a manifest-declared
      multi-gene promoter/upstream comparison bundle through GENtle's existing
      promoter/region extraction, canonical region-set import, and `SaveFile`
      paths. It accepts transcript/TSS windows and typed region sets while
      preserving comparison classes such as `ensembl_regulation` and
      `self_defined`. It collapses exact-identical windows into explicit
      sequence-equivalence classes before building a candidate-pool BLAST
      database and lower-stringency canonical-k-mer signatures.
      `scripts/prepare_ensembl_promoter_region_set.py` first captures exact
      promoter rows from content-verified, non-truncated locus reports through
      GENtle. The candidate database is not presented as a whole-genome
      uniqueness screen; see `docs/regulatory_region_comparison_preparation.md`
  - same command family for restriction-site details:
    - `inspect-feature-expert SEQ_ID restriction CUT_POS_1BASED [--enzyme NAME] [--start START_1BASED] [--end END_1BASED]`
    - JSON output includes `tooltip_lines[]` with the same concise
      caret-marked cut-site summary used by GUI hover tooltips
  - same command family for splicing:
    - `inspect-feature-expert SEQ_ID splicing FEATURE_ID`
    - `render-feature-expert-svg SEQ_ID splicing FEATURE_ID OUTPUT.svg`
  - same command family for transcript-first protein comparison, with or
    without stored external protein evidence:
    - `inspect-feature-expert SEQ_ID protein-comparison [--transcript TRANSCRIPT_ID] [--ensembl-entry ENTRY_ID] [--feature-key KEY]... [--feature-key-not KEY]...`
    - `render-feature-expert-svg SEQ_ID protein-comparison [--transcript TRANSCRIPT_ID] [--ensembl-entry ENTRY_ID] [--feature-key KEY]... [--feature-key-not KEY]... OUTPUT.svg`
    - this route opens the same transcript-first Protein Expert payload used by
      the GUI `Open Derived Protein Expert` action
    - `--ensembl-entry` layers one stored Ensembl protein entry onto that same
      transcript-first compare surface as optional external evidence without
      redefining which transcript-derived products exist
    - `--feature-key` / `--feature-key-not` can trim imported Ensembl or other
      external-opinion feature classes when the comparison view gets crowded
  - same command family for persisted UniProt protein mappings:
    - `inspect-feature-expert SEQ_ID uniprot-projection PROJECTION_ID [--feature-key KEY]... [--feature-key-not KEY]...`
    - `render-feature-expert-svg SEQ_ID uniprot-projection PROJECTION_ID [--feature-key KEY]... [--feature-key-not KEY]... OUTPUT.svg`
    - `CONFLICT` is hidden by default; explicit include/exclude keys let you
      trim noisy UniProt classes such as `VARIANT`, `STRAND`, or `HELIX`
  - direct Ensembl protein metadata routes:
    - `ensembl-protein fetch QUERY [--entry-id ID]`
    - `ensembl-protein list`
    - `ensembl-protein show ENTRY_ID`
    - `ensembl-protein import-sequence ENTRY_ID [--output-id ID]`
    - fetch/import keep Ensembl as optional external evidence plus an ordinary
      first-class protein import path; they do not replace transcript-native
      derivation as the authoritative product model
  - direct Ensembl gene metadata routes:
    - `ensembl-gene fetch QUERY [--species NAME] [--assembly NAME] [--flank-bp N|--flank-5p-bp N --flank-3p-bp N] [--entry-id ID]`
      - `QUERY` is an HGNC-approved human gene symbol such as `FUS` or `TP53`
        when `--species homo_sapiens` is used, or a stable Ensembl gene id such
        as `ENSG00000089280`; `--entry-id` names the stored GENtle metadata
        entry. `--assembly` verifies the assembly returned by Ensembl;
        `--flank-bp` requests the same strand-aware genomic context at both
        ends, while the directional options allow asymmetric context.
      - An HGNC ID such as `HGNC:11998` identifies an HGNC nomenclature record,
        not a GenBank, UniProt, or Ensembl accession. Resolve it to an approved
        symbol or linked database id before choosing a database-specific fetch.
    - `ensembl-gene list`
    - `ensembl-gene show ENTRY_ID`
    - `ensembl-gene import-sequence ENTRY_ID [--output-id ID]`
    - these are useful for one-off live Ensembl gene retrieval without whole
      reference preparation; prepared `genomes extract-gene|region|promoter`
      remains the preferred route for reproducible locus-context work
    - `scripts/fetch_ensembl_cdna_fixtures.sh` wraps this route for RNA-read
      preflight controls: it fetches the Ensembl gene entry, imports the locus,
      derives all transcript cDNAs with `transcripts derive`, exports them to
      `ensembl_<species>_<gene-token>_all.fasta`, and records sidecar provenance
      suitable for `--must-pass-transcript-fasta` /
      `--control-transcript-fasta`
    - `scripts/pancreas_gene_rna_screen.sh run` can call the same fixture
      helper when `--auto-fetch-fixtures` is provided; `--must-pass-gene`,
      `--positive-gene`, and `--control-gene` then resolve to Ensembl cDNA

      FASTA fixtures automatically before the RNA-read preflight
    - `scripts/list_ensembl_paralogs.py` queries Ensembl homology records for
      same-species paralogues and prints one target Ensembl gene id per line by
      default; `scripts/pancreas_gene_rna_screen.sh run --control-human-paralogs`
      uses that helper before BLAST/near-neighbor exploration, fetches the
      resulting cDNA fixtures, and adds available paralog fixtures as explicit
      preflight controls
  - direct Ensembl region/ROI route:
    - `ensembl-region fetch SPECIES CHR START END [--strand +|-] [--output-id ID] [--coord-system-version VERSION]`
      - `SPECIES` is an Ensembl species name such as `homo_sapiens`; `CHR`,
        `START`, and `END` are assembly coordinates for that Ensembl species
        and coordinate-system version; `--output-id` names the local GENtle
        sequence created from the fetched interval.
    - compact coordinates are accepted as
      `ensembl-region fetch SPECIES CHR:START..END[:STRAND]`
    - this imports the live REST sequence slice immediately as an anchored DNA
      sequence with provenance, which lets ClawBio answer "show me this
      interval" without preparing a whole local reference first
  - direct coding-DNA query for one projected UniProt feature:
    - `uniprot feature-coding-dna PROJECTION_ID FEATURE_QUERY [--transcript ID] [--mode genomic_as_encoded|translation_speed_optimized|both] [--speed-profile human|mouse|yeast|ecoli]`
    - returns one structured report per matching transcript feature span, including:
      - amino-acid interval
      - genomic coding DNA as encoded by the genome (spliced, coding-strand orientation)
      - optional preferred-codon translation-speed DNA
      - exon attribution and ordered exon pairs when the feature crosses a splice junction
  - reusable UniProt/Ensembl audit primitives now also exist for external AI orchestration:
    - `uniprot resolve-ensembl-links ...`
    - `uniprot transcript-accounting ...`
    - `uniprot compare-ensembl-exons ...`
    - `uniprot compare-ensembl-peptide ...`
  - the integrated audit persists a human-facing report plus a local unsent maintainer-email draft:
    - `uniprot audit-projection ...`
    - `uniprot audit-parity ...`
    - `uniprot audit-list|show|export ...`
    - `uniprot audit-parity-list|show|export ...`

Introspection models imported isoform panels as sequence-bound project facts:

- `isoform_panel.exists(PANEL_ID)` proves that a panel id is present.
- `isoform_panel.seq_id(PANEL_ID) == SEQ_ID` proves that the panel was imported
  for the sequence currently being inspected or rendered.
- `panels inspect-isoform`, `panels render-isoform-svg`, and
  `RenderIsoformArchitectureSvg` require both facts; `panels validate-isoform`
  validates an external file without mutating project state.

<!-- manual-page-nav:start -->

---

[CLI Manual](../cli.md) | [Previous: Cryptic Splicing](splicing.md) | [Next: Protein and Structure Evidence](protein-evidence.md)
<!-- manual-page-nav:end -->
