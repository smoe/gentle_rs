# Protein and Structure Evidence

<!-- manual-counterparts:start -->
Related GUI reference: [Loading and creating sequences](../gui/loading-and-creating-sequences.md).
These topic links do not imply identical feature coverage across interfaces.
<!-- manual-counterparts:end -->

UniProt mapping capability status:

- shared shell (`gentle_cli shell`, GUI shell): supported via `uniprot ...` commands
  - `uniprot fetch QUERY [--entry-id ID]`
    - `QUERY` is a UniProtKB/Swiss-Prot accession or entry name, for example
      `P04637` or `P53_HUMAN`; `--entry-id` names the stored GENtle metadata
      entry and is not the UniProt accession unless you deliberately reuse it.
  - `uniprot import-swissprot PATH [--entry-id ID]`
  - `uniprot list` / `uniprot show ENTRY_ID`
  - `uniprot map ENTRY_ID SEQ_ID [--projection-id ID] [--transcript ID]`
  - `uniprot projection-list [--seq SEQ_ID]`
  - `uniprot projection-show PROJECTION_ID`
  - `uniprot feature-coding-dna PROJECTION_ID FEATURE_QUERY [--transcript ID] [--mode genomic_as_encoded|translation_speed_optimized|both] [--speed-profile human|mouse|yeast|ecoli]`
  - `uniprot resolve-ensembl-links PROJECTION_ID [--transcript ID]`
  - `uniprot transcript-accounting PROJECTION_ID [--transcript ID]`
  - `uniprot compare-ensembl-exons PROJECTION_ID [--transcript ID] [--ensembl-entry ID]`
  - `uniprot compare-ensembl-peptide PROJECTION_ID [--transcript ID] [--ensembl-entry ID]`
  - `uniprot build-linked-transcript-inventory REQUEST_JSON_OR_@FILE`
    - builds `gentle.uniprot_linked_transcript_inventory.v1` from explicit
      UniProt entry/isoform/transcript links and one or more transcript FASTA
      resources; the request declares assembly, annotation release, resource
      identity, and output path;
    - every resolved record carries exact mature-cDNA SHA-256 and length;
      missing or conflicting transcript sequences remain explicit denominator
      records and are never inferred from a gene name or protein identifier.
  - `uniprot audit-projection PROJECTION_ID [--transcript ID] [--ensembl-entry ID] [--report-id ID]`
  - `uniprot audit-parity PROJECTION_ID [--transcript ID] [--ensembl-entry ID] [--report-id ID]`
  - `uniprot audit-list [--seq SEQ_ID]`
  - `uniprot audit-show REPORT_ID`
  - `uniprot audit-export REPORT_ID OUTPUT.json`
  - `uniprot audit-parity-list [--seq SEQ_ID]`
  - `uniprot audit-parity-show REPORT_ID`
  - `uniprot audit-parity-export REPORT_ID OUTPUT.json`
- shared feature-expert route: supported for stored UniProt genome projections
  via the same expert command family used by splicing/isoform inspection
  - `inspect-feature-expert SEQ_ID uniprot-projection PROJECTION_ID [--feature-key KEY]... [--feature-key-not KEY]...`
  - `render-feature-expert-svg SEQ_ID uniprot-projection PROJECTION_ID [--feature-key KEY]... [--feature-key-not KEY]... OUTPUT.svg`
  - default behavior hides UniProt `CONFLICT` annotations unless they are
    explicitly re-included with `--feature-key CONFLICT`
  - each transcript now gets its own projected protein-coverage rail, so
    missing or partial domains stay transcript-specific instead of every
    isoform inheriting the full reference-protein annotation set
  - when transcript features do not carry `cds_ranges_1based`, GENtle now
    prefers compatible `CDS` features before falling back to exon spans, which
    keeps RefSeq-style records closer to the actual encoded product
  - exported UniProt-projection SVGs now combine:
    - a coordinate-true genomic transcript/exon panel that preserves exon
      positions and introns
    - a shared genomic-exon-column transcript/product panel that reuses one
      exon-family color across aligned transcript segments and translated
      peptide segments
    - isoform-local protein axes, which make skipped-exon domain loss easier to
      read without discarding the genomic panel above
  - those CDS/exon colors are now keyed by genomic exon family/position, not by
    local segment order inside each transcript/product lane, so vertically
    matching locus exons keep the same color even when transcript boundaries
    differ slightly
  - membrane/topology-style features (`SIGNAL`, `TRANSIT`, `TOPO_DOM`,
    `TRANSMEM`, `INTRAMEM`) render in a dedicated lower band beneath the main
    protein rail so dense labels stay readable
- shared shell (`gentle_cli shell`, GUI shell): GenBank accession import
  - `genbank fetch ACCESSION [--as-id ID]`
    - `ACCESSION` is an NCBI GenBank nucleotide accession; `--as-id` names the
      local GENtle sequence id created from the fetched record.
- shared shell (`gentle_cli shell`, GUI shell): dbSNP-guided annotated region extraction
  - `dbsnp fetch RS_ID GENOME_ID [--flank-bp N] [--output-id ID] [--annotation-scope none|core|full] [--max-annotation-features N] [--catalog PATH] [--cache-dir PATH]`
    - `RS_ID` is a dbSNP identifier such as `rs9923231`; `GENOME_ID` is a
      prepared GENtle genome catalog id, and `--output-id` names the local
      sequence extracted around that variant.
    - A custom catalog whose reference identity yields no recognizable assembly
      family produces `dbsnp_assembly_check=unavailable`. The fetched slice remains
      inspectable, but allele materialization is refused even if its base matches;
      a matching base alone is not assembly or coordinate authority.
- engine operations behind those commands:
  - `FetchUniprotSwissProt`, `ImportUniprotSwissProt`, `ProjectUniprotToGenome`,
    `FetchGenBankAccession`, `FetchDbSnpRegion`
- reverse-translation is now also available through dedicated shared-shell /
  direct CLI commands:
  - `reverse-translate run PROTEIN_SEQ_ID [--output-id ID] [--speed-profile human|mouse|yeast|ecoli] [--speed-mark fast|slow] [--translation-table N] [--target-anneal-tm-c F] [--anneal-window-bp N]`
  - `reverse-translate list-reports [PROTEIN_SEQ_ID]`
  - `reverse-translate show-report REPORT_ID`
  - `reverse-translate export-report REPORT_ID OUTPUT.json`
  backed by `ReverseTranslateProteinSequence` plus persisted
  reverse-translation report inspection/export routes.
- construct-reasoning protein-to-DNA handoff is now also available through the
  shared shell:
  - `construct-reasoning build-protein-dna-handoff SEQ_ID PROTEIN_SEQ_ID [--transcript TRANSCRIPT_ID] [--projection-id ID] [--ensembl-entry ID] [--feature-query TEXT] [--ranking-goal balanced_provenance|native_fidelity|expression_optimized] [--speed-profile human|mouse|yeast|ecoli] [--speed-mark fast|slow] [--translation-table N] [--target-anneal-tm-c F] [--anneal-window-bp N] [--objective-id ID] [--graph-id ID]`
  - `construct-reasoning list-graphs [SEQ_ID]`
  - `construct-reasoning show-graph GRAPH_ID`
  - `construct-reasoning list-inspection-actions GRAPH_ID [--fact-id ID] [--annotation-id ID] [--candidate-id ID] [--evidence-id ID] [--seq-id ID] [--action-kind KIND] [--summary-id ID]`
  - `construct-reasoning run-inspection-action GRAPH_ID ACTION_ID [--word-size N] [--step N] [--max-mismatches N] [--tile-bp N] [--id DOTPLOT_ID] [--render-svg OUTPUT.svg]`
  - `construct-reasoning set-annotation-status GRAPH_ID ANNOTATION_ID draft|accepted|rejected|locked`
  - `construct-reasoning write-annotation GRAPH_ID ANNOTATION_ID`
  - `construct-reasoning export-graph GRAPH_ID OUTPUT.json`
  backed by `BuildProteinToDnaHandoffReasoning` plus the existing persisted
  construct-reasoning graph store.
  - `construct-reasoning show-graph` now returns both the full portable graph
    and a compact CLI-facing `summary` block:
    - `summary_lines`
    - `warning_lines`
    - `fact_summaries`
    - graph-level `inspection_actions[]` entries for deterministic
      repeat/similarity dotplot handoffs, including `action_id`,
      `rationale`, `driving_evidence_ids[]`, source fact/annotation/summary
      ids, mode, focus range, and optional RepeatMasker/UCSC `rmsk`-backed
      repeat-family provenance
    - current fact summaries now include adapter-capture review plus the new
      similarity-derived predictor rows for PCR/amplification,
      nanopore/direct-sequencing, repeat-driven mapping, and cloning stability
    - repeat/similarity fact summaries expose `task_severities[]` rows with
      task (`pcr`, `nanopore_sequencing`, `read_mapping`,
      `cloning_stability`, `construct_maintenance`), explicit objective
      applicability/basis, intrinsic `base_severity`/`base_score`, transparent
      `objective_adjustment`, effective severity/score, rationale, and
      supporting evidence ids, plus compact `task_severity: ...` detail lines
    - construct objectives may declare authoritative `intended_tasks[]`;
      missing intent retains conservative legacy positive inference while an
      explicit empty list means none of the current task vocabulary applies
    - repeat/mobile-element facts summarize overlapping materialized
      RepeatMasker/UCSC `rmsk`-style annotations as
      `curated_repeat_support[]` rows with repeat name/class/family and
      supporting evidence ids, while non-overlapping rmsk rows do not upgrade
      unrelated internal predictions
  - `construct-reasoning list-inspection-actions` returns the same graph-level
    `inspection_actions[]` objects, filtered by source fact/annotation/summary
    ids, candidate ids, driving evidence ids, sequence id, or action kind when
    requested, so CLI/agent/ClawBio layers do not rediscover GUI dotplot
    affordances from labels.
    Responses include `snapshot_status` (`current|stale|unknown` plus reasons),
    and action records can carry typed `repeat_family_provenances[]` for all
    overlapping curated families rather than one parsed note string.
  - `construct-reasoning run-inspection-action` resolves one `action_id`,
    feeds its mode and focus range to the shared `ComputeDotplot` operation,
    stores the resulting dotplot payload, and optionally invokes
    `RenderDotplotSvg` for an evidence artifact. The structured response also
    includes the resolved `compute_parameters` block used for the shared
    operation. The stored dotplot and its list summary retain `op_id`, `run_id`,
    and a typed `gentle.dotplot_inspection_provenance_citation.v1` citation to
    the graph, action, rationale, source facts/annotations/candidates/summaries,
    driving evidence, fingerprint, and resolved dotplot parameters. A
    fingerprinted stale graph is rejected until the reasoning build route that
    created it is rerun, preventing old action ranges from being run against
    changed sequence/features, objective intent, or reasoning rules. Request or
    source mismatches persist as `status=fail`; unverifiable legacy snapshots
    remain `unknown`. Manual dotplots omit the citation rather than recording a
    failed one.
  - `construct-reasoning set-annotation-status` updates one persisted
    annotation candidate in place and returns:
    - the updated portable graph
    - the updated `annotation_candidate`
    - the same shared `summary` block exposed by `construct-reasoning show-graph`
  - `construct-reasoning write-annotation` materializes one accepted or locked
    generated annotation candidate as an ordinary sequence feature and returns:
    - the refreshed portable graph
    - the refreshed or already-backed `annotation_candidate`
    - a `writeback` report (`gentle.annotation_candidate_writeback.v1`) with
      `created`, `already_present`, and the resulting `feature_id`
    - the same shared `summary` block exposed by `construct-reasoning show-graph`
  - for adapter/linker restriction-capture facts, that summary now surfaces the
    same shared-engine review as the GUI inspector:
    - per-motif conflicts across capture plus retrieval sites
    - the stronger `all named adapter motifs already occur on the insert` case
    - `methylation may rescue this, but the enzyme-specific methylation
      knowledge base is still incomplete`
- additional protein-sequence operations are still available through the shared
  engine contract (`apply_operation` / workflow JSON):
  - `ImportUniprotEntrySequence`
  - `DeriveProteinSequences`
    - transcript-first and self-sufficient; UniProt is optional comparison
      evidence, not the source of truth for what peptide products exist
  - protein and peptide molecule labels are treated as one first-class
    `protein` family for shared import/export and FASTA metadata purposes

## Protease Catalogs and Digestion

- `proteases list [--filter TEXT] [--output OUTPUT.json]`
  - Lists the built-in protease catalog through one deterministic JSON report.
  - Search matches the curated biotech-facing metadata now carried by the
    protease rows:
    - names and aliases
    - cleavage-pattern notation
    - specificity summaries
    - common applications such as `proteomics` or `tag_removal`
  - Catalog rows also drive the shared protein-digest operation below.
- `proteases show QUERY [--output OUTPUT.json]`
  - Resolves one protease by exact name or alias first, then by unique
    catalog-search match.
  - Useful for sequence-specific tag-cleavage proteases such as `TEV protease`,
    `HRV 3C protease`, `Thrombin`, `Factor Xa`, or residue-specific entries
    such as `Trypsin`.
- `proteases digest SEQ_ID PROTEASE[,PROTEASE...] [--output-prefix PREFIX] [--min-length-aa N] [--predict-only]`
  - Predicts cleavage sites on a first-class protein/peptide sequence using the
    built-in protease catalog and returns
    `gentle.protease_digest_report.v1`.
  - By default the predicted peptide products are materialized as first-class
    `peptide` sequences with source-protein and transcript-derived provenance
    qualifiers when the source protein carries them.
  - `--predict-only` keeps the operation report-only and does not create
    sequence entries; `--min-length-aa` filters short peptide products.
- `proteases digest-gel-svg SEQ_ID PROTEASE[,PROTEASE...] OUTPUT.svg [--min-length-aa N] [--ladder NAME] [--ladders CSV]`
  - Renders the same prediction-only protease digest as a protein-gel SVG using
    the shared protein ladder renderer.
  - The operation still returns `gentle.protease_digest_report.v1`, so ClawBio
    can combine a compact structured digest report with the PNG raster that the
    wrapper derives from the declared SVG artifact.

<!-- manual-page-nav:start -->

---

[CLI Manual](../cli.md) | [Previous: Isoform Architecture and Genomic Regions](isoform-architecture.md) | [Next: Primer Design and Handoff](primers.md)
<!-- manual-page-nav:end -->
