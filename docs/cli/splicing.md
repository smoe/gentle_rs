# Cryptic Splicing

<!-- manual-counterparts:start -->
Related GUI reference: [Splicing Expert Workspace](../gui/splicing-expert-workspace.md) | [Splicing Structure And Evidence Views](../gui/splicing-structure-and-evidence-views.md) | [Cryptic Splicing Screen](../gui/cryptic-splicing-screen.md).
These topic links do not imply identical feature coverage across interfaces.
<!-- manual-counterparts:end -->

Cryptic-splicing inspection:

- `splicing cryptic-screen REQUEST_JSON_OR_@FILE`
- `splicing cryptic-render REQUEST_JSON_OR_@FILE OUTPUT.svg`
- `splicing cryptic-export REQUEST_JSON_OR_@FILE OUTPUT.json`
- `splicing cryptic-overlay REQUEST_JSON_OR_@FILE [--output OUTPUT.json]`
- `splicing cryptic-protein REQUEST_JSON_OR_@FILE [--output OUTPUT.json]`
- `resources sync-maxent INPUT.directory_or_zip [OUTPUT.json] [--source-url URL] [--retrieved-on YYYY-MM-DD] [--redistribution-status STATUS]`
- `resources status` includes the active optional MaxEnt model's validation and
  provenance status.
- The request identifies exactly one loaded-sequence span using 1-based
  inclusive coordinates and a forward/reverse orientation. Optional
  `insert_span` labels vector-versus-insert boundaries, while optional
  `cds_feature_id` asks GENtle to project sequence-removal consequences against
  one CDS feature.
- The v1 implementation deterministically enumerates bounded `GT...AG`
  configurations and applies `min_pseudo_intron_bp`,
  `max_pseudo_intron_bp`, and `max_candidate_pairs`. The report records both
  the total admissible count, evaluated count, reported count, ranking
  completeness, and the exact truncation rule.
- `model_policy` defaults to `structural_only`. With a valid user-supplied
  snapshot, `use_active_maxent` scores donor and acceptor sites separately and
  uses their sum only as a prioritization heuristic. A capped scored run retains
  the deterministic global top-K after evaluating every admissible pair. This
  is not a probability or an observed-junction claim.
- GENtle does not download, bundle, or redistribute MaxEntScan tables. Supply a
  reviewed local directory/ZIP to `resources sync-maxent`; the normalized
  runtime snapshot records source URL, retrieval date, redistribution status,
  per-table SHA-256 values, and one content fingerprint. Missing, corrupt, or
  fingerprint-mismatched models leave the structural report available with
  model evidence absent/failed.
- `cryptic-overlay` accepts
  `{"screen_request": {...}, "rna_read_report_ids": [...], "nearby_tolerance_bp": 2}`.
  Only RNA reports bound to the same source-sequence digest and coordinate
  convention participate. Exact/nearby/no-retained/not-assessable are distinct;
  no retained support must not be read as biological absence. Stored RNA
  junction boundaries use ascending loaded-sequence coordinates; GENtle
  normalizes candidate geometry before matching reverse-strand screens.
- `cryptic-protein` accepts
  `{"screen_request": {...}, "uniprot_projection_id": "..."}` and reports
  overlap with digest-bound projected UniProt feature geometry as a separate
  evidence layer. It does not assert that a candidate is spliced.
- Example request file:

  ```json
  {
    "seq_id": "tap73alpha_cdna",
    "start_1based": 1,
    "end_1based": 1908,
    "strand": "forward",
    "insert_span": {"start_1based": 1, "end_1based": 1908},
    "cds_feature_id": 0,
    "model_policy": "structural_only",
    "min_pseudo_intron_bp": 50,
    "max_pseudo_intron_bp": 1500,
    "max_candidate_pairs": 1000
  }
  ```

  Then run `gentle_cli shell 'splicing cryptic-screen @request.json'`,
  `gentle_cli shell 'splicing cryptic-export @request.json report.json'`, or
  `gentle_cli shell 'splicing cryptic-render @request.json output.svg'`.
  These reports are read-only and computed on demand; there is deliberately no
  list/show store for candidate rows.

<!-- manual-page-nav:start -->

---

[CLI Manual](../cli.md) | [Previous: Source-Coherent Transcripts](transcripts.md) | [Next: Isoform Architecture and Genomic Regions](isoform-architecture.md)
<!-- manual-page-nav:end -->
