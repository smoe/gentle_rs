# TFBS and JASPAR Scans

<!-- manual-counterparts:start -->
Related GUI reference: [TFBS Annotation (Engine Ops)](../gui/tfbs-annotation-engine-ops.md).
These topic links do not imply identical feature coverage across interfaces.
<!-- manual-counterparts:end -->

- Shared shell:
  - Commands:
    - TFBS/JASPAR direct scan helper notes (`features tfbs-scan`):
      - non-mutating structured result schema:
        `gentle.tfbs_hit_scan.v1`
      - accepts either stored `SEQ_ID` or inline ASCII DNA via
        `--sequence-text`
      - uses the same shared local motif registry/IUPAC scoring helpers as the
        current `AnnotateTfbs` and TFBS score-track routes, so CLI, shell, and
        future GUI selection actions can stay in parity
      - `--min-llr-bits` and `--min-llr-quantile` set the default scan gate;
        `--per-tf-min-llr-bits TF=VALUE` and
        `--per-tf-min-llr-quantile TF=VALUE` override those defaults for one
        motif token / TF id / TF name
      - `--max-hits` caps the total retained row count across the scan and the
        report marks that truncation explicitly
      - `--motif ALL` or `--motif *` expands to the whole local motif registry
    - Optional precomputed genome-scan notes
      (`features genomic-motif-evidence`):
      - The rules below describe the original full-scan provider. The separate
        `genome_regulatory_tfbs_subset` provider uses delivered `file_inventory`,
        not original-atlas paths in `scan_file_inventory`. It preserves negative
        source-retained scores by default, exact annotation scope/releases,
        known-empty intersections and separate TSS/transcript ownership. Its
        catalog/gene/region commands work from an empty project; no TP73 anchor,
        `complete.json`, implicit download or DNA annotation is required.
        [Walkthrough and limits](../regulatory_motif_subset.md).
      - returns `gentle.genomic_motif_evidence.v1` without mutating project
        state
      - accepts one genome-anchored `SEQ_ID` (optionally narrowed by local
        0-based half-open `--range`), one saved genomic `--region-set ID`, or
        one or more explicit BED-style genomic `--region` values; the target
        forms are mutually exclusive
      - explicit intervals may be labelled, for example
        `--region TP73=1:3652515..3653715 --region PATZ1=22:31320000..31322000`
      - saved region sets retain their assembly reference in every report row.
        Before any hit query, a saved assembly accession must exactly match the
        package accession; without one, the saved assembly name must match the
        package name or accession. Missing identity or a mismatch returns
        `incompatible_package` with no hits. `--genome-id` cannot override this
        check, and mixed-assembly saved sets are rejected by the region store
      - requires 1-64 explicit JASPAR matrix IDs; `ALL` and `*` are refused
        because a sparse all-motif viewport query would open thousands of
        Parquet partitions
      - discovers the package through `--package` or
        `GENTLE_JASPAR_GENOME_SCAN_PACKAGE`, and DuckDB through `--duckdb`,
        `GENTLE_DUCKDB_BIN`, or `duckdb` on `PATH`
      - neither dependency is probed at startup. A missing package or DuckDB
        executable returns a typed unavailable report and leaves local
        `features tfbs-scan` and score-track commands fully usable
      - queries the package catalog read-only, selects exact payload paths from
        `scan_file_inventory`, and applies hard file, row, output-size, and
        timeout bounds; `--path` writes the same machine-readable report
      - coordinates are BED 0-based half-open. Minus-strand package hits stay
        in forward-reference coordinates; an anchored reverse sequence gets
        correctly projected local coordinates and orientation
      - per-motif coverage rows state the package retention floor, density
        limiting, requested threshold completeness, missing motifs, and row
        truncation. A storage floor is not a biological binding threshold
      - package-native scores use the package's recorded pseudocount and score
        policy and are not numerically interchangeable with GENtle's local
        score tracks; the two sources remain separate evidence lanes
      - the package's Ensembl release describes source FASTA provenance, not a
        gene-annotation model. Saved sets additionally verify declared assembly
        identity, but neither that check nor contig aliases, containment, and
        anchor geometry establish matching per-contig sequence hashes
      - examples:
        - `gentle_cli shell 'features genomic-motif-evidence tp73_context --motif MA0861.2 --range 15564..16764 --path /tmp/tp73_precomputed_tfbs.json'`
        - `gentle_cli shell 'features genomic-motif-evidence --motif MA0861.2 --motif MA0525.2 --region TP73=1:3652515..3653715 --region PATZ1=22:31320000..31322000 --path /tmp/tp73_patz1_precomputed_tfbs.json'`
        - `gentle_cli shell 'features genomic-motif-evidence --region-set conserved_promoters --motif MA0861.2 --path /tmp/conserved_promoter_tfbs.json'`

<!-- manual-page-nav:start -->

---

[CLI Manual](../cli.md) | [Previous: Accession-Pinned TSS Profiles](tss-profiles.md) | [Next: Repeats and Regulatory Annotations](repeats-regulation.md)
<!-- manual-page-nav:end -->
