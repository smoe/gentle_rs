# Collections

<!-- manual-counterparts:start -->
Related GUI reference: [Project Overview And Collections](../gui/project-overview-and-collections.md) | [Gene Set Inspector And Collection Operations](../gui/gene-set-inspector-and-collection-operations.md).
These topic links do not imply identical feature coverage across interfaces.
<!-- manual-counterparts:end -->

- Shared shell:
  - Commands:
    - `collections run tfbs-scan --tss-collection COLLECTION_ID --motif TOKEN [--path OUTPUT.json]` (validates all members before scanning; exclusive with gene-set/sequence/member-binding arguments)
    - `collections run restriction-scan --tss-collection COLLECTION_ID --enzyme NAME` (same validated membership; no member-sequence override)
    - `collections run digest --tss-collection COLLECTION_ID --enzyme NAME --dry-run` (application still needs `--apply --expected-plan-fingerprint-sha256 SHA256`)
    - `collections run primer-specificity --tss-collection COLLECTION_ID --pair-index N --target-genome GENOME_ID [--member-report MEMBER_ID=PRIMER_REPORT_ID]...` (validates membership before specificity; existing primer-report and backend requirements still apply)

- Shared shell:
  - Commands:
    - `collections run primer-specificity GENE_SET_REPORT_ID --member-report MEMBER_ID=PRIMER_REPORT_ID ... --pair-rank N --target-genome GENOME_ID [--policy JSON_OR_@FILE] [--catalog PATH] [--cache-dir DIR] [--path OUTPUT.json]`
    - `collections run primer-specificity --seq-ids SEQ_ID,... (--pair-rank N | --pair-index N) --target-genome GENOME_ID [--member-report MEMBER_ID=PRIMER_REPORT_ID]... [--policy JSON_OR_@FILE] [--catalog PATH] [--cache-dir DIR] [--path OUTPUT.json]`
    - `collections run restriction-scan GENE_SET_REPORT_ID --member-sequence MEMBER_ID=SEQ_ID ... [--enzyme NAME]... [--max-sites-per-enzyme N] [--no-cut-geometry] [--path OUTPUT.json]`
    - `collections run restriction-scan --seq-ids SEQ_ID,... [--enzyme NAME]... [--max-sites-per-enzyme N] [--no-cut-geometry] [--path OUTPUT.json]`
    - `collections run construct-reasoning-inspection GENE_SET_REPORT_ID --member-sequence MEMBER_ID=SEQ_ID ... [--path OUTPUT.json]`
    - `collections run construct-reasoning-inspection --seq-ids SEQ_ID,... [--path OUTPUT.json]`
    - `collections run tfbs-scan GENE_SET_REPORT_ID --member-sequence MEMBER_ID=SEQ_ID ... --motif TOKEN [--motif TOKEN ...] [--motifs CSV] [--min-llr-bits VALUE] [--min-llr-quantile VALUE] [--per-tf-min-llr-bits TF=VALUE] [--per-tf-min-llr-quantile TF=VALUE] [--max-hits N] [--path OUTPUT.json]`
    - `collections run tfbs-scan --seq-ids SEQ_ID,... --motif TOKEN [same threshold/cap/output options]`
    - `digest SEQ_ID --enzyme NAME [--enzyme NAME ...] [--enzymes CSV] [--output-prefix PREFIX]`
    - `collections run digest GENE_SET_REPORT_ID --member-sequence MEMBER_ID=SEQ_ID ... --enzyme NAME [--enzyme NAME ...] [--enzymes CSV] [--output-prefix PREFIX] [--preview] [--path REPORT.json]`
    - `collections run digest --seq-ids SEQ_ID,... --enzyme NAME [same enzyme/prefix/report options]`
    - `collections run digest SUBJECT [same inputs] --apply --expected-plan-fingerprint-sha256 SHA`
    - `collections run export-pool CONTAINER_ID --path OUTPUT.pool.gentle.json [--pool-id POOL_ID] [--human-id HUMAN_ID]`

<!-- manual-page-nav:start -->

---

[CLI Manual](../cli.md) | [Previous: Project and Sequence State](projects.md) | [Next: Operation History](history.md)
<!-- manual-page-nav:end -->
