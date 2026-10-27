# Variant and Reporter Evidence

<!-- manual-counterparts:start -->
Related GUI reference: [Promoter Cofactors](../gui/promoter-cofactors.md).
These topic links do not imply identical feature coverage across interfaces.
<!-- manual-counterparts:end -->

- Shared shell:
  - Commands:
    - `variant annotate-promoters SEQ_ID [--gene-label LABEL] [--transcript-id ID] [--upstream-bp N] [--downstream-bp N] [--collapse transcript|gene|tss-cluster[:TOLERANCE_BP]]`
    - `variant promoter-context SEQ_ID [--variant ID] [--gene-label LABEL] [--transcript-id ID] [--promoter-upstream-bp N] [--promoter-downstream-bp N] [--tfbs-focus-half-window-bp N] [--path FILE.json]`
    - `variant reporter-fragments SEQ_ID [--variant ID] [--gene-label LABEL] [--transcript-id ID] [--retain-downstream-from-tss-bp N] [--retain-upstream-beyond-variant-bp N] [--max-candidates N] [--path FILE.json]`
    - `reporters list [--catalog PATH] [--filter TEXT] [--limit N] [--output FILE.json]`
    - `reporters recommend [--catalog PATH] [--assay NAME] [--chassis HOST] [--live true|false] [--color COLOR] [--class CLASS] [--excitation-nm NM] [--emission-nm NM] [--fusion MODE] [--max-length-bp N] [--forbid-motif IUPAC] [--substrate-allowed true|false] [--limit N] [--output FILE.json]`
      - recommendation JSON includes `biological_intent` so agent-facing
        synthetic-biology flows can route promoter, luciferase, fusion, live,
        and endpoint reporter requests without guessing from prose
      - local reporter catalogs should provide `sequence_sha256`; legacy
        `sequence_sha1` values are accepted for older catalogs but are not
        recomputed by GENtle
    - `reporters export-corpus OUTPUT.json|OUTPUT.jsonl [--catalog PATH] [--format json|jsonl]`
    - `reporters plan-handoff CANDIDATE_SET.json [--candidate-id ID] [--catalog PATH] [--backbone-seq-id ID] [--backbone-path PATH] [--backbone-catalog-id ID] [--helper-catalog PATH] [--reference-fragment-seq-id ID] [--alternate-fragment-seq-id ID] [--output-prefix PREFIX] [--output FILE.json]`
      - handoff JSON includes `biological_intent =
        allele_paired_promoter_luciferase_reporter_handoff` for the V1
        promoter-luciferase bridge
      - `--backbone-catalog-id` requests an exact helper-vector identity. GENtle
        compares the loaded or path-supplied sequence with catalog-owned
        accession/version, length, topology, and required annotation bounds;
        failed or unavailable validation emits no construct commands
      - the built-in pGL4.10[luc2] entry identifies Promega E6651 as
        `AY738222.1` (4,242 bp, circular) but bundles no commercial-vector
        sequence. Retrieval remains an explicit `FetchGenBankAccession`
        operation. KpnI/Acc65I and SacI/EcoICRI equivalence is reported without
        pretending the absent names are loaded in the active enzyme catalog
    - `variant materialize-allele SEQ_ID --allele reference|alternate [--alternate-base BASE] [--variant ID] [--output-id ID]`
      - `--alternate-base` (alias `--alternate-allele`) selects one reported
        A/C/G/T alternative; required for multiallelic markers and invalid for
        `--allele reference`. Selection uses the loaded sequence's forward
        orientation, not clinical/gene-strand allele nomenclature

<!-- manual-page-nav:start -->

---

[CLI Manual](../cli.md) | [Previous: Construct Reasoning](construct-reasoning.md) | [Next: Synthetic Sequence Design](sequence-design.md)
<!-- manual-page-nav:end -->
