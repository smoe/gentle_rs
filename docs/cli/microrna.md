# MicroRNA Target Scans

<!-- manual-counterparts:start -->
Related GUI reference: [microRNA Target Scan](../gui/microrna-target-scan.md).
These topic links do not imply identical feature coverage across interfaces.
<!-- manual-counterparts:end -->

microRNA target-site scan:

- `gentle_cli mirna explain-seed hsa-miR-96-5p`
  - resolves the mature sequence from the built-in seed catalog and reports
    canonical target motifs such as `8mer`, `7mer-m8`, `7mer-A1`, and `6mer`.
- `gentle_cli mirna catalog-show hsa-miR-96-5p`
  - shows the catalog-backed identifiers (`MIMAT0000095`, NCBI Gene `407053`)
    and source notes.
- `gentle_cli --state PROJECT.json mirna scan-target hsa-miR-96-5p TP73 --regions 3utr,exon,intron,boundary --seed-classes 8mer,7mer-m8,7mer-A1,6mer --boundary-flank-bp 25 --format json`
  - scans annotated transcript partitions in the loaded TP73 sequence/project,
    returning `gentle.mirna_target_scan.v1` JSON.
  - uses a 25 bp boundary flank by default when splice-boundary regions are
    requested; adjust with `--boundary-flank-bp` for wider or narrower
    junction context.
  - reports every requested exact seed class independently, so stronger sites
    can also appear in shorter seed-class counts.
  - reports exact seed candidates as sequence evidence only; orthologous rat
    Tp73 PMID 37099528 context is tagged separately from direct human
    validation.

<!-- manual-page-nav:start -->

---

[CLI Manual](../cli.md) | [Previous: RNA Secondary Structure](rna-structure.md) | [Next: Measured Gel Images](gel-images.md)
<!-- manual-page-nav:end -->
