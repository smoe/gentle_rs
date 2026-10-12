# TP53 sparse-origin mapping screenshots

These native GENtle captures were recorded on 2026-10-01 from the canonical
online workflow after reusing a fully prepared Human GRCh38 Ensembl 116 cache.
The workflow processed the committed 10,000-read cDNA fixture and stored report
`tp53_family_sparse_template` against TP53 mRNA feature `1`.

The run is a boundary demonstration, not a three-locus family comparison. The
extracted sequence contains TP53 annotation only. GENtle therefore matched
TP53, reported TP63 and TP73 as absent from local annotation, added no extra
transcript lanes, seed-passed 11/10,000 reads, and deferred alignment.

| File | SHA-256 |
| --- | --- |
| `01-project-graph.png` | `d3f0383dd3dbed334ae12c47a182f7c7ebf24bcc0307c56d194b19730947e21e` |
| `02-sparse-report-warnings.png` | `8454fa514ff6eb00f481ed5d2fae2a92a013dcc689b2cb1080e062a0415c9ed1` |
