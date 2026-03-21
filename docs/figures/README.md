# Documentation Figures

`gibson_two_fragment_protocol_cartoon.svg` is a deterministic render of the
built-in protocol cartoon `gibson.two_fragment`.

Regenerate it from the repository root with:

```sh
cargo run --quiet --bin gentle_cli -- \
  protocol-cartoon render-svg \
  gibson.two_fragment \
  docs/figures/gibson_two_fragment_protocol_cartoon.svg
```

`tp53_ensembl116_panel_source.gb` is a synthetic TP53 locus slice whose
feature geometry comes from the Ensembl 116 GRCh38 TP53 annotation. Its
sequence bases are placeholder `N`s because the isoform-architecture render
consumes transcript/CDS coordinates plus `cds_ranges_1based`, not genomic base
content.

`tp53_isoform_architecture.svg` is rendered from that source plus the curated
panel resource `assets/panels/tp53_isoforms_v1.json`.

Regenerate it from the repository root with:

```sh
cargo run --quiet --bin gentle_cli -- \
  workflow @docs/figures/tp53_isoform_architecture.workflow.json
```
