# Documentation Figures

`gentle_system_overview.svg` is the hand-authored homepage schematic used near
the top of the README. Unlike the cloning and analysis showcases below, it is
not rendered from an engine operation; it is a maintained explanatory figure
that summarizes how interfaces, imported context, the shared engine, and
provenance interact. The Mermaid block in `README.md` mirrors the same
structure in text-native form.

`gibson_two_fragment_protocol_cartoon.svg` is a deterministic render of the
built-in protocol cartoon `gibson.two_fragment`.

Regenerate it from the repository root with:

```sh
cargo run --quiet --bin gentle_cli -- \
  protocol-cartoon render-svg \
  gibson.two_fragment \
  docs/figures/gibson_two_fragment_protocol_cartoon.svg
```

`gibson_single_insert_protocol_cartoon.svg` is a deterministic render of the
built-in protocol cartoon `gibson.single_insert_dual_junction`. It is the
canonical single-insert Gibson mechanism strip used in the README to show both
destination-insert junctions explicitly.

Regenerate it from the repository root with:

```sh
cargo run --quiet --bin gentle_cli -- \
  protocol-cartoon render-svg \
  gibson.single_insert_dual_junction \
  docs/figures/gibson_single_insert_protocol_cartoon.svg
```

`gibson_single_insert_readme.plan.json` is the deterministic pGEX + insert demo
Gibson plan used to generate the README lineage figure. It is not hand-drawn;
it reuses the same tutorial destination/insert pair and the same `SmaI`
cutpoint (`941..941`) used in the Gibson specialist walkthrough.

`gibson_single_insert_lineage.svg` and `gibson_single_insert_lineage.png` are
the corresponding project-level lineage graph exports after applying that
single-insert Gibson plan. They show one Gibson operation with two inputs and
three concrete outputs (left primer, right primer, assembled product).

Regenerate the lineage assets from the repository root with:

```sh
cargo run --quiet --bin gentle_cli -- \
  --state /tmp/gibson_readme_lineage.state.json \
  workflow @docs/examples/workflows/gibson_specialist_testing_baseline.json

cargo run --quiet --bin gentle_cli -- \
  --state /tmp/gibson_readme_lineage.state.json \
  gibson apply @docs/figures/gibson_single_insert_readme.plan.json

cargo run --quiet --bin gentle_cli -- \
  --state /tmp/gibson_readme_lineage.state.json \
  render-lineage-svg docs/figures/gibson_single_insert_lineage.svg

cargo run --quiet --bin gentle_examples_docs -- \
  svg-png \
  docs/figures/gibson_single_insert_lineage.svg \
  docs/figures/gibson_single_insert_lineage.png \
  --scale 2
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

`tp73_cdna_genomic_dotplot.workflow.json` renders an offline TP73 cDNA-vs-
genomic dotplot from the local `test_files/tp73.ncbi.gb` fixture. The workflow
loads the TP73 genomic locus, derives transcript feature `n-100`
(`feature_ids=[99]`, `NM_001126241.3`), computes a pair-forward dotplot against
the same genomic source, and writes the auditable SVG intermediate
`tp73_cdna_genomic_dotplot.svg`.

`tp73_cdna_genomic_dotplot.png` is the README-friendly raster derivative. It is
rendered from the SVG with `resvg` through `gentle_examples_docs svg-png` while
dropping the export metadata/title text and keeping the basepair axis labels.

Regenerate the TP73 dotplot assets from the repository root with:

```sh
cargo run --quiet --bin gentle_cli -- \
  --state /tmp/tp73_readme_dotplot.state.json \
  workflow @docs/figures/tp73_cdna_genomic_dotplot.workflow.json

cargo run --quiet --bin gentle_examples_docs -- \
  svg-png \
  docs/figures/tp73_cdna_genomic_dotplot.svg \
  docs/figures/tp73_cdna_genomic_dotplot.png \
  --drop-dotplot-metadata
```
