# Gene-set publication reports

GENtle can assemble a responsive analysis site and a printable companion from
one portable JSON manifest. This is the reusable replacement for analysis-local
HTML builders and separately maintained supplement scripts.

## Contract

`gentle.gene_set_publication_request.v1` identifies the title, conditions,
metrics, interpretation and provenance; a delimited primer table with explicit
column mapping; an ordered gene and figure inventory; supporting downloads;
and which figures are appended to print.

The generator resolves that request into
`gentle.gene_set_publication_report.v1`. Both renderers consume that exact
record. The complete primer list, gene order, captions, and figure inventory
therefore cannot drift between HTML and print. `resolved-report.json` remains
in the output bundle for audit and reuse.

The HTML is the primary interactive presentation. It provides a responsive
overview with metric cards, methods/rationale, a complete primer table, and
per-gene figure galleries. The PDF is a printable companion with the same
narrative and primer table followed by every selected locus, primer/exon-map,
and virtual-gel figure.

## Run

```sh
cargo run --bin gentle_publication_report -- \
  analysis/gene_set_publication_request.json \
  analysis/publication-report \
  --pdf
```

Without `--pdf`, GENtle produces the web/Markdown bundle only. PDF generation
uses `pandoc`, LuaLaTeX, `pdfunite`, and Ghostscript; override executable paths
with `GENTLE_PANDOC_BIN`, `GENTLE_PDFUNITE_BIN`, and
`GENTLE_GHOSTSCRIPT_BIN` when required.

Relative inputs resolve from the manifest directory. Assets are copied under
`figures/` and `data/`, making the result independently publishable.

Minimal shape (real reports normally add metrics, narrative, downloads, and
three figures per gene):

```json
{
  "schema": "gentle.gene_set_publication_request.v1",
  "report_id": "example-panel",
  "title": "Example assay panel",
  "primer_table": { "path": "final_primer_panel.tsv", "delimiter": "\\t" },
  "genes": [{
    "gene_symbol": "GENE1",
    "figures": [{
      "figure_id": "gene1-locus",
      "kind": "locus_evidence",
      "source_path": "GENE1_locus.svg",
      "include_in_pdf": true
    }]
  }]
}
```

## Scientific boundaries

This layer presents persisted evidence; it does not invent or re-score it.
Clariom abundance, Clariom activity differences, formal gene-level
differential-expression statistics, isoform evidence, and primer readiness
must remain separately labelled. Abundance near a primer region supports
assay-region prioritisation but does not prove primer binding. Differential
activity is not formal significance, and neither lane alone proves an isoform.

## Adapting another gene set

Change the manifest's gene records, primer table, figure inventory, conditions,
narrative, and provenance. No renderer edit is required. A PARK7-style overview
is therefore routine report structure rather than a one-off design.
