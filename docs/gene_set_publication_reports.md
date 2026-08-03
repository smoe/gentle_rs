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

Each bundle also contains `bundle-manifest.json` using
`gentle.gene_set_publication_bundle_manifest.v1`. It binds the exact request
bytes and every finalized HTML, Markdown, resolved-report, data, figure, and
optional PDF artifact by SHA-256 and media type. The manifest deliberately
does not list itself or `generation-report.json`; the generation receipt binds
the manifest path and digest and keeps `copied_files[]` for compatibility with
older consumers.

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
Configured HTML, Markdown, and PDF names must be plain filenames and cannot
replace GENtle's internal receipt paths. If two different source assets would
land at the same case-insensitive bundle path, generation stops with a clear
collision error rather than silently replacing one asset.

With `--pdf`, each external tool must both exit successfully and produce a
non-empty expected output. Errors identify the executable and corresponding
override variable. The HTML, Markdown, and resolved record are still written
without a PDF link when PDF generation cannot complete, so a failed optional
print step cannot masquerade as a complete PDF bundle.
The pre-PDF web manifest remains valid when an optional PDF tool fails; after
a successful PDF build GENtle rewrites the manifest to bind the final PDF and
the updated HTML/resolved report containing its link.

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

## Isoform-assay dossier

The additive `gentle.gene_isoform_assay_publication_request.v1` contract is for
reviewed gene-centred PCR work. Each gene supplies a content-addressed
`gentle.gene_isoform_assay_study_plan.v1`, zero or more
`gentle.experimental_assay_handoff.v1` reports, and zero or more
`gentle.oligo_order_form.v1` records. Order forms are accepted only when every
line retains a matching included handoff digest, readiness-policy id and
digest, assay/pair/oligo identities, and the complete `order_ready` row. The
row must also match the included handoff. This is stricter than copying a broad
`order_ready` label; the same bindings are preserved in each order-sheet TSV.
Each included handoff must in turn name a panel report id, source sequence, and
source feature from one exact `DesignTranscriptAssayPanel` payload in the
approved study plan.

GENtle writes `canonical-report.json`, `index.html`, one `gene-SYMBOL.html`
page per gene, `print.html`, per-gene order-sheet TSVs, and
`projection-report.json`. `--pdf` asks a Chromium-compatible browser to print
that exact `print.html`; it does not maintain a second PDF narrative. GENtle
searches `PATH` and conventional Linux, macOS, and Windows locations for
Chromium, Google Chrome, and Microsoft Edge, while `GENTLE_BROWSER_BIN` remains
the authoritative override. The
canonical report declares `content_blocks[]` and named `profiles[]` (`full`,
`review`, and `ordering`). An outer OpenClaw/ClawBio presentation may select
those ids or host the generated files, but cannot add a block, rewrite a GENtle
value, or treat layout text as scientific evidence.

Planner policy values shared by every gene appear once in the meta-page.
Gene-specific effective values are listed with a source JSON pointer and a
reason. The automatic recommendation and any explicit profile override remain
separate, so an override never erases the recommendation it replaced.

## Adapting another gene set

Change the manifest's gene records, primer table, figure inventory, conditions,
narrative, and provenance. No renderer edit is required. A PARK7-style overview
is therefore routine report structure rather than a one-off design.
