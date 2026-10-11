# Compare Sequence Environments of a Repeat Cohort

> Instead of choosing one gene first, can I select repeated elements and compare
> their surroundings without merging different genomic locations?

This 15-minute offline **software exercise** uses three synthetic rmsk rows on
the public bacterial `LocalProject` reference. The LINE/L1 labels test filtering;
they are not evidence of mammalian LINE biology. Read the
[fixture provenance](../../test_files/fixtures/repeat_cohort/README.md).
Then apply the same operations to a genuinely prepared, assembly-matched genome.

## Before You Start

Run from the repository root. Create `/tmp/repeat-lesson` and use a disposable
state and isolated cache. Do not overwrite an earlier cohort or its inputs.
GENtle owns the selection, coordinate projection and motif summaries; there is
no separate tutorial ranking script. The Shell commands below are exercised by
`repeat_cohort_local_project_shell_pipeline_keeps_distinct_loci`. Execution of
that new Rust test remains pending with Glen/CI, not a claimed tutorial replay.

## 1. Prepare and Select

```sh
gentle_cli --state /tmp/repeat-lesson/state.json genomes prepare LocalProject --catalog assets/genomes.json --cache-dir /tmp/repeat-lesson/cache
gentle_cli --state /tmp/repeat-lesson/state.json shell 'features repeat-query LocalProject --rmsk test_files/fixtures/repeat_cohort/local_project.rmsk.txt --alias LINE/L1'
gentle_cli --state /tmp/repeat-lesson/state.json shell 'features repeat-cohort LocalProject --rmsk test_files/fixtures/repeat_cohort/local_project.rmsk.txt --alias LINE/L1 --geometry repeat_midpoint --upstream-bp 10 --downstream-bp 10 --catalog assets/genomes.json --cache /tmp/repeat-lesson/cache --path /tmp/repeat-lesson/cohort.json'
```

Preparation uses retained local FASTA/GenBank inputs, not a remote reference.
BLAST indexing is not required for this motif exercise; retain any preparation
warning rather than interpreting it as a passed specificity check.
The raw file has three rows; the alias selects two. Keep `repName`, `repClass`,
`repFamily`, raw coordinates and strand alongside the normalized label.

The two forward-genomic repeat substrings are both `AGCCAG`, at BED intervals
`[1843,1849)` and `[3877,3883)`. One repeat is plus, one minus. They remain two
rows: identical sequence is not identical locus, and repeat orientation is not
the gene's orientation. Windows include the midpoint and the requested flanks;
read their returned coordinates instead of treating the repeat span as a window.

**Ask the inner agent:** "Prepare LocalProject in an isolated cache, preview
the committed teaching rmsk filter LINE/L1, then propose the repeat-cohort
command. Keep separate loci and explain missing transcript context. Do not
fetch human references or claim these bacterial rows are real LINEs."

## 2. Summarize the Stored Windows

```sh
gentle_cli --state /tmp/repeat-lesson/state.json shell 'features window-cohort-tfbs /tmp/repeat-lesson/cohort.json --motif AAC --score-kind llr_bits --allow-negative --catalog assets/genomes.json --cache /tmp/repeat-lesson/cache --path /tmp/repeat-lesson/tfbs.json'
```

`AAC` is a literal sequence-model teaching motif, not a claimed transcription
factor or JASPAR accession. The expected software checkpoint is **two returned
windows and two summary rows**, with actual score arrays and genomic peak
positions. Zero returned windows with loading warnings is not successful motif
comparison. These scores are model preferences, not occupancy or activity.

## 3. Change Geometry, Not the Selected Repeat Set

Repeat the cohort command into fresh files, changing `--geometry` to
`transcript_5utr_start`, `pol2_promoter_upstream`, then `cds_stop_context`.
All should preserve the two selected rows. This reference's GenBank annotation
does not provide the prepared transcript-index projection: those geometries
remain unavailable, with a reason and no invented anchor. A subsequent summary
skips them with explicit warnings. Missing transcript projection is not absence
of CDS in the source record and not evidence of absent regulation.

In a prepared GTF/GFF-backed genome, inspect associated transcript IDs, strand,
5' boundary, inferred UTR and annotated stop context before selecting a
promoter-oriented frame. Alternative starts are separate biological hypotheses;
the nearest gene is an assignment rule, not a demonstrated regulatory target.

## 4. Move to Real rmsk Data

Select the genome and **the matching assembly's** UCSC RepeatMasker `rmsk`
download. Keep the original file, source URL/date, assembly and SHA-256. Do not
use hg19 repeats with GRCh38 or relabel contigs to hide incompatibility. For
example, with GRCh38 resources independently prepared and reviewed:

```text
features repeat-cohort "Human GRCh38 Ensembl 116" --rmsk data/rmsk.txt.gz --alias LINE/L1 --geometry pol2_promoter_upstream --upstream-bp 2000 --downstream-bp 2000 --limit 50 --path /tmp/repeat-lesson/human-l1.json
features window-cohort-tfbs /tmp/repeat-lesson/human-l1.json --motif SP1 --motif TP73 --score-kind llr_background_tail_log10 --path /tmp/repeat-lesson/human-l1-tfbs.json
```

These are optional resource-dependent commands, not an offline acceptance pass.
Replace paths and catalog/cache with your prepared resources. A 50-row cap is
not a genome-wide survey; record truncation, unavailable contexts and warnings.
Use class/family/name filters for other repeat types rather than assuming LINE
is a privileged class. Defaults without explicit flanks are 2,000 bp each side.

**RNA-read ranking is not implemented for repeat cohorts.** The report's empty
RNA fields are placeholders, not measured zero and not a ranking. Separate RNA
analysis must remain assembly/interval-bound external evidence until a supported
join lands; it must not silently change repeat selection.

## GUI and Interpretation

The deterministic cohort and window summary are the shared Shell/CLI report
path. Do not claim a dedicated cohort GUI simply because JSON was produced.
Open a genuinely anchored sequence and its existing Promoter design or Splicing
expert for a chosen locus; do not invent a row-to-window command or TSS from this
bacterial exercise. [08.16](08-16_tss_regulatory_view_gui.md) shows a validated
per-TSS view; [08.13](08-13_motif_logo_to_promoter_trace.md) explains score units.

Before drawing a common-motif conclusion, ask whether windows use the same
geometry, strand convention, model and units, and whether repeated copies or
duplicated gene assignments dominate the sample. Save inputs and reports with
the exact producer SHA. Native/resource screenshots and scientific approval
remain separate Glen/owner acceptance, not inferred from source checks.
