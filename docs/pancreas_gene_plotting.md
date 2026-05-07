# Pancreas Gene-Screen Plotting Scripts

This document describes the small Python plotting layer used for pancreatic
Nanopore cDNA cohort summaries. The scripts are intentionally dependency-free:
they write SVG directly and do not require matplotlib, R, or a graphical
session. They do not re-run RNA-read interpretation or alignment; they only
consume already-produced GENtle report tables.

The current plotting split is deliberately simple:

- `scripts/plot_pancreas_gene_screen.py`: canonical single-gene plotter.
- `scripts/plot_pancreas_gene_family.py`: canonical gene-family plotter.
- `scripts/plot_tp73_pancreas_cohort.py`: TP73 compatibility wrapper.
- `scripts/plot_p53_family_pancreas.py`: p53-family compatibility wrapper.

The wrappers exist so older runbooks and command histories keep working. New
commands should prefer the two canonical scripts.

## Single-Gene Plotter

Use `scripts/plot_pancreas_gene_screen.py` after a gene screen has been
summarized with `scripts/pancreas_gene_rna_screen.sh summarize RUN_ROOT`.
The input is the generated figure-source TSV:

```text
RUN_ROOT/figures/<gene>_pancreas_figure_source.tsv
```

Typical command:

```bash
python3 scripts/plot_pancreas_gene_screen.py \
  "$RUN_ROOT/figures/e2f1_pancreas_figure_source.tsv" \
  --gene E2F1 \
  --bar-column strict_seed_passed_reads \
  --output "$RUN_ROOT/figures/e2f1_pancreas_overview.svg"
```

Useful options:

- `--gene GENE`: title/label gene name. If omitted, the script infers a leading
  gene token from the input file name.
- `--bar-column strict_seed_passed_reads`: conservative abundance view and the
  default, preferred setting for release comparisons. This is the number of
  reads that pass the first seed phase before sequence-alignment confirmation.
- `--bar-column strict_seed_accepted_target_reads`: strict seed-passed rows that
  were accepted as target-supporting by the downstream report contract after
  the confirmation stage.
- `--bar-column accepted_target_reads`: accepted target rows when this column is
  present in a generic figure-source TSV. This is a downstream accepted view,
  not the raw first-pass seed count.
- `--bar-column accepted_tp73_reads`: legacy TP73 figure-source alias retained
  for older cohort tables; it is likewise a downstream accepted view.
- `--label-column sample_id|sample_name|run_accession`: x-axis labels.
- `--title TEXT`: override the generated title.

The output SVG contains two panels:

- whole-library read-length context from all-read quantiles/mean
  (`all_q0_bp`, `all_q25_bp`, `all_q50_bp`, `all_q75_bp`, `all_q90_bp`,
  `all_q95_bp`, `all_q99_bp`, `all_q100_bp`, `all_mean_bp`), drawn on a
  logarithmic read-length axis;
- the selected single-gene support count as bars, with supporting-read maximum
  length overlaid when the relevant length column is available.

The script warns on stderr if the all-read quantiles are non-monotone, because
that usually indicates a TSV extraction or column-mapping problem.

## Gene-Family Plotter

Use `scripts/plot_pancreas_gene_family.py` when several gene screens should be
shown as neighboring bars for each sample. It first writes a canonical merged
family TSV, then renders an SVG. This makes the figure inputs auditable and
lets later scripts reuse the normalized table.

Inputs are supplied as repeatable `GENE=PATH` pairs. The script accepts three
source styles:

- `--batch-report GENE=PATH`: an RNA-read batch `out/batch_report.json`.
- `--batch-summary GENE=PATH`: an RNA-read batch `out/batch_summary.tsv`.
- `--figure-source GENE=PATH`: a `pancreas_gene_rna_screen.sh summarize`
  figure-source TSV.

Example for the p53 family:

```bash
python3 scripts/plot_pancreas_gene_family.py \
  --batch-report TP53="$TP53_ROOT/out/batch_report.json" \
  --figure-source TP63="$TP63_ROOT/figures/tp63_pancreas_figure_source.tsv" \
  --figure-source TP73="$TP73_ROOT/figures/tp73_pancreas_figure_source.tsv" \
  --genes TP53,TP63,TP73 \
  --metric per_million \
  --output "$WORK/figures/p53_family_seed_passed_grouped.svg" \
  --output-tsv "$WORK/figures/p53_family_seed_passed_grouped.tsv"
```

Useful options:

- `--genes TP53,TP63,TP73`: explicit display order. If omitted, input order is
  used and any extra genes found in the rows are appended.
- `--metric per_million`: normalize bars as seed-passed reads per million total
  reads. Raw read counts are still printed above bars.
- `--metric raw`: use raw seed-passed read counts as bar height.
- `--support-length-stat max|mean`: choose the gene-supporting read-length
  statistic drawn as one colored line per gene in the lower support panel. The
  default is `max`, matching the single-gene overview convention. Whole-library
  q90 remains library context only and is drawn in the upper panel when present.
- `--support-length-source strict_seed_passed|any`: choose which support-length
  provenance can be drawn. The default is `strict_seed_passed`, so the lower
  read-length lines describe the same seed-passed read population as the bars.
  `any` is a diagnostic mode for older reports and can draw accepted-target
  fallback lengths, which are not directly comparable to strict seed-passed bars.
- `--support-length-scale log|linear`: choose the lower read-length axis scale.
  The default is `log`, matching the whole-library read-length panel and making
  TP63/TP73-style length differences easier to compare without outlier
  domination.
- `--support-length-display points|lines`: choose whether support-read lengths
  are drawn as independent symbols or connected lines. The default is `points`,
  because sparse seed-passed samples should not imply a continuous trend across
  missing samples.
- `--support-length-genes TP63,TP73`: optionally draw support-read length
  symbols for only a subset of genes while keeping all gene bars. This is useful
  when one older input lacks strict seed-passed length summaries or only has
  downstream accepted-target fallback lengths.
- `--show-all-read-max`: include the whole-library `all_q99_bp` and
  `all_q100_bp`/maximum-read lines in the upper read-length panel. They are
  hidden by default in family plots because extreme read lengths are outlier
  statistics and tend to dominate the visual impression without helping the
  gene-support comparison. `all_q95_bp` is collected and displayed by default as
  the calmer upper-tail context line.
- `--label-column sample_id|sample_name|run_accession`: x-axis labels.
- `--output-tsv PATH`: canonical merged family table. Defaults to the SVG path
  with `.tsv` extension.
- `--title TEXT`: override the generated title.

The canonical family TSV uses these main fields:

- `gene`, `run_accession`, `sample_id`, `sample_name`
- `total_reads`
- all-read length fields (`all_q0_bp`, `all_q25_bp`, `all_q50_bp`,
  `all_q75_bp`, `all_q90_bp`, `all_q95_bp`, `all_q99_bp`, `all_q100_bp`, plus
  `all_mean_bp`)
- `seed_passed_reads` and `seed_passed_per_million`
- `accepted_target_count` and `accepted_target_per_million`
- optional supporting-read length fields
- `source`, pointing back to the input table/report

The grouped SVG keeps `seed_passed_reads` as the primary support measure. For a
three-gene p53-family comparison, every sample therefore has three neighboring
bars. `accepted_target_count` values are retained in the TSV for auditability but
are not used as the main family-support scale, because those counts can mean a
later, less conservative downstream contract than the strict seed gate.
Gene-supporting read lengths are drawn as colored symbols over the corresponding
gene bars in the same lower panel. By default, this overlay uses the selected
support statistic (`max` by default) only from strict seed-passed rows, so it
describes the same read population as the bars. Accepted-target fallback lengths
are retained in the TSV for older batch reports, but are drawn only when
`--support-length-source any` is selected explicitly. Connecting symbols with
`--support-length-display lines` is available as a diagnostic view, but the
release-facing default is unconnected points because missing samples should stay
visibly missing. The all-read q90 line belongs only to the upper whole-library
read-length panel.

## Compatibility Wrappers

### TP73 Wrapper

`scripts/plot_tp73_pancreas_cohort.py` preserves the old TP73-specific command
line and `$COHORT_ROOT` default input. It forwards to
`scripts/plot_pancreas_gene_screen.py` with `--gene TP73`.

Legacy command:

```bash
python3 scripts/plot_tp73_pancreas_cohort.py \
  "$COHORT_ROOT/figures/tp73_pancreas_figure_source.tsv" \
  --bar-column strict_seed_passed_reads \
  --output "$COHORT_ROOT/figures/tp73_pancreas_cohort_overview.svg"
```

Equivalent canonical command:

```bash
python3 scripts/plot_pancreas_gene_screen.py \
  "$COHORT_ROOT/figures/tp73_pancreas_figure_source.tsv" \
  --gene TP73 \
  --bar-column strict_seed_passed_reads \
  --output "$COHORT_ROOT/figures/tp73_pancreas_cohort_overview.svg"
```

### p53-Family Wrapper

`scripts/plot_p53_family_pancreas.py` preserves the old TP53/TP63/TP73 command
line and forwards to `scripts/plot_pancreas_gene_family.py` with
`--genes TP53,TP63,TP73`.

Legacy command:

```bash
python3 scripts/plot_p53_family_pancreas.py \
  --tp53-batch-report "$TP53_ROOT/out/batch_report.json" \
  --tp63-figure-source "$TP63_ROOT/figures/tp63_pancreas_figure_source.tsv" \
  --tp73-figure-source "$TP73_ROOT/figures/tp73_pancreas_figure_source.tsv" \
  --output "$WORK/figures/p53_family_seed_passed_grouped.svg"
```

Equivalent canonical command:

```bash
python3 scripts/plot_pancreas_gene_family.py \
  --batch-report TP53="$TP53_ROOT/out/batch_report.json" \
  --figure-source TP63="$TP63_ROOT/figures/tp63_pancreas_figure_source.tsv" \
  --figure-source TP73="$TP73_ROOT/figures/tp73_pancreas_figure_source.tsv" \
  --genes TP53,TP63,TP73 \
  --metric per_million \
  --support-length-stat max \
  --support-length-source strict_seed_passed \
  --support-length-scale log \
  --support-length-display points \
  --support-length-genes TP63,TP73 \
  --output "$WORK/figures/p53_family_seed_passed_grouped.svg"
```

## Column Expectations

For single-gene figure-source TSVs, the most useful columns are:

```text
run_accession
sample_id
sample_name
total_reads
all_q0_bp
all_q25_bp
all_q50_bp
all_q75_bp
all_q90_bp
all_q95_bp
all_q99_bp
all_q100_bp
all_mean_bp
strict_seed_passed_reads
strict_seed_passed_max_read_bp
accepted_target_reads or accepted_tp73_reads
accepted_target_max_read_bp
```

The gene-family plotter normalizes differences between batch JSON, batch TSV,
and figure-source TSV inputs. When multiple gene rows share one SRA accession,
all-read length context is copied from the row with the richest read-length
fields so every gene can be plotted on the same library axis.

## Recommended Release Figure Practice

For release-facing comparisons, prefer:

- `strict_seed_passed_reads` for single-gene abundance bars;
- `--metric per_million` for cross-sample family plots;
- explicit `--genes ...` order for family plots;
- committing or archiving the generated canonical family TSV beside the SVG;
- treating wrapper scripts as compatibility only, not as places for new plotting
  behavior.
