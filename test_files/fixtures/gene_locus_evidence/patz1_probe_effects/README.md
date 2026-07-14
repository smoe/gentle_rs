# PATZ1 Clariom probe-effect fixture

This directory contains a compact PATZ1-only fixture derived from the local
E-MTAB-14704 Clariom D probe-set activity summary used during the PATZ1
CUT&RUN evidence-display discussion. It exists so agents and tests can develop
GENtle's gene-locus probe-effect overlay without depending on Glen's untracked
`analysis/` directory or on bulky CEL-derived intermediate tables.

The fixture is intentionally small and keeps only fields needed for coordinate
projection, rendering, and sanity checks:

- `patz1_clariom_probe_effects.tsv`: one row per PATZ1 PSR or JUC probe set.
- `patz1_clariom_probe_effects_summary.tsv`: gene-level paired-contrast
  summary for PATZ1.

## Provenance

Source table on Glen's workstation:

`analysis/e_mtab_14704_tp73_microarray/gene_panel_probe_set_activity/probe_set_activity_summary.tsv`

Gene-level source table:

`analysis/e_mtab_14704_tp73_microarray/gene_panel_probe_set_activity/paired_gene_level_summary.tsv`

The source tables summarize raw PM-probe activity from SK-MEL-29 cells
overexpressing DNp73beta, TAp73alpha, or GFP. This fixture is for visualization
and software tests. It must not be interpreted as a formal differential
expression model, isoform-support call, or statistical significance result.

## Biological sanity checks

For PATZ1 in this compact fixture:

- probe/junction probe sets: 24
- PM probes: 130
- median DNp73beta-GFP across paired contrasts: -0.0671
- median TAp73alpha-GFP across paired contrasts: -0.3833
- mean DNp73beta-GFP across probe-set rows: about +0.1384
- mean TAp73alpha-GFP across probe-set rows: about -0.5729
- TAp73alpha-GFP is negative for all 24 rows
- DNp73beta-GFP is positive for 18 of 24 rows

PATZ1 is on the negative strand. Rows are listed in genomic ascending order, not
transcript 5'-to-3' order. Implementations that render transcript-relative
tracks should handle the strand explicitly.
