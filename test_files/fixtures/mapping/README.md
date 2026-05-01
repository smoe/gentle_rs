# `test_files/fixtures/mapping`

Compact RNA-seed mapping benchmark fixtures used by deterministic tests.

## Files

- `ensembl_chimp_tp73_all.fasta`
- `ensembl_human_tp73_all.fasta`
- `ensembl_human_tp53_all.fasta`
- `ensembl_human_tp63_all.fasta`
- `ensembl_mouse_trp73_all.fasta`

## Why this folder exists

- Keep committed mapping benchmarks small and CI-friendly.
- Decouple regression tests from larger exploratory TP73 cDNA datasets.
- Include TP53 and TP63 as practical near-family benchmarks for TP73
  specificity and isoform-preflight optimizer controls.

## Provenance

- These files are copied from
  `test_files/mapping/True_TP73/` (legacy mapping corpus).
- Headers and naming indicate Ensembl transcript FASTA exports for
  TP73/TP53/TP63 family genes across species.
- Contributing Ensembl release for these cDNA fixtures: version `115`.
- Exact upstream export URLs and parameter snapshots are not yet captured in
  this repository; this is tracked as a documentation/provenance follow-up in
  `docs/roadmap.md`.

## Current usage in GENtle

- Engine regression tests in `src/engine/tests.rs`:
  - `test_seed_filter_cross_species_and_close_family_specificity_sets`
  - `test_rna_read_preflight_uses_tp53_tp63_variant_controls`
- TP53 and TP63 fixtures are actively loaded as p53-family control transcript
  variants for RNA isoform preflight/optimizer specificity checks.

## Recreation guidance (deterministic target)

Until exact historical export parameters are frozen, regenerate by:

1. Export all transcript FASTA entries for the named genes/species from one
   pinned Ensembl release.
2. Preserve transcript headers (no post-processing).
3. Save with the exact filenames above.
4. Re-run the seed-filter regression to verify behavior parity.
