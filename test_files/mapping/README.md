# `test_files/mapping`

Exploratory RNA mapping datasets and benchmark corpora.

## Layout

- `True_TP73/`
  - positive-control TP73-family transcript corpus (legacy benchmark set).
- `False_TP73/`
  - decoy/control transcript corpus for TP73 specificity experiments.
- `SRR32957124*.fasta.gz`
  - larger read datasets used for local mapping/optimization runs.

## Manual provisioning (local-only corpora)

- Files in this directory are treated as local exploratory inputs and may not
  be present in every checkout/environment.
- Provide/copy these files manually when running large-corpus mapping
  benchmarks.
- Because of size/runtime/provenance constraints, these datasets are currently
  **not part of CI**.

## How to test (CI vs local)

1. CI-safe deterministic regression (committed compact fixtures):

   ```bash
   cargo test -q test_tp73_seed_filter_cross_species_and_tp53_specificity_sets
   ```

   This test uses `test_files/fixtures/mapping/` (not this folder).

2. Local exploratory mapping run (manual datasets in this folder):

   - Start from the workflow template
     `docs/examples/workflows/rna_reads_interpret_cdna_tp73_template.json`.
   - Update at least:
     - `input_path` -> one of:
       - `test_files/mapping/SRR32957124_10000.fasta.gz`
       - `test_files/mapping/SRR32957124_100000.fasta.gz`
       - `test_files/mapping/SRR32957124.fasta.gz`
     - `report_id` -> a unique local ID.
     - `seq_id` and `seed_feature_id` -> values valid for your local project
       state.
   - Run:

   ```bash
   cargo run --bin gentle_cli -- shell 'workflow @docs/examples/workflows/rna_reads_interpret_cdna_tp73_template.json'
   ```

   - Optional report exports:

   ```bash
   cargo run --bin gentle_cli -- shell 'rna-reads list-reports'
   cargo run --bin gentle_cli -- shell 'rna-reads export-score-density-svg REPORT_ID test_files/mapping/REPORT_ID_score_density.svg --scale log'
   cargo run --bin gentle_cli -- shell 'rna-reads export-abundance-tsv REPORT_ID test_files/mapping/REPORT_ID_exon_abundance.tsv'
   cargo run --bin gentle_cli -- shell 'rna-reads export-paths-tsv REPORT_ID test_files/mapping/REPORT_ID_exon_paths.tsv'
   ```
