# Probe-Region Output Fixtures

## `clariom_pm_probe_interpretation/`

- Origin: hand-crafted synthetic helper-output directory for the ClawBio
  Clariom PM-probe interpretation smoke path.
- Deterministic recreation:
  - load `test_files/tp73.ncbi.gb`, whose sequence is anchored to the TP73
    GRCh38.p14 locus on chromosome 1
  - choose tiny intervals inside the first TP73 exon span represented by the
    existing Clariom D TP73 vendor-subset fixture
  - write one probeset/region row and two PM probe rows with fixed synthetic
    sample, condition, and log2 fold-change values
  - mark PM probe rows as `probe_level_input` so `--level pm_probe` projection
    materializes only true probe-level features
- Primary usage:
  - ClawBio gentle-cloning workflow example
    `request_workflow_clariom_pm_probe_interpretation.json`
  - agent-visible smoke coverage for project-then-interpret behavior without
    running APT, R, or using CEL/vendor binary payloads
- Runtime relevance: mirrors a completed
  `arrays import-apt-probe-region-output ... --probe-intensity ...` helper
  directory closely enough for projection and interpretation routing tests; it
  is not biological evidence and is not a substitute for full Clariom
  annotation or CEL-derived analysis.

## `clariom_e_mtab_14704_tp73_validation/`

- Origin: adapter-generated validation fixture for the public E-MTAB-14704
  Clariom D TP73 workflow. It uses a compact Glen-style input derived from
  Glen's exploratory probe-region analysis shape, with TP73-compatible
  hg38/GRCh38.p14 coordinates and public sample labels.
- Deterministic recreation:
  - start from the publication-resource context in
    `data/publication_resources/rostock_p73_clariomd_e_mtab_14704/`
  - do not commit or require raw CEL files
  - regenerate from
    `test_files/fixtures/probe_region_adapter_inputs/clariom_e_mtab_14704_tp73_glen_style.csv`
  - write three parent probeset/region rows and fourteen PM-probe rows with
    nine public sample columns plus condition means and `log2FC_*` values
  - mark every PM-probe row as `probe_level_input`
- Primary usage:
  - release-grade validation of `arrays project-probe-region-output --level
    pm_probe`
  - `arrays interpret-probe-region-evidence --level pm_probe`
  - GUI action-to-capability parity tests
- Runtime relevance: stresses shared-vs-unique transcript geometry,
  junction-spanning review rows, and multiple PM probes per parent region
  without running R/APT or using CEL/vendor binary payloads. It is not
  biological evidence and does not support isoform-support claims.
