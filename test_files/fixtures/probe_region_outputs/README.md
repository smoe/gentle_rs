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
