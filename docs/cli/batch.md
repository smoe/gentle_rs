# Batch Workflows

<!-- manual-counterparts:start -->
No direct GUI counterpart: This chapter documents an automation interface, example or command-only workflow; no dedicated GUI chapter is paired with it.
<!-- manual-counterparts:end -->

Batch helpers:

- `batch plan MANIFEST.tsv --template WORKFLOW.json [--out-dir DIR] ...`
  expands a manifest into one workflow JSON file per row plus a runnable local
  shell script or Slurm array script. The batch layer only substitutes
  placeholders such as `${sample_id}`, `${input_path}`, `${row_id}`, and
  `${out_dir}` into an existing workflow template; all biology still runs
  through normal engine workflows.
- Use `--state-template STATE.json --state-mode copy-per-row` for cluster
  jobs. Each row then receives its own copied project state, avoiding
  concurrent writes to one mutable `.gentle_state.json`.
- `batch run MANIFEST.tsv --template WORKFLOW.json [--limit N]` runs the same
  expansion locally and sequentially against the current state. It is intended
  for smoke tests before submitting the generated plan to a cluster.

<!-- manual-page-nav:start -->

---

[CLI Manual](../cli.md) | [Previous: Window and Display Intents](display.md) | [Next: Macros and Routines](macros-routines.md)
<!-- manual-page-nav:end -->
