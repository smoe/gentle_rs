# Construct Reasoning

<!-- manual-counterparts:start -->
No direct GUI counterpart: This chapter documents an automation interface, example or command-only workflow; no dedicated GUI chapter is paired with it.
<!-- manual-counterparts:end -->

Construct-reasoning inspection capability status:

- `gentle_cli`: supported via shared shell/direct commands (`construct-reasoning list-graphs`, `construct-reasoning show-graph`, `construct-reasoning list-inspection-actions`, `construct-reasoning run-inspection-action`, `construct-reasoning set-annotation-status`, `construct-reasoning write-annotation`, `construct-reasoning export-graph`)
- `gentle_js`: supported via dedicated shared-shell-backed helpers (`list_construct_reasoning_graphs`, `show_construct_reasoning_graph`, `list_construct_reasoning_inspection_actions`, `run_construct_reasoning_inspection_action`, `set_construct_reasoning_annotation_status`, `write_back_construct_reasoning_annotation`)
- `gentle_lua`: supported via dedicated shared-shell-backed helpers (`list_construct_reasoning_graphs`, `show_construct_reasoning_graph`, `list_construct_reasoning_inspection_actions`, `run_construct_reasoning_inspection_action`, `set_construct_reasoning_annotation_status`, `write_back_construct_reasoning_annotation`)
- `gentle_mcp`: supported via thin tool wrappers over the same shared shell contracts (`construct_reasoning_graphs`, `construct_reasoning_graph`, `construct_reasoning_inspection_actions`, `construct_reasoning_run_inspection_action`, `construct_reasoning_set_annotation_status`, `construct_reasoning_write_annotation`)

<!-- manual-page-nav:start -->

---

[CLI Manual](../cli.md) | [Previous: Figure and Layout Exports](rendering.md) | [Next: Variant and Reporter Evidence](variants-reporters.md)
<!-- manual-page-nav:end -->
