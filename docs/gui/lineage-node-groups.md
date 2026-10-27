### Lineage Node Groups

<!-- manual-counterparts:start -->
No direct CLI counterpart: This chapter describes desktop presentation or interaction; no dedicated command reference is paired with it.
<!-- manual-counterparts:end -->

Node groups in lineage view:

- A `Node Groups` section is available above the lineage table/graph views.
- Groups are disjoint (a node can belong to only one group).
- Each group has one representative node and one or more member nodes.
- Table view:
  - representative rows are shown at top level
  - member rows are indented beneath the representative
  - collapse/expand toggles are available on representative rows
- Graph view:
  - grouped nodes are enclosed with a labeled outline
  - when a group is collapsed, only the representative node is shown
  - edges from/to hidden members are projected to the representative node
- Context menu actions (table and graph):
  - right-click a node to `Mark for node-group` / `Unmark for node-group`
  - right-click a node to `Use as draft representative` (uses currently marked nodes)
  - right-click a node to `Create group now (representative)` for one-step creation

<!-- manual-page-nav:start -->

---

[GUI Manual](../gui.md) | [Previous: Lineage Node Actions](lineage-node-actions.md) | [Next: Open windows and focus](open-windows-and-focus.md)
<!-- manual-page-nav:end -->
