## Engine Ops State Persistence

<!-- manual-counterparts:start -->
Related CLI reference: [Project and Sequence State](../cli/projects.md).
These topic links do not imply identical feature coverage across interfaces.
<!-- manual-counterparts:end -->

Engine Ops panel input state is persisted in project metadata per active
sequence id. This includes panel visibility and operation-form inputs.

The same persistence payload now also stores:

- shell panel visibility
- last shell command text

Metadata key format:

- `gui.engine_ops.<seq_id>`

<!-- manual-page-nav:start -->

---

[GUI Manual](../gui.md) | [Previous: Command Palette and History/Jobs Panels](command-palette-and-historyjobs-panels.md) | [Next: Loading and creating sequences](loading-and-creating-sequences.md)
<!-- manual-page-nav:end -->
