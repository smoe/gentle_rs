### Status History And Background Jobs

<!-- manual-counterparts:start -->
Related CLI reference: [Operation History](../cli/history.md).
These topic links do not imply identical feature coverage across interfaces.
<!-- manual-counterparts:end -->

Global productivity controls:

- Status bar (bottom of main window) shows:
  - hovered control stable name (when available)
  - undo/redo availability counters
  - latest app/job status message
    (including slow-open timing diagnostics for help/configuration windows:
    native Help-menu dispatch, configuration runtime sync, help payload load,
    focus-acquisition latency, first-frame render, total open latency, and
    native window-menu sync when those exceed threshold)
- `Edit -> Undo` / `Redo` and `Window -> Show Operation History` expose
  operation-level, multi-level history controls.
  - The history stack is session-local: it is not written into project files
    and is reset after restart/reload.
  - Undo/redo remains disabled while background jobs are active, so long-running
    imports or analyses cannot race against state restoration.
- `Window -> Show Background Jobs` opens a centralized progress panel for
  long-running tasks (prepare/import/BLAST).
  - Engine-backed dbSNP fetch, track import, primer design, TFBS, and RNA-read
    workers run against detached snapshots, so the sequence windows remain
    inspectable while computation or network work is active.
  - If the project changes before a mutating worker finishes, GENtle reports
    the result as stale and asks for a rerun instead of overwriting the newer
    project state.
- Sequence layer badges, custom-bin GC tracks, lineage graph projection/layout,
  JASPAR filtering, and sequence-window titles are cached by their relevant
  revisions. Large JASPAR and feature-tree lists render only visible rows.
- `Cmd/Ctrl+K` opens the Command Palette.

<!-- manual-page-nav:start -->

---

[GUI Manual](../gui.md) | [Previous: Project And Sequence Window Examples](project-and-sequence-window-examples.md) | [Next: Lineage Node Actions](lineage-node-actions.md)
<!-- manual-page-nav:end -->
