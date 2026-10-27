# Operation History

<!-- manual-counterparts:start -->
Related GUI reference: [Status History And Background Jobs](../gui/status-history-and-background-jobs.md) | [Command Palette and History/Jobs Panels](../gui/command-palette-and-historyjobs-panels.md).
These topic links do not imply identical feature coverage across interfaces.
<!-- manual-counterparts:end -->

Session-local operation history:

- `gentle_cli shell 'history status'`
  - returns the engine-owned multi-level undo/redo summary:
    `undo_count`, `redo_count`, `history_limit`, `operation_log_count`, and
    optional next undo/redo transition labels.
- `gentle_cli shell 'history undo'`
  - restores the previous operation-level engine checkpoint.
- `gentle_cli shell 'history redo'`
  - reapplies the most recently undone operation-level checkpoint.
- The shared parser also accepts `/history`, `/undo`, and `/redo` as concise
  aliases for the three commands above.
- Agent-returned mutating suggestions, including undo/redo, are never executed
  solely by `--allow-auto-exec`; use an explicit execution index or
  `--execute-all` after reviewing the suggestion.
- Undo/redo history is intentionally session-local and is not persisted into
  `.gentle.json` project files.

<!-- manual-page-nav:start -->

---

[CLI Manual](../cli.md) | [Previous: Collections](collections.md) | [Next: Window and Display Intents](display.md)
<!-- manual-page-nav:end -->
