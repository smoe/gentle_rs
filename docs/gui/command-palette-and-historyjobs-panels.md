## Command Palette and History/Jobs Panels

<!-- manual-counterparts:start -->
Related CLI reference: [Operation History](../cli/history.md).
These topic links do not imply identical feature coverage across interfaces.
<!-- manual-counterparts:end -->

Subject-bound DNA tools use the focused sequence viewer or an explicitly
selected sequence node in the main project graph, never the first sequence by
name. PCR Designer, Sequencing Confirmation, TATA evidence, cryptic-splicing and
saved-region/conservation launches reject protein/RNA targets and explain how to
select DNA. Feature Location Editor requires an explicit sequence but is not
restricted to DNA. Gibson keeps its input choosers: the selected DNA may prefill
the destination, but an insert is never picked automatically.

The palette captures its initiating subject before taking focus. Closing that
viewer, changing the project or removing the sequence invalidates the subject
rather than choosing another one. Migrated actions are dimmed when unavailable
and explain their prerequisites; broader action coverage and native macOS menu
presentation remain separate follow-ups. Availability means the advertised
launcher can be used, not that its analysis passed scientific validation.

- Command Palette:
  - open via `Cmd/Ctrl+K` or `Edit -> Command Palette...`
  - searchable action list for project, genome, planning, help, and window actions
  - UI-intent entries for reference genomes, helper genomes, PCR, sequencing,
    and agent assistant are generated from the shared `UiIntentTarget`
    metadata used by shell/agent capability discovery
  - includes `Planning` action (`Patterns -> Planning...` equivalent)
  - compact result rows show action names; a separate detail area explains the
    selected action's purpose and any unavailable reason/recovery guidance
  - hovering a disabled row or selecting it with `Up`/`Down` shows its explanation
    without enabling it; a stationary pointer does not override keyboard selection
  - details remain outside the result-list scroll, with their own scroll for long
    guidance in small windows; search still matches descriptions and keywords
  - supports keyboard navigation (`Up`/`Down`, `Enter`, `Esc`), with the selected
    result brought into view on arrow navigation
  - Undo/Redo use the same current session-history and background-job checks as
    the egui Edit menu, including a fresh check when invoked
  - Edit-menu checkpoint counts use the same nonblocking live snapshot as its
    buttons, not the status-bar cache. A busy project reports unavailable
    history instead of stale counts; the history limit remains in Operation History
- Operation History panel:
  - open via `Edit -> Operation History...` or `Window -> Show Operation History`
  - includes undo/redo buttons, available counts, history limit, next transition
    labels, and recent operation summaries
  - uses the same engine-owned history summary surfaced by shell/CLI
    `history status`
  - standalone no-path design previews in Agent Assistant retain bounded command
    results/receipts instead of adding project History rows; workflows, file-writing
    previews and raw CLI/MCP/scripting keep their existing journal behavior
- Background Jobs panel:
  - open via `Window -> Show Background Jobs`
  - central place for progress, cancel/retry actions, and recent completion/error events
  - retry actions capture structured argument snapshots (shown in-panel) and
    persist them with project metadata for restart-safe debugging
  - retry snapshot list supports kind/text filtering and filtered JSON export
    for triage handoff or offline analysis
  - retry snapshot retention controls support "retain newest N", explicit prune,
    and clear-all cleanup for long-running sessions
  - filtered maintenance actions now include direct filtered delete plus
    archive-and-delete export workflow for per-kind/per-origin cleanup
  - destructive filtered cleanup actions are staged and require explicit
    confirm with a match-summary preview before execution
  - staged destructive actions now include a dry-run diff panel showing
    "would remove" and "would remain" snapshot previews before confirm
  - destructive cleanup confirm now requires typing an action-specific
    phrase (for example `delete 3`) before the confirm button is enabled
  - successful destructive cleanup actions append persisted audit entries
    (action/filter/counts/archive path) visible in-panel for traceability
  - cleanup audit section provides filtered JSON report export using current
    action/text filters; report export is read-only and does not append
    additional audit entries
  - cleanup audit history supports action/text filtering plus independent
    retention tuning (`retain newest N`, `prune oldest`, `clear all`) for
    bounded, searchable long-running session traces
  - cleanup audit `clear all` now uses staged type-to-confirm before final
    removal

<!-- manual-page-nav:start -->

---

[GUI Manual](../gui.md) | [Previous: Help manuals](help-manuals.md) | [Next: Engine Ops State Persistence](engine-ops-state-persistence.md)
<!-- manual-page-nav:end -->
