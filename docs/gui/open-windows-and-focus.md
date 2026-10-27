## Open windows and focus

<!-- manual-counterparts:start -->
Related CLI reference: [Window and Display Intents](../cli/display.md).
These topic links do not imply identical feature coverage across interfaces.
<!-- manual-counterparts:end -->

GENtle tracks open native windows and can raise a selected one to front.

- Main project window menu: `Windows -> <window name>` jumps directly to that window
- macOS native mirrors:
  - `Window -> GENtle Open Windows… -> <window name>`
  - `GENtle -> GENtle Windows… -> <window name>`
- Native macOS window-menu synchronization is deferred while a new Help or
  Configuration open probe is active to reduce first-frame contention.
- `Windows` includes project, sequence/pool, and auxiliary windows
  (Help, Configuration, Prepare Genome, Retrieve, BLAST, Track Import,
  Planning, Agent Assistant, Protein Evidence, Operation History)
- In hosted/macOS mode, selecting a listed window now also raises the embedded
  egui child window to the front instead of only requesting viewport focus.
- For DNA-owned auxiliary windows, the sequence owner is kept visible but the
  auxiliary shell itself receives the foreground layer so it does not remain
  hidden behind the DNA viewer.
- Shortcut: `Cmd+Backtick` focuses the main project window
- Specialist windows (including DNA sequence windows) include a top-left nav
  strip with:
  - `Help`: opens the GUI manual
  - `Main`: raises the main project window
- `Cmd+W` closes the focused native window (sequence and auxiliary viewports,
  including Configuration/Help/BLAST/Track Import/Planning/Agent Assistant/History).
- Closing a sequence window no longer force-closes an already-open
  `RNA-read Mapping` auxiliary workspace; close that workspace directly when
  you are done with it.
- Help shortcuts:
  - `F1` (Windows/Linux)
  - `Ctrl+F1` (fallback in function-key reserved environments)
  - `Cmd+Shift+/` (macOS)

<!-- manual-page-nav:start -->

---

[GUI Manual](../gui.md) | [Previous: Lineage Node Groups](lineage-node-groups.md) | [Next: Help manuals](help-manuals.md)
<!-- manual-page-nav:end -->
