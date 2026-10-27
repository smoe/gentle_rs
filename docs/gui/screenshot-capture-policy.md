### Screenshot Capture Policy

<!-- manual-counterparts:start -->
Related CLI reference: [Screenshot Capture Policy](../cli/screenshot-capture.md).
These topic links do not imply identical feature coverage across interfaces.
<!-- manual-counterparts:end -->

- Desktop builds on Linux and macOS support egui/eframe viewport capture,
  subject to rendering-backend support. This captures GENtle's application
  content, not the desktop or another application.
- `Agent help` is an explicit, window-bound diagnostic action. A normal click
  captures only the GENtle viewport containing that button, then opens Agent
  Assistant with a local preview and a prompt for the user's problem report.
- An image-capable inner agent may instead return one typed screenshot request
  with a reason. GENtle captures nothing until the user selects an open,
  registered GENtle content window and clicks `Allow one screenshot`. The
  resulting image is still previewed before the separate `Ask Agent` action.
- Right-clicking `Agent help` on macOS 14 or newer offers an optional complete-
  window capture through ScreenCaptureKit. This requires Screen Recording
  permission and the `screenshot-capture` build feature (enabled by default).
  Normal viewport capture does not require that permission; declining it leaves
  only the optional native full-window path unavailable. GENtle may request a
  restart after permission is granted.
- No screenshot is sent until the user clicks `Ask Agent`. Auto-run cannot
  approve a screenshot; desktop capture and capture of non-GENtle windows are
  prohibited.
- The headless `--allow-screenshots` option and `screenshot-window` shared-shell
  command remain disabled; this narrow GUI help flow does not enable general
  screenshot automation.

macOS auxiliary-window stability note:

- On macOS, GENtle defaults to native OS child viewports with `egui/eframe`
  `0.34.3` because manual testing showed substantially smoother window
  movement, resizing, and DNA-map interaction than the hosted/nested fallback.
- The root GENtle window intentionally does not restore persisted fullscreen or
  maximized state on macOS. Starting the root in macOS fullscreen can trigger
  unstable Split View/tab-like flicker when the next native child window opens.
  If a sequence window is requested while the root is still fullscreen or
  maximized, GENtle first asks macOS to return the root to a regular window and
  opens the sequence window on a later frame.
- The dedicated investigation pack for earlier native-viewport regressions
  lives under `investigations/egui_macos_windowing/` so manual matrices and
  upstream-report notes stay separate from the main GUI manual.
- If native child windows regress, start GENtle with
  `GENTLE_MACOS_HOSTED_CHILD_VIEWPORTS=1 cargo run --bin gentle` to force the
  hosted/nested fallback. If both macOS window-mode environment variables are
  set, `GENTLE_MACOS_NATIVE_CHILD_VIEWPORTS=1` wins and keeps native child
  viewports active.
- In hosted fallback mode, the main/root GENtle window remains native, but it
  acts as a neutral workspace host instead of directly being the project
  biology UI.
- What used to be the privileged main/project surface is now the first hosted
  internal window inside that root workspace, and sequence/help/configuration
  plus other auxiliary workspaces are hosted as sibling internal windows.
- In that fallback, the project/lineage surface should therefore open as one
  draggable/resizable hosted peer window again, rather than filling the root
  host as if maximized.
- Hosted and specialist working windows are expected to route through the shared
  `src/egui_compat.rs` `HostedWindowSpec` / `show_hosted_window` wrapper. The
  wrapper owns stable egui ids, safe-area clamping, foreground-on-focus ordering,
  viewport-builder defaults, and cleanup of stale title-derived layers.
- In hosted mode, the project and root status bars may append short
  hosted-window activity notes while a window is being focused, moved/resized,
  or stale hosted layers are being cleaned up.
- With egui/eframe `0.34.3`, the hosted-window wrapper keeps movement
  restricted to GENtle's title-bar band and disables body-surface movement so
  drags inside one hosted window do not pass through to sibling windows.
- Hosted windows are kept inside a small safety inset from the root workspace
  edges so resizing/moving them does not immediately hit the physical screen
  border and trigger the macOS menu bar or Dock.
- Hosted project and sequence windows also clamp their initial size to that
  inset-safe workspace so resize handles remain reachable and large panes such
  as arrangements and feature details are not born partly off-screen.
- Hosted sequence windows reserve explicit drag room inside the root workspace
  and reset stale over-wide shell state on redraw, so older or over-expanded
  DNA viewers should not become almost immovable after a title-bar drag.
- Native sequence windows start at the same `1200 x 860` working size as hosted
  sequence windows, with an `820 x 520` minimum. This keeps wrapped toolbar
  rows, including `TFBS scan`, reachable on first open instead of clipping them
  behind the sequence client area. Window managers may still clamp the initial
  size to the available monitor.
- The DNA toolbar panel opens at 30% of the available window height (bounded to
  `180..420` px where space permits), may be resized down to 100 px, and keeps
  its own vertical scroll. This gives ordinary first-open access to action rows
  while preserving user control over map space.
- The hosted Configuration window now uses the same safe-area clamping, and
  its Graphics-tab control rows wrap/bound their width instead of expanding
  horizontally without limit.
- In hosted/macOS mode, the Configuration window's `Close` and `Cancel`
  buttons now dismiss the dialog correctly instead of being overridden by the
  embedded window shell state.
- The Help/Tutorial viewport now uses the same regular hosted working-window
  shell as other macOS embedded windows instead of taking over the full child
  viewport surface.
- When older legacy help-window layers are still visible, GENtle now resets
  those stale area states both in the hosted help viewport and in the main/root
  workspace so stray detached help title bars self-heal on reopen.
- Hosted DNA sequence windows and hosted auxiliary workspaces such as
  Splicing Expert and RNA-read Mapping now also clear older legacy
  title-derived embedded layers, so reopening those windows no longer leaves
  detached non-functional title bars behind.
- If eframe itself reports a sequence child viewport as already being an
  embedded window, GENtle now renders the DNA-window content directly into that
  outer shell instead of wrapping it in a second hosted window. This prevents
  the inner sequence client area from moving independently of the visible outer
  frame.
- The root hosted workspace also checks open auxiliary workspace titles and
  clears stale title-derived layers there, covering cases where RNA-read
  Mapping or Splicing Expert is opened or re-owned from inside a DNA sequence
  viewer.
- The hosted RNA-read Mapping shell now renders its workspace intro only once,
  and short saved-report preview tables expand a little further before opening
  an inner vertical scrollbar so wheel scrolling is less likely to get trapped
  in a tiny nested pane.
- The hosted help body now reserves the full remaining window height before its
  vertical scroll area is built so mouse-wheel scrolling stays with the help
  content instead of falling through to the background workspace.
- Hosted dialog/project/help content now paints its backdrop against the
  visible shell clip and reserves an explicit client-area frame inside the
  hosted window, reducing background spill and making resize/scroll ownership
  less ambiguous.
- While this hosted-window workaround is active, native macOS `Main` / `Windows`
  menu behavior can be imperfect because those menus still assume native child
  windows. The menu path is intentionally not being deepened further; the
  intended exit is to return to native child viewports once the upstream egui
  bug is fixed.
- Newly opened sequence windows are queued for foreground focus after their
  viewport is registered, so opening a sequence from a project/lineage view
  should raise that DNA viewer instead of leaving it behind older windows.
- Opening a saved sequence does not recompute or persist derived restriction
  caches in shared project state. Restriction scans remain explicit shared
  operations, so merely viewing a sequence does not mark its project modified.
- Hosted auxiliary workspaces opened from a sequence window, including
  Splicing Expert, RNA-read Mapping, Dotplot, Promoter design, and Isoform /
  Protein Expert, render in a root-level hosted pass after the owning sequence
  body. This avoids a second embedded eframe title shell and keeps the
  auxiliary workspace above the large DNA viewer instead of trapped behind or
  inside it.
- In hosted/macOS mode, `Windows` menu focus for a DNA-owned auxiliary
  workspace keeps the owner visible but does not queue that owner as the
  foreground shell; the auxiliary target itself gets the foreground layer.
- Simple one-shot modal prompts such as unsaved-changes, removal confirmation,
  cache cleanup, operation failure, and prepared-genome choice dialogs route
  through `src/egui_compat.rs` `ModalWindowSpec` / `show_modal_window`. The modal
  wrapper owns stable egui ids, foreground ordering, centered placement, and
  safe-area clamping for prompt-sized shells.

Manual GUI stability checklist for macOS window modes:

- Start with default macOS native child viewports. Keep the root GENtle window
  out of macOS fullscreen before opening the first project/DNA child window.
- For the hosted fallback checklist, start with
  `GENTLE_MACOS_HOSTED_CHILD_VIEWPORTS=1 cargo run --bin gentle`.
- Resize the root window down to a small usable viewport and verify project,
  sequence, Help, Configuration, Background Jobs, and specialist windows keep
  reachable title bars and resize handles.
- Open Help, Configuration, a sequence window, and one specialist window; reopen
  each from its menu route and verify the existing window comes forward without
  resetting the current tab/scroll/input state.
- Drag a large DNA sequence window by its title bar after resizing it close to
  the root workspace bounds; it should move as a whole shell and still leave
  visible workspace room on either side.
- Exercise modal prompts: unsaved-changes, About, removal confirmation,
  operation failure, prepared-genome choice, and cache cleanup close/cancel
  buttons should dismiss the prompt they own.
- Use `Window -> GENtle Open Windows...` or `GENtle -> GENtle Windows...` to
  focus a background hosted window and verify it rises above the project
  workspace.

<!-- manual-page-nav:start -->

---

[GUI Manual](../gui.md) | [Previous: Catalog Inspection And Provider Imports](catalog-inspection-and-provider-imports.md) | [Next: External Services Window](external-services-window.md)
<!-- manual-page-nav:end -->
