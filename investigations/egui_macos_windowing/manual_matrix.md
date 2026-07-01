# macOS egui Manual Matrix

Use `cargo run --bin gentle_egui_window_repro` and switch the child-window mode
from the combo box inside the repro app.

## Native child viewport modes

| Case | Mode | Steps | What to watch for |
| --- | --- | --- | --- |
| Immediate open/close | `immediate native viewport` | Open one child, close it from the button, reopen it. | stale twins, reappearing closed windows, focus jumps |
| Immediate resize | `immediate native viewport` | Open one child, resize root and child repeatedly. | duplicated surfaces, lost repaint, wrong title shell |
| Immediate maximize | `immediate native viewport` | Open one child and maximize/restore it several times. | respawned child, blank content, stuck maximize state |
| Deferred open/close | `deferred native viewport` | Repeat the immediate open/close case. | mode-specific lifecycle drift |
| Deferred resize/maximize | `deferred native viewport` | Repeat the immediate resize/maximize cases. | close/focus regressions unique to deferred mode |
| Nested child chain | `immediate native viewport` or `deferred native viewport` | Open child from root, then open another child from that child. | parent/child confusion, serial reuse, stale viewport IDs |

## Embedded hosted baseline

| Case | Mode | Steps | What to watch for |
| --- | --- | --- | --- |
| Embedded baseline | `embedded egui::Window` | Open several children, drag them, resize them, close them. | whether the bug disappears when native child viewports are removed |
| Embedded nested chain | `embedded egui::Window` | Open child from child a few times. | ordering/focus behavior in the fallback path |

## Main-app comparison

Run the same general interactions in the full app when you need to compare the
minimal harness with GENtle's real window stack:

- Hosted fallback on macOS:
  `cargo run --bin gentle`
- Native child viewport override:
  `GENTLE_MACOS_NATIVE_CHILD_VIEWPORTS=1 cargo run --bin gentle`

Compare:

- whether the issue reproduces in the minimal repro
- whether it only reproduces in GENtle after hosted-window adaptation
- whether it only appears under maximize/restore instead of plain resize
- whether focus ordering bugs exist in both the repro and the main app
