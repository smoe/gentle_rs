# Screenshot Capture Policy

<!-- manual-counterparts:start -->
Related GUI reference: [Screenshot Capture Policy](../gui/screenshot-capture-policy.md).
These topic links do not imply identical feature coverage across interfaces.
<!-- manual-counterparts:end -->

**Screenshot Capture: CLI Versus GUI**

- **CLI/shared-shell export remains disabled.** `screenshot-window OUTPUT.png`
  returns the deterministic `screenshot-window is disabled by security policy`
  error. `--allow-screenshots` is rejected by argument parsing; enabling the
  `screenshot-capture` build feature does not enable this command.
- **GUI screenshot help is available.** In desktop builds on Linux and macOS,
  `Agent help` requests an egui/eframe capture of that GENtle viewport's content,
  subject to rendering-backend support. This is not desktop capture or the
  reserved CLI command. The image is previewed locally and sent to the selected
  image-capable agent only when the user clicks `Ask Agent`.
- An agent may request a screenshot, but cannot capture one itself. The user
  must select a registered GENtle content window and click `Allow one screenshot`;
  preview and the separate send action still follow. See
  [GUI capture controls](../gui.md#screenshot-capture-policy).
- On macOS 14 or newer, right-clicking `Agent help` also offers complete-window
  capture through ScreenCaptureKit, including native window chrome. This optional
  path requires the `screenshot-capture` feature (enabled in default builds) and
  macOS Screen Recording permission. Normal egui viewport capture does not need
  that permission.
- Automatic/background, headless, whole-desktop and other-application capture
  remain prohibited. For reproducible documentation figures, prefer SVG exports;
  reviewed tutorial screenshots use a separate
  [evidence-and-staging process](../testing.md#62-candidate-bound-tutorial-gate-for-external-auditors),
  not automatic publication of Agent Assistant attachments.

<!-- manual-page-nav:start -->

---

[CLI Manual](../cli.md) | [Previous: Latency Diagnostics](diagnostics.md) | [Next: Agent Bridge](agents.md)
<!-- manual-page-nav:end -->
