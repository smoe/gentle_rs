## About GENtle

<!-- manual-counterparts:start -->
No direct CLI counterpart: This chapter describes desktop presentation or interaction; no dedicated command reference is paired with it.
<!-- manual-counterparts:end -->

The in-app `Help -> About GENtle` entry and the native macOS application-menu
entry open the same GENtle-owned About experience. The native label follows the
selected presentation language, including `GENtle -> Über GENtle` in German;
it no longer opens AppKit's sparse standard metadata panel.

The compact, Finder-inspired About card shows the application icon, purpose,
internal-preview status, semantic package version, build stamp, authorship, and
licence. `More Info …` / `Weitere Infos …` opens a larger mascot view with:

- the shared-engine and provenance model
- the internal-preview bench-use caution
- package version, build identity, and licence
- direct links to the GitHub project, documentation, and credits
- a copy action for reproducible build information

The copied payload and CLI `--version` draw from the same compile-time About
metadata. Windows and Linux use the same custom GUI panels; only the native
application-menu bridge is macOS-specific.

<!-- manual-page-nav:start -->

---

[GUI Manual](../gui.md) | [Previous: Agent Assistant](agent-assistant.md) | [Next: Engine Settings (Engine Ops)](engine-settings-engine-ops.md)
<!-- manual-page-nav:end -->
