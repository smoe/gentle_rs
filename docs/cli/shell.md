# Shared Shell

<!-- manual-counterparts:start -->
Related GUI reference: [GENtle Shell (GUI)](../gui/gentle-shell-gui.md).
These topic links do not imply identical feature coverage across interfaces.
<!-- manual-counterparts:end -->

Shared shell command:

- `shell '<command line>'`
  - Parses and executes commands through the same shared parser/executor used by
    the GUI `Shell` panel (`src/engine_shell.rs`).
  - Shell-supported commands:

- Shared shell:
  - Commands:
    - `op <operation-json-or-@file>`
    - `workflow <workflow-json-or-@file>`
    - `screenshot-window OUTPUT.png` (reserved CLI/shared-shell export; disabled by security policy)
    - `help [COMMAND ...] [--format text|json|markdown] [--interface all|cli-direct|cli-shell|gui-shell|gui-menu|js|lua|mcp]`
  - Use single quotes around JSON payloads to preserve whitespace:
    - `gentle_cli shell 'workflow {"run_id":"r1","ops":[]}'`
  - Structured help export for automation:
    - `gentle_cli help --format json`
    - `gentle_cli help --format markdown`
    - JSON output includes a `capability` object for glossary commands that are
      present in the shared protocol capability registry.

<!-- manual-page-nav:start -->

---

[CLI Manual](../cli.md) | [Previous: File Formats](file-formats.md) | [Next: Introspection, Facts and Readiness](discovery.md)
<!-- manual-page-nav:end -->
