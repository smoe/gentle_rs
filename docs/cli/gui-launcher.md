# GUI Launcher

<!-- manual-counterparts:start -->
Related GUI reference: [Start the GUI](../gui/start-the-gui.md).
These topic links do not imply identical feature coverage across interfaces.
<!-- manual-counterparts:end -->

## `gentle` (GUI launcher)

`gentle` starts the graphical application.

```bash
cargo run --bin gentle
cargo run --bin gentle -- path/to/project.gentle.json
cargo run --bin gentle -- --version
```

If a project path is provided, GENtle opens that project at startup.

Current behavior:

- opens GUI windows
- opens an empty project by default
- optionally opens a project path passed on startup
- supports project/sequence actions from the File menu

<!-- manual-page-nav:start -->

---

[CLI Manual](../cli.md) | [Previous: Build and Run](build.md) | [Next: CLI Options and Progress](cli-adapter.md)
<!-- manual-page-nav:end -->
