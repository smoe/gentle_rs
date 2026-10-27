# Build and Run

<!-- manual-counterparts:start -->
Related GUI reference: [Start the GUI](../gui/start-the-gui.md).
These topic links do not imply identical feature coverage across interfaces.
<!-- manual-counterparts:end -->

## Build and run

From the repository root:

```bash
cargo run --bin gentle
cargo run --features js-interface --bin gentle_js
cargo run --features lua-interface --bin gentle_lua
cargo run --bin gentle_cli -- capabilities
cargo run --bin gentle_examples_docs -- --check
cargo run --bin gentle_examples_docs -- tutorial-check
cargo run --bin gentle_mcp -- --help
cargo run --bin gentle -- --version
cargo run --bin gentle_cli -- --version
```

For optimized builds:

```bash
cargo run --release --bin gentle
cargo run --release --features js-interface --bin gentle_js
cargo run --release --features lua-interface --bin gentle_lua
cargo run --release --bin gentle_cli -- capabilities
cargo run --release --bin gentle_examples_docs -- --check
cargo run --release --bin gentle_examples_docs -- tutorial-check
cargo run --release --bin gentle_mcp -- --help
cargo run --release --bin gentle -- --version
cargo run --release --bin gentle_cli -- --version
```

Important: with `cargo run`, arguments for GENtle binaries must come after `--`.
Example: `cargo run --bin gentle -- --version`.

JavaScript and Lua shells are compile-time optional so default `cargo check`
and `cargo build` do not pull in their heavier runtime dependencies. The Python
wrapper in `integrations/python/gentle_py` remains separate from Cargo
features.

Release-build policy:

- default local builds stay lean unless JS/Lua support is requested explicitly
- native release packaging uses Cargo's default release profile and desktop
  features, without embedded JS/Lua; source builds can still opt into
  `--features script-interfaces`
- to build only the five native release entrypoints, use
  `cargo build --locked --release -j1 --bin gentle --bin gentle_cli --bin gentle_mcp --bin gentle_examples_docs --bin gentle_publication_report`

<!-- manual-page-nav:start -->

---

[CLI Manual](../cli.md) | [Previous: Interfaces and Command Discovery](interfaces.md) | [Next: GUI Launcher](gui-launcher.md)
<!-- manual-page-nav:end -->
