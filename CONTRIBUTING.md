# Contributing to GENtle

Thank you for your interest in contributing to GENtle.

This guide is intentionally short and practical. It focuses on getting a new
contributor from a clean macOS machine to a first local build and a minimal
validation loop before opening a pull request.

## macOS Setup

The recommended way to install local build dependencies on macOS is via
Homebrew.

### 1. Install Apple's command line tools

```sh
xcode-select --install
```

### 2. Install Homebrew

If Homebrew is not installed yet:

```sh
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```

### 3. Install Rust and common build helpers

```sh
brew install rust pkg-config
```

Check that the toolchain is available:

```sh
rustc --version
cargo --version
```

## Clone and Build

Clone the repository and run a first local build check:

```sh
git clone https://github.com/smoe/gentle_rs.git
cd gentle_rs
cargo check -q
```

## Run GENtle

Start the CLI:

```sh
cargo run --bin gentle_cli -- capabilities
```

Start the GUI:

```sh
cargo run --bin gentle
```

Important: when using `cargo run`, arguments for GENtle binaries must come
after `--`.

Example:

```sh
cargo run --bin gentle_cli -- --version
```

## Optional External Tools

Some GENtle features can use external helper applications when available, such
as:

- `blastn`
- `makeblastdb`
- `bigWigToBedGraph`
- `rnapkin`

These tools are optional for a first build. They can be configured in the GUI
under the external applications settings.

## Minimal Validation Loop

Before opening a pull request, a useful baseline is:

```sh
cargo check -q
cargo test -q workflow_examples -- --test-threads=1
cargo run --bin gentle_examples_docs -- --check
cargo run --bin gentle_examples_docs -- tutorial-check
```

## Project Documentation

Depending on the area you want to change, these documents are the most useful
entry points:

- `docs/architecture.md`
- `docs/roadmap.md`
- `docs/protocol.md`
- `docs/gui.md`
- `docs/cli.md`
- `docs/testing.md`

## Contribution Notes

- Keep behavior aligned across GUI, CLI, and other interfaces when possible.
- Update documentation when user-visible behavior changes.
- Prefer deterministic checks and reproducible examples when validating a
  change.
