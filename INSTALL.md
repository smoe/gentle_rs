# Installing and Running GENtle

GENtle is still under development. This guide gets a new contributor
from a clean machine to a first local build and a minimal runtime
check. The examples below focus on macOS because that is where setup
questions tend to come up most often, but the same general flow works
on Linux and Windows as well.

## Install `git` for source code management

GENtle is developed in `git`, so please make sure that command is
available first.

### macOS

Apple ships `git` together with its Xcode command line tools. For many
contributors that is the easiest place to start:

```sh
xcode-select --install
```

If you prefer a Homebrew-managed `git`, that is also perfectly fine.

### Linux

Most Linux systems already have `git` installed or make it available
through the standard package manager.

### Windows

Install `git` through Git for Windows or use the version that comes
with your preferred development environment.

## Install the Rust toolchain

The canonical installation route for Rust is `rustup`, documented at
https://rust-lang.org/learn/get-started/. That is the most portable
choice across macOS, Linux, and Windows.

On macOS or Linux you may also use Homebrew if that better matches your
local setup:

```sh
brew install rust pkg-config
```

Check that the toolchain is available:

```sh
rustc --version
cargo --version
```

GENtle has been compiled successfully with Rust 1.93 and 1.94.

## Optional package manager setup on macOS

If Homebrew is not installed yet, follow the instructions on
https://brew.sh. They currently look like this:

```sh
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```

If you want Homebrew to provide both `git` and Rust, install them with:

```sh
brew install git rust pkg-config
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

Optional scripting shells:

```sh
cargo run --features js-interface --bin gentle_js
cargo run --features lua-interface --bin gentle_lua
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
