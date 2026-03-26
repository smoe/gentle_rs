# Installing and Running GENtle

GENtle is still under development. This guide is meant to get a new
contributor from a clean machine to:

- a first local build,
- a first successful CLI invocation,
- a first GUI launch.

The examples below focus on practical contributor setup. macOS tends to
raise the most setup questions, but the same general flow works on Linux
and Windows as well.

## Prerequisites

Before you start, please make sure you have:

- about `10 GB` of free disk space in practice
  - Rust toolchains, build artifacts, cloned repositories, and optional
    genome/helper resources add up quickly
- a reasonable internet connection
  - needed for cloning the repository, downloading Rust crates, and, when
    used, downloading genome/helper resources
- a terminal or shell you are comfortable using
- permission to install developer tools on your machine

For a first local build you only need:

- `git`
- a Rust toolchain (`rustc`, `cargo`)

The external biology helper tools described later are optional.

## Install `git`

GENtle is developed in `git`, so please make sure that command is
available first.

### macOS

Apple ships `git` together with its Xcode command line tools. For many
contributors that is the easiest place to start:

```sh
xcode-select --install
```

If you prefer a Homebrew-managed `git`, that is also perfectly fine:

```sh
brew install git
```

### Linux

Most Linux systems either already have `git` installed or make it
available through the standard package manager.

Examples:

```sh
sudo apt install git
sudo dnf install git
sudo pacman -S git
```

Use whichever command matches your distribution.

### Windows

Install `git` through Git for Windows, or use the version that comes
with your preferred development environment.

If you already use WSL for development, install `git` inside that WSL
environment and run GENtle from there.

## Install the Rust toolchain

The canonical installation route for Rust is `rustup`, documented at:

- <https://rust-lang.org/learn/get-started/>

That is the most portable choice across macOS, Linux, and Windows.

After installation, check that the toolchain is available:

```sh
rustc --version
cargo --version
```

GENtle has been compiled successfully with Rust `1.93` and `1.94`.

### Alternative on macOS or Linux: Homebrew

If you prefer to use Homebrew for Rust on macOS or Linux:

```sh
brew install rust pkg-config
```

`pkg-config` is not always strictly required for every first build, but
it is a useful piece of general build infrastructure to have around.

## Optional package manager setup on macOS

If Homebrew is not installed yet, follow the instructions on
<https://brew.sh>. They currently look like this:

```sh
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```

If you want Homebrew to provide both `git` and Rust, install them with:

```sh
brew install git rust pkg-config
```

## Clone the repository

Clone the repository and move into the project directory:

```sh
git clone https://github.com/smoe/gentle_rs.git
cd gentle_rs
```

## First local build check

Run a first local build check:

```sh
cargo check
```

If this succeeds, your core contributor toolchain is in place.

## Alternative: run GENtle through Docker

If you do not want to prepare a local Rust toolchain first, GENtle can also be
run through the Debian-first container image.

This route is currently aimed at macOS and Linux hosts.

### Build the image locally

From the repository root:

```sh
docker build -t gentle:local .
```

The image includes:

- the GUI
- `gentle_cli`
- the MCP server
- the embedded JavaScript and Lua shells
- the Python wrapper path
- external helper-tool support already wired into GENtle

### Start the GUI in a browser

```sh
docker run --rm -it \
  -p 6080:6080 \
  -v "$(pwd)":/work \
  gentle:local
```

Then open:

- <http://localhost:6080/vnc.html?autoconnect=1&resize=scale>

### Run the CLI from the image

```sh
docker run --rm -it \
  -v "$(pwd)":/work \
  gentle:local cli capabilities
```

### Linux users who prefer Apptainer / Singularity

GENtle currently keeps the Dockerfile as the single maintained image
definition. Linux users can consume the same OCI image from Apptainer instead
of maintaining a separate `.def` recipe.

See:

- `docs/container.md`

## Run GENtle

### Start the graphical user interface (GUI)

This is the interface most users will want to use.

```sh
cargo run --bin gentle
```

The GUI opens an empty project unless you pass an existing project file.

### Access from the command line or shell programming languages

GENtle is designed so GUI actions correspond to shared engine operations,
and much of the functionality is already available from the command line.
In short: many "OK" actions in GENtle windows correspond to an
"operation" that can also be started without clicking through the same
windows again. This comes in handy when only changing a small parameter
(like the mutation induced), when redoing the same operations on another
genomic background, or when exchanging cloning concepts with an AI.

```sh
cargo run --bin gentle_cli -- capabilities
```

GENtle can also be programmed not from the command line but from
regular programming languages. Besides JavaScript and Lua, also an
interface for Python is offered.

```sh
cargo run --features js-interface --bin gentle_js
cargo run --features lua-interface --bin gentle_lua
```

The Python wrapper is documented separately in:

- `integrations/python/README.md`

and can be installed locally with:

```sh
python3 -m pip install -e integrations/python
```

Important: when using `cargo run`, arguments for GENtle binaries must
come after `--`.

Example:

```sh
cargo run --bin gentle_cli -- --version
```

Container note:

- the Docker image provides the same shared-entrypoint routes through
  `gui-web`, `cli`, `mcp`, `js`, and `lua`
- full details live in `docs/container.md`

## First successful run checklist

You are in a good state when all of the following work:

1. `git --version`
2. `rustc --version`
3. `cargo --version`
4. `cargo check -q`
5. `cargo run --bin gentle_cli -- capabilities`
6. `cargo run --bin gentle`

## Optional external helper tools

Some GENtle features can use external helper applications when
available, such as:

- `blastn`
- `makeblastdb`
- `bigWigToBedGraph`
- `rnapkin`

These tools are optional for a first build and for a first GUI launch.

They become useful for more advanced workflows:

- `blastn` and `makeblastdb`
  - local BLAST-backed genome/helper searches and index preparation
- `bigWigToBedGraph`
  - importing `BigWig` genome tracks through conversion
- `rnapkin`
  - RNA secondary-structure inspection/export
- `primer3`
  - primer design backend for primer-pair and qPCR workflows


Common ways to obtain them:

- package manager on your platform
- Homebrew
- Conda / Mamba environment
- manual installation into a directory on `PATH`

The following tools can be configured globally in:

- `Settings -> Configuration... -> External Applications`

This applies to:

- `blastn`
- `makeblastdb`
- `bigWigToBedGraph`
- `rnapkin`

`primer3` is configured through the primer-design workflow itself
(for example in `Engine Ops -> Primer and qPCR design reports`) or via
CLI/shared-shell options such as `--backend primer3 --primer3-exec PATH`.

Environment-variable overrides are also supported:

- `GENTLE_RNAPKIN_BIN`
- `GENTLE_MAKEBLASTDB_BIN`
- `GENTLE_BLASTN_BIN`
- `GENTLE_BIGWIG_TO_BEDGRAPH_BIN`

### Notes

- Add the optional biology helper tools only when you actually start
  using the workflows that require them.
- If you are on Windows and run into avoidable tool-availability issues,
  WSL can be a practical development environment for GENtle.
- If you prefer containerized execution on macOS or Linux, use
  `docs/container.md` as the primary runtime guide.
