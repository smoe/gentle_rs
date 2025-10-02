# GENtle

A modern re-implementation of [GENtle](https://github.com/GENtle-persons/gentle-m)
 * in Rust,
 * with an agentic interface
 * and cloning-steps represented as a directed graph for perfect
   - provenance,
   - documentation and 
   - replication.

## Documentation

- GUI manual: [`docs/gui.md`](docs/gui.md)
- CLI manual: [`docs/cli.md`](docs/cli.md)
- Architecture guide: [`docs/architecture.md`](docs/architecture.md)
- Engine protocol draft: [`docs/protocol.md`](docs/protocol.md)
- Testing strategy: [`docs/testing.md`](docs/testing.md)

## GUI toolbar buttons

![GENtle Sequence Display](assets/gentle_screenshot.png)

In each DNA window, the top toolbar contains icon buttons that control the map
and sequence views. Hovering a button in the app shows a tooltip with the same
description.

From left to right:

1. Circular/Linear map toggle
   - Switches DNA topology visualization between circular and linear map mode.
2. Show/Hide sequence panel
   - Toggles the sequence text panel.
3. Show/Hide map panel
   - Toggles the graphical DNA map panel.
4. Show/Hide annotated features
   - Toggles feature rendering (e.g., genes/CDS and other annotations).
5. Show/Hide restriction enzyme sites
   - Toggles displayed restriction cut-site markers and labels.
6. Show/Hide GC content
   - Toggles GC-content visualization on the map.
7. Show/Hide open reading frames (ORFs)
   - Toggles predicted ORF overlays.
8. Show/Hide methylation sites
   - Toggles methylation-site markers on the map.

## JavaScript interactive shell
You can run GENtle as an interactive shell using JavaScript.

### Example
```
> cargo run --release --bin gentle_js
Interactive JavaScript Shell (type 'exit' to quit)
GENtle> pgex = load_dna("test_files/pGEX-3X.gb");
GENtle> console.log(pgex.seq.seq.length)
4952
GENtle> results = digest(pgex,"BamHI,EcoRI"); // Digest pGex-3X with BamHI and EcoRI
GENtle> console.log(results.length) // Number of sequences
2
GENtle> console.log(results[0].seq.seq.length) // Length of first sequence
6
GENtle> console.log(results[1].seq.seq.length) // Length of second sequence
4938
```

## Lua interactive shell
You can run GENtle as an interactive shell using the [Lua programming language](https://www.lua.org/).

### Example
```
> cargo run --release --bin gentle_lua
Interactive Lua Shell (type 'exit' to quit)
(...)
GENtle> pgex = load_dna("test_files/pGEX-3X.gb") -- loads a GenBank sequence and performs some computations on it
GENtle> #pgex.seq.seq -- prints the length of the sequence
4952
GENtle> pgex.restriction_enzyme_sites -- shows the precomputed restriction enzyme sites
(...)
GENtle> pgex.methylation_sites -- shows the precomputed methylation_sites
(...)
GENtle> write_gb("output.gb",pgex) -- writes the sequence to a new GenBank file
```

## Install

- [Install Rust](https://www.rust-lang.org/tools/install)
```bash
git clone https://github.com/magnusmanske/gentle_rs/
cd gentle_rs
cargo run --release --bin gentle
```
Note: Cargo will compile >360 dependencies.
On non-Intel platforms it may be beneficial to explicitly specify
your build architecture to avoid an error in some packages.
E.g., for the MacBook M1 with conda, run
```[bash]
CFLAGS="-march=armv8-a" cargo run --release --bin gentle
```
Note: Currently only loads two test sequences, hardcoded.

### Bundling GENtle for MacOS

To present all files to the operating system, including the application's icon,
that are required for GENtle to run, all application data needs to be "bundled"
as follows:

```[bash]
cargo install cargo-bundle
CFLAGS="-march=armv8-a" cargo bundle --release
```
To install the .app just copy it to an "Application" folder:
```[bash]
cp -r target/release/bundle/osx/gentle.app ~/Applications/
```
