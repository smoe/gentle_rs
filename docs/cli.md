# GENtle CLI Manual

This page documents command-line entry points for GENtle.

## Overview

GENtle currently provides three binaries:

- `gentle`: graphical desktop app
- `gentle_js`: interactive JavaScript shell
- `gentle_lua`: interactive Lua shell
- `gentle_cli`: JSON operation/workflow CLI for automation and AI tools

## Build and run

From the repository root:

```bash
cargo run --bin gentle
cargo run --bin gentle_js
cargo run --bin gentle_lua
cargo run --bin gentle_cli -- capabilities
```

For optimized builds:

```bash
cargo run --release --bin gentle
cargo run --release --bin gentle_js
cargo run --release --bin gentle_lua
cargo run --release --bin gentle_cli -- capabilities
```

## `gentle` (GUI launcher)

`gentle` starts the graphical application.

```bash
cargo run --bin gentle
```

Current behavior:

- opens GUI windows
- loads default test sequences at startup
- supports opening additional files from `File -> Open`

## `gentle_js` (JavaScript shell)

`gentle_js` starts an interactive JavaScript REPL backed by GENtle data structures.

```bash
cargo run --bin gentle_js
```

Exit methods:

- type `exit`
- `Ctrl+C`
- `Ctrl+D`

### JavaScript shell functions

1. `load_dna(path)`
   - Loads a DNA sequence from file.
2. `write_gb(seq, path)`
   - Writes a sequence object to a GenBank file.
3. `digest(seq, "EnzymeA,EnzymeB,...")`
   - Performs full digest with comma-separated restriction enzyme names.

### JavaScript example

```javascript
pgex = load_dna("test_files/pGEX-3X.gb");
console.log(pgex.seq.seq.length);
parts = digest(pgex, "BamHI,EcoRI");
console.log(parts.length);
write_gb(parts[0], "output_part0.gb");
```

## `gentle_lua` (Lua shell)

`gentle_lua` starts an interactive Lua REPL backed by GENtle data structures.

```bash
cargo run --bin gentle_lua
```

Exit methods:

- type `exit`
- `Ctrl+C`
- `Ctrl+D`

### Lua shell functions

1. `load_dna(filename)`
   - Loads a DNA sequence from file.
2. `write_gb(filename, seq)`
   - Writes a sequence object to a GenBank file.

Note: In the current implementation, `digest` is available in JavaScript shell but not yet exposed in Lua shell.

### Lua example

```lua
pgex = load_dna("test_files/pGEX-3X.gb")
print(#pgex.seq.seq)
write_gb("output.gb", pgex)
```

## File format expectations

Current CLI workflows rely on sequence files supported by internal loaders:

- GenBank
- FASTA

## `gentle_cli` (automation/agent CLI)

`gentle_cli` exposes operation-based control with JSON input and output.
It is intended for scripted workflows and AI agents.

State is persisted to `.gentle_state.json` by default (override with `--state PATH`).

### Commands

```bash
cargo run --bin gentle_cli -- capabilities
cargo run --bin gentle_cli -- state-summary
cargo run --bin gentle_cli -- op '<operation-json>'
cargo run --bin gentle_cli -- workflow '<workflow-json>'
cargo run --bin gentle_cli -- export-state state.json
cargo run --bin gentle_cli -- import-state state.json
cargo run --bin gentle_cli -- render-svg pgex linear pgex.linear.svg
cargo run --bin gentle_cli -- render-svg pgex circular pgex.circular.svg
```

You can pass JSON from a file with `@file.json`.

### Example operations

Load a file:

```json
{"LoadFile":{"path":"test_files/pGEX-3X.gb","as_id":"pgex"}}
```

Digest:

```json
{"Digest":{"input":"pgex","enzymes":["BamHI","EcoRI"],"output_prefix":"pgex_frag"}}
```

Ligation:

```json
{"Ligation":{"inputs":["pgex_frag_1","pgex_frag_2"],"circularize_if_possible":true,"output_id":"re_ligated"}}
```

PCR:

```json
{"Pcr":{"template":"pgex","forward_primer":"ATGGCT","reverse_primer":"CGTACC","output_id":"amplicon1"}}
```

Extract region (`from` inclusive, `to` exclusive):

```json
{"ExtractRegion":{"input":"pgex","from":100,"to":900,"output_id":"insert_candidate"}}
```

Set visibility of a GUI-equivalent display toggle (example: features):

```json
{"SetDisplayVisibility":{"target":"Features","visible":false}}
```

Available `target` values:

- `SequencePanel`
- `MapPanel`
- `Features`
- `RestrictionEnzymes`
- `GcContents`
- `OpenReadingFrames`
- `MethylationSites`

Save as GenBank:

```json
{"SaveFile":{"seq_id":"pgex_frag_1","path":"frag1.gb","format":"GenBank"}}
```

### Current limitations in the new operation layer

- Ligation currently concatenates input sequences in order without sticky-end compatibility checks.
- PCR currently supports linear templates and exact primer matching only.

## Error behavior

- REPLs print runtime errors without exiting automatically.
- Invalid file paths or unsupported content produce load/write errors.
- Some functions assume computed annotations are available after `load_dna`.
