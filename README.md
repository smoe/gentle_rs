A Rust implementation of [GENtle](https://github.com/GENtle-persons/gentle-m).

# Install
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

# JavaScript interactive shell
You can run GENtle as an interactive shell using JavaScript.

## Example
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

# Lua interactive shell
You can run GENtle as an interactive shell using the [Lua programming language](https://www.lua.org/).

## Example
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
