A Rust implementation of [GENtle](https://github.com/GENtle-persons/gentle-m).

# Install
- [Install Rust](https://www.rust-lang.org/tools/install)
```bash
git clone https://github.com/magnusmanske/gentle_rs/
cd gentle_rs
cargo run --release
```
Note: Cargo will compile >360 dependencies.
On non-Intel platforms it may be beneficial to explicitly specify
your build architecture to avoid an error in some packages.
E.g., for the MacBook M1 with conda, run
```[bash]
CFLAGS="-march=armv8-a" cargo run --release --binary gentle
```
Note: Currently only loads two test sequences, hardcoded.

# Lua interactive interface
You can run GENtle as an interactive interface using the [Lua programming language](https://www.lua.org/).



## Example
```
> cargo run --release --binary gentle_cli
Interactive Lua Shell (type 'exit' to quit)
(...)
> pgex = load_dna("test_files/pGEX-3X.gb") -- loads a GenBank sequence and performs some computations on it
> #pgex.seq.seq -- prints the length of the sequence
4952
> pgex.restriction_enzyme_sites -- shows the precomputed restriction enzyme sites
(...)
> pgex.methylation_sites -- shows the precomputed methylation_sites
(...)
> write_gb("output.gb",pgex) -- writes the sequence to a new GenBank file
```
