# GENtle

A Rust implementation of [GENtle](https://github.com/GENtle-persons/gentle-m).

## Shell Interface

One of the upcoming features of that new implmentation shall be an extension
of the graphical user interface with a shell.

### JavaScript interactive shell
You can run GENtle as an interactive shell using JavaScript.

#### Example
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

### Lua interactive shell
You can run GENtle as an interactive shell using the [Lua programming language](https://www.lua.org/).

#### Example
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

### Vision: Declarative Cloning

The interaction with a cloning too is a bit ambivalent. If you know exactly what you want then you can just look up the sequences in a
database - and frankly this is how many are still doing it. The cloning tool then only comes in afterwards to generate graphics to
explain what you have done to your peers.

If however you do not know what to do, the a graphical representation is of little help since there are too many options for your single
mouse pointer and your single brain to process. And filters alone may not come at your rescue. This especially
holds for more complicated serial cloning events, when one change in the first step (say adjusting the coding frequencies to another species)
may downstream affect the primers or endonucleases to select.

What we are anticipating, without yet an exact plan how to implement it, is to have GENtle offering a series of core operations, and then derive cloning templates from them.
The exact protocol is then constituted by the order of these operations (which could be a template), and what would traditionally be
offered as a manual filter in a dialog box, may be automated as a set of user-driven constraints.
For instance, a core feature could be to derive a subsequence from a given sequence or its inverse complement.
And there could be a function to determine the annealing temperature of any sequence.
Then a PCR can be derived as getting two of these subsequences, but with constraints, roughly like the following:
 * different strands
 * lengths of oligos between 16 and 22
 * annealing temperature within 2 degrees of each other
 * region of interest between primer positions
and depending on the experimental setup it may be required to also demand
 * Genomic uniqueness of primers / primer pairing

Essentially, the idea is to express in writing (and likely with some support by the GUI) what you want to be cloned and then GENtle should figure out how it is to be done. 
The above mentioned PCR scenario is implemented in src/constraints.rs. For an early look at how JSON-formatted constraint could be used to direct the cloning, and to think along, please run

```bash
cargo test constraints
```



## Installation

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


