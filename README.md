A Rust implementation of [GENtle](https://github.com/GENtle-persons/gentle-m).

## Install
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
CFLAGS="-march=armv8-a" cargo run --release
```
Note: Currently only loads two test sequences, hardcoded.
