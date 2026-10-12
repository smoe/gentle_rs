# RNAPKIN Container Dependency Lock

The container still builds upstream **rnapkin 0.3.9**, without source changes,
optimization or disabled development assertions. Its published lockfile selects
`plotters-bitmap 0.3.2`, which performs aligned `u64` reads/writes on byte-aligned
RGB storage. Container run `36561534544`, job `109383399310`, at GENtle `8a87118f`
aborted in `rgb.rs:215` while rendering PNG with the unoptimized helper. SVG
rendering and GENtle compilation had succeeded.

## Sources And Scope

- Published source: <https://static.crates.io/crates/rnapkin/rnapkin-0.3.9.crate>
- Archive SHA-256, verified against the crates.io registry index:
  `4495690197e1cced9b16234d6a66b40ebf190e9613b9c1c6aea837adeed00f17`.
- The adjacent `Cargo.lock` is Cargo-generated from that archive's lockfile.
  Only three package versions/checksums differ: `plotters-bitmap 0.3.2 -> 0.3.3`,
  `plotters-backend 0.3.4 -> 0.3.7`, and `gif 0.11.4 -> 0.12.0`. Cargo also
  upgrades the lockfile format from version 3 to 4. Other resolutions are intact;
  GENtle's root lockfile is not involved.
- Plotters-bitmap 0.3.3 includes the upstream alignment corrections:
  [rectangle writes](https://github.com/plotters-rs/plotters/commit/9abcb2f6dfb3de3df225b815cdc4edda3098ae99),
  [blend reads](https://github.com/plotters-rs/plotters/commit/a0239c32a5caaa17d627b42918ed24ece2c55b18),
  and [blend writes](https://github.com/plotters-rs/plotters/commit/94e22a3605195b6f74ad4f14efed1579661679e3).
  The backend and GIF versions satisfy its newer dependency requirements.

## Recreate Or Update

In a disposable directory, download the archive, verify the SHA-256 above, and
extract it. Starting from its published lockfile:

```sh
cargo update --manifest-path rnapkin-0.3.9/Cargo.toml -p plotters-bitmap --precise 0.3.3
cargo update --manifest-path rnapkin-0.3.9/Cargo.toml -p plotters-backend --precise 0.3.7
cargo update --manifest-path rnapkin-0.3.9/Cargo.toml -p gif --precise 0.12.0
```

Review the diff against the archived lockfile before replacing this one. Do not
run an unconstrained update, remove `--locked`, or change GENtle's build profile
to work around a helper failure. Docker verifies the source archive before
extraction, substitutes this lockfile, and installs with `--locked --debug -j1`.

## Verification Boundary

The Dockerfile renders the same hand-crafted nine-base hairpin to SVG and PNG
before compiling GENtle. Container CI repeats both renders without networking
as the final image's unprivileged user and checks dynamic linking and fonts.
`scripts.test_container` guards the pinned build recipe and lock resolutions;
`scripts.test_tutorial_checkouts` guards lockfile bytes under LF/CRLF checkouts.
The container receipt records this helper lockfile's digest separately.

Dependency resolution and source inspection are not execution acceptance. No
local build or test was run for this repair; CI must confirm the corrected
helper and final image. Rendering this synthetic hairpin is a packaging check,
not validation of RNA folding or scientific output equivalence.
