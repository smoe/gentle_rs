# Synthetic Detailed TSS Context

All files are hand-crafted synthetic software fixtures, created 2026-09-11.
They contain no private gene data, experimental reads, downloaded sequences,
EPD records or motif weights. They extend the parent fixture's 11-base
`synthetic-minus` window, not a biological promoter.

The canonical source is the committed JSON/FASTA. Recreate `locus.fa` as
`ACG` + reverse-complement(`TTACGACGTAA`) + `TAC`, one 17-base record on a
positive loaded-sequence anchor at synthetic chromosome positions 290..306.
The TSS window is negative-strand, 293..303. The two synthetic exons, CDS
intervals and start/stop *markers* test geometry only; they are not a valid
protein-coding model. Two nonadjacent bedGraph-like intervals demonstrate
clipping and missing signal between supplied intervals. The control is
explicitly `not_prepared`.

`tata.json` illustrates all three distinct evidence kinds: source annotation,
motif prediction and EPD TSS classification. Its numbers and dummy annotation,
anchor and matrix hashes are synthetic, not measurements or validated model
outputs. To recreate its content hash, deserialize as `TataBoxScreenReport`,
clear `content_sha256`, serialize compactly with Rust serde_json, and SHA-256
those bytes. Rust float formatting and struct field order are part of that
existing report contract. Prefix this internal hash with `sha256:`. The context
manifest then binds the exact file bytes of the locus, FASTA and TATA files,
using `shasum -a 256` (unprefixed lowercase digests).

Used by the shared TSS export integration regression in
`src/engine/analysis/tss_profiles_context_tests.rs` and the offline commands in
`docs/tss_tfbs_profiles.md`. Tests verify the real resolver and renderer, not
manually injected projection rows. No source files are modified by the tests.
The positive-strand parent window intentionally has no context entry and keeps
its historical score-only display. The separate inline fixtures in the Rust
test module cover all four combinations of genomic and loaded-sequence strand.
