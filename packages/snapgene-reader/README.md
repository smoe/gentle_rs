# snapgene-reader

`snapgene-reader` reads SnapGene `.dna` files into reusable Rust data
structures.

The crate is intentionally parser-first:

- no GENtle dependency
- preserves structured SnapGene packets that map cleanly onto stable records
- keeps unknown packet payloads available for downstream inspection
- exposes an optional `gb-io` adapter for projects that already use
  GenBank-like sequence models

## Current packet support

Parsed directly:

- `0x09` cookie
- `0x00` DNA
- `0x05` primers
- `0x06` notes
- `0x08` additional sequence properties
- `0x0A` features

Retained as raw packets in v1:

- any unrecognized packet type, including currently unparsed `0x03` and `0x0D`

## Example

```rust
use snapgene_reader::parse_path;

let parsed = parse_path("vector.dna")?;
println!("{}", parsed.dna.sequence.len());
println!("{:?}", parsed.dna.topology);
# Ok::<(), snapgene_reader::SnapGeneError>(())
```

With the optional `gb-io` feature:

```rust
use gb_io::seq::Seq;
use snapgene_reader::parse_path;

let parsed = parse_path("vector.dna")?;
let seq = Seq::try_from(&parsed)?;
println!("{}", seq.features.len());
# Ok::<(), Box<dyn std::error::Error>>(())
```

## Status

V1 is read-only import support. It does not yet write `.dna` files and it does
not attempt to interpret SnapGene-only display/editor state as biological
behavior.

