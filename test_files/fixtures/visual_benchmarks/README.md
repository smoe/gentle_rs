# Visual Benchmark Fixtures

Small synthetic fixtures for deterministic SVG readability and structure tests.

## Provenance

- `dense_plasmid_map.json` is synthetic and hand-crafted for GENtle visual
  benchmark tests. It is not derived from a biological construct. The sequence
  deliberately places common cloning restriction motifs near annotated vector
  parts to stress dense label placement.
- `antisense_repeat_dotplot_context.json` is synthetic and hand-crafted for
  GENtle visual benchmark tests. It is not derived from a genome assembly. The
  reference span deliberately includes same-strand exons, antisense exons, and
  an rmsk-style `repeat_region` feature so dotplot genome-context rails can be
  tested without external resources.

## Deterministic Recreation

Both fixtures are plain JSON with explicit literal/repeat sequence segments and
0-based half-open feature ranges. To recreate them, expand each
`sequence_segments` entry in order:

- `{"literal": "ACGT"}` appends that literal sequence.
- `{"repeat": "ACGT", "count": 3}` appends `ACGTACGTACGT`.

Then apply each feature range exactly as listed. A one-range feature becomes a
single GENtle/GenBank location; a multi-range feature becomes a joined location;
`strand = "-"` wraps the location in a complement.

## Usage

- `dense_plasmid_map.json` is used by render-export and engine SVG benchmark
  tests for dense feature labels plus common restriction-enzyme labels.
- `antisense_repeat_dotplot_context.json` is used by engine dotplot SVG
  benchmark tests for genome-context side rails with exon, antisense, and
  RepeatMasker-style repeat intervals.
