# Public 26R266 Priming Comparisons Without Subtraction

These are **real image-display results**, not synthetic data, raw detector
traces, GUI screenshots or scientific acceptance. The owner approved public
use of the original data on 2026-10-18 and requested this no-subtraction
comparison on 2026-10-22. Use the original
[fixture/provenance](../../../test_files/fixtures/gel_images/tapestation_26r266/README.md)
and its `region_summaries.json` for sample associations. No cell-line names,
gene assignments or recovered primer identities are added.

## Sources And Execution

- After-cleanup PDF SHA-256:
  `859ce71d2ec5f6031c1fcd635adee639ec39a545d306de3a5612d8806526df7f`.
- Retained, unscaled embedded gel PNG SHA-256:
  `9bf36a8c8afa5a33c5b2be0aa6eb51e740e61f700a12abd685c7978117e44c7b`.
- Existing macOS arm64 GENtle CLI source:
  `40aa261d647af83b1ed3b6d23e8d9deac72d67f3`.
- Binary SHA-256:
  `5459d9cd6587c9e8a6ea51dbe04578420c9faf200926b97849ef863ae04eea2e`.
- Checkout at execution: `5ab0b02aef2c05edc15e8434be2cf1429ebf1ce5` plus
  the replay/tutorial changes. The nine sizing/workbench/protocol/renderer
  files listed in `observations.json` matched the binary's canonical Git blobs.
  No new Rust build or performance/profile acceptance is claimed.

The opt-in replay in `scripts/gel_tapestation_acceptance.py` imports the
retained PNGs, reuses saved annotation-assisted calibrations and invokes only
existing GENtle extraction/comparison/export operations. Cell 1 is C1/E1/G1;
cell 2 is D1/F1/H1. C1/D1/E1 use L1, F1/G1/H1 use L2. Both are retained explicit
calibrations, not a claim to recover TapeStation's internal calibration.

Polarity is dark bands, background is explicitly zero and additional masks
are empty. The requested 100-5000 bp window excludes the alignment-marker
centres, not necessarily every marker tail. Unscored edges stay excluded, not
zero. There is no lane subtraction, fitted multiplier, smoothing, normalization
of raw density heights or peak detector in this replay. The engine's CDF/QQ
and distances explicitly normalize each profile's weight on shared valid
support for shape comparison only. No p-value or enrichment estimate is made.

## Reproduce

Use an existing CLI with matching canonical implementation files and a **new**
output directory. This does not compile, fetch, invoke an agent or capture UI:

```sh
python3 scripts/gel_tapestation_acceptance.py \
  --gentle-cli /absolute/path/to/gentle_cli \
  --output-dir /absolute/path/to/new-priming-replay \
  --compare-priming
```

The output contains the retained before/after sizing results, a saved project,
six profile requests, six comparisons and two four-page native vector PDFs
with matching JSON. All exports reload the saved state in separate processes.
`receipt.json` hashes inputs, binary and outputs. No output directory or report
is overwritten. On the recorded host these full native report hashes were:

| File | SHA-256 |
| --- | --- |
| `cell_1-priming.pdf` | `eac573c27f5d2757d2cfd9793b31b578229a48e7f72f360df0bf98cc8c242f71` |
| `cell_2-priming.pdf` | `c4ee6e5105c706b6883c1760977649de78e89c79f77881697c91a59b33c6657b` |
| `cell_1-priming.json` | `87848f49ad5c8497198e97c306cf6cae2177eec4038f0a7615bb150a9476e1e5` |
| `cell_2-priming.json` | `ba8ed76ca10d785039823e232128ac66ec8ae2f3fb028ed3606cb1e2a45dcfc6` |

These are return identities, **not cross-platform golden outputs**. Fonts,
Poppler or host changes can alter rendered bytes. Retain new receipts rather
than changing these historical identities to conceal a difference.

## Retained Teaching Assets

`observations.json` is the replay's `priming-summary.json`, copied without
alteration. It excerpts the native requests, summaries, distances, exclusions
and null p-values; it does not calculate them independently. Full canonical
comparison reports remain reproducible in the replay's per-cell JSON exports.

The four PNGs are whole, unannotated native PDF pages rendered by Poppler
`pdftoppm` 26.06.0, at 96 dpi, without cropping, retouching or contrast changes:

```sh
pdftoppm -f 3 -l 3 -singlefile -r 96 -png cell_1-priming.pdf cell_1_custom_vs_combined
pdftoppm -f 4 -l 4 -singlefile -r 96 -png cell_1-priming.pdf cell_1_manufacturer_vs_combined
pdftoppm -f 3 -l 3 -singlefile -r 96 -png cell_2-priming.pdf cell_2_custom_vs_combined
pdftoppm -f 4 -l 4 -singlefile -r 96 -png cell_2-priming.pdf cell_2_manufacturer_vs_combined
```

`SHA256SUMS` binds these five assets. Recreate the checksum list with
`shasum -a 256 *.png observations.json > SHA256SUMS`. Scoped Git attributes and
the fast LF/CRLF checkout test preserve their bytes. They are used by tutorial
10.03 and `scripts.test_gel_tapestation` for reproducible pairing, return-value
and teaching-content checks; PNGs also provide offline Help illustrations.

## Interpretation Limits

The combined condition is closer to manufacturer-only by both reported
descriptive shape distances in each cell. This is conditional on this extraction,
coverage and display processing; it does not establish additive chemistry or
gene-specific enrichment. Background/marker policy, cross-ladder acceptance,
numerical detector comparisons and native GUI/platform/scientific acceptance
remain open. No human scientific sign-off is implied by retaining the figures.
