# Glen Audit: Native Linux GUI At Current `.11` Main

Prepared: 2026-09-27. This report separates native correctness observations
from interaction timing. It is an external Linux diagnostic, not package,
macOS, biological or `.12` release acceptance.

## Candidate And Environment

- Candidate source: `0467867cce153ab4c89c4093c32c003b9f70e513`, the
  fetched `origin/main` selected before the audit.
- Package workflow: Actions run
  [`36273325748`](https://github.com/smoe/gentle_rs/actions/runs/36273325748)
  built the exact-SHA Linux and Windows packages; macOS failed and collection
  and publication were skipped. The Linux artifact exists, but Glen's local
  `gh` credential could not download it. No extracted-package claim follows.
- Tested binary: local locked/offline rebuild of the current internal-installer
  `dev` recipe, `-j1`, `CARGO_INCREMENTAL=0`,
  `CARGO_PROFILE_DEV_DEBUG=0`, default features. It reports
  `0.1.0-internal.11+git.0467867cce153ab4c89c4093c32c003b9f70e513`.
- GUI SHA-256:
  `d5ab09b20ee27c5e461b19d2d457d9631a90d516e7b92247dfe771f8e2e4a525`.
  `Cargo.lock` SHA-256:
  `576ee59fbadd5075e12570be62365e39245a9be94b3e324276a2527108653c6c`.
- Toolchain: `rustc 1.95.0 (59807616e 2026-04-14)`,
  `x86_64-unknown-linux-gnu`.
- Display: Xvfb `1600x1000x24`, fresh HOME/XDG/runtime/temp directories for
  every subject. No EWMH window manager was installed, so window geometry and
  screenshots are useful diagnostics but do not constitute the full packaged
  Xvfb/Openbox contract.
- Inputs were public repository fixtures only. No private project, network
  query, external agent model or database query was used.

The complete from-empty-cache product build was first run while the worktree
contained only Glen's documentation commit `b89b2135`; that commit changes no
product source. It completed in 16:56.16 with 8,737,820 KiB maximum RSS, exit
zero and no swap. Identity verification correctly rejected that binary as the
candidate because it embedded `b89b2135`. Rebinding the root crate and five
binaries at detached `0467867c` completed in 6:38.85 with 8,859,824 KiB maximum
RSS, exit zero and no swap. Only the latter exact-identity binary was exercised.
This two-step local provenance is not equivalent to extracting the Actions
artifact.

## Correctness Observations

### Authentic PATZ1

The pinned preparation helper produced the same public 20,802-bp minus-strand
PATZ1 project and 13 Ensembl transcript rows used by tutorial 04.08. The project
SHA-256 is
`51e065e695e5e5b345a432f680b69bae2d768dee761f6c0448565b63549c0121`;
the locus report retained the established hash
`4f84705b05432cda7699d6fbe027f4c8fa99775b9cac904d241fd0583e0cac27`.

Ordinary X11 input and raster captures confirmed:

- the verified GRCh38 minus-strand anchor and 20,802-bp range;
- thirteen mRNA rows, including transcript IDs and product labels;
- transcript selection and its exact exon/product detail;
- mRNA layer off/on changing the map without losing the selected project;
- rendered content after `820x520`, `1200x800`, `1600x1000` and `1200x860`
  resizes;
- zoom from the full locus to `9752..11051` and explicit range navigation to
  `12000..13000`.

No coordinate or transcript-label discrepancy was observed. The observed
window geometry changed in 7--8 ms, but that measures X11 geometry acknowledgement,
not completed content paint. Selection, layer, zoom and navigation captures
used fixed 250--500 ms observation waits and therefore provide upper bounds,
not render-time estimates.

### TSS Regulatory Viewer

The public 11-base synthetic minus-strand annotated TSS fixture and bound
profile hashes were respectively
`8acf0b32a52d8f6f16cbdce57cefb52372fab9cee57e6e276b2fc63e8c197692`
and
`4478630c0ec4dded1605618e4d37bc45a4ac38547ec60aaf1827aa4769ca134f`.

The native viewer showed the annotated minus-strand TSS at genomic coordinate
300/local base 4, unavailable control as unavailable rather than zero, and
separate structure, signal and stored-motif lanes. Agent Assistant recognized
the literal reviewed command
`ui open tss-view --report .../synthetic-minus.profile-report.json` as a direct
shared command. Its result named the bound subject and stated that reference,
TSS geometry and sequence hash would be validated and that no scoring or
database query was started. The viewer then identified the exact producer and
displayed the report-provided `MA0004.1` curve lane.

This teaching profile deliberately contains no imported DuckDB hit rows. Its
`DuckDB peaks` control and the no-query/no-hit state were visible, but live
rendering of nonempty imported hits is **not accepted by this run**. A retained
profile with nonempty, hash-bound imported evidence is still required for that
checkpoint.

### Circular Map Regression

The public circular `pGEX-3X` GenBank fixture retained its 4,952-bp circular
topology. The shared `RenderSequenceSvg` operation produced a circular SVG with
one circle and fourteen feature paths, SHA-256
`659ae23c58db974bb4ba7faf94f3b4f6e31a8f9a13cc3b0ed8433d3ab31a396a`.
This is a correctness regression check only; no circular timing target or
native circular-interaction acceptance is claimed.

## Startup And First-Content Attribution

Three fresh-profile runs reproduced the same dominant startup span:

| Subject | Help preparation | Root workspace checkpoint |
| --- | ---: | ---: |
| Empty project | 37.461 s | 37.909 s |
| PATZ1 project | 37.807 s | 38.360 s |
| Synthetic TSS project | 37.698 s | 38.189 s |

`help_preparation` is therefore the next measured product-side startup slice on
this unoptimized package recipe. It is not project decode, credential refresh,
the feature-tree rebuild or compositor presentation. PATZ1 project read/decode
was 13.053 ms and install 0.531 ms; the TSS values were 2.474 ms and 0.503 ms.

For PATZ1, DNA open dispatch took 272.001 ms, deferred construction 271.915 ms,
hydration 295.271 ms and dispatch-to-first-native-content 833.346 ms. For the
tiny TSS fixture the corresponding values were 5.594 ms, 5.526 ms, 5.287 ms and
217.016 ms. These are CPU checkpoints; they do not establish pixels presented
to the user. The independently detected PATZ1 child window appeared about
410 ms after the double-click, before its content was confirmed.

## Verdict And Next Step

The native Linux run accepts the observed PATZ1 and TSS correctness paths at
the selected source revision, with the stated Xvfb/no-window-manager limit. It
does not accept nonempty DuckDB-hit rendering, macOS behavior or the exact
downloaded package. Correctness and speed remain separate: the displayed
states were sound, while startup was dominated by a repeatable 37.5--37.8 s
help-preparation span.

The next implementation slice should instrument and bound the work inside
`help_preparation`, then remove or defer its largest proven component. Do not
select another feature-rendering optimization from this evidence. Re-run the
same fresh-profile traces on the repaired build and on macOS. Separately,
download and smoke the exact Linux Actions artifact once authenticated, and
add one public nonempty imported-motif profile before claiming DuckDB GUI
acceptance.
