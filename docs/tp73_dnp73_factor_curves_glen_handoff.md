# Handoff: test, benchmark and illustrate the TP73 ΔNp73 factor-curve tutorial

Addressed to Glen. Scope: acceptance-test and benchmark the new native TSS
local-scoring and export routes, then add screenshots to the new tutorial.

- Tutorial: [`docs/tutorial/08-17_tp73_dnp73_factor_curves.md`](tutorial/08-17_tp73_dnp73_factor_curves.md)
- Reviewed revision: record the exact `git rev-parse HEAD` you tested. Do not
  assume this document's revision is still current — `main` moves.
- Input: `test_files/tp73.ncbi.gb` only. No network, no prepared genome, no
  private data. If a step seems to need any of those, stop and report it.

## What is new and therefore unproven

Two routes were added so the tutorial is executable by the inner/outer agent,
not only by hand in the GUI:

1. `Operation::ExportTssViewSvg`, shell `promoters tss-view-svg SEQ_ID OUT.svg
   [--report R] [--motif ACC]... [--score-kind K] [--keep-negative]
   [--span A..B] [--width PX]`. Fully headless: it decodes the annotated window,
   optionally attaches a profile report, optionally computes local curves, and
   writes the same SVG the GUI's **Export View SVG** produces.
2. `ui open|focus tss-view --local-score ACC[,ACC...] [--score-kind K]
   [--keep-negative]`, which mirrors the existing `--report` intent and is
   applied by the GUI host to the active viewer.

The implementation author verified the headless chain end-to-end on macOS and
added four regression tests for export, shell parsing, hosted settings adoption
and budget refusal.
What I did **not** do: run the GUI, capture any screenshot, measure performance,
or test on Windows. That is the gap this handoff asks you to close.

## 1. Reproduce the tutorial exactly as written

Work through [08-17](tutorial/08-17_tp73_dnp73_factor_curves.md) three times —
once via the GUI, once via Agent Assistant, once by asking the inner agent — and
confirm all three reach the same lanes. Please report each divergence verbatim
rather than smoothing it over in prose.

Expected values on the retained record, which you can check without the GUI:

| Check | Expected |
| --- | --- |
| Recovered anchor | `1:3652516-3736201 (GRCh38.p14, strand +)` |
| Annotated starts found | 5 |
| P1 / TAp73 start `1:3,652,516` | `missing_flanks` |
| ΔNp73 start `1:3,690,672` | `available`, 6 transcripts |
| Lanes in the exported view | 11 (3 local curves, 6 transcript models, gene, TSS marker) |
| Evaluated strand-windows | TP73 1372/1372, E2F1 1380/1380, PATZ1 1382/1382 |
| Coordinate alignment | `local 501 \| TSS +0 bp \| 1:3690672 (+)` |
| Amber unavailable bands | none in this window |

The three window counts differ only because the matrices are 16, 12 and 11 bp
long. If they come out equal, something is wrong.

## 2. Acceptance checks I could not run

Please confirm or refute each, with the revision and platform recorded:

1. **GUI/headless equivalence.** Produce the view once through **Export View
   SVG** and once through `promoters tss-view-svg` with the same motifs, score
   kind and span. The lane set and coordinates must agree. The SVG bytes may
   differ; say so explicitly if they do, and whether anything scientific differs.
2. **Hosted intent really applies.** From a focused TP73 window, run
   `ui open tss-view --local-score MA0861.2,MA0024.3,MA1961.2 --score-kind
   llr_background_tail_log10` in Agent Assistant. The viewer's own **Local
   scoring** controls should end up showing those three accessions and that score
   kind, and the computed lanes should match them. This adoption step is new.
3. **Independence.** Attach a profile report and compute local curves in the same
   viewer. Then detach the report; the local curves must survive. Then clear
   local scores; the report lanes must survive. Confirm the lanes never share a
   scale.
4. **Cancellation.** Start local scoring and cancel it. The previous complete
   result must remain, and no partial curve may appear.
5. **Budget refusals are honest.** Try 33 accessions, and try a window larger
   than 50,000 bp. Both must be refused with a message naming the limit, and
   must point at the headless route rather than silently truncating.
6. **Alias refusal.** `--motif TP73` (a factor name, not an accession) must fail.
   Confirm the error says why.
7. **Windows.** Run at least the headless export and the four new regression tests on
   native Windows. The export writes atomically through a temporary file in the
   target directory; that path is untested there. Per our cross-platform rules a
   macOS pass is not evidence for Windows.

## 3. Benchmark

The tutorial's window is 701 bp, which says nothing about behaviour at realistic
sizes. Two numbers matter, and they are different numbers:

**Scoring cost** — scales with window length × matrices × 2 strands. Measure
`promoters tss-view-svg` wall time with 1, 3, 8 and 16 matrices on windows of
roughly 700, 2,000, 10,000 and 50,000 bp. The 50,000 bp / 16 matrix corner is at
the documented admission ceiling, so it is the honest worst case for the route.
Note that a fresh run re-resolves the registry and re-derives the background
calibration per matrix; do not report a warm second run as the cost of the first.

**Render and paint cost** — the current exporter emits one hover title per
scored window start. The 701 bp / 3 matrix export is already **1.2 MB** of SVG
with 4,134 such titles. Please measure how file size and render time grow at
10,000 and 50,000 bp, and check whether a viewer or a downstream converter
becomes unusable before the scoring does. If SVG size turns out to be the real
ceiling rather than scoring, that is a finding worth its own roadmap entry — I
have deliberately not changed that behaviour.

For the GUI, please also measure interactive pan/zoom on a 10,000 bp window with
8 local curves attached, using the actual packaged profile and recipe: currently
unoptimized `dev` with `incremental=false` and `debug=0`. Do not substitute an
optimized `release` or `bench-audit` build. Bind the result to the exact revision
and effective build settings as in our other GUI latency audits.

## 4. Screenshots for the tutorial

Once you are satisfied the workflow is correct, please capture screenshots and
register them in `docs/tutorial/catalog.json` under the `08-17` entry's
`graphics` array, following the shape already used by the `08-16` entry
(`kind`, `path`, `caption`, `illustrates_step`, `capture_date`).

Most useful, in order:

1. **Step 2 — the inventory** showing all five starts with the P1 row's
   `missing_flanks` visible. This is the tutorial's main honesty point: GENtle
   declines a window it cannot build instead of padding one.
2. **Step 4 — the window before scoring**, with the three aligned rulers and the
   six transcript lanes. Establishes that curves are absent, not zero.
3. **Step 5 — the three curves**, ideally with one lane's Y-axis and printed
   score kind legible, to show the scales are separate and not cross-calibrated.
4. Optional, and more instructive than it sounds: **a window containing an `N`
   run**, so the amber unavailable band and its "upper half local +, lower half
   local −" convention are visible. This window has no ambiguous bases, so the
   tutorial cannot show it from TP73.

Captions should describe what is on screen and stay clear of implying occupancy.
A screenshot of a tall TP73 curve is not evidence that p73 binds its own
internal promoter.

## 5. What not to conclude

The tutorial is deliberately careful, and the screenshots should not undo that.
The three factors were chosen because the TAp73 → ΔNp73 feedback and E2F1's role
at this locus make them a coherent teaching set, with PATZ1 as an additional
candidate. That is motivation for the figure, not a result from it. Nothing in
this workflow measures occupancy, ranks the factors against each other, or says
anything about the unassessed P1 promoter.

If a reviewer reads the finished page as "GENtle shows that E2F1 and p73 bind the
ΔNp73 promoter", the page has failed and I would rather fix the wording than ship
the figures.
