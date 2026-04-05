# Release Notes: `v0.1.0-internal.4` (draft)

This draft internal release is centered on a new sequencing-confirmation track,
substantial RNA-read mapping refinement, broader exploratory visualization via
multi-sequence dotplots, and the first concrete step from a monolithic crate
into a structured `gentle-*` workspace. Compared with
`v0.1.0-internal.3`, this cut is less about polishing one existing window and
more about adding an adjacent validation workflow: from sequence design and
variant retrieval toward trace-backed confirmation and reportable evidence.

As with the earlier internal releases, these notes intentionally describe both
what is already usable and what is visibly becoming stable enough for broader
internal review.

## Highlights

- Sequencing confirmation is now a real cross-interface feature track:
  - raw Sanger trace intake landed
  - trace-aware confirmation logic is in place
  - sequencing-primer suggestion/proposal paths were added
  - GUI review now includes chromatogram-backed inspection and saved reports
- RNA-read mapping gained a more explicit second-stage review story:
  - phase-2 alignment can now be rerun over retained reports
  - individual retained rows can be selected and aligned on demand
  - histogram / score-density exploration and retained-hit inspection are more
    explicit and better instrumented
  - gene-level support summaries are now available from the shared shell/CLI
- GENtle now has a multi-sequence dotplot workflow rather than only one-query
  pairwise views.
- Variant / dbSNP retrieval became much more cloning- and review-friendly:
  - rsID lookup can retrieve genomic context directly
  - the resulting extracted locus is represented in the DNA window with a
    proper single-position variation feature and chromosome-aware interpretation
- The repository has started the planned workspace split:
  - `gentle-protocol`, `gentle-render`, `gentle-shell`, `gentle-gui`, and
    scaffolded `gentle-engine` now exist
  - stable public contracts and some rendering/help/backdrop layers have moved
    into those crates without forcing a big-bang import rewrite
- Help/tutorial behavior and hosted-window behavior were hardened across the
  same period, reducing friction during internal demonstrations.

## Notable Changes by Area

### 1) Sequencing Confirmation and Sanger Trace Review

- The sequence-confirmation track is now visibly established rather than only
  planned:
  - initial GUI for sequence confirmation landed
  - sequencing-primer suggestion/proposal paths were added
  - trace-aware confirmation became operational
  - final report-oriented sequencing-confirmation records were added
- Raw Sanger trace intake now exists:
  - AB1-style trace input support was introduced
  - GUI metadata and trace handling paths were added
  - trace fixtures and provenance notes now exist under
    `tests/fixtures/sequencing_confirmation/`
- The GUI now includes a chromatogram-backed inspection surface:
  - selected loci can display chromatogram curves
  - review is variant-focused and confirmation-oriented rather than a generic
    full trace editor
- Cross-adapter parity improved around this track:
  - CLI/tutorial coverage for sequencing confirmation was added
  - sequencing-primer overlays now have parity-oriented work behind them

### 2) RNA-Read Mapping and cDNA Review

- The RNA-read mapping workflow moved closer to a report-centric specialist
  surface:
  - retained reports can now enter an explicit alignment phase after phase 1
  - individual retained reads can be re-aligned selectively
  - score-density views and retained-hit review are more explicit
  - histogram-based read exploration is now available
- The shared shell/CLI gained a gene-level summary operation:
  - `rna-reads summarize-gene-support` can now summarize support for one or
    more genes from a saved report
- RNA-read UX was clarified further:
  - selection hover text and subset-alignment behavior were tightened
  - later internal review can now inspect retained rows with clearer intent
    rather than treating the RNA mapping run as one opaque pass
- This release continues to strengthen the “retain first, inspect later” story
  for RNA-read workflows, even though long-read alignment interpretation still
  needs more biological tuning.

### 3) Multi-Sequence Dotplots and Exploratory Visualization

- GENtle now supports multi-sequence dotplots, extending beyond the prior
  single query/reference emphasis.
- This broadens exploratory comparison use cases for transcript fragments,
  retained read subsets, and validation-oriented sequence comparisons.
- The feature aligns with the larger trend in this release: more specialist
  review surfaces that sit beside design/planning rather than only initial
  construction workflows.

### 4) Variant Retrieval and dbSNP Context

- A dedicated rsID-oriented retrieval path landed:
  - genomic context of dbSNP rsIDs can now be retrieved
  - chromosome alias handling was fixed/improved, including RefSeq-style
    accessions
  - extracted loci now include a concrete variation feature at the SNP
    position
- The GUI `Fetch GenBank / dbSNP...` path became more robust and easier to use
  during internal demonstrations.
- Tutorial/figure work also advanced around this area, including a first more
  concrete ClawBio-side artifact and a dedicated dbSNP presentation figure.

### 5) Help, Tutorial, and Hosted-Window Reliability

- Hosted Help/Tutorial windows received several rounds of fixes:
  - stale layer issues and problematic viewport behavior were addressed
  - Help no longer tends to open full screen by default in the same way as
    before
  - hosted project/help window behavior in the root workspace was restored and
    hardened
- The release also carries more tutorial material:
  - separate arrangements tutorial
  - sequencing confirmation CLI tutorial
  - continued figure/tutorial work around genomic variant retrieval and Gibson
    arrangements
- These changes improve internal onboarding and make it easier to demonstrate
  partially mature workflows without fighting the windowing system.

### 6) Workspace Split and Internal Architecture Progress

- The planned workspace decomposition is now visibly underway:
  - new workspace scaffold
  - initial separation into crates
  - subsequent extraction phases for protocol, render, shell, and GUI layers
- Concretely:
  - `gentle-protocol` now owns more stable analysis/report/display contracts
  - `gentle-render` owns feature-expert SVG, protocol-cartoon, and pool-gel
    rendering layers
  - `gentle-shell` owns shell help rendering
  - `gentle-gui` owns window-backdrop/resource helpers
- This is mainly an internal-engineering milestone, but it matters for release
  readers because it reduces future refactor risk and makes cross-interface
  parity work easier to sustain.
- The monolith is not gone; this release simply marks the point where the split
  became concrete rather than aspirational.

### 7) Packaging, CI, and Release Hygiene

- GitHub CI packaging saw additional fixes, including a `.dmg` target-path fix.
- Docker/build context handling was updated so the new workspace crates are
  copied correctly into image builds.
- `gb-io 0.9.0` transition work was completed, which is relevant both for build
  hygiene and for keeping sequence/annotation import behavior current.
- Shell command dispatch was also split to prevent recursive-shell stack
  overflows, which is primarily a safety/stability improvement for agent and
  workflow execution.

## Interim Release Readiness

- This draft release tells a much stronger validation story than
  `v0.1.0-internal.3`:
  GENtle is no longer only about sequence planning and feature visualization,
  but now also about bringing confirmation evidence back into the same shared
  environment.
- The strongest new narrative thread is:
  1. fetch or derive a sequence/locus,
  2. inspect/edit/plan,
  3. retrieve variants or supporting genomic context,
  4. import sequencing evidence,
  5. review confirmation state via reports and traces.
- RNA-read mapping also continues to move from “experimental filter pane”
  toward a reusable report workflow with explicit follow-up analysis steps.
- The workspace split is not release-marketing material by itself, but it is a
  meaningful readiness indicator: the codebase is being prepared for growth
  without abandoning the existing shared-engine contract.

## Release-Facing Known Limitations

- Sequencing confirmation is still strongest for:
  - called-sequence confirmation
  - primer suggestion/proposal
  - variant-focused chromatogram review
  It is not yet a general chromatogram editor or full base-calling workbench.
- RNA-read mapping still has important interpretation limits for long reads:
  - retained-report alignment is now explicit, but biological tuning for long
    composite/partial reads is still incomplete
  - the current seed/alignment model is useful for structured internal review,
    but not yet a final isoform-resolution endpoint
- Alternative cloning plans are still not ranked by practical complexity,
  elapsed time, reagent cost, failure risk, or local feasibility.
- The workspace split is intentionally partial; most production logic still
  lives in the root crate and is only beginning to migrate.
- Internal tutorial/figure coverage improved, but some new specialist features
  still outpace polished end-user walkthrough coverage.

## Install / Package Notes

- Version in `Cargo.toml` has not yet been bumped beyond
  `0.1.0-internal.3`; these notes are therefore a draft for the next internal
  cut rather than a statement that the repository is already retagged.
- Expected artifact story remains:
  - macOS `.dmg`
  - Windows `.zip` containing `gentle.exe`
  - Linux release metadata / container path as documented elsewhere in the
    release process
- Before tagging, verify that release automation still behaves correctly under
  the new workspace layout.

## Release-Shaped Smoke Check Results

Pending for this draft. Before tagging, run the local/internal release-shaped
matrix from [`docs/release.md`](docs/release.md) and record pass/fail here.
