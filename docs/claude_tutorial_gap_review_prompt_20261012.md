# Claude Review Prompt: Five Tutorial Gaps

Historical prompt prepared 2026-10-12 at `150a339074e7…`. Steffen supplied
Claude's read-only reply on that same baseline and authorized implementation.
The original prompt below is retained as historical planning, not current
pending work. The reconciliation at the end records the revised scope; Claude
has not reviewed the resulting code or new tutorial prose. No direct contact
with Claude or Glen is claimed.

---

Please review the following bounded GENtle tutorial plan against the actual
repository. Read-only: do not edit files or run Rust builds/tests. Give a concise
critique and a revised minimal plan, with verified source references. Report
HEAD drift before interpreting old claims.

## Repository And Authority

Repository: `/Users/u005069/.codex/worktrees/47dd/gentle_rs`.
Branch: `gentle_rs_2_main`.
Reviewed preparation baseline:
`150a339074e7a2d282567d791a8f50834b70ff9c`.
The clean worktree was fast-forwarded from `b0b2703a` to local `main`; no fetch
was performed in this preparation session. Do not describe this as a freshly
verified remote tip. The primary checkout's `paper` and `output/` are unrelated
and must remain untouched.

Read AGENTS.md, the relevant sections of docs/architecture.md,
docs/decisions.md, docs/roadmap.md, docs/testing.md and current interface docs.
The owner leaves Rust verification to Glen/CI. Do not request a detached
multi-gigabyte build merely to review teaching prose. Fast Python source checks
may be suggested, but do not confuse them with executed GENtle workflows.

## Owner's Request

Address, as well as practicable, the five tutorial gaps listed below, and
prepare explicit prompts for Glen when native screenshots, public data
preparation or independent scientific acceptance would be cumbersome locally.
Use existing shared engine contracts. Do not add biological algorithms or
adjacent features just to make a tutorial possible. If a required capability
is genuinely absent, name the gap rather than inventing a command or making
the tutorial imply it is implemented.

The owner has requested your second opinion before implementation. Codex has
not yet edited tutorials, source metadata, fixtures or catalog projections.

## Current Inventory

The non-executing inventory on the baseline reports:

- 62 discovery-catalog entries and 30 executable-manifest chapters.
- Six typed GUI contracts: five `gui_scientific`, one `gui_view_only`.
- 29 manual guides without a typed GUI contract, 24 workflow-only entries and
  three reference entries.
- Review metadata: 46 Codex-reviewed, two historically human-reviewed and
  14 unreviewed. A review badge is not a current acceptance pass.

Recreate the inventory without executing GENtle:

```sh
python3 scripts/tutorial_acceptance.py inventory --output /tmp/gentle-tutorial-inventory.json
```

Catalog absence and manifest absence mean different things. Existing manual
guides should not be made executable by relabeling metadata, and an unexecuted
GUI contract must remain `not_run`. Do not re-propose already landed PATZ1,
reference-bound specificity, VKORC1, motif-score, TSS-collection or local TP73
curve tutorials as missing.

## Original Codex Plan

1. Extend 08.16 into a coherent per-TSS investigation using a new companion
   input with nonempty imported motif evidence, treatment/control signal,
   separate local scoring and an explicit reporter-planning continuation.
2. Add a repeat-family cohort lesson using UCSC rmsk, transcript-aware geometry
   and motif comparison. Explain RNA support as a separate ranking layer only
   where a real supported join exists.
3. Add a Microarrays lesson: configuration, file roles and retrieval sources,
   a public TP73 projection exercise, nearest-feature navigation and cautious
   interpretation of expression/probe evidence.
4. Add a catalogued gel-image sizing lesson using a deterministic drawn image,
   manual ladder assignments and independently calculable bp/kDa answers.
5. Add a real FLNA splicing-interpretation lesson, grounded in retained public
   annotation and exact loaded UniProt links rather than a synthetic gene
   misleadingly named FLNA.
6. Update task navigation and genuinely stale discovery descriptions only
   where required for these lessons. Preserve historical reviews and receipts.
7. Prepare separate Glen handoffs for fixture preparation, exact-source replay,
   native screenshots and scientific/readability acceptance.

Each lesson should start with the wet-lab question, prerequisites, available
inputs, expected observations, a small interpretation exercise and recovery
instructions. Put GUI, shared Shell/CLI and agent proposal examples together
where those routes exist. A headless UI-intent receipt is not a GUI action;
agent prose is not mutation permission. Prefer existing chapters and fixtures
over a parallel tutorial framework.

## Evidence And Questions By Lesson

### A. Integrated TSS Investigation

Read:

- `docs/tutorial/08-16_tss_regulatory_view_gui.md` and its matching source JSON.
- `docs/tutorial/08-17_tp73_dnp73_factor_curves.md`.
- `docs/tutorial/reproducibility/tss_regulatory_view_gui/`.
- `src/tss_sequence_view/`, `src/main_area_dna/tss_view.rs`.
- `crates/gentle-protocol/src/tss_profiles.rs` and `tss_motif_evidence.rs`.
- `docs/integrated_locus_tss_profiles.md` and existing reporter tutorials.
- `docs/regulatory_motif_subset.md` and
  `test_files/fixtures/regulatory_motif_subset/{README.md,make_fixture.py}`.

08.16 already loads a hash-bound report through the real project-hosted
`ui open tss-view --report ...` route. Its 11-base synthetic minus-strand record
has local TSS base 4 / genomic base 300, stored MA0004.1 peaks, a supplied
synthetic treatment lane and an unavailable control. Its raw negative scores
produce an explicitly clipped flat-zero display. It intentionally contains no
imported motif hits; this makes it insufficient for learning imported-hit
presentation, not a defective fixture.

`TssProfileReport.imported_motif_evidence` carries source and canonical-report
hashes. `tss_motif_evidence::validate/project/covers_window/score_range` already
enforces identity, projection and coverage rules. Complete curves retain
unavailable windows independently of zero; imported sparse hits are not curves.

Review the smallest way to add a nonempty, reproducible companion example.
Prefer an actual small supported query/export path over hand-editing generated
reports. Never manufacture a DuckDB execution receipt, producer revision,
binary hash or experimental control. If a hand-crafted viewer-only teaching
record is warranted, label it explicitly and keep real query acceptance separate.
Retain the old report, receipt and screenshots unchanged.

A real, tiny DuckDB/Parquet exercise already exists in the regulatory-subset
walkthrough. Its seven synthetic orientation records include negative and zero
scores, partial retention and shared TSS ownership. It uses MA0861.2/MA1961.2,
not 08.16's MA0004.1, and a different synthetic reference. It must not simply be
attached to 08.16's report. Review reuse of this preparation machinery and guide
before proposing a second DuckDB fixture generator or runner. Its scores are
hand-crafted identifiers/values, not outputs calculated from those real PFMs.

Check both genomic strands, local strand inversion on the minus window,
complete versus partial package retention and a footprint clipped at a boundary.
Use exact matrix accessions and distinct score units/scales. Treat attached
reports, local computations and imported hits as separate evidence sources.
Use 08.17 for the existing local-scoring path where possible. A reporter
continuation must name a real existing request and any missing control/approval;
do not imply that a selected TSS window is already a validated construct.

### B. Repeat-Family Investigation

Read `src/engine/analysis/repeat_cohort.rs`, the operations in src/engine.rs,
the corresponding parser arms in src/engine_shell.rs, and docs/cli.md's
`features repeat-query`, `repeat-cohort` and `window-cohort-tfbs` sections.
Inspect existing rmsk fixtures and their provenance before adding another.

The shared CLI already accepts class/family/name/alias filters and all four
geometry modes. Missing transcript/UTR/CDS context is reported explicitly.
Prepared genome/catalog/cache and rmsk assembly must be mutually compatible;
a small source fixture must not be presented as a genome-wide biological sample.

Important newly confirmed limit: BuildRepeatEnvironmentCohort has no RNA-report
argument, and row construction currently sets `rna_rank_score: None` and
`rna_support_labels: vec![]` (repeat_cohort.rs around 1093). Please verify whether
another existing route performs the join. If not, teach selection/geometry/motif
analysis now and mark RNA ranking as unimplemented, with an external
assembly/interval-bound evidence comparison only if honestly supported.
Do not introduce a hidden Python biological ranking engine or fake a CLI flag.

Advise whether a tiny prepared public/synthetic genome is the best offline
exercise, or whether this should initially be a prepared-resource walkthrough.
Preserve distinct repeat loci even when their sequences are identical.

### C. Microarray Evidence

Read docs/gui.md's Microarrays configuration and Array sections,
src/app/configuration_ui.rs, the sequence Array setup handler and docs/cli.md's
`arrays inspect-microarray-track` / `project-microarray-track` contracts.

The dedicated `Configuration -> Microarrays` tab and local resource discovery
already exist. The Array context action can focus the nearest projected feature.
Do not describe old scrolling behavior as a missing tab without checking source.

Reuse `test_files/fixtures/microarray_tracks/clariomd.tp73_vendor_subset.*`
and `test_files/tp73.ncbi.gb` if compatible. The nearby README states that
probeset IDs and hg38 coordinates are vendor-derived, but the effect statistics
are synthetic. Those are useful teaching inputs, not Rostock expression results.
Full CEL/APT/private publication resources must remain optional.

Explain the separate roles of vendor annotation, CEL measurements, sample
metadata, analyzed contrasts and a coordinate-bearing manifest. Point to the
actual retrieval instructions already documented; do not invent a vendor URL
or imply GENtle downloads/preprocesses everything automatically. Include an
assembly mismatch and 'no visible array feature' recovery checkpoint.

### D. Measured Gel Images

Read docs/gel_image_analysis.md, src/gel_image/tests.rs,
src/app/gel_image_ui_tests.rs and the shared gel-image commands.
An independent editor and CLI already exist; automatic band detection does not.

The documented 160 x 220 synthetic coordinate example has ladder centers at
y=20/80/140/200 with 2000/1000/500/250 bp. A sample at y=110 should estimate
sqrt(1000*500), approximately 707.1 bp. An outside-range mark must remain null,
not extrapolated. A small reproducibly drawn image would make this teachable.
Choose PNG for the basic lesson; keep 16-bit/TIFF and damaged-image checks as
existing contributor regressions rather than complicating the first exercise.

Explain original zero-based pixel coordinates, image digest, ladder identity,
confirmed assignments, localization bounds (not confidence intervals), no-clobber
exports and saved-project reuse. A separate SDS exercise reports apparent kDa,
not molecular identity or amount. Prestained ladders require the correct gel
system and supplier values; synthetic kDa labels are not a commercial calibration.

### E. Real FLNA Splicing Interpretation

Read docs/glen_flna_splicing_readability_handoff_20261007.md, the Splicing
Expert presentation contract and existing 07.03/06.02/04.08 lessons.

The historical public smoke returned FLNA ENSG00000196924, canonical
ENST00000369850.10, and loaded reviewed UniProt P21333. It reported 36
transcripts and genomic minus strand with gene-oriented local plus coordinates.
Those are properties of that retained run, not timeless expected counts.
The cited temporary state/view files no longer exist locally, so do not write
an offline lesson pretending they are bundled or reuse feature ID 300 as a
universal identity. The REST response did not pin an annotation release.

Recommend a minimal public pinned input bundle or a clearly optional online
preparation step, with exact retrieval/version/digest provenance and honest
current availability. Glen can prepare and review the real input if necessary.
Teach boundary motifs, exceptional/canonical flags, exact loaded-xref grouping,
unknown joins, exon selection and resulting protein consequences. Distinguish
loaded evidence from global UniProt absence and prediction from observed splicing.
Keep live 2400-px/native acceptance and human interpretation pending.

## Navigation, Generation And Validation

Review docs/tutorial/README.md and landscape_overview.md. The README's long
list has older literal review badges despite newer review_manifest metadata;
the landscape still calls an rs9923231 tutorial planned although 08.04 exists.
Decide the smallest scoped reconciliation, not a new graph/discovery project.

Use the existing source JSON/catalog generator, review manifest and parity
contracts. Do not hand-fix generated catalog/manifest outputs or manufacture
review dates. New manual guides need discovery metadata, not a fake executable
manifest entry. Add fast LF/CRLF coverage for any new byte-bound fixture or
generated-file freshness check, following scripts/test_tutorial_checkouts.py
and .gitattributes. Update the roadmap only for genuine remaining gates; it
has a 400-line limit and does not accept 'Done' entries.

Suggest exact focused Python/Rust/catalog/parity commands for the eventual
implementer and Glen. Keep Rust execution remote-only here. No build, benchmark,
release, scientific or native-Windows acceptance is claimed by this review.

Glen handoffs should request exact final source/binary/lockfile/profile/input
identities, replay results, untouched raw screenshots plus separate semantic
snapshots/annotated derivatives, hashes and typed failure classifications.
Existing screenshots are history, not replaceable placeholders. Ask Glen to
return a branch/full SHA, artifact locations and a concise per-check verdict
through the owner; do not assume direct agent-to-agent messaging is authorized.

## Requested Answer

1. Which proposed gaps are real, already covered, or blocked by missing inputs
   or capability? Cite source, not only prior commit summaries.
2. What is risky or unnecessarily broad in the original Codex plan?
3. Give a revised minimal sequence with files, fixture provenance, validation
   and one clear completion criterion per slice.
4. Separate immediately implementable tutorial work from data preparation,
   screenshots and scientific checks that Glen should perform.

Please do not implement. The owner will forward your review to Codex before
the tutorial changes begin.

## Review Reconciliation (2026-10-12)

**Original Codex plan:** cover all five gaps using shared contracts, but bundle
imported TSS hits, local scoring, controls and reporter continuation in A;
use a repeat-cohort exercise in B; add microarray/gel teaching and an FLNA lesson.

**Claude feedback supplied by Steffen:** A's attach route already exists, but
the default subset fixture cannot attach to 08.16 because reference identities
and matrix widths differ. The 11-base window cannot teach control biology.
B's RNA join is unimplemented and cohort/summary execution needs a real test.
C has existing vendor-coordinate inputs with synthetic statistics; D needs a
drawn input, not a virtual-gel screenshot. E lacks retained real FLNA inputs and
must remain a draft. Remove duplicated review badges instead of resyncing them.
Minimal order: navigation, D, C, compatible A, tested B, then input-gated E.

**Revised Codex implementation:** navigation/C/D are independently committed in
`d95286e31caa96a59e711c7504347c0a21009b58` with 79 Python checks and the separate
`glen_microarray_gel_tutorial_review_20261012.md` scrutiny request. A extends the
existing generator with a bounded `tss-view` preset rather than arbitrary
reference/width parameters: the two supported presets are auditable and cannot
silently create incompatible combinations. It adds real query/export regression
source and an optional 08.16 companion, not a second reader, scorer or control.
B adds distinct-locus fixture checks, shared-shell execution regression source
and a manual 08.18 guide with unavailable geometry and no RNA-ranking claim.
E remains outside the catalog as an input-pending draft.

Glen's current tasks and exact test commands live in
`glen_tutorial_gap_preparation_20261012.md`. Local Python/DuckDB fixture I/O is
not executed GENtle query/export or Rust acceptance. Old reports, screenshots,
human review dates and release verdicts remain unchanged. No push is authorized.
