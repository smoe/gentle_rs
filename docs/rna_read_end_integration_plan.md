# General Read-End, Coverage And Pileup Integration

Status: owner-authorized implementation, 2026-11-08. No local Cargo build,
check or test; no push, tag, dispatch or private-data processing. Rust/runtime,
native and scientific acceptance are separate evidence gates, never inferred
from Python/source checks.

## Reconciliation

The original five-step plan separated shared UTR geometry, end evidence,
indexed input, adapters/exports and a thin GUI. The owner-supplied Claude review
adds a mandatory per-query-unit pileup and identifies three missing contracts:
model/normalized-annotation agreement, uncollapsed hypothesis/candidate records
and additive clip evidence. Accept these additions without changing censoring.

The reviewed draft is our own work, committed as `1df2db13` on main over
`1b73f9dc`; it is not another session's uncommitted code. The first core's
17 Rust tests remain unexecuted, with 31 shared literal cases checked only by
Python. The refinement must not relabel that as Rust acceptance.

The normalized-annotation hash prerequisite is established in
`engine/analysis/rna_panel/primer_payload/observed_sam.rs::annotations`: normalized
bindings use the verified annotation file's hash; prepared bindings compare
the file against the preparation artifact and also use that file's hash.
`comparison::sources` carries it into `AnnotationBinding`. It is not the GTF
hash. Pure functions accept already authenticated inputs; the service must
verify original bytes before calling them.

Keep existing SAM v1, Phase 2/S1 and Python wires unchanged. No implicit gene,
transcript, sample, barcode, reference, assembly, biological threshold,
smoothing, imputation, representative isoform or fractional assignment.
Identity remains `sam_declared_query_unit`, not an original read or molecule.

## Commit Ledger And Dependencies

Dependency order: 1 -> 2 -> 3 -> 5 -> 6. Step 4 is independent of 1-3;
step 5 also requires 4. Step 7 extends 1/3/5/6 without changing their UTR rules.
Each implementation commit updates this ledger,
changelog, roadmap and consolidated Glen checklist with actual checks and
remaining acceptance, not an inherited pass.

| Step | Deliverable | Implementation / Acceptance |
|---|---|---|
| 0 | Decisions, reconciliation, OD1-OD6 | Implemented; 19 documentation/checkout tests, status and whitespace pass; session-close 4 OK / 2 expected warnings |
| 1 | Annotation-checked core, uncollapsed projections, clip descriptors | Implemented; included in Glen's scoped Linux read-end pass at `7c5214ec`; no biological/platform-wide acceptance |
| 2 | Coverage/denominators, unavailable A/T, imported raw-signal evidence | Implemented and service-connected; included in the scoped Linux pass; no inferred basecalls or real-signal acceptance |
| 3 | Engine-owned pileup/view and deterministic headless renderer | Glen reports renderer 6/6 at `7fa6b15a`; retain the earlier `01d42b71` budget receipt; native GUI separate |
| 4 | Shared supervised process and bounded configured samtools ingress | Linux real-samtools group 9/9 plus overlapping exact CRAM/FAI and BAM/CRAM reruns pass at `1c4a4a5a`; Windows separate |
| 5 | Authenticated shared commands, immutable exports and discovery | Retain earlier generator/parity/VM/MCP/export receipts; oversized refusal 2/2 and actual small conversion 1/1 pass at `01d42b71` |
| 6 | Thin virtualized GUI and migrated synthetic tutorial | Glen reports all-format actual CLI/export replay at `7fa6b15a`; native GUI/scrolling not_run; Windows CI passes, not a native GUI verdict |
| 7 | Original-read orientation and literal A/T runs (owner 2026-11-09) | Glen's scoped `7fa6b15a` return passes both checks, renderer 6/6, root/VM/MCP 47/47 and rehearsal 8/8 (both actual-process cases, zero skips); native and overall CI separate |

### External Compile-Gate Return, 2026-11-08

Glen's first default-feature generator build at
`1d5eed33bcd58460f93dd45f90d85588e9fb5859` failed with exit 101 on Debian/Rust
1.97.1. Nine diagnostics reduce to a missing GUI helper import, coverage's
inferred unsized accumulator and two moves of a still-borrowed publication
identity. The generator did not run; the empty matrix diff is not a freshness
pass. All later execution gates remain not_run.

The owner/Claude repair already exists in `6adee466`; source inspection is not
a successful rebuild. Restart the same default-feature canonical generator on
the repaired pushed source before the remaining checklist. Archive integrity
and the exact retry boundary are recorded in the
[consolidated handoff](glen_next_verification_handoff.md#read-end-compile-gate-return-and-retry-2026-11-08).
No implementation scope, biological contract or prior acceptance is changed.

### Local macOS Execution, 2026-11-08

Claude compiled and ran `7c5214ec` on macOS arm64 (rustc 1.101.0-nightly,
`--locked`, lockfile with bio 4.2.2) with two test-only fixes to
`indexed_locus_tests.rs`: the fixture stub timeout (2 s was shorter than the
first launch of a freshly compiled stub) and a 1 KiB hash cap that refused the
2016-byte synthetic FASTA. The canonical generator then produced the matrix diff
(two new rows plus counts). Counts are in the
[handoff](glen_next_verification_handoff.md#local-macos-execution-2026-11-08).
This is one platform, not CI, Windows, real-samtools or native acceptance.

At `1c4a4a5a` the same macOS host compiled all targets and passed the new
selection (9/9) and renderer layout regressions after two further test-only
fixes: the renderer test's attribute parser matched `y="` inside
`font-family="` and the supervisor fixture's 5 s default timeout was shorter
than the first launch of the changed stub. Counts are in the
[handoff](glen_next_verification_handoff.md#local-macos-check-at-1c4a4a5a-2026-11-08).

### Glen Linux Return And Remaining Repairs, 2026-11-08

At `7c5214ec`, Glen's canonical generation and generated-tree registry/freshness,
default/headless checks, read-end core/services, renderer, supervisor, glossary,
actual JS/Lua/MCP and CLI tutorial 07.01 scopes pass. The imported matrix diff
matches his patch exactly. Preserve those receipts rather than repeating their
campaigns. The [active handoff](glen_next_verification_handoff.md#glen-linux-read-end-return-and-bounded-repairs-2026-11-08)
records counts, hashes and the narrowed retry.

Step 4 remains blocked by the synthetic FASTA fingerprint cap and the real
samtools test's version-output UTF-8 failure. Preserve the existing test cap/
timeout fixes; the new lossless diagnostic escaping does not weaken SAM or CRAM
validation. Step 6/export presentation also needs fresh acceptance after shared
text wrapping and cached virtual-lane height changes. These bounded repairs
have new regressions but no Codex Rust/native execution. Windows, native GUI,
real BAM/CRAM, normal CI and scientific acceptance are not closed by the passes.

### Glen Linux Follow-Up At 1c4a4a5a, 2026-11-08

The follow-up closes Linux real BAM/CRAM selection, including the 9/9 configured
samtools group and overlapping exact CRAM/FAI/fallback and BAM-to-CRAM reruns.
Actual CLI export passes 1/1 with fresh JSON/TSV/SVG/PNG/PDF, bounds/no-clobber
and no project state. Glen's visual review is limited to those exported pages.

The renderer group is still 2/3: its test parser matches `y="` inside
`font-family="` and attempts to parse `monospace`. Preserve the existing
attribute-separator fix, add explicit attribute-suffix controls, and rerun only
the three renderer tests externally. This is test-only; no new BAM/CRAM or CLI
campaign is required. Archive hashes, exact counts and remaining native/Windows
boundaries are in the [active handoff](glen_next_verification_handoff.md#glen-linux-return-at-1c4a4a5a-and-renderer-test-follow-up-2026-11-08).
Codex ran no local Rust or native verification. Scientific acceptance stays open.

### Actual-Layout Pixel Guard Follow-Up, 2026-11-09

Glen's `a5496266` focused bundle passes the three renderer tests and bio policy
checks. The broader source scan found a distinct guard drift: the export's fixed
54-pixel row estimate omitted wrapped header/lane height. The repair exposes
dimensions from the same checked layout used for SVG, checks actual width times
height before PNG/PDF conversion and retains the 32-million-pixel threshold.
Overflow and long-header/multiline-lane refusals plus actual small PNG dimensions
are authored regressions, not locally executed Rust evidence. SVG/JSON/TSV stay
outside the pixel cap under their existing output limits. Exact source, archive
hash and narrowed external commands are in the
[active handoff](glen_next_verification_handoff.md#read-end-actual-layout-pixel-guard-2026-11-09).
Earlier passing scopes, native/Windows boundaries and OD1-OD6 remain unchanged.

### Pixel-Guard Linux Return, 2026-11-09

Glen's verified `01d42b71` package passes renderer 4/4, budget 2/2, small actual
PNG/PDF conversion 1/1 and default test-target Cargo check, without failures or
skips. This closes the focused Linux pixel-guard execution gate. Its carried
PDF visual reference remains attributed to `1c4a4a5a`, not a freshly retained
`01d42b71` PDF. The existing same-source normal CI additionally passes Windows
and full Linux, headless and release safety; macOS was skipped. Native GUI and
release acceptance remain separate. Full evidence identity and the unselected
point-release inventory are in the
[current handoff](glen_next_verification_handoff.md#glen-pixel-guard-linux-return-and-release-inventory-2026-11-09).
No local Rust execution, candidate freeze, dispatch or publication by Codex.

### Oriented Original-Read Extension, 2026-11-09

The owner's forwarded product contract authorizes Step 7 and replaces only
OD2's unavailable-only default with an explicit hash-bound original-base input.
No automatic SAM/FASTQ capture/join, clipping reinterpretation, tail caller,
private analysis or new registry path is added. V2 pileups preserve v1 reading,
original-query literal A/T identity, all qualifying runs with terminal/internal
classification, direction alternatives, covered intervals and unfiltered unit
counters. Shared scenes group mRNA-direction lanes above/below a divider and
keep conflicts in Alternatives. The same renderer supplies native/export row
heights, glyphs and actual-layout pixel checks. The
[measurement contract](rna_read_end_metrics.md#original-bases-and-oriented-single-read-lanes)
specifies supplied orientation/provenance versus verified byte bindings.

The earlier `01d42b71` pixel/platform receipts remain historical, not acceptance
of the extended renderer/service. New source, VM/MCP and actual-process gates are
in the [active handoff](glen_next_verification_handoff.md#oriented-original-read-pileup-2026-11-09).
Implementation is committed on main at the owner's request. No local Rust/CLI
execution, push, merge of another branch, candidate selection or scientific
approval is claimed by source implementation.

### Glen Oriented-Pileup Return, 2026-11-09

Glen reports the complete focused checklist at exactly `7fa6b15a`: both locked
checks, renderer 6/6, root/real-VM/MCP group 47/47 and rehearsal 8/8 with both
CLI process cases and no skips. No local product defect was found. The
[source-bound return](glen_next_verification_handoff.md#glen-oriented-pileup-return-2026-11-09)
records the reported binary hash, locally confirmed lock hash, available CI
snapshot and absent raw-receipt details. Native GUI is not_run; Windows,
headless and release-safety CI pass while macOS is still running and full Linux
was skipped. This closes the scoped execution request, not release acceptance.

The separate exon-locus/shared-genome-track/scored-Assignments mode is now
source-authored below, not covered by this historical result. Its prompt/tests
11-21 remain general, not gene-specific, and do not extend Step 7's acceptance.

### Exon-Locus Extension Draft, 2026-11-09

Status: implemented in source at the owner's request; Rust/runtime/native
verification pending, not a new accepted slice or release gate.
Canonical [owner-forwarded Glen prompt](https://gentle.functional.domains/openclaw/canvas/codex-exon-locus-prompt-20261109-d72152e0/report.md),
SHA-256 `d72152e0197392e2f820312e5a338e860e56cc7a5737f94ef057fa610d79df6d`;
the Markdown bytes and its manifest entry were checked through normal host
HTTPS, without bypassing the allowlist. The owner supplied Claude's source-only
review on 2026-11-10; the resulting bounded repairs are recorded below. No local
Claude CLI invocation or Rust/runtime acceptance is claimed. The original
[review prompt](exon_locus_claude_review_prompt.md) remains historical.

The current `coverage::gene_assignment` intentionally excludes multi-gene units.
Do not weaken it to admit a sense/antisense locus: new locus membership must
start from authenticated saved genomic observations, before transcript/gene
assignment. Keep the verified v1/v2 transcript path and its reports unchanged.

1. **Pure locus and assignment records.** Explicit exon interval, contig and
   bounded flanks; versioned annotation/model bindings; one record per declared
   unit, factoring the existing sample/unit read-key recipe into one shared
   helper for both projections. Keep
   every legal hypothesis and off-locus connection, and distinguish supplementary
   flags from actual chimeric geometry. Add all seven requested compatibility
   states plus machine-readable evidence/conflicts and a versioned score tuple.
   Order lexicographically by splice compatibility, discriminating exonic
   structure, supported ends, then weak proximity. Lower evidence must never
   compensate for a structural conflict. Stable-ID ordering is presentation,
   not a winner; tied or uninformative models remain alternatives. Raw BAM
   direction and uncalibrated A/T cannot decide RNA origin. Unknown motif
   evidence remains unknown; do not label a novel junction canonical from
   annotation proximity. Author the independent fixtures for tests 11-17/19/20.
2. **One authenticated service and two projections.** Add versioned locus
   actions under `RnaReadEnds`, reusing original-base intake, bounded loading,
   cancellation/approval, owned publication and converters. Do not synthesize
   an old S2 compatible-model claim to reuse the gene-filtered path. Build
   genomic evidence once; a selected transcript changes only its explicit
   projection, never read identity, locus geometry, candidates or input counts.
   Declare each coordinate system. Missing annotation leaves genomic viewing
   available with unavailable assignments/UTR; failed raw hashes still refuse.
   Changed annotation gets a new binding/recomputed scores, never a stale cache
   or silent upgrade of historical reports. Cover tests 18/20/21.
3. **Shared scene and thin native controls.** Fixed genomic locus, both-strand
   annotation tracks with arrows and gene/model/version IDs, one displayed
   read row with separately identifiable hypothesis geometry, and edge markers
   for connections leaving the locus. Genome and transcript modes consume the
   same evidence. Preserve all candidates in JSON/TSV and inspectable reasons in
   the viewer; no coverage-only substitution. Reuse virtualization and checked
   actual dimensions. Keep full/partial/no reference overlap explicit rather
   than silently changing the meaning of v2's full-coverage boolean. Native
   acceptance requires a real GUI environment; tests 8/18/19/21 remain scoped.
4. **Process/adapter acceptance and documentation.** Source-matched synthetic
   CLI, real JS/Lua/MCP, both projections and all five formats; deterministic
   candidate ordering/exports, raw LF/CRLF bindings, quoted/Unicode paths,
   no-clobber/racing destinations and unchanged source/project state. Extend
   per-sample declared/locus/visible counters and UTR/alignment/run summaries
   without counting hypotheses/models twice. Recheck the affected tests 1-10
   alongside 11-21, not unrelated completed campaigns. Glen/CI owns Rust and
   native Windows/GUI execution. Update discovery and canonical projections
   only where their contracts actually change; never hand-edit generated rows.

No numeric biological probability, read-to-molecule conversion, implicit
representative model, TP73-specific label heuristic, calibrated adapter rule,
raw-signal tail caller or private-data run is part of this extension. Do not add the
new mode to `.13` acceptance merely because its prompt is available; candidate
scope and publication still require the owner.

Implementation follows those four parts: pure `read_ends/locus.rs` and
`locus/presentation.rs`; nested `rna_read_ends/locus.rs` actions; existing shared
renderer/native viewer; synthetic Rust and opt-in real-process/VM fixtures.
The old gene-assignment rule and v1/v2 wires remain unchanged. Preferred flags
use structural compatibility/discrimination, not weak endpoint proximity;
the full tuple orders presentation without forcing a unique biological winner.
Uncalibrated A/T has zero strand-assignment weight. A novel junction is not
called canonical without motif evidence; incompatible hypotheses remain unresolved.
Missing explicit annotation/CDS leaves genomic evidence and unavailable facts;
present-but-mismatched files still fail closed. Scores are freshly computed
against the supplied annotation binding, never taken from old comparisons.
Bounded retained geometry, candidates, runs, work and page sublanes refuse
rather than silently truncate. All saved units are per SAM, not molecule IDs.
See the [active verification scope](glen_next_verification_handoff.md#exon-locus-and-scored-alternatives-2026-11-09).

### Exon-Locus Review Repair Ledger, 2026-11-10

Source review against `7fa6b15a` plus this thread's uncommitted implementation;
the named structures/counterexamples were traced before changing them. These
are source-derived defects, not reproduced Rust/runtime failures. Repairs are
authored with Rust tests **unrun at authoring**. The later
[attributed macOS return](glen_next_verification_handoff.md#claude-exon-locus-return-at-b95aae07-2026-11-10)
records scoped execution and a test-only anchor correction, not current-lock,
default-feature, VM, actual-CLI, Windows or native acceptance.
No schema/RULE_VERSION bump for this unreleased v1, new dependency, additional
calibrated cutoff, mask policy, queue, release gate or legacy gene-only behavior.

1. Outside-window run anchors remain unmodified data/off-axis text, not invalid
   or clamped glyphs; JSON-only export has a service regression too.
2. Junction-backed conditional strand support; mirrored identical chains remain
   competitors. Preferred flags/per-model classes and overhang refusal stay.
3. Full-chain classification follows complete conflict-free annotated junctions,
   independent of deletion gaps and endpoints. Unavailable exclusions use serde
   names; raw-only locus membership is separately marked/counted.
4. Cache other-model exon unions once per sample/model, adding consumed `work`
   without tuning caps/charge rules or introducing a sweep algorithm. A
   30-model/12-exon/20-read/~40-block fixture requests the 100M cap and prints
   actual consumed work when externally run; no local measurement is invented.
5. Distinct declared-not-found diagnostics, de-duplicated run counters,
   alternative UTR statuses/reasons, unique annotation tracks and selected-track
   placement. Add realistic query-disjoint chimera, novel, CDS, two-sample,
   cancellation and named-budget regressions. Replace Python source grep with
   fixture behavior; actual GENtle LF/CRLF evidence remains the opt-in CLI test.

The two withdrawn claims are not acted on: query-overlapping segments are never
joined by S2's DAG, and admitted units already have known axes. Mask-aware
counts, terminal tolerance/calibration, run-storage redesign and capacity tuning
are outside this repair. The baseline diff in `read_ends.rs`, `oriented.rs` and
`pileup.rs` was reviewed as shared helper extraction/default-empty scene fields,
with no changed legacy charge, hash recipe or v1/v2 report field. This is a
source-level preservation check, not a runtime byte-equivalence receipt.

| Decision | Options (Owner Choice Pending) | Implemented Default |
|---|---|---|
| OD-L1 | Require junctions for conditional strand support, or later permit independently justified containment-only support; containment-only `AntisenseCompatible` class policy remains pending. | Every preferred assignment needs a matched junction for `model_supported`; otherwise strand unresolved with the containment-only reason. No score/preferred/class change. |
| OD-L2 | Tolerate declared-but-absent annotation/CDS with diagnostics, or refuse those declarations. | Preserve tolerance with distinct hash-bearing diagnostics; present bad hashes and absent declared originals still refuse. |
| OD-L3 | Keep terminal overhang as a structural conflict, or separately justify tolerance/calibration. | Unchanged; no new threshold or accepted relation. |
| OD-L4 | Measure current capacity, or later authorize cap/charge/algorithm changes. | Reported synthetic 30-model/20-read work is 5,543,670 at `b95aae07`; distinct-isoform and real capacity remain unmeasured. No limit/charge/algorithm change. |

Follow the [narrowed repair checks](glen_next_verification_handoff.md#exon-locus-review-repairs-2026-11-10),
not unrelated completed reader/BLAST/barcode campaigns. Native, Windows and
scientific acceptance remain separate. All four owner options remain pending
beyond these expressly selected defaults.

## Open Owner Decisions

Each remains **pending beyond the implemented default**. Sending the supplied
review authorizes only these defaults, not an implicit expansion later.

| ID | Options / Recommendation | Authorized Default |
|---|---|---|
| OD1 | Keep current clipping; later consider separately justified tail-only clipping policy. Recommend inspecting descriptors first. | No reinterpretation; every current status remains unchanged. |
| OD2 | Owner-authorized original-base document (2026-11-09); automatic SAM/FASTQ joins remain deferred. | Optional file binds exact saved analysis/unit/QNAME/full query axis; unavailable otherwise. Literal 60/10 terminal defaults; no inferred original or tail length. |
| OD3 | Separate imported file (recommended initially), or future caller BAM tags bound to provenance. | File only, explicit caller/version/model/kit/configuration; no aux-tag reading. |
| OD4 | Existing per-analysis key, or future BAM-bound identity (recommended for cross-locus aggregation). | Deduplicate within one analysis only; refuse multiple analyses per sample/gene. |
| OD5 | Freeze Python as reference oracle (recommended), or explicitly migrate/deprecate later. | Python unchanged; ranking/contrasts stay there, no new features ported. |
| OD6 | Keep V2 JSON/preparer (recommended initially), or version the GTF preparation kernel to retain CDS/stops. | Rust reads V2 JSON; no third GTF parser. |

## Core And Report Invariants

Every candidate used for a metric must match its authenticated normalized
annotation transcript's contig, strand and genomic exons, and its model gene
must be an annotation gene. Mismatch refuses; missing CDS geometry remains
unavailable. Saved annotation gene IDs and complete-model gene IDs are separate.
The core emits every optimal-hypothesis x compatible-candidate projection,
including intervals, saved relations, both metrics and clip length/kind.
Collapsed metrics derive from those records, never a frontend reconstruction.

Coverage and pileup retain one lane group per declared unit, with every optimal
hypothesis shown separately. Gene-assigned lanes require nonempty compatible
sets in every hypothesis and exactly one unioned annotation gene. Exclusions
are explicit: multi-gene, no-compatible-transcript, no-hypothesis or unavailable.
The selected-transcript subset splits unique/shared/hypothesis-dependent.
Filters never change either denominator. Explicit deterministic ordering,
tie-breaks, ranges and lane limits prevent silent sampling/truncation.

Saved blocks, actual N skips, deletion gaps and links between observations
remain distinct. Clips and end-evidence glyphs are off-axis, not aligned bases.
No new strand or molecule claim comes from query/model orientation.

The samtools adapter is explicitly configured and local only. Missing tool or
index refuses; index warnings refuse. Region selection never declares complete
search/views. No whole-input pass occurs unless separately requested. CRAM
requires its bound reference/FAI and an isolated reference lookup environment.
The existing declared conversion variant and its bytes keep their old meaning;
executed-selection provenance is a separate new contract.

All surfaces consume the engine-owned report/view. JSON/TSV/SVG share values,
hashes and selection rules. PNG/PDF derive from SVG with existing renderers.
Use owned staging, closed/synced handles, no-clobber publication and report-last
completion. Authenticate saved bindings/artifacts before interpretation.
Confirmation, nested-agent guards, stale/cancel handling and project separation
remain mandatory. Canonical projections must be generated, never hand-edited;
if a source-matched generator cannot run under the build restriction, record
that gate as pending rather than inventing freshness.

The tutorial catalog was regenerated with the existing canonical docs binary
reporting source `2d2b8b42`; `src/workflow_examples.rs`, its child sources and
`src/bin/gentle_examples_docs.rs` are byte-identical between that source and the
integration. Catalog/manifest checks pass (70 entries / 30 executable chapters).
That narrowly permits metadata generation, **not** acceptance of the newer core
or registry. The parity generator embeds registry/operation/glossary code, so
its source-matched regeneration remains an external prerequisite. No matrix was
hand-edited. All Steps 1-6 now have source implementations; runtime/VM/native,
actual-tool and scientific acceptance remain open, not satisfied by that ledger.

## Verification

Author focused Rust tests per step; Glen/CI compiles and executes with an
explicit run-owned target. Local checks are Python, formatting and documents.
Shared corpus literals cover both strands, CDS/stop conventions, missing
models, zero/clip/ambiguity and geometry refusals. Report/pileup tests add
denominators, stacked hypotheses, connectors, sort/filter/paging, unavailable
evidence and deterministic bytes. Process fixtures cover cancellation, caps,
timeouts, indexes, CRAM references and forbidden fallback without real data.
Byte-bound fixtures/projections require scoped attributes and actual LF/CRLF
checkout coverage. Native responsiveness, Windows handles, actual configured
samtools and private/biological acceptance remain independent, explicit gates.

Do not repeat completed BLAST/reader, containment or private 15/16 campaigns.
Retain their original receipts. Finish each handoff with
`python3 scripts/maintenance_chore.py session-close --plan docs/rna_read_end_integration_plan.md`.
