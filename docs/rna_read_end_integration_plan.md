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
step 5 also requires 4. Each implementation commit updates this ledger,
changelog, roadmap and consolidated Glen checklist with actual checks and
remaining acceptance, not an inherited pass.

| Step | Deliverable | Implementation / Acceptance |
|---|---|---|
| 0 | Decisions, reconciliation, OD1-OD6 | Implemented; 19 documentation/checkout tests, status and whitespace pass; session-close 4 OK / 2 expected warnings |
| 1 | Annotation-checked core, uncollapsed projections, clip descriptors | Implemented over `1df2db13`; 20 authored Rust tests unrun, 31+6 shared cases checked by Python in LF/CRLF |
| 2 | Coverage/denominators, unavailable A/T, imported raw-signal evidence | Pure records/functions implemented; 8 authored Rust tests unrun; no file/service route yet |
| 3 | Engine-owned pileup/view and deterministic headless renderer | Pending |
| 4 | Shared supervised process and bounded configured samtools ingress | Pending |
| 5 | Authenticated shared commands, immutable exports and discovery | Pending |
| 6 | Thin virtualized GUI and migrated synthetic tutorial | Pending |

## Open Owner Decisions

Each remains **pending beyond the implemented default**. Sending the supplied
review authorizes only these defaults, not an implicit expansion later.

| ID | Options / Recommendation | Authorized Default |
|---|---|---|
| OD1 | Keep current clipping; later consider separately justified tail-only clipping policy. Recommend inspecting descriptors first. | No reinterpretation; every current status remains unchanged. |
| OD2 | New versioned SAM intake retaining bounded terminal bases (recommended), or FASTQ join with a new identity contract. | A/T channel unavailable: `no_bound_basecalls`; no scanner or guessed joins. |
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
