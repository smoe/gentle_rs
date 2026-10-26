# Primer Payload Partner Evidence: Phase 2 Implementation Plan

Status: experimental implementation authored; Rust/runtime and scientific
acceptance are **pending**, not completed.
Prepared 2026-10-26 against local main
`b544a21b81b0d98896599f4ae21cdf60eea59b62`. This is a source-inspection baseline,
not a new candidate freeze or a claim about current cross-platform CI.

Input: Glen's owner-supplied
`GENtle_PRIMER_PAYLOAD_PARTNER_PHASE2_PROMPT_2026-10-26.md`, SHA-256
`501e107849f8a548786a9b17bca2fc5cc1c413f10b144d1762fbd8c61de8ed2d`.
The original Downloads document is unchanged. A read-only Claude consultation
was offered; this plan has not been Claude-reviewed.

## Implementation Boundary, 2026-10-26

The separate schemas, verified retained-primary/explicit FASTQ-replay preparation,
query-grouped PAF/annotation/streamed-reference import, best-equivalent segment
evidence, partner profiles, separate edge clusters, named artifact publication,
shared-shell/MCP/agent routes and synthetic rehearsal are authored. See the
[actual contract and claim prerequisites](primer_payload_evidence.md) and
[Glen's incremental execution checklist](glen_next_verification_handoff.md#primer-payload-phase-2-2026-10-26).
No hard-coded current-panel runtime policy or old Phase 1 schema change.

Bounded first-version choices: one batch per explicit reference per request;
query-grouped PAF rather than an implicit external sort; a nonoverlapping
query-span DAG with union support, not a second aligner; JSON-cell TSVs for
lossless nested ambiguity; encoded-data ceilings in addition to caller budgets.
Overlapping segments are alternatives, not silently merged into a biological
chain. Inspection verifies content hashes only when explicitly requested.

Twenty new Rust regressions are authored (16 engine, four adapter/lifecycle/
feedback), with an opt-in real-CLI LF/CRLF/spaced-path replay. They have **not run
locally**. The wider mandatory matrix below, full CI, canonical Rust-generated
parity freshness, Windows handles, live-agent responsiveness and private-cohort
accuracy/throughput remain external acceptance work. No item is accepted merely
because source/Python checks pass. The existing goal stays paused.

## 1. Outcome And Limits

Interpret externally aligned, primer-free payloads as a separate engine-owned
evidence record. Preserve all supported partner alternatives, target-base
unions, ordered segments, concrete primer-adjacent edges and recurrence counts.
Do not replace the Phase 1 local core/probe observation or RNA-panel nomination.

The main operation is `SummarizePrimerPayloadAlignments`, with
strict versioned request/report/result-binding/inspection records. Query
preparation needs its own small shared operation; it must not be a hidden
adapter script. Inspection and request-bound cancellation reuse shared routes.

No aligner execution/download, SAM/BAM, GUI, primer redesign, UMI inference,
fusion/trans-splicing confirmation, private-data run or release acceptance is
included. External alignments are evidence supplied to GENtle, not independently
reproduced alignment searches. Hash agreement is not authenticated provenance.

## 2. Generality Is An Acceptance Condition

- No runtime branches for a gene, species, primer label, barcode, contig name,
  fixed panel size, commercial kit, ordered-oligo length or current cohort.
- Primer sequences/components/direction come from the bound Phase 1 panel;
  inherit its declared limits rather than imposing a new 44-mer assumption.
- Expected targets require explicit stable-ID mappings between Phase 1
  locus/member identities and the new reference catalog. Do not equate a panel
  locus ID with a reference transcript or infer equality from gene symbols.
- Reference, annotation, cohort selection and producer identity are inputs.
  No automatic chromosome aliases, transcript-version stripping, reference
  downloads or resolution from labels. Ambiguous mappings remain sets.
- Biological thresholds and equivalence/clustering policies are visible,
  content-bound request fields. Seed lengths are a list, with 6/8/10/12 as the
  requested test profile, not the only admissible lengths. Short cores or
  unavailable context yield explicit unassessed results, never padding.
- Tests must vary primer count, oligo/core lengths, orientation and reference
  layout. Consistent opaque-ID renaming must preserve biological outcomes after
  inverse renaming; input-order permutations must not choose a different winner.

## 3. Source-Grounded Contract Corrections

Existing seams: [Phase 1 plan](primer_trace_phase01_plan_20261026.md),
[protocol](protocol.md#primer-anchored-read-trace-phase-01),
`src/engine/protocol/primer_trace.rs`,
`src/engine/analysis/rna_panel/primer_trace.rs`, and its `matching.rs`/tests.

- The Phase 1 inspector reads a bounded TSV prefix and reports
  `artifact_content_verified=false`. It is not sufficient input verification for
  Phase 2. Fully verify the selected report's allowlisted sidecars and lengths
  before admitting their records.
- Its optional payload FASTA contains the selected primary-anchor payload,
  normalized into primer direction. It does not export every retained hit.
  Missing FASTA or secondary-anchor sequence cannot be reconstructed from the
  report's coordinates alone. The query-preparation slice below closes this
  gap without modifying any Phase 1 schema, historical artifact or hash.
- Query identity needs trace digest, file ID, record index, read ID and stable
  anchor identity. Repeated read names, repeated primers within a read and
  competing primer assignments must not collide. Multiple PAF rows for one
  query are normal; duplicate manifest identity is a different error.
- Distinguish finalized report status from complete input coverage. Partial
  Phase 1 reports may supply only their EOF/hash-verified units, with retained
  outcome and missing-input diagnostics. Never silently promote them to a full
  cohort. The current Phase 1 fixture repairs still need Rust execution.

## 4. Bound Inputs And Coordinates

Create new strict schemas, leaving existing ones byte-compatible:

- Full `RnaPrimerTraceResultBinding`, selected unit receipts, query manifest
  and query FASTA hashes. Retain the exact source outcome and artifact receipts.
- Explicit cohort identity, full/selected coverage, selection-rule provenance,
  selected query identities and exclusions with reasons. Distinguish cohort
  scope from reference scope; neither is inferred from filenames.
- Per alignment batch: PAF path/raw SHA-256, query-batch digest, reference ID,
  transcriptome/genome scope, producer name/version/argument vector, exit and
  completion declaration, and coverage accounting. A canonical argument digest
  may bind arguments but is not a command to execute.
- Reference catalog: opaque reference ID, assembly/release, exact FASTA path
  and raw digest, sequence-name/length bindings, normalized transcript-to-gene
  associations and their source/raw/semantic digests. Genome evidence also
  requires declared exon geometry, genomic strand and transcript projections.
- Explicit expected-target crosswalk, policy, resource budgets, fresh output
  directory and optional request-bound cancellation marker.

Use zero-based, half-open intervals internally for raw-read, payload-query and
reference spans. Interbase edge coordinates are separate typed values. Retain
strand and coordinate-space identity; never mix transcript positions with
genomic positions or read orientation with transcript strand. Project genomic
donor/acceptor edges using the bound annotation, not a gene-name heuristic.

The query manifest records its exact sequence length/hash, raw-read interval,
anchor/core intervals, orientation transform and any excluded primer intervals.
Query coordinate zero is the anchor-adjacent payload edge in primer direction.
Reference sequence context is required for terminal-seed checks; PAF coordinates
or approximate identity alone cannot supply missing reference bases.

## 5. Query Completeness And Aligner Evidence

A manifest defines the intended query universe. Every selected query/reference
batch must be accounted for as mapped, producer-reported unmapped, or explicitly
unassessed. Missing accounting is an error in strict mode, not zero support.
Foreign identities, length conflicts, duplicate query definitions or inconsistent
mapped/unmapped declarations are hard failures. Known exclusions stay named and
outside assessed denominators. Valid partial-cohort analysis requires an explicit
selection scope, not ignoring missing queries.

Allow explicit PAF no-hit records or a checksum-bound accounting companion for
producers that omit them. Absence of an alignment row is not proof of absence of
a match. A producer's secondary-hit/output limits constrain the interpretation;
never claim exhaustive competitors or whole-reference uniqueness.

PAF implementation follows the [PAF definition](https://github.com/lh3/miniasm/blob/master/PAF.md)
and [minimap2 output specification](https://github.com/lh3/minimap2/blob/master/minimap2.1):
validate mandatory spans, strand, sequence lengths and typed tags; preserve
MAPQ 0 as low/ambiguous evidence and 255 as unavailable. `tp` is an aligner role,
not a generic SAM supplementary flag. Keep unsupported roles unknown. Use
validated `cg`/`cs` geometry for base/edge claims; reject contradictory tags.
The [minimap2 manual](https://github.com/lh3/minimap2/blob/master/minimap2.1)
documents `-c`, `--cs` and `--paf-no-hit`; these are external preparation options,
not new GENtle subprocess behavior. Mapping-only identity may be approximate.

Missing detail cannot produce exact base counts or splice claims. Strict mode
refuses a requested claim lacking its required tags/annotation. An explicitly
bounded evidence-only mode may retain coarse segments with unavailable fields
and reasons, never the same verdict with invented values. Publish the claim
prerequisite matrix with the schemas before implementation.

## 6. Seven Narrow Implementation Slices

### A. Versioned Contracts And Claim Prerequisites

Add `src/engine/protocol/primer_payload.rs` and export it from the shared protocol
module. Specify query/anchor identity, expected-target crosswalk, reference and
annotation schema, batch accounting, policies/budgets, outcomes and result binding.
Use `deny_unknown_fields`, typed scopes/statuses and explicit units. Document
required versus unavailable fields for each scientific claim. No live capability
is advertised before its complete engine path exists.

Author serialization, invalid-input and historical Phase 1/RNA wire-byte tests.
Ensure a Phase 2 field cannot be smuggled into an old strict schema.

### B. Verified Query Preparation And Input Admission

Implement `PreparePrimerPayloadAlignmentQueries` in a separate engine module.
It verifies report/binding, selected TSV receipts and source row identities.
Prefer existing verified payload FASTA when it satisfies the explicit anchor
selection. For absent FASTA or additional anchors, require explicit opt-in replay
of the exact completed FASTQ units and EOF hashes; never reopen private inputs
implicitly. Reuse existing bounded FASTQ intake, not a new parser or matcher.

Generate opaque query IDs plus a manifest; never parse biological identity from
read names. Anchor IDs bind observed primer/core intervals, direction and primer
alternatives, not incidental array order. Preserve multiple possible primer IDs
rather than treating Phase 1's deterministic primary anchor as unique attribution.

Each query has an explicit primer-free coordinate view. Other observed primer
intervals inside a candidate window must be excluded/masked in the new query
copy and all support calculations, or leave the candidate unassessed with a
reason. Never silently count a second primer as biological payload, change the
original FASTQ, or present these masks as inferred molecule segmentation.
Bind the selected policy and mask transform. Primary-only preparation remains
explicitly primary-only; do not label it all-anchor coverage.

No aligner is invoked. Publish fresh query FASTA/manifest receipts atomically.
If complete evidence for an anchor is unavailable, retain its coverage deficit
instead of fabricating a query or an off-target zero.

### C. Bounded PAF And Reference/Annotation Import

Implement a small PAF parser under
`src/engine/analysis/rna_panel/primer_payload/`; do not add a whole-genome aligner.
Stream/hash raw input bytes, enforce maximum line/tag/operation/query sizes,
check integer overflow, and validate CIGAR/difference-string consumption against
both coordinate spans. Cross-check target identity/length and batch provenance.
No version stripping or alias guesses.

Require explicit query-grouped input for bounded processing, or implement a
separately budgeted temporary spool; do not hide an unlimited in-memory sort.
Preserve all alternatives or fail at a budget, never keep only the first N.
Reference access must be bounded too: use explicitly bound existing indexes or
stream the FASTA while retaining only admitted context ranges. No silent cache
preparation or whole-genome materialization in RAM.

Missing/malformed annotation or incompatible assemblies cannot become empty
gene lists. Keep valid multi-gene annotation overlaps as sets. Require actual
alignment rows for additional segments; a supplementary tag is not another
independently verified alignment.

### D. Read/Anchor Evidence Without A Synthetic Winner

Normalize segments into payload order; retain reference/query coordinates,
strand, identity/edit provenance, aligner role and all equivalent transcripts.
Compute expected-target coverage from unions of aligned query-base blocks,
excluding primer masks and reference skips. Keep aligned support distinct from
exact matching bases. The target-fraction denominator is the number of unmasked
primer-free query bases, not just the aligned subset; retain raw candidate
length, excluded bases and this denominator separately. Zero usable bases stay
unassessed. Insertions/deletions are explicit.
Overlapping HSPs and compatible isoforms never multiply a base count.

Apply declared length, coverage, identity, mapping-quality and equivalence
policies. Different-reference scopes are not pooled or ranked on incomparable
raw aligner scores. Retain best-equivalent alternative segment chains; ambiguous
overlapping mappings are not simultaneous ordered biological segments. Chain
enumeration is budgeted and cannot truncate alternatives into a unique call.

Return expected/off-target dominant, mixed, ambiguous or unresolved evidence
only when its prerequisites hold. Bind thresholds; they are not calibrated
probabilities. Lack of expected-target mapping is not evidence of off-target
dominance. Preserve unassessed/unmapped/ambiguous distinctions.

Find the primer-adjacent mapped edge only when alignment geometry reaches it
under explicit policy. An unaligned query prefix or indel leaves an interval or
unavailable edge, not an extrapolated exact breakpoint. Separately report
internal segment boundaries and annotated exon-edge compatibility. Collapse
splice observations by genomic donor/acceptor/strand, retaining isoform citations.

Compute requested terminal-core seed comparisons against verified, oriented
reference context. Short/missing/ambiguous context stays unavailable. Seed
agreement is not a mispriming probability or proof that the primer acted in vivo.
`cross_locus` requires two sufficiently long, ordered payload segments with
compatible evidence, never a primer plus one foreign payload segment. Keep
technical artifact and biological-candidate observations on independent axes.

### E. Partner Profiles And Separate Edge Clusters

Partner profile keys include primer ID, cohort, sample/barcode/priming,
reference scope and the full partner-equivalence group. Count a source read once
per key even with multiple anchors, isoforms or duplicate alignment rows.
Expose anchor observations separately. Ambiguous group counts are not exclusive
member counts and must not be summed as independent reads or molecules.

Edge clusters additionally bind reference/contig identity, strand, edge type,
partner group and explicit coordinate tolerance. Use a specified deterministic
sorted rule with maximum cluster span bounded by tolerance; avoid transitive
nearest-neighbor chaining that joins distant edges. Retain actual min/max edges,
a documented representative and distinct-read counts. Internal segment-boundary
clusters stay separate from primer-edge clusters.

Export both rankings with eligible, assessed, unmapped, ambiguous, unassessed and
excluded denominators. Full-cohort and selected-sentinel scopes remain separate
even when their leading partners differ. Never silently add transcriptome and
genome observations of the same read as independent support.

### F. Publication, Managed Execution And Thin Adapters

Publish a versioned JSON report, read/anchor TSV, partner-profile TSV and
breakpoint-cluster TSV with artifact hashes, inputs, policy and non-claims.
Keep large rows in sidecars, not an unbounded JSON or agent response. Bounded
inspection reports exactly what it verified; full sidecar verification is an
explicit operation, not a display side effect.

Reuse fresh-path reservation, staging, handle closure, request-bound control and
existing managed execution lanes. A completed unit is a declared alignment
batch whose dependent inputs have all reached EOF/hash validation. Cancellation
or post-start control error retains only such units with explicit missing
coverage. Data/hash/resource/publication failures remain hard errors without a
final report; orphan staging is not evidence. Partial aggregates are labelled
and recomputed from retained complete units, not a half-read PAF.

Add narrow operation handlers, shared-shell routes, registry descriptors and
generic MCP/agent access. Use a named artifact-only execution classification;
do not generalize the proposed lifecycle framework. Writes and cancellation
require confirmation, project state/history/files remain unchanged, and late
host cancellation cannot relabel already published evidence. Agent feedback
must preserve completed/cancelled/control-error outcomes and remain responsive
to inspection/control while processing. No new GUI or OS-shell requirement.

### G. Synthetic Tutorial, Regression Matrix And External Acceptance

Use small hand-crafted fixtures with provenance/recreation instructions and
raw-byte manifests. Prefer arbitrary synthetic IDs over the current panel's
names. Update tutorial 10.04 from its canonical source where generated; document
the explicit external-alignment boundary without fetching data or running an
aligner automatically. Include all mandatory examples from Glen's prompt as
generic relationships, not biological shortcuts.

Required cases:

- Foreign-only versus expected payload; equivalent isoforms, overlapping HSP
  unions, MAPQ 0/255, primary/secondary alternatives and incomplete search scope.
- Different primers sharing one expected gene; same partner at different edges;
  same read with repeated names/anchors; ambiguous primer attribution and masks.
- Separate full-cohort/sentinel rankings, distinct partner/edge denominators,
  two genuine payload loci versus primer-plus-one-partner negative control.
- Positive/negative strands, transcript/genome projections, adjacent bases,
  insertions/deletions/skips, clipped query prefixes, both cs forms, malformed or
  contradictory tags, and isoform-independent genomic splice-edge collapse.
- Missing/foreign/duplicate identities, explicit unmapped versus omitted queries,
  reference/annotation mismatches, changed files, hashes and honestly rehashed
  structural tampering. Partial source units never acquire full-cohort status.
- Every resource cap, cooperative cancellation/control faults before and after
  verified units, existing outputs, traversal/symlink/nonregular inputs, Windows
  handle closure, spaced/apostrophe paths, LF/CRLF and input immutability.
- Old Phase 1/RNA wire/artifact compatibility; equivalent CLI/MCP/agent payloads,
  confirmation refusal, stale workers, late Cancel and read-only inspection.
- Renamed/reordered IDs, panel cardinalities and primer lengths outside the
  motivating examples; counts must remain covariant, not name-dependent.

## 7. Verification And Completion Gates

During implementation, run scoped rustfmt and Python/source/document checks
only. No local Rust build/test/check, push, CI dispatch, private-file replay or
live-model invocation is authorized by this plan. Keep the existing goal paused.
Each slice should be a narrow commit with tests and an honest execution boundary;
do not regenerate unrelated catalogs or repair unrelated pending features.

Glen/CI must execute focused new `primer_payload_` tests, old `primer_trace` and
RNA regression filters, protocol/parity/glossary and actual CLI/MCP replays,
all-features `--tests`, headless and locked checks. Record exact source,
toolchain/platform/profile, command counts, primary failures versus cascades,
skips and artifact locations. Canonical generators, LF/CRLF, Windows paths and
handle behavior need actual execution; source checks are not substitutes.

A live inner-agent run separately checks discovery, reviewed request creation,
write approval, progress, cancel/inspection responsiveness and explanation of
ambiguous/partial outcomes. It is not established by a deterministic mock.

Completion requires passing synthetic engine and process-level adapter paths,
documented schema/coordinate/count contracts, fresh required projections and
retained evidence for the claimed platforms. Private cohort accuracy/throughput
and scientific interpretation remain separate external gates. No new release
candidate or recurring bundle/SHA synchronization is required during development.

## 8. First Action On Implementation Approval

Implement slice A and the smallest end-to-end synthetic query-preparation case
from B first. Resolve the claim-prerequisite matrix and complete query accounting
before building partner rankings. This prevents an attractive summary from
depending on unverified source rows, incomplete query populations or hard-wired
assay labels. Do not begin with the private current panel.
