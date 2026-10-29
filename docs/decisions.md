# GENtle Decisions

## Repository Line Endings

Owner decision, 2026-10-21: source and ordinary repository text use Unix LF
line endings in both Git and working copies on every platform. The root
`.gitattributes` default `* text=auto eol=lf` overrides ambient Windows
`core.autocrlf=true`; `.editorconfig` requests the same editor behavior without
changing contributors' global Git settings. Binary files remain byte-exact.
Existing explicit `-text` rules for immutable audit originals also remain
byte-exact, including their historical mixed line endings; do not normalize
them or rewrite their recorded hashes. Scoped LF rules for hash-bound inputs
and generated outputs remain explicit documentation of those byte contracts.

`scripts.test_tutorial_checkouts` exercises real LF/CRLF-policy Git checkouts
with the committed attributes, source text, binaries and verbatim records.
Missing-scoped-rule negative controls also remove the general LF default so
they still reproduce the original unprotected CRLF boundary.

## Observed Junction Structures: S0 Technical Boundary

Owner-authorized 2026-10-27: implement the pure Rust S0 kernel from the
[reviewed proposal](observed_transcript_structures_plan.md), not the later
operation, adapters, file publication, SAM/BAM ingress or graphical slices.
`reference_gene` is mandatory and source-bound by stable gene ID. No gene,
species, assembly, strand, canonical transcript or protein label is a default.
Reuse the payload annotation shape and explicit source-specific gene bindings;
multiple providers do not require multiple alignment searches.

The `observed_junction_chain` digest identifies the declared genomic reference,
assembly, contig and exact ascending skips, not a full RNA or a strand call.
Selected gene, model strand, endpoints and support counts cannot alter this
physical key. Preserve sourced model strands, query-to-reference orientation,
payload-end constraints, masks and alignment receipts separately. Unspliced
observations retain record identity but receive no structural-family ID.
Known-junction recombinations, unannotated skips, coverage divergence and
unassessed detail remain distinct. Comparison is exact; no boundary snapping.

Factor the operation walker/projection and best-equivalent query DAG once for
Phase 2 and S0. Preserve Phase 2 wire bytes, decisions and work accounting.
S0 must use legal paths, not the union of optimal nodes, and must not trust
whole-span Phase 2 gene labels as transcript compatibility. Reconcile anchors
by original source record before support counting; retain conflicting/partial
views, incomplete producer declarations and unknown primer/flank qualification.
Neither recurrence nor a chain implies independent molecules, natural RNA ends,
retention mechanism, fusion, expression or protein identity.

The exported S0 Rust API remains experimental and pure. It checks declared
hashes/geometry in memory, not corresponding files or aligner search. The
separately owner-authorized S1 boundary below owns verified ingress/publication;
main's independent lifecycle Slice 3 does not implicitly cover new operations.
Native/platform replay and scientific acceptance remain source-bound. No assay
thresholds are approved.

### Observed Structures S2.0b Private Boundary

Owner-authorized 2026-10-28 after S2.0a's commit/rebase: use one private
format-neutral observation/path/comparison core with typed identity policies.
Project legacy results back to the unchanged v1 DTOs, retaining full input
serialization, digest tuples, charged work, sort keys, errors and stored-support
assertions. Never update the frozen preservation expectations to bless a refactor.

Decoded SAM uses declared query-unit keys, not QNAME-derived original records;
separate private context/observation domains share only the physical chain
recipe. Normalize clips/padding before the unchanged walker. Unknown alignment
identity stays absent, not a synthetic PAF statistic. Unsupported units supply
zero interpreted segments; reference/dictionary/search declarations cannot
become verification, and unknown masks remain independent from search coverage.
No SAM text/file parser, public schema, operation, annotation preparation,
publication, GUI or biological acceptance is authorized by this private slice.
See the [implementation and verification scope](observed_transcript_structures_s2_plan.md#s20b-implementation-2026-10-28).

### Observed Structures S2.1a Pure Preparation Boundary

Owner-authorized 2026-10-28: implement the reconciled strict GTF subset as a
private bounded incremental consumer, not annotation file/JSON ingress,
publication, operations or adapters. Keep gene/source/reference selectors and
gap/transcript-row policies explicit; no biology or provider defaults. Hash
supplied raw bytes without claiming original-file or provider verification.
Preparation failures/absence have typed outcomes and no annotation, never an
empty v1 annotation used as unavailable evidence. Preserve scoped diagnostics,
unknown strand, versions and excerpt/completeness limits.

Retain the existing annotation geometry, legacy enum/DTO meanings and all
Phase 2/S1 byte/work guards. Touching exons are refused only by preparation;
do not tighten legacy validators to enforce a receipt they do not consume.
Public dictionary/receipt contracts, safe publication and named lifecycle/
adapter migration belong to separately authorized S2.1b/c. See the
[implementation ledger](observed_transcript_structures_s2_plan.md#s21a-implementation-2026-10-28).

### Observed Structures S2.1b File Boundary

Owner-authorized after the S2.1a commit on 2026-10-28: expose strict versioned
preparation request/receipt/dictionary records and one direct engine method.
Reuse the canonical annotation geometry and the same typed preparation
vocabulary, not adapter-local biology or a new transcript model. No Operation,
CLI/MCP/GUI route or implicit lifecycle migration belongs to this slice.

Hash the raw supplied GTF stream that is parsed, bind every dictionary byte and
retain source/gap/version policies in the request/receipt digest. Keep EOF/data
verification separate from execution status and complete gene analysis. Partial
input never becomes absence or measured zero. Refusals publish receipt only;
completed geometry uses raw-snapshot release and remains structurally readable
by legacy validators without claiming they enforce the new receipt.

Use preparation-only limits and safe fresh-directory publication, with the
receipt as last completion marker. External writes remain separate from the
unchanged project/journal. Validate raw paths before normalization and existing
identities before publication. Never recursively remove potentially unowned
output after failure; preserve/disclose orphaned incomplete directories.
No data/hash/resource error is converted into a recoverable biological outcome.
The [protocol](protocol.md#s21b-annotation-preparation-direct-engine-only) records
the new digest recipe and conservative JSON/memory admission model; neither
changes existing S1/Phase 2 bytes, limits or errors. Native/platform/provider/
scientific acceptance and S2.1c named adapters remain separately gated.

### Observed Structures S2.1c Operation Boundary

Owner-requested continuation after the S2.1b local commit on 2026-10-28:
register only named annotation preparation and request-bound cancellation.
Keep one strict file/biology boundary. New operation-derived approval covers
every external receipt/marker write attempt, independently of model/candidate
flags and generic MCP confirmation. Explicitly migrate constructor/no-save,
no-checkpoint and managed policies; do not infer them from `ReadOnly` or output.

Typed standalone admission separates preparation execution and marker control
from the live project. Retain correlated artifact metadata without merging a
private journal or vetoing an already-published result on late host cancel.
Raw-text admission and workflows retain their conservative project-bound
host policies. Compact transport retains a separate gene outcome: complete
source refusal is unavailable evidence, not zero biology or successful models.
Do not re-query sources to validate a worker result; correlate request/receipt
digests, serialized receipt bytes and its resolved publication path.

The new optional result/worker-binding fields are omitted for old operations.
Phase 2/S1 wire bytes, existing schemas, helper limits and frozen assertions
remain unchanged. SAM ingress, S2.3 receipt consumption, GUI, provider/gap
policy acceptance and native/scientific verdicts remain separate. See the
[S2.1c ledger](observed_transcript_structures_s2_plan.md#s21c-operation-boundary-2026-10-28)
for authored versus executed scope.

### Observed Structures S2.2a Private Text Boundary

Owner-authorized continuation, 2026-10-28, from the reviewed plan: decode only
caller-supplied raw SAM lines in a private, bounded incremental consumer.
Keep raw source hashing/framing separate from actual file verification and
retain one row per supplied alignment ordinal. Equal/unavailable QNAME and
RG/SA declarations do not establish original read identity, unit membership,
search completeness, masks or RNA strand. Never synthesize alignments from SA.

Reuse decoded-SAM geometry/target checks without changing old inputs, public
DTOs, error strings, charged work, identities or Phase 2/S1 wire bytes. Validate
but discard sequence/quality strings, unused tags and program/comment text.
Return no successful prefix after data/resource/control failure. Explicit
budgets include temporary normalization and all rows/optional arrays, even
ones that cannot later contribute support. This is a documented subset, not a
claim of full SAM validation, provider verification or biological identity.

Verified file ingress, explicit row-to-unit/mask witnesses, annotation/receipt
joins, saved SAM schemas, operations/adapters and GUI remain later work. The
[S2.2a ledger](observed_transcript_structures_s2_plan.md#s22a-private-text-decoder-2026-10-28)
records execution separately from authored guards and the prior owner report.

### Observed Structures S1 Intake And Lifecycle

Owner-authorized 2026-10-27: consume one explicitly selected, bound genomic
Phase 2 batch through the shared engine, parser/CLI and generic MCP. Keep
`reference_gene` configurable by source/stable ID, preserve all old Phase 2
schemas and bytes, and publish a separate versioned result. Source trace,
prepared-query witnesses, full evidence accounting/hashes and independently
bound annotations are verified; original FASTQ, PAF and reference FASTA are not
reopened. Recorded producer receipts are not an independently repeated search.

Use annotation-bound target lengths. Unknown competitor lengths make the whole
query view unavailable, not a selectively filtered mapping. Excluded/missing
anchor views also suppress definite support from other anchors of that source
record. Source statuses/cohort gaps, per-skip source comparison and nonadditive
sample/barcode/primer support remain visible. Annotation uncertainty is not zero.

Summarize/inspect/cancel join the named request-aware preserve-content selector:
no constructor repair or implicit CLI/MCP project save. Raw `apply` retains the
normal execution audit without checkpoints or project mutation. Direct pure
analysis methods and managed standalone publication do not change the live
project/journal; the worker's private audit is not merged. Workflow publication
remains the existing conservative project-bound path, not a new artifact-only
workflow service. Explicit report/marker writes require approval independently
of this no-save fact; generic MCP confirmation stays mandatory.

The selected batch is an atomic analysis unit. Cancel/fault before completion
leaves analysis/counts absent, not zero; a complete unit can survive a final
cancel/fault with that outcome retained. Data/hash/resource errors publish no
final report. Bound writers stage fresh output and publish the report last.
No source schema rewrite, BAM/SAM ingress, aligner, GUI, protein assignment,
new dependency, push or scientific acceptance belongs to S1. See the
[contract and pending execution scope](observed_transcript_structures.md).

## Primer-Anchored Read Trace

Phase 2 technical implementation authorized 2026-10-26: keep a separate
[primer-free payload evidence contract](primer_payload_evidence.md). Source/query,
reference/annotation, producer accounting and expected-target mappings are
explicit inputs, never current-panel rules. Preserve primer ambiguity, mask
other observed primers in new query copies, separate partner/edge counts and
retain only verified complete units on cancellation. Named managed artifact
publication does not commit project history. No implicit original FASTQ replay,
aligner, biological calibration or mechanism approval. Runtime acceptance is
still external; this does not upgrade Phase 1 scientific status below.

Proposed 2026-10-26, pending owner confirmation of biological axis/default policy.
The owner authorized implementation of the revised technical Phase 0/1, not
kit identity, assay validation or calibrated thresholds. See the
[revised plan](primer_trace_phase01_plan_20261026.md).

A separate `TracePrimerReads` operation observes explicitly supplied oligos
and optionally checks core plus boundary-adjacent payload on the same mature
template interval. It does not alter RNA-panel nomination, legacy classes,
original reads or their orientation. Oligos are component-defined, not assumed
44-mers. Declared primer direction is independent of chromosome/read strand
and does not establish whether the oligo acted in RT or PCR.

Keep independent primer, position, multiplicity, locus-relation and artifact
axes; reserve junction/segment axes without emitting placeholders. Handle
alone never attributes a target. A seed only nominates a candidate; retain
all validated template alternatives. Boundary mismatch is an observation,
not a mechanism assertion, and missing/short/ambiguous evidence is unresolved.
Internal primers are segmentation candidates for a later phase, not inferred
junctions. A candidate payload containing another anchor remains unsegmented
and unresolved in Phase 1. Do not emit mechanism-confirmation labels including
"trans-splicing", "fusion" or "chimera confirmed" as result classes.

Policy/budgets, exact references and input bytes are content-bound. Independent
file effects require approval. Standalone managed trace publication is a
named artifact-only exception with its own result binding, never a generic
ReadOnly/lifecycle inference. Validation/show have bounded non-project work;
marker control is request-bound. Generic MCP `op` retains confirmation but
does not rewrite its state file for these named operations. No general lifecycle
proposal is implemented by this exception.

Reuse the existing strict FASTQ/hash/control primitives. Quality is retained
for this consumer without changing existing nomination calculations or output
schemas. Only EOF-verified completed units contribute; cancelled/control-error
outcomes keep missing inputs explicit and data/resource failures fail closed.
No private-data fixtures, new GUI, segmentation graph, aligner adapter, molecule
count, expression inference or scientific acceptance follows. The configurable
handle/poly(T) observation does not authenticate an ONT kit. Payload FASTA
remains explicitly opt-in. Real protocol roles and numerical defaults require
owner/Glen review before this decision becomes biologically accepted.

## Independent-Locus FASTQ Panel Nominations

Owner-authorized 2026-10-22: a portable RNA panel binds independent loaded
loci, source-record and mature-sequence hashes, declared assembly/release/gene
labels, local exon geometry and namespaced `(locus_id, transcript_feature_id)`
members. Labels and content hashes are not authenticated reference provenance.
Potential overlapping source views are disclosed and refused for screening,
not double-counted or silently turned into a consensus annotation.

One packed-word index retains all template postings. Exact seeds are chained
with the existing Rust-Bio LCSk++ implementation. Every qualifying member and
orientation is retained; partial seed co-nominations are not whole-read-compatible
transcript sets or exclusive isoform calls,
expression, full-length capture or whole-transcriptome uniqueness. Budgets fail
explicitly, never by annotation-order truncation. Junction support requires
declared flanks; substitutions/indels require separate alignment inspection.
The legacy single-locus winner/tie semantics are not a panel assignment oracle.

Parse four-line FASTQ directly, including concatenated gzip members. Compute
compressed-byte SHA-256 and supplier MD5 in that same pass, verify expected
digests at EOF and reject changing/incomplete copies. No upfront full-file
checksum pass or full FASTA conversion. MD5 checks transfer integrity only.
Publish a final report only after artifact handles close; failure is distinct
from zero support. Cancellation retains verified completed-file units, marks
the run cancelled and discards the partial unit. Continue with remaining
explicit files in a fresh run; there is no in-file or changed-input resume.

Demo preparation, build, screen, bounded inspection and cancellation are shared
engine commands available to the inner agent and other adapters. Artifact
writes require explicit approval; no recursive agent execution, automatic
trimming, oligo(dT) flip, implicit resource preparation, POD5 or scheduler is
introduced. A bounded TSV inspection exposes that it has not reverified the
full sidecar. See [tutorial 10.04](tutorial/10-04_targeted_nanopore_fastq_panel.md).
Rust/native/real-library calibration and throughput acceptance remain external.

Owner-authorized 2026-10-23: report-only count exports separate single-candidate
panel support, overlapping per-member ambiguous support and exact member-set
counts. Each read contributes once to its set, never to multiple exclusive
isoform counts. No fractional assignment, original molecule/UMI count or
expression inference follows from enrichment. Use the identical panel/full
policy, preserve sample/barcode/priming bindings and refuse repeated compressed
input hashes across reports. Different hashes do not prove read/molecule
independence. Only completed verified units contribute; cancelled sources
remain partial and samples without such units remain numerically unassessed.
Export revalidates retained report/panel content, not original reads or large
nomination sidecars. This shared, confirmed external export and subsequent
read-only inspection remain available to the inner agent without OS tools.

The narrowed 2026-10-23 repair preserves `(class, sorted members)` jointly in
v2 receipts and summaries. Never reconstruct that association for member-only
v1 reports: label it `legacy_class_unresolved`. A distinguishing junction in A
does not exclude B nominated through shared exons; there is no whole-read
compatibility filter or isoform resolver in this slice. Retain named requested
but uncompleted inputs and numerically unassessed samples. Do not merge a
delivery FASTQ with its constituent chunks: different hashes cannot detect
shared reads. Keep fail-closed resource limits, request-bound cancellation
(including matching pre-start requests), throttled marker I/O and honest
compressed-byte/read-ahead progress. Table verification is explicit, bounded
and digest-bound, not an implicit expensive display operation or a revalidation
of original reads/nomination sidecars. V1 fields remain readable without
changing their semantic digest or promoting their historical names into claims.

Owner-authorized 2026-10-23 follow-up: a marker error after startup preserves
checksum-verified file units under an explicit `control_error` outcome with
a source-bound diagnostic, never as an authorized cancellation or successful
analysis. Discard the unfinished unit, retain named missing inputs, and keep
control-error source summaries partial even if all requested files happened
to verify before finalization. Preflight marker errors and unrelated data,
resource and publication failures keep their existing hard-error behavior.
No automatic retry, overwrite, reconciliation or budget-policy change follows.

Owner-authorized 2026-10-24: native and command-line panel preparation share the
source inventory and builder. Default to the full annotated transcript group;
changing a competitor set is scientific input, not a display filter. Extension
binds an explicit saved base and verifies every original loaded source hash,
rejecting missing/changed originals and duplicate locus identities. Publish a
new panel only; edit/removal is explicit full recomposition. Old panel/report
files are never overwritten or reinterpreted, and different panel digests remain
incompatible summary inputs. GUI requests/drafts are session-local; source
provenance remains caller-reviewed, not authenticated by editable labels.
Reuse stored gene-set resolutions for target membership, with exact report
hashes and annotated stable gene IDs for each source. No name/alias inference
or parallel panel gene-list primitive is introduced. A resolved expansion is a
new operand and a new panel, not an in-place mutation of historical evidence.
Explicit `register_in_project` retains the panel snapshot and its set/source/
parent relationships for the GUI graph and shared SVG projection. Registration
is undoable project state; the external panel file is not. Graph nodes describe
reference definitions, never read counts or proof of isoform discrimination.

Owner-authorized 2026-10-24 follow-up: managed standalone `ScreenRnaReadPanel`
and `ExportRnaReadPanelSummary` are a named artifact-only allowlist, not an
inference from checkpoint kind, `state_changed=false`, or arbitrary output JSON.
They retain bound published results without merging their detached journals into
the live project. Other operations, workflows and registered panel builds keep
their normal optimistic mutation commit/history rules. The synchronous engine
API still returns its ordinary operation receipt; this does not authorize a
host to append an artifact worker's private journal to a concurrently edited
project. Later graph attachment is a separate short explicit undoable operation.

Submitted and effective requests are independently digest-bound: manifest MD5
resolution/normalization is not tampering. A retained artifact receipt also binds
panel identity, report content, location and scientific outcome. Optional reviewed
panel digests are checked before indexing/streaming/publication. New optional
request/report fields are omitted when absent so legacy semantic hashes and
cancel identities do not change. Hashes bind content, not authenticated origin.
Never reopen a presumed result solely because a report filename exists.

After engine publication, host cancellation cannot relabel completed/cancelled/
control-error evidence. Engine cancellation before publication is unchanged.
Cancellation requests and actual scientific outcomes remain distinct. A control
error is failed screening with verified units retained, not authorized cancel.
Metadata observation and cancellation use independent bounded host lanes, never
ordinary execution slots or a live project lock; expensive table verification
still consumes bounded execution capacity. Process-local receipts/job IDs are
not cross-process worker handles or restart recovery. Reopening published bound
evidence and request-bound marker cancellation remain shared headless operations.

Owner-authorized 2026-10-25 clarification: draft `--request @FILE` loading is
bounded observation work; admission captures a path, not file bytes or a live
project snapshot. Keep source I/O off the GUI thread. Model Stop, explicit
same-process command Cancel and request-bound marker control are distinct.
Stopping a model or replacing a project does not implicitly cancel an
artifact-only screen; ordinary project-bound work retains its owner-change
cancellation. Marker failures never trigger an implicit service-cancel fallback.

Owner-authorized historical-review follow-up, 2026-10-25: new locus records
add an optional uppercase-forward-DNA hash independently of the DNA-plus-feature
record hash. Equal present hashes flag potential duplicate sources even with
different source IDs and absent anchors, never biological identity. Historical
absence remains unknown and omitted, preserving old digests. Keep all supported
members/orientations and restricted-competitor warnings in reports. Resource
failures remain hard failures, with no new over-budget biological class or
implicit change to the competitor universe. A locus-position index requires
measurement and per-template candidate-equivalence evidence before replacing
current chaining; the best genomic chain cannot silently discard alternatives.

Owner-authorized duplicate-source/ambiguity follow-up, 2026-10-25: keep physical
source admission separate from exact mature-sequence sharing. Unequal extracts
with a repeated declared assembly/gene identity remain refused. Sequence-equal
selected transcripts across otherwise independent loci become a digest-bound
`shared_transcript_sequences` diagnostic, not inferred gene identity or an
automatic rejection/merge. Missing source relations require human provenance
review; labels, matching RNA or a hash do not resolve them. Preserve every
qualifying member and count the cross-gene co-nomination set once. Historical
diagnostic absence stays omitted, not an empty-check verdict. No changes to
legacy winner selection, posting caps or thresholds follow. Legacy regressions
compare complete reports before/after file-only panel building and equivalent
no-op expansion, excluding only creation time and explicit mode notices.

RNA residual follow-up, 2026-10-26: optional deterministic work profiles observe
the existing nomination path and verified file units, never select candidates
or change budgets. Keep profiles out of panel/request identity, omit historical
absence and keep timing outside scientific records. Carry selected sequence-
sharing evidence into summary JSON/TSV without inferring physical source or
causality: all nominated members may share one mature sequence while their
biological origin remains unresolved. Historical absence stays unknown, not
negative. Synthetic substitution/indel tests are software characterizations,
not calibrated Nanopore sensitivity. Zero-length FASTQ remains a hard error;
any future different intake class requires a separate owner decision.

RNA source-overlap correction, 2026-10-26: new panels explicitly bind
`source_overlap_assessment.rule=genomic_exons_v1`. Decide coordinate duplication
from selected exons mapped into the same declared assembly/chromosome and the
composition of record-anchor and template strands, not record orientation or
flanking/intronic extents. Same-strand exonic overlap remains a conservative
refusal. Opposite-strand overlap is retained annotated antisense evidence,
never read orientation, a gene-identity verdict or a reason to collapse counts.
Unknown/unmappable orientation or geometry retains disclosed conservative
whole-record admission; invalid bounds cannot prove disjointness. Same loaded
record, equal present forward-DNA hashes and repeated declared assembly/gene
identity remain independent safeguards. Historical absent assessments retain
the exact former rule and hash. Rebuilding/upgrading evidence is explicit and
changes panel identity. Keep diagnostics engine-owned in panels and summary
JSON for every adapter; summary inspection does not authenticate or freshly
re-map the source geometry. No nomination, budget or intake-policy change follows.

## Development And Release Verification Cadence

Owner decision, 2026-10-15: continue development on `main` without intermediate
candidate freezes, recurring bundle handoffs or waits for synchronized acceptance
at every commit. Run normal focused checks and available CI independently;
record each result's actual revision/toolchain/profile as evidence metadata,
not as a recurring coordination request. Fix observed regressions without
requiring the full acceptance campaign before the next unrelated change.
Consolidate the still-open cross-platform/native/agent/scroll/package gates at
the actual release preparation, using one selected release candidate and
rechecking affected gates if that candidate changes. Earlier passes never
silently certify newer code. No pending gate is waived or marked complete.
Existing no-local-Rust, no-push-without-request, authentication and publication
restrictions remain unchanged. Do not issue a new bundle or SHA-sync request
for routine development; an exceptional transfer requires an owner request.

## Workspace Compiler Compatibility

The root `Cargo.toml` owns workspace-wide edition, resolver, any declared Rust
minimum and shared dependency policy. First-party workspace members inherit
the edition; any future member Rust minimum must inherit `workspace.package`
rather than define a separate floor. No explicit Rust minimum is selected now.
Private component versions can remain independent. Third-party dependency
requirements and vendored upstream metadata are not compatibility promises for
the whole application and must not be erased to force dependency updates.

A successful build establishes compatibility only with its recorded toolchain,
features, dependencies and target; moving `stable` CI does not establish an older
minimum. The private sequence-design core's standalone extraction/publication
gate needs an explicitly verified compiler contract, not restoration of the
unverified Rust 1.85 declaration. A future pinned older-toolchain job building
and testing that dependency-free core alone would establish core-only evidence.
It cannot by itself justify a workspace-wide `workspace.package.rust-version`:
that broader declaration also needs matching workspace/dependency, feature and
target verification. Selecting a number and adding those jobs are deferred;
no new blanket publication gate for other members follows from this finding.

## Explicit Variant Allele Selection

`MaterializeVariantAllele` is the shared engine boundary for single-base
allele materialization. A multiallelic marker never authorizes an arbitrary
alternate: GUI, shell and agents must provide a reviewed reported A/C/G/T
choice, distinct from the reference. Reference materialization needs no
alternate; invalid/unlisted or contradictory selections fail before sequence
creation. Exact selected bases remain in output qualifiers and handoffs.

Alleles refer to the loaded sequence's forward orientation. dbSNP position,
reference and alternatives stay bound to the selected genomic placement;
gene strand or clinical nomenclature does not reinterpret them. A base/ref
mismatch is a refusal, not permission to guess a complement. Reporter preview
geometry, insert orientation and topology remain separate from SNP selection.
Historical SVG presence and synthetic backbones never establish assay validity.

New dbSNP markers retain separate reference and assembly-family checks. An
assembly fallback or unspecified family remains unavailable even when the
extracted base matches; this is inspection evidence, not a liftover or authority
to materialize alleles. Unavailable/mismatching checks block materialization.
Only unambiguous single-base placements receive `vcf_ref`/`vcf_alt`; every
retained SPDI allele must have canonical single-base deleted/inserted strings,
including alleles excluded from the normalized same-reference alternate list.
Non-SNV alleles stay raw SPDI with their original coordinates and empty strings,
not invented VCF anchors. Legacy markers without these checks still require the
existing fresh reference/base validation and establish no assembly-check claim.
The GUI validates both inserts through the shared engine helper before applying
either. The two applications share one write lock and the existing rollback
guard; a later refusal restores project/journal/undo/redo state before GUI
success feedback or output-name changes. Successful pairs retain the two
canonical operation/history entries; rollback keeps revision invalidation and
does not reuse consumed operation identities. Read-only handoffs emit an
explicit alternate only when the loaded source passes that same validator
with one alternate; otherwise the command requires review. No candidate report
silently supplies an allele or assay proof.

## Splicing UniProt Presentation Order

Status: active for the Splicing Expert canvas, transcript/exon matrix and SVG.

`SplicingExpertView::uniprot_presentation_layout()` is the sole display-order
policy. It joins lanes to saved matrix rows by exact feature ID, never ordinal.
Referenced, other evaluated and unevaluated groups sort by target first,
normalized transcript ID, feature ID and original lane index. The normalization
is for sorting only; existing engine-owned exact-xref classification is unchanged.
Headers appear only if a displayed lane was evaluated. Badges count sources
with actual matched transcript xrefs, in accession order, with review status,
overflow and diagnostic indicators; full evidence stays in JSON/hovers. `?`
indicates a loaded-evidence diagnostic, not necessarily transcript ambiguity.

Missing or duplicate joins are unevaluated, with diagnostics; unused matrix
records remain in original order at the end. Missing cells are unknown, not
absent exons. One shared half-open row-coordinate table drives GUI painting and
hit-testing; headers cannot select adjacent transcripts. Boundary summaries
deduplicate oriented presentation text only, retaining every source marker.
GUI caches summary labels/hovers and motif-detail rows with the existing
validated content identity; SVG consumes the same readable source-note wording.
For views without a valid fingerprint, the fallback cache identity includes
every ordered boundary record, including source annotations and pair flags;
changing only a boundary must invalidate cached summaries and motif details.
Boundary and paired-signature exception flags share one presentation predicate.

Grouping never reorders/mutates engine payloads or changes their fingerprints,
sequence derivation, biological coordinates or evidence. Other isoform-evidence
and RNA-read views deliberately retain their existing payload order. SVG labels
use deterministic collision suppression at readable sizes, with all suppressed
facts retained in shared summaries, hovers and structured records. No renderer
performs a UniProt fetch or claims transcript validity/global absence.

## Explicit Gene-Set Storage Projection

A logical gene set does not imply physical samples or placement permission.
`CreateGeneSetArrangement` requires a complete persisted resolution, one
explicit exclusive-singleton binding per distinct member and a separate true
physical-placement flag. The existing outer apply owns one serial arrangement,
its default rack and content-bound source/projection metadata; do not nest
operations or create a pool/aliquot implicitly. The minimal slice is 1..24 tubes.
Source-set/member order and sequence snapshots remain inspectable; later rack
edits do not rewrite the creation evidence. Digital placement is neither a
verified freezer inventory nor qualified fabrication. Generic reviewed `op`
routes supply parity; no GUI-only mapping or new automatic planner is implied.

## Read-Only Genome-Edit Intent

SNV intent validation is separate from editor/guide design, synthetic synonymous
redesign and experimental approval. `EvaluateGenomeEditIntent` is one shared
read-only engine path with strict reference/anchor/digest binding; a caller's
serialized projection status cannot certify those facts. Genomic and synthetic
identities remain disjoint. Retained unversioned names do not authorize guessed
accessions, inline reference association is caller-declared, and absent
annotations stay unavailable. Desired sequence strings in a transient report
are neither materialized project DNA nor observed alleles. For this operation,
in-process apply skips persistent reconciliation, and the named CLI/MCP direct
operation routes suppress their implicit state-file save. MCP's independent
tool confirmation does not authorize a state-file write for this inspection.
The owner-authorized [Lifecycle Slice 3](operation_lifecycle_slice3_plan.md)
extends this to decoded raw CLI/MCP and generic JS/Lua Op/Workflow construction:
Plan/S1 and nonempty all-covered workflows preserve post-loader/supplied state
and suppress implicit project saves. Ordinary `from_state`, empty/mixed
workflows and opaque wrappers remain legacy. Unchanged disk bytes and unchanged
decoded state are distinct guarantees; candidate hydration/warnings still occur.
Future
editor feasibility, persisted handoff/apply and laboratory review need their
own explicit contracts, not a widened S1 or transferred sequence-design approval.

S1's no-reconcile and no-state-file-save guarantees are explicit operation
policy, not implied by every registry `CapabilityMutation::ReadOnly` label.
That class also includes file writers such as `SaveFile`; legacy previews can
reconcile state, and `PlanDnaSequenceDesign` has an optional report path. Do not
replace confirmation, history, reconciliation or persistence rules with that
label alone. Before S2, separately review a shared lifecycle policy that
distinguishes project effects, external writes, reconciliation and confirmation;
preserve existing adapter gates. Slice 3's authored regressions await external
Rust/VM/native execution. The separately authorized Plan-only Slice 4 writer
below has its own pending verification; S2 remains separate.

## Degraded Candidate-Store Persistence

Owner-authorized [Lifecycle Slice 6](operation_lifecycle_slice6_plan.md): a
failed candidate sidecar load is unavailable evidence, not an empty collection.
Permissive loading retains the original reference and a structured unresolved
diagnostic; strict loading still refuses. Historical warning-only states remain
unresolved. Only full hydration of an original reference on reopening clears
that diagnostic; no inferred reference, silent repair or discard is allowed.

The shared project saver refuses before any save-owned filesystem change while
that dependency is unresolved, including Save As and implicit saves. Candidate
reads/writes also refuse; unrelated inspection stays available. Project facts
carry unavailable families so missing candidate evidence cannot prove absence.
Healthy empty-store cleanup remains intentional. This is not rollback of prior
in-memory commands/external writes, a hostile-parent/crash-durability guarantee,
or authorization for exporter/wrapper/queue migration. Runtime/native checks
remain source-bound and external under the current no-local-Rust restriction.

## Sequence-Design Preview Report Approval

Owner-authorized 2026-10-26: migrate only the conversation-approval facet of
`PlanDnaSequenceDesign`. The engine-owned typed selector adds approval for every
present report path, including empty/malformed paths; an absent path adds no
operation-derived requirement. Unmigrated operations retain legacy policy.
Explicit candidate and MCP/tool gates remain additive. The shared compatibility
flag `is_state_mutating()` includes external file effects, not just project
mutation. Static Plan metadata is conservatively External with confirmation,
consistent with TSS SVG export; it is not the runtime permission oracle.

Approval covers the exact write attempt, including a cancelled search report.
Successful file publication is not scientific/search completion; retain actual
status and report bytes. No constructor, reconciliation, checkpoint/journal,
autosave, alias/collision or overwrite semantics change follows from that approval
facet alone. The separately owner-authorized
[Slice 3](operation_lifecycle_slice3_plan.md) now supplies covered construction/
state preservation, with external verification pending. The Plan-only Slice 4
destination contract below is separately owner-authorized and awaits its own
verification; S2 remains gated. Focused macOS approval boundary and
parity checks have a [source/feature-bound execution record](operation_lifecycle_policy_proposal.md#executed-local-verification-2026-10-26).
They do not establish full-policy security, native or scientific acceptance.

## Sequence-Design Preview Destination Safety

Only Plan's optional report writer uses this policy; other exporters and project
saving retain their existing contracts. Every present path still needs the exact
write-attempt approval. Reports require a fresh filename in an existing parent:
no truncation, implicit directories, overwrite flag, delete/retry or new name
without review. Identical existing reports, dangling links and special files are
not fresh destinations. Stage the existing JSON bytes, sync a write-capable
temporary handle and publish with no-clobber. Preserve original publication errors.

A nonserialized engine `host_project_path` binds CLI/MCP defaults/arguments and
GUI current-project semantics, including absent files. Reuse ProjectState's
derived sidecar root to exclude the whole save-owned tree. Do not retain loader-
lost metadata or infer filenames in stateless JS/Lua. Python remains CLI-backed.
Explicit history-free clones and imported checkpoints preserve the binding;
guarded commit keeps a newer live Save As binding. The field changes no state,
history, dirty revision, command/request/report digest or approval policy.

Use native directory identity, with conservative ASCII ignore-case for protected
leaf names. Ambiguous same-parent non-ASCII comparisons fail closed rather than
claiming universal Unicode normalization. Native alias resolution does not
defeat adversarial parent swaps. No universal atomicity, crash durability or
perfect temporary-link cleanup is promised. A later explicit project save may
replace a report; an already published report may survive a later command failure,
cancel or stale commit. Scientific completion and file publication stay separate.
Authored regressions and [external checks](glen_next_verification_handoff.md#lifecycle-slice-4-2026-10-27)
still require Glen/CI execution; this is not native acceptance or S2 authorization.

## Synthetic Sequence-Design Boundaries

Status: active for experimental synonymous motif removal and explicit GC constraints;
Glen independently executed both strategies with motif/CDS-wide/window-GC checks
at `de9df8af0a56923df01a7492f61005b1fe079d8c`, scoped `prototype_pass`.
See the [source-bound acceptance record](sequence_design_acceptance.md#executed-audit-at-de9df8af).
This changes no invariant or release, usefulness, biology or performance claim.

Only explicit synthetic coding-insert requests authorize preview, never labels
or natural assay templates. The pure independently versioned core has no root,
GUI, serialization, Python, file or network dependency. Root context resolution
binds explicit codon/residue mapping rather than undeclared asset ordering.
Its private component version does not follow GENtle's application-version
rollover; only workspace-version-inheriting lockfile entries follow that bump.
It stays private and in-workspace until its standalone compiler contract is
verified under [workspace compiler policy](#workspace-compiler-compatibility).
Glen's Rust 1.96.1 audit and CI's moving `stable` selection do not select an MSRV.
Frozen bases, exact protein/start/stop and full-sequence motif evaluation are
invariants. Bounded-search failure is unresolved, not infeasible; cancellation
cannot publish an applicable proposal, including a valid intermediate result.
Feasible incomplete optimization does not establish minimum nucleotide edits.
Reverse translation and constrained redesign deliberately have different
contracts. `ReverseTranslateProteinSequence` generates candidate DNA from a
protein using sequential codon preferences and advisory local Tm/GC steering;
its preferred/alternative-codon diagnostics establish neither motif avoidance
nor whole-CDS/window-GC compliance. `PlanDnaSequenceDesign` takes an explicitly
authorized synthetic DNA insert with declared protein/CDS/protected bases and
jointly evaluates its requested constraints. Its shared output validator is the
authority for those constrained-output claims, not for search provenance or
expression suitability. A reverse-translated candidate may therefore need a
separate admitted design request; no automatic redesign, approval transfer or
replacement of the reverse-translation algorithm is implied. Keep both report
families and their metrics distinct; preferred codons are not minimum edits.
Solver selection is explicit and version-bound; default enumeration and its
serialized approvals stay compatible. Conflict search branches over all legal
nonoriginal assignments in a current violated window, fixes one codon per level,
and re-evaluates the full sequence. A complete tree, including sound edit-cost
pruning, proves feasibility/optimality or infeasibility; a budget stop does not.
The original input must fit the match-diagnostic cap. An edited candidate that
exceeds it is skipped, including its conflict-search subtree, and makes the
search incomplete; retain a freshly validated best if available, otherwise
report unresolved exhaustion. Never discard a valid best as invalid input or
claim a proof after that skip. Other validation errors and apply checks remain
strict, and cancellation still clears any best result.
No local-score substitution or performance claim follows from fewer visited
candidates. Application checks the selected strategy/algorithm/policy binding
without rerunning search.
Optional GC constraints cover the complete declared CDS, including frozen
start/stop, never flanks or inferred organism defaults. Inclusive basis-point
bounds become integer counts once; all search and apply predicates check both
motifs and GC. Missing GC/output facts are absent, not zero. Disjoint synonymous
count extrema prove impossibility only; overlap does not establish feasibility.
Windowed GC is separately explicit: every full one-base-step window contained
in the CDS, including frozen bases, with no flanks, shortened windows or stride.
Each fixed-length window uses integer count bounds; whole-CDS GC, motifs and
windowed GC must all pass. Conflict branching considers a codon's contribution
inside the conflicting window, not its whole-codon GC count, and retains all
synonyms of direction-capable unassigned codons. Reports retain every exact
window interval/count; apply rederives them. Unrequested/unadmitted facts and
missing output are absent, never zeros. Prefix counting checks complete windows
without sampling or substituting local checks for full validation.
GC presence is explicitly algorithm/policy-bound without changing legacy
GC/window-absent serialization. Recomputed approval hashes cannot bypass output
constraint or GC-fact verification; they still do not authenticate search claims.
Approval binds exact input/output/protection/mapping/configuration/evaluations;
apply revalidates without rerunning search and creates one undoable derivation.
The digest is a content binding, not an authenticated engine signature or proof
of search history. Apply verifies output constraints and exact edits only;
completeness, minimum-edit claims, search counts and reasons remain unverified.
Retain the exact approved report content under
`submitted_proposal={verification:"unverified_portable_preview",report:...}`,
separate from the engine receipt (`output_constraints_verified=true`,
`search_claims_verified=false`). The mandatory typed tag has no verified variant;
it does not change preview or approval bytes. Missing legacy verification fields
establish no claim; historical `proposal` or unlabelled `submitted_proposal`
metadata must not be interpreted as engine-verified search evidence. One apply
owns one new sequence, one named singleton container and one lineage node, all
with the outer operation ID. Inline input has `ImportedSynthetic` provenance;
loaded-source input has `Derived` provenance and a parent edge. All fallible
materialization work precedes state mutation, with no nested public creation
operation or inner detached commit. Host detached execution still checks the
live project instance, structure and journal before committing the outer apply.
Stale source annotations are omitted rather than inherited. No expression,
splicing, folding or suitability claim follows from translation preservation.
See [the scoped contract](dna_sequence_optimization_plan.md). This is outside
the `.12` release gate and is not full DNA Chisel compatibility.

## Transcript Capture Discovery Boundaries

Capture binding coverage, retained isoform information, specificity and
experimental success are separate claims. Target sets, role/window policy,
primer budget and permission to share are independent inputs. Never silently
move an unknown UTR window into coding sequence, add a primer beyond budget,
infer oligo availability from a design, or borrow paired-assay/unique-genomic-
hit approval for a multi-target single-primer pool. Fixed protocol oligos
participate only in declared stages; recorded heuristic Tm is not reaction Ta.
See [the portable capture contract](transcript_capture_pools.md).

## Contextual Launcher Readiness

The egui menu/palette pilot reuses the app's collection-launcher readiness
presentation, including typed collection rejection reasons. It is not a second
protocol capability registry or permission system. Compute bounded selection,
molecule-kind and stored-input facts once per surface; do not call full project
introspection, filesystem probes or scientific preflight from menu paint.
Host catalog validation has its own background snapshot and explicit refresh,
independent of project revision. Invocation rechecks its explicit subject;
shared macro/engine validation remains authoritative for execution. Static
capability descriptor lookup may cache immutable metadata, never host availability.
Importing a template and using it with bound inputs are separate actions.

Object context menus are a primary discovery surface: identify the clicked item
and lead with inspection/opening and applicable setup actions before bookkeeping.
Capture that item for the popup's lifetime, without retargeting on pointer
movement or borrowing an unrelated global selection. Explain disabled
choices and their recovery steps; opening a form is not permission to execute.
Reuse existing shared operations and viewer routes rather than context-menu-only
biology or analysis logic.

## Accession-Pinned TSS Document Invariants

Status: active

Strict TSS profile documents resolve exact, versioned full PFMs once, retaining
the actual input/matrix/source digests. Alias expansion and consensus fallback
remain available only to older permissive paths. Transcript-oriented FASTA
windows use one checked transform for local window starts, TSS offsets and
genomic coordinates; a negative genomic strand never triggers a second sequence
reversal. Missing windows are not measured zero scores. Raw computation and
display clipping are distinct, and report-only rendering cannot rescore data.
Within-factor correlation compares model outputs, not binding or co-regulation.
Bundle consistency is not independent reference authenticity. Existing target
export manifests bind their checksum inventory, not the reverse. Whole-string
genome/dataset expectations and the input source revision remain separate from
the scoring producer. FASTA/manifest transcript sets must be equal; integrated
selection memberships need only be a subset, joined by promoter ID and checked
TSS identity, not by different-length sequence hashes or intervals. Repeated DNA
at distinct loci is retained with a warning. Saved score JSON must round-trip
floating-point bits exactly so export-only replay does not change results. See the
[portable contract and explicit limits](tss_tfbs_profiles.md).

Explicit local rescoring of an annotated TSS window is a separate lane class,
not a refresh of attached report or imported evidence. Use the validated bases
as an `InlineSequence`, retain the actual matrix-count hashes and score kind,
and own settings/jobs/cache independently of other TFBS panels. Valid zero and
unevaluable windows remain distinct; legacy absence of evaluability/provenance
cannot authorize local lanes. Cancellation and source/request changes must not
publish partial or stale results. Navigation never silently triggers rescoring.

An explicit native-view export with local scoring and a subspan scans only
complete motif footprints contained in that span, on both local strands.
Boundary-crossing footprints are excluded, not padded or reported as zeros;
counts retain the same evaluated/possible strand-window meaning. Admission
still uses the full annotated window. Retain full-window identity plus the
exact scored span and its byte digest; identical substrings at different local
positions are not interchangeable. Empirical quantiles use the scanned span.
Cropping existing report/GUI curves remains presentation-only and retains their
original scope, calibration and scales.

Native TSS SVG hover annotations may be bounded only with a disclosed
per-lane/per-local-strand emitted/total/omitted count; the default cap is 128,
retaining strongest and sampled evaluated starts. Every curve point and validity
band remains. Highest raw-score starts (up to three) are a within-matrix,
displayed-span position summary, not a factor or binding-likelihood ranking.
Collection export reuses validated TSS membership and preserves each window's
supplied lane scales, including report-declared scaling; it introduces no new
cross-window calibration. Staged artifact receipts do not authorize changing
approved membership or treating missing flanks as biological absence.

Presentation-size reductions must not discard scores, change sampling or weaken
receipt checks. Raster-backed PDF image compression is lossless; encoding and
output hashes are recorded. The separate `vector_pdf` representation preserves
vector geometry and embeds selectable text, while remaining a static PDF; SVG
title-hover interaction belongs to the SVG/HTML peer and is not promised across
PDF viewers. Receipts bind the exact PDF backend, lockfile, geometry and used
font-source identities. CLI summaries are adapter-only projections of
successful exports, with an explicit full-JSON opt-in. Compute-only results
remain complete when there is no persisted report. GUI/MCP and the shared
engine contract are unchanged, and verification non-claims remain visible.

Resolved TF query tokens remain matrix accessions through shared scoring, even
when several PFMs share one factor name. Integrated reports compare the bound
locus JSON's resolved accessions with the TSS detail panel before composition.
Multi-page raster PDFs audit actual used fonts per page, just like single-page
exports. Compositor outputs are staged before publication, with a receipt-last
completion marker and rollback on caught failures; this does not claim a single
atomic transaction across multiple output paths or protection from power loss.

## 2026-08-01: derive web and print gene-set reports from one record

Publication-style multi-gene analyses use one portable manifest and one
resolved semantic record. Responsive HTML remains the main browsable report;
the PDF is a printable companion containing the full primer list and all
selected figures. This retains the useful PARK7-style overview while removing
hard-coded gene lists and independently maintained HTML/PDF assembly scripts.

This file records durable implementation constraints and architectural
decisions. It is intentionally short: if an entry becomes a backlog item, move
that work to [`roadmap.md`](roadmap.md); if it describes completed history,
move the outcome to [`CHANGELOG.md`](CHANGELOG.md).

Status values:

- `active`: currently governs implementation.
- `superseded`: kept for history; no longer governs new work.
- `deferred`: accepted direction, but not currently active work.

## DEC-001: Single Shared Engine

Status: active

GENtle uses one deterministic engine contract across GUI, CLI, shared shell,
JS, Lua, Python wrapper, MCP, and agent-facing routes. Adapters must not fork
biology or business logic.

## DEC-002: Protocol-First Machine Contracts

Status: active

Machine-facing routes use structured, versioned request/result records. Durable
schemas and wire contracts live in [`protocol.md`](protocol.md), not in the
roadmap.

## DEC-003: Thin Adapter Boundary

Status: active

Frontend and integration adapters should parse, route, present, and collect
artifacts. Core decisions, reports, findings, and explanation records should be
engine-owned and portable whenever practical.

## DEC-004: Roadmap/Changelog/Decision Split

Status: active

`roadmap.md` is for next work only. Completed work goes to
[`CHANGELOG.md`](CHANGELOG.md). Durable constraints go to this file. A session
should be able to start by reading only `roadmap.md` in under two minutes.

## DEC-005: No Roadmap `Done` Entries

Status: active

Do not add a `Done` entry to `roadmap.md`. Move completed work directly to
`CHANGELOG.md` in the same session. This prevents the roadmap from growing back
into a mixed history/planning document.

## DEC-006: Workspace Extraction Order

Status: active

The intended crate extraction order remains:
`gentle-protocol -> gentle-engine -> gentle-render -> gentle-shell -> gentle-gui`.
First-wave extraction avoids per-feature micro-crates; related analysis logic
should remain together until the engine boundary is stable.

## DEC-007: Root-Crate Compatibility During Extraction

Status: active

While production code still partly lives in the root crate, extracted crates
must expose stable contracts without forcing all callers through a one-shot
import rewrite. Root modules may re-export extracted contracts during the
transition.

Extraction does not require runtime dynamic libraries. Ordinary workspace
crates remain statically linked unless a separately reviewed plugin/FFI contract
justifies loader, ABI and packaging complexity. Measure compiler and linker RSS
separately before proposing dynamic linking as a memory repair.

## DEC-008: Shell Dispatch Split Pattern

Status: active

Shared shell parsing and execution should remain one behavior path for GUI
Shell, `gentle_cli shell`, MCP, JS/Lua wrapper helpers, and agent execution.
Large command families may split into helper dispatch functions, but those
helpers must preserve shared parser/executor semantics.

## DEC-009: Stack-Hardening Policy

Status: active

Stack-overflow fixes should reduce dispatcher frame depth without forking
behavior. Prefer narrow helper dispatch and expanded-stack workers for
confirmed small-stack failures in shell/engine routes.

The MCP stdio entry point runs its complete read/decode/dispatch/write loop on
one explicitly sized 16-MiB worker stack and synchronously joins it. This covers
typed operation decoding before shared engine execution, not just the inner
dispatcher. It neither detaches execution nor creates a worker per request.
Transport errors still reach the binary's stderr/nonzero-exit boundary, and
panics are propagated rather than converted to successful responses. Do not
rely on `RUST_MIN_STACK` to enlarge the Windows process main thread.

The standalone CLI likewise runs its complete argument parsing, direct-command
and shared-shell dispatch on one explicitly sized 16-MiB worker, synchronously
joined by main. SIGUSR1 diagnostic setup and the structured stderr/nonzero-exit
boundary remain in main; worker panics propagate. Do not rely on shell/engine
inner workers alone: direct CLI operations can exhaust the main stack before
reaching such a boundary. Real-binary tests must remove inherited
`RUST_MIN_STACK` for both successful operation output and error handling.

GUI dbSNP fetch uses an explicitly sized 16-MiB worker, matching GenBank
sequence ingress, because its detached snapshot enters the same typed engine
dispatcher. A large caller/test stack does not enlarge a nested worker; do not
depend on an ambient `RUST_MIN_STACK`. Worker creation errors clear the pending
fetch task and remain visible rather than reporting a started fetch.

## DEC-010: `#[inline(never)]` Helper Dispatch

Status: active

Until the monolithic root dispatchers are fully decomposed, stack-sensitive
command or operation families may dispatch through dedicated `#[inline(never)]`
helpers before entering broader match frames. Inner match branches should
delegate back to the same helper so there remains one implementation path.

## DEC-011: UI-Intent Shared Catalog

Status: active

UI-intent targets and discoverability metadata are shared catalog data.
Menus, command palette, shell/agent routes, MCP discovery, and ClawBio/OpenClaw
handoff should consume the same target catalog instead of hard-coding parallel
target lists.

## DEC-012: GUI Intent Handlers Stay Thin

Status: active

GUI `ui open|focus|close|selection ...` handlers may open/focus/close host
windows, dialogs, and active viewer selections, but must route target
resolution through shared `ui ...` shell command contracts where possible.

## DEC-013: Mutating Agent/MCP Safety

Status: active

Mutating `op` and `workflow` execution through agent/MCP paths requires
explicit confirmation. Agent-suggested commands must not recursively invoke
`agents ask`, `agents plan`, or `agents execute-plan`, including through macro
expansion paths.

The stored-plan executor checks the shared parsed command's mutation policy,
not just model-supplied `mutating` / `requires_confirmation` flags. Shell,
UI-intent, raw-operation and workflow candidates cannot relabel project changes
or classified external artifact writes as read-only to bypass `--confirm`.
Explicit candidate confirmation requirements still apply. Named read-only
inspections keep the shared policy's exemptions; direct MCP `op` retains its
independent confirmation requirement.

## DEC-014: Screenshot Capture Policy

Status: active

Project policy permits only user-approved, window-bound screenshots for Agent
Assistant help. A direct request may originate from an `Agent help` control in
the GENtle viewport being diagnosed. An agent response may alternatively carry
one typed, path-free request id and human-readable reason. That response emits
no capture command: GENtle shows a consent card, the user selects one currently
registered GENtle content viewport, and `Allow one screenshot` consumes a
turn/provider/project/viewport-bound approval before requesting one egui
capture. The agent cannot choose an arbitrary target or grant approval.

Both structured project-state context and screenshots disclose information to
the selected provider; this policy does not assume that a screenshot is more
sensitive. A project summary can be scientifically richer, while a visual
snapshot is often more useful for layout, direction, colour, and enabled-state
guidance. They remain separate review controls because their scope and
diagnostic purpose differ. Captured images are previewed locally and sent only
when the user explicitly submits the follow-up prompt. Conversation history
stores path-free request/attachment provenance, not live approval, image bytes,
or the temporary local image path.

Native macOS full-window capture is an optional same-process enhancement. It
must use ScreenCaptureKit, refuse windows not owned by GENtle, and remain
unavailable when the user does not grant macOS Screen Recording permission.
Automatic/background capture, desktop capture, capture of another application,
agent-selected capture, and headless capture remain prohibited. The reserved
`screenshot-window` CLI/shared-shell route therefore stays disabled;
deterministic SVG/export routes remain preferred for reproducible artifacts.

## DEC-015: ClawBio/OpenClaw Boundary

Status: active

The ClawBio/OpenClaw skill scaffold wraps deterministic `gentle_cli` routes and
writes reproducibility bundles. It should surface shared capabilities,
suggested actions, and UI-intent handoffs without creating ClawBio-only biology
logic.

Every wrapper run carries a topic-neutral, content-bound provider execution
receipt. The receipt binds the resolved request, immutable inputs, mutable
state before and after execution, wrapper/runtime identity, command outcomes,
native result identifiers/status fields, and emitted artifacts. It is not a
discussion ledger, permission record, participant claim, freshness model, or
scientific interpretation. External systems may retain the receipt, but own
their own request identity, policy, and interpretation records.

Capability-maturity discovery is a generated projection of explicit intent and
plan-step bindings plus a separate public evidence archive. Stable IDs refer to
existing registry identities; missing or ambiguous identities remain unassessed.
Parity/reachability, tutorial review and release prose never become execution,
native, performance or biological passes automatically. Discovery does not
probe a runtime and always remains unassessed. Scoped assessment requires the
exact capability definition, tested source/binary, platform/architecture,
profile/features, input class and resource hashes for one evidence kind.
Historical mismatched evidence remains inspectable, not applicable; conflicting
fail/incomplete/withdrawn evidence cannot be overridden by a passing record.
Content hashes bind evidence but do not authenticate its author or correctness.
Availability, maturity, execution outcome and existing approvals remain separate.
No maturity field grants mutation, network, order or release permission, and no
feedback collection/sharing or telemetry is introduced. Optional metadata must
not disable ordinary discovery/execution when missing or malformed. Existing
descriptor/request digest algorithms stay unchanged; an updated descriptor must
be co-deployed and approved afresh where its bytes participate in approval.

When a caller needs a policy-aware analysis-run record, the provider receipt
and caller-owned run remain two records joined by content hashes rather than
one shared mutable object. A post-run skill adapter may prepare a typed intake
packet, but it must not discover or mutate the caller's store. The caller must
revalidate current references and permissions, recompute the request digest,
and assign replay, dependency, and freshness state through its owning commands.

Descriptor-only skills may delegate to the generic runtime only through an
explicit versioned contract whose co-deployed catalog, descriptor, selected
route, and plan step verify before execution. This records the route that was
selected; it does not claim deterministic replay of natural-language routing.
Wrapper execution outcome and native biological status remain separate, so a
successfully computed scientific `fail` is neither an execution failure nor a
pass.

When a ClawBio presentation combines scientific artifacts, the raw GENtle
result remains authoritative. Strict-attribution presentations must label
executable-derived statements separately from wrapper presentation, caller
inputs, named external tools, and unverified prose. Every executable-derived
summary requires an exact native-result pointer or scope, content hash, and
named deterministic projection in a machine-readable claim/artifact ledger.
Wrapper-authored graphics may arrange GENtle artifacts but must not recalculate
or upgrade their scientific content.

Adapter execution timeouts are approval-bound operational ceilings, not
scientific parameters or expected runtimes. Long primer-design routes must
allow the ceiling to be reviewed explicitly; adapters must not present a fixed
wrapper timeout as a Primer3 limitation. Runtime expectations should consider
transcript-equivalence classes, target ROIs, requested junctions, assay mode,
and operation count in addition to sequence span.

## DEC-016: External AI/Automation Deployment

Status: active

For tool-driven external AI deployment, prefer the published GHCR image in
headless `mcp` mode with explicit project/state mounts and stdio communication
instead of depending on browser GUI containers.

The maintained OCI distribution builds only `gentle_cli`, `gentle_mcp` and
`gentle_examples_docs`, with Cargo default features disabled and no embedded
JS/Lua features. It retains the Python CLI wrapper and scientific helper tools.
Native desktop packages retain the GUI; embedded scripting is an optional
source build, not part of the current package inventory. The container dependency
guard must reject desktop and embedded scripting runtimes before compilation.
New bare release tags and `latest` refer to this headless image; historical
`gui` tags are not refreshed or presented as current releases.

## DEC-017: Fixture Provenance

Status: active

Every committed fixture under `test_files/` or `tests/` must document origin,
deterministic recreation/retrieval steps, and where GENtle uses it.

## DEC-018: Source Documentation Policy

Status: active

Public modules and key public records should carry concise rustdoc explaining
purpose, invariants, and non-obvious behavior. Prefer why/invariant/edge-case
documentation over line-by-line narration.

## DEC-019: Primer/Oligo Material Identity

Status: active

A planned primer/probe or in-silico design artifact is not physical stock.
GUI, CLI, JS, Lua, Python, MCP, and agent routes must preserve the distinction
between design, reviewed order, received material, local availability, and
consumption.

## DEC-020: Inline/Stateless Sequence Inspection

Status: active

Read-only operations that need only sequence letters, optional topology, and an
optional span should be state-optional. Promoting such results into project
state remains an explicit second step.

## DEC-021: Helper-Construct Terminology Migration

Status: deferred

The legacy "helper genome" vocabulary should eventually move toward
"helper construct" semantics. The migration should be atomic across contracts,
docs, GUI, and examples to avoid mixed terminology.

## DEC-022: Primer3 Backend Parity

Status: deferred

Primer3 is the intended external backend parity layer for primer design, but
deeper constraint mapping, fixture coverage, and backend equivalence checks
should land through shared engine/report contracts rather than GUI-only code.

## DEC-023: Presentation i18n Does Not Localise Machine Contracts

Status: active

GUI labels and dialog text may be translated at runtime, but shared shell
commands, protocol schema fields, saved project records, adapter payloads,
agent response schemas, and scientific identifiers remain deterministic
English. Localised text is presentation only.

## DEC-024: Introspection Composes Facts And Capabilities

Status: active

Words-only and agent-facing introspection uses additive fact domains rather
than renaming existing fact vocabulary. Existing project facts keep their
`sequence.*`, `report.*`, and related names with `domain: "project"`;
view/host/config facts use the same registry shape. Headless contexts must
project `ui.host_available=false` explicitly so GUI/view-intent readiness fails
with a named unmet fact rather than ambiguous absence. Capability readiness is
computed by instantiating templated fact atoms from declared arguments and then
using the shared fact evaluator; unbound catalog readiness reports missing
argument bindings as `unknown`. Host-domain availability facts must project
from existing deterministic preflight/probe logic rather than creating a
parallel detector.

Partial bindings preserve expression structure and the shared three-valued
truth rules: false OR unknown is unknown, true OR unknown is true, and NOT
unknown is unknown. Binding diagnostics remain inspectable even if another
alternative decides the result; a flat list of missing arguments is not a
readiness veto. Mutually exclusive TSS SVG targets are an explicit exception
at binding selection, not Boolean evaluation: supplying exactly one `SEQ_ID`
or `COLLECTION_ID` selects its prerequisite and omits the unselected target.
Neither target supplied remains unknown; operation admission still rejects
both targets supplied. No placeholder project fact may stand in for an
unbound argument.

Content-addressed approval workflows whose truth depends on external files,
catalog bytes, or an exact proposal digest use an operation-specific,
context-bound readiness report rather than a durable Boolean project fact. The
readiness operation must rerun the canonical validator/planner against current
state and preserve `not_supplied` separately from `blocked`; otherwise an old
fact could incorrectly authorize a changed proposal.

Agent-context projection must retain complete fact counts and explicit omission
metadata, without altering the canonical fact graph or treating omitted facts
as false. Request relevance may reorder rows within fair per-type allocations,
but must expose its deterministic reasons and never imply biological importance,
new facts, or readiness. Execution feedback is host-observed evidence, not model-authored
conversation or authorization. Receipt history and current applicability are
separate: a completed command can need rechecking after edits; equal project
revisions do not validate changed external files. Session receipts must not be
persisted through metadata writes that would clear redo or change dirty state.
Model-facing receipts contain bounded typed outcomes and hashes, not raw
commands, result bodies, or error strings. Successful dispatch is not background
job completion, GUI presentation, or scientific acceptance.
Host projection omits receipts without both a visible turn and a positive
suggestion index, retaining them locally rather than exposing anonymous outcomes.
Helper-catalog retrieval ranks current-prompt matches before the latest user
message and its assistant reply; older dialogue cannot displace current intent.
This bounded lexical ranking does not establish biological relevance or readiness.

## DEC-025: Legacy SHA-1 Download Verification Uses External Tools

Status: active

GENtle must not use an in-process SHA-1 implementation to verify external
downloads. If an upstream source or compatibility manifest still offers only a
SHA-1-style integrity field, GENtle may delegate that legacy digest to platform
tools such as `sha1sum`, `shasum -a 1`, or Windows `certutil`. If no suitable
tool is available, the file is not SHA-1 verified; preparation may still rely on
basic corruption checks such as successful gzip decompression, size checks, and
parser validation. `GENTLE_SHA1_TOOL` may point to an explicit tool, and
`GENTLE_DISABLE_LEGACY_SHA1=1` intentionally disables the legacy check. Legacy
SHA-1-shaped local provenance fields are compatibility fields; regenerated
local fingerprints should use algorithm-prefixed stronger digests such as
`sha256:<hex>` without requiring an in-process SHA-1 dependency.

## DEC-026: Long GUI Jobs Use Optimistic Engine Snapshots

Status: active

Long-running GUI network and computation jobs must not hold the shared
`GentleEngine` read or write lock for their full run. A mutating worker forks
the engine state, journal, operation counter, revisions, and history limit,
executes against that detached snapshot, and commits with a short write-lock
swap only if structural project data and journal history are unchanged.
Inherited undo/redo checkpoints stay in the live engine; commit appends only
checkpoints created by detached operations. Concurrent viewport/display edits
are preserved across the commit; disjoint UI/session metadata is merged into
the result and its undo checkpoints, while conflicting metadata makes the
result stale. Other stale results never overwrite intervening edits.
Persisted-state mutation revision is separate from structural and read-only
execution identity so GUI dirty checks remain constant-time without marking
inspections as edits.

Observation-only workers use a history-free engine clone and do not commit it.
Managed generic Op/Workflow commands are normally journal-publishing workers,
including covered Plan/S1 inspections except for Slice 5 below: they retain a
project snapshot and the
same guarded commit, owner/history/structural checks and cancellation boundary.
An inspection commit preserves persisted content, mutation/structural revisions
and undo/redo, but advances the journal/execution identity. Concurrent journal
publication can still make another job stale; Slice 3 does not waive that guard
or move these commands into the independent observation lane.

Owner-authorized Slice 5 exception (2026-10-28): only a standalone managed
`PlanDnaSequenceDesign` with `path=None` publishes a bounded command result
instead of its private journal. It still requires the project snapshot and an
Execution slot; publication checks owner/structural freshness and cancellation
under engine-read then job locks. It does not compare journal length there or
change the generic commit/admission guards. The receipt records `preview`, the
snapshot owner/revision and zero committed operations. Such previews add no
project History row; command receipts are bounded, local and nonpersistent.
Same-project structural edits after publication make the retained preview historical;
display and explicit model-sharing recheck its host-local binding. Replaced-
project results are not reattributed/shared. Scientific report hashes/statuses,
raw adapter journals, exact Apply validation and approval remain unchanged.
With-path previews, S1, workflows, opaque wrappers and arbitrary trusted work
remain journal-publishing controls. This is not RNA's project-independent lane,
a generic ReadOnly exemption, journal merge or authorization of editor biology.
See the [scope and pending checks](operation_lifecycle_slice5_plan.md).

A raw
`GentleEngine::clone` in a background worker is prohibited because it retains
unneeded history (even with shared checkpoints). Persisted metadata mutations outside
`apply` must advance the mutation revision and clear redo; otherwise a redo
checkpoint captured before the side-edit can overwrite that metadata.

Transactional macro rollback restores the complete pre-run engine baseline,
including the journal and both history stacks, on error or unwind. It advances
execution/mutation/structural revisions beyond the baseline and final failed
state, retaining the larger baseline/final operation counter. Private history
checkpoints may be shared immutably; rebasing them requires copy-on-write.
The public macro failure receipt is recorded after rollback. This guarantee
does not undo external files, resource preparation or running jobs.

[DEC-048](#dec-048-command-submission-and-observation-never-wait-for-execution)
extends the nonblocking requirement to command submission and observation across
all adapters; this snapshot/commit rule still governs project-changing workers.

## DEC-027: RNA Allele Evidence Does Not Invent Phase

Status: active

RNA-read allele screens may use reviewed local variants to extend reference
hashes only when the variant can be projected to an explicit transcript
coordinate and its reference allele matches that transcript. Unphased
genotypes remain allele-level evidence, and disconnected phase sets remain
independent blocks; neither may be combined into a synthetic gene-wide
haplotype. Reports must distinguish physical reads, fragments, and evaluated
evidence observations. Structural rescue matches are non-blocking candidate
evidence and must not be labelled as novel exons, retained introns, fusions, or
resolved isoforms without a separate alignment and junction-validation layer.
Annotated-junction and exon-intron-boundary catalogs index only k-mers that
cross an explicitly declared boundary. Retained-intron and intronic candidates
remain reportable, but same-read or same-fragment support from an annotated
splice junction for that gene is recorded separately as an RNA anchor. Such an
anchor strengthens an RNA-derived interpretation without excluding pre-mRNA or
genomic-DNA contamination and is not itself a retained-intron call.

## DEC-028: Construct Reasoning Outputs Are Fingerprinted Snapshots

Status: active

Construct-reasoning graphs persist their emitted inspection actions as part of
the graph snapshot. Loading or normalizing a non-empty action list must not
silently rerun current recommendation rules. New graphs fingerprint the full
source sequence/feature snapshot and normalized objective with SHA-256 and
record the reasoning rule-set version; readers compare those identities to
report `current`, `stale`, or `unknown`. A fingerprinted stale graph remains
inspectable but its executable inspection actions must be refused until its
reasoning build route is rerun to produce a new snapshot.

Task-specific reasoning keeps intrinsic evidence severity separate from
objective applicability and effective priority. Typed
`ConstructObjective.intended_tasks` is authoritative when present; legacy
free text may conservatively infer positive applicability but must not infer
non-applicability. Curated repeat provenance belongs in typed repeat annotation
and repeat-family records, with legacy note parsing retained only for backward
compatibility.

## DEC-029: Primer Specificity Is Target-Space And Database-Content Bound

Status: active

Primer-specificity evidence must name its biological target space. Genomic-DNA
screening assesses genomic carryover and unintended genomic products;
transcriptome/cDNA screening assesses cross-amplification among spliced
sequences. A junction-spanning primer can therefore pass a genomic screen
without a contiguous intended genomic product. GENtle must identify an intended
genomic product by prepared-FASTA subject plus genomic binding geometry, never
by equality to the expected cDNA amplicon length.

Specificity searches must not use BLAST `-max_target_seqs` as a substitute for
a post-search hit-review threshold, because doing so can silently exclude
contigs. Handoffs and finalization bind the prepared BLAST database's index
kind and deterministic content fingerprint in addition to its prefix and genome
label. Replacing database content at the same prefix invalidates a pending
handoff. Whole-prepare activity state does not override successful
component-level `blastdbcmd` validation.

Intended-target geometry carries its original reference/assembly/release. A
matching chromosome label cannot transfer coordinates to another assembly or
release. Explicit mappings bind the searched database, are labelled
caller-provided, and are not an automatic liftover or orthology claim. Missing
or incompatible binding stays `unknown`/`not_assessed`.

Existing panel readiness selects one exact reference independently for genomic
carryover and transcriptome specificity through companion selection receipts,
not report timestamps or communication-summary labels. Receipts bind the panel
design, full pair and annealing/tail interpretation, policy, sequence and
annotation fingerprints, and validated acceptance/report bytes. Redundant
communication summaries are excluded from selection design identity; existing
approval/handoff digest rules are unchanged. Explicit complete finalization can
replace a selection; superseded evidence is retained and incomplete finalization
cannot activate it. Legacy evidence without a receipt is inspectable, not
passing. Explicit new handoff construction revalidates current references;
standalone multi-reference summaries do not attach to these panel selections.
Their genomic/transcriptome verdicts remain independent: required failures persist,
unavailable/incomplete evidence is never a zero/pass, optional-only dimensions are
not required, and partial external execution cannot create an aggregate pass.
The parent binds expected commands/paths and retained-byte output receipts;
scientific summary identity excludes machine paths and execution timestamps.
Historical display does not probe resources or assert present applicability.
Import performs fresh source/reference checks, but a hash is not authentication
of an external process or a claim of laboratory/isoform/order validity.

Import parses and hashes the same retained output bytes. File hashes bind
content, not the authenticity of an external runner's execution claim.

Read-only panel-specificity inspection projects retained receipt snapshots,
not whichever report currently shares their ID in the global store. Each
selection keeps its reference, policy and design scope. Binding diagnostics,
historical verdicts and active-selection metadata remain distinct; no combined
panel pass or requirement policy is invented. Historical design evidence and
unlinked same-oligo standalone summaries cannot fill current-design evidence
gaps. This projection never records a selection, probes references or changes
the existing readiness contract.

Opt-in primer-import tracing is observational, operation-local and bounded.
It uses fixed phases/numeric counts only, never private identifiers, report
fields, content fingerprints or persistence. Wall-clock spans overlap and do
not establish CPU attribution or biological pass; early exits/truncation stay
explicit. Trace I/O failure must not change validation, error or state results.

## DEC-030: Collection Lifting Is Subject-Specific And Engine-Owned

Status: active

Collection-visible capabilities declare their lifting behavior in the typed
`gentle.collection_lift_policy_registry.v1` registry. Policies are keyed by
capability and subject kind because a logical gene set, physical container, and
ordered arrangement do not have interchangeable semantics. Supported modes are
`map`, `combine`, `compare`, `arrange`, and `derive`; incompatibility is an
explicit typed rejection rather than an adapter-local conditional.

Dynamic readiness failures remain operation/fact-graph results and must not be
misclassified as permanent policy rejection. Generic collection reports carry
typed member outcomes/errors and links to produced domain reports. When a
domain report has no persisted store/getter, a portable domain wrapper may own
the child reports directly; in that case `produced_report_ids` stays empty so
GENtle never emits dangling report references. Membership locks are set-like
except for explicitly ordered subjects such as arrangements, whose numeric
member order is part of the fingerprint.

TSS-derived collections are validated references to project-sequence sets, not
raw ID aliases. Their approvals/member snapshots bind authoritative sequence
records, overhangs and genomic anchor geometry, never recomputable caches or
unordered cache maps. Legacy cache-sensitive snapshots require explicit new
approval, not silent migration. Typed collection maps validate all members
before scanning; metadata-only forgetting preserves sequences and lineage.
Registry discovery is metadata-only and never implies member validation.
Legacy/malformed entries remain discoverable with explicit diagnostics, and
unavailable metadata counts are distinct from zero. A validating inspection is
bound to the current project instance and structural revision in the GUI.
Transcript 5-prime truncation must survive extraction, including dropped first
exons. A surviving exon boundary cannot become an exact annotated TSS merely
because a file library clips locations without marking them partial.
Inventory and extraction share uniform-strand endpoint interpretation, including
both complement/join encodings; neither may use a majority-strand guess for an
exact TSS. Missing gene IDs may be associated only through an unambiguous
explicit label/ID link in the same loaded locus, source and strand. Unlinked
annotations remain locus-level unassigned records, not gene-specific exclusions.

Collection wrappers over non-persisted child scans own those reports directly
and must describe aggregate completeness. A retained-row cap is not evidence
that later motifs or candidates were scanned: aggregate count fields must be
named as retained/observed values and accompanied by an explicit completeness
flag plus affected member ids or warnings.

Versioned collection result wrappers are additive output contracts: they
default missing fields and tolerate unknown future fields. Collection request
and member-binding records remain strict and may reject unknown fields so a
misspelled mutation or binding parameter cannot be silently ignored.

A mutating collection lift must expose a non-mutating preview before apply.
Its apply lock binds the exact source snapshots, effective parameters, and the
deterministically reserved output namespace in addition to the collection
membership fingerprint; the membership fingerprint alone is not a mutation
guard. Apply must reject a changed plan before inserting any output. Derived
sequence lineage is recorded from each source member to its own products, and
logical collection products do not imply that those products occupy one
physical container.

Shared digest convergence guards must be deterministic. Direct and lifted
digests use fragment-count, round-count, repeated-state, and no-progress
limits; wall-clock load must not decide whether an otherwise identical digest
succeeds.

A collection `combine` that writes one external artifact is atomic rather than
partially successful: all members validate before the single write, and a
successful report links every member to the same wrapper report. Pool export
requires an explicit physical container whose declared contents are exclusive.
In-silico `selection` containers are not physical pool sources and are rejected
before export; physical `singleton` and `pool` containers remain eligible.
`gentle.pool.v1` has no incompleteness/provenance slot, so exporting a
non-exclusive known subset is rejected before writing rather than silently
losing that meaning at the artifact boundary.

Gel-render lift modes follow the resulting lane geometry. An explicit
container is a `combine` into one physical lane; a serial arrangement is an
`arrange` that preserves its declared lane and role-label order. A bare project
sequence selection is not a valid collection lift because the direct input
form co-migrates every supplied sequence in one lane without proving that the
sequences share a physical sample. `ArrangementMode::Plate` incompatibility is
a dynamic readiness/execution error, not a permanent arrangement-subject
policy rejection. The shared `render-pool-gel-svg --arrangement` route remains
executable even while the generic collection launcher reports that it has no
typed arrangement adapter.

Explicit gene-set pooling is a physical aliquoting operation, not an inference
from logical membership. Every resolved member must bind to a distinct,
exclusive singleton source container, and incomplete resolutions reject the
whole plan. Apply creates one lineage-linked derived aliquot sequence per
source, places only those aliquots in one new exclusive `pool` container, and
retains the source containers and their latest-membership mappings. The plan
lock covers the resolution membership, source-container bindings, exact source
sequence snapshots, output ids, and reserved pool-container id.

Collection members bind biological interpretation through a report-owned
context registry. Context-sensitive policies are fail-closed: an undeclared
policy is `not_reviewed`, while reviewed consumers may require one homogeneous
context and reject missing, conflicting, or target-mismatched contexts before
coordinate lookup, sequence search, or other biological work. Generic portable
collection reports copy the registry because their source may otherwise be
available only by report id; domain reports embedding the source resolution do
not create a second registry. The canonical collection membership fingerprint
continues to lock membership/order only and must not be interpreted as a
context or biological-input fingerprint.

## DEC-031: Ortholog Context Binding Extends The Specialized Resource

Status: active

The existing offline `gentle.ortholog_resource.v1` remains the authoritative
ortholog mapping contract. Optional source/target context ids bind each
directional endpoint to the resource's biological-context registry, and
resolved cohort reports copy referenced contexts so they remain portable.
Explicit organism or genome conflicts fail before genome-catalog work.

Orthology type and confidence are open typed string vocabularies: GENtle
recognizes canonical cardinality/confidence values while preserving unknown
provider-specific text exactly. Legacy symbol-only rows remain accepted, but a
matching symbol is never treated as relationship or functional-equivalence
evidence. The shipped `reject` and deterministic warning-bearing `first`
ambiguity policies remain unchanged. The additive `preserve` policy leaves an
ambiguous target unresolved while carrying every ordered candidate mapping and
its portable context/provenance; candidates do not become resolved cohort
members.

Unlike provider-owned orthology type and confidence text, ambiguity policy is
a closed GENtle operation-control vocabulary. The v1 report schema is retained
because `preserve` is opt-in and does not change the meaning or wire shape of
existing `reject` and `first` reports. Readers must reject an unknown policy
rather than silently defaulting it; consequently, a pre-`preserve` v1 reader
cannot consume a report that explicitly requested `preserve`. Candidate ranks
and labels are presentation references within one report, not durable
cross-run identities.

A generic relation hierarchy or paralog contract still requires a concrete
resolver/consumer and is not inferred from the ortholog contract.

## DEC-032: Conversational Scientific Delegation Is Proposal-Then-Execute

Status: active

Natural-language routing may draft a structured GENtle request, but it is not
execution authority for a delegated route that selects biological material, a
pair or backend, mutates project state, or writes/replaces an artifact. Such a
route first emits a content-addressed proposal binding the normalized request,
exact command, delegated descriptor/route, runtime, state, scientific inputs,
resolved paths, reference context, and pinned backend. Execution consumes that
stored request only after a caller-supplied approval assertion names the exact
proposal digest; it must not reconstruct slots from prose a second time.
For this gated path, known local Cargo launchers are reduced to the already
built, hash-bound `gentle_cli` executable before the proposal is issued, and
OCI execution requires an immutable image digest. Approval must never authorize
a later source rebuild or mutable image-tag resolution.

GENtle verifies the digest and fails on bound-context drift. The caller or
OpenClaw control plane owns approver identity, authentication, authorization,
and proof that a human actually reviewed the proposal. Read-only list, show,
and preflight routes may remain automatic. Direct structured GENtle and
`gentle_cli` use remains independent of this conversational delegation gate.
Material cache/tool environment overrides are part of the approved context;
GENtle does not claim to bind arbitrary unlisted operating-system state.

## DEC-033: Isoform-Assay Planning Has Two Scientific Approval Digests

Status: active

Conversational gene isoform-assay planning separates recommendation from
execution. The first approval binds the complete normalized planner request,
including policy schema/version/hash, evidence report ids and hashes, prior
plan, observations, override, and every effective default. It authorizes
planning from that exact basis; it does not validate the biology. The second
approval binds the exact ordered `DesignTranscriptAssayPanel` payloads emitted
by the approved plan. Those payloads are not regenerated or altered after
approval, and one batch approval is valid only because its digest covers the
complete ordered set.

An optional pre-approved informative fallback is part of both approval bases,
not a runtime convenience flag. The normalized request digest binds its typed
policy and whole-panel budget. The emitted operation digest binds the exact
wrapper, unchanged strict child, deterministic fallback template, and distinct
report identity. Execution may derive that child only after a matching typed
strict coverage-infeasibility record; all other failures remain fail-closed.
The strict failure and partial fallback are separate durable records, and a
bounded-greedy selector certificate must never be presented as a proven global
optimum or as experimental/order readiness.

For a multi-study second stage, independently approved plan/workflow pairs are
first composed into one engine-owned batch record. That record binds each
plan/workflow file, their existing workflow and operation digests, their order,
and the combined ordered operation set. Execution revalidates the complete
batch before applying its first operation and runs it through one state-bound
command. An adapter must not create several execution proposals against one
pre-mutation project state and then run them serially; the first mutation would
make the remaining approvals stale.

An oligo order form may be projected from this work only through a specific
`gentle.experimental_assay_handoff.v1` whose bytes and named readiness policy
are retained. Every order line keeps the handoff hash, readiness-policy id and
hash, readiness row/card, and assay/pair/oligo identities. Presentation
profiles and content-block selection may project an already immutable
canonical report without another scientific approval; they still produce an
execution receipt and remain subject to ordinary file-overwrite safety. An
outer agent may host or shorten those declared projections, but cannot author
new scientific blocks.

## DEC-034: Approved Batch Reuse Is Explicit, Exact-Prefix, and Detached

Status: active

Long-running approved isoform-assay batches execute against a detached engine;
by default the live project is committed only after the complete ordered batch
succeeds. DEC-035 adds an explicit opt-in that lowers this unit from the batch
to the gene without removing it.
An explicitly configured private checkpoint directory may retain the detached
state and exact operation results after each successful approved operation.
These artifacts survive a later rollback, but they are not project mutations
and are never discovered or imported during ordinary execution.

Reuse is a separate read-only proposal followed by explicit approval of that
proposal digest. The proposal binds the target batch, baseline project hash,
exact approved-operation prefix, checkpoint manifest and engine bytes, and the
running GENtle executable identity. Execution revalidates every binding and
imports only into another detached run. Adding work after an unchanged prefix
is supported; arbitrary branch-state merging, cross-build reuse, payload drift,
and approval-based waivers of provenance mismatches are not. Approval
authorizes transfer of prior computation and does not validate the biology.

## DEC-035: Per-Gene Continuation Lowers Atomicity, It Does Not Remove It

Status: active

`--on-gene-failure continue` is opt-in, and the default stays `abort` so an
existing command keeps its exact behaviour and output. Under `continue` the
unit of atomicity moves from the batch to the gene: a failing gene is restored
to an in-memory boundary snapshot taken when that gene started, so a gene is
never committed half-applied, and the remaining precomputed genes still run. A
batch in which no gene completed commits nothing.

This rollback guarantee is limited to engine state. Therefore `continue`
refuses the complete batch before its first operation if any approved operation
declares an external output path; a file already written cannot be made atomic
by restoring an in-memory engine snapshot. Callers must omit per-operation
paths or retain the default batch-level `abort` policy. The boundary snapshot
also includes the reported checkpoint pointer, so a failed gene cannot leave
`last_checkpoint_manifest` naming a partial-gene checkpoint that was rolled
back. More strongly, `continue` persists checkpoints only at complete-gene
boundaries; it never writes a partial-gene checkpoint that could remain
discoverable after rollback. The default `abort` policy intentionally retains
per-operation checkpoints because those are its recovery units.

The flag governs execution failures only. Approval and verification failures
still refuse the whole batch under either policy, because they say the request
was not the approved one rather than that a result was not obtained. The
requested coverage policy is likewise untouched: coverage is a biological
statement about one assay, this flag is an execution statement about a batch,
and a gene failing `require-all` is skipped rather than downgraded. A report
from a run with failures sets `batch_complete: false` and never presents itself
as a finished study.

Checkpoints keep describing an exact contiguous ordered prefix. Rather than
teaching the checkpoint schema about gaps — which would cost the invariant that
makes reuse trustworthy — checkpoint writing freezes at the last contiguous
prefix once a gene has been skipped. Because such a run also commits the genes
that succeeded, its own frozen proposal is correctly refused afterwards, and
recovery is a new batch of only the failed genes against the committed
baseline. For the same reason `continue` refuses to resume from a reuse prefix
that ends inside a gene, where no boundary snapshot can exist.

## DEC-037: Primer3 Is Selected By Capability, Then Version, Then Recency

Status: active

Where several Primer3 builds are installed, taking the first `PATH` match makes
GENtle's behaviour depend on shell configuration. A bare configured name is
therefore resolved by ranking every candidate: a build advertising `--progress`
beats one that does not, because bounded search reporting and SIGUSR1 status
are built on it and a legacy build silently degrades them; then the highest
version, compared component by component rather than lexically; then the newest
file, which resolves the remaining ties without reintroducing `PATH` order.

A configured value naming a path is authoritative and is never overridden by
discovery, because an operator pinning a build has made exactly the decision
this ranking would otherwise take back. The selection is reported rather than
implicit: preflight lists every candidate with its capability, version and
creation time, and states why the winner won, so a result can be traced to the
build that produced it. It is also cached per candidate set and file identity,
because checkpoint reuse binds the running executable and must not see a
different choice mid-run.

## DEC-038: Collection Construct Inspection Does Not Persist Temporary Graphs

Status: active

Collection construct-reasoning inspection is a `map` over project sequences or
homogeneous gene-set resolutions. It builds each member's current reasoning
graph on an engine clone and returns graph identity, input fingerprint, and
inspection actions in the collection wrapper, without writing those temporary
graphs into project metadata. Consequently `graph_persisted` is false and the
generic member status does not advertise a `produced_report_id` that cannot be
retrieved later.

Project-sequence inspection is context-agnostic because the operation has no
reference-genome parameter; a gene-set resolution must nevertheless be
homogeneous so one report does not silently mix declared biological contexts.
Containers and arrangements are explicitly rejected: physical occupancy and
member order do not change independent per-sequence reasoning, and accepting
them would imply unsupported physical semantics.

## DEC-036: A Dossier States What Is Missing Rather Than Waiting For It

Status: active

Experimental work starts from whichever genes are ready, so an isoform-assay
dossier is publishable while a study is incomplete. Each gene carries a
`resolved`, `pending`, or `unresolved` status; `unresolved` requires a reason
stating why the established automatism could not address that gene. The index,
each incomplete gene page, and the print document say so visibly, and the
canonical report carries per-status counts and a `complete` flag.

An omitted status is derived from the gene's own contents rather than assumed
to be `resolved`, so a request written before this field existed reports
honestly instead of rendering a gene without results as finished. A declared
status is validated against those contents for the same reason: the label is
only worth having if a reader can rely on it, which is why `resolved` without
any assay handoff is an error rather than a silently empty section.

## DEC-049: Primer Similarity Guides Search But Does Not Replace Specificity

Status: active (formerly duplicated as DEC-038)

Background: Ye et al., *Primer-BLAST: A tool to design target-specific primers
for polymerase chain reaction* (BMC Bioinformatics 2012),
https://doi.org/10.1186/1471-2105-13-134.

A broad cDNA or genomic BLAST search may identify loci and intervals where a
primer is likely to cross-amplify, and GENtle may use that evidence to exclude
or deprioritize Primer3 search windows. The planning artifact must bind the
database content, query templates, coordinates, and exact treatment of each
interval. Similarity alone never proves paralogy, target identity, or assay
specificity.

Final specificity remains a separate complete target-space gate. BLAST finds
candidate loci; an explicit full-primer alignment may adjudicate their complete
binding geometry and terminal differences. A reviewed off-target allowance is
exact and evidence-bound, never a family-wide shortcut and never a way to turn
an incomplete search into a pass. Replacement primers proposed after failure
must re-enter the same cDNA and genomic gates.

For groups of caller-supplied related sequences, exact common primer placement
is derived before Primer3 from recorded global alignments. The longest member
is merely the default representative; it is not declared canonical. Alignment
and candidate budgets are explicit and content-bound, and `require_all` cannot
be weakened by runtime fallback. The resulting full pair-by-member product
matrix remains distinct from later genomic and transcriptome specificity.

Primer reports may additionally project accepted candidates onto a bounded
Pareto frontier. This is an explanation of non-dominated tradeoffs after hard
constraints, not a second acceptance rule. The declared design objective and
final specificity/readiness gates remain authoritative.

## DEC-039: External Auditor Owns The Performance Verdict

Status: active

GENtle developers provide deterministic benchmark targets, fixture provenance,
source-revision binding, and compile/smoke evidence. They do not certify the
performance of their own implementation from a development machine. A named
external auditor runs timed benchmarks on a stable host, retains the raw
Criterion artifacts and environment metadata, compares only compatible
baselines, and owns the release-facing performance verdict.

Criterion is a CPU-side characterization tool, not proof that the GUI is
responsive. GUI acceptance also requires the real application, representative
project state, window creation, repaint and resize behavior, and profiler
evidence where needed. Shared CI runners may compile or smoke benchmark targets
but do not establish strict timing gates without repeated stable-host evidence
and an auditor-approved threshold.

Routine thin-LTO `bench-audit` binaries and native release-profile binaries
are distinct evidence classes. Historical fat-LTO/single-codegen-unit and
Cargo-default (`lto=false`) release measurements retain their original profile
identity; neither the 2026-09-22 simplification nor the subsequent native
`lto="off"` change relabels those results. Baseline repeats are
comparable only when profile, host, toolchain, fixture hash, and GENtle revision
match; retained metadata and baseline labels must identify the effective profile
explicitly.

Native latency acceptance uses the actual packaged GUI and its effective
build recipe, not a Criterion harness, even if that harness inherits `release`.
The published `.11` baseline at `a51bbc06` used `dev` with incremental
compilation/debug information disabled. After its successful native/container
artifact cycle, the owner's 2026-10-01 first optimization step introduces
`package-opt1` for `.12` packaging, independent of version label. It inherits
`dev`, changes optimization to level 1 with all LTO off, and preserves 256
codegen units, assertions, overflow checks, unwinding and no stripping.
Receipts must bind those settings; neither this recipe nor higher optimization
is certified by the successful `.11` build. Local `dev`, pinned build helpers
and `bench-audit` remain unchanged.

On 2026-10-13, after the `.12` artifact cycle, the owner requests a bounded
increase in optimization. The next package recipe selects `package-opt2`,
changing only opt-level 1 to 2. Keep `package-opt1` unchanged, all LTO off,
256 codegen units, assertions/overflow checks, unwinding, no stripping, cold
build settings, one Cargo job and unoptimized helpers. Native/container paths,
guards and receipts must name the new profile; old opt1 evidence is not opt2
acceptance. This authorizes implementation, not a push, dispatch or publication.
Build reliability and same-SHA opt1/opt2 runtime/resource comparison remain
external gates; no import speedup or compiler-memory reduction is inferred.

Keep those verdicts separate: an optimized measurement cannot certify an
unoptimized package, nor can that package certify a future optimized release's
latency. Historical "release-like" measurements still mean the release profile;
do not relabel them or silently relax budgets for development-profile packages.
Retain extracted-package identity separately from local counterpart builds and
instrumented diagnostics. Keep `bench-audit` for historical ladder continuity.
Native package acceptance requires all five extracted entrypoints to report the
selected package version and full source SHA, with observed binary hashes bound
in build receipts. Plaintext package revision files and successful help commands
are insufficient; source identity observations do not authenticate build
provenance or replace independent archive verification and full-platform gates.
A controlled same-SHA, same-host native profile comparison may quantify the
recipe difference with repeated, condition-matched samples; it is not a
code-change baseline, cannot isolate LTO from other recipe differences, and
must not become a universal timing correction. Ordinary before/after code
comparisons name both revisions and hold the effective profile constant.
The [measurement-profile contract](dna_feature_rendering_latency_plan.md#measurement-profiles)
assigns this comparison to Glen without adding benchmark work to installer CI.

## DEC-040: Tutorial GUI Contracts Cannot Invent Execution Authority

Status: active

Tutorial GUI acceptance contracts may select only controls in a closed shared
semantic-control catalog. Each catalog row fixes the owning window, allowed
interaction, whether the action persists project state, whether it may
establish a scientific result, and any replacement-text policy. Merely knowing
a semantic snapshot id is not authority to type into or activate that control;
command-capable and snapshot-only controls are unavailable unless deliberately
added with a narrow typed policy.

Project persistence and scientific proof are independent. Metadata edits may
make the project dirty without claiming a scientific result. A scientific
effect requires before/after facts plus an engine/report/artifact/state
verifier; visual state never substitutes for that proof. GENtle validates the
contract and exposes read-only semantic state, while an isolated external
runner owns input delivery, evidence collection, and the final acceptance
verdict.

Candidate tutorial acceptance binds an explicitly frozen revision, helper/GUI/CLI
binary hashes, source contracts and fixture inputs. A derived coverage inventory
is not an execution verdict. Reference/manual and view-only chapters cannot
silently count as scientific GUI passes. External auditors may run the gate
outside hosted CI; private biological acceptance remains separately authorized.
Historical screenshot integrity is distinct from current-candidate acceptance.
Repairs produce new immutable run receipts, optionally referencing the failed
report by hash; they never combine passes across revisions. Publishing captures
remains an explicit reviewed action, not a side effect of a green test.

Semantic rectangles use egui logical client-area coordinates, not native
root-screen coordinates. Before emitting X11 input or annotating a retained
screenshot, the external runner must bind the owning semantic viewport to one
exact native client, sample a stable root-screen client origin, and record the
translation. Missing or ambiguous native geometry is a harness gap; it must
never be repaired with guessed title-bar offsets or classified as a product
failure.

Embedded semantic windows may share an egui/native viewport with their owning
tool window. That association is carried by the snapshot's egui viewport
identity; one semantic-window name does not imply one operating-system window.
Metadata saves may be deferred across dependent text-field/button interactions,
but must be flushed before the next scientific action. The scientific action
must then produce a fresh dirty transition before its result is saved and
verified; a pre-existing dirty flag is never evidence of completion. Child
viewports publish their own semantic frames, including exact pointer-hover
state, so the external runner can prove that the intended egui control observed
the pointer before emitting a button event.

## DEC-041: Inner-Agent Web Research Is Explicit, Public, And Audited

Status: active

GENtle's own structured catalogs are consulted before external research. Every
Agent Assistant request carries a bounded prompt-matched helper/vector catalog
projection, including relevant recent-conversation terms, so known product and
sequence identities are not needlessly guessed or rediscovered. Catalog
knowledge remains distinct from loaded project state.

Public-web research is an optional capability of a selected inner-agent
transport, not a general browser or command-execution privilege. It requires a
per-request user opt-in and may expose only narrowly allowlisted public search
and page-reading tools. Local/private/reserved network targets, credentials,
shell, filesystem, and project-state access remain unavailable. Redirects and
DNS results are checked under the same public-network policy. Because
prompt-derived search terms may leave the machine, the capability remains off
for sensitive projects; the agent contract forbids transmitting sequences,
local paths, credentials, personal data, or confidential project details.

Web source provenance is attached by the trusted adapter after actual tool
execution and is displayed separately from model prose. A model cannot certify
its own browsing by emitting provenance fields, and web content remains
untrusted evidence that cannot alter the response contract, parser validation,
confirmation policy, or nested-agent guardrail.

## DEC-042: Precomputed Genome Motif Evidence Is Optional And Scale-Separate

Status: active

Finalized sparse whole-genome motif packages are optional evidence providers,
not GENtle runtime resources. GENtle discovers and invokes an external DuckDB
executable only for an explicit query, opens only payload paths selected from
the package's exact inventory, and imposes motif, interval, file, row, output,
and wall-clock bounds. An absent package or executable therefore has no startup
cost and cannot disable local JASPAR scoring or any unrelated workflow.

Annotation-selected regulatory/TSS packages are a distinct provider, not a
weakened full-scan or TP73-conditioned cofactor contract. Their delivered
`file_inventory` owns payload paths and counts; the copied original-atlas
inventory supplies source policy only. Validate bounded metadata and selected
payloads, never scan all genome files on opening. Keep typed subset coverage,
annotation releases, full physical footprints and n:m transcript ownership in
portable reports and downstream views. A known-empty regulatory intersection,
a valid zero-row file and unsupported coverage must remain distinguishable.
Catalog inspection cannot invent unavailable score configuration. Package
hashes prove consistency, not independent reference authentication.

Package hits remain on the package-declared score and retention policy. They
must not share a numerical axis with GENtle's local scorer unless a future
contract proves cross-source calibration. In the current providers this is a
material estimator difference, not merely a naming concern: sparse JASPAR
packages declare additive per-base pseudocounts while GENtle's local smoother
uses its own finite-sample estimator, so zero-count cells can diverge sharply.
Detailed TSS exports may attach saved query reports, retaining source-file and
canonical-report hashes. Triangle orientation follows the displayed sequence,
not the sign of the score; original genomic strand remains explicit. Footprints
are clipped only for display, raw imported scores use their own labelled scale,
and an unqueried or sparsely retained interval never becomes measured zero.
Every report states per-motif
retention-floor completeness, density limiting, row truncation, and genome
compatibility. Contig aliases and geometry alone may support a
`contig_geometry_matched_only` result; they do not establish per-contig
sequence identity. Row queries require an explicit bounded motif set because
the sparse package is not an all-motif density index. Saved genomic region sets
are resolved to their canonical intervals before entering this same bounded
query contract; they do not define a parallel DuckDB route.
Their canonical references remain attached to resolved report rows and bound
into report identity. Before any motif payload query, a saved assembly
accession must exactly match the package accession; it takes precedence over
display names. Without a saved accession, its assembly name must exactly match
the package name or accession. No aliases, patch equivalence, or liftover are
inferred. Missing required identity or a mismatch rejects the complete query
as incompatible, even with a matching explicit `expected_genome_id`. This
validates declared assembly identity, not sequence bytes or annotation release.

## DEC-043: SCREEN cCRE Evidence Is Optional, Assembly-Bound, And Non-Causal

Status: active

ENCODE SCREEN Registry files are optional external evidence, not bundled
project data and not startup dependencies. GENtle accesses them only through an
explicit install, index, overlap-query, or materialization request. An absent
BED/index therefore cannot slow startup or disable local sequence, annotation,
TFBS, promoter, or cloning work.

Every supported source descriptor fixes provider release, species, taxon,
assembly aliases, subset, expected classes, URL, and BED field order. The
expected accession namespace is validated while indexing. The derived sparse
index binds the exact BED SHA-256 and byte size, and every query checks both
that content identity and the loaded sequence's genome anchor. Human and mouse
resources are never inferred from filenames or silently mixed.

A cCRE/ELS overlap is biochemical candidate-regulatory evidence. It does not
prove cell-type activity, identify a target gene, or establish causal
regulation. Materialized rows remain ordinary strandless
`regulatory_region` features carrying the source identity and these limits, so
all adapters and renderers consume the same engine-owned interpretation.

## DEC-044: Regulatory Annotation, Activity, And Primary Signal Stay Distinct

Status: active

Ensembl Regulation resources are optional external evidence and never startup
dependencies. Annotation snapshots are pinned by provider release, API and
pipeline versions, gene-annotation release, species, taxon, assembly name and
accession, then normalized and indexed under exact content hashes. Queries fail
closed on release or assembly mismatch; notably Ensembl GRCm39/mm39 and SCREEN
GRCm38/mm10 mouse resources are not interchangeable.

A regulatory-feature overlap, an associated-gene annotation, a regulatory
activity call, a primary peak, and a quantitative signal are separate claims.
GENtle must retain their independent source identities and cannot derive a
biosample-activity or causal target-gene claim from annotation overlap alone.
Future primary-signal display should use bounded region reads, explicit assay,
biosample, condition, replicate, output-type and unit provenance, and a
declared derived-activity rule rather than silently reducing signal to peaks.

## DEC-045: Portable Genomic Regions Are Assembly-Bound Evidence Ledgers

Status: active

A saved genomic region of interest has one canonical machine geometry:
explicit species when known, assembly, contig, genomic strand, and 0-based
half-open bounds. Human 1-based inclusive text, BED6, sequence-local
coordinates, and locus-figure overlays are projections of that record rather
than independent coordinate authorities. Imports never infer assembly from a
chromosome label or filename, and assembly conversion requires a future
explicit, provenance-bearing liftover operation.

Region geometry, selection purpose, and supporting evidence remain separate.
Provider annotations and CUT&RUN enrichment may explain why a span was saved,
but neither establishes affinity, causal regulation, or a target-gene verdict.
Combining intervals requires an explicit union, intersection, or hull that
retains content-bound parent identities. Canonical JSON is lossless; BED
round-trips require a sidecar binding the complete region set to the exact BED
bytes and coordinate contract.

Promoter-cofactor handoff reuses this ledger, not a parallel annotation store.
An explicit `promoter_cofactor` capture preserves the selected span and a bounded
original query `source_record`, checked against its canonical-JSON SHA-256.
This binds copied evidence to the displayed result; it does not authenticate
the producer or re-query/revalidate an unavailable source package. Raw motif
scores, selected-site presence, cohort statistics and sample/control depth
summaries stay distinct. Maxima cannot reconstruct a coverage track or establish
cofactor occupancy. Unspecified anchor/promoter orientation stays unstranded.

Saved cofactor motif hits can be attached explicitly as ordinary DNA features
through preview and digest-approved materialization. The engine verifies the
saved source row, structured catalog assembly/taxon metadata, exact contig,
full containment, verified anchor orientation and whole anchored DNA against
the exact locally prepared entry. Assembly metadata comes from the NCBI name
or the assembly suffix of the structured Ensembl file stem, never a display
label. Patch names are not silently equated. No compatible-genome fallback,
implicit download or inferred CUT&RUN coverage is allowed.
The approval binds ROI content, catalog entry, oriented reference DNA,
projection and the ordered annotation state. Mutation rechecks these inputs,
rejects duplicate attachment and uses the existing undoable feature writer.
Unstranded hits, anchors and promoters remain regions: the current simple
feature editor cannot express unknown orientation without implying a strand.
Features retain the original ROI/source report as a base64-encoded JSON qualifier
(`gentle_roi_json_base64`; remove flat-file wrapping whitespace to decode); those
checksums bind evidence but do not authenticate the producer's claims.

## DEC-046: Regulatory-Fragment Panels Are Exact Contrast Plans

Status: active

Upstream evidence-guided boundary selection is a separate read-only proposal.
Its fixed signal-measurement envelope, historical selection provenance, display
window and proposed insert are distinct geometries. Required biological context
precedes compactness; explained per-end changes retain alternatives and cannot
silently waive context. Vector/MCS constraints may motivate new proposals, not
modify approved downstream products. Raw coverage is not an enrichment claim;
optimizing insert boundaries cannot retrospectively validate them with the
same broad-window signal.

Fixed-window similarity searches use the shared core geometry operation, not
an adapter-owned window/union/intersection algorithm or the adaptive selector.
Called-peak seeds require explicit file and source-lane hashes plus declared
caller/version/parameters and sample/control provenance. Coverage labels or
overlap counts cannot promote a source to that tier. Peak retention is a
presence flag, not an independent-source count, and is not independent of the
source coverage. Neither a valid hash nor declared settings certify peak QA.

A regulatory-fragment panel consumes exact persisted genomic regions rather
than rediscovering or substituting candidate biology. Requests bind region-set
content, ROI identity and content, current sequence projection, assembly,
release, and catalog-validated vector context. Any stale projection, changed
content, missing source, or mixed assembly/release fails before candidate
generation; liftover and partner inference are separate future operations.

"Smallest panel" means a deterministic bounded cover of exact member pairs for
the requested questions, not a qualitative preference and not a Cartesian
product of geometries. Every retained member and every generated omission has
a typed reason. Sequence uniqueness, repeats, regulatory annotations, TFBS,
occupancy/chromatin, junction evidence, and cloning risk remain independent
evidence lanes. Missing or deferred evaluation is `not_evaluated`, never a
pass or a biological conclusion.

Planning is read-only and content-addressed. Approval binds the normalized
request, source/evidence/vector state, exact member order and operation
payloads, but does not validate regulatory activity. The legacy promoter-panel
materializer remains byte-compatible and single-fragment. Ordered multi-fragment
designs use a separate exact-product proposal and explicit digest approval.
The proposal binds full product sequences, annotations, output IDs and ordered
instances. Approval never authorizes silently regenerating different products.
Current inputs are rechecked; every product is committed atomically or none is.
Designed insertion-context replacement is not restriction/Gibson simulation:
final-product audits remain `not_evaluated`.
GUI product approval is distinct from panel review and is cleared on input
changes, replanning and execution attempts. The GUI submits the cached exact
proposal, never a freshly regenerated replacement under an older approval.

Sequence-derived evidence is computed inside the same read-only operation and
retained as separate content-addressed lane observations and findings. A bound
external report id and digest is provenance only: without typed report content
or a typed engine-owned report resolver, Ensembl Regulation, TFBS, CUT&RUN, and
chromatin lanes remain `not_evaluated`. The typed locus-document resolver hashes
and parses the same bounded bytes, verifies current sequence/anchor and
assembly/release, and retains separate annotation/model/occupancy payloads.
Unavailable payloads never become evaluated merely by being cited.
Rendering is a passive projection of
the authoritative plan and may not reinterpret those states.

## DEC-047: Region Homology Uses Query Coordinates And Explicit Orthology

Status: active

Genomic-region homology is an optional, read-only analysis over already
validated local genomic-DNA indexes. GENtle does not download or prepare an
index as a side effect of inspection, does not accept a cDNA database in its
place, and binds cache/report identity to database content rather than only a
path or genome label. An unavailable optional genome is unassessed evidence,
not biological absence.

Module uniqueness requires a completed, explicitly requested same-genome
comparison. A zero retained repetition percentage without that comparison is
unassessed, not a passed uniqueness check.

The query is the only display coordinate authority. Target insertions are
retained as structured evidence but never add alignment columns. Sequence
similarity does not establish orthology: only an explicit, provenance-bearing
expected locus licenses the `expected_ortholog` class. Same-genome non-self
similarity remains a repetition/ambiguity observation rather than a functional
off-target verdict.

Conserved blocks and their composition with selected CUT&RUN, motif,
annotation, or reporter spans produce inspectable experimental hypotheses.
They may motivate testing a block alone, a partner alone, their ordered
combination, and motif-disrupted controls, but they never prove that a fragment
is autonomously functional.

Paired-context conservation requires an actual shared ortholog locus with
strand-consistent block order and measured query/target gaps within the stated
bounds and gap-difference tolerance. Reuse the alignment's winning-HSP
projection for those coordinates; query proximity alone is insufficient.

## DEC-048: Command Submission And Observation Never Wait For Execution

Status: active

Every interface uses the same separation of command submission, execution and
observation. Submission acknowledges a stable identity and lifecycle state
promptly, or rejects with a structured reason. Expensive preflight, resource
resolution and computation happen after admission, not before acknowledgment.
Status/result-availability queries return the latest timestamped snapshot
without waiting for the worker or holding its execution lock. Cached status
must expose staleness; unknown progress/ETA is never invented completion.

The lifecycle must distinguish queued, validating, blocked/waiting, running,
cancellation/termination in progress, and terminal outcomes. The exact versioned
wire vocabulary belongs in the protocol implementation. A blocked command states
its reason and dependency/approval needed; it does not block the application,
agent conversation or unrelated work. Conflicting mutations retain explicit
ordering, short-lock commit and stale-result protection under DEC-026.

Asynchrony is a system responsibility, never a requirement for an agent or human
to select a special background command. Existing synchronous CLI/script behavior
may remain as a documented caller-side wait policy over the same lifecycle;
it must not block shared dispatch or keep the engine locked. An immediately
joined worker thread does not satisfy the asynchronous submission contract.

Acceptance authorizes no unapproved execution. The exact scientific request,
approval digests, result validation and atomicity remain unchanged by scheduling.
Job receipts and progress updates are execution evidence, not scientific project
edits. Cancellation is requested before it is confirmed. Worker ownership,
shutdown behavior and interruption must be explicit; no host-independent
execution or restart recovery is implied merely by issuing a job ID.

This decision records the invariant, not completed runtime migration. Existing
blocking paths and deterministic acceptance tests remain tracked in the roadmap.

## Measured-Gel Evidence Invariants

Status: active

Measured-image sizing and virtual gel simulation are distinct evidence paths.
Original image bytes remain unchanged and content-bound; display previews and
future contrast/geometry transforms never redefine original measurement pixels
silently. Confirmed band assignments and calibration settings belong to the
engine report, not a GUI session. DNA fragment length (bp), apparent protein
mass (kDa), and eventual intensity/abundance measurements are different quantities.
Reference/sample geometry uses one explicit image coordinate system and
migration direction. Missing reference bands, absent registration, and values
outside calibration must not become guessed measurements. A size match is not
proof of molecular identity, specificity, concentration or biological function.

Incomplete manual assignments are engine-owned, digest-bound **drafts**, separate
from validated immutable reports. Persisting a draft is not a sizing claim; analysis
must validate it again. GUI zoom/pan maps pixel centers, never thumbnail pixels.
Heavy image work runs on detached snapshots; only the live host commits a current,
uncanceled result. File exports remain explicit side effects and are not undone
by closing the editor or switching projects.
Draft-only undo checkpoints retain the draft store and journal, not a copy of
unrelated DNA state; these edits advance mutation but not structural revision.
Draft request allocations are shared across journal copies as well as project
snapshots; the serialized contract remains ordinary self-contained JSON.

Measured profiles retain an explicit evidence kind, integration measure, units,
coverage and source digest. Image display intensity is not calibrated detector
fluorescence, mass or molecule count. Record fixed stored-channel luminance and
native-alpha compositing explicitly,
without silent color-space/transfer conversion or clipping invalid channels.
Rebin integrated weights conservatively; null intervals/markers never become
zero or interpolation bridges. Shape
comparisons normalize independently on explicitly shared valid support and
record that conditioning. Weighted QQ/KS-distance/Wasserstein summaries do not
justify significance p-values from correlated pixels/bins. Inferential tests
require a separate sampling/replicate contract. Linear algebra keeps explicit
coefficients, missing operands and signed residuals; it establishes neither
gene identity nor additive reaction chemistry. Native encrypted measurement
members remain unsupported, even when a thumbnail is explicitly admitted.
Exports project frozen results; they never silently rescore or re-extract data.

Signed gel algebra v1 freezes its calculation, profile/request digest field
order and receipt bytes. Future algorithm/serialization changes require an
explicit new version rather than reinterpreting saved evidence. Verification
replays saved operands with bounded, memoized dependency traversal and compares
the entire derived profile, not only a supplied digest. This proves internal
consistency, not source authenticity. Unverifiable legacy results stay readable
but cannot be report-exported; new structured imports cannot manufacture signed
residuals. Derived identity uses the output ID and algebra input digest, not
inherited run/lane metadata. Show all operand calibrations and warnings.

The signed-density display divides each saved bin integral by its size width,
labels the resulting units per bp/kDa and preserves its signed area and gaps.
It does not change the integration measure or establish fluorescence, mass or
molar density. No implicit normalization, fitted/reference scaling, clipping
negative values, residual probability summaries or biological mixture claim is
permitted. Algebra renderers consume the engine-prepared presentation. Ordinary
comparison renderers consume the replay-verified frozen `GelComparisonReport`;
density (integral / bin width), like the CDF, is a pure display projection, not a
new result model. GUI replay/rasterization occurs only in a background inspection
bound to project revision, selected result and display mode. Omitted export mode
keeps legacy bin-integral pages; the new GUI session default is density. A common
numeric axis does not establish matched loading/exposure/calibration; no curve
is silently fitted to its own maximum. Display mode is PDF-only and requires a
selected comparison; it is rejected on JSON rather than ignored. Comparison
replay currently compares the entire v1 report, including warning prose; changing
that algorithm/prose requires an explicit compatibility policy, not silent
acceptance of older results under changed semantics.
Coefficient-weighted operand overlays are engine-prepared on the residual grid
and retain each nonzero term's identity, order and own gaps. They share one
axis with the residual and never imply that missing residual coverage is known.
Zero terms remain provenance/grid-only; repeated IDs are not collapsed.
This additive presentation changes no frozen algebra-v1 bytes, coefficients,
comparison pages or biological interpretation. Owner approval, 2026-10-20:
an explicitly selected non-residual custom-only reference may occupy a separate,
unscaled native-bin density panel with its own labelled Y axis and provenance.
It is not an operand or proven priming identity. Size kind/range and retained
sources must validate; no units conversion, rebinning, amplitude matching or
display multiplier is implicit. Pair export with exactly one selected residual;
reference selection participates in worker/cache identity. Real-data subtraction
coefficients and background/marker policy remain undecided. Fitted/manual
display multipliers still require separate approval.

### Cross-Ladder Residual Policy

Owner decision, 2026-10-25: residuals across different ladders are permitted.
Operands may retain distinct ladder calibrations; matching ladder IDs are not
required. Use each profile's recorded size axis and retain both calibrations,
source identities, coverage gaps and warnings. Existing size/evidence/unit/
integration-measure compatibility and source-verification rules still apply.
This permits cross-ladder residuals; it does not assert calibration equivalence,
align uncertain axes, choose a scale factor or validate a biological mixture.
Existing calibration/registration warnings remain applicable, not a prohibition
on different ladder IDs.

Real-data coefficients, background and marker masks remain undecided. In
particular, this decision does not approve the synthetic lesson's `1/-1` for
the cell-line data. Native/platform and scientific acceptance remain separate;
the retained public priming comparisons still perform no subtraction. No
automatic fitting, new quantitative yield claim or reprocessing is authorized.

## Advisory Cloning Time And Cost Boundaries

Status: active

Bench-planning allowances are separate from computational cloning macros.
`gentle.cloning_wetlab_steps.v1` initially admits only hand-authored synthetic
defaults for `restriction.digest_ligate_extract_sticky`, not laboratory
measurements or enzyme-specific protocols. The output-only estimate v2 keeps
unknown breakdowns as optional ranges; a missing model never means zero work.
Other routines retain explicitly labelled legacy elapsed scalars.

Hands-on and unattended work sum separately; elapsed minimum/typical/maximum
scenarios use their own dependency-only critical paths. Existing summed queue
and purchasing delays form a disclosed prefix. This is idealized concurrency,
not staffing, equipment-capacity or calendar scheduling, and no success,
yield or actual completion date follows. Invalid IDs, ranges and dependency
cycles fail explicitly rather than dropping steps or fabricating a breakdown.

The retained cost total preserves its old accumulation and is decomposed into
routine base, stock-drawn material, missing-material purchasing proxy and
machine run. Machine charges retain the legacy base-duration approximation.
Material values are known-unit subtotals, not package/quantity-aware stock
accounting or quotes. Unknown/unpriced classes remain visible and never
establish zero cost. All surfaces consume the same estimate and explanation;
ranking and guardrails remain advisory, not execution or purchase approval.

Bench-model presence and duration must not change policy outcomes. Preserve
the original family/base scalar plus summed queues/purchasing delays as the
time input to `composite_meta_score` and `max_time_hours`, even when the
displayed elapsed estimate has a critical-path basis. Expose this separately
as `explanation.time_policy` with hours, `legacy_family_scalar` basis and
the two consumers. Historical reports missing this record remain unknown;
do not reconstruct their scoring basis from a newer model or their displayed
elapsed value. Guardrail passes do not bound the displayed scenario ranges
or actual laboratory completion. Changing calibrated cross-method ranking is
a separately approved policy change, not a reward for richer documentation.

Protocol Overview is a read-only projection of the recorded bench schedule,
not a scheduler, executed lineage, validated recipe or material-flow model.
Its shared JSON/SVG binds the exact routine identity and estimate; consumers
verify that binding rather than trusting a supplied content digest. Columns
encode dependency depth, not time. Critical-path highlighting uses only the
recorded typical path and stays unavailable if absent. A diagram neither
changes ranking/guardrails nor establishes staffing, equipment capacity,
reagent/sequence identity, physical execution or readiness. Do not generate
bench steps from computational macro scripts or populate molecular protocol
cartoons without explicit molecule bindings.

Supplier discovery is independent of both these synthetic estimates and
implemented service adapters. A country-coded readable HTML/PDF protocol is a
dated document reference, not structured executable instructions, commercial
eligibility, sequence identity or validated laboratory performance. Global
`catalog_only` supplier rows remain ineligible for quote/submission even if
capability flags are supplied; listing does not contact suppliers or download
protocols. The named company's country is not inferred shipping coverage or
ultimate ownership. Unknown metadata remains absent. Promega's worked-example
status is an owner preference, never a ranking or automatic routing rule.

## Catalog Evidence Lanes And Local Update Authority

Status: active

Extend existing catalogs, not a second planning system. Biological product
identity, commercial offers, immutable protocol documents, reviewed normalized
procedures and private local observations are separate evidence lanes. Names,
aliases and readable protocol links never establish cross-lane applicability,
source authority, orderability or review. Machine-extracted procedures remain
drafts until explicit review against exact retained document passages.

The first guarded import supports existing provider configs only. Preview
binds the exact source set/bytes and scope/destination; apply requires explicit
approval and fresh validation under a nonblocking writer lock. Full-record
replacement or explicit withdrawal tombstones prevent ambiguous merges and
accidental fallback to lower overlays. Prior bytes and conditional recovery
receipts precede a single atomic catalog replacement; never claim multi-file
atomicity or replay an inbox automatically. Hashes bind content, not vendor
authenticity or signatures. Local stock/prices/equipment and saved biological
reports are not overwritten, and catalog imports do not alter planning ranking,
readiness or commercial authorization. See [the contract](catalog_updates.md).

## Labware Reference And Resource Boundaries

Status: active

Rack grid admission is shared across engine and GUI processing. Dimensions
must be positive, multiply without overflow, and fit the 1,536-position
engineering processing budget. Apply admission before grid expansion,
placement/reflow or export, never as a saved-project deserialization rule.
Unsupported historical racks remain intact and visible as unsupported; they
must not be silently truncated, processed as an empty valid grid, or exported
over an existing artifact. Increasing this budget requires resource evidence,
not a claim that a physical carrier holds that many tubes.

Tube reference profiles are a static, source-backed shared contract, not new
physical inventory records. Retain raw drawing digest/revision/page and mark
each dimension as drawing-backed, approximate or unknown. Side-view cap width
does not establish the plan-view hinge/cap outline. Visualization-only tube
models must disclose omitted cap/hinge/cavity and interpolated geometry; they
never imply fit, manufacturability, material or biological qualification.
Existing PCR/plate templates are not a compatibility authority for Safe-Lock
tubes. Owner decisions recorded 2026-10-19: offer both collar-seated and
floor-seated 1.5/2.0 ml carriers; use a conservative 20-mm approximate
side-view pitch without claiming measured cap clearance. The 96-position
intent is small 0.2-ml PCR tubes at 9-mm pitch, not 0.5-ml Safe-Lock tubes;
new PCR dimensions require an independently retained official drawing reference.
That gate was met on 2026-10-20: the [PCR drawing reference](../data/labware/README.md#pcr-02-ml-drawing-and-96-position-carrier)
binds raw PDF bytes, revision and page 4. The explicit floor-seated 8x12 carrier
uses nominal plate-format coordinates without claiming a plate, thermocycler
block, thermal contact or commercial-holder fit. Missing collar/cap dimensions
remain unknown; the visualization renderer must refuse incomplete dimensions
without panic or file overwrite. Standard coordinates do not qualify hardware.
Owner clarification, 2026-10-20: offer adjustable compact spacing below 20 mm
as a separate extension, preserving the existing templates and their numbers.
Storage is a distinct use from treatment/bench handling: prioritize dense
transport and freezer organization when samples are not being treated in the
rack. Commercial paper/cardboard boxes are possible storage destinations, not
necessarily objects to fabricate. A future box profile needs source/measured
interior, grid/divider, closed-cap-height and lid constraints; a nominal tube
volume or cardboard material never establishes fit or temperature suitability.
Compact pitch is now an explicit rack-scoped whole-micrometre override, omitted
when absent. A shared engine operation and pure geometry admission preserve
the original 20-mm defaults and outer margins. The setter admits collar-seated
spacing; each export revalidates its selected seating mode. Unsupported templates
refuse the override, never silently ignore it. Named compact web/collar allowances
are engineering assumptions, not cap-clearance, strength or material evidence.
Storage-box formats must remain variable: independent row/column grids, interior
dimensions, dividers, height and lid, not a hardcoded 9x9 box. Existing generic
grid/label operations support different formats but do not constitute a physical
box model or fit qualification. Owner clarification, 2026-10-20: 1.5/2.0 ml
tubes rest on the box bottom; the openings are square **when viewed from above**,
with snug tube-to-divider contact. Adjacent cells are separated only by divider
wall thickness, not a requirement for literal tube-to-tube contact. This does
not establish cap/lid clearance, dimensional tolerance or freezer suitability.
The [Semadeni example and missing fit evidence](../data/labware/README.md#different-storage-box-formats)
are recorded as a source lead, not a qualified compatibility profile.
Owner-authorized derived-box slice, 2026-10-25: `safe_lock_storage_box_derived`
is a GENtle-derived design, not a commercial-box profile. Preserve variable
grids, floor-supported tubes and square top-view cells. Body-snug is a design
target: 9.85-mm body plus 0.15-mm radial clearance gives a 10.15-mm cell.
Choose 1.2-mm dividers as an engineering assumption above the 1.15-mm minimum
needed for the existing 0.5-mm collar gap; pitch is 11.35 mm. Dividers stop at
the existing floor-seated depth, below the collar. Walls/floor, approximate
closed height, lid clearances and skirt are explicit assumptions, not fit or
temperature evidence. Caps/hinges may collide. The SCAD projects box, insert
and separate lid; SVG markers remain round with an explicit visible disclosure.
Commercial-box profiles, square SVG cells and measured fit remain unimplemented.
Saved compact-carrier pitch is refused for this separate template, never ignored.
New template metadata is
optional and omitted for legacy templates so their serialized bytes and
geometry remain unchanged. No new template establishes verified fit.

Owner request, 2026-10-22: model contiguous cut 0.2-ml PCR strips of one through
eight tubes. BRAND 781320 is a separate manufacturer-bound reference, not a
renaming of the Eppendorf tube or PCR96 carrier. Its additive strip layout and
optional export count use the existing tube inspection/export operations.
The dedicated strip silhouette consumes the documented envelope and pitch;
undimensioned taper transition, collar position/thickness and bridge details
are explicit visualization assumptions, not filled-in drawing dimensions.
Keep those unknown profile fields null and retain individual-tube refusal
rules. Cut count does not change per-tube nominal capacity, create inventory,
qualify fit, or establish manufacturer support for cut-end geometry.
