# Observed Transcript Structures S2: SAM And Annotation Preparation

Status: **S2.0a test-only guards implemented with scoped macOS verification**.
S2.0b and later steps remain proposals, not authorized for implementation.
Prepared 2026-10-27 at local `main` = `gentle_rs_2_main` =
`1804605f1ba42369f943213414d79fd9a40013cb`. No fetch or remote-tip claim.
Tracked files were clean; unrelated untracked `output/` is untouched.
The [roadmap](roadmap.md) owns priority/release scope. This is a bounded expansion
of [the umbrella plan's S2](observed_transcript_structures_plan.md#slice-boundaries-and-remaining-work),
not a new release gate. The [original Claude prompt](observed_transcript_structures_s2_claude_review_prompt.md)
was reviewed; the reconciliation below is Codex's, not a fresh Claude sign-off.
The later S2.0a authorization does not authorize a production refactor or SAM
ingress. Focused local Rust tests/check were separately authorized by the owner.

Integration context, 2026-10-28: the owner requested committing this proposal
and rebasing onto local `main`, now `7329f00f607c82c195ea13d931f83b0f1bc93150`.
Its two additional commits introduce the Plan-only destination writer and
reorganize Help/manual documentation. They do not implement observed-transcript
S2 or change its kernel/Phase 2 paths. Preserve the original inspection baseline
above; recheck the actual checkout before review. No fetch or push is requested.

Review reconciliation, 2026-10-28: Claude inspected main
`550e6e0ad67fbb301e9a10bb4cc996924db3388f`; Codex checked the relevant symbols
at worktree `6878d6083536d5a260feeff97913a41614654b1f`. Git comparison confirms
the S2 documents, observed kernel and payload paths are byte-identical between
those revisions. Behavioral verification is source inspection, not Rust
execution, ZIP/BAM or biological replay. Fast documentation validation is
recorded separately below; prior execution receipts remain historical.

## Purpose In Biological Terms

Given an explicitly selected gene, a saved alignment and independently supplied
annotation snapshots, ask which splice patterns the supplied query units show
and which reference models remain compatible. The new input is **SAM text
prepared outside GENtle**, not a BAM opened or an aligner launched by GENtle.
An explicit preparation step converts source annotations to the existing exon
model. Neither step chooses a gene, canonical transcript or protein for the user.

`reference_gene = {source_id, gene_id}` stays mandatory. Additional sources use
explicit `comparison_genes`; equal symbols are not cross-source gene identity.
There is no TP53, human, assembly, chromosome, strand or MAPQ60 default. TP53 is
one optional historical replay, not the target of the feature.

## Baseline And Constraints Found In Source

Use symbols rather than line-number anchors when reviewing a newer checkout:

| Existing seam | Consequence for S2 |
| --- | --- |
| `src/observed_transcript_structures.rs`: `Input`, `ReadRecord`, `QueryEvidence`, `AnnotationSource` | Original-record identity and payload-specific inputs must not silently become SAM QNAME identity. The pure Rust API is experimental, but S1 serializes its result. |
| `src/observed_transcript_structures/kernel.rs`: `observe`, `summarize_with_control` | The context digest binds serialized `Input`; observations contain `payload_query_length`. Preserve S1 digests and field meanings, not just comparable final counts. |
| `src/engine/protocol/primer_payload.rs`: `PrimerPayloadSegment` | `reported_identity_per_mille` is mandatory. Ordinary SAM `M` without supporting tags does not provide that value. Do not invent a PAF identity to enter the kernel. |
| `src/observed_transcript_structures/{geometry,dag}.rs`: `walk`, `best_equivalent`, `Dag::paths` | One existing `MIDN=X` walker and legal-path hypothesis calculation. SAM clipping/padding normalization belongs before them; do not broaden PAF's grammar. |
| `src/observed_transcript_structures/comparison.rs`: `sources` | Requires nonempty annotation release, unique transcript IDs, exact reference binding and source-specific gene selectors. Missing metadata needs an explicit preparation policy. |
| `src/engine/analysis/rna_panel/primer_payload/reference.rs`: `catalog` | Reusable annotation shape/validation; target lengths are required and cannot be inferred from the largest exon coordinate. |
| `src/genomes.rs`: `parse_annotation_gene_records_with_progress`, `parse_annotation_attributes` | Existing permissive gene indexing can skip malformed rows and synthesize a missing gene ID. It is not an exact transcript-normalization boundary. Reuse only suitable lexical helpers, without changing its current consumers. |
| `src/engine/protocol/observed_transcripts.rs`, `src/engine/analysis/rna_panel/primer_payload/observed.rs` | S1 consumes bound Phase 2 artifacts, not SAM. Keep its strict v1 requests, reports, TSVs, accounting and inspection intact. |
| `src/engine/operation_lifecycle.rs`, CLI/MCP/managed-operation paths | Named constructor/no-save/approval/publication policies do not automatically cover new operations. Migrate and test new names explicitly. |

S0/S1's scoped macOS execution remains recorded in the
[existing ledger](observed_transcript_structures_plan.md#s1-commit-and-rebase).
Those results are not evidence that S2 exists, or native/platform acceptance.

## Review Reconciliation

The original plan at `6878d608` proposed one input-side private seam, annotation
preparation, SAM intake and adapter integration. Claude validated its central
source claims, but exposed output-identity coupling and overbroad alternative/
unavailability wording. The current sections below incorporate these changes;
the original prompt is historical, not a second competing implementation plan.

| Review item | Codex finding and revised contract |
| --- | --- |
| F1: input-only seam is insufficient | Accept. `observe`, its digest tuple and `comparison::{sources,qualify,compare}` also couple to legacy types. Separate private segment/record/context/observation views and fixed identity policies, with old and SAM public projections. Do not commit to a public generic API merely because generics were suggested. |
| F2: default output is not definite | Accept with an axis distinction. Unknown search completeness suppresses definite mapping support/reconciliation; missing masks independently make flank qualification unassessed. Masks alone do not control mapping reconciliation. State the combined defaults' result explicitly. |
| F3: unavailable unit must lose all support segments | Accept. Drop every segment of an unsupported SAM unit from interpretation and register its private key as unassessed. S1 does this for unknown target length; its other missing-view cases can retain possible support. Do not claim S1 universally clears every unassessed view. |
| F4: duplicate accessions | Accept a typed `duplicate_transcript_accession` source refusal within the declared annotation scope. No automatic ID qualification, contig pruning or source merge. The reference dictionary must still cover alignment competitors outside annotation scope. |
| F5: touching exons and small reference gaps | Accept the structural risk: `compare` uses every exon boundary gap as a model skip, including zero-length gaps. Refuse touching/overlapping models in preparation. Do not infer that every short positive gap is an artifact or select a universal biological cutoff. Require an explicit gap-admission policy before real-data preparation. |
| F6: identifier/version mismatch | Accept. Keep source IDs exactly as written; retain separate version attributes in receipts. Require exact selector agreement with the prepared receipt; a spelling/version mismatch is not `missing_gene`. |
| F7: only maximal-support paths contribute | Accept as inherited S0 policy, not a new scoring defect. Retain all admitted observations and mark optimal-path membership; only best-equivalent supported-base paths contribute families. No change to weights or legacy bytes. |
| F8: raw family counts ignore flank qualification | Accept for new SAM presentation only. Split unique-unit support by qualification and retain independent mapping-certainty fields. An unqualified, unannotated skip is not a candidate novel biological junction. |
| F9: identical chain keys do not prove reference verification | Accept. Carry `coordinate_authority` plus its evidence; equal keys or joins cannot promote declared authority to verified. No implicit FASTA verification or replay. |
| Low: unsupported CIGAR/axis/scope | Admit genome scope only; make valid-but-unsupported geometry and overlong queries whole-unit unavailable before kernel admission. Unmapped records without an axis remain accounted/unverified, not length zero. |
| Low: first-slice tag scope | Narrow to consumed fields and `CG` detection. No ignored-tag inventory or new NM/MD/AS interpretation/retention. Validate bounded optional-field syntax without claiming to understand every tag. |

Retain the correct clip example, `D`/`N` distinction and MAPQ 0/255 separation.
Reject placeholder PAF identity. The review's pseudogene and tiny-gap examples
motivate synthetic software tests, not verified claims about the owner's reads
or particular annotation releases. No new assay threshold is approved.

## Reconciled Implementation Order

### S2.0a: Test-Only Preservation Guards First

The owner-authorized first implementation task is tests only. Freeze exact
S1 report/data/TSV bytes, `Summary`/`Observation`/
`RecordSupport`/`ChainFamily` field order, context and observation/chain digests,
`comparison_work`, covered error strings, query/node sort keys and the PAF
stored-support check. Retain all ten existing Phase 2 serialized-row/work cases.
Include the serialized effect of every PAF segment field, even fields the
geometry algorithm does not read, and the existing defaulted
`Input.unassessed_records`. No production changes or serde-default additions to
S1 v1. Bound any deterministic output paths before pinning file bytes; no
expected-byte regeneration to conceal changes during S2.0b.

### S2.0a Implementation And Verification, 2026-10-28

The reconciled plan was committed as `36702cce` and rebased onto local main
`550e6e0a` as `e84b3f73` before test work. The owner separately authorized local
Rust tests/check and the minimal missing test-path correction. Tests were
verified as a working diff; the owner subsequently requested committing and
rebasing it. Main first advanced to documentation-only `962b68a2`;
the results below belong to `e84b3f73` plus the fingerprinted Rust diff, not a
new-source/full-suite acceptance claim. The later integration target
`b0521123447d3636258cd6fca1702968c1ddb190` adds lifecycle-preview publication,
not observed-transcript kernel/payload changes. Recheck the integrated source;
do not silently relabel the pre-commit binary receipt below.

- `src/observed_transcript_structures/preservation_tests.rs` freezes five
  in-memory scenarios: exact forward, reverse/masked, unassessed/masks unknown,
  competing loci and distinct records sharing a name. Pin input/result compact
  bytes, pretty result bytes, individual observation/record/family arrays,
  context/observation/chain identities and work counts `54,51,24,71,142`.
  Separate tests preserve query/node ordering, all 14 segment fields in the
  serialized context, the existing defaulted `unassessed_records`, ten error
  strings and the PAF stored-support check on both orientations.
- `src/engine/analysis/rna_panel/primer_payload/observed_preservation_tests.rs`
  replays the old v1 request/report/data and five JSON-cell table recipes
  independently against actual S1 publication. Completed and incomplete cases
  retain actual temporary paths and every derived binding, not scrubbed hashes.
  Geometry uses the separately frozen kernel, not a second interpreter. Exact
  compact/pretty JSON, row quoting/order, receipts, nonclaims, absent counts,
  mandatory fields and unknown-field refusal are checked. Three request/source
  errors are pinned too. Fixtures remain inline and explicitly synthetic.
- Only test files/registration and the owner-approved one-line
  `#[path = "sequence_design_preview/tests.rs"]` under `#[cfg(test)]` change
  Rust. No runtime code, production serde/schema, algorithm, dependency or old
  ten-case Phase 2 oracle changes. There is no expectation-regeneration mode.

Keep the primary failures separate: initial test compilation and whole-project
fmt stopped at main's missing test-module path. After the approved path repair,
the first baseline run passed four tests and deliberately failed two unsealed
capture assertions. Their literal values were then pinned without runtime
changes; the final six pass. Slice 4 initially passed 13/14 in the sandbox:
`UnixListener::bind` at its special-file fixture returned `EPERM`. The unchanged
suite passes 14/14 on the approved host rerun. No assertion, timeout or safety
check was weakened.

All Rust commands below use `CARGO_INCREMENTAL=0`, `--locked --offline
--no-default-features --lib -j1`; tests use `-- --test-threads=1`:

| Command/filter | Result |
| --- | --- |
| `cargo test ... s2a_` | 6/6, no ignored tests |
| `cargo test ... observed_transcript_structures` | 28/28 |
| `cargo test ... observed_transcripts_` | 16/16 |
| `cargo test ... primer_payload` | 48/48, including the ten-case serialized-row/work oracle |
| `cargo test ... sequence_design_preview::tests` | 14/14 on host; retain the sandbox failure above |
| `cargo check -q ...` | Pass; existing headless dead-code warnings |

Filters overlap and must not be summed. The earlier `primer_payload_` filter
passed 35 tests but is only a subset; use the full-module 48-test result.
After confirming no active Cargo/rustc process, the owner-approved cleanup
removed 9.1 GiB of inactive incremental cache, preserving sources, test binaries
and `output/`. No new cache/profile setting is committed. This is not a build-
performance or runtime-performance verdict.

Execution identity: macOS/aarch64, Cargo `1.100.0-beta.5` (`7552a885e`), rustc
`1.100.0-beta.5` (`8ea38001b9da1d6683ac5ebbc993cac89d3a09df`), LLVM 23.1.3.
Unoptimized test profile: root debuginfo 1, dependencies debuginfo 0, 256 codegen
units, incremental disabled by environment, default assertion/unwind behavior;
`.cargo/config.toml` retains the 16-MiB test-thread stack.
Cargo.lock SHA-256 remains
`40ae5b6c9ba74ba080d39be11b4e57ad2dbc020dba5ecad7222e22f639f84918`.
The 20-path Rust/Cargo/config source manifest at
`/tmp/gentle-s2a-source-manifest.sha256` has SHA-256
`58e1d6d328af6b23b1da3f96f31c33db366fb2bea172219cb8d037d1c0b3d5f3`.
Test binary `target/debug/build/GENtle/002e669310863d66/out/gentle-002e669310863d66`
has SHA-256
`e79cffc89258bad13bce505498d808d94c396bfe046b70951315d5264952cf60`.
Logs are local `/tmp/gentle-s2a-*.log`, not committed fixtures or release assets.

Scoped rustfmt for all five changed Rust paths passes. Whole-workspace
`cargo fmt --check` now resolves the test module but fails on 12 unrelated
pre-existing Help/manual-adaptation files; leave them untouched. Fast checks
include documentation status, 18 LF/CRLF documentation cases, and 16 observed/
payload Python contract cases (15 passed, one binary-dependent skip).
Full suite, default-feature GUI, actual CLI/VM processes, Windows/Linux/native,
real-data/scientific and package acceptance are unrun. S2.0b remains separately
gated, not authorized by passing these guards. No push or Glen dispatch.

### S2.0b: Private Core And Pure SAM Geometry

After those guards execute successfully under separately authorized verification,
introduce a private format-neutral seam across both inputs **and outputs**:

- A minimal segment view with an optional stored-support assertion. PAF supplies
  its existing assertion; SAM derives support from normalized operations rather
  than fabricating a persisted PAF support value.
- A private typed record key and observation/support representation, allowing
  legacy original-record and SAM declared-unit identities without conflation.
  A private generic or discriminated key is an implementation choice, not a new
  public API or wire type. Preserve legacy ordering exactly.
- Fixed private identity policies: preserve PAF context and observation tuple
  recipes verbatim; use a distinct SAM context/observation domain. Callers must
  not select arbitrary digest domains through request JSON. The physical chain
  recipe stays shared and unchanged.
- A minimal comparison context for reference, selectors, sources and budgets,
  and shared geometry/qualification views rather than mandatory legacy `Input`
  and `Observation`. Preserve comparison call order, errors and charged work.

Project old PAF results back to their existing public records unchanged; the
SAM public projection uses `unit` and `query_axis_length`. Reuse the **same**
observation, path, reconciliation and comparison logic. Add a pure private
normalizer of already decoded in-memory SAM fields with a basewise oracle and
synthetic end-to-end kernel cases. No SAM file parser/I/O or new operation yet.
Do not create a second biological classifier or convert SAM into pretend PAF.
The private view carries only required geometry/support/context, not mandatory
PAF identity statistics. Unknown identity stays absent in the SAM presentation.

Keep old S0 input serialization and the legacy context-digest calculation at
their current boundary. Do not add serde defaults to Phase 2 or reinterpret
existing public fields. Preserve S1 v1 report/data/table bytes and identifiers,
as well as Phase 2 decisions, serialized rows, work counters and outer results.
SAM gets a separate versioned request/report/data envelope and digest domain.
Its public records use `query_axis_length` and a typed identity scope, not
`payload_query_length` or "original read" for undeclared provenance.

Physical `gentle.observed_junction_chain.v1` keys keep the same recipe: declared
reference FASTA digest, assembly, exact contig and ascending skipped intervals.
Equivalent SAM/PAF geometry with identical coordinate authority must agree on
that physical key; ingress/context/observation identities may legitimately differ.
Changing the reference gene must not change the same physical chain key.

### S2.1: Explicit Annotation Preparation

Propose `PrepareObservedTranscriptAnnotation` with a strict versioned request
and typed preparation receipt. It writes a new, bound
`gentle.primer_payload_annotation.v1` plus provenance/diagnostics; **no third
transcript geometry schema**. One request prepares one named source; subsequent
comparison consumes independently prepared sources without rerunning alignments.

Require raw input path/hash/byte count, explicit format, source/provider,
assembly, declared reference FASTA hash, bound contig dictionary, stable gene-ID
selection, input-scope declaration, budgets and fresh output directory. The
dictionary is a hashed exact-name/length document; reconcile it with SAM `@SQ`.
Do not guess aliases, liftover, assembly or contig lengths. Whole-FASTA SHA-256,
per-contig digests and SAM `M5` are different bindings, not interchangeable.
Default authority is **producer-declared reference identity**: verifying the
dictionary's bytes does not verify the original genome FASTA or aligned bases.
Missing reference identity refuses interpretable output; missing original FASTA
bytes remains a disclosed verification limit, as in S1, not a fabricated replay.

Recommended order is plain-text **GTF first**, then UCSC genePred in a separate
commit after its gap/crosswalk policy and tabular companion are ready. This is
a Codex recommendation, not an owner-approved format decision. Select the
format/column layout explicitly (including any genePred `bin` column); no GFF3,
BED12 or JSON sniffing. The archive's UCSC JSON needs an owner-prepared tabular
companion with extraction/hash provenance; do not promise automatic loading of
every archive member.

- For GTF, normalize exon rows from one-based inclusive to zero-based half-open
  geometry; require consistent source gene/transcript IDs, contig and strand.
  Preserve gene/transcript IDs exactly as written, without adding/removing
  version suffixes or changing case. Record separate version attributes in the
  receipt and reject contradictory claims. Use strict GTF attribute handling:
  reject conflicting/duplicate consumed keys rather than lowercasing them or
  choosing the first; do not reuse the permissive GFF/GTF attribute parser as-is.
- For genePred, validate declared columns, exon count, ordered arrays and bounds.
  A `name2` label is not a stable gene-ID crosswalk. Require an explicit,
  snapshot-bound transcript-to-gene-ID mapping when the source lacks gene IDs;
  never assign every row to the requested gene from a symbol or interval.
- Scan/hash the full supplied input to EOF. Selection retains **whole transcript
  exon lists**, not exons clipped to a locus window. Record supplied-file and
  selection coverage; a locus excerpt is not a whole-genome annotation search.
  Non-exon GTF features are counted as intentionally unused, not silently lost.
- Malformed rows, contradictory repeated IDs, mixed loci/strands for the same
  transcript, unsupported overlaps and missing selected-gene links refuse a
  completed annotation rather than silently dropping models. Identical duplicate
  rows may be deduplicated only under a recorded deterministic rule and count.
  Bind the annotation's included contig set separately from the full reference
  dictionary. A selected gene's repeated transcript accession inside that scope
  produces typed `duplicate_transcript_accession` source refusal: no merged
  models, fabricated qualified IDs or reactive contig pruning. This rule is not
  limited to alternate contigs. Missing competing-source gene selection remains
  a typed resolution gap at comparison; keep other alignment contigs accounted.
- A genuinely absent selected gene is a typed `missing_gene` preparation result,
  not a malformed row or an empty completed model set. Retain its source receipt
  as unavailable comparison evidence; do not silently drop that requested source.
- Preserve provider/source separation. A UCSC `ncbiRefSeq` track snapshot is not
  a direct NCBI GCF annotation release. Bind retrieval/source metadata and notices
  when supplied; no online retrieval, source replacement or inferred release.

The preparation receipt binds each exact source gene selector; S2 ingress
requires its `reference_gene`/`comparison_genes` to agree. A mismatched ID form
is `selector_binding_mismatch`, not a biologically absent gene or an opportunity
for alias fallback. A version-attribute map is provenance, not selector rewriting.

Refuse touching or overlapping exons: the current comparer would turn touching
exons into a zero-length model skip. For positive gaps, require an explicit
source-bound `minimum_reference_gap_bases` admission policy in the request/
receipt, with no default claimed as a biological minimum. A retained model below
that bound makes the selected source/gene typed `reference_gap_policy_refused`;
do not drop the model, snap boundaries or merge exons. Suspected alignment-gap
semantics also stay unresolved rather than being reinterpreted. Synthetic tests
use clearly labelled hand-chosen values. The real-data minimum and acceptability
of refusal remain owner decisions; no new RNA mechanism follows from a short gap.

**Reconciled snapshot rule:** for every S2-prepared annotation, use
`release = snapshot-sha256:<64-hex-raw-input-digest>`, even when a provider release
is supplied. Put that verbatim provider claim in nullable
`declared_provider_release` with its provenance in the receipt. Thus the reused
field has one explicit snapshot meaning in S2; old S1 annotations/semantics are
not rewritten. Require the receipt at S2 ingress and label the version kind as
snapshot, not provider release. Selection/crosswalk/gap-policy hashes also bind
the prepared artifact. Missing provider release limits claims without creating
a permanent activation blocker or a guessed annotation release.

### S2.2: Strict, Bounded SAM Intake And Query Accounting

Propose a separate `SummarizeObservedTranscriptStructuresFromSam` request, not
optional SAM fields grafted onto S1 v1. Required inputs are:

| Binding | Required meaning |
| --- | --- |
| SAM | Raw file path/hash/bytes; explicit SAM-text/unpadded-reference admission policy. No compressed/binary format autodetection. |
| Conversion provenance | Direct-SAM origin, or input BAM digest, producer tool/version, exact conversion/filter/selection command and scope. Producer declarations are not commands to execute or proof of the original BAM's contents. |
| Query-unit manifest | Hash-bound unit IDs plus every alignment-row ordinal assigned exactly once, matching QNAME/read-group metadata and declared unclipped axis length. Headers are not alignment ordinals. Missing units are explicit. |
| Coordinate authority | Assembly/reference digest and exact contig dictionary; discrepancies with supplied header names/lengths refuse the input. Missing `@SQ` refuses this initial subset, although other SAM consumers may allow it. |
| Annotation sources | Prepared annotation files and receipts from S2.1, or an explicitly verified equivalent normalized binding. Preserve the exact snapshot and selector of each source. |
| Gene/policy/control | Required `reference_gene`, explicit `comparison_genes`, flank policy, intake/kernel budgets, fresh output directory and request-bound cancellation path. |

Only genomic reference scope is admitted; transcriptome SAM is unsupported in S2.

**Counting default:** `identity_scope = sam_declared_query_unit`. A manifest may
group several secondary/supplementary rows as one producer-declared query unit;
different units may have the same QNAME. Never group by QNAME alone, count rows
as independent reads, or manufacture a FASTQ record index. Sample/barcode/primer
labels are optional producer metadata, not measured biological identity. Original
FASTQ identity verification and independent-molecule counting stay out of S2.
Public tables/counts explicitly say query units, not unique original reads.

All supplied rows, including unmapped/secondary/supplementary/QC-fail/duplicate
rows, are accounted for. Default: no MAPQ, locus, primary-only or duplicate
filtering. Retain admitted mapped rows across contigs, including off-target
ones; the selected gene must not filter them. This is observation retention,
not a promise that every mapping contributes to family support (see S2.3).
If a valid but unsupported row occurs, remove **all** segments of its unit from
kernel input and register the private unit key as unassessed, using S1's
unknown-target-length exclusion principle. Remaining good rows of that unit
must not generate possible/definite chains or increment family support.
Retain row-level diagnostics outside biological observations. Malformed data/
binding errors fail the request. A unit listed without alignment rows is
`not_provided`, not an observed unmapped read.

Accept unsorted input only within global record/operation/memory caps; group by
the manifest's exact row mapping, not an unbounded sort/spool or QNAME map.
Require consistent axis length and query orientation normalization across a
unit's rows. First version admits unpaired query units only: paired/mate flags
make the affected unit explicitly unsupported, not a merged molecule.

Keep file/accounting completion, selection scope and producer-declared search
coverage as separate facts. Default search completeness is unknown/incomplete;
reaching SAM EOF, MAPQ, an `NH` tag or absence of secondary rows cannot prove
exhaustive alignment. Missing `SA` companions or contradictory completeness
declarations must not strengthen support. Do not construct segments from `SA`.

Only an explicit, bound producer search-completeness declaration may request
`search_complete = true`, and only for units with no missing/unsupported views
or contradicted producer limitations. It remains **declared**, not independently
repeated alignment. By default, missing that declaration means reconciliation
`Unassessed` and zero definite-chain contributions; raw observations, possible
best-path chains and source/model geometry relations may still be available.
Separately, masks `None` make junction-flank qualification `Unassessed`, not
zero support. Even with a complete-search declaration, unknown masks do not
become primer-free qualified evidence. Report these defaults and axes explicitly.

### S2.3: Shared Interpretation, Saved Evidence And Adapters

Normalize admitted SAM geometry into the private S2.0 seam; use the existing
legal-path/reconciliation/model comparer. Do not infer a splice across separate
alignment rows, even if supplementary segments touch. Strand switching and
cross-contig hypotheses stay alternatives, not a fusion or single RNA call.

Publish a separate SAM report/data schema with request/raw/semantic/policy
digests, query identity scope, annotation receipts, coordinate authority,
all intake/coverage statuses and explicit nonclaims. Reuse portable table
projections where meanings agree; name the SAM table `query_structures.tsv`.
Catalog support counts are declared query units; junction rows are observations,
not independently summable support. No sequence/quality strings in output.

Preserve S0's hypothesis objective: maximize masked aligned supported-base
coverage, keeping all best-equivalent legal paths. `M` mismatch positions count
as supported under this inherited policy; there is no SAM NM/AS reranking.
Mark each observation's membership in retained optimal paths using their IDs.
Nonoptimal observations remain inspectable but do not contribute to `possible`
families. A longer pseudogene alignment can therefore displace a shorter spliced
alignment in support selection; this is a disclosed policy limitation, not
biological disproof of the observed shorter alignment.

In the **new SAM** catalog projection, separate mapping certainty from flank
qualification. Deduplicate per chain-or-junction and unit: qualification groups
are qualified-only, insufficient-only, unassessed-only or mixed across the
relevant retained occurrences, with mutually exclusive unit sets. Multiple rows,
hypotheses or skips must not multiply one unit's support. Preserve raw per-
observation statuses and exact count definitions; retain mapping-definite versus
possible support independently. This is not a change to old `ChainFamily` bytes.
Unqualified unannotated skips remain descriptive geometry, never candidate novel
biological junctions; qualified evidence still does not prove novelty/mechanism.

Carry `coordinate_authority = declared | verified` and the actual verification
receipt/scope. A hashed dictionary alone yields `declared`; a verified value
requires evidence of actual reference-byte verification, not just a supplied
digest. No implicit FASTA verification step in this slice. Equal physical keys,
provider models or joins across artifacts never upgrade source authority.

Read-only inspection consumes saved results and verifies requested bindings,
never reopens SAM to reinterpret it. Saved detailed analysis must remain usable
when its original SAM is offline, with artifact verification distinguished from
source revalidation. Optional bounded source-verification is explicit, not idle
work. Do not rewrite/import SAM observations as Phase 2 evidence or project DNA.

Use shared engine operations/parser/executor for CLI, generic MCP and agent
discovery, with thin scripting access. Tentative producer commands are
`rna-reads transcript-annotation-prepare` and
`rna-reads transcript-structures-from-sam`, both with strict JSON requests.
Specify typed saved-result inspection and request-bound cancellation companions
in the protocol before coding; reuse S1 helpers, not its strict v1 request type.
These names are **proposed, not currently available commands**. No GUI in S2.

## SAM Geometry Admission Contract

The parser must implement the documented SAM subset, not claim full SAM support.
Use the official [SAM v1.6 specification](https://samtools.github.io/hts-specs/SAMv1.pdf)
(printing `b5341fb`, 2025-08-12) and
[optional-tag specification](https://samtools.github.io/hts-specs/SAMtags.pdf)
(`a71164e`, 2024-09-09), inspected 2026-10-27. Proposed normalization:

- Check mandatory fields, header/order/dictionary consistency, flag ranges,
  integer overflow, positive operation lengths and legal terminal clipping.
  Support LF/CRLF parsing while hashing the **unaltered raw bytes**.
- Let `C` be query-consuming aligned operations, `a/b` terminal hard+soft clip
  lengths in stored orientation, and `L=a+C+b`. Forward query interval is
  `[a,a+C)`; reverse interval is `[L-a-C,L-a)`. This is a SAM-declared unclipped
  axis, not a verified original read length. Example `2H3S4M2I5M7S1H` gives
  `L=24`, forward `[5,16)`, reverse `[8,19)`. Test with an independent basewise
  oracle; do not copy the same formula into both implementations.
- Normalize terminal clips before `geometry::walk`; retain soft/hard clip and
  unknown-end evidence. `SEQ`/`QUAL` lengths must agree with declared operations
  when present; missing sequence does not imply missing geometry or zero length.
  Neither unclipped SAM ends nor complete splice chains establish RNA ends.
- Convert reference POS to zero-based; retain `N` as skip and `D` as deletion.
  Preserve `M` versus `=`/`X`; no invented identity or mismatch count. Do not
  interpret/retain NM/MD/AS as new evidence fields or ranking rules in this slice.
  MAPQ 255 is unavailable; zero remains reported zero.
- Admit `P` only as zero-consumption padding under the explicit unpadded policy,
  with a normalization diagnostic. Padded-reference SAM is unsupported. Refuse
  `CG` long-CIGAR/placeholder records as unsupported units rather than reading a
  placeholder as an intron. Normal textual long CIGAR is allowed within budgets.
- Validate bounded optional-field syntax and types/duplicates for consumed tags,
  including `CG` detection; define that small consumed set before implementation.
  No inventory of ignored tags. Preserve raw-source identity, and never use tags
  to synthesize absent alignments or RNA strand. FLAG reverse orientation is
  query/reference orientation, not transcription direction.
- A mapped row with absent CIGAR, no `M`/`=`/`X` support, or an unsupported
  placeholder is whole-unit unsupported before `walk`, not a lower-level span
  error that accidentally rejects unrelated units. A syntactically invalid
  CIGAR is still a data error. The same distinction applies to a valid declared
  query axis beyond S0's 2,000,000-base admission ceiling: retain a typed
  `query_axis_unsupported` unit, without bypassing global resource limits.
- An unmapped row with CIGAR/SEQ `*` remains accounted as unmapped with axis
  unverified. Its manifest length may be absent, never manufactured as zero or
  admitted to the kernel as a supported query. Contradictory producer length
  declarations are binding errors, not grounds to guess an axis.

Primer exclusions require a separately bound interval witness on the same
unclipped query axis, including its origin/transform. Project masks through
supported blocks with the shared helper, never through the aligned-span union.
Missing masks remain `None`/unassessed; a primer-ID list, absence of soft clips
or missing handle sequence cannot authorize an empty-mask assertion. An explicit
empty mask is a producer declaration with provenance, not an independent primer
search. Hard-clipped bases are unavailable, not observed primer-free bases.

## Bounds, Effects And Publication

Propose initial ceilings of 128 MiB total raw intake, 1 MiB per line,
65,536 SAM rows, 4,096 query units, 256 segments per unit, 100,000 dictionary
entries, 100,000 prepared transcript rows and 1,000,000 operations/exon entries.
Keep an explicit 128-MiB retained-data cap, bounded diagnostic examples plus
exact totals, the existing S0 work/path/model ceilings and a 128-MiB output cap.
Requests may lower limits. These are engineering admission defaults for review,
not measured throughput promises or whole-transcriptome support. Count all
files/tags/headers/crosswalks and allocation/work before accepting them; no
unbounded `read_to_end`, Cartesian materialization or silent truncation.
Hash and parse the same input stream/file snapshot, not one path opened once
for verification and again for interpretation. Diagnostics name row/field/reason
without copying raw SAM lines, read sequence or qualities into the report.

Each prepare request and each selected SAM cohort is one atomic analysis unit.
Cancel/control failure before completion leaves analysis/counts absent, not
zero. A completed unit retained after a late control failure keeps that outcome.
Hash/parse/resource errors publish no final report. Reuse fresh-directory,
bounded staged writers, close handles and publish the report last. Never
overwrite inputs, reports or project state. Destination safety must use the
applicable existing writer contract. Local main's separately implemented
Plan-only writer is not a general exporter contract or automatic protection for
new S2 operations; inspect its scope before considering any shared helper reuse.

Every report/marker write attempt requires confirmation. Name all new operations
in request-aware construction/no-save/approval policy and managed response/
result-binding handling. Retain generic MCP `confirm:true`. Raw execution keeps
ordinary auditing; standalone artifact-only publication does not merge a private
worker journal. Workflows remain conservative. Test project bytes, undo/redo,
structural revision and journals at each boundary; "read-only biology" is not
proof of no lifecycle effects.

## Required Verification Before Claiming Implementation

| Area | Deterministic acceptance path |
| --- | --- |
| Legacy freeze | Test-only S2.0a: ten existing Phase 2 cases plus exact S1 report/data/TSV, full Summary/Observation/RecordSupport/ChainFamily bytes, context/observation/chain digests, comparison work, covered errors, sort order and PAF stored-support assertion. No regenerated baselines to hide seam differences. |
| Gene variability | Two synthetic genes at distinct loci/opposite model strands; explicit selectors in engine/parser/CLI/MCP. Same physical chain, changed gene context: same chain key, different context and relations. |
| Geometry | Independently hand-calculated SAM/PAF equivalents; asymmetric clips, reverse flag, insertions/deletions/skips, `=`/`X`, `P`, hard clips and absent `SEQ`; exact bounds and path alternatives. Unknown identity never becomes zero or perfect. |
| SAM refusal | Malformed fields/clips/overflow, bad dictionary/SEQ/QUAL, unknown target, duplicate/bad consumed tags, padded mode, paired flags, CIGAR `*`, no supported bases, overlong axis and `CG`/untagged placeholders. Unsupported units send zero segments plus an unassessed key; no possible-family contribution. |
| Accounting | Repeated QNAME across distinct units; many rows per unit; arbitrary row order; missing/duplicate row assignments; zero-row units; unmapped, secondary/supplementary/off-target/QC-fail/duplicate rows; missing companions and incomplete search. |
| Primers/claims | Independently vary masks and declared search completeness; defaults give unassessed reconciliation/zero definite support and unassessed flanks. Assert absent masks alone do not change mapping certainty. Deduplicate qualification-split unit catalogs; no unqualified novel-junction label. |
| Hypothesis selection | A longer unspliced versus shorter spliced mapping: both observations remain, nonoptimal membership is explicit and only retained best-equivalent paths contribute support. No new score/filter or old-byte change. |
| Annotation | GTF `[101,110]` becomes `[100,110)`; exact ID forms, separate versions and selector mismatch; genePred crosswalk/columns; touching/overlapping/small-positive-gap policies and duplicate accessions; declared contig scope cannot drop an alignment competitor. Typed source refusal, no qualified-ID invention or silent model filtering. |
| Authority | Always snapshot-tagged prepared release plus nullable declared provider release. Same chain key with different verification evidence stays same key/different authority; joins do not upgrade declared to verified. |
| Resources/control | Exact-boundary/one-over cases for every cap, full EOF hash, cancellation/control/data errors, output failure, rollback/no final report, idempotent refusal of existing output and bounded saved inspection. |
| Portability | LF and CRLF raw-hash bindings through real file-backed boundary, paths with spaces through shared shell quoting, Windows handle-close/rename behavior. New committed byte-bound fixtures need provenance plus scoped `.gitattributes` and fast checkout guards. |
| Adapters/lifecycle | Parser/descriptor parity, forged planner flags still needing confirmation, raw/managed results, constructor/workflow/autosave and journal/redo preservation, actual CLI/MCP processes separately from unit tests. |

After separately authorized implementation, propose focused locked/offline Rust
tests for new paths, S0/S1/payload/lifecycle/parity regressions, headless plus
default-feature Cargo check, fmt/whitespace and fast Python/checkout guards.
Use the repository's existing test symbols and exact executable paths; finalize
new filters after tests exist. Do not present zero-matching filters or old S1
passes as S2 verification. No builds/tests of unimplemented code in this task.

Annotation format references: [Ensembl GTF fields](https://static.ensembl.org/info/website/upload/gff.html)
and [UCSC genePred definitions](https://genome.ucsc.edu/FAQ/FAQformat.html#format9),
inspected 2026-10-27. The live Ensembl `www` endpoint timed out; its official
static documentation was accessible. Source documents specify syntax, not
evidence that the historical archive has a complete reference or cohort.

## Glen Inputs And Deferred Scope

Through the GitHub-pushing thread's **single consolidated request**, after an
authorized committed candidate, ask Glen for a hash-bound SAM companion and
conversion/selection receipt for the existing BAMs, explicit query-unit mapping,
contig dictionary/reference authority, annotation companions/crosswalks and
available primer-mask evidence. Preserve originals; no script execution or
private-data publication is authorized by this plan. Missing inputs remain
named gaps. The archive's 21 rows and aggregate classes are cross-checks, not
expected answers to fit. Report ten exact versus one near canonical chain
separately; no boundary snapping to recover an old total.

Return one report with exact source/binary/lock/toolchain/profile/input hashes,
individual commands/counts, unavailable inputs and typed software/platform/
scientific verdicts. Do not require an interim main freeze or another full
audit per branch. No Glen request is sent by this planning task.

BAM/CRAM parsers, implicit `samtools`, aligners, remote data, full FASTQ/UMI
identity resolution, automatic annotation releases/aliases, source merging,
GUI/SVG/PDF/S3, protein consequences and real-data publication stay separate.
The supplied review challenged the private kernel seam, snapshot-version rule,
query-unit accounting and SAM normalization; it is reconciled above, not a new
sign-off on this revised text.

## Recommended Defaults And Remaining Owner Decisions

- Start with S2.0a tests only, followed by separately authorized S2.0b. No real
  SAM or annotation policy is needed to write synthetic preservation guards.
- Recommend GTF first, genePred later. Keep the exact row-to-unit manifest in
  this slice. A compact QNAME-uniqueness declaration could be a later reviewed
  producer contract, but must not silently replace current accounting or prove
  original FASTQ identity. Neither recommendation is an owner decision yet.
- Do not add contig aliases or qualified transcript IDs now. A future explicit
  alias table would need hash/authority/collision policy and its own review;
  it is not equivalent to guessing names, but is additional scope.
- Ask the owner to select the real-data minimum positive reference gap and
  accept/refine source refusal before real annotation admission. Choose no
  universal number here; keep the policy in the request/receipt digest.
- Glen's SAM/conversion, search-completeness, mask and annotation companion
  questions stay in the single consolidated handoff after a committed candidate.
  Missing those inputs does not block test-only or pure synthetic work, and
  receiving the review authorizes none of that implementation automatically.

## Verification Of This Planning Change

The initial planning session's documentation status and 11
`DocumentationStatusTests` pass, including their
LF/CRLF checkout coverage. Four local Markdown link targets in the two new
documents resolve; whitespace checks pass. Session-close reports 4 OK, 2 warnings
and no failures: intentional uncommitted documentation plus untouched `output/`,
and the manual plan-fidelity reminder. Scope review confirms five documentation
files only, no Rust/Cargo/fixture/generated-catalog changes, no implementation,
build, model invocation, commit, fetch, push or acceptance dispatch in that
initial session. The subsequent owner-authorized commit/rebase is a separate
integration action, not implementation or a new Rust verification result.

The 2026-10-28 review-reconciliation diff is documentation only. Source checks
confirm the input/output/digest coupling, independent completeness/qualification
axes, optimal-path support and annotation-boundary behavior described above.
They are not executed evidence for S2.0a, the planned normalizer or real data.
Documentation status, 18 LF/CRLF-aware documentation tests, five local plan/
prompt link targets and whitespace checks pass. Session-close has 4 OK,
2 warnings (intentional documentation diff/untouched `output/`, manual plan
fidelity), no failures. No Rust build/test, implementation, commit or push.
That paragraph records the earlier review-reconciliation session only; the
separately authorized S2.0a implementation/execution ledger above supersedes its
pre-implementation status, not its source-inspection limits.
