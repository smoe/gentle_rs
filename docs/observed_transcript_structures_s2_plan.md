# Observed Transcript Structures S2: SAM And Annotation Preparation

Status: **proposal for read-only review; not implemented or authorized for
implementation**. Prepared 2026-10-27 at local `main` = `gentle_rs_2_main` =
`1804605f1ba42369f943213414d79fd9a40013cb`. No fetch or remote-tip claim.
Tracked files were clean; unrelated untracked `output/` is untouched.
The [roadmap](roadmap.md) owns priority/release scope. This is a bounded expansion
of [the umbrella plan's S2](observed_transcript_structures_plan.md#slice-boundaries-and-remaining-work),
not a new release gate. The [Claude prompt](observed_transcript_structures_s2_claude_review_prompt.md)
has not been sent or reviewed. Implementation needs a separate owner instruction.

Integration context, 2026-10-28: the owner requested committing this proposal
and rebasing onto local `main`, now `7329f00f607c82c195ea13d931f83b0f1bc93150`.
Its two additional commits introduce the Plan-only destination writer and
reorganize Help/manual documentation. They do not implement observed-transcript
S2 or change its kernel/Phase 2 paths. Preserve the original inspection baseline
above; recheck the actual checkout before review. No fetch or push is requested.

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
| `src/observed_transcript_structures/kernel.rs`: `observation`, `summarize_with_control` | The context digest binds serialized `Input`; observations contain `payload_query_length`. Preserve S1 digests and field meanings, not just comparable final counts. |
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

## Original Codex Plan: Four Reviewable Steps

### S2.0: Freeze Old Bytes And Introduce Only The Necessary Private Seam

First add preservation guards, then a private format-neutral geometry/query
view consumed by the **same** observation, path, reconciliation and comparison
logic. Project existing PAF inputs losslessly; provide a separate SAM projection.
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

Initial formats are explicitly selected plain-text GTF and UCSC genePred
tabular variants with an explicit column layout (including any `bin` column).
No sniffing between GFF3, BED12 or JSON. The archive's UCSC JSON needs an
owner-prepared tabular companion with extraction/hash provenance; this task
does not promise automatic loading of every archive member.

- For GTF, normalize exon rows from one-based inclusive to zero-based half-open
  geometry; require consistent source gene/transcript IDs, contig and strand.
  Retain identifier versions; if version attributes accompany an unversioned ID,
  record the exact normalization policy and reject conflicting version claims.
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
  Accessions repeated across primary/alternate loci must not be merged. Missing
  competing-source gene selection remains a typed resolution gap at comparison.
- A genuinely absent selected gene is a typed `missing_gene` preparation result,
  not a malformed row or an empty completed model set. Retain its source receipt
  as unavailable comparison evidence; do not silently drop that requested source.
- Preserve provider/source separation. A UCSC `ncbiRefSeq` track snapshot is not
  a direct NCBI GCF annotation release. Bind retrieval/source metadata and notices
  when supplied; no online retrieval, source replacement or inferred release.

**Proposed missing-release rule, for Claude to challenge:** record nullable
`provider_release` and `version_kind` in the preparation receipt. With no known
provider release, set the reused annotation's required `release` to the explicit
content identity `snapshot-sha256:<64-hex-raw-input-digest>`, never "Ensembl 116" or a
date-derived guess. Require this receipt at S2 ingress and label it as a snapshot,
not a provider release. Selection/crosswalk/policy hashes also bind the prepared
artifact. A missing release therefore limits claims without creating a permanent
activation blocker or changing S1 validation semantics.

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

**Counting default:** `identity_scope = sam_declared_query_unit`. A manifest may
group several secondary/supplementary rows as one producer-declared query unit;
different units may have the same QNAME. Never group by QNAME alone, count rows
as independent reads, or manufacture a FASTQ record index. Sample/barcode/primer
labels are optional producer metadata, not measured biological identity. Original
FASTQ identity verification and independent-molecule counting stay out of S2.
Public tables/counts explicitly say query units, not unique original reads.

All supplied rows, including unmapped/secondary/supplementary/QC-fail/duplicate
rows, are accounted for. Default: no MAPQ, locus, primary-only or duplicate
filtering. Keep mapped alternatives across contigs, including off-target ones;
the selected gene must not filter them. If a valid but unsupported row occurs,
mark its **whole query unit** unavailable so supported competitors cannot be
promoted by omission. Malformed data/binding errors fail the request. A unit
listed without alignment rows is `not_provided`, not an observed unmapped read.

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
  Preserve `M` versus `=`/`X`; no invented identity or mismatch count. Retain
  producer `NM`/`MD`/`AS` evidence separately if present, not as new ranking rules.
  MAPQ 255 is unavailable; zero remains reported zero.
- Admit `P` only as zero-consumption padding under the explicit unpadded policy,
  with a normalization diagnostic. Padded-reference SAM is unsupported. Refuse
  `CG` long-CIGAR/placeholder records as unsupported units rather than reading a
  placeholder as an intron. Normal textual long CIGAR is allowed within budgets.
- Validate relevant optional-tag types/duplicates and bounded tag lengths;
  preserve raw-source identity and a diagnostic inventory for ignored tags.
  Never use tags to synthesize absent alignments or RNA strand. FLAG reverse
  orientation is query/reference orientation, not transcription direction.

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
| Legacy freeze | Run the ten old Phase 2 serialized-row/work cases and all affected payload tests. Add exact S1 report/data/TSV/context-digest guards, preserving optional-field serialization. No regenerated baselines to hide differences. |
| Gene variability | Two synthetic genes at distinct loci/opposite model strands; explicit selectors in engine/parser/CLI/MCP. Same physical chain, changed gene context: same chain key, different context and relations. |
| Geometry | Independently hand-calculated SAM/PAF equivalents; asymmetric clips, reverse flag, insertions/deletions/skips, `=`/`X`, `P`, hard clips and absent `SEQ`; exact bounds and path alternatives. Unknown identity never becomes zero or perfect. |
| SAM refusal | Malformed fields/clips/overflow, bad dictionary/SEQ/QUAL, unknown target, duplicate/bad relevant tags, padded mode, paired flags and unsupported `CG`; whole-unit unavailability never promotes another mapping. |
| Accounting | Repeated QNAME across distinct units; many rows per unit; arbitrary row order; missing/duplicate row assignments; zero-row units; unmapped, secondary/supplementary/off-target/QC-fail/duplicate rows; missing companions and incomplete search. |
| Primers/claims | Mask witnesses on both orientations, clipped/missing witnesses, missing versus explicit empty masks; no artificial junction, natural-end, RNA-strand, independent-molecule or protein assertion. |
| Annotation | GTF boundary `[101,110]` becomes `[100,110)`; genePred exon arrays; version/crosswalk conflicts; same accession on different loci/strands; snapshot-only version provenance; source disagreement, missing genes and complete-transcript selection. |
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
Before implementation, Claude should challenge especially the private kernel
seam, snapshot-version rule, query-unit accounting and SAM normalization here;
the owner then chooses the first bounded implementation step.

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
