# Primer Payload Evidence (Experimental)

This separate Phase 2 consumer does not change Phase 1/RNA nomination. Runtime
acceptance is pending. See the [plan](primer_payload_partner_phase2_plan_20261026.md)
and [synthetic rehearsal](tutorial/10-04_targeted_nanopore_fastq_panel.md#phase-2-primer-free-payload-evidence-synthetic-only).

## Shared Contract

| Operation | Shared shell | Effects |
|---|---|---|
| `PreparePrimerPayloadAlignmentQueries` | `rna-reads primer-payload-prepare REQUEST_JSON_OR_@FILE` | Fresh query directory |
| `SummarizePrimerPayloadAlignments` | `rna-reads primer-payload-summarize REQUEST_JSON_OR_@FILE` | Fresh evidence directory |
| `InspectPrimerPayload` | `rna-reads primer-payload-show REPORT [--rows N] [--binding JSON_OR_@FILE] [--verify-artifacts]` | Read-only, 0..1000 rows |
| `CancelPrimerPayload` | `rna-reads primer-payload-cancel CONTROL_JSON_OR_@FILE` | Request-bound marker |

CLI, shared shell, generic MCP `op`, scripting and the inner agent share these
engine operations. No aligner, download, FASTQ discovery or provider/model call.
Artifact writes/cancellation require agent approval; MCP `op` still requires
`confirm:true`, including inspection. Named artifact-only managed execution does
not commit the detached project journal. Observation/control lanes are separate;
late cancellation cannot veto/relabel published evidence. There is no new GUI.

Strict versioned types: `src/engine/protocol/primer_payload.rs`. Preparation uses
`gentle.primer_payload_prepare_request.v1` → `gentle.primer_payload_queries.v1`.
Interpretation uses `gentle.primer_payload_alignment_request.v1` →
`gentle.primer_payload_alignment_report.v1`. Annotation is
`gentle.primer_payload_annotation.v1`; publication is
`gentle.primer_payload_result_binding.v1` (`kind:queries|alignments`). Unknown
fields are refused. Cancel wraps the **exact** submitted request as
`{"kind":"queries|alignments","request":{...}}`; the inner request is hashed.

Raw bindings contain `path`, prefixed `sha256`, and `bytes`. Semantic digests
clear `content_sha256`, serialize recursively sorted keys without whitespace,
then SHA-256. All budgets/policies participate in request identity. Capture actual
bindings, not guessed paths. Content consistency is not a signature, trusted
producer, alignment truth or exhaustive search. Interpretation does not repeat
FASTQ replay or alignment search. Retained-primary bases are additionally
compared to verified Phase 1 FASTA; all-anchor replay remains a content-bound
preparation witness, not authenticated provenance.

Preparation verifies the Phase 1 report/binding, selected sidecars, row counts
and anchor identities. `retained_primary_fasta` never reopens FASTQ and permits
only `primary` anchors. `replay_verified_fastq` explicitly rereads exact completed
original units and EOF hashes; only this mode permits `all` anchors. Missing
sequence is not inferred from coordinates. IDs bind source/file/record, read ID,
anchor geometry/direction and competing primer IDs; names need not be unique.
Queries run in observed primer direction, query zero adjacent to the anchor.
Other observed primer intervals are masked with N in **new copies** and excluded
from support, or explicitly unassessed. The denominator is unmasked query length,
not aligned length. Original reads and Phase 1 artifacts remain unchanged.

References explicitly bind ID, transcriptome/genome scope, assembly/release,
FASTA and normalized annotation. Transcript records associate exact target
name/length with sorted unique gene IDs. Genome scope requires exact contig,
strand and ascending nonoverlapping zero-based half-open exons. Optional
transcriptome genomic projection must have exon lengths summing to transcript
length. No aliases or version stripping. Expected-target crosswalks bind the
actual Phase 1 primer's `(locus_id,feature_id)` to reference gene/transcript IDs.

## Claim Prerequisites

| Observation | Required evidence | If missing |
|---|---|---|
| Query population | Source anchor accounting, unique IDs/lengths/hashes | Hard refusal |
| Producer-unmapped | Explicit accounting or consistent PAF no-hit row | Never inferred from absent rows |
| Aligned query support | Validated cg/cs consumption, target identity, masks | Strict refusal; coarse mode has unavailable exact fields |
| Expected fraction | Explicit source/reference crosswalk, block union, unmasked denominator | Null/unresolved |
| Primer edge | Best-equivalent detailed alignment reaches zero without boundary indel/mask | No extrapolated edge |
| Annotated splice compatibility | Genomic skip and exact annotated donor/acceptor/strand geometry | No inferred splice |
| Terminal-core comparison | Requested length, sufficient core, verified oriented A/C/G/T flank | Null mismatches/reason |
| Cross-locus candidate | Two ordered long detailed payload segments with disjoint nonempty gene sets | False; primer plus foreign payload is insufficient |

PAF requires contiguous query groups; a repeated group is a hard error. One
batch per reference per request is supported. Every query has `mapped`,
`producer_reported_unmapped` or `unassessed` accounting. Exit must be zero;
search completion/limitations are producer declarations, not uniqueness proof.
cg/short or long cs must consume both spans consistently; contradictions, tag
types, edit counts and length conflicts fail closed. MAPQ 0 is low/ambiguous,
255 unavailable. tp is an aligner role, not a supplementary flag; SA text creates
no extra segment. `alignment_detail_verified` means **geometry consistency**,
not independently reproduced base matching. Identity is explicitly reported.

Segments form a maximum unmasked aligned-support DAG of nonoverlapping query
spans. Retain all best-equivalent nodes/transitions without exponential chain
enumeration. Conflicting overlapping mappings, competing primers and low/missing
MAPQ remain ambiguous. Expected support is a union, never an HSP/isoform sum.
Overlapping HSPs are alternatives, not automatic biological concatenation.
Separate references are not ranked on incomparable scores. Thresholds are
explicit uncalibrated policy, not probabilities.

## Counts, Budgets And Publication

Partner keys include reference, sample/barcode/priming, full primer-ID group and
partner-gene group; cohort is report-bound. Distinct `(file_id,record_index)`
counts avoid multiplying anchors/HSPs. Anchor observations are separate.
Ambiguous groups are not exclusive members. Eligible/assessed/ambiguous/unmapped/
unassessed/excluded counts are read sets, not mutually exclusive across a read's
different anchors. Missing source/query units remain explicit deficits. Do not
pool selected/full cohorts or transcriptome/genome counts.

Primer/internal-left/internal-right clusters have separate edge-type keys.
Ambiguous, unresolved and unassessed mappings remain in read evidence but not
exclusive edge recurrence. Unresolved/unassessed rows do not increase profile
support; their anchor observations and coverage counters remain explicit.
Sort actual coordinates; start a new cluster when maximum span would exceed
tolerance. Representative = minimum observed edge; no transitive bridging.
Retain distinct-read support and denominators independently from partner ranks.

Fresh `report.json` plus hashed JSON-cell TSVs retain read/anchor evidence,
profiles and clusters. One named JSON column preserves nested sets without
delimiter guessing. Inspection returns bounded metadata, not the entire request.
Default prefixes are unverified; explicit `--verify-artifacts` hashes declared
outputs within budget. Neither mode reopens original FASTQ.

Input/output bytes, records, queries, query/line length, alignments, CIGAR ops,
annotations, contexts and chain/annotation comparisons have hard budgets.
An additional 128 MB ceiling bounds JSON, retained sidecars and memory evidence.
FASTA streams once per reference, retaining only requested seed flanks. No top-N
truncation. Data/hash/resource errors publish no final report. Cancel/control
faults retain only EOF/hash-verified source units/reference batches. Stage, close
and reserve fresh paths; orphan files without a final report are not evidence.

Rust/CLI/MCP, Windows paths/handles, live-agent responsiveness, real accuracy/
throughput and scientific acceptance remain external gates. No molecule,
expression, exclusive isoform, fusion/trans-splicing or mechanism approval.
