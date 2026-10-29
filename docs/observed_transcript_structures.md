# Observed Transcript Structures (Experimental S1)

Source-relative near offsets, motif context and explicit sample recurrence use
the separate [J1 assessment](observed_splice_junction_assessment.md). It consumes
saved S1 evidence without changing this v1 contract.

Interpret **saved genomic Phase 2 evidence**, not raw reads. Ask which exact
splice patterns were observed and which supplied gene/transcript models remain
compatible. A compatible RNA pattern does not establish a complete transcript,
an exclusive isoform, independent molecules, expression or a protein.

Original implementation baseline: local `main`/`gentle_rs_2_main` at `a2637079`
plus the S1 working diff. Committed S1 was rebased as `929c54e4` onto `b8e0118f`;
the [post-rebase receipt](observed_transcript_structures_plan.md#s1-commit-and-rebase)
records the repeated checks separately. Focused macOS verification is not
a native, cross-platform, scientific or release verdict.
The [reviewed plan and ledger](observed_transcript_structures_plan.md) retain
the earlier S0 results separately. No new aligner, dependency, GUI or SAM/BAM
ingress is added, and existing Phase 2 report/TSV output is not rewritten.

## Shared Commands

```text
rna-reads transcript-structures REQUEST_JSON_OR_@FILE
rna-reads transcript-structures-show REPORT [--table TABLE.tsv] [--rows N] [--binding JSON_OR_@FILE] [--verify-artifacts]
rna-reads transcript-structures-cancel REQUEST_JSON_OR_@FILE
```

These are shared-shell commands, also available through CLI shell execution.
The raw operations are `SummarizeObservedTranscriptStructures`,
`InspectObservedTranscriptStructures` and `CancelObservedTranscriptStructures`.
Generic MCP `op` requires `confirm:true`, even for inspection; generic scripting
and agent routes use the same typed engine operation, not another analyser.
Capability discovery provides nonexecutable placeholders, not evidence or a
default gene. Report/marker write attempts require approval independently of
project-state mutation. There is no dedicated graphical viewer in this slice.

## Request

`src/engine/protocol/observed_transcripts.rs` defines strict versioned records.
Summarization takes `gentle.observed_transcript_structures_request.v1`:

| Field | Meaning |
|---|---|
| `source` | Exact `gentle.primer_payload_result_binding.v1`, kind `alignments`; preserve its actual outcome. |
| `source_report` | Raw `path`, prefixed `sha256`, `bytes` for that report. The path must agree with the binding. Raw and semantic hashes are distinct. |
| `batch_id` | Exactly one explicitly selected saved genomic batch. Transcriptome-only batches are refused, not projected speculatively. |
| `reference_gene` | Mandatory `{source_id,gene_id}` selecting a stable gene ID in a named source. No symbol, species, locus, strand, canonical transcript or TP53 default. |
| `comparison_genes` | Explicit stable gene IDs in additional sources; equal names are not automatic cross-source identity. |
| `annotation_sources` | `{source_id,provider,reference}` entries. `reference` reuses `PrimerPayloadReference` and the existing normalized annotation schema. |
| `required_flank_bases` | Explicit uncalibrated support requirement, included in the digest; not a biological universal. |
| `budgets` | Existing `PrimerPayloadBudgets` for intake/publication, all required and positive. |
| `kernel_budgets` | Existing S0 structural ceilings, including legal-path alternatives, operations, model comparisons and total work. |
| `output_dir`, `cancel_path` | A fresh output directory and optional exact-request cancellation marker. |

Each annotation reference must declare genomic scope and the selected batch's
exact assembly and reference FASTA digest. Its own normalized annotation ID,
release, raw hash, target lengths, gene IDs, strand and exon geometry are
validated independently. Providers may use different source-specific gene IDs.
Unknown/ambiguous gene resolution remains typed, without falling back to another
gene. Off-target observations are retained instead of filtered by the selector.

S1 validates the source report's semantic identity and actual raw bytes; source
trace and query manifests; retained query/trace dependencies and exclusions;
all published Phase 2 artifact hashes to EOF; and exact unique completed
`(batch_id,query_id)` row accounting. Read/source/sample/barcode/primer context
must agree with the bound queries. Rehashed contradictory geometry is still an
error. Old classification, expected-target fractions and broad transcript/gene
labels do not decide structure compatibility.

No original FASTQ, PAF or reference FASTA is reopened. Those source digests are
recorded producer evidence, not a new search or base-matching verification.
The saved prepared query FASTA and retained Phase 1 payload witnesses **are**
verified. A previously approved all-anchor replay is not run again here.

Phase 2 segments do not retain every target length. S1 uses bound annotation
lengths; if any competitor's length is unavailable, the **whole query view**
is unassessed, rather than promoting the remaining mappings. Excluded or
unassessed anchor views also prevent definite support from another anchor of
the same original record. Original source digest and record index, not QNAME
alone, govern counting. Missing input units remain disclosed source coverage
gaps; they are never extrapolated into records or biological zeros.

## Saved Results

The compact `gentle.observed_transcript_structures_report.v1` contains the exact
request, semantic digest, stage outcome, source/cohort/completed/missing-unit
metadata, optional counts, artifact hashes and nonclaims. MCP and oversized
managed results return `gentle.observed_transcript_structures_result_binding.v1`,
not an unbounded observation population. Small managed results retain the typed
metadata report; the job receipt always carries its result binding. Bindings
name the actual report path, request/content
digests and outcome. They establish internal consistency, not authentication.

The fresh directory contains `report.json` and these hashed `artifacts/` files:

| Artifact | Contents |
|---|---|
| `structures.json` | `gentle.observed_transcript_structures_data.v1`: full typed kernel analysis, if completed, plus unavailable views. |
| `read_structures.tsv` | Original-record support, legal hypotheses, query/primer context, endpoint constraints and observations. |
| `structure_catalog.tsv` | Exact physical chain families, record support, and nonadditive sample/barcode/priming/primer groups. |
| `reference_compatibility.tsv` | Every source-specific model relation for each observation; no winner or protein label. |
| `junction_catalog.tsv` | Exact skipped intervals by observation/record, flank qualification, best-path membership and source-specific annotated-model matches. Rows are observations, not independent counts. |
| `unavailable.tsv` | Excluded/missing views and reasons, including unknown target lengths and uncompleted batches. |

TSVs use one quoted `row_json` column, preserving nested geometry and alternatives
without inventing a second table analyser. Read sequences are not exported.
S0's Rust API remains experimental, but S1 v1's published serialized shape is
versioned: later changes need an explicit compatibility/version policy, not
regeneration of old evidence.
Different anchors may appear in multiple rows/groups; do not add their counts.
Junction annotation is `null` when comparison is unavailable, not an empty
annotation match. Novelty means relative to bound source snapshots only.

The physical chain digest covers exact genomic skips and coordinate authority,
not the selected gene, model strand or read endpoints. Query orientation is not
RNA transcription direction. Unspliced records retain observations but no
empty-chain family. Detailed `N` skips, supported `M/=/X` blocks, masks, `D`
deletions and coarse/unavailable geometry remain distinct. Tolerance is zero.

## Control, Limits And Lifecycle

The selected batch is one atomic analysis unit. Cancellation is cooperative
between bounded intake/query/segment/record units. It is not immediate worker
termination. A cancelled/control-error report has absent counts if analysis
did not complete; a completed unit may be retained with a later cancel/fault
outcome. An uncompleted source batch produces `unavailable` with absent counts,
not an empty measured distribution. `completed` describes the selected analysis,
not independent validation or completeness of the entire source experiment.

Hash, data, geometry and resource errors publish **no final report**. Output
must be a fresh directory; all artifacts are staged, writers are byte-bounded
and closed before publication, and the final report is written last. Orphan
artifacts after a filesystem publication failure are not completed evidence.
The report's semantic digest clears `content_sha256` and hashes sorted-key JSON.

S1 admits at most 32 annotation sources, 4096 prepared queries/evidence rows,
and 128 MB input/output budgets. Lower caller limits still apply globally,
including retained trace/query dependencies and nonselected-batch evidence.
Individual records are bounded by `max_line_bytes`; inspection returns at most
1000 rows and 8 MB of selected JSON. Over-budget requests fail, never truncate.

Inspection reads saved S1 output only; it does not reopen source annotations or
reports. Default zero-row inspection validates report identity, not sidecar
bytes. Requesting rows verifies the entire selected table at EOF; explicit
`--verify-artifacts` verifies all saved sidecars, not a new biological analysis.

These three named operations join the request-aware preserve-content lifecycle:
constructors do not repair the project and CLI/MCP suppress implicit saves.
Raw `apply` still appends the existing execution audit, without checkpoint,
undo/redo or project-content changes; it does **not** promise a journal-free
raw API. Managed standalone publication discards its detached journal and does
not merge into the live project. Shell inspection/control use independent lanes.
Mixed workflows keep the existing conservative lifecycle; no broad migration.

## Verification Scope

Synthetic Rust regressions cover actual file-backed Phase 2 intake, opposite
annotation strands and variable gene contexts, preserved source bytes, deletion
of original FASTQ/PAF/FASTA before consumption, LF/CRLF bytes, tampering,
accounting, resource errors, cancellation/faults, saved inspection, shared shell,
MCP and managed publication. The owner subsequently authorized focused local
S1 execution on 2026-10-27, following the existing reviewed plan without another
Claude consultation. Final macOS/aarch64 runs pass S1 **14/14**, S0 **24/24**,
primer-payload **46/46**, lifecycle Slice 3 **11/11**, full-registry introspection
**1/1**, and canonical parity producer/document equality **1/1**, with no failed
or ignored tests in those final groups. Filters overlap; these are not summed
as distinct tests. The Phase 2 group retains ten old serialized-row/work cases.

Locked/offline Cargo check passes for the headless library, CLI and MCP, and
with default features and `--tests` (including GUI result constructors).
The first build exposed omitted optional fields in three headless result
constructors and an invalid `OpResult::default()` in a new test. All ten explicit
headless/GUI constructors now set the fields to `None`; the serialization test
uses a real engine result. The first S1 execution passed 13/14: its sole failure
incorrectly expected the compact envelope for a small managed result. The test
now checks the typed small report **and** its exact receipt binding, without
changing publication behavior or weakening live journal/redo assertions.

The exact working-tree source/binary fingerprints and commands are in the
[executed ledger](observed_transcript_structures_plan.md#s1-executed-verification-2026-10-27).
Existing headless dead-code and large `__eh_frame` linker warnings remain;
none were suppressed. Actual CLI-process, scripting VM, native, Linux/Windows,
full-suite and scientific acceptance remain for the consolidated Glen/CI pass.

Local fast checks: 27 Python tests (26 passed, one expected binary-dependent
skip), documentation-status checks, formatting and whitespace pass. The Python
suite includes LF/CRLF checkout/documentation checks and six S1 source-contract
guards, not execution of the Rust analyser. Session-close reports four OK, two
expected warnings (intentional working diff/generated-looking glossary and
manual plan-fidelity review), and zero failures. Unrelated `output/` is untouched.

```bash
CARGO_BUILD_JOBS=1 cargo test --locked --offline --no-default-features --lib observed_transcripts -- --test-threads=1
CARGO_BUILD_JOBS=1 cargo test --locked --offline --no-default-features --lib observed_transcript_structures -- --test-threads=1
CARGO_BUILD_JOBS=1 cargo test --locked --offline --no-default-features --lib primer_payload -- --test-threads=1
CARGO_BUILD_JOBS=1 cargo test --locked --offline --no-default-features --lib lifecycle_slice3 -- --test-threads=1
CARGO_BUILD_JOBS=1 cargo test --locked --offline --no-default-features --lib execute_introspect_capabilities_projects_full_registry_with_fact_annotations -- --test-threads=1
CARGO_BUILD_JOBS=1 cargo test --locked --offline -p gentle-protocol shared_workbench_route_keeps_canonical_parity_projection_current -- --test-threads=1
CARGO_BUILD_JOBS=1 cargo check -q --locked --offline --no-default-features --lib --bin gentle_cli --bin gentle_mcp
CARGO_BUILD_JOBS=1 cargo check -q --locked --offline --tests
```

The parity test calls the same canonical producer as `parity-matrix-check`;
the standalone documentation-helper CLI was not built/run. Recheck affected
guards after integration or source changes, without repeating unchanged passes
as new evidence. Native/platform/scientific review remains separate;
SAM-text/BAM ingress, observed-lane GUI/rendering and protein consequences are
later slices, not implied by this delivery.
