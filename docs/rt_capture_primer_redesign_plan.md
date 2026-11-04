# RT Capture Primer Assessment And Redesign

Status: original source-grounded plan with the implementation ledgers below,
2026-11-04; no scientific acceptance. Prepared after fetching `origin/main` at
`89df315a80e0d3ae02e5685efd7a94e3895fbaa7`; the working branch already matched
that commit. A read-only Claude consultation was offered; this plan has not
been Claude-reviewed. At the owner's request, a
[ready-to-pass review prompt](rt_capture_primer_claude_review_prompt.md) is
prepared separately; no Claude invocation or waiting is required.

## Objective And Evidence

Assess an adapter-bearing RT primer as an initiation site on RNA, separately
from its full-length genomic matches and subsequent enrichment. First assess
supplied candidates; compose automatic redesign only after that assessment
has a deterministic contract. Keep this within the official 0.1
[specificity/redesign aim](official_0_1_release_aims.md#3-flagship-global-rt-pcr-specificity-and-redesign),
without adding an internal.13 release gate.

Source: Glen's [2026-11-03 report](https://gentle.functional.domains/openclaw/canvas/tp73-rt-primer-redesign-20261103/REPORT.md)
and [candidate table](https://gentle.functional.domains/openclaw/canvas/tp73-rt-primer-redesign-20261103/candidates.tsv),
read over HTTPS on 2026-11-04. These are external reported results, not a
GENtle replay. The source identifies GRCh38, Ensembl 116 and Primer3 2.6.1;
the underlying reference bytes, private reads and search completeness have
not been independently verified here.

| Reported candidate | Tail, 5′ to 3′ | Model coverage | Old-junction rows with ≥6/8 matches |
| --- | --- | --- | --- |
| Old | `TGCTGTAGTCGTGGCCCTGCTT` | 17/20 | 70,396/147,185 |
| D_shared | `GAAGTCCGCCCACTCGTCAG` | 17/20 | 109/147,185 |
| L14_long_utr | `GCTCGGAGAAGCACGATCTTAC` | 13/20 | 539/147,185 |

All three reportedly have one exact full-tail genomic hit. D_shared is a
coding-terminal-exon pilot; L14 is a long-ending 3′-UTR-group pilot. Neither
replacement is experimentally validated. The old tail's terminal `CTGCTT`
is complementary to `AAGCAG` at many observed MT-RNR2 transitions. This
supports an RT-mispriming hypothesis, not proof of reaction history.

The old-junction cohort was selected by the old primer. Its fractions cannot
estimate new candidates' prospective off-target rates. Retain that cohort as
a separate retrospective comparison. Glen's alternative of regular cDNA
followed by target-specific enrichment is a distinct workflow choice; the
software must not silently substitute it for a requested RT-primer design.

## Existing Owners And Actual Gaps

| Requested behavior | Existing shared owner | Remaining work |
| --- | --- | --- |
| Separate handle and target tail | `RnaLibraryPrimer` / `RnaPrimerComponent` in `src/engine/protocol/primer_trace.rs`; terminal-exon design records in `src/engine/protocol.rs` | Explicit RT assay meaning and validated conversion between report candidates and component-defined trace inputs |
| Terminal 6/8/10/12 comparisons | `PrimerPayloadPolicy.terminal_seed_lengths`; `primer_payload/evidence.rs` computes positional mismatches at verified mapped edges | Prospective whole-background counts for each new candidate, with explicit scope, completeness and denominators |
| Oriented primer–payload inspection | `rna_panel/primer_trace` and `rna_panel/primer_payload` | Joined inspection of retained anchors, seed evidence and partner/edge profiles; preserve current ambiguity and read-count units |
| Per-isoform coverage and CDS/UTR | `transcript_capture.rs`, shared transcript-template geometry and coverage-universe resolution | Per-binding CDS/5′-UTR/3′-UTR/boundary/unknown classification, separately from site absence |
| Candidate comparison | Terminal-exon and transcript-capture reports, shared oligo metrics and reference-bound specificity receipts | One portable assessment tying those records to prospective and retrospective evidence without borrowing pair-assay approval |

Paths in the last three rows are under `src/engine/analysis/`, with primer
trace/payload under `rna_panel/`; capture records are in
`src/engine/protocol/transcript_capture.rs`. Existing source is not evidence
of current-source runtime, native GUI or biological acceptance. Reuse the
[biological extension seams](biological_extension_guide.md),
[capture boundaries](transcript_capture_pools.md) and
[payload evidence contract](primer_payload_evidence.md).

## Proposed First Slice: Assess Supplied RT Candidates

This is the original Codex plan. It proposes a single shared assessment
request/report, with thin shared-shell and generic operation access. Names,
numeric admission budgets and exact wire schemas must be finalized with the
implementation; no command in this document is advertised as executable.

### Primer Identity And Binding Contract

Use a typed RT-specific view of the existing oligo/component representation,
not another independently editable sequence store. Require an explicit
5′-handle/3′-target-tail partition, complete ordered oligo, candidate identity,
source/revision citation, assay stage and intended transcript group. Validate
that the complete sequence equals handle plus tail. A core-only primer must
be an explicit permitted form; a shared handle cannot identify a gene.

Preserve original report identities and hashes when converting a design;
refuse disagreement rather than reconstructing an allegedly ordered oligo
from a gene label. Do not change old v1 digest semantics, component roles or
serialized output to add RT meaning. Prefer an assessment-owned wrapper over
a mandatory field on every legacy primer. Candidate identity also binds the
handle even when two candidates share the same target tail.

All oligos are supplied 5′ to 3′. RNA references are interpreted in biological
5′-to-3′ orientation; genomic sequences require explicit strand handling.
Coordinates are 0-based, half-open internally, with labelled conversion for
external reports. RT binding compares the tail against complementary RNA;
read orientation, genomic strand and assay role remain independent.

### Prospective Terminal-Risk Matrix

Assess terminal 6, 8, 10 and 12 bases of the target tail, never of the handle.
Report exact counts separately from a positional mismatch histogram for each
length. “≥6/8” means at least six equal positions in an oriented eight-base
comparison; it does not mean six contiguous matches. Also retain the exact
contiguous run from the primer's 3′ end, so mismatches at the terminal base
cannot masquerade as a compatible extension end. No indel search is implied.

Every cell binds candidate, reference content, target space, searched strand
policy, seed length and mismatch rule. Genome, transcriptome, mitochondrial,
ribosomal and user-supplied frequent-RNA scopes remain distinguishable.
Overlapping scopes and transcript isoforms can contain the same sequence;
do not sum their hits as independent molecules or independent risks.

Report both matching positions and distinct reference records with a hit;
state whether palindromic genomic occurrences count once physically or once
per orientation. Prefer one physical site with explicit strand alternatives.
Keep all hits and intended-target-excluded counts separate, using explicit
reference-bound intended-site mappings rather than gene-name exclusions.
Bind circular reference topology where supplied; linear scans cannot claim
to cover origin-spanning mitochondrial sites.

Unknown sequence windows, short references, missing inputs and unsupported
contexts have explicit assessed/unassessed counts and reasons. Ambiguous
bases cannot become exact matches or clean negative evidence. An omitted
reference is unassessed, and an incomplete or budget-limited scan cannot be
labelled zero-risk. Fail on exceeded hard limits rather than silently retain
only the first hits. Count exhaustively within the declared scope while
retaining only a separately bounded example-site list.

Scan each verified reference once for the supplied candidate set, reusing
existing reference/index access where semantics agree. Bind maximum candidates,
reference bytes/bases, evaluated windows/comparisons and retained sites in
the effective request. Choose numerical caps from a bounded synthetic work
estimate before coding. Hash/parse the same source bytes; include cancellation
and source-change handling at the established verified-unit boundary.

For the first slice, report unweighted burden. Abundance weighting is
`not_assessed` unless a later explicit, content-bound weight table supplies
reference identities, sample/context, units and missing-weight policy. Never
invent weights from RNA labels or equate PCR read counts with RNA abundance.

### Retrospective Evidence And Report

Consume an explicitly supplied saved junction-evidence binding only when its
reference flank and source identities can be verified. Reuse oriented trace
geometry and Phase 2 seed comparisons; do not replay private FASTQ or build a
new aligner. Report selected-cohort identity, unique source-read counts,
anchor-observation counts and unavailable-flank denominators separately.
If saved evidence lacks the necessary flanks, mark comparison unavailable.

The portable assessment contains normalized inputs and hashes, all reference
statuses, per-candidate matrix cells, diagnostics, completeness and explicit
limitations. Inspection/export read this record. New references or revised
candidates require a new assessment; inspection cannot refresh them silently.
Assessment is advisory and does not change existing order-readiness receipts.

Reuse the existing artifact publication or persisted-report policy after
choosing the actual operation effect. Do not infer write authorization from
a read-only label. All supported adapters reach the same engine operation;
GUI work is deferred until the report is usable through shared-shell/MCP.

## Subsequent Slices

1. **Isoform context and comparison:** classify every exact candidate binding
   against each transcript's annotated CDS and UTR intervals. A boundary-crossing
   site is mixed; absent binding and unknown annotation are different states.
   Preserve the requested transcript universe, including noncoding and missing
   models. Display tail/full-oligo Tm/GC and interaction diagnostics with their
   model and conditions. Reuse Primer3 where appropriate; heuristic runs do
   not become thermodynamic hairpin energies. Concentrations are recorded
   metadata until a calibrated competition model exists.
2. **Junction inspector projection:** join existing full-primer orientation,
   observed trimmed interval, masked payload, partner-equivalence profiles and
   edge evidence. Approximate primers trim at observed alignment boundaries;
   handle-only hits remain unattributed. Display recurrent partners with
   distinct-read denominators, never inferred molecule counts. Keep genuine
   biological alternatives unresolved and retain primer-free payload evidence.
3. **Constrained redesign:** compose existing candidate enumeration with the
   assessment. Explicit shared-transcript and long-ending-3′-UTR objectives
   have different coverage contracts. Never fall back from unknown UTR to CDS
   or drop excluded models from a denominator. Test candidate-specific new
   backgrounds, then report a transparent multi-criterion comparison with
   stable tie breaking; avoid an uncalibrated universal specificity score.
4. **Abundance and native presentation:** add explicit weight-table semantics
   and native transcript lanes/comparison controls after portable outputs and
   missing-data behavior are verified. Keep GUI handlers thin and measure
   large-reference work outside the paint thread.

## Deterministic Acceptance And Handoff

First-slice tests use hand-crafted sequences with documented provenance:

- handle-only matches versus terminal-tail hits; malformed partitions and
  legacy digest preservation;
- complementary RNA and both genomic orientations, reverse-strand targets,
  terminal-base mismatch versus internal mismatch, and ≥6/8 positional matches
  versus a contiguous terminal run;
- repeated sites versus distinct records, duplicate/overlapping reference
  scopes, palindromic sites, circular-origin windows, short tails/references
  and ambiguous bases;
- prospective results differing from the old-selected cohort, with no transfer
  of retrospective percentages into predicted specificity;
- explicit intended-site exclusions, changed reference/candidate hashes,
  unavailable flanks, zero assessed windows and absent abundance weights;
- deterministic counts and order, bounded retained examples, hard work-limit
  refusal, cancellation/source failures and unchanged legacy trace/payload bytes;
- typed operation/shared-shell/MCP reachability, write approval, serialization,
  inspection/export and no automatic readiness upgrade.

For later slices add mixed CDS/UTR and unknown-CDS oracles on both strands,
multiple binding sites, absent sites, retained isoform denominators, full-primer
indels and ambiguous partner grouping. Native inspection, real-reference
performance and wet-lab interpretation need separate evidence.

Run targeted Rust tests and `cargo check -q` for implementation, plus affected
registry/parity/tutorial and documentation checks. Follow the repository's
LF/CRLF and fixture-provenance rules if introducing file ingress or fixtures.
Do not add hypothetical commands to tutorials. Record exact tested source,
platform and unrun Windows/native checks; no earlier pass accepts later code.

No underlying reads or reference cache were fetched, no assay was ordered,
and no private analysis was rerun.

## First-Slice Implementation Ledger

The owner subsequently authorized supplied-candidate assessment against
explicit background references, without waiting for Claude. This selects the
prospective part of the original first slice; saved retrospective comparison
stays deferred with `retrospective_status: not_assessed`. No automatic design,
isoform/CDS/UTR classification, abundance weighting, thermodynamic comparison,
native form, private replay or order-readiness attachment is added.

Implemented `AssessRtPrimerCandidates` / `primers assess-rt-candidates`, strict
versioned request and portable result in `OpResult.rt_primer_assessment`.
Use the [executable contract](cli/rt-primer-assessment.md), rather than the
proposed text above, for exact schemas, counting units and admission limits.
The operation reuses `RnaLibraryPrimer`, existing panel validation and
`PrimerPayloadFile` with the same verified byte reader and bounded line parser.
Independent per-primer validation allows different handles to share a tail
without changing the legacy trace panel's duplicate-core rule or digest.

Genome counts use both template orientations; RNA scopes use the declared
biological RNA orientation. Physical-window mismatch and terminal-run
histograms are separate marginal summaries. Verified full-tail intended sites
remove only their own template orientation. Circular-origin seed windows are
supported without multi-lap matching; wrapping full-tail intended exclusions
remain unsupported and fail explicitly. Missing scopes, ambiguous windows,
short references and short tails retain their own unassessed accounting.

Limits are bounded subset support, not whole-human-genome acceptance: default
one million bases and 100 million base comparisons; one candidate on a
million-base both-strand reference requires at most 72 million comparisons.
Maximum admitted comparison budget is one billion. Counts stay exhaustive
when bounded exact-example lists fill; excess work or source/control failure
returns no partial assessment. The operation has no marker path or external
write. Existing callback polling supports cooperative cancellation, including
at the final return boundary. Managed execution is an expensive observation,
without a project snapshot; generic engine apply retains its normal journal,
with no content revision or undo checkpoint. Shared-shell invocation preserves
project/history. MCP's generic-op confirmation remains unchanged.

New files own protocol records and the assessment kernel/synthetic Rust
regressions. Existing dispatch, shared parser/introspection, registry, result
constructors, input-path summaries and managed observation classification are
changed only for this operation. Glossary-derived maturity digest projections
must be refreshed without upgrading scientific evidence. No committed fixture
or dependency/lockfile change is required; all FASTA tests generate documented
hand-crafted bytes in private temporary directories.

Verification on macOS arm64, 2026-11-04, in the uncommitted worktree based on
`89df315a80e0d3ae02e5685efd7a94e3895fbaa7`, using
`rustc 1.100.0-beta.7 (fd8e211ae 2026-10-17)`:

```sh
cargo test -q --offline --no-default-features --lib rt_primer_assessment_ -- --test-threads=1
cargo test -q --offline --no-default-features --lib primer_trace_ -- --test-threads=1
cargo test -q --offline --no-default-features --lib primer_payload_ -- --test-threads=1
python3 -m unittest scripts.test_primer_trace_contract scripts.test_capability_maturity scripts.test_tutorial_checkouts -q
python3 scripts/check_documentation_status.py
python3 scripts/manual_docs.py
python3 scripts/generate_capability_maturity.py --check
git diff --check
```

The Rust filters pass 14 new assessment/discovery/managed-policy tests, 29
trace regressions and 38 payload regressions. The 105 Python tests pass,
including real temporary LF/CRLF checkouts and manual/projection freshness.
Scoped rustfmt passes for the three new Rust files. A temporary standalone
Cargo program depending only on the current `gentle-protocol` calls
`render_gui_cli_mcp_parity_matrix_markdown()` and matches the projection
byte-for-byte; it does not depend on an older prebuilt docs executable.

The required default-feature `cargo check -q --offline` passes with temporary
Swift/Clang caches and elevated tool access:

```sh
CLANG_MODULE_CACHE_PATH=/private/tmp/gentle-rt-primer-swift-cache \
SWIFT_MODULECACHE_PATH=/private/tmp/gentle-rt-primer-swift-cache \
SWIFTPM_MODULECACHE_OVERRIDE=/private/tmp/gentle-rt-primer-swift-cache \
cargo check -q --offline
```

The initial outer-sandbox attempts failed in Apple dependency Swift builds,
first on the default user module cache, then on nested `sandbox-exec` after
redirecting caches. These were build-environment failures, not GENtle compiler
errors; no dependency, build-script or lockfile changes were made to bypass
them. Headless test linking warns about the large existing `__eh_frame`
compact-unwind limit; all selected tests still pass. Session-close reports
four OK, two reminders (intentional dirty/generated files and manual fidelity),
zero failures. All changes match this slice and its required integration/docs
projections; prior `outputs/` is untouched.

Native Windows/Linux, actual GUI/MCP/JS/Lua process/VM execution, full-suite,
authentic-reference performance and scientific acceptance remain unrun and
unclaimed. Local callback and managed-policy tests are not native worker or
responsiveness acceptance. No push, CI dispatch or publication occurred.


## Second-Slice Implementation Ledger

The owner authorized per-isoform CDS/UTR/absent classification and candidate
comparison using the existing shared transcript geometry, without waiting
for Claude. The scope is a separate `CompareRtPrimerCandidates` /
`primers compare-rt-candidates` request/report, not a changed first-slice v1
format, a new designer or a GUI-only calculator. Exact request and result
semantics are in the [comparison contract](cli/rt-primer-assessment.md#compare-candidates-across-annotated-isoforms).

The request embeds the existing prospective assessment request and explicit
loaded transcript sources. The comparison composes the background assessment
once through its existing method, preserving its own request/report hashes
and named missing-reference evidence. It does not accept supplied saved
counts as verified current analysis. Existing `TranscriptCaptureSource`,
coverage-universe resolution, mature transcript template/CDS construction and
source-coordinate mapping provide the transcript interpretation. Source
bindings retain complete DNA/annotation hashes and selected coverage
resolution; each member retains its cDNA hash and annotated CDS intervals.
Optional expected source-record hashes reject stale annotations as well as
changed bases.

Each candidate/member pair keeps independent presence and region evidence:
all canonical exact full-tail complement occurrences, ambiguity-compatible
possible occurrences and explicit absence within the selected scope. Members
without a site remain in the denominator. Multiple occurrences do not become
multiple covered isoforms. CDS uses its annotated interval union, not its
outer envelope; empty CDS and interior gaps remain unknown. Unknown bases
make a binding's classification unknown, while mixed means a crossing among
known CDS/UTR classes. Biological transcript orientation and shared mapping
preserve exon-junction and reverse-strand source projections. Ambiguous or
conflicting annotation resolution fails rather than choosing a gene-wide CDS.
Comparison-only guards refuse empty mature transcripts, multiple explicit
CDS qualifiers, more than 4,096 nonempty comma tokens, and mixed valid/malformed
explicit CDS tokens; entirely malformed/empty CDS remains unknown via the
existing resolver. No shared parser policy changes. Each binding maps only
the binary-selected overlapping exon-segment slice through the shared source
mapper, keeping repeated exact sites from rescanning unrelated exons.

Comparison rows retain input order, exact/possible-only/absent member counts
and overlapping per-region member counts, alongside existing tail/full-oligo
length, GC, descriptive Tm, homopolymer and complementary-run heuristics.
The embedded assessment's terminal matrix remains reference-scoped and joined
by primer identity. There is no combined specificity score, automatic winner,
coverage-universe shrinkage, readiness attachment or extrapolation from the
old-primer-selected cohort. The shared Tm model is disclosed;
`secondary_structure_status: not_assessed` keeps hairpin energies and new
Primer3 execution separate from inexpensive complementarity heuristics.

Admission bounds are 1–16 sources, 256 transcripts, one million mature cDNA
bases, 16 million loaded DNA bases and 100,000 source features. A cumulative
annotation preflight runs before source hashing or shared view construction:
at most 16,384 location nodes plus nonempty exon/CDS qualifier tokens, 16 MB
of feature-qualifier/external-location text and location depth 32. Its
conservative geometry estimate sums `U² × max(S, 1)` per source binding, capped
at 25 million units; `U` is the number of unique source/exon/CDS interval pairs
from splicing-seed features and `S` is the seed-feature count. These are
source-wide admission bounds, including repeated bindings, so the later
256-transcript/one-million-cDNA checks are not relied on to bound earlier
Splicing Expert construction. A narrow transcript selection does not bypass
these conservative loaded-annotation limits.
Comparison work defaults to 100 million base comparisons and 10,000 bindings,
with hard ceilings of one billion comparisons and 50,000 bindings. These limits do not
relax the embedded prospective assessment's independent budgets. Cancellation
and failures return no partial report. Managed comparison needs a project
snapshot for its loaded sources; it must not inherit the earlier file-only
assessment's project-free observation policy. Shared-shell execution remains
read-only; ordinary generic engine apply retains its journal without a new
content revision or undo checkpoint. No new lifecycle/publication exception
is needed.

### Second-Slice Local Verification

On macOS arm64, 2026-11-04, in the uncommitted worktree based on
`89df315a80e0d3ae02e5685efd7a94e3895fbaa7`, with the same beta Rust toolchain
recorded above:

```sh
cargo test -q --offline --no-default-features --lib rt_primer_comparison_ -- --test-threads=1
python3 -m unittest scripts.test_primer_trace_contract scripts.test_capability_maturity scripts.test_tutorial_checkouts -q
```

All 14 new comparison tests pass: both strands, CDS/UTR endpoints and crossings,
unknown/gapped CDS, junction-spanning source projection, possible/absent/exact
sites, repeated-site/member counts, explicit universe refusal/denominators,
annotation-only source hash changes, metrics/order, identical embedded
assessment bytes, work/output/cancellation/schema limits, malformed CDS,
conflicting/clipped geometry, pre-view complexity and typed/shared-shell parity.
The 105 Python checks pass, including temporary LF/CRLF checkouts. Current
`gentle-protocol::render_gui_cli_mcp_parity_matrix_markdown()` from the temporary
standalone renderer matches the updated parity document byte-for-byte.
Documentation/manual/maturity freshness, new-file rustfmt and whitespace pass.

Reuse the freshly compiled headless libtest executable for the existing capture
and assessment filters, with Cargo's configured test-thread stack:

```sh
RUST_MIN_STACK=16777216 target/debug/build/GENtle/40cefbf0f8c3f2c0/out/gentle-40cefbf0f8c3f2c0 transcript_capture_ --test-threads=1
RUST_MIN_STACK=16777216 target/debug/build/GENtle/40cefbf0f8c3f2c0/out/gentle-40cefbf0f8c3f2c0 rt_primer_assessment_ --test-threads=1
```

These pass 9 capture and 14 first-slice assessment tests. The executable's
SHA-256 is `236fb186d15c06af29d10966604027c20a0121f357c64aedabcfdc8bbb83f79d`.
A redundant Cargo rebuild after formatting-only edits was interrupted before
using this already-built binary; no functional change followed its compilation.
Initial direct libtest launches omitted `RUST_MIN_STACK` and aborted with stack
overflow. Supplying the existing `.cargo/config.toml` value fixes the launcher
boundary; this is not a source regression or a waived assertion. The existing
large-`__eh_frame` linker warning remains visible.

One broader headless regression was run on that same executable:
`engine_shell::tests::execute_introspect_capabilities_projects_full_registry_with_fact_annotations`.
It fails solely on `("primers_preflight", "registry_only")`; the new RT
descriptors are fact-annotated. The `primers_preflight` registry entry exists
at HEAD and its descriptor wiring is untouched by either RT slice. This is a
separate follow-up, not a green full-registry verdict or a baseline test run.
Do not repair unrelated discovery metadata in this diff.

The required default-feature `cargo check -q --offline` passes on the latest
source with the same temporary Swift/Clang caches and elevated tool access
shown in the first-slice receipt. No dependency/toolchain changes or build-script
bypasses were needed. Session-close has four OK, two expected reminders (intentional dirty/generated files and
manual plan fidelity), zero failures. Root manifest/lockfile and prior
`outputs/` are unchanged. The helper extraction and required adapter/projection
updates are intentional parts of this slice, without unrelated refactoring.
Native Windows/Linux, actual GUI/process/VM execution, full-suite, performance
and scientific acceptance remain unrun; no ordinary `.13` release gate changes.

Retrospective junction inspection, abundance/expression, thermodynamic
secondary structure, constrained redesign, native forms and order readiness
remain later work. No private reference/FASTQ replay, provider call, order,
push or publication belongs to this slice.


## Third-Slice Implementation Ledger

The owner authorized the primer–payload junction inspector after the isoform
comparison slice, without waiting for Claude. This slice projects existing
saved evidence; it does not rerun original FASTQ, PAF, references/annotations,
trimming or alignment, infer an RT mechanism, or revise earlier report schemas.

`InspectPrimerPayloadJunctions` / `rna-reads primer-junctions` returns the
separate portable `gentle.primer_junction_inspection.v1` through
`OpResult.primer_payload_junctions`. Its request supplies the exact saved Phase 2
result binding and raw-byte-bound report, optional preview filters/pagination
and existing `PrimerPayloadBudgets`. See the
[inspection contract](cli/rt-primer-assessment.md#inspect-saved-primerpayload-junction-evidence)
for the request, fields and counting units.

The inspector reuses saved trace metadata validation, prepared-query/FASTA and
trace-row binding, original anchor geometry and the existing Phase 2 evidence
population checks. A retained-primary preparation still requires its saved
retained payload FASTA; a historical replay mode does not cause a fresh FASTQ
read. Every Phase 2 sidecar is verified to EOF before a capped preview returns.
Existing `evidence::aggregate` recomputes partner profiles and edge clusters
and checks saved arithmetic instead of adding another grouping engine.

Joined rows retain original-read trace intervals, saved primer-oriented query
geometry and masks, all primer hit alternatives, payload segments, seed
mismatch/reason evidence and unassessed dispositions. Each joined row also
exposes `payload_previews`: query ID to the first at most 64 saved, masked,
already primer-oriented query bases. This reuses the verified query bytes,
without another orientation transform, new trimming or reference-flank
reconstruction. Approximate full anchors
use observed alignment boundaries. Core-only or handle-only observations do
not acquire an invented full-primer trim, handle observation or target.
A query's representative hit does not classify every primer sharing its
geometry. Partner and primer ambiguity sets remain intact; missing literal
reference flanks are not reconstructed from mismatch counts.

Distinct reads use `(file_id, record_index)`, with source-file hashes retained;
read names and anchor observations are separate units. Duplicate source-file
hashes in the selected verified cohort are refused. The cohort count includes
handle-only and excluded/unassessed trace records. Saved partner-profile and
edge-cluster scopes, overlapping denominator buckets and ambiguity restrictions
are preserved; no cross-reference/group sum is a molecule count. Preview
filters and pagination do not alter whole-cohort profiles or denominators.
Source trace/query/alignment outcomes, control errors and named missing units
remain visible independently of successful inspection.

Defaults are 64 MB cumulative saved input, 8 MB response, 4,096 source
records/queries and 1 MB lines, with a 128-alignment/query ceiling. Inspector
hard limits include 128 MB input, 8 MB output, 4,096 records/queries, 16 batches
and references, and preview rows/offset at most 1,000/4,096. Each evidence,
partner-profile and breakpoint-cluster sidecar is independently bounded by
`max_records` rows. Known report/sidecar/query-manifest bytes are reserved
before trace validation receives its residual metadata allowance; the query
validator checks the full retained dependency total before reading its data
artifacts. Exceeded evidence/output budgets refuse rather than silently
truncate counts. Callback cancellation is checked before/after work and during
streamed verification; saved-query validation remains a bounded synchronous
phase. Managed execution is an expensive, project-free observation using the
existing policy. No external file, content revision, undo checkpoint, order
readiness or new artifact-publication exception is introduced; ordinary generic
engine apply retains its journal, and shared-shell inspection does not.

### Third-Slice Local Verification

On macOS arm64, 2026-11-04, in the uncommitted worktree based on
`89df315a80e0d3ae02e5685efd7a94e3895fbaa7`, with the previously recorded
`rustc 1.100.0-beta.7`:

```sh
cargo test -q --offline --no-default-features --lib primer_junction_inspector_ -- --test-threads=1
python3 -m unittest scripts.test_primer_trace_contract scripts.test_capability_maturity scripts.test_tutorial_checkouts -q
```

All 15 new tests pass: forward/reverse previews, same-geometry exact/approximate
primer alternatives, core-only boundaries, distinct record versus repeated
anchor/read-name counting, whole partner-equivalence sets, no-anchor/excluded
records, page/filter-independent denominators, deleted original FASTQ/reference/
annotation/PAF, required retained-artifact refusal, hashes and complete evidence
population, synchronized rehashed aggregate/accounting/internal-edge tampering,
saved cancelled/control-error outcomes, budgets/cancellation/schema and typed/
shared-shell/@file parity. The 105 Python checks pass, including temporary
LF/CRLF checkouts. New-file and changed query-validator rustfmt, documentation,
manual/maturity and whitespace checks pass.

The same freshly compiled headless libtest executable has SHA-256
`6c7d4964e426148da0adfe82e22bdd6fd1148518861bc97f3eed426a9732bccd`.
With `RUST_MIN_STACK=16777216` (the existing Cargo configuration), direct
`primer_payload_` and `primer_trace_` filters pass 38 and 29 existing tests,
respectively. No functional source changed after the build. The existing
large-`__eh_frame` linker warning remains visible, without a performance claim.
The current protocol crate's canonical Rust parity renderer matches the
projection byte-for-byte using the temporary standalone renderer described
above, not an older prebuilt documentation executable.

The required default-feature `cargo check -q --offline` passes with the same
Swift/Clang temporary-cache and elevated-access workaround recorded above.
Root manifest/lockfile and the unrelated pre-existing `outputs/` are unchanged.
Session-close reports four OK, two expected reminders (intentional dirty/
generated files and manual plan fidelity), zero failures. The new files,
legacy-validator extraction and required operation/adapter/projection wiring
are scoped to this inspector. The broader `primers_preflight` discovery issue
from the previous ledger is unchanged; its failing gate was not rerun or
represented as green here.

Native GUI/process/VM execution, Windows/Linux, large real cohorts,
candidate-specific reference-flank rescoring, abundance, redesign and RT
mechanism/experimental acceptance remain separate. No private-data replay,
provider call, order, push or publication belongs to this slice.


## Fourth-Slice Implementation Ledger

The owner requested committing the previous three slices, then constrained
redesign for explicit shared-transcript or 3′-UTR groups without waiting for
Claude. Commit `fee72185` contains those three completed slices. The redesign
is a separate subsequent worktree change, not a claim that the earlier commit
contains it. The original plan and uninvoked read-only Claude prompt remain
historical planning evidence.

This is biological-extension class B: one small shared orchestration operation
composes canonical transcript window enumeration, existing tail/full-oligo
metrics, source/CDS/member preparation, per-isoform comparison and prospective
background assessment. `RedesignRtPrimers` / `primers redesign-rt` returns
`gentle.rt_primer_redesign.v1` in `OpResult.rt_primer_redesign` from the strict
`gentle.rt_primer_redesign_request.v1`. See the [complete request and semantics](cli/rt-primer-assessment.md#redesign-rt-primers-for-an-explicit-transcript-group).
Existing assessment/comparison schemas and their evidence remain unchanged.

Each source has an explicit transcript selection. The nonempty required-member
subset is independent of the complete comparison universe; comparison-only
members retain their denominators without expanding the generation obligation.
A caller-selected required anchor and half-open mature cDNA interval define
where antisense tails may originate. Shared mode requires a canonical exact
full-tail site per required member. Strict 3′-UTR mode requires every exact or
ambiguity-compatible site in each required member to be wholly in its known
annotated 3′-UTR. Unknown CDS, annotation gaps, boundary crossings and ambiguous
possible sites outside that UTR cannot satisfy the contract. Both source strands
and splice-spanning sites use the existing shared geometry.

Generated full oligos keep the fixed common handle and terminal target core as
separate components. Required inclusive tail length/GC/Tm constraints and
homopolymer/self/general-3′ complementary-run ceilings reuse shared heuristics;
run ceilings constrain tail and full oligo. No independent structure or
thermodynamic engine, native form, reference acquisition, aligner or RT-mechanism
claim is added.

Generation exhaustively enumerates admitted anchor windows in start-coordinate,
then length order, skips ambiguous words and deduplicates tails. The first
qualifying distinct tails fill the capped shortlist; all qualifying tails are
counted, and additional tails are explicitly background-unassessed. Candidate
background failures do not refill or rerank that shortlist. Declared seed limits
use existing all-site physical-window marginal mismatch histograms including
intended sites. Transcript geometry never infers an intended-site exclusion.
Proven over-limit counts fail; otherwise absent/ambiguous evidence remains
unresolved. Both failed and unresolved candidates stay visible with the embedded
unchanged comparison and assessment. No seed limits means a vacuous pass, not
specificity or order readiness.

Request admission validates metadata before generation. Defaults are 50,000
windows, 100 million conservative generation work units and 16 candidates;
hard ceilings are 100,000 windows, one billion work units and 16 candidates.
The normalized serialized request/response ceilings are one/eight MB, respectively. At most 64 seed
limits and 256 required members are accepted; the complete handle-plus-tail is
at most 64 bases, with a tail of at least 10. Shared source/annotation/materialized
transcript preflight, comparison and background budgets retain their independent
existing limits. Work includes canonical intake, conservative metric cost
`2 × L³ + 8 × L²` per distinct tail for complete oligo length `L`, and
full group-member matching. The default work budget can exhaust near 500
44-base oligos before member search; callers should narrow the declared anchor
or increase the explicit work budget within the ceiling. The enumeration
ceiling is not a throughput promise. Excess refuses rather than silently truncating
counts or becoming an empty biological result. An exhaustive empty result keeps
resolved source/member evidence but has no comparison and does not open
background files. Callback cancellation returns no partial redesign report.

The operation writes no files or project content. Shared-shell execution keeps
history unchanged; generic apply retains its normal journal without a content
revision or undo checkpoint. Managed execution uses the ordinary explicit
project snapshot/execution lane, not a project-free observation or publication
exception. Required routing, discovery, glossary and generated parity/maturity
updates belong to this slice. Maturity digest refresh does not upgrade evidence.
The unrelated pre-existing `outputs/` remains outside scope.

### Fourth-Slice Local Verification

On macOS arm64, 2026-11-04, with `rustc 1.100.0-beta.7
(fd8e211ae 2026-10-17)`, in the worktree based on
`fee7218507bf9a54be38275761e6c90e3994d5f3`:

```sh
cargo test -q --offline --no-default-features --lib rt_primer_redesign_ -- --test-threads=1
python3 -m unittest scripts.test_primer_trace_contract scripts.test_capability_maturity scripts.test_tutorial_checkouts -q
```

All 17 focused checks pass: 16 new synthetic redesign tests and the project-
snapshot execution-policy test. They exercise both strands, different shared
site coordinates, splice-junction source projections, required versus comparison-
only members, strict UTR and ambiguity rules, fixed-handle metric constraints,
sequence deduplication and shortlist ordering/no refill, positional seed limits,
missing/all-ambiguous/short-seed evidence, empty unopened backgrounds, annotation-
only stale hashes, resource/cancellation refusal, determinism and typed/shell/
quoted `@file` parity without content/history changes.

The initial focused run passed 16 of 17: an all-ambiguous verified reference
correctly remained unresolved but its saved ambiguity count was omitted from
the redesign projection. The fix retains that count independently of the
assessment's completed-cell status, while missing references and unsupported
seed lengths retain `None`. The final rerun passes all 17; no earlier assessment
schema or cell-completion rule changed. Read-only source audits covered the
composition, strict UTR rules, bounds and this correction; no Claude was invoked.

The freshly compiled headless libtest executable has SHA-256
`47bb9c07955961b9cda234ba92d4a68ecee69934dce271b1e40746662639e5f5`.
With `RUST_MIN_STACK=16777216` from the existing Cargo configuration, its
`rt_primer_comparison_`, `rt_primer_assessment_`, `transcript_capture_` and
`terminal_exon_rt_primer_` filters pass 14, 14, 9 and 5 tests, respectively.
These exercise the unchanged contracts and prior callers of the extracted
shared geometry/window helpers. No functional source changes follow this build.
The existing large-`__eh_frame` linker warning remains; it is not a new
performance acceptance result.

All 105 Python checks pass, including actual temporary LF/CRLF checkout
boundaries. The current local protocol crate's canonical Rust renderer produced
the updated parity projection, rather than an older documentation executable.
Scoped source rustfmt, documentation/manual/maturity freshness and whitespace
checks pass. Root manifest/lockfile and the unrelated `outputs/` are unchanged.
The generated maturity changes refresh only glossary input digests.

The required default-feature `cargo check -q --offline` passes with the same
temporary Swift/Clang caches and elevated-access workaround recorded in the
first ledger. No dependency, lockfile, build-script or toolchain change was
needed. Session-close reports four OK, two expected reminders (intentional dirty/
generated files and manual plan fidelity), zero failures. The scope review
confirms the new contract/kernel/tests, necessary shared helper extractions,
adapter/discovery wiring and required projections/documentation are intentional;
there is no unrelated feature or catalog expansion. The broader unchanged
`primers_preflight` discovery issue recorded in the second ledger was not rerun
or represented as green; this is not a full-registry/full-suite verdict.

Native GUI/process/VM execution, Windows/Linux, real-cohort/reference performance,
abundance, candidate-specific retrospective rescoring, calibrated secondary
structure, experimental validation and ordering acceptance remain separate.
No private-data replay, provider/Claude call, order, push or publication is
part of this implementation slice.
