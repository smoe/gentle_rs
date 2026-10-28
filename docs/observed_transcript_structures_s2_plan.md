# Observed Transcript Structures S2: SAM And Annotation Preparation

Status: **S2.0a/b retained; private S2.1a implemented with scoped verification**.
S2.0b adds a private shared core and decoded, in-memory SAM normalization only.
The owner subsequently authorized only S2.1a's pure GTF normalizer and focused
local Rust tests/check; its working-diff verification is recorded below.
Annotation file ingress/publication, operations/adapters and later slices remain
proposals, not authorized by either private step.
Prepared 2026-10-27 at local `main` = `gentle_rs_2_main` =
`1804605f1ba42369f943213414d79fd9a40013cb`. No fetch or remote-tip claim.
Tracked files were clean; unrelated untracked `output/` is untouched.
The [roadmap](roadmap.md) owns priority/release scope. This is a bounded expansion
of [the umbrella plan's S2](observed_transcript_structures_plan.md#slice-boundaries-and-remaining-work),
not a new release gate. The [original Claude prompt](observed_transcript_structures_s2_claude_review_prompt.md)
was reviewed; the reconciliation below is Codex's, not a fresh Claude sign-off.
S2.0a did not authorize a production refactor or SAM ingress. After its commit
and rebase, the owner separately authorized bounded S2.0b and focused local
Rust tests/check. No file ingress, new operation, GUI, dependency or push.

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

### S2.0b Implementation, 2026-10-28

Base: `a52c7b7473ef5e175f928c9504b98aadd435d99b` on `gentle_rs_2_main`,
after S2.0a's commit/rebase onto local main `b0521123`. No fetch/remote-tip claim.
The implementation was verified as a working diff, not a new public SAM capability:

- `core.rs` owns the one observation/DAG/reconciliation/comparison path over
  private typed keys, borrowed query/segment views and a minimal context.
  `kernel.rs` preserves legacy admission and projects back to the unchanged
  public v1 records. Legacy context serialization and observation/physical-chain
  digest recipes stay at their existing boundaries. The public DTOs, Phase 2
  walker, DAG, assertions and S2.0a expectations are not rewritten.
- `sam.rs` accepts already decoded fields and explicitly grouped units, not
  SAM text, files or optional tags. A future parser must validate text/tag
  syntax and provide exact row accounting; that boundary is still unimplemented.
  No QNAME grouping, original-read identity, invented PAF identity, NM/MD/AS
  inference or aligner is introduced. It has fixed private context/observation
  domains and reuses the unchanged physical-chain key. Private normalized
  results are not a versioned SAM request/report/publication contract.
- Normalization checks terminal `H/S`, unclipped axes, original-query reverse
  orientation, reference dictionary/POS and optional SEQ/QUAL lengths. Preserve
  `M/=/X/I/D/N`; only `N` creates skips. Under explicit unpadded admission,
  `P` consumes no coordinates and its operation count remains a diagnostic.
  Hard/soft counts remain stored-orientation counts; projected intervals use
  the original-query axis. MAPQ 255 is absent, while 0 remains reported zero.
- Paired/mate, padded-reference, CG/placeholder, missing CIGAR, no supported
  bases and overlong-axis rows make their entire unit unavailable. Unmapped
  rows also make the unit unassessed; zero-row units are `NotProvided`. No
  retained good row may contribute possible support for such a unit. Malformed
  geometry or contradictory declarations fail the entire call. Missing axes
  remain absent, not manufactured zero-length queries.
- All units/rows/operations are admitted before normalized copies, even for
  unavailable units; exact caps and cooperative control errors are tested.
  Dictionary/reference authority and search completeness remain producer
  declarations. Missing masks and incomplete search stay independent axes.
  Reference genes/sources remain explicit, with two opposite-strand synthetic
  loci. There is no source selection, transcript canonicalization, protein,
  biological novelty, molecule-count or natural-RNA-end assertion.
- `sam_tests.rs` contains inline synthetic provenance and an independent
  basewise geometry oracle, asymmetric and Cartesian clipping cases, SAM/PAF
  physical-key parity, configurable genes, repeated names/row roles, whole-unit
  refusal, missing axes/masks, alternatives, nonoptimal observations, strict
  data errors and resource/control tests. No committed external data/fixture
  or newline policy is added.

Verification of this working diff on macOS/aarch64 uses `CARGO_INCREMENTAL=0`
and `--locked --offline --no-default-features --lib -j1`; every test command
uses `-- --test-threads=1`:

| Command/filter | Result |
| --- | --- |
| `cargo test ... observed_transcript_structures` | 41/41, including 13 new SAM tests and all 28 existing tests |
| `cargo test ... s2a_` | 6/6 exact byte/identity/work/error and S1 publication guards |
| `cargo test ... observed_transcripts_` | 16/16 engine/shell/MCP-unit/lifecycle tests |
| `cargo test ... primer_payload` | 48/48, including the unchanged ten-case Phase 2 serialized-row/work oracle |
| `cargo check -q ...` | Pass; existing headless dead-code warnings |

Filters overlap; no ignored tests or regenerated preservation expectations in
these Rust runs. The test build retains the existing large `__eh_frame` linker
warning, not a measured exception/runtime-performance verdict. Cargo/rustc,
effective unoptimized profile, incremental setting and 16-MiB thread-stack
configuration match the S2.0a execution identity above. Cargo.lock remains
`40ae5b6c9ba74ba080d39be11b4e57ad2dbc020dba5ecad7222e22f639f84918`.
The 18-path Rust/Cargo/config manifest at
`/tmp/gentle-s2b-source-manifest.sha256` has SHA-256
`2e63cea03130af71dcf51a634c671ec3808d958163f6cc358fc5770b92bd2b89`;
the test binary at
`target/debug/build/GENtle/002e669310863d66/out/gentle-002e669310863d66`
has SHA-256
`f099940ff3cc411014feac07f512ea4d6a0a28563eca30d3838a22eac21aa028`.
These are local working-diff receipts, not committed fixtures or release assets.
Local main subsequently points at the same S2.0a base `a52c7b74`; no later
source change, fetch, merge or push is inferred.

The owner subsequently requested committing S2.0b and rebasing onto local main
`98f69995`, whose added commit changes lifecycle-preview tests and documentation,
not these RNA paths. Preserve the pre-commit source/binary receipts above;
integration does not relabel them as execution of the later tree. No fetch,
push or later-slice implementation is requested by this integration action.

Documentation status passes; 34 fast documentation/contract cases finish with
33 passes and one expected binary-dependent skip, including LF/CRLF documentation
coverage. Scoped rustfmt on the six changed Rust paths and `git diff --check`
pass. Whole-workspace `cargo fmt --check` still fails in the same 12 unrelated
baseline files (local `/tmp/gentle-s2b-fmt.log`); they are untouched. Session-close
reports 4 OK, 2 warnings (intentional working diff/untouched `output/`, manual
plan fidelity), no failures. Scope matches the separately authorized private
S2.0b contract; no adjacent behavior change or extra cleanup is included.

Full suite, default-feature/GUI check, actual CLI/VM processes, Linux/Windows,
native, real-data/scientific and release acceptance are unrun. Retain S2.0a's
historical execution identity separately. S2.1/S2.2/S2.3, text/file/optional-tag ingress, annotation gap and
duplicate-accession policies, qualification-split public SAM catalogs, receipts,
operations/lifecycle/adapters and saved publication stay separate.

### S2.1: Explicit Annotation Preparation

#### S2.1 Review Reconciliation, 2026-10-28

The owner supplied Claude's read-only S2.1 review against `c5a9d63f`.
Codex rechecked that same source: `main` and `gentle_rs_2_main` both point at
`c5a9d63fcddaa78cc02a791b6001081b39af80b6`; only unrelated untracked `output/`
is present. No fetch, new Claude invocation or runtime verification. The
following **revised Codex proposal replaces the earlier S2.1 outline**; it is
not implementation authorization or a fresh Claude/scientific sign-off.

Original recommendation: GTF-first preparation into the existing annotation
geometry plus a receipt, then SAM intake and saved interpretation separately.
The supplied review identified missing outcome/dictionary contracts and
underspecified grammar, source scope and lifecycle details.

| Finding | Source check and revised proposal |
| --- | --- |
| H1: empty/refused sources | Confirmed: `comparison::sources` permits zero transcripts; Phase 2 `reference::catalog` refuses them. Neither legacy enum nor `core::Context` represents a refused source. Publish a receipt only for gene-level refusal/absence, never an empty v1 annotation. Later S2.3 needs a SAM-only resolution projection; do not extend legacy `GeneResolution`. |
| H2: dictionary absent on disk | Confirmed: `sam::SamInput.dictionary` is only an in-memory map. Specify a strict ordered dictionary document below, rather than infer competitor lengths from one selected gene. |
| H3: transcript-ID reuse | Accept bounded whole-file indexing and required/not-required transcript-row policy. Qualify the review's stronger claim: transcript rows provide evidence of multiple declarations, not proof of biological reuse or complete exon coverage. Without them, consistent same-locus reuse may be undetectable; with them, missing internal exons still need a coverage nonclaim. |
| M1-M4: grammar, locus and scope | Specify the strict subset, an `ambiguous_gene_locus` receipt-only outcome and exact-ID missing-gene semantics. EOF verifies supplied bytes, not completeness of a provider annotation or correctness/stability of its IDs. |
| M5-M6: legacy compatibility | Accept mandatory annotation/receipt pairing for future S2 intake. Bare v1 geometry may still pass old consumers; do not claim those consumers enforce new policies. Keep touching-exon validation unchanged in both legacy paths and add preservation tests. |
| L1-L2: budgets/helpers | Use preparation-specific budgets, not irrelevant read/PAF fields. File helpers are reusable from a child beside `observed`, but not all are drop-in: `read_json` and `publish` take `PrimerPayloadBudgets`. Any narrow limit-helper extraction must preserve the old wrapper's limits, bytes and errors. No blanket visibility/refactor change. |
| L3-L4: paths/lifecycle | Accept explicit input/output/control collision checks and named prepare/cancel approval/no-save/managed policies. Existing `staging`/`publish` are not a complete collision or rollback contract. Do not import Plan-only writer assumptions. |

Codex also verified the genome-index parser beyond the review's module-level
claim: `parse_tabular_annotation_gene_records_with_progress` skips malformed
rows/UTF-8, while `parse_annotation_attributes` lowercases keys and keeps the
first value. Do not change or directly reuse that permissive parser.

#### Source, Dictionary And Identity Contract

Propose one strict versioned `PrepareObservedTranscriptAnnotation` request for
one source and configurable `{source_id, gene_id}`, with explicit GTF format,
`reference_id`, assembly, declared FASTA digest, input/dictionary path-hash-byte
bindings, policies, budgets and fresh output directory. Input scope is required:
`whole_source_file | producer_excerpt`, with optional declared derivation/source/
notice metadata. Neither value authenticates or proves whole-provider coverage.
No species, gene, strand, assembly, canonical model or release default.

Specify `gentle.reference_contig_dictionary.v1`: `schema`, `reference_id`,
`assembly`, `reference_fasta_sha256` and ordered `contigs: [{name, length}]`.
Refuse unknown/duplicate fields, empty dictionaries and duplicate exact names.
Names are nonempty UTF-8, at most 4096 bytes, not `*`, without whitespace or
control characters; do not rewrite case/aliases. Lengths are positive unsigned
integers with checked conversion to the consuming geometry's integer range.
Preserve list order in the bound document; use an exact-name lookup internally.
Check reference metadata equality across request/dictionary/output. Bound its
full bytes and entries independently of the selected annotation contig set.
Per-contig digests, `M5`, FASTA-index imports and inferred lengths stay outside v1.

For prepared geometry, retain `gentle.primer_payload_annotation.v1` unchanged.
Set `release = snapshot-sha256:<64-hex-raw-GTF-digest>`. Carry the request's
nullable `declared_provider_release` verbatim in the receipt, labelled as a
declaration, not a header-derived release. Provider labels are not interchangeable:
a UCSC track snapshot is not silently relabelled a direct NCBI GCF release.
Coordinate authority remains `declared`; hashing a dictionary does not verify
the genome FASTA. A different selection/policy can yield identical annotation
bytes but a different preparation-receipt identity.

#### Strict GTF Subset And Selection

Propose grammar `gtf_strict.v1`, not a claim to accept every GTF dialect:

- Count column-one `#` comments, including release-looking headers, without
  treating them as metadata authority. Reject actual blank/whitespace-only
  records in this subset; an ordinary final newline is not a blank record.
  Accept a final unterminated data line. Strip only LF/one terminal CR for
  parsing; hash unmodified bytes. Reject embedded controls and invalid UTF-8.
- Require exactly nine tab-separated columns. Coordinates are ASCII decimal,
  overflow-checked, with `1 <= start <= end <= dictionary length`; every data
  contig must exist in the dictionary, including unselected rows. Validate
  source/feature tokens (nonempty, at most 4096 bytes, no whitespace/controls),
  score as `.` or finite number, frame as
  `.|0|1|2`, and strand as `+|-|.`. Store exons as `[start-1,end)`.
- Parse case-sensitive ASCII attribute keys (`[A-Za-z_][A-Za-z0-9_]*`), separated
  from values by one or more ASCII spaces, with quote/escape-aware semicolon
  splitting. Require a terminating semicolon per attribute and nonempty key/
  value. Consumed IDs/versions are quoted strings; recognise escaped quote
  and backslash explicitly and reject unsupported escapes rather than silently
  normalizing them. Unknown attributes may contain quoted values or bare
  tokens, may repeat and are intentionally unused, not biological evidence.
- Consume `gene_id` on every data row; `transcript_id` must be absent on
  `gene` rows and present on other features in this declared subset. Consume
  optional `gene_version`/`transcript_version`; reject duplicate consumed
  keys even if equal. IDs remain exact, nonempty, at most 512 UTF-8 bytes,
  without controls/whitespace. Separate version attributes never alter IDs.
- Lexical, required-field, dictionary and coordinate failures anywhere refuse
  the request with no finalized report. For otherwise parseable unrelated
  models, recommend counting semantic conflicts with bounded examples, not
  refusing the selected gene. This proposed policy must be explicit in the
  grammar/receipt. It cannot hide a conflict linked to the selected gene.
- Keep a budgeted **whole-supplied-file** index from transcript ID to declared
  gene/contig/strand associations and transcript-feature declarations. Ordinary
  exon rows sharing one consistent transcript ID are not duplicates. Detect
  cross-gene/cross-locus collisions before selection discards other rows.
- For the selected gene, every row must have consistent locus/strand/linkage;
  require `+` or `-`, retain every supplied exon and sort ascending independently
  of row order or `exon_number`. Never clip a model to the requested locus or
  silently discard a failing transcript. Multiple loci/strands produce
  `ambiguous_gene_locus`.
- Require explicit `transcript_feature_rows: required | not_required`.
  Under `required`, every selected transcript needs one declaration; absence
  is `transcript_model_conflict` with a missing-declaration reason. A supplied
  declaration must be unique and its span equal the supplied exon hull.
  Repeated declarations are observational collision evidence, not a mechanism.
  Under `not_required`, disclose undetectable consistent same-locus reuse.
  Neither setting proves full transcript coverage in an excerpt.
- For selected rows, require consistent version values and all-or-none presence
  in each applicable scope: gene version across that gene, transcript version
  across rows bearing that transcript ID, not across gene-only rows. Missing/
  contradictory applicable declarations are model conflicts, not ID rewrites.
  Recommend refusing identical selected-model duplicate rows in v1 rather than
  deduplication; retain their count/reason. These technical policy recommendations
  are not real-provider or biological approval.
- Refuse overlapping/touching selected exons in preparation only. Require
  `minimum_reference_gap_bases >= 1` explicitly; gap equal to the minimum
  passes, a smaller positive gap refuses the whole selected gene. Never merge
  exons, snap boundaries or infer a mechanism from gap size. Real-data threshold/
  refusal acceptance remains an owner decision; synthetic tests use labelled
  hand-chosen values.

#### Outcomes, Receipts And Limits Of Reuse

Define a new preparation outcome, separate from execution status:
`completed`, `missing_gene`, `ambiguous_gene_locus`,
`duplicate_transcript_accession`, `transcript_model_conflict`,
`identical_duplicate_row`, `reference_gap_policy_refused`.
For collisions retain observed bases `multiple_transcript_rows`, `multiple_loci`
and/or `cross_gene_link`, not an inferred biological cause. Retain all bounded
reason counts and choose the primary refusal deterministically: collision,
ambiguous locus, identical duplicate, model conflict, then gap policy.
`missing_gene` means no exact `gene_id` match in supplied bytes, including
case/version-spelling differences; it is not proof of biological absence.
A matched gene without usable exon models is a model conflict, not missing.

Only `completed` publishes a nonempty annotation artifact. Every completed
scan with a gene-level absence/refusal publishes **receipt only**, never an empty
v1 annotation. A preparation receipt in `report.json` binds request/schema/
grammar version, raw GTF/dictionary bindings and verification status, reference
metadata, exact selector, source/scope/provider claims, policies, feature/index/
diagnostic counts, outcome/reasons, and the optional annotation artifact hash.
The receipt's canonical digest excludes only its own digest field; record its
raw file hash separately. Dictionary and all source bytes must be verified to
EOF before any completed gene outcome is published.

Cancellation/control failure before complete verification leaves gene outcome
and complete analysis/counts absent, not zero or `missing_gene`; any processed
diagnostics are explicitly partial. After a completed verified unit, retain its
analysis with the late interrupted execution status, following S1's distinction.
Data/hash/resource errors publish no final report. Define/test these states in
the request/receipt contract before implementing file publication.

S2.3 must consume a matching annotation/receipt pair or a typed receipt-only
unavailable source. Add a SAM-only resolution wrapper there; do not fabricate
an empty `AnnotationSource`, extend legacy `GeneResolution`/`AnnotationBinding`
or tighten old `sources`/`catalog`. Exact receipt-selector mismatch is
`selector_binding_mismatch`, not missing. Bare completed v1 geometry can still
be accepted by old validators: that is structural readability, **not** proof
that Phase 2/S1 enforce the receipt, gap policy or whole-reference coverage.

#### Bounded Implementation Order And Verification

The owner separately authorized S2.1a and focused local Rust verification on
2026-10-28, after reviewing the reconciliation. Its private working diff is
described below; S2.1b/c still require separate authorization:

1. **S2.1a -- pure private kernel and synthetic guards.** Implement private typed
   dictionary/policy/outcome views, a bounded incremental line consumer and the
   whole-file identity index under `observed_transcript_structures`. No file
   access, public schema, operation or adapter. Test normal exon repetition,
   two configurable opposite-strand genes, collisions, source scope, strict
   grammar, version/gap/refusal semantics and exact resource/control boundaries.
   Add direct guards that both unchanged legacy validators still admit touching
   exons; preparation must refuse them. Preserve every S2.0a/Phase 2/S1 guard.
2. **S2.1b -- versioned request/receipt/dictionary and file boundary.** Place
   I/O below `primer_payload`, beside `observed`; reuse suitable verified-reader,
   line, control and staging patterns. The existing `line()` strips terminators;
   S2.1a needs raw line framing for hashing, so do not reconstruct LF after
   decoding. Any raw-line helper must preserve legacy `line()` behavior.
   Use preparation-only raw/line/index/exon/
   work/memory/diagnostic/output limits, including unselected rows and comments.
   Raw GTF streaming and bounded dictionary JSON are separate limits: the
   existing S1/read-JSON cap is 128,000,000 bytes, not proof that the generic
   verified reader cannot stream larger sources. Do not raise any legacy cap.
   Initial preparation ceilings remain the explicit engineering proposals in
   Bounds below; larger real-provider files need a separately approved budget.
   Narrow helper extraction, if necessary, must preserve old signatures/behavior
   and avoid dummy query/PAF budget fields.
   Check raw path components before normalization; refuse input/output/control
   aliases and a cancellation marker within output or equal to either input.
   Check existing-file identities/resolved parents for aliases, regular-file/
   symlink and no-clobber rules; never treat lexical comparison as complete
   cross-platform collision protection. Close handles before publication.
   Reuse report-last publication only after checking its failure paths: existing
   `publish` can leave reserved/artifact directories after a later failure.
   No final receipt means no completed artifact set; clean only owned staging/
   reservations where safe or disclose orphaned output, never overwrite/reuse it.
3. **S2.1c -- named operation and thin adapters.** Register prepare plus a
   request-bound Cancel companion through shared engine/parser/shell, generic
   MCP, discovery/parity and named construction/no-save/write-approval/managed
   policies. Approval covers every receipt/marker write attempt, including
   refusal or cancellation. Preserve conservative workflows and project/journal
   boundaries. Do not add SAM interpretation routes, GUI or a second registry.

File-backed acceptance must cover LF/CRLF raw hashes and releases with identical
geometry, quoted paths with spaces, symlink/input changes, all output/control
collisions, hash mismatch, budget/cancellation/write failures and report-last
visibility. Check final directory contents, but do not infer write chronology
from directory iteration order. Same annotation bytes under two admitted policies
must have distinct receipt identities. Native Windows handles, actual CLI/MCP/
VM, real-data/scientific and release acceptance stay separately unverified.
GTF fixtures remain inline/synthetic initially; any later committed byte-bound
input requires provenance, scoped attributes and fast checkout coverage.

GTF-first is the recommendation, not authorization to add genePred, GFF3, BED12,
JSON sniffing, implicit conversion/network retrieval or guessed aliases. Later
genePred requires declared columns and a snapshot-bound stable gene-ID crosswalk;
a `name2` label alone is not that mapping. The archive needs explicit compatible
companions, not a promise that every retained file is already loadable.

#### S2.1a Implementation, 2026-10-28

Owner-authorized base: `c5a9d63fcddaa78cc02a791b6001081b39af80b6` on
`gentle_rs_2_main`, matching local main at start. The four-document review
reconciliation was already an intentional working diff; unrelated `output/`
remains untouched. No fetch, new model review, commit or push.

- `annotation.rs` and `annotation/gtf.rs` add only a private, incremental raw-line
  consumer. Private context/dictionary/policy/limit/outcome views are not a
  frozen schema or file parser. The caller supplies raw LF/CRLF or final
  unterminated lines; the kernel hashes exactly those bytes and requires an
  explicit finish. It cannot verify an original file, provider coverage or
  declaration. Failed/control-interrupted input poisons the consumer, preventing
  a later partial successful finish. No stream/file/control-marker is opened.
- Strict row/attribute admission and a bounded whole-input gene/transcript index
  preserve decoded IDs and separate version strings without aliases or suffix
  construction. Selected-linked cross-gene/locus collisions survive selection;
  normal exon repetition is accepted. Unknown strand is refused as missing
  strand evidence, not manufactured into a second known locus. Two declared
  opposing strands or different contigs remain ambiguity/collision evidence.
- Selected models retain all supplied exons in ascending order and use the
  unchanged `PrimerPayloadAnnotation` geometry with raw-snapshot release.
  Refusal/absence returns `annotation: None`, never an empty v1 annotation;
  unused conflicting models are counted, not silently repaired. Required versus
  optional transcript rows, duplicate-row refusal, version scope and gap policy
  follow the reconciled technical contract. No real-provider/gap policy or
  biological identity/completeness is approved.
- Results retain exact scope/selector/provider declarations, raw hash/bytes,
  feature/selection/index counts, version maps, entity/reason diagnostics and
  nonclaims. Diagnostic examples are capped separately from exact totals;
  truncating examples cannot improve the outcome. Counts across gene/transcript
  entities are not independently summable biological events. There is no public
  preparation receipt, persisted digest recipe or adapter projection yet.
- Preparation-specific limits cover raw/line bytes, all lines including
  comments, dictionary entries, indexed genes/transcripts, all exon rows,
  charged work, retained-state accounting and diagnostic examples. Accounting
  includes row digests, index/object/string charges, vector reserve allowances,
  output geometry and bounded parser scratch. The reported peak is this
  deterministic representation-accounting model, not allocator/RSS measurement.
  No legacy cap/profile/dependency is changed; input-owned borrowed context and
  dictionary storage are outside the consumer's retained-state accounting.
- Inline synthetic tests cover variable genes/opposite strands, source/scope
  nonclaims, LF/CRLF/unterminated input, quoting/Unicode/byte limits, strict data
  errors, whole-file collisions, versions, model refusals, exact/one-over limits
  and every bounded control poll. Direct tests of the unchanged legacy source
  validator and real bound-JSON Phase 2 catalog preserve touching-exon admission;
  only preparation refuses touching models. The catalog fixture is synthetic
  and does not replay FASTA/alignment or assert native-platform behavior.

Final working-diff verification on macOS/aarch64 uses `CARGO_INCREMENTAL=0`
and `--locked --offline --no-default-features --lib -j1`; every test command
uses `-- --test-threads=1`:

| Command/filter | Result |
| --- | --- |
| `cargo test ... observed_transcript_structures` | 58/58, including 16 new pure GTF tests, the new source-validator touching guard, all S2.0b/S0 tests and four kernel preservation guards |
| `cargo test ... s2a_` | 6/6, including exact S1 report/data/table bytes and kernel identity/work/errors |
| `cargo test ... observed_transcripts_` | 16/16 engine/shell/MCP-unit/lifecycle tests |
| `cargo test ... primer_payload` | 49/49, including the new real bound-JSON catalog touching guard and unchanged ten-case Phase 2 serialized-row/work oracle |
| `cargo check -q ...` | Pass; existing headless dead-code warnings |

Filters overlap; none of these Rust tests is ignored. There are 18 new tests in
total, including both legacy guards. Earlier runs passed the initial 15 pure
cases and then the 58-test filter; the final rerun above additionally covers
the final whitespace-only metadata/source-token admission checks. No frozen
expectation is regenerated. Existing dead-code and large `__eh_frame` linker
warnings remain separate from the passing results, with no performance verdict.

Execution identity: Cargo `1.100.0-beta.5` (`7552a885e`), rustc
`1.100.0-beta.5` (`8ea38001b9da1d6683ac5ebbc993cac89d3a09df`), LLVM 23.1.3,
`aarch64-apple-darwin`. Effective unoptimized profile, root debuginfo 1,
dependency debuginfo 0, 256 codegen units, incremental disabled by environment
and the existing 16-MiB test-thread stack remain unchanged. Cargo.lock SHA-256:
`40ae5b6c9ba74ba080d39be11b4e57ad2dbc020dba5ecad7222e22f639f84918`.
The 22-path Rust/Cargo/config manifest at
`/tmp/gentle-s21a-final-source-manifest.sha256` has SHA-256
`0c4368e8129ea25f3e36a53f666017b823fd8d4b88a92895dc7c8ad622322c62`.
Test binary `target/debug/build/GENtle/002e669310863d66/out/gentle-002e669310863d66`
has SHA-256
`c9ac125e35405ca1df4e002a789027ba564309a369c24be055717514c3ad0ef6`.
These local working-diff receipts are not committed fixtures or release assets.

Documentation status, 34 fast documentation/contract cases (33 passed, one
expected binary-dependent skip), local Markdown targets and whitespace pass.
Scoped rustfmt for all six changed Rust paths passes; whole-workspace fmt still
reports the 12 unrelated baseline files in `/tmp/gentle-s21a-fmt.log`, untouched.
Session-close reports 4 OK, 2 warnings for the intentional diff/untouched
`output/` and manual plan fidelity, no failures. Manual scope review confirms
the authorized pure consumer/tests plus documentation only. No other cache
cleanup, profile, dependency, public DTO, production validator or fixture change.

S2.0a expectations, public v1 DTOs and ten-case Phase 2 oracle remain unchanged;
old execution receipts are not reattributed to S2.1a. Full suite, default-GUI
check, actual CLI/VM processes, native Windows/Linux, real-data/scientific and
release acceptance are unrun. S2.1b/c file/JSON/tag ingress, public receipts,
publication/lifecycle, adapters and GUI remain separate.

While verification runs, local main advances to `54da62df`, changing Help
visibility/test code and documentation only. The working base and RNA source
remain the above `c5a9d63f` plus this diff; no rebase or later-tree pass is inferred.

#### S2.1b File Boundary, 2026-10-28

After committing S2.1a, the owner authorized S2.1b from the existing reviewed
plan and separately allowed focused locked/offline local Rust tests/checks.
S2.1a was committed as `18116bf7`, then rebased onto local main `54da62df`,
becoming `f78a592b`; both changelog entries were retained. All 22 recorded
Rust/Cargo/config hashes still match the S2.1a execution receipt. This is not
execution of the later integrated tree. No fetch or push; unrelated `output/`
remains untouched.

The S2.1b working diff adds the strict versioned dictionary/request/receipt
and one direct `&self` engine file method below `primer_payload`, beside S1's
`observed` boundary. It reuses the incremental consumer, canonical geometry,
verified-reader and control/writer patterns. The small preparation vocabulary
is re-exported once through protocol aliases; context, limits and parser stay
crate-private. No legacy public DTO, grammar, limit or validation change.

The [protocol](protocol.md#s21b-annotation-preparation-direct-engine-only)
specifies the raw/semantic digest recipes, complete-analysis versus execution
status, receipt-only refusals, conservative JSON/memory accounting and path/
publication limits. Exact raw GTF line framing is preserved, dictionary bytes
are independently bounded and verified, and existing-file identities supplement
the legacy reader's hash/size/mtime checks. The preparation-specific publisher
closes/synchronizes staged files, reserves fresh output and links the report
last with no clobber. Fault cleanup removes only identity-matching owned files
and empty directories; unowned additions survive with an explicit orphan error.
It does not strengthen or refactor the old publisher under cover of this slice.

Inline synthetic regressions cover raw LF/CRLF/final-line geometry and release
bindings, exact configurable IDs, policy-dependent receipts, missing/refused
sources, strict dictionary fields/metadata/competitors, hash/parse/resource
failures, request-bound early/late control, raw paths/hardlinks, Unix symlinks/
input replacement, report-last visibility and owned/orphan cleanup. No new
committed fixture or generated file requires checkout regeneration.
S2.1c named operations, cancellation command, lifecycle/approval/no-save
policies, shell/CLI/MCP/discovery, SAM ingress and GUI remain unimplemented.

Native Windows handles/case/rename behavior, real disk/sync failure, default
GUI/full suite, real-provider/gap policies and scientific/release acceptance
remain separate. Same-filesystem hard-link support is required for the final
completion marker. No hostile concurrent-filesystem or crash-durability
guarantee is claimed.

##### S2.1b Verification

On macOS/aarch64, the final working diff above `f78a592b` passes:

```text
env CARGO_INCREMENTAL=0 cargo test --locked --offline --no-default-features --lib -j1 primer_payload -- --test-threads=1
64 passed (15 new S2.1b file cases; unchanged legacy serialized-row/work oracle)
env CARGO_INCREMENTAL=0 cargo test --locked --offline --no-default-features --lib -j1 observed_transcript_structures -- --test-threads=1
58 passed
env CARGO_INCREMENTAL=0 cargo test --locked --offline --no-default-features --lib -j1 s2a_ -- --test-threads=1
6 passed
env CARGO_INCREMENTAL=0 cargo test --locked --offline --no-default-features --lib -j1 observed_transcripts_ -- --test-threads=1
16 passed
env CARGO_INCREMENTAL=0 cargo check -q --locked --offline --no-default-features --lib -j1
exit 0
```

Filters overlap; do not sum them. No ignored Rust tests in these filters.
The first check found two new helper type-inference errors, corrected before
successful execution. A first 13-case selected run passed; ownership review
then added the stale-temp-name guard and three final file tests. The final
payload rebuild took 3m56s; no timing/performance acceptance is inferred.
Baseline headless dead-code and linker `__eh_frame` warnings remain visible.

The 28-path Rust/Cargo/config manifest at
`/tmp/gentle-s21b-final-source-manifest.sha256` has SHA-256
`9794dd19c7105eec3c51723b2575cc5f6c97ea9d13771e3ed715a7a827a83cdb`;
all paths recheck unchanged. The test binary
`target/debug/build/GENtle/002e669310863d66/out/gentle-002e669310863d66`
has SHA-256
`f30a22a9f3bd57213b5e037f63699349f4970da5327cb1f17c50cdcd5e8d83b1`.
Cargo.lock remains
`40ae5b6c9ba74ba080d39be11b4e57ad2dbc020dba5ecad7222e22f639f84918`.
Toolchain: Cargo 1.100.0-beta.5 `7552a885e`, rustc 1.100.0-beta.5
`8ea38001b9da1d6683ac5ebbc993cac89d3a09df`, LLVM 23.1.3,
`aarch64-apple-darwin`. Effective unoptimized test/dev recipe: root debug 1,
dependency debug 0, 256 codegen units, `CARGO_INCREMENTAL=0`, unchanged
16-MiB configured test thread stack. These are local working-diff receipts,
not committed fixtures, native acceptance or a future merged-SHA claim.

Fast source/documentation checks run 34 cases: 33 pass, one explicitly skipped
for lack of a supplied actual-process payload CLI binary. Documentation status,
scoped rustfmt and whitespace pass. Whole-workspace `cargo fmt --check` still
reports formatting in 12 unrelated baseline files; preserve them untouched
(`/tmp/gentle-s21b-fmt.log`). Session-close reports 4 OK, 2 warnings for the
intentional diff/untouched `output/` and manual plan fidelity, no failures.
The diff stays within the authorized S2.1b contract; no adjacent cleanup,
dependency/profile/fixture changes or later slice is included. Native/default-
GUI/full-suite, real-provider, scientific and release checks remain unrun.

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
| S2.1a grammar/index | Comments, blank records, unterminated final line, field counts, consumed/unknown repeated attributes, quoted semicolons/escapes, invalid UTF-8/controls and bounds; interleaved genes, normal exon repetition, selected-linked collisions in otherwise unselected rows, required/missing/multiple transcript declarations, version scope and excerpt completeness nonclaims. Legacy validators still accept touching exons; preparation alone refuses them. |
| S2.1b dictionary/receipt | Duplicate names/fields, metadata mismatch, positive lengths and checked integer conversion; missing/refused gene produces receipt only, never an empty annotation; exact-ID absence versus matched gene without exons; identical geometry with changed admitted policy gives a different receipt identity; LF/CRLF changes raw hashes/snapshot release but not transcript geometry. |
| S2.1b publication | Partial/cancelled versus fully verified outcomes, early/late control failures, input/output/marker aliases and marker containment, source change/symlink/hash errors, output exists/write/sync failure, no-clobber and report-last visibility. Exact contents are checked without treating directory order as write chronology. |
| Authority | Always snapshot-tagged prepared release plus nullable declared provider release. Same chain key with different verification evidence stays same key/different authority; joins do not upgrade declared to verified. |
| Resources/control | Exact-boundary/one-over cases for every cap, including whole-file index/unselected data/comments; full EOF hash, cancellation/control/data errors, output failure and safe cleanup/orphan disclosure; no final report on data/hash/resource errors, no reuse of uncommitted output, idempotent refusal of existing output and bounded saved inspection. |
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

- S2.0a/b are implemented with scoped results above. S2.1a is subsequently
  owner-authorized with focused local verification; preserve all older guards
  and verify its working diff without inferring S2.1b/c authorization. The
  supplied review alone was not implementation approval.
- Recommend GTF first, genePred later, with explicit real-provider
  `transcript_feature_rows` policy and input scope. Counting unrelated semantic
  inconsistencies and refusing identical selected-model duplicate rows are
  Codex's v1 technical recommendations, not owner-approved provider policy.
- Retain the proposed finite preparation budgets; larger raw-stream/index
  limits for real whole-provider files need separate approval. Do not confuse
  S1's 128,000,000-byte JSON/input limits with the proposal's MiB units or treat
  a budget refusal as gene absence.
- Keep the exact row-to-unit manifest in later S2.2. A compact QNAME-uniqueness
  declaration could be a later reviewed producer contract, but must not silently
  replace current accounting or prove original FASTQ identity.
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

The later S2.1 review reconciliation was source-checked at `c5a9d63f` after the
interrupted read-only session resumed with the same clean tracked baseline.
It changes only this canonical plan, the umbrella plan's pointer, the existing
roadmap entry and changelog. Documentation status, all 18 LF/CRLF-aware
`DocumentationStatusTests`, local Markdown targets and whitespace pass.
Session-close reports 4 OK, 2 warnings (intentional four-document working diff
and untouched `output/`, manual plan fidelity), no failures. This is a scoped
technical recommendation, not owner authorization for S2.1a/b/c or a fresh
Claude verdict; no Rust code, schema, generated artifact, fixture, build/test,
commit, fetch, push, real-data or scientific acceptance is added.
That paragraph records the review-only reconciliation. The subsequent owner
authorization and executed S2.1a working diff are recorded in the implementation
ledger above; they do not authorize S2.1b/c or change that historical receipt.
That paragraph records the earlier review-reconciliation session only; the
separately authorized S2.0a implementation/execution ledger above supersedes its
pre-implementation status, not its source-inspection limits.
