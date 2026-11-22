# Bacterial and eukaryotic transcription termination and 3′-end annotation plan

Date: 2026-11-16. Status: bounded software implementation; historical local verification and later commit-review results are distinguished below. Scientific/native follow-ups remain open. No release commitment.
Repository inspected at `b62d21c971c207704b96396523eec0ce61659806`, with pre-existing local changes.


## Original Codex plan, review result and reconciliation

The original seven-step plan proposed: (1) explicit host/polymerase and search
context, (2) known-reference recognition with exact/bounded approximate matches,
(3) separate bacterial intrinsic, source-backed Rho, Pol II and Pol III modes,
(4) conservative transcript/UTR geometry, (5) portable evidence reports,
(6) shared-engine integration with selected undoable writeback, and
(7) a first milestone combining recognition, bacterial intrinsic candidates and
a mammalian Pol II profile, followed by deterministic and scientific validation.
The complete review prompt, including that original plan and repository evidence,
is saved in [the Claude review prompt](transcription_termination_claude_review_prompt.md).

**Review status:** owner-supplied Claude feedback received on 2026-11-16,
self-reported at `b62d21c9`, read-only with no code execution. The complete supplied
text is retained in [the review record](transcription_termination_claude_review.md).
Codex checked the cited feature-library intake, RNA input guard, cryptic-splicing
wrap refusal, helper-component fields and TATA materialization path against the
working tree. This is a supplied second opinion, not an independently reproduced
Claude execution receipt or a scientific sign-off.

The earlier local CLI attempts remain historical: one silent attempt was stopped
after more than six minutes; the retry failed with `401 OAuth access token is
invalid`. They produced no feedback. The newly supplied review supersedes the
previous unavailable-review status; no new Claude call was made.

**Reconciliation and revised Codex decisions:**

- Accept a firm ownership decision: Phase 1 supplies a small reusable exact
  reference-recognition primitive, with termination as its first consumer. This
  is the exact-recognition foundation for the parked feature-library idea, not
  a second termination-private matcher or approval for a general annotation suite.
- Accept catalog provenance as a prerequisite. No usable termination sequence
  catalog was identified in the review; helper components contain names/attributes,
  not matchable sequence records. Use a dedicated data catalog, not sequences
  embedded in helper semantics or Rust protocol definitions. Reuse the repository's
  redistribution vocabulary. Agents may gather sources and document terms;
  unresolved redistribution and scientific suitability require owner disposition.
- Accept TATA as the primary screen/materialize precedent, including selected IDs,
  digest binding, per-row qualifiers and `prepare_sequence`; do not reuse
  construct-reasoning eligibility as a prediction writeback API.
- Accept non-wrapping initial searches and a separate circular-support slice.
  Qualify “linear only”: explicit non-wrapping intervals on circular plasmids are
  useful, but whole-circle boundary-complete scanning remains unsupported and
  must never be reported as complete. Later wraps need segmented locations and
  topology-aware bounds; restriction-site bounds alone are not a full solution.
- Accept pure heuristic scoring, declared mismatch/wobble policy, and no external
  folding in v1. Existing RNAfold rejects DNA input and lacks a windowed API and
  tool-version receipt. Future folding needs a separately versioned supplement
  outside the core identity, or explicit version-pinned identity and refusal.
- Accept explicit subject binding, unassigned whole-sequence candidates, typed
  unsupported outcomes, uncalibrated Pol II defaults and absent-motif non-claims.
- Accept parity in every executable phase, catalog provenance and scoped checkout
  coverage from its introduction, and existing-terminator deduplication tests.
  “No hits” is a completed result outcome, not a fabricated biological hit row.
- Keep independent biological benchmarking required before accuracy claims.
  Qualify “owner-blocked”: no source/curator is selected yet; identifying sources
  is a prerequisite task, not evidence that all agent progress is blocked.
- Use proposed vocabulary `terminators screen` / `terminators materialize`, with
  `ScreenTranscriptionTermination` / `MaterializeTranscriptionTerminationFeatures`.
  These are proposed names, not implemented commands. Pol III stays deferred.

**Continuation reconciliation (sections 6–11):** accept the concrete protocol
module placement, first-slice parity inventory, generated matrix freshness,
confirmation refusal/state-byte test and catalog-readiness decision. Source
inspection confirms the dedicated TATA records, screen/materialize digest fields
and four command/operation capability registrations. Adopt the TATA contract
pattern narrowly rather than literally copying its biology or entire module.
Do not add approximate-match controls to the exact-only v1 request schema;
reporting exact matched extent is still legitimate geometry, not approximate
matching support. Defer a dedicated GUI honestly in the generated matrix; generic
MCP `op` needs no dedicated tool. Preserve the earlier qualifications about agents
helping collect sources and independent benchmark prerequisites. Accept a typed
catalog contract in the protocol crate; catalog sequence data remain a dedicated
external data document, not Rust constants or helper-component fields. Choose
caller-supplied catalog for v1, with explicit readiness/unavailability semantics;
no curated biological collection is assumed bundled. The proposed guide-table
correction is recorded as follow-up, not included in this scoped plan change.

The revised phases below replace the original combined first milestone. No
implementation or release commitment follows from accepting review suggestions.

## Objective and biological contract

Provide explainable candidate annotations for plasmid expression units and
incompletely annotated transcript 3′ regions, followed by explicit, undoable
application of selected features. A missing 3′ UTR is an uncertainty to describe,
not permission to infer an exact mature transcript end from a DNA motif.

- Subjects: exactly one `seq_id` plus explicit interval, a whole linear sequence,
  or one resolved transcript/CDS feature ID. Reject conflicting selectors and
  ambiguous identities. Whole-sequence hits are unassigned; transcript association
  requires an explicitly named transcript. A named CDS alone is CDS context, not
  a resolved transcript. Circular origin-crossing searches are deferred. Record intended expression host and polymerase for
  each expression unit; plasmid propagation host is not expression context.
- Coordinates: source intervals use zero-based half-open coordinates with explicit
  strand. Preserve ordered segments for spliced transcripts; later circular
  support must also preserve ordered origin-crossing segments;
  display conversions must not alter stored geometry. Downstream means biological
  5′-to-3′ direction, including on the reverse strand.
- Fixed semantics: distinguish a known element match, a predicted motif/structure,
  a cleavage window, an annotated mature RNA end, and a polymerase termination
  region. None automatically establishes another or proves efficiency.
- Adjustable policy: explicit versioned profile, scan bounds/flanks, structural
  thresholds, reference identity/coverage and output/work limits. No universal
  species-independent termination score or probability.
- Ambiguity: retain multiple transcript/end candidates; refuse unresolved required
  identities rather than choose an isoform. Report unknown strand/context and
  unsupported mechanisms as typed mechanism outcomes (including a requested Rho
  de novo analysis). Both-strand non-wrapping plasmid scans return unassigned hits,
  without inferring expression units or associating hits to the nearest CDS.
- Output: portable report with optional selected feature writeback. Existing
  annotations remain separate source evidence. Searching does not alter sequence
  or annotations.
- Provenance: sequence/topology digest, source/transcript/feature identities,
  search policy and actual bounds, profile/reference versions and hashes, evidence
  origin, component scores, exclusions and truncation/completeness status.
- Non-claims: no measured termination efficiency, expression, transcript abundance,
  exact cleavage site or completed transcript from motif presence alone. Literal
  read-end A/T runs are not physical tails or polyadenylation-site proof.

For Pol II, mature RNA cleavage and downstream polymerase termination must stay
separate. Primary mechanistic context:
[DNA-directed termination of RNA polymerase II transcription](https://pmc.ncbi.nlm.nih.gov/articles/PMC7615648/).
Bacterial intrinsic candidates require structure and sequence context, informed by
[Structural basis for intrinsic transcription termination](https://pmc.ncbi.nlm.nih.gov/articles/PMC10091898/).
These sources motivate the contract; they do not calibrate a new GENtle predictor.

## Existing implementation seams

This is a missing biological primitive plus orchestration, not an existing
parameter adjustment. Existing feature roles and overlap annotations do not supply
termination prediction.

- `crates/gentle-protocol/src/tata_boxes.rs` and
  `src/engine/analysis/tata_boxes.rs`: primary precedent for schema/non-claim
  constants, typed evidence, `content_sha256`, `expected_report_sha256`, selected
  `row_ids`, idempotence qualifiers and `prepare_sequence` after feature changes.
  Reuse conventions narrowly, not TATA biology or unchecked wholesale copies.
- `src/engine.rs`: operation registration, read-only/history classification and
  undo policy. Construct roles describe existing feature conventions only;
  generated-annotation/SequenceSpan eligibility is graph-specific, not this API.
- New records belong in `crates/gentle-protocol/src/transcription_termination.rs`,
  registered with a documented `pub mod` in that crate's `lib.rs`. Follow its
  operation-name registry pattern. `src/engine/protocol.rs` needs only the report
  persistence field using the dedicated protocol types, not new request/report
  definitions. The biological extension guide's central-file table is older than
  this extraction; its broader correction is separate follow-up work.
- `src/engine/analysis/feature_expert_ops.rs`: transcript groups/exon geometry;
  `src/engine/ops/operation_handlers.rs`: transcript-local projection and shared
  dispatch. Reuse these interpretations rather than reconstruct in a frontend.
- `src/genomes.rs`: `HelperConstructComponent` has no sequence field; its
  terminator classification is a naming layer, not a matchable catalog. Reuse
  `redistribution_status` terminology, not unrelated helper-container ownership. Existing promoter terminator-overlap evidence stays intact.
- `src/engine_shell.rs` and its parsers: shared commands, discovery and readiness;
  GUI, MCP, JS/Lua and Python project the same operation/report.

## Revised implementation sequence

### Phase 0: freeze the bounded contract

Freeze report/request fields, feature types/qualifiers, supported host profile,
non-wrapping interval rules and scan budgets. Exact recognition is one minimal
shared primitive serving the termination consumer and future feature-library work
(`docs/legacy_gentle_m_intake.md`); no parallel private matcher. The dedicated
versioned data catalog contains sequence, element identity, host/mechanism,
source/version/digest, retrieval/recreation instructions, redistribution status
and approval provenance. Its home is a scoped data directory with a nearby
manifest/README, finalized with the first implementation patch. Protocol records
validate catalog data; helper semantics may link identities without owning bases.

V1 catalog delivery is caller-supplied: the request identifies one catalog file
and expected content binding. No biological catalog is silently bundled or
fetched. Readiness advertises the catalog requirement and selected-source
availability when known; execution validates readability, schema, provenance and
hash binding. Missing catalog is `catalog_unavailable`, malformed input is an
input error, and changed content is a freshness refusal; none is a no-hit result.
A future bundled catalog can change readiness only with documented delivery,
source and redistribution decisions. Synthetic fixtures are not a default
biological catalog.

Schema conventions: named schema and non-claim constants; screen request with
`#[serde(default, deny_unknown_fields)]` and explicit visible `Default` policy;
report `content_sha256` computed with that field empty; materialize request with
`screen`, `expected_report_sha256` and `row_ids`, rejecting unknown fields. Put
catalog selection in the screen request so recomputation binds the same inputs.
Do not add identity-threshold, approximate coverage or indel knobs to v1.

Author hand-crafted non-biological test catalogs first. A shipped biological
catalog requires selected sources and documented terms; unresolved clearance is
an explicit owner decision. This phase adds no predictor or claim of usable
curated biological coverage. Reject unknown fields and test schema round trips.
Enumerate wrap refusal, ambiguous transcript, unsupported mechanism and missing
required context/flanks as structured outcomes. Invalid requests remain errors;
valid completed searches with zero hits remain distinguishable.

Acceptance: schema examples distinguish signal/cleavage/termination geometry,
missing-CDS versus incomplete-transcript cases, unsupported context, zero hits,
insufficient flanks and truncated searches. Profile thresholds require explicit
rationale and no calibrated-accuracy claim before benchmarking.

### Phase 1: known-element recognition and report-first application

Recognize a small explicitly versioned, redistribution-cleared set of termination
or cleavage/polyadenylation elements. Begin with exact sequence matching;
approximate matching is a separate extension with explicit identity, coverage,
indel and ambiguity policies. A complete known element and its short motif are
separate evidence objects. References retain documented host and mechanism.

Create one engine-owned screen report and a separate materialize action. Require
selected row IDs and the reviewed report/input binding, reject stale inputs, skip
existing source annotations, prevent duplicate features and support undo/redo.
Feature kinds must preserve distinctions such as `terminator`, `polyA_signal`,
`polyA_site` and candidate `misc_feature` as supported by actual format semantics;
never export an uncertain cleavage window as an exact observed polyA site.

Acceptance: deterministic forward/reverse matches within non-wrapping spans,
ambiguous-base policy, no-hit outcomes, provenance export, stale-report refusal,
idempotence, existing-annotation preservation and undo/redo. Include an existing real feature with
`regulatory_class=terminator` and matching geometry: retain it as source evidence
without a duplicate generated annotation. Define overlap/strand equivalence
explicitly; a nearby or opposite-strand feature is not automatically a duplicate.

Reject origin-wrapping requests with a structured unsupported-span outcome.
Non-wrapping intervals on circular plasmids must report boundary limitations;
whole-circle searches cannot claim complete coverage without wrap support.
Add registry, read-only/mutating classification, protocol docs, glossary,
capability/readiness and MCP/generic-op reachability plus shell/CLI/JS/Lua/Python
parity in this executable slice. GUI can use the existing shared operation surface
before dedicated refinement; no new frontend biology. Bind saved report identity
and catalog identity. This surface work is part of the phase estimate, not polish.

#### Phase 1 parity inventory and acceptance gate

Budget the following with the matcher and materializer, not as Phase 4 polish:

| Surface | Planned change or check |
| --- | --- |
| Protocol | Dedicated `transcription_termination.rs`, documented `pub mod`, operation-name registry in `crates/gentle-protocol/src/lib.rs` |
| Engine | Two variants and correct read-only/history classification in `src/engine.rs`; shared dispatcher in `src/engine/ops/operation_handlers.rs` |
| Analysis | One termination analysis module; minimal reusable exact-recognition helper; no construct-graph writeback dependency |
| Persistence | Typed optional report field in `src/engine/protocol.rs`, backward-compatible project deserialization and round-trip tests |
| Discovery | Four descriptors in `src/engine_shell.rs`: each shell command and operation name; consistent `mutating`/`requires_confirmation` and readiness facts |
| Parser/glossary | Working `terminators screen` and `terminators materialize` routes, help/error text and `docs/glossary.json` entries |
| MCP/adapters | Generic typed `op`; confirmation enforced for materialization; CLI, JS/Lua and Python shared-route reachability |
| Generated matrix | Regenerate `docs/gui_cli_mcp_parity.md` using `scripts/regenerate_parity_matrix.sh`; check actual generated command/operation rows and freshness |
| User contract | Protocol/CLI documentation and accurate GUI reachability description; dedicated GUI form deferred |

Do not add `docs/parity_matrix_overrides.json` entries merely to hide missing
reachability. Let the generated matrix record actual prominence, including
`shell-only` where appropriate. A dedicated GUI module, palette entry or subject
selection integration is not a Phase 1 gate.

Required deterministic checks include:

- MCP materialize without `confirm` returns an error and leaves the state file
  byte-identical; a confirmed valid request applies the reviewed rows.
- Screen's read-only classification and any explicit export/persistence route
  retain the existing write-confirmation contract; no hidden write in inspection.
- Generic glossary parser smoke tests
  `glossary_cli_usage_smoke_commands_parse` and
  `glossary_cli_usage_flags_parse_one_by_one` pass with the new entries.
- Generated parity matrix freshness passes via the existing generator check.
  Inspect existing `.gitattributes` and fast `scripts.test_tutorial_checkouts`
  coverage and extend only the missing scoped LF/CRLF cases; do not introduce a
  redundant freshness mechanism or regenerate unrelated projections.
- Catalog absent, malformed, changed and valid-zero-hit outcomes remain distinct.
- Existing kind `terminator` plus `regulatory_class=terminator` survives source
  recognition and reapplication without duplication, including export/reimport
  through the supported annotated format.
- Core identity remains independent of RNAfold installation. Test with controlled
  tool resolution rather than requiring ViennaRNA on CI or invoking it.

### Phase 2: bacterial intrinsic candidate analysis

Add a bounded, versioned structural heuristic for RNA-oriented hairpin plus
U-rich tract, with explicit stem/loop, mismatches, spacing and tract components.
Declare a versioned mismatch/wobble allowance as a heuristic policy, not proof
of RNA folding. Do not present a score as thermodynamic free energy or probability.
Keep the v1 detector pure engine code with no ViennaRNA subprocess or installation
dependency. Core report digests must be identical with and without RNAfold
installed; test this boundary. A future external-folding supplement must not
silently change the reviewed core identity or materialization behavior.

Rho-dependent termination remains source/reference evidence only. Absence of an
intrinsic candidate does not establish absent termination. Archaeal transcription
is explicitly outside the bacterial profile rather than silently included under
"prokaryotic".

Acceptance: synthetic positives and motif/structure negatives, both strands,
origin-wrap refusal, imperfect/ambiguous stems, bounds/work caps, stable ordering
and adapter parity. Independent biological benchmarks need a named source,
curation/rights review and held-out positives/negatives before sensitivity or
specificity claims. No source is selected yet; synthetic software acceptance can
proceed while scientific performance acceptance remains pending.

### Phase 3: narrowly scoped mammalian Pol II 3′-processing candidates

Implement one explicitly named host profile with a documented polyadenylation
signal/variant set, surrounding sequence context, spacing and candidate cleavage
windows. Report alternatives without collapsing motif occurrence, cleavage and
polymerase termination into one coordinate. Do not generalize this profile to
yeast, plants, all eukaryotes or non-polyadenylated transcript classes. Thresholds are authored defaults with no
measured sensitivity/specificity until independently benchmarked. Signal/context
matches do not establish actual cleavage-site usage; expression/3′-end evidence
is a separate requirement, and absence of a canonical hexamer does not imply
absent 3′ processing.

A CDS stop alone supports a downstream search region, not an exact UTR. A partial
transcript supplies only the bounds it actually contains. Known mature endpoints
and spliced transcript geometry can delimit an annotated UTR; predicted windows
remain candidates. Downstream-of-transcript flanks can be inspected only when
available in the supplied DNA and explicitly requested. Do not silently scan a
circular plasmid indefinitely or assign a hit to the nearest CDS.

Acceptance: motif-only negatives, multiple candidate signals, reverse/spliced
transcripts, uncertain CDS/stop/end annotations, missing flanks, explicit wrap
refusal, adapter parity and correct distinction between an annotated UTR and a candidate 3′ region.
No exact cleavage or polymerase-release claim from sequence-only screening.

### Phase 4: presentation refinement and separately scoped extensions

Use one evidence/geometry model for report tables, DNA-map overlays and exports.
Show host/profile, evidence type, uncertainty and alternatives in the annotation
preview. All interfaces invoke the same engine operation; dedicated forms are
optional, reachability through the shared contract is required.

Circular support is a separate deliverable: inspect
`src/restriction_enzyme.rs` topology-aware bounds, reuse suitable geometry and add
segmented feature locations, origin-crossing both-strand tests and no-hit boundary
completeness checks. Cryptic-splicing currently refuses wraps and TATA does not
supply this contract; neither is evidence that origin handling already works.

Later work: calibrated additional host profiles, bounded approximate reference
matching, Pol III oligo-T rules, imported cleavage/termination experiments and
source-bound read-end support. Pol I, archaeal, organellar and specialized
non-polyadenylated mechanisms require their own contracts. Read-end pileup
integration must preserve existing saved evidence and original-read non-claims.

## Verification and delivery requirements

Implementation requires targeted deterministic tests, adapter parser/executor
and capability parity, persistence/export/undo checks, and `cargo check -q`.
Cover single-base/empty boundaries, both strands, topology and origin splits,
ambiguous bases, truncated searches and resource limits, duplicate annotations,
stale reports, alternative endpoints and source freshness. Document all committed
fixtures with origin, deterministic recreation and runtime/test usage. Keep
training/reference material distinct from independent predictive evaluation;
software test success is not scientific validation.

From Phase 1, a hashed catalog is byte-bound input: define scoped checkout
policy and LF/CRLF tests in the fast `scripts.test_tutorial_checkouts` gate for
its fixtures/projections as applicable. Preserve digest semantics; do not weaken
hashes to accommodate checkout conversion. Follow `docs/testing.md` portability
rules for every added file boundary. Update protocol, CLI/GUI,
discovery/glossary, tutorial prose, roadmap and changelog when executable behavior
is implemented. Record native Windows/Linux, GUI and real-biological acceptance
separately; do not inherit old-source receipts.

The original planning-only session authorized no implementation. The owner
subsequently requested implementation of as much of this plan as possible.
That authorizes the scoped software work below, without provider-data processing,
publication, release-scope changes or scientific sign-off.

## Implementation ledger — 2026-11-16

Working tree based on `b62d21c971c207704b96396523eec0ce61659806`; unrelated
pre-existing edits and untracked work are preserved. This is new-source local
verification, not a receipt for a committed release candidate.

- Phase 0/1 sources: strict dedicated protocol module, caller-supplied byte-bound
  catalog, shared exact matcher, separate engine screen/materialize, report
  identity, explicit context/strand, non-wrapping bounds, source annotations,
  provenance retention, atomic selection and duplicate prevention. Later selected
  support merges into generated annotations with undo; original source annotations
  remain unchanged. Catalog inputs require regular files.
- Phase 2/3 sources: opt-in pure bacterial intrinsic and mammalian Pol II
  signal/context heuristics with authored versioned policies, component evidence,
  separate cleavage-window geometry (explicit interbase-boundary sets with
  correct reverse projection), resource bounds and missing-context outcomes.
  These are uncalibrated software capabilities, not validated predictors.
- Parity sources: direct CLI forwarding, shared shell/parser, four descriptors,
  glossary, generic MCP op, generic JS/Lua/Python operation routes and existing GUI
  shell/workflow access. Screen migrates to content-preserving lifecycle behavior;
  explicit export attempts retain confirmation. Generic MCP op always retains
  its existing confirm gate. Dedicated GUI form remains deferred.
- Inline synthetic fixtures document their origin/recreation/use in Rust test
  modules. No reference sequence catalog or biological fixture was newly bundled;
  the existing generated matrix LF/CRLF contract is reused.

Local verification (macOS arm64, `rustc 1.102.0-nightly (ac8310ba6
2026-11-09)`, working tree above):

| Check | Result |
| --- | --- |
| `cargo test -q --no-default-features --lib termination_` | 21 passed, including prediction geometry, provenance/GenBank, undo and shell/MCP/lifecycle boundaries |
| `cargo test -q --no-default-features --lib engine::exact_reference` | 2 passed |
| `cargo test -q --no-default-features --lib glossary_cli_usage` | 2 passed |
| `cargo test -q --no-default-features --bin gentle_cli test_termination_cli_forwarding` | 1 passed, preserving JSON, Windows-style paths, apostrophes and output paths |
| `cargo test -q -p gentle-protocol transcription_termination --lib` | 1 passed with strict schema/default tests |
| `cargo test -q --no-default-features --features js-interface --lib termination_js_adapter_matches_engine_without_state_changes` | 1 passed through the registered operation implementation |
| `cargo check -q` | Passed with default desktop features |
| `cargo test -q --no-default-features --features lua-interface --lib termination_lua_adapter_matches_engine_without_state_changes` | 1 passed through actual Lua VM; typed project content and report identity retained |
| `python3 scripts/test_transcription_termination.py --cli target/debug/gentle_cli` | Passed with rebuilt headless CLI: Python/direct parity, export, unchanged project bytes, apply/idempotence, GenBank reimport and raw LF/CRLF freshness |
| `target/debug/gentle_examples_docs parity-matrix-generate` then `parity-matrix-check` | Generated four new command/operation rows; freshness passed, no overrides added |
| `python3 -m unittest scripts.test_tutorial_checkouts` | 99 passed after generation, including real matrix LF/CRLF and negative-control coverage |
| `python3 scripts/manual_docs.py`; `python3 scripts/generate_capability_maturity.py --check`; `git diff --check` | Passed |

All checks in the table passed on the final source. The first JavaScript
attempt failed before GENtle tests because the restricted sandbox could not
resolve/download the V8 archive; the permitted retry passed. Earlier focused
failures exposed missing content preservation, Lua nil-vs-JSON spelling in the
test, and GenBank qualifier folding. The fixes retain the real boundaries rather
than weaken them. Existing headless dead-code and large-test-binary compact-unwind
linker warnings remain; no test failure is waived.

Session-close hygiene passes with two reminders: pre-existing dirty/untracked
work remains, and plan fidelity requires manual review. No full suite, native
Windows/Linux, dedicated GUI smoke or scientific validation is claimed.

Files changed for implementation:

- New contracts/analysis: `crates/gentle-protocol/src/transcription_termination.rs`,
  its `lib.rs` registration, `src/engine/analysis/exact_reference.rs`,
  `termination_predictions.rs` and `transcription_termination.rs`.
- Engine glue: `src/engine.rs`, `src/engine/protocol.rs`,
  `src/engine/operation_lifecycle.rs`, `src/engine/ops/operation_handlers.rs`.
- Adapter/discovery/tests: `src/engine_shell.rs`, `src/engine_shell/tests.rs`,
  `src/bin/gentle_cli.rs`, `src/mcp_server.rs`, `src/js_interface.rs`,
  `src/lua_interface.rs`, `scripts/test_transcription_termination.py`.
- Required result-field initialization only: `src/app.rs`, `src/app/tests.rs`,
  `src/main_area_dna/tests.rs`, `src/engine/analysis/rna_reads.rs`.
- Documentation: this ledger, `docs/protocol.md`, `docs/cli/features.md`,
  `docs/decisions.md`, `docs/CHANGELOG.md`, `docs/roadmap.md`,
  `docs/glossary.json` and generated `docs/gui_cli_mcp_parity.md`.
- Required glossary source-hash projections only: `assets/capability_maturity.json`
  and `integrations/clawbio/skills/{gentle-cloning,gentle-pcr-primer-design,gentle-tss-collection}/CAPABILITY_MATURITY.json`.
  No maturity status, dependency/lockfile or unrelated source change is included.

Remaining scope: independent scientific calibration and a redistribution-reviewed
biological catalog; origin-crossing/circular-complete analysis; dedicated GUI
presentation; approximate matching, additional hosts, Pol III and read-end
integration. None is implied by sequence candidate detection or existing source
annotation. A broader guide-table correction remains separate follow-up.

## Commit Review, 2026-11-20

Review began on `main` at `ccd839dc` with this implementation still uncommitted.
The preceding macOS test table belongs to the earlier implementation session;
it is retained as reported working-tree evidence, without an immutable candidate
or binary receipt and without claiming that Rust was rerun during this review.

One concrete discovery defect was reproduced with
`scripts.test_capability_maturity.DiscoveryCatalogTests.test_mcp_op_glossary_commands_are_advertised_by_the_operation_tool`:
`terminators screen` and `terminators materialize` were absent from the generic
MCP operation tool's `commandPaths`. Add only those two mappings and assert their
unique presence in the existing MCP inspection/confirmation test. The Python
guard then passes; the extended Rust test remains unexecuted here. Scoring,
catalog validation, report identity and writeback semantics are unchanged.

After the repair, 128 Python source/checkout tests pass, including LF/CRLF and
discovery coverage. Documentation status, manual/maturity freshness, new-module
rustfmt, MCP source parsing and whitespace checks pass. Session-close reports
four OK, two expected dirty-work/scope reminders and no failures. No additional
product defect was established by this source review; no Rust/native acceptance
is inferred from the Python results.

The [new-source handoff](glen_next_verification_handoff.md#termination-evidence-commit-review-2026-11-20)
requires compile, actual adapter/process and canonical parity checks after the
owner's push. No local Rust build/test, provider call, push or publication is
authorized or performed by this commit review.

## Updated-Lock Process And Native Follow-Up, 2026-11-22

With unchanged production source and the already built binaries described in
the [source/lock receipt](glen_next_verification_handoff.md#updated-lock-macos-pcr-verification-2026-11-22),
the actual Python/direct CLI companion passes. The existing default-feature
test executable also passes the selected MCP discovery/confirmation regression
(1 passed, 5178 filtered); this does not execute an external MCP server or JS/Lua
VMs. No further Rust build occurred.

Native macOS GUI Shell confirms exact-reference report parity, stale-digest
refusal, explicitly selected single-feature materialization and feature removal
by Undo. It also exposes a concrete history defect: Undo reports one redo entry,
but the immediate Redo fails. The shell automatically persists its changed UI
command text after success, reaching the engine's intentional side-edit redo
invalidation. Do not weaken that engine safeguard or claim native history green.
The GUI persistence boundary needs a separate repair and regressions covering
immediate Redo and read-only inspection between Undo/Redo. No repair is included
in this evidence-only update.

The [current handoff](glen_next_verification_handoff.md#termination-process-and-native-shell-return-2026-11-22)
records exact binary identities, unchanged synthetic project bytes and remaining
VM/stdio/parity/platform scope. No private inputs, scientific acceptance, push
or publication is involved; Windows remains untested natively.

## GUI Shell History Repair, 2026-11-22

The native follow-up above identified an adapter persistence defect, not a
termination-kernel or engine-checkpoint defect. Remove the unconditional UI
settings save after successful shared-shell execution. Preserve the command
result, display synchronization, legacy saved UI record and the engine's
intentional Redo invalidation for genuine persisted side edits. This is not a
special case for termination or the history command family.

Five synthetic controller regressions cover direct and alias Undo/Redo,
intervening inspection, failed commands, real new operations and explicit UI
side edits. They exercise `run_shell_command`, not just a new helper. No local
Rust build/test or fresh native acceptance is claimed; old binaries still
contain the reproduced failure. Follow the [focused compile/test/native request](glen_next_verification_handoff.md#gui-shell-redo-repair-2026-11-22).
The owner keeps native Windows separate as an external gate. Existing
termination process/MCP-handler passes remain source-bound and are not reissued.
