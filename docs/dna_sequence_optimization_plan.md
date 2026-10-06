# Bounded synthetic sequence design

Updated 2026-10-06, implementation started from `fc4458c3`.
Rebased during implementation onto fetched `main` at
`4a71d5bd15f6b0da60709ee79a431668756d110f`, then onto local `main` at
`564ea09be48149649748d620f57973d5f87cfff7` after its dependency/Help update.
The latter was ahead of the fetched remote at integration time.
Owner authorized work from the previously Claude-reviewed proposal. This is
experimental development outside the `.12` release gate, not release acceptance,
full DNA Chisel compatibility, or a measured performance improvement.
Owner-requested follow-up on 2026-10-05 adds the explicit conflict-search slice
below, based on integrated local `main` `d221aa1a`. Glen's later first-slice
audit reported `needs_fix` at `85d23ddb`; independent repair/new-slice acceptance
remains pending. Local tests do not substitute for his usefulness, performance
or independent correctness verdict.
The earlier Claude review covers the original bounded proposal; no additional
Claude review of this conflict-search implementation is claimed.
The owner approved the CDS-wide GC slice after Claude's read-only review of
`dc8c5252`, received 2026-10-06. Implementation rebased early onto fetched
`main` `2f3a60f0bbf9a9ef61a02fbae6978305f5ac7df3`, preserving its separation of
output validation from submitted search claims and its single apply/undo boundary.
The GC proposal, review and reconciliation are recorded below; no runtime,
native-platform, Glen or release acceptance follows from this authorization.
The owner authorized the next bounded slice on 2026-10-06: explicit sliding-window
GC through the same shared contract. Implementation starts from fetched `main`
`afbef44451ed61e5f347c3aeb586e7124201cf22`, retaining its typed unverified-preview
wrapper and single materialization boundary. An optional read-only Claude review
was offered, not performed; this windowed slice is not Claude-reviewed.

Review rationale is preserved in the [dated review history](dna_sequence_optimization_review_20261005.md)
and [original consultation prompt](dna_sequence_optimization_claude_prompt.md).
Those historical proposal stages do not override this owner-authorized contract.

## Adopted first slice

Remove explicitly specified finite IUPAC motifs from a synthetic coding insert
using synonymous substitutions. Preserve the supplied protein sequence, literal
ATG, literal terminal stop, flanks and protected bases. Changed nucleotides are
the secondary objective; do not infer minimum edits from incomplete search.
Natural PATZ1 assay templates, reverse translation and existing primer design
are unchanged. A label never implies permission to recode a sequence.

The first reviewed proposal reduced four constraint families to motif removal.
Translation/protection are mutation invariants, not optional objectives. The
later owner-approved GC slices add CDS-wide and explicit sliding-window constraints
below; codon adaptation, repeats/hairpins, folding and a general specification/plugin
framework remain deferred. The independent component option is adopted as
`crates/gentle-sequence-design`, version 0.1.0, private/unpublished, original MIT
code, no dependencies. It supports Rust 1.85 / edition 2024. No Python runtime,
new dynamic ABI, GUI dependency or published separate repository is introduced.
The core does not receive the root codon asset's positional 64-character string:
the root resolves every explicit codon/residue pair through the existing table
and binds the actual complete mapping and identity in the report.

Upstream [DNA Chisel](https://github.com/Edinburgh-Genome-Foundry/DnaChisel)
offers broader Python specifications and a GenBank-annotation CLI. GENtle's
initial parity is across its own shared operation, shell/CLI and inner agent;
its request format does not claim compatibility with DNA Chisel's annotation
language or algorithm. A future pinned comparison may use upstream CLI as
optional development tooling, comparing semantics, not identical chosen DNA.

## Admission and coordinates

- Required request schema: `gentle.dna_sequence_design_request.v1`.
- Explicit `purpose=synthetic_coding_insert`; exact inline DNA or loaded sequence.
- Uppercase A/C/G/T only; linear DNA, one explicit forward contiguous CDS,
  length divisible by three, complete ATG-to-stop, no internal stops.
- Required protein includes initial M but excludes stop; it must match input.
- Standard genetic code 1 only. The frozen ATG is not an initiation-efficiency
  model or proof of an authentic CDS. No inferred translation exceptions.
- All intervals are local zero-based half-open. Flanks and the start/stop triplets
  are frozen; intersect synonymous choices with every protected base.
- Patterns are uppercase finite IUPAC strings, not regexes/enzyme aliases.
  Explicit strand is `forward`, `reverse` or `both`. Reverse scanning uses the
  reverse-complement pattern on the unchanged supplied sequence.
- All full-sequence overlapping windows and CDS/flank boundary matches count.
  A palindromic mask pattern requested on both strands yields one `both` match.
  Other patterns can yield separate forward/reverse matches at the same span.
- Ambiguous DNA, incomplete/reverse/joined CDS requests, other codes, indels and
  circular topology are unsupported or invalid, never silently converted.
  Source annotations do not determine CDS geometry or synthetic intent.

Numeric limits, checked before search allocation: 12,000 bases, 16 motifs,
32 bases/motif, 256 protected intervals, 1..100,000 requested candidate
evaluations, 4,096 recorded matches. A conservative 50,000,000 full-evaluation
work-unit ceiling reduces the effective budget when necessary; reports disclose
both requested and effective budgets. Work includes full DNA traversal plus all
possible motif/base/strand checks. Excessive matches are a typed invalid resource
outcome, not truncated evidence. Search-space counting uses checked `u64`;
overflow is unknown, never infeasible.
At the engine approval boundary, loaded source IDs are capped at 4,096 bytes
and digest fields at 128 bytes before preview copying/hashing. Oversized inline
targets and mismatched request schema/purpose are also rejected before that work;
an attacker-supplied recomputed approval hash does not bypass these limits.

## Algorithm and outcomes

`synonymous_full_enumeration_v1` is a correctness-first prototype, deliberately
not an efficient general large-insert solver. Original codons precede synonyms
ordered by nucleotide edit count then lexicographic DNA. Codons intersecting
initial violations are enumerated first, then remaining codons by local position;
the first codon is the least significant mixed-radix digit. Feasible outputs
are compared by changed nucleotide count, then lexicographic whole DNA. There is
no random seed, cache, localized evaluator or hidden source mutation.

Each candidate receives complete translation/protection/motif evaluation.
Any chosen output receives a fresh full validation. Complete tiny-space
enumeration establishes feasibility/optimality and is checked against explicit
four-variant fixtures. A coupled-edit fixture has one violation before and after
either single edit, and none only after the two edits: it falsifies strict greedy
repair rather than merely containing a two-codon motif.

### Explicit Conflict-Search Slice

The optional `search_strategy` request field selects `conflict_directed`,
recorded as `synonymous_conflict_search_v1`. Omitted/default `full_enumeration`
retains the original implementation, output ordering and serialized approval
basis: the default field is omitted on write, including in old-report replay.
Unknown modes are invalid; no automatic strategy selection is introduced.

Conflict search starts from original codons. At each fully evaluated violating
candidate, take the first match in motif-input/local-position/strand order.
Try every nonoriginal legal synonym of every unassigned codon intersecting that
window, in local codon order and existing edit-count/lexicographic synonym order.
Each branch fixes that codon for its descendants. Backtracking restores it to
the original. The explicit depth-first stack is bounded by admitted CDS length;
it stores branch indices, not whole sequences, and uses no recursion or visited
sequence cache. Candidate and motif-work budgets remain unchanged. Counts mean
fully evaluated tree nodes, not necessarily distinct full-space variants.

Completeness follows from a feasible target needing a different unassigned
codon somewhere in every current violated window: otherwise that match persists.
The branches include its synonymous assignment. Follow that path until a
feasible candidate is reached; it changes a subset of the target's codons.
Thus a fully resolved tree includes a minimum-edit feasible candidate, including
the lexicographically preferred equal-cost output. Once a feasible best exists,
an infeasible node with at least that many edits may be pruned: descendants fix
additional nonoriginal codons, so they cannot attain equal or lower edit cost.
No assumed motif-local validity replaces the full validator, including for new
conflicts outside the original match and across CDS/flank boundaries.

Complete-tree proofs differ from complete full-space enumeration. A budget stop
still means unresolved or feasible/incomplete, never infeasible/optimal.
Search-space count overflow is separately unknown and is not a proof either.
Cancellation discards candidates and proof flags, including after feasibility.
Apply binds the explicit strategy, algorithm and its non-claims, and freshly
validates the exact approved output without executing either search.

Independent literal four-state and eight-state oracles check both solvers,
including degenerate motifs, requested strands, protected bases and tie breaks.
A CTG/GAA/TTC synthetic case requires repairing a newly introduced motif outside
the original site. A repeated-GCT insert demonstrates fewer visited candidates
and completed proof despite an overflowed full search-space count; this is not
a timing, large-insert suitability or general performance verdict. No GC,
adaptation, folding, Python dependency or new GUI editor is part of this slice.

### CDS-Wide GC Slice (2026-10-06)

Original Codex proposal: optional explicit minimum/maximum GC on the declared
CDS, jointly enforced with motifs by both existing strategies, exact approval,
independent exhaustive examples and shared CLI/GUI Shell/inner-agent guidance.
Claude's read-only critique correctly requires a shared feasibility predicate,
GC conflict branches, fresh apply-time GC validation and explicit work accounting.
Reconciliation: extend the existing decision record rather than creating another;
21 bp with 4762..5238 basis points is an empty count window, not an endpoint
acceptance test; check parity freshness without gratuitous regeneration.

Request `gc_content={min_basis_points,max_basis_points}` is optional, inclusive
0..10000, applies to exactly the complete supplied CDS, and includes frozen
ATG/stop and protected bases. Flanks are excluded. G+C count is invariant under
reverse complement; no strand knob, organism default, floating-point tolerance,
arbitrary interval or sliding-window policy is inferred. One explicit GC
constraint permits `avoid_motifs=[]`; without GC the old nonempty-motif rule stays.

For fixed length L, prepare converts once to minimum `ceil(min*L/10000)` and
maximum `floor(max*L/10000)`. Bad bounds are invalid; an empty integer interval
is proven infeasible, not malformed. Sum admissible per-codon minimum/maximum
GC contributions, retaining fixed contributions, to bound possible GC. Disjoint
intervals prove GC impossibility before candidate evaluation, after the first
cancellation poll. Intersection is only a relaxation and proves no feasibility
when integer gaps or motif conflicts remain. Full-evaluation work admission adds
the CDS GC pass; candidate/progress/cancellation limits stay bounded.

One full validator produces motif and GC violations. Enumeration, conflict
search, final output checks and `validate_output` use its satisfaction predicate.
Conflict search keeps motif-first deterministic selection. For GC below/above
bounds with no motif conflict, select all unassigned codons capable of changing
their GC contribution in that direction, then branch on **all** their nonoriginal
synonyms, including neutral/worsening assignments. No per-edit GC monotonicity
is required. A feasible target consistent with assigned codons must net-change
GC in the needed direction, so at least one differing unassigned codon is
direction-capable and its target assignment is included. That preserves a path
to every minimum-edit target; motif branching's argument and edit-cost pruning
remain valid. Unrestricted literal eight-state and 36-state oracles check this
restriction, including coupled motif/GC changes, mixed-direction synonyms and
budget-limited proofs. Each depth-first frame holds only a conflict and branch
cursor, not a copied whole-CDS branch list; stack storage grows with assigned
codon depth rather than the square of that depth. This is a resource bound, not
a runtime performance verdict.

GC-enabled algorithm IDs are `synonymous_full_enumeration_gc_v1` and
`synonymous_conflict_search_gc_v1`, with a GC-specific nonclaim. Absence retains
old algorithm/policy and serialized report/request approval bytes. The optional
report GC block carries integer denominator/bounds, relaxed extrema and input
measurement. Output measurement derives through the pure core from actual
validated output DNA; no output DNA means no output facts. Unrequested or
unadmitted GC facts are absent, never zeros. Apply freshly checks all constraints
and rederives all GC facts, rejecting inconsistencies even in client-rehashed
reports. Its receipt verifies output constraints, **not** submitted search history
or minimum-edit claims. No optimizer is rerun and no source gene is recoded.

The [synthetic MFK walkthrough](cli.md#cds-wide-gc-walkthrough) uses the same
typed JSON and preview/apply commands in CLI, GUI Shell and inner-agent guidance;
it is not a new tutorial catalog or GUI setup path. No expression, folding,
synthesis, large-insert performance or upstream DNA Chisel compatibility claim.

| Outcome | Meaning |
| --- | --- |
| `invalid` / `unsupported` | Admission/resource rule failed; no applicable preview. |
| `feasible` | Fully validated exact output; only complete enumeration/conflict search or unchanged feasible input proves minimum edits. A budget stop can retain a feasible candidate with optimization incomplete. |
| `search_exhausted` | Budget ended without a feasible candidate; unresolved, not infeasible. |
| `proven_infeasible` | An initial match lies wholly in frozen bases, GC count bounds are empty/unreachable, or complete enumeration/conflict search found no admissible variant. |
| `cancelled` | External cancellation; diagnostics only, never output DNA/approval, even if a feasible candidate was already found. |

Core cancellation polls before each candidate and before result publication.
The engine exposes typed sequence-design progress, with candidate counts rather
than a wall-clock estimate, through the existing shared callback. Timing of an
external cancellation is not claimed deterministic.

## Shared engine and approval

`PlanDnaSequenceDesign` / `sequence-design plan` prepares a read-only portable
`gentle.dna_sequence_design_report.v1`. It binds exact source/output DNA hashes,
complete codon mapping/hash, request/protection/motif policies, algorithm,
budget/outcome, matches and nucleotide edits. Preferred/alternative codon counts
from reverse-translation are not reused or relabeled as nucleotide edit counts.
Only `feasible` reports carry `approval_digest`.

`ApplyDnaSequenceDesign` / `sequence-design apply --approve DIGEST` is a separate
explicit mutation. It rechecks digest, source bytes, source annotations/topology,
current mapping and policies, and freshly validates approved output/edits.
It never reruns optimization. All fallible validation, materialization and
metadata serialization finish before state insertion. One outer apply owns the
journal/undo entry, lineage node and singleton container named
`Synthetic sequence design`; their operation IDs agree, without nested public
creation hooks or an inner detached commit. Inline input has `ImportedSynthetic`
provenance, while loaded-source input has `Derived` provenance and a parent edge.
Host detached execution still checks live project/structure/journal identity
before publishing. The retained `gentle.dna_sequence_design_receipt.v1` does not
change the original DNA.

The approval digest identifies reviewed content, not authenticated search
provenance. A caller can change a valid output and recompute that digest, so
apply must not certify the preview's completeness/minimum, counts or reason.
The engine receipt records `output_constraints_verified=true` and
`search_claims_verified=false`. The exact approved preview is retained only as
`submitted_proposal.report`, with the mandatory typed
`verification="unverified_portable_preview"` tag, separate from apply-time
checks; no false proof is copied into an engine-verified report. Missing legacy
receipt fields and historical `proposal` or unlabelled `submitted_proposal`
metadata establish no search-verification claim. This does not reject
valid non-minimum DNA or rerun either solver. Preview/approval bytes remain
unchanged; serialized receipts themselves are not authenticated signatures.

Source features are intentionally omitted, even though coordinates did not
change. The derived sequence carries only the explicitly supplied synthetic CDS,
translation/table, approval and non-claims. Translation preservation establishes
neither expression, splicing, regulatory function, folding nor experimental
suitability. Approval is not a laboratory/order decision.

No GUI redesign: GUI Shell and the inner agent use the same shared commands.
Solver selection is a field of that same request, not a second adapter route.
The agent must disclose an explicit conflict-search choice and obtain a fresh
preview/approval after a strategy change, never infer one from a receipt.
CLI direct forwarding, JSON `op`, MCP `op`, JS/Lua shared shell and Python `op`
reach the same operation. The agent must ask for unspecified biological inputs,
inspect actual previews and request explicit application approval; it must not
invent codons, output strings or digest values. See [CLI usage](cli.md#synthetic-sequence-design).
The existing reviewed-result draft handoff remains opt-in; execution receipts
alone do not disclose DNA or an approval digest to the model. No additional
automatic project-content forwarding or approval bypass is introduced.

## Verification

### First-Slice Verification History

On macOS with Rust/Cargo 1.100.0-beta.1, against the uncommitted implementation
on local `main` at `564ea09b`, the following checks passed after the rebase:

- `cargo test -p gentle-sequence-design --locked --offline`: 7 pure-core tests,
  including the independent exhaustive oracle and cancellation after feasibility.
- `cargo test --lib --no-default-features --locked --offline -j1 sequence_design
  -- --test-threads=1`: 6 engine/MCP tests covering exact apply, undo/redo,
  stale/tampered inputs, admission, shared shell discovery and denied confirmation.
- The same headless library test invocation with `glossary_cli_usage`: both
  glossary usage/parser checks passed, also rerun directly on the rebuilt test
  binary after the final admission guard.
- `cargo test --bin gentle_cli --no-default-features --locked --offline -j1
  sequence_design -- --test-threads=1`: direct CLI preview/apply parity passed.
- `cargo check -q --locked --offline -j1` with default features.
- `cargo check -q --locked --offline --lib --tests -j1` compiles desktop/test
  code; it does not execute native GUI workflows.
- `cargo test --no-default-features --locked --offline -j1
  --test mcp_capability_surface --test parity_matrix_freshness
  --test release_version_consistency`: 6 MCP surface, 1 canonical generated
  matrix and 4 version-consistency checks passed.
- `python3 -m unittest scripts.test_tutorial_checkouts
  scripts.test_tutorial_walkthroughs scripts.test_publish_tutorial_gui_screenshots`:
  27 tests, with 3 expected platform skips.
- `cargo fmt --all --check` and `git diff --check`.

The final admission regression re-hashes oversized previews to ensure early
rejection does not merely rely on a stale approval digest. Existing headless
unused-code and macOS large-unwind linker warnings remain visible and unchanged.
The temporary canonical-renderer helper's own Cargo target was cleaned; the
repository target and other sessions were not cleaned or interrupted.

Full-workspace tests, native GUI/Pi, native Windows, script-adapter runtime,
the declared Rust 1.85 MSRV, upstream DNA Chisel comparison and Glen's
independent review/measurements have not been run. Compile time and
focused test success are not runtime-performance or release acceptance evidence.
Session-close maintenance reported 4 OK, 2 warnings and no failures: all dirty
files are intentional task edits, and manual plan-fidelity review remains a
reminder rather than independent acceptance. The release checklist explicitly
excludes this independently versioned private crate from application rollovers.

### Conflict-Slice Verification

Verified the patched working tree based on `d221aa1accd18b8f83f686279081885b83ade3db`
on Darwin 27 arm64, Rust/Cargo 1.100.0-beta.1, default development/test features.
This is local pre-commit verification, not acceptance of a frozen/pushed package.

```sh
cargo test -p gentle-sequence-design --locked --offline -j 1
cargo test -p gentle-protocol --locked --offline -j 1 sequence_design
cargo test --lib --locked --offline -j 1 sequence_design -- --test-threads=1
cargo test --bin gentle_cli --locked --offline -j 1 sequence_design -- --test-threads=1
cargo check -q --locked --offline -j 1
cargo fmt --all --check
git diff --check
python3 -m unittest scripts.test_tutorial_checkouts scripts.test_tutorial_walkthroughs scripts.test_publish_tutorial_gui_screenshots -q
```

All 24 targeted Rust tests passed: 11 core, one protocol, ten engine/shell/MCP
and two direct-CLI/shared-parser tests. Python passed 25 tests with three expected
platform skips. Locked Cargo check, formatting and whitespace checks passed.
Root-test and CLI linking retain the known non-blocking large `__eh_frame`
warning. These CLI tests execute the forwarding/parser/engine paths in-process,
not newly packaged native executables. No full workspace, native GUI/Pi,
Windows/Linux, Rust 1.85, script-runtime or upstream DNA Chisel comparison was
run. Candidate-count assertions are functional evidence, not benchmarks.

Session-close reported four OK, two warnings and no failures: all 18 changed/new
files are intentional scoped implementation/test/contract edits, and the manual
plan-fidelity reminder is not an independent review. No generated catalog,
scientific fixture/export, Cargo dependency/profile or workflow changed. No
push, tag, dispatch or release/experimental acceptance is implied.

### Approval-Boundary Audit Follow-Up

Glen's 2026-10-05 Linux audit of `85d23ddb` found correct constrained outputs on
its tested cases, but a `needs_fix` reporting/provenance verdict. Both defects
are visible by source inspection at current `dc8c5252`: client-rehashed valid
two-edit `ATGGAGTTTTAA` could retain fabricated minimum/completeness claims,
and nested creation produced duplicate containers with split operation IDs.
Its enumeration-only acceptance does not certify the later conflict solver,
native GUI/Windows/macOS or `.12` packages; its independent evidence remains
separate from local verification.

The scoped repair separates verified output constraints from submitted search
claims, and materializes once under the outer apply. Added regressions cover
that two-edit example and rehashed frozen-stop rejection in both strategies,
legacy receipt defaults, one journal/checkpoint/container/node, loaded-source
edges, unchanged source provenance, full state undo/redo and shared-shell
receipts. No solver, dependency, budget or performance change is included.
Local builds/tests are deliberately unrun at the owner's request. Repeat the
focused protocol/engine/shell tests and direct CLI/confirmed MCP replay at the
committed repair SHA before calling these findings independently closed.

### CDS-Wide GC Verification

Verified the patched, uncommitted tree based on fetched `main`
`2f3a60f0bbf9a9ef61a02fbae6978305f5ac7df3`, Darwin 27 arm64,
Rust/Cargo 1.100.0-beta.1, default development/test features. These are local
functional checks, not a frozen binary, native GUI or independent audit verdict.

```sh
cargo test -p gentle-sequence-design --locked --offline -j 1
cargo test -p gentle-protocol --locked --offline -j 1 sequence_design
cargo test --lib --locked --offline -j 1 sequence_design -- --test-threads=1
cargo test --bin gentle_cli --locked --offline -j 1 sequence_design -- --test-threads=1
cargo test --lib --locked --offline -j 1 glossary_cli_usage -- --test-threads=1
cargo test --lib --locked --offline -j 1 cli_docs_mention_all_cli_glossary_paths -- --test-threads=1
cargo run --bin gentle_examples_docs --locked --offline -j 1 -- parity-matrix-check
python3 -m unittest scripts.test_tutorial_checkouts scripts.test_tutorial_walkthroughs scripts.test_publish_tutorial_gui_screenshots -q
cargo check -q --locked --offline -j 1
cargo fmt --all --check
git diff --check
python3 scripts/maintenance_chore.py session-close --plan docs/dna_sequence_optimization_plan.md
```

All 43 targeted Rust tests pass: 17 core, three protocol, 17 engine/shell/MCP,
three CLI and three glossary/documentation checks, with none ignored. Python
passes 25 tests with three expected platform skips, including the fast LF/CRLF
checkout gate. The canonical parity checker passes without regeneration; the
queued integration-target invocation was stopped before building to avoid
unrelated application binaries. Locked Cargo check, formatting and whitespace
checks pass. Retain the known non-blocking large `__eh_frame` linker warning.

The first root build exposed `2f3a60f0`'s `set_name(Some(String))` E0277; use the
setter's required `String` without changing materialization/provenance semantics.
The first new root run then exposed two test mistakes, not solver failures:
preview records its own journal entry, so apply needs a delta assertion; the
coupled MEF parity fixture needs exact 25% GC rather than a range admitting a
one-edit solution. The corrected core/root/CLI tests were rerun successfully.

Task changes are confined to the core GC/search/tests/README, shared protocol,
root sequence-design boundary/tests, shell/agent discoverability, CLI/MCP tests
and the existing interface/decision/plan/roadmap/changelog docs (19 files).
No dependencies, profiles, workflows, generated catalogs or scientific fixtures
change. The unrelated pre-existing deletions of
`docs/glen_sequence_design_handoff_20261005.md` and
`docs/internal_12_candidate_preparation.md` remain untouched. Session-close:
four OK, two warnings, no failures; warnings cover intentional uncommitted work
and the manual plan-fidelity reminder, not independent acceptance.

No full workspace, native Windows/Linux/GUI/Pi, Rust 1.85, JS/Lua/Python runtime,
upstream DNA Chisel comparison or runtime/RSS measurement was run. CLI tests
exercise parser/forwarding/engine paths in-process, not packaged executables.
Glen must audit the exact authorized committed candidate separately; earlier
motif-only evidence does not certify the GC slice or its search claims.

## Sliding-Window GC Slice (2026-10-06)

`gc_window={window_bp,min_basis_points,max_basis_points}` is optional. Check every
complete window in the declared CDS at a fixed one-base step, ascending local
coordinates, including frozen start/stop/protected bases. Exclude flanks and
shortened edge windows. Length must be 1..CDS length; bounds are inclusive
0..10000 basis points. No stride, arbitrary interval, strand, organism default
or percentage tolerance is inferred. Window-only design permits empty motifs;
whole-CDS GC may coexist and is never replaced by window checks.

Convert bounds once using the fixed window denominator. Malformed geometry or
bounds are invalid; an empty integer count interval is infeasible. A violating
fully frozen window also proves impossibility, after the first cancellation poll.
No per-window reachable-extrema shortcut or general infeasibility heuristic is
added. All other failures remain unresolved until the bounded search proves
infeasibility. Candidate validation constructs complete CDS prefix counts and
checks all windows; work admission adds CDS length plus one and window count to
the previous conservative per-evaluation cost. Workspace/report rows are linear
in the admitted 12,000-bp sequence limit, not a runtime-performance claim.

The common predicate covers motifs, whole-CDS GC and windows in that order.
For the earliest violating window, conflict search selects unassigned codons
capable of changing **their overlapping bases'** GC in the needed direction.
It branches over all their nonoriginal synonyms, including neutral/worsening
assignments. A feasible target must net-change GC inside that window in the
needed direction, so at least one differing unassigned codon is direction-capable.
Whole-codon GC would be an unsound substitute when a window cuts through a codon.
Full candidate/final/apply validation and edit-cost pruning remain unchanged.

Window-enabled algorithms are `synonymous_full_enumeration_window_gc_v1` and
`synonymous_conflict_search_window_gc_v1`, independent of whole-CDS GC presence.
Omitting windows preserves prior algorithm and serialized approval bytes.
`report.gc_window` retains length/count, integer bounds and every input window's
exact interval/count/satisfaction in order. Output facts exist only with actual
validated DNA. Apply caps rows before hashing, rederives all facts and checks the
approved DNA without rerunning search; false rehashed facts are refused. Receipts
continue to verify constraints, not submitted minimum/completeness claims.

Deterministic coverage includes literal synthetic MFK global-vs-local divergence,
flanks, non-codon edges, complete-CDS-sized windows, frozen/empty/malformed cases,
budget/cancellation, partial-codon GC-neutral synonyms, and a literal 36-state
L/R oracle with coupled motifs/global GC/protection. Shared shell/file, direct
CLI and confirmed MCP tests use both strategies with/without global GC.
The [CLI walkthrough](cli.md#sliding-window-gc-walkthrough) is shared with GUI
Shell/inner-agent guidance, not a generated tutorial chapter or dedicated editor.
Independent Glen, native GUI/platform, MSRV, upstream comparison and runtime
acceptance remain open; no `.12` gate is added.

### Sliding-Window Verification

Local verification on Darwin 27.0.0 arm64, Rust/Cargo 1.100.0-beta.1, used
`afbef444` plus this slice's working changes before commit, with default features
and Cargo's unoptimized test/dev profiles, not a performance baseline:

```sh
cargo test -p gentle-sequence-design --locked --offline -j 1
cargo test -p gentle-protocol sequence_design --locked --offline -j 1
cargo test --lib sequence_design --locked --offline -j 1 -- --test-threads=1
cargo test --bin gentle_cli sequence_design --locked --offline -j 1 -- --test-threads=1
python3 -m unittest scripts.test_tutorial_checkouts
cargo run --bin gentle_examples_docs --locked --offline -j 1 -- parity-matrix-check
cargo check -q --locked --offline -j 1
cargo fmt --all --check
git diff --check
```

Core **23/23**, protocol **5/5**, engine/MCP **26/26**, CLI **4/4** and checkout
**18/18** pass, with no ignored tests in these focused sets. The canonical parity
matrix is current. Locked offline check, formatting and whitespace checks pass.
The newly merged persistence test initially failed on widened `f32` JSON values
and, after fixing that, nondeterministic row order in the restriction-group
`HashMap` serialized as pairs. Compare the same JSON-text path and sort only that
map's pairs for comparison, preserving all values and every other state field;
do not change approval bytes, caches or production serialization. Final focused
tests include this correction and the executable Markdown JSON walkthrough.
The corrected persistence test passes ten separate process repetitions.
macOS debug links emit the large-unwind-table warning, without a profile change.

Full workspace, live GUI/Pi, native Windows/Linux, MSRV, upstream DNA Chisel
comparison and runtime/RSS measurements remain unrun. Earlier local-build
deferrals above describe prior patches; these focused results cover the current
integrated code, not independent scientific, usefulness or release acceptance.
Glen should rerun these commands at one committed authorized SHA and retain
source, toolchain/profile, inputs, results and limitations before issuing a verdict.
Session-close has four OK, two warnings and no failures: intentional dirty files
(including an unrelated handoff deletion) and the manual plan-fidelity reminder.
The persistence-test repair is the sole named prerequisite beyond windowed GC;
no generated catalog, fixture, dependency, profile or release changes are included.

## Remaining stages

1. Glen rechecks the approval/provenance repair at its exact integrated SHA,
   including the two-edit/non-minimum example, unverified search receipt and
   one coherent derivation through direct CLI, shared shell and confirmed MCP.
   Then independently reviews the candidates and measures representative and
   adversarial synthetic inserts, now including both explicit strategies at one
   authorized exact revision. Keep the [first-slice pinned audit](glen_sequence_design_handoff_20261005.md)
   historical rather than relabeling its SHA as new-search acceptance. Compare
   feasibility/edits/proof status separately from candidate counts and wall/RSS;
   retain coupled/new-conflict cases, budget stops and cancellation. Keep
   maintenance value, performance and biological suitability as separate claims;
   no `.12` delay or gate added.
2. Glen should now include the explicitly requested CDS-wide and sliding-window GC contracts and
   coupled constraints at one authorized revision, without treating prior
   motif-only proofs as GC acceptance. Decide usefulness before further search,
   publication or specification expansion; never substitute local checks for
   global constraints.
3. Codon harmonization/adaptation, distant repeat/hairpin interactions,
   RNA folding, other initiation/codes and DNA Chisel comparison remain explicit
   follow-ups with their own contracts, not unfinished obligations of this slice.

## Claude Follow-Up: Apply Metadata And Execution Ownership

The follow-up to `2f3a60f0` starts from local `main` `164b44ce`, preserving the
new CDS-wide GC constraints. Persist a typed, explicitly unverified preview
wrapper; reject unknown receipt/wrapper fields without changing preview digests
or the legacy default-false receipt policy. Remove the inner materialization
fork: validation and metadata serialization finish before mutation, and the
outer operation/host retain journaling, undo and stale-result checks. Name the
single container and distinguish inline synthetic imports from loaded-source
derivations. No migration of historical metadata, solver or adapter change.

Added deterministic regressions cover persisted JSON/reload counts and origins,
wrapper refusal of unlabelled/verified/unknown-field variants, unchanged approval
content, successful detached commit/undo/redo, stale structural and reopened-
project rejection without leaked output/metadata/provenance, and unknown receipt
fields. The forged two-edit and current GC tests remain in the same focused set.

No local Rust builds/tests/checks are run under the owner's restriction. The
previous GC slice's recorded passes do not validate this patch. Glen/CI should
run the following on one frozen committed revision and report its SHA/platform;
these are acceptance instructions, not local pass claims:

```sh
cargo test -p gentle-sequence-design --locked --offline -j 1
cargo test -p gentle-protocol --locked --offline -j 1 sequence_design
cargo test --lib --no-default-features --locked --offline -j 1 sequence_design -- --test-threads=1
cargo test --bin gentle_cli --no-default-features --locked --offline -j 1 sequence_design -- --test-threads=1
cargo check -q --locked --offline -j 1
cargo run --bin gentle_examples_docs --locked --offline -j 1 -- parity-matrix-check
```

Live GUI, native-platform packaging, independent search/usefulness/performance
and `.12` acceptance remain separate. No private reports are regenerated.

## GC Review Reconciliation (2026-10-06)

Claude's extra GC review examined `dc8c5252`, before the implementation in
`164b44ce`; apply ownership was subsequently repaired in `13377b12`. Recheck the
actual shared paths instead of rebuilding the proposed slice:

| Finding | Current implementation and scoped follow-up |
| --- | --- |
| F-1: false GC-only infeasibility | `conflict_search.rs` selects motif or global GC conflicts and all synonyms of direction-capable codons. Literal unrestricted eight/36-state tests already check the restriction. |
| F-2: apply missing GC validation | Core `validate_output` uses the full predicate; root apply rederives GC facts. Extend the honestly rehashed output/facts/stripped-policy test to both strategies. |
| F-3: drifting feasibility predicates | `Prepared::validate` and `Violations::satisfied` are already shared by both solvers, final checks and fresh apply validation. |
| F-4: GC-only admission | Explicit bounds already permit empty motifs. Clarify missing-both admission as `at_least_one_motif_or_explicit_gc_bounds_required`; keep invalid status, the motif cap and zero evaluated candidates. |
| F-5: approval/policy identity | GC algorithms, the GC-specific nonclaim and complete optional request/report facts are already bound. `decisions.md` owns the CDS/integer/approval invariants; legacy GC-absent feasible previews remain unchanged. |
| F-6: cheap negative proof | `PreparedGc` already separates malformed bounds, empty count windows and disjoint relaxed extrema before candidate evaluation. Extend malformed-bound status assertions to both strategies. |
| F-7: work admission | Preparation already adds one CDS-only GC pass. Replace the loose budget comparison with exact assertions for motif-only, motif+GC and GC-only requests on a flanked synthetic CDS, under both solvers. This is a work-budget contract, not timing evidence. |
| F-8: discovery/parity | Capability request text, inner-agent guidance and the CLI GC walkthrough already describe the policy and nonclaims. No new routes/progress fields or generated rows are needed; Glen/CI still checks the canonical matrix at the exact follow-up SHA. |

The original D-1 example needs care: for a 21-bp CDS, bounds 4762..5238 give
minimum 11 and maximum 10, an empty count window. Bounds 4761..5239 permit 10..11.
The existing integer-only endpoint test checks both examples, accepted endpoints
and one-base-outside refusal. Existing frozen/extrema, coupled motif/GC, both-
strategy oracle, absent/unknown output, legacy-byte and cancellation tests stay
in scope; no new solver, constraint framework or heuristics are introduced.

The follow-up command set above includes the pure core and canonical parity
check. Rust execution remains deferred under the owner's no-local-build rule;
neither static review nor `164b44ce`'s recorded passes certify this patch.
Keep independent usefulness, native, performance and release acceptance pending.
