# Bounded synthetic sequence design

Updated 2026-10-05, implementation started from `fc4458c3`.
Rebased during implementation onto fetched `main` at
`4a71d5bd15f6b0da60709ee79a431668756d110f`, then onto local `main` at
`564ea09be48149649748d620f57973d5f87cfff7` after its dependency/Help update.
The latter was ahead of the fetched remote at integration time.
Owner authorized work from the previously Claude-reviewed proposal. This is
experimental development outside the `.12` release gate, not release acceptance,
full DNA Chisel compatibility, or a measured performance improvement.
Owner-requested follow-up on 2026-10-05 adds the explicit conflict-search slice
below, based on integrated local `main` `d221aa1a`. Glen's first-slice audit is
still pending; this addition does not substitute local tests for his usefulness,
performance or independent correctness verdict.
The earlier Claude review covers the original bounded proposal; no additional
Claude review of this conflict-search implementation is claimed.

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

The reviewed proposal reduced four constraint families to motif removal.
Translation/protection are mutation invariants, not optional objectives. GC,
codon adaptation, repeats/hairpins, folding and a general specification/plugin
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

| Outcome | Meaning |
| --- | --- |
| `invalid` / `unsupported` | Admission/resource rule failed; no applicable preview. |
| `feasible` | Fully validated exact output; only complete enumeration/conflict search or unchanged feasible input proves minimum edits. A budget stop can retain a feasible candidate with optimization incomplete. |
| `search_exhausted` | Budget ended without a feasible candidate; unresolved, not infeasible. |
| `proven_infeasible` | An initial match lies wholly in frozen bases, or complete enumeration/conflict search found no admissible variant. |
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
It never reruns optimization. The exact unused output ID is created in detached
execution, checked before one undoable commit, with source lineage and retained
`gentle.dna_sequence_design_receipt.v1`. Original DNA remains unchanged.
One outer apply owns the journal/undo entry, derived lineage node and singleton
container; their operation IDs agree, without nested public creation hooks.

The approval digest identifies reviewed content, not authenticated search
provenance. A caller can change a valid output and recompute that digest, so
apply must not certify the preview's completeness/minimum, counts or reason.
The engine receipt records `output_constraints_verified=true` and
`search_claims_verified=false`. The exact approved preview is retained only as
`submitted_proposal`, separate from apply-time checks; no false proof is copied
into an engine-verified report. Missing legacy receipt fields and historical
`proposal` metadata establish no search-verification claim. This does not reject
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
2. Use that evidence to decide whether further search changes, independent
   publication or richer specifications merit their cost. GC constraints need
   their own whole-sequence/window and threshold contract before implementation;
   do not replace global constraints with assumed local checks.
3. GC, codon harmonization/adaptation, distant repeat/hairpin interactions,
   RNA folding, other initiation/codes and DNA Chisel comparison remain explicit
   follow-ups with their own contracts, not unfinished obligations of this slice.
