# Rust DNA sequence optimization review prompt

Prepared 2026-10-04 against GENtle `main` at
`f523a757cd2a2861196b1a5c0b36df67983ce0e6`.
Status: original Codex proposal, retained unchanged below for comparison.
The owner supplied Claude's review on 2026-10-05; see the
[revised proposal and review reconciliation](dna_sequence_optimization_plan.md).
Neither document is approval for implementation or a release commitment.

## Review request

Please perform a read-only architectural review of the proposal below. The
repository is `/Users/u005069/GitHub/gentle_rs`. Read its `AGENTS.md` and relevant
guidance before inspecting code. Do not edit files, run builds or tests, update
dependencies, invoke other agents, commit, push, or change project status.

The owner asked us to inspect DNA Chisel, identify concepts missing from GENtle,
and consider either integration or a Rust reimplementation maintained as a
separate crate. Codex recommends a bounded Rust implementation of the useful
concepts, not a promise of complete Python API or solver-output compatibility.
Please challenge that recommendation, especially its maintenance cost and
whether the first useful capability justifies a new solver.

GENtle's `.12` priorities remain responsiveness and release acceptance. This
proposal is future work, not an additional `.12` blocker. Existing primer-panel
and biological-context behavior must remain unchanged.

## Evidence already gathered

DNA Chisel distinguishes hard constraints from optimization objectives and
represents legal sequence variants with a mutation space. Its specifications
support regional evaluation and localization. These are relevant concepts for
joint sequence redesign, rather than merely filtering existing candidates.
See the [core classes](https://edinburgh-genome-foundry.github.io/DnaChisel/ref/core_classes.html).

Its built-ins include translation preservation, protected sequence, motif
avoidance, GC bounds, and several codon-use objectives. Those establish a useful
comparison vocabulary, not an obligation to implement every specification.
See the [built-in specifications](https://edinburgh-genome-foundry.github.io/DnaChisel/ref/builtin_specifications.html).

The [upstream project](https://github.com/Edinburgh-Genome-Foundry/DnaChisel)
provides Python source and an MIT license. Preserve applicable copyright and
license notices if source or test material is adapted; verify the exact pinned
revision and distribution before importing anything. Do not assume permission
to reuse unrelated biological datasets or imply upstream endorsement.

Relevant local evidence, with paths relative to the repository above:

| File or symbol | Observation and implication |
| --- | --- |
| `src/engine/ops/operation_handlers.rs`: `choose_back_translation_codon`, `build_reverse_translated_coding_sequence` | Existing reverse translation selects codons sequentially using preferences and local Tm/GC penalties. Extend neither its semantics nor its callers silently. |
| Same file: `FilterByDesignConstraints`, `ScoreCandidateSetWeightedObjective`, `ParetoFrontierCandidateSet` | Candidate filtering and ranking already exist. A redesign solver should not duplicate them or be presented as a new primer designer. |
| `src/amino_acids.rs`, `assets/codon_tables.json` | Codon interpretation is currently rooted in GENtle. The independent crate must not depend on the root crate merely to obtain a table. |
| `src/engine/protocol.rs`: `PrimerDesignSideConstraint`, `PrimerDesignScoreTerm` | Primer constraints and score explanations provide integration precedent, not the mutation-space implementation. |
| `src/engine/ops/reporter_ops.rs` | Existing reporter recommendation has hard rejection and soft ranking. Check reuse opportunities without expanding this task into a reporter refactor. |
| `Cargo.toml`, `crates/gentle-engine/Cargo.toml`, `crates/gentle-engine/src/lib.rs` | Preserve workspace dependency direction and headless reuse. |
| `docs/protocol.md` | Inspect existing preview, approval, materialization and provenance contracts before proposing a new operation family. Not all existing workflows use identical approval mechanics. |
| `docs/architecture.md`, `docs/decisions.md`, `docs/biological_extension_guide.md`, `docs/testing.md`, `docs/roadmap.md` | Shared engine ownership, biological non-claims, portable reports, fixture provenance, byte-exact checkout checks, and scoped future work remain binding. |

These observations do not prove that no suitable Rust crate already exists.
Check alternatives before committing to ownership of a new one.

## Original Codex plan

### 1. Define the first biological contract

The first use case is redesigning an explicitly selected synthetic linear coding
insert, not rewriting natural assay targets such as PATZ1 transcripts.

Accept unambiguous A/C/G/T DNA with one contiguous, complete, forward-oriented
CDS and an explicitly resolved genetic code. Initially permit synonymous codon
substitutions only; freeze sequence outside the CDS. Freeze the original start
and terminal stop codons and any user-protected bases. Reject partial CDS,
ambiguous sequence, internal stops, unsupported translation exceptions,
overlapping CDS, joined locations, circular topology, and indels explicitly.
The engine may later normalize other orientations, but must not imply initial
support for them.

Provide four hard constraints: unchanged amino-acid sequence, unchanged protected
bases, absence of specified finite-length IUPAC motifs on explicitly selected
strands, and requested global or sliding-window GC bounds. Motif checks cover
the whole supplied sequence, including junctions with frozen flanks. If an
unavoidable site is entirely frozen, report that limitation rather than editing
the flank or pretending the input is solvable.

Define zero-based half-open coordinates, complete-window GC edge behavior,
inclusive GC thresholds, motif overlap reporting, reverse-complement matching,
and palindromic duplicate handling before implementing search. Use integer GC
counts and rational bounds where practical to avoid threshold rounding drift.

The sole initial soft objective is minimizing nucleotide substitutions from the
original input. A heuristic result is not automatically the minimum-edit result.
Translation preservation does not establish preservation of expression,
splicing, regulatory activity, folding, or experimental performance.

### 2. Create an independent core crate

Use a provisional neutral crate name, subject to an availability and dependency
audit. Keep its API and release version independent of GENtle's internal release
number. Developing it initially inside the workspace does not require creating
or publishing a separate repository now.

The crate owns validated sequence input, regions, built-in specifications,
evaluation records, legal codon choices, bounded search, and an edit report.
Start with a small typed set of built-ins rather than a general plugin DSL or
arbitrary user-code hooks. Keep any serialization support optional.

GENtle supplies an explicit canonical 64-codon translation mapping and the
resolved table identity; the core validates and uses that mapping. Standalone
callers must supply equivalent explicit inputs. This avoids a dependency cycle
and an unrelated extraction of all amino-acid logic. Check whether this boundary
is sufficient for genetic-code and initiation semantics.

Do not add GUI, filesystem, network, Python, external aligner, or project-state
dependencies to the core. Bound input size, motif count/length, window work,
search evaluations, and memory before allocation or enumeration. Choose concrete
initial limits from the intended insert workload during the first implementation
stage; do not leave production admission unbounded.

### 3. Establish correctness before local optimization

Construct each codon's legal synonymous choices, intersecting them with every
protected base. Preserve the caller's original input immutably; preparation must
not silently rewrite it. Treat choices for whole codons as units rather than
independent nucleotides that could break translation.

First implement complete evaluators and exhaustive enumeration for small legal
spaces. Use those as an independent oracle for later search and edit-count
optimality. Then add bounded deterministic search over larger spaces, including
multi-codon moves where single-codon repair gets trapped. Start with full
evaluation; introduce local invalidation/caching only after profiling and tests
prove equivalence, including motifs and GC windows crossing edit boundaries.

Search feasibility first, then reduce edit count without sacrificing any hard
constraint. Fix ordering, tie-breaks, solver version, evaluation budget and, if
used, the RNG algorithm and seed. Ordinary completion should use work budgets,
not wall-clock timing, as its reproducible stopping rule.

Return distinct outcomes for invalid input, unsupported requests, a feasible
result, resource/search exhaustion, and cancellation. Claim proven infeasibility
only after complete enumeration of the admitted legal space or a specific sound
contradiction proof, never from failure of a local/random search. A feasible
result may have incomplete optimization; report that separately. Cancelled or
partial runs must not silently become applicable proposals.

Before accepting a result, re-evaluate every hard constraint against the complete
final sequence without relying on cached or assumed passes. Keep the original
state unchanged on every failure.

### 4. Integrate through the GENtle engine

After the core passes its correctness gates, add a narrow engine-owned planning
operation, a portable report, and a shared-shell route. Proposed names and schema
versions must follow the repository's existing conventions; do not finalize a
new protocol during this review.

Planning binds source sequence bytes, CDS geometry, frozen bases, normalized
specifications, resolved genetic-code/catalog inputs, solver configuration and
version, evaluations, the edit script, and exact proposed sequence bytes. Treat
labels as explanations, not mutation authority. Preserve relevant source
annotations without claiming that every annotation remains biologically valid.

Materialization requires the exact approved proposal digest, checks currency and
all necessary inputs, validates the proposed output, and creates a derived
sequence with lineage in one undoable transaction. It does not rerun search at
approval time, overwrite the source, or invoke a new agent. The same approved
bytes must reach every adapter. No automatic change to reverse translation,
primer design, or assay readiness is part of this integration.

### 5. Validate in stages

| Stage | Deliverable and admission condition |
| --- | --- |
| A | Specify semantics, alternatives, translation boundary, licensing, minimum supported Rust version and resource limits. Agree the minimal useful scope before implementation. |
| B | Standalone evaluators and mutation-space construction. Hand-crafted cases plus generated tiny spaces prove translation/protection and enumeration correctness. |
| C | Bounded solver. Compare against exhaustive tiny-space results; test infeasibility, exhaustion, cancellation, traps needing multi-codon edits, and final full validation. |
| D | GENtle proposal/report/materialization integration. Test source drift, tampered proposals, input changes, undo/lineage, and shell/protocol parity before GUI work. |
| E | Glen measures representative inserts and adversarial cases at an exact revision. Admit performance/cache changes only with semantic equivalence evidence. Consider independent publication after API review. |

Use a pinned DNA Chisel version as an optional development comparison, not a
required packaged runtime or the only correctness oracle. Compare agreed
constraint semantics, translations, protected bases and reported feasibility;
do not demand identical optimized sequences when several valid solutions exist.
Resolve disagreements against the written contract and independent evaluators.

Every fixture needs source/recreation/use documentation. Any committed
hash-bound or byte-exact generated artifact needs scoped Git attributes and the
existing LF/CRLF checkout gate. CI must exercise the core on Linux, macOS and
Windows. Standalone/core-only builds must not pull in the desktop application.

This task prepares a review plan only. No local builds, tests, benchmarks,
implementation, publication or commits are authorized by the plan itself.

### 6. Defer adjacent features deliberately

Defer CAI maximization, target codon distributions and harmonization, repeat or
hairpin objectives, RNA folding, BLAST/Bowtie-backed specifications, arbitrary
regular expressions, indels, overlapping CDS, circular support, GenBank design
annotation import/export, graphical design controls, and full upstream API
compatibility. Add each only for a demonstrated use case and a separate semantic
and validation contract. Rust does not by itself establish speed gains or lower
GENtle build memory.

## Requested critique

Please return concise findings in risk order, grounded in current repository
paths/symbols, followed by a revised minimal plan. Specifically assess:

1. Is a Rust crate justified now, or would a pinned optional Python adapter be a
   better first experiment? Separate short-term delivery from long-term ownership.
2. Is the first use case still too broad? Identify the smallest slice that is
   useful without making false biological or solver-completeness claims.
3. Does supplying a 64-codon mapping plus table identity avoid duplication while
   correctly handling translation, starts/stops and unsupported exceptions?
4. Are the search, mutation-space, boundary, failure-state and reproducibility
   contracts sufficient? Identify contradictions rather than assuming Rust
   reproduces DNA Chisel's semantics automatically.
5. Can the existing engine proposal/history/report machinery support exact-byte
   approval without a broad refactor or new root-crate dependency?
6. Which test or acceptance gate would most cheaply falsify this design? Which
   deferred concept is actually necessary for the first useful result?

Explicitly distinguish confirmed code observations from recommendations and
unknowns. Do not mark the proposal approved, implemented, tested, scientifically
validated, or part of the `.12` release gate.
