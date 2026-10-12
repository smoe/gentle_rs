# Sequence-design review history (2026-10-05)

Archived from `docs/dna_sequence_optimization_plan.md` in commit
`7669ed1c`. The status/stages below describe the proposal before owner-authorized
implementation; they are retained verbatim as review history, not current tasks.
The [active bounded contract](dna_sequence_optimization_plan.md) supersedes them.

---

# Rust sequence design proposal after Claude review

Updated 2026-10-05 against `f523a757cd2a2861196b1a5c0b36df67983ce0e6`.
Status: proposed future work, not approved implementation or a `.12` release gate.
The owner supplied the Claude review; Codex checked the cited local behavior
without invoking another consultation, running tests, or changing runtime code.

The recommendation is a small Rust component for synonymous motif removal from
synthetic coding inserts, not a full DNA Chisel clone. GC constraints and a
general specification framework are deferred. Independent maintenance remains
the design aim; Rust alone establishes neither a speedup nor a reduction in
GENtle's build memory.

## Original proposal and review reconciliation

The [original Codex proposal](dna_sequence_optimization_claude_prompt.md) called
for an independent crate with translation preservation, protected bases, motif
avoidance and GC constraints, followed by bounded search and engine integration.
Claude recommended reducing the first slice to motif removal, clarifying
initiation/cancellation/negative-result semantics, and postponing crate
extraction until an engine-internal prototype demonstrates value.

| Review finding | Codex disposition |
| --- | --- |
| Four initial constraint families are too much. | Accept the scope reduction. Translation and protected bases remain mandatory mutation invariants; motif removal is the sole initial design goal. Defer global/window GC. |
| The codon mapping does not describe initiation. | Confirmed. `assets/codon_tables.json` has elongation mappings and stop markers, not start-codon sets. State the limitation and narrow initial admission below. |
| Ambiguous codons always return unknown. | Correct this claim: `codon2aa` resolves degenerate codons when every expansion has the same translation. Existing tests assert `GCN -> A` and `GAN -> UNKNOWN_CODON`. Reject ambiguous input here to simplify redesign, not because GENtle cannot translate any ambiguous codon. |
| An engine-internal module is sufficient initially. | Valid alternative, not a correctness finding. Prefer a minimal independent crate to serve the requested standalone-maintenance aim; defer publication, a separate repository and general extensibility. Use an internal pure module as a fallback if Stage A identifies a concrete ownership or integration obstacle. |
| Cancellation is outside deterministic completion. | Accept. External cancellation can occur at different work counts and must never produce an applicable proposal, even if a valid intermediate candidate exists. |
| Exhaustive infeasibility proof is rarely practical. | Accept the clarification. Exhaustive enumeration is a tiny-space correctness oracle, not the production strategy for realistic inserts. Bounded-search failure means unresolved, not infeasible. |
| Reverse-translation diagnostics were omitted. | Add `reverse_translation_choice_diagnostics` as reporting precedent. Preferred/alternative codon counts are not nucleotide edit counts and must not be relabeled as such. |
| Stage A lacks an owner and checkable output. | Assign proposed ownership and require a written semantic contract with resolved admission and resource limits before implementation. |

The existing exact-product approval mechanism is confirmed in
`src/engine/analysis/regulatory_fragment_panel/materialization.rs`. It binds
products in a digest, rejects altered or stale proposals, creates the approved
sequence strings in detached execution, and checks their hashes before commit.
It is a useful precedent, not evidence that the new solver integration already
exists. Its source re-planning must not become optimization re-execution during
approval of a redesigned sequence.

## First useful capability

Given an explicitly selected synthetic insert, remove specified finite-length
IUPAC motifs using synonymous substitutions while preserving the supplied
protein translation, frozen flanks and protected bases. Never recode natural
PATZ1 assay targets or change current reverse-translation/primer behavior
implicitly. Minimize changed nucleotides as a secondary objective; report when
that objective has not been proven optimal.

Initial admission is unambiguous A/C/G/T DNA with one forward, contiguous CDS,
length divisible by three, an ATG start, a terminal stop under explicitly
resolved standard genetic code 1, and no internal stops or translation
exceptions. Freeze the literal start and terminal stop codons. This does not
infer a biologically authentic CDS or model initiation efficiency. Alternative
initiators and other genetic codes remain unsupported in this first slice,
rather than being guessed from an elongation table.

The root engine resolves context and converts the existing table into explicit
codon-to-residue entries; do not pass the asset's positional 64-character string
without declaring its ordering. Bind the resolved identity and mapping content.
The core consumes that validated mapping without depending on the root crate,
but the first supported contract remains standard code 1 with the frozen ATG
start. Do not extract all amino-acid code merely to supply this input.

Freeze everything outside the CDS and intersect synonymous choices with every
protected base inside it. Motif evaluation covers the full supplied sequence,
including CDS/flank boundaries, overlapping matches and explicitly selected
strands. Define reverse-complement and palindromic duplicate semantics. A match
entirely within frozen bases is a specific contradiction proof. Partial CDS,
ambiguous bases, joined/overlapping CDS, reverse-oriented requests, circular
topology and indels receive explicit unsupported/invalid outcomes.

Preserving translation does not establish preservation of expression,
splicing, regulatory function, folding or experimental suitability. Report
these non-claims alongside the proposed sequence, not only in documentation.

## Core and search boundaries

Prefer a tiny independent workspace crate with a provisional neutral name and
independent version. Keep it private/unpublished until the first useful API has
been reviewed. Its initial surface is validated input, legal codon choices,
motif evaluation, bounded search, outcome records and a nucleotide edit script.
Do not build a plugin system, GenBank design language or generic weighted
optimizer in this slice. Keep GENtle state, GUI, files, network and external
tools outside it; optional serialization is a separate interface concern.

Implement full evaluators and complete enumeration for tiny spaces first.
Validate coupled edits with independently enumerated fixtures; merely spanning
two codons does not prove that a motif needs a multi-codon repair. Add bounded
deterministic search for larger inputs only after those fixtures establish the
correctness baseline. Freeze candidate order, tie-breaks, algorithm version,
work budget and any future RNG algorithm/seed. Admit resource limits before
allocation, including sequence length, motif work and search-space counting.

Do not cache or localize evaluation initially. Every accepted output receives a
fresh full-sequence check of translation, protected bases and every requested
motif. Future GC bounds are a second slice with their own threshold/window-edge
contract. Repeat/hairpin specifications may depend on pairs of distant regions;
do not assume all later evaluators can be reduced to one local window.

| Outcome | Meaning and proposal eligibility |
| --- | --- |
| Invalid or unsupported | Admission failed; no proposal. |
| Feasible | Full validation passed. A proposal may be produced with edit count and explicit optimization-completion/optimality status. |
| Search exhausted | Budget ended without a feasible result; no proposal and no infeasibility claim. |
| Proven infeasible | A recorded sound contradiction or completed tiny-space enumeration excluded every admitted variant; no proposal. |
| Cancelled | Externally interrupted, potentially at a nondeterministic work count; diagnostics only, never an applicable proposal. |

A normal work-budget stop with a fully validated best candidate is `feasible`
with optimization incomplete, not cancellation or proof of minimum edits.
Preparation/search never mutates the caller's original input or project state.

## Delivery stages

Proposed owners below describe future work, not dispatches or acceptance.

| Stage | Owner and checkable output |
| --- | --- |
| A: semantic decision | Codex drafts the bounded contract; the owner decides scope/component ownership, with biological review for the claims. Resolve numeric admission/work limits, mapping/ordering, motif/coordinate semantics, dependency/license alternatives and supported Rust version in this scoped design document. Record active invariants in `docs/decisions.md` only once adopted. |
| B: smallest proof | Implementer builds the pure core and tiny-space oracle. Independent tests cover translation, start/stop freezing, protected bases, both strands, overlapping sites, flanking-boundary sites, unchanged feasible input and frozen-site contradiction. A crafted coupled-edit case must falsify a naive greedy repair. |
| C: bounded useful search | Implementer adds larger-space search and records all outcomes above. Compare feasibility and edit optimality against exhaustive tiny cases; test budget stops and cancellation after a valid candidate is found. No global-optimum or complete-solver claim for bounded runs. |
| D: GENtle integration | Codex adds only the engine request/report/planning/materialization and shared-shell route required for this capability. Test tampering, source/table/protection drift, exact approved output, undo/lineage and adapter reachability. No GUI redesign. |
| E: ownership and measurements | Glen measures exact-revision representative and adversarial inserts. Assess maintenance and standalone API readiness before publication or further specifications; performance results are separate from biological validity. |

GENtle proposals bind exact source/output DNA, CDS geometry, frozen bases, motif
specifications, mapping identity/content, solver version/configuration, final
evaluations and edits. Materialization checks current inputs and validates the
approved output without rerunning search, then creates a derived sequence in
one undoable transaction. Labels do not authorize edits. Avoid carrying stale
biological annotation claims forward simply because coordinates are unchanged.

Python/DNA Chisel comparison is optional development tooling only. It is not a
dependency of the Rust core, normal CI, GENtle packages or runtime. Any future
separate comparison job needs its own pinned environment and explicit scope;
do not silently add a Python platform matrix. Compare specified semantics, not
identical final DNA among multiple valid solutions, and retain independent
evaluators as the correctness oracle.

CI/core tests should cover Linux, macOS and Windows without pulling in the GUI.
Fixtures require origin/recreation/use documentation. Any committed byte-exact
projection or hash-bound artifact also requires the scoped LF/CRLF checkout
policy. No local test, build, benchmark, commit or implementation has been run
as part of this planning revision.

## Deferred work

GC constraints, codon adaptation/harmonization, repeats/hairpins, RNA folding,
external aligners, broader genetic codes/initiation models, circular/joined or
reverse-strand CDS, indels, arbitrary plugins, GenBank design annotations and
full DNA Chisel compatibility remain separate future decisions. The `.12`
roadmap and release gates are unchanged.
