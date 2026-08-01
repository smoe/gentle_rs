# GENtle Decisions

## 2026-08-01: derive web and print gene-set reports from one record

Publication-style multi-gene analyses use one portable manifest and one
resolved semantic record. Responsive HTML remains the main browsable report;
the PDF is a printable companion containing the full primer list and all
selected figures. This retains the useful PARK7-style overview while removing
hard-coded gene lists and independently maintained HTML/PDF assembly scripts.

This file records durable implementation constraints and architectural
decisions. It is intentionally short: if an entry becomes a backlog item, move
that work to [`roadmap.md`](roadmap.md); if it describes completed history,
move the outcome to [`CHANGELOG.md`](CHANGELOG.md).

Status values:

- `active`: currently governs implementation.
- `superseded`: kept for history; no longer governs new work.
- `deferred`: accepted direction, but not currently active work.

## DEC-001: Single Shared Engine

Status: active

GENtle uses one deterministic engine contract across GUI, CLI, shared shell,
JS, Lua, Python wrapper, MCP, and agent-facing routes. Adapters must not fork
biology or business logic.

## DEC-002: Protocol-First Machine Contracts

Status: active

Machine-facing routes use structured, versioned request/result records. Durable
schemas and wire contracts live in [`protocol.md`](protocol.md), not in the
roadmap.

## DEC-003: Thin Adapter Boundary

Status: active

Frontend and integration adapters should parse, route, present, and collect
artifacts. Core decisions, reports, findings, and explanation records should be
engine-owned and portable whenever practical.

## DEC-004: Roadmap/Changelog/Decision Split

Status: active

`roadmap.md` is for next work only. Completed work goes to
[`CHANGELOG.md`](CHANGELOG.md). Durable constraints go to this file. A session
should be able to start by reading only `roadmap.md` in under two minutes.

## DEC-005: No Roadmap `Done` Entries

Status: active

Do not add a `Done` entry to `roadmap.md`. Move completed work directly to
`CHANGELOG.md` in the same session. This prevents the roadmap from growing back
into a mixed history/planning document.

## DEC-006: Workspace Extraction Order

Status: active

The intended crate extraction order remains:
`gentle-protocol -> gentle-engine -> gentle-render -> gentle-shell -> gentle-gui`.
First-wave extraction avoids per-feature micro-crates; related analysis logic
should remain together until the engine boundary is stable.

## DEC-007: Root-Crate Compatibility During Extraction

Status: active

While production code still partly lives in the root crate, extracted crates
must expose stable contracts without forcing all callers through a one-shot
import rewrite. Root modules may re-export extracted contracts during the
transition.

## DEC-008: Shell Dispatch Split Pattern

Status: active

Shared shell parsing and execution should remain one behavior path for GUI
Shell, `gentle_cli shell`, MCP, JS/Lua wrapper helpers, and agent execution.
Large command families may split into helper dispatch functions, but those
helpers must preserve shared parser/executor semantics.

## DEC-009: Stack-Hardening Policy

Status: active

Stack-overflow fixes should reduce dispatcher frame depth without forking
behavior. Prefer narrow helper dispatch and expanded-stack workers for
confirmed small-stack failures in shell/engine routes.

## DEC-010: `#[inline(never)]` Helper Dispatch

Status: active

Until the monolithic root dispatchers are fully decomposed, stack-sensitive
command or operation families may dispatch through dedicated `#[inline(never)]`
helpers before entering broader match frames. Inner match branches should
delegate back to the same helper so there remains one implementation path.

## DEC-011: UI-Intent Shared Catalog

Status: active

UI-intent targets and discoverability metadata are shared catalog data.
Menus, command palette, shell/agent routes, MCP discovery, and ClawBio/OpenClaw
handoff should consume the same target catalog instead of hard-coding parallel
target lists.

## DEC-012: GUI Intent Handlers Stay Thin

Status: active

GUI `ui open|focus|close|selection ...` handlers may open/focus/close host
windows, dialogs, and active viewer selections, but must route target
resolution through shared `ui ...` shell command contracts where possible.

## DEC-013: Mutating Agent/MCP Safety

Status: active

Mutating `op` and `workflow` execution through agent/MCP paths requires
explicit confirmation. Agent-suggested commands must not recursively invoke
`agents ask`, `agents plan`, or `agents execute-plan`, including through macro
expansion paths.

## DEC-014: Screenshot Capture Policy

Status: active

Screenshot capture paths stay disabled unless project policy explicitly
approves them. Documentation should prefer deterministic SVG/export routes
while screenshot execution is policy-disabled.

## DEC-015: ClawBio/OpenClaw Boundary

Status: active

The ClawBio/OpenClaw skill scaffold wraps deterministic `gentle_cli` routes and
writes reproducibility bundles. It should surface shared capabilities,
suggested actions, and UI-intent handoffs without creating ClawBio-only biology
logic.

## DEC-016: External AI/Automation Deployment

Status: active

For tool-driven external AI deployment, prefer the published GHCR image in
headless `mcp` mode with explicit project/state mounts and stdio communication
instead of depending on browser GUI containers.

## DEC-017: Fixture Provenance

Status: active

Every committed fixture under `test_files/` or `tests/` must document origin,
deterministic recreation/retrieval steps, and where GENtle uses it.

## DEC-018: Source Documentation Policy

Status: active

Public modules and key public records should carry concise rustdoc explaining
purpose, invariants, and non-obvious behavior. Prefer why/invariant/edge-case
documentation over line-by-line narration.

## DEC-019: Primer/Oligo Material Identity

Status: active

A planned primer/probe or in-silico design artifact is not physical stock.
GUI, CLI, JS, Lua, Python, MCP, and agent routes must preserve the distinction
between design, reviewed order, received material, local availability, and
consumption.

## DEC-020: Inline/Stateless Sequence Inspection

Status: active

Read-only operations that need only sequence letters, optional topology, and an
optional span should be state-optional. Promoting such results into project
state remains an explicit second step.

## DEC-021: Helper-Construct Terminology Migration

Status: deferred

The legacy "helper genome" vocabulary should eventually move toward
"helper construct" semantics. The migration should be atomic across contracts,
docs, GUI, and examples to avoid mixed terminology.

## DEC-022: Primer3 Backend Parity

Status: deferred

Primer3 is the intended external backend parity layer for primer design, but
deeper constraint mapping, fixture coverage, and backend equivalence checks
should land through shared engine/report contracts rather than GUI-only code.

## DEC-023: Presentation i18n Does Not Localise Machine Contracts

Status: active

GUI labels and dialog text may be translated at runtime, but shared shell
commands, protocol schema fields, saved project records, adapter payloads,
agent response schemas, and scientific identifiers remain deterministic
English. Localised text is presentation only.

## DEC-024: Introspection Composes Facts And Capabilities

Status: active

Words-only and agent-facing introspection uses additive fact domains rather
than renaming existing fact vocabulary. Existing project facts keep their
`sequence.*`, `report.*`, and related names with `domain: "project"`;
view/host/config facts use the same registry shape. Headless contexts must
project `ui.host_available=false` explicitly so GUI/view-intent readiness fails
with a named unmet fact rather than ambiguous absence. Capability readiness is
computed by instantiating templated fact atoms from declared arguments and then
using the shared fact evaluator; unbound catalog readiness reports missing
argument bindings as `unknown`. Host-domain availability facts must project
from existing deterministic preflight/probe logic rather than creating a
parallel detector.

## DEC-025: Legacy SHA-1 Download Verification Uses External Tools

Status: active

GENtle must not use an in-process SHA-1 implementation to verify external
downloads. If an upstream source or compatibility manifest still offers only a
SHA-1-style integrity field, GENtle may delegate that legacy digest to platform
tools such as `sha1sum`, `shasum -a 1`, or Windows `certutil`. If no suitable
tool is available, the file is not SHA-1 verified; preparation may still rely on
basic corruption checks such as successful gzip decompression, size checks, and
parser validation. `GENTLE_SHA1_TOOL` may point to an explicit tool, and
`GENTLE_DISABLE_LEGACY_SHA1=1` intentionally disables the legacy check. Legacy
SHA-1-shaped local provenance fields are compatibility fields; regenerated
local fingerprints should use algorithm-prefixed stronger digests such as
`sha256:<hex>` without requiring an in-process SHA-1 dependency.

## DEC-026: Long GUI Jobs Use Optimistic Engine Snapshots

Status: active

Long-running GUI network and computation jobs must not hold the shared
`GentleEngine` read or write lock for their full run. A mutating worker forks
the engine state, journal, operation counter, revisions, and history limit,
executes against that detached snapshot, and commits with a short write-lock
swap only if structural project data and journal history are unchanged.
Inherited undo/redo checkpoints stay in the live engine; commit appends only
checkpoints created by detached operations. Concurrent viewport/display edits
are preserved across the commit; disjoint UI/session metadata is merged into
the result and its undo checkpoints, while conflicting metadata makes the
result stale. Other stale results never overwrite intervening edits.
Persisted-state mutation revision is separate from structural and read-only
execution identity so GUI dirty checks remain constant-time without marking
inspections as edits.

Read-only workers use a history-free engine clone and do not commit it. A raw
`GentleEngine::clone` in a background worker is prohibited because it can clone
hundreds of full-state history checkpoints. Persisted metadata mutations outside
`apply` must advance the mutation revision and clear redo; otherwise a redo
checkpoint captured before the side-edit can overwrite that metadata.

## DEC-027: RNA Allele Evidence Does Not Invent Phase

Status: active

RNA-read allele screens may use reviewed local variants to extend reference
hashes only when the variant can be projected to an explicit transcript
coordinate and its reference allele matches that transcript. Unphased
genotypes remain allele-level evidence, and disconnected phase sets remain
independent blocks; neither may be combined into a synthetic gene-wide
haplotype. Reports must distinguish physical reads, fragments, and evaluated
evidence observations. Structural rescue matches are non-blocking candidate
evidence and must not be labelled as novel exons, retained introns, fusions, or
resolved isoforms without a separate alignment and junction-validation layer.
Annotated-junction and exon-intron-boundary catalogs index only k-mers that
cross an explicitly declared boundary. Retained-intron and intronic candidates
remain reportable, but same-read or same-fragment support from an annotated
splice junction for that gene is recorded separately as an RNA anchor. Such an
anchor strengthens an RNA-derived interpretation without excluding pre-mRNA or
genomic-DNA contamination and is not itself a retained-intron call.

## DEC-028: Construct Reasoning Outputs Are Fingerprinted Snapshots

Status: active

Construct-reasoning graphs persist their emitted inspection actions as part of
the graph snapshot. Loading or normalizing a non-empty action list must not
silently rerun current recommendation rules. New graphs fingerprint the full
source sequence/feature snapshot and normalized objective with SHA-256 and
record the reasoning rule-set version; readers compare those identities to
report `current`, `stale`, or `unknown`. A fingerprinted stale graph remains
inspectable but its executable inspection actions must be refused until its
reasoning build route is rerun to produce a new snapshot.

Task-specific reasoning keeps intrinsic evidence severity separate from
objective applicability and effective priority. Typed
`ConstructObjective.intended_tasks` is authoritative when present; legacy
free text may conservatively infer positive applicability but must not infer
non-applicability. Curated repeat provenance belongs in typed repeat annotation
and repeat-family records, with legacy note parsing retained only for backward
compatibility.

## DEC-029: Primer Specificity Is Target-Space And Database-Content Bound

Status: active

Primer-specificity evidence must name its biological target space. Genomic-DNA
screening assesses genomic carryover and unintended genomic products;
transcriptome/cDNA screening assesses cross-amplification among spliced
sequences. A junction-spanning primer can therefore pass a genomic screen
without a contiguous intended genomic product. GENtle must identify an intended
genomic product by prepared-FASTA subject plus genomic binding geometry, never
by equality to the expected cDNA amplicon length.

Specificity searches must not use BLAST `-max_target_seqs` as a substitute for
a post-search hit-review threshold, because doing so can silently exclude
contigs. Handoffs and finalization bind the prepared BLAST database's index
kind and deterministic content fingerprint in addition to its prefix and genome
label. Replacing database content at the same prefix invalidates a pending
handoff. Whole-prepare activity state does not override successful
component-level `blastdbcmd` validation.

## DEC-030: Collection Lifting Is Subject-Specific And Engine-Owned

Status: active

Collection-visible capabilities declare their lifting behavior in the typed
`gentle.collection_lift_policy_registry.v1` registry. Policies are keyed by
capability and subject kind because a logical gene set, physical container, and
ordered arrangement do not have interchangeable semantics. Supported modes are
`map`, `combine`, `compare`, `arrange`, and `derive`; incompatibility is an
explicit typed rejection rather than an adapter-local conditional.

Dynamic readiness failures remain operation/fact-graph results and must not be
misclassified as permanent policy rejection. Generic collection reports carry
typed member outcomes/errors and links to produced domain reports. Membership
locks are set-like except for explicitly ordered subjects such as
arrangements, whose numeric member order is part of the fingerprint.

Collection members bind biological interpretation through a report-owned
context registry. Context-sensitive policies are fail-closed: an undeclared
policy is `not_reviewed`, while reviewed consumers may require one homogeneous
context and reject missing, conflicting, or target-mismatched contexts before
coordinate lookup, sequence search, or other biological work. Generic portable
collection reports copy the registry because their source may otherwise be
available only by report id; domain reports embedding the source resolution do
not create a second registry. The canonical collection membership fingerprint
continues to lock membership/order only and must not be interpreted as a
context or biological-input fingerprint.

## DEC-031: Ortholog Context Binding Extends The Specialized Resource

Status: active

The existing offline `gentle.ortholog_resource.v1` remains the authoritative
ortholog mapping contract. Optional source/target context ids bind each
directional endpoint to the resource's biological-context registry, and
resolved cohort reports copy referenced contexts so they remain portable.
Explicit organism or genome conflicts fail before genome-catalog work.

Orthology type and confidence are open typed string vocabularies: GENtle
recognizes canonical cardinality/confidence values while preserving unknown
provider-specific text exactly. Legacy symbol-only rows remain accepted, but a
matching symbol is never treated as relationship or functional-equivalence
evidence. The shipped `reject` and deterministic warning-bearing `first`
ambiguity policies remain unchanged. The additive `preserve` policy leaves an
ambiguous target unresolved while carrying every ordered candidate mapping and
its portable context/provenance; candidates do not become resolved cohort
members.

Unlike provider-owned orthology type and confidence text, ambiguity policy is
a closed GENtle operation-control vocabulary. The v1 report schema is retained
because `preserve` is opt-in and does not change the meaning or wire shape of
existing `reject` and `first` reports. Readers must reject an unknown policy
rather than silently defaulting it; consequently, a pre-`preserve` v1 reader
cannot consume a report that explicitly requested `preserve`. Candidate ranks
and labels are presentation references within one report, not durable
cross-run identities.

A generic relation hierarchy or paralog contract still requires a concrete
resolver/consumer and is not inferred from the ortholog contract.
