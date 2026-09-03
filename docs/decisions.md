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

Project policy permits only user-approved, window-bound screenshots for Agent
Assistant help. A direct request may originate from an `Agent help` control in
the GENtle viewport being diagnosed. An agent response may alternatively carry
one typed, path-free request id and human-readable reason. That response emits
no capture command: GENtle shows a consent card, the user selects one currently
registered GENtle content viewport, and `Allow one screenshot` consumes a
turn/provider/project/viewport-bound approval before requesting one egui
capture. The agent cannot choose an arbitrary target or grant approval.

Both structured project-state context and screenshots disclose information to
the selected provider; this policy does not assume that a screenshot is more
sensitive. A project summary can be scientifically richer, while a visual
snapshot is often more useful for layout, direction, colour, and enabled-state
guidance. They remain separate review controls because their scope and
diagnostic purpose differ. Captured images are previewed locally and sent only
when the user explicitly submits the follow-up prompt. Conversation history
stores path-free request/attachment provenance, not live approval, image bytes,
or the temporary local image path.

Native macOS full-window capture is an optional same-process enhancement. It
must use ScreenCaptureKit, refuse windows not owned by GENtle, and remain
unavailable when the user does not grant macOS Screen Recording permission.
Automatic/background capture, desktop capture, capture of another application,
agent-selected capture, and headless capture remain prohibited. The reserved
`screenshot-window` CLI/shared-shell route therefore stays disabled;
deterministic SVG/export routes remain preferred for reproducible artifacts.

## DEC-015: ClawBio/OpenClaw Boundary

Status: active

The ClawBio/OpenClaw skill scaffold wraps deterministic `gentle_cli` routes and
writes reproducibility bundles. It should surface shared capabilities,
suggested actions, and UI-intent handoffs without creating ClawBio-only biology
logic.

Every wrapper run carries a topic-neutral, content-bound provider execution
receipt. The receipt binds the resolved request, immutable inputs, mutable
state before and after execution, wrapper/runtime identity, command outcomes,
native result identifiers/status fields, and emitted artifacts. It is not a
discussion ledger, permission record, participant claim, freshness model, or
scientific interpretation. External systems may retain the receipt, but own
their own request identity, policy, and interpretation records.

When a caller needs a policy-aware analysis-run record, the provider receipt
and caller-owned run remain two records joined by content hashes rather than
one shared mutable object. A post-run skill adapter may prepare a typed intake
packet, but it must not discover or mutate the caller's store. The caller must
revalidate current references and permissions, recompute the request digest,
and assign replay, dependency, and freshness state through its owning commands.

Descriptor-only skills may delegate to the generic runtime only through an
explicit versioned contract whose co-deployed catalog, descriptor, selected
route, and plan step verify before execution. This records the route that was
selected; it does not claim deterministic replay of natural-language routing.
Wrapper execution outcome and native biological status remain separate, so a
successfully computed scientific `fail` is neither an execution failure nor a
pass.

When a ClawBio presentation combines scientific artifacts, the raw GENtle
result remains authoritative. Strict-attribution presentations must label
executable-derived statements separately from wrapper presentation, caller
inputs, named external tools, and unverified prose. Every executable-derived
summary requires an exact native-result pointer or scope, content hash, and
named deterministic projection in a machine-readable claim/artifact ledger.
Wrapper-authored graphics may arrange GENtle artifacts but must not recalculate
or upgrade their scientific content.

Adapter execution timeouts are approval-bound operational ceilings, not
scientific parameters or expected runtimes. Long primer-design routes must
allow the ceiling to be reviewed explicitly; adapters must not present a fixed
wrapper timeout as a Primer3 limitation. Runtime expectations should consider
transcript-equivalence classes, target ROIs, requested junctions, assay mode,
and operation count in addition to sequence span.

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
typed member outcomes/errors and links to produced domain reports. When a
domain report has no persisted store/getter, a portable domain wrapper may own
the child reports directly; in that case `produced_report_ids` stays empty so
GENtle never emits dangling report references. Membership locks are set-like
except for explicitly ordered subjects such as arrangements, whose numeric
member order is part of the fingerprint.

Collection wrappers over non-persisted child scans own those reports directly
and must describe aggregate completeness. A retained-row cap is not evidence
that later motifs or candidates were scanned: aggregate count fields must be
named as retained/observed values and accompanied by an explicit completeness
flag plus affected member ids or warnings.

Versioned collection result wrappers are additive output contracts: they
default missing fields and tolerate unknown future fields. Collection request
and member-binding records remain strict and may reject unknown fields so a
misspelled mutation or binding parameter cannot be silently ignored.

A mutating collection lift must expose a non-mutating preview before apply.
Its apply lock binds the exact source snapshots, effective parameters, and the
deterministically reserved output namespace in addition to the collection
membership fingerprint; the membership fingerprint alone is not a mutation
guard. Apply must reject a changed plan before inserting any output. Derived
sequence lineage is recorded from each source member to its own products, and
logical collection products do not imply that those products occupy one
physical container.

Shared digest convergence guards must be deterministic. Direct and lifted
digests use fragment-count, round-count, repeated-state, and no-progress
limits; wall-clock load must not decide whether an otherwise identical digest
succeeds.

A collection `combine` that writes one external artifact is atomic rather than
partially successful: all members validate before the single write, and a
successful report links every member to the same wrapper report. Pool export
requires an explicit physical container whose declared contents are exclusive.
In-silico `selection` containers are not physical pool sources and are rejected
before export; physical `singleton` and `pool` containers remain eligible.
`gentle.pool.v1` has no incompleteness/provenance slot, so exporting a
non-exclusive known subset is rejected before writing rather than silently
losing that meaning at the artifact boundary.

Gel-render lift modes follow the resulting lane geometry. An explicit
container is a `combine` into one physical lane; a serial arrangement is an
`arrange` that preserves its declared lane and role-label order. A bare project
sequence selection is not a valid collection lift because the direct input
form co-migrates every supplied sequence in one lane without proving that the
sequences share a physical sample. `ArrangementMode::Plate` incompatibility is
a dynamic readiness/execution error, not a permanent arrangement-subject
policy rejection. The shared `render-pool-gel-svg --arrangement` route remains
executable even while the generic collection launcher reports that it has no
typed arrangement adapter.

Explicit gene-set pooling is a physical aliquoting operation, not an inference
from logical membership. Every resolved member must bind to a distinct,
exclusive singleton source container, and incomplete resolutions reject the
whole plan. Apply creates one lineage-linked derived aliquot sequence per
source, places only those aliquots in one new exclusive `pool` container, and
retains the source containers and their latest-membership mappings. The plan
lock covers the resolution membership, source-container bindings, exact source
sequence snapshots, output ids, and reserved pool-container id.

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

## DEC-032: Conversational Scientific Delegation Is Proposal-Then-Execute

Status: active

Natural-language routing may draft a structured GENtle request, but it is not
execution authority for a delegated route that selects biological material, a
pair or backend, mutates project state, or writes/replaces an artifact. Such a
route first emits a content-addressed proposal binding the normalized request,
exact command, delegated descriptor/route, runtime, state, scientific inputs,
resolved paths, reference context, and pinned backend. Execution consumes that
stored request only after a caller-supplied approval assertion names the exact
proposal digest; it must not reconstruct slots from prose a second time.
For this gated path, known local Cargo launchers are reduced to the already
built, hash-bound `gentle_cli` executable before the proposal is issued, and
OCI execution requires an immutable image digest. Approval must never authorize
a later source rebuild or mutable image-tag resolution.

GENtle verifies the digest and fails on bound-context drift. The caller or
OpenClaw control plane owns approver identity, authentication, authorization,
and proof that a human actually reviewed the proposal. Read-only list, show,
and preflight routes may remain automatic. Direct structured GENtle and
`gentle_cli` use remains independent of this conversational delegation gate.
Material cache/tool environment overrides are part of the approved context;
GENtle does not claim to bind arbitrary unlisted operating-system state.

## DEC-033: Isoform-Assay Planning Has Two Scientific Approval Digests

Status: active

Conversational gene isoform-assay planning separates recommendation from
execution. The first approval binds the complete normalized planner request,
including policy schema/version/hash, evidence report ids and hashes, prior
plan, observations, override, and every effective default. It authorizes
planning from that exact basis; it does not validate the biology. The second
approval binds the exact ordered `DesignTranscriptAssayPanel` payloads emitted
by the approved plan. Those payloads are not regenerated or altered after
approval, and one batch approval is valid only because its digest covers the
complete ordered set.

An optional pre-approved informative fallback is part of both approval bases,
not a runtime convenience flag. The normalized request digest binds its typed
policy and whole-panel budget. The emitted operation digest binds the exact
wrapper, unchanged strict child, deterministic fallback template, and distinct
report identity. Execution may derive that child only after a matching typed
strict coverage-infeasibility record; all other failures remain fail-closed.
The strict failure and partial fallback are separate durable records, and a
bounded-greedy selector certificate must never be presented as a proven global
optimum or as experimental/order readiness.

For a multi-study second stage, independently approved plan/workflow pairs are
first composed into one engine-owned batch record. That record binds each
plan/workflow file, their existing workflow and operation digests, their order,
and the combined ordered operation set. Execution revalidates the complete
batch before applying its first operation and runs it through one state-bound
command. An adapter must not create several execution proposals against one
pre-mutation project state and then run them serially; the first mutation would
make the remaining approvals stale.

An oligo order form may be projected from this work only through a specific
`gentle.experimental_assay_handoff.v1` whose bytes and named readiness policy
are retained. Every order line keeps the handoff hash, readiness-policy id and
hash, readiness row/card, and assay/pair/oligo identities. Presentation
profiles and content-block selection may project an already immutable
canonical report without another scientific approval; they still produce an
execution receipt and remain subject to ordinary file-overwrite safety. An
outer agent may host or shorten those declared projections, but cannot author
new scientific blocks.

## DEC-034: Approved Batch Reuse Is Explicit, Exact-Prefix, and Detached

Status: active

Long-running approved isoform-assay batches execute against a detached engine;
by default the live project is committed only after the complete ordered batch
succeeds. DEC-035 adds an explicit opt-in that lowers this unit from the batch
to the gene without removing it.
An explicitly configured private checkpoint directory may retain the detached
state and exact operation results after each successful approved operation.
These artifacts survive a later rollback, but they are not project mutations
and are never discovered or imported during ordinary execution.

Reuse is a separate read-only proposal followed by explicit approval of that
proposal digest. The proposal binds the target batch, baseline project hash,
exact approved-operation prefix, checkpoint manifest and engine bytes, and the
running GENtle executable identity. Execution revalidates every binding and
imports only into another detached run. Adding work after an unchanged prefix
is supported; arbitrary branch-state merging, cross-build reuse, payload drift,
and approval-based waivers of provenance mismatches are not. Approval
authorizes transfer of prior computation and does not validate the biology.

## DEC-035: Per-Gene Continuation Lowers Atomicity, It Does Not Remove It

Status: active

`--on-gene-failure continue` is opt-in, and the default stays `abort` so an
existing command keeps its exact behaviour and output. Under `continue` the
unit of atomicity moves from the batch to the gene: a failing gene is restored
to an in-memory boundary snapshot taken when that gene started, so a gene is
never committed half-applied, and the remaining precomputed genes still run. A
batch in which no gene completed commits nothing.

This rollback guarantee is limited to engine state. Therefore `continue`
refuses the complete batch before its first operation if any approved operation
declares an external output path; a file already written cannot be made atomic
by restoring an in-memory engine snapshot. Callers must omit per-operation
paths or retain the default batch-level `abort` policy. The boundary snapshot
also includes the reported checkpoint pointer, so a failed gene cannot leave
`last_checkpoint_manifest` naming a partial-gene checkpoint that was rolled
back. More strongly, `continue` persists checkpoints only at complete-gene
boundaries; it never writes a partial-gene checkpoint that could remain
discoverable after rollback. The default `abort` policy intentionally retains
per-operation checkpoints because those are its recovery units.

The flag governs execution failures only. Approval and verification failures
still refuse the whole batch under either policy, because they say the request
was not the approved one rather than that a result was not obtained. The
requested coverage policy is likewise untouched: coverage is a biological
statement about one assay, this flag is an execution statement about a batch,
and a gene failing `require-all` is skipped rather than downgraded. A report
from a run with failures sets `batch_complete: false` and never presents itself
as a finished study.

Checkpoints keep describing an exact contiguous ordered prefix. Rather than
teaching the checkpoint schema about gaps — which would cost the invariant that
makes reuse trustworthy — checkpoint writing freezes at the last contiguous
prefix once a gene has been skipped. Because such a run also commits the genes
that succeeded, its own frozen proposal is correctly refused afterwards, and
recovery is a new batch of only the failed genes against the committed
baseline. For the same reason `continue` refuses to resume from a reuse prefix
that ends inside a gene, where no boundary snapshot can exist.

## DEC-037: Primer3 Is Selected By Capability, Then Version, Then Recency

Status: active

Where several Primer3 builds are installed, taking the first `PATH` match makes
GENtle's behaviour depend on shell configuration. A bare configured name is
therefore resolved by ranking every candidate: a build advertising `--progress`
beats one that does not, because bounded search reporting and SIGUSR1 status
are built on it and a legacy build silently degrades them; then the highest
version, compared component by component rather than lexically; then the newest
file, which resolves the remaining ties without reintroducing `PATH` order.

A configured value naming a path is authoritative and is never overridden by
discovery, because an operator pinning a build has made exactly the decision
this ranking would otherwise take back. The selection is reported rather than
implicit: preflight lists every candidate with its capability, version and
creation time, and states why the winner won, so a result can be traced to the
build that produced it. It is also cached per candidate set and file identity,
because checkpoint reuse binds the running executable and must not see a
different choice mid-run.

## DEC-038: Collection Construct Inspection Does Not Persist Temporary Graphs

Status: active

Collection construct-reasoning inspection is a `map` over project sequences or
homogeneous gene-set resolutions. It builds each member's current reasoning
graph on an engine clone and returns graph identity, input fingerprint, and
inspection actions in the collection wrapper, without writing those temporary
graphs into project metadata. Consequently `graph_persisted` is false and the
generic member status does not advertise a `produced_report_id` that cannot be
retrieved later.

Project-sequence inspection is context-agnostic because the operation has no
reference-genome parameter; a gene-set resolution must nevertheless be
homogeneous so one report does not silently mix declared biological contexts.
Containers and arrangements are explicitly rejected: physical occupancy and
member order do not change independent per-sequence reasoning, and accepting
them would imply unsupported physical semantics.

## DEC-036: A Dossier States What Is Missing Rather Than Waiting For It

Status: active

Experimental work starts from whichever genes are ready, so an isoform-assay
dossier is publishable while a study is incomplete. Each gene carries a
`resolved`, `pending`, or `unresolved` status; `unresolved` requires a reason
stating why the established automatism could not address that gene. The index,
each incomplete gene page, and the print document say so visibly, and the
canonical report carries per-status counts and a `complete` flag.

An omitted status is derived from the gene's own contents rather than assumed
to be `resolved`, so a request written before this field existed reports
honestly instead of rendering a gene without results as finished. A declared
status is validated against those contents for the same reason: the label is
only worth having if a reader can rely on it, which is why `resolved` without
any assay handoff is an error rather than a silently empty section.

## DEC-038: Primer Similarity Guides Search But Does Not Replace Specificity

Status: active

Background: Ye et al., *Primer-BLAST: A tool to design target-specific primers
for polymerase chain reaction* (BMC Bioinformatics 2012),
https://doi.org/10.1186/1471-2105-13-134.

A broad cDNA or genomic BLAST search may identify loci and intervals where a
primer is likely to cross-amplify, and GENtle may use that evidence to exclude
or deprioritize Primer3 search windows. The planning artifact must bind the
database content, query templates, coordinates, and exact treatment of each
interval. Similarity alone never proves paralogy, target identity, or assay
specificity.

Final specificity remains a separate complete target-space gate. BLAST finds
candidate loci; an explicit full-primer alignment may adjudicate their complete
binding geometry and terminal differences. A reviewed off-target allowance is
exact and evidence-bound, never a family-wide shortcut and never a way to turn
an incomplete search into a pass. Replacement primers proposed after failure
must re-enter the same cDNA and genomic gates.

For groups of caller-supplied related sequences, exact common primer placement
is derived before Primer3 from recorded global alignments. The longest member
is merely the default representative; it is not declared canonical. Alignment
and candidate budgets are explicit and content-bound, and `require_all` cannot
be weakened by runtime fallback. The resulting full pair-by-member product
matrix remains distinct from later genomic and transcriptome specificity.

Primer reports may additionally project accepted candidates onto a bounded
Pareto frontier. This is an explanation of non-dominated tradeoffs after hard
constraints, not a second acceptance rule. The declared design objective and
final specificity/readiness gates remain authoritative.

## DEC-039: External Auditor Owns The Performance Verdict

Status: active

GENtle developers provide deterministic benchmark targets, fixture provenance,
source-revision binding, and compile/smoke evidence. They do not certify the
performance of their own implementation from a development machine. A named
external auditor runs timed benchmarks on a stable host, retains the raw
Criterion artifacts and environment metadata, compares only compatible
baselines, and owns the release-facing performance verdict.

Criterion is a CPU-side characterization tool, not proof that the GUI is
responsive. GUI acceptance also requires the real application, representative
project state, window creation, repaint and resize behavior, and profiler
evidence where needed. Shared CI runners may compile or smoke benchmark targets
but do not establish strict timing gates without repeated stable-host evidence
and an auditor-approved threshold.

Routine thin-LTO, unwind-panic audit binaries and exact release-like fat-LTO,
abort-panic audit binaries are distinct evidence classes. Baselines are
comparable only when profile, host, toolchain, fixture hash, and GENtle revision
match; retained metadata and baseline labels must identify the profile
explicitly.

## DEC-040: Tutorial GUI Contracts Cannot Invent Execution Authority

Status: active

Tutorial GUI acceptance contracts may select only controls in a closed shared
semantic-control catalog. Each catalog row fixes the owning window, allowed
interaction, whether the action persists project state, whether it may
establish a scientific result, and any replacement-text policy. Merely knowing
a semantic snapshot id is not authority to type into or activate that control;
command-capable and snapshot-only controls are unavailable unless deliberately
added with a narrow typed policy.

Project persistence and scientific proof are independent. Metadata edits may
make the project dirty without claiming a scientific result. A scientific
effect requires before/after facts plus an engine/report/artifact/state
verifier; visual state never substitutes for that proof. GENtle validates the
contract and exposes read-only semantic state, while an isolated external
runner owns input delivery, evidence collection, and the final acceptance
verdict.

## DEC-041: Precomputed Genome Motif Evidence Is Optional And Scale-Separate

Status: active

Finalized sparse whole-genome motif packages are optional evidence providers,
not GENtle runtime resources. GENtle discovers and invokes an external DuckDB
executable only for an explicit query, opens only payload paths selected from
the package's exact inventory, and imposes motif, interval, file, row, output,
and wall-clock bounds. An absent package or executable therefore has no startup
cost and cannot disable local JASPAR scoring or any unrelated workflow.

Package hits remain on the package-declared score and retention policy. They
must not share a numerical axis with GENtle's local scorer unless a future
contract proves cross-source calibration. Every report states per-motif
retention-floor completeness, density limiting, row truncation, and genome
compatibility. Contig aliases and geometry alone may support a
`contig_geometry_matched_only` result; they do not establish per-contig
sequence identity. Row queries require an explicit bounded motif set because
the sparse package is not an all-motif density index.

## DEC-042: SCREEN cCRE Evidence Is Optional, Assembly-Bound, And Non-Causal

Status: active

ENCODE SCREEN Registry files are optional external evidence, not bundled
project data and not startup dependencies. GENtle accesses them only through an
explicit install, index, overlap-query, or materialization request. An absent
BED/index therefore cannot slow startup or disable local sequence, annotation,
TFBS, promoter, or cloning work.

Every supported source descriptor fixes provider release, species, taxon,
assembly aliases, subset, expected classes, URL, and BED field order. The
expected accession namespace is validated while indexing. The derived sparse
index binds the exact BED SHA-256 and byte size, and every query checks both
that content identity and the loaded sequence's genome anchor. Human and mouse
resources are never inferred from filenames or silently mixed.

A cCRE/ELS overlap is biochemical candidate-regulatory evidence. It does not
prove cell-type activity, identify a target gene, or establish causal
regulation. Materialized rows remain ordinary strandless
`regulatory_region` features carrying the source identity and these limits, so
all adapters and renderers consume the same engine-owned interpretation.

## DEC-043: Regulatory Annotation, Activity, And Primary Signal Stay Distinct

Status: active

Ensembl Regulation resources are optional external evidence and never startup
dependencies. Annotation snapshots are pinned by provider release, API and
pipeline versions, gene-annotation release, species, taxon, assembly name and
accession, then normalized and indexed under exact content hashes. Queries fail
closed on release or assembly mismatch; notably Ensembl GRCm39/mm39 and SCREEN
GRCm38/mm10 mouse resources are not interchangeable.

A regulatory-feature overlap, an associated-gene annotation, a regulatory
activity call, a primary peak, and a quantitative signal are separate claims.
GENtle must retain their independent source identities and cannot derive a
biosample-activity or causal target-gene claim from annotation overlap alone.
Future primary-signal display should use bounded region reads, explicit assay,
biosample, condition, replicate, output-type and unit provenance, and a
declared derived-activity rule rather than silently reducing signal to peaks.
