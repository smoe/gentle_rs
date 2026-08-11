# Extending GENtle With a Biological Adjustment

GENtle is intended to be a vehicle for turning a precise biological request
into a reproducible operation, not merely a collection of finished menu items.
This guide explains how to decide whether a request is already configurable,
can be composed from existing GENtle capabilities, or needs one new shared
engine operation.

The central rule is:

> Reuse biological interpretation and calculation inside the engine. Add one
> versioned request/report only when the new biological decision needs a stable
> identity of its own. Keep GUI, CLI, MCP, JavaScript, Lua, Python, and agent
> adapters as thin ways to invoke and inspect that same decision.

This document is for contributors and coding agents. GENtle's inner Agent
Assistant can classify a request, inspect current capabilities/readiness, and
prepare an implementation brief. It cannot edit the source tree and must not
invent a command for a capability that does not exist.

## 1. Classify Before Implementing

Every requested adjustment starts in one of three classes.

### A. Existing parameter

The operation already expresses the biological decision and only needs a
different value, constraint, target, or rendering choice.

Examples:

- use an exact primer length rather than a range
- move a search window closer to a feature boundary
- retain more ranked candidates
- change a display filter without changing the underlying report

Use `capabilities`, `introspect capabilities`, `introspect readiness`, the
glossary, and the protocol documentation before changing code. Do not add a
second operation merely to provide a convenient preset.

### B. Composition of existing capabilities

GENtle already has the required biological interpretations and calculations,
but no operation yet joins them into the requested result.

Examples:

- resolve one annotated transcript, obtain its mature 5'-to-3' sequence, select
  windows in its terminal exon, attach an adapter, and score a multiplex pool
- derive candidate regions between local sequence anchors, attach scores, and
  apply a top-k or Pareto selector
- combine a genome anchor with imported track evidence before ranking local
  candidates

This normally deserves one small engine-owned orchestration operation and one
portable report. Reuse existing private helpers rather than reconstructing
transcripts, coordinates, scores, or features in an adapter.

### C. Missing biological primitive

GENtle lacks an interpretation or calculation needed by the request.

Examples include a genuinely new evidence model, a new coordinate transform,
or an assay rule that cannot be stated with current constraints. Add the
smallest deterministic primitive first, test its biological boundary, and then
compose it into an operation. Do not hide the primitive in a GUI callback or
one-off script.

## 2. Write the Biological Contract First

Before choosing files or types, describe the request under these headings:

1. **Objective**: what decision or artifact is needed?
2. **Biological subject**: sequence, transcript, exon, feature, collection,
   container, genome assembly, track, or report?
3. **Orientation and coordinates**: genomic ascending, feature strand,
   transcript 5'-to-3', local sequence, or genome anchored?
4. **Fixed requirements**: exact motifs, adapters, lengths, roles, evidence,
   or target priority that must not drift?
5. **Adjustable policy**: windows, thresholds, scoring, ranking, tie breaks,
   or presentation choices?
6. **Ambiguity policy**: when must GENtle ask instead of choosing?
7. **Output**: selected candidates, all evaluated candidates, explanation,
   features, sequence artifacts, or an external handoff?
8. **Provenance**: which inputs, versions, coordinates, and policies must make
   the result reproducible?
9. **Non-claims**: what does the result explicitly not establish?
10. **Edge cases**: reverse strand, multiple transcripts, ambiguous bases,
    missing annotations, empty results, stale evidence, and repeated execution?

This contract should be understandable by a wet-lab biologist before it is
translated into Rust types.

## 3. Find the Existing Composition Seams

Start from concepts, not frontend widgets. The following source locations are
the usual route through GENtle:

| Need | Preferred source seam |
| --- | --- |
| Stable request/result records | `src/engine/protocol.rs` |
| Public operation identity and history policy | `src/engine.rs` |
| Shared operation orchestration | `src/engine/ops/operation_handlers.rs` |
| Transcript groups, exons, splice geometry | `src/engine/analysis/feature_expert_ops.rs` |
| Primer metrics, complementarity, Primer3 integration | `src/engine/state/sequence_ops.rs` and operation helpers |
| Genome assemblies, extraction, BLAST | `src/genomes.rs` and engine genome operations |
| Candidate sets, scores, filters, optimization | engine candidate/score operations documented in `docs/protocol.md` |
| Shared shell parsing/execution | `src/engine_shell/command_parsers.rs` and `src/engine_shell.rs` |
| GUI presentation/background invocation | relevant `src/app*` or `src/main_area_dna*` module |
| Runtime discovery and readiness | capability descriptors in `src/engine_shell.rs` plus project facts |
| User-facing command vocabulary | `docs/glossary.json` |

Useful private transcript/oligo seams currently include:

- `build_splicing_expert_view`: one engine interpretation of transcript groups,
  exon order, strand, and splice evidence
- `build_qpcr_transcript_design_templates`: mature-transcript sequence plus
  transcript-local/source exon geometry
- `resolve_single_transcript_design_template`: explicit isoform resolution that
  refuses ambiguous groups
- `enumerate_canonical_transcript_segment_windows`: fixed-length canonical
  windows in transcript-oriented coordinates, without embedding assay policy
- `map_transcript_local_interval`: strand-aware projection back to source DNA
- `compute_primer_heuristic_metrics`,
  `compute_primer_self_3prime_complementary_run`, and
  `compute_primer_pair_dimer_metrics`: shared oligo diagnostics

These are source-level extension seams, not public APIs. They may be renamed or
refined while the public operation/report remains compatible. Their private
Rust documentation can be inspected with:

```sh
cargo doc --document-private-items --no-deps
```
If a related request needs the same biological meaning, reuse or generalize a
seam narrowly. Do not copy its implementation into another handler.

## 4. Promote Only the Stable Decision

When a composition has a meaningful reusable output, give the decision a
versioned request and report:

1. Add strict request/report records to `src/engine/protocol.rs`.
2. Add one `Operation` variant in `src/engine.rs`.
3. Implement it in the shared operation handler.
4. Decide whether it is read-only, display-only, or persisted/undoable.
5. Record input digests, resolved identities, coordinate basis, policy, and
   warnings in the report.
6. Store the report in engine/project state when later inspection matters.
7. Add explicit export as a read-only operation when users need a file.

The public report should describe biological meaning, not private helper calls.
For example, report "terminal exon in transcript-oriented coordinates" rather
than "the result of `local_exon_segments.last()`".

## 5. Project One Operation Onto Every Interface

After the engine behavior is deterministic:

- shared shell/CLI parses input and prints the report
- GUI collects inputs, starts the operation in the background, and renders the
  report
- JavaScript and Lua call the shared shell or operation adapter
- Python calls `GentleClient.op(...)`
- MCP uses a dedicated stable tool or the generic typed `op` route
- the inner agent proposes only parser-valid commands and respects readiness
- ClawBio/OpenClaw wrappers remain deterministic wrappers around the same route

Prominent buttons or dedicated MCP tools are optional. Reachability through the
shared operation is not optional. No adapter may infer a transcript, transform
coordinates, rerank candidates, or reinterpret evidence differently.

## 6. Make the New Capability Explain Itself

A completed biological extension should leave enough information for a human,
small local model, or external coding agent to understand it without reading a
large dispatcher:

- capability description says what biological decision is made
- arguments identify subjects and coordinate conventions
- readiness lists required project facts
- expected effects identify the persisted report/artifact
- protocol docs define schemas and non-claims
- GUI/CLI docs show one exact use path
- glossary maps the command to its engine operation and reachable interfaces
- changelog records the completed behavior
- source comments explain non-obvious private invariants and composition seams

Introspection describes stable executable capabilities. It should not expose a
catalog of private Rust function names: that would turn implementation details
into an accidental wire contract.

## 7. Deterministic Validation Checklist

At minimum, test:

- request validation and unknown-field rejection
- deterministic selection and tie breaking
- correct forward- and reverse-strand coordinate mapping
- explicit rejection of unresolved biological ambiguity
- canonical/ambiguous sequence behavior
- persistence, export, undo/redo classification, and repeated execution
- shell parser/executor parity
- capability/glossary/MCP reachability parity
- GUI state parsing separately from biology calculations
- explicit non-claims in report/docs when a later evidence step is required

Use synthetic fixtures for edge geometry unless a real fixture is necessary.
Document fixture provenance as required by `AGENTS.md`.

## 8. Worked Example: Adapter-Bearing Terminal-Exon RT Primers

The request was: keep a fixed 5' adapter, derive an exact 22 bp
sequence-specific 3' segment near the start of each target's terminal exon,
respect reverse-strand transcripts, and reduce interference within an ordered
multiplex pool.

GENtle could already provide:

- annotated transcript grouping and explicit isoform identity
- mature-transcript sequence in biological 5'-to-3' order
- exon/source coordinate projection on both strands
- reverse complementation and primer diagnostics
- operation history, report persistence, exports, background GUI jobs, and
  adapter reachability

The missing piece was therefore class B: one composition operation, not a new
transcript parser and not a Python reimplementation. The resulting
`DesignTerminalExonRtPrimerPool` operation owns only the new decisions:

- target order is explicit biological priority
- candidate windows start from the transcript-oriented terminal-exon boundary
- the variable primer is the reverse complement of the target window
- the complete oligo is fixed adapter followed by variable segment
- adapter/self/prior-pool complementarity precedes distance in ranking
- Tm is descriptive only
- global transcriptome/genome specificity is not claimed

This is the preferred pattern for related requests: identify what GENtle
already knows, add only the missing biological decision, and preserve one
explainable result across every interface.

## 9. Prompt for a Coding Agent

```text
Objective:
Add or adapt <BIOLOGICAL WORKFLOW> in GENtle.

Repository:
<PATH TO gentle_rs>

Biological contract:
- subject and purpose: ...
- orientation/coordinates: ...
- fixed requirements: ...
- adjustable policy: ...
- ambiguity behavior: ...
- output and provenance: ...
- explicit non-claims: ...

First classify the request as:
1. existing parameter;
2. composition of existing engine services;
3. missing biological primitive.

Inspect docs/biological_extension_guide.md and the relevant protocol,
operation-handler, analysis, and adapter paths. Reuse existing private engine
seams where their biological meaning matches. Do not add biology to a GUI,
shell, MCP, Python, JavaScript, Lua, ClawBio, or OpenClaw adapter.

Implement the smallest shared-engine change, strict versioned records when a
new stable decision is required, provenance, parity, and deterministic edge
tests. State what remains unsupported rather than broadening the claim.
```
