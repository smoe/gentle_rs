# GENtle Introspection Contract (Proposal / Working Draft)

Last updated: 2026-06-26

Status: design proposal. This document is intentionally a starting point, not a
frozen spec. Sections marked "Degrees of freedom" call out where the
implementing agent (Codex) is expected to adjust shapes to its own view of the
code.

Revision note: incorporates the first Codex review (argument binding, fact
domains, effect modality, subroutes, narrowed "lossless", explicit host fact).
Where this draft keeps a different mechanism than the review literally proposed,
the divergence is called out inline as "Reconciliation".

Related shared documents:

- `docs/architecture.md`: single-engine invariant, adapter parity levels,
  `ui ...` intent plane, voice-path note, `data-gentle-role` SVG hooks.
- `docs/protocol.md`: `gentle.project_fact_graph.v1`, `gentle.fact_expression.v1`,
  `gentle.fact_evaluation.v1`, agent response schema.
- `docs/agent_interface.md`: agent-facing routes.
- `docs/glossary.json`: current source of truth for shell command semantics.

## 1. Why this exists

GENtle is moving toward a words-only interaction model: every operation,
including classically mouse-driven visualization controls, must be reachable and
inspectable through language alone. A words-only client (text agent, voice
transport, or accessibility frontend) must be able to ask, at any moment:

- **"What is true right now?"** — structured state read-back.
- **"What can I do right now?"** — capability/intent discovery, annotated with
  whether each action is currently ready.

Introspection answers both over **one controlled vocabulary** (the project fact
registry, `PROJECT_FACT_TYPE_SPECS`). The verbs (operations and intents) and the
nouns (facts) share that vocabulary so an operation's declared *effect* is stated
in the same language as another operation's *precondition*. That is what lets a
planner chain steps and check readiness deterministically instead of guessing.

### In scope

- Deterministic read-back of state as typed facts, **lossless for discrete,
  command-addressable state** (see section 7).
- A unified capability catalog where every verb self-describes its arguments,
  mutation/confirmation posture, and its preconditions/effects in fact
  vocabulary, **bound to its own arguments** (section 4).
- Per-capability readiness computed against the current fact graph, in both an
  unbound and a bound mode (section 4).

### Explicitly out of scope (downstream)

- **Interpretation of a visualization.** Narrating the *meaning* of a rendered
  map, gel, or dotplot is an interpretation layer built on top of introspection.
  Introspection reports facts ("selection is 100..200; linear view; three EcoRI
  sites shown"); it does not read pixels. Full blind/accessibility support
  depends on this layer and is deferred. See section 8.
- **Continuous/transient GUI state.** Live drag-in-progress, hover state,
  pixel-level layout, and transient focus quirks are not introspection
  obligations (section 7).

## 2. The two questions, one vocabulary

```
            verbs (capabilities)                 nouns (facts, by domain)
        +-------------------------+        +-----------------------+
        | operation               |        | project.* (sequence,  |
        | view_intent             |  reads |   report, restriction)|
        | host_config             | -----> | view.*                |
        |   precondition_expr ----+--------> host.*                |
        |     (bound to $args)    |        | config.*              |
        |   effects + effect_kind +--------> (post-run verification)|
        +-------------------------+        +-----------------------+
                     \                                /
                      \___ shared registry vocabulary _/
                           (PROJECT_FACT_TYPE_SPECS)
```

Both halves already partially exist. Introspection binds them and fills the gaps
without forking any execution path: it composes the existing `capabilities`,
`ui ...` discovery, `facts graph`, and the fact evaluator. **Composition, not
duplication** — it must not become a parallel command surface.

## 3. Routes and payload

**Reconciliation (review pt 6):** one root with several subroutes, not a single
scoped route.

- `introspect facts` — state read-back (the fact graph, grouped by domain).
- `introspect capabilities` — the verb catalog with self-descriptions.
- `introspect readiness` — capability readiness against current facts.
- `introspect all` — the aggregate, for clients that want one round trip.

Scoping applies within each subroute (filter by `domain`, `kind`, `readiness`,
or subject `--seq-id`), because the catalog is large (hundreds of shell
commands). A words-only client must be able to ask "what can I do with `insert`
right now?" and get a short, ready-annotated list.

`introspect capabilities` descriptor shape:

```json
{
  "schema": "gentle.introspection.v1",
  "capabilities": [
    {
      "id": "features restriction-scan",
      "kind": "operation",
      "mutating": false,
      "requires_confirmation": false,
      "args": [
        { "name": "SEQ_ID", "required": true, "detail": "loaded sequence id" },
        { "name": "--enzyme", "required": false, "detail": "narrow to one enzyme" }
      ],
      "reads":   [ { "fact": "sequence.exists", "subject": { "arg": "SEQ_ID" } } ],
      "effects": [ { "fact": "report.exists", "report_kind": "restriction_scan",
                     "subject": { "arg": "SEQ_ID" }, "effect_kind": "must_on_success" } ],
      "precondition_expr": { "all": [ { "fact": "sequence.exists", "subject": { "arg": "SEQ_ID" } } ] },
      "description": "Scan a loaded sequence for restriction-enzyme sites."
    },
    {
      "id": "ui select",
      "kind": "view_intent",
      "mutating": false,
      "requires_confirmation": false,
      "args": [
        { "name": "SEQ_ID", "required": true, "detail": "subject to select within" },
        { "name": "--range", "required": false, "detail": "START..END (0-based, end-exclusive)" }
      ],
      "reads":   [ { "fact": "sequence.exists", "subject": { "arg": "SEQ_ID" } },
                   { "fact": "ui.host_available" } ],
      "effects": [ { "fact": "view.selection", "subject": { "arg": "SEQ_ID" },
                     "effect_kind": "view_session" } ],
      "precondition_expr": { "all": [
        { "fact": "ui.host_available" },
        { "fact": "sequence.exists", "subject": { "arg": "SEQ_ID" } }
      ] },
      "description": "Select a region so later view/analysis intents can act on it."
    }
  ]
}
```

`introspect facts` payload groups by domain and stays
**deterministic and timestamp-free** (like `state-summary` / `facts graph`):

```json
{
  "schema": "gentle.introspection.v1",
  "facts": {
    "project": { "schema": "gentle.project_fact_graph.v1", "facts": [ "..." ] },
    "view":    { "host_attached": true,  "facts": [ "..." ] },
    "host":    { "facts": [ "..." ] },
    "config":  { "facts": [ "..." ] }
  }
}
```

## 4. Capability self-description, argument binding, and readiness

**Review pt 2 (adopted, with the consequence made explicit).** A fact name alone
is not enough: `sequence.exists` must bind to a command argument, or readiness
can only say "*some* sequence exists," never "this command is ready for *this*
sequence." So a capability's precondition/effect atoms are **templates** that
reference the capability's own arguments.

This splits fact atoms into two forms:

- **Ground atom** — a concrete subject, e.g.
  `{"fact":"sequence.exists","subject":{"kind":"sequence","id":"insert"}}`.
  This is what `facts graph` projects and what the evaluator already consumes.
- **Templated atom** — references an argument, e.g.
  `{"fact":"sequence.exists","subject":{"arg":"SEQ_ID"}}`. This appears only in
  capability descriptors and must be **instantiated** before evaluation.

Consequently there are two readiness modes:

- **Unbound (catalog) readiness** — no argument values supplied. Templated atoms
  resolve to `unknown` with reason "unbound argument"; only argument-free atoms
  (`ui.host_available`, host/view availability) resolve concretely. Answers
  "what *could* I do?" and is cheap to compute for the whole catalog.
- **Bound readiness** — caller supplies argument values; the engine substitutes
  each `{"arg":...}` with the concrete subject, then runs the existing evaluator.
  Answers "is this command ready for *these* inputs?" This is the EDITtoTrEMBL
  discipline: the precondition binds to the actual input, not to the catalog.

Instantiation is a deterministic pre-pass over `FactExpression`; the tri-state
evaluator and `FactTruth` are reused unchanged.

### Effect modality (review pt 5, adopted)

Not all effects are guaranteed. Each effect carries an `effect_kind`:

- `must_on_success` — the effect must hold if the command succeeds. **Only these
  are asserted by post-run verification** (diff the fact graph before/after).
- `may_on_success` — conditional/optional; not asserted.
- `external_handoff` — an external/provider action, not a local state change.
  Aligns with the existing `mutating: "external"` posture.
- `view_session` — view/session-only, never recorded in project provenance.
  Aligns with `kind: view_intent`.

### Capability kinds

- `operation` — engine operation; may mutate; may require confirmation; recorded
  in provenance/undo when mutating.
- `view_intent` — non-mutating presentation control (select, set view mode,
  zoom/fit, pan, toggle tracks). Routes through the `ui ...` plane, which already
  rejects `state_changed=true`. Not recorded as a project mutation.
- `host_config` — host/session configuration (external tool paths, language,
  graphics) through a shared host-config contract. Not project state.

Reachability is identical for all three kinds (section 6); they differ only in
routing and in whether the action is recorded.

### Definition of "introspection-complete"

A capability is introspection-complete when it declares `args` (reuse glossary
operand metadata), `mutating`/`requires_confirmation` (consistent with MCP
`tools/list`), and `reads[]`/`effects[]` as registry-vocabulary atoms bound to
its arguments. Unknown/forward facts remain legal and evaluate to `unknown`,
never an error.

## 5. Fact domains and new namespaces

**Review pt 3 (adopted as fact domains) — Reconciliation: do not rename the
landed facts.** The shipped `gentle.project_fact_graph.v1` facts
(`sequence.*`, `report.*`, `restriction_site.*`) stay exactly as they are. Add an
**additive `domain` discriminator** (default `project`) to the registry spec and
to projected facts, then group by domain in introspection. This gets the clean
separation the review asks for without breaking the schema, the generated prompt,
or its tests.

Domains:

- `project` — saved project state (the existing facts, unchanged).
- `view` — current GUI/session projection.
- `host` — local machine/tool availability.
- `config` — GENtle parameters/configuration.

Suggested additions:

| fact | domain | world | requires_basis | notes |
|---|---|---|---|---|
| `ui.host_available` | view | closed | no | **already exists**; reclassified into `view`. The canonical "is a GUI host attached" fact — see section below |
| `view.mode` | view | closed | no | linear / circular / cdna / protein; `equals` |
| `view.selection` | view | closed | no | subject + `range` |
| `view.viewport` | view | closed | no | visible bp window / lane window |
| `view.visible_tracks` | view | closed | no | which tracks/lanes are shown |
| `config.param.<name>` | config | closed | no | `SetParameter` values; `equals`/`compare` |
| `host.tool_available` | host | closed | no | e.g. primer3/blastn resolvable; see Reconciliation below |

These are closed-world: the host/adapter resolves them deterministically when
asked, so a missing fact means `unsatisfied`, not `unknown` (unlike open-world
proof facts such as `restriction_site.absent`).

### Headless behavior (review pt 4, adopted with reconciliation)

For CLI/MCP without a GUI, omitting `view.*` facts is ambiguous (no selection vs
no host). Fix: **always project the host-presence fact explicitly.**

**Reconciliation:** rather than add `view.host_attached` *and* `view.available`
as new spellings, reuse the existing `ui.host_available` as the single
host-presence fact (now in the `view` domain). In headless contexts it projects
`false` and all other `view.*` facts are omitted by design. Then bound readiness
for `ui select` resolves to `unsatisfied` with `ui.host_available` as the named
unmet atom — surfaced to the user as "unavailable: no view host," **not** as a
new fourth `FactTruth` value. The three-valued evaluator is preserved; the
"unavailable" distinction is a presentation of unsatisfied-with-reason.

### Reconciliation: where does `host.tool_available` live?

Tool availability already has adjacent contracts (`agents preflight`,
helper/genome `status`, primer3 preflight). `host.tool_available` should
**project from those existing probes**, not introduce a parallel detector. Open
for Codex: whether it belongs in the introspection `host` domain or in a separate
host-capability projection that introspection merely references.

### Degrees of freedom

Namespace spelling (`config.param.foo` vs `param.foo`), whether `view.*` facts
are projected lazily (only with a host attached) vs always with a marker, and
whether `domain` is a field on each fact vs a wrapper per group.

## 6. Reachability and parity

Success criterion for the words-only model is twofold:

1. **No control without a word path.** Every state-mutating control and every
   visualization/mouse-driven control resolves to some agent-reachable contract
   (operation, view_intent, host_config). A GUI-only mutation path with no shared
   route is a reachability-parity gap (a bug). Continuous gestures are not
   commands, but their discrete forms are (`ui select`, `ui zoom-to`,
   `ui fit-features`), and those must exist.
2. **No capability without a self-description** (section 4).

First-class surfacing stays per-adapter and may differ; reachability and
introspection-completeness are strict.

## 7. "Lossless", narrowed (review pt 7)

Read-back is **lossless for discrete, command-addressable state** only:
selection range, view mode, visible tracks, parameter values, loaded sequences,
reports. It is explicitly **not** an obligation for continuous GUI gestures,
hover state, drag-in-progress, transient focus, or pixel-level layout. Those
either have a discrete command form (which is what gets introspected) or are
presentation-only and out of scope.

## 8. The accessibility boundary

Introspection is the structured substrate for blind/voice support but stops at
facts. The bridge to the downstream interpretation layer already exists: the
engine-owned SVG exports expose `data-gentle-role` / `data-gentle-feature-kind`
on semantic glyphs.

- **Now (introspection):** `view.*` facts + capability catalog answer "what am I
  looking at and what can I do?" with no interpretation.
- **Downstream (accessibility):** narrating the meaning of a rendered figure
  consumes the same `data-gentle-role` structure plus the fact graph. Separate
  capability, out of scope here.

## 9. Constraints

- **Single engine.** No business/biology logic in adapters; introspection
  composes engine-owned contracts shared across GUI/CLI/MCP/JS/Lua.
- **Determinism.** Same state -> byte-identical payload; no volatile timestamps.
- **Safety.** Reachable is not auto-executed. Mutating capabilities keep their
  confirmation gate; view_intents may run freely.
- **Reuse, don't fork.** `facts`, vocabulary, `readiness`, and capability
  metadata come from existing sources.

## 10. Suggested first slice (review tweak adopted)

1. Add the `domain` discriminator and the `view` domain, including
   `ui.host_available` projected as `false` in headless contexts (read side
   first — lowest risk, immediately useful to a words-only client).
2. Add `introspect facts` and `introspect capabilities` as separate subroutes.
3. Attach engine-owned bound self-description (`args`, `reads`, `effects` with
   `effect_kind`, `precondition_expr`) to **three** capabilities only:
   one mutating operation, one `ui` view intent, one host/config capability.
4. Add readiness tests proving the **same descriptor** resolves correctly in both
   a GUI-attached and a headless context (the binding + host-domain payoff).
5. Defer: full catalog annotation, `introspect readiness`/`all`, post-run effect
   verification, glossary-as-projection inversion, and any
   visualization-interpretation work.

## 11. Resolved decisions (pending first-slice validation)

These were open questions in the prior draft; the review resolved them. Promote
the non-negotiables into `docs/decisions.md` **only after** the first slice
validates the shape.

- **Routes:** subroutes + aggregate, not one scoped route (section 3).
- **Domains:** additive `domain` field; landed facts unchanged (section 5).
- **Headless:** explicit `ui.host_available` fact; "unavailable" is
  unsatisfied-with-reason, not a new truth value (section 5).
- **Effects:** carry `effect_kind`; only `must_on_success` is verified
  (section 4).
- **Args source (long term):** an engine/protocol-side descriptor becomes the
  root truth and `docs/glossary.json` becomes a projection of it. Direction
  accepted; the inversion is **deferred** past the first slice.

## 12. Still open for Codex

- Exact representation of `domain` (per-fact field vs per-group wrapper).
- Whether `host.tool_available` lives in the introspection `host` domain or a
  separate host-capability projection it references.
- Whether bound readiness should be exposed as part of `introspect capabilities`
  (when args are supplied) or only via `introspect readiness`.
