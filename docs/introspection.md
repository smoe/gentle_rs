# GENtle Introspection Contract (Proposal / Working Draft)

Last updated: 2026-06-27

Status: implemented for `introspect facts`, `introspect capabilities`,
`introspect readiness`, and `introspect all`. The capability route now projects
the shared protocol capability registry and overlays fact-aware annotations for
the validated first slice. Full fact annotation for every catalog row remains a
working design. Sections marked "Degrees of freedom" call out where later
implementation should adjust shapes to validated code.

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
- `introspect verify-effects` — post-run verification for hard effects.
- `introspect all` — the aggregate, for clients that want one round trip.

Scoping applies within each subroute (filter by `domain`, `kind`, `readiness`,
or subject `--seq-id`), because the catalog is large (hundreds of shell
commands). A words-only client must be able to ask "what can I do with `insert`
right now?" and get a short, ready-annotated list.

Implemented scoping in the first slice:

- `introspect facts --domain ... --seq-id ...`
- `introspect capabilities --kind ...`
- `introspect readiness --seq-id ... --readiness ...`
- `introspect verify-effects CAPABILITY_ID --seq-id ... --evidence ...`

`introspect capabilities` descriptor shape:

```json
{
  "schema": "gentle.introspection.v1",
  "capabilities": [
    {
      "id": "sequence create",
      "kind": "operation",
      "mutating": "true",
      "requires_confirmation": false,
      "args": [
        { "name": "--sequence-text", "required": true, "detail": "inline DNA sequence text" },
        { "name": "OUTPUT_ID", "required": true, "detail": "explicit id supplied with --output-id" }
      ],
      "reads": [],
      "effects": [ { "fact": "sequence.exists", "subject": { "arg": "OUTPUT_ID" },
                     "effect_kind": "must_on_success" } ],
      "precondition_expr": { "all": [] },
      "description": "Create a persistent project sequence from inline sequence text."
    },
    {
      "id": "features restriction-scan",
      "kind": "operation",
      "mutating": "false",
      "requires_confirmation": false,
      "args": [
        { "name": "SEQ_ID", "required": true, "detail": "loaded sequence id" },
        { "name": "--enzyme", "required": false, "detail": "narrow to one enzyme" }
      ],
      "reads":   [ { "fact": "sequence.exists", "subject": { "arg": "SEQ_ID" } } ],
      "effects": [ { "fact": "report.exists", "report_kind": "restriction_scan",
                     "equals": "restriction_scan", "effect_kind": "must_on_success" } ],
      "precondition_expr": { "all": [ { "fact": "sequence.exists", "subject": { "arg": "SEQ_ID" } } ] },
      "description": "Scan a loaded sequence for restriction-enzyme sites."
    },
    {
      "id": "variant materialize-allele",
      "kind": "operation",
      "mutating": "true",
      "requires_confirmation": false,
      "args": [
        { "name": "SEQ_ID", "required": true, "detail": "loaded variant-bearing sequence id" },
        { "name": "--allele", "required": true, "detail": "reference or alternate" },
        { "name": "OUTPUT_ID", "required": false, "detail": "explicit id supplied with --output-id; required for deterministic effect verification" }
      ],
      "reads": [ { "fact": "sequence.exists", "subject": { "arg": "SEQ_ID" } } ],
      "effects": [ { "fact": "sequence.exists", "subject": { "arg": "OUTPUT_ID" },
                     "effect_kind": "must_on_success" } ],
      "precondition_expr": { "all": [ { "fact": "sequence.exists", "subject": { "arg": "SEQ_ID" } } ] },
      "description": "Create one reference or alternate single-allele sequence from a variant-bearing sequence."
    },
    {
      "id": "align compute",
      "kind": "operation",
      "mutating": "false",
      "requires_confirmation": false,
      "args": [
        { "name": "QUERY_SEQ_ID", "required": true, "detail": "loaded query sequence id" },
        { "name": "TARGET_SEQ_ID", "required": true, "detail": "loaded target sequence id" }
      ],
      "reads": [
        { "fact": "sequence.exists", "subject": { "arg": "QUERY_SEQ_ID" } },
        { "fact": "sequence.exists", "subject": { "arg": "TARGET_SEQ_ID" } }
      ],
      "effects": [],
      "precondition_expr": { "all": [
        { "fact": "sequence.exists", "subject": { "arg": "QUERY_SEQ_ID" } },
        { "fact": "sequence.exists", "subject": { "arg": "TARGET_SEQ_ID" } }
      ] },
      "description": "Compute a pairwise alignment between two loaded project sequences."
    },
    {
      "id": "render-svg",
      "kind": "operation",
      "mutating": "false",
      "requires_confirmation": false,
      "args": [
        { "name": "SEQ_ID", "required": true, "detail": "loaded sequence id to render" },
        { "name": "MODE", "required": true, "detail": "linear or circular render mode" },
        { "name": "OUTPUT_PATH", "required": true, "detail": "external SVG output path" }
      ],
      "reads": [
        { "fact": "sequence.exists", "subject": { "arg": "SEQ_ID" } }
      ],
      "effects": [ { "fact": "artifact.written", "subject": { "arg": "OUTPUT_PATH" },
                     "effect_kind": "external_handoff" } ],
      "precondition_expr": { "all": [
        { "fact": "sequence.exists", "subject": { "arg": "SEQ_ID" } }
      ] },
      "description": "Render a loaded sequence map to an external SVG file."
    },
    {
      "id": "render-rna-svg",
      "kind": "operation",
      "mutating": "false",
      "requires_confirmation": false,
      "args": [
        { "name": "SEQ_ID", "required": true, "detail": "loaded RNA sequence id to render" },
        { "name": "OUTPUT_PATH", "required": true, "detail": "external SVG output path" }
      ],
      "reads": [
        { "fact": "sequence.exists", "subject": { "arg": "SEQ_ID" } }
      ],
      "effects": [ { "fact": "artifact.written", "subject": { "arg": "OUTPUT_PATH" },
                     "effect_kind": "external_handoff" } ],
      "precondition_expr": { "all": [
        { "fact": "sequence.exists", "subject": { "arg": "SEQ_ID" } }
      ] },
      "description": "Render a loaded RNA secondary-structure SVG file."
    },
    {
      "id": "render-lineage-svg",
      "kind": "operation",
      "mutating": "false",
      "requires_confirmation": false,
      "args": [
        { "name": "OUTPUT_PATH", "required": true, "detail": "external SVG output path" }
      ],
      "reads": [],
      "effects": [ { "fact": "artifact.written", "subject": { "arg": "OUTPUT_PATH" },
                     "effect_kind": "external_handoff" } ],
      "precondition_expr": { "all": [] },
      "description": "Render the current project lineage graph to an external SVG file."
    },
    {
      "id": "rna-info",
      "kind": "operation",
      "mutating": "false",
      "requires_confirmation": false,
      "args": [
        { "name": "SEQ_ID", "required": true, "detail": "loaded sequence id to inspect with the RNA structure text reporter" }
      ],
      "reads": [
        { "fact": "sequence.exists", "subject": { "arg": "SEQ_ID" } }
      ],
      "effects": [],
      "precondition_expr": { "all": [
        { "fact": "sequence.exists", "subject": { "arg": "SEQ_ID" } }
      ] },
      "description": "Inspect one loaded sequence with the RNA structure text reporter."
    },
    {
      "id": "reverse-translate run",
      "kind": "operation",
      "mutating": "true",
      "requires_confirmation": false,
      "args": [
        { "name": "PROTEIN_SEQ_ID", "required": true, "detail": "loaded protein sequence id" },
        { "name": "OUTPUT_ID", "required": false, "detail": "explicit coding-DNA sequence id supplied with --output-id; required for deterministic effect verification" }
      ],
      "reads": [
        { "fact": "sequence.exists", "subject": { "arg": "PROTEIN_SEQ_ID" } },
        { "fact": "sequence.kind", "subject": { "arg": "PROTEIN_SEQ_ID" }, "equals": "protein" }
      ],
      "effects": [ { "fact": "sequence.exists", "subject": { "arg": "OUTPUT_ID" },
                     "effect_kind": "must_on_success" } ],
      "precondition_expr": { "all": [
        { "fact": "sequence.exists", "subject": { "arg": "PROTEIN_SEQ_ID" } },
        { "fact": "sequence.kind", "subject": { "arg": "PROTEIN_SEQ_ID" }, "equals": "protein" }
      ] },
      "description": "Reverse-translate one protein sequence into a synthetic coding-DNA product."
    },
    {
      "id": "reverse-translate list-reports",
      "kind": "operation",
      "mutating": "false",
      "requires_confirmation": false,
      "args": [
        { "name": "PROTEIN_SEQ_ID", "required": false, "detail": "optional protein sequence id filter" }
      ],
      "reads": [],
      "effects": [],
      "precondition_expr": { "all": [] },
      "description": "List persisted reverse-translation reports, optionally filtered by protein sequence id."
    },
    {
      "id": "reverse-translate show-report",
      "kind": "operation",
      "mutating": "false",
      "requires_confirmation": false,
      "args": [
        { "name": "REPORT_ID", "required": true, "detail": "persisted reverse-translation report id" }
      ],
      "reads": [
        { "fact": "report.exists", "subject": { "arg": "REPORT_ID" }, "equals": "reverse_translation" }
      ],
      "effects": [],
      "precondition_expr": { "all": [
        { "fact": "report.exists", "subject": { "arg": "REPORT_ID" }, "equals": "reverse_translation" }
      ] },
      "description": "Inspect one persisted reverse-translation provenance report."
    },
    {
      "id": "primers design",
      "kind": "operation",
      "mutating": "true",
      "requires_confirmation": false,
      "args": [
        { "name": "REQUEST_JSON", "required": true, "detail": "DesignPrimerPairs operation payload or @file" },
        { "name": "TEMPLATE_SEQ_ID", "required": true, "detail": "template sequence id carried by the payload" },
        { "name": "REPORT_ID", "required": false, "detail": "explicit primer-design report id carried by the payload; required for deterministic effect verification" }
      ],
      "reads": [
        { "fact": "sequence.exists", "subject": { "arg": "TEMPLATE_SEQ_ID" } }
      ],
      "effects": [ { "fact": "report.exists", "subject": { "arg": "REPORT_ID" },
                     "equals": "primer_design", "effect_kind": "must_on_success" } ],
      "precondition_expr": { "all": [
        { "fact": "sequence.exists", "subject": { "arg": "TEMPLATE_SEQ_ID" } }
      ] },
      "description": "Generate and persist a ranked primer-pair design report."
    },
    {
      "id": "primers list-reports",
      "kind": "operation",
      "mutating": "false",
      "requires_confirmation": false,
      "args": [],
      "reads": [],
      "effects": [],
      "precondition_expr": { "all": [] },
      "description": "List persisted primer-pair design reports."
    },
    {
      "id": "primers show-report",
      "kind": "operation",
      "mutating": "false",
      "requires_confirmation": false,
      "args": [
        { "name": "REPORT_ID", "required": true, "detail": "persisted primer-design report id" }
      ],
      "reads": [
        { "fact": "report.exists", "subject": { "arg": "REPORT_ID" }, "equals": "primer_design" }
      ],
      "effects": [],
      "precondition_expr": { "all": [
        { "fact": "report.exists", "subject": { "arg": "REPORT_ID" }, "equals": "primer_design" }
      ] },
      "description": "Inspect one persisted primer-pair design report."
    },
    {
      "id": "primers design-qpcr",
      "kind": "operation",
      "mutating": "true",
      "requires_confirmation": false,
      "args": [
        { "name": "REQUEST_JSON", "required": true, "detail": "DesignQpcrAssays operation payload or seed payload" },
        { "name": "TEMPLATE_SEQ_ID", "required": true, "detail": "template sequence id carried by the payload" },
        { "name": "REPORT_ID", "required": false, "detail": "explicit qPCR design report id carried by the payload; required for deterministic effect verification" }
      ],
      "reads": [
        { "fact": "sequence.exists", "subject": { "arg": "TEMPLATE_SEQ_ID" } }
      ],
      "effects": [ { "fact": "report.exists", "subject": { "arg": "REPORT_ID" },
                     "equals": "qpcr_design", "effect_kind": "must_on_success" } ],
      "precondition_expr": { "all": [
        { "fact": "sequence.exists", "subject": { "arg": "TEMPLATE_SEQ_ID" } }
      ] },
      "description": "Generate and persist a ranked qPCR assay design report."
    },
    {
      "id": "primers list-qpcr-reports",
      "kind": "operation",
      "mutating": "false",
      "requires_confirmation": false,
      "args": [],
      "reads": [],
      "effects": [],
      "precondition_expr": { "all": [] },
      "description": "List persisted qPCR assay design reports."
    },
    {
      "id": "primers show-qpcr-report",
      "kind": "operation",
      "mutating": "false",
      "requires_confirmation": false,
      "args": [
        { "name": "REPORT_ID", "required": true, "detail": "persisted qPCR design report id" }
      ],
      "reads": [
        { "fact": "report.exists", "subject": { "arg": "REPORT_ID" }, "equals": "qpcr_design" }
      ],
      "effects": [],
      "precondition_expr": { "all": [
        { "fact": "report.exists", "subject": { "arg": "REPORT_ID" }, "equals": "qpcr_design" }
      ] },
      "description": "Inspect one persisted qPCR assay design report."
    },
    {
      "id": "primers prepare-restriction-cloning",
      "kind": "operation",
      "mutating": "true",
      "requires_confirmation": false,
      "args": [
        { "name": "REQUEST_JSON", "required": true, "detail": "PrepareRestrictionCloningPcrHandoff operation payload or @file" },
        { "name": "TEMPLATE_SEQ_ID", "required": true, "detail": "template sequence id carried by the payload" },
        { "name": "PRIMER_REPORT_ID", "required": true, "detail": "persisted primer-design report id carried by the payload" },
        { "name": "DESTINATION_VECTOR_SEQ_ID", "required": true, "detail": "destination vector sequence id carried by the payload" },
        { "name": "REPORT_ID", "required": false, "detail": "restriction-cloning handoff report id returned by execution; required for deterministic effect verification" }
      ],
      "reads": [
        { "fact": "sequence.exists", "subject": { "arg": "TEMPLATE_SEQ_ID" } },
        { "fact": "report.exists", "subject": { "arg": "PRIMER_REPORT_ID" }, "equals": "primer_design" },
        { "fact": "sequence.exists", "subject": { "arg": "DESTINATION_VECTOR_SEQ_ID" } }
      ],
      "effects": [ { "fact": "report.exists", "subject": { "arg": "REPORT_ID" },
                     "equals": "restriction_cloning_pcr_handoff", "effect_kind": "must_on_success" } ],
      "precondition_expr": { "all": [
        { "fact": "sequence.exists", "subject": { "arg": "TEMPLATE_SEQ_ID" } },
        { "fact": "report.exists", "subject": { "arg": "PRIMER_REPORT_ID" }, "equals": "primer_design" },
        { "fact": "sequence.exists", "subject": { "arg": "DESTINATION_VECTOR_SEQ_ID" } }
      ] },
      "description": "Create a restriction-site cloning PCR handoff from one persisted primer-pair report and a destination vector."
    },
    {
      "id": "primers list-restriction-cloning-handoffs",
      "kind": "operation",
      "mutating": "false",
      "requires_confirmation": false,
      "args": [],
      "reads": [],
      "effects": [],
      "precondition_expr": { "all": [] },
      "description": "List persisted restriction-site cloning PCR handoff reports."
    },
    {
      "id": "primers show-restriction-cloning-handoff",
      "kind": "operation",
      "mutating": "false",
      "requires_confirmation": false,
      "args": [
        { "name": "REPORT_ID", "required": true, "detail": "persisted restriction-cloning handoff report id" }
      ],
      "reads": [
        { "fact": "report.exists", "subject": { "arg": "REPORT_ID" }, "equals": "restriction_cloning_pcr_handoff" }
      ],
      "effects": [],
      "precondition_expr": { "all": [
        { "fact": "report.exists", "subject": { "arg": "REPORT_ID" }, "equals": "restriction_cloning_pcr_handoff" }
      ] },
      "description": "Inspect one persisted restriction-site cloning PCR handoff report."
    },
    {
      "id": "ui selection",
      "kind": "view_intent",
      "mutating": "false",
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
  are asserted by post-run verification** through
  `introspect verify-effects`. The verifier evaluates declared hard effects
  against the current fact graph plus supplied evidence.
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
| `dotplot.exists` | project | closed | no | stored dotplot analysis payload; subject id is the dotplot id |
| `flexibility_track.exists` | project | closed | no | stored sequence-flexibility track payload; subject id is the track id |
| `candidate_set.exists` | project | closed | no | persisted candidate-window set; subject id is the set name |
| `container.exists` | project | closed | no | persisted wet-lab container/pool/selection; subject id is the container id |
| `arrangement.exists` | project | closed | no | persisted lane/plate arrangement; subject id is the arrangement id |
| `rack.exists` | project | closed | no | persisted physical rack layout; subject id is the rack id |
| `workflow_macro_template.exists` | project | closed | no | persisted workflow macro template; subject id is the template name |
| `candidate_macro_template.exists` | project | closed | no | persisted candidate-window macro template; subject id is the template name |
| `macro_instance.exists` | project | closed | no | recorded workflow macro-instance lineage row; subject id is the macro-instance id |
| `construct_reasoning_graph.exists` | project | closed | no | persisted construct-reasoning graph; subject id is the graph id |
| `sequencing_trace.exists` | project | closed | no | persisted sequencing-trace evidence record; subject id is the trace id |
| `uniprot_entry.exists` | project | closed | no | stored UniProt metadata entry; subject id is the entry id |
| `ensembl_gene_entry.exists` | project | closed | no | stored Ensembl gene metadata entry; subject id is the entry id |
| `ensembl_protein_entry.exists` | project | closed | no | stored Ensembl protein metadata entry; subject id is the entry id |
| `guide_set.exists` | project | closed | no | persisted guide-design set; subject id is the guide-set id |
| `guide_filter_report.exists` | project | closed | no | persisted practical guide-filter report; subject id is the guide-set id |
| `guide_oligo_set.exists` | project | closed | no | persisted guide-oligo output set; subject id is the oligo-set id |
| `config.param` | config | closed | no | engine-owned `SetParameter` values; subject id is the parameter name, `value` is the JSON value |
| `host.tool_available` | host | closed | no | e.g. primer3/blastn resolvable; see Reconciliation below |
| `artifact.written` | host | open | yes | external file/handoff artifact written outside saved project state |

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
for `ui selection` resolves to `unsatisfied` with `ui.host_available` as the named
unmet atom — surfaced to the user as "unavailable: no view host," **not** as a
new fourth `FactTruth` value. The three-valued evaluator is preserved; the
"unavailable" distinction is a presentation of unsatisfied-with-reason.

### Reconciliation: where does `host.tool_available` live?

Tool availability already has adjacent contracts (`agents preflight`,
helper/genome `status`, primer3 preflight). `host.tool_available` should
**project from those existing probes**, not introduce a parallel detector. The
first implementation keeps it directly in the introspection `host` domain and
currently projects deterministic agent-system availability from the configured
agent catalog. This keeps `agents preflight` readiness inside the same evaluator
as project and view facts.

### Config parameters

The implemented config fact keeps the fact registry finite:
`config.param` uses `subject.kind = "other"` and `subject.id = PARAM_NAME`
instead of minting a new fact name for every parameter. The first projection is
limited to `EngineParameters` (`max_fragments_per_container`,
`require_verified_genome_anchor_for_extension`,
`genome_anchor_prepared_fallback_policy`, `primer_design_backend`, and
`primer3_executable`). Display-specific `set-param` aliases can be folded into
the same pattern later when their canonical names are normalized.

`set-param` is fact-annotated as a `host_config` capability with no project
preconditions and a `must_on_success` effect:
`config.param(PARAM_NAME) == PARAM_VALUE`, where `PARAM_VALUE` is parsed as
JSON before verification.

### Degrees of freedom

Whether `view.*` facts are projected lazily (only with a host attached) vs
always with a marker, and whether `domain` is a field on each fact vs a wrapper
per group.

## 6. Reachability and parity

Success criterion for the words-only model is twofold:

1. **No control without a word path.** Every state-mutating control and every
   visualization/mouse-driven control resolves to some agent-reachable contract
   (operation, view_intent, host_config). A GUI-only mutation path with no shared
   route is a reachability-parity gap (a bug). Continuous gestures are not
   commands, but their discrete forms are (`ui selection`, `ui zoom-to`,
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

## 10. First slice (review tweak adopted)

1. Add the `domain` discriminator and the `view` domain, including
   `ui.host_available` projected as `false` in headless contexts (implemented).
2. Add `introspect facts` and `introspect capabilities` as separate subroutes
   (implemented; `introspect readiness` and `introspect all` are also present).
3. Attach engine-owned bound self-description (`args`, `reads`, `effects` with
   `effect_kind`, `precondition_expr`) to the first fact-annotated capabilities:
   `help`, `capabilities`, `state-summary`, `state_summary`,
   `history status`, `history undo`, `history redo`,
   `load-project`, `load_project`, `save-project`, `save_project`,
   `facts graph`,
   `facts eval`, `sequence create`,
   `LoadFile`, `load_dna`, `genbank fetch`, `FetchGenBankAccession`,
   `ensembl-region fetch`, `FetchEnsemblRegion`,
   `dbsnp fetch`, `FetchDbSnpRegion`, `FetchUniprotLinkedGenBank`,
   `ImportUniprotEntrySequence`, `ensembl-gene import-sequence`,
   `ImportEnsemblGeneSequence`, `ensembl-protein import-sequence`,
   `ImportEnsemblProteinSequence`,
   `SaveFile`, `Digest`, `Pcr`, `PcrAdvanced`, `PcrMutagenesis`,
   `PcrOverlapExtensionMutagenesis`,
   `uniprot fetch`, `FetchUniprotSwissProt`,
   `uniprot import-swissprot`, `ImportUniprotSwissProt`,
   `uniprot list`, `uniprot show`,
   `ensembl-gene fetch`, `FetchEnsemblGene`,
   `ensembl-gene list`, `ensembl-gene show`,
   `ensembl-protein fetch`, `FetchEnsemblProtein`,
   `ensembl-protein list`, `ensembl-protein show`,
   `introspect facts`, `introspect capabilities`, `introspect readiness`,
   `introspect verify-effects`, `introspect all`,
   `Reverse`, `Complement`, `ReverseComplement`, `Branch`, `ExtractRegion`,
   `features restriction-scan`, `features query`, `features export-bed`,
   `ExportFeaturesBed`,
   `InspectSequenceContextView`, `ExportSequenceContextBundle`,
   `inspect-feature-expert`, `render-feature-expert-svg`,
   `features tfbs-summary`, `SummarizeTfbsRegion`,
   `AnnotateTfbs`, `QueryProteinResidueGenomicCoordinates`,
   `variant annotate-promoters`, `AnnotatePromoterWindows`,
   `variant promoter-context`, `SummarizeVariantPromoterContext`,
   `variant reporter-fragments`, `SuggestPromoterReporterFragments`,
   `variant materialize-allele`, `MaterializeVariantAllele`,
   `align compute`, `AlignSequences`,
   `dotplot compute`, `ComputeDotplot`,
   `dotplot overlay-compute`, `ComputeDotplotOverlay`,
   `dotplot show`, `RenderDotplotSvg`,
   `flex compute`, `ComputeFlexibilityTrack`, `flex show`,
   `candidates generate`, `GenerateCandidateSet`,
   `candidates generate-between-anchors`, `GenerateCandidateSetBetweenAnchors`,
   `candidates show`, `candidates metrics`,
   `candidates score`, `ScoreCandidateSetExpression`,
   `candidates score-distance`, `ScoreCandidateSetDistance`,
   `candidates score-weighted`, `ScoreCandidateSetWeightedObjective`,
   `candidates filter`, `FilterCandidateSet`,
   `candidates top-k`, `TopKCandidateSet`,
   `candidates pareto`, `ParetoFrontierCandidateSet`,
   `candidates set-op`, `CandidateSetOp`,
   `candidates delete`, `DeleteCandidateSet`,
   `candidates template-show`, `candidates template-put`,
   `UpsertCandidateMacroTemplate`, `candidates template-delete`,
   `DeleteCandidateMacroTemplate`, `candidates template-run`,
   `guides list`, `guides put`, `UpsertGuideSet`,
   `guides show`, `guides delete`, `DeleteGuideSet`,
   `guides filter`, `FilterGuidesPractical`, `guides filter-show`,
   `guides oligos-generate`, `GenerateGuideOligos`,
   `guides oligos-list`, `guides oligos-show`,
   `guides oligos-export`, `ExportGuideOligos`,
   `guides protocol-export`, `ExportGuideProtocolText`,
   `render-svg`, `render-rna-svg`, `render-lineage-svg`, `RenderSequenceSvg`,
   `RenderFeatureExpertSvg`, `RenderTfbsScoreTrackCorrelationSvg`,
   `RenderRnaStructureSvg`, `RenderLineageSvg`,
   `rna-info`, `reverse-translate run`, `ReverseTranslateProteinSequence`,
   `DeriveTranscriptSequences`, `DeriveProteinSequences`,
   `DeriveSplicingReferences`,
   `reverse-translate list-reports`,
   `reverse-translate show-report`, `reverse-translate export-report`,
   `proteases list`, `proteases show`, `proteases digest`,
   `ProteaseDigestProteinSequence`, `proteases digest-gel-svg`,
   `RenderProteaseDigestGelSvg`, `RenderProteinGelSvg`,
   `RenderProteinGelReportsSvg`, `RenderProtein2dGelSvg`,
   `primers design`, `DesignPrimerPairs`, `DesignInsertionPrimerPairs`,
   `primers list-reports`,
   `primers show-report`, `primers export-report`, `primers design-qpcr`,
   `DesignQpcrAssays`, `primers list-qpcr-reports`, `primers show-qpcr-report`,
   `primers export-qpcr-report`, `ExportPrimerDesignReport`,
   `primers prepare-restriction-cloning`,
   `primers list-restriction-cloning-handoffs`,
   `primers show-restriction-cloning-handoff`,
   `primers export-restriction-cloning-handoff`,
   `primers preflight`, `primers seed-from-feature`,
   `primers seed-from-splicing`,
   `primers restriction-cloning-vector-suggestions`,
   `primers seed-restriction-cloning-handoff`,
   `seq-trace import`, `ImportSequencingTrace`,
   `seq-trace list`, `ListSequencingTraces`,
   `seq-trace show`, `ShowSequencingTrace`,
   `seq-confirm list-reports`, `seq-confirm run`, `seq-confirm show-report`,
   `seq-confirm export-report`, `seq-confirm export-support-tsv`,
   `ListSequencingConfirmationReports`,
   `ConfirmConstructReads`,
   `ShowSequencingConfirmationReport`,
   `ExportSequencingConfirmationReport`,
   `ExportSequencingConfirmationSupportTsv`,
   `seq-primer suggest`, `SuggestSequencingPrimers`,
   `cutrun list-read-reports`, `cutrun show-read-report`,
   `cutrun export-coverage`, `ListCutRunReadReports`, `ShowCutRunReadReport`,
   `ExportCutRunReadCoverage`, `rna-reads list-reports`,
   `rna-reads show-report`, `rna-reads export-report`,
   `rna-reads export-hits-fasta`, `rna-reads export-target-quality`,
   `rna-reads export-paths-tsv`, `rna-reads export-abundance-tsv`,
   `rna-reads export-score-density-svg`, `rna-reads export-alignments-tsv`,
   `rna-reads export-isoform-triage-tsv`,
   `rna-reads export-alignment-dotplot-svg`,
   `ListRnaReadReports`, `ShowRnaReadReport`, `ExportRnaReadReport`,
   `ExportRnaReadHitsFasta`, `ExportRnaReadTargetQuality`,
   `ExportRnaReadExonPathsTsv`, `ExportRnaReadExonAbundanceTsv`,
   `ExportRnaReadScoreDensitySvg`, `ExportRnaReadAlignmentsTsv`,
   `ExportRnaReadIsoformTriageTsv`, `ExportRnaReadAlignmentDotplotSvg`,
   `rna-reads align-report`, `AlignRnaReadReport`,
   `rna-reads preflight-isoforms`, `PreflightRnaReadIsoforms`,
   `rna-reads materialize-hits`, `MaterializeRnaReadHitSequences`,
   `SummarizeRnaReadGeneSupport`, `InspectRnaReadGeneSupport`,
   `SummarizeJasparEntries`, `BenchmarkJasparRegistry`,
   `ListJasparCatalog`, `ResolveTfQueries`, `ListReporterCatalog`,
   `RecommendReporters`,
   `reporters list`, `reporters recommend`, `reporters export-corpus`,
   `ExportReporterCorpus`,
   `services status`, `services providers list`,
   `services providers doctor`,
   `planning consult cloning`, `planning protein-expression-handoff`,
   `planning profile show`, `planning profile set`, `planning profile clear`,
   `planning objective show`, `planning objective set`,
   `planning objective clear`, `planning suggestions list`,
   `planning suggestions accept`, `planning suggestions reject`,
   `planning sync status`, `planning sync pull`, `planning sync push`,
   `resources summarize-jaspar`, `resources status`, `cache clear`,
   `resources sync-rebase`, `sync_rebase`,
   `resources sync-jaspar`, `sync_jaspar`,
   `resources sync-ucsc-rmsk`,
   `resources import-gene-list-cache`,
   `resources import-ontology-assignment-cache`,
   `resources import-co-regulated-cache`,
   `resources install-ucsc-rmsk`,
   `resources prepare-ucsc-rmsk-index`,
   `resources sync-jaspar-remote-metadata`,
   `SyncJasparRemoteMetadata`, `resources benchmark-jaspar`,
   `resources prepare-publication-dataset`,
   `resources suggest-ucsc-rmsk-index`, `resources list-jaspar`,
   `resources inspect-jaspar`, `resources resolve-tf-query`,
   `resources list-publication-datasets`,
   `resources status-publication-dataset`,
   `genomes validate-catalog`, `helpers validate-catalog`,
   `helpers vocabulary list`, `helpers vocabulary doctor`,
   `services delivery-route`, `services project-preflight`,
   `services project-quote`, `services handoff`, `services guide`,
   `hosts list`, `list_host_profile_catalog_entries`,
   `host_profile_catalog_entries`, `list_helper_catalog_entries`,
   `helper_catalog_entries`, `helper_semantics_vocabulary`,
   `helper_interpretation`, `proteases list`, `proteases show`,
   `mirna explain-seed`, `mirna catalog-show`,
   `ladders list`, `inspect_dna_ladders`, `inspect_rna_ladders`,
   `list_dna_ladders`, `list_rna_ladders`, `ladders export`,
   `export_dna_ladders`, `export_rna_ladders`,
   `ExportDnaLadders`, `ExportRnaLadders`,
   `agents list`, `agent_systems`, `list_agent_systems`,
   `agents ask`, `ask_agent_system`, `agents plan`, `agent_plan`,
   `agents execute-plan`, `agent_execute_plan`,
   `reads acquire status`, `reads acquire prepare`,
   `reads acquire inspect`, `reads acquire cancel`,
   `ReadAcquireStatus`, `ReadAcquirePrepare`,
   `ReadAcquireInspect`, `ReadAcquireCancel`,
   `protocol-cartoon list`, `protocol-cartoon render-svg`,
   `protocol-cartoon render-template-svg`,
   `protocol-cartoon render-with-bindings`,
   `protocol-cartoon template-export`, `protocol-cartoon template-validate`,
   `RenderProtocolCartoonSvg`, `RenderProtocolCartoonTemplateSvg`,
   `RenderProtocolCartoonTemplateWithBindingsSvg`,
   `ExportProtocolCartoonTemplateJson`, `ValidateProtocolCartoonTemplate`,
   `cache inspect`, `cutrun list`, `ListCutRunDatasets`, `cutrun status`,
   `ShowCutRunDatasetStatus`, `arrays inspect-microarray-track`,
   `arrays inspect-probe-region-output`,
   `arrays render-probe-region-output-svg`,
   `tracks import-bed`, `import_genome_bed_track`, `ImportGenomeBedTrack`,
   `tracks import-bigwig`, `import_genome_bigwig_track`,
   `ImportGenomeBigWigTrack`,
   `tracks import-vcf`, `import_genome_vcf_track`, `ImportGenomeVcfTrack`,
   `ImportBlastHitsTrack`, `genomes blast-track`, `helpers blast-track`,
   `tracks tracked list`, `tracks tracked add`, `tracks tracked remove`,
   `tracks tracked clear`, `tracks tracked apply`,
   `arrays project-microarray-track`, `ProjectMicroarrayTrack`,
   `arrays project-probe-region-output`,
   `FindRestrictionSites`, `features tfbs-score-tracks-svg`,
   `RenderTfbsScoreTracksSvg`, `SummarizeTfbsScoreTracks`,
   `features tfbs-track-similarity`, `SummarizeTfbsTrackSimilarity`,
   `candidates list`, `candidates template-list`, `guides list`,
   `macros instance-list`, `macros template-list`, `routines list`,
   `routines explain`, `routines compare`,
   `construct-reasoning list-graphs`, `construct_reasoning_graphs`,
   `dotplot list`, `flex list`,
   `genomes list`, `helpers list`, `genomes status`, `helpers status`,
   `genomes genes`, `helpers genes`,
   `genomes ensembl-available`, `helpers ensembl-available`,
   `ensembl_installable_genomes`, `list_ensembl_installable_genomes`,
   `genomes preview-ensembl-specs`, `helpers preview-ensembl-specs`,
   `genomes update-ensembl-specs`, `helpers update-ensembl-specs`,
   `list_reference_genomes`, `list_reference_catalog_entries`,
   `reference_catalog_entries`,
   `is_reference_genome_prepared`, `list_reference_genome_genes`,
   `genomes blast-status`, `helpers blast-status`,
   `genomes blast-list`, `helpers blast-list`,
   `blast_async_status`, `blast_async_list`,
   `macros instance-show`,
   `macros template-show`, `macros template-put`,
   `UpsertWorkflowMacroTemplate`, `macros template-delete`,
   `DeleteWorkflowMacroTemplate`, `macros template-import`,
   `macros template-run`,
   `construct-reasoning show-graph`, `construct_reasoning_graph`,
   `construct-reasoning list-inspection-actions`,
   `construct_reasoning_inspection_actions`,
   `construct-reasoning run-inspection-action`,
   `construct_reasoning_run_inspection_action`,
   `construct-reasoning set-annotation-status`,
   `construct_reasoning_set_annotation_status`,
   `construct-reasoning write-annotation`,
   `construct_reasoning_write_annotation`,
   `construct-reasoning export-graph`,
   `ensembl-gene list`,
   `ensembl-protein list`, `gene-groups list`, `gene-groups show`,
   `gene-groups resolve`, `gene-groups doctor`,
   `display`,
   `SetLinearViewport`, `ui intents`, `ui_intents`,
   `ui open`, `ui focus`, `ui close`, `ui_intent`,
   `ui_prepared_genomes`, `ui_latest_prepared`, `ui selection`, `set-param`,
   `SetTopology`, `RecomputeFeatures`, `agents preflight`,
   `agents discover-models`, `agent_preflight`, and `agent_models`
   (implemented).
4. Add readiness tests proving the **same descriptor** resolves correctly in both
   a GUI-attached and a headless context (implemented through explicit
   `--ui-host true|false` projection).
5. Defer: full catalog fact annotation, glossary-as-projection inversion, and
   any visualization-interpretation work.

## 11. Resolved decisions

These were open questions in the prior draft; the review and implemented first
slice resolved them. The non-negotiables are promoted into `docs/decisions.md`.

- **Routes:** subroutes + aggregate, not one scoped route (section 3).
- **Domains:** additive `domain` field; landed facts unchanged (section 5).
- **Headless:** explicit `ui.host_available` fact; "unavailable" is
  unsatisfied-with-reason, not a new truth value (section 5).
- **Effects:** carry `effect_kind`; only `must_on_success` is verified by
  `introspect verify-effects` (section 4).
- **Args source (long term):** an engine/protocol-side descriptor becomes the
  root truth and `docs/glossary.json` becomes a projection of it. Direction
  accepted; the inversion is **deferred** past the first slice.
- **Domain representation:** projected facts carry a per-fact `domain` field
  and `introspect facts` additionally groups rows by domain for client
  convenience.
- **Bound readiness route:** bound readiness lives in `introspect readiness`.
  `introspect capabilities` remains catalog/descriptor discovery.

## 12. Remaining implementation work

- Add fact-aware annotations for more glossary/engine capability rows.
- Migrate glossary/help generation toward projection from an engine/protocol-side
  descriptor once enough annotated rows prove the shape.
