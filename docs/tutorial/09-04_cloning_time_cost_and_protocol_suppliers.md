# 09.04 - Plan Cloning Work, Waiting And Cost With Protocol-Publishing Suppliers

This offline walkthrough is for a wet-lab biologist deciding how much bench
attention, waiting and purchasing a cloning idea might require. Promega's
pGL4.10[luc2] is the worked product example because the owner prefers it, not
because GENtle has ranked it above other suppliers.

**Outcome:** inspect an advisory cloning estimate, explain its four cost
components, distinguish missing data from zero, and find suppliers with readable
protocols using explicit country and document-format metadata. No construct is
made, supplier contacted, sequence downloaded or order prepared.

## 1. Keep Four Kinds Of Evidence Apart

| Question | Evidence to inspect | What it does not establish |
| --- | --- | --- |
| Which product? | Promega pGL4.10[luc2], E6651; official product page and TM259 manual | That another pGL4 vector is interchangeable |
| Which actual DNA? | Existing helper identity AY738222.1, 4,242 bp, circular; later verify the acquired sequence/hash and annotations | That a product name proves local sequence identity |
| Which laboratory procedure? | Current supplier protocol, exact enzymes/insert/vector, and local review | That a readable PDF is an executable GENtle workflow |
| How much effort and cost? | Shared planning estimate and the declared local profile below | That synthetic defaults are measured Promega protocol durations or a quotation |

pGL4.10[luc2] is a promoterless luciferase reporter vector: an inserted region
can be investigated in a reporter experiment. A vector listing or planning
score does not predict luciferase activity, autonomous promoter function or
successful cloning. This is not a generic high-yield protein-expression vector;
use [09.03](./09-03_protein_expression_handoff_cli.md) for that separate question.

Official references: [Promega product E6651](https://www.promega.com/products/luciferase-assays/genetic-reporter-vectors-and-cell-lines/promoterless-firefly-luciferase-basic-vectors/?catNum=E6651),
[pGL4 technical manual TM259](https://worldwide.promega.com/-/media/files/resources/protocols/technical-manuals/0/pgl4-luciferase-reporter-vectors-protocol.pdf?rev=395cb76a0eed4807ae0bbd6f67b7bb77&sc_lang=en).
Opening these is an optional explicit network action; the exercise does not
retrieve them. The commercial vector sequence is not bundled.

## 2. Start With A Disposable Offline Project

From the repository root, use an already-built `gentle_cli` on your PATH:

```bash
ROOT="$PWD"
GENTLE_CLI="$(command -v gentle_cli)"
test -n "$GENTLE_CLI"
RUN="$(mktemp -d "${TMPDIR:-/tmp}/gentle-planning-09-04.XXXXXX")"
STATE="$RUN/tutorial.gentle.json"
export GENTLE_ASSET_ROOT="$ROOT"
export GENTLE_SYSTEM_CONFIG_ROOT="$RUN/system"
export XDG_CONFIG_HOME="$RUN/user"
export GENTLE_PROJECT_ROOT="$RUN/project"
g() { "$GENTLE_CLI" --state "$STATE" --full-report "$@"; }
g planning profile show --scope effective
g planning profile set @docs/examples/planning/cloning_in_stock.json --scope project_override
```

The profile-setting step is an explicit project mutation. Use this disposable
state, not your study. Read the [synthetic input provenance](../examples/planning/README.md).
Its EUR amounts and availability are teaching values, not Promega prices,
physical stock checks or an assertion about what your laboratory owns.

## 3. Separate Bench Attention From Waiting

```bash
g routines explain restriction.digest_ligate_extract_sticky > "$RUN/in-stock.json"
g planning consult cloning --format text
g planning consult cloning --objective '{"schema":"gentle.planning_objective.v1","helper_profile_id":"Reporter Promega Luciferase AY738222 (online)"}' --format json > "$RUN/promega-context.json"
g routines compare restriction.digest_ligate_extract_sticky gibson.two_fragment_overlap_preview
```

In `in-stock.json`, find `planning.estimate`, schema
`gentle.planning_estimate.v2`. The same estimate is used by consultation,
routine listing, explanation, comparison and Routine Assistant.
The second consultation names the existing Promega helper explicitly; inspect
`objective_summary.objective.helper_profile_id` and `vector_candidates` to see
that choice. It changes planning context, not supplier selection or sequence
validation. It does not prepare/download the helper or prove enzyme/insert
compatibility. The unbound consultation is retained for comparison.

| Typical scenario, no delays | Expected value | Why you care |
| --- | --- | --- |
| `hands_on_hours.typical` | 0.8 h | Bench attention is not the whole time until the next decision |
| `unattended_hours.typical` | 2.6 h | Incubation still occupies elapsed time without continuous human work |
| `elapsed_hours.typical` | 3.2 h | Dependency-compatible preparation overlaps; do not add every step serially |
| `step_model` / `estimated_time_hours_basis` | `present` / `step_model_critical_path` | These identify the origin of the displayed scalar |

Inspect `schedule[]`: step IDs, dependencies, earliest start/finish ranges and
source notes explain the calculation. `explanation.critical_paths` records a
separate path for minimum, typical and maximum scenarios. These are scenarios,
not probabilities or confidence intervals. Human work plus unattended work is
3.4 h, but the dependency critical path is 3.2 h.

The branch preparing later ligation materials can overlap the digest branch.
This is idealized dependency-only overlap, **not** a staff, equipment-capacity,
calendar or collision-free booking schedule. Transformation, colony screening,
sequencing and experimental success are not silently included. The figures are
hand-authored synthetic allowances, not a Promega protocol or timings measured
in a laboratory. The in-silico digest/ligation/extraction macro is a different
thing and is not executed by `routines explain`.

In the Gibson comparison, the three breakdown ranges are `null`/`unknown`,
`step_model=absent`, and the scalar has basis `legacy_family_scalar`.
Only the restriction/ligation routine currently has a step model. A missing
breakdown does not mean zero effort, nor does comparing differently based
scalars prove which method is faster in your laboratory.

## 4. Explain Cost Before Comparing It

The synthetic in-stock estimate has a **known subtotal** of 39 EUR:

| Field | In stock | Meaning |
| --- | --- | --- |
| `routine_base_cost` | 0 | Existing routine base; not proof of no overhead |
| `stock_drawn_material_cost` | 9 | Per-class allowances 2 + 7; not quantity-aware consumption |
| `missing_material_purchase_proxy_cost` | 0 | No required profiled material is marked unavailable |
| `machine_run_cost` | 30 | Synthetic rate 5 x legacy 6 h occupancy basis |
| `estimated_cost` | 39 | Sum of the four components, not an orderable price |

The machine approximation deliberately still uses the legacy base duration,
not the new 3.2 h critical path or 0.8 h hands-on time. Inspect
`explanation.machine_run_cost_basis_hours`, `material_cost_rows`,
`machine_cost_rows`, `unknown_material_classes`, `unpriced_material_classes`,
`unpriced_machine_classes` and `cost_assumptions` before treating a subtotal
as useful. Unpriced/unknown is not free. Package size, tax, shipping, minimum
orders, expiry, remaining aliquots and repeat allowances are not modelled.

## 5. Introduce Procurement And An Equipment Queue

```bash
g planning profile set @docs/examples/planning/cloning_waiting.json --scope project_override
g routines explain restriction.digest_ligate_extract_sticky > "$RUN/waiting.json"
```

Only two facts change: ligase is unavailable with a two-business-day delay,
and the thermocycler queue is half a business day. Hands-on and unattended
ranges stay unchanged. The disclosed delay prefix adds
`2.5 x 24 x 7/5 = 84 h`, making typical elapsed time 87.2 h.
This is the existing business-day approximation, not a promised completion
date: no start date, holidays, weekend calendar or vendor delivery is known.

The known subtotal stays 39 EUR, but stock-draw value becomes 2 and the
missing-material purchasing proxy becomes 7. The distinction helps explain
why a method can consume stock yet require no new purchase, or wait for a
purchase without adding human bench work. It does not trigger procurement.

## 6. Find Suppliers By Protocol And Country, Not Brand Guessing

```bash
g services providers doctor
g services providers list > "$RUN/suppliers.json"
python3 -c 'import json,sys; c=json.load(open(sys.argv[1])); print(json.dumps([p for p in c["providers"] if p.get("country_code")=="DE" and p.get("protocols")],indent=2))' "$RUN/suppliers.json"
python3 -c 'import json,sys; c=json.load(open(sys.argv[1])); print(json.dumps([p for p in c["providers"] if p["provider"]=="promega"],indent=2))' "$RUN/suppliers.json"
```

Replace `DE` with another ISO alpha-2 country code. The new directory is global,
not Europe-only; it is source-backed but **non-exhaustive**. An unknown country
or absent protocol list is unassessed for this filter, not evidence that a
supplier lacks protocols. `country_code` identifies the named company/operation,
not shipping coverage, ownership or manufacturing. Promega's named company is
`US`, even though it has European offices.

Each `protocols[]` item has `title`, `url`, `format`, `checked_on` and scope
notes. Both `html` and `text_pdf` are readable documents, not structured JSON/XML
instructions or an executable GENtle workflow. Filter a provider's protocols
by `format` if needed. A product page or order form alone is not a protocol.
The [directory and official sources](../../assets/external_service_providers.d/README.md)
explain inclusion, country sources, older manuals and how to add further
qualifying suppliers. Reviewing current license, product and protocol revision
remains a human responsibility; no live stock or commercial status is queried.

Try the deliberate no-adapter boundary:

```bash
g services project-preflight '{"schema":"gentle.external_service_request.v1","provider":"promega","service_kind":"vector_supply","source_target":{"kind":"catalog_vector","catalog_number":"E6651"}}'
```

Expected: `capability_status=catalog_only`, `eligible=false`, both availability
flags false, and turnaround/cost hint `null`. A published protocol does not
create a GENtle adapter. The catalog number also supplies no actual DNA, so
source validation can report a missing sequence as well; do not invent vector
DNA to clear that message. Listing Promega must not inherit a GeneArt turnaround
claim. GeneArt and metabion keep their separate, existing local handoff
contracts; none of these commands submits an order or contacts a provider.

## 7. Inspect The Same Facts In The GUI Or Through An Agent

1. Open a disposable project, then `Patterns -> Routine Assistant...`.
2. Select the restriction/ligation routine and inspect its shared planning
   estimate, labels and comparison with the unmodelled Gibson alternative.
3. In the GUI Shell, inspect/set the synthetic project profile only after
   explicit approval; do not assume the GUI and a separate CLI process share
   unsaved state. Compare the same effective profile before comparing numbers.
4. Open `Services -> External Services...`, refresh the catalog and inspect
   Promega's discovery-only status and raw shared catalog metadata. The
   country/format filtering above operates on JSON; this chapter does not
   claim a new graphical country-filter control.
5. Use local preflight to see the blocked discovery-only result. Do not use
   submission, external browsing or account credentials for this exercise.

The inner agent and an outer CLI/MCP agent should request these shared shell
commands and explain the returned fields, not invent another estimate or
promote a supplier listing to validated ordering. Profile edits retain mutation
approval. These are existing shared-shell routes; the chapter does not claim
a new named MCP tool or a live inner-model/native-GUI acceptance run.

## Acceptance And Limits

- The optional real-CLI companion in `scripts/test_specialist_tutorials.py`
  exercises both synthetic profiles, cost reconciliation, unknown legacy
  breakdowns, filterable supplier metadata and blocked Promega preflight.
- Source checks cover catalog registration, local links and fixture intent;
  the existing LF/CRLF tutorial-checkout gate covers the generated catalog.
- No screenshots are fabricated. Native GUI review/capture and human
  scientific/procurement acceptance remain pending; this chapter is a
  `manual/hybrid` guide, not another completed executable workflow replay.
- No preferred supplier, fastest method, real completion date, quote,
  experimental success or reporter-function claim follows from this exercise.
