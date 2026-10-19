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
RUN="$(cd "$RUN" && pwd -P)"
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
| `explanation.time_policy.hours` | 6 h | Preserved legacy proxy for ranking and the time guardrail; not measured bench time |

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

Giving one procedure a richer description must not, by itself, make it win a
comparison. `composite_meta_score` and the `max_time_hours` guardrail therefore
retain the **old family/base scalar plus delays** for every routine, separately
recorded in `explanation.time_policy`. Here that proxy is 6 h, even though the
displayed synthetic critical path is 3.2 h. Changing the bench model cannot
change ordering or a guardrail verdict. A guardrail pass is not a promise that
the experiment, or even the displayed maximum scenario, fits within that time.
Historical reports missing the policy record display its basis as unknown;
do not guess it from their elapsed value. A calibrated timing comparison would
need a separately reviewed policy and suitable evidence for every method.

Try the policy limit explicitly while the in-stock profile is still active:

```bash
g planning consult cloning --objective '{"schema":"gentle.planning_objective.v1","max_time_hours":5}' --format json > "$RUN/time-guardrail.json"
```

Find the `restriction` row in `strategy_candidates`. Its estimate has
`passes_guardrails=false` and a `max_time_hours_exceeded (6.0 > 5.0)` failure,
despite displaying 3.2 h. The two numbers answer different questions. This
advisory rejection neither changes the protocol nor executes another method.

### 3.1 Read The Protocol Overview Without Running A Routine

In **Patterns -> Routine Assistant**, find and select
`restriction.digest_ligate_extract_sticky`. In the comparison stage, switch
from **Explanation** to **Protocol Overview**. This is a dependency diagram,
not an execution record or a validated Promega protocol.

- Locate the two independent roots: digest-material preparation and preparation
  of later ligation supplies. Follow both arrows into the ligation setup path.
  This explains why preparation need not be added serially to the critical
  path; it does not promise that one operator or shared equipment can do both.
- Select or hover a step to inspect its source, prerequisites and
  minimum/typical/maximum ranges. Magenta marks the **recorded typical** critical
  path, not the minimum/maximum scenario or a completed laboratory action.
- Zoom, fit or scroll the overview. Horizontal spacing is dependency depth,
  **not an elapsed-time axis**. The 6-hour policy proxy remains separate from
  the displayed 3.2-hour scenario; neither is changed by opening this tab.
- Use **Copy SVG** for the portable figure with hover details. **Copy JSON**
  retains the full explanation and the estimate needed to check its binding.
- Select the Gibson alternative: its bench-step overview must say unavailable,
  not display invented steps or an empty successful procedure. Legacy reports
  without the new record also cannot establish a verified diagram.

The graph currently has **seven synthetic bench steps**, not molecular inputs
or outputs. It deliberately does not draw your vector, insert or assembled
product: those identities and transfers are not bound by this model. GENtle's
existing molecular `protocol-cartoon` illustrations remain separate, as does
the lineage graph of actual project operations. Do not run a routine merely
to obtain a diagram.

### 3.2 Retain The Same Figure From Shared JSON

The CLI, GUI and agent-facing shared shell return the same
`planning.protocol_overview`, schema `gentle.planning_protocol_overview.v1`.
MCP/script callers use the existing agent-plan executor with a read-only shell
candidate `routines explain restriction.digest_ligate_extract_sticky`; MCP's
tool is `agent_execute_plan`, not a new protocol-specific tool. This exercises
the deterministic command path and requires no inner model or supplier call.
Using the `in-stock.json` retained above:

```bash
python3 - "$RUN/in-stock.json" "$RUN/protocol-overview.svg" <<'PY'
import json
import sys
with open(sys.argv[1], encoding="utf-8") as handle:
    explanation = json.load(handle)
overview = explanation["planning"]["protocol_overview"]
if overview["status"] != "available":
    raise SystemExit(overview["reason"])
with open(sys.argv[2], "w", encoding="utf-8", newline="\n") as handle:
    handle.write(overview["svg"])
print(overview["source_estimate_sha256"], overview["content_sha256"])
PY
```

Keep the full JSON alongside its SVG. The content bindings identify this
routine/estimate/projection; they are not authenticated search, supplier or
laboratory evidence. The GUI verifies the diagram against the actual estimate,
including changed payload bytes, rather than trusting a claimed digest alone.
Opening the view and copying these records do not change project state,
ranking, guardrails, stock, lineage or execution approvals.

**Native acceptance remains pending.** From a matching verified GUI build,
retain real screenshots of: the complete seven-step overview with its timing
proxy; a selected step's source/ranges; and the Gibson unavailable state.
Record the exact revision, explanation/content hashes and active routine.
Compare copied SVG node/edge IDs and hover details with JSON, and verify that
viewing/copying leaves the project unchanged. No screenshots are fabricated;
these checkpoints are not established by the authored headless tests.

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
The ranking/time-guardrail proxy separately becomes `6 + 84 = 90 h`, not
87.2 h. Machine charges still use the original 6 h occupancy approximation,
not either delay-inclusive figure. This keeps the three purposes inspectable.

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

### 6.1 Review A Disposable Supplier Update Without Changing Vendor Facts

Why separate five lanes? A vector's biological identity, a seller's package,
the protocol PDF, a reviewed procedure and your lab's stock observations answer
different questions. A protocol link is not proof of DNA identity, pricing or
an executable recipe. The new catalog routes expose those differences:

```bash
g catalog doctor > "$RUN/catalog-health.json"
g catalog search --entity product --query luciferas > "$RUN/product-identities.json"
g catalog search --entity provider --query promega > "$RUN/promega-record.json"
g catalog search --entity protocol --query luciferase --limit 20
```

Product rows are existing helper/vector/reference identities, not automatically
orderable products. Provider rows include full config and `record_sha256`,
distinct from the raw `source_sha256`. Protocol rows have a linked document,
`document_sha256=null` and `procedure_review=not_evaluated`. Commercial offers
and reviewed vendor procedures are unavailable; private profiles are not
inspected. None of this changes the earlier time, cost, scores or guardrails.

Use only the disposable roots established in section 2 for this optional
exercise. The following complete record is hand-crafted synthetic data, not
Promega, not a supplier offer and not a reviewed protocol:

```bash
python3 - "$RUN/provider-update.json" <<'PY'
import json, sys
record = dict(provider="synthetic_tutorial_supplier", display_name="Synthetic tutorial only",
              support_status="catalog_only", website_url="https://example.invalid",
              dashboard_url="https://example.invalid", api_documentation_url=None)
for key in ("channels", "product_templates", "validation_rules", "default_delivery_hints",
            "default_purification_hints", "default_qc_hints", "required_followup",
            "account_enablement_notes", "warnings"):
    record[key] = []
record["capabilities"] = [dict(service_kind="synthetic_protocol", track="dna",
    display_name="Synthetic", quote_handoff_supported=False,
    direct_api_documented=False, direct_api_implemented=False,
    supported_submission_modes=[], status_tracking="none", artifact_kinds=[], notes=[])]
bundle = dict(schema="gentle.provider_update_bundle.v1", changes=[dict(
    action="replace", provider=record["provider"], base_record_sha256=None, record=record)])
with open(sys.argv[1], "w", encoding="utf-8", newline="\n") as handle:
    json.dump(bundle, handle, indent=2)
    handle.write("\n")
PY
g catalog updates preview "$RUN/provider-update.json" --scope project --output "$RUN/provider-preview.json"
```

Read the preview **before** approving it. Expect `add`, null previous hash,
`protocol_document_evidence_unavailable`, and the exact project destination
under `$GENTLE_PROJECT_ROOT/.gentle/catalogs/`. Confirm `applicable=true` and
review source bindings and impact limitations. Preview writes only its fresh
artifact, not the active catalog. Do not substitute a real supplier casually.
For this synthetic, disposable example only, after review:

```bash
DIGEST="$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["preview_sha256"])' "$RUN/provider-preview.json")"
g catalog updates apply "$RUN/provider-preview.json" --approve-digest "$DIGEST" > "$RUN/provider-receipt.json"
g catalog search --entity provider --query synthetic_tutorial_supplier
g catalog doctor
```

The receipt identifies retained history and the active raw-byte hash. Check
that `proposed.json` matches that hash and that no `previous.bin` exists for
the previously absent project catalog. Existing built-in suppliers remain.
Repeating the old approval fails: current sources changed. A new preview of
the old null-base addition reports a conflict, not permission to overwrite.
For a real replacement, obtain the current **full** config and semantic base
hash from provider search, then review a new complete bundle.

Withdrawal uses `action="withdraw"` with that current base hash, **without** a
record. It installs an explicit tombstone, preventing a lower overlay from
silently reappearing. Do not delete a catalog row to imply withdrawal. Both
writing commands require agent confirmation, even if a model calls them
read-only. They do not alter project undo history or contact a supplier.
Source changes, masking by later `.d` files or tampering require a new preview.
No-op-only previews are not applicable and do not rewrite files.

See [bundle limits and interrupted-apply recovery](../catalog_updates.md):
never auto-replay prepared history, delete receipts to make a retry pass, or
interpret a content hash as authenticated vendor approval. Joined offers,
retained/reviewed documents, quantities/taxes/shipping, procedure extraction
and inbox automation remain later work, not features of this example.

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
