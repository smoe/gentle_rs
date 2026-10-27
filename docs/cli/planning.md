# Experiment Planning

<!-- manual-counterparts:start -->
Related GUI reference: [Protein Expression Handoff](../gui/protein-expression-handoff.md) | [Experiment Planning](../gui/experiment-planning.md).
These topic links do not imply identical feature coverage across interfaces.
<!-- manual-counterparts:end -->

Planning meta-layer commands (`gentle_cli planning ...` or `gentle_cli shell 'planning ...'`):

- Planning schemas:
  - `gentle.planning_profile.v1`
  - `gentle.planning_objective.v1`
  - `gentle.planning_estimate.v2` (output-only; profile/objective stay v1)
  - `gentle.planning_suggestion.v1`
  - `gentle.planning_sync_status.v1`
  - `gentle.planning_cloning_consultation.v1`
  - `gentle.protein_expression_handoff.v2` (legacy v1 reports remain readable)
- Effective profile merge precedence:
  - `global_profile -> confirmed_agent_overlay -> project_override`
- Schema compatibility:
  - payloads declaring a non-matching planning schema id are rejected
    (not auto-upgraded silently).
- `planning consult cloning [--seq-id SEQ_ID] [--objective JSON_OR_@FILE] [--profile-scope effective] [--format json|text]`
  - Returns a read-only cloning strategy/vector consultation from the effective
    profile, current or supplied objective, host/helper catalogs, and existing
    routine estimate logic.
  - `PlanningObjective.biological_intent` may be set to
    `protein_expression_max_yield`; phrase-like values such as "give me the
    maximal amount of protein" are normalized to the same intent.
  - For `protein_expression_max_yield`, the consultation adds explicit
    missing questions for total/soluble/active/purified/secreted yield,
    expression chassis, folding/cofactor requirements, toxicity/induction, and
    scale/purification endpoint before any construct route is accepted.
  - `strategy_candidates[]` are one best routine per the 11 catalogued routine
    families; reporter handoffs are not treated as peer strategy families.
  - Each candidate's estimate separates hands-on work, unattended work and
    elapsed time as minimum/typical/maximum scenarios. Only
    `restriction.digest_ligate_extract_sticky` has a synthetic bench-step
    model. Other routines retain labelled legacy elapsed scalars, with all
    three breakdown ranges unknown, not zero. The text report uses the same
    labels as Routine Assistant; JSON retains dependencies and provenance.
  - Ranking and `max_time_hours` still use the original family/base scalar
    plus summed delays, not the displayed critical path. Inspect
    `planning_estimate.explanation.time_policy` (or
    `planning.estimate.explanation.time_policy` in `routines explain`) for its
    hours, legacy basis and consumers. Text and
    comparison explicitly label this policy proxy; old reports without it
    stay unknown. Bench-model coverage/duration cannot change scores, family
    winners or guardrails. Neither the proxy nor a pass is validated lab time.
  - Cost is split into stock-drawn material value, one-unit missing-material
    purchasing proxy, machine-run cost and routine-base cost. The four sum to
    the retained total, with unknown/unpriced classes explicitly listed.
    These are not package prices, orderable quotes or quantity-aware stock
    accounting. Machine occupancy retains the old base-duration approximation.
  - `routines list` exposes the same fields and full `planning_estimate`;
    `routines explain restriction.digest_ligate_extract_sticky` shows each
    step, delay and critical-path scenario. The model permits dependency-only
    overlap, not staffing/resource scheduling; times are synthetic defaults,
    not validated wet-lab performance or a completion date. See the
    [estimate contract](../protocol.md) and
    [asset provenance](../../assets/cloning_wetlab_steps.README.md).
  - `vector_candidates[]` use structured helper/vector catalog fields only;
    marker, promoter/expression, host-target, and MCS/site constraints remain
    explicit `missing_questions[]` until structured profile fields exist.
  - `--seq-id` is traceability-only in v1 and is not consumed as a
    construct-candidate graph input.
- `planning protein-expression-handoff [--seq-id SEQ_ID] [--objective JSON_OR_@FILE] [--requirements JSON_OR_@FILE] [--provider PROVIDER] [--profile-scope effective] [--format json|text]`
  - Returns a read-only protein-expression planning handoff for underspecified
    requests such as "give me the maximal amount of protein."
  - Emits `product_definition`, `product_readiness`, `sequence_context`,
    `cds_assessment`, `tag_assessment`, optional reviewed `requirements`, ranked `host_chassis_candidates[]`,
    ranked `vector_route_candidates[]`, explicit `missing_questions[]`,
    `neutral_service_draft`, optional `provider_projection`, typed
    `provider_binding`, and readiness-driven `suggested_next_actions[]`.
    The draft uses `gentle.external_service_delivery_route_request.v1`, with
    exact project `seq_id` or `protein_seq_id`, bp/aa length and readiness.
    Missing, unusable or ambiguous products have no draft; no bundled example
    replaces them. Draft source context includes CDS evidence and explicitly
    unresolved coding DNA for protein-only targets.
  - `--requirements` accepts inline JSON, a JSON path, or `@FILE` using schema
    `gentle.protein_expression_requirements.v1`. Reviewed topic records resolve
    only their matching missing questions; partial requirements remain partial.
    Supplied topics populate neutral optimization, tag-policy and delivery
    fields (including yield in the optimization target); metadata retains
    service intent, review provenance, nonclaims and unanswered question IDs.
    Missing fields stay absent; supplied false values and reviewed empty lists
    remain intact. Annotated tags are evidence, not a manufactured vector.
  - Without `--provider`, `provider_binding.status = not_selected`,
    `provider_projection` is null and `service_handoff_candidates[]` is empty.
    No provider-specific preflight/quote command appears. Explicit selection
    produces an existing `gentle.external_service_request.v1` projection only
    if outsourcing is allowed, the product is usable/unambiguous and the active
    catalog supports `protein_expression`; a non-empty `approved_providers`
    must include the selection. Otherwise binding is `withheld`, with typed
    reasons, and no provider actions. In-house actions remain available.
    A usable neutral draft with outsourcing permission suggests reviewing
    `services providers list` before explicitly selecting a provider; no default
    provider is chosen. Quote suggestions additionally require an eligible,
    quote-capable canonical local preflight, not merely complete answers.
  - If `--seq-id` is supplied, the route now performs read-only product-context
    analysis: annotated CDS features are summarized first, whole-sequence CDS
    sanity is used only as an explicit fallback, and the report records
    nucleotide/protein length, GC percent/range, ambiguous bases,
    translation possibility, start/stop sanity, internal stops, and annotated
    tag/signal context.
  - If no usable CDS/protein context is found, `missing_questions[]` asks for
    coding sequence, ORF, CDS annotation, or target-protein boundaries. If a
    product is inferable, the questions shift toward total vs soluble/active/
    purified/secreted yield, purification endpoint, tag preference,
    host/chassis, toxicity/induction, PTMs/cofactors, localization/secretion,
    scale, and delivery endpoint.
  - `product_definition.readiness.status` drives the next actions:
    `whole_sequence_cds_candidate` and `annotated_cds_review_required` point to
    cloning consultation and, only after permitted provider selection, local
    preflight and question-resolved quote-packet preparation;
    `protein_sequence_review_required` points to reverse
    translation or provider protein-target handoff review; `needs_cds_boundary`
    points to CDS/ORF boundary inspection.
  - Text output derives provider lines from the same structured action rows,
    and product-specific preflight/quote actions serialize the same request
    preview instead of pointing at the tutorial protein.
  - Partial requirements leave their questions visible even when a selected
    provider permits a draft/local preflight; they withhold quote suggestions.
    The bundled GeneArt catalog keeps direct submission unavailable.
  - The route never designs, optimizes, submits, orders, writes files, changes
    project state, queries live providers, or promises wet-lab yield.
  - `ui open|focus|close protein-expression-handoff` addresses the read-only
    Services inspector (no arguments). Headless hosts record the intent only;
    use this `planning` command to obtain the actual report. The inspector runs
    the same report on a background snapshot, preserves explicit product/provider
    selection and exposes suggested commands for copying, not execution.
  - External agents can call the named MCP `protein_expression_handoff` tool
    with the same explicit inputs; `requirements`/`objective` are inline JSON
    objects instead of `@FILE` operands. Save the project first when comparing
    with MCP's `state_path`. The inner assistant proposes this shared command;
    neither agent path selects a provider or runs follow-up actions implicitly.
- `planning profile show [--scope global|project_override|confirmed_agent_overlay|effective]`
  - Shows one profile scope or merged effective profile.
- `planning profile set JSON_OR_@FILE [--scope global|project_override|confirmed_agent_overlay]`
  - Sets/replaces one profile scope.
- `planning profile clear [--scope global|project_override|confirmed_agent_overlay]`
  - Clears one profile scope.
- `planning objective show`
  - Shows current objective weights/guardrails.
  - Objective payloads may also include optional `helper_profile_id` and
    `preferred_routine_families[]` fields for helper-aware routine ranking.
- `planning objective set JSON_OR_@FILE`
  - Sets/replaces current planning objective.
- `planning objective clear`
  - Clears objective and falls back to engine defaults.
- `planning suggestions list [--status pending|accepted|rejected]`
  - Lists planning sync suggestions.
- `planning suggestions accept SUGGESTION_ID`
  - Accepts a pending suggestion and applies its patches.
- `planning suggestions reject SUGGESTION_ID [--reason TEXT]`
  - Rejects a pending suggestion.
- `planning sync status`
  - Shows sync metadata and pending count.
- `planning sync pull JSON_OR_@FILE [--source ID] [--confidence N] [--snapshot-id ID]`
  - Registers a pending inbound advisory suggestion.
- `planning sync push JSON_OR_@FILE [--source ID] [--confidence N] [--snapshot-id ID]`
  - Registers a pending outbound advisory suggestion.
- Sync payload shape (`pull`/`push`):
  - optional `profile_patch` (`gentle.planning_profile.v1`)
  - optional `objective_patch` (`gentle.planning_objective.v1`)
  - optional `message`

Planning estimate rule (v1 purchasing simplification):

- Missing required local material adds procurement delay using
  `procurement_business_days_default` (default `10`) unless overridden on the
  item (`inventory.<class>.procurement_business_days`).
- Delay is applied once per missing required material class (deduplicated).
- Business-day model in v1: Monday-Friday only (no holiday calendar).
- Business-day delays are converted to `estimated_time_hours` with a
  deterministic weekend-aware factor (`24h * 7/5` per business day).
- When helper/routine preferences are present in the planning objective,
  `routines list` and `routines compare` emit a synthesized
  `routine_preference_context` and a transparent
  `routine_family_alignment_bonus` inside each planning estimate explanation.
- When `--seq-id` is present and the sequence has construct reasoning with
  variant assay suggestions, that same context also contributes:
  - `construct_reasoning_seq_id`
  - `variant_effect_tags`
  - `variant_suggested_assay_ids`
  - `variant_derived_preferred_routine_families`
  - `routine_family_alignment_sources += variant_derived`

Shipped starter assets:

- Legacy pack:
  - `assets/cloning_patterns.json` (`gentle.cloning_patterns.v1`)
- Hierarchical catalog (one template file per macro):
  - `assets/cloning_patterns_catalog/**/*.json`
  - each file schema: `gentle.cloning_pattern_template.v1`
  - folder hierarchy is used by GUI `Patterns` menu hierarchy
- Typed routine manifest:
  - `assets/cloning_routines.json`
  - schema: `gentle.cloning_routines.v1`
  - includes `gibson.two_fragment_overlap_preview`
    -> `assets/cloning_patterns_catalog/gibson/overlap_assembly/gibson_two_fragment_overlap_preview.json`
- Import commands:
  - `gentle_cli shell 'macros template-import assets/cloning_patterns.json'`
  - `gentle_cli shell 'macros template-import assets/cloning_patterns_catalog'`

<!-- manual-page-nav:start -->

---

[CLI Manual](../cli.md) | [Previous: Restriction Cloning](restriction-cloning.md) | [Next: Candidate Sets](candidates.md)
<!-- manual-page-nav:end -->
