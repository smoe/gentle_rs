# Construct Reasoning Graph Plan

Last updated: 2026-04-09

Purpose: define an engine-owned, sequence-linked construct-planning layer that
captures why GENtle proposes one genomic window, transcript/CDS span, fusion
boundary, or cloning route over another.

This plan is deliberately architecture-first:

- keep reasoning inspectable and editable by humans,
- keep sequence-linked evidence visible on the map,
- keep routines/protocol execution deterministic,
- keep ClawBio/OpenClaw as thin orchestration over GENtle rather than a second
  biology engine.

## 1. Current gap

GENtle already knows many isolated facts and can already execute many
deterministic transforms:

- restriction-site geometry and exact cut positions,
- transcript/exon/CDS derivation,
- cDNA/RNA-read support summaries,
- UniProt projection and isoform architecture inspection,
- TFBS hits,
- Gibson/Golden Gate/restriction planning primitives,
- routine comparison and planning metadata.

What is missing is a first-class shared layer that answers:

- which exact region should be taken for a target construct,
- how much surrounding context should be preserved,
- which neighboring parts should be fused next,
- which observations are structural facts versus context-sensitive evidence,
- how competing alternatives should be compared and explained.

## 2. Core idea

Introduce an engine-owned, sequence-linked reasoning graph:

- sequence spans are evidence-bearing nodes,
- decision nodes consume evidence and produce derived facts,
- derived facts support construct candidates,
- construct candidates feed routine selection and deterministic execution.

The graph is not a black box. It is intended to be:

- human-inspectable,
- human-editable,
- portable across GUI/CLI/JS/Lua/MCP/ClawBio,
- storable in project metadata and run bundles,
- usable offline.

## 3. Evidence model

Do not collapse everything into one confidence score.

Separate at least four classes:

1. `hard_fact`
   - exact sequence-true observations or deterministic geometry.
   - examples:
     - restriction-site positions,
     - exact cut geometry,
     - strand,
     - start/stop codon positions,
     - cDNA-confirmed exon boundaries,
     - explicitly derived CDS frame geometry.
2. `reliable_annotation`
   - high-trust sequence annotation that is not universally context-true.
   - examples:
     - annotated exons without direct tissue-usage evidence,
     - imported CDS/transcript annotations,
     - manually curated feature intervals.
3. `context_evidence`
   - context-dependent support that can strengthen or weaken a role under a
     species/cell/tissue/assay objective.
   - examples:
     - exon usage frequency in one tissue,
     - host-expression suitability,
     - organelle/translation-table context,
     - conservation or assay support.
4. `soft_hypothesis`
   - plausible but intrinsically fuzzy signals.
   - examples:
     - TFBS hits,
     - promoter-like windows,
     - learned regulatory-role predictions,
     - “compact but likely still functional” span suggestions.

This distinction is central:

- a cDNA-supported exon boundary can be a hard fact,
- the same exon’s relevance in one tissue is context evidence,
- TFBS hits remain soft hypotheses unless independently validated.

## 4. Proposed contract surface

Add new portable records in `crates/gentle-protocol`:

- `gentle.construct_objective.v1`
- `gentle.design_evidence.v1`
- `gentle.design_fact.v1`
- `gentle.design_decision_node.v1`
- `gentle.construct_candidate.v1`
- `gentle.construct_reasoning_graph.v1`
- `gentle.construct_reasoning_store.v1`

Suggested enums:

- `ConstructRole`
  - `promoter`
  - `enhancer`
  - `utr_5prime`
  - `cds`
  - `utr_3prime`
  - `terminator`
  - `linker`
  - `tag`
  - `signal_peptide`
  - `localization_signal`
  - `homology_arm`
  - `fusion_boundary`
  - `restriction_site`
  - `splice_boundary`
  - `tfbs`
  - `context_baggage`
  - `other`
- `EvidenceClass`
  - `hard_fact`
  - `reliable_annotation`
  - `context_evidence`
  - `soft_hypothesis`
  - `user_override`
- `DecisionMethod`
  - `hard_rule`
  - `weighted_rule`
  - `fuzzy_rule`
  - `model_summary`
  - `pareto_rank`
  - `user_override`
- `EditableStatus`
  - `draft`
  - `accepted`
  - `rejected`
  - `locked`

### Rust skeletons

```rust
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct ConstructObjective {
    pub schema: String,
    pub objective_id: String,
    pub title: String,
    pub goal: String,
    pub host_species: Option<String>,
    pub cell_type: Option<String>,
    pub tissue: Option<String>,
    pub organelle: Option<String>,
    pub expression_intent: Option<String>,
    pub required_roles: Vec<ConstructRole>,
    pub forbidden_roles: Vec<ConstructRole>,
    pub preferred_routine_families: Vec<String>,
    pub notes: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct DesignEvidence {
    pub evidence_id: String,
    pub seq_id: SeqId,
    pub start_0based: usize,
    pub end_0based_exclusive: usize,
    pub strand: Option<String>,
    pub role: ConstructRole,
    pub evidence_class: EvidenceClass,
    pub label: String,
    pub rationale: String,
    pub score: Option<f64>,
    pub confidence: Option<f64>,
    pub specificity_bias: Option<f64>,
    pub sensitivity_bias: Option<f64>,
    pub context_tags: Vec<String>,
    pub provenance_kind: String,
    pub provenance_refs: Vec<String>,
    pub editable_status: EditableStatus,
    pub warnings: Vec<String>,
    pub notes: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct DesignFact {
    pub fact_id: String,
    pub fact_type: String,
    pub label: String,
    pub rationale: String,
    pub based_on_evidence_ids: Vec<String>,
    pub value_json: serde_json::Value,
    pub confidence: Option<f64>,
    pub editable_status: EditableStatus,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct DesignDecisionNode {
    pub decision_id: String,
    pub decision_type: String,
    pub method: DecisionMethod,
    pub title: String,
    pub rationale: String,
    pub input_evidence_ids: Vec<String>,
    pub input_fact_ids: Vec<String>,
    pub output_fact_ids: Vec<String>,
    pub output_candidate_ids: Vec<String>,
    pub parameters_json: serde_json::Value,
    pub editable_status: EditableStatus,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct ConstructCandidate {
    pub candidate_id: String,
    pub objective_id: String,
    pub title: String,
    pub component_ids: Vec<String>,
    pub derived_from_fact_ids: Vec<String>,
    pub suggested_routine_ids: Vec<String>,
    pub compactness_score: Option<f64>,
    pub confidence_score: Option<f64>,
    pub cloning_complexity_score: Option<f64>,
    pub host_fit_score: Option<f64>,
    pub warnings: Vec<String>,
    pub notes: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct ConstructReasoningGraph {
    pub schema: String,
    pub graph_id: String,
    pub objective: ConstructObjective,
    pub evidence: Vec<DesignEvidence>,
    pub facts: Vec<DesignFact>,
    pub decisions: Vec<DesignDecisionNode>,
    pub candidates: Vec<ConstructCandidate>,
}
```

## 5. File placement

### Protocol crate

- `crates/gentle-protocol/src/construct_reasoning.rs`
- `crates/gentle-protocol/src/lib.rs`

Responsibility:

- portable schemas only,
- no biology execution logic,
- no GUI-specific rendering assumptions.

### Engine

- `src/engine/analysis/construct_reasoning.rs`
- `src/engine/state/construct_reasoning.rs`
- `src/engine/protocol.rs`
- `src/engine.rs`

Responsibility:

- normalize/store graph payloads,
- derive evidence from existing operations/reports/features,
- build decision nodes and construct candidates,
- export stable reports and graph payloads.

### GUI

- `src/main_area_dna.rs`
- follow-up helper module if needed:
  `src/main_area_dna/construct_reasoning.rs`
- optional specialist window in `src/app.rs`

Responsibility:

- sequence-linked overlays,
- graph inspector/editor surface,
- accept/reject/lock and boundary-edit actions,
- no GUI-local biology scoring forks.

### Documentation

- this file:
  `docs/construct_reasoning_plan.md`
- roadmap pointer in `docs/roadmap.md`
- short durable architecture constraints in `docs/architecture.md`

## 6. First evidence extractors

V1 should reuse existing GENtle outputs before introducing any external model:

1. restriction enzymes
   - source: existing restriction-site computation
   - classification: `hard_fact`
2. transcript/exon/CDS derivation
   - source: transcript derivation and splicing/isoform helpers
   - classification:
     - cDNA-confirmed exon bounds: `hard_fact`
     - imported/derived exon/CDS annotations: `reliable_annotation`
3. cDNA/RNA-read support
   - source: RNA-read reports and transcript support summaries
   - classification:
     - confirmed exon path support: `hard_fact` or strong `context_evidence`
     - tissue/abundance usage support: `context_evidence`
4. UniProt projection / isoform architecture
   - source: UniProt mapping and protein-architecture views
   - classification: `context_evidence`
5. TFBS
   - source: existing TFBS hits
   - classification: `soft_hypothesis`

## 7. Decision-node model

Not every node needs fuzzy logic.

Decision nodes should initially support:

- `choose_transcript_or_isoform`
- `choose_cds_span`
- `choose_promoter_window`
- `evaluate_context_baggage`
- `evaluate_fusion_boundary`
- `derive_construct_candidate`
- `suggest_routine_family`

Fuzzy logic belongs only where soft thresholds are meaningful, for example:

- promoter-window quality,
- compactness vs baggage tradeoff,
- fusion compatibility,
- host/cell expression suitability.

Do not use fuzzy logic for:

- exact restriction-site positions,
- cDNA-confirmed exon boundaries,
- strand,
- CDS frame geometry.

## 8. Human inspectability and editability

This layer must remain editable by humans.

Required edit operations:

- accept one evidence node,
- reject one evidence node,
- lock one fact,
- adjust one span boundary,
- add one manual evidence/fact node,
- attach a user note,
- promote one candidate as preferred.

All overrides should remain explicit and portable.

Suggested additive operations:

- `UpsertConstructObjective`
- `BuildConstructReasoningGraph`
- `AcceptDesignEvidence`
- `RejectDesignEvidence`
- `LockDesignFact`
- `AdjustDesignEvidenceSpan`
- `DeriveConstructCandidates`
- `ExportConstructReasoningJson`
- `RenderConstructReasoningSvg`

## 9. Sequence-linked UI model

The DNA window should become one surface of the graph.

Sequence overlay rules:

- color hue = role,
- border style = evidence class,
- fill opacity = confidence/weight,
- solid border = accepted or hard fact,
- dashed border = soft hypothesis,
- muted/struck style = rejected,
- lock badge = human-locked fact.

Interaction rules:

- clicking a sequence span highlights connected decision nodes,
- clicking a decision node highlights all contributing sequence spans,
- selecting a candidate highlights its contributing spans and decisions,
- accepted/rejected/locked edits immediately update the same shared graph
  payload, not a GUI-only shadow state.

## 10. Routine Assistant relationship

Do not force protocol choice before construct reasoning.

Preferred order:

1. define objective,
2. build reasoning graph,
3. derive construct candidates,
4. compare/suggest routines,
5. preflight and execute chosen routine.

This keeps construct intent separate from protocol execution while allowing
Routine Assistant to consume richer candidate context.

## 11. Offline-first policy

V1 should not require internet or external model calls.

- graphs, evidence, and overrides live in project metadata,
- export/import should work offline,
- model-derived evidence is additive and optional later,
- ClawBio/OpenClaw can inspect or edit the same graph payloads offline from
  saved state/run-bundle artifacts.

## 12. Delivery plan

### Phase 1: protocol contracts + metadata store

- add portable schemas in `crates/gentle-protocol`
- add engine metadata store + normalization helpers
- re-export through `src/engine/protocol.rs`

Acceptance:

- empty graph/objective round-trips deterministically
- metadata load/save survives project persistence

### Phase 2: read-only deterministic evidence extraction

- build graph from existing sources:
  restriction sites, transcript/CDS/exon derivation, RNA-read support,
  UniProt projection, TFBS
- no user edits yet

Acceptance:

- deterministic graph build for one test construct
- evidence classes are stable and reproducible

### Phase 3: sequence-linked visualization

- add DNA-window overlay
- add read-only graph/detail inspector

Acceptance:

- selecting evidence on sequence highlights connected nodes
- one SVG/JSON export path exists

### Phase 4: editable reasoning graph

- accept/reject/lock
- manual note support
- boundary adjustment for span-backed evidence

Acceptance:

- edits persist in metadata
- reload reproduces the same graph state

### Phase 5: construct candidates + routine handoff

- derive `construct_candidate` rows
- feed candidates into Routine Assistant compare/rank stages

Acceptance:

- at least one objective produces multiple candidates with explainable
  tradeoffs
- routine suggestion is derived from construct candidates, not guessed ad hoc

### Phase 6: optional learned/model-backed evidence

- add optional sidecar/model adapters
- keep outputs normalized into the same graph contracts

Acceptance:

- model-backed evidence remains additive and inspectable
- offline deterministic baseline remains intact without the model

## 13. Minimal v1 milestone

The smallest useful milestone is:

- one `construct_objective`,
- one `construct_reasoning_graph`,
- evidence extractors for:
  - restriction sites,
  - exon/CDS/transcript spans,
  - TFBS,
- one read-only DNA overlay,
- one JSON export/import path.

That is enough to prove:

- sequence-linked reasoning is useful,
- hard facts and fuzzy hypotheses can coexist without being conflated,
- the graph can later drive construct-candidate derivation and routine choice.

## 14. Test matrix

- protocol/schema tests:
  - deterministic serde round-trip for all new records
- engine tests:
  - graph normalization/order stability
  - evidence-class assignment from existing sources
  - metadata persistence and reload
- GUI tests:
  - overlay mapping from graph payload to visible sequence regions
  - selection/highlight linkage between span and decision/fact panes
- export tests:
  - deterministic JSON export
  - deterministic SVG export once rendering exists

## 15. Non-goals for v1

- no hidden autonomous construct designer,
- no model-required online execution path,
- no replacement of existing deterministic cloning operations,
- no adapter-only reasoning layers in GUI/ClawBio,
- no premature single-score ranking that destroys tradeoff visibility.
