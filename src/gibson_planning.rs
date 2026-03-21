//! Gibson assembly planning contracts and deterministic preview derivation.
//!
//! This module promotes `gentle.gibson_assembly_plan.v1` from a docs-only draft
//! into an engine-consumed planning artifact for a restricted v1 scope:
//! destination-first, single-insert Gibson planning.
//! The preview result is non-mutating and is intended to drive GUI, shell/CLI,
//! and protocol-cartoon rendering from one shared derivation path.

use crate::{
    enzymes::active_restriction_enzymes,
    engine::{EngineError, ErrorCode, GentleEngine},
    feature_location::collect_location_ranges_usize,
    protocol_cartoon::{
        DnaEndStyle, OverhangPolarity, ProtocolCartoonKind, ProtocolCartoonTemplateBindings,
        ProtocolCartoonTemplateEventBinding, ProtocolCartoonTemplateFeatureBinding,
        ProtocolCartoonTemplateMoleculeBinding,
    },
    restriction_enzyme::RestrictionEnzymeSite,
};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, HashMap, HashSet};

pub const GIBSON_ASSEMBLY_PLAN_SCHEMA: &str = "gentle.gibson_assembly_plan.v1";
pub const GIBSON_ASSEMBLY_PREVIEW_SCHEMA: &str = "gentle.gibson_assembly_preview.v1";

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct GibsonAssemblyPlan {
    pub schema: String,
    pub id: String,
    pub title: String,
    pub summary: String,
    pub destination: GibsonPlanDestination,
    pub product: GibsonPlanProduct,
    pub fragments: Vec<GibsonPlanFragment>,
    pub assembly_order: Vec<GibsonPlanAssemblyMember>,
    pub junctions: Vec<GibsonPlanJunction>,
    pub validation_policy: GibsonPlanValidationPolicy,
    #[serde(default)]
    pub derived_design: Option<serde_json::Value>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct GibsonPlanDestination {
    pub seq_id: String,
    pub topology_before_opening: String,
    pub opening: GibsonPlanOpening,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct GibsonPlanOpening {
    pub mode: String,
    pub label: String,
    pub start_0based: Option<usize>,
    pub end_0based_exclusive: Option<usize>,
    pub left_end_id: String,
    pub right_end_id: String,
    pub uniqueness_requirement: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct GibsonPlanProduct {
    pub topology: String,
    pub output_id_hint: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct GibsonPlanFragment {
    pub id: String,
    pub seq_id: String,
    pub role: String,
    pub orientation: String,
    pub left_end_strategy: Option<GibsonPlanEndStrategy>,
    pub right_end_strategy: Option<GibsonPlanEndStrategy>,
    #[serde(default)]
    pub source_span_1based: Option<GibsonPlanSourceSpan>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct GibsonPlanSourceSpan {
    pub source_seq_id: String,
    pub start: usize,
    pub end: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct GibsonPlanEndStrategy {
    pub mode: String,
    pub target_junction_id: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct GibsonPlanAssemblyMember {
    pub kind: String,
    pub id: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct GibsonPlanJunction {
    pub id: String,
    pub left_member: GibsonPlanAssemblyMember,
    pub right_member: GibsonPlanAssemblyMember,
    pub required_overlap_bp: Option<usize>,
    pub overlap_partition: Option<GibsonPlanOverlapPartition>,
    pub overlap_source: String,
    #[serde(default)]
    pub distinct_from: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct GibsonPlanOverlapPartition {
    pub left_member_bp: usize,
    pub right_member_bp: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct GibsonPlanValidationPolicy {
    pub require_unambiguous_destination_opening: bool,
    pub require_distinct_terminal_junctions: bool,
    pub adjacency_overlap_mismatch: String,
    pub design_targets: GibsonPlanDesignTargets,
    pub uniqueness_checks: GibsonPlanUniquenessChecks,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default)]
pub struct GibsonPlanDesignTargets {
    pub overlap_bp_min: usize,
    pub overlap_bp_max: usize,
    pub minimum_overlap_tm_celsius: f64,
    pub priming_segment_tm_min_celsius: f64,
    pub priming_segment_tm_max_celsius: f64,
    pub priming_segment_min_length_bp: usize,
    pub priming_segment_max_length_bp: usize,
    pub max_anneal_hits: usize,
}

impl Default for GibsonPlanDesignTargets {
    fn default() -> Self {
        Self {
            overlap_bp_min: 20,
            overlap_bp_max: 40,
            minimum_overlap_tm_celsius: 60.0,
            priming_segment_tm_min_celsius: 58.0,
            priming_segment_tm_max_celsius: 68.0,
            priming_segment_min_length_bp: 18,
            priming_segment_max_length_bp: 35,
            max_anneal_hits: 1,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct GibsonPlanUniquenessChecks {
    pub destination_context: String,
    pub participating_fragments: String,
    #[serde(default)]
    pub reference_contexts: Vec<GibsonPlanReferenceContext>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct GibsonPlanReferenceContext {
    pub reference_id: String,
    pub severity: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct GibsonDestinationOpeningSuggestion {
    pub kind: String,
    pub label: String,
    pub summary: String,
    pub start_0based: usize,
    pub end_0based_exclusive: usize,
    pub enzyme_name: Option<String>,
    pub end_geometry: String,
    pub overhang_bp: Option<usize>,
    pub in_mcs_context: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct GibsonAssemblyPreview {
    pub schema: String,
    pub plan_id: String,
    pub title: String,
    pub summary: String,
    pub can_execute: bool,
    pub destination: GibsonPreviewDestination,
    pub insert: GibsonPreviewInsert,
    #[serde(default)]
    pub resolved_junctions: Vec<GibsonResolvedJunctionPreview>,
    #[serde(default)]
    pub primer_suggestions: Vec<GibsonPrimerSuggestion>,
    #[serde(default)]
    pub warnings: Vec<String>,
    #[serde(default)]
    pub errors: Vec<String>,
    #[serde(default)]
    pub notes: Vec<String>,
    pub cartoon: GibsonCartoonPreview,
    pub routine_handoff: GibsonRoutineHandoffPreview,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct GibsonPreviewDestination {
    pub seq_id: String,
    pub topology_before_opening: String,
    pub actual_topology: String,
    pub length_bp: usize,
    pub opening_mode: String,
    pub opening_label: String,
    pub opening_start_0based: Option<usize>,
    pub opening_end_0based_exclusive: Option<usize>,
    pub removed_span_bp: Option<usize>,
    pub left_end_id: String,
    pub right_end_id: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct GibsonPreviewInsert {
    pub fragment_id: String,
    pub seq_id: String,
    pub role: String,
    pub orientation: String,
    pub length_bp: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct GibsonResolvedJunctionPreview {
    pub junction_id: String,
    pub left_member_id: String,
    pub right_member_id: String,
    pub overlap_bp: usize,
    pub overlap_tm_celsius: f64,
    pub overlap_sequence: String,
    pub overlap_source: String,
    pub distinct_from: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct GibsonPrimerSuggestion {
    pub primer_id: String,
    pub side: String,
    pub fragment_id: String,
    pub template_seq_id: String,
    pub full_sequence: String,
    pub overlap_5prime: GibsonPrimerSegment,
    pub priming_3prime: GibsonPrimerSegment,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct GibsonPrimerSegment {
    pub sequence: String,
    pub length_bp: usize,
    pub tm_celsius: f64,
    pub gc_fraction: f64,
    pub anneal_hits: usize,
    pub note: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default)]
pub struct GibsonCartoonPreview {
    pub protocol_id: String,
    pub representative_overlap_bp: usize,
    pub bindings: ProtocolCartoonTemplateBindings,
}

impl Default for GibsonCartoonPreview {
    fn default() -> Self {
        Self {
            protocol_id: String::new(),
            representative_overlap_bp: 0,
            bindings: empty_protocol_cartoon_bindings(),
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct GibsonRoutineHandoffPreview {
    pub supported: bool,
    pub routine_id: String,
    pub reason: Option<String>,
    #[serde(default)]
    pub bindings: BTreeMap<String, String>,
}

#[derive(Debug, Clone)]
struct GibsonResolvedOpening {
    mode: String,
    label: String,
    start_0based: Option<usize>,
    end_0based_exclusive: Option<usize>,
    removed_span_bp: Option<usize>,
    left_end_id: String,
    right_end_id: String,
    destination_is_circular: bool,
}

#[derive(Debug, Clone)]
struct GibsonPrimerCandidate {
    sequence: String,
    tm_celsius: f64,
    gc_fraction: f64,
    anneal_hits: usize,
}

#[derive(Debug, Clone)]
struct GibsonMcsFeatureHint {
    label: String,
    start_0based: usize,
    end_0based_exclusive: usize,
    candidate_enzymes: Vec<String>,
}

impl GentleEngine {
    /// Preview one restricted v1 Gibson assembly plan without mutating engine state.
    pub fn preview_gibson_assembly_plan(
        &self,
        plan: &GibsonAssemblyPlan,
    ) -> Result<GibsonAssemblyPreview, EngineError> {
        preview_gibson_assembly_plan(self, plan)
    }

    /// Suggest destination openings for the Gibson specialist UI.
    pub fn suggest_gibson_destination_openings(
        &self,
        seq_id: &str,
    ) -> Result<Vec<GibsonDestinationOpeningSuggestion>, EngineError> {
        suggest_gibson_destination_openings(self, seq_id)
    }
}

pub fn suggest_gibson_destination_openings(
    engine: &GentleEngine,
    seq_id: &str,
) -> Result<Vec<GibsonDestinationOpeningSuggestion>, EngineError> {
    let seq_id = seq_id.trim();
    if seq_id.is_empty() {
        return Ok(vec![]);
    }
    let dna = engine
        .state()
        .sequences
        .get(seq_id)
        .ok_or_else(|| EngineError {
            code: ErrorCode::NotFound,
            message: format!("Destination sequence '{}' not found", seq_id),
        })?;
    let seq_len = dna.len();
    if seq_len == 0 {
        return Ok(vec![]);
    }

    let mut suggestions = vec![];
    let mcs_hints = collect_mcs_feature_hints(dna);
    let unique_sites = unique_forward_sites(dna);

    for hint in &mcs_hints {
        suggestions.push(GibsonDestinationOpeningSuggestion {
            kind: "mcs_feature".to_string(),
            label: format!(
                "{} ({}..{})",
                hint.label,
                hint.start_0based + 1,
                hint.end_0based_exclusive
            ),
            summary: format!(
                "Use the annotated MCS feature span as the destination opening window; fills defined-site {}..{} (0-based).",
                hint.start_0based, hint.end_0based_exclusive
            ),
            start_0based: hint.start_0based,
            end_0based_exclusive: hint.end_0based_exclusive,
            enzyme_name: None,
            end_geometry: "feature_span".to_string(),
            overhang_bp: None,
            in_mcs_context: true,
        });
    }

    let mut seen_enzymes = HashSet::new();
    let mut mcs_site_suggestions = vec![];
    for hint in &mcs_hints {
        for enzyme_name in &hint.candidate_enzymes {
            let Some(site) = unique_sites.get(enzyme_name) else {
                continue;
            };
            if !seen_enzymes.insert(enzyme_name.clone()) {
                continue;
            }
            if let Some(suggestion) = restriction_site_suggestion(site, true, seq_len) {
                mcs_site_suggestions.push(suggestion);
            }
        }
    }
    mcs_site_suggestions.sort_by(|left, right| {
        geometry_priority(&left.end_geometry)
            .cmp(&geometry_priority(&right.end_geometry))
            .then(left.start_0based.cmp(&right.start_0based))
            .then(left.label.cmp(&right.label))
    });
    if !mcs_site_suggestions.is_empty() {
        suggestions.extend(mcs_site_suggestions);
        return Ok(suggestions);
    }

    let mut fallback_unique_sites = unique_sites
        .values()
        .filter_map(|site| restriction_site_suggestion(site, false, seq_len))
        .collect::<Vec<_>>();
    fallback_unique_sites.sort_by(|left, right| {
        geometry_priority(&left.end_geometry)
            .cmp(&geometry_priority(&right.end_geometry))
            .then(left.start_0based.cmp(&right.start_0based))
            .then(left.label.cmp(&right.label))
    });
    suggestions.extend(fallback_unique_sites);
    Ok(suggestions)
}

fn unique_forward_sites<'a>(
    dna: &'a crate::dna_sequence::DNAsequence,
) -> HashMap<String, &'a RestrictionEnzymeSite> {
    let mut counts: HashMap<String, usize> = HashMap::new();
    let mut first_sites: HashMap<String, &RestrictionEnzymeSite> = HashMap::new();
    for site in dna.restriction_enzyme_sites().iter().filter(|site| site.forward_strand) {
        let name = site.enzyme.name.clone();
        *counts.entry(name.clone()).or_insert(0) += 1;
        first_sites.entry(name).or_insert(site);
    }
    first_sites
        .into_iter()
        .filter_map(|(name, site)| (counts.get(&name).copied() == Some(1)).then_some((name, site)))
        .collect()
}

fn restriction_site_suggestion(
    site: &RestrictionEnzymeSite,
    in_mcs_context: bool,
    seq_len: usize,
) -> Option<GibsonDestinationOpeningSuggestion> {
    let start_0based = usize::try_from(site.offset).ok()?;
    let end_0based_exclusive = start_0based.checked_add(site.enzyme.sequence.len())?;
    if end_0based_exclusive > seq_len {
        return None;
    }
    let (end_geometry, overhang_bp, geometry_label) =
        restriction_end_geometry_descriptor(site.enzyme.overlap);
    let context_label = if in_mcs_context {
        "MCS-linked unique cutter"
    } else {
        "Unique cutter"
    };
    Some(GibsonDestinationOpeningSuggestion {
        kind: "unique_restriction_site".to_string(),
        label: format!(
            "{} | {} | {}..{}",
            site.enzyme.name,
            geometry_label,
            start_0based + 1,
            end_0based_exclusive
        ),
        summary: format!(
            "{} {} using recognition span {}..{} (1-based); fills defined-site {}..{} (0-based).",
            context_label,
            site.enzyme.name,
            start_0based + 1,
            end_0based_exclusive,
            start_0based,
            end_0based_exclusive
        ),
        start_0based,
        end_0based_exclusive,
        enzyme_name: Some(site.enzyme.name.clone()),
        end_geometry: end_geometry.to_string(),
        overhang_bp,
        in_mcs_context,
    })
}

fn restriction_end_geometry_descriptor(overlap: isize) -> (&'static str, Option<usize>, String) {
    if overlap == 0 {
        ("blunt", None, "blunt".to_string())
    } else if overlap > 0 {
        let bp = overlap as usize;
        (
            "5prime_overhang",
            Some(bp),
            format!("5' overhang ({} bp)", bp),
        )
    } else {
        let bp = overlap.unsigned_abs();
        (
            "3prime_overhang",
            Some(bp),
            format!("3' overhang ({} bp)", bp),
        )
    }
}

fn geometry_priority(end_geometry: &str) -> usize {
    match end_geometry {
        "blunt" => 0,
        "5prime_overhang" => 1,
        "3prime_overhang" => 2,
        "feature_span" => 3,
        _ => 4,
    }
}

fn collect_mcs_feature_hints(dna: &crate::dna_sequence::DNAsequence) -> Vec<GibsonMcsFeatureHint> {
    let lookup = rebase_name_lookup_by_normalized();
    let available_sites: HashSet<String> = dna
        .restriction_enzyme_sites()
        .iter()
        .filter(|site| site.forward_strand)
        .map(|site| site.enzyme.name.clone())
        .collect();
    let mut hints = vec![];
    for feature in dna.features() {
        let mcs_hint = feature
            .qualifier_values("mcs_expected_sites".into())
            .next()
            .is_some()
            || feature
                .qualifier_values("mcs_preset".into())
                .next()
                .is_some()
            || feature_first_nonempty_qualifier(
                feature,
                &["label", "note", "gene", "name", "standard_name"],
            )
            .map(|text| text_mentions_mcs(&text))
            .unwrap_or(false);
        if !mcs_hint {
            continue;
        }

        let mut ranges = vec![];
        collect_location_ranges_usize(&feature.location, &mut ranges);
        if ranges.is_empty() {
            continue;
        }
        let start_0based = ranges.iter().map(|(start, _)| *start).min().unwrap_or(0);
        let end_0based_exclusive = ranges
            .iter()
            .map(|(_, end)| *end)
            .max()
            .unwrap_or(start_0based);
        if end_0based_exclusive <= start_0based {
            continue;
        }

        let mut candidate_enzymes = vec![];
        let mut seen_candidates = HashSet::new();
        for raw in feature.qualifier_values("mcs_expected_sites".into()) {
            for token in raw
                .split(',')
                .map(str::trim)
                .filter(|value| !value.is_empty())
            {
                if let Some(name) = canonicalize_rebase_enzyme_name(token, &lookup) {
                    if seen_candidates.insert(name.clone()) && available_sites.contains(&name) {
                        candidate_enzymes.push(name);
                    }
                } else {
                    for name in extract_rebase_enzyme_names_from_text(token, &lookup) {
                        if seen_candidates.insert(name.clone()) && available_sites.contains(&name) {
                            candidate_enzymes.push(name);
                        }
                    }
                }
            }
        }
        for key in ["note", "label", "gene", "name", "standard_name"] {
            for raw in feature.qualifier_values(key.into()) {
                for name in extract_rebase_enzyme_names_from_text(raw, &lookup) {
                    if seen_candidates.insert(name.clone()) && available_sites.contains(&name) {
                        candidate_enzymes.push(name);
                    }
                }
            }
        }

        let label = feature_first_nonempty_qualifier(
            feature,
            &["label", "standard_name", "name", "gene", "note"],
        )
        .unwrap_or_else(|| "Multiple Cloning Site (MCS)".to_string());
        hints.push(GibsonMcsFeatureHint {
            label,
            start_0based,
            end_0based_exclusive,
            candidate_enzymes,
        });
    }
    hints.sort_by(|left, right| {
        left.start_0based
            .cmp(&right.start_0based)
            .then(left.end_0based_exclusive.cmp(&right.end_0based_exclusive))
            .then(left.label.cmp(&right.label))
    });
    hints
}

fn normalize_enzyme_match_token(raw: &str) -> String {
    raw.chars()
        .filter(|c| c.is_ascii_alphanumeric())
        .map(|c| c.to_ascii_uppercase())
        .collect()
}

fn rebase_name_lookup_by_normalized() -> HashMap<String, String> {
    let mut lookup = HashMap::new();
    for enzyme in active_restriction_enzymes() {
        let normalized = normalize_enzyme_match_token(&enzyme.name);
        if !normalized.is_empty() {
            lookup.entry(normalized).or_insert(enzyme.name);
        }
    }
    lookup
}

fn canonicalize_rebase_enzyme_name(
    raw: &str,
    lookup: &HashMap<String, String>,
) -> Option<String> {
    let normalized = normalize_enzyme_match_token(raw);
    lookup.get(&normalized).cloned()
}

fn extract_rebase_enzyme_names_from_text(
    text: &str,
    lookup: &HashMap<String, String>,
) -> Vec<String> {
    let tokens = text
        .split(|c: char| !c.is_ascii_alphanumeric())
        .map(str::trim)
        .filter(|value| !value.is_empty())
        .map(str::to_string)
        .collect::<Vec<_>>();
    let mut out = vec![];
    let mut seen = HashSet::new();
    let mut push_candidate = |candidate: String| {
        if let Some(name) = canonicalize_rebase_enzyme_name(&candidate, lookup) {
            if seen.insert(name.clone()) {
                out.push(name);
            }
        }
    };
    for idx in 0..tokens.len() {
        push_candidate(tokens[idx].clone());
        if idx + 1 < tokens.len() {
            push_candidate(format!("{}{}", tokens[idx], tokens[idx + 1]));
        }
        if idx + 2 < tokens.len() {
            push_candidate(format!(
                "{}{}{}",
                tokens[idx],
                tokens[idx + 1],
                tokens[idx + 2]
            ));
        }
    }
    out
}

fn text_mentions_mcs(text: &str) -> bool {
    let lower = text.to_ascii_lowercase();
    lower.contains("multiple cloning site")
        || lower == "mcs"
        || lower.contains(" mcs ")
        || lower.starts_with("mcs ")
        || lower.ends_with(" mcs")
        || lower.contains("(mcs)")
}

fn feature_first_nonempty_qualifier(feature: &gb_io::seq::Feature, keys: &[&str]) -> Option<String> {
    for key in keys {
        for value in feature.qualifier_values((*key).into()) {
            let normalized = value.split_whitespace().collect::<Vec<_>>().join(" ");
            let normalized = normalized.trim().to_string();
            if !normalized.is_empty() {
                return Some(normalized);
            }
        }
    }
    None
}

pub fn preview_gibson_assembly_plan(
    engine: &GentleEngine,
    plan: &GibsonAssemblyPlan,
) -> Result<GibsonAssemblyPreview, EngineError> {
    if plan.schema.trim() != GIBSON_ASSEMBLY_PLAN_SCHEMA {
        return Err(EngineError {
            code: ErrorCode::InvalidInput,
            message: format!(
                "Unsupported Gibson assembly plan schema '{}' (expected '{}')",
                plan.schema, GIBSON_ASSEMBLY_PLAN_SCHEMA
            ),
        });
    }

    let destination_id = plan.destination.seq_id.trim();
    if destination_id.is_empty() {
        return Err(EngineError {
            code: ErrorCode::InvalidInput,
            message: "Gibson plan requires destination.seq_id".to_string(),
        });
    }
    let destination = engine
        .state()
        .sequences
        .get(destination_id)
        .ok_or_else(|| EngineError {
            code: ErrorCode::NotFound,
            message: format!("Destination sequence '{}' not found", destination_id),
        })?
        .clone();

    let mut preview = GibsonAssemblyPreview {
        schema: GIBSON_ASSEMBLY_PREVIEW_SCHEMA.to_string(),
        plan_id: plan.id.trim().to_string(),
        title: if plan.title.trim().is_empty() {
            format!("Gibson preview for '{}'", destination_id)
        } else {
            plan.title.trim().to_string()
        },
        summary: plan.summary.trim().to_string(),
        can_execute: false,
        destination: GibsonPreviewDestination {
            seq_id: destination_id.to_string(),
            topology_before_opening: plan.destination.topology_before_opening.trim().to_string(),
            actual_topology: if destination.is_circular() {
                "circular".to_string()
            } else {
                "linear".to_string()
            },
            length_bp: destination.len(),
            opening_mode: plan.destination.opening.mode.trim().to_string(),
            opening_label: plan.destination.opening.label.trim().to_string(),
            opening_start_0based: plan.destination.opening.start_0based,
            opening_end_0based_exclusive: plan.destination.opening.end_0based_exclusive,
            removed_span_bp: None,
            left_end_id: plan.destination.opening.left_end_id.trim().to_string(),
            right_end_id: plan.destination.opening.right_end_id.trim().to_string(),
        },
        insert: GibsonPreviewInsert::default(),
        resolved_junctions: vec![],
        primer_suggestions: vec![],
        warnings: vec![],
        errors: vec![],
        notes: vec![],
        cartoon: GibsonCartoonPreview {
            protocol_id: ProtocolCartoonKind::GibsonTwoFragment.id().to_string(),
            representative_overlap_bp: 30,
            bindings: empty_protocol_cartoon_bindings(),
        },
        routine_handoff: GibsonRoutineHandoffPreview::default(),
    };

    if plan.fragments.len() != 1 {
        preview.errors.push(format!(
            "Gibson specialist v1 currently supports exactly one insert fragment, but plan declares {} fragment(s)",
            plan.fragments.len()
        ));
        return finalize_preview(preview);
    }
    let insert_plan = &plan.fragments[0];
    let insert_seq_id = insert_plan.seq_id.trim();
    if insert_seq_id.is_empty() {
        preview
            .errors
            .push("Gibson specialist v1 requires fragments[0].seq_id".to_string());
        return finalize_preview(preview);
    }
    let insert = engine
        .state()
        .sequences
        .get(insert_seq_id)
        .ok_or_else(|| EngineError {
            code: ErrorCode::NotFound,
            message: format!("Insert sequence '{}' not found", insert_seq_id),
        })?
        .clone();
    preview.insert = GibsonPreviewInsert {
        fragment_id: if insert_plan.id.trim().is_empty() {
            "insert".to_string()
        } else {
            insert_plan.id.trim().to_string()
        },
        seq_id: insert_seq_id.to_string(),
        role: if insert_plan.role.trim().is_empty() {
            "insert".to_string()
        } else {
            insert_plan.role.trim().to_string()
        },
        orientation: normalized_orientation(&insert_plan.orientation),
        length_bp: insert.len(),
    };

    validate_destination_topology(&plan.destination, &destination, &mut preview);
    let preview_fragment_id = preview.insert.fragment_id.clone();
    validate_single_insert_shape(plan, &preview_fragment_id, &mut preview);

    let resolved_opening = resolve_opening(plan, &destination, &mut preview);
    let Some(resolved_opening) = resolved_opening else {
        return finalize_preview(preview);
    };
    preview.destination.opening_mode = resolved_opening.mode.clone();
    preview.destination.opening_label = resolved_opening.label.clone();
    preview.destination.opening_start_0based = resolved_opening.start_0based;
    preview.destination.opening_end_0based_exclusive = resolved_opening.end_0based_exclusive;
    preview.destination.removed_span_bp = resolved_opening.removed_span_bp;
    preview.destination.left_end_id = resolved_opening.left_end_id.clone();
    preview.destination.right_end_id = resolved_opening.right_end_id.clone();

    let oriented_insert_seq = if preview.insert.orientation.eq_ignore_ascii_case("reverse") {
        GentleEngine::reverse_complement(&insert.get_forward_string())
    } else {
        insert.get_forward_string().to_ascii_uppercase()
    };
    let insert_template_forward = insert.get_forward_string().to_ascii_uppercase();

    let left_junction = plan
        .junctions
        .iter()
        .find(|junction| junction.id.eq_ignore_ascii_case("junction_left"))
        .or_else(|| {
            plan.junctions.iter().find(|junction| {
                junction
                    .left_member
                    .id
                    .eq_ignore_ascii_case(&resolved_opening.left_end_id)
            })
        })
        .cloned();
    let right_junction = plan
        .junctions
        .iter()
        .find(|junction| junction.id.eq_ignore_ascii_case("junction_right"))
        .or_else(|| {
            plan.junctions.iter().find(|junction| {
                junction
                    .right_member
                    .id
                    .eq_ignore_ascii_case(&resolved_opening.right_end_id)
            })
        })
        .cloned();

    let Some(left_junction) = left_junction else {
        preview
            .errors
            .push("Could not resolve a left terminal junction in the Gibson plan".to_string());
        return finalize_preview(preview);
    };
    let Some(right_junction) = right_junction else {
        preview.errors.push(
            "Could not resolve a right terminal junction in the Gibson plan".to_string(),
        );
        return finalize_preview(preview);
    };

    let left_overlap = derive_terminal_overlap(
        &destination.get_forward_string().to_ascii_uppercase(),
        &resolved_opening,
        TerminalSide::Left,
        &left_junction,
        &plan.validation_policy.design_targets,
        &mut preview,
    );
    let right_overlap = derive_terminal_overlap(
        &destination.get_forward_string().to_ascii_uppercase(),
        &resolved_opening,
        TerminalSide::Right,
        &right_junction,
        &plan.validation_policy.design_targets,
        &mut preview,
    );
    let (Some(left_overlap), Some(right_overlap)) = (left_overlap, right_overlap) else {
        return finalize_preview(preview);
    };

    if plan.validation_policy.require_distinct_terminal_junctions {
        let right_overlap_rc = GentleEngine::reverse_complement(&right_overlap.2);
        if left_overlap.2.eq_ignore_ascii_case(&right_overlap.2)
            || left_overlap.2.eq_ignore_ascii_case(&right_overlap_rc)
        {
            preview.errors.push(format!(
                "Terminal overlap regions are not distinct enough for a destination-first Gibson plan (left='{}', right='{}')",
                left_overlap.2, right_overlap.2
            ));
        }
    }

    preview.resolved_junctions.push(GibsonResolvedJunctionPreview {
        junction_id: left_junction.id.clone(),
        left_member_id: left_junction.left_member.id.clone(),
        right_member_id: left_junction.right_member.id.clone(),
        overlap_bp: left_overlap.1,
        overlap_tm_celsius: GentleEngine::estimate_primer_tm_c(left_overlap.2.as_bytes()),
        overlap_sequence: left_overlap.2.clone(),
        overlap_source: left_junction.overlap_source.clone(),
        distinct_from: left_junction.distinct_from.clone(),
    });
    preview.resolved_junctions.push(GibsonResolvedJunctionPreview {
        junction_id: right_junction.id.clone(),
        left_member_id: right_junction.left_member.id.clone(),
        right_member_id: right_junction.right_member.id.clone(),
        overlap_bp: right_overlap.1,
        overlap_tm_celsius: GentleEngine::estimate_primer_tm_c(right_overlap.2.as_bytes()),
        overlap_sequence: right_overlap.2.clone(),
        overlap_source: right_junction.overlap_source.clone(),
        distinct_from: right_junction.distinct_from.clone(),
    });

    if preview.errors.is_empty() {
        let left_priming = choose_insert_priming_segment(
            &insert_template_forward,
            &oriented_insert_seq,
            TerminalSide::Left,
            &plan.validation_policy.design_targets,
            &mut preview,
        );
        let right_priming = choose_insert_priming_segment(
            &insert_template_forward,
            &oriented_insert_seq,
            TerminalSide::Right,
            &plan.validation_policy.design_targets,
            &mut preview,
        );
        if let (Some(left_priming), Some(right_priming)) = (left_priming, right_priming) {
            let left_full = format!("{}{}", left_overlap.2, left_priming.sequence);
            let right_overlap_primer = GentleEngine::reverse_complement(&right_overlap.2);
            let right_full = format!("{}{}", right_overlap_primer, right_priming.sequence);
            preview.primer_suggestions.push(GibsonPrimerSuggestion {
                primer_id: format!("{}_left_insert_primer", preview.insert.fragment_id),
                side: "left_insert_primer".to_string(),
                fragment_id: preview.insert.fragment_id.clone(),
                template_seq_id: preview.insert.seq_id.clone(),
                full_sequence: left_full,
                overlap_5prime: GibsonPrimerSegment {
                    sequence: left_overlap.2.clone(),
                    length_bp: left_overlap.1,
                    tm_celsius: GentleEngine::estimate_primer_tm_c(left_overlap.2.as_bytes()),
                    gc_fraction: GentleEngine::sequence_gc_fraction(left_overlap.2.as_bytes())
                        .unwrap_or(0.0),
                    anneal_hits: count_hits_both_strands(
                        &destination.get_forward_string().to_ascii_uppercase(),
                        &left_overlap.2,
                    ),
                    note: "5' non-priming overlap segment derived from the left destination flank"
                        .to_string(),
                },
                priming_3prime: GibsonPrimerSegment {
                    note: "3' gene-specific priming segment chosen from the insert terminus to satisfy the Gibson PCR target window".to_string(),
                    ..left_priming
                },
            });
            preview.primer_suggestions.push(GibsonPrimerSuggestion {
                primer_id: format!("{}_right_insert_primer", preview.insert.fragment_id),
                side: "right_insert_primer".to_string(),
                fragment_id: preview.insert.fragment_id.clone(),
                template_seq_id: preview.insert.seq_id.clone(),
                full_sequence: right_full,
                overlap_5prime: GibsonPrimerSegment {
                    sequence: right_overlap_primer,
                    length_bp: right_overlap.1,
                    tm_celsius: GentleEngine::estimate_primer_tm_c(right_overlap.2.as_bytes()),
                    gc_fraction: GentleEngine::sequence_gc_fraction(right_overlap.2.as_bytes())
                        .unwrap_or(0.0),
                    anneal_hits: count_hits_both_strands(
                        &destination.get_forward_string().to_ascii_uppercase(),
                        &right_overlap.2,
                    ),
                    note: "5' non-priming overlap segment derived from the right destination flank and reverse-complemented for the right insert primer".to_string(),
                },
                priming_3prime: GibsonPrimerSegment {
                    note: "3' gene-specific priming segment chosen from the opposite insert terminus to satisfy the Gibson PCR target window".to_string(),
                    ..right_priming
                },
            });
        }
    }

    apply_uniqueness_advisories(
        &plan.validation_policy.uniqueness_checks.destination_context,
        "destination overlap",
        &destination.get_forward_string().to_ascii_uppercase(),
        &[left_overlap.2.clone(), right_overlap.2.clone()],
        &mut preview,
    );
    apply_uniqueness_advisories(
        &plan.validation_policy.uniqueness_checks.participating_fragments,
        "insert priming segment",
        &insert_template_forward,
        &preview
            .primer_suggestions
            .iter()
            .map(|primer| primer.priming_3prime.sequence.clone())
            .collect::<Vec<_>>(),
        &mut preview,
    );

    preview.notes.push(format!(
        "Destination opening resolves '{}' and '{}' as the two terminal Gibson junctions.",
        resolved_opening.left_end_id, resolved_opening.right_end_id
    ));
    preview.notes.push(
        "Primer suggestions stay Gibson-specific: 5' overlap plus 3' gene-specific priming segment. Full generic PCR controls remain outside this specialist flow.".to_string(),
    );
    preview.notes.push(
        "The protocol-cartoon payload remains conceptual for the two-fragment Gibson mechanism while captions and labels are rebound from the resolved plan.".to_string(),
    );

    let representative_overlap_bp = if left_overlap.1 == right_overlap.1 {
        left_overlap.1
    } else {
        left_overlap.1.min(right_overlap.1)
    };
    preview.cartoon = GibsonCartoonPreview {
        protocol_id: ProtocolCartoonKind::GibsonTwoFragment.id().to_string(),
        representative_overlap_bp,
        bindings: build_cartoon_bindings(
            &preview,
            &left_overlap.2,
            &right_overlap.2,
            representative_overlap_bp,
        ),
    };
    preview.routine_handoff = build_routine_handoff(&preview, representative_overlap_bp);

    finalize_preview(preview)
}

fn empty_protocol_cartoon_bindings() -> ProtocolCartoonTemplateBindings {
    ProtocolCartoonTemplateBindings {
        schema: "gentle.protocol_cartoon_template_bindings.v1".to_string(),
        template_id: None,
        defaults: Default::default(),
        event_overrides: vec![],
        molecule_overrides: vec![],
        feature_overrides: vec![],
    }
}

fn finalize_preview(mut preview: GibsonAssemblyPreview) -> Result<GibsonAssemblyPreview, EngineError> {
    preview.can_execute = preview.errors.is_empty();
    Ok(preview)
}

fn validate_destination_topology(
    destination: &GibsonPlanDestination,
    dna: &crate::dna_sequence::DNAsequence,
    preview: &mut GibsonAssemblyPreview,
) {
    let requested = destination.topology_before_opening.trim();
    if requested.is_empty() {
        return;
    }
    let actual = if dna.is_circular() { "circular" } else { "linear" };
    if !requested.eq_ignore_ascii_case(actual) {
        preview.errors.push(format!(
            "destination.topology_before_opening='{}' does not match the actual topology '{}' of '{}'",
            requested, actual, destination.seq_id
        ));
    }
}

fn validate_single_insert_shape(
    plan: &GibsonAssemblyPlan,
    fragment_id: &str,
    preview: &mut GibsonAssemblyPreview,
) {
    if !plan.product.topology.trim().is_empty()
        && !matches!(
            plan.product.topology.trim().to_ascii_lowercase().as_str(),
            "linear" | "circular"
        )
    {
        preview.errors.push(format!(
            "Unsupported product.topology '{}' for Gibson specialist v1",
            plan.product.topology.trim()
        ));
    }

    if plan.assembly_order.len() != 3 {
        preview.errors.push(format!(
            "Gibson specialist v1 expects assembly_order to contain exactly 3 members (left end -> insert -> right end), but found {}",
            plan.assembly_order.len()
        ));
        return;
    }
    let left = &plan.assembly_order[0];
    let middle = &plan.assembly_order[1];
    let right = &plan.assembly_order[2];
    if !left.kind.eq_ignore_ascii_case("destination_end")
        || !middle.kind.eq_ignore_ascii_case("fragment")
        || !right.kind.eq_ignore_ascii_case("destination_end")
    {
        preview.errors.push(
            "Gibson specialist v1 expects assembly_order kinds [destination_end, fragment, destination_end]"
                .to_string(),
        );
    }
    if !middle.id.eq_ignore_ascii_case(fragment_id) {
        preview.errors.push(format!(
            "assembly_order middle fragment '{}' does not match the single insert fragment '{}'",
            middle.id, fragment_id
        ));
    }
}

fn resolve_opening(
    plan: &GibsonAssemblyPlan,
    destination: &crate::dna_sequence::DNAsequence,
    preview: &mut GibsonAssemblyPreview,
) -> Option<GibsonResolvedOpening> {
    let mode = normalized_opening_mode(&plan.destination.opening.mode, destination.is_circular());
    let left_end_id = if plan.destination.opening.left_end_id.trim().is_empty() {
        "dest_left".to_string()
    } else {
        plan.destination.opening.left_end_id.trim().to_string()
    };
    let right_end_id = if plan.destination.opening.right_end_id.trim().is_empty() {
        "dest_right".to_string()
    } else {
        plan.destination.opening.right_end_id.trim().to_string()
    };
    let label = if plan.destination.opening.label.trim().is_empty() {
        match mode.as_str() {
            "existing_termini" => "existing linear termini".to_string(),
            _ => "defined opening span".to_string(),
        }
    } else {
        plan.destination.opening.label.trim().to_string()
    };
    let len = destination.len();
    if len == 0 {
        preview
            .errors
            .push("Destination sequence must be non-empty for Gibson preview".to_string());
        return None;
    }
    match mode.as_str() {
        "existing_termini" => {
            if destination.is_circular() {
                preview.errors.push(
                    "opening.mode='existing_termini' requires a linear destination sequence"
                        .to_string(),
                );
                return None;
            }
            Some(GibsonResolvedOpening {
                mode,
                label,
                start_0based: Some(0),
                end_0based_exclusive: Some(len),
                removed_span_bp: Some(0),
                left_end_id,
                right_end_id,
                destination_is_circular: false,
            })
        }
        _ => {
            let Some(start) = plan.destination.opening.start_0based else {
                preview.errors.push(
                    "opening.mode='defined_site' requires destination.opening.start_0based"
                        .to_string(),
                );
                return None;
            };
            let Some(end) = plan.destination.opening.end_0based_exclusive else {
                preview.errors.push(
                    "opening.mode='defined_site' requires destination.opening.end_0based_exclusive"
                        .to_string(),
                );
                return None;
            };
            if start >= end || end > len {
                preview.errors.push(format!(
                    "Destination opening span {}..{} is invalid for sequence length {}",
                    start, end, len
                ));
                return None;
            }
            Some(GibsonResolvedOpening {
                mode: "defined_site".to_string(),
                label,
                start_0based: Some(start),
                end_0based_exclusive: Some(end),
                removed_span_bp: Some(end.saturating_sub(start)),
                left_end_id,
                right_end_id,
                destination_is_circular: destination.is_circular(),
            })
        }
    }
}

#[derive(Clone, Copy)]
enum TerminalSide {
    Left,
    Right,
}

fn derive_terminal_overlap(
    destination_seq: &str,
    opening: &GibsonResolvedOpening,
    side: TerminalSide,
    junction: &GibsonPlanJunction,
    targets: &GibsonPlanDesignTargets,
    preview: &mut GibsonAssemblyPreview,
) -> Option<(String, usize, String)> {
    let Some(partition) = &junction.overlap_partition else {
        preview.errors.push(format!(
            "Junction '{}' requires overlap_partition for Gibson specialist v1",
            junction.id
        ));
        return None;
    };
    let required_overlap_bp = junction.required_overlap_bp.unwrap_or(0);
    if required_overlap_bp > 0
        && partition.left_member_bp.saturating_add(partition.right_member_bp) != required_overlap_bp
    {
        preview.errors.push(format!(
            "Junction '{}' overlap_partition {}+{} does not match required_overlap_bp {}",
            junction.id,
            partition.left_member_bp,
            partition.right_member_bp,
            required_overlap_bp
        ));
        return None;
    }
    match side {
        TerminalSide::Left => {
            if partition.right_member_bp != 0
                || (required_overlap_bp > 0 && partition.left_member_bp == 0)
            {
                preview.errors.push(format!(
                    "Junction '{}' is not destination-derived on the left side; Gibson specialist v1 currently expects left terminal overlap_partition to be N+0",
                    junction.id
                ));
                return None;
            }
        }
        TerminalSide::Right => {
            if partition.left_member_bp != 0
                || (required_overlap_bp > 0 && partition.right_member_bp == 0)
            {
                preview.errors.push(format!(
                    "Junction '{}' is not destination-derived on the right side; Gibson specialist v1 currently expects right terminal overlap_partition to be 0+N",
                    junction.id
                ));
                return None;
            }
        }
    }

    let preferred_len = if required_overlap_bp > 0 {
        required_overlap_bp
    } else {
        0
    };
    let actual_len = if preferred_len > 0 {
        preferred_len
    } else {
        match choose_overlap_len_for_side(destination_seq, opening, side, targets) {
            Some(value) => value,
            None => {
                preview.errors.push(format!(
                    "Could not derive an overlap length for junction '{}' within {}..{} bp at minimum overlap Tm {:.1} C",
                    junction.id,
                    targets.overlap_bp_min,
                    targets.overlap_bp_max,
                    targets.minimum_overlap_tm_celsius
                ));
                return None;
            }
        }
    };

    if actual_len < targets.overlap_bp_min || actual_len > targets.overlap_bp_max {
        preview.warnings.push(format!(
            "Junction '{}' uses {} bp overlap outside the preferred Gibson target range {}..{} bp",
            junction.id, actual_len, targets.overlap_bp_min, targets.overlap_bp_max
        ));
    }
    let sequence = match overlap_sequence_for_side(destination_seq, opening, side, actual_len) {
        Some(seq) => seq,
        None => {
            preview.errors.push(format!(
                "Destination does not have enough flank sequence to derive {} bp overlap for junction '{}'",
                actual_len, junction.id
            ));
            return None;
        }
    };
    let tm = GentleEngine::estimate_primer_tm_c(sequence.as_bytes());
    if tm < targets.minimum_overlap_tm_celsius {
        preview.errors.push(format!(
            "Junction '{}' overlap Tm {:.1} C is below the configured minimum {:.1} C",
            junction.id, tm, targets.minimum_overlap_tm_celsius
        ));
    }
    Some((junction.id.clone(), actual_len, sequence))
}

fn choose_overlap_len_for_side(
    destination_seq: &str,
    opening: &GibsonResolvedOpening,
    side: TerminalSide,
    targets: &GibsonPlanDesignTargets,
) -> Option<usize> {
    let mut best: Option<usize> = None;
    let target_mid = (targets.overlap_bp_min + targets.overlap_bp_max) / 2;
    for len in targets.overlap_bp_min..=targets.overlap_bp_max {
        let Some(sequence) = overlap_sequence_for_side(destination_seq, opening, side, len) else {
            continue;
        };
        let tm = GentleEngine::estimate_primer_tm_c(sequence.as_bytes());
        if tm < targets.minimum_overlap_tm_celsius {
            continue;
        }
        match best {
            None => best = Some(len),
            Some(previous) => {
                let better = len.abs_diff(target_mid) < previous.abs_diff(target_mid)
                    || (len.abs_diff(target_mid) == previous.abs_diff(target_mid)
                        && len < previous);
                if better {
                    best = Some(len);
                }
            }
        }
    }
    best
}

fn overlap_sequence_for_side(
    destination_seq: &str,
    opening: &GibsonResolvedOpening,
    side: TerminalSide,
    len: usize,
) -> Option<String> {
    if len == 0 || destination_seq.is_empty() {
        return None;
    }
    match opening.mode.as_str() {
        "existing_termini" => match side {
            TerminalSide::Left => destination_seq.get(0..len).map(str::to_string),
            TerminalSide::Right => {
                let total = destination_seq.len();
                total.checked_sub(len)
                    .and_then(|start| destination_seq.get(start..total))
                    .map(str::to_string)
            }
        },
        _ => {
            let start = opening.start_0based?;
            let end = opening.end_0based_exclusive?;
            if opening.destination_is_circular {
                match side {
                    TerminalSide::Left => circular_window_suffix(destination_seq, start, len),
                    TerminalSide::Right => circular_window_prefix(destination_seq, end, len),
                }
            } else {
                match side {
                    TerminalSide::Left => start
                        .checked_sub(len)
                        .and_then(|left_start| destination_seq.get(left_start..start))
                        .map(str::to_string),
                    TerminalSide::Right => end
                        .checked_add(len)
                        .and_then(|right_end| destination_seq.get(end..right_end))
                        .map(str::to_string),
                }
            }
        }
    }
}

fn circular_window_suffix(sequence: &str, end_exclusive: usize, len: usize) -> Option<String> {
    if sequence.is_empty() || len == 0 || end_exclusive > sequence.len() {
        return None;
    }
    let total = sequence.len();
    let start = (end_exclusive + total - (len % total)) % total;
    circular_window(sequence, start, len)
}

fn circular_window_prefix(sequence: &str, start: usize, len: usize) -> Option<String> {
    if sequence.is_empty() || len == 0 || start > sequence.len() {
        return None;
    }
    circular_window(sequence, start % sequence.len(), len)
}

fn circular_window(sequence: &str, start: usize, len: usize) -> Option<String> {
    if sequence.is_empty() || len == 0 {
        return None;
    }
    let total = sequence.len();
    let mut out = String::with_capacity(len);
    for idx in 0..len {
        let pos = (start + idx) % total;
        let ch = sequence.as_bytes().get(pos).copied()?;
        out.push(ch as char);
    }
    Some(out)
}

fn choose_insert_priming_segment(
    insert_template_forward: &str,
    oriented_insert_seq: &str,
    side: TerminalSide,
    targets: &GibsonPlanDesignTargets,
    preview: &mut GibsonAssemblyPreview,
) -> Option<GibsonPrimerSegment> {
    if oriented_insert_seq.is_empty() {
        preview
            .errors
            .push("Insert sequence must be non-empty for Gibson primer derivation".to_string());
        return None;
    }
    let target_mid = (targets.priming_segment_tm_min_celsius
        + targets.priming_segment_tm_max_celsius)
        / 2.0;
    let mut best: Option<GibsonPrimerCandidate> = None;
    for len in targets.priming_segment_min_length_bp..=targets.priming_segment_max_length_bp {
        let sequence = match side {
            TerminalSide::Left => oriented_insert_seq.get(0..len).map(str::to_string),
            TerminalSide::Right => oriented_insert_seq
                .get(oriented_insert_seq.len().saturating_sub(len)..oriented_insert_seq.len())
                .map(GentleEngine::reverse_complement),
        };
        let Some(sequence) = sequence else {
            continue;
        };
        let tm = GentleEngine::estimate_primer_tm_c(sequence.as_bytes());
        let gc_fraction = GentleEngine::sequence_gc_fraction(sequence.as_bytes()).unwrap_or(0.0);
        let anneal_hits = count_hits_both_strands(insert_template_forward, &sequence);
        if tm < targets.priming_segment_tm_min_celsius
            || tm > targets.priming_segment_tm_max_celsius
            || anneal_hits == 0
            || anneal_hits > targets.max_anneal_hits
        {
            continue;
        }
        let candidate = GibsonPrimerCandidate {
            sequence,
            tm_celsius: tm,
            gc_fraction,
            anneal_hits,
        };
        let replace = match &best {
            None => true,
            Some(previous) => {
                let prev_dist = (previous.tm_celsius - target_mid).abs();
                let next_dist = (candidate.tm_celsius - target_mid).abs();
                next_dist < prev_dist
                    || (next_dist - prev_dist).abs() < f64::EPSILON
                        && candidate.anneal_hits < previous.anneal_hits
                    || (next_dist - prev_dist).abs() < f64::EPSILON
                        && candidate.anneal_hits == previous.anneal_hits
                        && candidate.sequence.len() < previous.sequence.len()
            }
        };
        if replace {
            best = Some(candidate);
        }
    }

    let Some(best) = best else {
        preview.errors.push(format!(
            "Could not derive a Gibson priming segment within {}..{} bp and {:.1}..{:.1} C for the {} insert end",
            targets.priming_segment_min_length_bp,
            targets.priming_segment_max_length_bp,
            targets.priming_segment_tm_min_celsius,
            targets.priming_segment_tm_max_celsius,
            match side {
                TerminalSide::Left => "left",
                TerminalSide::Right => "right",
            }
        ));
        return None;
    };
    Some(GibsonPrimerSegment {
        length_bp: best.sequence.len(),
        sequence: best.sequence,
        tm_celsius: best.tm_celsius,
        gc_fraction: best.gc_fraction,
        anneal_hits: best.anneal_hits,
        note: String::new(),
    })
}

fn count_hits_both_strands(template_forward: &str, query: &str) -> usize {
    if query.is_empty() || template_forward.is_empty() {
        return 0;
    }
    let forward = GentleEngine::find_all_subsequences(template_forward.as_bytes(), query.as_bytes());
    let reverse_template = GentleEngine::reverse_complement(template_forward);
    let reverse = GentleEngine::find_all_subsequences(reverse_template.as_bytes(), query.as_bytes());
    forward.len().saturating_add(reverse.len())
}

fn apply_uniqueness_advisories(
    severity: &str,
    label: &str,
    context_sequence: &str,
    sequences: &[String],
    preview: &mut GibsonAssemblyPreview,
) {
    let level = severity.trim().to_ascii_lowercase();
    if level.is_empty() || level == "off" {
        return;
    }
    for sequence in sequences {
        if sequence.trim().is_empty() {
            continue;
        }
        let hits = count_hits_both_strands(context_sequence, sequence);
        if hits <= 1 {
            continue;
        }
        let message = format!(
            "{} '{}' is not unique in its checked context ({} exact hits across both strands)",
            label, sequence, hits
        );
        if level == "error" {
            preview.errors.push(message);
        } else {
            preview.warnings.push(message);
        }
    }
}

fn build_cartoon_bindings(
    preview: &GibsonAssemblyPreview,
    left_overlap_seq: &str,
    right_overlap_seq: &str,
    representative_overlap_bp: usize,
) -> ProtocolCartoonTemplateBindings {
    let destination_label = if preview.destination.seq_id.trim().is_empty() {
        "Destination".to_string()
    } else {
        preview.destination.seq_id.clone()
    };
    let insert_label = if preview.insert.seq_id.trim().is_empty() {
        "Insert".to_string()
    } else {
        preview.insert.seq_id.clone()
    };
    let left_bp = preview
        .resolved_junctions
        .first()
        .map(|row| row.overlap_bp)
        .unwrap_or(representative_overlap_bp);
    let right_bp = preview
        .resolved_junctions
        .get(1)
        .map(|row| row.overlap_bp)
        .unwrap_or(representative_overlap_bp);
    let joined_label = format!("{} + {}", destination_label, insert_label);
    let overlap_note = if left_bp == right_bp {
        format!("{} bp terminal overlap", left_bp)
    } else {
        format!("{} bp / {} bp terminal overlaps", left_bp, right_bp)
    };
    ProtocolCartoonTemplateBindings {
        schema: "gentle.protocol_cartoon_template_bindings.v1".to_string(),
        template_id: Some(ProtocolCartoonKind::GibsonTwoFragment.id().to_string()),
        defaults: Default::default(),
        event_overrides: vec![
            ProtocolCartoonTemplateEventBinding {
                event_id: "context".to_string(),
                title: Some("Context".to_string()),
                caption: Some(format!(
                    "The opened destination and insert are configured with {} and {} bp terminal homology regions (left='{}', right='{}').",
                    left_bp, right_bp, left_overlap_seq, right_overlap_seq
                )),
                action: None,
            },
            ProtocolCartoonTemplateEventBinding {
                event_id: "chew_back".to_string(),
                title: Some("Chew-back".to_string()),
                caption: Some(format!(
                    "A 5' exonuclease exposes single-stranded ends so the {} can anneal during Gibson assembly.",
                    overlap_note
                )),
                action: None,
            },
            ProtocolCartoonTemplateEventBinding {
                event_id: "anneal".to_string(),
                title: Some("Anneal".to_string()),
                caption: Some(format!(
                    "Complementary terminal homology anneals while destination- and insert-derived unique tails remain to be filled." 
                )),
                action: None,
            },
            ProtocolCartoonTemplateEventBinding {
                event_id: "extend".to_string(),
                title: Some("Extend".to_string()),
                caption: Some(
                    "DNA polymerase fills the remaining single-stranded tails but leaves one nick on each strand.".to_string(),
                ),
                action: None,
            },
            ProtocolCartoonTemplateEventBinding {
                event_id: "seal".to_string(),
                title: Some("Seal".to_string()),
                caption: Some(
                    "DNA ligase seals the nicks, leaving one continuous assembled duplex with preserved origin colors.".to_string(),
                ),
                action: None,
            },
        ],
        molecule_overrides: vec![
            ProtocolCartoonTemplateMoleculeBinding {
                event_id: "context".to_string(),
                molecule_id: "fragment_a_context".to_string(),
                label: Some(destination_label.clone()),
                topology: None,
                left_end: Some(DnaEndStyle::Continuation),
                right_end: Some(DnaEndStyle::Blunt),
            },
            ProtocolCartoonTemplateMoleculeBinding {
                event_id: "context".to_string(),
                molecule_id: "fragment_b_context".to_string(),
                label: Some(insert_label.clone()),
                topology: None,
                left_end: Some(DnaEndStyle::Blunt),
                right_end: Some(DnaEndStyle::Continuation),
            },
            ProtocolCartoonTemplateMoleculeBinding {
                event_id: "chew_back".to_string(),
                molecule_id: "fragment_a_chewed".to_string(),
                label: Some(format!("{} (chewed)", destination_label)),
                topology: None,
                left_end: Some(DnaEndStyle::Continuation),
                right_end: Some(DnaEndStyle::Sticky {
                    polarity: OverhangPolarity::ThreePrime,
                    nt: representative_overlap_bp.max(1),
                }),
            },
            ProtocolCartoonTemplateMoleculeBinding {
                event_id: "chew_back".to_string(),
                molecule_id: "fragment_b_chewed".to_string(),
                label: Some(format!("{} (chewed)", insert_label)),
                topology: None,
                left_end: Some(DnaEndStyle::Sticky {
                    polarity: OverhangPolarity::ThreePrime,
                    nt: representative_overlap_bp.max(1),
                }),
                right_end: Some(DnaEndStyle::Continuation),
            },
            ProtocolCartoonTemplateMoleculeBinding {
                event_id: "anneal".to_string(),
                molecule_id: "annealed_intermediate".to_string(),
                label: Some(format!("{} annealed", joined_label)),
                topology: None,
                left_end: Some(DnaEndStyle::Continuation),
                right_end: Some(DnaEndStyle::Continuation),
            },
            ProtocolCartoonTemplateMoleculeBinding {
                event_id: "extend".to_string(),
                molecule_id: "extended_intermediate".to_string(),
                label: Some(format!("{} extended duplex", joined_label)),
                topology: None,
                left_end: Some(DnaEndStyle::Continuation),
                right_end: Some(DnaEndStyle::Continuation),
            },
            ProtocolCartoonTemplateMoleculeBinding {
                event_id: "seal".to_string(),
                molecule_id: "sealed_product".to_string(),
                label: Some(format!("{} sealed product", joined_label)),
                topology: None,
                left_end: Some(DnaEndStyle::Continuation),
                right_end: Some(DnaEndStyle::Continuation),
            },
        ],
        feature_overrides: vec![
            ProtocolCartoonTemplateFeatureBinding {
                event_id: "context".to_string(),
                molecule_id: "fragment_a_context".to_string(),
                feature_id: "a_overlap".to_string(),
                label: Some(format!("Left {} bp overlap", left_bp)),
                length_bp: Some(representative_overlap_bp.max(1)),
                top_length_bp: Some(representative_overlap_bp.max(1)),
                bottom_length_bp: Some(representative_overlap_bp.max(1)),
                color_hex: None,
                bottom_color_hex: None,
                top_nick_after: None,
                bottom_nick_after: None,
            },
            ProtocolCartoonTemplateFeatureBinding {
                event_id: "context".to_string(),
                molecule_id: "fragment_b_context".to_string(),
                feature_id: "b_overlap".to_string(),
                label: Some(format!("Right {} bp overlap", right_bp)),
                length_bp: Some(representative_overlap_bp.max(1)),
                top_length_bp: Some(representative_overlap_bp.max(1)),
                bottom_length_bp: Some(representative_overlap_bp.max(1)),
                color_hex: None,
                bottom_color_hex: None,
                top_nick_after: None,
                bottom_nick_after: None,
            },
            ProtocolCartoonTemplateFeatureBinding {
                event_id: "anneal".to_string(),
                molecule_id: "annealed_intermediate".to_string(),
                feature_id: "hybrid_overlap".to_string(),
                label: Some(overlap_note.clone()),
                length_bp: Some(representative_overlap_bp.max(1)),
                top_length_bp: Some(representative_overlap_bp.max(1)),
                bottom_length_bp: Some(representative_overlap_bp.max(1)),
                color_hex: None,
                bottom_color_hex: None,
                top_nick_after: None,
                bottom_nick_after: None,
            },
        ],
    }
}

fn build_routine_handoff(
    preview: &GibsonAssemblyPreview,
    representative_overlap_bp: usize,
) -> GibsonRoutineHandoffPreview {
    if preview.destination.opening_mode.eq_ignore_ascii_case("existing_termini")
        && preview.destination.actual_topology.eq_ignore_ascii_case("linear")
    {
        let mut bindings = BTreeMap::new();
        bindings.insert(
            "left_seq_id".to_string(),
            preview.destination.seq_id.clone(),
        );
        bindings.insert("right_seq_id".to_string(), preview.insert.seq_id.clone());
        bindings.insert(
            "overlap_bp".to_string(),
            representative_overlap_bp.to_string(),
        );
        if !preview.destination.seq_id.trim().is_empty() {
            bindings.insert(
                "assembly_prefix".to_string(),
                format!("{}_gibson", preview.destination.seq_id),
            );
        }
        if !preview.destination.seq_id.trim().is_empty() && !preview.insert.seq_id.trim().is_empty() {
            bindings.insert(
                "output_id".to_string(),
                format!("{}_with_{}", preview.destination.seq_id, preview.insert.seq_id),
            );
        }
        GibsonRoutineHandoffPreview {
            supported: true,
            routine_id: "gibson.two_fragment_overlap_preview".to_string(),
            reason: None,
            bindings,
        }
    } else {
        GibsonRoutineHandoffPreview {
            supported: false,
            routine_id: "gibson.two_fragment_overlap_preview".to_string(),
            reason: Some(
                "Routine Assistant handoff is currently automatic only when the destination is already linear and uses existing termini. Defined opening spans remain preview-only in this v1 specialist.".to_string(),
            ),
            bindings: BTreeMap::new(),
        }
    }
}

fn normalized_opening_mode(raw: &str, is_circular: bool) -> String {
    match raw.trim().to_ascii_lowercase().as_str() {
        "existing_termini" | "existing-termini" => "existing_termini".to_string(),
        "defined_span" | "defined_site" | "selection" | "window" => "defined_site".to_string(),
        "" if !is_circular => "existing_termini".to_string(),
        _ => "defined_site".to_string(),
    }
}

fn normalized_orientation(raw: &str) -> String {
    match raw.trim().to_ascii_lowercase().as_str() {
        "reverse" | "rev" | "reverse_complement" => "reverse".to_string(),
        _ => "forward".to_string(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{dna_sequence::DNAsequence, engine::ProjectState, enzymes::active_restriction_enzymes};
    use gb_io::seq::{Feature, FeatureKind, Location};

    fn test_engine_with_sequences() -> GentleEngine {
        let mut engine = GentleEngine::new();
        let mut destination = DNAsequence::from_sequence(
            "AAACCCGGGTTTAAACCCGGGTTTAAACCCGGGTTTAAACCCGGGTTT",
        )
        .expect("destination sequence");
        destination.set_name("destination_vector".to_string());
        destination.set_circular(true);
        let mut insert = DNAsequence::from_sequence(
            "ATGCGTACGTTAGCGTACGATCGTACGTAGCTAGCTAGCATCGATCGA",
        )
        .expect("insert sequence");
        insert.set_name("insert_x_amplicon".to_string());
        insert.set_circular(false);
        engine.state_mut().sequences.insert("destination_vector".to_string(), destination);
        engine.state_mut().sequences.insert("insert_x_amplicon".to_string(), insert);
        engine.state_mut().display = ProjectState::default().display;
        engine
    }

    fn restriction_ready_dna(sequence: &str) -> DNAsequence {
        let mut dna = DNAsequence::from_sequence(sequence).expect("sequence");
        *dna.restriction_enzymes_mut() = active_restriction_enzymes();
        dna.set_max_restriction_enzyme_sites(Some(2));
        dna.update_computed_features();
        dna
    }

    fn single_insert_plan() -> GibsonAssemblyPlan {
        serde_json::from_str(
            r#"{
  "schema": "gentle.gibson_assembly_plan.v1",
  "id": "preview_test",
  "title": "Preview test",
  "summary": "single insert",
  "destination": {
    "seq_id": "destination_vector",
    "topology_before_opening": "circular",
    "opening": {
      "mode": "defined_site",
      "label": "selected window",
      "start_0based": 12,
      "end_0based_exclusive": 18,
      "left_end_id": "dest_left",
      "right_end_id": "dest_right",
      "uniqueness_requirement": "must_be_unambiguous"
    }
  },
  "product": {"topology": "circular", "output_id_hint": "out"},
  "fragments": [
    {
      "id": "insert_x",
      "seq_id": "insert_x_amplicon",
      "role": "insert",
      "orientation": "forward",
      "left_end_strategy": {"mode": "primer_added_overlap", "target_junction_id": "junction_left"},
      "right_end_strategy": {"mode": "primer_added_overlap", "target_junction_id": "junction_right"}
    }
  ],
  "assembly_order": [
    {"kind": "destination_end", "id": "dest_left"},
    {"kind": "fragment", "id": "insert_x"},
    {"kind": "destination_end", "id": "dest_right"}
  ],
  "junctions": [
    {
      "id": "junction_left",
      "left_member": {"kind": "destination_end", "id": "dest_left"},
      "right_member": {"kind": "fragment", "id": "insert_x"},
      "required_overlap_bp": 20,
      "overlap_partition": {"left_member_bp": 20, "right_member_bp": 0},
      "overlap_source": "derive_from_destination_left_flank",
      "distinct_from": ["junction_right"]
    },
    {
      "id": "junction_right",
      "left_member": {"kind": "fragment", "id": "insert_x"},
      "right_member": {"kind": "destination_end", "id": "dest_right"},
      "required_overlap_bp": 20,
      "overlap_partition": {"left_member_bp": 0, "right_member_bp": 20},
      "overlap_source": "derive_from_destination_right_flank",
      "distinct_from": ["junction_left"]
    }
  ],
  "validation_policy": {
    "require_unambiguous_destination_opening": true,
    "require_distinct_terminal_junctions": true,
    "adjacency_overlap_mismatch": "error",
    "design_targets": {
      "overlap_bp_min": 18,
      "overlap_bp_max": 24,
      "minimum_overlap_tm_celsius": 48.0,
      "priming_segment_tm_min_celsius": 48.0,
      "priming_segment_tm_max_celsius": 70.0,
      "priming_segment_min_length_bp": 18,
      "priming_segment_max_length_bp": 30,
      "max_anneal_hits": 4
    },
    "uniqueness_checks": {
      "destination_context": "warn",
      "participating_fragments": "warn",
      "reference_contexts": []
    }
  }
}"#,
        )
        .expect("plan json")
    }

    #[test]
    fn preview_single_insert_plan_returns_junctions_and_primers() {
        let engine = test_engine_with_sequences();
        let preview = preview_gibson_assembly_plan(&engine, &single_insert_plan())
            .expect("preview output");
        assert!(preview.can_execute);
        assert_eq!(preview.resolved_junctions.len(), 2);
        assert_eq!(preview.primer_suggestions.len(), 2);
        assert_eq!(preview.cartoon.protocol_id, "gibson.two_fragment");
    }

    #[test]
    fn suggest_gibson_destination_openings_prefers_mcs_blunt_unique_cutter() {
        let mut dna = restriction_ready_dna("TTTGGATCCAAACCCGGGTTTGAATTCTTT");
        dna.set_circular(true);
        dna.features_mut().push(Feature {
            kind: FeatureKind::from("misc_feature"),
            location: Location::simple_range(3, 21),
            qualifiers: vec![
                (
                    "label".into(),
                    Some("Multiple Cloning Site (MCS)".to_string()),
                ),
                (
                    "note".into(),
                    Some("contains BamHI, SmaI and EcoR I".to_string()),
                ),
            ],
        });
        let mut engine = GentleEngine::new();
        engine.state_mut().sequences.insert("vector".to_string(), dna);

        let suggestions =
            suggest_gibson_destination_openings(&engine, "vector").expect("suggestions");

        assert!(!suggestions.is_empty());
        assert_eq!(suggestions[0].kind, "mcs_feature");
        let site_order = suggestions
            .iter()
            .filter_map(|row| row.enzyme_name.as_deref())
            .collect::<Vec<_>>();
        assert_eq!(site_order.first().copied(), Some("SmaI"));
        let smai = suggestions
            .iter()
            .find(|row| row.enzyme_name.as_deref() == Some("SmaI"))
            .expect("SmaI suggestion");
        assert_eq!(smai.end_geometry, "blunt");
        assert!(smai.in_mcs_context);
    }

    #[test]
    fn preview_requires_valid_opening_span() {
        let engine = test_engine_with_sequences();
        let mut plan = single_insert_plan();
        plan.destination.opening.start_0based = Some(25);
        plan.destination.opening.end_0based_exclusive = Some(10);
        let preview = preview_gibson_assembly_plan(&engine, &plan).expect("preview output");
        assert!(!preview.can_execute);
        assert!(preview
            .errors
            .iter()
            .any(|row| row.contains("Destination opening span")));
    }

    #[test]
    fn preview_flags_non_distinct_terminal_overlaps() {
        let engine = test_engine_with_sequences();
        let mut plan = single_insert_plan();
        plan.destination.opening.start_0based = Some(3);
        plan.destination.opening.end_0based_exclusive = Some(9);
        let preview = preview_gibson_assembly_plan(&engine, &plan).expect("preview output");
        assert!(!preview.can_execute);
        assert!(preview
            .errors
            .iter()
            .any(|row| row.contains("Terminal overlap regions")));
    }

    #[test]
    fn preview_reports_impossible_priming_targets() {
        let engine = test_engine_with_sequences();
        let mut plan = single_insert_plan();
        plan.validation_policy.design_targets.priming_segment_tm_min_celsius = 200.0;
        plan.validation_policy.design_targets.priming_segment_tm_max_celsius = 210.0;
        let preview = preview_gibson_assembly_plan(&engine, &plan).expect("preview output");
        assert!(!preview.can_execute);
        assert!(preview
            .errors
            .iter()
            .any(|row| row.contains("Could not derive a Gibson priming segment")));
    }

    #[test]
    fn preview_existing_linear_termini_enables_routine_handoff() {
        let mut engine = GentleEngine::new();
        let destination = DNAsequence::from_sequence(
            "ATGCGTACGATTCGATGCACTGATCGTACGTAGCTAACGTAGGCTAAC",
        )
        .expect("destination sequence");
        let insert = DNAsequence::from_sequence(
            "ATGCGTACGTTAGCGTACGATCGTACGTAGCTAGCTAGCATCGATCGA",
        )
        .expect("insert sequence");
        engine
            .state_mut()
            .sequences
            .insert("destination_vector".to_string(), destination);
        engine
            .state_mut()
            .sequences
            .insert("insert_x_amplicon".to_string(), insert);
        let mut plan = single_insert_plan();
        plan.destination.topology_before_opening = "linear".to_string();
        plan.destination.opening.mode = "existing_termini".to_string();
        plan.destination.opening.start_0based = None;
        plan.destination.opening.end_0based_exclusive = None;
        let preview = preview_gibson_assembly_plan(&engine, &plan).expect("preview output");
        assert!(preview.can_execute);
        assert!(preview.routine_handoff.supported);
        assert_eq!(
            preview.routine_handoff.routine_id,
            "gibson.two_fragment_overlap_preview"
        );
    }
}
