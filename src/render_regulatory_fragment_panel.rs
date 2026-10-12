//! Deterministic SVG presentation for regulatory-fragment contrast plans.
//!
//! The renderer is deliberately passive: it displays the exact engine-owned
//! plan, including evidence states and non-claims, without recomputing biology.

use crate::engine::{
    RegulatoryFragmentEvidenceDimensionKind, RegulatoryFragmentEvidenceState,
    RegulatoryFragmentOrientation, RegulatoryFragmentPanelPlan, RegulatoryFragmentRole,
};
use std::collections::BTreeMap;
use std::fmt::Write;

const WIDTH: f32 = 1200.0;
const LEFT: f32 = 190.0;
const RIGHT: f32 = 32.0;
const TRACK_WIDTH: f32 = WIDTH - LEFT - RIGHT;
const ROW_HEIGHT: f32 = 28.0;

fn escape(raw: &str) -> String {
    raw.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
        .replace('\'', "&apos;")
}

fn role_token(role: RegulatoryFragmentRole) -> &'static str {
    match role {
        RegulatoryFragmentRole::Candidate => "candidate",
        RegulatoryFragmentRole::Partner => "partner",
        RegulatoryFragmentRole::MinimalPromoter => "minimal promoter",
        RegulatoryFragmentRole::ReferenceControl => "reference control",
    }
}

fn role_color(role: RegulatoryFragmentRole) -> &'static str {
    match role {
        RegulatoryFragmentRole::Candidate => "#15803d",
        RegulatoryFragmentRole::Partner => "#ca8a04",
        RegulatoryFragmentRole::MinimalPromoter => "#2563eb",
        RegulatoryFragmentRole::ReferenceControl => "#dc2626",
    }
}

fn evidence_kind_token(kind: RegulatoryFragmentEvidenceDimensionKind) -> &'static str {
    match kind {
        RegulatoryFragmentEvidenceDimensionKind::ReferenceGenomicUniqueness => {
            "Reference-context matches"
        }
        RegulatoryFragmentEvidenceDimensionKind::PanelSequenceSimilarity => {
            "Panel / vector similarity"
        }
        RegulatoryFragmentEvidenceDimensionKind::RepeatsAndLowComplexity => {
            "Repeats / low complexity"
        }
        RegulatoryFragmentEvidenceDimensionKind::PairSpecificJunctionUniqueness => {
            "Junction exact matches"
        }
        RegulatoryFragmentEvidenceDimensionKind::RestrictionAndCloningRisk => {
            "Restriction / cloning"
        }
        RegulatoryFragmentEvidenceDimensionKind::EnsemblRegulatoryOverlap => "Ensembl regulation",
        RegulatoryFragmentEvidenceDimensionKind::TfbsModelScoreContext => "TFBS model context",
        RegulatoryFragmentEvidenceDimensionKind::CutrunAndChromatinContext => "CUT&RUN / chromatin",
    }
}

fn evidence_state_token(state: RegulatoryFragmentEvidenceState) -> &'static str {
    match state {
        RegulatoryFragmentEvidenceState::NotEvaluated => "not evaluated",
        RegulatoryFragmentEvidenceState::Evaluated => "evaluated",
        RegulatoryFragmentEvidenceState::Unavailable => "unavailable",
        RegulatoryFragmentEvidenceState::Stale => "stale",
    }
}

fn evidence_state_color(state: RegulatoryFragmentEvidenceState) -> &'static str {
    match state {
        RegulatoryFragmentEvidenceState::Evaluated => "#15803d",
        RegulatoryFragmentEvidenceState::NotEvaluated => "#64748b",
        RegulatoryFragmentEvidenceState::Unavailable => "#b45309",
        RegulatoryFragmentEvidenceState::Stale => "#b91c1c",
    }
}

fn section_title(svg: &mut String, y: f32, title: &str) {
    let _ = write!(
        svg,
        "<text x=\"24\" y=\"{y:.1}\" font-family=\"sans-serif\" font-size=\"16\" font-weight=\"700\" fill=\"#111827\">{}</text>",
        escape(title)
    );
}

fn text(svg: &mut String, x: f32, y: f32, size: f32, fill: &str, value: &str) {
    let _ = write!(
        svg,
        "<text x=\"{x:.1}\" y=\"{y:.1}\" font-family=\"sans-serif\" font-size=\"{size:.1}\" fill=\"{fill}\">{}</text>",
        escape(value)
    );
}

fn clipped(raw: &str, max_chars: usize) -> String {
    if raw.chars().count() <= max_chars {
        return raw.to_string();
    }
    let mut out = raw
        .chars()
        .take(max_chars.saturating_sub(3))
        .collect::<String>();
    out.push_str("...");
    out
}

fn row_label(svg: &mut String, baseline: f32, label: &str, full_label: &str) {
    let top = baseline - 11.0;
    let width = LEFT - 36.0;
    let _ = write!(
        svg,
        "<svg x=\"24\" y=\"{top:.1}\" width=\"{width:.1}\" height=\"16\" overflow=\"hidden\"><title>{}</title>",
        escape(full_label)
    );
    text(svg, 0.0, 11.0, 11.0, "#334155", &clipped(label, 25));
    svg.push_str("</svg>");
}

/// Render one plan as a compact, deterministic inspection figure.
pub fn render_regulatory_fragment_panel_svg(plan: &RegulatoryFragmentPanelPlan) -> String {
    let genomic_rows = plan.request.fragments.len().max(1);
    let construct_rows = plan.members.len().max(1);
    let contrast_rows = plan.contrasts.len().max(1);
    let evidence_rows = plan.evidence_dimensions.len().max(1);
    let finding_rows = (plan.blockers.len() + plan.omitted_variants.len()).clamp(1, 8);
    let nonclaim_rows = plan.nonclaims.len().clamp(1, 4);
    let height = 170.0
        + ROW_HEIGHT
            * (genomic_rows
                + construct_rows
                + contrast_rows
                + evidence_rows
                + finding_rows
                + nonclaim_rows) as f32
        + 180.0;
    let mut svg = String::new();
    let _ = write!(
        svg,
        "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{WIDTH:.0}\" height=\"{height:.0}\" viewBox=\"0 0 {WIDTH:.0} {height:.0}\" data-gentle-schema=\"gentle.regulatory_fragment_panel_plan.v1\">"
    );
    svg.push_str("<rect width=\"100%\" height=\"100%\" fill=\"#ffffff\"/>");
    text(
        &mut svg,
        24.0,
        34.0,
        22.0,
        "#111827",
        &format!("Regulatory-fragment panel: {}", plan.plan_id),
    );
    text(
        &mut svg,
        24.0,
        58.0,
        12.0,
        "#475569",
        &format!(
            "Planning label: {:?} | {} constructs | {} contrasts | proposal {}",
            plan.planning_label,
            plan.members.len(),
            plan.contrasts.len(),
            clipped(&plan.proposal_digest, 28)
        ),
    );

    let mut y = 96.0;
    section_title(&mut svg, y, "Genome-anchored source fragments");
    y += 18.0;
    let mut contig_bounds = BTreeMap::<String, (u64, u64)>::new();
    for fragment in &plan.request.fragments {
        let interval = &fragment.region.interval;
        contig_bounds
            .entry(interval.reference.contig_name.clone())
            .and_modify(|bounds| {
                bounds.0 = bounds.0.min(interval.start_0based);
                bounds.1 = bounds.1.max(interval.end_0based_exclusive);
            })
            .or_insert((interval.start_0based, interval.end_0based_exclusive));
    }
    for fragment in &plan.request.fragments {
        let interval = &fragment.region.interval;
        let (min, max) = contig_bounds
            .get(&interval.reference.contig_name)
            .copied()
            .unwrap_or((interval.start_0based, interval.end_0based_exclusive));
        let span = max.saturating_sub(min).max(1) as f32;
        let x0 = LEFT + (interval.start_0based.saturating_sub(min) as f32 / span) * TRACK_WIDTH;
        let x1 =
            LEFT + (interval.end_0based_exclusive.saturating_sub(min) as f32 / span) * TRACK_WIDTH;
        let label = format!("{} ({})", fragment.fragment_id, role_token(fragment.role));
        row_label(&mut svg, y + 16.0, &label, &label);
        let rect_y = y + 6.0;
        let rect_width = (x1 - x0).max(3.0);
        let _ = write!(
            svg,
            "<line x1=\"{LEFT:.1}\" y1=\"{:.1}\" x2=\"{:.1}\" y2=\"{:.1}\" stroke=\"#cbd5e1\"/><rect data-gentle-role=\"source-fragment\" x=\"{x0:.1}\" y=\"{rect_y:.1}\" width=\"{rect_width:.1}\" height=\"12\" fill=\"{}\" rx=\"2\"><title>{}</title></rect>",
            y + 12.0,
            LEFT + TRACK_WIDTH,
            y + 12.0,
            role_color(fragment.role),
            escape(&format!(
                "{}:{}..{} {:?}",
                interval.reference.contig_name,
                interval.start_0based,
                interval.end_0based_exclusive,
                interval.strand
            ))
        );
        text(
            &mut svg,
            LEFT,
            y + 27.0,
            9.0,
            "#64748b",
            &format!("{}:{}..{}", interval.reference.contig_name, min, max),
        );
        y += ROW_HEIGHT;
    }
    if plan.request.fragments.is_empty() {
        text(
            &mut svg,
            24.0,
            y + 14.0,
            11.0,
            "#64748b",
            "No source fragments",
        );
        y += ROW_HEIGHT;
    }

    y += 20.0;
    section_title(&mut svg, y, "Selected construct geometry");
    y += 18.0;
    let max_insert = plan
        .members
        .iter()
        .map(|member| member.insert_length_bp)
        .max()
        .unwrap_or(1)
        .max(1) as f32;
    let member_aliases: BTreeMap<_, _> = plan
        .members
        .iter()
        .enumerate()
        .map(|(index, member)| (member.member_id.as_str(), format!("C{}", index + 1)))
        .collect();
    let member_prefix = format!("{}_", plan.plan_id);
    for member in &plan.members {
        let alias = &member_aliases[member.member_id.as_str()];
        let suffix = member
            .member_id
            .strip_prefix(&member_prefix)
            .unwrap_or(&member.member_id);
        row_label(
            &mut svg,
            y + 17.0,
            &format!("{alias}: {suffix}"),
            &member.member_id,
        );
        let _ = write!(
            svg,
            "<line x1=\"{LEFT:.1}\" y1=\"{:.1}\" x2=\"{:.1}\" y2=\"{:.1}\" stroke=\"#e2e8f0\"/>",
            y + 13.0,
            LEFT + TRACK_WIDTH,
            y + 13.0
        );
        for instance in &member.instances {
            let x = LEFT + instance.assembled_start_0based as f32 / max_insert * TRACK_WIDTH;
            let length_bp = instance
                .assembled_end_0based_exclusive
                .saturating_sub(instance.assembled_start_0based);
            let width = (length_bp as f32 / max_insert * TRACK_WIDTH).max(3.0);
            let arrow = match instance.orientation {
                RegulatoryFragmentOrientation::Forward => ">",
                RegulatoryFragmentOrientation::ReverseComplement => "<",
            };
            let _ = write!(
                svg,
                "<rect data-gentle-role=\"construct-fragment\" x=\"{x:.1}\" y=\"{:.1}\" width=\"{width:.1}\" height=\"16\" fill=\"{}\" rx=\"2\"><title>{}</title></rect>",
                y + 5.0,
                role_color(instance.role),
                escape(&format!(
                    "{} {} {} bp",
                    instance.fragment_id, arrow, length_bp
                ))
            );
            let label_chars = ((width - 6.0) / 9.0).floor() as usize;
            if label_chars >= 5 {
                text(
                    &mut svg,
                    x + 3.0,
                    y + 17.0,
                    9.0,
                    "#ffffff",
                    &format!(
                        "{} {arrow}",
                        clipped(&instance.fragment_id, label_chars.saturating_sub(2).min(14))
                    ),
                );
            }
        }
        text(
            &mut svg,
            LEFT + TRACK_WIDTH - 74.0,
            y + 27.0,
            9.0,
            "#64748b",
            &format!("{} bp", member.insert_length_bp),
        );
        y += ROW_HEIGHT;
    }
    if plan.members.is_empty() {
        text(
            &mut svg,
            24.0,
            y + 14.0,
            11.0,
            "#64748b",
            "No selected constructs",
        );
        y += ROW_HEIGHT;
    }

    y += 20.0;
    section_title(&mut svg, y, "Contrast matrix");
    y += 20.0;
    for contrast in &plan.contrasts {
        let alias = |id: &str| {
            member_aliases
                .get(id)
                .map(String::as_str)
                .unwrap_or("unlisted")
        };
        let _ = write!(
            svg,
            "<g data-gentle-role=\"contrast\"><title>{}</title>",
            escape(&format!(
                "{} vs {}: {}",
                contrast.left_member_id, contrast.right_member_id, contrast.interpretation
            ))
        );
        text(
            &mut svg,
            24.0,
            y + 14.0,
            11.0,
            "#334155",
            &format!(
                "{:?}: {} vs {}",
                contrast.question,
                alias(&contrast.left_member_id),
                alias(&contrast.right_member_id)
            ),
        );
        text(
            &mut svg,
            610.0,
            y + 14.0,
            10.0,
            "#64748b",
            &clipped(&contrast.interpretation, 86),
        );
        svg.push_str("</g>");
        y += ROW_HEIGHT;
    }
    if plan.contrasts.is_empty() {
        text(
            &mut svg,
            24.0,
            y + 14.0,
            11.0,
            "#64748b",
            "No induced contrasts",
        );
        y += ROW_HEIGHT;
    }

    y += 20.0;
    section_title(&mut svg, y, "Independent evidence lanes");
    y += 20.0;
    for dimension in &plan.evidence_dimensions {
        let color = evidence_state_color(dimension.state);
        let _ = write!(
            svg,
            "<circle data-gentle-role=\"evidence-state\" cx=\"34\" cy=\"{:.1}\" r=\"5\" fill=\"{color}\"/>",
            y + 10.0
        );
        text(
            &mut svg,
            48.0,
            y + 14.0,
            11.0,
            "#334155",
            evidence_kind_token(dimension.kind),
        );
        text(
            &mut svg,
            310.0,
            y + 14.0,
            10.0,
            color,
            evidence_state_token(dimension.state),
        );
        text(
            &mut svg,
            420.0,
            y + 14.0,
            10.0,
            "#64748b",
            &format!(
                "{} observation(s), {} warning(s), {} blocker(s){} | {}",
                dimension.observations.len(),
                dimension.warnings.len(),
                dimension.blockers.len(),
                if dimension.truncated {
                    " (truncated)"
                } else {
                    ""
                },
                clipped(&dimension.detail, 102)
            ),
        );
        y += ROW_HEIGHT;
    }

    y += 20.0;
    section_title(&mut svg, y, "Blockers and omitted variants");
    y += 20.0;
    let mut findings = plan
        .blockers
        .iter()
        .map(|finding| format!("BLOCKER {}: {}", finding.code, finding.detail))
        .chain(plan.omitted_variants.iter().map(|omitted| {
            format!(
                "OMITTED {} ({:?}): {}",
                omitted.member_id, omitted.reason, omitted.detail
            )
        }))
        .take(8)
        .collect::<Vec<_>>();
    if findings.is_empty() {
        findings.push("No planner blockers or generated-variant omissions.".to_string());
    }
    for finding in findings {
        text(
            &mut svg,
            24.0,
            y + 14.0,
            10.0,
            "#7c2d12",
            &clipped(&finding, 175),
        );
        y += ROW_HEIGHT;
    }

    y += 20.0;
    section_title(&mut svg, y, "Scientific non-claims");
    y += 20.0;
    for nonclaim in plan.nonclaims.iter().take(4) {
        text(
            &mut svg,
            24.0,
            y + 14.0,
            10.0,
            "#475569",
            &format!("- {}", clipped(nonclaim, 178)),
        );
        y += ROW_HEIGHT;
    }
    if plan.nonclaims.is_empty() {
        text(&mut svg, 24.0, y + 14.0, 10.0, "#475569", "- None recorded");
    }
    text(
        &mut svg,
        WIDTH - 220.0,
        height - 16.0,
        9.0,
        "#94a3b8",
        "GENtle regulatory-fragment plan",
    );
    svg.push_str("</svg>");
    svg
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::engine::{
        RegulatoryFragmentBinding, RegulatoryFragmentContrast, RegulatoryFragmentEvidenceDimension,
        RegulatoryFragmentPanelMember, RegulatoryFragmentPanelRequest, RegulatoryFragmentQuestion,
        RegulatoryFragmentResolvedInstance,
    };
    use gentle_protocol as gp;

    #[test]
    fn regulatory_fragment_panel_svg_exposes_geometry_contrasts_and_evidence_states() {
        let mut candidate = RegulatoryFragmentBinding {
            fragment_id: "candidate_a".to_string(),
            role: RegulatoryFragmentRole::Candidate,
            ..RegulatoryFragmentBinding::default()
        };
        candidate.region.interval.reference.contig_name = "chr1".to_string();
        candidate.region.interval.start_0based = 100;
        candidate.region.interval.end_0based_exclusive = 220;
        candidate.region.interval.strand = gp::GenomicRegionStrand::Plus;
        let plan = RegulatoryFragmentPanelPlan {
            plan_id: "tp73_partner_test".to_string(),
            proposal_digest: "sha256:0123456789abcdef".to_string(),
            request: RegulatoryFragmentPanelRequest {
                fragments: vec![candidate],
                ..RegulatoryFragmentPanelRequest::default()
            },
            members: vec![RegulatoryFragmentPanelMember {
                member_id: "candidate_alone".to_string(),
                insert_length_bp: 120,
                instances: vec![RegulatoryFragmentResolvedInstance {
                    fragment_id: "candidate_a".to_string(),
                    role: RegulatoryFragmentRole::Candidate,
                    orientation: RegulatoryFragmentOrientation::Forward,
                    assembled_end_0based_exclusive: 120,
                    ..RegulatoryFragmentResolvedInstance::default()
                }],
                ..RegulatoryFragmentPanelMember::default()
            }],
            contrasts: vec![RegulatoryFragmentContrast {
                contrast_id: "standalone".to_string(),
                question: RegulatoryFragmentQuestion::StandaloneCandidate,
                left_member_id: "promoterless_control".to_string(),
                right_member_id: "candidate_alone".to_string(),
                interpretation: "Tests a candidate-associated reporter difference.".to_string(),
            }],
            evidence_dimensions: vec![RegulatoryFragmentEvidenceDimension {
                kind: RegulatoryFragmentEvidenceDimensionKind::PanelSequenceSimilarity,
                state: RegulatoryFragmentEvidenceState::Evaluated,
                assessment_id: "panel:similarity".to_string(),
                assessment_sha256: "sha256:evidence".to_string(),
                detail: "Similarity is context evidence.".to_string(),
                ..RegulatoryFragmentEvidenceDimension::default()
            }],
            nonclaims: vec!["This figure does not establish sufficiency.".to_string()],
            ..RegulatoryFragmentPanelPlan::default()
        };
        let svg = render_regulatory_fragment_panel_svg(&plan);
        assert!(svg.contains("data-gentle-schema=\"gentle.regulatory_fragment_panel_plan.v1\""));
        assert!(svg.contains("Selected construct geometry"));
        assert!(svg.contains("data-gentle-role=\"source-fragment\""));
        assert!(svg.contains(
            "data-gentle-role=\"source-fragment\" x=\"190.0\" y=\"120.0\" width=\"978.0\""
        ));
        assert!(svg.contains("data-gentle-role=\"construct-fragment\""));
        assert!(svg.contains("<title>candidate_a &gt; 120 bp</title>"));
        assert!(svg.contains("StandaloneCandidate"));
        assert!(svg.contains("Panel / vector similarity"));
        assert!(svg.contains("does not establish sufficiency"));
    }

    #[test]
    fn regulatory_fragment_panel_svg_uses_traceable_short_labels_for_long_member_ids() {
        let prefix = "a_long_tutorial_panel_name_with_shared_member_prefix";
        let left = format!("{prefix}_candidate_alone");
        let right = format!("{prefix}_reference_combination");
        let plan = RegulatoryFragmentPanelPlan {
            plan_id: prefix.to_string(),
            members: vec![
                RegulatoryFragmentPanelMember {
                    member_id: left.clone(),
                    ..Default::default()
                },
                RegulatoryFragmentPanelMember {
                    member_id: right.clone(),
                    ..Default::default()
                },
            ],
            contrasts: vec![RegulatoryFragmentContrast {
                question: RegulatoryFragmentQuestion::PartnerDependence,
                left_member_id: left.clone(),
                right_member_id: right.clone(),
                ..Default::default()
            }],
            ..Default::default()
        };
        let svg = render_regulatory_fragment_panel_svg(&plan);
        assert!(svg.contains("PartnerDependence: C1 vs C2</text>"));
        assert!(svg.contains("C1: candidate_alone</text>"));
        assert!(svg.contains("C2: reference_combination</text>"));
        assert!(svg.contains(&format!("<title>{left}</title>")));
        assert!(svg.contains(&format!("<title>{left} vs {right}: </title>")));
        assert!(svg.contains("width=\"154.0\" height=\"16\" overflow=\"hidden\""));
    }

    #[test]
    fn regulatory_fragment_panel_svg_is_byte_stable() {
        let plan = RegulatoryFragmentPanelPlan {
            plan_id: "stable".to_string(),
            ..RegulatoryFragmentPanelPlan::default()
        };
        assert_eq!(
            render_regulatory_fragment_panel_svg(&plan),
            render_regulatory_fragment_panel_svg(&plan)
        );
    }
}
