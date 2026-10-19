//! Read-only, content-bound bench dependency diagrams, not executed lineage.

use crate::digest_utils::{canonical_json_string, sha256_prefixed_bytes};
use crate::engine::{PlanningEstimate, PlanningStepModel, ScheduledStep, TimeRange};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};
use std::fmt::Write;

/// Wire schema for the read-only dependency projection.
pub const SCHEMA: &str = "gentle.planning_protocol_overview.v1";
#[cfg(test)]
mod tests;
const MAX_STEPS: usize = 128;
const NODE_WIDTH: u32 = 240;
const NODE_HEIGHT: u32 = 116;

/// Missing bench structure is unavailable, not a valid empty procedure.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum ProtocolOverviewStatus {
    Available,
    Unavailable,
}

/// Dependency geometry is a projection; it never creates operation/history nodes.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ProtocolNode {
    pub step: ScheduledStep,
    pub label_lines: Vec<String>,
    pub x: u32,
    pub y: u32,
    pub width: u32,
    pub height: u32,
    pub typical_critical_path: bool,
}

/// An explicitly recorded prerequisite, not inferred material transfer.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ProtocolEdge {
    pub from: String,
    pub to: String,
    pub typical_critical_path: bool,
    pub points: Vec<[u32; 2]>,
}

/// Portable diagram and SVG shared by GUI, shell, CLI and agent consumers.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ProtocolOverview {
    pub schema: String,
    pub routine_id: String,
    pub status: ProtocolOverviewStatus,
    pub reason: String,
    /// Canonical JSON digest of `[routine_id, estimate]`, not a signature.
    pub source_estimate_sha256: String,
    /// Canonical complete report digest with this field empty during hashing.
    pub content_sha256: String,
    pub width: u32,
    pub height: u32,
    pub nodes: Vec<ProtocolNode>,
    pub edges: Vec<ProtocolEdge>,
    pub typical_critical_path_ids: Option<Vec<String>>,
    pub nonclaims: Vec<String>,
    pub svg: String,
}

fn error(code: &str, detail: &str) -> String {
    serde_json::json!({"code":code,"detail":detail}).to_string()
}

fn valid_range(range: TimeRange) -> bool {
    [range.min, range.typical, range.max]
        .iter()
        .all(|value| value.is_finite() && *value >= 0.0)
        && range.min <= range.typical
        && range.typical <= range.max
}

fn label_lines(label: &str) -> Vec<String> {
    let mut lines = vec![String::new()];
    for word in label.split_whitespace() {
        let line = lines.last_mut().unwrap();
        if !line.is_empty() && line.chars().count() + word.chars().count() + 1 > 29 {
            if lines.len() == 2 {
                lines[1].push_str(" ...");
                break;
            }
            lines.push(String::new());
        }
        let line = lines.last_mut().unwrap();
        if !line.is_empty() {
            line.push(' ');
        }
        line.extend(word.chars().take(29));
    }
    lines
}

/// Construct from the exact shared estimate, without recomputing timings/policy.
pub fn build(routine_id: &str, estimate: &PlanningEstimate) -> Result<ProtocolOverview, String> {
    if [
        estimate.estimated_time_hours,
        estimate.estimated_cost,
        estimate.local_fit_score,
        estimate.composite_meta_score,
        estimate.stock_drawn_material_cost,
        estimate.missing_material_purchase_proxy_cost,
        estimate.machine_run_cost,
        estimate.routine_base_cost,
    ]
    .iter()
    .any(|value| !value.is_finite())
        || [
            estimate.hands_on_hours,
            estimate.unattended_hours,
            estimate.elapsed_hours,
        ]
        .into_iter()
        .flatten()
        .any(|range| !valid_range(range))
    {
        return Err(error(
            "invalid_protocol_overview_estimate",
            "Non-finite estimate or invalid scenario range cannot be content-bound",
        ));
    }
    let bytes = canonical_json_string(&(routine_id, estimate)).map_err(|err| err.to_string())?;
    let mut report = ProtocolOverview {
        schema: SCHEMA.into(),
        routine_id: routine_id.into(),
        status: ProtocolOverviewStatus::Unavailable,
        reason: "No bench-step model recorded; a legacy elapsed scalar is not a protocol.".into(),
        source_estimate_sha256: sha256_prefixed_bytes(bytes.as_bytes()),
        content_sha256: String::new(),
        width: 900,
        height: 200,
        nodes: vec![],
        edges: vec![],
        typical_critical_path_ids: None,
        nonclaims: vec![
            "Planned synthetic bench allowances; not executed work or validated laboratory instructions.".into(),
            "Arrows are dependencies, not verified material transfers or sequence compatibility.".into(),
            "Independent branches do not establish staffing or shared-equipment feasibility.".into(),
            "Material identities, amounts, equipment assignments and verification success are not bound by this model.".into(),
            "No scheduling, execution, supplier contact, ranking change or readiness approval.".into(),
        ],
        svg: String::new(),
    };
    if estimate.step_model == PlanningStepModel::Present {
        if estimate.schedule.is_empty() || estimate.schedule.len() > MAX_STEPS {
            return Err(error(
                "protocol_overview_step_budget",
                "Require 1..=128 recorded steps",
            ));
        }
        let mut pending = BTreeMap::new();
        for step in &estimate.schedule {
            if step.step_id.is_empty()
                || step.label.trim().is_empty()
                || !valid_range(step.hands_on_hours)
                || !valid_range(step.unattended_hours)
                || !valid_range(step.earliest_start_hours)
                || !valid_range(step.earliest_finish_hours)
                || step.source.kind != "synthetic_planning_default"
                || step.source.note.trim().is_empty()
                || pending.insert(step.step_id.clone(), step).is_some()
            {
                return Err(error("invalid_protocol_overview_step", &step.step_id));
            }
        }
        for step in pending.values() {
            let unique = step.depends_on.iter().collect::<BTreeSet<_>>();
            if unique.len() != step.depends_on.len()
                || unique.iter().any(|id| !pending.contains_key(*id))
            {
                return Err(error("invalid_protocol_overview_dependency", &step.step_id));
            }
        }
        let critical: Vec<String> = estimate
            .explanation
            .get("critical_paths")
            .and_then(|value| value.get(1))
            .cloned()
            .map(serde_json::from_value)
            .transpose()
            .map_err(|_| {
                error(
                    "invalid_protocol_overview_critical_path",
                    "Expected step IDs",
                )
            })?
            .unwrap_or_default();
        if critical.iter().any(|id| !pending.contains_key(id))
            || critical
                .windows(2)
                .any(|pair| !pending[&pair[1]].depends_on.contains(&pair[0]))
        {
            return Err(error(
                "invalid_protocol_overview_critical_path",
                "Unknown or disconnected step",
            ));
        }
        report.typical_critical_path_ids = (!critical.is_empty()).then(|| critical.clone());
        let mut levels = BTreeMap::<String, usize>::new();
        let mut rows = BTreeMap::<usize, u32>::new();
        while !pending.is_empty() {
            let Some(id) = pending
                .iter()
                .find(|(_, step)| step.depends_on.iter().all(|id| levels.contains_key(id)))
                .map(|(id, _)| id.clone())
            else {
                return Err(error(
                    "protocol_overview_dependency_cycle",
                    "No admissible next step",
                ));
            };
            let step = pending.remove(&id).unwrap();
            let level = step
                .depends_on
                .iter()
                .map(|id| levels[id] + 1)
                .max()
                .unwrap_or(0);
            let row = rows.entry(level).or_default();
            report.nodes.push(ProtocolNode {
                step: step.clone(),
                label_lines: label_lines(&step.label),
                x: 24 + level as u32 * 304,
                y: 100 + *row * 148,
                width: NODE_WIDTH,
                height: NODE_HEIGHT,
                typical_critical_path: critical.contains(&id),
            });
            *row += 1;
            levels.insert(id.clone(), level);
            let mut dependencies = step.depends_on.clone();
            dependencies.sort();
            for from in dependencies {
                report.edges.push(ProtocolEdge {
                    points: vec![],
                    typical_critical_path: critical
                        .windows(2)
                        .any(|pair| pair[0] == from && pair[1] == id),
                    from,
                    to: id.clone(),
                });
            }
        }
        report
            .nodes
            .sort_by(|a, b| a.step.step_id.cmp(&b.step.step_id));
        report
            .edges
            .sort_by(|a, b| (&a.from, &a.to).cmp(&(&b.from, &b.to)));
        report.width = report
            .nodes
            .iter()
            .map(|node| node.x + node.width + 24)
            .max()
            .unwrap()
            .max(900);
        report.height = report
            .nodes
            .iter()
            .map(|node| node.y + node.height + 80)
            .max()
            .unwrap();
        for (index, edge) in report.edges.iter_mut().enumerate() {
            let from = report
                .nodes
                .iter()
                .find(|node| node.step.step_id == edge.from)
                .unwrap();
            let to = report
                .nodes
                .iter()
                .find(|node| node.step.step_id == edge.to)
                .unwrap();
            let start = [from.x + from.width, from.y + from.height / 2];
            let end = [to.x, to.y + to.height / 2];
            edge.points = if to.x == from.x + 304 {
                let middle = (start[0] + end[0]) / 2;
                vec![start, [middle, start[1]], [middle, end[1]], end]
            } else {
                // Skip-column edges use gutters and a shared bottom channel,
                // never a vertical segment through an intervening step box.
                let channel = report.height - 64 + (index as u32 % 4) * 8;
                vec![
                    start,
                    [start[0] + 16, start[1]],
                    [start[0] + 16, channel],
                    [end[0] - 16, channel],
                    [end[0] - 16, end[1]],
                    end,
                ]
            };
        }
        report.status = ProtocolOverviewStatus::Available;
        report.reason = "Synthetic bench dependency overview; not executed lineage.".into();
    } else if !estimate.schedule.is_empty() {
        return Err(error(
            "protocol_overview_model_mismatch",
            "Absent model has recorded steps",
        ));
    }
    report.svg = render_svg(&report);
    report.content_sha256 = sha256_prefixed_bytes(
        canonical_json_string(&report)
            .map_err(|err| err.to_string())?
            .as_bytes(),
    );
    Ok(report)
}

/// Rebind portable presentation to its estimate, rather than trusting a wire digest.
pub fn from_explanation(output: &serde_json::Value) -> Result<ProtocolOverview, String> {
    let id = output["routine"]["routine_id"]
        .as_str()
        .ok_or_else(|| error("protocol_overview_unavailable", "Routine identity absent"))?;
    let estimate: PlanningEstimate = serde_json::from_value(output["planning"]["estimate"].clone())
        .map_err(|err| error("protocol_overview_unavailable", &err.to_string()))?;
    let recorded: ProtocolOverview =
        serde_json::from_value(output["planning"]["protocol_overview"].clone()).map_err(|_| {
            error(
                "protocol_overview_unavailable",
                "No readable recorded overview",
            )
        })?;
    let expected = build(id, &estimate)?;
    if serde_json::to_value(&recorded).unwrap() != serde_json::to_value(&expected).unwrap() {
        return Err(error(
            "protocol_overview_binding_mismatch",
            "Overview does not match routine and estimate",
        ));
    }
    Ok(recorded)
}

/// Identical tooltip/details content for the graph and portable SVG hover text.
pub fn node_details(node: &ProtocolNode) -> String {
    let step = &node.step;
    format!(
        "{}\n{}\nHands-on (min/typical/max): {:.2}/{:.2}/{:.2} h\nUnattended (min/typical/max): {:.2}/{:.2}/{:.2} h\nEarliest start (min/typical/max): {:.2}/{:.2}/{:.2} h\nEarliest finish (min/typical/max): {:.2}/{:.2}/{:.2} h\nSource: {}\n{}\nPrerequisites: {}",
        step.step_id,
        step.label,
        step.hands_on_hours.min,
        step.hands_on_hours.typical,
        step.hands_on_hours.max,
        step.unattended_hours.min,
        step.unattended_hours.typical,
        step.unattended_hours.max,
        step.earliest_start_hours.min,
        step.earliest_start_hours.typical,
        step.earliest_start_hours.max,
        step.earliest_finish_hours.min,
        step.earliest_finish_hours.typical,
        step.earliest_finish_hours.max,
        step.source.kind,
        step.source.note,
        if step.depends_on.is_empty() {
            "none recorded".into()
        } else {
            step.depends_on.join(", ")
        }
    )
}

fn xml(text: &str) -> String {
    text.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
        .replace('\'', "&apos;")
}

/// Scientific legend is identical in the GUI and the printable projection.
pub fn legend(report: &ProtocolOverview) -> &'static str {
    if report.typical_critical_path_ids.is_some() {
        "Magenta: recorded typical critical path. Arrows: prerequisites, not material transfer. Columns: dependency depth, not time."
    } else {
        "Typical critical path unavailable. Arrows: prerequisites, not material transfer. Columns: dependency depth, not time."
    }
}

fn render_svg(report: &ProtocolOverview) -> String {
    let mut svg = format!(
        r##"<svg xmlns="http://www.w3.org/2000/svg" width="{}" height="{}" viewBox="0 0 {} {}" role="img" data-schema="{}" data-source-estimate-sha256="{}">
<title>Planned bench dependencies: {}</title>
<desc>{}</desc>
<defs><marker id="arrow" markerWidth="8" markerHeight="8" refX="7" refY="4" orient="auto"><path d="M0,0 L8,4 L0,8 Z" fill="#555555"/></marker></defs>
<rect width="100%" height="100%" fill="white"/>
<g font-family="sans-serif" fill="#222222"><text x="24" y="28" font-size="18">Planned bench dependencies</text><text x="24" y="50" font-size="12">{}</text><text x="24" y="72" font-size="12">{}</text>
"##,
        report.width,
        report.height,
        report.width,
        report.height,
        SCHEMA,
        report.source_estimate_sha256,
        xml(&report.routine_id),
        xml(&report.reason),
        xml(&report.reason),
        xml(legend(report))
    );
    let metadata = serde_json::json!({"schema":SCHEMA,"routine_id":report.routine_id,
        "source_estimate_sha256":report.source_estimate_sha256,"nonclaims":report.nonclaims});
    let _ = writeln!(
        svg,
        "<metadata>{}</metadata>",
        xml(&canonical_json_string(&metadata).unwrap())
    );
    for edge in &report.edges {
        let path = edge
            .points
            .iter()
            .enumerate()
            .map(|(index, point)| {
                format!(
                    "{}{},{}",
                    if index == 0 { "M" } else { "L" },
                    point[0],
                    point[1]
                )
            })
            .collect::<Vec<_>>()
            .join(" ");
        let color = if edge.typical_critical_path {
            "#b12b62"
        } else {
            "#777777"
        };
        let _ = writeln!(
            svg,
            r#"<path data-from="{}" data-to="{}" d="{path}" fill="none" stroke="{color}" stroke-width="2" marker-end="url(#arrow)"/>"#,
            xml(&edge.from),
            xml(&edge.to)
        );
    }
    for node in &report.nodes {
        let color = if node.typical_critical_path {
            "#b12b62"
        } else {
            "#777777"
        };
        let _ = writeln!(
            svg,
            r##"<g data-step-id="{}"><title>{}</title><rect x="{}" y="{}" width="{}" height="{}" rx="4" fill="#fafafa" stroke="{color}" stroke-width="2"/>"##,
            xml(&node.step.step_id),
            xml(&node_details(node)),
            node.x,
            node.y,
            node.width,
            node.height
        );
        for (index, line) in node.label_lines.iter().enumerate() {
            let _ = writeln!(
                svg,
                r#"<text x="{}" y="{}" font-size="13">{}</text>"#,
                node.x + 12,
                node.y + 22 + index as u32 * 18,
                xml(line)
            );
        }
        let _ = writeln!(
            svg,
            r##"<text x="{}" y="{}" font-size="12" fill="#24734b">Hands-on {:.2} h</text><text x="{}" y="{}" font-size="12" fill="#906400">Unattended {:.2} h</text></g>"##,
            node.x + 12,
            node.y + 78,
            node.step.hands_on_hours.typical,
            node.x + 12,
            node.y + 98,
            node.step.unattended_hours.typical
        );
    }
    if report.nodes.is_empty() {
        svg.push_str(
            "<text x=\"24\" y=\"120\" font-size=\"14\">Bench-step overview unavailable</text>\n",
        );
    }
    let _ = writeln!(
        svg,
        r#"<text x="24" y="{}" font-size="12">Synthetic allowances; no validated protocol, execution, staffing or material-flow claim.</text></g></svg>"#,
        report.height - 24
    );
    svg
}
