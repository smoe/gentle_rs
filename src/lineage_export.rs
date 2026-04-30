//! Lineage graph export and serialization utilities.

use crate::{
    engine::{GentleEngine, LineageMacroPortBinding, Operation, OperationRecord, ProjectState},
    gibson_planning::GibsonAssemblyPlan,
};
use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use svg::Document;
use svg::node::element::{Circle, Line, Polygon, Rectangle, Text};

const MIN_H: f32 = 180.0;
const MIN_W: f32 = 960.0;
const TOP_CONTENT_Y: f32 = 120.0;
const LEFT_CONTENT_X: f32 = 110.0;
const RIGHT_PADDING: f32 = 340.0;
const BOTTOM_PADDING: f32 = 72.0;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum LineageSvgNodeKind {
    Sequence,
    Pool,
    Arrangement,
    Macro,
    Analysis,
    OperationHub,
}

#[derive(Clone, Debug, PartialEq)]
pub struct LineageSvgNode {
    pub node_id: String,
    pub title: String,
    pub subtitle: String,
    pub kind: LineageSvgNodeKind,
    pub x: f32,
    pub y: f32,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct LineageSvgEdge {
    pub from_node_id: String,
    pub to_node_id: String,
    pub label: String,
}

fn lineage_canvas_width(pos_by_node: &HashMap<String, (f32, f32)>) -> f32 {
    let max_content_x = pos_by_node
        .values()
        .map(|(x, _)| *x)
        .fold(LEFT_CONTENT_X, f32::max);
    (max_content_x + RIGHT_PADDING).max(MIN_W)
}

fn lineage_canvas_height(pos_by_node: &HashMap<String, (f32, f32)>) -> f32 {
    let max_content_y = pos_by_node
        .values()
        .map(|(_, y)| *y)
        .fold(TOP_CONTENT_Y, f32::max);
    (max_content_y + BOTTOM_PADDING).max(MIN_H)
}

pub fn export_projected_lineage_svg(
    title: &str,
    nodes: &[LineageSvgNode],
    edges: &[LineageSvgEdge],
) -> String {
    let node_by_id: HashMap<String, &LineageSvgNode> = nodes
        .iter()
        .map(|node| (node.node_id.clone(), node))
        .collect();
    let pos_by_node: HashMap<String, (f32, f32)> = nodes
        .iter()
        .map(|node| (node.node_id.clone(), (node.x, node.y)))
        .collect();
    let canvas_width = lineage_canvas_width(&pos_by_node);
    let canvas_height = lineage_canvas_height(&pos_by_node);

    let mut doc = Document::new()
        .set("viewBox", (0, 0, canvas_width, canvas_height))
        .set("width", canvas_width)
        .set("height", canvas_height)
        .set("style", "background:#ffffff");

    doc = doc.add(
        Text::new(title.to_string())
            .set("x", 24)
            .set("y", 34)
            .set("font-family", "Helvetica, Arial, sans-serif")
            .set("font-size", 24)
            .set("fill", "#202020"),
    );

    for edge in edges {
        let Some((fx, fy)) = pos_by_node.get(&edge.from_node_id).copied() else {
            continue;
        };
        let Some((tx, ty)) = pos_by_node.get(&edge.to_node_id).copied() else {
            continue;
        };
        doc = doc.add(
            Line::new()
                .set("x1", fx)
                .set("y1", fy)
                .set("x2", tx)
                .set("y2", ty)
                .set("stroke", "#8a8a8a")
                .set("stroke-width", 1.2),
        );
        let suppress_edge_label = node_by_id
            .get(&edge.from_node_id)
            .into_iter()
            .chain(node_by_id.get(&edge.to_node_id))
            .any(|node| node.kind == LineageSvgNodeKind::OperationHub);
        if !suppress_edge_label && !edge.label.trim().is_empty() {
            let mx = (fx + tx) * 0.5;
            let my = (fy + ty) * 0.5 - 6.0;
            doc = doc
                .add(
                    Rectangle::new()
                        .set("x", mx - 52.0)
                        .set("y", my - 12.0)
                        .set("width", 104)
                        .set("height", 16)
                        .set("fill", "#f5f5f5")
                        .set("stroke", "#e0e0e0")
                        .set("rx", 2),
                )
                .add(
                    Text::new(edge.label.clone())
                        .set("x", mx)
                        .set("y", my)
                        .set("text-anchor", "middle")
                        .set("dominant-baseline", "middle")
                        .set("font-family", "Helvetica, Arial, sans-serif")
                        .set("font-size", 10)
                        .set("fill", "#222222"),
                );
        }
    }

    for node in nodes {
        let x = node.x;
        let y = node.y;
        match node.kind {
            LineageSvgNodeKind::Sequence => {
                doc = doc.add(
                    Circle::new()
                        .set("cx", x)
                        .set("cy", y)
                        .set("r", 16)
                        .set("fill", "#5a8cd2")
                        .set("stroke", "#3b6aaa")
                        .set("stroke-width", 1),
                );
            }
            LineageSvgNodeKind::Pool => {
                let points = format!(
                    "{},{} {},{} {},{} {},{}",
                    x,
                    y - 18.0,
                    x + 18.0,
                    y,
                    x,
                    y + 18.0,
                    x - 18.0,
                    y
                );
                doc = doc.add(
                    Polygon::new()
                        .set("points", points)
                        .set("fill", "#b47846")
                        .set("stroke", "#a05f2b")
                        .set("stroke-width", 1),
                );
            }
            LineageSvgNodeKind::Arrangement => {
                doc = doc.add(
                    Rectangle::new()
                        .set("x", x - 36.0)
                        .set("y", y - 14.0)
                        .set("width", 72)
                        .set("height", 28)
                        .set("rx", 5)
                        .set("fill", "#6c9a7a")
                        .set("stroke", "#537563")
                        .set("stroke-width", 1),
                );
            }
            LineageSvgNodeKind::Macro => {
                doc = doc.add(
                    Rectangle::new()
                        .set("x", x - 40.0)
                        .set("y", y - 15.0)
                        .set("width", 80)
                        .set("height", 30)
                        .set("rx", 4)
                        .set("fill", "#62626c")
                        .set("stroke", "#4f4f58")
                        .set("stroke-width", 1),
                );
            }
            LineageSvgNodeKind::Analysis => {
                doc = doc.add(
                    Rectangle::new()
                        .set("x", x - 31.0)
                        .set("y", y - 14.0)
                        .set("width", 62)
                        .set("height", 28)
                        .set("rx", 4)
                        .set("fill", "#8c62ac")
                        .set("stroke", "#704a8d")
                        .set("stroke-width", 1),
                );
            }
            LineageSvgNodeKind::OperationHub => {
                doc = doc.add(
                    Rectangle::new()
                        .set("x", x - 37.0)
                        .set("y", y - 15.0)
                        .set("width", 74)
                        .set("height", 30)
                        .set("rx", 4)
                        .set("fill", "#3e8e7c")
                        .set("stroke", "#2d6f61")
                        .set("stroke-width", 1),
                );
                doc = doc.add(
                    Text::new(node.title.clone())
                        .set("x", x)
                        .set("y", y)
                        .set("text-anchor", "middle")
                        .set("dominant-baseline", "middle")
                        .set("font-family", "Helvetica, Arial, sans-serif")
                        .set("font-size", 12)
                        .set("fill", "#ffffff"),
                );
                continue;
            }
        }

        doc = doc
            .add(
                Text::new(node.title.clone())
                    .set("x", x + 24.0)
                    .set("y", y - 2.0)
                    .set("font-family", "Helvetica, Arial, sans-serif")
                    .set("font-size", 12)
                    .set("fill", "#101010"),
            )
            .add(
                Text::new(node.subtitle.clone())
                    .set("x", x + 24.0)
                    .set("y", y + 12.0)
                    .set("font-family", "Helvetica, Arial, sans-serif")
                    .set("font-size", 10)
                    .set("fill", "#222222"),
            );
    }

    doc.to_string()
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum EngineLineageRenderRowKind {
    Sequence,
    Arrangement,
    Macro,
    Analysis,
    OperationHub,
}

#[derive(Clone, Debug)]
struct EngineLineageRenderRow {
    node_id: String,
    seq_id: String,
    display_name: String,
    created_by_op: String,
    created_at: u128,
    kind: EngineLineageRenderRowKind,
    length: usize,
    circular: bool,
    pool_size: usize,
    arrangement_id: Option<String>,
    arrangement_mode: Option<String>,
    lane_count: usize,
    analysis_kind: Option<String>,
    analysis_status: Option<String>,
    analysis_reference_seq_id: Option<String>,
    analysis_read_count: Option<usize>,
    analysis_trace_count: Option<usize>,
    analysis_target_count: Option<usize>,
    analysis_variant_count: Option<usize>,
    macro_instance_id: Option<String>,
    macro_op_count: usize,
}

fn operation_variant_name(op: &Operation) -> String {
    serde_json::to_value(op)
        .ok()
        .and_then(|value| value.as_object().and_then(|obj| obj.keys().next().cloned()))
        .unwrap_or_else(|| "UnknownOperation".to_string())
}

fn humanize_operation_variant_name(raw: &str) -> String {
    let mut out = String::with_capacity(raw.len() + 8);
    let mut prev_is_lower_or_digit = false;
    for ch in raw.chars() {
        let is_upper = ch.is_ascii_uppercase();
        if is_upper && prev_is_lower_or_digit && !out.is_empty() {
            out.push(' ');
        }
        out.push(ch);
        prev_is_lower_or_digit = ch.is_ascii_lowercase() || ch.is_ascii_digit();
    }
    out
}

fn summarize_gibson_operation(plan_json: &str) -> String {
    let parsed: Result<GibsonAssemblyPlan, _> = serde_json::from_str(plan_json);
    match parsed {
        Ok(plan) => {
            let insert_seq_ids = if plan.fragments.is_empty() {
                "-".to_string()
            } else {
                plan.fragments
                    .iter()
                    .map(|fragment| fragment.seq_id.clone())
                    .collect::<Vec<_>>()
                    .join(",")
            };
            let opening_mode = if plan.destination.opening.mode.trim().is_empty() {
                "-".to_string()
            } else {
                plan.destination.opening.mode.trim().to_string()
            };
            let output = if plan.product.output_id_hint.trim().is_empty() {
                "-".to_string()
            } else {
                plan.product.output_id_hint.trim().to_string()
            };
            format!(
                "Gibson cloning: destination={}, inserts={}, opening={}, output={}",
                plan.destination.seq_id, insert_seq_ids, opening_mode, output
            )
        }
        Err(_) => "Gibson cloning".to_string(),
    }
}

fn summarize_operation_for_lineage(op: &Operation) -> String {
    match op {
        Operation::ApplyGibsonAssemblyPlan { plan_json } => summarize_gibson_operation(plan_json),
        Operation::ConfirmConstructReads {
            expected_seq_id,
            baseline_seq_id,
            read_seq_ids,
            trace_ids,
            report_id,
            ..
        } => format!(
            "Sequencing confirmation: expected={}, baseline={}, reads={}, traces={}, report_id={}",
            expected_seq_id,
            baseline_seq_id.as_deref().unwrap_or("-"),
            read_seq_ids.len(),
            trace_ids.len(),
            report_id.as_deref().unwrap_or("-"),
        ),
        _ => humanize_operation_variant_name(&operation_variant_name(op)),
    }
}

fn infer_gibson_like_operation_ids_from_state(state: &ProjectState) -> HashSet<String> {
    let mut output_seq_ids_by_op: HashMap<String, Vec<String>> = HashMap::new();
    for node in state.lineage.nodes.values() {
        let Some(op_id) = node.created_by_op.as_ref() else {
            continue;
        };
        output_seq_ids_by_op
            .entry(op_id.clone())
            .or_default()
            .push(node.seq_id.clone());
    }

    let mut parent_nodes_by_op: HashMap<String, HashSet<String>> = HashMap::new();
    let mut child_nodes_by_op: HashMap<String, HashSet<String>> = HashMap::new();
    for edge in &state.lineage.edges {
        parent_nodes_by_op
            .entry(edge.op_id.clone())
            .or_default()
            .insert(edge.from_node_id.clone());
        child_nodes_by_op
            .entry(edge.op_id.clone())
            .or_default()
            .insert(edge.to_node_id.clone());
    }

    output_seq_ids_by_op
        .into_iter()
        .filter_map(|(op_id, output_seq_ids)| {
            let has_left_primer = output_seq_ids
                .iter()
                .any(|seq_id| seq_id.ends_with("_left_insert_primer"));
            let has_right_primer = output_seq_ids
                .iter()
                .any(|seq_id| seq_id.ends_with("_right_insert_primer"));
            if !has_left_primer || !has_right_primer {
                return None;
            }
            let parent_count = parent_nodes_by_op
                .get(&op_id)
                .map(|nodes| nodes.len())
                .unwrap_or(0);
            let child_count = child_nodes_by_op
                .get(&op_id)
                .map(|nodes| nodes.len())
                .unwrap_or(0);
            (parent_count >= 2 && child_count >= 2).then_some(op_id)
        })
        .collect()
}

fn lineage_layout_positions(
    order_by_layer: &BTreeMap<usize, Vec<String>>,
) -> HashMap<String, usize> {
    let mut positions = HashMap::new();
    for nodes in order_by_layer.values() {
        for (index, node_id) in nodes.iter().enumerate() {
            positions.insert(node_id.clone(), index);
        }
    }
    positions
}

fn compute_lineage_dag_layout(
    rows: &[EngineLineageRenderRow],
    edges: &[(String, String, String)],
) -> HashMap<String, (usize, usize)> {
    if rows.is_empty() {
        return HashMap::new();
    }

    let mut row_index: HashMap<String, usize> = HashMap::new();
    for (index, row) in rows.iter().enumerate() {
        row_index.insert(row.node_id.clone(), index);
    }

    let mut parents_by_node: HashMap<String, Vec<String>> = HashMap::new();
    let mut children_by_node: HashMap<String, Vec<String>> = HashMap::new();
    let mut indegree_by_node: HashMap<String, usize> = HashMap::new();
    for row in rows {
        parents_by_node.insert(row.node_id.clone(), Vec::new());
        children_by_node.insert(row.node_id.clone(), Vec::new());
        indegree_by_node.insert(row.node_id.clone(), 0);
    }

    let mut seen_edges: HashSet<(String, String)> = HashSet::new();
    for (from_node, to_node, _op_id) in edges {
        if !row_index.contains_key(from_node) || !row_index.contains_key(to_node) {
            continue;
        }
        if !seen_edges.insert((from_node.clone(), to_node.clone())) {
            continue;
        }
        if let Some(children) = children_by_node.get_mut(from_node) {
            children.push(to_node.clone());
        }
        if let Some(parents) = parents_by_node.get_mut(to_node) {
            parents.push(from_node.clone());
        }
        if let Some(indegree) = indegree_by_node.get_mut(to_node) {
            *indegree = indegree.saturating_add(1);
        }
    }

    for parents in parents_by_node.values_mut() {
        parents.sort_by_key(|node_id| row_index.get(node_id).copied().unwrap_or(usize::MAX));
    }
    for children in children_by_node.values_mut() {
        children.sort_by_key(|node_id| row_index.get(node_id).copied().unwrap_or(usize::MAX));
    }

    let mut ready: Vec<String> = indegree_by_node
        .iter()
        .filter_map(|(node_id, indegree)| (*indegree == 0).then(|| node_id.clone()))
        .collect();
    ready.sort_by_key(|node_id| row_index.get(node_id).copied().unwrap_or(usize::MAX));

    let mut topo_order: Vec<String> = Vec::with_capacity(rows.len());
    let mut topo_seen: HashSet<String> = HashSet::with_capacity(rows.len());
    while !ready.is_empty() {
        let node_id = ready.remove(0);
        if !topo_seen.insert(node_id.clone()) {
            continue;
        }
        topo_order.push(node_id.clone());
        if let Some(children) = children_by_node.get(&node_id) {
            for child_id in children {
                if let Some(indegree) = indegree_by_node.get_mut(child_id)
                    && *indegree > 0
                {
                    *indegree -= 1;
                    if *indegree == 0 {
                        ready.push(child_id.clone());
                    }
                }
            }
        }
        ready.sort_by_key(|candidate| row_index.get(candidate).copied().unwrap_or(usize::MAX));
    }

    for row in rows {
        if topo_seen.insert(row.node_id.clone()) {
            topo_order.push(row.node_id.clone());
        }
    }

    let mut layer_by_node: HashMap<String, usize> = HashMap::new();
    for node_id in &topo_order {
        let layer = parents_by_node
            .get(node_id)
            .map(|parents| {
                parents
                    .iter()
                    .filter_map(|parent_id| layer_by_node.get(parent_id).copied())
                    .max()
                    .map(|max_parent_layer| max_parent_layer + 1)
                    .unwrap_or(0)
            })
            .unwrap_or(0);
        layer_by_node.insert(node_id.clone(), layer);
    }

    let mut order_by_layer: BTreeMap<usize, Vec<String>> = BTreeMap::new();
    for node_id in &topo_order {
        let layer = layer_by_node.get(node_id).copied().unwrap_or(0);
        order_by_layer
            .entry(layer)
            .or_default()
            .push(node_id.clone());
    }
    for nodes in order_by_layer.values_mut() {
        nodes.sort_by_key(|node_id| row_index.get(node_id).copied().unwrap_or(usize::MAX));
    }

    let max_layer = order_by_layer.keys().copied().max().unwrap_or(0);
    let barycenter = |neighbors: &[String], positions: &HashMap<String, usize>| -> Option<f32> {
        let mut sum = 0.0f32;
        let mut count = 0usize;
        for node_id in neighbors {
            if let Some(pos) = positions.get(node_id) {
                sum += *pos as f32;
                count += 1;
            }
        }
        (count > 0).then(|| sum / count as f32)
    };

    for _ in 0..6 {
        for layer in 1..=max_layer {
            let positions = lineage_layout_positions(&order_by_layer);
            let Some(mut nodes) = order_by_layer.remove(&layer) else {
                continue;
            };
            nodes.sort_by(|left, right| {
                let left_score = parents_by_node
                    .get(left)
                    .and_then(|parents| barycenter(parents, &positions));
                let right_score = parents_by_node
                    .get(right)
                    .and_then(|parents| barycenter(parents, &positions));
                let fallback = row_index
                    .get(left)
                    .copied()
                    .unwrap_or(usize::MAX)
                    .cmp(&row_index.get(right).copied().unwrap_or(usize::MAX));
                match (left_score, right_score) {
                    (Some(l), Some(r)) => l
                        .partial_cmp(&r)
                        .unwrap_or(std::cmp::Ordering::Equal)
                        .then(fallback),
                    (Some(_), None) => std::cmp::Ordering::Less,
                    (None, Some(_)) => std::cmp::Ordering::Greater,
                    (None, None) => fallback,
                }
            });
            order_by_layer.insert(layer, nodes);
        }

        if max_layer == 0 {
            break;
        }

        for layer in (0..max_layer).rev() {
            let positions = lineage_layout_positions(&order_by_layer);
            let Some(mut nodes) = order_by_layer.remove(&layer) else {
                continue;
            };
            nodes.sort_by(|left, right| {
                let left_score = children_by_node
                    .get(left)
                    .and_then(|children| barycenter(children, &positions));
                let right_score = children_by_node
                    .get(right)
                    .and_then(|children| barycenter(children, &positions));
                let fallback = row_index
                    .get(left)
                    .copied()
                    .unwrap_or(usize::MAX)
                    .cmp(&row_index.get(right).copied().unwrap_or(usize::MAX));
                match (left_score, right_score) {
                    (Some(l), Some(r)) => l
                        .partial_cmp(&r)
                        .unwrap_or(std::cmp::Ordering::Equal)
                        .then(fallback),
                    (Some(_), None) => std::cmp::Ordering::Less,
                    (None, Some(_)) => std::cmp::Ordering::Greater,
                    (None, None) => fallback,
                }
            });
            order_by_layer.insert(layer, nodes);
        }
    }

    let mut layout_by_node: HashMap<String, (usize, usize)> = HashMap::new();
    for (layer, nodes) in &order_by_layer {
        for (rank, node_id) in nodes.iter().enumerate() {
            layout_by_node.insert(node_id.clone(), (*layer, rank));
        }
    }

    for row in rows {
        let fallback_rank = row_index.get(&row.node_id).copied().unwrap_or(0);
        layout_by_node
            .entry(row.node_id.clone())
            .or_insert((0, fallback_rank));
    }

    layout_by_node
}

fn sequence_nodes_for_binding(
    state: &ProjectState,
    binding: &LineageMacroPortBinding,
) -> Vec<String> {
    let mut node_ids: Vec<String> = vec![];
    let mut seen: HashSet<String> = HashSet::new();
    match binding.kind.trim().to_ascii_lowercase().as_str() {
        "sequence" => {
            for seq_id in &binding.values {
                if let Some(node_id) = state.lineage.seq_to_node.get(seq_id)
                    && seen.insert(node_id.clone())
                {
                    node_ids.push(node_id.clone());
                }
            }
        }
        "container" => {
            for container_id in &binding.values {
                let Some(container) = state.container_state.containers.get(container_id) else {
                    continue;
                };
                for seq_id in &container.members {
                    if let Some(node_id) = state.lineage.seq_to_node.get(seq_id)
                        && seen.insert(node_id.clone())
                    {
                        node_ids.push(node_id.clone());
                    }
                }
            }
        }
        _ => {}
    }
    node_ids
}

fn project_lineage_operation_hubs(
    rows: &[EngineLineageRenderRow],
    edges: &[(String, String, String)],
    hub_op_ids: &HashSet<String>,
) -> (Vec<EngineLineageRenderRow>, Vec<(String, String, String)>) {
    let valid_node_ids: HashSet<String> = rows.iter().map(|row| row.node_id.clone()).collect();
    let row_by_node: HashMap<String, EngineLineageRenderRow> = rows
        .iter()
        .map(|row| (row.node_id.clone(), row.clone()))
        .collect();
    let mut op_parents: HashMap<String, BTreeSet<String>> = HashMap::new();
    let mut op_children: HashMap<String, BTreeSet<String>> = HashMap::new();
    for (from_node, to_node, op_id) in edges {
        if !hub_op_ids.contains(op_id) {
            continue;
        }
        if !valid_node_ids.contains(from_node) || !valid_node_ids.contains(to_node) {
            continue;
        }
        op_parents
            .entry(op_id.clone())
            .or_default()
            .insert(from_node.clone());
        op_children
            .entry(op_id.clone())
            .or_default()
            .insert(to_node.clone());
    }

    let hubbed_ops: Vec<String> = hub_op_ids
        .iter()
        .filter(|op_id| {
            op_parents.get(*op_id).map(|rows| rows.len()).unwrap_or(0) >= 2
                && op_children.get(*op_id).map(|rows| rows.len()).unwrap_or(0) >= 2
        })
        .cloned()
        .collect();
    if hubbed_ops.is_empty() {
        return (rows.to_vec(), edges.to_vec());
    }

    let hubbed_op_ids: HashSet<String> = hubbed_ops.iter().cloned().collect();
    let mut out_rows = rows.to_vec();
    let mut out_edges: Vec<(String, String, String)> = vec![];
    let mut seen_edges: HashSet<(String, String, String)> = HashSet::new();

    for op_id in hubbed_ops {
        let parents = op_parents
            .get(&op_id)
            .cloned()
            .unwrap_or_default()
            .into_iter()
            .collect::<Vec<_>>();
        let children = op_children
            .get(&op_id)
            .cloned()
            .unwrap_or_default()
            .into_iter()
            .collect::<Vec<_>>();
        let hub_node_id = format!("operation:{op_id}");
        let created_at = parents
            .iter()
            .filter_map(|node_id| row_by_node.get(node_id))
            .map(|row| row.created_at)
            .max()
            .unwrap_or(0);
        out_rows.push(EngineLineageRenderRow {
            node_id: hub_node_id.clone(),
            seq_id: op_id.clone(),
            display_name: "Gibson cloning".to_string(),
            created_by_op: op_id.clone(),
            created_at,
            kind: EngineLineageRenderRowKind::OperationHub,
            length: 0,
            circular: false,
            pool_size: 0,
            arrangement_id: None,
            arrangement_mode: None,
            lane_count: 0,
            analysis_kind: None,
            analysis_status: None,
            analysis_reference_seq_id: None,
            analysis_read_count: None,
            analysis_trace_count: None,
            analysis_target_count: None,
            analysis_variant_count: None,
            macro_instance_id: None,
            macro_op_count: 0,
        });

        let inbound_op_id = format!("{op_id}::hub_in");
        let outbound_op_id = format!("{op_id}::hub_out");
        for parent in parents {
            let key = (parent, hub_node_id.clone(), inbound_op_id.clone());
            if seen_edges.insert(key.clone()) {
                out_edges.push(key);
            }
        }
        for child in children {
            let key = (hub_node_id.clone(), child, outbound_op_id.clone());
            if seen_edges.insert(key.clone()) {
                out_edges.push(key);
            }
        }
    }

    for edge in edges {
        if hubbed_op_ids.contains(&edge.2) {
            continue;
        }
        if seen_edges.insert(edge.clone()) {
            out_edges.push(edge.clone());
        }
    }

    out_rows.sort_by(|a, b| {
        a.created_at
            .cmp(&b.created_at)
            .then(a.node_id.cmp(&b.node_id))
    });
    (out_rows, out_edges)
}

fn lineage_svg_node_kind(row: &EngineLineageRenderRow) -> LineageSvgNodeKind {
    match row.kind {
        EngineLineageRenderRowKind::Sequence if row.pool_size > 1 => LineageSvgNodeKind::Pool,
        EngineLineageRenderRowKind::Sequence => LineageSvgNodeKind::Sequence,
        EngineLineageRenderRowKind::Arrangement => LineageSvgNodeKind::Arrangement,
        EngineLineageRenderRowKind::Macro => LineageSvgNodeKind::Macro,
        EngineLineageRenderRowKind::Analysis => LineageSvgNodeKind::Analysis,
        EngineLineageRenderRowKind::OperationHub => LineageSvgNodeKind::OperationHub,
    }
}

fn lineage_svg_node_title(row: &EngineLineageRenderRow) -> String {
    let display = row.display_name.trim();
    if !display.is_empty() {
        display.to_string()
    } else {
        row.seq_id.clone()
    }
}

fn lineage_svg_node_subtitle(row: &EngineLineageRenderRow) -> String {
    match row.kind {
        EngineLineageRenderRowKind::Sequence if row.pool_size > 1 => {
            format!(
                "{} | pool n={} | {} bp",
                row.seq_id, row.pool_size, row.length
            )
        }
        EngineLineageRenderRowKind::Sequence => {
            let topology = if row.circular { "circular" } else { "linear" };
            format!("{} ({} bp, {topology})", row.seq_id, row.length)
        }
        EngineLineageRenderRowKind::Arrangement => format!(
            "{} | {} | lanes={}",
            row.arrangement_id.as_deref().unwrap_or(&row.seq_id),
            row.arrangement_mode.as_deref().unwrap_or("-"),
            row.lane_count
        ),
        EngineLineageRenderRowKind::Macro => format!(
            "{} | ops={}",
            row.macro_instance_id.as_deref().unwrap_or(&row.seq_id),
            row.macro_op_count
        ),
        EngineLineageRenderRowKind::Analysis => match row.analysis_kind.as_deref() {
            Some("sequencing_confirmation") => {
                let baseline = row.analysis_reference_seq_id.as_deref().unwrap_or("-");
                format!(
                    "status={} | baseline={} | reads={} | traces={} | targets={} | variants={}",
                    row.analysis_status.as_deref().unwrap_or("-"),
                    baseline,
                    row.analysis_read_count.unwrap_or(0),
                    row.analysis_trace_count.unwrap_or(0),
                    row.analysis_target_count.unwrap_or(0),
                    row.analysis_variant_count.unwrap_or(0)
                )
            }
            _ => "analysis".to_string(),
        },
        EngineLineageRenderRowKind::OperationHub => format!("op={}", row.created_by_op),
    }
}

pub fn build_lineage_svg_graph(
    state: &ProjectState,
    operation_log: &[OperationRecord],
) -> (Vec<LineageSvgNode>, Vec<LineageSvgEdge>) {
    let mut op_created_count: HashMap<String, usize> = HashMap::new();
    let mut op_label_by_id: HashMap<String, String> = HashMap::new();
    let mut sequencing_confirmation_op_by_report_id: HashMap<
        String,
        (String, String, Option<String>),
    > = HashMap::new();
    let mut hub_op_ids = infer_gibson_like_operation_ids_from_state(state);
    let mut individually_rendered_multi_output_ops = hub_op_ids.clone();
    for record in operation_log {
        op_created_count.insert(
            record.result.op_id.clone(),
            record.result.created_seq_ids.len(),
        );
        op_label_by_id.insert(
            record.result.op_id.clone(),
            summarize_operation_for_lineage(&record.op),
        );
        if matches!(record.op, Operation::ApplyGibsonAssemblyPlan { .. }) {
            hub_op_ids.insert(record.result.op_id.clone());
            individually_rendered_multi_output_ops.insert(record.result.op_id.clone());
        }
        if let Some(report) = record.result.sequencing_confirmation_report.as_ref() {
            sequencing_confirmation_op_by_report_id.insert(
                report.report_id.clone(),
                (
                    record.result.op_id.clone(),
                    report.expected_seq_id.clone(),
                    report.baseline_seq_id.clone(),
                ),
            );
        }
    }
    for op_id in &hub_op_ids {
        op_label_by_id
            .entry(op_id.clone())
            .or_insert_with(|| "Gibson cloning".to_string());
    }

    let mut rows: Vec<EngineLineageRenderRow> = state
        .lineage
        .nodes
        .values()
        .map(|node| {
            let display_name = state
                .sequences
                .get(&node.seq_id)
                .and_then(|dna| dna.name().clone())
                .unwrap_or_else(|| node.seq_id.clone());
            let (length, circular) = state
                .sequences
                .get(&node.seq_id)
                .map(|dna| (dna.len(), dna.is_circular()))
                .unwrap_or((0, false));
            let pool_size = node
                .created_by_op
                .as_ref()
                .and_then(|op_id| {
                    (!individually_rendered_multi_output_ops.contains(op_id))
                        .then(|| op_created_count.get(op_id).copied())
                        .flatten()
                })
                .unwrap_or(1);
            EngineLineageRenderRow {
                node_id: node.node_id.clone(),
                seq_id: node.seq_id.clone(),
                display_name,
                created_by_op: node
                    .created_by_op
                    .clone()
                    .unwrap_or_else(|| "-".to_string()),
                created_at: node.created_at_unix_ms,
                kind: EngineLineageRenderRowKind::Sequence,
                length,
                circular,
                pool_size,
                arrangement_id: None,
                arrangement_mode: None,
                lane_count: 0,
                analysis_kind: None,
                analysis_status: None,
                analysis_reference_seq_id: None,
                analysis_read_count: None,
                analysis_trace_count: None,
                analysis_target_count: None,
                analysis_variant_count: None,
                macro_instance_id: None,
                macro_op_count: 0,
            }
        })
        .collect();
    rows.sort_by(|a, b| {
        a.created_at
            .cmp(&b.created_at)
            .then(a.node_id.cmp(&b.node_id))
    });

    let raw_lineage_edges: Vec<(String, String, String)> = state
        .lineage
        .edges
        .iter()
        .filter_map(|edge| {
            let from = state.lineage.nodes.get(&edge.from_node_id)?.node_id.clone();
            let to = state.lineage.nodes.get(&edge.to_node_id)?.node_id.clone();
            Some((from, to, edge.op_id.clone()))
        })
        .collect();

    let (mut projected_rows, mut projected_edges) =
        project_lineage_operation_hubs(&rows, &raw_lineage_edges, &hub_op_ids);
    for op_id in &hub_op_ids {
        op_label_by_id.insert(op_id.clone(), "Gibson cloning".to_string());
        op_label_by_id.insert(format!("{op_id}::hub_in"), "Gibson cloning".to_string());
        op_label_by_id.insert(format!("{op_id}::hub_out"), "Gibson cloning".to_string());
    }

    for (id, arrangement) in &state.container_state.arrangements {
        let mut source_node_ids: Vec<String> = vec![];
        let mut seen: HashSet<String> = HashSet::new();
        for container_id in &arrangement.lane_container_ids {
            if let Some(container) = state.container_state.containers.get(container_id) {
                for seq_id in &container.members {
                    if let Some(node_id) = state.lineage.seq_to_node.get(seq_id)
                        && seen.insert(node_id.clone())
                    {
                        source_node_ids.push(node_id.clone());
                    }
                }
            }
        }
        let arrangement_node_id = format!("arr:{id}");
        let created_by_op = arrangement
            .created_by_op
            .clone()
            .unwrap_or_else(|| "CreateArrangementSerial".to_string());
        projected_rows.push(EngineLineageRenderRow {
            node_id: arrangement_node_id.clone(),
            seq_id: id.clone(),
            display_name: arrangement.name.clone().unwrap_or_else(|| id.clone()),
            created_by_op: created_by_op.clone(),
            created_at: arrangement.created_at_unix_ms,
            kind: EngineLineageRenderRowKind::Arrangement,
            length: 0,
            circular: false,
            pool_size: 0,
            arrangement_id: Some(id.clone()),
            arrangement_mode: Some(format!("{:?}", arrangement.mode)),
            lane_count: arrangement.lane_container_ids.len(),
            analysis_kind: None,
            analysis_status: None,
            analysis_reference_seq_id: None,
            analysis_read_count: None,
            analysis_trace_count: None,
            analysis_target_count: None,
            analysis_variant_count: None,
            macro_instance_id: None,
            macro_op_count: 0,
        });
        for from_node_id in source_node_ids {
            projected_edges.push((
                from_node_id,
                arrangement_node_id.clone(),
                created_by_op.clone(),
            ));
        }
    }

    for instance in &state.lineage.macro_instances {
        let macro_node_id = format!("macro:{}", instance.macro_instance_id);
        projected_rows.push(EngineLineageRenderRow {
            node_id: macro_node_id.clone(),
            seq_id: instance.macro_instance_id.clone(),
            display_name: instance
                .routine_title
                .clone()
                .or_else(|| instance.routine_id.clone())
                .or_else(|| instance.template_name.clone())
                .unwrap_or_else(|| "Macro".to_string()),
            created_by_op: instance
                .expanded_op_ids
                .first()
                .cloned()
                .unwrap_or_else(|| "-".to_string()),
            created_at: instance.created_at_unix_ms,
            kind: EngineLineageRenderRowKind::Macro,
            length: 0,
            circular: false,
            pool_size: 0,
            arrangement_id: None,
            arrangement_mode: None,
            lane_count: 0,
            analysis_kind: None,
            analysis_status: None,
            analysis_reference_seq_id: None,
            analysis_read_count: None,
            analysis_trace_count: None,
            analysis_target_count: None,
            analysis_variant_count: None,
            macro_instance_id: Some(instance.macro_instance_id.clone()),
            macro_op_count: instance.expanded_op_ids.len(),
        });
        for binding in &instance.bound_inputs {
            let edge_label_id = format!(
                "macro:{}:in:{}",
                instance.macro_instance_id, binding.port_id
            );
            op_label_by_id
                .entry(edge_label_id.clone())
                .or_insert_with(|| format!("in:{}", binding.port_id));
            for from_node_id in sequence_nodes_for_binding(state, binding) {
                projected_edges.push((from_node_id, macro_node_id.clone(), edge_label_id.clone()));
            }
        }
        for binding in &instance.bound_outputs {
            let edge_label_id = format!(
                "macro:{}:out:{}",
                instance.macro_instance_id, binding.port_id
            );
            op_label_by_id
                .entry(edge_label_id.clone())
                .or_insert_with(|| format!("out:{}", binding.port_id));
            for to_node_id in sequence_nodes_for_binding(state, binding) {
                projected_edges.push((macro_node_id.clone(), to_node_id, edge_label_id.clone()));
            }
        }
    }

    for report in GentleEngine::sequencing_confirmation_reports_from_state(state) {
        let node_id = format!("analysis:seq_confirm:{}", report.report_id);
        let op_binding = sequencing_confirmation_op_by_report_id
            .get(&report.report_id)
            .cloned();
        let created_by_op = op_binding
            .as_ref()
            .map(|(op_id, _, _)| op_id.clone())
            .unwrap_or_else(|| "-".to_string());
        let edge_op_id = if created_by_op == "-" {
            format!("analysis:seq_confirm:{}", report.report_id)
        } else {
            created_by_op.clone()
        };
        op_label_by_id.entry(edge_op_id.clone()).or_insert_with(|| {
            format!(
                "Sequencing confirmation: expected={}, baseline={}, report_id={}",
                report.expected_seq_id,
                report.baseline_seq_id.as_deref().unwrap_or("-"),
                report.report_id
            )
        });
        let expected_seq_id = report.expected_seq_id.clone();
        let baseline_seq_id = report.baseline_seq_id.clone().or_else(|| {
            op_binding
                .as_ref()
                .and_then(|(_, _, baseline_seq_id)| baseline_seq_id.clone())
        });
        projected_rows.push(EngineLineageRenderRow {
            node_id: node_id.clone(),
            seq_id: expected_seq_id.clone(),
            display_name: report.report_id.clone(),
            created_by_op,
            created_at: report.generated_at_unix_ms,
            kind: EngineLineageRenderRowKind::Analysis,
            length: 0,
            circular: false,
            pool_size: 0,
            arrangement_id: None,
            arrangement_mode: None,
            lane_count: 0,
            analysis_kind: Some("sequencing_confirmation".to_string()),
            analysis_status: Some(report.overall_status.as_str().to_string()),
            analysis_reference_seq_id: baseline_seq_id.clone(),
            analysis_read_count: Some(report.read_seq_ids.len()),
            analysis_trace_count: Some(report.trace_ids.len()),
            analysis_target_count: Some(report.targets.len()),
            analysis_variant_count: Some(report.variants.len()),
            macro_instance_id: None,
            macro_op_count: 0,
        });
        let mut seen_sources: HashSet<String> = HashSet::new();
        for source_seq_id in std::iter::once(expected_seq_id)
            .chain(baseline_seq_id.into_iter())
            .collect::<Vec<_>>()
        {
            let Some(source_node_id) = state.lineage.seq_to_node.get(&source_seq_id) else {
                continue;
            };
            if seen_sources.insert(source_node_id.clone()) {
                projected_edges.push((source_node_id.clone(), node_id.clone(), edge_op_id.clone()));
            }
        }
    }

    projected_rows.sort_by(|a, b| {
        a.created_at
            .cmp(&b.created_at)
            .then(a.node_id.cmp(&b.node_id))
    });
    let layout_by_node = compute_lineage_dag_layout(&projected_rows, &projected_edges);
    let nodes: Vec<LineageSvgNode> = projected_rows
        .iter()
        .enumerate()
        .map(|(fallback_rank, row)| {
            let (layer, rank) = layout_by_node
                .get(&row.node_id)
                .copied()
                .unwrap_or((0, fallback_rank));
            LineageSvgNode {
                node_id: row.node_id.clone(),
                title: lineage_svg_node_title(row),
                subtitle: lineage_svg_node_subtitle(row),
                kind: lineage_svg_node_kind(row),
                x: 120.0 + layer as f32 * 220.0,
                y: 120.0 + rank as f32 * 110.0,
            }
        })
        .collect();
    let edges: Vec<LineageSvgEdge> = projected_edges
        .into_iter()
        .map(|(from_node_id, to_node_id, op_id)| LineageSvgEdge {
            from_node_id,
            to_node_id,
            label: op_label_by_id.get(&op_id).cloned().unwrap_or(op_id),
        })
        .collect();
    (nodes, edges)
}

pub fn export_lineage_svg(state: &ProjectState, operation_log: &[OperationRecord]) -> String {
    let (nodes, edges) = build_lineage_svg_graph(state, operation_log);
    export_projected_lineage_svg("GENtle Lineage (DALG)", &nodes, &edges)
}
