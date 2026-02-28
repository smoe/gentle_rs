//! Lineage graph export and serialization utilities.

use crate::engine::ProjectState;
use std::collections::{HashMap, HashSet};
use svg::Document;
use svg::node::element::{Circle, Line, Polygon, Rectangle, Text};

const W: f32 = 1600.0;
const H: f32 = 900.0;

pub fn export_lineage_svg(state: &ProjectState) -> String {
    let mut rows: Vec<(String, String, String, usize, usize)> = state
        .lineage
        .nodes
        .values()
        .map(|node| {
            let display_name = state
                .sequences
                .get(&node.seq_id)
                .and_then(|dna| dna.name().clone())
                .unwrap_or_else(|| node.seq_id.clone());
            let length = state
                .sequences
                .get(&node.seq_id)
                .map(|dna| dna.len())
                .unwrap_or(0);
            (
                node.node_id.clone(),
                node.seq_id.clone(),
                display_name,
                length,
                node.created_at_unix_ms as usize,
            )
        })
        .collect();
    rows.sort_by(|a, b| a.4.cmp(&b.4).then(a.0.cmp(&b.0)));

    let mut op_created_count: HashMap<String, usize> = HashMap::new();
    for node in state.lineage.nodes.values() {
        if let Some(op_id) = &node.created_by_op {
            *op_created_count.entry(op_id.clone()).or_insert(0) += 1;
        }
    }

    let mut pos_by_node: HashMap<String, (f32, f32)> = HashMap::new();
    for (idx, row) in rows.iter().enumerate() {
        let mut h = std::collections::hash_map::DefaultHasher::new();
        use std::hash::{Hash, Hasher};
        row.0.hash(&mut h);
        let lane = (h.finish() % 5) as f32;
        let x = 110.0 + idx as f32 * 170.0;
        let y = 120.0 + lane * 130.0;
        pos_by_node.insert(row.0.clone(), (x, y));
    }

    struct ArrangementRenderRow {
        node_id: String,
        arrangement_id: String,
        display_name: String,
        created_by_op: String,
        created_at: usize,
        mode: String,
        lane_container_ids: Vec<String>,
        ladders: Vec<String>,
        source_node_ids: Vec<String>,
    }

    let mut arrangement_rows: Vec<ArrangementRenderRow> = state
        .container_state
        .arrangements
        .iter()
        .map(|(id, arrangement)| {
            let mut source_node_ids: Vec<String> = vec![];
            let mut seen: HashSet<String> = HashSet::new();
            for container_id in &arrangement.lane_container_ids {
                if let Some(container) = state.container_state.containers.get(container_id) {
                    for seq_id in &container.members {
                        if let Some(node_id) = state.lineage.seq_to_node.get(seq_id) {
                            if seen.insert(node_id.clone()) {
                                source_node_ids.push(node_id.clone());
                            }
                        }
                    }
                }
            }
            ArrangementRenderRow {
                node_id: format!("arr:{id}"),
                arrangement_id: id.clone(),
                display_name: arrangement.name.clone().unwrap_or_else(|| id.clone()),
                created_by_op: arrangement
                    .created_by_op
                    .clone()
                    .unwrap_or_else(|| "CreateArrangementSerial".to_string()),
                created_at: arrangement.created_at_unix_ms as usize,
                mode: format!("{:?}", arrangement.mode),
                lane_container_ids: arrangement.lane_container_ids.clone(),
                ladders: arrangement.ladders.clone(),
                source_node_ids,
            }
        })
        .collect();
    arrangement_rows.sort_by(|a, b| {
        a.created_at
            .cmp(&b.created_at)
            .then(a.arrangement_id.cmp(&b.arrangement_id))
    });

    let max_sequence_x = pos_by_node.values().map(|(x, _)| *x).fold(110.0, f32::max);
    let mut arrangement_layer_counts: HashMap<usize, usize> = HashMap::new();
    let mut arrangement_positions: HashMap<String, (f32, f32)> = HashMap::new();
    for arrangement in &arrangement_rows {
        let mut source_x_max = 110.0f32;
        let mut source_y_sum = 0.0f32;
        let mut source_y_count = 0usize;
        for source_node_id in &arrangement.source_node_ids {
            if let Some((sx, sy)) = pos_by_node.get(source_node_id) {
                source_x_max = source_x_max.max(*sx);
                source_y_sum += *sy;
                source_y_count += 1;
            }
        }
        let lane = (((source_x_max - 110.0) / 170.0).round() as usize).saturating_add(1);
        let lane_idx = arrangement_layer_counts
            .entry(lane)
            .and_modify(|idx| *idx += 1)
            .or_insert(0usize);
        let x = (source_x_max + 190.0).max(max_sequence_x + 120.0);
        let y = if source_y_count > 0 {
            source_y_sum / source_y_count as f32
        } else {
            120.0 + (*lane_idx as f32) * 120.0
        };
        arrangement_positions.insert(arrangement.node_id.clone(), (x, y));
        pos_by_node.insert(arrangement.node_id.clone(), (x, y));
    }

    struct MacroRenderRow {
        node_id: String,
        macro_instance_id: String,
        display_name: String,
        created_at: usize,
        template_name: Option<String>,
        routine_id: Option<String>,
        status: String,
        status_message: Option<String>,
        input_bindings: Vec<crate::engine::LineageMacroPortBinding>,
        output_bindings: Vec<crate::engine::LineageMacroPortBinding>,
        expanded_op_ids: Vec<String>,
    }

    let mut macro_rows: Vec<MacroRenderRow> = state
        .lineage
        .macro_instances
        .iter()
        .map(|instance| MacroRenderRow {
            node_id: format!("macro:{}", instance.macro_instance_id),
            macro_instance_id: instance.macro_instance_id.clone(),
            display_name: instance
                .routine_title
                .clone()
                .or_else(|| instance.routine_id.clone())
                .or_else(|| instance.template_name.clone())
                .unwrap_or_else(|| "Macro".to_string()),
            created_at: instance.created_at_unix_ms as usize,
            template_name: instance.template_name.clone(),
            routine_id: instance.routine_id.clone(),
            status: format!("{:?}", instance.status).to_ascii_lowercase(),
            status_message: instance.status_message.clone(),
            input_bindings: instance.bound_inputs.clone(),
            output_bindings: instance.bound_outputs.clone(),
            expanded_op_ids: instance.expanded_op_ids.clone(),
        })
        .collect();
    macro_rows.sort_by(|a, b| {
        a.created_at
            .cmp(&b.created_at)
            .then(a.macro_instance_id.cmp(&b.macro_instance_id))
    });
    for (idx, macro_row) in macro_rows.iter().enumerate() {
        let mut h = std::collections::hash_map::DefaultHasher::new();
        use std::hash::{Hash, Hasher};
        macro_row.node_id.hash(&mut h);
        let lane = (h.finish() % 5) as f32;
        let x = 140.0 + idx as f32 * 170.0;
        let y = 160.0 + lane * 130.0;
        pos_by_node.insert(macro_row.node_id.clone(), (x, y));
    }

    let mut doc = Document::new()
        .set("viewBox", (0, 0, W, H))
        .set("width", W)
        .set("height", H)
        .set("style", "background:#ffffff");

    doc = doc.add(
        Text::new("GENtle Lineage (DALG)")
            .set("x", 24)
            .set("y", 34)
            .set("font-family", "Helvetica, Arial, sans-serif")
            .set("font-size", 24)
            .set("fill", "#202020"),
    );

    for edge in &state.lineage.edges {
        let Some((fx, fy)) = state
            .lineage
            .nodes
            .get(&edge.from_node_id)
            .and_then(|n| pos_by_node.get(&n.node_id))
            .cloned()
        else {
            continue;
        };
        let Some((tx, ty)) = state
            .lineage
            .nodes
            .get(&edge.to_node_id)
            .and_then(|n| pos_by_node.get(&n.node_id))
            .cloned()
        else {
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
                Text::new(edge.op_id.clone())
                    .set("x", mx)
                    .set("y", my)
                    .set("text-anchor", "middle")
                    .set("dominant-baseline", "middle")
                    .set("font-family", "Helvetica, Arial, sans-serif")
                    .set("font-size", 10)
                    .set("fill", "#222222"),
            );
    }

    for arrangement in &arrangement_rows {
        let Some((tx, ty)) = arrangement_positions.get(&arrangement.node_id).copied() else {
            continue;
        };
        for from_node_id in &arrangement.source_node_ids {
            let Some((fx, fy)) = pos_by_node.get(from_node_id).copied() else {
                continue;
            };
            doc = doc.add(
                Line::new()
                    .set("x1", fx)
                    .set("y1", fy)
                    .set("x2", tx)
                    .set("y2", ty)
                    .set("stroke", "#5b6f65")
                    .set("stroke-width", 1.1)
                    .set("stroke-dasharray", "4 3"),
            );
            let mx = (fx + tx) * 0.5;
            let my = (fy + ty) * 0.5 - 7.0;
            doc = doc
                .add(
                    Rectangle::new()
                        .set("x", mx - 58.0)
                        .set("y", my - 12.0)
                        .set("width", 116)
                        .set("height", 16)
                        .set("fill", "#eef6f1")
                        .set("stroke", "#d4e6dd")
                        .set("rx", 2),
                )
                .add(
                    Text::new(arrangement.created_by_op.clone())
                        .set("x", mx)
                        .set("y", my)
                        .set("text-anchor", "middle")
                        .set("dominant-baseline", "middle")
                        .set("font-family", "Helvetica, Arial, sans-serif")
                        .set("font-size", 9)
                        .set("fill", "#1f2b25"),
                );
        }
    }

    let sequence_nodes_for_binding =
        |binding: &crate::engine::LineageMacroPortBinding| -> Vec<String> {
            let mut node_ids: Vec<String> = vec![];
            let mut seen: HashSet<String> = HashSet::new();
            match binding.kind.trim().to_ascii_lowercase().as_str() {
                "sequence" => {
                    for seq_id in &binding.values {
                        if let Some(node_id) = state.lineage.seq_to_node.get(seq_id) {
                            if seen.insert(node_id.clone()) {
                                node_ids.push(node_id.clone());
                            }
                        }
                    }
                }
                "container" => {
                    for container_id in &binding.values {
                        let Some(container) = state.container_state.containers.get(container_id)
                        else {
                            continue;
                        };
                        for seq_id in &container.members {
                            if let Some(node_id) = state.lineage.seq_to_node.get(seq_id) {
                                if seen.insert(node_id.clone()) {
                                    node_ids.push(node_id.clone());
                                }
                            }
                        }
                    }
                }
                _ => {}
            }
            node_ids
        };

    for macro_row in &macro_rows {
        let Some((mx, my)) = pos_by_node.get(&macro_row.node_id).copied() else {
            continue;
        };
        for binding in &macro_row.input_bindings {
            let source_nodes = sequence_nodes_for_binding(binding);
            for source_node_id in source_nodes {
                let Some((sx, sy)) = pos_by_node.get(&source_node_id).copied() else {
                    continue;
                };
                doc = doc.add(
                    Line::new()
                        .set("x1", sx)
                        .set("y1", sy)
                        .set("x2", mx)
                        .set("y2", my)
                        .set("stroke", "#6f6f7b")
                        .set("stroke-width", 1.0)
                        .set("stroke-dasharray", "5 3"),
                );
                let lx = (sx + mx) * 0.5;
                let ly = (sy + my) * 0.5 - 6.0;
                doc = doc
                    .add(
                        Rectangle::new()
                            .set("x", lx - 44.0)
                            .set("y", ly - 10.0)
                            .set("width", 88)
                            .set("height", 14)
                            .set("fill", "#f3f3f7")
                            .set("stroke", "#e2e2ea")
                            .set("rx", 2),
                    )
                    .add(
                        Text::new(format!("in:{}", binding.port_id))
                            .set("x", lx)
                            .set("y", ly)
                            .set("text-anchor", "middle")
                            .set("dominant-baseline", "middle")
                            .set("font-family", "Helvetica, Arial, sans-serif")
                            .set("font-size", 8.5)
                            .set("fill", "#30303a"),
                    );
            }
        }
        for binding in &macro_row.output_bindings {
            let target_nodes = sequence_nodes_for_binding(binding);
            for target_node_id in target_nodes {
                let Some((tx, ty)) = pos_by_node.get(&target_node_id).copied() else {
                    continue;
                };
                doc = doc.add(
                    Line::new()
                        .set("x1", mx)
                        .set("y1", my)
                        .set("x2", tx)
                        .set("y2", ty)
                        .set("stroke", "#6f6f7b")
                        .set("stroke-width", 1.0)
                        .set("stroke-dasharray", "5 3"),
                );
                let lx = (tx + mx) * 0.5;
                let ly = (ty + my) * 0.5 - 6.0;
                doc = doc
                    .add(
                        Rectangle::new()
                            .set("x", lx - 44.0)
                            .set("y", ly - 10.0)
                            .set("width", 88)
                            .set("height", 14)
                            .set("fill", "#f3f3f7")
                            .set("stroke", "#e2e2ea")
                            .set("rx", 2),
                    )
                    .add(
                        Text::new(format!("out:{}", binding.port_id))
                            .set("x", lx)
                            .set("y", ly)
                            .set("text-anchor", "middle")
                            .set("dominant-baseline", "middle")
                            .set("font-family", "Helvetica, Arial, sans-serif")
                            .set("font-size", 8.5)
                            .set("fill", "#30303a"),
                    );
            }
        }
    }

    for row in &rows {
        let Some((x, y)) = pos_by_node.get(&row.0).cloned() else {
            continue;
        };
        let pool_size = state
            .lineage
            .nodes
            .get(&row.0)
            .and_then(|n| n.created_by_op.clone())
            .and_then(|op| op_created_count.get(&op).cloned())
            .unwrap_or(1);
        if pool_size > 1 {
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
            doc = doc
                .add(
                    Polygon::new()
                        .set("points", points)
                        .set("fill", "#b47846")
                        .set("stroke", "#a05f2b")
                        .set("stroke-width", 1),
                )
                .add(
                    Text::new(format!("n={pool_size}"))
                        .set("x", x + 22.0)
                        .set("y", y - 12.0)
                        .set("font-family", "Helvetica, Arial, sans-serif")
                        .set("font-size", 10)
                        .set("fill", "#5c4300"),
                );
        } else {
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

        doc = doc
            .add(
                Text::new(row.2.clone())
                    .set("x", x + 24.0)
                    .set("y", y - 2.0)
                    .set("font-family", "Helvetica, Arial, sans-serif")
                    .set("font-size", 12)
                    .set("fill", "#101010"),
            )
            .add(
                Text::new(format!("{} ({} bp)", row.1, row.3))
                    .set("x", x + 24.0)
                    .set("y", y + 12.0)
                    .set("font-family", "Helvetica, Arial, sans-serif")
                    .set("font-size", 10)
                    .set("fill", "#222222"),
            );
    }

    for arrangement in &arrangement_rows {
        let Some((x, y)) = arrangement_positions.get(&arrangement.node_id).copied() else {
            continue;
        };
        let ladders = if arrangement.ladders.is_empty() {
            "auto".to_string()
        } else {
            arrangement.ladders.join(" + ")
        };
        doc = doc
            .add(
                Rectangle::new()
                    .set("x", x - 36.0)
                    .set("y", y - 14.0)
                    .set("width", 72)
                    .set("height", 28)
                    .set("rx", 5)
                    .set("fill", "#6c9a7a")
                    .set("stroke", "#537563")
                    .set("stroke-width", 1),
            )
            .add(
                Text::new(arrangement.arrangement_id.clone())
                    .set("x", x)
                    .set("y", y)
                    .set("text-anchor", "middle")
                    .set("dominant-baseline", "middle")
                    .set("font-family", "Helvetica, Arial, sans-serif")
                    .set("font-size", 10)
                    .set("fill", "#ffffff"),
            )
            .add(
                Text::new(arrangement.display_name.clone())
                    .set("x", x + 42.0)
                    .set("y", y - 2.0)
                    .set("font-family", "Helvetica, Arial, sans-serif")
                    .set("font-size", 12)
                    .set("fill", "#101010"),
            )
            .add(
                Text::new(format!(
                    "{} | lanes={} | ladders={}",
                    arrangement.mode.to_lowercase(),
                    arrangement.lane_container_ids.len(),
                    ladders
                ))
                .set("x", x + 42.0)
                .set("y", y + 12.0)
                .set("font-family", "Helvetica, Arial, sans-serif")
                .set("font-size", 10)
                .set("fill", "#1f2b25"),
            );
    }

    for macro_row in &macro_rows {
        let Some((x, y)) = pos_by_node.get(&macro_row.node_id).copied() else {
            continue;
        };
        let op_count = macro_row.expanded_op_ids.len();
        doc = doc
            .add(
                Rectangle::new()
                    .set("x", x - 40.0)
                    .set("y", y - 15.0)
                    .set("width", 80)
                    .set("height", 30)
                    .set("rx", 4)
                    .set("fill", "#62626c")
                    .set("stroke", "#4f4f58")
                    .set("stroke-width", 1),
            )
            .add(
                Text::new(macro_row.macro_instance_id.clone())
                    .set("x", x)
                    .set("y", y)
                    .set("text-anchor", "middle")
                    .set("dominant-baseline", "middle")
                    .set("font-family", "Helvetica, Arial, sans-serif")
                    .set("font-size", 10)
                    .set("fill", "#ffffff"),
            )
            .add(
                Text::new(macro_row.display_name.clone())
                    .set("x", x + 46.0)
                    .set("y", y - 2.0)
                    .set("font-family", "Helvetica, Arial, sans-serif")
                    .set("font-size", 12)
                    .set("fill", "#101010"),
            )
            .add(
                Text::new(format!(
                    "status={} | template={} | routine={} | ops={}",
                    macro_row.status,
                    macro_row.template_name.as_deref().unwrap_or("-"),
                    macro_row.routine_id.as_deref().unwrap_or("-"),
                    op_count
                ))
                .set("x", x + 46.0)
                .set("y", y + 12.0)
                .set("font-family", "Helvetica, Arial, sans-serif")
                .set("font-size", 10)
                .set("fill", "#30303a"),
            );
        if let Some(status_message) = macro_row.status_message.as_deref() {
            let compact = if status_message.chars().count() > 72 {
                let mut truncated = status_message.chars().take(72).collect::<String>();
                truncated.push_str("...");
                truncated
            } else {
                status_message.to_string()
            };
            doc = doc.add(
                Text::new(format!("note={compact}"))
                    .set("x", x + 46.0)
                    .set("y", y + 24.0)
                    .set("font-family", "Helvetica, Arial, sans-serif")
                    .set("font-size", 9)
                    .set("fill", "#7a1c1c"),
            );
        }
    }

    doc.to_string()
}
