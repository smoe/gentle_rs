use crate::engine::ProjectState;
use std::collections::HashMap;
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

    doc.to_string()
}
