//! Virtual pool gel model and rendering primitives.

use gentle_protocol::{
    GelBufferModel, GelRunConditions, GelTopologyForm, LadderCatalog, default_dna_ladders,
};
use std::collections::BTreeSet;
use std::sync::LazyLock;
use svg::Document;
use svg::node::element::{Line, Rectangle, Text};

const SVG_WIDTH: f32 = 1320.0;
const SVG_HEIGHT: f32 = 760.0;
const GEL_LEFT: f32 = 90.0;
const GEL_RIGHT: f32 = SVG_WIDTH - 420.0;
const GEL_TOP: f32 = 90.0;
const GEL_BOTTOM: f32 = SVG_HEIGHT - 110.0;
const DETAIL_PANEL_LEFT: f32 = GEL_RIGHT + 90.0;
const DETAIL_PANEL_TOP: f32 = GEL_TOP + 8.0;
const DETAIL_PANEL_WIDTH: f32 = SVG_WIDTH - DETAIL_PANEL_LEFT - 36.0;

static DNA_LADDERS: LazyLock<LadderCatalog> = LazyLock::new(default_dna_ladders);

#[derive(Clone, Debug)]
pub struct GelSampleMember {
    pub seq_id: String,
    pub bp: usize,
    pub topology_form: GelTopologyForm,
}

#[derive(Clone, Debug)]
pub struct PoolGelBand {
    pub bp: usize,
    pub min_bp: usize,
    pub apparent_bp: usize,
    pub intensity: f32,
    pub count: usize,
    pub estimated_mass_units: f32,
    pub topology_label: String,
    pub labels: Vec<String>,
}

#[derive(Clone, Debug)]
pub struct PoolGelLane {
    pub name: String,
    pub role_label: Option<String>,
    pub is_ladder: bool,
    pub bands: Vec<PoolGelBand>,
}

#[derive(Clone, Debug)]
pub struct PoolGelLayout {
    pub lanes: Vec<PoolGelLane>,
    pub selected_ladders: Vec<String>,
    pub sample_count: usize,
    pub pool_member_count: usize,
    pub range_min_bp: usize,
    pub range_max_bp: usize,
    pub conditions: GelRunConditions,
}

#[derive(Clone, Debug)]
pub struct GelSampleInput {
    pub name: String,
    pub role_label: Option<String>,
    pub members: Vec<GelSampleMember>,
}

impl PoolGelLayout {
    pub fn y_for_bp(&self, bp: usize, top: f32, bottom: f32) -> f32 {
        let min_bp = self.range_min_bp.max(1) as f64;
        let max_bp = self.range_max_bp.max(self.range_min_bp + 1) as f64;
        let bp = bp.clamp(self.range_min_bp.max(1), self.range_max_bp.max(2)) as f64;
        let log_min = min_bp.log10();
        let log_max = max_bp.log10();
        let denom = (log_max - log_min).max(1e-6);
        let f = ((log_max - bp.log10()) / denom) as f32;
        let conditions = self.conditions.normalized();
        let mut exponent = 1.0 + (conditions.agarose_percent - 1.0) * 0.18;
        if matches!(conditions.buffer_model, GelBufferModel::Tbe) {
            exponent += 0.05;
        }
        let curved = f.clamp(0.0, 1.0).powf(exponent.clamp(0.72, 1.55));
        top + curved * (bottom - top)
    }
}

fn normalize_ladder_name(name: &str) -> String {
    name.trim().to_ascii_lowercase()
}

fn resolve_ladder_names(requested: &[String], min_bp: usize, max_bp: usize) -> Vec<String> {
    let names = DNA_LADDERS.names_sorted();
    if names.is_empty() {
        return vec![];
    }

    if !requested.is_empty() {
        let mut picked: Vec<String> = vec![];
        for requested_name in requested.iter().take(2) {
            let needle = normalize_ladder_name(requested_name);
            if needle.is_empty() {
                continue;
            }
            if let Some(found) = names
                .iter()
                .find(|name| normalize_ladder_name(name) == needle)
                .cloned()
            {
                if !picked.iter().any(|n| n == &found) {
                    picked.push(found);
                }
                continue;
            }
            if let Some(found) = names
                .iter()
                .find(|name| normalize_ladder_name(name).contains(&needle))
                .cloned()
            {
                if !picked.iter().any(|n| n == &found) {
                    picked.push(found);
                }
            }
        }
        if !picked.is_empty() {
            return picked;
        }
    }

    DNA_LADDERS.choose_for_range(min_bp, max_bp, 2)
}

fn apparent_bp_for_member(member: &GelSampleMember, conditions: &GelRunConditions) -> usize {
    let normalized = conditions.normalized();
    let mut apparent_bp = member.bp as f32;
    if normalized.topology_aware {
        let log_span = (member.bp.max(10) as f32).log10().clamp(2.0, 4.2);
        let circular_baseline = 0.80 + (log_span - 2.0) * 0.04;
        let agarose_adjust = 1.0 - (normalized.agarose_percent - 1.0) * 0.05;
        let factor = match member.topology_form {
            GelTopologyForm::Linear => 1.0,
            GelTopologyForm::Circular => (circular_baseline * agarose_adjust).clamp(0.72, 0.92),
            GelTopologyForm::Supercoiled => {
                (circular_baseline * 0.90 * agarose_adjust).clamp(0.62, 0.86)
            }
            GelTopologyForm::RelaxedCircular => {
                (1.08 + (normalized.agarose_percent - 1.0) * 0.04).clamp(1.02, 1.18)
            }
            GelTopologyForm::NickedCircular => {
                (1.18 + (normalized.agarose_percent - 1.0) * 0.06).clamp(1.08, 1.32)
            }
        };
        apparent_bp *= factor;
    }
    apparent_bp.round().max(1.0) as usize
}

fn estimate_member_mass_units(member: &GelSampleMember) -> f32 {
    member.bp.max(1) as f32
}

fn co_migration_log_tolerance(conditions: &GelRunConditions) -> f32 {
    let normalized = conditions.normalized();
    let mut tolerance = 0.022 - (normalized.agarose_percent - 1.0) * 0.004;
    if matches!(normalized.buffer_model, GelBufferModel::Tbe) {
        tolerance -= 0.002;
    }
    tolerance.clamp(0.012, 0.028)
}

fn normalize_sample_lane_intensities(lanes: &mut [PoolGelLane]) {
    let max_mass = lanes
        .iter()
        .filter(|lane| !lane.is_ladder)
        .flat_map(|lane| lane.bands.iter())
        .map(|band| band.estimated_mass_units)
        .fold(0.0_f32, f32::max)
        .max(1.0);
    for lane in lanes.iter_mut().filter(|lane| !lane.is_ladder) {
        for band in &mut lane.bands {
            let scaled = (band.estimated_mass_units / max_mass).sqrt();
            band.intensity = (0.28 + 0.72 * scaled).clamp(0.24, 1.0);
        }
    }
}

fn sample_bands(members: &[GelSampleMember], conditions: &GelRunConditions) -> Vec<PoolGelBand> {
    let mut projected = members
        .iter()
        .map(|member| (member, apparent_bp_for_member(member, conditions)))
        .collect::<Vec<_>>();
    projected.sort_by(|a, b| {
        b.1.cmp(&a.1)
            .then(b.0.bp.cmp(&a.0.bp))
            .then(a.0.seq_id.cmp(&b.0.seq_id))
    });

    let tolerance = co_migration_log_tolerance(conditions);
    let mut groups: Vec<Vec<(&GelSampleMember, usize)>> = vec![];
    for (member, apparent_bp) in projected {
        let apparent_bp_f = apparent_bp.max(1) as f32;
        let apparent_log = apparent_bp_f.log10();
        let can_merge = groups.last().is_some_and(|group| {
            let representative_log = group
                .iter()
                .map(|(_, bp)| (*bp).max(1) as f32)
                .map(f32::log10)
                .sum::<f32>()
                / group.len() as f32;
            (representative_log - apparent_log).abs() <= tolerance
        });
        if can_merge {
            groups.last_mut().unwrap().push((member, apparent_bp));
        } else {
            groups.push(vec![(member, apparent_bp)]);
        }
    }

    groups
        .into_iter()
        .map(|grouped| {
            let count = grouped.len();
            let bp = grouped
                .iter()
                .map(|(member, _)| member.bp)
                .max()
                .unwrap_or(1);
            let min_bp = grouped
                .iter()
                .map(|(member, _)| member.bp)
                .min()
                .unwrap_or(bp);
            let estimated_mass_units = grouped
                .iter()
                .map(|(member, _)| estimate_member_mass_units(member))
                .sum::<f32>();
            let apparent_bp = (grouped
                .iter()
                .map(|(member, apparent_bp)| {
                    *apparent_bp as f32 * estimate_member_mass_units(member)
                })
                .sum::<f32>()
                / estimated_mass_units.max(1.0))
            .round()
            .max(1.0) as usize;
            let topology_label = {
                let mut forms = grouped.iter().map(|(member, _)| member.topology_form);
                let first = forms.next().unwrap_or(GelTopologyForm::Linear);
                if forms.all(|form| form == first) {
                    first.display_label().to_string()
                } else {
                    "mixed".to_string()
                }
            };
            let mut labels = grouped
                .iter()
                .map(|(member, _)| {
                    if member.topology_form.is_circular() {
                        format!(
                            "{} ({} bp, {})",
                            member.seq_id,
                            member.bp,
                            member.topology_form.display_label()
                        )
                    } else {
                        format!("{} ({} bp)", member.seq_id, member.bp)
                    }
                })
                .collect::<Vec<_>>();
            labels.sort();
            PoolGelBand {
                bp,
                min_bp,
                apparent_bp,
                intensity: 0.5,
                count,
                estimated_mass_units,
                topology_label,
                labels,
            }
        })
        .collect::<Vec<_>>()
}

fn ladder_lane_names_for_display(selected_ladders: &[String]) -> Vec<String> {
    match selected_ladders {
        [] => vec![],
        [name] => vec![name.clone(), name.clone()],
        [left, right] => vec![left.clone(), right.clone()],
        more => {
            let mut names = vec![more[0].clone()];
            names.extend(more[1..more.len() - 1].iter().cloned());
            names.push(more[more.len() - 1].clone());
            names
        }
    }
}

pub fn build_serial_gel_layout(
    samples: &[GelSampleInput],
    requested_ladders: &[String],
    conditions: Option<&GelRunConditions>,
) -> Result<PoolGelLayout, String> {
    if samples.is_empty() {
        return Err("Serial gel needs at least one sample lane".to_string());
    }
    let conditions = conditions.cloned().unwrap_or_default().normalized();
    let mut normalized_samples: Vec<GelSampleInput> = vec![];
    let mut all_members: Vec<GelSampleMember> = vec![];
    for (sample_idx, sample) in samples.iter().enumerate() {
        let mut valid_members = sample
            .members
            .iter()
            .filter(|member| member.bp > 0)
            .cloned()
            .collect::<Vec<_>>();
        if valid_members.is_empty() {
            return Err(format!(
                "Sample lane {} has no sequence lengths > 0 bp",
                sample_idx + 1
            ));
        }
        valid_members.sort_by(|a, b| b.bp.cmp(&a.bp).then(a.seq_id.cmp(&b.seq_id)));
        all_members.extend(valid_members.clone());
        let lane_name = if sample.name.trim().is_empty() {
            format!("Sample {} (n={})", sample_idx + 1, valid_members.len())
        } else {
            sample.name.clone()
        };
        normalized_samples.push(GelSampleInput {
            name: lane_name,
            role_label: sample.role_label.clone(),
            members: valid_members,
        });
    }

    let pool_min = all_members
        .iter()
        .map(|member| member.bp)
        .min()
        .unwrap_or(1);
    let pool_max = all_members
        .iter()
        .map(|member| member.bp)
        .max()
        .unwrap_or(pool_min);
    let selected_ladders = resolve_ladder_names(requested_ladders, pool_min, pool_max);
    if selected_ladders.is_empty() {
        return Err("No DNA ladders available for pool-gel rendering".to_string());
    }

    let mut lanes: Vec<PoolGelLane> = vec![];
    let mut all_band_bps: Vec<usize> = all_members
        .iter()
        .map(|member| apparent_bp_for_member(member, &conditions))
        .collect();

    let display_ladders = ladder_lane_names_for_display(&selected_ladders);
    let left_ladders = display_ladders
        .iter()
        .take(display_ladders.len().div_ceil(2))
        .cloned()
        .collect::<Vec<_>>();
    let right_ladders = display_ladders
        .iter()
        .skip(display_ladders.len().div_ceil(2))
        .cloned()
        .collect::<Vec<_>>();

    for ladder_name in &left_ladders {
        let Some(ladder) = DNA_LADDERS.get(ladder_name) else {
            continue;
        };
        let mut parts = ladder.bands().clone();
        parts.sort_by(|a, b| b.length_bp().total_cmp(&a.length_bp()));
        let max_strength = parts
            .iter()
            .filter_map(|p| p.relative_strength)
            .fold(0.0_f64, f64::max)
            .max(1.0);
        let bands = parts
            .into_iter()
            .map(|p| {
                let bp = p.length_bp().round().max(1.0) as usize;
                all_band_bps.push(bp);
                let raw = p.relative_strength.unwrap_or(1.0).max(0.1);
                let intensity = (raw / max_strength).clamp(0.18, 1.0) as f32;
                PoolGelBand {
                    bp,
                    min_bp: bp,
                    apparent_bp: bp,
                    intensity,
                    count: 1,
                    estimated_mass_units: bp as f32,
                    topology_label: "ladder".to_string(),
                    labels: vec![format!("{bp} bp")],
                }
            })
            .collect::<Vec<_>>();
        lanes.push(PoolGelLane {
            name: ladder.name().to_string(),
            role_label: None,
            is_ladder: true,
            bands,
        });
    }

    if lanes.is_empty() {
        return Err("None of the selected DNA ladders could be resolved".to_string());
    }

    for sample in &normalized_samples {
        lanes.push(PoolGelLane {
            name: sample.name.clone(),
            role_label: sample.role_label.clone(),
            is_ladder: false,
            bands: sample_bands(&sample.members, &conditions),
        });
    }
    normalize_sample_lane_intensities(&mut lanes);

    for ladder_name in &right_ladders {
        let Some(ladder) = DNA_LADDERS.get(ladder_name) else {
            continue;
        };
        let mut parts = ladder.bands().clone();
        parts.sort_by(|a, b| b.length_bp().total_cmp(&a.length_bp()));
        let max_strength = parts
            .iter()
            .filter_map(|p| p.relative_strength)
            .fold(0.0_f64, f64::max)
            .max(1.0);
        let bands = parts
            .into_iter()
            .map(|p| {
                let bp = p.length_bp().round().max(1.0) as usize;
                all_band_bps.push(bp);
                let raw = p.relative_strength.unwrap_or(1.0).max(0.1);
                let intensity = (raw / max_strength).clamp(0.18, 1.0) as f32;
                PoolGelBand {
                    bp,
                    min_bp: bp,
                    apparent_bp: bp,
                    intensity,
                    count: 1,
                    estimated_mass_units: bp as f32,
                    topology_label: "ladder".to_string(),
                    labels: vec![format!("{bp} bp")],
                }
            })
            .collect::<Vec<_>>();
        lanes.push(PoolGelLane {
            name: ladder.name().to_string(),
            role_label: None,
            is_ladder: true,
            bands,
        });
    }

    let min_band = all_band_bps.iter().copied().min().unwrap_or(pool_min);
    let max_band = all_band_bps.iter().copied().max().unwrap_or(pool_max);
    let range_min_bp = ((min_band as f64) * 0.72).floor().max(1.0) as usize;
    let mut range_max_bp = ((max_band as f64) * 1.30).ceil().max(2.0) as usize;
    if range_max_bp <= range_min_bp {
        range_max_bp = range_min_bp + 1;
    }

    Ok(PoolGelLayout {
        lanes,
        selected_ladders,
        sample_count: normalized_samples.len(),
        pool_member_count: all_members.len(),
        range_min_bp,
        range_max_bp,
        conditions,
    })
}

pub fn build_pool_gel_layout(
    pool_members: &[GelSampleMember],
    requested_ladders: &[String],
    conditions: Option<&GelRunConditions>,
) -> Result<PoolGelLayout, String> {
    let members = pool_members
        .iter()
        .filter(|member| member.bp > 0)
        .cloned()
        .collect::<Vec<_>>();
    if members.is_empty() {
        return Err("Pool gel needs sequence lengths > 0 bp".to_string());
    }
    let sample = GelSampleInput {
        name: format!("Pool (n={})", members.len()),
        role_label: None,
        members,
    };
    build_serial_gel_layout(&[sample], requested_ladders, conditions)
}

fn format_bp_label(bp: usize) -> String {
    if bp >= 1_000 {
        if bp % 1_000 == 0 {
            format!("{} kb", bp / 1_000)
        } else {
            format!("{:.1} kb", bp as f32 / 1_000.0)
        }
    } else {
        format!("{bp} bp")
    }
}

fn merged_band_note_lines(layout: &PoolGelLayout) -> Vec<String> {
    let mut lines = vec![];
    for lane in layout.lanes.iter().filter(|lane| !lane.is_ladder) {
        for band in lane.bands.iter().filter(|band| band.count > 1) {
            let actual_label = if band.min_bp == band.bp {
                format_bp_label(band.bp)
            } else {
                format!(
                    "{}..{}",
                    format_bp_label(band.min_bp),
                    format_bp_label(band.bp)
                )
            };
            lines.push(format!(
                "{}: {} fragments -> one observed {} band from {}",
                lane.name,
                band.count,
                format_bp_label(band.apparent_bp),
                actual_label
            ));
        }
    }
    lines
}

fn normalized_lane_role(lane: &PoolGelLane) -> String {
    let role = lane
        .role_label
        .as_deref()
        .unwrap_or(lane.name.as_str())
        .trim()
        .to_ascii_lowercase();
    if role.starts_with("insert") {
        "insert".to_string()
    } else if role.starts_with("vector") || role.starts_with("destination") {
        "vector".to_string()
    } else if role.starts_with("product") || role.starts_with("assembled") {
        "product".to_string()
    } else {
        role
    }
}

fn canonical_role_display(role: &str) -> &'static str {
    match role {
        "insert" => "Insert",
        "vector" => "Vector",
        "product" => "Product",
        _ => "Sample",
    }
}

fn lane_hint_prefix(lane: &PoolGelLane, role: &str, plural: bool) -> String {
    let canonical = canonical_role_display(role);
    if plural {
        return format!("{canonical} lanes");
    }
    let lane_name = lane.name.trim();
    let normalized_name = lane_name.to_ascii_lowercase();
    let normalized_role = role.to_ascii_lowercase();
    if lane_name.is_empty() || normalized_name.starts_with(&normalized_role) {
        format!("{canonical} lane")
    } else {
        format!("{canonical} lane ({lane_name})")
    }
}

fn lane_role_badge(role: &str) -> Option<(&'static str, &'static str, &'static str)> {
    match role {
        "insert" => Some(("INSERT", "#d1fae5", "#065f46")),
        "vector" => Some(("VECTOR", "#dbeafe", "#1d4ed8")),
        "product" => Some(("PRODUCT", "#fee2e2", "#b45309")),
        _ => None,
    }
}

fn singleton_actual_bp(lane: &PoolGelLane) -> Option<usize> {
    if lane.is_ladder || lane.bands.len() != 1 {
        return None;
    }
    let band = lane.bands.first()?;
    (band.count == 1).then_some(band.bp)
}

fn comparison_hint_lines(layout: &PoolGelLayout) -> Vec<String> {
    let sample_lanes = layout
        .lanes
        .iter()
        .filter(|lane| !lane.is_ladder)
        .collect::<Vec<_>>();
    let insert_lanes = sample_lanes
        .iter()
        .copied()
        .filter(|lane| normalized_lane_role(lane) == "insert")
        .collect::<Vec<_>>();
    let vector_lane = sample_lanes
        .iter()
        .copied()
        .find(|lane| normalized_lane_role(lane) == "vector");
    let product_lane = sample_lanes
        .iter()
        .copied()
        .find(|lane| normalized_lane_role(lane) == "product");

    let mut lines = vec![];
    if !insert_lanes.is_empty() {
        if insert_lanes.len() == 1 {
            let prefix = lane_hint_prefix(insert_lanes[0], "insert", false);
            if let Some(insert_bp) = singleton_actual_bp(insert_lanes[0]) {
                lines.push(format!(
                    "{prefix}: compare against the fine ladder and confirm the expected {} insert band.",
                    format_bp_label(insert_bp)
                ));
            } else {
                lines.push(
                    format!(
                        "{prefix}: compare against the fine ladder and confirm the small-fragment insert readout."
                    ),
                );
            }
        } else {
            let prefix = lane_hint_prefix(insert_lanes[0], "insert", true);
            let insert_total_bp = insert_lanes
                .iter()
                .filter_map(|lane| singleton_actual_bp(lane))
                .sum::<usize>();
            if insert_total_bp > 0 {
                lines.push(format!(
                    "{prefix}: compare each insert to the fine ladder; combined expected added payload is {}.",
                    format_bp_label(insert_total_bp)
                ));
            } else {
                lines.push(
                    format!(
                        "{prefix}: compare each insert to the fine ladder before reading the product shift."
                    ),
                );
            }
        }
    }

    if let (Some(vector_lane), Some(product_lane)) = (vector_lane, product_lane) {
        let vector_label = if vector_lane.name.trim().is_empty()
            || vector_lane.name.trim().eq_ignore_ascii_case("vector")
        {
            "Vector".to_string()
        } else {
            format!("Vector ({})", vector_lane.name.trim())
        };
        let product_label = if product_lane.name.trim().is_empty()
            || product_lane.name.trim().eq_ignore_ascii_case("product")
        {
            "product".to_string()
        } else {
            format!("product ({})", product_lane.name.trim())
        };
        if let (Some(vector_bp), Some(product_bp)) = (
            singleton_actual_bp(vector_lane),
            singleton_actual_bp(product_lane),
        ) {
            let delta_bp = product_bp.saturating_sub(vector_bp);
            if delta_bp > 0 {
                lines.push(format!(
                    "{vector_label} vs {product_label}: product should run as the larger construct, about {} above the vector backbone.",
                    format_bp_label(delta_bp)
                ));
            } else if delta_bp == 0 {
                lines.push(
                    format!(
                        "{vector_label} vs {product_label}: backbone-sized lanes match closely, so rely on topology label and insert lane confirmation."
                    ),
                );
            }
            let insert_total_bp = insert_lanes
                .iter()
                .filter_map(|lane| singleton_actual_bp(lane))
                .sum::<usize>();
            if insert_total_bp > 0 {
                if delta_bp == insert_total_bp {
                    lines.push(format!(
                        "Consistency check: product-vector delta matches the summed insert payload ({}).",
                        format_bp_label(insert_total_bp)
                    ));
                } else if delta_bp > 0 {
                    lines.push(format!(
                        "Consistency check: product-vector delta is {} while the summed insert payload is {}.",
                        format_bp_label(delta_bp),
                        format_bp_label(insert_total_bp)
                    ));
                }
            }
        } else {
            lines.push(
                "Vector vs product: compare the large-fragment lanes directly; the product should remain the larger assembly."
                    .to_string(),
            );
        }
    }

    lines
}

pub fn export_pool_gel_svg(layout: &PoolGelLayout) -> String {
    let lane_count = layout.lanes.len().max(1);
    let lane_gap = (GEL_RIGHT - GEL_LEFT) / (lane_count as f32 + 1.0);
    let gel_width = GEL_RIGHT - GEL_LEFT;
    let gel_height = GEL_BOTTOM - GEL_TOP;
    let sample_lane_indices = layout
        .lanes
        .iter()
        .enumerate()
        .filter_map(|(idx, lane)| (!lane.is_ladder).then_some(idx))
        .collect::<Vec<_>>();

    let mut doc = Document::new()
        .set("viewBox", (0, 0, SVG_WIDTH, SVG_HEIGHT))
        .set("width", SVG_WIDTH)
        .set("height", SVG_HEIGHT)
        .add(
            Rectangle::new()
                .set("x", 0)
                .set("y", 0)
                .set("width", SVG_WIDTH)
                .set("height", SVG_HEIGHT)
                .set("fill", "#f9fafb"),
        )
        .add(
            Rectangle::new()
                .set("x", GEL_LEFT)
                .set("y", GEL_TOP)
                .set("width", gel_width)
                .set("height", gel_height)
                .set("rx", 10)
                .set("ry", 10)
                .set("fill", "#111315"),
        )
        .add(
            Rectangle::new()
                .set("x", DETAIL_PANEL_LEFT - 18.0)
                .set("y", GEL_TOP)
                .set("width", DETAIL_PANEL_WIDTH + 18.0)
                .set("height", gel_height)
                .set("rx", 10)
                .set("ry", 10)
                .set("fill", "#f3f4f6"),
        );

    let mut tick_bps = BTreeSet::new();
    for lane in layout.lanes.iter().filter(|l| l.is_ladder) {
        for band in &lane.bands {
            tick_bps.insert(band.bp);
        }
    }
    if tick_bps.is_empty() {
        tick_bps.insert(layout.range_min_bp);
        tick_bps.insert(layout.range_max_bp);
    }
    let mut accepted_ticks: Vec<usize> = vec![];
    let mut last_y: Option<f32> = None;
    for bp in tick_bps.iter().rev() {
        let y = layout.y_for_bp(*bp, GEL_TOP, GEL_BOTTOM);
        if last_y.map(|v| (v - y).abs() >= 16.0).unwrap_or(true) {
            accepted_ticks.push(*bp);
            last_y = Some(y);
        }
        if accepted_ticks.len() >= 20 {
            break;
        }
    }
    for bp in accepted_ticks {
        let y = layout.y_for_bp(bp, GEL_TOP, GEL_BOTTOM);
        doc = doc
            .add(
                Line::new()
                    .set("x1", GEL_LEFT)
                    .set("y1", y)
                    .set("x2", GEL_RIGHT)
                    .set("y2", y)
                    .set("stroke", "#2d3238")
                    .set("stroke-width", 1),
            )
            .add(
                Text::new(format_bp_label(bp))
                    .set("x", GEL_RIGHT + 12.0)
                    .set("y", y + 4.0)
                    .set("font-family", "monospace")
                    .set("font-size", 12)
                    .set("fill", "#374151"),
            );
    }

    for (lane_idx, lane) in layout.lanes.iter().enumerate() {
        let x = GEL_LEFT + lane_gap * (lane_idx as f32 + 1.0);
        let lane_fill = if lane.is_ladder { "#1a2028" } else { "#1f252e" };
        doc = doc
            .add(
                Rectangle::new()
                    .set("x", x - 34.0)
                    .set("y", GEL_TOP + 10.0)
                    .set("width", 68.0)
                    .set("height", gel_height - 20.0)
                    .set("rx", 6)
                    .set("ry", 6)
                    .set("fill", lane_fill),
            )
            .add(
                Text::new(lane.name.clone())
                    .set("x", x)
                    .set("y", GEL_BOTTOM + 26.0)
                    .set("text-anchor", "middle")
                    .set("font-family", "monospace")
                    .set("font-size", 13)
                    .set("fill", "#0f172a"),
            );

        if !lane.is_ladder {
            if let Some((badge_text, badge_fill, badge_text_fill)) =
                lane_role_badge(&normalized_lane_role(lane))
            {
                doc = doc
                    .add(
                        Rectangle::new()
                            .set("x", x - 28.0)
                            .set("y", GEL_BOTTOM + 34.0)
                            .set("width", 56.0)
                            .set("height", 16.0)
                            .set("rx", 8)
                            .set("ry", 8)
                            .set("fill", badge_fill),
                    )
                    .add(
                        Text::new(badge_text)
                            .set("x", x)
                            .set("y", GEL_BOTTOM + 45.0)
                            .set("text-anchor", "middle")
                            .set("font-family", "monospace")
                            .set("font-size", 9)
                            .set("font-weight", 700)
                            .set("fill", badge_text_fill),
                    );
            }
        }

        for band in &lane.bands {
            let y = layout.y_for_bp(band.apparent_bp, GEL_TOP + 14.0, GEL_BOTTOM - 14.0);
            let width = if lane.is_ladder {
                32.0 + 18.0 * band.intensity
            } else {
                38.0 + 24.0 * band.intensity
            };
            let height = if lane.is_ladder {
                2.5 + 2.5 * band.intensity
            } else {
                3.0 + 3.0 * band.intensity
            };
            let fill = if lane.is_ladder { "#e5e7eb" } else { "#f59e0b" };
            doc = doc.add(
                Rectangle::new()
                    .set("x", x - width * 0.5)
                    .set("y", y - height * 0.5)
                    .set("width", width)
                    .set("height", height)
                    .set("rx", 2)
                    .set("ry", 2)
                    .set("fill", fill)
                    .set("opacity", (0.42 + 0.58 * band.intensity).clamp(0.35, 1.0)),
            );

            if !lane.is_ladder && band.count > 0 {
                let mut label = if band.min_bp == band.bp {
                    format_bp_label(band.bp)
                } else {
                    format!(
                        "{}..{}",
                        format_bp_label(band.min_bp),
                        format_bp_label(band.bp)
                    )
                };
                if band.apparent_bp != band.bp {
                    label.push_str(&format!(" -> {}", format_bp_label(band.apparent_bp)));
                }
                if band.count > 1 {
                    label.push_str(&format!(" | merged x{}", band.count));
                }
                if !band.topology_label.trim().is_empty() {
                    label.push_str(&format!(" | {}", band.topology_label));
                }
                doc = doc.add(
                    Text::new(label)
                        .set("x", x + 44.0)
                        .set("y", y + 4.0)
                        .set("font-family", "monospace")
                        .set("font-size", 11)
                        .set("fill", "#111827"),
                );
            }
        }
    }

    let ladder_caption = if layout.selected_ladders.is_empty() {
        "auto ladder".to_string()
    } else {
        layout.selected_ladders.join(" + ")
    };
    let title = format!(
        "Serial Gel Preview ({} sample lane(s), {} member(s)) | ladders: {}",
        layout.sample_count, layout.pool_member_count, ladder_caption
    );
    doc = doc
        .add(
            Text::new(title)
                .set("x", GEL_LEFT)
                .set("y", 42.0)
                .set("font-family", "monospace")
                .set("font-size", 16)
                .set("fill", "#0f172a"),
        )
        .add(
            Text::new(format!(
                "range {}..{} | lanes: {} | conditions: {}",
                format_bp_label(layout.range_min_bp),
                format_bp_label(layout.range_max_bp),
                lane_count,
                layout.conditions.describe()
            ))
            .set("x", GEL_LEFT)
            .set("y", 62.0)
            .set("font-family", "monospace")
            .set("font-size", 12)
            .set("fill", "#334155"),
        );

    for sample_idx in sample_lane_indices {
        let x = GEL_LEFT + lane_gap * (sample_idx as f32 + 1.0);
        doc = doc.add(
            Line::new()
                .set("x1", x)
                .set("y1", GEL_TOP)
                .set("x2", x)
                .set("y2", GEL_BOTTOM)
                .set("stroke", "#f59e0b")
                .set("stroke-width", 1.5)
                .set("opacity", 0.25),
        );
    }

    let comparison_hints = comparison_hint_lines(layout);
    let merged_notes = merged_band_note_lines(layout);
    let mut header_y = DETAIL_PANEL_TOP + 24.0;
    if !comparison_hints.is_empty() {
        doc = doc
            .add(
                Text::new("Comparison hints")
                    .set("x", DETAIL_PANEL_LEFT)
                    .set("y", header_y)
                    .set("font-family", "monospace")
                    .set("font-size", 12)
                    .set("font-weight", 700)
                    .set("fill", "#0f172a"),
            )
            .add(
                Text::new(
                    "read these lanes together: insert vs fine ladder, then vector vs product",
                )
                .set("x", DETAIL_PANEL_LEFT + 4.0)
                .set("y", header_y + 16.0)
                .set("font-family", "monospace")
                .set("font-size", 10)
                .set("fill", "#475569"),
            );
        header_y += 32.0;
        for hint in comparison_hints.iter().take(4) {
            doc = doc.add(
                Text::new(hint.clone())
                    .set("x", DETAIL_PANEL_LEFT + 8.0)
                    .set("y", header_y)
                    .set("font-family", "monospace")
                    .set("font-size", 10)
                    .set("fill", "#334155"),
            );
            header_y += 12.0;
        }
        if comparison_hints.len() > 4 {
            doc = doc.add(
                Text::new(format!(
                    "+{} more comparison hint(s)",
                    comparison_hints.len() - 4
                ))
                .set("x", DETAIL_PANEL_LEFT + 8.0)
                .set("y", header_y)
                .set("font-family", "monospace")
                .set("font-size", 10)
                .set("fill", "#64748b"),
            );
            header_y += 12.0;
        }
        header_y += 10.0;
    }

    if !merged_notes.is_empty() {
        doc = doc
            .add(
                Text::new("Merged-band notes")
                    .set("x", DETAIL_PANEL_LEFT)
                    .set("y", header_y)
                    .set("font-family", "monospace")
                    .set("font-size", 12)
                    .set("font-weight", 700)
                    .set("fill", "#0f172a"),
            )
            .add(
                Text::new(
                    "merged xN means nearby fragments co-migrate under current gel conditions",
                )
                .set("x", DETAIL_PANEL_LEFT + 4.0)
                .set("y", header_y + 16.0)
                .set("font-family", "monospace")
                .set("font-size", 10)
                .set("fill", "#475569"),
            )
            .add(
                Text::new("observed = apparent band position | actual = source-size span")
                    .set("x", DETAIL_PANEL_LEFT + 4.0)
                    .set("y", header_y + 30.0)
                    .set("font-family", "monospace")
                    .set("font-size", 10)
                    .set("fill", "#64748b"),
            );
        header_y += 46.0;
        for note in merged_notes.iter().take(4) {
            doc = doc.add(
                Text::new(note.clone())
                    .set("x", DETAIL_PANEL_LEFT + 8.0)
                    .set("y", header_y)
                    .set("font-family", "monospace")
                    .set("font-size", 10)
                    .set("fill", "#334155"),
            );
            header_y += 12.0;
        }
        if merged_notes.len() > 4 {
            doc = doc.add(
                Text::new(format!(
                    "+{} more merged-band note(s)",
                    merged_notes.len() - 4
                ))
                .set("x", DETAIL_PANEL_LEFT + 8.0)
                .set("y", header_y)
                .set("font-family", "monospace")
                .set("font-size", 10)
                .set("fill", "#64748b"),
            );
            header_y += 12.0;
        }
        header_y += 10.0;
    }

    doc = doc
        .add(
            Text::new("Fragment table")
                .set("x", DETAIL_PANEL_LEFT)
                .set("y", header_y + 4.0)
                .set("font-family", "monospace")
                .set("font-size", 14)
                .set("font-weight", 700)
                .set("fill", "#0f172a"),
        )
        .add(
            Text::new("lane | observed | actual | topology | mass")
                .set("x", DETAIL_PANEL_LEFT)
                .set("y", header_y + 24.0)
                .set("font-family", "monospace")
                .set("font-size", 11)
                .set("fill", "#64748b"),
        );

    let mut detail_y = header_y + 44.0;
    for lane in layout.lanes.iter().filter(|lane| !lane.is_ladder) {
        doc = doc.add(
            Text::new(lane.name.clone())
                .set("x", DETAIL_PANEL_LEFT)
                .set("y", detail_y)
                .set("font-family", "monospace")
                .set("font-size", 12)
                .set("font-weight", 700)
                .set("fill", "#111827"),
        );
        detail_y += 16.0;
        for band in &lane.bands {
            let actual_label = if band.min_bp == band.bp {
                format_bp_label(band.bp)
            } else {
                format!(
                    "{}..{}",
                    format_bp_label(band.min_bp),
                    format_bp_label(band.bp)
                )
            };
            let mut row = format!(
                "{} | {} | {} | {:.0} au",
                format_bp_label(band.apparent_bp),
                actual_label,
                band.topology_label,
                band.estimated_mass_units
            );
            if band.count > 1 {
                row.push_str(&format!(" | merged x{}", band.count));
            }
            doc = doc.add(
                Text::new(row)
                    .set("x", DETAIL_PANEL_LEFT + 4.0)
                    .set("y", detail_y)
                    .set("font-family", "monospace")
                    .set("font-size", 11)
                    .set("fill", "#334155"),
            );
            detail_y += 14.0;
            for label in band.labels.iter().take(3) {
                doc = doc.add(
                    Text::new(label.clone())
                        .set("x", DETAIL_PANEL_LEFT + 14.0)
                        .set("y", detail_y)
                        .set("font-family", "monospace")
                        .set("font-size", 10)
                        .set("fill", "#475569"),
                );
                detail_y += 12.0;
            }
            if band.labels.len() > 3 {
                doc = doc.add(
                    Text::new(format!("+{} more", band.labels.len() - 3))
                        .set("x", DETAIL_PANEL_LEFT + 14.0)
                        .set("y", detail_y)
                        .set("font-family", "monospace")
                        .set("font-size", 10)
                        .set("fill", "#64748b"),
                );
                detail_y += 12.0;
            }
            detail_y += 6.0;
            if detail_y > GEL_BOTTOM - 24.0 {
                break;
            }
        }
        detail_y += 8.0;
        if detail_y > GEL_BOTTOM - 24.0 {
            break;
        }
    }

    doc.to_string()
}

#[cfg(test)]
mod tests {
    use super::*;
    #[cfg(feature = "snapshot-tests")]
    use std::fs;

    #[test]
    fn test_build_pool_gel_layout_auto_ladders() {
        let members = vec![
            GelSampleMember {
                seq_id: "frag_a".to_string(),
                bp: 420,
                topology_form: GelTopologyForm::Linear,
            },
            GelSampleMember {
                seq_id: "frag_b".to_string(),
                bp: 950,
                topology_form: GelTopologyForm::Linear,
            },
            GelSampleMember {
                seq_id: "frag_c".to_string(),
                bp: 1210,
                topology_form: GelTopologyForm::Linear,
            },
            GelSampleMember {
                seq_id: "frag_d".to_string(),
                bp: 1210,
                topology_form: GelTopologyForm::Linear,
            },
        ];
        let layout = build_pool_gel_layout(&members, &[], None).unwrap();
        assert!(!layout.selected_ladders.is_empty());
        assert!(layout.lanes.iter().any(|l| l.is_ladder));
        assert!(layout.lanes.iter().any(|l| !l.is_ladder));
        assert_eq!(layout.sample_count, 1);
        assert!(layout.range_max_bp > layout.range_min_bp);
    }

    #[test]
    fn test_export_pool_gel_svg() {
        let members = vec![
            GelSampleMember {
                seq_id: "frag_a".to_string(),
                bp: 350,
                topology_form: GelTopologyForm::Linear,
            },
            GelSampleMember {
                seq_id: "frag_b".to_string(),
                bp: 820,
                topology_form: GelTopologyForm::Linear,
            },
            GelSampleMember {
                seq_id: "frag_c".to_string(),
                bp: 1650,
                topology_form: GelTopologyForm::Linear,
            },
        ];
        let layout = build_pool_gel_layout(&members, &[], None).unwrap();
        let svg = export_pool_gel_svg(&layout);
        assert!(svg.contains("<svg"));
        assert!(svg.contains("Serial Gel Preview"));
        assert!(svg.contains("Fragment table"));
        assert!(svg.contains("conditions:"));
    }

    #[test]
    fn test_build_serial_gel_layout_flanks_samples_with_ladders() {
        let samples = vec![
            GelSampleInput {
                name: "Vector".to_string(),
                role_label: Some("vector".to_string()),
                members: vec![GelSampleMember {
                    seq_id: "vector".to_string(),
                    bp: 4952,
                    topology_form: GelTopologyForm::Circular,
                }],
            },
            GelSampleInput {
                name: "Insert".to_string(),
                role_label: Some("insert_1".to_string()),
                members: vec![GelSampleMember {
                    seq_id: "insert".to_string(),
                    bp: 314,
                    topology_form: GelTopologyForm::Linear,
                }],
            },
            GelSampleInput {
                name: "Product".to_string(),
                role_label: Some("product".to_string()),
                members: vec![GelSampleMember {
                    seq_id: "product".to_string(),
                    bp: 5266,
                    topology_form: GelTopologyForm::Circular,
                }],
            },
        ];
        let layout = build_serial_gel_layout(
            &samples,
            &[
                "NEB 100bp DNA Ladder".to_string(),
                "NEB 1kb DNA Ladder".to_string(),
            ],
            None,
        )
        .expect("layout");
        assert!(layout.lanes.first().is_some_and(|lane| lane.is_ladder));
        assert!(layout.lanes.last().is_some_and(|lane| lane.is_ladder));
        assert_eq!(
            layout
                .lanes
                .iter()
                .filter(|lane| !lane.is_ladder)
                .map(|lane| lane.name.clone())
                .collect::<Vec<_>>(),
            vec![
                "Vector".to_string(),
                "Insert".to_string(),
                "Product".to_string()
            ]
        );
    }

    #[test]
    fn test_export_pool_gel_svg_adds_comparison_hints_for_vector_insert_product() {
        let layout = build_serial_gel_layout(
            &[
                GelSampleInput {
                    name: "Vector".to_string(),
                    role_label: Some("vector".to_string()),
                    members: vec![GelSampleMember {
                        seq_id: "vector".to_string(),
                        bp: 4952,
                        topology_form: GelTopologyForm::Circular,
                    }],
                },
                GelSampleInput {
                    name: "Insert".to_string(),
                    role_label: Some("insert_1".to_string()),
                    members: vec![GelSampleMember {
                        seq_id: "insert".to_string(),
                        bp: 314,
                        topology_form: GelTopologyForm::Linear,
                    }],
                },
                GelSampleInput {
                    name: "Product".to_string(),
                    role_label: Some("product".to_string()),
                    members: vec![GelSampleMember {
                        seq_id: "product".to_string(),
                        bp: 5266,
                        topology_form: GelTopologyForm::Circular,
                    }],
                },
            ],
            &[
                "NEB 100bp DNA Ladder".to_string(),
                "NEB 1kb DNA Ladder".to_string(),
            ],
            None,
        )
        .expect("layout");
        let svg = export_pool_gel_svg(&layout);
        assert!(svg.contains("Comparison hints"));
        assert!(svg.contains("Insert lane: compare against the fine ladder"));
        assert!(svg.contains("Vector vs product: product should run as the larger construct"));
        assert!(
            svg.contains(
                "Consistency check: product-vector delta matches the summed insert payload"
            )
        );
        assert!(svg.contains("314 bp"));
    }

    #[test]
    #[cfg(feature = "snapshot-tests")]
    fn snapshot_pool_gel_svg() {
        let members = vec![
            GelSampleMember {
                seq_id: "frag_a".to_string(),
                bp: 350,
                topology_form: GelTopologyForm::Linear,
            },
            GelSampleMember {
                seq_id: "frag_b".to_string(),
                bp: 820,
                topology_form: GelTopologyForm::Linear,
            },
            GelSampleMember {
                seq_id: "frag_c".to_string(),
                bp: 1650,
                topology_form: GelTopologyForm::Linear,
            },
        ];
        let ladders = vec![
            "NEB 100bp DNA Ladder".to_string(),
            "NEB 1kb DNA Ladder".to_string(),
        ];
        let layout = build_pool_gel_layout(&members, &ladders, None).unwrap();
        let svg = export_pool_gel_svg(&layout);
        let expected = include_str!("../tests/snapshots/pool_gel/minimal.svg");
        assert_eq!(svg, expected);
    }

    #[test]
    #[cfg(feature = "snapshot-tests")]
    #[ignore]
    fn write_pool_gel_snapshot() {
        let members = vec![
            GelSampleMember {
                seq_id: "frag_a".to_string(),
                bp: 350,
                topology_form: GelTopologyForm::Linear,
            },
            GelSampleMember {
                seq_id: "frag_b".to_string(),
                bp: 820,
                topology_form: GelTopologyForm::Linear,
            },
            GelSampleMember {
                seq_id: "frag_c".to_string(),
                bp: 1650,
                topology_form: GelTopologyForm::Linear,
            },
        ];
        let ladders = vec![
            "NEB 100bp DNA Ladder".to_string(),
            "NEB 1kb DNA Ladder".to_string(),
        ];
        let layout = build_pool_gel_layout(&members, &ladders, None).unwrap();
        let svg = export_pool_gel_svg(&layout);
        fs::create_dir_all("tests/snapshots/pool_gel").unwrap();
        fs::write("tests/snapshots/pool_gel/minimal.svg", svg).unwrap();
    }

    #[test]
    fn test_topology_aware_circular_dna_runs_lower_than_linear_same_bp() {
        let members = vec![
            GelSampleMember {
                seq_id: "linear".to_string(),
                bp: 5000,
                topology_form: GelTopologyForm::Linear,
            },
            GelSampleMember {
                seq_id: "circular".to_string(),
                bp: 5000,
                topology_form: GelTopologyForm::Circular,
            },
        ];
        let layout = build_pool_gel_layout(&members, &[], None).expect("layout");
        let sample_lane = layout
            .lanes
            .iter()
            .find(|lane| !lane.is_ladder)
            .expect("sample lane");
        assert_eq!(sample_lane.bands.len(), 2);
        let circular_band = sample_lane
            .bands
            .iter()
            .find(|band| band.topology_label == "circular")
            .expect("circular band");
        let linear_band = sample_lane
            .bands
            .iter()
            .find(|band| band.topology_label == "linear")
            .expect("linear band");
        assert!(circular_band.apparent_bp < linear_band.apparent_bp);
    }

    #[test]
    fn test_topology_aware_explicit_circular_forms_span_supercoiled_to_nicked() {
        let members = vec![
            GelSampleMember {
                seq_id: "supercoiled".to_string(),
                bp: 5000,
                topology_form: GelTopologyForm::Supercoiled,
            },
            GelSampleMember {
                seq_id: "circular".to_string(),
                bp: 5000,
                topology_form: GelTopologyForm::Circular,
            },
            GelSampleMember {
                seq_id: "linear".to_string(),
                bp: 5000,
                topology_form: GelTopologyForm::Linear,
            },
            GelSampleMember {
                seq_id: "relaxed".to_string(),
                bp: 5000,
                topology_form: GelTopologyForm::RelaxedCircular,
            },
            GelSampleMember {
                seq_id: "nicked".to_string(),
                bp: 5000,
                topology_form: GelTopologyForm::NickedCircular,
            },
        ];
        let layout = build_pool_gel_layout(&members, &[], None).expect("layout");
        let sample_lane = layout
            .lanes
            .iter()
            .find(|lane| !lane.is_ladder)
            .expect("sample lane");
        let band = |label: &str| {
            sample_lane
                .bands
                .iter()
                .find(|band| band.topology_label == label)
                .expect("topology band")
        };
        assert!(band("supercoiled").apparent_bp < band("circular").apparent_bp);
        assert!(band("circular").apparent_bp < band("linear").apparent_bp);
        assert!(band("linear").apparent_bp < band("relaxed circular").apparent_bp);
        assert!(band("relaxed circular").apparent_bp < band("nicked circular").apparent_bp);
    }

    #[test]
    fn test_mass_based_intensity_favors_larger_single_band() {
        let members = vec![
            GelSampleMember {
                seq_id: "small".to_string(),
                bp: 300,
                topology_form: GelTopologyForm::Linear,
            },
            GelSampleMember {
                seq_id: "large".to_string(),
                bp: 5000,
                topology_form: GelTopologyForm::Linear,
            },
        ];
        let layout = build_serial_gel_layout(
            &[
                GelSampleInput {
                    name: "Small".to_string(),
                    role_label: None,
                    members: vec![members[0].clone()],
                },
                GelSampleInput {
                    name: "Large".to_string(),
                    role_label: None,
                    members: vec![members[1].clone()],
                },
            ],
            &[],
            None,
        )
        .expect("layout");
        let sample_lanes = layout
            .lanes
            .iter()
            .filter(|lane| !lane.is_ladder)
            .collect::<Vec<_>>();
        assert!(sample_lanes[1].bands[0].intensity > sample_lanes[0].bands[0].intensity);
    }

    #[test]
    fn test_co_migration_groups_nearby_fragments_into_one_band() {
        let members = vec![
            GelSampleMember {
                seq_id: "frag_a".to_string(),
                bp: 1000,
                topology_form: GelTopologyForm::Linear,
            },
            GelSampleMember {
                seq_id: "frag_b".to_string(),
                bp: 1035,
                topology_form: GelTopologyForm::Linear,
            },
            GelSampleMember {
                seq_id: "frag_c".to_string(),
                bp: 1400,
                topology_form: GelTopologyForm::Linear,
            },
        ];
        let layout = build_pool_gel_layout(&members, &[], None).expect("layout");
        let sample_lane = layout
            .lanes
            .iter()
            .find(|lane| !lane.is_ladder)
            .expect("sample lane");
        assert_eq!(sample_lane.bands.len(), 2);
        let merged_band = sample_lane
            .bands
            .iter()
            .find(|band| band.count == 2)
            .expect("merged band");
        assert_eq!(merged_band.min_bp, 1000);
        assert_eq!(merged_band.bp, 1035);
        assert_eq!(merged_band.topology_label, "linear");
    }

    #[test]
    fn test_export_pool_gel_svg_marks_merged_band_annotation() {
        let layout = build_pool_gel_layout(
            &[
                GelSampleMember {
                    seq_id: "frag_a".to_string(),
                    bp: 1000,
                    topology_form: GelTopologyForm::Linear,
                },
                GelSampleMember {
                    seq_id: "frag_b".to_string(),
                    bp: 1035,
                    topology_form: GelTopologyForm::Linear,
                },
            ],
            &[],
            None,
        )
        .expect("layout");
        let svg = export_pool_gel_svg(&layout);
        assert!(svg.contains("merged x2"));
        assert!(svg.contains("frag_a (1000 bp)"));
        assert!(svg.contains("frag_b (1035 bp)"));
        assert!(svg.contains("Merged-band notes"));
        assert!(
            svg.contains(
                "merged xN means nearby fragments co-migrate under current gel conditions"
            )
        );
    }

    #[test]
    fn test_export_pool_gel_svg_labels_explicit_circular_forms() {
        let layout = build_pool_gel_layout(
            &[
                GelSampleMember {
                    seq_id: "supercoiled_plasmid".to_string(),
                    bp: 5000,
                    topology_form: GelTopologyForm::Supercoiled,
                },
                GelSampleMember {
                    seq_id: "nicked_plasmid".to_string(),
                    bp: 5000,
                    topology_form: GelTopologyForm::NickedCircular,
                },
            ],
            &[],
            None,
        )
        .expect("layout");
        let svg = export_pool_gel_svg(&layout);
        assert!(svg.contains("supercoiled"));
        assert!(svg.contains("nicked circular"));
        assert!(svg.contains("supercoiled_plasmid (5000 bp, supercoiled)"));
        assert!(svg.contains("nicked_plasmid (5000 bp, nicked circular)"));
    }
}
