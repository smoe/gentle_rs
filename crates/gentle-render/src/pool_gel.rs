//! Virtual pool gel model and rendering primitives.

use gentle_protocol::{
    GelBandLabelLayout, GelBufferModel, GelIsoformMarkerMode, GelLaneLabelLayout, GelRunConditions,
    GelTopologyForm, LadderCatalog, PoolGelRenderOptions, default_dna_ladders,
};
use std::collections::{BTreeMap, BTreeSet};
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
const LANE_LABEL_BASELINE: f32 = GEL_BOTTOM + 26.0;
const LANE_LABEL_FONT_SIZE: f32 = 13.0;
const LANE_LABEL_LINE_HEIGHT: f32 = 14.0;
const MIN_ANGLED_LANE_LABEL_DEGREES: f32 = 55.0;
const MAX_ANGLED_LANE_LABEL_DEGREES: f32 = 88.5;
const ISOFORM_LEGEND_ROW_HEIGHT: f32 = 16.0;
const ISOFORM_COLORS: [&str; 12] = [
    "#0072b2", "#d55e00", "#009e73", "#cc79a7", "#56b4e9", "#e69f00", "#f0e442", "#332288",
    "#88ccee", "#aa4499", "#44aa99", "#882255",
];

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
    pub isoform_ids: Vec<String>,
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

#[derive(Clone, Debug, PartialEq)]
enum ResolvedLaneLabelPlacement {
    Horizontal { row: usize },
    Wrapped(Vec<String>),
    Angled { degrees: f32 },
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct GelIsoformIdentity {
    id: String,
    ordinal: usize,
    color: &'static str,
    binary_code: String,
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
                && !picked.iter().any(|n| n == &found)
            {
                picked.push(found);
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

fn normalize_transcript_accession(raw: &str) -> String {
    let accession = raw.trim().to_ascii_uppercase();
    accession
        .rsplit_once('.')
        .filter(|(_, version)| !version.is_empty() && version.chars().all(|ch| ch.is_ascii_digit()))
        .map(|(base, _)| base.to_string())
        .unwrap_or(accession)
}

fn ensembl_transcript_accession(raw: &str) -> Option<String> {
    raw.split(|ch: char| !(ch.is_ascii_alphanumeric() || ch == '.'))
        .filter_map(|token| {
            let upper = token.to_ascii_uppercase();
            if !upper.starts_with("ENS") {
                return None;
            }
            let transcript_marker = upper.char_indices().find_map(|(idx, ch)| {
                (ch == 'T'
                    && upper[idx + ch.len_utf8()..]
                        .chars()
                        .next()
                        .is_some_and(|next| next.is_ascii_digit()))
                .then_some(idx)
            })?;
            (transcript_marker >= 3).then(|| normalize_transcript_accession(&upper))
        })
        .last()
}

fn refseq_transcript_accession(raw: &str) -> Option<String> {
    let upper = raw.to_ascii_uppercase();
    let mut best: Option<(usize, String)> = None;
    for prefix in ["NM_", "NR_", "XM_", "XR_"] {
        let mut search_from = 0usize;
        while let Some(relative_start) = upper[search_from..].find(prefix) {
            let start = search_from + relative_start;
            let prefix_is_bounded =
                start == 0 || !upper.as_bytes()[start - 1].is_ascii_alphanumeric();
            let digits_start = start + prefix.len();
            let mut end = digits_start;
            while end < upper.len() && upper.as_bytes()[end].is_ascii_digit() {
                end += 1;
            }
            if prefix_is_bounded && end > digits_start {
                if end < upper.len() && upper.as_bytes()[end] == b'.' {
                    let version_start = end + 1;
                    let mut version_end = version_start;
                    while version_end < upper.len()
                        && upper.as_bytes()[version_end].is_ascii_digit()
                    {
                        version_end += 1;
                    }
                    if version_end > version_start {
                        end = version_end;
                    }
                }
                let accession = normalize_transcript_accession(&upper[start..end]);
                if best
                    .as_ref()
                    .is_none_or(|(best_start, _)| start > *best_start)
                {
                    best = Some((start, accession));
                }
            }
            search_from = digits_start.min(upper.len());
        }
    }
    best.map(|(_, accession)| accession)
}

fn isoform_identity_from_product_id(raw: &str) -> Option<String> {
    ensembl_transcript_accession(raw).or_else(|| refseq_transcript_accession(raw))
}

fn binary_identity_code(ordinal: usize, width: usize) -> String {
    (0..width.max(1))
        .rev()
        .map(|bit| {
            if ordinal & (1usize << bit) == 0 {
                'O'
            } else {
                'I'
            }
        })
        .collect()
}

fn collect_isoform_identities(
    layout: &PoolGelLayout,
    mode: GelIsoformMarkerMode,
) -> Vec<GelIsoformIdentity> {
    if mode == GelIsoformMarkerMode::Off {
        return vec![];
    }
    let ids = layout
        .lanes
        .iter()
        .filter(|lane| !lane.is_ladder)
        .flat_map(|lane| lane.bands.iter())
        .flat_map(|band| band.isoform_ids.iter().cloned())
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect::<Vec<_>>();
    let code_width = if ids.len() <= 2 {
        1
    } else {
        (usize::BITS - (ids.len() - 1).leading_zeros()) as usize
    };
    ids.into_iter()
        .enumerate()
        .map(|(ordinal, id)| GelIsoformIdentity {
            id,
            ordinal,
            color: ISOFORM_COLORS[ordinal % ISOFORM_COLORS.len()],
            binary_code: binary_identity_code(ordinal, code_width),
        })
        .collect()
}

fn isoform_marker_x(center_x: f32, band_width: f32, ordinal: usize, count: usize) -> f32 {
    center_x - band_width * 0.5
        + band_width * (ordinal.saturating_add(1) as f32 / count.saturating_add(1) as f32)
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
            let isoform_ids = grouped
                .iter()
                .filter_map(|(member, _)| isoform_identity_from_product_id(&member.seq_id))
                .collect::<BTreeSet<_>>()
                .into_iter()
                .collect::<Vec<_>>();
            PoolGelBand {
                bp,
                min_bp,
                apparent_bp,
                intensity: 0.5,
                count,
                estimated_mass_units,
                topology_label,
                labels,
                isoform_ids,
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
                    isoform_ids: vec![],
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
                    isoform_ids: vec![],
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
        if bp.is_multiple_of(1_000) {
            format!("{} kb", bp / 1_000)
        } else {
            format!("{:.1} kb", bp as f32 / 1_000.0)
        }
    } else {
        format!("{bp} bp")
    }
}

fn estimated_monospace_width(text: &str, font_size: f32) -> f32 {
    text.chars().count() as f32 * font_size * 0.62
}

fn split_long_label_token(token: &str, max_chars: usize) -> Vec<String> {
    let chars = token.chars().collect::<Vec<_>>();
    chars
        .chunks(max_chars.max(1))
        .map(|chunk| chunk.iter().collect::<String>())
        .collect()
}

fn wrap_lane_label(label: &str, max_chars: usize) -> Vec<String> {
    let max_chars = max_chars.max(1);
    let mut lines = vec![];
    let mut current = String::new();
    for word in label.split_whitespace() {
        let word_len = word.chars().count();
        if word_len > max_chars {
            if !current.is_empty() {
                lines.push(std::mem::take(&mut current));
            }
            lines.extend(split_long_label_token(word, max_chars));
            continue;
        }
        let candidate_len = if current.is_empty() {
            word_len
        } else {
            current.chars().count() + 1 + word_len
        };
        if candidate_len <= max_chars {
            if !current.is_empty() {
                current.push(' ');
            }
            current.push_str(word);
        } else {
            lines.push(std::mem::take(&mut current));
            current.push_str(word);
        }
    }
    if !current.is_empty() {
        lines.push(current);
    }
    if lines.is_empty() {
        lines.push(label.to_string());
    }
    lines
}

fn angled_lane_label_degrees(label: &str, lane_gap: f32) -> f32 {
    let label_width = estimated_monospace_width(label, LANE_LABEL_FONT_SIZE);
    let available_width = (lane_gap - 8.0).max(1.0);
    if label_width <= available_width {
        return MIN_ANGLED_LANE_LABEL_DEGREES;
    }
    (available_width / label_width)
        .clamp(0.0, 1.0)
        .acos()
        .to_degrees()
        .clamp(MIN_ANGLED_LANE_LABEL_DEGREES, MAX_ANGLED_LANE_LABEL_DEGREES)
}

fn resolve_lane_label_placements(
    layout: &PoolGelLayout,
    lane_gap: f32,
    requested: GelLaneLabelLayout,
) -> Vec<ResolvedLaneLabelPlacement> {
    let horizontal_width = (lane_gap - 8.0).max(24.0);
    let overlong = layout
        .lanes
        .iter()
        .map(|lane| estimated_monospace_width(&lane.name, LANE_LABEL_FONT_SIZE) > horizontal_width)
        .collect::<Vec<_>>();

    layout
        .lanes
        .iter()
        .enumerate()
        .map(|(lane_idx, lane)| match requested {
            GelLaneLabelLayout::Horizontal => ResolvedLaneLabelPlacement::Horizontal { row: 0 },
            GelLaneLabelLayout::Staggered => {
                ResolvedLaneLabelPlacement::Horizontal { row: lane_idx % 2 }
            }
            GelLaneLabelLayout::Angled => ResolvedLaneLabelPlacement::Angled {
                degrees: angled_lane_label_degrees(&lane.name, lane_gap),
            },
            GelLaneLabelLayout::Wrapped => {
                let max_chars = ((lane_gap - 8.0).max(1.0) / (LANE_LABEL_FONT_SIZE * 0.62))
                    .floor()
                    .max(1.0) as usize;
                ResolvedLaneLabelPlacement::Wrapped(wrap_lane_label(&lane.name, max_chars))
            }
            GelLaneLabelLayout::Auto => {
                if !overlong[lane_idx] {
                    return ResolvedLaneLabelPlacement::Horizontal { row: 0 };
                }
                let has_long_neighbor = lane_idx.checked_sub(1).is_some_and(|idx| overlong[idx])
                    || overlong.get(lane_idx + 1).copied().unwrap_or(false);
                if !has_long_neighbor {
                    let max_chars = ((lane_gap - 8.0).max(1.0) / (LANE_LABEL_FONT_SIZE * 0.62))
                        .floor()
                        .max(1.0) as usize;
                    let lines = wrap_lane_label(&lane.name, max_chars);
                    if lines.len() <= 8 {
                        return ResolvedLaneLabelPlacement::Wrapped(lines);
                    }
                }
                ResolvedLaneLabelPlacement::Angled {
                    degrees: angled_lane_label_degrees(&lane.name, lane_gap),
                }
            }
        })
        .collect()
}

fn lane_label_bottom(label: &str, placement: &ResolvedLaneLabelPlacement) -> f32 {
    match placement {
        ResolvedLaneLabelPlacement::Horizontal { row } => {
            LANE_LABEL_BASELINE + *row as f32 * 18.0 + 3.0
        }
        ResolvedLaneLabelPlacement::Wrapped(lines) => {
            LANE_LABEL_BASELINE
                + lines.len().saturating_sub(1) as f32 * LANE_LABEL_LINE_HEIGHT
                + 3.0
        }
        ResolvedLaneLabelPlacement::Angled { degrees } => {
            let radians = degrees.to_radians();
            LANE_LABEL_BASELINE
                + estimated_monospace_width(label, LANE_LABEL_FONT_SIZE) * radians.sin()
                + 3.0
        }
    }
}

fn band_label_fits_inline(label: &str, lane_idx: usize, lane_count: usize, lane_gap: f32) -> bool {
    let x = GEL_LEFT + lane_gap * (lane_idx as f32 + 1.0);
    let text_left = x + 44.0;
    let right_boundary = if lane_idx + 1 < lane_count {
        let next_x = GEL_LEFT + lane_gap * (lane_idx as f32 + 2.0);
        next_x - 42.0
    } else {
        GEL_RIGHT - 8.0
    };
    estimated_monospace_width(label, 11.0) <= (right_boundary - text_left).max(0.0)
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
    export_pool_gel_svg_with_options(layout, &PoolGelRenderOptions::default())
}

pub fn export_pool_gel_svg_with_options(
    layout: &PoolGelLayout,
    render_options: &PoolGelRenderOptions,
) -> String {
    let lane_count = layout.lanes.len().max(1);
    let lane_gap = (GEL_RIGHT - GEL_LEFT) / (lane_count as f32 + 1.0);
    let gel_width = GEL_RIGHT - GEL_LEFT;
    let gel_height = GEL_BOTTOM - GEL_TOP;
    let lane_label_placements =
        resolve_lane_label_placements(layout, lane_gap, render_options.lane_label_layout);
    let label_bottom = layout
        .lanes
        .iter()
        .zip(&lane_label_placements)
        .map(|(lane, placement)| lane_label_bottom(&lane.name, placement))
        .fold(LANE_LABEL_BASELINE, f32::max);
    let role_badge_y = (label_bottom + 5.0).max(GEL_BOTTOM + 34.0);
    let has_role_badges = layout
        .lanes
        .iter()
        .any(|lane| !lane.is_ladder && lane_role_badge(&normalized_lane_role(lane)).is_some());
    let isoform_identities = collect_isoform_identities(layout, render_options.isoform_marker_mode);
    let isoform_identity_indices = isoform_identities
        .iter()
        .enumerate()
        .map(|(idx, identity)| (identity.id.clone(), idx))
        .collect::<BTreeMap<_, _>>();
    let isoform_legend_height = if isoform_identities.is_empty() {
        0.0
    } else {
        52.0 + isoform_identities.len() as f32 * ISOFORM_LEGEND_ROW_HEIGHT
    };
    let content_bottom = if has_role_badges {
        role_badge_y + 16.0
    } else {
        label_bottom
    };
    let svg_height = (SVG_HEIGHT + isoform_legend_height).max(content_bottom + 36.0);
    let detail_panel_height = gel_height + isoform_legend_height;
    let sample_lane_indices = layout
        .lanes
        .iter()
        .enumerate()
        .filter_map(|(idx, lane)| (!lane.is_ladder).then_some(idx))
        .collect::<Vec<_>>();

    let mut doc = Document::new()
        .set("viewBox", (0, 0, SVG_WIDTH, svg_height))
        .set("width", SVG_WIDTH)
        .set("height", svg_height)
        .add(
            Rectangle::new()
                .set("x", 0)
                .set("y", 0)
                .set("width", SVG_WIDTH)
                .set("height", svg_height)
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
                .set("height", detail_panel_height)
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

    for (lane_idx, (lane, label_placement)) in
        layout.lanes.iter().zip(&lane_label_placements).enumerate()
    {
        let x = GEL_LEFT + lane_gap * (lane_idx as f32 + 1.0);
        let lane_fill = if lane.is_ladder { "#1a2028" } else { "#1f252e" };
        doc = doc.add(
            Rectangle::new()
                .set("x", x - 34.0)
                .set("y", GEL_TOP + 10.0)
                .set("width", 68.0)
                .set("height", gel_height - 20.0)
                .set("rx", 6)
                .set("ry", 6)
                .set("fill", lane_fill),
        );
        match label_placement {
            ResolvedLaneLabelPlacement::Horizontal { row } => {
                doc = doc.add(
                    Text::new(lane.name.clone())
                        .set("x", x)
                        .set("y", LANE_LABEL_BASELINE + *row as f32 * 18.0)
                        .set("text-anchor", "middle")
                        .set("font-family", "monospace")
                        .set("font-size", LANE_LABEL_FONT_SIZE)
                        .set("fill", "#0f172a"),
                );
            }
            ResolvedLaneLabelPlacement::Wrapped(lines) => {
                for (line_idx, line) in lines.iter().enumerate() {
                    doc = doc.add(
                        Text::new(line.clone())
                            .set("x", x)
                            .set(
                                "y",
                                LANE_LABEL_BASELINE + line_idx as f32 * LANE_LABEL_LINE_HEIGHT,
                            )
                            .set("text-anchor", "middle")
                            .set("font-family", "monospace")
                            .set("font-size", LANE_LABEL_FONT_SIZE)
                            .set("fill", "#0f172a")
                            .set("data-label-layout", "wrapped")
                            .set("data-full-label", lane.name.clone()),
                    );
                }
            }
            ResolvedLaneLabelPlacement::Angled { degrees } => {
                doc = doc.add(
                    Text::new(lane.name.clone())
                        .set("x", x)
                        .set("y", LANE_LABEL_BASELINE)
                        .set("text-anchor", "start")
                        .set(
                            "transform",
                            format!("rotate({degrees:.1} {x} {LANE_LABEL_BASELINE})"),
                        )
                        .set("font-family", "monospace")
                        .set("font-size", LANE_LABEL_FONT_SIZE)
                        .set("fill", "#0f172a")
                        .set("data-label-layout", "angled"),
                );
            }
        }

        if !lane.is_ladder
            && let Some((badge_text, badge_fill, badge_text_fill)) =
                lane_role_badge(&normalized_lane_role(lane))
        {
            doc = doc
                .add(
                    Rectangle::new()
                        .set("x", x - 28.0)
                        .set("y", role_badge_y)
                        .set("width", 56.0)
                        .set("height", 16.0)
                        .set("rx", 8)
                        .set("ry", 8)
                        .set("fill", badge_fill),
                )
                .add(
                    Text::new(badge_text)
                        .set("x", x)
                        .set("y", role_badge_y + 11.0)
                        .set("text-anchor", "middle")
                        .set("font-family", "monospace")
                        .set("font-size", 9)
                        .set("font-weight", 700)
                        .set("fill", badge_text_fill),
                );
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

            if !lane.is_ladder && !isoform_identities.is_empty() {
                let marker_outer_size =
                    (width / (isoform_identities.len() as f32 + 1.0) * 0.72).clamp(3.0, 6.0);
                let marker_inner_size = (marker_outer_size - 1.8).max(1.4);
                for isoform_id in &band.isoform_ids {
                    let Some(identity_idx) = isoform_identity_indices.get(isoform_id).copied()
                    else {
                        continue;
                    };
                    let identity = &isoform_identities[identity_idx];
                    let marker_x =
                        isoform_marker_x(x, width, identity.ordinal, isoform_identities.len());
                    doc = doc
                        .add(
                            Rectangle::new()
                                .set("x", marker_x - marker_outer_size * 0.5)
                                .set("y", y - marker_outer_size * 0.5)
                                .set("width", marker_outer_size)
                                .set("height", marker_outer_size)
                                .set("rx", 0.8)
                                .set("ry", 0.8)
                                .set("fill", "#f8fafc")
                                .set("stroke", "#111827")
                                .set("stroke-width", 0.6),
                        )
                        .add(
                            Rectangle::new()
                                .set("x", marker_x - marker_inner_size * 0.5)
                                .set("y", y - marker_inner_size * 0.5)
                                .set("width", marker_inner_size)
                                .set("height", marker_inner_size)
                                .set("rx", 0.4)
                                .set("ry", 0.4)
                                .set("fill", identity.color)
                                .set("data-gentle-role", "isoform-identity-marker")
                                .set("data-isoform-id", identity.id.clone())
                                .set("data-isoform-code", identity.binary_code.clone())
                                .set("data-isoform-ordinal", identity.ordinal),
                        );
                }
            }

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
                let show_inline = match render_options.band_label_layout {
                    GelBandLabelLayout::Inline => true,
                    GelBandLabelLayout::Panel => false,
                    GelBandLabelLayout::Auto => {
                        band_label_fits_inline(&label, lane_idx, lane_count, lane_gap)
                    }
                };
                if show_inline {
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
    if !isoform_identities.is_empty() {
        doc = doc
            .add(
                Text::new("Isoform identity key")
                    .set("x", DETAIL_PANEL_LEFT)
                    .set("y", header_y)
                    .set("font-family", "monospace")
                    .set("font-size", 12)
                    .set("font-weight", 700)
                    .set("fill", "#0f172a"),
            )
            .add(
                Text::new("colour + marker position repeat across bands")
                    .set("x", DETAIL_PANEL_LEFT + 4.0)
                    .set("y", header_y + 16.0)
                    .set("font-family", "monospace")
                    .set("font-size", 9)
                    .set("fill", "#475569"),
            )
            .add(
                Text::new("code: O=0, I=1 (colour-independent fallback)")
                    .set("x", DETAIL_PANEL_LEFT + 4.0)
                    .set("y", header_y + 29.0)
                    .set("font-family", "monospace")
                    .set("font-size", 9)
                    .set("fill", "#64748b"),
            );
        header_y += 44.0;
        let legend_band_x = DETAIL_PANEL_LEFT + 8.0;
        let legend_band_width = 54.0;
        for identity in &isoform_identities {
            let row_y = header_y;
            let marker_x = isoform_marker_x(
                legend_band_x + legend_band_width * 0.5,
                legend_band_width,
                identity.ordinal,
                isoform_identities.len(),
            );
            let code_x = legend_band_x + legend_band_width + 10.0;
            let id_x = code_x + estimated_monospace_width(&identity.binary_code, 10.0) + 10.0;
            doc = doc
                .add(
                    Rectangle::new()
                        .set("x", legend_band_x)
                        .set("y", row_y - 5.0)
                        .set("width", legend_band_width)
                        .set("height", 6.0)
                        .set("rx", 2)
                        .set("ry", 2)
                        .set("fill", "#f59e0b")
                        .set("opacity", 0.78),
                )
                .add(
                    Rectangle::new()
                        .set("x", marker_x - 3.0)
                        .set("y", row_y - 5.0)
                        .set("width", 6.0)
                        .set("height", 6.0)
                        .set("rx", 0.8)
                        .set("ry", 0.8)
                        .set("fill", "#f8fafc")
                        .set("stroke", "#111827")
                        .set("stroke-width", 0.6),
                )
                .add(
                    Rectangle::new()
                        .set("x", marker_x - 2.0)
                        .set("y", row_y - 4.0)
                        .set("width", 4.0)
                        .set("height", 4.0)
                        .set("rx", 0.4)
                        .set("ry", 0.4)
                        .set("fill", identity.color)
                        .set("data-gentle-role", "isoform-legend-marker")
                        .set("data-isoform-id", identity.id.clone())
                        .set("data-isoform-code", identity.binary_code.clone()),
                )
                .add(
                    Text::new(identity.binary_code.clone())
                        .set("x", code_x)
                        .set("y", row_y)
                        .set("font-family", "monospace")
                        .set("font-size", 10)
                        .set("font-weight", 700)
                        .set("fill", identity.color),
                )
                .add(
                    Text::new(identity.id.clone())
                        .set("x", id_x)
                        .set("y", row_y)
                        .set("font-family", "monospace")
                        .set("font-size", 10)
                        .set("fill", "#334155"),
                );
            header_y += ISOFORM_LEGEND_ROW_HEIGHT;
        }
        header_y += 8.0;
    }
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
    let detail_bottom = GEL_BOTTOM + isoform_legend_height - 24.0;
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
                let identity = isoform_identity_from_product_id(label)
                    .and_then(|id| isoform_identity_indices.get(&id).copied())
                    .map(|idx| &isoform_identities[idx]);
                let display_label = identity
                    .map(|identity| format!("[{}] {label}", identity.binary_code))
                    .unwrap_or_else(|| label.clone());
                let mut label_text = Text::new(display_label)
                    .set("x", DETAIL_PANEL_LEFT + 14.0)
                    .set("y", detail_y)
                    .set("font-family", "monospace")
                    .set("font-size", 10)
                    .set("fill", "#475569");
                if let Some(identity) = identity {
                    label_text = label_text
                        .set("data-isoform-id", identity.id.clone())
                        .set("data-isoform-code", identity.binary_code.clone());
                }
                doc = doc.add(label_text);
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
            if detail_y > detail_bottom {
                break;
            }
        }
        detail_y += 8.0;
        if detail_y > detail_bottom {
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

    fn dense_labeled_serial_layout() -> PoolGelLayout {
        let samples = (1..=14)
            .map(|index| GelSampleInput {
                name: format!("A{index}"),
                role_label: None,
                members: vec![GelSampleMember {
                    seq_id: format!("amplicon_{index}"),
                    bp: 313 + index * 100,
                    topology_form: GelTopologyForm::Linear,
                }],
            })
            .collect::<Vec<_>>();
        build_serial_gel_layout(
            &samples,
            &["Plasmid Factory 1kb DNA Ladder".to_string()],
            None,
        )
        .expect("dense serial layout")
    }

    fn pou2f2_isoform_layout() -> PoolGelLayout {
        let samples = vec![
            GelSampleInput {
                name: "A1".to_string(),
                role_label: None,
                members: vec![
                    GelSampleMember {
                        seq_id: "POU2F2_A1_ENST00000342301.8_p1_5075bp".to_string(),
                        bp: 5_075,
                        topology_form: GelTopologyForm::Linear,
                    },
                    GelSampleMember {
                        seq_id: "POU2F2_A1_ENST00000389341_p1_5001bp".to_string(),
                        bp: 5_001,
                        topology_form: GelTopologyForm::Linear,
                    },
                ],
            },
            GelSampleInput {
                name: "A2".to_string(),
                role_label: None,
                members: vec![GelSampleMember {
                    seq_id: "POU2F2_A2_ENST00000534559_p1_352bp".to_string(),
                    bp: 352,
                    topology_form: GelTopologyForm::Linear,
                }],
            },
            GelSampleInput {
                name: "A3".to_string(),
                role_label: None,
                members: vec![
                    GelSampleMember {
                        seq_id: "POU2F2_A3_ENST00000342301_p1_832bp".to_string(),
                        bp: 832,
                        topology_form: GelTopologyForm::Linear,
                    },
                    GelSampleMember {
                        seq_id: "POU2F2_A3_ENST00000389341_p1_758bp".to_string(),
                        bp: 758,
                        topology_form: GelTopologyForm::Linear,
                    },
                    GelSampleMember {
                        seq_id: "POU2F2_A3_ENST00000534559_p1_758bp".to_string(),
                        bp: 758,
                        topology_form: GelTopologyForm::Linear,
                    },
                ],
            },
        ];
        build_serial_gel_layout(
            &samples,
            &["Plasmid Factory 1kb DNA Ladder".to_string()],
            None,
        )
        .expect("POU2F2 isoform layout")
    }

    #[test]
    fn dense_gel_auto_wraps_isolated_long_lane_labels_and_hides_crowded_band_text() {
        let layout = dense_labeled_serial_layout();
        let svg = export_pool_gel_svg(&layout);

        assert!(svg.contains("data-label-layout=\"wrapped\""));
        assert!(svg.contains("data-full-label=\"Plasmid Factory 1kb DNA Ladder\""));
        assert!(svg.contains("\nA1\n"));
        assert_eq!(svg.matches("413 bp | linear").count(), 1);
        assert!(svg.contains("413 bp | 413 bp | linear"));
    }

    #[test]
    fn explicit_gel_label_layouts_override_auto_placement() {
        let layout = dense_labeled_serial_layout();
        let angled_panel = export_pool_gel_svg_with_options(
            &layout,
            &PoolGelRenderOptions {
                lane_label_layout: GelLaneLabelLayout::Angled,
                band_label_layout: GelBandLabelLayout::Panel,
                isoform_marker_mode: GelIsoformMarkerMode::Auto,
            },
        );
        assert!(angled_panel.contains("data-label-layout=\"angled\""));
        assert_eq!(angled_panel.matches("413 bp | linear").count(), 1);

        let horizontal_inline = export_pool_gel_svg_with_options(
            &layout,
            &PoolGelRenderOptions {
                lane_label_layout: GelLaneLabelLayout::Horizontal,
                band_label_layout: GelBandLabelLayout::Inline,
                isoform_marker_mode: GelIsoformMarkerMode::Auto,
            },
        );
        assert!(!horizontal_inline.contains("data-label-layout=\"wrapped\""));
        assert_eq!(horizontal_inline.matches("413 bp | linear").count(), 2);
    }

    #[test]
    fn band_label_auto_fit_tracks_lane_spacing() {
        assert!(!band_label_fits_inline(
            "413 bp | linear",
            1,
            16,
            810.0 / 17.0
        ));
        assert!(band_label_fits_inline("413 bp | linear", 1, 3, 810.0 / 4.0));
    }

    #[test]
    fn adaptive_angle_limits_long_label_horizontal_projection_to_lane_width() {
        let lane_gap = 50.0;
        let label = "POU2F2 assay with a deliberately long adjacent lane label";
        let degrees = angled_lane_label_degrees(label, lane_gap);
        let projected_width =
            estimated_monospace_width(label, LANE_LABEL_FONT_SIZE) * degrees.to_radians().cos();

        assert!(degrees > MIN_ANGLED_LANE_LABEL_DEGREES);
        assert!(projected_width <= lane_gap - 8.0 + 0.01);
    }

    #[test]
    fn transcript_identity_detection_normalizes_ensembl_and_refseq_versions() {
        assert_eq!(
            isoform_identity_from_product_id("POU2F2_A1_ENST00000342301.8_p1_5075bp").as_deref(),
            Some("ENST00000342301")
        );
        assert_eq!(
            isoform_identity_from_product_id("assay_NM_001256789.3_product").as_deref(),
            Some("NM_001256789")
        );
        assert_eq!(
            isoform_identity_from_product_id(
                "source_ENST00000011111_cdna_product_ENST00000342301_p1"
            )
            .as_deref(),
            Some("ENST00000342301")
        );
        assert_eq!(isoform_identity_from_product_id("ordinary_plasmid"), None);
    }

    #[test]
    fn pou2f2_gel_repeats_color_position_and_binary_identity_across_bands() {
        let layout = pou2f2_isoform_layout();
        let a1 = layout
            .lanes
            .iter()
            .find(|lane| lane.name == "A1")
            .expect("A1 lane");
        assert!(a1.bands.iter().any(|band| {
            band.isoform_ids == vec!["ENST00000342301".to_string(), "ENST00000389341".to_string()]
        }));

        let svg = export_pool_gel_svg(&layout);
        assert!(svg.contains("Isoform identity key"));
        assert!(svg.contains("code: O=0, I=1"));
        assert_eq!(
            svg.matches("data-gentle-role=\"isoform-identity-marker\"")
                .count(),
            6
        );
        assert!(svg.contains("data-isoform-id=\"ENST00000342301\""));
        assert!(svg.contains("data-isoform-code=\"OO\""));
        assert!(svg.contains("data-isoform-code=\"OI\""));
        assert!(svg.contains("data-isoform-code=\"IO\""));
        assert!(svg.contains("[OO] POU2F2_A1_ENST00000342301.8_p1_5075bp"));
    }

    #[test]
    fn isoform_marker_mode_off_restores_unmarked_gel() {
        let layout = pou2f2_isoform_layout();
        let svg = export_pool_gel_svg_with_options(
            &layout,
            &PoolGelRenderOptions {
                isoform_marker_mode: GelIsoformMarkerMode::Off,
                ..PoolGelRenderOptions::default()
            },
        );
        assert!(!svg.contains("Isoform identity key"));
        assert!(!svg.contains("data-gentle-role=\"isoform-identity-marker\""));
        assert!(!svg.contains("data-isoform-code"));
    }

    #[test]
    fn isoform_marker_position_is_relative_and_stable_across_band_widths() {
        let narrow = isoform_marker_x(100.0, 40.0, 1, 4);
        let wide = isoform_marker_x(200.0, 80.0, 1, 4);
        let narrow_fraction = (narrow - 80.0) / 40.0;
        let wide_fraction = (wide - 160.0) / 80.0;
        assert!((narrow_fraction - wide_fraction).abs() < 1e-6);
        assert!((narrow_fraction - 0.4).abs() < 1e-6);
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
    #[ignore = "regenerates pool-gel SVG snapshot; run manually"]
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
        let members = [
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
