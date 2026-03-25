//! Core operation dispatch handler for `GentleEngine::apply_internal`.
//!
//! This file is the mutation hub that maps one `Operation` variant onto the
//! appropriate helper routines. Keep adapter entry points thin and route new
//! engine behavior through this shared dispatch path.

use super::*;
use crate::{
    gibson_planning::{GibsonAssemblyPlan, derive_gibson_execution_plan},
    uniprot::UniprotNucleotideXref,
};

impl GentleEngine {
    fn dotplot_svg_xml_escape(raw: &str) -> String {
        raw.replace('&', "&amp;")
            .replace('<', "&lt;")
            .replace('>', "&gt;")
            .replace('"', "&quot;")
            .replace('\'', "&apos;")
    }

    fn dotplot_svg_mix_rgb(from: [u8; 3], to: [u8; 3], t: f32) -> [u8; 3] {
        let t = t.clamp(0.0, 1.0);
        let lerp = |a: u8, b: u8| -> u8 {
            ((a as f32) + ((b as f32) - (a as f32)) * t)
                .round()
                .clamp(0.0, 255.0) as u8
        };
        [
            lerp(from[0], to[0]),
            lerp(from[1], to[1]),
            lerp(from[2], to[2]),
        ]
    }

    fn dotplot_sampling_overlap_summary(window_len: usize, step_bp: usize) -> String {
        let safe_step = step_bp.max(1);
        let overlap_bp = window_len.saturating_sub(safe_step);
        let windows_per_base = (window_len.saturating_add(safe_step).saturating_sub(1)) / safe_step;
        if safe_step == 1 {
            format!(
                "adjacent windows overlap by {overlap_bp} bp; each interior base participates in {windows_per_base} consecutive ordered windows (dense sliding)"
            )
        } else if safe_step < window_len {
            format!(
                "adjacent windows overlap by {overlap_bp} bp; each interior base participates in up to {windows_per_base} consecutive ordered windows (subsampled sliding)"
            )
        } else if safe_step == window_len {
            "adjacent windows do not overlap; each interior base participates in at most 1 ordered window (edge-touching sampling)".to_string()
        } else {
            format!(
                "adjacent windows do not overlap and can leave up to {} bp unsampled between starts; each interior base participates in at most 1 ordered window (sparse sampling)",
                safe_step - window_len
            )
        }
    }

    fn build_dotplot_svg_document_for_export(
        view: &DotplotView,
        flex_track: Option<&FlexibilityTrack>,
        density_threshold: f32,
        intensity_gain: f32,
    ) -> String {
        const RENDER_MAX_POINTS: usize = 120_000;
        const CONNECT_DIAGONALS_MAX_CELLS: usize = 80_000;
        let density_threshold = density_threshold.clamp(0.0, 0.99);
        let intensity_gain = intensity_gain.clamp(0.1, 16.0);
        let query_span = view
            .span_end_0based
            .saturating_sub(view.span_start_0based)
            .max(1);
        let query_span_max = query_span.saturating_sub(1).max(1);
        let reference_span = view
            .reference_span_end_0based
            .saturating_sub(view.reference_span_start_0based)
            .max(1);
        let reference_span_max = reference_span.saturating_sub(1).max(1);
        let reference_seq_label = view
            .reference_seq_id
            .as_deref()
            .unwrap_or(view.seq_id.as_str());

        let outer_margin = 18.0_f32;
        let left_margin = 56.0_f32;
        let right_margin = 16.0_f32;
        let top_margin = 26.0_f32;
        let bottom_margin = 24.0_f32;
        let dotplot_width = 1600.0_f32;
        let dotplot_height =
            (dotplot_width * (reference_span as f32 / query_span as f32)).clamp(420.0, 1280.0);
        let flex_height = if flex_track.is_some() { 108.0 } else { 0.0 };
        let gap = if flex_track.is_some() { 10.0 } else { 0.0 };
        let canvas_width = outer_margin + left_margin + dotplot_width + right_margin + outer_margin;
        let canvas_height = outer_margin
            + top_margin
            + dotplot_height
            + gap
            + flex_height
            + bottom_margin
            + outer_margin;
        let dotplot_left = outer_margin + left_margin;
        let dotplot_top = outer_margin + top_margin;
        let dotplot_right = dotplot_left + dotplot_width;
        let dotplot_bottom = dotplot_top + dotplot_height;
        let cols = dotplot_width.max(2.0).round() as i32;
        let rows = dotplot_height.max(2.0).round() as i32;
        let sample_stride = (view.points.len() / RENDER_MAX_POINTS).max(1);

        let mut cells: HashMap<(i32, i32), (usize, usize)> = HashMap::new();
        for point in view.points.iter().step_by(sample_stride) {
            let x_local = point
                .x_0based
                .saturating_sub(view.span_start_0based)
                .min(query_span_max);
            let y_local = point
                .y_0based
                .saturating_sub(view.reference_span_start_0based)
                .min(reference_span_max);
            let x_frac = (x_local as f32 / query_span_max as f32).clamp(0.0, 1.0);
            let y_frac = (y_local as f32 / reference_span_max as f32).clamp(0.0, 1.0);
            let x_cell = ((x_frac * (cols - 1) as f32).round() as i32).clamp(0, cols - 1);
            let y_cell = ((y_frac * (rows - 1) as f32).round() as i32).clamp(0, rows - 1);
            let entry = cells
                .entry((x_cell, y_cell))
                .or_insert((0, point.mismatches));
            entry.0 = entry.0.saturating_add(1);
            entry.1 = entry.1.min(point.mismatches);
        }
        let max_cell_count = cells.values().map(|(count, _)| *count).max().unwrap_or(1) as f32;
        let mut visible_cells: Vec<((i32, i32), [u8; 3], f32)> = vec![];
        for ((x_cell, y_cell), (count, min_mismatch)) in &cells {
            let density_raw = (*count as f32 / max_cell_count).clamp(0.0, 1.0);
            if density_raw < density_threshold {
                continue;
            }
            let normalized = if density_threshold <= 0.0 {
                density_raw
            } else {
                ((density_raw - density_threshold) / (1.0 - density_threshold)).clamp(0.0, 1.0)
            };
            let density = (normalized * intensity_gain).clamp(0.0, 1.0).sqrt();
            let mismatch_frac = if view.max_mismatches == 0 {
                0.0
            } else {
                (*min_mismatch as f32 / view.max_mismatches as f32).clamp(0.0, 1.0)
            };
            let rgb = Self::dotplot_svg_mix_rgb([29, 78, 216], [180, 83, 9], mismatch_frac);
            let opacity = ((90.0 + 165.0 * density) / 255.0).clamp(0.0, 1.0);
            visible_cells.push(((*x_cell, *y_cell), rgb, opacity));
        }
        visible_cells.sort_by_key(|((x, y), _, _)| (*y, *x));
        let visible_lookup: HashMap<(i32, i32), ([u8; 3], f32)> = visible_cells
            .iter()
            .map(|(cell, rgb, opacity)| (*cell, (*rgb, *opacity)))
            .collect();

        let mut svg = String::new();
        svg.push_str(&format!(
            "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{:.0}\" height=\"{:.0}\" viewBox=\"0 0 {:.0} {:.0}\">",
            canvas_width, canvas_height, canvas_width, canvas_height
        ));
        svg.push_str(&format!(
            "<rect x=\"0\" y=\"0\" width=\"{:.0}\" height=\"{:.0}\" fill=\"#ffffff\"/>",
            canvas_width, canvas_height
        ));
        svg.push_str(&format!(
            "<rect x=\"{:.1}\" y=\"{:.1}\" width=\"{:.1}\" height=\"{:.1}\" fill=\"#f8fafc\" stroke=\"#cbd5e1\" stroke-width=\"1\" rx=\"6\" ry=\"6\"/>",
            outer_margin,
            outer_margin,
            canvas_width - outer_margin * 2.0,
            canvas_height - outer_margin * 2.0
        ));
        svg.push_str(&format!(
            "<rect x=\"{:.1}\" y=\"{:.1}\" width=\"{:.1}\" height=\"{:.1}\" fill=\"#ffffff\" stroke=\"#cbd5e1\" stroke-width=\"1\" rx=\"4\" ry=\"4\"/>",
            dotplot_left, dotplot_top, dotplot_width, dotplot_height
        ));

        for tick_idx in 1..10 {
            let frac = tick_idx as f32 / 10.0;
            let gx = dotplot_left + frac * dotplot_width;
            let gy = dotplot_top + frac * dotplot_height;
            svg.push_str(&format!(
                "<line x1=\"{:.2}\" y1=\"{:.2}\" x2=\"{:.2}\" y2=\"{:.2}\" stroke=\"#f1f5f9\" stroke-width=\"1\"/>",
                gx, dotplot_top, gx, dotplot_bottom
            ));
            svg.push_str(&format!(
                "<line x1=\"{:.2}\" y1=\"{:.2}\" x2=\"{:.2}\" y2=\"{:.2}\" stroke=\"#f1f5f9\" stroke-width=\"1\"/>",
                dotplot_left, gy, dotplot_right, gy
            ));
        }

        let cell_w = dotplot_width / cols as f32;
        let cell_h = dotplot_height / rows as f32;
        for ((x_cell, y_cell), rgb, opacity) in &visible_cells {
            let x0 = dotplot_left + (*x_cell as f32 / cols as f32) * dotplot_width;
            let y0 = dotplot_top + (*y_cell as f32 / rows as f32) * dotplot_height;
            let x1 = dotplot_left + ((*x_cell + 1) as f32 / cols as f32) * dotplot_width;
            let y1 = dotplot_top + ((*y_cell + 1) as f32 / rows as f32) * dotplot_height;
            let width = (x1 - x0).max(0.8);
            let height = (y1 - y0).max(0.8);
            svg.push_str(&format!(
                "<rect x=\"{:.2}\" y=\"{:.2}\" width=\"{:.2}\" height=\"{:.2}\" fill=\"rgb({},{},{})\" fill-opacity=\"{:.4}\"/>",
                x0, y0, width, height, rgb[0], rgb[1], rgb[2], opacity
            ));
        }

        if !visible_cells.is_empty() && visible_cells.len() <= CONNECT_DIAGONALS_MAX_CELLS {
            let stroke_width = (cell_w.min(cell_h) * 0.45).clamp(0.8, 2.0);
            for ((x_cell, y_cell), rgb, opacity) in &visible_cells {
                for dy in [-1, 0, 1] {
                    let neighbor = (*x_cell + 1, *y_cell + dy);
                    if visible_lookup.contains_key(&neighbor) {
                        let x0 = dotplot_left + (*x_cell as f32 + 0.5) * cell_w;
                        let y0 = dotplot_top + (*y_cell as f32 + 0.5) * cell_h;
                        let x1 = dotplot_left + (neighbor.0 as f32 + 0.5) * cell_w;
                        let y1 = dotplot_top + (neighbor.1 as f32 + 0.5) * cell_h;
                        svg.push_str(&format!(
                            "<line x1=\"{:.2}\" y1=\"{:.2}\" x2=\"{:.2}\" y2=\"{:.2}\" stroke=\"rgb({},{},{})\" stroke-opacity=\"{:.4}\" stroke-width=\"{:.2}\"/>",
                            x0, y0, x1, y1, rgb[0], rgb[1], rgb[2], (opacity * 0.55).clamp(0.0, 1.0), stroke_width
                        ));
                    }
                }
            }
        }

        if visible_cells.is_empty() {
            svg.push_str(&format!(
                "<text x=\"{:.2}\" y=\"{:.2}\" text-anchor=\"middle\" dominant-baseline=\"middle\" font-family=\"monospace\" font-size=\"14\" fill=\"#64748b\">{}</text>",
                (dotplot_left + dotplot_right) * 0.5,
                (dotplot_top + dotplot_bottom) * 0.5,
                Self::dotplot_svg_xml_escape(&format!(
                    "No visible cells (threshold={:.2}, points={})",
                    density_threshold, view.point_count
                ))
            ));
        }

        let header = format!(
            "Dotplot workspace export: {} | mode={} | query={} [{}..{}] | reference={} [{}..{}] | word={} step={} mismatches={} | points={} | threshold={:.2} gain={:.2}",
            view.dotplot_id,
            view.mode.as_str(),
            view.seq_id,
            view.span_start_0based.saturating_add(1),
            view.span_end_0based,
            reference_seq_label,
            view.reference_span_start_0based.saturating_add(1),
            view.reference_span_end_0based,
            view.word_size,
            view.step_bp,
            view.max_mismatches,
            view.point_count,
            density_threshold,
            intensity_gain
        );
        let sampling_line = format!(
            "{} | rendered_cells={} sampled_points={} sample_stride={}",
            Self::dotplot_sampling_overlap_summary(view.word_size.max(1), view.step_bp.max(1)),
            visible_cells.len(),
            view.points.len().div_ceil(sample_stride),
            sample_stride
        );
        svg.push_str(&format!(
            "<text x=\"{:.1}\" y=\"{:.1}\" font-family=\"monospace\" font-size=\"13\" fill=\"#0f172a\">{}</text>",
            dotplot_left,
            outer_margin + 16.0,
            Self::dotplot_svg_xml_escape(&header)
        ));
        svg.push_str(&format!(
            "<text x=\"{:.1}\" y=\"{:.1}\" font-family=\"monospace\" font-size=\"10\" fill=\"#64748b\">{}</text>",
            dotplot_left,
            outer_margin + 31.0,
            Self::dotplot_svg_xml_escape(&sampling_line)
        ));
        svg.push_str(&format!(
            "<text x=\"{:.1}\" y=\"{:.1}\" font-family=\"monospace\" font-size=\"10\" fill=\"#334155\">x: {}</text>",
            dotplot_left,
            outer_margin + 45.0,
            Self::dotplot_svg_xml_escape(view.seq_id.as_str())
        ));
        svg.push_str(&format!(
            "<text x=\"{:.1}\" y=\"{:.1}\" text-anchor=\"end\" font-family=\"monospace\" font-size=\"10\" fill=\"#334155\">y: {}</text>",
            dotplot_right,
            outer_margin + 45.0,
            Self::dotplot_svg_xml_escape(reference_seq_label)
        ));
        svg.push_str(&format!(
            "<text x=\"{:.1}\" y=\"{:.1}\" font-family=\"monospace\" font-size=\"10\" fill=\"#475569\">{}</text>",
            dotplot_left,
            dotplot_bottom + 14.0,
            view.span_start_0based.saturating_add(1)
        ));
        svg.push_str(&format!(
            "<text x=\"{:.1}\" y=\"{:.1}\" text-anchor=\"end\" font-family=\"monospace\" font-size=\"10\" fill=\"#475569\">{}</text>",
            dotplot_right,
            dotplot_bottom + 14.0,
            view.span_end_0based
        ));
        svg.push_str(&format!(
            "<text x=\"{:.1}\" y=\"{:.1}\" text-anchor=\"end\" font-family=\"monospace\" font-size=\"10\" fill=\"#475569\">{}</text>",
            dotplot_left - 6.0,
            dotplot_top + 10.0,
            view.reference_span_start_0based.saturating_add(1)
        ));
        svg.push_str(&format!(
            "<text x=\"{:.1}\" y=\"{:.1}\" text-anchor=\"end\" font-family=\"monospace\" font-size=\"10\" fill=\"#475569\">{}</text>",
            dotplot_left - 6.0,
            dotplot_bottom,
            view.reference_span_end_0based
        ));

        if let Some(track) = flex_track {
            let flex_top = dotplot_bottom + gap;
            let flex_bottom = canvas_height - (outer_margin + bottom_margin);
            let flex_height_px = (flex_bottom - flex_top).max(10.0);
            svg.push_str(&format!(
                "<rect x=\"{:.1}\" y=\"{:.1}\" width=\"{:.1}\" height=\"{:.1}\" fill=\"#ffffff\" stroke=\"#cbd5e1\" stroke-width=\"1\" rx=\"4\" ry=\"4\"/>",
                dotplot_left, flex_top, dotplot_width, flex_height_px
            ));
            let score_span = (track.max_score - track.min_score).abs().max(1e-12);
            let span_bp = track
                .span_end_0based
                .saturating_sub(track.span_start_0based)
                .max(1) as f32;
            let mut bins: Vec<&_> = track.bins.iter().collect();
            bins.sort_by_key(|bin| bin.start_0based);
            let mut points: Vec<(f32, f32)> = Vec::with_capacity(bins.len());
            for bin in bins {
                let x_fraction = (bin.start_0based.saturating_sub(track.span_start_0based) as f32
                    / span_bp)
                    .clamp(0.0, 1.0);
                let y_fraction =
                    ((bin.score - track.min_score) / score_span).clamp(0.0, 1.0) as f32;
                let x = dotplot_left + x_fraction * dotplot_width;
                let y = flex_bottom - y_fraction * flex_height_px;
                points.push((x, y));
            }
            for pair in points.windows(2) {
                let (left_x, left_y) = pair[0];
                let (right_x, right_y) = pair[1];
                svg.push_str(&format!(
                    "<line x1=\"{:.2}\" y1=\"{:.2}\" x2=\"{:.2}\" y2=\"{:.2}\" stroke=\"#166534\" stroke-width=\"1.5\"/>",
                    left_x, left_y, right_x, right_y
                ));
            }
            svg.push_str(&format!(
                "<text x=\"{:.1}\" y=\"{:.1}\" font-family=\"monospace\" font-size=\"10\" fill=\"#0f172a\">{}</text>",
                dotplot_left + 6.0,
                flex_top + 14.0,
                Self::dotplot_svg_xml_escape(&format!(
                    "{} [{}..{}] min={:.3} max={:.3}",
                    track.model.as_str(),
                    track.span_start_0based.saturating_add(1),
                    track.span_end_0based,
                    track.min_score,
                    track.max_score
                ))
            ));
        }

        svg.push_str(&format!(
            "<text x=\"{:.1}\" y=\"{:.1}\" text-anchor=\"end\" font-family=\"monospace\" font-size=\"10\" fill=\"#64748b\">GENtle dotplot SVG export</text>",
            canvas_width - outer_margin,
            canvas_height - 6.0
        ));
        svg.push_str("</svg>");
        svg
    }

    fn render_dotplot_svg_to_path(
        &self,
        seq_id: &str,
        dotplot_id: &str,
        path: &str,
        flex_track_id: Option<&str>,
        display_density_threshold: Option<f32>,
        display_intensity_gain: Option<f32>,
    ) -> Result<(), EngineError> {
        let normalized_dotplot_id = dotplot_id.trim();
        if normalized_dotplot_id.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "RenderDotplotSvg requires dotplot_id".to_string(),
            });
        }
        let view = self.get_dotplot_view(normalized_dotplot_id)?;
        if !view.seq_id.eq_ignore_ascii_case(seq_id) {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Dotplot '{}' belongs to seq_id='{}', not '{}'",
                    normalized_dotplot_id, view.seq_id, seq_id
                ),
            });
        }
        let flex_track = if let Some(flex_track_id) = flex_track_id
            .map(str::trim)
            .filter(|value| !value.is_empty())
        {
            let track = self.get_flexibility_track(flex_track_id)?;
            if !track.seq_id.eq_ignore_ascii_case(seq_id) {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!(
                        "Flexibility track '{}' belongs to seq_id='{}', not '{}'",
                        flex_track_id, track.seq_id, seq_id
                    ),
                });
            }
            Some(track)
        } else {
            None
        };
        let density_threshold = display_density_threshold.unwrap_or(0.0);
        let intensity_gain = display_intensity_gain.unwrap_or(1.0);
        let svg = Self::build_dotplot_svg_document_for_export(
            &view,
            flex_track.as_ref(),
            density_threshold,
            intensity_gain,
        );
        std::fs::write(path, svg).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!("Could not write SVG output '{path}': {e}"),
        })?;
        Ok(())
    }

    fn is_transcript_feature_for_derivation(feature: &gb_io::seq::Feature) -> bool {
        let kind = feature.kind.to_string();
        kind.eq_ignore_ascii_case("mRNA") || kind.eq_ignore_ascii_case("transcript")
    }

    fn qualifier_text_for_derivation(feature: &gb_io::seq::Feature, key: &str) -> Option<String> {
        feature
            .qualifier_values(key.into())
            .map(|value| value.split_whitespace().collect::<Vec<_>>().join(" "))
            .map(|value| value.trim().to_string())
            .find(|value| !value.is_empty())
    }

    fn first_nonempty_qualifier_for_derivation(
        feature: &gb_io::seq::Feature,
        keys: &[&str],
    ) -> Option<String> {
        for key in keys {
            if let Some(value) = Self::qualifier_text_for_derivation(feature, key) {
                return Some(value);
            }
        }
        None
    }

    fn transcript_feature_ids_for_derivation(
        &self,
        seq_id: &str,
        feature_ids: &[usize],
        scope: Option<SplicingScopePreset>,
    ) -> Result<Vec<usize>, EngineError> {
        let dna = self
            .state
            .sequences
            .get(seq_id)
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!("Sequence '{seq_id}' not found"),
            })?;
        let features = dna.features();

        if let Some(scope) = scope {
            if feature_ids.len() != 1 {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message:
                        "DeriveTranscriptSequences with scope requires exactly one seed feature id"
                            .to_string(),
                });
            }
            let seed_feature_id = feature_ids[0];
            if seed_feature_id >= features.len() {
                return Err(EngineError {
                    code: ErrorCode::NotFound,
                    message: format!(
                        "Feature id '{}' was not found in sequence '{}'",
                        seed_feature_id, seq_id
                    ),
                });
            }
            let expert = self.inspect_feature_expert(
                seq_id,
                &FeatureExpertTarget::SplicingFeature {
                    feature_id: seed_feature_id,
                    scope,
                },
            )?;
            let splicing = match expert {
                FeatureExpertView::Splicing(view) => view,
                _ => {
                    return Err(EngineError {
                        code: ErrorCode::Internal,
                        message:
                            "Unexpected expert-view payload while deriving transcript sequences"
                                .to_string(),
                    });
                }
            };
            let mut ids = splicing
                .transcripts
                .iter()
                .map(|row| row.transcript_feature_id)
                .collect::<Vec<_>>();
            ids.sort_unstable();
            ids.dedup();
            if ids.is_empty() {
                return Err(EngineError {
                    code: ErrorCode::NotFound,
                    message: format!(
                        "No mRNA transcript features matched splicing scope '{}' for seed feature {}",
                        scope.as_str(),
                        seed_feature_id + 1
                    ),
                });
            }
            return Ok(ids);
        }

        if feature_ids.is_empty() {
            let mut ids = features
                .iter()
                .enumerate()
                .filter_map(|(idx, feature)| {
                    Self::is_transcript_feature_for_derivation(feature).then_some(idx)
                })
                .collect::<Vec<_>>();
            ids.sort_unstable();
            ids.dedup();
            if ids.is_empty() {
                return Err(EngineError {
                    code: ErrorCode::NotFound,
                    message: format!("Sequence '{}' has no mRNA/transcript features", seq_id),
                });
            }
            return Ok(ids);
        }

        let mut ids = feature_ids.to_vec();
        ids.sort_unstable();
        ids.dedup();
        for feature_id in &ids {
            let feature = features.get(*feature_id).ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!(
                    "Feature id '{}' was not found in sequence '{}'",
                    feature_id, seq_id
                ),
            })?;
            if !Self::is_transcript_feature_for_derivation(feature) {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!(
                        "Feature n-{} in '{}' is not an mRNA/transcript feature",
                        feature_id + 1,
                        seq_id
                    ),
                });
            }
        }
        Ok(ids)
    }

    fn derive_transcript_sequence_from_feature(
        source_sequence_upper: &[u8],
        source_feature: &gb_io::seq::Feature,
        source_feature_id: usize,
        source_seq_id: &str,
    ) -> Result<(DNAsequence, String, String, bool, usize), EngineError> {
        let mut exon_ranges: Vec<(usize, usize)> = vec![];
        collect_location_ranges_usize(&source_feature.location, &mut exon_ranges);
        if exon_ranges.is_empty() {
            let (from, to) = source_feature
                .location
                .find_bounds()
                .map_err(|e| EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!(
                        "Could not parse transcript feature n-{} location: {e}",
                        source_feature_id + 1
                    ),
                })?;
            if from >= 0 && to >= 0 {
                exon_ranges.push((from as usize, to as usize));
            }
        }
        exon_ranges.sort_unstable_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
        exon_ranges.dedup();
        exon_ranges.retain(|(start, end)| *end > *start);
        if exon_ranges.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Transcript feature n-{} has no usable exon ranges",
                    source_feature_id + 1
                ),
            });
        }

        let source_len = source_sequence_upper.len();
        let mut assembled: Vec<u8> = vec![];
        let mut exon_segments: Vec<(usize, usize, usize, usize)> = vec![];
        for (start_0based, end_0based_exclusive) in &exon_ranges {
            if *end_0based_exclusive > source_len {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!(
                        "Transcript feature n-{} exon range {}..{} exceeds source sequence length {}",
                        source_feature_id + 1,
                        start_0based + 1,
                        end_0based_exclusive,
                        source_len
                    ),
                });
            }
            let local_start = assembled.len();
            assembled
                .extend_from_slice(&source_sequence_upper[*start_0based..*end_0based_exclusive]);
            let local_end = assembled.len();
            exon_segments.push((*start_0based, *end_0based_exclusive, local_start, local_end));
        }
        if assembled.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Transcript feature n-{} produced an empty spliced sequence",
                    source_feature_id + 1
                ),
            });
        }

        let is_reverse = feature_is_reverse(source_feature);
        let mut derived_sequence = String::from_utf8_lossy(&assembled).to_ascii_uppercase();
        if is_reverse {
            derived_sequence = Self::reverse_complement(&derived_sequence);
            let total_len = derived_sequence.len();
            for segment in &mut exon_segments {
                let new_start = total_len.saturating_sub(segment.3);
                let new_end = total_len.saturating_sub(segment.2);
                segment.2 = new_start;
                segment.3 = new_end;
            }
        }
        exon_segments.sort_unstable_by(|a, b| a.2.cmp(&b.2).then(a.3.cmp(&b.3)));

        let transcript_id = Self::first_nonempty_qualifier_for_derivation(
            source_feature,
            &[
                "transcript_id",
                "standard_name",
                "label",
                "name",
                "product",
                "gene",
            ],
        )
        .unwrap_or_else(|| format!("transcript_{}", source_feature_id + 1));
        let transcript_label = Self::first_nonempty_qualifier_for_derivation(
            source_feature,
            &[
                "label",
                "name",
                "standard_name",
                "product",
                "transcript_id",
                "gene",
            ],
        )
        .unwrap_or_else(|| transcript_id.clone());

        let mut derived =
            DNAsequence::from_sequence(&derived_sequence).map_err(|e| EngineError {
                code: ErrorCode::Internal,
                message: format!(
                    "Could not construct derived transcript sequence for feature n-{}: {e}",
                    source_feature_id + 1
                ),
            })?;
        derived.set_name(transcript_label.clone());
        let total_len = derived_sequence.len();
        let strand_text = if is_reverse { "-" } else { "+" }.to_string();

        let mut transcript_qualifiers = vec![
            ("transcript_id".into(), Some(transcript_id.clone())),
            ("label".into(), Some(transcript_label.clone())),
            ("source_seq_id".into(), Some(source_seq_id.to_string())),
            (
                "source_feature_id".into(),
                Some((source_feature_id + 1).to_string()),
            ),
            (
                "synthetic_origin".into(),
                Some("mrna_transcript_derived".to_string()),
            ),
            ("strand".into(), Some(strand_text.clone())),
        ];
        for key in ["gene", "gene_id", "locus_tag", "note"] {
            if let Some(value) = Self::qualifier_text_for_derivation(source_feature, key) {
                transcript_qualifiers.push((key.into(), Some(value)));
            }
        }
        derived.features_mut().push(gb_io::seq::Feature {
            kind: gb_io::seq::FeatureKind::from("mRNA"),
            location: gb_io::seq::Location::simple_range(0, total_len as i64),
            qualifiers: transcript_qualifiers,
        });

        for (idx, segment) in exon_segments.iter().enumerate() {
            let mut exon_qualifiers = vec![
                (
                    "exon_number".into(),
                    Some((idx.saturating_add(1)).to_string()),
                ),
                ("transcript_id".into(), Some(transcript_id.clone())),
                ("source_seq_id".into(), Some(source_seq_id.to_string())),
                (
                    "source_feature_id".into(),
                    Some((source_feature_id + 1).to_string()),
                ),
                (
                    "genomic_start_1based".into(),
                    Some((segment.0 + 1).to_string()),
                ),
                ("genomic_end_1based".into(), Some(segment.1.to_string())),
                (
                    "synthetic_origin".into(),
                    Some("mrna_transcript_derived".to_string()),
                ),
                ("strand".into(), Some(strand_text.clone())),
            ];
            for key in ["gene", "gene_id", "locus_tag"] {
                if let Some(value) = Self::qualifier_text_for_derivation(source_feature, key) {
                    exon_qualifiers.push((key.into(), Some(value)));
                }
            }
            derived.features_mut().push(gb_io::seq::Feature {
                kind: gb_io::seq::FeatureKind::from("exon"),
                location: gb_io::seq::Location::simple_range(segment.2 as i64, segment.3 as i64),
                qualifiers: exon_qualifiers,
            });
        }

        Self::prepare_sequence(&mut derived);
        Ok((
            derived,
            transcript_id,
            transcript_label,
            is_reverse,
            exon_segments.len(),
        ))
    }

    fn choose_uniprot_nucleotide_xref(
        entry: &UniprotEntry,
        accession_override: Option<&str>,
    ) -> Result<UniprotNucleotideXref, EngineError> {
        if entry.nucleotide_xrefs.is_empty() {
            return Err(EngineError {
                code: ErrorCode::NotFound,
                message: format!(
                    "UniProt entry '{}' has no EMBL/GenBank nucleotide cross-reference",
                    entry.entry_id
                ),
            });
        }

        if let Some(accession_override) = accession_override
            .map(str::trim)
            .filter(|value| !value.is_empty())
        {
            if let Some(found) = entry
                .nucleotide_xrefs
                .iter()
                .find(|xref| xref.accession.eq_ignore_ascii_case(accession_override))
            {
                return Ok(found.clone());
            }
            let available = entry
                .nucleotide_xrefs
                .iter()
                .map(|xref| xref.accession.as_str())
                .collect::<Vec<_>>()
                .join(", ");
            return Err(EngineError {
                code: ErrorCode::NotFound,
                message: format!(
                    "UniProt entry '{}' does not contain nucleotide accession '{}' (available: {})",
                    entry.entry_id, accession_override, available
                ),
            });
        }

        entry
            .nucleotide_xrefs
            .iter()
            .find(|xref| {
                xref.molecule_type
                    .as_deref()
                    .map(|value| value.eq_ignore_ascii_case("mRNA"))
                    .unwrap_or(false)
            })
            .or_else(|| {
                entry.nucleotide_xrefs.iter().find(|xref| {
                    xref.molecule_type
                        .as_deref()
                        .map(|value| value.eq_ignore_ascii_case("Genomic_DNA"))
                        .unwrap_or(false)
                })
            })
            .or_else(|| entry.nucleotide_xrefs.first())
            .cloned()
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!(
                    "UniProt entry '{}' has no usable EMBL/GenBank nucleotide cross-reference",
                    entry.entry_id
                ),
            })
    }

    fn fetch_genbank_accession_into_state(
        &mut self,
        result: &mut OpResult,
        accession_trimmed: &str,
        as_id: Option<SeqId>,
        provenance_operation: &str,
        source_note: Option<&str>,
    ) -> Result<SeqId, EngineError> {
        let (source_url, text) = Self::fetch_genbank_accession_text(accession_trimmed)?;
        let mut tmp = tempfile::Builder::new()
            .prefix("gentle_genbank_fetch_")
            .suffix(".gb")
            .tempfile()
            .map_err(|e| EngineError {
                code: ErrorCode::Io,
                message: format!(
                    "Could not create temporary GenBank file for '{}': {e}",
                    accession_trimmed
                ),
            })?;
        tmp.write_all(text.as_bytes()).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!(
                "Could not write temporary GenBank file for '{}': {e}",
                accession_trimmed
            ),
        })?;
        let path = tmp.path().to_string_lossy().to_string();
        let mut dna = GENtleApp::load_from_file(&path).map_err(|e| EngineError {
            code: ErrorCode::InvalidInput,
            message: format!(
                "Could not parse fetched GenBank accession '{}' from '{}': {e}",
                accession_trimmed, source_url
            ),
        })?;
        Self::prepare_sequence(&mut dna);
        let base = as_id.unwrap_or_else(|| accession_trimmed.to_string());
        let seq_id = self.unique_seq_id(&base);
        self.state.sequences.insert(seq_id.clone(), dna);
        self.add_lineage_node(
            &seq_id,
            SequenceOrigin::ImportedGenomic,
            Some(&result.op_id),
        );

        let imported_anchor = self
            .state
            .sequences
            .get(&seq_id)
            .and_then(|loaded| Self::infer_imported_genbank_anchor(&path, loaded));
        if let Some(anchor) = imported_anchor {
            let mut anchor_verified: Option<bool> = None;
            if let Some(loaded) = self.state.sequences.get(&seq_id) {
                match Self::verify_anchor_sequence_against_catalog(
                    loaded,
                    &anchor,
                    DEFAULT_GENOME_CATALOG_PATH,
                    None,
                ) {
                    Ok(is_match) => {
                        anchor_verified = Some(is_match);
                        if is_match {
                            result.messages.push(format!(
                                "Verified imported GenBank anchor '{}' against catalog '{}' ({}:{}-{})",
                                seq_id,
                                DEFAULT_GENOME_CATALOG_PATH,
                                anchor.genome_id,
                                anchor.chromosome,
                                anchor.start_1based
                            ));
                        } else {
                            result.warnings.push(format!(
                                "Imported GenBank anchor '{}' does not match catalog sequence at {}:{}:{}-{} (catalog='{}')",
                                seq_id,
                                anchor.genome_id,
                                anchor.chromosome,
                                anchor.start_1based,
                                anchor.end_1based,
                                DEFAULT_GENOME_CATALOG_PATH
                            ));
                        }
                    }
                    Err(err) => {
                        result.warnings.push(format!(
                            "Could not verify imported GenBank anchor '{}' against catalog '{}': {}",
                            seq_id, DEFAULT_GENOME_CATALOG_PATH, err
                        ));
                    }
                }
            }
            self.append_genome_extraction_provenance(GenomeExtractionProvenance {
                seq_id: seq_id.clone(),
                recorded_at_unix_ms: Self::now_unix_ms(),
                operation: provenance_operation.to_string(),
                genome_id: anchor.genome_id.clone(),
                catalog_path: DEFAULT_GENOME_CATALOG_PATH.to_string(),
                cache_dir: None,
                chromosome: Some(anchor.chromosome.clone()),
                start_1based: Some(anchor.start_1based),
                end_1based: Some(anchor.end_1based),
                gene_query: None,
                occurrence: None,
                gene_id: None,
                gene_name: None,
                strand: None,
                anchor_strand: anchor.strand,
                anchor_verified,
                sequence_source_type: Some("genbank_accession".to_string()),
                annotation_source_type: Some("genbank_accession".to_string()),
                sequence_source: Some(source_url.clone()),
                annotation_source: Some(source_url.clone()),
                sequence_sha1: None,
                annotation_sha1: None,
            });
            let strand = anchor.strand.unwrap_or('+');
            let verification_label = match anchor_verified {
                Some(true) => "verified",
                Some(false) => "unverified",
                None => "verification n/a",
            };
            result.messages.push(format!(
                "Detected GenBank genome anchor for '{}': {}:{}-{} ({}, strand {}, {})",
                seq_id,
                anchor.chromosome,
                anchor.start_1based,
                anchor.end_1based,
                anchor.genome_id,
                strand,
                verification_label
            ));
        }

        result.created_seq_ids.push(seq_id.clone());
        result.messages.push(format!(
            "Fetched GenBank accession '{}' from '{}' as '{}'",
            accession_trimmed, source_url, seq_id
        ));
        if let Some(source_note) = source_note.map(str::trim).filter(|value| !value.is_empty()) {
            result
                .messages
                .push(format!("Nucleotide source chain: {}", source_note));
        }
        Ok(seq_id)
    }

    pub(super) fn apply_internal(
        &mut self,
        op: Operation,
        run_id: &str,
        on_progress: &mut dyn FnMut(OperationProgress) -> bool,
    ) -> Result<OpResult, EngineError> {
        self.reconcile_lineage_nodes();
        self.reconcile_containers();
        let op_for_containers = op.clone();
        let op = match op {
            Operation::MergeContainersById {
                container_ids,
                output_prefix,
            } => Operation::MergeContainers {
                inputs: self.flatten_container_members(&container_ids)?,
                output_prefix,
            },
            Operation::LigationContainer {
                container_id,
                circularize_if_possible,
                output_id,
                protocol,
                output_prefix,
                unique,
            } => Operation::Ligation {
                inputs: self.container_members(&container_id)?,
                circularize_if_possible,
                output_id,
                protocol,
                output_prefix,
                unique,
            },
            Operation::FilterContainerByMolecularWeight {
                container_id,
                min_bp,
                max_bp,
                error,
                unique,
                output_prefix,
            } => Operation::FilterByMolecularWeight {
                inputs: self.container_members(&container_id)?,
                min_bp,
                max_bp,
                error,
                unique,
                output_prefix,
            },
            other => other,
        };
        let op_id = self.next_op_id();
        let mut parent_seq_ids: Vec<SeqId> = vec![];
        let mut result = OpResult {
            op_id,
            created_seq_ids: vec![],
            changed_seq_ids: vec![],
            warnings: vec![],
            messages: vec![],
            genome_annotation_projection: None,
            sequence_alignment: None,
        };

        match op {
            Operation::LoadFile { path, as_id } => {
                let mut dna = GENtleApp::load_from_file(&path).map_err(|e| EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!("Could not load sequence file '{path}': {e}"),
                })?;
                Self::prepare_sequence(&mut dna);

                let base = as_id.unwrap_or_else(|| Self::derive_seq_id(&path));
                let seq_id = self.unique_seq_id(&base);
                self.state.sequences.insert(seq_id.clone(), dna);
                self.add_lineage_node(
                    &seq_id,
                    Self::classify_import_origin(
                        &path,
                        self.state
                            .sequences
                            .get(&seq_id)
                            .expect("sequence just inserted"),
                    ),
                    Some(&result.op_id),
                );
                let imported_anchor = self
                    .state
                    .sequences
                    .get(&seq_id)
                    .and_then(|loaded| Self::infer_imported_genbank_anchor(&path, loaded));
                if let Some(anchor) = imported_anchor {
                    let mut anchor_verified: Option<bool> = None;
                    if let Some(loaded) = self.state.sequences.get(&seq_id) {
                        match Self::verify_anchor_sequence_against_catalog(
                            loaded,
                            &anchor,
                            DEFAULT_GENOME_CATALOG_PATH,
                            None,
                        ) {
                            Ok(is_match) => {
                                anchor_verified = Some(is_match);
                                if is_match {
                                    result.messages.push(format!(
                                        "Verified imported GenBank anchor '{}' against catalog '{}' ({}:{}-{})",
                                        seq_id,
                                        DEFAULT_GENOME_CATALOG_PATH,
                                        anchor.genome_id,
                                        anchor.chromosome,
                                        anchor.start_1based
                                    ));
                                } else {
                                    result.warnings.push(format!(
                                        "Imported GenBank anchor '{}' does not match catalog sequence at {}:{}:{}-{} (catalog='{}')",
                                        seq_id,
                                        anchor.genome_id,
                                        anchor.chromosome,
                                        anchor.start_1based,
                                        anchor.end_1based,
                                        DEFAULT_GENOME_CATALOG_PATH
                                    ));
                                }
                            }
                            Err(err) => {
                                result.warnings.push(format!(
                                    "Could not verify imported GenBank anchor '{}' against catalog '{}': {}",
                                    seq_id, DEFAULT_GENOME_CATALOG_PATH, err
                                ));
                            }
                        }
                    }
                    self.append_genome_extraction_provenance(GenomeExtractionProvenance {
                        seq_id: seq_id.clone(),
                        recorded_at_unix_ms: Self::now_unix_ms(),
                        operation: "LoadFileGenBankRegion".to_string(),
                        genome_id: anchor.genome_id.clone(),
                        // Imported GenBank files are sequence sources, not catalog JSON.
                        // Keep the default catalog path so later anchor-extension flows
                        // resolve against real genome catalogs instead of the .gb file.
                        catalog_path: DEFAULT_GENOME_CATALOG_PATH.to_string(),
                        cache_dir: None,
                        chromosome: Some(anchor.chromosome.clone()),
                        start_1based: Some(anchor.start_1based),
                        end_1based: Some(anchor.end_1based),
                        gene_query: None,
                        occurrence: None,
                        gene_id: None,
                        gene_name: None,
                        strand: None,
                        anchor_strand: anchor.strand,
                        anchor_verified,
                        sequence_source_type: Some("genbank_file".to_string()),
                        annotation_source_type: Some("genbank_file".to_string()),
                        sequence_source: Some(path.clone()),
                        annotation_source: Some(path.clone()),
                        sequence_sha1: None,
                        annotation_sha1: None,
                    });
                    let strand = anchor.strand.unwrap_or('+');
                    let verification_label = match anchor_verified {
                        Some(true) => "verified",
                        Some(false) => "unverified",
                        None => "verification n/a",
                    };
                    result.messages.push(format!(
                        "Detected GenBank genome anchor for '{}': {}:{}-{} ({}, strand {}, {})",
                        seq_id,
                        anchor.chromosome,
                        anchor.start_1based,
                        anchor.end_1based,
                        anchor.genome_id,
                        strand,
                        verification_label
                    ));
                }
                result.created_seq_ids.push(seq_id.clone());
                result
                    .messages
                    .push(format!("Loaded '{path}' as '{seq_id}'"));
            }
            Operation::SaveFile {
                seq_id,
                path,
                format,
            } => {
                let dna = self
                    .state
                    .sequences
                    .get(&seq_id)
                    .ok_or_else(|| EngineError {
                        code: ErrorCode::NotFound,
                        message: format!("Sequence '{seq_id}' not found"),
                    })?;

                match format {
                    ExportFormat::GenBank => {
                        dna.write_genbank_file(&path).map_err(|e| EngineError {
                            code: ErrorCode::Io,
                            message: format!("Could not write GenBank file '{path}': {e}"),
                        })?;
                    }
                    ExportFormat::Fasta => Self::save_as_fasta(&seq_id, dna, &path)?,
                }

                result.changed_seq_ids.push(seq_id.clone());
                result
                    .messages
                    .push(format!("Wrote '{seq_id}' to '{path}'"));
            }
            Operation::RenderSequenceSvg { seq_id, mode, path } => {
                let dna = self
                    .state
                    .sequences
                    .get(&seq_id)
                    .ok_or_else(|| EngineError {
                        code: ErrorCode::NotFound,
                        message: format!("Sequence '{seq_id}' not found"),
                    })?;
                let svg = match mode {
                    RenderSvgMode::Linear => export_linear_svg(dna, &self.state.display),
                    RenderSvgMode::Circular => export_circular_svg(dna, &self.state.display),
                };
                std::fs::write(&path, svg).map_err(|e| EngineError {
                    code: ErrorCode::Io,
                    message: format!("Could not write SVG output '{path}': {e}"),
                })?;
                result.messages.push(format!(
                    "Wrote {:?} SVG for '{}' to '{}'",
                    mode, seq_id, path
                ));
            }
            Operation::RenderDotplotSvg {
                seq_id,
                dotplot_id,
                path,
                flex_track_id,
                display_density_threshold,
                display_intensity_gain,
            } => {
                self.render_dotplot_svg_to_path(
                    &seq_id,
                    &dotplot_id,
                    &path,
                    flex_track_id.as_deref(),
                    display_density_threshold,
                    display_intensity_gain,
                )?;
                result.messages.push(format!(
                    "Wrote dotplot SVG for '{}' dotplot='{}' to '{}'",
                    seq_id, dotplot_id, path
                ));
            }
            Operation::RenderFeatureExpertSvg {
                seq_id,
                target,
                path,
            } => {
                self.render_feature_expert_svg_to_path(&seq_id, &target, &path)?;
                result.messages.push(format!(
                    "Wrote feature expert SVG for '{}' target={} to '{}'",
                    seq_id,
                    target.describe(),
                    path
                ));
            }
            Operation::RenderIsoformArchitectureSvg {
                seq_id,
                panel_id,
                path,
            } => {
                self.render_isoform_architecture_svg_to_path(&seq_id, &panel_id, &path)?;
                result.messages.push(format!(
                    "Wrote isoform architecture SVG for '{}' panel='{}' to '{}'",
                    seq_id, panel_id, path
                ));
            }
            Operation::RenderRnaStructureSvg { seq_id, path } => {
                let report = self.render_rna_structure_svg_to_path(&seq_id, &path)?;
                result.messages.push(format!(
                    "Wrote RNA structure SVG for '{}' to '{}' using {}",
                    seq_id, path, report.tool
                ));
            }
            Operation::RenderLineageSvg { path } => {
                let svg = export_lineage_svg(&self.state);
                std::fs::write(&path, svg).map_err(|e| EngineError {
                    code: ErrorCode::Io,
                    message: format!("Could not write SVG output '{path}': {e}"),
                })?;
                result
                    .messages
                    .push(format!("Wrote lineage SVG to '{}'", path));
            }
            Operation::RenderProtocolCartoonSvg { protocol, path } => {
                let svg = crate::protocol_cartoon::render_protocol_cartoon_svg(&protocol);
                std::fs::write(&path, svg).map_err(|e| EngineError {
                    code: ErrorCode::Io,
                    message: format!("Could not write SVG output '{path}': {e}"),
                })?;
                result.messages.push(format!(
                    "Wrote protocol cartoon SVG for '{}' to '{}'",
                    protocol.id(),
                    path
                ));
            }
            Operation::RenderProtocolCartoonTemplateSvg {
                template_path,
                path,
            } => {
                let template_json =
                    std::fs::read_to_string(&template_path).map_err(|e| EngineError {
                        code: ErrorCode::Io,
                        message: format!(
                            "Could not read protocol cartoon template '{}': {e}",
                            template_path
                        ),
                    })?;
                let template: crate::protocol_cartoon::ProtocolCartoonTemplate =
                    serde_json::from_str(&template_json).map_err(|e| EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "Could not parse protocol cartoon template JSON from '{}': {e}",
                            template_path
                        ),
                    })?;
                let spec = crate::protocol_cartoon::resolve_protocol_cartoon_template(&template)
                    .map_err(|e| EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "Invalid protocol cartoon template '{}': {}",
                            template_path, e
                        ),
                    })?;
                let svg = crate::protocol_cartoon::render_protocol_cartoon_spec_svg(&spec);
                std::fs::write(&path, svg).map_err(|e| EngineError {
                    code: ErrorCode::Io,
                    message: format!("Could not write SVG output '{path}': {e}"),
                })?;
                result.messages.push(format!(
                    "Wrote protocol cartoon template SVG for '{}' from '{}' to '{}'",
                    spec.id, template_path, path
                ));
            }
            Operation::ValidateProtocolCartoonTemplate { template_path } => {
                let template_json =
                    std::fs::read_to_string(&template_path).map_err(|e| EngineError {
                        code: ErrorCode::Io,
                        message: format!(
                            "Could not read protocol cartoon template '{}': {e}",
                            template_path
                        ),
                    })?;
                let template: crate::protocol_cartoon::ProtocolCartoonTemplate =
                    serde_json::from_str(&template_json).map_err(|e| EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "Could not parse protocol cartoon template JSON from '{}': {e}",
                            template_path
                        ),
                    })?;
                let spec = crate::protocol_cartoon::resolve_protocol_cartoon_template(&template)
                    .map_err(|e| EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "Invalid protocol cartoon template '{}': {}",
                            template_path, e
                        ),
                    })?;
                result.messages.push(format!(
                    "Validated protocol cartoon template '{}' from '{}' (events={})",
                    spec.id,
                    template_path,
                    spec.events.len()
                ));
            }
            Operation::RenderProtocolCartoonTemplateWithBindingsSvg {
                template_path,
                bindings_path,
                path,
            } => {
                let template_json =
                    std::fs::read_to_string(&template_path).map_err(|e| EngineError {
                        code: ErrorCode::Io,
                        message: format!(
                            "Could not read protocol cartoon template '{}': {e}",
                            template_path
                        ),
                    })?;
                let template: crate::protocol_cartoon::ProtocolCartoonTemplate =
                    serde_json::from_str(&template_json).map_err(|e| EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "Could not parse protocol cartoon template JSON from '{}': {e}",
                            template_path
                        ),
                    })?;
                let bindings_json =
                    std::fs::read_to_string(&bindings_path).map_err(|e| EngineError {
                        code: ErrorCode::Io,
                        message: format!(
                            "Could not read protocol cartoon bindings '{}': {e}",
                            bindings_path
                        ),
                    })?;
                let bindings: crate::protocol_cartoon::ProtocolCartoonTemplateBindings =
                    serde_json::from_str(&bindings_json).map_err(|e| EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "Could not parse protocol cartoon bindings JSON from '{}': {e}",
                            bindings_path
                        ),
                    })?;
                let spec =
                    crate::protocol_cartoon::resolve_protocol_cartoon_template_with_bindings(
                        &template, &bindings,
                    )
                    .map_err(|e| EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "Invalid protocol cartoon template/bindings ('{}', '{}'): {}",
                            template_path, bindings_path, e
                        ),
                    })?;
                let svg = crate::protocol_cartoon::render_protocol_cartoon_spec_svg(&spec);
                std::fs::write(&path, svg).map_err(|e| EngineError {
                    code: ErrorCode::Io,
                    message: format!("Could not write SVG output '{path}': {e}"),
                })?;
                result.messages.push(format!(
                    "Wrote protocol cartoon template SVG for '{}' from '{}' with bindings '{}' to '{}'",
                    spec.id, template_path, bindings_path, path
                ));
            }
            Operation::ExportProtocolCartoonTemplateJson { protocol, path } => {
                let template =
                    crate::protocol_cartoon::protocol_cartoon_template_for_kind(&protocol);
                let json = serde_json::to_string_pretty(&template).map_err(|e| EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!(
                        "Could not serialize protocol cartoon template '{}' as JSON: {e}",
                        protocol.id()
                    ),
                })?;
                std::fs::write(&path, format!("{json}\n")).map_err(|e| EngineError {
                    code: ErrorCode::Io,
                    message: format!("Could not write JSON output '{path}': {e}"),
                })?;
                result.messages.push(format!(
                    "Exported protocol cartoon template '{}' to '{}'",
                    protocol.id(),
                    path
                ));
            }
            Operation::ApplyGibsonAssemblyPlan { plan_json } => {
                let plan: GibsonAssemblyPlan =
                    serde_json::from_str(&plan_json).map_err(|err| EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!("Could not parse Gibson assembly plan JSON: {err}"),
                    })?;
                let execution = derive_gibson_execution_plan(self, &plan)?;
                parent_seq_ids = execution.parent_seq_ids.clone();
                result.warnings.extend(execution.preview.warnings.clone());
                for output in execution.outputs {
                    let seq_id = self.unique_seq_id(&output.base_seq_id);
                    let mut dna = output.dna;
                    Self::prepare_sequence(&mut dna);
                    let length_bp = dna.len();
                    let topology_label = if dna.is_circular() {
                        "circular"
                    } else {
                        "linear"
                    };
                    self.state.sequences.insert(seq_id.clone(), dna);
                    self.add_lineage_node(&seq_id, SequenceOrigin::Derived, Some(&result.op_id));
                    result.created_seq_ids.push(seq_id.clone());
                    result.messages.push(format!(
                        "Created Gibson {} '{}' ({} bp, {})",
                        output.kind.replace('_', " "),
                        seq_id,
                        length_bp,
                        topology_label
                    ));
                }
                let plan_label = execution.preview.plan_id.trim();
                let insert_summary = if execution.preview.fragments.is_empty() {
                    "-".to_string()
                } else {
                    execution
                        .preview
                        .fragments
                        .iter()
                        .map(|fragment| fragment.seq_id.clone())
                        .collect::<Vec<_>>()
                        .join(", ")
                };
                result.messages.push(format!(
                    "Applied Gibson assembly plan '{}' from destination '{}' and inserts '{}'",
                    if plan_label.is_empty() {
                        execution.preview.title.as_str()
                    } else {
                        plan_label
                    },
                    execution.preview.destination.seq_id,
                    insert_summary
                ));
            }
            Operation::CreateArrangementSerial {
                container_ids,
                arrangement_id,
                name,
                ladders,
            } => {
                let created_id = self.add_serial_arrangement(
                    &container_ids,
                    arrangement_id,
                    name,
                    ladders,
                    Some(&result.op_id),
                )?;
                result.messages.push(format!(
                    "Created serial arrangement '{}' with {} lane container(s)",
                    created_id,
                    self.state
                        .container_state
                        .arrangements
                        .get(&created_id)
                        .map(|a| a.lane_container_ids.len())
                        .unwrap_or(0)
                ));
            }
            Operation::RenderPoolGelSvg {
                inputs,
                path,
                ladders,
                container_ids,
                arrangement_id,
            } => {
                let mut ladder_names = ladders
                    .unwrap_or_default()
                    .into_iter()
                    .map(|v| v.trim().to_string())
                    .filter(|v| !v.is_empty())
                    .collect::<Vec<_>>();
                let samples: Vec<GelSampleInput> = if let Some(arrangement_id) =
                    arrangement_id.as_deref().map(str::trim)
                {
                    if arrangement_id.is_empty() {
                        return Err(EngineError {
                            code: ErrorCode::InvalidInput,
                            message: "arrangement_id cannot be empty".to_string(),
                        });
                    }
                    let (arrangement_samples, arrangement_ladders) =
                        self.gel_samples_from_arrangement(arrangement_id)?;
                    if ladder_names.is_empty() {
                        ladder_names = arrangement_ladders;
                    }
                    if !inputs.is_empty() {
                        result.warnings.push(
                            "RenderPoolGelSvg ignored 'inputs' because arrangement_id was provided"
                                .to_string(),
                        );
                    }
                    if container_ids.as_ref().is_some_and(|ids| !ids.is_empty()) {
                        result.warnings.push(
                                "RenderPoolGelSvg ignored 'container_ids' because arrangement_id was provided"
                                    .to_string(),
                            );
                    }
                    arrangement_samples
                } else if let Some(container_ids) = container_ids {
                    if container_ids.is_empty() {
                        return Err(EngineError {
                            code: ErrorCode::InvalidInput,
                            message: "container_ids was provided but empty".to_string(),
                        });
                    }
                    if !inputs.is_empty() {
                        result.warnings.push(
                            "RenderPoolGelSvg ignored 'inputs' because container_ids were provided"
                                .to_string(),
                        );
                    }
                    self.gel_samples_from_container_ids(&container_ids)?
                } else {
                    if inputs.is_empty() {
                        return Err(EngineError {
                                code: ErrorCode::InvalidInput,
                                message: "RenderPoolGelSvg requires either inputs, container_ids, or arrangement_id"
                                    .to_string(),
                            });
                    }
                    let mut members: Vec<(String, usize)> = Vec::with_capacity(inputs.len());
                    for seq_id in &inputs {
                        let dna = self
                            .state
                            .sequences
                            .get(seq_id)
                            .ok_or_else(|| EngineError {
                                code: ErrorCode::NotFound,
                                message: format!("Sequence '{seq_id}' not found"),
                            })?;
                        members.push((seq_id.clone(), dna.len()));
                    }
                    vec![GelSampleInput {
                        name: format!("Input tube (n={})", members.len()),
                        members,
                    }]
                };
                let layout =
                    build_serial_gel_layout(&samples, &ladder_names).map_err(|e| EngineError {
                        code: ErrorCode::InvalidInput,
                        message: e,
                    })?;
                let svg = export_pool_gel_svg(&layout);
                std::fs::write(&path, svg).map_err(|e| EngineError {
                    code: ErrorCode::Io,
                    message: format!("Could not write SVG output '{path}': {e}"),
                })?;
                let ladders_used = if layout.selected_ladders.is_empty() {
                    "auto".to_string()
                } else {
                    layout.selected_ladders.join(" + ")
                };
                result.messages.push(format!(
                    "Wrote serial gel SVG for {} sample lane(s), {} sequence(s) to '{}' (ladders: {})",
                    layout.sample_count,
                    layout.pool_member_count,
                    path,
                    ladders_used
                ));
            }
            Operation::ExportDnaLadders { path, name_filter } => {
                let report = Self::export_dna_ladders(&path, name_filter.as_deref())?;
                let filter_text = name_filter
                    .as_deref()
                    .map(str::trim)
                    .filter(|v| !v.is_empty())
                    .unwrap_or("-");
                result.messages.push(format!(
                    "Wrote DNA ladders catalog ({} ladder(s), filter={}) to '{}'",
                    report.ladder_count, filter_text, report.path
                ));
            }
            Operation::ExportRnaLadders { path, name_filter } => {
                let report = Self::export_rna_ladders(&path, name_filter.as_deref())?;
                let filter_text = name_filter
                    .as_deref()
                    .map(str::trim)
                    .filter(|v| !v.is_empty())
                    .unwrap_or("-");
                result.messages.push(format!(
                    "Wrote RNA ladders catalog ({} ladder(s), filter={}) to '{}'",
                    report.ladder_count, filter_text, report.path
                ));
            }
            Operation::ExportPool {
                inputs,
                path,
                pool_id,
                human_id,
            } => {
                let (pool_id, human_id, count) =
                    self.export_pool_file(&inputs, &path, pool_id, human_id)?;
                result.messages.push(format!(
                    "Wrote pool export '{}' ({}) with {} members to '{}'",
                    pool_id, human_id, count, path
                ));
            }
            Operation::ExportProcessRunBundle { path, run_id } => {
                let bundle = self.export_process_run_bundle_file(&path, run_id.as_deref())?;
                let run_scope = bundle
                    .run_id_filter
                    .as_deref()
                    .map(str::trim)
                    .filter(|value| !value.is_empty())
                    .unwrap_or("all");
                result.messages.push(format!(
                    "Wrote process run bundle '{}' with {} record(s) (run_id={}) to '{}'",
                    bundle.schema, bundle.selected_record_count, run_scope, path
                ));
            }
            Operation::PrepareGenome {
                genome_id,
                catalog_path,
                cache_dir,
                timeout_seconds,
            } => {
                let report = Self::prepare_reference_genome_once(
                    &genome_id,
                    catalog_path.as_deref(),
                    cache_dir.as_deref(),
                    timeout_seconds,
                    &mut |p| on_progress(OperationProgress::GenomePrepare(p)),
                )?;
                result.messages.push(Self::format_prepare_genome_message(
                    &genome_id,
                    cache_dir.as_deref(),
                    &report,
                ));
                if !report.warnings.is_empty() {
                    result.warnings.extend(report.warnings.clone());
                }
            }
            Operation::ExtractGenomeRegion {
                genome_id,
                chromosome,
                start_1based,
                end_1based,
                output_id,
                annotation_scope,
                max_annotation_features,
                include_genomic_annotation,
                catalog_path,
                cache_dir,
            } => {
                let catalog_path =
                    catalog_path.unwrap_or_else(|| DEFAULT_GENOME_CATALOG_PATH.to_string());
                let catalog =
                    GenomeCatalog::from_json_file(&catalog_path).map_err(|e| EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!("Could not open genome catalog '{catalog_path}': {e}"),
                    })?;
                let sequence = catalog
                    .get_sequence_region_with_cache(
                        &genome_id,
                        &chromosome,
                        start_1based,
                        end_1based,
                        cache_dir.as_deref(),
                    )
                    .map_err(|e| EngineError {
                        code: ErrorCode::NotFound,
                        message: format!(
                            "Could not load genome region {}:{}-{} from '{}': {}",
                            chromosome, start_1based, end_1based, genome_id, e
                        ),
                    })?;
                let default_id = format!(
                    "{}_{}_{}_{}",
                    Self::normalize_id_token(&genome_id),
                    Self::normalize_id_token(&chromosome),
                    start_1based,
                    end_1based
                );
                let base = output_id.unwrap_or(default_id);
                let seq_id = self.import_genome_slice_sequence(&mut result, sequence, base)?;
                let requested_scope = Self::resolve_extract_region_annotation_scope(
                    annotation_scope,
                    include_genomic_annotation,
                );
                let effective_cap =
                    max_annotation_features
                        .filter(|value| *value > 0)
                        .or_else(|| {
                            matches!(requested_scope, GenomeAnnotationScope::Full)
                                .then_some(DEFAULT_EXTRACT_REGION_ANNOTATION_FEATURE_CAP)
                        });
                let mut genes: Vec<GenomeGeneRecord> = vec![];
                let mut transcripts: Vec<GenomeTranscriptRecord> = vec![];
                if !matches!(requested_scope, GenomeAnnotationScope::None) {
                    match catalog.list_gene_regions(&genome_id, cache_dir.as_deref()) {
                        Ok(records) => {
                            genes = records
                                .into_iter()
                                .filter(|record| {
                                    Self::genome_chromosome_matches(&record.chromosome, &chromosome)
                                        && record.end_1based >= start_1based
                                        && record.start_1based <= end_1based
                                })
                                .collect();
                        }
                        Err(e) => result.warnings.push(format!(
                            "Could not inspect gene records for extracted region '{}': {}",
                            seq_id, e
                        )),
                    }
                    match catalog.list_gene_transcript_records(
                        &genome_id,
                        &chromosome,
                        start_1based,
                        end_1based,
                        None,
                        None,
                        cache_dir.as_deref(),
                    ) {
                        Ok(records) => {
                            transcripts = records;
                        }
                        Err(e) => result.warnings.push(format!(
                            "Could not inspect transcript/exon annotation for extracted region '{}': {}",
                            seq_id, e
                        )),
                    }
                }
                let mut effective_scope = requested_scope;
                let mut fallback_reason: Option<String> = None;
                let mut projection = Self::build_extract_region_annotation_projection(
                    &genes,
                    &transcripts,
                    start_1based,
                    end_1based,
                    requested_scope,
                );
                let candidate_before_fallback = projection.feature_count();
                if let Some(cap) = effective_cap {
                    if projection.feature_count() > cap {
                        if matches!(requested_scope, GenomeAnnotationScope::Full) {
                            let core_projection = Self::build_extract_region_annotation_projection(
                                &genes,
                                &transcripts,
                                start_1based,
                                end_1based,
                                GenomeAnnotationScope::Core,
                            );
                            if core_projection.feature_count() <= cap {
                                effective_scope = GenomeAnnotationScope::Core;
                                projection = core_projection;
                                fallback_reason = Some(format!(
                                    "Projected full annotation would attach {} feature(s), exceeding max_annotation_features={cap}; fell back to core projection. Re-run with annotation_scope=full --max-annotation-features 0 to force full transfer.",
                                    candidate_before_fallback
                                ));
                            } else {
                                effective_scope = GenomeAnnotationScope::None;
                                projection = ExtractRegionAnnotationProjectionBatch::default();
                                fallback_reason = Some(format!(
                                    "Projected annotation exceeded max_annotation_features={cap} even after core fallback; annotation transfer was disabled for this extraction. Re-run with --max-annotation-features 0 for unrestricted transfer."
                                ));
                            }
                        } else {
                            effective_scope = GenomeAnnotationScope::None;
                            projection = ExtractRegionAnnotationProjectionBatch::default();
                            fallback_reason = Some(format!(
                                "Projected annotation exceeded max_annotation_features={cap}; annotation transfer was disabled for this extraction."
                            ));
                        }
                    }
                }

                let attached_feature_count = projection.feature_count();
                let genes_attached = projection.gene_count;
                let transcripts_attached = projection.transcript_count;
                let exons_attached = projection.exon_count;
                let cds_attached = projection.cds_count;
                if attached_feature_count > 0 {
                    if let Some(dna) = self.state.sequences.get_mut(&seq_id) {
                        dna.features_mut().extend(projection.features);
                        Self::prepare_sequence(dna);
                    }
                    result.messages.push(format!(
                        "Attached {} genomic annotation feature(s) to extracted region '{}' (scope={}, genes={}, transcripts={}, exons={}, cds={})",
                        attached_feature_count,
                        seq_id,
                        effective_scope.as_str(),
                        genes_attached,
                        transcripts_attached,
                        exons_attached,
                        cds_attached
                    ));
                } else if matches!(effective_scope, GenomeAnnotationScope::None) {
                    result.messages.push(format!(
                        "Extracted region '{}' without projected genomic annotation (scope={})",
                        seq_id,
                        effective_scope.as_str()
                    ));
                } else {
                    result.messages.push(format!(
                        "No overlapping genomic annotation features were found for extracted region '{}'",
                        seq_id
                    ));
                }
                if let Some(reason) = fallback_reason.as_ref() {
                    result.warnings.push(reason.clone());
                }
                result.genome_annotation_projection = Some(GenomeAnnotationProjectionTelemetry {
                    requested_scope: requested_scope.as_str().to_string(),
                    effective_scope: effective_scope.as_str().to_string(),
                    max_features_cap: effective_cap,
                    candidate_feature_count: candidate_before_fallback,
                    attached_feature_count,
                    dropped_feature_count: candidate_before_fallback
                        .saturating_sub(attached_feature_count),
                    genes_attached,
                    transcripts_attached,
                    exons_attached,
                    cds_attached,
                    fallback_applied: fallback_reason.is_some(),
                    fallback_reason,
                });
                self.maybe_attach_known_helper_mcs_annotation(&seq_id, &genome_id, &mut result);
                let source_plan = catalog.source_plan(&genome_id, cache_dir.as_deref()).ok();
                let inspection = catalog
                    .inspect_prepared_genome(&genome_id, cache_dir.as_deref())
                    .ok()
                    .flatten();
                let (
                    sequence_source_type,
                    annotation_source_type,
                    sequence_source,
                    annotation_source,
                    sequence_sha1,
                    annotation_sha1,
                ) = Self::genome_source_snapshot(source_plan.as_ref(), inspection.as_ref());
                self.append_genome_extraction_provenance(GenomeExtractionProvenance {
                    seq_id: seq_id.clone(),
                    recorded_at_unix_ms: Self::now_unix_ms(),
                    operation: "ExtractGenomeRegion".to_string(),
                    genome_id: genome_id.clone(),
                    catalog_path: catalog_path.clone(),
                    cache_dir: cache_dir.clone(),
                    chromosome: Some(chromosome.clone()),
                    start_1based: Some(start_1based),
                    end_1based: Some(end_1based),
                    gene_query: None,
                    occurrence: None,
                    gene_id: None,
                    gene_name: None,
                    strand: None,
                    anchor_strand: Some('+'),
                    anchor_verified: Some(true),
                    sequence_source_type,
                    annotation_source_type,
                    sequence_source,
                    annotation_source,
                    sequence_sha1,
                    annotation_sha1,
                });
                result.messages.push(format!(
                    "Extracted genome region {}:{}-{} from '{}' as '{}'",
                    chromosome, start_1based, end_1based, genome_id, seq_id
                ));
            }
            Operation::ExtractGenomeGene {
                genome_id,
                gene_query,
                occurrence,
                output_id,
                annotation_scope,
                max_annotation_features,
                include_genomic_annotation,
                catalog_path,
                cache_dir,
            } => {
                let query = gene_query.trim();
                if query.is_empty() {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: "Gene query cannot be empty".to_string(),
                    });
                }
                let catalog_path =
                    catalog_path.unwrap_or_else(|| DEFAULT_GENOME_CATALOG_PATH.to_string());
                let catalog =
                    GenomeCatalog::from_json_file(&catalog_path).map_err(|e| EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!("Could not open genome catalog '{catalog_path}': {e}"),
                    })?;
                let genes = catalog
                    .list_gene_regions(&genome_id, cache_dir.as_deref())
                    .map_err(|e| EngineError {
                        code: ErrorCode::NotFound,
                        message: format!(
                            "Could not load gene index for genome '{}': {}",
                            genome_id, e
                        ),
                    })?;
                let mut exact_matches: Vec<&GenomeGeneRecord> = genes
                    .iter()
                    .filter(|record| Self::genome_gene_matches_exact(record, query))
                    .collect();
                let used_fuzzy = if exact_matches.is_empty() {
                    let query_lower = query.to_ascii_lowercase();
                    exact_matches = genes
                        .iter()
                        .filter(|record| Self::genome_gene_matches_contains(record, &query_lower))
                        .collect();
                    true
                } else {
                    false
                };
                if exact_matches.is_empty() {
                    return Err(EngineError {
                        code: ErrorCode::NotFound,
                        message: format!("No genes in '{}' match query '{}'", genome_id, query),
                    });
                }
                let requested_occurrence = occurrence;
                let occurrence = requested_occurrence.unwrap_or(1);
                if occurrence == 0 {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: "Gene occurrence must be >= 1".to_string(),
                    });
                }
                if exact_matches.len() > 1 && requested_occurrence.is_none() {
                    result.warnings.push(format!(
                        "Gene query '{}' matched {} records in '{}'; using first match (set occurrence for another match).",
                        query,
                        exact_matches.len(),
                        genome_id
                    ));
                }
                let Some(selected_gene) = exact_matches.get(occurrence - 1) else {
                    return Err(EngineError {
                        code: ErrorCode::NotFound,
                        message: format!(
                            "Gene query '{}' matched {} records, but occurrence {} was requested",
                            query,
                            exact_matches.len(),
                            occurrence
                        ),
                    });
                };
                let sequence = catalog
                    .get_sequence_region_with_cache(
                        &genome_id,
                        &selected_gene.chromosome,
                        selected_gene.start_1based,
                        selected_gene.end_1based,
                        cache_dir.as_deref(),
                    )
                    .map_err(|e| EngineError {
                        code: ErrorCode::NotFound,
                        message: format!(
                            "Could not load gene region {}:{}-{} from '{}': {}",
                            selected_gene.chromosome,
                            selected_gene.start_1based,
                            selected_gene.end_1based,
                            genome_id,
                            e
                        ),
                    })?;
                let gene_label = selected_gene
                    .gene_name
                    .as_ref()
                    .or(selected_gene.gene_id.as_ref())
                    .cloned()
                    .unwrap_or_else(|| query.to_string());
                let default_id = format!(
                    "{}_{}_{}_{}",
                    Self::normalize_id_token(&genome_id),
                    Self::normalize_id_token(&gene_label),
                    selected_gene.start_1based,
                    selected_gene.end_1based
                );
                let base = output_id.unwrap_or(default_id);
                let seq_id = self.import_genome_slice_sequence(&mut result, sequence, base)?;
                let preferred_name = self
                    .state
                    .sequences
                    .get(&seq_id)
                    .and_then(|dna| {
                        dna.name()
                            .as_deref()
                            .map(str::trim)
                            .filter(|v| !v.is_empty())
                            .map(str::to_string)
                    })
                    .unwrap_or_else(|| seq_id.clone());
                if preferred_name.eq_ignore_ascii_case("<unnamed sequence>")
                    || preferred_name.eq_ignore_ascii_case("<no name>")
                {
                    if let Some(dna) = self.state.sequences.get_mut(&seq_id) {
                        dna.set_name(seq_id.clone());
                        Self::prepare_sequence(dna);
                    }
                }

                let mut transcript_records = match catalog.list_gene_transcript_records(
                    &genome_id,
                    &selected_gene.chromosome,
                    selected_gene.start_1based,
                    selected_gene.end_1based,
                    selected_gene.gene_id.as_deref(),
                    selected_gene.gene_name.as_deref(),
                    cache_dir.as_deref(),
                ) {
                    Ok(records) => records,
                    Err(e) => {
                        result.warnings.push(format!(
                            "Could not inspect transcript/exon annotation for extracted gene '{}': {}",
                            seq_id, e
                        ));
                        vec![]
                    }
                };
                if transcript_records.is_empty() {
                    match catalog.list_gene_transcript_records(
                        &genome_id,
                        &selected_gene.chromosome,
                        selected_gene.start_1based,
                        selected_gene.end_1based,
                        None,
                        None,
                        cache_dir.as_deref(),
                    ) {
                        Ok(mut records) => {
                            let selected_gene_id_norm = selected_gene
                                .gene_id
                                .as_deref()
                                .map(Self::normalize_id_token)
                                .filter(|v| !v.is_empty());
                            let selected_gene_name_norm = selected_gene
                                .gene_name
                                .as_deref()
                                .map(Self::normalize_id_token)
                                .filter(|v| !v.is_empty());
                            if selected_gene_id_norm.is_some() || selected_gene_name_norm.is_some()
                            {
                                records.retain(|record| {
                                    let transcript_gene_id_norm = record
                                        .gene_id
                                        .as_deref()
                                        .map(Self::normalize_id_token)
                                        .filter(|v| !v.is_empty());
                                    let transcript_gene_name_norm = record
                                        .gene_name
                                        .as_deref()
                                        .map(Self::normalize_id_token)
                                        .filter(|v| !v.is_empty());
                                    selected_gene_id_norm
                                        .as_ref()
                                        .zip(transcript_gene_id_norm.as_ref())
                                        .map(|(left, right)| left == right)
                                        .unwrap_or(false)
                                        || selected_gene_name_norm
                                            .as_ref()
                                            .zip(transcript_gene_name_norm.as_ref())
                                            .map(|(left, right)| left == right)
                                            .unwrap_or(false)
                                });
                            }
                            if !records.is_empty() {
                                result.warnings.push(format!(
                                    "Gene-scoped transcript filter returned no records for '{}'; using overlap fallback with {} transcript candidate(s).",
                                    seq_id,
                                    records.len()
                                ));
                                transcript_records = records;
                            }
                        }
                        Err(e) => {
                            result.warnings.push(format!(
                                "Could not run transcript fallback for extracted gene '{}': {}",
                                seq_id, e
                            ));
                        }
                    }
                }

                if let Some(extracted_sequence) = self
                    .state
                    .sequences
                    .get(&seq_id)
                    .map(|dna| dna.get_forward_string())
                    && let Some(exon_projection) = Self::build_exon_concatenated_projection(
                        &extracted_sequence,
                        selected_gene.start_1based,
                        selected_gene.end_1based,
                        selected_gene.strand,
                        &transcript_records,
                        EXON_CONCAT_SPACER_BP,
                    )
                {
                    let exon_default_id = format!("{seq_id}__exons");
                    let exon_seq_id = self.unique_seq_id(&exon_default_id);
                    let mut exon_dna = DNAsequence::from_sequence(&exon_projection.sequence)
                        .map_err(|e| EngineError {
                            code: ErrorCode::Internal,
                            message: format!(
                                "Could not construct exon-concatenated sequence for '{}': {e}",
                                seq_id
                            ),
                        })?;
                    exon_dna.set_name(exon_seq_id.clone());
                    let sequence_len = exon_projection.sequence.len();
                    if sequence_len > 0 {
                        let mut transcript_qualifiers = vec![
                            ("chromosome".into(), Some(selected_gene.chromosome.clone())),
                            (
                                "genomic_start_1based".into(),
                                Some(selected_gene.start_1based.to_string()),
                            ),
                            (
                                "genomic_end_1based".into(),
                                Some(selected_gene.end_1based.to_string()),
                            ),
                            (
                                "synthetic_origin".into(),
                                Some("exon_concatenated".to_string()),
                            ),
                            (
                                "synthetic_spacer_bp".into(),
                                Some(EXON_CONCAT_SPACER_BP.to_string()),
                            ),
                        ];
                        if let Some(gene_name) = selected_gene.gene_name.as_ref() {
                            transcript_qualifiers.push(("gene".into(), Some(gene_name.clone())));
                            transcript_qualifiers.push(("label".into(), Some(gene_name.clone())));
                        }
                        if let Some(gene_id) = selected_gene.gene_id.as_ref() {
                            transcript_qualifiers.push(("gene_id".into(), Some(gene_id.clone())));
                            if selected_gene.gene_name.is_none() {
                                transcript_qualifiers.push(("label".into(), Some(gene_id.clone())));
                            }
                        }
                        if let Some(strand) = selected_gene.strand {
                            transcript_qualifiers.push(("strand".into(), Some(strand.to_string())));
                        }
                        exon_dna.features_mut().push(gb_io::seq::Feature {
                            kind: gb_io::seq::FeatureKind::from("mRNA"),
                            location: gb_io::seq::Location::simple_range(0, sequence_len as i64),
                            qualifiers: transcript_qualifiers,
                        });
                    }
                    for (index, block) in exon_projection.blocks.iter().enumerate() {
                        if block.local_end_0based_exclusive <= block.local_start_0based {
                            continue;
                        }
                        let mut qualifiers = vec![
                            (
                                "exon_number".into(),
                                Some((index.saturating_add(1)).to_string()),
                            ),
                            ("chromosome".into(), Some(selected_gene.chromosome.clone())),
                            (
                                "genomic_start_1based".into(),
                                Some(block.genomic_start_1based.to_string()),
                            ),
                            (
                                "genomic_end_1based".into(),
                                Some(block.genomic_end_1based.to_string()),
                            ),
                            (
                                "synthetic_origin".into(),
                                Some("exon_concatenated".to_string()),
                            ),
                        ];
                        if let Some(gene_name) = selected_gene.gene_name.as_ref() {
                            qualifiers.push(("gene".into(), Some(gene_name.clone())));
                            qualifiers.push(("label".into(), Some(gene_name.clone())));
                        }
                        if let Some(gene_id) = selected_gene.gene_id.as_ref() {
                            qualifiers.push(("gene_id".into(), Some(gene_id.clone())));
                            if selected_gene.gene_name.is_none() {
                                qualifiers.push(("label".into(), Some(gene_id.clone())));
                            }
                        }
                        if let Some(strand) = selected_gene.strand {
                            qualifiers.push(("strand".into(), Some(strand.to_string())));
                        }
                        exon_dna.features_mut().push(gb_io::seq::Feature {
                            kind: gb_io::seq::FeatureKind::from("exon"),
                            location: gb_io::seq::Location::simple_range(
                                block.local_start_0based as i64,
                                block.local_end_0based_exclusive as i64,
                            ),
                            qualifiers,
                        });
                    }
                    Self::prepare_sequence(&mut exon_dna);
                    self.state.sequences.insert(exon_seq_id.clone(), exon_dna);
                    self.add_lineage_node(
                        &exon_seq_id,
                        SequenceOrigin::Derived,
                        Some(&result.op_id),
                    );
                    self.add_lineage_edges(
                        &[seq_id.clone()],
                        std::slice::from_ref(&exon_seq_id),
                        &result.op_id,
                        run_id,
                    );
                    result.created_seq_ids.push(exon_seq_id.clone());
                    result.messages.push(format!(
                        "Created exon-concatenated synthetic sequence '{}' from '{}' ({} merged exon block(s), spacer={} bp).",
                        exon_seq_id,
                        seq_id,
                        exon_projection.blocks.len(),
                        EXON_CONCAT_SPACER_BP
                    ));
                }

                let requested_scope = Self::resolve_extract_region_annotation_scope(
                    annotation_scope,
                    include_genomic_annotation,
                );
                let max_annotation_features = match max_annotation_features {
                    Some(0) | None => None,
                    Some(value) => Some(value),
                };
                let mut effective_scope = requested_scope;
                let mut fallback_reason: Option<String> = None;
                let gene_records = vec![(*selected_gene).clone()];
                let mut projection = Self::build_extract_region_annotation_projection(
                    &gene_records,
                    &transcript_records,
                    selected_gene.start_1based,
                    selected_gene.end_1based,
                    requested_scope,
                );
                let candidate_before_fallback = projection.feature_count();
                if let Some(cap) = max_annotation_features {
                    if projection.feature_count() > cap {
                        if matches!(requested_scope, GenomeAnnotationScope::Full) {
                            let core_projection = Self::build_extract_region_annotation_projection(
                                &gene_records,
                                &transcript_records,
                                selected_gene.start_1based,
                                selected_gene.end_1based,
                                GenomeAnnotationScope::Core,
                            );
                            if core_projection.feature_count() <= cap {
                                effective_scope = GenomeAnnotationScope::Core;
                                projection = core_projection;
                                fallback_reason = Some(format!(
                                    "Projected full annotation would attach {} feature(s), exceeding max_annotation_features={cap}; fell back to core projection. Re-run with annotation_scope=full --max-annotation-features 0 to force full transfer.",
                                    candidate_before_fallback
                                ));
                            } else {
                                effective_scope = GenomeAnnotationScope::None;
                                projection = ExtractRegionAnnotationProjectionBatch::default();
                                fallback_reason = Some(format!(
                                    "Projected annotation exceeded max_annotation_features={cap} even after core fallback; annotation transfer was disabled for this extraction. Re-run with --max-annotation-features 0 for unrestricted transfer."
                                ));
                            }
                        } else {
                            effective_scope = GenomeAnnotationScope::None;
                            projection = ExtractRegionAnnotationProjectionBatch::default();
                            fallback_reason = Some(format!(
                                "Projected annotation exceeded max_annotation_features={cap}; annotation transfer was disabled for this extraction."
                            ));
                        }
                    }
                }
                let attached_feature_count = projection.feature_count();
                let genes_attached = projection.gene_count;
                let transcripts_attached = projection.transcript_count;
                let exons_attached = projection.exon_count;
                let cds_attached = projection.cds_count;
                if attached_feature_count > 0 {
                    if let Some(dna) = self.state.sequences.get_mut(&seq_id) {
                        dna.features_mut().extend(projection.features);
                        Self::prepare_sequence(dna);
                    }
                    result.messages.push(format!(
                        "Attached {} genomic annotation feature(s) to extracted gene '{}' (scope={} -> {}, genes={}, transcripts={}, exons={}, cds={})",
                        attached_feature_count,
                        seq_id,
                        requested_scope.as_str(),
                        effective_scope.as_str(),
                        genes_attached,
                        transcripts_attached,
                        exons_attached,
                        cds_attached
                    ));
                } else if !matches!(requested_scope, GenomeAnnotationScope::None) {
                    result.warnings.push(format!(
                        "No overlapping genomic annotation features were attached to extracted gene '{}'",
                        seq_id
                    ));
                }
                if let Some(reason) = fallback_reason.as_ref() {
                    result.warnings.push(reason.clone());
                }
                result.genome_annotation_projection = Some(GenomeAnnotationProjectionTelemetry {
                    requested_scope: requested_scope.as_str().to_string(),
                    effective_scope: effective_scope.as_str().to_string(),
                    max_features_cap: max_annotation_features,
                    candidate_feature_count: candidate_before_fallback,
                    attached_feature_count,
                    dropped_feature_count: candidate_before_fallback
                        .saturating_sub(attached_feature_count),
                    genes_attached,
                    transcripts_attached,
                    exons_attached,
                    cds_attached,
                    fallback_applied: fallback_reason.is_some(),
                    fallback_reason,
                });
                self.maybe_attach_known_helper_mcs_annotation(&seq_id, &genome_id, &mut result);
                let source_plan = catalog.source_plan(&genome_id, cache_dir.as_deref()).ok();
                let inspection = catalog
                    .inspect_prepared_genome(&genome_id, cache_dir.as_deref())
                    .ok()
                    .flatten();
                let (
                    sequence_source_type,
                    annotation_source_type,
                    sequence_source,
                    annotation_source,
                    sequence_sha1,
                    annotation_sha1,
                ) = Self::genome_source_snapshot(source_plan.as_ref(), inspection.as_ref());
                self.append_genome_extraction_provenance(GenomeExtractionProvenance {
                    seq_id: seq_id.clone(),
                    recorded_at_unix_ms: Self::now_unix_ms(),
                    operation: "ExtractGenomeGene".to_string(),
                    genome_id: genome_id.clone(),
                    catalog_path: catalog_path.clone(),
                    cache_dir: cache_dir.clone(),
                    chromosome: Some(selected_gene.chromosome.clone()),
                    start_1based: Some(selected_gene.start_1based),
                    end_1based: Some(selected_gene.end_1based),
                    gene_query: Some(query.to_string()),
                    occurrence: Some(occurrence),
                    gene_id: selected_gene.gene_id.clone(),
                    gene_name: selected_gene.gene_name.clone(),
                    strand: selected_gene.strand,
                    anchor_strand: Some('+'),
                    anchor_verified: Some(true),
                    sequence_source_type,
                    annotation_source_type,
                    sequence_source,
                    annotation_source,
                    sequence_sha1,
                    annotation_sha1,
                });
                let match_mode = if used_fuzzy { "fuzzy" } else { "exact" };
                result.messages.push(format!(
                    "Extracted genome gene '{}' [{} match, occurrence {}] as '{}' from '{}' ({})",
                    query,
                    match_mode,
                    occurrence,
                    seq_id,
                    genome_id,
                    Self::genome_gene_display_label(selected_gene)
                ));
            }
            Operation::ExtendGenomeAnchor {
                seq_id,
                side,
                length_bp,
                output_id,
                catalog_path,
                cache_dir,
                prepared_genome_id,
            } => {
                if length_bp == 0 {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: "ExtendGenomeAnchor requires length_bp >= 1".to_string(),
                    });
                }
                if !self.state.sequences.contains_key(&seq_id) {
                    return Err(EngineError {
                        code: ErrorCode::NotFound,
                        message: format!("Sequence '{seq_id}' not found"),
                    });
                }
                let anchor = self.latest_genome_anchor_for_seq(&seq_id)?;
                if self
                    .state
                    .parameters
                    .require_verified_genome_anchor_for_extension
                    && anchor.anchor_verified != Some(true)
                {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "ExtendGenomeAnchor requires a verified genome anchor when parameter 'require_verified_genome_anchor_for_extension' is true (seq_id='{}')",
                            seq_id
                        ),
                    });
                }
                let explicit_catalog_requested = catalog_path
                    .as_deref()
                    .map(str::trim)
                    .map(|v| !v.is_empty())
                    .unwrap_or(false);
                let requested_catalog_path = catalog_path
                    .or(anchor.catalog_path.clone())
                    .map(|v| v.trim().to_string())
                    .filter(|v| !v.is_empty())
                    .unwrap_or_else(|| DEFAULT_GENOME_CATALOG_PATH.to_string());
                let mut resolved_catalog_path = requested_catalog_path.clone();
                let resolved_cache_dir = cache_dir.or(anchor.cache_dir.clone());
                let mut catalog_fallback_warning: Option<String> = None;
                let catalog = match GenomeCatalog::from_json_file(&requested_catalog_path) {
                    Ok(catalog) => catalog,
                    Err(primary_err) => {
                        let default_catalog_path = DEFAULT_GENOME_CATALOG_PATH.to_string();
                        let should_fallback_to_default = !explicit_catalog_requested
                            && requested_catalog_path != default_catalog_path;
                        if should_fallback_to_default {
                            match GenomeCatalog::from_json_file(&default_catalog_path) {
                                Ok(catalog) => {
                                    resolved_catalog_path = default_catalog_path.clone();
                                    catalog_fallback_warning = Some(format!(
                                        "Could not open genome catalog '{}' from anchor provenance ({}). Falling back to default '{}'.",
                                        requested_catalog_path, primary_err, default_catalog_path
                                    ));
                                    catalog
                                }
                                Err(default_err) => {
                                    return Err(EngineError {
                                        code: ErrorCode::InvalidInput,
                                        message: format!(
                                            "Could not open genome catalog '{}' ({}) and fallback '{}' ({})",
                                            requested_catalog_path,
                                            primary_err,
                                            default_catalog_path,
                                            default_err
                                        ),
                                    });
                                }
                            }
                        } else {
                            return Err(EngineError {
                                code: ErrorCode::InvalidInput,
                                message: format!(
                                    "Could not open genome catalog '{}': {}",
                                    requested_catalog_path, primary_err
                                ),
                            });
                        }
                    }
                };
                if let Some(warning) = catalog_fallback_warning {
                    result.warnings.push(warning);
                }

                let anchor_is_reverse = anchor.strand == Some('-');
                let (new_start_1based, new_end_1based) = match (anchor_is_reverse, side) {
                    (false, GenomeAnchorSide::FivePrime) => (
                        anchor.start_1based.saturating_sub(length_bp).max(1),
                        anchor.end_1based,
                    ),
                    (false, GenomeAnchorSide::ThreePrime) => (
                        anchor.start_1based,
                        anchor.end_1based.saturating_add(length_bp),
                    ),
                    (true, GenomeAnchorSide::FivePrime) => (
                        anchor.start_1based,
                        anchor.end_1based.saturating_add(length_bp),
                    ),
                    (true, GenomeAnchorSide::ThreePrime) => (
                        anchor.start_1based.saturating_sub(length_bp).max(1),
                        anchor.end_1based,
                    ),
                };
                let preferred_prepared = prepared_genome_id
                    .as_deref()
                    .map(str::trim)
                    .filter(|v| !v.is_empty())
                    .map(|v| v.to_string());
                let resolution_query = preferred_prepared
                    .as_deref()
                    .unwrap_or(&anchor.genome_id)
                    .to_string();
                let fallback_policy = if preferred_prepared.is_some() {
                    PreparedGenomeFallbackPolicy::Off
                } else {
                    self.genome_anchor_fallback_policy_for_extension()
                };
                let prepared_resolution = catalog
                    .resolve_prepared_genome_id_with_policy(
                        &resolution_query,
                        resolved_cache_dir.as_deref(),
                        fallback_policy,
                    )
                    .map_err(|e| EngineError {
                        code: ErrorCode::NotFound,
                        message: format!(
                            "Could not resolve prepared genome '{}' for extension: {}",
                            resolution_query, e
                        ),
                    })?;
                if let Some(warning) = prepared_resolution.fallback_warning {
                    result.warnings.push(warning);
                }
                let effective_genome_id = prepared_resolution.resolved_genome_id;
                let mut sequence = catalog
                    .get_sequence_region_with_cache(
                        &effective_genome_id,
                        &anchor.chromosome,
                        new_start_1based,
                        new_end_1based,
                        resolved_cache_dir.as_deref(),
                    )
                    .map_err(|e| EngineError {
                        code: ErrorCode::NotFound,
                        message: format!(
                            "Could not load extended genome region {}:{}-{} from '{}': {}",
                            anchor.chromosome,
                            new_start_1based,
                            new_end_1based,
                            effective_genome_id,
                            e
                        ),
                    })?;
                if anchor_is_reverse {
                    sequence = Self::reverse_complement(&sequence);
                }

                let side_token = side.as_str();
                let default_id = format!("{seq_id}_ext_{side_token}_{length_bp}");
                let base = output_id.unwrap_or(default_id);
                let mut dna = DNAsequence::from_sequence(&sequence).map_err(|e| EngineError {
                    code: ErrorCode::Internal,
                    message: format!("Could not construct DNA sequence from extended anchor: {e}"),
                })?;
                Self::prepare_sequence(&mut dna);
                let extended_seq_id = self.unique_seq_id(&base);
                dna.set_name(extended_seq_id.clone());
                self.state.sequences.insert(extended_seq_id.clone(), dna);
                self.add_lineage_node(
                    &extended_seq_id,
                    SequenceOrigin::Derived,
                    Some(&result.op_id),
                );
                result.created_seq_ids.push(extended_seq_id.clone());
                parent_seq_ids.push(seq_id.clone());

                let source_plan = catalog
                    .source_plan(&effective_genome_id, resolved_cache_dir.as_deref())
                    .ok();
                let inspection = catalog
                    .inspect_prepared_genome(&effective_genome_id, resolved_cache_dir.as_deref())
                    .ok()
                    .flatten();
                let (
                    sequence_source_type,
                    annotation_source_type,
                    sequence_source,
                    annotation_source,
                    sequence_sha1,
                    annotation_sha1,
                ) = Self::genome_source_snapshot(source_plan.as_ref(), inspection.as_ref());
                self.append_genome_extraction_provenance(GenomeExtractionProvenance {
                    seq_id: extended_seq_id.clone(),
                    recorded_at_unix_ms: Self::now_unix_ms(),
                    operation: "ExtendGenomeAnchor".to_string(),
                    genome_id: effective_genome_id.clone(),
                    catalog_path: resolved_catalog_path.clone(),
                    cache_dir: resolved_cache_dir.clone(),
                    chromosome: Some(anchor.chromosome.clone()),
                    start_1based: Some(new_start_1based),
                    end_1based: Some(new_end_1based),
                    gene_query: None,
                    occurrence: None,
                    gene_id: None,
                    gene_name: None,
                    strand: anchor.strand,
                    anchor_strand: Some(anchor.strand.unwrap_or('+')),
                    anchor_verified: Some(true),
                    sequence_source_type,
                    annotation_source_type,
                    sequence_source,
                    annotation_source,
                    sequence_sha1,
                    annotation_sha1,
                });
                let side_label = match side {
                    GenomeAnchorSide::FivePrime => "5'",
                    GenomeAnchorSide::ThreePrime => "3'",
                };
                let anchor_strand = anchor.strand.unwrap_or('+');
                result.messages.push(format!(
                    "Extended genome anchor '{}' on {} by {} bp (anchor strand {}) via '{}' => {}:{}-{} as '{}'",
                    seq_id,
                    side_label,
                    length_bp,
                    anchor_strand,
                    effective_genome_id,
                    anchor.chromosome,
                    new_start_1based,
                    new_end_1based,
                    extended_seq_id
                ));
                let lower_bound_clipped = matches!(
                    (anchor_is_reverse, side),
                    (false, GenomeAnchorSide::FivePrime) | (true, GenomeAnchorSide::ThreePrime)
                ) && new_start_1based == 1
                    && anchor.start_1based <= length_bp;
                if lower_bound_clipped {
                    result.warnings.push(format!(
                        "Requested {} bp {} extension for '{}' clipped at chromosome start position 1",
                        length_bp, side_label, seq_id
                    ));
                }
            }
            Operation::VerifyGenomeAnchor {
                seq_id,
                catalog_path,
                cache_dir,
                prepared_genome_id,
            } => {
                let anchor = self.latest_genome_anchor_for_seq(&seq_id)?;
                let dna = self
                    .state
                    .sequences
                    .get(&seq_id)
                    .ok_or_else(|| EngineError {
                        code: ErrorCode::NotFound,
                        message: format!("Sequence '{seq_id}' not found"),
                    })?
                    .clone();
                let explicit_catalog_requested = catalog_path
                    .as_deref()
                    .map(str::trim)
                    .map(|v| !v.is_empty())
                    .unwrap_or(false);
                let requested_catalog_path = catalog_path
                    .or(anchor.catalog_path.clone())
                    .map(|v| v.trim().to_string())
                    .filter(|v| !v.is_empty())
                    .unwrap_or_else(|| DEFAULT_GENOME_CATALOG_PATH.to_string());
                let mut resolved_catalog_path = requested_catalog_path.clone();
                let resolved_cache_dir = cache_dir.or(anchor.cache_dir.clone());
                let mut catalog_fallback_warning: Option<String> = None;
                let catalog = match GenomeCatalog::from_json_file(&requested_catalog_path) {
                    Ok(catalog) => catalog,
                    Err(primary_err) => {
                        let default_catalog_path = DEFAULT_GENOME_CATALOG_PATH.to_string();
                        let should_fallback_to_default = !explicit_catalog_requested
                            && requested_catalog_path != default_catalog_path;
                        if should_fallback_to_default {
                            match GenomeCatalog::from_json_file(&default_catalog_path) {
                                Ok(catalog) => {
                                    resolved_catalog_path = default_catalog_path.clone();
                                    catalog_fallback_warning = Some(format!(
                                        "Could not open genome catalog '{}' from anchor provenance ({}). Falling back to default '{}'.",
                                        requested_catalog_path, primary_err, default_catalog_path
                                    ));
                                    catalog
                                }
                                Err(default_err) => {
                                    return Err(EngineError {
                                        code: ErrorCode::InvalidInput,
                                        message: format!(
                                            "Could not open genome catalog '{}' ({}) and fallback '{}' ({})",
                                            requested_catalog_path,
                                            primary_err,
                                            default_catalog_path,
                                            default_err
                                        ),
                                    });
                                }
                            }
                        } else {
                            return Err(EngineError {
                                code: ErrorCode::InvalidInput,
                                message: format!(
                                    "Could not open genome catalog '{}': {}",
                                    requested_catalog_path, primary_err
                                ),
                            });
                        }
                    }
                };
                if let Some(warning) = catalog_fallback_warning {
                    result.warnings.push(warning);
                }

                let preferred_prepared = prepared_genome_id
                    .as_deref()
                    .map(str::trim)
                    .filter(|v| !v.is_empty())
                    .map(|v| v.to_string());
                let resolution_query = preferred_prepared
                    .as_deref()
                    .unwrap_or(&anchor.genome_id)
                    .to_string();
                let fallback_policy = if preferred_prepared.is_some() {
                    PreparedGenomeFallbackPolicy::Off
                } else {
                    self.genome_anchor_fallback_policy_for_extension()
                };
                let prepared_resolution = catalog
                    .resolve_prepared_genome_id_with_policy(
                        &resolution_query,
                        resolved_cache_dir.as_deref(),
                        fallback_policy,
                    )
                    .map_err(|e| EngineError {
                        code: ErrorCode::NotFound,
                        message: format!(
                            "Could not resolve prepared genome '{}' for verification: {}",
                            resolution_query, e
                        ),
                    })?;
                if let Some(warning) = prepared_resolution.fallback_warning {
                    result.warnings.push(warning);
                }
                let effective_genome_id = prepared_resolution.resolved_genome_id;
                let anchor_for_verification = GenomeSequenceAnchor {
                    genome_id: effective_genome_id.clone(),
                    chromosome: anchor.chromosome.clone(),
                    start_1based: anchor.start_1based,
                    end_1based: anchor.end_1based,
                    strand: anchor.strand,
                    anchor_verified: anchor.anchor_verified,
                    catalog_path: anchor.catalog_path.clone(),
                    cache_dir: anchor.cache_dir.clone(),
                };
                let is_match = Self::verify_anchor_sequence_against_catalog(
                    &dna,
                    &anchor_for_verification,
                    &resolved_catalog_path,
                    resolved_cache_dir.as_deref(),
                )
                .map_err(|e| EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!(
                        "Could not verify genome anchor '{}' against catalog '{}': {}",
                        seq_id, resolved_catalog_path, e
                    ),
                })?;
                let source_plan = catalog
                    .source_plan(&effective_genome_id, resolved_cache_dir.as_deref())
                    .ok();
                let inspection = catalog
                    .inspect_prepared_genome(&effective_genome_id, resolved_cache_dir.as_deref())
                    .ok()
                    .flatten();
                let (
                    sequence_source_type,
                    annotation_source_type,
                    sequence_source,
                    annotation_source,
                    sequence_sha1,
                    annotation_sha1,
                ) = Self::genome_source_snapshot(source_plan.as_ref(), inspection.as_ref());
                self.append_genome_extraction_provenance(GenomeExtractionProvenance {
                    seq_id: seq_id.clone(),
                    recorded_at_unix_ms: Self::now_unix_ms(),
                    operation: "VerifyGenomeAnchor".to_string(),
                    genome_id: effective_genome_id.clone(),
                    catalog_path: resolved_catalog_path.clone(),
                    cache_dir: resolved_cache_dir.clone(),
                    chromosome: Some(anchor.chromosome.clone()),
                    start_1based: Some(anchor.start_1based),
                    end_1based: Some(anchor.end_1based),
                    gene_query: None,
                    occurrence: None,
                    gene_id: None,
                    gene_name: None,
                    strand: anchor.strand,
                    anchor_strand: Some(anchor.strand.unwrap_or('+')),
                    anchor_verified: Some(is_match),
                    sequence_source_type,
                    annotation_source_type,
                    sequence_source,
                    annotation_source,
                    sequence_sha1,
                    annotation_sha1,
                });
                result.changed_seq_ids.push(seq_id.clone());
                if is_match {
                    result.messages.push(format!(
                        "Genome anchor for '{}' verified against '{}' via prepared '{}'",
                        seq_id, resolved_catalog_path, effective_genome_id
                    ));
                } else {
                    result.warnings.push(format!(
                        "Genome anchor for '{}' is unverified against '{}' via prepared '{}'",
                        seq_id, resolved_catalog_path, effective_genome_id
                    ));
                    result.messages.push(format!(
                        "Recorded anchor verification status for '{}' as unverified",
                        seq_id
                    ));
                }
            }
            Operation::ImportGenomeBedTrack {
                seq_id,
                path,
                track_name,
                min_score,
                max_score,
                clear_existing,
            } => {
                if path.trim().is_empty() {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: "ImportGenomeBedTrack requires a non-empty BED path".to_string(),
                    });
                }
                if min_score
                    .zip(max_score)
                    .map(|(min, max)| min > max)
                    .unwrap_or(false)
                {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: "ImportGenomeBedTrack requires min_score <= max_score".to_string(),
                    });
                }

                let anchor = self.latest_genome_anchor_for_seq(&seq_id)?;
                let _ = self.ensure_lineage_node(&seq_id);
                let dna = self
                    .state
                    .sequences
                    .get_mut(&seq_id)
                    .ok_or_else(|| EngineError {
                        code: ErrorCode::NotFound,
                        message: format!("Sequence '{seq_id}' not found"),
                    })?;
                let seq_id_for_progress = seq_id.clone();
                let path_for_progress = path.clone();
                let mut progress_cb = |parsed_records: usize,
                                       imported_features: usize,
                                       skipped_records: usize,
                                       done: bool| {
                    on_progress(OperationProgress::GenomeTrackImport(
                        GenomeTrackImportProgress {
                            seq_id: seq_id_for_progress.clone(),
                            source: "BED".to_string(),
                            path: path_for_progress.clone(),
                            parsed_records,
                            imported_features,
                            skipped_records,
                            done,
                        },
                    ))
                };

                let report = Self::import_genome_bed_track(
                    dna,
                    &anchor,
                    &path,
                    track_name.as_deref(),
                    min_score,
                    max_score,
                    clear_existing.unwrap_or(false),
                    Some(&mut progress_cb),
                )?;

                result.changed_seq_ids.push(seq_id.clone());
                result.warnings.extend(report.warnings);
                let anchor_strand = anchor.strand.unwrap_or('+');
                result.messages.push(format!(
                    "Imported {} BED feature(s) into '{}' from '{}' as track '{}' (anchor={} {}:{}-{} strand {}, parsed={}, skipped={})",
                    report.imported_features,
                    seq_id,
                    path,
                    report.track_name,
                    anchor.genome_id,
                    anchor.chromosome,
                    anchor.start_1based,
                    anchor.end_1based,
                    anchor_strand,
                    report.parsed_records,
                    report.skipped_records
                ));
                if report.skipped_missing_score > 0 {
                    result.warnings.push(format!(
                        "{} BED record(s) were skipped because score filters were set but the BED score column was missing",
                        report.skipped_missing_score
                    ));
                }
                if report.skipped_outside_score_range > 0 {
                    result.messages.push(format!(
                        "{} BED record(s) were outside score filter bounds",
                        report.skipped_outside_score_range
                    ));
                }
                if report.truncated_at_limit {
                    result.warnings.push(format!(
                        "BED import was truncated after {} features (limit={})",
                        report.imported_features, MAX_IMPORTED_SIGNAL_FEATURES
                    ));
                }
            }
            Operation::ImportGenomeBigWigTrack {
                seq_id,
                path,
                track_name,
                min_score,
                max_score,
                clear_existing,
            } => {
                if path.trim().is_empty() {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: "ImportGenomeBigWigTrack requires a non-empty BigWig path"
                            .to_string(),
                    });
                }
                if min_score
                    .zip(max_score)
                    .map(|(min, max)| min > max)
                    .unwrap_or(false)
                {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: "ImportGenomeBigWigTrack requires min_score <= max_score"
                            .to_string(),
                    });
                }

                let anchor = self.latest_genome_anchor_for_seq(&seq_id)?;
                let _ = self.ensure_lineage_node(&seq_id);
                let dna = self
                    .state
                    .sequences
                    .get_mut(&seq_id)
                    .ok_or_else(|| EngineError {
                        code: ErrorCode::NotFound,
                        message: format!("Sequence '{seq_id}' not found"),
                    })?;
                let seq_id_for_progress = seq_id.clone();
                let path_for_progress = path.clone();
                let mut progress_cb = |parsed_records: usize,
                                       imported_features: usize,
                                       skipped_records: usize,
                                       done: bool| {
                    on_progress(OperationProgress::GenomeTrackImport(
                        GenomeTrackImportProgress {
                            seq_id: seq_id_for_progress.clone(),
                            source: "BigWig".to_string(),
                            path: path_for_progress.clone(),
                            parsed_records,
                            imported_features,
                            skipped_records,
                            done,
                        },
                    ))
                };

                let report = Self::import_genome_bigwig_track(
                    dna,
                    &anchor,
                    &path,
                    track_name.as_deref(),
                    min_score,
                    max_score,
                    clear_existing.unwrap_or(false),
                    Some(&mut progress_cb),
                )?;

                result.changed_seq_ids.push(seq_id.clone());
                result.warnings.extend(report.warnings);
                let anchor_strand = anchor.strand.unwrap_or('+');
                result.messages.push(format!(
                    "Imported {} BigWig feature(s) into '{}' from '{}' as track '{}' (anchor={} {}:{}-{} strand {}, parsed={}, skipped={})",
                    report.imported_features,
                    seq_id,
                    path,
                    report.track_name,
                    anchor.genome_id,
                    anchor.chromosome,
                    anchor.start_1based,
                    anchor.end_1based,
                    anchor_strand,
                    report.parsed_records,
                    report.skipped_records
                ));
                if report.skipped_missing_score > 0 {
                    result.warnings.push(format!(
                        "{} converted bedGraph record(s) were skipped because score filters were set but no value was available",
                        report.skipped_missing_score
                    ));
                }
                if report.skipped_outside_score_range > 0 {
                    result.messages.push(format!(
                        "{} converted bedGraph record(s) were outside score filter bounds",
                        report.skipped_outside_score_range
                    ));
                }
                if report.truncated_at_limit {
                    result.warnings.push(format!(
                        "BigWig import was truncated after {} features (limit={})",
                        report.imported_features, MAX_IMPORTED_SIGNAL_FEATURES
                    ));
                }
            }
            Operation::ImportGenomeVcfTrack {
                seq_id,
                path,
                track_name,
                min_score,
                max_score,
                clear_existing,
            } => {
                if path.trim().is_empty() {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: "ImportGenomeVcfTrack requires a non-empty VCF path".to_string(),
                    });
                }
                if min_score
                    .zip(max_score)
                    .map(|(min, max)| min > max)
                    .unwrap_or(false)
                {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: "ImportGenomeVcfTrack requires min_score <= max_score".to_string(),
                    });
                }

                let anchor = self.latest_genome_anchor_for_seq(&seq_id)?;
                let _ = self.ensure_lineage_node(&seq_id);
                let dna = self
                    .state
                    .sequences
                    .get_mut(&seq_id)
                    .ok_or_else(|| EngineError {
                        code: ErrorCode::NotFound,
                        message: format!("Sequence '{seq_id}' not found"),
                    })?;
                let seq_id_for_progress = seq_id.clone();
                let path_for_progress = path.clone();
                let mut progress_cb = |parsed_records: usize,
                                       imported_features: usize,
                                       skipped_records: usize,
                                       done: bool| {
                    on_progress(OperationProgress::GenomeTrackImport(
                        GenomeTrackImportProgress {
                            seq_id: seq_id_for_progress.clone(),
                            source: "VCF".to_string(),
                            path: path_for_progress.clone(),
                            parsed_records,
                            imported_features,
                            skipped_records,
                            done,
                        },
                    ))
                };

                let report = Self::import_genome_vcf_track(
                    dna,
                    &anchor,
                    &path,
                    track_name.as_deref(),
                    min_score,
                    max_score,
                    clear_existing.unwrap_or(false),
                    Some(&mut progress_cb),
                )?;

                result.changed_seq_ids.push(seq_id.clone());
                result.warnings.extend(report.warnings);
                let anchor_strand = anchor.strand.unwrap_or('+');
                result.messages.push(format!(
                    "Imported {} VCF feature(s) into '{}' from '{}' as track '{}' (anchor={} {}:{}-{} strand {}, parsed={}, skipped={})",
                    report.imported_features,
                    seq_id,
                    path,
                    report.track_name,
                    anchor.genome_id,
                    anchor.chromosome,
                    anchor.start_1based,
                    anchor.end_1based,
                    anchor_strand,
                    report.parsed_records,
                    report.skipped_records
                ));
                if report.skipped_missing_score > 0 {
                    result.warnings.push(format!(
                        "{} VCF record(s) were skipped because QUAL-based score filters were set but QUAL was missing",
                        report.skipped_missing_score
                    ));
                }
                if report.skipped_outside_score_range > 0 {
                    result.messages.push(format!(
                        "{} VCF record(s) were outside QUAL score filter bounds",
                        report.skipped_outside_score_range
                    ));
                }
                if report.truncated_at_limit {
                    result.warnings.push(format!(
                        "VCF import was truncated after {} features (limit={})",
                        report.imported_features, MAX_IMPORTED_SIGNAL_FEATURES
                    ));
                }
            }
            Operation::ImportIsoformPanel {
                seq_id,
                panel_path,
                panel_id,
                strict,
            } => {
                if panel_path.trim().is_empty() {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: "ImportIsoformPanel requires a non-empty panel_path".to_string(),
                    });
                }
                if !self.state.sequences.contains_key(&seq_id) {
                    return Err(EngineError {
                        code: ErrorCode::NotFound,
                        message: format!("Sequence '{seq_id}' not found"),
                    });
                }
                let resource = Self::load_isoform_panel_resource(&panel_path, panel_id.as_deref())?;
                let panel_id = resource.panel_id.clone();
                let preview = self.build_isoform_architecture_expert_view_from_resource(
                    &seq_id, &panel_id, &resource, strict,
                )?;
                self.upsert_isoform_panel_record(IsoformPanelRecord {
                    seq_id: seq_id.clone(),
                    panel_id: panel_id.clone(),
                    imported_at_unix_ms: Self::now_unix_ms(),
                    source_path: panel_path.clone(),
                    strict,
                    resource: resource.clone(),
                })?;
                result.changed_seq_ids.push(seq_id.clone());
                result.warnings.extend(preview.warnings.clone());
                result.messages.push(format!(
                    "Imported isoform panel '{}' for '{}' (isoforms={}, strict={}) from '{}'",
                    panel_id,
                    seq_id,
                    resource.isoforms.len(),
                    strict,
                    panel_path
                ));
            }
            Operation::ImportUniprotSwissProt { path, entry_id } => {
                if path.trim().is_empty() {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: "ImportUniprotSwissProt requires a non-empty path".to_string(),
                    });
                }
                let text = std::fs::read_to_string(&path).map_err(|e| EngineError {
                    code: ErrorCode::Io,
                    message: format!("Could not read SWISS-PROT text from '{}': {e}", path),
                })?;
                let source = format!("file://{}", path);
                let entry =
                    Self::parse_uniprot_entry_text(&text, &source, None, entry_id.as_deref())?;
                let resolved_entry_id = entry.entry_id.clone();
                self.upsert_uniprot_entry(entry)?;
                result.messages.push(format!(
                    "Imported UniProt SWISS-PROT entry '{}' from '{}'",
                    resolved_entry_id, path
                ));
            }
            Operation::FetchUniprotSwissProt { query, entry_id } => {
                let query_trimmed = query.trim();
                if query_trimmed.is_empty() {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: "FetchUniprotSwissProt requires a non-empty query".to_string(),
                    });
                }
                let (source_url, text) = Self::fetch_uniprot_swiss_prot_text(query_trimmed)?;
                let entry = Self::parse_uniprot_entry_text(
                    &text,
                    &source_url,
                    Some(query_trimmed),
                    entry_id.as_deref(),
                )?;
                let resolved_entry_id = entry.entry_id.clone();
                let accession = entry.accession.clone();
                self.upsert_uniprot_entry(entry)?;
                result.messages.push(format!(
                    "Fetched UniProt entry '{}' (accession '{}') from '{}'",
                    resolved_entry_id, accession, source_url
                ));
            }
            Operation::FetchGenBankAccession { accession, as_id } => {
                let accession_trimmed = accession.trim();
                if accession_trimmed.is_empty() {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: "FetchGenBankAccession requires a non-empty accession".to_string(),
                    });
                }
                let _ = self.fetch_genbank_accession_into_state(
                    &mut result,
                    accession_trimmed,
                    as_id,
                    "FetchGenBankAccession",
                    None,
                )?;
            }
            Operation::FetchUniprotLinkedGenBank {
                entry_id,
                accession,
                as_id,
            } => {
                let entry_id_trimmed = entry_id.trim();
                if entry_id_trimmed.is_empty() {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: "FetchUniprotLinkedGenBank requires a non-empty entry_id"
                            .to_string(),
                    });
                }
                let entry = self.get_uniprot_entry(entry_id_trimmed)?;
                let selected = Self::choose_uniprot_nucleotide_xref(&entry, accession.as_deref())?;
                let source_note = format!(
                    "UniProt '{}' -> {} accession '{}'{}",
                    entry.entry_id,
                    selected.database,
                    selected.accession,
                    selected
                        .molecule_type
                        .as_deref()
                        .map(|kind| format!(" ({})", kind))
                        .unwrap_or_default()
                );
                let _ = self.fetch_genbank_accession_into_state(
                    &mut result,
                    &selected.accession,
                    as_id,
                    "FetchUniprotLinkedGenBank",
                    Some(&source_note),
                )?;
            }
            Operation::ImportUniprotEntrySequence {
                entry_id,
                output_id,
            } => {
                let entry_id = entry_id.trim();
                if entry_id.is_empty() {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: "ImportUniprotEntrySequence requires a non-empty entry_id"
                            .to_string(),
                    });
                }
                let _ = self.import_uniprot_entry_sequence(
                    &mut result,
                    entry_id,
                    output_id.as_deref(),
                )?;
            }
            Operation::ProjectUniprotToGenome {
                seq_id,
                entry_id,
                projection_id,
                transcript_id,
            } => {
                let projection = self.project_uniprot_to_genome(
                    &seq_id,
                    &entry_id,
                    projection_id.as_deref(),
                    transcript_id.as_deref(),
                )?;
                let projection_id = projection.projection_id.clone();
                let transcript_count = projection.transcript_projections.len();
                result.warnings.extend(projection.warnings.clone());
                for transcript in &projection.transcript_projections {
                    result.warnings.extend(transcript.warnings.clone());
                }
                self.upsert_uniprot_projection(projection)?;
                result.messages.push(format!(
                    "Projected UniProt entry '{}' onto '{}' as '{}' (transcripts={})",
                    entry_id, seq_id, projection_id, transcript_count
                ));
            }
            Operation::ImportBlastHitsTrack {
                seq_id,
                hits,
                track_name,
                clear_existing,
                blast_provenance,
            } => {
                if hits.is_empty() {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: "ImportBlastHitsTrack requires at least one hit".to_string(),
                    });
                }
                let selected_track_name = track_name
                    .as_deref()
                    .map(str::trim)
                    .filter(|v| !v.is_empty())
                    .unwrap_or("blast_hits")
                    .to_string();
                let clear_existing = clear_existing.unwrap_or(false);

                let _ = self.ensure_lineage_node(&seq_id);
                let dna = self
                    .state
                    .sequences
                    .get_mut(&seq_id)
                    .ok_or_else(|| EngineError {
                        code: ErrorCode::NotFound,
                        message: format!("Sequence '{seq_id}' not found"),
                    })?;
                let seq_len = dna.len();
                if seq_len == 0 {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "Sequence '{}' is empty and cannot receive BLAST hit features",
                            seq_id
                        ),
                    });
                }

                let mut removed_count = 0usize;
                if clear_existing {
                    let before = dna.features().len();
                    Self::remove_generated_blast_hit_features(dna.features_mut());
                    removed_count = before.saturating_sub(dna.features().len());
                }

                let mut imported_count = 0usize;
                let mut skipped_count = 0usize;
                for (idx, hit) in hits.iter().enumerate() {
                    let start = hit.query_start_1based.min(hit.query_end_1based);
                    let end = hit.query_start_1based.max(hit.query_end_1based);
                    if start == 0 || end == 0 {
                        skipped_count += 1;
                        if result.warnings.len() < 20 {
                            result.warnings.push(format!(
                                "BLAST hit {} skipped because query coordinates are invalid: {}..{}",
                                idx + 1,
                                hit.query_start_1based,
                                hit.query_end_1based
                            ));
                        }
                        continue;
                    }
                    if start > seq_len {
                        skipped_count += 1;
                        if result.warnings.len() < 20 {
                            result.warnings.push(format!(
                                "BLAST hit {} skipped because query start {} is outside sequence length {}",
                                idx + 1,
                                start,
                                seq_len
                            ));
                        }
                        continue;
                    }
                    let end_clamped = end.min(seq_len);
                    if end_clamped < start {
                        skipped_count += 1;
                        continue;
                    }
                    if end_clamped < end && result.warnings.len() < 20 {
                        result.warnings.push(format!(
                            "BLAST hit {} query range {}..{} was clamped to {}..{} for sequence length {}",
                            idx + 1,
                            start,
                            end,
                            start,
                            end_clamped,
                            seq_len
                        ));
                    }
                    let local_start_0based = start.saturating_sub(1);
                    let local_end_0based_exclusive = end_clamped;
                    let local_strand =
                        if hit.subject_start_1based == 0 || hit.subject_end_1based == 0 {
                            None
                        } else if hit.subject_start_1based <= hit.subject_end_1based {
                            Some('+')
                        } else {
                            Some('-')
                        };
                    let feature = Self::build_blast_hit_feature(
                        hit,
                        &selected_track_name,
                        local_start_0based,
                        local_end_0based_exclusive,
                        local_strand,
                    );
                    dna.features_mut().push(feature);
                    imported_count += 1;
                }

                if imported_count > 0 || removed_count > 0 {
                    result.changed_seq_ids.push(seq_id.clone());
                }
                result.messages.push(format!(
                    "Imported {} BLAST hit feature(s) into '{}' as track '{}' (input_hits={}, skipped={}, cleared_existing={})",
                    imported_count,
                    seq_id,
                    selected_track_name,
                    hits.len(),
                    skipped_count,
                    removed_count
                ));
                if let Some(provenance) = blast_provenance {
                    let command_line = if provenance.command_line.trim().is_empty() {
                        if provenance.command.is_empty() {
                            provenance.blastn_executable.clone()
                        } else {
                            format!(
                                "{} {}",
                                provenance.blastn_executable,
                                provenance.command.join(" ")
                            )
                        }
                    } else {
                        provenance.command_line.clone()
                    };
                    result.messages.push(format!(
                        "BLAST provenance: genome='{}', query='{}' ({} bp), task='{}', max_hits={}, command={}",
                        provenance.genome_id,
                        provenance.query_label,
                        provenance.query_length,
                        provenance.task,
                        provenance.max_hits,
                        command_line
                    ));
                }
            }
            Operation::Digest {
                input,
                enzymes,
                output_prefix,
            } => {
                parent_seq_ids.push(input.clone());
                let dna = self
                    .state
                    .sequences
                    .get(&input)
                    .ok_or_else(|| EngineError {
                        code: ErrorCode::NotFound,
                        message: format!("Sequence '{input}' not found"),
                    })?
                    .clone();

                let (found, missing) = self.resolve_enzymes(&enzymes)?;
                if !missing.is_empty() {
                    result
                        .warnings
                        .push(format!("Unknown enzymes ignored: {}", missing.join(",")));
                }

                let fragments =
                    Self::digest_with_guard(&dna, found, self.max_fragments_per_container())?;
                let prefix = output_prefix.unwrap_or_else(|| format!("{input}_digest"));

                for (i, mut fragment) in fragments.into_iter().enumerate() {
                    // Keep digest interactive by deferring expensive feature recomputation.
                    Self::prepare_sequence_light(&mut fragment);
                    let candidate = format!("{}_{}", prefix, i + 1);
                    let seq_id = self.unique_seq_id(&candidate);
                    self.state.sequences.insert(seq_id.clone(), fragment);
                    self.add_lineage_node(&seq_id, SequenceOrigin::Derived, Some(&result.op_id));
                    result.created_seq_ids.push(seq_id);
                }

                result.messages.push(format!(
                    "Digest created {} fragment(s); feature recomputation deferred",
                    result.created_seq_ids.len()
                ));
            }
            Operation::DigestContainer {
                container_id,
                enzymes,
                output_prefix,
            } => {
                let inputs = self.container_members(&container_id)?;
                parent_seq_ids.extend(inputs.clone());
                let (found, missing) = self.resolve_enzymes(&enzymes)?;
                if !missing.is_empty() {
                    result
                        .warnings
                        .push(format!("Unknown enzymes ignored: {}", missing.join(",")));
                }
                let prefix = output_prefix.unwrap_or_else(|| format!("{container_id}_digest"));
                for input in inputs {
                    let dna = self
                        .state
                        .sequences
                        .get(&input)
                        .ok_or_else(|| EngineError {
                            code: ErrorCode::NotFound,
                            message: format!("Sequence '{input}' not found"),
                        })?
                        .clone();
                    let fragments = Self::digest_with_guard(
                        &dna,
                        found.clone(),
                        self.max_fragments_per_container(),
                    )?;
                    for (i, mut fragment) in fragments.into_iter().enumerate() {
                        // Keep digest interactive by deferring expensive feature recomputation.
                        Self::prepare_sequence_light(&mut fragment);
                        let candidate = format!("{}_{}_{}", prefix, input, i + 1);
                        let seq_id = self.unique_seq_id(&candidate);
                        self.state.sequences.insert(seq_id.clone(), fragment);
                        self.add_lineage_node(
                            &seq_id,
                            SequenceOrigin::Derived,
                            Some(&result.op_id),
                        );
                        result.created_seq_ids.push(seq_id);
                        if result.created_seq_ids.len() > self.max_fragments_per_container() {
                            return Err(EngineError {
                                code: ErrorCode::InvalidInput,
                                message: format!(
                                    "Digest produced more than max_fragments_per_container={}",
                                    self.max_fragments_per_container()
                                ),
                            });
                        }
                    }
                }

                result.messages.push(format!(
                    "Digest container '{}' created {} fragment(s); feature recomputation deferred",
                    container_id,
                    result.created_seq_ids.len()
                ));
            }
            Operation::MergeContainersById { .. }
            | Operation::LigationContainer { .. }
            | Operation::FilterContainerByMolecularWeight { .. } => {
                unreachable!("container operation variants are normalized before execution")
            }
            Operation::MergeContainers {
                inputs,
                output_prefix,
            } => {
                if inputs.is_empty() {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: "MergeContainers requires at least one input sequence".to_string(),
                    });
                }
                if inputs.len() > self.max_fragments_per_container() {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "MergeContainers input count {} exceeds max_fragments_per_container={}",
                            inputs.len(),
                            self.max_fragments_per_container()
                        ),
                    });
                }
                parent_seq_ids.extend(inputs.clone());
                let prefix = output_prefix.unwrap_or_else(|| "merged".to_string());
                for (i, input) in inputs.iter().enumerate() {
                    let dna = self
                        .state
                        .sequences
                        .get(input)
                        .ok_or_else(|| EngineError {
                            code: ErrorCode::NotFound,
                            message: format!("Sequence '{input}' not found"),
                        })?
                        .clone();
                    let seq_id = self.unique_seq_id(&format!("{}_{}", prefix, i + 1));
                    self.state.sequences.insert(seq_id.clone(), dna);
                    self.add_lineage_node(&seq_id, SequenceOrigin::Derived, Some(&result.op_id));
                    result.created_seq_ids.push(seq_id);
                }
                result.messages.push(format!(
                    "Merged {} input sequence(s) into container prefix '{}'",
                    inputs.len(),
                    prefix
                ));
            }
            Operation::Ligation {
                inputs,
                circularize_if_possible,
                output_id,
                protocol,
                output_prefix,
                unique,
            } => {
                parent_seq_ids.extend(inputs.clone());
                if inputs.len() < 2 {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: "Ligation requires at least two input sequences".to_string(),
                    });
                }
                if 1 > self.max_fragments_per_container() {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "Ligation product count exceeds max_fragments_per_container={}",
                            self.max_fragments_per_container()
                        ),
                    });
                }
                let mut accepted: Vec<(String, String, String)> = vec![];
                for (i, left_id) in inputs.iter().enumerate() {
                    for (j, right_id) in inputs.iter().enumerate() {
                        if i == j {
                            continue;
                        }
                        let left =
                            self.state
                                .sequences
                                .get(left_id)
                                .ok_or_else(|| EngineError {
                                    code: ErrorCode::NotFound,
                                    message: format!("Sequence '{left_id}' not found"),
                                })?;
                        let right =
                            self.state
                                .sequences
                                .get(right_id)
                                .ok_or_else(|| EngineError {
                                    code: ErrorCode::NotFound,
                                    message: format!("Sequence '{right_id}' not found"),
                                })?;

                        let ok = match protocol {
                            LigationProtocol::Sticky => Self::sticky_compatible(left, right),
                            LigationProtocol::Blunt => {
                                Self::right_end_is_blunt(left) && Self::left_end_is_blunt(right)
                            }
                        };
                        if !ok {
                            continue;
                        }

                        let product = format!(
                            "{}{}",
                            left.get_forward_string(),
                            right.get_forward_string()
                        );
                        accepted.push((left_id.clone(), right_id.clone(), product));
                        if accepted.len() > self.max_fragments_per_container() {
                            return Err(EngineError {
                                code: ErrorCode::InvalidInput,
                                message: format!(
                                    "Ligation produced more than max_fragments_per_container={}",
                                    self.max_fragments_per_container()
                                ),
                            });
                        }
                    }
                }

                if accepted.is_empty() {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "No ligation products found for protocol '{:?}'",
                            protocol
                        ),
                    });
                }

                let require_unique = unique.unwrap_or(false);
                if require_unique && accepted.len() != 1 {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "Ligation unique=true requires exactly one product, found {}",
                            accepted.len()
                        ),
                    });
                }
                if output_id.is_some() && accepted.len() != 1 {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: "Ligation output_id can only be used when exactly one product is produced"
                            .to_string(),
                    });
                }

                let prefix = output_prefix.unwrap_or_else(|| "ligation".to_string());
                for (idx, (left_id, right_id, merged)) in accepted.into_iter().enumerate() {
                    let mut product =
                        DNAsequence::from_sequence(&merged).map_err(|e| EngineError {
                            code: ErrorCode::Internal,
                            message: format!("Could not create ligation product: {e}"),
                        })?;
                    product.set_circular(circularize_if_possible);
                    Self::prepare_sequence(&mut product);

                    let seq_id = if idx == 0 {
                        if let Some(id) = output_id.clone() {
                            self.unique_seq_id(&id)
                        } else {
                            self.unique_seq_id(&format!("{}_{}", prefix, idx + 1))
                        }
                    } else {
                        self.unique_seq_id(&format!("{}_{}", prefix, idx + 1))
                    };
                    self.state.sequences.insert(seq_id.clone(), product);
                    self.add_lineage_node(&seq_id, SequenceOrigin::Derived, Some(&result.op_id));
                    result.created_seq_ids.push(seq_id.clone());
                    result.messages.push(format!(
                        "Ligation product '{}' from '{}' + '{}'",
                        seq_id, left_id, right_id
                    ));
                }
            }
            Operation::Pcr {
                template,
                forward_primer,
                reverse_primer,
                output_id,
                unique,
            } => {
                parent_seq_ids.push(template.clone());
                let dna = self
                    .state
                    .sequences
                    .get(&template)
                    .ok_or_else(|| EngineError {
                        code: ErrorCode::NotFound,
                        message: format!("Sequence '{template}' not found"),
                    })?
                    .clone();

                if dna.is_circular() {
                    return Err(EngineError {
                        code: ErrorCode::Unsupported,
                        message: "PCR on circular templates is not implemented yet".to_string(),
                    });
                }

                let template_seq = dna.get_forward_string().to_ascii_uppercase();
                let template_bytes = template_seq.as_bytes();
                let fwd = Self::normalize_dna_text(&forward_primer);
                let rev = Self::normalize_dna_text(&reverse_primer);
                if fwd.is_empty() || rev.is_empty() {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: "PCR primers must not be empty".to_string(),
                    });
                }

                let rev_binding = Self::reverse_complement(&rev);
                let fwd_sites = Self::find_all_subsequences(template_bytes, fwd.as_bytes());
                if fwd_sites.is_empty() {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: "Forward primer not found on template".to_string(),
                    });
                }
                let rev_sites = Self::find_all_subsequences(template_bytes, rev_binding.as_bytes());
                if rev_sites.is_empty() {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: "Reverse primer binding site not found on template".to_string(),
                    });
                }

                let mut amplicon_ranges: Vec<(usize, usize)> = vec![];
                for fwd_pos in &fwd_sites {
                    for rev_pos in &rev_sites {
                        if *rev_pos < *fwd_pos {
                            continue;
                        }
                        let amplicon_end = rev_pos + rev_binding.len();
                        if amplicon_end <= template_seq.len() {
                            amplicon_ranges.push((*fwd_pos, amplicon_end));
                        }
                    }
                }
                amplicon_ranges.sort_unstable();
                amplicon_ranges.dedup();

                if amplicon_ranges.is_empty() {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: "No valid forward/reverse primer pair produced an amplicon"
                            .to_string(),
                    });
                }
                if amplicon_ranges.len() > self.max_fragments_per_container() {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "PCR produced {} amplicons, exceeding max_fragments_per_container={}",
                            amplicon_ranges.len(),
                            self.max_fragments_per_container()
                        ),
                    });
                }

                let require_unique = unique.unwrap_or(false);
                if require_unique && amplicon_ranges.len() != 1 {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "PCR unique=true requires exactly one amplicon, found {}",
                            amplicon_ranges.len()
                        ),
                    });
                }
                if output_id.is_some() && amplicon_ranges.len() != 1 {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message:
                            "PCR output_id can only be used when exactly one amplicon is produced"
                                .to_string(),
                    });
                }

                let default_base = format!("{template}_pcr");
                for (i, (start, end)) in amplicon_ranges.iter().enumerate() {
                    let amplicon = &template_seq[*start..*end];
                    let mut pcr_product =
                        DNAsequence::from_sequence(amplicon).map_err(|e| EngineError {
                            code: ErrorCode::Internal,
                            message: format!("Could not create PCR product: {e}"),
                        })?;
                    pcr_product.set_circular(false);
                    Self::prepare_sequence(&mut pcr_product);

                    let requested = if i == 0 { output_id.clone() } else { None };
                    let seq_id = if let Some(id) = requested {
                        self.unique_seq_id(&id)
                    } else if amplicon_ranges.len() == 1 {
                        self.unique_seq_id(&default_base)
                    } else {
                        self.unique_seq_id(&format!("{}_{}", default_base, i + 1))
                    };
                    self.state.sequences.insert(seq_id.clone(), pcr_product);
                    self.add_lineage_node(&seq_id, SequenceOrigin::Derived, Some(&result.op_id));
                    result.created_seq_ids.push(seq_id.clone());
                    result.messages.push(format!(
                        "PCR product '{}' created: {}..{} (len {})",
                        seq_id,
                        start,
                        end,
                        end - start
                    ));
                }
            }
            Operation::PcrAdvanced {
                template,
                forward_primer,
                reverse_primer,
                output_id,
                unique,
            } => {
                parent_seq_ids.push(template.clone());
                let dna = self
                    .state
                    .sequences
                    .get(&template)
                    .ok_or_else(|| EngineError {
                        code: ErrorCode::NotFound,
                        message: format!("Sequence '{template}' not found"),
                    })?
                    .clone();

                if dna.is_circular() {
                    return Err(EngineError {
                        code: ErrorCode::Unsupported,
                        message: "PCR on circular templates is not implemented yet".to_string(),
                    });
                }

                let template_seq = dna.get_forward_string().to_ascii_uppercase();
                let template_bytes = template_seq.as_bytes();

                let fwd_variants = Self::expand_primer_variants(
                    &forward_primer,
                    self.max_fragments_per_container(),
                )?;
                let rev_variants = Self::expand_primer_variants(
                    &reverse_primer,
                    self.max_fragments_per_container(),
                )?;
                if fwd_variants.is_empty() || rev_variants.is_empty() {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: "PCR primers must not be empty".to_string(),
                    });
                }

                let mut candidates: Vec<(usize, usize, String)> = vec![];
                let mut seen_amplicons: HashSet<String> = HashSet::new();
                for fwd_full in &fwd_variants {
                    let fwd_anneal_len = forward_primer.anneal_len.unwrap_or(fwd_full.len());
                    if fwd_anneal_len == 0 || fwd_anneal_len > fwd_full.len() {
                        return Err(EngineError {
                            code: ErrorCode::InvalidInput,
                            message: "Invalid anneal_len in PCR primer spec".to_string(),
                        });
                    }
                    let fwd_anneal = &fwd_full[fwd_full.len() - fwd_anneal_len..];
                    let fwd_sites = Self::find_anneal_sites(
                        template_bytes,
                        fwd_anneal.as_bytes(),
                        forward_primer.max_mismatches.unwrap_or(0),
                        forward_primer.require_3prime_exact_bases.unwrap_or(0),
                        true,
                    );
                    if fwd_sites.is_empty() {
                        continue;
                    }

                    for rev_full in &rev_variants {
                        let rev_anneal_len = reverse_primer.anneal_len.unwrap_or(rev_full.len());
                        if rev_anneal_len == 0 || rev_anneal_len > rev_full.len() {
                            return Err(EngineError {
                                code: ErrorCode::InvalidInput,
                                message: "Invalid anneal_len in PCR primer spec".to_string(),
                            });
                        }
                        let rev_anneal = &rev_full[rev_full.len() - rev_anneal_len..];
                        let rev_binding = Self::reverse_complement(rev_anneal);
                        let rev_sites = Self::find_anneal_sites(
                            template_bytes,
                            rev_binding.as_bytes(),
                            reverse_primer.max_mismatches.unwrap_or(0),
                            reverse_primer.require_3prime_exact_bases.unwrap_or(0),
                            false,
                        );
                        if rev_sites.is_empty() {
                            continue;
                        }
                        let rev_full_rc = Self::reverse_complement(rev_full);
                        for fwd_pos in &fwd_sites {
                            for rev_pos in &rev_sites {
                                let interior_start = *fwd_pos + fwd_anneal_len;
                                let interior_end = *rev_pos;
                                if interior_start > interior_end {
                                    continue;
                                }
                                let interior = &template_seq[interior_start..interior_end];
                                let amplicon = format!("{fwd_full}{interior}{rev_full_rc}");
                                if seen_amplicons.insert(amplicon.clone()) {
                                    candidates.push((*fwd_pos, *rev_pos, amplicon));
                                    if candidates.len() > self.max_fragments_per_container() {
                                        return Err(EngineError {
                                            code: ErrorCode::InvalidInput,
                                            message: format!(
                                                "PCR produced more than max_fragments_per_container={}",
                                                self.max_fragments_per_container()
                                            ),
                                        });
                                    }
                                }
                            }
                        }
                    }
                }
                if candidates.is_empty() {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: "No valid advanced PCR amplicon could be formed".to_string(),
                    });
                }

                let require_unique = unique.unwrap_or(false);
                if require_unique && candidates.len() != 1 {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "PCR unique=true requires exactly one amplicon, found {}",
                            candidates.len()
                        ),
                    });
                }
                if output_id.is_some() && candidates.len() != 1 {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message:
                            "PCR output_id can only be used when exactly one amplicon is produced"
                                .to_string(),
                    });
                }

                let default_base = format!("{template}_pcr");
                for (i, (fwd_pos, rev_pos, amplicon)) in candidates.into_iter().enumerate() {
                    let mut pcr_product =
                        DNAsequence::from_sequence(&amplicon).map_err(|e| EngineError {
                            code: ErrorCode::Internal,
                            message: format!("Could not create PCR product: {e}"),
                        })?;
                    pcr_product.set_circular(false);
                    Self::prepare_sequence(&mut pcr_product);

                    let requested = if i == 0 { output_id.clone() } else { None };
                    let seq_id = if let Some(id) = requested {
                        self.unique_seq_id(&id)
                    } else if seen_amplicons.len() == 1 {
                        self.unique_seq_id(&default_base)
                    } else {
                        self.unique_seq_id(&format!("{}_{}", default_base, i + 1))
                    };
                    self.state.sequences.insert(seq_id.clone(), pcr_product);
                    self.add_lineage_node(&seq_id, SequenceOrigin::Derived, Some(&result.op_id));
                    result.created_seq_ids.push(seq_id.clone());
                    result.messages.push(format!(
                        "Advanced PCR product '{}' created from fwd@{} rev@{} (len {})",
                        seq_id,
                        fwd_pos,
                        rev_pos,
                        amplicon.len()
                    ));
                }
            }
            Operation::PcrMutagenesis {
                template,
                forward_primer,
                reverse_primer,
                mutations,
                output_id,
                unique,
                require_all_mutations,
            } => {
                parent_seq_ids.push(template.clone());
                let dna = self
                    .state
                    .sequences
                    .get(&template)
                    .ok_or_else(|| EngineError {
                        code: ErrorCode::NotFound,
                        message: format!("Sequence '{template}' not found"),
                    })?
                    .clone();

                if dna.is_circular() {
                    return Err(EngineError {
                        code: ErrorCode::Unsupported,
                        message: "PCR on circular templates is not implemented yet".to_string(),
                    });
                }
                if mutations.is_empty() {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: "PcrMutagenesis requires at least one mutation".to_string(),
                    });
                }

                let template_seq = dna.get_forward_string().to_ascii_uppercase();
                let template_bytes = template_seq.as_bytes();

                let mut normalized_mutations: Vec<(usize, u8, u8)> = vec![];
                for m in &mutations {
                    if m.zero_based_position >= template_bytes.len() {
                        return Err(EngineError {
                            code: ErrorCode::InvalidInput,
                            message: format!(
                                "Mutation position {} is out of bounds for template length {}",
                                m.zero_based_position,
                                template_bytes.len()
                            ),
                        });
                    }
                    let ref_nt = Self::normalize_dna_text(&m.reference).to_ascii_uppercase();
                    let alt_nt = Self::normalize_dna_text(&m.alternate).to_ascii_uppercase();
                    if ref_nt.len() != 1 || alt_nt.len() != 1 {
                        return Err(EngineError {
                            code: ErrorCode::InvalidInput,
                            message: "Mutation reference/alternate must be single nucleotides"
                                .to_string(),
                        });
                    }
                    let ref_b = ref_nt.as_bytes()[0];
                    let alt_b = alt_nt.as_bytes()[0];
                    if template_bytes[m.zero_based_position] != ref_b {
                        return Err(EngineError {
                            code: ErrorCode::InvalidInput,
                            message: format!(
                                "Mutation reference mismatch at position {}: template has '{}', expected '{}'",
                                m.zero_based_position,
                                template_bytes[m.zero_based_position] as char,
                                ref_b as char
                            ),
                        });
                    }
                    normalized_mutations.push((m.zero_based_position, ref_b, alt_b));
                }

                let fwd_variants = Self::expand_primer_variants(
                    &forward_primer,
                    self.max_fragments_per_container(),
                )?;
                let rev_variants = Self::expand_primer_variants(
                    &reverse_primer,
                    self.max_fragments_per_container(),
                )?;
                if fwd_variants.is_empty() || rev_variants.is_empty() {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: "PCR primers must not be empty".to_string(),
                    });
                }

                let require_all = require_all_mutations.unwrap_or(true);
                let mut selected: Vec<((usize, usize), String)> = vec![];
                let mut seen_amplicons: HashSet<String> = HashSet::new();
                for fwd_full in &fwd_variants {
                    let fwd_anneal_len = forward_primer.anneal_len.unwrap_or(fwd_full.len());
                    if fwd_anneal_len == 0 || fwd_anneal_len > fwd_full.len() {
                        return Err(EngineError {
                            code: ErrorCode::InvalidInput,
                            message: "Invalid anneal_len in PCR primer spec".to_string(),
                        });
                    }
                    let fwd_anneal = &fwd_full[fwd_full.len() - fwd_anneal_len..];
                    let fwd_sites = Self::find_anneal_sites(
                        template_bytes,
                        fwd_anneal.as_bytes(),
                        forward_primer.max_mismatches.unwrap_or(0),
                        forward_primer.require_3prime_exact_bases.unwrap_or(0),
                        true,
                    );
                    if fwd_sites.is_empty() {
                        continue;
                    }

                    for rev_full in &rev_variants {
                        let rev_anneal_len = reverse_primer.anneal_len.unwrap_or(rev_full.len());
                        if rev_anneal_len == 0 || rev_anneal_len > rev_full.len() {
                            return Err(EngineError {
                                code: ErrorCode::InvalidInput,
                                message: "Invalid anneal_len in PCR primer spec".to_string(),
                            });
                        }
                        let rev_anneal = &rev_full[rev_full.len() - rev_anneal_len..];
                        let rev_binding = Self::reverse_complement(rev_anneal);
                        let rev_sites = Self::find_anneal_sites(
                            template_bytes,
                            rev_binding.as_bytes(),
                            reverse_primer.max_mismatches.unwrap_or(0),
                            reverse_primer.require_3prime_exact_bases.unwrap_or(0),
                            false,
                        );
                        if rev_sites.is_empty() {
                            continue;
                        }
                        let rev_full_rc = Self::reverse_complement(rev_full);

                        for fwd_pos in &fwd_sites {
                            for rev_pos in &rev_sites {
                                let interior_start = *fwd_pos + fwd_anneal_len;
                                let interior_end = *rev_pos;
                                if interior_start > interior_end {
                                    continue;
                                }
                                let interior = &template_seq[interior_start..interior_end];
                                let amplicon = format!("{fwd_full}{interior}{rev_full_rc}");
                                let mut introduced_count = 0usize;
                                let mut valid = true;
                                for (pos, _ref_b, alt_b) in &normalized_mutations {
                                    if *pos < *fwd_pos || *pos >= (*rev_pos + rev_anneal_len) {
                                        if require_all {
                                            valid = false;
                                            break;
                                        }
                                        continue;
                                    }

                                    let observed = if *pos < (*fwd_pos + fwd_anneal_len) {
                                        let offset =
                                            fwd_full.len() - fwd_anneal_len + (*pos - *fwd_pos);
                                        fwd_full.as_bytes()[offset]
                                    } else if *pos < *rev_pos {
                                        template_bytes[*pos]
                                    } else {
                                        let offset = *pos - *rev_pos;
                                        rev_full_rc.as_bytes()[offset]
                                    };

                                    if observed == *alt_b {
                                        introduced_count += 1;
                                    } else if require_all {
                                        valid = false;
                                        break;
                                    }
                                }

                                let keep = if require_all {
                                    valid && introduced_count == normalized_mutations.len()
                                } else {
                                    introduced_count > 0
                                };
                                if keep && seen_amplicons.insert(amplicon.clone()) {
                                    selected.push(((*fwd_pos, *rev_pos), amplicon));
                                    if selected.len() > self.max_fragments_per_container() {
                                        return Err(EngineError {
                                            code: ErrorCode::InvalidInput,
                                            message: format!(
                                                "Mutagenesis PCR produced more than max_fragments_per_container={}",
                                                self.max_fragments_per_container()
                                            ),
                                        });
                                    }
                                }
                            }
                        }
                    }
                }

                if selected.is_empty() {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: if require_all {
                            "No amplicon introduced all requested mutations".to_string()
                        } else {
                            "No amplicon introduced any requested mutation".to_string()
                        },
                    });
                }
                if selected.len() > self.max_fragments_per_container() {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "Mutagenesis PCR produced {} amplicons, exceeding max_fragments_per_container={}",
                            selected.len(),
                            self.max_fragments_per_container()
                        ),
                    });
                }

                let require_unique = unique.unwrap_or(false);
                if require_unique && selected.len() != 1 {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "PCR unique=true requires exactly one amplicon, found {}",
                            selected.len()
                        ),
                    });
                }
                if output_id.is_some() && selected.len() != 1 {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message:
                            "PCR output_id can only be used when exactly one amplicon is produced"
                                .to_string(),
                    });
                }

                let default_base = format!("{template}_pcr_mut");
                for (i, ((fwd_pos, rev_pos), amplicon)) in selected.into_iter().enumerate() {
                    let mut pcr_product =
                        DNAsequence::from_sequence(&amplicon).map_err(|e| EngineError {
                            code: ErrorCode::Internal,
                            message: format!("Could not create PCR product: {e}"),
                        })?;
                    pcr_product.set_circular(false);
                    Self::prepare_sequence(&mut pcr_product);

                    let requested = if i == 0 { output_id.clone() } else { None };
                    let seq_id = if let Some(id) = requested {
                        self.unique_seq_id(&id)
                    } else if seen_amplicons.len() == 1 {
                        self.unique_seq_id(&default_base)
                    } else {
                        self.unique_seq_id(&format!("{}_{}", default_base, i + 1))
                    };
                    self.state.sequences.insert(seq_id.clone(), pcr_product);
                    self.add_lineage_node(&seq_id, SequenceOrigin::Derived, Some(&result.op_id));
                    result.created_seq_ids.push(seq_id.clone());
                    result.messages.push(format!(
                        "Mutagenesis PCR product '{}' created from fwd@{} rev@{}",
                        seq_id, fwd_pos, rev_pos
                    ));
                }
            }
            Operation::DesignPrimerPairs {
                template,
                roi_start_0based,
                roi_end_0based,
                forward,
                reverse,
                pair_constraints,
                min_amplicon_bp,
                max_amplicon_bp,
                max_tm_delta_c,
                max_pairs,
                report_id,
            } => {
                let dna = self
                    .state
                    .sequences
                    .get(&template)
                    .ok_or_else(|| EngineError {
                        code: ErrorCode::NotFound,
                        message: format!("Sequence '{template}' not found"),
                    })?
                    .clone();
                if dna.is_circular() {
                    return Err(EngineError {
                        code: ErrorCode::Unsupported,
                        message: "DesignPrimerPairs currently supports linear templates only"
                            .to_string(),
                    });
                }

                let template_seq = dna.get_forward_string().to_ascii_uppercase();
                let template_bytes = template_seq.as_bytes();
                if template_bytes.is_empty() {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: "DesignPrimerPairs requires a non-empty template sequence"
                            .to_string(),
                    });
                }
                if roi_start_0based >= roi_end_0based || roi_end_0based > template_bytes.len() {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "DesignPrimerPairs ROI {}..{} is invalid for template length {}",
                            roi_start_0based,
                            roi_end_0based,
                            template_bytes.len()
                        ),
                    });
                }
                if min_amplicon_bp == 0 {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: "DesignPrimerPairs min_amplicon_bp must be >= 1".to_string(),
                    });
                }
                if min_amplicon_bp > max_amplicon_bp {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "DesignPrimerPairs min_amplicon_bp ({min_amplicon_bp}) must be <= max_amplicon_bp ({max_amplicon_bp})"
                        ),
                    });
                }
                let max_tm_delta_c = max_tm_delta_c.unwrap_or(2.0);
                if max_tm_delta_c < 0.0 {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "DesignPrimerPairs max_tm_delta_c ({max_tm_delta_c}) must be >= 0.0"
                        ),
                    });
                }
                let max_pairs = max_pairs.unwrap_or(200);
                if max_pairs == 0 {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: "DesignPrimerPairs max_pairs must be >= 1".to_string(),
                    });
                }

                Self::validate_primer_design_side_constraints("forward", &forward)?;
                Self::validate_primer_design_side_constraints("reverse", &reverse)?;
                let forward_sequence_constraints =
                    Self::normalize_primer_side_sequence_constraints(&forward)?;
                let reverse_sequence_constraints =
                    Self::normalize_primer_side_sequence_constraints(&reverse)?;
                let pair_constraints_normalized =
                    Self::normalize_primer_pair_constraints(&pair_constraints)?;
                for (label, side) in [("forward", &forward), ("reverse", &reverse)] {
                    if let Some(location) = side.location_0based {
                        if location >= template_bytes.len() {
                            return Err(EngineError {
                                code: ErrorCode::InvalidInput,
                                message: format!(
                                    "{label}.location_0based ({location}) is outside template length {}",
                                    template_bytes.len()
                                ),
                            });
                        }
                    }
                    if let Some(start) = side.start_0based {
                        if start >= template_bytes.len() {
                            return Err(EngineError {
                                code: ErrorCode::InvalidInput,
                                message: format!(
                                    "{label}.start_0based ({start}) is outside template length {}",
                                    template_bytes.len()
                                ),
                            });
                        }
                    }
                    if let Some(end) = side.end_0based {
                        if end == 0 || end > template_bytes.len() {
                            return Err(EngineError {
                                code: ErrorCode::InvalidInput,
                                message: format!(
                                    "{label}.end_0based ({end}) must be in 1..={} for this template",
                                    template_bytes.len()
                                ),
                            });
                        }
                    }
                }
                if let Some(start) = pair_constraints.fixed_amplicon_start_0based {
                    if start >= template_bytes.len() {
                        return Err(EngineError {
                            code: ErrorCode::InvalidInput,
                            message: format!(
                                "pair_constraints.fixed_amplicon_start_0based ({start}) is outside template length {}",
                                template_bytes.len()
                            ),
                        });
                    }
                }
                if let Some(end) = pair_constraints.fixed_amplicon_end_0based_exclusive {
                    if end == 0 || end > template_bytes.len() {
                        return Err(EngineError {
                            code: ErrorCode::InvalidInput,
                            message: format!(
                                "pair_constraints.fixed_amplicon_end_0based_exclusive ({end}) must be in 1..={} for this template",
                                template_bytes.len()
                            ),
                        });
                    }
                }

                let requested_backend = self.state.parameters.primer_design_backend;
                let primer3_executable = {
                    let raw = self.state.parameters.primer3_executable.trim();
                    if raw.is_empty() { "primer3_core" } else { raw }
                };
                let mut backend = PrimerDesignBackendInfo {
                    requested: requested_backend.as_str().to_string(),
                    ..PrimerDesignBackendInfo::default()
                };
                let (pairs, rejection_summary) = match requested_backend {
                    PrimerDesignBackend::Internal => {
                        backend.used = PrimerDesignBackend::Internal.as_str().to_string();
                        Self::design_primer_pairs_internal(
                            template_bytes,
                            roi_start_0based,
                            roi_end_0based,
                            &forward,
                            &forward_sequence_constraints,
                            &reverse,
                            &reverse_sequence_constraints,
                            &pair_constraints_normalized,
                            min_amplicon_bp,
                            max_amplicon_bp,
                            max_tm_delta_c,
                            max_pairs,
                        )
                    }
                    PrimerDesignBackend::Primer3 => {
                        let (pairs, rejection_summary, version, explain, request_boulder_io) =
                            Self::design_primer_pairs_primer3(
                                &template_seq,
                                roi_start_0based,
                                roi_end_0based,
                                &forward,
                                &forward_sequence_constraints,
                                &reverse,
                                &reverse_sequence_constraints,
                                &pair_constraints_normalized,
                                min_amplicon_bp,
                                max_amplicon_bp,
                                max_tm_delta_c,
                                max_pairs,
                                primer3_executable,
                            )?;
                        backend.used = PrimerDesignBackend::Primer3.as_str().to_string();
                        backend.primer3_executable = Some(primer3_executable.to_string());
                        backend.primer3_version = version;
                        backend.primer3_explain = explain;
                        backend.primer3_request_boulder_io = Some(request_boulder_io);
                        (pairs, rejection_summary)
                    }
                    PrimerDesignBackend::Auto => {
                        match Self::design_primer_pairs_primer3(
                            &template_seq,
                            roi_start_0based,
                            roi_end_0based,
                            &forward,
                            &forward_sequence_constraints,
                            &reverse,
                            &reverse_sequence_constraints,
                            &pair_constraints_normalized,
                            min_amplicon_bp,
                            max_amplicon_bp,
                            max_tm_delta_c,
                            max_pairs,
                            primer3_executable,
                        ) {
                            Ok((
                                pairs,
                                rejection_summary,
                                version,
                                explain,
                                request_boulder_io,
                            )) => {
                                backend.used = PrimerDesignBackend::Primer3.as_str().to_string();
                                backend.primer3_executable = Some(primer3_executable.to_string());
                                backend.primer3_version = version;
                                backend.primer3_explain = explain;
                                backend.primer3_request_boulder_io = Some(request_boulder_io);
                                (pairs, rejection_summary)
                            }
                            Err(err) => {
                                backend.used = PrimerDesignBackend::Internal.as_str().to_string();
                                backend.primer3_executable = Some(primer3_executable.to_string());
                                backend.fallback_reason = Some(err.message.clone());
                                result.warnings.push(format!(
                                    "Primer3 backend unavailable in auto mode: {}. Falling back to internal primer design backend.",
                                    err.message
                                ));
                                Self::design_primer_pairs_internal(
                                    template_bytes,
                                    roi_start_0based,
                                    roi_end_0based,
                                    &forward,
                                    &forward_sequence_constraints,
                                    &reverse,
                                    &reverse_sequence_constraints,
                                    &pair_constraints_normalized,
                                    min_amplicon_bp,
                                    max_amplicon_bp,
                                    max_tm_delta_c,
                                    max_pairs,
                                )
                            }
                        }
                    }
                };

                let report_id = Self::render_primer_design_report_id(report_id, &template);
                let report = PrimerDesignReport {
                    schema: PRIMER_DESIGN_REPORT_SCHEMA.to_string(),
                    report_id: report_id.clone(),
                    template: template.clone(),
                    generated_at_unix_ms: Self::now_unix_ms(),
                    roi_start_0based,
                    roi_end_0based,
                    forward,
                    reverse,
                    pair_constraints,
                    min_amplicon_bp,
                    max_amplicon_bp,
                    max_tm_delta_c,
                    max_pairs,
                    pair_count: pairs.len(),
                    pairs,
                    rejection_summary,
                    backend,
                };
                if report.rejection_summary.pair_evaluation_limit_skipped > 0 {
                    result.warnings.push(format!(
                        "Internal primer-pair search reached its evaluation limit and skipped {} candidate combinations; narrow ROI/constraints for a more exhaustive run",
                        report.rejection_summary.pair_evaluation_limit_skipped
                    ));
                }
                let mut store = self.read_primer_design_store();
                let replaced = store
                    .reports
                    .insert(report.report_id.clone(), report.clone())
                    .is_some();
                self.write_primer_design_store(store)?;
                result.messages.push(format!(
                    "{} primer-design report '{}' for template '{}' (pairs={})",
                    if replaced { "Updated" } else { "Created" },
                    report.report_id,
                    report.template,
                    report.pair_count
                ));
                if !report.pairs.is_empty() {
                    parent_seq_ids.push(template.clone());
                    let report_token = Self::normalize_id_token(&report.report_id);
                    let mut pair_container_ids: Vec<String> =
                        Vec::with_capacity(report.pairs.len());
                    for (pair_index, pair) in report.pairs.iter().enumerate() {
                        let pair_rank = if pair.rank == 0 {
                            pair_index + 1
                        } else {
                            pair.rank
                        };
                        let pair_token = format!("r{pair_rank:02}");
                        let base = format!("{report_token}_{pair_token}");
                        let forward_seq_id = self.unique_seq_id(&format!("{base}_fwd"));
                        let reverse_seq_id = self.unique_seq_id(&format!("{base}_rev"));

                        let mut forward_primer = DNAsequence::from_sequence(&pair.forward.sequence)
                            .map_err(|e| EngineError {
                                code: ErrorCode::Internal,
                                message: format!(
                                    "Could not materialize forward primer sequence for report '{}' pair {}: {e}",
                                    report.report_id, pair_rank
                                ),
                            })?;
                        forward_primer.set_circular(false);
                        forward_primer.set_name(forward_seq_id.clone());
                        Self::prepare_sequence(&mut forward_primer);
                        self.state
                            .sequences
                            .insert(forward_seq_id.clone(), forward_primer);
                        self.add_lineage_node(
                            &forward_seq_id,
                            SequenceOrigin::Derived,
                            Some(&result.op_id),
                        );
                        result.created_seq_ids.push(forward_seq_id.clone());

                        let mut reverse_primer = DNAsequence::from_sequence(&pair.reverse.sequence)
                            .map_err(|e| EngineError {
                                code: ErrorCode::Internal,
                                message: format!(
                                    "Could not materialize reverse primer sequence for report '{}' pair {}: {e}",
                                    report.report_id, pair_rank
                                ),
                            })?;
                        reverse_primer.set_circular(false);
                        reverse_primer.set_name(reverse_seq_id.clone());
                        Self::prepare_sequence(&mut reverse_primer);
                        self.state
                            .sequences
                            .insert(reverse_seq_id.clone(), reverse_primer);
                        self.add_lineage_node(
                            &reverse_seq_id,
                            SequenceOrigin::Derived,
                            Some(&result.op_id),
                        );
                        result.created_seq_ids.push(reverse_seq_id.clone());

                        let pair_members = vec![forward_seq_id.clone(), reverse_seq_id.clone()];
                        if let Some(container_id) = self.add_container(
                            &pair_members,
                            ContainerKind::Pool,
                            Some(format!("Primer pair {} {}", report.report_id, pair_token)),
                            Some(&result.op_id),
                        ) {
                            pair_container_ids.push(container_id.clone());
                            result.messages.push(format!(
                                "Created primer-pair container '{}' for {} ({}, {})",
                                container_id, pair_token, forward_seq_id, reverse_seq_id
                            ));
                        }
                    }
                    result.messages.push(format!(
                        "Materialized {} primer sequence(s) and {} primer-pair container(s) for report '{}'",
                        report.pairs.len().saturating_mul(2),
                        pair_container_ids.len(),
                        report.report_id
                    ));
                }
                if let Some(top_pair) = report.pairs.first() {
                    let advisories = Self::primer_pair_heuristic_advisories(top_pair);
                    if !advisories.is_empty() {
                        result.warnings.push(format!(
                            "Top primer pair in report '{}' has heuristic advisories: {}",
                            report.report_id,
                            advisories.join("; ")
                        ));
                    }
                }
                if report.pair_count == 0 {
                    result.warnings.push(format!(
                        "No primer pairs satisfied constraints for report '{}'",
                        report.report_id
                    ));
                    if let Some(explain) = report.backend.primer3_explain.as_deref() {
                        result.warnings.push(format!(
                            "Primer3 explain for report '{}': {}",
                            report.report_id, explain
                        ));
                    }
                    if report.backend.primer3_request_boulder_io.is_some() {
                        result.warnings.push(format!(
                            "Primer3 request payload was captured for report '{}' and can be exported for local reruns",
                            report.report_id
                        ));
                    }
                }
            }
            Operation::DesignQpcrAssays {
                template,
                roi_start_0based,
                roi_end_0based,
                forward,
                reverse,
                probe,
                pair_constraints,
                min_amplicon_bp,
                max_amplicon_bp,
                max_tm_delta_c,
                max_probe_tm_delta_c,
                max_assays,
                report_id,
            } => {
                let dna = self
                    .state
                    .sequences
                    .get(&template)
                    .ok_or_else(|| EngineError {
                        code: ErrorCode::NotFound,
                        message: format!("Sequence '{template}' not found"),
                    })?
                    .clone();
                if dna.is_circular() {
                    return Err(EngineError {
                        code: ErrorCode::Unsupported,
                        message: "DesignQpcrAssays currently supports linear templates only"
                            .to_string(),
                    });
                }

                let template_seq = dna.get_forward_string().to_ascii_uppercase();
                let template_bytes = template_seq.as_bytes();
                if template_bytes.is_empty() {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: "DesignQpcrAssays requires a non-empty template sequence"
                            .to_string(),
                    });
                }
                if roi_start_0based >= roi_end_0based || roi_end_0based > template_bytes.len() {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "DesignQpcrAssays ROI {}..{} is invalid for template length {}",
                            roi_start_0based,
                            roi_end_0based,
                            template_bytes.len()
                        ),
                    });
                }
                if min_amplicon_bp == 0 {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: "DesignQpcrAssays min_amplicon_bp must be >= 1".to_string(),
                    });
                }
                if min_amplicon_bp > max_amplicon_bp {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "DesignQpcrAssays min_amplicon_bp ({min_amplicon_bp}) must be <= max_amplicon_bp ({max_amplicon_bp})"
                        ),
                    });
                }
                let max_tm_delta_c = max_tm_delta_c.unwrap_or(2.0);
                if max_tm_delta_c < 0.0 {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "DesignQpcrAssays max_tm_delta_c ({max_tm_delta_c}) must be >= 0.0"
                        ),
                    });
                }
                let max_probe_tm_delta_c = max_probe_tm_delta_c.unwrap_or(10.0);
                if max_probe_tm_delta_c < 0.0 {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "DesignQpcrAssays max_probe_tm_delta_c ({max_probe_tm_delta_c}) must be >= 0.0"
                        ),
                    });
                }
                let max_assays = max_assays.unwrap_or(200);
                if max_assays == 0 {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: "DesignQpcrAssays max_assays must be >= 1".to_string(),
                    });
                }

                Self::validate_primer_design_side_constraints("forward", &forward)?;
                Self::validate_primer_design_side_constraints("reverse", &reverse)?;
                Self::validate_primer_design_side_constraints("probe", &probe)?;
                let forward_sequence_constraints =
                    Self::normalize_primer_side_sequence_constraints(&forward)?;
                let reverse_sequence_constraints =
                    Self::normalize_primer_side_sequence_constraints(&reverse)?;
                let probe_sequence_constraints =
                    Self::normalize_primer_side_sequence_constraints(&probe)?;
                let pair_constraints_normalized =
                    Self::normalize_primer_pair_constraints(&pair_constraints)?;
                for (label, side) in [
                    ("forward", &forward),
                    ("reverse", &reverse),
                    ("probe", &probe),
                ] {
                    if let Some(location) = side.location_0based {
                        if location >= template_bytes.len() {
                            return Err(EngineError {
                                code: ErrorCode::InvalidInput,
                                message: format!(
                                    "{label}.location_0based ({location}) is outside template length {}",
                                    template_bytes.len()
                                ),
                            });
                        }
                    }
                    if let Some(start) = side.start_0based {
                        if start >= template_bytes.len() {
                            return Err(EngineError {
                                code: ErrorCode::InvalidInput,
                                message: format!(
                                    "{label}.start_0based ({start}) is outside template length {}",
                                    template_bytes.len()
                                ),
                            });
                        }
                    }
                    if let Some(end) = side.end_0based {
                        if end == 0 || end > template_bytes.len() {
                            return Err(EngineError {
                                code: ErrorCode::InvalidInput,
                                message: format!(
                                    "{label}.end_0based ({end}) must be in 1..={} for this template",
                                    template_bytes.len()
                                ),
                            });
                        }
                    }
                }
                if let Some(start) = pair_constraints.fixed_amplicon_start_0based {
                    if start >= template_bytes.len() {
                        return Err(EngineError {
                            code: ErrorCode::InvalidInput,
                            message: format!(
                                "pair_constraints.fixed_amplicon_start_0based ({start}) is outside template length {}",
                                template_bytes.len()
                            ),
                        });
                    }
                }
                if let Some(end) = pair_constraints.fixed_amplicon_end_0based_exclusive {
                    if end == 0 || end > template_bytes.len() {
                        return Err(EngineError {
                            code: ErrorCode::InvalidInput,
                            message: format!(
                                "pair_constraints.fixed_amplicon_end_0based_exclusive ({end}) must be in 1..={} for this template",
                                template_bytes.len()
                            ),
                        });
                    }
                }

                let pair_generation_limit = max_assays.saturating_mul(25).clamp(max_assays, 5000);
                let requested_backend = self.state.parameters.primer_design_backend;
                let primer3_executable = {
                    let raw = self.state.parameters.primer3_executable.trim();
                    if raw.is_empty() { "primer3_core" } else { raw }
                };
                let mut backend = PrimerDesignBackendInfo {
                    requested: requested_backend.as_str().to_string(),
                    ..PrimerDesignBackendInfo::default()
                };
                let (pair_candidates, pair_rejections) = match requested_backend {
                    PrimerDesignBackend::Internal => {
                        backend.used = PrimerDesignBackend::Internal.as_str().to_string();
                        Self::design_primer_pairs_internal(
                            template_bytes,
                            roi_start_0based,
                            roi_end_0based,
                            &forward,
                            &forward_sequence_constraints,
                            &reverse,
                            &reverse_sequence_constraints,
                            &pair_constraints_normalized,
                            min_amplicon_bp,
                            max_amplicon_bp,
                            max_tm_delta_c,
                            pair_generation_limit,
                        )
                    }
                    PrimerDesignBackend::Primer3 => {
                        let (pairs, rejection_summary, version, explain, request_boulder_io) =
                            Self::design_primer_pairs_primer3(
                                &template_seq,
                                roi_start_0based,
                                roi_end_0based,
                                &forward,
                                &forward_sequence_constraints,
                                &reverse,
                                &reverse_sequence_constraints,
                                &pair_constraints_normalized,
                                min_amplicon_bp,
                                max_amplicon_bp,
                                max_tm_delta_c,
                                pair_generation_limit,
                                primer3_executable,
                            )?;
                        backend.used = PrimerDesignBackend::Primer3.as_str().to_string();
                        backend.primer3_executable = Some(primer3_executable.to_string());
                        backend.primer3_version = version;
                        backend.primer3_explain = explain;
                        backend.primer3_request_boulder_io = Some(request_boulder_io);
                        (pairs, rejection_summary)
                    }
                    PrimerDesignBackend::Auto => {
                        match Self::design_primer_pairs_primer3(
                            &template_seq,
                            roi_start_0based,
                            roi_end_0based,
                            &forward,
                            &forward_sequence_constraints,
                            &reverse,
                            &reverse_sequence_constraints,
                            &pair_constraints_normalized,
                            min_amplicon_bp,
                            max_amplicon_bp,
                            max_tm_delta_c,
                            pair_generation_limit,
                            primer3_executable,
                        ) {
                            Ok((
                                pairs,
                                rejection_summary,
                                version,
                                explain,
                                request_boulder_io,
                            )) => {
                                backend.used = PrimerDesignBackend::Primer3.as_str().to_string();
                                backend.primer3_executable = Some(primer3_executable.to_string());
                                backend.primer3_version = version;
                                backend.primer3_explain = explain;
                                backend.primer3_request_boulder_io = Some(request_boulder_io);
                                (pairs, rejection_summary)
                            }
                            Err(err) => {
                                backend.used = PrimerDesignBackend::Internal.as_str().to_string();
                                backend.primer3_executable = Some(primer3_executable.to_string());
                                backend.fallback_reason = Some(err.message.clone());
                                result.warnings.push(format!(
                                    "Primer3 backend unavailable in auto mode: {}. Falling back to internal qPCR design backend.",
                                    err.message
                                ));
                                Self::design_primer_pairs_internal(
                                    template_bytes,
                                    roi_start_0based,
                                    roi_end_0based,
                                    &forward,
                                    &forward_sequence_constraints,
                                    &reverse,
                                    &reverse_sequence_constraints,
                                    &pair_constraints_normalized,
                                    min_amplicon_bp,
                                    max_amplicon_bp,
                                    max_tm_delta_c,
                                    pair_generation_limit,
                                )
                            }
                        }
                    }
                };

                let (assays, rejection_summary) = Self::design_qpcr_assays_from_pairs(
                    template_bytes,
                    roi_start_0based,
                    roi_end_0based,
                    &probe,
                    &probe_sequence_constraints,
                    max_probe_tm_delta_c,
                    max_assays,
                    pair_candidates,
                    pair_rejections,
                );

                let report_id = Self::render_primer_design_report_id(report_id, &template);
                let report = QpcrDesignReport {
                    schema: QPCR_DESIGN_REPORT_SCHEMA.to_string(),
                    report_id: report_id.clone(),
                    template: template.clone(),
                    generated_at_unix_ms: Self::now_unix_ms(),
                    roi_start_0based,
                    roi_end_0based,
                    forward,
                    reverse,
                    probe,
                    pair_constraints,
                    min_amplicon_bp,
                    max_amplicon_bp,
                    max_tm_delta_c,
                    max_probe_tm_delta_c,
                    max_assays,
                    assay_count: assays.len(),
                    assays,
                    rejection_summary,
                    backend,
                };
                if report
                    .rejection_summary
                    .primer_pair
                    .pair_evaluation_limit_skipped
                    > 0
                {
                    result.warnings.push(format!(
                        "Internal primer-pair candidate generation for qPCR reached its evaluation limit and skipped {} candidate combinations; narrow ROI/constraints for a more exhaustive run",
                        report
                            .rejection_summary
                            .primer_pair
                            .pair_evaluation_limit_skipped
                    ));
                }
                let mut store = self.read_primer_design_store();
                let replaced = store
                    .qpcr_reports
                    .insert(report.report_id.clone(), report.clone())
                    .is_some();
                self.write_primer_design_store(store)?;
                result.messages.push(format!(
                    "{} qPCR-design report '{}' for template '{}' (assays={})",
                    if replaced { "Updated" } else { "Created" },
                    report.report_id,
                    report.template,
                    report.assay_count
                ));
                if let Some(top_assay) = report.assays.first() {
                    let dimer_metrics = Self::compute_primer_pair_dimer_metrics(
                        top_assay.forward.sequence.as_bytes(),
                        top_assay.reverse.sequence.as_bytes(),
                    );
                    let pair_like = PrimerDesignPairRecord {
                        rank: top_assay.rank,
                        score: top_assay.score,
                        forward: top_assay.forward.clone(),
                        reverse: top_assay.reverse.clone(),
                        amplicon_start_0based: top_assay.amplicon_start_0based,
                        amplicon_end_0based_exclusive: top_assay.amplicon_end_0based_exclusive,
                        amplicon_length_bp: top_assay.amplicon_length_bp,
                        tm_delta_c: top_assay.primer_tm_delta_c,
                        primer_pair_complementary_run_bp: dimer_metrics.max_complementary_run_bp,
                        primer_pair_3prime_complementary_run_bp: dimer_metrics
                            .max_3prime_complementary_run_bp,
                        rule_flags: PrimerDesignPairRuleFlags {
                            roi_covered: top_assay.rule_flags.roi_covered,
                            amplicon_size_in_range: top_assay.rule_flags.amplicon_size_in_range,
                            tm_delta_in_range: top_assay.rule_flags.primer_tm_delta_in_range,
                            forward_secondary_structure_ok: top_assay
                                .forward
                                .longest_homopolymer_run_bp
                                <= PRIMER_RECOMMENDED_MAX_HOMOPOLYMER_RUN_BP
                                && top_assay.forward.self_complementary_run_bp
                                    <= PRIMER_RECOMMENDED_MAX_SELF_COMPLEMENTARY_RUN_BP,
                            reverse_secondary_structure_ok: top_assay
                                .reverse
                                .longest_homopolymer_run_bp
                                <= PRIMER_RECOMMENDED_MAX_HOMOPOLYMER_RUN_BP
                                && top_assay.reverse.self_complementary_run_bp
                                    <= PRIMER_RECOMMENDED_MAX_SELF_COMPLEMENTARY_RUN_BP,
                            primer_pair_dimer_risk_low: dimer_metrics.max_complementary_run_bp
                                <= PRIMER_RECOMMENDED_MAX_PAIR_DIMER_RUN_BP
                                && dimer_metrics.max_3prime_complementary_run_bp
                                    <= PRIMER_RECOMMENDED_MAX_PAIR_3PRIME_DIMER_RUN_BP,
                            forward_three_prime_gc_clamp: top_assay.forward.three_prime_gc_clamp,
                            reverse_three_prime_gc_clamp: top_assay.reverse.three_prime_gc_clamp,
                        },
                    };
                    let advisories = Self::primer_pair_heuristic_advisories(&pair_like);
                    if !advisories.is_empty() {
                        result.warnings.push(format!(
                            "Top qPCR primer pair in report '{}' has heuristic advisories: {}",
                            report.report_id,
                            advisories.join("; ")
                        ));
                    }
                }
                if report.assay_count == 0 {
                    result.warnings.push(format!(
                        "No qPCR assays satisfied constraints for report '{}'",
                        report.report_id
                    ));
                }
            }
            Operation::DeriveTranscriptSequences {
                seq_id,
                feature_ids,
                scope,
                output_prefix,
            } => {
                let source_sequence_upper = self
                    .state
                    .sequences
                    .get(&seq_id)
                    .ok_or_else(|| EngineError {
                        code: ErrorCode::NotFound,
                        message: format!("Sequence '{seq_id}' not found"),
                    })?
                    .get_forward_string()
                    .to_ascii_uppercase()
                    .into_bytes();
                let source_features = self
                    .state
                    .sequences
                    .get(&seq_id)
                    .ok_or_else(|| EngineError {
                        code: ErrorCode::NotFound,
                        message: format!("Sequence '{seq_id}' not found"),
                    })?
                    .features()
                    .to_vec();
                let transcript_feature_ids =
                    self.transcript_feature_ids_for_derivation(&seq_id, &feature_ids, scope)?;
                parent_seq_ids.push(seq_id.clone());
                let normalized_prefix = output_prefix
                    .as_deref()
                    .map(str::trim)
                    .filter(|value| !value.is_empty())
                    .map(|value| value.to_string())
                    .unwrap_or_else(|| format!("{seq_id}__mrna"));
                for transcript_feature_id in transcript_feature_ids {
                    let source_feature =
                        source_features
                            .get(transcript_feature_id)
                            .ok_or_else(|| EngineError {
                                code: ErrorCode::NotFound,
                                message: format!(
                                    "Feature id '{}' was not found in sequence '{}'",
                                    transcript_feature_id, seq_id
                                ),
                            })?;
                    let (derived_dna, transcript_id, transcript_label, is_reverse, exon_count) =
                        Self::derive_transcript_sequence_from_feature(
                            &source_sequence_upper,
                            source_feature,
                            transcript_feature_id,
                            &seq_id,
                        )?;
                    let transcript_token = Self::normalize_id_token(&transcript_id);
                    let transcript_token = if transcript_token.is_empty() {
                        format!("feature_{}", transcript_feature_id + 1)
                    } else {
                        transcript_token
                    };
                    let base_seq_id = format!(
                        "{normalized_prefix}__f{}__{}",
                        transcript_feature_id + 1,
                        transcript_token
                    );
                    let derived_seq_id = self.unique_seq_id(&base_seq_id);
                    self.state
                        .sequences
                        .insert(derived_seq_id.clone(), derived_dna);
                    self.add_lineage_node(
                        &derived_seq_id,
                        SequenceOrigin::Derived,
                        Some(&result.op_id),
                    );
                    result.created_seq_ids.push(derived_seq_id.clone());
                    let strand = if is_reverse { "-" } else { "+" };
                    let derived_len = self
                        .state
                        .sequences
                        .get(&derived_seq_id)
                        .map(|dna| dna.len())
                        .unwrap_or(0);
                    result.messages.push(format!(
                        "Derived transcript '{}' from mRNA feature n-{} ('{}', label='{}', strand {}, exons={}, length={} bp).",
                        derived_seq_id,
                        transcript_feature_id + 1,
                        transcript_id,
                        transcript_label,
                        strand,
                        exon_count,
                        derived_len
                    ));
                }
                if result.created_seq_ids.is_empty() {
                    return Err(EngineError {
                        code: ErrorCode::NotFound,
                        message: format!(
                            "DeriveTranscriptSequences did not produce transcripts for '{}'",
                            seq_id
                        ),
                    });
                }
            }
            Operation::ComputeDotplot {
                seq_id,
                reference_seq_id,
                span_start_0based,
                span_end_0based,
                reference_span_start_0based,
                reference_span_end_0based,
                mode,
                word_size,
                step_bp,
                max_mismatches,
                tile_bp,
                store_as,
            } => {
                let dna = self
                    .state
                    .sequences
                    .get(&seq_id)
                    .ok_or_else(|| EngineError {
                        code: ErrorCode::NotFound,
                        message: format!("Sequence '{seq_id}' not found"),
                    })?;
                let query_text = dna.get_forward_string().to_ascii_uppercase();
                let query_bytes = query_text.as_bytes();
                let (span_start_0based, span_end_0based) = Self::resolve_analysis_span(
                    query_bytes.len(),
                    span_start_0based,
                    span_end_0based,
                )?;
                let query_span = &query_bytes[span_start_0based..span_end_0based];
                let (
                    reference_label,
                    reference_seq_id_for_view,
                    reference_text,
                    reference_span_start_0based,
                    reference_span_end_0based,
                ) = match mode {
                    DotplotMode::SelfForward | DotplotMode::SelfReverseComplement => (
                        seq_id.clone(),
                        None,
                        query_text.clone(),
                        span_start_0based,
                        span_end_0based,
                    ),
                    DotplotMode::PairForward | DotplotMode::PairReverseComplement => {
                        let ref_seq_id = reference_seq_id
                            .as_ref()
                            .map(|value| value.trim())
                            .filter(|value| !value.is_empty())
                            .ok_or_else(|| EngineError {
                                code: ErrorCode::InvalidInput,
                                message: format!(
                                    "ComputeDotplot mode '{}' requires reference_seq_id",
                                    mode.as_str()
                                ),
                            })?;
                        let reference_dna =
                            self.state
                                .sequences
                                .get(ref_seq_id)
                                .ok_or_else(|| EngineError {
                                    code: ErrorCode::NotFound,
                                    message: format!(
                                        "Reference sequence '{}' not found",
                                        ref_seq_id
                                    ),
                                })?;
                        let reference_text =
                            reference_dna.get_forward_string().to_ascii_uppercase();
                        let reference_bytes = reference_text.as_bytes();
                        let (ref_start, ref_end) = Self::resolve_analysis_span(
                            reference_bytes.len(),
                            reference_span_start_0based,
                            reference_span_end_0based,
                        )?;
                        (
                            ref_seq_id.to_string(),
                            Some(ref_seq_id.to_string()),
                            reference_text,
                            ref_start,
                            ref_end,
                        )
                    }
                };
                let reference_bytes = reference_text.as_bytes();
                let reference_span =
                    &reference_bytes[reference_span_start_0based..reference_span_end_0based];
                let (points, truncated) = Self::compute_dotplot_points(
                    query_span,
                    reference_span,
                    span_start_0based,
                    reference_span_start_0based,
                    mode,
                    word_size,
                    step_bp,
                    max_mismatches,
                    MAX_DOTPLOT_POINTS,
                )?;
                let boxplot_bins = Self::compute_dotplot_boxplot_bins(
                    &points,
                    span_start_0based,
                    span_end_0based,
                    DOTPLOT_BOXPLOT_DEFAULT_BINS,
                );
                let dotplot_id = if let Some(raw_id) = store_as.as_deref() {
                    Self::normalize_analysis_id(raw_id, "dotplot")?
                } else {
                    format!("dotplot_{}", result.op_id)
                };
                let view = DotplotView {
                    schema: DOTPLOT_VIEW_SCHEMA.to_string(),
                    dotplot_id: dotplot_id.clone(),
                    seq_id: seq_id.clone(),
                    reference_seq_id: reference_seq_id_for_view,
                    generated_at_unix_ms: Self::now_unix_ms(),
                    span_start_0based,
                    span_end_0based,
                    reference_span_start_0based,
                    reference_span_end_0based,
                    mode,
                    word_size,
                    step_bp,
                    max_mismatches,
                    tile_bp,
                    point_count: points.len(),
                    points,
                    boxplot_bin_count: boxplot_bins.len(),
                    boxplot_bins,
                };
                let replaced = self
                    .read_dotplot_analysis_store()
                    .dotplots
                    .contains_key(dotplot_id.as_str());
                self.upsert_dotplot_view(view.clone())?;
                result.messages.push(format!(
                    "{} dotplot '{}' for '{}' vs '{}' (mode={}, query_span={}..{}, reference_span={}..{}, points={})",
                    if replaced { "Updated" } else { "Created" },
                    dotplot_id,
                    seq_id,
                    reference_label,
                    mode.as_str(),
                    span_start_0based,
                    span_end_0based,
                    reference_span_start_0based,
                    reference_span_end_0based,
                    view.point_count
                ));
                if truncated {
                    result.warnings.push(format!(
                        "Dotplot '{}' reached MAX_DOTPLOT_POINTS ({MAX_DOTPLOT_POINTS}); result was truncated",
                        dotplot_id
                    ));
                }
            }
            Operation::ComputeFlexibilityTrack {
                seq_id,
                span_start_0based,
                span_end_0based,
                model,
                bin_bp,
                smoothing_bp,
                store_as,
            } => {
                let dna = self
                    .state
                    .sequences
                    .get(&seq_id)
                    .ok_or_else(|| EngineError {
                        code: ErrorCode::NotFound,
                        message: format!("Sequence '{seq_id}' not found"),
                    })?;
                let seq_text = dna.get_forward_string().to_ascii_uppercase();
                let seq_bytes = seq_text.as_bytes();
                let (span_start_0based, span_end_0based) = Self::resolve_analysis_span(
                    seq_bytes.len(),
                    span_start_0based,
                    span_end_0based,
                )?;
                let span = &seq_bytes[span_start_0based..span_end_0based];
                let bins = Self::compute_flexibility_track_bins(
                    span,
                    span_start_0based,
                    model,
                    bin_bp,
                    smoothing_bp,
                )?;
                let min_score = bins
                    .iter()
                    .map(|bin| bin.score)
                    .fold(f64::INFINITY, f64::min);
                let max_score = bins
                    .iter()
                    .map(|bin| bin.score)
                    .fold(f64::NEG_INFINITY, f64::max);
                let (min_score, max_score) = if bins.is_empty() {
                    (0.0, 0.0)
                } else {
                    (min_score, max_score)
                };
                let track_id = if let Some(raw_id) = store_as.as_deref() {
                    Self::normalize_analysis_id(raw_id, "track")?
                } else {
                    format!("flex_{}", result.op_id)
                };
                let track = FlexibilityTrack {
                    schema: FLEXIBILITY_TRACK_SCHEMA.to_string(),
                    track_id: track_id.clone(),
                    seq_id: seq_id.clone(),
                    generated_at_unix_ms: Self::now_unix_ms(),
                    span_start_0based,
                    span_end_0based,
                    model,
                    bin_bp,
                    smoothing_bp,
                    min_score,
                    max_score,
                    bins,
                };
                let replaced = self
                    .read_dotplot_analysis_store()
                    .flexibility_tracks
                    .contains_key(track_id.as_str());
                self.upsert_flexibility_track(track.clone())?;
                result.messages.push(format!(
                    "{} flexibility track '{}' for '{}' (model={}, span={}..{}, bins={})",
                    if replaced { "Updated" } else { "Created" },
                    track_id,
                    seq_id,
                    model.as_str(),
                    span_start_0based,
                    span_end_0based,
                    track.bins.len()
                ));
            }
            Operation::DeriveSplicingReferences {
                seq_id,
                span_start_0based,
                span_end_0based,
                seed_feature_id,
                scope,
                output_prefix,
            } => {
                parent_seq_ids.push(seq_id.clone());
                let dna = self
                    .state
                    .sequences
                    .get(&seq_id)
                    .ok_or_else(|| EngineError {
                        code: ErrorCode::NotFound,
                        message: format!("Sequence '{seq_id}' not found"),
                    })?
                    .clone();
                let forward = dna.forward_bytes();
                let (span_start_0based, span_end_0based) = Self::resolve_analysis_span(
                    forward.len(),
                    Some(span_start_0based),
                    Some(span_end_0based),
                )?;
                let resolved_seed_feature_id = if let Some(feature_id) = seed_feature_id {
                    feature_id
                } else {
                    let mut best: Option<(usize, usize, usize)> = None;
                    for (feature_idx, feature) in dna.features().iter().enumerate() {
                        if !Self::is_mrna_feature(feature) {
                            continue;
                        }
                        let mut ranges = vec![];
                        collect_location_ranges_usize(&feature.location, &mut ranges);
                        if ranges.is_empty() {
                            if let Ok((from, to)) = feature.location.find_bounds() {
                                if from >= 0 && to >= 0 {
                                    ranges.push((from as usize, to as usize));
                                }
                            }
                        }
                        if ranges.is_empty() {
                            continue;
                        }
                        let mut overlap_bp = 0usize;
                        let mut min_start = usize::MAX;
                        for range in ranges {
                            min_start = min_start.min(range.0);
                            if let Some((overlap_start, overlap_end)) =
                                Self::range_intersection_0based(
                                    range,
                                    (span_start_0based, span_end_0based),
                                )
                            {
                                overlap_bp = overlap_bp
                                    .saturating_add(overlap_end.saturating_sub(overlap_start));
                            }
                        }
                        if overlap_bp == 0 {
                            continue;
                        }
                        match best {
                            None => best = Some((overlap_bp, min_start, feature_idx)),
                            Some((best_overlap, best_start, best_idx)) => {
                                if overlap_bp > best_overlap
                                    || (overlap_bp == best_overlap
                                        && (min_start < best_start
                                            || (min_start == best_start && feature_idx < best_idx)))
                                {
                                    best = Some((overlap_bp, min_start, feature_idx));
                                }
                            }
                        }
                    }
                    best.map(|(_, _, feature_idx)| feature_idx)
                        .ok_or_else(|| EngineError {
                            code: ErrorCode::NotFound,
                            message: format!(
                                "DeriveSplicingReferences requires a seed mRNA feature id or at least one mRNA feature overlapping span {}..{} in '{}'",
                                span_start_0based, span_end_0based, seq_id
                            ),
                        })?
                };
                let splicing =
                    self.build_splicing_expert_view(&seq_id, resolved_seed_feature_id, scope)?;
                let prefix = output_prefix
                    .as_deref()
                    .map(str::trim)
                    .filter(|value| !value.is_empty())
                    .map(ToString::to_string)
                    .unwrap_or_else(|| format!("{seq_id}_splicing_refs"));

                let mut dna_window = dna
                    .extract_region_preserving_features(span_start_0based, span_end_0based)
                    .ok_or_else(|| EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "Could not derive DNA window {}..{} from '{}'",
                            span_start_0based, span_end_0based, seq_id
                        ),
                    })?;
                if dna_window.is_empty() {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "Could not derive DNA window {}..{} from '{}'",
                            span_start_0based, span_end_0based, seq_id
                        ),
                    });
                }
                dna_window.set_circular(false);
                let dna_window_seq_id = self.unique_seq_id(&format!("{prefix}_dna"));
                dna_window.set_name(dna_window_seq_id.clone());
                Self::prepare_sequence(&mut dna_window);
                self.state
                    .sequences
                    .insert(dna_window_seq_id.clone(), dna_window);
                self.add_lineage_node(
                    &dna_window_seq_id,
                    SequenceOrigin::Derived,
                    Some(&result.op_id),
                );
                result.created_seq_ids.push(dna_window_seq_id.clone());
                result.messages.push(format!(
                    "Derived DNA window '{}' from '{}' span {}..{}",
                    dna_window_seq_id, seq_id, span_start_0based, span_end_0based
                ));
                if let Some(container_id) = self.add_container(
                    std::slice::from_ref(&dna_window_seq_id),
                    ContainerKind::Singleton,
                    Some("Derived splicing DNA window".to_string()),
                    Some(&result.op_id),
                ) {
                    result
                        .messages
                        .push(format!("Created DNA-window container '{}'", container_id));
                }

                let mut mrna_seq_ids: Vec<SeqId> = vec![];
                for (lane_index, lane) in splicing.transcripts.iter().enumerate() {
                    let template = Self::make_transcript_template(&dna, lane, 0);
                    if template.sequence.is_empty() {
                        continue;
                    }
                    let rna_bytes = template
                        .sequence
                        .iter()
                        .map(|base| match base.to_ascii_uppercase() {
                            b'T' => b'U',
                            b'A' | b'C' | b'G' | b'U' | b'N' => base.to_ascii_uppercase(),
                            _ => b'N',
                        })
                        .collect::<Vec<_>>();
                    let rna_sequence = String::from_utf8(rna_bytes).map_err(|e| EngineError {
                        code: ErrorCode::Internal,
                        message: format!("Could not render mRNA sequence bytes: {e}"),
                    })?;
                    let mut mrna =
                        DNAsequence::from_sequence(&rna_sequence).map_err(|e| EngineError {
                            code: ErrorCode::Internal,
                            message: format!(
                                "Could not create mRNA sequence '{}' from transcript '{}': {e}",
                                lane.transcript_id, lane.label
                            ),
                        })?;
                    mrna.set_circular(false);
                    let transcript_token = lane
                        .transcript_id
                        .chars()
                        .map(|ch| {
                            if ch.is_ascii_alphanumeric() || matches!(ch, '_' | '-' | '.') {
                                ch
                            } else {
                                '_'
                            }
                        })
                        .collect::<String>()
                        .trim_matches('_')
                        .to_string();
                    let transcript_token = if transcript_token.is_empty() {
                        format!("tx{}", lane_index + 1)
                    } else {
                        transcript_token
                    };
                    let mrna_seq_id =
                        self.unique_seq_id(&format!("{prefix}_mrna_{transcript_token}"));
                    mrna.set_name(mrna_seq_id.clone());
                    Self::prepare_sequence(&mut mrna);
                    self.state.sequences.insert(mrna_seq_id.clone(), mrna);
                    self.add_lineage_node(
                        &mrna_seq_id,
                        SequenceOrigin::Derived,
                        Some(&result.op_id),
                    );
                    result.created_seq_ids.push(mrna_seq_id.clone());
                    mrna_seq_ids.push(mrna_seq_id);
                }
                if mrna_seq_ids.is_empty() {
                    result.warnings.push(format!(
                        "No transcript-derived mRNA sequences were generated for '{}' (seed feature {}, scope={})",
                        seq_id,
                        resolved_seed_feature_id,
                        scope.as_str()
                    ));
                } else {
                    result.messages.push(format!(
                        "Derived {} mRNA sequence(s): {}",
                        mrna_seq_ids.len(),
                        mrna_seq_ids.join(", ")
                    ));
                    if let Some(container_id) = self.add_container(
                        &mrna_seq_ids,
                        ContainerKind::Pool,
                        Some("Derived splicing mRNA isoforms".to_string()),
                        Some(&result.op_id),
                    ) {
                        result
                            .messages
                            .push(format!("Created mRNA-isoform container '{}'", container_id));
                    }
                }

                let mut unique_exons: Vec<(usize, usize)> = splicing
                    .unique_exons
                    .iter()
                    .map(|exon| (exon.start_1based.saturating_sub(1), exon.end_1based))
                    .collect();
                if splicing.strand.trim() == "-" {
                    unique_exons.reverse();
                }
                let mut exon_reference_bytes: Vec<u8> = vec![];
                for (start_0based, end_1based) in unique_exons {
                    if start_0based >= forward.len() {
                        continue;
                    }
                    let end_0based_exclusive = end_1based.min(forward.len());
                    if end_0based_exclusive <= start_0based {
                        continue;
                    }
                    let segment = &forward[start_0based..end_0based_exclusive];
                    if splicing.strand.trim() == "-" {
                        exon_reference_bytes.extend(
                            Self::reverse_complement_bytes(segment)
                                .into_iter()
                                .map(|base| base.to_ascii_uppercase()),
                        );
                    } else {
                        exon_reference_bytes
                            .extend(segment.iter().copied().map(Self::normalize_nucleotide_base));
                    }
                }
                if exon_reference_bytes.is_empty() {
                    result.warnings.push(format!(
                        "Could not derive an exon-reference sequence for '{}' (seed feature {}, scope={})",
                        seq_id,
                        resolved_seed_feature_id,
                        scope.as_str()
                    ));
                } else {
                    let exon_reference_sequence =
                        String::from_utf8(exon_reference_bytes).map_err(|e| EngineError {
                            code: ErrorCode::Internal,
                            message: format!("Could not render exon-reference sequence bytes: {e}"),
                        })?;
                    let mut exon_reference = DNAsequence::from_sequence(&exon_reference_sequence)
                        .map_err(|e| EngineError {
                        code: ErrorCode::Internal,
                        message: format!(
                            "Could not create exon-reference sequence from '{}' transcripts: {e}",
                            seq_id
                        ),
                    })?;
                    exon_reference.set_circular(false);
                    let exon_reference_seq_id =
                        self.unique_seq_id(&format!("{prefix}_exon_reference"));
                    exon_reference.set_name(exon_reference_seq_id.clone());
                    Self::prepare_sequence(&mut exon_reference);
                    self.state
                        .sequences
                        .insert(exon_reference_seq_id.clone(), exon_reference);
                    self.add_lineage_node(
                        &exon_reference_seq_id,
                        SequenceOrigin::Derived,
                        Some(&result.op_id),
                    );
                    result.created_seq_ids.push(exon_reference_seq_id.clone());
                    result.messages.push(format!(
                        "Derived exon-reference sequence '{}' (seed feature {}, scope={}, unique_exons={})",
                        exon_reference_seq_id,
                        resolved_seed_feature_id,
                        scope.as_str(),
                        splicing.unique_exons.len()
                    ));
                    if let Some(container_id) = self.add_container(
                        std::slice::from_ref(&exon_reference_seq_id),
                        ContainerKind::Singleton,
                        Some("Derived exon-reference sequence".to_string()),
                        Some(&result.op_id),
                    ) {
                        result.messages.push(format!(
                            "Created exon-reference container '{}'",
                            container_id
                        ));
                    }
                }
            }
            Operation::AlignSequences {
                query_seq_id,
                target_seq_id,
                query_span_start_0based,
                query_span_end_0based,
                target_span_start_0based,
                target_span_end_0based,
                mode,
                match_score,
                mismatch_score,
                gap_open,
                gap_extend,
            } => {
                if match_score <= 0 {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: "AlignSequences requires match_score > 0".to_string(),
                    });
                }
                if gap_open > 0 || gap_extend > 0 {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: "AlignSequences requires gap_open and gap_extend to be <= 0"
                            .to_string(),
                    });
                }
                let query_dna =
                    self.state
                        .sequences
                        .get(&query_seq_id)
                        .ok_or_else(|| EngineError {
                            code: ErrorCode::NotFound,
                            message: format!("Sequence '{query_seq_id}' not found"),
                        })?;
                let target_dna =
                    self.state
                        .sequences
                        .get(&target_seq_id)
                        .ok_or_else(|| EngineError {
                            code: ErrorCode::NotFound,
                            message: format!("Sequence '{target_seq_id}' not found"),
                        })?;
                let query_text = query_dna.get_forward_string().to_ascii_uppercase();
                let target_text = target_dna.get_forward_string().to_ascii_uppercase();
                let query_bytes = query_text.as_bytes();
                let target_bytes = target_text.as_bytes();
                let (query_span_start_0based, query_span_end_0based) = Self::resolve_analysis_span(
                    query_bytes.len(),
                    query_span_start_0based,
                    query_span_end_0based,
                )?;
                let (target_span_start_0based, target_span_end_0based) =
                    Self::resolve_analysis_span(
                        target_bytes.len(),
                        target_span_start_0based,
                        target_span_end_0based,
                    )?;
                let query_span = &query_bytes[query_span_start_0based..query_span_end_0based];
                let target_span = &target_bytes[target_span_start_0based..target_span_end_0based];

                let score = |a: u8, b: u8| {
                    if a.eq_ignore_ascii_case(&b) {
                        match_score
                    } else {
                        mismatch_score
                    }
                };
                let mut aligner =
                    bio::alignment::pairwise::Aligner::new(gap_open, gap_extend, &score);
                let alignment = match mode {
                    PairwiseAlignmentMode::Global => aligner.global(query_span, target_span),
                    PairwiseAlignmentMode::Local => aligner.local(query_span, target_span),
                };

                let mut matches = 0usize;
                let mut mismatches = 0usize;
                let mut insertions = 0usize;
                let mut deletions = 0usize;
                let mut cigar = String::new();
                let mut run_len = 0usize;
                let mut run_code = 'M';
                for op in &alignment.operations {
                    let code = match op {
                        bio::alignment::AlignmentOperation::Match => {
                            matches += 1;
                            Some('=')
                        }
                        bio::alignment::AlignmentOperation::Subst => {
                            mismatches += 1;
                            Some('X')
                        }
                        bio::alignment::AlignmentOperation::Ins => {
                            insertions += 1;
                            Some('I')
                        }
                        bio::alignment::AlignmentOperation::Del => {
                            deletions += 1;
                            Some('D')
                        }
                        bio::alignment::AlignmentOperation::Xclip(_) => Some('S'),
                        bio::alignment::AlignmentOperation::Yclip(_) => None,
                    };
                    let Some(code) = code else {
                        continue;
                    };
                    if run_len == 0 {
                        run_code = code;
                        run_len = 1;
                    } else if run_code == code {
                        run_len += 1;
                    } else {
                        cigar.push_str(&format!("{run_len}{run_code}"));
                        run_code = code;
                        run_len = 1;
                    }
                }
                if run_len > 0 {
                    cigar.push_str(&format!("{run_len}{run_code}"));
                }
                if cigar.is_empty() {
                    cigar = "*".to_string();
                }

                let aligned_columns = matches + mismatches + insertions + deletions;
                let query_covered = alignment.xend.saturating_sub(alignment.xstart);
                let target_covered = alignment.yend.saturating_sub(alignment.ystart);
                let identity_fraction = if aligned_columns == 0 {
                    0.0
                } else {
                    matches as f64 / aligned_columns as f64
                };
                let query_coverage_fraction = if query_span.is_empty() {
                    0.0
                } else {
                    query_covered as f64 / query_span.len() as f64
                };
                let target_coverage_fraction = if target_span.is_empty() {
                    0.0
                } else {
                    target_covered as f64 / target_span.len() as f64
                };
                let report = SequenceAlignmentReport {
                    schema: SEQUENCE_ALIGNMENT_REPORT_SCHEMA.to_string(),
                    mode,
                    query_seq_id: query_seq_id.clone(),
                    target_seq_id: target_seq_id.clone(),
                    query_span_start_0based,
                    query_span_end_0based,
                    target_span_start_0based,
                    target_span_end_0based,
                    aligned_query_start_0based: query_span_start_0based + alignment.xstart,
                    aligned_query_end_0based_exclusive: query_span_start_0based + alignment.xend,
                    aligned_target_start_0based: target_span_start_0based + alignment.ystart,
                    aligned_target_end_0based_exclusive: target_span_start_0based + alignment.yend,
                    score: alignment.score,
                    match_score,
                    mismatch_score,
                    gap_open,
                    gap_extend,
                    aligned_columns,
                    matches,
                    mismatches,
                    insertions,
                    deletions,
                    identity_fraction,
                    query_coverage_fraction,
                    target_coverage_fraction,
                    cigar,
                };
                result.sequence_alignment = Some(report.clone());
                result.messages.push(format!(
                    "Computed {} alignment '{}' vs '{}' (score={}, identity={:.3}, query_cov={:.3}, target_cov={:.3}, cigar={})",
                    mode.as_str(),
                    query_seq_id,
                    target_seq_id,
                    report.score,
                    report.identity_fraction,
                    report.query_coverage_fraction,
                    report.target_coverage_fraction,
                    report.cigar
                ));
            }
            Operation::InterpretRnaReads {
                seq_id,
                seed_feature_id,
                profile,
                input_path,
                input_format,
                scope,
                origin_mode,
                target_gene_ids,
                roi_seed_capture_enabled,
                seed_filter,
                align_config,
                report_id,
                report_mode,
                checkpoint_path,
                checkpoint_every_reads,
                resume_from_checkpoint,
            } => {
                let mut keep_running = || true;
                let report = self.compute_rna_read_report_with_options_and_progress_and_cancel(
                    &seq_id,
                    seed_feature_id,
                    profile,
                    &input_path,
                    input_format,
                    scope,
                    origin_mode,
                    &target_gene_ids,
                    roi_seed_capture_enabled,
                    &seed_filter,
                    &align_config,
                    report_id.as_deref(),
                    &RnaReadInterpretOptions {
                        report_mode,
                        checkpoint_path: checkpoint_path.clone(),
                        checkpoint_every_reads,
                        resume_from_checkpoint,
                    },
                    on_progress,
                    &mut keep_running,
                )?;
                self.push_rna_read_report_result_message(report, &mut result)?;
            }
            Operation::AlignRnaReadReport {
                report_id,
                selection,
                align_config_override,
                selected_record_indices,
            } => {
                let mut keep_running = || true;
                let report = self.align_rna_read_report_with_progress_and_cancel(
                    &report_id,
                    selection,
                    align_config_override.clone(),
                    &selected_record_indices,
                    on_progress,
                    &mut keep_running,
                )?;
                self.push_rna_read_report_result_message(report.clone(), &mut result)?;
                result.messages.push(format!(
                    "Alignment phase updated report '{}' (selection={}{} aligned={}, msa_eligible(retained)={})",
                    report.report_id,
                    selection.as_str(),
                    if selected_record_indices.is_empty() {
                        ", ".to_string()
                    } else {
                        format!(
                            ", selected_record_indices={}, ",
                            selected_record_indices.len()
                        )
                    },
                    report.read_count_aligned,
                    report.retained_count_msa_eligible
                ));
            }
            Operation::ListRnaReadReports { seq_id } => {
                let rows = self.list_rna_read_reports(seq_id.as_deref());
                result.messages.push(format!(
                    "RNA-read reports: {} row(s){}",
                    rows.len(),
                    seq_id
                        .as_deref()
                        .map(|s| format!(" (seq_id='{}')", s))
                        .unwrap_or_default()
                ));
                for row in rows.iter().take(8) {
                    result.messages.push(format!(
                        "  - {}",
                        Self::format_rna_read_report_summary_row(row)
                    ));
                }
                if rows.len() > 8 {
                    result.messages.push(format!(
                        "  ... {} additional report row(s) omitted",
                        rows.len() - 8
                    ));
                }
            }
            Operation::ShowRnaReadReport { report_id } => {
                let report = self.get_rna_read_report(&report_id)?;
                result.messages.push(format!(
                    "RNA-read report summary: {}",
                    Self::format_rna_read_report_detail_summary(&report)
                ));
                if !report.target_gene_ids.is_empty() {
                    let preview = report
                        .target_gene_ids
                        .iter()
                        .take(8)
                        .cloned()
                        .collect::<Vec<_>>()
                        .join(",");
                    let suffix = if report.target_gene_ids.len() > 8 {
                        format!(" (+{} more)", report.target_gene_ids.len() - 8)
                    } else {
                        String::new()
                    };
                    result
                        .messages
                        .push(format!("  target_genes={}{}", preview, suffix));
                }
                if !report.origin_class_counts.is_empty() {
                    let class_summary = report
                        .origin_class_counts
                        .iter()
                        .map(|(class, count)| format!("{class}:{count}"))
                        .collect::<Vec<_>>()
                        .join(",");
                    result
                        .messages
                        .push(format!("  origin_classes={class_summary}"));
                }
            }
            Operation::ExportRnaReadReport { report_id, path } => {
                let report = self.export_rna_read_report(&report_id, &path)?;
                result.messages.push(format!(
                    "Exported RNA-read report '{}' to '{}'",
                    report.report_id, path
                ));
            }
            Operation::ExportRnaReadHitsFasta {
                report_id,
                path,
                selection,
                selected_record_indices,
            } => {
                let count = self.export_rna_read_hits_fasta(
                    &report_id,
                    &path,
                    selection,
                    &selected_record_indices,
                )?;
                result.messages.push(format!(
                    "Exported {} RNA-read hit sequence(s) from '{}' to '{}' (selection={}, selected_record_indices={})",
                    count,
                    report_id,
                    path,
                    selection.as_str(),
                    selected_record_indices.len()
                ));
            }
            Operation::ExportRnaReadSampleSheet {
                path,
                seq_id,
                report_ids,
                append,
            } => {
                let export = self.export_rna_read_sample_sheet(
                    &path,
                    seq_id.as_deref(),
                    &report_ids,
                    append,
                )?;
                result.messages.push(format!(
                    "Exported RNA-read sample sheet '{}' with {} report row(s) (append={})",
                    export.path, export.report_count, export.appended
                ));
            }
            Operation::ExportRnaReadExonPathsTsv {
                report_id,
                path,
                selection,
            } => {
                let export = self.export_rna_read_exon_paths_tsv(&report_id, &path, selection)?;
                result.messages.push(format!(
                    "Exported RNA-read exon paths '{}' to '{}' (selection={}, rows={})",
                    export.report_id,
                    export.path,
                    export.selection.as_str(),
                    export.row_count
                ));
            }
            Operation::ExportRnaReadExonAbundanceTsv {
                report_id,
                path,
                selection,
            } => {
                let export =
                    self.export_rna_read_exon_abundance_tsv(&report_id, &path, selection)?;
                result.messages.push(format!(
                    "Exported RNA-read exon abundance '{}' to '{}' (selection={}, selected_reads={}, exon_rows={}, transition_rows={})",
                    export.report_id,
                    export.path,
                    export.selection.as_str(),
                    export.selected_read_count,
                    export.exon_row_count,
                    export.transition_row_count
                ));
            }
            Operation::ExportRnaReadScoreDensitySvg {
                report_id,
                path,
                scale,
            } => {
                let export = self.export_rna_read_score_density_svg(&report_id, &path, scale)?;
                result.messages.push(format!(
                    "Exported RNA-read score-density SVG '{}' to '{}' (scale={}, bins={}, max_bin_count={})",
                    export.report_id,
                    export.path,
                    export.scale.as_str(),
                    export.bin_count,
                    export.max_bin_count
                ));
                if export.derived_from_report_hits_only {
                    result.warnings.push(format!(
                        "Report '{}' did not persist score_density_bins; derived bins from retained hits only",
                        export.report_id
                    ));
                }
            }
            Operation::ExportRnaReadAlignmentsTsv {
                report_id,
                path,
                selection,
                limit,
                selected_record_indices,
            } => {
                let export = self.export_rna_read_alignments_tsv(
                    &report_id,
                    &path,
                    selection,
                    limit,
                    &selected_record_indices,
                )?;
                let limit_text = export
                    .limit
                    .map(|value| value.to_string())
                    .unwrap_or_else(|| "all".to_string());
                result.messages.push(format!(
                    "Exported RNA-read alignment TSV '{}' to '{}' (selection={}, rows={}, aligned_total={}, limit={}, selected_record_indices={})",
                    export.report_id,
                    export.path,
                    export.selection.as_str(),
                    export.row_count,
                    export.aligned_count,
                    limit_text,
                    selected_record_indices.len()
                ));
            }
            Operation::ExportRnaReadAlignmentDotplotSvg {
                report_id,
                path,
                selection,
                max_points,
            } => {
                let export = self.export_rna_read_alignment_dotplot_svg(
                    &report_id, &path, selection, max_points,
                )?;
                result.messages.push(format!(
                    "Exported RNA-read alignment dotplot SVG '{}' to '{}' (selection={}, rendered_points={}, total_points={}, max_points={})",
                    export.report_id,
                    export.path,
                    export.selection.as_str(),
                    export.rendered_point_count,
                    export.point_count,
                    export.max_points
                ));
            }
            Operation::MaterializeRnaReadHitSequences {
                report_id,
                selection,
                selected_record_indices,
                output_prefix,
            } => {
                let report = self.get_rna_read_report(&report_id)?;
                parent_seq_ids.push(report.seq_id.clone());
                let mut explicit_record_indices = selected_record_indices.clone();
                explicit_record_indices.sort_unstable();
                explicit_record_indices.dedup();
                let explicit_selection = !explicit_record_indices.is_empty();
                let prefix = output_prefix
                    .as_deref()
                    .map(Self::normalize_id_token)
                    .filter(|value| !value.is_empty())
                    .unwrap_or_else(|| {
                        format!(
                            "{}_rna_read_{}",
                            Self::normalize_id_token(&report.seq_id),
                            Self::normalize_id_token(&report.report_id)
                        )
                    });
                let selected_hits = if explicit_selection {
                    report
                        .hits
                        .iter()
                        .filter(|hit| {
                            explicit_record_indices
                                .binary_search(&hit.record_index)
                                .is_ok()
                        })
                        .collect::<Vec<_>>()
                } else {
                    report
                        .hits
                        .iter()
                        .filter(|hit| Self::include_rna_read_hit_by_selection(hit, selection))
                        .collect::<Vec<_>>()
                };
                if selected_hits.is_empty() {
                    return Err(EngineError {
                        code: ErrorCode::NotFound,
                        message: if explicit_selection {
                            format!(
                                "MaterializeRnaReadHitSequences did not find any hits for report '{}' and record_indices={:?}",
                                report.report_id, explicit_record_indices
                            )
                        } else {
                            format!(
                                "MaterializeRnaReadHitSequences found no hits for report '{}' with selection '{}'",
                                report.report_id,
                                selection.as_str()
                            )
                        },
                    });
                }
                let mut created = Vec::<String>::with_capacity(selected_hits.len());
                for hit in selected_hits {
                    let header_token = {
                        let normalized = Self::normalize_id_token(&hit.header_id);
                        if normalized.is_empty() {
                            "read".to_string()
                        } else {
                            normalized
                        }
                    };
                    let base = format!("{}_r{}_{}", prefix, hit.record_index + 1, header_token);
                    let seq_id = self.unique_seq_id(&base);
                    let mut dna =
                        DNAsequence::from_sequence(&hit.sequence).map_err(|e| EngineError {
                            code: ErrorCode::Internal,
                            message: format!(
                                "Could not create RNA-read hit sequence '{}' from report '{}': {e}",
                                hit.header_id, report.report_id
                            ),
                        })?;
                    dna.set_name(seq_id.clone());
                    dna.set_circular(false);
                    Self::prepare_sequence(&mut dna);
                    self.state.sequences.insert(seq_id.clone(), dna);
                    self.add_lineage_node(
                        &seq_id,
                        SequenceOrigin::InSilicoSelection,
                        Some(&result.op_id),
                    );
                    result.created_seq_ids.push(seq_id.clone());
                    created.push(seq_id);
                }
                if explicit_selection && created.len() < explicit_record_indices.len() {
                    result.warnings.push(format!(
                        "Requested {} explicit record_index values but materialized {} retained hits from report '{}'",
                        explicit_record_indices.len(),
                        created.len(),
                        report.report_id
                    ));
                }
                let preview = created
                    .iter()
                    .take(4)
                    .cloned()
                    .collect::<Vec<_>>()
                    .join(", ");
                let suffix = if created.len() > 4 {
                    format!(" (+{} more)", created.len() - 4)
                } else {
                    String::new()
                };
                result.messages.push(format!(
                    "Materialized {} RNA-read hit sequence(s) from report '{}' into {}{}",
                    created.len(),
                    report.report_id,
                    preview,
                    suffix
                ));
            }
            Operation::ExtractRegion {
                input,
                from,
                to,
                output_id,
            } => {
                parent_seq_ids.push(input.clone());
                let dna = self
                    .state
                    .sequences
                    .get(&input)
                    .ok_or_else(|| EngineError {
                        code: ErrorCode::NotFound,
                        message: format!("Sequence '{input}' not found"),
                    })?
                    .clone();

                if from == to {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: "ExtractRegion requires from != to".to_string(),
                    });
                }
                let mut out = dna
                    .extract_region_preserving_features(from, to)
                    .ok_or_else(|| EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "Could not extract region {}..{} from sequence '{}'",
                            from, to, input
                        ),
                    })?;
                if out.is_empty() {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "Could not extract region {}..{} from sequence '{}'",
                            from, to, input
                        ),
                    });
                }
                out.set_circular(false);
                Self::prepare_sequence(&mut out);

                let base = output_id.unwrap_or_else(|| format!("{input}_region"));
                let seq_id = self.unique_seq_id(&base);
                self.state.sequences.insert(seq_id.clone(), out);
                self.add_lineage_node(&seq_id, SequenceOrigin::Derived, Some(&result.op_id));
                result.created_seq_ids.push(seq_id.clone());
                result.messages.push(format!(
                    "Extracted region {}..{} from '{}' into '{}'",
                    from, to, input, seq_id
                ));
            }
            Operation::ExtractAnchoredRegion {
                input,
                anchor,
                direction,
                target_length_bp,
                length_tolerance_bp,
                required_re_sites,
                required_tf_motifs,
                forward_primer,
                reverse_primer,
                output_prefix,
                unique,
                max_candidates,
            } => {
                parent_seq_ids.push(input.clone());
                let dna = self
                    .state
                    .sequences
                    .get(&input)
                    .ok_or_else(|| EngineError {
                        code: ErrorCode::NotFound,
                        message: format!("Sequence '{input}' not found"),
                    })?
                    .clone();

                if target_length_bp == 0 {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: "ExtractAnchoredRegion requires target_length_bp >= 1".to_string(),
                    });
                }

                let anchor_pos = Self::resolve_sequence_anchor_position(&dna, &anchor, "anchor")?;
                let min_len = target_length_bp.saturating_sub(length_tolerance_bp).max(1);
                let max_len = target_length_bp
                    .saturating_add(length_tolerance_bp)
                    .max(min_len);

                let forward_primer = match forward_primer {
                    Some(p) => {
                        let v = Self::normalize_iupac_text(&p)?;
                        if v.is_empty() { None } else { Some(v) }
                    }
                    None => None,
                };
                let reverse_primer = match reverse_primer {
                    Some(p) => {
                        let v = Self::normalize_iupac_text(&p)?;
                        if v.is_empty() { None } else { Some(v) }
                    }
                    None => None,
                };
                let reverse_primer_rc = reverse_primer
                    .as_ref()
                    .map(|p| Self::reverse_complement_iupac(p))
                    .transpose()?;

                let mut tf_motifs = Vec::new();
                for motif in required_tf_motifs {
                    let m = Self::resolve_tf_motif_or_iupac(&motif)?;
                    if !m.is_empty() {
                        tf_motifs.push(m);
                    }
                }

                let (required_enzymes, missing_enzymes) = if required_re_sites.is_empty() {
                    (Vec::new(), Vec::new())
                } else {
                    self.resolve_enzymes(&required_re_sites)?
                };
                if !missing_enzymes.is_empty() {
                    result.warnings.push(format!(
                        "Unknown anchored-region enzymes ignored: {}",
                        missing_enzymes.join(",")
                    ));
                }

                #[derive(Clone)]
                struct Candidate {
                    start: usize,
                    end: usize,
                    len: usize,
                    sequence: String,
                    score: usize,
                }

                let mut candidates = Vec::new();
                for len in min_len..=max_len {
                    let Some((start, end)) = Self::anchored_range(
                        anchor_pos,
                        len,
                        &direction,
                        dna.len(),
                        dna.is_circular(),
                    ) else {
                        continue;
                    };
                    let Some(fragment) = dna.get_range_safe(start..end) else {
                        continue;
                    };
                    if fragment.is_empty() {
                        continue;
                    }
                    let fragment_text = String::from_utf8_lossy(&fragment).to_string();
                    let fragment_bytes = fragment_text.as_bytes();

                    if let Some(fwd) = &forward_primer {
                        if !Self::iupac_match_at(fragment_bytes, fwd.as_bytes(), 0) {
                            continue;
                        }
                    }
                    if let Some(rev_rc) = &reverse_primer_rc {
                        if rev_rc.len() > fragment_bytes.len() {
                            continue;
                        }
                        let start_idx = fragment_bytes.len() - rev_rc.len();
                        if !Self::iupac_match_at(fragment_bytes, rev_rc.as_bytes(), start_idx) {
                            continue;
                        }
                    }

                    let mut motif_ok = true;
                    for motif in &tf_motifs {
                        if !Self::contains_motif_any_strand(fragment_bytes, motif)? {
                            motif_ok = false;
                            break;
                        }
                    }
                    if !motif_ok {
                        continue;
                    }

                    if !required_enzymes.is_empty() {
                        let frag_dna = DNAsequence::from_sequence(&fragment_text).map_err(|e| {
                            EngineError {
                                code: ErrorCode::Internal,
                                message: format!(
                                    "Could not evaluate anchored-region enzyme constraints: {e}"
                                ),
                            }
                        })?;
                        let mut enzymes_ok = true;
                        for enzyme in &required_enzymes {
                            if enzyme.get_sites(&frag_dna, None).is_empty() {
                                enzymes_ok = false;
                                break;
                            }
                        }
                        if !enzymes_ok {
                            continue;
                        }
                    }

                    candidates.push(Candidate {
                        start,
                        end,
                        len,
                        sequence: fragment_text,
                        score: len.abs_diff(target_length_bp),
                    });
                }

                if candidates.is_empty() {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message:
                            "No anchored-region candidate satisfied the configured constraints"
                                .to_string(),
                    });
                }

                candidates.sort_by(|a, b| {
                    a.score
                        .cmp(&b.score)
                        .then(a.start.cmp(&b.start))
                        .then(a.end.cmp(&b.end))
                });

                let require_unique = unique.unwrap_or(false);
                if require_unique && candidates.len() != 1 {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "ExtractAnchoredRegion unique=true requires exactly one candidate, found {}",
                            candidates.len()
                        ),
                    });
                }

                let limit = max_candidates
                    .unwrap_or(candidates.len())
                    .min(self.max_fragments_per_container());
                if limit == 0 {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: "ExtractAnchoredRegion max_candidates must be >= 1".to_string(),
                    });
                }
                if candidates.len() > limit {
                    result.warnings.push(format!(
                        "Anchored-region candidates truncated from {} to {} by max_candidates/max_fragments_per_container",
                        candidates.len(),
                        limit
                    ));
                    candidates.truncate(limit);
                }

                let prefix = output_prefix.unwrap_or_else(|| format!("{input}_anchored"));
                for (idx, cand) in candidates.into_iter().enumerate() {
                    let mut out =
                        DNAsequence::from_sequence(&cand.sequence).map_err(|e| EngineError {
                            code: ErrorCode::Internal,
                            message: format!("Could not create anchored-region sequence: {e}"),
                        })?;
                    out.set_circular(false);
                    Self::prepare_sequence(&mut out);
                    let seq_id = self.unique_seq_id(&format!("{}_{}", prefix, idx + 1));
                    self.state.sequences.insert(seq_id.clone(), out);
                    self.add_lineage_node(&seq_id, SequenceOrigin::Derived, Some(&result.op_id));
                    result.created_seq_ids.push(seq_id.clone());
                    result.messages.push(format!(
                        "Anchored-region candidate '{}' [{}..{}, {} bp] score={}",
                        seq_id, cand.start, cand.end, cand.len, cand.score
                    ));
                }
            }
            Operation::SelectCandidate {
                input,
                criterion,
                output_id,
            } => {
                parent_seq_ids.push(input.clone());
                let dna = self
                    .state
                    .sequences
                    .get(&input)
                    .ok_or_else(|| EngineError {
                        code: ErrorCode::NotFound,
                        message: format!("Sequence '{input}' not found"),
                    })?
                    .clone();

                let base = output_id.unwrap_or_else(|| format!("{input}_selected"));
                let seq_id = self.unique_seq_id(&base);
                self.state.sequences.insert(seq_id.clone(), dna);
                self.add_lineage_node(
                    &seq_id,
                    SequenceOrigin::InSilicoSelection,
                    Some(&result.op_id),
                );
                result.created_seq_ids.push(seq_id.clone());
                result.warnings.push(
                    "Selection operation is in-silico and may not directly correspond to a unique wet-lab product"
                        .to_string(),
                );
                result.messages.push(format!(
                    "Selected candidate '{}' from '{}' using criterion '{}'",
                    seq_id, input, criterion
                ));
            }
            Operation::FilterByMolecularWeight {
                inputs,
                min_bp,
                max_bp,
                error,
                unique,
                output_prefix,
            } => {
                if inputs.is_empty() {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: "FilterByMolecularWeight requires at least one input sequence"
                            .to_string(),
                    });
                }
                if min_bp > max_bp {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!("min_bp ({min_bp}) must be <= max_bp ({max_bp})"),
                    });
                }
                if !(0.0..=1.0).contains(&error) {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: "error must be between 0.0 and 1.0".to_string(),
                    });
                }

                let min_allowed = ((min_bp as f64) * (1.0 - error)).floor() as usize;
                let max_allowed = ((max_bp as f64) * (1.0 + error)).ceil() as usize;

                let mut matches: Vec<(SeqId, DNAsequence)> = vec![];
                for input in &inputs {
                    parent_seq_ids.push(input.clone());
                    let dna = self
                        .state
                        .sequences
                        .get(input)
                        .ok_or_else(|| EngineError {
                            code: ErrorCode::NotFound,
                            message: format!("Sequence '{input}' not found"),
                        })?
                        .clone();
                    let bp = dna.len();
                    if bp >= min_allowed && bp <= max_allowed {
                        matches.push((input.clone(), dna));
                    }
                }

                if matches.len() > self.max_fragments_per_container() {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "FilterByMolecularWeight produced {} candidates, exceeding max_fragments_per_container={}",
                            matches.len(),
                            self.max_fragments_per_container()
                        ),
                    });
                }

                if unique && matches.len() != 1 {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "unique=true requires exactly one match, found {}",
                            matches.len()
                        ),
                    });
                }

                let prefix = output_prefix.unwrap_or_else(|| "mw_filter".to_string());
                for (i, (_source_id, mut dna)) in matches.into_iter().enumerate() {
                    Self::prepare_sequence(&mut dna);
                    let candidate = format!("{}_{}", prefix, i + 1);
                    let seq_id = self.unique_seq_id(&candidate);
                    self.state.sequences.insert(seq_id.clone(), dna);
                    self.add_lineage_node(
                        &seq_id,
                        SequenceOrigin::InSilicoSelection,
                        Some(&result.op_id),
                    );
                    result.created_seq_ids.push(seq_id);
                }

                result.messages.push(format!(
                    "Molecular-weight filter kept {} sequence(s) in effective range {}-{} bp (requested {}-{} bp, error={:.3})",
                    result.created_seq_ids.len(),
                    min_allowed,
                    max_allowed,
                    min_bp,
                    max_bp,
                    error
                ));
            }
            Operation::FilterByDesignConstraints {
                inputs,
                gc_min,
                gc_max,
                max_homopolymer_run,
                reject_ambiguous_bases,
                avoid_u6_terminator_tttt,
                forbidden_motifs,
                unique,
                output_prefix,
            } => {
                if inputs.is_empty() {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: "FilterByDesignConstraints requires at least one input sequence"
                            .to_string(),
                    });
                }

                if let Some(min) = gc_min {
                    if !(0.0..=1.0).contains(&min) {
                        return Err(EngineError {
                            code: ErrorCode::InvalidInput,
                            message: format!("gc_min ({min}) must be between 0.0 and 1.0"),
                        });
                    }
                }
                if let Some(max) = gc_max {
                    if !(0.0..=1.0).contains(&max) {
                        return Err(EngineError {
                            code: ErrorCode::InvalidInput,
                            message: format!("gc_max ({max}) must be between 0.0 and 1.0"),
                        });
                    }
                }
                if let (Some(min), Some(max)) = (gc_min, gc_max) {
                    if min > max {
                        return Err(EngineError {
                            code: ErrorCode::InvalidInput,
                            message: format!("gc_min ({min}) must be <= gc_max ({max})"),
                        });
                    }
                }
                if let Some(max_run) = max_homopolymer_run {
                    if max_run == 0 {
                        return Err(EngineError {
                            code: ErrorCode::InvalidInput,
                            message: "max_homopolymer_run must be >= 1".to_string(),
                        });
                    }
                }

                let reject_ambiguous_bases = reject_ambiguous_bases.unwrap_or(true);
                let avoid_u6_terminator_tttt = avoid_u6_terminator_tttt.unwrap_or(true);
                let mut forbidden_motifs_normalized: Vec<String> = vec![];
                for motif in forbidden_motifs {
                    let normalized = Self::normalize_iupac_text(&motif)?;
                    if !normalized.is_empty() {
                        forbidden_motifs_normalized.push(normalized);
                    }
                }

                let mut matches: Vec<(SeqId, DNAsequence)> = vec![];
                let mut rejected = 0usize;
                let mut rejection_warnings_left = 32usize;
                let input_count = inputs.len();

                for input in &inputs {
                    parent_seq_ids.push(input.clone());
                    let dna = self
                        .state
                        .sequences
                        .get(input)
                        .ok_or_else(|| EngineError {
                            code: ErrorCode::NotFound,
                            message: format!("Sequence '{input}' not found"),
                        })?
                        .clone();
                    let sequence = Self::normalized_sequence_for_quality(&dna);
                    let mut reasons: Vec<String> = vec![];

                    if reject_ambiguous_bases && Self::has_ambiguous_bases(&sequence) {
                        reasons.push("contains_ambiguous_base".to_string());
                    }

                    if gc_min.is_some() || gc_max.is_some() {
                        match Self::sequence_gc_fraction(&sequence) {
                            Some(gc) => {
                                if let Some(min) = gc_min {
                                    if gc < min {
                                        reasons.push(format!("gc_too_low({gc:.3}<{min:.3})"));
                                    }
                                }
                                if let Some(max) = gc_max {
                                    if gc > max {
                                        reasons.push(format!("gc_too_high({gc:.3}>{max:.3})"));
                                    }
                                }
                            }
                            None => reasons.push("gc_not_computable".to_string()),
                        }
                    }

                    if let Some(max_run) = max_homopolymer_run {
                        let observed = Self::max_homopolymer_run(&sequence);
                        if observed > max_run {
                            reasons.push(format!("homopolymer_run_exceeded({observed}>{max_run})"));
                        }
                    }

                    if avoid_u6_terminator_tttt && Self::contains_u6_terminator_t4(&sequence) {
                        reasons.push("u6_terminator_t4".to_string());
                    }

                    for motif in &forbidden_motifs_normalized {
                        if Self::contains_motif_any_strand(&sequence, motif)? {
                            reasons.push(format!("forbidden_motif_present({motif})"));
                        }
                    }

                    if reasons.is_empty() {
                        matches.push((input.clone(), dna));
                    } else {
                        rejected += 1;
                        if rejection_warnings_left > 0 {
                            result.warnings.push(format!(
                                "Sequence '{}' rejected by design constraints: {}",
                                input,
                                reasons.join(", ")
                            ));
                            rejection_warnings_left -= 1;
                        }
                    }
                }

                if matches.len() > self.max_fragments_per_container() {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "FilterByDesignConstraints produced {} candidates, exceeding max_fragments_per_container={}",
                            matches.len(),
                            self.max_fragments_per_container()
                        ),
                    });
                }

                if unique && matches.len() != 1 {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "unique=true requires exactly one match, found {}",
                            matches.len()
                        ),
                    });
                }

                let prefix = output_prefix.unwrap_or_else(|| "design_filter".to_string());
                for (i, (_source_id, mut dna)) in matches.into_iter().enumerate() {
                    Self::prepare_sequence(&mut dna);
                    let candidate = format!("{}_{}", prefix, i + 1);
                    let seq_id = self.unique_seq_id(&candidate);
                    self.state.sequences.insert(seq_id.clone(), dna);
                    self.add_lineage_node(
                        &seq_id,
                        SequenceOrigin::InSilicoSelection,
                        Some(&result.op_id),
                    );
                    result.created_seq_ids.push(seq_id);
                }

                if rejected > 32 {
                    result.warnings.push(format!(
                        "{} additional sequence(s) were rejected (warning output truncated)",
                        rejected - 32
                    ));
                }

                result.messages.push(format!(
                    "Design-constraint filter kept {} of {} sequence(s)",
                    result.created_seq_ids.len(),
                    input_count
                ));
                result.messages.push(format!(
                    "Applied filters: gc_min={:?}, gc_max={:?}, max_homopolymer_run={:?}, reject_ambiguous_bases={}, avoid_u6_terminator_tttt={}, forbidden_motifs={}",
                    gc_min,
                    gc_max,
                    max_homopolymer_run,
                    reject_ambiguous_bases,
                    avoid_u6_terminator_tttt,
                    forbidden_motifs_normalized.len()
                ));
            }
            Operation::GenerateCandidateSet {
                set_name,
                seq_id,
                length_bp,
                step_bp,
                feature_kinds,
                feature_label_regex,
                max_distance_bp,
                feature_geometry_mode,
                feature_boundary_mode,
                feature_strand_relation,
                limit,
            } => {
                self.op_generate_candidate_set(
                    set_name,
                    seq_id,
                    length_bp,
                    step_bp,
                    feature_kinds,
                    feature_label_regex,
                    max_distance_bp,
                    feature_geometry_mode,
                    feature_boundary_mode,
                    feature_strand_relation,
                    limit,
                    &mut result,
                )?;
            }
            Operation::GenerateCandidateSetBetweenAnchors {
                set_name,
                seq_id,
                anchor_a,
                anchor_b,
                length_bp,
                step_bp,
                limit,
            } => {
                self.op_generate_candidate_set_between_anchors(
                    set_name,
                    seq_id,
                    anchor_a,
                    anchor_b,
                    length_bp,
                    step_bp,
                    limit,
                    &mut result,
                )?;
            }
            Operation::DeleteCandidateSet { set_name } => {
                self.op_delete_candidate_set(set_name, &mut result)?;
            }
            Operation::UpsertGuideSet {
                guide_set_id,
                guides,
            } => {
                self.op_upsert_guide_set(guide_set_id, guides, &mut result)?;
            }
            Operation::DeleteGuideSet { guide_set_id } => {
                self.op_delete_guide_set(guide_set_id, &mut result)?;
            }
            Operation::FilterGuidesPractical {
                guide_set_id,
                config,
                output_guide_set_id,
            } => {
                self.op_filter_guides_practical(
                    guide_set_id,
                    config,
                    output_guide_set_id,
                    &mut result,
                )?;
            }
            Operation::GenerateGuideOligos {
                guide_set_id,
                template_id,
                apply_5prime_g_extension,
                output_oligo_set_id,
                passed_only,
            } => {
                self.op_generate_guide_oligos(
                    guide_set_id,
                    template_id,
                    apply_5prime_g_extension,
                    output_oligo_set_id,
                    passed_only,
                    &mut result,
                )?;
            }
            Operation::ExportGuideOligos {
                guide_set_id,
                oligo_set_id,
                format,
                path,
                plate_format,
            } => {
                self.op_export_guide_oligos(
                    guide_set_id,
                    oligo_set_id,
                    format,
                    path,
                    plate_format,
                    &mut result,
                )?;
            }
            Operation::ExportGuideProtocolText {
                guide_set_id,
                oligo_set_id,
                path,
                include_qc_checklist,
            } => {
                self.op_export_guide_protocol_text(
                    guide_set_id,
                    oligo_set_id,
                    path,
                    include_qc_checklist,
                    &mut result,
                )?;
            }
            Operation::ScoreCandidateSetExpression {
                set_name,
                metric,
                expression,
            } => {
                self.op_score_candidate_set_expression(set_name, metric, expression, &mut result)?;
            }
            Operation::ScoreCandidateSetDistance {
                set_name,
                metric,
                feature_kinds,
                feature_label_regex,
                feature_geometry_mode,
                feature_boundary_mode,
                feature_strand_relation,
            } => {
                self.op_score_candidate_set_distance(
                    set_name,
                    metric,
                    feature_kinds,
                    feature_label_regex,
                    feature_geometry_mode,
                    feature_boundary_mode,
                    feature_strand_relation,
                    &mut result,
                )?;
            }
            Operation::FilterCandidateSet {
                input_set,
                output_set,
                metric,
                min,
                max,
                min_quantile,
                max_quantile,
            } => {
                self.op_filter_candidate_set(
                    input_set,
                    output_set,
                    metric,
                    min,
                    max,
                    min_quantile,
                    max_quantile,
                    &mut result,
                )?;
            }
            Operation::CandidateSetOp {
                op,
                left_set,
                right_set,
                output_set,
            } => {
                self.op_candidate_set_op(op, left_set, right_set, output_set, &mut result)?;
            }
            Operation::ScoreCandidateSetWeightedObjective {
                set_name,
                metric,
                objectives,
                normalize_metrics,
            } => {
                self.op_score_candidate_set_weighted_objective(
                    set_name,
                    metric,
                    objectives,
                    normalize_metrics,
                    &mut result,
                )?;
            }
            Operation::TopKCandidateSet {
                input_set,
                output_set,
                metric,
                k,
                direction,
                tie_break,
            } => {
                self.op_top_k_candidate_set(
                    input_set,
                    output_set,
                    metric,
                    k,
                    direction,
                    tie_break,
                    &mut result,
                )?;
            }
            Operation::ParetoFrontierCandidateSet {
                input_set,
                output_set,
                objectives,
                max_candidates,
                tie_break,
            } => {
                self.op_pareto_frontier_candidate_set(
                    input_set,
                    output_set,
                    objectives,
                    max_candidates,
                    tie_break,
                    &mut result,
                )?;
            }
            Operation::UpsertWorkflowMacroTemplate {
                name,
                description,
                details_url,
                parameters,
                input_ports,
                output_ports,
                script,
            } => {
                self.op_upsert_workflow_macro_template(
                    name,
                    description,
                    details_url,
                    parameters,
                    input_ports,
                    output_ports,
                    script,
                    &mut result,
                )?;
            }
            Operation::DeleteWorkflowMacroTemplate { name } => {
                self.op_delete_workflow_macro_template(name, &mut result)?;
            }
            Operation::UpsertCandidateMacroTemplate {
                name,
                description,
                details_url,
                parameters,
                script,
            } => {
                self.op_upsert_candidate_macro_template(
                    name,
                    description,
                    details_url,
                    parameters,
                    script,
                    &mut result,
                )?;
            }
            Operation::DeleteCandidateMacroTemplate { name } => {
                self.op_delete_candidate_macro_template(name, &mut result)?;
            }
            Operation::Reverse { input, output_id } => {
                parent_seq_ids.push(input.clone());
                let dna = self
                    .state
                    .sequences
                    .get(&input)
                    .ok_or_else(|| EngineError {
                        code: ErrorCode::NotFound,
                        message: format!("Sequence '{input}' not found"),
                    })?
                    .clone();
                let mut text = dna.get_forward_string();
                text = text.chars().rev().collect();
                let mut out = DNAsequence::from_sequence(&text).map_err(|e| EngineError {
                    code: ErrorCode::Internal,
                    message: format!("Could not create reverse sequence: {e}"),
                })?;
                out.set_circular(dna.is_circular());
                Self::prepare_sequence(&mut out);

                let base = output_id.unwrap_or_else(|| format!("{input}_rev"));
                let seq_id = self.unique_seq_id(&base);
                self.state.sequences.insert(seq_id.clone(), out);
                self.add_lineage_node(&seq_id, SequenceOrigin::Derived, Some(&result.op_id));
                result.created_seq_ids.push(seq_id.clone());
                result.messages.push(format!(
                    "Created reverse sequence '{}' from '{}'",
                    seq_id, input
                ));
            }
            Operation::Complement { input, output_id } => {
                parent_seq_ids.push(input.clone());
                let dna = self
                    .state
                    .sequences
                    .get(&input)
                    .ok_or_else(|| EngineError {
                        code: ErrorCode::NotFound,
                        message: format!("Sequence '{input}' not found"),
                    })?
                    .clone();
                let text: String = dna
                    .get_forward_string()
                    .as_bytes()
                    .iter()
                    .map(|c| IupacCode::letter_complement(*c))
                    .map(char::from)
                    .collect();
                let mut out = DNAsequence::from_sequence(&text).map_err(|e| EngineError {
                    code: ErrorCode::Internal,
                    message: format!("Could not create complement sequence: {e}"),
                })?;
                out.set_circular(dna.is_circular());
                Self::prepare_sequence(&mut out);

                let base = output_id.unwrap_or_else(|| format!("{input}_comp"));
                let seq_id = self.unique_seq_id(&base);
                self.state.sequences.insert(seq_id.clone(), out);
                self.add_lineage_node(&seq_id, SequenceOrigin::Derived, Some(&result.op_id));
                result.created_seq_ids.push(seq_id.clone());
                result.messages.push(format!(
                    "Created complement sequence '{}' from '{}'",
                    seq_id, input
                ));
            }
            Operation::ReverseComplement { input, output_id } => {
                parent_seq_ids.push(input.clone());
                let dna = self
                    .state
                    .sequences
                    .get(&input)
                    .ok_or_else(|| EngineError {
                        code: ErrorCode::NotFound,
                        message: format!("Sequence '{input}' not found"),
                    })?
                    .clone();
                let text = Self::reverse_complement(&dna.get_forward_string());
                let mut out = DNAsequence::from_sequence(&text).map_err(|e| EngineError {
                    code: ErrorCode::Internal,
                    message: format!("Could not create reverse-complement sequence: {e}"),
                })?;
                out.set_circular(dna.is_circular());
                Self::prepare_sequence(&mut out);

                let base = output_id.unwrap_or_else(|| format!("{input}_revcomp"));
                let seq_id = self.unique_seq_id(&base);
                self.state.sequences.insert(seq_id.clone(), out);
                self.add_lineage_node(&seq_id, SequenceOrigin::Derived, Some(&result.op_id));
                result.created_seq_ids.push(seq_id.clone());
                result.messages.push(format!(
                    "Created reverse-complement sequence '{}' from '{}'",
                    seq_id, input
                ));
            }
            Operation::Branch { input, output_id } => {
                parent_seq_ids.push(input.clone());
                let dna = self
                    .state
                    .sequences
                    .get(&input)
                    .ok_or_else(|| EngineError {
                        code: ErrorCode::NotFound,
                        message: format!("Sequence '{input}' not found"),
                    })?
                    .clone();

                let base = output_id.unwrap_or_else(|| format!("{input}_branch"));
                let seq_id = self.unique_seq_id(&base);
                self.state.sequences.insert(seq_id.clone(), dna);
                self.add_lineage_node(&seq_id, SequenceOrigin::Branch, Some(&result.op_id));
                result.created_seq_ids.push(seq_id.clone());
                result
                    .messages
                    .push(format!("Branched '{}' into '{}'", input, seq_id));
            }
            Operation::SetDisplayVisibility { target, visible } => {
                let (name, slot): (&str, &mut bool) = match target {
                    DisplayTarget::SequencePanel => (
                        "sequence_panel",
                        &mut self.state.display.show_sequence_panel,
                    ),
                    DisplayTarget::MapPanel => {
                        ("map_panel", &mut self.state.display.show_map_panel)
                    }
                    DisplayTarget::Features => ("features", &mut self.state.display.show_features),
                    DisplayTarget::CdsFeatures => {
                        ("cds_features", &mut self.state.display.show_cds_features)
                    }
                    DisplayTarget::GeneFeatures => {
                        ("gene_features", &mut self.state.display.show_gene_features)
                    }
                    DisplayTarget::MrnaFeatures => {
                        ("mrna_features", &mut self.state.display.show_mrna_features)
                    }
                    DisplayTarget::Tfbs => ("tfbs", &mut self.state.display.show_tfbs),
                    DisplayTarget::RestrictionEnzymes => (
                        "restriction_enzymes",
                        &mut self.state.display.show_restriction_enzymes,
                    ),
                    DisplayTarget::GcContents => {
                        ("gc_contents", &mut self.state.display.show_gc_contents)
                    }
                    DisplayTarget::OpenReadingFrames => (
                        "open_reading_frames",
                        &mut self.state.display.show_open_reading_frames,
                    ),
                    DisplayTarget::MethylationSites => (
                        "methylation_sites",
                        &mut self.state.display.show_methylation_sites,
                    ),
                };
                *slot = visible;
                result
                    .messages
                    .push(format!("Set display target '{name}' to {visible}"));
            }
            Operation::SetLinearViewport { start_bp, span_bp } => {
                self.state.display.linear_view_start_bp = start_bp;
                self.state.display.linear_view_span_bp = span_bp;
                result.messages.push(format!(
                    "Set linear viewport start_bp={start_bp}, span_bp={span_bp}"
                ));
            }
            Operation::SetTopology { seq_id, circular } => {
                let _ = self.ensure_lineage_node(&seq_id);
                let dna = self
                    .state
                    .sequences
                    .get_mut(&seq_id)
                    .ok_or_else(|| EngineError {
                        code: ErrorCode::NotFound,
                        message: format!("Sequence '{seq_id}' not found"),
                    })?;
                dna.set_circular(circular);
                dna.update_computed_features();
                result.changed_seq_ids.push(seq_id.clone());
                result.messages.push(format!(
                    "Set topology of '{seq_id}' to {}",
                    if circular { "circular" } else { "linear" }
                ));
            }
            Operation::RecomputeFeatures { seq_id } => {
                let _ = self.ensure_lineage_node(&seq_id);
                let dna = self
                    .state
                    .sequences
                    .get_mut(&seq_id)
                    .ok_or_else(|| EngineError {
                        code: ErrorCode::NotFound,
                        message: format!("Sequence '{seq_id}' not found"),
                    })?;
                dna.update_computed_features();
                result.changed_seq_ids.push(seq_id.clone());
                result
                    .messages
                    .push(format!("Recomputed features for '{seq_id}'"));
            }
            Operation::AnnotateTfbs {
                seq_id,
                motifs,
                min_llr_bits,
                min_llr_quantile,
                per_tf_thresholds,
                clear_existing,
                max_hits,
            } => {
                const DEFAULT_MAX_TFBS_HITS: usize = 500;
                let motifs = if motifs.len() == 1
                    && matches!(motifs[0].trim().to_ascii_uppercase().as_str(), "ALL" | "*")
                {
                    tf_motifs::all_motif_ids()
                } else {
                    motifs
                };

                if motifs.is_empty() {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: "AnnotateTfbs requires at least one motif".to_string(),
                    });
                }
                let default_min_llr_bits = min_llr_bits.unwrap_or(f64::NEG_INFINITY);
                let default_min_llr_quantile = min_llr_quantile.unwrap_or(0.0);
                Self::validate_tf_thresholds(default_min_llr_quantile)?;

                let mut override_map: HashMap<String, (Option<f64>, Option<f64>)> = HashMap::new();
                for o in &per_tf_thresholds {
                    let key = o.tf.trim().to_ascii_uppercase();
                    if key.is_empty() {
                        continue;
                    }
                    if let Some(q) = o.min_llr_quantile {
                        Self::validate_tf_thresholds(q)?;
                    }
                    override_map.insert(key, (o.min_llr_bits, o.min_llr_quantile));
                }

                let _ = self.ensure_lineage_node(&seq_id);
                let dna = self
                    .state
                    .sequences
                    .get_mut(&seq_id)
                    .ok_or_else(|| EngineError {
                        code: ErrorCode::NotFound,
                        message: format!("Sequence '{seq_id}' not found"),
                    })?;
                let seq_text = dna.get_forward_string();
                let seq_bytes = seq_text.as_bytes();

                if clear_existing.unwrap_or(true) {
                    Self::remove_generated_tfbs_features(dna.features_mut());
                }

                let mut added = 0usize;
                let max_hits = match max_hits {
                    Some(0) => None,
                    Some(v) => Some(v),
                    None => Some(DEFAULT_MAX_TFBS_HITS),
                };
                let motif_count = motifs.len();
                let mut motifs_scanned = 0usize;
                let mut cap_reached = false;
                'motif_loop: for (motif_idx, token) in motifs.into_iter().enumerate() {
                    let token_key = token.trim().to_ascii_uppercase();
                    if token_key.is_empty() {
                        continue;
                    }
                    let (tf_id, tf_name, _consensus, matrix_counts) =
                        Self::resolve_tf_motif_for_scoring(&token)?;
                    let (llr_matrix, true_log_odds_matrix) =
                        Self::prepare_scoring_matrices(&matrix_counts);
                    if llr_matrix.is_empty() || llr_matrix.len() > seq_bytes.len() {
                        result.warnings.push(format!(
                            "TF '{}' skipped: motif length {} exceeds sequence length {}",
                            tf_id,
                            llr_matrix.len(),
                            seq_bytes.len()
                        ));
                        continue;
                    }

                    let mut eff_bits = default_min_llr_bits;
                    let mut eff_quantile = default_min_llr_quantile;
                    let id_key = tf_id.to_ascii_uppercase();
                    let name_key = tf_name
                        .as_ref()
                        .map(|n| n.trim().to_ascii_uppercase())
                        .unwrap_or_default();
                    for key in [token_key.as_str(), id_key.as_str(), name_key.as_str()] {
                        if key.is_empty() {
                            continue;
                        }
                        if let Some((b, q)) = override_map.get(key) {
                            if let Some(v) = b {
                                eff_bits = *v;
                            }
                            if let Some(v) = q {
                                eff_quantile = *v;
                            }
                            break;
                        }
                    }

                    motifs_scanned += 1;
                    let seq_id_for_progress = seq_id.clone();
                    let tf_id_for_progress = tf_id.clone();
                    let motif_index = motif_idx + 1;
                    let hits = Self::scan_tf_scores(
                        seq_bytes,
                        &llr_matrix,
                        &true_log_odds_matrix,
                        |scanned_steps, total_steps| {
                            let motif_fraction = if total_steps == 0 {
                                1.0
                            } else {
                                (scanned_steps as f64 / total_steps as f64).clamp(0.0, 1.0)
                            };
                            let total_fraction = if motif_count == 0 {
                                1.0
                            } else {
                                ((motif_index - 1) as f64 + motif_fraction) / motif_count as f64
                            }
                            .clamp(0.0, 1.0);
                            on_progress(OperationProgress::Tfbs(TfbsProgress {
                                seq_id: seq_id_for_progress.clone(),
                                motif_id: tf_id_for_progress.clone(),
                                motif_index,
                                motif_count,
                                scanned_steps,
                                total_steps,
                                motif_percent: motif_fraction * 100.0,
                                total_percent: total_fraction * 100.0,
                            }));
                        },
                    );
                    let mut kept = 0usize;
                    for (
                        start,
                        reverse,
                        llr_bits,
                        llr_quantile,
                        true_log_odds_bits,
                        true_log_odds_quantile,
                    ) in hits
                    {
                        if llr_bits < eff_bits || llr_quantile < eff_quantile {
                            continue;
                        }
                        let end = start + llr_matrix.len();
                        dna.features_mut().push(Self::build_tfbs_feature(
                            start,
                            end,
                            reverse,
                            llr_matrix.len(),
                            &tf_id,
                            tf_name.as_deref(),
                            llr_bits,
                            llr_quantile,
                            true_log_odds_bits,
                            true_log_odds_quantile,
                        ));
                        kept += 1;
                        added += 1;
                        if let Some(limit) = max_hits {
                            if added >= limit {
                                cap_reached = true;
                                break;
                            }
                        }
                    }
                    result.messages.push(format!(
                        "TF '{}' annotated {} hit(s){}",
                        tf_id,
                        kept,
                        Self::format_tf_threshold_summary(eff_bits, eff_quantile)
                    ));
                    on_progress(OperationProgress::Tfbs(TfbsProgress {
                        seq_id: seq_id.clone(),
                        motif_id: tf_id,
                        motif_index,
                        motif_count,
                        scanned_steps: 1,
                        total_steps: 1,
                        motif_percent: 100.0,
                        total_percent: (motif_index as f64 / motif_count.max(1) as f64) * 100.0,
                    }));
                    if cap_reached {
                        if let Some(limit) = max_hits {
                            result.warnings.push(format!(
                                "TFBS hit cap ({limit}) reached after scanning {motifs_scanned}/{motif_count} motif(s); skipping remaining motif scans"
                            ));
                        }
                        break 'motif_loop;
                    }
                }
                result.messages.push(format!(
                    "TFBS motif scan coverage: {motifs_scanned}/{motif_count} motif(s)"
                ));

                result.changed_seq_ids.push(seq_id.clone());
                result.messages.push(format!(
                    "Annotated {} TFBS feature(s) on '{}'",
                    added, seq_id
                ));
            }
            Operation::SetParameter { name, value } => match name.as_str() {
                "max_fragments_per_container" => {
                    let raw = value.as_u64().ok_or_else(|| EngineError {
                        code: ErrorCode::InvalidInput,
                        message:
                            "SetParameter max_fragments_per_container requires a positive integer"
                                .to_string(),
                    })?;
                    if raw == 0 {
                        return Err(EngineError {
                            code: ErrorCode::InvalidInput,
                            message: "max_fragments_per_container must be >= 1".to_string(),
                        });
                    }
                    self.state.parameters.max_fragments_per_container = raw as usize;
                    result.messages.push(format!(
                        "Set parameter '{}' to {}",
                        name, self.state.parameters.max_fragments_per_container
                    ));
                }
                "genome_anchor_prepared_fallback_policy"
                | "genome_anchor_fallback_mode"
                | "genome_anchor_prepared_mode" => {
                    let raw = value.as_str().ok_or_else(|| EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!("SetParameter {name} requires a string value"),
                    })?;
                    let normalized = raw.trim().to_ascii_lowercase().replace('-', "_");
                    let policy = match normalized.as_str() {
                        "off" | "disabled" | "none" => GenomeAnchorPreparedFallbackPolicy::Off,
                        "single_compatible" | "single" | "auto" => {
                            GenomeAnchorPreparedFallbackPolicy::SingleCompatible
                        }
                        "always_explicit" | "explicit" => {
                            GenomeAnchorPreparedFallbackPolicy::AlwaysExplicit
                        }
                        other => {
                            return Err(EngineError {
                                code: ErrorCode::InvalidInput,
                                message: format!(
                                    "Unsupported {name} '{other}' (expected off|single_compatible|always_explicit)"
                                ),
                            });
                        }
                    };
                    self.state.parameters.genome_anchor_prepared_fallback_policy = policy;
                    result.messages.push(format!(
                        "Set parameter 'genome_anchor_prepared_fallback_policy' to {}",
                        policy.as_str()
                    ));
                }
                "require_verified_genome_anchor_for_extension"
                | "strict_genome_anchor_verification"
                | "strict_anchor_verification" => {
                    let raw = value.as_bool().ok_or_else(|| EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!("SetParameter {name} requires a boolean"),
                    })?;
                    self.state
                        .parameters
                        .require_verified_genome_anchor_for_extension = raw;
                    result.messages.push(format!(
                        "Set parameter 'require_verified_genome_anchor_for_extension' to {}",
                        self.state
                            .parameters
                            .require_verified_genome_anchor_for_extension
                    ));
                }
                "primer_design_backend" | "primers_design_backend" => {
                    let raw = value.as_str().ok_or_else(|| EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!("SetParameter {name} requires a string value"),
                    })?;
                    let normalized = raw.trim().to_ascii_lowercase();
                    let backend = match normalized.as_str() {
                        "auto" => PrimerDesignBackend::Auto,
                        "internal" | "baseline" => PrimerDesignBackend::Internal,
                        "primer3" | "primer3_core" | "external_primer3" | "external-primer3" => {
                            PrimerDesignBackend::Primer3
                        }
                        other => {
                            return Err(EngineError {
                                code: ErrorCode::InvalidInput,
                                message: format!(
                                    "Unsupported primer_design_backend '{other}' (expected auto|internal|primer3)"
                                ),
                            });
                        }
                    };
                    self.state.parameters.primer_design_backend = backend;
                    result.messages.push(format!(
                        "Set parameter 'primer_design_backend' to {}",
                        backend.as_str()
                    ));
                }
                "primer3_executable" | "primer3_backend_executable" | "primer3_path" => {
                    if value.is_null() {
                        self.state.parameters.primer3_executable = "primer3_core".to_string();
                        result.messages.push(
                            "Reset parameter 'primer3_executable' to default 'primer3_core'"
                                .to_string(),
                        );
                    } else {
                        let raw = value.as_str().ok_or_else(|| EngineError {
                            code: ErrorCode::InvalidInput,
                            message: format!("SetParameter {name} requires a string or null value"),
                        })?;
                        let trimmed = raw.trim();
                        if trimmed.is_empty() {
                            self.state.parameters.primer3_executable = "primer3_core".to_string();
                            result.messages.push(
                                "Reset parameter 'primer3_executable' to default 'primer3_core'"
                                    .to_string(),
                            );
                        } else {
                            self.state.parameters.primer3_executable = trimmed.to_string();
                            result.messages.push(format!(
                                "Set parameter 'primer3_executable' to '{}'",
                                trimmed
                            ));
                        }
                    }
                }
                "feature_details_font_size" | "feature_detail_font_size" => {
                    let raw = value.as_f64().ok_or_else(|| EngineError {
                        code: ErrorCode::InvalidInput,
                        message: "SetParameter feature_details_font_size requires a number"
                            .to_string(),
                    })?;
                    if !raw.is_finite() {
                        return Err(EngineError {
                            code: ErrorCode::InvalidInput,
                            message: "feature_details_font_size must be a finite number"
                                .to_string(),
                        });
                    }
                    if !(8.0..=24.0).contains(&raw) {
                        return Err(EngineError {
                            code: ErrorCode::InvalidInput,
                            message: "feature_details_font_size must be between 8.0 and 24.0"
                                .to_string(),
                        });
                    }
                    self.state.display.feature_details_font_size = raw as f32;
                    result.messages.push(format!(
                        "Set parameter 'feature_details_font_size' to {:.2}",
                        self.state.display.feature_details_font_size
                    ));
                }
                "linear_external_feature_label_font_size"
                | "linear_feature_label_font_size"
                | "feature_label_font_size" => {
                    let raw = value.as_f64().ok_or_else(|| EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!("SetParameter {name} requires a number"),
                    })?;
                    if !raw.is_finite() {
                        return Err(EngineError {
                            code: ErrorCode::InvalidInput,
                            message:
                                "linear_external_feature_label_font_size must be a finite number"
                                    .to_string(),
                        });
                    }
                    if !(8.0..=24.0).contains(&raw) {
                        return Err(EngineError {
                            code: ErrorCode::InvalidInput,
                            message: "linear_external_feature_label_font_size must be between 8.0 and 24.0".to_string(),
                        });
                    }
                    self.state.display.linear_external_feature_label_font_size = raw as f32;
                    result.messages.push(format!(
                        "Set parameter 'linear_external_feature_label_font_size' to {:.2}",
                        self.state.display.linear_external_feature_label_font_size
                    ));
                }
                "linear_external_feature_label_background_opacity"
                | "linear_feature_label_background_opacity"
                | "feature_label_background_opacity" => {
                    let raw = value.as_f64().ok_or_else(|| EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!("SetParameter {name} requires a number"),
                    })?;
                    if !raw.is_finite() {
                        return Err(EngineError {
                            code: ErrorCode::InvalidInput,
                            message: "linear_external_feature_label_background_opacity must be a finite number".to_string(),
                        });
                    }
                    if !(0.0..=1.0).contains(&raw) {
                        return Err(EngineError {
                            code: ErrorCode::InvalidInput,
                            message: "linear_external_feature_label_background_opacity must be between 0.0 and 1.0".to_string(),
                        });
                    }
                    self.state
                        .display
                        .linear_external_feature_label_background_opacity = raw as f32;
                    result.messages.push(format!(
                        "Set parameter 'linear_external_feature_label_background_opacity' to {:.3}",
                        self.state
                            .display
                            .linear_external_feature_label_background_opacity
                    ));
                }
                "reverse_strand_visual_opacity"
                | "linear_reverse_strand_visual_opacity"
                | "linear_reverse_strand_letter_opacity"
                | "reverse_strand_letter_opacity" => {
                    let raw = value.as_f64().ok_or_else(|| EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!("SetParameter {name} requires a number"),
                    })?;
                    if !raw.is_finite() {
                        return Err(EngineError {
                            code: ErrorCode::InvalidInput,
                            message: "reverse_strand_visual_opacity must be a finite number"
                                .to_string(),
                        });
                    }
                    if !(0.2..=1.0).contains(&raw) {
                        return Err(EngineError {
                            code: ErrorCode::InvalidInput,
                            message: "reverse_strand_visual_opacity must be between 0.2 and 1.0"
                                .to_string(),
                        });
                    }
                    self.state.display.reverse_strand_visual_opacity = raw as f32;
                    result.messages.push(format!(
                        "Set parameter 'reverse_strand_visual_opacity' to {:.3}",
                        self.state.display.reverse_strand_visual_opacity
                    ));
                }
                "regulatory_feature_max_view_span_bp"
                | "regulatory_max_view_span_bp"
                | "regulatory_max_span_bp" => {
                    let raw = value.as_u64().ok_or_else(|| EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!("SetParameter {name} requires a non-negative integer"),
                    })?;
                    self.state.display.regulatory_feature_max_view_span_bp = raw as usize;
                    result.messages.push(format!(
                        "Set parameter 'regulatory_feature_max_view_span_bp' to {}",
                        self.state.display.regulatory_feature_max_view_span_bp
                    ));
                }
                "sequence_panel_max_text_length_bp"
                | "sequence_text_panel_max_length_bp"
                | "sequence_panel_max_length_bp" => {
                    let raw = value.as_u64().ok_or_else(|| EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!("SetParameter {name} requires a non-negative integer"),
                    })?;
                    self.state.display.sequence_panel_max_text_length_bp = raw as usize;
                    let mode = if self.state.display.sequence_panel_max_text_length_bp == 0 {
                        "unlimited".to_string()
                    } else {
                        format!(
                            "{} bp",
                            self.state.display.sequence_panel_max_text_length_bp
                        )
                    };
                    result.messages.push(format!(
                        "Set parameter 'sequence_panel_max_text_length_bp' to {mode}"
                    ));
                }
                "linear_sequence_base_text_max_view_span_bp"
                | "linear_base_text_max_view_span_bp"
                | "sequence_base_text_max_view_span_bp" => {
                    // See docs/architecture.md (single shared rendering contract):
                    // fixed bp thresholds are compatibility-only while adaptive
                    // viewport-density routing is authoritative.
                    let _raw = value.as_u64().ok_or_else(|| EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!("SetParameter {name} requires a non-negative integer"),
                    })?;
                    result.messages.push(format!(
                        "Set parameter '{}' accepted as deprecated no-op under adaptive linear DNA letter routing",
                        name
                    ));
                }
                "gc_content_bin_size_bp" | "gc_bin_size_bp" => {
                    let raw = value.as_u64().ok_or_else(|| EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!("SetParameter {name} requires a positive integer"),
                    })?;
                    if raw == 0 {
                        return Err(EngineError {
                            code: ErrorCode::InvalidInput,
                            message: "gc_content_bin_size_bp must be >= 1".to_string(),
                        });
                    }
                    self.state.display.gc_content_bin_size_bp = raw as usize;
                    result.messages.push(format!(
                        "Set parameter 'gc_content_bin_size_bp' to {}",
                        self.state.display.gc_content_bin_size_bp
                    ));
                }
                "linear_sequence_helical_max_view_span_bp" | "linear_helical_max_view_span_bp" => {
                    let _raw = value.as_u64().ok_or_else(|| EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!("SetParameter {name} requires a non-negative integer"),
                    })?;
                    result.messages.push(format!(
                        "Set parameter '{}' accepted as deprecated no-op under adaptive linear DNA letter routing",
                        name
                    ));
                }
                "linear_sequence_condensed_max_view_span_bp"
                | "linear_condensed_max_view_span_bp" => {
                    let _raw = value.as_u64().ok_or_else(|| EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!("SetParameter {name} requires a non-negative integer"),
                    })?;
                    result.messages.push(format!(
                        "Set parameter '{}' accepted as deprecated no-op under adaptive linear DNA letter routing",
                        name
                    ));
                }
                "linear_sequence_letter_layout_mode" | "linear_helical_letter_layout_mode" => {
                    let raw = value.as_str().ok_or_else(|| EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!("SetParameter {name} requires a string value"),
                    })?;
                    let normalized = raw.trim().to_ascii_lowercase();
                    let layout_mode = match normalized.as_str() {
                        "auto_adaptive" | "auto-adaptive" | "auto" | "adaptive" => {
                            LinearSequenceLetterLayoutMode::AutoAdaptive
                        }
                        "standard_linear" | "standard-linear" | "standard" => {
                            LinearSequenceLetterLayoutMode::StandardLinear
                        }
                        "helical" | "continuous_helical" | "continuous-helical" | "continuous" => {
                            LinearSequenceLetterLayoutMode::ContinuousHelical
                        }
                        "condensed_10_row" | "condensed-10-row" | "condensed10row"
                        | "condensed" | "10_row" | "10-row" => {
                            LinearSequenceLetterLayoutMode::Condensed10Row
                        }
                        _ => {
                            return Err(EngineError {
                                code: ErrorCode::InvalidInput,
                                message: format!(
                                    "Unsupported linear sequence letter layout mode '{raw}' (expected auto|adaptive|standard|helical|condensed_10_row; legacy aliases accepted)"
                                ),
                            });
                        }
                    };
                    self.state.display.linear_sequence_letter_layout_mode = layout_mode;
                    result.messages.push(format!(
                        "Set parameter 'linear_sequence_letter_layout_mode' to {:?}",
                        self.state.display.linear_sequence_letter_layout_mode
                    ));
                }
                "linear_sequence_helical_phase_offset_bp" | "linear_helical_phase_offset_bp" => {
                    let raw = value.as_u64().ok_or_else(|| EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!("SetParameter {name} requires a non-negative integer"),
                    })?;
                    if raw > 9 {
                        return Err(EngineError {
                            code: ErrorCode::InvalidInput,
                            message:
                                "linear_sequence_helical_phase_offset_bp must be between 0 and 9"
                                    .to_string(),
                        });
                    }
                    self.state.display.linear_sequence_helical_phase_offset_bp = raw as usize;
                    result.messages.push(format!(
                        "Set parameter 'linear_sequence_helical_phase_offset_bp' to {}",
                        self.state.display.linear_sequence_helical_phase_offset_bp
                    ));
                }
                "blast_options_override" => {
                    if value.is_null() {
                        self.state
                            .metadata
                            .remove(BLAST_OPTIONS_OVERRIDE_METADATA_KEY);
                        result
                            .messages
                            .push("Cleared parameter 'blast_options_override'".to_string());
                    } else {
                        let parsed =
                            Self::parse_blast_run_options_from_value(&value, "SetParameter")?;
                        self.state.metadata.insert(
                            BLAST_OPTIONS_OVERRIDE_METADATA_KEY.to_string(),
                            serde_json::to_value(parsed).map_err(|e| EngineError {
                                code: ErrorCode::Internal,
                                message: format!(
                                    "Could not serialize blast_options_override metadata: {e}"
                                ),
                            })?,
                        );
                        result
                            .messages
                            .push("Set parameter 'blast_options_override'".to_string());
                    }
                }
                "blast_options_defaults_path" => {
                    if value.is_null() {
                        self.state
                            .metadata
                            .remove(BLAST_OPTIONS_DEFAULTS_PATH_METADATA_KEY);
                        result
                            .messages
                            .push("Cleared parameter 'blast_options_defaults_path'".to_string());
                    } else {
                        let raw = value.as_str().ok_or_else(|| EngineError {
                            code: ErrorCode::InvalidInput,
                            message:
                                "SetParameter blast_options_defaults_path requires a string or null"
                                    .to_string(),
                        })?;
                        let trimmed = raw.trim();
                        if trimmed.is_empty() {
                            self.state
                                .metadata
                                .remove(BLAST_OPTIONS_DEFAULTS_PATH_METADATA_KEY);
                            result.messages.push(
                                "Cleared parameter 'blast_options_defaults_path'".to_string(),
                            );
                        } else {
                            self.state.metadata.insert(
                                BLAST_OPTIONS_DEFAULTS_PATH_METADATA_KEY.to_string(),
                                serde_json::Value::String(trimmed.to_string()),
                            );
                            result.messages.push(format!(
                                "Set parameter 'blast_options_defaults_path' to '{}'",
                                trimmed
                            ));
                        }
                    }
                }
                "vcf_display_show_snp"
                | "vcf_display_show_ins"
                | "vcf_display_show_del"
                | "vcf_display_show_sv"
                | "vcf_display_show_other"
                | "vcf_display_pass_only"
                | "vcf_display_use_min_qual"
                | "vcf_display_use_max_qual"
                | "show_tfbs"
                | "tfbs_display_use_llr_bits"
                | "tfbs_display_use_llr_quantile"
                | "tfbs_display_use_true_log_odds_bits"
                | "tfbs_display_use_true_log_odds_quantile"
                | "auto_hide_sequence_panel_when_linear_bases_visible"
                | "linear_show_double_strand_bases"
                | "linear_show_reverse_strand_bases"
                | "show_linear_double_strand_bases"
                | "show_linear_reverse_strand_bases"
                | "linear_helical_parallel_strands"
                | "linear_sequence_helical_parallel_strands"
                | "linear_hide_backbone_when_sequence_bases_visible"
                | "hide_linear_backbone_when_bases_visible"
                | "linear_sequence_helical_letters_enabled"
                | "linear_reverse_strand_use_upside_down_letters"
                | "linear_reverse_strand_upside_down_letters" => {
                    let raw = value.as_bool().ok_or_else(|| EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!("SetParameter {} requires a boolean", name),
                    })?;
                    match name.as_str() {
                        "vcf_display_show_snp" => self.state.display.vcf_display_show_snp = raw,
                        "vcf_display_show_ins" => self.state.display.vcf_display_show_ins = raw,
                        "vcf_display_show_del" => self.state.display.vcf_display_show_del = raw,
                        "vcf_display_show_sv" => self.state.display.vcf_display_show_sv = raw,
                        "vcf_display_show_other" => self.state.display.vcf_display_show_other = raw,
                        "vcf_display_pass_only" => self.state.display.vcf_display_pass_only = raw,
                        "vcf_display_use_min_qual" => {
                            self.state.display.vcf_display_use_min_qual = raw
                        }
                        "vcf_display_use_max_qual" => {
                            self.state.display.vcf_display_use_max_qual = raw
                        }
                        "show_tfbs" => self.state.display.show_tfbs = raw,
                        "tfbs_display_use_llr_bits" => {
                            self.state.display.tfbs_display_use_llr_bits = raw
                        }
                        "tfbs_display_use_llr_quantile" => {
                            self.state.display.tfbs_display_use_llr_quantile = raw
                        }
                        "tfbs_display_use_true_log_odds_bits" => {
                            self.state.display.tfbs_display_use_true_log_odds_bits = raw
                        }
                        "tfbs_display_use_true_log_odds_quantile" => {
                            self.state.display.tfbs_display_use_true_log_odds_quantile = raw
                        }
                        "auto_hide_sequence_panel_when_linear_bases_visible" => {
                            self.state
                                .display
                                .auto_hide_sequence_panel_when_linear_bases_visible = raw
                        }
                        "linear_show_double_strand_bases"
                        | "linear_show_reverse_strand_bases"
                        | "show_linear_double_strand_bases"
                        | "show_linear_reverse_strand_bases" => {
                            self.state.display.linear_show_double_strand_bases = raw
                        }
                        "linear_helical_parallel_strands"
                        | "linear_sequence_helical_parallel_strands" => {
                            self.state.display.linear_helical_parallel_strands = raw
                        }
                        "linear_hide_backbone_when_sequence_bases_visible"
                        | "hide_linear_backbone_when_bases_visible" => {
                            self.state
                                .display
                                .linear_hide_backbone_when_sequence_bases_visible = raw
                        }
                        "linear_sequence_helical_letters_enabled" => {
                            self.state.display.linear_sequence_helical_letters_enabled = raw
                        }
                        "linear_reverse_strand_use_upside_down_letters"
                        | "linear_reverse_strand_upside_down_letters" => {
                            self.state
                                .display
                                .linear_reverse_strand_use_upside_down_letters = raw
                        }
                        _ => unreachable!(),
                    }
                    result
                        .messages
                        .push(format!("Set parameter '{}' to {}", name, raw));
                }
                "tfbs_display_min_llr_bits" | "tfbs_display_min_true_log_odds_bits" => {
                    let raw = value.as_f64().ok_or_else(|| EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!("SetParameter {} requires a number", name),
                    })?;
                    if !raw.is_finite() {
                        return Err(EngineError {
                            code: ErrorCode::InvalidInput,
                            message: format!("{} must be a finite number", name),
                        });
                    }
                    if name == "tfbs_display_min_llr_bits" {
                        self.state.display.tfbs_display_min_llr_bits = raw;
                    } else {
                        self.state.display.tfbs_display_min_true_log_odds_bits = raw;
                    }
                    result
                        .messages
                        .push(format!("Set parameter '{}' to {:.6}", name, raw));
                }
                "tfbs_display_min_llr_quantile" | "tfbs_display_min_true_log_odds_quantile" => {
                    let raw = value.as_f64().ok_or_else(|| EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!("SetParameter {} requires a number", name),
                    })?;
                    if !raw.is_finite() {
                        return Err(EngineError {
                            code: ErrorCode::InvalidInput,
                            message: format!("{} must be a finite number", name),
                        });
                    }
                    Self::validate_tf_thresholds(raw)?;
                    if name == "tfbs_display_min_llr_quantile" {
                        self.state.display.tfbs_display_min_llr_quantile = raw;
                    } else {
                        self.state.display.tfbs_display_min_true_log_odds_quantile = raw;
                    }
                    result
                        .messages
                        .push(format!("Set parameter '{}' to {:.6}", name, raw));
                }
                "vcf_display_min_qual" | "vcf_display_max_qual" => {
                    let raw = value.as_f64().ok_or_else(|| EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!("SetParameter {} requires a number", name),
                    })?;
                    if !raw.is_finite() {
                        return Err(EngineError {
                            code: ErrorCode::InvalidInput,
                            message: format!("{} must be a finite number", name),
                        });
                    }
                    if name == "vcf_display_min_qual" {
                        self.state.display.vcf_display_min_qual = raw;
                    } else {
                        self.state.display.vcf_display_max_qual = raw;
                    }
                    result
                        .messages
                        .push(format!("Set parameter '{}' to {:.6}", name, raw));
                }
                "vcf_display_required_info_keys"
                | "vcf_display_required_info_keys_csv"
                | "vcf_display_required_info" => {
                    let mut keys = if let Some(array) = value.as_array() {
                        let mut values = Vec::with_capacity(array.len());
                        for entry in array {
                            let Some(raw) = entry.as_str() else {
                                return Err(EngineError {
                                        code: ErrorCode::InvalidInput,
                                        message: "vcf_display_required_info_keys array entries must be strings".to_string(),
                                    });
                            };
                            values.push(raw.to_string());
                        }
                        values
                    } else if let Some(raw) = value.as_str() {
                        raw.split(',')
                            .map(str::trim)
                            .filter(|v| !v.is_empty())
                            .map(|v| v.to_string())
                            .collect::<Vec<_>>()
                    } else {
                        return Err(EngineError {
                                code: ErrorCode::InvalidInput,
                                message: "SetParameter vcf_display_required_info_keys requires a string (CSV) or string array".to_string(),
                            });
                    };
                    for key in &mut keys {
                        *key = key.trim().to_ascii_uppercase();
                    }
                    keys.retain(|key| !key.is_empty());
                    keys.sort();
                    keys.dedup();
                    self.state.display.vcf_display_required_info_keys = keys.clone();
                    result.messages.push(format!(
                        "Set parameter 'vcf_display_required_info_keys' to [{}]",
                        keys.join(",")
                    ));
                }
                _ => {
                    return Err(EngineError {
                        code: ErrorCode::Unsupported,
                        message: format!("Unknown parameter '{}'", name),
                    });
                }
            },
        }

        self.add_lineage_edges(
            &parent_seq_ids,
            &result.created_seq_ids,
            &result.op_id,
            run_id,
        );
        self.add_container_from_result(&op_for_containers, &result);

        Ok(result)
    }
}
