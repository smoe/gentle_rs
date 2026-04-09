//! Core operation dispatch handler for `GentleEngine::apply_internal`.
//!
//! This file is the mutation hub that maps one `Operation` variant onto the
//! appropriate helper routines. Keep adapter entry points thin and route new
//! engine behavior through this shared dispatch path.
//!
//! Start here when tracing:
//! - which helper actually executes a given `Operation` variant
//! - state mutations, `OpResult` messages, and journal side effects
//! - cross-cutting operation glue that does not belong in adapter code

use super::*;
use crate::{
    AMINO_ACIDS,
    amino_acids::{STOP_CODON, UNKNOWN_CODON},
    dna_ladder::default_dna_ladders,
    gibson_planning::{GibsonAssemblyPlan, derive_gibson_execution_plan},
    uniprot::UniprotNucleotideXref,
};

impl GentleEngine {
    fn gibson_arrangement_insert_seq_ids(plan: &GibsonAssemblyPlan) -> Vec<String> {
        let fragments_by_id = plan
            .fragments
            .iter()
            .map(|fragment| {
                (
                    fragment.id.trim().to_string(),
                    fragment.seq_id.trim().to_string(),
                )
            })
            .collect::<HashMap<_, _>>();
        plan.assembly_order
            .iter()
            .filter(|member| member.kind.trim().eq_ignore_ascii_case("fragment"))
            .filter_map(|member| fragments_by_id.get(member.id.trim()).cloned())
            .filter(|seq_id| !seq_id.trim().is_empty())
            .collect()
    }

    fn choose_gibson_arrangement_ladders(&self, lane_seq_ids: &[String]) -> Vec<String> {
        let mut lengths = lane_seq_ids
            .iter()
            .filter_map(|seq_id| self.state.sequences.get(seq_id).map(|dna| dna.len()))
            .filter(|bp| *bp > 0)
            .collect::<Vec<_>>();
        if lengths.is_empty() {
            return vec![];
        }
        lengths.sort_unstable();
        let min_bp = *lengths.first().unwrap_or(&1);
        let max_bp = *lengths.last().unwrap_or(&min_bp);
        default_dna_ladders().choose_for_range(min_bp, max_bp, 2)
    }

    fn maybe_create_gibson_arrangement(
        &mut self,
        plan: &GibsonAssemblyPlan,
        product_seq_id: &str,
        op_id: &str,
    ) -> Result<Option<String>, EngineError> {
        let destination_seq_id = plan.destination.seq_id.trim();
        if destination_seq_id.is_empty() || product_seq_id.trim().is_empty() {
            return Ok(None);
        }

        let mut lane_seq_ids = vec![destination_seq_id.to_string()];
        let insert_seq_ids = Self::gibson_arrangement_insert_seq_ids(plan);
        lane_seq_ids.extend(insert_seq_ids.clone());
        lane_seq_ids.push(product_seq_id.trim().to_string());

        let mut lane_container_ids: Vec<String> = vec![];
        for seq_id in &lane_seq_ids {
            let container_id =
                self.ensure_exact_singleton_container_for_seq(seq_id, Some(op_id))?;
            lane_container_ids.push(container_id);
        }

        let ladders = self.choose_gibson_arrangement_ladders(&lane_seq_ids);
        let arrangement_name = Some(format!("Gibson lane set: {}", product_seq_id.trim()));
        let mut role_labels = vec!["vector".to_string()];
        role_labels.extend(
            insert_seq_ids
                .iter()
                .enumerate()
                .map(|(idx, _)| format!("insert_{}", idx + 1)),
        );
        role_labels.push("product".to_string());
        let arrangement_id = self.add_serial_arrangement(
            &lane_container_ids,
            None,
            arrangement_name,
            Some(ladders),
            Some(role_labels),
            None,
        )?;
        Ok(Some(arrangement_id))
    }

    fn dotplot_svg_xml_escape(raw: &str) -> String {
        raw.replace('&', "&amp;")
            .replace('<', "&lt;")
            .replace('>', "&gt;")
            .replace('"', "&quot;")
            .replace('\'', "&apos;")
    }

    fn derive_extract_region_default_base(input: &str, from: usize, to: usize) -> String {
        let stem =
            Self::strip_dbsnp_auto_region_id_to_rsid(input).unwrap_or_else(|| input.to_string());
        let start_1based = from.saturating_add(1);
        let end_1based = to;
        format!("{stem}_{start_1based}_{end_1based}_region")
    }

    fn strip_dbsnp_auto_region_id_to_rsid(input: &str) -> Option<String> {
        let mut parts = input.split('_').collect::<Vec<_>>();
        if parts.len() < 4 {
            return None;
        }
        let rs_id = parts.first().copied()?.trim();
        if !(rs_id.starts_with("rs") && rs_id[2..].chars().all(|ch| ch.is_ascii_digit())) {
            return None;
        }
        let last = parts.pop()?.trim();
        let penultimate = parts.pop()?.trim();
        if !last.chars().all(|ch| ch.is_ascii_digit())
            || !penultimate.chars().all(|ch| ch.is_ascii_digit())
        {
            return None;
        }
        Some(rs_id.to_string())
    }

    fn extract_genome_region_into_state(
        &mut self,
        result: &mut OpResult,
        genome_id: &str,
        chromosome: &str,
        start_1based: usize,
        end_1based: usize,
        output_id: Option<SeqId>,
        annotation_scope: Option<GenomeAnnotationScope>,
        max_annotation_features: Option<usize>,
        include_genomic_annotation: Option<bool>,
        catalog_path: Option<String>,
        cache_dir: Option<String>,
        provenance_operation: &str,
    ) -> Result<SeqId, EngineError> {
        let catalog_path = catalog_path.unwrap_or_else(|| DEFAULT_GENOME_CATALOG_PATH.to_string());
        let catalog = GenomeCatalog::from_json_file(&catalog_path).map_err(|e| EngineError {
            code: ErrorCode::InvalidInput,
            message: format!("Could not open genome catalog '{catalog_path}': {e}"),
        })?;
        let sequence = catalog
            .get_sequence_region_with_cache(
                genome_id,
                chromosome,
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
            Self::normalize_id_token(genome_id),
            Self::normalize_id_token(chromosome),
            start_1based,
            end_1based
        );
        let base = output_id.unwrap_or(default_id);
        let seq_id = self.import_genome_slice_sequence(result, sequence, base)?;
        let requested_scope = Self::resolve_extract_region_annotation_scope(
            annotation_scope,
            include_genomic_annotation,
        );
        let effective_cap = match max_annotation_features {
            Some(0) => None,
            Some(value) => Some(value),
            None => matches!(requested_scope, GenomeAnnotationScope::Full)
                .then_some(DEFAULT_EXTRACT_REGION_ANNOTATION_FEATURE_CAP),
        };
        let mut genes: Vec<GenomeGeneRecord> = vec![];
        let mut transcripts: Vec<GenomeTranscriptRecord> = vec![];
        if !matches!(requested_scope, GenomeAnnotationScope::None) {
            match catalog.list_gene_regions(genome_id, cache_dir.as_deref()) {
                Ok(records) => {
                    genes = records
                        .into_iter()
                        .filter(|record| {
                            Self::genome_chromosome_matches(&record.chromosome, chromosome)
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
                genome_id,
                chromosome,
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
            dropped_feature_count: candidate_before_fallback.saturating_sub(attached_feature_count),
            genes_attached,
            transcripts_attached,
            exons_attached,
            cds_attached,
            fallback_applied: fallback_reason.is_some(),
            fallback_reason,
        });
        self.maybe_attach_known_helper_mcs_annotation(&seq_id, genome_id, result);
        let source_plan = catalog.source_plan(genome_id, cache_dir.as_deref()).ok();
        let inspection = catalog
            .inspect_prepared_genome(genome_id, cache_dir.as_deref())
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
            operation: provenance_operation.to_string(),
            genome_id: genome_id.to_string(),
            catalog_path: catalog_path.clone(),
            cache_dir: cache_dir.clone(),
            chromosome: Some(chromosome.to_string()),
            start_1based: Some(start_1based),
            end_1based: Some(end_1based),
            gene_query: None,
            occurrence: None,
            gene_extract_mode: None,
            promoter_upstream_bp: None,
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
        Ok(seq_id)
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
        #[derive(Clone)]
        struct OverlaySeriesSvgRenderData {
            label: String,
            color_rgb: [u8; 3],
            visible_cells: Vec<((i32, i32), f32)>,
        }

        const RENDER_MAX_POINTS: usize = 120_000;
        const CONNECT_DIAGONALS_MAX_CELLS: usize = 80_000;
        let overlay_mode = view.series_count > 1 || view.query_series.len() > 1;
        let flex_track = if overlay_mode { None } else { flex_track };
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

        if overlay_mode {
            let total_points: usize = view
                .query_series
                .iter()
                .map(|series| series.point_count)
                .sum();
            let average_query_span = view
                .query_series
                .iter()
                .map(|series| {
                    series
                        .span_end_0based
                        .saturating_sub(series.span_start_0based)
                        .max(1)
                })
                .sum::<usize>()
                .checked_div(view.query_series.len().max(1))
                .unwrap_or(1)
                .max(1);
            let reference_annotation = view
                .reference_annotation
                .as_ref()
                .filter(|track| !track.intervals.is_empty());

            let outer_margin = 18.0_f32;
            let left_margin = if reference_annotation.is_some() {
                78.0_f32
            } else {
                56.0_f32
            };
            let right_margin = 16.0_f32;
            let top_margin = 48.0_f32;
            let bottom_margin = 24.0_f32;
            let dotplot_width = 1600.0_f32;
            let dotplot_height = (dotplot_width
                * (reference_span as f32 / average_query_span as f32))
                .clamp(420.0, 1280.0);
            let canvas_width =
                outer_margin + left_margin + dotplot_width + right_margin + outer_margin;
            let canvas_height =
                outer_margin + top_margin + dotplot_height + bottom_margin + outer_margin;
            let dotplot_left = outer_margin + left_margin;
            let dotplot_top = outer_margin + top_margin;
            let dotplot_right = dotplot_left + dotplot_width;
            let dotplot_bottom = dotplot_top + dotplot_height;
            let cols = dotplot_width.max(2.0).round() as i32;
            let rows = dotplot_height.max(2.0).round() as i32;
            let reference_span_max = reference_span.saturating_sub(1).max(1);
            let series_sample_cap = (RENDER_MAX_POINTS / view.query_series.len().max(1)).max(1);
            let mut rendered_series: Vec<OverlaySeriesSvgRenderData> = vec![];
            let mut total_visible_cells = 0usize;

            for series in &view.query_series {
                let query_series_span = series
                    .span_end_0based
                    .saturating_sub(series.span_start_0based)
                    .max(1);
                let query_series_span_max = query_series_span.saturating_sub(1).max(1);
                let sample_stride = (series.points.len() / series_sample_cap).max(1);
                let mut cells: HashMap<(i32, i32), usize> = HashMap::new();
                for point in series.points.iter().step_by(sample_stride) {
                    let x_local = point
                        .x_0based
                        .saturating_sub(series.span_start_0based)
                        .min(query_series_span_max);
                    let y_local = point
                        .y_0based
                        .saturating_sub(view.reference_span_start_0based)
                        .min(reference_span_max);
                    let x_frac = (x_local as f32 / query_series_span_max as f32).clamp(0.0, 1.0);
                    let y_frac = (y_local as f32 / reference_span_max as f32).clamp(0.0, 1.0);
                    let x_cell = ((x_frac * (cols - 1) as f32).round() as i32).clamp(0, cols - 1);
                    let y_cell = ((y_frac * (rows - 1) as f32).round() as i32).clamp(0, rows - 1);
                    let entry = cells.entry((x_cell, y_cell)).or_insert(0);
                    *entry = entry.saturating_add(1);
                }
                let max_cell_count = cells.values().copied().max().unwrap_or(1).max(1) as f32;
                let mut visible_cells: Vec<((i32, i32), f32)> = vec![];
                for ((x_cell, y_cell), count) in &cells {
                    let density_raw = (*count as f32 / max_cell_count).clamp(0.0, 1.0);
                    if density_raw < density_threshold {
                        continue;
                    }
                    let normalized = if density_threshold <= 0.0 {
                        density_raw
                    } else {
                        ((density_raw - density_threshold) / (1.0 - density_threshold))
                            .clamp(0.0, 1.0)
                    };
                    let density = (normalized * intensity_gain).clamp(0.0, 1.0).sqrt();
                    let opacity = ((72.0 + 178.0 * density) / 255.0).clamp(0.0, 1.0);
                    visible_cells.push(((*x_cell, *y_cell), opacity));
                    total_visible_cells = total_visible_cells.saturating_add(1);
                }
                visible_cells.sort_by_key(|((x, y), _)| (*y, *x));
                rendered_series.push(OverlaySeriesSvgRenderData {
                    label: series.label.clone(),
                    color_rgb: series.color_rgb,
                    visible_cells,
                });
            }

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
            for series in &rendered_series {
                for ((x_cell, y_cell), opacity) in &series.visible_cells {
                    let x0 = dotplot_left + (*x_cell as f32 / cols as f32) * dotplot_width;
                    let y0 = dotplot_top + (*y_cell as f32 / rows as f32) * dotplot_height;
                    let x1 = dotplot_left + ((*x_cell + 1) as f32 / cols as f32) * dotplot_width;
                    let y1 = dotplot_top + ((*y_cell + 1) as f32 / rows as f32) * dotplot_height;
                    svg.push_str(&format!(
                        "<rect x=\"{:.2}\" y=\"{:.2}\" width=\"{:.2}\" height=\"{:.2}\" fill=\"rgb({},{},{})\" fill-opacity=\"{:.4}\"/>",
                        x0,
                        y0,
                        (x1 - x0).max(0.8),
                        (y1 - y0).max(0.8),
                        series.color_rgb[0],
                        series.color_rgb[1],
                        series.color_rgb[2],
                        opacity
                    ));
                }
                if !series.visible_cells.is_empty()
                    && series.visible_cells.len()
                        <= CONNECT_DIAGONALS_MAX_CELLS
                            .checked_div(rendered_series.len().max(1))
                            .unwrap_or(CONNECT_DIAGONALS_MAX_CELLS)
                            .max(1)
                {
                    let visible_lookup = series
                        .visible_cells
                        .iter()
                        .map(|(cell, opacity)| (*cell, *opacity))
                        .collect::<HashMap<_, _>>();
                    let stroke_width = (cell_w.min(cell_h) * 0.45).clamp(0.8, 2.0);
                    for ((x_cell, y_cell), opacity) in &series.visible_cells {
                        for dy in [-1, 0, 1] {
                            let neighbor = (*x_cell + 1, *y_cell + dy);
                            if visible_lookup.contains_key(&neighbor) {
                                let x0 = dotplot_left + (*x_cell as f32 + 0.5) * cell_w;
                                let y0 = dotplot_top + (*y_cell as f32 + 0.5) * cell_h;
                                let x1 = dotplot_left + (neighbor.0 as f32 + 0.5) * cell_w;
                                let y1 = dotplot_top + (neighbor.1 as f32 + 0.5) * cell_h;
                                svg.push_str(&format!(
                                    "<line x1=\"{:.2}\" y1=\"{:.2}\" x2=\"{:.2}\" y2=\"{:.2}\" stroke=\"rgb({},{},{})\" stroke-opacity=\"{:.4}\" stroke-width=\"{:.2}\"/>",
                                    x0,
                                    y0,
                                    x1,
                                    y1,
                                    series.color_rgb[0],
                                    series.color_rgb[1],
                                    series.color_rgb[2],
                                    (opacity * 0.5).clamp(0.0, 1.0),
                                    stroke_width
                                ));
                            }
                        }
                    }
                }
            }

            if total_visible_cells == 0 {
                svg.push_str(&format!(
                    "<text x=\"{:.2}\" y=\"{:.2}\" text-anchor=\"middle\" dominant-baseline=\"middle\" font-family=\"monospace\" font-size=\"14\" fill=\"#64748b\">{}</text>",
                    (dotplot_left + dotplot_right) * 0.5,
                    (dotplot_top + dotplot_bottom) * 0.5,
                    Self::dotplot_svg_xml_escape(&format!(
                        "No visible cells (threshold={:.2}, total_points={}, series={})",
                        density_threshold,
                        total_points,
                        rendered_series.len()
                    ))
                ));
            }

            if let Some(track) = reference_annotation {
                let annotation_left = dotplot_left - 16.0;
                let annotation_width = 10.0;
                svg.push_str(&format!(
                    "<rect x=\"{:.1}\" y=\"{:.1}\" width=\"{:.1}\" height=\"{:.1}\" fill=\"#f8fafc\" stroke=\"#cbd5e1\" stroke-width=\"1\" rx=\"2\" ry=\"2\"/>",
                    annotation_left,
                    dotplot_top,
                    annotation_width,
                    dotplot_height
                ));
                let reference_span_f32 = reference_span.max(1) as f32;
                for interval in &track.intervals {
                    let local_start = interval
                        .start_0based
                        .saturating_sub(view.reference_span_start_0based)
                        .min(reference_span);
                    let local_end = interval
                        .end_0based_exclusive
                        .saturating_sub(view.reference_span_start_0based)
                        .min(reference_span);
                    let y0 =
                        dotplot_top + (local_start as f32 / reference_span_f32) * dotplot_height;
                    let y1 = dotplot_top + (local_end as f32 / reference_span_f32) * dotplot_height;
                    svg.push_str(&format!(
                        "<rect x=\"{:.2}\" y=\"{:.2}\" width=\"{:.2}\" height=\"{:.2}\" fill=\"#22c55e\" rx=\"1\" ry=\"1\"/>",
                        annotation_left + 1.0,
                        y0,
                        annotation_width - 2.0,
                        (y1 - y0).max(1.0)
                    ));
                }
                svg.push_str(&format!(
                    "<text x=\"{:.1}\" y=\"{:.1}\" text-anchor=\"middle\" font-family=\"monospace\" font-size=\"9\" fill=\"#334155\">{}</text>",
                    annotation_left + annotation_width * 0.5,
                    outer_margin + 14.0,
                    Self::dotplot_svg_xml_escape(track.label.as_str())
                ));
            }

            let header = format!(
                "Dotplot workspace export: {} | overlay owner={} reference={} [{}..{}] | series={} total_points={} | word={} step={} mismatches={} | threshold={:.2} gain={:.2}",
                view.dotplot_id,
                view.owner_seq_id,
                reference_seq_label,
                view.reference_span_start_0based.saturating_add(1),
                view.reference_span_end_0based,
                rendered_series.len(),
                total_points,
                view.word_size,
                view.step_bp,
                view.max_mismatches,
                density_threshold,
                intensity_gain
            );
            let sampling_line = format!(
                "{} | rendered_cells={} | merged_exons={}",
                Self::dotplot_sampling_overlap_summary(view.word_size.max(1), view.step_bp.max(1)),
                total_visible_cells,
                reference_annotation
                    .map(|track| track.interval_count)
                    .unwrap_or(0)
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

            let mut legend_x = dotplot_left;
            let legend_y = outer_margin + 42.0;
            for series in &rendered_series {
                svg.push_str(&format!(
                    "<rect x=\"{:.2}\" y=\"{:.2}\" width=\"10\" height=\"10\" fill=\"rgb({},{},{})\" rx=\"1\" ry=\"1\"/>",
                    legend_x,
                    legend_y - 8.0,
                    series.color_rgb[0],
                    series.color_rgb[1],
                    series.color_rgb[2]
                ));
                svg.push_str(&format!(
                    "<text x=\"{:.2}\" y=\"{:.2}\" font-family=\"monospace\" font-size=\"9.5\" fill=\"#334155\">{}</text>",
                    legend_x + 14.0,
                    legend_y,
                    Self::dotplot_svg_xml_escape(series.label.as_str())
                ));
                legend_x += 24.0 + series.label.len() as f32 * 6.1;
            }

            svg.push_str(&format!(
                "<text x=\"{:.1}\" y=\"{:.1}\" font-family=\"monospace\" font-size=\"10\" fill=\"#334155\">x: normalized isoform queries</text>",
                dotplot_left,
                outer_margin + 56.0
            ));
            svg.push_str(&format!(
                "<text x=\"{:.1}\" y=\"{:.1}\" text-anchor=\"end\" font-family=\"monospace\" font-size=\"10\" fill=\"#334155\">y: {}</text>",
                dotplot_right,
                outer_margin + 56.0,
                Self::dotplot_svg_xml_escape(reference_seq_label)
            ));
            svg.push_str(&format!(
                "<text x=\"{:.1}\" y=\"{:.1}\" font-family=\"monospace\" font-size=\"10\" fill=\"#475569\">0%</text>",
                dotplot_left,
                dotplot_bottom + 14.0
            ));
            svg.push_str(&format!(
                "<text x=\"{:.1}\" y=\"{:.1}\" text-anchor=\"end\" font-family=\"monospace\" font-size=\"10\" fill=\"#475569\">100%</text>",
                dotplot_right,
                dotplot_bottom + 14.0
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
            svg.push_str(&format!(
                "<text x=\"{:.1}\" y=\"{:.1}\" text-anchor=\"end\" font-family=\"monospace\" font-size=\"10\" fill=\"#64748b\">GENtle dotplot SVG export</text>",
                canvas_width - outer_margin,
                canvas_height - 6.0
            ));
            svg.push_str("</svg>");
            return svg;
        }

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
        if !view.owner_seq_id.eq_ignore_ascii_case(seq_id) {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Dotplot '{}' belongs to seq_id='{}', not '{}'",
                    normalized_dotplot_id, view.owner_seq_id, seq_id
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
            .qualifier_values(key)
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

    fn qualifier_usize_for_derivation(feature: &gb_io::seq::Feature, key: &str) -> Option<usize> {
        Self::qualifier_text_for_derivation(feature, key).and_then(|value| {
            value
                .trim()
                .parse::<usize>()
                .ok()
                .filter(|parsed| *parsed > 0)
        })
    }

    fn translation_table_display_name(table_id: usize) -> String {
        match table_id {
            11 => "Bacterial, Archaeal and Plant Plastid".to_string(),
            _ => AMINO_ACIDS
                .codon_tables
                .get(&table_id)
                .map(|table| table.organism.trim().to_string())
                .filter(|label| !label.is_empty())
                .unwrap_or_else(|| format!("translation table {table_id}")),
        }
    }

    fn infer_translation_speed_profile_hint(organism: Option<&str>) -> Option<String> {
        let normalized = organism?.trim().to_ascii_lowercase();
        if normalized.contains("homo sapiens") {
            return Some("human".to_string());
        }
        if normalized.contains("mus musculus") {
            return Some("mouse".to_string());
        }
        if normalized.contains("saccharomyces cerevisiae") {
            return Some("yeast".to_string());
        }
        if normalized.contains("escherichia coli") || normalized.contains("e. coli") {
            return Some("ecoli".to_string());
        }
        None
    }

    fn first_source_feature_for_derivation(
        source_features: &[gb_io::seq::Feature],
    ) -> Option<&gb_io::seq::Feature> {
        source_features
            .iter()
            .find(|feature| feature.kind.to_string().eq_ignore_ascii_case("source"))
    }

    fn feature_ranges_0based_for_derivation(feature: &gb_io::seq::Feature) -> Vec<(usize, usize)> {
        let mut ranges: Vec<(usize, usize)> = vec![];
        collect_location_ranges_usize(&feature.location, &mut ranges);
        if ranges.is_empty()
            && let Ok((from, to)) = feature.location.find_bounds()
            && from >= 0
            && to >= 0
            && (to as usize) > (from as usize)
        {
            ranges.push((from as usize, to as usize));
        }
        ranges.sort_unstable_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
        ranges.dedup();
        ranges.retain(|(start, end)| *end > *start);
        ranges
    }

    fn ranges_fit_within_exons(
        ranges_0based: &[(usize, usize)],
        exon_ranges_0based: &[(usize, usize)],
    ) -> bool {
        !ranges_0based.is_empty()
            && ranges_0based.iter().all(|(start, end)| {
                exon_ranges_0based
                    .iter()
                    .any(|(exon_start, exon_end)| *start >= *exon_start && *end <= *exon_end)
            })
    }

    fn merge_adjacent_ranges_0based(ranges_0based: &[(usize, usize)]) -> Vec<(usize, usize)> {
        let mut sorted = ranges_0based.to_vec();
        sorted.sort_unstable_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
        let mut merged: Vec<(usize, usize)> = vec![];
        for (start, end) in sorted {
            if end <= start {
                continue;
            }
            if let Some(last) = merged.last_mut()
                && start <= last.1
            {
                last.1 = last.1.max(end);
            } else {
                merged.push((start, end));
            }
        }
        merged
    }

    fn trim_leading_bases_from_ranges_0based(
        ranges_0based: &[(usize, usize)],
        mut trim_bp: usize,
    ) -> Vec<(usize, usize)> {
        let mut trimmed: Vec<(usize, usize)> = vec![];
        for (start, end) in ranges_0based {
            let span = end.saturating_sub(*start);
            if span == 0 {
                continue;
            }
            if trim_bp >= span {
                trim_bp -= span;
                continue;
            }
            trimmed.push((start.saturating_add(trim_bp), *end));
            trim_bp = 0;
        }
        trimmed
    }

    fn collect_matching_cds_features_for_derivation<'a>(
        source_features: &'a [gb_io::seq::Feature],
        source_feature: &gb_io::seq::Feature,
        exon_ranges_0based: &[(usize, usize)],
    ) -> Vec<&'a gb_io::seq::Feature> {
        let transcript_id = Self::qualifier_text_for_derivation(source_feature, "transcript_id");
        let transcript_gene = Self::first_nonempty_qualifier_for_derivation(
            source_feature,
            &["gene_id", "gene", "locus_tag"],
        );
        let transcript_reverse = feature_is_reverse(source_feature);
        let mut matches = source_features
            .iter()
            .filter(|feature| feature.kind.to_string().eq_ignore_ascii_case("CDS"))
            .filter(|feature| feature_is_reverse(feature) == transcript_reverse)
            .filter(|feature| {
                let feature_ranges = Self::feature_ranges_0based_for_derivation(feature);
                Self::ranges_fit_within_exons(&feature_ranges, exon_ranges_0based)
            })
            .filter(|feature| {
                if let Some(transcript_id) = transcript_id.as_deref() {
                    return feature
                        .qualifier_values("transcript_id")
                        .any(|value| value.trim().eq_ignore_ascii_case(transcript_id));
                }
                if let Some(transcript_gene) = transcript_gene.as_deref() {
                    return ["gene_id", "gene", "locus_tag"].iter().any(|key| {
                        feature
                            .qualifier_values(key)
                            .any(|value| value.trim().eq_ignore_ascii_case(transcript_gene))
                    });
                }
                false
            })
            .collect::<Vec<_>>();
        matches.sort_by(|left, right| {
            let mut left_ranges = Self::feature_ranges_0based_for_derivation(left);
            let mut right_ranges = Self::feature_ranges_0based_for_derivation(right);
            left_ranges.sort_unstable_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
            right_ranges.sort_unstable_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
            left_ranges.cmp(&right_ranges)
        });
        matches
    }

    fn map_source_ranges_to_transcript_local_ranges_0based(
        source_ranges_0based: &[(usize, usize)],
        exon_segments_forward: &[(usize, usize, usize, usize)],
        is_reverse: bool,
        total_len: usize,
    ) -> Vec<(usize, usize)> {
        let mut local_ranges: Vec<(usize, usize)> = vec![];
        for (source_start, source_end) in source_ranges_0based {
            for (exon_start, exon_end, local_start, local_end) in exon_segments_forward {
                let overlap_start = (*source_start).max(*exon_start);
                let overlap_end = (*source_end).min(*exon_end);
                if overlap_end <= overlap_start {
                    continue;
                }
                let offset_start = overlap_start.saturating_sub(*exon_start);
                let offset_end = overlap_end.saturating_sub(*exon_start);
                let mapped_start = local_start.saturating_add(offset_start);
                let mapped_end = local_start.saturating_add(offset_end);
                if mapped_end > mapped_start && mapped_end <= *local_end {
                    local_ranges.push((mapped_start, mapped_end));
                }
            }
        }
        let mut merged = Self::merge_adjacent_ranges_0based(&local_ranges);
        if is_reverse {
            merged = merged
                .into_iter()
                .map(|(start, end)| {
                    (
                        total_len.saturating_sub(end),
                        total_len.saturating_sub(start),
                    )
                })
                .collect();
            merged.sort_unstable_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
            merged = Self::merge_adjacent_ranges_0based(&merged);
        }
        merged
    }

    fn infer_organelle_for_derivation(
        source_feature: &gb_io::seq::Feature,
        cds_feature: Option<&gb_io::seq::Feature>,
        source_context_feature: Option<&gb_io::seq::Feature>,
    ) -> Option<String> {
        let explicit = cds_feature
            .and_then(|feature| Self::qualifier_text_for_derivation(feature, "organelle"))
            .or_else(|| Self::qualifier_text_for_derivation(source_feature, "organelle"))
            .or_else(|| {
                source_context_feature
                    .and_then(|feature| Self::qualifier_text_for_derivation(feature, "organelle"))
            });
        if explicit.is_some() {
            return explicit;
        }
        let chromosome_hint = cds_feature
            .and_then(|feature| Self::qualifier_text_for_derivation(feature, "chromosome"))
            .or_else(|| Self::qualifier_text_for_derivation(source_feature, "chromosome"))
            .or_else(|| {
                source_context_feature
                    .and_then(|feature| Self::qualifier_text_for_derivation(feature, "chromosome"))
            })?;
        (Self::normalize_chromosome_alias(&chromosome_hint) == "mt")
            .then_some("mitochondrion".to_string())
    }

    fn resolve_translation_table_for_derivation(
        source_feature: &gb_io::seq::Feature,
        cds_feature: Option<&gb_io::seq::Feature>,
        source_context_feature: Option<&gb_io::seq::Feature>,
        organelle: Option<&str>,
    ) -> (usize, TranscriptProteinTranslationTableSource, Vec<String>) {
        let mut warnings: Vec<String> = vec![];
        if let Some(value) = cds_feature
            .and_then(|feature| Self::qualifier_usize_for_derivation(feature, "transl_table"))
        {
            return (
                value,
                TranscriptProteinTranslationTableSource::ExplicitCdsQualifier,
                warnings,
            );
        }
        if let Some(value) = Self::qualifier_usize_for_derivation(source_feature, "transl_table") {
            return (
                value,
                TranscriptProteinTranslationTableSource::ExplicitTranscriptQualifier,
                warnings,
            );
        }
        if let Some(value) = source_context_feature
            .and_then(|feature| Self::qualifier_usize_for_derivation(feature, "transl_table"))
        {
            return (
                value,
                TranscriptProteinTranslationTableSource::ExplicitSourceQualifier,
                warnings,
            );
        }

        let organelle_normalized = organelle
            .map(|value| value.trim().to_ascii_lowercase())
            .unwrap_or_default();
        if organelle_normalized.contains("plastid")
            || organelle_normalized.contains("chloroplast")
            || organelle_normalized.contains("chromoplast")
            || organelle_normalized.contains("leucoplast")
            || organelle_normalized.contains("apicoplast")
        {
            return (
                11,
                TranscriptProteinTranslationTableSource::OrganellePlastidDefault,
                warnings,
            );
        }
        if organelle_normalized.contains("mitochond")
            || organelle_normalized.contains("kinetoplast")
        {
            warnings.push(
                "Mitochondrial context was detected but no explicit /transl_table qualifier was present; defaulted to translation table 1 because lineage-specific mitochondrial code inference is not implemented yet.".to_string(),
            );
            return (
                1,
                TranscriptProteinTranslationTableSource::AmbiguousMitochondrialDefault,
                warnings,
            );
        }
        (
            1,
            TranscriptProteinTranslationTableSource::StandardDefault,
            warnings,
        )
    }

    fn translate_transcript_cds_to_protein(
        derived_sequence: &str,
        cds_ranges_0based: &[(usize, usize)],
        translation_table: usize,
    ) -> (String, bool, Vec<String>) {
        let mut warnings: Vec<String> = vec![];
        let bytes = derived_sequence.as_bytes();
        let mut cds_bytes: Vec<u8> = vec![];
        for (start, end) in cds_ranges_0based {
            if *end > *start && *end <= bytes.len() {
                cds_bytes.extend_from_slice(&bytes[*start..*end]);
            }
        }
        if cds_bytes.is_empty() {
            return (String::new(), false, warnings);
        }
        if !cds_bytes.len().is_multiple_of(3) {
            warnings.push(format!(
                "CDS length {} bp is not divisible by 3; trailing partial codon was ignored.",
                cds_bytes.len()
            ));
        }
        let mut protein = String::new();
        for codon in cds_bytes.chunks(3) {
            if codon.len() < 3 {
                break;
            }
            let aa = AMINO_ACIDS.codon2aa([codon[0], codon[1], codon[2]], Some(translation_table));
            protein.push(match aa {
                STOP_CODON => '*',
                UNKNOWN_CODON => 'X',
                other => other,
            });
        }
        let terminal_stop_trimmed = protein.ends_with('*');
        if terminal_stop_trimmed {
            protein.pop();
        }
        if protein.contains('*') {
            warnings.push(
                "Translated protein contains an internal stop codon after trimming any terminal stop."
                    .to_string(),
            );
        }
        if protein.contains('X') {
            warnings.push(
                "Translated protein contains ambiguous codon(s) that were rendered as 'X'."
                    .to_string(),
            );
        }
        (protein, terminal_stop_trimmed, warnings)
    }

    fn build_transcript_protein_derivation(
        derived_sequence: &str,
        source_feature: &gb_io::seq::Feature,
        source_feature_id: usize,
        source_seq_id: &str,
        source_features: &[gb_io::seq::Feature],
        exon_ranges_0based: &[(usize, usize)],
        exon_segments_forward: &[(usize, usize, usize, usize)],
        is_reverse: bool,
        transcript_id: &str,
        transcript_label: &str,
    ) -> Result<Option<TranscriptProteinDerivation>, EngineError> {
        let transcript_cds_ranges_0based =
            Self::feature_qualifier_ranges_0based(source_feature, "cds_ranges_1based");
        let matching_cds_features = Self::collect_matching_cds_features_for_derivation(
            source_features,
            source_feature,
            exon_ranges_0based,
        );
        let representative_cds_feature = matching_cds_features.first().copied();
        let source_cds_ranges_0based = if !transcript_cds_ranges_0based.is_empty() {
            transcript_cds_ranges_0based
        } else {
            Self::merge_adjacent_ranges_0based(
                &matching_cds_features
                    .iter()
                    .flat_map(|feature| Self::feature_ranges_0based_for_derivation(feature))
                    .collect::<Vec<_>>(),
            )
        };
        if source_cds_ranges_0based.is_empty() {
            return Ok(None);
        }

        let source_context_feature = Self::first_source_feature_for_derivation(source_features);
        let organism = representative_cds_feature
            .and_then(|feature| Self::qualifier_text_for_derivation(feature, "organism"))
            .or_else(|| Self::qualifier_text_for_derivation(source_feature, "organism"))
            .or_else(|| {
                source_context_feature
                    .and_then(|feature| Self::qualifier_text_for_derivation(feature, "organism"))
            });
        let organelle = Self::infer_organelle_for_derivation(
            source_feature,
            representative_cds_feature,
            source_context_feature,
        );
        let speed_profile_hint = Self::infer_translation_speed_profile_hint(organism.as_deref());

        let (translation_table, translation_table_source, mut warnings) =
            Self::resolve_translation_table_for_derivation(
                source_feature,
                representative_cds_feature,
                source_context_feature,
                organelle.as_deref(),
            );
        let codon_start = Self::qualifier_usize_for_derivation(source_feature, "codon_start")
            .or_else(|| {
                representative_cds_feature.and_then(|feature| {
                    Self::qualifier_usize_for_derivation(feature, "codon_start")
                })
            })
            .or_else(|| {
                Self::qualifier_usize_for_derivation(source_feature, "phase")
                    .filter(|phase| *phase <= 2)
                    .map(|phase| phase + 1)
            })
            .or_else(|| {
                representative_cds_feature.and_then(|feature| {
                    Self::qualifier_usize_for_derivation(feature, "phase")
                        .filter(|phase| *phase <= 2)
                        .map(|phase| phase + 1)
                })
            })
            .unwrap_or(1);

        let local_cds_ranges_0based = Self::map_source_ranges_to_transcript_local_ranges_0based(
            &source_cds_ranges_0based,
            exon_segments_forward,
            is_reverse,
            derived_sequence.len(),
        );
        let trimmed_local_cds_ranges_0based = Self::trim_leading_bases_from_ranges_0based(
            &local_cds_ranges_0based,
            codon_start.saturating_sub(1),
        );
        let trimmed_local_cds_ranges_0based =
            Self::merge_adjacent_ranges_0based(&trimmed_local_cds_ranges_0based);
        if trimmed_local_cds_ranges_0based.is_empty() {
            warnings.push(
                "CDS ranges resolved, but applying codon_start/phase removed the entire coding span."
                    .to_string(),
            );
            return Ok(Some(TranscriptProteinDerivation {
                transcript_id: transcript_id.to_string(),
                transcript_label: transcript_label.to_string(),
                source_seq_id: source_seq_id.to_string(),
                source_feature_id: source_feature_id + 1,
                derivation_mode: TranscriptProteinDerivationMode::AnnotatedCds,
                cds_ranges_1based: vec![],
                cds_length_bp: 0,
                protein_sequence: String::new(),
                protein_length_aa: 0,
                translation_table,
                translation_table_label: Self::translation_table_display_name(translation_table),
                translation_table_source,
                codon_start,
                organism,
                organelle,
                translation_speed_profile_hint: speed_profile_hint,
                terminal_stop_trimmed: false,
                warnings,
            }));
        }
        let cds_ranges_1based = trimmed_local_cds_ranges_0based
            .iter()
            .map(|(start, end)| (start.saturating_add(1), *end))
            .collect::<Vec<_>>();
        let cds_length_bp = trimmed_local_cds_ranges_0based
            .iter()
            .map(|(start, end)| end.saturating_sub(*start))
            .sum::<usize>();
        let (protein_sequence, terminal_stop_trimmed, translation_warnings) =
            Self::translate_transcript_cds_to_protein(
                derived_sequence,
                &trimmed_local_cds_ranges_0based,
                translation_table,
            );
        warnings.extend(translation_warnings);
        Ok(Some(TranscriptProteinDerivation {
            transcript_id: transcript_id.to_string(),
            transcript_label: transcript_label.to_string(),
            source_seq_id: source_seq_id.to_string(),
            source_feature_id: source_feature_id + 1,
            derivation_mode: TranscriptProteinDerivationMode::AnnotatedCds,
            cds_ranges_1based,
            cds_length_bp,
            protein_sequence: protein_sequence.clone(),
            protein_length_aa: protein_sequence.len(),
            translation_table,
            translation_table_label: Self::translation_table_display_name(translation_table),
            translation_table_source,
            codon_start,
            organism,
            organelle,
            translation_speed_profile_hint: speed_profile_hint,
            terminal_stop_trimmed,
            warnings,
        }))
    }

    fn infer_transcript_protein_derivation_without_annotation(
        derived_sequence: &str,
        source_feature: &gb_io::seq::Feature,
        source_feature_id: usize,
        source_seq_id: &str,
        source_features: &[gb_io::seq::Feature],
        transcript_id: &str,
        transcript_label: &str,
    ) -> Result<Option<TranscriptProteinDerivation>, EngineError> {
        if derived_sequence.len() < 3 {
            return Ok(None);
        }
        let source_context_feature = Self::first_source_feature_for_derivation(source_features);
        let organism = Self::qualifier_text_for_derivation(source_feature, "organism").or_else(|| {
            source_context_feature
                .and_then(|feature| Self::qualifier_text_for_derivation(feature, "organism"))
        });
        let organelle =
            Self::infer_organelle_for_derivation(source_feature, None, source_context_feature);
        let speed_profile_hint = Self::infer_translation_speed_profile_hint(organism.as_deref());
        let (translation_table, translation_table_source, mut warnings) =
            Self::resolve_translation_table_for_derivation(
                source_feature,
                None,
                source_context_feature,
                organelle.as_deref(),
            );

        #[derive(Clone, Copy)]
        struct Candidate {
            start_0based: usize,
            end_0based_exclusive: usize,
            aa_len: usize,
            has_terminal_stop: bool,
            derivation_mode: TranscriptProteinDerivationMode,
        }

        let bytes = derived_sequence.as_bytes();
        let mut best: Option<Candidate> = None;
        let adopt_candidate = |best: &mut Option<Candidate>, candidate: Candidate| {
            if candidate.aa_len == 0 || candidate.end_0based_exclusive <= candidate.start_0based {
                return;
            }
            let should_replace = match *best {
                None => true,
                Some(current) => {
                    if current.derivation_mode != TranscriptProteinDerivationMode::InferredOrf
                        && candidate.derivation_mode == TranscriptProteinDerivationMode::InferredOrf
                    {
                        true
                    } else if candidate.has_terminal_stop && !current.has_terminal_stop {
                        true
                    } else if candidate.aa_len > current.aa_len {
                        true
                    } else {
                        candidate.aa_len == current.aa_len
                            && candidate.start_0based < current.start_0based
                    }
                }
            };
            if should_replace {
                *best = Some(candidate);
            }
        };

        for start_0based in 0..=derived_sequence.len().saturating_sub(3) {
            if !bytes[start_0based..].starts_with(b"ATG") {
                continue;
            }
            let mut pos = start_0based;
            let mut has_terminal_stop = false;
            while pos + 3 <= bytes.len() {
                let aa = AMINO_ACIDS.codon2aa(
                    [bytes[pos], bytes[pos + 1], bytes[pos + 2]],
                    Some(translation_table),
                );
                pos += 3;
                if aa == STOP_CODON {
                    has_terminal_stop = true;
                    break;
                }
            }
            let translated_end = if has_terminal_stop {
                pos
            } else {
                start_0based + ((bytes.len().saturating_sub(start_0based)) / 3) * 3
            };
            let aa_len = translated_end.saturating_sub(start_0based) / 3
                - usize::from(has_terminal_stop);
            adopt_candidate(&mut best, Candidate {
                start_0based,
                end_0based_exclusive: translated_end,
                aa_len,
                has_terminal_stop,
                derivation_mode: TranscriptProteinDerivationMode::InferredOrf,
            });
        }

        if best.is_none() {
            for frame in 0..3usize {
                let mut segment_start = frame;
                let mut pos = frame;
                while pos + 3 <= bytes.len() {
                    let aa = AMINO_ACIDS.codon2aa(
                        [bytes[pos], bytes[pos + 1], bytes[pos + 2]],
                        Some(translation_table),
                    );
                    if aa == STOP_CODON {
                        let aa_len = pos.saturating_sub(segment_start) / 3;
                        adopt_candidate(&mut best, Candidate {
                            start_0based: segment_start,
                            end_0based_exclusive: pos,
                            aa_len,
                            has_terminal_stop: false,
                            derivation_mode: TranscriptProteinDerivationMode::HeuristicLongestFrame,
                        });
                        segment_start = pos + 3;
                    }
                    pos += 3;
                }
                let frame_end = frame + ((bytes.len().saturating_sub(frame)) / 3) * 3;
                let aa_len = frame_end.saturating_sub(segment_start) / 3;
                adopt_candidate(&mut best, Candidate {
                    start_0based: segment_start,
                    end_0based_exclusive: frame_end,
                    aa_len,
                    has_terminal_stop: false,
                    derivation_mode: TranscriptProteinDerivationMode::HeuristicLongestFrame,
                });
            }
        }

        let Some(best) = best else {
            return Ok(None);
        };
        let cds_ranges_0based = vec![(best.start_0based, best.end_0based_exclusive)];
        let cds_ranges_1based = cds_ranges_0based
            .iter()
            .map(|(start, end)| (start.saturating_add(1), *end))
            .collect::<Vec<_>>();
        let cds_length_bp = cds_ranges_0based
            .iter()
            .map(|(start, end)| end.saturating_sub(*start))
            .sum::<usize>();
        let (protein_sequence, terminal_stop_trimmed, translation_warnings) =
            Self::translate_transcript_cds_to_protein(
                derived_sequence,
                &cds_ranges_0based,
                translation_table,
            );
        warnings.extend(translation_warnings);
        warnings.push(match best.derivation_mode {
            TranscriptProteinDerivationMode::InferredOrf => format!(
                "CDS annotation was absent; inferred a forward ORF starting at transcript position {}{}.",
                best.start_0based.saturating_add(1),
                if best.has_terminal_stop {
                    " and ending at the first in-frame stop codon"
                } else {
                    " without a terminating in-frame stop codon"
                }
            ),
            TranscriptProteinDerivationMode::HeuristicLongestFrame => format!(
                "CDS annotation and ATG-start ORFs were absent; used the longest stop-free reading frame segment starting at transcript position {}.",
                best.start_0based.saturating_add(1)
            ),
            TranscriptProteinDerivationMode::AnnotatedCds => unreachable!(),
        });

        Ok(Some(TranscriptProteinDerivation {
            transcript_id: transcript_id.to_string(),
            transcript_label: transcript_label.to_string(),
            source_seq_id: source_seq_id.to_string(),
            source_feature_id: source_feature_id + 1,
            derivation_mode: best.derivation_mode,
            cds_ranges_1based,
            cds_length_bp,
            protein_sequence: protein_sequence.clone(),
            protein_length_aa: protein_sequence.len(),
            translation_table,
            translation_table_label: Self::translation_table_display_name(translation_table),
            translation_table_source,
            codon_start: 1,
            organism,
            organelle,
            translation_speed_profile_hint: speed_profile_hint,
            terminal_stop_trimmed,
            warnings,
        }))
    }

    fn infer_translation_speed_profile_enum(raw: Option<&str>) -> Option<TranslationSpeedProfile> {
        match raw?.trim().to_ascii_lowercase().as_str() {
            "human" => Some(TranslationSpeedProfile::Human),
            "mouse" => Some(TranslationSpeedProfile::Mouse),
            "yeast" => Some(TranslationSpeedProfile::Yeast),
            "ecoli" | "e_coli" | "e-coli" => Some(TranslationSpeedProfile::Ecoli),
            _ => None,
        }
    }

    fn codon_profile_species_label(
        profile: TranslationSpeedProfile,
    ) -> (&'static str, Option<&'static str>) {
        match profile {
            TranslationSpeedProfile::Human => ("Human", None),
            TranslationSpeedProfile::Mouse => (
                "Rattus norvegicus",
                Some("Mouse codon-speed bias currently uses the bundled rat codon-preference proxy because a dedicated Mus musculus table is not bundled yet."),
            ),
            TranslationSpeedProfile::Yeast => ("Saccharomyces cerevisiae", None),
            TranslationSpeedProfile::Ecoli => ("E. coli", None),
        }
    }

    fn sequence_feature_translation_speed_hint(
        dna: &DNAsequence,
    ) -> Option<TranslationSpeedProfile> {
        dna.features()
            .iter()
            .find_map(|feature| {
                Self::feature_qualifier_text(feature, "translation_speed_profile_hint")
                    .or_else(|| Self::feature_qualifier_text(feature, "speed_profile"))
            })
            .and_then(|raw| Self::infer_translation_speed_profile_enum(Some(raw.as_str())))
    }

    fn choose_back_translation_codon(
        protein_char: char,
        translation_table: usize,
        preferred_species_label: Option<&str>,
        speed_mark: Option<TranslationSpeedMark>,
        target_anneal_tm_c: Option<f64>,
        anneal_window_bp: usize,
        current_dna: &str,
    ) -> Option<String> {
        let aa = if protein_char == '*' {
            STOP_CODON
        } else {
            protein_char.to_ascii_uppercase()
        };
        let mut candidates = AMINO_ACIDS
            .aa2codons(aa, Some(translation_table))
            .into_iter()
            .map(|codon| String::from_utf8_lossy(&codon).to_string())
            .collect::<Vec<_>>();
        candidates.sort();
        if candidates.is_empty() {
            return None;
        }
        let preferred = preferred_species_label
            .and_then(|label| AMINO_ACIDS.preferred_species_codon(aa, label))
            .or_else(|| AMINO_ACIDS.preferred_species_codon(aa, "Default"));
        let multiple_candidates = candidates.len() > 1;
        candidates.sort_by(|left, right| {
            let score = |candidate: &str| {
                let is_preferred = preferred.as_deref() == Some(candidate);
                let preference_penalty = match speed_mark {
                    Some(TranslationSpeedMark::Fast) => {
                        if is_preferred { 0.0 } else { 10.0 }
                    }
                    Some(TranslationSpeedMark::Slow) => {
                        if !is_preferred && multiple_candidates { 0.0 } else { 10.0 }
                    }
                    None => {
                        if is_preferred { 0.0 } else { 1.0 }
                    }
                };
                let tm_penalty = if let Some(target_tm) = target_anneal_tm_c {
                    let appended = format!("{current_dna}{candidate}");
                    let window_len = anneal_window_bp.max(6).min(appended.len());
                    let suffix = &appended.as_bytes()[appended.len().saturating_sub(window_len)..];
                    (Self::estimate_primer_tm_c(suffix) - target_tm).abs() / 10.0
                } else {
                    0.0
                };
                let gc_penalty = if target_anneal_tm_c.is_some() {
                    let gc = Self::sequence_gc_fraction(candidate.as_bytes()).unwrap_or(0.0);
                    (gc - 0.5).abs()
                } else {
                    0.0
                };
                (
                    preference_penalty + tm_penalty + gc_penalty,
                    candidate.to_string(),
                )
            };
            score(left)
                .partial_cmp(&score(right))
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        candidates.into_iter().next()
    }

    fn build_reverse_translated_coding_sequence(
        protein_sequence: &str,
        translation_table: usize,
        preferred_species_label: Option<&str>,
        speed_mark: Option<TranslationSpeedMark>,
        target_anneal_tm_c: Option<f64>,
        anneal_window_bp: usize,
    ) -> (String, Vec<String>) {
        let mut dna = String::with_capacity(protein_sequence.len().saturating_mul(3));
        let mut warnings = vec![];
        for (aa_index, aa) in protein_sequence.chars().enumerate() {
            match Self::choose_back_translation_codon(
                aa,
                translation_table,
                preferred_species_label,
                speed_mark,
                target_anneal_tm_c,
                anneal_window_bp,
                &dna,
            ) {
                Some(codon) => dna.push_str(&codon),
                None => {
                    dna.push_str("NNN");
                    warnings.push(format!(
                        "Protein residue {} ('{}') has no bundled codon-choice rule for translation table {}; used 'NNN'.",
                        aa_index + 1,
                        aa,
                        translation_table
                    ));
                }
            }
        }
        (dna, warnings)
    }

    fn build_transcript_derived_protein_sequence(
        protein_sequence: &str,
        seq_name: &str,
        source_seq_id: &str,
        source_feature_id: usize,
        transcript_feature: &gb_io::seq::Feature,
        representative_cds_feature: Option<&gb_io::seq::Feature>,
        derivation: &TranscriptProteinDerivation,
    ) -> Result<DNAsequence, EngineError> {
        let mut protein = DNAsequence::from_sequence(protein_sequence).map_err(|e| EngineError {
            code: ErrorCode::Internal,
            message: format!("Could not construct protein sequence '{seq_name}': {e}"),
        })?;
        protein.set_name(seq_name);
        protein.set_molecule_type("protein");
        let mut qualifiers = vec![
            ("label".into(), Some(seq_name.to_string())),
            ("source_seq_id".into(), Some(source_seq_id.to_string())),
            (
                "source_feature_id".into(),
                Some((source_feature_id + 1).to_string()),
            ),
            ("transcript_id".into(), Some(derivation.transcript_id.clone())),
            (
                "protein_derivation_mode".into(),
                Some(derivation.derivation_mode.as_str().to_string()),
            ),
            (
                "translation_table".into(),
                Some(derivation.translation_table.to_string()),
            ),
            (
                "translation_table_label".into(),
                Some(derivation.translation_table_label.clone()),
            ),
            (
                "translation_table_source".into(),
                Some(derivation.translation_table_source.as_str().to_string()),
            ),
            (
                "synthetic_origin".into(),
                Some("transcript_protein_derived".to_string()),
            ),
        ];
        if let Some(speed_hint) = derivation.translation_speed_profile_hint.as_ref() {
            qualifiers.push((
                "translation_speed_profile_hint".into(),
                Some(speed_hint.clone()),
            ));
        }
        if let Some(organism) = derivation.organism.as_ref() {
            qualifiers.push(("translation_context_organism".into(), Some(organism.clone())));
        }
        if let Some(organelle) = derivation.organelle.as_ref() {
            qualifiers.push(("translation_context_organelle".into(), Some(organelle.clone())));
        }
        for key in ["gene", "gene_id", "locus_tag", "product", "protein_id", "note"] {
            if let Some(value) = representative_cds_feature
                .and_then(|feature| Self::qualifier_text_for_derivation(feature, key))
                .or_else(|| Self::qualifier_text_for_derivation(transcript_feature, key))
            {
                qualifiers.push((key.into(), Some(value)));
            }
        }
        let protein_len = protein.len();
        if protein_len > 0 {
            protein.features_mut().push(gb_io::seq::Feature {
                kind: "Protein".into(),
                location: gb_io::seq::Location::simple_range(0, protein_len as i64),
                qualifiers,
            });
        }
        Self::prepare_sequence(&mut protein);
        Ok(protein)
    }

    fn build_reverse_translated_coding_dna(
        coding_sequence: &str,
        seq_name: &str,
        protein_seq_id: &str,
        protein: &DNAsequence,
        translation_table: usize,
        translation_table_label: &str,
        speed_profile: Option<TranslationSpeedProfile>,
        preferred_species_label: Option<&str>,
        speed_mark: Option<TranslationSpeedMark>,
        target_anneal_tm_c: Option<f64>,
        anneal_window_bp: usize,
        warnings: &[String],
    ) -> Result<DNAsequence, EngineError> {
        let mut dna = DNAsequence::from_sequence(coding_sequence).map_err(|e| EngineError {
            code: ErrorCode::Internal,
            message: format!("Could not construct coding DNA sequence '{seq_name}': {e}"),
        })?;
        dna.set_name(seq_name);
        dna.set_molecule_type("dsDNA");
        let mut qualifiers = vec![
            ("label".into(), Some(seq_name.to_string())),
            ("protein_seq_id".into(), Some(protein_seq_id.to_string())),
            ("translation".into(), Some(protein.get_forward_string())),
            (
                "translation_table".into(),
                Some(translation_table.to_string()),
            ),
            (
                "translation_table_label".into(),
                Some(translation_table_label.to_string()),
            ),
            (
                "synthetic_origin".into(),
                Some("protein_reverse_translated".to_string()),
            ),
        ];
        if let Some(profile) = speed_profile {
            qualifiers.push((
                "translation_speed_profile_hint".into(),
                Some(profile.as_str().to_string()),
            ));
        }
        if let Some(species) = preferred_species_label {
            qualifiers.push((
                "translation_speed_reference_species".into(),
                Some(species.to_string()),
            ));
        }
        if let Some(mark) = speed_mark {
            qualifiers.push(("translation_speed_mark".into(), Some(mark.as_str().to_string())));
        }
        if let Some(target_tm) = target_anneal_tm_c {
            qualifiers.push((
                "target_anneal_tm_c".into(),
                Some(format!("{target_tm:.2}")),
            ));
            qualifiers.push((
                "target_anneal_window_bp".into(),
                Some(anneal_window_bp.to_string()),
            ));
        }
        for warning in warnings {
            qualifiers.push(("reverse_translation_warning".into(), Some(warning.clone())));
        }
        let dna_len = dna.len();
        if dna_len > 0 {
            dna.features_mut().push(gb_io::seq::Feature {
                kind: "CDS".into(),
                location: gb_io::seq::Location::simple_range(0, dna_len as i64),
                qualifiers,
            });
        }
        Self::prepare_sequence(&mut dna);
        Ok(dna)
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
        source_features: &[gb_io::seq::Feature],
        source_feature_id: usize,
        source_seq_id: &str,
    ) -> Result<
        (
            DNAsequence,
            String,
            String,
            bool,
            usize,
            Option<TranscriptProteinDerivation>,
        ),
        EngineError,
    > {
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
        let mut exon_segments_forward: Vec<(usize, usize, usize, usize)> = vec![];
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
            exon_segments_forward.push((
                *start_0based,
                *end_0based_exclusive,
                local_start,
                local_end,
            ));
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
        let mut exon_segments = exon_segments_forward.clone();
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
        let protein_derivation = Self::build_transcript_protein_derivation(
            &derived_sequence,
            source_feature,
            source_feature_id,
            source_seq_id,
            source_features,
            &exon_ranges,
            &exon_segments_forward,
            is_reverse,
            &transcript_id,
            &transcript_label,
        )?;
        let representative_cds_feature = Self::collect_matching_cds_features_for_derivation(
            source_features,
            source_feature,
            &exon_ranges,
        )
        .into_iter()
        .next();

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
        if let Some(derivation) = protein_derivation.as_ref() {
            if let Some(cds_encoded) = Self::serialize_ranges_1based(&derivation.cds_ranges_1based)
            {
                transcript_qualifiers.push(("cds_ranges_1based".into(), Some(cds_encoded)));
            }
            transcript_qualifiers.push((
                "protein_length_aa".into(),
                Some(derivation.protein_length_aa.to_string()),
            ));
            transcript_qualifiers.push((
                "translation_table".into(),
                Some(derivation.translation_table.to_string()),
            ));
            transcript_qualifiers.push((
                "translation_table_label".into(),
                Some(derivation.translation_table_label.clone()),
            ));
            transcript_qualifiers.push((
                "translation_table_source".into(),
                Some(derivation.translation_table_source.as_str().to_string()),
            ));
            transcript_qualifiers.push((
                "protein_derivation_mode".into(),
                Some(derivation.derivation_mode.as_str().to_string()),
            ));
            if let Some(organism) = derivation.organism.as_ref() {
                transcript_qualifiers.push((
                    "translation_context_organism".into(),
                    Some(organism.clone()),
                ));
            }
            if let Some(organelle) = derivation.organelle.as_ref() {
                transcript_qualifiers.push((
                    "translation_context_organelle".into(),
                    Some(organelle.clone()),
                ));
            }
            if let Some(speed_hint) = derivation.translation_speed_profile_hint.as_ref() {
                transcript_qualifiers.push((
                    "translation_speed_profile_hint".into(),
                    Some(speed_hint.clone()),
                ));
            }
            if !derivation.protein_sequence.is_empty() {
                transcript_qualifiers.push((
                    "derived_protein_translation".into(),
                    Some(derivation.protein_sequence.clone()),
                ));
            }
        }
        derived.features_mut().push(gb_io::seq::Feature {
            kind: "mRNA".into(),
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
                kind: "exon".into(),
                location: gb_io::seq::Location::simple_range(segment.2 as i64, segment.3 as i64),
                qualifiers: exon_qualifiers,
            });
        }

        if let Some(derivation) = protein_derivation.as_ref()
            && !derivation.cds_ranges_1based.is_empty()
        {
            let mut cds_parts = derivation
                .cds_ranges_1based
                .iter()
                .filter_map(|(start_1based, end_1based)| {
                    let start_0based = start_1based.saturating_sub(1);
                    (*end_1based > start_0based).then_some(gb_io::seq::Location::simple_range(
                        start_0based as i64,
                        *end_1based as i64,
                    ))
                })
                .collect::<Vec<_>>();
            let cds_location = if cds_parts.len() == 1 {
                cds_parts.remove(0)
            } else {
                gb_io::seq::Location::Join(cds_parts)
            };
            let mut cds_qualifiers = vec![
                ("transcript_id".into(), Some(transcript_id.clone())),
                ("label".into(), Some(format!("{transcript_label} CDS"))),
                ("source_seq_id".into(), Some(source_seq_id.to_string())),
                (
                    "source_feature_id".into(),
                    Some((source_feature_id + 1).to_string()),
                ),
                (
                    "synthetic_origin".into(),
                    Some("mrna_transcript_derived".to_string()),
                ),
                (
                    "codon_start".into(),
                    Some(derivation.codon_start.to_string()),
                ),
                (
                    "transl_table".into(),
                    Some(derivation.translation_table.to_string()),
                ),
                (
                    "translation_table_label".into(),
                    Some(derivation.translation_table_label.clone()),
                ),
                (
                    "translation_table_source".into(),
                    Some(derivation.translation_table_source.as_str().to_string()),
                ),
                (
                    "protein_derivation_mode".into(),
                    Some(derivation.derivation_mode.as_str().to_string()),
                ),
                (
                    "protein_length_aa".into(),
                    Some(derivation.protein_length_aa.to_string()),
                ),
                (
                    "translation".into(),
                    Some(derivation.protein_sequence.clone()),
                ),
            ];
            if let Some(cds_encoded) = Self::serialize_ranges_1based(&derivation.cds_ranges_1based)
            {
                cds_qualifiers.push(("cds_ranges_1based".into(), Some(cds_encoded)));
            }
            if derivation.terminal_stop_trimmed {
                cds_qualifiers.push(("terminal_stop_trimmed".into(), Some("true".to_string())));
            }
            if let Some(organism) = derivation.organism.as_ref() {
                cds_qualifiers.push((
                    "translation_context_organism".into(),
                    Some(organism.clone()),
                ));
            }
            if let Some(organelle) = derivation.organelle.as_ref() {
                cds_qualifiers.push((
                    "translation_context_organelle".into(),
                    Some(organelle.clone()),
                ));
            }
            if let Some(speed_hint) = derivation.translation_speed_profile_hint.as_ref() {
                cds_qualifiers.push((
                    "translation_speed_profile_hint".into(),
                    Some(speed_hint.clone()),
                ));
            }
            for key in [
                "gene",
                "gene_id",
                "locus_tag",
                "product",
                "protein_id",
                "note",
            ] {
                if let Some(value) = representative_cds_feature
                    .and_then(|feature| Self::qualifier_text_for_derivation(feature, key))
                    .or_else(|| Self::qualifier_text_for_derivation(source_feature, key))
                {
                    cds_qualifiers.push((key.into(), Some(value)));
                }
            }
            for warning in &derivation.warnings {
                cds_qualifiers.push(("translation_warning".into(), Some(warning.clone())));
            }
            derived.features_mut().push(gb_io::seq::Feature {
                kind: "CDS".into(),
                location: cds_location,
                qualifiers: cds_qualifiers,
            });
        }

        Self::prepare_sequence(&mut derived);
        Ok((
            derived,
            transcript_id,
            transcript_label,
            is_reverse,
            exon_segments.len(),
            protein_derivation,
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
                gene_extract_mode: None,
                promoter_upstream_bp: None,
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

    fn build_dbsnp_variant_marker_feature(
        rs_id: &str,
        chromosome: &str,
        chromosome_display: &str,
        position_1based: usize,
        assembly_name: Option<&str>,
        gene_symbols: &[String],
        local_start_0based: usize,
    ) -> gb_io::seq::Feature {
        let chromosome_label = if chromosome_display.trim().is_empty() {
            chromosome.trim()
        } else {
            chromosome_display.trim()
        };
        let assembly_label = assembly_name
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .unwrap_or("unknown");
        let accession_note = if chromosome.trim().is_empty()
            || chromosome.trim().eq_ignore_ascii_case(chromosome_label)
        {
            String::new()
        } else {
            format!(" [accession={}]", chromosome.trim())
        };
        let mut note = format!(
            "dbSNP variant {rs_id} resolved from NCBI Variation to {chromosome_label}:{position_1based}{accession_note} (assembly={assembly_label})"
        );
        if !gene_symbols.is_empty() {
            note.push_str(&format!(" [genes={}]", gene_symbols.join(", ")));
        }
        let mut qualifiers = vec![
            ("label".into(), Some(rs_id.to_string())),
            ("note".into(), Some(note)),
            ("db_xref".into(), Some(format!("dbSNP:{rs_id}"))),
            ("chromosome".into(), Some(chromosome.to_string())),
            (
                "genomic_position_1based".into(),
                Some(position_1based.to_string()),
            ),
            ("assembly_name".into(), Some(assembly_label.to_string())),
            (
                "gentle_generated".into(),
                Some(DBSNP_VARIANT_MARKER_GENERATED_TAG.to_string()),
            ),
        ];
        if !chromosome.trim().is_empty()
            && !chromosome.trim().eq_ignore_ascii_case(chromosome_label)
        {
            qualifiers.push(("refseq_accession".into(), Some(chromosome.to_string())));
        }
        gb_io::seq::Feature {
            kind: "variation".into(),
            location: gb_io::seq::Location::simple_range(
                local_start_0based as i64,
                local_start_0based.saturating_add(1) as i64,
            ),
            qualifiers,
        }
    }

    fn normalize_primer_insertion_intent(
        insertion: &PrimerInsertionIntent,
        template_len: usize,
    ) -> Result<PrimerInsertionIntent, EngineError> {
        if template_len == 0 {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "Insertion intent requires a non-empty template sequence".to_string(),
            });
        }
        if insertion.forward_window_start_0based >= insertion.forward_window_end_0based_exclusive
            || insertion.forward_window_end_0based_exclusive > template_len
        {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "insertion.forward window {}..{} is invalid for template length {}",
                    insertion.forward_window_start_0based,
                    insertion.forward_window_end_0based_exclusive,
                    template_len
                ),
            });
        }
        if insertion.reverse_window_start_0based >= insertion.reverse_window_end_0based_exclusive
            || insertion.reverse_window_end_0based_exclusive > template_len
        {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "insertion.reverse window {}..{} is invalid for template length {}",
                    insertion.reverse_window_start_0based,
                    insertion.reverse_window_end_0based_exclusive,
                    template_len
                ),
            });
        }
        if insertion.requested_forward_3prime_end_0based_exclusive == 0
            || insertion.requested_forward_3prime_end_0based_exclusive > template_len
        {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "insertion.requested_forward_3prime_end_0based_exclusive ({}) must be in 1..={}",
                    insertion.requested_forward_3prime_end_0based_exclusive, template_len
                ),
            });
        }
        if insertion.requested_reverse_3prime_start_0based >= template_len {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "insertion.requested_reverse_3prime_start_0based ({}) must be < {}",
                    insertion.requested_reverse_3prime_start_0based, template_len
                ),
            });
        }
        if insertion.requested_forward_3prime_end_0based_exclusive
            > insertion.requested_reverse_3prime_start_0based
        {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "Insertion anchor geometry is invalid: requested forward 3' end must be <= requested reverse 3' start".to_string(),
            });
        }
        let mut normalized = insertion.clone();
        normalized.forward_extension_5prime =
            Self::normalize_iupac_text(&normalized.forward_extension_5prime)?;
        normalized.reverse_extension_5prime =
            Self::normalize_iupac_text(&normalized.reverse_extension_5prime)?;
        normalized.max_anchor_shift_bp = Some(normalized.max_anchor_shift_bp.unwrap_or(12));
        Ok(normalized)
    }

    fn derive_insertion_intent_roi_bounds(
        insertion: &PrimerInsertionIntent,
        template_len: usize,
    ) -> (usize, usize) {
        let left = insertion.requested_forward_3prime_end_0based_exclusive;
        let right = insertion.requested_reverse_3prime_start_0based;
        if left < right {
            (left, right)
        } else if left < template_len {
            (left, (left + 1).min(template_len))
        } else {
            (template_len.saturating_sub(1), template_len)
        }
    }

    fn build_primer_insertion_context_report(
        template_seq: &str,
        insertion: &PrimerInsertionIntent,
        pairs: &[PrimerDesignPairRecord],
    ) -> PrimerInsertionContextReport {
        let max_shift_bp = insertion.max_anchor_shift_bp.unwrap_or(12);
        let mut rows: Vec<PrimerInsertionPairCompensation> = Vec::with_capacity(pairs.len());
        let mut uncompensable_pair_count = 0usize;
        let mut out_of_shift_budget_pair_count = 0usize;
        for (pair_index, pair) in pairs.iter().enumerate() {
            let rank = if pair.rank == 0 {
                pair_index + 1
            } else {
                pair.rank
            };
            let actual_forward_end = pair.forward.end_0based_exclusive;
            let actual_reverse_start = pair.reverse.start_0based;
            let forward_shift_bp = insertion.requested_forward_3prime_end_0based_exclusive as isize
                - actual_forward_end as isize;
            let reverse_shift_bp = actual_reverse_start as isize
                - insertion.requested_reverse_3prime_start_0based as isize;
            let within_shift_budget = forward_shift_bp.unsigned_abs() <= max_shift_bp
                && reverse_shift_bp.unsigned_abs() <= max_shift_bp;
            let compensable = forward_shift_bp >= 0 && reverse_shift_bp >= 0;
            if !compensable {
                uncompensable_pair_count = uncompensable_pair_count.saturating_add(1);
            }
            if !within_shift_budget {
                out_of_shift_budget_pair_count = out_of_shift_budget_pair_count.saturating_add(1);
            }

            let forward_compensation_5prime = if forward_shift_bp > 0 {
                template_seq
                    [actual_forward_end..insertion.requested_forward_3prime_end_0based_exclusive]
                    .to_string()
            } else {
                String::new()
            };
            let reverse_compensation_5prime = if reverse_shift_bp > 0 {
                Self::reverse_complement(
                    &template_seq
                        [insertion.requested_reverse_3prime_start_0based..actual_reverse_start],
                )
            } else {
                String::new()
            };

            let compensated_forward_5prime_tail = format!(
                "{}{}",
                insertion.forward_extension_5prime, forward_compensation_5prime
            );
            let compensated_reverse_5prime_tail = format!(
                "{}{}",
                insertion.reverse_extension_5prime, reverse_compensation_5prime
            );

            let forward_seq = pair.forward.sequence.to_ascii_uppercase();
            let reverse_seq = pair.reverse.sequence.to_ascii_uppercase();
            let forward_anneal_len = pair.forward.anneal_length_bp.min(forward_seq.len());
            let reverse_anneal_len = pair.reverse.anneal_length_bp.min(reverse_seq.len());
            let forward_anneal =
                forward_seq[forward_seq.len().saturating_sub(forward_anneal_len)..].to_string();
            let reverse_anneal =
                reverse_seq[reverse_seq.len().saturating_sub(reverse_anneal_len)..].to_string();

            rows.push(PrimerInsertionPairCompensation {
                rank,
                forward_anchor_shift_bp: forward_shift_bp,
                reverse_anchor_shift_bp: reverse_shift_bp,
                within_shift_budget,
                compensable,
                forward_compensation_5prime,
                reverse_compensation_5prime,
                compensated_forward_5prime_tail: compensated_forward_5prime_tail.clone(),
                compensated_reverse_5prime_tail: compensated_reverse_5prime_tail.clone(),
                compensated_forward_sequence: format!(
                    "{compensated_forward_5prime_tail}{forward_anneal}"
                ),
                compensated_reverse_sequence: format!(
                    "{compensated_reverse_5prime_tail}{reverse_anneal}"
                ),
            });
        }
        PrimerInsertionContextReport {
            requested_forward_3prime_end_0based_exclusive: insertion
                .requested_forward_3prime_end_0based_exclusive,
            requested_reverse_3prime_start_0based: insertion.requested_reverse_3prime_start_0based,
            forward_extension_5prime: insertion.forward_extension_5prime.clone(),
            reverse_extension_5prime: insertion.reverse_extension_5prime.clone(),
            forward_window_start_0based: insertion.forward_window_start_0based,
            forward_window_end_0based_exclusive: insertion.forward_window_end_0based_exclusive,
            reverse_window_start_0based: insertion.reverse_window_start_0based,
            reverse_window_end_0based_exclusive: insertion.reverse_window_end_0based_exclusive,
            max_anchor_shift_bp: max_shift_bp,
            uncompensable_pair_count,
            out_of_shift_budget_pair_count,
            pairs: rows,
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn execute_design_primer_pairs(
        &mut self,
        result: &mut OpResult,
        parent_seq_ids: &mut Vec<SeqId>,
        template: SeqId,
        roi_start_0based: usize,
        roi_end_0based: usize,
        forward: PrimerDesignSideConstraint,
        reverse: PrimerDesignSideConstraint,
        pair_constraints: PrimerDesignPairConstraint,
        min_amplicon_bp: usize,
        max_amplicon_bp: usize,
        max_tm_delta_c: Option<f64>,
        max_pairs: Option<usize>,
        report_id: Option<String>,
        insertion_intent: Option<PrimerInsertionIntent>,
    ) -> Result<(), EngineError> {
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
                message: "DesignPrimerPairs currently supports linear templates only".to_string(),
            });
        }

        let template_seq = dna.get_forward_string().to_ascii_uppercase();
        let template_bytes = template_seq.as_bytes();
        if template_bytes.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "DesignPrimerPairs requires a non-empty template sequence".to_string(),
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
                    Ok((pairs, rejection_summary, version, explain, request_boulder_io)) => {
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

        let insertion_context = insertion_intent.as_ref().map(|insertion| {
            Self::build_primer_insertion_context_report(&template_seq, insertion, &pairs)
        });
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
            insertion_context,
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
        if let Some(context) = report.insertion_context.as_ref() {
            result.messages.push(format!(
                "Insertion intent annotated for report '{}': requested anchors fwd_end={} rev_start={} with max shift {} bp",
                report.report_id,
                context.requested_forward_3prime_end_0based_exclusive,
                context.requested_reverse_3prime_start_0based,
                context.max_anchor_shift_bp
            ));
            if context.uncompensable_pair_count > 0 {
                result.warnings.push(format!(
                    "{} pair(s) in report '{}' require non-compensable downstream anchor shifts; review insertion compensation rows",
                    context.uncompensable_pair_count,
                    report.report_id
                ));
            }
            if context.out_of_shift_budget_pair_count > 0 {
                result.warnings.push(format!(
                    "{} pair(s) in report '{}' exceed max_anchor_shift_bp={}; review insertion compensation rows",
                    context.out_of_shift_budget_pair_count,
                    report.report_id,
                    context.max_anchor_shift_bp
                ));
            }
        }
        if !report.pairs.is_empty() {
            parent_seq_ids.push(template.clone());
            let report_token = Self::normalize_id_token(&report.report_id);
            let mut pair_container_ids: Vec<String> = Vec::with_capacity(report.pairs.len());
            let mut amplicon_seq_count = 0usize;
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
                let amplicon_seq_id = self.unique_seq_id(&format!("{base}_amplicon"));

                let mut forward_primer = DNAsequence::from_sequence(&pair.forward.sequence).map_err(|e| EngineError {
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

                let mut reverse_primer = DNAsequence::from_sequence(&pair.reverse.sequence).map_err(|e| EngineError {
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

                let amplicon_sequence = Self::build_tailed_amplicon_sequence_from_primer_pair(
                    &template_seq,
                    pair,
                )
                .map_err(|e| EngineError {
                    code: e.code,
                    message: format!(
                        "Could not derive predicted amplicon sequence for report '{}' pair {}: {}",
                        report.report_id, pair_rank, e.message
                    ),
                })?;
                let mut amplicon =
                    DNAsequence::from_sequence(&amplicon_sequence).map_err(|e| EngineError {
                        code: ErrorCode::Internal,
                        message: format!(
                            "Could not materialize amplicon sequence for report '{}' pair {}: {e}",
                            report.report_id, pair_rank
                        ),
                    })?;
                amplicon.set_circular(false);
                amplicon.set_name(amplicon_seq_id.clone());
                Self::prepare_sequence(&mut amplicon);
                self.state
                    .sequences
                    .insert(amplicon_seq_id.clone(), amplicon);
                self.add_lineage_node(
                    &amplicon_seq_id,
                    SequenceOrigin::Derived,
                    Some(&result.op_id),
                );
                result.created_seq_ids.push(amplicon_seq_id.clone());
                amplicon_seq_count = amplicon_seq_count.saturating_add(1);

                let pair_members = vec![
                    forward_seq_id.clone(),
                    reverse_seq_id.clone(),
                    amplicon_seq_id.clone(),
                ];
                if let Some(container_id) = self.add_container(
                    &pair_members,
                    ContainerKind::Pool,
                    Some(format!("Primer pair {} {}", report.report_id, pair_token)),
                    Some(&result.op_id),
                ) {
                    pair_container_ids.push(container_id.clone());
                    result.messages.push(format!(
                        "Created primer-pair container '{}' for {} ({}, {}, {})",
                        container_id, pair_token, forward_seq_id, reverse_seq_id, amplicon_seq_id
                    ));
                }
            }
            result.messages.push(format!(
                "Materialized {} primer sequence(s), {} predicted amplicon sequence(s), and {} primer-pair container(s) for report '{}'",
                report.pairs.len().saturating_mul(2),
                amplicon_seq_count,
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
        Ok(())
    }

    fn materialize_derived_linear_sequence(
        &mut self,
        result: &mut OpResult,
        base_id: &str,
        sequence: &str,
    ) -> Result<String, EngineError> {
        let mut dna = DNAsequence::from_sequence(sequence).map_err(|e| EngineError {
            code: ErrorCode::Internal,
            message: format!("Could not materialize derived sequence '{}': {e}", base_id),
        })?;
        dna.set_circular(false);
        let seq_id = self.unique_seq_id(base_id);
        dna.set_name(seq_id.clone());
        Self::prepare_sequence(&mut dna);
        self.state.sequences.insert(seq_id.clone(), dna);
        self.add_lineage_node(&seq_id, SequenceOrigin::Derived, Some(&result.op_id));
        result.created_seq_ids.push(seq_id.clone());
        Ok(seq_id)
    }

    #[allow(clippy::too_many_arguments)]
    fn execute_overlap_extension_mutagenesis(
        &mut self,
        result: &mut OpResult,
        parent_seq_ids: &mut Vec<SeqId>,
        template: SeqId,
        edit_start_0based: usize,
        edit_end_0based_exclusive: usize,
        insert_sequence: String,
        constraints: OverlapExtensionMutagenesisConstraints,
        output_prefix: Option<String>,
    ) -> Result<(), EngineError> {
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
                message: "PcrOverlapExtensionMutagenesis currently supports linear templates only"
                    .to_string(),
            });
        }
        let template_seq = dna.get_forward_string().to_ascii_uppercase();
        let template_bytes = template_seq.as_bytes();
        let template_len = template_bytes.len();
        if template_len == 0 {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "PcrOverlapExtensionMutagenesis requires a non-empty template sequence"
                    .to_string(),
            });
        }
        if edit_start_0based > edit_end_0based_exclusive || edit_end_0based_exclusive > template_len
        {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "PcrOverlapExtensionMutagenesis edit window {}..{} is invalid for template length {}",
                    edit_start_0based, edit_end_0based_exclusive, template_len
                ),
            });
        }
        if constraints.overlap_bp == 0 {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "PcrOverlapExtensionMutagenesis requires overlap_bp >= 1".to_string(),
            });
        }

        let insert_sequence = if insert_sequence.trim().is_empty() {
            String::new()
        } else {
            Self::normalize_iupac_text(&insert_sequence)?
        };
        let deleted_bp = edit_end_0based_exclusive.saturating_sub(edit_start_0based);
        let inserted_bp = insert_sequence.len();
        if deleted_bp == 0 && inserted_bp == 0 {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "PcrOverlapExtensionMutagenesis requires an insertion/deletion/replacement (received a no-op edit)".to_string(),
            });
        }

        let mut mutated_template = String::with_capacity(template_len + inserted_bp);
        mutated_template.push_str(&template_seq[..edit_start_0based]);
        mutated_template.push_str(&insert_sequence);
        mutated_template.push_str(&template_seq[edit_end_0based_exclusive..]);
        if mutated_template.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "PcrOverlapExtensionMutagenesis cannot produce an empty mutant template"
                    .to_string(),
            });
        }
        let deleted_bp_i64 = i64::try_from(deleted_bp).map_err(|_| EngineError {
            code: ErrorCode::InvalidInput,
            message: "Deletion span is too large".to_string(),
        })?;
        let inserted_bp_i64 = i64::try_from(inserted_bp).map_err(|_| EngineError {
            code: ErrorCode::InvalidInput,
            message: "Insertion span is too large".to_string(),
        })?;
        let delta_bp_i64 = inserted_bp_i64 - deleted_bp_i64;
        let map_original_to_mutant_pos = |pos: usize| -> Option<usize> {
            if pos <= edit_start_0based {
                Some(pos)
            } else if pos >= edit_end_0based_exclusive {
                if delta_bp_i64 >= 0 {
                    pos.checked_add(delta_bp_i64 as usize)
                } else {
                    pos.checked_sub((-delta_bp_i64) as usize)
                }
            } else {
                None
            }
        };

        let mut outer_forward = constraints.outer_forward.clone();
        outer_forward.end_0based = Some(
            outer_forward
                .end_0based
                .unwrap_or(edit_start_0based)
                .min(edit_start_0based),
        );
        let mut outer_reverse = constraints.outer_reverse.clone();
        outer_reverse.start_0based = Some(
            outer_reverse
                .start_0based
                .unwrap_or(edit_end_0based_exclusive)
                .max(edit_end_0based_exclusive),
        );
        let mut inner_forward = constraints.inner_forward.clone();
        inner_forward.start_0based = Some(
            inner_forward
                .start_0based
                .unwrap_or(edit_end_0based_exclusive)
                .max(edit_end_0based_exclusive),
        );
        let mut inner_reverse = constraints.inner_reverse.clone();
        inner_reverse.end_0based = Some(
            inner_reverse
                .end_0based
                .unwrap_or(edit_start_0based)
                .min(edit_start_0based),
        );

        for (label, side) in [
            ("outer_forward", &outer_forward),
            ("outer_reverse", &outer_reverse),
            ("inner_forward", &inner_forward),
            ("inner_reverse", &inner_reverse),
        ] {
            Self::validate_primer_design_side_constraints(label, side)?;
            if let Some(location) = side.location_0based {
                if location >= template_len {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "{label}.location_0based ({location}) is outside template length {template_len}",
                        ),
                    });
                }
            }
            if let Some(start) = side.start_0based {
                if start >= template_len {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "{label}.start_0based ({start}) is outside template length {template_len}",
                        ),
                    });
                }
            }
            if let Some(end) = side.end_0based {
                if end == 0 || end > template_len {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "{label}.end_0based ({end}) must be in 1..={template_len} for this template",
                        ),
                    });
                }
            }
        }

        let outer_forward_seq_constraints =
            Self::normalize_primer_side_sequence_constraints(&outer_forward)?;
        let outer_reverse_seq_constraints =
            Self::normalize_primer_side_sequence_constraints(&outer_reverse)?;
        let inner_forward_seq_constraints =
            Self::normalize_primer_side_sequence_constraints(&inner_forward)?;
        let inner_reverse_seq_constraints =
            Self::normalize_primer_side_sequence_constraints(&inner_reverse)?;
        let inner_forward_base_tail = inner_forward_seq_constraints
            .non_annealing_5prime_tail
            .clone()
            .unwrap_or_default();
        let inner_reverse_base_tail = inner_reverse_seq_constraints
            .non_annealing_5prime_tail
            .clone()
            .unwrap_or_default();

        let mut inner_forward_generation = inner_forward.clone();
        inner_forward_generation.non_annealing_5prime_tail = None;
        let mut inner_reverse_generation = inner_reverse.clone();
        inner_reverse_generation.non_annealing_5prime_tail = None;
        let relaxed_inner_constraints = NormalizedPrimerSideSequenceConstraints::default();

        let mut rejection_summary = PrimerDesignRejectionSummary::default();
        let outer_forward_candidates = Self::generate_primer_side_candidates(
            template_bytes,
            &outer_forward,
            &outer_forward_seq_constraints,
            false,
            &mut rejection_summary,
        );
        let outer_reverse_candidates = Self::generate_primer_side_candidates(
            template_bytes,
            &outer_reverse,
            &outer_reverse_seq_constraints,
            true,
            &mut rejection_summary,
        );
        let inner_forward_candidates = Self::generate_primer_side_candidates(
            template_bytes,
            &inner_forward_generation,
            &relaxed_inner_constraints,
            false,
            &mut rejection_summary,
        );
        let inner_reverse_candidates = Self::generate_primer_side_candidates(
            template_bytes,
            &inner_reverse_generation,
            &relaxed_inner_constraints,
            true,
            &mut rejection_summary,
        );

        if outer_forward_candidates.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "No outer_forward primer candidates satisfied constraints".to_string(),
            });
        }
        if outer_reverse_candidates.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "No outer_reverse primer candidates satisfied constraints".to_string(),
            });
        }
        if inner_forward_candidates.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "No inner_forward primer candidates satisfied constraints".to_string(),
            });
        }
        if inner_reverse_candidates.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "No inner_reverse primer candidates satisfied constraints".to_string(),
            });
        }

        #[derive(Clone)]
        struct OverlapExtensionDesign {
            score: f64,
            overlap_sequence: String,
            outer_forward: PrimerDesignCandidate,
            outer_reverse: PrimerDesignCandidate,
            inner_forward: PrimerDesignCandidate,
            inner_reverse: PrimerDesignCandidate,
            inner_forward_full: String,
            inner_reverse_full: String,
            stage1_left_product: String,
            stage1_right_product: String,
            final_product: String,
        }

        let build_tailed_amplicon = |source: &str,
                                     forward_start_0based: usize,
                                     forward_anneal_len: usize,
                                     forward_full: &str,
                                     reverse_start_0based: usize,
                                     reverse_full: &str|
         -> Option<String> {
            let forward_anneal_end = forward_start_0based.checked_add(forward_anneal_len)?;
            if forward_anneal_end > reverse_start_0based {
                return None;
            }
            let interior = source.get(forward_anneal_end..reverse_start_0based)?;
            let reverse_full_rc = Self::reverse_complement(reverse_full);
            Some(format!("{forward_full}{interior}{reverse_full_rc}"))
        };
        let primer_quality =
            |sequence: &str, anneal_len: usize, tm_c: f64, gc_fraction: f64, anneal_hits: usize| {
                let metrics = Self::compute_primer_heuristic_metrics(sequence.as_bytes());
                let mut score = 1000.0;
                score -= (tm_c - 62.0).abs() * 5.0;
                score -= (gc_fraction - 0.5).abs() * 50.0;
                score -= Self::preferred_primer_length_penalty(anneal_len) * 6.0;
                score -= Self::primer_secondary_structure_penalty(metrics);
                score -= (anneal_hits.saturating_sub(1) as f64) * 25.0;
                score
            };
        let better_than = |left: &OverlapExtensionDesign, right: &OverlapExtensionDesign| -> bool {
            if (left.score - right.score).abs() > f64::EPSILON {
                return left.score > right.score;
            }
            (
                left.outer_forward.start_0based,
                left.inner_reverse.start_0based,
                left.inner_forward.start_0based,
                left.outer_reverse.start_0based,
                left.outer_forward.sequence.as_str(),
                left.inner_reverse_full.as_str(),
                left.inner_forward_full.as_str(),
                left.outer_reverse.sequence.as_str(),
            ) < (
                right.outer_forward.start_0based,
                right.inner_reverse.start_0based,
                right.inner_forward.start_0based,
                right.outer_reverse.start_0based,
                right.outer_forward.sequence.as_str(),
                right.inner_reverse_full.as_str(),
                right.inner_forward_full.as_str(),
                right.outer_reverse.sequence.as_str(),
            )
        };

        let max_evaluations = 250_000usize;
        let mut evaluated = 0usize;
        let mut evaluation_limited = false;
        let mut best_design: Option<OverlapExtensionDesign> = None;

        'inner_search: for inner_reverse_candidate in &inner_reverse_candidates {
            if inner_reverse_candidate.end_0based_exclusive > edit_start_0based {
                continue;
            }
            for inner_forward_candidate in &inner_forward_candidates {
                if inner_forward_candidate.start_0based < edit_end_0based_exclusive {
                    continue;
                }
                let Some(inner_forward_start_mut) =
                    map_original_to_mutant_pos(inner_forward_candidate.start_0based)
                else {
                    continue;
                };
                let inner_reverse_end_mut = inner_reverse_candidate.end_0based_exclusive;
                if inner_reverse_end_mut > inner_forward_start_mut {
                    continue;
                }
                let overlap_sequence =
                    &mutated_template[inner_reverse_end_mut..inner_forward_start_mut];
                if overlap_sequence.len() < constraints.overlap_bp {
                    continue;
                }

                let inner_forward_full = format!(
                    "{}{}{}",
                    inner_forward_base_tail, overlap_sequence, inner_forward_candidate.sequence
                );
                let inner_reverse_full = format!(
                    "{}{}{}",
                    inner_reverse_base_tail,
                    Self::reverse_complement(overlap_sequence),
                    inner_reverse_candidate.sequence
                );
                if !Self::primer_sequence_matches_side_constraints(
                    inner_forward_full.as_bytes(),
                    &inner_forward_seq_constraints,
                ) {
                    continue;
                }
                if !Self::primer_sequence_matches_side_constraints(
                    inner_reverse_full.as_bytes(),
                    &inner_reverse_seq_constraints,
                ) {
                    continue;
                }

                for outer_forward_candidate in &outer_forward_candidates {
                    if outer_forward_candidate.start_0based >= inner_reverse_candidate.start_0based
                    {
                        continue;
                    }
                    let outer_forward_anneal_len = outer_forward_candidate
                        .end_0based_exclusive
                        .saturating_sub(outer_forward_candidate.start_0based);
                    for outer_reverse_candidate in &outer_reverse_candidates {
                        if evaluated >= max_evaluations {
                            evaluation_limited = true;
                            break 'inner_search;
                        }
                        evaluated = evaluated.saturating_add(1);

                        if outer_reverse_candidate.start_0based
                            <= inner_forward_candidate.end_0based_exclusive
                        {
                            continue;
                        }
                        let Some(outer_reverse_start_mut) =
                            map_original_to_mutant_pos(outer_reverse_candidate.start_0based)
                        else {
                            continue;
                        };
                        if outer_forward_candidate.end_0based_exclusive > outer_reverse_start_mut {
                            continue;
                        }

                        let inner_forward_anneal_len = inner_forward_candidate
                            .end_0based_exclusive
                            .saturating_sub(inner_forward_candidate.start_0based);
                        let stage1_left_product = match build_tailed_amplicon(
                            &template_seq,
                            outer_forward_candidate.start_0based,
                            outer_forward_anneal_len,
                            &outer_forward_candidate.sequence,
                            inner_reverse_candidate.start_0based,
                            &inner_reverse_full,
                        ) {
                            Some(value) => value,
                            None => continue,
                        };
                        let stage1_right_product = match build_tailed_amplicon(
                            &template_seq,
                            inner_forward_candidate.start_0based,
                            inner_forward_anneal_len,
                            &inner_forward_full,
                            outer_reverse_candidate.start_0based,
                            &outer_reverse_candidate.sequence,
                        ) {
                            Some(value) => value,
                            None => continue,
                        };
                        let final_product = match build_tailed_amplicon(
                            &mutated_template,
                            outer_forward_candidate.start_0based,
                            outer_forward_anneal_len,
                            &outer_forward_candidate.sequence,
                            outer_reverse_start_mut,
                            &outer_reverse_candidate.sequence,
                        ) {
                            Some(value) => value,
                            None => continue,
                        };

                        let score = primer_quality(
                            &outer_forward_candidate.sequence,
                            outer_forward_anneal_len,
                            outer_forward_candidate.tm_c,
                            outer_forward_candidate.gc_fraction,
                            outer_forward_candidate.anneal_hits,
                        ) + primer_quality(
                            &outer_reverse_candidate.sequence,
                            outer_reverse_candidate
                                .end_0based_exclusive
                                .saturating_sub(outer_reverse_candidate.start_0based),
                            outer_reverse_candidate.tm_c,
                            outer_reverse_candidate.gc_fraction,
                            outer_reverse_candidate.anneal_hits,
                        ) + primer_quality(
                            &inner_forward_full,
                            inner_forward_anneal_len,
                            inner_forward_candidate.tm_c,
                            inner_forward_candidate.gc_fraction,
                            inner_forward_candidate.anneal_hits,
                        ) + primer_quality(
                            &inner_reverse_full,
                            inner_reverse_candidate
                                .end_0based_exclusive
                                .saturating_sub(inner_reverse_candidate.start_0based),
                            inner_reverse_candidate.tm_c,
                            inner_reverse_candidate.gc_fraction,
                            inner_reverse_candidate.anneal_hits,
                        ) - (overlap_sequence
                            .len()
                            .saturating_sub(constraints.overlap_bp)
                            as f64)
                            * 0.35
                            - ((outer_forward_candidate.tm_c - outer_reverse_candidate.tm_c).abs()
                                * 4.0)
                            - ((inner_forward_candidate.tm_c - inner_reverse_candidate.tm_c).abs()
                                * 4.0);

                        let candidate = OverlapExtensionDesign {
                            score,
                            overlap_sequence: overlap_sequence.to_string(),
                            outer_forward: outer_forward_candidate.clone(),
                            outer_reverse: outer_reverse_candidate.clone(),
                            inner_forward: inner_forward_candidate.clone(),
                            inner_reverse: inner_reverse_candidate.clone(),
                            inner_forward_full: inner_forward_full.clone(),
                            inner_reverse_full: inner_reverse_full.clone(),
                            stage1_left_product,
                            stage1_right_product,
                            final_product,
                        };
                        match &best_design {
                            Some(current) if !better_than(&candidate, current) => {}
                            _ => best_design = Some(candidate),
                        }
                    }
                }
            }
        }

        if evaluation_limited {
            result.warnings.push(format!(
                "Overlap-extension mutagenesis candidate search reached evaluation limit ({max_evaluations}); narrow primer windows or use tighter constraints for a more exhaustive search",
            ));
        }
        let Some(best_design) = best_design else {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "No overlap-extension insertion/deletion mutagenesis primer set satisfied constraints".to_string(),
            });
        };

        parent_seq_ids.push(template.clone());
        let base_prefix = output_prefix
            .as_deref()
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .map(|value| value.to_string())
            .unwrap_or_else(|| format!("{template}_oe_mut"));

        let outer_forward_id = self.materialize_derived_linear_sequence(
            result,
            &format!("{base_prefix}_outer_fwd"),
            &best_design.outer_forward.sequence,
        )?;
        let outer_reverse_id = self.materialize_derived_linear_sequence(
            result,
            &format!("{base_prefix}_outer_rev"),
            &best_design.outer_reverse.sequence,
        )?;
        let inner_forward_id = self.materialize_derived_linear_sequence(
            result,
            &format!("{base_prefix}_inner_fwd"),
            &best_design.inner_forward_full,
        )?;
        let inner_reverse_id = self.materialize_derived_linear_sequence(
            result,
            &format!("{base_prefix}_inner_rev"),
            &best_design.inner_reverse_full,
        )?;
        let stage1_left_id = self.materialize_derived_linear_sequence(
            result,
            &format!("{base_prefix}_stage1_left"),
            &best_design.stage1_left_product,
        )?;
        let stage1_right_id = self.materialize_derived_linear_sequence(
            result,
            &format!("{base_prefix}_stage1_right"),
            &best_design.stage1_right_product,
        )?;
        let mutant_product_id = self.materialize_derived_linear_sequence(
            result,
            &format!("{base_prefix}_mutant"),
            &best_design.final_product,
        )?;

        if let Some(container_id) = self.add_container(
            &[
                outer_forward_id.clone(),
                inner_reverse_id.clone(),
                stage1_left_id.clone(),
            ],
            ContainerKind::Pool,
            Some(format!(
                "OE stage1 left ({})",
                Self::normalize_id_token(&base_prefix)
            )),
            Some(&result.op_id),
        ) {
            result.messages.push(format!(
                "Created overlap-extension stage-1 left container '{}' ({}, {}, {})",
                container_id, outer_forward_id, inner_reverse_id, stage1_left_id
            ));
        }
        if let Some(container_id) = self.add_container(
            &[
                inner_forward_id.clone(),
                outer_reverse_id.clone(),
                stage1_right_id.clone(),
            ],
            ContainerKind::Pool,
            Some(format!(
                "OE stage1 right ({})",
                Self::normalize_id_token(&base_prefix)
            )),
            Some(&result.op_id),
        ) {
            result.messages.push(format!(
                "Created overlap-extension stage-1 right container '{}' ({}, {}, {})",
                container_id, inner_forward_id, outer_reverse_id, stage1_right_id
            ));
        }
        if let Some(container_id) = self.add_container(
            &[
                outer_forward_id.clone(),
                outer_reverse_id.clone(),
                mutant_product_id.clone(),
            ],
            ContainerKind::Pool,
            Some(format!(
                "OE stage2 final ({})",
                Self::normalize_id_token(&base_prefix)
            )),
            Some(&result.op_id),
        ) {
            result.messages.push(format!(
                "Created overlap-extension stage-2 container '{}' ({}, {}, {})",
                container_id, outer_forward_id, outer_reverse_id, mutant_product_id
            ));
        }

        if inserted_bp > 0 {
            // Keep the context lane proportional to the configured overlap so
            // the strip remains readable across short and long overlap designs.
            let flank_bp = constraints.overlap_bp.saturating_mul(4).max(18);
            let overlap_bp = best_design.overlap_sequence.len().max(1);
            let insert_bp = inserted_bp;
            let protocol = ProtocolCartoonKind::PcrOeSubstitution;
            result.protocol_cartoon_preview = Some(ProtocolCartoonPreviewTelemetry {
                protocol: protocol.id().to_string(),
                flank_bp,
                overlap_bp,
                insert_bp,
                bindings: crate::protocol_cartoon::pcr_oe_substitution_geometry_bindings(
                    flank_bp, overlap_bp, insert_bp,
                ),
            });
            result.messages.push(format!(
                "Protocol cartoon preview attached: protocol='{}' flank_bp={} overlap_bp={} insert_bp={}",
                protocol.id(),
                flank_bp,
                overlap_bp,
                insert_bp
            ));
        }

        let mode = if deleted_bp == 0 {
            "insertion"
        } else if inserted_bp == 0 {
            "deletion"
        } else {
            "replacement"
        };
        result.messages.push(format!(
            "Overlap-extension {} mutagenesis from '{}' edit {}..{} (+{} bp / -{} bp) selected overlap={} bp and produced staged artifacts: left='{}', right='{}', final='{}'",
            mode,
            template,
            edit_start_0based,
            edit_end_0based_exclusive,
            inserted_bp,
            deleted_bp,
            best_design.overlap_sequence.len(),
            stage1_left_id,
            stage1_right_id,
            mutant_product_id
        ));
        result.messages.push(format!(
            "Selected primer coordinates: outer_fwd={}..{}, inner_rev={}..{}, inner_fwd={}..{}, outer_rev={}..{}",
            best_design.outer_forward.start_0based,
            best_design.outer_forward.end_0based_exclusive,
            best_design.inner_reverse.start_0based,
            best_design.inner_reverse.end_0based_exclusive,
            best_design.inner_forward.start_0based,
            best_design.inner_forward.end_0based_exclusive,
            best_design.outer_reverse.start_0based,
            best_design.outer_reverse.end_0based_exclusive
        ));
        Ok(())
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
            protocol_cartoon_preview: None,
            genome_annotation_projection: None,
            sequence_alignment: None,
            sequencing_confirmation_report: None,
            sequencing_primer_overlay_report: None,
            sequencing_trace_import_report: None,
            sequencing_trace_record: None,
            sequencing_trace_summaries: None,
            rna_read_gene_support_summary: None,
            rna_read_gene_support_audit: None,
            tfbs_region_summary: None,
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
                        gene_extract_mode: None,
                        promoter_upstream_bp: None,
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
                let svg = export_lineage_svg(&self.state, self.operation_log());
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
                let mut assembled_product_seq_id: Option<String> = None;
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
                    let container_name = self
                        .state
                        .sequences
                        .get(&seq_id)
                        .and_then(|dna| dna.name().clone())
                        .filter(|name| !name.trim().is_empty())
                        .unwrap_or_else(|| seq_id.clone());
                    let _ = self.add_container(
                        std::slice::from_ref(&seq_id),
                        ContainerKind::Singleton,
                        Some(container_name),
                        Some(&result.op_id),
                    );
                    if output.kind == "assembled_product" {
                        assembled_product_seq_id = Some(seq_id.clone());
                    }
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
                if let Some(product_seq_id) = assembled_product_seq_id.as_deref() {
                    match self.maybe_create_gibson_arrangement(&plan, product_seq_id, &result.op_id)
                    {
                        Ok(Some(arrangement_id)) => {
                            let arrangement = self
                                .state
                                .container_state
                                .arrangements
                                .get(&arrangement_id)
                                .cloned();
                            let lane_count = arrangement
                                .as_ref()
                                .map(|arrangement| arrangement.lane_container_ids.len())
                                .unwrap_or(0);
                            let ladders = arrangement
                                .as_ref()
                                .map(|arrangement| {
                                    if arrangement.ladders.is_empty() {
                                        "auto".to_string()
                                    } else {
                                        arrangement.ladders.join(" + ")
                                    }
                                })
                                .unwrap_or_else(|| "auto".to_string());
                            let rack_id = arrangement
                                .and_then(|arrangement| arrangement.default_rack_id)
                                .unwrap_or_else(|| "draft-on-demand".to_string());
                            result.messages.push(format!(
                                "Created Gibson serial arrangement '{}' with {} sample lane(s), ladders '{}', and rack draft '{}'",
                                arrangement_id, lane_count, ladders, rack_id
                            ));
                        }
                        Ok(None) => {}
                        Err(err) => {
                            result.warnings.push(format!(
                                "Gibson arrangement was not created automatically: {}",
                                err.message
                            ));
                        }
                    }
                }
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
                    None,
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
            Operation::SetArrangementLadders {
                arrangement_id,
                ladders,
            } => {
                let normalized = self.set_arrangement_ladders(&arrangement_id, ladders)?;
                let ladder_text = if normalized.is_empty() {
                    "auto".to_string()
                } else {
                    normalized.join(" + ")
                };
                result.messages.push(format!(
                    "Updated arrangement '{}' ladders to {}",
                    arrangement_id.trim(),
                    ladder_text
                ));
            }
            Operation::CreateRackFromArrangement {
                arrangement_id,
                rack_id,
                name,
                profile,
            } => {
                let created_id = self.create_rack_from_arrangement(
                    &arrangement_id,
                    rack_id,
                    name,
                    profile,
                    Some(&result.op_id),
                    true,
                )?;
                result.messages.push(format!(
                    "Created rack '{}' from arrangement '{}'",
                    created_id,
                    arrangement_id.trim()
                ));
            }
            Operation::PlaceArrangementOnRack {
                arrangement_id,
                rack_id,
            } => {
                let start_index = self.place_arrangement_on_rack(&arrangement_id, &rack_id)?;
                result.messages.push(format!(
                    "Placed arrangement '{}' onto rack '{}' starting at slot {}",
                    arrangement_id.trim(),
                    rack_id.trim(),
                    start_index + 1
                ));
            }
            Operation::MoveRackPlacement {
                rack_id,
                from_coordinate,
                to_coordinate,
                move_block,
            } => {
                self.move_rack_placement(&rack_id, &from_coordinate, &to_coordinate, move_block)?;
                result.messages.push(if move_block {
                    format!(
                        "Moved arrangement block on rack '{}' from '{}' to '{}'",
                        rack_id.trim(),
                        from_coordinate.trim(),
                        to_coordinate.trim()
                    )
                } else {
                    format!(
                        "Moved rack placement on '{}' from '{}' to '{}'",
                        rack_id.trim(),
                        from_coordinate.trim(),
                        to_coordinate.trim()
                    )
                });
            }
            Operation::MoveRackSamples {
                rack_id,
                from_coordinates,
                to_coordinate,
            } => {
                self.move_rack_samples(&rack_id, &from_coordinates, &to_coordinate)?;
                result.messages.push(format!(
                    "Moved {} rack sample(s) on '{}' to '{}'",
                    from_coordinates.len(),
                    rack_id.trim(),
                    to_coordinate.trim()
                ));
            }
            Operation::MoveRackArrangementBlocks {
                rack_id,
                arrangement_ids,
                to_coordinate,
            } => {
                self.move_rack_arrangement_blocks(&rack_id, &arrangement_ids, &to_coordinate)?;
                result.messages.push(format!(
                    "Moved {} arrangement block(s) on rack '{}' to '{}'",
                    arrangement_ids.len(),
                    rack_id.trim(),
                    to_coordinate.trim()
                ));
            }
            Operation::SetRackProfile { rack_id, profile } => {
                self.set_rack_profile(&rack_id, profile)?;
                result.messages.push(format!(
                    "Updated rack '{}' profile to '{}'",
                    rack_id.trim(),
                    profile.as_str()
                ));
            }
            Operation::ApplyRackTemplate { rack_id, template } => {
                self.apply_rack_template(&rack_id, template)?;
                result.messages.push(format!(
                    "Applied rack template '{}' to '{}'",
                    template.as_str(),
                    rack_id.trim()
                ));
            }
            Operation::SetRackFillDirection {
                rack_id,
                fill_direction,
            } => {
                self.set_rack_fill_direction(&rack_id, fill_direction)?;
                result.messages.push(format!(
                    "Updated rack '{}' fill direction to '{}'",
                    rack_id.trim(),
                    fill_direction.as_str()
                ));
            }
            Operation::SetRackProfileCustom {
                rack_id,
                rows,
                columns,
            } => {
                self.set_rack_profile_custom(&rack_id, rows, columns)?;
                result.messages.push(format!(
                    "Updated rack '{}' profile to custom {}x{}",
                    rack_id.trim(),
                    rows,
                    columns
                ));
            }
            Operation::SetRackBlockedCoordinates {
                rack_id,
                blocked_coordinates,
            } => {
                self.set_rack_blocked_coordinates(&rack_id, blocked_coordinates.clone())?;
                result.messages.push(if blocked_coordinates.is_empty() {
                    format!("Cleared blocked rack positions on '{}'", rack_id.trim())
                } else {
                    format!(
                        "Updated blocked rack positions on '{}' ({})",
                        rack_id.trim(),
                        blocked_coordinates.join(", ")
                    )
                });
            }
            Operation::ExportRackLabelsSvg {
                rack_id,
                path,
                arrangement_id,
                preset,
            } => {
                let count = self.export_rack_labels_svg(
                    &rack_id,
                    arrangement_id.as_deref(),
                    preset,
                    &path,
                )?;
                result.messages.push(format!(
                    "Wrote {} rack label(s) for '{}' to '{}' using preset '{}'",
                    count,
                    rack_id.trim(),
                    path,
                    preset.as_str()
                ));
            }
            Operation::ExportRackFabricationSvg {
                rack_id,
                path,
                template,
            } => {
                let spec = self.export_rack_fabrication_svg(&rack_id, template, &path)?;
                result.messages.push(format!(
                    "Wrote rack fabrication SVG for '{}' to '{}' using template '{}' ({:.1} x {:.1} mm)",
                    rack_id.trim(),
                    path,
                    template.as_str(),
                    spec.overall_width_mm,
                    spec.overall_depth_mm
                ));
            }
            Operation::ExportRackIsometricSvg {
                rack_id,
                path,
                template,
            } => {
                let spec = self.export_rack_isometric_svg(&rack_id, template, &path)?;
                result.messages.push(format!(
                    "Wrote rack isometric SVG for '{}' to '{}' using template '{}' ({:.1} x {:.1} mm)",
                    rack_id.trim(),
                    path,
                    template.as_str(),
                    spec.overall_width_mm,
                    spec.overall_depth_mm
                ));
            }
            Operation::ExportRackOpenScad {
                rack_id,
                path,
                template,
            } => {
                let spec = self.export_rack_openscad(&rack_id, template, &path)?;
                result.messages.push(format!(
                    "Wrote rack OpenSCAD for '{}' to '{}' using template '{}' ({:.1} x {:.1} mm)",
                    rack_id.trim(),
                    path,
                    template.as_str(),
                    spec.overall_width_mm,
                    spec.overall_depth_mm
                ));
            }
            Operation::ExportRackCarrierLabelsSvg {
                rack_id,
                path,
                arrangement_id,
                template,
                preset,
            } => {
                let (label_count, spec) = self.export_rack_carrier_labels_svg(
                    &rack_id,
                    arrangement_id.as_deref(),
                    template,
                    preset,
                    &path,
                )?;
                result.messages.push(format!(
                    "Wrote {} rack carrier label artifact(s) for '{}' to '{}' using template '{}' and preset '{}' ({:.1} x {:.1} mm)",
                    label_count,
                    rack_id.trim(),
                    path,
                    template.as_str(),
                    preset.as_str(),
                    spec.overall_width_mm,
                    spec.overall_depth_mm
                ));
            }
            Operation::ExportRackSimulationJson {
                rack_id,
                path,
                template,
            } => {
                let spec = self.export_rack_simulation_json(&rack_id, template, &path)?;
                result.messages.push(format!(
                    "Wrote rack simulation JSON for '{}' to '{}' using template '{}' ({:.1} x {:.1} mm)",
                    rack_id.trim(),
                    path,
                    template.as_str(),
                    spec.overall_width_mm,
                    spec.overall_depth_mm
                ));
            }
            Operation::RenderPoolGelSvg {
                inputs,
                path,
                ladders,
                container_ids,
                arrangement_id,
                conditions,
            } => {
                if let Some(arrangement_id) = arrangement_id.as_deref().map(str::trim) {
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
                    if arrangement_id.is_empty() {
                        return Err(EngineError {
                            code: ErrorCode::InvalidInput,
                            message: "arrangement_id cannot be empty".to_string(),
                        });
                    }
                } else if let Some(container_ids) = container_ids.as_ref() {
                    if !inputs.is_empty() {
                        result.warnings.push(
                            "RenderPoolGelSvg ignored 'inputs' because container_ids were provided"
                                .to_string(),
                        );
                    }
                    if container_ids.is_empty() {
                        return Err(EngineError {
                            code: ErrorCode::InvalidInput,
                            message: "container_ids was provided but empty".to_string(),
                        });
                    }
                }
                let layout = self.build_serial_gel_layout_for_render(
                    &inputs,
                    container_ids.as_deref(),
                    arrangement_id.as_deref(),
                    ladders.as_deref(),
                    conditions.as_ref(),
                )?;
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
                    "Wrote serial gel SVG for {} sample lane(s), {} sequence(s) to '{}' (ladders: {}, conditions: {})",
                    layout.sample_count,
                    layout.pool_member_count,
                    path,
                    ladders_used,
                    layout.conditions.describe()
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
                let _ = self.extract_genome_region_into_state(
                    &mut result,
                    &genome_id,
                    &chromosome,
                    start_1based,
                    end_1based,
                    output_id,
                    annotation_scope,
                    max_annotation_features,
                    include_genomic_annotation,
                    catalog_path,
                    cache_dir,
                    "ExtractGenomeRegion",
                )?;
            }
            Operation::FetchDbSnpRegion {
                rs_id,
                genome_id,
                flank_bp,
                output_id,
                annotation_scope,
                max_annotation_features,
                catalog_path,
                cache_dir,
            } => {
                let genome_id = genome_id.trim().to_string();
                if genome_id.is_empty() {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: "Genome id cannot be empty for dbSNP fetch".to_string(),
                    });
                }
                let (display_rs_id, refsnp_id) = Self::normalize_dbsnp_rs_id(&rs_id)?;
                let mut emit_progress = |stage: DbSnpFetchStage, detail: String| {
                    let _ = on_progress(OperationProgress::DbSnpFetch(DbSnpFetchProgress {
                        rs_id: display_rs_id.clone(),
                        genome_id: genome_id.clone(),
                        stage,
                        detail,
                    }));
                };
                emit_progress(
                    DbSnpFetchStage::ValidateInput,
                    format!(
                        "Validated dbSNP request '{}' for prepared genome '{}'",
                        display_rs_id, genome_id
                    ),
                );
                let resolved_catalog_path =
                    catalog_path.unwrap_or_else(|| DEFAULT_GENOME_CATALOG_PATH.to_string());
                let catalog =
                    GenomeCatalog::from_json_file(&resolved_catalog_path).map_err(|e| {
                        EngineError {
                            code: ErrorCode::InvalidInput,
                            message: format!(
                                "Could not open genome catalog '{}': {}",
                                resolved_catalog_path, e
                            ),
                        }
                    })?;
                emit_progress(
                    DbSnpFetchStage::InspectPreparedGenome,
                    format!(
                        "Inspecting prepared genome '{}' via catalog '{}'",
                        genome_id, resolved_catalog_path
                    ),
                );
                let compatibility = catalog
                    .inspect_prepared_genome_compatibility(&genome_id, cache_dir.as_deref())
                    .map_err(|e| EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "Could not inspect prepared compatibility for '{}': {}",
                            genome_id, e
                        ),
                    })?;
                let source_url = Self::dbsnp_refsnp_url(&refsnp_id);
                emit_progress(
                    DbSnpFetchStage::ContactServer,
                    format!(
                        "Contacting NCBI Variation refSNP service at '{}' for '{}'",
                        source_url, display_rs_id
                    ),
                );
                emit_progress(
                    DbSnpFetchStage::WaitResponse,
                    format!("Waiting for dbSNP response from '{}'", source_url),
                );
                let (_source_url_from_fetch, raw_document) =
                    Self::fetch_dbsnp_refsnp_text(&refsnp_id)?;
                emit_progress(
                    DbSnpFetchStage::ParseResponse,
                    format!("Parsing dbSNP response for '{}'", display_rs_id),
                );
                let document = Self::parse_dbsnp_refsnp_json(&refsnp_id, &raw_document)?;
                emit_progress(
                    DbSnpFetchStage::ResolvePlacement,
                    format!(
                        "Resolving genomic placement for '{}' against genome '{}'",
                        display_rs_id, genome_id
                    ),
                );
                let placement = Self::resolve_dbsnp_primary_placement(
                    &document,
                    &display_rs_id,
                    compatibility.requested_family.as_deref(),
                )?;
                let flank_bp = flank_bp.unwrap_or(3000);
                let start_1based = placement.position_1based.saturating_sub(flank_bp).max(1);
                let end_1based = placement.position_1based.saturating_add(flank_bp);
                let chromosome_label = if placement.chromosome_display.trim().is_empty() {
                    placement.chromosome.clone()
                } else {
                    placement.chromosome_display.clone()
                };
                let chromosome_message_label = if placement.chromosome.trim().is_empty()
                    || placement
                        .chromosome
                        .trim()
                        .eq_ignore_ascii_case(chromosome_label.trim())
                {
                    chromosome_label.clone()
                } else {
                    format!("{} [{}]", chromosome_label, placement.chromosome)
                };
                let output_id = output_id.or_else(|| {
                    Some(format!(
                        "{}_{}_{}_{}",
                        Self::normalize_id_token(&display_rs_id),
                        Self::normalize_id_token(&chromosome_label),
                        start_1based,
                        end_1based
                    ))
                });
                let annotation_scope =
                    Some(annotation_scope.unwrap_or(GenomeAnnotationScope::Full));
                let include_genomic_annotation = Some(!matches!(
                    annotation_scope,
                    Some(GenomeAnnotationScope::None)
                ));
                let max_annotation_features = Some(max_annotation_features.unwrap_or(0));
                let assembly_name = placement
                    .assembly_name
                    .as_deref()
                    .map(str::trim)
                    .filter(|value| !value.is_empty())
                    .unwrap_or("unknown");
                let genes = if placement.gene_symbols.is_empty() {
                    String::new()
                } else {
                    format!(" [genes={}]", placement.gene_symbols.join(", "))
                };
                emit_progress(
                    DbSnpFetchStage::ExtractRegion,
                    format!(
                        "Extracting annotated slice {}:{}-{} from '{}'",
                        chromosome_message_label, start_1based, end_1based, genome_id
                    ),
                );
                result.messages.push(format!(
                    "Resolved dbSNP '{}' from '{}' to {}:{} (assembly={}) and extracted +/-{} bp on '{}'{}",
                    placement.rs_id,
                    source_url,
                    chromosome_message_label,
                    placement.position_1based,
                    assembly_name,
                    flank_bp,
                    genome_id,
                    genes
                ));
                let seq_id = self.extract_genome_region_into_state(
                    &mut result,
                    &genome_id,
                    &placement.chromosome,
                    start_1based,
                    end_1based,
                    output_id,
                    annotation_scope,
                    max_annotation_features,
                    include_genomic_annotation,
                    Some(resolved_catalog_path),
                    cache_dir,
                    "FetchDbSnpRegion",
                )?;
                let local_start_0based = placement.position_1based.saturating_sub(start_1based);
                if let Some(dna) = self.state.sequences.get_mut(&seq_id) {
                    if local_start_0based < dna.len() {
                        emit_progress(
                            DbSnpFetchStage::AttachVariantMarker,
                            format!(
                                "Attaching local rs marker '{}' at sequence position {}",
                                display_rs_id,
                                local_start_0based.saturating_add(1)
                            ),
                        );
                        dna.features_mut()
                            .push(Self::build_dbsnp_variant_marker_feature(
                                &display_rs_id,
                                &placement.chromosome,
                                &placement.chromosome_display,
                                placement.position_1based,
                                placement.assembly_name.as_deref(),
                                &placement.gene_symbols,
                                local_start_0based,
                            ));
                        Self::prepare_sequence(dna);
                    } else {
                        result.warnings.push(format!(
                            "Resolved dbSNP '{}' fell outside extracted sequence '{}' (local position {}, length {})",
                            display_rs_id,
                            seq_id,
                            local_start_0based + 1,
                            dna.len()
                        ));
                    }
                }
            }
            Operation::ExtractGenomeGene {
                genome_id,
                gene_query,
                occurrence,
                output_id,
                extract_mode,
                promoter_upstream_bp,
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
                let extract_mode = extract_mode.unwrap_or(GenomeGeneExtractMode::Gene);
                let promoter_upstream_bp = promoter_upstream_bp.unwrap_or(0);
                if matches!(extract_mode, GenomeGeneExtractMode::Gene) && promoter_upstream_bp > 0 {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: "promoter_upstream_bp requires extract_mode=coding_with_promoter"
                            .to_string(),
                    });
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
                            "Could not inspect transcript/exon annotation for gene '{}': {}",
                            Self::genome_gene_display_label(selected_gene),
                            e
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
                                    Self::genome_gene_display_label(selected_gene),
                                    records.len()
                                ));
                                transcript_records = records;
                            }
                        }
                        Err(e) => {
                            result.warnings.push(format!(
                                "Could not run transcript fallback for gene '{}': {}",
                                Self::genome_gene_display_label(selected_gene),
                                e
                            ));
                        }
                    }
                }

                let (extract_start_1based, extract_end_1based) =
                    Self::resolve_extract_genome_gene_interval(
                        selected_gene,
                        &transcript_records,
                        extract_mode,
                        promoter_upstream_bp,
                    )
                    .map_err(|message| EngineError {
                        code: ErrorCode::NotFound,
                        message,
                    })?;
                let sequence = catalog
                    .get_sequence_region_with_cache(
                        &genome_id,
                        &selected_gene.chromosome,
                        extract_start_1based,
                        extract_end_1based,
                        cache_dir.as_deref(),
                    )
                    .map_err(|e| EngineError {
                        code: ErrorCode::NotFound,
                        message: format!(
                            "Could not load gene extraction interval {}:{}-{} from '{}': {}",
                            selected_gene.chromosome,
                            extract_start_1based,
                            extract_end_1based,
                            genome_id,
                            e
                        ),
                    })?;
                let default_id = Self::default_extract_genome_gene_output_id(
                    &genome_id,
                    selected_gene,
                    extract_mode,
                    promoter_upstream_bp,
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

                if let Some(extracted_sequence) = self
                    .state
                    .sequences
                    .get(&seq_id)
                    .map(|dna| dna.get_forward_string())
                    && let Some(exon_projection) = Self::build_exon_concatenated_projection(
                        &extracted_sequence,
                        extract_start_1based,
                        extract_end_1based,
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
                                Some(extract_start_1based.to_string()),
                            ),
                            (
                                "genomic_end_1based".into(),
                                Some(extract_end_1based.to_string()),
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
                            kind: "mRNA".into(),
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
                            kind: "exon".into(),
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
                    extract_start_1based,
                    extract_end_1based,
                    requested_scope,
                );
                let candidate_before_fallback = projection.feature_count();
                if let Some(cap) = max_annotation_features {
                    if projection.feature_count() > cap {
                        if matches!(requested_scope, GenomeAnnotationScope::Full) {
                            let core_projection = Self::build_extract_region_annotation_projection(
                                &gene_records,
                                &transcript_records,
                                extract_start_1based,
                                extract_end_1based,
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
                    start_1based: Some(extract_start_1based),
                    end_1based: Some(extract_end_1based),
                    gene_query: Some(query.to_string()),
                    occurrence: Some(occurrence),
                    gene_extract_mode: Some(extract_mode.as_str().to_string()),
                    promoter_upstream_bp: matches!(
                        extract_mode,
                        GenomeGeneExtractMode::CodingWithPromoter
                    )
                    .then_some(promoter_upstream_bp),
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
                    "Extracted genome gene '{}' [{} match, occurrence {}] as '{}' from '{}' ({}, extract_mode={}, promoter_upstream_bp={})",
                    query,
                    match_mode,
                    occurrence,
                    seq_id,
                    genome_id,
                    Self::genome_gene_display_label(selected_gene),
                    extract_mode.as_str(),
                    promoter_upstream_bp
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
                    gene_extract_mode: None,
                    promoter_upstream_bp: None,
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
                    gene_extract_mode: None,
                    promoter_upstream_bp: None,
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
                self.execute_design_primer_pairs(
                    &mut result,
                    &mut parent_seq_ids,
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
                    None,
                )?;
            }
            Operation::DesignInsertionPrimerPairs {
                template,
                insertion,
                mut forward,
                mut reverse,
                pair_constraints,
                min_amplicon_bp,
                max_amplicon_bp,
                max_tm_delta_c,
                max_pairs,
                report_id,
            } => {
                let template_len = self
                    .state
                    .sequences
                    .get(&template)
                    .ok_or_else(|| EngineError {
                        code: ErrorCode::NotFound,
                        message: format!("Sequence '{template}' not found"),
                    })?
                    .get_forward_string()
                    .len();
                let insertion = Self::normalize_primer_insertion_intent(&insertion, template_len)?;
                let (roi_start_0based, roi_end_0based) =
                    Self::derive_insertion_intent_roi_bounds(&insertion, template_len);

                forward.start_0based = Some(insertion.forward_window_start_0based);
                forward.end_0based = Some(insertion.forward_window_end_0based_exclusive);
                forward.non_annealing_5prime_tail =
                    Some(insertion.forward_extension_5prime.clone());
                reverse.start_0based = Some(insertion.reverse_window_start_0based);
                reverse.end_0based = Some(insertion.reverse_window_end_0based_exclusive);
                reverse.non_annealing_5prime_tail =
                    Some(insertion.reverse_extension_5prime.clone());

                self.execute_design_primer_pairs(
                    &mut result,
                    &mut parent_seq_ids,
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
                    Some(insertion),
                )?;
            }
            Operation::PcrOverlapExtensionMutagenesis {
                template,
                edit_start_0based,
                edit_end_0based_exclusive,
                insert_sequence,
                constraints,
                output_prefix,
            } => {
                self.execute_overlap_extension_mutagenesis(
                    &mut result,
                    &mut parent_seq_ids,
                    template,
                    edit_start_0based,
                    edit_end_0based_exclusive,
                    insert_sequence,
                    constraints,
                    output_prefix,
                )?;
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
                    let (
                        derived_dna,
                        transcript_id,
                        transcript_label,
                        is_reverse,
                        exon_count,
                        protein_derivation,
                    ) = Self::derive_transcript_sequence_from_feature(
                        &source_sequence_upper,
                        source_feature,
                        &source_features,
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
                    if let Some(derivation) = protein_derivation {
                        if derivation.protein_length_aa > 0 {
                            result.messages.push(format!(
                                "Derived protein for transcript '{}' using translation table {} ('{}', source={}, {} aa).",
                                derived_seq_id,
                                derivation.translation_table,
                                derivation.translation_table_label,
                                derivation.translation_table_source.as_str(),
                                derivation.protein_length_aa
                            ));
                        } else {
                            result.warnings.push(format!(
                                "Transcript '{}' resolved CDS context but did not yield a non-empty protein sequence.",
                                derived_seq_id
                            ));
                        }
                        for warning in derivation.warnings {
                            result.warnings.push(format!(
                                "Transcript '{}': {}",
                                derivation.transcript_id, warning
                            ));
                        }
                    } else {
                        result.warnings.push(format!(
                            "Transcript '{}' has no CDS annotation; protein translation was not derived.",
                            derived_seq_id
                        ));
                    }
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
            Operation::DeriveProteinSequences {
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
                    .unwrap_or_else(|| format!("{seq_id}__protein"));
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
                    let (
                        derived_transcript,
                        transcript_id,
                        transcript_label,
                        _is_reverse,
                        _exon_count,
                        annotated_derivation,
                    ) = Self::derive_transcript_sequence_from_feature(
                        &source_sequence_upper,
                        source_feature,
                        &source_features,
                        transcript_feature_id,
                        &seq_id,
                    )?;
                    let representative_cds_feature = Self::collect_matching_cds_features_for_derivation(
                        &source_features,
                        source_feature,
                        &Self::feature_ranges_0based_for_derivation(source_feature),
                    )
                    .into_iter()
                    .next();
                    let derivation = match annotated_derivation {
                        Some(derivation) => Some(derivation),
                        None => Self::infer_transcript_protein_derivation_without_annotation(
                            &derived_transcript.get_forward_string(),
                            source_feature,
                            transcript_feature_id,
                            &seq_id,
                            &source_features,
                            &transcript_id,
                            &transcript_label,
                        )?,
                    };
                    let Some(derivation) = derivation else {
                        result.warnings.push(format!(
                            "Transcript feature n-{} in '{}' could not yield a protein sequence.",
                            transcript_feature_id + 1,
                            seq_id
                        ));
                        continue;
                    };
                    if derivation.protein_length_aa == 0 || derivation.protein_sequence.is_empty() {
                        result.warnings.push(format!(
                            "Transcript '{}' did not yield a non-empty protein sequence.",
                            derivation.transcript_id
                        ));
                        continue;
                    }
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
                    let protein_seq_id = self.unique_seq_id(&base_seq_id);
                    let protein_name = format!("{transcript_label} protein");
                    let protein = Self::build_transcript_derived_protein_sequence(
                        &derivation.protein_sequence,
                        &protein_name,
                        &seq_id,
                        transcript_feature_id,
                        source_feature,
                        representative_cds_feature,
                        &derivation,
                    )?;
                    self.state
                        .sequences
                        .insert(protein_seq_id.clone(), protein);
                    self.add_lineage_node(
                        &protein_seq_id,
                        SequenceOrigin::Derived,
                        Some(&result.op_id),
                    );
                    result.created_seq_ids.push(protein_seq_id.clone());
                    result.messages.push(format!(
                        "Derived protein '{}' from transcript '{}' using translation table {} ('{}', mode={}, {} aa).",
                        protein_seq_id,
                        derivation.transcript_id,
                        derivation.translation_table,
                        derivation.translation_table_label,
                        derivation.derivation_mode.as_str(),
                        derivation.protein_length_aa
                    ));
                    for warning in derivation.warnings {
                        result.warnings.push(format!(
                            "Protein '{}' (transcript '{}'): {}",
                            protein_seq_id, derivation.transcript_id, warning
                        ));
                    }
                }
                if result.created_seq_ids.is_empty() {
                    return Err(EngineError {
                        code: ErrorCode::NotFound,
                        message: format!(
                            "DeriveProteinSequences did not produce proteins for '{}'",
                            seq_id
                        ),
                    });
                }
            }
            Operation::ReverseTranslateProteinSequence {
                seq_id,
                output_id,
                speed_profile,
                speed_mark,
                translation_table,
                target_anneal_tm_c,
                anneal_window_bp,
            } => {
                let protein = self
                    .state
                    .sequences
                    .get(&seq_id)
                    .ok_or_else(|| EngineError {
                        code: ErrorCode::NotFound,
                        message: format!("Sequence '{seq_id}' not found"),
                    })?
                    .clone();
                if !protein.is_protein_sequence() {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "ReverseTranslateProteinSequence requires a protein sequence; '{}' has molecule_type={:?}",
                            seq_id,
                            protein.molecule_type()
                        ),
                    });
                }
                let effective_speed_profile = speed_profile
                    .or_else(|| Self::sequence_feature_translation_speed_hint(&protein));
                let (preferred_species_label, mut warnings) = if let Some(profile) = effective_speed_profile
                {
                    let (species, warning) = Self::codon_profile_species_label(profile);
                    let mut warnings = vec![];
                    if let Some(warning) = warning {
                        warnings.push(warning.to_string());
                    }
                    (Some(species), warnings)
                } else {
                    (None, vec![])
                };
                let effective_translation_table = translation_table
                    .or_else(|| {
                        protein.features().iter().find_map(|feature| {
                            Self::feature_qualifier_text(feature, "translation_table")
                                .or_else(|| Self::feature_qualifier_text(feature, "transl_table"))
                                .and_then(|raw| raw.parse::<usize>().ok())
                        })
                    })
                    .unwrap_or(1);
                let anneal_window_bp = anneal_window_bp.unwrap_or(20).max(6);
                let protein_sequence = protein.get_forward_string().to_ascii_uppercase();
                let (coding_sequence, reverse_translation_warnings) =
                    Self::build_reverse_translated_coding_sequence(
                        &protein_sequence,
                        effective_translation_table,
                        preferred_species_label,
                        speed_mark,
                        target_anneal_tm_c,
                        anneal_window_bp,
                    );
                warnings.extend(reverse_translation_warnings);
                let base_seq_id = output_id
                    .as_deref()
                    .map(str::trim)
                    .filter(|value| !value.is_empty())
                    .map(|value| value.to_string())
                    .unwrap_or_else(|| format!("{seq_id}__coding"));
                let coding_seq_id = self.unique_seq_id(&base_seq_id);
                let coding_name = format!(
                    "{} coding sequence",
                    protein.name().as_deref().unwrap_or(seq_id.as_str())
                );
                let coding = Self::build_reverse_translated_coding_dna(
                    &coding_sequence,
                    &coding_name,
                    &seq_id,
                    &protein,
                    effective_translation_table,
                    &Self::translation_table_display_name(effective_translation_table),
                    effective_speed_profile,
                    preferred_species_label,
                    speed_mark,
                    target_anneal_tm_c,
                    anneal_window_bp,
                    &warnings,
                )?;
                self.state.sequences.insert(coding_seq_id.clone(), coding);
                self.add_lineage_node(
                    &coding_seq_id,
                    SequenceOrigin::Derived,
                    Some(&result.op_id),
                );
                parent_seq_ids.push(seq_id.clone());
                result.created_seq_ids.push(coding_seq_id.clone());
                result.messages.push(format!(
                    "Reverse-translated protein '{}' into coding sequence '{}' using translation table {} ('{}'){}{}.",
                    seq_id,
                    coding_seq_id,
                    effective_translation_table,
                    Self::translation_table_display_name(effective_translation_table),
                    effective_speed_profile
                        .map(|profile| format!(", speed_profile={}", profile.as_str()))
                        .unwrap_or_default(),
                    speed_mark
                        .map(|mark| format!(", speed_mark={}", mark.as_str()))
                        .unwrap_or_default()
                ));
                if let Some(target_tm) = target_anneal_tm_c {
                    result.messages.push(format!(
                        "Applied local reverse-translation Tm heuristic target {:.1} °C over {} bp windows.",
                        target_tm, anneal_window_bp
                    ));
                }
                result.warnings.extend(warnings);
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
                let primary_series = Self::build_dotplot_query_series(
                    format!("{dotplot_id}_series_1"),
                    seq_id.clone(),
                    seq_id.clone(),
                    Self::default_dotplot_series_color(0),
                    mode,
                    span_start_0based,
                    span_end_0based,
                    points,
                    boxplot_bins,
                );
                let view = DotplotView {
                    schema: DOTPLOT_VIEW_SCHEMA.to_string(),
                    dotplot_id: dotplot_id.clone(),
                    owner_seq_id: seq_id.clone(),
                    seq_id: seq_id.clone(),
                    reference_seq_id: reference_seq_id_for_view.clone(),
                    generated_at_unix_ms: Self::now_unix_ms(),
                    span_start_0based: primary_series.span_start_0based,
                    span_end_0based: primary_series.span_end_0based,
                    reference_span_start_0based,
                    reference_span_end_0based,
                    mode,
                    word_size,
                    step_bp,
                    max_mismatches,
                    tile_bp,
                    point_count: primary_series.point_count,
                    points: primary_series.points.clone(),
                    boxplot_bin_count: primary_series.boxplot_bin_count,
                    boxplot_bins: primary_series.boxplot_bins.clone(),
                    series_count: 1,
                    query_series: vec![primary_series],
                    reference_annotation: reference_seq_id_for_view.as_deref().and_then(|ref_id| {
                        self.build_dotplot_reference_annotation_track(
                            ref_id,
                            reference_span_start_0based,
                            reference_span_end_0based,
                        )
                    }),
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
            Operation::ComputeDotplotOverlay {
                owner_seq_id,
                reference_seq_id,
                reference_span_start_0based,
                reference_span_end_0based,
                queries,
                word_size,
                step_bp,
                max_mismatches,
                tile_bp,
                store_as,
            } => {
                let owner_seq_id = owner_seq_id.trim().to_string();
                if owner_seq_id.is_empty() {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: "ComputeDotplotOverlay requires owner_seq_id".to_string(),
                    });
                }
                if !self.state.sequences.contains_key(owner_seq_id.as_str()) {
                    return Err(EngineError {
                        code: ErrorCode::NotFound,
                        message: format!("Owner sequence '{}' not found", owner_seq_id),
                    });
                }
                let reference_seq_id = reference_seq_id.trim().to_string();
                if reference_seq_id.is_empty() {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: "ComputeDotplotOverlay requires reference_seq_id".to_string(),
                    });
                }
                let reference_dna =
                    self.state
                        .sequences
                        .get(&reference_seq_id)
                        .ok_or_else(|| EngineError {
                            code: ErrorCode::NotFound,
                            message: format!("Reference sequence '{}' not found", reference_seq_id),
                        })?;
                if queries.is_empty() {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: "ComputeDotplotOverlay requires at least one query spec"
                            .to_string(),
                    });
                }
                let reference_text = reference_dna.get_forward_string().to_ascii_uppercase();
                let reference_bytes = reference_text.as_bytes();
                let (reference_span_start_0based, reference_span_end_0based) =
                    Self::resolve_analysis_span(
                        reference_bytes.len(),
                        reference_span_start_0based,
                        reference_span_end_0based,
                    )?;
                let reference_span =
                    &reference_bytes[reference_span_start_0based..reference_span_end_0based];
                let dotplot_id = if let Some(raw_id) = store_as.as_deref() {
                    Self::normalize_analysis_id(raw_id, "dotplot")?
                } else {
                    format!("dotplot_{}", result.op_id)
                };
                let reference_annotation = self.build_dotplot_reference_annotation_track(
                    &reference_seq_id,
                    reference_span_start_0based,
                    reference_span_end_0based,
                );
                let mut query_series: Vec<DotplotQuerySeries> = vec![];
                let mut truncated_series_labels: Vec<String> = vec![];
                for (index, query) in queries.iter().enumerate() {
                    let query_seq_id = query.seq_id.trim();
                    if query_seq_id.is_empty() {
                        return Err(EngineError {
                            code: ErrorCode::InvalidInput,
                            message: format!(
                                "ComputeDotplotOverlay query spec #{} requires seq_id",
                                index + 1
                            ),
                        });
                    }
                    if !matches!(
                        query.mode,
                        DotplotMode::PairForward | DotplotMode::PairReverseComplement
                    ) {
                        return Err(EngineError {
                            code: ErrorCode::InvalidInput,
                            message: format!(
                                "ComputeDotplotOverlay query '{}' requires pair mode, got '{}'",
                                query_seq_id,
                                query.mode.as_str()
                            ),
                        });
                    }
                    let query_dna =
                        self.state
                            .sequences
                            .get(query_seq_id)
                            .ok_or_else(|| EngineError {
                                code: ErrorCode::NotFound,
                                message: format!("Query sequence '{}' not found", query_seq_id),
                            })?;
                    let query_text = query_dna.get_forward_string().to_ascii_uppercase();
                    let query_bytes = query_text.as_bytes();
                    let (query_span_start_0based, query_span_end_0based) =
                        Self::resolve_analysis_span(
                            query_bytes.len(),
                            query.span_start_0based,
                            query.span_end_0based,
                        )?;
                    let query_span = &query_bytes[query_span_start_0based..query_span_end_0based];
                    let (points, truncated) = Self::compute_dotplot_points(
                        query_span,
                        reference_span,
                        query_span_start_0based,
                        reference_span_start_0based,
                        query.mode,
                        word_size,
                        step_bp,
                        max_mismatches,
                        MAX_DOTPLOT_POINTS,
                    )?;
                    let boxplot_bins = Self::compute_dotplot_boxplot_bins(
                        &points,
                        query_span_start_0based,
                        query_span_end_0based,
                        DOTPLOT_BOXPLOT_DEFAULT_BINS,
                    );
                    let label = if query.label.trim().is_empty() {
                        query_seq_id.to_string()
                    } else {
                        query.label.trim().to_string()
                    };
                    if truncated {
                        truncated_series_labels.push(label.clone());
                    }
                    query_series.push(Self::build_dotplot_query_series(
                        format!("{dotplot_id}_series_{}", index + 1),
                        query_seq_id.to_string(),
                        label,
                        query
                            .color_rgb
                            .unwrap_or_else(|| Self::default_dotplot_series_color(index)),
                        query.mode,
                        query_span_start_0based,
                        query_span_end_0based,
                        points,
                        boxplot_bins,
                    ));
                }
                let primary_series = query_series.first().cloned().ok_or_else(|| EngineError {
                    code: ErrorCode::Internal,
                    message: "ComputeDotplotOverlay did not produce a primary query series"
                        .to_string(),
                })?;
                let view = DotplotView {
                    schema: DOTPLOT_VIEW_SCHEMA.to_string(),
                    dotplot_id: dotplot_id.clone(),
                    owner_seq_id: owner_seq_id.clone(),
                    seq_id: primary_series.seq_id.clone(),
                    reference_seq_id: Some(reference_seq_id.clone()),
                    generated_at_unix_ms: Self::now_unix_ms(),
                    span_start_0based: primary_series.span_start_0based,
                    span_end_0based: primary_series.span_end_0based,
                    reference_span_start_0based,
                    reference_span_end_0based,
                    mode: primary_series.mode,
                    word_size,
                    step_bp,
                    max_mismatches,
                    tile_bp,
                    point_count: primary_series.point_count,
                    points: primary_series.points.clone(),
                    boxplot_bin_count: primary_series.boxplot_bin_count,
                    boxplot_bins: primary_series.boxplot_bins.clone(),
                    series_count: query_series.len(),
                    query_series,
                    reference_annotation,
                };
                let replaced = self
                    .read_dotplot_analysis_store()
                    .dotplots
                    .contains_key(dotplot_id.as_str());
                self.upsert_dotplot_view(view.clone())?;
                result.messages.push(format!(
                    "{} overlay dotplot '{}' for owner '{}' against '{}' (series={}, reference_span={}..{}, word={}, step={}, mismatches={})",
                    if replaced { "Updated" } else { "Created" },
                    dotplot_id,
                    owner_seq_id,
                    reference_seq_id,
                    view.series_count,
                    reference_span_start_0based,
                    reference_span_end_0based,
                    word_size,
                    step_bp,
                    max_mismatches
                ));
                if !truncated_series_labels.is_empty() {
                    result.warnings.push(format!(
                        "Overlay dotplot '{}' truncated {} series at MAX_DOTPLOT_POINTS ({MAX_DOTPLOT_POINTS}): {}",
                        dotplot_id,
                        truncated_series_labels.len(),
                        truncated_series_labels.join(", ")
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
                let query_text = query_dna.get_forward_string();
                let target_text = target_dna.get_forward_string();
                let report = Self::compute_pairwise_alignment_report(
                    &query_seq_id,
                    query_text.as_str(),
                    query_span_start_0based,
                    query_span_end_0based,
                    &target_seq_id,
                    target_text.as_str(),
                    target_span_start_0based,
                    target_span_end_0based,
                    mode,
                    match_score,
                    mismatch_score,
                    gap_open,
                    gap_extend,
                )?
                .report;
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
            Operation::ConfirmConstructReads {
                expected_seq_id,
                baseline_seq_id,
                read_seq_ids,
                trace_ids,
                targets,
                alignment_mode,
                match_score,
                mismatch_score,
                gap_open,
                gap_extend,
                min_identity_fraction,
                min_target_coverage_fraction,
                allow_reverse_complement,
                report_id,
            } => {
                let report = self.confirm_construct_reads(
                    &expected_seq_id,
                    baseline_seq_id.as_deref(),
                    &read_seq_ids,
                    &trace_ids,
                    &targets,
                    alignment_mode,
                    match_score,
                    mismatch_score,
                    gap_open,
                    gap_extend,
                    min_identity_fraction,
                    min_target_coverage_fraction,
                    allow_reverse_complement,
                    report_id.as_deref(),
                )?;
                result.sequencing_confirmation_report = Some(report.clone());
                result.messages.push(format!(
                    "Sequencing confirmation stored report '{}': {}",
                    report.report_id,
                    Self::format_sequencing_confirmation_report_detail_summary(&report)
                ));
                for target in report.targets.iter().take(8) {
                    result.messages.push(format!(
                        "  - {} [{}] {}",
                        target.label,
                        target.target_id,
                        target.status.as_str()
                    ));
                }
                if report.targets.len() > 8 {
                    result.messages.push(format!(
                        "  ... {} additional target row(s) omitted",
                        report.targets.len() - 8
                    ));
                }
            }
            Operation::SuggestSequencingPrimers {
                expected_seq_id,
                primer_seq_ids,
                confirmation_report_id,
                min_3prime_anneal_bp,
                predicted_read_length_bp,
            } => {
                let report = self.suggest_sequencing_primers(
                    &expected_seq_id,
                    &primer_seq_ids,
                    confirmation_report_id.as_deref(),
                    min_3prime_anneal_bp,
                    predicted_read_length_bp,
                )?;
                result.sequencing_primer_overlay_report = Some(report.clone());
                result.messages.push(format!(
                    "Sequencing-primer overlay report for '{}' (primers={}, suggestions={}, guidance_rows={}, proposal_rows={}, min_3prime_anneal_bp={}, predicted_read_length_bp={})",
                    report.expected_seq_id,
                    report.primer_seq_ids.len(),
                    report.suggestion_count,
                    report.problem_guidance_count,
                    report.proposal_count,
                    report.min_3prime_anneal_bp,
                    report.predicted_read_length_bp
                ));
                for row in report.suggestions.iter().take(8) {
                    result.messages.push(format!(
                        "  - {} [{}] {} {}..{} read_span={}..{} flagged_targets={} flagged_variants={}",
                        row.primer_label,
                        row.primer_seq_id,
                        row.orientation.as_str(),
                        row.anneal_start_0based,
                        row.anneal_end_0based_exclusive,
                        row.predicted_read_span_start_0based,
                        row.predicted_read_span_end_0based_exclusive,
                        row.covered_problem_target_ids.len(),
                        row.covered_problem_variant_ids.len()
                    ));
                }
                if report.suggestions.len() > 8 {
                    result.messages.push(format!(
                        "  ... {} additional sequencing-primer suggestion row(s) omitted",
                        report.suggestions.len() - 8
                    ));
                }
            }
            Operation::ImportSequencingTrace {
                path,
                trace_id,
                seq_id,
            } => {
                let report =
                    self.import_sequencing_trace(&path, trace_id.as_deref(), seq_id.as_deref())?;
                let record = self.get_sequencing_trace(&report.trace_id)?;
                result.sequencing_trace_import_report = Some(report.clone());
                result.sequencing_trace_record = Some(record.clone());
                result.messages.push(format!(
                    "Imported sequencing trace '{}': {}",
                    report.trace_id,
                    Self::format_sequencing_trace_detail_summary(&record)
                ));
                if !report.warnings.is_empty() {
                    for warning in report.warnings.iter().take(8) {
                        result.warnings.push(warning.clone());
                    }
                    if report.warnings.len() > 8 {
                        result.warnings.push(format!(
                            "... {} additional sequencing-trace warning(s) omitted",
                            report.warnings.len() - 8
                        ));
                    }
                }
            }
            Operation::ListSequencingTraces { seq_id } => {
                let rows = self.list_sequencing_traces(seq_id.as_deref());
                result.sequencing_trace_summaries = Some(rows.clone());
                result.messages.push(format!(
                    "Sequencing traces: {} row(s){}",
                    rows.len(),
                    seq_id
                        .as_deref()
                        .map(|s| format!(" (seq_id='{}')", s))
                        .unwrap_or_default()
                ));
                for row in rows.iter().take(8) {
                    result.messages.push(format!(
                        "  - {}",
                        Self::format_sequencing_trace_summary_row(row)
                    ));
                }
                if rows.len() > 8 {
                    result.messages.push(format!(
                        "  ... {} additional trace row(s) omitted",
                        rows.len() - 8
                    ));
                }
            }
            Operation::ShowSequencingTrace { trace_id } => {
                let record = self.get_sequencing_trace(&trace_id)?;
                result.sequencing_trace_record = Some(record.clone());
                result.messages.push(format!(
                    "Sequencing trace summary: {}",
                    Self::format_sequencing_trace_detail_summary(&record)
                ));
            }
            Operation::ListSequencingConfirmationReports { expected_seq_id } => {
                let rows = self.list_sequencing_confirmation_reports(expected_seq_id.as_deref());
                result.messages.push(format!(
                    "Sequencing-confirmation reports: {} row(s){}",
                    rows.len(),
                    expected_seq_id
                        .as_deref()
                        .map(|s| format!(" (expected_seq_id='{}')", s))
                        .unwrap_or_default()
                ));
                for row in rows.iter().take(8) {
                    result.messages.push(format!(
                        "  - {}",
                        Self::format_sequencing_confirmation_report_summary_row(row)
                    ));
                }
                if rows.len() > 8 {
                    result.messages.push(format!(
                        "  ... {} additional report row(s) omitted",
                        rows.len() - 8
                    ));
                }
            }
            Operation::ShowSequencingConfirmationReport { report_id } => {
                let report = self.get_sequencing_confirmation_report(&report_id)?;
                result.messages.push(format!(
                    "Sequencing-confirmation report summary: {}",
                    Self::format_sequencing_confirmation_report_detail_summary(&report)
                ));
                for target in report.targets.iter().take(8) {
                    result.messages.push(format!(
                        "  {} [{}] status={} coverage={}/{} support={} contradict={}",
                        target.label,
                        target.target_id,
                        target.status.as_str(),
                        target.covered_bp,
                        target.target_length_bp,
                        target.support_read_ids.len(),
                        target.contradicting_read_ids.len()
                    ));
                }
                if report.targets.len() > 8 {
                    result.messages.push(format!(
                        "  ... {} additional target row(s) omitted",
                        report.targets.len() - 8
                    ));
                }
            }
            Operation::ExportSequencingConfirmationReport { report_id, path } => {
                let report = self.export_sequencing_confirmation_report(&report_id, &path)?;
                result.messages.push(format!(
                    "Exported sequencing-confirmation report '{}' to '{}'",
                    report.report_id, path
                ));
            }
            Operation::ExportSequencingConfirmationSupportTsv { report_id, path } => {
                let report = self.export_sequencing_confirmation_support_tsv(&report_id, &path)?;
                result.messages.push(format!(
                    "Exported sequencing-confirmation support TSV '{}' to '{}'",
                    report.report_id, path
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
            Operation::SummarizeRnaReadGeneSupport {
                report_id,
                gene_ids,
                selected_record_indices,
                complete_rule,
                path,
            } => {
                let summary = self.summarize_rna_read_gene_support(
                    &report_id,
                    &gene_ids,
                    &selected_record_indices,
                    complete_rule,
                )?;
                if let Some(path) = path.as_deref() {
                    self.write_rna_read_gene_support_summary_json(&summary, path)?;
                    result.messages.push(format!(
                        "Wrote RNA-read gene-support summary '{}' to '{}'",
                        summary.report_id, path
                    ));
                }
                result.messages.push(format!(
                    "RNA-read gene-support summary for '{}' matched {} gene(s), missing {}, accepted_target_reads={}, complete_rule={}",
                    summary.report_id,
                    summary.matched_gene_ids.len(),
                    summary.missing_gene_ids.len(),
                    summary.accepted_target_count,
                    summary.complete_rule.as_str()
                ));
                result.rna_read_gene_support_summary = Some(summary);
            }
            Operation::InspectRnaReadGeneSupport {
                report_id,
                gene_ids,
                selected_record_indices,
                complete_rule,
                cohort_filter,
                path,
            } => {
                let audit = self.inspect_rna_read_gene_support(
                    &report_id,
                    &gene_ids,
                    &selected_record_indices,
                    complete_rule,
                    cohort_filter,
                )?;
                if let Some(path) = path.as_deref() {
                    self.write_rna_read_gene_support_audit_json(&audit, path)?;
                    result.messages.push(format!(
                        "Wrote RNA-read gene-support audit '{}' to '{}'",
                        audit.report_id, path
                    ));
                }
                result.messages.push(format!(
                    "RNA-read gene-support audit for '{}' matched {} gene(s), missing {}, rows={}, cohort={}",
                    audit.report_id,
                    audit.matched_gene_ids.len(),
                    audit.missing_gene_ids.len(),
                    audit.row_count,
                    audit.cohort_filter.as_str()
                ));
                result.rna_read_gene_support_audit = Some(audit);
            }
            Operation::SummarizeTfbsRegion {
                seq_id,
                focus_start_0based,
                focus_end_0based_exclusive,
                context_start_0based,
                context_end_0based_exclusive,
                min_focus_occurrences,
                min_context_occurrences,
                limit,
                path,
            } => {
                let summary = self.summarize_tfbs_region(TfbsRegionSummaryRequest {
                    seq_id: seq_id.clone(),
                    focus_start_0based,
                    focus_end_0based_exclusive,
                    context_start_0based,
                    context_end_0based_exclusive,
                    min_focus_occurrences,
                    min_context_occurrences,
                    limit,
                })?;
                if let Some(path) = path.as_deref() {
                    self.write_tfbs_region_summary_json(&summary, path)?;
                    result.messages.push(format!(
                        "Wrote TFBS region summary for '{}' to '{}'",
                        summary.seq_id, path
                    ));
                }
                result.messages.push(format!(
                    "TFBS region summary for '{}' matched {} factor(s) in focus {}..{} against context {}..{}",
                    summary.seq_id,
                    summary.matched_tf_count,
                    summary.focus_start_0based,
                    summary.focus_end_0based_exclusive,
                    summary.context_start_0based,
                    summary.context_end_0based_exclusive,
                ));
                result.tfbs_region_summary = Some(summary);
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
                subset_spec,
            } => {
                let count = self.export_rna_read_hits_fasta(
                    &report_id,
                    &path,
                    selection,
                    &selected_record_indices,
                    subset_spec.as_deref(),
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
                gene_ids,
                complete_rule,
                append,
            } => {
                let export = self.export_rna_read_sample_sheet(
                    &path,
                    seq_id.as_deref(),
                    &report_ids,
                    &gene_ids,
                    complete_rule,
                    append,
                )?;
                result.messages.push(format!(
                    "Exported RNA-read sample sheet '{}' with {} report row(s) (genes={}, complete_rule={}, append={})",
                    export.path,
                    export.report_count,
                    export.gene_ids.len(),
                    export.complete_rule.as_str(),
                    export.appended
                ));
            }
            Operation::ExportRnaReadExonPathsTsv {
                report_id,
                path,
                selection,
                selected_record_indices,
                subset_spec,
            } => {
                let export = self.export_rna_read_exon_paths_tsv(
                    &report_id,
                    &path,
                    selection,
                    &selected_record_indices,
                    subset_spec.as_deref(),
                )?;
                result.messages.push(format!(
                    "Exported RNA-read exon paths '{}' to '{}' (selection={}, rows={}, selected_record_indices={})",
                    export.report_id,
                    export.path,
                    export.selection.as_str(),
                    export.row_count,
                    selected_record_indices.len()
                ));
            }
            Operation::ExportRnaReadExonAbundanceTsv {
                report_id,
                path,
                selection,
                selected_record_indices,
                subset_spec,
            } => {
                let export = self.export_rna_read_exon_abundance_tsv(
                    &report_id,
                    &path,
                    selection,
                    &selected_record_indices,
                    subset_spec.as_deref(),
                )?;
                result.messages.push(format!(
                    "Exported RNA-read exon abundance '{}' to '{}' (selection={}, selected_reads={}, exon_rows={}, transition_rows={}, selected_record_indices={})",
                    export.report_id,
                    export.path,
                    export.selection.as_str(),
                    export.selected_read_count,
                    export.exon_row_count,
                    export.transition_row_count,
                    selected_record_indices.len()
                ));
            }
            Operation::ExportRnaReadScoreDensitySvg {
                report_id,
                path,
                scale,
                variant,
                seed_filter_override,
            } => {
                let export = self.export_rna_read_score_density_svg(
                    &report_id,
                    &path,
                    scale,
                    variant,
                    seed_filter_override.as_ref(),
                )?;
                result.messages.push(format!(
                    "Exported RNA-read score-density SVG '{}' to '{}' (scale={}, variant={}, bins={}, max_bin_count={})",
                    export.report_id,
                    export.path,
                    export.scale.as_str(),
                    export.variant.as_str(),
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
                subset_spec,
            } => {
                let export = self.export_rna_read_alignments_tsv(
                    &report_id,
                    &path,
                    selection,
                    limit,
                    &selected_record_indices,
                    subset_spec.as_deref(),
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

                let base = output_id
                    .unwrap_or_else(|| Self::derive_extract_region_default_base(&input, from, to));
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
