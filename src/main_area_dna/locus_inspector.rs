//! Interactive projection of portable locus reports and validated DNA selections.

use super::*;

const EXON: egui::Color32 = egui::Color32::from_rgb(120, 165, 220);
const CDS: egui::Color32 = egui::Color32::from_rgb(37, 99, 175);
const FORWARD: egui::Color32 = egui::Color32::from_rgb(190, 24, 93);
const REVERSE: egui::Color32 = egui::Color32::from_rgb(37, 99, 235);

#[derive(Debug)]
struct EvidenceInterval {
    start: usize,
    end: usize,
    height: f32,
    color: egui::Color32,
}

#[cfg(test)]
mod tests {
    use super::*;

    fn anchored_area() -> (MainAreaDna, GeneLocusEvidenceDisplayReport) {
        let dna = DNAsequence::from_sequence(&"ACGT".repeat(25)).unwrap();
        let mut state = crate::engine::ProjectState::default();
        state.sequences.insert("roi_locus".into(), dna.clone());
        state.metadata.insert(
            "provenance".into(),
            serde_json::json!({"genome_extractions": [{
                "seq_id": "roi_locus", "genome_id": "GRCh38", "chromosome": "chr7",
                "start_1based": 1001, "end_1based": 1100, "anchor_strand": "-",
                "anchor_verified": true, "recorded_at_unix_ms": 123
            }]}),
        );
        let area = MainAreaDna::new(
            dna,
            Some("roi_locus".into()),
            Some(Arc::new(RwLock::new(GentleEngine::from_state(state)))),
        );
        let report = GeneLocusEvidenceDisplayReport {
            schema: gentle_protocol::GENE_LOCUS_EVIDENCE_DISPLAY_SCHEMA.into(),
            seq_id: "roi_locus".into(),
            gene_symbol: "SYNTHETIC".into(),
            sequence_binding: Some(crate::locus_report::sequence_binding(
                &area.dna.read().unwrap(),
                area.active_sequence_anchor_summary().as_ref(),
            )),
            locus_local_start_1based: 1,
            locus_local_end_1based: 100,
            gene_strand: "-".into(),
            local_axis_direction:
                gentle_protocol::GeneLocusLocalAxisDirection::IncreasingLeftToRight,
            ..Default::default()
        };
        (area, report)
    }

    #[test]
    fn locus_import_rejects_changed_bases_and_save_rechecks_staged_binding() {
        let (mut area, report) = anchored_area();
        area.open_genomic_region_manager(Some((10, 21)));
        area.genomic_region_pending_locus_report = Some(Arc::new(report.clone()));
        *area.dna.write().unwrap() = DNAsequence::from_sequence(&"TGCA".repeat(25)).unwrap();
        let document =
            crate::locus_report::LocusDocument::from_json(&serde_json::to_vec(&report).unwrap())
                .unwrap();
        assert!(
            area.load_splicing_locus_document(document)
                .unwrap_err()
                .contains("digest/length")
        );
        area.save_pending_genomic_region_selection();
        assert!(area.genomic_region_status.contains("digest/length"));
        assert!(
            area.engine
                .as_ref()
                .unwrap()
                .read()
                .unwrap()
                .genomic_region_store_snapshot()
                .unwrap()
                .sets
                .is_empty()
        );
        assert!(area.splicing_locus_report.is_none());
    }

    #[test]
    fn locus_presentation_preserves_transcripts_unavailable_lanes_sites_and_reporters() {
        let (_, mut report) = anchored_area();
        report.isoform_evidence.splicing = Some(serde_json::from_value(serde_json::json!({
            "seq_id":"roi_locus", "target_feature_id":0,"group_label":"SYNTHETIC","strand":"+",
            "region_start_1based":1,"region_end_1based":100,"transcript_count":2,"unique_exon_count":2,
            "instruction":"Synthetic test", "unique_exons":[],"matrix_rows":[],"boundaries":[],"junctions":[],"events":[],
            "transcripts": [
                {"transcript_id":"tx1", "transcript_feature_id":1, "label":"tx1","strand":"+","introns":[],"has_target_feature":true,"exons":[{"start_1based":1,"end_1based":20}]},
                {"transcript_id":"tx2", "transcript_feature_id":2, "label":"tx2","strand":"+","introns":[],"has_target_feature":true,"exons":[{"start_1based":5,"end_1based":20}]}]
        })).unwrap());
        report
            .transcript_metrics
            .push(gentle_protocol::GeneLocusTranscriptMetrics {
                transcript_feature_id: 1,
                cds_ranges_local_1based: vec![(10, 20)],
                ..Default::default()
            });
        report.occupancy_groups = serde_json::from_value(serde_json::json!([{"group_id":"cutrun", "lanes":[{"source_id":"missing", "state":"not_prepared"}]}])).unwrap();
        report.regulatory_score_tracks = serde_json::from_value(serde_json::json!([{"track_id":"tf", "forward_scores":[-2,0,4], "score_kind":"llr_log2", "sites":[{"local_start_0based":8,"local_end_0based_exclusive":12,"local_strand":"-"}]}])).unwrap();
        let overlay = gentle_render::GeneLocusEvidenceOverlay {
            legend_items: vec![gentle_render::GeneLocusEvidenceOverlayLegendItem {
                item_id: "spliced".into(),
                label: "Spliced cDNA".into(),
                fill: "#d58a3a".into(),
            }],
            rows: vec![gentle_render::GeneLocusEvidenceOverlayRow {
                row_id: "reporter".into(),
                label: "Spliced reporter".into(),
                segments: vec![gentle_render::GeneLocusEvidenceOverlaySegment {
                    local_start_1based: 5,
                    local_end_1based: 20,
                    fill: "#d58a3a".into(),
                    ..Default::default()
                }],
                ..Default::default()
            }],
            ..Default::default()
        };
        let prepared = LocusPresentation::from_report(&report, Some(overlay));
        assert_eq!(
            prepared
                .lanes
                .iter()
                .filter(|lane| lane.id.starts_with("transcript:"))
                .count(),
            2
        );
        assert_eq!(prepared.lanes[0].intervals[1].color, CDS);
        assert!(
            prepared
                .lanes
                .iter()
                .find(|lane| lane.id.starts_with("occupancy:"))
                .unwrap()
                .detail
                .contains("NotPrepared")
        );
        assert!(
            prepared
                .lanes
                .iter()
                .find(|lane| lane.id.starts_with("occupancy:"))
                .unwrap()
                .intervals
                .is_empty()
        );
        assert_eq!(
            prepared
                .lanes
                .iter()
                .find(|lane| lane.id.starts_with("sites:"))
                .unwrap()
                .intervals[0]
                .start,
            9
        );
        assert_eq!(
            prepared
                .lanes
                .iter()
                .find(|lane| lane.id.starts_with("sites:"))
                .unwrap()
                .intervals[0]
                .color,
            REVERSE
        );
        assert!(
            prepared
                .lanes
                .iter()
                .any(|lane| lane.id == "reporter:reporter")
        );
        assert!(
            prepared
                .legend
                .iter()
                .any(|(label, _)| label == "Spliced cDNA")
        );
        assert_eq!(prepared.score_bounds, vec![(0.0, 4.0)]);
        assert_eq!(MainAreaDna::locus_inspector_color(Some("#é1234"), CDS), CDS);
    }

    #[cfg(feature = "gui-test-support")]
    #[test]
    fn locus_gui_drag_stage_save_reopen_export_uses_negative_anchor() {
        use crate::gui_test_support::{self as semantics, GuiTestSnapshot};
        let (mut area, report) = anchored_area();
        area.load_splicing_locus_document(
            crate::locus_report::LocusDocument::from_json(&serde_json::to_vec(&report).unwrap())
                .unwrap(),
        )
        .unwrap();
        let ctx = egui::Context::default();
        let frame =
            |area: &mut MainAreaDna, events: Vec<egui::Event>, manager: bool| -> GuiTestSnapshot {
                ctx.begin_pass(egui::RawInput {
                    screen_rect: Some(egui::Rect::from_min_size(
                        egui::Pos2::ZERO,
                        egui::vec2(1400.0, 1200.0),
                    )),
                    events,
                    ..Default::default()
                });
                semantics::begin_frame(&ctx);
                if manager {
                    area.render_genomic_region_manager(&ctx);
                } else {
                    crate::egui_compat::show_central_panel_for_test_context(
                        &ctx,
                        egui::CentralPanel::default(),
                        |ui| area.render_splicing_locus_interactive_inspector(ui, &report),
                    );
                }
                crate::egui_compat::discard_test_pass_output(&ctx);
                semantics::snapshot(&ctx)
            };
        let find = |snapshot: &GuiTestSnapshot, id: &str| {
            let item = snapshot
                .items
                .iter()
                .find(|item| item.semantic_id == id)
                .unwrap_or_else(|| panic!("missing {id}: {snapshot:?}"));
            let r = item.rect_logical_points;
            egui::Rect::from_min_max(egui::pos2(r.min_x, r.min_y), egui::pos2(r.max_x, r.max_y))
        };
        let button = |pos, pressed| {
            vec![
                egui::Event::PointerMoved(pos),
                egui::Event::PointerButton {
                    pos,
                    button: egui::PointerButton::Primary,
                    pressed,
                    modifiers: egui::Modifiers::default(),
                },
            ]
        };
        frame(&mut area, vec![], false);
        let snapshot = frame(&mut area, vec![], false);
        let canvas = find(&snapshot, "splicing.locus.inspector.canvas");
        let left = canvas.left() + 250.0;
        let width = canvas.right() - 12.0 - left;
        let from = egui::pos2(left + width * 10.0 / 99.0, canvas.center().y);
        let to = egui::pos2(left + width * 20.0 / 99.0, canvas.center().y);
        frame(&mut area, button(from, true), false);
        frame(&mut area, vec![egui::Event::PointerMoved(to)], false);
        let snapshot = frame(&mut area, button(to, false), false);
        assert_eq!(
            area.splicing_locus_inspector_selection,
            Some((11, 21)),
            "drag must retain the press origin, not the drag-threshold position"
        );
        let stage = find(&snapshot, "splicing.locus.inspector.stage_reporter_region").center();
        frame(&mut area, button(stage, true), false);
        frame(&mut area, button(stage, false), false);
        assert_eq!(area.genomic_region_pending_selection, Some((10, 21)));
        assert!(area.genomic_region_pending_locus_report.is_some());
        area.genomic_region_new_label = "Reviewed synthetic reporter".into();
        area.genomic_region_new_color_hex = "#336699".into();
        frame(&mut area, vec![], true);
        let snapshot = frame(&mut area, vec![], true);
        let save = find(
            &snapshot,
            crate::tutorial_gui_semantics::GENOMIC_REGION_SAVE_PENDING,
        )
        .center();
        frame(&mut area, button(save, true), true);
        frame(&mut area, button(save, false), true);
        assert!(
            area.genomic_region_pending_selection.is_none(),
            "{}",
            area.genomic_region_status
        );
        let engine = area.engine.as_ref().unwrap().read().unwrap();
        let store = engine.genomic_region_store_snapshot().unwrap();
        assert_eq!(store.sets.len(), 1);
        let region = &store.sets[0].regions[0];
        assert_eq!(region.label.as_deref(), Some("Reviewed synthetic reporter"));
        assert_eq!(region.display_color_hex.as_deref(), Some("#336699"));
        assert_eq!(
            region.purpose,
            gentle_protocol::GenomicRegionPurpose::ReporterCandidate
        );
        let temp = tempfile::tempdir().unwrap();
        let project = temp.path().join("reopen.gentle.json");
        engine
            .state()
            .save_to_path(project.to_str().unwrap())
            .unwrap();
        let state = crate::engine::ProjectState::load_from_path(project.to_str().unwrap()).unwrap();
        let mut reopened = GentleEngine::from_state(state);
        assert_eq!(reopened.genomic_region_store_snapshot().unwrap(), store);
        let bed = temp.path().join("reporter.bed");
        reopened
            .apply(Operation::ExportGenomicRegionSet {
                request: gentle_protocol::GenomicRegionExportRequest {
                    set_id: store.sets[0].set_id.clone(),
                    bed_path: Some(bed.display().to_string()),
                    manifest_path: Some(
                        temp.path()
                            .join("reporter.manifest.json")
                            .display()
                            .to_string(),
                    ),
                    json_path: Some(temp.path().join("reporter.json").display().to_string()),
                    include_local_paths: false,
                },
            })
            .unwrap();
        assert!(
            std::fs::read_to_string(bed)
                .unwrap()
                .contains("chr7\t1079\t1090\t"),
            "local 11..21 on the reverse anchor must export genomic 1080..1090 (BED 1079..1090)"
        );
    }
}

#[derive(Debug)]
struct EvidenceLane {
    id: String,
    label: String,
    detail: String,
    intervals: Vec<EvidenceInterval>,
    markers: Vec<(usize, String)>,
    signal_scale: Option<f64>,
}

/// Presentation geometry prepared once per report, not per repaint.
#[derive(Debug, Default)]
pub(super) struct LocusPresentation {
    lanes: Vec<EvidenceLane>,
    score_bounds: Vec<(f64, f64)>,
    legend: Vec<(String, egui::Color32)>,
}

impl LocusPresentation {
    pub(super) fn from_report(
        report: &GeneLocusEvidenceDisplayReport,
        overlay: Option<gentle_render::GeneLocusEvidenceOverlay>,
    ) -> Self {
        let mut result = Self {
            score_bounds: report
                .regulatory_score_tracks
                .iter()
                .map(gentle_render::locus_regulatory_display_bounds)
                .collect(),
            ..Default::default()
        };
        if let Some(splicing) = &report.isoform_evidence.splicing {
            for lane in &splicing.transcripts {
                let mut intervals = lane
                    .exons
                    .iter()
                    .map(|exon| EvidenceInterval {
                        start: exon.start_1based,
                        end: exon.end_1based,
                        height: 0.55,
                        color: EXON,
                    })
                    .collect::<Vec<_>>();
                if let Some(metrics) = report
                    .transcript_metrics
                    .iter()
                    .find(|metrics| metrics.transcript_feature_id == lane.transcript_feature_id)
                {
                    intervals.extend(metrics.cds_ranges_local_1based.iter().map(
                        |&(start, end)| EvidenceInterval {
                            start,
                            end,
                            height: 1.0,
                            color: CDS,
                        },
                    ));
                }
                result.lanes.push(EvidenceLane {
                    id: format!("transcript:{}", lane.transcript_feature_id),
                    label: lane.transcript_id.clone(),
                    detail: format!(
                        "{} | local strand {} | thin blue: exon; thick blue: annotated CDS",
                        lane.label, lane.strand
                    ),
                    intervals,
                    markers: report
                        .codon_markers
                        .iter()
                        .filter(|marker| marker.transcript_id == lane.transcript_id)
                        .map(|marker| {
                            (
                                marker.local_position_1based,
                                format!("{:?}: {}", marker.kind, marker.basis),
                            )
                        })
                        .collect(),
                    signal_scale: None,
                });
            }
        }
        if result.lanes.is_empty() {
            result.lanes.push(EvidenceLane {
                id: "merged_exons".into(), label: "Merged exons".into(),
                detail: "No per-transcript model was supplied; these are merged exon intervals, not CDS assignments".into(),
                intervals: report.isoform_evidence.exon_families.iter().map(|exon| EvidenceInterval { start: exon.local_start_1based, end: exon.local_end_1based, height: 0.55, color: EXON }).collect(),
                markers: vec![], signal_scale: None,
            });
        }
        result
            .legend
            .extend([("Exon".into(), EXON), ("CDS".into(), CDS)]);
        for group in &report.occupancy_groups {
            for (index, lane) in group.lanes.iter().enumerate() {
                let color = match lane.role {
                    gentle_protocol::GeneLocusOccupancyLaneRole::ChromatinContext => {
                        egui::Color32::from_rgb(180, 90, 30)
                    }
                    _ => egui::Color32::from_rgb(30, 130, 105),
                };
                let scale = lane.display_abs_max_score.max(f64::EPSILON);
                let available =
                    lane.state == gentle_protocol::GeneLocusOccupancyLaneState::Available;
                result.lanes.push(EvidenceLane {
                    id: format!("occupancy:{}:{index}", group.group_id),
                    label: format!(
                        "{} [{}]",
                        lane.display_label.as_deref().unwrap_or(&lane.source_id),
                        lane.state.as_str()
                    ),
                    detail: format!(
                        "{} | {:?} | {:?} | assay {} / mark {} / factor {} | source {}",
                        group.label,
                        lane.role,
                        lane.state,
                        lane.assay.as_deref().unwrap_or("unspecified"),
                        lane.mark.as_deref().unwrap_or("unspecified"),
                        lane.factor.as_deref().unwrap_or("unspecified"),
                        lane.source_id
                    ),
                    intervals: if available {
                        lane.lane
                            .intervals
                            .iter()
                            .filter_map(|interval| {
                                interval
                                    .score
                                    .filter(|score| score.is_finite())
                                    .map(|score| EvidenceInterval {
                                        start: interval.local_start_1based,
                                        end: interval.local_end_1based,
                                        height: (score / scale).clamp(-1.0, 1.0) as f32,
                                        color,
                                    })
                            })
                            .collect()
                    } else {
                        vec![]
                    },
                    markers: if available {
                        lane.lane
                            .intervals
                            .iter()
                            .filter(|interval| interval.score.is_none())
                            .map(|interval| {
                                (
                                    interval.local_start_1based,
                                    format!(
                                        "Unscored occupancy interval {}..{}",
                                        interval.local_start_1based, interval.local_end_1based
                                    ),
                                )
                            })
                            .collect()
                    } else {
                        vec![]
                    },
                    signal_scale: (available
                        && lane
                            .lane
                            .intervals
                            .iter()
                            .any(|interval| interval.score.is_some_and(f64::is_finite)))
                    .then_some(scale),
                });
            }
        }
        for track in &report.regulatory_score_tracks {
            if !track.sites.is_empty() {
                result.lanes.push(EvidenceLane {
                    id: format!("sites:{}", track.track_id),
                    label: format!("{} sites", track.label),
                    detail: format!(
                        "Retained predicted site calls | {} | {} | not occupancy evidence",
                        track.score_kind, track.score_units
                    ),
                    intervals: track
                        .sites
                        .iter()
                        .filter(|site| site.score.is_finite() && site.score >= 0.0)
                        .map(|site| EvidenceInterval {
                            start: site.local_start_0based.saturating_add(1),
                            end: site.local_end_0based_exclusive,
                            height: 0.7,
                            color: if site.local_strand == "-"
                                || (site.local_strand.is_empty() && site.strand == "-")
                            {
                                REVERSE
                            } else {
                                FORWARD
                            },
                        })
                        .collect(),
                    markers: vec![],
                    signal_scale: None,
                });
            }
        }
        for track in &report.motif_tracks {
            if !report
                .regulatory_score_tracks
                .iter()
                .any(|score| score.source_ids.contains(&track.motif_id))
            {
                result.lanes.push(EvidenceLane {
                    id: format!("motif:{}", track.motif_id),
                    label: format!("{} hits", track.motif_id),
                    detail: format!(
                        "Legacy retained motif hits | {} | not occupancy evidence",
                        track.score_kind
                    ),
                    intervals: track
                        .top_hits
                        .iter()
                        .filter(|hit| hit.score.is_finite() && hit.score >= 0.0)
                        .map(|hit| EvidenceInterval {
                            start: hit.local_start_0based.saturating_add(1),
                            end: hit.local_end_0based_exclusive,
                            height: 0.7,
                            color: if hit.strand == "-" { REVERSE } else { FORWARD },
                        })
                        .collect(),
                    markers: vec![],
                    signal_scale: None,
                });
            }
        }
        if let Some(overlay) = overlay {
            result
                .legend
                .extend(overlay.legend_items.into_iter().map(|item| {
                    (
                        item.label,
                        MainAreaDna::locus_inspector_color(Some(&item.fill), CDS),
                    )
                }));
            for row in overlay.rows {
                let mut detail = row.detail;
                if let Some(tail) = &row.schematic_tail {
                    detail.push_str(&format!(
                        " | {}: {} (schematic, not genomic scale)",
                        tail.label, tail.detail
                    ));
                }
                let mut markers = row
                    .marker_local_1based
                    .map(|position| (position, row.marker_label.unwrap_or_default()))
                    .into_iter()
                    .collect::<Vec<_>>();
                if let Some(tail) = row.schematic_tail {
                    markers.push((
                        tail.anchor_local_1based,
                        format!("{} (schematic)", tail.label),
                    ));
                }
                result.lanes.push(EvidenceLane {
                    id: format!("reporter:{}", row.row_id),
                    label: row.label,
                    detail,
                    intervals: row
                        .segments
                        .into_iter()
                        .map(|segment| EvidenceInterval {
                            start: segment.local_start_1based,
                            end: segment.local_end_1based,
                            height: 0.8,
                            color: MainAreaDna::locus_inspector_color(Some(&segment.fill), CDS),
                        })
                        .collect(),
                    markers,
                    signal_scale: None,
                });
            }
        }
        result
    }
}

#[derive(Debug, Clone)]
pub(super) struct LocusBindingCache {
    key: (u64, u64, u64, usize),
    result: Result<(), String>,
}

impl MainAreaDna {
    pub(super) fn splicing_locus_svg(&self, report: &GeneLocusEvidenceDisplayReport) -> String {
        if let Some(document) = &self.splicing_locus_imported_document {
            document.render_svg()
        } else {
            gentle_render::render_gene_locus_evidence_with_overlay_svg(report, None)
        }
    }

    pub(super) fn splicing_locus_json(
        &self,
        report: &GeneLocusEvidenceDisplayReport,
    ) -> Result<Vec<u8>, String> {
        if let Some(document) = &self.splicing_locus_imported_document {
            document.to_json()
        } else {
            serde_json::to_vec_pretty(report).map_err(|error| error.to_string())
        }
    }

    pub(super) fn load_splicing_locus_document(
        &mut self,
        document: crate::locus_report::LocusDocument,
    ) -> Result<(), String> {
        let report = document.locus();
        if report.seq_id != self.seq_id.as_deref().unwrap_or_default() {
            return Err("Report sequence identity does not match the active sequence".into());
        }
        if report.sequence_binding.is_some() {
            self.verify_splicing_locus_binding(report)?;
        }
        let report = report.clone();
        self.splicing_locus_imported_document = Some(Arc::new(document));
        self.reset_splicing_locus_inspector(&report);
        self.splicing_locus_report = Some(Arc::new(report.clone()));
        self.splicing_locus_preview_png = None;
        self.cache_splicing_locus_preview(&report)
    }

    pub(super) fn verify_splicing_locus_binding(
        &self,
        report: &GeneLocusEvidenceDisplayReport,
    ) -> Result<(), String> {
        let anchor = self.active_sequence_anchor_summary();
        let dna = self.dna.try_read().map_err(|_| "Sequence is busy")?;
        let binding = crate::locus_report::sequence_binding(&dna, anchor.as_ref());
        crate::locus_report::verify_live_binding(
            report,
            self.seq_id.as_deref().unwrap_or_default(),
            &binding,
        )?;
        if let Some(engine) = &self.engine {
            let engine = engine.try_read().map_err(|_| "Engine is busy")?;
            let source = engine
                .state()
                .sequences
                .get(&report.seq_id)
                .ok_or("Report sequence is no longer present in the engine")?;
            let engine_binding = crate::locus_report::sequence_binding(source, anchor.as_ref());
            crate::locus_report::verify_live_binding(report, &report.seq_id, &engine_binding)?;
        }
        Ok(())
    }

    fn cached_locus_binding(
        &mut self,
        report: &GeneLocusEvidenceDisplayReport,
    ) -> Result<(), String> {
        let engine_revision = match self.engine.as_ref() {
            Some(engine) => engine
                .try_read()
                .map_err(|_| "Engine is busy")?
                .mutation_revision(),
            None => 0,
        };
        let display_revision = self
            .dna_display
            .try_read()
            .map_err(|_| "DNA display is busy")?
            .revision();
        let (generation, len) = {
            let dna = self.dna.try_read().map_err(|_| "Sequence is busy")?;
            (dna.feature_generation(), dna.len())
        };
        let key = (engine_revision, display_revision, generation, len);
        if let Some(cache) = &self.splicing_locus_binding_cache
            && cache.key == key
        {
            return cache.result.clone();
        }
        let result = self.verify_splicing_locus_binding(report);
        self.splicing_locus_binding_cache = Some(LocusBindingCache {
            key,
            result: result.clone(),
        });
        result
    }

    pub(super) fn locus_inspector_color(
        value: Option<&str>,
        fallback: egui::Color32,
    ) -> egui::Color32 {
        let Some(value) = value.and_then(|value| value.strip_prefix('#')) else {
            return fallback;
        };
        if value.len() != 6 || !value.bytes().all(|value| value.is_ascii_hexdigit()) {
            return fallback;
        }
        let parse = |range| u8::from_str_radix(&value[range], 16).ok();
        match (parse(0..2), parse(2..4), parse(4..6)) {
            (Some(r), Some(g), Some(b)) => egui::Color32::from_rgb(r, g, b),
            _ => fallback,
        }
    }

    pub(super) fn render_splicing_locus_interactive_inspector(
        &mut self,
        ui: &mut egui::Ui,
        report: &GeneLocusEvidenceDisplayReport,
    ) {
        let binding = self.cached_locus_binding(report);
        let live = binding.is_ok();
        #[cfg(feature = "gui-test-support")]
        let subject_scope = crate::gui_test_support::pseudonymous_subject_scope(&[self
            .seq_id
            .as_deref()
            .unwrap_or("unnamed")]);
        if let Err(error) = &binding {
            ui.colored_label(egui::Color32::DARK_RED, error);
            ui.small("Historical evidence remains inspectable/exportable. Live bases and DNA selection are disabled.");
        }
        let presentation = self.splicing_locus_presentation.clone();
        let locus_start = report.locus_local_start_1based.max(1);
        let locus_end = report.locus_local_end_1based.max(locus_start);
        self.splicing_locus_inspector_start_1based = self
            .splicing_locus_inspector_start_1based
            .clamp(locus_start, locus_end);
        self.splicing_locus_inspector_end_1based = self
            .splicing_locus_inspector_end_1based
            .clamp(self.splicing_locus_inspector_start_1based, locus_end);
        let view_span = self
            .splicing_locus_inspector_end_1based
            .saturating_sub(self.splicing_locus_inspector_start_1based)
            .saturating_add(1);

        ui.horizontal_wrapped(|ui| {
            ui.strong("Interactive locus inspection");
            let fit = ui.small_button("Fit locus");
            #[cfg(feature = "gui-test-support")]
            crate::gui_test_support::register_response(
                &fit,
                "splicing.locus.inspector.fit",
                crate::tutorial_gui_semantics::WINDOW_SPLICING_EXPERT,
                Some(&subject_scope),
                crate::gui_test_support::GuiTestWidgetKind::Button,
                false,
            );
            if fit.clicked() {
                self.splicing_locus_inspector_start_1based = locus_start;
                self.splicing_locus_inspector_end_1based = locus_end;
            }
            let center = self
                .splicing_locus_inspector_start_1based
                .saturating_add(view_span / 2);
            if ui.small_button("−").on_hover_text("Zoom out").clicked() {
                let span = view_span.saturating_mul(2).min(locus_end - locus_start + 1);
                let start = center.saturating_sub(span / 2).max(locus_start);
                self.splicing_locus_inspector_start_1based = start.min(locus_end + 1 - span);
                self.splicing_locus_inspector_end_1based =
                    self.splicing_locus_inspector_start_1based + span - 1;
            }
            if ui.small_button("+").on_hover_text("Zoom in").clicked() {
                let span = (view_span / 2).max(10).min(view_span);
                let start = center.saturating_sub(span / 2).max(locus_start);
                self.splicing_locus_inspector_start_1based = start.min(locus_end + 1 - span);
                self.splicing_locus_inspector_end_1based =
                    self.splicing_locus_inspector_start_1based + span - 1;
            }
            ui.checkbox(
                &mut self.splicing_locus_inspector_show_sequence,
                "Sequence at nucleotide zoom",
            );
            let open_map =
                ui.add_enabled(live, egui::Button::new("Open this span in DNA map").small());
            #[cfg(feature = "gui-test-support")]
            crate::gui_test_support::register_response(
                &open_map,
                "splicing.locus.inspector.open_in_dna_map",
                crate::tutorial_gui_semantics::WINDOW_SPLICING_EXPERT,
                Some(&subject_scope),
                crate::gui_test_support::GuiTestWidgetKind::Button,
                false,
            );
            if open_map.clicked() {
                match self.verify_splicing_locus_binding(report) {
                    Ok(()) => {
                        self.set_linear_viewport(
                            self.splicing_locus_inspector_start_1based.saturating_sub(1),
                            self.splicing_locus_inspector_end_1based
                                .saturating_sub(self.splicing_locus_inspector_start_1based)
                                + 1,
                        );
                        self.show_map = true;
                        self.show_sequence = true;
                    }
                    Err(error) => self.splicing_locus_status = error,
                }
            }
            ui.monospace(format!(
                "local {}..{} ({} bp)",
                self.splicing_locus_inspector_start_1based,
                self.splicing_locus_inspector_end_1based,
                view_span
            ));
            ui.label("start");
            ui.add(
                egui::DragValue::new(&mut self.splicing_locus_inspector_start_1based)
                    .range(locus_start..=locus_end),
            );
            ui.label("end");
            ui.add(
                egui::DragValue::new(&mut self.splicing_locus_inspector_end_1based)
                    .range(locus_start..=locus_end),
            );
        });
        self.splicing_locus_inspector_start_1based = self
            .splicing_locus_inspector_start_1based
            .clamp(locus_start, locus_end);
        self.splicing_locus_inspector_end_1based = self
            .splicing_locus_inspector_end_1based
            .clamp(self.splicing_locus_inspector_start_1based, locus_end);
        let view_span = self
            .splicing_locus_inspector_end_1based
            .saturating_sub(self.splicing_locus_inspector_start_1based)
            .saturating_add(1);

        egui::CollapsingHeader::new("Visible evidence lanes")
            .default_open(true)
            .show(ui, |ui| {
                ui.horizontal_wrapped(|ui| {
                    for lane in &presentation.lanes {
                        let mut visible = !self.splicing_locus_inspector_hidden_scores.contains(&lane.id);
                        if ui.checkbox(&mut visible, &lane.label).on_hover_text(&lane.detail).changed() {
                            if visible { self.splicing_locus_inspector_hidden_scores.remove(&lane.id); }
                            else { self.splicing_locus_inspector_hidden_scores.insert(lane.id.clone()); }
                        }
                    }
                    for track in &report.regulatory_score_tracks {
                        let mut visible = !self
                            .splicing_locus_inspector_hidden_scores
                            .contains(&track.track_id);
                        if ui.checkbox(&mut visible, &track.label).changed() {
                            if visible {
                                self.splicing_locus_inspector_hidden_scores.remove(&track.track_id);
                            } else {
                                self.splicing_locus_inspector_hidden_scores
                                    .insert(track.track_id.clone());
                            }
                        }
                    }
                    if let Some(ensembl) = report.ensembl_regulation.as_ref() {
                        for row in &ensembl.rows {
                            let mut visible = !self
                                .splicing_locus_inspector_hidden_ensembl
                                .contains(&row.feature_id);
                            if ui
                                .checkbox(
                                    &mut visible,
                                    format!("{} ({})", row.feature_id, row.feature_type),
                                )
                                .changed()
                            {
                                if visible {
                                    self.splicing_locus_inspector_hidden_ensembl
                                        .remove(&row.feature_id);
                                } else {
                                    self.splicing_locus_inspector_hidden_ensembl
                                        .insert(row.feature_id.clone());
                                }
                            }
                        }
                    }
                });
                ui.small(
                    "Visibility is an inspection choice only. Add or remove score sources in the composition request, then recompute to change the scientific report.",
                );
            });

        let visible_scores = report
            .regulatory_score_tracks
            .iter()
            .enumerate()
            .filter(|(_, track)| {
                !self
                    .splicing_locus_inspector_hidden_scores
                    .contains(&track.track_id)
            })
            .collect::<Vec<_>>();
        let visible_lanes = presentation
            .lanes
            .iter()
            .filter(|lane| {
                !self
                    .splicing_locus_inspector_hidden_scores
                    .contains(&lane.id)
            })
            .collect::<Vec<_>>();
        ui.horizontal_wrapped(|ui| {
            for (label, color) in &presentation.legend { ui.colored_label(*color, label); }
            ui.colored_label(egui::Color32::from_rgb(194, 65, 12), "Selection / saved ROI (custom colour)");
            if report.ensembl_regulation.is_some() { ui.colored_label(egui::Color32::from_rgb(124, 58, 237), "Ensembl annotation"); }
            ui.colored_label(FORWARD, "Local +");
            ui.colored_label(REVERSE, "Local -");
            ui.small("Green: occupancy; brown: chromatin context. Hover a lane for source/state; markers show reported translation or reporter boundaries.");
        });
        let ensembl_rows = report
            .ensembl_regulation
            .iter()
            .flat_map(|evidence| evidence.rows.iter())
            .filter(|row| {
                !self
                    .splicing_locus_inspector_hidden_ensembl
                    .contains(&row.feature_id)
                    && row.displayed_local_end_1based >= self.splicing_locus_inspector_start_1based
                    && row.displayed_local_start_1based <= self.splicing_locus_inspector_end_1based
            })
            .collect::<Vec<_>>();
        let max_sequence_bases = ((ui.available_width() - 262.0).max(90.0) / 9.0) as usize;
        let sequence_lane =
            live && self.splicing_locus_inspector_show_sequence && view_span <= max_sequence_bases;
        let lane_height = 62.0;
        let canvas_height = 64.0
            + visible_lanes.len() as f32 * 42.0
            + visible_scores.len() as f32 * lane_height
            + if sequence_lane { 34.0 } else { 0.0 };
        let (response, painter) = ui.allocate_painter(
            egui::vec2(ui.available_width().max(360.0), canvas_height),
            egui::Sense::click_and_drag(),
        );
        let rect = response.rect;
        #[cfg(feature = "gui-test-support")]
        crate::gui_test_support::register_rect(
            ui.ctx().clone(),
            "splicing.locus.inspector.canvas",
            crate::tutorial_gui_semantics::WINDOW_SPLICING_EXPERT,
            Some(&subject_scope),
            crate::gui_test_support::GuiTestWidgetKind::Row,
            rect,
            true,
            true,
            response.hovered(),
            Some(if live { "ready" } else { "historical" }),
        );
        painter.rect_filled(rect, 4.0, egui::Color32::from_rgb(248, 250, 252));
        let plot_left = rect.left() + 250.0;
        let plot_right = rect.right() - 12.0;
        let increasing = !gentle_render::gene_locus_local_axis_decreases(report);
        let view_start = self.splicing_locus_inspector_start_1based;
        let view_end = self.splicing_locus_inspector_end_1based;
        let x_for = |position: usize| {
            let fraction = position.saturating_sub(view_start) as f32
                / view_end.saturating_sub(view_start).max(1) as f32;
            if increasing {
                egui::lerp(plot_left..=plot_right, fraction.clamp(0.0, 1.0))
            } else {
                egui::lerp(plot_right..=plot_left, fraction.clamp(0.0, 1.0))
            }
        };
        let position_for_x = |x: f32| {
            let raw = ((x - plot_left) / (plot_right - plot_left).max(1.0)).clamp(0.0, 1.0);
            let fraction = if increasing { raw } else { 1.0 - raw };
            view_start + (fraction * view_end.saturating_sub(view_start) as f32).round() as usize
        };

        for row in &report.saved_region_overlays {
            if row.local_end_1based < view_start || row.local_start_1based > view_end {
                continue;
            }
            let x1 = x_for(row.local_start_1based.max(view_start));
            let x2 = x_for(row.local_end_1based.min(view_end));
            let color = Self::locus_inspector_color(
                row.display_color_hex.as_deref(),
                egui::Color32::from_rgb(194, 65, 12),
            );
            painter.rect_filled(
                egui::Rect::from_x_y_ranges(x1.min(x2)..=x1.max(x2), rect.top()..=rect.bottom()),
                0.0,
                color.gamma_multiply(0.14),
            );
        }

        let axis_y = rect.top() + 25.0;
        painter.line_segment(
            [
                egui::pos2(plot_left, axis_y),
                egui::pos2(plot_right, axis_y),
            ],
            egui::Stroke::new(1.0, egui::Color32::GRAY),
        );
        painter.text(
            egui::pos2(rect.left() + 4.0, axis_y),
            egui::Align2::LEFT_CENTER,
            "shared local DNA axis (bp)",
            egui::FontId::monospace(10.0),
            egui::Color32::DARK_GRAY,
        );
        for position in [
            view_start,
            view_start + (view_end - view_start) / 2,
            view_end,
        ] {
            painter.text(
                egui::pos2(x_for(position), axis_y - 10.0),
                egui::Align2::CENTER_CENTER,
                position.to_string(),
                egui::FontId::monospace(10.0),
                egui::Color32::DARK_GRAY,
            );
        }
        for row in ensembl_rows {
            let x1 = x_for(row.displayed_local_start_1based.max(view_start));
            let x2 = x_for(row.displayed_local_end_1based.min(view_end));
            painter.rect_filled(
                egui::Rect::from_min_max(
                    egui::pos2(x1.min(x2), axis_y + 11.0),
                    egui::pos2(x1.max(x2).max(x1.min(x2) + 2.0), axis_y + 18.0),
                ),
                1.0,
                egui::Color32::from_rgb(124, 58, 237),
            );
        }

        let mut lane_y = rect.top() + 61.0;
        for lane in visible_lanes {
            let lane_rect = egui::Rect::from_min_max(
                egui::pos2(rect.left(), lane_y - 18.0),
                egui::pos2(rect.right(), lane_y + 20.0),
            );
            if !ui.is_rect_visible(lane_rect) {
                lane_y += 42.0;
                continue;
            }
            let label_painter = painter.with_clip_rect(lane_rect.intersect(
                egui::Rect::from_x_y_ranges(rect.left()..=plot_left - 40.0, rect.y_range()),
            ));
            label_painter.text(
                egui::pos2(rect.left() + 4.0, lane_y - 6.0),
                egui::Align2::LEFT_CENTER,
                &lane.label,
                egui::FontId::monospace(10.0),
                egui::Color32::DARK_GRAY,
            );
            // Empty requested lanes retain their state instead of suggesting zero signal.
            if lane.intervals.is_empty() {
                label_painter.text(
                    egui::pos2(rect.left() + 4.0, lane_y + 8.0),
                    egui::Align2::LEFT_CENTER,
                    &lane.detail,
                    egui::FontId::monospace(9.0),
                    egui::Color32::DARK_RED,
                );
            }
            let hover = ui.interact(lane_rect, response.id.with(&lane.id), egui::Sense::hover());
            hover.on_hover_text(&lane.detail);
            if let Some(scale) = lane.signal_scale {
                for (value, y) in [
                    (scale, lane_y - 13.0),
                    (0.0, lane_y),
                    (-scale, lane_y + 13.0),
                ] {
                    painter.text(
                        egui::pos2(plot_left - 5.0, y),
                        egui::Align2::RIGHT_CENTER,
                        format!("{value:.1}"),
                        egui::FontId::monospace(8.0),
                        egui::Color32::GRAY,
                    );
                }
            }
            let in_view = lane
                .intervals
                .iter()
                .filter(|interval| interval.end >= view_start && interval.start <= view_end)
                .collect::<Vec<_>>();
            if let (Some(start), Some(end)) = (
                in_view.iter().map(|interval| interval.start).min(),
                in_view.iter().map(|interval| interval.end).max(),
            ) {
                painter.line_segment(
                    [
                        egui::pos2(x_for(start), lane_y),
                        egui::pos2(x_for(end), lane_y),
                    ],
                    egui::Stroke::new(0.5, egui::Color32::LIGHT_GRAY),
                );
            }
            for interval in in_view {
                let x1 = x_for(interval.start.max(view_start));
                let x2 = x_for(interval.end.min(view_end));
                let (top, bottom) = if lane.signal_scale.is_some() {
                    let tip = lane_y - interval.height * 13.0;
                    (tip.min(lane_y), tip.max(lane_y))
                } else {
                    (
                        lane_y - interval.height * 8.0,
                        lane_y + interval.height * 8.0,
                    )
                };
                painter.rect_filled(
                    egui::Rect::from_min_max(
                        egui::pos2(x1.min(x2), top),
                        egui::pos2((x1.max(x2)).max(x1.min(x2) + 1.0), bottom),
                    ),
                    0.0,
                    interval.color,
                );
            }
            for (position, label) in &lane.markers {
                if !(view_start..=view_end).contains(position) {
                    continue;
                }
                let x = x_for(*position);
                if label.starts_with("LUC") {
                    painter.text(
                        egui::pos2(x + 3.0, lane_y - 13.0),
                        egui::Align2::LEFT_BOTTOM,
                        "LUC (schematic)",
                        egui::FontId::monospace(8.0),
                        egui::Color32::DARK_RED,
                    );
                }
                painter.line_segment(
                    [egui::pos2(x, lane_y - 13.0), egui::pos2(x, lane_y + 13.0)],
                    egui::Stroke::new(1.5, egui::Color32::BLACK),
                );
                ui.interact(
                    egui::Rect::from_center_size(egui::pos2(x, lane_y), egui::vec2(6.0, 26.0)),
                    response.id.with((&lane.id, position, label)),
                    egui::Sense::hover(),
                )
                .on_hover_text(label);
            }
            lane_y += 42.0;
        }
        for (index, track) in visible_scores {
            let (min, max) = presentation
                .score_bounds
                .get(index)
                .copied()
                .unwrap_or((0.0, track.display_scale_max.max(1.0)));
            let score_rect = egui::Rect::from_min_max(
                egui::pos2(rect.left(), lane_y - 20.0),
                egui::pos2(rect.right(), lane_y + 35.0),
            );
            if !ui.is_rect_visible(score_rect) {
                lane_y += lane_height;
                continue;
            }
            let label_painter = painter.with_clip_rect(score_rect.intersect(
                egui::Rect::from_x_y_ranges(rect.left()..=plot_left - 40.0, rect.y_range()),
            ));
            label_painter.text(
                egui::pos2(rect.left() + 4.0, lane_y + 12.0),
                egui::Align2::LEFT_CENTER,
                gentle_render::locus_regulatory_axis_label(&track.score_kind, &track.score_units),
                egui::FontId::monospace(9.0),
                egui::Color32::DARK_GRAY,
            );
            label_painter.text(
                egui::pos2(rect.left() + 4.0, lane_y),
                egui::Align2::LEFT_CENTER,
                &track.label,
                egui::FontId::monospace(9.0),
                egui::Color32::DARK_GRAY,
            );
            ui.interact(
                score_rect,
                response.id.with((&track.track_id, "score")),
                egui::Sense::hover(),
            )
            .on_hover_text(format!(
                "{} | {} | {} | {:?} | {}",
                track.label, track.score_kind, track.score_units, track.state, track.provenance
            ));
            for (value, y) in [
                (min, lane_y + 14.0),
                ((min + max) / 2.0, lane_y + 1.5),
                (max, lane_y - 11.0),
            ] {
                painter.text(
                    egui::pos2(plot_left - 5.0, y),
                    egui::Align2::RIGHT_CENTER,
                    format!("{value:.1}"),
                    egui::FontId::monospace(8.0),
                    egui::Color32::GRAY,
                );
            }
            painter.line_segment(
                [
                    egui::pos2(plot_left, lane_y + 14.0),
                    egui::pos2(plot_right, lane_y + 14.0),
                ],
                egui::Stroke::new(0.5, egui::Color32::LIGHT_GRAY),
            );
            let scale = (max - min).max(f64::EPSILON);
            let stride = track.stride_bp.max(1);
            let score_points = |scores: &[f64]| {
                let first = view_start
                    .saturating_sub(track.track_start_0based.saturating_add(1))
                    .div_ceil(stride)
                    .min(scores.len());
                let last = (view_end.saturating_sub(track.track_start_0based.saturating_add(1))
                    / stride)
                    .saturating_add(1)
                    .min(scores.len());
                let mut segments = Vec::new();
                let mut segment = Vec::new();
                for (index, score) in scores.iter().enumerate().take(last).skip(first) {
                    let local = track
                        .track_start_0based
                        .saturating_add(index.saturating_mul(stride))
                        .saturating_add(1);
                    if !score.is_finite() || *score < 0.0 || local > view_end {
                        if !segment.is_empty() {
                            segments.push(std::mem::take(&mut segment));
                        }
                        continue;
                    }
                    let fraction = ((*score - min) / scale).clamp(0.0, 1.0);
                    segment.push(egui::pos2(
                        x_for(local),
                        lane_y + 14.0 - fraction as f32 * 25.0,
                    ));
                }
                if !segment.is_empty() {
                    segments.push(segment);
                }
                segments
            };
            for (scores, color) in [
                (&track.forward_scores, FORWARD),
                (&track.reverse_scores, REVERSE),
            ] {
                for points in score_points(scores) {
                    if points.len() == 1 {
                        painter.circle_filled(points[0], 1.0, color);
                    } else {
                        painter.add(egui::Shape::line(points, egui::Stroke::new(1.0, color)));
                    }
                }
            }
            lane_y += lane_height;
        }

        if sequence_lane {
            let sequence = self
                .dna
                .read()
                .ok()
                .and_then(|dna| dna.get_range_safe(view_start.saturating_sub(1)..view_end));
            if let Some(sequence) = sequence {
                painter.text(
                    egui::pos2(rect.left() + 4.0, lane_y),
                    egui::Align2::LEFT_CENTER,
                    "loaded DNA (local + bases)",
                    egui::FontId::monospace(9.0),
                    egui::Color32::DARK_GRAY,
                );
                for (offset, base) in sequence.iter().enumerate() {
                    let local = view_start + offset;
                    painter.text(
                        egui::pos2(x_for(local), lane_y),
                        egui::Align2::CENTER_CENTER,
                        (*base as char).to_ascii_uppercase(),
                        egui::FontId::monospace(10.0),
                        egui::Color32::BLACK,
                    );
                }
            }
        }

        if live
            && response.drag_started()
            && let Some(origin) = ui.input(|input| input.pointer.press_origin())
            && origin.x >= plot_left
            && origin.x <= plot_right
            && rect.intersect(ui.clip_rect()).contains(origin)
        {
            self.splicing_locus_inspector_drag_anchor = Some(position_for_x(origin.x));
        }
        if live
            && response.dragged()
            && let Some(anchor) = self.splicing_locus_inspector_drag_anchor
            && let Some(origin) = response.interact_pointer_pos()
        {
            let current = position_for_x(origin.x);
            self.splicing_locus_inspector_selection =
                Some((anchor.min(current), anchor.max(current)));
        }
        if response.drag_stopped() {
            self.splicing_locus_inspector_drag_anchor = None;
        }
        if let Some(pointer) = response.hover_pos()
            && pointer.x >= plot_left
        {
            response.clone().on_hover_text(format!(
                "Shared local DNA coordinate: {} bp",
                position_for_x(pointer.x)
            ));
        }
        if let Some((start, end)) = self.splicing_locus_inspector_selection {
            let x1 = x_for(start);
            let x2 = x_for(end);
            painter.rect_filled(
                egui::Rect::from_x_y_ranges(x1.min(x2)..=x1.max(x2), rect.top()..=rect.bottom()),
                0.0,
                egui::Color32::from_rgba_unmultiplied(234, 88, 12, 42),
            );
            ui.horizontal_wrapped(|ui| {
                ui.label(format!(
                    "Candidate cloning region: local {start}..{end} ({} bp)",
                    end - start + 1
                ));
                let stage = ui.add_enabled(
                    live,
                    egui::Button::new("Stage as persistent reporter region"),
                );
                #[cfg(feature = "gui-test-support")]
                crate::gui_test_support::register_response(
                    &stage,
                    "splicing.locus.inspector.stage_reporter_region",
                    crate::tutorial_gui_semantics::WINDOW_SPLICING_EXPERT,
                    Some(&subject_scope),
                    crate::gui_test_support::GuiTestWidgetKind::Button,
                    false,
                );
                if stage.clicked() {
                    match self.verify_splicing_locus_binding(report) {
                        Ok(()) => {
                            self.genomic_region_new_purpose =
                                gentle_protocol::GenomicRegionPurpose::ReporterCandidate;
                            self.genomic_region_new_label = format!(
                                "{} reporter candidate {}-{}",
                                report.gene_symbol, start, end
                            );
                            self.genomic_region_new_color_hex = "#C2410C".to_string();
                            self.open_genomic_region_manager(Some((start - 1, end)));
                            self.genomic_region_pending_locus_report =
                                Some(Arc::new(report.clone()));
                        }
                        Err(error) => self.splicing_locus_status = error,
                    }
                }
                if ui.small_button("Clear selection").clicked() {
                    self.splicing_locus_inspector_selection = None;
                }
            });
        } else {
            ui.small("Drag across the shared locus axis to define a candidate cloning interval.");
        }
    }
}
