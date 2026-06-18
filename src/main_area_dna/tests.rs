use super::{
    DnaPresentationMode, FeatureCopyPayloadKind, LINEAR_TOPOLOGY_SWITCH_MAX_INITIAL_SPAN_BP,
    MainAreaDna, PcrDesignerMode, PcrPaintRole, PrimaryMapMode, QpcrTranscriptIntentUiMode,
    RnaReadConcatemerSubsetMode, RnaReadInterpretOpsUiState, RnaReadTask, RnaReadTaskMessage,
    RnaReadTaskOutcome, SPLICING_ATTRACT_EAGER_BOUNDARY_THRESHOLD,
    SequencingConfirmationOverviewSelection, SequencingConfirmationReviewFocusKind,
    SplicingIntronSignalKey, SplicingIntronSignalRow, ViewSvgExportProfile,
};
use crate::{
    dna_display::{ConstructReasoningOverlay, ConstructReasoningOverlaySpan, Selection},
    dna_sequence::DNAsequence,
    engine::{
        AlternativePromoterComparisonRow, AnnotationCandidateSummary, AttractRegionClass,
        AttractSpeciesMatchMode, AttractSplicingEvidenceHitRow, AttractSplicingEvidenceSettings,
        AttractSplicingEvidenceView, ConstructRole, CutRunAlignConfig, CutRunInputFormat,
        CutRunReadLayout, CutRunRegulatoryTfbsConfirmationStatus, CutRunSeedFilterConfig,
        DotplotMode, DotplotOverlayAnchorExonRef, DotplotOverlayXAxisMode, DotplotView,
        EditableStatus, Engine, EvidenceClass, FlexibilityModel, FlexibilityTrack, GentleEngine,
        LinearSequenceLetterLayoutMode, OpResult, Operation, PairwiseAlignmentMode,
        PrimerDesignBackend, PrimerDesignPairConstraint, PrimerDesignProgress,
        PrimerDesignSideConstraint, ProbeRegionEvidenceInterpretationReport, ProjectState,
        PromoterExpressionEvidenceInput, PromoterReporterCandidateSet,
        ProtocolCartoonPreviewTelemetry, QpcrTranscriptSpecificityEvidence,
        QpcrTranscriptTargeting, QpcrTranscriptTargetingMode, RestrictionCloningPcrHandoffMode,
        RestrictionEnzymeDisplayMode, RnaReadAlignmentEffect, RnaReadAlignmentInspection,
        RnaReadAlignmentInspectionRow, RnaReadHitSelection, RnaReadInputFormat,
        RnaReadInterpretProgress, RnaReadInterpretationHit, RnaReadInterpretationProfile,
        RnaReadInterpretationReport, RnaReadInterpretationReportSummary, RnaReadIsoformSupportRow,
        RnaReadIsoformTriageBin, RnaReadMappingHit, RnaReadOriginMode, RnaReadReportMode,
        RnaReadScoreDensityVariant, RnaReadSeedFilterConfig, SequenceAlignmentReport,
        SequenceGenomeAnchorSummary, SequencingConfirmationReadResult,
        SequencingConfirmationReport, SequencingConfirmationStatus,
        SequencingConfirmationTargetKind, SequencingConfirmationTargetResult,
        SequencingConfirmationVariantRow, SequencingPrimerOrientation,
        SequencingPrimerOverlayReport, SequencingPrimerOverlaySuggestion,
        SequencingPrimerProblemKind, SequencingPrimerProposalRow, SequencingTraceChannelData,
        SequencingTraceFormat, SequencingTraceImportReport, SequencingTraceRecord,
        SplicingScopePreset, TfbsScoreTrackReport, TfbsScoreTrackValueKind,
        TfbsTrackSimilarityRankingMetric, TfbsTrackSimilarityReport, VariantPromoterContextReport,
        parse_required_usize_or_formula_text_on_sequence,
    },
    engine_shell::ShellCommand,
    enzymes::active_restriction_enzymes,
    feature_expert::{
        FeatureExpertView, IsoformArchitectureExpertView, RestrictionSiteExpertView,
        SplicingBoundaryMarker, SplicingExonSummary, SplicingExpertView, SplicingIntronSignal,
        SplicingJunctionArc, SplicingRange, SplicingTranscriptLane,
    },
    linear_base_routing::{LinearBaseRenderMode, LinearBaseRoutePolicy},
    protocol_cartoon::pcr_oe_substitution_geometry_bindings,
};
use eframe::egui;
use gb_io::seq::{Feature, Location};
use serde_json::json;
use std::{
    collections::{BTreeMap, BTreeSet},
    fs,
    path::Path,
    sync::{Arc, Mutex, RwLock, atomic::AtomicBool, mpsc},
    time::{Duration, Instant},
};
use tempfile::{TempDir, tempdir};

#[test]
fn cutrun_tfbs_confirmation_label_covers_occupancy_support_statuses() {
    assert_eq!(
        MainAreaDna::cutrun_tfbs_confirmation_label(
            CutRunRegulatoryTfbsConfirmationStatus::Confirmed
        ),
        "confirmed"
    );
    assert_eq!(
        MainAreaDna::cutrun_tfbs_confirmation_label(CutRunRegulatoryTfbsConfirmationStatus::Nearby),
        "nearby"
    );
    assert_eq!(
        MainAreaDna::cutrun_tfbs_confirmation_label(CutRunRegulatoryTfbsConfirmationStatus::Absent),
        "absent"
    );
    assert_eq!(
        MainAreaDna::cutrun_tfbs_confirmation_label(
            CutRunRegulatoryTfbsConfirmationStatus::MotifPoor
        ),
        "motif-poor"
    );
    assert_eq!(
        MainAreaDna::cutrun_tfbs_confirmation_label(
            CutRunRegulatoryTfbsConfirmationStatus::Unconfirmed
        ),
        "unconfirmed"
    );
}

fn make_feature(kind: &str, qualifiers: Vec<(&str, &str)>) -> Feature {
    Feature {
        kind: kind.to_string().into(),
        location: Location::simple_range(0, 10),
        qualifiers: qualifiers
            .into_iter()
            .map(|(key, value)| (key.to_string().into(), Some(value.to_string())))
            .collect(),
    }
}

fn primary_press_input(pos: egui::Pos2, pressed: bool) -> egui::RawInput {
    egui::RawInput {
        events: vec![
            egui::Event::PointerMoved(pos),
            egui::Event::PointerButton {
                pos,
                button: egui::PointerButton::Primary,
                pressed,
                modifiers: egui::Modifiers::default(),
            },
        ],
        ..Default::default()
    }
}

fn with_linear_map_response(
    ctx: &egui::Context,
    raw_input: egui::RawInput,
    rect: egui::Rect,
    inspect: impl FnOnce(&egui::Response),
) {
    ctx.begin_pass(raw_input);
    crate::egui_compat::show_central_panel_for_test_context(
        ctx,
        egui::CentralPanel::default(),
        |ui| {
            let response = ui.interact(
                rect,
                ui.make_persistent_id("linear_map_drag_origin_guard_test"),
                egui::Sense::click_and_drag(),
            );
            inspect(&response);
        },
    );
    let _ = ctx.end_pass();
}

fn probe_region_validation_fixture_dir() -> String {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("test_files/fixtures/probe_region_outputs/clariom_e_mtab_14704_tp73_validation")
        .to_string_lossy()
        .to_string()
}

fn add_tp73_probe_region_validation_transcripts(dna: &mut DNAsequence) {
    dna.features_mut().push(Feature {
        kind: "mRNA".into(),
        location: Location::Join(vec![
            Location::simple_range(4, 43),
            Location::simple_range(54, 80),
            Location::simple_range(104, 135),
        ]),
        qualifiers: vec![
            ("gene".into(), Some("TP73".to_string())),
            ("transcript_id".into(), Some("TP73-201".to_string())),
            ("label".into(), Some("TP73-201".to_string())),
            ("strand".into(), Some("+".to_string())),
        ],
    });
    dna.features_mut().push(Feature {
        kind: "mRNA".into(),
        location: Location::Join(vec![
            Location::simple_range(4, 43),
            Location::simple_range(74, 100),
            Location::simple_range(144, 175),
        ]),
        qualifiers: vec![
            ("gene".into(), Some("TP73".to_string())),
            ("transcript_id".into(), Some("TP73-202".to_string())),
            ("label".into(), Some("TP73-202".to_string())),
            ("strand".into(), Some("+".to_string())),
        ],
    });
}

fn make_tp73_probe_region_validation_engine() -> GentleEngine {
    let mut dna =
        DNAsequence::from_sequence(&"A".repeat(256)).expect("valid TP73 validation sequence");
    add_tp73_probe_region_validation_transcripts(&mut dna);
    let mut state = ProjectState::default();
    state.sequences.insert("array_slice".to_string(), dna);
    state.metadata.insert(
        "provenance".to_string(),
        json!({
            "genome_extractions": [
                {
                    "seq_id": "array_slice",
                    "genome_id": "hg38",
                    "chromosome": "chr1",
                    "start_1based": 3652516,
                    "end_1based": 3736201,
                    "anchor_strand": "+",
                    "anchor_verified": true,
                    "recorded_at_unix_ms": 123
                }
            ]
        }),
    );
    GentleEngine::from_state(state)
}

fn make_tp73_probe_region_validation_area() -> MainAreaDna {
    let engine = make_tp73_probe_region_validation_engine();
    let dna = engine
        .state()
        .sequences
        .get("array_slice")
        .expect("array slice")
        .clone();
    let mut area = MainAreaDna::new(dna, Some("array_slice".to_string()), None);
    area.engine = Some(Arc::new(RwLock::new(engine)));
    area.probe_region_output_dir = probe_region_validation_fixture_dir();
    area.probe_region_projection_seq_id = "array_slice".to_string();
    area
}

#[test]
fn linear_map_drag_handler_rejects_press_origin_outside_map() {
    let mut dna = DNAsequence::from_sequence(&"ACGT".repeat(64)).expect("sequence");
    dna.set_circular(false);
    let mut area = MainAreaDna::new(dna, None, None);
    let ctx = egui::Context::default();
    let map_rect = egui::Rect::from_min_size(egui::pos2(40.0, 40.0), egui::vec2(240.0, 96.0));
    let outside_press = egui::pos2(12.0, 70.0);
    let inside_drag = egui::pos2(80.0, 70.0);

    with_linear_map_response(
        &ctx,
        primary_press_input(outside_press, true),
        map_rect,
        |_| {},
    );
    area.linear_drag_selection_anchor_bp = Some(3);
    area.linear_pan_drag_origin_bp = Some((5, 7.0));
    area.pcr_paint_drag_interval = Some((PcrPaintRole::Roi, 3, 4));

    with_linear_map_response(
        &ctx,
        egui::RawInput {
            events: vec![egui::Event::PointerMoved(inside_drag)],
            ..Default::default()
        },
        map_rect,
        |response| {
            assert!(response.rect.contains(inside_drag));
            area.handle_linear_map_drag_for_pcr_paint(response, &ctx);
        },
    );

    assert_eq!(area.linear_drag_selection_anchor_bp, None);
    assert_eq!(area.linear_pan_drag_origin_bp, None);
    assert_eq!(area.pcr_paint_drag_interval, None);
}

#[test]
fn linear_map_drag_handler_rejects_drag_owned_by_other_widget() {
    let mut dna = DNAsequence::from_sequence(&"ACGT".repeat(64)).expect("sequence");
    dna.set_circular(false);
    let mut area = MainAreaDna::new(dna, None, None);
    let ctx = egui::Context::default();
    let map_rect = egui::Rect::from_min_size(egui::pos2(40.0, 40.0), egui::vec2(240.0, 96.0));
    let press_inside_map = egui::pos2(80.0, 70.0);

    with_linear_map_response(
        &ctx,
        primary_press_input(press_inside_map, true),
        map_rect,
        |response| {
            assert!(response.rect.contains(press_inside_map));
            ctx.set_dragged_id(egui::Id::new("hosted_window_resize_drag_owner"));
            area.linear_drag_selection_anchor_bp = Some(3);
            area.linear_pan_drag_origin_bp = Some((5, 7.0));
            area.pcr_paint_drag_interval = Some((PcrPaintRole::Roi, 3, 4));
            area.handle_linear_map_drag_for_pcr_paint(response, &ctx);
        },
    );

    assert_eq!(area.linear_drag_selection_anchor_bp, None);
    assert_eq!(area.linear_pan_drag_origin_bp, None);
    assert_eq!(area.pcr_paint_drag_interval, None);
}

#[test]
fn sequence_description_panel_text_unwraps_prose_but_keeps_metadata_rows() {
    let lines = vec![
        "REFSEQ INFORMATION: The reference sequence is identical to".to_string(),
        "CM000663.2.".to_string(),
        String::new(),
        "The DNA sequence is composed of genomic sequence, primarily".to_string(),
        "finished clones that were sequenced as part of the Human Genome".to_string(),
        "Project.".to_string(),
        String::new(),
        "##Genome-Annotation-Data-START##".to_string(),
        "Annotation Pipeline         :: NCBI eukaryotic genome annotation".to_string(),
        "pipeline".to_string(),
        "Annotation Method           :: Best-placed RefSeq; Gnomon;".to_string(),
        "RefSeqFE; cmsearch; tRNAscan-SE".to_string(),
        "##Genome-Annotation-Data-END##".to_string(),
    ];

    let text = MainAreaDna::sequence_description_panel_text(&lines);

    assert!(text.contains("identical to CM000663.2."));
    assert!(
        text.contains(
            "The DNA sequence is composed of genomic sequence, primarily finished clones"
        )
    );
    assert!(
        text.contains("Annotation Pipeline         :: NCBI eukaryotic genome annotation pipeline")
    );
    assert!(text.contains(
        "Annotation Method           :: Best-placed RefSeq; Gnomon; RefSeqFE; cmsearch; tRNAscan-SE"
    ));
    assert!(!text.contains("identical to\nCM000663.2."));
}

#[test]
fn feature_copy_identifier_prefers_transcript_id() {
    let feature = make_feature(
        "mRNA",
        vec![
            ("gene", "TP73"),
            ("transcript_id", "NM_005427.4"),
            ("product", "tumor protein p73, transcript variant 1"),
        ],
    );

    assert_eq!(
        MainAreaDna::feature_copy_payload(&feature, FeatureCopyPayloadKind::Identifier),
        "NM_005427.4"
    );
}

#[test]
fn feature_copy_description_prefers_product_text() {
    let feature = make_feature(
        "mRNA",
        vec![
            ("transcript_id", "NM_005427.4"),
            ("product", "tumor protein p73, transcript variant 1"),
        ],
    );

    assert_eq!(
        MainAreaDna::feature_copy_payload(&feature, FeatureCopyPayloadKind::Description),
        "tumor protein p73, transcript variant 1"
    );
}

#[test]
fn feature_copy_popup_text_includes_title_range_and_details() {
    let feature = make_feature(
        "mRNA",
        vec![
            ("transcript_id", "NM_005427.4"),
            ("product", "tumor protein p73, transcript variant 1"),
            ("gene", "TP73"),
        ],
    );

    let text = MainAreaDna::feature_copy_payload(&feature, FeatureCopyPayloadKind::PopupText);

    assert!(text.contains("NM_005427.4"));
    assert!(text.contains("1..10 (10 bp)"));
    assert!(text.contains("product: tumor protein p73, transcript variant 1"));
    assert!(text.contains("gene: TP73"));
}

#[test]
fn restriction_site_copy_payloads_use_tooltip_lines_and_shared_json() {
    let view = RestrictionSiteExpertView {
        seq_id: "pBR322".to_string(),
        cut_pos_1based: 410,
        paired_cut_pos_1based: 416,
        recognition_start_1based: 410,
        recognition_end_1based: 417,
        cut_index_0based: 1,
        paired_cut_index_0based: 7,
        end_geometry: "5prime_overhang".to_string(),
        number_of_cuts_for_enzyme: 1,
        selected_enzyme: Some("SgrAI".to_string()),
        enzyme_names: vec!["SgrAI".to_string()],
        recognition_iupac: Some("CRCCGGYG".to_string()),
        site_sequence: "CACCGGTG".to_string(),
        site_sequence_complement: "GTGGCCAC".to_string(),
        enzyme_cut_offset_0based: Some(1),
        overlap_bp: Some(6),
        enzyme_note: Some("Efficient cleavage requires two copies.".to_string()),
        rebase_url: Some("https://rebase.neb.com/rebase/enz/SgrAI.html".to_string()),
        tooltip_lines: vec![
            "SgrAI | 1 site | 5' overhang (6 bp)".to_string(),
            "pBR322:410..417 | cuts 410|416".to_string(),
            "5' C^ACCGGTG 3'".to_string(),
            "3' GTGGCCA^C 5'".to_string(),
        ],
        instruction: "Inspect restriction-site cleavage geometry.".to_string(),
    };

    let summary = MainAreaDna::restriction_site_copy_summary_text(&view);
    assert!(summary.contains("SgrAI | 1 site | 5' overhang (6 bp)"));
    assert!(summary.contains("5' C^ACCGGTG 3'"));

    let json = MainAreaDna::restriction_site_copy_json_text(&view);
    assert!(json.contains("\"kind\": \"restriction_site\""));
    assert!(json.contains("\"tooltip_lines\""));
    assert!(json.contains("\"selected_enzyme\": \"SgrAI\""));

    assert_eq!(
        MainAreaDna::restriction_site_context_menu_title(&view),
        "Restriction site: SgrAI at 410 bp"
    );
}

fn sanitize_for_path_for_test(s: &str) -> String {
    let mut out = String::new();
    for c in s.chars() {
        if c.is_ascii_alphanumeric() {
            out.push(c.to_ascii_lowercase());
        } else if matches!(c, ' ' | '-' | '_' | '.') && !out.ends_with('_') {
            out.push('_');
        }
    }
    let trimmed = out.trim_matches('_');
    if trimmed.is_empty() {
        "genome".to_string()
    } else {
        trimmed.to_string()
    }
}

fn restriction_ready_dna(sequence: &str) -> DNAsequence {
    let mut dna = DNAsequence::from_sequence(sequence).expect("sequence");
    *dna.restriction_enzymes_mut() = active_restriction_enzymes();
    dna.set_max_restriction_enzyme_sites(None);
    dna.update_computed_features();
    dna
}

fn transcript_derivation_test_sequence() -> DNAsequence {
    let mut bases = vec![b'A'; 40];
    for (idx, base) in [
        (8usize, b'G'),
        (9, b'T'),
        (10, b'A'),
        (11, b'G'),
        (14, b'A'),
        (15, b'G'),
        (20, b'G'),
        (21, b'T'),
        (24, b'A'),
        (25, b'G'),
    ] {
        bases[idx] = base;
    }
    let sequence = String::from_utf8(bases).expect("valid ASCII DNA");
    let mut dna = DNAsequence::from_sequence(&sequence).expect("valid DNA sequence");
    dna.features_mut().push(Feature {
        kind: "mRNA".into(),
        location: Location::Join(vec![
            Location::simple_range(2, 8),
            Location::simple_range(12, 20),
            Location::simple_range(26, 34),
        ]),
        qualifiers: vec![
            ("gene".into(), Some("GENE1".to_string())),
            ("transcript_id".into(), Some("NM_TEST_1".to_string())),
            ("label".into(), Some("NM_TEST_1".to_string())),
        ],
    });
    dna.features_mut().push(Feature {
        kind: "mRNA".into(),
        location: Location::Join(vec![
            Location::simple_range(2, 8),
            Location::simple_range(16, 20),
            Location::simple_range(26, 34),
        ]),
        qualifiers: vec![
            ("gene".into(), Some("GENE1".to_string())),
            ("transcript_id".into(), Some("NM_TEST_2".to_string())),
            ("label".into(), Some("NM_TEST_2".to_string())),
        ],
    });
    dna
}

fn write_single_line_fasta_index(path: &Path, sequence_len: usize) {
    let line = format!(
        "chr1\t{}\t6\t{}\t{}\n",
        sequence_len,
        sequence_len,
        sequence_len + 1
    );
    fs::write(path, line).expect("write fasta index");
}

fn collect_rendered_text_from_shape(shape: &egui::epaint::Shape, out: &mut Vec<String>) {
    match shape {
        egui::epaint::Shape::Text(text) => {
            let raw = text.galley.job.text.trim();
            if !raw.is_empty() {
                out.push(raw.to_string());
            }
        }
        egui::epaint::Shape::Vec(shapes) => {
            for shape in shapes {
                collect_rendered_text_from_shape(shape, out);
            }
        }
        _ => {}
    }
}

fn collect_rendered_text_rects_from_shape(
    shape: &egui::epaint::Shape,
    out: &mut Vec<(String, egui::Rect)>,
) {
    match shape {
        egui::epaint::Shape::Text(text) => {
            let raw = text.galley.job.text.trim();
            if !raw.is_empty() {
                out.push((raw.to_string(), text.visual_bounding_rect()));
            }
        }
        egui::epaint::Shape::Vec(shapes) => {
            for shape in shapes {
                collect_rendered_text_rects_from_shape(shape, out);
            }
        }
        _ => {}
    }
}

fn collect_pass_texts(ctx: &egui::Context) -> Vec<String> {
    let full_output = ctx.end_pass();
    let mut texts = Vec::new();
    for clipped in full_output.shapes {
        collect_rendered_text_from_shape(&clipped.shape, &mut texts);
    }
    texts
}

fn render_selection_formula_control_pass(
    ctx: &egui::Context,
    area: &mut MainAreaDna,
    raw_input: egui::RawInput,
) -> Vec<(String, egui::Rect)> {
    ctx.begin_pass(raw_input);
    crate::egui_compat::show_central_panel_for_test_context(
        ctx,
        egui::CentralPanel::default(),
        |ui| {
            area.render_selection_formula_inline_controls(ui, 420.0);
        },
    );
    let full_output = ctx.end_pass();
    let mut rects = Vec::new();
    for clipped in full_output.shapes {
        collect_rendered_text_rects_from_shape(&clipped.shape, &mut rects);
    }
    rects
}

fn render_feature_tree_pass(
    ctx: &egui::Context,
    area: &mut MainAreaDna,
    raw_input: egui::RawInput,
) -> Vec<(String, egui::Rect)> {
    ctx.begin_pass(raw_input);
    crate::egui_compat::show_central_panel_for_test_context(
        ctx,
        egui::CentralPanel::default(),
        |ui| {
            area.render_features(ui);
        },
    );
    let full_output = ctx.end_pass();
    let mut rects = Vec::new();
    for clipped in full_output.shapes {
        collect_rendered_text_rects_from_shape(&clipped.shape, &mut rects);
    }
    rects
}

fn center_of_rendered_text(rects: &[(String, egui::Rect)], needle: &str) -> Option<egui::Pos2> {
    rects
        .iter()
        .find_map(|(text, rect)| (text == needle).then_some(rect.center()))
}

fn make_area_with_unique_compatible_anchor(
    fallback_policy: &str,
) -> (MainAreaDna, TempDir, String, String) {
    let td = tempdir().expect("tempdir");
    let root = td.path();
    let fasta = root.join("grch38.fa");
    let gtf = root.join("grch38.gtf");
    fs::write(&fasta, ">chr1\nACGTACGTACGT\n").expect("write fasta");
    fs::write(
        &gtf,
        "chr1\tsrc\tgene\t1\t12\t.\t+\t.\tgene_id \"GENE1\"; gene_name \"GENE1\";\n",
    )
    .expect("write gtf");

    let cache_dir = root.join("cache");
    fs::create_dir_all(&cache_dir).expect("create cache dir");

    let prepared_genome_id = "Human GRCh38 Ensembl 116".to_string();
    let alias_catalog_key = "Human GRCh38 NCBI RefSeq GCF_000001405.40";
    let requested_alias = "GRCh38.p14".to_string();
    let seq_id = "alias_slice".to_string();

    // Seed one prepared install deterministically (without requiring external tools).
    let install_dir = cache_dir.join(sanitize_for_path_for_test(&prepared_genome_id));
    fs::create_dir_all(&install_dir).expect("create install dir");
    let installed_sequence = install_dir.join("sequence.fa");
    let installed_annotation = install_dir.join("annotation.gtf");
    let installed_fai = install_dir.join("sequence.fa.fai");
    fs::copy(&fasta, &installed_sequence).expect("copy sequence");
    fs::copy(&gtf, &installed_annotation).expect("copy annotation");
    write_single_line_fasta_index(&installed_fai, 12);
    let manifest = json!({
        "genome_id": prepared_genome_id,
        "sequence_source": installed_sequence.to_string_lossy(),
        "annotation_source": installed_annotation.to_string_lossy(),
        "sequence_source_type": "local",
        "annotation_source_type": "local",
        "sequence_sha1": null,
        "annotation_sha1": null,
        "sequence_path": installed_sequence.to_string_lossy(),
        "annotation_path": installed_annotation.to_string_lossy(),
        "fasta_index_path": installed_fai.to_string_lossy(),
        "gene_index_path": null,
        "blast_db_prefix": null,
        "blast_index_executable": null,
        "blast_indexed_at_unix_ms": null,
        "installed_at_unix_ms": 1
    });
    fs::write(
        install_dir.join("manifest.json"),
        serde_json::to_string_pretty(&manifest).expect("serialize manifest"),
    )
    .expect("write manifest");

    let catalog_path = root.join("catalog.json");
    let catalog_json = json!({
        prepared_genome_id.clone(): {
            "ncbi_taxonomy_id": 9606,
            "ncbi_assembly_accession": "GCF_000001405.40",
            "ncbi_assembly_name": "GRCh38",
            "sequence_local": fasta.to_string_lossy(),
            "annotations_local": gtf.to_string_lossy(),
            "cache_dir": cache_dir.to_string_lossy()
        },
        alias_catalog_key: {
            "ncbi_taxonomy_id": 9606,
            "ncbi_assembly_accession": "GCF_000001405.40",
            "ncbi_assembly_name": "GRCh38.p14",
            "sequence_local": fasta.to_string_lossy(),
            "annotations_local": gtf.to_string_lossy(),
            "cache_dir": cache_dir.to_string_lossy()
        }
    });
    fs::write(
        &catalog_path,
        serde_json::to_string_pretty(&catalog_json).expect("serialize catalog"),
    )
    .expect("write catalog");
    let catalog_path_str = catalog_path.to_string_lossy().to_string();

    let mut engine = GentleEngine::new();
    engine
        .apply(Operation::ExtractGenomeRegion {
            genome_id: requested_alias.clone(),
            chromosome: "chr1".to_string(),
            start_1based: 3,
            end_1based: 10,
            output_id: Some(seq_id.clone()),
            annotation_scope: None,
            max_annotation_features: None,
            include_genomic_annotation: None,
            catalog_path: Some(catalog_path_str),
            cache_dir: None,
        })
        .expect("extract region through compatibility fallback");
    engine
        .apply(Operation::SetParameter {
            name: "genome_anchor_prepared_fallback_policy".to_string(),
            value: json!(fallback_policy),
        })
        .expect("set fallback policy");
    let dna = engine
        .state()
        .sequences
        .get(&seq_id)
        .cloned()
        .expect("sequence created by extract");
    let engine = Arc::new(RwLock::new(engine));
    (
        MainAreaDna::new(dna, Some(seq_id.clone()), Some(engine)),
        td,
        seq_id,
        prepared_genome_id,
    )
}

#[test]
fn feature_tree_count_label_shows_visible_over_total_for_linear_view() {
    assert_eq!(
        MainAreaDna::format_feature_tree_count_label("mRNA", 3, 7, true),
        "mRNA (3/7)"
    );
}

#[test]
fn feature_tree_count_label_uses_total_when_all_visible() {
    assert_eq!(
        MainAreaDna::format_feature_tree_count_label("mRNA", 4, 4, true),
        "mRNA (4)"
    );
}

#[test]
fn feature_tree_count_label_uses_total_for_non_linear_scope() {
    assert_eq!(
        MainAreaDna::format_feature_tree_count_label("mRNA", 3, 7, false),
        "mRNA (7)"
    );
}

#[test]
fn feature_tree_filter_terms_parse_scoped_and_free_terms() {
    assert_eq!(
        MainAreaDna::parse_feature_tree_filter_terms("kind:mrna TP73 range:1..10"),
        vec![
            (Some("kind".to_string()), "mrna".to_string()),
            (None, "tp73".to_string()),
            (Some("range".to_string()), "1..10".to_string()),
        ]
    );
}

#[test]
fn feature_tree_filter_term_presence_is_case_insensitive() {
    assert!(MainAreaDna::filter_term_present(
        "kind:tracks source:bed",
        "SOURCE:BED"
    ));
    assert!(!MainAreaDna::filter_term_present(
        "kind:tracks source:bed",
        "source:vcf"
    ));
}

#[test]
fn feature_tree_filter_term_add_remove_and_toggle_roundtrip() {
    let mut filter = "kind:tracks source:bed".to_string();
    MainAreaDna::set_filter_term_enabled(&mut filter, "source:vcf", true);
    assert!(MainAreaDna::filter_term_present(&filter, "source:vcf"));

    MainAreaDna::set_filter_term_enabled(&mut filter, "SOURCE:BED", false);
    assert!(!MainAreaDna::filter_term_present(&filter, "source:bed"));
    assert!(MainAreaDna::filter_term_present(&filter, "kind:tracks"));
    assert!(MainAreaDna::filter_term_present(&filter, "source:vcf"));
}

#[test]
fn feature_tree_filter_matches_array_scopes() {
    let feature = make_feature(
        "track",
        vec![
            ("gentle_track_source", "Array"),
            ("gentle_array_platform", "Clariom D human"),
            ("gentle_array_contrast", "AdTAp73alpha-AdGFP"),
            ("gene", "TP73"),
            ("feature_id", "PSR0001"),
        ],
    );
    assert!(MainAreaDna::feature_tree_matches_filter(
        &feature,
        "source:array track:Clariom contrast:AdTAp73alpha-AdGFP gene:TP73",
        "array",
        "AdTAp73alpha-AdGFP PSR0001",
        "10..20",
    ));
    assert!(!MainAreaDna::feature_tree_matches_filter(
        &feature,
        "contrast:missing",
        "array",
        "AdTAp73alpha-AdGFP PSR0001",
        "10..20",
    ));
}

#[test]
fn feature_tree_filter_matches_focus_terms_for_dense_layers() {
    let regulatory = make_feature(
        "regulatory",
        vec![
            ("regulatory_class", "enhancer"),
            ("note", "H3K27ac active region"),
        ],
    );
    assert!(MainAreaDna::feature_tree_matches_filter(
        &regulatory,
        "regulatory",
        "regulatory",
        "enhancer: H3K27ac active region",
        "1..10",
    ));
    assert!(MainAreaDna::feature_tree_matches_filter(
        &regulatory,
        "regulatory:enhancer",
        "regulatory",
        "enhancer: H3K27ac active region",
        "1..10",
    ));

    let repeat = make_feature(
        "repeat_region",
        vec![
            ("repName", "L1PA2"),
            ("repClass", "LINE"),
            ("repFamily", "L1"),
        ],
    );
    assert!(MainAreaDna::feature_tree_matches_filter(
        &repeat,
        "repeat",
        "repeat_region",
        "L1PA2 (LINE / L1)",
        "1..10",
    ));

    let track = make_feature(
        "track",
        vec![
            ("gentle_generated", "genome_bed_track"),
            ("gentle_track_source", "BED"),
            ("gentle_track_name", "H3K27ac"),
            ("gentle_track_file", "/tmp/peaks.bed"),
            ("label", "peak 1"),
        ],
    );
    assert!(MainAreaDna::feature_tree_matches_filter(
        &track, "track", "tracks", "peak 1", "1..10",
    ));
}

#[test]
fn feature_tree_focus_preset_sets_filter_and_grouping() {
    let dna = DNAsequence::from_sequence("AAAAAAAAAA").expect("sequence");
    let mut area = MainAreaDna::new(dna, Some("feature_tree_focus".to_string()), None);

    area.apply_feature_tree_focus_preset(super::FeatureTreeFocusPreset::Regulatory);
    assert_eq!(area.feature_tree_filter, "regulatory");
    assert_eq!(
        area.feature_tree_grouping_mode,
        super::FeatureTreeGroupingMode::Always
    );

    area.apply_feature_tree_focus_preset(super::FeatureTreeFocusPreset::Cloning);
    assert!(area.feature_tree_filter.is_empty());
    assert_eq!(
        area.feature_tree_grouping_mode,
        super::FeatureTreeGroupingMode::Auto
    );
}

#[test]
fn parse_positive_usize_text_accepts_positive_integer() {
    assert_eq!(
        MainAreaDna::parse_positive_usize_text("250", "extension length").unwrap(),
        250
    );
}

#[test]
fn build_sequencing_confirmation_targets_includes_full_span_and_junctions() {
    let dna = DNAsequence::from_sequence(&"ACGT".repeat(50)).expect("sequence");
    let mut area = MainAreaDna::new(dna, Some("seq1".to_string()), None);
    area.sequencing_confirmation_ui.include_full_span_target = true;
    area.sequencing_confirmation_ui.junction_positions_0based = "25,60".to_string();
    area.sequencing_confirmation_ui.junction_flank_bp = "12".to_string();

    let targets = area
        .build_sequencing_confirmation_targets()
        .expect("targets should build");

    assert_eq!(targets.len(), 3);
    assert_eq!(targets[0].kind, SequencingConfirmationTargetKind::FullSpan);
    assert_eq!(targets[0].start_0based, 0);
    assert_eq!(targets[0].end_0based_exclusive, 200);
    assert_eq!(targets[1].kind, SequencingConfirmationTargetKind::Junction);
    assert_eq!(targets[1].junction_left_end_0based, Some(25));
    assert_eq!(targets[1].start_0based, 13);
    assert_eq!(targets[1].end_0based_exclusive, 37);
    assert_eq!(targets[2].junction_left_end_0based, Some(60));
}

#[test]
fn build_confirm_construct_reads_operation_uses_gui_state() {
    let dna = DNAsequence::from_sequence(&"ACGT".repeat(40)).expect("sequence");
    let mut area = MainAreaDna::new(dna, Some("expected_seq".to_string()), None);
    area.sequencing_confirmation_ui.baseline_seq_id_text = "baseline_seq".to_string();
    area.sequencing_confirmation_ui.read_seq_ids_text = "read_a,read_b".to_string();
    area.sequencing_confirmation_ui.trace_ids_text = "trace_a,trace_b".to_string();
    area.sequencing_confirmation_ui.report_id = String::new();
    area.sequencing_confirmation_ui.include_full_span_target = false;
    area.sequencing_confirmation_ui.junction_positions_0based = "40".to_string();
    area.sequencing_confirmation_ui.junction_flank_bp = "10".to_string();
    area.sequencing_confirmation_ui.allow_reverse_complement = false;
    area.sequencing_confirmation_ui.alignment_mode = PairwiseAlignmentMode::Global;
    area.sequencing_confirmation_ui.match_score = "3".to_string();
    area.sequencing_confirmation_ui.mismatch_score = "-4".to_string();
    area.sequencing_confirmation_ui.gap_open = "-6".to_string();
    area.sequencing_confirmation_ui.gap_extend = "-2".to_string();
    area.sequencing_confirmation_ui.min_identity_fraction = "0.90".to_string();
    area.sequencing_confirmation_ui.min_target_coverage_fraction = "0.75".to_string();

    let (report_id, op) = area
        .build_confirm_construct_reads_operation()
        .expect("operation should build");

    assert_eq!(report_id, "expected_seq_seq_confirm");
    match op {
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
            assert_eq!(expected_seq_id, "expected_seq");
            assert_eq!(baseline_seq_id, Some("baseline_seq".to_string()));
            assert_eq!(
                read_seq_ids,
                vec!["read_a".to_string(), "read_b".to_string()]
            );
            assert_eq!(
                trace_ids,
                vec!["trace_a".to_string(), "trace_b".to_string()]
            );
            assert_eq!(targets.len(), 1);
            assert_eq!(targets[0].junction_left_end_0based, Some(40));
            assert_eq!(alignment_mode, PairwiseAlignmentMode::Global);
            assert_eq!(match_score, 3);
            assert_eq!(mismatch_score, -4);
            assert_eq!(gap_open, -6);
            assert_eq!(gap_extend, -2);
            assert_eq!(min_identity_fraction, 0.90);
            assert_eq!(min_target_coverage_fraction, 0.75);
            assert!(!allow_reverse_complement);
            assert_eq!(report_id, Some("expected_seq_seq_confirm".to_string()));
        }
        other => panic!("unexpected operation: {other:?}"),
    }
}

#[test]
fn build_confirm_construct_reads_operation_allows_trace_only_gui_state() {
    let dna = DNAsequence::from_sequence(&"ACGT".repeat(24)).expect("sequence");
    let mut area = MainAreaDna::new(dna, Some("expected_seq".to_string()), None);
    area.sequencing_confirmation_ui.read_seq_ids_text.clear();
    area.sequencing_confirmation_ui.trace_ids_text = "trace_only".to_string();
    area.sequencing_confirmation_ui.include_full_span_target = true;

    let (_report_id, op) = area
        .build_confirm_construct_reads_operation()
        .expect("trace-only GUI operation should build");

    match op {
        Operation::ConfirmConstructReads {
            baseline_seq_id,
            read_seq_ids,
            trace_ids,
            ..
        } => {
            assert_eq!(baseline_seq_id, None);
            assert!(read_seq_ids.is_empty());
            assert_eq!(trace_ids, vec!["trace_only".to_string()]);
        }
        other => panic!("unexpected operation: {other:?}"),
    }
}

#[test]
fn build_import_sequencing_trace_operation_uses_gui_state() {
    let dna = DNAsequence::from_sequence(&"ACGT".repeat(24)).expect("sequence");
    let mut area = MainAreaDna::new(dna, Some("expected_seq".to_string()), None);
    area.sequencing_confirmation_ui.trace_import_path = "/tmp/demo_trace.ab1".to_string();
    area.sequencing_confirmation_ui.trace_import_id = "trace_demo".to_string();
    area.sequencing_confirmation_ui
        .trace_import_associate_with_expected_seq = true;

    let op = area
        .build_import_sequencing_trace_operation()
        .expect("trace import operation should build");

    match op {
        Operation::ImportSequencingTrace {
            path,
            trace_id,
            seq_id,
        } => {
            assert_eq!(path, "/tmp/demo_trace.ab1");
            assert_eq!(trace_id, Some("trace_demo".to_string()));
            assert_eq!(seq_id, Some("expected_seq".to_string()));
        }
        other => panic!("unexpected operation: {other:?}"),
    }
}

#[test]
fn handle_imported_sequencing_trace_result_selects_trace_and_appends_to_run() {
    let dna = DNAsequence::from_sequence(&"ACGT".repeat(24)).expect("sequence");
    let mut area = MainAreaDna::new(dna, Some("expected_seq".to_string()), None);
    area.sequencing_confirmation_ui.trace_ids_text = "trace_a".to_string();
    area.sequencing_confirmation_ui.trace_import_add_to_run = true;

    area.handle_imported_sequencing_trace_result(&OpResult {
        op_id: "op-import-trace".to_string(),
        created_seq_ids: vec![],
        changed_seq_ids: vec![],
        warnings: vec![],
        messages: vec![],
        protocol_cartoon_preview: None,
        genome_annotation_projection: None,
        sequence_alignment: None,
        protein_derivation_report: None,
        reverse_translation_report: None,
        protease_digest_report: None,
        protein_residue_genomic_coordinates: None,
        exon_skip_selection_plan: None,
        exon_skip_materialization: None,
        cdna_assay_test_report: None,
        cdna_assay_product_materialization: None,
        transcript_qpcr_panel: None,
        primer_specificity_report: None,
        construct_reasoning_graph: None,
        sequencing_confirmation_report: None,
        sequencing_trace_import_report: Some(SequencingTraceImportReport {
            schema: "gentle.sequencing_trace_import_report.v2".to_string(),
            trace_id: "trace_b".to_string(),
            format: SequencingTraceFormat::AbiAb1,
            source_path: "/tmp/demo_trace.ab1".to_string(),
            imported_at_unix_ms: 1,
            seq_id: Some("expected_seq".to_string()),
            sample_name: None,
            run_name: None,
            called_base_count: 4,
            confidence_value_count: 4,
            peak_location_count: 4,
            has_curve_data: true,
            channel_count: 4,
            warnings: vec![],
        }),
        sequencing_trace_record: None,
        sequencing_trace_summaries: None,
        sequencing_primer_overlay_report: None,
        cutrun_dataset_list: None,
        cutrun_dataset_status: None,
        cutrun_dataset_projection: None,
        cutrun_read_report: None,
        cutrun_read_report_summaries: None,
        cutrun_read_coverage_export: None,
        cutrun_regulatory_support: None,
        gene_set_resolution: None,
        gene_set_promoter_cohort: None,
        gene_set_cutrun_regulatory_support: None,
        read_acquisition_report: None,
        microarray_projection: None,
        probe_region_evidence_interpretation: None,
        genome_coordinate_projection: None,
        rna_read_gene_support_summary: None,
        rna_read_gene_support_audit: None,
        rna_read_target_quality_export: None,
        rna_read_batch_map_report: None,
        rna_read_isoform_preflight: None,
        tfbs_region_summary: None,
        tfbs_score_tracks: None,
        tfbs_track_similarity: None,
        multi_gene_promoter_tfbs: None,
        promoter_cohort_comparison: None,
        repeat_annotation_query: None,
        sequence_repeat_overlaps: None,
        repeat_feature_materialization: None,
        repeat_environment_cohort: None,
        window_cohort_tfbs: None,
        tfbs_hit_scan: None,
        restriction_site_scan: None,
        jaspar_remote_metadata_snapshot: None,
        jaspar_catalog_report: None,
        tf_query_resolution_report: None,
        jaspar_entry_expert_view: None,
        jaspar_registry_benchmark: None,
        jaspar_entry_presentation: None,
        sequence_context_view: None,
        sequence_context_bundle: None,
        alternative_promoter_comparison: None,
        variant_promoter_context: None,
        promoter_evidence_matrix: None,
        isoform_promoter_comparison: None,
        promoter_expression_evidence: None,
        promoter_artifact_manifest: None,
        promoter_reporter_candidates: None,
        reporter_catalog: None,
        reporter_recommendation: None,
        reporter_corpus_export: None,
        reporter_construct_handoff: None,
        uniprot_projection_audit: None,
        uniprot_projection_audit_parity: None,
        lab_assistant_instructions: None,
    });

    assert_eq!(area.sequencing_confirmation_ui.selected_trace_id, "trace_b");
    assert_eq!(
        area.sequencing_confirmation_ui.trace_ids_text,
        "trace_a,trace_b"
    );
}

#[test]
fn build_suggest_sequencing_primers_operation_uses_gui_state() {
    let dna = DNAsequence::from_sequence(&"ACGT".repeat(24)).expect("sequence");
    let mut area = MainAreaDna::new(dna, Some("expected_seq".to_string()), None);
    area.sequencing_confirmation_ui.primer_seq_ids_text = "primer_a,primer_b".to_string();
    area.sequencing_confirmation_ui.selected_report_id = "confirm_report".to_string();
    area.sequencing_confirmation_ui.primer_min_3prime_anneal_bp = "20".to_string();
    area.sequencing_confirmation_ui
        .primer_predicted_read_length_bp = "650".to_string();

    let op = area
        .build_suggest_sequencing_primers_operation()
        .expect("sequencing-primer overlay operation should build");

    match op {
        Operation::SuggestSequencingPrimers {
            expected_seq_id,
            primer_seq_ids,
            confirmation_report_id,
            min_3prime_anneal_bp,
            predicted_read_length_bp,
        } => {
            assert_eq!(expected_seq_id, "expected_seq");
            assert_eq!(
                primer_seq_ids,
                vec!["primer_a".to_string(), "primer_b".to_string()]
            );
            assert_eq!(confirmation_report_id, Some("confirm_report".to_string()));
            assert_eq!(min_3prime_anneal_bp, 20);
            assert_eq!(predicted_read_length_bp, 650);
        }
        other => panic!("unexpected operation: {other:?}"),
    }
}

#[test]
fn build_suggest_sequencing_primers_operation_allows_report_only_gui_state() {
    let dna = DNAsequence::from_sequence(&"ACGT".repeat(24)).expect("sequence");
    let mut area = MainAreaDna::new(dna, Some("expected_seq".to_string()), None);
    area.sequencing_confirmation_ui.primer_seq_ids_text.clear();
    area.sequencing_confirmation_ui.selected_report_id = "confirm_report".to_string();

    let op = area
        .build_suggest_sequencing_primers_operation()
        .expect("report-only sequencing-primer overlay operation should build");

    match op {
        Operation::SuggestSequencingPrimers {
            expected_seq_id,
            primer_seq_ids,
            confirmation_report_id,
            ..
        } => {
            assert_eq!(expected_seq_id, "expected_seq");
            assert!(primer_seq_ids.is_empty());
            assert_eq!(confirmation_report_id, Some("confirm_report".to_string()));
        }
        other => panic!("unexpected operation: {other:?}"),
    }
}

#[test]
fn sequencing_confirmation_unresolved_summary_markdown_uses_matching_primer_overlay() {
    let dna = DNAsequence::from_sequence(&"ACGT".repeat(24)).expect("sequence");
    let mut area = MainAreaDna::new(dna, Some("expected_seq".to_string()), None);
    let report = SequencingConfirmationReport {
        report_id: "confirm_report".to_string(),
        expected_seq_id: "expected_seq".to_string(),
        overall_status: SequencingConfirmationStatus::InsufficientEvidence,
        targets: vec![SequencingConfirmationTargetResult {
            target_id: "junction_1".to_string(),
            label: "Insert junction".to_string(),
            kind: SequencingConfirmationTargetKind::Junction,
            start_0based: 8,
            end_0based_exclusive: 12,
            required: true,
            status: SequencingConfirmationStatus::InsufficientEvidence,
            covered_bp: 0,
            target_length_bp: 4,
            reason: "Coverage stops before the junction.".to_string(),
            ..Default::default()
        }],
        ..Default::default()
    };
    let mut state = ProjectState::default();
    state.metadata.insert(
        crate::engine::SEQUENCING_CONFIRMATION_REPORTS_METADATA_KEY.to_string(),
        serde_json::json!({
            "schema": "gentle.sequencing_confirmation_reports.v1",
            "updated_at_unix_ms": 0,
            "reports": {
                "confirm_report": serde_json::to_value(&report)
                    .expect("serialize confirmation report")
            }
        }),
    );
    area.engine = Some(Arc::new(RwLock::new(GentleEngine::from_state(state))));
    area.sequencing_confirmation_ui.primer_overlay_report = Some(SequencingPrimerOverlayReport {
        schema: "gentle.sequencing_primer_overlay_report.v1".to_string(),
        expected_seq_id: "expected_seq".to_string(),
        confirmation_report_id: Some("confirm_report".to_string()),
        min_3prime_anneal_bp: 20,
        predicted_read_length_bp: 650,
        primer_seq_ids: vec!["primer_a".to_string()],
        problem_guidance_count: 1,
        problem_guidance: vec![crate::engine::SequencingPrimerProblemGuidanceRow {
            problem_id: "junction_1".to_string(),
            problem_kind: SequencingPrimerProblemKind::Target,
            problem_label: "Insert junction".to_string(),
            problem_summary: "No retained evidence row spans the junction yet.".to_string(),
            recommended_primer_seq_id: Some("primer_a".to_string()),
            recommended_orientation: Some(SequencingPrimerOrientation::ForwardRead),
            recommended_three_prime_distance_bp: Some(42),
            candidate_count: 1,
            reason: "Closest existing sequencing primer covers the unresolved junction."
                .to_string(),
            ..Default::default()
        }],
        proposal_count: 1,
        proposals: vec![SequencingPrimerProposalRow {
            proposal_id: "proposal_a".to_string(),
            problem_id: "junction_1".to_string(),
            problem_kind: SequencingPrimerProblemKind::Target,
            problem_label: "Insert junction".to_string(),
            problem_summary: "No retained evidence row spans the junction yet.".to_string(),
            orientation: SequencingPrimerOrientation::ForwardRead,
            primer_sequence: "ACGTACGTACGTACGTACGT".to_string(),
            anneal_sequence: "ACGTACGTACGTACGTACGT".to_string(),
            anneal_start_0based: 0,
            anneal_end_0based_exclusive: 20,
            three_prime_position_0based: 19,
            predicted_read_span_start_0based: 0,
            predicted_read_span_end_0based_exclusive: 120,
            three_prime_distance_bp: 19,
            tm_c: 61.5,
            gc_fraction: 0.5,
            anneal_hits: 1,
            reason: "Fresh primer would cover the unresolved junction directly.".to_string(),
            ..Default::default()
        }],
        ..Default::default()
    });

    let (_report, text) = area
        .build_sequencing_confirmation_unresolved_summary_markdown("confirm_report")
        .expect("build unresolved summary markdown");

    assert!(text.contains("## Primer Guidance"));
    assert!(text.contains("primer=`primer_a`"));
    assert!(text.contains("Fresh Primer Proposals"));
    assert!(text.contains("ACGTACGTACGTACGTACGT"));
}

#[test]
fn sequencing_trace_curve_unavailable_message_guides_reimport_for_legacy_trace() {
    let trace = SequencingTraceRecord {
        schema: "gentle.sequencing_trace_record.v1".to_string(),
        trace_id: "legacy_trace".to_string(),
        format: SequencingTraceFormat::AbiAb1,
        source_path: "legacy.ab1".to_string(),
        imported_at_unix_ms: 1,
        seq_id: None,
        sample_name: None,
        sample_well: None,
        run_name: None,
        machine_name: None,
        machine_model: None,
        called_bases: "ACGT".to_string(),
        called_base_confidence_values: vec![40, 40, 40, 40],
        peak_locations: vec![10, 20, 30, 40],
        channel_data: vec![],
        channel_summaries: vec![],
        clip_start_base_index: None,
        clip_end_base_index_exclusive: None,
        comments_text: None,
    };

    assert_eq!(
        MainAreaDna::sequencing_trace_curve_unavailable_message(&trace),
        Some("curve data unavailable; re-import this trace to inspect chromatogram curves")
    );
}

#[test]
fn sequencing_trace_base_index_for_variant_uses_nearest_peak() {
    let trace = SequencingTraceRecord {
        schema: "gentle.sequencing_trace_record.v2".to_string(),
        trace_id: "trace_demo".to_string(),
        format: SequencingTraceFormat::AbiAb1,
        source_path: "trace.ab1".to_string(),
        imported_at_unix_ms: 1,
        seq_id: None,
        sample_name: None,
        sample_well: None,
        run_name: None,
        machine_name: None,
        machine_model: None,
        called_bases: "ACGT".to_string(),
        called_base_confidence_values: vec![40, 41, 42, 43],
        peak_locations: vec![10, 20, 30, 40],
        channel_data: vec![SequencingTraceChannelData {
            channel: "A".to_string(),
            trace_set: "DATA9".to_string(),
            points: vec![0; 64],
        }],
        channel_summaries: vec![],
        clip_start_base_index: Some(1),
        clip_end_base_index_exclusive: Some(4),
        comments_text: None,
    };
    let variant = SequencingConfirmationVariantRow {
        variant_id: "var-1".to_string(),
        label: "edit".to_string(),
        peak_center: Some(27),
        ..Default::default()
    };

    assert_eq!(
        MainAreaDna::sequencing_trace_base_index_for_variant(&trace, &variant),
        Some(2)
    );
}

#[test]
fn sequencing_trace_sample_window_uses_trace_browser_base_index() {
    let trace = SequencingTraceRecord {
        schema: "gentle.sequencing_trace_record.v2".to_string(),
        trace_id: "trace_demo".to_string(),
        format: SequencingTraceFormat::AbiAb1,
        source_path: "trace.ab1".to_string(),
        imported_at_unix_ms: 1,
        seq_id: None,
        sample_name: None,
        sample_well: None,
        run_name: None,
        machine_name: None,
        machine_model: None,
        called_bases: "ACGTN".to_string(),
        called_base_confidence_values: vec![35, 36, 37, 38, 39],
        peak_locations: vec![10, 20, 30, 40, 50],
        channel_data: vec![SequencingTraceChannelData {
            channel: "A".to_string(),
            trace_set: "DATA9".to_string(),
            points: vec![0; 80],
        }],
        channel_summaries: vec![],
        clip_start_base_index: None,
        clip_end_base_index_exclusive: None,
        comments_text: None,
    };

    assert_eq!(
        MainAreaDna::sequencing_trace_sample_window(&trace, None, Some(3), 1),
        Some((12, 68))
    );
}

#[test]
fn splicing_boundary_motif_rows_pair_donor_and_acceptor_markers() {
    let view = SplicingExpertView {
        seq_id: "seq".to_string(),
        target_feature_id: 1,
        scope: SplicingScopePreset::AllOverlappingAnyStrand,
        group_label: "GENE1".to_string(),
        strand: "+".to_string(),
        region_start_1based: 1,
        region_end_1based: 40,
        transcript_count: 1,
        unique_exon_count: 2,
        instruction: String::new(),
        transcripts: vec![SplicingTranscriptLane {
            transcript_feature_id: 11,
            transcript_id: "tx1".to_string(),
            label: "tx1".to_string(),
            strand: "+".to_string(),
            exons: vec![],
            exon_cds_phases: vec![],
            introns: vec![],
            has_target_feature: true,
        }],
        unique_exons: vec![],
        matrix_rows: vec![],
        boundaries: vec![
            SplicingBoundaryMarker {
                transcript_feature_id: 11,
                transcript_id: "tx1".to_string(),
                side: "donor".to_string(),
                position_1based: 10,
                motif_2bp: "GC".to_string(),
                canonical: false,
                canonical_pair: false,
                partner_position_1based: 20,
                paired_motif_signature: "GC-AG".to_string(),
                motif_class: "gc_ag_major_noncanonical".to_string(),
                annotation:
                    "GC-AG intron; a known non-canonical major/U2-type splice-site motif class."
                        .to_string(),
            },
            SplicingBoundaryMarker {
                transcript_feature_id: 11,
                transcript_id: "tx1".to_string(),
                side: "acceptor".to_string(),
                position_1based: 20,
                motif_2bp: "AG".to_string(),
                canonical: true,
                canonical_pair: false,
                partner_position_1based: 10,
                paired_motif_signature: "GC-AG".to_string(),
                motif_class: "gc_ag_major_noncanonical".to_string(),
                annotation:
                    "GC-AG intron; a known non-canonical major/U2-type splice-site motif class."
                        .to_string(),
            },
        ],
        intron_signals: vec![],
        junctions: vec![],
        events: vec![],
    };

    let rows = MainAreaDna::splicing_boundary_motif_rows(&view);
    assert_eq!(rows.len(), 1);
    assert_eq!(rows[0].transcript_feature_id, 11);
    assert_eq!(rows[0].donor_position_1based, 10);
    assert_eq!(rows[0].donor_motif_2bp, "GC");
    assert_eq!(rows[0].acceptor_position_1based, 20);
    assert_eq!(rows[0].acceptor_motif_2bp, "AG");
    assert_eq!(rows[0].paired_motif_signature, "GC-AG");
    assert_eq!(rows[0].motif_class, "gc_ag_major_noncanonical");
}

#[test]
fn splicing_intron_signal_rows_preserve_branchpoint_and_polyy_metadata() {
    let view = SplicingExpertView {
        seq_id: "seq".to_string(),
        target_feature_id: 1,
        scope: SplicingScopePreset::AllOverlappingAnyStrand,
        group_label: "GENE1".to_string(),
        strand: "+".to_string(),
        region_start_1based: 1,
        region_end_1based: 80,
        transcript_count: 1,
        unique_exon_count: 2,
        instruction: String::new(),
        transcripts: vec![],
        unique_exons: vec![],
        matrix_rows: vec![],
        boundaries: vec![],
        intron_signals: vec![SplicingIntronSignal {
            transcript_feature_id: 11,
            transcript_id: "tx1".to_string(),
            donor_position_1based: 10,
            acceptor_position_1based: 40,
            intron_length_bp: 31,
            branchpoint_position_1based: Some(18),
            branchpoint_motif: "CTAAC".to_string(),
            branchpoint_score: 4.1,
            branchpoint_annotation:
                "Strong branchpoint-like adenine in the usual 18-40 nt acceptor-proximal window (heuristic, not a splice predictor)."
                    .to_string(),
            polypyrimidine_start_1based: Some(28),
            polypyrimidine_end_1based: Some(37),
            polypyrimidine_fraction: 0.9,
            polypyrimidine_annotation:
                "Strong acceptor-proximal polypyrimidine tract by simple pyrimidine-density heuristic."
                    .to_string(),
        }],
        junctions: vec![],
        events: vec![],
    };

    let rows = MainAreaDna::splicing_intron_signal_rows(&view);
    assert_eq!(rows.len(), 1);
    assert_eq!(rows[0].branchpoint_position_1based, Some(18));
    assert_eq!(rows[0].branchpoint_motif, "CTAAC");
    assert_eq!(rows[0].polypyrimidine_start_1based, Some(28));
    assert_eq!(rows[0].polypyrimidine_end_1based, Some(37));
    assert_eq!(rows[0].polypyrimidine_fraction_percent, 90);
}

#[test]
fn splicing_selected_intron_signal_key_falls_back_and_preserves_valid_selection() {
    let dna = DNAsequence::from_sequence("ACGTACGTACGT").expect("dna");
    let mut area = MainAreaDna::new(dna, Some("seq".to_string()), None);
    let view = SplicingExpertView {
        seq_id: "seq".to_string(),
        target_feature_id: 1,
        scope: SplicingScopePreset::AllOverlappingAnyStrand,
        group_label: "GENE1".to_string(),
        strand: "+".to_string(),
        region_start_1based: 1,
        region_end_1based: 80,
        transcript_count: 1,
        unique_exon_count: 2,
        instruction: String::new(),
        transcripts: vec![],
        unique_exons: vec![],
        matrix_rows: vec![],
        boundaries: vec![],
        intron_signals: vec![
            SplicingIntronSignal {
                transcript_feature_id: 11,
                transcript_id: "tx1".to_string(),
                donor_position_1based: 10,
                acceptor_position_1based: 40,
                intron_length_bp: 31,
                branchpoint_position_1based: Some(18),
                branchpoint_motif: "CTAAC".to_string(),
                branchpoint_score: 4.1,
                branchpoint_annotation: "branchpoint note".to_string(),
                polypyrimidine_start_1based: Some(28),
                polypyrimidine_end_1based: Some(37),
                polypyrimidine_fraction: 0.9,
                polypyrimidine_annotation: "polyY note".to_string(),
            },
            SplicingIntronSignal {
                transcript_feature_id: 12,
                transcript_id: "tx2".to_string(),
                donor_position_1based: 20,
                acceptor_position_1based: 60,
                intron_length_bp: 41,
                branchpoint_position_1based: None,
                branchpoint_motif: String::new(),
                branchpoint_score: 0.0,
                branchpoint_annotation: "none".to_string(),
                polypyrimidine_start_1based: None,
                polypyrimidine_end_1based: None,
                polypyrimidine_fraction: 0.0,
                polypyrimidine_annotation: "none".to_string(),
            },
        ],
        junctions: vec![],
        events: vec![],
    };

    let fallback = area
        .splicing_selected_intron_signal_key(&view)
        .expect("fallback signal");
    assert_eq!(fallback.transcript_feature_id, 11);
    assert_eq!(fallback.donor_position_1based, 10);
    assert_eq!(fallback.acceptor_position_1based, 40);

    area.splicing_expert_selected_intron_signal_key = Some(SplicingIntronSignalKey {
        transcript_feature_id: 12,
        donor_position_1based: 20,
        acceptor_position_1based: 60,
    });
    let preserved = area
        .splicing_selected_intron_signal_key(&view)
        .expect("preserved selection");
    assert_eq!(preserved.transcript_feature_id, 12);
    assert_eq!(preserved.donor_position_1based, 20);
    assert_eq!(preserved.acceptor_position_1based, 60);
}

#[test]
fn splicing_selected_intron_signal_row_recovers_when_previous_selection_is_missing() {
    let dna = DNAsequence::from_sequence("ACGTACGTACGT").expect("dna");
    let mut area = MainAreaDna::new(dna, Some("seq".to_string()), None);
    area.splicing_expert_selected_intron_signal_key = Some(SplicingIntronSignalKey {
        transcript_feature_id: 99,
        donor_position_1based: 1,
        acceptor_position_1based: 2,
    });
    let view = SplicingExpertView {
        seq_id: "seq".to_string(),
        target_feature_id: 1,
        scope: SplicingScopePreset::AllOverlappingAnyStrand,
        group_label: "GENE1".to_string(),
        strand: "+".to_string(),
        region_start_1based: 1,
        region_end_1based: 80,
        transcript_count: 1,
        unique_exon_count: 2,
        instruction: String::new(),
        transcripts: vec![],
        unique_exons: vec![],
        matrix_rows: vec![],
        boundaries: vec![],
        intron_signals: vec![SplicingIntronSignal {
            transcript_feature_id: 11,
            transcript_id: "tx1".to_string(),
            donor_position_1based: 10,
            acceptor_position_1based: 40,
            intron_length_bp: 31,
            branchpoint_position_1based: Some(18),
            branchpoint_motif: "CTAAC".to_string(),
            branchpoint_score: 4.1,
            branchpoint_annotation: "branchpoint note".to_string(),
            polypyrimidine_start_1based: Some(28),
            polypyrimidine_end_1based: Some(37),
            polypyrimidine_fraction: 0.9,
            polypyrimidine_annotation: "polyY note".to_string(),
        }],
        junctions: vec![],
        events: vec![],
    };

    let selected = area
        .splicing_selected_intron_signal_row(&view)
        .expect("selected row");
    assert_eq!(selected.transcript_feature_id, 11);
    assert_eq!(selected.branchpoint_motif, "CTAAC");
    assert_eq!(
        area.splicing_expert_selected_intron_signal_key,
        Some(SplicingIntronSignalKey {
            transcript_feature_id: 11,
            donor_position_1based: 10,
            acceptor_position_1based: 40,
        })
    );
}

#[test]
fn splicing_intron_signal_description_reports_detected_sites_and_fallbacks() {
    let detected = SplicingIntronSignalRow {
        transcript_feature_id: 11,
        transcript_id: "tx1".to_string(),
        donor_position_1based: 10,
        acceptor_position_1based: 40,
        intron_length_bp: 31,
        branchpoint_position_1based: Some(18),
        branchpoint_motif: "CTAAC".to_string(),
        branchpoint_annotation: "branchpoint note".to_string(),
        polypyrimidine_start_1based: Some(28),
        polypyrimidine_end_1based: Some(37),
        polypyrimidine_fraction_percent: 90,
        polypyrimidine_annotation: "polyY note".to_string(),
    };
    assert_eq!(
        MainAreaDna::splicing_intron_signal_branchpoint_summary(&detected),
        "18 (CTAAC)"
    );
    assert_eq!(
        MainAreaDna::splicing_intron_signal_polypyrimidine_summary(&detected),
        "28..37 (90% pyrimidines)"
    );
    let detected_description = MainAreaDna::splicing_intron_signal_description_text(&detected);
    assert!(detected_description.contains("31 bp intron is now highlighted"));
    assert!(detected_description.contains("18 (CTAAC)"));
    assert!(detected_description.contains("28..37 (90% pyrimidines)"));

    let missing = SplicingIntronSignalRow {
        branchpoint_position_1based: None,
        polypyrimidine_start_1based: None,
        polypyrimidine_end_1based: None,
        ..detected
    };
    assert_eq!(
        MainAreaDna::splicing_intron_signal_branchpoint_summary(&missing),
        "not detected above the current heuristic"
    );
    assert_eq!(
        MainAreaDna::splicing_intron_signal_polypyrimidine_summary(&missing),
        "not detected above the current heuristic"
    );
    let missing_description = MainAreaDna::splicing_intron_signal_description_text(&missing);
    assert!(missing_description.contains("not detected above the current heuristic"));
}

#[test]
fn sequencing_trace_selected_base_summary_reports_clip_status() {
    let trace = SequencingTraceRecord {
        schema: "gentle.sequencing_trace_record.v2".to_string(),
        trace_id: "trace_demo".to_string(),
        format: SequencingTraceFormat::AbiAb1,
        source_path: "trace.ab1".to_string(),
        imported_at_unix_ms: 1,
        seq_id: None,
        sample_name: None,
        sample_well: None,
        run_name: None,
        machine_name: None,
        machine_model: None,
        called_bases: "ACGT".to_string(),
        called_base_confidence_values: vec![40, 41, 42, 43],
        peak_locations: vec![10, 20, 30, 40],
        channel_data: vec![SequencingTraceChannelData {
            channel: "A".to_string(),
            trace_set: "DATA9".to_string(),
            points: vec![0; 64],
        }],
        channel_summaries: vec![],
        clip_start_base_index: Some(1),
        clip_end_base_index_exclusive: Some(3),
        comments_text: None,
    };

    let summary =
        MainAreaDna::sequencing_trace_selected_base_summary(&trace, 2).expect("base summary");
    assert!(summary.contains("base 3 / 4"));
    assert!(summary.contains("'G'"));
    assert!(summary.contains("inside clip"));
}

#[test]
fn sequencing_confirmation_selected_read_falls_back_to_first_row() {
    let report = SequencingConfirmationReport {
        reads: vec![
            SequencingConfirmationReadResult {
                evidence_id: "evidence_a".to_string(),
                read_seq_id: "read_a".to_string(),
                ..Default::default()
            },
            SequencingConfirmationReadResult {
                evidence_id: "evidence_b".to_string(),
                read_seq_id: "read_b".to_string(),
                ..Default::default()
            },
        ],
        ..Default::default()
    };

    let selected = MainAreaDna::sequencing_confirmation_selected_read(&report, "missing_evidence")
        .expect("fallback read");
    assert_eq!(selected.evidence_id, "evidence_a");
}

#[test]
fn sequencing_confirmation_selected_evidence_id_prefers_matching_row() {
    let report = SequencingConfirmationReport {
        reads: vec![
            SequencingConfirmationReadResult {
                evidence_id: "evidence_a".to_string(),
                read_seq_id: "read_a".to_string(),
                ..Default::default()
            },
            SequencingConfirmationReadResult {
                evidence_id: "evidence_b".to_string(),
                read_seq_id: "read_b".to_string(),
                ..Default::default()
            },
        ],
        ..Default::default()
    };

    let selected = MainAreaDna::sequencing_confirmation_selected_evidence_id(&report, "EVIDENCE_B")
        .expect("selected evidence id");
    assert_eq!(selected, "evidence_b");
}

#[test]
fn sequencing_confirmation_selected_target_id_prefers_unresolved_rows_when_requested() {
    let report = SequencingConfirmationReport {
        targets: vec![
            SequencingConfirmationTargetResult {
                target_id: "confirmed_target".to_string(),
                label: "Confirmed".to_string(),
                status: SequencingConfirmationStatus::Confirmed,
                start_0based: 0,
                end_0based_exclusive: 10,
                ..Default::default()
            },
            SequencingConfirmationTargetResult {
                target_id: "insufficient_target".to_string(),
                label: "Insufficient".to_string(),
                status: SequencingConfirmationStatus::InsufficientEvidence,
                start_0based: 10,
                end_0based_exclusive: 20,
                ..Default::default()
            },
        ],
        ..Default::default()
    };

    let selected = MainAreaDna::sequencing_confirmation_selected_target_id(&report, "", true)
        .expect("selected target id");
    assert_eq!(selected, "insufficient_target");
}

#[test]
fn sequencing_confirmation_sync_target_selection_updates_linked_rows() {
    let dna = DNAsequence::from_sequence("ACGTACGTACGTACGTACGT").expect("dna");
    let mut area = MainAreaDna::new(dna, Some("expected_seq".to_string()), None);
    let report = SequencingConfirmationReport {
        targets: vec![SequencingConfirmationTargetResult {
            target_id: "target_b".to_string(),
            label: "Target B".to_string(),
            status: SequencingConfirmationStatus::Contradicted,
            start_0based: 4,
            end_0based_exclusive: 12,
            contradicting_read_ids: vec!["evidence_b".to_string()],
            ..Default::default()
        }],
        reads: vec![SequencingConfirmationReadResult {
            evidence_id: "evidence_b".to_string(),
            read_seq_id: "read_b".to_string(),
            trace_id: Some("trace_b".to_string()),
            contradicted_target_ids: vec!["target_b".to_string()],
            covered_target_ids: vec!["target_b".to_string()],
            ..Default::default()
        }],
        variants: vec![SequencingConfirmationVariantRow {
            variant_id: "variant_b".to_string(),
            label: "Variant B".to_string(),
            target_id: Some("target_b".to_string()),
            evidence_id: "evidence_b".to_string(),
            trace_id: Some("trace_b".to_string()),
            start_0based: 6,
            end_0based_exclusive: 7,
            ..Default::default()
        }],
        ..Default::default()
    };

    area.sequencing_confirmation_sync_target_selection(&report, "target_b");

    assert_eq!(
        area.sequencing_confirmation_ui.selected_target_id,
        "target_b"
    );
    assert_eq!(
        area.sequencing_confirmation_ui.selected_evidence_id,
        "evidence_b"
    );
    assert_eq!(
        area.sequencing_confirmation_ui.selected_variant_id,
        "variant_b"
    );
    assert_eq!(area.sequencing_confirmation_ui.selected_trace_id, "trace_b");
}

#[test]
fn sequencing_confirmation_apply_overview_selection_coverage_gap_selects_gap_target() {
    let dna = DNAsequence::from_sequence("ACGTACGTACGTACGTACGT").expect("dna");
    let mut area = MainAreaDna::new(dna, Some("expected_seq".to_string()), None);
    let report = SequencingConfirmationReport {
        targets: vec![
            SequencingConfirmationTargetResult {
                target_id: "confirmed_target".to_string(),
                label: "Confirmed".to_string(),
                status: SequencingConfirmationStatus::Confirmed,
                start_0based: 0,
                end_0based_exclusive: 8,
                ..Default::default()
            },
            SequencingConfirmationTargetResult {
                target_id: "gap_target".to_string(),
                label: "Gap target".to_string(),
                status: SequencingConfirmationStatus::InsufficientEvidence,
                start_0based: 10,
                end_0based_exclusive: 18,
                ..Default::default()
            },
        ],
        ..Default::default()
    };

    let applied = area.sequencing_confirmation_apply_overview_selection(
        &report,
        SequencingConfirmationOverviewSelection::CoverageGap(11, 16),
    );

    assert!(applied);
    assert_eq!(
        area.sequencing_confirmation_ui.selected_target_id,
        "gap_target"
    );
    assert_eq!(
        area.sequencing_confirmation_ui.selected_gap_start_0based,
        Some(11)
    );
    assert_eq!(
        area.sequencing_confirmation_ui
            .selected_gap_end_0based_exclusive,
        Some(16)
    );
}

#[test]
fn sequencing_confirmation_apply_overview_selection_coverage_gap_without_target_keeps_gap_focus() {
    let dna = DNAsequence::from_sequence("ACGTACGTACGTACGTACGT").expect("dna");
    let mut area = MainAreaDna::new(dna, Some("expected_seq".to_string()), None);
    let report = SequencingConfirmationReport::default();

    let applied = area.sequencing_confirmation_apply_overview_selection(
        &report,
        SequencingConfirmationOverviewSelection::CoverageGap(4, 12),
    );

    assert!(applied);
    assert!(
        area.sequencing_confirmation_ui
            .selected_target_id
            .is_empty()
    );
    assert_eq!(
        area.sequencing_confirmation_ui.selected_gap_start_0based,
        Some(4)
    );
    assert_eq!(
        area.sequencing_confirmation_ui
            .selected_gap_end_0based_exclusive,
        Some(12)
    );
}

#[test]
fn sequencing_confirmation_gap_flanking_reads_pick_nearest_edges() {
    let report = SequencingConfirmationReport {
        reads: vec![
            SequencingConfirmationReadResult {
                evidence_id: "left_far".to_string(),
                best_alignment: SequenceAlignmentReport {
                    aligned_target_start_0based: 0,
                    aligned_target_end_0based_exclusive: 8,
                    identity_fraction: 0.90,
                    ..Default::default()
                },
                ..Default::default()
            },
            SequencingConfirmationReadResult {
                evidence_id: "left_near".to_string(),
                best_alignment: SequenceAlignmentReport {
                    aligned_target_start_0based: 3,
                    aligned_target_end_0based_exclusive: 11,
                    identity_fraction: 0.95,
                    ..Default::default()
                },
                ..Default::default()
            },
            SequencingConfirmationReadResult {
                evidence_id: "right_near".to_string(),
                best_alignment: SequenceAlignmentReport {
                    aligned_target_start_0based: 18,
                    aligned_target_end_0based_exclusive: 24,
                    identity_fraction: 0.91,
                    ..Default::default()
                },
                ..Default::default()
            },
            SequencingConfirmationReadResult {
                evidence_id: "right_far".to_string(),
                best_alignment: SequenceAlignmentReport {
                    aligned_target_start_0based: 24,
                    aligned_target_end_0based_exclusive: 30,
                    identity_fraction: 0.99,
                    ..Default::default()
                },
                ..Default::default()
            },
        ],
        ..Default::default()
    };

    let (left, right) = MainAreaDna::sequencing_confirmation_gap_flanking_reads(&report, 12, 18);

    assert_eq!(left.expect("left flank").evidence_id, "left_near");
    assert_eq!(right.expect("right flank").evidence_id, "right_near");
}

#[test]
fn sequencing_confirmation_gap_primer_suggestions_prioritize_covering_spans() {
    let overlay = SequencingPrimerOverlayReport {
        suggestions: vec![
            SequencingPrimerOverlaySuggestion {
                primer_seq_id: "primer_overlap".to_string(),
                orientation: SequencingPrimerOrientation::ForwardRead,
                three_prime_position_0based: 9,
                predicted_read_span_start_0based: 8,
                predicted_read_span_end_0based_exclusive: 20,
                ..Default::default()
            },
            SequencingPrimerOverlaySuggestion {
                primer_seq_id: "primer_miss".to_string(),
                orientation: SequencingPrimerOrientation::ReverseRead,
                three_prime_position_0based: 2,
                predicted_read_span_start_0based: 0,
                predicted_read_span_end_0based_exclusive: 6,
                ..Default::default()
            },
        ],
        ..Default::default()
    };

    let ranked = MainAreaDna::sequencing_confirmation_gap_primer_suggestions(&overlay, 10, 14);

    assert_eq!(ranked[0].primer_seq_id, "primer_overlap");
}

#[test]
fn sequencing_confirmation_unresolved_review_queue_tracks_targets_variants_and_gaps() {
    let report = SequencingConfirmationReport {
        targets: vec![SequencingConfirmationTargetResult {
            target_id: "target_a".to_string(),
            label: "Target A".to_string(),
            status: SequencingConfirmationStatus::Contradicted,
            start_0based: 2,
            end_0based_exclusive: 8,
            ..Default::default()
        }],
        reads: vec![SequencingConfirmationReadResult {
            evidence_id: "evidence_a".to_string(),
            best_alignment: SequenceAlignmentReport {
                aligned_target_start_0based: 0,
                aligned_target_end_0based_exclusive: 6,
                ..Default::default()
            },
            ..Default::default()
        }],
        variants: vec![SequencingConfirmationVariantRow {
            variant_id: "variant_a".to_string(),
            label: "Variant A".to_string(),
            status: SequencingConfirmationStatus::InsufficientEvidence,
            start_0based: 10,
            end_0based_exclusive: 11,
            ..Default::default()
        }],
        ..Default::default()
    };

    let queue = MainAreaDna::sequencing_confirmation_unresolved_review_queue(&report, 20);

    assert_eq!(
        queue,
        vec![
            SequencingConfirmationOverviewSelection::Target("target_a".to_string()),
            SequencingConfirmationOverviewSelection::Variant("variant_a".to_string()),
            SequencingConfirmationOverviewSelection::CoverageGap(6, 20),
        ]
    );
}

#[test]
fn sequencing_confirmation_step_unresolved_focus_advances_from_target_to_variant() {
    let dna = DNAsequence::from_sequence("ACGTACGTACGTACGTACGT").expect("dna");
    let mut area = MainAreaDna::new(dna, Some("expected_seq".to_string()), None);
    let report = SequencingConfirmationReport {
        targets: vec![SequencingConfirmationTargetResult {
            target_id: "target_a".to_string(),
            label: "Target A".to_string(),
            status: SequencingConfirmationStatus::Contradicted,
            start_0based: 2,
            end_0based_exclusive: 8,
            contradicting_read_ids: vec!["evidence_a".to_string()],
            ..Default::default()
        }],
        reads: vec![SequencingConfirmationReadResult {
            evidence_id: "evidence_a".to_string(),
            covered_target_ids: vec!["target_a".to_string()],
            best_alignment: SequenceAlignmentReport {
                aligned_target_start_0based: 0,
                aligned_target_end_0based_exclusive: 6,
                ..Default::default()
            },
            ..Default::default()
        }],
        variants: vec![SequencingConfirmationVariantRow {
            variant_id: "variant_a".to_string(),
            label: "Variant A".to_string(),
            status: SequencingConfirmationStatus::InsufficientEvidence,
            target_id: Some("target_a".to_string()),
            evidence_id: "evidence_a".to_string(),
            start_0based: 10,
            end_0based_exclusive: 11,
            ..Default::default()
        }],
        ..Default::default()
    };

    area.sequencing_confirmation_sync_target_selection(&report, "target_a");

    let moved = area.sequencing_confirmation_step_unresolved_focus(&report, 20, 1);

    assert!(moved);
    assert_eq!(
        area.sequencing_confirmation_ui.selected_review_focus_kind,
        Some(SequencingConfirmationReviewFocusKind::Variant)
    );
    assert_eq!(
        area.sequencing_confirmation_ui.selected_variant_id,
        "variant_a"
    );
}

#[test]
fn extract_rebase_enzyme_names_from_text_handles_ecor_i_spacing() {
    let lookup = MainAreaDna::rebase_name_lookup_by_normalized();
    let names = MainAreaDna::extract_rebase_enzyme_names_from_text(
        "Multiple Cloning Site (MCS); contains BamHI, SmaI and EcoR I",
        &lookup,
    );
    assert!(names.iter().any(|name| name == "BamHI"));
    assert!(names.iter().any(|name| name == "SmaI"));
    assert!(names.iter().any(|name| name == "EcoRI"));
}

#[test]
fn use_mcs_enzymes_prefills_digest_from_mcs_feature() {
    let mut dna = restriction_ready_dna("TTTGAATTCTTTGGATCCTTT");
    dna.features_mut().push(Feature {
        kind: "misc_feature".into(),
        location: Location::simple_range(1, 5),
        qualifiers: vec![
            (
                "label".into(),
                Some("Multiple Cloning Site (MCS)".to_string()),
            ),
            (
                "note".into(),
                Some("contains BamHI, SmaI and EcoR I".to_string()),
            ),
        ],
    });
    let mut area = MainAreaDna::new(dna, Some("s".to_string()), None);
    area.set_digest_enzymes_from_mcs();
    let selected = MainAreaDna::parse_ids(&area.digest_enzymes_text);
    assert!(selected.iter().any(|name| name == "BamHI"));
    assert!(selected.iter().any(|name| name == "EcoRI"));
    assert!(!selected.iter().any(|name| name == "SmaI"));
    assert!(area.op_status.contains("MCS-linked enzyme"));
}

#[test]
fn single_cutter_cds_filter_excludes_non_cds_single_cutters() {
    let mut dna = restriction_ready_dna(&format!(
        "{}{}{}{}",
        "GAATTC",
        "A".repeat(16),
        "GGATCC",
        "A".repeat(20)
    ));
    dna.features_mut().push(Feature {
        kind: "CDS".into(),
        location: Location::simple_range(20, 30),
        qualifiers: vec![("label".into(), Some("coding".to_string()))],
    });
    let mut area = MainAreaDna::new(dna, Some("s".to_string()), None);
    area.set_digest_enzymes_single_cutters(false);
    let all_single = MainAreaDna::parse_ids(&area.digest_enzymes_text);
    assert!(all_single.iter().any(|name| name == "EcoRI"));
    assert!(all_single.iter().any(|name| name == "BamHI"));

    area.set_digest_enzymes_single_cutters(true);
    let cds_single = MainAreaDna::parse_ids(&area.digest_enzymes_text);
    assert!(cds_single.iter().any(|name| name == "BamHI"));
    assert!(!cds_single.iter().any(|name| name == "EcoRI"));
}

#[test]
fn restriction_layer_counts_respect_display_mode_and_total_sites() {
    let mut dna = DNAsequence::from_sequence("GAATTCAAAAGAATTCAAAAGGATCCAAA").expect("seq");
    *dna.restriction_enzymes_mut() = active_restriction_enzymes()
        .into_iter()
        .filter(|enzyme| matches!(enzyme.name.as_str(), "EcoRI" | "BamHI"))
        .collect();
    dna.set_max_restriction_enzyme_sites(None);
    dna.update_computed_features();
    let area = MainAreaDna::new(dna, Some("s".to_string()), None);
    {
        let mut display = area.dna_display.write().expect("display");
        display.set_restriction_enzyme_display_mode(RestrictionEnzymeDisplayMode::PreferredOnly);
        display.set_preferred_restriction_enzymes(vec!["EcoRI".to_string()]);
    }
    let preferred_counts = area.compute_layer_visibility_counts();
    assert_eq!(preferred_counts.restriction_site_total_count, 3);
    assert_eq!(preferred_counts.restriction_site_preferred_count, 2);
    assert_eq!(preferred_counts.restriction_site_unique_count, 1);
    assert_eq!(preferred_counts.restriction_site_count, 2);

    area.dna_display
        .write()
        .expect("display")
        .set_restriction_enzyme_display_mode(RestrictionEnzymeDisplayMode::PreferredAndUnique);
    let preferred_and_unique_counts = area.compute_layer_visibility_counts();
    assert_eq!(preferred_and_unique_counts.restriction_site_count, 3);

    area.dna_display
        .write()
        .expect("display")
        .set_restriction_enzyme_display_mode(RestrictionEnzymeDisplayMode::AllInView);
    let all_counts = area.compute_layer_visibility_counts();
    assert_eq!(all_counts.restriction_site_count, 3);
}

#[test]
fn collect_tfbs_scan_settings_merges_and_sorts_per_tf_thresholds() {
    let dna = DNAsequence::from_sequence("TTTTATAAAGGGTATAAATTT").expect("sequence");
    let mut area = MainAreaDna::new(dna, Some("seq1".to_string()), None);
    area.tfbs_motifs = "TP73,SP1".to_string();
    area.tfbs_min_llr_bits = "1.5".to_string();
    area.tfbs_min_llr_quantile = "0.8".to_string();
    area.tfbs_per_tf_min_llr_bits = "TP73=0.5,SP1=1.25".to_string();
    area.tfbs_per_tf_min_llr_quantile = "TP73=0.95".to_string();

    let (motifs, min_bits, min_quantile, overrides) = area
        .collect_tfbs_scan_settings()
        .expect("tfbs scan settings");

    assert_eq!(motifs, vec!["TP73".to_string(), "SP1".to_string()]);
    assert_eq!(min_bits, Some(1.5));
    assert_eq!(min_quantile, Some(0.8));
    assert_eq!(overrides.len(), 2);
    assert_eq!(overrides[0].tf, "SP1");
    assert_eq!(overrides[0].min_llr_bits, Some(1.25));
    assert_eq!(overrides[0].min_llr_quantile, None);
    assert_eq!(overrides[1].tf, "TP73");
    assert_eq!(overrides[1].min_llr_bits, Some(0.5));
    assert_eq!(overrides[1].min_llr_quantile, Some(0.95));
}

#[test]
fn collect_tfbs_scan_settings_requires_motif_selection() {
    let dna = DNAsequence::from_sequence("TTTTATAAAGGGTATAAATTT").expect("sequence");
    let area = MainAreaDna::new(dna, Some("seq1".to_string()), None);

    let err = area
        .collect_tfbs_scan_settings()
        .expect_err("missing motif selection should fail");

    assert!(err.contains("Provide at least one TF motif"));
}

#[test]
fn current_engine_ops_state_records_tfbs_score_track_controls() {
    let dna = DNAsequence::from_sequence("TTTTATAAAGGGTATAAATTT").expect("sequence");
    let mut area = MainAreaDna::new(dna, Some("seq1".to_string()), None);
    area.tfbs_score_track_value_kind = TfbsScoreTrackValueKind::TrueLogOddsBackgroundTailLog10;
    area.tfbs_score_track_clip_negative = false;

    let state = area.current_engine_ops_state();

    assert_eq!(
        state.tfbs_score_track_value_kind,
        TfbsScoreTrackValueKind::TrueLogOddsBackgroundTailLog10
    );
    assert!(!state.tfbs_score_track_clip_negative);
}

#[test]
fn current_engine_ops_state_records_tfbs_track_similarity_controls() {
    let dna = DNAsequence::from_sequence("TTTTATAAAGGGTATAAATTT").expect("sequence");
    let mut area = MainAreaDna::new(dna, Some("seq1".to_string()), None);
    area.tfbs_track_similarity_anchor_motif = "TP53".to_string();
    area.tfbs_track_similarity_all_candidates = false;
    area.tfbs_track_similarity_candidate_motifs = "TP63,TP73".to_string();
    area.tfbs_track_similarity_species_filters = "human,mouse".to_string();
    area.tfbs_track_similarity_include_remote_metadata = true;
    area.tfbs_track_similarity_limit = "15".to_string();
    area.tfbs_track_similarity_ranking_metric = TfbsTrackSimilarityRankingMetric::SmoothedPearson;

    let state = area.current_engine_ops_state();

    assert_eq!(state.tfbs_track_similarity_anchor_motif, "TP53");
    assert!(!state.tfbs_track_similarity_all_candidates);
    assert_eq!(state.tfbs_track_similarity_candidate_motifs, "TP63,TP73");
    assert_eq!(state.tfbs_track_similarity_species_filters, "human,mouse");
    assert!(state.tfbs_track_similarity_include_remote_metadata);
    assert_eq!(state.tfbs_track_similarity_limit, "15");
    assert_eq!(
        state.tfbs_track_similarity_ranking_metric,
        TfbsTrackSimilarityRankingMetric::SmoothedPearson
    );
}

#[test]
fn tfbs_score_track_controls_default_when_missing_in_serialized_engine_ops_state() {
    let dna = DNAsequence::from_sequence("ACGT").unwrap();
    let area = MainAreaDna::new(dna, None, None);
    let mut value = serde_json::to_value(area.current_engine_ops_state()).unwrap();
    value
        .as_object_mut()
        .unwrap()
        .remove("tfbs_score_track_value_kind");
    value
        .as_object_mut()
        .unwrap()
        .remove("tfbs_score_track_clip_negative");

    let decoded: super::EngineOpsUiState = serde_json::from_value(value).unwrap();

    assert_eq!(
        decoded.tfbs_score_track_value_kind,
        TfbsScoreTrackValueKind::LlrBits
    );
    assert!(decoded.tfbs_score_track_clip_negative);
}

#[test]
fn tfbs_track_similarity_controls_default_when_missing_in_serialized_engine_ops_state() {
    let dna = DNAsequence::from_sequence("ACGT").unwrap();
    let area = MainAreaDna::new(dna, None, None);
    let mut value = serde_json::to_value(area.current_engine_ops_state()).unwrap();
    let object = value.as_object_mut().unwrap();
    object.remove("tfbs_track_similarity_anchor_motif");
    object.remove("tfbs_track_similarity_all_candidates");
    object.remove("tfbs_track_similarity_candidate_motifs");
    object.remove("tfbs_track_similarity_species_filters");
    object.remove("tfbs_track_similarity_include_remote_metadata");
    object.remove("tfbs_track_similarity_limit");
    object.remove("tfbs_track_similarity_ranking_metric");

    let decoded: super::EngineOpsUiState = serde_json::from_value(value).unwrap();

    assert!(decoded.tfbs_track_similarity_anchor_motif.is_empty());
    assert!(decoded.tfbs_track_similarity_all_candidates);
    assert!(decoded.tfbs_track_similarity_candidate_motifs.is_empty());
    assert!(decoded.tfbs_track_similarity_species_filters.is_empty());
    assert!(!decoded.tfbs_track_similarity_include_remote_metadata);
    assert!(decoded.tfbs_track_similarity_limit.is_empty());
    assert_eq!(
        decoded.tfbs_track_similarity_ranking_metric,
        TfbsTrackSimilarityRankingMetric::SmoothedSpearman
    );
}

#[test]
fn whole_sequence_tfbs_scan_uses_shared_engine_report_path() {
    let dna = DNAsequence::from_sequence("TTTTATAAAGGGTATAAATTT").expect("sequence");
    let mut state = ProjectState::default();
    state.sequences.insert("seq1".to_string(), dna.clone());
    let engine = Arc::new(RwLock::new(GentleEngine::from_state(state)));
    let mut area = MainAreaDna::new(dna, Some("seq1".to_string()), Some(engine));
    area.tfbs_motifs = "TATAAA".to_string();

    area.scan_whole_sequence_for_tfbs_hits();

    assert!(area.op_status.contains("TFBS scan for whole sequence"));
    assert!(area.op_status.contains("matched"));
    assert!(area.op_status.contains("TATAAA@3+"));
    let report = area
        .cached_tfbs_hit_scan
        .as_ref()
        .expect("cached TFBS hit scan");
    assert_eq!(report.scan_start_0based, 0);
    assert_eq!(report.scan_end_0based_exclusive, 21);
    assert_eq!(report.motifs_requested, vec!["TATAAA".to_string()]);
    assert!(!report.rows.is_empty());
    assert!(
        report
            .rows
            .iter()
            .any(|row| row.source_match_start_0based == 3 && row.forward_strand)
    );
}

#[test]
fn whole_sequence_restriction_site_scan_caches_shared_engine_report() {
    let dna = DNAsequence::from_sequence("GAATTCCCGGGATCC").expect("sequence");
    let mut state = ProjectState::default();
    state.sequences.insert("seq1".to_string(), dna.clone());
    let engine = Arc::new(RwLock::new(GentleEngine::from_state(state)));
    let mut area = MainAreaDna::new(dna, Some("seq1".to_string()), Some(engine));

    area.scan_whole_sequence_for_restriction_sites();

    assert!(
        area.op_status
            .contains("Restriction-site scan for whole sequence")
    );
    let report = area
        .cached_restriction_site_scan
        .as_ref()
        .expect("cached restriction-site scan");
    assert_eq!(report.scan_start_0based, 0);
    assert_eq!(report.scan_end_0based_exclusive, 15);
    assert!(
        report.rows.iter().any(|row| row.enzyme_name == "EcoRI"),
        "EcoRI should be present in the cached scan report"
    );
    assert!(
        report.rows.iter().any(|row| row.enzyme_name == "SmaI"),
        "SmaI should be present in the cached scan report"
    );
}

#[test]
fn whole_sequence_tfbs_score_tracks_use_shared_engine_report_path() {
    let dna = DNAsequence::from_sequence("TTTTATAAAGGGTATAAATTT").expect("sequence");
    let mut state = ProjectState::default();
    state.sequences.insert("seq1".to_string(), dna.clone());
    let engine = Arc::new(RwLock::new(GentleEngine::from_state(state)));
    let mut area = MainAreaDna::new(dna, Some("seq1".to_string()), Some(engine));
    area.tfbs_motifs = "TATAAA".to_string();

    area.show_whole_sequence_tfbs_score_tracks();
    assert!(area.tfbs_task.is_some());

    let ctx = eframe::egui::Context::default();
    let wait_started = Instant::now();
    while area.tfbs_task.is_some() && wait_started.elapsed() < Duration::from_secs(15) {
        area.poll_tfbs_task(&ctx);
        std::thread::sleep(Duration::from_millis(2));
    }
    if area.tfbs_task.is_some() {
        area.poll_tfbs_task(&ctx);
    }
    assert!(
        area.tfbs_task.is_none(),
        "tfbs score-track task should finish"
    );
    assert!(
        area.op_status
            .contains("TFBS score tracks for whole sequence")
    );
    assert!(area.op_status.contains("0..21"));
    let report = area
        .cached_tfbs_score_tracks
        .as_ref()
        .expect("cached TFBS score tracks");
    assert_eq!(report.view_start_0based, 0);
    assert_eq!(report.view_end_0based_exclusive, 21);
    assert_eq!(report.motifs_requested, vec!["TATAAA".to_string()]);
}

#[test]
fn whole_sequence_tfbs_track_similarity_uses_shared_engine_report_path() {
    let dna = DNAsequence::from_sequence(&"ACGT".repeat(16)).expect("sequence");
    let mut state = ProjectState::default();
    state.sequences.insert("seq1".to_string(), dna.clone());
    let engine = Arc::new(RwLock::new(GentleEngine::from_state(state)));
    let mut area = MainAreaDna::new(dna, Some("seq1".to_string()), Some(engine));
    area.tfbs_track_similarity_anchor_motif = "ACGT".to_string();
    area.tfbs_track_similarity_all_candidates = false;
    area.tfbs_track_similarity_candidate_motifs = "CGTA,TACG".to_string();

    area.show_whole_sequence_tfbs_track_similarity();
    assert!(area.tfbs_task.is_some());

    let ctx = eframe::egui::Context::default();
    let wait_started = Instant::now();
    while area.tfbs_task.is_some() && wait_started.elapsed() < Duration::from_secs(15) {
        area.poll_tfbs_task(&ctx);
        std::thread::sleep(Duration::from_millis(2));
    }
    if area.tfbs_task.is_some() {
        area.poll_tfbs_task(&ctx);
    }
    assert!(
        area.tfbs_task.is_none(),
        "tfbs similarity task should finish"
    );
    assert!(
        area.op_status
            .contains("TFBS similarity for whole sequence")
    );
    let report = area
        .cached_tfbs_track_similarity
        .as_ref()
        .expect("cached TFBS similarity report");
    assert_eq!(report.view_start_0based, 0);
    assert_eq!(report.view_end_0based_exclusive, 64);
    assert_eq!(report.anchor_requested, "ACGT");
    assert_eq!(
        report.candidates_requested,
        vec!["CGTA".to_string(), "TACG".to_string()]
    );
    assert_eq!(
        report.ranking_metric,
        TfbsTrackSimilarityRankingMetric::SmoothedSpearman
    );
    assert!(!report.rows.is_empty());
}

#[test]
fn cutrun_regulatory_support_gui_uses_shared_engine_report_path() {
    let td = tempdir().expect("tempdir");
    let root = td.path();
    let sp1_consensus = "TATAAA";
    let supported_prefix = "TTGACCAA";
    let sequence = format!("{supported_prefix}{sp1_consensus}AACCGGTT");
    let fasta = root.join("toy.fa");
    let gtf = root.join("toy.gtf");
    fs::write(&fasta, format!(">chr1\n{sequence}\n")).expect("write FASTA");
    fs::write(
        &gtf,
        format!(
            "chr1\tsrc\tgene\t1\t{}\t.\t+\t.\tgene_id \"GENE1\"; gene_name \"GENE1\";\n",
            sequence.len()
        ),
    )
    .expect("write GTF");
    let catalog_path = root.join("catalog.json");
    let cache_dir = root.join("cache");
    fs::write(
        &catalog_path,
        serde_json::to_string_pretty(&json!({
            "ToyGenome": {
                "sequence_local": fasta.to_string_lossy(),
                "annotations_local": gtf.to_string_lossy(),
                "cache_dir": cache_dir.to_string_lossy()
            }
        }))
        .expect("serialize catalog"),
    )
    .expect("write catalog");
    let mut engine_value = GentleEngine::new();
    engine_value
        .apply(Operation::PrepareGenome {
            genome_id: "ToyGenome".to_string(),
            catalog_path: Some(catalog_path.to_string_lossy().to_string()),
            cache_dir: None,
            timeout_seconds: None,
        })
        .expect("prepare toy genome");
    engine_value
        .apply(Operation::ExtractGenomeRegion {
            genome_id: "ToyGenome".to_string(),
            chromosome: "chr1".to_string(),
            start_1based: 1,
            end_1based: sequence.len(),
            output_id: Some("seq1".to_string()),
            annotation_scope: None,
            max_annotation_features: None,
            include_genomic_annotation: None,
            catalog_path: Some(catalog_path.to_string_lossy().to_string()),
            cache_dir: None,
        })
        .expect("extract anchored toy ROI");
    let dna = {
        let dna = engine_value
            .state_mut()
            .sequences
            .get_mut("seq1")
            .expect("extracted sequence");
        dna.features_mut().push(Feature {
            kind: "TFBS".into(),
            location: Location::simple_range(
                supported_prefix.len() as i64,
                (supported_prefix.len() + sp1_consensus.len()) as i64,
            ),
            qualifiers: vec![
                ("label".into(), Some("TATAAA".to_string())),
                ("bound_moiety".into(), Some("TATAAA".to_string())),
            ],
        });
        dna.clone()
    };
    let engine = Arc::new(RwLock::new(engine_value));
    let reads_path = root.join("cutrun_reads.fa");
    fs::write(
        &reads_path,
        format!(
            ">support_a\n{supported_prefix}{sp1_consensus}\n>support_b\n{supported_prefix}{sp1_consensus}\n"
        ),
    )
    .expect("write CUT&RUN reads");
    engine
        .write()
        .expect("engine lock")
        .apply(Operation::InterpretCutRunReads {
            seq_id: "seq1".to_string(),
            input_r1_path: Some(reads_path.to_string_lossy().to_string()),
            input_r2_path: None,
            dataset_id: None,
            catalog_path: None,
            cache_dir: None,
            input_format: CutRunInputFormat::Fasta,
            read_layout: CutRunReadLayout::SingleEnd,
            roi_flank_bp: 0,
            seed_filter: CutRunSeedFilterConfig {
                kmer_len: 4,
                min_seed_matches: 1,
            },
            align_config: CutRunAlignConfig {
                max_mismatches: 0,
                min_identity_fraction: 1.0,
                max_fragment_span_bp: 64,
            },
            deduplicate_fragments: false,
            report_id: Some("gui_cutrun_reads".to_string()),
            checkpoint_path: None,
            checkpoint_every_reads: 10,
        })
        .expect("interpret CUT&RUN reads");

    let mut area = MainAreaDna::new(dna, Some("seq1".to_string()), Some(engine));
    area.cutrun_regulatory_read_report_ids = "gui_cutrun_reads".to_string();

    area.inspect_cutrun_regulatory_support_for_active_sequence();

    assert!(
        area.op_status
            .contains("CUT&RUN regulatory support for 'seq1'"),
        "expected CUT&RUN summary status, got: {}",
        area.op_status
    );
    assert!(
        area.op_status.contains("merged 1 source(s)"),
        "expected merged-source count in status, got: {}",
        area.op_status
    );
    let report = area
        .cached_cutrun_regulatory_support
        .as_ref()
        .expect("cached CUT&RUN regulatory-support report");
    assert_eq!(report.seq_id, "seq1");
    assert_eq!(report.evidence_sources.len(), 1);
    assert!(!report.support_windows.is_empty());
    assert_eq!(report.confirmed_tfbs_rows.len(), 1);
    assert_eq!(
        report.confirmed_tfbs_rows[0].motif_label.as_deref(),
        Some("TATAAA")
    );
}

#[test]
fn whole_sequence_tfbs_score_tracks_can_be_cancelled() {
    let dna = DNAsequence::from_sequence("TTTTATAAAGGGTATAAATTT").expect("sequence");
    let mut state = ProjectState::default();
    state.sequences.insert("seq1".to_string(), dna.clone());
    let engine = Arc::new(RwLock::new(GentleEngine::from_state(state)));
    let mut area = MainAreaDna::new(dna, Some("seq1".to_string()), Some(engine));
    area.tfbs_motifs = "TATAAA".to_string();

    area.show_whole_sequence_tfbs_score_tracks();
    let task = area.tfbs_task.as_ref().expect("tfbs task started");
    task.cancel_requested
        .store(true, std::sync::atomic::Ordering::Relaxed);

    let ctx = eframe::egui::Context::default();
    for _ in 0..300 {
        area.poll_tfbs_task(&ctx);
        if area.tfbs_task.is_none() {
            break;
        }
        std::thread::sleep(Duration::from_millis(2));
    }
    assert!(
        area.tfbs_task.is_none(),
        "cancelled tfbs score-track task should stop"
    );
    assert!(area.cached_tfbs_score_tracks.is_none());
    assert!(area.op_status.to_ascii_lowercase().contains("cancelled"));
}

#[test]
fn parse_positive_usize_text_rejects_zero_and_non_integer() {
    assert!(
        MainAreaDna::parse_positive_usize_text("0", "extension length")
            .unwrap_err()
            .contains(">= 1")
    );
    assert!(
        MainAreaDna::parse_positive_usize_text("abc", "extension length")
            .unwrap_err()
            .contains("expected an integer")
    );
}

#[test]
fn dotplot_window_count_respects_word_and_step() {
    assert_eq!(MainAreaDna::dotplot_window_count(1001, 11, 11), 91);
    assert_eq!(MainAreaDna::dotplot_window_count(1001, 9, 4), 249);
    assert_eq!(MainAreaDna::dotplot_window_count(8, 9, 1), 0);
    assert_eq!(MainAreaDna::dotplot_window_count(100, 0, 1), 0);
}

#[test]
fn dotplot_svg_default_filename_mentions_parameters() {
    let dna = DNAsequence::from_sequence("ACGTACGTACGTACGT").expect("sequence");
    let area = MainAreaDna::new(dna, Some("seq1".to_string()), None);
    let mut view = DotplotView::default();
    view.dotplot_id = "tp73.dotplot.primary".to_string();
    view.seq_id = "tp73.ncbi".to_string();
    view.reference_seq_id = Some("NC_000001".to_string());
    view.span_start_0based = 3_652_306;
    view.span_end_0based = 3_736_201;
    view.reference_span_start_0based = 3_652_306;
    view.reference_span_end_0based = 3_736_201;
    view.mode = DotplotMode::PairForward;
    view.word_size = 9;
    view.step_bp = 4;
    view.max_mismatches = 1;
    view.tile_bp = Some(250);
    let mut track = FlexibilityTrack::default();
    track.track_id = "tp73.flex.track".to_string();
    track.model = FlexibilityModel::AtRichness;
    track.bin_bp = 50;
    track.smoothing_bp = Some(150);
    let file_name = area.default_dotplot_svg_file_name(
        &view,
        Some(&track),
        0.35,
        1.75,
        DotplotOverlayXAxisMode::PercentLength,
        None,
    );
    assert!(file_name.ends_with(".svg"));
    assert!(file_name.contains("pair_forward"));
    assert!(file_name.contains("_w9_s4_mm1_tile250_"));
    assert!(file_name.contains("_th35_gain175"));
    assert!(file_name.contains("_fx-tp73_flex_track_at_richness_b50_sm150"));
}

#[test]
fn dotplot_svg_export_document_contains_summary_and_cells() {
    let mut view = DotplotView::default();
    view.dotplot_id = "toy_plot".to_string();
    view.seq_id = "q1".to_string();
    view.reference_seq_id = Some("r1".to_string());
    view.mode = DotplotMode::SelfForward;
    view.span_start_0based = 0;
    view.span_end_0based = 200;
    view.reference_span_start_0based = 0;
    view.reference_span_end_0based = 200;
    view.word_size = 7;
    view.step_bp = 1;
    view.max_mismatches = 1;
    view.point_count = 3;
    view.points = vec![
        crate::engine::DotplotMatchPoint {
            x_0based: 10,
            y_0based: 10,
            mismatches: 0,
        },
        crate::engine::DotplotMatchPoint {
            x_0based: 30,
            y_0based: 35,
            mismatches: 1,
        },
        crate::engine::DotplotMatchPoint {
            x_0based: 150,
            y_0based: 151,
            mismatches: 0,
        },
    ];
    let svg = MainAreaDna::build_dotplot_svg_document(&view, None, 0.10, 1.50, Some((25, 26)));
    assert!(svg.starts_with("<svg "));
    assert!(svg.contains("Dotplot workspace export: toy_plot"));
    assert!(svg.contains("mode=self_forward"));
    assert!(svg.contains("threshold=0.10 gain=1.50"));
    assert!(svg.contains("overlap by 6 bp"));
    assert!(svg.contains("7 consecutive ordered windows"));
    assert!(svg.contains("GENtle dotplot SVG export"));
    assert!(svg.contains("<rect "));
}

#[test]
fn dotplot_svg_default_filename_mentions_overlay_owner_and_series_count() {
    let dna = DNAsequence::from_sequence("ACGTACGTACGTACGT").expect("sequence");
    let area = MainAreaDna::new(dna, Some("ref".to_string()), None);
    let mut view = DotplotView::default();
    view.dotplot_id = "tp53.overlay".to_string();
    view.owner_seq_id = "tp53_genomic".to_string();
    view.seq_id = "tp53_201".to_string();
    view.reference_seq_id = Some("tp53_genomic".to_string());
    view.reference_span_start_0based = 7_661_778;
    view.reference_span_end_0based = 7_687_546;
    view.mode = DotplotMode::PairForward;
    view.word_size = 7;
    view.step_bp = 1;
    view.max_mismatches = 0;
    view.series_count = 2;
    view.query_series = vec![
        crate::engine::DotplotQuerySeries {
            series_id: "tp53.overlay.series1".to_string(),
            seq_id: "tp53_201".to_string(),
            label: "TP53-201".to_string(),
            transcript_feature_id: None,
            query_anchor_0based: None,
            query_anchor_label: None,
            color_rgb: [29, 78, 216],
            mode: DotplotMode::PairForward,
            span_start_0based: 0,
            span_end_0based: 1_234,
            point_count: 0,
            points: vec![],
            boxplot_bin_count: 0,
            boxplot_bins: vec![],
        },
        crate::engine::DotplotQuerySeries {
            series_id: "tp53.overlay.series2".to_string(),
            seq_id: "tp53_202".to_string(),
            label: "TP53-202".to_string(),
            transcript_feature_id: None,
            query_anchor_0based: None,
            query_anchor_label: None,
            color_rgb: [220, 38, 38],
            mode: DotplotMode::PairForward,
            span_start_0based: 0,
            span_end_0based: 1_118,
            point_count: 0,
            points: vec![],
            boxplot_bin_count: 0,
            boxplot_bins: vec![],
        },
    ];
    let file_name = area.default_dotplot_svg_file_name(
        &view,
        None,
        0.20,
        1.25,
        DotplotOverlayXAxisMode::RightAlignedBp,
        None,
    );
    assert!(file_name.ends_with(".svg"));
    assert!(file_name.contains("owner-tp53_genomic"));
    assert!(file_name.contains("overlay2_series"));
    assert!(file_name.contains("pair_forward"));
    assert!(file_name.contains("x-right_aligned_bp"));
}

#[test]
fn dotplot_svg_default_filename_mentions_shared_exon_anchor() {
    let dna = DNAsequence::from_sequence("ACGTACGTACGTACGT").expect("sequence");
    let area = MainAreaDna::new(dna, Some("ref".to_string()), None);
    let mut view = DotplotView::default();
    view.dotplot_id = "toy.anchor.overlay".to_string();
    view.owner_seq_id = "toy_anchor_locus".to_string();
    view.seq_id = "toy_anchor_tx1".to_string();
    view.reference_seq_id = Some("toy_anchor_locus".to_string());
    view.reference_span_start_0based = 0;
    view.reference_span_end_0based = 240;
    view.mode = DotplotMode::PairForward;
    view.word_size = 7;
    view.step_bp = 1;
    view.max_mismatches = 0;
    view.series_count = 3;
    view.query_series = vec![
        crate::engine::DotplotQuerySeries {
            series_id: "s1".to_string(),
            seq_id: "toy_anchor_tx1".to_string(),
            label: "toy-1".to_string(),
            color_rgb: [29, 78, 216],
            transcript_feature_id: Some(0),
            query_anchor_0based: None,
            query_anchor_label: None,
            mode: DotplotMode::PairForward,
            span_start_0based: 0,
            span_end_0based: 120,
            point_count: 0,
            points: vec![],
            boxplot_bin_count: 0,
            boxplot_bins: vec![],
        },
        crate::engine::DotplotQuerySeries {
            series_id: "s2".to_string(),
            seq_id: "toy_anchor_tx2".to_string(),
            label: "toy-2".to_string(),
            color_rgb: [220, 38, 38],
            transcript_feature_id: Some(1),
            query_anchor_0based: None,
            query_anchor_label: None,
            mode: DotplotMode::PairForward,
            span_start_0based: 0,
            span_end_0based: 132,
            point_count: 0,
            points: vec![],
            boxplot_bin_count: 0,
            boxplot_bins: vec![],
        },
        crate::engine::DotplotQuerySeries {
            series_id: "s3".to_string(),
            seq_id: "toy_anchor_tx3".to_string(),
            label: "toy-3".to_string(),
            color_rgb: [5, 150, 105],
            transcript_feature_id: Some(2),
            query_anchor_0based: None,
            query_anchor_label: None,
            mode: DotplotMode::PairForward,
            span_start_0based: 0,
            span_end_0based: 140,
            point_count: 0,
            points: vec![],
            boxplot_bin_count: 0,
            boxplot_bins: vec![],
        },
    ];
    let file_name = area.default_dotplot_svg_file_name(
        &view,
        None,
        0.08,
        2.80,
        DotplotOverlayXAxisMode::SharedExonAnchor,
        Some(&DotplotOverlayAnchorExonRef {
            start_1based: 181,
            end_1based: 215,
        }),
    );
    assert!(file_name.contains("x-shared_exon_anchor"));
    assert!(file_name.contains("anchor-181_215"));
}

#[test]
fn dotplot_svg_default_filename_mentions_query_anchor_label() {
    let dna = DNAsequence::from_sequence("ACGTACGTACGTACGT").expect("sequence");
    let area = MainAreaDna::new(dna, Some("ref".to_string()), None);
    let mut view = DotplotView::default();
    view.dotplot_id = "p53.family.overlay".to_string();
    view.owner_seq_id = "tp73_family_ref".to_string();
    view.seq_id = "tp73_family_ref".to_string();
    view.reference_seq_id = Some("tp73_family_ref".to_string());
    view.reference_span_start_0based = 0;
    view.reference_span_end_0based = 5192;
    view.mode = DotplotMode::PairForward;
    view.word_size = 9;
    view.step_bp = 2;
    view.max_mismatches = 0;
    view.series_count = 3;
    view.query_series = vec![
        crate::engine::DotplotQuerySeries {
            series_id: "tp73".to_string(),
            seq_id: "tp73_family_ref".to_string(),
            label: "TP73".to_string(),
            color_rgb: [29, 78, 216],
            transcript_feature_id: None,
            query_anchor_0based: Some(926),
            query_anchor_label: Some("shared core motif".to_string()),
            mode: DotplotMode::PairForward,
            span_start_0based: 0,
            span_end_0based: 5192,
            point_count: 0,
            points: vec![],
            boxplot_bin_count: 0,
            boxplot_bins: vec![],
        },
        crate::engine::DotplotQuerySeries {
            series_id: "tp63".to_string(),
            seq_id: "tp63_family_query".to_string(),
            label: "TP63".to_string(),
            color_rgb: [220, 38, 38],
            transcript_feature_id: None,
            query_anchor_0based: Some(1044),
            query_anchor_label: Some("shared core motif".to_string()),
            mode: DotplotMode::PairForward,
            span_start_0based: 0,
            span_end_0based: 4944,
            point_count: 0,
            points: vec![],
            boxplot_bin_count: 0,
            boxplot_bins: vec![],
        },
        crate::engine::DotplotQuerySeries {
            series_id: "tp53".to_string(),
            seq_id: "tp53_family_query".to_string(),
            label: "TP53".to_string(),
            color_rgb: [5, 150, 105],
            transcript_feature_id: None,
            query_anchor_0based: Some(707),
            query_anchor_label: Some("shared core motif".to_string()),
            mode: DotplotMode::PairForward,
            span_start_0based: 0,
            span_end_0based: 1018,
            point_count: 0,
            points: vec![],
            boxplot_bin_count: 0,
            boxplot_bins: vec![],
        },
    ];
    let file_name = area.default_dotplot_svg_file_name(
        &view,
        None,
        0.08,
        2.20,
        DotplotOverlayXAxisMode::QueryAnchorBp,
        None,
    );
    assert!(file_name.contains("x-query_anchor_bp"));
    assert!(file_name.contains("anchor-shared_core_motif"));
}

#[test]
fn dotplot_query_context_prefers_override_sequence() {
    let dna = DNAsequence::from_sequence("ACGTACGT").expect("sequence");
    let mut area = MainAreaDna::new(dna, Some("seq1".to_string()), None);
    assert_eq!(area.current_dotplot_query_seq_id().as_deref(), Some("seq1"));
    assert!(!area.using_dotplot_query_override());
    area.dotplot_query_override_seq_id = "read_query".to_string();
    area.dotplot_query_override_source_label = "RNA-read #7".to_string();
    assert_eq!(
        area.current_dotplot_query_seq_id().as_deref(),
        Some("read_query")
    );
    assert!(area.using_dotplot_query_override());
    assert!(area.current_dotplot_query_label().contains("RNA-read #7"));
}

#[test]
fn dotplot_selection_sync_is_disabled_for_query_override_views() {
    let dna = DNAsequence::from_sequence("ACGTACGT").expect("sequence");
    let mut area = MainAreaDna::new(dna, Some("seq1".to_string()), None);
    let mut view = DotplotView::default();
    view.seq_id = "seq1".to_string();
    assert!(area.dotplot_selection_sync_enabled_for_view(&view));
    area.dotplot_query_override_seq_id = "read_query".to_string();
    view.seq_id = "read_query".to_string();
    assert!(!area.dotplot_selection_sync_enabled_for_view(&view));
}

#[test]
fn dotplot_selection_sync_is_disabled_for_overlay_views() {
    let dna = DNAsequence::from_sequence("ACGTACGT").expect("sequence");
    let area = MainAreaDna::new(dna, Some("ref".to_string()), None);
    let mut view = DotplotView::default();
    view.seq_id = "iso_a".to_string();
    view.series_count = 2;
    view.query_series = vec![
        crate::engine::DotplotQuerySeries {
            series_id: "overlay.series1".to_string(),
            seq_id: "iso_a".to_string(),
            label: "Isoform A".to_string(),
            transcript_feature_id: None,
            query_anchor_0based: None,
            query_anchor_label: None,
            color_rgb: [29, 78, 216],
            mode: DotplotMode::PairForward,
            span_start_0based: 0,
            span_end_0based: 10,
            point_count: 0,
            points: vec![],
            boxplot_bin_count: 0,
            boxplot_bins: vec![],
        },
        crate::engine::DotplotQuerySeries {
            series_id: "overlay.series2".to_string(),
            seq_id: "iso_b".to_string(),
            label: "Isoform B".to_string(),
            transcript_feature_id: None,
            query_anchor_0based: None,
            query_anchor_label: None,
            color_rgb: [220, 38, 38],
            mode: DotplotMode::PairForward,
            span_start_0based: 0,
            span_end_0based: 10,
            point_count: 0,
            points: vec![],
            boxplot_bin_count: 0,
            boxplot_bins: vec![],
        },
    ];
    assert!(!area.dotplot_selection_sync_enabled_for_view(&view));
}

#[test]
fn dotplot_reference_normalization_replaces_stale_genomic_span_for_short_query_reference() {
    let dna = DNAsequence::from_sequence("ACGTACGTACGT").expect("sequence");
    let mut state = ProjectState::default();
    state.sequences.insert("seq1".to_string(), dna.clone());
    state.sequences.insert(
        "read_query".to_string(),
        DNAsequence::from_sequence("ACGTTTAA").expect("read"),
    );
    let engine = Arc::new(RwLock::new(GentleEngine::from_state(state)));
    let mut area = MainAreaDna::new(dna, Some("seq1".to_string()), Some(engine));
    area.dotplot_query_override_seq_id = "read_query".to_string();
    area.dotplot_ui.reference_seq_id = "read_query".to_string();
    area.dotplot_ui.reference_span_start_0based = "0".to_string();
    area.dotplot_ui.reference_span_end_0based = "25768".to_string();

    assert!(area.maybe_normalize_dotplot_reference_from_current_input());
    assert_eq!(area.dotplot_ui.reference_span_start_0based, "0");
    assert_eq!(area.dotplot_ui.reference_span_end_0based, "8");
}

#[test]
fn dotplot_transcript_reference_choices_follow_active_mapping_locus() {
    let dna = DNAsequence::from_sequence("ACGTACGT").expect("sequence");
    let mut area = MainAreaDna::new(dna, Some("seq1".to_string()), None);
    area.rna_read_mapping_window_view = Some(Arc::new(SplicingExpertView {
        seq_id: "seq1".to_string(),
        target_feature_id: 17,
        scope: SplicingScopePreset::AllOverlappingAnyStrand,
        group_label: "TP73".to_string(),
        strand: "+".to_string(),
        region_start_1based: 1,
        region_end_1based: 8,
        transcript_count: 1,
        unique_exon_count: 0,
        instruction: "mapping".to_string(),
        transcripts: vec![SplicingTranscriptLane {
            transcript_feature_id: 4,
            transcript_id: "NM_005427.4".to_string(),
            label: "TP73".to_string(),
            strand: "+".to_string(),
            exons: vec![],
            exon_cds_phases: vec![],
            introns: vec![],
            has_target_feature: true,
        }],
        unique_exons: vec![],
        matrix_rows: vec![],
        boundaries: vec![],
        intron_signals: vec![],
        junctions: vec![],
        events: vec![],
    }));

    let choices = area.dotplot_transcript_reference_choices();
    assert_eq!(choices.len(), 1);
    assert_eq!(choices[0].0, 4);
    assert!(choices[0].1.contains("TP73"));
}

#[test]
fn dotplot_can_derive_annotated_transcript_reference_without_switching_query_window() {
    let dna = transcript_derivation_test_sequence();
    let mut state = ProjectState::default();
    state.sequences.insert("seq1".to_string(), dna.clone());
    let engine = Arc::new(RwLock::new(GentleEngine::from_state(state)));
    let mut area = MainAreaDna::new(dna, Some("seq1".to_string()), Some(engine.clone()));
    area.rna_read_mapping_window_view = Some(Arc::new(SplicingExpertView {
        seq_id: "seq1".to_string(),
        target_feature_id: 0,
        scope: SplicingScopePreset::AllOverlappingAnyStrand,
        group_label: "GENE1".to_string(),
        strand: "+".to_string(),
        region_start_1based: 1,
        region_end_1based: 40,
        transcript_count: 2,
        unique_exon_count: 0,
        instruction: "mapping".to_string(),
        transcripts: vec![
            SplicingTranscriptLane {
                transcript_feature_id: 0,
                transcript_id: "NM_TEST_1".to_string(),
                label: "NM_TEST_1".to_string(),
                strand: "+".to_string(),
                exons: vec![],
                exon_cds_phases: vec![],
                introns: vec![],
                has_target_feature: true,
            },
            SplicingTranscriptLane {
                transcript_feature_id: 1,
                transcript_id: "NM_TEST_2".to_string(),
                label: "NM_TEST_2".to_string(),
                strand: "+".to_string(),
                exons: vec![],
                exon_cds_phases: vec![],
                introns: vec![],
                has_target_feature: false,
            },
        ],
        unique_exons: vec![],
        matrix_rows: vec![],
        boundaries: vec![],
        intron_signals: vec![],
        junctions: vec![],
        events: vec![],
    }));

    area.use_annotated_transcript_as_dotplot_reference(0)
        .expect("derive transcript reference");

    assert_eq!(area.seq_id.as_deref(), Some("seq1"));
    assert_eq!(area.dotplot_ui.mode, DotplotMode::PairForward);
    assert_eq!(
        area.dotplot_ui.reference_seq_id,
        "seq1__mrna__f1__nm_test_1"
    );
    assert_eq!(area.dotplot_ui.reference_span_start_0based, "0");
    let guard = engine.read().expect("engine");
    let derived_len = guard
        .state()
        .sequences
        .get("seq1__mrna__f1__nm_test_1")
        .map(|dna| dna.len())
        .expect("derived transcript sequence");
    assert_eq!(
        area.dotplot_ui.reference_span_end_0based,
        derived_len.to_string()
    );
}

#[test]
fn rna_read_dotplot_svg_default_filename_mentions_parameters() {
    let hit = RnaReadInterpretationHit {
        record_index: 6,
        header_id: "read alpha".to_string(),
        ..RnaReadInterpretationHit::default()
    };
    let file_name = MainAreaDna::default_rna_read_sequence_dotplot_svg_file_name(
        "cdna_report",
        &hit,
        9,
        1,
        1,
        None,
    );
    assert!(file_name.ends_with(".svg"));
    assert!(file_name.contains("r7"));
    assert!(file_name.contains("_w9_s1_mm1_tileauto"));
}

#[test]
fn async_task_repaint_delay_scales_with_queue_pressure() {
    assert_eq!(
        MainAreaDna::async_task_repaint_delay(0, false),
        Duration::from_millis(700)
    );
    assert_eq!(
        MainAreaDna::async_task_repaint_delay(3, false),
        Duration::from_millis(120)
    );
    assert_eq!(
        MainAreaDna::async_task_repaint_delay(512, true),
        Duration::from_millis(20)
    );
}

#[test]
fn poll_rna_read_task_keeps_live_progress_out_of_dna_window_status() {
    let dna = DNAsequence::from_sequence("ACGTACGT").expect("sequence");
    let mut area = MainAreaDna::new(dna, Some("seq1".to_string()), None);
    let (_tx, rx) = mpsc::channel::<RnaReadTaskMessage>();
    area.show_rna_read_mapping_window = true;
    area.op_status = "keep this status in the DNA window".to_string();
    area.rna_read_progress = Some(RnaReadInterpretProgress {
        seq_id: "seq1".to_string(),
        reads_processed: 8,
        reads_total: 16,
        read_bases_processed: 3200,
        mean_read_length_bp: 400.0,
        median_read_length_bp: 400,
        p95_read_length_bp: 500,
        input_bytes_processed: 1024,
        input_bytes_total: 4096,
        seed_passed: 5,
        aligned: 1,
        tested_kmers: 64,
        matched_kmers: 18,
        seed_compute_ms: 0.0,
        align_compute_ms: 0.0,
        io_read_ms: 0.0,
        fasta_parse_ms: 0.0,
        normalize_compute_ms: 0.0,
        inference_compute_ms: 0.0,
        progress_emit_ms: 0.0,
        update_every_reads: 1000,
        done: false,
        bins: vec![],
        score_density_bins: vec![],
        seed_pass_score_density_bins: vec![],
        top_hits_preview: vec![],
        transition_support_rows: vec![],
        isoform_support_rows: vec![],
        mapped_exon_support_frequencies: vec![],
        mapped_junction_support_frequencies: vec![],
        mapped_isoform_support_rows: vec![],
        reads_with_transition_support: 0,
        transition_confirmations: 0,
        junction_crossing_seed_bits_indexed: 0,
        origin_class_counts: BTreeMap::new(),
    });
    area.rna_read_task = Some(RnaReadTask {
        started: Instant::now(),
        seq_id: "seq1".to_string(),
        seed_feature_id: 0,
        report_id_hint: Some("rna_report".to_string()),
        input_path: "reads.fa.gz".to_string(),
        operation_label: "Nanopore interpretation".to_string(),
        cancel_requested: Arc::new(AtomicBool::new(false)),
        receiver: Arc::new(Mutex::new(rx)),
    });

    area.poll_rna_read_task(&egui::Context::default());

    assert_eq!(area.op_status, "keep this status in the DNA window");
}

#[test]
fn completed_rna_read_task_routes_final_status_to_mapping_workspace() {
    let dna = DNAsequence::from_sequence("ACGTACGT").expect("sequence");
    let mut state = ProjectState::default();
    state.sequences.insert("seq1".to_string(), dna.clone());
    let engine = Arc::new(RwLock::new(GentleEngine::from_state(state)));
    let mut area = MainAreaDna::new(dna, Some("seq1".to_string()), Some(engine));
    area.show_rna_read_mapping_window = true;
    area.rna_read_mapping_window_view = Some(Arc::new(SplicingExpertView {
        seq_id: "seq1".to_string(),
        target_feature_id: 9,
        scope: SplicingScopePreset::AllOverlappingAnyStrand,
        group_label: "TP73".to_string(),
        strand: "+".to_string(),
        region_start_1based: 1,
        region_end_1based: 8,
        transcript_count: 1,
        unique_exon_count: 0,
        instruction: String::new(),
        transcripts: vec![],
        unique_exons: vec![],
        matrix_rows: vec![],
        boundaries: vec![],
        intron_signals: vec![],
        junctions: vec![],
        events: vec![],
    }));
    area.op_status = "keep this status in the DNA window".to_string();
    let (tx, rx) = mpsc::channel::<RnaReadTaskMessage>();
    area.rna_read_task = Some(RnaReadTask {
        started: Instant::now(),
        seq_id: "seq1".to_string(),
        seed_feature_id: 9,
        report_id_hint: Some("rna_report".to_string()),
        input_path: "reads.fa.gz".to_string(),
        operation_label: "Nanopore interpretation".to_string(),
        cancel_requested: Arc::new(AtomicBool::new(false)),
        receiver: Arc::new(Mutex::new(rx)),
    });
    tx.send(RnaReadTaskMessage::Done(Ok(RnaReadTaskOutcome::Interpret(
        RnaReadInterpretationReport {
            report_id: "rna_report".to_string(),
            seq_id: "seq1".to_string(),
            input_path: "reads.fa.gz".to_string(),
            read_count_total: 1,
            ..RnaReadInterpretationReport::default()
        },
    ))))
    .expect("send outcome");

    area.poll_rna_read_task(&egui::Context::default());

    assert_eq!(area.op_status, "keep this status in the DNA window");
    assert!(
        area.rna_read_mapping_status
            .contains("RNA-read report 'rna_report'")
    );
}

#[test]
fn commit_completed_rna_read_task_outcome_persists_report_on_ui_thread() {
    let dna = DNAsequence::from_sequence("ACGTACGT").expect("sequence");
    let mut state = ProjectState::default();
    state.sequences.insert("seq1".to_string(), dna.clone());
    let engine = Arc::new(RwLock::new(GentleEngine::from_state(state)));
    let mut area = MainAreaDna::new(dna, Some("seq1".to_string()), Some(engine.clone()));
    let report = RnaReadInterpretationReport {
        report_id: "rna_ui_commit".to_string(),
        seq_id: "seq1".to_string(),
        input_path: "reads.fa".to_string(),
        read_count_total: 1,
        ..RnaReadInterpretationReport::default()
    };
    let result = area
        .commit_completed_rna_read_task_outcome(RnaReadTaskOutcome::Interpret(report.clone()))
        .expect("commit outcome");
    assert!(
        result
            .messages
            .iter()
            .any(|msg| msg.contains("rna_ui_commit"))
    );
    let stored = engine
        .read()
        .expect("engine")
        .get_rna_read_report("rna_ui_commit")
        .expect("stored report");
    assert_eq!(stored.report_id, report.report_id);
    assert_eq!(stored.seq_id, "seq1");
}

#[test]
fn defer_feature_tree_until_interaction_sets_runtime_flag() {
    let dna = DNAsequence::from_sequence("ACGT").expect("sequence");
    let mut area = MainAreaDna::new(dna, Some("seq1".to_string()), None);
    assert!(!area.feature_tree_deferred_until_interaction);
    area.defer_feature_tree_until_interaction();
    assert!(area.feature_tree_deferred_until_interaction);
}

#[test]
fn auto_enable_deferred_feature_tree_flips_flag_once() {
    let dna = DNAsequence::from_sequence("ACGT").expect("sequence");
    let mut area = MainAreaDna::new(dna, Some("seq1".to_string()), None);
    area.defer_feature_tree_until_interaction();
    assert!(area.maybe_auto_enable_deferred_feature_tree());
    assert!(!area.feature_tree_deferred_until_interaction);
    assert!(!area.maybe_auto_enable_deferred_feature_tree());
}

#[test]
fn replace_loaded_sequence_replaces_active_sequence_content() {
    let dna = DNAsequence::from_sequence("ACGT").expect("sequence");
    let mut area = MainAreaDna::new(dna, Some("seq1".to_string()), None);
    let replacement = DNAsequence::from_sequence("TTTT").expect("replacement sequence");
    area.replace_loaded_sequence(replacement);
    let loaded = area.dna().read().expect("dna lock").get_forward_string();
    assert_eq!(loaded, "TTTT");
}

#[test]
fn support_ratio_percent_handles_zero_total() {
    assert_eq!(MainAreaDna::support_ratio_percent(5, 0), 0.0);
}

#[test]
fn format_support_fraction_formats_percent() {
    assert_eq!(MainAreaDna::format_support_fraction(3, 8), "3/8 (37.5%)");
}

#[test]
fn prepared_genome_choice_dialog_opens_for_verify_when_always_explicit() {
    let (mut area, _td, seq_id, prepared_genome_id) =
        make_area_with_unique_compatible_anchor("always_explicit");
    let opened = area
        .maybe_prompt_prepared_genome_choice(super::PendingGenomeAnchorAction::Verify { seq_id });
    assert!(
        opened,
        "always_explicit should force explicit prepared choice"
    );
    let dialog = area
        .anchor_prepared_choice_dialog
        .as_ref()
        .expect("prepared-choice dialog should exist");
    assert_eq!(dialog.options, vec![prepared_genome_id]);
    assert_eq!(dialog.selected_index, 0);
}

#[test]
fn prepared_genome_choice_dialog_stays_closed_for_verify_when_single_compatible() {
    let (mut area, _td, seq_id, _prepared_genome_id) =
        make_area_with_unique_compatible_anchor("single_compatible");
    let opened = area
        .maybe_prompt_prepared_genome_choice(super::PendingGenomeAnchorAction::Verify { seq_id });
    assert!(
        !opened,
        "single_compatible should auto-resolve unique compatible prepared genome"
    );
    assert!(area.anchor_prepared_choice_dialog.is_none());
}

#[test]
fn parse_primer_locked_positions_accepts_offset_base_entries() {
    let locks = MainAreaDna::parse_primer_locked_positions("0:A,3:R,10:T", "locks")
        .expect("locked positions");
    assert_eq!(locks.len(), 3);
    assert_eq!(locks[0].offset_0based, 0);
    assert_eq!(locks[0].base, "A");
    assert_eq!(locks[1].offset_0based, 3);
    assert_eq!(locks[1].base, "R");
    assert_eq!(locks[2].offset_0based, 10);
    assert_eq!(locks[2].base, "T");
}

#[test]
fn parse_primer_locked_positions_rejects_duplicate_offsets() {
    let err = MainAreaDna::parse_primer_locked_positions("2:A,2:G", "locks")
        .expect_err("duplicate offsets should fail");
    assert!(err.contains("duplicate offset"));
}

#[test]
fn seed_primer_roi_updates_both_primer_and_qpcr_forms() {
    let dna = DNAsequence::from_sequence(&"ACGT".repeat(100)).expect("sequence");
    let mut area = MainAreaDna::new(dna, Some("seq1".to_string()), None);
    assert!(!area.show_engine_ops);
    area.seed_primer_design_roi_0based(25, 150, "test");
    assert_eq!(area.primer_design_ui.roi_start_0based, "25");
    assert_eq!(area.primer_design_ui.roi_end_0based, "150");
    assert_eq!(area.qpcr_design_ui.roi_start_0based, "25");
    assert_eq!(area.qpcr_design_ui.roi_end_0based, "150");
    assert!(!area.show_engine_ops);
}

#[test]
fn queue_current_selection_adds_pcr_region() {
    let dna = DNAsequence::from_sequence(&"ACGT".repeat(100)).expect("sequence");
    let mut area = MainAreaDna::new(dna, Some("seq1".to_string()), None);
    area.dna_display
        .write()
        .expect("display lock")
        .select(Selection::new(15, 75, 400));

    area.queue_current_selection_for_pcr();

    assert_eq!(area.pcr_queued_regions_ui.len(), 1);
    let queued = &area.pcr_queued_regions_ui[0];
    assert_eq!(queued.template, "seq1");
    assert_eq!(queued.source_label, "current sequence selection");
    assert_eq!(queued.start_0based, 15);
    assert_eq!(queued.end_0based_exclusive, 75);
    assert_eq!(area.primer_design_ui.roi_start_0based, "0");
    assert_eq!(area.primer_design_ui.roi_end_0based, "0");
    assert!(!area.show_engine_ops);
}

#[test]
fn set_primer_roi_from_current_selection_updates_forms() {
    let dna = DNAsequence::from_sequence(&"ACGT".repeat(100)).expect("sequence");
    let mut area = MainAreaDna::new(dna, Some("seq1".to_string()), None);
    area.dna_display
        .write()
        .expect("display lock")
        .select(Selection::new(25, 145, 400));

    area.set_primer_design_roi_from_current_selection("test selection")
        .expect("selection should seed ROI");

    assert_eq!(area.primer_design_ui.roi_start_0based, "25");
    assert_eq!(area.primer_design_ui.roi_end_0based, "145");
    assert_eq!(area.qpcr_design_ui.roi_start_0based, "25");
    assert_eq!(area.qpcr_design_ui.roi_end_0based, "145");
}

#[test]
fn seed_primer_design_roi_from_current_selection_updates_forms_without_engine_ops() {
    let dna = DNAsequence::from_sequence(&"ACGT".repeat(100)).expect("sequence");
    let mut area = MainAreaDna::new(dna, Some("seq1".to_string()), None);
    area.dna_display
        .write()
        .expect("display lock")
        .select(Selection::new(25, 145, 400));

    area.seed_primer_design_roi_from_active_region()
        .expect("selection should seed PCR ROI");

    assert_eq!(area.primer_design_ui.roi_start_0based, "25");
    assert_eq!(area.primer_design_ui.roi_end_0based, "145");
    assert_eq!(area.qpcr_design_ui.roi_start_0based, "25");
    assert_eq!(area.qpcr_design_ui.roi_end_0based, "145");
    assert!(!area.show_engine_ops);
    assert!(
        area.op_status
            .contains("Seeded primer/qPCR ROI from current sequence selection")
    );
    assert!(area.op_status.contains("open PCR Designer"));
}

#[test]
fn seed_primer_design_roi_from_painted_roi_interval_updates_forms_without_engine_ops() {
    let dna = DNAsequence::from_sequence(&"ACGT".repeat(100)).expect("sequence");
    let mut area = MainAreaDna::new(dna, Some("seq1".to_string()), None);
    area.apply_pcr_painted_interval_0based(PcrPaintRole::Roi, 41, 92, false);

    area.seed_primer_design_roi_from_active_region()
        .expect("painted ROI should seed PCR ROI");

    assert_eq!(area.primer_design_ui.roi_start_0based, "41");
    assert_eq!(area.primer_design_ui.roi_end_0based, "92");
    assert_eq!(area.qpcr_design_ui.roi_start_0based, "41");
    assert_eq!(area.qpcr_design_ui.roi_end_0based, "92");
    assert!(!area.show_engine_ops);
    assert!(
        area.op_status
            .contains("Seeded primer/qPCR ROI from painted ROI interval")
    );
    assert!(area.op_status.contains("open PCR Designer"));
}

#[test]
fn apply_restriction_cloning_single_site_recommendation_sets_mode_and_clears_saved_handoff() {
    let dna = DNAsequence::from_sequence(&"ACGT".repeat(100)).expect("sequence");
    let mut area = MainAreaDna::new(dna, Some("seq1".to_string()), None);
    area.primer_design_ui
        .restriction_cloning
        .destination_vector_seq_id = "vec1".to_string();
    area.primer_design_ui
        .restriction_cloning
        .selected_saved_report_id = "handoff_old".to_string();
    area.primer_design_ui.restriction_cloning.mode = RestrictionCloningPcrHandoffMode::DirectedPair;
    area.primer_design_ui.restriction_cloning.reverse_enzyme = "HindIII".to_string();

    area.apply_restriction_cloning_single_site_recommendation("EcoRI");

    assert_eq!(
        area.primer_design_ui.restriction_cloning.mode,
        RestrictionCloningPcrHandoffMode::SingleSite
    );
    assert_eq!(
        area.primer_design_ui.restriction_cloning.forward_enzyme,
        "EcoRI"
    );
    assert_eq!(
        area.primer_design_ui.restriction_cloning.reverse_enzyme,
        "EcoRI"
    );
    assert!(
        area.primer_design_ui
            .restriction_cloning
            .selected_saved_report_id
            .is_empty()
    );
    assert!(
        area.op_status
            .contains("single-site recommendation 'EcoRI'")
    );
    assert!(area.op_status.contains("vec1"));
}

#[test]
fn apply_restriction_cloning_directed_pair_recommendation_sets_mode_and_order() {
    let dna = DNAsequence::from_sequence(&"ACGT".repeat(100)).expect("sequence");
    let mut area = MainAreaDna::new(dna, Some("seq1".to_string()), None);
    area.primer_design_ui
        .restriction_cloning
        .destination_vector_seq_id = "vec1".to_string();
    area.primer_design_ui
        .restriction_cloning
        .selected_saved_report_id = "handoff_old".to_string();

    area.apply_restriction_cloning_directed_pair_recommendation(
        "EcoRI",
        "HindIII",
        "mcs_expected_sites",
    );

    assert_eq!(
        area.primer_design_ui.restriction_cloning.mode,
        RestrictionCloningPcrHandoffMode::DirectedPair
    );
    assert_eq!(
        area.primer_design_ui.restriction_cloning.forward_enzyme,
        "EcoRI"
    );
    assert_eq!(
        area.primer_design_ui.restriction_cloning.reverse_enzyme,
        "HindIII"
    );
    assert!(
        area.primer_design_ui
            .restriction_cloning
            .selected_saved_report_id
            .is_empty()
    );
    assert!(area.op_status.contains("EcoRI -> HindIII"));
    assert!(area.op_status.contains("mcs_expected_sites"));
    assert!(area.op_status.contains("vec1"));
}

#[test]
fn seed_simple_pcr_roi_sets_flanking_and_min_amplicon() {
    let dna = DNAsequence::from_sequence(&"ACGT".repeat(100)).expect("sequence");
    let mut area = MainAreaDna::new(dna, Some("seq1".to_string()), None);

    let seeded = area
        .seed_simple_pcr_roi_0based(25, 145, "test selection")
        .expect("seed simple pcr");

    assert_eq!(seeded, (25, 145));
    assert_eq!(area.primer_design_ui.roi_start_0based, "25");
    assert_eq!(area.primer_design_ui.roi_end_0based, "145");
    assert!(area.primer_design_ui.pair_constraints.require_roi_flanking);
    assert_eq!(area.primer_design_ui.min_amplicon_bp, "120");
}

#[test]
fn apply_simple_pcr_flank_windows_uses_roi_and_distance() {
    let dna = DNAsequence::from_sequence(&"ACGT".repeat(100)).expect("sequence");
    let mut area = MainAreaDna::new(dna, Some("seq1".to_string()), None);
    area.primer_design_ui.roi_start_0based = "25".to_string();
    area.primer_design_ui.roi_end_0based = "145".to_string();
    area.simple_pcr_max_primer_distance_bp = "40".to_string();

    area.apply_simple_pcr_flank_windows_from_roi()
        .expect("apply simple pcr windows");

    assert_eq!(area.primer_design_ui.forward.start_0based, "0");
    assert_eq!(area.primer_design_ui.forward.end_0based, "25");
    assert_eq!(area.primer_design_ui.reverse.start_0based, "145");
    assert_eq!(area.primer_design_ui.reverse.end_0based, "185");
    assert!(area.primer_design_ui.pair_constraints.require_roi_flanking);
}

#[test]
fn open_simple_pcr_from_current_selection_seeds_roi_and_windows() {
    let dna = DNAsequence::from_sequence(&"ACGT".repeat(100)).expect("sequence");
    let mut area = MainAreaDna::new(dna, Some("seq1".to_string()), None);
    area.simple_pcr_max_primer_distance_bp = "30".to_string();
    area.dna_display
        .write()
        .expect("display lock")
        .select(Selection::new(25, 145, 400));

    area.open_simple_pcr_designer_from_current_selection();

    assert_eq!(area.primer_design_ui.roi_start_0based, "25");
    assert_eq!(area.primer_design_ui.roi_end_0based, "145");
    assert_eq!(area.primer_design_ui.min_amplicon_bp, "120");
    assert_eq!(area.primer_design_ui.forward.start_0based, "0");
    assert_eq!(area.primer_design_ui.forward.end_0based, "25");
    assert_eq!(area.primer_design_ui.reverse.start_0based, "145");
    assert_eq!(area.primer_design_ui.reverse.end_0based, "175");
    assert!(area.primer_design_ui.pair_constraints.require_roi_flanking);
}

#[test]
fn roi_coordinate_formula_resolves_feature_boundary_and_offset() {
    let mut dna = DNAsequence::from_sequence(&"ACGT".repeat(100)).expect("sequence");
    dna.features_mut().push(Feature {
        kind: "CDS".into(),
        location: Location::simple_range(20, 80),
        qualifiers: vec![("label".into(), Some("CDS_A".to_string()))],
    });
    let area = MainAreaDna::new(dna, Some("seq1".to_string()), None);

    let dna = area.dna.read().expect("dna");
    let start = parse_required_usize_or_formula_text_on_sequence(
        &dna,
        "=CDS.start+10",
        "primer_design.roi_start_0based",
    )
    .expect("formula start");
    let end = parse_required_usize_or_formula_text_on_sequence(
        &dna,
        "=CDS.end-5",
        "primer_design.roi_end_0based",
    )
    .expect("formula end");

    assert_eq!(start, 30);
    assert_eq!(end, 75);
}

#[test]
fn roi_range_formula_in_start_field_resolves_both_coordinates() {
    let mut dna = DNAsequence::from_sequence(&"ACGT".repeat(100)).expect("sequence");
    dna.features_mut().push(Feature {
        kind: "CDS".into(),
        location: Location::simple_range(20, 80),
        qualifiers: vec![("label".into(), Some("CDS_A".to_string()))],
    });
    let area = MainAreaDna::new(dna, Some("seq1".to_string()), None);

    let (start, end) = area
        .resolve_roi_range_inputs_0based("=CDS.start+10 .. CDS.end-5", "0", "primer_design")
        .expect("range formula");

    assert_eq!((start, end), (30, 75));
}

#[test]
fn queue_current_primer_roi_fields_supports_range_formula() {
    let mut dna = DNAsequence::from_sequence(&"ACGT".repeat(100)).expect("sequence");
    dna.features_mut().push(Feature {
        kind: "CDS".into(),
        location: Location::simple_range(20, 80),
        qualifiers: vec![("label".into(), Some("CDS_A".to_string()))],
    });
    let mut area = MainAreaDna::new(dna, Some("seq1".to_string()), None);
    area.primer_design_ui.roi_start_0based = "=CDS.start+10 .. CDS.end-5".to_string();
    area.primer_design_ui.roi_end_0based = "0".to_string();

    area.queue_current_primer_roi_fields_for_pcr();

    assert_eq!(area.pcr_queued_regions_ui.len(), 1);
    let queued = &area.pcr_queued_regions_ui[0];
    assert_eq!(queued.template, "seq1");
    assert_eq!(queued.start_0based, 30);
    assert_eq!(queued.end_0based_exclusive, 75);
    assert_eq!(area.primer_design_ui.roi_start_0based, "30");
    assert_eq!(area.primer_design_ui.roi_end_0based, "75");
}

#[test]
fn selection_formula_applies_current_selection_range() {
    let mut dna = DNAsequence::from_sequence(&"ACGT".repeat(120)).expect("sequence");
    dna.features_mut().push(Feature {
        kind: "CDS".into(),
        location: Location::simple_range(20, 100),
        qualifiers: vec![("label".into(), Some("CDS_A".to_string()))],
    });
    let mut area = MainAreaDna::new(dna, Some("seq1".to_string()), None);
    area.selection_formula_text = "=CDS.start+10 .. CDS.end-5".to_string();

    area.apply_selection_formula();

    assert_eq!(area.current_selection_range_0based(), Some((30, 95)));
}

#[test]
fn selection_formula_inline_controls_apply_button_resolves_formula() {
    let mut dna = DNAsequence::from_sequence(&"ACGT".repeat(120)).expect("sequence");
    dna.features_mut().push(Feature {
        kind: "CDS".into(),
        location: Location::simple_range(20, 100),
        qualifiers: vec![("label".into(), Some("CDS_A".to_string()))],
    });
    let mut area = MainAreaDna::new(dna, Some("seq1".to_string()), None);
    area.selection_formula_text = "=CDS.start+10 .. CDS.end-5".to_string();

    let ctx = egui::Context::default();
    let screen_rect = egui::Rect::from_min_size(egui::Pos2::ZERO, egui::vec2(720.0, 180.0));
    let rects = render_selection_formula_control_pass(
        &ctx,
        &mut area,
        egui::RawInput {
            screen_rect: Some(screen_rect),
            ..Default::default()
        },
    );
    assert!(
        rects.iter().any(|(text, _)| text == "Selection formula"),
        "inline control should render its label: {rects:?}"
    );
    let apply_center = rects
        .iter()
        .find_map(|(text, rect)| (text == "Apply Sel").then_some(rect.center()))
        .expect("Apply Sel button text should render");

    render_selection_formula_control_pass(
        &ctx,
        &mut area,
        egui::RawInput {
            screen_rect: Some(screen_rect),
            events: vec![
                egui::Event::PointerMoved(apply_center),
                egui::Event::PointerButton {
                    pos: apply_center,
                    button: egui::PointerButton::Primary,
                    pressed: true,
                    modifiers: egui::Modifiers::default(),
                },
            ],
            ..Default::default()
        },
    );
    render_selection_formula_control_pass(
        &ctx,
        &mut area,
        egui::RawInput {
            screen_rect: Some(screen_rect),
            events: vec![
                egui::Event::PointerMoved(apply_center),
                egui::Event::PointerButton {
                    pos: apply_center,
                    button: egui::PointerButton::Primary,
                    pressed: false,
                    modifiers: egui::Modifiers::default(),
                },
            ],
            ..Default::default()
        },
    );

    assert_eq!(area.current_selection_range_0based(), Some((30, 95)));
    assert!(area.op_status.contains("Selected region from formula"));
}

#[test]
fn selection_formula_inline_controls_enter_key_resolves_typed_formula() {
    let mut dna = DNAsequence::from_sequence(&"ACGT".repeat(120)).expect("sequence");
    dna.features_mut().push(Feature {
        kind: "CDS".into(),
        location: Location::simple_range(20, 100),
        qualifiers: vec![("label".into(), Some("CDS_A".to_string()))],
    });
    let mut area = MainAreaDna::new(dna, Some("seq1".to_string()), None);

    let ctx = egui::Context::default();
    let screen_rect = egui::Rect::from_min_size(egui::Pos2::ZERO, egui::vec2(720.0, 180.0));
    let rects = render_selection_formula_control_pass(
        &ctx,
        &mut area,
        egui::RawInput {
            screen_rect: Some(screen_rect),
            ..Default::default()
        },
    );
    let formula_field_center = rects
        .iter()
        .find_map(|(text, rect)| {
            text.starts_with("=CDS.start+10")
                .then_some(rect.center() + egui::vec2(8.0, 0.0))
        })
        .expect("formula hint text should render inside the text field");

    for pressed in [true, false] {
        render_selection_formula_control_pass(
            &ctx,
            &mut area,
            egui::RawInput {
                screen_rect: Some(screen_rect),
                events: vec![
                    egui::Event::PointerMoved(formula_field_center),
                    egui::Event::PointerButton {
                        pos: formula_field_center,
                        button: egui::PointerButton::Primary,
                        pressed,
                        modifiers: egui::Modifiers::default(),
                    },
                ],
                ..Default::default()
            },
        );
    }

    render_selection_formula_control_pass(
        &ctx,
        &mut area,
        egui::RawInput {
            screen_rect: Some(screen_rect),
            events: vec![egui::Event::Text("=CDS.start+10 .. CDS.end-5".to_string())],
            ..Default::default()
        },
    );
    assert_eq!(area.selection_formula_text, "=CDS.start+10 .. CDS.end-5");

    render_selection_formula_control_pass(
        &ctx,
        &mut area,
        egui::RawInput {
            screen_rect: Some(screen_rect),
            events: vec![egui::Event::Key {
                key: egui::Key::Enter,
                physical_key: Some(egui::Key::Enter),
                pressed: true,
                repeat: false,
                modifiers: egui::Modifiers::default(),
            }],
            ..Default::default()
        },
    );

    assert_eq!(area.current_selection_range_0based(), Some((30, 95)));
    assert!(area.op_status.contains("Selected region from formula"));
}

#[test]
fn selection_formula_inline_controls_enter_key_preserves_selection_on_invalid_formula() {
    let mut dna = DNAsequence::from_sequence(&"ACGT".repeat(120)).expect("sequence");
    dna.features_mut().push(Feature {
        kind: "CDS".into(),
        location: Location::simple_range(20, 100),
        qualifiers: vec![("label".into(), Some("CDS_A".to_string()))],
    });
    let mut area = MainAreaDna::new(dna, Some("seq1".to_string()), None);
    area.set_selection_range_0based(5, 15)
        .expect("initial selection");

    let ctx = egui::Context::default();
    let screen_rect = egui::Rect::from_min_size(egui::Pos2::ZERO, egui::vec2(720.0, 180.0));
    let rects = render_selection_formula_control_pass(
        &ctx,
        &mut area,
        egui::RawInput {
            screen_rect: Some(screen_rect),
            ..Default::default()
        },
    );
    let formula_field_center = rects
        .iter()
        .find_map(|(text, rect)| {
            text.starts_with("=CDS.start+10")
                .then_some(rect.center() + egui::vec2(8.0, 0.0))
        })
        .expect("formula hint text should render inside the text field");

    for pressed in [true, false] {
        render_selection_formula_control_pass(
            &ctx,
            &mut area,
            egui::RawInput {
                screen_rect: Some(screen_rect),
                events: vec![
                    egui::Event::PointerMoved(formula_field_center),
                    egui::Event::PointerButton {
                        pos: formula_field_center,
                        button: egui::PointerButton::Primary,
                        pressed,
                        modifiers: egui::Modifiers::default(),
                    },
                ],
                ..Default::default()
            },
        );
    }

    render_selection_formula_control_pass(
        &ctx,
        &mut area,
        egui::RawInput {
            screen_rect: Some(screen_rect),
            events: vec![egui::Event::Text("=CDS.sideways .. CDS.end".to_string())],
            ..Default::default()
        },
    );
    assert_eq!(area.selection_formula_text, "=CDS.sideways .. CDS.end");

    render_selection_formula_control_pass(
        &ctx,
        &mut area,
        egui::RawInput {
            screen_rect: Some(screen_rect),
            events: vec![egui::Event::Key {
                key: egui::Key::Enter,
                physical_key: Some(egui::Key::Enter),
                pressed: true,
                repeat: false,
                modifiers: egui::Modifiers::default(),
            }],
            ..Default::default()
        },
    );

    assert_eq!(area.current_selection_range_0based(), Some((5, 15)));
    assert!(area.op_status.contains("unknown boundary 'sideways'"));
}

#[test]
fn selection_formula_supports_label_filter_and_occurrence() {
    let mut dna = DNAsequence::from_sequence(&"ACGT".repeat(200)).expect("sequence");
    dna.features_mut().push(Feature {
        kind: "gene".into(),
        location: Location::simple_range(10, 40),
        qualifiers: vec![("label".into(), Some("TP73".to_string()))],
    });
    dna.features_mut().push(Feature {
        kind: "gene".into(),
        location: Location::simple_range(80, 140),
        qualifiers: vec![("label".into(), Some("TP73".to_string()))],
    });
    let area = MainAreaDna::new(dna, Some("seq1".to_string()), None);

    let dna = area.dna.read().expect("dna");
    let start = parse_required_usize_or_formula_text_on_sequence(
        &dna,
        "=gene[label=TP73].start",
        "primer_design.roi_start_0based",
    )
    .expect("label formula");
    let second_end = parse_required_usize_or_formula_text_on_sequence(
        &dna,
        "=gene[2].end",
        "primer_design.roi_end_0based",
    )
    .expect("occurrence formula");

    assert_eq!(start, 10);
    assert_eq!(second_end, 140);
}

#[test]
fn painted_role_interval_updates_selected_role_only() {
    let dna = DNAsequence::from_sequence(&"ACGT".repeat(100)).expect("sequence");
    let mut area = MainAreaDna::new(dna, Some("seq1".to_string()), None);

    area.apply_pcr_painted_interval_0based(PcrPaintRole::Roi, 20, 90, false);
    area.apply_pcr_painted_interval_0based(PcrPaintRole::UpstreamPrimerWindow, 5, 30, false);
    area.apply_pcr_painted_interval_0based(PcrPaintRole::DownstreamPrimerWindow, 90, 130, false);

    assert_eq!(area.pcr_paint_intervals.roi, Some((20, 90)));
    assert_eq!(area.pcr_paint_intervals.upstream_window, Some((5, 30)));
    assert_eq!(area.pcr_paint_intervals.downstream_window, Some((90, 130)));
}

#[test]
fn shift_drag_roi_behavior_queues_region_immediately() {
    let dna = DNAsequence::from_sequence(&"ACGT".repeat(100)).expect("sequence");
    let mut area = MainAreaDna::new(dna, Some("seq1".to_string()), None);
    area.pcr_paint_role = PcrPaintRole::Roi;

    area.finalize_pcr_map_drag_paint(PcrPaintRole::Roi, 25, 80, true);

    assert_eq!(area.pcr_queued_regions_ui.len(), 1);
    let queued = &area.pcr_queued_regions_ui[0];
    assert_eq!(queued.template, "seq1");
    assert_eq!(queued.start_0based, 25);
    assert_eq!(queued.end_0based_exclusive, 81);
    assert_eq!(queued.source_label, "painted ROI (shift+drag)");
}

#[test]
fn painted_interval_editor_applies_typed_coordinates() {
    let dna = DNAsequence::from_sequence(&"ACGT".repeat(100)).expect("sequence");
    let mut area = MainAreaDna::new(dna, Some("seq1".to_string()), None);
    area.apply_pcr_painted_interval_0based(PcrPaintRole::Roi, 20, 90, false);
    area.pcr_paint_last_drag_start_text = "31".to_string();
    area.pcr_paint_last_drag_end_text = "101".to_string();

    area.apply_pcr_post_drag_interval_editor_fields(PcrPaintRole::Roi);

    assert_eq!(area.pcr_paint_intervals.roi, Some((31, 101)));
    assert_eq!(
        area.pcr_paint_last_drag_interval,
        Some((PcrPaintRole::Roi, 31, 101))
    );
}

#[test]
fn painted_pair_geometry_requires_all_three_regions() {
    let dna = DNAsequence::from_sequence(&"ACGT".repeat(100)).expect("sequence");
    let mut area = MainAreaDna::new(dna, Some("seq1".to_string()), None);
    area.apply_pcr_painted_interval_0based(PcrPaintRole::Roi, 20, 90, false);
    area.apply_pcr_painted_interval_0based(PcrPaintRole::UpstreamPrimerWindow, 5, 30, false);
    assert!(area.painted_pair_geometry_bp().is_none());
    area.apply_pcr_painted_interval_0based(PcrPaintRole::DownstreamPrimerWindow, 90, 120, false);
    assert_eq!(area.painted_pair_geometry_bp(), Some((70, 25, 30)));
}

#[test]
fn queue_current_primer_roi_fields_adds_region_spec() {
    let dna = DNAsequence::from_sequence(&"ACGT".repeat(100)).expect("sequence");
    let mut area = MainAreaDna::new(dna, Some("seq1".to_string()), None);
    area.primer_design_ui.roi_start_0based = "25".to_string();
    area.primer_design_ui.roi_end_0based = "125".to_string();

    area.queue_current_primer_roi_fields_for_pcr();

    assert_eq!(area.pcr_queued_regions_ui.len(), 1);
    let queued = &area.pcr_queued_regions_ui[0];
    assert_eq!(queued.template, "seq1");
    assert_eq!(queued.source_label, "primer ROI form");
    assert_eq!(queued.start_0based, 25);
    assert_eq!(queued.end_0based_exclusive, 125);
    assert!(area.op_status.contains("queue stores region specs"));
}

#[test]
fn build_design_primer_pairs_operation_rejects_min_amplicon_larger_than_roi() {
    let dna = DNAsequence::from_sequence(&"ACGT".repeat(100)).expect("sequence");
    let mut area = MainAreaDna::new(dna, Some("seq1".to_string()), None);
    area.primer_design_ui.roi_start_0based = "25".to_string();
    area.primer_design_ui.roi_end_0based = "75".to_string();
    area.primer_design_ui.min_amplicon_bp = "60".to_string();
    area.primer_design_ui.max_amplicon_bp = "120".to_string();

    let err = area
        .build_design_primer_pairs_operation("seq1")
        .expect_err("min amplicon larger than ROI should be rejected");

    assert!(err.contains("must be <= ROI length (50)"));
}

#[test]
fn prepare_primer_pair_design_batch_inputs_rejects_min_amplicon_larger_than_queued_roi() {
    let dna = DNAsequence::from_sequence(&"ACGT".repeat(100)).expect("sequence");
    let mut area = MainAreaDna::new(dna, Some("seq1".to_string()), None);
    area.pcr_queued_regions_ui = vec![super::PcrQueuedRegionUiState {
        template: "seq1".to_string(),
        source_label: "short_roi".to_string(),
        start_0based: 30,
        end_0based_exclusive: 70,
    }];
    area.primer_design_ui.min_amplicon_bp = "50".to_string();
    area.primer_design_ui.max_amplicon_bp = "120".to_string();

    let err = area
        .prepare_primer_pair_design_batch_inputs()
        .expect_err("queued ROI shorter than min amplicon should fail");

    assert!(err.contains("Queued PCR region #01 'short_roi' is 40 bp long"));
    assert!(err.contains("must be <= ROI length"));
}

#[test]
fn queue_selected_features_adds_multiple_regions() {
    let mut dna = DNAsequence::from_sequence(&"A".repeat(500)).expect("sequence");
    dna.features_mut().push(Feature {
        kind: "gene".into(),
        location: Location::simple_range(25, 120),
        qualifiers: vec![("label".into(), Some("f1".to_string()))],
    });
    dna.features_mut().push(Feature {
        kind: "gene".into(),
        location: Location::simple_range(180, 260),
        qualifiers: vec![("label".into(), Some("f2".to_string()))],
    });
    let mut area = MainAreaDna::new(dna, Some("seq1".to_string()), None);
    area.focus_feature(0);
    area.toggle_feature_multi_select(1);

    area.queue_selected_features_for_pcr();

    assert_eq!(area.pcr_queued_regions_ui.len(), 2);
    assert_eq!(area.pcr_queued_regions_ui[0].template, "seq1");
    assert_eq!(
        area.pcr_queued_regions_ui[0].source_label,
        "feature n-0".to_string()
    );
    assert_eq!(area.pcr_queued_regions_ui[0].start_0based, 25);
    assert_eq!(area.pcr_queued_regions_ui[0].end_0based_exclusive, 120);
    assert_eq!(
        area.pcr_queued_regions_ui[1].source_label,
        "feature n-1".to_string()
    );
    assert_eq!(area.pcr_queued_regions_ui[1].start_0based, 180);
    assert_eq!(area.pcr_queued_regions_ui[1].end_0based_exclusive, 260);
    assert_eq!(area.primer_design_ui.roi_start_0based, "0");
    assert_eq!(area.primer_design_ui.roi_end_0based, "0");
}

#[test]
fn apply_linear_viewport_from_input_fields_sets_requested_range() {
    let dna = DNAsequence::from_sequence(&"ACGT".repeat(100)).expect("sequence");
    let mut area = MainAreaDna::new(dna, Some("seq1".to_string()), None);
    area.linear_view_start_1based_input = "51".to_string();
    area.linear_view_end_1based_input = "150".to_string();

    area.apply_linear_viewport_from_input_fields();

    let (start, span, seq_len) = area.current_linear_viewport();
    assert_eq!(seq_len, 400);
    assert_eq!(start, 50);
    assert_eq!(span, 100);
    assert!(area.op_status.contains("Set linear map view to 51..150"));
}

#[test]
fn select_current_visible_span_promotes_viewport_to_selection() {
    let dna = DNAsequence::from_sequence(&"ACGT".repeat(100)).expect("sequence");
    let mut area = MainAreaDna::new(dna, Some("seq1".to_string()), None);
    area.set_linear_viewport(40, 75);

    area.select_current_visible_span();

    assert_eq!(area.current_selection_range_0based(), Some((40, 115)));
    assert!(area.op_status.contains("Selected visible span: 40..115"));
}

#[test]
fn inspect_sequence_span_selects_and_centers_linear_view() {
    let dna = DNAsequence::from_sequence(&"ACGT".repeat(100)).expect("sequence");
    let mut area = MainAreaDna::new(dna, Some("seq1".to_string()), None);
    area.set_linear_viewport(0, 100);

    area.inspect_sequence_span_0based(220, 240, "demo site")
        .expect("inspect span");

    assert_eq!(area.current_selection_range_0based(), Some((220, 240)));
    let (start, span, seq_len) = area.current_linear_viewport();
    assert_eq!(seq_len, 400);
    assert_eq!(span, 100);
    assert_eq!(start, 180);
    assert!(area.op_status.contains("Focused demo site: 220..240"));
}

#[test]
fn inspect_sequence_span_supports_wrapped_circular_selection() {
    let mut dna = DNAsequence::from_sequence(&"ACGT".repeat(25)).expect("sequence");
    dna.set_circular(true);
    let mut area = MainAreaDna::new(dna, Some("seq1".to_string()), None);

    area.inspect_sequence_span_0based(90, 10, "wrapped site")
        .expect("inspect wrapped span");

    let selection = area
        .dna_display
        .read()
        .expect("display lock")
        .selection()
        .expect("selection");
    assert_eq!(selection.from(), 90);
    assert_eq!(selection.to(), 10);
    assert!(area.op_status.contains("wrapped span 90..10"));
}

#[test]
fn linear_drag_selection_updates_selection_formula_text() {
    let dna = DNAsequence::from_sequence(&"ACGT".repeat(100)).expect("sequence");
    let mut area = MainAreaDna::new(dna, Some("seq1".to_string()), None);

    area.update_linear_drag_selection(40, 74);

    assert_eq!(area.current_selection_range_0based(), Some((40, 75)));
    assert_eq!(area.selection_formula_text, "=40 .. 75");
}

#[test]
fn linear_selection_resize_edge_detection_tracks_visible_edges() {
    let dna = DNAsequence::from_sequence(&"ACGT".repeat(100)).expect("sequence");
    let area = MainAreaDna::new(dna, Some("seq1".to_string()), None);
    area.set_linear_viewport(0, 100);
    area.dna_display
        .write()
        .expect("display lock")
        .select(Selection::new(20, 60, 400));
    let rect = egui::Rect::from_min_size(egui::pos2(0.0, 0.0), egui::vec2(1000.0, 60.0));

    assert_eq!(
        area.linear_selection_resize_edge_at_pos(rect, egui::pos2(202.0, 20.0)),
        Some(super::LinearSelectionResizeEdge::Start)
    );
    assert_eq!(
        area.linear_selection_resize_edge_at_pos(rect, egui::pos2(598.0, 20.0)),
        Some(super::LinearSelectionResizeEdge::End)
    );
    assert_eq!(
        area.linear_selection_resize_edge_at_pos(rect, egui::pos2(400.0, 20.0)),
        None
    );
}

#[test]
fn linear_selection_resize_drag_updates_selection_formula_text() {
    let dna = DNAsequence::from_sequence(&"ACGT".repeat(100)).expect("sequence");
    let mut area = MainAreaDna::new(dna, Some("seq1".to_string()), None);
    area.dna_display
        .write()
        .expect("display lock")
        .select(Selection::new(40, 80, 400));

    area.apply_linear_selection_resize_drag_bp(
        super::LinearSelectionResizeDrag {
            edge: super::LinearSelectionResizeEdge::Start,
            fixed_start: 40,
            fixed_end_exclusive: 80,
        },
        55,
    );
    assert_eq!(area.current_selection_range_0based(), Some((55, 80)));
    assert_eq!(area.selection_formula_text, "=55 .. 80");

    area.dna_display
        .write()
        .expect("display lock")
        .select(Selection::new(40, 80, 400));
    area.apply_linear_selection_resize_drag_bp(
        super::LinearSelectionResizeDrag {
            edge: super::LinearSelectionResizeEdge::End,
            fixed_start: 40,
            fixed_end_exclusive: 80,
        },
        89,
    );
    assert_eq!(area.current_selection_range_0based(), Some((40, 90)));
    assert_eq!(area.selection_formula_text, "=40 .. 90");
}

#[test]
fn primer_design_async_worker_completes_and_applies_operation() {
    let mut area = make_primer_batch_area();
    let op = Operation::SetParameter {
        name: "max_fragments_per_container".to_string(),
        value: json!(12345),
    };
    area.start_primer_design_operation(op, "Primer-pair design");
    assert!(area.primer_design_task.is_some());

    let ctx = eframe::egui::Context::default();
    for _ in 0..300 {
        area.poll_primer_design_task(&ctx);
        if area.primer_design_task.is_none() {
            break;
        }
        std::thread::sleep(Duration::from_millis(2));
    }
    assert!(
        area.primer_design_task.is_none(),
        "primer task should finish"
    );
    let engine = area.engine.as_ref().expect("engine");
    let max_fragments = engine
        .read()
        .expect("engine lock")
        .state()
        .parameters
        .max_fragments_per_container;
    assert_eq!(max_fragments, 12345);
}

#[test]
fn protocol_cartoon_preview_svg_render_resolves_bound_oe_geometry() {
    let preview = ProtocolCartoonPreviewTelemetry {
        protocol: "pcr.oe.substitution".to_string(),
        flank_bp: 80,
        overlap_bp: 11,
        insert_bp: 72,
        bindings: pcr_oe_substitution_geometry_bindings(80, 11, 72),
    };
    let svg = MainAreaDna::render_protocol_cartoon_preview_svg(&preview)
        .expect("render protocol cartoon preview svg");
    assert!(svg.starts_with("<svg"));
    assert!(!svg.contains("Invalid protocol cartoon"));
}

#[test]
fn handle_operation_success_captures_protocol_cartoon_preview_payload() {
    let mut area = make_primer_batch_area();
    let preview = ProtocolCartoonPreviewTelemetry {
        protocol: "pcr.oe.substitution".to_string(),
        flank_bp: 64,
        overlap_bp: 29,
        insert_bp: 6,
        bindings: pcr_oe_substitution_geometry_bindings(64, 29, 6),
    };
    let interpretation = ProbeRegionEvidenceInterpretationReport {
        schema: "gentle.probe_region_evidence_interpretation.v1".to_string(),
        seq_id: "array_slice".to_string(),
        gene_label: Some("PATZ1".to_string()),
        level: "pm_probe".to_string(),
        array_feature_count: 1,
        transcript_count: 2,
        ..Default::default()
    };
    area.handle_operation_success(
        super::OpResult {
            op_id: "op-preview".to_string(),
            created_seq_ids: vec![],
            changed_seq_ids: vec![],
            warnings: vec![],
            messages: vec!["preview attached".to_string()],
            protocol_cartoon_preview: Some(preview.clone()),
            genome_annotation_projection: None,
            sequence_alignment: None,
            protein_derivation_report: None,
            reverse_translation_report: None,
            protease_digest_report: None,
            protein_residue_genomic_coordinates: None,
            exon_skip_selection_plan: None,
            exon_skip_materialization: None,
            cdna_assay_test_report: None,
            cdna_assay_product_materialization: None,
            transcript_qpcr_panel: None,
            primer_specificity_report: None,
            construct_reasoning_graph: None,
            sequencing_confirmation_report: None,
            sequencing_primer_overlay_report: None,
            sequencing_trace_import_report: None,
            sequencing_trace_record: None,
            sequencing_trace_summaries: None,
            cutrun_dataset_list: None,
            cutrun_dataset_status: None,
            cutrun_dataset_projection: None,
            cutrun_read_report: None,
            cutrun_read_report_summaries: None,
            cutrun_read_coverage_export: None,
            cutrun_regulatory_support: None,
            gene_set_resolution: None,
            gene_set_promoter_cohort: None,
            gene_set_cutrun_regulatory_support: None,
            read_acquisition_report: None,
            microarray_projection: None,
            probe_region_evidence_interpretation: Some(interpretation.clone()),
            genome_coordinate_projection: None,
            rna_read_gene_support_summary: None,
            rna_read_gene_support_audit: None,
            rna_read_target_quality_export: None,
            rna_read_batch_map_report: None,
            rna_read_isoform_preflight: None,
            tfbs_region_summary: None,
            tfbs_score_tracks: None,
            tfbs_track_similarity: None,
            multi_gene_promoter_tfbs: None,
            promoter_cohort_comparison: None,
            repeat_annotation_query: None,
            sequence_repeat_overlaps: None,
            repeat_feature_materialization: None,
            repeat_environment_cohort: None,
            window_cohort_tfbs: None,
            tfbs_hit_scan: None,
            restriction_site_scan: None,
            jaspar_remote_metadata_snapshot: None,
            jaspar_catalog_report: None,
            tf_query_resolution_report: None,
            jaspar_entry_expert_view: None,
            jaspar_registry_benchmark: None,
            jaspar_entry_presentation: None,
            sequence_context_view: None,
            sequence_context_bundle: None,
            alternative_promoter_comparison: None,
            variant_promoter_context: None,
            promoter_evidence_matrix: None,
            isoform_promoter_comparison: None,
            promoter_expression_evidence: None,
            promoter_artifact_manifest: None,
            promoter_reporter_candidates: None,
            reporter_catalog: None,
            reporter_recommendation: None,
            reporter_corpus_export: None,
            reporter_construct_handoff: None,
            uniprot_projection_audit: None,
            uniprot_projection_audit_parity: None,
            lab_assistant_instructions: None,
        },
        Instant::now(),
    );
    assert_eq!(area.last_protocol_cartoon_preview_op_id, "op-preview");
    assert_eq!(
        area.last_protocol_cartoon_preview
            .as_ref()
            .map(|row| row.protocol.as_str()),
        Some("pcr.oe.substitution")
    );
    let captured = area
        .last_protocol_cartoon_preview
        .as_ref()
        .expect("captured preview");
    assert_eq!(captured.flank_bp, 64);
    assert_eq!(captured.overlap_bp, 29);
    assert_eq!(captured.insert_bp, 6);
    let cached = area
        .cached_probe_region_interpretation
        .as_ref()
        .expect("cached probe-region interpretation");
    assert_eq!(cached.seq_id, "array_slice");
    assert_eq!(cached.gene_label.as_deref(), Some("PATZ1"));
    assert_eq!(cached.level, "pm_probe");
}

#[test]
fn main_area_dna_projects_probe_region_output_through_shared_shell_capability() {
    let mut area = make_tp73_probe_region_validation_area();
    area.probe_region_projection_contrasts = "AdTAp73alpha-AdGFP".to_string();
    area.probe_region_projection_level = "pm_probe".to_string();
    area.probe_region_projection_min_abs_logfc = "0.5".to_string();
    area.probe_region_projection_max_features = "20".to_string();
    area.probe_region_projection_clear_existing = true;

    area.project_probe_region_output_for_current_path();

    let report = area
        .cached_probe_region_projection
        .as_ref()
        .expect("cached projection report");
    assert_eq!(report.schema, "gentle.microarray_projection_report.v1");
    assert_eq!(report.seq_id, "array_slice");
    assert_eq!(report.level, "pm_probe");
    assert_eq!(
        report.projected_contrasts.as_slice(),
        ["AdTAp73alpha-AdGFP"]
    );
    assert_eq!(report.parsed_rows, 14);
    assert_eq!(report.imported_features, 14);
    assert!(area.cached_probe_region_interpretation.is_none());
    assert!(
        area.op_status
            .contains("Probe-region projection imported 14 feature(s)")
    );
    assert!(area.op_status.contains("contrasts=AdTAp73alpha-AdGFP"));
}

#[test]
fn main_area_dna_interprets_probe_region_evidence_through_shared_shell_capability() {
    let temp = tempdir().expect("tempdir");
    let report_path = temp.path().join("tp73_probe_region_interpretation.json");
    let mut area = make_tp73_probe_region_validation_area();
    area.probe_region_projection_contrasts = "AdTAp73alpha-AdGFP".to_string();
    area.probe_region_projection_level = "pm_probe".to_string();
    area.probe_region_projection_min_abs_logfc = "0.5".to_string();
    area.probe_region_projection_max_features = "20".to_string();
    area.probe_region_projection_clear_existing = true;
    area.project_probe_region_output_for_current_path();
    assert!(
        area.cached_probe_region_projection.is_some(),
        "projection must populate features before interpretation"
    );

    area.probe_region_interpretation_gene_label = "TP73".to_string();
    area.probe_region_interpretation_level = "pm_probe".to_string();
    area.probe_region_interpretation_min_abs_logfc = "0.5".to_string();
    area.probe_region_interpretation_output_path = report_path.to_string_lossy().to_string();

    area.interpret_probe_region_evidence_for_current_sequence();

    assert!(report_path.exists());
    let report = area
        .cached_probe_region_interpretation
        .as_ref()
        .expect("cached interpretation report");
    assert_eq!(
        report.schema,
        "gentle.probe_region_evidence_interpretation.v1"
    );
    assert_eq!(report.seq_id, "array_slice");
    assert_eq!(report.gene_label.as_deref(), Some("TP73"));
    assert_eq!(report.level, "pm_probe");
    assert_eq!(report.array_feature_count, 14);
    assert_eq!(report.transcript_count, 2);
    assert_eq!(report.evidence_rows.len(), 14);
    assert!(report.evidence_rows.iter().any(|row| {
        row.feature_id == "342828"
            && row.parent_feature_id.as_deref() == Some("PSR0100145780.hg.1")
            && row.start_1based == Some(3652576)
            && row.end_1based == Some(3652586)
    }));
    assert!(report.evidence_rows.iter().all(|row| {
        row.ambiguity_tags
            .iter()
            .any(|tag| tag == "multi_hit_not_assessed")
            && row
                .ambiguity_tags
                .iter()
                .any(|tag| tag == "isoform_support_not_inferred")
    }));
    assert!(area.op_status.contains("Probe-region evidence interpreted"));
    assert!(
        area.op_status
            .contains("14 array feature(s), 2 transcript model(s), 14 evidence row(s)")
    );
}

#[test]
fn main_area_dna_exports_probe_region_evidence_svg_through_shared_shell_capability() {
    let temp = tempdir().expect("tempdir");
    let report_path = temp.path().join("tp73_probe_region_interpretation.json");
    let svg_path = temp.path().join("tp73_probe_region_evidence.svg");
    let mut area = make_tp73_probe_region_validation_area();
    area.probe_region_projection_contrasts = "AdTAp73alpha-AdGFP".to_string();
    area.probe_region_projection_level = "pm_probe".to_string();
    area.probe_region_projection_min_abs_logfc = "0.5".to_string();
    area.probe_region_projection_max_features = "20".to_string();
    area.probe_region_projection_clear_existing = true;
    area.project_probe_region_output_for_current_path();

    area.probe_region_interpretation_gene_label = "TP73".to_string();
    area.probe_region_interpretation_level = "pm_probe".to_string();
    area.probe_region_interpretation_min_abs_logfc = "0.5".to_string();
    area.probe_region_interpretation_output_path = report_path.to_string_lossy().to_string();
    area.interpret_probe_region_evidence_for_current_sequence();
    assert!(report_path.exists());
    assert!(area.cached_probe_region_interpretation.is_some());
    let expected_junction_span_count = area
        .cached_probe_region_interpretation
        .as_ref()
        .expect("cached interpretation")
        .evidence_rows
        .iter()
        .flat_map(|row| &row.transcript_mappings)
        .map(|mapping| mapping.junction_spans.len())
        .sum::<usize>();

    area.probe_region_evidence_svg_output_path = svg_path.to_string_lossy().to_string();
    let command = area
        .build_probe_region_evidence_svg_command_for_current_path()
        .expect("build evidence SVG shell command");
    match command {
        ShellCommand::ArraysRenderProbeRegionEvidenceSvg { report, output } => {
            assert_eq!(report, report_path.to_string_lossy().to_string());
            assert_eq!(output, svg_path.to_string_lossy().to_string());
        }
        other => panic!("unexpected command: {other:?}"),
    }

    area.export_probe_region_evidence_svg_for_current_path();

    assert!(svg_path.exists());
    let svg = fs::read_to_string(&svg_path).expect("read evidence svg");
    assert!(svg.contains("gentle.probe_region_evidence_svg_export.v1"));
    assert!(svg.contains("class=\"junction-span\""));
    assert!(svg.contains("class=\"transcript\""));
    assert!(area.op_status.contains("Probe-region evidence SVG exported"));
    assert!(area.op_status.contains("14 evidence row(s)"));
    assert!(
        area.op_status
            .contains(&format!("{expected_junction_span_count} junction span(s)"))
    );
    assert!(
        area.op_status
            .contains("report_local_geometry_aligned_to_evidence_axis_without_full_gene_model")
    );
}

#[test]
fn queued_primer_batch_async_worker_completes_and_populates_results() {
    let mut area = make_primer_batch_area();
    area.primer_design_ui.report_id = "batch_async".to_string();
    area.pcr_queued_regions_ui = vec![
        super::PcrQueuedRegionUiState {
            template: "missing_template_a".to_string(),
            source_label: "missing_a".to_string(),
            start_0based: 20,
            end_0based_exclusive: 90,
        },
        super::PcrQueuedRegionUiState {
            template: "missing_template_b".to_string(),
            source_label: "missing_b".to_string(),
            start_0based: 30,
            end_0based_exclusive: 100,
        },
    ];

    area.start_queued_primer_pair_design_batch();
    assert!(area.primer_design_task.is_some());

    let ctx = eframe::egui::Context::default();
    for _ in 0..1_500 {
        area.poll_primer_design_task(&ctx);
        if area.primer_design_task.is_none() {
            break;
        }
        std::thread::sleep(Duration::from_millis(2));
    }

    assert!(
        area.primer_design_task.is_none(),
        "queued primer batch task should finish"
    );
    assert_eq!(area.pcr_batch_results_ui.len(), 2);
    assert!(
        area.op_status
            .contains("2 queued region(s), 0 succeeded, 2 failed")
    );
    let engine = area.engine.as_ref().expect("engine");
    let reports = engine
        .read()
        .expect("engine lock")
        .list_primer_design_reports()
        .into_iter()
        .map(|row| row.report_id)
        .collect::<Vec<_>>();
    assert!(reports.is_empty());
}

#[test]
fn execute_primer_batch_reports_progress_for_each_region() {
    let mut area = make_primer_batch_area();
    area.primer_design_ui.report_id = "batch_progress".to_string();
    area.pcr_queued_regions_ui = vec![
        super::PcrQueuedRegionUiState {
            template: "missing_template_a".to_string(),
            source_label: "missing_a".to_string(),
            start_0based: 20,
            end_0based_exclusive: 90,
        },
        super::PcrQueuedRegionUiState {
            template: "missing_template_b".to_string(),
            source_label: "missing_b".to_string(),
            start_0based: 30,
            end_0based_exclusive: 100,
        },
    ];
    let prepared = area
        .prepare_primer_pair_design_batch_inputs()
        .expect("prepared batch");
    let engine = area.engine.as_ref().expect("engine");
    let mut progress_rows: Vec<(usize, usize, String)> = vec![];
    {
        let mut guard = engine.write().expect("engine lock");
        let _outcome = MainAreaDna::execute_primer_pair_design_batch(
            &mut guard,
            &prepared.queued_regions,
            &prepared.spec,
            &prepared.report_base,
            prepared.create_copies,
            |progress| {
                progress_rows.push((
                    progress.region_index_1based,
                    progress.region_count,
                    progress.template,
                ));
            },
        );
    }
    assert_eq!(progress_rows.len(), 2);
    assert_eq!(progress_rows[0], (1, 2, "missing_template_a".to_string()));
    assert_eq!(progress_rows[1], (2, 2, "missing_template_b".to_string()));
}

fn make_primer_batch_area() -> MainAreaDna {
    let mut state = ProjectState::default();
    state.parameters.primer_design_backend = PrimerDesignBackend::Internal;
    state.sequences.insert(
        "tpl".to_string(),
        DNAsequence::from_sequence(
            "ACGTTGCATGTCAGTACGATCGTACGTAGCTAGTCGATCGTACGATCGTAGCTAGCATCGATGCTAGCTAGTACGTAGCATCGATCGTAGCTAGCATGCTAGCTAGTCGATCGATCGTACGATCG",
        )
        .expect("sequence"),
    );
    let engine = Arc::new(RwLock::new(GentleEngine::from_state(state)));
    let dna = engine
        .read()
        .expect("engine lock")
        .state()
        .sequences
        .get("tpl")
        .cloned()
        .expect("template sequence");
    let mut area = MainAreaDna::new(dna, Some("tpl".to_string()), Some(engine));
    area.primer_design_ui.forward.min_tm_c = "0".to_string();
    area.primer_design_ui.forward.max_tm_c = "100".to_string();
    area.primer_design_ui.forward.min_gc_fraction = "0.0".to_string();
    area.primer_design_ui.forward.max_gc_fraction = "1.0".to_string();
    area.primer_design_ui.forward.max_anneal_hits = "1000".to_string();
    area.primer_design_ui.reverse.min_tm_c = "0".to_string();
    area.primer_design_ui.reverse.max_tm_c = "100".to_string();
    area.primer_design_ui.reverse.min_gc_fraction = "0.0".to_string();
    area.primer_design_ui.reverse.max_gc_fraction = "1.0".to_string();
    area.primer_design_ui.reverse.max_anneal_hits = "1000".to_string();
    area.primer_design_ui.min_amplicon_bp = "40".to_string();
    area.primer_design_ui.max_amplicon_bp = "140".to_string();
    area.primer_design_ui.max_tm_delta_c = "100.0".to_string();
    area.primer_design_ui.max_pairs = "5".to_string();
    area
}

#[test]
fn queued_primer_batch_creates_one_report_per_region_with_deterministic_suffixes() {
    let mut area = make_primer_batch_area();
    area.primer_design_ui.report_id = "batch_gui".to_string();
    area.pcr_queued_regions_ui = vec![
        super::PcrQueuedRegionUiState {
            template: "tpl".to_string(),
            source_label: "feature n-0".to_string(),
            start_0based: 20,
            end_0based_exclusive: 90,
        },
        super::PcrQueuedRegionUiState {
            template: "tpl".to_string(),
            source_label: "feature n-1".to_string(),
            start_0based: 30,
            end_0based_exclusive: 100,
        },
    ];

    area.run_queued_primer_pair_design_batch();

    let engine = area.engine.as_ref().expect("engine");
    let reports = engine
        .read()
        .expect("engine lock")
        .list_primer_design_reports()
        .into_iter()
        .map(|row| row.report_id)
        .collect::<Vec<_>>();
    assert_eq!(reports, vec!["batch_gui_r01", "batch_gui_r02"]);
    assert!(
        area.op_status
            .contains("2 queued region(s), 2 succeeded, 0 failed")
    );
}

#[test]
fn queued_primer_batch_copy_mode_emits_extract_region_artifacts() {
    let mut area = make_primer_batch_area();
    area.primer_design_ui.report_id = "batch_copy".to_string();
    area.pcr_batch_create_extract_copies = true;
    area.pcr_queued_regions_ui = vec![super::PcrQueuedRegionUiState {
        template: "tpl".to_string(),
        source_label: "selection".to_string(),
        start_0based: 25,
        end_0based_exclusive: 95,
    }];

    area.run_queued_primer_pair_design_batch();

    let engine = area.engine.as_ref().expect("engine");
    let guard = engine.read().expect("engine lock");
    assert!(guard.state().sequences.contains_key("tpl_pcr_roi_1"));
    assert_eq!(
        guard
            .list_primer_design_reports()
            .into_iter()
            .map(|row| row.report_id)
            .collect::<Vec<_>>(),
        vec!["batch_copy_r01".to_string()]
    );
    assert!(
        area.op_status
            .contains("Extracted region copies: 1 succeeded, 0 failed")
    );
    assert_eq!(area.pcr_batch_results_ui.len(), 1);
    assert_eq!(area.pcr_batch_results_ui[0].report_id, "batch_copy_r01");
    assert_eq!(
        area.pcr_batch_results_ui[0].copy_seq_id.as_deref(),
        Some("tpl_pcr_roi_1")
    );
}

#[test]
fn queued_primer_batch_populates_result_rows_with_failure_details() {
    let mut area = make_primer_batch_area();
    area.primer_design_ui.report_id = "batch_result".to_string();
    area.pcr_batch_create_extract_copies = true;
    area.pcr_queued_regions_ui = vec![
        super::PcrQueuedRegionUiState {
            template: "tpl".to_string(),
            source_label: "selection".to_string(),
            start_0based: 20,
            end_0based_exclusive: 90,
        },
        super::PcrQueuedRegionUiState {
            template: "missing_template".to_string(),
            source_label: "missing".to_string(),
            start_0based: 30,
            end_0based_exclusive: 80,
        },
    ];

    area.run_queued_primer_pair_design_batch();

    assert_eq!(area.pcr_batch_results_ui.len(), 2);
    let first = &area.pcr_batch_results_ui[0];
    assert_eq!(first.region_index_1based, 1);
    assert_eq!(first.report_id, "batch_result_r01");
    assert!(first.report_error.is_none());
    assert_eq!(first.copy_seq_id.as_deref(), Some("tpl_pcr_roi_1"));
    let second = &area.pcr_batch_results_ui[1];
    assert_eq!(second.region_index_1based, 2);
    assert_eq!(second.report_id, "batch_result_r02");
    assert!(second.report_error.is_some());
    assert!(second.copy_error.is_some());
}

#[test]
fn queued_primer_batch_marks_zero_pair_reports_as_empty() {
    let mut area = make_primer_batch_area();
    area.primer_design_ui.report_id = "batch_empty".to_string();
    area.primer_design_ui.forward.fixed_5prime = "AAAAAAAAAAAAAAAAAAAA".to_string();
    area.primer_design_ui.reverse.fixed_5prime = "AAAAAAAAAAAAAAAAAAAA".to_string();
    area.pcr_queued_regions_ui = vec![super::PcrQueuedRegionUiState {
        template: "tpl".to_string(),
        source_label: "selection".to_string(),
        start_0based: 20,
        end_0based_exclusive: 90,
    }];

    area.run_queued_primer_pair_design_batch();

    assert_eq!(area.pcr_batch_results_ui.len(), 1);
    let row = &area.pcr_batch_results_ui[0];
    assert_eq!(row.report_status_label(), "report:empty");
    assert_eq!(row.report_pair_count, Some(0));
    assert!(row.report_error.is_none());
    assert!(
        row.report_note
            .as_deref()
            .unwrap_or_default()
            .contains("no accepted primer pairs")
    );
    assert!(
        area.op_status
            .contains("Reports with no accepted primer pairs: 1")
    );
}

#[test]
fn open_sequence_by_id_switches_active_sequence_from_batch_copy() {
    let mut area = make_primer_batch_area();
    area.primer_design_ui.report_id = "batch_open".to_string();
    area.pcr_batch_create_extract_copies = true;
    area.pcr_queued_regions_ui = vec![super::PcrQueuedRegionUiState {
        template: "tpl".to_string(),
        source_label: "selection".to_string(),
        start_0based: 25,
        end_0based_exclusive: 95,
    }];
    area.run_queued_primer_pair_design_batch();
    assert!(area.open_sequence_by_id("tpl_pcr_roi_1").is_ok());
    assert_eq!(area.seq_id.as_deref(), Some("tpl_pcr_roi_1"));
    assert_eq!(area.dna.read().expect("dna lock").len(), 70);
}

#[test]
fn queued_primer_batch_reports_partial_failures_deterministically() {
    let mut area = make_primer_batch_area();
    area.primer_design_ui.report_id = "batch_partial".to_string();
    area.pcr_queued_regions_ui = vec![
        super::PcrQueuedRegionUiState {
            template: "tpl".to_string(),
            source_label: "selection".to_string(),
            start_0based: 20,
            end_0based_exclusive: 90,
        },
        super::PcrQueuedRegionUiState {
            template: "missing_template".to_string(),
            source_label: "bad".to_string(),
            start_0based: 30,
            end_0based_exclusive: 80,
        },
    ];

    area.run_queued_primer_pair_design_batch();

    let engine = area.engine.as_ref().expect("engine");
    let reports = engine
        .read()
        .expect("engine lock")
        .list_primer_design_reports()
        .into_iter()
        .map(|row| row.report_id)
        .collect::<Vec<_>>();
    assert_eq!(reports, vec!["batch_partial_r01"]);
    assert!(
        area.op_status
            .contains("2 queued region(s), 1 succeeded, 1 failed")
    );
    assert!(area.op_status.contains("r02 missing_template"));
}

#[test]
fn splicing_group_roi_range_converts_to_zero_based_end_exclusive() {
    let view = SplicingExpertView {
        seq_id: "seq1".to_string(),
        target_feature_id: 7,
        scope: SplicingScopePreset::AllOverlappingAnyStrand,
        group_label: "TP73".to_string(),
        strand: "+".to_string(),
        region_start_1based: 101,
        region_end_1based: 640,
        transcript_count: 0,
        unique_exon_count: 0,
        instruction: String::new(),
        transcripts: vec![],
        unique_exons: vec![],
        matrix_rows: vec![],
        boundaries: vec![],
        intron_signals: vec![],
        junctions: vec![],
        events: vec![],
    };
    assert_eq!(
        MainAreaDna::splicing_group_roi_range_0based(&view),
        Some((100, 640))
    );
}

#[test]
fn transcript_derivation_compact_feedback_summarizes_multi_transcript_results() {
    let dna = transcript_derivation_test_sequence();
    let mut state = ProjectState::default();
    state.sequences.insert("seq1".to_string(), dna.clone());
    let engine = Arc::new(RwLock::new(GentleEngine::from_state(state)));
    let mut area = MainAreaDna::new(dna, Some("seq1".to_string()), Some(engine));

    let result = area
        .derive_transcript_sequences_with_compact_feedback(
            "seq1".to_string(),
            vec![],
            None,
            "all transcript features on 'seq1'",
        )
        .expect("derive transcripts");

    assert_eq!(result.created_seq_ids.len(), 2);
    assert_eq!(area.last_created_seq_ids.len(), 2);
    assert!(area.op_status.contains("Derived 2 transcript sequence(s)"));
    assert!(
        area.op_status
            .contains("Engine Ops pool/export targets now point to this set")
    );
    assert!(
        area.op_status
            .contains("created: seq1__mrna__f1__nm_test_1")
    );
    assert!(!area.op_status.contains("messages:"));
    assert!(
        !area
            .op_status
            .contains("Derived transcript 'seq1__mrna__f1__nm_test_1' from mRNA feature")
    );
}

#[test]
fn transcript_derivation_compact_feedback_preserves_single_result_window_switch() {
    let dna = transcript_derivation_test_sequence();
    let mut state = ProjectState::default();
    state.sequences.insert("seq1".to_string(), dna.clone());
    let engine = Arc::new(RwLock::new(GentleEngine::from_state(state)));
    let mut area = MainAreaDna::new(dna, Some("seq1".to_string()), Some(engine));

    let result = area
        .derive_transcript_sequences_with_compact_feedback(
            "seq1".to_string(),
            vec![0],
            None,
            "selected transcript lane n-1",
        )
        .expect("derive one transcript");

    let derived_seq_id = result
        .created_seq_ids
        .first()
        .cloned()
        .expect("created sequence");
    assert_eq!(area.seq_id.as_deref(), Some(derived_seq_id.as_str()));
    assert!(
        area.op_status
            .contains(&format!("Derived transcript '{derived_seq_id}'"))
    );
    assert!(
        area.op_status
            .contains("switched this window to the new cDNA sequence")
    );
    assert!(!area.op_status.contains("messages:"));
}

#[test]
fn recommend_pair_dotplot_step_respects_pair_evaluation_limit() {
    let step = MainAreaDna::recommend_pair_dotplot_step(5000, 83686, 7, 1);
    let query_windows = MainAreaDna::dotplot_window_count(5000, 7, step);
    let reference_windows = MainAreaDna::dotplot_window_count(83686, 7, step);
    assert!(
        query_windows.saturating_mul(reference_windows)
            <= crate::engine::MAX_DOTPLOT_PAIR_EVALUATIONS
    );
    assert!(step >= 1);
}

#[test]
fn splicing_transcript_dotplot_mode_uses_transcript_lane_strand() {
    let view = SplicingExpertView {
        seq_id: "seq1".to_string(),
        target_feature_id: 3,
        scope: SplicingScopePreset::AllOverlappingAnyStrand,
        group_label: "GENE1".to_string(),
        strand: "+".to_string(),
        region_start_1based: 10,
        region_end_1based: 90,
        transcript_count: 2,
        unique_exon_count: 0,
        instruction: String::new(),
        transcripts: vec![
            SplicingTranscriptLane {
                transcript_feature_id: 11,
                transcript_id: "TX_PLUS".to_string(),
                label: "TX_PLUS".to_string(),
                strand: "+".to_string(),
                exons: vec![SplicingRange {
                    start_1based: 10,
                    end_1based: 20,
                }],
                exon_cds_phases: vec![],
                introns: vec![],
                has_target_feature: false,
            },
            SplicingTranscriptLane {
                transcript_feature_id: 12,
                transcript_id: "TX_MINUS".to_string(),
                label: "TX_MINUS".to_string(),
                strand: "-".to_string(),
                exons: vec![SplicingRange {
                    start_1based: 30,
                    end_1based: 40,
                }],
                exon_cds_phases: vec![],
                introns: vec![],
                has_target_feature: false,
            },
        ],
        unique_exons: vec![],
        matrix_rows: vec![],
        boundaries: vec![],
        intron_signals: vec![],
        junctions: vec![],
        events: vec![],
    };
    assert_eq!(
        MainAreaDna::splicing_transcript_dotplot_mode(&view, 11),
        DotplotMode::PairForward
    );
    assert_eq!(
        MainAreaDna::splicing_transcript_dotplot_mode(&view, 12),
        DotplotMode::PairReverseComplement
    );
}

#[test]
fn list_primer_design_reports_updates_status_line() {
    let mut state = ProjectState::default();
    state.sequences.insert(
        "tpl".to_string(),
        DNAsequence::from_sequence(
            "ACGTTGCATGTCAGTACGATCGTACGTAGCTAGTCGATCGTACGATCGTAGCTAGCATCGATGCTAGCTAGTACGTAGCATCGATCGTAGCTAGCATGCTAGCTAGTCGATCGATCGTACGATCG",
        )
        .expect("sequence"),
    );
    let mut engine = GentleEngine::from_state(state);
    engine.state_mut().parameters.primer_design_backend = PrimerDesignBackend::Internal;
    engine
        .apply(Operation::DesignPrimerPairs {
            template: "tpl".to_string(),
            roi_start_0based: 30,
            roi_end_0based: 70,
            forward: PrimerDesignSideConstraint {
                min_length: 20,
                max_length: 20,
                location_0based: Some(5),
                start_0based: None,
                end_0based: None,
                min_tm_c: 40.0,
                max_tm_c: 90.0,
                min_gc_fraction: 0.0,
                max_gc_fraction: 1.0,
                max_anneal_hits: 10,
                ..Default::default()
            },
            reverse: PrimerDesignSideConstraint {
                min_length: 20,
                max_length: 20,
                location_0based: Some(90),
                start_0based: None,
                end_0based: None,
                min_tm_c: 40.0,
                max_tm_c: 90.0,
                min_gc_fraction: 0.0,
                max_gc_fraction: 1.0,
                max_anneal_hits: 10,
                ..Default::default()
            },
            pair_constraints: PrimerDesignPairConstraint::default(),
            min_amplicon_bp: 40,
            max_amplicon_bp: 130,
            max_tm_delta_c: Some(50.0),
            max_pairs: Some(10),
            report_id: Some("primer_ui_test".to_string()),
        })
        .expect("design primer pairs");
    let dna = engine
        .state()
        .sequences
        .get("tpl")
        .cloned()
        .expect("template sequence");
    let engine = Arc::new(RwLock::new(engine));
    let mut area = MainAreaDna::new(dna, Some("tpl".to_string()), Some(engine));
    area.list_primer_design_reports();
    assert!(area.op_status.contains("Primer reports: 1 total"));
    assert!(area.op_status.contains("primer_ui_test"));
}

#[test]
fn primer_pair_core_geometry_reports_distances_and_overlap() {
    let report = crate::engine::PrimerDesignReport {
        roi_start_0based: 30,
        roi_end_0based: 70,
        ..Default::default()
    };
    let flanking_pair = crate::engine::PrimerDesignPairRecord {
        forward: crate::engine::PrimerDesignPrimerRecord {
            end_0based_exclusive: 24,
            ..Default::default()
        },
        reverse: crate::engine::PrimerDesignPrimerRecord {
            start_0based: 79,
            ..Default::default()
        },
        ..Default::default()
    };
    let overlapping_pair = crate::engine::PrimerDesignPairRecord {
        forward: crate::engine::PrimerDesignPrimerRecord {
            end_0based_exclusive: 35,
            ..Default::default()
        },
        reverse: crate::engine::PrimerDesignPrimerRecord {
            start_0based: 63,
            ..Default::default()
        },
        ..Default::default()
    };

    let flanking = report.pair_core_geometry(&flanking_pair);
    assert_eq!(flanking.left_distance_from_core_bp, 6);
    assert_eq!(flanking.right_distance_from_core_bp, 9);
    assert_eq!(flanking.left_overlap_into_core_bp, 0);
    assert_eq!(flanking.right_overlap_into_core_bp, 0);
    assert!(flanking.flanks_core_cleanly());
    assert_eq!(flanking.left_label(), "6 bp");
    assert_eq!(flanking.right_label(), "9 bp");

    let overlapping = report.pair_core_geometry(&overlapping_pair);
    assert_eq!(overlapping.left_distance_from_core_bp, 0);
    assert_eq!(overlapping.right_distance_from_core_bp, 0);
    assert_eq!(overlapping.left_overlap_into_core_bp, 5);
    assert_eq!(overlapping.right_overlap_into_core_bp, 7);
    assert!(!overlapping.flanks_core_cleanly());
    assert_eq!(overlapping.left_label(), "overlap 5 bp");
    assert_eq!(overlapping.right_label(), "overlap 7 bp");
}

#[test]
fn show_primer_design_report_includes_core_distances() {
    let mut state = ProjectState::default();
    state.sequences.insert(
        "tpl".to_string(),
        DNAsequence::from_sequence(
            "ACGTTGCATGTCAGTACGATCGTACGTAGCTAGTCGATCGTACGATCGTAGCTAGCATCGATGCTAGCTAGTACGTAGCATCGATCGTAGCTAGCATGCTAGCTAGTCGATCGATCGTACGATCG",
        )
        .expect("sequence"),
    );
    let mut engine = GentleEngine::from_state(state);
    engine.state_mut().parameters.primer_design_backend = PrimerDesignBackend::Internal;
    engine
        .apply(Operation::DesignPrimerPairs {
            template: "tpl".to_string(),
            roi_start_0based: 30,
            roi_end_0based: 70,
            forward: PrimerDesignSideConstraint {
                min_length: 20,
                max_length: 20,
                location_0based: Some(5),
                start_0based: None,
                end_0based: None,
                min_tm_c: 40.0,
                max_tm_c: 90.0,
                min_gc_fraction: 0.0,
                max_gc_fraction: 1.0,
                max_anneal_hits: 10,
                ..Default::default()
            },
            reverse: PrimerDesignSideConstraint {
                min_length: 20,
                max_length: 20,
                location_0based: Some(90),
                start_0based: None,
                end_0based: None,
                min_tm_c: 40.0,
                max_tm_c: 90.0,
                min_gc_fraction: 0.0,
                max_gc_fraction: 1.0,
                max_anneal_hits: 10,
                ..Default::default()
            },
            pair_constraints: PrimerDesignPairConstraint::default(),
            min_amplicon_bp: 40,
            max_amplicon_bp: 150,
            max_tm_delta_c: Some(100.0),
            max_pairs: Some(10),
            report_id: Some("primer_ui_distances".to_string()),
        })
        .expect("design primer pairs");
    let dna = engine
        .state()
        .sequences
        .get("tpl")
        .cloned()
        .expect("template sequence");
    let engine = Arc::new(RwLock::new(engine));
    let mut area = MainAreaDna::new(dna, Some("tpl".to_string()), Some(engine));

    area.show_primer_design_report("primer_ui_distances");

    assert!(area.op_status.contains("left_to_core="));
    assert!(area.op_status.contains("right_to_core="));
    assert!(area.op_status.contains("flanks_core=true"));
}

#[test]
fn show_qpcr_design_report_updates_status_line() {
    let mut state = ProjectState::default();
    state.sequences.insert(
        "tpl".to_string(),
        DNAsequence::from_sequence(
            "GGGGGGGGGGGGGGGGGGGGCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAAAAATTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT",
        )
        .expect("sequence"),
    );
    let mut engine = GentleEngine::from_state(state);
    engine
        .apply(Operation::DesignQpcrAssays {
            template: "tpl".to_string(),
            roi_start_0based: 30,
            roi_end_0based: 70,
            forward: PrimerDesignSideConstraint {
                min_length: 20,
                max_length: 20,
                location_0based: Some(5),
                start_0based: None,
                end_0based: None,
                min_tm_c: 40.0,
                max_tm_c: 90.0,
                min_gc_fraction: 0.0,
                max_gc_fraction: 1.0,
                max_anneal_hits: 100,
                ..Default::default()
            },
            reverse: PrimerDesignSideConstraint {
                min_length: 20,
                max_length: 20,
                location_0based: Some(60),
                start_0based: None,
                end_0based: None,
                min_tm_c: 40.0,
                max_tm_c: 90.0,
                min_gc_fraction: 0.0,
                max_gc_fraction: 1.0,
                max_anneal_hits: 100,
                ..Default::default()
            },
            probe: PrimerDesignSideConstraint {
                min_length: 20,
                max_length: 20,
                location_0based: Some(35),
                start_0based: None,
                end_0based: None,
                min_tm_c: 40.0,
                max_tm_c: 90.0,
                min_gc_fraction: 0.0,
                max_gc_fraction: 1.0,
                max_anneal_hits: 100,
                ..Default::default()
            },
            pair_constraints: PrimerDesignPairConstraint::default(),
            min_amplicon_bp: 40,
            max_amplicon_bp: 130,
            max_tm_delta_c: Some(50.0),
            max_probe_tm_delta_c: Some(50.0),
            max_assays: Some(10),
            transcript_targeting: None,
            report_id: Some("qpcr_ui_test".to_string()),
        })
        .expect("design qpcr assays");
    let dna = engine
        .state()
        .sequences
        .get("tpl")
        .cloned()
        .expect("template sequence");
    let engine = Arc::new(RwLock::new(engine));
    let mut area = MainAreaDna::new(dna, Some("tpl".to_string()), Some(engine));
    area.show_qpcr_design_report("qpcr_ui_test");
    assert!(area.op_status.contains("qPCR report 'qpcr_ui_test'"));
    assert!(area.op_status.contains("assays="));
}

#[test]
fn build_design_qpcr_operation_genomic_mode_emits_no_transcript_targeting() {
    let dna = DNAsequence::from_sequence(&"ACGT".repeat(80)).expect("sequence");
    let mut area = MainAreaDna::new(dna, Some("tpl".to_string()), None);
    area.qpcr_design_ui.roi_start_0based = "20".to_string();
    area.qpcr_design_ui.roi_end_0based = "120".to_string();

    let op = area
        .build_design_qpcr_operation("tpl")
        .expect("build genomic qPCR op");

    match op {
        Operation::DesignQpcrAssays {
            transcript_targeting,
            ..
        } => assert!(transcript_targeting.is_none()),
        other => panic!("unexpected operation: {other:?}"),
    }
}

#[test]
fn build_design_qpcr_operation_shared_mode_uses_attached_splicing_context() {
    let dna = DNAsequence::from_sequence(&"ACGT".repeat(80)).expect("sequence");
    let mut area = MainAreaDna::new(dna, Some("tpl".to_string()), None);
    area.qpcr_design_ui.roi_start_0based = "20".to_string();
    area.qpcr_design_ui.roi_end_0based = "120".to_string();
    area.qpcr_design_ui.transcript_targeting.assay_intent =
        QpcrTranscriptIntentUiMode::SharedAcrossTranscripts;
    area.qpcr_design_ui.transcript_targeting.source_seq_id = "tpl".to_string();
    area.qpcr_design_ui.transcript_targeting.source_feature_id = Some(17);
    area.qpcr_design_ui.transcript_targeting.group_label = "TP73-AS3".to_string();
    area.qpcr_design_ui.transcript_targeting.transcript_count = 3;
    area.qpcr_design_ui.transcript_targeting.strand = "-".to_string();

    let op = area
        .build_design_qpcr_operation("tpl")
        .expect("build shared-transcript qPCR op");

    match op {
        Operation::DesignQpcrAssays {
            transcript_targeting,
            ..
        } => {
            let targeting = transcript_targeting.expect("shared transcript targeting");
            assert_eq!(targeting.source_feature_id, 17);
            assert_eq!(targeting.mode, QpcrTranscriptTargetingMode::SharedGene);
            assert!(targeting.transcript_id.is_none());
            assert!(targeting.specificity_evidence.is_none());
        }
        other => panic!("unexpected operation: {other:?}"),
    }
}

#[test]
fn build_design_qpcr_operation_specific_mode_emits_specificity_evidence() {
    let dna = DNAsequence::from_sequence(&"ACGT".repeat(80)).expect("sequence");
    let mut area = MainAreaDna::new(dna, Some("tpl".to_string()), None);
    area.qpcr_design_ui.roi_start_0based = "20".to_string();
    area.qpcr_design_ui.roi_end_0based = "120".to_string();
    area.qpcr_design_ui.transcript_targeting.assay_intent =
        QpcrTranscriptIntentUiMode::SpecificTranscript;
    area.qpcr_design_ui.transcript_targeting.source_seq_id = "tpl".to_string();
    area.qpcr_design_ui.transcript_targeting.source_feature_id = Some(17);
    area.qpcr_design_ui.transcript_targeting.group_label = "TP73-AS3".to_string();
    area.qpcr_design_ui.transcript_targeting.transcript_count = 3;
    area.qpcr_design_ui.transcript_targeting.strand = "-".to_string();
    area.qpcr_design_ui
        .transcript_targeting
        .selected_transcript_id = "NR_187362.1".to_string();
    area.qpcr_design_ui
        .transcript_targeting
        .selected_transcript_label = "TP73 antisense RNA 3, transcript variant 1".to_string();
    area.qpcr_design_ui
        .transcript_targeting
        .specificity_evidence = QpcrTranscriptSpecificityEvidence::UniqueExonOrChain;

    let op = area
        .build_design_qpcr_operation("tpl")
        .expect("build transcript-specific qPCR op");

    match op {
        Operation::DesignQpcrAssays {
            transcript_targeting,
            ..
        } => {
            let targeting = transcript_targeting.expect("specific transcript targeting");
            assert_eq!(targeting.source_feature_id, 17);
            assert_eq!(
                targeting.mode,
                QpcrTranscriptTargetingMode::DistinguishTranscript
            );
            assert_eq!(targeting.transcript_id.as_deref(), Some("NR_187362.1"));
            assert_eq!(
                targeting.specificity_evidence,
                Some(QpcrTranscriptSpecificityEvidence::UniqueExonOrChain)
            );
        }
        other => panic!("unexpected operation: {other:?}"),
    }
}

#[test]
fn validate_qpcr_transcript_targeting_requires_attached_splicing_context() {
    let dna = DNAsequence::from_sequence(&"ACGT".repeat(80)).expect("sequence");
    let mut area = MainAreaDna::new(dna, Some("tpl".to_string()), None);
    area.qpcr_design_ui.transcript_targeting.assay_intent =
        QpcrTranscriptIntentUiMode::SharedAcrossTranscripts;

    let err = area
        .validate_qpcr_transcript_targeting_ui("tpl")
        .expect_err("missing splicing context should fail");

    assert!(err.contains("splicing overview"));
}

#[test]
fn validate_qpcr_transcript_targeting_requires_selected_transcript() {
    let dna = DNAsequence::from_sequence(&"ACGT".repeat(80)).expect("sequence");
    let mut area = MainAreaDna::new(dna, Some("tpl".to_string()), None);
    area.qpcr_design_ui.transcript_targeting.assay_intent =
        QpcrTranscriptIntentUiMode::SpecificTranscript;
    area.qpcr_design_ui.transcript_targeting.source_seq_id = "tpl".to_string();
    area.qpcr_design_ui.transcript_targeting.source_feature_id = Some(17);
    area.qpcr_design_ui.transcript_targeting.group_label = "TP73-AS3".to_string();
    area.qpcr_design_ui.transcript_targeting.transcript_count = 3;
    area.qpcr_design_ui.transcript_targeting.strand = "-".to_string();

    let err = area
        .validate_qpcr_transcript_targeting_ui("tpl")
        .expect_err("missing transcript selection should fail");

    assert!(err.contains("Select a transcript"));
}

#[test]
fn qpcr_report_targeting_summary_and_ui_sync_preserve_specificity_details() {
    let dna = DNAsequence::from_sequence("ACGT").expect("dna");
    let mut area = MainAreaDna::new(dna, Some("tp73".to_string()), None);
    let report = crate::engine::QpcrDesignReport {
        report_id: "qpcr_targeting".to_string(),
        template: "tp73".to_string(),
        transcript_targeting: Some(QpcrTranscriptTargeting {
            source_feature_id: 17,
            mode: QpcrTranscriptTargetingMode::DistinguishTranscript,
            transcript_id: Some("NR_187362.1".to_string()),
            specificity_evidence: Some(QpcrTranscriptSpecificityEvidence::EitherPreferJunction),
        }),
        transcript_targeting_result: Some(crate::engine::QpcrTranscriptTargetingResult {
            source_feature_id: 17,
            mode: QpcrTranscriptTargetingMode::DistinguishTranscript,
            group_label: "TP73-AS3".to_string(),
            strand: "-".to_string(),
            transcript_count_considered: 3,
            transcript_id: Some("NR_187362.1".to_string()),
            transcript_label: Some("TP73 antisense RNA 3, transcript variant 1".to_string()),
            realized_specificity_evidence: Some(QpcrTranscriptSpecificityEvidence::JunctionOnly),
            selected_support_transcript_count: 1,
            selected_support_transcript_fraction: 1.0 / 3.0,
            used_shared_support_fallback: false,
            warnings: vec![],
        }),
        ..Default::default()
    };

    let summary = MainAreaDna::qpcr_report_targeting_summary(&report).expect("targeting summary");
    assert!(summary.contains("requested=specific transcript"));
    assert!(summary.contains("requested evidence=Either (prefer junction)"));
    assert!(summary.contains("realized evidence=Junction only"));

    area.sync_qpcr_transcript_targeting_ui_from_report(&report);
    assert_eq!(
        area.qpcr_design_ui.transcript_targeting.assay_intent,
        QpcrTranscriptIntentUiMode::SpecificTranscript
    );
    assert_eq!(
        area.qpcr_design_ui
            .transcript_targeting
            .selected_transcript_id,
        "NR_187362.1"
    );
    assert_eq!(
        area.qpcr_design_ui
            .transcript_targeting
            .specificity_evidence,
        QpcrTranscriptSpecificityEvidence::EitherPreferJunction
    );
    assert_eq!(
        area.qpcr_design_ui.transcript_targeting.group_label,
        "TP73-AS3"
    );
}

#[test]
fn qpcr_preview_geometry_from_assay_tracks_top_assay_lengths() {
    let report = crate::engine::QpcrDesignReport {
        report_id: "qpcr_preview_geometry".to_string(),
        template: "tpl".to_string(),
        roi_start_0based: 30,
        roi_end_0based: 70,
        assay_count: 1,
        assays: vec![crate::engine::QpcrAssayRecord {
            rank: 1,
            score: 97.2,
            forward: crate::engine::PrimerDesignPrimerRecord {
                start_0based: 8,
                end_0based_exclusive: 28,
                length_bp: 20,
                ..Default::default()
            },
            reverse: crate::engine::PrimerDesignPrimerRecord {
                start_0based: 76,
                end_0based_exclusive: 98,
                length_bp: 22,
                ..Default::default()
            },
            probe: crate::engine::PrimerDesignPrimerRecord {
                start_0based: 42,
                end_0based_exclusive: 60,
                length_bp: 18,
                ..Default::default()
            },
            amplicon_start_0based: 8,
            amplicon_end_0based_exclusive: 98,
            amplicon_length_bp: 90,
            primer_tm_delta_c: 0.8,
            probe_tm_delta_c: 4.3,
            transcript_context: None,
            rule_flags: Default::default(),
        }],
        ..Default::default()
    };
    let top = report.assays.first().expect("top assay");

    let geometry = MainAreaDna::qpcr_preview_geometry_from_assay(&report, top)
        .expect("preview geometry from assay");

    assert!(geometry.source_label.contains("qpcr_preview_geometry"));
    assert_eq!(geometry.forward_primer_site_bp, top.forward.length_bp);
    assert_eq!(geometry.reverse_primer_site_bp, top.reverse.length_bp);
    assert_eq!(geometry.probe_site_bp, top.probe.length_bp);
    assert!(geometry.bindings_feature_override_count() > 0);
}

#[test]
fn selected_qpcr_assay_rank_defaults_to_first_for_invalid_rank() {
    let dna = DNAsequence::from_sequence("ACGTACGTACGT").expect("dna");
    let mut area = MainAreaDna::new(dna, Some("tpl".to_string()), None);
    area.qpcr_design_ui.selected_assay_rank_1based = "99".to_string();
    let report = crate::engine::QpcrDesignReport {
        report_id: "qpcr_preview_geometry".to_string(),
        template: "tpl".to_string(),
        assays: vec![
            crate::engine::QpcrAssayRecord {
                rank: 1,
                ..Default::default()
            },
            crate::engine::QpcrAssayRecord {
                rank: 2,
                ..Default::default()
            },
        ],
        ..Default::default()
    };

    let rank = area
        .selected_qpcr_assay_rank_1based_for_report(&report)
        .expect("normalized selected assay rank");

    assert_eq!(rank, 1);
    assert_eq!(area.qpcr_design_ui.selected_assay_rank_1based, "1");
}

#[test]
fn qpcr_splicing_context_summary_flags_junction_crossing_probe() {
    let view = SplicingExpertView {
        seq_id: "seq".to_string(),
        target_feature_id: 1,
        scope: SplicingScopePreset::AllOverlappingAnyStrand,
        group_label: "GENE1".to_string(),
        strand: "+".to_string(),
        region_start_1based: 1,
        region_end_1based: 80,
        transcript_count: 1,
        unique_exon_count: 2,
        instruction: String::new(),
        transcripts: vec![],
        unique_exons: vec![
            SplicingExonSummary {
                start_1based: 1,
                end_1based: 20,
                support_transcript_count: 1,
                constitutive: true,
            },
            SplicingExonSummary {
                start_1based: 31,
                end_1based: 60,
                support_transcript_count: 1,
                constitutive: true,
            },
        ],
        matrix_rows: vec![],
        boundaries: vec![],
        intron_signals: vec![],
        junctions: vec![SplicingJunctionArc {
            donor_1based: 20,
            acceptor_1based: 31,
            support_transcript_count: 1,
            transcript_feature_ids: vec![11],
        }],
        events: vec![],
    };
    let assay = crate::engine::QpcrAssayRecord {
        rank: 1,
        forward: crate::engine::PrimerDesignPrimerRecord {
            start_0based: 5,
            end_0based_exclusive: 15,
            length_bp: 10,
            ..Default::default()
        },
        probe: crate::engine::PrimerDesignPrimerRecord {
            start_0based: 17,
            end_0based_exclusive: 31,
            length_bp: 14,
            ..Default::default()
        },
        reverse: crate::engine::PrimerDesignPrimerRecord {
            start_0based: 40,
            end_0based_exclusive: 50,
            length_bp: 10,
            ..Default::default()
        },
        amplicon_start_0based: 5,
        amplicon_end_0based_exclusive: 50,
        amplicon_length_bp: 45,
        ..Default::default()
    };

    let summary = MainAreaDna::qpcr_splicing_context_summary(&view, &assay)
        .expect("splicing-aware qPCR summary");

    assert_eq!(summary.assay_class_label, "junction-crossing probe");
    assert!(
        summary
            .explanation
            .contains("probe overlaps a modeled exon junction")
    );
    assert_eq!(summary.covered_junction_labels, vec!["20→31".to_string()]);
}

#[test]
fn qpcr_persisted_context_summary_preserves_transcript_aware_details() {
    let assay = crate::engine::QpcrAssayRecord {
        rank: 1,
        transcript_context: Some(crate::engine::QpcrTranscriptAssayContext {
            assay_class_label: "distinguishing junction primer".to_string(),
            explanation: "Forward primer spans a transcript-specific exon junction.".to_string(),
            probe_placement: "inside_shared_exon".to_string(),
            design_transcript_feature_id: 11,
            design_transcript_id: "NR_187362.1".to_string(),
            design_transcript_label: "TP73-AS3 variant 1".to_string(),
            support_transcript_count: 1,
            support_transcript_fraction: 1.0 / 3.0,
            supported_transcript_ids: vec!["NR_187362.1".to_string()],
            forward_source_ranges_0based: vec![],
            reverse_source_ranges_0based: vec![],
            probe_source_ranges_0based: vec![],
            amplicon_source_ranges_0based: vec![],
            covered_junction_labels: vec!["120→241".to_string()],
            forward_spans_junction: true,
            reverse_spans_junction: false,
            probe_spans_junction: false,
            genomic_equivalent_start_0based: None,
            genomic_equivalent_end_0based_exclusive: None,
            genomic_equivalent_length_bp: None,
            genomic_carryover_risk: "low".to_string(),
            genomic_carryover_rationale: "Forward primer spans a junction.".to_string(),
            transcript_distinguishing_primer: Some("forward".to_string()),
            realized_specificity_evidence: Some(QpcrTranscriptSpecificityEvidence::JunctionOnly),
            satisfies_requested_targeting: true,
        }),
        ..Default::default()
    };

    let summary = MainAreaDna::qpcr_persisted_context_summary(&assay)
        .expect("persisted transcript-aware summary");

    assert_eq!(summary.assay_class_label, "distinguishing junction primer");
    assert!(
        summary
            .explanation
            .contains("Forward primer spans a transcript-specific exon junction.")
    );
    assert!(
        summary
            .explanation
            .contains("Genomic-DNA carryover risk: low.")
    );
    assert!(
        summary
            .explanation
            .contains("Supported transcripts: 1 (33% of considered set).")
    );
    assert!(
        summary
            .explanation
            .contains("Design transcript: NR_187362.1.")
    );
    assert!(
        summary
            .explanation
            .contains("Distinguishing primer(s): forward.")
    );
    assert!(
        summary
            .explanation
            .contains("Realized specificity evidence: Junction only.")
    );
    assert_eq!(summary.covered_junction_labels, vec!["120→241".to_string()]);
}

#[test]
fn qpcr_preview_geometry_from_ui_reflects_current_lengths() {
    let dna = DNAsequence::from_sequence(
        "GGGGGGGGGGGGGGGGGGGGCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAAAAATTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT",
    )
    .expect("sequence");
    let mut area = MainAreaDna::new(dna, Some("tpl".to_string()), None);
    area.qpcr_design_ui.roi_start_0based = "25".to_string();
    area.qpcr_design_ui.roi_end_0based = "95".to_string();
    area.qpcr_design_ui.forward.min_length = "19".to_string();
    area.qpcr_design_ui.reverse.min_length = "21".to_string();
    area.qpcr_design_ui.probe.min_length = "17".to_string();

    let geometry = area
        .qpcr_preview_geometry_from_ui()
        .expect("preview geometry from ui");

    assert_eq!(geometry.forward_primer_site_bp, 19);
    assert_eq!(geometry.reverse_primer_site_bp, 21);
    assert_eq!(geometry.probe_site_bp, 17);
    assert!(geometry.total_roi_bp() >= 17);
}

#[test]
fn focus_primer_design_report_selects_report_and_shows_engine_ops() {
    let mut state = ProjectState::default();
    state.sequences.insert(
        "tpl".to_string(),
        DNAsequence::from_sequence(
            "GGGGGGGGGGGGGGGGGGGGCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAAAAATTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT",
        )
        .expect("sequence"),
    );
    let mut engine = GentleEngine::from_state(state);
    engine
        .apply(Operation::DesignPrimerPairs {
            template: "tpl".to_string(),
            roi_start_0based: 30,
            roi_end_0based: 70,
            forward: PrimerDesignSideConstraint {
                min_length: 20,
                max_length: 20,
                location_0based: Some(5),
                start_0based: None,
                end_0based: None,
                min_tm_c: 40.0,
                max_tm_c: 90.0,
                min_gc_fraction: 0.0,
                max_gc_fraction: 1.0,
                max_anneal_hits: 10,
                ..Default::default()
            },
            reverse: PrimerDesignSideConstraint {
                min_length: 20,
                max_length: 20,
                location_0based: Some(60),
                start_0based: None,
                end_0based: None,
                min_tm_c: 40.0,
                max_tm_c: 90.0,
                min_gc_fraction: 0.0,
                max_gc_fraction: 1.0,
                max_anneal_hits: 10,
                ..Default::default()
            },
            pair_constraints: PrimerDesignPairConstraint::default(),
            min_amplicon_bp: 40,
            max_amplicon_bp: 130,
            max_tm_delta_c: Some(50.0),
            max_pairs: Some(10),
            report_id: Some("primer_ui_focus".to_string()),
        })
        .expect("design primer pairs");
    let dna = engine
        .state()
        .sequences
        .get("tpl")
        .cloned()
        .expect("template sequence");
    let engine = Arc::new(RwLock::new(engine));
    let mut area = MainAreaDna::new(dna, Some("tpl".to_string()), Some(engine));

    area.focus_primer_design_report("primer_ui_focus");

    assert!(area.show_engine_ops);
    assert_eq!(area.primer_design_ui.report_id, "primer_ui_focus");
    assert_eq!(area.pcr_designer_mode, PcrDesignerMode::PrimerPairs);
    assert!(area.op_status.contains("Primer report 'primer_ui_focus'"));
}

#[test]
fn focus_qpcr_design_report_selects_report_and_shows_engine_ops() {
    let mut state = ProjectState::default();
    state.sequences.insert(
        "tpl".to_string(),
        DNAsequence::from_sequence(
            "GGGGGGGGGGGGGGGGGGGGCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAAAAATTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT",
        )
        .expect("sequence"),
    );
    let mut engine = GentleEngine::from_state(state);
    engine
        .apply(Operation::DesignQpcrAssays {
            template: "tpl".to_string(),
            roi_start_0based: 30,
            roi_end_0based: 70,
            forward: PrimerDesignSideConstraint {
                min_length: 20,
                max_length: 20,
                location_0based: Some(5),
                start_0based: None,
                end_0based: None,
                min_tm_c: 40.0,
                max_tm_c: 90.0,
                min_gc_fraction: 0.0,
                max_gc_fraction: 1.0,
                max_anneal_hits: 100,
                ..Default::default()
            },
            reverse: PrimerDesignSideConstraint {
                min_length: 20,
                max_length: 20,
                location_0based: Some(60),
                start_0based: None,
                end_0based: None,
                min_tm_c: 40.0,
                max_tm_c: 90.0,
                min_gc_fraction: 0.0,
                max_gc_fraction: 1.0,
                max_anneal_hits: 100,
                ..Default::default()
            },
            probe: PrimerDesignSideConstraint {
                min_length: 20,
                max_length: 20,
                location_0based: Some(35),
                start_0based: None,
                end_0based: None,
                min_tm_c: 40.0,
                max_tm_c: 90.0,
                min_gc_fraction: 0.0,
                max_gc_fraction: 1.0,
                max_anneal_hits: 100,
                ..Default::default()
            },
            pair_constraints: PrimerDesignPairConstraint::default(),
            min_amplicon_bp: 40,
            max_amplicon_bp: 130,
            max_tm_delta_c: Some(50.0),
            max_probe_tm_delta_c: Some(50.0),
            max_assays: Some(10),
            transcript_targeting: None,
            report_id: Some("qpcr_ui_focus".to_string()),
        })
        .expect("design qpcr assays");
    let dna = engine
        .state()
        .sequences
        .get("tpl")
        .cloned()
        .expect("template sequence");
    let engine = Arc::new(RwLock::new(engine));
    let mut area = MainAreaDna::new(dna, Some("tpl".to_string()), Some(engine));
    area.qpcr_design_ui.selected_assay_rank_1based = "7".to_string();

    area.focus_qpcr_design_report("qpcr_ui_focus");

    assert!(area.show_engine_ops);
    assert_eq!(area.qpcr_design_ui.report_id, "qpcr_ui_focus");
    assert_eq!(area.qpcr_design_ui.selected_assay_rank_1based, "1");
    assert_eq!(area.pcr_designer_mode, PcrDesignerMode::QpcrAssays);
    assert!(area.op_status.contains("qPCR report 'qpcr_ui_focus'"));
}

#[test]
fn focus_restriction_cloning_handoff_report_selects_saved_handoff_and_updates_pcr_designer() {
    let mut state = ProjectState::default();
    state.sequences.insert(
        "tpl".to_string(),
        DNAsequence::from_sequence(
            "ACGTTGCATGTCAGTACGATCGTACGTAGCTAGTCGATCGTACGATCGTAGCTAGCATCGATGCTAGCTAGTACGTAGCATCGATCGTAGCTAGCATGCTAGCTAGTCGATCGATCGTACGATCG",
        )
        .expect("sequence"),
    );
    let mut vector = DNAsequence::from_sequence("AAAAGAATTCGGGGGAAGCTTTTTT").expect("vector");
    *vector.restriction_enzymes_mut() = active_restriction_enzymes();
    let vector_len_i64 = vector.len().try_into().unwrap();
    vector.features_mut().push(gb_io::seq::Feature {
        kind: "misc_feature".into(),
        location: gb_io::seq::Location::simple_range(0, vector_len_i64),
        qualifiers: vec![
            ("label".into(), Some("MCS".to_string())),
            (
                "mcs_expected_sites".into(),
                Some("EcoRI,HindIII".to_string()),
            ),
        ],
    });
    vector.update_computed_features();
    state.sequences.insert("vec".to_string(), vector);
    let mut engine = GentleEngine::from_state(state);
    engine.state_mut().parameters.primer_design_backend = PrimerDesignBackend::Internal;
    engine
        .apply(Operation::DesignPrimerPairs {
            template: "tpl".to_string(),
            roi_start_0based: 40,
            roi_end_0based: 80,
            forward: PrimerDesignSideConstraint {
                min_length: 20,
                max_length: 20,
                location_0based: Some(5),
                start_0based: None,
                end_0based: None,
                min_tm_c: 0.0,
                max_tm_c: 100.0,
                min_gc_fraction: 0.0,
                max_gc_fraction: 1.0,
                max_anneal_hits: 1000,
                ..Default::default()
            },
            reverse: PrimerDesignSideConstraint {
                min_length: 20,
                max_length: 20,
                location_0based: Some(90),
                start_0based: None,
                end_0based: None,
                min_tm_c: 0.0,
                max_tm_c: 100.0,
                min_gc_fraction: 0.0,
                max_gc_fraction: 1.0,
                max_anneal_hits: 1000,
                ..Default::default()
            },
            pair_constraints: PrimerDesignPairConstraint::default(),
            min_amplicon_bp: 40,
            max_amplicon_bp: 150,
            max_tm_delta_c: Some(100.0),
            max_pairs: Some(10),
            report_id: Some("primer_ui_handoff".to_string()),
        })
        .expect("design primer pairs");
    engine
        .apply(Operation::PrepareRestrictionCloningPcrHandoff {
            template: "tpl".to_string(),
            primer_report_id: "primer_ui_handoff".to_string(),
            pair_index: 0,
            destination_vector_seq_id: "vec".to_string(),
            mode: RestrictionCloningPcrHandoffMode::DirectedPair,
            forward_enzyme: "EcoRI".to_string(),
            reverse_enzyme: Some("HindIII".to_string()),
            forward_leader_5prime: Some("GC".to_string()),
            reverse_leader_5prime: Some("AT".to_string()),
        })
        .expect("prepare restriction-cloning handoff");
    let handoff_id = engine
        .list_restriction_cloning_pcr_handoffs()
        .into_iter()
        .find(|row| row.primer_report_id == "primer_ui_handoff")
        .expect("handoff summary")
        .report_id;
    let dna = engine
        .state()
        .sequences
        .get("tpl")
        .cloned()
        .expect("template sequence");
    let engine = Arc::new(RwLock::new(engine));
    let mut area = MainAreaDna::new(dna, Some("tpl".to_string()), Some(engine));

    area.focus_restriction_cloning_handoff_report(&handoff_id);

    assert!(area.show_engine_ops);
    assert_eq!(area.primer_design_ui.report_id, "primer_ui_handoff");
    assert_eq!(
        area.primer_design_ui
            .restriction_cloning
            .destination_vector_seq_id,
        "vec"
    );
    assert_eq!(
        area.primer_design_ui.restriction_cloning.mode,
        RestrictionCloningPcrHandoffMode::DirectedPair
    );
    assert_eq!(
        area.primer_design_ui.restriction_cloning.forward_enzyme,
        "EcoRI"
    );
    assert_eq!(
        area.primer_design_ui.restriction_cloning.reverse_enzyme,
        "HindIII"
    );
    assert_eq!(
        area.primer_design_ui
            .restriction_cloning
            .selected_saved_report_id,
        handoff_id
    );
    assert!(area.op_status.contains("Restriction-cloning handoff"));
}

#[test]
fn seed_restriction_cloning_handoff_ui_uses_engine_recommendation_when_enzymes_blank() {
    let mut state = ProjectState::default();
    state.sequences.insert(
        "tpl".to_string(),
        DNAsequence::from_sequence(
            "ACGTTGCATGTCAGTACGATCGTACGTAGCTAGTCGATCGTACGATCGTAGCTAGCATCGATGCTAGCTAGTACGTAGCATCGATCGTAGCTAGCATGCTAGCTAGTCGATCGATCGTACGATCG",
        )
        .expect("sequence"),
    );
    let mut vector =
        DNAsequence::from_sequence("AAAAGAATTCGGGGGAAGCTTTTTTGCGGCCGCTTTT").expect("vector");
    *vector.restriction_enzymes_mut() = active_restriction_enzymes();
    let vector_len_i64 = vector.len().try_into().unwrap();
    vector.features_mut().push(gb_io::seq::Feature {
        kind: "misc_feature".into(),
        location: gb_io::seq::Location::simple_range(0, vector_len_i64),
        qualifiers: vec![
            ("label".into(), Some("MCS".to_string())),
            (
                "mcs_expected_sites".into(),
                Some("EcoRI,HindIII".to_string()),
            ),
        ],
    });
    vector.update_computed_features();
    state.sequences.insert("vec".to_string(), vector);
    let mut engine = GentleEngine::from_state(state);
    engine.state_mut().parameters.primer_design_backend = PrimerDesignBackend::Internal;
    engine
        .apply(Operation::DesignPrimerPairs {
            template: "tpl".to_string(),
            roi_start_0based: 40,
            roi_end_0based: 80,
            forward: PrimerDesignSideConstraint {
                min_length: 20,
                max_length: 20,
                location_0based: Some(5),
                start_0based: None,
                end_0based: None,
                min_tm_c: 0.0,
                max_tm_c: 100.0,
                min_gc_fraction: 0.0,
                max_gc_fraction: 1.0,
                max_anneal_hits: 1000,
                ..Default::default()
            },
            reverse: PrimerDesignSideConstraint {
                min_length: 20,
                max_length: 20,
                location_0based: Some(90),
                start_0based: None,
                end_0based: None,
                min_tm_c: 0.0,
                max_tm_c: 100.0,
                min_gc_fraction: 0.0,
                max_gc_fraction: 1.0,
                max_anneal_hits: 1000,
                ..Default::default()
            },
            pair_constraints: PrimerDesignPairConstraint::default(),
            min_amplicon_bp: 40,
            max_amplicon_bp: 150,
            max_tm_delta_c: Some(100.0),
            max_pairs: Some(10),
            report_id: Some("primer_ui_seed".to_string()),
        })
        .expect("design primer pairs");
    let dna = engine
        .state()
        .sequences
        .get("tpl")
        .cloned()
        .expect("template sequence");
    let engine = Arc::new(RwLock::new(engine));
    let mut area = MainAreaDna::new(dna, Some("tpl".to_string()), Some(engine));
    area.primer_design_ui.report_id = "primer_ui_seed".to_string();
    area.primer_design_ui
        .restriction_cloning
        .selected_pair_rank_1based = "1".to_string();
    area.primer_design_ui
        .restriction_cloning
        .destination_vector_seq_id = "vec".to_string();
    area.primer_design_ui.restriction_cloning.mode = RestrictionCloningPcrHandoffMode::SingleSite;
    area.primer_design_ui
        .restriction_cloning
        .forward_enzyme
        .clear();
    area.primer_design_ui
        .restriction_cloning
        .reverse_enzyme
        .clear();

    let seeded = area
        .seed_restriction_cloning_pcr_handoff_request()
        .expect("seed restriction-cloning request");

    assert_eq!(seeded.selection_source, "recommended_single_site");
    assert_eq!(seeded.forward_enzyme, "EcoRI");
    assert_eq!(seeded.reverse_enzyme, "EcoRI");
    match seeded.operation {
        Operation::PrepareRestrictionCloningPcrHandoff {
            forward_enzyme,
            reverse_enzyme,
            ..
        } => {
            assert_eq!(forward_enzyme, "EcoRI");
            assert_eq!(reverse_enzyme.as_deref(), Some("EcoRI"));
        }
        other => panic!("unexpected seeded operation: {other:?}"),
    }
}

#[test]
fn seed_restriction_cloning_handoff_ui_surfaces_engine_order_validation() {
    let mut state = ProjectState::default();
    state.sequences.insert(
        "tpl".to_string(),
        DNAsequence::from_sequence(
            "ACGTTGCATGTCAGTACGATCGTACGTAGCTAGTCGATCGTACGATCGTAGCTAGCATCGATGCTAGCTAGTACGTAGCATCGATCGTAGCTAGCATGCTAGCTAGTCGATCGATCGTACGATCG",
        )
        .expect("sequence"),
    );
    let mut vector =
        DNAsequence::from_sequence("AAAAGAATTCGGGGGAAGCTTTTTTGCGGCCGCTTTT").expect("vector");
    *vector.restriction_enzymes_mut() = active_restriction_enzymes();
    let vector_len_i64 = vector.len().try_into().unwrap();
    vector.features_mut().push(gb_io::seq::Feature {
        kind: "misc_feature".into(),
        location: gb_io::seq::Location::simple_range(0, vector_len_i64),
        qualifiers: vec![
            ("label".into(), Some("MCS".to_string())),
            (
                "mcs_expected_sites".into(),
                Some("EcoRI,HindIII".to_string()),
            ),
        ],
    });
    vector.update_computed_features();
    state.sequences.insert("vec".to_string(), vector);
    let mut engine = GentleEngine::from_state(state);
    engine.state_mut().parameters.primer_design_backend = PrimerDesignBackend::Internal;
    engine
        .apply(Operation::DesignPrimerPairs {
            template: "tpl".to_string(),
            roi_start_0based: 40,
            roi_end_0based: 80,
            forward: PrimerDesignSideConstraint {
                min_length: 20,
                max_length: 20,
                location_0based: Some(5),
                start_0based: None,
                end_0based: None,
                min_tm_c: 0.0,
                max_tm_c: 100.0,
                min_gc_fraction: 0.0,
                max_gc_fraction: 1.0,
                max_anneal_hits: 1000,
                ..Default::default()
            },
            reverse: PrimerDesignSideConstraint {
                min_length: 20,
                max_length: 20,
                location_0based: Some(90),
                start_0based: None,
                end_0based: None,
                min_tm_c: 0.0,
                max_tm_c: 100.0,
                min_gc_fraction: 0.0,
                max_gc_fraction: 1.0,
                max_anneal_hits: 1000,
                ..Default::default()
            },
            pair_constraints: PrimerDesignPairConstraint::default(),
            min_amplicon_bp: 40,
            max_amplicon_bp: 150,
            max_tm_delta_c: Some(100.0),
            max_pairs: Some(10),
            report_id: Some("primer_ui_seed_bad_order".to_string()),
        })
        .expect("design primer pairs");
    let dna = engine
        .state()
        .sequences
        .get("tpl")
        .cloned()
        .expect("template sequence");
    let engine = Arc::new(RwLock::new(engine));
    let mut area = MainAreaDna::new(dna, Some("tpl".to_string()), Some(engine));
    area.primer_design_ui.report_id = "primer_ui_seed_bad_order".to_string();
    area.primer_design_ui
        .restriction_cloning
        .selected_pair_rank_1based = "1".to_string();
    area.primer_design_ui
        .restriction_cloning
        .destination_vector_seq_id = "vec".to_string();
    area.primer_design_ui.restriction_cloning.mode = RestrictionCloningPcrHandoffMode::DirectedPair;
    area.primer_design_ui.restriction_cloning.forward_enzyme = "HindIII".to_string();
    area.primer_design_ui.restriction_cloning.reverse_enzyme = "EcoRI".to_string();

    let err = area
        .seed_restriction_cloning_pcr_handoff_request()
        .expect_err("reverse order should be rejected by engine seed helper");

    assert!(err.contains("does not match the valid order"));
    assert!(err.contains("HindIII -> EcoRI"));
}

#[test]
fn feature_tree_matches_filter_scoped_terms() {
    let feature = make_feature("gene", vec![("label", "TP53")]);
    assert!(MainAreaDna::feature_tree_matches_filter(
        &feature,
        "kind:gene label:tp53",
        "gene",
        "TP53",
        "1..10"
    ));
    assert!(!MainAreaDna::feature_tree_matches_filter(
        &feature,
        "kind:mrna",
        "gene",
        "TP53",
        "1..10"
    ));
}

#[test]
fn feature_tree_matches_filter_scoped_track_source_and_path_terms() {
    let feature = make_feature(
        "track",
        vec![
            ("gentle_track_source", "BED"),
            ("gentle_track_name", "chip_signal_tp73"),
            ("gentle_track_file", "/tmp/tracks/tp73_peaks.bed"),
            ("label", "TP73 peaks"),
            ("note", "enhancer-like signal"),
        ],
    );
    assert!(MainAreaDna::feature_tree_matches_filter(
        &feature,
        "source:bed track:chip path:tp73_peaks.bed",
        "tracks",
        "TP73 peaks",
        "1..10"
    ));
    assert!(MainAreaDna::feature_tree_matches_filter(
        &feature,
        "file:tp73_peaks.bed note:enhancer",
        "tracks",
        "TP73 peaks",
        "1..10"
    ));
    assert!(!MainAreaDna::feature_tree_matches_filter(
        &feature,
        "path:missing.bed",
        "tracks",
        "TP73 peaks",
        "1..10"
    ));
}

#[test]
fn feature_tree_matches_filter_scoped_range_terms() {
    let feature = make_feature("gene", vec![("label", "TP73")]);
    assert!(MainAreaDna::feature_tree_matches_filter(
        &feature,
        "range:6128..16430",
        "gene",
        "TP73",
        "6128..16430"
    ));
    assert!(!MainAreaDna::feature_tree_matches_filter(
        &feature,
        "range:7000..8000",
        "gene",
        "TP73",
        "6128..16430"
    ));
}

#[test]
fn feature_tree_should_group_respects_grouping_mode() {
    assert!(!MainAreaDna::feature_tree_should_group(
        super::FeatureTreeGroupingMode::Off,
        8
    ));
    assert!(!MainAreaDna::feature_tree_should_group(
        super::FeatureTreeGroupingMode::Auto,
        1
    ));
    assert!(MainAreaDna::feature_tree_should_group(
        super::FeatureTreeGroupingMode::Auto,
        2
    ));
    assert!(MainAreaDna::feature_tree_should_group(
        super::FeatureTreeGroupingMode::Always,
        1
    ));
}

#[test]
fn feature_tree_display_label_disambiguates_rna_by_transcript_id() {
    let feature = make_feature(
        "mRNA",
        vec![("label", "TP73-AS3"), ("transcript_id", "NM_001126112.3")],
    );
    assert_eq!(
        MainAreaDna::feature_tree_display_label(&feature, "TP73-AS3", "6128..16430"),
        "TP73-AS3 [NM_001126112.3]"
    );
}

#[test]
fn feature_tree_display_label_disambiguates_rna_by_range_fallback() {
    let feature = make_feature("ncRNA", vec![("label", "TP73-AS3")]);
    assert_eq!(
        MainAreaDna::feature_tree_display_label(&feature, "TP73-AS3", "6128..16430"),
        "TP73-AS3 [6128..16430]"
    );
}

#[test]
fn feature_tree_subgroup_label_does_not_group_gene_features() {
    let feature = make_feature("gene", vec![("gene", "TP73")]);
    assert_eq!(
        MainAreaDna::feature_tree_subgroup_label(
            &feature,
            "TP73",
            super::FeatureTreeGroupingMode::Always
        ),
        None
    );
}

#[test]
fn feature_tree_model_keeps_gene_rows_splicing_action_capable() {
    let mut dna = DNAsequence::from_sequence("AAAAAAAAAAAAAAAAAAAA").expect("sequence");
    dna.features_mut().push(make_feature(
        "gene",
        vec![("gene", "TP73"), ("label", "TP73")],
    ));
    let mut area = MainAreaDna::new(dna, Some("seq_gene".to_string()), None);
    area.feature_tree_grouping_mode = super::FeatureTreeGroupingMode::Always;

    let key = area.current_feature_tree_cache_key(None);
    let model = area.build_feature_tree_model(&key);
    let gene_group = model
        .groups
        .iter()
        .find(|group| group.kind.eq_ignore_ascii_case("gene"))
        .expect("gene group");

    assert_eq!(gene_group.ungrouped_entry_indices, vec![0]);
    let entry = &gene_group.entries[gene_group.ungrouped_entry_indices[0]];
    assert!(entry.disable_grouping);
    assert!(entry.supports_splicing_expert);
}

#[test]
fn feature_tree_model_groups_dense_regulatory_repeat_track_and_array_surfaces() {
    let mut dna = DNAsequence::from_sequence("A".repeat(240).as_str()).expect("sequence");
    dna.features_mut().push(Feature {
        kind: "regulatory".into(),
        location: Location::simple_range(10, 30),
        qualifiers: vec![
            ("regulatory_class".into(), Some("enhancer".to_string())),
            (
                "note".into(),
                Some("H3K4me1 active enhancer one".to_string()),
            ),
        ],
    });
    dna.features_mut().push(Feature {
        kind: "regulatory".into(),
        location: Location::simple_range(40, 60),
        qualifiers: vec![
            ("regulatory_class".into(), Some("enhancer".to_string())),
            (
                "note".into(),
                Some("H3K4me1 active enhancer two".to_string()),
            ),
        ],
    });
    dna.features_mut().push(Feature {
        kind: "repeat_region".into(),
        location: Location::simple_range(70, 90),
        qualifiers: vec![
            ("repName".into(), Some("L1PA2".to_string())),
            ("repClass".into(), Some("LINE".to_string())),
            ("repFamily".into(), Some("L1".to_string())),
        ],
    });
    dna.features_mut().push(Feature {
        kind: "track".into(),
        location: Location::simple_range(100, 120),
        qualifiers: vec![
            (
                "gentle_generated".into(),
                Some("genome_bed_track".to_string()),
            ),
            ("gentle_track_source".into(), Some("BED".to_string())),
            ("gentle_track_name".into(), Some("H3K27ac".to_string())),
            (
                "gentle_track_file".into(),
                Some("/tmp/peaks.bed".to_string()),
            ),
            ("label".into(), Some("peak one".to_string())),
        ],
    });
    dna.features_mut().push(Feature {
        kind: "track".into(),
        location: Location::simple_range(130, 150),
        qualifiers: vec![
            (
                "gentle_generated".into(),
                Some("microarray_track_projection".to_string()),
            ),
            ("gentle_track_source".into(), Some("Array".to_string())),
            (
                "gentle_array_platform".into(),
                Some("Clariom D".to_string()),
            ),
            (
                "gentle_array_contrast".into(),
                Some("AdTAp73alpha-AdGFP".to_string()),
            ),
            ("feature_id".into(), Some("PSR0001".to_string())),
        ],
    });
    let mut area = MainAreaDna::new(dna, Some("dense_tree".to_string()), None);
    area.feature_tree_grouping_mode = super::FeatureTreeGroupingMode::Always;

    let key = area.current_feature_tree_cache_key(None);
    let model = area.build_feature_tree_model(&key);

    let regulatory_group = model
        .groups
        .iter()
        .find(|group| group.kind.eq_ignore_ascii_case("regulatory"))
        .expect("regulatory group");
    let enhancer_group = regulatory_group
        .regulatory_primary_groups
        .iter()
        .find(|group| group.label == "enhancer")
        .expect("enhancer primary group");
    assert_eq!(enhancer_group.entry_indices.len(), 2);
    assert_eq!(enhancer_group.secondary_groups[0].label, "H3K4me1");

    let repeat_group = model
        .groups
        .iter()
        .find(|group| group.kind.eq_ignore_ascii_case("repeat_region"))
        .expect("repeat group");
    assert_eq!(repeat_group.grouped_entries[0].label, "LINE / L1");

    let track_group = model
        .groups
        .iter()
        .find(|group| group.kind.eq_ignore_ascii_case("tracks"))
        .expect("generic track group");
    assert_eq!(
        track_group.grouped_entries[0].label,
        "BED: peaks.bed (H3K27ac)"
    );

    let array_group = model
        .groups
        .iter()
        .find(|group| group.kind.eq_ignore_ascii_case("array"))
        .expect("array group");
    assert_eq!(
        array_group.grouped_entries[0].label,
        "Array: Clariom D (AdTAp73alpha-AdGFP)"
    );
}

#[test]
fn feature_tree_layer_summary_labels_count_current_model_categories() {
    let mut dna = DNAsequence::from_sequence("A".repeat(240).as_str()).expect("sequence");
    dna.features_mut()
        .push(make_feature("gene", vec![("gene", "TP73")]));
    dna.features_mut().push(make_feature(
        "mRNA",
        vec![("gene", "TP73"), ("transcript_id", "TX1")],
    ));
    dna.features_mut().push(make_feature(
        "regulatory",
        vec![("regulatory_class", "enhancer"), ("note", "H3K27ac")],
    ));
    dna.features_mut().push(make_feature(
        "repeat_region",
        vec![
            ("repName", "AluY"),
            ("repClass", "SINE"),
            ("repFamily", "Alu"),
        ],
    ));
    dna.features_mut().push(make_feature(
        "track",
        vec![
            ("gentle_generated", "genome_bed_track"),
            ("gentle_track_source", "BED"),
            ("gentle_track_name", "H3K27ac"),
            ("gentle_track_file", "/tmp/peaks.bed"),
        ],
    ));
    dna.features_mut().push(make_feature(
        "track",
        vec![
            ("gentle_generated", "microarray_track_projection"),
            ("gentle_track_source", "Array"),
            ("gentle_array_platform", "Clariom D"),
            ("gentle_array_contrast", "AdTAp73alpha-AdGFP"),
        ],
    ));
    let area = MainAreaDna::new(dna, Some("summary_tree".to_string()), None);
    let key = area.current_feature_tree_cache_key(None);
    let model = area.build_feature_tree_model(&key);

    let labels = MainAreaDna::feature_tree_layer_summary_labels(&model, false);

    assert!(labels.contains(&"core (2)".to_string()));
    assert!(labels.contains(&"regulatory (1)".to_string()));
    assert!(labels.contains(&"repeats (1)".to_string()));
    assert!(labels.contains(&"tracks (1)".to_string()));
    assert!(labels.contains(&"array (1)".to_string()));
}

#[test]
fn feature_tree_ui_opens_group_and_subgroup_before_rendering_rows() {
    let mut dna = DNAsequence::from_sequence("A".repeat(200).as_str()).expect("sequence");
    dna.features_mut().push(Feature {
        kind: "mRNA".into(),
        location: Location::simple_range(10, 70),
        qualifiers: vec![
            ("gene".into(), Some("TP73".to_string())),
            ("label".into(), Some("TP73-201".to_string())),
            ("transcript_id".into(), Some("TX1".to_string())),
        ],
    });
    dna.features_mut().push(Feature {
        kind: "mRNA".into(),
        location: Location::simple_range(90, 150),
        qualifiers: vec![
            ("gene".into(), Some("TP73".to_string())),
            ("label".into(), Some("TP73-202".to_string())),
            ("transcript_id".into(), Some("TX2".to_string())),
        ],
    });
    let mut area = MainAreaDna::new(dna, Some("feature_tree_ui".to_string()), None);
    area.feature_tree_grouping_mode = super::FeatureTreeGroupingMode::Always;

    let ctx = egui::Context::default();
    let screen_rect = egui::Rect::from_min_size(egui::Pos2::ZERO, egui::vec2(720.0, 520.0));
    let rects = render_feature_tree_pass(
        &ctx,
        &mut area,
        egui::RawInput {
            screen_rect: Some(screen_rect),
            ..Default::default()
        },
    );
    let mrna_center =
        center_of_rendered_text(&rects, "mRNA (2)").expect("mRNA group header should render");
    assert!(
        center_of_rendered_text(&rects, "TP73 (2)").is_none(),
        "subgroup should stay hidden while the kind group is collapsed: {rects:?}"
    );

    for pressed in [true, false] {
        render_feature_tree_pass(
            &ctx,
            &mut area,
            egui::RawInput {
                screen_rect: Some(screen_rect),
                events: vec![
                    egui::Event::PointerMoved(mrna_center),
                    egui::Event::PointerButton {
                        pos: mrna_center,
                        button: egui::PointerButton::Primary,
                        pressed,
                        modifiers: egui::Modifiers::default(),
                    },
                ],
                ..Default::default()
            },
        );
    }

    // Advance test time so the first collapsing animation has a drawable body
    // before the nested subgroup header is targeted.
    let rects = render_feature_tree_pass(
        &ctx,
        &mut area,
        egui::RawInput {
            screen_rect: Some(screen_rect),
            time: Some(1.0),
            ..Default::default()
        },
    );
    let subgroup_center = center_of_rendered_text(&rects, "TP73 (2)")
        .expect("TP73 subgroup header should render after opening mRNA group");
    let subgroup_toggle_center = subgroup_center - egui::vec2(20.0, 0.0);
    assert!(
        center_of_rendered_text(&rects, "TX1").is_none(),
        "transcript rows should stay hidden while the subgroup is collapsed: {rects:?}"
    );

    render_feature_tree_pass(
        &ctx,
        &mut area,
        egui::RawInput {
            screen_rect: Some(screen_rect),
            events: vec![egui::Event::PointerMoved(subgroup_toggle_center)],
            ..Default::default()
        },
    );
    render_feature_tree_pass(
        &ctx,
        &mut area,
        egui::RawInput {
            screen_rect: Some(screen_rect),
            events: vec![
                egui::Event::PointerMoved(subgroup_toggle_center),
                egui::Event::PointerButton {
                    pos: subgroup_toggle_center,
                    button: egui::PointerButton::Primary,
                    pressed: true,
                    modifiers: egui::Modifiers::default(),
                },
                egui::Event::PointerButton {
                    pos: subgroup_toggle_center,
                    button: egui::PointerButton::Primary,
                    pressed: false,
                    modifiers: egui::Modifiers::default(),
                },
            ],
            ..Default::default()
        },
    );

    let rects = render_feature_tree_pass(
        &ctx,
        &mut area,
        egui::RawInput {
            screen_rect: Some(screen_rect),
            ..Default::default()
        },
    );
    assert!(
        center_of_rendered_text(&rects, "TX1").is_some(),
        "first transcript row should render after opening the subgroup: {rects:?}"
    );
}

#[test]
fn feature_tree_ui_auto_opens_selected_regulatory_nested_group() {
    let mut dna = DNAsequence::from_sequence("A".repeat(200).as_str()).expect("sequence");
    dna.features_mut().push(Feature {
        kind: "regulatory".into(),
        location: Location::simple_range(10, 30),
        qualifiers: vec![
            ("regulatory_class".into(), Some("enhancer".to_string())),
            (
                "note".into(),
                Some("H3K4me1 active enhancer one".to_string()),
            ),
        ],
    });
    dna.features_mut().push(Feature {
        kind: "regulatory".into(),
        location: Location::simple_range(40, 60),
        qualifiers: vec![
            ("regulatory_class".into(), Some("enhancer".to_string())),
            (
                "note".into(),
                Some("H3K4me1 active enhancer two".to_string()),
            ),
        ],
    });
    let mut area = MainAreaDna::new(dna, Some("regulatory_tree_ui".to_string()), None);
    area.feature_tree_grouping_mode = super::FeatureTreeGroupingMode::Always;
    area.focus_feature(0);

    let ctx = egui::Context::default();
    let rects = render_feature_tree_pass(
        &ctx,
        &mut area,
        egui::RawInput {
            screen_rect: Some(egui::Rect::from_min_size(
                egui::Pos2::ZERO,
                egui::vec2(720.0, 520.0),
            )),
            ..Default::default()
        },
    );

    assert!(center_of_rendered_text(&rects, "enhancer (2)").is_some());
    assert!(center_of_rendered_text(&rects, "H3K4me1 (2)").is_some());
    assert!(
        center_of_rendered_text(&rects, "active enhancer one").is_some(),
        "selected nested regulatory row should render after parent groups auto-open: {rects:?}"
    );
}

#[test]
fn feature_tree_ui_auto_opens_selected_repeat_group() {
    let mut dna = DNAsequence::from_sequence("A".repeat(200).as_str()).expect("sequence");
    dna.features_mut().push(Feature {
        kind: "repeat_region".into(),
        location: Location::simple_range(10, 20),
        qualifiers: vec![
            ("repName".into(), Some("L1PA2".to_string())),
            ("repClass".into(), Some("LINE".to_string())),
            ("repFamily".into(), Some("L1".to_string())),
        ],
    });
    dna.features_mut().push(Feature {
        kind: "repeat_region".into(),
        location: Location::simple_range(40, 55),
        qualifiers: vec![
            ("repName".into(), Some("L1PA3".to_string())),
            ("repClass".into(), Some("LINE".to_string())),
            ("repFamily".into(), Some("L1".to_string())),
        ],
    });
    let mut area = MainAreaDna::new(dna, Some("repeat_tree_ui".to_string()), None);
    area.feature_tree_grouping_mode = super::FeatureTreeGroupingMode::Always;
    area.focus_feature(0);

    let ctx = egui::Context::default();
    let rects = render_feature_tree_pass(
        &ctx,
        &mut area,
        egui::RawInput {
            screen_rect: Some(egui::Rect::from_min_size(
                egui::Pos2::ZERO,
                egui::vec2(720.0, 520.0),
            )),
            ..Default::default()
        },
    );

    assert!(center_of_rendered_text(&rects, "repeat_region (2)").is_some());
    assert!(center_of_rendered_text(&rects, "LINE / L1 (2)").is_some());
    assert!(
        center_of_rendered_text(&rects, "11..20 (10 bp)").is_some(),
        "selected repeat row should render after parent group auto-opens: {rects:?}"
    );
}

#[test]
fn derive_regulatory_feature_grouping_groups_enhancers_by_marker_and_strips_prefix() {
    let feature = make_feature(
        "regulatory",
        vec![
            ("regulatory_class", "enhancer"),
            (
                "note",
                "H3K4me1 hESC enhancer chr1:3600129-3600889 (GRCh37/hg19 assembly coordinates)",
            ),
        ],
    );
    let grouped = MainAreaDna::derive_regulatory_feature_grouping(
        &feature,
        "enhancer: H3K4me1 hESC enhancer chr1:3600129-3600889 (GRCh37/hg19 assembly coordinates)",
    )
    .expect("expected regulatory grouping");
    assert_eq!(grouped.primary_key, "enhancer");
    assert_eq!(grouped.secondary_label.as_deref(), Some("H3K4me1"));
    assert!(
        !grouped
            .display_label
            .to_ascii_lowercase()
            .starts_with("enhancer")
    );
    assert!(grouped.display_label.contains("chr1:"));
}

#[test]
fn derive_regulatory_feature_grouping_uses_active_region_secondary_group() {
    let feature = make_feature(
        "regulatory",
        vec![
            ("regulatory_class", "silencer"),
            ("note", "active region_65"),
        ],
    );
    let grouped =
        MainAreaDna::derive_regulatory_feature_grouping(&feature, "silencer: active region_65")
            .expect("expected regulatory grouping");
    assert_eq!(grouped.primary_key, "silencer");
    assert_eq!(grouped.secondary_key.as_deref(), Some("active_region"));
    assert_eq!(grouped.secondary_label.as_deref(), Some("active region"));
    assert!(
        !grouped
            .display_label
            .to_ascii_lowercase()
            .starts_with("active region")
    );
}

#[test]
fn derive_regulatory_feature_grouping_accepts_active_region_colon_prefix_variant() {
    let feature = make_feature(
        "regulatory",
        vec![
            ("regulatory_class", "silencer"),
            ("note", "active region: 122"),
        ],
    );
    let grouped =
        MainAreaDna::derive_regulatory_feature_grouping(&feature, "silencer: active region: 122")
            .expect("expected active-region grouping");
    assert_eq!(grouped.primary_key, "silencer");
    assert_eq!(grouped.secondary_key.as_deref(), Some("active_region"));
    assert_eq!(grouped.secondary_label.as_deref(), Some("active region"));
    assert_eq!(grouped.display_label, "122");
}

#[test]
fn derive_regulatory_feature_grouping_groups_enhancers_by_non_histone_marker() {
    let feature = make_feature(
        "regulatory",
        vec![
            ("regulatory_class", "enhancer"),
            (
                "note",
                "CTCF-bound enhancer chr1:3650019-3650807 (GRCh37/hg19 assembly coordinates)",
            ),
        ],
    );
    let grouped = MainAreaDna::derive_regulatory_feature_grouping(
        &feature,
        "enhancer: CTCF-bound enhancer chr1:3650019-3650807 (GRCh37/hg19 assembly coordinates)",
    )
    .expect("expected marker grouping");
    assert_eq!(grouped.primary_key, "enhancer");
    assert_eq!(grouped.secondary_label.as_deref(), Some("CTCF"));
    assert_eq!(grouped.secondary_key.as_deref(), Some("marker:ctcf"));
    assert!(
        grouped
            .display_label
            .to_ascii_lowercase()
            .starts_with("bound enhancer chr1:")
    );
}

#[test]
fn derive_regulatory_feature_grouping_falls_back_to_label_prefix_when_class_missing() {
    let feature = make_feature(
        "regulatory",
        vec![(
            "note",
            "fragment chr1:3596796-3596963 (GRCh37/hg19 assembly coordinates)",
        )],
    );
    let grouped = MainAreaDna::derive_regulatory_feature_grouping(
        &feature,
        "silencer: fragment chr1:3596796-3596963 (GRCh37/hg19 assembly coordinates)",
    )
    .expect("expected fallback grouping");
    assert_eq!(grouped.primary_key, "silencer");
    assert!(grouped.secondary_key.is_none());
    assert!(
        !grouped
            .display_label
            .to_ascii_lowercase()
            .starts_with("silencer")
    );
}

#[test]
fn derive_regulatory_feature_grouping_routes_non_enhancer_silencer_classes_to_other() {
    let feature = make_feature(
        "regulatory",
        vec![
            ("regulatory_class", "promoter"),
            ("note", "promoter: TATA-rich region chr1:3650019-3650101"),
        ],
    );
    let grouped = MainAreaDna::derive_regulatory_feature_grouping(
        &feature,
        "promoter: TATA-rich region chr1:3650019-3650101",
    )
    .expect("expected other-group routing");
    assert_eq!(grouped.primary_key, "other");
    assert_eq!(grouped.primary_label, "other");
    assert!(grouped.secondary_key.is_none());
    assert!(
        grouped
            .display_label
            .to_ascii_lowercase()
            .starts_with("promoter")
    );
}

#[test]
fn derive_regulatory_feature_grouping_routes_unclassified_labels_to_other() {
    let feature = make_feature("regulatory", vec![("note", "candidate boundary region")]);
    let grouped =
        MainAreaDna::derive_regulatory_feature_grouping(&feature, "candidate boundary region")
            .expect("expected other-group routing");
    assert_eq!(grouped.primary_key, "other");
    assert_eq!(grouped.primary_label, "other");
    assert!(grouped.secondary_key.is_none());
}

#[test]
fn feature_tree_subgroup_label_groups_mrna_by_gene() {
    let feature = make_feature(
        "mRNA",
        vec![
            ("label", "NM_001126112.3"),
            ("gene", "TP73"),
            ("transcript_id", "NM_001126112.3"),
        ],
    );
    assert_eq!(
        MainAreaDna::feature_tree_subgroup_label(
            &feature,
            "NM_001126112.3",
            super::FeatureTreeGroupingMode::Auto
        ),
        Some("TP73".to_string())
    );
}

#[test]
fn feature_tree_subgroup_label_groups_rmsk_repeats_by_class_and_family() {
    let feature = make_feature(
        "repeat_region",
        vec![
            ("repName", "L1PA2"),
            ("repClass", "LINE"),
            ("repFamily", "L1"),
        ],
    );

    assert_eq!(
        MainAreaDna::feature_tree_subgroup_label(
            &feature,
            "L1PA2 (LINE / L1)",
            super::FeatureTreeGroupingMode::Always,
        ),
        Some("LINE / L1".to_string())
    );
    assert!(MainAreaDna::feature_tree_matches_filter(
        &feature,
        "repeat_class:line repfamily:l1",
        "repeat_region",
        "L1PA2 (LINE / L1)",
        "1..10"
    ));
}

#[test]
fn seed_anchored_promoter_from_forward_mrna_feature_uses_start_and_upstream() {
    let mut dna = DNAsequence::from_sequence("A".repeat(200).as_str()).expect("sequence");
    dna.features_mut().push(Feature {
        kind: "mRNA".into(),
        location: Location::simple_range(10, 40),
        qualifiers: vec![("gene".into(), Some("TP73".to_string()))],
    });
    let mut area = MainAreaDna::new(dna, None, None);

    area.seed_anchored_promoter_from_feature_id(0);

    assert!(area.anchored_mode_feature);
    assert_eq!(area.anchored_feature_kind, "mRNA");
    assert_eq!(area.anchored_feature_label, "TP73");
    assert!(area.anchored_feature_boundary_start);
    assert!(area.anchored_direction_upstream);
    assert_eq!(area.anchored_output_prefix, "promoter_tp73");
}

#[test]
fn seed_anchored_promoter_from_reverse_mrna_feature_uses_end_and_downstream() {
    let mut dna = DNAsequence::from_sequence("A".repeat(200).as_str()).expect("sequence");
    dna.features_mut().push(Feature {
        kind: "mRNA".into(),
        location: Location::Complement(Box::new(Location::simple_range(10, 40))),
        qualifiers: vec![("gene".into(), Some("TP73".to_string()))],
    });
    let mut area = MainAreaDna::new(dna, None, None);

    area.seed_anchored_promoter_from_feature_id(0);

    assert_eq!(area.anchored_feature_kind, "mRNA");
    assert_eq!(area.anchored_feature_label, "TP73");
    assert!(!area.anchored_feature_boundary_start);
    assert!(!area.anchored_direction_upstream);
    assert_eq!(area.anchored_output_prefix, "promoter_tp73");
}

#[test]
fn refresh_from_engine_sequence_state_reloads_open_sequence() {
    let initial_dna = DNAsequence::from_sequence("AAAA").expect("sequence");
    let mut updated_dna = DNAsequence::from_sequence("AAAA").expect("sequence");
    updated_dna
        .features_mut()
        .push(make_feature("gene", vec![("label", "TP53")]));

    let mut state = ProjectState::default();
    state.sequences.insert("seq1".to_string(), updated_dna);
    let engine = Arc::new(RwLock::new(GentleEngine::from_state(state)));
    let mut area = MainAreaDna::new(initial_dna, Some("seq1".to_string()), Some(engine));

    area.refresh_from_engine_sequence_state();

    let dna = area.dna().read().expect("DNA lock poisoned");
    assert_eq!(dna.features().len(), 1);
}

#[test]
fn linear_sequence_window_hides_text_panel_by_default() {
    let dna = DNAsequence::from_sequence("AAAA").expect("sequence");
    let area = MainAreaDna::new(dna, None, None);

    assert!(!area.show_sequence);
}

#[test]
fn circular_sequence_window_keeps_text_panel_by_default() {
    let mut dna = DNAsequence::from_sequence("AAAA").expect("sequence");
    dna.set_circular(true);
    let area = MainAreaDna::new(dna, None, None);

    assert!(area.show_sequence);
}

#[test]
fn engine_backed_linear_sequence_window_uses_linear_text_panel_default() {
    let dna = DNAsequence::from_sequence("AAAA").expect("sequence");
    let mut state = ProjectState::default();
    state.sequences.insert("seq1".to_string(), dna.clone());
    let engine = Arc::new(RwLock::new(GentleEngine::from_state(state)));
    let area = MainAreaDna::new(dna, Some("seq1".to_string()), Some(engine));

    assert!(!area.show_sequence);
}

#[test]
fn engine_backed_circular_sequence_window_uses_circular_text_panel_default() {
    let mut dna = DNAsequence::from_sequence("AAAA").expect("sequence");
    dna.set_circular(true);
    let mut state = ProjectState::default();
    state.sequences.insert("seq1".to_string(), dna.clone());
    let engine = Arc::new(RwLock::new(GentleEngine::from_state(state)));
    let area = MainAreaDna::new(dna, Some("seq1".to_string()), Some(engine));

    assert!(area.show_sequence);
}

#[test]
fn compact_lane_layout_hides_sequence_panel_without_changing_shared_display_default() {
    let dna = DNAsequence::from_sequence("AAAA").expect("sequence");
    let mut state = ProjectState::default();
    state.sequences.insert("seq1".to_string(), dna.clone());
    state.display.show_sequence_panel = true;
    state.display.show_map_panel = true;
    let engine = Arc::new(RwLock::new(GentleEngine::from_state(state)));
    let mut area = MainAreaDna::new(dna, Some("seq1".to_string()), Some(engine.clone()));

    area.enable_compact_lane_layout();
    area.refresh_from_engine_settings();

    assert!(!area.show_sequence);
    assert!(area.show_map);
    assert!(
        engine
            .read()
            .expect("engine lock")
            .state()
            .display
            .show_sequence_panel,
        "compact lane windows should not flip the shared display default"
    );
}

#[test]
fn replace_active_dna_marks_layout_dirty_and_expands_full_view_on_growth() {
    let initial_dna = DNAsequence::from_sequence("A".repeat(100).as_str()).expect("sequence");
    let mut area = MainAreaDna::new(initial_dna, None, None);
    area.set_linear_viewport(0, 100);

    let grown_dna = DNAsequence::from_sequence("A".repeat(180).as_str()).expect("sequence");
    area.replace_active_dna(grown_dna, true);

    let (start, span, seq_len) = area.current_linear_viewport();
    assert_eq!(start, 0);
    assert_eq!(span, 180);
    assert_eq!(seq_len, 180);
    let display = area.dna_display.read().expect("display lock");
    assert!(display.update_layout().update_map_dna());
    assert!(display.update_layout().update_map_sequence());
}

#[test]
fn replace_active_dna_preserves_non_full_linear_viewport_on_growth() {
    let initial_dna = DNAsequence::from_sequence("A".repeat(200).as_str()).expect("sequence");
    let mut area = MainAreaDna::new(initial_dna, None, None);
    area.set_linear_viewport(40, 60);

    let grown_dna = DNAsequence::from_sequence("A".repeat(320).as_str()).expect("sequence");
    area.replace_active_dna(grown_dna, true);

    let (start, span, seq_len) = area.current_linear_viewport();
    assert_eq!(start, 40);
    assert_eq!(span, 60);
    assert_eq!(seq_len, 320);
}

#[test]
fn replace_active_dna_limits_initial_linear_view_after_circular_topology_switch() {
    let mut initial_dna =
        DNAsequence::from_sequence("A".repeat(120_000).as_str()).expect("sequence");
    initial_dna.set_circular(true);
    let mut area = MainAreaDna::new(initial_dna, None, None);

    let mut linear_dna =
        DNAsequence::from_sequence("A".repeat(120_000).as_str()).expect("sequence");
    linear_dna.set_circular(false);
    area.replace_active_dna(linear_dna, true);

    let (start, span, seq_len) = area.current_linear_viewport();
    assert_eq!(start, 0);
    assert_eq!(span, LINEAR_TOPOLOGY_SWITCH_MAX_INITIAL_SPAN_BP);
    assert_eq!(seq_len, 120_000);
    assert!(!area.map_dna.is_circular());
}

#[test]
fn replace_active_dna_invalidates_feature_tree_cache() {
    let mut initial_dna = DNAsequence::from_sequence("A".repeat(100).as_str()).expect("sequence");
    initial_dna
        .features_mut()
        .push(make_feature("gene", vec![("label", "TP53")]));
    let mut area = MainAreaDna::new(initial_dna, None, None);

    area.ensure_feature_tree_cache_current(None);
    assert!(area.feature_tree_cache.is_some());

    let mut replacement_dna =
        DNAsequence::from_sequence("A".repeat(120).as_str()).expect("sequence");
    replacement_dna
        .features_mut()
        .push(make_feature("gene", vec![("label", "TP73")]));
    area.replace_active_dna(replacement_dna, true);

    assert!(area.feature_tree_cache.is_none());
}

#[test]
fn feature_tree_cache_rebuilds_when_filter_changes() {
    let mut dna = DNAsequence::from_sequence("A".repeat(100).as_str()).expect("sequence");
    dna.features_mut()
        .push(make_feature("gene", vec![("label", "TP53")]));
    dna.features_mut().push(make_feature(
        "mRNA",
        vec![("gene", "TP53"), ("label", "TP53-201")],
    ));
    let mut area = MainAreaDna::new(dna, None, None);

    area.ensure_feature_tree_cache_current(None);
    let initial_entry_count: usize = area
        .feature_tree_cache
        .as_ref()
        .expect("feature tree cache")
        .model
        .groups
        .iter()
        .map(|group| group.entries.len())
        .sum();
    assert_eq!(initial_entry_count, 2);

    area.feature_tree_filter = "label:missing".to_string();
    area.ensure_feature_tree_cache_current(None);

    let cache = area
        .feature_tree_cache
        .as_ref()
        .expect("feature tree cache after filter change");
    assert_eq!(cache.key.feature_filter_text, "label:missing");
    assert_eq!(cache.model.filter_total_count, 2);
    assert_eq!(cache.model.filter_matched_count, 0);
    assert!(cache.model.groups.is_empty());
}

#[test]
fn ensure_full_restriction_site_catalog_current_invalidates_feature_tree_cache() {
    let dna = restriction_ready_dna("GAATTCGAATTC");
    let mut area = MainAreaDna::new(dna, None, None);

    area.ensure_feature_tree_cache_current(None);
    assert!(area.feature_tree_cache.is_some());

    area.ensure_full_restriction_site_catalog_current();

    assert!(area.feature_tree_cache.is_none());
}

#[test]
fn feature_multi_select_toggle_adds_and_removes_ids() {
    let mut dna = DNAsequence::from_sequence("A".repeat(200).as_str()).expect("sequence");
    dna.features_mut()
        .push(make_feature("gene", vec![("label", "gene_1")]));
    dna.features_mut()
        .push(make_feature("gene", vec![("label", "gene_2")]));
    let mut area = MainAreaDna::new(dna, None, None);

    area.focus_feature(0);
    area.toggle_feature_multi_select(1);
    assert_eq!(
        area.multi_selected_feature_ids,
        BTreeSet::from([0usize, 1usize])
    );
    assert_eq!(area.get_selected_feature_id(), Some(1));

    area.toggle_feature_multi_select(1);
    assert_eq!(area.multi_selected_feature_ids, BTreeSet::from([0usize]));
    assert_eq!(area.get_selected_feature_id(), Some(0));
}

#[test]
fn feature_multi_select_toggle_clears_when_last_feature_is_removed() {
    let mut dna = DNAsequence::from_sequence("A".repeat(120).as_str()).expect("sequence");
    dna.features_mut()
        .push(make_feature("gene", vec![("label", "gene_1")]));
    let mut area = MainAreaDna::new(dna, None, None);

    area.focus_feature(0);
    area.toggle_feature_multi_select(0);
    assert!(area.multi_selected_feature_ids.is_empty());
    assert_eq!(area.get_selected_feature_id(), None);
    assert_eq!(area.focused_feature_id, None);
}

#[test]
fn clear_multi_select_keeps_primary_feature_selected() {
    let mut dna = DNAsequence::from_sequence("A".repeat(240).as_str()).expect("sequence");
    dna.features_mut()
        .push(make_feature("gene", vec![("label", "gene_1")]));
    dna.features_mut()
        .push(make_feature("gene", vec![("label", "gene_2")]));
    dna.features_mut()
        .push(make_feature("gene", vec![("label", "gene_3")]));
    let mut area = MainAreaDna::new(dna, None, None);

    area.focus_feature(0);
    area.toggle_feature_multi_select(1);
    area.toggle_feature_multi_select(2);
    assert_eq!(
        area.multi_selected_feature_ids,
        BTreeSet::from([0usize, 1usize, 2usize])
    );
    assert_eq!(area.focused_feature_id, Some(2));

    area.clear_multi_select_keep_primary();
    assert_eq!(area.multi_selected_feature_ids, BTreeSet::from([2usize]));
    assert_eq!(area.focused_feature_id, Some(2));
    assert_eq!(area.get_selected_feature_id(), Some(2));
}

#[test]
fn linear_base_status_and_auto_hide_follow_shared_adaptive_routing() {
    let dna = DNAsequence::from_sequence("A".repeat(5000).as_str()).expect("sequence");
    let mut area = MainAreaDna::new(dna, None, None);
    area.last_linear_map_width_px = 1000.0;
    area.show_map = true;
    area.show_sequence = true;
    {
        let mut display = area.dna_display.write().expect("display lock");
        display.set_linear_viewport(0, 300);
        display.set_linear_sequence_helical_letters_enabled(true);
        display
            .set_linear_sequence_letter_layout_mode(LinearSequenceLetterLayoutMode::AutoAdaptive);
        display.set_auto_hide_sequence_panel_when_linear_bases_visible(true);
    }

    let decision = area.linear_base_routing_decision();
    assert_eq!(decision.route_policy, LinearBaseRoutePolicy::Auto);
    assert_ne!(decision.active_mode, LinearBaseRenderMode::Off);
    assert!(area.linear_map_sequence_bases_visible());
    assert!(area.should_auto_hide_sequence_panel());
}

#[test]
fn sequence_panel_layout_config_keeps_visible_panel_resizable() {
    let config = MainAreaDna::sequence_panel_layout_config(800.0, true);
    assert!(config.resizable);
    assert_eq!(config.exact_height_px, None);
    assert_eq!(config.default_height_px, 400.0);
    assert_eq!(config.min_height_px, 200.0);
    assert_eq!(config.max_height_px, 400.0);
}

#[test]
fn sequence_panel_layout_config_collapses_hidden_panel() {
    let config = MainAreaDna::sequence_panel_layout_config(800.0, false);
    assert!(!config.resizable);
    assert_eq!(config.exact_height_px, Some(0.0));
    assert_eq!(config.default_height_px, 0.0);
    assert_eq!(config.min_height_px, 0.0);
    assert_eq!(config.max_height_px, 0.0);
}

#[test]
fn linear_base_auto_hide_uses_width_cache_fallback() {
    let dna = DNAsequence::from_sequence("A".repeat(5000).as_str()).expect("sequence");
    let mut area = MainAreaDna::new(dna, None, None);
    area.last_linear_map_width_px = 0.0;
    area.show_map = true;
    area.show_sequence = true;
    {
        let mut display = area.dna_display.write().expect("display lock");
        display.set_linear_viewport(0, 3000);
        display.set_linear_sequence_helical_letters_enabled(true);
        display
            .set_linear_sequence_letter_layout_mode(LinearSequenceLetterLayoutMode::AutoAdaptive);
        display.set_auto_hide_sequence_panel_when_linear_bases_visible(true);
    }

    let decision = area.linear_base_routing_decision();
    assert_eq!(decision.active_mode, LinearBaseRenderMode::Off);
    assert!(!area.linear_map_sequence_bases_visible());
    assert!(!area.should_auto_hide_sequence_panel());
}

#[test]
fn linear_vertical_pan_is_independent_from_horizontal_viewport() {
    let dna = DNAsequence::from_sequence("A".repeat(5000).as_str()).expect("sequence");
    let area = MainAreaDna::new(dna, None, None);
    area.set_linear_viewport(120, 640);
    area.pan_linear_vertical_viewport(72.0);

    let (start, span, _) = area.current_linear_viewport();
    assert_eq!(start, 120);
    assert_eq!(span, 640);
    assert!((area.current_linear_vertical_offset_px() - 72.0).abs() < 0.001);
}

#[test]
fn linear_vertical_pan_offset_is_clamped() {
    let dna = DNAsequence::from_sequence("A".repeat(400).as_str()).expect("sequence");
    let area = MainAreaDna::new(dna, None, None);
    area.set_linear_vertical_offset_px(50_000.0);
    assert_eq!(area.current_linear_vertical_offset_px(), 10_000.0);
    area.set_linear_vertical_offset_px(-50_000.0);
    assert_eq!(area.current_linear_vertical_offset_px(), -10_000.0);
}

#[test]
fn linear_vertical_fit_delta_centers_content_with_margin() {
    let delta = MainAreaDna::linear_vertical_fit_delta(0.0, 400.0, 80.0, 240.0, 20.0);
    // Viewport center is 200.0, content center is 160.0.
    assert!((delta - 40.0).abs() < 0.001);
}

#[test]
fn linear_vertical_fit_delta_rejects_invalid_ranges() {
    assert_eq!(
        MainAreaDna::linear_vertical_fit_delta(0.0, 0.0, 10.0, 20.0, 10.0),
        0.0
    );
    assert_eq!(
        MainAreaDna::linear_vertical_fit_delta(0.0, 200.0, 20.0, 10.0, 10.0),
        0.0
    );
}

#[test]
fn linear_vertical_fit_delta_clamp_limits_extreme_pan() {
    assert_eq!(
        MainAreaDna::clamp_linear_vertical_fit_delta(12_000.0, 0.0, 400.0),
        800.0
    );
    assert_eq!(
        MainAreaDna::clamp_linear_vertical_fit_delta(-12_000.0, 0.0, 400.0),
        -800.0
    );
}

#[test]
fn view_svg_export_profiles_scale_canvas_and_context() {
    let screen =
        MainAreaDna::view_svg_export_layout(ViewSvgExportProfile::Screen, 1280.0, true, true);
    let wide =
        MainAreaDna::view_svg_export_layout(ViewSvgExportProfile::WideContext, 1280.0, true, true);
    let print = MainAreaDna::view_svg_export_layout(
        ViewSvgExportProfile::PrintA3Landscape,
        1280.0,
        true,
        true,
    );

    assert!(wide.canvas_width_px > screen.canvas_width_px);
    assert!(print.canvas_width_px > wide.canvas_width_px);
    assert!(wide.viewport_span_multiplier > screen.viewport_span_multiplier);
    assert!(print.viewport_span_multiplier > wide.viewport_span_multiplier);
    assert_eq!(print.print_size_mm, Some((420.0, 297.0)));
}

#[test]
fn expand_linear_export_window_expands_and_clamps_to_bounds() {
    assert_eq!(
        MainAreaDna::expand_linear_export_window(100, 200, 1000, 2.0),
        (0, 400)
    );
    assert_eq!(
        MainAreaDna::expand_linear_export_window(850, 200, 1000, 2.0),
        (600, 400)
    );
    assert_eq!(
        MainAreaDna::expand_linear_export_window(10, 0, 0, 4.0),
        (0, 0)
    );
}

#[test]
fn splicing_expert_viewport_id_is_stable_and_feature_scoped() {
    let baseline = MainAreaDna::splicing_expert_viewport_id("seq1", 7);
    assert_eq!(
        baseline,
        MainAreaDna::splicing_expert_viewport_id("seq1", 7)
    );
    assert_ne!(
        baseline,
        MainAreaDna::splicing_expert_viewport_id("seq1", 8)
    );
    assert_ne!(
        baseline,
        MainAreaDna::splicing_expert_viewport_id("seq2", 7)
    );
}

#[test]
fn splicing_expert_window_content_is_wider_than_default_viewport() {
    let default_size = MainAreaDna::splicing_expert_window_default_size();
    let min_size = MainAreaDna::splicing_expert_window_min_size();
    let content_min_size = MainAreaDna::splicing_expert_window_content_min_size();
    assert!(default_size.x < content_min_size.x);
    assert!(min_size.x < content_min_size.x);
    assert!(default_size.x >= min_size.x);
}

#[test]
fn collect_open_auxiliary_window_entries_includes_splicing_and_isoform_windows() {
    let dna = DNAsequence::from_sequence("ACGT").expect("sequence");
    let mut area = MainAreaDna::new(dna, Some("seq1".to_string()), None);
    area.show_splicing_expert_window = true;
    area.splicing_expert_window_view = Some(Arc::new(SplicingExpertView {
        seq_id: "seq1".to_string(),
        target_feature_id: 17,
        scope: SplicingScopePreset::AllOverlappingAnyStrand,
        group_label: "TP73".to_string(),
        strand: "+".to_string(),
        region_start_1based: 1,
        region_end_1based: 4,
        transcript_count: 0,
        unique_exon_count: 0,
        instruction: "splicing".to_string(),
        transcripts: vec![],
        unique_exons: vec![],
        matrix_rows: vec![],
        boundaries: vec![],
        intron_signals: vec![],
        junctions: vec![],
        events: vec![],
    }));
    area.show_isoform_expert_window = true;
    area.isoform_expert_window_panel_id = Some("TA".to_string());
    area.isoform_expert_window_view = Some(Arc::new(IsoformArchitectureExpertView {
        seq_id: "seq1".to_string(),
        panel_id: "TA".to_string(),
        gene_symbol: "TP73".to_string(),
        transcript_geometry_mode: "exon".to_string(),
        panel_source: None,
        region_start_1based: 1,
        region_end_1based: 4,
        instruction: "isoform".to_string(),
        transcript_lanes: vec![],
        protein_lanes: vec![],
        expression_matrix: None,
        warnings: vec![],
    }));

    let entries = area.collect_open_auxiliary_window_entries();
    assert_eq!(entries.len(), 2);
    assert!(
        entries
            .iter()
            .any(|(_, title, _)| title.contains("Splicing Expert - TP73 (seq1)"))
    );
    assert!(
        entries
            .iter()
            .any(|(_, title, _)| title.contains("Isoform Expert - TA (seq1)"))
    );
}

#[test]
fn isoform_expert_window_title_uses_protein_expert_for_transcript_native_source() {
    let view = IsoformArchitectureExpertView {
        seq_id: "seq1".to_string(),
        panel_id: "protein_compare@seq1".to_string(),
        gene_symbol: "TP73".to_string(),
        transcript_geometry_mode: "cds".to_string(),
        panel_source: Some(
            "Transcript-native protein derivation (external protein opinions optional)".to_string(),
        ),
        region_start_1based: 1,
        region_end_1based: 4,
        instruction: "protein".to_string(),
        transcript_lanes: vec![],
        protein_lanes: vec![],
        expression_matrix: None,
        warnings: vec![],
    };
    assert_eq!(
        MainAreaDna::isoform_expert_window_title("protein_compare@seq1", "seq1", &view),
        "Protein Expert - TP73 (seq1)"
    );
}

#[test]
fn isoform_expert_window_title_uses_protein_expert_for_ensembl_backed_source() {
    let view = IsoformArchitectureExpertView {
        seq_id: "seq1".to_string(),
        panel_id: "ENSP00000362111".to_string(),
        gene_symbol: "TP73".to_string(),
        transcript_geometry_mode: "cds".to_string(),
        panel_source: Some(
            "Transcript-native proteins with optional Ensembl opinion ENSP00000362111 (ENST00000378295)"
                .to_string(),
        ),
        region_start_1based: 1,
        region_end_1based: 4,
        instruction: "protein".to_string(),
        transcript_lanes: vec![],
        protein_lanes: vec![],
        expression_matrix: None,
        warnings: vec![],
    };
    assert_eq!(
        MainAreaDna::isoform_expert_window_title("ENSP00000362111", "seq1", &view),
        "Protein Expert - TP73 (seq1)"
    );
}

#[test]
fn collect_open_auxiliary_window_entries_describes_transcript_first_panel_as_protein_expert() {
    let dna = DNAsequence::from_sequence("ACGT").expect("sequence");
    let mut area = MainAreaDna::new(dna, Some("seq1".to_string()), None);
    area.show_isoform_expert_window = true;
    area.isoform_expert_window_panel_id = Some("ENSP00000362111".to_string());
    area.isoform_expert_window_view = Some(Arc::new(IsoformArchitectureExpertView {
        seq_id: "seq1".to_string(),
        panel_id: "ENSP00000362111".to_string(),
        gene_symbol: "TP73".to_string(),
        transcript_geometry_mode: "cds".to_string(),
        panel_source: Some(
            "Transcript-native proteins with optional Ensembl opinion ENSP00000362111 (ENST00000378295)"
                .to_string(),
        ),
        region_start_1based: 1,
        region_end_1based: 4,
        instruction: "protein".to_string(),
        transcript_lanes: vec![],
        protein_lanes: vec![],
        expression_matrix: None,
        warnings: vec![],
    }));

    let entries = area.collect_open_auxiliary_window_entries();
    assert!(
        entries.iter().any(|(_, title, description)| {
            title == "Protein Expert - TP73 (seq1)"
                && description == "Protein Expert 'ENSP00000362111' on 'seq1'"
        }),
        "expected transcript-first protein panel to be described as Protein Expert"
    );
}

#[test]
fn isoform_expert_embed_mode_uses_single_hosted_shell_without_viewport_title_layer() {
    let ctx = egui::Context::default();
    ctx.set_embed_viewports(true);
    let dna = DNAsequence::from_sequence("ACGT").expect("sequence");
    let mut area = MainAreaDna::new(dna, Some("seq1".to_string()), None);
    let view = IsoformArchitectureExpertView {
        seq_id: "seq1".to_string(),
        panel_id: "ENSP00000362111".to_string(),
        gene_symbol: "TP73".to_string(),
        transcript_geometry_mode: "cds".to_string(),
        panel_source: Some(
            "Transcript-native proteins with optional Ensembl opinion ENSP00000362111".to_string(),
        ),
        region_start_1based: 1,
        region_end_1based: 4,
        instruction: "protein".to_string(),
        transcript_lanes: vec![],
        protein_lanes: vec![],
        expression_matrix: None,
        warnings: vec![],
    };
    area.show_isoform_expert_window = true;
    area.isoform_expert_window_panel_id = Some("ENSP00000362111".to_string());
    area.isoform_expert_window_view = Some(Arc::new(view));
    let hosted_layer_id = egui::LayerId::new(
        egui::Order::Middle,
        egui::Id::new("isoform_expert_window_embedded_seq1_ENSP00000362111"),
    );
    let stale_viewport_layer_id = egui::LayerId::new(
        egui::Order::Middle,
        egui::Id::new(MainAreaDna::isoform_expert_viewport_id(
            "seq1",
            "ENSP00000362111",
        )),
    );

    ctx.begin_pass(egui::RawInput::default());
    area.render_isoform_expert_window(&ctx);

    assert!(ctx.memory(|mem| mem.areas().is_visible(&hosted_layer_id)));
    assert!(!ctx.memory(|mem| mem.areas().is_visible(&stale_viewport_layer_id)));
    let _ = ctx.end_pass();
}

#[test]
fn collect_open_auxiliary_window_entries_includes_rna_read_mapping_window() {
    let dna = DNAsequence::from_sequence("ACGT").expect("sequence");
    let mut area = MainAreaDna::new(dna, Some("seq1".to_string()), None);
    area.show_rna_read_mapping_window = true;
    area.rna_read_mapping_window_view = Some(Arc::new(SplicingExpertView {
        seq_id: "seq1".to_string(),
        target_feature_id: 17,
        scope: SplicingScopePreset::AllOverlappingAnyStrand,
        group_label: "TP73".to_string(),
        strand: "+".to_string(),
        region_start_1based: 1,
        region_end_1based: 4,
        transcript_count: 0,
        unique_exon_count: 0,
        instruction: "mapping".to_string(),
        transcripts: vec![],
        unique_exons: vec![],
        matrix_rows: vec![],
        boundaries: vec![],
        intron_signals: vec![],
        junctions: vec![],
        events: vec![],
    }));

    let entries = area.collect_open_auxiliary_window_entries();
    assert_eq!(entries.len(), 1);
    assert!(
        entries
            .iter()
            .any(|(_, title, _)| title.contains("RNA-read Mapping - TP73 (seq1)"))
    );
}

#[test]
fn rna_read_mapping_embedded_window_registers_visible_window_area() {
    let ctx = egui::Context::default();
    let dna = DNAsequence::from_sequence("ACGT").expect("sequence");
    let mut area = MainAreaDna::new(dna, Some("seq1".to_string()), None);
    let view = SplicingExpertView {
        seq_id: "seq1".to_string(),
        target_feature_id: 17,
        scope: SplicingScopePreset::AllOverlappingAnyStrand,
        group_label: "TP73".to_string(),
        strand: "+".to_string(),
        region_start_1based: 1,
        region_end_1based: 4,
        transcript_count: 0,
        unique_exon_count: 0,
        instruction: "mapping".to_string(),
        transcripts: vec![],
        unique_exons: vec![],
        matrix_rows: vec![],
        boundaries: vec![],
        intron_signals: vec![],
        junctions: vec![],
        events: vec![],
    };
    let layer_id = egui::LayerId::new(
        egui::Order::Middle,
        MainAreaDna::rna_read_mapping_embedded_window_id(&view),
    );

    ctx.begin_pass(egui::RawInput::default());
    area.show_rna_read_mapping_window = true;
    area.render_rna_read_mapping_embedded_window_shell(
        &ctx,
        &MainAreaDna::rna_read_mapping_window_title(&view),
        &view,
        MainAreaDna::rna_read_mapping_window_default_size(),
        MainAreaDna::rna_read_mapping_window_content_min_size(),
        false,
        false,
    );
    assert!(ctx.memory(|mem| mem.areas().is_visible(&layer_id)));
    let _ = ctx.end_pass();
}

#[test]
fn rna_read_mapping_embedded_window_focus_request_uses_foreground_order() {
    let ctx = egui::Context::default();
    let dna = DNAsequence::from_sequence("ACGT").expect("sequence");
    let mut area = MainAreaDna::new(dna, Some("seq1".to_string()), None);
    let view = SplicingExpertView {
        seq_id: "seq1".to_string(),
        target_feature_id: 17,
        scope: SplicingScopePreset::AllOverlappingAnyStrand,
        group_label: "TP73".to_string(),
        strand: "+".to_string(),
        region_start_1based: 1,
        region_end_1based: 4,
        transcript_count: 0,
        unique_exon_count: 0,
        instruction: "mapping".to_string(),
        transcripts: vec![],
        unique_exons: vec![],
        matrix_rows: vec![],
        boundaries: vec![],
        intron_signals: vec![],
        junctions: vec![],
        events: vec![],
    };
    let layer_id = egui::LayerId::new(
        egui::Order::Foreground,
        MainAreaDna::rna_read_mapping_embedded_window_id(&view),
    );

    ctx.begin_pass(egui::RawInput::default());
    area.show_rna_read_mapping_window = true;
    area.render_rna_read_mapping_embedded_window_shell(
        &ctx,
        &MainAreaDna::rna_read_mapping_window_title(&view),
        &view,
        MainAreaDna::rna_read_mapping_window_default_size(),
        MainAreaDna::rna_read_mapping_window_content_min_size(),
        false,
        true,
    );
    assert!(ctx.memory(|mem| mem.areas().is_visible(&layer_id)));
    let _ = ctx.end_pass();
}

#[test]
fn rna_read_mapping_embedded_window_clears_legacy_title_bar_area() {
    let ctx = egui::Context::default();
    let dna = DNAsequence::from_sequence("ACGT").expect("sequence");
    let mut area = MainAreaDna::new(dna, Some("seq1".to_string()), None);
    let view = SplicingExpertView {
        seq_id: "seq1".to_string(),
        target_feature_id: 17,
        scope: SplicingScopePreset::AllOverlappingAnyStrand,
        group_label: "TP73".to_string(),
        strand: "+".to_string(),
        region_start_1based: 1,
        region_end_1based: 4,
        transcript_count: 0,
        unique_exon_count: 0,
        instruction: "mapping".to_string(),
        transcripts: vec![],
        unique_exons: vec![],
        matrix_rows: vec![],
        boundaries: vec![],
        intron_signals: vec![],
        junctions: vec![],
        events: vec![],
    };
    let title = MainAreaDna::rna_read_mapping_window_title(&view);
    let hosted_layer_id = egui::LayerId::new(
        egui::Order::Middle,
        MainAreaDna::rna_read_mapping_embedded_window_id(&view),
    );
    let stale_title_layer_id = MainAreaDna::stale_auxiliary_window_title_layer_id(&title);

    ctx.begin_pass(egui::RawInput::default());
    crate::egui_compat::show_legacy_layer_for_tests(&ctx, stale_title_layer_id, |ui| {
        ui.label("legacy mapping title shell");
    });
    assert!(ctx.memory(|mem| mem.areas().is_visible(&stale_title_layer_id)));
    let _ = ctx.end_pass();

    ctx.begin_pass(egui::RawInput::default());
    area.show_rna_read_mapping_window = true;
    area.render_rna_read_mapping_embedded_window_shell(
        &ctx,
        &title,
        &view,
        MainAreaDna::rna_read_mapping_window_default_size(),
        MainAreaDna::rna_read_mapping_window_content_min_size(),
        false,
        false,
    );
    assert!(ctx.memory(|mem| mem.areas().is_visible(&hosted_layer_id)));
    let _ = ctx.end_pass();
    assert!(!ctx.memory(|mem| mem.areas().is_visible(&stale_title_layer_id)));
}

#[test]
fn rna_read_mapping_embed_mode_uses_single_hosted_shell_without_viewport_title_layer() {
    let ctx = egui::Context::default();
    ctx.set_embed_viewports(true);
    let dna = DNAsequence::from_sequence("ACGT").expect("sequence");
    let mut area = MainAreaDna::new(dna, Some("seq1".to_string()), None);
    area.seed_rna_read_mapping_window_for_tests("seq1", 17, "TP73");
    let view = area
        .rna_read_mapping_window_view
        .as_deref()
        .expect("seeded RNA-read Mapping view")
        .clone();
    let title = MainAreaDna::rna_read_mapping_window_title(&view);
    let hosted_layer_id = egui::LayerId::new(
        egui::Order::Middle,
        MainAreaDna::rna_read_mapping_embedded_window_id(&view),
    );
    let stale_viewport_layer_id = egui::LayerId::new(
        egui::Order::Middle,
        egui::Id::new(MainAreaDna::rna_read_mapping_viewport_id(
            &view.seq_id,
            view.target_feature_id,
        )),
    );
    let stale_title_layer_id = MainAreaDna::stale_auxiliary_window_title_layer_id(&title);

    ctx.begin_pass(egui::RawInput::default());
    area.render_rna_read_mapping_window(&ctx);

    assert!(ctx.memory(|mem| mem.areas().is_visible(&hosted_layer_id)));
    assert!(!ctx.memory(|mem| mem.areas().is_visible(&stale_viewport_layer_id)));
    assert!(!ctx.memory(|mem| mem.areas().is_visible(&stale_title_layer_id)));
    let _ = ctx.end_pass();
}

#[test]
fn splicing_expert_embed_mode_uses_single_hosted_shell_without_viewport_title_layer() {
    let ctx = egui::Context::default();
    ctx.set_embed_viewports(true);
    let dna = DNAsequence::from_sequence("ACGT").expect("sequence");
    let mut area = MainAreaDna::new(dna, Some("seq1".to_string()), None);
    area.seed_splicing_expert_window_for_tests("seq1", 17, "TP73");
    let view = area
        .splicing_expert_window_view
        .as_deref()
        .expect("seeded Splicing Expert view")
        .clone();
    let title = MainAreaDna::splicing_expert_window_title(&view);
    let hosted_layer_id = egui::LayerId::new(
        egui::Order::Middle,
        MainAreaDna::splicing_expert_embedded_window_id(&view),
    );
    let stale_viewport_layer_id = egui::LayerId::new(
        egui::Order::Middle,
        egui::Id::new(MainAreaDna::splicing_expert_viewport_id(
            &view.seq_id,
            view.target_feature_id,
        )),
    );
    let stale_title_layer_id = MainAreaDna::stale_auxiliary_window_title_layer_id(&title);

    ctx.begin_pass(egui::RawInput::default());
    area.render_splicing_expert_window(&ctx);

    assert!(ctx.memory(|mem| mem.areas().is_visible(&hosted_layer_id)));
    assert!(!ctx.memory(|mem| mem.areas().is_visible(&stale_viewport_layer_id)));
    assert!(!ctx.memory(|mem| mem.areas().is_visible(&stale_title_layer_id)));
    let _ = ctx.end_pass();
}

#[test]
fn rna_read_mapping_embedded_window_renders_intro_only_once() {
    let ctx = egui::Context::default();
    let dna = DNAsequence::from_sequence("ACGT").expect("sequence");
    let mut area = MainAreaDna::new(dna, Some("seq1".to_string()), None);
    let view = SplicingExpertView {
        seq_id: "seq1".to_string(),
        target_feature_id: 17,
        scope: SplicingScopePreset::AllOverlappingAnyStrand,
        group_label: "TP73".to_string(),
        strand: "+".to_string(),
        region_start_1based: 1,
        region_end_1based: 4,
        transcript_count: 0,
        unique_exon_count: 0,
        instruction: "mapping".to_string(),
        transcripts: vec![],
        unique_exons: vec![],
        matrix_rows: vec![],
        boundaries: vec![],
        intron_signals: vec![],
        junctions: vec![],
        events: vec![],
    };
    let intro_summary = format!(
        "Locus: {} | feature n-{} | group '{}' | transcripts={}",
        view.seq_id, view.target_feature_id, view.group_label, view.transcript_count
    );

    ctx.begin_pass(egui::RawInput::default());
    crate::egui_compat::show_central_panel_for_test_context(
        &ctx,
        egui::CentralPanel::default(),
        |ui| {
            area.render_rna_read_mapping_window_body(&ctx, ui, &view, false);
        },
    );
    let texts = collect_pass_texts(&ctx);
    assert_eq!(
        texts
            .iter()
            .filter(|text| text.contains("Workspace guide"))
            .count(),
        1
    );
    assert_eq!(
        texts
            .iter()
            .filter(|text| text.contains(&intro_summary))
            .count(),
        1
    );
}

#[test]
fn rna_read_preview_table_row_budget_grows_for_short_lists() {
    assert_eq!(MainAreaDna::target_rna_read_preview_body_rows(0), 12);
    assert_eq!(MainAreaDna::target_rna_read_preview_body_rows(12), 12);
    assert_eq!(MainAreaDna::target_rna_read_preview_body_rows(13), 13);
    assert_eq!(MainAreaDna::target_rna_read_preview_body_rows(16), 16);
    assert_eq!(MainAreaDna::target_rna_read_preview_body_rows(24), 16);
}

#[test]
fn collect_open_auxiliary_window_entries_includes_dotplot_window() {
    let dna = DNAsequence::from_sequence("ACGT").expect("sequence");
    let mut area = MainAreaDna::new(dna, Some("seq1".to_string()), None);
    area.show_dotplot_window = true;
    let entries = area.collect_open_auxiliary_window_entries();
    assert_eq!(entries.len(), 1);
    assert_eq!(
        entries[0].0,
        MainAreaDna::dotplot_viewport_id("seq1"),
        "dotplot viewport id should be deterministic for sequence scope"
    );
    assert!(
        entries[0].1.contains("Dotplot - seq1"),
        "dotplot title should be listed in auxiliary windows"
    );
}

#[test]
fn dotplot_embed_mode_uses_single_hosted_shell_without_viewport_title_layer() {
    let ctx = egui::Context::default();
    ctx.set_embed_viewports(true);
    let dna = DNAsequence::from_sequence("ACGT").expect("sequence");
    let mut area = MainAreaDna::new(dna, Some("seq1".to_string()), None);
    area.show_dotplot_window = true;
    let hosted_layer_id = egui::LayerId::new(
        egui::Order::Middle,
        egui::Id::new("dotplot_window_embedded_seq1"),
    );
    let stale_viewport_layer_id = egui::LayerId::new(
        egui::Order::Middle,
        egui::Id::new(MainAreaDna::dotplot_viewport_id("seq1")),
    );

    ctx.begin_pass(egui::RawInput::default());
    area.render_dotplot_window(&ctx);

    assert!(ctx.memory(|mem| mem.areas().is_visible(&hosted_layer_id)));
    assert!(!ctx.memory(|mem| mem.areas().is_visible(&stale_viewport_layer_id)));
    let _ = ctx.end_pass();
}

#[test]
fn variant_followup_embed_mode_uses_single_hosted_shell_without_viewport_title_layer() {
    let ctx = egui::Context::default();
    ctx.set_embed_viewports(true);
    let dna = DNAsequence::from_sequence("ACGT").expect("sequence");
    let mut area = MainAreaDna::new(dna, Some("seq1".to_string()), None);
    area.show_variant_followup_window = true;
    area.variant_followup_window_pending_initial_render = false;
    area.variant_followup_ui.source_seq_id = "seq1".to_string();
    area.variant_followup_ui.source_feature_id = Some(17);
    area.variant_followup_ui.gene_label = "TP73".to_string();
    let hosted_layer_id = egui::LayerId::new(
        egui::Order::Middle,
        egui::Id::new("variant_followup_window_embedded_seq1_17"),
    );
    let stale_viewport_layer_id = egui::LayerId::new(
        egui::Order::Middle,
        egui::Id::new(egui::ViewportId::from_hash_of((
            "variant_followup_viewport",
            "seq1",
            17usize,
        ))),
    );

    ctx.begin_pass(egui::RawInput::default());
    area.render_variant_followup_window(&ctx);

    assert!(ctx.memory(|mem| mem.areas().is_visible(&hosted_layer_id)));
    assert!(!ctx.memory(|mem| mem.areas().is_visible(&stale_viewport_layer_id)));
    let _ = ctx.end_pass();
}

#[test]
fn dotplot_window_uses_compact_layout_without_loaded_payload() {
    let compact_default = MainAreaDna::dotplot_window_default_size(false);
    let loaded_default = MainAreaDna::dotplot_window_default_size(true);
    let compact_min = MainAreaDna::dotplot_window_min_size(false);
    let loaded_min = MainAreaDna::dotplot_window_min_size(true);
    let compact_content = MainAreaDna::dotplot_window_content_min_size(false);
    let loaded_content = MainAreaDna::dotplot_window_content_min_size(true);

    assert_eq!(compact_default.x, loaded_default.x);
    assert!(compact_default.y < loaded_default.y);
    assert!(compact_min.y < loaded_min.y);
    assert_eq!(compact_content.x, loaded_content.x);
    assert_eq!(compact_content.y, 0.0);
    assert!(compact_content.y < loaded_content.y);
}

#[test]
fn current_engine_ops_state_records_primary_map_mode() {
    let dna = DNAsequence::from_sequence("ACGT").unwrap();
    let mut area = MainAreaDna::new(dna, None, None);
    area.primary_map_mode = PrimaryMapMode::Splicing;
    let state = area.current_engine_ops_state();
    assert_eq!(state.primary_map_mode, PrimaryMapMode::Splicing);
}

#[test]
fn current_engine_ops_state_records_dna_presentation_mode() {
    let dna = DNAsequence::from_sequence("ACGT").unwrap();
    let mut area = MainAreaDna::new(dna, None, None);
    area.dna_presentation_mode = DnaPresentationMode::Cdna;
    let state = area.current_engine_ops_state();
    assert_eq!(state.dna_presentation_mode, DnaPresentationMode::Cdna);
}

#[test]
fn current_engine_ops_state_records_contextual_transcript_toggle() {
    let dna = DNAsequence::from_sequence("ACGT").unwrap();
    let mut area = MainAreaDna::new(dna, None, None);
    area.show_all_contextual_transcripts = true;
    let state = area.current_engine_ops_state();
    assert!(state.show_all_contextual_transcripts);
}

#[test]
fn current_engine_ops_state_records_extended_top_panel_height() {
    let dna = DNAsequence::from_sequence("ACGT").unwrap();
    let mut area = MainAreaDna::new(dna, None, None);
    area.extended_top_panel_height_px = 512.0;
    let state = area.current_engine_ops_state();
    assert!((state.extended_top_panel_height_px - 512.0).abs() <= f32::EPSILON);
}

#[test]
fn current_engine_ops_state_records_pcr_designer_mode() {
    let dna = DNAsequence::from_sequence("ACGT").unwrap();
    let mut area = MainAreaDna::new(dna, None, None);
    area.pcr_designer_mode = PcrDesignerMode::QpcrAssays;
    let state = area.current_engine_ops_state();
    assert_eq!(state.pcr_designer_mode, PcrDesignerMode::QpcrAssays);
}

#[test]
fn primary_map_mode_defaults_when_missing_in_serialized_engine_ops_state() {
    let dna = DNAsequence::from_sequence("ACGT").unwrap();
    let area = MainAreaDna::new(dna, None, None);
    let mut value = serde_json::to_value(area.current_engine_ops_state()).unwrap();
    value.as_object_mut().unwrap().remove("primary_map_mode");
    let decoded: super::EngineOpsUiState = serde_json::from_value(value).unwrap();
    assert_eq!(decoded.primary_map_mode, PrimaryMapMode::Standard);
}

#[test]
fn dna_presentation_mode_defaults_when_missing_in_serialized_engine_ops_state() {
    let dna = DNAsequence::from_sequence("ACGT").unwrap();
    let area = MainAreaDna::new(dna, None, None);
    let mut value = serde_json::to_value(area.current_engine_ops_state()).unwrap();
    value
        .as_object_mut()
        .unwrap()
        .remove("dna_presentation_mode");
    let decoded: super::EngineOpsUiState = serde_json::from_value(value).unwrap();
    assert_eq!(decoded.dna_presentation_mode, DnaPresentationMode::Region);
}

#[test]
fn pcr_designer_mode_defaults_when_missing_in_serialized_engine_ops_state() {
    let dna = DNAsequence::from_sequence("ACGT").unwrap();
    let area = MainAreaDna::new(dna, None, None);
    let mut value = serde_json::to_value(area.current_engine_ops_state()).unwrap();
    value.as_object_mut().unwrap().remove("pcr_designer_mode");
    let decoded: super::EngineOpsUiState = serde_json::from_value(value).unwrap();
    assert_eq!(decoded.pcr_designer_mode, PcrDesignerMode::PrimerPairs);
}

#[test]
fn contextual_transcript_toggle_defaults_when_missing_in_serialized_engine_ops_state() {
    let dna = DNAsequence::from_sequence("ACGT").unwrap();
    let area = MainAreaDna::new(dna, None, None);
    let mut value = serde_json::to_value(area.current_engine_ops_state()).unwrap();
    value
        .as_object_mut()
        .unwrap()
        .remove("show_all_contextual_transcripts");
    let decoded: super::EngineOpsUiState = serde_json::from_value(value).unwrap();
    assert!(!decoded.show_all_contextual_transcripts);
}

#[test]
fn extended_top_panel_height_defaults_when_missing_in_serialized_engine_ops_state() {
    let dna = DNAsequence::from_sequence("ACGT").unwrap();
    let area = MainAreaDna::new(dna, None, None);
    let mut value = serde_json::to_value(area.current_engine_ops_state()).unwrap();
    value
        .as_object_mut()
        .unwrap()
        .remove("extended_top_panel_height_px");
    let decoded: super::EngineOpsUiState = serde_json::from_value(value).unwrap();
    assert!(
        (decoded.extended_top_panel_height_px - super::default_extended_top_panel_height_px())
            .abs()
            <= f32::EPSILON
    );
}

#[test]
fn simple_pcr_distance_defaults_when_missing_in_serialized_engine_ops_state() {
    let dna = DNAsequence::from_sequence("ACGT").unwrap();
    let area = MainAreaDna::new(dna, None, None);
    let mut value = serde_json::to_value(area.current_engine_ops_state()).unwrap();
    value
        .as_object_mut()
        .unwrap()
        .remove("simple_pcr_max_primer_distance_bp");
    let decoded: super::EngineOpsUiState = serde_json::from_value(value).unwrap();
    assert_eq!(
        decoded.simple_pcr_max_primer_distance_bp,
        super::default_simple_pcr_max_primer_distance_bp()
    );
}

#[test]
fn rmsk_materialization_controls_default_when_missing_in_serialized_engine_ops_state() {
    let dna = DNAsequence::from_sequence("ACGT").unwrap();
    let area = MainAreaDna::new(dna, None, None);
    let mut value = serde_json::to_value(area.current_engine_ops_state()).unwrap();
    let object = value.as_object_mut().unwrap();
    object.remove("rmsk_index_path");
    object.remove("rmsk_max_features");
    object.remove("rmsk_append_features");
    let decoded: super::EngineOpsUiState = serde_json::from_value(value).unwrap();
    assert_eq!(
        decoded.rmsk_index_path,
        crate::ucsc_rmsk::DEFAULT_UCSC_RMSK_INDEX_PATH
    );
    assert_eq!(decoded.rmsk_max_features, "5000");
    assert!(!decoded.rmsk_append_features);
}

#[test]
fn clamp_extended_top_panel_height_respects_bounds() {
    let max_height = MainAreaDna::max_extended_top_panel_height(900.0, true);
    assert!(
        (MainAreaDna::clamp_extended_top_panel_height(0.0, 900.0, true)
            - super::EXTENDED_TOP_PANEL_MIN_HEIGHT_PX)
            .abs()
            <= f32::EPSILON
    );
    assert!(
        (MainAreaDna::clamp_extended_top_panel_height(9_999.0, 900.0, true) - max_height).abs()
            <= f32::EPSILON
    );
}

#[test]
fn sync_contextual_transcript_visibility_filter_is_mode_aware() {
    let dna = DNAsequence::from_sequence("ACGT").unwrap();
    let mut area = MainAreaDna::new(dna, None, None);

    area.dna_presentation_mode = DnaPresentationMode::Region;
    area.show_all_contextual_transcripts = false;
    area.sync_contextual_transcript_visibility_filter();
    assert!(
        area.dna_display
            .read()
            .expect("DNA display lock poisoned")
            .show_contextual_transcript_features()
    );

    area.dna_presentation_mode = DnaPresentationMode::Cdna;
    area.show_all_contextual_transcripts = false;
    area.sync_contextual_transcript_visibility_filter();
    assert!(
        !area
            .dna_display
            .read()
            .expect("DNA display lock poisoned")
            .show_contextual_transcript_features()
    );

    area.show_all_contextual_transcripts = true;
    area.sync_contextual_transcript_visibility_filter();
    assert!(
        area.dna_display
            .read()
            .expect("DNA display lock poisoned")
            .show_contextual_transcript_features()
    );
}

#[test]
fn dna_presentation_mode_policy_gates_toolbar_groups() {
    assert!(!DnaPresentationMode::Genome.allows_roi_tools());
    assert!(!DnaPresentationMode::Genome.allows_derivation_tools());
    assert!(!DnaPresentationMode::Genome.allows_engine_shell_panels());

    assert!(DnaPresentationMode::Chromosomal.allows_roi_tools());
    assert!(!DnaPresentationMode::Chromosomal.allows_derivation_tools());
    assert!(!DnaPresentationMode::Chromosomal.allows_engine_shell_panels());

    assert!(DnaPresentationMode::Region.allows_roi_tools());
    assert!(DnaPresentationMode::Region.allows_derivation_tools());
    assert!(DnaPresentationMode::Region.allows_engine_shell_panels());

    assert!(DnaPresentationMode::Gene.allows_roi_tools());
    assert!(DnaPresentationMode::Gene.allows_derivation_tools());
    assert!(DnaPresentationMode::Gene.allows_engine_shell_panels());

    assert!(DnaPresentationMode::Cdna.allows_roi_tools());
    assert!(DnaPresentationMode::Cdna.allows_derivation_tools());
    assert!(DnaPresentationMode::Cdna.allows_engine_shell_panels());
}

#[test]
fn dotplot_ui_defaults_when_missing_in_serialized_engine_ops_state() {
    let dna = DNAsequence::from_sequence("ACGT").unwrap();
    let area = MainAreaDna::new(dna, None, None);
    let mut value = serde_json::to_value(area.current_engine_ops_state()).unwrap();
    value.as_object_mut().unwrap().remove("dotplot_ui");
    let decoded: super::EngineOpsUiState = serde_json::from_value(value).unwrap();
    assert_eq!(decoded.dotplot_ui.half_window_bp, "");
    assert_eq!(decoded.dotplot_ui.word_size, "7");
    assert_eq!(decoded.dotplot_ui.step_bp, "1");
    assert_eq!(decoded.dotplot_ui.max_mismatches, "0");
    assert_eq!(
        decoded.dotplot_ui.mode,
        crate::engine::DotplotMode::SelfReverseComplement
    );
}

#[test]
fn dotplot_half_window_default_normalizes_to_active_sequence_length() {
    let dna = DNAsequence::from_sequence("ACGTAC").unwrap();
    let mut area = MainAreaDna::new(dna, Some("seq1".to_string()), None);
    area.dotplot_ui.half_window_bp.clear();
    assert_eq!(area.dotplot_ui.half_window_bp, "");
    assert!(area.normalize_dotplot_half_window_default_if_needed());
    assert_eq!(area.dotplot_ui.half_window_bp, "6");
    assert!(!area.normalize_dotplot_half_window_default_if_needed());
}

#[test]
fn dotplot_half_window_manual_legacy_like_value_is_preserved() {
    let dna = DNAsequence::from_sequence("ACGTAC").unwrap();
    let mut area = MainAreaDna::new(dna, None, None);
    area.dotplot_ui.half_window_bp = "500".to_string();
    assert!(!area.normalize_dotplot_half_window_default_if_needed());
    assert_eq!(area.dotplot_ui.half_window_bp, "500");
}

#[test]
fn rna_reads_ui_defaults_when_missing_in_serialized_engine_ops_state() {
    let dna = DNAsequence::from_sequence("ACGT").unwrap();
    let area = MainAreaDna::new(dna, None, None);
    let mut value = serde_json::to_value(area.current_engine_ops_state()).unwrap();
    value.as_object_mut().unwrap().remove("rna_reads_ui");
    let decoded: super::EngineOpsUiState = serde_json::from_value(value).unwrap();
    assert_eq!(decoded.rna_reads_ui.kmer_len, "10");
    assert_eq!(decoded.rna_reads_ui.seed_stride_bp, "1");
    assert_eq!(decoded.rna_reads_ui.min_seed_hit_fraction, "0.30");
    assert_eq!(decoded.rna_reads_ui.min_unique_matched_kmers, "12");
    assert_eq!(decoded.rna_reads_ui.min_chain_consistency_fraction, "0.40");
    assert_eq!(decoded.rna_reads_ui.max_median_transcript_gap, "4.0");
    assert_eq!(decoded.rna_reads_ui.min_confirmed_exon_transitions, "1");
    assert_eq!(decoded.rna_reads_ui.min_transition_support_fraction, "0.05");
    assert_eq!(
        decoded.rna_reads_ui.origin_mode,
        crate::engine::RnaReadOriginMode::SingleGene
    );
    assert!(decoded.rna_reads_ui.target_gene_ids.is_empty());
    assert!(!decoded.rna_reads_ui.roi_seed_capture_enabled);
    assert_eq!(
        decoded.rna_reads_ui.scope,
        SplicingScopePreset::AllOverlappingAnyStrand
    );
    assert_eq!(
        decoded.rna_reads_ui.report_mode,
        crate::engine::RnaReadReportMode::Full
    );
    assert!(decoded.rna_reads_ui.checkpoint_path.is_empty());
    assert_eq!(decoded.rna_reads_ui.checkpoint_every_reads, "10000");
    assert!(!decoded.rna_reads_ui.resume_from_checkpoint);
}

#[test]
fn rna_read_evidence_ui_defaults_when_missing_in_serialized_engine_ops_state() {
    let dna = DNAsequence::from_sequence("ACGT").unwrap();
    let area = MainAreaDna::new(dna, None, None);
    let mut value = serde_json::to_value(area.current_engine_ops_state()).unwrap();
    value
        .as_object_mut()
        .unwrap()
        .remove("rna_read_evidence_ui");
    let decoded: super::EngineOpsUiState = serde_json::from_value(value).unwrap();
    assert!(decoded.rna_read_evidence_ui.selected_report_id.is_empty());
    assert!(decoded.rna_read_evidence_ui.score_density_use_log_scale);
}

#[test]
fn rna_read_concatemer_ui_defaults_when_missing_in_serialized_engine_ops_state() {
    let dna = DNAsequence::from_sequence("ACGT").unwrap();
    let area = MainAreaDna::new(dna, None, None);
    let mut value = serde_json::to_value(area.current_engine_ops_state()).unwrap();
    value
        .as_object_mut()
        .unwrap()
        .remove("rna_read_concatemer_ui");
    let decoded: super::EngineOpsUiState = serde_json::from_value(value).unwrap();
    assert_eq!(
        decoded.rna_read_concatemer_ui.selection,
        RnaReadHitSelection::Aligned
    );
    assert_eq!(
        decoded.rna_read_concatemer_ui.subset_mode,
        RnaReadConcatemerSubsetMode::AllMatchingRows
    );
    assert_eq!(decoded.rna_read_concatemer_ui.limit, "50");
}

#[test]
fn current_rna_read_concatemer_request_uses_selected_rows_and_multiline_paths() {
    let dna = DNAsequence::from_sequence("ACGT").unwrap();
    let mut area = MainAreaDna::new(dna, None, None);
    area.rna_read_concatemer_ui.selection = RnaReadHitSelection::Aligned;
    area.rna_read_concatemer_ui.subset_mode = RnaReadConcatemerSubsetMode::SelectedRowsOnly;
    area.rna_read_concatemer_ui.limit = "25".to_string();
    area.rna_read_concatemer_ui.adapter_fasta_path = "/tmp/adapters.fa".to_string();
    area.rna_read_concatemer_ui.transcript_fasta_paths_text =
        "/tmp/cdna.fa.gz\n\n/tmp/ncrna.fa.gz".to_string();
    area.rna_read_concatemer_ui.transcript_index_paths_text =
        "/tmp/human.index.json\n/tmp/extra.index.json".to_string();
    area.rna_read_concatemer_ui.fragment_max_parts = "0".to_string();
    area.rna_seed_selected_record_indices = [8usize, 2usize].into_iter().collect();

    let (selection, limit, selected_record_indices, settings) =
        area.current_rna_read_concatemer_request().unwrap();

    assert_eq!(selection, RnaReadHitSelection::Aligned);
    assert_eq!(limit, 25);
    assert_eq!(selected_record_indices, vec![2, 8]);
    assert_eq!(
        settings.adapter_fasta_path.as_deref(),
        Some("/tmp/adapters.fa")
    );
    assert_eq!(
        settings.transcript_fasta_paths,
        vec![
            "/tmp/cdna.fa.gz".to_string(),
            "/tmp/ncrna.fa.gz".to_string()
        ]
    );
    assert_eq!(
        settings.transcript_index_paths,
        vec![
            "/tmp/human.index.json".to_string(),
            "/tmp/extra.index.json".to_string()
        ]
    );
    assert_eq!(settings.fragment_max_parts, 0);
}

#[test]
fn rna_reads_ui_checkpoint_defaults_when_missing_from_partial_state() {
    let dna = DNAsequence::from_sequence("ACGT").unwrap();
    let area = MainAreaDna::new(dna, None, None);
    let mut value = serde_json::to_value(area.current_engine_ops_state()).unwrap();
    let rna = value
        .as_object_mut()
        .unwrap()
        .get_mut("rna_reads_ui")
        .and_then(|v| v.as_object_mut())
        .expect("rna_reads_ui object");
    rna.remove("seed_stride_bp");
    rna.remove("report_mode");
    rna.remove("checkpoint_path");
    rna.remove("checkpoint_every_reads");
    rna.remove("resume_from_checkpoint");
    let decoded: super::EngineOpsUiState = serde_json::from_value(value).unwrap();
    assert_eq!(decoded.rna_reads_ui.seed_stride_bp, "1");
    assert_eq!(
        decoded.rna_reads_ui.report_mode,
        crate::engine::RnaReadReportMode::Full
    );
    assert!(decoded.rna_reads_ui.checkpoint_path.is_empty());
    assert_eq!(decoded.rna_reads_ui.checkpoint_every_reads, "10000");
    assert!(!decoded.rna_reads_ui.resume_from_checkpoint);
}

#[test]
fn parse_rna_target_gene_ids_splits_and_deduplicates_case_insensitive() {
    let genes = MainAreaDna::parse_rna_target_gene_ids("TP73, tp53 TP73;  tp63");
    assert_eq!(
        genes,
        vec!["TP73".to_string(), "tp53".to_string(), "tp63".to_string()]
    );
}

#[test]
fn rna_read_scope_selection_label_mentions_target_gene_when_applicable() {
    let view = SplicingExpertView {
        seq_id: "seq1".to_string(),
        target_feature_id: 17,
        scope: SplicingScopePreset::AllOverlappingAnyStrand,
        group_label: "TP73".to_string(),
        strand: "+".to_string(),
        region_start_1based: 1,
        region_end_1based: 4,
        transcript_count: 0,
        unique_exon_count: 0,
        instruction: String::new(),
        transcripts: vec![],
        unique_exons: vec![],
        matrix_rows: vec![],
        boundaries: vec![],
        intron_signals: vec![],
        junctions: vec![],
        events: vec![],
    };
    assert_eq!(
        MainAreaDna::rna_read_scope_selection_label(
            &view,
            SplicingScopePreset::TargetGroupTargetStrand
        ),
        "Target gene/group TP73 (target-gene strand)"
    );
    assert_eq!(
        MainAreaDna::rna_read_scope_selection_label(
            &view,
            SplicingScopePreset::AllOverlappingAnyStrand
        ),
        "All overlapping genes incl. antisense (any strand)"
    );
}

#[test]
fn rna_read_mapping_button_hover_text_guides_when_no_feature_is_selected() {
    let hover = MainAreaDna::rna_read_mapping_button_hover_text(None);
    assert!(hover.contains("Select an mRNA"));
    assert!(hover.contains("feature first"));
}

#[test]
fn rna_read_scope_cell_included_matches_scope_flags() {
    assert!(MainAreaDna::rna_read_scope_cell_included(
        SplicingScopePreset::AllOverlappingAnyStrand,
        false,
        false,
    ));
    assert!(!MainAreaDna::rna_read_scope_cell_included(
        SplicingScopePreset::TargetGroupAnyStrand,
        false,
        true,
    ));
    assert!(!MainAreaDna::rna_read_scope_cell_included(
        SplicingScopePreset::AllOverlappingTargetStrand,
        true,
        false,
    ));
    assert!(MainAreaDna::rna_read_scope_cell_included(
        SplicingScopePreset::TargetGroupTargetStrand,
        true,
        true,
    ));
}

#[test]
fn rna_read_gene_scope_summary_mentions_explicit_sparse_gene_list() {
    let view = SplicingExpertView {
        seq_id: "seq1".to_string(),
        target_feature_id: 17,
        scope: SplicingScopePreset::AllOverlappingAnyStrand,
        group_label: "TP73".to_string(),
        strand: "+".to_string(),
        region_start_1based: 10,
        region_end_1based: 50,
        transcript_count: 0,
        unique_exon_count: 0,
        instruction: String::new(),
        transcripts: vec![],
        unique_exons: vec![],
        matrix_rows: vec![],
        boundaries: vec![],
        intron_signals: vec![],
        junctions: vec![],
        events: vec![],
    };
    let (summary, _color) = MainAreaDna::rna_read_gene_scope_summary(
        &view,
        SplicingScopePreset::AllOverlappingTargetStrand,
        RnaReadOriginMode::MultiGeneSparse,
        "TP73, TP53",
    );
    assert!(
        summary
            .contains("all genes overlapping seq1:10..50 on the selected target gene/group strand")
    );
    assert!(summary.contains("TP73, TP53"));
}

#[test]
fn latest_matching_rna_read_report_id_prefers_most_recent_summary() {
    let summaries = vec![
        RnaReadInterpretationReportSummary {
            report_id: "older".to_string(),
            generated_at_unix_ms: 10,
            ..RnaReadInterpretationReportSummary {
                report_id: String::new(),
                generated_at_unix_ms: 0,
                report_mode: RnaReadReportMode::Full,
                seq_id: "seq1".to_string(),
                op_id: None,
                run_id: None,
                profile: RnaReadInterpretationProfile::NanoporeCdnaV1,
                input_path: String::new(),
                input_format: RnaReadInputFormat::Fasta,
                seed_feature_id: 7,
                scope: SplicingScopePreset::TargetGroupTargetStrand,
                origin_mode: RnaReadOriginMode::SingleGene,
                input_orientation_mode: "cdna_oriented".to_string(),
                input_orientation_label: "cDNA-oriented".to_string(),
                target_gene_count: 0,
                roi_seed_capture_enabled: false,
                read_count_total: 0,
                read_count_seed_passed: 0,
                read_count_aligned: 0,
                retained_count_msa_eligible: 0,
            }
        },
        RnaReadInterpretationReportSummary {
            report_id: "newer".to_string(),
            generated_at_unix_ms: 20,
            ..RnaReadInterpretationReportSummary {
                report_id: String::new(),
                generated_at_unix_ms: 0,
                report_mode: RnaReadReportMode::Full,
                seq_id: "seq1".to_string(),
                op_id: None,
                run_id: None,
                profile: RnaReadInterpretationProfile::NanoporeCdnaV1,
                input_path: String::new(),
                input_format: RnaReadInputFormat::Fasta,
                seed_feature_id: 7,
                scope: SplicingScopePreset::TargetGroupTargetStrand,
                origin_mode: RnaReadOriginMode::SingleGene,
                input_orientation_mode: "cdna_oriented".to_string(),
                input_orientation_label: "cDNA-oriented".to_string(),
                target_gene_count: 0,
                roi_seed_capture_enabled: false,
                read_count_total: 0,
                read_count_seed_passed: 0,
                read_count_aligned: 0,
                retained_count_msa_eligible: 0,
            }
        },
    ];
    assert_eq!(
        MainAreaDna::latest_matching_rna_read_report_id(&summaries).as_deref(),
        Some("older"),
        "summaries are expected newest-first before selection"
    );
}

#[test]
fn open_rna_read_mapping_workspace_seeds_report_id_from_selected_evidence_report() {
    let dna = DNAsequence::from_sequence("ACGT").unwrap();
    let mut area = MainAreaDna::new(dna, Some("seq1".to_string()), None);
    area.rna_read_evidence_ui.selected_report_id = "tp73_saved".to_string();
    let view = SplicingExpertView {
        seq_id: "seq1".to_string(),
        target_feature_id: 17,
        scope: SplicingScopePreset::AllOverlappingAnyStrand,
        group_label: "TP73".to_string(),
        strand: "+".to_string(),
        region_start_1based: 1,
        region_end_1based: 4,
        transcript_count: 0,
        unique_exon_count: 0,
        instruction: String::new(),
        transcripts: vec![],
        unique_exons: vec![],
        matrix_rows: vec![],
        boundaries: vec![],
        intron_signals: vec![],
        junctions: vec![],
        events: vec![],
    };

    area.open_rna_read_mapping_workspace_for_view(&view);
    assert!(area.show_rna_read_mapping_window);
    assert!(area.rna_read_mapping_window_pending_initial_render);
    assert!(area.rna_read_mapping_window_focus_requested);
    assert_eq!(area.rna_read_mapping_window_feature_id, Some(17));
    assert_eq!(area.rna_reads_ui.report_id, "tp73_saved");
}

#[test]
fn open_rna_read_mapping_workspace_replaces_stale_report_id_with_matching_report() {
    let dna = DNAsequence::from_sequence("ACGT").unwrap();
    let mut state = ProjectState::default();
    state.sequences.insert("seq1".to_string(), dna.clone());
    let engine = Arc::new(RwLock::new(GentleEngine::from_state(state)));
    engine
        .write()
        .expect("engine")
        .commit_rna_read_report(RnaReadInterpretationReport {
            report_id: "other_locus".to_string(),
            seq_id: "seq1".to_string(),
            seed_feature_id: 9,
            generated_at_unix_ms: 2,
            ..RnaReadInterpretationReport::default()
        })
        .expect("commit stale report");
    engine
        .write()
        .expect("engine")
        .commit_rna_read_report(RnaReadInterpretationReport {
            report_id: "tp73_saved".to_string(),
            seq_id: "seq1".to_string(),
            seed_feature_id: 17,
            generated_at_unix_ms: 1,
            ..RnaReadInterpretationReport::default()
        })
        .expect("commit matching report");
    let mut area = MainAreaDna::new(dna, Some("seq1".to_string()), Some(engine));
    area.rna_read_evidence_ui.selected_report_id = "other_locus".to_string();
    area.rna_reads_ui.report_id = "other_locus".to_string();
    let view = SplicingExpertView {
        seq_id: "seq1".to_string(),
        target_feature_id: 17,
        scope: SplicingScopePreset::AllOverlappingAnyStrand,
        group_label: "TP73".to_string(),
        strand: "+".to_string(),
        region_start_1based: 1,
        region_end_1based: 4,
        transcript_count: 0,
        unique_exon_count: 0,
        instruction: String::new(),
        transcripts: vec![],
        unique_exons: vec![],
        matrix_rows: vec![],
        boundaries: vec![],
        intron_signals: vec![],
        junctions: vec![],
        events: vec![],
    };

    area.open_rna_read_mapping_workspace_for_view(&view);

    assert_eq!(area.rna_reads_ui.report_id, "tp73_saved");
}

#[test]
fn open_current_rna_read_mapping_workspace_uses_selected_splicing_feature() {
    let mut dna =
        DNAsequence::from_sequence("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA").expect("sequence");
    dna.features_mut().push(Feature {
        kind: "mRNA".into(),
        location: Location::Join(vec![
            Location::simple_range(2, 8),
            Location::simple_range(12, 20),
            Location::simple_range(26, 34),
        ]),
        qualifiers: vec![
            ("gene".into(), Some("GENE1".to_string())),
            ("transcript_id".into(), Some("NM_TEST_1".to_string())),
            ("label".into(), Some("NM_TEST_1".to_string())),
        ],
    });
    let mut state = ProjectState::default();
    state.sequences.insert("seq_gene".to_string(), dna.clone());
    let engine = Arc::new(RwLock::new(GentleEngine::from_state(state)));
    let mut area = MainAreaDna::new(dna, Some("seq_gene".to_string()), Some(engine));

    area.focus_feature(0);

    let view = area
        .open_current_rna_read_mapping_workspace()
        .expect("open RNA-read Mapping from current selection");

    assert!(area.show_rna_read_mapping_window);
    assert!(area.rna_read_mapping_window_focus_requested);
    assert_eq!(view.seq_id, "seq_gene");
    assert_eq!(view.target_feature_id, 0);
    assert_eq!(view.group_label, "GENE1");
    assert_eq!(area.rna_read_mapping_window_feature_id, Some(0));
}

#[test]
fn open_current_rna_read_mapping_workspace_falls_back_when_description_cache_is_empty() {
    let mut dna =
        DNAsequence::from_sequence("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA").expect("sequence");
    dna.features_mut().push(Feature {
        kind: "mRNA".into(),
        location: Location::Join(vec![
            Location::simple_range(2, 8),
            Location::simple_range(12, 20),
            Location::simple_range(26, 34),
        ]),
        qualifiers: vec![
            ("gene".into(), Some("GENE1".to_string())),
            ("transcript_id".into(), Some("NM_TEST_1".to_string())),
            ("label".into(), Some("NM_TEST_1".to_string())),
        ],
    });
    let mut state = ProjectState::default();
    state.sequences.insert("seq_gene".to_string(), dna.clone());
    let engine = Arc::new(RwLock::new(GentleEngine::from_state(state)));
    let mut area = MainAreaDna::new(dna, Some("seq_gene".to_string()), Some(engine));

    area.focus_feature(0);
    area.description_cache_expert_view = None;
    area.description_cache_expert_error = None;

    let view = area
        .open_current_rna_read_mapping_workspace()
        .expect("open RNA-read Mapping from current selection without cached expert view");

    assert!(area.show_rna_read_mapping_window);
    assert!(area.rna_read_mapping_window_focus_requested);
    assert_eq!(view.seq_id, "seq_gene");
    assert_eq!(view.target_feature_id, 0);
    assert_eq!(view.group_label, "GENE1");
}

#[test]
fn open_rna_read_mapping_for_feature_opens_workspace_on_explicit_request() {
    let mut dna =
        DNAsequence::from_sequence("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA").expect("sequence");
    dna.features_mut().push(Feature {
        kind: "gene".into(),
        location: Location::simple_range(2, 34),
        qualifiers: vec![
            ("gene".into(), Some("GENE1".to_string())),
            ("label".into(), Some("GENE1".to_string())),
        ],
    });
    dna.features_mut().push(Feature {
        kind: "mRNA".into(),
        location: Location::Join(vec![
            Location::simple_range(2, 8),
            Location::simple_range(12, 20),
            Location::simple_range(26, 34),
        ]),
        qualifiers: vec![
            ("gene".into(), Some("GENE1".to_string())),
            ("transcript_id".into(), Some("NM_TEST_1".to_string())),
            ("label".into(), Some("NM_TEST_1".to_string())),
        ],
    });
    let mut state = ProjectState::default();
    state.sequences.insert("seq_gene".to_string(), dna.clone());
    let engine = Arc::new(RwLock::new(GentleEngine::from_state(state)));
    let mut area = MainAreaDna::new(dna, Some("seq_gene".to_string()), Some(engine));

    let opened = area.open_rna_read_mapping_for_feature(0, "test");

    assert!(opened);
    assert!(area.show_rna_read_mapping_window);
    assert_eq!(area.rna_read_mapping_window_feature_id, Some(0));
    match area.description_cache_expert_view.as_ref() {
        Some(FeatureExpertView::Splicing(view)) => {
            assert_eq!(view.target_feature_id, 0);
            assert_eq!(view.group_label, "GENE1");
        }
        other => panic!("expected splicing expert view for mapping action, got {other:?}"),
    }
}

#[test]
fn open_variant_followup_for_feature_seeds_window_from_variation_feature() {
    let mut dna =
        DNAsequence::from_sequence("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA").expect("sequence");
    dna.features_mut().push(Feature {
        kind: "variation".into(),
        location: Location::simple_range(10, 11),
        qualifiers: vec![
            ("label".into(), Some("rs9923231".to_string())),
            ("db_xref".into(), Some("dbSNP:rs9923231".to_string())),
            ("gene".into(), Some("VKORC1".to_string())),
        ],
    });
    let mut area = MainAreaDna::new(dna, Some("vkorc1_context".to_string()), None);

    let opened = area.open_variant_followup_for_feature(0, "test");

    assert!(opened);
    assert!(area.show_variant_followup_window);
    assert!(area.variant_followup_window_pending_initial_render);
    assert!(area.variant_followup_window_focus_requested);
    assert_eq!(area.variant_followup_ui.source_seq_id, "vkorc1_context");
    assert_eq!(area.variant_followup_ui.source_feature_id, Some(0));
    assert_eq!(area.variant_followup_ui.variant_label_or_id, "rs9923231");
    assert_eq!(area.variant_followup_ui.gene_label, "VKORC1");
    assert_eq!(
        area.variant_followup_ui.fragment_output_id,
        "rs9923231_promoter_fragment"
    );
    assert_eq!(
        area.variant_followup_ui.reference_output_id,
        "rs9923231_promoter_reference"
    );
    assert_eq!(
        area.variant_followup_ui.alternate_output_id,
        "rs9923231_promoter_alternate"
    );
    assert_eq!(area.variant_followup_ui.score_track_motifs, "SP1");
    assert_eq!(
        area.variant_followup_ui.score_track_value_kind,
        TfbsScoreTrackValueKind::LlrBits
    );
    assert!(area.variant_followup_ui.score_track_clip_negative);
    assert_eq!(
        area.variant_followup_ui.tfbs_track_similarity_anchor_motif,
        ""
    );
    assert!(
        area.variant_followup_ui
            .tfbs_track_similarity_all_candidates
    );
    assert_eq!(area.variant_followup_ui.tfbs_track_similarity_limit, "25");
    assert_eq!(
        area.variant_followup_ui
            .tfbs_track_similarity_ranking_metric,
        TfbsTrackSimilarityRankingMetric::SmoothedSpearman
    );
    assert!(
        !area
            .variant_followup_ui
            .tfbs_track_similarity_include_remote_metadata
    );
    assert!(
        area.variant_followup_ui
            .cached_tfbs_track_similarity
            .is_none()
    );
    assert!(
        area.variant_followup_ui
            .cached_alternative_promoter_comparison
            .is_none()
    );
}

#[test]
fn feature_supports_variant_followup_accepts_gene_mrna_promoter_and_variation() {
    let mut dna =
        DNAsequence::from_sequence("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA").expect("sequence");
    dna.features_mut().push(Feature {
        kind: "variation".into(),
        location: Location::simple_range(10, 11),
        qualifiers: vec![("label".into(), Some("rs9923231".to_string()))],
    });
    dna.features_mut().push(Feature {
        kind: "gene".into(),
        location: Location::simple_range(2, 20),
        qualifiers: vec![("label".into(), Some("VKORC1".to_string()))],
    });
    dna.features_mut().push(Feature {
        kind: "mRNA".into(),
        location: Location::simple_range(2, 20),
        qualifiers: vec![("label".into(), Some("NM_VKORC1_1".to_string()))],
    });
    dna.features_mut().push(Feature {
        kind: "promoter".into(),
        location: Location::simple_range(0, 12),
        qualifiers: vec![("label".into(), Some("VKORC1 promoter".to_string()))],
    });
    let area = MainAreaDna::new(dna, Some("seq1".to_string()), None);

    assert!(area.feature_supports_variant_followup(0));
    assert!(area.feature_supports_variant_followup(1));
    assert!(area.feature_supports_variant_followup(2));
    assert!(area.feature_supports_variant_followup(3));
}

#[test]
fn open_variant_followup_for_gene_feature_seeds_promoter_design_presets() {
    let mut dna =
        DNAsequence::from_sequence("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA").expect("sequence");
    dna.features_mut().push(Feature {
        kind: "gene".into(),
        location: Location::simple_range(2, 20),
        qualifiers: vec![("label".into(), Some("TP73".to_string()))],
    });
    let mut area = MainAreaDna::new(dna, Some("tp73_context".to_string()), None);

    let opened = area.open_variant_followup_for_feature(0, "test");

    assert!(opened);
    assert!(area.show_variant_followup_window);
    assert!(area.variant_followup_window_pending_initial_render);
    assert_eq!(area.variant_followup_ui.source_seq_id, "tp73_context");
    assert_eq!(area.variant_followup_ui.source_feature_id, Some(0));
    assert_eq!(area.variant_followup_ui.variant_label_or_id, "");
    assert_eq!(area.variant_followup_ui.gene_label, "TP73");
    assert_eq!(
        area.variant_followup_ui.score_track_motifs,
        "TP73,SP1,BACH2,PATZ1"
    );
    assert_eq!(area.variant_followup_ui.score_track_start_0based, "0");
    assert_eq!(
        area.variant_followup_ui.score_track_end_0based_exclusive,
        "40"
    );
    assert_eq!(
        area.variant_followup_ui.score_track_value_kind,
        TfbsScoreTrackValueKind::LlrBits
    );
    assert!(area.variant_followup_ui.score_track_clip_negative);
}

#[test]
fn open_variant_followup_for_reasoning_promoter_span_resolves_transcript_feature() {
    let mut dna = DNAsequence::from_sequence(&"A".repeat(4000)).expect("sequence");
    dna.features_mut().push(Feature {
        kind: "mRNA".into(),
        location: Location::simple_range(1200, 1800),
        qualifiers: vec![
            ("gene".into(), Some("TP73".to_string())),
            ("transcript_id".into(), Some("ENSTTP73A".to_string())),
            ("label".into(), Some("TP73-201".to_string())),
        ],
    });
    let mut state = ProjectState::default();
    state
        .sequences
        .insert("tp73_reasoning_ui".to_string(), dna.clone());
    let mut engine = GentleEngine::from_state(state);
    let graph = engine
        .build_construct_reasoning_graph("tp73_reasoning_ui", None, Some("graph_tp73_ui"))
        .expect("graph");
    let promoter_evidence_id = graph
        .evidence
        .iter()
        .find(|row| row.role == ConstructRole::Promoter)
        .map(|row| row.evidence_id.clone())
        .expect("promoter evidence");
    let engine = Arc::new(RwLock::new(engine));
    let mut area = MainAreaDna::new(dna, Some("tp73_reasoning_ui".to_string()), Some(engine));
    area.focus_construct_reasoning_graph("graph_tp73_ui");

    let opened =
        area.open_variant_followup_for_reasoning_span(promoter_evidence_id.as_str(), "test");

    assert!(opened);
    assert!(area.show_variant_followup_window);
    assert_eq!(area.variant_followup_ui.source_feature_id, Some(0));
    assert_eq!(area.variant_followup_ui.gene_label, "TP73");
    assert_eq!(area.variant_followup_ui.transcript_id, "ENSTTP73A");
}

#[test]
fn use_variant_followup_alternative_promoter_row_retargets_promoter_design_span() {
    let dna = DNAsequence::from_sequence(&"A".repeat(4000)).expect("sequence");
    let mut area = MainAreaDna::new(dna, Some("tp73_context".to_string()), None);
    area.variant_followup_ui.cached_score_tracks = Some(TfbsScoreTrackReport::default());
    area.variant_followup_ui.cached_tfbs_track_similarity =
        Some(TfbsTrackSimilarityReport::default());
    area.variant_followup_ui.cached_report = Some(VariantPromoterContextReport::default());
    area.variant_followup_ui.cached_candidates = Some(PromoterReporterCandidateSet::default());

    let row = AlternativePromoterComparisonRow {
        label: "TP73 promoter window (2 tx)".to_string(),
        gene_label: Some("TP73".to_string()),
        gene_id: None,
        representative_transcript_id: Some("ENSTTP73A".to_string()),
        representative_transcript_label: Some("TP73-201".to_string()),
        representative_transcript_feature_id: Some(7),
        transcript_count: 2,
        transcript_ids: vec!["ENSTTP73A".to_string(), "ENSTTP73B".to_string()],
        transcript_labels: vec!["TP73-201".to_string(), "TP73-202".to_string()],
        strand: "+".to_string(),
        representative_tss_local_0based: 1200,
        start_0based: 200,
        end_0based_exclusive: 1400,
        upstream_bp: 1000,
        downstream_bp: 200,
        source: "transcript_tss".to_string(),
    };

    area.use_variant_followup_alternative_promoter_row(&row);

    assert_eq!(area.variant_followup_ui.transcript_id, "ENSTTP73A");
    assert_eq!(area.variant_followup_ui.gene_label, "TP73");
    assert_eq!(area.variant_followup_ui.promoter_upstream_bp, "1000");
    assert_eq!(area.variant_followup_ui.promoter_downstream_bp, "200");
    assert_eq!(area.variant_followup_ui.score_track_start_0based, "200");
    assert_eq!(
        area.variant_followup_ui.score_track_end_0based_exclusive,
        "1400"
    );
    assert!(area.variant_followup_ui.cached_score_tracks.is_none());
    assert!(
        area.variant_followup_ui
            .cached_tfbs_track_similarity
            .is_none()
    );
    assert!(area.variant_followup_ui.cached_report.is_none());
    assert!(area.variant_followup_ui.cached_candidates.is_none());
    assert!(
        area.op_status.contains("ENSTTP73A"),
        "status was: {}",
        area.op_status
    );
}

#[test]
fn use_variant_followup_isoform_promoter_group_retargets_promoter_design_span() {
    let dna = DNAsequence::from_sequence(&"A".repeat(4000)).expect("sequence");
    let mut area = MainAreaDna::new(dna, Some("tp73_context".to_string()), None);
    area.variant_followup_ui.cached_score_tracks = Some(TfbsScoreTrackReport::default());
    area.variant_followup_ui.cached_tfbs_track_similarity =
        Some(TfbsTrackSimilarityReport::default());
    area.variant_followup_ui.cached_report = Some(VariantPromoterContextReport::default());
    area.variant_followup_ui.cached_candidates = Some(PromoterReporterCandidateSet::default());
    area.variant_followup_ui.cached_promoter_evidence_matrix =
        Some(crate::engine::PromoterEvidenceMatrixReport::default());
    area.variant_followup_ui.cached_isoform_promoter_comparison =
        Some(crate::engine::IsoformPromoterComparisonReport::default());

    let group = crate::engine::IsoformPromoterComparisonGroup {
        group_id: "promoter_group_1".to_string(),
        label: "TP73 promoter group 1".to_string(),
        gene_label: Some("TP73".to_string()),
        strand: "+".to_string(),
        start_0based: 500,
        end_0based_exclusive: 1700,
        transcript_count: 2,
        transcript_ids: vec!["ENSTTP73A".to_string(), "ENSTTP73B".to_string()],
        transcript_labels: vec!["TP73-201".to_string(), "TP73-202".to_string()],
        ..crate::engine::IsoformPromoterComparisonGroup::default()
    };

    area.use_variant_followup_isoform_promoter_group(&group, 1000, 200);

    assert_eq!(area.variant_followup_ui.transcript_id, "ENSTTP73A");
    assert_eq!(area.variant_followup_ui.gene_label, "TP73");
    assert_eq!(area.variant_followup_ui.promoter_upstream_bp, "1000");
    assert_eq!(area.variant_followup_ui.promoter_downstream_bp, "200");
    assert_eq!(area.variant_followup_ui.score_track_start_0based, "500");
    assert_eq!(
        area.variant_followup_ui.score_track_end_0based_exclusive,
        "1700"
    );
    assert!(area.variant_followup_ui.cached_score_tracks.is_none());
    assert!(
        area.variant_followup_ui
            .cached_tfbs_track_similarity
            .is_none()
    );
    assert!(area.variant_followup_ui.cached_report.is_none());
    assert!(area.variant_followup_ui.cached_candidates.is_none());
    assert!(
        area.variant_followup_ui
            .cached_promoter_evidence_matrix
            .is_none()
    );
    assert!(
        area.variant_followup_ui
            .cached_isoform_promoter_comparison
            .is_some()
    );
    assert!(
        area.op_status.contains("TP73 promoter group 1"),
        "status was: {}",
        area.op_status
    );
}

#[test]
fn variant_followup_promoter_expression_evidence_runs_shared_op_and_caches_report() {
    let mut dna = DNAsequence::from_sequence(&"A".repeat(2500)).expect("sequence");
    for (transcript_id, transcript_label, start) in [
        ("ENSTTP73A", "TP73-201", 1200_i64),
        ("ENSTTP73B", "TP73-202", 1200_i64),
    ] {
        dna.features_mut().push(Feature {
            kind: "mRNA".into(),
            location: Location::simple_range(start, start + 500),
            qualifiers: vec![
                ("gene".into(), Some("TP73".to_string())),
                ("transcript_id".into(), Some(transcript_id.to_string())),
                ("label".into(), Some(transcript_label.to_string())),
            ],
        });
    }

    let mut state = ProjectState::default();
    state
        .sequences
        .insert("tp73_expression_gui".to_string(), dna.clone());
    let engine = Arc::new(RwLock::new(GentleEngine::from_state(state)));
    let mut area = MainAreaDna::new(dna, Some("tp73_expression_gui".to_string()), Some(engine));
    area.variant_followup_ui.source_seq_id = "tp73_expression_gui".to_string();
    area.variant_followup_ui.gene_label = "TP73".to_string();
    area.variant_followup_ui.promoter_upstream_bp = "200".to_string();
    area.variant_followup_ui.promoter_downstream_bp = "50".to_string();
    area.variant_followup_ui.promoter_expression_source_label =
        "synthetic expression table".to_string();
    area.variant_followup_ui.promoter_expression_rows_json = serde_json::to_string(&vec![
        PromoterExpressionEvidenceInput {
            transcript_id: Some("ENSTTP73A".to_string()),
            sample_id: Some("case_1".to_string()),
            condition: Some("case".to_string()),
            value: 10.0,
            unit: Some("TPM".to_string()),
            source: Some("synthetic RNA-seq".to_string()),
            ..PromoterExpressionEvidenceInput::default()
        },
        PromoterExpressionEvidenceInput {
            transcript_id: Some("ENSTUNRELATED".to_string()),
            value: 99.0,
            unit: Some("TPM".to_string()),
            ..PromoterExpressionEvidenceInput::default()
        },
    ])
    .expect("serialize expression rows");

    area.summarize_variant_followup_promoter_expression_evidence();

    let report = area
        .variant_followup_ui
        .cached_promoter_expression_evidence
        .as_ref()
        .expect("cached promoter expression evidence report");
    assert_eq!(report.schema, "gentle.promoter_expression_evidence.v1");
    assert_eq!(report.seq_id, "tp73_expression_gui");
    assert_eq!(report.promoter_group_count, 1);
    assert_eq!(report.supplied_expression_record_count, 2);
    assert_eq!(report.assigned_expression_record_count, 1);
    assert_eq!(report.unassigned_expression_records.len(), 1);
    assert_eq!(report.expression_source_label, "synthetic expression table");
    assert!(
        report
            .warnings
            .iter()
            .any(|warning| warning.contains("could not be assigned")),
        "warnings were: {:?}",
        report.warnings
    );
    let row = report.rows.first().expect("promoter expression row");
    assert_eq!(row.expression_record_count, 1);
    assert_eq!(row.mean_value, Some(10.0));
    assert_eq!(row.max_value, Some(10.0));
    assert_eq!(row.unit.as_deref(), Some("TPM"));
    assert_eq!(row.records[0].matched_by, vec!["transcript_id"]);
    assert_eq!(row.records[0].source, "synthetic RNA-seq");
}

#[test]
fn splicing_intron_regulatory_rows_merge_cached_attract_hits_with_intron_signals() {
    let view = SplicingExpertView {
        seq_id: "tp73".to_string(),
        target_feature_id: 1,
        scope: SplicingScopePreset::TargetGroupTargetStrand,
        group_label: "TP73".to_string(),
        strand: "+".to_string(),
        region_start_1based: 1,
        region_end_1based: 30,
        transcript_count: 1,
        unique_exon_count: 2,
        instruction: String::new(),
        transcripts: vec![SplicingTranscriptLane {
            transcript_feature_id: 1,
            transcript_id: "ENSTTP73A".to_string(),
            label: "TP73-201".to_string(),
            strand: "+".to_string(),
            exons: vec![
                SplicingRange {
                    start_1based: 1,
                    end_1based: 10,
                },
                SplicingRange {
                    start_1based: 21,
                    end_1based: 30,
                },
            ],
            exon_cds_phases: vec![],
            introns: vec![SplicingRange {
                start_1based: 11,
                end_1based: 20,
            }],
            has_target_feature: true,
        }],
        unique_exons: vec![
            SplicingExonSummary {
                start_1based: 1,
                end_1based: 10,
                support_transcript_count: 1,
                constitutive: true,
            },
            SplicingExonSummary {
                start_1based: 21,
                end_1based: 30,
                support_transcript_count: 1,
                constitutive: true,
            },
        ],
        matrix_rows: vec![],
        boundaries: vec![],
        intron_signals: vec![SplicingIntronSignal {
            transcript_feature_id: 1,
            transcript_id: "ENSTTP73A".to_string(),
            donor_position_1based: 10,
            acceptor_position_1based: 21,
            intron_length_bp: 10,
            branchpoint_position_1based: Some(17),
            branchpoint_motif: "TACTAAC".to_string(),
            branchpoint_score: 9.5,
            branchpoint_annotation: "heuristic".to_string(),
            polypyrimidine_start_1based: Some(18),
            polypyrimidine_end_1based: Some(20),
            polypyrimidine_fraction: 0.91,
            polypyrimidine_annotation: "tract".to_string(),
        }],
        junctions: vec![],
        events: vec![],
    };
    let evidence = AttractSplicingEvidenceView {
        schema: "gentle.attract_splicing_evidence.v1".to_string(),
        seq_id: "tp73".to_string(),
        target_feature_id: 1,
        scope: SplicingScopePreset::TargetGroupTargetStrand,
        group_label: "TP73".to_string(),
        target_strand: "+".to_string(),
        settings: AttractSplicingEvidenceSettings::default(),
        requested_organism: Some("Homo sapiens".to_string()),
        resolved_organism: Some("Homo sapiens".to_string()),
        species_match_mode: AttractSpeciesMatchMode::ExactOrganism,
        scanned_transcript_count: 1,
        scanned_window_count: 3,
        unique_rbp_count: 3,
        hit_count: 4,
        pwm_scored_hit_count: 0,
        exact_length_pwm_hit_count: 0,
        windowed_pwm_hit_count: 0,
        consensus_hit_count: 4,
        active_resource_source: "synthetic".to_string(),
        active_resource_item_count: 4,
        active_resource_pwm_row_count: 0,
        active_resource_consensus_only_row_count: 4,
        active_resource_fingerprint: None,
        alternate_policy_summary: None,
        summary_rows: vec![],
        hit_rows: vec![
            AttractSplicingEvidenceHitRow {
                transcript_feature_id: 1,
                transcript_id: "ENSTTP73A".to_string(),
                transcript_label: "TP73-201".to_string(),
                transcript_strand: "+".to_string(),
                gene_name: "SRSF1".to_string(),
                organism: "Homo sapiens".to_string(),
                matrix_id: "M001".to_string(),
                motif_iupac: "GAAGAA".to_string(),
                model_kind: "consensus_iupac".to_string(),
                region_class: AttractRegionClass::DonorFlank,
                region_start_1based: 9,
                region_end_1based: 11,
                region_local_start_1based: 1,
                region_local_end_1based: 3,
                matched_sequence: "GAA".to_string(),
                match_score: 1.0,
                match_score_kind: "llr_bits".to_string(),
                match_score_quantile: None,
                quality_score: 3.0,
                exact_species_match: true,
                pwm_mapping_status: "none".to_string(),
                mapping_policy_used: "strict".to_string(),
                pfm_subwindow_start_1based: None,
                pfm_subwindow_end_1based: None,
                warnings: vec![],
            },
            AttractSplicingEvidenceHitRow {
                transcript_feature_id: 1,
                transcript_id: "ENSTTP73A".to_string(),
                transcript_label: "TP73-201".to_string(),
                transcript_strand: "+".to_string(),
                gene_name: "PTBP1".to_string(),
                organism: "Homo sapiens".to_string(),
                matrix_id: "M002".to_string(),
                motif_iupac: "TCTT".to_string(),
                model_kind: "consensus_iupac".to_string(),
                region_class: AttractRegionClass::AcceptorFlank,
                region_start_1based: 20,
                region_end_1based: 22,
                region_local_start_1based: 1,
                region_local_end_1based: 3,
                matched_sequence: "TCT".to_string(),
                match_score: 1.0,
                match_score_kind: "llr_bits".to_string(),
                match_score_quantile: None,
                quality_score: 3.0,
                exact_species_match: true,
                pwm_mapping_status: "none".to_string(),
                mapping_policy_used: "strict".to_string(),
                pfm_subwindow_start_1based: None,
                pfm_subwindow_end_1based: None,
                warnings: vec![],
            },
            AttractSplicingEvidenceHitRow {
                transcript_feature_id: 1,
                transcript_id: "ENSTTP73A".to_string(),
                transcript_label: "TP73-201".to_string(),
                transcript_strand: "+".to_string(),
                gene_name: "RBM5".to_string(),
                organism: "Homo sapiens".to_string(),
                matrix_id: "M003".to_string(),
                motif_iupac: "UGCAUG".to_string(),
                model_kind: "consensus_iupac".to_string(),
                region_class: AttractRegionClass::IntronBody,
                region_start_1based: 12,
                region_end_1based: 18,
                region_local_start_1based: 1,
                region_local_end_1based: 7,
                matched_sequence: "TGCATG".to_string(),
                match_score: 1.0,
                match_score_kind: "llr_bits".to_string(),
                match_score_quantile: None,
                quality_score: 3.0,
                exact_species_match: true,
                pwm_mapping_status: "none".to_string(),
                mapping_policy_used: "strict".to_string(),
                pfm_subwindow_start_1based: None,
                pfm_subwindow_end_1based: None,
                warnings: vec![],
            },
            AttractSplicingEvidenceHitRow {
                transcript_feature_id: 99,
                transcript_id: "OTHER".to_string(),
                transcript_label: "OTHER".to_string(),
                transcript_strand: "+".to_string(),
                gene_name: "IGNORE".to_string(),
                organism: "Homo sapiens".to_string(),
                matrix_id: "M004".to_string(),
                motif_iupac: "AAAA".to_string(),
                model_kind: "consensus_iupac".to_string(),
                region_class: AttractRegionClass::IntronBody,
                region_start_1based: 12,
                region_end_1based: 18,
                region_local_start_1based: 1,
                region_local_end_1based: 7,
                matched_sequence: "AAAA".to_string(),
                match_score: 1.0,
                match_score_kind: "llr_bits".to_string(),
                match_score_quantile: None,
                quality_score: 3.0,
                exact_species_match: true,
                pwm_mapping_status: "none".to_string(),
                mapping_policy_used: "strict".to_string(),
                pfm_subwindow_start_1based: None,
                pfm_subwindow_end_1based: None,
                warnings: vec![],
            },
        ],
        warnings: vec![],
    };

    let rows = MainAreaDna::splicing_intron_regulatory_rows(&view, Some(&evidence));

    assert_eq!(rows.len(), 1);
    assert_eq!(rows[0].donor_flank_hit_count, 1);
    assert_eq!(rows[0].acceptor_flank_hit_count, 1);
    assert_eq!(rows[0].intron_body_hit_count, 1);
    assert_eq!(rows[0].unique_rbp_count, 3);
    assert_eq!(rows[0].top_rbps, vec!["SRSF1", "PTBP1", "RBM5"]);
    let summary = MainAreaDna::splicing_intron_regulatory_summary_text(&rows[0]);
    assert!(summary.contains("branchpoint-like @ 17"));
    assert!(summary.contains("polyY tract 18..20"));
    assert!(summary.contains("3 unique RBP(s): SRSF1, PTBP1, RBM5"));
}

#[test]
fn variant_followup_handoff_result_uses_bundle_relative_artifacts() {
    let dna =
        DNAsequence::from_sequence("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA").expect("sequence");
    let mut area = MainAreaDna::new(dna, Some("vkorc1_context".to_string()), None);
    area.variant_followup_ui.source_seq_id = "vkorc1_context".to_string();
    area.variant_followup_ui.variant_label_or_id = "rs9923231".to_string();
    area.variant_followup_ui.gene_label = "VKORC1".to_string();
    area.variant_followup_ui.fragment_output_id = "vkorc1_fragment".to_string();
    area.variant_followup_ui.reference_output_id = "vkorc1_reference".to_string();
    area.variant_followup_ui.alternate_output_id = "vkorc1_alternate".to_string();
    area.variant_followup_ui.reporter_backbone_seq_id =
        "gentle_mammalian_luciferase_backbone_v1".to_string();
    area.variant_followup_ui.reporter_backbone_path =
        "data/tutorial_inputs/gentle_mammalian_luciferase_backbone_v1.gb".to_string();
    let artifacts = area.variant_followup_bundle_artifacts(
        "vkorc1_reporter_reference",
        "vkorc1_reporter_alternate",
    );
    let report = VariantPromoterContextReport {
        schema: "gentle.variant_promoter_context.v1".to_string(),
        seq_id: "vkorc1_context".to_string(),
        sequence_length_bp: 6001,
        generated_at_unix_ms: 0,
        op_id: None,
        run_id: None,
        variant_label: "rs9923231".to_string(),
        variant_feature_id: Some(0),
        variant_start_0based: 3000,
        variant_end_0based_exclusive: 3001,
        variant_class: Some("snv".to_string()),
        genomic_ref: Some("C".to_string()),
        genomic_alt: Some("T".to_string()),
        genome_anchor: Some(SequenceGenomeAnchorSummary {
            seq_id: "vkorc1_context".to_string(),
            genome_id: "Human GRCh38 Ensembl 116".to_string(),
            chromosome: "16".to_string(),
            start_1based: 31093368,
            end_1based: 31099368,
            strand: Some('+'),
            anchor_verified: Some(true),
        }),
        chosen_gene_label: Some("VKORC1".to_string()),
        chosen_transcript_id: Some("ENST_TEST".to_string()),
        chosen_transcript_label: Some("VKORC1-201".to_string()),
        transcript_ambiguity_status: "resolved".to_string(),
        promoter_upstream_bp: 1000,
        promoter_downstream_bp: 200,
        promoter_overlap: true,
        signed_tss_distance_bp: Some(-88),
        overlapping_gene_labels: vec!["VKORC1".to_string()],
        overlapping_transcript_labels: vec!["VKORC1-201".to_string()],
        overlapping_promoter_labels: vec!["VKORC1 promoter".to_string()],
        overlapping_tfbs_labels: vec!["SP1".to_string()],
        overlapping_evidence: vec![],
        promoter_windows_considered: vec![],
        effect_tags: vec!["promoter_variant_candidate".to_string()],
        suggested_assay_ids: vec!["allele_paired_promoter_luciferase_reporter".to_string()],
        tfbs_focus_half_window_bp: 100,
        tfbs_near_variant_status: "annotations_present".to_string(),
        tfbs_region_summary: None,
        rationale: "Variant falls into the derived promoter window.".to_string(),
    };
    let candidates = PromoterReporterCandidateSet {
        schema: "gentle.promoter_reporter_candidates.v1".to_string(),
        seq_id: "vkorc1_context".to_string(),
        sequence_length_bp: 6001,
        generated_at_unix_ms: 0,
        op_id: None,
        run_id: None,
        variant_label: "rs9923231".to_string(),
        chosen_gene_label: Some("VKORC1".to_string()),
        chosen_transcript_id: Some("ENST_TEST".to_string()),
        chosen_transcript_label: Some("VKORC1-201".to_string()),
        transcript_ambiguity_status: "resolved".to_string(),
        retain_downstream_from_tss_bp: 200,
        retain_upstream_beyond_variant_bp: 500,
        max_candidates: 5,
        recommended_candidate_id: "cand-1".to_string(),
        suggested_assay_ids: vec!["allele_paired_promoter_luciferase_reporter".to_string()],
        candidates: vec![],
    };
    let recommended = crate::engine::PromoterReporterFragmentCandidate {
        candidate_id: "cand-1".to_string(),
        gene_label: Some("VKORC1".to_string()),
        transcript_id: "ENST_TEST".to_string(),
        transcript_label: "VKORC1-201".to_string(),
        strand: "-".to_string(),
        tss_local_0based: 2412,
        variant_start_0based: 3000,
        variant_end_0based_exclusive: 3001,
        start_0based: 2412,
        end_0based_exclusive: 3501,
        length_bp: 1089,
        retain_downstream_from_tss_bp: 200,
        retain_upstream_beyond_variant_bp: 500,
        promoter_overlap: true,
        signed_tss_distance_bp: -88,
        rank: 1,
        recommended: true,
        rationale: "Top-ranked reverse-strand promoter fragment".to_string(),
    };

    let result = area.build_variant_followup_handoff_result_json(
        &artifacts,
        &report,
        &candidates,
        &recommended,
        "vkorc1_reporter_reference",
        "vkorc1_reporter_alternate",
    );

    assert_eq!(result["schema"].as_str(), Some("gentle.clawbio_handoff.v2"));
    assert_eq!(
        result["artifacts"]["promoter_context_svg"].as_str(),
        Some(artifacts.promoter_context_svg.as_str())
    );
    assert_eq!(
        result["design"]["construct_previews"]["reference"]["svg_path"].as_str(),
        Some(artifacts.reference_reporter_svg.as_str())
    );
    assert_eq!(
        result["design"]["construct_previews"]["alternate"]["svg_path"].as_str(),
        Some(artifacts.alternate_reporter_svg.as_str())
    );
}

#[test]
fn variant_followup_handoff_commands_render_svg_paths_expand_bundle_dir() {
    let dna =
        DNAsequence::from_sequence("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA").expect("sequence");
    let mut area = MainAreaDna::new(dna, Some("vkorc1_context".to_string()), None);
    area.variant_followup_ui.source_seq_id = "vkorc1_context".to_string();
    area.variant_followup_ui.variant_label_or_id = "rs9923231".to_string();
    area.variant_followup_ui.gene_label = "VKORC1".to_string();
    area.variant_followup_ui.fragment_output_id = "vkorc1_fragment".to_string();
    area.variant_followup_ui.reference_output_id = "vkorc1_reference".to_string();
    area.variant_followup_ui.alternate_output_id = "vkorc1_alternate".to_string();
    area.variant_followup_ui.reporter_backbone_seq_id =
        "gentle_mammalian_luciferase_backbone_v1".to_string();
    area.variant_followup_ui.reporter_backbone_path =
        "data/tutorial_inputs/gentle_mammalian_luciferase_backbone_v1.gb".to_string();
    let artifacts =
        area.variant_followup_bundle_artifacts("rs9923231_reference", "rs9923231_alternate");
    let report = VariantPromoterContextReport {
        schema: "gentle.variant_promoter_context.v1".to_string(),
        seq_id: "vkorc1_context".to_string(),
        sequence_length_bp: 6001,
        generated_at_unix_ms: 0,
        op_id: None,
        run_id: None,
        variant_label: "rs9923231".to_string(),
        variant_feature_id: Some(0),
        variant_start_0based: 3000,
        variant_end_0based_exclusive: 3001,
        variant_class: Some("snv".to_string()),
        genomic_ref: Some("C".to_string()),
        genomic_alt: Some("T".to_string()),
        genome_anchor: Some(SequenceGenomeAnchorSummary {
            seq_id: "vkorc1_context".to_string(),
            genome_id: "Human GRCh38 Ensembl 116".to_string(),
            chromosome: "16".to_string(),
            start_1based: 31093368,
            end_1based: 31099368,
            strand: Some('+'),
            anchor_verified: Some(true),
        }),
        chosen_gene_label: Some("VKORC1".to_string()),
        chosen_transcript_id: Some("ENST_TEST".to_string()),
        chosen_transcript_label: Some("VKORC1-201".to_string()),
        transcript_ambiguity_status: "resolved".to_string(),
        promoter_upstream_bp: 1000,
        promoter_downstream_bp: 200,
        promoter_overlap: true,
        signed_tss_distance_bp: Some(-88),
        overlapping_gene_labels: vec!["VKORC1".to_string()],
        overlapping_transcript_labels: vec!["VKORC1-201".to_string()],
        overlapping_promoter_labels: vec!["VKORC1 promoter".to_string()],
        overlapping_tfbs_labels: vec!["SP1".to_string()],
        overlapping_evidence: vec![],
        promoter_windows_considered: vec![],
        effect_tags: vec!["promoter_variant_candidate".to_string()],
        suggested_assay_ids: vec!["allele_paired_promoter_luciferase_reporter".to_string()],
        tfbs_focus_half_window_bp: 100,
        tfbs_near_variant_status: "annotations_present".to_string(),
        tfbs_region_summary: None,
        rationale: "Variant falls into the derived promoter window.".to_string(),
    };
    let recommended = crate::engine::PromoterReporterFragmentCandidate {
        candidate_id: "cand-1".to_string(),
        gene_label: Some("VKORC1".to_string()),
        transcript_id: "ENST_TEST".to_string(),
        transcript_label: "VKORC1-201".to_string(),
        strand: "-".to_string(),
        tss_local_0based: 2412,
        variant_start_0based: 3000,
        variant_end_0based_exclusive: 3001,
        start_0based: 2412,
        end_0based_exclusive: 3501,
        length_bp: 1089,
        retain_downstream_from_tss_bp: 200,
        retain_upstream_beyond_variant_bp: 500,
        promoter_overlap: true,
        signed_tss_distance_bp: -88,
        rank: 1,
        recommended: true,
        rationale: "Top-ranked reverse-strand promoter fragment".to_string(),
    };

    let commands = area.build_variant_followup_handoff_commands(&artifacts, &report, &recommended);

    assert!(commands.contains("render-svg vkorc1_context linear \"$BUNDLE_DIR/"));
    assert!(commands.contains("render-svg rs9923231_reference circular \"$BUNDLE_DIR/"));
    assert!(commands.contains("render-svg rs9923231_alternate circular \"$BUNDLE_DIR/"));
    assert!(!commands.contains(
        "RenderSequenceSvg\":{\"seq_id\":\"vkorc1_context\",\"mode\":\"Linear\",\"path\":\"$BUNDLE_DIR/"
    ));
}

#[test]
fn open_dotplot_for_feature_opens_transcript_dotplot_on_explicit_request() {
    let mut dna =
        DNAsequence::from_sequence("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA").expect("sequence");
    dna.features_mut().push(Feature {
        kind: "gene".into(),
        location: Location::simple_range(2, 34),
        qualifiers: vec![
            ("gene".into(), Some("GENE1".to_string())),
            ("label".into(), Some("GENE1".to_string())),
        ],
    });
    dna.features_mut().push(Feature {
        kind: "mRNA".into(),
        location: Location::Join(vec![
            Location::simple_range(2, 8),
            Location::simple_range(12, 20),
            Location::simple_range(26, 34),
        ]),
        qualifiers: vec![
            ("gene".into(), Some("GENE1".to_string())),
            ("transcript_id".into(), Some("NM_TEST_1".to_string())),
            ("label".into(), Some("NM_TEST_1".to_string())),
        ],
    });
    let mut state = ProjectState::default();
    state.sequences.insert("seq_gene".to_string(), dna.clone());
    let engine = Arc::new(RwLock::new(GentleEngine::from_state(state)));
    let mut area = MainAreaDna::new(dna, Some("seq_gene".to_string()), Some(engine));

    let opened = area.open_dotplot_for_feature(0, "test");

    assert!(opened);
    assert!(area.show_dotplot_window);
    assert_eq!(area.primary_map_mode, PrimaryMapMode::Dotplot);
    assert_eq!(area.dotplot_ui.reference_seq_id, "seq_gene");
    match area.description_cache_expert_view.as_ref() {
        Some(FeatureExpertView::Splicing(view)) => {
            assert_eq!(view.target_feature_id, 0);
            assert_eq!(view.group_label, "GENE1");
        }
        other => panic!("expected splicing expert view for dotplot action, got {other:?}"),
    }
}

#[test]
fn show_splicing_expert_for_rna_read_mapping_view_selects_current_mapping_report() {
    let dna = DNAsequence::from_sequence("ACGT").unwrap();
    let mut area = MainAreaDna::new(dna, Some("seq1".to_string()), None);
    area.rna_reads_ui.report_id = "tp73_run".to_string();
    area.rna_read_mapping_window_view = Some(Arc::new(SplicingExpertView {
        seq_id: "seq1".to_string(),
        target_feature_id: 17,
        scope: SplicingScopePreset::AllOverlappingAnyStrand,
        group_label: "TP73".to_string(),
        strand: "+".to_string(),
        region_start_1based: 1,
        region_end_1based: 4,
        transcript_count: 0,
        unique_exon_count: 0,
        instruction: String::new(),
        transcripts: vec![],
        unique_exons: vec![],
        matrix_rows: vec![],
        boundaries: vec![],
        intron_signals: vec![],
        junctions: vec![],
        events: vec![],
    }));

    area.show_splicing_expert_for_rna_read_mapping_view();
    assert!(area.show_splicing_expert_window);
    assert_eq!(area.rna_read_evidence_ui.selected_report_id, "tp73_run");
}

#[test]
fn sanitize_workflow_run_id_component_normalizes_symbols_and_case() {
    assert_eq!(
        MainAreaDna::sanitize_workflow_run_id_component("TP73 (tp73.ncbi)/RNA"),
        "tp73_tp73_ncbi_rna"
    );
    assert_eq!(MainAreaDna::sanitize_workflow_run_id_component("___"), "");
}

#[test]
fn default_rna_read_workflow_run_id_uses_sanitized_components() {
    assert_eq!(
        MainAreaDna::default_rna_read_workflow_run_id("TP73.ncbi", "cdna SRR32957124"),
        "workflow_rna_reads_tp73_ncbi_cdna_srr32957124"
    );
    assert_eq!(
        MainAreaDna::default_rna_read_workflow_run_id("", ""),
        "workflow_rna_reads_seq_cdna"
    );
}

#[test]
fn default_rna_read_report_id_tracks_scope_origin_and_target_genes() {
    let view = SplicingExpertView {
        seq_id: "TP73.ncbi".to_string(),
        target_feature_id: 42,
        scope: SplicingScopePreset::TargetGroupTargetStrand,
        group_label: "TP73".to_string(),
        strand: "+".to_string(),
        region_start_1based: 1,
        region_end_1based: 8,
        transcript_count: 0,
        unique_exon_count: 0,
        instruction: String::new(),
        transcripts: vec![],
        unique_exons: vec![],
        matrix_rows: vec![],
        boundaries: vec![],
        intron_signals: vec![],
        junctions: vec![],
        events: vec![],
    };
    let mut ui_state = RnaReadInterpretOpsUiState::default();
    ui_state.input_path = "/tmp/SRR32957124.fasta.gz".to_string();
    ui_state.scope = SplicingScopePreset::TargetGroupTargetStrand;
    ui_state.origin_mode = RnaReadOriginMode::MultiGeneSparse;
    ui_state.target_gene_ids = "TP53 TP63".to_string();

    assert_eq!(
        MainAreaDna::default_rna_read_report_id(&view, &ui_state),
        "rna_tp73_srr32957124_ncdna_tgts_mgs_tg_tp53_tp63"
    );
}

#[test]
fn build_splicing_rna_read_workflow_uses_interpret_payload_and_autorunid() {
    let dna = DNAsequence::from_sequence("ACGTACGT").expect("sequence");
    let mut area = MainAreaDna::new(dna, Some("TP73.ncbi".to_string()), None);
    area.rna_reads_ui.input_path = "/tmp/SRR32957124.fasta.gz".to_string();
    area.rna_reads_ui.seed_stride_bp = "3".to_string();
    area.rna_reads_ui.report_id.clear();
    area.workflow_run_id.clear();
    let view = SplicingExpertView {
        seq_id: "TP73.ncbi".to_string(),
        target_feature_id: 42,
        scope: SplicingScopePreset::AllOverlappingAnyStrand,
        group_label: "TP73".to_string(),
        strand: "+".to_string(),
        region_start_1based: 1,
        region_end_1based: 8,
        transcript_count: 0,
        unique_exon_count: 0,
        instruction: String::new(),
        transcripts: vec![],
        unique_exons: vec![],
        matrix_rows: vec![],
        boundaries: vec![],
        intron_signals: vec![],
        junctions: vec![],
        events: vec![],
    };

    let workflow = area
        .build_splicing_rna_read_workflow(&view)
        .expect("workflow");
    assert_eq!(
        workflow.run_id,
        "workflow_rna_reads_tp73_ncbi_rna_tp73_srr32957124_ncdna_aoas_sg"
    );
    assert_eq!(workflow.ops.len(), 1);
    match &workflow.ops[0] {
        Operation::InterpretRnaReads {
            seq_id,
            seed_feature_id,
            input_path,
            seed_filter,
            report_id,
            ..
        } => {
            assert_eq!(seq_id, "TP73.ncbi");
            assert_eq!(*seed_feature_id, 42);
            assert_eq!(input_path, "/tmp/SRR32957124.fasta.gz");
            assert_eq!(seed_filter.seed_stride_bp, 3);
            assert_eq!(
                report_id.as_deref(),
                Some("rna_tp73_srr32957124_ncdna_aoas_sg")
            );
        }
        other => panic!("unexpected workflow op: {other:?}"),
    }
    assert_eq!(
        area.rna_reads_ui.report_id,
        "rna_tp73_srr32957124_ncdna_aoas_sg"
    );
}

#[test]
fn build_splicing_rna_read_workflow_keeps_manual_report_id_override() {
    let dna = DNAsequence::from_sequence("ACGTACGT").expect("sequence");
    let mut area = MainAreaDna::new(dna, Some("TP73.ncbi".to_string()), None);
    area.rna_reads_ui.input_path = "/tmp/SRR32957124.fasta.gz".to_string();
    area.rna_reads_ui.report_id = "tp73_manual_review".to_string();
    area.rna_reads_ui.report_id_auto_sync = false;
    area.workflow_run_id.clear();
    let view = SplicingExpertView {
        seq_id: "TP73.ncbi".to_string(),
        target_feature_id: 42,
        scope: SplicingScopePreset::AllOverlappingAnyStrand,
        group_label: "TP73".to_string(),
        strand: "+".to_string(),
        region_start_1based: 1,
        region_end_1based: 8,
        transcript_count: 0,
        unique_exon_count: 0,
        instruction: String::new(),
        transcripts: vec![],
        unique_exons: vec![],
        matrix_rows: vec![],
        boundaries: vec![],
        intron_signals: vec![],
        junctions: vec![],
        events: vec![],
    };

    let workflow = area
        .build_splicing_rna_read_workflow(&view)
        .expect("workflow");
    assert_eq!(
        workflow.run_id,
        "workflow_rna_reads_tp73_ncbi_tp73_manual_review"
    );
    match &workflow.ops[0] {
        Operation::InterpretRnaReads { report_id, .. } => {
            assert_eq!(report_id.as_deref(), Some("tp73_manual_review"));
        }
        other => panic!("unexpected workflow op: {other:?}"),
    }
    assert_eq!(area.rna_reads_ui.report_id, "tp73_manual_review");
}

#[test]
fn legacy_rna_read_ui_state_keeps_manual_report_id_after_load() {
    let legacy_json = serde_json::json!({
        "input_path": "/tmp/SRR32957124.fasta.gz",
        "report_id": "tp73_manual_review"
    });
    let mut loaded =
        serde_json::from_value::<RnaReadInterpretOpsUiState>(legacy_json).expect("legacy ui");
    assert!(
        loaded.report_id_auto_sync,
        "legacy default still deserializes as auto"
    );
    assert!(
        !loaded.report_id_auto_sync_explicit,
        "legacy state should not carry the explicit auto/manual marker"
    );
    loaded.normalize_loaded_state();
    assert!(
        !loaded.report_id_auto_sync,
        "legacy manual report IDs should stay manual after normalization"
    );
    assert!(loaded.report_id_auto_sync_explicit);

    let dna = DNAsequence::from_sequence("ACGTACGT").expect("sequence");
    let mut area = MainAreaDna::new(dna, Some("TP73.ncbi".to_string()), None);
    area.rna_reads_ui = loaded;
    area.workflow_run_id.clear();
    let view = SplicingExpertView {
        seq_id: "TP73.ncbi".to_string(),
        target_feature_id: 42,
        scope: SplicingScopePreset::AllOverlappingAnyStrand,
        group_label: "TP73".to_string(),
        strand: "+".to_string(),
        region_start_1based: 1,
        region_end_1based: 8,
        transcript_count: 0,
        unique_exon_count: 0,
        instruction: String::new(),
        transcripts: vec![],
        unique_exons: vec![],
        matrix_rows: vec![],
        boundaries: vec![],
        intron_signals: vec![],
        junctions: vec![],
        events: vec![],
    };

    let workflow = area
        .build_splicing_rna_read_workflow(&view)
        .expect("workflow");
    assert_eq!(
        workflow.run_id,
        "workflow_rna_reads_tp73_ncbi_tp73_manual_review"
    );
    match &workflow.ops[0] {
        Operation::InterpretRnaReads { report_id, .. } => {
            assert_eq!(report_id.as_deref(), Some("tp73_manual_review"));
        }
        other => panic!("unexpected workflow op: {other:?}"),
    }
    assert_eq!(area.rna_reads_ui.report_id, "tp73_manual_review");
}

#[test]
fn explicit_auto_rna_read_ui_state_stays_auto_after_load() {
    let mut ui_state = RnaReadInterpretOpsUiState::default();
    ui_state.input_path = "/tmp/SRR32957124.fasta.gz".to_string();
    ui_state.report_id = "rna_tp73_srr32957124_ncdna_aoas_sg".to_string();
    let serialized = serde_json::to_value(&ui_state).expect("serialize ui");
    let mut loaded =
        serde_json::from_value::<RnaReadInterpretOpsUiState>(serialized).expect("reload ui");
    loaded.normalize_loaded_state();
    assert!(loaded.report_id_auto_sync);
    assert!(loaded.report_id_auto_sync_explicit);
}

#[test]
fn apply_rna_reads_demo_specificity_preset_sets_expected_seed_gates() {
    let dna = DNAsequence::from_sequence("ACGT").unwrap();
    let mut area = MainAreaDna::new(dna, None, None);
    area.apply_rna_reads_demo_specificity_preset();
    assert_eq!(
        area.rna_reads_ui.scope,
        SplicingScopePreset::TargetGroupTargetStrand
    );
    assert_eq!(area.rna_reads_ui.kmer_len, "10");
    assert_eq!(area.rna_reads_ui.seed_stride_bp, "1");
    assert!(area.rna_reads_ui.cdna_poly_t_flip_enabled);
    assert_eq!(area.rna_reads_ui.min_seed_hit_fraction, "0.30");
    assert_eq!(area.rna_reads_ui.min_weighted_seed_hit_fraction, "0.07");
    assert_eq!(area.rna_reads_ui.min_unique_matched_kmers, "16");
    assert_eq!(area.rna_reads_ui.min_chain_consistency_fraction, "0.60");
    assert_eq!(area.rna_reads_ui.max_median_transcript_gap, "3.0");
    assert_eq!(area.rna_reads_ui.min_confirmed_exon_transitions, "1");
    assert_eq!(area.rna_reads_ui.min_transition_support_fraction, "0.10");
    assert!(area.rna_reads_ui.show_advanced);
}

#[test]
fn apply_rna_read_dense_similarity_preset_sets_hash_and_dotplot_knobs() {
    let dna = DNAsequence::from_sequence("ACGT").unwrap();
    let mut area = MainAreaDna::new(dna, None, None);
    area.apply_rna_read_dense_similarity_preset();
    assert_eq!(area.rna_reads_ui.kmer_len, "9");
    assert_eq!(area.rna_reads_ui.seed_stride_bp, "1");
    assert_eq!(area.dotplot_ui.word_size, "9");
    assert_eq!(area.dotplot_ui.step_bp, "1");
    assert_eq!(area.dotplot_ui.max_mismatches, "0");
    assert!(area.dotplot_ui.tile_bp.is_empty());
    assert!(area.rna_reads_ui.show_advanced);
}

#[test]
fn rna_read_overlap_summaries_state_overlap_and_order_depth() {
    let dna = DNAsequence::from_sequence("ACGT").unwrap();
    let mut area = MainAreaDna::new(dna, None, None);
    area.rna_reads_ui.kmer_len = "9".to_string();
    area.rna_reads_ui.seed_stride_bp = "1".to_string();
    area.dotplot_ui.word_size = "9".to_string();
    area.dotplot_ui.step_bp = "3".to_string();
    let hash = area.rna_read_hash_parameter_summary();
    let dotplot = area.rna_read_dotplot_parameter_summary();
    assert!(hash.contains("overlap by 8 bp"));
    assert!(hash.contains("9 consecutive ordered windows"));
    assert!(dotplot.contains("overlap by 6 bp"));
    assert!(dotplot.contains("up to 3 consecutive ordered windows"));
}

#[test]
fn rna_read_statistics_tab_defaults_to_thresholded_cdna() {
    let dna = DNAsequence::from_sequence("ACGT").unwrap();
    let area = MainAreaDna::new(dna, None, None);
    assert_eq!(
        area.rna_read_statistics_tab,
        super::RnaReadEvidenceSourceTab::ThresholdedCdna
    );
}

#[test]
fn rna_read_mapped_cdna_subview_defaults_to_read_effects() {
    let dna = DNAsequence::from_sequence("ACGT").unwrap();
    let area = MainAreaDna::new(dna, None, None);
    assert_eq!(
        area.rna_read_mapped_cdna_subview,
        super::RnaReadMappedCdnaSubview::ReadEffects
    );
}

#[test]
fn rna_read_alignment_effect_controls_default_to_all_rank_no_search() {
    let dna = DNAsequence::from_sequence("ACGT").unwrap();
    let area = MainAreaDna::new(dna, None, None);
    assert_eq!(
        area.rna_read_alignment_effect_filter,
        super::RnaReadAlignmentEffectFilter::AllAligned
    );
    assert_eq!(
        area.rna_read_alignment_effect_sort_key,
        super::RnaReadAlignmentEffectSortKey::Rank
    );
    assert!(area.rna_read_alignment_effect_search.is_empty());
    assert_eq!(area.rna_read_alignment_effect_score_bin_index, None);
}

#[test]
fn saved_rna_read_report_payloads_are_cached_and_invalidated_on_commit() {
    let dna = DNAsequence::from_sequence("ACGTACGT").expect("sequence");
    let mut state = ProjectState::default();
    state.sequences.insert("seq1".to_string(), dna.clone());
    let engine = Arc::new(RwLock::new(GentleEngine::from_state(state)));
    let report = RnaReadInterpretationReport {
        report_id: "rna_cache_demo".to_string(),
        seq_id: "seq1".to_string(),
        seed_feature_id: 9,
        read_count_total: 1,
        read_count_seed_passed: 1,
        read_count_aligned: 1,
        score_density_bins: vec![0; 40],
        hits: vec![RnaReadInterpretationHit {
            record_index: 0,
            header_id: "read0".to_string(),
            sequence: "ACGTACGT".to_string(),
            read_length_bp: 8,
            tested_kmers: 8,
            matched_kmers: 8,
            seed_hit_fraction: 1.0,
            passed_seed_filter: true,
            best_mapping: Some(RnaReadMappingHit {
                transcript_feature_id: 9,
                transcript_id: "TX1".to_string(),
                transcript_label: "TX1".to_string(),
                strand: "+".to_string(),
                query_end_0based_exclusive: 8,
                target_start_1based: 1,
                target_end_1based: 8,
                target_end_offset_0based_exclusive: 8,
                matches: 8,
                score: 16,
                identity_fraction: 1.0,
                query_coverage_fraction: 1.0,
                ..RnaReadMappingHit::default()
            }),
            ..RnaReadInterpretationHit::default()
        }],
        ..RnaReadInterpretationReport::default()
    };
    engine
        .write()
        .expect("engine")
        .commit_rna_read_report(report.clone())
        .expect("commit report");
    let mut area = MainAreaDna::new(dna, Some("seq1".to_string()), Some(engine));
    area.rna_read_evidence_ui.selected_report_id = report.report_id.clone();
    let view = SplicingExpertView {
        seq_id: "seq1".to_string(),
        target_feature_id: 9,
        scope: SplicingScopePreset::AllOverlappingAnyStrand,
        group_label: "demo".to_string(),
        strand: "+".to_string(),
        region_start_1based: 1,
        region_end_1based: 8,
        transcript_count: 1,
        unique_exon_count: 0,
        instruction: String::new(),
        transcripts: vec![],
        unique_exons: vec![],
        matrix_rows: vec![],
        boundaries: vec![],
        intron_signals: vec![],
        junctions: vec![],
        events: vec![],
    };

    let summaries_a = area.matching_rna_read_report_summaries_for_splicing_view(&view);
    let summaries_b = area.matching_rna_read_report_summaries_for_splicing_view(&view);
    assert!(Arc::ptr_eq(&summaries_a, &summaries_b));

    let saved_report_a = area.current_saved_rna_read_report().expect("saved report");
    let saved_report_b = area
        .current_saved_rna_read_report()
        .expect("saved report cache");
    assert!(Arc::ptr_eq(&saved_report_a, &saved_report_b));

    let progress_a = area
        .current_rna_read_evidence_progress_for_view(&view, Some(saved_report_a.as_ref()))
        .expect("cached progress");
    let progress_b = area
        .current_rna_read_evidence_progress_for_view(&view, Some(saved_report_a.as_ref()))
        .expect("cached progress");
    assert!(Arc::ptr_eq(&progress_a, &progress_b));

    let inspection_a = area
        .current_saved_rna_read_alignment_inspection(saved_report_a.hits.len(), None)
        .expect("inspection");
    let inspection_b = area
        .current_saved_rna_read_alignment_inspection(saved_report_a.hits.len(), None)
        .expect("inspection cache");
    assert!(Arc::ptr_eq(&inspection_a, &inspection_b));

    let mut updated_report = report.clone();
    updated_report.read_count_total = 2;
    updated_report.read_count_seed_passed = 2;
    updated_report.hits.push(RnaReadInterpretationHit {
        record_index: 1,
        header_id: "read1".to_string(),
        sequence: "ACGTACGA".to_string(),
        read_length_bp: 8,
        tested_kmers: 8,
        matched_kmers: 7,
        seed_hit_fraction: 0.875,
        passed_seed_filter: true,
        ..RnaReadInterpretationHit::default()
    });
    area.commit_completed_rna_read_task_outcome(RnaReadTaskOutcome::Interpret(updated_report))
        .expect("commit updated report");

    let summaries_c = area.matching_rna_read_report_summaries_for_splicing_view(&view);
    assert!(!Arc::ptr_eq(&summaries_a, &summaries_c));
    let saved_report_c = area
        .current_saved_rna_read_report()
        .expect("refreshed saved report");
    assert!(!Arc::ptr_eq(&saved_report_a, &saved_report_c));
    let progress_c = area
        .current_rna_read_evidence_progress_for_view(&view, Some(saved_report_c.as_ref()))
        .expect("refreshed progress");
    assert!(!Arc::ptr_eq(&progress_a, &progress_c));
    assert_eq!(progress_c.reads_processed, 2);
}

#[test]
fn mapping_workspace_uses_saved_report_progress_without_evidence_selection() {
    let dna = DNAsequence::from_sequence("ACGTACGT").expect("sequence");
    let mut state = ProjectState::default();
    state.sequences.insert("seq1".to_string(), dna.clone());
    let engine = Arc::new(RwLock::new(GentleEngine::from_state(state)));
    let report = RnaReadInterpretationReport {
        report_id: "rna_workspace_demo".to_string(),
        seq_id: "seq1".to_string(),
        seed_feature_id: 9,
        read_count_total: 1,
        read_count_seed_passed: 1,
        read_count_aligned: 1,
        score_density_bins: {
            let mut bins = vec![0; 40];
            bins[39] = 1;
            bins
        },
        hits: vec![RnaReadInterpretationHit {
            record_index: 0,
            header_id: "read0".to_string(),
            sequence: "ACGTACGT".to_string(),
            read_length_bp: 8,
            tested_kmers: 8,
            matched_kmers: 8,
            seed_hit_fraction: 1.0,
            weighted_seed_hit_fraction: 1.0,
            passed_seed_filter: true,
            best_mapping: Some(RnaReadMappingHit {
                transcript_feature_id: 9,
                transcript_id: "TX1".to_string(),
                transcript_label: "TX1".to_string(),
                strand: "+".to_string(),
                query_end_0based_exclusive: 8,
                target_start_1based: 1,
                target_end_1based: 8,
                target_end_offset_0based_exclusive: 8,
                matches: 8,
                score: 16,
                identity_fraction: 1.0,
                query_coverage_fraction: 1.0,
                ..RnaReadMappingHit::default()
            }),
            ..RnaReadInterpretationHit::default()
        }],
        ..RnaReadInterpretationReport::default()
    };
    engine
        .write()
        .expect("engine")
        .commit_rna_read_report(report)
        .expect("commit report");
    let mut area = MainAreaDna::new(dna, Some("seq1".to_string()), Some(engine));
    area.rna_reads_ui.report_id = "rna_workspace_demo".to_string();
    area.rna_read_evidence_ui.selected_report_id.clear();
    let view = SplicingExpertView {
        seq_id: "seq1".to_string(),
        target_feature_id: 9,
        scope: SplicingScopePreset::AllOverlappingAnyStrand,
        group_label: "demo".to_string(),
        strand: "+".to_string(),
        region_start_1based: 1,
        region_end_1based: 8,
        transcript_count: 1,
        unique_exon_count: 0,
        instruction: String::new(),
        transcripts: vec![],
        unique_exons: vec![],
        matrix_rows: vec![],
        boundaries: vec![],
        intron_signals: vec![],
        junctions: vec![],
        events: vec![],
    };

    area.sync_rna_read_evidence_selection_to_mapping_report();
    assert_eq!(
        area.rna_read_evidence_ui.selected_report_id,
        "rna_workspace_demo"
    );

    let progress = area
        .current_rna_read_mapping_progress_for_view(&view)
        .expect("workspace progress");
    assert_eq!(progress.reads_processed, 1);
    assert_eq!(progress.aligned, 1);
    assert_eq!(progress.score_density_bins[39], 1);
    assert_eq!(progress.top_hits_preview.len(), 1);
}

#[test]
fn sync_rna_read_evidence_selection_to_mapping_report_clears_report_scoped_row_state() {
    let dna = DNAsequence::from_sequence("ACGTACGT").expect("sequence");
    let mut area = MainAreaDna::new(dna, Some("seq1".to_string()), None);
    area.rna_read_evidence_ui.selected_report_id = "old_report".to_string();
    area.rna_reads_ui.report_id = "new_report".to_string();
    area.rna_seed_selected_record_indices = BTreeSet::from([2usize, 7usize]);
    area.rna_seed_highlight_record_index = Some(7);
    area.rna_read_alignment_effect_score_bin_index = Some(39);
    area.rna_read_alignment_detail_visible_key = Some("old_report|7".to_string());

    area.sync_rna_read_evidence_selection_to_mapping_report();

    assert_eq!(area.rna_read_evidence_ui.selected_report_id, "new_report");
    assert!(area.rna_seed_selected_record_indices.is_empty());
    assert_eq!(area.rna_seed_highlight_record_index, None);
    assert_eq!(area.rna_read_alignment_effect_score_bin_index, None);
    assert_eq!(area.rna_read_alignment_detail_visible_key, None);
}

#[test]
fn mapping_workspace_ignores_saved_report_progress_for_other_splicing_view() {
    let dna = DNAsequence::from_sequence("ACGTACGT").expect("sequence");
    let mut state = ProjectState::default();
    state.sequences.insert("seq1".to_string(), dna.clone());
    let engine = Arc::new(RwLock::new(GentleEngine::from_state(state)));
    engine
        .write()
        .expect("engine")
        .commit_rna_read_report(RnaReadInterpretationReport {
            report_id: "other_locus".to_string(),
            seq_id: "seq1".to_string(),
            seed_feature_id: 7,
            read_count_total: 1,
            read_count_seed_passed: 1,
            hits: vec![RnaReadInterpretationHit {
                record_index: 0,
                header_id: "read0".to_string(),
                sequence: "ACGTACGT".to_string(),
                read_length_bp: 8,
                tested_kmers: 8,
                matched_kmers: 8,
                seed_hit_fraction: 1.0,
                passed_seed_filter: true,
                ..RnaReadInterpretationHit::default()
            }],
            ..RnaReadInterpretationReport::default()
        })
        .expect("commit report");
    let mut area = MainAreaDna::new(dna, Some("seq1".to_string()), Some(engine));
    area.rna_reads_ui.report_id = "other_locus".to_string();
    area.rna_reads_ui.report_id_auto_sync = false;
    area.rna_reads_ui.input_path = "/tmp/reads.fa".to_string();
    let view = SplicingExpertView {
        seq_id: "seq1".to_string(),
        target_feature_id: 9,
        scope: SplicingScopePreset::AllOverlappingAnyStrand,
        group_label: "demo".to_string(),
        strand: "+".to_string(),
        region_start_1based: 1,
        region_end_1based: 8,
        transcript_count: 1,
        unique_exon_count: 0,
        instruction: String::new(),
        transcripts: vec![],
        unique_exons: vec![],
        matrix_rows: vec![],
        boundaries: vec![],
        intron_signals: vec![],
        junctions: vec![],
        events: vec![],
    };

    assert!(
        area.current_rna_read_mapping_progress_for_view(&view)
            .is_none()
    );
    let align_error = area
        .build_splicing_rna_read_align_operation(&view, None)
        .expect_err("stale report id should not align from this workspace");
    assert!(align_error.contains("belongs to seq1 feature n-7"));
    let interpret_error = area
        .build_splicing_rna_read_interpret_operation(&view)
        .expect_err("stale report id should not be overwritten from this workspace");
    assert!(interpret_error.contains("belongs to seq1 feature n-7"));
}

#[test]
fn unaligned_saved_report_progress_skips_alignment_inspection() {
    let dna = DNAsequence::from_sequence("ACGTACGT").expect("sequence");
    let mut state = ProjectState::default();
    state.sequences.insert("seq1".to_string(), dna.clone());
    let engine = Arc::new(RwLock::new(GentleEngine::from_state(state)));
    let report = RnaReadInterpretationReport {
        report_id: "rna_unaligned_demo".to_string(),
        seq_id: "seq1".to_string(),
        seed_feature_id: 9,
        read_count_total: 1,
        read_count_seed_passed: 1,
        hits: vec![RnaReadInterpretationHit {
            record_index: 0,
            header_id: "read0".to_string(),
            sequence: "ACGTACGT".to_string(),
            read_length_bp: 8,
            tested_kmers: 8,
            matched_kmers: 8,
            seed_hit_fraction: 1.0,
            passed_seed_filter: true,
            ..RnaReadInterpretationHit::default()
        }],
        ..RnaReadInterpretationReport::default()
    };
    engine
        .write()
        .expect("engine")
        .commit_rna_read_report(report.clone())
        .expect("commit report");
    let mut area = MainAreaDna::new(dna, Some("seq1".to_string()), Some(engine));
    area.rna_read_evidence_ui.selected_report_id = report.report_id.clone();
    let view = SplicingExpertView {
        seq_id: "seq1".to_string(),
        target_feature_id: 9,
        scope: SplicingScopePreset::AllOverlappingAnyStrand,
        group_label: "demo".to_string(),
        strand: "+".to_string(),
        region_start_1based: 1,
        region_end_1based: 8,
        transcript_count: 1,
        unique_exon_count: 0,
        instruction: String::new(),
        transcripts: vec![],
        unique_exons: vec![],
        matrix_rows: vec![],
        boundaries: vec![],
        intron_signals: vec![],
        junctions: vec![],
        events: vec![],
    };

    let saved_report = area.current_saved_rna_read_report().expect("saved report");
    let progress = area
        .current_rna_read_evidence_progress_for_view(&view, Some(saved_report.as_ref()))
        .expect("saved progress");

    assert!(progress.mapped_exon_support_frequencies.is_empty());
    assert!(progress.mapped_junction_support_frequencies.is_empty());
    assert!(area.cached_rna_read_alignment_inspections.is_empty());
}

#[test]
fn rna_read_alignment_inspection_errors_are_cached() {
    let dna = DNAsequence::from_sequence("ACGT").expect("sequence");
    let mut state = ProjectState::default();
    state.sequences.insert("seq1".to_string(), dna.clone());
    let engine = Arc::new(RwLock::new(GentleEngine::from_state(state)));
    let mut area = MainAreaDna::new(dna, Some("seq1".to_string()), Some(engine));

    let err_a = area
        .saved_rna_read_alignment_inspection_for_report_id("missing_report", 1, None)
        .expect_err("missing report should fail");
    assert_eq!(area.cached_rna_read_alignment_inspections.len(), 1);
    match &area.cached_rna_read_alignment_inspections[0].result {
        Err(cached) => assert_eq!(cached, &err_a),
        Ok(_) => panic!("expected cached alignment inspection error"),
    }

    let err_b = area
        .saved_rna_read_alignment_inspection_for_report_id("missing_report", 1, None)
        .expect_err("missing report should still fail");
    assert_eq!(err_a, err_b);
    assert_eq!(area.cached_rna_read_alignment_inspections.len(), 1);
}

#[test]
fn rna_read_alignment_display_errors_are_cached() {
    let dna = DNAsequence::from_sequence("ACGTACGT").expect("sequence");
    let mut state = ProjectState::default();
    state.sequences.insert("seq1".to_string(), dna.clone());
    let engine = Arc::new(RwLock::new(GentleEngine::from_state(state)));
    let report = RnaReadInterpretationReport {
        report_id: "rna_alignment_display_err".to_string(),
        seq_id: "seq1".to_string(),
        seed_feature_id: 9,
        read_count_total: 1,
        read_count_seed_passed: 1,
        hits: vec![RnaReadInterpretationHit {
            record_index: 0,
            header_id: "read0".to_string(),
            sequence: "ACGTACGT".to_string(),
            read_length_bp: 8,
            tested_kmers: 8,
            matched_kmers: 8,
            seed_hit_fraction: 1.0,
            passed_seed_filter: true,
            ..RnaReadInterpretationHit::default()
        }],
        ..RnaReadInterpretationReport::default()
    };
    engine
        .write()
        .expect("engine")
        .commit_rna_read_report(report.clone())
        .expect("commit report");
    let mut area = MainAreaDna::new(dna, Some("seq1".to_string()), Some(engine));
    area.rna_read_evidence_ui.selected_report_id = report.report_id.clone();

    let err_a = area
        .current_rna_read_alignment_display(0)
        .expect_err("unaligned row should not have a pairwise display");
    let cached = area
        .cached_rna_read_alignment_display
        .as_ref()
        .expect("alignment display error should be cached");
    assert_eq!(
        cached.cache_key,
        MainAreaDna::rna_read_alignment_display_cache_key(&report.report_id, 0)
    );
    match &cached.result {
        Err(cached_err) => assert_eq!(cached_err, &err_a),
        Ok(_) => panic!("expected cached alignment display error"),
    }

    let err_b = area
        .current_rna_read_alignment_display(0)
        .expect_err("unaligned row should still not have a pairwise display");
    assert_eq!(err_a, err_b);
}

#[test]
fn select_rna_read_report_score_bin_record_indices_matches_requested_bin() {
    let report = RnaReadInterpretationReport {
        hits: vec![
            RnaReadInterpretationHit {
                record_index: 1,
                seed_hit_fraction: 0.12,
                ..RnaReadInterpretationHit::default()
            },
            RnaReadInterpretationHit {
                record_index: 4,
                seed_hit_fraction: 0.88,
                ..RnaReadInterpretationHit::default()
            },
            RnaReadInterpretationHit {
                record_index: 7,
                seed_hit_fraction: 0.89,
                ..RnaReadInterpretationHit::default()
            },
        ],
        score_density_bins: vec![0; 40],
        ..RnaReadInterpretationReport::default()
    };
    assert_eq!(
        MainAreaDna::select_rna_read_report_score_bin_record_indices(
            &report,
            35,
            40,
            RnaReadScoreDensityVariant::AllScored,
            None,
        ),
        vec![4, 7]
    );
}

#[test]
fn select_rna_read_report_score_bin_record_indices_composite_gate_respects_passed_seed_filter() {
    let report = RnaReadInterpretationReport {
        hits: vec![
            RnaReadInterpretationHit {
                record_index: 4,
                seed_hit_fraction: 0.88,
                passed_seed_filter: false,
                ..RnaReadInterpretationHit::default()
            },
            RnaReadInterpretationHit {
                record_index: 7,
                seed_hit_fraction: 0.89,
                passed_seed_filter: true,
                ..RnaReadInterpretationHit::default()
            },
        ],
        score_density_bins: vec![0; 40],
        seed_pass_score_density_bins: vec![0; 40],
        ..RnaReadInterpretationReport::default()
    };
    assert_eq!(
        MainAreaDna::select_rna_read_report_score_bin_record_indices(
            &report,
            35,
            40,
            RnaReadScoreDensityVariant::CompositeSeedGate,
            None,
        ),
        vec![7]
    );
}

#[test]
fn select_rna_read_report_score_bin_record_indices_retained_replay_respects_override() {
    let replay_filter = RnaReadSeedFilterConfig {
        min_seed_hit_fraction: 0.80,
        min_weighted_seed_hit_fraction: 0.80,
        min_unique_matched_kmers: 10,
        min_chain_consistency_fraction: 0.50,
        max_median_transcript_gap: 4.0,
        min_confirmed_exon_transitions: 1,
        min_transition_support_fraction: 0.50,
        ..RnaReadSeedFilterConfig::default()
    };
    let report = RnaReadInterpretationReport {
        hits: vec![
            RnaReadInterpretationHit {
                record_index: 4,
                seed_hit_fraction: 0.88,
                weighted_seed_hit_fraction: 0.88,
                tested_kmers: 12,
                matched_kmers: 12,
                seed_chain_support_fraction: 1.0,
                seed_median_transcript_gap: 1.0,
                seed_transcript_gap_count: 1,
                exon_transitions_confirmed: 1,
                exon_transitions_total: 1,
                ..RnaReadInterpretationHit::default()
            },
            RnaReadInterpretationHit {
                record_index: 7,
                seed_hit_fraction: 0.89,
                weighted_seed_hit_fraction: 0.89,
                tested_kmers: 12,
                matched_kmers: 12,
                seed_chain_support_fraction: 0.05,
                seed_median_transcript_gap: 1.0,
                seed_transcript_gap_count: 1,
                exon_transitions_confirmed: 1,
                exon_transitions_total: 1,
                ..RnaReadInterpretationHit::default()
            },
        ],
        score_density_bins: vec![0; 40],
        ..RnaReadInterpretationReport::default()
    };
    assert_eq!(
        MainAreaDna::select_rna_read_report_score_bin_record_indices(
            &report,
            35,
            40,
            RnaReadScoreDensityVariant::RetainedReplayCurrentControls,
            Some(&replay_filter),
        ),
        vec![4]
    );
}

#[test]
fn rna_read_alignment_selection_summary_explains_all_vs_seed_passed_counts() {
    let report = RnaReadInterpretationReport {
        hits: vec![
            RnaReadInterpretationHit {
                record_index: 1,
                seed_hit_fraction: 0.12,
                passed_seed_filter: false,
                ..RnaReadInterpretationHit::default()
            },
            RnaReadInterpretationHit {
                record_index: 4,
                seed_hit_fraction: 0.31,
                passed_seed_filter: false,
                ..RnaReadInterpretationHit::default()
            },
            RnaReadInterpretationHit {
                record_index: 7,
                seed_hit_fraction: 0.88,
                passed_seed_filter: true,
                ..RnaReadInterpretationHit::default()
            },
        ],
        seed_filter: RnaReadSeedFilterConfig {
            min_seed_hit_fraction: 0.30,
            ..RnaReadSeedFilterConfig::default()
        },
        ..RnaReadInterpretationReport::default()
    };

    let summary_all =
        MainAreaDna::format_rna_read_alignment_selection_summary(&report, RnaReadHitSelection::All);
    assert!(summary_all.contains("aligns all retained saved-report rows (3)"));
    assert!(
        summary_all.contains("1 currently pass the composite seed gate recorded for this report")
    );
    assert!(summary_all.contains("2 are at or above raw min_hit=0.30"));
    assert!(summary_all.contains("The red histogram line marks only that raw min_hit component"));

    let summary_seed_passed = MainAreaDna::format_rna_read_alignment_selection_summary(
        &report,
        RnaReadHitSelection::SeedPassed,
    );
    assert!(
        summary_seed_passed.contains("aligns the retained saved-report rows that currently pass")
    );
    assert!(summary_seed_passed.contains("(1)"));
    assert!(summary_seed_passed.contains("passed_seed_filter=yes"));
}

#[test]
fn rna_read_support_and_preview_tables_reserve_practical_minimum_heights() {
    let ctx = egui::Context::default();
    ctx.begin_pass(egui::RawInput::default());
    egui::Area::new("rna_read_table_height_probe".into()).show(&ctx, |ui| {
        let support_height = MainAreaDna::default_rna_read_support_table_height(ui);
        let preview_height = MainAreaDna::default_rna_read_preview_table_height_for_rows(ui, 12);

        assert!(
            support_height >= 170.0,
            "support tables should reserve about eight visible rows"
        );
        assert!(
            preview_height > support_height,
            "preview table should reserve at least as much space as the support tables"
        );
    });
    let _ = ctx.end_pass();
}

#[test]
fn collect_rna_read_top_hit_previews_for_score_bin_uses_all_saved_report_rows() {
    let report = RnaReadInterpretationReport {
        hits: vec![
            RnaReadInterpretationHit {
                record_index: 1,
                header_id: "r1".to_string(),
                seed_hit_fraction: 0.12,
                ..RnaReadInterpretationHit::default()
            },
            RnaReadInterpretationHit {
                record_index: 4,
                header_id: "r4".to_string(),
                seed_hit_fraction: 0.88,
                ..RnaReadInterpretationHit::default()
            },
            RnaReadInterpretationHit {
                record_index: 7,
                header_id: "r7".to_string(),
                seed_hit_fraction: 0.89,
                ..RnaReadInterpretationHit::default()
            },
        ],
        score_density_bins: vec![0; 40],
        ..RnaReadInterpretationReport::default()
    };
    let rows = MainAreaDna::collect_rna_read_top_hit_previews_for_score_bin(
        &report,
        35,
        40,
        RnaReadScoreDensityVariant::AllScored,
        None,
    );
    assert_eq!(
        rows.iter().map(|row| row.record_index).collect::<Vec<_>>(),
        vec![7, 4]
    );
}

#[test]
fn rna_read_selected_export_operations_preserve_exact_record_indices() {
    let selected_record_indices = vec![2usize, 5usize, 9usize];

    let exon_paths = super::RnaReadSelectedExportKind::ExonPathsTsv.build_operation(
        "report1".to_string(),
        "/tmp/paths.tsv".to_string(),
        selected_record_indices.clone(),
        Some("filter=disagreement only | sort=score | search=tp53".to_string()),
    );
    match exon_paths {
        Operation::ExportRnaReadExonPathsTsv {
            report_id,
            path,
            selection,
            selected_record_indices,
            subset_spec,
        } => {
            assert_eq!(report_id, "report1");
            assert_eq!(path, "/tmp/paths.tsv");
            assert_eq!(selection, RnaReadHitSelection::Aligned);
            assert_eq!(selected_record_indices, vec![2, 5, 9]);
            assert_eq!(
                subset_spec.as_deref(),
                Some("filter=disagreement only | sort=score | search=tp53")
            );
        }
        other => panic!("unexpected operation: {other:?}"),
    }

    let abundance = super::RnaReadSelectedExportKind::ExonAbundanceTsv.build_operation(
        "report2".to_string(),
        "/tmp/abundance.tsv".to_string(),
        selected_record_indices,
        None,
    );
    match abundance {
        Operation::ExportRnaReadExonAbundanceTsv {
            report_id,
            path,
            selection,
            selected_record_indices,
            subset_spec,
        } => {
            assert_eq!(report_id, "report2");
            assert_eq!(path, "/tmp/abundance.tsv");
            assert_eq!(selection, RnaReadHitSelection::Aligned);
            assert_eq!(selected_record_indices, vec![2, 5, 9]);
            assert!(subset_spec.is_none());
        }
        other => panic!("unexpected operation: {other:?}"),
    }
}

#[test]
fn build_splicing_rna_read_align_operation_uses_explicit_sorted_subset() {
    let dna = DNAsequence::from_sequence("ACGT").expect("sequence");
    let mut area = MainAreaDna::new(dna, None, None);
    area.rna_reads_ui.report_id = "report1".to_string();
    area.rna_reads_ui.align_phase_selection = RnaReadHitSelection::All;
    let view = SplicingExpertView {
        seq_id: "seq1".to_string(),
        target_feature_id: 9,
        scope: SplicingScopePreset::AllOverlappingAnyStrand,
        group_label: "demo".to_string(),
        strand: "+".to_string(),
        region_start_1based: 1,
        region_end_1based: 4,
        transcript_count: 1,
        unique_exon_count: 0,
        instruction: String::new(),
        transcripts: vec![],
        unique_exons: vec![],
        matrix_rows: vec![],
        boundaries: vec![],
        intron_signals: vec![],
        junctions: vec![],
        events: vec![],
    };

    let op = area
        .build_splicing_rna_read_align_operation(&view, Some(vec![9, 2, 9, 5]))
        .expect("build align op");

    match op {
        Operation::AlignRnaReadReport {
            report_id,
            selection,
            align_config_override,
            selected_record_indices,
        } => {
            assert_eq!(report_id, "report1");
            assert_eq!(selection, RnaReadHitSelection::All);
            assert_eq!(selected_record_indices, vec![2, 5, 9]);
            assert_eq!(
                align_config_override.expect("align config").band_width_bp,
                24
            );
        }
        other => panic!("unexpected operation: {other:?}"),
    }
}

#[test]
fn summarize_rna_read_alignment_effects_counts_each_category() {
    let report = RnaReadInterpretationReport {
        report_id: "rna_effects".to_string(),
        hits: vec![
            RnaReadInterpretationHit {
                record_index: 0,
                passed_seed_filter: true,
                best_mapping: Some(RnaReadMappingHit::default()),
                ..RnaReadInterpretationHit::default()
            },
            RnaReadInterpretationHit {
                record_index: 1,
                passed_seed_filter: true,
                best_mapping: Some(RnaReadMappingHit::default()),
                ..RnaReadInterpretationHit::default()
            },
            RnaReadInterpretationHit {
                record_index: 2,
                passed_seed_filter: true,
                best_mapping: Some(RnaReadMappingHit::default()),
                ..RnaReadInterpretationHit::default()
            },
            RnaReadInterpretationHit {
                record_index: 3,
                passed_seed_filter: true,
                best_mapping: None,
                ..RnaReadInterpretationHit::default()
            },
        ],
        ..RnaReadInterpretationReport::default()
    };
    let inspection = RnaReadAlignmentInspection {
        aligned_count: 3,
        rows: vec![
            RnaReadAlignmentInspectionRow {
                alignment_effect: RnaReadAlignmentEffect::ConfirmedAssignment,
                ..RnaReadAlignmentInspectionRow::default()
            },
            RnaReadAlignmentInspectionRow {
                alignment_effect: RnaReadAlignmentEffect::ReassignedTranscript,
                ..RnaReadAlignmentInspectionRow::default()
            },
            RnaReadAlignmentInspectionRow {
                alignment_effect: RnaReadAlignmentEffect::AlignedWithoutPhase1Assignment,
                ..RnaReadAlignmentInspectionRow::default()
            },
        ],
        ..RnaReadAlignmentInspection::default()
    };

    let summary = MainAreaDna::summarize_rna_read_alignment_effects(&report, &inspection);
    assert_eq!(summary.aligned_rows, 3);
    assert_eq!(summary.confirmed_assignments, 1);
    assert_eq!(summary.reassigned_transcripts, 1);
    assert_eq!(summary.aligned_without_phase1_assignment, 1);
    assert_eq!(summary.seed_passed_but_unaligned, 1);
}

#[test]
fn collect_visible_rna_read_alignment_effect_rows_filters_sorts_and_searches() {
    let inspection = RnaReadAlignmentInspection {
        rows: vec![
            RnaReadAlignmentInspectionRow {
                record_index: 0,
                rank: 3,
                header_id: "read_confirmed".to_string(),
                phase1_primary_transcript_id: "TXA".to_string(),
                transcript_id: "TXA".to_string(),
                transcript_label: "TXA".to_string(),
                identity_fraction: 0.81,
                query_coverage_fraction: 0.91,
                score: 120,
                alignment_effect: RnaReadAlignmentEffect::ConfirmedAssignment,
                ..RnaReadAlignmentInspectionRow::default()
            },
            RnaReadAlignmentInspectionRow {
                record_index: 1,
                rank: 1,
                header_id: "read_reassigned".to_string(),
                phase1_primary_transcript_id: "TX_seed".to_string(),
                transcript_id: "TXB".to_string(),
                transcript_label: "TP53-201".to_string(),
                identity_fraction: 0.95,
                query_coverage_fraction: 0.62,
                score: 310,
                alignment_effect: RnaReadAlignmentEffect::ReassignedTranscript,
                ..RnaReadAlignmentInspectionRow::default()
            },
            RnaReadAlignmentInspectionRow {
                record_index: 2,
                rank: 2,
                header_id: "orphan_aligned".to_string(),
                transcript_id: "TXC".to_string(),
                transcript_label: "GFOD3P".to_string(),
                identity_fraction: 0.89,
                query_coverage_fraction: 0.99,
                score: 205,
                alignment_effect: RnaReadAlignmentEffect::AlignedWithoutPhase1Assignment,
                ..RnaReadAlignmentInspectionRow::default()
            },
        ],
        ..RnaReadAlignmentInspection::default()
    };
    let selected = BTreeSet::from([2usize]);

    let by_score = MainAreaDna::collect_visible_rna_read_alignment_effect_rows(
        &inspection,
        super::RnaReadAlignmentEffectFilter::AllAligned,
        "",
        &selected,
        super::RnaReadAlignmentEffectSortKey::Score,
    );
    assert_eq!(
        by_score
            .iter()
            .map(|row| row.record_index)
            .collect::<Vec<_>>(),
        vec![1, 2, 0]
    );

    let confirmed_only = MainAreaDna::collect_visible_rna_read_alignment_effect_rows(
        &inspection,
        super::RnaReadAlignmentEffectFilter::ConfirmedOnly,
        "",
        &selected,
        super::RnaReadAlignmentEffectSortKey::Rank,
    );
    assert_eq!(
        confirmed_only
            .iter()
            .map(|row| row.record_index)
            .collect::<Vec<_>>(),
        vec![0]
    );

    let disagreements = MainAreaDna::collect_visible_rna_read_alignment_effect_rows(
        &inspection,
        super::RnaReadAlignmentEffectFilter::DisagreementOnly,
        "",
        &selected,
        super::RnaReadAlignmentEffectSortKey::Score,
    );
    assert_eq!(
        disagreements
            .iter()
            .map(|row| row.record_index)
            .collect::<Vec<_>>(),
        vec![1, 2]
    );

    let selected_only = MainAreaDna::collect_visible_rna_read_alignment_effect_rows(
        &inspection,
        super::RnaReadAlignmentEffectFilter::SelectedOnly,
        "",
        &selected,
        super::RnaReadAlignmentEffectSortKey::Rank,
    );
    assert_eq!(
        selected_only
            .iter()
            .map(|row| row.record_index)
            .collect::<Vec<_>>(),
        vec![2]
    );

    let transcript_search = MainAreaDna::collect_visible_rna_read_alignment_effect_rows(
        &inspection,
        super::RnaReadAlignmentEffectFilter::AllAligned,
        "tp53",
        &selected,
        super::RnaReadAlignmentEffectSortKey::Rank,
    );
    assert_eq!(
        transcript_search
            .iter()
            .map(|row| row.record_index)
            .collect::<Vec<_>>(),
        vec![1]
    );
}

#[test]
fn collect_rna_read_alignment_effect_record_indices_preserves_display_order() {
    let inspection = RnaReadAlignmentInspection {
        rows: vec![
            RnaReadAlignmentInspectionRow {
                record_index: 7,
                ..RnaReadAlignmentInspectionRow::default()
            },
            RnaReadAlignmentInspectionRow {
                record_index: 2,
                ..RnaReadAlignmentInspectionRow::default()
            },
            RnaReadAlignmentInspectionRow {
                record_index: 11,
                ..RnaReadAlignmentInspectionRow::default()
            },
        ],
        ..RnaReadAlignmentInspection::default()
    };
    let rows = inspection.rows.iter().collect::<Vec<_>>();
    assert_eq!(
        MainAreaDna::collect_rna_read_alignment_effect_record_indices(&rows),
        vec![7, 2, 11]
    );
}

#[test]
fn rna_read_record_indices_all_selected_requires_full_non_empty_coverage() {
    let dna = DNAsequence::from_sequence("ACGT").expect("sequence");
    let mut area = MainAreaDna::new(dna, None, None);
    area.rna_seed_selected_record_indices = [2usize, 7usize].into_iter().collect();

    assert!(area.rna_read_record_indices_all_selected([2, 7]));
    assert!(!area.rna_read_record_indices_all_selected([2, 7, 11]));
    assert!(!area.rna_read_record_indices_all_selected(std::iter::empty()));
}

#[test]
fn set_rna_read_record_indices_selected_only_changes_requested_rows() {
    let dna = DNAsequence::from_sequence("ACGT").expect("sequence");
    let mut area = MainAreaDna::new(dna, None, None);
    area.rna_seed_selected_record_indices = [2usize, 11usize].into_iter().collect();

    let affected = area.set_rna_read_record_indices_selected([7, 13], true);
    assert_eq!(affected, 2);
    assert_eq!(
        area.selected_rna_record_indices(),
        vec![2, 7, 11, 13],
        "selecting listed rows should preserve hidden selections",
    );

    let affected = area.set_rna_read_record_indices_selected([7, 13], false);
    assert_eq!(affected, 2);
    assert_eq!(
        area.selected_rna_record_indices(),
        vec![2, 11],
        "clearing listed rows should leave hidden selections intact",
    );
}

#[test]
fn rna_read_alignment_effect_subset_spec_is_formal_and_stable() {
    assert_eq!(
        MainAreaDna::rna_read_alignment_effect_subset_spec(
            super::RnaReadAlignmentEffectFilter::DisagreementOnly,
            "tp53",
            super::RnaReadAlignmentEffectSortKey::Score,
            RnaReadScoreDensityVariant::AllScored,
            Some(39),
            40,
        ),
        "filter=disagreement only | histogram=all scored | score_bin=39/40 [0.975,1.000) | sort=score | search=tp53"
    );
    assert_eq!(
        MainAreaDna::rna_read_alignment_effect_subset_spec(
            super::RnaReadAlignmentEffectFilter::AllAligned,
            "   ",
            super::RnaReadAlignmentEffectSortKey::Rank,
            RnaReadScoreDensityVariant::CompositeSeedGate,
            None,
            40,
        ),
        "filter=all aligned | histogram=composite gate | score_bin=<none> | sort=rank | search=<none>"
    );
}

#[test]
fn rna_read_alignment_effect_subset_struct_maps_to_engine_contract() {
    let selected = BTreeSet::from([9usize, 2usize, 2usize]);
    let subset = MainAreaDna::rna_read_alignment_effect_subset_struct(
        super::RnaReadAlignmentEffectFilter::SelectedOnly,
        " tp53 ",
        super::RnaReadAlignmentEffectSortKey::Score,
        &selected,
        RnaReadScoreDensityVariant::CompositeSeedGate,
        None,
        Some(7),
        40,
    );
    assert_eq!(
        subset.effect_filter,
        crate::engine::RnaReadAlignmentInspectionEffectFilter::SelectedOnly
    );
    assert_eq!(
        subset.sort_key,
        crate::engine::RnaReadAlignmentInspectionSortKey::Score
    );
    assert_eq!(subset.search, "tp53");
    assert_eq!(subset.selected_record_indices, vec![2, 9]);
    assert_eq!(
        subset.score_density_variant,
        crate::engine::RnaReadScoreDensityVariant::CompositeSeedGate
    );
    assert!(subset.score_density_seed_filter_override.is_none());
    assert_eq!(subset.score_bin_index, Some(7));
    assert_eq!(subset.score_bin_count, 40);
}

#[test]
fn collect_rna_read_mapped_support_contributors_tracks_record_indices() {
    let inspection = RnaReadAlignmentInspection {
        rows: vec![
            RnaReadAlignmentInspectionRow {
                record_index: 0,
                transcript_id: "TXA".to_string(),
                mapped_exon_support: vec![crate::engine::RnaReadMappedSupportExonAttribution {
                    start_1based: 10,
                    end_1based: 20,
                }],
                mapped_junction_support: vec![
                    crate::engine::RnaReadMappedSupportJunctionAttribution {
                        donor_1based: 20,
                        acceptor_1based: 30,
                    },
                ],
                ..RnaReadAlignmentInspectionRow::default()
            },
            RnaReadAlignmentInspectionRow {
                record_index: 2,
                transcript_id: "TXA".to_string(),
                mapped_exon_support: vec![crate::engine::RnaReadMappedSupportExonAttribution {
                    start_1based: 10,
                    end_1based: 20,
                }],
                mapped_junction_support: vec![
                    crate::engine::RnaReadMappedSupportJunctionAttribution {
                        donor_1based: 20,
                        acceptor_1based: 30,
                    },
                ],
                ..RnaReadAlignmentInspectionRow::default()
            },
        ],
        ..RnaReadAlignmentInspection::default()
    };

    let exon = MainAreaDna::collect_rna_read_mapped_exon_contributors(&inspection);
    let junction = MainAreaDna::collect_rna_read_mapped_junction_contributors(&inspection);
    let isoform = MainAreaDna::collect_rna_read_mapped_isoform_contributors(&inspection);

    assert_eq!(exon.get(&(10, 20)).cloned().unwrap_or_default(), vec![0, 2]);
    assert_eq!(
        junction.get(&(20, 30)).cloned().unwrap_or_default(),
        vec![0, 2]
    );
    assert_eq!(isoform.get("TXA").cloned().unwrap_or_default(), vec![0, 2]);

    let dna = DNAsequence::from_sequence("ACGT").expect("sequence");
    let mut area = MainAreaDna::new(dna, None, None);
    area.focus_rna_read_alignment_effect_record_indices(
        isoform.get("TXA").cloned().unwrap_or_default(),
        "mapped isoform TXA",
    );
    assert_eq!(area.selected_rna_record_indices(), vec![0, 2]);
    assert_eq!(
        area.rna_read_alignment_effect_filter,
        super::RnaReadAlignmentEffectFilter::SelectedOnly
    );
    assert_eq!(
        area.rna_read_mapped_cdna_subview,
        super::RnaReadMappedCdnaSubview::ReadEffects
    );
}

#[test]
fn rna_read_isoform_triage_labels_and_colors_are_categorical() {
    assert_eq!(
        MainAreaDna::rna_read_isoform_triage_bin_label(Some(
            RnaReadIsoformTriageBin::KnownIsoformConfirmed
        )),
        "confirmed"
    );
    assert_eq!(
        MainAreaDna::rna_read_isoform_triage_bin_label(Some(
            RnaReadIsoformTriageBin::KnownIsoformAmbiguous
        )),
        "ambiguous"
    );
    assert_eq!(
        MainAreaDna::rna_read_isoform_triage_bin_label(None),
        "no triage"
    );
    assert_ne!(
        MainAreaDna::rna_read_isoform_triage_bin_color(Some(
            RnaReadIsoformTriageBin::KnownIsoformConfirmed
        )),
        MainAreaDna::rna_read_isoform_triage_bin_color(Some(
            RnaReadIsoformTriageBin::KnownIsoformAmbiguous
        ))
    );
}

#[test]
fn rna_read_isoform_triage_segment_fractions_follow_priority_order() {
    let counts = BTreeMap::from([
        (
            RnaReadIsoformTriageBin::KnownIsoformAmbiguous
                .as_str()
                .to_string(),
            2usize,
        ),
        (
            RnaReadIsoformTriageBin::KnownIsoformConfirmed
                .as_str()
                .to_string(),
            1usize,
        ),
        (
            RnaReadIsoformTriageBin::OffTargetOrBadSeed
                .as_str()
                .to_string(),
            1usize,
        ),
    ]);
    let segments = MainAreaDna::rna_read_isoform_triage_segment_fractions(&counts, 4);
    assert_eq!(
        segments
            .iter()
            .map(|(bin, count, _fraction)| (*bin, *count))
            .collect::<Vec<_>>(),
        vec![
            (RnaReadIsoformTriageBin::KnownIsoformConfirmed, 1),
            (RnaReadIsoformTriageBin::KnownIsoformAmbiguous, 2),
            (RnaReadIsoformTriageBin::OffTargetOrBadSeed, 1),
        ]
    );
    assert!((segments[1].2 - 0.5).abs() < 0.001);
}

#[test]
fn splicing_isoform_read_support_empty_text_distinguishes_missing_inputs() {
    assert_eq!(
        MainAreaDna::splicing_isoform_read_support_empty_text(false, false),
        "Select a saved RNA-read report to inspect isoform-level read support."
    );
    assert_eq!(
        MainAreaDna::splicing_isoform_read_support_empty_text(true, false),
        "No mapped isoform support rows are available yet; run or rerun phase-2 alignment in RNA-read Mapping."
    );
    assert_eq!(
        MainAreaDna::splicing_isoform_read_support_empty_text(true, true),
        ""
    );
}

#[test]
fn thresholded_cdna_exon_support_rows_follow_seed_passed_isoform_assignments() {
    let view = SplicingExpertView {
        seq_id: "seq1".to_string(),
        target_feature_id: 0,
        scope: SplicingScopePreset::AllOverlappingAnyStrand,
        group_label: "GENE1".to_string(),
        strand: "+".to_string(),
        region_start_1based: 1,
        region_end_1based: 60,
        transcript_count: 2,
        unique_exon_count: 3,
        instruction: String::new(),
        transcripts: vec![
            SplicingTranscriptLane {
                transcript_feature_id: 11,
                transcript_id: "TX1".to_string(),
                label: "TX1".to_string(),
                strand: "+".to_string(),
                exons: vec![
                    SplicingRange {
                        start_1based: 1,
                        end_1based: 10,
                    },
                    SplicingRange {
                        start_1based: 21,
                        end_1based: 30,
                    },
                ],
                exon_cds_phases: vec![],
                introns: vec![],
                has_target_feature: true,
            },
            SplicingTranscriptLane {
                transcript_feature_id: 12,
                transcript_id: "TX2".to_string(),
                label: "TX2".to_string(),
                strand: "+".to_string(),
                exons: vec![
                    SplicingRange {
                        start_1based: 1,
                        end_1based: 10,
                    },
                    SplicingRange {
                        start_1based: 41,
                        end_1based: 50,
                    },
                ],
                exon_cds_phases: vec![],
                introns: vec![],
                has_target_feature: false,
            },
        ],
        unique_exons: vec![
            SplicingExonSummary {
                start_1based: 1,
                end_1based: 10,
                support_transcript_count: 2,
                constitutive: true,
            },
            SplicingExonSummary {
                start_1based: 21,
                end_1based: 30,
                support_transcript_count: 1,
                constitutive: false,
            },
            SplicingExonSummary {
                start_1based: 41,
                end_1based: 50,
                support_transcript_count: 1,
                constitutive: false,
            },
        ],
        matrix_rows: vec![],
        boundaries: vec![],
        intron_signals: vec![],
        junctions: vec![SplicingJunctionArc {
            donor_1based: 10,
            acceptor_1based: 21,
            support_transcript_count: 1,
            transcript_feature_ids: vec![11],
        }],
        events: vec![],
    };
    let progress = RnaReadInterpretProgress {
        seq_id: "seq1".to_string(),
        reads_processed: 5,
        reads_total: 5,
        read_bases_processed: 0,
        mean_read_length_bp: 0.0,
        median_read_length_bp: 0,
        p95_read_length_bp: 0,
        input_bytes_processed: 0,
        input_bytes_total: 0,
        seed_passed: 5,
        aligned: 0,
        tested_kmers: 0,
        matched_kmers: 0,
        seed_compute_ms: 0.0,
        align_compute_ms: 0.0,
        io_read_ms: 0.0,
        fasta_parse_ms: 0.0,
        normalize_compute_ms: 0.0,
        inference_compute_ms: 0.0,
        progress_emit_ms: 0.0,
        update_every_reads: 1000,
        done: true,
        bins: vec![],
        score_density_bins: vec![],
        seed_pass_score_density_bins: vec![],
        top_hits_preview: vec![],
        transition_support_rows: vec![],
        isoform_support_rows: vec![
            RnaReadIsoformSupportRow {
                transcript_feature_id: 11,
                transcript_id: "TX1".to_string(),
                transcript_label: "TX1".to_string(),
                strand: "+".to_string(),
                exon_count: 2,
                expected_transition_count: 1,
                reads_assigned: 3,
                reads_seed_passed: 3,
                ..Default::default()
            },
            RnaReadIsoformSupportRow {
                transcript_feature_id: 12,
                transcript_id: "TX2".to_string(),
                transcript_label: "TX2".to_string(),
                strand: "+".to_string(),
                exon_count: 2,
                expected_transition_count: 1,
                reads_assigned: 2,
                reads_seed_passed: 2,
                ..Default::default()
            },
        ],
        mapped_exon_support_frequencies: vec![],
        mapped_junction_support_frequencies: vec![],
        mapped_isoform_support_rows: vec![],
        reads_with_transition_support: 0,
        transition_confirmations: 0,
        junction_crossing_seed_bits_indexed: 0,
        origin_class_counts: std::collections::BTreeMap::new(),
    };

    let rows = MainAreaDna::collect_thresholded_cdna_exon_support_rows(
        &view,
        &progress.isoform_support_rows,
        progress.seed_passed,
    );
    assert_eq!(rows.len(), 3);
    assert_eq!(rows[0].support_read_count, 5);
    assert_eq!(rows[1].support_read_count, 3);
    assert_eq!(rows[2].support_read_count, 2);
    assert!((rows[0].support_fraction - 1.0).abs() < f64::EPSILON);
    assert!((rows[1].support_fraction - 0.6).abs() < 1e-9);
    assert!((rows[2].support_fraction - 0.4).abs() < 1e-9);
}

#[test]
fn select_rna_read_report_max_score_record_indices_returns_all_ties() {
    let report = RnaReadInterpretationReport {
        report_id: "rna_reads_top".to_string(),
        hits: vec![
            RnaReadInterpretationHit {
                record_index: 0,
                seed_hit_fraction: 0.95,
                ..RnaReadInterpretationHit::default()
            },
            RnaReadInterpretationHit {
                record_index: 3,
                seed_hit_fraction: 1.0,
                ..RnaReadInterpretationHit::default()
            },
            RnaReadInterpretationHit {
                record_index: 8,
                seed_hit_fraction: 1.0,
                ..RnaReadInterpretationHit::default()
            },
        ],
        ..RnaReadInterpretationReport::default()
    };
    let indices = MainAreaDna::select_rna_read_report_max_score_record_indices(&report);
    assert_eq!(indices, vec![3, 8]);
}

#[test]
fn select_rna_read_report_rightmost_score_bin_record_indices_uses_saved_report_hits() {
    let report = RnaReadInterpretationReport {
        report_id: "rna_reads_bins".to_string(),
        score_density_bins: vec![0; 40],
        hits: vec![
            RnaReadInterpretationHit {
                record_index: 1,
                seed_hit_fraction: 0.32,
                ..RnaReadInterpretationHit::default()
            },
            RnaReadInterpretationHit {
                record_index: 5,
                seed_hit_fraction: 0.98,
                ..RnaReadInterpretationHit::default()
            },
            RnaReadInterpretationHit {
                record_index: 7,
                seed_hit_fraction: 0.99,
                ..RnaReadInterpretationHit::default()
            },
        ],
        ..RnaReadInterpretationReport::default()
    };
    let indices = MainAreaDna::select_rna_read_report_rightmost_score_bin_record_indices(
        &report,
        RnaReadScoreDensityVariant::AllScored,
        None,
    );
    assert_eq!(indices, vec![5, 7]);
}

#[test]
fn primer_backend_controls_default_when_missing_in_serialized_engine_ops_state() {
    let dna = DNAsequence::from_sequence("ACGT").unwrap();
    let area = MainAreaDna::new(dna, None, None);
    let mut value = serde_json::to_value(area.current_engine_ops_state()).unwrap();
    value.as_object_mut().unwrap().remove("primer_backend");
    value.as_object_mut().unwrap().remove("primer3_executable");
    let decoded: super::EngineOpsUiState = serde_json::from_value(value).unwrap();
    assert_eq!(decoded.primer_backend, PrimerDesignBackend::Auto);
    assert_eq!(decoded.primer3_executable, "primer3_core");
}

#[test]
fn apply_engine_ops_state_migrates_legacy_off_grouping_to_auto() {
    let dna = DNAsequence::from_sequence("ACGT").unwrap();
    let mut area = MainAreaDna::new(dna, None, None);
    let mut value = serde_json::to_value(area.current_engine_ops_state()).unwrap();
    value.as_object_mut().unwrap().insert(
        "feature_tree_grouping_mode".to_string(),
        serde_json::json!("off"),
    );
    value.as_object_mut().unwrap().insert(
        "feature_tree_second_level_grouping".to_string(),
        serde_json::json!(false),
    );
    value
        .as_object_mut()
        .unwrap()
        .remove("feature_tree_grouping_mode_explicit");
    let decoded: super::EngineOpsUiState = serde_json::from_value(value).unwrap();
    area.apply_engine_ops_state(decoded);
    assert_eq!(
        area.feature_tree_grouping_mode,
        super::FeatureTreeGroupingMode::Auto
    );
}

#[test]
fn apply_engine_ops_state_preserves_explicit_off_grouping() {
    let dna = DNAsequence::from_sequence("ACGT").unwrap();
    let mut area = MainAreaDna::new(dna, None, None);
    let mut value = serde_json::to_value(area.current_engine_ops_state()).unwrap();
    value.as_object_mut().unwrap().insert(
        "feature_tree_grouping_mode".to_string(),
        serde_json::json!("off"),
    );
    value.as_object_mut().unwrap().insert(
        "feature_tree_grouping_mode_explicit".to_string(),
        serde_json::json!(true),
    );
    let decoded: super::EngineOpsUiState = serde_json::from_value(value).unwrap();
    area.apply_engine_ops_state(decoded);
    assert_eq!(
        area.feature_tree_grouping_mode,
        super::FeatureTreeGroupingMode::Off
    );
}

#[test]
fn new_area_syncs_primer_backend_controls_from_engine_parameters() {
    let mut state = ProjectState::default();
    state.sequences.insert(
        "tpl".to_string(),
        DNAsequence::from_sequence("ACGTACGT").unwrap(),
    );
    let mut engine = GentleEngine::from_state(state);
    engine.state_mut().parameters.primer_design_backend = PrimerDesignBackend::Primer3;
    engine.state_mut().parameters.primer3_executable = "/tmp/fake_primer3_core".to_string();
    let engine = Arc::new(RwLock::new(engine));

    let area = MainAreaDna::new(
        DNAsequence::from_sequence("ACGTACGT").unwrap(),
        Some("tpl".to_string()),
        Some(engine),
    );
    assert_eq!(area.primer_backend, PrimerDesignBackend::Primer3);
    assert_eq!(area.primer3_executable, "/tmp/fake_primer3_core");
}

#[test]
fn primer_backend_help_summary_explains_auto_fallback_and_primer3_requirement() {
    let dna = DNAsequence::from_sequence("ACGT").unwrap();
    let mut area = MainAreaDna::new(dna, None, None);

    area.primer_backend = PrimerDesignBackend::Auto;
    assert!(area.primer_backend_help_summary().contains("falls back"));

    area.primer_backend = PrimerDesignBackend::Internal;
    assert!(
        area.primer_backend_help_summary()
            .contains("does not require Primer3")
    );

    area.primer_backend = PrimerDesignBackend::Primer3;
    assert!(
        area.primer_backend_help_summary()
            .contains("design requests fail")
    );
}

#[test]
fn primer_design_progress_summary_mentions_pair_evaluation_counts() {
    let summary = MainAreaDna::primer_design_progress_summary(&PrimerDesignProgress {
        seq_id: "tpl".to_string(),
        design_kind: "primer_pairs".to_string(),
        backend_requested: "internal".to_string(),
        backend_used: "internal".to_string(),
        stage: "pair_search".to_string(),
        detail: "Evaluated 42 pair candidate combination(s)".to_string(),
        roi_start_0based: 40,
        roi_end_0based_exclusive: 80,
        forward_candidate_count: Some(8),
        reverse_candidate_count: Some(7),
        probe_candidate_count: None,
        pair_candidate_combinations: Some(56),
        pair_evaluated: Some(42),
        pair_evaluation_limit: Some(100),
        pair_evaluation_limited: Some(false),
        accepted_pair_count: Some(3),
        assay_candidate_combinations: None,
        assays_evaluated: None,
        accepted_assay_count: None,
        max_output: 10,
        done: false,
    });

    assert!(summary.contains("checked 42/100 pair combinations"));
    assert!(summary.contains("accepted 3 candidate pairs"));
}

#[test]
fn primer_design_progress_summary_mentions_qpcr_probe_evaluation_counts() {
    let summary = MainAreaDna::primer_design_progress_summary(&PrimerDesignProgress {
        seq_id: "tpl".to_string(),
        design_kind: "qpcr_assays".to_string(),
        backend_requested: "internal".to_string(),
        backend_used: "internal".to_string(),
        stage: "assay_search".to_string(),
        detail: "Evaluated 12 assay combination(s)".to_string(),
        roi_start_0based: 40,
        roi_end_0based_exclusive: 120,
        forward_candidate_count: Some(6),
        reverse_candidate_count: Some(4),
        probe_candidate_count: Some(4),
        pair_candidate_combinations: Some(24),
        pair_evaluated: Some(6),
        pair_evaluation_limit: Some(40),
        pair_evaluation_limited: Some(false),
        accepted_pair_count: Some(6),
        assay_candidate_combinations: Some(24),
        assays_evaluated: Some(12),
        accepted_assay_count: Some(2),
        max_output: 5,
        done: false,
    });

    assert!(summary.contains("screened 12/24 probe placements"));
    assert!(summary.contains("accepted 2 assays"));
}

#[test]
fn bounded_center_window_keeps_requested_span_near_edges() {
    let window = MainAreaDna::bounded_center_window(80_000, 120, 500).expect("window");
    assert_eq!(window, (0, 1001));
    let window = MainAreaDna::bounded_center_window(80_000, 79_900, 500).expect("window");
    assert_eq!(window, (78_999, 80_000));
}

#[test]
fn bounded_center_window_returns_none_for_empty_sequence() {
    assert!(MainAreaDna::bounded_center_window(0, 0, 500).is_none());
}

#[test]
fn dotplot_crosshair_selection_bounds_spans_inclusive_interval() {
    assert_eq!(
        MainAreaDna::dotplot_crosshair_selection_bounds(12, 20, 100),
        Some((12, 21))
    );
    assert_eq!(
        MainAreaDna::dotplot_crosshair_selection_bounds(20, 12, 100),
        Some((12, 21))
    );
    assert_eq!(
        MainAreaDna::dotplot_crosshair_selection_bounds(99, 99, 100),
        Some((99, 100))
    );
}

#[test]
fn dotplot_crosshair_selection_bounds_handles_empty_sequence() {
    assert_eq!(
        MainAreaDna::dotplot_crosshair_selection_bounds(0, 0, 0),
        None
    );
}

#[test]
fn dotplot_crosshair_selection_bounds_for_pair_mode_uses_query_axis() {
    assert_eq!(
        MainAreaDna::dotplot_crosshair_selection_bounds_for_mode(
            DotplotMode::PairForward,
            12,
            80,
            100,
        ),
        Some((12, 13))
    );
    assert_eq!(
        MainAreaDna::dotplot_crosshair_selection_bounds_for_mode(
            DotplotMode::PairReverseComplement,
            99,
            1,
            100,
        ),
        Some((99, 100))
    );
}

#[test]
fn dotplot_reference_hit_envelope_applies_padding_and_clamp() {
    let mut view = DotplotView::default();
    view.reference_span_start_0based = 1_000;
    view.reference_span_end_0based = 2_000;
    view.points = vec![
        crate::engine::DotplotMatchPoint {
            x_0based: 20,
            y_0based: 1_120,
            mismatches: 0,
        },
        crate::engine::DotplotMatchPoint {
            x_0based: 40,
            y_0based: 1_360,
            mismatches: 1,
        },
    ];
    assert_eq!(
        MainAreaDna::dotplot_reference_hit_envelope(&view, 50),
        Some((1_070, 1_411))
    );
}

#[test]
fn dotplot_reference_hit_envelope_uses_all_overlay_series_points() {
    let mut view = DotplotView::default();
    view.reference_span_start_0based = 1_000;
    view.reference_span_end_0based = 2_000;
    view.series_count = 2;
    view.query_series = vec![
        crate::engine::DotplotQuerySeries {
            series_id: "overlay.series1".to_string(),
            seq_id: "iso_a".to_string(),
            label: "Isoform A".to_string(),
            transcript_feature_id: None,
            query_anchor_0based: None,
            query_anchor_label: None,
            color_rgb: [29, 78, 216],
            mode: DotplotMode::PairForward,
            span_start_0based: 0,
            span_end_0based: 200,
            point_count: 1,
            points: vec![crate::engine::DotplotMatchPoint {
                x_0based: 20,
                y_0based: 1_120,
                mismatches: 0,
            }],
            boxplot_bin_count: 0,
            boxplot_bins: vec![],
        },
        crate::engine::DotplotQuerySeries {
            series_id: "overlay.series2".to_string(),
            seq_id: "iso_b".to_string(),
            label: "Isoform B".to_string(),
            transcript_feature_id: None,
            query_anchor_0based: None,
            query_anchor_label: None,
            color_rgb: [220, 38, 38],
            mode: DotplotMode::PairForward,
            span_start_0based: 0,
            span_end_0based: 200,
            point_count: 1,
            points: vec![crate::engine::DotplotMatchPoint {
                x_0based: 40,
                y_0based: 1_360,
                mismatches: 0,
            }],
            boxplot_bin_count: 0,
            boxplot_bins: vec![],
        },
    ];
    assert_eq!(
        MainAreaDna::dotplot_reference_hit_envelope(&view, 50),
        Some((1_070, 1_411))
    );
}

#[test]
fn dotplot_reference_hit_envelope_returns_none_without_points() {
    let mut view = DotplotView::default();
    view.reference_span_start_0based = 1_000;
    view.reference_span_end_0based = 2_000;
    assert_eq!(
        MainAreaDna::dotplot_reference_hit_envelope(&view, 100),
        None
    );
}

#[test]
fn auto_fit_reference_span_for_view_if_requested_returns_fit_span() {
    let mut view = DotplotView::default();
    view.mode = DotplotMode::PairForward;
    view.reference_span_start_0based = 1_000;
    view.reference_span_end_0based = 2_000;
    view.points = vec![
        crate::engine::DotplotMatchPoint {
            x_0based: 20,
            y_0based: 1_120,
            mismatches: 0,
        },
        crate::engine::DotplotMatchPoint {
            x_0based: 40,
            y_0based: 1_360,
            mismatches: 0,
        },
    ];
    assert_eq!(
        MainAreaDna::auto_fit_reference_span_for_view_if_requested(true, &view, 50),
        Some((1_070, 1_411))
    );
}

#[test]
fn auto_fit_reference_span_for_view_if_requested_skips_without_request() {
    let mut view = DotplotView::default();
    view.mode = DotplotMode::PairForward;
    view.reference_span_start_0based = 1_000;
    view.reference_span_end_0based = 2_000;
    view.points = vec![crate::engine::DotplotMatchPoint {
        x_0based: 20,
        y_0based: 1_120,
        mismatches: 0,
    }];
    assert_eq!(
        MainAreaDna::auto_fit_reference_span_for_view_if_requested(false, &view, 50),
        None
    );
}

#[test]
fn splicing_lane_index_at_y_returns_expected_lane() {
    assert_eq!(
        MainAreaDna::splicing_lane_index_at_y(100.0, 80.0, 20.0, 4),
        Some(1)
    );
    assert_eq!(
        MainAreaDna::splicing_lane_index_at_y(159.9, 80.0, 20.0, 4),
        Some(3)
    );
}

#[test]
fn splicing_lane_index_at_y_rejects_out_of_range_positions() {
    assert_eq!(
        MainAreaDna::splicing_lane_index_at_y(79.0, 80.0, 20.0, 4),
        None
    );
    assert_eq!(
        MainAreaDna::splicing_lane_index_at_y(160.0, 80.0, 20.0, 4),
        None
    );
    assert_eq!(
        MainAreaDna::splicing_lane_index_at_y(100.0, 80.0, 0.0, 4),
        None
    );
}

#[test]
fn large_splicing_matrix_defaults_to_collapsed_section() {
    assert!(!MainAreaDna::splicing_matrix_should_default_collapsed(
        20, 200
    ));
    assert!(MainAreaDna::splicing_matrix_should_default_collapsed(
        100, 200
    ));
}

#[test]
fn large_splicing_transition_matrix_defaults_to_collapsed_section() {
    assert!(!MainAreaDna::splicing_transition_should_default_collapsed(
        80
    ));
    assert!(MainAreaDna::splicing_transition_should_default_collapsed(
        81
    ));
}

#[test]
fn large_splicing_attract_section_defaults_to_collapsed() {
    let small_view = SplicingExpertView {
        seq_id: "seq1".to_string(),
        target_feature_id: 1,
        scope: SplicingScopePreset::AllOverlappingAnyStrand,
        group_label: "TP73".to_string(),
        strand: "+".to_string(),
        region_start_1based: 1,
        region_end_1based: 100,
        transcript_count: 0,
        unique_exon_count: 0,
        instruction: String::new(),
        transcripts: vec![],
        unique_exons: vec![],
        matrix_rows: vec![],
        boundaries: vec![],
        intron_signals: vec![],
        junctions: vec![],
        events: vec![],
    };
    let mut large_view = small_view.clone();
    large_view.boundaries = (0..(SPLICING_ATTRACT_EAGER_BOUNDARY_THRESHOLD + 1))
        .map(|idx| SplicingBoundaryMarker {
            transcript_feature_id: idx,
            transcript_id: format!("tx{idx}"),
            side: "donor".to_string(),
            position_1based: idx + 1,
            motif_2bp: "GT".to_string(),
            canonical: true,
            canonical_pair: true,
            partner_position_1based: idx + 2,
            paired_motif_signature: "GT-AG".to_string(),
            motif_class: "gt_ag_major_canonical".to_string(),
            annotation: String::new(),
        })
        .collect();

    assert!(!MainAreaDna::splicing_attract_should_default_collapsed(
        &small_view
    ));
    assert!(MainAreaDna::splicing_attract_should_default_collapsed(
        &large_view
    ));
}

#[test]
fn splicing_map_placeholder_svg_mentions_selected_sequence() {
    let svg = MainAreaDna::render_splicing_map_placeholder_svg("tp53.seq");
    assert!(svg.contains("tp53.seq"));
    assert!(svg.contains("No splicing-seed feature is selected"));
    assert!(svg.contains("<svg"));
}

#[test]
fn refresh_description_cache_builds_splicing_view_for_gene_seed_feature() {
    let mut dna =
        DNAsequence::from_sequence("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA").expect("sequence");
    dna.features_mut().push(Feature {
        kind: "mRNA".into(),
        location: Location::Join(vec![
            Location::simple_range(2, 8),
            Location::simple_range(12, 20),
            Location::simple_range(26, 34),
        ]),
        qualifiers: vec![
            ("gene".into(), Some("GENE1".to_string())),
            ("transcript_id".into(), Some("NM_TEST_1".to_string())),
            ("label".into(), Some("NM_TEST_1".to_string())),
        ],
    });
    dna.features_mut().push(Feature {
        kind: "gene".into(),
        location: Location::simple_range(2, 34),
        qualifiers: vec![
            ("gene".into(), Some("GENE1".to_string())),
            ("label".into(), Some("GENE1".to_string())),
        ],
    });
    let mut state = ProjectState::default();
    state.sequences.insert("seq_gene".to_string(), dna.clone());
    let engine = Arc::new(RwLock::new(GentleEngine::from_state(state)));
    let mut area = MainAreaDna::new(dna, Some("seq_gene".to_string()), Some(engine));

    area.focus_feature(1);
    area.refresh_description_cache();

    match area.description_cache_expert_view.as_ref() {
        Some(FeatureExpertView::Splicing(view)) => {
            assert_eq!(view.target_feature_id, 1);
            assert_eq!(view.group_label, "GENE1");
            assert_eq!(view.transcript_count, 1);
        }
        other => panic!("expected splicing expert view for gene seed, got {other:?}"),
    }
}

#[test]
fn refresh_description_cache_builds_splicing_view_for_misc_rna_seed_feature() {
    let mut dna =
        DNAsequence::from_sequence("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA").expect("sequence");
    dna.features_mut().push(Feature {
        kind: "misc_RNA".into(),
        location: Location::Join(vec![
            Location::simple_range(2, 8),
            Location::simple_range(12, 20),
            Location::simple_range(26, 34),
        ]),
        qualifiers: vec![
            ("gene".into(), Some("GFOD3P".to_string())),
            ("transcript_id".into(), Some("NR_TEST_MISC_1".to_string())),
            ("label".into(), Some("GFOD3P misc_RNA".to_string())),
        ],
    });
    let mut state = ProjectState::default();
    state.sequences.insert("seq_misc".to_string(), dna.clone());
    let engine = Arc::new(RwLock::new(GentleEngine::from_state(state)));
    let mut area = MainAreaDna::new(dna, Some("seq_misc".to_string()), Some(engine));

    area.focus_feature(0);
    area.refresh_description_cache();

    match area.description_cache_expert_view.as_ref() {
        Some(FeatureExpertView::Splicing(view)) => {
            assert_eq!(view.target_feature_id, 0);
            assert_eq!(view.group_label, "GFOD3P");
            assert_eq!(view.transcript_count, 1);
        }
        other => panic!("expected splicing expert view for misc_RNA seed, got {other:?}"),
    }
}

#[test]
fn refresh_description_cache_does_not_auto_open_splicing_window() {
    let mut dna =
        DNAsequence::from_sequence("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA").expect("sequence");
    dna.features_mut().push(Feature {
        kind: "mRNA".into(),
        location: Location::Join(vec![
            Location::simple_range(2, 8),
            Location::simple_range(12, 20),
            Location::simple_range(26, 34),
        ]),
        qualifiers: vec![
            ("gene".into(), Some("GENE1".to_string())),
            ("transcript_id".into(), Some("NM_TEST_1".to_string())),
            ("label".into(), Some("NM_TEST_1".to_string())),
        ],
    });
    let mut state = ProjectState::default();
    state.sequences.insert("seq_gene".to_string(), dna.clone());
    let engine = Arc::new(RwLock::new(GentleEngine::from_state(state)));
    let mut area = MainAreaDna::new(dna, Some("seq_gene".to_string()), Some(engine));
    area.focus_feature(0);

    area.refresh_description_cache();

    assert!(
        !area.show_splicing_expert_window,
        "single-click/selection should not auto-open Splicing Expert"
    );
}

#[test]
fn refresh_description_cache_uses_construct_reasoning_selection_details() {
    let dna = DNAsequence::from_sequence("ACGTACGTACGTACGTACGT").expect("sequence");
    let mut area = MainAreaDna::new(dna, Some("seq_reasoning".to_string()), None);
    area.dna_display
        .write()
        .expect("display")
        .set_construct_reasoning_overlay(Some(ConstructReasoningOverlay {
            graph_id: "graph_reasoning".to_string(),
            seq_id: "seq_reasoning".to_string(),
            objective_title: "Assemble promoter cassette".to_string(),
            objective_goal: "Inspect promoter evidence".to_string(),
            evidence: vec![ConstructReasoningOverlaySpan {
                annotation_id: "annotation_promoter".to_string(),
                evidence_id: "ev_promoter".to_string(),
                summary_id: "summary_promoter".to_string(),
                summary_title: "Promoter candidate".to_string(),
                summary_subtitle: "Transcript interpretations disagree".to_string(),
                summary_candidate_count: 2,
                summary_review_status: "draft".to_string(),
                start_0based: 3,
                end_0based_exclusive: 12,
                strand: Some("+".to_string()),
                role: ConstructRole::Promoter,
                evidence_class: EvidenceClass::ReliableAnnotation,
                label: "Promoter candidate".to_string(),
                rationale: "Annotation and restriction context agree".to_string(),
                score: Some(0.75),
                confidence: Some(0.91),
                context_tags: vec!["cassette".to_string(), "5prime".to_string()],
                provenance_kind: "annotation_projection".to_string(),
                provenance_refs: vec!["gene:TP73".to_string()],
                source_kind: "supporting_annotation".to_string(),
                supporting_fact_labels: vec!["Variant effect candidates derived".to_string()],
                supporting_decision_titles: vec!["Evaluate Variant Effect Context".to_string()],
                transcript_context_status: Some("multi_transcript_ambiguous".to_string()),
                effect_tags: vec![
                    "promoter_variant_candidate".to_string(),
                    "transcript_context_ambiguous".to_string(),
                ],
                editable_status: EditableStatus::Draft,
                warnings: vec![],
                notes: vec!["inspect upstream fusion".to_string()],
            }],
        }));
    area.map_dna
        .select_reasoning_evidence(Some("ev_promoter".to_string()));

    area.refresh_description_cache();

    assert_eq!(area.description_cache_title, "Promoter candidate");
    assert_eq!(
        area.description_cache_range.as_deref(),
        Some("4..12 (9 bp, strand +)")
    );
    assert!(
        area.description_cache_details
            .iter()
            .any(|line| line.contains("role: Promoter"))
    );
    assert!(
        area.description_cache_details
            .iter()
            .any(|line| line.contains("why: Annotation and restriction context agree"))
    );
    assert!(
        area.description_cache_details
            .iter()
            .any(|line| line.contains("supports_facts: Variant effect candidates derived"))
    );
    assert!(
        area.description_cache_details
            .iter()
            .any(|line| line.contains("transcript_context: multi_transcript_ambiguous"))
    );
    assert!(
        area.description_cache_details
            .iter()
            .any(|line| line.contains("summary: Transcript interpretations disagree"))
    );
    assert!(area.description_cache_expert_view.is_none());
}

#[test]
fn refresh_description_cache_includes_non_sequence_construct_reasoning_context() {
    let dna = DNAsequence::from_sequence("ATGCGTATGCGTATGCGTATGCGT").expect("sequence");
    let mut state = ProjectState::default();
    state
        .sequences
        .insert("seq_reasoning_context".to_string(), dna.clone());
    let mut engine = GentleEngine::from_state(state);
    let objective = engine
        .upsert_construct_objective(crate::engine::ConstructObjective {
            title: "GUI reasoning".to_string(),
            goal: "Show host/helper/growth context in the inspector".to_string(),
            propagation_host_profile_id: Some("ecoli_k12".to_string()),
            expression_host_profile_id: Some("ecoli_k12".to_string()),
            helper_profile_id: Some("pUC19".to_string()),
            medium_conditions: vec!["ampicillin".to_string()],
            ..crate::engine::ConstructObjective::default()
        })
        .expect("objective");
    engine
        .build_construct_reasoning_graph(
            "seq_reasoning_context",
            Some(&objective.objective_id),
            None,
        )
        .expect("build graph");
    let engine = Arc::new(RwLock::new(engine));
    let mut area = MainAreaDna::new(dna, Some("seq_reasoning_context".to_string()), Some(engine));

    area.refresh_description_cache();

    let reasoning = area
        .description_cache_construct_reasoning
        .as_ref()
        .expect("construct reasoning cache");
    assert_eq!(reasoning.objective_title, "GUI reasoning");
    assert!(
        reasoning
            .fact_entries
            .iter()
            .any(|entry| entry.title == "Selection/complementation context supported")
    );
    assert!(reasoning.fact_entries.iter().any(|entry| {
        entry
            .detail_lines
            .iter()
            .any(|line| line.contains("ampicillin_selection"))
    }));
    assert!(
        reasoning
            .decision_entries
            .iter()
            .any(|entry| entry.title == "Evaluate Selection/Complementation Fit")
    );
}

#[test]
fn refresh_description_cache_includes_annotation_summary_entries() {
    let mut dna = DNAsequence::from_sequence("ATGCGTATGCGT").expect("sequence");
    dna.features_mut().push(gb_io::seq::Feature {
        kind: "exon".into(),
        location: gb_io::seq::Location::simple_range(0, 10),
        qualifiers: vec![
            ("label".into(), Some("TP73 exon 1".to_string())),
            ("note".into(), Some("confirmed by cDNA mapping".to_string())),
        ],
    });
    dna.features_mut().push(gb_io::seq::Feature {
        kind: "exon".into(),
        location: gb_io::seq::Location::simple_range(2, 12),
        qualifiers: vec![
            ("label".into(), Some("TP73 exon 1 alt".to_string())),
            ("evidence".into(), Some("supported by cDNA".to_string())),
        ],
    });
    dna.update_computed_features();
    let mut state = ProjectState::default();
    state
        .sequences
        .insert("seq_reasoning_summary".to_string(), dna.clone());
    let mut engine = GentleEngine::from_state(state);
    let graph = engine
        .build_construct_reasoning_graph("seq_reasoning_summary", None, None)
        .expect("build graph");
    assert!(
        graph.annotation_candidate_summaries.iter().any(|summary| {
            summary.role == ConstructRole::Exon
                && summary.candidate_count >= 2
                && summary.subtitle.contains("overlapping")
        }),
        "graph summaries: {:?}",
        graph.annotation_candidate_summaries
    );
    let engine = Arc::new(RwLock::new(engine));
    let mut area = MainAreaDna::new(dna, Some("seq_reasoning_summary".to_string()), Some(engine));

    area.focus_construct_reasoning_graph(&graph.graph_id);
    area.refresh_description_cache();

    let reasoning = area
        .description_cache_construct_reasoning
        .as_ref()
        .expect("construct reasoning cache");
    assert!(
        reasoning.annotation_summary_entries.iter().any(|entry| {
            entry.title.contains("Exon")
                && entry
                    .detail_lines
                    .iter()
                    .any(|line| line.contains("collapsed_candidates: 2"))
                && entry
                    .detail_lines
                    .iter()
                    .any(|line| line.contains("overlapping exon candidates"))
        }),
        "summary entries: {:?}",
        reasoning.annotation_summary_entries
    );
}

#[test]
fn construct_reasoning_summary_entries_use_stable_ids_not_titles() {
    let graph = crate::engine::ConstructReasoningGraph::default();
    let first = MainAreaDna::construct_reasoning_inspector_entry_for_annotation_summary(
        &graph,
        &AnnotationCandidateSummary {
            summary_id: "summary-1".to_string(),
            title: "Inverted repeat cluster".to_string(),
            ..AnnotationCandidateSummary::default()
        },
    );
    let second = MainAreaDna::construct_reasoning_inspector_entry_for_annotation_summary(
        &graph,
        &AnnotationCandidateSummary {
            summary_id: "summary-2".to_string(),
            title: "Inverted repeat cluster".to_string(),
            ..AnnotationCandidateSummary::default()
        },
    );

    assert_eq!(first.title, second.title);
    assert_ne!(first.stable_id, second.stable_id);
    assert_eq!(first.stable_id, "summary-1");
    assert_eq!(second.stable_id, "summary-2");
}

#[test]
fn refresh_description_cache_includes_similarity_reasoning_facts() {
    let sequence = format!(
        "{}{}{}{}{}",
        "ACGT".repeat(12),
        "AAAAAAAAAAAAAA",
        "ATATATATATATATATATAT",
        "GATTACAGATTACCCGGGGATTACAGATTA",
        "GCGTACGCTATTTTTAGCGTACGC"
    );
    let dna = DNAsequence::from_sequence(&sequence).expect("sequence");
    let mut state = ProjectState::default();
    state
        .sequences
        .insert("seq_reasoning_similarity".to_string(), dna.clone());
    let mut engine = GentleEngine::from_state(state);
    let graph = engine
        .build_construct_reasoning_graph("seq_reasoning_similarity", None, None)
        .expect("build graph");
    let engine = Arc::new(RwLock::new(engine));
    let mut area = MainAreaDna::new(
        dna,
        Some("seq_reasoning_similarity".to_string()),
        Some(engine),
    );

    area.focus_construct_reasoning_graph(&graph.graph_id);
    area.refresh_description_cache();

    let reasoning = area
        .description_cache_construct_reasoning
        .as_ref()
        .expect("construct reasoning cache");
    assert!(reasoning.fact_entries.iter().any(|entry| {
        entry.title == "Low-complexity / tandem-repeat context detected"
            || entry.title == "Low-complexity context detected"
    }));
    assert!(
        reasoning
            .fact_entries
            .iter()
            .any(|entry| entry.title == "Similarity/low-complexity operational risks detected")
    );
    assert!(
        reasoning
            .fact_entries
            .iter()
            .any(|entry| entry.title == "PCR/amplification review suggested")
    );
    assert!(
        reasoning
            .fact_entries
            .iter()
            .any(|entry| entry.title == "Nanopore/direct-sequencing review suggested")
    );
    assert!(
        reasoning
            .fact_entries
            .iter()
            .any(|entry| entry.title == "Repeat-driven mapping review suggested")
    );
    assert!(
        reasoning
            .fact_entries
            .iter()
            .any(|entry| entry.title == "Cloning stability review suggested")
    );
}

#[test]
fn construct_reasoning_similarity_fact_entries_offer_dotplot_actions() {
    let sequence = format!(
        "{}{}{}{}{}",
        "ACGT".repeat(12),
        "AAAAAAAAAAAAAA",
        "ATATATATATATATATATAT",
        "GATTACAGATTACCCGGGGATTACAGATTA",
        "GCGTACGCTATTTTTAGCGTACGC"
    );
    let dna = DNAsequence::from_sequence(&sequence).expect("sequence");
    let mut state = ProjectState::default();
    state
        .sequences
        .insert("seq_reasoning_similarity".to_string(), dna.clone());
    let mut engine = GentleEngine::from_state(state);
    let graph = engine
        .build_construct_reasoning_graph("seq_reasoning_similarity", None, None)
        .expect("build graph");

    let cache = MainAreaDna::build_construct_reasoning_inspector_cache(&graph);

    let pcr_entry = cache
        .fact_entries
        .iter()
        .find(|entry| entry.title == "PCR/amplification review suggested")
        .expect("pcr reasoning entry");
    assert!(
        pcr_entry
            .dotplot_actions
            .iter()
            .any(|action| action.mode == DotplotMode::SelfForward),
        "pcr entry actions: {:?}",
        pcr_entry.dotplot_actions
    );

    let cloning_entry = cache
        .fact_entries
        .iter()
        .find(|entry| entry.title == "Cloning stability review suggested")
        .expect("cloning reasoning entry");
    assert!(
        cloning_entry
            .dotplot_actions
            .iter()
            .any(|action| action.mode == DotplotMode::SelfReverseComplement),
        "cloning entry actions: {:?}",
        cloning_entry.dotplot_actions
    );
}

#[test]
fn open_construct_reasoning_dotplot_action_focuses_repeat_region() {
    let sequence = format!(
        "{}{}{}{}{}",
        "ACGT".repeat(12),
        "AAAAAAAAAAAAAA",
        "ATATATATATATATATATAT",
        "GATTACAGATTACCCGGGGATTACAGATTA",
        "GCGTACGCTATTTTTAGCGTACGC"
    );
    let dna = DNAsequence::from_sequence(&sequence).expect("sequence");
    let mut state = ProjectState::default();
    state
        .sequences
        .insert("seq_reasoning_similarity".to_string(), dna.clone());
    let mut engine = GentleEngine::from_state(state);
    let graph = engine
        .build_construct_reasoning_graph("seq_reasoning_similarity", None, None)
        .expect("build graph");
    let engine = Arc::new(RwLock::new(engine));
    let mut area = MainAreaDna::new(
        dna,
        Some("seq_reasoning_similarity".to_string()),
        Some(engine),
    );
    let cache = MainAreaDna::build_construct_reasoning_inspector_cache(&graph);
    let action = cache
        .fact_entries
        .iter()
        .find(|entry| entry.title == "Cloning stability review suggested")
        .and_then(|entry| {
            entry
                .dotplot_actions
                .iter()
                .find(|action| action.mode == DotplotMode::SelfReverseComplement)
        })
        .cloned()
        .expect("cloning revcomp dotplot action");

    area.open_construct_reasoning_dotplot_action(&action)
        .expect("open reasoning-guided dotplot");

    assert!(area.show_dotplot_window);
    assert_eq!(area.primary_map_mode, PrimaryMapMode::Dotplot);
    assert_eq!(area.dotplot_ui.mode, DotplotMode::SelfReverseComplement);
    assert!(!area.dotplot_ui.dotplot_id.trim().is_empty());
    let view = area
        .dotplot_cached_view
        .as_ref()
        .expect("cached reasoning-guided dotplot");
    let focus_midpoint_0based = action.focus_start_0based.saturating_add(
        action
            .focus_end_0based_exclusive
            .saturating_sub(action.focus_start_0based)
            / 2,
    );
    assert!(view.span_start_0based <= focus_midpoint_0based);
    assert!(view.span_end_0based > focus_midpoint_0based);
}

#[test]
fn focus_construct_reasoning_graph_prefers_requested_graph() {
    let dna = DNAsequence::from_sequence("ATGCGTATGCGTATGCGTATGCGT").expect("sequence");
    let mut state = ProjectState::default();
    state
        .sequences
        .insert("seq_reasoning_focus".to_string(), dna.clone());
    let mut engine = GentleEngine::from_state(state);
    let objective_one = engine
        .upsert_construct_objective(crate::engine::ConstructObjective {
            objective_id: "obj_reasoning_one".to_string(),
            title: "Reasoning One".to_string(),
            goal: "Inspect first graph".to_string(),
            ..crate::engine::ConstructObjective::default()
        })
        .expect("objective one");
    let objective_two = engine
        .upsert_construct_objective(crate::engine::ConstructObjective {
            objective_id: "obj_reasoning_two".to_string(),
            title: "Reasoning Two".to_string(),
            goal: "Inspect second graph".to_string(),
            ..crate::engine::ConstructObjective::default()
        })
        .expect("objective two");
    engine
        .build_construct_reasoning_graph(
            "seq_reasoning_focus",
            Some(&objective_one.objective_id),
            Some("graph_reasoning_one"),
        )
        .expect("build graph one");
    engine
        .build_construct_reasoning_graph(
            "seq_reasoning_focus",
            Some(&objective_two.objective_id),
            Some("graph_reasoning_two"),
        )
        .expect("build graph two");
    let engine = Arc::new(RwLock::new(engine));
    let mut area = MainAreaDna::new(dna, Some("seq_reasoning_focus".to_string()), Some(engine));

    area.focus_construct_reasoning_graph("graph_reasoning_one");
    area.refresh_description_cache();

    let reasoning = area
        .description_cache_construct_reasoning
        .as_ref()
        .expect("construct reasoning cache");
    assert_eq!(reasoning.graph_id, "graph_reasoning_one");
    assert_eq!(reasoning.objective_title, "Reasoning One");

    area.focus_construct_reasoning_graph("graph_reasoning_two");
    area.refresh_description_cache();

    let reasoning = area
        .description_cache_construct_reasoning
        .as_ref()
        .expect("construct reasoning cache");
    assert_eq!(reasoning.graph_id, "graph_reasoning_two");
    assert_eq!(reasoning.objective_title, "Reasoning Two");
}

#[test]
fn set_construct_reasoning_annotation_candidate_status_refreshes_overlay_and_details() {
    let mut dna = DNAsequence::from_sequence("ATGCGTATGCGT").expect("sequence");
    dna.features_mut().push(Feature {
        kind: "exon".into(),
        location: Location::simple_range(2, 10),
        qualifiers: vec![
            ("label".into(), Some("Confirmed exon".to_string())),
            ("evidence".into(), Some("supported by cDNA".to_string())),
        ],
    });
    dna.update_computed_features();

    let mut state = ProjectState::default();
    state
        .sequences
        .insert("seq_reasoning_status".to_string(), dna.clone());
    let mut engine = GentleEngine::from_state(state);
    let graph = engine
        .build_construct_reasoning_graph(
            "seq_reasoning_status",
            None,
            Some("graph_reasoning_status"),
        )
        .expect("build graph");
    let candidate = graph
        .annotation_candidates
        .iter()
        .find(|candidate| candidate.role == ConstructRole::Exon)
        .expect("annotation candidate");
    let engine = Arc::new(RwLock::new(engine));
    let mut area = MainAreaDna::new(dna, Some("seq_reasoning_status".to_string()), Some(engine));

    area.focus_construct_reasoning_graph("graph_reasoning_status");
    area.map_dna
        .select_reasoning_evidence(Some(candidate.evidence_id.clone()));
    area.refresh_description_cache();

    area.set_construct_reasoning_annotation_candidate_status(
        "graph_reasoning_status",
        &candidate.annotation_id,
        EditableStatus::Accepted,
    );

    let (_, span) = area
        .find_construct_reasoning_span(&candidate.evidence_id)
        .expect("updated reasoning span");
    assert_eq!(span.editable_status, EditableStatus::Accepted);
    assert!(
        area.description_cache_details
            .iter()
            .any(|line| line.contains("editable_status: accepted"))
    );
    let reasoning = area
        .description_cache_construct_reasoning
        .as_ref()
        .expect("construct reasoning cache");
    assert!(reasoning.annotation_entries.iter().any(|entry| {
        entry.annotation_id.as_deref() == Some(candidate.annotation_id.as_str())
            && entry
                .detail_lines
                .iter()
                .any(|line| line.contains("status: accepted"))
    }));
}

#[test]
fn write_back_construct_reasoning_annotation_candidate_refreshes_sequence_and_status() {
    let mut dna = DNAsequence::from_sequence(&"A".repeat(6001)).expect("sequence");
    dna.features_mut().push(Feature {
        kind: "mRNA".into(),
        location: Location::Complement(Box::new(Location::simple_range(2000, 2613))),
        qualifiers: vec![
            ("gene".into(), Some("VKORC1".to_string())),
            ("transcript_id".into(), Some("ENSTVKORC1".to_string())),
            ("label".into(), Some("VKORC1-201".to_string())),
        ],
    });
    dna.update_computed_features();

    let mut state = ProjectState::default();
    state
        .sequences
        .insert("seq_reasoning_writeback".to_string(), dna.clone());
    let mut engine = GentleEngine::from_state(state);
    let graph = engine
        .build_construct_reasoning_graph(
            "seq_reasoning_writeback",
            None,
            Some("graph_reasoning_writeback"),
        )
        .expect("build graph");
    let candidate = graph
        .annotation_candidates
        .iter()
        .find(|candidate| candidate.role == ConstructRole::Promoter)
        .cloned()
        .expect("generated promoter annotation candidate");
    engine
        .set_construct_reasoning_annotation_candidate_status(
            &graph.graph_id,
            &candidate.annotation_id,
            EditableStatus::Accepted,
        )
        .expect("accept annotation candidate");
    let engine = Arc::new(RwLock::new(engine));
    let mut area = MainAreaDna::new(
        dna,
        Some("seq_reasoning_writeback".to_string()),
        Some(engine.clone()),
    );

    area.focus_construct_reasoning_graph("graph_reasoning_writeback");
    area.refresh_description_cache();
    area.write_back_construct_reasoning_annotation_candidate(
        "graph_reasoning_writeback",
        &candidate.annotation_id,
    );

    let stored_dna = engine
        .read()
        .expect("engine")
        .state()
        .sequences
        .get("seq_reasoning_writeback")
        .cloned()
        .expect("stored sequence");
    assert!(stored_dna.features().iter().any(|feature| {
        feature
            .qualifier_values("construct_reasoning_annotation_id")
            .any(|value| value == candidate.annotation_id.as_str())
    }));
    assert!(area.op_status.contains("Wrote back annotation candidate"));
    let reasoning = area
        .description_cache_construct_reasoning
        .as_ref()
        .expect("construct reasoning cache");
    assert!(reasoning.annotation_entries.iter().any(|entry| {
        entry
            .detail_lines
            .iter()
            .any(|line| line.contains("status: accepted"))
    }));
}

#[test]
fn sync_from_engine_display_does_not_rebuild_construct_reasoning_overlay_each_frame() {
    let dna = DNAsequence::from_sequence("ATGCGTATGCGTATGCGTATGCGT").expect("sequence");
    let mut state = ProjectState::default();
    state
        .sequences
        .insert("seq_reasoning_sync".to_string(), dna.clone());
    let engine = Arc::new(RwLock::new(GentleEngine::from_state(state)));
    let mut area = MainAreaDna::new(
        dna,
        Some("seq_reasoning_sync".to_string()),
        Some(engine.clone()),
    );

    let initial_generated_at = engine
        .read()
        .expect("engine")
        .construct_reasoning_graph_for_seq_id("seq_reasoning_sync")
        .expect("graph generated during initial sync")
        .generated_at_unix_ms;
    std::thread::sleep(Duration::from_millis(2));

    area.sync_from_engine_display();

    let repeated_generated_at = engine
        .read()
        .expect("engine")
        .construct_reasoning_graph_for_seq_id("seq_reasoning_sync")
        .expect("graph still present")
        .generated_at_unix_ms;

    assert_eq!(
        repeated_generated_at, initial_generated_at,
        "ordinary repaint sync should reuse the current construct-reasoning overlay instead of rebuilding it"
    );
}

#[test]
fn refresh_description_cache_includes_variant_reasoning_context() {
    let mut dna = DNAsequence::from_sequence("ATGGAATTTACGTACGT").expect("sequence");
    dna.features_mut().push(Feature {
        kind: "transcript".into(),
        location: Location::simple_range(0, 9),
        qualifiers: vec![("label".into(), Some("Demo transcript".to_string()))],
    });
    dna.features_mut().push(Feature {
        kind: "CDS".into(),
        location: Location::simple_range(0, 9),
        qualifiers: vec![("label".into(), Some("Demo CDS".to_string()))],
    });
    dna.features_mut().push(Feature {
        kind: "variation".into(),
        location: Location::simple_range(3, 4),
        qualifiers: vec![
            ("label".into(), Some("rsGui".to_string())),
            (
                "gentle_generated".into(),
                Some("genome_vcf_track".to_string()),
            ),
            ("vcf_ref".into(), Some("G".to_string())),
            ("vcf_alt".into(), Some("A".to_string())),
        ],
    });
    dna.update_computed_features();

    let mut state = ProjectState::default();
    state
        .sequences
        .insert("seq_reasoning_variant".to_string(), dna.clone());
    let mut engine = GentleEngine::from_state(state);
    engine
        .build_construct_reasoning_graph("seq_reasoning_variant", None, None)
        .expect("build graph");
    let engine = Arc::new(RwLock::new(engine));
    let mut area = MainAreaDna::new(dna, Some("seq_reasoning_variant".to_string()), Some(engine));

    area.refresh_description_cache();

    let reasoning = area
        .description_cache_construct_reasoning
        .as_ref()
        .expect("construct reasoning cache");
    assert!(reasoning.fact_entries.iter().any(|entry| {
        entry.title == "Variant effect candidates derived"
            && entry
                .detail_lines
                .iter()
                .any(|line| line.contains("coding_variant_candidate"))
            && entry
                .detail_lines
                .iter()
                .any(|line| line.contains("transcript_context"))
    }));
    assert!(reasoning.fact_entries.iter().any(|entry| {
        entry.title == "Promoter-design assays suggested"
            && entry
                .detail_lines
                .iter()
                .any(|line| line.contains("allele_paired_expression_compare"))
            && entry
                .detail_lines
                .iter()
                .any(|line| line.contains("per_variant"))
    }));
}

#[test]
fn open_splicing_expert_for_feature_opens_window_on_explicit_request() {
    let mut dna =
        DNAsequence::from_sequence("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA").expect("sequence");
    dna.features_mut().push(Feature {
        kind: "mRNA".into(),
        location: Location::Join(vec![
            Location::simple_range(2, 8),
            Location::simple_range(12, 20),
            Location::simple_range(26, 34),
        ]),
        qualifiers: vec![
            ("gene".into(), Some("GENE1".to_string())),
            ("transcript_id".into(), Some("NM_TEST_1".to_string())),
            ("label".into(), Some("NM_TEST_1".to_string())),
        ],
    });
    let mut state = ProjectState::default();
    state.sequences.insert("seq_gene".to_string(), dna.clone());
    let engine = Arc::new(RwLock::new(GentleEngine::from_state(state)));
    let mut area = MainAreaDna::new(dna, Some("seq_gene".to_string()), Some(engine));

    let opened = area.open_splicing_expert_for_feature(0, "test");

    assert!(opened);
    assert!(area.show_splicing_expert_window);
    assert!(area.splicing_expert_window_pending_initial_render);
    assert_eq!(area.splicing_expert_window_feature_id, Some(0));
}

#[test]
fn open_splicing_expert_window_marks_initial_render_pending() {
    let dna = DNAsequence::from_sequence("ACGT").expect("sequence");
    let mut area = MainAreaDna::new(dna, Some("seq1".to_string()), None);
    let view = SplicingExpertView {
        seq_id: "seq1".to_string(),
        target_feature_id: 17,
        scope: SplicingScopePreset::AllOverlappingAnyStrand,
        group_label: "TP73".to_string(),
        strand: "+".to_string(),
        region_start_1based: 1,
        region_end_1based: 4,
        transcript_count: 0,
        unique_exon_count: 0,
        instruction: String::new(),
        transcripts: vec![],
        unique_exons: vec![],
        matrix_rows: vec![],
        boundaries: vec![],
        intron_signals: vec![],
        junctions: vec![],
        events: vec![],
    };

    area.open_splicing_expert_window_for_view(&view);

    assert!(area.show_splicing_expert_window);
    assert!(area.splicing_expert_window_pending_initial_render);
    assert_eq!(area.splicing_expert_window_feature_id, Some(17));
}

#[test]
fn tooltip_help_splicing_expert_window_mentions_transcripts_and_rna_reads() {
    let help = MainAreaDna::splicing_expert_window_help_text();
    assert!(help.contains("transcript"));
    assert!(help.contains("primer/qPCR"));
    assert!(help.contains("RNA-read evidence"));
}

#[test]
fn tooltip_help_nanopore_cdna_panel_mentions_two_phase_report_flow() {
    let help = MainAreaDna::splicing_nanopore_cdna_panel_help_text();
    assert!(help.contains("two-phase"));
    assert!(help.contains("shared across mapping profiles"));
    assert!(help.contains("InterpretRnaReads"));
    assert!(help.contains("AlignRnaReadReport"));
    assert!(help.contains("Report ID"));
    assert!(help.contains("whole-genome mapper"));
}

#[test]
fn rna_read_mapping_labels_are_profile_generic() {
    assert_eq!(
        MainAreaDna::rna_read_mapping_parameter_section_title(),
        "RNA-read mapping parameters"
    );
    assert_eq!(
        MainAreaDna::rna_read_mapping_run_button_label(),
        "Run RNA-read interpretation"
    );
    assert!(!MainAreaDna::rna_read_mapping_parameter_section_title().contains("Nanopore"));
    assert!(!MainAreaDna::rna_read_mapping_run_button_label().contains("Nanopore"));
}

#[test]
fn default_rna_align_selection_prefers_seed_passed_rows() {
    assert_eq!(
        super::default_rna_align_selection(),
        RnaReadHitSelection::SeedPassed
    );
    assert_eq!(
        super::RnaReadInterpretOpsUiState::default().align_phase_selection,
        RnaReadHitSelection::SeedPassed
    );
}

#[test]
fn rna_read_align_selection_ui_labels_are_distinct_and_explicit() {
    assert_eq!(
        super::MainAreaDna::rna_read_align_selection_ui_label(RnaReadHitSelection::SeedPassed),
        "seed_passed"
    );
    assert_eq!(
        super::MainAreaDna::rna_read_align_selection_ui_label(RnaReadHitSelection::All),
        "all retained"
    );
    assert_eq!(
        super::MainAreaDna::rna_read_align_selection_ui_label(RnaReadHitSelection::Aligned),
        "already_aligned"
    );
}

#[test]
fn rna_read_progress_eta_text_reports_remaining_time_for_known_working_set() {
    assert_eq!(
        super::MainAreaDna::format_rna_read_progress_eta(25, 100, 10.0).as_deref(),
        Some("ETA: 30s")
    );
    assert_eq!(
        super::MainAreaDna::format_rna_read_progress_eta(0, 100, 10.0),
        None
    );
    assert_eq!(
        super::MainAreaDna::format_rna_read_progress_eta(100, 100, 10.0),
        Some("ETA: 0s".to_string())
    );
}

#[cfg(test)]
mod embedded_scope_tests {
    use super::MainAreaDna;
    use crate::dna_sequence::DNAsequence;

    #[test]
    fn panel_scope_key_uses_window_scope_id_when_present() {
        let dna = DNAsequence::from_sequence("ACGT").expect("valid DNA");
        let mut area = MainAreaDna::new(dna, Some("seq1".to_string()), None);
        assert_eq!(area.panel_scope_key(), "seq1");
        area.set_window_scope_id("viewport-42");
        assert_eq!(area.panel_scope_key(), "seq1::viewport-42");
    }
}
