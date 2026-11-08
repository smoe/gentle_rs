//! Synthetic core cases shared with the unchanged Python projector; no I/O at runtime.

use super::*;
use serde_json::{Value, json};

const CORPUS: &str = include_str!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/test_files/rna_read_end_core/cases.json"
));

pub(super) fn normalized(
    document: &ModelDocument,
    analysis: &ObservedSamAnalysis,
) -> crate::engine::PrimerPayloadAnnotation {
    use crate::engine::{
        PRIMER_PAYLOAD_ANNOTATION_SCHEMA, PrimerPayloadAnnotation, PrimerPayloadTranscript,
    };
    let mut transcripts = BTreeMap::new();
    for model in &document.transcripts {
        transcripts
            .entry(model.transcript_id.clone())
            .or_insert(PrimerPayloadTranscript {
                transcript_id: model.transcript_id.clone(),
                gene_ids: vec![model.gene_id.clone()],
                target_name: model.contig.clone(),
                target_length: usize::MAX,
                genomic_contig: Some(model.contig.clone()),
                genomic_strand: Some(model.strand),
                genomic_exons: model.genomic_exons_0based_half_open.clone(),
            });
    }
    for observation in &analysis.observations {
        for source in &observation.comparisons {
            for model in source.comparison.iter().flat_map(|c| &c.models) {
                let fallback = document
                    .transcripts
                    .first()
                    .map(|m| m.genomic_exons_0based_half_open.clone())
                    .unwrap_or_default();
                let value = transcripts.entry(model.transcript_id.clone()).or_insert(
                    PrimerPayloadTranscript {
                        transcript_id: model.transcript_id.clone(),
                        gene_ids: model.gene_ids.clone(),
                        target_name: observation.contig.clone(),
                        target_length: usize::MAX,
                        genomic_contig: Some(observation.contig.clone()),
                        genomic_strand: Some(model.model_strand),
                        genomic_exons: fallback,
                    },
                );
                value.gene_ids.extend(model.gene_ids.clone());
                value.gene_ids.sort();
                value.gene_ids.dedup();
            }
        }
    }
    PrimerPayloadAnnotation {
        schema: PRIMER_PAYLOAD_ANNOTATION_SCHEMA.into(),
        reference_id: "synthetic_reference".into(),
        assembly: document.assembly.clone(),
        release: "synthetic_snapshot".into(),
        reference_fasta_sha256: document.reference_fasta_sha256.clone(),
        coordinate_convention: document.coordinate_convention.clone(),
        transcripts: transcripts.into_values().collect(),
    }
}

pub(super) fn project_read_ends(
    document: &ModelDocument,
    analysis: &ObservedSamAnalysis,
    limits: Limits,
) -> Result<Projection, EngineError> {
    super::project_read_ends(document, analysis, &normalized(document, analysis), limits)
}

fn inputs(corpus: &Value, case: &Value) -> (ModelDocument, ObservedSamAnalysis) {
    let mut document = corpus["binding"].clone();
    document["transcripts"] = Value::Array(
        case["models"]
            .as_array()
            .unwrap()
            .iter()
            .map(|name| corpus["models"][name.as_str().unwrap()].clone())
            .collect(),
    );
    let unit = json!({
        "source_sha256": format!("sha256:{}", "d".repeat(64)),
        "unit_id": "synthetic_unit"
    });
    let observations: Vec<Value> = case["observations"]
        .as_array()
        .unwrap()
        .iter()
        .enumerate()
        .map(|(index, row)| {
            let raw = row["blocks"].as_array().unwrap();
            let blocks: Vec<Value> = raw
                .iter()
                .map(|b| json!({"query": [b[0], b[1]], "target": [b[2], b[3]]}))
                .collect();
            let models: Vec<Value> = row["candidates"]
                .as_array()
                .unwrap()
                .iter()
                .map(|id| {
                    let strand = document["transcripts"]
                        .as_array()
                        .unwrap()
                        .iter()
                        .find(|m| m["transcript_id"] == *id)
                        .map(|m| m["strand"].clone())
                        .unwrap_or(json!("+"));
                    json!({
                        "transcript_id": id, "gene_ids": ["GENE_SYN_A"],
                        "model_strand": strand,
                        "relation": row.get("relation").cloned().unwrap_or(json!("exact_partial_chain")),
                        "retained_intron_supported_bases": null
                    })
                })
                .collect();
            json!({
                "observation_id": row["id"], "unit": unit, "source_ordinal": index + 1,
                "on_optimal_path": true, "contig": "synthetic_contig",
                "observed_span": [raw.iter().map(|b| b[2].as_u64().unwrap()).min().unwrap(),
                                  raw.iter().map(|b| b[3].as_u64().unwrap()).max().unwrap()],
                "aligned_query_interval": [raw.iter().map(|b| b[0].as_u64().unwrap()).min().unwrap(),
                                           raw.iter().map(|b| b[1].as_u64().unwrap()).max().unwrap()],
                "query_axis_length": row["axis"], "alignment_clipped": [false, false],
                "query_to_reference_strand": row["strand"], "reported_mapq": null,
                "ends": [[],[]], "detail": "detailed", "topology": "unassessed",
                "blocks": blocks, "excluded_query_intervals": null, "skips": [],
                "flank_qualification": [], "observed_junction_chain": null, "clipping": null,
                "comparisons": [{
                    "source_id": document["source_id"], "availability": "available",
                    "comparison": {
                        "source_id": document["source_id"], "gene_id": "GENE_SYN_A",
                        "resolution": "resolved", "junction_relation": null,
                        "junctions": null, "models": models
                    }
                }]
            })
        })
        .collect();
    let analysis = json!({
        "context_sha256": format!("sha256:{}", "c".repeat(64)),
        "reference": {"assembly": document["assembly"], "fasta_sha256": document["reference_fasta_sha256"]},
        "reference_gene": {"source_id": document["source_id"], "gene_id": "GENE_SYN_A"},
        "annotation_bindings": [{
            "source_id": document["source_id"], "provider": null,
            "annotation_sha256": document["annotation_sha256"], "reference_id": "synthetic_reference",
            "assembly": document["assembly"], "release": "synthetic_snapshot",
            "reference_fasta_sha256": document["reference_fasta_sha256"]
        }],
        "gene_resolutions": [], "observations": observations,
        "units": [{
            "unit": unit, "unclipped_axis_length": case["axis"], "unavailable": [], "rows": [],
            "search_declared_complete": false, "views_declared_complete": false,
            "all_query_views_assessed": false, "hypotheses": case["hypotheses"],
            "reconciliation": "unassessed", "topologies": [],
            "definite_chains": [], "possible_chains": []
        }],
        "chain_families": [], "junctions": [], "comparison_work": 0, "projection_work": 0
    });
    (
        serde_json::from_value(document).unwrap(),
        serde_json::from_value(analysis).unwrap(),
    )
}

pub(super) fn fixture(name: &str) -> (ModelDocument, ObservedSamAnalysis) {
    let corpus: Value = serde_json::from_str(CORPUS).unwrap();
    let case = corpus["cases"]
        .as_array()
        .unwrap()
        .iter()
        .find(|case| case["name"] == name)
        .unwrap();
    inputs(&corpus, case)
}

fn triple(metric: &Metric) -> Value {
    match metric {
        Metric::Observed { value_bp } => json!([value_bp, "observed", "resolved"]),
        Metric::RightCensored { lower_bound_bp } => {
            json!([lower_bound_bp, "right_censored", "clipped_transcript_end"])
        }
        Metric::Missing { reason } => json!([null, "missing", reason]),
    }
}

fn semantic(unit: &UnitProjection) -> Value {
    json!({
        "five": triple(&unit.utr5), "three": triple(&unit.utr3),
        "ids": unit.candidate_transcript_ids, "genes": unit.model_gene_ids
    })
}

fn error(document: &ModelDocument, analysis: &ObservedSamAnalysis) -> String {
    project_read_ends(document, analysis, Limits::default())
        .unwrap_err()
        .message
}

#[test]
fn utr_core_shared_python_corpus_matches_literals_and_preserves_inputs() {
    let corpus: Value = serde_json::from_str(CORPUS).unwrap();
    let cases = corpus["cases"].as_array().unwrap();
    assert_eq!(cases.len(), 31);
    for case in cases {
        let (document, analysis) = inputs(&corpus, case);
        let before = serde_json::to_value((&document, &analysis)).unwrap();
        let out = project_read_ends(&document, &analysis, Limits::default()).unwrap();
        assert_eq!(
            semantic(&out.units[0]),
            case["expected"],
            "{}",
            case["name"]
        );
        assert_eq!(out.units[0].unit, analysis.units[0].unit);
        assert_eq!(
            serde_json::to_value((&document, &analysis)).unwrap(),
            before
        );
        assert!(out.work > 0);
    }
}

#[test]
fn utr_core_v2_requires_explicit_stop_field_and_strict_shapes() {
    let (document, _) = fixture("plus");
    let raw = serde_json::to_value(document).unwrap();
    for field in ["genomic_stop_codon_0based_half_open", "strand"] {
        let mut altered = raw.clone();
        altered["transcripts"][0]
            .as_object_mut()
            .unwrap()
            .remove(field);
        assert!(serde_json::from_value::<ModelDocument>(altered).is_err());
    }
    for value in [Value::Null, json!("unknown")] {
        let mut altered = raw.clone();
        altered["transcripts"][0]["genomic_stop_codon_0based_half_open"] = value;
        assert!(serde_json::from_value::<ModelDocument>(altered).is_err());
    }
    for field in ["search_complete", "polya_tail_length"] {
        let mut altered = raw.clone();
        altered[field] = json!(true);
        assert!(serde_json::from_value::<ModelDocument>(altered).is_err());
    }
    let mut altered = raw;
    altered["transcripts"][0]["extra"] = json!(true);
    assert!(serde_json::from_value::<ModelDocument>(altered).is_err());
}

#[test]
fn utr_core_unknown_schema_and_legacy_are_not_guessed_as_v2() {
    let (mut document, analysis) = fixture("plus");
    for schema in ["gentle.rna_read_end_transcript_models.v1", "unknown"] {
        document.schema = schema.into();
        assert_eq!(error(&document, &analysis), "rna_read_end_model_document");
    }
}

#[test]
fn utr_core_reference_and_exact_source_binding_are_required() {
    let (document, analysis) = fixture("plus");
    let mut other = analysis.clone();
    other.reference.assembly.push_str("_other");
    assert_eq!(error(&document, &other), "rna_read_end_reference_binding");
    other = analysis.clone();
    other.reference.fasta_sha256 = format!("sha256:{}", "e".repeat(64));
    assert_eq!(error(&document, &other), "rna_read_end_reference_binding");
    other = analysis.clone();
    other.annotation_bindings.clear();
    assert_eq!(error(&document, &other), "rna_read_end_annotation_binding");
    other = analysis.clone();
    other
        .annotation_bindings
        .push(other.annotation_bindings[0].clone());
    assert_eq!(error(&document, &other), "rna_read_end_annotation_binding");
    other = analysis;
    other.annotation_bindings[0].annotation_sha256 = format!("sha256:{}", "e".repeat(64));
    assert_eq!(error(&document, &other), "rna_read_end_annotation_binding");
}

#[test]
fn utr_core_invalid_model_geometry_and_duplicate_ids_refuse() {
    let (document, analysis) = fixture("plus");
    for exons in [
        vec![(100, 110), (110, 140)],
        vec![(100, 110), (105, 140)],
        vec![(120, 140), (100, 110)],
    ] {
        let mut changed = document.clone();
        changed.transcripts[0].genomic_exons_0based_half_open = exons;
        assert_eq!(error(&changed, &analysis), "rna_read_end_model_intervals");
    }
    for stops in [
        vec![(110, 113)],
        vec![(122, 125)],
        vec![(130, 134)],
        vec![(125, 126), (130, 132)],
    ] {
        let mut changed = document.clone();
        changed.transcripts[0].genomic_stop_codon_0based_half_open = stops;
        assert!(error(&changed, &analysis).starts_with("rna_read_end_stop_"));
    }
    let mut changed = document.clone();
    changed.transcripts[0].genomic_cds_0based_half_open = vec![(109, 121)];
    assert_eq!(error(&changed, &analysis), "rna_read_end_cds_outside_exons");
    changed = document;
    changed.transcripts.push(changed.transcripts[0].clone());
    assert_eq!(error(&changed, &analysis), "rna_read_end_duplicate_model");
}

#[test]
fn utr_core_complete_partial_and_split_stop_boundaries_are_explicit() {
    for (name, expected) in [
        ("plus", StopEvidence::Included),
        ("separate_stop", StopEvidence::Separate),
        ("overlapping_stop", StopEvidence::Overlapping),
        ("split_plus_stop", StopEvidence::Separate),
        ("split_minus_stop", StopEvidence::Separate),
        ("partial_stop", StopEvidence::Partial),
        ("unknown_stop", StopEvidence::Unknown),
    ] {
        let (document, analysis) = fixture(name);
        let out = project_read_ends(&document, &analysis, Limits::default()).unwrap();
        assert_eq!(out.models[0].stop_evidence, expected, "{name}");
        assert_eq!(
            out.models[0].utr3_start.is_none(),
            matches!(expected, StopEvidence::Unknown | StopEvidence::Partial)
        );
    }
}

#[test]
fn utr_core_wrong_unit_axis_nonoptimal_or_repeated_path_nodes_refuse() {
    let (document, analysis) = fixture("plus");
    let mut changed = analysis.clone();
    changed.observations[0].unit.unit_id = "a_different_unit".into();
    assert_eq!(error(&document, &changed), "rna_read_end_unit_binding");
    changed = analysis.clone();
    changed.units[0].unclipped_axis_length = Some(26);
    assert_eq!(error(&document, &changed), "rna_read_end_unit_binding");
    changed = analysis.clone();
    changed.observations[0].on_optimal_path = false;
    assert_eq!(error(&document, &changed), "rna_read_end_unit_binding");
    changed = analysis;
    changed.units[0].hypotheses[0].push("o1".into());
    assert_eq!(
        error(&document, &changed),
        "rna_read_end_hypothesis_identity"
    );
}

#[test]
fn utr_core_duplicate_observations_units_and_invalid_axes_refuse() {
    let (document, analysis) = fixture("plus");
    let mut changed = analysis.clone();
    changed.observations.push(changed.observations[0].clone());
    assert_eq!(
        error(&document, &changed),
        "rna_read_end_observation_identity"
    );
    changed = analysis.clone();
    changed.units.push(changed.units[0].clone());
    assert_eq!(error(&document, &changed), "rna_read_end_unit_identity");
    changed = analysis;
    changed.units[0].unclipped_axis_length = Some(0);
    assert_eq!(error(&document, &changed), "rna_read_end_unit_axis");
}

#[test]
fn utr_core_bad_blocks_refuse_instead_of_wrapping_or_extrapolating() {
    let (document, mut analysis) = fixture("plus");
    for range in [(10, 10), (20, 10), (0, 26)] {
        analysis.observations[0].blocks[0].query = range;
        assert_eq!(
            error(&document, &analysis),
            "rna_read_end_observation_geometry"
        );
    }
}

#[test]
fn utr_core_unknown_axis_does_not_become_zero_or_change_end_geometry() {
    let (document, mut analysis) = fixture("plus");
    analysis.units[0].unclipped_axis_length = None;
    let out = project_read_ends(&document, &analysis, Limits::default()).unwrap();
    assert_eq!(
        out.units[0].read_length,
        missing(Reason::AxisLengthUnavailable)
    );
    assert_eq!(out.units[0].utr5, Metric::Observed { value_bp: 5 });
    assert_eq!(out.units[0].utr3, Metric::Observed { value_bp: 5 });
}

#[test]
fn utr_core_reordered_candidates_models_and_observations_keep_semantics() {
    let (mut document, mut analysis) = fixture("agreeing_candidates");
    let before = project_read_ends(&document, &analysis, Limits::default()).unwrap();
    document.transcripts.reverse();
    analysis.observations[0].comparisons[0]
        .comparison
        .as_mut()
        .unwrap()
        .models
        .reverse();
    let after = project_read_ends(&document, &analysis, Limits::default()).unwrap();
    assert_eq!(after.units, before.units);
    assert_eq!(after.models, before.models);
    let (document, mut analysis) = fixture("intersection_before_availability");
    let before = project_read_ends(&document, &analysis, Limits::default()).unwrap();
    analysis.observations.reverse();
    assert_eq!(
        project_read_ends(&document, &analysis, Limits::default())
            .unwrap()
            .units,
        before.units
    );
}

#[test]
fn utr_core_work_and_size_limits_refuse_without_partial_output() {
    let (document, analysis) = fixture("plus");
    let out = project_read_ends(&document, &analysis, Limits::default()).unwrap();
    let exact = Limits {
        max_work: out.work,
        ..Limits::default()
    };
    assert_eq!(project_read_ends(&document, &analysis, exact).unwrap(), out);
    assert_eq!(
        project_read_ends(
            &document,
            &analysis,
            Limits {
                max_work: out.work - 1,
                ..exact
            }
        )
        .unwrap_err()
        .message,
        "rna_read_end_work_budget"
    );
    for limits in [
        Limits {
            max_intervals: 1,
            ..Limits::default()
        },
        Limits {
            max_models: 0,
            ..Limits::default()
        },
        Limits {
            max_units: 0,
            ..Limits::default()
        },
        Limits {
            max_observations: 0,
            ..Limits::default()
        },
    ] {
        assert!(project_read_ends(&document, &analysis, limits).is_err());
    }
}

#[test]
fn utr_core_names_are_data_not_current_panel_defaults() {
    let (mut document, mut analysis) = fixture("plus");
    document.source_id = "another_annotation".into();
    document.transcripts[0].transcript_id = "another_transcript".into();
    document.transcripts[0].gene_id = "another_gene".into();
    document.transcripts[0].contig = "another_contig".into();
    analysis.annotation_bindings[0].source_id = document.source_id.clone();
    analysis.reference_gene.source_id = document.source_id.clone();
    analysis.reference_gene.gene_id = "another_gene".into();
    analysis.observations[0].contig = "another_contig".into();
    let comparison = &mut analysis.observations[0].comparisons[0];
    comparison.source_id = document.source_id.clone();
    let data = comparison.comparison.as_mut().unwrap();
    data.source_id = document.source_id.clone();
    data.gene_id = Some("another_gene".into());
    data.models[0].transcript_id = "another_transcript".into();
    data.models[0].gene_ids = vec!["another_gene".into()];
    let out = project_read_ends(&document, &analysis, Limits::default()).unwrap();
    assert_eq!(
        out.units[0].candidate_transcript_ids,
        ["another_transcript"]
    );
    assert_eq!(out.units[0].model_gene_ids, ["another_gene"]);
    assert_eq!(out.units[0].utr3, Metric::Observed { value_bp: 5 });
    let (mut document, analysis) = fixture("agreeing_candidates");
    document.transcripts[1].gene_id = "OTHER_SYN_GENE".into();
    let out = project_read_ends(&document, &analysis, Limits::default()).unwrap();
    assert_eq!(
        out.units[0].model_gene_ids,
        ["GENE_SYN_A", "OTHER_SYN_GENE"]
    );
    assert_eq!(out.units[0].candidate_transcript_ids.len(), 2);
}

#[test]
fn utr_core_existing_sam_kernel_assessments_feed_the_same_projection() {
    let (_root, request, _manifest) = super::super::sam_intake::tests::fixture();
    let (analysis, _) =
        super::super::summarize_sam_bound(&request, &BTreeMap::new(), 1_000_000, &mut || Ok(()))
            .unwrap();
    let (mut document, _) = fixture("plus");
    let source = &request.annotation_sources[0];
    document.source_id = source.source_id.clone();
    document.assembly = request.reference.assembly.clone();
    document.reference_fasta_sha256 = request.reference.fasta_sha256.clone();
    document.annotation_sha256 = source.annotation_sha256.clone();
    let transcript = &source.annotation.transcripts[0];
    document.transcripts[0].transcript_id = transcript.transcript_id.clone();
    document.transcripts[0].gene_id = transcript.gene_ids[0].clone();
    document.transcripts[0].contig = transcript.target_name.clone();
    let before = serde_json::to_value(&analysis).unwrap();
    let out = project_read_ends(&document, &analysis, Limits::default()).unwrap();
    assert_eq!(out.units.len(), 4);
    let projected = out
        .units
        .iter()
        .find(|u| u.unit.unit_id == "unit_a")
        .unwrap();
    assert_eq!(projected.utr5, Metric::Observed { value_bp: 5 });
    assert_eq!(projected.utr3, Metric::Observed { value_bp: 10 });
    assert_eq!(serde_json::to_value(analysis).unwrap(), before);
}

#[test]
fn utr_core_alternative_paths_do_not_choose_the_first_end() {
    let (document, mut analysis) = fixture("plus");
    let mut alternative = analysis.observations[0].clone();
    alternative.observation_id = "other".into();
    alternative.blocks[1].query = (10, 20);
    alternative.blocks[1].target = (120, 130);
    alternative.aligned_query_interval = (0, 20);
    alternative.observed_span = (100, 130);
    analysis.observations.push(alternative);
    analysis.units[0].hypotheses.push(vec!["other".into()]);
    let before = project_read_ends(&document, &analysis, Limits::default()).unwrap();
    assert_eq!(before.units[0].utr5, Metric::Observed { value_bp: 5 });
    assert_eq!(
        before.units[0].utr3,
        missing(Reason::AmbiguousTranscriptEndGeometry)
    );
    analysis.units[0].hypotheses.reverse();
    let after = project_read_ends(&document, &analysis, Limits::default()).unwrap();
    assert_eq!(semantic(&after.units[0]), semantic(&before.units[0]));
    assert_eq!(after.units[0].hypotheses.len(), 2);
    assert_eq!(after.units[0].hypotheses[0].observation_ids, ["other"]);
}

#[test]
fn utr_core_unavailable_or_unresolved_comparison_does_not_assign_a_model() {
    let (document, analysis) = fixture("plus");
    for change in 0..4 {
        let mut changed = analysis.clone();
        let row = &mut changed.observations[0].comparisons[0];
        match change {
            0 => row.availability = ObservedSamAnnotationAvailability::Refused,
            1 => row.comparison.as_mut().unwrap().resolution = GeneResolution::AmbiguousLocus,
            2 => row.comparison.as_mut().unwrap().models[0].relation = ModelRelation::Intronic,
            _ => row.source_id = "unrelated_annotation".into(),
        }
        let out = project_read_ends(&document, &changed, Limits::default()).unwrap();
        assert_eq!(out.units[0].utr3, missing(Reason::NoCompatibleCdsModel));
        assert!(out.units[0].candidate_transcript_ids.is_empty());
        assert!(out.units[0].model_gene_ids.is_empty());
    }
}

#[test]
fn utr_core_large_genomic_coordinates_do_not_wrap_transcript_offsets() {
    let (mut document, mut analysis) = fixture("plus");
    let end = usize::MAX;
    document.transcripts[0].genomic_exons_0based_half_open = vec![(end - 100, end)];
    document.transcripts[0].genomic_cds_0based_half_open = vec![(end - 90, end - 20)];
    document.transcripts[0].genomic_stop_codon_0based_half_open = vec![(end - 23, end - 20)];
    analysis.units[0].unclipped_axis_length = Some(100);
    analysis.observations[0].query_axis_length = 100;
    analysis.observations[0].aligned_query_interval = (0, 100);
    analysis.observations[0].observed_span = (end - 100, end);
    analysis.observations[0].blocks = vec![SupportedBlock {
        query: (0, 100),
        target: (end - 100, end),
    }];
    let out = project_read_ends(&document, &analysis, Limits::default()).unwrap();
    assert_eq!(out.models[0].transcript_length, 100);
    assert_eq!(out.units[0].utr5, Metric::Observed { value_bp: 10 });
    assert_eq!(out.units[0].utr3, Metric::Observed { value_bp: 20 });
}

#[test]
fn utr_core_model_must_match_actual_normalized_annotation_geometry() {
    let (document, analysis) = fixture("plus");
    let annotation = normalized(&document, &analysis);
    for change in 0..4 {
        let mut changed = annotation.clone();
        let model = &mut changed.transcripts[0];
        match change {
            0 => model.genomic_contig = Some("other_contig".into()),
            1 => model.genomic_strand = Some('-'),
            2 => model.genomic_exons[0].0 += 1,
            _ => model.gene_ids = vec!["unrelated_gene".into()],
        }
        assert_eq!(
            super::project_read_ends(&document, &analysis, &changed, Limits::default())
                .unwrap_err()
                .message,
            "rna_read_end_model_annotation_geometry"
        );
    }
}

#[test]
fn utr_core_records_are_uncollapsed_and_retain_noncoding_annotation_genes() {
    let (document, analysis) = fixture("missing_candidate_model");
    let out = project_read_ends(&document, &analysis, Limits::default()).unwrap();
    let unit = &out.units[0];
    assert_eq!(unit.annotation_gene_ids, ["GENE_SYN_A"]);
    assert!(unit.model_gene_ids.is_empty());
    assert_eq!(unit.hypotheses[0].candidates.len(), 2);
    let known = unit.hypotheses[0]
        .candidates
        .iter()
        .find(|c| c.transcript_id == "tx_plus")
        .unwrap();
    assert_eq!(known.utr3, Metric::Observed { value_bp: 5 });
    assert_eq!(known.transcript_intervals, [(0, 25)]);
    assert_eq!(unit.utr3, missing(Reason::CandidateModelUnavailable));
    assert_eq!(known.relations[0].observation_id, "o1");
}

#[test]
fn utr_core_clip_kinds_are_additive_without_reinterpreting_censoring() {
    let corpus: Value = serde_json::from_str(CORPUS).unwrap();
    let cases = corpus["clip_cases"].as_array().unwrap();
    assert_eq!(cases.len(), 6);
    for entry in cases {
        let (document, mut analysis) = fixture(entry["case"].as_str().unwrap());
        if let Some(axis) = entry.get("axis").and_then(Value::as_u64) {
            analysis.units[0].unclipped_axis_length = Some(axis as usize);
            for observation in &mut analysis.observations {
                observation.query_axis_length = axis as usize;
            }
        }
        analysis.observations[0].clipping = if entry["soft"].is_null() {
            None
        } else {
            Some(ObservedSamClipping {
                soft: serde_json::from_value(entry["soft"].clone()).unwrap(),
                hard: serde_json::from_value(entry["hard"].clone()).unwrap(),
                padding_operations: 0,
            })
        };
        let out = project_read_ends(&document, &analysis, Limits::default()).unwrap();
        let candidate = &out.units[0].hypotheses[0].candidates[0];
        assert_eq!(
            serde_json::to_value(&candidate.clip3).unwrap(),
            entry["end3"],
            "{entry}"
        );
        assert_eq!(triple(&candidate.utr3), entry["utr3"]);
        assert_eq!(triple(&out.units[0].utr3), entry["utr3"]);
        assert_eq!(candidate.clip5.kind, ClipKind::None);
    }
}
