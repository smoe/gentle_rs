//! Inline invented SAM, locus, annotation and base/quality fixtures; no biological samples.
//! Helpers recreate all raw-byte-bound inputs and genuine saved S2 artifacts in private tempdirs.
//! The tests exercise comparison ingress and shared routes, not external aligners or read acquisition.

use super::*;
use crate::engine::{
    Engine, GentleEngine, OBSERVED_SAM_QUERY_UNIT_MANIFEST_SCHEMA, ObservedSamAnnotationBinding,
    ObservedSamMaskWitness, ObservedSamQueryUnit, ObservedSamQueryUnitManifest,
    ObservedSamRowWitness, ObservedSamSearchDeclaration, Operation,
};
use crate::engine_shell::{execute_shell_command, parse_shell_line, quote_shell_arg};
use std::{
    fs,
    path::{Path, PathBuf},
};

struct Fixture {
    _root: tempfile::TempDir,
    request: CdnaFileComparisonRequest,
    original_inputs: Vec<PathBuf>,
}

fn bound(path: &Path) -> PrimerPayloadFile {
    let bytes = fs::read(path).unwrap();
    PrimerPayloadFile {
        path: path.to_str().unwrap().into(),
        sha256: crate::digest_utils::sha256_prefixed_bytes(&bytes),
        bytes: bytes.len() as u64,
    }
}

fn save<T: Serialize>(path: &Path, value: &T) -> PrimerPayloadFile {
    fs::write(path, serde_json::to_vec_pretty(value).unwrap()).unwrap();
    bound(path)
}

fn reverse_complement(s: &str) -> String {
    s.bytes()
        .rev()
        .map(|c| match c {
            b'A' => 'T',
            b'C' => 'G',
            b'G' => 'C',
            b'T' => 'A',
            _ => 'N',
        })
        .collect()
}

fn fixture(reverse: bool, two_samples: bool) -> Fixture {
    let (root, mut template) = GentleEngine::observed_sam_test_fixture();
    let directory = root.path().canonicalize().unwrap();
    let annotation_file = match &template.annotation_sources[0] {
        ObservedSamAnnotationBinding::Normalized { annotation, .. } => annotation.clone(),
        _ => panic!("literal normalized fixture"),
    };
    let mut annotation: PrimerPayloadAnnotation =
        serde_json::from_slice(&fs::read(&annotation_file.path).unwrap()).unwrap();
    annotation.transcripts[0].genomic_strand = Some(if reverse { '-' } else { '+' });
    let annotation_file = save(
        &directory.join("comparison's annotation α.json"),
        &annotation,
    );
    match &mut template.annotation_sources[0] {
        ObservedSamAnnotationBinding::Normalized { annotation, .. } => {
            *annotation = annotation_file.clone()
        }
        _ => unreachable!(),
    }
    let locus = directory.join("comparison's locus α.fa");
    fs::write(
        &locus,
        format!(">ctg_a synthetic locus\n{}\n", "A".repeat(50)),
    )
    .unwrap();
    let mut samples = Vec::new();
    let mut original_inputs = vec![PathBuf::from(&template.dictionary.path)];
    for index in 0..if two_samples { 2 } else { 1 } {
        let mut source = template.clone();
        let name = format!("synthetic_{index}");
        let qname = format!("invented_read_{index}");
        let sam = directory.join(format!("{name}'s input.sam"));
        let flag = if reverse { 16 } else { 0 };
        let raw = format!(
            "@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:ctg_a\tLN:1000\n@SQ\tSN:ctg_b\tLN:1000\n@RG\tID:library\n{qname}\t{flag}\tctg_a\t101\t60\t10M10N20M\t*\t0\t0\t*\t*\tRG:Z:library\n"
        );
        fs::write(&sam, raw).unwrap();
        source.sam = bound(&sam);
        let manifest = ObservedSamQueryUnitManifest {
            schema: OBSERVED_SAM_QUERY_UNIT_MANIFEST_SCHEMA.into(),
            sam_sha256: source.sam.sha256.clone(),
            units: vec![ObservedSamQueryUnit {
                unit_id: name.clone(),
                unclipped_axis_length: Some(30),
                rows: vec![ObservedSamRowWitness {
                    ordinal: 1,
                    qname: qname.clone(),
                    read_group: Some("library".into()),
                }],
                masks: Some(ObservedSamMaskWitness {
                    coordinate_convention: "zero_based_half_open_unclipped_query".into(),
                    unclipped_axis_length: 30,
                    intervals: vec![],
                    source_sha256: format!("sha256:{}", "c".repeat(64)),
                    derivation: "Invented fixture explicitly declares no excluded bases".into(),
                }),
                ends: [vec![], vec![]],
                search: Some(ObservedSamSearchDeclaration {
                    search_complete: true,
                    views_complete: true,
                    limitations: vec![],
                    producer: "Literal fixture, not an aligner".into(),
                    version: "synthetic-v1".into(),
                    arguments: vec![],
                    selection_scope: "All invented competing rows supplied".into(),
                }),
            }],
        };
        source.manifest = save(&directory.join(format!("{name}'s units.json")), &manifest);
        source.output_dir = directory
            .join(format!("{name}'s saved analysis"))
            .to_str()
            .unwrap()
            .into();
        let engine = GentleEngine::new();
        let report = engine
            .summarize_observed_transcript_structures_from_sam(source.clone(), &mut |_| true)
            .unwrap();
        let analysis_binding = GentleEngine::observed_sam_result_binding(&source, &report).unwrap();
        let (_, data_file) = crate::rna_read_ends::analysis_file(&analysis_binding).unwrap();
        let data: ObservedSamData =
            crate::rna_read_ends::bound_json(&data_file, &mut || Ok(())).unwrap();
        let analysis = data.analysis.unwrap();
        let mut sequence = "A".repeat(30);
        if index == 0 {
            sequence.replace_range(3..4, "C");
            sequence.replace_range(15..16, "G");
        }
        if reverse {
            sequence = reverse_complement(&sequence);
        }
        let document = CdnaBasecallDocument {
            schema: CDNA_BASECALL_SCHEMA.into(),
            analysis_context_sha256: analysis.context_sha256,
            orientation: "original_sequencing_5_to_3".into(),
            provenance: "Hand-crafted original sequence/quality; no physical sequencing claim"
                .into(),
            reads: vec![CdnaOriginalBasecall {
                unit: analysis.units[0].unit.clone(),
                read_id: qname,
                sequence,
                qualities_phred33: Some("I".repeat(30)),
            }],
        };
        samples.push(CdnaComparisonSampleFiles {
            sample_id: name,
            label: format!("Invented cell-line-like sample {index}; synthetic"),
            analysis: analysis_binding,
            basecalls: save(
                &directory.join(format!("sample {index}'s basecalls α.json")),
                &document,
            ),
            additional_partitions: vec![],
            partition_source_sha256: None,
        });
        original_inputs.extend([
            PathBuf::from(source.sam.path),
            PathBuf::from(source.manifest.path),
        ]);
    }
    let mut policy = core::CdnaComparisonPolicy::default();
    policy.minimum_alternate_units = 1;
    policy.minimum_haplotype_units = 1;
    Fixture {
        _root: root,
        request: CdnaFileComparisonRequest {
            schema: CDNA_FILE_COMPARISON_REQUEST_SCHEMA.into(),
            gene: "invented_gene_display".into(),
            source_id: template.reference_gene.source_id,
            gene_id: template.reference_gene.gene_id,
            reference: CdnaLocusReference {
                assembly: template.reference.assembly,
                contig: "ctg_a".into(),
                start_0based: 100,
                genomic_reference_sha256: template.reference.fasta_sha256,
                locus_fasta: bound(&locus),
            },
            normalized_annotation: annotation_file,
            samples,
            supplied_sites: vec![],
            phase_marker_positions: None,
            policy,
            budgets: core::CdnaComparisonBudgets::default(),
            max_input_bytes: default_input_bytes(),
            max_basecall_bases: default_basecall_bases(),
        },
        original_inputs,
    }
}

fn run(request: CdnaFileComparisonRequest) -> CdnaFileComparisonReport {
    compare_cdna_files(request, &mut |_| true).unwrap()
}

fn rewrite_calls(
    request: &mut CdnaFileComparisonRequest,
    change: impl FnOnce(&mut CdnaBasecallDocument),
) {
    let file = &request.samples[0].basecalls;
    let mut calls: CdnaBasecallDocument =
        serde_json::from_slice(&fs::read(&file.path).unwrap()).unwrap();
    change(&mut calls);
    request.samples[0].basecalls = save(Path::new(&file.path), &calls);
}

const PARTITION_BAM_SHA: &str =
    "sha256:aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa";

/// Recreate one invented original-BAM record through the real saved-SAM producer.
/// Local SAM row ordinals restart at one; the producer unit ID names the original record.
fn partition(
    fixture: &Fixture,
    index: usize,
    original_unit_id: &str,
    alternate: bool,
    bam_sha: &str,
    annotation_override: Option<PrimerPayloadFile>,
) -> CdnaComparisonPartitionFiles {
    let directory = fixture._root.path().canonicalize().unwrap();
    let original: crate::engine::ObservedSamReport = serde_json::from_slice(
        &fs::read(&fixture.request.samples[0].analysis.report_path).unwrap(),
    )
    .unwrap();
    let mut source = original.request;
    let qname = format!("invented_partition_read_{index}");
    let sam = directory.join(format!("partition {index}'s primary α.sam"));
    fs::write(
        &sam,
        format!(
            "@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:ctg_a\tLN:1000\n@SQ\tSN:ctg_b\tLN:1000\n@RG\tID:library\n{qname}\t0\tctg_a\t101\t60\t10M10N20M\t*\t0\t0\t*\t*\tRG:Z:library\n"
        ),
    )
    .unwrap();
    source.sam = bound(&sam);
    source.conversion = crate::engine::ObservedSamConversion::BamConversion {
        bam_sha256: bam_sha.into(),
        producer: "Invented partition conversion".into(),
        version: "synthetic-v1".into(),
        arguments: vec![],
        selection_scope: "One explicitly named original record in a synthetic partition".into(),
    };
    if let Some(annotation) = annotation_override {
        match &mut source.annotation_sources[0] {
            ObservedSamAnnotationBinding::Normalized {
                annotation: file, ..
            } => *file = annotation,
            _ => unreachable!(),
        }
    }
    let mut manifest: ObservedSamQueryUnitManifest =
        serde_json::from_slice(&fs::read(&source.manifest.path).unwrap()).unwrap();
    manifest.sam_sha256 = source.sam.sha256.clone();
    manifest.units[0].unit_id = original_unit_id.into();
    manifest.units[0].rows[0].qname = qname.clone();
    let declaration = manifest.units[0].search.as_mut().unwrap();
    declaration.search_complete = false;
    declaration.views_complete = false;
    declaration.limitations = vec!["Only supplied primary records are available".into()];
    source.manifest = save(
        &directory.join(format!("partition {index}'s units.json")),
        &manifest,
    );
    source.output_dir = directory
        .join(format!("partition {index}'s saved evidence"))
        .to_str()
        .unwrap()
        .into();
    let report = GentleEngine::new()
        .summarize_observed_transcript_structures_from_sam(source.clone(), &mut |_| true)
        .unwrap();
    let binding = GentleEngine::observed_sam_result_binding(&source, &report).unwrap();
    let (_, data_file) = crate::rna_read_ends::analysis_file(&binding).unwrap();
    let data: ObservedSamData =
        crate::rna_read_ends::bound_json(&data_file, &mut || Ok(())).unwrap();
    let analysis = data.analysis.unwrap();
    let mut sequence = "A".repeat(30);
    if alternate {
        sequence.replace_range(3..4, "C");
        sequence.replace_range(15..16, "G");
    }
    CdnaComparisonPartitionFiles {
        analysis: binding,
        basecalls: save(
            &directory.join(format!("partition {index}'s basecalls.json")),
            &CdnaBasecallDocument {
                schema: CDNA_BASECALL_SCHEMA.into(),
                analysis_context_sha256: analysis.context_sha256,
                orientation: "original_sequencing_5_to_3".into(),
                provenance: "Invented original BAM ordinal and full bases/Phred40; no real sample"
                    .into(),
                reads: vec![CdnaOriginalBasecall {
                    unit: analysis.units[0].unit.clone(),
                    read_id: qname,
                    sequence,
                    qualities_phred33: Some("I".repeat(30)),
                }],
            },
        ),
    }
}

fn partitioned_fixture() -> Fixture {
    let mut fixture = fixture(false, false);
    let mut parts = (0..4)
        .map(|index| {
            partition(
                &fixture,
                index,
                &format!("bam_record_{:08}", index + 1),
                index < 2,
                PARTITION_BAM_SHA,
                None,
            )
        })
        .collect::<Vec<_>>();
    let first = parts.remove(0);
    fixture.request.samples[0].analysis = first.analysis;
    fixture.request.samples[0].basecalls = first.basecalls;
    fixture.request.samples[0].additional_partitions = parts;
    fixture.request.samples[0].partition_source_sha256 = Some(PARTITION_BAM_SHA.into());
    fixture.request.policy.minimum_alternate_units = 2;
    fixture.request.policy.minimum_haplotype_units = 2;
    fixture
}

/// The same six invented original records, emitted either together or in disjoint chunks.
/// Record identities, bases, Phred bytes and skip patterns are independent of the chunk layout.
fn equivalence_evidence(
    fixture: &Fixture,
    label: &str,
    original_indices: &[usize],
) -> CdnaComparisonPartitionFiles {
    let directory = fixture._root.path().canonicalize().unwrap();
    let original: crate::engine::ObservedSamReport = serde_json::from_slice(
        &fs::read(&fixture.request.samples[0].analysis.report_path).unwrap(),
    )
    .unwrap();
    let mut source = original.request;
    let template: ObservedSamQueryUnitManifest =
        serde_json::from_slice(&fs::read(&source.manifest.path).unwrap()).unwrap();
    let mut raw = "@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:ctg_a\tLN:1000\n@SQ\tSN:ctg_b\tLN:1000\n@RG\tID:library\n".to_owned();
    let mut units = Vec::new();
    let mut originals = Vec::new();
    for (local_index, &original_index) in original_indices.iter().enumerate() {
        let qname = format!("equivalence_read_{original_index:04}");
        let unit_id = format!("bam_record_{:08}", original_index + 1);
        let alternate = original_index < 3;
        let reverse = original_index % 2 == 1;
        let cigar = if alternate { "10M10N20M" } else { "15M10N15M" };
        let mut aligned_sequence = "A".repeat(30);
        if alternate {
            aligned_sequence.replace_range(3..4, "C");
            aligned_sequence.replace_range(20..21, "G");
        }
        let mut aligned_quality = "I".repeat(30);
        aligned_quality.replace_range(3..4, "@"); // Phred31 at genomic103.
        aligned_quality.replace_range(20..21, "C"); // Phred34 at genomic130.
        raw.push_str(&format!(
            "{qname}\t{}\tctg_a\t101\t60\t{cigar}\t*\t0\t0\t{aligned_sequence}\t{aligned_quality}\tRG:Z:library\n",
            if reverse { 16 } else { 0 }
        ));
        let mut unit = template.units[0].clone();
        unit.unit_id = unit_id.clone();
        unit.rows[0].ordinal = (local_index + 1) as u64;
        unit.rows[0].qname = qname.clone();
        let search = unit.search.as_mut().unwrap();
        search.search_complete = false;
        search.views_complete = false;
        search.limitations =
            vec!["Only the invented supplied primary records are available".into()];
        units.push(unit);
        originals.push((
            unit_id,
            qname,
            if reverse {
                reverse_complement(&aligned_sequence)
            } else {
                aligned_sequence
            },
            if reverse {
                aligned_quality.chars().rev().collect::<String>()
            } else {
                aligned_quality
            },
        ));
    }
    let sam_path = directory.join(format!("equivalence {label}'s records.sam"));
    fs::write(&sam_path, raw).unwrap();
    source.sam = bound(&sam_path);
    source.conversion = crate::engine::ObservedSamConversion::BamConversion {
        bam_sha256: PARTITION_BAM_SHA.into(),
        producer: "Hand-crafted original-record equivalence fixture".into(),
        version: "synthetic-v1".into(),
        arguments: vec![],
        selection_scope:
            "Every named original record exactly once in this declared disjoint subset".into(),
    };
    source.manifest = save(
        &directory.join(format!("equivalence {label}'s units.json")),
        &ObservedSamQueryUnitManifest {
            schema: OBSERVED_SAM_QUERY_UNIT_MANIFEST_SCHEMA.into(),
            sam_sha256: source.sam.sha256.clone(),
            units,
        },
    );
    source.output_dir = directory
        .join(format!("equivalence {label}'s saved evidence"))
        .to_str()
        .unwrap()
        .into();
    let report = GentleEngine::new()
        .summarize_observed_transcript_structures_from_sam(source.clone(), &mut |_| true)
        .unwrap();
    let binding = GentleEngine::observed_sam_result_binding(&source, &report).unwrap();
    let (_, data_file) = crate::rna_read_ends::analysis_file(&binding).unwrap();
    let data: ObservedSamData =
        crate::rna_read_ends::bound_json(&data_file, &mut || Ok(())).unwrap();
    let analysis = data.analysis.unwrap();
    let reads = originals
        .into_iter()
        .map(
            |(unit_id, read_id, sequence, qualities)| CdnaOriginalBasecall {
                unit: analysis
                    .units
                    .iter()
                    .find(|u| u.unit.unit_id == unit_id)
                    .unwrap()
                    .unit
                    .clone(),
                read_id,
                sequence,
                qualities_phred33: Some(qualities),
            },
        )
        .collect();
    CdnaComparisonPartitionFiles {
        analysis: binding,
        basecalls: save(
            &directory.join(format!("equivalence {label}'s basecalls.json")),
            &CdnaBasecallDocument {
                schema: CDNA_BASECALL_SCHEMA.into(),
                analysis_context_sha256: analysis.context_sha256,
                orientation: "original_sequencing_5_to_3".into(),
                provenance:
                    "Literal synthetic original BAM records, identical across chunk layouts".into(),
                reads,
            },
        ),
    }
}

#[test]
fn cdna_comparison_partitioned_and_unpartitioned_original_records_are_semantically_identical() {
    let fixture = fixture(false, false);
    let whole = equivalence_evidence(&fixture, "whole", &[0, 1, 2, 3, 4, 5]);
    let parts = [
        equivalence_evidence(&fixture, "first", &[0, 1]),
        equivalence_evidence(&fixture, "second", &[2, 3]),
        equivalence_evidence(&fixture, "third", &[4, 5]),
    ];
    let original_calls = |parts: &[CdnaComparisonPartitionFiles]| {
        let mut calls = parts
            .iter()
            .flat_map(|p| {
                serde_json::from_slice::<CdnaBasecallDocument>(
                    &fs::read(&p.basecalls.path).unwrap(),
                )
                .unwrap()
                .reads
            })
            .map(|r| (r.read_id, r.sequence, r.qualities_phred33))
            .collect::<Vec<_>>();
        calls.sort();
        calls
    };
    assert_eq!(
        original_calls(std::slice::from_ref(&whole)),
        original_calls(&parts)
    );
    let mut unpartitioned = fixture.request.clone();
    unpartitioned.samples[0].analysis = whole.analysis;
    unpartitioned.samples[0].basecalls = whole.basecalls;
    unpartitioned.samples[0].partition_source_sha256 = Some(PARTITION_BAM_SHA.into());
    unpartitioned.policy.minimum_alternate_units = 3;
    unpartitioned.policy.minimum_haplotype_units = 3;
    let mut partitioned = unpartitioned.clone();
    partitioned.samples[0].analysis = parts[0].analysis.clone();
    partitioned.samples[0].basecalls = parts[0].basecalls.clone();
    partitioned.samples[0].additional_partitions = parts[1..].to_vec();
    let whole_report = run(unpartitioned);
    let chunk_report = run(partitioned);
    assert_eq!(whole_report.inputs.len(), 1);
    assert_eq!(chunk_report.inputs.len(), 3);
    // Input/request/context digests intentionally identify the different immutable containers.
    // Every biological result, normalized original base/quality and same-read skip set must agree.
    let semantic = |report: &CdnaFileComparisonReport| {
        let mut core = serde_json::to_value(&report.comparison).unwrap();
        core.as_object_mut().unwrap().remove("content_sha256");
        core.as_object_mut().unwrap().remove("request_sha256");
        for sample in core["samples"].as_array_mut().unwrap() {
            sample.as_object_mut().unwrap().remove("source_sha256");
            sample
                .as_object_mut()
                .unwrap()
                .remove("analysis_context_sha256");
        }
        core
    };
    assert_eq!(semantic(&whole_report), semantic(&chunk_report));
    let sample = &chunk_report.comparison.samples[0];
    assert_eq!(sample.total_units, 6);
    assert_eq!(sample.site_calls.len(), 12);
    assert_eq!(sample.site_observations.len(), 12);
    assert!(sample.site_observations.iter().all(
        |o| match o.observation.reference_position_0based {
            103 => o.observation.phred == Some(31),
            130 => o.observation.phred == Some(34),
            _ => false,
        }
    ));
    assert!(
        sample
            .site_observations
            .iter()
            .any(|o| o.observation.query_to_reference_strand == '-')
    );
    assert_eq!(sample.phase_blocks.len(), 1);
    assert_eq!(sample.phase_blocks[0].candidate_haplotype_count, Some(2));
    assert_eq!(sample.haplotype_junction_support.len(), 2);
    assert!(sample.haplotype_junction_support.iter().all(|j| {
        (j.junctions == vec![(110, 120)] || j.junctions == vec![(115, 125)])
            && j.possible_units == 3
            && j.exclusive_units == 3
    }));
}

#[test]
fn cdna_comparison_partitions_combine_support_before_one_sample_phase_run() {
    let mut fixture = partitioned_fixture();
    fixture.request.supplied_sites = vec![103, 125];
    let unscoped = run(fixture.request.clone());
    fixture.request.phase_marker_positions = Some(vec![125, 103]);
    let report = run(fixture.request);
    assert_eq!(report.request.phase_marker_positions, Some(vec![125, 103]));
    assert_eq!(
        report.comparison.phase_marker_positions,
        Some(vec![103, 125])
    );
    assert_eq!(
        serde_json::to_value(&report.comparison.sites).unwrap(),
        serde_json::to_value(&unscoped.comparison.sites).unwrap()
    );
    assert_eq!(
        serde_json::to_value(&report.comparison.samples[0].site_calls).unwrap(),
        serde_json::to_value(&unscoped.comparison.samples[0].site_calls).unwrap()
    );
    assert_eq!(
        serde_json::to_value(&report.comparison.samples[0].site_observations).unwrap(),
        serde_json::to_value(&unscoped.comparison.samples[0].site_observations).unwrap()
    );
    assert_eq!(report.inputs.len(), 4);
    assert_eq!(
        report
            .inputs
            .iter()
            .map(|p| p.partition_index)
            .collect::<Vec<_>>(),
        [0, 1, 2, 3]
    );
    assert!(report.inputs.iter().all(|p| p.sample_id == "synthetic_0"
        && p.units_with_incomplete_search_declaration == 1
        && p.units_with_incomplete_view_declaration == 1
        && p.units_with_unassessed_supplied_views == 0));
    assert_eq!(report.comparison.samples.len(), 1);
    assert_eq!(report.comparison.samples[0].total_units, 4);
    assert_eq!(report.comparison.sites.len(), 2);
    for site in &report.comparison.sites {
        assert_eq!(site.samples[0].total_units, 4);
        assert_eq!(site.samples[0].callable_units, 4);
        assert_eq!(site.samples[0].base_counts[&'A'], 2);
        assert_eq!(site.samples[0].candidate_alternate_bases.len(), 1);
    }
    // No partition alone meets the two-unit alternate or seed threshold.
    // The single combined sample has both two-marker templates, with two original units each.
    assert_eq!(report.comparison.samples[0].phase_blocks.len(), 1);
    let block = &report.comparison.samples[0].phase_blocks[0];
    assert_eq!(block.sites_0based, [103, 125]);
    assert_eq!(block.candidate_haplotype_count, Some(2));
    assert!(
        block
            .haplotypes
            .iter()
            .all(|h| h.exact_seed_units.len() == 2 && h.assigned_units == 2)
    );
    assert_eq!(block.assignments.len(), 4);
    assert_eq!(block.copy_number_status, "not_assessed");
}

#[test]
fn cdna_comparison_partitions_reject_reused_original_units_despite_distinct_sam_hashes() {
    for across_samples in [false, true] {
        let mut fixture = fixture(false, false);
        let first = partition(
            &fixture,
            0,
            "original_record_1",
            true,
            PARTITION_BAM_SHA,
            None,
        );
        let second = partition(
            &fixture,
            1,
            "original_record_1",
            true,
            PARTITION_BAM_SHA,
            None,
        );
        fixture.request.samples[0].analysis = first.analysis;
        fixture.request.samples[0].basecalls = first.basecalls;
        fixture.request.samples[0].partition_source_sha256 = Some(PARTITION_BAM_SHA.into());
        if across_samples {
            fixture.request.samples.push(CdnaComparisonSampleFiles {
                sample_id: "another_sample_label".into(),
                label: "Same original unit under another sample label".into(),
                analysis: second.analysis,
                basecalls: second.basecalls,
                additional_partitions: vec![],
                partition_source_sha256: Some(PARTITION_BAM_SHA.into()),
            });
        } else {
            fixture.request.samples[0]
                .additional_partitions
                .push(second);
        }
        let error = compare_cdna_files(fixture.request, &mut |_| true).unwrap_err();
        assert!(error.message.contains("partition_unit_reused"), "{error:?}");
    }
}

#[test]
fn cdna_comparison_partitions_require_one_matching_source_and_bounded_partition_count() {
    for which in 0..3 {
        let mut fixture = partitioned_fixture();
        match which {
            0 => fixture.request.samples[0].partition_source_sha256 = None,
            1 => {
                fixture.request.samples[0].partition_source_sha256 =
                    Some(format!("sha256:{}", "b".repeat(64)))
            }
            _ => {
                let duplicate = fixture.request.samples[0].additional_partitions[0].clone();
                fixture.request.samples[0].additional_partitions = vec![duplicate; 64];
            }
        }
        let error = compare_cdna_files(fixture.request, &mut |_| true).unwrap_err();
        let reason = if which == 1 {
            "partition_source_binding"
        } else {
            "partition_schema_or_source"
        };
        assert!(error.message.contains(reason), "case {which}: {error:?}");
    }
}

#[test]
fn cdna_comparison_partitions_preserve_one_annotation_authority_and_cumulative_budgets() {
    let mut fixture = partitioned_fixture();
    let mut annotation: PrimerPayloadAnnotation =
        serde_json::from_slice(&fs::read(&fixture.request.normalized_annotation.path).unwrap())
            .unwrap();
    annotation.release = "different invented annotation revision".into();
    let other = save(
        &fixture
            ._root
            .path()
            .canonicalize()
            .unwrap()
            .join("other annotation.json"),
        &annotation,
    );
    let extra = partition(
        &fixture,
        4,
        "original_record_5",
        false,
        PARTITION_BAM_SHA,
        Some(other),
    );
    fixture.request.samples[0].additional_partitions.push(extra);
    let error = compare_cdna_files(fixture.request, &mut |_| true).unwrap_err();
    assert!(
        error.message.contains("analysis_binding_or_duplicate"),
        "{error:?}"
    );

    for which in 0..3 {
        let mut fixture = partitioned_fixture();
        match which {
            0 => fixture.request.max_basecall_bases = 90,
            1 => fixture.request.budgets.max_units = 3,
            _ => {
                let file = &fixture.request.samples[0].additional_partitions[2].basecalls;
                let mut bytes = fs::read(&file.path).unwrap();
                bytes.push(b' ');
                fs::write(&file.path, bytes).unwrap();
            }
        }
        let error = compare_cdna_files(fixture.request, &mut |_| true).unwrap_err();
        if which == 0 {
            assert!(error.message.contains("basecall_base_budget"), "{error:?}");
        } else if which == 1 {
            assert!(error.message.contains("unit_budget"), "{error:?}");
        }
    }
}

#[test]
fn cdna_comparison_ingress_discovers_two_sample_snvs_and_annotates_without_forcing_two_alleles() {
    let fixture = fixture(false, true);
    let report = run(fixture.request.clone());
    assert_eq!(report.schema, CDNA_FILE_COMPARISON_SCHEMA);
    assert_eq!(report.inputs.len(), 2);
    assert!(report.inputs.iter().all(|s| s.source_eof_verified));
    assert_eq!(
        report
            .comparison
            .sites
            .iter()
            .map(|s| s.position_0based)
            .collect::<Vec<_>>(),
        [103, 125]
    );
    let first = &report.comparison.sites[0];
    assert_eq!(first.reference_base, 'A');
    assert_eq!(first.candidate_alternate_bases, ['C']);
    assert_eq!(first.samples[0].base_counts.get(&'C'), Some(&1));
    assert_eq!(first.samples[1].base_counts.get(&'A'), Some(&1));
    assert_eq!(first.transcripts[0].cdna_position_1based, Some(4));
    assert_eq!(first.transcripts[0].transcript_alternate_bases, ['C']);
    assert_eq!(report.comparison.samples.len(), 2);
    assert!(
        report
            .comparison
            .samples
            .iter()
            .flat_map(|s| &s.phase_blocks)
            .all(|b| b.haplotypes.len() <= 1)
    );
    assert_eq!(report.content_sha256, run(fixture.request).content_sha256);
}

#[test]
fn cdna_comparison_ingress_reverse_original_bases_and_qualities_use_shared_projection() {
    let mut fixture = fixture(true, false);
    fixture.request.supplied_sites = vec![103, 125];
    rewrite_calls(&mut fixture.request, |d| {
        // Reference-forward query position 15 maps to original reverse-query 14.
        d.reads[0]
            .qualities_phred33
            .as_mut()
            .unwrap()
            .replace_range(14..15, "&");
    });
    let report = run(fixture.request);
    let sample = &report.comparison.samples[0];
    let high = sample
        .site_calls
        .iter()
        .find(|s| s.position_0based == 103)
        .unwrap();
    let low = sample
        .site_calls
        .iter()
        .find(|s| s.position_0based == 125)
        .unwrap();
    assert_eq!(high.base, Some('C'));
    assert_eq!(high.minimum_phred, Some(40));
    assert_eq!(high.status, core::CdnaCallStatus::AlternateSupported);
    assert_eq!(low.status, core::CdnaCallStatus::Unassessed);
    let annotation = &report
        .comparison
        .sites
        .iter()
        .find(|s| s.position_0based == 103)
        .unwrap()
        .transcripts[0];
    assert_eq!(annotation.cdna_position_1based, Some(27));
    assert_eq!(annotation.transcript_reference_base, Some('T'));
    assert_eq!(annotation.transcript_alternate_bases, ['G']);
}

#[test]
fn cdna_comparison_ingress_reverse_iupac_original_base_remains_unassessed() {
    let mut fixture = fixture(true, false);
    fixture.request.supplied_sites = vec![103];
    rewrite_calls(&mut fixture.request, |d| {
        // Original reverse-query 26 projects to genomic 103; R complements to Y.
        d.reads[0].sequence.replace_range(26..27, "R");
    });
    let report = run(fixture.request);
    let site = report
        .comparison
        .sites
        .iter()
        .find(|s| s.position_0based == 103)
        .unwrap();
    assert!(site.candidate_alternate_bases.is_empty());
    assert_eq!(site.samples[0].callable_units, 0);
    assert_eq!(site.samples[0].unassessed_units, 1);
    let call = report.comparison.samples[0]
        .site_calls
        .iter()
        .find(|s| s.position_0based == 103)
        .unwrap();
    assert_eq!(call.status, core::CdnaCallStatus::Unassessed);
    assert_eq!(call.base, None);
}

#[test]
fn cdna_comparison_ingress_missing_quality_or_base_evidence_never_becomes_reference_support() {
    for which in 0..3 {
        let mut fixture = fixture(false, false);
        fixture.request.supplied_sites = vec![103];
        rewrite_calls(&mut fixture.request, |d| match which {
            0 => d.reads[0].qualities_phred33 = None,
            1 => d.reads.clear(),
            _ => d.reads[0].qualities_phred33 = Some("!".repeat(30)),
        });
        let report = run(fixture.request);
        let counts = &report.comparison.sites[0].samples[0];
        assert_eq!(counts.total_units, 1);
        assert_eq!(counts.callable_units, 0);
        assert_eq!(counts.unassessed_units, 1);
        assert_eq!(
            report.comparison.samples[0].site_calls[0].status,
            core::CdnaCallStatus::Unassessed
        );
    }
}

#[test]
fn cdna_comparison_ingress_rejects_rehashed_context_axis_name_quality_and_foreign_units() {
    for which in 0..8 {
        let mut fixture = fixture(false, false);
        rewrite_calls(&mut fixture.request, |d| match which {
            0 => d.analysis_context_sha256 = format!("sha256:{}", "f".repeat(64)),
            1 => d.orientation = "reference_forward".into(),
            2 => d.reads[0].read_id = "different_read".into(),
            3 => d.reads[0].sequence.push('A'),
            4 => d.reads[0].qualities_phred33 = Some("I".repeat(29)),
            5 => d.reads[0]
                .qualities_phred33
                .as_mut()
                .unwrap()
                .replace_range(0..1, " "),
            6 => d.reads[0].unit.unit_id = "foreign_unit".into(),
            _ => d.reads.push(d.reads[0].clone()),
        });
        assert!(
            compare_cdna_files(fixture.request, &mut |_| true).is_err(),
            "case {which}"
        );
    }
}

#[test]
fn cdna_comparison_ingress_rejects_raw_tampering_in_basecalls_reference_and_saved_sidecars() {
    for which in 0..3 {
        let fixture = fixture(false, false);
        let path = match which {
            0 => PathBuf::from(&fixture.request.samples[0].basecalls.path),
            1 => PathBuf::from(&fixture.request.reference.locus_fasta.path),
            _ => Path::new(&fixture.request.samples[0].analysis.report_path)
                .parent()
                .unwrap()
                .join("artifacts/analysis.json"),
        };
        let mut raw = fs::read(&path).unwrap();
        raw.push(b' ');
        fs::write(path, raw).unwrap();
        assert!(
            compare_cdna_files(fixture.request, &mut |_| true).is_err(),
            "case {which}"
        );
    }
}

#[test]
fn cdna_comparison_ingress_preserves_exact_reference_source_gene_and_transcript_binding() {
    for which in 0..5 {
        let mut fixture = fixture(false, false);
        match which {
            0 => {
                fixture.request.reference.genomic_reference_sha256 =
                    format!("sha256:{}", "f".repeat(64))
            }
            1 => fixture.request.reference.start_0based = 101,
            2 => fixture.request.gene_id = "missing_gene".into(),
            3 => fixture.request.source_id = "other_source".into(),
            _ => fixture.request.reference.contig = "other_contig".into(),
        }
        assert!(
            compare_cdna_files(fixture.request, &mut |_| true).is_err(),
            "case {which}"
        );
    }
}

#[test]
fn cdna_comparison_ingress_rejects_duplicate_sample_or_analysis_instead_of_multiplying_reads() {
    for distinct_name in [false, true] {
        let mut fixture = fixture(false, false);
        let mut duplicate = fixture.request.samples[0].clone();
        if distinct_name {
            duplicate.sample_id = "another_label_for_same_analysis".into();
        }
        fixture.request.samples.push(duplicate);
        assert!(compare_cdna_files(fixture.request, &mut |_| true).is_err());
    }
}

#[test]
fn cdna_comparison_ingress_distinguishes_assessed_supplied_views_from_complete_declarations() {
    let mut fixture = fixture(false, false);
    let directory = fixture._root.path().canonicalize().unwrap();
    let original = &fixture.request.samples[0];
    let report: crate::engine::ObservedSamReport =
        serde_json::from_slice(&fs::read(&original.analysis.report_path).unwrap()).unwrap();
    let mut source = report.request;
    let mut manifest: ObservedSamQueryUnitManifest =
        serde_json::from_slice(&fs::read(&source.manifest.path).unwrap()).unwrap();
    let declaration = manifest.units[0].search.as_mut().unwrap();
    declaration.search_complete = false;
    declaration.views_complete = false;
    declaration.limitations = vec!["Only the supplied primary view is available".into()];
    source.manifest = save(&directory.join("primary scoped units.json"), &manifest);
    source.output_dir = directory
        .join("primary scoped analysis")
        .to_str()
        .unwrap()
        .into();
    let mut calls: CdnaBasecallDocument =
        serde_json::from_slice(&fs::read(&original.basecalls.path).unwrap()).unwrap();
    let report = GentleEngine::new()
        .summarize_observed_transcript_structures_from_sam(source.clone(), &mut |_| true)
        .unwrap();
    let binding = GentleEngine::observed_sam_result_binding(&source, &report).unwrap();
    let (_, data_file) = crate::rna_read_ends::analysis_file(&binding).unwrap();
    let data: ObservedSamData =
        crate::rna_read_ends::bound_json(&data_file, &mut || Ok(())).unwrap();
    calls.analysis_context_sha256 = data.analysis.unwrap().context_sha256;
    fixture.request.samples[0].analysis = binding;
    fixture.request.samples[0].basecalls =
        save(&directory.join("primary scoped bases.json"), &calls);
    let report = run(fixture.request);
    let scope = &report.inputs[0];
    assert_eq!(
        scope.inference_scope,
        "supplied_eligible_mapping_hypotheses"
    );
    assert_eq!(scope.units_with_incomplete_search_declaration, 1);
    assert_eq!(scope.units_with_incomplete_view_declaration, 1);
    assert_eq!(scope.units_with_unassessed_supplied_views, 0);
    assert_eq!(report.comparison.sites[0].samples[0].callable_units, 1);
    assert_eq!(
        report.comparison.sites[0].samples[0].base_counts.get(&'C'),
        Some(&1)
    );
}

#[test]
fn cdna_comparison_ingress_rejects_same_declared_unit_in_independently_saved_reports() {
    let mut fixture = fixture(false, false);
    let directory = fixture._root.path().canonicalize().unwrap();
    let first = &fixture.request.samples[0];
    let report: crate::engine::ObservedSamReport =
        serde_json::from_slice(&fs::read(&first.analysis.report_path).unwrap()).unwrap();
    let mut source = report.request;
    let mut manifest: ObservedSamQueryUnitManifest =
        serde_json::from_slice(&fs::read(&source.manifest.path).unwrap()).unwrap();
    // Different qualification settings change the analysis context, not the physical SAM unit.
    source.required_flank_bases += 1;
    manifest.units[0].search.as_mut().unwrap().producer = "Other synthetic declaration".into();
    source.manifest = save(&directory.join("other producer's units.json"), &manifest);
    source.output_dir = directory
        .join("independently saved same unit")
        .to_str()
        .unwrap()
        .into();
    let report = GentleEngine::new()
        .summarize_observed_transcript_structures_from_sam(source.clone(), &mut |_| true)
        .unwrap();
    let binding = GentleEngine::observed_sam_result_binding(&source, &report).unwrap();
    let (_, data_file) = crate::rna_read_ends::analysis_file(&binding).unwrap();
    let data: ObservedSamData =
        crate::rna_read_ends::bound_json(&data_file, &mut || Ok(())).unwrap();
    let analysis = data.analysis.unwrap();
    let mut calls: CdnaBasecallDocument =
        serde_json::from_slice(&fs::read(&first.basecalls.path).unwrap()).unwrap();
    assert_eq!(calls.reads[0].unit, analysis.units[0].unit);
    assert_ne!(calls.analysis_context_sha256, analysis.context_sha256);
    calls.analysis_context_sha256 = analysis.context_sha256;
    fixture.request.samples.push(CdnaComparisonSampleFiles {
        sample_id: "relabelled_same_synthetic_unit".into(),
        label: "A different label for an already counted unit".into(),
        analysis: binding,
        basecalls: save(&directory.join("relabelled basecalls.json"), &calls),
        additional_partitions: vec![],
        partition_source_sha256: None,
    });
    let error = compare_cdna_files(fixture.request, &mut |_| true).unwrap_err();
    assert!(error.message.contains("unit_reused_across_samples"));
}

#[test]
fn cdna_comparison_ingress_saved_geometry_does_not_reopen_original_sam_manifest_or_dictionary() {
    let fixture = fixture(false, true);
    for path in &fixture.original_inputs {
        if path.exists() {
            fs::remove_file(path).unwrap();
        }
    }
    let report = run(fixture.request);
    assert_eq!(report.comparison.sites.len(), 2);
    assert_eq!(report.comparison.samples.len(), 2);
}

#[test]
fn cdna_comparison_ingress_budgets_schema_and_cancellation_fail_explicitly() {
    for which in 0..5 {
        let mut fixture = fixture(false, false);
        match which {
            0 => fixture.request.max_input_bytes = 1,
            1 => fixture.request.max_basecall_bases = 1,
            2 => fixture.request.budgets.max_comparisons = 1,
            3 => fixture.request.schema = "unknown_schema".into(),
            _ => fixture.request.budgets.max_output_bytes = 1,
        }
        assert!(
            compare_cdna_files(fixture.request, &mut |_| true).is_err(),
            "case {which}"
        );
    }
    let fixture = fixture(false, false);
    let error = compare_cdna_files(fixture.request.clone(), &mut |_| false).unwrap_err();
    assert!(error.message.contains("cancelled"));
    let mut polls = 0;
    assert!(
        compare_cdna_files(fixture.request.clone(), &mut |_| {
            polls += 1;
            polls < 4
        })
        .is_err()
    );
    assert!(polls >= 4);
    let mut value = serde_json::to_value(fixture.request).unwrap();
    value["impute_uncovered_reference"] = serde_json::json!(true);
    assert!(serde_json::from_value::<CdnaFileComparisonRequest>(value).is_err());
}

#[test]
fn cdna_comparison_ingress_typed_inline_and_quoted_file_routes_preserve_content_and_shell_history()
{
    let fixture = fixture(false, true);
    let request = fixture.request;
    let mut engine = GentleEngine::new();
    let before = serde_json::to_value(engine.snapshot()).unwrap();
    let revision = engine.mutation_revision();
    let undo = engine.undo_available();
    let typed = engine
        .apply(Operation::CompareCdnaAlleles {
            request: Box::new(request.clone()),
        })
        .unwrap();
    let expected = serde_json::to_value(typed.cdna_allele_comparison.unwrap()).unwrap();
    assert_eq!(serde_json::to_value(engine.snapshot()).unwrap(), before);
    assert_eq!(engine.mutation_revision(), revision);
    assert_eq!(engine.undo_available(), undo);
    let journal = engine.journal_len();
    let path = fixture._root.path().join("comparison request's α.json");
    fs::write(&path, serde_json::to_vec_pretty(&request).unwrap()).unwrap();
    for payload in [
        serde_json::to_string(&request).unwrap(),
        format!("@{}", path.to_str().unwrap()),
    ] {
        let parsed = parse_shell_line(&format!(
            "rna-reads allele-hash-screen --compare-cdna {}",
            quote_shell_arg(&payload)
        ))
        .unwrap();
        assert!(!parsed.is_state_mutating());
        let actual = execute_shell_command(&mut engine, &parsed).unwrap();
        assert!(!actual.state_changed);
        assert_eq!(actual.output["result"]["cdna_allele_comparison"], expected);
        assert_eq!(engine.journal_len(), journal);
        assert_eq!(serde_json::to_value(engine.snapshot()).unwrap(), before);
    }
}
