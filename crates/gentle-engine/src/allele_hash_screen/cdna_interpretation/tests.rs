//! Invented two-exon loci and read units, recreated in memory; no biological sample fixtures.
//! Translation fixtures read the existing standard-code asset rather than a second code table.

use super::*;
use ring::digest::{SHA256, digest};

fn code() -> BTreeMap<String, char> {
    let tables: serde_json::Value =
        serde_json::from_str(include_str!("../../../../../assets/codon_tables.json")).unwrap();
    let table = tables
        .as_array()
        .unwrap()
        .iter()
        .find(|t| t["id"] == 1)
        .unwrap()["sequence"]
        .as_str()
        .unwrap()
        .as_bytes();
    let mut out = BTreeMap::new();
    for (i, a) in b"TCAG".iter().enumerate() {
        for (j, b) in b"TCAG".iter().enumerate() {
            for (k, c) in b"TCAG".iter().enumerate() {
                let aa = table[i * 16 + j * 4 + k] as char;
                out.insert(
                    String::from_utf8(vec![*a, *b, *c]).unwrap(),
                    if aa == '|' { '*' } else { aa },
                );
            }
        }
    }
    out
}

fn fixture(reverse: bool) -> GeneInput {
    let transcript = "CCATGTCCTAAGG";
    let mut genomic = format!("{}{}{}", &transcript[..6], "N".repeat(14), &transcript[6..]);
    let mut exons = vec![(100, 106), (120, 127)];
    let mut positions = vec![120, 121];
    if reverse {
        genomic = genomic.chars().rev().map(|c| complement(c, '-')).collect();
        exons = exons
            .iter()
            .rev()
            .map(|&(a, z)| (227 - z, 227 - a))
            .collect();
        positions = positions.into_iter().map(|p| 226 - p).collect();
    }
    let strand = if reverse { '-' } else { '+' };
    let mut units = Vec::new();
    for (i, alleles) in [['C', 'C'], ['A', 'C'], ['C', 'A'], ['A', 'A']]
        .into_iter()
        .enumerate()
    {
        let mut bases = Vec::new();
        for &(a, z) in &exons {
            for p in a..z {
                let mut observed = genomic.as_bytes()[p - 100] as char;
                if let Some(index) = positions.iter().position(|v| *v == p) {
                    observed = complement(alleles[index], strand);
                }
                bases.push(evidence::CdnaBaseObservation {
                    reference_position_0based: p,
                    query_position_0based: p - 100,
                    observed_base: observed,
                    original_base: complement(observed, strand),
                    query_to_reference_strand: strand,
                    phred: Some(40),
                    mapq: Some(60),
                });
            }
        }
        units.push(evidence::CdnaUnit {
            unit_id: format!("synthetic-unit-{i}"),
            all_query_views_assessed: true,
            search_declared_complete: true,
            hypotheses: vec![evidence::CdnaHypothesis {
                hypothesis_id: "0".into(),
                eligible: true,
                bases,
                junctions: vec![(exons[0].1, exons[1].0)],
            }],
        });
    }
    let hash = format!(
        "sha256:{}",
        digest(&SHA256, genomic.as_bytes())
            .as_ref()
            .iter()
            .map(|b| format!("{b:02x}"))
            .collect::<String>()
    );
    let t = evidence::CdnaTranscript {
        transcript_id: "synthetic-transcript".into(),
        strand,
        exons: exons.clone(),
    };
    let mut policy = evidence::CdnaComparisonPolicy::default();
    policy.minimum_haplotype_units = 1;
    policy.minimum_alternate_units = 1;
    let sample = evidence::CdnaSample {
        sample_id: "synthetic-sample".into(),
        label: "invented".into(),
        source_sha256: hash.clone(),
        analysis_context_sha256: hash.clone(),
        units,
    };
    let request = evidence::CdnaComparisonRequest {
        schema: evidence::CDNA_COMPARISON_REQUEST_SCHEMA.into(),
        gene: "invented-gene".into(),
        reference: evidence::CdnaReference {
            assembly: "synthetic-assembly".into(),
            contig: "chrSynthetic".into(),
            start_0based: 100,
            sequence: genomic.clone(),
            sha256: hash,
        },
        transcripts: vec![t.clone()],
        samples: vec![sample.clone()],
        supplied_sites: positions,
        phase_marker_positions: None,
        policy,
        budgets: Default::default(),
    };
    let report = evidence::compare_cdna_alleles(request, &mut |_, _| true).unwrap();
    let geometry = sample
        .units
        .iter()
        .flat_map(|u| {
            u.hypotheses.iter().map(|h| {
                (
                    (
                        sample.sample_id.clone(),
                        u.unit_id.clone(),
                        h.hypothesis_id.clone(),
                    ),
                    AlignmentGeometry {
                        blocks: exons
                            .iter()
                            .map(|&(a, z)| ("chrSynthetic".into(), a, z))
                            .collect(),
                        junctions: vec![("chrSynthetic".into(), exons[0].1, exons[1].0)],
                    },
                )
            })
        })
        .collect();
    GeneInput {
        source_id: "synthetic-source".into(),
        gene_id: "synthetic-gene".into(),
        report,
        models: vec![CodingModel {
            transcript: t,
            accession: Some("SYN_T.2".into()),
            protein_accession: Some("SYN_P.3".into()),
            sequence: transcript.into(),
            cds: Some((2, 11)),
            translation_status: "complete_verified".into(),
            translation_reasons: Vec::new(),
        }],
        genomic_sequence: genomic,
        linkage: Some(LinkageEvidence {
            samples: vec![sample],
            geometry,
        }),
    }
}

fn run(input: GeneInput) -> GeneReport {
    interpret(
        &[input],
        &Policy::default(),
        &code(),
        None,
        &Budgets::default(),
        &mut || true,
    )
    .unwrap()
    .0
    .remove(0)
}

#[test]
fn cdna_interpretation_all_four_same_unit_codon_patterns_and_mirrored_coordinates() {
    for reverse in [false, true] {
        let report = run(fixture(reverse));
        assert_eq!(report.codon_combinations.len(), 4);
        assert_eq!(report.codon_groups[0].complete_units, 4);
        assert!(report.codon_groups[0].incomplete_units.is_empty());
        let texts = report
            .codon_combinations
            .iter()
            .filter_map(|c| c.consequence.protein_hgvs.clone())
            .collect::<BTreeSet<_>>();
        assert_eq!(
            texts,
            BTreeSet::from([
                "SYN_P.3:p.(Ser2=)".into(),
                "SYN_P.3:p.(Ser2Tyr)".into(),
                "SYN_P.3:p.(Ser2Ter)".into()
            ])
        );
        let stop = report
            .codon_combinations
            .iter()
            .find(|c| {
                c.consequence
                    .consequence_terms
                    .contains(&"stop_gained".into())
            })
            .unwrap();
        assert_eq!(
            stop.consequence.nucleotide_hgvs.as_deref(),
            Some("SYN_T.2:c.5_6delinsAA")
        );
        assert_eq!(stop.unit_ids, vec!["synthetic-unit-3"]);
        assert_eq!(report.sites.len(), 2);
        assert!(
            report
                .sites
                .iter()
                .all(|r| r.alternate_units == 2 && r.callable_units == 4)
        );
        assert_eq!(
            report
                .sites
                .iter()
                .filter_map(|r| r.consequence.nucleotide_hgvs.as_deref())
                .collect::<BTreeSet<_>>(),
            BTreeSet::from(["SYN_T.2:c.5C>A", "SYN_T.2:c.6C>A"])
        );
    }
}

#[test]
fn cdna_interpretation_hgvs_uses_transcript_reference_across_introns_and_one_base_gap() {
    let input = fixture(false);
    let model = &input.models[0];
    let changes = vec![
        change(&input, model, 105, 'C').unwrap(),
        change(&input, model, 121, 'A').unwrap(),
    ];
    let mut progress = || true;
    let mut meter = Meter {
        work: 0,
        limit: 10000,
        progress: &mut progress,
        remaining_output: 10000,
    };
    let result = hgvs::consequence(model, &changes, &code(), &mut meter).unwrap();
    assert_eq!(
        result.nucleotide_hgvs.as_deref(),
        Some("SYN_T.2:c.4_6delinsCCA")
    );
    assert_eq!(result.protein_hgvs.as_deref(), Some("SYN_P.3:p.(Ser2Pro)"));
}

#[test]
fn cdna_interpretation_disjoint_calls_never_make_a_codon_combination() {
    let mut input = fixture(false);
    for (i, call) in input.report.samples[0].site_calls.iter_mut().enumerate() {
        if (i % 2 == 0 && call.position_0based == 120)
            || (i % 2 == 1 && call.position_0based == 121)
        {
            call.base = None;
            call.status = evidence::CdnaCallStatus::Unassessed;
        }
    }
    // Explicitly retain only one assessed position per unit, independent of report ordering.
    for call in &mut input.report.samples[0].site_calls {
        let unit = call
            .unit_id
            .rsplit('-')
            .next()
            .unwrap()
            .parse::<usize>()
            .unwrap();
        if call.position_0based != if unit % 2 == 0 { 120 } else { 121 } {
            call.base = None;
            call.status = evidence::CdnaCallStatus::Unassessed;
        }
    }
    let report = run(input);
    assert!(report.codon_combinations.is_empty());
    assert_eq!(report.codon_groups[0].incomplete_units.len(), 4);
    assert!(
        report
            .sites
            .iter()
            .all(|r| r.codon_context == "codon_context_unresolved")
    );
}

#[test]
fn cdna_interpretation_complete_reference_patterns_cannot_resolve_partial_alternate_context() {
    let mut input = fixture(false);
    input.linkage = None;
    for call in &mut input.report.samples[0].site_calls {
        if call.position_0based == 121
            && matches!(
                call.unit_id.as_str(),
                "synthetic-unit-1" | "synthetic-unit-3"
            )
        {
            call.base = None;
            call.status = evidence::CdnaCallStatus::Unassessed;
        }
    }
    let report = run(input);
    assert_eq!(report.codon_combinations.len(), 2);
    assert_eq!(report.codon_groups[0].incomplete_units.len(), 2);
    let unlinked = report
        .sites
        .iter()
        .find(|r| r.genomic_position_0based == 120)
        .unwrap();
    assert_eq!(unlinked.codon_context, "codon_context_unresolved");
    let linked = report
        .sites
        .iter()
        .find(|r| r.genomic_position_0based == 121)
        .unwrap();
    assert_eq!(linked.codon_context, "see_codon_haplotype");
}

#[test]
fn cdna_interpretation_utr_noncoding_not_present_and_unavailable_models() {
    let input = fixture(false);
    let mut model = input.models[0].clone();
    let mut progress = || true;
    let mut meter = Meter {
        work: 0,
        limit: 10000,
        progress: &mut progress,
        remaining_output: 10000,
    };
    for (p, alt, expected) in [(100, 'A', "SYN_T.2:c.-2C>A"), (126, 'A', "SYN_T.2:c.*2G>A")] {
        let c = change(&input, &model, p, alt).unwrap();
        assert_eq!(
            hgvs::consequence(&model, &[c], &code(), &mut meter)
                .unwrap()
                .nucleotide_hgvs
                .as_deref(),
            Some(expected)
        );
    }
    model.cds = None;
    model.translation_status = "noncoding".into();
    assert_eq!(
        hgvs::consequence(
            &model,
            &[change(&input, &model, 120, 'A').unwrap()],
            &code(),
            &mut meter
        )
        .unwrap()
        .nucleotide_hgvs
        .as_deref(),
        Some("SYN_T.2:n.7C>A")
    );
    assert_eq!(position(&model, 110), None);
    let absent = hgvs::consequence(&model, &[], &code(), &mut meter).unwrap();
    assert_eq!(absent.region, "not_present");
    assert!(absent.nucleotide_hgvs.is_none());
    for status in [
        "incomplete_start",
        "incomplete_end",
        "phase_inconsistent",
        "length_not_multiple_of_3",
        "internal_stop",
        "supplied_protein_mismatch",
    ] {
        let mut model = input.models[0].clone();
        model.translation_status = status.into();
        let result = hgvs::consequence(
            &model,
            &[change(&input, &model, 120, 'A').unwrap()],
            &code(),
            &mut meter,
        )
        .unwrap();
        assert!(!result.assessable);
        assert!(result.nucleotide_hgvs.is_none());
        assert!(result.protein_hgvs.is_none());
    }
    let mut model = input.models[0].clone();
    model.accession = None;
    let result = hgvs::consequence(
        &model,
        &[change(&input, &model, 120, 'A').unwrap()],
        &code(),
        &mut meter,
    )
    .unwrap();
    assert!(result.assessable);
    assert!(result.protein_hgvs.is_none());
    assert!(result.nucleotide_hgvs.is_none());
}

#[test]
fn cdna_interpretation_start_stop_retained_lost_and_bounded_extension() {
    let input = fixture(false);
    let model = &input.models[0];
    let code = code();
    let mut progress = || true;
    let mut meter = Meter {
        work: 0,
        limit: 10000,
        progress: &mut progress,
        remaining_output: 10000,
    };
    for (p, alt, text) in [
        (102, 'C', "SYN_P.3:p.(Met1?)"),
        (124, 'G', "SYN_P.3:p.(Ter3=)"),
        (122, 'C', "SYN_P.3:p.(Ter3GlnextTer?)"),
    ] {
        let result = hgvs::consequence(
            model,
            &[change(&input, model, p, alt).unwrap()],
            &code,
            &mut meter,
        )
        .unwrap();
        assert_eq!(result.protein_hgvs.as_deref(), Some(text));
    }
    let mut extended = model.clone();
    extended.sequence.push_str("ATAA");
    let result = hgvs::consequence(
        &extended,
        &[change(&input, model, 122, 'C').unwrap()],
        &code,
        &mut meter,
    )
    .unwrap();
    // Added sequence is Gln (former stop), Gly (GGA), Ter (TAA): new stop at 3.
    assert_eq!(
        result.protein_hgvs.as_deref(),
        Some("SYN_P.3:p.(Ter3GlnextTer3)")
    );
}

#[test]
fn cdna_interpretation_linkage_distinguishes_compatible_discriminating_and_intron_conflicts() {
    let mut input = fixture(false);
    let mut competitor = input.models[0].clone();
    competitor.transcript.transcript_id = "identical-alternative".into();
    input.models.push(competitor.clone());
    input.report.transcripts.push(competitor.transcript.clone());
    for site in &mut input.report.sites {
        let mut annotation = site.transcripts[0].clone();
        annotation.transcript_id = "identical-alternative".into();
        site.transcripts.push(annotation);
    }
    let report = run(input.clone());
    assert!(
        report
            .sites
            .iter()
            .all(|r| r.linkage.as_ref().unwrap().status == "same_unit_compatible_only")
    );
    input.models[1].transcript.exons = vec![(100, 127)];
    input.models[1].sequence = input.genomic_sequence.clone();
    input.models[1].cds = None;
    input.models[1].translation_status = "noncoding".into();
    input.report.transcripts[1] = input.models[1].transcript.clone();
    for site in &mut input.report.sites {
        site.transcripts[1].cdna_position_1based = Some(site.position_0based - 100 + 1);
    }
    let report = run(input.clone());
    let row = report
        .sites
        .iter()
        .find(|r| r.transcript_id == "synthetic-transcript")
        .unwrap();
    assert_eq!(
        row.linkage.as_ref().unwrap().status,
        "same_unit_discriminating_support"
    );
    assert_eq!(row.linkage.as_ref().unwrap().unique_model_units.len(), 2);
    for geometry in input.linkage.as_mut().unwrap().geometry.values_mut() {
        geometry.blocks = vec![("chrSynthetic".into(), 100, 127)];
    }
    let report = run(input);
    assert!(
        report
            .sites
            .iter()
            .filter(|r| r.transcript_id == "synthetic-transcript")
            .all(|r| r.linkage.as_ref().unwrap().compatible_units.is_empty())
    );
}

#[test]
fn cdna_interpretation_never_combines_allele_and_junction_across_hypotheses() {
    let mut input = fixture(false);
    let link = input.linkage.as_mut().unwrap();
    for unit in &mut link.samples[0].units {
        let mut alternative = unit.hypotheses[0].clone();
        alternative.hypothesis_id = "other".into();
        alternative.bases.clear();
        let key: (String, String, String) =
            ("synthetic-sample".into(), unit.unit_id.clone(), "0".into());
        link.geometry.insert(
            (key.0.clone(), key.1.clone(), "other".into()),
            link.geometry[&key].clone(),
        );
        unit.hypotheses.push(alternative);
    }
    assert!(
        run(input)
            .sites
            .iter()
            .all(|r| r.linkage.as_ref().unwrap().compatible_units.is_empty())
    );
}

#[test]
fn cdna_interpretation_frozen_known_labels_match_exact_supported_composite_only() {
    let input = fixture(false);
    let manifest = KnownVariantManifest {
        schema: "gentle.known_variant_manifest.v1".into(),
        source_name: "synthetic-source".into(),
        source_release: "invented-v1".into(),
        retrieved_utc: "2026-11-12T00:00:00Z".into(),
        license: "synthetic-redistributable".into(),
        attribution: "hand-crafted test only".into(),
        assembly: "synthetic-assembly".into(),
        entries: vec![KnownVariant {
            contig: "chrSynthetic".into(),
            position_1based: 121,
            r#ref: "CC".into(),
            alt: "AA".into(),
            identifier: "invented-composite".into(),
        }],
    };
    let out = interpret(
        &[input.clone()],
        &Policy::default(),
        &code(),
        Some(&manifest),
        &Budgets::default(),
        &mut || true,
    )
    .unwrap()
    .0;
    assert_eq!(
        out[0]
            .codon_combinations
            .iter()
            .filter(|c| !c.known_variants.is_empty())
            .count(),
        1
    );
    assert!(out[0].sites.iter().all(|r| r.known_variants.is_empty()));
    let mut invalid = manifest.clone();
    invalid.assembly = "wrong".into();
    assert!(
        interpret(
            &[input],
            &Policy::default(),
            &code(),
            Some(&invalid),
            &Budgets::default(),
            &mut || true
        )
        .is_err()
    );
    let mut value = serde_json::to_value(manifest).unwrap();
    value["clinical_significance"] = serde_json::json!("not allowed");
    assert!(serde_json::from_value::<KnownVariantManifest>(value).is_err());
}

#[test]
fn cdna_interpretation_inherited_support_disconnected_blocks_and_sample_boundaries() {
    let mut input = fixture(false);
    input.report.sites[0].samples[0]
        .candidate_alternate_bases
        .clear();
    input.report.samples[0].phase_blocks = vec![
        evidence::CdnaPhaseBlock {
            block_id: "a".into(),
            sites_0based: vec![120],
            status: "unresolved".into(),
            candidate_haplotype_count: None,
            excess_expressed_haplotype_diversity: false,
            copy_number_status: "not_assessed".into(),
            haplotypes: Vec::new(),
            assignments: Vec::new(),
        },
        evidence::CdnaPhaseBlock {
            block_id: "b".into(),
            sites_0based: vec![121],
            status: "unresolved".into(),
            candidate_haplotype_count: None,
            excess_expressed_haplotype_diversity: false,
            copy_number_status: "not_assessed".into(),
            haplotypes: Vec::new(),
            assignments: Vec::new(),
        },
    ];
    let report = run(input);
    assert_eq!(report.local_phase_relations[0].relation, "different_blocks");
    assert!(
        report
            .sites
            .iter()
            .filter(|r| r.genomic_position_0based == 120)
            .all(|r| r.support_status == "below_inherited_threshold")
    );
}

#[test]
fn cdna_interpretation_cancellation_budgets_policy_and_binding_refuse_without_partial_results() {
    let input = fixture(false);
    assert!(
        interpret(
            &[input.clone()],
            &Policy::default(),
            &code(),
            None,
            &Budgets::default(),
            &mut || false
        )
        .is_err()
    );
    let mut tiny = Budgets::default();
    tiny.max_work = 1;
    assert!(
        interpret(
            &[input.clone()],
            &Policy::default(),
            &code(),
            None,
            &tiny,
            &mut || true
        )
        .is_err()
    );
    tiny = Budgets::default();
    tiny.max_output_bytes = 1;
    assert!(
        interpret(
            &[input.clone()],
            &Policy::default(),
            &code(),
            None,
            &tiny,
            &mut || true
        )
        .is_err()
    );
    let mut invalid = input;
    invalid.models[0].sequence.replace_range(0..1, "G");
    assert!(
        interpret(
            &[invalid],
            &Policy::default(),
            &code(),
            None,
            &Budgets::default(),
            &mut || true
        )
        .is_err()
    );
}

#[test]
fn cdna_interpretation_low_depth_is_inherited_and_cannot_lower_the_composite_threshold() {
    let mut input = fixture(false);
    let mut policy = input.report.policy.clone();
    policy.minimum_alternate_units = 3;
    policy.minimum_haplotype_units = 3;
    let request = evidence::CdnaComparisonRequest {
        schema: evidence::CDNA_COMPARISON_REQUEST_SCHEMA.into(),
        gene: input.report.gene.clone(),
        reference: evidence::CdnaReference {
            assembly: input.report.reference.assembly.clone(),
            contig: input.report.reference.contig.clone(),
            start_0based: 100,
            sequence: input.genomic_sequence.clone(),
            sha256: input.report.reference.sequence_sha256.clone(),
        },
        transcripts: input.report.transcripts.clone(),
        samples: input.linkage.as_ref().unwrap().samples.clone(),
        supplied_sites: vec![120, 121],
        phase_marker_positions: None,
        policy,
        budgets: Default::default(),
    };
    input.report = evidence::compare_cdna_alleles(request, &mut |_, _| true).unwrap();
    let result = run(input.clone());
    assert_eq!(result.sites.len(), 2);
    assert!(result.sites.iter().all(|r| r.alternate_units == 2
        && r.callable_units == 4
        && r.support_status == "below_inherited_threshold"));
    assert!(
        result
            .codon_combinations
            .iter()
            .all(|c| !c.supported && c.minimum_support == 3)
    );
    let mut override_policy = Policy::default();
    override_policy.minimum_codon_haplotype_units = Some(1);
    assert!(
        interpret(
            &[input],
            &override_policy,
            &code(),
            None,
            &Budgets::default(),
            &mut || true
        )
        .is_err()
    );
}

#[test]
fn cdna_interpretation_adjacent_multi_codon_changes_keep_per_codon_predictions_only() {
    let input = fixture(false);
    let model = &input.models[0];
    let changes = vec![
        change(&input, model, 104, 'A').unwrap(),
        change(&input, model, 120, 'A').unwrap(),
    ];
    let mut progress = || true;
    let mut meter = Meter {
        work: 0,
        limit: 10000,
        progress: &mut progress,
        remaining_output: 10000,
    };
    let consequence = hgvs::consequence(model, &changes, &code(), &mut meter).unwrap();
    assert_eq!(
        consequence.nucleotide_hgvs.as_deref(),
        Some("SYN_T.2:c.[3G>A;5C>A]")
    );
    assert_eq!(consequence.coding_changes.len(), 2);
    assert!(consequence.protein_hgvs.is_none());
    assert!(
        consequence
            .reasons
            .contains(&"multi_codon_composite_not_rendered".into())
    );
}

#[test]
fn cdna_interpretation_samples_with_disjoint_partial_calls_never_share_compositions() {
    let mut input = fixture(false);
    input.linkage = None;
    let mut second = input.report.samples[0].clone();
    second.sample_id = "another-synthetic-sample".into();
    for call in &mut second.site_calls {
        call.unit_id = format!("other:{}", call.unit_id);
        if call.position_0based == 120 {
            call.base = None;
            call.status = evidence::CdnaCallStatus::Unassessed;
        }
    }
    for call in &mut input.report.samples[0].site_calls {
        if call.position_0based == 121 {
            call.base = None;
            call.status = evidence::CdnaCallStatus::Unassessed;
        }
    }
    input.report.samples.push(second);
    for site in &mut input.report.sites {
        let mut counts = site.samples[0].clone();
        counts.sample_id = "another-synthetic-sample".into();
        site.samples.push(counts);
    }
    let report = run(input);
    assert!(report.codon_combinations.is_empty());
    assert_eq!(report.codon_groups.len(), 2);
    assert!(
        report
            .codon_groups
            .iter()
            .all(|g| g.complete_units == 0 && g.incomplete_units.len() == 4)
    );
}

#[test]
fn cdna_interpretation_manifest_timestamp_and_input_indels_are_explicit() {
    assert!(valid_utc("2024-02-29T23:59:59Z"));
    for invalid in [
        "2023-02-29T00:00:00Z",
        "2026-13-01T00:00:00Z",
        "2026-01-01T24:00:00Z",
        "2026-01-01T00:60:00Z",
    ] {
        assert!(!valid_utc(invalid));
    }
    let input = fixture(false);
    let known = KnownVariantManifest {
        schema: "gentle.known_variant_manifest.v1".into(),
        source_name: "synthetic".into(),
        source_release: "v1".into(),
        retrieved_utc: "2026-11-12T00:00:00Z".into(),
        license: "synthetic".into(),
        attribution: "invented only".into(),
        assembly: input.report.reference.assembly.clone(),
        entries: vec![KnownVariant {
            contig: "chrSynthetic".into(),
            position_1based: 121,
            r#ref: "C".into(),
            alt: "CA".into(),
            identifier: "unrelated-indel".into(),
        }],
    };
    let report = interpret(
        &[input],
        &Policy::default(),
        &code(),
        Some(&known),
        &Budgets::default(),
        &mut || true,
    )
    .unwrap()
    .0;
    assert!(report[0].sites.iter().all(|s| s.known_variants.is_empty()));
    assert!(
        report[0]
            .codon_combinations
            .iter()
            .all(|s| s.known_variants.is_empty())
    );
}

#[test]
fn cdna_interpretation_batch_keeps_gene_identity_counts_statuses_and_policy_boundaries() {
    let a = fixture(false);
    let mut b = fixture(true);
    b.gene_id = "other-invented-gene".into();
    b.report.gene = "other label, not a selector".into();
    let (rows, _) = interpret(
        &[a.clone(), b.clone()],
        &Policy::default(),
        &code(),
        None,
        &Budgets::default(),
        &mut || true,
    )
    .unwrap();
    assert_eq!(rows.len(), 2);
    for (row, input) in rows.iter().zip([&a, &b]) {
        assert_eq!(row.gene_id, input.gene_id);
        assert_eq!(
            serde_json::to_value(&row.inherited_samples).unwrap(),
            serde_json::to_value(&input.report.samples).unwrap()
        );
        assert_eq!(row.inherited_warnings, input.report.warnings);
    }
    b.report.policy.minimum_base_quality += 1;
    assert!(
        interpret(
            &[a.clone(), b],
            &Policy::default(),
            &code(),
            None,
            &Budgets::default(),
            &mut || true
        )
        .unwrap_err()
        .contains("policy_mismatch")
    );
    assert!(
        interpret(
            &[a.clone(), a],
            &Policy::default(),
            &code(),
            None,
            &Budgets::default(),
            &mut || true
        )
        .unwrap_err()
        .contains("duplicate_gene")
    );
}

#[test]
fn cdna_interpretation_incomplete_query_views_and_foreign_contigs_cannot_identify_models() {
    for foreign in [false, true] {
        let mut input = fixture(false);
        let linkage = input.linkage.as_mut().unwrap();
        if foreign {
            for geometry in linkage.geometry.values_mut() {
                geometry.blocks.push(("another-contig".into(), 100, 101));
            }
        } else {
            for unit in &mut linkage.samples[0].units {
                unit.all_query_views_assessed = false;
            }
        }
        assert!(
            run(input).sites.iter().all(|r| r
                .linkage
                .as_ref()
                .unwrap()
                .compatible_units
                .is_empty())
        );
    }
}
