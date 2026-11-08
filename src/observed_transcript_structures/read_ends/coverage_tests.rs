//! Synthetic saved-projection controls; no BAM, files, private reads or scientific claims.

use super::super::tests::{fixture, project_read_ends};
use super::*;

fn bindings() -> InputBindings {
    let hash = format!("sha256:{}", "a".repeat(64));
    InputBindings {
        report_content_sha256: hash.clone(),
        data_sha256: hash.clone(),
        model_sha256: hash.clone(),
        projection_sha256: hash,
    }
}

fn run(analysis: &ObservedSamAnalysis, projection: &Projection) -> CoverageReport {
    summarize(
        &projection.source_id,
        "GENE_SYN_A",
        Some("tx_plus"),
        &[0.0, 0.5, 1.0],
        &[SampleInput {
            sample_id: "synthetic_sample",
            gene_id: "GENE_SYN_A",
            bindings: &bindings(),
            analysis,
            projection,
        }],
        Limits::default(),
    )
    .unwrap()
}

#[test]
fn read_end_coverage_keeps_noncoding_gene_units_and_distinct_denominators() {
    let (document, analysis) = fixture("missing_candidate_model");
    let projection = project_read_ends(&document, &analysis, Limits::default()).unwrap();
    let report = run(&analysis, &projection);
    let sample = &report.samples[0];
    assert_eq!(sample.declared_units, 1);
    assert_eq!(sample.assessed_units, 0);
    assert_eq!(sample.gene_assigned_units, 1);
    assert_eq!(sample.selected_transcript_compatible_units, 1);
    assert_eq!(sample.compatibility_counts[&Compatibility::Shared], 1);
    assert_eq!(sample.utr3.missing, 1);
    assert_eq!(sample.utr3.exact_quantiles[1], (0.5, None));
}

#[test]
fn read_end_coverage_exclusions_and_compatibility_classes_are_explicit() {
    let (document, analysis) = fixture("plus");
    let projection = project_read_ends(&document, &analysis, Limits::default()).unwrap();
    let baseline = &projection.units[0];
    assert_eq!(
        compatibility(baseline, "tx_plus"),
        Some(Compatibility::Unique)
    );
    assert_eq!(
        gene_assignment(baseline, true, "GENE_SYN_A").unwrap_err().0,
        ExclusionBucket::UnavailableUnit
    );
    let mut unit = baseline.clone();
    unit.hypotheses.clear();
    assert_eq!(
        gene_assignment(&unit, false, "GENE_SYN_A").unwrap_err().0,
        ExclusionBucket::NoOptimalHypothesis
    );
    unit = baseline.clone();
    unit.hypotheses[0].candidates.clear();
    assert_eq!(
        gene_assignment(&unit, false, "GENE_SYN_A").unwrap_err().0,
        ExclusionBucket::NoCompatibleTranscript
    );
    unit = baseline.clone();
    unit.annotation_gene_ids.push("OTHER_SYN_GENE".into());
    assert_eq!(
        gene_assignment(&unit, false, "GENE_SYN_A").unwrap_err().0,
        ExclusionBucket::MultiGene
    );
    unit = baseline.clone();
    let mut other = unit.hypotheses[0].clone();
    other.hypothesis_index = 1;
    other.candidates[0].transcript_id = "other_transcript".into();
    unit.hypotheses.push(other);
    assert_eq!(
        compatibility(&unit, "tx_plus"),
        Some(Compatibility::HypothesisDependent)
    );
}

#[test]
fn read_end_coverage_fraction_denominators_and_exact_quantiles_are_not_interchangeable() {
    let values = [
        Metric::Observed { value_bp: 0 },
        Metric::Observed { value_bp: 10 },
        Metric::RightCensored {
            lower_bound_bp: 900,
        },
        missing(Reason::StopCodonUnavailable),
    ];
    let summary = summary(values.iter().collect(), &[0.5]);
    assert_eq!(summary.exact_quantiles, [(0.5, Some(5.0))]);
    assert_eq!(summary.censored_among_observed_and_censored.denominator, 3);
    assert_eq!(summary.censored_among_all.denominator, 4);
    assert_eq!(summary.missing_by_reason[&Reason::StopCodonUnavailable], 1);
}

#[test]
fn read_end_coverage_alternatives_never_double_count_and_empty_is_not_zero_quantile() {
    let (document, mut analysis) = fixture("plus");
    let repeated = analysis.units[0].hypotheses[0].clone();
    analysis.units[0].hypotheses.push(repeated);
    let mut projection = project_read_ends(&document, &analysis, Limits::default()).unwrap();
    let report = run(&analysis, &projection);
    assert_eq!(report.samples[0].coverage[0].transcript_interval, (0, 25));
    assert_eq!(report.samples[0].coverage[0].maximum_conditional_depth, 1);
    assert_eq!(report.samples[0].coverage[0].minimum_conditional_depth, 1);
    projection.units.clear();
    analysis.units.clear();
    analysis.observations.clear();
    let empty = run(&analysis, &projection);
    assert_eq!(empty.samples[0].gene_assigned_units, 0);
    assert_eq!(empty.samples[0].coverage[0].minimum_fraction.value, None);
    assert_eq!(empty.samples[0].utr3.exact_quantiles[1], (0.5, None));
}

#[test]
fn read_end_coverage_duplicate_sample_gene_and_mismatched_bindings_refuse() {
    let (document, analysis) = fixture("plus");
    let projection = project_read_ends(&document, &analysis, Limits::default()).unwrap();
    let hash = bindings();
    let inputs = [
        SampleInput {
            sample_id: "s",
            gene_id: "GENE_SYN_A",
            bindings: &hash,
            analysis: &analysis,
            projection: &projection,
        },
        SampleInput {
            sample_id: "s",
            gene_id: "GENE_SYN_A",
            bindings: &hash,
            analysis: &analysis,
            projection: &projection,
        },
    ];
    assert_eq!(
        summarize(
            &projection.source_id,
            "GENE_SYN_A",
            Some("tx_plus"),
            &[0.5],
            &inputs,
            Limits::default()
        )
        .unwrap_err()
        .message,
        "rna_read_end_coverage_sample_gene_analysis"
    );
    let before = serde_json::to_vec(&run(&analysis, &projection)).unwrap();
    assert_eq!(
        serde_json::to_vec(&run(&analysis, &projection)).unwrap(),
        before
    );
    let mut changed = analysis.clone();
    changed.reference_gene.gene_id = "OTHER".into();
    assert!(
        summarize(
            &projection.source_id,
            "GENE_SYN_A",
            None,
            &[],
            &[SampleInput {
                sample_id: "s",
                gene_id: "GENE_SYN_A",
                bindings: &hash,
                analysis: &changed,
                projection: &projection
            }],
            Limits::default()
        )
        .is_err()
    );
}
