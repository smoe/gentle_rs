//! Inline synthetic SAM projection regressions; no real reads, tools or provider claims.

use super::super::sam_intake::tests::{bound, complete, fixture, manifest_file};
use super::*;
use std::fs;
use std::path::Path;

fn rewrite(r: &mut sam_intake::Request, m: &mut ObservedSamQueryUnitManifest, raw: &str) {
    fs::write(&r.sam.path, raw).unwrap();
    r.sam = bound(Path::new(&r.sam.path));
    m.sam_sha256 = r.sam.sha256.clone();
    manifest_file(r, m);
}

fn run(r: &sam_intake::Request) -> ObservedSamAnalysis {
    summarize_bound(r, &BTreeMap::new(), 1_000_000, &mut || Ok(()))
        .unwrap()
        .0
}

#[test]
fn sam_report_roundtrip_retains_query_units_and_unknown_axes_without_read_identity() {
    let (_dir, r, _m) = fixture();
    let result = run(&r);
    assert_eq!(result.units.len(), 4);
    assert_eq!(result.observations.len(), 3);
    assert!(result.observations.iter().all(|o| o.on_optimal_path));
    assert_eq!(result.chain_families.len(), 1);
    assert_eq!(result.chain_families[0].units.len(), 1);
    assert_eq!(
        result.chain_families[0]
            .counts
            .possible_only
            .unassessed_only,
        1
    );
    assert!(result.units.iter().all(|u| u.definite_chains.is_empty()));
    let json = serde_json::to_value(&result).unwrap();
    let back: ObservedSamAnalysis = serde_json::from_value(json.clone()).unwrap();
    assert_eq!(serde_json::to_value(&back).unwrap(), json);
    assert!(json.get("records").is_none());
    assert!(json["observations"][0].get("record_index").is_none());
    let absent = result
        .units
        .iter()
        .find(|u| u.unit.unit_id == "unit_missing")
        .unwrap();
    assert_eq!(absent.unclipped_axis_length, None);
    assert!(
        absent
            .unavailable
            .contains(&ObservedSamUnavailableReason::NotProvided)
    );
    assert!(!absent.all_query_views_assessed);
}

#[test]
fn sam_report_nonoptimal_spliced_observations_do_not_contribute_to_catalogs() {
    let (_dir, mut r, mut m) = fixture();
    let raw = fs::read_to_string(&r.sam.path)
        .unwrap()
        .replace("10M10N20M", "10M10N10M10S");
    let other = m.units[1].rows.pop().unwrap();
    m.units[0].rows.push(other);
    m.units[0].search = Some(complete());
    rewrite(&mut r, &mut m, &raw);
    let data = run(&r);
    assert_eq!(data.observations.len(), 3);
    assert_eq!(
        data.observations
            .iter()
            .filter(|o| o.on_optimal_path)
            .count(),
        1
    );
    assert!(
        data.observations
            .iter()
            .filter(|o| !o.skips.is_empty())
            .all(|o| !o.on_optimal_path)
    );
    assert!(data.chain_families.is_empty());
    assert!(data.junctions.is_empty());
    let unit = data
        .units
        .iter()
        .find(|u| u.unit.unit_id == "unit_a")
        .unwrap();
    assert_eq!(unit.topologies, [Topology::Unspliced]);
    assert_eq!(unit.hypotheses.len(), 1);
}

#[test]
fn sam_report_qualification_is_exclusive_and_deduplicated_across_optimal_rows() {
    let (_dir, mut r, mut m) = fixture();
    let raw = fs::read_to_string(&r.sam.path)
        .unwrap()
        .replace("same_name\t2048", "same_name\t2064");
    m.units[0].search = Some(complete());
    m.units[0].masks = Some(ObservedSamMaskWitness {
        coordinate_convention: "zero_based_half_open_unclipped_query".into(),
        unclipped_axis_length: 30,
        intervals: vec![(8, 10)],
        source_sha256: format!("sha256:{}", "a".repeat(64)),
        derivation: "synthetic primer interval".into(),
    });
    rewrite(&mut r, &mut m, &raw);
    let data = run(&r);
    let family = &data.chain_families[0];
    assert_eq!(family.units.len(), 1);
    assert_eq!(
        family.units[0].certainty,
        ObservedSamMappingCertainty::Definite
    );
    assert_eq!(
        family.units[0].qualification,
        ObservedSamQualificationGroup::Mixed
    );
    assert_eq!(family.counts.definite.mixed, 1);
    assert_eq!(family.counts.definite.qualified_only, 0);
    assert_eq!(data.junctions[0].counts.definite.mixed, 1);
    assert_eq!(data.junctions[0].units.len(), 1);
    assert!(
        data.observations
            .iter()
            .any(|o| o.flank_qualification == [Qualification::Qualified])
    );
    assert!(
        data.observations
            .iter()
            .any(|o| o.flank_qualification == [Qualification::InsufficientSupport])
    );
}

#[test]
fn sam_report_shared_junction_can_be_definite_when_full_chains_are_alternatives() {
    let (_dir, mut r, mut m) = fixture();
    let raw = fs::read_to_string(&r.sam.path)
        .unwrap()
        .replacen("10M10N20M", "10M10N10M10N10M", 1);
    m.units[0].search = Some(complete());
    rewrite(&mut r, &mut m, &raw);
    let data = run(&r);
    let unit = data
        .units
        .iter()
        .find(|u| u.unit.unit_id == "unit_a")
        .unwrap();
    assert_eq!(unit.possible_chains.len(), 2);
    assert!(unit.definite_chains.is_empty());
    let shared = data
        .junctions
        .iter()
        .find(|j| j.skipped_interval == (110, 120))
        .unwrap();
    let other = data
        .junctions
        .iter()
        .find(|j| j.skipped_interval == (130, 140))
        .unwrap();
    assert_eq!(shared.counts.definite.unassessed_only, 1);
    assert_eq!(other.counts.possible_only.unassessed_only, 1);
    assert_eq!(shared.units.len(), 1);
    assert_eq!(other.units.len(), 1);
}

#[test]
fn sam_report_refused_sources_are_not_relabelled_as_missing_genes_or_empty_models() {
    let (_dir, mut r, _m) = fixture();
    r.annotation_sources.clear();
    let states = BTreeMap::from([(
        "source_a".into(),
        ObservedSamAnnotationAvailability::Refused,
    )]);
    let (data, _) = summarize_bound(&r, &states, 1_000_000, &mut || Ok(())).unwrap();
    assert!(data.annotation_bindings.is_empty());
    assert_eq!(
        data.gene_resolutions[0].availability,
        ObservedSamAnnotationAvailability::Refused
    );
    assert!(data.gene_resolutions[0].comparison.is_none());
    assert!(
        data.observations
            .iter()
            .all(|o| o.comparisons.iter().all(|c| c.comparison.is_none()))
    );
}

#[test]
fn sam_report_retains_length_but_not_sam_sequence_or_quality_text() {
    let (_dir, mut r, mut m) = fixture();
    let sequence = format!("{}AC", "ACGT".repeat(7));
    let quality = "~".repeat(30);
    let raw = fs::read_to_string(&r.sam.path).unwrap().replacen(
        "\t*\t*\tRG",
        &format!("\t{sequence}\t{quality}\tRG"),
        1,
    );
    rewrite(&mut r, &mut m, &raw);
    let data = run(&r);
    let serialized = serde_json::to_string(&data).unwrap();
    assert!(!serialized.contains(&sequence));
    assert!(!serialized.contains(&quality));
    let row = data
        .units
        .iter()
        .flat_map(|u| &u.rows)
        .find(|r| r.ordinal == 1)
        .unwrap();
    assert_eq!(row.sequence_length, Some(30));
    assert_eq!(row.quality_length, Some(30));
}

#[test]
fn sam_report_projection_budget_and_cancel_refuse_without_publication() {
    let (dir, mut r, _m) = fixture();
    assert!(summarize_bound(&r, &BTreeMap::new(), 1, &mut || Ok(())).is_err());
    r.kernel_budgets.max_alternatives_per_record = 1_000_000;
    let error = summarize_bound(&r, &BTreeMap::new(), 1_000_000, &mut || Ok(())).unwrap_err();
    assert_eq!(error.message, "observed_sam_projection_memory_budget");
    r.kernel_budgets.max_alternatives_per_record = 100;
    let mut polls = 0;
    summarize_bound(&r, &BTreeMap::new(), 1_000_000, &mut || {
        polls += 1;
        Ok(())
    })
    .unwrap();
    for stop in 1..=polls {
        let mut count = 0;
        let error = summarize_bound(&r, &BTreeMap::new(), 1_000_000, &mut || {
            count += 1;
            if count == stop {
                Err(bad("synthetic_cancel"))
            } else {
                Ok(())
            }
        })
        .unwrap_err();
        assert_eq!(error.message, "synthetic_cancel");
    }
    assert_eq!(fs::read_dir(dir.path()).unwrap().count(), 3);
}
