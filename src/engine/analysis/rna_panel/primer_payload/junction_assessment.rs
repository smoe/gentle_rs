//! J1 consumes verified saved S1 evidence, never realigns or changes its report.

use super::*;
use crate::observed_transcript_structures as s1;
use serde_json::json;

#[path = "junction_assessment/geometry.rs"]
mod geometry;
use geometry::*;
#[path = "junction_assessment/publication.rs"]
mod publication;
#[path = "junction_assessment/source.rs"]
mod source;
#[path = "junction_assessment/support.rs"]
mod support;

const TABLES: [&str; 4] = [
    "junction_assessment.tsv",
    "chain_assessment.tsv",
    "sample_support.tsv",
    "unavailable.tsv",
];

fn admit(r: &AssessObservedSpliceJunctionsRequest) -> Result<(), EngineError> {
    let b = &r.budgets;
    if r.schema != OBSERVED_JUNCTION_REQUEST_SCHEMA
        || !(1..=50).contains(&r.near_window_bp)
        || !(1..=1_000_000_000_000).contains(&b.max_input_bytes)
        || !(1..=64_000_000_000).contains(&b.max_reference_bytes)
        || !(1..=8_000_000).contains(&b.max_manifest_bytes)
        || !(1..=100_000).contains(&b.max_distinct_junctions)
        || !(1..=1_000_000).contains(&b.max_annotated_junctions_compared)
        || !(1..=100_000).contains(&b.max_near_hits_per_junction)
        || !(1..=1_000_000_000_000).contains(&b.max_total_work)
        || !(1..=128_000_000).contains(&b.max_output_bytes)
        || r.recurrence_policy.as_ref().is_some_and(|p| {
            p.min_units < 2 || p.min_units > 1_000_000 || !p.require_qualified_flanks
        })
    {
        return Err(bad("junction_assessment_admission"));
    }
    let ObservedJunctionSource::S1Report { report, result } = &r.source;
    if !is_hash(&report.sha256) || Path::new(&report.path) != Path::new(&result.report_path) {
        return Err(bad("junction_assessment_source_binding"));
    }
    for p in [&report.path, &r.output_dir]
        .into_iter()
        .chain(r.cancel_path.iter())
        .chain(r.reference_sequence.iter().map(|f| &f.file.path))
        .chain(r.sample_manifest.iter().map(|f| &f.path))
    {
        publication::raw_path(p)?;
    }
    if let Some(marker) = &r.cancel_path {
        let out = Path::new(&r.output_dir);
        let marker = Path::new(marker);
        if marker == out || marker.starts_with(out) {
            return Err(bad("junction_assessment_control_collision"));
        }
    }
    if let Some(reference) = &r.reference_sequence {
        if !is_hash(&reference.file.sha256) || reference.file.bytes > b.max_reference_bytes {
            return Err(bad("junction_assessment_reference_budget"));
        }
        if let ObservedJunctionReferenceAuthority::DeclaredExcerpt {
            contig,
            offset_0based,
            length,
        } = &reference.authority
        {
            if contig.is_empty()
                || contig.len() > 4096
                || contig.chars().any(char::is_whitespace)
                || *length == 0
                || offset_0based
                    .checked_add(*length)
                    .is_none_or(|z| z > i64::MAX as usize)
            {
                return Err(bad("junction_assessment_excerpt_binding"));
            }
        }
    }
    if r.sample_manifest
        .as_ref()
        .is_some_and(|f| !is_hash(&f.sha256) || f.bytes > b.max_manifest_bytes)
    {
        return Err(bad("junction_assessment_manifest_budget"));
    }
    Ok(())
}

struct Meter<'a> {
    b: &'a ObservedJunctionBudgets,
    input: u64,
    work: u64,
    comparisons: usize,
}
impl Meter<'_> {
    fn bytes(&mut self, n: u64) -> Result<(), EngineError> {
        self.input = self
            .input
            .checked_add(n)
            .filter(|n| *n <= self.b.max_input_bytes)
            .ok_or_else(|| bad("junction_assessment_input_budget"))?;
        self.charge(n)
    }
    fn charge(&mut self, n: u64) -> Result<(), EngineError> {
        self.work = self
            .work
            .checked_add(n)
            .filter(|n| *n <= self.b.max_total_work)
            .ok_or_else(|| bad("junction_assessment_work_budget"))?;
        Ok(())
    }
    fn comparisons(&mut self, n: usize) -> Result<(), EngineError> {
        self.comparisons = self
            .comparisons
            .checked_add(n)
            .filter(|n| *n <= self.b.max_annotated_junctions_compared)
            .ok_or_else(|| bad("junction_assessment_comparison_budget"))?;
        self.charge(n as u64)
    }
}

fn read_bound<T: serde::de::DeserializeOwned>(
    f: &PrimerPayloadFile,
    limit: u64,
    m: &mut Meter<'_>,
    ctl: &mut Control<'_>,
) -> Result<T, EngineError> {
    let mut reader = VerifiedReader::open(f, limit.min(128_000_000))?;
    m.bytes(f.bytes)?;
    let mut bytes = Vec::new();
    let mut buf = [0u8; 16384];
    loop {
        ctl.check()?;
        let n = reader.read(&mut buf).map_err(|e| bad(e.to_string()))?;
        if n == 0 {
            break;
        }
        bytes.extend_from_slice(&buf[..n]);
    }
    reader.finish()?;
    serde_json::from_slice(&bytes).map_err(|e| bad(format!("junction_assessment_json: {e}")))
}

fn empty_data() -> ObservedJunctionAssessmentData {
    ObservedJunctionAssessmentData {
        schema: OBSERVED_JUNCTION_DATA_SCHEMA.into(),
        junctions: Vec::new(),
        chains: Vec::new(),
        libraries: Vec::new(),
        unassessed_contexts: Vec::new(),
        unused_manifest_rows: Vec::new(),
        source_unavailable: Vec::new(),
    }
}

impl GentleEngine {
    /// Assess saved S1 evidence with separately bound bases and declared sample hierarchy.
    pub fn assess_observed_splice_junctions(
        &self,
        r: AssessObservedSpliceJunctionsRequest,
        progress: &mut dyn FnMut(OperationProgress) -> bool,
    ) -> Result<ObservedSpliceJunctionAssessmentReport, EngineError> {
        admit(&r)?;
        publication::execution_paths(&r)?;
        let identity = hash(&r)?;
        let mut ctl = Control::new(r.cancel_path.clone(), identity.clone(), 1, progress)?;
        let stage = staging(&r.output_dir)?;
        let mut meter = Meter {
            b: &r.budgets,
            input: 0,
            work: 0,
            comparisons: 0,
        };
        let mut data = empty_data();
        let mut source_status = "not_read".to_string();
        let mut counts = None;
        let work = (|| {
            let loaded = source::load(&r, &mut meter, &mut ctl)?;
            source_status = loaded.report.status.clone();
            data.source_unavailable = loaded.data.unavailable.clone();
            if source_status != "completed" {
                return Ok(());
            }
            let summary = loaded
                .data
                .analysis
                .as_ref()
                .ok_or_else(|| bad("junction_assessment_completed_source_without_analysis"))?;
            let mut physical = BTreeSet::new();
            for o in &summary.observations {
                ctl.check()?;
                meter.charge(o.skips.len() as u64 + 1)?;
                for &skip in &o.skips {
                    if skip.0 >= skip.1 || skip.1 > i64::MAX as usize {
                        return Err(bad("junction_assessment_coordinates"));
                    }
                    physical.insert((o.contig.clone(), skip));
                    if physical.len() > r.budgets.max_distinct_junctions {
                        return Err(bad("junction_assessment_junction_budget"));
                    }
                }
            }
            let bases = source::reference_bases(
                &r,
                summary,
                &physical,
                &loaded.catalogs,
                &mut meter,
                &mut ctl,
            )?;
            let manifest = source::manifest(&r, &mut meter, &mut ctl)?;
            meter.charge(
                (summary.records.len() as u64)
                    .checked_mul(loaded.data.unavailable.len() as u64 + 1)
                    .ok_or_else(|| bad("junction_assessment_work_budget"))?,
            )?;
            let metadata = support::metadata(
                summary,
                &loaded.contexts,
                &manifest,
                &loaded.report,
                &loaded.data.unavailable,
            );
            data.libraries = metadata.libraries;
            data.unassessed_contexts = metadata.unassessed;
            data.unused_manifest_rows = metadata.unused;
            let record_index =
                support::record_index(summary, &loaded.contexts, &mut meter, &mut ctl)?;
            for (contig, skip) in physical {
                ctl.check()?;
                let mut comparisons = Vec::new();
                let mut hits = 0usize;
                for c in &loaded.catalogs {
                    meter.charge(
                        c.models
                            .iter()
                            .map(|m| m.genomic_exons.len() as u64 + 1)
                            .sum(),
                    )?;
                    meter.comparisons(
                        c.models
                            .iter()
                            .map(|m| m.genomic_exons.len().saturating_sub(1))
                            .sum(),
                    )?;
                    let result = source::compare_source(&contig, skip, c, r.near_window_bp);
                    hits = hits
                        .checked_add(result.near_hits.len())
                        .filter(|n| *n <= r.budgets.max_near_hits_per_junction)
                        .ok_or_else(|| bad("junction_assessment_near_hit_budget"))?;
                    comparisons.push(result);
                }
                let orientation = assigned_orientation(&comparisons);
                let length = loaded
                    .catalogs
                    .iter()
                    .flat_map(|c| &c.models)
                    .find(|m| m.target_name == contig)
                    .map(|m| m.target_length);
                let right = skip
                    .1
                    .checked_add(1)
                    .filter(|n| length.is_some_and(|len| *n <= len));
                let left = (skip.0 > 0).then_some(skip.0);
                let (donor, acceptor) = match orientation {
                    Some('+') => (left, right),
                    Some('-') => (right, left),
                    _ => (None, None),
                };
                let first = bases
                    .values
                    .get(&(contig.clone(), skip.0))
                    .and_then(|b| std::str::from_utf8(b).ok());
                let last = skip
                    .1
                    .checked_sub(2)
                    .and_then(|p| bases.values.get(&(contig.clone(), p)))
                    .and_then(|b| std::str::from_utf8(b).ok());
                let motif = motif_from_dinucleotides(first, last, orientation, &bases.authority);
                let records = record_index
                    .get(&(contig.clone(), skip))
                    .cloned()
                    .unwrap_or_default();
                let support = support::counts(&records);
                let grouped_support = support::groups(&records, &manifest);
                let recurrence = r.recurrence_policy.as_ref().map(|p| {
                    support::recurrence(
                        p,
                        &records,
                        &manifest,
                        loaded.report.source_inputs_missing > 0
                            || loaded.report.source_query_units_missing > 0,
                    )
                });
                let junction_id = junction_id(&summary.reference, &contig, skip)?;
                data.junctions.push(ObservedSpliceJunctionAssessment {
                    display_id: junction_id.strip_prefix("sha256:").unwrap_or(&junction_id)[..8]
                        .into(),
                    junction_id,
                    contig,
                    skipped_interval: skip,
                    intron_first_base_1based: skip.0 + 1,
                    intron_last_base_1based: skip.1,
                    donor_exon_last_base_1based: donor,
                    acceptor_exon_first_base_1based: acceptor,
                    orientation,
                    orientation_reason: if orientation.is_some() {
                        "all_bound_models_agree"
                    } else {
                        "unassigned_model_context"
                    }
                    .into(),
                    motif,
                    sources: comparisons,
                    support,
                    records,
                    grouped_support,
                    recurrence,
                });
            }
            data.chains = summary
                .observations
                .iter()
                .map(|o| ObservedChainAssessment {
                    observation_id: o.observation_id.clone(),
                    record: o.record.clone(),
                    chain_id: o.observed_junction_chain.clone(),
                    comparisons: o.comparisons.clone(),
                    end_constraints: o.ends.clone(),
                    alignment_clipped: o.alignment_clipped,
                    rna_end_status: "not_verified".into(),
                })
                .collect();
            data.chains
                .sort_by(|a, b| a.observation_id.cmp(&b.observation_id));
            counts = Some(ObservedJunctionAssessmentCounts {
                source_records: summary.records.len(),
                distinct_junctions: data.junctions.len(),
                unassessed_contexts: data.unassessed_contexts.len(),
                annotated_comparisons: meter.comparisons,
                charged_input_bytes: meter.input,
                charged_work: meter.work,
            });
            ctl.completed = 1;
            Ok::<_, EngineError>(())
        })();
        match work {
            Err(e) if e.message == INTERRUPTED => {}
            Err(e) => return Err(e),
            Ok(()) => {}
        }
        if counts.is_none() {
            let unavailable = std::mem::take(&mut data.source_unavailable);
            data = empty_data();
            data.source_unavailable = unavailable;
        }
        ctl.stop(true);
        let mut report = ObservedSpliceJunctionAssessmentReport {
            schema: OBSERVED_JUNCTION_REPORT_SCHEMA.into(),
            content_sha256: String::new(),
            request_sha256: identity,
            request: r.clone(),
            status: if ctl.status() == "completed" && counts.is_none() {
                "unavailable"
            } else {
                ctl.status()
            }.into(),
            control_error: ctl.error.clone(),
            source_status,
            counts,
            artifacts: Vec::new(),
            nonclaims: vec![
                "Not in these bound references never means biologically novel or absent.".into(),
                "No protein, translation, expression, molecule, secretion or isoform exclusivity claim.".into(),
                "Sample hierarchy and independence are producer-declared, not authenticated.".into(),
                "Motif authority is verified only for a full assembly hashed to EOF; excerpt-to-assembly identity remains declared.".into(),
                "Possible-only and definite support are disjoint; grouping axes and contexts are not additive.".into(),
                "Payload-end constraints and a complete junction chain do not establish complete RNA molecules.".into(),
            ],
        };
        report.artifacts = publication::artifacts(stage.path(), &data, &report)?;
        report.content_sha256 = hash(&report)?;
        validate_report(&report)?;
        publication::publish(stage, &report, &mut |_| Ok(()))?;
        Ok(report)
    }

    pub fn observed_splice_junction_result_binding(
        request: &AssessObservedSpliceJunctionsRequest,
        report: &ObservedSpliceJunctionAssessmentReport,
    ) -> Result<ObservedSpliceJunctionResultBinding, EngineError> {
        if hash(request)? != report.request_sha256 {
            return Err(bad("junction_assessment_executed_request_mismatch"));
        }
        validate_report(report)
    }
}

fn validate_report(
    r: &ObservedSpliceJunctionAssessmentReport,
) -> Result<ObservedSpliceJunctionResultBinding, EngineError> {
    admit(&r.request)?;
    let mut basis = r.clone();
    basis.content_sha256.clear();
    let expected = TABLES
        .into_iter()
        .chain(["assessment.json"])
        .collect::<BTreeSet<_>>();
    let names = r
        .artifacts
        .iter()
        .map(|a| a.path.as_str())
        .collect::<BTreeSet<_>>();
    if r.schema != OBSERVED_JUNCTION_REPORT_SCHEMA
        || hash(&basis)? != r.content_sha256
        || hash(&r.request)? != r.request_sha256
        || !(binding::outcome(&r.status, &r.control_error)
            || (r.status == "unavailable" && r.control_error.is_none()))
        || (r.status == "completed" && r.counts.is_none())
        || (r.status == "unavailable" && r.counts.is_some())
        || names != expected
        || r.artifacts.len() != 5
        || r.artifacts
            .iter()
            .any(|a| !is_hash(&a.sha256) || a.bytes > r.request.budgets.max_output_bytes)
    {
        return Err(bad("junction_assessment_report_binding"));
    }
    Ok(ObservedSpliceJunctionResultBinding {
        schema: OBSERVED_JUNCTION_BINDING_SCHEMA.into(),
        request_sha256: r.request_sha256.clone(),
        report_content_sha256: r.content_sha256.clone(),
        report_path: Path::new(&r.request.output_dir)
            .join("report.json")
            .display()
            .to_string(),
        outcome: r.status.clone(),
    })
}

#[cfg(test)]
#[path = "junction_assessment/tests.rs"]
mod tests;
