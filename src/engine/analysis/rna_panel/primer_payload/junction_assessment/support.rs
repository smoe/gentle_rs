//! Original-record reconciliation and declared hierarchy counts; no molecule/independence inference.

use super::source::Manifest;
use super::*;

pub(super) type Contexts = BTreeMap<s1::ReadRecord, Vec<ObservedSampleContext>>;
pub(super) fn record_index(
    summary: &s1::Summary,
    contexts: &Contexts,
    m: &mut Meter<'_>,
    ctl: &mut Control<'_>,
) -> Result<BTreeMap<(String, s1::Interval), Vec<ObservedJunctionRecordSupport>>, EngineError> {
    let observations = summary
        .observations
        .iter()
        .map(|o| (o.observation_id.as_str(), o))
        .collect::<BTreeMap<_, _>>();
    let families = summary
        .families
        .iter()
        .map(|f| (f.chain_id.as_str(), f))
        .collect::<BTreeMap<_, _>>();
    if observations.len() != summary.observations.len() || families.len() != summary.families.len()
    {
        return Err(bad("junction_assessment_source_duplicate_identity"));
    }
    let mut out = BTreeMap::<(String, s1::Interval), Vec<ObservedJunctionRecordSupport>>::new();
    for r in &summary.records {
        ctl.check()?;
        m.charge(1)?;
        if !r
            .definite_chains
            .iter()
            .all(|id| r.possible_chains.contains(id))
        {
            return Err(bad("junction_assessment_source_certainty_binding"));
        }
        let mut keys = BTreeMap::<(String, s1::Interval), bool>::new();
        for id in &r.possible_chains {
            let f = families
                .get(id.as_str())
                .ok_or_else(|| bad("junction_assessment_source_chain_binding"))?;
            m.charge(f.skips.len() as u64 + 1)?;
            let definite = r.search_declared_complete
                && r.all_query_views_assessed
                && r.definite_chains.contains(id);
            for skip in &f.skips {
                *keys.entry((f.contig.clone(), *skip)).or_default() |= definite;
            }
        }
        let mut qualification = BTreeMap::<(String, s1::Interval), BTreeSet<u8>>::new();
        for id in r.hypotheses.iter().flatten().collect::<BTreeSet<_>>() {
            m.charge(1)?;
            let o = observations
                .get(id.as_str())
                .ok_or_else(|| bad("junction_assessment_source_path_binding"))?;
            if o.record != r.record || o.skips.len() != o.flank_qualification.len() {
                return Err(bad("junction_assessment_source_flank_binding"));
            }
            for (skip, q) in o.skips.iter().zip(&o.flank_qualification) {
                m.charge(1)?;
                qualification
                    .entry((o.contig.clone(), *skip))
                    .or_default()
                    .insert(match q {
                        s1::Qualification::Qualified => 0,
                        s1::Qualification::InsufficientSupport => 1,
                        s1::Qualification::Unassessed => 2,
                    });
            }
        }
        for (key, definite) in keys {
            let q = qualification
                .get(&key)
                .ok_or_else(|| bad("junction_assessment_source_skip_without_witness"))?;
            let qualification = if q.len() > 1 {
                "mixed"
            } else {
                match q.first() {
                    Some(0) => "qualified_only",
                    Some(1) => "insufficient_only",
                    _ => "unassessed_only",
                }
            };
            out.entry(key)
                .or_default()
                .push(ObservedJunctionRecordSupport {
                    record: r.record.clone(),
                    certainty: if definite { "definite" } else { "possible" }.into(),
                    qualification: qualification.into(),
                    contexts: contexts.get(&r.record).cloned().unwrap_or_default(),
                });
        }
    }
    for rows in out.values_mut() {
        rows.sort_by(|a, b| a.record.cmp(&b.record));
    }
    Ok(out)
}

fn add(counts: &mut ObservedJunctionSupportCounts, row: &ObservedJunctionRecordSupport) {
    let q = if row.certainty == "definite" {
        &mut counts.definite
    } else {
        &mut counts.possible_only
    };
    match row.qualification.as_str() {
        "qualified_only" => q.qualified_only += 1,
        "insufficient_only" => q.insufficient_only += 1,
        "unassessed_only" => q.unassessed_only += 1,
        _ => q.mixed += 1,
    }
}
pub(super) fn counts(records: &[ObservedJunctionRecordSupport]) -> ObservedJunctionSupportCounts {
    let mut c = ObservedJunctionSupportCounts::default();
    for r in records {
        add(&mut c, r);
    }
    c
}
pub(super) fn groups(
    records: &[ObservedJunctionRecordSupport],
    manifest: &Manifest,
) -> Vec<ObservedGroupedJunctionSupport> {
    let mut groups = BTreeMap::<(String, String), ObservedJunctionSupportCounts>::new();
    for r in records {
        let mut keys = BTreeSet::new();
        for context in &r.contexts {
            if let Some(row) = manifest.get(context) {
                keys.extend([
                    ("library".into(), row.library_id.clone()),
                    ("cdna_prep".into(), row.cdna_prep_id.clone()),
                    ("biological_sample".into(), row.biological_sample_id.clone()),
                ]);
                if let Some(line) = &row.cell_line {
                    keys.insert(("cell_line".into(), line.clone()));
                }
            }
        }
        for key in keys {
            add(groups.entry(key).or_default(), r);
        }
    }
    groups
        .into_iter()
        .map(|((dimension, id), counts)| ObservedGroupedJunctionSupport {
            dimension,
            id,
            counts,
        })
        .collect()
}

pub(super) fn recurrence(
    policy: &ObservedRecurrencePolicy,
    records: &[ObservedJunctionRecordSupport],
    manifest: &Manifest,
    source_gaps: bool,
) -> ObservedJunctionRecurrence {
    let mut units = BTreeSet::new();
    let mut unknown = source_gaps;
    for record in records {
        if record.qualification != "qualified_only" {
            unknown |= record.qualification != "insufficient_only";
            continue;
        }
        if policy.certainty == ObservedRecurrenceCertainty::Definite
            && record.certainty != "definite"
        {
            unknown = true;
            continue;
        }
        let mut known = BTreeSet::new();
        let mut missing = record.contexts.is_empty();
        for context in &record.contexts {
            match manifest.get(context) {
                Some(row) => {
                    known.insert(match policy.unit {
                        ObservedRecurrenceUnit::Library => row.library_id.clone(),
                        ObservedRecurrenceUnit::CdnaPrep => row.cdna_prep_id.clone(),
                        ObservedRecurrenceUnit::BiologicalSample => {
                            row.biological_sample_id.clone()
                        }
                    });
                }
                None => missing = true,
            }
        }
        // One physical record in multiple claimed units cannot establish independent recurrence.
        if !missing && known.len() == 1 {
            units.extend(known);
        } else {
            unknown = true;
        }
    }
    let (status, reason) = if units.len() >= policy.min_units {
        ("recurrent_candidate", "explicit_policy_qualified_units")
    } else if unknown {
        ("not_assessable", "incomplete_context_or_required_evidence")
    } else {
        (
            "not_recurrent_under_policy",
            "qualified_units_below_explicit_threshold",
        )
    };
    ObservedJunctionRecurrence {
        status: status.into(),
        reason: reason.into(),
        independent_units: units.into_iter().collect(),
    }
}

pub(super) struct Metadata {
    pub libraries: Vec<ObservedLibraryDenominator>,
    pub unassessed: Vec<ObservedJunctionUnassessedContext>,
    pub unused: Vec<ObservedSampleBinding>,
}
pub(super) fn metadata(
    summary: &s1::Summary,
    contexts: &Contexts,
    manifest: &Manifest,
    source: &ObservedTranscriptStructuresReport,
    unavailable: &[ObservedTranscriptUnavailable],
) -> Metadata {
    #[derive(Default)]
    struct Sets {
        all: BTreeSet<s1::ReadRecord>,
        gene: BTreeSet<s1::ReadRecord>,
        excluded: BTreeSet<s1::ReadRecord>,
        unassessed: BTreeSet<s1::ReadRecord>,
    }
    let mut used = BTreeSet::new();
    let mut libraries = BTreeMap::<String, Sets>::new();
    let mut unassessed = Vec::new();
    let coherent = summary
        .observations
        .iter()
        .filter(|o| {
            o.comparisons.iter().any(|c| {
                c.source_id == summary.reference_gene.source_id
                    && c.resolution == s1::GeneResolution::Resolved
                    && c.models
                        .iter()
                        .any(|m| m.relation != s1::ModelRelation::OutsideModel)
            })
        })
        .map(|o| o.record.clone())
        .collect::<BTreeSet<_>>();
    let gene_assessed = summary.gene_resolutions.iter().any(|c| {
        c.source_id == summary.reference_gene.source_id
            && c.resolution == s1::GeneResolution::Resolved
    });
    for record in &summary.records {
        let rows = contexts
            .get(&record.record)
            .map(Vec::as_slice)
            .unwrap_or(&[]);
        if rows.is_empty() {
            unassessed.push(ObservedJunctionUnassessedContext {
                record: record.record.clone(),
                context: None,
                reason: "context_not_retained_by_s1".into(),
            });
        }
        for c in rows {
            if let Some(row) = manifest.get(c) {
                used.insert(c.clone());
                let sets = libraries.entry(row.library_id.clone()).or_default();
                sets.all.insert(record.record.clone());
                if coherent.contains(&record.record) {
                    sets.gene.insert(record.record.clone());
                }
                if unavailable.iter().any(|u| {
                    u.query_id.is_none()
                        && u.file_id == c.file_id
                        && u.record_index == record.record.record_index
                        && u.read_id == record.record.read_id
                }) {
                    sets.excluded.insert(record.record.clone());
                }
                if !record.all_query_views_assessed
                    || !record.search_declared_complete
                    || !gene_assessed
                {
                    sets.unassessed.insert(record.record.clone());
                }
            } else {
                unassessed.push(ObservedJunctionUnassessedContext {
                    record: record.record.clone(),
                    context: Some(c.clone()),
                    reason: "unassessed_context".into(),
                });
            }
        }
    }
    unassessed.sort_by(|a, b| (&a.record, &a.context).cmp(&(&b.record, &b.context)));
    let complete = unassessed.is_empty()
        && source.source_inputs_missing == 0
        && source.source_query_units_missing == 0;
    Metadata {
        libraries: libraries
            .into_iter()
            .map(|(library_id, s)| ObservedLibraryDenominator {
                library_id,
                total_records: s.all.len(),
                gene_coherent_records: s.gene.len(),
                excluded_records: s.excluded.len(),
                unassessed_records: s.unassessed.len(),
                context_coverage_complete: complete,
            })
            .collect(),
        unassessed,
        unused: manifest
            .iter()
            .filter(|(c, _)| !used.contains(*c))
            .map(|(_, v)| v.clone())
            .collect(),
    }
}
