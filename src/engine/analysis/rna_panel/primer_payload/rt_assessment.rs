//! Bounded prospective RT screening, reusing oligo validation and verified FASTA byte intake.

use super::*;

const SEEDS: [usize; 4] = [6, 8, 10, 12];
const LINE_LIMIT: usize = 1_048_576;

fn counts(n: usize) -> RtPrimerSeedCounts {
    RtPrimerSeedCounts {
        assessed_windows: 0,
        ambiguous_windows: 0,
        mismatch_histogram: vec![0; n + 1],
        terminal_run_histogram: vec![0; n + 1],
        exact_matching_records: 0,
    }
}

fn admission(request: &AssessRtPrimerCandidatesRequest) -> Result<Vec<Vec<u8>>, EngineError> {
    let b = &request.budgets;
    if request.schema != RT_PRIMER_ASSESSMENT_REQUEST_SCHEMA
        || !(1..=16).contains(&request.primers.len())
        || !(1..=16).contains(&request.references.len())
        || request.intended_sites.len() > 1024
        || !(1..=256_000_000).contains(&b.max_input_bytes)
        || !(1..=16_000_000).contains(&b.max_reference_bases)
        || !(1..=16_000_000).contains(&b.max_record_bases)
        || !(1..=100_000).contains(&b.max_records)
        || !(1..=1_000_000_000).contains(&b.max_base_comparisons)
        || b.max_examples_per_cell > 32
    {
        return Err(bad("rt_assessment_schema_or_budget"));
    }
    let mut ids = BTreeSet::new();
    let mut tails = Vec::new();
    for primer in &request.primers {
        if !ids.insert(&primer.primer_id)
            || primer.core_orientation != PrimerCoreOrientation::Antisense3prime
            || primer.components.len() > 2
        {
            return Err(bad("rt_assessment_requires_unique_rt_primers"));
        }
        // Validate independently: different ordered handles may legitimately share a tail.
        let panel = RnaPrimerPanel {
            schema: RNA_PRIMER_PANEL_SCHEMA.into(),
            content_sha256: String::new(),
            primers: vec![primer.clone()],
            policy: RnaPrimerTracePolicy::default(),
        };
        GentleEngine::validate_rna_primer_panel(
            &RnaPrimerPanelBinding::Inline {
                panel: Box::new(panel),
            },
            true,
        )?;
        let core = primer
            .components
            .iter()
            .find(|c| c.role == "target_core")
            .ok_or_else(|| bad("rt_assessment_core_missing"))?;
        tails.push(primer.sequence.as_bytes()[core.start..core.end].to_vec());
    }
    let mut refs = BTreeSet::new();
    let mut bytes = 0u64;
    for r in &request.references {
        let circular = r.circular_record_ids.iter().collect::<BTreeSet<_>>();
        if !safe_id(&r.reference_id)
            || !refs.insert(&r.reference_id)
            || r.source.trim().is_empty()
            || r.source.len() > 4096
            || r.fasta.is_some() == r.unavailable_reason.is_some()
            || r.unavailable_reason
                .as_ref()
                .is_some_and(|s| s.trim().is_empty() || s.len() > 4096)
            || circular.len() != r.circular_record_ids.len()
            || circular.len() > b.max_records
            || circular.iter().any(|s| {
                s.is_empty()
                    || s.len() > 512
                    || s.chars().any(|c| c.is_whitespace() || c.is_control())
            })
        {
            return Err(bad("rt_assessment_reference_binding"));
        }
        if let Some(f) = &r.fasta {
            if !is_hash(&f.sha256)
                || f.path.trim().is_empty()
                || f.path.len() > 4096
                || f.bytes == 0
            {
                return Err(bad("rt_assessment_fasta_binding"));
            }
            bytes = bytes
                .checked_add(f.bytes)
                .filter(|v| *v <= b.max_input_bytes)
                .ok_or_else(|| bad("rt_assessment_input_byte_budget"))?;
        } else if !circular.is_empty() {
            return Err(bad("rt_assessment_missing_reference_topology"));
        }
    }
    let mut sites = BTreeSet::new();
    for s in &request.intended_sites {
        let primer = request
            .primers
            .iter()
            .position(|p| p.primer_id == s.primer_id)
            .ok_or_else(|| bad("rt_assessment_intended_primer"))?;
        let r = request
            .references
            .iter()
            .find(|r| r.reference_id == s.reference_id)
            .ok_or_else(|| bad("rt_assessment_intended_reference"))?;
        if r.fasta.is_none()
            || s.record_id.is_empty()
            || s.record_id.len() > 512
            || s.record_id
                .chars()
                .any(|c| c.is_whitespace() || c.is_control())
            || s.tail_end_0based_exclusive.checked_sub(s.tail_start_0based)
                != Some(tails[primer].len())
            || !matches!(s.template_strand, '+' | '-')
            || (r.scope != RtPrimerBackgroundScope::Genome && s.template_strand != '+')
            || !sites.insert((
                &s.primer_id,
                &s.reference_id,
                &s.record_id,
                s.tail_start_0based,
                s.template_strand,
            ))
        {
            return Err(bad("rt_assessment_intended_site_geometry"));
        }
    }
    Ok(tails)
}

fn record(
    request: &AssessRtPrimerCandidatesRequest,
    reference: &RtPrimerBackground,
    name: &str,
    sequence: &[u8],
    tails: &[Vec<u8>],
    cells: &mut [RtPrimerSeedAssessment],
    work: &mut u64,
    validated_sites: &mut BTreeSet<usize>,
    ctl: &mut Control<'_>,
) -> Result<(), EngineError> {
    if sequence.is_empty() {
        return Err(bad("rt_assessment_empty_fasta_record"));
    }
    let circular = reference.circular_record_ids.iter().any(|id| id == name);
    let strands: &[char] = if reference.scope == RtPrimerBackgroundScope::Genome {
        &['+', '-']
    } else {
        &['+']
    };
    let sites = request
        .intended_sites
        .iter()
        .enumerate()
        .filter(|(_, s)| s.reference_id == reference.reference_id && s.record_id == name)
        .collect::<Vec<_>>();
    for &(index, s) in &sites {
        let p = request
            .primers
            .iter()
            .position(|p| p.primer_id == s.primer_id)
            .expect("admitted");
        let context = sequence
            .get(s.tail_start_0based..s.tail_end_0based_exclusive)
            .ok_or_else(|| bad("rt_assessment_intended_site_out_of_bounds"))?;
        let expected = if s.template_strand == '+' {
            GentleEngine::reverse_complement_bytes(&tails[p])
        } else {
            tails[p].clone()
        };
        if context != expected {
            return Err(bad("rt_assessment_intended_site_not_exact"));
        }
        validated_sites.insert(index);
    }
    for (pi, tail) in tails.iter().enumerate() {
        for (ni, &n) in SEEDS.iter().enumerate() {
            let cell = &mut cells[pi * SEEDS.len() + ni];
            if tail.len() < n {
                continue;
            }
            if sequence.len() < n {
                cell.records_shorter_than_seed += 1;
                continue;
            }
            let windows = if circular {
                sequence.len()
            } else {
                sequence.len() - n + 1
            };
            // Bound work before scanning a record/candidate/seed, including ambiguous windows.
            let cost = (windows as u64)
                .checked_mul(n as u64)
                .and_then(|v| v.checked_mul(strands.len() as u64))
                .ok_or_else(|| bad("rt_assessment_comparison_budget"))?;
            *work = work
                .checked_add(cost)
                .filter(|v| *v <= request.budgets.max_base_comparisons)
                .ok_or_else(|| bad("rt_assessment_comparison_budget"))?;
            let mut exact_record = false;
            let mut background_exact_record = false;
            let excluded = sites
                .iter()
                .filter(|(_, s)| s.primer_id == cell.primer_id)
                .map(|(_, s)| {
                    (
                        if s.template_strand == '+' {
                            s.tail_start_0based
                        } else {
                            s.tail_end_0based_exclusive - n
                        },
                        s.template_strand,
                    )
                })
                .collect::<BTreeSet<_>>();
            for start in 0..windows {
                if start % 1024 == 0 {
                    ctl.check()?;
                }
                let mut window = [0u8; 12];
                for j in 0..n {
                    window[j] = sequence[(start + j) % sequence.len()];
                }
                let remaining = strands
                    .iter()
                    .copied()
                    .filter(|strand| !excluded.contains(&(start, *strand)))
                    .collect::<Vec<_>>();
                if !window[..n].iter().all(|x| b"ACGT".contains(x)) {
                    cell.all_sites.ambiguous_windows += 1;
                    if !remaining.is_empty() {
                        cell.excluding_intended_orientations.ambiguous_windows += 1;
                    }
                    continue;
                }
                let mut values = Vec::with_capacity(2);
                for &strand in strands {
                    // '+' RNA template is reverse-complemented into primer 5′–3′ coordinates.
                    let oriented = if strand == '+' {
                        GentleEngine::reverse_complement_bytes(&window[..n])
                    } else {
                        window[..n].to_vec()
                    };
                    let seed = &tail[tail.len() - n..];
                    let mismatches = seed.iter().zip(&oriented).filter(|(a, b)| a != b).count();
                    let run = seed
                        .iter()
                        .rev()
                        .zip(oriented.iter().rev())
                        .take_while(|(a, b)| a == b)
                        .count();
                    values.push((strand, mismatches, run));
                }
                let update = |c: &mut RtPrimerSeedCounts, vals: &[(char, usize, usize)]| {
                    if vals.is_empty() {
                        return false;
                    }
                    c.assessed_windows += 1;
                    let mismatches = vals.iter().map(|x| x.1).min().expect("nonempty");
                    c.mismatch_histogram[mismatches] += 1;
                    c.terminal_run_histogram[vals.iter().map(|x| x.2).max().expect("nonempty")] +=
                        1;
                    mismatches == 0
                };
                let exact = update(&mut cell.all_sites, &values);
                let background = values
                    .iter()
                    .copied()
                    .filter(|x| remaining.contains(&x.0))
                    .collect::<Vec<_>>();
                background_exact_record |=
                    update(&mut cell.excluding_intended_orientations, &background);
                exact_record |= exact;
                if exact {
                    if cell.exact_examples.len() < request.budgets.max_examples_per_cell {
                        cell.exact_examples.push(RtPrimerSeedExample {
                            record_id: name.into(),
                            start_0based: start,
                            end_0based_exclusive: start + n,
                            topology: if circular {
                                InlineSequenceTopology::Circular
                            } else {
                                InlineSequenceTopology::Linear
                            },
                            template_strands: values
                                .iter()
                                .filter(|x| x.1 == 0)
                                .map(|x| x.0)
                                .collect(),
                        });
                    } else {
                        cell.exact_examples_truncated = true;
                    }
                }
            }
            cell.all_sites.exact_matching_records += usize::from(exact_record);
            cell.excluding_intended_orientations.exact_matching_records +=
                usize::from(background_exact_record);
        }
    }
    Ok(())
}

impl GentleEngine {
    /// Validate the existing request contract without opening backgrounds, for empty redesigns.
    pub(in crate::engine) fn validate_rt_primer_assessment_request(
        request: &AssessRtPrimerCandidatesRequest,
    ) -> Result<(), EngineError> {
        admission(request).map(|_| ())
    }

    /// Read-only prospective assessment. Errors/cancellation return no partial report or project effects.
    pub fn assess_rt_primer_candidates(
        &self,
        request: AssessRtPrimerCandidatesRequest,
        progress: &mut dyn FnMut(OperationProgress) -> bool,
    ) -> Result<RtPrimerAssessmentReport, EngineError> {
        let tails = admission(&request)?;
        let identity = hash(&request)?;
        let mut ctl = Control::new(None, identity.clone(), request.references.len(), progress)?;
        if ctl.stop(true) {
            return Err(bad("rt_assessment_cancelled"));
        }
        let mut report = RtPrimerAssessmentReport {
            schema: RT_PRIMER_ASSESSMENT_SCHEMA.into(), content_sha256: String::new(),
            request_sha256: identity, request: request.clone(), references: Vec::new(), cells: Vec::new(),
            base_comparisons: 0, abundance_status: "not_assessed".into(), retrospective_status: "not_assessed".into(),
            warnings: vec![
                "Terminal compatibility counts are not RT initiation probabilities, experimental validation, full-tail specificity or order approval.".into(),
                "Reference scopes may overlap; do not sum cells as independent sites or RNA molecules. Missing scopes and ambiguous windows remain unassessed.".into(),
                "Genome searches both template orientations; RNA inputs use biological 5-to-3 DNA-alphabet sequence and search only its complementary primer orientation.".into(),
                "Histograms count each physical window once: minimum positional mismatches and independently longest terminal run across eligible orientations. Bins are separate marginal summaries.".into(),
                "Intended exclusions remove only the verified intended template orientation. Exact example lists describe all sites and may be capped independently of exhaustive counts.".into(),
                "This bounded scan is not a whole-genome scalability claim. No abundance weights, retrospective old-primer cohort, isoform/CDS classification or thermodynamic assessment is performed.".into(),
            ],
        };
        let mut total_bases = 0usize;
        let mut total_records = 0usize;
        let mut validated_sites = BTreeSet::new();
        for reference in &request.references {
            let mut cells = request
                .primers
                .iter()
                .flat_map(|p| {
                    SEEDS.map(|n| RtPrimerSeedAssessment {
                        primer_id: p.primer_id.clone(),
                        reference_id: reference.reference_id.clone(),
                        seed_length: n,
                        status: "not_assessed".into(),
                        reason: reference.unavailable_reason.clone(),
                        records_shorter_than_seed: 0,
                        all_sites: counts(n),
                        excluding_intended_orientations: counts(n),
                        exact_examples: Vec::new(),
                        exact_examples_truncated: false,
                    })
                })
                .collect::<Vec<_>>();
            let mut assessment = RtPrimerBackgroundAssessment {
                reference_id: reference.reference_id.clone(),
                status: "not_assessed".into(),
                reason: reference.unavailable_reason.clone(),
                records: 0,
                bases: 0,
                template_strands: if reference.scope == RtPrimerBackgroundScope::Genome {
                    vec!['+', '-']
                } else {
                    vec!['+']
                },
            };
            if let Some(file) = &reference.fasta {
                let mut reader =
                    BufReader::new(VerifiedReader::open(file, request.budgets.max_input_bytes)?);
                let mut seen = BTreeSet::new();
                let mut name = String::new();
                let mut sequence = Vec::new();
                while let Some(text) = line(&mut reader, LINE_LIMIT)? {
                    ctl.check()?;
                    if let Some(header) = text.strip_prefix('>') {
                        if !name.is_empty() {
                            record(
                                &request,
                                reference,
                                &name,
                                &sequence,
                                &tails,
                                &mut cells,
                                &mut report.base_comparisons,
                                &mut validated_sites,
                                &mut ctl,
                            )?;
                        }
                        name = header.split_ascii_whitespace().next().unwrap_or("").into();
                        total_records += 1;
                        if name.is_empty()
                            || name.chars().any(|c| c.is_whitespace() || c.is_control())
                            || name.len() > 512
                            || !seen.insert(name.clone())
                            || total_records > request.budgets.max_records
                        {
                            return Err(bad("rt_assessment_fasta_record_identity_or_budget"));
                        }
                        assessment.records += 1;
                        sequence.clear();
                    } else {
                        if name.is_empty()
                            || text.is_empty()
                            || !text
                                .bytes()
                                .all(|x| b"ACGTRYSWKMBDHVNacgtryswkmbdhvn".contains(&x))
                        {
                            return Err(bad("rt_assessment_fasta_sequence"));
                        }
                        total_bases = total_bases
                            .checked_add(text.len())
                            .filter(|v| *v <= request.budgets.max_reference_bases)
                            .ok_or_else(|| bad("rt_assessment_reference_base_budget"))?;
                        if sequence
                            .len()
                            .checked_add(text.len())
                            .is_none_or(|v| v > request.budgets.max_record_bases)
                        {
                            return Err(bad("rt_assessment_record_base_budget"));
                        }
                        assessment.bases += text.len();
                        sequence.extend(text.bytes().map(|x| x.to_ascii_uppercase()));
                    }
                }
                if name.is_empty() {
                    return Err(bad("rt_assessment_empty_fasta"));
                }
                record(
                    &request,
                    reference,
                    &name,
                    &sequence,
                    &tails,
                    &mut cells,
                    &mut report.base_comparisons,
                    &mut validated_sites,
                    &mut ctl,
                )?;
                reader.into_inner().finish()?;
                if reference
                    .circular_record_ids
                    .iter()
                    .any(|id| !seen.contains(id))
                {
                    return Err(bad("rt_assessment_unknown_circular_record"));
                }
                assessment.status = "verified".into();
                for (pi, tail) in tails.iter().enumerate() {
                    for (ni, &n) in SEEDS.iter().enumerate() {
                        let c = &mut cells[pi * SEEDS.len() + ni];
                        if tail.len() < n {
                            c.reason = Some("target_tail_shorter_than_seed".into());
                        } else if c.all_sites.assessed_windows == 0 {
                            c.reason = Some("no_unambiguous_windows_long_enough".into());
                        } else {
                            c.status = "completed".into();
                        }
                    }
                }
            }
            report.references.push(assessment);
            report.cells.extend(cells);
            ctl.completed += 1;
        }
        if validated_sites.len() != request.intended_sites.len() {
            return Err(bad("rt_assessment_intended_record_missing"));
        }
        if ctl.stop(true) {
            return Err(bad("rt_assessment_cancelled"));
        }
        report.content_sha256 = hash(&report)?;
        Ok(report)
    }
}

#[cfg(test)]
#[path = "rt_assessment_tests.rs"]
mod tests;
