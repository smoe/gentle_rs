//! Shared deterministic analyses for multi-promoter reporter panels.
//!
//! This module composes existing vector-site ranking and restriction-site
//! scanning. It does not introduce a second enzyme catalog or cloning model.

use super::*;

impl GentleEngine {
    /// Select one directional restriction pair that is usable across every
    /// panel insert, or return an explicit Gibson fallback.
    pub fn restriction_cloning_panel_strategy(
        &self,
        destination_vector_seq_id: &str,
        insert_seq_ids: &[String],
    ) -> Result<PromoterReporterPanelCloningStrategyReport, EngineError> {
        let destination_vector_seq_id = destination_vector_seq_id.trim();
        if destination_vector_seq_id.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "Panel cloning strategy requires a destination vector seq_id".to_string(),
                cause_chain: vec![],
            });
        }
        if insert_seq_ids.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "Panel cloning strategy requires at least one insert seq_id".to_string(),
                cause_chain: vec![],
            });
        }

        let mut normalized_insert_ids = Vec::with_capacity(insert_seq_ids.len());
        let mut seen = HashSet::new();
        for raw in insert_seq_ids {
            let insert_seq_id = raw.trim();
            if insert_seq_id.is_empty() {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: "Panel cloning strategy insert seq_ids must be non-empty".to_string(),
                    cause_chain: vec![],
                });
            }
            if insert_seq_id == destination_vector_seq_id {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!(
                        "Destination vector '{}' cannot also be a panel insert",
                        destination_vector_seq_id
                    ),
                    cause_chain: vec![],
                });
            }
            if !seen.insert(insert_seq_id.to_string()) {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!("Duplicate panel insert seq_id '{}'", insert_seq_id),
                    cause_chain: vec![],
                });
            }
            if !self.state.sequences.contains_key(insert_seq_id) {
                return Err(EngineError {
                    code: ErrorCode::NotFound,
                    message: format!("Panel insert sequence '{}' not found", insert_seq_id),
                    cause_chain: vec![],
                });
            }
            normalized_insert_ids.push(insert_seq_id.to_string());
        }

        let vector_suggestions =
            self.restriction_cloning_vector_enzyme_suggestions(destination_vector_seq_id)?;
        let candidate_enzymes = vector_suggestions
            .recommended_directed_pairs
            .iter()
            .flat_map(|pair| [&pair.forward_enzyme, &pair.reverse_enzyme])
            .cloned()
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect::<Vec<_>>();

        let mut insert_site_summaries = Vec::with_capacity(normalized_insert_ids.len());
        for insert_seq_id in &normalized_insert_ids {
            let insert_length_bp = self
                .state
                .sequences
                .get(insert_seq_id)
                .map(DNAsequence::len)
                .unwrap_or_default();
            let mut physical_sites = BTreeSet::new();
            if !candidate_enzymes.is_empty() {
                let scan = self.find_restriction_sites(
                    SequenceScanTarget::SeqId {
                        seq_id: insert_seq_id.clone(),
                        span_start_0based: None,
                        span_end_0based_exclusive: None,
                    },
                    &candidate_enzymes,
                    None,
                    false,
                    None,
                    None,
                )?;
                for row in scan.rows {
                    physical_sites.insert((
                        row.enzyme_name,
                        row.recognition_start_0based,
                        row.recognition_end_0based_exclusive,
                    ));
                }
            }
            let mut site_count_by_enzyme = candidate_enzymes
                .iter()
                .cloned()
                .map(|enzyme| (enzyme, 0usize))
                .collect::<BTreeMap<_, _>>();
            for (enzyme, _, _) in physical_sites {
                *site_count_by_enzyme.entry(enzyme).or_default() += 1;
            }
            insert_site_summaries.push(PromoterReporterPanelInsertSiteSummary {
                insert_seq_id: insert_seq_id.clone(),
                insert_length_bp,
                site_count_by_enzyme,
            });
        }

        let mut pair_evaluations =
            Vec::with_capacity(vector_suggestions.recommended_directed_pairs.len());
        let mut selected_pair = None;
        for (pair_index, pair) in vector_suggestions
            .recommended_directed_pairs
            .iter()
            .enumerate()
        {
            let mut blockers = vec![];
            for insert in &insert_site_summaries {
                for enzyme in [&pair.forward_enzyme, &pair.reverse_enzyme] {
                    let site_count = insert
                        .site_count_by_enzyme
                        .get(enzyme)
                        .copied()
                        .unwrap_or_default();
                    if site_count > 0 {
                        blockers.push(PromoterReporterPanelCloningBlocker {
                            insert_seq_id: insert.insert_seq_id.clone(),
                            enzyme: enzyme.clone(),
                            site_count,
                        });
                    }
                }
            }
            let compatible = blockers.is_empty();
            if compatible && selected_pair.is_none() {
                selected_pair = Some(pair.clone());
            }
            pair_evaluations.push(PromoterReporterPanelPairEvaluation {
                rank_1based: pair_index + 1,
                pair: pair.clone(),
                compatible,
                blockers,
            });
        }

        let (strategy, fallback_reason) = if selected_pair.is_some() {
            (
                PromoterReporterPanelCloningStrategy::DirectionalRestriction,
                None,
            )
        } else if vector_suggestions.recommended_directed_pairs.is_empty() {
            (
                PromoterReporterPanelCloningStrategy::Gibson,
                Some(format!(
                    "Destination vector '{}' does not expose two ordered unique restriction sites; use Gibson assembly",
                    destination_vector_seq_id
                )),
            )
        } else {
            (
                PromoterReporterPanelCloningStrategy::Gibson,
                Some(
                    "Every vector-ranked directional restriction pair is blocked by at least one internal site across the panel; use Gibson assembly"
                        .to_string(),
                ),
            )
        };

        Ok(PromoterReporterPanelCloningStrategyReport {
            schema: PROMOTER_REPORTER_PANEL_CLONING_STRATEGY_SCHEMA.to_string(),
            destination_vector_seq_id: destination_vector_seq_id.to_string(),
            insert_seq_ids: normalized_insert_ids,
            vector_suggestions,
            insert_site_summaries,
            pair_evaluations,
            strategy,
            selected_pair,
            fallback_reason,
            warnings: vec![
                "Restriction compatibility is a sequence-derived planning screen; verify enzyme availability and reaction conditions before use."
                    .to_string(),
            ],
        })
    }

    /// Apply the fixed p53-family core-disruption rule and audit its local PWM
    /// and selected restriction-site consequences.
    pub fn design_p53_family_motif_disruption(
        &self,
        sequence_text: &str,
        motif_start_0based: usize,
        motif_end_0based_exclusive: usize,
        motif_forward_strand: bool,
        selected_enzymes: &[String],
    ) -> Result<P53FamilyMotifDisruptionReport, EngineError> {
        const STATED_RULE: &str = "For each oriented 10-bp p53-family half-site, replace position 4 C with A and position 7 G with T (1-based within the half-site)";
        const P53_FAMILY_MOTIFS: [&str; 3] = ["MA0861.1", "MA0106.3", "MA0525.2"];

        let normalized_sequence = sequence_text
            .chars()
            .filter(|base| !base.is_ascii_whitespace())
            .map(|base| base.to_ascii_uppercase())
            .collect::<String>();
        if normalized_sequence.is_empty()
            || normalized_sequence
                .bytes()
                .any(|base| !matches!(base, b'A' | b'C' | b'G' | b'T'))
        {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "p53-family motif disruption requires a non-empty A/C/G/T sequence"
                    .to_string(),
                cause_chain: vec![],
            });
        }
        if motif_start_0based >= motif_end_0based_exclusive
            || motif_end_0based_exclusive > normalized_sequence.len()
        {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "p53-family motif interval {}..{} is outside a {} bp sequence",
                    motif_start_0based,
                    motif_end_0based_exclusive,
                    normalized_sequence.len()
                ),
                cause_chain: vec![],
            });
        }
        let motif_length = motif_end_0based_exclusive - motif_start_0based;
        if !matches!(motif_length, 10 | 20) {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "p53-family motif disruption v1 requires one or two contiguous 10-bp half-sites, found {} bp",
                    motif_length
                ),
                cause_chain: vec![],
            });
        }

        let original_motif_source =
            normalized_sequence[motif_start_0based..motif_end_0based_exclusive].to_string();
        let original_motif_oriented = if motif_forward_strand {
            original_motif_source.clone()
        } else {
            Self::reverse_complement(&original_motif_source)
        };
        let oriented_bytes = original_motif_oriented.as_bytes();
        for half_site_index in 0..(motif_length / 10) {
            let offset = half_site_index * 10;
            if oriented_bytes[offset + 3] != b'C' || oriented_bytes[offset + 6] != b'G' {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!(
                        "p53-family half-site {} does not match the stated-rule core: expected C at oriented position {} and G at position {}",
                        half_site_index + 1,
                        offset + 4,
                        offset + 7
                    ),
                    cause_chain: vec![],
                });
            }
        }

        let mut mutant_bytes = normalized_sequence.as_bytes().to_vec();
        let mut changes = vec![];
        for half_site_index in 0..(motif_length / 10) {
            for (within_half_0based, oriented_replacement) in [(3usize, b'A'), (6usize, b'T')] {
                let oriented_position_0based = half_site_index * 10 + within_half_0based;
                let source_position_0based = if motif_forward_strand {
                    motif_start_0based + oriented_position_0based
                } else {
                    motif_end_0based_exclusive - 1 - oriented_position_0based
                };
                let source_replacement = if motif_forward_strand {
                    oriented_replacement
                } else {
                    IupacCode::letter_complement(oriented_replacement).to_ascii_uppercase()
                };
                let original_base = mutant_bytes[source_position_0based];
                mutant_bytes[source_position_0based] = source_replacement;
                changes.push(P53FamilyMotifBaseChange {
                    source_position_0based,
                    oriented_position_0based,
                    half_site_index_1based: half_site_index + 1,
                    original_base: original_base as char,
                    replacement_base: source_replacement as char,
                    rule: if within_half_0based == 3 {
                        "oriented_core_C_to_A".to_string()
                    } else {
                        "oriented_core_G_to_T".to_string()
                    },
                });
            }
        }
        changes.sort_by_key(|change| change.source_position_0based);
        let mutant_sequence = String::from_utf8(mutant_bytes).map_err(|error| EngineError {
            code: ErrorCode::Internal,
            message: format!("Could not encode p53-family motif mutant: {error}"),
            cause_chain: vec![],
        })?;
        let mutant_motif_source =
            mutant_sequence[motif_start_0based..motif_end_0based_exclusive].to_string();
        let mutant_motif_oriented = if motif_forward_strand {
            mutant_motif_source.clone()
        } else {
            Self::reverse_complement(&mutant_motif_source)
        };
        let changed_positions = changes
            .iter()
            .map(|change| change.source_position_0based)
            .collect::<BTreeSet<_>>();

        let motif_tokens = P53_FAMILY_MOTIFS
            .iter()
            .map(|motif| motif.to_string())
            .collect::<Vec<_>>();
        let wild_type_pwm = self.scan_tfbs_hits(
            SequenceScanTarget::InlineSequence {
                sequence_text: normalized_sequence.clone(),
                topology: InlineSequenceTopology::Linear,
                id_hint: Some("p53_family_wild_type".to_string()),
                span_start_0based: None,
                span_end_0based_exclusive: None,
            },
            &motif_tokens,
            None,
            Some(0.0),
            &[],
            None,
            None,
            None,
        )?;
        let mutant_pwm = self.scan_tfbs_hits(
            SequenceScanTarget::InlineSequence {
                sequence_text: mutant_sequence.clone(),
                topology: InlineSequenceTopology::Linear,
                id_hint: Some("p53_family_mutant".to_string()),
                span_start_0based: None,
                span_end_0based_exclusive: None,
            },
            &motif_tokens,
            None,
            Some(0.0),
            &[],
            None,
            None,
            None,
        )?;
        let mut pwm_audits = vec![];
        for tf_id in P53_FAMILY_MOTIFS {
            let wild_type_scores = Self::p53_family_affected_pwm_scores(
                &wild_type_pwm.rows,
                tf_id,
                &changed_positions,
            );
            let mutant_scores =
                Self::p53_family_affected_pwm_scores(&mutant_pwm.rows, tf_id, &changed_positions);
            let keys = wild_type_scores
                .keys()
                .chain(mutant_scores.keys())
                .copied()
                .collect::<BTreeSet<_>>();
            let stronger_mutant_window_count = keys
                .iter()
                .filter(|key| {
                    mutant_scores.get(key).copied().unwrap_or(f64::NEG_INFINITY)
                        > wild_type_scores
                            .get(key)
                            .copied()
                            .unwrap_or(f64::NEG_INFINITY)
                            + 1e-9
                })
                .count();
            let wild_type_max_llr_bits = Self::finite_max(wild_type_scores.values().copied());
            let mutant_max_llr_bits = Self::finite_max(mutant_scores.values().copied());
            let max_llr_delta_bits = wild_type_max_llr_bits
                .zip(mutant_max_llr_bits)
                .map(|(wild_type, mutant)| mutant - wild_type);
            let tf_name = wild_type_pwm
                .rows
                .iter()
                .chain(mutant_pwm.rows.iter())
                .find(|row| row.tf_id == tf_id)
                .and_then(|row| row.tf_name.clone());
            pwm_audits.push(P53FamilyMotifPwmAudit {
                tf_id: tf_id.to_string(),
                tf_name,
                windows_compared: keys.len(),
                wild_type_max_llr_bits,
                mutant_max_llr_bits,
                max_llr_delta_bits,
                stronger_mutant_window_count,
            });
        }
        let no_stronger_p53_family_hit_near_edit = pwm_audits
            .iter()
            .all(|audit| audit.stronger_mutant_window_count == 0);

        let (restriction_site_audits, creates_selected_restriction_site) = self
            .p53_family_restriction_site_audits(
                &normalized_sequence,
                &mutant_sequence,
                selected_enzymes,
            )?;

        Ok(P53FamilyMotifDisruptionReport {
            schema: P53_FAMILY_MOTIF_DISRUPTION_SCHEMA.to_string(),
            stated_rule: STATED_RULE.to_string(),
            source_sequence_sha256: sha256_prefixed_str(&normalized_sequence),
            motif_start_0based,
            motif_end_0based_exclusive,
            motif_forward_strand,
            original_motif_source,
            original_motif_oriented,
            mutant_motif_source,
            mutant_motif_oriented,
            length_preserved: mutant_sequence.len() == normalized_sequence.len(),
            mutant_sequence,
            changes,
            pwm_audits,
            no_stronger_p53_family_hit_near_edit,
            restriction_site_audits,
            creates_selected_restriction_site,
            warnings: vec![
                "This is a sequence-only disruption proposal, not evidence of lost binding, occupancy, transcriptional effect, or reporter activity."
                    .to_string(),
                "PWM comparisons are post-edit audits of the stated rule; they do not choose or optimize the substitutions."
                    .to_string(),
            ],
        })
    }

    fn p53_family_affected_pwm_scores(
        rows: &[TfbsHitScanRow],
        tf_id: &str,
        changed_positions: &BTreeSet<usize>,
    ) -> BTreeMap<(usize, bool), f64> {
        rows.iter()
            .filter(|row| row.tf_id == tf_id && !row.wraps_origin)
            .filter(|row| {
                changed_positions.iter().any(|position| {
                    row.match_start_0based <= *position
                        && *position < row.match_end_0based_exclusive
                })
            })
            .map(|row| ((row.match_start_0based, row.forward_strand), row.llr_bits))
            .collect()
    }

    fn finite_max(values: impl Iterator<Item = f64>) -> Option<f64> {
        values
            .filter(|value| value.is_finite())
            .reduce(|left, right| left.max(right))
    }

    fn p53_family_restriction_site_audits(
        &self,
        wild_type_sequence: &str,
        mutant_sequence: &str,
        selected_enzymes: &[String],
    ) -> Result<(Vec<P53FamilyRestrictionSiteAudit>, bool), EngineError> {
        let mut enzymes = vec![];
        let mut seen = HashSet::new();
        for raw in selected_enzymes {
            let trimmed = raw.trim();
            if !trimmed.is_empty() && seen.insert(trimmed.to_ascii_uppercase()) {
                enzymes.push(trimmed.to_string());
            }
        }
        if enzymes.is_empty() {
            return Ok((vec![], false));
        }
        let scan = |sequence_text: &str, id_hint: &str| {
            self.find_restriction_sites(
                SequenceScanTarget::InlineSequence {
                    sequence_text: sequence_text.to_string(),
                    topology: InlineSequenceTopology::Linear,
                    id_hint: Some(id_hint.to_string()),
                    span_start_0based: None,
                    span_end_0based_exclusive: None,
                },
                &enzymes,
                None,
                false,
                None,
                None,
            )
        };
        let wild_type = scan(wild_type_sequence, "p53_family_wild_type")?;
        let mutant = scan(mutant_sequence, "p53_family_mutant")?;
        let to_sites = |report: RestrictionSiteScanReport| {
            let mut by_enzyme: BTreeMap<String, BTreeSet<P53FamilyRestrictionSiteInterval>> =
                report
                    .enzymes_scanned
                    .iter()
                    .cloned()
                    .map(|enzyme| (enzyme, BTreeSet::new()))
                    .collect();
            for row in report.rows {
                by_enzyme.entry(row.enzyme_name).or_default().insert(
                    P53FamilyRestrictionSiteInterval {
                        start_0based: row.recognition_start_0based,
                        end_0based_exclusive: row.recognition_end_0based_exclusive,
                    },
                );
            }
            by_enzyme
        };
        let wild_type_sites = to_sites(wild_type);
        let mutant_sites = to_sites(mutant);
        let enzyme_names = wild_type_sites
            .keys()
            .chain(mutant_sites.keys())
            .cloned()
            .collect::<BTreeSet<_>>();
        let mut audits = vec![];
        for enzyme in enzyme_names {
            let before = wild_type_sites.get(&enzyme).cloned().unwrap_or_default();
            let after = mutant_sites.get(&enzyme).cloned().unwrap_or_default();
            audits.push(P53FamilyRestrictionSiteAudit {
                enzyme,
                wild_type_site_count: before.len(),
                mutant_site_count: after.len(),
                newly_created_sites: after.difference(&before).cloned().collect(),
                lost_sites: before.difference(&after).cloned().collect(),
            });
        }
        let creates_selected_restriction_site = audits
            .iter()
            .any(|audit| !audit.newly_created_sites.is_empty());
        Ok((audits, creates_selected_restriction_site))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn promoter_reporter_panel_vector(sequence: &str) -> DNAsequence {
        let mut dna = DNAsequence::from_sequence(sequence).expect("vector sequence");
        *dna.restriction_enzymes_mut() = crate::enzymes::active_restriction_enzymes();
        dna.features_mut().push(gb_io::seq::Feature {
            kind: "misc_feature".into(),
            location: gb_io::seq::Location::simple_range(
                0,
                sequence.len().try_into().expect("sequence length"),
            ),
            qualifiers: vec![
                ("label".into(), Some("MCS".to_string())),
                ("note".into(), Some("Multiple cloning site".to_string())),
                (
                    "mcs_expected_sites".into(),
                    Some("EcoRI,HindIII,KpnI".to_string()),
                ),
            ],
        });
        dna.update_computed_features();
        dna
    }

    fn promoter_reporter_panel_engine(insert_a: &str, insert_b: &str) -> GentleEngine {
        let mut state = ProjectState::default();
        state.sequences.insert(
            "vector".to_string(),
            promoter_reporter_panel_vector("AAAAGAATTCGGGGGAAGCTTTTTTGGTACCAAAA"),
        );
        state.sequences.insert(
            "insert_a".to_string(),
            DNAsequence::from_sequence(insert_a).expect("insert A"),
        );
        state.sequences.insert(
            "insert_b".to_string(),
            DNAsequence::from_sequence(insert_b).expect("insert B"),
        );
        GentleEngine::from_state(state)
    }

    #[test]
    fn promoter_reporter_panel_selects_first_pair_clear_across_every_insert() {
        let engine = promoter_reporter_panel_engine(
            "ACGTACGTACGTGAATTCACGTACGTACGT",
            "TGCATGCATGCATGCATGCATGCA",
        );
        let report = engine
            .restriction_cloning_panel_strategy(
                "vector",
                &["insert_a".to_string(), "insert_b".to_string()],
            )
            .expect("panel strategy");

        assert_eq!(
            report.strategy,
            PromoterReporterPanelCloningStrategy::DirectionalRestriction
        );
        let selected = report.selected_pair.expect("selected pair");
        assert_eq!(selected.forward_enzyme, "HindIII");
        assert_eq!(selected.reverse_enzyme, "KpnI");
        assert_eq!(report.pair_evaluations.len(), 3);
        assert!(!report.pair_evaluations[0].compatible);
        assert!(!report.pair_evaluations[1].compatible);
        assert!(report.pair_evaluations[2].compatible);
        assert_eq!(
            report.insert_site_summaries[0]
                .site_count_by_enzyme
                .get("EcoRI"),
            Some(&1)
        );
        assert_eq!(
            report.insert_site_summaries[1]
                .site_count_by_enzyme
                .get("EcoRI"),
            Some(&0)
        );
    }

    #[test]
    fn promoter_reporter_panel_falls_back_to_gibson_when_all_pairs_are_blocked() {
        let engine = promoter_reporter_panel_engine(
            "ACGTGAATTCACGTAAGCTTACGTGGTACCACGT",
            "TGCATGCATGCATGCATGCA",
        );
        let report = engine
            .restriction_cloning_panel_strategy(
                "vector",
                &["insert_a".to_string(), "insert_b".to_string()],
            )
            .expect("panel strategy");

        assert_eq!(
            report.strategy,
            PromoterReporterPanelCloningStrategy::Gibson
        );
        assert!(report.selected_pair.is_none());
        assert!(report.pair_evaluations.iter().all(|row| !row.compatible));
        assert!(
            report
                .fallback_reason
                .as_deref()
                .is_some_and(|reason| reason.contains("internal site"))
        );
    }

    #[test]
    fn promoter_reporter_panel_p53_rule_is_length_preserving_and_audited() {
        let engine = GentleEngine::default();
        let report = engine
            .design_p53_family_motif_disruption(
                "GGGCATGCCCGGGCATGCCC",
                0,
                20,
                true,
                &["KpnI".to_string(), "HindIII".to_string()],
            )
            .expect("p53-family disruption");

        assert_eq!(report.mutant_motif_oriented, "GGGAATTCCCGGGAATTCCC");
        assert_eq!(report.changes.len(), 4);
        assert!(report.length_preserved);
        assert_eq!(report.pwm_audits.len(), 3);
        assert!(report.no_stronger_p53_family_hit_near_edit);
        assert!(!report.creates_selected_restriction_site);
        assert!(
            report
                .restriction_site_audits
                .iter()
                .all(|audit| audit.newly_created_sites.is_empty())
        );
    }

    #[test]
    fn promoter_reporter_panel_p53_rule_maps_reverse_strand_changes_to_source() {
        let engine = GentleEngine::default();
        let source = GentleEngine::reverse_complement("GGGCATGCCCGGGCATGCCC");
        let report = engine
            .design_p53_family_motif_disruption(&source, 0, 20, false, &[])
            .expect("reverse p53-family disruption");

        assert_eq!(report.original_motif_oriented, "GGGCATGCCCGGGCATGCCC");
        assert_eq!(report.mutant_motif_oriented, "GGGAATTCCCGGGAATTCCC");
        assert_eq!(
            report.mutant_sequence,
            GentleEngine::reverse_complement("GGGAATTCCCGGGAATTCCC")
        );
        assert_eq!(
            report
                .changes
                .iter()
                .map(|change| change.source_position_0based)
                .collect::<Vec<_>>(),
            vec![3, 6, 13, 16]
        );
    }

    #[test]
    fn promoter_reporter_panel_p53_rule_rejects_noncanonical_core() {
        let engine = GentleEngine::default();
        let error = engine
            .design_p53_family_motif_disruption("GGGAATGCCC", 0, 10, true, &[])
            .expect_err("noncanonical core must fail");
        assert_eq!(error.code, ErrorCode::InvalidInput);
        assert!(error.message.contains("expected C"));
    }
}
