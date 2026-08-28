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
}
