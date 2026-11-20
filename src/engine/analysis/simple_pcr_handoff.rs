//! Read-only validation of a saved, tail-free primer pair for PCR form preparation.
//!
//! This does not execute PCR, modify a report or renew specificity acceptance.
//! Tests construct synthetic sequences and reports in memory, without fixtures/I/O.

use super::{EngineError, ErrorCode, GentleEngine, PrimerDesignReport};
use crate::digest_utils::{sha256_prefixed_bytes, sha256_prefixed_str};

/// Point-in-time source binding for an explicit copy into the simple PCR form.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SimplePcrPairPreparation {
    pub report_id: String,
    pub pair_rank: usize,
    pub template: String,
    pub report_sha256: String,
    pub template_sha256: String,
    pub forward_primer: String,
    pub reverse_primer: String,
    pub backend_requested: String,
    pub backend_used: String,
}

fn invalid(message: &str) -> EngineError {
    EngineError {
        code: ErrorCode::InvalidInput,
        message: format!("simple_pcr_handoff:{message}"),
        cause_chain: vec![],
    }
}

impl GentleEngine {
    /// Re-read the saved report and check the exact displayed snapshot, rank and
    /// current template geometry. Legacy reports have no original template digest:
    /// matching current binding sites is not historical whole-template freshness.
    pub fn prepare_simple_pcr_pair(
        &self,
        displayed_report: &PrimerDesignReport,
        pair_rank: usize,
        template: &str,
    ) -> Result<SimplePcrPairPreparation, EngineError> {
        let report = self.get_primer_design_report(&displayed_report.report_id)?;
        let encoded = |report: &PrimerDesignReport| {
            serde_json::to_vec(report).map_err(|_| invalid("invalid_report"))
        };
        let report_bytes = encoded(&report)?;
        if report_bytes != encoded(displayed_report)? {
            return Err(invalid("report_changed; reopen the saved report"));
        }
        if report.schema != super::PRIMER_DESIGN_REPORT_SCHEMA || report.template != template {
            return Err(invalid("report_template_mismatch"));
        }
        let mut candidates = report.pairs.iter().filter(|pair| pair.rank == pair_rank);
        let pair = candidates
            .next()
            .filter(|_| pair_rank > 0)
            .ok_or_else(|| invalid("pair_rank_missing"))?;
        if candidates.next().is_some() {
            return Err(invalid("pair_rank_ambiguous"));
        }
        let dna = self
            .state
            .sequences
            .get(template)
            .ok_or_else(|| invalid("template_missing"))?;
        if dna.is_circular() {
            return Err(invalid("circular_template_unsupported"));
        }
        let sequence = dna.get_forward_string().to_ascii_uppercase();
        for (primer, reverse) in [(&pair.forward, false), (&pair.reverse, true)] {
            if primer.non_annealing_5prime_tail_bp != 0
                || primer.anneal_length_bp != primer.sequence.len()
            {
                return Err(invalid("5prime_tails_require_advanced_pcr"));
            }
            if primer.sequence.is_empty()
                || !primer.sequence.is_ascii()
                || primer.length_bp != primer.sequence.len()
                || primer.end_0based_exclusive.checked_sub(primer.start_0based)
                    != Some(primer.sequence.len())
            {
                return Err(invalid("invalid_primer_geometry"));
            }
            let site = sequence
                .as_bytes()
                .get(primer.start_0based..primer.end_0based_exclusive)
                .ok_or_else(|| invalid("template_binding_changed"))?;
            let binding = if reverse {
                Self::reverse_complement(&primer.sequence.to_ascii_uppercase())
            } else {
                primer.sequence.to_ascii_uppercase()
            };
            if site != binding.as_bytes() {
                return Err(invalid("template_binding_changed"));
            }
        }
        if pair.forward.start_0based > pair.reverse.start_0based {
            return Err(invalid("incompatible_pair_geometry"));
        }
        Ok(SimplePcrPairPreparation {
            report_id: report.report_id.clone(),
            pair_rank,
            template: template.to_string(),
            report_sha256: sha256_prefixed_bytes(&report_bytes),
            template_sha256: sha256_prefixed_str(&sequence),
            forward_primer: pair.forward.sequence.clone(),
            reverse_primer: pair.reverse.sequence.clone(),
            backend_requested: report.backend.requested.clone(),
            backend_used: report.backend.used.clone(),
        })
    }
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;
    use crate::{dna_sequence::DNAsequence, engine::*};

    /// Hand-crafted, non-palindromic reverse primer, unequal lengths and rank 7.
    pub(crate) fn fixture() -> (GentleEngine, PrimerDesignReport) {
        let report = PrimerDesignReport {
            schema: PRIMER_DESIGN_REPORT_SCHEMA.to_string(),
            report_id: "synthetic_pcr_report".to_string(),
            template: "synthetic_template".to_string(),
            pair_count: 1,
            pairs: vec![PrimerDesignPairRecord {
                rank: 7,
                forward: PrimerDesignPrimerRecord {
                    sequence: "AACCG".to_string(),
                    start_0based: 0,
                    end_0based_exclusive: 5,
                    length_bp: 5,
                    anneal_length_bp: 5,
                    ..Default::default()
                },
                reverse: PrimerDesignPrimerRecord {
                    sequence: "TAGTGC".to_string(),
                    start_0based: 14,
                    end_0based_exclusive: 20,
                    length_bp: 6,
                    anneal_length_bp: 6,
                    ..Default::default()
                },
                ..Default::default()
            }],
            backend: PrimerDesignBackendInfo {
                requested: "primer3".to_string(),
                used: "internal".to_string(),
                ..Default::default()
            },
            ..Default::default()
        };
        let mut state = ProjectState::default();
        state.sequences.insert(
            report.template.clone(),
            DNAsequence::from_sequence("AACCGGATTTGCATGCACTA").unwrap(),
        );
        let mut engine = GentleEngine::from_state(state);
        install(&mut engine, &report);
        (engine, report)
    }

    pub(crate) fn install(engine: &mut GentleEngine, report: &PrimerDesignReport) {
        engine.state_mut().metadata.insert(
            PRIMER_DESIGN_REPORTS_METADATA_KEY.to_string(),
            serde_json::json!({
                "schema": "gentle.primer_design_reports.v1",
                "reports": {"synthetic_pcr_report": report}
            }),
        );
    }

    #[test]
    fn simple_pcr_handoff_engine_copies_original_oligos_without_state_changes() {
        let (engine, report) = fixture();
        let before = serde_json::to_value(engine.state()).unwrap();
        let history = serde_json::to_value(engine.history_summary()).unwrap();
        let prepared = engine
            .prepare_simple_pcr_pair(&report, 7, &report.template)
            .unwrap();
        assert_eq!(prepared.forward_primer, "AACCG");
        assert_eq!(prepared.reverse_primer, "TAGTGC");
        assert_eq!(prepared.backend_requested, "primer3");
        assert_eq!(prepared.backend_used, "internal");
        assert_eq!(
            prepared,
            engine
                .prepare_simple_pcr_pair(&report, 7, &report.template)
                .unwrap()
        );
        assert_eq!(before, serde_json::to_value(engine.state()).unwrap());
        assert_eq!(
            history,
            serde_json::to_value(engine.history_summary()).unwrap()
        );
    }

    #[test]
    fn simple_pcr_handoff_engine_refuses_invalid_rank_template_and_deleted_report() {
        let (mut engine, report) = fixture();
        for rank in [0, 1, 8] {
            assert!(
                engine
                    .prepare_simple_pcr_pair(&report, rank, &report.template)
                    .is_err()
            );
        }
        assert!(engine.prepare_simple_pcr_pair(&report, 7, "other").is_err());
        engine.state_mut().sequences.remove(&report.template);
        assert!(
            engine
                .prepare_simple_pcr_pair(&report, 7, &report.template)
                .is_err()
        );
        engine.state_mut().metadata.clear();
        assert!(
            engine
                .prepare_simple_pcr_pair(&report, 7, &report.template)
                .is_err()
        );
    }

    #[test]
    fn simple_pcr_handoff_engine_refuses_tails_malformed_and_duplicate_pairs() {
        for case in 0..5 {
            let (mut engine, mut report) = fixture();
            match case {
                0 => report.pairs[0].forward.non_annealing_5prime_tail_bp = 1,
                1 => report.pairs[0].reverse.anneal_length_bp = 5,
                2 => report.pairs[0].forward.length_bp = 9,
                3 => report.pairs[0].reverse.end_0based_exclusive = usize::MAX,
                _ => report.pairs.push(report.pairs[0].clone()),
            }
            install(&mut engine, &report);
            assert!(
                engine
                    .prepare_simple_pcr_pair(&report, 7, &report.template)
                    .is_err()
            );
        }
    }

    #[test]
    fn simple_pcr_handoff_engine_refuses_changed_report_and_binding() {
        let (mut engine, report) = fixture();
        let mut changed = report.clone();
        changed.backend.used = "different".to_string();
        install(&mut engine, &changed);
        assert!(
            engine
                .prepare_simple_pcr_pair(&report, 7, &report.template)
                .is_err()
        );
        install(&mut engine, &report);
        engine.state_mut().sequences.insert(
            report.template.clone(),
            DNAsequence::from_sequence("TACCGGATTTGCATGCACTA").unwrap(),
        );
        assert!(
            engine
                .prepare_simple_pcr_pair(&report, 7, &report.template)
                .is_err()
        );
    }

    #[test]
    fn simple_pcr_handoff_engine_refuses_circular_templates_and_wrong_report_schema() {
        let (mut engine, mut report) = fixture();
        engine
            .state_mut()
            .sequences
            .get_mut(&report.template)
            .unwrap()
            .set_circular(true);
        assert!(
            engine
                .prepare_simple_pcr_pair(&report, 7, &report.template)
                .is_err()
        );
        engine
            .state_mut()
            .sequences
            .get_mut(&report.template)
            .unwrap()
            .set_circular(false);
        report.schema = "unrelated.report.v1".to_string();
        install(&mut engine, &report);
        assert!(
            engine
                .prepare_simple_pcr_pair(&report, 7, &report.template)
                .is_err()
        );
    }

    #[test]
    fn simple_pcr_handoff_engine_keeps_basic_pcr_support_for_overlapping_binding_sites() {
        let (mut engine, mut report) = fixture();
        report.pairs[0].reverse.sequence = "AATCCG".to_string();
        report.pairs[0].reverse.start_0based = 3;
        report.pairs[0].reverse.end_0based_exclusive = 9;
        install(&mut engine, &report);
        let prepared = engine
            .prepare_simple_pcr_pair(&report, 7, &report.template)
            .unwrap();
        assert_eq!(prepared.reverse_primer, "AATCCG");
    }
}
