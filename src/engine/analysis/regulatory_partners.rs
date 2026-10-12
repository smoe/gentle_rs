//! Gene-set regulatory-partner motif tuples and replayable evidence trees.
//!
//! The implementation reuses the shared promoter resolver and exact TFBS hit
//! scanner. It deliberately makes no enrichment or causal claim: v1 is an
//! inspectable promoter-scoped evidence ledger.

use super::*;
use serde::Serialize;
use std::collections::{BTreeMap, BTreeSet};

impl GentleEngine {
    fn regulatory_partner_source_ref<T: Serialize>(
        report: &T,
        report_kind: &str,
        schema: &str,
        op_id: Option<&str>,
        run_id: Option<&str>,
    ) -> Result<RegulatoryPartnerSourceReportRef, EngineError> {
        let bytes = serde_json::to_vec(report).map_err(|error| {
            EngineError::internal(format!(
                "Could not serialize {report_kind} for regulatory-partner provenance: {error}"
            ))
        })?;
        Ok(RegulatoryPartnerSourceReportRef {
            report_kind: report_kind.to_string(),
            schema: schema.to_string(),
            content_sha256: crate::digest_utils::sha256_prefixed_bytes(&bytes),
            op_id: op_id.map(str::to_string),
            run_id: run_id.map(str::to_string),
        })
    }

    fn regulatory_partner_interval_relation(
        anchor: &RegulatoryPartnerMotifHit,
        partner: &RegulatoryPartnerMotifHit,
    ) -> RegulatoryPartnerIntervalRelation {
        let a_start = anchor.promoter_start_0based;
        let a_end = anchor.promoter_end_0based_exclusive;
        let p_start = partner.promoter_start_0based;
        let p_end = partner.promoter_end_0based_exclusive;
        if a_end == p_start || p_end == a_start {
            RegulatoryPartnerIntervalRelation::Abutting
        } else if a_start < p_end && p_start < a_end {
            if (a_start <= p_start && a_end >= p_end) || (p_start <= a_start && p_end >= a_end) {
                RegulatoryPartnerIntervalRelation::Contained
            } else {
                RegulatoryPartnerIntervalRelation::Overlap
            }
        } else {
            RegulatoryPartnerIntervalRelation::Disjoint
        }
    }

    fn regulatory_partner_orientation(
        anchor: &RegulatoryPartnerMotifHit,
        partner: &RegulatoryPartnerMotifHit,
    ) -> RegulatoryPartnerOrientation {
        if anchor.promoter_forward_strand == partner.promoter_forward_strand {
            return RegulatoryPartnerOrientation::Same;
        }
        let anchor_center_x2 = anchor.promoter_start_0based + anchor.promoter_end_0based_exclusive;
        let partner_center_x2 =
            partner.promoter_start_0based + partner.promoter_end_0based_exclusive;
        if anchor_center_x2 == partner_center_x2 {
            return RegulatoryPartnerOrientation::Ambiguous;
        }
        let (left_forward, right_forward) = if anchor_center_x2 < partner_center_x2 {
            (
                anchor.promoter_forward_strand,
                partner.promoter_forward_strand,
            )
        } else {
            (
                partner.promoter_forward_strand,
                anchor.promoter_forward_strand,
            )
        };
        match (left_forward, right_forward) {
            (true, false) => RegulatoryPartnerOrientation::Convergent,
            (false, true) => RegulatoryPartnerOrientation::Divergent,
            _ => RegulatoryPartnerOrientation::Ambiguous,
        }
    }

    fn regulatory_partner_decision_tree() -> RegulatoryPartnerDecisionTree {
        let node = |node_id: &str, label: &str, question: &str, terminal: bool| {
            RegulatoryPartnerDecisionNode {
                node_id: node_id.to_string(),
                label: label.to_string(),
                question: question.to_string(),
                terminal,
            }
        };
        let edge = |from: &str, to: &str, outcome| RegulatoryPartnerDecisionEdge {
            from_node_id: from.to_string(),
            to_node_id: to.to_string(),
            outcome,
        };
        RegulatoryPartnerDecisionTree {
            schema: REGULATORY_PARTNER_DECISION_TREE_SCHEMA.to_string(),
            root_node_id: "promoter_resolved".to_string(),
            nodes: vec![
                node(
                    "promoter_resolved",
                    "Promoter resolved",
                    "Was a strand-aware promoter window and sequence resolved?",
                    false,
                ),
                node(
                    "anchor_motif_present",
                    "Anchor motif",
                    "Does the promoter contain a threshold-passing anchor motif?",
                    false,
                ),
                node(
                    "partner_motif_present",
                    "Partner motif",
                    "Does the promoter contain a threshold-passing partner motif?",
                    false,
                ),
                node(
                    "proximal_tuple_present",
                    "Proximal tuple",
                    "Is at least one anchor-partner tuple within the declared distance?",
                    false,
                ),
                node(
                    "promoter_occupancy",
                    "Promoter occupancy",
                    "Was promoter-level CUT&RUN evidence evaluated and supported?",
                    false,
                ),
                node(
                    "unresolved",
                    "Unresolved",
                    "Promoter could not be evaluated.",
                    true,
                ),
                node(
                    "no_anchor",
                    "No anchor",
                    "No anchor motif passed the threshold.",
                    true,
                ),
                node(
                    "anchor_only",
                    "Anchor only",
                    "Anchor evidence was present without a partner motif.",
                    true,
                ),
                node(
                    "distant_partner",
                    "Distant partner",
                    "Anchor and partner motifs were present, but no tuple met the distance rule.",
                    true,
                ),
                node(
                    "occupancy_supported_tuple",
                    "Occupancy-supported tuple",
                    "A proximal motif tuple and promoter-level occupancy support were present.",
                    true,
                ),
                node(
                    "tuple_without_occupancy_support",
                    "Tuple, zero occupancy support",
                    "A proximal motif tuple was present in an evaluated promoter with zero support.",
                    true,
                ),
                node(
                    "tuple_occupancy_unknown",
                    "Tuple, occupancy unknown",
                    "A proximal motif tuple was present, but occupancy was not evaluated or not required.",
                    true,
                ),
            ],
            edges: vec![
                edge(
                    "promoter_resolved",
                    "anchor_motif_present",
                    RegulatoryPartnerDecisionOutcome::Pass,
                ),
                edge(
                    "promoter_resolved",
                    "unresolved",
                    RegulatoryPartnerDecisionOutcome::Fail,
                ),
                edge(
                    "anchor_motif_present",
                    "partner_motif_present",
                    RegulatoryPartnerDecisionOutcome::Pass,
                ),
                edge(
                    "anchor_motif_present",
                    "no_anchor",
                    RegulatoryPartnerDecisionOutcome::Fail,
                ),
                edge(
                    "partner_motif_present",
                    "proximal_tuple_present",
                    RegulatoryPartnerDecisionOutcome::Pass,
                ),
                edge(
                    "partner_motif_present",
                    "anchor_only",
                    RegulatoryPartnerDecisionOutcome::Fail,
                ),
                edge(
                    "proximal_tuple_present",
                    "promoter_occupancy",
                    RegulatoryPartnerDecisionOutcome::Pass,
                ),
                edge(
                    "proximal_tuple_present",
                    "distant_partner",
                    RegulatoryPartnerDecisionOutcome::Fail,
                ),
                edge(
                    "promoter_occupancy",
                    "occupancy_supported_tuple",
                    RegulatoryPartnerDecisionOutcome::Pass,
                ),
                edge(
                    "promoter_occupancy",
                    "tuple_without_occupancy_support",
                    RegulatoryPartnerDecisionOutcome::Fail,
                ),
                edge(
                    "promoter_occupancy",
                    "tuple_occupancy_unknown",
                    RegulatoryPartnerDecisionOutcome::Unknown,
                ),
            ],
        }
    }

    fn regulatory_partner_terminal_step(
        node_id: &str,
        detail: impl Into<String>,
        evidence_ids: Vec<String>,
    ) -> RegulatoryPartnerDecisionTraceStep {
        RegulatoryPartnerDecisionTraceStep {
            node_id: node_id.to_string(),
            outcome: RegulatoryPartnerDecisionOutcome::Pass,
            detail: detail.into(),
            evidence_ids,
        }
    }

    fn regulatory_partner_trace(
        promoter_resolved: bool,
        unresolved_reason: Option<&str>,
        anchor_hits: &[RegulatoryPartnerMotifHit],
        partner_hits: &[RegulatoryPartnerMotifHit],
        proximal_tuples: &[RegulatoryPartnerTupleRow],
        anchor_mode: RegulatoryPartnerAnchorMode,
        occupancy_state: GeneSetCutRunEvaluationState,
        occupancy_supported: bool,
    ) -> (String, Vec<RegulatoryPartnerDecisionTraceStep>) {
        let mut trace = vec![RegulatoryPartnerDecisionTraceStep {
            node_id: "promoter_resolved".to_string(),
            outcome: if promoter_resolved {
                RegulatoryPartnerDecisionOutcome::Pass
            } else {
                RegulatoryPartnerDecisionOutcome::Fail
            },
            detail: unresolved_reason
                .map(str::to_string)
                .unwrap_or_else(|| "Strand-aware promoter sequence resolved.".to_string()),
            evidence_ids: vec![],
        }];
        if !promoter_resolved {
            trace.push(Self::regulatory_partner_terminal_step(
                "unresolved",
                "No sequence-level motif decision was made.",
                vec![],
            ));
            return ("unresolved".to_string(), trace);
        }

        let anchor_ids = anchor_hits
            .iter()
            .map(|hit| hit.hit_id.clone())
            .collect::<Vec<_>>();
        trace.push(RegulatoryPartnerDecisionTraceStep {
            node_id: "anchor_motif_present".to_string(),
            outcome: if anchor_hits.is_empty() {
                RegulatoryPartnerDecisionOutcome::Fail
            } else {
                RegulatoryPartnerDecisionOutcome::Pass
            },
            detail: format!(
                "{} anchor motif hit(s) passed absolute thresholds.",
                anchor_hits.len()
            ),
            evidence_ids: anchor_ids.clone(),
        });
        if anchor_hits.is_empty() {
            trace.push(Self::regulatory_partner_terminal_step(
                "no_anchor",
                "Partner evidence was not promoted to a tuple without an anchor motif.",
                vec![],
            ));
            return ("no_anchor".to_string(), trace);
        }

        let partner_ids = partner_hits
            .iter()
            .map(|hit| hit.hit_id.clone())
            .collect::<Vec<_>>();
        trace.push(RegulatoryPartnerDecisionTraceStep {
            node_id: "partner_motif_present".to_string(),
            outcome: if partner_hits.is_empty() {
                RegulatoryPartnerDecisionOutcome::Fail
            } else {
                RegulatoryPartnerDecisionOutcome::Pass
            },
            detail: format!(
                "{} partner motif hit(s) passed absolute thresholds.",
                partner_hits.len()
            ),
            evidence_ids: partner_ids,
        });
        if partner_hits.is_empty() {
            trace.push(Self::regulatory_partner_terminal_step(
                "anchor_only",
                "Anchor motif evidence was retained without a partner tuple.",
                anchor_ids,
            ));
            return ("anchor_only".to_string(), trace);
        }

        let tuple_ids = proximal_tuples
            .iter()
            .map(|row| row.tuple_id.clone())
            .collect::<Vec<_>>();
        trace.push(RegulatoryPartnerDecisionTraceStep {
            node_id: "proximal_tuple_present".to_string(),
            outcome: if proximal_tuples.is_empty() {
                RegulatoryPartnerDecisionOutcome::Fail
            } else {
                RegulatoryPartnerDecisionOutcome::Pass
            },
            detail: format!(
                "{} tuple(s) met the declared anchor-partner distance.",
                proximal_tuples.len()
            ),
            evidence_ids: tuple_ids.clone(),
        });
        if proximal_tuples.is_empty() {
            trace.push(Self::regulatory_partner_terminal_step(
                "distant_partner",
                "All exact motif intervals remain available in the ledger.",
                vec![],
            ));
            return ("distant_partner".to_string(), trace);
        }

        let (outcome, terminal_class, detail) =
            if anchor_mode == RegulatoryPartnerAnchorMode::MotifOnly {
                (
                    RegulatoryPartnerDecisionOutcome::Unknown,
                    "tuple_occupancy_unknown",
                    "Motif-only mode did not require occupancy evidence.",
                )
            } else if occupancy_state == GeneSetCutRunEvaluationState::Unevaluated {
                (
                    RegulatoryPartnerDecisionOutcome::Unknown,
                    "tuple_occupancy_unknown",
                    "No overlapping prepared CUT&RUN evidence evaluated this promoter.",
                )
            } else if occupancy_supported {
                (
                    RegulatoryPartnerDecisionOutcome::Pass,
                    "occupancy_supported_tuple",
                    "The promoter had at least one CUT&RUN support window.",
                )
            } else {
                (
                    RegulatoryPartnerDecisionOutcome::Fail,
                    "tuple_without_occupancy_support",
                    "The promoter was evaluated and had zero CUT&RUN support windows.",
                )
            };
        trace.push(RegulatoryPartnerDecisionTraceStep {
            node_id: "promoter_occupancy".to_string(),
            outcome,
            detail: detail.to_string(),
            evidence_ids: vec![],
        });
        trace.push(Self::regulatory_partner_terminal_step(
            terminal_class,
            detail,
            tuple_ids,
        ));
        (terminal_class.to_string(), trace)
    }

    #[allow(clippy::too_many_arguments)]
    pub(crate) fn inspect_regulatory_partner_screen(
        &self,
        genome_id: &str,
        resolution: GeneSetResolutionReport,
        cutrun_support: Option<GeneSetCutRunRegulatorySupportReport>,
        anchor_motifs: &[String],
        partner_motifs: &[String],
        min_llr_bits: f64,
        per_motif_thresholds: &[RegulatoryPartnerMotifThreshold],
        max_distance_bp: usize,
        anchor_mode: RegulatoryPartnerAnchorMode,
        upstream_bp: usize,
        downstream_bp: usize,
        genome_catalog_path: Option<&str>,
        cache_dir: Option<&str>,
    ) -> Result<RegulatoryPartnerScreenReport, EngineError> {
        if anchor_motifs.is_empty() {
            return Err(EngineError::invalid_input(
                "InspectRegulatoryPartnerScreen requires at least one anchor motif",
            ));
        }
        if partner_motifs.is_empty() {
            return Err(EngineError::invalid_input(
                "InspectRegulatoryPartnerScreen requires at least one partner motif",
            ));
        }
        if !min_llr_bits.is_finite() {
            return Err(EngineError::invalid_input(
                "InspectRegulatoryPartnerScreen min_llr_bits must be finite",
            ));
        }
        if max_distance_bp == 0 {
            return Err(EngineError::invalid_input(
                "InspectRegulatoryPartnerScreen max_distance_bp must be >= 1",
            ));
        }
        if let Some(row) = per_motif_thresholds
            .iter()
            .find(|row| row.motif.trim().is_empty() || !row.min_llr_bits.is_finite())
        {
            return Err(EngineError::invalid_input(format!(
                "Invalid per-motif threshold for '{}': motif must be non-empty and min_llr_bits finite",
                row.motif
            )));
        }

        let promoter_cohort = self.build_gene_set_promoter_cohort(
            genome_id,
            resolution,
            GeneSetCohortRelationship::Unspecified,
            upstream_bp,
            downstream_bp,
            genome_catalog_path,
            cache_dir,
        )?;
        let effective_catalog_path =
            genome_catalog_path.unwrap_or(crate::genomes::default_catalog_discovery_token(false));
        let (catalog, _) = Self::open_reference_genome_catalog(Some(effective_catalog_path))?;
        let scan_thresholds = per_motif_thresholds
            .iter()
            .map(|row| TfThresholdOverride {
                tf: row.motif.clone(),
                min_llr_bits: Some(row.min_llr_bits),
                min_llr_quantile: None,
            })
            .collect::<Vec<_>>();

        let mut source_reports = vec![Self::regulatory_partner_source_ref(
            &promoter_cohort.gene_set_resolution,
            "gene_set_resolution",
            &promoter_cohort.gene_set_resolution.schema,
            promoter_cohort.gene_set_resolution.op_id.as_deref(),
            promoter_cohort.gene_set_resolution.run_id.as_deref(),
        )?];
        if let Some(report) = cutrun_support.as_ref() {
            if report.genome_id != genome_id {
                return Err(EngineError::invalid_input(format!(
                    "CUT&RUN gene-set report genome '{}' does not match requested genome '{}'",
                    report.genome_id, genome_id
                )));
            }
            source_reports.push(Self::regulatory_partner_source_ref(
                report,
                "gene_set_cutrun_regulatory_support",
                &report.schema,
                report.op_id.as_deref(),
                report.run_id.as_deref(),
            )?);
        }
        let occupancy_by_member = cutrun_support
            .as_ref()
            .map(|report| {
                report
                    .member_support
                    .iter()
                    .map(|row| (row.member_dedup_key.clone(), row))
                    .collect::<BTreeMap<_, _>>()
            })
            .unwrap_or_default();

        let mut warnings = promoter_cohort.warnings.clone();
        if cutrun_support.is_none()
            && anchor_mode == RegulatoryPartnerAnchorMode::OccupancyPreferred
        {
            warnings.push(
                "Occupancy-preferred mode received no gene-set CUT&RUN support report; occupancy branches remain unknown."
                    .to_string(),
            );
        }
        if cutrun_support.is_some() {
            warnings.push(
                "CUT&RUN support is promoter-level. It does not prove that an occupancy interval overlaps a particular motif hit, and the source report's dataset target identity must be reviewed from its listed dataset ids."
                    .to_string(),
            );
        }

        let mut motif_hits = vec![];
        let mut tuples = vec![];
        let mut genes = vec![];
        let mut window_keys = BTreeSet::new();

        for window in &promoter_cohort.windows {
            window_keys.insert(window.member_dedup_key.clone());
            let raw_sequence = catalog
                .get_sequence_region_with_cache(
                    genome_id,
                    &window.chromosome,
                    window.promoter_start_1based,
                    window.promoter_end_1based,
                    cache_dir,
                )
                .map_err(|error| {
                    EngineError::new(
                        ErrorCode::NotFound,
                        format!(
                            "Could not load promoter {}:{}-{} for '{}': {error}",
                            window.chromosome,
                            window.promoter_start_1based,
                            window.promoter_end_1based,
                            window.symbol
                        ),
                    )
                })?;
            let strand = window.strand.chars().next();
            let promoter_sequence = Self::promoter_aligned_sequence(&raw_sequence, strand);
            let tss_position_0based = Self::promoter_oriented_tss_position_0based(
                promoter_sequence.len(),
                window.promoter_start_1based,
                window.tss_1based,
                strand,
            )
            .ok_or_else(|| {
                EngineError::internal(format!(
                    "Could not place TSS for '{}' into its promoter sequence",
                    window.symbol
                ))
            })?;
            let target = |role: &str| SequenceScanTarget::InlineSequence {
                sequence_text: promoter_sequence.clone(),
                topology: InlineSequenceTopology::Linear,
                id_hint: Some(format!("{}:{role}", window.member_dedup_key)),
                span_start_0based: None,
                span_end_0based_exclusive: None,
            };
            let anchor_scan = self.scan_tfbs_hits(
                target("anchor"),
                anchor_motifs,
                Some(min_llr_bits),
                None,
                &scan_thresholds,
                None,
                None,
                None,
            )?;
            let partner_scan = self.scan_tfbs_hits(
                target("partner"),
                partner_motifs,
                Some(min_llr_bits),
                None,
                &scan_thresholds,
                None,
                None,
                None,
            )?;
            if anchor_scan.truncated_at_max_hits || partner_scan.truncated_at_max_hits {
                return Err(EngineError::internal(
                    "Regulatory-partner authoritative hit scan was unexpectedly truncated",
                ));
            }
            let overlapping_tf_ids = anchor_scan
                .motifs_scanned
                .iter()
                .filter(|tf_id| partner_scan.motifs_scanned.contains(tf_id))
                .cloned()
                .collect::<BTreeSet<_>>();
            if !overlapping_tf_ids.is_empty() {
                return Err(EngineError::invalid_input(format!(
                    "Anchor and partner motif sets resolve to the same motif id(s): {}",
                    overlapping_tf_ids
                        .into_iter()
                        .collect::<Vec<_>>()
                        .join(", ")
                )));
            }

            let occupancy = occupancy_by_member.get(&window.member_dedup_key).copied();
            let occupancy_state = occupancy
                .map(|row| row.evaluation_state)
                .unwrap_or(GeneSetCutRunEvaluationState::Unevaluated);
            let occupancy_supported = occupancy.is_some_and(|row| row.support_window_count > 0);

            let mut gene_hit_ids = vec![];
            let mut gene_anchor_hits = vec![];
            let mut gene_partner_hits = vec![];
            let mut seen_hits = BTreeSet::new();
            for (role, rows) in [
                (RegulatoryPartnerMotifRole::Anchor, anchor_scan.rows),
                (RegulatoryPartnerMotifRole::Partner, partner_scan.rows),
            ] {
                for row in rows {
                    let identity = (
                        role,
                        row.tf_id.clone(),
                        row.match_start_0based,
                        row.match_end_0based_exclusive,
                        row.forward_strand,
                    );
                    if !seen_hits.insert(identity) {
                        continue;
                    }
                    let local_end_base = row.match_end_0based_exclusive.saturating_sub(1);
                    let genomic_left = Self::promoter_local_position_to_genomic_1based(
                        strand,
                        window.promoter_start_1based,
                        window.promoter_end_1based,
                        promoter_sequence.len(),
                        row.match_start_0based,
                    )
                    .ok_or_else(|| EngineError::internal("Could not map motif start to genome"))?;
                    let genomic_right = Self::promoter_local_position_to_genomic_1based(
                        strand,
                        window.promoter_start_1based,
                        window.promoter_end_1based,
                        promoter_sequence.len(),
                        local_end_base,
                    )
                    .ok_or_else(|| EngineError::internal("Could not map motif end to genome"))?;
                    let role_label = match role {
                        RegulatoryPartnerMotifRole::Anchor => "anchor",
                        RegulatoryPartnerMotifRole::Partner => "partner",
                    };
                    let strand_label = if row.forward_strand { "f" } else { "r" };
                    let hit_id = format!(
                        "hit:{}:{}:{}:{}:{}",
                        window.member_dedup_key,
                        role_label,
                        row.tf_id,
                        row.match_start_0based,
                        strand_label
                    );
                    let hit = RegulatoryPartnerMotifHit {
                        hit_id: hit_id.clone(),
                        member_dedup_key: window.member_dedup_key.clone(),
                        gene_symbol: window.symbol.clone(),
                        role,
                        tf_id: row.tf_id,
                        tf_name: row.tf_name,
                        motif_consensus_iupac: row.motif_consensus_iupac,
                        motif_length_bp: row.motif_length_bp,
                        promoter_start_0based: row.match_start_0based,
                        promoter_end_0based_exclusive: row.match_end_0based_exclusive,
                        genomic_start_1based: genomic_left.min(genomic_right),
                        genomic_end_1based: genomic_left.max(genomic_right),
                        promoter_forward_strand: row.forward_strand,
                        genomic_forward_strand: if strand == Some('-') {
                            !row.forward_strand
                        } else {
                            row.forward_strand
                        },
                        matched_sequence: row.matched_sequence,
                        llr_bits: row.llr_bits,
                        llr_quantile_within_promoter: row.llr_quantile,
                        signed_start_to_tss_bp: row.match_start_0based as i64
                            - tss_position_0based as i64,
                        signed_center_to_tss_bp: ((row.match_start_0based
                            + row.match_end_0based_exclusive)
                            as i64
                            - (2 * tss_position_0based) as i64)
                            / 2,
                    };
                    gene_hit_ids.push(hit_id);
                    match role {
                        RegulatoryPartnerMotifRole::Anchor => gene_anchor_hits.push(hit.clone()),
                        RegulatoryPartnerMotifRole::Partner => gene_partner_hits.push(hit.clone()),
                    }
                    motif_hits.push(hit);
                }
            }

            let mut gene_tuple_ids = vec![];
            let mut gene_proximal_tuples = vec![];
            for anchor in &gene_anchor_hits {
                for partner in &gene_partner_hits {
                    let center_delta_x2 = (partner.promoter_start_0based
                        + partner.promoter_end_0based_exclusive)
                        as i64
                        - (anchor.promoter_start_0based + anchor.promoter_end_0based_exclusive)
                            as i64;
                    let center_distance = center_delta_x2 as f64 / 2.0;
                    let absolute_distance = (center_distance.abs().round()) as usize;
                    let edge_distance = if partner.promoter_start_0based
                        >= anchor.promoter_end_0based_exclusive
                    {
                        (partner.promoter_start_0based - anchor.promoter_end_0based_exclusive)
                            as i64
                    } else if partner.promoter_end_0based_exclusive <= anchor.promoter_start_0based
                    {
                        -((anchor.promoter_start_0based - partner.promoter_end_0based_exclusive)
                            as i64)
                    } else {
                        0
                    };
                    let tuple_id = format!(
                        "tuple:{}:{}:{}",
                        window.member_dedup_key,
                        anchor.hit_id.trim_start_matches("hit:"),
                        partner.hit_id.trim_start_matches("hit:")
                    );
                    let within_requested_distance = absolute_distance <= max_distance_bp;
                    let tuple = RegulatoryPartnerTupleRow {
                        tuple_id: tuple_id.clone(),
                        member_dedup_key: window.member_dedup_key.clone(),
                        gene_symbol: window.symbol.clone(),
                        anchor_hit_id: anchor.hit_id.clone(),
                        partner_hit_id: partner.hit_id.clone(),
                        signed_anchor_to_partner_center_distance_bp: center_distance,
                        signed_anchor_to_partner_edge_distance_bp: edge_distance,
                        absolute_anchor_to_partner_distance_bp: absolute_distance,
                        interval_relation: Self::regulatory_partner_interval_relation(
                            anchor, partner,
                        ),
                        orientation: Self::regulatory_partner_orientation(anchor, partner),
                        within_requested_distance,
                        occupancy_evaluation_state: occupancy_state,
                        promoter_occupancy_supported: occupancy_supported,
                        occupancy_dataset_ids: occupancy
                            .map(|row| row.contributing_dataset_ids.clone())
                            .unwrap_or_default(),
                        occupancy_read_report_ids: occupancy
                            .map(|row| row.contributing_read_report_ids.clone())
                            .unwrap_or_default(),
                    };
                    gene_tuple_ids.push(tuple_id);
                    if within_requested_distance {
                        gene_proximal_tuples.push(tuple.clone());
                    }
                    tuples.push(tuple);
                }
            }
            let (terminal_class, trace) = Self::regulatory_partner_trace(
                true,
                None,
                &gene_anchor_hits,
                &gene_partner_hits,
                &gene_proximal_tuples,
                anchor_mode,
                occupancy_state,
                occupancy_supported,
            );
            genes.push(RegulatoryPartnerGeneRow {
                member_dedup_key: window.member_dedup_key.clone(),
                gene_symbol: window.symbol.clone(),
                gene_id: window.gene_id.clone(),
                promoter_resolved: true,
                unresolved_reason: None,
                transcript_id: Some(window.transcript_id.clone()),
                chromosome: Some(window.chromosome.clone()),
                strand: Some(window.strand.clone()),
                promoter_start_1based: Some(window.promoter_start_1based),
                promoter_end_1based: Some(window.promoter_end_1based),
                tss_1based: Some(window.tss_1based),
                tss_position_0based: Some(tss_position_0based),
                promoter_sequence_5to3: Some(promoter_sequence),
                motif_hit_ids: gene_hit_ids,
                tuple_ids: gene_tuple_ids,
                occupancy_evaluation_state: occupancy_state,
                promoter_occupancy_supported: occupancy_supported,
                terminal_class,
                trace,
            });
        }

        for member in &promoter_cohort.gene_set_resolution.resolved_members {
            if window_keys.contains(&member.dedup_key) {
                continue;
            }
            let unresolved_reason = promoter_cohort
                .unresolved_members
                .iter()
                .find(|row| row.source_id.as_deref() == Some(member.dedup_key.as_str()))
                .map(|row| row.reason.clone())
                .unwrap_or_else(|| "Promoter window was not returned.".to_string());
            let (terminal_class, trace) = Self::regulatory_partner_trace(
                false,
                Some(&unresolved_reason),
                &[],
                &[],
                &[],
                anchor_mode,
                GeneSetCutRunEvaluationState::Unevaluated,
                false,
            );
            genes.push(RegulatoryPartnerGeneRow {
                member_dedup_key: member.dedup_key.clone(),
                gene_symbol: member.symbol.clone(),
                gene_id: member.gene_id.clone(),
                unresolved_reason: Some(unresolved_reason),
                terminal_class,
                trace,
                ..RegulatoryPartnerGeneRow::default()
            });
        }
        for (index, member) in promoter_cohort
            .gene_set_resolution
            .unresolved_members
            .iter()
            .enumerate()
        {
            let dedup_key = format!("unresolved:{index}:{}", member.query);
            let (terminal_class, trace) = Self::regulatory_partner_trace(
                false,
                Some(&member.reason),
                &[],
                &[],
                &[],
                anchor_mode,
                GeneSetCutRunEvaluationState::Unevaluated,
                false,
            );
            genes.push(RegulatoryPartnerGeneRow {
                member_dedup_key: dedup_key,
                gene_symbol: member.query.clone(),
                unresolved_reason: Some(member.reason.clone()),
                terminal_class,
                trace,
                ..RegulatoryPartnerGeneRow::default()
            });
        }
        genes.sort_by(|left, right| {
            left.gene_symbol
                .to_ascii_lowercase()
                .cmp(&right.gene_symbol.to_ascii_lowercase())
                .then(left.member_dedup_key.cmp(&right.member_dedup_key))
        });
        motif_hits.sort_by(|left, right| {
            left.gene_symbol
                .to_ascii_lowercase()
                .cmp(&right.gene_symbol.to_ascii_lowercase())
                .then(left.promoter_start_0based.cmp(&right.promoter_start_0based))
                .then(left.role.cmp(&right.role))
                .then(left.tf_id.cmp(&right.tf_id))
                .then(
                    left.promoter_forward_strand
                        .cmp(&right.promoter_forward_strand),
                )
        });
        tuples.sort_by(|left, right| {
            left.gene_symbol
                .to_ascii_lowercase()
                .cmp(&right.gene_symbol.to_ascii_lowercase())
                .then(
                    left.absolute_anchor_to_partner_distance_bp
                        .cmp(&right.absolute_anchor_to_partner_distance_bp),
                )
                .then(left.tuple_id.cmp(&right.tuple_id))
        });

        let request = RegulatoryPartnerScreenRequest {
            genome_id: genome_id.to_string(),
            anchor_motifs: anchor_motifs.to_vec(),
            partner_motifs: partner_motifs.to_vec(),
            default_min_llr_bits: min_llr_bits,
            per_motif_thresholds: per_motif_thresholds.to_vec(),
            max_anchor_partner_distance_bp: max_distance_bp,
            upstream_bp,
            downstream_bp,
            anchor_mode,
        };
        let ledger = RegulatoryPartnerTupleLedger {
            schema: REGULATORY_PARTNER_TUPLE_LEDGER_SCHEMA.to_string(),
            request,
            promoter_cohort: Box::new(promoter_cohort.clone()),
            source_reports,
            requested_member_count: promoter_cohort.gene_set_resolution.requested_member_count,
            resolved_promoter_count: promoter_cohort.returned_window_count,
            motif_hit_count: motif_hits.len(),
            tuple_count: tuples.len(),
            genes,
            motif_hits,
            tuples,
            unresolved_members: promoter_cohort.unresolved_members.clone(),
        };
        Ok(RegulatoryPartnerScreenReport {
            schema: REGULATORY_PARTNER_SCREEN_SCHEMA.to_string(),
            generated_at_unix_ms: Self::now_unix_ms(),
            op_id: None,
            run_id: None,
            ledger,
            decision_tree: Self::regulatory_partner_decision_tree(),
            warnings,
            non_claims: vec![
                "Motif hits nominate DNA-binding-family evidence; shared PWMs do not identify a unique protein or cofactor."
                    .to_string(),
                "This v1 screen is promoter-scoped and does not assess distal regulatory elements."
                    .to_string(),
                "Co-occurrence, distance, and promoter occupancy are association evidence, not proof of direct or causal co-regulation."
                    .to_string(),
                "Early versus secondary response timing is not assessed without time-resolved expression evidence."
                    .to_string(),
            ],
        })
    }
}
