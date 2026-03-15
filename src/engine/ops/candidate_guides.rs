//! Candidate-set generation/scoring and guide-design operation handlers.

use super::*;

impl GentleEngine {
    pub(super) fn op_generate_candidate_set(
        &mut self,
        set_name: String,
        seq_id: SeqId,
        length_bp: usize,
        step_bp: usize,
        feature_kinds: Vec<String>,
        feature_label_regex: Option<String>,
        max_distance_bp: Option<usize>,
        feature_geometry_mode: Option<CandidateFeatureGeometryMode>,
        feature_boundary_mode: Option<CandidateFeatureBoundaryMode>,
        feature_strand_relation: Option<CandidateFeatureStrandRelation>,
        limit: Option<usize>,
        result: &mut OpResult,
    ) -> Result<(), EngineError> {
        let set_name = Self::normalize_candidate_set_name(&set_name)?;
        if length_bp == 0 {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "GenerateCandidateSet requires length_bp >= 1".to_string(),
            });
        }
        if step_bp == 0 {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "GenerateCandidateSet requires step_bp >= 1".to_string(),
            });
        }
        let limit = limit.unwrap_or(self.max_fragments_per_container());
        if limit == 0 {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "GenerateCandidateSet requires limit >= 1".to_string(),
            });
        }
        let dna = self
            .state
            .sequences
            .get(&seq_id)
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!("Sequence '{}' not found", seq_id),
            })?;
        if dna.len() < length_bp {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "GenerateCandidateSet length_bp {} exceeds sequence '{}' length {}",
                    length_bp,
                    seq_id,
                    dna.len()
                ),
            });
        }
        let label_regex =
            Self::compile_optional_regex(&feature_label_regex, "feature_label_regex")?;
        let mut kind_filter_upper = feature_kinds
            .iter()
            .map(|kind| kind.trim().to_ascii_uppercase())
            .filter(|kind| !kind.is_empty())
            .collect::<Vec<_>>();
        kind_filter_upper.sort();
        kind_filter_upper.dedup();
        let feature_geometry_mode = feature_geometry_mode.unwrap_or_default();
        let requested_boundary_mode = feature_boundary_mode.unwrap_or_default();
        let feature_strand_relation = feature_strand_relation.unwrap_or_default();
        let effective_boundary_mode =
            if feature_geometry_mode == CandidateFeatureGeometryMode::FeatureBoundaries {
                requested_boundary_mode
            } else {
                CandidateFeatureBoundaryMode::Any
            };
        if feature_geometry_mode != CandidateFeatureGeometryMode::FeatureBoundaries
            && feature_boundary_mode.is_some()
        {
            result.warnings.push(
                "feature_boundary_mode is ignored unless feature_geometry_mode=feature_boundaries"
                    .to_string(),
            );
        }
        let feature_targets = Self::collect_feature_distance_targets(
            dna,
            feature_geometry_mode,
            effective_boundary_mode,
        );
        let matching_feature_count = Self::matching_feature_count(
            &feature_targets,
            &kind_filter_upper,
            label_regex.as_ref(),
            feature_strand_relation,
        );
        let has_feature_filter = !kind_filter_upper.is_empty()
            || label_regex.is_some()
            || max_distance_bp.is_some()
            || feature_strand_relation != CandidateFeatureStrandRelation::Any;
        if has_feature_filter && matching_feature_count == 0 {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "No features matched feature filters while generating candidate set"
                    .to_string(),
            });
        }

        let mut candidates = vec![];
        let mut considered = 0usize;
        let mut truncated = false;
        let upper = dna.len().saturating_sub(length_bp);
        let mut start = 0usize;
        while start <= upper {
            let end = start + length_bp;
            considered += 1;
            let distance_any = Self::nearest_feature_distance(
                start,
                end,
                &feature_targets,
                &[],
                None,
                CandidateFeatureStrandRelation::Any,
            );
            let distance_filtered = Self::nearest_feature_distance(
                start,
                end,
                &feature_targets,
                &kind_filter_upper,
                label_regex.as_ref(),
                feature_strand_relation,
            );
            let selected_distance = if !kind_filter_upper.is_empty()
                || label_regex.is_some()
                || feature_strand_relation != CandidateFeatureStrandRelation::Any
            {
                distance_filtered
            } else {
                distance_any
            };
            if let Some(max_distance) = max_distance_bp {
                let Some(distance) = selected_distance else {
                    start = start.saturating_add(step_bp);
                    continue;
                };
                if distance > max_distance {
                    start = start.saturating_add(step_bp);
                    continue;
                }
            }

            let Some(fragment) = dna.get_range_safe(start..end) else {
                start = start.saturating_add(step_bp);
                continue;
            };
            let sequence = String::from_utf8_lossy(&fragment).to_string();
            let mut metrics =
                Self::compute_candidate_metrics(&fragment, start, end, dna.len(), distance_any);
            if let Some(distance) = selected_distance {
                metrics.insert(
                    "distance_to_filtered_feature_bp".to_string(),
                    distance as f64,
                );
            }
            candidates.push(CandidateRecord {
                seq_id: seq_id.clone(),
                start_0based: start,
                end_0based: end,
                sequence,
                metrics,
            });
            if candidates.len() >= limit {
                truncated = true;
                break;
            }
            start = start.saturating_add(step_bp);
        }
        if candidates.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "No candidates matched generation constraints".to_string(),
            });
        }
        let generated = candidates.len();
        let mut store = self.read_candidate_store();
        let replaced_existing = store
            .sets
            .insert(
                set_name.clone(),
                CandidateSet {
                    name: set_name.clone(),
                    created_at_unix_ms: Self::now_unix_ms(),
                    source_seq_ids: vec![seq_id.clone()],
                    candidates,
                },
            )
            .is_some();
        let metric_names = store
            .sets
            .get(&set_name)
            .map(Self::metric_names_for_candidate_set)
            .unwrap_or_default();
        self.write_candidate_store(store)?;
        result.messages.push(format!(
            "Generated candidate set '{}' from '{}' ({} candidates, {} windows considered)",
            set_name, seq_id, generated, considered
        ));
        if truncated {
            result.warnings.push(format!(
                "Candidate generation for '{}' was truncated at limit={}",
                set_name, limit
            ));
        }
        if replaced_existing {
            result.warnings.push(format!(
                "Candidate set '{}' replaced existing set",
                set_name
            ));
        }
        if !metric_names.is_empty() {
            result.messages.push(format!(
                "Candidate set '{}' metrics: {}",
                set_name,
                metric_names.join(", ")
            ));
        }
        if matching_feature_count > 0 {
            result.messages.push(format!(
                "Candidate set '{}' matching feature count: {}",
                set_name, matching_feature_count
            ));
        }
        if !feature_targets.is_empty() {
            result.messages.push(format!(
                "Candidate set '{}' feature distance mode: geometry='{}', boundary='{}', strand_relation='{}'",
                set_name,
                feature_geometry_mode.as_str(),
                effective_boundary_mode.as_str(),
                feature_strand_relation.as_str()
            ));
        }
        Ok(())
    }

    pub(super) fn op_generate_candidate_set_between_anchors(
        &mut self,
        set_name: String,
        seq_id: SeqId,
        anchor_a: SequenceAnchor,
        anchor_b: SequenceAnchor,
        length_bp: usize,
        step_bp: usize,
        limit: Option<usize>,
        result: &mut OpResult,
    ) -> Result<(), EngineError> {
        let set_name = Self::normalize_candidate_set_name(&set_name)?;
        if length_bp == 0 {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "GenerateCandidateSetBetweenAnchors requires length_bp >= 1".to_string(),
            });
        }
        if step_bp == 0 {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "GenerateCandidateSetBetweenAnchors requires step_bp >= 1".to_string(),
            });
        }
        let limit = limit.unwrap_or(self.max_fragments_per_container());
        if limit == 0 {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "GenerateCandidateSetBetweenAnchors requires limit >= 1".to_string(),
            });
        }

        let dna = self
            .state
            .sequences
            .get(&seq_id)
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!("Sequence '{}' not found", seq_id),
            })?;
        let is_circular = dna.is_circular();

        let anchor_a_pos = Self::resolve_sequence_anchor_position(dna, &anchor_a, "anchor_a")?;
        let anchor_b_pos = Self::resolve_sequence_anchor_position(dna, &anchor_b, "anchor_b")?;
        if anchor_a_pos == anchor_b_pos {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "GenerateCandidateSetBetweenAnchors requires distinct anchor positions"
                    .to_string(),
            });
        }

        let left = anchor_a_pos.min(anchor_b_pos);
        let right = anchor_a_pos.max(anchor_b_pos);
        let span = right.saturating_sub(left);
        if span < length_bp {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Anchor span {} bp is shorter than requested candidate length {} bp",
                    span, length_bp
                ),
            });
        }

        let mut candidates = vec![];
        let mut considered = 0usize;
        let mut truncated = false;
        let upper = right.saturating_sub(length_bp);
        let mut start = left;
        while start <= upper {
            let end = start + length_bp;
            considered += 1;
            let Some(fragment) = dna.get_range_safe(start..end) else {
                start = start.saturating_add(step_bp);
                continue;
            };
            let sequence = String::from_utf8_lossy(&fragment).to_string();
            let mut metrics =
                Self::compute_candidate_metrics(&fragment, start, end, dna.len(), None);
            let dist_a = Self::boundary_distance(start, end, anchor_a_pos);
            let dist_b = Self::boundary_distance(start, end, anchor_b_pos);
            metrics.insert("anchor_a_position_bp".to_string(), anchor_a_pos as f64);
            metrics.insert("anchor_b_position_bp".to_string(), anchor_b_pos as f64);
            metrics.insert("distance_to_anchor_a_bp".to_string(), dist_a as f64);
            metrics.insert("distance_to_anchor_b_bp".to_string(), dist_b as f64);
            metrics.insert(
                "distance_to_nearest_anchor_bp".to_string(),
                dist_a.min(dist_b) as f64,
            );
            metrics.insert("anchor_interval_start_bp".to_string(), left as f64);
            metrics.insert("anchor_interval_end_bp".to_string(), right as f64);
            metrics.insert("anchor_interval_span_bp".to_string(), span as f64);
            metrics.insert(
                "anchor_order_sign".to_string(),
                if anchor_a_pos <= anchor_b_pos {
                    1.0
                } else {
                    -1.0
                },
            );
            metrics.insert(
                "candidate_offset_from_anchor_a_bp".to_string(),
                if anchor_a_pos <= anchor_b_pos {
                    start.saturating_sub(anchor_a_pos) as f64
                } else {
                    anchor_a_pos.saturating_sub(end) as f64
                },
            );
            candidates.push(CandidateRecord {
                seq_id: seq_id.clone(),
                start_0based: start,
                end_0based: end,
                sequence,
                metrics,
            });
            if candidates.len() >= limit {
                truncated = true;
                break;
            }
            start = start.saturating_add(step_bp);
        }

        if candidates.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "No candidates generated between the two anchors".to_string(),
            });
        }

        let generated = candidates.len();
        let mut store = self.read_candidate_store();
        let replaced_existing = store
            .sets
            .insert(
                set_name.clone(),
                CandidateSet {
                    name: set_name.clone(),
                    created_at_unix_ms: Self::now_unix_ms(),
                    source_seq_ids: vec![seq_id.clone()],
                    candidates,
                },
            )
            .is_some();
        let metric_names = store
            .sets
            .get(&set_name)
            .map(Self::metric_names_for_candidate_set)
            .unwrap_or_default();
        self.write_candidate_store(store)?;

        result.messages.push(format!(
            "Generated candidate set '{}' between anchors on '{}' ({} candidates, {} windows considered, span={} bp)",
            set_name, seq_id, generated, considered, span
        ));
        result.messages.push(format!(
            "Anchor positions on '{}': anchor_a={}, anchor_b={}",
            seq_id, anchor_a_pos, anchor_b_pos
        ));
        if truncated {
            result.warnings.push(format!(
                "Candidate generation for '{}' was truncated at limit={}",
                set_name, limit
            ));
        }
        if replaced_existing {
            result.warnings.push(format!(
                "Candidate set '{}' replaced existing set",
                set_name
            ));
        }
        if !metric_names.is_empty() {
            result.messages.push(format!(
                "Candidate set '{}' metrics: {}",
                set_name,
                metric_names.join(", ")
            ));
        }
        if is_circular {
            result.warnings.push(
                "GenerateCandidateSetBetweenAnchors currently uses linearized anchor interval semantics on circular sequences"
                    .to_string(),
            );
        }
        Ok(())
    }

    pub(super) fn append_guide_design_audit(
        store: &mut GuideDesignStore,
        operation: &str,
        guide_set_id: &str,
        payload: serde_json::Value,
    ) {
        store.audit_log.push(GuideDesignAuditEntry {
            unix_ms: Self::now_unix_ms(),
            operation: operation.to_string(),
            guide_set_id: guide_set_id.to_string(),
            payload,
        });
        if store.audit_log.len() > 512 {
            let drain = store.audit_log.len() - 512;
            store.audit_log.drain(0..drain);
        }
    }

    pub(super) fn op_upsert_guide_set(
        &mut self,
        guide_set_id: String,
        guides: Vec<GuideCandidate>,
        result: &mut OpResult,
    ) -> Result<(), EngineError> {
        let guide_set_id = Self::normalize_guide_set_id(&guide_set_id)?;
        let guides = Self::normalize_guide_candidates(guides)?;
        let mut store = self.read_guide_design_store();
        let now = Self::now_unix_ms();
        let created_at_unix_ms = store
            .guide_sets
            .get(&guide_set_id)
            .map(|set| set.created_at_unix_ms)
            .unwrap_or(now);
        let replaced_existing = store
            .guide_sets
            .insert(
                guide_set_id.clone(),
                GuideSet {
                    guide_set_id: guide_set_id.clone(),
                    guides: guides.clone(),
                    created_at_unix_ms,
                    updated_at_unix_ms: now,
                },
            )
            .is_some();
        Self::append_guide_design_audit(
            &mut store,
            "UpsertGuideSet",
            &guide_set_id,
            json!({
                "guide_count": guides.len(),
                "replaced_existing": replaced_existing
            }),
        );
        self.write_guide_design_store(store)?;
        result.messages.push(format!(
            "Upserted guide set '{}' with {} guide(s)",
            guide_set_id,
            guides.len()
        ));
        if replaced_existing {
            result.warnings.push(format!(
                "Guide set '{}' replaced existing content",
                guide_set_id
            ));
        }
        Ok(())
    }

    pub(super) fn op_delete_guide_set(
        &mut self,
        guide_set_id: String,
        result: &mut OpResult,
    ) -> Result<(), EngineError> {
        let guide_set_id = Self::normalize_guide_set_id(&guide_set_id)?;
        let mut store = self.read_guide_design_store();
        let removed = store.guide_sets.remove(&guide_set_id).is_some();
        let _ = store.practical_filter_reports.remove(&guide_set_id);
        let oligo_ids = store
            .oligo_sets
            .iter()
            .filter(|(_, set)| set.guide_set_id == guide_set_id)
            .map(|(id, _)| id.clone())
            .collect::<Vec<_>>();
        for oligo_id in &oligo_ids {
            let _ = store.oligo_sets.remove(oligo_id);
        }
        let _ = store.latest_oligo_set_by_guide_set.remove(&guide_set_id);
        if removed {
            Self::append_guide_design_audit(
                &mut store,
                "DeleteGuideSet",
                &guide_set_id,
                json!({
                    "removed_oligo_sets": oligo_ids.len()
                }),
            );
        }
        self.write_guide_design_store(store)?;
        if removed {
            result.messages.push(format!(
                "Deleted guide set '{}' and {} associated oligo set(s)",
                guide_set_id,
                oligo_ids.len()
            ));
        } else {
            result
                .warnings
                .push(format!("Guide set '{}' was not present", guide_set_id));
        }
        Ok(())
    }

    pub(super) fn op_filter_guides_practical(
        &mut self,
        guide_set_id: String,
        config: GuidePracticalFilterConfig,
        output_guide_set_id: Option<String>,
        result: &mut OpResult,
    ) -> Result<(), EngineError> {
        let guide_set_id = Self::normalize_guide_set_id(&guide_set_id)?;
        let output_guide_set_id = output_guide_set_id
            .as_deref()
            .map(Self::normalize_guide_set_id)
            .transpose()?;
        let config = Self::normalize_practical_filter_config(config)?;
        let mut store = self.read_guide_design_store();
        let guide_set = store
            .guide_sets
            .get(&guide_set_id)
            .cloned()
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!("Guide set '{}' not found", guide_set_id),
            })?;
        if guide_set.guides.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!("Guide set '{}' is empty", guide_set_id),
            });
        }

        let mut passed_guides: Vec<GuideCandidate> = Vec::new();
        let mut rows: Vec<GuidePracticalFilterResult> = Vec::with_capacity(guide_set.guides.len());

        for guide in &guide_set.guides {
            let spacer = guide
                .protospacer
                .as_bytes()
                .iter()
                .map(|b| match b.to_ascii_uppercase() {
                    b'U' => b'T',
                    other => other,
                })
                .collect::<Vec<_>>();
            let tail = guide
                .pam
                .as_bytes()
                .iter()
                .map(|b| match b.to_ascii_uppercase() {
                    b'U' => b'T',
                    other => other,
                })
                .collect::<Vec<_>>();

            let mut row = GuidePracticalFilterResult {
                guide_id: guide.guide_id.clone(),
                passed: true,
                reasons: vec![],
                warnings: vec![],
                metrics: HashMap::new(),
            };

            if config.reject_ambiguous_bases && Self::has_ambiguous_bases(&spacer) {
                row.reasons.push(GuideFilterReason {
                    code: "contains_ambiguous_base".to_string(),
                    message: "Protospacer contains non-ACGT bases".to_string(),
                });
            }

            if config.gc_min.is_some() || config.gc_max.is_some() {
                if let Some(gc) = Self::sequence_gc_fraction(&spacer) {
                    row.metrics.insert("gc_fraction".to_string(), gc);
                    if let Some(min) = config.gc_min {
                        if gc < min {
                            row.reasons.push(GuideFilterReason {
                                code: "gc_too_low".to_string(),
                                message: format!(
                                    "GC fraction {:.3} is below minimum {:.3}",
                                    gc, min
                                ),
                            });
                        }
                    }
                    if let Some(max) = config.gc_max {
                        if gc > max {
                            row.reasons.push(GuideFilterReason {
                                code: "gc_too_high".to_string(),
                                message: format!(
                                    "GC fraction {:.3} is above maximum {:.3}",
                                    gc, max
                                ),
                            });
                        }
                    }
                }
            }

            let max_run = Self::max_homopolymer_run(&spacer);
            row.metrics
                .insert("max_homopolymer_run".to_string(), max_run as f64);
            if let Some(limit) = config.max_homopolymer_run {
                if max_run > limit {
                    row.reasons.push(GuideFilterReason {
                        code: "homopolymer_run_exceeded".to_string(),
                        message: format!(
                            "Max homopolymer run {} exceeds configured limit {}",
                            max_run, limit
                        ),
                    });
                }
            }
            for (base, limit) in &config.max_homopolymer_run_per_base {
                let probe = base.as_bytes()[0];
                let observed = Self::max_homopolymer_run_for_base(&spacer, probe);
                row.metrics.insert(
                    format!("max_homopolymer_run_{}", base.to_ascii_lowercase()),
                    observed as f64,
                );
                if observed > *limit {
                    row.reasons.push(GuideFilterReason {
                        code: "homopolymer_run_exceeded".to_string(),
                        message: format!(
                            "Homopolymer run for base '{}' ({}) exceeds configured limit {}",
                            base, observed, limit
                        ),
                    });
                }
            }

            if config.avoid_u6_terminator_tttt {
                let mut window = spacer.clone();
                if config.u6_terminator_window == GuideU6TerminatorWindow::SpacerPlusTail {
                    window.extend_from_slice(&tail);
                }
                if Self::contains_u6_terminator_t4(&window) {
                    row.reasons.push(GuideFilterReason {
                        code: "u6_terminator_t4".to_string(),
                        message: "Contains TTTT in configured U6 terminator scan window"
                            .to_string(),
                    });
                }
            }

            let max_repeat = Self::max_dinucleotide_repeat_units(&spacer);
            row.metrics.insert(
                "max_dinucleotide_repeat_units".to_string(),
                max_repeat as f64,
            );
            if let Some(limit) = config.max_dinucleotide_repeat_units {
                if max_repeat > limit {
                    row.reasons.push(GuideFilterReason {
                        code: "dinucleotide_repeat_exceeded".to_string(),
                        message: format!(
                            "Max dinucleotide repeat units {} exceeds configured limit {}",
                            max_repeat, limit
                        ),
                    });
                }
            }

            for motif in &config.forbidden_motifs {
                if Self::contains_motif_any_strand(&spacer, motif)? {
                    row.reasons.push(GuideFilterReason {
                        code: "forbidden_motif_present".to_string(),
                        message: format!("Protospacer contains forbidden motif '{}'", motif),
                    });
                }
            }

            if let Some(required_base) = config.required_5prime_base.as_ref() {
                if let Some(actual) = spacer.first().map(|b| (*b as char).to_string()) {
                    if actual != *required_base {
                        if config.allow_5prime_g_extension && required_base == "G" {
                            row.warnings.push(
                                "Missing required 5' G but can be rescued by 5' G extension"
                                    .to_string(),
                            );
                        } else {
                            row.reasons.push(GuideFilterReason {
                                code: "required_5prime_base_missing".to_string(),
                                message: format!(
                                    "Protospacer 5' base '{}' does not match required '{}'",
                                    actual, required_base
                                ),
                            });
                        }
                    }
                }
            }

            row.passed = row.reasons.is_empty();
            if row.passed {
                passed_guides.push(guide.clone());
            }
            rows.push(row);
        }

        let report = GuidePracticalFilterReport {
            guide_set_id: guide_set_id.clone(),
            generated_at_unix_ms: Self::now_unix_ms(),
            config: config.clone(),
            passed_count: passed_guides.len(),
            rejected_count: rows.len().saturating_sub(passed_guides.len()),
            results: rows,
        };
        store
            .practical_filter_reports
            .insert(guide_set_id.clone(), report.clone());

        if let Some(output_guide_set_id) = output_guide_set_id {
            let now = Self::now_unix_ms();
            let created_at_unix_ms = store
                .guide_sets
                .get(&output_guide_set_id)
                .map(|set| set.created_at_unix_ms)
                .unwrap_or(now);
            let replaced = store
                .guide_sets
                .insert(
                    output_guide_set_id.clone(),
                    GuideSet {
                        guide_set_id: output_guide_set_id.clone(),
                        guides: passed_guides.clone(),
                        created_at_unix_ms,
                        updated_at_unix_ms: now,
                    },
                )
                .is_some();
            result.messages.push(format!(
                "Wrote passed guides to set '{}' ({} guide(s))",
                output_guide_set_id,
                passed_guides.len()
            ));
            if replaced {
                result.warnings.push(format!(
                    "Guide set '{}' replaced existing content",
                    output_guide_set_id
                ));
            }
        }

        Self::append_guide_design_audit(
            &mut store,
            "FilterGuidesPractical",
            &guide_set_id,
            json!({
                "passed_count": report.passed_count,
                "rejected_count": report.rejected_count
            }),
        );
        self.write_guide_design_store(store)?;
        result.messages.push(format!(
            "Filtered guide set '{}' (passed {}, rejected {})",
            guide_set_id, report.passed_count, report.rejected_count
        ));
        Ok(())
    }

    pub(super) fn op_generate_guide_oligos(
        &mut self,
        guide_set_id: String,
        template_id: String,
        apply_5prime_g_extension: Option<bool>,
        output_oligo_set_id: Option<String>,
        passed_only: Option<bool>,
        result: &mut OpResult,
    ) -> Result<(), EngineError> {
        let guide_set_id = Self::normalize_guide_set_id(&guide_set_id)?;
        let template = Self::built_in_guide_oligo_template(&template_id).ok_or_else(|| {
            EngineError {
                code: ErrorCode::NotFound,
                message: format!(
                    "Unknown guide oligo template '{}'; supported templates: lenti_bsmbi_u6_default, plain_forward_reverse",
                    template_id
                ),
            }
        })?;
        let apply_5prime_g_extension = apply_5prime_g_extension.unwrap_or(false);
        let passed_only = passed_only.unwrap_or(false);

        let mut store = self.read_guide_design_store();
        let guide_set = store
            .guide_sets
            .get(&guide_set_id)
            .cloned()
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!("Guide set '{}' not found", guide_set_id),
            })?;
        if guide_set.guides.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!("Guide set '{}' is empty", guide_set_id),
            });
        }

        let passed_lookup = if passed_only {
            let report = store
                .practical_filter_reports
                .get(&guide_set_id)
                .ok_or_else(|| EngineError {
                    code: ErrorCode::NotFound,
                    message: format!(
                        "No practical filter report found for '{}' (required because passed_only=true)",
                        guide_set_id
                    ),
                })?;
            Some(
                report
                    .results
                    .iter()
                    .filter(|r| r.passed)
                    .map(|r| r.guide_id.clone())
                    .collect::<HashSet<_>>(),
            )
        } else {
            None
        };

        let mut guides = guide_set.guides.clone();
        guides.sort_by(|a, b| {
            a.rank
                .unwrap_or(usize::MAX)
                .cmp(&b.rank.unwrap_or(usize::MAX))
                .then(a.guide_id.cmp(&b.guide_id))
        });
        let mut records = vec![];
        for guide in &guides {
            if let Some(lookup) = &passed_lookup {
                if !lookup.contains(&guide.guide_id) {
                    continue;
                }
            }
            let mut notes = vec![];
            let mut spacer = guide
                .protospacer
                .as_bytes()
                .iter()
                .map(|b| match b.to_ascii_uppercase() {
                    b'U' => b'T',
                    other => other,
                } as char)
                .collect::<String>();
            if apply_5prime_g_extension && !spacer.starts_with('G') {
                spacer = format!("G{spacer}");
                notes.push("5' G extension applied".to_string());
            }
            let mut forward = format!(
                "{}{}{}",
                template.forward_prefix, spacer, template.forward_suffix
            );
            let reverse_spacer = if template.reverse_uses_reverse_complement_of_spacer {
                Self::reverse_complement(&spacer)
            } else {
                spacer
            };
            let mut reverse = format!(
                "{}{}{}",
                template.reverse_prefix, reverse_spacer, template.reverse_suffix
            );
            if template.uppercase_output {
                forward = forward.to_ascii_uppercase();
                reverse = reverse.to_ascii_uppercase();
            }
            records.push(GuideOligoRecord {
                guide_id: guide.guide_id.clone(),
                rank: guide.rank,
                forward_oligo: forward,
                reverse_oligo: reverse,
                notes,
            });
        }
        if records.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "No guides selected for oligo generation in set '{}'",
                    guide_set_id
                ),
            });
        }

        let now = Self::now_unix_ms();
        let requested = output_oligo_set_id
            .as_deref()
            .map(Self::normalize_oligo_set_id)
            .transpose()?;
        let default_id = format!("{}_{}_{}", guide_set_id, template.template_id, now);
        let oligo_set_id =
            Self::unique_oligo_set_id(&store, requested.as_deref().unwrap_or(default_id.as_str()));
        store
            .latest_oligo_set_by_guide_set
            .insert(guide_set_id.clone(), oligo_set_id.clone());
        store.oligo_sets.insert(
            oligo_set_id.clone(),
            GuideOligoSet {
                oligo_set_id: oligo_set_id.clone(),
                guide_set_id: guide_set_id.clone(),
                generated_at_unix_ms: now,
                template: template.clone(),
                apply_5prime_g_extension,
                records: records.clone(),
            },
        );
        Self::append_guide_design_audit(
            &mut store,
            "GenerateGuideOligos",
            &guide_set_id,
            json!({
                "oligo_set_id": oligo_set_id,
                "record_count": records.len(),
                "template_id": template.template_id,
                "passed_only": passed_only,
                "apply_5prime_g_extension": apply_5prime_g_extension
            }),
        );
        self.write_guide_design_store(store)?;
        result.messages.push(format!(
            "Generated oligos for guide set '{}' ({} guide(s), template '{}')",
            guide_set_id,
            records.len(),
            template.template_id
        ));
        if passed_only {
            result
                .messages
                .push("Oligo generation used only practical-filter passing guides".to_string());
        }
        Ok(())
    }

    pub(super) fn csv_escape(value: &str) -> String {
        if value.contains(',') || value.contains('"') || value.contains('\n') {
            format!("\"{}\"", value.replace('"', "\"\""))
        } else {
            value.to_string()
        }
    }

    pub(super) fn ensure_output_parent_dir(path: &str) -> Result<(), EngineError> {
        let parent = Path::new(path)
            .parent()
            .map(Path::to_path_buf)
            .unwrap_or_else(|| PathBuf::from("."));
        std::fs::create_dir_all(&parent).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!(
                "Could not create output directory '{}' for export '{}': {e}",
                parent.display(),
                path
            ),
        })
    }

    pub(super) fn op_export_guide_oligos(
        &mut self,
        guide_set_id: String,
        oligo_set_id: Option<String>,
        format: GuideOligoExportFormat,
        path: String,
        plate_format: Option<GuideOligoPlateFormat>,
        result: &mut OpResult,
    ) -> Result<(), EngineError> {
        let guide_set_id = Self::normalize_guide_set_id(&guide_set_id)?;
        if path.trim().is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "ExportGuideOligos requires non-empty path".to_string(),
            });
        }
        let mut store = self.read_guide_design_store();
        if !store.guide_sets.contains_key(&guide_set_id) {
            return Err(EngineError {
                code: ErrorCode::NotFound,
                message: format!("Guide set '{}' not found", guide_set_id),
            });
        }
        let oligo_set =
            Self::resolve_oligo_set_for_export(&store, &guide_set_id, oligo_set_id.as_deref())?;

        let text = match format {
            GuideOligoExportFormat::CsvTable => {
                let mut rows = vec!["guide_id,rank,forward_oligo,reverse_oligo,notes".to_string()];
                for record in &oligo_set.records {
                    let rank = record.rank.map(|v| v.to_string()).unwrap_or_default();
                    let notes = record.notes.join("; ");
                    rows.push(format!(
                        "{},{},{},{},{}",
                        Self::csv_escape(&record.guide_id),
                        Self::csv_escape(&rank),
                        Self::csv_escape(&record.forward_oligo),
                        Self::csv_escape(&record.reverse_oligo),
                        Self::csv_escape(&notes),
                    ));
                }
                rows.join("\n")
            }
            GuideOligoExportFormat::PlateCsv => {
                let plate_format = plate_format.unwrap_or_default();
                let (rows_per_plate, cols_per_plate) = plate_format.dimensions();
                let capacity = rows_per_plate * cols_per_plate;
                let mut rows =
                    vec!["plate,well,guide_id,rank,forward_oligo,reverse_oligo,notes".to_string()];
                for (idx, record) in oligo_set.records.iter().enumerate() {
                    let plate_index = idx / capacity + 1;
                    let within_plate = idx % capacity;
                    let row = within_plate / cols_per_plate;
                    let col = within_plate % cols_per_plate + 1;
                    let row_char = (b'A' + row as u8) as char;
                    let well = format!("{row_char}{:02}", col);
                    let rank = record.rank.map(|v| v.to_string()).unwrap_or_default();
                    let notes = record.notes.join("; ");
                    rows.push(format!(
                        "{},{},{},{},{},{},{}",
                        plate_index,
                        Self::csv_escape(&well),
                        Self::csv_escape(&record.guide_id),
                        Self::csv_escape(&rank),
                        Self::csv_escape(&record.forward_oligo),
                        Self::csv_escape(&record.reverse_oligo),
                        Self::csv_escape(&notes),
                    ));
                }
                rows.join("\n")
            }
            GuideOligoExportFormat::Fasta => {
                let mut out = String::new();
                for record in &oligo_set.records {
                    out.push_str(&format!(
                        ">{}|forward\n{}\n",
                        record.guide_id, record.forward_oligo
                    ));
                    out.push_str(&format!(
                        ">{}|reverse\n{}\n",
                        record.guide_id, record.reverse_oligo
                    ));
                }
                out
            }
        };

        Self::ensure_output_parent_dir(&path)?;
        std::fs::write(&path, text).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!("Could not write oligo export '{}': {e}", path),
        })?;

        let report = GuideOligoExportReport {
            guide_set_id: guide_set_id.clone(),
            oligo_set_id: oligo_set.oligo_set_id.clone(),
            format: format.as_str().to_string(),
            path: path.clone(),
            exported_records: oligo_set.records.len(),
        };
        Self::append_guide_design_audit(
            &mut store,
            "ExportGuideOligos",
            &guide_set_id,
            serde_json::to_value(&report).unwrap_or_else(|_| json!({ "path": path })),
        );
        self.write_guide_design_store(store)?;
        result.messages.push(format!(
            "Exported {} oligo records for '{}' to '{}' as {}",
            report.exported_records, guide_set_id, report.path, report.format
        ));
        Ok(())
    }

    pub(super) fn op_export_guide_protocol_text(
        &mut self,
        guide_set_id: String,
        oligo_set_id: Option<String>,
        path: String,
        include_qc_checklist: Option<bool>,
        result: &mut OpResult,
    ) -> Result<(), EngineError> {
        let guide_set_id = Self::normalize_guide_set_id(&guide_set_id)?;
        if path.trim().is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "ExportGuideProtocolText requires non-empty path".to_string(),
            });
        }
        let mut store = self.read_guide_design_store();
        if !store.guide_sets.contains_key(&guide_set_id) {
            return Err(EngineError {
                code: ErrorCode::NotFound,
                message: format!("Guide set '{}' not found", guide_set_id),
            });
        }
        let oligo_set =
            Self::resolve_oligo_set_for_export(&store, &guide_set_id, oligo_set_id.as_deref())?;
        let include_qc = include_qc_checklist.unwrap_or(true);

        let mut text = String::new();
        text.push_str("GENtle Guide Oligo Protocol\n");
        text.push_str("==========================\n\n");
        text.push_str(&format!("Guide set: {}\n", guide_set_id));
        text.push_str(&format!("Oligo set: {}\n", oligo_set.oligo_set_id));
        text.push_str(&format!(
            "Template: {} ({})\n",
            oligo_set.template.template_id, oligo_set.template.description
        ));
        text.push_str(&format!(
            "Generated guides: {}\n\n",
            oligo_set.records.len()
        ));
        text.push_str("Suggested steps:\n");
        text.push_str("1. Prepare oligo stocks according to ordering sheet.\n");
        text.push_str("2. Anneal forward/reverse oligo pairs.\n");
        text.push_str("3. Clone into vector backbone using template-defined overhang strategy.\n");
        text.push_str("4. Transform, recover colonies, and verify insertion by sequencing.\n\n");
        text.push_str("Guide oligos:\n");
        for record in &oligo_set.records {
            let rank = record
                .rank
                .map(|v| format!("rank={v}"))
                .unwrap_or_else(|| "rank=-".to_string());
            text.push_str(&format!(
                "- {} ({})\n  forward: {}\n  reverse: {}\n",
                record.guide_id, rank, record.forward_oligo, record.reverse_oligo
            ));
            if !record.notes.is_empty() {
                text.push_str(&format!("  notes: {}\n", record.notes.join("; ")));
            }
        }
        if include_qc {
            text.push_str("\nQC checklist:\n");
            text.push_str("- Confirm oligo lengths and overhang sequences.\n");
            text.push_str(
                "- Confirm no guide contains forbidden motifs for your cloning strategy.\n",
            );
            text.push_str("- Confirm expected insert size by colony PCR or digest.\n");
            text.push_str("- Confirm sequence identity by Sanger/NGS.\n");
        }

        Self::ensure_output_parent_dir(&path)?;
        std::fs::write(&path, text).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!("Could not write protocol export '{}': {e}", path),
        })?;

        let report = GuideProtocolExportReport {
            guide_set_id: guide_set_id.clone(),
            oligo_set_id: oligo_set.oligo_set_id.clone(),
            path: path.clone(),
            guide_count: oligo_set.records.len(),
        };
        Self::append_guide_design_audit(
            &mut store,
            "ExportGuideProtocolText",
            &guide_set_id,
            serde_json::to_value(&report).unwrap_or_else(|_| json!({ "path": path })),
        );
        self.write_guide_design_store(store)?;
        result.messages.push(format!(
            "Exported guide protocol text for '{}' to '{}'",
            guide_set_id, report.path
        ));
        Ok(())
    }

    pub(super) fn op_delete_candidate_set(
        &mut self,
        set_name: String,
        result: &mut OpResult,
    ) -> Result<(), EngineError> {
        let set_name = Self::normalize_candidate_set_name(&set_name)?;
        let mut store = self.read_candidate_store();
        let removed = store.sets.remove(&set_name).is_some();
        self.write_candidate_store(store)?;
        if removed {
            result
                .messages
                .push(format!("Deleted candidate set '{}'", set_name));
        } else {
            result
                .warnings
                .push(format!("Candidate set '{}' was not present", set_name));
        }
        Ok(())
    }

    pub(super) fn op_score_candidate_set_expression(
        &mut self,
        set_name: String,
        metric: String,
        expression: String,
        result: &mut OpResult,
    ) -> Result<(), EngineError> {
        let set_name = Self::normalize_candidate_set_name(&set_name)?;
        let metric_name = Self::normalize_metric_name(&metric);
        let expr = Self::parse_metric_expression(&expression)?;
        let mut store = self.read_candidate_store();
        let set = store.sets.get_mut(&set_name).ok_or_else(|| EngineError {
            code: ErrorCode::NotFound,
            message: format!("Candidate set '{}' not found", set_name),
        })?;
        if set.candidates.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!("Candidate set '{}' is empty", set_name),
            });
        }
        let mut values = Vec::with_capacity(set.candidates.len());
        for (idx, candidate) in set.candidates.iter().enumerate() {
            let value =
                Self::evaluate_metric_expression(&expr, &candidate.metrics).map_err(|e| {
                    EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "Could not evaluate expression for candidate {} in '{}': {}",
                            idx, set_name, e.message
                        ),
                    }
                })?;
            values.push(value);
        }
        for (candidate, value) in set.candidates.iter_mut().zip(values.iter()) {
            candidate.metrics.insert(metric_name.clone(), *value);
        }
        let min_value = values.iter().copied().fold(f64::INFINITY, f64::min);
        let max_value = values.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        self.write_candidate_store(store)?;
        result.messages.push(format!(
            "Scored candidate set '{}' with metric '{}' from expression '{}'",
            set_name, metric_name, expression
        ));
        result.messages.push(format!(
            "Metric '{}' range in '{}': [{:.6}, {:.6}]",
            metric_name, set_name, min_value, max_value
        ));
        Ok(())
    }

    pub(super) fn op_score_candidate_set_distance(
        &mut self,
        set_name: String,
        metric: String,
        feature_kinds: Vec<String>,
        feature_label_regex: Option<String>,
        feature_geometry_mode: Option<CandidateFeatureGeometryMode>,
        feature_boundary_mode: Option<CandidateFeatureBoundaryMode>,
        feature_strand_relation: Option<CandidateFeatureStrandRelation>,
        result: &mut OpResult,
    ) -> Result<(), EngineError> {
        let set_name = Self::normalize_candidate_set_name(&set_name)?;
        let metric_name = Self::normalize_metric_name(&metric);
        let label_regex =
            Self::compile_optional_regex(&feature_label_regex, "feature_label_regex")?;
        let mut kind_filter_upper = feature_kinds
            .iter()
            .map(|kind| kind.trim().to_ascii_uppercase())
            .filter(|kind| !kind.is_empty())
            .collect::<Vec<_>>();
        kind_filter_upper.sort();
        kind_filter_upper.dedup();
        let feature_geometry_mode = feature_geometry_mode.unwrap_or_default();
        let requested_boundary_mode = feature_boundary_mode.unwrap_or_default();
        let feature_strand_relation = feature_strand_relation.unwrap_or_default();
        let effective_boundary_mode =
            if feature_geometry_mode == CandidateFeatureGeometryMode::FeatureBoundaries {
                requested_boundary_mode
            } else {
                CandidateFeatureBoundaryMode::Any
            };
        if feature_geometry_mode != CandidateFeatureGeometryMode::FeatureBoundaries
            && feature_boundary_mode.is_some()
        {
            result.warnings.push(
                "feature_boundary_mode is ignored unless feature_geometry_mode=feature_boundaries"
                    .to_string(),
            );
        }

        let mut store = self.read_candidate_store();
        let set = store.sets.get_mut(&set_name).ok_or_else(|| EngineError {
            code: ErrorCode::NotFound,
            message: format!("Candidate set '{}' not found", set_name),
        })?;
        if set.candidates.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!("Candidate set '{}' is empty", set_name),
            });
        }

        let mut feature_cache: HashMap<String, Vec<FeatureDistanceTarget>> = HashMap::new();
        for seq_id in set
            .candidates
            .iter()
            .map(|c| c.seq_id.clone())
            .collect::<BTreeSet<_>>()
        {
            let Some(dna) = self.state.sequences.get(&seq_id) else {
                return Err(EngineError {
                    code: ErrorCode::NotFound,
                    message: format!(
                        "Candidate set '{}' references missing sequence '{}'",
                        set_name, seq_id
                    ),
                });
            };
            feature_cache.insert(
                seq_id.clone(),
                Self::collect_feature_distance_targets(
                    dna,
                    feature_geometry_mode,
                    effective_boundary_mode,
                ),
            );
        }

        let mut values = Vec::with_capacity(set.candidates.len());
        for (idx, candidate) in set.candidates.iter().enumerate() {
            let features = feature_cache
                .get(&candidate.seq_id)
                .ok_or_else(|| EngineError {
                    code: ErrorCode::Internal,
                    message: format!(
                        "Missing feature cache for sequence '{}' while scoring candidate set '{}'",
                        candidate.seq_id, set_name
                    ),
                })?;
            let distance = Self::nearest_feature_distance(
                candidate.start_0based,
                candidate.end_0based,
                features,
                &kind_filter_upper,
                label_regex.as_ref(),
                feature_strand_relation,
            )
            .ok_or_else(|| EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "No matching features found for candidate {} (seq='{}')",
                    idx, candidate.seq_id
                ),
            })?;
            values.push(distance as f64);
        }

        for (candidate, value) in set.candidates.iter_mut().zip(values.iter()) {
            candidate.metrics.insert(metric_name.clone(), *value);
        }
        let min_value = values.iter().copied().fold(f64::INFINITY, f64::min);
        let max_value = values.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        self.write_candidate_store(store)?;
        result.messages.push(format!(
            "Scored candidate set '{}' with distance metric '{}'",
            set_name, metric_name
        ));
        result.messages.push(format!(
            "Metric '{}' range in '{}': [{:.6}, {:.6}]",
            metric_name, set_name, min_value, max_value
        ));
        result.messages.push(format!(
            "Distance scoring mode for '{}': geometry='{}', boundary='{}', strand_relation='{}'",
            set_name,
            feature_geometry_mode.as_str(),
            effective_boundary_mode.as_str(),
            feature_strand_relation.as_str()
        ));
        Ok(())
    }

    pub(super) fn op_filter_candidate_set(
        &mut self,
        input_set: String,
        output_set: String,
        metric: String,
        min: Option<f64>,
        max: Option<f64>,
        min_quantile: Option<f64>,
        max_quantile: Option<f64>,
        result: &mut OpResult,
    ) -> Result<(), EngineError> {
        let input_set = Self::normalize_candidate_set_name(&input_set)?;
        let output_set = Self::normalize_candidate_set_name(&output_set)?;
        let metric_name = Self::normalize_metric_name(&metric);
        if min.zip(max).map(|(a, b)| a > b).unwrap_or(false) {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "FilterCandidateSet requires min <= max".to_string(),
            });
        }
        if min_quantile
            .zip(max_quantile)
            .map(|(a, b)| a > b)
            .unwrap_or(false)
        {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "FilterCandidateSet requires min_quantile <= max_quantile".to_string(),
            });
        }
        for (name, value) in [
            ("min_quantile", min_quantile),
            ("max_quantile", max_quantile),
        ] {
            if let Some(q) = value {
                if !(0.0..=1.0).contains(&q) {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!("FilterCandidateSet {} must be between 0 and 1", name),
                    });
                }
            }
        }
        if min.is_none() && max.is_none() && min_quantile.is_none() && max_quantile.is_none() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "FilterCandidateSet requires at least one bound or quantile".to_string(),
            });
        }

        let mut store = self.read_candidate_store();
        let input = store
            .sets
            .get(&input_set)
            .cloned()
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!("Candidate set '{}' not found", input_set),
            })?;
        if input.candidates.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!("Candidate set '{}' is empty", input_set),
            });
        }

        let mut metric_values = Vec::with_capacity(input.candidates.len());
        for (idx, candidate) in input.candidates.iter().enumerate() {
            let value = candidate
                .metrics
                .get(&metric_name)
                .copied()
                .ok_or_else(|| EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!(
                        "Candidate {} in '{}' is missing metric '{}'",
                        idx, input_set, metric_name
                    ),
                })?;
            if !value.is_finite() {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!(
                        "Candidate {} in '{}' has non-finite metric '{}'",
                        idx, input_set, metric_name
                    ),
                });
            }
            metric_values.push(value);
        }

        let min_q_threshold =
            min_quantile.and_then(|q| Self::quantile_threshold(&metric_values, q));
        let max_q_threshold =
            max_quantile.and_then(|q| Self::quantile_threshold(&metric_values, q));
        let lower_bound = match (min, min_q_threshold) {
            (Some(a), Some(b)) => Some(a.max(b)),
            (Some(a), None) => Some(a),
            (None, Some(b)) => Some(b),
            (None, None) => None,
        };
        let upper_bound = match (max, max_q_threshold) {
            (Some(a), Some(b)) => Some(a.min(b)),
            (Some(a), None) => Some(a),
            (None, Some(b)) => Some(b),
            (None, None) => None,
        };

        let filtered_candidates = input
            .candidates
            .into_iter()
            .filter(|candidate| {
                let Some(value) = candidate.metrics.get(&metric_name).copied() else {
                    return false;
                };
                if let Some(lo) = lower_bound {
                    if value < lo {
                        return false;
                    }
                }
                if let Some(hi) = upper_bound {
                    if value > hi {
                        return false;
                    }
                }
                true
            })
            .collect::<Vec<_>>();
        let kept = filtered_candidates.len();
        let dropped = metric_values.len().saturating_sub(kept);
        let replaced_existing = store
            .sets
            .insert(
                output_set.clone(),
                CandidateSet {
                    name: output_set.clone(),
                    created_at_unix_ms: Self::now_unix_ms(),
                    source_seq_ids: input.source_seq_ids.clone(),
                    candidates: filtered_candidates,
                },
            )
            .is_some();
        self.write_candidate_store(store)?;
        result.messages.push(format!(
            "Filtered candidate set '{}' by metric '{}' into '{}' (kept {}, dropped {})",
            input_set, metric_name, output_set, kept, dropped
        ));
        if replaced_existing {
            result.warnings.push(format!(
                "FilterCandidateSet output '{}' replaced existing candidate set",
                output_set
            ));
        }
        Ok(())
    }

    pub(super) fn op_candidate_set_op(
        &mut self,
        op: CandidateSetOperator,
        left_set: String,
        right_set: String,
        output_set: String,
        result: &mut OpResult,
    ) -> Result<(), EngineError> {
        let left_set = Self::normalize_candidate_set_name(&left_set)?;
        let right_set = Self::normalize_candidate_set_name(&right_set)?;
        let output_set = Self::normalize_candidate_set_name(&output_set)?;
        let mut store = self.read_candidate_store();
        let left = store
            .sets
            .get(&left_set)
            .cloned()
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!("Candidate set '{}' not found", left_set),
            })?;
        let right = store
            .sets
            .get(&right_set)
            .cloned()
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!("Candidate set '{}' not found", right_set),
            })?;
        let mut output_candidates = vec![];
        match op {
            CandidateSetOperator::Union => {
                let mut seen = HashSet::new();
                for candidate in left.candidates.iter().chain(right.candidates.iter()) {
                    let key = Self::candidate_key(candidate);
                    if seen.insert(key) {
                        output_candidates.push(candidate.clone());
                    }
                }
            }
            CandidateSetOperator::Intersect => {
                let right_keys = right
                    .candidates
                    .iter()
                    .map(Self::candidate_key)
                    .collect::<HashSet<_>>();
                for candidate in &left.candidates {
                    if right_keys.contains(&Self::candidate_key(candidate)) {
                        output_candidates.push(candidate.clone());
                    }
                }
            }
            CandidateSetOperator::Subtract => {
                let right_keys = right
                    .candidates
                    .iter()
                    .map(Self::candidate_key)
                    .collect::<HashSet<_>>();
                for candidate in &left.candidates {
                    if !right_keys.contains(&Self::candidate_key(candidate)) {
                        output_candidates.push(candidate.clone());
                    }
                }
            }
        }
        let mut source_seq_ids = left.source_seq_ids.clone();
        source_seq_ids.extend(right.source_seq_ids.clone());
        source_seq_ids.sort();
        source_seq_ids.dedup();
        let output_count = output_candidates.len();
        let replaced_existing = store
            .sets
            .insert(
                output_set.clone(),
                CandidateSet {
                    name: output_set.clone(),
                    created_at_unix_ms: Self::now_unix_ms(),
                    source_seq_ids,
                    candidates: output_candidates,
                },
            )
            .is_some();
        self.write_candidate_store(store)?;
        result.messages.push(format!(
            "Candidate set {} '{}' and '{}' into '{}' ({} candidates)",
            op.as_str(),
            left_set,
            right_set,
            output_set,
            output_count
        ));
        if replaced_existing {
            result.warnings.push(format!(
                "CandidateSetOp output '{}' replaced existing candidate set",
                output_set
            ));
        }
        Ok(())
    }

    pub(super) fn compare_candidates_by_tie_break(
        a: &CandidateRecord,
        b: &CandidateRecord,
        policy: CandidateTieBreakPolicy,
    ) -> Ordering {
        let a_len = a.end_0based.saturating_sub(a.start_0based);
        let b_len = b.end_0based.saturating_sub(b.start_0based);
        match policy {
            CandidateTieBreakPolicy::SeqStartEnd => a
                .seq_id
                .cmp(&b.seq_id)
                .then(a.start_0based.cmp(&b.start_0based))
                .then(a.end_0based.cmp(&b.end_0based))
                .then(a.sequence.cmp(&b.sequence)),
            CandidateTieBreakPolicy::SeqEndStart => a
                .seq_id
                .cmp(&b.seq_id)
                .then(a.end_0based.cmp(&b.end_0based))
                .then(a.start_0based.cmp(&b.start_0based))
                .then(a.sequence.cmp(&b.sequence)),
            CandidateTieBreakPolicy::LengthAscending => a_len
                .cmp(&b_len)
                .then(a.seq_id.cmp(&b.seq_id))
                .then(a.start_0based.cmp(&b.start_0based))
                .then(a.end_0based.cmp(&b.end_0based))
                .then(a.sequence.cmp(&b.sequence)),
            CandidateTieBreakPolicy::LengthDescending => b_len
                .cmp(&a_len)
                .then(a.seq_id.cmp(&b.seq_id))
                .then(a.start_0based.cmp(&b.start_0based))
                .then(a.end_0based.cmp(&b.end_0based))
                .then(a.sequence.cmp(&b.sequence)),
            CandidateTieBreakPolicy::SequenceLexicographic => a
                .sequence
                .cmp(&b.sequence)
                .then(a.seq_id.cmp(&b.seq_id))
                .then(a.start_0based.cmp(&b.start_0based))
                .then(a.end_0based.cmp(&b.end_0based)),
        }
    }

    pub(super) fn op_score_candidate_set_weighted_objective(
        &mut self,
        set_name: String,
        metric: String,
        objectives: Vec<CandidateWeightedObjectiveTerm>,
        normalize_metrics: Option<bool>,
        result: &mut OpResult,
    ) -> Result<(), EngineError> {
        let set_name = Self::normalize_candidate_set_name(&set_name)?;
        let metric_name = Self::normalize_metric_name(&metric);
        if objectives.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "ScoreCandidateSetWeightedObjective requires at least one objective term"
                    .to_string(),
            });
        }
        let normalize_metrics = normalize_metrics.unwrap_or(true);

        let mut compiled = Vec::with_capacity(objectives.len());
        for term in objectives {
            let metric = Self::normalize_metric_name(&term.metric);
            if !term.weight.is_finite() || term.weight <= 0.0 {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!(
                        "Invalid weighted objective term for metric '{}': weight must be finite and > 0",
                        metric
                    ),
                });
            }
            compiled.push((metric, term.weight, term.direction));
        }

        let mut store = self.read_candidate_store();
        let set = store.sets.get_mut(&set_name).ok_or_else(|| EngineError {
            code: ErrorCode::NotFound,
            message: format!("Candidate set '{}' not found", set_name),
        })?;
        if set.candidates.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!("Candidate set '{}' is empty", set_name),
            });
        }

        let mut bounds = Vec::with_capacity(compiled.len());
        for (metric_name, _, _) in &compiled {
            let mut min_value = f64::INFINITY;
            let mut max_value = f64::NEG_INFINITY;
            for (idx, candidate) in set.candidates.iter().enumerate() {
                let value =
                    candidate
                        .metrics
                        .get(metric_name)
                        .copied()
                        .ok_or_else(|| EngineError {
                            code: ErrorCode::InvalidInput,
                            message: format!(
                                "Candidate {} in '{}' is missing metric '{}'",
                                idx, set_name, metric_name
                            ),
                        })?;
                if !value.is_finite() {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "Candidate {} in '{}' has non-finite metric '{}'",
                            idx, set_name, metric_name
                        ),
                    });
                }
                min_value = min_value.min(value);
                max_value = max_value.max(value);
            }
            bounds.push((min_value, max_value));
        }

        let mut combined_values = Vec::with_capacity(set.candidates.len());
        for candidate in &mut set.candidates {
            let mut combined = 0.0f64;
            for ((metric_name, weight, direction), (min_value, max_value)) in
                compiled.iter().zip(bounds.iter())
            {
                let raw_value =
                    candidate
                        .metrics
                        .get(metric_name)
                        .copied()
                        .ok_or_else(|| EngineError {
                            code: ErrorCode::InvalidInput,
                            message: format!(
                                "Candidate in '{}' is missing metric '{}'",
                                set_name, metric_name
                            ),
                        })?;
                let objective_value = if normalize_metrics {
                    let scaled = if *max_value > *min_value {
                        (raw_value - *min_value) / (*max_value - *min_value)
                    } else {
                        0.5
                    };
                    match direction {
                        CandidateObjectiveDirection::Maximize => scaled,
                        CandidateObjectiveDirection::Minimize => 1.0 - scaled,
                    }
                } else {
                    match direction {
                        CandidateObjectiveDirection::Maximize => raw_value,
                        CandidateObjectiveDirection::Minimize => -raw_value,
                    }
                };
                combined += *weight * objective_value;
            }
            candidate.metrics.insert(metric_name.clone(), combined);
            combined_values.push(combined);
        }

        let min_combined = combined_values
            .iter()
            .copied()
            .fold(f64::INFINITY, f64::min);
        let max_combined = combined_values
            .iter()
            .copied()
            .fold(f64::NEG_INFINITY, f64::max);
        self.write_candidate_store(store)?;
        result.messages.push(format!(
            "Scored candidate set '{}' with weighted objective metric '{}'",
            set_name, metric_name
        ));
        result.messages.push(format!(
            "Weighted objective mode for '{}': normalize_metrics={}",
            set_name, normalize_metrics
        ));
        result.messages.push(format!(
            "Metric '{}' range in '{}': [{:.6}, {:.6}]",
            metric_name, set_name, min_combined, max_combined
        ));
        Ok(())
    }

    pub(super) fn op_top_k_candidate_set(
        &mut self,
        input_set: String,
        output_set: String,
        metric: String,
        k: usize,
        direction: Option<CandidateObjectiveDirection>,
        tie_break: Option<CandidateTieBreakPolicy>,
        result: &mut OpResult,
    ) -> Result<(), EngineError> {
        let input_set = Self::normalize_candidate_set_name(&input_set)?;
        let output_set = Self::normalize_candidate_set_name(&output_set)?;
        let metric_name = Self::normalize_metric_name(&metric);
        if k == 0 {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "TopKCandidateSet requires k >= 1".to_string(),
            });
        }
        let direction = direction.unwrap_or_default();
        let tie_break = tie_break.unwrap_or_default();

        let mut store = self.read_candidate_store();
        let input = store
            .sets
            .get(&input_set)
            .cloned()
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!("Candidate set '{}' not found", input_set),
            })?;
        if input.candidates.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!("Candidate set '{}' is empty", input_set),
            });
        }

        let mut scored = Vec::with_capacity(input.candidates.len());
        for (idx, candidate) in input.candidates.into_iter().enumerate() {
            let value = candidate
                .metrics
                .get(&metric_name)
                .copied()
                .ok_or_else(|| EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!(
                        "Candidate {} in '{}' is missing metric '{}'",
                        idx, input_set, metric_name
                    ),
                })?;
            if !value.is_finite() {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!(
                        "Candidate {} in '{}' has non-finite metric '{}'",
                        idx, input_set, metric_name
                    ),
                });
            }
            scored.push((candidate, value));
        }
        scored.sort_by(|(left, left_value), (right, right_value)| {
            let primary = match direction {
                CandidateObjectiveDirection::Maximize => right_value
                    .partial_cmp(left_value)
                    .unwrap_or(Ordering::Equal),
                CandidateObjectiveDirection::Minimize => left_value
                    .partial_cmp(right_value)
                    .unwrap_or(Ordering::Equal),
            };
            if primary == Ordering::Equal {
                Self::compare_candidates_by_tie_break(left, right, tie_break)
            } else {
                primary
            }
        });

        let selected = scored
            .into_iter()
            .take(k)
            .map(|(candidate, _)| candidate)
            .collect::<Vec<_>>();
        let selected_count = selected.len();
        let replaced_existing = store
            .sets
            .insert(
                output_set.clone(),
                CandidateSet {
                    name: output_set.clone(),
                    created_at_unix_ms: Self::now_unix_ms(),
                    source_seq_ids: input.source_seq_ids,
                    candidates: selected,
                },
            )
            .is_some();
        self.write_candidate_store(store)?;
        result.messages.push(format!(
            "Selected top {} candidate(s) from '{}' into '{}' by metric '{}' ({}, tie_break={})",
            selected_count,
            input_set,
            output_set,
            metric_name,
            direction.as_str(),
            tie_break.as_str()
        ));
        if replaced_existing {
            result.warnings.push(format!(
                "TopKCandidateSet output '{}' replaced existing candidate set",
                output_set
            ));
        }
        Ok(())
    }

    pub(super) fn candidate_dominates(
        left_values: &[f64],
        right_values: &[f64],
        objectives: &[CandidateObjectiveSpec],
    ) -> bool {
        let mut strictly_better = false;
        for (idx, objective) in objectives.iter().enumerate() {
            let left = left_values.get(idx).copied().unwrap_or(0.0);
            let right = right_values.get(idx).copied().unwrap_or(0.0);
            match objective.direction {
                CandidateObjectiveDirection::Maximize => {
                    if left < right {
                        return false;
                    }
                    if left > right {
                        strictly_better = true;
                    }
                }
                CandidateObjectiveDirection::Minimize => {
                    if left > right {
                        return false;
                    }
                    if left < right {
                        strictly_better = true;
                    }
                }
            }
        }
        strictly_better
    }

    pub(super) fn op_pareto_frontier_candidate_set(
        &mut self,
        input_set: String,
        output_set: String,
        objectives: Vec<CandidateObjectiveSpec>,
        max_candidates: Option<usize>,
        tie_break: Option<CandidateTieBreakPolicy>,
        result: &mut OpResult,
    ) -> Result<(), EngineError> {
        let input_set = Self::normalize_candidate_set_name(&input_set)?;
        let output_set = Self::normalize_candidate_set_name(&output_set)?;
        if objectives.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "ParetoFrontierCandidateSet requires at least one objective".to_string(),
            });
        }
        if max_candidates == Some(0) {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "ParetoFrontierCandidateSet max_candidates must be >= 1".to_string(),
            });
        }
        let tie_break = tie_break.unwrap_or_default();

        let mut compiled = Vec::with_capacity(objectives.len());
        for objective in objectives {
            compiled.push(CandidateObjectiveSpec {
                metric: Self::normalize_metric_name(&objective.metric),
                direction: objective.direction,
            });
        }

        let mut store = self.read_candidate_store();
        let input = store
            .sets
            .get(&input_set)
            .cloned()
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!("Candidate set '{}' not found", input_set),
            })?;
        if input.candidates.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!("Candidate set '{}' is empty", input_set),
            });
        }

        let mut objective_values: Vec<Vec<f64>> = vec![];
        for (candidate_idx, candidate) in input.candidates.iter().enumerate() {
            let mut row = Vec::with_capacity(compiled.len());
            for objective in &compiled {
                let value = candidate
                    .metrics
                    .get(&objective.metric)
                    .copied()
                    .ok_or_else(|| EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "Candidate {} in '{}' is missing metric '{}'",
                            candidate_idx, input_set, objective.metric
                        ),
                    })?;
                if !value.is_finite() {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "Candidate {} in '{}' has non-finite metric '{}'",
                            candidate_idx, input_set, objective.metric
                        ),
                    });
                }
                row.push(value);
            }
            objective_values.push(row);
        }

        let mut dominated = vec![false; input.candidates.len()];
        for i in 0..input.candidates.len() {
            if dominated[i] {
                continue;
            }
            for j in 0..input.candidates.len() {
                if i == j {
                    continue;
                }
                if Self::candidate_dominates(&objective_values[j], &objective_values[i], &compiled)
                {
                    dominated[i] = true;
                    break;
                }
            }
        }

        let mut frontier = input
            .candidates
            .into_iter()
            .zip(dominated.into_iter())
            .filter_map(
                |(candidate, is_dominated)| if is_dominated { None } else { Some(candidate) },
            )
            .collect::<Vec<_>>();
        let raw_frontier_count = frontier.len();
        if let Some(limit) = max_candidates {
            if frontier.len() > limit {
                frontier.sort_by(|a, b| Self::compare_candidates_by_tie_break(a, b, tie_break));
                frontier.truncate(limit);
                result.warnings.push(format!(
                    "Pareto frontier for '{}' had {} candidates and was truncated to {} by tie-break policy '{}'",
                    input_set,
                    raw_frontier_count,
                    limit,
                    tie_break.as_str()
                ));
            }
        }

        let output_count = frontier.len();
        let replaced_existing = store
            .sets
            .insert(
                output_set.clone(),
                CandidateSet {
                    name: output_set.clone(),
                    created_at_unix_ms: Self::now_unix_ms(),
                    source_seq_ids: input.source_seq_ids,
                    candidates: frontier,
                },
            )
            .is_some();
        self.write_candidate_store(store)?;
        result.messages.push(format!(
            "Computed Pareto frontier from '{}' into '{}' ({} candidate(s), objectives={})",
            input_set,
            output_set,
            output_count,
            compiled
                .iter()
                .map(|objective| format!("{}:{}", objective.metric, objective.direction.as_str()))
                .collect::<Vec<_>>()
                .join(", ")
        ));
        if replaced_existing {
            result.warnings.push(format!(
                "ParetoFrontierCandidateSet output '{}' replaced existing candidate set",
                output_set
            ));
        }
        Ok(())
    }

    pub(super) fn op_upsert_workflow_macro_template(
        &mut self,
        name: String,
        description: Option<String>,
        details_url: Option<String>,
        parameters: Vec<WorkflowMacroTemplateParam>,
        input_ports: Vec<WorkflowMacroTemplatePort>,
        output_ports: Vec<WorkflowMacroTemplatePort>,
        script: String,
        result: &mut OpResult,
    ) -> Result<(), EngineError> {
        let name = Self::normalize_workflow_macro_template_name(&name)?;
        let script = script.trim();
        if script.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "Workflow macro template script cannot be empty".to_string(),
            });
        }
        let details_url = Self::normalize_workflow_macro_template_details_url(details_url)?;

        let mut normalized_parameters = Vec::with_capacity(parameters.len());
        let mut seen = HashSet::new();
        for parameter in parameters {
            let param_name = Self::normalize_workflow_macro_param_name(&parameter.name)?;
            if !seen.insert(param_name.clone()) {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!(
                        "Workflow macro template '{}' contains duplicate parameter '{}'",
                        name, param_name
                    ),
                });
            }
            let default_value = parameter
                .default_value
                .map(|value| value.trim().to_string())
                .filter(|value| !value.is_empty());
            let required = if default_value.is_some() {
                false
            } else {
                parameter.required
            };
            normalized_parameters.push(WorkflowMacroTemplateParam {
                name: param_name,
                default_value,
                required,
            });
        }

        let declared = normalized_parameters
            .iter()
            .map(|parameter| parameter.name.clone())
            .collect::<HashSet<_>>();
        let placeholder_regex =
            Regex::new(r"\$\{([A-Za-z_][A-Za-z0-9_]*)\}").map_err(|e| EngineError {
                code: ErrorCode::Internal,
                message: format!("Could not compile workflow macro placeholder regex: {e}"),
            })?;
        for captures in placeholder_regex.captures_iter(script) {
            if let Some(param_name) = captures.get(1).map(|m| m.as_str()) {
                if !declared.contains(param_name) {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "Workflow macro template '{}' references undeclared parameter '{}' in script",
                            name, param_name
                        ),
                    });
                }
            }
        }

        let normalize_port = |port: WorkflowMacroTemplatePort,
                              direction_label: &str|
         -> Result<WorkflowMacroTemplatePort, EngineError> {
            let port_id = Self::normalize_workflow_macro_param_name(&port.port_id)?;
            let kind = port.kind.trim().to_ascii_lowercase();
            if kind.is_empty() {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!(
                        "Workflow macro template '{}' has {} port '{}' with empty kind",
                        name, direction_label, port_id
                    ),
                });
            }
            let cardinality = match port.cardinality.trim().to_ascii_lowercase().as_str() {
                "" | "one" => "one".to_string(),
                "many" => "many".to_string(),
                other => {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "Workflow macro template '{}' has {} port '{}' with unsupported cardinality '{}' (expected one|many)",
                            name, direction_label, port_id, other
                        ),
                    });
                }
            };
            Ok(WorkflowMacroTemplatePort {
                port_id,
                kind,
                required: port.required,
                cardinality,
                description: port
                    .description
                    .map(|value| value.trim().to_string())
                    .filter(|value| !value.is_empty()),
            })
        };

        let mut normalized_input_ports: Vec<WorkflowMacroTemplatePort> = vec![];
        let mut seen_input_ports: HashSet<String> = HashSet::new();
        for port in input_ports {
            let port = normalize_port(port, "input")?;
            if !seen_input_ports.insert(port.port_id.clone()) {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!(
                        "Workflow macro template '{}' contains duplicate input port '{}'",
                        name, port.port_id
                    ),
                });
            }
            normalized_input_ports.push(port);
        }

        let mut normalized_output_ports: Vec<WorkflowMacroTemplatePort> = vec![];
        let mut seen_output_ports: HashSet<String> = HashSet::new();
        for port in output_ports {
            let port = normalize_port(port, "output")?;
            if !seen_output_ports.insert(port.port_id.clone()) {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!(
                        "Workflow macro template '{}' contains duplicate output port '{}'",
                        name, port.port_id
                    ),
                });
            }
            normalized_output_ports.push(port);
        }

        let now = Self::now_unix_ms();
        let mut store = self.read_workflow_macro_template_store();
        let created_at_unix_ms = store
            .templates
            .get(&name)
            .map(|template| template.created_at_unix_ms)
            .unwrap_or(now);
        let replaced = store
            .templates
            .insert(
                name.clone(),
                WorkflowMacroTemplate {
                    name: name.clone(),
                    description: description
                        .map(|text| text.trim().to_string())
                        .filter(|text| !text.is_empty()),
                    details_url,
                    parameters: normalized_parameters,
                    input_ports: normalized_input_ports,
                    output_ports: normalized_output_ports,
                    template_schema: CLONING_MACRO_TEMPLATE_SCHEMA.to_string(),
                    script: script.to_string(),
                    created_at_unix_ms,
                    updated_at_unix_ms: now,
                },
            )
            .is_some();
        self.write_workflow_macro_template_store(store)?;
        if replaced {
            result
                .messages
                .push(format!("Updated workflow macro template '{}'", name));
        } else {
            result
                .messages
                .push(format!("Added workflow macro template '{}'", name));
        }
        Ok(())
    }

    pub(super) fn op_delete_workflow_macro_template(
        &mut self,
        name: String,
        result: &mut OpResult,
    ) -> Result<(), EngineError> {
        let name = Self::normalize_workflow_macro_template_name(&name)?;
        let mut store = self.read_workflow_macro_template_store();
        let removed = store.templates.remove(&name).is_some();
        self.write_workflow_macro_template_store(store)?;
        if removed {
            result
                .messages
                .push(format!("Deleted workflow macro template '{}'", name));
        } else {
            result.warnings.push(format!(
                "Workflow macro template '{}' was not present",
                name
            ));
        }
        Ok(())
    }

    pub(super) fn op_upsert_candidate_macro_template(
        &mut self,
        name: String,
        description: Option<String>,
        details_url: Option<String>,
        parameters: Vec<CandidateMacroTemplateParam>,
        script: String,
        result: &mut OpResult,
    ) -> Result<(), EngineError> {
        let name = Self::normalize_candidate_macro_template_name(&name)?;
        let script = script.trim();
        if script.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "Candidate macro template script cannot be empty".to_string(),
            });
        }
        let details_url = Self::normalize_candidate_macro_template_details_url(details_url)?;

        let mut normalized_parameters = Vec::with_capacity(parameters.len());
        let mut seen = HashSet::new();
        for parameter in parameters {
            let param_name = Self::normalize_candidate_macro_param_name(&parameter.name)?;
            if !seen.insert(param_name.clone()) {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!(
                        "Candidate macro template '{}' contains duplicate parameter '{}'",
                        name, param_name
                    ),
                });
            }
            let default_value = parameter
                .default_value
                .map(|value| value.trim().to_string())
                .filter(|value| !value.is_empty());
            let required = if default_value.is_some() {
                false
            } else {
                parameter.required
            };
            normalized_parameters.push(CandidateMacroTemplateParam {
                name: param_name,
                default_value,
                required,
            });
        }

        let declared = normalized_parameters
            .iter()
            .map(|parameter| parameter.name.clone())
            .collect::<HashSet<_>>();
        let placeholder_regex =
            Regex::new(r"\$\{([A-Za-z_][A-Za-z0-9_]*)\}").map_err(|e| EngineError {
                code: ErrorCode::Internal,
                message: format!("Could not compile candidate macro placeholder regex: {e}"),
            })?;
        for captures in placeholder_regex.captures_iter(script) {
            if let Some(param_name) = captures.get(1).map(|m| m.as_str()) {
                if !declared.contains(param_name) {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "Candidate macro template '{}' references undeclared parameter '{}' in script",
                            name, param_name
                        ),
                    });
                }
            }
        }

        let now = Self::now_unix_ms();
        let mut store = self.read_candidate_macro_template_store();
        let created_at_unix_ms = store
            .templates
            .get(&name)
            .map(|template| template.created_at_unix_ms)
            .unwrap_or(now);
        let replaced = store
            .templates
            .insert(
                name.clone(),
                CandidateMacroTemplate {
                    name: name.clone(),
                    description: description
                        .map(|text| text.trim().to_string())
                        .filter(|text| !text.is_empty()),
                    details_url,
                    parameters: normalized_parameters,
                    script: script.to_string(),
                    created_at_unix_ms,
                    updated_at_unix_ms: now,
                },
            )
            .is_some();
        self.write_candidate_macro_template_store(store)?;
        if replaced {
            result
                .messages
                .push(format!("Updated candidate macro template '{}'", name));
        } else {
            result
                .messages
                .push(format!("Added candidate macro template '{}'", name));
        }
        Ok(())
    }

    pub(super) fn op_delete_candidate_macro_template(
        &mut self,
        name: String,
        result: &mut OpResult,
    ) -> Result<(), EngineError> {
        let name = Self::normalize_candidate_macro_template_name(&name)?;
        let mut store = self.read_candidate_macro_template_store();
        let removed = store.templates.remove(&name).is_some();
        self.write_candidate_macro_template_store(store)?;
        if removed {
            result
                .messages
                .push(format!("Deleted candidate macro template '{}'", name));
        } else {
            result.warnings.push(format!(
                "Candidate macro template '{}' was not present",
                name
            ));
        }
        Ok(())
    }

}
