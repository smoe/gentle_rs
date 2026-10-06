//! Root context/approval boundary for the dependency-free synonymous design core.
//! Search is read-only. Apply validates the approved bytes, never reruns search,
//! and deliberately drops source biological annotations rather than copying claims.
//! Submitted search claims are retained as unverified evidence, not engine proof.

use super::*;
use gentle_protocol::sequence_design::*;
use gentle_sequence_design as core;

const MAX_SOURCE_ID_BYTES: usize = 4_096;
const MAX_DIGEST_BYTES: usize = 128;
const APPLY_SEARCH_NONCLAIM: &str = "Apply verified output constraints and exact edits only; preview search history, completeness and minimum edits remain unverified. A content digest is not authenticated search provenance.";

fn interval(value: &DesignInterval) -> core::Interval {
    core::Interval {
        start: value.start_0based,
        end: value.end_0based_exclusive,
    }
}

fn strand(value: DesignStrand) -> core::Strand {
    match value {
        DesignStrand::Forward => core::Strand::Forward,
        DesignStrand::Reverse => core::Strand::Reverse,
        DesignStrand::Both => core::Strand::Both,
    }
}

fn wire_strand(value: core::Strand) -> DesignStrand {
    match value {
        core::Strand::Forward => DesignStrand::Forward,
        core::Strand::Reverse => DesignStrand::Reverse,
        core::Strand::Both => DesignStrand::Both,
    }
}

fn wire_edits(edits: Vec<core::Edit>) -> Vec<DesignNucleotideEdit> {
    edits
        .into_iter()
        .map(|edit| DesignNucleotideEdit {
            position_0based: edit.position,
            before: (edit.before as char).to_string(),
            after: (edit.after as char).to_string(),
        })
        .collect()
}

fn search_strategy(value: DnaDesignSearchStrategy) -> core::SearchStrategy {
    match value {
        DnaDesignSearchStrategy::FullEnumeration => core::SearchStrategy::FullEnumeration,
        DnaDesignSearchStrategy::ConflictDirected => core::SearchStrategy::ConflictDirected,
    }
}

fn nonclaims(strategy: DnaDesignSearchStrategy) -> Vec<String> {
    let mut statements = vec![
        "Synthetic coding-insert redesign only; no natural assay target is implicitly recoded.".into(),
        "Translation preservation does not establish expression, splicing, regulatory function, folding or experimental suitability.".into(),
        "Standard code 1 elongation and literal frozen ATG/stop only; no inferred initiation efficiency or authentic CDS claim.".into(),
        "Bounded enumeration is not full DNA Chisel compatibility or an efficient general optimizer; incomplete results do not prove minimum edits.".into(),
        "Source annotations are omitted from the derived sequence; only the explicitly supplied synthetic CDS is recorded.".into(),
    ];
    if strategy == DnaDesignSearchStrategy::ConflictDirected {
        statements[3] = "Bounded conflict-directed search is not full DNA Chisel compatibility; only completed search establishes minimum edits and no runtime-performance claim is made.".into();
    }
    statements
}

impl GentleEngine {
    fn dna_design_request_limits(request: &DnaSequenceDesignRequest) -> Result<(), EngineError> {
        if request.schema != REQUEST_SCHEMA || request.purpose != "synthetic_coding_insert" {
            return Err(EngineError::invalid_input(
                "Declare gentle.dna_sequence_design_request.v1 and purpose=synthetic_coding_insert; labels do not authorize recoding",
            ));
        }
        let oversized_target = match &request.target {
            DnaDesignTarget::InlineSequence { sequence } => sequence.len() > core::MAX_SEQUENCE_BP,
            DnaDesignTarget::LoadedSequence { seq_id } => seq_id.len() > MAX_SOURCE_ID_BYTES,
        };
        if oversized_target
            || request.protected_intervals.len() > core::MAX_PROTECTED_INTERVALS
            || request.avoid_motifs.len() > core::MAX_MOTIFS
            || request
                .avoid_motifs
                .iter()
                .any(|m| m.pattern.len() > core::MAX_MOTIF_BP)
            || request.protein_sequence.len() > core::MAX_SEQUENCE_BP / 3
            || request.output_seq_id.len() > 128
        {
            return Err(EngineError::invalid_input(
                "Sequence-design request exceeds bounded input limits",
            ));
        }
        Ok(())
    }

    fn dna_design_hash(value: &impl Serialize) -> Result<String, EngineError> {
        serde_json::to_vec(value)
            .map(|bytes| sha256_prefixed_bytes(&bytes))
            .map_err(|error| {
                EngineError::invalid_input(format!("Sequence-design serialization failed: {error}"))
            })
    }

    fn dna_design_source(
        &self,
        request: &DnaSequenceDesignRequest,
    ) -> Result<(String, Option<String>, usize), EngineError> {
        Self::dna_design_request_limits(request)?;
        let id = &request.output_seq_id;
        if id.is_empty()
            || id.len() > 128
            || !id
                .bytes()
                .all(|b| b.is_ascii_alphanumeric() || b"_-.".contains(&b))
        {
            return Err(EngineError::invalid_input(
                "Design output ID requires 1..128 ASCII letters, digits, '_', '-' or '.'",
            ));
        }
        if self.state.sequences.contains_key(id) {
            return Err(EngineError::invalid_input(
                "Sequence design never overwrites an existing output ID",
            ));
        }
        match &request.target {
            DnaDesignTarget::InlineSequence { sequence } => {
                if sequence.len() > core::MAX_SEQUENCE_BP {
                    return Err(EngineError::invalid_input(
                        "Sequence design admits at most 12000 bases",
                    ));
                }
                Ok((sequence.clone(), None, 0))
            }
            DnaDesignTarget::LoadedSequence { seq_id } => {
                let dna = self.state.sequences.get(seq_id).ok_or_else(|| {
                    EngineError::invalid_input(format!(
                        "Design source sequence '{seq_id}' is not loaded"
                    ))
                })?;
                if dna.is_circular()
                    || dna.len() > core::MAX_SEQUENCE_BP
                    || dna
                        .molecule_type()
                        .is_some_and(|kind| !kind.eq_ignore_ascii_case("DNA"))
                {
                    return Err(EngineError::invalid_input(
                        "Design source must be linear DNA of at most 12000 bases; circular/RNA/protein sources are unsupported",
                    ));
                }
                Ok((
                    dna.get_forward_string(),
                    Some(Self::dna_design_hash(&dna.features())?),
                    dna.features().len(),
                ))
            }
        }
    }

    fn dna_design_input(
        request: &DnaSequenceDesignRequest,
        sequence: String,
    ) -> Result<(core::Input, Vec<DesignCodon>), EngineError> {
        let mut codons = BTreeMap::new();
        let mut mapping = vec![];
        // Explicit codon keys remove dependence on the asset's TCAG positional order.
        for a in b"ACGT" {
            for b in b"ACGT" {
                for c in b"ACGT" {
                    let codon = [*a, *b, *c];
                    let residue = crate::AMINO_ACIDS.codon2aa(codon, Some(request.genetic_code));
                    let residue = if residue == crate::amino_acids::STOP_CODON {
                        '*'
                    } else {
                        residue
                    };
                    codons.insert(codon, residue as u8);
                    mapping.push(DesignCodon {
                        codon: String::from_utf8(codon.to_vec()).unwrap(),
                        residue: residue.to_string(),
                    });
                }
            }
        }
        Ok((
            core::Input {
                sequence,
                cds: interval(&request.cds),
                protein: request.protein_sequence.clone(),
                code: core::GeneticCode {
                    id: request.genetic_code,
                    codons,
                },
                protected: request.protected_intervals.iter().map(interval).collect(),
                motifs: request
                    .avoid_motifs
                    .iter()
                    .map(|m| core::Motif {
                        pattern: m.pattern.clone(),
                        strand: strand(m.strand),
                    })
                    .collect(),
                max_evaluations: request.max_evaluations,
            },
            mapping,
        ))
    }

    fn dna_design_approval(report: &DnaSequenceDesignReport) -> Result<String, EngineError> {
        let mut basis = report.clone();
        basis.approval_digest = None;
        Self::dna_design_hash(&basis)
    }

    /// Preview without altering source/project state. Cancellation discards any candidate.
    pub fn plan_dna_sequence_design(
        &self,
        request: DnaSequenceDesignRequest,
        cancel: impl FnMut(u64) -> bool,
    ) -> Result<DnaSequenceDesignReport, EngineError> {
        let (sequence, source_features_sha256, omitted_source_feature_count) =
            self.dna_design_source(&request)?;
        let (input, mapping) = Self::dna_design_input(&request, sequence.clone())?;
        let strategy = request.search_strategy;
        let result = core::solve_with_strategy(&input, search_strategy(strategy), cancel);
        let status = match result.status {
            core::Status::Invalid => DnaDesignStatus::Invalid,
            core::Status::Unsupported => DnaDesignStatus::Unsupported,
            core::Status::Feasible => DnaDesignStatus::Feasible,
            core::Status::SearchExhausted => DnaDesignStatus::SearchExhausted,
            core::Status::ProvenInfeasible => DnaDesignStatus::ProvenInfeasible,
            core::Status::Cancelled => DnaDesignStatus::Cancelled,
        };
        let mut report = DnaSequenceDesignReport {
            schema: REPORT_SCHEMA.into(),
            request,
            source_sha256: sha256_prefixed_str(&sequence),
            source_sequence: sequence,
            source_features_sha256,
            omitted_source_feature_count,
            genetic_code_mapping_sha256: Self::dna_design_hash(&mapping)?,
            genetic_code_mapping: mapping,
            algorithm: search_strategy(strategy).algorithm().into(),
            status,
            reason: result.reason,
            output_sha256: result.sequence.as_deref().map(sha256_prefixed_str),
            output_sequence: result.sequence,
            edits: wire_edits(result.edits),
            initial_matches: result
                .initial_matches
                .into_iter()
                .map(|h| DesignMotifMatch {
                    motif_index: h.motif_index,
                    interval: DesignInterval {
                        start_0based: h.interval.start,
                        end_0based_exclusive: h.interval.end,
                    },
                    strand: wire_strand(h.strand),
                })
                .collect(),
            evaluated_candidates: result.evaluated_candidates,
            effective_evaluation_budget: result.effective_evaluation_budget,
            search_space: result.search_space,
            optimization_complete: result.optimization_complete,
            minimum_edits_proven: result.minimum_edits_proven,
            approval_digest: None,
            nonclaims: nonclaims(strategy),
        };
        if status == DnaDesignStatus::Feasible {
            report.approval_digest = Some(Self::dna_design_approval(&report)?);
        }
        Ok(report)
    }

    /// Execute the shared apply operation with one provenance/undo boundary.
    /// Output constraints are revalidated; search claims are not re-proven.
    pub fn apply_dna_sequence_design(
        &mut self,
        proposal: DnaSequenceDesignReport,
        approval: &str,
    ) -> Result<DnaSequenceDesignReceipt, EngineError> {
        self.apply(Operation::ApplyDnaSequenceDesign {
            proposal: Box::new(proposal),
            approval_digest: approval.into(),
        })?
        .dna_sequence_design_receipt
        .ok_or_else(|| EngineError {
            code: ErrorCode::Internal,
            message: "Sequence-design apply returned no receipt".into(),
            cause_chain: vec![],
        })
    }

    /// Materialize validated bytes without nested public-operation hooks.
    /// The outer apply owns lineage, container creation, journaling and undo.
    pub(super) fn materialize_approved_dna_design(
        &mut self,
        proposal: DnaSequenceDesignReport,
        approval: &str,
    ) -> Result<DnaSequenceDesignReceipt, EngineError> {
        Self::dna_design_request_limits(&proposal.request)?;
        if proposal.source_sequence.len() > core::MAX_SEQUENCE_BP
            || proposal.source_sha256.len() > MAX_DIGEST_BYTES
            || proposal.genetic_code_mapping_sha256.len() > MAX_DIGEST_BYTES
            || proposal
                .source_features_sha256
                .as_ref()
                .is_some_and(|s| s.len() > MAX_DIGEST_BYTES)
            || proposal
                .output_sha256
                .as_ref()
                .is_some_and(|s| s.len() > MAX_DIGEST_BYTES)
            || proposal
                .output_sequence
                .as_ref()
                .is_some_and(|s| s.len() > core::MAX_SEQUENCE_BP)
            || proposal.genetic_code_mapping.len() != 64
            || proposal
                .genetic_code_mapping
                .iter()
                .any(|c| c.codon.len() != 3 || c.residue.len() != 1)
            || proposal.edits.len() > core::MAX_SEQUENCE_BP
            || proposal
                .edits
                .iter()
                .any(|e| e.before.len() != 1 || e.after.len() != 1)
            || proposal.initial_matches.len() > core::MAX_REPORTED_MATCHES
            || proposal.nonclaims != nonclaims(proposal.request.search_strategy)
            || proposal.reason.len() > 1024
            || proposal.algorithm.len() > 128
            || approval.len() > MAX_DIGEST_BYTES
        {
            return Err(EngineError::invalid_input(
                "Sequence-design preview exceeds bounded input limits or changes policy",
            ));
        }
        if proposal.schema != REPORT_SCHEMA
            || proposal.algorithm != search_strategy(proposal.request.search_strategy).algorithm()
            || proposal.status != DnaDesignStatus::Feasible
            || proposal.approval_digest.as_deref() != Some(approval)
            || Self::dna_design_approval(&proposal)? != approval
        {
            return Err(EngineError::invalid_input(
                "Sequence-design approval must match an exact feasible preview; cancelled/exhausted/altered proposals are not applicable",
            ));
        }
        let (source, features_digest, feature_count) = self.dna_design_source(&proposal.request)?;
        let (input, mapping) = Self::dna_design_input(&proposal.request, source.clone())?;
        if source != proposal.source_sequence
            || sha256_prefixed_str(&source) != proposal.source_sha256
            || features_digest != proposal.source_features_sha256
            || feature_count != proposal.omitted_source_feature_count
            || mapping != proposal.genetic_code_mapping
            || Self::dna_design_hash(&mapping)? != proposal.genetic_code_mapping_sha256
            || proposal.nonclaims != nonclaims(proposal.request.search_strategy)
        {
            return Err(EngineError::invalid_input(
                "Sequence-design source, annotations, genetic-code mapping or policy changed; prepare and review a new preview",
            ));
        }
        let output = proposal.output_sequence.as_deref().ok_or_else(|| {
            EngineError::invalid_input("Feasible proposal has no exact output DNA")
        })?;
        if output.len() > core::MAX_SEQUENCE_BP
            || proposal.output_sha256.as_deref() != Some(sha256_prefixed_str(output).as_str())
        {
            return Err(EngineError::invalid_input(
                "Approved sequence-design output digest/size mismatch",
            ));
        }
        let validated_edits =
            core::validate_output(&input, output).map_err(EngineError::invalid_input)?;
        if wire_edits(validated_edits) != proposal.edits {
            return Err(EngineError::invalid_input(
                "Approved nucleotide edit script disagrees with full output validation",
            ));
        }
        let mut receipt_nonclaims = proposal.nonclaims.clone();
        receipt_nonclaims.push(APPLY_SEARCH_NONCLAIM.into());
        let receipt = DnaSequenceDesignReceipt {
            schema: RECEIPT_SCHEMA.into(),
            approval_digest: approval.into(),
            created_seq_id: proposal.request.output_seq_id.clone(),
            output_sha256: sha256_prefixed_str(output),
            output_constraints_verified: true,
            search_claims_verified: false,
            nonclaims: receipt_nonclaims,
        };
        let mut dna = DNAsequence::from_sequence(output).map_err(|error| {
            EngineError::invalid_input(format!("Could not materialize approved DNA: {error}"))
        })?;
        dna.set_circular(false);
        dna.set_name(Some(receipt.created_seq_id.clone()));
        Self::prepare_sequence(&mut dna);
        if sha256_prefixed_str(&dna.get_forward_string()) != receipt.output_sha256 {
            return Err(EngineError::invalid_input(
                "Materialized DNA differs from the approved output",
            ));
        }
        *dna.features_mut() = vec![gb_io::seq::Feature {
            kind: "CDS".into(),
            location: gb_io::seq::Location::simple_range(
                proposal.request.cds.start_0based as i64,
                proposal.request.cds.end_0based_exclusive as i64,
            ),
            qualifiers: vec![
                (
                    "label".into(),
                    Some("Synthetic CDS; supplied translation preserved".into()),
                ),
                (
                    "translation".into(),
                    Some(proposal.request.protein_sequence.clone()),
                ),
                ("transl_table".into(), Some("1".into())),
                ("gentle_design_approval".into(), Some(approval.into())),
                ("note".into(), Some(receipt.nonclaims.join(" "))),
            ],
        }];
        let mut detached = self.fork_detached_execution();
        let engine = detached.engine_mut();
        engine
            .state_mut()
            .sequences
            .insert(receipt.created_seq_id.clone(), dna);
        engine.state.metadata.insert(
            format!("dna_sequence_design:{}", receipt.created_seq_id),
            serde_json::json!({"submitted_proposal": proposal, "receipt": receipt}),
        );
        self.commit_detached_execution(&mut detached)?;
        Ok(receipt)
    }
}

#[cfg(test)]
#[path = "sequence_design/tests.rs"]
mod tests;
