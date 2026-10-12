//! Separately approved exact design products for ordered regulatory-fragment panels.

use super::*;
use gb_io::seq::{Feature, Location, Seq};

const PROPOSAL_SCHEMA: &str = "gentle.regulatory_fragment_materialization_proposal.v1";
const RECEIPT_SCHEMA: &str = "gentle.regulatory_fragment_materialization_receipt.v1";
const METHOD: &str = "exact_insertion_context_replacement_v1";

impl GentleEngine {
    /// Prepare exact products without mutating state. This is sequence engineering,
    /// not a claim that a particular cloning reaction will produce these molecules.
    pub fn plan_regulatory_fragment_materialization(
        &self,
        plan: RegulatoryFragmentPanelPlan,
        output_prefix: String,
    ) -> Result<RegulatoryFragmentMaterializationProposal, EngineError> {
        self.validate_regulatory_fragment_panel_approval(&plan, &plan.proposal_digest)?;
        let prefix = Self::regulatory_fragment_slug(&output_prefix);
        if prefix.is_empty() || prefix != output_prefix || !plan.blockers.is_empty() {
            return Err(Self::regulatory_fragment_error(
                "materialization_not_ready",
                [&plan.plan_id],
                "Use a non-empty normalized output prefix and resolve plan blockers first.",
            ));
        }
        let context = &plan.vector_context;
        let vector = self
            .state
            .sequences
            .get(&context.vector_seq_id)
            .ok_or_else(|| {
                Self::regulatory_fragment_error(
                    "vector_missing",
                    [&context.vector_seq_id],
                    "Vector is not loaded",
                )
            })?;
        let sequence = vector.get_forward_string();
        let start = context.insertion_start_0based;
        let end = context.insertion_end_0based_exclusive;
        if start >= end || end > sequence.len() {
            return Err(Self::regulatory_fragment_error(
                "invalid_insertion_context",
                [&plan.plan_id],
                "Expected an exact, non-wrapping vector insertion interval",
            ));
        }
        let mut products = vec![];
        let mut ids = BTreeSet::new();
        for member in &plan.members {
            let output_seq_id = format!("{prefix}_{}", member.member_id);
            if !ids.insert(output_seq_id.clone())
                || self.state.sequences.contains_key(&output_seq_id)
            {
                return Err(Self::regulatory_fragment_error(
                    "output_id_exists",
                    [&output_seq_id],
                    "Materialization never overwrites or renames an output.",
                ));
            }
            let insert = &member.insert_sequence_5prime_to_3prime;
            // CreateSequenceFromText stores canonical uppercase DNA. Bind those
            // exact product bytes at proposal time, not only after approval.
            let full = format!("{}{}{}", &sequence[..start], insert, &sequence[end..])
                .to_ascii_uppercase();
            let mut features = vec![];
            let mut omitted = vec![];
            let mut coordinate_frame = Seq::empty();
            coordinate_frame.seq = full.as_bytes().to_vec();
            // Features crossing the replaced interval must not silently acquire a new
            // biological meaning. Preserve only wholly unaffected records and list omissions.
            for (index, feature) in vector.features().iter().enumerate() {
                let (low, high) = feature.location.find_bounds().map_err(|e| {
                    Self::regulatory_fragment_error(
                        "unsupported_vector_annotation",
                        [&output_seq_id],
                        e.to_string(),
                    )
                })?;
                if low < 0 || high < low || high as usize > sequence.len() {
                    return Err(Self::regulatory_fragment_error(
                        "invalid_vector_annotation",
                        [&output_seq_id],
                        "Vector feature lies outside its source sequence",
                    ));
                }
                if high <= start as i64 {
                    features.push(feature.clone());
                } else if low >= end as i64 {
                    features.push(
                        coordinate_frame
                            .relocate_feature(
                                feature.clone(),
                                insert.len() as i64 - (end - start) as i64,
                            )
                            .map_err(|e| {
                                Self::regulatory_fragment_error(
                                    "feature_projection_failed",
                                    [&output_seq_id],
                                    e.to_string(),
                                )
                            })?,
                    );
                } else {
                    omitted.push(index);
                }
            }
            for (ordinal, instance) in member.instances.iter().enumerate() {
                let lo = start + instance.assembled_start_0based;
                let hi = start + instance.assembled_end_0based_exclusive;
                let qualifier = |key: &'static str, value: String| (key.into(), Some(value));
                features.push(Feature {
                    kind: "misc_feature".into(),
                    location: Location::simple_range(lo as i64, hi as i64),
                    qualifiers: vec![
                        qualifier("label", format!("{} instance {}", instance.fragment_id, ordinal + 1)),
                        qualifier("gentle_fragment_instance", serde_json::to_string(instance).map_err(|e|
                            Self::regulatory_fragment_error("instance_serialization_failed", [&output_seq_id], e.to_string()))?),
                        qualifier("gentle_plan_digest", plan.proposal_digest.clone()),
                        qualifier("note", "Engineered fragment placement; orientation is relative to its canonical ROI, not a newly inferred gene strand.".into()),
                    ],
                });
                if !instance.spacer_before.is_empty() {
                    features.push(Feature {
                        kind: "misc_feature".into(),
                        location: Location::simple_range(
                            (lo - instance.spacer_before.len()) as i64,
                            lo as i64,
                        ),
                        qualifiers: vec![
                            qualifier("label", format!("spacer before instance {}", ordinal + 1)),
                            qualifier(
                                "note",
                                "Explicit designed spacer; no regulatory function asserted".into(),
                            ),
                        ],
                    });
                }
            }
            products.push(RegulatoryFragmentDesignedProduct {
                member_id: member.member_id.clone(),
                output_seq_id,
                sequence_sha256: sha256_prefixed_str(&full),
                sequence_5prime_to_3prime: full,
                circular: vector.is_circular(),
                features,
                omitted_vector_feature_indices: omitted,
                instances: member.instances.clone(),
            });
        }
        let mut proposal = RegulatoryFragmentMaterializationProposal {
            schema: PROPOSAL_SCHEMA.into(), proposal_digest: String::new(), plan: Box::new(plan),
            output_prefix, method: METHOD.into(),
            vector_features_sha256: Self::regulatory_fragment_value_sha256(
                &vector.features(), "source vector annotations")?,
            products,
            nonclaims: vec![
                "Exact designed vector-context replacement, not a simulated restriction/Gibson reaction, validated product, or order-ready panel.".into(),
                "Features crossing the replaced vector context are explicitly omitted; source fragment annotations are not inferred across engineered junctions.".into(),
                "Functional activity, sufficiency and partner dependence require experimental comparison.".into(),
            ],
        };
        proposal.proposal_digest = Self::regulatory_fragment_materialization_digest(&proposal)?;
        Ok(proposal)
    }

    fn regulatory_fragment_materialization_digest(
        proposal: &RegulatoryFragmentMaterializationProposal,
    ) -> Result<String, EngineError> {
        let mut basis = proposal.clone();
        basis.proposal_digest.clear();
        Self::regulatory_fragment_value_sha256(&basis, "exact regulatory products")
    }

    /// Commit the exact ordered products atomically, with one undoable outer operation.
    pub fn materialize_regulatory_fragment_panel(
        &mut self,
        proposal: RegulatoryFragmentMaterializationProposal,
        approval_digest: &str,
    ) -> Result<RegulatoryFragmentMaterializationReceipt, EngineError> {
        if proposal.schema != PROPOSAL_SCHEMA
            || proposal.method != METHOD
            || approval_digest != proposal.proposal_digest
            || Self::regulatory_fragment_materialization_digest(&proposal)?
                != proposal.proposal_digest
        {
            return Err(Self::regulatory_fragment_error(
                "materialization_approval_mismatch",
                [&proposal.plan.plan_id],
                "Approval must bind the exact ordered products, features and source plan.",
            ));
        }
        let fresh = self.plan_regulatory_fragment_materialization(
            *proposal.plan.clone(),
            proposal.output_prefix.clone(),
        )?;
        if fresh.proposal_digest != proposal.proposal_digest {
            return Err(Self::regulatory_fragment_error(
                "materialization_source_stale",
                [&proposal.plan.plan_id],
                "Bound evidence, vector annotation or planned products changed. Re-plan and review.",
            ));
        }
        let mut detached = self.fork_detached_execution();
        let engine = detached.engine_mut();
        let mut created = vec![];
        for product in &proposal.products {
            let result = engine.apply(Operation::CreateSequenceFromText {
                sequence_text: product.sequence_5prime_to_3prime.clone(),
                output_id: Some(product.output_seq_id.clone()),
                name: Some(product.member_id.clone()),
                circular: product.circular,
            })?;
            if result.created_seq_ids != [product.output_seq_id.clone()] {
                return Err(Self::regulatory_fragment_error(
                    "unexpected_product_identity",
                    [&product.output_seq_id],
                    "Exact output identity was not preserved",
                ));
            }
            let dna = engine
                .state
                .sequences
                .get_mut(&product.output_seq_id)
                .ok_or_else(|| {
                    Self::regulatory_fragment_error(
                        "product_missing",
                        [&product.output_seq_id],
                        "Product was not created",
                    )
                })?;
            *dna.features_mut() = product.features.clone();
            if sha256_prefixed_str(&dna.get_forward_string()) != product.sequence_sha256 {
                return Err(Self::regulatory_fragment_error(
                    "product_sequence_mismatch",
                    [&product.output_seq_id],
                    "Product sequence differs from the approved bytes",
                ));
            }
            created.push(product.output_seq_id.clone());
        }
        let receipt = RegulatoryFragmentMaterializationReceipt {
            schema: RECEIPT_SCHEMA.into(),
            approved_proposal_digest: proposal.proposal_digest.clone(),
            plan_digest: proposal.plan.proposal_digest.clone(),
            created_seq_ids: created,
            product_sequence_sha256: proposal
                .products
                .iter()
                .map(|p| (p.output_seq_id.clone(), p.sequence_sha256.clone()))
                .collect(),
            final_product_audit_state: RegulatoryFragmentFinalProductAuditState::NotEvaluated,
            nonclaims: proposal.nonclaims.clone(),
        };
        engine.state.metadata.insert(
            format!(
                "regulatory_fragment_materialization:{}",
                proposal.output_prefix
            ),
            serde_json::json!({"proposal": proposal, "receipt": receipt}),
        );
        self.commit_detached_execution(&mut detached)?;
        Ok(receipt)
    }
}
