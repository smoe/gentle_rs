//! Provider-neutral expression drafts and explicit read-only provider binding.
//! Existing product assessment supplies evidence; this module neither invents
//! requirements nor loads vendor examples, changes sequences, or executes services.

use super::*;
use crate::engine::{
    ProteinExpressionProviderBinding, ProteinExpressionProviderBindingReason as Reason,
    ProteinExpressionProviderBindingStatus as BindingStatus,
};

fn omit_unanswered_fields(value: &mut Value) {
    match value {
        Value::Object(fields) => {
            fields.retain(|_, value| !value.is_null());
            for value in fields.values_mut() {
                omit_unanswered_fields(value);
            }
        }
        Value::Array(values) => {
            for value in values {
                omit_unanswered_fields(value);
            }
        }
        _ => {}
    }
}

fn section(fields: serde_json::Map<String, Value>) -> Option<Value> {
    (!fields.is_empty()).then_some(Value::Object(fields))
}

pub(super) fn protein_expression_neutral_service_draft(
    product: &ProteinExpressionProductDefinition,
    sequence: Option<&DNAsequence>,
    cds: &ProteinExpressionCdsAssessment,
    requirements: Option<&ProteinExpressionRequirements>,
    tags: &ProteinExpressionTagAssessment,
    missing_questions: &[PlanningCloningMissingQuestion],
) -> Result<Option<ExternalServiceDeliveryRouteRequest>, String> {
    let (Some(sequence), Some(seq_id)) = (
        sequence,
        product
            .seq_id
            .as_deref()
            .map(str::trim)
            .filter(|id| !id.is_empty()),
    ) else {
        return Ok(None);
    };
    if !product.sequence_present
        || sequence.is_empty()
        || cds.annotated_cds_features.len() > 1
        || !matches!(
            product.readiness.status.as_str(),
            "whole_sequence_cds_candidate"
                | "annotated_cds_review_required"
                | "protein_sequence_review_required"
        )
    {
        return Ok(None);
    }
    let mut source = json!({
        "kind": "project_sequence",
        "readiness_status": product.readiness.status,
        "readiness": product.readiness,
        "cds_assessment": cds,
        "sequence_content": "project_reference_only",
    });
    source["molecule_type"] = json!(sequence.molecule_type().unwrap_or(
        if sequence.is_protein_sequence() {
            "protein"
        } else {
            "DNA"
        }
    ));
    if let Some(name) = sequence.name().as_ref() {
        source["name"] = json!(name);
    }
    if sequence.is_protein_sequence() {
        source["protein_seq_id"] = json!(seq_id);
        source["length_aa"] = json!(sequence.len());
        source["coding_dna_status"] = json!("unresolved");
        source["reverse_translation_status"] = json!("not_performed");
    } else {
        source["seq_id"] = json!(seq_id);
        source["length_bp"] = json!(sequence.len());
    }
    let mut reviewed = requirements
        .map(serde_json::to_value)
        .transpose()
        .map_err(|error| format!("Could not serialize expression requirements: {error}"))?
        .unwrap_or(Value::Null);
    omit_unanswered_fields(&mut reviewed);
    let mut optimization = serde_json::Map::new();
    for key in [
        "yield_goal",
        "chassis",
        "localization",
        "folding",
        "toxicity_induction",
    ] {
        if let Some(value) = reviewed.get(key) {
            optimization.insert(key.to_string(), value.clone());
        }
    }
    let mut vector = serde_json::Map::new();
    if let Some(policy) = reviewed.get("tag_policy") {
        vector.insert("tag_policy".to_string(), policy.clone());
    }
    if !tags.annotated_tags.is_empty() {
        vector.insert(
            "annotated_tag_evidence".to_string(),
            json!(tags.annotated_tags),
        );
        vector.insert("annotation_is_context_only".to_string(), json!(true));
    }
    let mut delivery = serde_json::Map::new();
    if let Some(goal) = reviewed.get("yield_goal") {
        delivery.insert("yield_goal".to_string(), goal.clone());
    }
    if let Some(Value::Object(endpoint)) = reviewed.get("scale_purification") {
        delivery.extend(endpoint.clone());
    }
    let mut metadata = json!({
        "service_intent": "protein_expression",
        "human_review_only": true,
        "review_status": "review_required",
        "generated_by": "planning protein-expression-handoff",
        "nonclaims": [
            "Project references are not inline sequence content or a designed construct.",
            "Reviewed yield and optimization requirements are requests, not expression, folding or synthesis predictions.",
            "Provider selection is a local review projection, not permission for network I/O, submission, purchase or provider contact."
        ],
        "product_readiness_status": product.readiness.status,
        "unresolved_question_ids": missing_questions.iter().map(|question| &question.question_id).collect::<Vec<_>>(),
    });
    if let Some(requirements) = requirements {
        metadata["requirements_schema"] = json!(requirements.schema);
        metadata["reviewed_requirement_topics"] =
            json!(protein_expression_reviewed_requirement_topics(requirements));
        for key in [
            "requirements_id",
            "reviewed_by",
            "reviewed_at_unix_ms",
            "outsourcing",
            "notes",
        ] {
            if let Some(value) = reviewed.get(key) {
                metadata[key] = value.clone();
            }
        }
    }
    Ok(Some(ExternalServiceDeliveryRouteRequest {
        source_target: source,
        optimization_target: section(optimization),
        vector_spec: section(vector),
        delivery_options: section(delivery),
        request_metadata: Some(metadata),
        ..ExternalServiceDeliveryRouteRequest::default()
    }))
}

pub(super) fn protein_expression_provider_projection(
    provider: Option<&str>,
    product: &ProteinExpressionProductDefinition,
    sequence: Option<&DNAsequence>,
    cds: &ProteinExpressionCdsAssessment,
    requirements: Option<&ProteinExpressionRequirements>,
    draft: Option<&ExternalServiceDeliveryRouteRequest>,
) -> (
    ProteinExpressionProviderBinding,
    Option<ExternalServiceRequest>,
) {
    let Some(provider) = provider else {
        return (ProteinExpressionProviderBinding::default(), None);
    };
    let provider = crate::external_service_providers::normalize_service_token(provider);
    let mut binding = ProteinExpressionProviderBinding {
        requested_provider: Some(provider.clone()),
        status: BindingStatus::Withheld,
        ..ProteinExpressionProviderBinding::default()
    };
    if !product.sequence_present
        || sequence.is_none_or(DNAsequence::is_empty)
        || cds.annotated_cds_features.len() > 1
    {
        binding.reasons.push(Reason::ProductIdentityUnresolved);
    }
    if !matches!(
        product.readiness.status.as_str(),
        "whole_sequence_cds_candidate"
            | "annotated_cds_review_required"
            | "protein_sequence_review_required"
    ) {
        binding.reasons.push(Reason::ProductNotReady);
    }
    let outsourcing = requirements.and_then(|requirements| requirements.outsourcing.as_ref());
    match outsourcing.and_then(|outsourcing| outsourcing.allowed) {
        Some(true) => {}
        Some(false) => binding.reasons.push(Reason::OutsourcingDisallowed),
        None => binding.reasons.push(Reason::OutsourcingPermissionMissing),
    }
    if let Some(outsourcing) = outsourcing
        && !outsourcing.approved_providers.is_empty()
        && !outsourcing.approved_providers.iter().any(|approved| {
            crate::external_service_providers::normalize_service_token(approved) == provider
        })
    {
        binding.reasons.push(Reason::ProviderNotApproved);
    }
    if !binding.reasons.is_empty() {
        return (binding, None);
    }
    let Some(draft) = draft else {
        binding.reasons.push(Reason::ProductNotReady);
        return (binding, None);
    };
    match service_readiness::project_external_service_delivery_request(
        draft,
        &provider,
        "protein_expression",
    ) {
        Ok(request) => {
            binding.status = BindingStatus::Projected;
            (binding, Some(request))
        }
        Err(error) => {
            binding.reasons.push(Reason::ProviderUnavailable);
            binding.detail = Some(error);
            (binding, None)
        }
    }
}
