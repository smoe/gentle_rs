//! Portable biological-reference context for collection members.
//!
//! Coordinates and biological identifiers are only interpretable against a
//! taxon, reference assembly, annotation release, and identifier namespace.
//! Reports own a context registry and members refer to entries by id so one
//! shared context is not duplicated across every row.

use serde::{Deserialize, Serialize};
use std::{collections::BTreeSet, error::Error, fmt};

/// Default context id generated when legacy report-level fields are promoted.
pub const DEFAULT_BIOLOGICAL_CONTEXT_ID: &str = "default";
/// Synthetic id used while resolving legacy fields that have no registry row.
pub const LEGACY_BIOLOGICAL_CONTEXT_ID: &str = "legacy_report_context";

/// One biological and reference-resource context.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct BiologicalContext {
    pub context_id: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub organism: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub taxon_id: Option<String>,
    /// GENtle genome-catalog identity; this is not implicitly an assembly id.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub genome_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub assembly_accession: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub annotation_source: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub annotation_release: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub symbol_namespace: Option<String>,
}

impl BiologicalContext {
    /// True when the row carries no biological/reference information.
    pub fn is_empty(&self) -> bool {
        self.organism.is_none()
            && self.taxon_id.is_none()
            && self.genome_id.is_none()
            && self.assembly_accession.is_none()
            && self.annotation_source.is_none()
            && self.annotation_release.is_none()
            && self.symbol_namespace.is_none()
    }

    /// Compare context meaning while ignoring report-local `context_id`.
    pub fn semantically_matches(&self, other: &Self) -> bool {
        self.organism == other.organism
            && self.taxon_id == other.taxon_id
            && self.genome_id == other.genome_id
            && self.assembly_accession == other.assembly_accession
            && self.annotation_source == other.annotation_source
            && self.annotation_release == other.annotation_release
            && self.symbol_namespace == other.symbol_namespace
    }
}

/// Context rows owned by one portable report.
///
/// This record is flattened into report JSON so `contexts` and
/// `default_context_id` remain directly discoverable machine-contract fields.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct BiologicalContextRegistry {
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub contexts: Vec<BiologicalContext>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub default_context_id: Option<String>,
}

/// Deterministic context-registry/reference failure.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(tag = "reason", rename_all = "snake_case")]
pub enum BiologicalContextResolutionError {
    EmptyContextId,
    DuplicateContextId {
        context_id: String,
    },
    UnknownContextId {
        context_id: String,
    },
    ConflictingField {
        field: String,
        primary_context_id: String,
        primary_value: String,
        fallback_context_id: String,
        fallback_value: String,
    },
}

impl fmt::Display for BiologicalContextResolutionError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptyContextId => write!(formatter, "biological context id must not be empty"),
            Self::DuplicateContextId { context_id } => {
                write!(formatter, "duplicate biological context id '{context_id}'")
            }
            Self::UnknownContextId { context_id } => {
                write!(formatter, "unknown biological context id '{context_id}'")
            }
            Self::ConflictingField {
                field,
                primary_context_id,
                primary_value,
                fallback_context_id,
                fallback_value,
            } => write!(
                formatter,
                "biological contexts '{primary_context_id}' and '{fallback_context_id}' conflict for {field}: '{primary_value}' versus '{fallback_value}'"
            ),
        }
    }
}

impl Error for BiologicalContextResolutionError {}

impl BiologicalContextRegistry {
    /// Validate context ids and the optional default reference.
    pub fn validate(&self) -> Result<(), BiologicalContextResolutionError> {
        let mut ids = BTreeSet::new();
        for context in &self.contexts {
            let context_id = context.context_id.trim();
            if context_id.is_empty() {
                return Err(BiologicalContextResolutionError::EmptyContextId);
            }
            if !ids.insert(context_id) {
                return Err(BiologicalContextResolutionError::DuplicateContextId {
                    context_id: context.context_id.clone(),
                });
            }
        }
        if let Some(default_context_id) = self.default_context_id.as_deref()
            && !ids.contains(default_context_id.trim())
        {
            return Err(BiologicalContextResolutionError::UnknownContextId {
                context_id: default_context_id.to_string(),
            });
        }
        Ok(())
    }

    /// Resolve one member context, then the report default, then legacy fields.
    ///
    /// Lower-precedence rows may fill absent fields but may not contradict a
    /// value already supplied by a higher-precedence row.
    pub fn resolve(
        &self,
        member_context_id: Option<&str>,
        legacy_context: Option<&BiologicalContext>,
    ) -> Result<Option<BiologicalContext>, BiologicalContextResolutionError> {
        self.validate()?;
        let member_context_id = member_context_id
            .map(str::trim)
            .filter(|context_id| !context_id.is_empty());
        let mut resolved = member_context_id
            .map(|context_id| self.context(context_id).cloned())
            .transpose()?;

        if let Some(default_context_id) = self.default_context_id.as_deref()
            && member_context_id != Some(default_context_id)
        {
            let fallback = self.context(default_context_id)?;
            merge_context(&mut resolved, fallback)?;
        }
        if let Some(legacy_context) = legacy_context.filter(|context| !context.is_empty()) {
            merge_context(&mut resolved, legacy_context)?;
        }
        Ok(resolved.filter(|context| !context.is_empty()))
    }

    /// Promote legacy report-level fields into one deterministic default row.
    pub fn ensure_default_from_legacy(
        &mut self,
        legacy_context: &BiologicalContext,
    ) -> Result<(), BiologicalContextResolutionError> {
        if legacy_context.is_empty() {
            return self.validate();
        }
        self.validate()?;
        let default_context_id = self
            .default_context_id
            .clone()
            .or_else(|| (self.contexts.len() == 1).then(|| self.contexts[0].context_id.clone()))
            .unwrap_or_else(|| DEFAULT_BIOLOGICAL_CONTEXT_ID.to_string());

        if let Some(existing) = self
            .contexts
            .iter_mut()
            .find(|context| context.context_id == default_context_id)
        {
            let mut resolved = Some(existing.clone());
            merge_context(&mut resolved, legacy_context)?;
            *existing = resolved.expect("existing biological context remains present");
            existing.context_id = default_context_id.clone();
        } else {
            let mut context = legacy_context.clone();
            context.context_id = default_context_id.clone();
            self.contexts.push(context);
        }
        self.default_context_id = Some(default_context_id);
        self.validate()
    }

    fn context(
        &self,
        context_id: &str,
    ) -> Result<&BiologicalContext, BiologicalContextResolutionError> {
        self.contexts
            .iter()
            .find(|context| context.context_id == context_id)
            .ok_or_else(|| BiologicalContextResolutionError::UnknownContextId {
                context_id: context_id.to_string(),
            })
    }
}

fn merge_context(
    resolved: &mut Option<BiologicalContext>,
    fallback: &BiologicalContext,
) -> Result<(), BiologicalContextResolutionError> {
    let Some(primary) = resolved.as_mut() else {
        *resolved = Some(fallback.clone());
        return Ok(());
    };
    let primary_context_id = primary.context_id.clone();
    let fallback_context_id = fallback.context_id.clone();
    merge_field(
        &mut primary.organism,
        &fallback.organism,
        "organism",
        &primary_context_id,
        &fallback_context_id,
    )?;
    merge_field(
        &mut primary.taxon_id,
        &fallback.taxon_id,
        "taxon_id",
        &primary_context_id,
        &fallback_context_id,
    )?;
    merge_field(
        &mut primary.genome_id,
        &fallback.genome_id,
        "genome_id",
        &primary_context_id,
        &fallback_context_id,
    )?;
    merge_field(
        &mut primary.assembly_accession,
        &fallback.assembly_accession,
        "assembly_accession",
        &primary_context_id,
        &fallback_context_id,
    )?;
    merge_field(
        &mut primary.annotation_source,
        &fallback.annotation_source,
        "annotation_source",
        &primary_context_id,
        &fallback_context_id,
    )?;
    merge_field(
        &mut primary.annotation_release,
        &fallback.annotation_release,
        "annotation_release",
        &primary_context_id,
        &fallback_context_id,
    )?;
    merge_field(
        &mut primary.symbol_namespace,
        &fallback.symbol_namespace,
        "symbol_namespace",
        &primary_context_id,
        &fallback_context_id,
    )
}

fn merge_field(
    primary: &mut Option<String>,
    fallback: &Option<String>,
    field: &str,
    primary_context_id: &str,
    fallback_context_id: &str,
) -> Result<(), BiologicalContextResolutionError> {
    match (primary.as_ref(), fallback.as_ref()) {
        (Some(primary_value), Some(fallback_value))
            if primary_value.trim() != fallback_value.trim() =>
        {
            Err(BiologicalContextResolutionError::ConflictingField {
                field: field.to_string(),
                primary_context_id: primary_context_id.to_string(),
                primary_value: primary_value.clone(),
                fallback_context_id: fallback_context_id.to_string(),
                fallback_value: fallback_value.clone(),
            })
        }
        (None, Some(fallback_value)) => {
            *primary = Some(fallback_value.clone());
            Ok(())
        }
        _ => Ok(()),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn member_context_inherits_compatible_default_and_legacy_fields() {
        let registry = BiologicalContextRegistry {
            contexts: vec![
                BiologicalContext {
                    context_id: "human".to_string(),
                    taxon_id: Some("9606".to_string()),
                    genome_id: Some("GRCh38".to_string()),
                    ..BiologicalContext::default()
                },
                BiologicalContext {
                    context_id: "member".to_string(),
                    organism: Some("Homo sapiens".to_string()),
                    ..BiologicalContext::default()
                },
            ],
            default_context_id: Some("human".to_string()),
        };
        let legacy = BiologicalContext {
            context_id: LEGACY_BIOLOGICAL_CONTEXT_ID.to_string(),
            annotation_release: Some("Ensembl 116".to_string()),
            ..BiologicalContext::default()
        };

        let resolved = registry
            .resolve(Some("member"), Some(&legacy))
            .expect("resolve compatible context")
            .expect("resolved context");
        assert_eq!(resolved.organism.as_deref(), Some("Homo sapiens"));
        assert_eq!(resolved.taxon_id.as_deref(), Some("9606"));
        assert_eq!(resolved.genome_id.as_deref(), Some("GRCh38"));
        assert_eq!(resolved.annotation_release.as_deref(), Some("Ensembl 116"));
    }

    #[test]
    fn conflicting_context_sources_fail_instead_of_silently_overriding() {
        let registry = BiologicalContextRegistry {
            contexts: vec![BiologicalContext {
                context_id: "mouse".to_string(),
                taxon_id: Some("10090".to_string()),
                ..BiologicalContext::default()
            }],
            default_context_id: Some("mouse".to_string()),
        };
        let legacy = BiologicalContext {
            context_id: LEGACY_BIOLOGICAL_CONTEXT_ID.to_string(),
            taxon_id: Some("9606".to_string()),
            ..BiologicalContext::default()
        };

        assert!(matches!(
            registry.resolve(None, Some(&legacy)),
            Err(BiologicalContextResolutionError::ConflictingField {
                field,
                ..
            }) if field == "taxon_id"
        ));
    }
}
