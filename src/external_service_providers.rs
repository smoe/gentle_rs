//! Overlay-discoverable external-service provider configuration.
//!
//! Provider-specific order-channel behavior lives in JSON catalogs so GUI,
//! CLI, MCP, and ClawBio keep one provider-neutral request contract while
//! operators can update vendor handoff details without recompiling GENtle.

use gentle_protocol::{
    EXTERNAL_SERVICE_PROVIDER_CATALOG_SCHEMA, EXTERNAL_SERVICE_PROVIDER_CONFIG_DOCTOR_SCHEMA,
    EXTERNAL_SERVICE_PROVIDER_CONFIG_SCHEMA, ExternalServiceCapability,
    ExternalServiceChannelConfig, ExternalServiceProductTemplateMapping,
    ExternalServiceProviderCatalog, ExternalServiceProviderConfigCatalog,
    ExternalServiceProviderConfigDoctorReport, ExternalServiceProviderConfigRecord,
    ExternalServiceProviderConfigSourceReport, ExternalServiceProviderRecord,
    ExternalServiceValidationRule,
};
use sha1::{Digest, Sha1};
use std::{
    collections::{BTreeMap, BTreeSet},
    fs,
    path::{Path, PathBuf},
};

pub const DEFAULT_EXTERNAL_SERVICE_PROVIDER_CONFIG_PATH: &str =
    "assets/external_service_providers.json";
pub const DEFAULT_EXTERNAL_SERVICE_PROVIDER_CONFIG_DISCOVERY_TOKEN: &str =
    "gentle://catalog/external-service-providers/default";

#[derive(Clone, Debug)]
struct ExternalServiceProviderConfigSourceCandidate {
    scope: &'static str,
    path: PathBuf,
}

#[derive(Clone, Debug)]
pub struct LoadedExternalServiceProviderConfig {
    pub record: ExternalServiceProviderConfigRecord,
    pub source_scope: String,
    pub source_path: String,
}

#[derive(Clone, Debug)]
pub struct ExternalServiceProviderConfigIndex {
    catalog_label: String,
    providers: Vec<LoadedExternalServiceProviderConfig>,
    summary_lines: Vec<String>,
    warnings: Vec<String>,
}

pub fn default_external_service_provider_config_discovery_label() -> &'static str {
    DEFAULT_EXTERNAL_SERVICE_PROVIDER_CONFIG_DISCOVERY_TOKEN
}

fn normalize_service_token(raw: &str) -> String {
    raw.trim()
        .to_ascii_lowercase()
        .replace([' ', '-'], "_")
        .split('_')
        .filter(|part| !part.is_empty())
        .collect::<Vec<_>>()
        .join("_")
}

fn external_service_provider_config_discovery_candidates()
-> Vec<ExternalServiceProviderConfigSourceCandidate> {
    let mut candidates = vec![];
    let built_in_root = crate::genomes::configured_builtin_asset_root();
    candidates.push(ExternalServiceProviderConfigSourceCandidate {
        scope: "built-in",
        path: built_in_root.join(DEFAULT_EXTERNAL_SERVICE_PROVIDER_CONFIG_PATH),
    });
    candidates.push(ExternalServiceProviderConfigSourceCandidate {
        scope: "built-in",
        path: built_in_root
            .join("assets")
            .join("external_service_providers.d"),
    });

    let system_root = crate::genomes::configured_system_config_root();
    candidates.push(ExternalServiceProviderConfigSourceCandidate {
        scope: "system",
        path: system_root
            .join("catalogs")
            .join("external_service_providers.json"),
    });
    candidates.push(ExternalServiceProviderConfigSourceCandidate {
        scope: "system",
        path: system_root
            .join("catalogs")
            .join("external_service_providers.d"),
    });

    if let Some(user_root) = crate::genomes::configured_user_config_root() {
        candidates.push(ExternalServiceProviderConfigSourceCandidate {
            scope: "user",
            path: user_root
                .join("catalogs")
                .join("external_service_providers.json"),
        });
        candidates.push(ExternalServiceProviderConfigSourceCandidate {
            scope: "user",
            path: user_root
                .join("catalogs")
                .join("external_service_providers.d"),
        });
    }

    if let Some(project_root) = crate::genomes::configured_project_root() {
        candidates.push(ExternalServiceProviderConfigSourceCandidate {
            scope: "project",
            path: project_root
                .join(".gentle")
                .join("catalogs")
                .join("external_service_providers.json"),
        });
        candidates.push(ExternalServiceProviderConfigSourceCandidate {
            scope: "project",
            path: project_root
                .join(".gentle")
                .join("catalogs")
                .join("external_service_providers.d"),
        });
    }

    let mut seen = BTreeSet::<String>::new();
    let mut out = vec![];
    for candidate in candidates {
        let key = candidate.path.to_string_lossy().to_string();
        if seen.insert(key) {
            out.push(candidate);
        }
    }
    out
}

fn discovered_external_service_provider_config_sources()
-> Vec<ExternalServiceProviderConfigSourceCandidate> {
    external_service_provider_config_discovery_candidates()
        .into_iter()
        .filter(|candidate| fs::metadata(&candidate.path).is_ok())
        .collect()
}

fn json_files_in_dir(path: &Path) -> Result<Vec<PathBuf>, String> {
    let mut files = vec![];
    for entry in fs::read_dir(path).map_err(|e| {
        format!(
            "Could not read provider config directory '{}': {e}",
            path.display()
        )
    })? {
        let entry = entry.map_err(|e| {
            format!(
                "Could not read provider config directory entry in '{}': {e}",
                path.display()
            )
        })?;
        let entry_path = entry.path();
        if entry_path
            .extension()
            .and_then(|ext| ext.to_str())
            .is_some_and(|ext| ext.eq_ignore_ascii_case("json"))
        {
            files.push(entry_path);
        }
    }
    files.sort();
    Ok(files)
}

fn source_sha1(path: &Path) -> Option<String> {
    fs::read(path)
        .ok()
        .map(|bytes| format!("{:x}", Sha1::digest(&bytes)))
}

fn validate_provider_config_catalog(
    catalog: &ExternalServiceProviderConfigCatalog,
    path: &Path,
) -> Vec<String> {
    let mut errors = vec![];
    if catalog.schema != EXTERNAL_SERVICE_PROVIDER_CONFIG_SCHEMA {
        errors.push(format!(
            "Provider config '{}' declares schema '{}' but expected '{}'",
            path.display(),
            catalog.schema,
            EXTERNAL_SERVICE_PROVIDER_CONFIG_SCHEMA
        ));
    }
    let mut seen = BTreeSet::<String>::new();
    for provider in &catalog.providers {
        let provider_id = normalize_service_token(&provider.provider);
        if provider_id.is_empty() {
            errors.push(format!(
                "Provider config '{}' contains a provider with empty provider id",
                path.display()
            ));
            continue;
        }
        if !seen.insert(provider_id.clone()) {
            errors.push(format!(
                "Provider config '{}' declares provider '{}' more than once",
                path.display(),
                provider.provider
            ));
        }
        if provider.display_name.trim().is_empty() {
            errors.push(format!(
                "Provider '{}' in '{}' requires display_name",
                provider.provider,
                path.display()
            ));
        }
        if provider.website_url.trim().is_empty() {
            errors.push(format!(
                "Provider '{}' in '{}' requires website_url",
                provider.provider,
                path.display()
            ));
        }
        if provider.dashboard_url.trim().is_empty() {
            errors.push(format!(
                "Provider '{}' in '{}' requires dashboard_url",
                provider.provider,
                path.display()
            ));
        }
        if provider.capabilities.is_empty() {
            errors.push(format!(
                "Provider '{}' in '{}' requires at least one capability",
                provider.provider,
                path.display()
            ));
        }
        for capability in &provider.capabilities {
            if capability.service_kind.trim().is_empty() {
                errors.push(format!(
                    "Provider '{}' in '{}' has a capability with empty service_kind",
                    provider.provider,
                    path.display()
                ));
            }
        }
    }
    errors
}

impl LoadedExternalServiceProviderConfig {
    pub fn provider_record(&self) -> ExternalServiceProviderRecord {
        ExternalServiceProviderRecord {
            provider: normalize_service_token(&self.record.provider),
            display_name: self.record.display_name.clone(),
            support_status: self.record.support_status.clone(),
            website_url: self.record.website_url.clone(),
            dashboard_url: self.record.dashboard_url.clone(),
            api_documentation_url: self.record.api_documentation_url.clone(),
            capabilities: self.record.capabilities.clone(),
            account_enablement_notes: self.record.account_enablement_notes.clone(),
            warnings: self.record.warnings.clone(),
        }
    }

    pub fn capability_for(&self, service_kind: &str) -> Option<ExternalServiceCapability> {
        let normalized = normalize_service_token(service_kind);
        self.record
            .capabilities
            .iter()
            .find(|capability| normalize_service_token(&capability.service_kind) == normalized)
            .cloned()
    }

    pub fn validation_rule_for(&self, service_kind: &str) -> Option<ExternalServiceValidationRule> {
        let normalized = normalize_service_token(service_kind);
        self.record
            .validation_rules
            .iter()
            .find(|rule| normalize_service_token(&rule.service_kind) == normalized)
            .cloned()
    }

    pub fn product_template_for(
        &self,
        service_kind: &str,
    ) -> Option<ExternalServiceProductTemplateMapping> {
        let normalized = normalize_service_token(service_kind);
        self.record
            .product_templates
            .iter()
            .find(|mapping| normalize_service_token(&mapping.service_kind) == normalized)
            .cloned()
    }

    pub fn channel_for(&self, channel: &str) -> Option<ExternalServiceChannelConfig> {
        let normalized = normalize_service_token(channel);
        self.record
            .channels
            .iter()
            .find(|row| normalize_service_token(&row.channel) == normalized)
            .cloned()
    }
}

impl ExternalServiceProviderConfigIndex {
    pub fn from_default_discovery() -> Result<Self, String> {
        let sources = discovered_external_service_provider_config_sources();
        if sources.is_empty() {
            return Err(format!(
                "No external-service provider config sources were found for {}",
                default_external_service_provider_config_discovery_label()
            ));
        }
        Self::from_sources(
            &sources,
            default_external_service_provider_config_discovery_label().to_string(),
        )
    }

    pub fn from_explicit_path(path: &str) -> Result<Self, String> {
        Self::from_sources(
            &[ExternalServiceProviderConfigSourceCandidate {
                scope: "explicit",
                path: PathBuf::from(path),
            }],
            path.to_string(),
        )
    }

    fn from_sources(
        sources: &[ExternalServiceProviderConfigSourceCandidate],
        catalog_label: String,
    ) -> Result<Self, String> {
        let mut provider_by_id = BTreeMap::<String, LoadedExternalServiceProviderConfig>::new();
        let mut summary_lines = vec![];
        let mut warnings = vec![];
        for source in sources {
            merge_config_source(
                source,
                &mut provider_by_id,
                &mut summary_lines,
                &mut warnings,
            )?;
        }
        Ok(Self {
            catalog_label,
            providers: provider_by_id.into_values().collect(),
            summary_lines,
            warnings,
        })
    }

    pub fn provider(&self, provider: &str) -> Option<&LoadedExternalServiceProviderConfig> {
        let normalized = normalize_service_token(provider);
        self.providers
            .iter()
            .find(|row| normalize_service_token(&row.record.provider) == normalized)
    }

    pub fn provider_ids(&self) -> Vec<String> {
        self.providers
            .iter()
            .map(|row| normalize_service_token(&row.record.provider))
            .collect()
    }

    pub fn to_provider_catalog(
        &self,
        generated_at_unix_ms: u128,
    ) -> ExternalServiceProviderCatalog {
        let mut summary_lines = self.summary_lines.clone();
        summary_lines.extend(
            self.warnings
                .iter()
                .map(|warning| format!("Warning: {warning}")),
        );
        summary_lines.push(format!(
            "External-service provider behavior loaded via {}.",
            self.catalog_label
        ));
        ExternalServiceProviderCatalog {
            schema: EXTERNAL_SERVICE_PROVIDER_CATALOG_SCHEMA.to_string(),
            generated_at_unix_ms,
            providers: self
                .providers
                .iter()
                .map(LoadedExternalServiceProviderConfig::provider_record)
                .collect(),
            summary_lines,
        }
    }
}

fn merge_config_source(
    source: &ExternalServiceProviderConfigSourceCandidate,
    provider_by_id: &mut BTreeMap<String, LoadedExternalServiceProviderConfig>,
    summary_lines: &mut Vec<String>,
    warnings: &mut Vec<String>,
) -> Result<(), String> {
    let metadata = fs::metadata(&source.path).map_err(|e| {
        format!(
            "Could not read external-service provider config '{}': {e}",
            source.path.display()
        )
    })?;
    if metadata.is_dir() {
        let json_files = json_files_in_dir(&source.path)?;
        for file_path in json_files {
            merge_config_file(
                source.scope,
                &file_path,
                provider_by_id,
                summary_lines,
                warnings,
            )?;
        }
        return Ok(());
    }
    merge_config_file(
        source.scope,
        &source.path,
        provider_by_id,
        summary_lines,
        warnings,
    )
}

fn merge_config_file(
    scope: &str,
    path: &Path,
    provider_by_id: &mut BTreeMap<String, LoadedExternalServiceProviderConfig>,
    summary_lines: &mut Vec<String>,
    warnings: &mut Vec<String>,
) -> Result<(), String> {
    let text = fs::read_to_string(path).map_err(|e| {
        format!(
            "Could not read external-service provider config '{}': {e}",
            path.display()
        )
    })?;
    let catalog: ExternalServiceProviderConfigCatalog =
        serde_json::from_str(&text).map_err(|e| {
            format!(
                "Could not parse external-service provider config '{}': {e}",
                path.display()
            )
        })?;
    let validation_errors = validate_provider_config_catalog(&catalog, path);
    if !validation_errors.is_empty() {
        return Err(validation_errors.join("; "));
    }
    summary_lines.extend(catalog.summary_lines);
    for provider in catalog.providers {
        let provider_id = normalize_service_token(&provider.provider);
        if let Some(previous) = provider_by_id.get(&provider_id) {
            warnings.push(format!(
                "Provider '{}' from '{}' overrides {} source '{}'",
                provider.provider,
                path.display(),
                previous.source_scope,
                previous.source_path
            ));
        }
        provider_by_id.insert(
            provider_id,
            LoadedExternalServiceProviderConfig {
                record: provider,
                source_scope: scope.to_string(),
                source_path: path.display().to_string(),
            },
        );
    }
    Ok(())
}

pub fn external_service_provider_config_index() -> Result<ExternalServiceProviderConfigIndex, String>
{
    ExternalServiceProviderConfigIndex::from_default_discovery()
}

pub fn doctor_external_service_provider_config(
    explicit_path: Option<&str>,
) -> ExternalServiceProviderConfigDoctorReport {
    let (sources, catalog_label) = if let Some(path) = explicit_path {
        (
            vec![ExternalServiceProviderConfigSourceCandidate {
                scope: "explicit",
                path: PathBuf::from(path),
            }],
            path.to_string(),
        )
    } else {
        (
            external_service_provider_config_discovery_candidates(),
            default_external_service_provider_config_discovery_label().to_string(),
        )
    };
    let mut report = ExternalServiceProviderConfigDoctorReport {
        schema: EXTERNAL_SERVICE_PROVIDER_CONFIG_DOCTOR_SCHEMA.to_string(),
        catalog_label,
        source_count: sources.len(),
        ..Default::default()
    };
    let mut provider_ids = BTreeSet::<String>::new();
    for source in sources {
        if fs::metadata(&source.path).is_err() {
            report
                .sources
                .push(ExternalServiceProviderConfigSourceReport {
                    scope: source.scope.to_string(),
                    path: source.path.display().to_string(),
                    status: "missing".to_string(),
                    ..Default::default()
                });
            continue;
        }
        let paths = if source.path.is_dir() {
            match json_files_in_dir(&source.path) {
                Ok(paths) => paths,
                Err(err) => {
                    report
                        .sources
                        .push(ExternalServiceProviderConfigSourceReport {
                            scope: source.scope.to_string(),
                            path: source.path.display().to_string(),
                            status: "error".to_string(),
                            errors: vec![err],
                            ..Default::default()
                        });
                    continue;
                }
            }
        } else {
            vec![source.path.clone()]
        };
        for path in paths {
            let text = match fs::read_to_string(&path) {
                Ok(text) => text,
                Err(err) => {
                    report
                        .sources
                        .push(ExternalServiceProviderConfigSourceReport {
                            scope: source.scope.to_string(),
                            path: path.display().to_string(),
                            status: "error".to_string(),
                            sha1: source_sha1(&path),
                            errors: vec![format!(
                                "Could not read external-service provider config '{}': {err}",
                                path.display()
                            )],
                            ..Default::default()
                        });
                    continue;
                }
            };
            match serde_json::from_str::<ExternalServiceProviderConfigCatalog>(&text) {
                Ok(catalog) => {
                    let errors = validate_provider_config_catalog(&catalog, &path);
                    let status = if errors.is_empty() { "ok" } else { "error" };
                    if errors.is_empty() {
                        report.parsed_source_count += 1;
                        for provider in &catalog.providers {
                            provider_ids.insert(normalize_service_token(&provider.provider));
                        }
                    }
                    report
                        .sources
                        .push(ExternalServiceProviderConfigSourceReport {
                            scope: source.scope.to_string(),
                            path: path.display().to_string(),
                            status: status.to_string(),
                            provider_count: catalog.providers.len(),
                            sha1: source_sha1(&path),
                            errors,
                            ..Default::default()
                        });
                }
                Err(err) => report
                    .sources
                    .push(ExternalServiceProviderConfigSourceReport {
                        scope: source.scope.to_string(),
                        path: path.display().to_string(),
                        status: "error".to_string(),
                        sha1: source_sha1(&path),
                        errors: vec![format!(
                            "Could not parse external-service provider config '{}': {err}",
                            path.display()
                        )],
                        ..Default::default()
                    }),
            }
        }
    }
    report.provider_count = provider_ids.len();
    report.warnings = report
        .sources
        .iter()
        .flat_map(|source| source.warnings.clone())
        .collect();
    report.errors = report
        .sources
        .iter()
        .flat_map(|source| source.errors.clone())
        .collect();
    report.warning_count = report.warnings.len();
    report.error_count = report.errors.len();
    report
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::tempdir;

    fn write_config(path: &Path, provider: &str, display_name: &str) {
        let text = format!(
            r#"{{
  "schema": "gentle.external_service_provider_config.v1",
  "providers": [
    {{
      "provider": "{provider}",
      "display_name": "{display_name}",
      "support_status": "quote_handoff_ready",
      "website_url": "https://example.invalid/{provider}",
      "dashboard_url": "https://example.invalid/{provider}/portal",
      "capabilities": [
        {{
          "service_kind": "dna_oligo_single_tube",
          "track": "dna",
          "display_name": "DNA oligo",
          "quote_handoff_supported": true,
          "direct_api_documented": false,
          "direct_api_implemented": false,
          "supported_submission_modes": ["quote_handoff"],
          "status_tracking": "manual_or_imported_status",
          "artifact_kinds": ["quote_metadata"]
        }}
      ]
    }}
  ]
}}"#
        );
        fs::write(path, text).expect("write config");
    }

    #[test]
    fn provider_config_overlay_later_source_overrides_provider() {
        let dir = tempdir().expect("tempdir");
        let built_in = dir.path().join("built_in.json");
        let project = dir.path().join("project.json");
        write_config(&built_in, "metabion", "metabion built-in");
        write_config(&project, "metabion", "metabion project");
        let index = ExternalServiceProviderConfigIndex::from_sources(
            &[
                ExternalServiceProviderConfigSourceCandidate {
                    scope: "built-in",
                    path: built_in,
                },
                ExternalServiceProviderConfigSourceCandidate {
                    scope: "project",
                    path: project,
                },
            ],
            "test".to_string(),
        )
        .expect("provider config index");
        let provider = index.provider("metabion").expect("metabion");
        assert_eq!(provider.record.display_name, "metabion project");
        assert_eq!(provider.source_scope, "project");
        assert!(
            index
                .warnings
                .iter()
                .any(|warning| warning.contains("overrides"))
        );
    }

    #[test]
    fn provider_config_doctor_reports_schema_errors() {
        let dir = tempdir().expect("tempdir");
        let path = dir.path().join("bad.json");
        fs::write(
            &path,
            r#"{"schema":"wrong","providers":[{"provider":"x"}]}"#,
        )
        .expect("write bad config");
        let report = doctor_external_service_provider_config(Some(&path.to_string_lossy()));
        assert_eq!(report.error_count, 5);
        assert_eq!(report.sources[0].status, "error");
    }
}
