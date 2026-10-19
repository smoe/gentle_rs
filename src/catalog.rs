//! Shared catalog inspection/search and explicitly approved provider-only updates.
//!
//! These routes reuse existing overlays. They neither join biological identity
//! to sellers nor expose private planning profiles or execute vendor protocols.

use crate::engine::GentleEngine;
use crate::external_service_providers::{doctor_external_service_provider_config, updates};
use serde_json::{Value, json};
use std::path::Path;

/// Shared-shell catalog actions. Both artifact-writing actions require confirmation.
#[derive(Debug, Clone)]
pub enum CatalogCommand {
    Doctor,
    Search {
        entity: String,
        query: String,
        limit: usize,
    },
    Preview {
        bundle_path: String,
        scope: updates::UpdateScope,
        output: String,
    },
    Apply {
        preview_path: String,
        approval: String,
    },
}

impl CatalogCommand {
    pub fn mutating(&self) -> bool {
        matches!(self, Self::Preview { .. } | Self::Apply { .. })
    }
}

/// Parse exact option sets; never infer destinations, providers or approvals.
pub fn parse(tokens: &[String]) -> Result<CatalogCommand, String> {
    let args: Vec<&str> = tokens.iter().map(String::as_str).collect();
    match args.as_slice() {
        ["catalog", "doctor"] => Ok(CatalogCommand::Doctor),
        ["catalog", "search", rest @ ..] => {
            let mut entity = None; let mut query = None; let mut limit = None;
            for pair in rest.chunks(2) {
                let [flag, value] = pair else { return Err("catalog: option requires a value".into()); };
                match *flag {
                    "--entity" if entity.is_none() => entity = Some((*value).to_string()),
                    "--query" if query.is_none() => query = Some((*value).to_string()),
                    "--limit" if limit.is_none() => limit = Some(value.parse::<usize>().map_err(|_| "catalog: invalid limit")?),
                    _ => return Err("catalog: unknown or duplicate option".into()),
                }
            }
            let entity = entity.ok_or("catalog: --entity required")?;
            if !["product", "provider", "protocol"].contains(&entity.as_str()) { return Err("catalog: entity must be product, provider or protocol".into()); }
            let query = query.ok_or("catalog: --query required")?;
            let limit = limit.unwrap_or(50);
            if query.trim().is_empty() || query.len() > 256 || !(1..=100).contains(&limit) { return Err("catalog: nonempty query <=256 bytes and limit 1..100 required".into()); }
            Ok(CatalogCommand::Search { entity, query, limit })
        }
        ["catalog", "updates", "preview", bundle, rest @ ..] => {
            let mut scope = None; let mut output = None;
            for pair in rest.chunks(2) {
                let [flag, value] = pair else { return Err("catalog: option requires a value".into()); };
                match *flag {
                    "--scope" if scope.is_none() => scope = Some(match *value {
                        "user" => updates::UpdateScope::User, "project" => updates::UpdateScope::Project,
                        _ => return Err("catalog: scope must be user or project".into()),
                    }),
                    "--output" if output.is_none() => output = Some((*value).to_string()),
                    _ => return Err("catalog: unknown or duplicate option".into()),
                }
            }
            Ok(CatalogCommand::Preview { bundle_path: (*bundle).into(), scope: scope.ok_or("catalog: --scope required")?, output: output.ok_or("catalog: --output required")? })
        }
        ["catalog", "updates", "apply", preview, "--approve-digest", approval] => Ok(CatalogCommand::Apply { preview_path: (*preview).into(), approval: (*approval).into() }),
        _ => Err("catalog: expected doctor, search --entity KIND --query TEXT, updates preview BUNDLE --scope user|project --output FILE, or updates apply PREVIEW --approve-digest SHA256".into()),
    }
}

fn lanes() -> Value {
    json!({
        "product_identity": "existing_helper_catalog_identity_only",
        "commercial_offer": "unavailable_no_joined_offer_contract",
        "protocol_document": "linked_references_only_not_retained_documents",
        "normalized_procedure": "unavailable_no_reviewed_vendor_procedure",
        "local_profile_and_observations": "private_not_inspected"
    })
}

fn search(entity: &str, query: &str, limit: usize) -> Result<Value, String> {
    if query.trim().is_empty() || query.len() > 256 || !(1..=100).contains(&limit) {
        return Err("catalog: nonempty query <=256 bytes and limit 1..100 required".into());
    }
    let candidates = match entity {
        "product" => {
            GentleEngine::list_helper_vector_cards(None, None).map_err(|e| e.to_string())?.cards.into_iter().map(|card| {
                let bytes = crate::digest_utils::canonical_json_string(&card).map_err(|e| e.to_string())?;
                Ok(json!({ "id": card.helper_id, "identity_kind": "helper_catalog_identity", "record": card, "record_sha256": crate::digest_utils::sha256_prefixed_str(&bytes), "commercial_offer_status": "not_evaluated" }))
            }).collect::<Result<Vec<_>,String>>()?
        }
        "provider" => updates::inspection_rows()?,
        "protocol" => {
            let mut rows = Vec::new();
            for provider in updates::inspection_rows()? {
                if let Some(protocols) = provider["record"]["protocols"].as_array() {
                    for (ordinal, protocol) in protocols.iter().enumerate() {
                        rows.push(json!({ "id": format!("{}:protocol:{}", provider["id"].as_str().unwrap_or_default(), ordinal), "provider": provider["id"], "country_code": provider["record"]["country_code"], "reference": protocol, "source_sha256": provider["source_sha256"], "document_sha256": null, "procedure_review": "not_evaluated" }));
                    }
                }
            }
            rows
        }
        _ => return Err("catalog: unsupported entity".into()),
    };
    // Literal case-insensitive substring, not biological identity resolution.
    let query_lower = query.to_lowercase();
    let mut matches = candidates
        .into_iter()
        .filter(|row| {
            let searchable = row
                .get("record")
                .or_else(|| row.get("reference"))
                .unwrap_or(row);
            searchable.to_string().to_lowercase().contains(&query_lower)
        })
        .collect::<Vec<_>>();
    matches.sort_by(|a, b| a["id"].as_str().cmp(&b["id"].as_str()));
    let total = matches.len();
    matches.truncate(limit);
    Ok(
        json!({ "schema": "gentle.catalog_search.v1", "entity": entity, "query": query, "match_rule": "literal_case_insensitive_substring", "total_matches": total, "truncated": total > limit, "rows": matches, "lanes": lanes(), "nonclaims": ["Helper rows include reference/host entries, not only orderable products. No seller, sequence, protocol applicability, availability or price is inferred."] }),
    )
}

/// One shared executor for CLI, GUI Shell and approved agent-plan pass-through.
pub fn execute(command: &CatalogCommand) -> Result<Value, String> {
    match command {
        CatalogCommand::Doctor => {
            let helper = match GentleEngine::doctor_helper_vector_catalog(None) {
                Ok(report) => serde_json::to_value(report).map_err(|e| e.to_string())?,
                Err(error) => json!({ "status": "unavailable", "error": error.to_string() }),
            };
            Ok(
                json!({ "schema": "gentle.catalog_doctor.v1", "providers": doctor_external_service_provider_config(None), "helper_identities": helper, "lanes": lanes(), "network_access": false, "private_profiles_inspected": false, "overall_acceptance": "not_evaluated" }),
            )
        }
        CatalogCommand::Search {
            entity,
            query,
            limit,
        } => search(entity, query, *limit),
        CatalogCommand::Preview {
            bundle_path,
            scope,
            output,
        } => {
            let (_, bytes) =
                updates::read_json::<updates::ProviderUpdateBundle>(Path::new(bundle_path))?;
            let report = updates::preview(&bytes, *scope)?;
            updates::write_preview(Path::new(output), &report)?;
            serde_json::to_value(report).map_err(|e| e.to_string())
        }
        CatalogCommand::Apply {
            preview_path,
            approval,
        } => {
            let (preview, _) = updates::read_json(Path::new(preview_path))?;
            serde_json::to_value(updates::apply(&preview, approval)?).map_err(|e| e.to_string())
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn catalog_parser_guards_scope_and_mutations() {
        for (line, mutation) in [
            ("catalog doctor", false),
            ("catalog search --entity product --query luciferas", false),
            (
                "catalog updates preview 'bundle with spaces.json' --scope project --output preview.json",
                true,
            ),
            (
                "catalog updates apply preview.json --approve-digest sha256:example",
                true,
            ),
        ] {
            let parsed = crate::engine_shell::parse_shell_line(line).unwrap();
            assert_eq!(parsed.is_state_mutating(), mutation);
        }
        for line in [
            "catalog updates preview x --scope system --output p",
            "catalog search --entity offer --query a",
            "catalog search --entity product --query a --limit 0",
            "catalog search --entity product --query a --limit 101",
            "catalog doctor --output x",
        ] {
            assert!(
                crate::engine_shell::parse_shell_line(line).is_err(),
                "{line}"
            );
        }
    }
    #[test]
    fn catalog_doctor_is_read_only_and_missing_lanes_are_not_zero() {
        let mut engine = GentleEngine::from_state(Default::default());
        let before = serde_json::to_value(engine.state()).unwrap();
        let command = crate::engine_shell::parse_shell_line("catalog doctor").unwrap();
        let result = crate::engine_shell::execute_shell_command(&mut engine, &command).unwrap();
        assert!(!result.state_changed);
        assert_eq!(
            result.output["lanes"]["local_profile_and_observations"],
            "private_not_inspected"
        );
        assert_eq!(result.output["overall_acceptance"], "not_evaluated");
        assert_eq!(serde_json::to_value(engine.state()).unwrap(), before);
    }
    #[test]
    fn catalog_import_plan_cannot_hide_mutation_or_bypass_confirmation() {
        let mut engine = GentleEngine::from_state(Default::default());
        for command in [
            "catalog updates preview missing.json --scope project --output preview.json",
            "catalog updates apply missing.json --approve-digest sha256:example",
        ] {
            let plan = serde_json::from_value(json!({
                "schema":"gentle.agent_plan_result.v1", "assistant_message":"Synthetic approval test", "questions":[],
                "candidates":[{"candidate_id":"update", "title":"Update", "rationale":"Synthetic", "kind":"shell", "mutating":false, "requires_confirmation":false, "execution_mode":"auto", "shell_command":command}]
            })).unwrap();
            let error = crate::agent_execution::execute_agent_plan_candidate(
                &mut engine,
                &plan,
                "update",
                false,
                &Default::default(),
            )
            .unwrap_err();
            assert!(error.contains("explicit --confirm"), "{error}");
        }
    }
}
