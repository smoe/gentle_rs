//! Local provider overlay updates. Preview is evidence, never vendor authentication.
//!
//! One atomic catalog replacement is the commit point. Before it, an immutable
//! transaction directory retains the old bytes, proposed bytes, preview and
//! conditional receipt. A crash is resolved by comparing the target hash with
//! the receipt; history is never replayed or rolled back automatically.

use super::*;
use crate::digest_utils::{canonical_json_string, sha256_prefixed_bytes};
use serde::{Deserialize, Serialize};
use serde_json::Value;
use std::io::Write;

const MAX_BYTES: u64 = 4 * 1024 * 1024;
const MAX_SOURCES: usize = 128;

fn error(code: &str, detail: impl std::fmt::Display) -> String {
    format!("catalog_update:{code}: {detail}")
}

/// Only local user/project overlays are writable; built-in/system sources are read-only.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum UpdateScope {
    User,
    Project,
}

/// Full replacement or explicit tombstone. Patches are deliberately not admitted.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "action", rename_all = "snake_case", deny_unknown_fields)]
pub enum ProviderChange {
    Replace {
        provider: String,
        base_record_sha256: Option<String>,
        record: Value,
    },
    Withdraw {
        provider: String,
        base_record_sha256: String,
    },
}

/// A retained, hash-bound bundle of at most 64 provider changes. No network reads.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ProviderUpdateBundle {
    pub schema: String,
    pub changes: Vec<ProviderChange>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct SourceBinding {
    pub scope: String,
    pub path: String,
    pub kind: String,
    pub sha256: Option<String>,
    pub files: Vec<SourceFileBinding>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct SourceFileBinding {
    pub path: String,
    pub sha256: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RecordDecision {
    pub provider: String,
    pub action: String,
    pub previous_record_sha256: Option<String>,
    pub proposed_record_sha256: Option<String>,
    pub diagnostics: Vec<String>,
}

/// All fields are recomputed before apply, including destination and overlay membership.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ProviderUpdatePreview {
    pub schema: String,
    pub scope: UpdateScope,
    pub target_path: String,
    pub history_directory: String,
    pub output_files: Vec<String>,
    pub bundle: ProviderUpdateBundle,
    pub bundle_source: String,
    pub bundle_bytes_sha256: String,
    pub sources: Vec<SourceBinding>,
    pub previous_target_sha256: Option<String>,
    pub proposed_target_sha256: String,
    pub proposed_catalog: ExternalServiceProviderConfigCatalog,
    pub decisions: Vec<RecordDecision>,
    pub applicable: bool,
    pub impacts: Vec<String>,
    pub warnings: Vec<String>,
    pub preview_sha256: String,
}

/// The target hash is the commit condition, not a multi-file transaction claim.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ProviderUpdateReceipt {
    pub schema: String,
    pub preview_sha256: String,
    pub target_path: String,
    pub previous_target_sha256: Option<String>,
    pub applied_target_sha256: String,
    pub history_directory: String,
    pub commit_condition: String,
}

#[derive(Clone)]
struct Context {
    scope: UpdateScope,
    target: PathBuf,
    candidates: Vec<ExternalServiceProviderConfigSourceCandidate>,
}

impl Context {
    fn discover(scope: UpdateScope) -> Result<Self, String> {
        let candidates = external_service_provider_config_discovery_candidates();
        let scope_name = match scope {
            UpdateScope::User => "user",
            UpdateScope::Project => "project",
        };
        let target = candidates
            .iter()
            .find(|source| source.scope == scope_name && source.path.extension().is_some())
            .ok_or_else(|| error("scope_unavailable", scope_name))?
            .path
            .clone();
        Ok(Self {
            scope,
            target: checked_path(&target)?,
            candidates,
        })
    }
}

/// Reject raw traversal and symlink ancestors before path normalization.
fn checked_path(path: &Path) -> Result<PathBuf, String> {
    if path
        .components()
        .any(|part| part == std::path::Component::ParentDir || part.as_os_str() == "..")
    {
        return Err(error("unsafe_path", path.display()));
    }
    let absolute = std::path::absolute(path).map_err(|e| error("path", e))?;
    let mut current = PathBuf::new();
    for component in absolute.components() {
        current.push(component.as_os_str());
        match fs::symlink_metadata(&current) {
            Ok(meta) if meta.file_type().is_symlink() => {
                return Err(error("symlink", current.display()));
            }
            Ok(_) => {}
            Err(e) if e.kind() == std::io::ErrorKind::NotFound => {}
            Err(e) => return Err(error("path", e)),
        }
    }
    Ok(absolute)
}

/// Parse and fingerprint exactly the same retained bytes, with bounded admission.
pub fn read_json<T: serde::de::DeserializeOwned>(path: &Path) -> Result<(T, Vec<u8>), String> {
    let bytes = read_bytes(path)?;
    let value = serde_json::from_slice(&bytes).map_err(|e| error("json", e))?;
    Ok((value, bytes))
}

fn read_bytes(path: &Path) -> Result<Vec<u8>, String> {
    use std::io::Read;
    let path = checked_path(path)?;
    let metadata = fs::metadata(&path).map_err(|e| error("read", e))?;
    if !metadata.is_file() {
        return Err(error("not_regular_file", path.display()));
    }
    let mut bytes = Vec::new();
    fs::File::open(&path)
        .map_err(|e| error("read", e))?
        .take(MAX_BYTES + 1)
        .read_to_end(&mut bytes)
        .map_err(|e| error("read", e))?;
    if bytes.len() as u64 > MAX_BYTES {
        return Err(error("byte_budget", path.display()));
    }
    Ok(bytes)
}

pub fn record_sha256(record: &ExternalServiceProviderConfigRecord) -> Result<String, String> {
    canonical_json_string(record)
        .map(|text| sha256_prefixed_bytes(text.as_bytes()))
        .map_err(|e| error("json", e))
}

fn catalog_bytes(catalog: &ExternalServiceProviderConfigCatalog) -> Result<Vec<u8>, String> {
    let mut bytes = serde_json::to_vec_pretty(catalog).map_err(|e| error("json", e))?;
    bytes.push(b'\n');
    if bytes.len() as u64 > MAX_BYTES {
        return Err(error("byte_budget", "proposed catalog"));
    }
    Ok(bytes)
}

struct Captured {
    bindings: Vec<SourceBinding>,
    files: Vec<(String, PathBuf, Vec<u8>)>,
}

fn capture(context: &Context) -> Result<Captured, String> {
    let mut bindings = Vec::new();
    let mut files = Vec::new();
    for source in &context.candidates {
        let path = checked_path(&source.path)?;
        let mut paths = Vec::new();
        let kind = match fs::metadata(&path) {
            Ok(meta) if meta.is_dir() => {
                paths = json_files_in_dir(&path)?;
                "directory"
            }
            Ok(meta) if meta.is_file() => {
                paths.push(path.clone());
                "file"
            }
            Ok(_) => return Err(error("source_kind", path.display())),
            Err(e) if e.kind() == std::io::ErrorKind::NotFound => "missing",
            Err(e) => return Err(error("source", e)),
        };
        let mut members = Vec::new();
        for file in paths {
            if files.len() >= MAX_SOURCES {
                return Err(error("source_budget", MAX_SOURCES));
            }
            let file = checked_path(&file)?;
            let bytes = read_bytes(&file)?;
            members.push(SourceFileBinding {
                path: file.display().to_string(),
                sha256: sha256_prefixed_bytes(&bytes),
            });
            files.push((source.scope.to_string(), file, bytes));
            if files.iter().map(|(_, _, b)| b.len()).sum::<usize>() > 32 * 1024 * 1024 {
                return Err(error("byte_budget", "combined catalog sources"));
            }
        }
        bindings.push(SourceBinding {
            scope: source.scope.to_string(),
            path: path.display().to_string(),
            kind: kind.into(),
            sha256: if kind == "missing" {
                None
            } else if kind == "file" {
                Some(members[0].sha256.clone())
            } else {
                Some(sha256_prefixed_bytes(
                    canonical_json_string(&members)
                        .map_err(|e| error("json", e))?
                        .as_bytes(),
                ))
            },
            files: members,
        });
    }
    Ok(Captured { bindings, files })
}

fn parse_catalog(
    bytes: &[u8],
    path: &Path,
) -> Result<ExternalServiceProviderConfigCatalog, String> {
    let original: Value = serde_json::from_slice(bytes).map_err(|e| error("json", e))?;
    let catalog: ExternalServiceProviderConfigCatalog =
        serde_json::from_value(original.clone()).map_err(|e| error("json", e))?;
    // Legacy defaults are allowed; unknown fields must not disappear during rewrite.
    fn retained(input: &Value, output: &Value) -> bool {
        match (input, output) {
            (Value::Object(a), Value::Object(b)) => a
                .iter()
                .all(|(k, v)| b.get(k).is_some_and(|o| retained(v, o))),
            (Value::Array(a), Value::Array(b)) => {
                a.len() == b.len() && a.iter().zip(b).all(|(a, b)| retained(a, b))
            }
            _ => input == output,
        }
    }
    let normalized = serde_json::to_value(&catalog).map_err(|e| error("json", e))?;
    // Empty optional fields may be omitted by serialization, so validate retained
    // keys after removing only these known, semantically empty fields.
    let mut original = original;
    if original.get("withdrawn_providers") == Some(&serde_json::json!([])) {
        original
            .as_object_mut()
            .unwrap()
            .remove("withdrawn_providers");
    }
    if let Some(rows) = original.get_mut("providers").and_then(Value::as_array_mut) {
        for row in rows {
            for key in ["country_code", "country_source_url", "protocols"] {
                if row
                    .get(key)
                    .is_some_and(|v| v.is_null() || v == &serde_json::json!([]))
                {
                    row.as_object_mut().unwrap().remove(key);
                }
            }
        }
    }
    if !retained(&original, &normalized) {
        return Err(error("unknown_fields", path.display()));
    }
    let errors = validate_provider_config_catalog(&catalog, path);
    if !errors.is_empty() {
        return Err(error("validation", errors.join("; ")));
    }
    Ok(catalog)
}

fn effective(
    files: &[(String, PathBuf, Vec<u8>)],
) -> Result<BTreeMap<String, (ExternalServiceProviderConfigRecord, PathBuf)>, String> {
    let mut rows = BTreeMap::new();
    for (_, path, bytes) in files {
        let catalog = parse_catalog(bytes, path)?;
        for id in catalog.withdrawn_providers {
            rows.remove(&normalize_service_token(&id));
        }
        for record in catalog.providers {
            rows.insert(
                normalize_service_token(&record.provider),
                (record, path.clone()),
            );
        }
    }
    Ok(rows)
}

/// Effective full records and exact source-byte fingerprints for inspection and bundle authorship.
pub fn inspection_rows() -> Result<Vec<Value>, String> {
    let context = Context {
        scope: UpdateScope::Project,
        target: PathBuf::new(),
        candidates: external_service_provider_config_discovery_candidates(),
    };
    let captured = capture(&context)?;
    if captured.files.is_empty() {
        return Err(error(
            "sources_unavailable",
            "no provider catalog files; not an empty provider inventory",
        ));
    }
    effective(&captured.files)?.into_iter().map(|(id, (record, path))| {
        let (scope, _, bytes) = captured.files.iter().find(|(_, p, _)| p == &path).ok_or_else(|| error("source", "missing binding"))?;
        Ok(serde_json::json!({ "id": id, "record": record, "record_sha256": record_sha256(&record)?, "source_scope": scope, "source_path": path, "source_sha256": sha256_prefixed_bytes(bytes) }))
    }).collect()
}

/// Preview never writes catalogs, contacts suppliers, fetches documents or reads lab profiles.
pub fn preview(bytes: &[u8], scope: UpdateScope) -> Result<ProviderUpdatePreview, String> {
    preview_in(bytes, &Context::discover(scope)?)
}

fn preview_in(bytes: &[u8], context: &Context) -> Result<ProviderUpdatePreview, String> {
    if bytes.len() > 512 * 1024 {
        return Err(error("byte_budget", "bundle"));
    }
    let bundle: ProviderUpdateBundle =
        serde_json::from_slice(bytes).map_err(|e| error("json", e))?;
    if bundle.schema != "gentle.provider_update_bundle.v1"
        || bundle.changes.is_empty()
        || bundle.changes.len() > 64
    {
        return Err(error("bundle_admission", "expected v1 and 1..64 changes"));
    }
    let captured = capture(context)?;
    let before = effective(&captured.files)?;
    let prior = captured
        .files
        .iter()
        .find(|(_, path, _)| path == &context.target)
        .map(|(_, _, b)| b.clone());
    let mut proposed = match &prior {
        Some(bytes) => parse_catalog(bytes, &context.target)?,
        None => ExternalServiceProviderConfigCatalog {
            schema: EXTERNAL_SERVICE_PROVIDER_CONFIG_SCHEMA.into(),
            ..Default::default()
        },
    };
    let mut decisions = Vec::new();
    let mut seen = BTreeSet::new();
    for change in &bundle.changes {
        let (id, base) = match change {
            ProviderChange::Replace {
                provider,
                base_record_sha256,
                ..
            } => (
                normalize_service_token(provider),
                base_record_sha256.as_ref(),
            ),
            ProviderChange::Withdraw {
                provider,
                base_record_sha256,
            } => (normalize_service_token(provider), Some(base_record_sha256)),
        };
        if id.is_empty() || !seen.insert(id.clone()) {
            return Err(error("duplicate_or_empty_id", id));
        }
        let previous = before.get(&id).map(|(r, _)| record_sha256(r)).transpose()?;
        let mut diagnostics = Vec::new();
        if base != previous.as_ref() {
            diagnostics.push("base_record_conflict".into());
        }
        let (action, next) = match change {
            ProviderChange::Replace { record, .. } => {
                let row: ExternalServiceProviderConfigRecord =
                    serde_json::from_value(record.clone()).map_err(|e| error("replacement", e))?;
                if normalize_service_token(&row.provider) != id
                    || serde_json::to_value(&row).map_err(|e| error("json", e))? != *record
                {
                    return Err(error("complete_record_required", id));
                }
                let check = ExternalServiceProviderConfigCatalog {
                    schema: EXTERNAL_SERVICE_PROVIDER_CONFIG_SCHEMA.into(),
                    providers: vec![row.clone()],
                    ..Default::default()
                };
                let errors = validate_provider_config_catalog(&check, &context.target);
                if !errors.is_empty() {
                    return Err(error("validation", errors.join("; ")));
                }
                for protocol in &row.protocols {
                    diagnostics.push(format!(
                        "linked_document_only:{}:checked_on={}:freshness_not_evaluated",
                        protocol.url, protocol.checked_on
                    ));
                }
                if row.protocols.is_empty() {
                    diagnostics.push("protocol_document_evidence_unavailable".into());
                }
                let next = record_sha256(&row)?;
                if previous.as_ref() != Some(&next) {
                    proposed
                        .providers
                        .retain(|r| normalize_service_token(&r.provider) != id);
                    proposed
                        .withdrawn_providers
                        .retain(|p| normalize_service_token(p) != id);
                    proposed.providers.push(row);
                }
                (
                    if previous.as_ref() == Some(&next) {
                        "no_op"
                    } else if previous.is_some() {
                        "replace"
                    } else {
                        "add"
                    },
                    Some(next),
                )
            }
            ProviderChange::Withdraw { .. } => {
                if previous.is_none() {
                    diagnostics.push("withdraw_requires_existing_record".into());
                }
                proposed
                    .providers
                    .retain(|r| normalize_service_token(&r.provider) != id);
                proposed
                    .withdrawn_providers
                    .retain(|p| normalize_service_token(p) != id);
                proposed.withdrawn_providers.push(id.clone());
                ("withdraw", None)
            }
        };
        decisions.push(RecordDecision {
            provider: id,
            action: action.into(),
            previous_record_sha256: previous,
            proposed_record_sha256: next,
            diagnostics,
        });
    }
    proposed
        .providers
        .sort_by_key(|r| normalize_service_token(&r.provider));
    proposed.withdrawn_providers.sort();
    if !proposed.withdrawn_providers.is_empty() {
        proposed.schema = EXTERNAL_SERVICE_PROVIDER_CONFIG_WITHDRAWALS_SCHEMA.into();
    }
    let proposed_bytes = catalog_bytes(&proposed)?;
    let mut after_files = Vec::new();
    for source in &context.candidates {
        let path = checked_path(&source.path)?;
        if path == context.target {
            after_files.push((source.scope.to_string(), path, proposed_bytes.clone()));
        } else {
            for row in &captured.files {
                if (row.1 == path || row.1.parent() == Some(path.as_path()))
                    && row.0 == source.scope
                {
                    after_files.push(row.clone());
                }
            }
        }
    }
    let after = effective(&after_files)?;
    for decision in &mut decisions {
        let actual = after
            .get(&decision.provider)
            .map(|(r, _)| record_sha256(r))
            .transpose()?;
        if actual != decision.proposed_record_sha256 {
            decision.diagnostics.push("masked_by_later_overlay".into());
        }
        if decision.diagnostics.iter().any(|d| {
            matches!(
                d.as_str(),
                "base_record_conflict"
                    | "masked_by_later_overlay"
                    | "withdraw_requires_existing_record"
            )
        }) {
            decision.action = "conflict".into();
        }
    }
    let history_key = sha256_prefixed_bytes(
        canonical_json_string(&(
            context.scope,
            &context.target,
            &captured.bindings,
            sha256_prefixed_bytes(bytes),
            sha256_prefixed_bytes(&proposed_bytes),
        ))
        .map_err(|e| error("json", e))?
        .as_bytes(),
    );
    let parent = context
        .target
        .parent()
        .ok_or_else(|| error("path", "missing parent"))?;
    let history = parent
        .join(".provider-update-history")
        .join(history_key.strip_prefix("sha256:").unwrap());
    let mut output_files = vec![
        context.target.display().to_string(),
        parent.join(".provider-update.lock").display().to_string(),
    ];
    for name in ["preview.json", "receipt.json", "proposed.json"] {
        output_files.push(history.join(name).display().to_string());
    }
    if prior.is_some() {
        output_files.push(history.join("previous.bin").display().to_string());
    }
    let mut report = ProviderUpdatePreview {
        schema: "gentle.provider_update_preview.v1".into(), scope: context.scope,
        target_path: context.target.display().to_string(), history_directory: history.display().to_string(), output_files, bundle,
        bundle_source: String::from_utf8(bytes.to_vec()).map_err(|e| error("utf8", e))?,
        bundle_bytes_sha256: sha256_prefixed_bytes(bytes), sources: captured.bindings,
        previous_target_sha256: prior.as_ref().map(|b| sha256_prefixed_bytes(b)),
        proposed_target_sha256: sha256_prefixed_bytes(&proposed_bytes), proposed_catalog: proposed,
        applicable: !decisions.iter().any(|d| d.action == "conflict") && decisions.iter().any(|d| d.action != "no_op"), decisions,
        impacts: vec!["Changed providers affect future discovery, capability checks and handoffs; regenerate those deliberately.".into(), "Products, commercial offers and reviewed procedures are not joined; affected planning results are not enumerated and saved reports are not rewritten.".into(), "Private stock, negotiated prices, equipment and observations are not read or modified.".into()],
        warnings: vec!["Hashes bind content, not supplier authenticity or signature trust. No document retrieval, freshness verdict, procedure review, quote or submission authorization.".into()],
        preview_sha256: String::new(),
    };
    report.preview_sha256 = sha256_prefixed_bytes(
        canonical_json_string(&report)
            .map_err(|e| error("json", e))?
            .as_bytes(),
    );
    Ok(report)
}

/// Save a fresh preview artifact; never overwrite any existing destination.
pub fn write_preview(path: &Path, preview: &ProviderUpdatePreview) -> Result<(), String> {
    let path = checked_path(path)?;
    for source in external_service_provider_config_discovery_candidates() {
        let source = checked_path(&source.path)?;
        if source == path || path.starts_with(&source) {
            return Err(error("preview_is_catalog_source", path.display()));
        }
    }
    let mut bytes = serde_json::to_vec_pretty(preview).map_err(|e| error("json", e))?;
    bytes.push(b'\n');
    if bytes.len() as u64 > MAX_BYTES {
        return Err(error("byte_budget", "preview artifact"));
    }
    let parent = path
        .parent()
        .ok_or_else(|| error("path", "missing parent"))?;
    let mut file = tempfile::NamedTempFile::new_in(parent).map_err(|e| error("write", e))?;
    file.write_all(&bytes)
        .and_then(|_| file.as_file().sync_all())
        .map_err(|e| error("write", e))?;
    file.persist_noclobber(path)
        .map_err(|e| error("preview_exists", e))?;
    Ok(())
}

/// Apply only an exactly approved, freshly reproducible preview. No project-state mutation.
pub fn apply(
    preview: &ProviderUpdatePreview,
    approval: &str,
) -> Result<ProviderUpdateReceipt, String> {
    apply_in(preview, approval, &Context::discover(preview.scope)?)
}

fn apply_in(
    preview: &ProviderUpdatePreview,
    approval: &str,
    context: &Context,
) -> Result<ProviderUpdateReceipt, String> {
    if approval != preview.preview_sha256 {
        return Err(error("approval_mismatch", "exact preview digest required"));
    }
    // Validate before creating any directories or lock file, then again under lock.
    let validate = || -> Result<(), String> {
        let current = preview_in(preview.bundle_source.as_bytes(), context)?;
        if !current.applicable
            || canonical_json_string(&current).map_err(|e| error("json", e))?
                != canonical_json_string(preview).map_err(|e| error("json", e))?
        {
            return Err(error("stale_or_tampered_preview", "regenerate preview"));
        }
        Ok(())
    };
    validate()?;
    let parent = context
        .target
        .parent()
        .ok_or_else(|| error("path", "missing parent"))?;
    fs::create_dir_all(parent).map_err(|e| error("directory", e))?;
    checked_path(parent)?;
    let lock_path = checked_path(&parent.join(".provider-update.lock"))?;
    let lock = fs::OpenOptions::new()
        .read(true)
        .write(true)
        .create(true)
        .truncate(false)
        .open(lock_path)
        .map_err(|e| error("lock", e))?;
    lock.try_lock().map_err(|e| error("busy", e))?;
    validate()?;
    let history = checked_path(&parent.join(".provider-update-history"))?;
    fs::create_dir_all(&history).map_err(|e| error("history", e))?;
    let destination = checked_path(Path::new(&preview.history_directory))?;
    if destination.exists() {
        return Err(error("recovery_required", destination.display()));
    }
    let receipt = ProviderUpdateReceipt {
        schema: "gentle.provider_update_receipt.v1".into(), preview_sha256: approval.into(),
        target_path: context.target.display().to_string(), previous_target_sha256: preview.previous_target_sha256.clone(),
        applied_target_sha256: preview.proposed_target_sha256.clone(), history_directory: destination.display().to_string(),
        commit_condition: "Compare target raw SHA-256: applied hash = committed; previous hash/absence = not committed; any other value = later change or manual recovery. Never auto-replay.".into(),
    };
    let staging = tempfile::tempdir_in(&history).map_err(|e| error("history", e))?;
    for (name, bytes) in [
        (
            "preview.json",
            serde_json::to_vec_pretty(preview).map_err(|e| error("json", e))?,
        ),
        (
            "receipt.json",
            serde_json::to_vec_pretty(&receipt).map_err(|e| error("json", e))?,
        ),
        ("proposed.json", catalog_bytes(&preview.proposed_catalog)?),
    ] {
        let mut file = fs::OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(staging.path().join(name))
            .map_err(|e| error("history", e))?;
        file.write_all(&bytes)
            .and_then(|_| file.sync_all())
            .map_err(|e| error("history", e))?;
    }
    if preview.previous_target_sha256.is_some() {
        let bytes = read_bytes(&context.target)?;
        if Some(sha256_prefixed_bytes(&bytes)) != preview.previous_target_sha256 {
            return Err(error("concurrent_writer", "target changed"));
        }
        let mut file = fs::OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(staging.path().join("previous.bin"))
            .map_err(|e| error("history", e))?;
        file.write_all(&bytes)
            .and_then(|_| file.sync_all())
            .map_err(|e| error("history", e))?;
    }
    fs::rename(staging.path(), &destination).map_err(|e| error("history", e))?;
    validate()?;
    let mut file = tempfile::NamedTempFile::new_in(parent).map_err(|e| error("write", e))?;
    file.write_all(&catalog_bytes(&preview.proposed_catalog)?)
        .and_then(|_| file.as_file().sync_all())
        .map_err(|e| error("write", e))?;
    file.persist(&context.target)
        .map_err(|e| error("commit_failed_recovery_retained", e))?;
    Ok(receipt)
}

#[cfg(test)]
mod tests;
