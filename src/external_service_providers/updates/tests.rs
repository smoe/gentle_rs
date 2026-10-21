//! Hand-crafted provider/overlay fixtures, recreated inline; no real supplier facts or contacts.

use super::*;
use serde_json::json;

fn record(id: &str, name: &str) -> ExternalServiceProviderConfigRecord {
    serde_json::from_value(json!({
        "provider":id,"display_name":name,"support_status":"catalog_only",
        "website_url":"https://example.invalid","dashboard_url":"https://example.invalid",
        "capabilities":[{"service_kind":"synthetic_protocol","track":"dna","display_name":"Synthetic","quote_handoff_supported":false,"direct_api_documented":false,"direct_api_implemented":false,"supported_submission_modes":[],"status_tracking":"none","artifact_kinds":[]}]
    })).unwrap()
}

fn catalog(rows: Vec<ExternalServiceProviderConfigRecord>) -> ExternalServiceProviderConfigCatalog {
    ExternalServiceProviderConfigCatalog {
        schema: EXTERNAL_SERVICE_PROVIDER_CONFIG_SCHEMA.into(),
        providers: rows,
        ..Default::default()
    }
}

fn fixture() -> (tempfile::TempDir, Context) {
    let temp = tempfile::tempdir().unwrap();
    let root = temp.path().canonicalize().unwrap();
    let built_in = root.join("built-in.json");
    fs::write(
        &built_in,
        catalog_bytes(&catalog(vec![record("synthetic", "Original")])).unwrap(),
    )
    .unwrap();
    let target = root
        .join("project")
        .join("catalogs")
        .join("external_service_providers.json");
    let candidates = vec![
        ExternalServiceProviderConfigSourceCandidate {
            scope: "built-in",
            path: built_in,
        },
        ExternalServiceProviderConfigSourceCandidate {
            scope: "project",
            path: target.clone(),
        },
        ExternalServiceProviderConfigSourceCandidate {
            scope: "project",
            path: target
                .parent()
                .unwrap()
                .join("external_service_providers.d"),
        },
    ];
    (
        temp,
        Context {
            scope: UpdateScope::Project,
            target,
            candidates,
        },
    )
}

fn replace() -> Vec<u8> {
    serde_json::to_vec_pretty(&ProviderUpdateBundle {
        schema: "gentle.provider_update_bundle.v1".into(),
        changes: vec![ProviderChange::Replace {
            provider: "synthetic".into(),
            base_record_sha256: Some(record_sha256(&record("synthetic", "Original")).unwrap()),
            record: serde_json::to_value(record("synthetic", "Replacement")).unwrap(),
        }],
    })
    .unwrap()
}

fn discovered_fixture_index(context: &Context) -> ExternalServiceProviderConfigIndex {
    let sources = context
        .candidates
        .iter()
        .filter(|source| source.path.try_exists().expect("fixture source metadata"))
        .cloned()
        .collect::<Vec<_>>();
    ExternalServiceProviderConfigIndex::from_sources(&sources, "synthetic".into()).unwrap()
}

#[cfg(windows)]
#[test]
fn provider_update_windows_drive_and_verbatim_paths_preserve_raw_traversal_checks() {
    let (temp, _context) = fixture();
    let plain = temp.path().join("built-in.json");
    let canonical = plain.canonicalize().unwrap();
    let expected = fs::read(&plain).unwrap();
    for path in [&plain, &canonical] {
        assert_eq!(
            checked_path(path).unwrap(),
            std::path::absolute(path).unwrap()
        );
        let (_, bytes) = read_json::<Value>(path).unwrap();
        assert_eq!(bytes, expected);
    }

    let mut raw = temp.path().as_os_str().to_os_string();
    raw.push("\\unused\\..\\built-in.json");
    let traversal = PathBuf::from(raw);
    assert!(
        traversal
            .components()
            .any(|part| matches!(part, std::path::Component::ParentDir))
    );
    assert!(
        checked_path(&traversal)
            .unwrap_err()
            .starts_with("catalog_update:unsafe_path:")
    );

    let mut raw = canonical.parent().unwrap().as_os_str().to_os_string();
    raw.push("\\unused\\..\\built-in.json");
    let traversal = PathBuf::from(raw);
    assert!(traversal.components().any(|part| part.as_os_str() == ".."));
    assert!(
        checked_path(&traversal)
            .unwrap_err()
            .starts_with("catalog_update:unsafe_path:")
    );
}

#[test]
fn provider_update_preview_apply_retains_raw_history_and_emits_receipt() {
    let (_temp, context) = fixture();
    let preview = preview_in(&replace(), &context).unwrap();
    assert!(preview.applicable);
    assert_eq!(preview.decisions[0].action, "replace");
    assert!(!context.target.exists());
    let receipt = apply_in(&preview, &preview.preview_sha256, &context).unwrap();
    let bytes = fs::read(&context.target).unwrap();
    assert_eq!(sha256_prefixed_bytes(&bytes), receipt.applied_target_sha256);
    let history = Path::new(&receipt.history_directory);
    assert!(!history.join("previous.bin").exists()); // previously absent, not empty
    assert_eq!(fs::read(history.join("proposed.json")).unwrap(), bytes);
    let index = discovered_fixture_index(&context);
    assert_eq!(
        index.provider("synthetic").unwrap().record.display_name,
        "Replacement"
    );
    assert!(apply_in(&preview, &preview.preview_sha256, &context).is_err());
}

#[test]
fn provider_update_withdrawal_masks_lower_record_and_legacy_overlay_still_works() {
    let (_temp, context) = fixture();
    let bytes = serde_json::to_vec(&ProviderUpdateBundle {
        schema: "gentle.provider_update_bundle.v1".into(),
        changes: vec![ProviderChange::Withdraw {
            provider: "synthetic".into(),
            base_record_sha256: record_sha256(&record("synthetic", "Original")).unwrap(),
        }],
    })
    .unwrap();
    let preview = preview_in(&bytes, &context).unwrap();
    assert_eq!(preview.decisions[0].action, "withdraw");
    assert_eq!(
        preview.proposed_catalog.schema,
        EXTERNAL_SERVICE_PROVIDER_CONFIG_WITHDRAWALS_SCHEMA
    );
    apply_in(&preview, &preview.preview_sha256, &context).unwrap();
    let index = discovered_fixture_index(&context);
    assert!(index.provider("synthetic").is_none());
    assert!(index.provider_ids().is_empty());
    assert!(index.to_provider_catalog(0).providers.is_empty());
    // Legacy full-record overlay can deliberately restore the provider later.
    let directory = &context.candidates[2].path;
    fs::create_dir(directory).unwrap();
    fs::write(
        directory.join("restore.json"),
        catalog_bytes(&catalog(vec![record("synthetic", "Restored")])).unwrap(),
    )
    .unwrap();
    let index = discovered_fixture_index(&context);
    assert_eq!(
        index.provider("synthetic").unwrap().record.display_name,
        "Restored"
    );
}

#[test]
fn provider_update_stale_overlay_raw_bytes_membership_and_target_are_rejected() {
    for changed in ["raw", "member", "target"] {
        let (_temp, context) = fixture();
        let preview = preview_in(&replace(), &context).unwrap();
        match changed {
            "raw" => {
                let source = &context.candidates[0].path;
                let bytes = fs::read_to_string(source).unwrap().replace('\n', "\r\n");
                fs::write(source, bytes).unwrap();
            }
            "member" => {
                let directory = &context.candidates[2].path;
                fs::create_dir_all(directory).unwrap();
                fs::write(
                    directory.join("new.json"),
                    catalog_bytes(&catalog(vec![])).unwrap(),
                )
                .unwrap();
            }
            _ => {
                fs::create_dir_all(context.target.parent().unwrap()).unwrap();
                fs::write(&context.target, catalog_bytes(&catalog(vec![])).unwrap()).unwrap();
            }
        }
        let before = fs::read(&context.target).ok();
        assert!(
            apply_in(&preview, &preview.preview_sha256, &context)
                .unwrap_err()
                .contains("stale_or_tampered_preview")
        );
        assert_eq!(fs::read(&context.target).ok(), before);
        assert!(
            !context
                .target
                .parent()
                .unwrap()
                .join(".provider-update-history")
                .exists()
        );
    }
}

#[test]
fn provider_update_detects_base_conflict_masking_and_rejects_incomplete_records() {
    let (_temp, context) = fixture();
    let mut bundle: Value = serde_json::from_slice(&replace()).unwrap();
    bundle["changes"][0]["base_record_sha256"] = Value::Null;
    let preview = preview_in(&serde_json::to_vec(&bundle).unwrap(), &context).unwrap();
    assert!(!preview.applicable);
    assert_eq!(preview.decisions[0].action, "conflict");
    let directory = &context.candidates[2].path;
    fs::create_dir_all(directory).unwrap();
    fs::write(
        directory.join("later.json"),
        catalog_bytes(&catalog(vec![record("synthetic", "Original")])).unwrap(),
    )
    .unwrap();
    let preview = preview_in(&replace(), &context).unwrap();
    assert!(!preview.applicable);
    assert!(
        preview.decisions[0]
            .diagnostics
            .contains(&"masked_by_later_overlay".into())
    );
    let mut bundle: Value = serde_json::from_slice(&replace()).unwrap();
    bundle["changes"][0]["record"]
        .as_object_mut()
        .unwrap()
        .remove("channels");
    assert!(
        preview_in(&serde_json::to_vec(&bundle).unwrap(), &context)
            .unwrap_err()
            .contains("complete_record_required")
    );
    bundle["changes"][0]["record"]["channels"] = json!([]);
    bundle["changes"][0]["record"]["invented"] = json!(true);
    assert!(preview_in(&serde_json::to_vec(&bundle).unwrap(), &context).is_err());
}

#[test]
fn provider_update_rejects_rehashed_tampering_wrong_scope_and_approval_without_writes() {
    for changed in ["target", "record", "source", "hash", "bundle"] {
        let (_temp, context) = fixture();
        let mut preview = preview_in(&replace(), &context).unwrap();
        match changed {
            "target" => preview.target_path.push_str(".other"),
            "record" => preview.proposed_catalog.providers[0].display_name = "Forged".into(),
            "source" => preview.sources.clear(),
            "hash" => preview.bundle_bytes_sha256 = "sha256:forged".into(),
            _ => preview.bundle.changes.clear(),
        }
        preview.preview_sha256.clear();
        preview.preview_sha256 =
            sha256_prefixed_bytes(canonical_json_string(&preview).unwrap().as_bytes());
        assert!(
            apply_in(&preview, &preview.preview_sha256, &context).is_err(),
            "{changed}"
        );
        assert!(!context.target.parent().unwrap().exists());
    }
    let (_temp, context) = fixture();
    let preview = preview_in(&replace(), &context).unwrap();
    assert!(apply_in(&preview, "sha256:wrong", &context).is_err());
    let mut other = context.clone();
    other.scope = UpdateScope::User;
    assert!(apply_in(&preview, &preview.preview_sha256, &other).is_err());
    assert!(!context.target.parent().unwrap().exists());
}

#[test]
fn provider_update_add_and_no_op_do_not_fabricate_procedure_evidence() {
    let (_temp, context) = fixture();
    let bundle = ProviderUpdateBundle {
        schema: "gentle.provider_update_bundle.v1".into(),
        changes: vec![ProviderChange::Replace {
            provider: "new".into(),
            base_record_sha256: None,
            record: serde_json::to_value(record("new", "New")).unwrap(),
        }],
    };
    let preview = preview_in(&serde_json::to_vec(&bundle).unwrap(), &context).unwrap();
    assert_eq!(preview.decisions[0].action, "add");
    assert!(
        preview.decisions[0]
            .diagnostics
            .contains(&"protocol_document_evidence_unavailable".into())
    );
    let mut bundle: ProviderUpdateBundle = serde_json::from_slice(&replace()).unwrap();
    if let ProviderChange::Replace { record: row, .. } = &mut bundle.changes[0] {
        *row = serde_json::to_value(record("synthetic", "Original")).unwrap();
    }
    let preview = preview_in(&serde_json::to_vec(&bundle).unwrap(), &context).unwrap();
    assert_eq!(preview.decisions[0].action, "no_op");
    assert!(!preview.applicable);
}

#[test]
fn provider_update_lf_crlf_preserves_prior_bytes_and_same_byte_parse_hash() {
    for newline in ["\n", "\r\n"] {
        let (_temp, context) = fixture();
        fs::create_dir_all(context.target.parent().unwrap()).unwrap();
        let prior = String::from_utf8(catalog_bytes(&catalog(vec![])).unwrap())
            .unwrap()
            .replace('\n', newline)
            .into_bytes();
        fs::write(&context.target, &prior).unwrap();
        let input = String::from_utf8(replace())
            .unwrap()
            .replace('\n', newline)
            .into_bytes();
        let preview = preview_in(&input, &context).unwrap();
        assert_eq!(preview.bundle_bytes_sha256, sha256_prefixed_bytes(&input));
        let receipt = apply_in(&preview, &preview.preview_sha256, &context).unwrap();
        assert_eq!(
            fs::read(Path::new(&receipt.history_directory).join("previous.bin")).unwrap(),
            prior
        );
        assert_eq!(
            receipt.previous_target_sha256,
            Some(sha256_prefixed_bytes(&prior))
        );
    }
}

#[test]
fn provider_update_raw_paths_and_shell_quoted_cross_platform_strings_round_trip() {
    let raw = "relative/../escape.json";
    assert!(
        Path::new(&raw)
            .components()
            .any(|p| p == std::path::Component::ParentDir)
    );
    assert!(
        read_json::<Value>(Path::new(&raw))
            .unwrap_err()
            .contains("unsafe_path")
    );
    for path in [
        "a directory/a 'quoted' bundle.json",
        r"C:\Users\test\bundle.json",
        r"\\server\share\bundle.json",
        r"\\?\C:\test\bundle.json",
    ] {
        let line = format!(
            "catalog updates preview {} --scope project --output preview.json",
            crate::engine_shell::quote_shell_arg(path)
        );
        let parsed = crate::engine_shell::parse_shell_line(&line).unwrap();
        match parsed {
            crate::engine_shell::ShellCommand::Catalog(
                crate::catalog::CatalogCommand::Preview { bundle_path, .. },
            ) => assert_eq!(bundle_path, path),
            _ => panic!("wrong command"),
        }
    }
}

#[cfg(windows)]
#[test]
fn provider_update_rejects_verbatim_parent_components_before_normalization() {
    let raw = Path::new(r"\\?\C:\synthetic\..\escape.json");
    assert!(raw.components().any(|part| part.as_os_str() == ".."));
    assert!(read_json::<Value>(raw).unwrap_err().contains("unsafe_path"));
}

#[test]
fn provider_update_preview_artifacts_are_fresh_and_hash_bound_on_reopen() {
    let (temp, context) = fixture();
    let preview = preview_in(&replace(), &context).unwrap();
    let path = temp
        .path()
        .canonicalize()
        .unwrap()
        .join("preview with spaces.json");
    write_preview(&path, &preview).unwrap();
    let saved = fs::read(&path).unwrap();
    assert!(write_preview(&path, &preview).is_err());
    assert_eq!(fs::read(&path).unwrap(), saved);
    let (reopened, bytes): (ProviderUpdatePreview, _) = read_json(&path).unwrap();
    assert_eq!(bytes, saved);
    let receipt = apply_in(&reopened, &reopened.preview_sha256, &context).unwrap();
    assert_eq!(receipt.history_directory, preview.history_directory);
    assert!(preview.output_files.contains(&receipt.target_path));
    assert!(
        preview.output_files.contains(
            &Path::new(&receipt.history_directory)
                .join("receipt.json")
                .display()
                .to_string()
        )
    );
}

#[test]
fn provider_update_recovery_directory_blocks_retry_without_overwrite() {
    let (_temp, context) = fixture();
    let preview = preview_in(&replace(), &context).unwrap();
    let history = PathBuf::from(&preview.history_directory);
    fs::create_dir_all(&history).unwrap();
    fs::write(history.join("retained.bin"), b"synthetic recovery evidence").unwrap();
    assert!(
        apply_in(&preview, &preview.preview_sha256, &context)
            .unwrap_err()
            .contains("recovery_required")
    );
    assert!(!context.target.exists());
    assert_eq!(
        fs::read(history.join("retained.bin")).unwrap(),
        b"synthetic recovery evidence"
    );
}

#[test]
fn provider_update_writer_lock_is_immediate_and_history_is_not_written_when_busy() {
    let (_temp, context) = fixture();
    let preview = preview_in(&replace(), &context).unwrap();
    let parent = context.target.parent().unwrap();
    fs::create_dir_all(parent).unwrap();
    let lock = fs::OpenOptions::new()
        .read(true)
        .write(true)
        .create(true)
        .truncate(false)
        .open(parent.join(".provider-update.lock"))
        .unwrap();
    lock.try_lock().unwrap();
    assert!(
        apply_in(&preview, &preview.preview_sha256, &context)
            .unwrap_err()
            .contains("busy")
    );
    assert!(!context.target.exists());
    assert!(!parent.join(".provider-update-history").exists());
}

#[test]
fn provider_update_admission_budgets_and_duplicate_tombstones_are_explicit() {
    let (_temp, context) = fixture();
    assert!(
        preview_in(&vec![b' '; 512 * 1024 + 1], &context)
            .unwrap_err()
            .contains("byte_budget")
    );
    let mut bundle: Value = serde_json::from_slice(&replace()).unwrap();
    let change = bundle["changes"][0].clone();
    bundle["changes"] = json!(vec![change; 65]);
    assert!(
        preview_in(&serde_json::to_vec(&bundle).unwrap(), &context)
            .unwrap_err()
            .contains("bundle_admission")
    );
    let mut invalid = catalog(vec![]);
    invalid.withdrawn_providers = vec!["synthetic".into(), "Synthetic".into()];
    assert!(!validate_provider_config_catalog(&invalid, &context.target).is_empty());
    invalid.withdrawn_providers = vec!["synthetic".into()];
    invalid.providers.push(record("synthetic", "Overlap"));
    assert!(!validate_provider_config_catalog(&invalid, &context.target).is_empty());
}

#[test]
fn provider_update_tombstones_require_v2_and_legacy_bytes_remain_unchanged() {
    let (_temp, context) = fixture();
    let legacy = catalog(vec![]);
    let output = serde_json::to_value(&legacy).unwrap();
    assert_eq!(
        serde_json::to_string(&legacy).unwrap(),
        r#"{"schema":"gentle.external_service_provider_config.v1","providers":[],"summary_lines":[]}"#
    );
    assert!(output.get("withdrawn_providers").is_none());
    assert_eq!(output["schema"], EXTERNAL_SERVICE_PROVIDER_CONFIG_SCHEMA);
    let mut withdrawal = legacy;
    withdrawal.withdrawn_providers.push("synthetic".into());
    assert!(!validate_provider_config_catalog(&withdrawal, &context.target).is_empty());
    withdrawal.schema = EXTERNAL_SERVICE_PROVIDER_CONFIG_WITHDRAWALS_SCHEMA.into();
    assert!(validate_provider_config_catalog(&withdrawal, &context.target).is_empty());
}

#[test]
fn provider_update_one_conflict_blocks_all_other_record_changes() {
    let (_temp, context) = fixture();
    let mut bundle: ProviderUpdateBundle = serde_json::from_slice(&replace()).unwrap();
    if let ProviderChange::Replace {
        base_record_sha256, ..
    } = &mut bundle.changes[0]
    {
        *base_record_sha256 = None;
    }
    bundle.changes.push(ProviderChange::Replace {
        provider: "another_synthetic".into(),
        base_record_sha256: None,
        record: serde_json::to_value(record("another_synthetic", "New")).unwrap(),
    });
    let preview = preview_in(&serde_json::to_vec(&bundle).unwrap(), &context).unwrap();
    assert_eq!(preview.decisions[0].action, "conflict");
    assert_eq!(preview.decisions[1].action, "add");
    assert!(apply_in(&preview, &preview.preview_sha256, &context).is_err());
    assert!(!context.target.parent().unwrap().exists());
}

#[test]
fn provider_update_unknown_bundle_duplicate_ids_and_semantic_errors_fail_closed() {
    let (_temp, context) = fixture();
    for change in ["unknown", "duplicate", "country", "empty"] {
        let mut bundle: Value = serde_json::from_slice(&replace()).unwrap();
        match change {
            "unknown" => bundle["untrusted_patch"] = json!({}),
            "duplicate" => {
                let row = bundle["changes"][0].clone();
                bundle["changes"].as_array_mut().unwrap().push(row);
            }
            "country" => {
                bundle["changes"][0]["record"]["country_code"] = json!("DE");
            }
            _ => bundle["changes"] = json!([]),
        }
        assert!(
            preview_in(&serde_json::to_vec(&bundle).unwrap(), &context).is_err(),
            "{change}"
        );
    }
    assert!(!context.target.parent().unwrap().exists());
}

#[cfg(unix)]
#[test]
fn provider_update_rejects_symlink_sources_and_destinations() {
    let (_temp, context) = fixture();
    let link = context.candidates[0].path.with_file_name("link.json");
    std::os::unix::fs::symlink(&context.candidates[0].path, &link).unwrap();
    assert!(read_json::<Value>(&link).unwrap_err().contains("symlink"));
    let mut linked = context.clone();
    linked.target = link;
    let preview = preview_in(&replace(), &context).unwrap();
    assert!(apply_in(&preview, &preview.preview_sha256, &linked).is_err());
}
