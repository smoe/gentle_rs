//! Thin shared shell parsing for strict TSS profile computation and report exports.

use super::*;
use gentle_protocol::tss_profiles::*;

pub(super) fn parse_tss_profiles_command(tokens: &[String]) -> Result<ShellCommand, String> {
    let export_only = tokens.get(1).map(String::as_str) == Some("tss-tfbs-profiles-export");
    let mut values = BTreeMap::new();
    let mut fasta = vec![];
    let mut index = 2;
    while index < tokens.len() {
        let flag = tokens[index].as_str();
        let allowed = match flag {
            "--output-dir" | "--formats" | "--scale-mode" | "--panels-per-page"
            | "--context-manifest" => true,
            "--report" => export_only,
            "--manifest"
            | "--panel"
            | "--fasta"
            | "--selection"
            | "--expected-genome-id"
            | "--expected-assembly"
            | "--expected-annotation-release"
            | "--expected-dataset-id" => !export_only,
            _ => false,
        };
        if !allowed {
            return Err(format!("Unknown option '{flag}' for {}", tokens[1]));
        }
        let value = tokens
            .get(index + 1)
            .filter(|s| !s.is_empty() && !s.starts_with("--"))
            .ok_or_else(|| format!("{flag} requires a value"))?
            .clone();
        if flag == "--fasta" {
            fasta.push(value);
        } else if values.insert(flag.to_string(), value).is_some() {
            return Err(format!("Duplicate option {flag}"));
        }
        index += 2;
    }
    let required = |flag: &str| {
        values
            .get(flag)
            .cloned()
            .ok_or_else(|| format!("{} requires {flag}", tokens[1]))
    };
    let scale_mode = values
        .get("--scale-mode")
        .map(|raw| match raw.as_str() {
            "shared" => Ok(TssScaleMode::Shared),
            "independent" => Ok(TssScaleMode::Independent),
            _ => Err("--scale-mode must be shared or independent".to_string()),
        })
        .transpose()?;
    let panels_per_page = values
        .get("--panels-per-page")
        .map(|s| s.parse::<usize>())
        .transpose()
        .map_err(|_| "--panels-per-page requires a positive integer")?
        .unwrap_or(1);
    if panels_per_page == 0 || panels_per_page > 32 {
        return Err("--panels-per-page must be 1..32".into());
    }
    let formats = values
        .get("--formats")
        .map(|s| s.as_str())
        .unwrap_or("svg")
        .split(',')
        .map(|s| match s {
            "svg" => Ok(TssExportFormat::Svg),
            "png" => Ok(TssExportFormat::Png),
            "pdf" => Ok(TssExportFormat::Pdf),
            _ => Err("--formats accepts svg,png,pdf".to_string()),
        })
        .collect::<Result<Vec<_>, _>>()?;
    for (i, format) in formats.iter().enumerate() {
        if formats[..i].contains(format) {
            return Err("Duplicate export format".into());
        }
    }
    let export = ExportTssProfilesRequest {
        context_manifest: values.get("--context-manifest").cloned(),
        output_dir: required("--output-dir")?,
        rendering: TssProfileRenderOptions {
            scale_mode,
            panels_per_page,
        },
        formats,
    };
    let op = if export_only {
        let path = required("--report")?;
        let file = crate::tss_fasta_bundle::open_regular_input(Path::new(&path))
            .map_err(|e| format!("Report file: {e}"))?;
        let mut bytes = vec![];
        std::io::Read::read_to_end(
            &mut std::io::Read::take(file, 256 * 1024 * 1024 + 1),
            &mut bytes,
        )
        .map_err(|e| format!("Report read: {e}"))?;
        if bytes.len() > 256 * 1024 * 1024 {
            return Err("Report exceeds 256 MiB".into());
        }
        let report = serde_json::from_slice(&bytes)
            .map_err(|e| format!("Invalid TSS profile report: {e}"))?;
        Operation::ExportTssTfbsProfiles {
            report: Box::new(report),
            request: export,
        }
    } else {
        let request = ComputeTssProfilesRequest {
            manifest: required("--manifest")?,
            panel: required("--panel")?,
            fasta,
            selection: values.get("--selection").cloned(),
            expected_genome_id: required("--expected-genome-id")?,
            expected_assembly: values.get("--expected-assembly").cloned(),
            expected_annotation_release: values.get("--expected-annotation-release").cloned(),
            expected_dataset_id: values.get("--expected-dataset-id").cloned(),
        };
        Operation::ComputeTssTfbsProfiles {
            request: Box::new(request),
            export: Some(export),
        }
    };
    Ok(ShellCommand::Op {
        payload: serde_json::to_string(&op).map_err(|e| e.to_string())?,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    fn parse(args: &str) -> Result<ShellCommand, String> {
        let tokens = args
            .split_whitespace()
            .map(str::to_string)
            .collect::<Vec<_>>();
        parse_tss_profiles_command(&tokens)
    }
    #[test]
    fn tss_profiles_parser_keeps_repeated_fasta_and_exact_reference() {
        let parsed=parse("features tss-tfbs-profiles --manifest m.json --panel p.json --expected-genome-id ref --fasta a.fa --fasta b.fa --output-dir out --formats svg,pdf").unwrap();
        let ShellCommand::Op { payload } = parsed else {
            panic!("typed operation")
        };
        let Operation::ComputeTssTfbsProfiles { request, export } =
            serde_json::from_str(&payload).unwrap()
        else {
            panic!("compute")
        };
        assert_eq!(request.fasta, vec!["a.fa", "b.fa"]);
        assert_eq!(request.expected_genome_id, "ref");
        assert_eq!(export.unwrap().formats.len(), 2);
    }

    #[test]
    fn tss_profiles_parser_preserves_whole_genome_and_dataset_labels() {
        let tokens = [
            "features",
            "tss-tfbs-profiles",
            "--manifest",
            "manifest.json",
            "--panel",
            "panel.json",
            "--output-dir",
            "output",
            "--expected-genome-id",
            "Synthetic assembly release 7",
            "--expected-assembly",
            "synthetic-assembly",
            "--expected-dataset-id",
            "synthetic_release7_promoterome",
        ]
        .map(str::to_string);
        let ShellCommand::Op { payload } = parse_tss_profiles_command(&tokens).unwrap() else {
            panic!("typed operation");
        };
        let Operation::ComputeTssTfbsProfiles { request, .. } =
            serde_json::from_str(&payload).unwrap()
        else {
            panic!("compute operation");
        };
        assert_eq!(request.expected_genome_id, "Synthetic assembly release 7");
        assert_eq!(
            request.expected_dataset_id.as_deref(),
            Some("synthetic_release7_promoterome")
        );
        assert_eq!(
            request.expected_assembly.as_deref(),
            Some("synthetic-assembly")
        );
        assert!(request.expected_annotation_release.is_none());
    }

    #[test]
    fn tss_profiles_parser_forwards_context_for_compute_and_report_only_export() {
        for id in [
            "features tss-tfbs-profiles",
            "ComputeTssTfbsProfiles",
            "features tss-tfbs-profiles-export",
            "ExportTssTfbsProfiles",
        ] {
            let descriptor = tss_profile_capability_descriptor(
                id,
                !id.contains("export") && id != "ExportTssTfbsProfiles",
            );
            assert!(
                descriptor["args"]
                    .as_array()
                    .unwrap()
                    .iter()
                    .any(|arg| arg["name"] == "CONTEXT_MANIFEST" && arg["required"] == false)
            );
        }
        let parsed = parse("features tss-tfbs-profiles --manifest m --panel p --expected-genome-id g --output-dir out --context-manifest context.json").unwrap();
        let ShellCommand::Op { payload } = parsed else {
            panic!("typed operation");
        };
        let Operation::ComputeTssTfbsProfiles { export, .. } =
            serde_json::from_str(&payload).unwrap()
        else {
            panic!("compute");
        };
        assert_eq!(
            export.unwrap().context_manifest.as_deref(),
            Some("context.json")
        );
        let temp = tempfile::tempdir().unwrap();
        let path = temp.path().join("report.json");
        std::fs::write(
            &path,
            serde_json::to_vec(&crate::tss_profile_export::tests::synthetic_report()).unwrap(),
        )
        .unwrap();
        let tokens = [
            "features",
            "tss-tfbs-profiles-export",
            "--report",
            path.to_str().unwrap(),
            "--output-dir",
            "out",
            "--context-manifest",
            "context with spaces.json",
        ]
        .map(str::to_string);
        let ShellCommand::Op { payload } = parse_tss_profiles_command(&tokens).unwrap() else {
            panic!("typed export");
        };
        let Operation::ExportTssTfbsProfiles { request, .. } =
            serde_json::from_str(&payload).unwrap()
        else {
            panic!("export");
        };
        assert_eq!(
            request.context_manifest.as_deref(),
            Some("context with spaces.json")
        );
        assert!(
            parse("features tss-tfbs-profiles --context-manifest a --context-manifest b")
                .unwrap_err()
                .contains("Duplicate")
        );
    }
    #[test]
    fn tss_profiles_parser_rejects_unknown_missing_and_duplicate_options() {
        for args in [
            "features tss-tfbs-profiles --wat value",
            "features tss-tfbs-profiles --panel",
            "features tss-tfbs-profiles --manifest m --output-dir out --expected-genome-id ref",
            "features tss-tfbs-profiles --panel a --panel b",
            "features tss-tfbs-profiles --panels-per-page 0",
        ] {
            assert!(parse(args).is_err(), "{args}");
        }
    }

    #[test]
    fn tss_profiles_shared_command_forwards_cancellation_without_writes() {
        let _registry_guard = crate::tf_motifs::test_registry_lock().lock().unwrap();
        let fixture = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("test_files/fixtures/tss_profiles");
        let scratch = tempfile::tempdir().unwrap();
        let output = std::fs::canonicalize(scratch.path())
            .unwrap()
            .join("cancelled");
        let tokens = vec![
            "features".into(),
            "tss-tfbs-profiles".into(),
            "--manifest".into(),
            fixture.join("manifest.json").to_string_lossy().into_owned(),
            "--panel".into(),
            fixture.join("panel.json").to_string_lossy().into_owned(),
            "--expected-genome-id".into(),
            "synthetic-genome-v1".into(),
            "--output-dir".into(),
            output.to_string_lossy().into_owned(),
        ];
        let command = parse_shell_tokens(&tokens).unwrap();
        let called = Arc::new(std::sync::atomic::AtomicBool::new(false));
        let capture = called.clone();
        let callback: ShellProgressCallback = Arc::new(Mutex::new(Box::new(move |_| {
            capture.store(true, std::sync::atomic::Ordering::SeqCst);
            false
        })));
        let mut engine = GentleEngine::from_state(ProjectState::default());
        let before = serde_json::to_value(engine.state()).unwrap();
        let result = execute_shell_command_with_options(
            &mut engine,
            &command,
            &ShellExecutionOptions {
                progress_callback: Some(callback),
                ..Default::default()
            },
        );
        assert!(result.unwrap_err().to_ascii_lowercase().contains("cancel"));
        assert!(called.load(std::sync::atomic::Ordering::SeqCst));
        assert!(!output.exists());
        assert_eq!(serde_json::to_value(engine.state()).unwrap(), before);
    }
}
