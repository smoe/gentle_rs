#![allow(dead_code)]
//! Phase-1 GENtle-side ClawBio adapter model.

use super::*;
use serde_json::Value;
use std::{
    fs,
    path::{Path, PathBuf},
    process::{Command, Stdio},
    sync::{
        Arc,
        atomic::{AtomicBool, Ordering},
    },
    thread,
    time::Duration,
};

pub(crate) const DEFAULT_CLAWBIO_SKILL_ALIAS: &str = "gentle-cloning";

#[derive(Clone, Default)]
pub(crate) struct CancelToken(Arc<AtomicBool>);

impl CancelToken {
    pub(crate) fn cancel(&self) {
        self.0.store(true, Ordering::Relaxed);
    }

    pub(crate) fn is_cancelled(&self) -> bool {
        self.0.load(Ordering::Relaxed)
    }
}

#[derive(Clone)]
pub(crate) struct ClawBioRunSpec {
    pub(crate) clawbio_bin: String,
    pub(crate) skill_alias: String,
    pub(crate) request: Value,
    pub(crate) run_root: PathBuf,
}

#[derive(Clone)]
pub(crate) struct ClawBioResultView {
    pub(crate) output_dir: String,
    pub(crate) result_json_path: String,
    pub(crate) report_md_path: Option<String>,
    pub(crate) workflow_state: Option<Value>,
    pub(crate) chat_summary_lines: Vec<String>,
    pub(crate) preferred_artifacts: Vec<Value>,
    pub(crate) suggested_actions: Vec<Value>,
    pub(crate) raw_result: Value,
}

pub(crate) trait ClawBioTransport {
    fn dispatch(
        &self,
        spec: ClawBioRunSpec,
        cancel: CancelToken,
    ) -> Result<ClawBioResultView, String>;
}

pub(crate) struct SubprocessTransport;

impl ClawBioTransport for SubprocessTransport {
    fn dispatch(
        &self,
        spec: ClawBioRunSpec,
        cancel: CancelToken,
    ) -> Result<ClawBioResultView, String> {
        let output_dir = next_clawbio_run_dir(&spec.run_root, &spec.skill_alias)?;
        fs::create_dir_all(&output_dir)
            .map_err(|e| format!("Could not create ClawBio output directory: {e}"))?;
        let request_path = output_dir.join("request.json");
        fs::write(
            &request_path,
            serde_json::to_string_pretty(&spec.request)
                .map_err(|e| format!("Could not serialize ClawBio request: {e}"))?
                + "\n",
        )
        .map_err(|e| format!("Could not write ClawBio request: {e}"))?;

        let mut child = Command::new(spec.clawbio_bin.trim())
            .arg("run")
            .arg(&spec.skill_alias)
            .arg("--input")
            .arg(&request_path)
            .arg("--output")
            .arg(&output_dir)
            .stdout(Stdio::null())
            .stderr(Stdio::null())
            .spawn()
            .map_err(|e| format!("Could not start ClawBio: {e}"))?;

        let status = loop {
            if cancel.is_cancelled() {
                let _ = child.kill();
                let _ = child.wait();
                return Err("ClawBio run cancelled.".to_string());
            }
            match child
                .try_wait()
                .map_err(|e| format!("Could not poll ClawBio process: {e}"))?
            {
                Some(status) => break status,
                None => thread::sleep(Duration::from_millis(100)),
            }
        };
        if !status.success() {
            return Err(format!("ClawBio exited with status {status}."));
        }
        ClawBioResultView::read_from_output_dir(&output_dir)
    }
}

impl ClawBioResultView {
    pub(crate) fn read_from_output_dir(output_dir: &Path) -> Result<Self, String> {
        let result_json_path = output_dir.join("result.json");
        let raw_text = fs::read_to_string(&result_json_path)
            .map_err(|e| format!("Could not read ClawBio result.json: {e}"))?;
        let raw_result: Value = serde_json::from_str(&raw_text)
            .map_err(|e| format!("Could not parse ClawBio result.json: {e}"))?;
        let stdout_json = raw_result.get("stdout_json");
        let workflow_state = raw_result
            .get("workflow_state")
            .or_else(|| stdout_json.and_then(|v| v.get("workflow_state")))
            .cloned();
        Ok(Self {
            output_dir: output_dir.display().to_string(),
            result_json_path: result_json_path.display().to_string(),
            report_md_path: raw_result
                .pointer("/artifacts/report_md")
                .and_then(Value::as_str)
                .or_else(|| raw_result.get("report_md").and_then(Value::as_str))
                .map(|path| resolve_result_path(output_dir, path)),
            workflow_state,
            chat_summary_lines: string_list_field(&raw_result, "chat_summary_lines")
                .or_else(|| stdout_json.and_then(|v| string_list_field(v, "chat_summary_lines")))
                .unwrap_or_default(),
            preferred_artifacts: value_array_field(&raw_result, "preferred_artifacts")
                .or_else(|| stdout_json.and_then(|v| value_array_field(v, "preferred_artifacts")))
                .unwrap_or_default(),
            suggested_actions: value_array_field(&raw_result, "suggested_actions")
                .or_else(|| stdout_json.and_then(|v| value_array_field(v, "suggested_actions")))
                .unwrap_or_default(),
            raw_result,
        })
    }
}

impl GENtleApp {
    pub(crate) fn clawbio_context_request(&self, skill_alias: Option<&str>) -> Value {
        serde_json::json!({
            "schema": "gentle.clawbio_skill_request.v1",
            "mode": "skill-info",
            "skill_alias": skill_alias.unwrap_or(DEFAULT_CLAWBIO_SKILL_ALIAS),
            "gentle_context": self.clawbio_context_payload(),
        })
    }

    pub(crate) fn clawbio_context_payload(&self) -> Value {
        let (seq_id, selected_region) = self
            .active_dna_window_context()
            .map(|(seq_id, selection)| {
                let region = selection.map(|(start, end)| {
                    serde_json::json!({
                        "start_0based": start,
                        "end_0based": end,
                        "strand": Value::Null
                    })
                });
                (Some(seq_id), region)
            })
            .unwrap_or((None, None));
        let resources = seq_id
            .as_ref()
            .and_then(|id| {
                self.engine
                    .read()
                    .ok()
                    .and_then(|engine| engine.sequence_genome_anchor_summary(id).ok())
            })
            .map(|anchor| self.clawbio_anchor_resource(anchor))
            .into_iter()
            .collect::<Vec<_>>();
        serde_json::json!({
            "current_sequence_id": seq_id,
            "selected_region": selected_region,
            "active_project_id": self.current_project_path.clone().or_else(|| Some("unsaved".to_string())),
            "resources": resources,
        })
    }

    fn clawbio_anchor_resource(&self, anchor: SequenceGenomeAnchorSummary) -> Value {
        serde_json::json!({
            "kind": "genome_sequence_anchor",
            "seq_id": anchor.seq_id,
            "genome_id": anchor.genome_id,
            "chromosome": anchor.chromosome,
            "start_1based": anchor.start_1based,
            "end_1based": anchor.end_1based,
            "strand": anchor.strand,
            "anchor_verified": anchor.anchor_verified,
            "catalog_path": self.genome_catalog_path_opt(),
            "cache_dir": self.genome_cache_dir_opt(),
        })
    }
}

pub(crate) fn default_clawbio_bin() -> String {
    env::var("CLAWBIO_BIN").unwrap_or_else(|_| "clawbio.py".to_string())
}

pub(crate) fn clawbio_run_root() -> PathBuf {
    if let Ok(value) = env::var("GENTLE_CLAWBIO_RUN_DIR") {
        return PathBuf::from(value);
    }
    env::var("HOME")
        .map(PathBuf::from)
        .unwrap_or_else(|_| env::temp_dir())
        .join(".gentle")
        .join("clawbio_runs")
}

fn next_clawbio_run_dir(root: &Path, skill_alias: &str) -> Result<PathBuf, String> {
    let safe_skill: String = skill_alias
        .chars()
        .map(|ch| if ch.is_ascii_alphanumeric() { ch } else { '_' })
        .collect();
    fs::create_dir_all(root).map_err(|e| format!("Could not create ClawBio run root: {e}"))?;
    Ok(root.join(format!("{}_{}", safe_skill, GENtleApp::now_unix_ms())))
}

fn resolve_result_path(output_dir: &Path, path: &str) -> String {
    let path = PathBuf::from(path);
    if path.is_absolute() {
        path.display().to_string()
    } else {
        output_dir.join(path).display().to_string()
    }
}

fn string_list_field(value: &Value, key: &str) -> Option<Vec<String>> {
    value.get(key)?.as_array().map(|rows| {
        rows.iter()
            .filter_map(Value::as_str)
            .map(str::to_string)
            .collect()
    })
}

fn value_array_field(value: &Value, key: &str) -> Option<Vec<Value>> {
    value.get(key)?.as_array().cloned()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn clawbio_result_view_extracts_summary_actions_and_artifacts() {
        let temp = tempfile::tempdir().expect("tempdir");
        fs::write(
            temp.path().join("result.json"),
            serde_json::json!({
                "status": "ok",
                "workflow_state": {"lifecycle": "ready", "state_label": "Demo"},
                "chat_summary_lines": ["one", "two"],
                "preferred_artifacts": [{"artifact_id": "fig", "path": "figure.svg"}],
                "suggested_actions": [{"label": "Continue", "request": {"mode": "skill-info"}}],
                "artifacts": {"report_md": "report.md"}
            })
            .to_string(),
        )
        .expect("write result");

        let view = ClawBioResultView::read_from_output_dir(temp.path()).expect("parse result");
        assert_eq!(view.chat_summary_lines, vec!["one", "two"]);
        assert_eq!(view.preferred_artifacts.len(), 1);
        assert_eq!(view.suggested_actions.len(), 1);
        assert!(view.report_md_path.as_ref().unwrap().ends_with("report.md"));
        assert_eq!(view.workflow_state.as_ref().unwrap()["lifecycle"], "ready");
    }

    #[test]
    fn empty_clawbio_context_payload_is_valid() {
        let app = GENtleApp::default();
        let context = app.clawbio_context_payload();
        assert!(context["current_sequence_id"].is_null());
        assert!(context["selected_region"].is_null());
        assert_eq!(context["active_project_id"], "unsaved");
        assert!(context["resources"].as_array().unwrap().is_empty());
    }

    #[test]
    fn cancel_token_reports_cancelled_state() {
        let token = CancelToken::default();
        assert!(!token.is_cancelled());
        token.cancel();
        assert!(token.is_cancelled());
    }
}
