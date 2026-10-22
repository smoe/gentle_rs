//! Explicit bounded verification of allowlisted summary tables, never original reads.

use super::*;
use ring::digest::{Context, SHA256};

const MAX_TABLE_BYTES: u64 = 256_000_000;

impl GentleEngine {
    /// Verify the exact inspected summary's table receipts; metadata-only display never calls this.
    pub fn verify_rna_read_panel_summary_artifacts(
        &self,
        path: &str,
        summary: &RnaReadPanelSummary,
    ) -> Result<RnaReadPanelArtifactVerification, EngineError> {
        let mut basis = summary.clone();
        basis.content_sha256.clear();
        let third = match summary.schema.as_str() {
            RNA_PANEL_SUMMARY_SCHEMA => "artifacts/co_nomination_sets.tsv",
            RNA_PANEL_SUMMARY_LEGACY_SCHEMA => "artifacts/ambiguous_sets.tsv",
            _ => return Err(bad("panel_verify_requires_summary")),
        };
        if summary.content_sha256 != sha256_prefixed_bytes(&encode(&basis)?) {
            return Err(bad("panel_summary_binding"));
        }
        let expected =
            BTreeSet::from(["artifacts/samples.tsv", "artifacts/transcripts.tsv", third]);
        if summary.artifacts.len() != 3
            || summary
                .artifacts
                .iter()
                .map(|a| a.path.as_str())
                .collect::<BTreeSet<_>>()
                != expected
            || summary.artifacts.iter().any(|a| a.locus_id.is_some())
        {
            return Err(bad("panel_verify_artifact_allowlist"));
        }
        let mut remaining = MAX_TABLE_BYTES;
        for artifact in &summary.artifacts {
            remaining = remaining
                .checked_sub(artifact.bytes)
                .ok_or_else(|| bad("panel_verify_byte_budget: at most 256 MB of summary tables"))?;
        }
        let parent = Path::new(path)
            .parent()
            .filter(|p| !p.as_os_str().is_empty())
            .unwrap_or(Path::new("."));
        let metadata = fs::symlink_metadata(parent.join("artifacts"))
            .map_err(|e| bad(format!("panel_verify_directory: {e}")))?;
        if !metadata.is_dir() || metadata.file_type().is_symlink() {
            return Err(bad("panel_verify_directory_not_regular"));
        }
        let root = fs::canonicalize(parent).map_err(|e| bad(e.to_string()))?;
        let mut verified_bytes = 0;
        for artifact in &summary.artifacts {
            let target = parent.join(&artifact.path);
            let before = fs::symlink_metadata(&target)
                .map_err(|e| bad(format!("panel_verify_open: {e}")))?;
            if !before.is_file()
                || before.file_type().is_symlink()
                || before.len() != artifact.bytes
                || !fs::canonicalize(&target)
                    .map_err(|e| bad(e.to_string()))?
                    .starts_with(&root)
            {
                return Err(bad("panel_verify_file_binding"));
            }
            let file = fs::File::open(&target).map_err(|e| bad(e.to_string()))?;
            let opened = file.metadata().map_err(|e| bad(e.to_string()))?;
            let mut reader = file.take(artifact.bytes + 1);
            let mut hasher = Context::new(&SHA256);
            let mut bytes = 0u64;
            let mut buffer = [0u8; 65_536];
            loop {
                let count = reader
                    .read(&mut buffer)
                    .map_err(|e| bad(format!("panel_verify_read: {e}")))?;
                if count == 0 {
                    break;
                }
                bytes += count as u64;
                hasher.update(&buffer[..count]);
            }
            let after = reader
                .get_ref()
                .metadata()
                .map_err(|e| bad(e.to_string()))?;
            drop(reader);
            let final_path = fs::symlink_metadata(&target).map_err(|e| bad(e.to_string()))?;
            if bytes != artifact.bytes
                || opened.len() != bytes
                || after.len() != bytes
                || opened.modified().ok() != after.modified().ok()
                || final_path.file_type().is_symlink()
                || !final_path.is_file()
                || final_path.len() != bytes
                || final_path.modified().ok() != after.modified().ok()
                || format!("sha256:{}", input::hex(hasher.finish().as_ref())) != artifact.sha256
            {
                return Err(bad("panel_verify_artifact_hash_or_size"));
            }
            verified_bytes += bytes;
        }
        Ok(RnaReadPanelArtifactVerification {
            schema: "gentle.rna_read_panel_artifact_verification.v1".into(),
            summary_sha256: summary.content_sha256.clone(), artifacts: summary.artifacts.clone(), verified_bytes,
            original_inputs_reverified: false,
            scope: "Exact bytes of three summary TSVs only; not original FASTQ, nomination sidecars, authenticated provenance or a scientific reassessment".into(),
        })
    }
}
