//! Request-bound cancellation markers and a deterministic, shared polling throttle.

use super::*;
use std::time::{Duration, Instant};

#[derive(Default)]
pub(super) struct PollGate(Option<Instant>);

impl PollGate {
    pub(super) fn due(&mut self, now: Instant, force: bool) -> bool {
        if !force
            && self
                .0
                .is_some_and(|last| now.duration_since(last) < Duration::from_millis(100))
        {
            return false;
        }
        self.0 = Some(now);
        true
    }
}

pub(super) fn cancel_identity(request: &ScreenRnaReadPanelRequest) -> Result<String, EngineError> {
    Ok(sha256_prefixed_bytes(&encode(request)?))
}

pub(super) fn diagnostic(
    error: EngineError,
    marker_path: Option<&str>,
    phase: RnaReadPanelControlErrorPhase,
    file_id: Option<&str>,
) -> RnaReadPanelControlError {
    RnaReadPanelControlError {
        code: "panel_control_marker_error".into(),
        message: error.message,
        // marker_requested can only fail when a path was supplied.
        marker_path: marker_path.unwrap_or_default().into(),
        phase,
        file_id: file_id.map(str::to_string),
    }
}

/// A matching pre-start request remains a valid cancellation; stale/malformed markers fail closed.
pub(super) fn marker_requested(
    path: Option<&str>,
    identities: &[String],
) -> Result<bool, EngineError> {
    let Some(path) = path else {
        return Ok(false);
    };
    let metadata = match fs::symlink_metadata(path) {
        Ok(metadata) => metadata,
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => return Ok(false),
        Err(error) => return Err(bad(format!("panel_cancel_marker: {error}"))),
    };
    if !metadata.is_file() || metadata.file_type().is_symlink() || metadata.len() > 128 {
        return Err(bad(
            "panel_cancel_marker_invalid: expected a bounded regular request-bound marker",
        ));
    }
    let mut marker = String::new();
    fs::File::open(path)
        .map_err(|e| bad(format!("panel_cancel_marker: {e}")))?
        .take(129)
        .read_to_string(&mut marker)
        .map_err(|e| bad(format!("panel_cancel_marker: {e}")))?;
    if marker.len() > 128 || !identities.iter().any(|id| id == marker.trim()) {
        return Err(bad(format!(
            "panel_cancel_marker_mismatch: {path}; use a fresh cancel path for this request"
        )));
    }
    Ok(true)
}
