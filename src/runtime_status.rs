//! Process-local runtime activity registry and SIGUSR1 diagnostics.
//!
//! The registry is intentionally in-memory only. It is a live activity stack,
//! not a persisted job log: significant routines push a frame while they run,
//! update phase/progress fields as work advances, and drop the frame when the
//! work leaves scope.

use serde::{Deserialize, Serialize};
use std::{
    cell::RefCell,
    collections::BTreeMap,
    process,
    sync::{
        LazyLock, Mutex,
        atomic::{AtomicBool, AtomicU64, Ordering},
    },
    thread,
    time::{SystemTime, UNIX_EPOCH},
};

/// Schema identifier for the process-local runtime status snapshot.
pub const RUNTIME_STATUS_SCHEMA: &str = "gentle.runtime_status.v1";

thread_local! {
    static RUNTIME_FRAME_STACK: RefCell<Vec<String>> = RefCell::new(Vec::new());
}

static RUNTIME_STATUS_REGISTRY: LazyLock<RuntimeStatusRegistry> =
    LazyLock::new(RuntimeStatusRegistry::default);

#[cfg(unix)]
static SIGUSR1_DUMP_THREAD_INSTALLED: AtomicBool = AtomicBool::new(false);

/// Why a runtime-status snapshot was requested.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum RuntimeStatusTrigger {
    Manual,
    Sigusr1,
    Agent,
    Mcp,
    Shell,
    Test,
}

impl Default for RuntimeStatusTrigger {
    fn default() -> Self {
        Self::Manual
    }
}

impl RuntimeStatusTrigger {
    fn as_str(self) -> &'static str {
        match self {
            Self::Manual => "manual",
            Self::Sigusr1 => "sigusr1",
            Self::Agent => "agent",
            Self::Mcp => "mcp",
            Self::Shell => "shell",
            Self::Test => "test",
        }
    }
}

/// Broad class of a live activity frame.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum RuntimeStatusFrameKind {
    ForegroundOperation,
    BackgroundJob,
    GuiAction,
    ShellCommand,
    AgentRequest,
    McpRequest,
    ResourcePrepare,
    ExternalTool,
    Other,
}

impl Default for RuntimeStatusFrameKind {
    fn default() -> Self {
        Self::Other
    }
}

/// Current state of a live activity frame.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum RuntimeStatusFrameState {
    Running,
    Waiting,
    Completed,
    Failed,
    Cancelled,
}

impl Default for RuntimeStatusFrameState {
    fn default() -> Self {
        Self::Running
    }
}

/// Source subsystem for an observed activity record.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum RuntimeStatusActivitySource {
    PrepareGenome,
    CutrunSharedAsset,
    BlastAsyncJob,
    Other,
}

impl Default for RuntimeStatusActivitySource {
    fn default() -> Self {
        Self::Other
    }
}

/// Where the activity observation came from.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum RuntimeStatusActivityScope {
    ProcessLocal,
    PersistedActivity,
    ProjectAsyncRegistry,
}

impl Default for RuntimeStatusActivityScope {
    fn default() -> Self {
        Self::ProcessLocal
    }
}

/// Freshness/state class for an observed activity.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum RuntimeStatusActivityObservation {
    Live,
    CrossProcess,
    Completed,
    Failed,
    Cancelled,
    Stale,
    Unknown,
}

impl Default for RuntimeStatusActivityObservation {
    fn default() -> Self {
        Self::Unknown
    }
}

/// One active activity frame in the process-local stack.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(default)]
pub struct RuntimeStatusFrame {
    pub frame_id: String,
    pub parent_frame_id: Option<String>,
    pub kind: RuntimeStatusFrameKind,
    pub label: String,
    pub phase: Option<String>,
    pub detail: Option<String>,
    pub state: RuntimeStatusFrameState,
    pub started_at_unix_ms: u128,
    pub updated_at_unix_ms: u128,
    pub progress_percent: Option<f64>,
    pub bytes_done: Option<u64>,
    pub bytes_total: Option<u64>,
    pub thread: Option<String>,
}

impl Default for RuntimeStatusFrame {
    fn default() -> Self {
        Self {
            frame_id: String::new(),
            parent_frame_id: None,
            kind: RuntimeStatusFrameKind::Other,
            label: String::new(),
            phase: None,
            detail: None,
            state: RuntimeStatusFrameState::Running,
            started_at_unix_ms: 0,
            updated_at_unix_ms: 0,
            progress_percent: None,
            bytes_done: None,
            bytes_total: None,
            thread: None,
        }
    }
}

/// One observed activity from an existing persisted or project-backed ledger.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(default)]
pub struct RuntimeStatusActivity {
    pub activity_id: String,
    pub source: RuntimeStatusActivitySource,
    pub scope: RuntimeStatusActivityScope,
    pub label: String,
    pub lifecycle_status: String,
    pub observation: RuntimeStatusActivityObservation,
    pub phase: Option<String>,
    pub detail: Option<String>,
    pub progress_percent: Option<f64>,
    pub bytes_done: Option<u64>,
    pub bytes_total: Option<u64>,
    pub started_at_unix_ms: Option<u128>,
    pub updated_at_unix_ms: Option<u128>,
    pub finished_at_unix_ms: Option<u128>,
    pub origin_process_id: Option<u32>,
    pub source_path: Option<String>,
    pub stale_reason: Option<String>,
}

impl RuntimeStatusActivity {
    pub fn new(
        source: RuntimeStatusActivitySource,
        scope: RuntimeStatusActivityScope,
        activity_id: impl Into<String>,
        label: impl Into<String>,
        lifecycle_status: impl Into<String>,
    ) -> Self {
        let lifecycle_status = lifecycle_status.into();
        Self {
            activity_id: activity_id.into(),
            source,
            scope,
            label: label.into(),
            observation: classify_activity_observation(&lifecycle_status),
            lifecycle_status,
            ..Self::default()
        }
    }
}

impl Default for RuntimeStatusActivity {
    fn default() -> Self {
        Self {
            activity_id: String::new(),
            source: RuntimeStatusActivitySource::Other,
            scope: RuntimeStatusActivityScope::ProcessLocal,
            label: String::new(),
            lifecycle_status: String::new(),
            observation: RuntimeStatusActivityObservation::Unknown,
            phase: None,
            detail: None,
            progress_percent: None,
            bytes_done: None,
            bytes_total: None,
            started_at_unix_ms: None,
            updated_at_unix_ms: None,
            finished_at_unix_ms: None,
            origin_process_id: None,
            source_path: None,
            stale_reason: None,
        }
    }
}

/// Snapshot of the process-local runtime registry.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(default)]
pub struct RuntimeStatusSnapshot {
    pub schema: String,
    pub generated_at_unix_ms: u128,
    pub process_id: u32,
    pub process_label: String,
    pub trigger: RuntimeStatusTrigger,
    pub frames: Vec<RuntimeStatusFrame>,
    pub activities: Vec<RuntimeStatusActivity>,
    pub observability_notes: Vec<String>,
    pub summary_lines: Vec<String>,
}

impl Default for RuntimeStatusSnapshot {
    fn default() -> Self {
        Self {
            schema: RUNTIME_STATUS_SCHEMA.to_string(),
            generated_at_unix_ms: 0,
            process_id: process::id(),
            process_label: String::new(),
            trigger: RuntimeStatusTrigger::Manual,
            frames: Vec::new(),
            activities: Vec::new(),
            observability_notes: Vec::new(),
            summary_lines: Vec::new(),
        }
    }
}

impl RuntimeStatusSnapshot {
    /// Attach observed activity records gathered from existing ledgers.
    pub fn with_activities(
        mut self,
        mut activities: Vec<RuntimeStatusActivity>,
        observability_notes: Vec<String>,
    ) -> Self {
        activities.sort_by(|left, right| {
            format!("{:?}", left.source)
                .cmp(&format!("{:?}", right.source))
                .then(left.activity_id.cmp(&right.activity_id))
                .then(left.label.cmp(&right.label))
        });
        self.activities = activities;
        self.observability_notes = observability_notes;
        self.summary_lines = runtime_summary_lines(&self.frames, &self.activities);
        self
    }
}

/// Process-local active activity registry.
#[derive(Debug)]
pub struct RuntimeStatusRegistry {
    frames: Mutex<BTreeMap<String, RuntimeStatusFrame>>,
    counter: AtomicU64,
}

impl Default for RuntimeStatusRegistry {
    fn default() -> Self {
        Self {
            frames: Mutex::new(BTreeMap::new()),
            counter: AtomicU64::new(0),
        }
    }
}

impl RuntimeStatusRegistry {
    /// Push one active frame and return a guard that removes it on drop.
    pub fn push(
        &'static self,
        kind: RuntimeStatusFrameKind,
        label: impl Into<String>,
    ) -> RuntimeStatusGuard {
        self.push_with_detail(kind, label, None)
    }

    /// Push one active frame with an initial detail string.
    pub fn push_with_detail(
        &'static self,
        kind: RuntimeStatusFrameKind,
        label: impl Into<String>,
        detail: Option<String>,
    ) -> RuntimeStatusGuard {
        let sequence = self.counter.fetch_add(1, Ordering::Relaxed) + 1;
        let frame_id = format!("runtime-frame-{sequence}");
        let parent_frame_id = current_parent_frame_id();
        let now = now_unix_ms();
        let frame = RuntimeStatusFrame {
            frame_id: frame_id.clone(),
            parent_frame_id,
            kind,
            label: label.into(),
            phase: None,
            detail,
            state: RuntimeStatusFrameState::Running,
            started_at_unix_ms: now,
            updated_at_unix_ms: now,
            progress_percent: None,
            bytes_done: None,
            bytes_total: None,
            thread: thread::current().name().map(ToString::to_string),
        };
        if let Ok(mut frames) = self.frames.lock() {
            frames.insert(frame_id.clone(), frame);
        }
        RUNTIME_FRAME_STACK.with(|stack| stack.borrow_mut().push(frame_id.clone()));
        RuntimeStatusGuard {
            registry: self,
            frame_id,
            active: true,
        }
    }

    fn update_frame(&self, frame_id: &str, update: impl FnOnce(&mut RuntimeStatusFrame)) {
        if let Ok(mut frames) = self.frames.lock()
            && let Some(frame) = frames.get_mut(frame_id)
        {
            update(frame);
            frame.updated_at_unix_ms = now_unix_ms();
        }
    }

    fn remove_frame(&self, frame_id: &str) {
        if let Ok(mut frames) = self.frames.lock() {
            frames.remove(frame_id);
        }
    }

    /// Return a structured snapshot of all active frames in this process.
    pub fn snapshot(&self, trigger: RuntimeStatusTrigger) -> RuntimeStatusSnapshot {
        let mut frames = self
            .frames
            .lock()
            .map(|frames| frames.values().cloned().collect::<Vec<_>>())
            .unwrap_or_default();
        frames.sort_by(|left, right| {
            left.started_at_unix_ms
                .cmp(&right.started_at_unix_ms)
                .then_with(|| left.frame_id.cmp(&right.frame_id))
        });
        RuntimeStatusSnapshot {
            schema: RUNTIME_STATUS_SCHEMA.to_string(),
            generated_at_unix_ms: now_unix_ms(),
            process_id: process::id(),
            process_label: process_label(),
            trigger,
            summary_lines: runtime_summary_lines(&frames, &[]),
            frames,
            activities: Vec::new(),
            observability_notes: Vec::new(),
        }
    }
}

/// RAII guard for one active runtime status frame.
#[derive(Debug)]
pub struct RuntimeStatusGuard {
    registry: &'static RuntimeStatusRegistry,
    frame_id: String,
    active: bool,
}

impl RuntimeStatusGuard {
    pub fn update_phase(&self, phase: impl Into<String>) {
        let phase = phase.into();
        self.registry.update_frame(&self.frame_id, |frame| {
            frame.phase = Some(phase);
        });
    }

    pub fn update_detail(&self, detail: impl Into<String>) {
        let detail = detail.into();
        self.registry.update_frame(&self.frame_id, |frame| {
            frame.detail = Some(detail);
        });
    }

    pub fn update_state(&self, state: RuntimeStatusFrameState) {
        self.registry.update_frame(&self.frame_id, |frame| {
            frame.state = state;
        });
    }

    pub fn update_progress(&self, progress_percent: Option<f64>) {
        let progress_percent = progress_percent.map(|value| value.clamp(0.0, 100.0));
        self.registry.update_frame(&self.frame_id, |frame| {
            frame.progress_percent = progress_percent;
        });
    }

    pub fn update_bytes(&self, bytes_done: u64, bytes_total: Option<u64>) {
        self.registry.update_frame(&self.frame_id, |frame| {
            frame.bytes_done = Some(bytes_done);
            frame.bytes_total = bytes_total;
        });
    }

    pub fn update_from_progress(
        &self,
        phase: impl Into<String>,
        detail: impl Into<String>,
        bytes_done: u64,
        bytes_total: Option<u64>,
        progress_percent: Option<f64>,
    ) {
        let phase = phase.into();
        let detail = detail.into();
        let progress_percent = progress_percent.map(|value| value.clamp(0.0, 100.0));
        self.registry.update_frame(&self.frame_id, |frame| {
            frame.phase = Some(phase);
            frame.detail = Some(detail);
            frame.bytes_done = Some(bytes_done);
            frame.bytes_total = bytes_total;
            frame.progress_percent = progress_percent;
            frame.state = RuntimeStatusFrameState::Running;
        });
    }

    pub fn fail(&self, message: impl Into<String>) {
        let message = message.into();
        self.registry.update_frame(&self.frame_id, |frame| {
            frame.state = RuntimeStatusFrameState::Failed;
            frame.detail = Some(message);
        });
    }

    pub fn cancel(&self, message: impl Into<String>) {
        let message = message.into();
        self.registry.update_frame(&self.frame_id, |frame| {
            frame.state = RuntimeStatusFrameState::Cancelled;
            frame.detail = Some(message);
        });
    }

    pub fn frame_id(&self) -> &str {
        &self.frame_id
    }
}

impl Drop for RuntimeStatusGuard {
    fn drop(&mut self) {
        if !self.active {
            return;
        }
        self.registry.remove_frame(&self.frame_id);
        RUNTIME_FRAME_STACK.with(|stack| {
            let mut stack = stack.borrow_mut();
            if let Some(pos) = stack.iter().rposition(|id| id == &self.frame_id) {
                stack.remove(pos);
            }
        });
        self.active = false;
    }
}

/// Return the global process-local runtime status registry.
pub fn runtime_status_registry() -> &'static RuntimeStatusRegistry {
    &RUNTIME_STATUS_REGISTRY
}

/// Return a structured snapshot of active runtime frames in this process.
pub fn runtime_status_snapshot(trigger: RuntimeStatusTrigger) -> RuntimeStatusSnapshot {
    runtime_status_registry().snapshot(trigger)
}

/// Render a human-readable runtime status snapshot for stderr diagnostics.
pub fn runtime_status_snapshot_text(trigger: RuntimeStatusTrigger) -> String {
    let snapshot = runtime_status_snapshot(trigger);
    let mut out = String::new();
    out.push_str(&format!(
        "GENtle runtime status ({}) pid={} process={} frames={}\n",
        snapshot.trigger.as_str(),
        snapshot.process_id,
        snapshot.process_label,
        snapshot.frames.len()
    ));
    for line in &snapshot.summary_lines {
        out.push_str("- ");
        out.push_str(line);
        out.push('\n');
    }
    for note in &snapshot.observability_notes {
        out.push_str("- note: ");
        out.push_str(note);
        out.push('\n');
    }
    if snapshot.frames.is_empty() {
        out.push_str("- no active runtime frames in this process\n");
    } else {
        for frame in &snapshot.frames {
            out.push_str(&format!(
                "- {} {:?} state={:?}",
                frame.frame_id, frame.kind, frame.state
            ));
            out.push_str(&format!(" label={}", frame.label));
            if let Some(phase) = &frame.phase {
                out.push_str(&format!(" phase={phase}"));
            }
            if let Some(progress) = frame.progress_percent {
                out.push_str(&format!(" progress={progress:.1}%"));
            }
            if let Some(bytes_done) = frame.bytes_done {
                out.push_str(&format!(" bytes_done={bytes_done}"));
            }
            if let Some(bytes_total) = frame.bytes_total {
                out.push_str(&format!(" bytes_total={bytes_total}"));
            }
            if let Some(detail) = &frame.detail {
                out.push_str(&format!(" detail={detail}"));
            }
            if let Some(parent) = &frame.parent_frame_id {
                out.push_str(&format!(" parent={parent}"));
            }
            if let Some(thread) = &frame.thread {
                out.push_str(&format!(" thread={thread}"));
            }
            out.push('\n');
        }
    }
    if snapshot.activities.is_empty() {
        out.push_str("- no observed persisted/project activity records\n");
    } else {
        for activity in &snapshot.activities {
            out.push_str(&format!(
                "- activity {} {:?}/{:?} observation={:?}",
                activity.activity_id, activity.source, activity.scope, activity.observation
            ));
            out.push_str(&format!(" status={}", activity.lifecycle_status));
            out.push_str(&format!(" label={}", activity.label));
            if let Some(phase) = &activity.phase {
                out.push_str(&format!(" phase={phase}"));
            }
            if let Some(progress) = activity.progress_percent {
                out.push_str(&format!(" progress={progress:.1}%"));
            }
            if let Some(bytes_done) = activity.bytes_done {
                out.push_str(&format!(" bytes_done={bytes_done}"));
            }
            if let Some(bytes_total) = activity.bytes_total {
                out.push_str(&format!(" bytes_total={bytes_total}"));
            }
            if let Some(detail) = &activity.detail {
                out.push_str(&format!(" detail={detail}"));
            }
            if let Some(pid) = activity.origin_process_id {
                out.push_str(&format!(" owner_pid={pid}"));
            }
            if let Some(path) = &activity.source_path {
                out.push_str(&format!(" source_path={path}"));
            }
            if let Some(reason) = &activity.stale_reason {
                out.push_str(&format!(" stale_reason={reason}"));
            }
            out.push('\n');
        }
    }
    out
}

/// Install a Unix SIGUSR1 waiter thread that prints the runtime snapshot to stderr.
///
/// The implementation uses `sigwait` instead of an async signal handler so it
/// can safely lock the registry and format JSON/prose without doing unsafe work
/// inside the signal handler itself.
#[cfg(unix)]
pub fn install_sigusr1_stderr_dump_thread() -> Result<bool, String> {
    if SIGUSR1_DUMP_THREAD_INSTALLED.swap(true, Ordering::SeqCst) {
        return Ok(false);
    }
    let set = make_sigusr1_set()?;
    let mask_rc = unsafe { libc::pthread_sigmask(libc::SIG_BLOCK, &set, std::ptr::null_mut()) };
    if mask_rc != 0 {
        SIGUSR1_DUMP_THREAD_INSTALLED.store(false, Ordering::SeqCst);
        return Err(format!(
            "Could not block SIGUSR1 before installing runtime-status waiter: {}",
            std::io::Error::from_raw_os_error(mask_rc)
        ));
    }
    thread::Builder::new()
        .name("gentle-sigusr1-runtime-status".to_string())
        .spawn(|| {
            let Ok(set) = make_sigusr1_set() else {
                return;
            };
            loop {
                let mut signal = 0;
                let rc = unsafe { libc::sigwait(&set, &mut signal) };
                if rc == 0 && signal == libc::SIGUSR1 {
                    eprintln!(
                        "{}",
                        runtime_status_snapshot_text(RuntimeStatusTrigger::Sigusr1)
                    );
                } else if rc != 0 {
                    eprintln!(
                        "GENtle runtime status SIGUSR1 waiter failed: {}",
                        std::io::Error::from_raw_os_error(rc)
                    );
                }
            }
        })
        .map_err(|error| {
            SIGUSR1_DUMP_THREAD_INSTALLED.store(false, Ordering::SeqCst);
            format!("Could not start SIGUSR1 runtime-status waiter thread: {error}")
        })?;
    Ok(true)
}

/// Non-Unix platforms do not have SIGUSR1; the read APIs remain available.
#[cfg(not(unix))]
pub fn install_sigusr1_stderr_dump_thread() -> Result<bool, String> {
    Ok(false)
}

#[cfg(unix)]
fn make_sigusr1_set() -> Result<libc::sigset_t, String> {
    let mut set = unsafe { std::mem::zeroed::<libc::sigset_t>() };
    let empty_rc = unsafe { libc::sigemptyset(&mut set) };
    if empty_rc != 0 {
        return Err(format!(
            "Could not initialize signal set: {}",
            std::io::Error::last_os_error()
        ));
    }
    let add_rc = unsafe { libc::sigaddset(&mut set, libc::SIGUSR1) };
    if add_rc != 0 {
        return Err(format!(
            "Could not add SIGUSR1 to signal set: {}",
            std::io::Error::last_os_error()
        ));
    }
    Ok(set)
}

fn current_parent_frame_id() -> Option<String> {
    RUNTIME_FRAME_STACK.with(|stack| stack.borrow().last().cloned())
}

fn runtime_summary_lines(
    frames: &[RuntimeStatusFrame],
    activities: &[RuntimeStatusActivity],
) -> Vec<String> {
    let mut lines = Vec::new();
    if frames.is_empty() {
        lines.push("No active runtime frames in this process.".to_string());
    } else {
        let mut by_kind = BTreeMap::<String, usize>::new();
        for frame in frames {
            *by_kind.entry(format!("{:?}", frame.kind)).or_default() += 1;
        }
        let kind_summary = by_kind
            .into_iter()
            .map(|(kind, count)| format!("{kind}={count}"))
            .collect::<Vec<_>>()
            .join(", ");
        lines.push(format!(
            "{} active runtime frame(s): {}",
            frames.len(),
            kind_summary
        ));
    }
    if activities.is_empty() {
        lines.push("No observed persisted/project activity records.".to_string());
    } else {
        let mut by_source = BTreeMap::<String, usize>::new();
        for activity in activities {
            *by_source
                .entry(format!("{:?}", activity.source))
                .or_default() += 1;
        }
        let activity_summary = by_source
            .into_iter()
            .map(|(source, count)| format!("{source}={count}"))
            .collect::<Vec<_>>()
            .join(", ");
        lines.push(format!(
            "{} observed activity record(s): {}",
            activities.len(),
            activity_summary
        ));
    }
    lines
}

fn classify_activity_observation(lifecycle_status: &str) -> RuntimeStatusActivityObservation {
    match lifecycle_status.trim().to_ascii_lowercase().as_str() {
        "running" | "queued" | "preparing" | "in_progress" => RuntimeStatusActivityObservation::Live,
        "ready" | "completed" | "done" | "ok" => RuntimeStatusActivityObservation::Completed,
        "failed" | "error" => RuntimeStatusActivityObservation::Failed,
        "cancelled" | "canceled" => RuntimeStatusActivityObservation::Cancelled,
        "stale" => RuntimeStatusActivityObservation::Stale,
        _ => RuntimeStatusActivityObservation::Unknown,
    }
}

fn process_label() -> String {
    std::env::args()
        .next()
        .filter(|value| !value.trim().is_empty())
        .unwrap_or_else(|| "gentle".to_string())
}

fn now_unix_ms() -> u128 {
    SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_millis())
        .unwrap_or(0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn runtime_status_snapshot_reports_nested_active_frames() {
        let outer = runtime_status_registry().push_with_detail(
            RuntimeStatusFrameKind::ShellCommand,
            "runtime-test outer shell",
            Some("state-summary".to_string()),
        );
        let inner = runtime_status_registry().push(
            RuntimeStatusFrameKind::ResourcePrepare,
            "runtime-test GRCh38",
        );
        inner.update_from_progress(
            "index_transcripts",
            "Transcript Index",
            40,
            Some(100),
            Some(40.0),
        );

        let snapshot = runtime_status_snapshot(RuntimeStatusTrigger::Test);

        assert_eq!(snapshot.schema, RUNTIME_STATUS_SCHEMA);
        let inner_frame = snapshot
            .frames
            .iter()
            .find(|frame| frame.label == "runtime-test GRCh38")
            .expect("inner runtime frame");
        assert_eq!(
            inner_frame.parent_frame_id.as_deref(),
            Some(outer.frame_id())
        );
        assert_eq!(inner_frame.phase.as_deref(), Some("index_transcripts"));
        assert_eq!(inner_frame.progress_percent, Some(40.0));
        drop(inner);
        drop(outer);
        let after = runtime_status_snapshot(RuntimeStatusTrigger::Test);
        assert!(
            !after
                .frames
                .iter()
                .any(|frame| frame.label == "runtime-test GRCh38")
        );
    }

    #[test]
    fn runtime_status_text_is_stderr_friendly() {
        let _guard = runtime_status_registry().push(
            RuntimeStatusFrameKind::AgentRequest,
            "runtime-test ask agent",
        );
        let text = runtime_status_snapshot_text(RuntimeStatusTrigger::Sigusr1);

        assert!(text.contains("GENtle runtime status"));
        assert!(text.contains("runtime-test ask agent"));
    }

    #[test]
    fn runtime_status_snapshot_can_attach_observed_activities_without_files() {
        let activity = RuntimeStatusActivity::new(
            RuntimeStatusActivitySource::PrepareGenome,
            RuntimeStatusActivityScope::PersistedActivity,
            "prepare-genome:reference:ToyGenome",
            "Reference genome prepare: ToyGenome",
            "running",
        );

        let snapshot = runtime_status_snapshot(RuntimeStatusTrigger::Test).with_activities(
            vec![activity],
            vec!["test observation note".to_string()],
        );

        assert_eq!(snapshot.activities.len(), 1);
        assert_eq!(
            snapshot.activities[0].observation,
            RuntimeStatusActivityObservation::Live
        );
        assert!(
            snapshot
                .summary_lines
                .iter()
                .any(|line| line.contains("observed activity record"))
        );
        let text = {
            let mut out = String::new();
            for line in &snapshot.summary_lines {
                out.push_str(line);
            }
            out
        };
        assert!(text.contains("PrepareGenome=1"));
    }
}
