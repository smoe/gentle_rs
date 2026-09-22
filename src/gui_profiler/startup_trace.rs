//! Opt-in startup attribution without per-frame file I/O or biological identifiers.
//!
//! Checkpoints measure CPU-side work from Rust `main`, not OS process launch,
//! GPU presentation, input responsiveness or independently confirmed content.
//! Events use a bounded, nonblocking recorder; losses remain explicit. Only the
//! session owner writes a new JSON file, after the application loop has returned.

use serde::Serialize;
use std::fs::OpenOptions;
use std::io::{self, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::sync::{Arc, Mutex, OnceLock};
use std::time::Instant;

/// Explicit output destination. Unset/empty means no recording and no writes.
pub const STARTUP_TRACE_ENV: &str = "GENTLE_GUI_STARTUP_TRACE";
const MAX_EVENTS: usize = 512;
static PROCESS_TRACE: OnceLock<TraceContext> = OnceLock::new();

/// Fixed vocabulary deliberately excludes user-controlled labels and paths.
#[derive(Clone, Copy, Debug, Serialize)]
#[serde(rename_all = "snake_case")]
#[repr(u8)]
pub enum Phase {
    MainEntered,
    NativeRun,
    AppInitialize,
    AppDefaults,
    HelpPreparation,
    ConfigurationLoad,
    CredentialRefresh,
    InitialStatePreparation,
    ProjectLoad,
    ProjectReadDecode,
    ProjectInstall,
    RootFirstFrame,
    RootWorkspaceFrame,
    DnaOpenDispatch,
    DnaConstruct,
    DnaPlaceholderConstruct,
    DnaWorkerScheduled,
    DnaEngineReadLock,
    DnaSequenceClone,
    DnaWorkerResult,
    DnaHydrate,
    DnaLoadFailed,
    DnaNativeContentFrame,
    DnaEmbeddedContentFrame,
}

#[derive(Clone, Copy, Debug, Serialize)]
#[serde(rename_all = "snake_case")]
enum EventKind {
    Begin,
    Completed,
    Failed,
    Interrupted,
    Checkpoint,
}

#[derive(Debug, Serialize)]
struct Event {
    elapsed_us: u64,
    subject: u64,
    phase: Phase,
    kind: EventKind,
    #[serde(skip_serializing_if = "Option::is_none")]
    span: Option<u64>,
}

#[derive(Debug)]
struct Recorder {
    origin: Instant,
    events: Mutex<Vec<Event>>,
    next_id: AtomicU64,
    dropped: AtomicU64,
    closed: AtomicBool,
}

impl Recorder {
    fn new(origin: Instant) -> Self {
        Self {
            origin,
            events: Mutex::new(Vec::with_capacity(MAX_EVENTS)),
            next_id: AtomicU64::new(1),
            dropped: AtomicU64::new(0),
            closed: AtomicBool::new(false),
        }
    }

    fn record(&self, subject: u64, phase: Phase, kind: EventKind, span: Option<u64>) {
        if self.closed.load(Ordering::Relaxed) {
            return;
        }
        // A diagnostic must never wait behind another thread on the paint path.
        if let Ok(mut events) = self.events.try_lock() {
            // The exit snapshot may have closed the recorder after the first check.
            if self.closed.load(Ordering::Relaxed) {
                return;
            }
            if events.len() < MAX_EVENTS {
                events.push(Event {
                    elapsed_us: self.origin.elapsed().as_micros().min(u64::MAX as u128) as u64,
                    subject,
                    phase,
                    kind,
                    span,
                });
                return;
            }
        }
        self.dropped.fetch_add(1, Ordering::Relaxed);
    }

    fn report(&self) -> io::Result<serde_json::Value> {
        self.closed.store(true, Ordering::Relaxed);
        let events = self
            .events
            .lock()
            .map_err(|_| io::Error::other("Startup trace recorder poisoned"))?;
        Ok(serde_json::json!({
            "schema": "gentle.gui_startup_trace.v1",
            "source_revision": crate::about::GENTLE_SOURCE_REVISION,
            "git_commit": crate::about::GENTLE_GIT_COMMIT,
            "os": std::env::consts::OS,
            "architecture": std::env::consts::ARCH,
            "debug_assertions": cfg!(debug_assertions),
            "gui_test_support": cfg!(feature = "gui-test-support"),
            "gui_profiler": cfg!(feature = "gui-profiler"),
            "clock": "monotonic microseconds since Rust main entry; excludes process loader",
            "scope": "CPU checkpoints only; not compositor presentation or native usability acceptance",
            "event_limit": MAX_EVENTS,
            "dropped_events": self.dropped.load(Ordering::Relaxed),
            "events": *events,
        }))
    }
}

#[derive(Debug)]
struct Subject {
    recorder: Arc<Recorder>,
    id: u64,
    checkpoints: AtomicU64,
}

/// Process-local ordinal identity, never a sequence ID, path or persistent hash.
#[derive(Clone, Debug, Default)]
pub(crate) struct TraceContext(Option<Arc<Subject>>);

impl TraceContext {
    fn process(recorder: Arc<Recorder>) -> Self {
        Self(Some(Arc::new(Subject {
            recorder,
            id: 0,
            checkpoints: AtomicU64::new(0),
        })))
    }

    pub(crate) fn child(&self) -> Self {
        Self(self.0.as_ref().map(|subject| {
            Arc::new(Subject {
                id: subject.recorder.next_id.fetch_add(1, Ordering::Relaxed),
                recorder: subject.recorder.clone(),
                checkpoints: AtomicU64::new(0),
            })
        }))
    }

    pub(crate) fn checkpoint(&self, phase: Phase) {
        if let Some(subject) = &self.0 {
            let mask = 1_u64 << phase as u8;
            if subject.checkpoints.fetch_or(mask, Ordering::Relaxed) & mask == 0 {
                subject
                    .recorder
                    .record(subject.id, phase, EventKind::Checkpoint, None);
            }
        }
    }

    pub(crate) fn span(&self, phase: Phase) -> Span {
        Span(self.0.as_ref().map(|subject| {
            let id = subject.recorder.next_id.fetch_add(1, Ordering::Relaxed);
            subject
                .recorder
                .record(subject.id, phase, EventKind::Begin, Some(id));
            (subject.clone(), phase, id)
        }))
    }
}

/// An unfinished or unwinding scope is recorded as interrupted, never successful.
#[must_use]
pub struct Span(Option<(Arc<Subject>, Phase, u64)>);

impl Span {
    pub fn finish(mut self, succeeded: bool) {
        if let Some((subject, phase, id)) = self.0.take() {
            subject.recorder.record(
                subject.id,
                phase,
                if succeeded {
                    EventKind::Completed
                } else {
                    EventKind::Failed
                },
                Some(id),
            );
        }
    }
}

impl Drop for Span {
    fn drop(&mut self) {
        if let Some((subject, phase, id)) = self.0.take() {
            subject
                .recorder
                .record(subject.id, phase, EventKind::Interrupted, Some(id));
        }
    }
}

pub(crate) fn context() -> TraceContext {
    #[cfg(test)]
    if let Some(context) = TEST_CONTEXT.with(|slot| slot.borrow().clone()) {
        return context;
    }
    PROCESS_TRACE.get().cloned().unwrap_or_default()
}

/// Time a named phase; inert unless recording was explicitly enabled.
pub fn span(phase: Phase) -> Span {
    context().span(phase)
}

/// Bounded process checkpoint, recorded at most once during this invocation.
pub fn checkpoint(phase: Phase) {
    context().checkpoint(phase);
}

pub(crate) fn measure<T>(phase: Phase, operation: impl FnOnce() -> T) -> T {
    let span = span(phase);
    let value = operation();
    span.finish(true);
    value
}

/// Owns the sole output write. Keep it alive until the native loop has returned.
pub struct Session {
    recorder: Arc<Recorder>,
    path: PathBuf,
}

/// Start explicitly requested recording, without opening a file on the GUI path.
///
/// The destination must not already exist; it is written only on orderly scope
/// exit (including unwinding). Forced termination may leave no trace file.
pub fn start_from_env(main_entered: Instant) -> Option<Session> {
    let path = std::env::var_os(STARTUP_TRACE_ENV).filter(|value| !value.is_empty())?;
    let recorder = Arc::new(Recorder::new(main_entered));
    let context = TraceContext::process(recorder.clone());
    if PROCESS_TRACE.set(context).is_err() {
        return None;
    }
    checkpoint(Phase::MainEntered);
    Some(Session {
        recorder,
        path: path.into(),
    })
}

fn write_report(path: &Path, report: &serde_json::Value) -> io::Result<()> {
    let bytes = serde_json::to_vec_pretty(report).map_err(io::Error::other)?;
    let mut file = OpenOptions::new().write(true).create_new(true).open(path)?;
    file.write_all(&bytes)?;
    file.write_all(b"\n")?;
    file.flush()
}

impl Drop for Session {
    fn drop(&mut self) {
        let result = self
            .recorder
            .report()
            .and_then(|report| write_report(&self.path, &report));
        if let Err(error) = result {
            eprintln!("W GENtle startup trace was not saved: {error}");
        }
    }
}

#[cfg(test)]
thread_local! {
    static TEST_CONTEXT: std::cell::RefCell<Option<TraceContext>> = const { std::cell::RefCell::new(None) };
}

/// Isolate test observations without process-environment mutation or shared state.
#[cfg(test)]
pub(crate) fn capture<T>(operation: impl FnOnce() -> T) -> (T, serde_json::Value) {
    struct Reset(Option<TraceContext>);
    impl Drop for Reset {
        fn drop(&mut self) {
            TEST_CONTEXT.with(|slot| *slot.borrow_mut() = self.0.take());
        }
    }
    let recorder = Arc::new(Recorder::new(Instant::now()));
    let reset = Reset(
        TEST_CONTEXT.with(|slot| slot.replace(Some(TraceContext::process(recorder.clone())))),
    );
    let value = operation();
    drop(reset);
    (value, recorder.report().unwrap())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn startup_trace_is_inert_without_a_context() {
        let context = TraceContext::default();
        context.checkpoint(Phase::RootFirstFrame);
        assert!(context.child().span(Phase::DnaConstruct).0.is_none());
    }

    #[test]
    fn startup_trace_binds_scopes_and_deduplicates_cloned_subjects() {
        let (_, report) = capture(|| {
            let subject = context().child();
            subject.span(Phase::DnaConstruct).finish(true);
            subject.checkpoint(Phase::DnaEmbeddedContentFrame);
            subject.clone().checkpoint(Phase::DnaEmbeddedContentFrame);
            subject.span(Phase::DnaHydrate).finish(false);
            drop(subject.span(Phase::DnaSequenceClone));
        });
        let events = report["events"].as_array().unwrap();
        assert_eq!(events.len(), 7);
        assert_eq!(events[0]["span"], events[1]["span"]);
        assert_eq!(events[1]["kind"], "completed");
        assert_eq!(events[4]["kind"], "failed");
        assert_eq!(events[6]["kind"], "interrupted");
        assert!(
            events
                .iter()
                .all(|event| event["subject"] == events[0]["subject"])
        );
        assert!(
            events
                .windows(2)
                .all(|pair| pair[0]["elapsed_us"].as_u64() <= pair[1]["elapsed_us"].as_u64())
        );
    }

    #[test]
    fn startup_trace_records_losses_instead_of_waiting_or_growing() {
        let recorder = Arc::new(Recorder::new(Instant::now()));
        let context = TraceContext::process(recorder.clone());
        let locked = recorder.events.lock().unwrap();
        context.checkpoint(Phase::RootFirstFrame);
        drop(locked);
        for _ in 0..MAX_EVENTS {
            context.span(Phase::ProjectLoad).finish(true);
        }
        let report = recorder.report().unwrap();
        assert_eq!(report["events"].as_array().unwrap().len(), MAX_EVENTS);
        assert_eq!(report["dropped_events"], MAX_EVENTS + 1);
        assert!((Phase::DnaEmbeddedContentFrame as u8) < 64);
    }

    #[test]
    fn startup_trace_transfers_subject_to_worker_without_private_names() {
        let (_, report) = capture(|| {
            let subject = context().child();
            let child = subject.clone();
            std::thread::spawn(move || child.span(Phase::DnaSequenceClone).finish(true))
                .join()
                .unwrap();
            subject.checkpoint(Phase::DnaEmbeddedContentFrame);
        });
        let events = report["events"].as_array().unwrap();
        assert!(
            events
                .iter()
                .all(|event| event["subject"] == events[0]["subject"])
        );
        for event in events {
            assert!(event.as_object().unwrap().keys().all(|key| {
                ["elapsed_us", "subject", "phase", "kind", "span"].contains(&key.as_str())
            }));
        }
    }

    #[test]
    fn startup_trace_exit_snapshot_stops_late_worker_events() {
        let recorder = Arc::new(Recorder::new(Instant::now()));
        let context = TraceContext::process(recorder.clone());
        let unfinished = context.span(Phase::DnaHydrate);
        let snapshot = recorder.report().unwrap();
        std::thread::spawn(move || {
            unfinished.finish(true);
            context.checkpoint(Phase::DnaWorkerResult);
        })
        .join()
        .unwrap();
        assert_eq!(recorder.report().unwrap(), snapshot);
        assert_eq!(snapshot["events"].as_array().unwrap().len(), 1);
        assert_eq!(snapshot["events"][0]["kind"], "begin");
    }

    #[test]
    fn startup_trace_output_is_new_json_and_never_overwrites_evidence() {
        let temp = tempfile::tempdir().unwrap();
        let path = temp.path().join("startup trace.json");
        let (_, report) = capture(|| checkpoint(Phase::MainEntered));
        write_report(&path, &report).unwrap();
        let bytes = std::fs::read(&path).unwrap();
        assert!(bytes.ends_with(b"\n"));
        assert_eq!(
            serde_json::from_slice::<serde_json::Value>(&bytes).unwrap(),
            report
        );
        assert!(write_report(&path, &serde_json::json!({"replacement":true})).is_err());
        assert_eq!(std::fs::read(&path).unwrap(), bytes);
        assert!(write_report(&temp.path().join("missing/trace.json"), &report).is_err());
    }
}
