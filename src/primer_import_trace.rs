//! Opt-in, bounded wall-clock attribution for synchronous primer imports.
//! Fixed phase names and numeric counts exclude paths, sequences and biological
//! identifiers. Events go only to stderr, never reports, fingerprints or state.

use serde::Serialize;
use std::cell::{Cell, RefCell};
use std::io::Write;
use std::rc::Rc;
use std::sync::atomic::{AtomicU64, Ordering};
use std::time::Instant;

const ENV: &str = "GENTLE_PRIMER_IMPORT_TRACE";
const MAX_EVENTS: usize = 512;
static NEXT_IMPORT: AtomicU64 = AtomicU64::new(1);

#[derive(Clone, Copy, Serialize)]
#[serde(rename_all = "snake_case")]
pub(crate) enum Phase {
    Import,
    Admission,
    Reference,
    DatabaseValidation,
    DatabaseInspection,
    AnnotationHash,
    AnnotationRead,
    AnnotationDecode,
    AnnotationLookup,
    OutputValidation,
    TabularParse,
    Interpretation,
    SubjectAliases,
    ForwardHits,
    ReverseHits,
    AmpliconJoin,
    Persistence,
}

#[derive(Clone, Copy, Serialize)]
#[serde(rename_all = "snake_case")]
enum Kind {
    Begin,
    Completed,
    Failed,
    Abandoned,
    Truncated,
}

#[derive(Serialize)]
struct Event {
    schema: &'static str,
    import: u64,
    span: u64,
    phase: Phase,
    kind: Kind,
    elapsed_us: u64,
    #[serde(skip_serializing_if = "Option::is_none")]
    duration_us: Option<u64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    count: Option<usize>,
}

struct Recorder {
    origin: Instant,
    import: u64,
    next_span: Cell<u64>,
    emitted: Cell<usize>,
    #[cfg(test)]
    captured: Option<Rc<RefCell<Vec<serde_json::Value>>>>,
}

thread_local! {
    static ACTIVE: RefCell<Option<Rc<Recorder>>> = const { RefCell::new(None) };
}

impl Recorder {
    fn new() -> Self {
        Self {
            origin: Instant::now(),
            import: NEXT_IMPORT.fetch_add(1, Ordering::Relaxed),
            next_span: Cell::new(1),
            emitted: Cell::new(0),
            #[cfg(test)]
            captured: None,
        }
    }

    fn emit(&self, mut event: Event) {
        let emitted = self.emitted.get();
        if emitted == MAX_EVENTS {
            return;
        }
        if emitted == MAX_EVENTS - 1 {
            event.kind = Kind::Truncated;
            event.duration_us = None;
            event.count = None;
        }
        self.emitted.set(emitted + 1);
        #[cfg(test)]
        if let Some(captured) = &self.captured {
            captured
                .borrow_mut()
                .push(serde_json::to_value(event).expect("fixed trace event serializes"));
            return;
        }
        // One lock per bounded event; sink failure must not fail the import.
        write_event(&mut std::io::stderr().lock(), &event);
    }

    fn event(
        &self,
        span: u64,
        phase: Phase,
        kind: Kind,
        started: Option<Instant>,
        count: Option<usize>,
    ) {
        self.emit(Event {
            schema: "gentle.primer_import_trace.v1",
            import: self.import,
            span,
            phase,
            kind,
            elapsed_us: micros(self.origin),
            duration_us: started.map(micros),
            count,
        });
    }
}

fn micros(started: Instant) -> u64 {
    started.elapsed().as_micros().min(u64::MAX as u128) as u64
}

fn write_event(writer: &mut impl Write, event: &Event) {
    if let Ok(bytes) = serde_json::to_vec(event) {
        let _ = writer.write_all(&bytes);
        let _ = writer.write_all(b"\n");
    }
}

pub(crate) struct Session {
    previous: Option<Rc<Recorder>>,
    active: bool,
}

impl Session {
    pub(crate) fn from_env() -> Self {
        Self::from_setting(std::env::var(ENV).ok().as_deref())
    }

    fn from_setting(setting: Option<&str>) -> Self {
        if ACTIVE.with(|active| active.borrow().is_some()) || setting != Some("1") {
            return Self {
                previous: None,
                active: false,
            };
        }
        Self::activate(Recorder::new())
    }

    fn activate(recorder: Recorder) -> Self {
        Self {
            previous: ACTIVE.with(|active| active.replace(Some(Rc::new(recorder)))),
            active: true,
        }
    }
}

impl Drop for Session {
    fn drop(&mut self) {
        if self.active {
            let _ = ACTIVE.with(|active| active.replace(self.previous.take()));
        }
    }
}

struct ActiveSpan {
    recorder: Rc<Recorder>,
    id: u64,
    phase: Phase,
    started: Instant,
}

pub(crate) struct Span(Option<ActiveSpan>);

pub(crate) fn span(phase: Phase) -> Span {
    Span(ACTIVE.with(|active| {
        let recorder = active.borrow().clone()?;
        let id = recorder.next_span.get();
        recorder.next_span.set(id.saturating_add(1));
        let started = Instant::now();
        recorder.event(id, phase, Kind::Begin, None, None);
        Some(ActiveSpan {
            recorder,
            id,
            phase,
            started,
        })
    }))
}

impl Span {
    pub(crate) fn finish(mut self, success: bool, count: Option<usize>) {
        self.end(
            if success {
                Kind::Completed
            } else {
                Kind::Failed
            },
            count,
        );
    }

    fn end(&mut self, kind: Kind, count: Option<usize>) {
        if let Some(active) = self.0.take() {
            active
                .recorder
                .event(active.id, active.phase, kind, Some(active.started), count);
        }
    }
}

impl Drop for Span {
    fn drop(&mut self) {
        self.end(Kind::Abandoned, None);
    }
}

pub(crate) fn measure_result<T, E>(
    phase: Phase,
    work: impl FnOnce() -> Result<T, E>,
) -> Result<T, E> {
    let trace = span(phase);
    let result = work();
    trace.finish(result.is_ok(), None);
    result
}

#[cfg(test)]
pub(crate) fn capture<T>(work: impl FnOnce() -> T) -> (T, Vec<serde_json::Value>) {
    let events = Rc::new(RefCell::new(Vec::new()));
    let mut recorder = Recorder::new();
    recorder.captured = Some(events.clone());
    let session = Session::activate(recorder);
    let result = work();
    drop(session);
    let captured = std::mem::take(&mut *events.borrow_mut());
    (result, captured)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn primer_import_trace_requires_exact_opt_in_without_changing_process_environment() {
        for setting in [
            None,
            Some(""),
            Some("0"),
            Some("true"),
            Some("01"),
            Some(" 1"),
        ] {
            let session = Session::from_setting(setting);
            assert!(!session.active);
            assert!(span(Phase::Import).0.is_none());
        }
        let session = Session::from_setting(Some("1"));
        assert!(session.active);
        drop(session);
        assert!(ACTIVE.with(|active| active.borrow().is_none()));
    }

    #[test]
    fn primer_import_trace_sink_failure_is_not_an_import_error() {
        struct FailedSink;
        impl Write for FailedSink {
            fn write(&mut self, _: &[u8]) -> std::io::Result<usize> {
                Err(std::io::ErrorKind::BrokenPipe.into())
            }
            fn flush(&mut self) -> std::io::Result<()> {
                Ok(())
            }
        }
        write_event(
            &mut FailedSink,
            &Event {
                schema: "gentle.primer_import_trace.v1",
                import: 1,
                span: 1,
                phase: Phase::Import,
                kind: Kind::Begin,
                elapsed_us: 0,
                duration_us: None,
                count: None,
            },
        );
    }

    #[test]
    fn primer_import_trace_is_inert_outside_a_session() {
        assert!(ACTIVE.with(|active| active.borrow().is_none()));
        assert!(span(Phase::AnnotationRead).0.is_none());
        assert_eq!(measure_result(Phase::Admission, || Ok::<_, ()>(7)), Ok(7));
        assert_eq!(
            measure_result(Phase::Admission, || Err::<(), _>("refusal")),
            Err("refusal")
        );
    }

    #[test]
    fn primer_import_trace_binds_spans_counts_errors_and_unfinished_work() {
        let (result, events) = capture(|| {
            let outer = span(Phase::Import);
            assert_eq!(
                measure_result(Phase::AnnotationRead, || Err::<(), _>("same error")),
                Err("same error")
            );
            drop(span(Phase::AnnotationDecode));
            outer.finish(true, Some(2));
            7
        });
        assert_eq!(result, 7);
        assert_eq!(events.len(), 6);
        assert_eq!(events[0]["span"], events[5]["span"]);
        assert_eq!(events[2]["kind"], "failed");
        assert_eq!(events[4]["kind"], "abandoned");
        assert_eq!(events[5]["count"], 2);
        assert!(
            events
                .iter()
                .all(|event| event["import"] == events[0]["import"])
        );
        assert!(
            events
                .iter()
                .all(|event| event.as_object().unwrap().keys().all(|key| {
                    [
                        "schema",
                        "import",
                        "span",
                        "phase",
                        "kind",
                        "elapsed_us",
                        "duration_us",
                        "count",
                    ]
                    .contains(&key.as_str())
                }))
        );
        assert!(
            events
                .windows(2)
                .all(|pair| pair[0]["elapsed_us"].as_u64() <= pair[1]["elapsed_us"].as_u64())
        );
    }

    #[test]
    fn primer_import_trace_is_bounded_and_restores_nested_and_panicked_contexts() {
        let (_, outer) = capture(|| {
            span(Phase::Import).finish(true, None);
            let (_, inner) = capture(|| {
                for _ in 0..MAX_EVENTS {
                    span(Phase::AnnotationRead).finish(true, None);
                }
            });
            assert_eq!(inner.len(), MAX_EVENTS);
            assert_eq!(inner.last().unwrap()["kind"], "truncated");
            assert!(
                std::panic::catch_unwind(|| capture(|| {
                    let _trace = span(Phase::AnnotationDecode);
                    panic!("diagnostic test");
                }))
                .is_err()
            );
            span(Phase::Persistence).finish(true, None);
        });
        assert_eq!(outer.len(), 4);
        assert!(ACTIVE.with(|active| active.borrow().is_none()));
    }
}
