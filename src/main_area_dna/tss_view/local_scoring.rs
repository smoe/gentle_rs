//! Independent, bounded and cancellable native TSS scoring jobs and result cache.

use super::*;
use crate::engine::{TfbsProgress, TfbsScoreTrackReport, TfbsScoreTrackValueKind};
use crate::tss_sequence_view::TssLocalScoreRequest;
use std::sync::atomic::{AtomicBool, Ordering};

#[derive(Debug)]
struct LocalScoreJob {
    source: Arc<TssSequenceView>,
    request: TssLocalScoreRequest,
    cancel: Arc<AtomicBool>,
    progress: Arc<Mutex<Option<TfbsProgress>>>,
    receiver: Mutex<Receiver<Result<TfbsScoreTrackReport, String>>>,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tss_local_scoring_job_acceptance_cancel_and_stale_source_are_independent() {
        let _guard = crate::tf_motifs::test_registry_lock().lock().unwrap();
        let dna = crate::tss_sequence_view::tests::fixture(false);
        let source = Arc::new(TssSequenceView::from_dna(&dna).unwrap());
        let mut area = MainAreaDna::new(dna.clone(), None, None);
        area.tss_view_available();
        area.tss_ui.document = Some(Ok(source.clone()));
        area.tss_ui.local_scores.matrix_text = "MA0004.1".into();
        let request = area.tss_ui.local_scores.request();
        let report = source
            .compute_local_scores(&dna.get_forward_string(), &request, &mut |_| true)
            .unwrap();
        let ctx = egui::Context::default();

        for stale in [false, true] {
            area.tss_ui.document = Some(Ok(source.clone()));
            let (send, receiver) = std::sync::mpsc::channel();
            let cancel = Arc::new(AtomicBool::new(false));
            area.tss_ui.local_scores.job = Some(Arc::new(LocalScoreJob {
                source: source.clone(),
                request: request.clone(),
                cancel: cancel.clone(),
                progress: Arc::new(Mutex::new(None)),
                receiver: Mutex::new(receiver),
            }));
            send.send(Ok(report.clone())).unwrap();
            if stale {
                area.tss_ui.document = Some(Ok(Arc::new((*source).clone())));
            }
            area.poll_tss_local_scores(&ctx);
            assert!(area.tss_ui.local_scores.job.is_none());
            assert!(cancel.load(Ordering::Relaxed));
            let displayed = area.tss_ui.document.as_ref().unwrap().as_ref().unwrap();
            assert_eq!(displayed.local_scoring.is_some(), !stale);
            assert!(area.tfbs_task.is_none());
            if !stale {
                assert!(area.tss_ui.local_scores.report.is_some());
                area.tss_ui.traces = false;
                let lanes = area.tss_ui.visible_lanes(displayed);
                assert!(
                    lanes
                        .iter()
                        .any(|&i| displayed.lanes[i].kind == TssLaneKind::LocalScoreTrace)
                );
                area.tss_ui.local_traces = false;
                assert!(
                    area.tss_ui
                        .visible_lanes(displayed)
                        .iter()
                        .all(|&i| displayed.lanes[i].kind != TssLaneKind::LocalScoreTrace)
                );
            }
        }
        // Form changes cannot relabel old results, even for the same source object.
        area.tss_ui.document = Some(Ok(source.clone()));
        let (send, receiver) = std::sync::mpsc::channel();
        area.tss_ui.local_scores.job = Some(Arc::new(LocalScoreJob {
            source: source.clone(),
            request: request.clone(),
            cancel: Arc::new(AtomicBool::new(false)),
            progress: Arc::new(Mutex::new(None)),
            receiver: Mutex::new(receiver),
        }));
        send.send(Ok(report)).unwrap();
        area.tss_ui.local_scores.clip_negative = !request.clip_negative;
        area.poll_tss_local_scores(&ctx);
        assert!(
            area.tss_ui
                .document
                .as_ref()
                .unwrap()
                .as_ref()
                .unwrap()
                .local_scoring
                .is_none()
        );
        assert!(
            area.tss_ui
                .local_scores
                .error
                .as_ref()
                .unwrap()
                .contains("request changed")
        );

        let (_, receiver) = std::sync::mpsc::channel();
        let cancel = Arc::new(AtomicBool::new(false));
        area.tss_ui.local_scores.job = Some(Arc::new(LocalScoreJob {
            source,
            request,
            cancel: cancel.clone(),
            progress: Arc::new(Mutex::new(None)),
            receiver: Mutex::new(receiver),
        }));
        assert!(area.tss_svg_snapshot(ViewSvgExportProfile::Screen).is_err());
        area.replace_loaded_sequence(crate::tss_sequence_view::tests::fixture(true));
        area.tss_view_available();
        assert!(
            cancel.load(Ordering::Relaxed),
            "replaced DNA cancels the old worker"
        );
        assert!(area.tss_ui.local_scores.report.is_none());
    }

    #[test]
    fn tss_local_scoring_background_worker_never_uses_tfbs_panel_slot() {
        let _guard = crate::tf_motifs::test_registry_lock().lock().unwrap();
        let dna = crate::tss_sequence_view::tests::fixture(false);
        let source = Arc::new(TssSequenceView::from_dna(&dna).unwrap());
        let mut area = MainAreaDna::new(dna, None, None);
        area.tss_view_available();
        area.tss_ui.document = Some(Ok(source.clone()));
        area.tss_ui.local_scores.matrix_text = "MA0004.1".into();
        let ctx = egui::Context::default();
        area.start_tss_local_scores(source, &ctx).unwrap();
        assert!(area.tfbs_task.is_none());
        let deadline = std::time::Instant::now() + std::time::Duration::from_secs(30);
        while area.tss_ui.local_scores.running() && std::time::Instant::now() < deadline {
            area.poll_tss_local_scores(&ctx);
            std::thread::sleep(std::time::Duration::from_millis(5));
        }
        assert!(!area.tss_ui.local_scores.running());
        assert!(
            area.tss_ui.local_scores.error.is_none(),
            "{:?}",
            area.tss_ui.local_scores.error
        );
        assert!(area.tss_ui.local_scores.report.is_some());
        assert!(area.tfbs_task.is_none());
    }
}

impl Drop for LocalScoreJob {
    fn drop(&mut self) {
        self.cancel.store(true, Ordering::Relaxed);
    }
}

#[derive(Clone, Debug, Default)]
pub(super) struct LocalScoreState {
    matrix_text: String,
    score_kind: TfbsScoreTrackValueKind,
    clip_negative: bool,
    job: Option<Arc<LocalScoreJob>>,
    report: Option<Arc<TfbsScoreTrackReport>>,
    error: Option<String>,
}

impl LocalScoreState {
    pub(super) fn running(&self) -> bool {
        self.job.is_some()
    }

    fn request(&self) -> TssLocalScoreRequest {
        TssLocalScoreRequest {
            matrix_ids: self
                .matrix_text
                .split(|c: char| c == ',' || c.is_whitespace())
                .filter(|s| !s.is_empty())
                .map(str::to_owned)
                .collect(),
            score_kind: self.score_kind,
            clip_negative: self.clip_negative,
        }
    }
}

impl MainAreaDna {
    fn start_tss_local_scores(
        &mut self,
        source: Arc<TssSequenceView>,
        ctx: &egui::Context,
    ) -> Result<(), String> {
        if self.tss_ui.local_scores.running() || self.tss_ui.profile_load.is_some() {
            return Err("A TSS scoring/attachment job is already running".into());
        }
        if !self
            .tss_ui
            .document
            .as_ref()
            .and_then(|d| d.as_ref().ok())
            .is_some_and(|current| Arc::ptr_eq(current, &source))
        {
            return Err("TSS document changed before local scoring admission".into());
        }
        let request = self.tss_ui.local_scores.request();
        request.validate_budget(source.geometry.length().ok_or("Invalid TSS geometry")?)?;
        let sequence = self
            .dna
            .read()
            .map_err(|_| "Could not read TSS DNA")?
            .get_forward_string();
        source.local_score_target(&sequence)?;
        let cancel = Arc::new(AtomicBool::new(false));
        let progress = Arc::new(Mutex::new(None));
        let (send, receiver) = std::sync::mpsc::channel();
        let worker_source = source.clone();
        let worker_request = request.clone();
        let worker_cancel = cancel.clone();
        let worker_progress = progress.clone();
        let ctx = ctx.clone();
        // Matrix resolution, calibration and scoring stay off the UI thread.
        std::thread::Builder::new()
            .name("tss-local-scores".into())
            .stack_size(8 * 1024 * 1024)
            .spawn(move || {
                let result =
                    worker_source.compute_local_scores(&sequence, &worker_request, &mut |event| {
                        if let OperationProgress::Tfbs(p) = event {
                            if let Ok(mut latest) = worker_progress.lock() {
                                *latest = Some(p);
                            }
                        }
                        ctx.request_repaint();
                        !worker_cancel.load(Ordering::Relaxed)
                    });
                if !worker_cancel.load(Ordering::Relaxed) {
                    let _ = send.send(result);
                }
                ctx.request_repaint();
            })
            .map_err(|e| format!("Cannot start TSS scoring worker: {e}"))?;
        self.tss_ui.local_scores.job = Some(Arc::new(LocalScoreJob {
            source,
            request,
            cancel,
            progress,
            receiver: Mutex::new(receiver),
        }));
        self.tss_ui.local_scores.error = None;
        Ok(())
    }

    pub(super) fn poll_tss_local_scores(&mut self, ctx: &egui::Context) {
        let Some(job) = self.tss_ui.local_scores.job.clone() else {
            return;
        };
        let received = job
            .receiver
            .lock()
            .map_err(|_| TryRecvError::Disconnected)
            .and_then(|r| r.try_recv());
        if matches!(received, Err(TryRecvError::Empty)) {
            ctx.request_repaint_after(std::time::Duration::from_millis(100));
            return;
        }
        self.tss_ui.local_scores.job = None;
        let current = self.tss_ui.document.as_ref().and_then(|d| d.as_ref().ok());
        if job.cancel.load(Ordering::Relaxed)
            || !current.is_some_and(|view| Arc::ptr_eq(view, &job.source))
            || self.tss_ui.local_scores.request() != job.request
        {
            self.tss_ui.local_scores.error =
                Some("Local scoring result discarded: document or request changed".into());
            return;
        }
        let result = received
            .map_err(|_| "Local scoring worker stopped".to_owned())
            .and_then(|r| r)
            .and_then(|report| {
                job.source
                    .with_local_scores(&job.request, &report)
                    .map(|view| (report, view))
            });
        match result {
            Ok((report, view)) => {
                self.tss_ui.local_scores.report = Some(Arc::new(report));
                self.tss_ui.document = Some(Ok(Arc::new(view)));
                self.tss_ui.selected = None;
            }
            Err(error) => self.tss_ui.local_scores.error = Some(error),
        }
    }

    pub(super) fn render_tss_local_scoring(
        &mut self,
        ui: &mut egui::Ui,
        view: &Arc<TssSequenceView>,
    ) {
        ui.collapsing("Local scoring", |ui| {
            ui.horizontal_wrapped(|ui| {
                ui.label("Matrix IDs");
                ui.add(egui::TextEdit::singleline(&mut self.tss_ui.local_scores.matrix_text)
                    .desired_width(260.0).hint_text("MA0861.2, MA0106.3"));
                egui::ComboBox::from_id_salt(("tss_local_score_kind", self.panel_scope_key()))
                    .selected_text(self.tss_ui.local_scores.score_kind.as_str())
                    .show_ui(ui, |ui| {
                        use TfbsScoreTrackValueKind::*;
                        for kind in [LlrBits, LlrQuantile, LlrBackgroundQuantile, LlrBackgroundTailLog10,
                            TrueLogOddsBits, TrueLogOddsQuantile, TrueLogOddsBackgroundQuantile, TrueLogOddsBackgroundTailLog10] {
                            ui.selectable_value(&mut self.tss_ui.local_scores.score_kind, kind, kind.as_str());
                        }
                    });
                ui.checkbox(&mut self.tss_ui.local_scores.clip_negative, "Clip negative scores");
            });
            ui.horizontal_wrapped(|ui| {
                if ui.add_enabled(!self.tss_ui.local_scores.running() && self.tss_ui.profile_load.is_none(),
                    egui::Button::new("Compute local scores")).clicked() {
                    if let Err(error) = self.start_tss_local_scores(view.clone(), ui.ctx()) {
                        self.tss_ui.local_scores.error = Some(error);
                    }
                }
                if let Some(job) = &self.tss_ui.local_scores.job {
                    ui.spinner();
                    if let Some(p) = job.progress.lock().ok().and_then(|v| v.clone()) {
                        ui.add(egui::ProgressBar::new((p.total_percent / 100.0) as f32)
                            .desired_width(150.0).text(format!("{}/{} {}", p.motif_index, p.motif_count,
                                p.stage_label.as_deref().unwrap_or("scoring"))));
                    }
                    if ui.button("Cancel scoring").clicked() {
                        // Dropping the job requests cooperative cancellation, retaining the last complete result.
                        self.tss_ui.local_scores.job = None;
                    }
                } else if view.local_scoring.is_some() && self.tss_ui.profile_load.is_none()
                    && ui.button("Clear local scores").clicked() {
                    let mut cleared = (**view).clone();
                    cleared.clear_local_scores();
                    self.tss_ui.document = Some(Ok(Arc::new(cleared)));
                    self.tss_ui.local_scores.report = None;
                    self.tss_ui.selected = None;
                }
            });
            if let Some(local) = &view.local_scoring {
                if local.request != self.tss_ui.local_scores.request() {
                    ui.label("Displayed local curves belong to the previous request; pending settings have not been computed.");
                }
            }
            if let Some(error) = &self.tss_ui.local_scores.error {
                ui.colored_label(ui.visuals().error_fg_color, error);
            }
        });
    }
}
