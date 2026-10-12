//! Refreshable, project-independent template-import menu snapshots.

use super::*;
use crate::runtime_assets::resolve_runtime_asset_path;

struct PatternCatalogSnapshot {
    root: PathBuf,
    catalog_path: PathBuf,
    entries: Result<Vec<CloningPatternCatalogEntry>, String>,
    routines: Result<Vec<CloningRoutineCatalogRow>, String>,
    grna: HashMap<grna_routine_ui::GrnaRoutine, Result<CloningRoutineCatalogRow, String>>,
}

impl PatternCatalogSnapshot {
    fn load(root: PathBuf, catalog_path: PathBuf) -> Self {
        let root = root.canonicalize().unwrap_or(root);
        let catalog_path = catalog_path.canonicalize().unwrap_or(catalog_path);
        let mut routines =
            crate::engine_shell::cloning_routine_catalog_rows(&catalog_path.to_string_lossy())
                .and_then(|rows| {
                    serde_json::from_value::<Vec<CloningRoutineCatalogRow>>(rows)
                        .map_err(|err| err.to_string())
                });
        let mut titles = HashMap::new();
        if let Ok(rows) = &mut routines {
            for row in rows.iter_mut() {
                if let Some(path) = &mut row.template_path {
                    let relative = Path::new(path);
                    let resolved = if relative.is_absolute() {
                        relative.to_path_buf()
                    } else {
                        catalog_path
                            .parent()
                            .unwrap_or(Path::new("."))
                            .join(relative.strip_prefix("assets").unwrap_or(relative))
                    };
                    let resolved = resolved.canonicalize().unwrap_or(resolved);
                    titles.insert(resolved.clone(), row.title.clone());
                    *path = resolved.to_string_lossy().into_owned();
                }
            }
            rows.sort_by(|a, b| a.family.cmp(&b.family).then(a.title.cmp(&b.title)));
        }
        let mut entries = GENtleApp::collect_cloning_pattern_catalog_entries(&root);
        if let Ok(entries) = &mut entries {
            Self::apply_titles(entries, &titles);
        }
        let grna = grna_routine_ui::GrnaRoutine::ALL
            .into_iter()
            .map(|routine| {
                let result = routines.as_ref().map_err(Clone::clone).and_then(|rows| {
                    let row = rows
                        .iter()
                        .find(|row| row.routine_id == routine.id())
                        .ok_or_else(|| format!("Routine '{}' is absent", routine.id()))?;
                    routine.validate(row)?;
                    Ok(row.clone())
                });
                (routine, result)
            })
            .collect();
        Self {
            root,
            catalog_path,
            entries,
            routines,
            grna,
        }
    }

    fn apply_titles(entries: &mut [CloningPatternCatalogEntry], titles: &HashMap<PathBuf, String>) {
        for entry in entries {
            if entry.is_file {
                let path = PathBuf::from(&entry.path);
                let path = path.canonicalize().unwrap_or(path);
                if let Some(title) = titles.get(&path) {
                    entry.label = title.clone();
                }
                entry.path = path.to_string_lossy().into_owned();
            } else {
                Self::apply_titles(&mut entry.children, titles);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::atomic::AtomicUsize;

    fn wait_for_load(cache: &mut PatternCatalogCache) {
        let deadline = Instant::now() + Duration::from_secs(30);
        while cache.receiver.is_some() {
            cache.poll();
            assert!(Instant::now() < deadline, "catalog worker timed out");
            std::thread::sleep(Duration::from_millis(1));
        }
        assert!(cache.error.is_none(), "{:?}", cache.error);
    }

    fn leaves(entries: &[CloningPatternCatalogEntry]) -> Vec<&CloningPatternCatalogEntry> {
        entries
            .iter()
            .flat_map(|entry| {
                if entry.is_file {
                    vec![entry]
                } else {
                    leaves(&entry.children)
                }
            })
            .collect()
    }

    #[test]
    fn canonical_catalog_titles_and_template_paths_survive_browsing() {
        let base = Path::new(env!("CARGO_MANIFEST_DIR"));
        let snapshot = PatternCatalogSnapshot::load(
            base.join(DEFAULT_CLONING_PATTERN_CATALOG_DIR),
            base.join(DEFAULT_CLONING_ROUTINE_CATALOG_PATH),
        );
        let entries = snapshot.entries.as_ref().unwrap();
        let files = leaves(entries);
        assert!(!files.is_empty());
        let raw = GENtleApp::collect_cloning_pattern_catalog_entries(&snapshot.root).unwrap();
        assert_eq!(files.len(), leaves(&raw).len());
        assert!(entries.iter().any(|entry| entry.label == "CRISPR"));
        assert!(
            files
                .iter()
                .any(|entry| entry.label == "gRNA Anchor Window Scan")
        );
        assert!(files.iter().all(|entry| !entry.label.contains("Grna")));
        for routine in snapshot.routines.as_ref().unwrap() {
            if let Some(path) = &routine.template_path {
                let entry = files
                    .iter()
                    .find(|entry| &entry.path == path)
                    .expect("catalog template in tree");
                assert_eq!(entry.label, routine.title);
            }
        }
        assert_eq!(
            GENtleApp::humanize_catalog_label("grna_candidate-priority_scan"),
            "gRNA Candidate Priority Scan"
        );
        assert_eq!(
            GENtleApp::humanize_catalog_label("crispr_pcr_dna_rna_tfbs"),
            "CRISPR PCR DNA RNA TFBS"
        );
    }

    #[test]
    fn steady_menu_frames_use_snapshot_and_explicit_refresh_reloads() {
        // Synthetic temporary import catalog; never executed or committed as a fixture.
        let dir = tempfile::tempdir().unwrap();
        let root = dir.path().join("patterns");
        fs::create_dir(&root).unwrap();
        fs::write(root.join("grna_toy.json"), "{}").unwrap();
        let catalog = dir.path().join("routines.json");
        fs::write(
            &catalog,
            r#"{"schema":"gentle.cloning_routines.v1","routines":[]}"#,
        )
        .unwrap();
        let ctx = egui::Context::default();
        let mut app = GENtleApp::default();
        let loads = Arc::new(AtomicUsize::new(0));
        let start = |app: &mut GENtleApp| {
            let (root, catalog, loads) = (root.clone(), catalog.clone(), loads.clone());
            app.pattern_catalog.start_with(&ctx, move || {
                loads.fetch_add(1, Ordering::SeqCst);
                PatternCatalogSnapshot::load(root, catalog)
            });
        };
        start(&mut app);
        wait_for_load(&mut app.pattern_catalog);
        let first = app.pattern_catalog.snapshot.clone().unwrap();
        assert_eq!(first.entries.as_ref().unwrap().len(), 1);
        fs::remove_dir_all(&root).unwrap();
        fs::remove_file(&catalog).unwrap();
        for _ in 0..20 {
            ctx.run_ui(egui::RawInput::default(), |ui| {
                app.render_pattern_template_catalog_menu(ui);
                app.render_pattern_routine_catalog_menu(ui);
            })
            .drop_without_applying_deltas();
        }
        assert_eq!(loads.load(Ordering::SeqCst), 1);
        assert!(Arc::ptr_eq(
            &first,
            app.pattern_catalog.snapshot.as_ref().unwrap()
        ));
        assert!(app.engine.read().unwrap().state().sequences.is_empty());
        start(&mut app);
        wait_for_load(&mut app.pattern_catalog);
        assert_eq!(loads.load(Ordering::SeqCst), 2);
        let refreshed = app.pattern_catalog.snapshot.as_ref().unwrap();
        assert!(refreshed.entries.is_err());
        assert!(refreshed.routines.is_err());
    }

    #[test]
    fn pending_worker_does_not_block_menu_or_spawn_again() {
        let ctx = egui::Context::default();
        let mut app = GENtleApp::default();
        let (release, wait) = mpsc::channel();
        app.pattern_catalog.start_with(&ctx, move || {
            wait.recv_timeout(Duration::from_secs(30)).unwrap();
            PatternCatalogSnapshot {
                root: "missing".into(),
                catalog_path: "missing".into(),
                entries: Err("missing".into()),
                routines: Err("missing".into()),
                grna: HashMap::new(),
            }
        });
        for _ in 0..5 {
            ctx.run_ui(egui::RawInput::default(), |ui| {
                app.render_pattern_template_catalog_menu(ui)
            })
            .drop_without_applying_deltas();
        }
        app.pattern_catalog
            .start_with(&ctx, || panic!("must not replace running loader"));
        assert!(app.pattern_catalog.snapshot.is_none());
        release.send(()).unwrap();
        wait_for_load(&mut app.pattern_catalog);
    }

    #[test]
    fn disconnected_worker_is_retryable_without_automatic_loop() {
        let (sender, receiver) = mpsc::channel();
        drop(sender);
        let mut cache = PatternCatalogCache {
            started: true,
            receiver: Some(receiver),
            ..Default::default()
        };
        cache.poll();
        assert!(cache.receiver.is_none());
        assert!(cache.error.as_ref().unwrap().contains("retry"));
        assert!(cache.started);
    }
}

#[cfg(test)]
pub(super) fn load_test_catalog() -> PatternCatalogCache {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    PatternCatalogCache {
        started: true,
        snapshot: Some(Arc::new(PatternCatalogSnapshot::load(
            root.join(DEFAULT_CLONING_PATTERN_CATALOG_DIR),
            root.join(DEFAULT_CLONING_ROUTINE_CATALOG_PATH),
        ))),
        ..Default::default()
    }
}

#[derive(Default)]
pub(super) struct PatternCatalogCache {
    started: bool,
    receiver: Option<mpsc::Receiver<PatternCatalogSnapshot>>,
    snapshot: Option<Arc<PatternCatalogSnapshot>>,
    error: Option<String>,
}

impl PatternCatalogCache {
    pub(super) fn ensure_started(&mut self, ctx: &egui::Context) {
        if !self.started {
            self.start(ctx);
        }
        self.poll();
    }

    pub(super) fn grna_routine(
        &self,
        routine: grna_routine_ui::GrnaRoutine,
    ) -> Result<CloningRoutineCatalogRow, ActionReadiness> {
        if let Some(error) = &self.error {
            return Err(ActionReadiness::AdapterUnavailable {
                detail: error.clone(),
            });
        }
        let Some(snapshot) = &self.snapshot else {
            return Err(ActionReadiness::Checking {
                detail: "Loading routine catalog; use Refresh template catalog to retry".into(),
            });
        };
        snapshot
            .grna
            .get(&routine)
            .cloned()
            .unwrap_or_else(|| Err("Routine not in this catalog".into()))
            .map_err(|error| ActionReadiness::AdapterUnavailable {
                detail: format!(
                    "{error}. Refresh template catalog or import a template explicitly."
                ),
            })
    }

    fn start(&mut self, ctx: &egui::Context) {
        self.start_with(ctx, || {
            PatternCatalogSnapshot::load(
                resolve_runtime_asset_path(DEFAULT_CLONING_PATTERN_CATALOG_DIR),
                resolve_runtime_asset_path(DEFAULT_CLONING_ROUTINE_CATALOG_PATH),
            )
        });
    }

    fn start_with(
        &mut self,
        ctx: &egui::Context,
        load: impl FnOnce() -> PatternCatalogSnapshot + Send + 'static,
    ) {
        if self.receiver.is_some() {
            return;
        }
        self.started = true;
        self.snapshot = None;
        self.error = None;
        let (sender, receiver) = mpsc::channel();
        let ctx = ctx.clone();
        match std::thread::Builder::new()
            .name("pattern-catalog".into())
            .spawn(move || {
                let snapshot = load();
                let _ = sender.send(snapshot);
                ctx.request_repaint();
            }) {
            Ok(_) => self.receiver = Some(receiver),
            Err(err) => self.error = Some(format!("Could not start catalog loading: {err}")),
        }
    }

    fn poll(&mut self) {
        let Some(receiver) = &self.receiver else {
            return;
        };
        match receiver.try_recv() {
            Ok(snapshot) => {
                self.snapshot = Some(Arc::new(snapshot));
                self.receiver = None;
            }
            Err(mpsc::TryRecvError::Empty) => {}
            Err(mpsc::TryRecvError::Disconnected) => {
                self.receiver = None;
                self.error =
                    Some("Catalog loading stopped; use Refresh template catalog to retry".into());
            }
        }
    }
}

impl GENtleApp {
    pub(super) fn render_pattern_template_catalog_menu(&mut self, ui: &mut Ui) {
        if !self.pattern_catalog.started {
            self.pattern_catalog.start(ui.ctx());
        }
        self.pattern_catalog.poll();
        if ui
            .add_enabled(
                self.pattern_catalog.receiver.is_none(),
                egui::Button::new("Refresh template catalog"),
            )
            .clicked()
        {
            self.pattern_catalog.start(ui.ctx());
        }
        let Some(snapshot) = self.pattern_catalog.snapshot.clone() else {
            if let Some(err) = &self.pattern_catalog.error {
                ui.label(err);
            } else {
                ui.label("Loading template catalog...");
                ui.ctx().request_repaint_after(Duration::from_millis(50));
            }
            return;
        };
        match &snapshot.entries {
            Ok(entries) if entries.is_empty() => {
                ui.label(self.tr("menu.patterns.catalog_empty"));
            }
            Ok(entries) => {
                if ui
                    .button(self.tr("menu.patterns.import_full_catalog"))
                    .on_hover_text(format!(
                        "Import all templates from {}; does not run them",
                        snapshot.root.display()
                    ))
                    .clicked()
                {
                    self.import_workflow_macro_templates_from_path(
                        &snapshot.root.to_string_lossy(),
                    );
                    ui.close();
                }
                ui.small("Import templates (does not run a workflow)");
                let mut selected_path = None;
                Self::render_cloning_pattern_catalog_menu_entries(ui, entries, &mut selected_path);
                if let Some(path) = selected_path {
                    self.import_workflow_macro_templates_from_path(&path);
                    ui.close();
                }
            }
            Err(err) => {
                ui.label(self.tr("menu.patterns.catalog_unavailable"));
                ui.small(err);
            }
        }
    }

    pub(super) fn render_pattern_routine_catalog_menu(&mut self, ui: &mut Ui) {
        let Some(snapshot) = self.pattern_catalog.snapshot.clone() else {
            return;
        };
        ui.small("Import routine templates");
        let routines = match &snapshot.routines {
            Ok(routines) => routines,
            Err(err) => {
                ui.label(self.tr("menu.patterns.routine_catalog_unavailable"));
                ui.small(err);
                return;
            }
        };
        if ui
            .button(self.tr("menu.patterns.show_routine_catalog_summary"))
            .clicked()
        {
            self.app_status = format!(
                "Loaded {} routine(s) from '{}' (template import only)",
                routines.len(),
                snapshot.catalog_path.display()
            );
            ui.close();
        }
        let mut by_family: BTreeMap<&str, Vec<CloningRoutineCatalogRow>> = BTreeMap::new();
        let mut by_status: BTreeMap<&str, Vec<CloningRoutineCatalogRow>> = BTreeMap::new();
        for routine in routines {
            by_family
                .entry(&routine.family)
                .or_default()
                .push(routine.clone());
            by_status
                .entry(&routine.status)
                .or_default()
                .push(routine.clone());
        }
        let mut selected_path = None;
        let mut status_message = None;
        for (label, groups) in [
            ("menu.patterns.browse_by_family", &by_family),
            ("menu.patterns.browse_by_status", &by_status),
        ] {
            ui.menu_button(self.tr(label), |ui| {
                for (group, rows) in groups {
                    ui.menu_button(format!("{group} ({})", rows.len()), |ui| {
                        Self::render_cloning_routine_menu_entries(
                            ui,
                            rows,
                            &mut selected_path,
                            &mut status_message,
                        );
                    });
                }
            });
        }
        if let Some(path) = selected_path {
            self.import_workflow_macro_templates_from_path(&path);
            ui.close();
        } else if let Some(message) = status_message {
            self.app_status = message;
        }
    }
}
