use std::{
    collections::HashMap,
    hash::{Hash, Hasher},
    path::Path,
    sync::{Arc, RwLock},
};

use crate::{
    about,
    dna_sequence::{self, DNAsequence},
    engine::{DisplayTarget, Engine, GentleEngine, Operation, ProjectState},
    enzymes,
    icons::APP_ICON,
    lineage_export::export_lineage_svg,
    resource_sync, tf_motifs,
    window::Window,
    TRANSLATIONS,
};
use anyhow::{anyhow, Result};
use eframe::egui::{self, menu, Key, KeyboardShortcut, Modifiers, Pos2, Ui, Vec2, ViewportId};

pub struct GENtleApp {
    engine: Arc<RwLock<GentleEngine>>,
    new_windows: Vec<Window>,
    windows: HashMap<ViewportId, Arc<RwLock<Window>>>,
    windows_to_close: Arc<RwLock<Vec<ViewportId>>>,
    viewport_id_counter: usize,
    update_has_run_before: bool,
    show_about_dialog: bool,
    current_project_path: Option<String>,
    lineage_graph_view: bool,
    clean_state_fingerprint: u64,
    pending_project_action: Option<ProjectAction>,
}

#[derive(Clone, Copy)]
enum ProjectAction {
    New,
    Open,
    Close,
}

impl Default for GENtleApp {
    fn default() -> Self {
        Self {
            engine: Arc::new(RwLock::new(GentleEngine::new())),
            new_windows: vec![],
            windows: HashMap::new(),
            windows_to_close: Arc::new(RwLock::new(vec![])),
            viewport_id_counter: 0,
            update_has_run_before: false,
            show_about_dialog: false,
            current_project_path: None,
            lineage_graph_view: false,
            clean_state_fingerprint: 0,
            pending_project_action: None,
        }
    }
}

impl GENtleApp {
    pub fn new() -> Self {
        let mut app = Self::default();
        app.mark_clean_snapshot();
        app
    }

    pub fn new_with_project(project_path: Option<&str>) -> Self {
        let mut ret = Self::new();
        if let Some(path) = project_path {
            let _ = ret.load_project_from_file(path);
        }
        ret
    }

    fn load_dna_from_genbank_file(filename: &str) -> Result<DNAsequence> {
        let dna = dna_sequence::DNAsequence::from_genbank_file(filename)?
            .pop()
            .ok_or_else(|| anyhow!("Could not read GenBank file {filename}"))?;
        Ok(dna)
    }

    fn load_dna_from_fasta_file(filename: &str) -> Result<DNAsequence> {
        let dna = dna_sequence::DNAsequence::from_fasta_file(filename)?
            .pop()
            .ok_or_else(|| anyhow!("Could not read fasta file {filename}"))?;
        Ok(dna)
    }

    fn new_dna_window(&mut self, seq_id: String, dna: DNAsequence) {
        self.new_windows
            .push(Window::new_dna(dna, seq_id, self.engine.clone()));
    }

    fn open_sequence_window(&mut self, seq_id: &str) {
        let dna = self
            .engine
            .read()
            .unwrap()
            .state()
            .sequences
            .get(seq_id)
            .cloned();
        if let Some(dna) = dna {
            self.new_dna_window(seq_id.to_string(), dna);
        }
    }

    fn open_pool_window(&mut self, representative_seq_id: &str, pool_seq_ids: Vec<String>) {
        let dna = self
            .engine
            .read()
            .unwrap()
            .state()
            .sequences
            .get(representative_seq_id)
            .cloned();
        if let Some(dna) = dna {
            let mut window =
                Window::new_dna(dna, representative_seq_id.to_string(), self.engine.clone());
            window.set_pool_context(pool_seq_ids);
            self.new_windows.push(window);
        }
    }

    pub fn load_from_file(path: &str) -> Result<DNAsequence> {
        if let Ok(dna) = Self::load_dna_from_genbank_file(path) {
            Ok(dna)
        } else if let Ok(dna) = Self::load_dna_from_fasta_file(path) {
            Ok(dna)
        } else {
            Err(anyhow!("Could not load file '{path}'"))
        }
    }

    fn open_new_window_from_file(&mut self, path: &str) {
        let op = Operation::LoadFile {
            path: path.to_string(),
            as_id: None,
        };
        let load_result = {
            let mut engine = self.engine.write().unwrap();
            engine.apply(op)
        };
        if let Ok(result) = load_result {
            if let Some(seq_id) = result.created_seq_ids.first() {
                let seq_id = seq_id.to_string();
                let dna = self
                    .engine
                    .read()
                    .unwrap()
                    .state()
                    .sequences
                    .get(&seq_id)
                    .cloned();
                if let Some(dna) = dna {
                    self.new_dna_window(seq_id, dna);
                    return;
                }
            }
            // TODO warning: load op ran but no created sequence
        } else {
            // TODO error, could not load file through engine
        }
    }

    fn save_project_to_file(&self, path: &str) -> Result<()> {
        self.engine
            .read()
            .unwrap()
            .state()
            .save_to_path(path)
            .map_err(|e| anyhow!(e.to_string()))
    }

    fn current_state_fingerprint(&self) -> u64 {
        let state_json = {
            let engine = self.engine.read().unwrap();
            serde_json::to_vec(engine.state()).unwrap_or_default()
        };
        let mut hasher = std::collections::hash_map::DefaultHasher::new();
        state_json.hash(&mut hasher);
        hasher.finish()
    }

    fn mark_clean_snapshot(&mut self) {
        self.clean_state_fingerprint = self.current_state_fingerprint();
    }

    fn is_project_dirty(&self) -> bool {
        self.current_state_fingerprint() != self.clean_state_fingerprint
    }

    fn reset_to_empty_project(&mut self) {
        self.engine = Arc::new(RwLock::new(GentleEngine::new()));
        self.current_project_path = None;
        self.new_windows.clear();
        self.windows.clear();
        self.windows_to_close.write().unwrap().clear();
        self.viewport_id_counter = 0;
        self.pending_project_action = None;
        self.mark_clean_snapshot();
    }

    fn prompt_open_project(&mut self) {
        if let Some(path) = rfd::FileDialog::new()
            .add_filter("GENtle project", &["json"])
            .pick_file()
        {
            let path = path.display().to_string();
            let _ = self.load_project_from_file(&path);
        }
    }

    fn prompt_open_sequence(&mut self) {
        if let Some(path) = rfd::FileDialog::new().pick_file() {
            let path = path.display().to_string();
            self.open_new_window_from_file(&path);
        }
    }

    fn prompt_save_project(&mut self) {
        if let Some(path) = rfd::FileDialog::new()
            .set_file_name("project.gentle.json")
            .add_filter("GENtle project", &["json"])
            .save_file()
        {
            let path = path.display().to_string();
            if self.save_project_to_file(&path).is_ok() {
                self.current_project_path = Some(path);
                self.mark_clean_snapshot();
            }
        }
    }

    fn prompt_export_lineage_svg(&mut self) {
        if let Some(path) = rfd::FileDialog::new()
            .set_file_name("lineage.svg")
            .add_filter("SVG", &["svg"])
            .save_file()
        {
            let path = path.display().to_string();
            let svg = {
                let engine = self.engine.read().unwrap();
                export_lineage_svg(engine.state())
            };
            let _ = std::fs::write(path, svg);
        }
    }

    fn refresh_project_restriction_enzymes(&mut self, resource_path: &str) -> Result<usize> {
        let enzymes = enzymes::load_restriction_enzymes_from_path(resource_path)?;
        let mut engine = self.engine.write().unwrap();
        for dna in engine.state_mut().sequences.values_mut() {
            enzymes.clone_into(dna.restriction_enzymes_mut());
            dna.update_computed_features();
        }
        Ok(enzymes.len())
    }

    fn prompt_import_rebase_resource(&mut self) {
        if let Some(path) = rfd::FileDialog::new()
            .add_filter("REBASE withrefm", &["txt", "withrefm", "f19", "f20", "f21"])
            .pick_file()
        {
            let input = path.display().to_string();
            match resource_sync::sync_rebase(
                &input,
                Some(resource_sync::DEFAULT_REBASE_RESOURCE_PATH),
                false,
            ) {
                Ok(report) => {
                    let seq_count = self.engine.read().unwrap().state().sequences.len();
                    let refresh_status = self.refresh_project_restriction_enzymes(&report.output);
                    match refresh_status {
                        Ok(loaded_count) => println!(
                            "Imported REBASE from '{}' ({} enzymes); active set now {} enzymes, refreshed {} sequence(s)",
                            report.source, report.item_count, loaded_count, seq_count
                        ),
                        Err(e) => eprintln!(
                            "Imported REBASE from '{}', but could not refresh loaded sequences: {}",
                            report.source, e
                        ),
                    }
                }
                Err(e) => eprintln!("Could not import REBASE from '{}': {}", input, e),
            }
        }
    }

    fn prompt_import_jaspar_resource(&mut self) {
        if let Some(path) = rfd::FileDialog::new()
            .add_filter("JASPAR motif matrix", &["txt", "jaspar", "pfm"])
            .pick_file()
        {
            let input = path.display().to_string();
            match resource_sync::sync_jaspar(
                &input,
                Some(resource_sync::DEFAULT_JASPAR_RESOURCE_PATH),
            ) {
                Ok(report) => {
                    tf_motifs::reload();
                    println!(
                        "Imported JASPAR from '{}' ({} motifs) to '{}'",
                        report.source, report.item_count, report.output
                    );
                }
                Err(e) => eprintln!("Could not import JASPAR from '{}': {}", input, e),
            }
        }
    }

    fn save_current_project(&mut self) -> bool {
        if let Some(path) = self.current_project_path.clone() {
            if self.save_project_to_file(&path).is_ok() {
                self.mark_clean_snapshot();
                return true;
            }
            return false;
        }

        if let Some(path) = rfd::FileDialog::new()
            .set_file_name("project.gentle.json")
            .add_filter("GENtle project", &["json"])
            .save_file()
        {
            let path = path.display().to_string();
            if self.save_project_to_file(&path).is_ok() {
                self.current_project_path = Some(path);
                self.mark_clean_snapshot();
                return true;
            }
        }
        false
    }

    fn execute_project_action(&mut self, action: ProjectAction) {
        match action {
            ProjectAction::New => self.new_project(),
            ProjectAction::Open => self.prompt_open_project(),
            ProjectAction::Close => self.close_project(),
        }
    }

    fn request_project_action(&mut self, action: ProjectAction) {
        if self.is_project_dirty() {
            self.pending_project_action = Some(action);
        } else {
            self.execute_project_action(action);
        }
    }

    fn new_project(&mut self) {
        self.reset_to_empty_project();
    }

    fn close_project(&mut self) {
        self.reset_to_empty_project();
    }

    fn load_project_from_file(&mut self, path: &str) -> Result<()> {
        let state = ProjectState::load_from_path(path).map_err(|e| anyhow!(e.to_string()))?;

        self.engine = Arc::new(RwLock::new(GentleEngine::from_state(state)));
        self.current_project_path = Some(path.to_string());
        self.new_windows.clear();
        self.windows.clear();
        self.windows_to_close.write().unwrap().clear();

        self.mark_clean_snapshot();
        Ok(())
    }

    fn current_project_name(&self) -> String {
        match &self.current_project_path {
            Some(path) => Path::new(path)
                .file_name()
                .map(|s| s.to_string_lossy().to_string())
                .unwrap_or_else(|| path.clone()),
            None => "Untitled Project".to_string(),
        }
    }

    pub fn render_menu_bar(&mut self, ui: &mut Ui) {
        menu::bar(ui, |ui| {
            ui.menu_button(TRANSLATIONS.get("m_file"), |ui| {
                if ui.button("New Project").clicked() {
                    self.request_project_action(ProjectAction::New);
                    ui.close_menu();
                }
                if ui.button("Open Project...").clicked() {
                    self.request_project_action(ProjectAction::Open);
                    ui.close_menu();
                }
                if ui.button("Close Project").clicked() {
                    self.request_project_action(ProjectAction::Close);
                    ui.close_menu();
                }
                if ui.button("Open Sequence...").clicked() {
                    self.prompt_open_sequence();
                    ui.close_menu();
                }
                if ui.button("Import REBASE Data...").clicked() {
                    self.prompt_import_rebase_resource();
                    ui.close_menu();
                }
                if ui.button("Import JASPAR Data...").clicked() {
                    self.prompt_import_jaspar_resource();
                    ui.close_menu();
                }
                if ui.button("Save Project...").clicked() {
                    self.prompt_save_project();
                    ui.close_menu();
                }
                if ui.button("Export DALG SVG...").clicked() {
                    self.prompt_export_lineage_svg();
                    ui.close_menu();
                }
            });
            ui.menu_button("Help", |ui| {
                if ui.button("About GENtle").clicked() {
                    self.show_about_dialog = true;
                    ui.close_menu();
                }
            });
        });
    }

    fn show_window(&self, ctx: &egui::Context, id: ViewportId, window: Arc<RwLock<Window>>) {
        let windows_to_close = self.windows_to_close.clone();
        let window_number = self.get_window_number_from_id(id);
        let window_pos = Pos2 {
            x: window_number as f32 * 200.0,
            y: window_number as f32 * 200.0,
        };
        let window_title = window.read().unwrap().name();
        ctx.show_viewport_deferred(
            id,
            egui::ViewportBuilder::default()
                .with_title(window_title)
                // .with_maximized(true),
                .with_position(window_pos),
            move |ctx, class| {
                assert!(
                    class == egui::ViewportClass::Deferred,
                    "This egui backend doesn't support multiple viewports"
                );

                // Draw the window
                window.write().unwrap().update(ctx);

                // "Close window" action
                if ctx.input(|i| i.viewport().close_requested()) {
                    windows_to_close.write().unwrap().push(id);
                }
            },
        );
    }

    fn get_window_number_from_id(&self, id: ViewportId) -> usize {
        let window_number = self
            .windows
            .keys()
            .enumerate()
            .find(|(_num, viewport_id)| **viewport_id == id)
            .unwrap()
            .0;
        window_number
    }

    fn render_main_lineage(&mut self, ui: &mut Ui) {
        #[derive(Clone)]
        struct LineageRow {
            node_id: String,
            seq_id: String,
            display_name: String,
            origin: String,
            created_by_op: String,
            created_at: u128,
            parents: Vec<String>,
            length: usize,
            circular: bool,
            pool_size: usize,
            pool_members: Vec<String>,
        }

        #[derive(Clone)]
        struct ContainerRow {
            container_id: String,
            kind: String,
            member_count: usize,
            representative: String,
            members: Vec<String>,
        }

        let (rows, lineage_edges, op_label_by_id, containers) = {
            let engine = self.engine.read().unwrap();
            let state = engine.state();
            let mut op_created_count: HashMap<String, usize> = HashMap::new();
            let mut op_created_ids: HashMap<String, Vec<String>> = HashMap::new();
            let mut op_label_by_id: HashMap<String, String> = HashMap::new();
            for rec in engine.operation_log() {
                op_created_count.insert(rec.result.op_id.clone(), rec.result.created_seq_ids.len());
                op_created_ids.insert(rec.result.op_id.clone(), rec.result.created_seq_ids.clone());
                op_label_by_id.insert(rec.result.op_id.clone(), Self::summarize_operation(&rec.op));
            }

            let mut out: Vec<LineageRow> = state
                .lineage
                .nodes
                .values()
                .map(|node| {
                    let parents: Vec<String> = state
                        .lineage
                        .edges
                        .iter()
                        .filter(|e| e.to_node_id == node.node_id)
                        .filter_map(|e| state.lineage.nodes.get(&e.from_node_id))
                        .map(|p| p.seq_id.clone())
                        .collect();
                    let (length, circular) = state
                        .sequences
                        .get(&node.seq_id)
                        .map(|dna| (dna.len(), dna.is_circular()))
                        .unwrap_or((0, false));
                    let display_name = state
                        .sequences
                        .get(&node.seq_id)
                        .and_then(|dna| dna.name().clone())
                        .unwrap_or_else(|| node.seq_id.clone());
                    let pool_size = node
                        .created_by_op
                        .as_ref()
                        .and_then(|op| op_created_count.get(op))
                        .cloned()
                        .unwrap_or(1);
                    let pool_members = node
                        .created_by_op
                        .as_ref()
                        .and_then(|op| op_created_ids.get(op))
                        .cloned()
                        .unwrap_or_else(|| vec![node.seq_id.clone()]);
                    LineageRow {
                        node_id: node.node_id.clone(),
                        seq_id: node.seq_id.clone(),
                        display_name,
                        origin: format!("{:?}", node.origin),
                        created_by_op: node
                            .created_by_op
                            .clone()
                            .unwrap_or_else(|| "-".to_string()),
                        created_at: node.created_at_unix_ms,
                        parents,
                        length,
                        circular,
                        pool_size,
                        pool_members,
                    }
                })
                .collect();
            out.sort_by(|a, b| {
                a.created_at
                    .cmp(&b.created_at)
                    .then(a.node_id.cmp(&b.node_id))
            });
            let lineage_edges: Vec<(String, String, String)> = state
                .lineage
                .edges
                .iter()
                .filter_map(|e| {
                    let from = state.lineage.nodes.get(&e.from_node_id)?.node_id.clone();
                    let to = state.lineage.nodes.get(&e.to_node_id)?.node_id.clone();
                    Some((from, to, e.op_id.clone()))
                })
                .collect();
            let mut containers: Vec<ContainerRow> = state
                .container_state
                .containers
                .iter()
                .map(|(id, c)| ContainerRow {
                    container_id: id.clone(),
                    kind: format!("{:?}", c.kind),
                    member_count: c.members.len(),
                    representative: c.members.first().cloned().unwrap_or_default(),
                    members: c.members.clone(),
                })
                .collect();
            containers.sort_by(|a, b| a.container_id.cmp(&b.container_id));
            (out, lineage_edges, op_label_by_id, containers)
        };

        ui.heading("Lineage Graph");
        ui.label(format!("Project: {}", self.current_project_name()));
        ui.label("Project-level sequence lineage (branch and merge aware)");
        ui.horizontal(|ui| {
            let table = ui
                .selectable_label(!self.lineage_graph_view, "Table")
                .on_hover_text("Show lineage as table");
            if table.clicked() {
                self.lineage_graph_view = false;
            }
            let graph = ui
                .selectable_label(self.lineage_graph_view, "Graph")
                .on_hover_text("Show lineage as node graph");
            if graph.clicked() {
                self.lineage_graph_view = true;
            }
        });
        ui.separator();
        ui.horizontal(|ui| {
            ui.label("Legend:");
            ui.colored_label(egui::Color32::from_rgb(90, 140, 210), "● single sequence");
            ui.colored_label(
                egui::Color32::from_rgb(180, 120, 70),
                "◆ pool (n = expected variants)",
            );
        });
        ui.separator();

        let mut open_seq: Option<String> = None;
        let mut open_pool: Option<(String, Vec<String>)> = None;
        let mut select_candidate_from: Option<String> = None;
        if self.lineage_graph_view {
            egui::ScrollArea::both().show(ui, |ui| {
                let width = (rows.len().max(1) as f32) * 180.0 + 120.0;
                let height = 440.0;
                let (resp, painter) =
                    ui.allocate_painter(Vec2::new(width, height), egui::Sense::hover());
                let rect = resp.rect;
                let mut pos_by_node: HashMap<String, Pos2> = HashMap::new();
                for (idx, row) in rows.iter().enumerate() {
                    let mut h = std::collections::hash_map::DefaultHasher::new();
                    row.node_id.hash(&mut h);
                    let lane = (h.finish() % 4) as f32;
                    let x = rect.left() + 80.0 + idx as f32 * 170.0;
                    let y = rect.top() + 70.0 + lane * 80.0;
                    pos_by_node.insert(row.node_id.clone(), Pos2::new(x, y));
                }
                let mut used_label_rects: Vec<egui::Rect> = Vec::new();
                for (from_node, to_node, op_id) in &lineage_edges {
                    let Some(from) = pos_by_node.get(from_node) else {
                        continue;
                    };
                    let Some(to) = pos_by_node.get(to_node) else {
                        continue;
                    };
                    painter.line_segment([*from, *to], egui::Stroke::new(1.0, egui::Color32::GRAY));
                    let mid = Pos2::new((from.x + to.x) * 0.5, (from.y + to.y) * 0.5);
                    let op_label = op_label_by_id
                        .get(op_id)
                        .cloned()
                        .unwrap_or_else(|| op_id.clone());
                    let display = op_label;
                    let galley = painter.layout_no_wrap(
                        display.clone(),
                        egui::FontId::proportional(10.0),
                        egui::Color32::BLACK,
                    );
                    let edge = *to - *from;
                    let edge_len = edge.length();
                    if edge_len < 0.1 {
                        continue;
                    }
                    let edge_dir = edge / edge_len;
                    let perp = Vec2::new(-edge_dir.y, edge_dir.x);
                    let mut placed = None;
                    for idx in 0..14 {
                        let sign = if idx % 2 == 0 { 1.0 } else { -1.0 };
                        let step = (idx / 2) as f32;
                        let candidate_center = mid + perp * (12.0 + step * 12.0) * sign;
                        let candidate_rect = egui::Rect::from_center_size(
                            candidate_center,
                            Vec2::new(galley.size().x + 10.0, galley.size().y + 6.0),
                        );
                        if !used_label_rects
                            .iter()
                            .any(|r| r.intersects(candidate_rect))
                        {
                            placed = Some((candidate_center, candidate_rect));
                            used_label_rects.push(candidate_rect);
                            break;
                        }
                    }
                    let (label_center, bg_rect) = placed.unwrap_or_else(|| {
                        let fallback_center = mid + perp * 12.0;
                        let rect = egui::Rect::from_center_size(
                            fallback_center,
                            Vec2::new(galley.size().x + 10.0, galley.size().y + 6.0),
                        );
                        (fallback_center, rect)
                    });
                    used_label_rects.push(bg_rect);
                    painter.rect_filled(
                        bg_rect,
                        3.0,
                        egui::Color32::from_rgba_premultiplied(245, 245, 245, 235),
                    );
                    painter.text(
                        label_center,
                        egui::Align2::CENTER_CENTER,
                        display,
                        egui::FontId::proportional(10.0),
                        egui::Color32::BLACK,
                    );
                }
                for row in &rows {
                    let Some(pos) = pos_by_node.get(&row.node_id).cloned() else {
                        continue;
                    };
                    if row.pool_size > 1 {
                        let points = vec![
                            pos + Vec2::new(0.0, -16.0),
                            pos + Vec2::new(16.0, 0.0),
                            pos + Vec2::new(0.0, 16.0),
                            pos + Vec2::new(-16.0, 0.0),
                        ];
                        painter.add(egui::Shape::convex_polygon(
                            points,
                            egui::Color32::from_rgb(180, 120, 70),
                            egui::Stroke::new(1.0, egui::Color32::from_rgb(235, 196, 150)),
                        ));
                        painter.text(
                            pos + Vec2::new(19.0, -14.0),
                            egui::Align2::LEFT_TOP,
                            format!("n={}", row.pool_size),
                            egui::FontId::proportional(10.0),
                            egui::Color32::YELLOW,
                        );
                    } else {
                        painter.circle_filled(pos, 16.0, egui::Color32::from_rgb(90, 140, 210));
                    }
                    painter.text(
                        pos + Vec2::new(22.0, -4.0),
                        egui::Align2::LEFT_BOTTOM,
                        &row.display_name,
                        egui::FontId::proportional(12.0),
                        egui::Color32::BLACK,
                    );
                    painter.text(
                        pos + Vec2::new(22.0, 10.0),
                        egui::Align2::LEFT_TOP,
                        format!("{} ({} bp)", row.seq_id, row.length),
                        egui::FontId::proportional(10.0),
                        egui::Color32::BLACK,
                    );
                }

                // Interactions: single-click opens sequence, double-click selects candidate.
                if let Some(pointer) = resp.interact_pointer_pos() {
                    let mut hit_row: Option<LineageRow> = None;
                    for row in &rows {
                        let Some(pos) = pos_by_node.get(&row.node_id).cloned() else {
                            continue;
                        };
                        let hit = if row.pool_size > 1 {
                            egui::Rect::from_center_size(pos, Vec2::new(30.0, 24.0))
                                .contains(pointer)
                        } else {
                            pointer.distance(pos) <= 18.0
                        };
                        if hit {
                            hit_row = Some(row.clone());
                            break;
                        }
                    }
                    if let Some(row) = hit_row {
                        if resp.double_clicked() {
                            select_candidate_from = Some(row.seq_id.clone());
                        } else if resp.clicked() {
                            if row.pool_size > 1 {
                                open_pool = Some((row.seq_id.clone(), row.pool_members.clone()));
                            } else {
                                open_seq = Some(row.seq_id.clone());
                            }
                        }
                    }
                }
            });
        } else {
            egui::ScrollArea::vertical().show(ui, |ui| {
                egui::Grid::new("lineage_grid")
                    .striped(true)
                    .show(ui, |ui| {
                        ui.strong("Node");
                        ui.strong("Sequence");
                        ui.strong("Parents");
                        ui.strong("Origin");
                        ui.strong("Op");
                        ui.strong("Length");
                        ui.strong("Topology");
                        ui.strong("Action");
                        ui.end_row();

                        for row in &rows {
                            ui.monospace(&row.node_id);
                            if ui.button(&row.seq_id).clicked() {
                                open_seq = Some(row.seq_id.clone());
                            }
                            ui.label(if row.parents.is_empty() {
                                "-".to_string()
                            } else {
                                row.parents.join(" + ")
                            });
                            ui.label(&row.origin);
                            ui.monospace(&row.created_by_op);
                            ui.monospace(format!("{} bp", row.length));
                            ui.label(if row.circular { "circular" } else { "linear" });
                            if ui.button("Select").clicked() {
                                select_candidate_from = Some(row.seq_id.clone());
                            }
                            ui.end_row();
                        }
                    });
            });
        }
        ui.separator();
        ui.heading("Containers");
        ui.label("Container-level view of candidate sequence sets");
        egui::ScrollArea::vertical()
            .max_height(180.0)
            .show(ui, |ui| {
                egui::Grid::new("container_grid")
                    .striped(true)
                    .show(ui, |ui| {
                        ui.strong("Container");
                        ui.strong("Kind");
                        ui.strong("Members");
                        ui.strong("Representative");
                        ui.strong("Action");
                        ui.end_row();
                        for c in &containers {
                            ui.monospace(&c.container_id);
                            ui.label(&c.kind);
                            ui.monospace(format!("{}", c.member_count));
                            ui.monospace(&c.representative);
                            if c.member_count > 1 {
                                if ui.button("Open Pool").clicked() {
                                    open_pool = Some((c.representative.clone(), c.members.clone()));
                                }
                            } else if !c.representative.is_empty() {
                                if ui.button("Open Seq").clicked() {
                                    open_seq = Some(c.representative.clone());
                                }
                            } else {
                                ui.label("-");
                            }
                            ui.end_row();
                        }
                    });
            });

        if let Some(input) = select_candidate_from {
            let criterion = format!("gui_lineage_select:{input}");
            let result = self
                .engine
                .write()
                .unwrap()
                .apply(Operation::SelectCandidate {
                    input: input.clone(),
                    criterion,
                    output_id: None,
                });
            if let Ok(op_result) = result {
                if let Some(seq_id) = op_result.created_seq_ids.first() {
                    open_seq = Some(seq_id.clone());
                }
            }
        }

        if let Some(seq_id) = open_seq {
            self.open_sequence_window(&seq_id);
        }
        if let Some((representative, pool_members)) = open_pool {
            self.open_pool_window(&representative, pool_members);
        }
    }

    fn render_unsaved_changes_dialog(&mut self, ctx: &egui::Context) {
        if self.pending_project_action.is_none() {
            return;
        }
        egui::Window::new("Unsaved Changes")
            .collapsible(false)
            .resizable(false)
            .show(ctx, |ui| {
                ui.label("Save changes to the current project before continuing?");
                ui.horizontal(|ui| {
                    if ui.button("Save").clicked() {
                        if self.save_current_project() {
                            if let Some(action) = self.pending_project_action.take() {
                                self.execute_project_action(action);
                            }
                        }
                    }
                    if ui.button("Don't Save").clicked() {
                        if let Some(action) = self.pending_project_action.take() {
                            self.execute_project_action(action);
                        }
                    }
                    if ui.button("Cancel").clicked() {
                        self.pending_project_action = None;
                    }
                });
            });
    }

    fn render_about_dialog(&mut self, ctx: &egui::Context) {
        if !self.show_about_dialog {
            return;
        }
        egui::Window::new("About GENtle")
            .open(&mut self.show_about_dialog)
            .collapsible(false)
            .resizable(false)
            .show(ctx, |ui| {
                ui.vertical_centered(|ui| {
                    ui.add(APP_ICON.clone().fit_to_exact_size(Vec2::new(96.0, 96.0)));
                    for line in about::version_cli_text().lines() {
                        if line.starts_with("GENtle ") {
                            ui.heading(line);
                        } else {
                            ui.label(line);
                        }
                    }
                });
            });
    }

    fn summarize_operation(op: &Operation) -> String {
        match op {
            Operation::LoadFile { path, as_id } => match as_id {
                Some(id) => format!("Load file: path={path}, as_id={id}"),
                None => format!("Load file: path={path}"),
            },
            Operation::DigestContainer {
                container_id,
                enzymes,
                output_prefix,
            } => format!(
                "Digest container: container_id={container_id}, enzymes=[{}], output_prefix={}",
                enzymes.join(", "),
                output_prefix.clone().unwrap_or_else(|| "-".to_string())
            ),
            Operation::MergeContainersById {
                container_ids,
                output_prefix,
            } => format!(
                "Merge containers by id: container_ids={}, output_prefix={}",
                container_ids.join(", "),
                output_prefix.clone().unwrap_or_else(|| "-".to_string())
            ),
            Operation::LigationContainer {
                container_id,
                circularize_if_possible,
                protocol,
                output_prefix,
                unique,
                ..
            } => format!(
                "Ligation container: container_id={container_id}, protocol={:?}, circularize_if_possible={}, output_prefix={}, unique={}",
                protocol,
                circularize_if_possible,
                output_prefix.clone().unwrap_or_else(|| "-".to_string()),
                unique.unwrap_or(false)
            ),
            Operation::FilterContainerByMolecularWeight {
                container_id,
                min_bp,
                max_bp,
                error,
                unique,
                ..
            } => format!(
                "Molecular weight filter (container): container_id={}, min_bp={}, max_bp={}, error={:.2}, unique={}",
                container_id, min_bp, max_bp, error, unique
            ),
            Operation::Digest {
                input,
                enzymes,
                output_prefix,
            } => {
                format!(
                    "Digest: input={input}, enzymes=[{}], output_prefix={}",
                    enzymes.join(", "),
                    output_prefix.clone().unwrap_or_else(|| "-".to_string())
                )
            }
            Operation::MergeContainers {
                inputs,
                output_prefix,
            } => format!(
                "Merge containers: inputs={}, output_prefix={}",
                inputs.join(", "),
                output_prefix.clone().unwrap_or_else(|| "-".to_string())
            ),
            Operation::Ligation {
                inputs,
                circularize_if_possible,
                protocol,
                output_prefix,
                unique,
                ..
            } => format!(
                "Ligation: inputs={}, protocol={:?}, circularize_if_possible={}, output_prefix={}, unique={}",
                inputs.join(", "),
                protocol,
                circularize_if_possible,
                output_prefix.clone().unwrap_or_else(|| "-".to_string()),
                unique.unwrap_or(false)
            ),
            Operation::Pcr {
                template,
                forward_primer,
                reverse_primer,
                unique,
                ..
            } => format!(
                "PCR: template={template}, forward_primer={}, reverse_primer={}, unique={}",
                forward_primer,
                reverse_primer,
                unique.unwrap_or(false)
            ),
            Operation::PcrAdvanced {
                template,
                forward_primer,
                reverse_primer,
                unique,
                ..
            } => format!(
                "PCR advanced: template={template}, fwd_len={}, fwd_anneal={:?}, fwd_mm={:?}, rev_len={}, rev_anneal={:?}, rev_mm={:?}, unique={}",
                forward_primer.sequence.len(),
                forward_primer.anneal_len,
                forward_primer.max_mismatches,
                reverse_primer.sequence.len(),
                reverse_primer.anneal_len,
                reverse_primer.max_mismatches,
                unique.unwrap_or(false)
            ),
            Operation::PcrMutagenesis {
                template,
                forward_primer,
                reverse_primer,
                mutations,
                require_all_mutations,
                unique,
                ..
            } => format!(
                "PCR mutagenesis: template={template}, fwd_len={}, rev_len={}, mutations={}, require_all_mutations={}, unique={}",
                forward_primer.sequence.len(),
                reverse_primer.sequence.len(),
                mutations.len(),
                require_all_mutations.unwrap_or(false),
                unique.unwrap_or(false)
            ),
            Operation::FilterByMolecularWeight {
                inputs,
                min_bp,
                max_bp,
                error,
                unique,
                ..
            } => format!(
                "Molecular weight filter: inputs={}, min_bp={}, max_bp={}, error={:.2}, unique={}",
                inputs.join(", "),
                min_bp,
                max_bp,
                error,
                unique
            ),
            Operation::SelectCandidate {
                input,
                criterion,
                output_id,
            } => format!(
                "Select candidate: input={input}, criterion={criterion}, output_id={}",
                output_id.clone().unwrap_or_else(|| "-".to_string())
            ),
            Operation::Reverse { input, output_id } => format!(
                "Reverse: input={input}, output_id={}",
                output_id.clone().unwrap_or_else(|| "-".to_string())
            ),
            Operation::Complement { input, output_id } => format!(
                "Complement: input={input}, output_id={}",
                output_id.clone().unwrap_or_else(|| "-".to_string())
            ),
            Operation::ReverseComplement { input, output_id } => format!(
                "Reverse complement: input={input}, output_id={}",
                output_id.clone().unwrap_or_else(|| "-".to_string())
            ),
            Operation::Branch { input, output_id } => format!(
                "Branch: input={input}, output_id={}",
                output_id.clone().unwrap_or_else(|| "-".to_string())
            ),
            Operation::SetDisplayVisibility { target, visible } => {
                let target_name = match target {
                    DisplayTarget::SequencePanel => "Sequence panel",
                    DisplayTarget::MapPanel => "Map panel",
                    DisplayTarget::Features => "Features",
                    DisplayTarget::Tfbs => "TFBS",
                    DisplayTarget::RestrictionEnzymes => "Restriction enzymes",
                    DisplayTarget::GcContents => "GC contents",
                    DisplayTarget::OpenReadingFrames => "Open reading frames",
                    DisplayTarget::MethylationSites => "Methylation sites",
                };
                format!("Set display visibility: target={target_name}, visible={visible}")
            }
            Operation::SetTopology { seq_id, circular } => {
                if *circular {
                    format!("Set topology: seq_id={seq_id}, topology=circular")
                } else {
                    format!("Set topology: seq_id={seq_id}, topology=linear")
                }
            }
            Operation::RecomputeFeatures { seq_id } => {
                format!("Recompute features: seq_id={seq_id}")
            }
            Operation::SetParameter { name, value } => {
                format!("Set parameter: name={name}, value={value}")
            }
            Operation::AnnotateTfbs {
                seq_id,
                motifs,
                min_llr_bits,
                min_llr_quantile,
                per_tf_thresholds,
                clear_existing,
                max_hits,
            } => format!(
                "Annotate TFBS: seq_id={seq_id}, motifs=[{}], min_llr_bits={:?}, min_llr_quantile={:?}, per_tf_overrides={}, clear_existing={}, max_hits={}",
                motifs.join(", "),
                min_llr_bits,
                min_llr_quantile,
                per_tf_thresholds.len(),
                clear_existing.unwrap_or(true),
                max_hits
                    .map(|v| v.to_string())
                    .unwrap_or_else(|| "default(500)".to_string())
            ),
            Operation::ExtractRegion {
                input,
                from,
                to,
                output_id,
            } => format!(
                "Extract region: input={input}, from={from}, to={to}, output_id={}",
                output_id.clone().unwrap_or_else(|| "-".to_string())
            ),
            Operation::ExtractAnchoredRegion {
                input,
                anchor,
                direction,
                target_length_bp,
                length_tolerance_bp,
                required_re_sites,
                required_tf_motifs,
                output_prefix,
                unique,
                max_candidates,
                ..
            } => format!(
                "Extract anchored region: input={input}, anchor={anchor:?}, direction={direction:?}, target_length_bp={}, tolerance_bp={}, re_sites=[{}], tf_motifs=[{}], output_prefix={}, unique={}, max_candidates={}",
                target_length_bp,
                length_tolerance_bp,
                required_re_sites.join(", "),
                required_tf_motifs.join(", "),
                output_prefix.clone().unwrap_or_else(|| "-".to_string()),
                unique.unwrap_or(false),
                max_candidates
                    .map(|v| v.to_string())
                    .unwrap_or_else(|| "-".to_string())
            ),
            Operation::SaveFile {
                seq_id,
                path,
                format,
            } => format!("Save file: seq_id={seq_id}, path={path}, format={format:?}"),
            Operation::RenderSequenceSvg { seq_id, mode, path } => format!(
                "Render sequence SVG: seq_id={seq_id}, mode={mode:?}, path={path}"
            ),
            Operation::RenderLineageSvg { path } => {
                format!("Render lineage SVG: path={path}")
            }
            Operation::ExportPool {
                inputs,
                path,
                pool_id,
                human_id,
            } => format!(
                "Export pool: inputs={}, path={}, pool_id={}, human_id={}",
                inputs.join(", "),
                path,
                pool_id.clone().unwrap_or_else(|| "-".to_string()),
                human_id.clone().unwrap_or_else(|| "-".to_string())
            ),
        }
    }
}

impl eframe::App for GENtleApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        if !self.update_has_run_before {
            egui_extras::install_image_loaders(ctx);
            self.update_has_run_before = true;
        }
        let dirty_marker = if self.is_project_dirty() { " *" } else { "" };
        ctx.send_viewport_cmd(egui::ViewportCommand::Title(format!(
            "GENtle - {}{} (v{})",
            self.current_project_name(),
            dirty_marker,
            about::GENTLE_DISPLAY_VERSION
        )));

        let open_project = KeyboardShortcut::new(Modifiers::COMMAND, Key::O);
        let new_project = KeyboardShortcut::new(Modifiers::COMMAND, Key::N);
        let open_sequence = KeyboardShortcut::new(Modifiers::COMMAND | Modifiers::SHIFT, Key::O);
        let save_project = KeyboardShortcut::new(Modifiers::COMMAND, Key::S);
        let close_project = KeyboardShortcut::new(Modifiers::COMMAND | Modifiers::SHIFT, Key::W);
        if ctx.input_mut(|i| i.consume_shortcut(&new_project)) {
            self.request_project_action(ProjectAction::New);
        }
        if ctx.input_mut(|i| i.consume_shortcut(&open_project)) {
            self.request_project_action(ProjectAction::Open);
        }
        if ctx.input_mut(|i| i.consume_shortcut(&open_sequence)) {
            self.prompt_open_sequence();
        }
        if ctx.input_mut(|i| i.consume_shortcut(&save_project)) {
            let _ = self.save_current_project();
        }
        if ctx.input_mut(|i| i.consume_shortcut(&close_project)) {
            self.request_project_action(ProjectAction::Close);
        }

        // Show menu bar
        egui::TopBottomPanel::top("top").show(ctx, |ui| {
            self.render_menu_bar(ui);
        });

        // Show main window
        egui::CentralPanel::default().show(ctx, |ui| {
            self.render_main_lineage(ui);
            if self.is_project_dirty() {
                ui.label("Status: unsaved changes");
            } else {
                ui.label("Status: saved");
            }
        });
        self.render_about_dialog(ctx);
        self.render_unsaved_changes_dialog(ctx);

        // Open new windows
        let mut new_windows: Vec<Window> = self.new_windows.drain(..).collect();
        for window in new_windows.drain(..) {
            let id = format!("Viewport {}", self.viewport_id_counter);
            let id = ViewportId::from_hash_of(id);
            self.windows.insert(id, Arc::new(RwLock::new(window)));
            self.viewport_id_counter += 1;
        }

        // Close windows
        for id in self.windows_to_close.write().unwrap().drain(..) {
            self.windows.remove(&id);
        }

        // Show windows
        for (id, window) in self.windows.iter() {
            let id = id.to_owned();
            self.show_window(ctx, id, window.clone());
        }
    }
}
