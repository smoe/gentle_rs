//! Help/tutorial markdown loading and rewriting helpers for `GENtleApp`.
//!
//! This submodule keeps the top-level app file focused on window/session state
//! while the help-system path owns its markdown parsing, tutorial discovery,
//! and relative-image rewriting in one place.
//!
//! Look here for:
//! - markdown cleanup/rewrite helpers used before rendering docs in-app
//! - tutorial/help discovery and path resolution
//! - relative asset/image rewriting and in-app chapter-link navigation

use super::*;

#[cfg(test)]
#[path = "help_docs/tests.rs"]
mod tests;

pub(super) const AGENT_INTERFACES_TUTORIAL_PATH: &str = "docs/tutorial/01-01_agent_interfaces.md";
pub(super) const AGENT_INTERFACES_TUTORIAL_TITLE: &str =
    "GENtle Agent Assistant and Agent Interfaces Tutorial";
pub(super) const AGENT_INTERFACES_TUTORIAL_SUMMARY: &str = "docs/tutorial/01-01_agent_interfaces.md\n\
type: operational_reference\n\
status: manual/reference\n\
Practical guide for the in-app Agent Assistant, provider quick starts, reviewed shared-shell suggestions, CLI/shared shell, MCP, and external coding agents.";

impl GENtleApp {
    pub(super) fn open_help_for_active_window(&mut self) {
        let topic = self
            .open_window_model_cache
            .model
            .entries
            .iter()
            .find(|entry| Some(entry.native_menu_key) == self.active_window_menu_key)
            .map_or(HelpTopic::MainWindow, |entry| entry.help_topic);
        self.open_help_topic(topic);
    }

    pub(super) fn open_help_topic(&mut self, topic: HelpTopic) {
        let (id, fragment) = topic.destination();
        self.ensure_help_docs_loaded();
        let result = super::help_manual_sources::GUI_SOURCES
            .iter()
            .find(|source| source.id == id)
            .ok_or_else(|| format!("No GUI manual topic for '{id}'."))
            .and_then(|source| self.open_gui_manual_target(source.path, fragment));
        if let Err(message) = result {
            self.open_help_doc(HelpDoc::Gui);
            self.help_link_notice = Some(message);
        }
    }

    pub(super) fn consume_manual_help_request(&mut self, ctx: &egui::Context) {
        if let Some(topic) = super::help_topics::take_request(ctx) {
            self.open_help_topic(topic);
        }
    }

    fn manual_book_at(doc: HelpDoc, root: Option<&Path>) -> (super::help_manual::ManualBook, Option<String>) {
        let sources = match doc {
            HelpDoc::Gui => super::help_manual_sources::GUI_SOURCES,
            HelpDoc::Cli => super::help_manual_sources::CLI_SOURCES,
            _ => unreachable!("Only GUI and CLI are catalogued manuals"),
        };
        let Some(root) = root else {
            return (super::help_manual::ManualBook::embedded(sources), None);
        };
        let load = || -> Option<super::help_manual::ManualBook> {
            let catalog = fs::read_to_string(root.join("docs/manuals.json")).ok()?;
            let actual: serde_json::Value = serde_json::from_str(&catalog).ok()?;
            let expected: serde_json::Value = serde_json::from_str(super::help_manual_sources::CATALOG_JSON).ok()?;
            if actual != expected {
                return None;
            }
            let texts = sources.iter().map(|source| fs::read_to_string(root.join(source.path)).ok())
                .collect::<Option<Vec<_>>>()?;
            Some(super::help_manual::ManualBook::from_texts(sources, texts, Some(root.to_path_buf())))
        };
        match load() {
            Some(book) => (book, None),
            None => (
                super::help_manual::ManualBook::embedded(sources),
                Some(format!("External {doc:?} manual is incomplete or has a different catalog. Using the complete embedded manual; no external and embedded chapters were mixed.")),
            ),
        }
    }

    pub(super) fn load_manual_books(&mut self) {
        let root = Self::resolve_runtime_doc_path("docs/manuals.json")
            .and_then(|path| path.canonicalize().ok())
            .and_then(|path| path.parent()?.parent().map(Path::to_path_buf));
        let (book, gui_notice) = Self::manual_book_at(HelpDoc::Gui, root.as_deref());
        self.help_gui_markdown.clone_from(&book.markdown);
        self.help_gui_book = Some(Arc::new(book));
        let (book, cli_notice) = Self::manual_book_at(HelpDoc::Cli, root.as_deref());
        self.help_cli_markdown.clone_from(&book.markdown);
        self.help_cli_book = Some(Arc::new(book));
        let notices: Vec<_> = [gui_notice, cli_notice].into_iter().flatten().collect();
        self.help_manual_load_notice = (!notices.is_empty()).then(|| notices.join("\n"));
        self.help_gui_pages = None;
        self.help_cli_pages = None;
        self.help_manual_back.clear();
        self.help_manual_forward.clear();
        self.help_rendered_markdown_cache.retain(|key, _| !key.starts_with("gui:") && !key.starts_with("cli:"));
        self.help_markdown_cache = CommonMarkCache::default();
    }

    fn render_manual_groups(
        ui: &mut Ui,
        pages: &super::help_manual::ManualPages,
        selected: &mut usize,
        reveal: bool,
    ) {
        let mut start = 0;
        while start < pages.pages.len() {
            let group = &pages.pages[start].group;
            let end = (start + 1..pages.pages.len())
                .find(|&index| pages.pages[index].group != *group)
                .unwrap_or(pages.pages.len());
            if group.is_empty() {
                Self::render_manual_tree(ui, pages, start..end, selected, reveal);
            } else {
                let id = ui.make_persistent_id(("manual_group", group));
                let mut state = egui::collapsing_header::CollapsingState::load_with_default_open(
                    ui.ctx(), id, (start..end).contains(selected),
                );
                if reveal && (start..end).contains(selected) {
                    state.set_open(true);
                }
                state.show_header(ui, |ui| { ui.label(group); }).body(|ui| {
                    Self::render_manual_tree(ui, pages, start..end, selected, reveal);
                });
            }
            start = end;
        }
    }

    fn render_manual_tree(
        ui: &mut Ui,
        pages: &super::help_manual::ManualPages,
        range: std::ops::Range<usize>,
        selected: &mut usize,
        reveal: bool,
    ) {
        let mut index = range.start;
        while index < range.end {
            let page = &pages.pages[index];
            let level = page.heading_level + u8::from(page.anchor.is_none());
            let end = (index + 1..range.end)
                .find(|&next| {
                    let next = &pages.pages[next];
                    next.heading_level + u8::from(next.anchor.is_none()) <= level
                })
                .unwrap_or(range.end);
            if end > index + 1 {
                let id = ui.make_persistent_id(("manual_tree", pages.source_hash, index));
                let mut state = egui::collapsing_header::CollapsingState::load_with_default_open(
                    ui.ctx(), id, level == 1,
                );
                if reveal && (index..end).contains(selected) {
                    state.set_open(true);
                }
                state.show_header(ui, |ui| {
                    if ui.add(egui::Button::new(&page.title).selected(*selected == index).wrap()).clicked() {
                        *selected = index;
                    }
                }).body(|ui| Self::render_manual_tree(ui, pages, index + 1..end, selected, reveal));
            } else if ui.add(egui::Button::new(&page.title).selected(*selected == index).wrap()).clicked() {
                *selected = index;
            }
            index = end;
        }
    }

    pub(super) fn render_help_navigation_and_page(&mut self, ui: &mut Ui, size: Vec2) {
        if matches!(self.help_doc, HelpDoc::Gui | HelpDoc::Cli) {
            let doc = self.help_doc;
            let pages = self.ensure_manual_pages(doc);
            let mut selected = self.manual_page(doc);
            let reveal = self.help_scroll_to_top;
            if size.x >= 720.0 {
                ui.horizontal_top(|ui| {
                    let tree_width = (size.x * 0.28).clamp(180.0, 280.0);
                    ui.allocate_ui_with_layout(Vec2::new(tree_width, size.y), egui::Layout::top_down(egui::Align::Min), |ui| {
                        ui.set_max_width(tree_width);
                        egui::ScrollArea::vertical().id_salt("manual_tree_scroll").max_height(size.y).show(ui, |ui| {
                            Self::render_manual_groups(ui, &pages, &mut selected, reveal);
                        });
                    });
                    ui.separator();
                    self.select_manual_location(doc, selected, true);
                    self.render_help_page(ui, Vec2::new(ui.available_width(), size.y));
                });
                return;
            }
            egui::CollapsingHeader::new("Contents").id_salt("manual_tree_compact").default_open(true).show(ui, |ui| {
                egui::ScrollArea::vertical().id_salt("manual_tree_compact_scroll").max_height(size.y * 0.3).show(ui, |ui| {
                    Self::render_manual_groups(ui, &pages, &mut selected, reveal);
                });
            });
            self.select_manual_location(doc, selected, true);
        }
        let height = ui.available_height().min(size.y).max(40.0);
        self.render_help_page(ui, Vec2::new(ui.available_width(), height));
    }

    fn render_help_page(&mut self, ui: &mut Ui, size: Vec2) {
        let content_width = size.x.max(HELP_BODY_MIN_WIDTH_PX);
        if Self::help_content_width_requires_relayout(self.help_last_content_width, content_width) {
            self.help_markdown_cache = CommonMarkCache::default();
        }
        self.help_last_content_width = content_width;
        let active_doc = self.help_doc;
        let source_hash = Self::help_markdown_hash(self.help_source_markdown(active_doc));
        let scroll_source_id = ("help_body_scroll", self.help_render_cache_key(active_doc), source_hash);
        let mut scroll = egui::ScrollArea::vertical().id_salt(scroll_source_id)
            .auto_shrink([false, false]).max_height(size.y);
        if std::mem::take(&mut self.help_scroll_to_top) {
            scroll = scroll.vertical_scroll_offset(0.0);
        }
        ui.set_width(content_width);
        ui.set_max_width(content_width);
        scroll.show(ui, |ui| {
            scroll_input_policy::apply_scrollarea_keyboard_navigation(
                ui, scroll_input_policy::DEFAULT_SCROLLAREA_KEYBOARD_STEP,
            );
            if self.help_selectable_text_mode {
                let mut text = self.help_source_markdown(active_doc).to_string();
                ui.add(egui::TextEdit::multiline(&mut text).font(egui::TextStyle::Monospace)
                    .desired_width(content_width).lock_focus(true));
            } else {
                let rendered = self.rendered_help_markdown_for(active_doc);
                let command_start = ui.ctx().output(|output| output.commands.len());
                CommonMarkViewer::new()
                    .max_image_width(Some(Self::help_markdown_max_image_width(content_width)))
                    .show(ui, &mut self.help_markdown_cache, rendered.as_ref());
                self.follow_help_markdown_links(ui.ctx(), command_start);
            }
        });
    }

    #[cfg(test)]
    pub(super) fn ensure_gui_manual_pages(&mut self) -> Arc<super::help_manual::ManualPages> {
        self.ensure_manual_pages(HelpDoc::Gui)
    }

    pub(super) fn ensure_manual_pages(&mut self, doc: HelpDoc) -> Arc<super::help_manual::ManualPages> {
        let (text, book, pages, selected) = match doc {
            HelpDoc::Gui => (&self.help_gui_markdown, &self.help_gui_book, &mut self.help_gui_pages, &mut self.help_gui_page),
            HelpDoc::Cli => (&self.help_cli_markdown, &self.help_cli_book, &mut self.help_cli_pages, &mut self.help_cli_page),
            _ => unreachable!("Only GUI and CLI are catalogued manuals"),
        };
        let source_hash = Self::help_markdown_hash(text);
        if pages.as_ref().is_none_or(|index| index.source_hash != source_hash) {
            let changed = pages.is_some();
            let index = match book.as_ref().filter(|book| book.source_hash == source_hash) {
                Some(book) => super::help_manual::ManualPages::from_book(book),
                None => {
                    let mut index = super::help_manual::ManualPages::new(text, source_hash);
                    if doc == HelpDoc::Cli {
                        for page in &mut index.pages {
                            page.document_path = "docs/cli.md".to_string();
                        }
                    }
                    index
                }
            };
            *selected = (*selected).min(index.pages.len() - 1);
            *pages = Some(Arc::new(index));
            if changed {
                self.help_manual_back.clear();
                self.help_manual_forward.clear();
            }
            self.help_rendered_markdown_cache
                .retain(|key, _| !key.starts_with(if doc == HelpDoc::Gui { "gui:" } else { "cli:" }));
            self.help_markdown_cache = CommonMarkCache::default();
        }
        pages.as_ref().unwrap().clone()
    }

    #[cfg(test)]
    pub(super) fn select_gui_manual_page(&mut self, index: usize) {
        self.select_manual_location(HelpDoc::Gui, index, true);
    }

    pub(super) fn manual_page(&self, doc: HelpDoc) -> usize {
        if doc == HelpDoc::Cli { self.help_cli_page } else { self.help_gui_page }
    }

    pub(super) fn select_manual_location(&mut self, doc: HelpDoc, index: usize, record_history: bool) {
        let pages = self.ensure_manual_pages(doc);
        let index = index.min(pages.pages.len() - 1);
        if self.help_doc != doc || self.manual_page(doc) != index {
            if record_history && matches!(self.help_doc, HelpDoc::Gui | HelpDoc::Cli) {
                self.help_manual_back.push((self.help_doc, self.manual_page(self.help_doc)));
                if self.help_manual_back.len() > 128 {
                    self.help_manual_back.remove(0);
                }
                self.help_manual_forward.clear();
            }
            let changed_doc = self.help_doc != doc;
            self.help_doc = doc;
            if doc == HelpDoc::Cli {
                self.help_cli_page = index;
            } else {
                self.help_gui_page = index;
            }
            if changed_doc {
                self.update_help_search_matches();
            }
            self.help_markdown_cache = CommonMarkCache::default();
            self.help_last_content_width = 0.0;
            self.help_scroll_to_top = true;
        }
    }

    fn navigate_manual_history(&mut self, forward: bool) {
        self.ensure_manual_pages(self.help_doc);
        let target = if forward {
            self.help_manual_forward.pop()
        } else {
            self.help_manual_back.pop()
        };
        if let Some((doc, target)) = target {
            let current = (self.help_doc, self.manual_page(self.help_doc));
            if forward {
                self.help_manual_back.push(current);
            } else {
                self.help_manual_forward.push(current);
            }
            self.select_manual_location(doc, target, false);
        }
    }

    fn open_help_fragment(&mut self, fragment: &str) -> Result<(), String> {
        if matches!(self.help_doc, HelpDoc::Gui | HelpDoc::Cli) {
            let pages = self.ensure_manual_pages(self.help_doc);
            let path = pages.pages[self.manual_page(self.help_doc)].document_path.clone();
            self.open_manual_target(self.help_doc, &path, fragment)?;
        } else {
            self.help_link_notice = Some(format!(
                "Page opened. Use Find for section '{fragment}'; automatic section scrolling is not supported."
            ));
        }
        Ok(())
    }

    fn open_gui_manual_target(&mut self, path: &str, fragment: &str) -> Result<(), String> {
        self.open_manual_target(HelpDoc::Gui, path, fragment)
    }

    fn open_manual_target(&mut self, doc: HelpDoc, path: &str, fragment: &str) -> Result<(), String> {
        let pages = self.ensure_manual_pages(doc);
        let (book, root, aliases) = if doc == HelpDoc::Gui {
            (&self.help_gui_book, "docs/gui.md", super::help_manual_sources::GUI_ALIASES)
        } else {
            (&self.help_cli_book, "docs/cli.md", super::help_manual_sources::CLI_ALIASES)
        };
        let has_book = book.as_ref().is_some_and(|book| book.source_hash == pages.source_hash);
        let alias = if path == root && has_book {
            let decoded = super::help_manual::decoded_fragment(fragment);
            aliases.iter().find(|alias| Some(alias.anchor) == decoded.as_deref())
        } else {
            None
        };
        let (path, fragment) = alias.map_or((path, fragment), |alias| (alias.path, alias.fragment));
        let target = pages.page_for_document(path, fragment)
            .ok_or_else(|| format!("Help section '{path}#{fragment}' could not be found."))?;
        self.select_manual_location(doc, target, true);
        self.open_help_doc(doc);
        self.help_scroll_to_top = true;
        self.help_link_notice = None;
        Ok(())
    }

    fn logical_help_document_path(&self) -> Option<String> {
        Some(match self.help_doc {
            HelpDoc::Gui => self.help_gui_pages.as_ref()
                .and_then(|pages| pages.pages.get(self.help_gui_page))
                .map(|page| page.document_path.clone()).unwrap_or_else(|| "docs/gui.md".to_string()),
            HelpDoc::Cli => self.help_cli_pages.as_ref()
                .and_then(|pages| pages.pages.get(self.help_cli_page))
                .map(|page| page.document_path.clone()).unwrap_or_else(|| "docs/cli.md".to_string()),
            HelpDoc::AgentInterface => "docs/agent_interface.md".to_string(),
            HelpDoc::ReviewerPreview => "docs/reviewer_preview.md".to_string(),
            HelpDoc::Shell => "docs/glossary.json".to_string(),
            HelpDoc::Tutorial => self.help_tutorial_entries.get(self.help_tutorial_selected)?.path.clone(),
        })
    }

    pub(super) fn show_manual_search_match(&mut self) {
        if matches!(self.help_doc, HelpDoc::Gui | HelpDoc::Cli)
            && let Some(hit) = self.help_search_matches.get(self.help_search_selected)
        {
            let line = hit.line_number;
            let pages = self.ensure_manual_pages(self.help_doc);
            self.select_manual_location(self.help_doc, pages.page_for_line(line), true);
        }
    }

    pub(super) fn help_document_path(&self) -> Option<PathBuf> {
        if matches!(self.help_doc, HelpDoc::Gui | HelpDoc::Cli) {
            let path = self.logical_help_document_path()?;
            let (book, text) = if self.help_doc == HelpDoc::Gui {
                (&self.help_gui_book, &self.help_gui_markdown)
            } else {
                (&self.help_cli_book, &self.help_cli_markdown)
            };
            if let Some(book) = book
                && book.source_hash == Self::help_markdown_hash(text)
            {
                return book.root.as_ref().map(|root| root.join(path));
            }
            return Self::resolve_runtime_doc_path(&path);
        }
        let path = match self.help_doc {
            HelpDoc::Gui => "docs/gui.md",
            HelpDoc::Cli => "docs/cli.md",
            HelpDoc::AgentInterface => "docs/agent_interface.md",
            HelpDoc::ReviewerPreview => "docs/reviewer_preview.md",
            HelpDoc::Shell => "docs/glossary.json",
            HelpDoc::Tutorial => {
                return self
                    .help_tutorial_entries
                    .get(self.help_tutorial_selected)
                    .map(|entry| PathBuf::from(&entry.path));
            }
        };
        Self::resolve_runtime_doc_path(path)
    }

    fn help_relative_markdown_link(destination: &str) -> Option<&str> {
        if destination.starts_with('#') {
            return Some("");
        }
        let path = destination.split(['#', '?']).next()?;
        if path.is_empty()
            || path.starts_with('/')
            || path.starts_with('\\')
            || Self::is_windows_drive_path(path)
            || Self::has_uri_scheme(path)
            || !Path::new(path)
                .extension()
                .is_some_and(|extension| extension.eq_ignore_ascii_case("md"))
        {
            return None;
        }
        Some(path)
    }

    fn open_help_markdown_link(&mut self, destination: &str) -> Result<(), String> {
        let relative = Self::help_relative_markdown_link(destination)
            .ok_or_else(|| "Not a relative Help Markdown link".to_string())?;
        if relative.is_empty() {
            return self.open_help_fragment(destination.trim_start_matches('#'));
        }
        if let Some(source) = self.logical_help_document_path()
            && let Some(target) = super::help_manual::relative_document_path(&source, relative)
        {
            for (doc, sources) in [(HelpDoc::Gui, super::help_manual_sources::GUI_SOURCES),
                                   (HelpDoc::Cli, super::help_manual_sources::CLI_SOURCES)] {
                if sources.iter().any(|source| source.path == target) {
                    let fragment = destination.split_once('#').map_or("", |(_, fragment)| fragment);
                    return self.open_manual_target(doc, &target, fragment);
                }
            }
            // These manuals have embedded single-document fallbacks as well.
            for (path, doc) in [("docs/agent_interface.md", HelpDoc::AgentInterface),
                                ("docs/reviewer_preview.md", HelpDoc::ReviewerPreview)] {
                if target == path {
                    self.open_help_doc(doc);
                    if let Some((_, fragment)) = destination.split_once('#') {
                        self.open_help_fragment(fragment)?;
                    }
                    return Ok(());
                }
            }
        }
        let source = self
            .help_document_path()
            .ok_or_else(|| "The current Help page has no source path".to_string())?;
        let target = source
            .parent()
            .ok_or_else(|| "The current Help page has no source directory".to_string())?
            .join(relative)
            .canonicalize()
            .map_err(|_| format!("Help link '{destination}' could not be found."))?;
        if !target.is_file() {
            return Err(format!("Help link '{destination}' is not a document."));
        }
        // Tutorial entries can carry absolute installed paths. Resolve those
        // back to a known book/file before applying legacy fragments, rather
        // than interpreting the fragment against the last selected page.
        let book_target = [
            (HelpDoc::Gui, super::help_manual_sources::GUI_SOURCES, &self.help_gui_book),
            (HelpDoc::Cli, super::help_manual_sources::CLI_SOURCES, &self.help_cli_book),
        ].into_iter().find_map(|(doc, sources, book)| {
            sources.iter().find(|source| {
                book.as_ref().and_then(|book| book.root.as_ref())
                    .map(|root| root.join(source.path))
                    .or_else(|| Self::resolve_runtime_doc_path(source.path))
                    .and_then(|path| path.canonicalize().ok()).as_ref() == Some(&target)
            }).map(|source| (doc, source.path))
        });
        if let Some((doc, path)) = book_target {
            let fragment = destination.split_once('#').map_or("", |(_, fragment)| fragment);
            return self.open_manual_target(doc, path, fragment);
        }
        let manual = [
            (HelpDoc::AgentInterface, "docs/agent_interface.md"),
            (HelpDoc::ReviewerPreview, "docs/reviewer_preview.md"),
        ]
        .into_iter()
        .find(|(_, path)| {
            Self::resolve_runtime_doc_path(path)
                .and_then(|path| path.canonicalize().ok())
                .as_ref()
                == Some(&target)
        });
        if let Some((doc, _)) = manual {
            self.open_help_doc(doc);
        } else if let Some(index) = self
            .help_tutorial_entries
            .iter()
            .position(|entry| Path::new(&entry.path).canonicalize().ok().as_ref() == Some(&target))
        {
            self.open_help_tutorial_doc(index);
        } else {
            let path = target
                .to_str()
                .ok_or_else(|| "The Help document path is not valid UTF-8".to_string())?;
            self.open_help_tutorial_path(path, "", "")?;
        }
        self.help_markdown_cache = CommonMarkCache::default();
        self.help_last_content_width = 0.0;
        self.help_link_notice = None;
        if let Some((_, fragment)) = destination.split_once('#')
            && !fragment.is_empty()
        {
            self.open_help_fragment(fragment)?;
        }
        Ok(())
    }

    pub(super) fn follow_help_markdown_links(&mut self, ctx: &egui::Context, command_start: usize) {
        // Only consume relative document links emitted by this Help render.
        // Browser URLs and output from other windows keep their native handling.
        let destinations = ctx.output_mut(|output| {
            let commands = output.commands.split_off(command_start);
            let mut destinations = Vec::new();
            for command in commands {
                if let egui::OutputCommand::OpenUrl(url) = &command
                    && Self::help_relative_markdown_link(&url.url).is_some()
                {
                    destinations.push(url.url.clone());
                } else {
                    output.commands.push(command);
                }
            }
            destinations
        });
        for destination in destinations {
            self.help_link_notice = None;
            if let Err(message) = self.open_help_markdown_link(&destination) {
                self.help_link_notice = Some(message);
            }
            ctx.request_repaint();
        }
    }

    fn help_tutorial_entry_from_catalog_entry(
        entry: crate::workflow_examples::TutorialCatalogEntry,
    ) -> Option<HelpTutorialDocEntry> {
        let resolved_path = Self::resolve_runtime_doc_path(&entry.path)?;
        let title = Self::load_help_tutorial_heading(&resolved_path)
            .ok()
            .flatten()
            .unwrap_or_else(|| {
                if entry.title.trim().is_empty() {
                    Self::markdown_title_from_path(&resolved_path)
                } else {
                    entry.title.clone()
                }
            });
        let mut summary = format!(
            "{}\ntype: {}\nstatus: {}",
            entry.path, entry.entry_type, entry.status
        );
        if let Some(minutes) = entry.estimated_minutes {
            summary.push_str(&format!("\ntime: about {minutes} minutes"));
        }
        if let Some(difficulty) = entry.difficulty.as_deref() {
            summary.push_str("\ndifficulty: ");
            summary.push_str(difficulty);
        }
        if let Some(network) = entry.network.as_deref() {
            summary.push_str("\nnetwork: ");
            summary.push_str(network);
        }
        if !entry.interfaces.is_empty() {
            summary.push_str("\ninterfaces: ");
            summary.push_str(&entry.interfaces.join(", "));
        }
        if let Some(starting_state) = entry.starting_state.as_deref() {
            summary.push_str("\nstart: ");
            summary.push_str(starting_state);
        }
        if !entry.produces.is_empty() {
            summary.push_str("\nproduces: ");
            summary.push_str(&entry.produces.join(", "));
        }
        summary.push('\n');
        summary.push_str(&crate::workflow_examples::tutorial_review_badge_label(
            entry.review_status.as_deref(),
            entry.review_stale,
            entry.codex_reviewed_at.as_deref(),
            entry.human_reviewed_at.as_deref(),
            entry.human_reviewer.as_deref(),
        ));
        if !entry.notes.trim().is_empty() {
            summary.push('\n');
            summary.push_str(entry.notes.trim());
        }
        Some(HelpTutorialDocEntry {
            tutorial_id: entry.id,
            title,
            path: resolved_path.to_string_lossy().to_string(),
            summary,
            audiences: entry.audiences,
            group_label: entry.group_label,
            group_order: entry.group_order,
            group_position: entry.group_position,
            decimal_id: entry.decimal_id,
            prerequisites: entry.prerequisites,
            next_steps: entry.next_steps,
            review_status: entry.review_status,
            codex_reviewed_at: entry.codex_reviewed_at,
            human_reviewed_at: entry.human_reviewed_at,
            human_reviewer: entry.human_reviewer,
            review_stale: entry.review_stale,
        })
    }

    pub(super) fn split_leading_markdown_front_matter(markdown: &str) -> Option<(&str, &str)> {
        let start = if markdown.starts_with("---\n") {
            "---\n".len()
        } else if markdown.starts_with("---\r\n") {
            "---\r\n".len()
        } else {
            return None;
        };
        let mut offset = start;
        while offset < markdown.len() {
            let line_start = offset;
            let line_end = markdown[offset..]
                .find('\n')
                .map(|relative| offset + relative)
                .unwrap_or(markdown.len());
            let line = markdown[line_start..line_end].trim_end_matches('\r');
            let next_line_start = if line_end < markdown.len() {
                line_end + 1
            } else {
                line_end
            };
            if line == "---" {
                return Some((&markdown[start..line_start], &markdown[next_line_start..]));
            }
            offset = next_line_start;
        }
        None
    }

    pub(super) fn markdown_front_matter_value(front_matter: &str, key: &str) -> Option<String> {
        let prefix = format!("{key}:");
        front_matter.lines().find_map(|line| {
            let value = line.trim().strip_prefix(&prefix)?.trim();
            if value.is_empty() {
                return None;
            }
            Some(value.trim_matches('"').to_string())
        })
    }

    pub(super) fn insert_markdown_note_after_first_heading(markdown: &str, note: &str) -> String {
        let trimmed = markdown.trim_start_matches(['\r', '\n']);
        if !trimmed.starts_with('#') {
            return format!("{note}{trimmed}");
        }
        let Some(first_line_end) = trimmed.find('\n') else {
            return format!("{trimmed}\n\n{note}");
        };
        let (heading, body) = trimmed.split_at(first_line_end + 1);
        format!("{heading}\n{note}{}", body.trim_start_matches(['\r', '\n']))
    }

    pub(super) fn summarize_markdown_front_matter_for_help(markdown: &str) -> String {
        let Some((front_matter, body)) = Self::split_leading_markdown_front_matter(markdown) else {
            return markdown.to_string();
        };
        let tier = Self::markdown_front_matter_value(front_matter, "tier")
            .unwrap_or_else(|| "unspecified tier".to_string());
        let automated = Self::markdown_front_matter_value(front_matter, "automated_status")
            .unwrap_or_else(|| "automation not recorded".to_string());
        let review = Self::markdown_front_matter_value(front_matter, "review_status")
            .unwrap_or_else(|| "review not recorded".to_string());
        let note = format!(
            "_Quick facts: `{tier}` · `{automated}` · review `{review}`. Full provenance is at the end._\n\n"
        );
        Self::insert_markdown_note_after_first_heading(body, &note)
    }

    pub(super) fn rewrite_markdown_inline_code_soft_breaks(markdown: &str) -> String {
        let mut out = String::with_capacity(markdown.len() + markdown.len() / 16);
        let mut idx = 0usize;
        let len = markdown.len();
        let mut in_fenced_block: Option<(u8, usize)> = None;
        let mut line_start = true;
        let bytes = markdown.as_bytes();

        while idx < len {
            if line_start {
                let mut look = idx;
                while look < len && bytes[look] == b' ' {
                    look += 1;
                }
                if look < len && (bytes[look] == b'`' || bytes[look] == b'~') {
                    let fence_byte = bytes[look];
                    let mut fence_end = look;
                    while fence_end < len && bytes[fence_end] == fence_byte {
                        fence_end += 1;
                    }
                    let fence_count = fence_end - look;
                    if fence_count >= 3 {
                        match in_fenced_block {
                            Some((current_byte, current_count))
                                if current_byte == fence_byte && fence_count >= current_count =>
                            {
                                in_fenced_block = None;
                            }
                            None => {
                                in_fenced_block = Some((fence_byte, fence_count));
                            }
                            _ => {}
                        }
                    }
                }
            }

            if in_fenced_block.is_none() && bytes[idx] == b'`' {
                let mut tick_end = idx;
                while tick_end < len && bytes[tick_end] == b'`' {
                    tick_end += 1;
                }
                let tick_count = tick_end - idx;
                if let Some(close_start) =
                    Self::find_matching_inline_backtick_run(markdown, tick_end, tick_count)
                {
                    out.push_str(&markdown[idx..tick_end]);
                    out.push_str(&Self::insert_soft_breaks_into_inline_code(
                        &markdown[tick_end..close_start],
                    ));
                    out.push_str(&markdown[close_start..close_start + tick_count]);
                    idx = close_start + tick_count;
                    line_start = false;
                    continue;
                }
            }

            let ch = markdown[idx..].chars().next().unwrap();
            out.push(ch);
            line_start = ch == '\n';
            idx += ch.len_utf8();
        }

        out
    }

    pub(super) fn find_matching_inline_backtick_run(
        markdown: &str,
        start: usize,
        tick_count: usize,
    ) -> Option<usize> {
        let bytes = markdown.as_bytes();
        let len = markdown.len();
        let mut idx = start;
        while idx < len {
            if bytes[idx] == b'\n' {
                return None;
            }
            if bytes[idx] == b'`' {
                let mut end = idx;
                while end < len && bytes[end] == b'`' {
                    end += 1;
                }
                if end - idx == tick_count {
                    return Some(idx);
                }
                idx = end;
                continue;
            }
            idx += markdown[idx..].chars().next().unwrap().len_utf8();
        }
        None
    }

    pub(super) fn insert_soft_breaks_into_inline_code(text: &str) -> String {
        const ZWSP: char = '\u{200B}';
        let mut out = String::with_capacity(text.len() + text.len() / 8);
        for ch in text.chars() {
            out.push(ch);
            if matches!(
                ch,
                ' ' | '/'
                    | '_'
                    | '-'
                    | '.'
                    | ':'
                    | ','
                    | ';'
                    | '('
                    | ')'
                    | '['
                    | ']'
                    | '{'
                    | '}'
                    | '='
                    | '+'
                    | '|'
                    | '>'
            ) {
                out.push(ZWSP);
            }
        }
        out
    }

    pub(super) fn rewrite_inline_image_destination(
        span: &str,
        absolute_dest: &Path,
    ) -> Option<String> {
        let (dest_start, dest_end) = Self::find_inline_image_destination(span)?;
        let mut rewritten = String::with_capacity(span.len() + 64);
        let absolute_uri = Self::file_uri_from_path(absolute_dest);
        rewritten.push_str(&span[..dest_start]);
        rewritten.push_str(&absolute_uri);
        rewritten.push_str(&span[dest_end..]);
        Some(rewritten)
    }

    pub(super) fn help_svg_png_cache_path(svg_path: &Path) -> Option<PathBuf> {
        if !svg_path
            .extension()
            .and_then(|ext| ext.to_str())
            .map(|ext| ext.eq_ignore_ascii_case("svg"))
            .unwrap_or(false)
        {
            return None;
        }
        let metadata = fs::metadata(svg_path).ok();
        let mut fingerprint_input = svg_path.to_string_lossy().to_string();
        if let Some(metadata) = metadata {
            fingerprint_input.push_str(&format!("\nlen={}", metadata.len()));
            if let Ok(modified) = metadata.modified()
                && let Ok(duration) = modified.duration_since(UNIX_EPOCH)
            {
                fingerprint_input.push_str(&format!(
                    "\nmtime={}.{}",
                    duration.as_secs(),
                    duration.subsec_nanos()
                ));
            }
        }
        let digest = crate::digest_utils::sha256_hex_str(&fingerprint_input);
        let stem = svg_path
            .file_stem()
            .and_then(|stem| stem.to_str())
            .unwrap_or("image")
            .chars()
            .map(|ch| {
                if ch.is_ascii_alphanumeric() || matches!(ch, '-' | '_') {
                    ch
                } else {
                    '_'
                }
            })
            .collect::<String>();
        Some(
            env::temp_dir()
                .join("gentle_help_svg_png")
                .join(format!("{stem}.{digest}.png")),
        )
    }

    pub(super) fn help_image_render_path(absolute_dest: &Path) -> PathBuf {
        use crate::gui_profiler::startup_trace::{
            HelpImageObservation, measure_help_rasterization, observe_help_image,
        };
        let Some(png_path) = Self::help_svg_png_cache_path(absolute_dest) else {
            return absolute_dest.to_path_buf();
        };
        observe_help_image(HelpImageObservation::SvgReference);
        if png_path.is_file() {
            observe_help_image(HelpImageObservation::CacheHit);
            return png_path;
        }
        if let Some(parent) = png_path.parent()
            && fs::create_dir_all(parent).is_err()
        {
            observe_help_image(HelpImageObservation::PreparationFailed);
            return absolute_dest.to_path_buf();
        }
        match measure_help_rasterization(|| {
            crate::svg_png::render_svg_file_to_png(
                absolute_dest,
                &png_path,
                crate::svg_png::SvgPngRenderOptions {
                    scale: 1.0,
                    drop_dotplot_metadata: false,
                },
            )
        }) {
            Ok(_) => png_path,
            Err(_) => absolute_dest.to_path_buf(),
        }
    }

    pub(super) fn find_inline_image_destination(span: &str) -> Option<(usize, usize)> {
        let bytes = span.as_bytes();
        if bytes.len() < 4 || bytes[0] != b'!' || bytes[1] != b'[' {
            return None;
        }

        let mut i = 2usize;
        let mut bracket_depth = 1usize;
        while i < bytes.len() {
            match bytes[i] {
                b'\\' => {
                    i = (i + 2).min(bytes.len());
                }
                b'[' => {
                    bracket_depth += 1;
                    i += 1;
                }
                b']' => {
                    bracket_depth = bracket_depth.saturating_sub(1);
                    i += 1;
                    if bracket_depth == 0 {
                        break;
                    }
                }
                _ => i += 1,
            }
        }
        if bracket_depth != 0 {
            return None;
        }

        while i < bytes.len() && bytes[i].is_ascii_whitespace() {
            i += 1;
        }
        if i >= bytes.len() || bytes[i] != b'(' {
            return None;
        }
        i += 1;

        while i < bytes.len() && bytes[i].is_ascii_whitespace() {
            i += 1;
        }
        if i >= bytes.len() {
            return None;
        }

        if bytes[i] == b'<' {
            let dest_start = i + 1;
            i += 1;
            while i < bytes.len() {
                match bytes[i] {
                    b'\\' => {
                        i = (i + 2).min(bytes.len());
                    }
                    b'>' => return Some((dest_start, i)),
                    _ => i += 1,
                }
            }
            return None;
        }

        let dest_start = i;
        let mut paren_depth = 0usize;
        while i < bytes.len() {
            match bytes[i] {
                b'\\' => {
                    i = (i + 2).min(bytes.len());
                }
                b'(' => {
                    paren_depth += 1;
                    i += 1;
                }
                b')' => {
                    if paren_depth == 0 {
                        break;
                    }
                    paren_depth -= 1;
                    i += 1;
                }
                b if b.is_ascii_whitespace() && paren_depth == 0 => break,
                _ => i += 1,
            }
        }
        if dest_start == i {
            return None;
        }
        Some((dest_start, i))
    }

    pub(super) fn resolve_relative_image_path(path: &str, base_dir: &Path) -> Option<PathBuf> {
        if path.is_empty()
            || path.starts_with('/')
            || path.starts_with('\\')
            || path.starts_with('#')
            || Self::is_windows_drive_path(path)
            || Self::has_uri_scheme(path)
        {
            return None;
        }
        let joined = base_dir.join(path);
        Some(joined.canonicalize().unwrap_or(joined))
    }

    pub(super) fn is_windows_drive_path(path: &str) -> bool {
        let bytes = path.as_bytes();
        bytes.len() >= 3
            && bytes[0].is_ascii_alphabetic()
            && bytes[1] == b':'
            && (bytes[2] == b'/' || bytes[2] == b'\\')
    }

    pub(super) fn has_uri_scheme(path: &str) -> bool {
        if Self::is_windows_drive_path(path) {
            return false;
        }
        let mut chars = path.chars();
        let Some(first) = chars.next() else {
            return false;
        };
        if !first.is_ascii_alphabetic() {
            return false;
        }
        for ch in chars {
            if ch == ':' {
                return true;
            }
            if !(ch.is_ascii_alphanumeric() || matches!(ch, '+' | '-' | '.')) {
                return false;
            }
        }
        false
    }

    pub(super) fn encode_hex_nibble(value: u8) -> char {
        match value {
            0..=9 => (b'0' + value) as char,
            _ => (b'A' + (value - 10)) as char,
        }
    }

    pub(super) fn percent_encode_file_uri_path(path: &str) -> String {
        let mut encoded = String::with_capacity(path.len());
        for byte in path.bytes() {
            let keep = byte.is_ascii_alphanumeric()
                || matches!(byte, b'-' | b'_' | b'.' | b'~' | b'/' | b':');
            if keep {
                encoded.push(byte as char);
                continue;
            }
            encoded.push('%');
            encoded.push(Self::encode_hex_nibble(byte >> 4));
            encoded.push(Self::encode_hex_nibble(byte & 0x0f));
        }
        encoded
    }

    pub(super) fn file_uri_from_path(path: &Path) -> String {
        let normalized = path.to_string_lossy().replace('\\', "/");
        let rooted = if normalized.starts_with('/') {
            normalized
        } else {
            format!("/{normalized}")
        };
        let encoded = Self::percent_encode_file_uri_path(&rooted);
        format!("file://{encoded}")
    }

    pub(super) fn load_help_doc(path: &str, fallback: &'static str) -> String {
        if let Some(runtime_path) = Self::resolve_runtime_doc_path(path)
            && let Ok(text) = fs::read_to_string(&runtime_path)
        {
            if let Some(base_dir) = runtime_path.parent() {
                return Self::rewrite_markdown_relative_image_links(&text, base_dir);
            }
            return text;
        }
        fallback.to_string()
    }

    pub(super) fn load_help_markdown_from_path(path: &Path) -> Option<String> {
        let text = fs::read_to_string(path).ok()?;
        if let Some(base_dir) = path.parent() {
            return Some(Self::rewrite_markdown_relative_image_links(&text, base_dir));
        }
        Some(text)
    }

    fn load_help_tutorial_heading(path: &Path) -> std::io::Result<Option<String>> {
        // Discovery preserves the document's heading without preparing images
        // (and their fonts). Only the selected page uses the rendering loader.
        let markdown = fs::read_to_string(path)?;
        Ok(Self::markdown_first_heading(&markdown))
    }

    pub(super) fn collect_markdown_paths_recursive(root: &Path, out: &mut Vec<PathBuf>) {
        let Ok(entries) = fs::read_dir(root) else {
            return;
        };
        for entry in entries.flatten() {
            let path = entry.path();
            if path.is_dir() {
                Self::collect_markdown_paths_recursive(&path, out);
                continue;
            }
            let Some(ext) = path.extension().and_then(|value| value.to_str()) else {
                continue;
            };
            if ext.eq_ignore_ascii_case("md") {
                out.push(path);
            }
        }
    }

    pub(super) fn markdown_first_heading(markdown: &str) -> Option<String> {
        markdown.lines().find_map(|line| {
            let trimmed = line.trim();
            let heading = trimmed.strip_prefix('#')?;
            let title = heading.trim_start_matches('#').trim();
            if title.is_empty() {
                return None;
            }
            Some(title.to_string())
        })
    }

    pub(super) fn markdown_title_from_path(path: &Path) -> String {
        let stem = path
            .file_stem()
            .and_then(|value| value.to_str())
            .unwrap_or("tutorial");
        let cleaned =
            stem.trim_start_matches(|ch: char| ch.is_ascii_digit() || ch == '_' || ch == '-');
        let normalized = cleaned
            .replace(['_', '-'], " ")
            .split_whitespace()
            .collect::<Vec<_>>()
            .join(" ");
        if normalized.is_empty() {
            "Tutorial".to_string()
        } else {
            normalized
        }
    }

    pub(super) fn discover_help_tutorial_entries() -> Vec<HelpTutorialDocEntry> {
        if let Some(catalog_path) = Self::resolve_runtime_doc_path(DEFAULT_TUTORIAL_CATALOG_PATH)
            && let Ok(catalog) = load_tutorial_catalog(&catalog_path)
        {
            let mut catalog_entries = catalog
                .entries
                .into_iter()
                .filter_map(Self::help_tutorial_entry_from_catalog_entry)
                .collect::<Vec<_>>();
            Self::ensure_agent_interfaces_tutorial_entry(&mut catalog_entries);
            Self::sort_help_tutorial_entries_by_audience_group(&mut catalog_entries);
            if !catalog_entries.is_empty() {
                return catalog_entries;
            }
        }
        let Some(tutorial_root) = Self::resolve_runtime_doc_path("docs/tutorial") else {
            return vec![];
        };
        let mut markdown_paths = Vec::new();
        Self::collect_markdown_paths_recursive(&tutorial_root, &mut markdown_paths);
        markdown_paths.sort_by(|left, right| {
            let left_rel = left.strip_prefix(&tutorial_root).ok();
            let right_rel = right.strip_prefix(&tutorial_root).ok();
            let left_depth = left_rel
                .map(|value| value.components().count())
                .unwrap_or(usize::MAX);
            let right_depth = right_rel
                .map(|value| value.components().count())
                .unwrap_or(usize::MAX);
            left_depth.cmp(&right_depth).then_with(|| left.cmp(right))
        });

        let mut entries = markdown_paths
            .into_iter()
            .map(|path| {
                let title = Self::load_help_tutorial_heading(&path)
                    .ok()
                    .flatten()
                    .unwrap_or_else(|| Self::markdown_title_from_path(&path));
                let relative = path
                    .strip_prefix(&tutorial_root)
                    .map(|value| value.to_string_lossy().replace('\\', "/"))
                    .unwrap_or_else(|_| path.to_string_lossy().to_string());
                HelpTutorialDocEntry {
                    tutorial_id: String::new(),
                    title,
                    path: path.to_string_lossy().to_string(),
                    summary: format!("docs/tutorial/{relative}"),
                    audiences: vec![],
                    group_label: None,
                    group_order: None,
                    group_position: None,
                    decimal_id: None,
                    prerequisites: vec![],
                    next_steps: vec![],
                    review_status: None,
                    codex_reviewed_at: None,
                    human_reviewed_at: None,
                    human_reviewer: None,
                    review_stale: false,
                }
            })
            .collect::<Vec<_>>();
        Self::ensure_agent_interfaces_tutorial_entry(&mut entries);
        Self::sort_help_tutorial_entries_by_audience_group(&mut entries);
        entries
    }

    pub(super) fn discover_guided_walkthrough_entries() -> Vec<HelpTutorialDocEntry> {
        if let Some(catalog_path) = Self::resolve_runtime_doc_path(DEFAULT_TUTORIAL_CATALOG_PATH)
            && let Ok(catalog) = load_tutorial_catalog(&catalog_path)
        {
            let mut entries = catalog
                .entries
                .into_iter()
                .filter(|entry| entry.status == "manual/reference")
                .filter_map(Self::help_tutorial_entry_from_catalog_entry)
                .collect::<Vec<_>>();
            Self::ensure_agent_interfaces_tutorial_entry(&mut entries);
            Self::sort_help_tutorial_entries_by_audience_group(&mut entries);
            if !entries.is_empty() {
                return entries;
            }
        }
        let mut entries = vec![];
        Self::ensure_agent_interfaces_tutorial_entry(&mut entries);
        Self::sort_help_tutorial_entries_by_audience_group(&mut entries);
        entries
    }

    pub(super) fn ensure_agent_interfaces_tutorial_entry(entries: &mut Vec<HelpTutorialDocEntry>) {
        let Some(resolved_path) = Self::resolve_runtime_doc_path(AGENT_INTERFACES_TUTORIAL_PATH)
        else {
            return;
        };
        let resolved_string = resolved_path.to_string_lossy().to_string();
        if entries.iter().any(|entry| entry.path == resolved_string) {
            return;
        }
        let title = Self::load_help_tutorial_heading(&resolved_path)
            .ok()
            .flatten()
            .unwrap_or_else(|| AGENT_INTERFACES_TUTORIAL_TITLE.to_string());
        entries.push(HelpTutorialDocEntry {
            tutorial_id: "agent_interfaces".to_string(),
            title,
            path: resolved_string,
            summary: AGENT_INTERFACES_TUTORIAL_SUMMARY.to_string(),
            audiences: vec![
                "orientation_interfaces".to_string(),
                "agent_users".to_string(),
                "mcp_users".to_string(),
            ],
            group_label: Some("Getting Started & Interfaces".to_string()),
            group_order: Some(1),
            group_position: Some(1),
            decimal_id: Some("01.01".to_string()),
            prerequisites: vec![],
            next_steps: vec![],
            review_status: Some("unreviewed".to_string()),
            codex_reviewed_at: None,
            human_reviewed_at: None,
            human_reviewer: None,
            review_stale: false,
        });
    }

    pub(super) fn set_help_tutorial_selected(&mut self, tutorial_index: usize) -> bool {
        if self.help_tutorial_entries.is_empty() {
            self.help_tutorial_selected = 0;
            self.help_tutorial_title = "Tutorial".to_string();
            self.help_tutorial_markdown =
                "# Tutorials\n\nNo tutorial entries were found from `docs/tutorial/catalog.json` or the `docs/tutorial` markdown tree."
                    .to_string();
            return false;
        }
        let clamped_index = tutorial_index.min(self.help_tutorial_entries.len() - 1);
        let changed = clamped_index != self.help_tutorial_selected;
        let trace = (self.show_help_dialog && changed).then(|| {
            crate::gui_profiler::startup_trace::span_once(
                crate::gui_profiler::startup_trace::Phase::HelpTutorialSwitch,
            )
        });
        self.help_tutorial_selected = clamped_index;
        if let Some(entry) = self.help_tutorial_entries.get(clamped_index) {
            self.help_tutorial_title = entry.title.clone();
            let tutorial_path = PathBuf::from(&entry.path);
            self.help_tutorial_markdown = Self::load_help_markdown_from_path(&tutorial_path)
                .unwrap_or_else(|| {
                    format!(
                        "# {}\n\nCould not load tutorial markdown from `{}`.",
                        entry.title, entry.summary
                    )
                });
        }
        if let Some(trace) = trace {
            trace.finish(true);
        }
        changed
    }

    pub(super) fn help_tutorial_navigation_targets(
        &self,
        tutorial_ids: &[String],
    ) -> Vec<(usize, String, String)> {
        tutorial_ids
            .iter()
            .filter_map(|tutorial_id| {
                self.help_tutorial_entries
                    .iter()
                    .enumerate()
                    .find(|(_, entry)| entry.tutorial_id == *tutorial_id)
                    .map(|(index, entry)| {
                        (
                            index,
                            entry
                                .decimal_id
                                .clone()
                                .unwrap_or_else(|| entry.title.clone()),
                            format!("{}\n{}", entry.title, entry.summary),
                        )
                    })
            })
            .collect()
    }

    pub(super) fn open_help_tutorial_doc(&mut self, tutorial_index: usize) {
        self.help_link_notice = None;
        let trace = crate::gui_profiler::startup_trace::span_once(
            crate::gui_profiler::startup_trace::Phase::HelpTutorialOpen,
        );
        if self.help_tutorial_entries.is_empty() {
            self.help_tutorial_entries = Self::discover_help_tutorial_entries();
        }
        let opening_from_closed_state = !self.show_help_dialog;
        if !self.show_help_dialog {
            self.mark_viewport_open_requested(Self::help_viewport_id());
        }
        let previous_doc = self.help_doc;
        let previous_tutorial_index = self.help_tutorial_selected;
        self.set_help_tutorial_selected(tutorial_index);
        self.help_doc = HelpDoc::Tutorial;
        if opening_from_closed_state
            || previous_doc != HelpDoc::Tutorial
            || previous_tutorial_index != self.help_tutorial_selected
        {
            self.refresh_help_search_matches();
        }
        self.help_focus_search_box = true;
        self.show_help_dialog = true;
        self.queue_focus_viewport(Self::help_viewport_id());
        trace.finish(true);
    }

    pub(super) fn open_help_tutorial_path(
        &mut self,
        tutorial_path: &str,
        fallback_title: &str,
        fallback_summary: &str,
    ) -> std::result::Result<(), String> {
        let resolved_path = if Path::new(tutorial_path).is_absolute() {
            PathBuf::from(tutorial_path)
        } else {
            Self::resolve_runtime_doc_path(tutorial_path).ok_or_else(|| {
                format!(
                    "Tutorial guide '{}' could not be resolved at runtime",
                    tutorial_path
                )
            })?
        };
        let heading = Self::load_help_tutorial_heading(&resolved_path).map_err(|_| {
            format!(
                "Tutorial guide '{}' could not be loaded",
                resolved_path.display()
            )
        })?;
        if self.help_tutorial_entries.is_empty() {
            self.help_tutorial_entries = Self::discover_help_tutorial_entries();
        }
        let resolved_string = resolved_path.to_string_lossy().to_string();
        let tutorial_index = match self
            .help_tutorial_entries
            .iter()
            .position(|entry| entry.path == resolved_string)
        {
            Some(index) => index,
            None => {
                let title = heading.unwrap_or_else(|| {
                    let trimmed = fallback_title.trim();
                    if trimmed.is_empty() {
                        Self::markdown_title_from_path(&resolved_path)
                    } else {
                        trimmed.to_string()
                    }
                });
                let summary = if fallback_summary.trim().is_empty() {
                    tutorial_path.to_string()
                } else {
                    fallback_summary.trim().to_string()
                };
                self.help_tutorial_entries.push(HelpTutorialDocEntry {
                    tutorial_id: String::new(),
                    title,
                    path: resolved_string,
                    summary,
                    audiences: vec![],
                    group_label: None,
                    group_order: None,
                    group_position: None,
                    decimal_id: None,
                    prerequisites: vec![],
                    next_steps: vec![],
                    review_status: None,
                    codex_reviewed_at: None,
                    human_reviewed_at: None,
                    human_reviewer: None,
                    review_stale: false,
                });
                self.help_tutorial_entries.len() - 1
            }
        };
        self.open_help_tutorial_doc(tutorial_index);
        Ok(())
    }
}
