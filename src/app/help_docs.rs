//! Help/tutorial markdown loading and rewriting helpers for `GENtleApp`.
//!
//! This submodule keeps the top-level app file focused on window/session state
//! while the help-system path owns its markdown parsing, tutorial discovery,
//! and relative-image rewriting in one place.
//!
//! Look here for:
//! - markdown cleanup/rewrite helpers used before rendering docs in-app
//! - tutorial/help discovery and path resolution
//! - relative asset/image rewriting for embedded help pages

use super::*;

pub(super) const AGENT_INTERFACES_TUTORIAL_PATH: &str = "docs/agent_interfaces_tutorial.md";
pub(super) const AGENT_INTERFACES_TUTORIAL_TITLE: &str =
    "GENtle Agent Assistant and Agent Interfaces Tutorial";
pub(super) const AGENT_INTERFACES_TUTORIAL_SUMMARY: &str = "docs/agent_interfaces_tutorial.md\n\
type: operational_reference\n\
status: manual/reference\n\
Practical guide for the in-app Agent Assistant, provider quick starts, reviewed shared-shell suggestions, CLI/shared shell, MCP, and external coding agents.";

impl GENtleApp {
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
        let chapter_id = Self::markdown_front_matter_value(front_matter, "chapter_id");
        let source_example = Self::markdown_front_matter_value(front_matter, "source_example");
        let mut note = String::from("_Provenance note: ");
        match (chapter_id, source_example) {
            (Some(chapter_id), Some(source_example)) => {
                note.push_str("this generated tutorial is tracked as chapter `");
                note.push_str(&chapter_id);
                note.push_str("` from workflow `");
                note.push_str(&source_example);
                note.push_str("`.");
            }
            (Some(chapter_id), None) => {
                note.push_str("this generated tutorial is tracked as chapter `");
                note.push_str(&chapter_id);
                note.push_str("`.");
            }
            (None, Some(source_example)) => {
                note.push_str("this generated tutorial is tracked from workflow `");
                note.push_str(&source_example);
                note.push_str("`.");
            }
            (None, None) => {
                note.push_str("this generated tutorial has machine-readable source metadata.");
            }
        }
        note.push_str(" The note is only there to preserve tutorial provenance; the hands-on walkthrough starts here, and full canonical source details are repeated at the end._\n\n");
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
        let mut hasher = Sha1::new();
        hasher.update(svg_path.to_string_lossy().as_bytes());
        if let Some(metadata) = metadata {
            hasher.update(metadata.len().to_le_bytes());
            if let Ok(modified) = metadata.modified() {
                if let Ok(duration) = modified.duration_since(UNIX_EPOCH) {
                    hasher.update(duration.as_secs().to_le_bytes());
                    hasher.update(duration.subsec_nanos().to_le_bytes());
                }
            }
        }
        let digest = format!("{:x}", hasher.finalize());
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
        let Some(png_path) = Self::help_svg_png_cache_path(absolute_dest) else {
            return absolute_dest.to_path_buf();
        };
        if png_path.is_file() {
            return png_path;
        }
        if let Some(parent) = png_path.parent() {
            if fs::create_dir_all(parent).is_err() {
                return absolute_dest.to_path_buf();
            }
        }
        match crate::svg_png::render_svg_file_to_png(
            absolute_dest,
            &png_path,
            crate::svg_png::SvgPngRenderOptions {
                scale: 1.0,
                drop_dotplot_metadata: false,
            },
        ) {
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
        if let Some(runtime_path) = Self::resolve_runtime_doc_path(path) {
            if let Ok(text) = fs::read_to_string(&runtime_path) {
                if let Some(base_dir) = runtime_path.parent() {
                    return Self::rewrite_markdown_relative_image_links(&text, base_dir);
                }
                return text;
            }
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
        if let Some(catalog_path) = Self::resolve_runtime_doc_path(DEFAULT_TUTORIAL_CATALOG_PATH) {
            if let Ok(catalog) = load_tutorial_catalog(&catalog_path) {
                let mut catalog_entries = catalog
                    .entries
                    .into_iter()
                    .filter_map(|entry| {
                        let resolved_path = Self::resolve_runtime_doc_path(&entry.path)?;
                        let markdown =
                            Self::load_help_markdown_from_path(&resolved_path).unwrap_or_default();
                        let title = Self::markdown_first_heading(&markdown).unwrap_or_else(|| {
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
                        if !entry.notes.trim().is_empty() {
                            summary.push('\n');
                            summary.push_str(entry.notes.trim());
                        }
                        Some(HelpTutorialDocEntry {
                            title,
                            path: resolved_path.to_string_lossy().to_string(),
                            summary,
                        })
                    })
                    .collect::<Vec<_>>();
                Self::ensure_agent_interfaces_tutorial_entry(&mut catalog_entries);
                if !catalog_entries.is_empty() {
                    return catalog_entries;
                }
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
                let markdown = Self::load_help_markdown_from_path(&path).unwrap_or_default();
                let title = Self::markdown_first_heading(&markdown)
                    .unwrap_or_else(|| Self::markdown_title_from_path(&path));
                let relative = path
                    .strip_prefix(&tutorial_root)
                    .map(|value| value.to_string_lossy().replace('\\', "/"))
                    .unwrap_or_else(|_| path.to_string_lossy().to_string());
                HelpTutorialDocEntry {
                    title,
                    path: path.to_string_lossy().to_string(),
                    summary: format!("docs/tutorial/{relative}"),
                }
            })
            .collect::<Vec<_>>();
        Self::ensure_agent_interfaces_tutorial_entry(&mut entries);
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
        let markdown = Self::load_help_markdown_from_path(&resolved_path).unwrap_or_default();
        let title = Self::markdown_first_heading(&markdown)
            .unwrap_or_else(|| AGENT_INTERFACES_TUTORIAL_TITLE.to_string());
        entries.push(HelpTutorialDocEntry {
            title,
            path: resolved_string,
            summary: AGENT_INTERFACES_TUTORIAL_SUMMARY.to_string(),
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
        changed
    }

    pub(super) fn open_help_tutorial_doc(&mut self, tutorial_index: usize) {
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
        let markdown = Self::load_help_markdown_from_path(&resolved_path).ok_or_else(|| {
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
                let title = Self::markdown_first_heading(&markdown).unwrap_or_else(|| {
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
                    title,
                    path: resolved_string,
                    summary,
                });
                self.help_tutorial_entries.len() - 1
            }
        };
        self.open_help_tutorial_doc(tutorial_index);
        Ok(())
    }
}
