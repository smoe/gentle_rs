//! Presentation-only manual books and pagination at complete Markdown boundaries.
//! Topic files remain authoritative; the joined text is a transient search/copy
//! projection, never a second maintained manual. Each page retains its source.

use pulldown_cmark::{Event, Options, Parser, Tag, TagEnd};
use std::collections::{HashMap, HashSet};
use std::hash::{Hash, Hasher};
use std::ops::Range;
use std::path::{Component, Path, PathBuf};

pub(super) const PAGE_TARGET_BYTES: usize = 12_000;

pub(super) struct ManualSource {
    pub id: &'static str,
    pub path: &'static str,
    pub title: &'static str,
    pub group: &'static str,
    pub markdown: &'static str,
}

pub(super) struct ManualAlias {
    pub anchor: &'static str,
    pub path: &'static str,
    pub fragment: &'static str,
}

pub(super) struct ManualDocument {
    pub path: String,
    pub title: String,
    pub group: String,
    pub source_range: Range<usize>,
}

pub(super) struct ManualBook {
    pub markdown: String,
    pub source_hash: u64,
    pub documents: Vec<ManualDocument>,
    pub root: Option<PathBuf>,
}

impl ManualBook {
    pub fn embedded(sources: &[ManualSource]) -> Self {
        Self::from_texts(sources, sources.iter().map(|source| source.markdown.to_owned()).collect(), None)
    }

    pub fn from_texts(sources: &[ManualSource], texts: Vec<String>, root: Option<PathBuf>) -> Self {
        assert_eq!(sources.len(), texts.len(), "one text per manual source");
        let mut markdown = String::new();
        let mut documents = Vec::new();
        for (source, text) in sources.iter().zip(texts) {
            // Compatibility stubs are for old web bookmarks. The native reader
            // routes those aliases directly, so they must not crowd its tree/search.
            let text = if source.id == "index" {
                text.split("<!-- manual-legacy:start -->").next().unwrap_or(&text)
            } else {
                &text
            };
            let start = markdown.len();
            markdown.push_str(text);
            markdown.push_str("\n\n");
            documents.push(ManualDocument {
                path: source.path.to_string(),
                title: source.title.to_string(),
                group: source.group.to_string(),
                source_range: start..markdown.len(),
            });
        }
        let mut hasher = std::collections::hash_map::DefaultHasher::new();
        markdown.hash(&mut hasher);
        Self { markdown, source_hash: hasher.finish(), documents, root }
    }
}

/// Resolve a known manual's logical URL without requiring files on disk.
/// Filesystem/tutorial navigation remains a separate, existing reader path.
pub(super) fn relative_document_path(source: &str, destination: &str) -> Option<String> {
    let decoded = decoded_fragment(destination)?;
    let target = Path::new(source).parent()?.join(decoded);
    let mut parts = Vec::new();
    for component in target.components() {
        match component {
            Component::Normal(part) => parts.push(part.to_str()?),
            Component::ParentDir => { parts.pop()?; }
            Component::CurDir => {}
            _ => return None,
        }
    }
    Some(parts.join("/"))
}

pub(super) struct ManualPage {
    pub title: String,
    pub heading_level: u8,
    pub anchor: Option<String>,
    pub document_path: String,
    pub group: String,
    pub source_range: Range<usize>,
    pub first_line: usize,
    pub markdown: String,
}

pub(super) struct ManualPages {
    pub source_hash: u64,
    pub pages: Vec<ManualPage>,
    anchors: HashMap<String, usize>,
}

/// GitHub-style heading fragments; punctuation is removed, spaces become dashes.
/// Duplicate headings are disambiguated document-wide by `ManualPages::new`.
fn heading_slug(title: &str) -> String {
    title
        .to_lowercase()
        .chars()
        .filter_map(|ch| {
            if ch.is_alphanumeric() || ch == '-' || ch == '_' {
                Some(ch)
            } else if ch.is_whitespace() {
                Some('-')
            } else {
                None
            }
        })
        .collect()
}

pub(super) fn decoded_fragment(fragment: &str) -> Option<String> {
    let mut bytes = Vec::with_capacity(fragment.len());
    let mut input = fragment.as_bytes().iter().copied();
    while let Some(byte) = input.next() {
        if byte == b'%' {
            let high = char::from(input.next()?).to_digit(16)?;
            let low = char::from(input.next()?).to_digit(16)?;
            bytes.push((high * 16 + low) as u8);
        } else {
            bytes.push(byte);
        }
    }
    String::from_utf8(bytes).ok()
}

impl ManualPages {
    pub fn new(source: &str, source_hash: u64) -> Self {
        let options = Options::ENABLE_TABLES
            | Options::ENABLE_TASKLISTS
            | Options::ENABLE_STRIKETHROUGH
            | Options::ENABLE_FOOTNOTES
            | Options::ENABLE_DEFINITION_LIST;
        let parser = Parser::new_ext(source, options);
        let mut definitions: Vec<_> = parser
            .reference_definitions()
            .iter()
            .map(|(_, definition)| definition.span.clone())
            .collect();
        definitions.sort_by_key(|range| range.start);
        let mut boundaries = Vec::new();
        let mut heading_levels = HashMap::new();
        let mut depth = 0usize;
        let mut root_bullet_list = false;
        let mut heading = None;
        let mut heading_text = String::new();
        let mut footnotes = Vec::new();
        let mut footnote_references = Vec::new();
        for (event, range) in parser.into_offset_iter() {
            match event {
                Event::Start(tag) => {
                    // Bullet items are safe boundaries, but keep numbered lists
                    // intact: restarting an implicit "1." list changes numbering.
                    if depth == 0 || (depth == 1 && root_bullet_list && matches!(tag, Tag::Item)) {
                        boundaries.push((range.start, None));
                        if let Tag::Heading { level, .. } = tag {
                            heading = Some(boundaries.len() - 1);
                            heading_levels.insert(range.start, level as u8);
                            heading_text.clear();
                        }
                    }
                    if depth == 0 {
                        root_bullet_list = matches!(tag, Tag::List(None));
                        if let Tag::FootnoteDefinition(name) = &tag {
                            footnotes.push((name.to_string(), range.clone()));
                        }
                    }
                    depth += 1;
                }
                Event::End(tag) => {
                    depth = depth.saturating_sub(1);
                    if matches!(tag, TagEnd::Heading(_))
                        && let Some(index) = heading.take()
                    {
                        boundaries[index].1 = Some(heading_text.trim().to_string());
                    }
                }
                Event::Text(text) | Event::Code(text) if heading.is_some() => {
                    heading_text.push_str(&text);
                }
                Event::SoftBreak | Event::HardBreak if heading.is_some() => {
                    heading_text.push(' ');
                }
                Event::Rule | Event::Html(_) if depth == 0 => {
                    boundaries.push((range.start, None));
                }
                Event::FootnoteReference(name) => {
                    footnote_references.push((name.to_string(), range));
                }
                _ => {}
            }
        }
        let mut ranges = Vec::new();
        let mut start = 0;
        let mut title = "Overview".to_string();
        let mut continuation = 1;
        let mut level = 1;
        let mut anchor = None;
        let mut used_anchors = HashSet::new();
        for (offset, next_heading) in boundaries {
            let next_anchor = next_heading.as_ref().map(|heading| {
                let base = heading_slug(heading);
                let mut candidate = base.clone();
                let mut suffix = 0;
                while !used_anchors.insert(candidate.clone()) {
                    suffix += 1;
                    candidate = format!("{base}-{suffix}");
                }
                candidate
            });
            if offset <= start {
                if let Some(next_heading) = next_heading {
                    title = next_heading;
                    level = heading_levels[&offset];
                    anchor = next_anchor;
                }
                continue;
            }
            if next_heading.is_some() || offset - start >= PAGE_TARGET_BYTES {
                ranges.push((start..offset, title.clone(), continuation, level, anchor.take()));
                start = offset;
                if let Some(next_heading) = next_heading {
                    title = next_heading;
                    continuation = 1;
                    level = heading_levels[&offset];
                    anchor = next_anchor;
                } else {
                    continuation += 1;
                }
            }
        }
        ranges.push((start..source.len(), title, continuation, level, anchor));
        let pages: Vec<_> = ranges
            .into_iter()
            .map(|(range, title, continuation, level, anchor)| {
                let body = &source[range.clone()];
                let mut markdown = body.to_string();
                for definition in &definitions {
                    if definition.start < range.start || definition.end > range.end {
                        markdown.push_str("\n\n");
                        markdown.push_str(&source[definition.clone()]);
                    }
                }
                // Preserve cross-page footnotes (including references from a
                // footnote to another) without showing unrelated definitions.
                let mut included = std::collections::HashSet::new();
                let mut pending: Vec<_> = footnote_references
                    .iter()
                    .filter(|(_, span)| range.contains(&span.start))
                    .map(|(name, _)| name.clone())
                    .collect();
                while let Some(name) = pending.pop() {
                    if !included.insert(name.clone()) {
                        continue;
                    }
                    if let Some((_, definition)) = footnotes.iter().find(|(id, _)| *id == name) {
                        pending.extend(
                            footnote_references
                                .iter()
                                .filter(|(_, span)| definition.contains(&span.start))
                                .map(|(name, _)| name.clone()),
                        );
                        if definition.start < range.start || definition.end > range.end {
                            markdown.push_str("\n\n");
                            markdown.push_str(&source[definition.clone()]);
                        }
                    }
                }
                ManualPage {
                    heading_level: level,
                    anchor,
                    document_path: "docs/gui.md".to_string(),
                    group: String::new(),
                    title: if continuation == 1 {
                        title
                    } else {
                        format!("{title} (continued {continuation})")
                    },
                    first_line: source[..range.start]
                        .bytes()
                        .filter(|b| *b == b'\n')
                        .count()
                        + 1,
                    source_range: range,
                    markdown,
                }
            })
            .collect();
        let anchors = pages
            .iter()
            .enumerate()
            .filter_map(|(index, page)| page.anchor.clone().map(|anchor| (anchor, index)))
            .collect();
        Self {
            source_hash,
            pages,
            anchors,
        }
    }

    pub fn page_for_anchor(&self, fragment: &str) -> Option<usize> {
        self.anchors.get(&decoded_fragment(fragment)?).copied()
    }

    pub fn from_book(book: &ManualBook) -> Self {
        let mut pages = Vec::new();
        for document in &book.documents {
            let local = Self::new(&book.markdown[document.source_range.clone()], book.source_hash);
            let first_level = local.pages.first().map_or(1, |page| page.heading_level);
            let line_offset = book.markdown[..document.source_range.start]
                .bytes().filter(|&byte| byte == b'\n').count();
            for (index, mut page) in local.pages.into_iter().enumerate() {
                page.document_path.clone_from(&document.path);
                page.group.clone_from(&document.group);
                page.source_range.start += document.source_range.start;
                page.source_range.end += document.source_range.start;
                page.first_line += line_offset;
                page.heading_level = 1 + page.heading_level.saturating_sub(first_level);
                if index == 0 {
                    page.title.clone_from(&document.title);
                }
                pages.push(page);
            }
        }
        Self { source_hash: book.source_hash, pages, anchors: HashMap::new() }
    }

    pub fn page_for_document(&self, path: &str, fragment: &str) -> Option<usize> {
        if path == "docs/gui.md" && !fragment.is_empty() && !self.anchors.is_empty() {
            return self.page_for_anchor(fragment);
        }
        let fragment = decoded_fragment(fragment)?;
        self.pages.iter().position(|page| {
            page.document_path == path && (fragment.is_empty() || page.anchor.as_deref() == Some(&fragment))
        })
    }

    pub fn breadcrumbs(&self, selected: usize) -> Vec<usize> {
        let mut ancestors: Vec<usize> = Vec::new();
        for (index, page) in self.pages.iter().enumerate().take(selected + 1) {
            if page.anchor.is_none() {
                continue;
            }
            while ancestors.last().is_some_and(|&ancestor| {
                self.pages[ancestor].heading_level >= page.heading_level
            }) {
                ancestors.pop();
            }
            ancestors.push(index);
        }
        ancestors
    }

    pub fn page_for_line(&self, line: usize) -> usize {
        self.pages
            .iter()
            .rposition(|page| page.first_line <= line)
            .unwrap_or(0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn manual_book_scopes_duplicate_fragments_and_references_to_their_file() {
        let sources = [
            ManualSource { id: "a", path: "docs/gui/a.md", title: "A", group: "Group", markdown: "## Same\n\n[A][ref]\n\n[ref]: a.md\n" },
            ManualSource { id: "b", path: "docs/gui/b.md", title: "B", group: "Group", markdown: "## Same\n\n[B][ref]\n\n[ref]: b.md\n\n### Same\n\nAgain.\n" },
        ];
        let book = ManualBook::embedded(&sources);
        let pages = ManualPages::from_book(&book);
        let a = pages.page_for_document("docs/gui/a.md", "same").unwrap();
        let b = pages.page_for_document("docs/gui/b.md", "same").unwrap();
        assert_ne!(a, b);
        assert!(pages.pages[a].markdown.contains("[ref]: a.md"));
        assert!(!pages.pages[a].markdown.contains("[ref]: b.md"));
        assert!(pages.page_for_document("docs/gui/b.md", "same-1").is_some());
        assert!(pages.page_for_document("docs/gui/a.md", "same-1").is_none());
        assert_eq!(pages.breadcrumbs(b), vec![b]);
        assert_eq!(relative_document_path("docs/gui/a.md", "../gui/b.md"), Some("docs/gui/b.md".into()));
        assert!(relative_document_path("docs/gui/a.md", "../../../outside.md").is_none());
    }

    #[test]
    fn every_embedded_gui_topic_and_legacy_anchor_has_a_page() {
        let book = ManualBook::embedded(super::super::help_manual_sources::GUI_SOURCES);
        let pages = ManualPages::from_book(&book);
        for source in super::super::help_manual_sources::GUI_SOURCES {
            assert!(pages.page_for_document(source.path, "").is_some(), "{}", source.path);
        }
        for alias in super::super::help_manual_sources::GUI_ALIASES {
            assert!(pages.page_for_document(alias.path, alias.fragment).is_some(), "{}", alias.anchor);
        }
    }

    #[test]
    fn manual_heading_anchors_locate_exact_sections_and_preserve_hierarchy() {
        let source = "# Manual\n\n## DNA & RNA\n\n### **Über** / CDS\n\nFirst.\n\n### **Über** / CDS\n\nSecond.\n\n## End\n\n```md\n## Not a heading\n```\n";
        for source in [source.to_string(), source.replace('\n', "\r\n")] {
            let index = ManualPages::new(&source, 1);
            assert_eq!(index.page_for_anchor("dna--rna"), Some(1));
            assert_eq!(index.page_for_anchor("%C3%BCber--cds"), Some(2));
            assert_eq!(index.page_for_anchor("über--cds-1"), Some(3));
            assert_eq!(index.page_for_anchor("not-a-heading"), None);
            assert_eq!(index.page_for_anchor("%ZZ"), None);
            assert_eq!(index.page_for_anchor("%FF"), None);
            assert_eq!(index.breadcrumbs(3), vec![0, 1, 3]);
            assert_eq!(index.breadcrumbs(4), vec![0, 4]);
            assert!(index.pages[3].markdown.starts_with("### **Über** / CDS"));
            assert!(index.pages[3].markdown.contains("Second."));
        }
    }

    #[test]
    fn gui_manual_pages_keep_source_complete_and_real_pages_small() {
        let book = ManualBook::embedded(super::super::help_manual_sources::GUI_SOURCES);
        let source = &book.markdown;
        let index = ManualPages::from_book(&book);
        assert!(book.documents.len() > 40);
        assert!(index.pages.len() > 40);
        let restored: String = index
            .pages
            .iter()
            .map(|page| &source[page.source_range.clone()])
            .collect();
        assert_eq!(&restored, source);
        let largest = index
            .pages
            .iter()
            .map(|page| page.markdown.len())
            .max()
            .unwrap();
        assert!(
            largest < 40_000,
            "largest indivisible manual page: {largest}"
        );
        assert!(index.pages[0].markdown.len() < PAGE_TARGET_BYTES);
        for (i, page) in index.pages.iter().enumerate() {
            assert_eq!(index.page_for_line(page.first_line), i);
        }
    }

    #[test]
    fn gui_manual_pages_preserve_code_tables_lists_references_unicode_and_crlf() {
        // Synthetic Markdown, recreated in memory; no external fixtures.
        let fixture = format!(
            "# Über GENtle\n\nIntro [jump][later] and footnote[^note].\n\n## First\n\n```md\n## Not a heading\n{}\n```\n\n| A | B |\n|---|---|\n| x | y |\n\n- First item\n  - Nested child\n\n- Second item\n\n## Second\n\nBody[^note].\n\n[later]: ./next.md \"Next chapter\"\n\n[^note]: Retained footnote.\n",
            "x".repeat(PAGE_TARGET_BYTES)
        );
        for source in [fixture.clone(), fixture.replace('\n', "\r\n")] {
            let index = ManualPages::new(&source, 2);
            assert!(!index.pages.iter().any(|p| p.title == "Not a heading"));
            let code = index
                .pages
                .iter()
                .find(|p| p.markdown.contains("```md"))
                .unwrap();
            assert_eq!(code.markdown.matches("```").count(), 2);
            let table = index
                .pages
                .iter()
                .find(|p| p.markdown.contains("| A |"))
                .unwrap();
            assert!(table.markdown.contains("| x | y |"));
            assert!(
                index
                    .pages
                    .iter()
                    .any(|p| p.markdown.contains("- First item")
                        && p.markdown.contains("- Nested child"))
            );
            assert!(index.pages[0].markdown.contains("[later]: ./next.md"));
            assert!(
                index.pages[0]
                    .markdown
                    .contains("[^note]: Retained footnote.")
            );
            assert!(Parser::new(&index.pages[0].markdown).any(|event| matches!(event, Event::Start(Tag::Link { dest_url, .. }) if dest_url.as_ref() == "./next.md")));
            let note = index
                .pages
                .iter()
                .find(|p| p.markdown.contains("Body[^note]"))
                .unwrap();
            assert!(note.markdown.contains("[^note]: Retained footnote."));
            assert!(index.pages[0].title.contains("Über"));
        }
    }

    #[test]
    fn gui_manual_pages_do_not_restart_numbered_lists_or_split_oversized_blocks() {
        let ordered = (0..500)
            .map(|i| format!("1. Numbered step {i}: {}\n", "x".repeat(40)))
            .collect::<String>();
        let source = format!("# Manual\n\n## Ordered\n\n{ordered}\n## End\n\nEnd\n");
        let pages = ManualPages::new(&source, 3);
        let numbered = pages
            .pages
            .iter()
            .find(|p| p.markdown.contains("Numbered step 0:"))
            .unwrap();
        assert!(numbered.markdown.contains("Numbered step 499:"));
        assert_eq!(
            Parser::new(&numbered.markdown)
                .filter(|e| matches!(e, Event::Start(Tag::List(Some(1)))))
                .count(),
            1
        );
        assert_eq!(numbered.markdown.matches("Numbered step").count(), 500);
        assert_eq!(ManualPages::new("", 4).pages.len(), 1);
    }
}
