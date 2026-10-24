//! Presentation-only pagination of a manual at complete Markdown boundaries.
//! The source remains one document; pages never split code, tables or nested items.

use pulldown_cmark::{Event, Options, Parser, Tag, TagEnd};
use std::ops::Range;

pub(super) const PAGE_TARGET_BYTES: usize = 12_000;

pub(super) struct ManualPage {
    pub title: String,
    pub source_range: Range<usize>,
    pub first_line: usize,
    pub markdown: String,
}

pub(super) struct ManualPages {
    pub source_hash: u64,
    pub pages: Vec<ManualPage>,
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
                        if let Tag::Heading { .. } = tag {
                            heading = Some(boundaries.len() - 1);
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
        for (offset, next_heading) in boundaries {
            if offset <= start {
                if let Some(next_heading) = next_heading {
                    title = next_heading;
                }
                continue;
            }
            if next_heading.is_some() || offset - start >= PAGE_TARGET_BYTES {
                ranges.push((start..offset, title.clone(), continuation));
                start = offset;
                if let Some(next_heading) = next_heading {
                    title = next_heading;
                    continuation = 1;
                } else {
                    continuation += 1;
                }
            }
        }
        ranges.push((start..source.len(), title, continuation));
        let pages = ranges
            .into_iter()
            .map(|(range, title, continuation)| {
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
        Self { source_hash, pages }
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
    fn gui_manual_pages_keep_source_complete_and_real_pages_small() {
        let source = include_str!("../../docs/gui.md");
        let index = ManualPages::new(source, 1);
        assert!(index.pages.len() > 40);
        let restored: String = index
            .pages
            .iter()
            .map(|page| &source[page.source_range.clone()])
            .collect();
        assert_eq!(restored, source);
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
