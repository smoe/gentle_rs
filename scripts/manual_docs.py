"""Generate and check the manual index links and embedded source inventory.

Topic Markdown is authoritative; this small catalog only owns ordering, groups
and legacy destinations. No Rust tools or network access are used. The scanner
recognizes the repository's ATX headings outside fenced examples, not arbitrary
CommonMark; actual rendering and anchor handling use pulldown-cmark in GENtle.
"""

from __future__ import annotations

import argparse
import json
import posixpath
from pathlib import Path, PurePosixPath
import re
from urllib.parse import unquote, urlsplit

ROOT = Path(__file__).resolve().parents[1]
CATALOG = "docs/manuals.json"
EMBEDDED = "src/app/help_manual_sources.rs"
START = "<!-- manual-catalog:start -->"
END = "<!-- manual-catalog:end -->"
LEGACY_START = "<!-- manual-legacy:start -->"
LEGACY_END = "<!-- manual-legacy:end -->"
NAV_START = "<!-- manual-page-nav:start -->"
NAV_END = "<!-- manual-page-nav:end -->"
CROSS_START = "<!-- manual-counterparts:start -->"
CROSS_END = "<!-- manual-counterparts:end -->"


def inline_targets(text):
    """Repository inline links outside fenced examples and inline code spans."""
    fence = None
    for line in text.splitlines():
        mark = re.match(r"^ {0,3}(`{3,}|~{3,})", line)
        if mark:
            token = mark[1]
            if fence is None:
                fence = token
            elif token[0] == fence[0] and len(token) >= len(fence):
                fence = None
            continue
        if fence is None:
            for part in re.split(r"(`+[^`]*`+)", line):
                if not part.startswith("`"):
                    yield from re.findall(r"\]\(([^\s)]+)\)", part)


def link_errors(root=ROOT):
    errors = []
    catalog = read_catalog(root)
    for book in books(catalog):
        for relative in paths(book):
            document = root / relative
            for target in inline_targets(document.read_text(encoding="utf-8")):
                url = urlsplit(target)
                if url.scheme or url.netloc:
                    continue
                destination = (document.parent / unquote(url.path)).resolve() if url.path else document
                if not destination.is_file():
                    errors.append(f"{relative}: missing resource {target}")
                elif url.fragment and destination.suffix == ".md":
                    fragments = {row[3] for row in heading_rows(destination.read_text(encoding="utf-8"))}
                    if unquote(url.fragment) not in fragments:
                        errors.append(f"{relative}: missing section {target}")
    return errors


def heading_rows(text):
    """Yield (line, level, raw title, GitHub-style fragment), ignoring fences."""
    fence = None
    used = set()
    for index, line in enumerate(text.splitlines()):
        mark = re.match(r"^ {0,3}(`{3,}|~{3,})", line)
        if mark:
            token = mark[1]
            if fence is None:
                fence = token
            elif token[0] == fence[0] and len(token) >= len(fence):
                fence = None
            continue
        match = re.match(r"^(#{1,6}) (.+?)\s*#*\s*$", line)
        if fence is not None or not match:
            continue
        title = match[2]
        plain = re.sub(r"[`*]", "", title)
        slug = "".join(ch if ch.isalnum() or ch in "_-" else "-" if ch.isspace() else ""
                       for ch in plain.lower())
        anchor = slug
        suffix = 0
        while anchor in used:
            suffix += 1
            anchor = f"{slug}-{suffix}"
        used.add(anchor)
        yield index, len(match[1]), title, anchor


def read_catalog(root=ROOT):
    catalog = json.loads((root / CATALOG).read_text(encoding="utf-8"))
    if catalog.get("version") != 1:
        raise ValueError("Unsupported manual catalog version")
    return catalog


def books(catalog):
    return catalog["books"]


def paths(book):
    return [book["index"], *(topic["path"] for topic in book["topics"])]


def book_text(book_id, root=ROOT):
    """Read the complete canonical book for documentation-contract checks."""
    book = next(book for book in books(read_catalog(root)) if book["id"] == book_id)
    return "\n\n".join((root / path).read_text(encoding="utf-8") for path in paths(book))


def replace_region(text, start, end, body):
    if text.count(start) != 1 or text.count(end) != 1:
        raise ValueError(f"Expected one generated region: {start}")
    before, rest = text.split(start)
    _, after = rest.split(end)
    return before + start + "\n" + body.rstrip() + "\n" + end + after


def validate(root, catalog):
    seen_books = set()
    for book in books(catalog):
        if book["id"] in seen_books:
            raise ValueError("Duplicate manual book id")
        seen_books.add(book["id"])
        if book["id"] not in ("gui", "cli"):
            raise ValueError("Unknown manual book")
        if book["index"] != f"docs/{book['id']}.md":
            raise ValueError("Manual index must retain its public path")
        ids = [topic["id"] for topic in book["topics"]]
        if len(ids) != len(set(ids)) or len(paths(book)) != len(set(paths(book))):
            raise ValueError("Duplicate topic identity/path")
        actual = {path.relative_to(root).as_posix() for path in (root / "docs" / book["id"]).rglob("*.md")}
        expected = {topic["path"] for topic in book["topics"]}
        if actual != expected:
            raise ValueError(f"Unlisted or missing {book['id']} topics: {sorted(actual ^ expected)}")
        anchors = {}
        for path in paths(book):
            relative = PurePosixPath(path)
            if relative.is_absolute() or ".." in relative.parts or not path.startswith("docs/"):
                raise ValueError(f"Unsafe manual path: {path}")
            anchors[path] = {row[3] for row in heading_rows((root / path).read_text(encoding="utf-8"))}
        legacy = [entry["anchor"] for entry in book["legacy"]]
        if len(legacy) != len(set(legacy)):
            raise ValueError("Duplicate legacy fragment")
        for alias in book["legacy"]:
            if alias["path"] not in anchors or alias["fragment"] not in anchors[alias["path"]]:
                raise ValueError(f"Missing legacy destination: {alias}")
    topics = {(book["id"], topic["id"]): topic for book in books(catalog) for topic in book["topics"]}
    seen_pairs = set()
    linked = set()
    for pair in catalog["counterparts"]:
        if len(pair) != 2 or tuple(pair) in seen_pairs:
            raise ValueError("Duplicate or invalid manual counterpart pair")
        seen_pairs.add(tuple(pair))
        for book, topic in zip(("gui", "cli"), pair):
            if (book, topic) not in topics:
                raise ValueError(f"Unknown counterpart: {book}/{topic}")
            linked.add((book, topic))
    for key, topic in topics.items():
        reason = topic.get("no_counterpart_reason", "").strip()
        if (key in linked) == bool(reason):
            raise ValueError(f"Topic needs either counterparts or an explicit absence reason: {key}")


def generated(root=ROOT):
    catalog = read_catalog(root)
    validate(root, catalog)
    rust = [
        "//! Generated manual source inventory. Edit docs/manuals.json and topic Markdown.",
        "//! Regenerate with python3 scripts/manual_docs.py --write.",
        "", "use super::help_manual::{ManualAlias, ManualSource};", "",
        'pub(super) const CATALOG_JSON: &str = include_str!("../../docs/manuals.json");', "",
    ]
    outputs = {}
    for book in books(catalog):
        name = book["id"].upper()
        rust.append(f"pub(super) const {name}_SOURCES: &[ManualSource] = &[")
        sources = [{"id": "index", "path": book["index"], "title": book["title"], "group": "Start here"}, *book["topics"]]
        for topic in sources:
            rust.extend(["    ManualSource {", *[
                f"        {key}: {json.dumps(topic[key], ensure_ascii=False)},"
                for key in ("id", "path", "title", "group")
            ], (f'        markdown: super::{name}_MANUAL_MD,' if topic["id"] == "index"
                else f'        markdown: include_str!("../../{topic["path"]}"),'), "    },"])
        rust.extend(["];", "", f"pub(super) const {name}_ALIASES: &[ManualAlias] = &["])
        for entry in book["legacy"]:
            rust.extend(["    ManualAlias {", *[
                f"        {key}: {json.dumps(entry[key], ensure_ascii=False)},"
                for key in ("anchor", "path", "fragment")
            ], "    },"])
        rust.extend(["];", ""])
        rows = ["## Topics", ""]
        last_group = None
        for topic in book["topics"]:
            if topic["group"] != last_group:
                rows.extend([f"### {topic['group']}", ""])
                last_group = topic["group"]
            relative = PurePosixPath(topic["path"]).relative_to("docs")
            rows.append(f"- [{topic['title']}]({relative})")
        old_links = ["## Earlier Section Links", "", "These entry points retain older bookmarks. Follow the link to the full topic.", ""]
        for alias in book["legacy"]:
            relative = PurePosixPath(alias["path"]).relative_to("docs")
            # The original heading preserves GitHub's public fragment as well
            # as a meaningful no-JavaScript fallback for external bookmarks.
            old_links.extend([f"{alias['heading']}", "", f"[Open this topic]({relative}#{alias['fragment']}).", ""])
        text = (root / book["index"]).read_text(encoding="utf-8")
        text = replace_region(text, START, END, "\n".join(rows))
        text = replace_region(text, LEGACY_START, LEGACY_END, "\n".join(old_links))
        outputs[book["index"]] = text.encode("utf-8")
        for index, topic in enumerate(book["topics"]):
            navigation = [f"[{book['title']}](../{book['id']}.md)"]
            if index:
                previous = book["topics"][index - 1]
                navigation.append(f"[Previous: {previous['title']}]({PurePosixPath(previous['path']).name})")
            if index + 1 < len(book["topics"]):
                following = book["topics"][index + 1]
                navigation.append(f"[Next: {following['title']}]({PurePosixPath(following['path']).name})")
            text = (root / topic["path"]).read_text(encoding="utf-8")
            if CROSS_START not in text:
                first, rest = text.split("\n", 1)
                text = first + f"\n\n{CROSS_START}\n{CROSS_END}\n" + rest
            other_id = "cli" if book["id"] == "gui" else "gui"
            position = 0 if book["id"] == "gui" else 1
            counterpart_ids = [pair[1-position] for pair in catalog["counterparts"] if pair[position] == topic["id"]]
            other = next(item for item in books(catalog) if item["id"] == other_id)
            counterparts = [item for item in other["topics"] if item["id"] in counterpart_ids]
            if counterparts:
                links = [f"[{item['title']}]({posixpath.relpath(item['path'], posixpath.dirname(topic['path']))})" for item in counterparts]
                cross = f"Related {other_id.upper()} reference: " + " | ".join(links) + ".\nThese topic links do not imply identical feature coverage across interfaces."
            else:
                cross = f"No direct {other_id.upper()} counterpart: {topic['no_counterpart_reason']}"
            text = replace_region(text, CROSS_START, CROSS_END, cross)
            if NAV_START not in text:
                text = text.rstrip() + f"\n\n{NAV_START}\n{NAV_END}\n"
            text = replace_region(text, NAV_START, NAV_END, "\n---\n\n" + " | ".join(navigation) + "\n")
            outputs[topic["path"]] = text.encode("utf-8")
    outputs[EMBEDDED] = ("\n".join(rust)).encode("utf-8")
    return outputs


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true", help="Regenerate index regions and embedded source inventory")
    args = parser.parse_args()
    failures = []
    for relative, expected in generated().items():
        path = ROOT / relative
        if args.write:
            path.write_bytes(expected)
        elif not path.is_file() or path.read_bytes() != expected:
            failures.append(relative)
    if failures:
        print("Stale manual projections: " + ", ".join(failures))
        return 1
    print("OK: manual catalog, index links and embedded source inventory")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
