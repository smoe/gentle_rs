#!/usr/bin/env python3
"""Generate the GENtle ClawBio SKILL.md trigger keyword list from INTENTS.json."""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[1]
SKILL_ROOT = REPO_ROOT / "integrations" / "clawbio" / "skills" / "gentle-cloning"
INTENTS_PATH = SKILL_ROOT / "INTENTS.json"
ADDITIONS_PATH = SKILL_ROOT / "trigger_keyword_additions.json"
SKILL_MD_PATH = SKILL_ROOT / "SKILL.md"


def _normalise_keyword(value: str) -> str:
    return re.sub(r"\s+", " ", value.strip()).lower()


def _append_keyword(keywords: list[str], seen: set[str], value: Any) -> None:
    if not isinstance(value, str):
        return
    keyword = _normalise_keyword(value)
    if keyword and keyword not in seen:
        seen.add(keyword)
        keywords.append(keyword)


def generated_keywords(
    intents_path: Path = INTENTS_PATH,
    additions_path: Path = ADDITIONS_PATH,
) -> list[str]:
    descriptor = json.loads(intents_path.read_text(encoding="utf-8"))
    additions = json.loads(additions_path.read_text(encoding="utf-8"))
    keywords: list[str] = []
    seen: set[str] = set()
    for alias in descriptor.get("aliases", []):
        _append_keyword(keywords, seen, alias)
    for route in descriptor.get("routes", []):
        if not isinstance(route, dict):
            continue
        for term in route.get("trigger_terms", []):
            _append_keyword(keywords, seen, term)
        for alias in route.get("aliases", []):
            _append_keyword(keywords, seen, alias)
    for addition in additions.get("keywords", []):
        _append_keyword(keywords, seen, addition)
    return keywords


def _extract_trigger_keyword_block(text: str) -> tuple[int, int, str]:
    lines = text.splitlines()
    for index, line in enumerate(lines):
        if line.strip() != "trigger_keywords:":
            continue
        indent = line[: len(line) - len(line.lstrip())]
        end = index + 1
        while end < len(lines):
            stripped = lines[end].strip()
            if not stripped:
                end += 1
                continue
            current_indent = lines[end][: len(lines[end]) - len(lines[end].lstrip())]
            if not current_indent.startswith(indent + "  "):
                break
            end += 1
        return index + 1, end, indent + "  "
    raise SystemExit("Could not find trigger_keywords block in SKILL.md")


def render_skill_md_with_keywords(skill_md: str, keywords: list[str]) -> str:
    lines = skill_md.splitlines()
    start, end, item_indent = _extract_trigger_keyword_block(skill_md)
    replacement = [f"{item_indent}- {keyword}" for keyword in keywords]
    trailing_newline = "\n" if skill_md.endswith("\n") else ""
    return "\n".join([*lines[:start], *replacement, *lines[end:]]) + trailing_newline


def current_keywords_from_skill_md(skill_md: str) -> list[str]:
    lines = skill_md.splitlines()
    start, end, _item_indent = _extract_trigger_keyword_block(skill_md)
    keywords: list[str] = []
    for line in lines[start:end]:
        match = re.match(r"\s*-\s+(.*)\s*$", line)
        if match:
            keywords.append(_normalise_keyword(match.group(1)))
    return keywords


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Regenerate or verify gentle-cloning SKILL.md trigger_keywords."
    )
    parser.add_argument("--check", action="store_true", help="Fail when SKILL.md is stale")
    parser.add_argument("--write", action="store_true", help="Rewrite SKILL.md in place")
    args = parser.parse_args()

    keywords = generated_keywords()
    current = SKILL_MD_PATH.read_text(encoding="utf-8")
    expected = render_skill_md_with_keywords(current, keywords)
    if args.check:
        if current != expected:
            current_keywords = current_keywords_from_skill_md(current)
            missing = [keyword for keyword in keywords if keyword not in current_keywords]
            extra = [keyword for keyword in current_keywords if keyword not in keywords]
            print(
                "SKILL.md trigger_keywords are stale; run "
                "`python3 scripts/generate_clawbio_trigger_keywords.py --write`.",
                file=sys.stderr,
            )
            if missing:
                print("Missing generated keywords:", file=sys.stderr)
                for keyword in missing:
                    print(f"  - {keyword}", file=sys.stderr)
            if extra:
                print("Extra SKILL.md keywords:", file=sys.stderr)
                for keyword in extra:
                    print(f"  - {keyword}", file=sys.stderr)
            return 1
        return 0
    if args.write:
        SKILL_MD_PATH.write_text(expected, encoding="utf-8")
        return 0
    for keyword in keywords:
        print(keyword)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
