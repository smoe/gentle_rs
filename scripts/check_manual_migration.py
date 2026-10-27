"""Read-only loss check for a documentation move against an explicit Git source.

This is a migration aid, not a Markdown parser or a permanent prose freeze.
Every nonblank source line must survive with its whitespace and multiplicity.
New headings/navigation may be added. LF and CRLF are equivalent; all other
whitespace (including code/list indentation) is significant. A successful
check proves text retention, not Markdown structure or hyperlink correctness.
Run the separate documentation/checkout checks for those contracts.
"""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import posixpath
import re
import subprocess
from urllib.parse import urlsplit


ROOT = Path(__file__).resolve().parents[1]


def retained_lines(text: str) -> Counter[str]:
    """Count exact content lines; ignore only empty lines and LF/CRLF choice."""
    return Counter(line for line in text.replace("\r\n", "\n").split("\n") if line.strip())


def canonical_links(text: str, document: str, repairs: dict | None = None) -> str:
    """Compare relocated inline resources by their repository target.

    Scope is the same inline-destination syntax used by the existing doc guard;
    fenced examples are kept byte-for-byte. This does not authorize arbitrary
    prose/format normalization or heading deletion.
    """
    result = []
    fence = None
    for line in text.splitlines(keepends=True):
        mark = re.match(r"^ {0,3}(`{3,}|~{3,})", line)
        if mark:
            token = mark[1]
            if fence is None:
                fence = token
            elif token[0] == fence[0] and len(token) >= len(fence):
                fence = None
            result.append(line)
            continue
        if fence is None:
            def canonical(match):
                destination = match[1]
                parts = urlsplit(destination)
                if parts.scheme or parts.netloc or destination.startswith("/"):
                    return match[0]
                path = posixpath.normpath(posixpath.join(posixpath.dirname(document), parts.path)) if parts.path else document
                path = (repairs or {}).get(path, path)
                suffix = ("?" + parts.query if parts.query else "") + ("#" + parts.fragment if parts.fragment else "")
                return f"](repo:{path}{suffix})"
            line = re.sub(r"\]\(([^\s)]+)\)", canonical, line)
        result.append(line)
    return "".join(result)


def compare(source: str, targets: list[str]) -> dict:
    before = retained_lines(source)
    after = Counter()
    for target in targets:
        after.update(retained_lines(target))
    missing = before - after
    return {
        "source_content_lines": sum(before.values()),
        "retained_content_lines": sum((before & after).values()),
        "missing": [{"text": line, "occurrences": count} for line, count in sorted(missing.items())],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-ref", required=True, help="Explicit source commit, before the move")
    parser.add_argument("--source", default="docs/gui.md", help="Source path at that commit")
    parser.add_argument("--relocated-links", action="store_true", help="Compare moved inline links by resolved repository target")
    parser.add_argument("--repair-map", type=Path, help="Explicit reviewed legacy_resources mapping for pre-existing broken links")
    parser.add_argument("targets", nargs="+", type=Path, help="Current files containing the relocated text")
    args = parser.parse_args()
    revision = subprocess.check_output(
        ["git", "rev-parse", "--verify", f"{args.base_ref}^{{commit}}"], cwd=ROOT, text=True,
    ).strip()
    source = subprocess.check_output(["git", "show", f"{revision}:{args.source}"], cwd=ROOT)
    original = source.decode("utf-8")
    targets = [path.read_bytes().decode("utf-8") for path in args.targets]
    repairs = json.loads(args.repair_map.read_text(encoding="utf-8"))["legacy_resources"] if args.repair_map else {}
    if repairs and not args.relocated_links:
        parser.error("--repair-map requires --relocated-links")
    if args.relocated_links:
        original = canonical_links(original, args.source, repairs)
        targets = [canonical_links(text, path.resolve().relative_to(ROOT).as_posix(), repairs) for path, text in zip(args.targets, targets)]
    result = compare(original, targets)
    result.update({
        "source_revision": revision,
        "source_path": args.source,
        "source_sha256": hashlib.sha256(source).hexdigest(),
        "targets": [str(path) for path in args.targets],
        "scope": "Nonblank lines, exact whitespace and multiplicity; not a semantic Markdown check",
        "relocated_links": args.relocated_links,
        "documented_link_repairs": repairs,
    })
    print(json.dumps(result, indent=2, ensure_ascii=False))
    return 1 if result["missing"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
