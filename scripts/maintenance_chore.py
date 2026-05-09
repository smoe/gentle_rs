#!/usr/bin/env python3
"""Run low-noise repository maintenance chore checks.

The first implemented chore is `session-close`, matching
docs/maintenance_chore_plan.md. The script intentionally reports concrete
evidence and avoids guessing whether a dirty file is related to the current
session.
"""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable, Optional


ROADMAP_MIN_LINES = 200
ROADMAP_MAX_LINES = 400
COMPLETED_HISTORY_RE = re.compile(
    r"^\s*[-*]\s*(Done|shipped|baseline is now available|already implemented)",
    re.IGNORECASE,
)
CONFLICT_MARKER_RE = re.compile(r"^(<<<<<<< .+|={7}|>>>>>>> .+)$")
GENERATED_CANDIDATE_EXTENSIONS = {
    ".csv",
    ".fasta",
    ".fastq",
    ".json",
    ".mmd",
    ".png",
    ".svg",
    ".tsv",
    ".txt",
}

SHARED_CONTRACT_PREFIXES = (
    "src/engine/",
    "src/engine_shell/",
    "crates/gentle-protocol/",
    "crates/gentle-shell/",
)
SHARED_CONTRACT_EXACT = {
    "src/engine.rs",
    "src/engine_shell.rs",
    "src/lib.rs",
}
TEST_PREFIXES = (
    "tests/",
    "integrations/python/tests/",
)
TEST_EXACT = {
    "src/test_support.rs",
}
LOCAL_MARKDOWN_LINK_RE = re.compile(r"\[[^\]]+\]\(([^)\s]+)(?:\s+\"[^\"]*\")?\)")


@dataclass
class Check:
    name: str
    status: str
    evidence: list[str]
    findings: list[str]
    minimal_fix: list[str]


@dataclass
class ChangedPath:
    status: str
    path: str


def repo_root_from_script() -> Path:
    return Path(__file__).resolve().parents[1]


def read_lines(path: Path) -> list[str]:
    return path.read_text(encoding="utf-8").splitlines()


def relative(repo_root: Path, path: Path) -> str:
    try:
        return str(path.relative_to(repo_root))
    except ValueError:
        return str(path)


def run_git(repo_root: Path, *args: str) -> tuple[int, str, str]:
    proc = subprocess.run(
        ["git", *args],
        cwd=repo_root,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    return proc.returncode, proc.stdout, proc.stderr


def all_doc_paths(repo_root: Path) -> Iterable[Path]:
    docs = repo_root / "docs"
    if docs.exists():
        for path in docs.rglob("*"):
            if path.is_file():
                yield path
    agents = repo_root / "AGENTS.md"
    if agents.exists():
        yield agents


def check_roadmap_line_count(repo_root: Path) -> Check:
    path = repo_root / "docs" / "roadmap.md"
    lines = read_lines(path)
    count = len(lines)
    ok = ROADMAP_MIN_LINES <= count <= ROADMAP_MAX_LINES
    return Check(
        name="roadmap line count",
        status="ok" if ok else "fail",
        evidence=[f"{relative(repo_root, path)} has {count} lines"],
        findings=[]
        if ok
        else [
            f"Expected {ROADMAP_MIN_LINES}-{ROADMAP_MAX_LINES} lines; found {count}."
        ],
        minimal_fix=[]
        if ok
        else [
            "Move completed history to docs/CHANGELOG.md, durable rules to "
            "docs/decisions.md, or condense roadmap-only context."
        ],
    )


def check_completed_history_bullets(repo_root: Path) -> Check:
    path = repo_root / "docs" / "roadmap.md"
    matches = [
        f"{relative(repo_root, path)}:{idx}: {line.strip()}"
        for idx, line in enumerate(read_lines(path), start=1)
        if COMPLETED_HISTORY_RE.search(line)
    ]
    return Check(
        name="roadmap completed-history bullets",
        status="ok" if not matches else "fail",
        evidence=[f"scanned {relative(repo_root, path)} with completed-history regex"],
        findings=matches,
        minimal_fix=[]
        if not matches
        else ["Move completed-history bullets to docs/CHANGELOG.md in this session."],
    )


def check_maintenance_links(repo_root: Path) -> Check:
    required = {
        "docs/roadmap.md": [
            "CHANGELOG.md",
            "decisions.md",
            "Do not add a `Done` entry",
            "Move completed work directly",
        ],
        "docs/architecture.md": [
            "docs/decisions.md",
            "docs/CHANGELOG.md",
        ],
        "AGENTS.md": [
            "docs/maintenance_chore_plan.md",
            "Do not add a `Done` entry",
            "docs/CHANGELOG.md",
        ],
    }
    must_exist = [
        "docs/roadmap.md",
        "docs/CHANGELOG.md",
        "docs/decisions.md",
        "docs/maintenance_chore_plan.md",
        "docs/architecture.md",
        "AGENTS.md",
    ]
    findings: list[str] = []
    evidence: list[str] = []
    for rel_path in must_exist:
        path = repo_root / rel_path
        if path.exists():
            evidence.append(f"found {rel_path}")
        else:
            findings.append(f"{rel_path} is missing.")
    for rel_path, snippets in required.items():
        path = repo_root / rel_path
        if not path.exists():
            continue
        text = path.read_text(encoding="utf-8")
        evidence.append(f"checked invariant links in {rel_path}")
        for snippet in snippets:
            if snippet not in text:
                findings.append(f"{rel_path} is missing required text: {snippet!r}")
    return Check(
        name="maintenance invariant links",
        status="ok" if not findings else "fail",
        evidence=evidence,
        findings=findings,
        minimal_fix=[]
        if not findings
        else [
            "Restore the roadmap/changelog/decisions maintenance invariant links "
            "in docs/roadmap.md, docs/architecture.md, and AGENTS.md."
        ],
    )


def check_conflict_markers(repo_root: Path) -> Check:
    findings: list[str] = []
    checked = 0
    for path in all_doc_paths(repo_root):
        try:
            lines = read_lines(path)
        except UnicodeDecodeError:
            continue
        checked += 1
        for idx, line in enumerate(lines, start=1):
            if CONFLICT_MARKER_RE.match(line):
                findings.append(f"{relative(repo_root, path)}:{idx}: {line}")
    return Check(
        name="conflict markers",
        status="ok" if not findings else "fail",
        evidence=[f"scanned {checked} docs/AGENTS text files"],
        findings=findings,
        minimal_fix=[] if not findings else ["Resolve conflict markers before handoff."],
    )


def parse_status_porcelain(output: str) -> list[str]:
    return [line for line in output.splitlines() if line.strip()]


def status_path(line: str) -> str:
    # Porcelain v1 status is two columns plus a space. Rename/copy entries may
    # contain "old -> new"; the final path is the relevant handoff candidate.
    path = line[3:] if len(line) > 3 else ""
    if " -> " in path:
        path = path.split(" -> ", 1)[1]
    return path.strip()


def parse_diff_name_status(output: str) -> list[ChangedPath]:
    changes: list[ChangedPath] = []
    for line in output.splitlines():
        if not line.strip():
            continue
        fields = line.split("\t")
        if len(fields) < 2:
            continue
        status = fields[0][0]
        path = fields[-1]
        changes.append(ChangedPath(status=status, path=path))
    return changes


def is_shared_contract_path(path: str) -> bool:
    return path in SHARED_CONTRACT_EXACT or path.startswith(SHARED_CONTRACT_PREFIXES)


def is_test_evidence_path(path: str) -> bool:
    if path in TEST_EXACT or path.startswith(TEST_PREFIXES):
        return True
    name = Path(path).name
    return (
        name == "tests.rs"
        or name.startswith("test_")
        or name.endswith("_test.rs")
        or name.endswith("_tests.rs")
    )


def markdown_section(lines: list[str], heading: str) -> list[str]:
    heading_line = f"## {heading}"
    start: Optional[int] = None
    for idx, line in enumerate(lines):
        if line.strip() == heading_line:
            start = idx
            break
    if start is None:
        return []

    end = len(lines)
    for idx in range(start + 1, len(lines)):
        line = lines[idx]
        if line.startswith("## "):
            end = idx
            break
    return lines[start:end]


def has_heading(lines: list[str], heading: str) -> bool:
    return any(line.strip() == f"## {heading}" for line in lines)


def release_gate_section(repo_root: Path) -> tuple[Path, list[str]]:
    roadmap = repo_root / "docs" / "roadmap.md"
    if not roadmap.exists():
        return roadmap, []
    return roadmap, markdown_section(read_lines(roadmap), "Release Gate")


def local_markdown_links(lines: list[str]) -> list[str]:
    links: list[str] = []
    for line in lines:
        for target in LOCAL_MARKDOWN_LINK_RE.findall(line):
            target = target.strip()
            if (
                target.startswith("http://")
                or target.startswith("https://")
                or target.startswith("#")
                or "://" in target
            ):
                continue
            path_part = target.split("#", 1)[0]
            if path_part:
                links.append(path_part)
    return links


def check_git_status(repo_root: Path) -> Check:
    code, stdout, stderr = run_git(repo_root, "status", "--short")
    if code != 0:
        return Check(
            name="git status",
            status="fail",
            evidence=["git status --short failed"],
            findings=[stderr.strip() or "unknown git status error"],
            minimal_fix=["Run git status manually and resolve repository access issues."],
        )

    lines = parse_status_porcelain(stdout)
    generated = []
    for line in lines:
        path_text = status_path(line)
        if Path(path_text).suffix.lower() in GENERATED_CANDIDATE_EXTENSIONS:
            generated.append(line)

    findings: list[str] = []
    status = "ok"
    if lines:
        status = "warn"
        findings.append("Working tree is not clean; confirm each entry is intentional.")
    if generated:
        findings.append(
            "Generated-looking dirty/untracked artifacts are present: "
            + "; ".join(generated[:12])
            + ("; ..." if len(generated) > 12 else "")
        )

    evidence = ["git status --short"]
    if lines:
        evidence.extend(lines[:20])
        if len(lines) > 20:
            evidence.append(f"... {len(lines) - 20} more status entries")
    else:
        evidence.append("working tree clean")

    return Check(
        name="worktree and generated artifacts",
        status=status,
        evidence=evidence,
        findings=findings,
        minimal_fix=[]
        if not findings
        else [
            "Call out unrelated dirty files in handoff; document reproduction "
            "commands or release-attachment intent for generated artifacts."
        ],
    )


def check_plan_fidelity(repo_root: Path, plan: Optional[str]) -> Check:
    if not plan:
        return Check(
            name="plan fidelity reminder",
            status="warn",
            evidence=["no --plan argument supplied"],
            findings=[
                "Plan fidelity is a manual check; compare this session's diff "
                "against the user request, issue, PR, or plan before handoff."
            ],
            minimal_fix=[
                "Pass --plan PATH_OR_LABEL when a concrete plan document or issue "
                "exists, and note any divergence in docs/CHANGELOG.md or handoff."
            ],
        )

    plan_path = Path(plan)
    if not plan_path.is_absolute():
        plan_path = repo_root / plan_path
    evidence = [f"plan reference: {plan}"]
    findings: list[str] = []
    if plan_path.exists():
        evidence.append(f"found plan file {relative(repo_root, plan_path)}")
    else:
        evidence.append("plan reference is treated as a label, not a file path")

    return Check(
        name="plan fidelity reminder",
        status="warn",
        evidence=evidence,
        findings=[
            "Manual confirmation still required: verify the diff matches the "
            "referenced plan/request and record intentional divergence."
        ],
        minimal_fix=[
            "If the diff diverged, update docs/CHANGELOG.md or handoff notes with "
            "the reason instead of silently broadening scope."
        ],
    )


def session_close(repo_root: Path, plan: Optional[str]) -> list[Check]:
    return [
        check_roadmap_line_count(repo_root),
        check_completed_history_bullets(repo_root),
        check_maintenance_links(repo_root),
        check_conflict_markers(repo_root),
        check_git_status(repo_root),
        check_plan_fidelity(repo_root, plan),
    ]


def collect_changed_paths(repo_root: Path, base_ref: str) -> tuple[Optional[Check], list[ChangedPath]]:
    code, stdout, stderr = run_git(repo_root, "diff", "--name-status", f"{base_ref}...HEAD")
    if code != 0:
        return (
            Check(
                name="changed path collection",
                status="fail",
                evidence=[f"git diff --name-status {base_ref}...HEAD failed"],
                findings=[stderr.strip() or "unknown git diff error"],
                minimal_fix=[
                    "Pass a valid --base-ref that is reachable from the current "
                    "branch, or fetch the base ref before running drift-scan."
                ],
            ),
            [],
        )
    return None, parse_diff_name_status(stdout)


def check_shared_contract_changes_have_tests(changes: list[ChangedPath]) -> Check:
    paths = [change.path for change in changes]
    shared_contract_paths = sorted(path for path in paths if is_shared_contract_path(path))
    test_paths = sorted(path for path in paths if is_test_evidence_path(path))

    evidence = [
        f"changed paths scanned: {len(paths)}",
        "shared-contract paths: "
        + (", ".join(shared_contract_paths) if shared_contract_paths else "none"),
        "test evidence paths: " + (", ".join(test_paths) if test_paths else "none"),
    ]

    if not shared_contract_paths:
        return Check(
            name="shared contract changes include separate test evidence",
            status="ok",
            evidence=evidence,
            findings=[],
            minimal_fix=[],
        )
    if test_paths:
        return Check(
            name="shared contract changes include separate test evidence",
            status="ok",
            evidence=evidence,
            findings=[],
            minimal_fix=[],
        )
    return Check(
        name="shared contract changes include separate test evidence",
        status="warn",
        evidence=evidence,
        findings=[
            "Shared engine/protocol/shell-contract paths changed, but this diff "
            "does not include obvious separate test-file evidence. This is a "
            "review prompt, not proof of a bug; inline Rust test changes may "
            "still need manual confirmation."
        ],
        minimal_fix=[
            "Add or update deterministic tests for the changed contract behavior, "
            "or document why existing or inline tests already cover these paths."
        ],
    )


def drift_scan(repo_root: Path, base_ref: str) -> list[Check]:
    collection_failure, changes = collect_changed_paths(repo_root, base_ref)
    checks: list[Check] = []
    if collection_failure is not None:
        checks.append(collection_failure)
        return checks

    checks.append(
        Check(
            name="changed path collection",
            status="ok",
            evidence=[
                f"git diff --name-status {base_ref}...HEAD",
                f"changed paths: {len(changes)}",
            ]
            + [f"{change.status}\t{change.path}" for change in changes[:30]]
            + ([f"... {len(changes) - 30} more changed paths"] if len(changes) > 30 else []),
            findings=[],
            minimal_fix=[],
        )
    )
    checks.append(check_shared_contract_changes_have_tests(changes))
    return checks


def check_release_gate_sections(repo_root: Path) -> Check:
    roadmap = repo_root / "docs" / "roadmap.md"
    if not roadmap.exists():
        return Check(
            name="release gate roadmap sections",
            status="fail",
            evidence=[f"{relative(repo_root, roadmap)} is missing"],
            findings=["docs/roadmap.md is required for release-gate readiness."],
            minimal_fix=["Restore docs/roadmap.md before running release-gate checks."],
        )
    lines = read_lines(roadmap)
    required = ["Release Gate", "Next Session Priorities"]
    missing = [heading for heading in required if not has_heading(lines, heading)]
    return Check(
        name="release gate roadmap sections",
        status="ok" if not missing else "fail",
        evidence=[f"checked headings in {relative(repo_root, roadmap)}"],
        findings=[f"Missing heading: ## {heading}" for heading in missing],
        minimal_fix=[] if not missing else ["Restore the compact roadmap release-gate structure."],
    )


def check_release_gate_links(repo_root: Path) -> Check:
    roadmap, section = release_gate_section(repo_root)
    if not section:
        return Check(
            name="release gate local links",
            status="fail",
            evidence=[f"Release Gate section not found in {relative(repo_root, roadmap)}"],
            findings=["Cannot verify release-gate links without a Release Gate section."],
            minimal_fix=["Restore ## Release Gate in docs/roadmap.md."],
        )

    docs_dir = repo_root / "docs"
    links = local_markdown_links(section)
    missing: list[str] = []
    resolved: list[str] = []
    for link in links:
        target = (docs_dir / link).resolve()
        resolved_label = f"{link} -> {relative(repo_root, target)}"
        if target.exists():
            resolved.append(resolved_label)
        else:
            missing.append(resolved_label)

    return Check(
        name="release gate local links",
        status="ok" if not missing else "fail",
        evidence=[f"local links in Release Gate: {len(links)}"]
        + (resolved if resolved else ["no local links found"]),
        findings=[f"Missing release-gate link target: {item}" for item in missing],
        minimal_fix=[] if not missing else ["Restore or retarget missing release-gate documents."],
    )


def check_release_docs_present(repo_root: Path) -> Check:
    required_paths = [
        repo_root / "INSTALL.md",
        repo_root / "README.md",
        repo_root / "docs" / "release.md",
    ]
    release_notes = sorted(repo_root.glob("release_notes_v*.md"))
    missing = [relative(repo_root, path) for path in required_paths if not path.exists()]
    evidence = [f"found {relative(repo_root, path)}" for path in required_paths if path.exists()]
    evidence.append(
        "release notes: "
        + (", ".join(relative(repo_root, path) for path in release_notes) if release_notes else "none")
    )
    if not release_notes:
        missing.append("release_notes_v*.md")

    return Check(
        name="release-facing docs present",
        status="ok" if not missing else "fail",
        evidence=evidence,
        findings=[f"Missing release-facing doc: {path}" for path in missing],
        minimal_fix=[] if not missing else ["Restore release-facing docs before release handoff."],
    )


def check_cargo_gate_alignment(repo_root: Path) -> Check:
    roadmap, section = release_gate_section(repo_root)
    release_doc = repo_root / "docs" / "release.md"
    roadmap_text = "\n".join(section)
    evidence = [f"checked {relative(repo_root, roadmap)} Release Gate"]

    if "cargo check" not in roadmap_text:
        return Check(
            name="cargo check release alignment",
            status="ok",
            evidence=evidence + ["roadmap Release Gate does not mention cargo check"],
            findings=[],
            minimal_fix=[],
        )
    if not release_doc.exists():
        return Check(
            name="cargo check release alignment",
            status="fail",
            evidence=evidence + [f"{relative(repo_root, release_doc)} is missing"],
            findings=["Roadmap requires cargo check, but docs/release.md is missing."],
            minimal_fix=["Restore docs/release.md or move release-signoff guidance to an existing doc."],
        )

    release_text = release_doc.read_text(encoding="utf-8")
    if "cargo check" in release_text:
        return Check(
            name="cargo check release alignment",
            status="ok",
            evidence=evidence + [f"{relative(repo_root, release_doc)} contains cargo check"],
            findings=[],
            minimal_fix=[],
        )
    return Check(
        name="cargo check release alignment",
        status="warn",
        evidence=evidence
        + [
            "roadmap Release Gate contains cargo check",
            f"{relative(repo_root, release_doc)} does not contain cargo check",
        ],
        findings=[
            "Roadmap release gate names `cargo check`, but docs/release.md does "
            "not mirror that release-signoff command. This is a review prompt, "
            "not proof that the release process is broken."
        ],
        minimal_fix=[
            "Add `cargo check -q` to the Local Pre-Tag Smoke Checklist in docs/release.md, "
            "or update the roadmap if the release gate changed."
        ],
    )


def release_gate(repo_root: Path) -> list[Check]:
    return [
        check_release_gate_sections(repo_root),
        check_release_gate_links(repo_root),
        check_release_docs_present(repo_root),
        check_cargo_gate_alignment(repo_root),
    ]


def status_counts(checks: list[Check]) -> dict[str, int]:
    counts = {"ok": 0, "warn": 0, "fail": 0}
    for check in checks:
        counts[check.status] = counts.get(check.status, 0) + 1
    return counts


def render_markdown(checks: list[Check], chore: str) -> str:
    counts = status_counts(checks)
    lines = [
        f"# Maintenance Chore: {chore}",
        "",
        f"Summary: {counts['ok']} ok, {counts['warn']} warn, {counts['fail']} fail.",
        "",
        "## Evidence And Findings",
    ]
    for check in checks:
        lines.extend(["", f"### {check.status.upper()}: {check.name}", "", "Evidence:"])
        lines.extend(f"- {item}" for item in (check.evidence or ["not recorded"]))
        if check.findings:
            lines.extend(["", "Findings:"])
            lines.extend(f"- {item}" for item in check.findings)
        else:
            lines.extend(["", "Findings:", "- No finding."])
        if check.minimal_fix:
            lines.extend(["", "Minimal fix:"])
            lines.extend(f"- {item}" for item in check.minimal_fix)
    if counts["fail"] == 0 and counts["warn"] == 0:
        lines.extend(["", "## No-Finding Note", "", f"- No warnings or failures in {chore}."])
    return "\n".join(lines) + "\n"


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run GENtle maintenance chores.")
    parser.add_argument(
        "--repo-root",
        type=Path,
        default=repo_root_from_script(),
        help="Repository root. Defaults to the parent of this script directory.",
    )
    parser.add_argument(
        "--format",
        choices=("markdown", "json"),
        default="markdown",
        help="Output format.",
    )
    subparsers = parser.add_subparsers(dest="chore", required=True)

    session = subparsers.add_parser(
        "session-close",
        help="Run session-close hygiene checks from docs/maintenance_chore_plan.md.",
    )
    session.add_argument(
        "--format",
        choices=("markdown", "json"),
        default=argparse.SUPPRESS,
        help="Output format.",
    )
    session.add_argument(
        "--plan",
        help="Optional plan file path or label used for manual plan-fidelity evidence.",
    )

    drift = subparsers.add_parser(
        "drift-scan",
        help="Run low-noise shared-contract drift checks for a commit range.",
    )
    drift.add_argument(
        "--format",
        choices=("markdown", "json"),
        default=argparse.SUPPRESS,
        help="Output format.",
    )
    drift.add_argument(
        "--base-ref",
        required=True,
        help="Base ref for git diff --name-status BASE...HEAD.",
    )

    release = subparsers.add_parser(
        "release-gate",
        help="Run release-gate readiness checks against roadmap and release docs.",
    )
    release.add_argument(
        "--format",
        choices=("markdown", "json"),
        default=argparse.SUPPRESS,
        help="Output format.",
    )

    return parser.parse_args(argv)


def main(argv: list[str]) -> int:
    args = parse_args(argv)
    repo_root = args.repo_root.resolve()
    if args.chore == "session-close":
        checks = session_close(repo_root, args.plan)
    elif args.chore == "drift-scan":
        checks = drift_scan(repo_root, args.base_ref)
    elif args.chore == "release-gate":
        checks = release_gate(repo_root)
    else:
        raise AssertionError(f"unknown chore {args.chore}")

    if args.format == "json":
        print(
            json.dumps(
                {
                    "chore": args.chore,
                    "summary": status_counts(checks),
                    "checks": [asdict(check) for check in checks],
                },
                indent=2,
                sort_keys=True,
            )
        )
    else:
        print(render_markdown(checks, args.chore), end="")

    return 1 if any(check.status == "fail" for check in checks) else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
