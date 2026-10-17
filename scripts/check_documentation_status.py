#!/usr/bin/env python3
"""Semantic guards for the documentation lifecycle reconciled on 2026-10-07.

Read-only, offline checks: no generated-byte comparison, shell execution, model,
build or scientific acceptance. Text is newline-neutral; fast checkout tests
exercise the preview and catalog contracts through actual LF/CRLF Git clones.
"""

from collections import Counter
import json
from pathlib import Path
import re
import shlex
import tomllib
from urllib.parse import unquote, urlsplit


ROOT = Path(__file__).resolve().parents[1]
RETIRED = (
    "construct_reasoning_inspection_cli_parity_plan.md",
    "contextual_action_availability_review_prompt.md",
    "tss_view_export_followup_codex_prompt.md",
    "internal_12_tutorial_regeneration_glen_prompt.md",
    "gene_assay_study_glen_acceptance.md",
)
STATUS_DOCS = (
    "reviewer_preview.md", "primer_design_specialist_plan.md",
    "cloning_pattern_catalog.md", "cloning_routine_catalog_plan.md",
    "rna_read_origin_sparse_index_plan.md", "rna_seq_nanopore_cloning_regions_plan.md",
    "legacy_gentle_m_intake.md", "sequencing_confirmation_plan.md",
    "gene_assay_study_gui_plan.md", "internal_12_tutorial_baseline_acceptance.md",
    "construct_reasoning_plan.md", "gene_set_retrieval_producers_plan.md",
    "gui_gene_set_collection_operations_plan.md", "roadmap.md",
    "maintenance_chore_plan.md", "clawbio_gentle_integration_onepager.md",
    "release_notes/release_notes_v0.1.0-internal.12.md",
)


def links(text):
    """Only Markdown inline targets, not bare historical filenames in prose."""
    return re.findall(r"\]\(([^\s)]+)\)", text)


def readme_version_errors(root):
    """Mirror the Rust release guard before CI spends time compiling/testing."""
    version = tomllib.loads((root / "Cargo.toml").read_text(encoding="utf-8"))["workspace"]["package"]["version"]
    readme = (root / "README.md").read_text(encoding="utf-8")
    if f"Current package version: `{version}`." not in readme:
        return [f"README project status must name Cargo package version {version}"]
    return []


def preview_errors(root):
    errors = []
    preview = (root / "docs/reviewer_preview.md").read_text(encoding="utf-8")
    version = tomllib.loads((root / "Cargo.toml").read_text(encoding="utf-8"))["workspace"]["package"]["version"]
    if f"GENtle `{version}`" not in preview:
        errors.append("embedded preview does not name the current workspace version")
    if "mktemp -d" not in preview or re.search(r"\brm\s+-", preview):
        errors.append("preview must use a fresh directory without destructive cleanup")
    if "not an acceptance result" not in preview or "acceptance remain pending" not in preview:
        errors.append("preview must distinguish inspection from pending acceptance")
    commands = [shlex.split(line) for line in preview.splitlines()
                if line.startswith("gentle_cli ") and " workflow " in line]
    states = []
    if len(commands) != 2:
        errors.append("preview must retain the two short offline workflow checks")
    for command in commands:
        state = command[command.index("--state") + 1]
        states.append(state)
        path = command[command.index("workflow") + 1]
        example = json.loads((root / path.removeprefix("@")).read_text(encoding="utf-8"))
        if not path.startswith("@") or example["test_mode"] != "always" or "offline" not in example["tags"]:
            errors.append(f"preview workflow is not an explicit offline example: {path}")
        for required in example["required_files"]:
            if not (root / required).is_file():
                errors.append(f"preview workflow input missing: {required}")
        # Compare stated outcomes with the declared workflow, not a guessed runtime result.
        for operation in example["workflow"]["ops"]:
            payload = next(iter(operation.values()))
            if output_id := payload.get("output_id"):
                if f"`{output_id}`" not in preview:
                    errors.append(f"preview omits declared output: {output_id}")
            if "ExtractRegion" in operation and payload["to"] - payload["from"] != 120:
                errors.append("preview extract length no longer matches its workflow")
    if len(set(states)) != len(states) or any(not state.startswith("$WORK/") for state in states):
        errors.append("preview examples must use separate fresh state files")
    for path, marker in [("build.rs", '"docs/reviewer_preview.md"'),
                         ("src/app.rs", 'include_str!("../docs/reviewer_preview.md")'),
                         ("crates/gentle-shell/src/ui_intent.rs", '"pcr-design"')]:
        if marker not in (root / path).read_text(encoding="utf-8"):
            errors.append(f"documented embedded-help/UI route changed: {path}")
    specialist = (root / "docs/primer_design_specialist_plan.md").read_text(encoding="utf-8")
    if "remaining proposal" not in specialist.lower() or "not registered" not in specialist:
        errors.append("persisted primer filtered views must be labelled as a proposal")
    if "PCR Designer" not in specialist or "no dedicated, report-centric" in specialist:
        errors.append("primer plan must recognize the existing PCR Designer")
    return errors


def catalog_errors(root):
    errors = []
    directory = root / "assets/cloning_patterns_catalog"
    actual = Counter(path.relative_to(directory).parts[0] for path in directory.rglob("*.json"))
    text = (root / "docs/cloning_pattern_catalog.md").read_text(encoding="utf-8")
    rows = re.findall(r"^\| ([a-z_]+) \| ([0-9]+) \|$", text, re.MULTILINE)
    if len(rows) != len(dict(rows)) or {name: int(n) for name, n in rows} != dict(actual):
        errors.append("documented template-family counts differ from the JSON inventory")
    routines = json.loads((root / "assets/cloning_routines.json").read_text(encoding="utf-8"))["routines"]
    plan = (root / "docs/cloning_routine_catalog_plan.md").read_text(encoding="utf-8")
    if f"{len(routines)} definitions" not in plan or f"{len(routines)}-entry" not in text:
        errors.append("documented routine count differs from the maintained catalog")
    if f"{sum(actual.values())} templates" not in text:
        errors.append("documented template total differs from the inventory")
    for routine in routines:
        if not (root / routine["template_path"]).is_file():
            errors.append(f"routine template missing: {routine['routine_id']}")
    return errors


def lifecycle_errors(root):
    errors = []
    for name in RETIRED:
        if (root / "docs" / name).exists():
            errors.append(f"completed request restored without lifecycle review: {name}")
    for doc in (root / "docs").rglob("*.md"):
        for target in links(doc.read_text(encoding="utf-8")):
            parts = urlsplit(target)
            # Historical external links do not restore a local retired request.
            if not parts.scheme and not parts.netloc and Path(unquote(parts.path)).name in RETIRED:
                errors.append(f"{doc.relative_to(root)} still links to retired request {target}")
    for relative in STATUS_DOCS:
        doc = root / "docs" / relative
        for target in links(doc.read_text(encoding="utf-8")):
            parts = urlsplit(target)
            if not parts.scheme and parts.path and not (doc.parent / unquote(parts.path)).exists():
                errors.append(f"{doc.relative_to(root)} has unavailable local target {target}")
    rna = (root / "src/engine/analysis/rna_reads.rs").read_text(encoding="utf-8")
    plan = (root / "docs/rna_read_origin_sparse_index_plan.md").read_text(encoding="utf-8")
    warning = "ROI capture layer is planned but not active in phase-1 scoring"
    if (warning in rna) != (warning in plan):
        errors.append("RNA capture implementation/status disagreement requires reconciliation")
    return errors


def check(root=ROOT):
    return readme_version_errors(root) + preview_errors(root) + catalog_errors(root) + lifecycle_errors(root)


if __name__ == "__main__":
    findings = check()
    for finding in findings:
        print(f"FAIL: {finding}")
    if not findings:
        print("OK: documentation status, preview and catalog contracts (not runtime acceptance)")
    raise SystemExit(bool(findings))
