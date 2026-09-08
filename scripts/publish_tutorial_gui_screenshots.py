#!/usr/bin/env python3
"""Publish selected, provenance-bound tutorial screenshots from a live GUI run."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from pathlib import Path


EVIDENCE_SCHEMA = "gentle.tutorial_gui_screenshot_evidence.v1"
SELECTION_SCHEMA = "gentle.tutorial_gui_screenshot_selection.v1"
PUBLICATION_SCHEMA = "gentle.tutorial_gui_screenshot_publication.v1"
CHECKPOINT_SUFFIXES = (
    "raw.png",
    "orientation.svg",
    "context.svg",
    "screenshot.json",
    "snapshot.json",
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def write_json(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    text = json.dumps(value, indent=2, sort_keys=True, ensure_ascii=False) + "\n"
    path.write_text(text, encoding="utf-8")


def require(condition: bool, message: str) -> None:
    if not condition:
        raise SystemExit(message)


def publish(evidence_root: Path, selection_path: Path, output_root: Path) -> None:
    selection = load_json(selection_path)
    require(selection.get("schema") == SELECTION_SCHEMA, "unsupported selection schema")
    report_path = evidence_root / "acceptance-report.json"
    report = load_json(report_path)
    require(report.get("status") == "pass", "source acceptance run is not green")

    expected_revision = selection["source_revision"]
    rows: list[dict] = []
    shared_binary_hash: str | None = None
    for chapter in selection["chapters"]:
        chapter_id = chapter["chapter_id"]
        for step_id in chapter["checkpoints"]:
            source_dir = evidence_root / chapter_id / "checkpoints"
            record_path = source_dir / f"{step_id}.screenshot.json"
            record = load_json(record_path)
            require(record.get("schema") == EVIDENCE_SCHEMA, f"bad record: {record_path}")
            require(record.get("source_revision") == expected_revision, f"revision mismatch: {record_path}")
            require(record.get("chapter_id") == chapter_id, f"chapter mismatch: {record_path}")
            require(record.get("step_id") == step_id, f"step mismatch: {record_path}")
            if shared_binary_hash is None:
                shared_binary_hash = record["gentle_binary_sha256"]
            require(
                record["gentle_binary_sha256"] == shared_binary_hash,
                f"mixed GUI binaries: {record_path}",
            )

            raw_path = source_dir / f"{step_id}.raw.png"
            require(
                sha256_file(raw_path) == record["capture"]["raw"]["sha256"],
                f"raw screenshot hash mismatch: {raw_path}",
            )
            for derived in record["derived_views"]:
                derived_path = source_dir / Path(derived["path"]).name
                require(
                    sha256_file(derived_path) == derived["sha256"],
                    f"derived screenshot hash mismatch: {derived_path}",
                )
            snapshot_path = source_dir / f"{step_id}.snapshot.json"
            require(
                sha256_file(snapshot_path)
                == record["semantic_snapshot"]["retained_file_sha256"],
                f"semantic snapshot hash mismatch: {snapshot_path}",
            )

            target_dir = output_root / chapter_id
            target_dir.mkdir(parents=True, exist_ok=True)
            files = []
            for suffix in CHECKPOINT_SUFFIXES:
                source = source_dir / f"{step_id}.{suffix}"
                target = target_dir / source.name
                shutil.copyfile(source, target)
                files.append(
                    {
                        "path": target.as_posix(),
                        "sha256": sha256_file(target),
                        "size_bytes": target.stat().st_size,
                    }
                )
            rows.append(
                {
                    "chapter_id": chapter_id,
                    "prose_step": record["prose_step"],
                    "step_id": step_id,
                    "semantic_target": record["requested_target"],
                    "source_record_sha256": sha256_file(record_path),
                    "semantic_snapshot_canonical_sha256": record["semantic_snapshot"][
                        "canonical_sha256"
                    ],
                    "captured_at_unix_ms": record["capture"]["captured_at_unix_ms"],
                    "files": files,
                }
            )

    manifest = {
        "schema": PUBLICATION_SCHEMA,
        "source_revision": expected_revision,
        "gentle_binary_sha256": shared_binary_hash,
        "capture_date": selection["capture_date"],
        "acceptance_report": {
            "source_sha256": sha256_file(report_path),
            "status": report["status"],
        },
        "selection_sha256": sha256_file(selection_path),
        "path_note": (
            "The copied screenshot sidecars are byte-identical capture records. "
            "Their absolute paths describe the original run; this manifest binds "
            "the portable repository copies by relative path and SHA-256."
        ),
        "checkpoints": rows,
    }
    write_json(output_root / "publication-manifest.json", manifest)


def check(output_root: Path) -> None:
    manifest_path = output_root / "publication-manifest.json"
    manifest = load_json(manifest_path)
    require(manifest.get("schema") == PUBLICATION_SCHEMA, "unsupported publication schema")
    for checkpoint in manifest["checkpoints"]:
        for file_row in checkpoint["files"]:
            path = Path(file_row["path"])
            require(path.is_file(), f"missing published screenshot file: {path}")
            require(sha256_file(path) == file_row["sha256"], f"hash mismatch: {path}")
            require(path.stat().st_size == file_row["size_bytes"], f"size mismatch: {path}")
        chapter_dir = output_root / checkpoint["chapter_id"]
        raw_name = f"{checkpoint['step_id']}.raw.png"
        for role in ("orientation", "context"):
            svg = chapter_dir / f"{checkpoint['step_id']}.{role}.svg"
            require(raw_name in svg.read_text(encoding="utf-8"), f"broken SVG raw-image link: {svg}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--evidence-root", type=Path)
    parser.add_argument(
        "--selection",
        type=Path,
        default=Path("docs/screenshots/tutorial_gui_acceptance/selection.json"),
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        default=Path("docs/screenshots/tutorial_gui_acceptance"),
    )
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if args.check:
        check(args.output_root)
        return
    require(args.evidence_root is not None, "--evidence-root is required when publishing")
    publish(args.evidence_root, args.selection, args.output_root)
    check(args.output_root)


if __name__ == "__main__":
    main()
