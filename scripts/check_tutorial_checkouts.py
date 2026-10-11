#!/usr/bin/env python3
"""Replay parity and tutorial checks under LF and CRLF Git checkout policies.

Uses local, temporary clones sharing read-only objects; never builds Rust,
fetches references, initializes submodules, or modifies the source checkout.
The binary must already be built from the revision being checked. A local
attributes-only overlay is explicit and is not exact-commit acceptance.
Visual-only checks inspect canonical declarations and retained archive bytes;
they do not capture the desktop or certify current-candidate GUI behavior.
"""

import argparse
import datetime
import hashlib
import json
import math
import os
from pathlib import Path
import re
import subprocess
import tempfile
import xml.etree.ElementTree as ET
import zlib


ROOT = Path(__file__).resolve().parents[1]
MODES = (("lf", "false", "lf"), ("crlf", "true", "crlf"))
VISUAL_ROOTS = ("docs/screenshots/", "docs/figures/",
                "docs/tutorial/generated/artifacts/", "docs/tutorial/reproducibility/")
MAX_SVG_BYTES = 16 * 1024 * 1024


def visual_require(condition, code, detail):
    if not condition:
        raise RuntimeError(f"tutorial visual {code}: {detail}")


def visual_json(path):
    try:
        value = json.loads(path.read_bytes())
    except (OSError, ValueError) as error:
        raise RuntimeError(f"tutorial visual invalid_json: {path}: {error}") from error
    visual_require(isinstance(value, dict), "invalid_json_object", str(path))
    return value


def visual_sha256(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def visual_hex(value, size):
    return isinstance(value, str) and bool(re.fullmatch(r"[0-9a-f]{%d}" % size, value))


def visual_path(root, name, relative_to=None):
    """Admit local documentation images, including resolved Markdown links."""
    visual_require(isinstance(name, str) and bool(name) and not re.search(r"[:\\\\?#]", name)
                   and not name.startswith("/"), "unsafe_path", str(name))
    if relative_to is None:
        visual_require(".." not in name.split("/"), "unsafe_path", name)
    root = root.resolve()
    path = ((relative_to or root) / name).resolve()
    visual_require(path.is_relative_to(root), "unsafe_path", name)
    relative = path.relative_to(root).as_posix()
    visual_require(any(relative.startswith(prefix) for prefix in VISUAL_ROOTS),
                   "unapproved_asset_root", name)
    visual_require(path.is_file(), "missing_asset", relative)
    return path


def inspect_visual_asset(root, path):
    """Bounded structural checks, not image decoding or screenshot acceptance."""
    if path.suffix.lower() == ".png":
        with path.open("rb") as stream:
            header = stream.read(33)
        visual_require(len(header) == 33 and header[:8] == b"\x89PNG\r\n\x1a\n"
                       and header[8:12] == b"\x00\x00\x00\x0d" and header[12:16] == b"IHDR"
                       and zlib.crc32(header[12:29]) == int.from_bytes(header[29:33], "big"),
                       "invalid_png", path.relative_to(root).as_posix())
        dimensions = [int.from_bytes(header[16:20], "big"), int.from_bytes(header[20:24], "big")]
        visual_require(all(value > 0 for value in dimensions), "invalid_png_dimensions", str(path))
        return {"format": "png", "dimensions": dimensions}
    visual_require(path.suffix.lower() == ".svg", "unsupported_format", str(path))
    visual_require(path.stat().st_size <= MAX_SVG_BYTES, "svg_size_budget", str(path))
    payload = path.read_bytes()
    visual_require(b"<!DOCTYPE" not in payload.upper() and b"<!ENTITY" not in payload.upper(),
                   "unsafe_svg_xml", str(path))
    try:
        document = ET.fromstring(payload)
    except ET.ParseError as error:
        raise RuntimeError(f"tutorial visual invalid_svg: {path}: {error}") from error
    visual_require(document.tag in ("svg", "{http://www.w3.org/2000/svg}svg"),
                   "invalid_svg_root", str(path))
    try:
        if "viewBox" in document.attrib:
            box = [float(value) for value in re.split(r"[ ,]+", document.attrib["viewBox"].strip())]
            valid_size = len(box) == 4 and all(math.isfinite(v) for v in box) and all(v > 0 for v in box[2:])
        else:
            pattern = r"(\d+(?:\.\d+)?)(?:px|pt|pc|mm|cm|in)?"
            sizes = [re.fullmatch(pattern, document.get(key, "")) for key in ("width", "height")]
            valid_size = all(match and math.isfinite(float(match[1])) and float(match[1]) > 0 for match in sizes)
    except ValueError:
        valid_size = False
    visual_require(valid_size, "invalid_svg_dimensions", str(path))
    images = []
    for node in document.iter():
        if node.tag.split("}")[-1] != "image":
            continue
        target = node.get("href", node.get("{http://www.w3.org/1999/xlink}href", ""))
        if target.startswith("data:image/"):
            images.append("embedded")
        else:
            image = visual_path(root, target, path.parent)
            visual_require(image.suffix.lower() == ".png", "unsupported_svg_image", target)
            inspect_visual_asset(root, image)
            images.append(image.relative_to(root).as_posix())
    return {"format": "svg", "raster_images": images}


def markdown_visuals(text):
    """Ignore fenced examples, so prose describing image syntax is not a link."""
    text = re.sub(r"(?ms)^(`{3,}|~{3,})[^\n]*\n.*?^\1[ \t]*$", "", text)
    return re.findall(r"!\[(?:\\.|[^\]])*\]\((?:<([^>]+)>|([^\s)]+))(?:\s+\"[^\"]*\")?\)", text)


def screenshot_provenance(root, path, graphic):
    """Inspect the existing archive records; never promote them to current passes."""
    stem = re.sub(r"\.(?:raw\.png|context\.svg|orientation\.svg)$", "", path.name)
    sidecar = path.with_name(stem + ".screenshot.json")
    if sidecar.is_file():
        record = visual_json(sidecar)
        visual_require(record.get("schema") == "gentle.tutorial_gui_screenshot_evidence.v1",
                       "invalid_capture_schema", str(sidecar))
        for field, size in (("source_revision", 40), ("gentle_binary_sha256", 64),
                            ("tutorial_manifest_sha256", 64), ("step_sha256", 64),
                            ("acceptance_contract_sha256", 64)):
            visual_require(visual_hex(record.get(field), size),
                           "missing_capture_binding", f"{sidecar}: {field}")
        visual_require(isinstance(record.get("chapter_id"), str) and bool(record["chapter_id"])
                       and isinstance(record.get("step_id"), str) and bool(record["step_id"])
                       and type(record.get("prose_step")) is int and record["prose_step"] > 0,
                       "missing_capture_step", str(sidecar))
        capture, derived = record.get("capture"), record.get("derived_views", [])
        visual_require(isinstance(capture, dict) and isinstance(capture.get("raw"), dict)
                       and isinstance(derived, list), "missing_capture_binding", str(sidecar))
        raw = capture["raw"]
        expected = [raw] + derived
        visual_require(all(isinstance(row, dict) and isinstance(row.get("path"), str)
                           and bool(row["path"]) and visual_hex(row.get("sha256"), 64)
                           for row in expected), "missing_capture_binding", str(sidecar))
        match = next((row for row in expected if Path(row.get("path", "")).name == path.name), None)
        visual_require(match is not None, "capture_path_mismatch", str(path))
        # A raw-only publication need not copy every optional private crop.
        # The selected asset and its authoritative raw PNG must remain bound.
        for row in (raw, match) if match != raw else (raw,):
            archive_path = visual_path(root, Path(row.get("path", "")).name, path.parent)
            inspect_visual_asset(root, archive_path)
            visual_require(visual_sha256(archive_path) == row.get("sha256"),
                           "capture_hash_mismatch", str(archive_path))
        snapshot = record.get("semantic_snapshot")
        visual_require(isinstance(snapshot, dict) and isinstance(snapshot.get("retained_path"), str)
                       and visual_hex(snapshot.get("retained_file_sha256"), 64),
                       "missing_capture_snapshot", str(sidecar))
        snapshot_name = Path(snapshot.get("retained_path", "")).name
        visual_require(bool(snapshot_name) and path.with_name(snapshot_name).is_file(),
                       "missing_capture_snapshot", str(sidecar))
        snapshot_path = visual_path(root, snapshot_name, path.parent)
        visual_require(visual_sha256(snapshot_path) == snapshot.get("retained_file_sha256"),
                       "capture_hash_mismatch", f"{sidecar}: semantic snapshot")
        timestamp = capture.get("captured_at_unix_ms")
        visual_require(type(timestamp) is int and timestamp > 0, "missing_capture_date", str(sidecar))
        try:
            date = datetime.datetime.fromtimestamp(timestamp / 1000, datetime.timezone.utc).date().isoformat()
        except (ValueError, OverflowError, OSError):
            raise RuntimeError(f"tutorial visual missing_capture_date: {sidecar}") from None
        visual_require(date == graphic["capture_date"], "capture_date_mismatch", str(sidecar))
        return {"status": "historical_provenance_bound", "revision": record["source_revision"],
                "captured_chapter": record["chapter_id"], "captured_step": record["prose_step"]}
    # Older manual archives have several documented schemas. Preserve that
    # distinction; matching bytes and an exact producer are not a GUI verdict.
    evidence = path.parent / "evidence.json"
    if evidence.is_file():
        record = visual_json(evidence)
        matches = []

        def binds_asset(name):
            # Manual archives use repository-relative or sibling-relative paths.
            # A matching basename in a different directory is not a binding.
            return isinstance(name, str) and not Path(name).is_absolute() and (
                (root / name).resolve() == path or (path.parent / name).resolve() == path)

        def find_rows(value):
            if isinstance(value, dict):
                if binds_asset(value.get("path")) and "sha256" in value:
                    matches.append(value["sha256"])
                for field in ("raw_png", "context_svg"):
                    if binds_asset(value.get(field)):
                        matches.append(value.get(field + "_sha256"))
                for child in value.values():
                    find_rows(child)
            elif isinstance(value, list):
                for child in value:
                    find_rows(child)

        find_rows(record)
        visual_require(bool(matches) and all(digest == visual_sha256(path)
                                            for digest in matches), "capture_hash_mismatch", str(path))
        producer = record.get("display", record.get("producer", {}))
        visual_require(isinstance(producer, dict), "missing_capture_binding", str(evidence))
        revision = producer.get("revision", producer.get("source_revision", ""))
        visual_require(visual_hex(revision, 40),
                       "missing_capture_binding", str(evidence))
        visual_require(visual_hex(producer.get("gentle_gui_sha256"), 64),
                       "missing_capture_binding", str(evidence))
        if "capture_date" in record:
            visual_require(record["capture_date"] == graphic["capture_date"],
                           "capture_date_mismatch", str(evidence))
        return {"status": "historical_manual_bound", "revision": revision}
    visual_require(not path.name.endswith((".raw.png", ".context.svg", ".orientation.svg")),
                   "missing_capture_provenance", str(path))
    return {"status": "legacy_date_only_not_acceptance", "revision": None}


def check_tutorial_visuals(root):
    """Derive coverage from canonical source/catalog/manifest, not another registry."""
    root = root.resolve()
    catalog = visual_json(root / "docs/tutorial/catalog.json")
    manifest = visual_json(root / "docs/tutorial/manifest.json")
    visual_require(catalog.get("schema") == "gentle.tutorial_catalog.v2"
                   and manifest.get("schema") == "gentle.tutorial_manifest.v2",
                   "unsupported_catalog_schema", str(root))
    entries = {row["id"]: row for row in catalog["entries"]}
    chapters = {row["id"]: row for row in manifest["chapters"]}
    generated_pages = {}
    for page in sorted((root / "docs/tutorial/generated/chapters").glob("*.md")):
        text = page.read_text(encoding="utf-8")
        identifier = re.search(r'(?m)^chapter_id: "([^"]+)"$', text)
        if identifier:
            generated_pages[identifier[1]] = (page, text)
    rows, inspected = [], {}

    def inspect(path):
        if path not in inspected:
            inspected[path] = inspect_visual_asset(root, path)
        return inspected[path]

    for source_path in sorted((root / "docs/tutorial/sources").glob("*.json")):
        source = visual_json(source_path)
        if "id" not in source:
            continue
        identifier = source["id"]
        visual_require(isinstance(identifier, str), "invalid_metadata", str(source_path))
        visual_require(identifier in entries, "missing_catalog_entry", identifier)
        entry, chapter = entries[identifier], chapters.get(identifier)
        graphics = source.get("graphics", [])
        visual_require(isinstance(graphics, list), "invalid_metadata", identifier)
        visual_require(entry.get("graphics", []) == graphics
                       and (chapter is None or chapter.get("graphics", []) == graphics),
                       "metadata_projection_mismatch", identifier)
        page = (root / entry["path"]).resolve()
        visual_require(page.is_relative_to(root) and page.is_file(), "missing_page", identifier)
        text = page.read_text(encoding="utf-8")
        page_images = [visual_path(root, a or b, page.parent) for a, b in markdown_visuals(text)]
        for image in page_images:
            inspect(image)
        details = []
        for graphic in graphics:
            visual_require(isinstance(graphic, dict), "invalid_metadata", identifier)
            kind = graphic.get("kind")
            visual_require(kind in ("screenshot", "generated"), "invalid_kind", identifier)
            visual_require(isinstance(graphic.get("caption"), str) and bool(graphic["caption"].strip()),
                           "missing_caption", identifier)
            step = graphic.get("illustrates_step")
            visual_require(type(step) is int and step > 0, "invalid_step", identifier)
            path = visual_path(root, graphic.get("path"))
            asset = inspect(path)
            visual_require(path in page_images, "missing_placement", f"{identifier}: {graphic['path']}")
            if chapter is not None:
                visual_require(step <= len(chapter["gui_steps"]), "invalid_step", identifier)
                visual_require(identifier in generated_pages, "missing_generated_page", identifier)
                generated_page, generated_text = generated_pages[identifier]
                walkthrough = re.search(r"(?ms)^## (?:Walkthrough[^\n]*|GUI Steps)\n(.*?)(?=^## |\Z)", generated_text)
                visual_require(walkthrough is not None, "missing_walkthrough", identifier)
                body = walkthrough[1]
                all_matches = [visual_path(root, a or b, generated_page.parent)
                               for a, b in markdown_visuals(body)]
                visual_require(all_matches.count(path) == 1, "duplicate_step_placement",
                               f"{identifier}: {graphic['path']}")
                headings = list(re.finditer(r"(?m)^(?:### Step (\d+):|([0-9]+)\. )", body))
                segment = next((body[m.end():headings[i + 1].start() if i + 1 < len(headings) else len(body)]
                                for i, m in enumerate(headings) if int(m[1] or m[2]) == step), "")
                matches = [visual_path(root, a or b, generated_page.parent) for a, b in markdown_visuals(segment)]
                visual_require(matches.count(path) == 1, "wrong_step_placement",
                               f"{identifier}: step {step}: {graphic['path']}")
            if kind == "screenshot":
                try:
                    datetime.date.fromisoformat(graphic.get("capture_date", ""))
                except (ValueError, TypeError):
                    raise RuntimeError(f"tutorial visual missing_capture_date: {identifier}: {path}") from None
                visual_require(path.relative_to(root).as_posix().startswith("docs/screenshots/")
                               and (asset["format"] != "svg" or bool(asset["raster_images"])),
                               "figure_as_screenshot", f"{identifier}: {path}")
                provenance = screenshot_provenance(root, path, graphic)
            else:
                visual_require(not graphic.get("capture_date") and not path.relative_to(root).as_posix().startswith("docs/screenshots/"),
                               "screenshot_as_figure", f"{identifier}: {path}")
                provenance = {"status": "deterministic_figure_not_capture"}
            details.append({"kind": kind, "path": graphic["path"], "step": step, **provenance})
        category = ("typed_gui_acceptance" if chapter and chapter.get("gui_acceptance") else
                    "generated_workflow_companion" if chapter else
                    "manual_gui_walkthrough" if "gui" in entry.get("type", "") else "cli_reference")
        rows.append({"id": identifier, "category": category, "visuals": details})
    return {"tutorial_count": len(rows), "declared_visual_count": sum(len(row["visuals"]) for row in rows),
            "inspected_asset_count": len(inspected), "tutorials": rows,
            "non_claim": "Archive integrity and visual placement are not current-candidate GUI acceptance."}


def git(root, *args):
    return subprocess.run(
        ["git", "-C", str(root), "-c", f"core.hooksPath={Path(root).resolve() / '.disabled-hooks'}",
         "-c", "submodule.recurse=false", *args],
        check=True, capture_output=True,
    ).stdout


def prepare_checkout(source, destination, mode, attributes=None, revision=None):
    """Preserve Git history for review badges, without copying ignored caches."""
    name, autocrlf, eol = mode
    revision = revision or git(source, "rev-parse", "HEAD").decode().strip()
    git(source, "clone", "--shared", "--no-checkout", "--quiet",
        str(source.resolve()), str(destination))
    git(destination, "config", "core.autocrlf", autocrlf)
    git(destination, "config", "core.eol", eol)
    git(destination, "config", "core.longpaths", "true")
    git(destination, "update-ref", "--no-deref", "HEAD", revision)
    git(destination, "read-tree", revision)
    if attributes is not None:
        (destination / ".gitattributes").write_bytes(attributes)
        git(destination, "add", "--", ".gitattributes")
    # Export once, after the policy is in the index; otherwise Git may retain
    # already-written CRLF files whose stat cache still says they are unchanged.
    git(destination, "checkout-index", "--all", "--force")
    print(f"Tutorial checkout: {name}, revision={revision}, "
          f"attributes={'overlay' if attributes is not None else 'committed'}", flush=True)
    return revision


def check_checkout(binary, checkout, mode, timeout):
    check_tutorial_visuals(checkout)
    env = os.environ.copy()
    # This gate checks committed offline teaching artifacts, not live services.
    env.pop("GENTLE_TEST_ONLINE", None)
    for argument in ("parity-matrix-check", "tutorial-catalog-check", "--check", "tutorial-check"):
        print(f"[{mode}] {binary.name} {argument}", flush=True)
        try:
            subprocess.run([str(binary), argument], cwd=checkout, env=env,
                           check=True, timeout=timeout)
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as error:
            raise RuntimeError(f"{mode} checkout failed: {argument}: {error}") from error


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=ROOT)
    parser.add_argument("--binary", type=Path)
    parser.add_argument("--visuals-only", action="store_true",
                        help="Check current source, visual links and archive metadata without a binary")
    parser.add_argument("--timeout-seconds", type=float, default=900)
    parser.add_argument("--attributes-from", type=Path,
                        help="Explicit local .gitattributes overlay; not exact-commit acceptance")
    args = parser.parse_args(argv)
    if args.visuals_only:
        print(json.dumps(check_tutorial_visuals(args.repo_root), indent=2, sort_keys=True))
        return 0
    if args.binary is None:
        parser.error("--binary is required unless --visuals-only is selected")
    binary = args.binary.resolve()
    if not binary.is_file():
        parser.error(f"Build gentle_examples_docs first: {binary} does not exist")
    if not math.isfinite(args.timeout_seconds) or args.timeout_seconds <= 0:
        parser.error("--timeout-seconds must be positive")
    attributes = args.attributes_from.read_bytes() if args.attributes_from else None
    revision = git(args.repo_root, "rev-parse", "HEAD").decode().strip()
    for mode in MODES:
        # Release each checkout before creating the next one on small machines.
        with tempfile.TemporaryDirectory(prefix=f"gentle-tutorial-{mode[0]}-",
                                         dir=os.environ.get("RUNNER_TEMP")) as tmp:
            checkout = Path(tmp) / "repo"
            prepare_checkout(args.repo_root, checkout, mode, attributes, revision)
            check_checkout(binary, checkout, mode[0], args.timeout_seconds)
    print("Tutorial checkout portability: LF and CRLF passed", flush=True)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except subprocess.CalledProcessError as error:
        detail = error.stderr.decode(errors="replace") if error.stderr else ""
        raise SystemExit(f"{error}\n{detail}") from error
    except (OSError, RuntimeError) as error:
        raise SystemExit(str(error)) from error
