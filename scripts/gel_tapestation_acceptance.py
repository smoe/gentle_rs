#!/usr/bin/env python3
"""Replay the public 26R266 gel fixtures through an explicitly supplied GENtle.

No build, network, band detector or sizing implementation lives here. Coordinates
are the annotation-assisted marks retained in the fixture. Only GENtle computes
sizes; this script compares its results with separately retained vendor values.
Optional priming comparisons also use only GENtle's saved display profiles.
Outputs are written to a new directory, never over an existing audit.
"""

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import re
import subprocess
import tempfile


ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "test_files/fixtures/gel_images/tapestation_26r266"
ALGORITHM_FILES = (
    "src/gel_image.rs", "src/engine/io/gel_images.rs",
    "crates/gentle-protocol/src/gel_image.rs", "crates/gentle-render/src/gel_image.rs",
)
PRIMING_FILES = (*ALGORITHM_FILES, "src/gel_workbench.rs",
                 "src/engine/io/gel_workbench.rs",
                 "crates/gentle-protocol/src/gel_workbench.rs",
                 "crates/gentle-render/src/gel_workbench.rs",
                 "crates/gentle-render/src/gel_workbench/comparison_display.rs")


def sha256(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def fixture_file(root, name):
    if not isinstance(name, str) or not re.fullmatch(r"[A-Za-z0-9_.-]+", name) or name in (".", ".."):
        raise ValueError(f"Unsafe fixture filename: {name!r}")
    path = root / name
    if path.is_symlink() or not path.is_file() or not path.resolve().is_relative_to(root.resolve()):
        raise ValueError(f"Missing or nonlocal fixture: {name}")
    return path


def verify_source_files(binary_source, root=ROOT, files=ALGORITHM_FILES):
    """Compare canonical Git blobs, not CRLF-converted working-tree bytes."""
    root = Path(root)
    subprocess.run(["git", "diff", "--quiet", "HEAD", "--", *files],
                   cwd=root, capture_output=True, timeout=30, check=True)
    matching = {}
    for name in files:
        def blob(revision):
            return subprocess.run(["git", "show", f"{revision}:{name}"], cwd=root,
                                  capture_output=True, timeout=30, check=True).stdout
        compiled = blob(binary_source)
        if compiled != blob("HEAD"):
            raise ValueError(f"Supplied CLI's gel source differs from the current source: {name}")
        matching[name] = hashlib.sha256(compiled).hexdigest()
    return matching


def load_fixture(root=FIXTURE):
    root = Path(root)
    manifest = json.loads(fixture_file(root, "manifest.json").read_bytes())
    if manifest.get("schema") != "gentle.tapestation_gel_regression.v1":
        raise ValueError("Unsupported TapeStation fixture")
    for group in ("sources", "cases"):
        ids = [row.get("id") for row in manifest[group]]
        if any(not isinstance(identifier, str) or not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]{0,79}", identifier)
               for identifier in ids) or len(set(ids)) != len(ids):
            raise ValueError(f"Invalid or duplicate {group} identity; filenames must stay inside the audit directory")
    for row in manifest["sources"]:
        for kind in ("pdf", "image"):
            if sha256(fixture_file(root, row[kind])) != row[f"{kind}_sha256"]:
                raise ValueError(f"{row['id']} {kind} digest mismatch")
    return manifest


def request_for_case(manifest, case):
    """Project retained geometry into the existing engine contract, without sizing."""
    source = next(row for row in manifest["sources"] if row["id"] == case["source"])
    lanes = {row["id"]: row for row in source["lanes"]}
    selected = [case["calibration_lane"], *case["sample_lanes"]]
    return {
        "report_id": case["id"], "image_id": source["id"],
        "image_sha256": source["image_sha256"], "size_kind": "linear_dna_bp",
        "migration": "down",
        "lanes": [{"id": identifier, "label": identifier,
                   "min": {"x": lanes[identifier]["center_x"] - 31.5, "y": 100.0},
                   "max": {"x": lanes[identifier]["center_x"] + 31.5, "y": 640.0}}
                  for identifier in selected],
        "ladder": {
            "lane_id": case["calibration_lane"], "label": "D5000 100-5000 bp",
            "source": case["ladder_source"], "prestained": False,
            "bands": [{"id": f"ref_{size}", "size": size,
                       "center": {"x": lanes[case["calibration_lane"]]["center_x"], "y": y}}
                      for size, y in case["references"]],
        },
        "sample_bands": [{"id": band["id"], "lane_id": band["lane"],
                          "label": band["id"],
                          "center": {"x": lanes[band["lane"]]["center_x"], "y": band["y"]},
                          "position_half_width_px": 1.0}
                         for band in case["bands"]],
    }


def compare_report(report, case, tolerance):
    estimates = report["estimates"]
    by_id = {row["band_id"]: row for row in estimates}
    if len(by_id) != len(estimates) or set(by_id) != {row["id"] for row in case["bands"]}:
        raise ValueError("GENtle omitted, duplicated or invented an estimate")
    rows = []
    sizes = [row[0] for row in case["references"]]
    for band in case["bands"]:
        measured = by_id[band["id"]]
        if measured["lane_id"] != band["lane"] or measured["center"]["y"] != band["y"]:
            raise ValueError("GENtle changed the retained lane or pixel coordinate")
        value = measured["estimated_size"]
        row = {"band_id": band["id"], "lane": band["lane"],
               "role": band.get("role", "sample_peak"),
               "estimated_bp": value, "vendor_peak_bp": band.get("vendor_peak_bp"),
               "relative_error": None, "agreement_pass": None, "check_pass": True}
        if "expected_status" in band:
            row["comparison"] = "outside_calibration_control"
            row["check_pass"] = (measured["status"] == band["expected_status"] and value is None
                                 and not measured["reference_band_ids"])
        else:
            valid = (measured["status"] == "interpolated" and isinstance(value, (int, float))
                     and not isinstance(value, bool) and math.isfinite(value)
                     and min(sizes) <= value <= max(sizes)
                     and len(measured["reference_band_ids"]) == 2)
            expected = band.get("vendor_peak_bp")
            if expected is None:
                row["comparison"] = "vendor_peak_unavailable"
                row["check_pass"] = valid
            else:
                row["comparison"] = "vendor_peak" if row["role"] == "sample_peak" else "cross_ladder"
                row["relative_error"] = abs(value - expected) / expected if valid else None
                row["agreement_pass"] = valid and row["relative_error"] <= tolerance
                row["check_pass"] = row["agreement_pass"]
        rows.append(row)
    return rows


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, ensure_ascii=True, allow_nan=False) + "\n",
                    encoding="utf-8", newline="\n")


def priming_requests(manifest, fixture=FIXTURE):
    """Compose explicit requests; no profile, density or statistic is calculated."""
    fixture = Path(fixture)
    path = fixture_file(fixture, "region_summaries.json")
    sums = dict(line.split("  ", 1)[::-1] for line in
                fixture_file(fixture, "SHA256SUMS").read_text(encoding="utf-8").splitlines())
    if sha256(path) != sums.get(path.name):
        raise ValueError("Priming pairing source digest mismatch")
    data = json.loads(path.read_bytes())
    source = next(row for row in manifest["sources"] if row["id"] == "after")
    if data["source_pdf_sha256"]["after"] != source["pdf_sha256"]:
        raise ValueError("Priming groups refer to a different source PDF")
    rows = {(r["source"], r["lane"]): r for r in data["rows"]}
    lanes = {r["id"]: r for r in source["lanes"]}
    groups = data["after_priming_groups"]
    if [g["cell_line"] for g in groups] != ["cell_1", "cell_2"]:
        raise ValueError("Expected the two explicit owner-associated cell groups")
    plan = []
    for group in groups:
        cell = group["cell_line"]
        profiles = []
        for role, priming in (("custom", "custom"), ("manufacturer", "manufacturer_crta"),
                              ("combined", "custom_and_manufacturer")):
            lane = group[f"{role}_lane"]
            if not isinstance(lane, str) or not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]{0,79}", lane):
                raise ValueError("Invalid priming lane identity for retained request filename")
            if (rows["after", lane]["cell_line"], rows["after", lane]["priming"]) != (cell, priming):
                raise ValueError("Priming lane role/cell differs from its retained association")
            if lanes[lane]["role"] != "sample":
                raise ValueError("A ladder cannot be a priming sample")
            cases = [c for c in manifest["cases"] if c["source"] == "after"
                     and lane in c["sample_lanes"] and c["calibration_lane"] != lane]
            if len(cases) != 1:
                raise ValueError("Priming sample needs one explicit retained calibration")
            case = cases[0]
            profiles.append({"role": role, "lane": lane, "ladder": case["calibration_lane"],
                             "request": {"action": "extract_profile",
                                         "profile_id": f"after_{cell}_{role}_{lane}",
                                         "analysis_id": case["id"], "lane_id": lane,
                                         "run_id": "2026-10-13_after_cleanup",
                                         "dark_bands": True, "background": 0.0,
                                         "masked_ranges": []}})
        by_role = {p["role"]: p["request"]["profile_id"] for p in profiles}
        comparisons = [{"action": "compare", "comparison_id": f"{cell}_{left}_vs_{right}",
                        "left_profile": by_role[left], "right_profile": by_role[right],
                        "range": [100.0, 5000.0], "short_threshold": None}
                       for left, right in (("custom", "manufacturer"), ("custom", "combined"),
                                           ("manufacturer", "combined"))]
        plan.append({"cell_line": cell, "profiles": profiles, "comparisons": comparisons})
    if len({p["lane"] for g in plan for p in g["profiles"]}) != 6:
        raise ValueError("Priming comparisons must retain six distinct sample lanes")
    return plan


def run_priming_comparisons(command, manifest, output, fixture):
    plan = priming_requests(manifest, fixture)
    summary = {"schema": "gentle.tapestation_priming_replay.v1",
               "pairing_source_sha256": sha256(fixture / "region_summaries.json"),
               "range_bp": [100.0, 5000.0], "background": 0.0, "masked_ranges": [],
               "evidence": "image_display; dark bands; unscaled signal per bp, not detector fluorescence",
               "nonclaims": ["No subtraction, fitted coefficients, marker/background fitting or gene assignment",
                             "CDF/QQ/KS-distance describe display-weighted shape, not fragment counts or significance",
                             "15/10000 bp alignment markers lie outside the chosen 100-5000 bp window",
                             "L1/L2 calibration, expiry warnings and vendor display processing remain limitations"],
               "cells": []}

    def request(name, value):
        path = output / f"{name}.request.json"
        write_json(path, value)
        return command("gel-image", "workbench", "@" + str(path))["result"]["gel_workbench"]

    for group in plan:
        cell = group["cell_line"]
        for p in group["profiles"]:
            request(p["request"]["profile_id"], p["request"])
        for comparison in group["comparisons"]:
            request(comparison["comparison_id"], comparison)
        ids = [r["comparison_id"] for r in group["comparisons"]]
        base = {"action": "export", "analysis_ids": [], "comparison_ids": ids,
                "title": f"26R266 {cell}: three priming conditions; no subtraction"}
        # Separate CLI invocations reload stored reports; export does not re-extract them.
        for extension in ("json", "pdf"):
            export = dict(base, path=str(output / f"{cell}-priming.{extension}"))
            if extension == "pdf":
                export["comparison_display"] = "signal_density_per_size"
            request(f"{cell}-priming-{extension}", export)
        saved = json.loads((output / f"{cell}-priming.json").read_bytes())["comparisons"]
        expected = [{k: v for k, v in r.items() if k != "action"} for r in group["comparisons"]]
        if [r["request"] for r in saved] != expected:
            raise ValueError("Saved comparison request/order changed")
        profiles = {p["request"]["profile_id"]: p for p in group["profiles"]}
        for report in saved:
            for p in (report["left"], report["right"]):
                if p["evidence_kind"] != "image_display" or p["lane_id"] != profiles[p["id"]]["lane"]:
                    raise ValueError("GENtle changed priming profile identity/evidence")
            if report["p_value"] is not None:
                raise ValueError("Display comparisons must not manufacture a significance p-value")
        summary["cells"].append({"cell_line": cell,
                                 "profiles": [{k: v for k, v in p.items() if k != "request"}
                                              | {"profile_id": p["request"]["profile_id"]}
                                              for p in group["profiles"]],
                                 "comparison_ids": ids,
                                 "reports": f"{cell}-priming.json",
                                 "descriptive_results": [{key: r[key] for key in (
                                     "request", "left_summary", "right_summary", "ks_distance",
                                     "wasserstein_distance", "p_value", "excluded_intervals")}
                                     for r in saved]})
    return summary


def run_audit(binary, output, fixture=FIXTURE, compare_priming=False):
    binary, output, fixture = Path(binary).resolve(), Path(output), Path(fixture).resolve()
    if not binary.is_file() or not os.access(binary, os.X_OK):
        raise ValueError(f"Supply an existing executable GENtle CLI: {binary}")
    manifest = load_fixture(fixture)
    output.mkdir(parents=True, exist_ok=False)
    output = output.resolve()
    binary_digest = sha256(binary)
    comparisons = []
    reports = {}
    with tempfile.TemporaryDirectory(prefix="gentle gel audit's ") as temporary:
        work = Path(temporary)
        env = {key: value for key, value in os.environ.items() if not key.startswith("GENTLE_")}
        for name, keys in (
            ("home", ("HOME", "USERPROFILE")), ("config", ("XDG_CONFIG_HOME",)),
            ("cache", ("XDG_CACHE_HOME",)), ("data", ("XDG_DATA_HOME",)),
            ("runtime", ("XDG_RUNTIME_DIR",)), ("tmp", ("TMPDIR", "TMP", "TEMP")),
            ("system", ("GENTLE_SYSTEM_CONFIG_ROOT",)), ("project", ("GENTLE_PROJECT_ROOT",)),
        ):
            directory = work / name
            directory.mkdir()
            env.update({key: str(directory) for key in keys})
        env["GENTLE_ASSET_ROOT"] = str(ROOT)

        def invoke(*args):
            completed = subprocess.run([str(binary), *args], cwd=work, env=env,
                                       text=True, capture_output=True, timeout=120)
            if completed.returncode:
                raise RuntimeError(f"GENtle failed ({completed.returncode}): {completed.stderr}")
            return completed.stdout

        version = invoke("--version")
        source = re.search(r"\+git\.([0-9a-f]{40})(?:\s|$)", version)
        if not source:
            raise ValueError("CLI does not identify a full source SHA")
        matching_files = verify_source_files(source[1], files=PRIMING_FILES if compare_priming else ALGORITHM_FILES)
        state = output / "project.gentle.json"

        def command(*args):
            return json.loads(invoke("--state", str(state), "--full-report", *args))

        for row in manifest["sources"]:
            result = command("gel-image", "import", json.dumps({
                "image_id": row["id"], "path": str(fixture / row["image"])}))
            if result["result"]["gel_image"]["sha256"] != row["image_sha256"]:
                raise ValueError("Imported source identity changed")

        for case in manifest["cases"]:
            request = request_for_case(manifest, case)
            request_path = output / f"{case['id']}.request.json"
            write_json(request_path, request)
            result = command("gel-image", "analyze", "@" + str(request_path))
            report = result["result"]["gel_image_analysis"]
            # A separate process reloads the project through the shared Shell.
            inspected = command("shell", f"gel-image inspect {case['id']}")
            if inspected["result"]["gel_image_analysis"] != report:
                raise ValueError("Persisted/shared-shell report differs")
            reports[case["id"]] = report
            comparisons.extend(dict(case=case["id"], **row) for row in compare_report(
                report, case, manifest["comparison_relative_tolerance"]))
            for format_name in ("json", "tsv", "svg"):
                command("gel-image", "export", json.dumps({
                    "report_id": case["id"], "format": format_name,
                    "path": str(output / f"{case['id']}.{format_name}")}))

        priming = run_priming_comparisons(command, manifest, output, fixture) if compare_priming else None
        if priming is not None:
            priming.update({"binary_source_sha": source[1], "binary_sha256": binary_digest,
                            "source_image_sha256": next(r["image_sha256"] for r in manifest["sources"]
                                                        if r["id"] == "after"),
                            "matching_gel_source_files": matching_files})
            write_json(output / "priming-summary.json", priming)
        if sha256(binary) != binary_digest:
            raise ValueError("CLI binary changed during replay; discard this audit")
        write_json(output / "comparison.json", comparisons)
        receipt = {
            "schema": "gentle.tapestation_gel_audit.v1", "binary_sha256": binary_digest,
            "binary_version": version, "binary_source_sha": source[1],
            "profile_evidence": "Not embedded by supplied binary; no timing/profile acceptance claimed",
            "fixture_manifest_sha256": sha256(fixture / "manifest.json"),
            "matching_gel_source_files": matching_files,
            "source_comparison_basis": "Canonical Git blobs with clean tracked implementation files; independent of checkout line endings",
            "source_context": manifest["sources"], "nonclaims": manifest["nonclaims"],
            "all_required_checks_pass": all(row["check_pass"] for row in comparisons),
            "sample_peak_comparisons": sum(row["comparison"] == "vendor_peak" for row in comparisons),
            "cross_ladder_comparisons": sum(row["comparison"] == "cross_ladder" for row in comparisons),
            "vendor_peak_unavailable": sum(row["comparison"] == "vendor_peak_unavailable" for row in comparisons),
            "outside_calibration_controls": sum(row["comparison"] == "outside_calibration_control" for row in comparisons),
            "artifacts": {path.name: sha256(path) for path in sorted(output.iterdir()) if path.is_file()},
        }
        if priming is not None:
            receipt["priming_comparisons"] = priming
        write_json(output / "receipt.json", receipt)
    return receipt, comparisons, reports


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gentle-cli", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--fixture-dir", type=Path, default=FIXTURE)
    parser.add_argument("--compare-priming", action="store_true",
                        help="Also compare all three after-cleanup priming conditions per cell, without subtraction")
    args = parser.parse_args()
    receipt, _, _ = run_audit(args.gentle_cli, args.output_dir, args.fixture_dir, args.compare_priming)
    print(json.dumps(receipt, indent=2))
    return 0 if receipt["all_required_checks_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
