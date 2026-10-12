#!/usr/bin/env python3
"""Offline tutorial orchestration; scientific outputs come only from GENtle.

Uses the checked-in synthetic workflow and evidence cases, retains both exact
operation payloads, and asserts report facts. No shell, sequence scoring, or
post-hoc editing of GENtle reports occurs here. BLAST+ must be installed.
"""

import argparse
import hashlib
import json
from pathlib import Path
import subprocess


ROOT = Path(__file__).resolve().parents[2]
ASSETS = ROOT / "docs/examples/assets/region_homology_demo"


def digest(path):
    with path.open("rb") as handle:
        return "sha256:" + hashlib.file_digest(handle, "sha256").hexdigest()


def write_json(path, payload):
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def require(condition, description):
    if not condition:
        raise RuntimeError(description)


def run(gentle, output):
    gentle = gentle.resolve(strict=True)
    output = output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    require(not any(output.iterdir()), "Use an empty output directory; old evidence must not be overwritten")
    state = output / "tutorial.project.json"
    receipts = []

    def execute(workflow, name, state_path):
        path = output / (name + ".workflow.json")
        write_json(path, workflow)
        command = [str(gentle), "--state", str(state_path), "workflow", "@" + str(path)]
        completed = subprocess.run(command, cwd=ROOT, capture_output=True, timeout=180)
        stdout, stderr = output / (name + ".stdout"), output / (name + ".stderr")
        stdout.write_bytes(completed.stdout)
        stderr.write_bytes(completed.stderr)
        receipts.append({"command": command, "exit_code": completed.returncode,
                         "workflow_sha256": digest(path), "stdout_sha256": digest(stdout),
                         "stderr_sha256": digest(stderr)})
        require(completed.returncode == 0, f"GENtle failed; inspect {stderr} and {stdout}")

    example = json.loads((ROOT / "docs/examples/workflows/region_homology_promoter_modules_offline.json").read_text())
    workflow = example["workflow"]
    for operation in workflow["ops"]:
        _, payload = next(iter(operation.items()))
        request = payload.get("request", payload)
        if "catalog_path" in request:
            request["catalog_path"] = str(ASSETS / "genomes.json")
        if "cache_dir" in request:
            request["cache_dir"] = str(output / "genome_cache")
        if "ScreenGenomicRegionHomology" in operation:
            payload["path"] = str(output / "homology_report.json")
            write_json(output / "homology_request.json", request)
    execute(workflow, "screen", state)
    report = json.loads((output / "homology_report.json").read_text())
    require(report["schema"] == "gentle.genomic_region_homology_screen.v1", "Wrong screen schema")
    require(report["omitted_insertions"], "Expected synthetic target insertion is missing")
    require(all(len(row["query_projection"]) == len(report["query"]["sequence"]) for row in report["alignment_rows"]), "Target insertion changed displayed query width")
    require(report["same_genome_nonself_locus_count"] > 0, "Synthetic competing copy was lost")
    cases_path = ASSETS / "module_cases.json"
    cases = json.loads(cases_path.read_text())
    operations = [{"RenderGenomicRegionHomologySvg": {"report": report, "path": str(output / "homology.svg")}}]
    for case in cases["cases"]:
        spans = [dict(cases["evidence_spans"][key], source_sha256=digest(cases_path)) for key in case["spans"]]
        request = {"homology_report": report, "selected_evidence_spans": spans,
                   "max_partner_gap_bp": 100, "max_partner_gap_difference_bp": 0,
                   "max_same_genome_query_coverage_percent": case.get("max_same_genome_query_coverage_percent", 80.0)}
        write_json(output / (case["id"] + ".request.json"), request)
        operations.append({"AssessPromoterConservedModules": {"request": request, "path": str(output / (case["id"] + ".json"))}})
    state_before = digest(state)
    # These operations consume self-contained reports, not the source project.
    execute({"run_id": "synthetic_module_assessments", "ops": operations},
            "assess_and_render", output / "report_only.project.json")
    require(digest(state) == state_before, "Read-only assessment/rendering changed the project")
    facts = []
    for case in cases["cases"]:
        result = json.loads((output / (case["id"] + ".json")).read_text())
        require(result["hypothesis"] == case["expected_hypothesis"], f"Unexpected decision for {case['id']}: {result['hypothesis']}")
        require(result["homology_report_sha256"] == report["content_sha256"], "Assessment bound another homology report")
        require(result["non_claims"], "Missing scientific interpretation boundary")
        rules = {rule["rule_id"]: rule for rule in result["decision_trace"]}
        require(rules["same_genome_evidence_available"]["satisfied"], "Same-genome comparison was not assessed")
        if case["id"] == "paired":
            require(any(context["passed"] for context in result["partner_contexts"]), "No target-coordinate evidence for the paired result")
        facts.append({"case": case["id"], "hypothesis": result["hypothesis"], "report_sha256": result["content_sha256"], "rules": rules})
    require((output / "homology.svg").stat().st_size > 100, "Empty SVG")
    write_json(output / "tutorial_receipt.json", {
        "schema": "tutorial_harness.conservation_acceptance.v1",
        "producer": "Repository tutorial harness; not a GENtle scientific report",
        "scientific_authority": "GENtle reports",
        "binary_sha256": digest(gentle), "input_cases_sha256": digest(cases_path),
        "commands": receipts, "verified_cases": facts,
        "artifacts": {path.name: digest(path) for path in sorted(output.iterdir()) if path.is_file()},
    })
    print(f"PASS: four GENtle module decisions, fixed-width alignment, insertion provenance, SVG, unchanged state.\n{output / 'tutorial_receipt.json'}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gentle", type=Path, default=ROOT / "target/debug/gentle_cli")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    run(args.gentle, args.output)
