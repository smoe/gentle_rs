#!/usr/bin/env python3
"""Prepare public synthetic PATZ1 study-inspection inputs using a supplied GENtle CLI.

No build, download, inner agent, BLAST, order or approved-study execution.
The three comparison panels are the unchanged scientific operations from 04.06;
only output paths change. The study uses 04.07's separate, stricter request.
"""

import argparse
import copy
import hashlib
import json
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[1]
WORKFLOW = ROOT / "docs/examples/workflows/patz1_endpoint_sybr_transcript_assay_panel_offline.json"
STUDY = ROOT / "docs/tutorial/inputs/transcript_assay_followup_study.json"


def digest(path):
    return "sha256:" + hashlib.sha256(path.read_bytes()).hexdigest()


def prepare(binary, output):
    binary = Path(binary).resolve(strict=True)
    output = Path(output).resolve()
    output.mkdir()  # A previous acceptance must never be overwritten.
    state = output / "patz1-study.gentle.json"
    calls = []

    def cli(*args):
        command = [str(binary), "--state", str(state), *map(str, args)]
        run = subprocess.run(command, cwd=ROOT, text=True, capture_output=True,
                             timeout=180, check=False)
        ordinal = len(calls) + 1
        stdout = output / f"call-{ordinal:02d}.stdout.json"
        stderr = output / f"call-{ordinal:02d}.stderr.txt"
        stdout.write_text(run.stdout, encoding="utf-8")
        stderr.write_text(run.stderr, encoding="utf-8")
        calls.append({"command": command, "exit_code": run.returncode,
                      "stdout": stdout.name, "stderr": stderr.name})
        if run.returncode:
            raise RuntimeError(f"GENtle call {ordinal} failed; inspect {stderr}")

    workflow = json.loads(WORKFLOW.read_text(encoding="utf-8"))
    operations = copy.deepcopy(workflow["workflow"]["ops"][:5])
    for ordinal, operation in enumerate(operations):
        if "DesignTranscriptAssayPanel" in operation:
            spec = operation["DesignTranscriptAssayPanel"]
            spec["path"] = str(output / f'{spec["report_id"]}.json')
        path = output / f"operation-{ordinal + 1:02d}.json"
        path.write_text(json.dumps(operation, indent=2), encoding="utf-8")
        cli("op", "@" + str(path))

    # Preserve the source request and its input bindings; do not turn the
    # successful relaxed panel recipe into an allegedly approved study result.
    request = json.loads(STUDY.read_text(encoding="utf-8"))
    request["isoform_evidence_path"] = str(ROOT / request["isoform_evidence_path"])
    for item in request["junction_evidence"]:
        item["path"] = str(ROOT / item["path"])
    request_path = output / "study.request.json"
    request_path.write_text(json.dumps(request, indent=2), encoding="utf-8")
    plan_path = output / "study.plan.json"
    cli("primers", "plan-gene-isoform-study", "@" + str(request_path),
        "--normalized-request", output / "study.normalized.json",
        "--path", plan_path, "--workflow", output / "study.workflow.json")
    plan = json.loads(plan_path.read_text(encoding="utf-8"))
    publication = {
        "schema": "gentle.gene_isoform_assay_publication_request.v1",
        "report_id": "patz1_gui_study_teaching",
        "title": "Synthetic PATZ1 planning review: assays not yet executed",
        "genes": [{
            "gene_symbol": plan["gene_symbol"],
            "study_plan": {"path": plan_path.name, "expected_sha256": digest(plan_path)},
            "status": "pending",
            "status_reason": "This study workflow has not been executed. The project also contains separate 04.06 comparison panels, not results of this study. No order approval or experimental validation.",
        }],
    }
    publication_path = output / "publication.request.json"
    publication_path.write_text(json.dumps(publication, indent=2), encoding="utf-8")
    cli("primers", "publish-gene-isoform-study", publication_path, output / "cli-dossier")
    inputs = [WORKFLOW, STUDY, ROOT / workflow["required_files"][0]]
    inputs.extend(ROOT / path for path in workflow["required_files"][1:])
    receipt = {
        "schema": "gentle.tutorial_gene_assay_study_preparation.v1",
        "synthetic": True, "study_executed": False, "gui_accepted": False,
        "binary": str(binary), "binary_sha256": digest(binary),
        "inputs": [{"path": str(p.relative_to(ROOT)), "sha256": digest(p)} for p in inputs],
        "calls": calls,
        "artifacts": [{"path": str(p.relative_to(output)), "sha256": digest(p)}
                      for p in sorted(output.rglob("*")) if p.is_file()],
    }
    (output / "preparation-receipt.json").write_text(json.dumps(receipt, indent=2), encoding="utf-8")
    return output


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gentle-cli", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path,
                        help="New directory; parent must exist")
    options = parser.parse_args()
    print(prepare(options.gentle_cli, options.output_dir))
