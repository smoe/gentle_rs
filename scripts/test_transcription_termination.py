"""Offline, actual-process termination smoke through Python and direct CLI.

All DNA/catalog data are hand-crafted synthetic literals below, recreated in a
private temporary directory. They exercise GENtle transport/provenance/annotation
boundaries and do not represent validated biological elements. No network or
external folding program is used; no repository project is read or modified.

Run after building gentle_cli:
  python3 scripts/test_transcription_termination.py --cli target/debug/gentle_cli
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "integrations" / "python"))
from gentle_py import GentleClient  # noqa: E402


def digest(raw: bytes) -> str:
    return "sha256:" + hashlib.sha256(raw).hexdigest()


def run(cli: Path) -> None:
    with tempfile.TemporaryDirectory(prefix="gentle-termination-") as scratch:
        directory = Path(scratch) / "O'Brien inputs"
        directory.mkdir()
        state = directory / "project.json"
        client = GentleClient(cli_cmd=[str(cli)], state_path=str(state), cwd=str(directory))
        element = "GCGCGCAAAAGCGCGCTTTTTTTT"
        client.op({"CreateSequenceFromText": {"sequence_text": "CCC" + element + "CCC",
                   "output_id": "toy", "name": "Synthetic transport fixture", "circular": False}})
        catalog = {
            "schema": "gentle.transcription_termination_catalog.v1", "catalog_id": "synthetic", "version": "1",
            "entries": [{"id": "invented", "label": "Synthetic element", "sequence": element,
                         "sequence_sha256": digest(element.encode()), "host": "bacteria",
                         "mechanism": "bacterial_intrinsic", "element_kind": "polymerase_termination_region",
                         "source_url": "synthetic:scripts/test_transcription_termination.py", "source_version": "1",
                         "recreation": "Literal synthetic element in this script", "redistribution_status": "synthetic",
                         "approval_provenance": "Hand-crafted transport fixture; no scientific approval"}],
        }
        catalog_path = directory / "catalog.json"
        catalog_bytes = (json.dumps(catalog, indent=2) + "\n").encode()
        catalog_path.write_bytes(catalog_bytes)
        request = {"seq_id": "toy", "host": "bacteria", "polymerase": "bacterial_rnap",
                   "catalog": {"path": str(catalog_path), "expected_sha256": digest(catalog_bytes)}}
        before = state.read_bytes()
        report = client.op({"ScreenTranscriptionTermination": {"request": request, "path": None}})["transcription_termination_screen"]
        assert state.read_bytes() == before, "Python inspection changed project bytes"
        assert report["status"] == "completed" and report["rows"]
        request_path = directory / "request.json"
        request_path.write_text(json.dumps(request), encoding="utf-8")
        exported = directory / "saved report.json"
        process = subprocess.run([str(cli), "--state", str(state), "terminators", "screen", "@" + str(request_path),
                                  "--output", str(exported)], capture_output=True, text=True, check=True)
        assert process.stdout.strip(), "Direct CLI did not return a report"
        assert json.loads(exported.read_bytes()) == report, "Python and direct CLI reports differ"
        assert state.read_bytes() == before, "Direct CLI inspection changed project bytes"
        selected = {"screen": request, "expected_report_sha256": report["content_sha256"],
                    "row_ids": [report["rows"][0]["row_id"]]}
        result = client.op({"MaterializeTranscriptionTerminationFeatures": {"request": selected}})
        assert result["changed_seq_ids"] == ["toy"]
        after = state.read_bytes()
        fresh = client.op({"ScreenTranscriptionTermination": {"request": request, "path": None}})["transcription_termination_screen"]
        selected["expected_report_sha256"] = fresh["content_sha256"]
        result = client.op({"MaterializeTranscriptionTerminationFeatures": {"request": selected}})
        assert not result["changed_seq_ids"], "Fresh reapply duplicated an annotation"
        annotated = directory / "round trip.gb"
        client.op({"SaveFile": {"seq_id": "toy", "path": str(annotated), "format": "GenBank"}})
        client.op({"LoadFile": {"path": str(annotated), "as_id": "imported_toy"}})
        imported_request = dict(request, seq_id="imported_toy")
        imported = client.op({"ScreenTranscriptionTermination": {"request": imported_request, "path": None}})["transcription_termination_screen"]
        imported_row = next(row for row in imported["rows"] if row["evidence_kind"] == "reference_match")
        imported_selection = {"screen": imported_request, "expected_report_sha256": imported["content_sha256"],
                              "row_ids": [imported_row["row_id"]]}
        reapplied = client.op({"MaterializeTranscriptionTerminationFeatures": {"request": imported_selection}})
        assert not reapplied["changed_seq_ids"], "GenBank round-trip lost portable evidence identity"
        # Rebinding rejects changed raw bytes even when JSON meaning is identical.
        catalog_path.write_bytes(catalog_bytes.replace(b"\n", b"\r\n"))
        refused = client.op({"ScreenTranscriptionTermination": {"request": request, "path": None}})["transcription_termination_screen"]
        assert refused["status"] == "freshness_refused"
        # A later operation may persist its journal; compare inspection to the
        # current project, not to the pre-apply file.
        final = state.read_bytes()
        client.op({"ScreenTranscriptionTermination": {"request": request, "path": None}})
        assert state.read_bytes() == final and after != before
        print("PASS: Python/direct CLI report parity, export, project preservation, confirmed apply, idempotence, GenBank round-trip and raw-byte freshness")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cli", type=Path, required=True)
    run(parser.parse_args().cli.resolve())
