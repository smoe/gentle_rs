#!/usr/bin/env python3
"""Prepare authentic PATZ1 annotation and primer-design inputs, offline.

GENtle owns sequence import, transcript construction, comparison and study
planning. This helper verifies pinned files and submits shared operations. It
does not execute primer design, BLAST, an inner agent or an order.
"""

import argparse
import hashlib
import json
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "test_files/fixtures/transcript_assay_panel/patz1_reference"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")


def prepare(binary, output):
    binary = Path(binary).resolve(strict=True)
    manifest = json.loads((FIXTURE / "manifest.json").read_text())
    for name, expected in manifest["files"].items():
        if sha(FIXTURE / name) != expected:
            raise ValueError(f"Pinned PATZ1 input hash mismatch: {name}")
    output = Path(output).resolve()
    output.mkdir()  # Never replace an earlier review or acceptance.
    state = output / "patz1.gentle.json"
    calls = []

    def cli(*args):
        command = [str(binary), "--state", str(state), *map(str, args)]
        run = subprocess.run(command, cwd=ROOT, text=True, capture_output=True,
                             timeout=300, check=False)
        ordinal = len(calls) + 1
        stdout = output / f"call-{ordinal:02d}.stdout.json"
        stderr = output / f"call-{ordinal:02d}.stderr.txt"
        stdout.write_text(run.stdout, encoding="utf-8")
        stderr.write_text(run.stderr, encoding="utf-8")
        calls.append({"command": command, "exit_code": run.returncode,
                      "stdout": stdout.name, "stderr": stderr.name})
        if run.returncode:
            raise RuntimeError(f"GENtle call {ordinal} failed; inspect {stderr}")

    entry = json.loads((FIXTURE / "ensembl_entry.json").read_text())
    # The imported entry is already gene-oriented by GENtle; hash it unchanged.
    sequence_sha = hashlib.sha256(entry["sequence"].upper().encode("ascii")).hexdigest()
    sources = []
    for provider, filename, chromosome, gene_id, release, accession in [
        ("ensembl", "ensembl.gff3", "22", "gene:ENSG00000100105", "Ensembl 116", "ENSG00000100105.20"),
        ("ref_seq", "refseq.gff3", "NC_000022.11", "gene-PATZ1", manifest["refseq_snapshot"], "NC_000022.11 / GeneID:23598"),
    ]:
        sources.append({"provider": provider, "format": "gff3",
                        "path": str(FIXTURE / filename), "sha256": manifest["files"][filename],
                        "assembly": "GRCh38", "release": release, "accession": accession,
                        "chromosome": chromosome, "gene_ids": [gene_id],
                        "locus_sequence_sha256": sequence_sha})
    write_json(output / "annotation-sources.json", sources)
    request = {
        "schema": "gentle.gene_locus_evidence_preparation_request.v1",
        "gene_query": "ENSG00000100105", "species": "homo_sapiens", "assembly": "GRCh38",
        "allow_ensembl_network": False, "ensembl_entry_path": str(FIXTURE / "ensembl_entry.json"),
        "flank_5prime_bp": 0, "flank_3prime_bp": 0,
        "display_upstream_bp": 0, "display_downstream_bp": 0,
        "output_seq_id": "patz1_ensembl_116", "output_panel_id": "patz1_ensembl_116",
        "annotation_release": "Ensembl 116", "transcript_annotation_sources": sources,
        "svg_path": str(output / "patz1-source-comparison.svg"),
        "display_report_path": str(output / "locus.report.json"),
        "receipt_path": str(output / "locus.receipt.json"),
    }
    write_json(output / "locus.request.json", request)
    cli("shell", f'gene-locus prepare @"{output / "locus.request.json"}"')
    report = json.loads((output / "locus.report.json").read_text())
    evidence = report["isoform_evidence"]
    write_json(output / "isoform-evidence.json", evidence)
    feature_id = evidence["splicing"]["target_feature_id"]
    design = {"DesignTranscriptAssayPanel": {
        "seq_id": report["seq_id"], "source_feature_id": feature_id,
        "assay_kind": "sybr_qpcr", "objective": "minimal_discrimination_panel",
        "coverage_policy": "best_effort", "assay_tier": "isoform_discrimination",
        "cdna_synthesis": "unspecified",
        "min_amplicon_bp": 80, "max_amplicon_bp": 250, "max_tm_delta_c": 3,
        "max_assays_per_class": 4, "max_mismatches": 0, "require_3prime_exact_bases": 8,
        "forward": {"min_length": 20, "max_length": 26, "min_tm_c": 58, "max_tm_c": 64,
                    "min_gc_fraction": 0.3, "max_gc_fraction": 0.7},
        "reverse": {"min_length": 20, "max_length": 26, "min_tm_c": 58, "max_tm_c": 64,
                    "min_gc_fraction": 0.3, "max_gc_fraction": 0.7},
        "annotation_release": "Ensembl 116", "report_id": "patz1_real_discrimination",
        "path": str(output / "primer-panel.json"),
    }}
    write_json(output / "primer-discrimination.operation.json", design)
    study = {
        "schema": "gentle.gene_isoform_assay_study_plan_request.v1",
        "plan_id": "patz1_real_study", "label": "Human PATZ1: annotation-led study, no expression evidence supplied",
        "isoform_evidence_path": str(output / "isoform-evidence.json"),
        "expected_isoform_evidence_sha256": "sha256:" + sha(output / "isoform-evidence.json"),
    }
    write_json(output / "study.request.json", study)
    cli("primers", "plan-gene-isoform-study", "@" + str(output / "study.request.json"),
        "--normalized-request", output / "study.normalized.json", "--path", output / "study.plan.json",
        "--workflow", output / "study.workflow.json")
    publication = {
        "schema": "gentle.gene_isoform_assay_publication_request.v1",
        "report_id": "patz1_real_tutorial", "title": "Human PATZ1 annotation-led primer study",
        "genes": [{"gene_symbol": "PATZ1", "status": "pending",
                   "study_plan": {"path": "study.plan.json", "expected_sha256": "sha256:" + sha(output / "study.plan.json")},
                   "status_reason": "Study not executed. The separate discrimination request is not a study result. Reference-wide specificity and laboratory validation are pending."}],
    }
    write_json(output / "publication.request.json", publication)
    receipt = {
        "schema": "gentle.tutorial_real_patz1_preparation.v1", "synthetic": False,
        "study_executed": False, "primer_design_executed": False, "gui_accepted": False,
        "binary": str(binary), "binary_sha256": sha(binary), "fixture_manifest_sha256": sha(FIXTURE / "manifest.json"),
        "calls": calls, "artifacts": [{"path": str(p.relative_to(output)), "sha256": sha(p)}
                                       for p in sorted(output.rglob("*")) if p.is_file()],
    }
    write_json(output / "preparation-receipt.json", receipt)
    return output


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gentle-cli", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    print(prepare(args.gentle_cli, args.output_dir))
