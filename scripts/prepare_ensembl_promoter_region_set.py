#!/usr/bin/env python3
"""Capture promoter features from bound GENtle locus reports as a region set.

This companion does not reinterpret Ensembl intervals. It passes the exact
normalized feature row and verified source binding back through GENtle's
``regions capture`` operation, then exports the resulting canonical
``gentle.genomic_region_set.v1`` document for sequence comparison.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import re
import subprocess
import sys
from typing import Any

try:
    from .prepare_regulatory_region_indexes import (
        checked_id, run_command, sha256_file, write_json, require,
    )
except ImportError:  # Direct ``python3 scripts/...`` execution.
    from prepare_regulatory_region_indexes import (  # type: ignore[no-redef]
        checked_id, run_command, sha256_file, write_json, require,
    )


REPORT_SCHEMA = "gentle.gene_locus_evidence_display.v1"
RECEIPT_SCHEMA = "gentle.ensembl_promoter_region_set_preparation_receipt.v1"
NON_ID = re.compile(r"[^A-Za-z0-9_.-]+")


def safe_component(value: str) -> str:
    cleaned = NON_ID.sub("_", value).strip("_.-")
    return checked_id(cleaned, "derived region id component")


def capture_requests(report_path: Path, set_id: str,
                     feature_types: set[str]) -> list[dict[str, Any]]:
    report = json.loads(report_path.read_text(encoding="utf-8"))
    require(report.get("schema") == REPORT_SCHEMA,
            f"{report_path}: expected schema '{REPORT_SCHEMA}'")
    ensembl = report.get("ensembl_regulation")
    require(isinstance(ensembl, dict), f"{report_path}: Ensembl Regulation evidence is absent")
    require(ensembl.get("availability") == "available",
            f"{report_path}: Ensembl Regulation evidence is not available")
    binding = ensembl.get("source_binding")
    require(isinstance(binding, dict), f"{report_path}: source_binding is absent")
    require(binding.get("content_identity_verified") is True,
            f"{report_path}: Ensembl source content identity is not verified")
    require(binding.get("truncated") is False,
            f"{report_path}: Ensembl overlap rows are truncated")
    source = binding.get("source")
    require(isinstance(source, dict), f"{report_path}: source descriptor is absent")
    chromosome = report.get("isoform_evidence", {}).get("chromosome")
    require(isinstance(chromosome, str) and chromosome.strip(),
            f"{report_path}: isoform_evidence.chromosome is absent")
    gene_symbol = report.get("gene_symbol")
    require(isinstance(gene_symbol, str) and gene_symbol.strip(),
            f"{report_path}: gene_symbol is absent")
    rows = ensembl.get("rows")
    require(isinstance(rows, list), f"{report_path}: Ensembl rows must be an array")

    reference = {
        "species_scientific_name": source.get("species_scientific_name"),
        "taxon_id": source.get("taxon_id"),
        "assembly_name": source.get("assembly_name"),
        "assembly_accession": source.get("assembly_accession"),
        "contig_name": chromosome,
        "contig_aliases": [f"chr{chromosome}"] if not chromosome.startswith("chr") else [],
    }
    report_digest = sha256_file(report_path)
    requests: list[dict[str, Any]] = []
    for row in rows:
        require(isinstance(row, dict), f"{report_path}: Ensembl row must be an object")
        feature_type = row.get("feature_type")
        if feature_type not in feature_types:
            continue
        feature_id = row.get("feature_id")
        require(isinstance(feature_id, str) and feature_id.strip(),
                f"{report_path}: selected row lacks feature_id")
        region_id = safe_component(f"{gene_symbol}_{feature_id}_{feature_type}")
        requests.append({
            "set_id": set_id,
            "set_label": "Ensembl Regulation promoter features",
            "region_id": region_id,
            "label": f"{gene_symbol} {feature_type} {feature_id}",
            "description": (
                "Release-bound Ensembl Regulation feature selected for sequence comparison; "
                "its annotation does not establish activity or regulation of the named gene."
            ),
            "purpose": "promoter_region",
            "display_color_hex": "#7C3AED",
            "source": {
                "source_kind": "gene_locus_ensembl_regulatory_feature",
                "row": row,
                "source_binding": binding,
                "reference": reference,
                "interval_kind": "core",
                "evidence_statement": ensembl.get("evidence_statement", ""),
                "non_claims": ensembl.get("non_claims", []),
            },
            "notes": [
                f"Source gene-locus report: {report_path.name}",
                f"Source gene-locus report digest: {report_digest}",
            ],
            "collision_policy": "reject",
        })
    return requests


def prepare(args: argparse.Namespace) -> None:
    repo_root = args.repo_root.resolve(strict=True)
    gentle = args.gentle.resolve(strict=True)
    output = args.output.resolve()
    require(not output.exists() or not any(output.iterdir()),
            "Output directory must be absent or empty")
    output.mkdir(parents=True, exist_ok=True)
    logs = output / "logs"
    requests_dir = output / "requests"
    logs.mkdir()
    requests_dir.mkdir()
    state_path = output / "ensembl_promoters.project.gentle.json"
    region_set_path = output / "ensembl_promoters.region_set.json"
    set_id = checked_id(args.set_id, "set_id")
    feature_types = set(args.feature_type or ["promoter"])
    require(feature_types, "at least one --feature-type is required")
    receipts: list[dict[str, Any]] = []
    version = run_command([str(gentle), "--version"], cwd=repo_root, timeout=30,
                          log_dir=logs, label="gentle_version", receipts=receipts)

    all_requests: list[dict[str, Any]] = []
    reports: list[dict[str, str]] = []
    for report_value in args.report:
        report_path = report_value.resolve(strict=True)
        requests = capture_requests(report_path, set_id, feature_types)
        require(requests, f"{report_path}: no requested feature types were found")
        all_requests.extend(requests)
        reports.append({"path": str(report_path), "sha256": sha256_file(report_path)})
    ids = [request["region_id"] for request in all_requests]
    require(len(ids) == len(set(ids)), "selected Ensembl rows produce duplicate region IDs")

    for index, request in enumerate(all_requests, start=1):
        request_path = requests_dir / f"{index:03d}_{request['region_id']}.capture.json"
        write_json(request_path, request)
        run_command([
            str(gentle), "--state", str(state_path), "shell",
            f"regions capture @{request_path}",
        ], cwd=repo_root, timeout=args.timeout, log_dir=logs,
            label=f"{index:03d}_{request['region_id']}_capture", receipts=receipts)

    export_request_path = requests_dir / "export.json"
    write_json(export_request_path, {"set_id": set_id, "json_path": str(region_set_path)})
    run_command([
        str(gentle), "--state", str(state_path), "shell",
        f"regions export @{export_request_path}",
    ], cwd=repo_root, timeout=args.timeout, log_dir=logs,
        label="export_region_set", receipts=receipts)
    exported = json.loads(region_set_path.read_text(encoding="utf-8"))
    require(exported.get("schema") == "gentle.genomic_region_set.v1",
            "GENtle export did not produce a canonical genomic region set")
    require(len(exported.get("regions", [])) == len(all_requests),
            "GENtle export region count does not match captured requests")

    write_json(output / "receipt.json", {
        "schema": RECEIPT_SCHEMA,
        "set_id": set_id,
        "feature_types": sorted(feature_types),
        "source_reports": reports,
        "gentle_binary": str(gentle),
        "gentle_binary_sha256": sha256_file(gentle),
        "gentle_version_stdout": version.stdout.decode(errors="replace").strip(),
        "commands": receipts,
        "output_region_set": str(region_set_path),
        "output_region_set_sha256": sha256_file(region_set_path),
        "region_count": len(all_requests),
        "non_claims": [
            "Ensembl feature type promoter is a provider annotation, not evidence of activity in a biosample.",
            "Geometric association with a gene does not establish causal regulation.",
            "Sequence similarity does not establish promoter function or functional equivalence.",
        ],
    })
    print(json.dumps({
        "status": "ok", "set_id": set_id, "regions": len(all_requests),
        "region_set": str(region_set_path), "receipt": str(output / "receipt.json"),
    }, indent=2))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", action="append", type=Path, required=True,
                        help="Repeat for each bound gentle.gene_locus_evidence_display.v1 report")
    parser.add_argument("--set-id", required=True)
    parser.add_argument("--feature-type", action="append",
                        help="Ensembl feature type to capture (default: promoter)")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--gentle", type=Path,
                        default=Path(__file__).resolve().parents[1] / "target/debug/gentle_cli")
    parser.add_argument("--timeout", type=int, default=3600)
    args = parser.parse_args()
    require(args.timeout > 0, "--timeout must be positive")
    prepare(args)


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, RuntimeError, subprocess.TimeoutExpired) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(2)
