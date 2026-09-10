#!/usr/bin/env python3
"""Compose one bound locus-context page with its selected TSS TFBS pages."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
from typing import Any
import xml.etree.ElementTree as ET


LOCUS_RECEIPT_SCHEMA = "gentle.tss_local_similarity_locus_report_receipt.v1"
LOCUS_SVG_SCHEMA = "gentle.gene_locus_evidence_display.v1"
TSS_REPORT_SCHEMA = "gentle.tss_tfbs_profiles.v1"
TSS_INDEX_SCHEMA = "gentle.tss_tfbs_profile_index.v1"
TSS_RECEIPT_SCHEMA = "gentle.tss_tfbs_profile_receipt.v1"
OUTPUT_SCHEMA = "gentle.integrated_locus_tss_tfbs_pdf_receipt.v1"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_bytes())
    require(isinstance(value, dict), f"{path} is not a JSON object")
    return value


def normalized_digest(value: Any) -> str:
    require(isinstance(value, str), "expected a SHA-256 string")
    digest = value.removeprefix("sha256:")
    require(len(digest) == 64 and all(c in "0123456789abcdef" for c in digest),
            "expected a lowercase SHA-256 digest")
    return digest


def bound_output(receipt: dict[str, Any], relative_name: str, path: Path) -> None:
    outputs = receipt.get("outputs")
    require(isinstance(outputs, dict), "receipt outputs must be an object")
    require(relative_name in outputs, f"receipt does not bind output {relative_name}")
    require(normalized_digest(outputs[relative_name]) == sha256(path),
            f"output hash mismatch for {relative_name}")


def locus_bands(svg_path: Path) -> tuple[dict[str, str], list[tuple[str, int, int]]]:
    root = ET.fromstring(svg_path.read_bytes())
    require(root.get("data-gentle-schema") == LOCUS_SVG_SCHEMA,
            "locus SVG has the wrong scientific schema")
    bands = []
    for item in root.iter():
        name = item.get("data-gentle-tss-stretch-band")
        if name is None:
            continue
        start = int(item.get("data-gentle-genomic-start", "0"))
        end = int(item.get("data-gentle-genomic-end", "0"))
        require(name and 0 < start <= end, "locus SVG contains an invalid TSS band")
        bands.append((name, start, end))
    require(bands, "locus SVG contains no TSS background bands")
    return root.attrib, bands


def select_pages(
    gene: str,
    report: dict[str, Any],
    index: dict[str, Any],
    tss_dir: Path,
    tss_receipt: dict[str, Any],
    bands: list[tuple[str, int, int]],
) -> tuple[list[Path], list[dict[str, Any]]]:
    selected = [window for window in report.get("windows", [])
                if window.get("selected") is True
                and window.get("record", {}).get("gene_symbol") == gene]
    require(selected, f"TSS report has no selected windows for {gene}")
    order = {row["record"]["promoter_id"]: i for i, row in enumerate(selected)}
    require(len(order) == len(selected), "selected promoter IDs are not unique")
    bindings = []
    for window in selected:
        record = window["record"]
        geometry = record.get("geometry", {})
        tss = geometry.get("tss_1based")
        require(isinstance(tss, int) and tss > 0, "selected TSS coordinate is invalid")
        covering = [name for name, start, end in bands if start <= tss <= end]
        require(covering, f"selected TSS {record['promoter_id']} is outside every locus band")
        require(geometry.get("strand") in {"+", "-"}, "selected TSS strand is invalid")
        require(isinstance(geometry.get("chromosome"), str), "selected chromosome is missing")
        evidence = window.get("selection_evidence")
        require(isinstance(evidence, dict) and evidence.get("label"),
                "selected TSS lacks bound selection evidence")
        bindings.append({
            "promoter_id": record["promoter_id"],
            "gene_id": record["gene_id"],
            "chromosome": geometry["chromosome"],
            "strand": geometry["strand"],
            "tss_1based": tss,
            "start_1based": geometry["start_1based"],
            "end_1based": geometry["end_1based"],
            "sequence_sha256": record["sequence_sha256"],
            "transcripts": record["transcripts"],
            "covering_tss_bands": covering,
            "selection_label": evidence["label"],
            "selection_factor": evidence.get("factor"),
            "selection_criterion": evidence.get("criterion"),
        })

    genes = [row for row in index.get("genes", []) if row.get("gene_symbol") == gene]
    require(len(genes) == 1, f"TSS index must contain exactly one {gene} entry")
    chosen: list[tuple[int, Path]] = []
    covered: set[str] = set()
    for page in genes[0].get("pages", []):
        ids = page.get("promoter_ids", [])
        overlap = [promoter_id for promoter_id in ids if promoter_id in order]
        if not overlap:
            continue
        require(len(overlap) == len(ids),
                "a selected detail page also contains unselected TSS panels")
        svgs = [name for name in page.get("files", []) if name.endswith(".svg")]
        require(len(svgs) == 1, "selected detail page must bind exactly one SVG")
        path = tss_dir / svgs[0]
        require(path.is_file(), f"selected detail SVG is missing: {path}")
        bound_output(tss_receipt, svgs[0], path)
        first = min(order[promoter_id] for promoter_id in overlap)
        chosen.append((first, path))
        for promoter_id in overlap:
            require(promoter_id not in covered, "selected TSS appears on multiple detail pages")
            covered.add(promoter_id)
    require(covered == set(order), "selected detail pages do not cover the exact selected TSS set")
    chosen.sort(key=lambda row: row[0])
    bindings.sort(key=lambda row: order[row["promoter_id"]])
    return [path for _, path in chosen], bindings


def compose(args: argparse.Namespace) -> dict[str, Any]:
    locus_svg = args.locus_svg.resolve()
    locus_receipt_path = args.locus_receipt.resolve()
    tss_report_path = args.tss_report.resolve()
    tss_index_path = args.tss_index.resolve()
    tss_receipt_path = args.tss_receipt.resolve()
    tss_dir = tss_report_path.parent
    output_pdf = args.output_pdf.resolve()
    output_receipt = args.output_receipt.resolve()
    gentle_cli = args.gentle_cli.resolve()
    require(output_pdf != output_receipt, "output PDF and receipt must be distinct")
    require(not output_pdf.exists() and not output_receipt.exists(),
            "outputs must not already exist")
    require(gentle_cli.is_file(), "gentle_cli does not exist")

    locus_receipt = read_json(locus_receipt_path)
    require(locus_receipt.get("schema") == LOCUS_RECEIPT_SCHEMA,
            "unexpected locus receipt schema")
    require(locus_receipt.get("gene") == args.gene, "locus receipt belongs to another gene")
    bound_output(locus_receipt, locus_svg.name, locus_svg)
    root, bands = locus_bands(locus_svg)

    tss_report = read_json(tss_report_path)
    tss_index = read_json(tss_index_path)
    tss_receipt = read_json(tss_receipt_path)
    require(tss_report.get("schema") == TSS_REPORT_SCHEMA, "unexpected TSS report schema")
    require(tss_index.get("schema") == TSS_INDEX_SCHEMA, "unexpected TSS index schema")
    require(tss_receipt.get("schema") == TSS_RECEIPT_SCHEMA, "unexpected TSS receipt schema")
    bound_output(tss_receipt, tss_report_path.name, tss_report_path)
    bound_output(tss_receipt, tss_index_path.name, tss_index_path)
    require(normalized_digest(tss_index.get("report_sha256")) == sha256(tss_report_path),
            "TSS index/report binding mismatch")

    detail_pages, bindings = select_pages(
        args.gene, tss_report, tss_index, tss_dir, tss_receipt, bands,
    )
    page_paths = [locus_svg, *detail_pages]
    output_pdf.parent.mkdir(parents=True, exist_ok=True)
    partial = output_pdf.with_name(output_pdf.name + ".partial")
    require(not partial.exists(), "stale partial output exists")
    command = [str(gentle_cli), "svg-pdf-set", str(partial),
               *[str(path) for path in page_paths]]
    try:
        result = subprocess.run(command, check=True, capture_output=True, text=True, timeout=900)
        summary = json.loads(result.stdout)
        require(summary.get("page_count") == len(page_paths),
                "multi-page renderer returned the wrong page count")
        require(partial.is_file() and partial.stat().st_size > 0,
                "multi-page renderer produced no PDF")
        os.replace(partial, output_pdf)
    except Exception:
        if partial.exists():
            partial.unlink()
        raise

    receipt = {
        "schema": OUTPUT_SCHEMA,
        "gene_symbol": args.gene,
        "gene_id": bindings[0]["gene_id"],
        "join_policy": "selected promoter_id plus gene/chromosome/TSS/strand; TSS must fall inside a bound locus background band",
        "page_order": ["locus_context", *["selected_tss_tfbs" for _ in detail_pages]],
        "page_count": len(page_paths),
        "selected_tss_count": len(bindings),
        "selected_tss_bindings": bindings,
        "inputs": {
            "locus_svg": {"path": str(locus_svg), "sha256": sha256(locus_svg)},
            "locus_receipt": {"path": str(locus_receipt_path), "sha256": sha256(locus_receipt_path)},
            "tss_report": {"path": str(tss_report_path), "sha256": sha256(tss_report_path)},
            "tss_index": {"path": str(tss_index_path), "sha256": sha256(tss_index_path)},
            "tss_receipt": {"path": str(tss_receipt_path), "sha256": sha256(tss_receipt_path)},
            "pages": [{"path": str(path), "sha256": sha256(path)} for path in page_paths],
        },
        "output": {
            "pdf": str(output_pdf),
            "sha256": sha256(output_pdf),
            "bytes": output_pdf.stat().st_size,
        },
        "producer": {
            "script": str(Path(__file__).resolve()),
            "script_sha256": sha256(Path(__file__).resolve()),
            "revision": args.producer_revision,
            "gentle_cli": str(gentle_cli),
            "gentle_cli_sha256": sha256(gentle_cli),
            "renderer_summary": summary,
        },
        "locus_panel_id": root.get("data-gentle-panel-id"),
        "score_policy": tss_report.get("score_policy"),
        "reference": tss_report.get("reference"),
        "verification": tss_report.get("verification"),
        "non_claims": tss_report.get("non_claims"),
    }
    payload = json.dumps(receipt, indent=2, sort_keys=True, ensure_ascii=False) + "\n"
    partial_receipt = output_receipt.with_name(output_receipt.name + ".partial")
    require(not partial_receipt.exists(), "stale partial receipt exists")
    partial_receipt.write_text(payload)
    os.replace(partial_receipt, output_receipt)
    return receipt


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gene", required=True)
    parser.add_argument("--locus-svg", type=Path, required=True)
    parser.add_argument("--locus-receipt", type=Path, required=True)
    parser.add_argument("--tss-report", type=Path, required=True)
    parser.add_argument("--tss-index", type=Path, required=True)
    parser.add_argument("--tss-receipt", type=Path, required=True)
    parser.add_argument("--gentle-cli", type=Path, required=True)
    parser.add_argument("--producer-revision", required=True)
    parser.add_argument("--output-pdf", type=Path, required=True)
    parser.add_argument("--output-receipt", type=Path, required=True)
    return parser.parse_args()


if __name__ == "__main__":
    try:
        result = compose(parse_args())
        print(json.dumps({
            "status": "ok",
            "schema": result["schema"],
            "gene_symbol": result["gene_symbol"],
            "page_count": result["page_count"],
            "selected_tss_count": result["selected_tss_count"],
            "output": result["output"],
        }, sort_keys=True))
    except Exception as error:
        print(f"error: {error}", file=sys.stderr)
        raise SystemExit(2)
