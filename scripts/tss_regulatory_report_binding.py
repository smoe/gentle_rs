"""Exact input bindings shared by TSS-feature preparation and both renderers."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any
import xml.etree.ElementTree as ET

try:
    from .compare_candidates_to_promoterome import require, same_digest
except ImportError:
    from compare_candidates_to_promoterome import require, same_digest

LOCUS_SCHEMA = "gentle.gene_locus_evidence_display.v1"


def read_locus_report(path: Path) -> tuple[dict[str, Any], str]:
    payload = path.read_bytes()
    report = json.loads(payload)
    require(report.get("schema") == LOCUS_SCHEMA, "expected a portable gene-locus display report")
    require(all(report.get(key) for key in ("gene_symbol", "panel_id", "seq_id")),
            "locus report lacks gene/panel/sequence identity")
    return report, hashlib.sha256(payload).hexdigest()


def validate_locus_svg(svg: str, report: dict[str, Any]) -> None:
    """Check declared panel identity, not biological truth of the drawing."""
    root = ET.fromstring(svg)
    require(root.tag in {"svg", "{http://www.w3.org/2000/svg}svg"}
            and root.get("data-gentle-schema") == LOCUS_SCHEMA
            and root.get("data-gentle-panel-id") == report["panel_id"],
            "base SVG does not match the locus report's schema/panel identity")
    require(not any(item.get("data-gentle-panel") == "tss-local-promoter-similarity"
                    or item.get("data-gentle-shifted-footer")
                    or item.get("data-gentle-tss-stretch-overview") for item in root.iter()),
            "base SVG already contains a TSS similarity section")


def load_bound_locus_report(
    candidates: dict[str, Any], path: Path, gene: str | None = None,
) -> tuple[dict[str, Any], str]:
    report, digest = read_locus_report(path)
    gene = gene or report["gene_symbol"]
    require(report["gene_symbol"] == gene, "locus report belongs to another gene")
    expected = candidates.get("source_bindings", {}).get("locus_reports", {}).get(gene)
    require(isinstance(expected, str) and same_digest(digest, expected),
            f"locus report hash mismatch for {gene}; use the preparation-bound report")
    regions = [row for row in candidates["regions"] if row["gene_query"] == gene]
    require(regions, f"no candidate regions for locus report {gene}")
    require(all(same_digest(row["source_region"]["source_report_sha256"], digest)
                for row in regions), "per-region locus report binding mismatch")
    return report, digest


def load_bound_locus_svg(
    candidates: dict[str, Any], path: Path, report: dict[str, Any],
) -> tuple[str, str]:
    payload = path.read_bytes()
    digest = hashlib.sha256(payload).hexdigest()
    gene = report["gene_symbol"]
    expected = candidates.get("source_bindings", {}).get("locus_svgs", {}).get(gene)
    require(isinstance(expected, str) and same_digest(digest, expected),
            f"base SVG is not preparation-bound for {gene}; reprepare with --locus-svg GENE=PATH")
    svg = payload.decode("utf-8")
    validate_locus_svg(svg, report)
    return svg, digest
