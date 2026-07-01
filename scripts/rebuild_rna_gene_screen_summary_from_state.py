#!/usr/bin/env python3
"""Rebuild a canonical RNA gene-screen summary from an older GENtle state.

This is a recovery/normalization helper for legacy `rna-reads batch-map` runs
whose `batch_summary.tsv` kept seed-passed counts but did not yet export the
strict seed-passed read-length quantiles. The persisted GENtle state may still
contain the full RNA-read reports with `read_length_counts_seed_passed`; this
script joins those reports back to the batch summary by report id and writes the
canonical `gentle.rna_read_gene_screen_summary.v1` TSV used by the family plotter.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import sys
from pathlib import Path
from typing import Any

SCHEMA = "gentle.rna_read_gene_screen_summary.v1"

OUTPUT_FIELDS = [
    "schema",
    "gene",
    "run_accession",
    "sample_id",
    "sample_name",
    "source_kind",
    "source_path",
    "analysis_phase",
    "report_id",
    "seq_id",
    "seed_feature_id",
    "input_path",
    "total_reads",
    "all_q0_bp",
    "all_q25_bp",
    "all_q50_bp",
    "all_q75_bp",
    "all_q90_bp",
    "all_q95_bp",
    "all_q99_bp",
    "all_q100_bp",
    "all_mean_bp",
    "seed_passed_reads",
    "seed_passed_per_million",
    "seed_passed_q90_bp",
    "seed_passed_q95_bp",
    "seed_passed_q99_bp",
    "seed_passed_max_bp",
    "seed_passed_mean_bp",
    "accepted_target_count",
    "accepted_target_per_million",
    "accepted_target_max_bp",
    "accepted_target_mean_bp",
    "align_selection",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Rebuild canonical RNA gene-screen summary TSV from legacy "
            "batch_summary.tsv plus persisted GENtle RNA-read reports."
        )
    )
    parser.add_argument("--gene", required=True, help="Gene label for output rows, e.g. TP53.")
    parser.add_argument("--state", required=True, type=Path, help="GENtle state JSON containing rna_read_reports metadata.")
    parser.add_argument("--batch-summary", required=True, type=Path, help="Legacy out/batch_summary.tsv.")
    parser.add_argument("--output", required=True, type=Path, help="Output canonical gene_screen_summary.tsv.")
    parser.add_argument(
        "--allow-missing-reports",
        action="store_true",
        help="Keep rows without a matching persisted report, with blank strict seed-length columns.",
    )
    return parser.parse_args()


def text(value: Any) -> str:
    if value is None:
        return ""
    return str(value)


def first(row: dict[str, Any], *keys: str) -> str:
    for key in keys:
        value = row.get(key)
        if value is not None and str(value).strip() != "":
            return str(value)
    return ""


def number(value: Any, default: float = 0.0) -> float:
    if value is None:
        return default
    raw = str(value).strip()
    if not raw:
        return default
    try:
        parsed = float(raw)
    except ValueError:
        return default
    if math.isnan(parsed) or math.isinf(parsed):
        return default
    return parsed


def number_text(value: Any) -> str:
    raw = str(value).strip() if value is not None else ""
    if not raw:
        return ""
    parsed = number(raw, math.nan)
    if math.isnan(parsed):
        return raw
    if parsed.is_integer():
        return str(int(parsed))
    return f"{parsed:.12g}"


def total_count(counts: Any) -> int:
    if not isinstance(counts, list):
        return 0
    return sum(int(c or 0) for c in counts)


def total_bases(counts: Any) -> int:
    if not isinstance(counts, list):
        return 0
    return sum(idx * int(c or 0) for idx, c in enumerate(counts))


def quantile_from_counts(counts: Any, q: float) -> str:
    if not isinstance(counts, list):
        return ""
    total = total_count(counts)
    if total <= 0:
        return ""
    rank = min(max(math.ceil(total * max(0.0, min(1.0, q))), 1), total)
    cumulative = 0
    for length_bp, count in enumerate(counts):
        cumulative += int(count or 0)
        if cumulative >= rank:
            return str(length_bp)
    return str(max(len(counts) - 1, 0))


def mean_from_counts(counts: Any) -> str:
    total = total_count(counts)
    if total <= 0:
        return ""
    return number_text(total_bases(counts) / total)


def load_reports(state_path: Path) -> dict[str, dict[str, Any]]:
    state = json.loads(state_path.read_text(encoding="utf-8"))
    metadata = state.get("metadata") if isinstance(state, dict) else None
    store = metadata.get("rna_read_reports") if isinstance(metadata, dict) else None
    reports = store.get("reports") if isinstance(store, dict) else None
    if isinstance(reports, dict):
        return {str(k): v for k, v in reports.items() if isinstance(v, dict)}
    raise SystemExit(f"No metadata.rna_read_reports.reports object found in {state_path}")


def extract_srr(*values: Any) -> str:
    for value in values:
        match = re.search(r"SRR\d+", str(value or ""))
        if match:
            return match.group(0)
    return ""


def main() -> int:
    args = parse_args()
    reports = load_reports(args.state)
    rows = list(csv.DictReader(args.batch_summary.open(newline="", encoding="utf-8"), delimiter="\t"))
    output_rows: list[dict[str, Any]] = []
    missing: list[str] = []

    for row in rows:
        report_id = first(row, "report_id")
        report = reports.get(report_id) if report_id else None
        if report is None:
            missing.append(report_id or first(row, "sample_id", "run_accession") or "<unknown>")
            if not args.allow_missing_reports:
                continue
            report = {}

        total_reads = number(first(row, "total_reads", "read_count_total") or report.get("read_count_total"))
        seed_reads = number(first(row, "seed_passed_reads", "read_count_seed_passed") or report.get("read_count_seed_passed"))
        accepted = number(first(row, "accepted_target_count") or row.get("gene_support_accepted_target_count"))
        run = first(row, "run_accession", "sra_accession") or extract_srr(report_id, row.get("input_path"), report.get("input_path"))
        all_counts = report.get("read_length_counts_all")
        seed_counts = report.get("read_length_counts_seed_passed")

        output_rows.append(
            {
                "schema": SCHEMA,
                "gene": args.gene.upper(),
                "run_accession": run,
                "sample_id": first(row, "sample_id"),
                "sample_name": first(row, "sample_name"),
                "source_kind": "recovered_legacy_state",
                "source_path": str(args.state),
                "analysis_phase": first(row, "analysis_phase") or "legacy_batch_recovered",
                "report_id": report_id,
                "seq_id": first(row, "seq_id") or text(report.get("seq_id")),
                "seed_feature_id": first(row, "seed_feature_id") or text(report.get("seed_feature_id")),
                "input_path": first(row, "input_path") or text(report.get("input_path")),
                "total_reads": number_text(total_reads),
                "all_q0_bp": quantile_from_counts(all_counts, 0.0) or first(row, "all_q0_bp"),
                "all_q25_bp": quantile_from_counts(all_counts, 0.25) or first(row, "all_q25_bp"),
                "all_q50_bp": quantile_from_counts(all_counts, 0.50) or first(row, "all_q50_bp"),
                "all_q75_bp": quantile_from_counts(all_counts, 0.75) or first(row, "all_q75_bp"),
                "all_q90_bp": quantile_from_counts(all_counts, 0.90) or first(row, "all_q90_bp"),
                "all_q95_bp": quantile_from_counts(all_counts, 0.95) or first(row, "all_q95_bp"),
                "all_q99_bp": quantile_from_counts(all_counts, 0.99) or first(row, "all_q99_bp"),
                "all_q100_bp": quantile_from_counts(all_counts, 1.0) or first(row, "all_q100_bp"),
                "all_mean_bp": mean_from_counts(all_counts) or first(row, "all_mean_bp", "mean_read_length_bp"),
                "seed_passed_reads": number_text(seed_reads),
                "seed_passed_per_million": number_text(seed_reads * 1_000_000.0 / total_reads) if total_reads > 0 else "0",
                "seed_passed_q90_bp": quantile_from_counts(seed_counts, 0.90),
                "seed_passed_q95_bp": quantile_from_counts(seed_counts, 0.95),
                "seed_passed_q99_bp": quantile_from_counts(seed_counts, 0.99),
                "seed_passed_max_bp": quantile_from_counts(seed_counts, 1.0),
                "seed_passed_mean_bp": mean_from_counts(seed_counts),
                "accepted_target_count": number_text(accepted),
                "accepted_target_per_million": number_text(accepted * 1_000_000.0 / total_reads) if total_reads > 0 else "0",
                "accepted_target_max_bp": first(row, "accepted_target_max_bp"),
                "accepted_target_mean_bp": first(row, "accepted_target_mean_bp"),
                "align_selection": first(row, "align_selection") or "legacy_batch",
            }
        )

    if missing:
        print(
            f"WARNING: {len(missing)} batch row(s) had no matching persisted RNA-read report: {', '.join(missing[:8])}",
            file=sys.stderr,
        )
        if not args.allow_missing_reports:
            print("Use --allow-missing-reports to keep those rows with blank strict seed lengths.", file=sys.stderr)
    if not output_rows:
        return 1

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=OUTPUT_FIELDS)
        writer.writeheader()
        writer.writerows(output_rows)
    print(f"Wrote {args.output} ({len(output_rows)} row(s))")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
