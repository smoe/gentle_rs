#!/usr/bin/env python3
"""Render pancreas RNA-read support plots for an arbitrary gene family.

The plotter deliberately consumes already-produced GENtle report tables instead
of re-running any biology. It makes one canonical family TSV first, then renders
an SVG with library read-length context and grouped per-gene support bars.
The primary support metric is always the conservative seed-passed read count.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import sys
from pathlib import Path
from xml.sax.saxutils import escape


DEFAULT_PALETTE = (
    "#2563eb",
    "#ea580c",
    "#16a34a",
    "#7c3aed",
    "#dc2626",
    "#0891b2",
    "#ca8a04",
    "#be185d",
)
READ_LENGTH_SERIES = [
    ("all_q25_bp", "q25", "#2563eb", ""),
    ("all_q50_bp", "q50", "#111827", ""),
    ("all_q75_bp", "q75", "#f59e0b", ""),
    ("all_q90_bp", "q90", "#7c3aed", "8,4"),
    ("all_q95_bp", "q95", "#0ea5e9", "4,3"),
    ("all_q99_bp", "q99", "#be123c", "3,3"),
    ("all_q100_bp", "q100", "#dc2626", ""),
    ("all_mean_bp", "mean", "#059669", "6,4"),
]
DEFAULT_READ_LENGTH_SERIES = [
    series for series in READ_LENGTH_SERIES if series[0] not in {"all_q99_bp", "all_q100_bp"}
]
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
    "support_length_source",
    "support_high_read_bp",
    "support_max_read_bp",
    "support_mean_read_bp",
    "source",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Render a grouped gene-family pancreas SVG from already-produced "
            "GENtle batch reports and figure-source TSVs."
        )
    )
    parser.add_argument(
        "--canonical-summary",
        action="append",
        default=[],
        metavar="PATH",
        help=(
            "Repeatable canonical gentle.rna_read_gene_screen_summary.v1 TSV. "
            "Prefer this when available; it already carries the gene column."
        ),
    )
    parser.add_argument(
        "--batch-report",
        action="append",
        default=[],
        metavar="GENE=PATH",
        help="Repeatable rna-reads batch_report.json source for one gene.",
    )
    parser.add_argument(
        "--batch-summary",
        action="append",
        default=[],
        metavar="GENE=PATH",
        help="Repeatable rna-reads batch_summary.tsv source for one gene.",
    )
    parser.add_argument(
        "--figure-source",
        action="append",
        default=[],
        metavar="GENE=PATH",
        help="Repeatable pancreas_gene_rna_screen.sh figure-source TSV for one gene.",
    )
    parser.add_argument(
        "--genes",
        default=None,
        help="Comma-separated display order. Defaults to input order.",
    )
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Output SVG path.",
    )
    parser.add_argument(
        "--output-tsv",
        default=None,
        help="Canonical family TSV path. Defaults to OUTPUT.with_suffix('.tsv').",
    )
    parser.add_argument(
        "--metric",
        choices=["per_million", "raw"],
        default="per_million",
        help="Grouped-bar height metric. Raw read counts are printed inside bars in both modes.",
    )
    parser.add_argument(
        "--support-length-stat",
        choices=["max", "mean"],
        default="max",
        help=(
            "Gene-supporting read-length statistic drawn as colored symbols in "
            "the lower panel. Whole-library q90 stays in the upper panel."
        ),
    )
    parser.add_argument(
        "--support-length-source",
        choices=["strict_seed_passed", "any"],
        default="strict_seed_passed",
        help=(
            "Which support-read length provenance may be drawn. The default only "
            "plots lengths measured on the same strict seed-passed population as "
            "the bars; 'any' also allows older accepted-target fallback lengths."
        ),
    )
    parser.add_argument(
        "--support-length-scale",
        choices=["log", "linear"],
        default="log",
        help="Scale for the lower support-read-length axis. Default: log.",
    )
    parser.add_argument(
        "--support-length-display",
        choices=["points", "lines"],
        default="points",
        help=(
            "How to draw support-read lengths in the lower panel. Default: "
            "points, because missing samples make connecting lines misleading."
        ),
    )
    parser.add_argument(
        "--support-length-genes",
        default=None,
        help=(
            "Optional comma-separated gene list for support-read length symbols. "
            "Bars are still drawn for all genes."
        ),
    )
    parser.add_argument(
        "--show-all-read-max",
        action="store_true",
        help=(
            "Draw the whole-library q99 and q100/max read-length lines in the "
            "upper panel. They are hidden by default because extreme read "
            "lengths are outlier statistics."
        ),
    )
    parser.add_argument(
        "--label-column",
        choices=["sample_id", "sample_name", "run_accession"],
        default="sample_id",
        help="Sample label shown on the x axis.",
    )
    parser.add_argument(
        "--title",
        default="Gene-family seed-passed Nanopore cDNA support",
        help="SVG title.",
    )
    return parser.parse_args()


def parse_gene_path(value: str, option: str) -> tuple[str, Path]:
    if "=" not in value:
        raise SystemExit(f"{option} expects GENE=PATH, got: {value}")
    gene, path = value.split("=", 1)
    gene = gene.strip().upper()
    path = path.strip()
    if not gene or not path:
        raise SystemExit(f"{option} expects non-empty GENE=PATH, got: {value}")
    return gene, Path(path)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def to_float(value: object, default: float = 0.0) -> float:
    if value is None:
        return default
    text = str(value).replace(",", "").strip()
    if text in {"", "null", "NA", "NaN", "nan"}:
        return default
    try:
        parsed = float(text)
    except ValueError:
        return default
    if math.isnan(parsed) or math.isinf(parsed):
        return default
    return parsed


def optional_number_text(value: object) -> str:
    numeric = to_float(value, math.nan)
    if math.isnan(numeric):
        return ""
    if numeric.is_integer():
        return str(int(numeric))
    return f"{numeric:.6g}"


def optional_float(row: dict[str, object], *keys: str) -> float | None:
    for key in keys:
        if key not in row:
            continue
        raw = row.get(key)
        if raw is None or str(raw).strip() == "":
            continue
        parsed = to_float(raw, math.nan)
        if not math.isnan(parsed):
            return parsed
    return None


def extract_srr(*values: object) -> str:
    for value in values:
        match = re.search(r"SRR\d+", str(value or ""))
        if match:
            return match.group(0)
    return ""


def sort_key_for_run(run: str) -> tuple[int, str]:
    match = re.search(r"SRR(\d+)", run or "")
    if match:
        return (int(match.group(1)), run)
    return (10**18, run or "")


def first_text(row: dict[str, object], *keys: str) -> str:
    for key in keys:
        value = row.get(key)
        if value is not None and str(value).strip() != "":
            return str(value).strip()
    return ""


def normalize_support_row(
    *,
    gene: str,
    source: Path,
    row: dict[str, object],
    seed_keys: tuple[str, ...],
    accepted_keys: tuple[str, ...],
) -> dict[str, object]:
    run = first_text(row, "run_accession", "sra_accession")
    if not run:
        run = extract_srr(row.get("report_id"), row.get("input_path"), row.get("source"))
    total = to_float(first_text(row, "total_reads", "read_count_total"))
    seed = 0.0
    for key in seed_keys:
        if key in row and str(row.get(key) or "").strip() != "":
            seed = to_float(row.get(key))
            break
    accepted = 0.0
    for key in accepted_keys:
        if key in row and str(row.get(key) or "").strip() != "":
            accepted = to_float(row.get(key))
            break
    out: dict[str, object] = {
        "schema": first_text(row, "schema") or "gentle.rna_read_gene_screen_summary.v1",
        "gene": gene,
        "run_accession": run,
        "sample_id": first_text(row, "sample_id"),
        "sample_name": first_text(row, "sample_name"),
        "source_kind": first_text(row, "source_kind") or "legacy_input",
        "source_path": first_text(row, "source_path") or str(source),
        "analysis_phase": first_text(row, "analysis_phase"),
        "report_id": first_text(row, "report_id"),
        "seq_id": first_text(row, "seq_id"),
        "seed_feature_id": first_text(row, "seed_feature_id"),
        "input_path": first_text(row, "input_path"),
        "total_reads": int(total) if float(total).is_integer() else total,
        "seed_passed_reads": int(seed) if float(seed).is_integer() else seed,
        "seed_passed_per_million": (seed * 1_000_000.0 / total) if total > 0 else 0.0,
        "seed_passed_q90_bp": first_text(row, "seed_passed_q90_bp"),
        "seed_passed_q95_bp": first_text(row, "seed_passed_q95_bp"),
        "seed_passed_q99_bp": first_text(row, "seed_passed_q99_bp"),
        "seed_passed_max_bp": first_text(row, "seed_passed_max_bp"),
        "seed_passed_mean_bp": first_text(row, "seed_passed_mean_bp"),
        "accepted_target_count": int(accepted) if float(accepted).is_integer() else accepted,
        "accepted_target_per_million": (accepted * 1_000_000.0 / total) if total > 0 else 0.0,
        "accepted_target_max_bp": first_text(row, "accepted_target_max_bp"),
        "accepted_target_mean_bp": first_text(row, "accepted_target_mean_bp"),
        "align_selection": first_text(row, "align_selection"),
        "support_length_source": "",
        "support_high_read_bp": "",
        "support_max_read_bp": "",
        "support_mean_read_bp": "",
        "source": str(source),
    }
    seed_high = optional_float(
        row,
        "strict_seed_passed_q90_read_bp",
        "strict_seed_passed_q90_read_length_bp",
        "strict_seed_passed_q90_bp",
        "strict_seed_passed_q95_read_bp",
        "strict_seed_passed_q95_read_length_bp",
        "strict_seed_passed_q95_bp",
        "strict_seed_passed_q99_read_bp",
        "strict_seed_passed_q99_read_length_bp",
        "strict_seed_passed_q99_bp",
        "seed_passed_q90_read_bp",
        "seed_passed_q90_read_length_bp",
        "seed_passed_q90_bp",
        "seed_passed_q95_read_bp",
        "seed_passed_q95_read_length_bp",
        "seed_passed_q95_bp",
        "seed_passed_q99_read_bp",
        "seed_passed_q99_read_length_bp",
        "seed_passed_q99_bp",
        "strict_seed_q90_bp",
        "strict_seed_q90_read_bp",
        "strict_seed_q90_read_length_bp",
        "strict_seed_q95_bp",
        "strict_seed_q95_read_bp",
        "strict_seed_q95_read_length_bp",
        "strict_seed_q99_bp",
        "strict_seed_q99_read_bp",
        "strict_seed_q99_read_length_bp",
    )
    seed_max = optional_float(
        row,
        "strict_seed_passed_max_read_bp",
        "strict_seed_passed_max_read_length_bp",
        "strict_seed_passed_q100_read_bp",
        "strict_seed_passed_q100_read_length_bp",
        "strict_seed_passed_q100_bp",
        "strict_seed_passed_max_bp",
        "seed_passed_max_read_bp",
        "seed_passed_max_read_length_bp",
        "seed_passed_q100_read_bp",
        "seed_passed_q100_read_length_bp",
        "seed_passed_q100_bp",
        "seed_passed_max_bp",
        "strict_seed_q100_bp",
        "strict_seed_q100_read_bp",
        "strict_seed_q100_read_length_bp",
        "strict_seed_max_read_bp",
        "strict_seed_max_read_length_bp",
        "strict_seed_max_bp",
    )
    seed_mean = optional_float(
        row,
        "strict_seed_passed_mean_read_bp",
        "strict_seed_passed_mean_read_length_bp",
        "strict_seed_passed_mean_bp",
        "seed_passed_mean_read_bp",
        "seed_passed_mean_read_length_bp",
        "seed_passed_mean_bp",
        "strict_seed_mean_bp",
        "strict_seed_mean_read_bp",
        "strict_seed_mean_read_length_bp",
    )
    if seed > 0 and (seed_high is not None or seed_max is not None or seed_mean is not None):
        out["support_length_source"] = "strict_seed_passed"
        out["support_high_read_bp"] = optional_number_text(seed_high)
        out["support_max_read_bp"] = optional_number_text(seed_max)
        out["support_mean_read_bp"] = optional_number_text(seed_mean)
        out["seed_passed_q90_bp"] = out["seed_passed_q90_bp"] or optional_number_text(seed_high)
        out["seed_passed_max_bp"] = out["seed_passed_max_bp"] or optional_number_text(seed_max)
        out["seed_passed_mean_bp"] = out["seed_passed_mean_bp"] or optional_number_text(seed_mean)
    else:
        accepted_high = optional_float(
            row,
            "accepted_target_q90_read_bp",
            "accepted_target_q90_read_length_bp",
            "accepted_target_p95_read_bp",
            "accepted_target_p95_read_length_bp",
            "accepted_tp73_q90_read_bp",
            "accepted_tp73_q90_read_length_bp",
        )
        accepted_max = optional_float(
            row,
            "accepted_target_max_read_bp",
            "accepted_target_max_read_length_bp",
            "accepted_target_q100_read_bp",
            "accepted_target_q100_read_length_bp",
            "accepted_target_max_bp",
            "accepted_tp73_max_read_bp",
            "accepted_tp73_max_read_length_bp",
            "accepted_tp73_q100_read_bp",
            "accepted_tp73_q100_read_length_bp",
            "accepted_max_bp",
        )
        accepted_mean = optional_float(
            row,
            "accepted_target_mean_read_bp",
            "accepted_target_mean_read_length_bp",
            "accepted_target_mean_bp",
            "accepted_tp73_mean_read_bp",
            "accepted_tp73_mean_read_length_bp",
            "accepted_tp73_mean_bp",
            "accepted_mean_bp",
        )
        if accepted > 0 and (accepted_high is not None or accepted_max is not None or accepted_mean is not None):
            out["support_length_source"] = "accepted_target"
            out["support_high_read_bp"] = optional_number_text(accepted_high)
            out["support_max_read_bp"] = optional_number_text(accepted_max)
            out["support_mean_read_bp"] = optional_number_text(accepted_mean)
            out["accepted_target_max_bp"] = out["accepted_target_max_bp"] or optional_number_text(accepted_max)
            out["accepted_target_mean_bp"] = out["accepted_target_mean_bp"] or optional_number_text(accepted_mean)
    for key in (
        "all_q0_bp",
        "all_q25_bp",
        "all_q50_bp",
        "all_q75_bp",
        "all_q90_bp",
        "all_q95_bp",
        "all_q99_bp",
        "all_q100_bp",
        "all_mean_bp",
    ):
        out[key] = optional_number_text(row.get(key, ""))
    return out


def load_canonical_rows(path: Path) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for row in read_tsv(path):
        gene = first_text(row, "gene").strip().upper()
        if not gene:
            gene = path.stem.split("_", 1)[0].upper()
        normalized = normalize_support_row(
            gene=gene,
            source=path,
            row=row,
            seed_keys=("seed_passed_reads", "strict_seed_passed_reads", "read_count_seed_passed"),
            accepted_keys=("accepted_target_count", "accepted_target_reads"),
        )
        normalized["source_kind"] = normalized.get("source_kind") or "canonical_summary"
        normalized["source_path"] = normalized.get("source_path") or str(path)
        rows.append(normalized)
    return rows


def maybe_fill_lengths_from_gene_support(
    row_out: dict[str, object],
    source_row: dict[str, object],
    source_path: Path,
) -> None:
    """Fill target-support read lengths from a gene-support summary JSON.

    Batch-map rows do not currently carry seed-passed read length histograms.
    They do, however, point at a gene-support summary with accepted-target read
    length distributions. Those are explicitly marked as accepted_target so the
    plot does not pretend they are whole-library or seed-passed quantiles.
    """
    if row_out.get("support_length_source"):
        return
    if to_float(row_out.get("accepted_target_count")) <= 0:
        return
    path_text = first_text(source_row, "gene_support_summary_json_path")
    if not path_text:
        return
    candidate = Path(path_text)
    if not candidate.is_absolute():
        candidate = source_path.parent / candidate
    if not candidate.exists():
        return
    try:
        summary = json.loads(candidate.read_text())
    except Exception:
        return
    lengths = summary.get("accepted_target_read_lengths") or {}
    if not isinstance(lengths, dict):
        return
    high = optional_float(lengths, "q90_length_bp", "p95_length_bp", "q75_length_bp")
    max_len = optional_float(lengths, "max_length_bp")
    mean = optional_float(lengths, "mean_length_bp")
    if high is None and max_len is None and mean is None:
        return
    row_out["support_length_source"] = "accepted_target"
    row_out["support_high_read_bp"] = optional_number_text(high)
    row_out["support_max_read_bp"] = optional_number_text(max_len)
    row_out["support_mean_read_bp"] = optional_number_text(mean)


def load_batch_rows(gene: str, path: Path, is_json: bool) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    if is_json:
        with path.open() as handle:
            report = json.load(handle)
        source_rows = report.get("rows")
        if not isinstance(source_rows, list):
            raise SystemExit(f"{gene} batch report has no rows[] array: {path}")
    else:
        source_rows = read_tsv(path)
    for row in source_rows:
        if not isinstance(row, dict):
            continue
        row_gene = str(row.get("gene") or "").strip().upper()
        if row_gene and row_gene != gene:
            continue
        normalized = normalize_support_row(
            gene=gene,
            source=path,
            row=row,
            seed_keys=("read_count_seed_passed", "seed_passed_reads"),
            accepted_keys=("accepted_target_count",),
        )
        maybe_fill_lengths_from_gene_support(normalized, row, path)
        rows.append(normalized)
    return rows


def load_figure_rows(gene: str, path: Path) -> list[dict[str, object]]:
    return [
        normalize_support_row(
            gene=gene,
            source=path,
            row=row,
            seed_keys=("strict_seed_passed_reads", "seed_passed_reads"),
            accepted_keys=("accepted_target_reads", "accepted_tp73_reads", "accepted_target_count"),
        )
        for row in read_tsv(path)
    ]


def merge_read_length_context(rows: list[dict[str, object]]) -> None:
    """Fill missing all-read length context from any same-SRR row that has it."""
    by_run: dict[str, dict[str, object]] = {}
    length_keys = [key for key in OUTPUT_FIELDS if key.startswith("all_")]
    for row in rows:
        run = str(row.get("run_accession") or "")
        if not run:
            continue
        current = by_run.get(run)
        current_score = sum(1 for key in length_keys if current and str(current.get(key) or ""))
        row_score = sum(1 for key in length_keys if str(row.get(key) or ""))
        if current is None or row_score > current_score:
            by_run[run] = row
    for row in rows:
        context = by_run.get(str(row.get("run_accession") or ""))
        if not context:
            continue
        for key in length_keys:
            if not str(row.get(key) or "") and str(context.get(key) or ""):
                row[key] = context[key]


def write_family_tsv(rows: list[dict[str, object]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=OUTPUT_FIELDS)
        writer.writeheader()
        writer.writerows(
            {field: row.get(field, "") for field in OUTPUT_FIELDS}
            for row in rows
        )


def validate_read_length_quantiles(rows: list[dict[str, object]], source_label: str) -> None:
    order = [
        "all_q0_bp",
        "all_q25_bp",
        "all_q50_bp",
        "all_q75_bp",
        "all_q90_bp",
        "all_q95_bp",
        "all_q99_bp",
        "all_q100_bp",
    ]
    for row in rows:
        previous_key = None
        previous_value = None
        for key in order:
            current = value(row, key)
            if current is None:
                continue
            if previous_value is not None and current < previous_value:
                print(
                    "WARNING: non-monotone all-read length quantiles in "
                    f"{source_label}: gene={row.get('gene')} run={row.get('run_accession')} "
                    f"{previous_key}={previous_value:g} > {key}={current:g}",
                    file=sys.stderr,
                )
            previous_key = key
            previous_value = current


def value(row: dict[str, object], key: str) -> float | None:
    if row is None or key not in row or row.get(key) is None:
        return None
    text = str(row.get(key)).strip()
    if not text:
        return None
    parsed = to_float(text, math.nan)
    return None if math.isnan(parsed) else parsed


def sample_label(
    run: str,
    rows_by_gene_run: dict[tuple[str, str], dict[str, object]],
    genes: list[str],
    label_column: str,
) -> str:
    for gene in genes:
        row = rows_by_gene_run.get((gene, run))
        if not row:
            continue
        label = str(row.get(label_column) or "").strip()
        if label:
            return label
        fallback = str(row.get("sample_id") or row.get("sample_name") or "").strip()
        if fallback:
            return fallback
    return run


def nice_linear_max(raw: float) -> float:
    if raw <= 0:
        return 1.0
    magnitude = 10 ** math.floor(math.log10(raw))
    scaled = raw / magnitude
    if scaled <= 1:
        nice = 1
    elif scaled <= 2:
        nice = 2
    elif scaled <= 5:
        nice = 5
    else:
        nice = 10
    return nice * magnitude


def format_bp(raw: float) -> str:
    if raw >= 1000:
        return f"{raw / 1000.0:.1f}k".replace(".0k", "k")
    return f"{raw:.0f}"


def format_metric(raw: float, metric: str) -> str:
    if metric == "raw":
        return f"{raw:.0f}"
    if raw >= 10:
        return f"{raw:.1f}"
    return f"{raw:.2f}"



def svg_text(
    x: float,
    y: float,
    text: object,
    size: int = 14,
    fill: str = "#111827",
    anchor: str = "start",
    weight: str = "400",
    extra: str = "",
) -> str:
    return (
        f'<text x="{x:.1f}" y="{y:.1f}" font-family="Avenir Next, Helvetica, Arial, sans-serif" '
        f'font-size="{size}" fill="{fill}" text-anchor="{anchor}" font-weight="{weight}" {extra}>'
        f"{escape(str(text))}</text>"
    )


def polyline(points: list[tuple[float, float]], color: str, dash: str = "") -> str:
    if not points:
        return ""
    coords = " ".join(f"{x:.1f},{y:.1f}" for x, y in points)
    dash_attr = f' stroke-dasharray="{dash}"' if dash else ""
    circles = "\n".join(
        f'<circle cx="{x:.1f}" cy="{y:.1f}" r="3.8" fill="{color}" stroke="white" stroke-width="1.4"/>'
        for x, y in points
    )
    return (
        f'<polyline points="{coords}" fill="none" stroke="{color}" '
        f'stroke-width="2.4" stroke-linejoin="round" stroke-linecap="round"{dash_attr}/>\n'
        f"{circles}"
    )


def render_svg(
    rows: list[dict[str, object]],
    genes: list[str],
    gene_colors: dict[str, str],
    output_tsv: Path,
    metric: str,
    label_column: str,
    title: str,
    support_length_stat: str,
    support_length_source: str,
    support_length_scale: str,
    support_length_display: str,
    support_length_genes: set[str],
    show_all_read_max: bool,
) -> str:
    rows_by_gene_run = {
        (str(row.get("gene")), str(row.get("run_accession"))): row
        for row in rows
        if row.get("gene") and row.get("run_accession")
    }
    runs = sorted({str(row.get("run_accession")) for row in rows if row.get("run_accession")}, key=sort_key_for_run)
    if not runs:
        raise SystemExit("No run accessions found for gene-family plot")

    n = len(runs)
    width = max(1420, 92 * n + 300)
    height = 950
    margin_left = 95
    margin_right = 130
    top_x = margin_left
    top_y = 102
    top_h = 300
    chart_w = width - margin_left - margin_right
    bottom_y = 535
    bottom_h = 280
    x_step = chart_w / max(n - 1, 1)
    x_positions = [top_x + idx * x_step for idx in range(n)]
    if n == 1:
        x_positions = [top_x + chart_w / 2.0]

    read_length_series = READ_LENGTH_SERIES if show_all_read_max else DEFAULT_READ_LENGTH_SERIES

    read_values = [
        number
        for run in runs
        for row in [next((rows_by_gene_run.get((gene, run)) for gene in genes if rows_by_gene_run.get((gene, run))), None)]
        if row
        for key, _label, _color, _dash in read_length_series
        if (number := value(row, key)) is not None and number > 0
    ]
    show_read_panel = bool(read_values)
    log_min = math.floor(math.log10(min(read_values))) if show_read_panel else 0
    log_max = math.ceil(math.log10(max(read_values))) if show_read_panel else 1
    if log_min == log_max:
        log_max += 1

    def y_read(raw: float) -> float:
        return top_y + top_h - ((math.log10(raw) - log_min) / (log_max - log_min)) * top_h

    metric_key = "seed_passed_per_million" if metric == "per_million" else "seed_passed_reads"
    metric_values = [
        value(row, metric_key) or 0.0
        for row in rows_by_gene_run.values()
    ]
    metric_max = nice_linear_max(max(metric_values) * 1.15 if metric_values else 1.0)
    support_length_key = {
        "max": "support_max_read_bp",
        "mean": "support_mean_read_bp",
    }[support_length_stat]
    support_length_label = {
        "max": "strict seed-passed max read length (bp)",
        "mean": "strict seed-passed mean read length (bp)",
    }[support_length_stat]
    if support_length_source == "any":
        support_length_label = support_length_label.replace("strict seed-passed", "gene-supporting")

    def support_source_allowed(row: dict[str, object]) -> bool:
        if support_length_source == "any":
            return True
        return str(row.get("support_length_source") or "") == "strict_seed_passed"

    def support_gene_allowed(gene: str) -> bool:
        return not support_length_genes or gene.upper() in support_length_genes

    support_length_values = [
        raw
        for row in rows_by_gene_run.values()
        if support_gene_allowed(str(row.get("gene") or ""))
        and support_source_allowed(row)
        and (raw := value(row, support_length_key)) is not None and raw > 0
    ]
    support_length_raw_max = max(support_length_values) if support_length_values else 0.0
    support_length_axis_max = (
        nice_linear_max(max(support_length_values) * 1.15)
        if support_length_values and support_length_scale == "linear"
        else 0.0
    )
    support_log_min = 0
    support_log_max = 1
    if support_length_values and support_length_scale == "log":
        support_log_min = math.floor(math.log10(min(support_length_values)))
        support_log_max = math.ceil(math.log10(max(support_length_values)))
        if support_log_min == support_log_max:
            support_log_max += 1
        support_length_label = support_length_label.replace("(bp)", "(bp, log scale)")

    def y_bar(raw: float) -> float:
        return bottom_y + bottom_h - (raw / metric_max) * bottom_h

    def y_support_len(raw: float) -> float:
        if support_length_raw_max <= 0:
            return bottom_y + bottom_h
        if support_length_scale == "log":
            return bottom_y + bottom_h - (
                (math.log10(raw) - support_log_min) / (support_log_max - support_log_min)
            ) * bottom_h
        return bottom_y + bottom_h - (raw / support_length_axis_max) * bottom_h

    metric_label = (
        "strict seed-passed reads per million total reads"
        if metric == "per_million"
        else "strict seed-passed reads"
    )

    parts: list[str] = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" role="img">',
        "<defs>",
        "<style>.axis{stroke:#334155;stroke-width:1.2}.grid{stroke:#e2e8f0;stroke-width:1}.note{fill:#64748b}</style>",
        "</defs>",
        f'<rect width="{width}" height="{height}" fill="#fbfaf7"/>',
        svg_text(38, 42, title, 24, weight="700"),
        svg_text(
            38,
            68,
            "Grouped bars use conservative seed-passed reads; raw counts are printed in each bar.",
            13,
            "#475569",
        ),
    ]

    parts.append(svg_text(top_x, top_y - 22, "Library read-length distributions", 17, weight="700"))
    if show_read_panel:
        for decade in range(log_min, log_max + 1):
            raw = 10**decade
            y = y_read(raw)
            parts.append(f'<line x1="{top_x}" y1="{y:.1f}" x2="{top_x + chart_w}" y2="{y:.1f}" class="grid"/>')
            parts.append(svg_text(top_x - 12, y + 4, format_bp(raw), 12, "#64748b", "end"))
        parts.append(f'<line x1="{top_x}" y1="{top_y}" x2="{top_x}" y2="{top_y + top_h}" class="axis"/>')
        parts.append(f'<line x1="{top_x}" y1="{top_y + top_h}" x2="{top_x + chart_w}" y2="{top_y + top_h}" class="axis"/>')
        parts.append(svg_text(26, top_y + top_h / 2, "read length (bp, log scale)", 12, "#475569", "middle", extra='transform="rotate(-90 26 %.1f)"' % (top_y + top_h / 2)))

        legend_x = top_x + chart_w - 430
        legend_y = top_y - 30
        for idx, (_key, label, color, dash) in enumerate(read_length_series):
            lx = legend_x + idx * 77
            parts.append(f'<line x1="{lx:.1f}" y1="{legend_y:.1f}" x2="{lx + 22:.1f}" y2="{legend_y:.1f}" stroke="{color}" stroke-width="2.5" stroke-dasharray="{dash}"/>')
            parts.append(svg_text(lx + 28, legend_y + 4, label, 12, "#334155"))

        for key, _label, color, dash in read_length_series:
            points = []
            for idx, run in enumerate(runs):
                row = next((rows_by_gene_run.get((gene, run)) for gene in genes if rows_by_gene_run.get((gene, run))), None)
                raw = value(row, key) if row else None
                if raw is not None and raw > 0:
                    points.append((x_positions[idx], y_read(raw)))
            parts.append(polyline(points, color, dash))
    else:
        parts.append(svg_text(top_x, top_y + 90, "No all-read length quantiles available in input tables.", 14, "#64748b"))

    family_label = " / ".join(genes)
    parts.append(svg_text(top_x, bottom_y - 24, f"{family_label} seed-passed support", 17, weight="700"))
    for tick in [0.0, 0.25, 0.5, 0.75, 1.0]:
        raw = metric_max * tick
        y = y_bar(raw)
        parts.append(f'<line x1="{top_x}" y1="{y:.1f}" x2="{top_x + chart_w}" y2="{y:.1f}" class="grid"/>')
        parts.append(svg_text(top_x - 12, y + 4, format_metric(raw, metric), 12, "#64748b", "end"))
    parts.append(f'<line x1="{top_x}" y1="{bottom_y}" x2="{top_x}" y2="{bottom_y + bottom_h}" class="axis"/>')
    if support_length_raw_max > 0:
        parts.append(f'<line x1="{top_x + chart_w}" y1="{bottom_y}" x2="{top_x + chart_w}" y2="{bottom_y + bottom_h}" class="axis"/>')
    parts.append(f'<line x1="{top_x}" y1="{bottom_y + bottom_h}" x2="{top_x + chart_w}" y2="{bottom_y + bottom_h}" class="axis"/>')
    parts.append(svg_text(28, bottom_y + bottom_h / 2, metric_label, 12, "#475569", "middle", extra='transform="rotate(-90 28 %.1f)"' % (bottom_y + bottom_h / 2)))
    if support_length_raw_max > 0:
        if support_length_scale == "log":
            for decade in range(support_log_min, support_log_max + 1):
                raw = 10**decade
                parts.append(svg_text(top_x + chart_w + 12, y_support_len(raw) + 4, format_bp(raw), 12, "#7c2d12", "start"))
        else:
            for tick in [0.25, 0.5, 0.75, 1.0]:
                raw = support_length_axis_max * tick
                parts.append(svg_text(top_x + chart_w + 12, y_support_len(raw) + 4, format_bp(raw), 12, "#7c2d12", "start"))
        parts.append(svg_text(top_x + chart_w + 84, bottom_y + bottom_h / 2, support_length_label, 12, "#7c2d12", "middle", extra='transform="rotate(90 %.1f %.1f)"' % (top_x + chart_w + 84, bottom_y + bottom_h / 2)))

    legend_gene_w = 86
    legend_extra_w = 140 if support_length_raw_max > 0 else 0
    legend_width = len(genes) * legend_gene_w + legend_extra_w + 22
    legend_x = top_x + chart_w - legend_width + 12
    legend_y = bottom_y - 43
    parts.append(f'<rect x="{legend_x - 12:.1f}" y="{legend_y - 17:.1f}" width="{legend_width}" height="32" rx="8" fill="#ffffff" stroke="#e2e8f0"/>')
    for idx, gene in enumerate(genes):
        lx = legend_x + idx * legend_gene_w
        parts.append(f'<rect x="{lx:.1f}" y="{legend_y - 9:.1f}" width="16" height="12" rx="2" fill="{gene_colors[gene]}"/>')
        parts.append(svg_text(lx + 22, legend_y + 2, gene, 12, "#334155", weight="700"))
    if support_length_raw_max > 0:
        lx = legend_x + len(genes) * legend_gene_w
        if support_length_display == "lines":
            parts.append(f'<line x1="{lx - 2:.1f}" y1="{legend_y - 3:.1f}" x2="{lx + 20:.1f}" y2="{legend_y - 3:.1f}" stroke="#111827" stroke-width="2.4"/>')
        parts.append(f'<circle cx="{lx + 9:.1f}" cy="{legend_y - 3:.1f}" r="3.8" fill="#111827" stroke="white" stroke-width="1"/>')
        parts.append(svg_text(lx + 28, legend_y + 2, f"{support_length_stat} read bp", 12, "#334155"))

    desired_group_w = max(58.0, len(genes) * 18.0 + max(0, len(genes) - 1) * 2.5)
    group_w = min(desired_group_w, chart_w / max(n, 1) * 0.72)
    bar_gap = 2.5
    bar_count = max(len(genes), 1)
    bar_w = max(4.0, (group_w - max(0, bar_count - 1) * bar_gap) / bar_count)
    for run_idx, run in enumerate(runs):
        base_x = x_positions[run_idx] - group_w / 2.0
        for gene_idx, gene in enumerate(genes):
            row = rows_by_gene_run.get((gene, run))
            metric_raw = value(row, metric_key) if row else 0.0
            count_raw = value(row, "seed_passed_reads") if row else 0.0
            x = base_x + gene_idx * (bar_w + bar_gap)
            y = y_bar(metric_raw or 0.0)
            h = bottom_y + bottom_h - y
            fill = gene_colors[gene]
            label_fill = fill
            parts.append(f'<rect x="{x:.1f}" y="{y:.1f}" width="{bar_w:.1f}" height="{max(h, 1.0):.1f}" rx="3" fill="{fill}" opacity="0.84">')
            parts.append(f'<title>{escape(gene)} {escape(run)}: {count_raw:.0f} seed-passed reads; {metric_raw or 0.0:.4g} {"per million" if metric == "per_million" else "raw"}</title>')
            parts.append("</rect>")
            if count_raw is not None and count_raw > 0:
                parts.append(svg_text(x + bar_w / 2.0, y - 5, f"{count_raw:.0f}", 10, label_fill, "middle", weight="700"))

    if support_length_raw_max > 0:
        for gene_idx, gene in enumerate(genes):
            if not support_gene_allowed(gene):
                continue
            points: list[tuple[float, float, str, float, str]] = []
            for idx, run in enumerate(runs):
                row = rows_by_gene_run.get((gene, run))
                raw = value(row, support_length_key) if row and support_source_allowed(row) else None
                if raw is None or raw <= 0:
                    continue
                source = str(row.get("support_length_source") or "support")
                base_x = x_positions[idx] - group_w / 2.0
                point_x = base_x + gene_idx * (bar_w + bar_gap) + bar_w / 2.0
                points.append((point_x, y_support_len(raw), run, raw, source))
            if support_length_display == "lines" and len(points) >= 2:
                coords = " ".join(f"{x:.1f},{y:.1f}" for x, y, _run, _raw, _source in points)
                parts.append(
                    f'<polyline points="{coords}" fill="none" stroke="{gene_colors[gene]}" '
                    'stroke-width="2.4" stroke-linejoin="round" stroke-linecap="round" opacity="0.95"/>'
                )
            for x, y, run, raw, source in points:
                parts.append(
                    f'<circle cx="{x:.1f}" cy="{y:.1f}" r="3.8" fill="{gene_colors[gene]}" '
                    'stroke="white" stroke-width="1.3">'
                )
                parts.append(
                    f'<title>{escape(gene)} {escape(run)} {escape(source)} '
                    f'{escape(support_length_stat)} read length: {raw:.0f} bp</title>'
                )
                parts.append("</circle>")

    for idx, run in enumerate(runs):
        x = x_positions[idx]
        label = sample_label(run, rows_by_gene_run, genes, label_column)
        parts.append(f'<line x1="{x:.1f}" y1="{top_y + top_h}" x2="{x:.1f}" y2="{top_y + top_h + 5}" stroke="#334155"/>')
        parts.append(f'<line x1="{x:.1f}" y1="{bottom_y + bottom_h}" x2="{x:.1f}" y2="{bottom_y + bottom_h + 5}" stroke="#334155"/>')
        parts.append(svg_text(x, bottom_y + bottom_h + 27, label, 12, "#334155", "end", extra=f'transform="rotate(-45 {x:.1f} {bottom_y + bottom_h + 27:.1f})"'))

    parts.append(svg_text(38, height - 42, f"Canonical source table: {output_tsv}", 10, "#64748b"))
    source_note = "strict seed-passed only" if support_length_source == "strict_seed_passed" else "any provenance, including accepted-target fallback"
    parts.append(svg_text(38, height - 24, f"accepted_target_count is retained in the TSV but not used as the main support scale; lower-panel symbols show {support_length_stat} support-read length ({source_note}).", 10, "#64748b"))
    parts.append("</svg>")
    return "\n".join(parts) + "\n"



def main() -> int:
    args = parse_args()
    output = Path(args.output)
    output_tsv = Path(args.output_tsv) if args.output_tsv else output.with_suffix(".tsv")

    rows = []
    input_genes: list[str] = []
    for value_text in args.canonical_summary:
        path = Path(value_text)
        if not path.exists():
            print(f"Input not found: {path}", file=sys.stderr)
            return 2
        canonical_rows = load_canonical_rows(path)
        for row in canonical_rows:
            gene = str(row.get("gene") or "").upper()
            if gene and gene not in input_genes:
                input_genes.append(gene)
        rows.extend(canonical_rows)
    for value_text in args.batch_report:
        gene, path = parse_gene_path(value_text, "--batch-report")
        if not path.exists():
            print(f"Input not found: {path}", file=sys.stderr)
            return 2
        input_genes.append(gene)
        rows.extend(load_batch_rows(gene, path, is_json=True))
    for value_text in args.batch_summary:
        gene, path = parse_gene_path(value_text, "--batch-summary")
        if not path.exists():
            print(f"Input not found: {path}", file=sys.stderr)
            return 2
        input_genes.append(gene)
        rows.extend(load_batch_rows(gene, path, is_json=False))
    for value_text in args.figure_source:
        gene, path = parse_gene_path(value_text, "--figure-source")
        if not path.exists():
            print(f"Input not found: {path}", file=sys.stderr)
            return 2
        input_genes.append(gene)
        rows.extend(load_figure_rows(gene, path))
    if not rows:
        print(
            "No input rows loaded; provide at least one --canonical-summary, "
            "--batch-report, --batch-summary, or --figure-source.",
            file=sys.stderr,
        )
        return 2

    if args.genes:
        genes = [gene.strip().upper() for gene in args.genes.split(",") if gene.strip()]
    else:
        genes = []
        for gene in input_genes:
            if gene not in genes:
                genes.append(gene)
    row_genes = {str(row.get("gene")) for row in rows if row.get("gene")}
    missing = [gene for gene in genes if gene not in row_genes]
    if missing:
        print(f"WARNING: requested gene(s) have no rows: {', '.join(missing)}", file=sys.stderr)
    extra = sorted(row_genes.difference(genes))
    genes.extend(extra)
    gene_colors = {
        gene: DEFAULT_PALETTE[idx % len(DEFAULT_PALETTE)]
        for idx, gene in enumerate(genes)
    }

    merge_read_length_context(rows)
    rows.sort(key=lambda row: (sort_key_for_run(str(row.get("run_accession") or "")), str(row.get("gene") or "")))
    validate_read_length_quantiles(rows, "gene-family canonical table")
    support_length_genes = {
        gene.strip().upper()
        for gene in (args.support_length_genes or "").split(",")
        if gene.strip()
    }

    write_family_tsv(rows, output_tsv)
    svg = render_svg(
        rows,
        genes,
        gene_colors,
        output_tsv,
        args.metric,
        args.label_column,
        args.title,
        args.support_length_stat,
        args.support_length_source,
        args.support_length_scale,
        args.support_length_display,
        support_length_genes,
        args.show_all_read_max,
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(svg, encoding="utf-8")
    print(f"Wrote {output_tsv}")
    print(f"Wrote {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
