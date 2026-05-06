#!/usr/bin/env python3
"""Render p53-family pancreas RNA-read support plots.

The plotter deliberately consumes already-produced GENtle report tables instead
of re-running any biology. It makes one canonical family TSV first, then renders
an SVG with library read-length context and grouped TP53/TP63/TP73 support bars.
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


GENES = ("TP53", "TP63", "TP73")
GENE_COLORS = {
    "TP53": "#2563eb",
    "TP63": "#ea580c",
    "TP73": "#16a34a",
}
READ_LENGTH_SERIES = [
    ("all_q25_bp", "q25", "#2563eb", ""),
    ("all_q50_bp", "q50", "#111827", ""),
    ("all_q75_bp", "q75", "#f59e0b", ""),
    ("all_q90_bp", "q90", "#7c3aed", "8,4"),
    ("all_q100_bp", "q100", "#dc2626", ""),
    ("all_mean_bp", "mean", "#059669", "6,4"),
]
OUTPUT_FIELDS = [
    "gene",
    "run_accession",
    "sample_id",
    "sample_name",
    "total_reads",
    "all_q0_bp",
    "all_q25_bp",
    "all_q50_bp",
    "all_q75_bp",
    "all_q90_bp",
    "all_q100_bp",
    "all_mean_bp",
    "seed_passed_reads",
    "seed_passed_per_million",
    "accepted_target_count",
    "accepted_target_per_million",
    "source",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Render a p53-family pancreas grouped-bar SVG from TP53 batch "
            "results plus TP63/TP73 figure-source TSVs."
        )
    )
    tp53 = parser.add_mutually_exclusive_group(required=True)
    tp53.add_argument(
        "--tp53-batch-report",
        help="TP53 out/batch_report.json from rna-reads batch-map.",
    )
    tp53.add_argument(
        "--tp53-batch-summary",
        help="TP53 out/batch_summary.tsv from rna-reads batch-map.",
    )
    parser.add_argument(
        "--tp63-figure-source",
        required=True,
        help="TP63 pancreas figure-source TSV from pancreas_gene_rna_screen.sh summarize.",
    )
    parser.add_argument(
        "--tp73-figure-source",
        required=True,
        help="TP73 pancreas figure-source TSV.",
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
        "--label-column",
        choices=["sample_id", "sample_name", "run_accession"],
        default="sample_id",
        help="Sample label shown on the x axis.",
    )
    parser.add_argument(
        "--title",
        default="p53-family seed-passed Nanopore cDNA support",
        help="SVG title.",
    )
    return parser.parse_args()


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
        "gene": gene,
        "run_accession": run,
        "sample_id": first_text(row, "sample_id"),
        "sample_name": first_text(row, "sample_name"),
        "total_reads": int(total) if float(total).is_integer() else total,
        "seed_passed_reads": int(seed) if float(seed).is_integer() else seed,
        "seed_passed_per_million": (seed * 1_000_000.0 / total) if total > 0 else 0.0,
        "accepted_target_count": int(accepted) if float(accepted).is_integer() else accepted,
        "accepted_target_per_million": (accepted * 1_000_000.0 / total) if total > 0 else 0.0,
        "source": str(source),
    }
    for key in (
        "all_q0_bp",
        "all_q25_bp",
        "all_q50_bp",
        "all_q75_bp",
        "all_q90_bp",
        "all_q100_bp",
        "all_mean_bp",
    ):
        out[key] = optional_number_text(row.get(key, ""))
    return out


def load_tp53_rows(path: Path, is_json: bool) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    if is_json:
        with path.open() as handle:
            report = json.load(handle)
        source_rows = report.get("rows")
        if not isinstance(source_rows, list):
            raise SystemExit(f"TP53 batch report has no rows[] array: {path}")
    else:
        source_rows = read_tsv(path)
    for row in source_rows:
        if not isinstance(row, dict):
            continue
        if str(row.get("gene") or "").strip().upper() not in {"", "TP53"}:
            continue
        rows.append(
            normalize_support_row(
                gene="TP53",
                source=path,
                row=row,
                seed_keys=("read_count_seed_passed", "seed_passed_reads"),
                accepted_keys=("accepted_target_count",),
            )
        )
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
        writer.writerows(rows)


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
    label_column: str,
) -> str:
    for gene in GENES:
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
    output_tsv: Path,
    metric: str,
    label_column: str,
    title: str,
) -> str:
    rows_by_gene_run = {
        (str(row.get("gene")), str(row.get("run_accession"))): row
        for row in rows
        if row.get("gene") and row.get("run_accession")
    }
    runs = sorted({str(row.get("run_accession")) for row in rows if row.get("run_accession")}, key=sort_key_for_run)
    if not runs:
        raise SystemExit("No run accessions found for p53-family plot")

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

    read_values = [
        number
        for run in runs
        for row in [next((rows_by_gene_run.get((gene, run)) for gene in GENES if rows_by_gene_run.get((gene, run))), None)]
        if row
        for key, _label, _color, _dash in READ_LENGTH_SERIES
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

    def y_bar(raw: float) -> float:
        return bottom_y + bottom_h - (raw / metric_max) * bottom_h

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

        legend_x = top_x + chart_w - 505
        legend_y = top_y - 30
        for idx, (_key, label, color, dash) in enumerate(READ_LENGTH_SERIES):
            lx = legend_x + idx * 77
            parts.append(f'<line x1="{lx:.1f}" y1="{legend_y:.1f}" x2="{lx + 22:.1f}" y2="{legend_y:.1f}" stroke="{color}" stroke-width="2.5" stroke-dasharray="{dash}"/>')
            parts.append(svg_text(lx + 28, legend_y + 4, label, 12, "#334155"))

        for key, _label, color, dash in READ_LENGTH_SERIES:
            points = []
            for idx, run in enumerate(runs):
                row = next((rows_by_gene_run.get((gene, run)) for gene in GENES if rows_by_gene_run.get((gene, run))), None)
                raw = value(row, key) if row else None
                if raw is not None and raw > 0:
                    points.append((x_positions[idx], y_read(raw)))
            parts.append(polyline(points, color, dash))
    else:
        parts.append(svg_text(top_x, top_y + 90, "No all-read length quantiles available in input tables.", 14, "#64748b"))

    parts.append(svg_text(top_x, bottom_y - 24, "TP53 / TP63 / TP73 seed-passed support", 17, weight="700"))
    for tick in [0.0, 0.25, 0.5, 0.75, 1.0]:
        raw = metric_max * tick
        y = y_bar(raw)
        parts.append(f'<line x1="{top_x}" y1="{y:.1f}" x2="{top_x + chart_w}" y2="{y:.1f}" class="grid"/>')
        parts.append(svg_text(top_x - 12, y + 4, format_metric(raw, metric), 12, "#64748b", "end"))
    parts.append(f'<line x1="{top_x}" y1="{bottom_y}" x2="{top_x}" y2="{bottom_y + bottom_h}" class="axis"/>')
    parts.append(f'<line x1="{top_x}" y1="{bottom_y + bottom_h}" x2="{top_x + chart_w}" y2="{bottom_y + bottom_h}" class="axis"/>')
    parts.append(svg_text(28, bottom_y + bottom_h / 2, metric_label, 12, "#475569", "middle", extra='transform="rotate(-90 28 %.1f)"' % (bottom_y + bottom_h / 2)))

    legend_x = top_x + chart_w - 265
    legend_y = bottom_y - 43
    parts.append(f'<rect x="{legend_x - 12:.1f}" y="{legend_y - 17:.1f}" width="275" height="32" rx="8" fill="#ffffff" stroke="#e2e8f0"/>')
    for idx, gene in enumerate(GENES):
        lx = legend_x + idx * 86
        parts.append(f'<rect x="{lx:.1f}" y="{legend_y - 9:.1f}" width="16" height="12" rx="2" fill="{GENE_COLORS[gene]}"/>')
        parts.append(svg_text(lx + 22, legend_y + 2, gene, 12, "#334155", weight="700"))

    group_w = min(58.0, chart_w / max(n, 1) * 0.58)
    bar_gap = 2.5
    bar_w = (group_w - 2 * bar_gap) / 3.0
    for run_idx, run in enumerate(runs):
        base_x = x_positions[run_idx] - group_w / 2.0
        for gene_idx, gene in enumerate(GENES):
            row = rows_by_gene_run.get((gene, run))
            metric_raw = value(row, metric_key) if row else 0.0
            count_raw = value(row, "seed_passed_reads") if row else 0.0
            x = base_x + gene_idx * (bar_w + bar_gap)
            y = y_bar(metric_raw or 0.0)
            h = bottom_y + bottom_h - y
            fill = GENE_COLORS[gene]
            label_fill = fill
            parts.append(f'<rect x="{x:.1f}" y="{y:.1f}" width="{bar_w:.1f}" height="{max(h, 1.0):.1f}" rx="3" fill="{fill}" opacity="0.84">')
            parts.append(f'<title>{escape(gene)} {escape(run)}: {count_raw:.0f} seed-passed reads; {metric_raw or 0.0:.4g} {"per million" if metric == "per_million" else "raw"}</title>')
            parts.append("</rect>")
            if count_raw is not None and count_raw > 0:
                parts.append(svg_text(x + bar_w / 2.0, y - 5, f"{count_raw:.0f}", 10, label_fill, "middle", weight="700"))

    for idx, run in enumerate(runs):
        x = x_positions[idx]
        label = sample_label(run, rows_by_gene_run, label_column)
        parts.append(f'<line x1="{x:.1f}" y1="{top_y + top_h}" x2="{x:.1f}" y2="{top_y + top_h + 5}" stroke="#334155"/>')
        parts.append(f'<line x1="{x:.1f}" y1="{bottom_y + bottom_h}" x2="{x:.1f}" y2="{bottom_y + bottom_h + 5}" stroke="#334155"/>')
        parts.append(svg_text(x, bottom_y + bottom_h + 27, label, 12, "#334155", "end", extra=f'transform="rotate(-45 {x:.1f} {bottom_y + bottom_h + 27:.1f})"'))

    parts.append(svg_text(38, height - 42, f"Canonical source table: {output_tsv}", 10, "#64748b"))
    parts.append(svg_text(38, height - 24, "accepted_target_count is retained in the TSV but deliberately not used as the main family-support scale.", 10, "#64748b"))
    parts.append("</svg>")
    return "\n".join(parts) + "\n"


def main() -> int:
    args = parse_args()
    output = Path(args.output)
    output_tsv = Path(args.output_tsv) if args.output_tsv else output.with_suffix(".tsv")
    tp53_path = Path(args.tp53_batch_report or args.tp53_batch_summary)
    tp63_path = Path(args.tp63_figure_source)
    tp73_path = Path(args.tp73_figure_source)
    for path in (tp53_path, tp63_path, tp73_path):
        if not path.exists():
            print(f"Input not found: {path}", file=sys.stderr)
            return 2

    rows = []
    rows.extend(load_tp53_rows(tp53_path, is_json=bool(args.tp53_batch_report)))
    rows.extend(load_figure_rows("TP63", tp63_path))
    rows.extend(load_figure_rows("TP73", tp73_path))
    merge_read_length_context(rows)
    rows.sort(key=lambda row: (sort_key_for_run(str(row.get("run_accession") or "")), str(row.get("gene") or "")))

    write_family_tsv(rows, output_tsv)
    svg = render_svg(rows, output_tsv, args.metric, args.label_column, args.title)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(svg, encoding="utf-8")
    print(f"Wrote {output_tsv}")
    print(f"Wrote {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
