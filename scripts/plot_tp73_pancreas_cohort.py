#!/usr/bin/env python3
"""Render TP73 pancreas cohort summary plots from the figure-source TSV.

The script intentionally writes SVG directly instead of depending on matplotlib.
That keeps it usable on minimal cluster installs where only Python is present.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import sys
from pathlib import Path
from xml.sax.saxutils import escape


READ_LENGTH_SERIES = [
    ("all_q0_bp", "q0", "#94a3b8", ""),
    ("all_q25_bp", "q25", "#2563eb", ""),
    ("all_q50_bp", "q50", "#111827", ""),
    ("all_q75_bp", "q75", "#f59e0b", ""),
    ("all_q100_bp", "q100", "#dc2626", ""),
    ("all_mean_bp", "mean", "#059669", "6,4"),
]


def parse_args() -> argparse.Namespace:
    cohort_root = os.environ.get("COHORT_ROOT")
    default_input = (
        Path(cohort_root) / "figures" / "tp73_pancreas_figure_source.tsv"
        if cohort_root
        else None
    )
    parser = argparse.ArgumentParser(
        description=(
            "Render TP73 pancreas cohort read-length and TP73-support SVG plots "
            "from tp73_pancreas_figure_source.tsv."
        )
    )
    parser.add_argument(
        "input_tsv",
        nargs="?",
        default=str(default_input) if default_input else None,
        help=(
            "Input TSV. Defaults to "
            "$COHORT_ROOT/figures/tp73_pancreas_figure_source.tsv when COHORT_ROOT is set."
        ),
    )
    parser.add_argument(
        "-o",
        "--output",
        default=None,
        help="Output SVG path. Defaults to INPUT.with_suffix('.svg').",
    )
    parser.add_argument(
        "--label-column",
        default="sample_id",
        choices=["sample_id", "sample_name", "run_accession"],
        help="Column used for x-axis sample labels.",
    )
    parser.add_argument(
        "--title",
        default="TP73 support across pancreatic cancer Nanopore cDNA samples",
        help="Figure title.",
    )
    return parser.parse_args()


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise SystemExit(f"No rows found in {path}")
    return rows


def number(row: dict[str, str], key: str) -> float | None:
    raw = (row.get(key) or "").replace(",", "").strip()
    if raw in {"", "null", "NA", "NaN", "nan"}:
        return None
    try:
        value = float(raw)
    except ValueError:
        return None
    if math.isnan(value) or math.isinf(value):
        return None
    return value


def sample_label(row: dict[str, str], label_column: str, index: int) -> str:
    label = (row.get(label_column) or "").strip()
    if not label:
        label = (row.get("sample_id") or row.get("run_accession") or f"sample_{index + 1}").strip()
    if label == row.get("run_accession", "") and row.get("sample_name", "").strip():
        label = row["run_accession"].strip()
    return label


def nice_linear_max(value: float) -> float:
    if value <= 0:
        return 1.0
    magnitude = 10 ** math.floor(math.log10(value))
    scaled = value / magnitude
    if scaled <= 1:
        nice = 1
    elif scaled <= 2:
        nice = 2
    elif scaled <= 5:
        nice = 5
    else:
        nice = 10
    return nice * magnitude


def format_bp(value: float) -> str:
    if value >= 1000:
        rounded = value / 1000.0
        return f"{rounded:.1f}k".replace(".0k", "k")
    return f"{value:.0f}"


def format_count(value: float) -> str:
    if value >= 1000:
        rounded = value / 1000.0
        return f"{rounded:.1f}k".replace(".0k", "k")
    return f"{value:.0f}"


def svg_text(
    x: float,
    y: float,
    text: str,
    size: int = 14,
    fill: str = "#111827",
    anchor: str = "start",
    weight: str = "400",
    extra: str = "",
) -> str:
    return (
        f'<text x="{x:.1f}" y="{y:.1f}" font-family="Avenir Next, Helvetica, Arial, sans-serif" '
        f'font-size="{size}" fill="{fill}" text-anchor="{anchor}" font-weight="{weight}" {extra}>'
        f"{escape(text)}</text>"
    )


def polyline(points: list[tuple[float, float]], color: str, dash: str = "") -> str:
    if not points:
        return ""
    coords = " ".join(f"{x:.1f},{y:.1f}" for x, y in points)
    dash_attr = f' stroke-dasharray="{dash}"' if dash else ""
    circles = "\n".join(
        f'<circle cx="{x:.1f}" cy="{y:.1f}" r="4.0" fill="{color}" stroke="white" stroke-width="1.5"/>'
        for x, y in points
    )
    return (
        f'<polyline points="{coords}" fill="none" stroke="{color}" '
        f'stroke-width="2.5" stroke-linejoin="round" stroke-linecap="round"{dash_attr}/>\n'
        f"{circles}"
    )


def render_svg(rows: list[dict[str, str]], args: argparse.Namespace, source_path: Path) -> str:
    labels = [sample_label(row, args.label_column, idx) for idx, row in enumerate(rows)]
    n = len(rows)
    width = max(1280, 90 * n + 260)
    height = 900
    margin_left = 95
    margin_right = 130
    title_y = 42
    top_x = margin_left
    top_y = 95
    chart_w = width - margin_left - margin_right
    top_h = 330
    bottom_y = 535
    bottom_h = 260
    x_step = chart_w / max(n - 1, 1)
    x_positions = [top_x + idx * x_step for idx in range(n)]
    if n == 1:
        x_positions = [top_x + chart_w / 2.0]

    read_values = [
        value
        for row in rows
        for key, _label, _color, _dash in READ_LENGTH_SERIES
        if (value := number(row, key)) is not None and value > 0
    ]
    if not read_values:
        raise SystemExit("No read-length values found in input TSV")
    log_min = math.floor(math.log10(min(read_values)))
    log_max = math.ceil(math.log10(max(read_values)))
    if log_min == log_max:
        log_max += 1

    def y_read(value: float) -> float:
        return top_y + top_h - ((math.log10(value) - log_min) / (log_max - log_min)) * top_h

    tp73_counts = [
        number(row, "accepted_tp73_reads")
        if number(row, "accepted_tp73_reads") is not None
        else number(row, "strict_seed_passed_reads") or 0.0
        for row in rows
    ]
    tp73_max_lengths = [number(row, "accepted_tp73_max_read_bp") for row in rows]
    count_max = nice_linear_max(max(tp73_counts) * 1.15 if tp73_counts else 1.0)
    length_positive = [value for value in tp73_max_lengths if value is not None and value > 0]
    length_max = nice_linear_max(max(length_positive) * 1.15 if length_positive else 1.0)

    def y_count(value: float) -> float:
        return bottom_y + bottom_h - (value / count_max) * bottom_h

    def y_tp73_len(value: float) -> float:
        return bottom_y + bottom_h - (value / length_max) * bottom_h

    parts: list[str] = []
    parts.append('<?xml version="1.0" encoding="UTF-8"?>')
    parts.append(
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}" role="img">'
    )
    parts.append("<defs>")
    parts.append(
        "<style>"
        ".axis{stroke:#334155;stroke-width:1.2}.grid{stroke:#e2e8f0;stroke-width:1}.soft{fill:#64748b}"
        "</style>"
    )
    parts.append("</defs>")
    parts.append(f'<rect width="{width}" height="{height}" fill="#fbfaf7"/>')
    parts.append(svg_text(38, title_y, args.title, 24, "#111827", weight="700"))
    parts.append(
        svg_text(
            38,
            title_y + 26,
            "Read-length distribution uses all reads; TP73 bars use accepted target reads, not retained phase-2 working rows.",
            13,
            "#475569",
        )
    )

    # Top chart frame and grid.
    parts.append(svg_text(top_x, top_y - 20, "Library read-length distributions", 17, "#111827", weight="700"))
    for decade in range(log_min, log_max + 1):
        value = 10**decade
        y = y_read(value)
        parts.append(f'<line x1="{top_x}" y1="{y:.1f}" x2="{top_x + chart_w}" y2="{y:.1f}" class="grid"/>')
        parts.append(svg_text(top_x - 12, y + 4, format_bp(value), 12, "#64748b", "end"))
    parts.append(f'<line x1="{top_x}" y1="{top_y}" x2="{top_x}" y2="{top_y + top_h}" class="axis"/>')
    parts.append(f'<line x1="{top_x}" y1="{top_y + top_h}" x2="{top_x + chart_w}" y2="{top_y + top_h}" class="axis"/>')
    parts.append(svg_text(26, top_y + top_h / 2, "read length (bp, log scale)", 12, "#475569", "middle", extra='transform="rotate(-90 26 %.1f)"' % (top_y + top_h / 2)))

    legend_x = top_x + chart_w - 470
    legend_y = top_y - 28
    for idx, (_key, label, color, dash) in enumerate(READ_LENGTH_SERIES):
        lx = legend_x + idx * 76
        parts.append(f'<line x1="{lx}" y1="{legend_y}" x2="{lx + 22}" y2="{legend_y}" stroke="{color}" stroke-width="2.5" stroke-dasharray="{dash}"/>')
        parts.append(svg_text(lx + 28, legend_y + 4, label, 12, "#334155"))

    for key, _label, color, dash in READ_LENGTH_SERIES:
        points = []
        for idx, row in enumerate(rows):
            value = number(row, key)
            if value is not None and value > 0:
                points.append((x_positions[idx], y_read(value)))
        parts.append(polyline(points, color, dash))

    # Bottom chart frame and axes.
    parts.append(svg_text(top_x, bottom_y - 24, "TP73-positive reads and maximum TP73 read length", 17, "#111827", weight="700"))
    for tick in [0.0, 0.25, 0.5, 0.75, 1.0]:
        value = count_max * tick
        y = y_count(value)
        parts.append(f'<line x1="{top_x}" y1="{y:.1f}" x2="{top_x + chart_w}" y2="{y:.1f}" class="grid"/>')
        parts.append(svg_text(top_x - 12, y + 4, format_count(value), 12, "#64748b", "end"))
        len_value = length_max * tick
        parts.append(svg_text(top_x + chart_w + 12, y + 4, format_bp(len_value), 12, "#c2410c", "start"))
    parts.append(f'<line x1="{top_x}" y1="{bottom_y}" x2="{top_x}" y2="{bottom_y + bottom_h}" class="axis"/>')
    parts.append(f'<line x1="{top_x + chart_w}" y1="{bottom_y}" x2="{top_x + chart_w}" y2="{bottom_y + bottom_h}" class="axis"/>')
    parts.append(f'<line x1="{top_x}" y1="{bottom_y + bottom_h}" x2="{top_x + chart_w}" y2="{bottom_y + bottom_h}" class="axis"/>')
    parts.append(svg_text(28, bottom_y + bottom_h / 2, "accepted TP73 reads", 12, "#475569", "middle", extra='transform="rotate(-90 28 %.1f)"' % (bottom_y + bottom_h / 2)))
    parts.append(svg_text(top_x + chart_w + 86, bottom_y + bottom_h / 2, "max TP73 read length (bp)", 12, "#c2410c", "middle", extra='transform="rotate(90 %.1f %.1f)"' % (top_x + chart_w + 86, bottom_y + bottom_h / 2)))

    bar_w = min(34.0, chart_w / max(n, 1) * 0.48)
    for idx, count in enumerate(tp73_counts):
        x = x_positions[idx] - bar_w / 2
        y = y_count(count)
        h = bottom_y + bottom_h - y
        parts.append(f'<rect x="{x:.1f}" y="{y:.1f}" width="{bar_w:.1f}" height="{max(h, 1):.1f}" rx="3" fill="#2563eb" opacity="0.82"/>')
        parts.append(svg_text(x_positions[idx], y - 6, format_count(count), 11, "#1e3a8a", "middle"))

    length_points = [
        (x_positions[idx], y_tp73_len(value))
        for idx, value in enumerate(tp73_max_lengths)
        if value is not None and value > 0
    ]
    parts.append(polyline(length_points, "#ea580c", ""))
    parts.append(f'<rect x="{top_x + chart_w - 245:.1f}" y="{bottom_y - 45:.1f}" width="236" height="30" rx="8" fill="#ffffff" stroke="#e2e8f0"/>')
    parts.append(f'<rect x="{top_x + chart_w - 232:.1f}" y="{bottom_y - 35:.1f}" width="18" height="12" rx="2" fill="#2563eb" opacity="0.82"/>')
    parts.append(svg_text(top_x + chart_w - 207, bottom_y - 25, "accepted TP73 reads", 12, "#334155"))
    parts.append(f'<line x1="{top_x + chart_w - 84:.1f}" y1="{bottom_y - 29:.1f}" x2="{top_x + chart_w - 58:.1f}" y2="{bottom_y - 29:.1f}" stroke="#ea580c" stroke-width="2.5"/>')
    parts.append(f'<circle cx="{top_x + chart_w - 71:.1f}" cy="{bottom_y - 29:.1f}" r="4" fill="#ea580c" stroke="white" stroke-width="1.5"/>')
    parts.append(svg_text(top_x + chart_w - 50, bottom_y - 25, "max bp", 12, "#334155"))

    # Shared x labels.
    for idx, label in enumerate(labels):
        x = x_positions[idx]
        parts.append(f'<line x1="{x:.1f}" y1="{top_y + top_h}" x2="{x:.1f}" y2="{top_y + top_h + 5}" stroke="#334155"/>')
        parts.append(f'<line x1="{x:.1f}" y1="{bottom_y + bottom_h}" x2="{x:.1f}" y2="{bottom_y + bottom_h + 5}" stroke="#334155"/>')
        parts.append(
            svg_text(
                x,
                bottom_y + bottom_h + 26,
                label,
                12,
                "#334155",
                "end",
                extra=f'transform="rotate(-45 {x:.1f} {bottom_y + bottom_h + 26:.1f})"',
            )
        )

    parts.append(
        svg_text(
            38,
            height - 34,
            f"Source: {source_path} | q0/q25/q50/q75/q100/mean from all reads; TP73 max from accepted target reads.",
            11,
            "#64748b",
        )
    )
    parts.append("</svg>")
    return "\n".join(parts) + "\n"


def main() -> int:
    args = parse_args()
    if not args.input_tsv:
        print(
            "No input TSV provided and COHORT_ROOT is not set. "
            "Pass the TSV path explicitly.",
            file=sys.stderr,
        )
        return 2
    input_path = Path(args.input_tsv)
    if not input_path.exists():
        print(f"Input TSV not found: {input_path}", file=sys.stderr)
        return 2
    output_path = Path(args.output) if args.output else input_path.with_suffix(".svg")
    rows = read_rows(input_path)
    svg = render_svg(rows, args, input_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(svg, encoding="utf-8")
    print(f"Wrote {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
