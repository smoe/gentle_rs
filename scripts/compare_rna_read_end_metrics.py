#!/usr/bin/env python3
"""Compare observed RNA-read length/end metrics across samples.

Adapted from Glen / smoe-bot commit a1cf936d0cbccbbc6f1e2cbed8474006d740bc25.
Biological projection (for example, transcript-bound 5' or 3' UTR reach) must
happen upstream; this script never infers transcript coordinates from lengths.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import html
import json
import math
import re
from collections import defaultdict
from contextlib import contextmanager
import io
from itertools import combinations
import os
from pathlib import Path
import tempfile


SCHEMA = "gentle.rna_read_metric_comparison.v1"
HEADER = ["sample_id", "metric", "value_bp", "observation_status"]
STATUSES = {"observed", "left_censored", "right_censored", "missing"}
TOKEN = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]{0,63}$")
MAX_VALUE_BP = 2**53 - 1


def quantile(values: list[int], probability: float) -> float:
    """Return the deterministic type-7 sample quantile."""
    if not values:
        raise ValueError("quantile requires at least one observed value")
    position = (len(values) - 1) * probability
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return float(values[lower])
    fraction = position - lower
    return values[lower] * (1.0 - fraction) + values[upper] * fraction


def parse_input(path: Path, positive_metrics: frozenset[str] = frozenset()):
    raw = path.read_bytes()
    try:
        text = raw.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise ValueError("input TSV must be UTF-8") from exc
    reader = csv.DictReader(io.StringIO(text, newline=""), delimiter="\t", strict=True)
    observed = defaultdict(list)
    counts = defaultdict(lambda: {status: 0 for status in sorted(STATUSES)})
    try:
        if reader.fieldnames != HEADER:
            raise ValueError(f"input TSV requires exactly: {','.join(HEADER)}")
        for row in reader:
            line_number = reader.line_num
            if None in row or any(row[field] is None for field in HEADER):
                raise ValueError(f"line {line_number}: expected exactly four TSV fields")
            sample = row["sample_id"].strip()
            metric = row["metric"].strip()
            status = row["observation_status"].strip()
            value = row["value_bp"].strip()
            if not TOKEN.fullmatch(sample) or not TOKEN.fullmatch(metric):
                raise ValueError(f"line {line_number}: invalid sample_id or metric token")
            if status not in STATUSES:
                raise ValueError(f"line {line_number}: unsupported observation_status {status!r}")
            key = (sample, metric)
            if status == "observed":
                if not re.fullmatch(r"[0-9]+", value):
                    raise ValueError(f"line {line_number}: observed value_bp must be a nonnegative integer")
                parsed = int(value)
                if parsed > MAX_VALUE_BP:
                    raise ValueError(f"line {line_number}: value_bp exceeds the exact integer limit {MAX_VALUE_BP}")
                if parsed == 0 and metric in positive_metrics:
                    raise ValueError(f"line {line_number}: {metric} requires a positive observed value")
                observed[key].append(parsed)
            elif value:
                raise ValueError(
                    f"line {line_number}: censored/missing values must be blank, never encoded as zero"
                )
            counts[key][status] += 1
    except csv.Error as exc:
        raise ValueError(f"line {reader.line_num}: invalid TSV quoting: {exc}") from exc
    if not counts:
        raise ValueError("input TSV contains no observations")
    for values in observed.values():
        values.sort()
    return observed, counts, raw


@contextmanager
def staged_output(path: Path, *, newline: str = ""):
    """Close/sync a private staging file before atomic no-clobber publication."""
    descriptor, name = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".partial", dir=path.parent)
    partial = Path(name)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline=newline) as handle:
            yield handle
            handle.flush()
            os.fsync(handle.fileno())
        os.link(partial, path)
    finally:
        partial.unlink(missing_ok=True)


def write_tsv(path: Path, header: list[str], rows: list[list[object]]) -> None:
    with staged_output(path) as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def svg_for_metric(metric: str, pairs: list[dict[str, object]], output: Path) -> None:
    columns, panel_w, panel_h = 3, 260, 230
    rows = math.ceil(len(pairs) / columns)
    width, height = columns * panel_w + 30, rows * panel_h + 60
    all_values = [float(value) for pair in pairs for point in pair["points"] for value in point[1:]]
    low, high = min(all_values), max(all_values)
    if low == high:
        high = low + 1.0
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        f'<text x="15" y="25" font-family="sans-serif" font-size="17" font-weight="700">All-pairs Q-Q: {html.escape(metric)}</text>',
        '<text x="15" y="45" font-family="sans-serif" font-size="11" fill="#475569">Equal-scale axes; gray diagonal is identity. Censored/missing observations are excluded, not zero.</text>',
    ]
    for index, pair in enumerate(pairs):
        row, column = divmod(index, columns)
        x0, y0, side = 45 + column * panel_w, 78 + row * panel_h, 155
        parts.extend([
            f'<text x="{x0}" y="{y0 - 12}" font-family="sans-serif" font-size="12">{html.escape(str(pair["left_sample"]))} vs {html.escape(str(pair["right_sample"]))}</text>',
            f'<line x1="{x0}" y1="{y0 + side}" x2="{x0 + side}" y2="{y0 + side}" stroke="#111827"/>',
            f'<line x1="{x0}" y1="{y0}" x2="{x0}" y2="{y0 + side}" stroke="#111827"/>',
            f'<line x1="{x0}" y1="{y0 + side}" x2="{x0 + side}" y2="{y0}" stroke="#94a3b8"/>',
        ])
        for _probability, left, right in pair["points"]:
            x = x0 + (float(left) - low) / (high - low) * side
            y = y0 + side - (float(right) - low) / (high - low) * side
            parts.append(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="1.45" fill="#2563eb"/>')
        parts.append(f'<text x="{x0}" y="{y0 + side + 17}" font-family="sans-serif" font-size="9" fill="#475569">{low:.0f}–{high:.0f} bp on both axes</text>')
    parts.append("</svg>\n")
    with staged_output(output) as handle:
        handle.write("\n".join(parts))


def run(input_path: Path, output_dir: Path, quantile_points: int,
        positive_metrics: tuple[str, ...] = ()) -> None:
    if quantile_points < 11 or quantile_points > 1001:
        raise ValueError("--quantile-points must be between 11 and 1001")
    positive_metrics = tuple(sorted(set(positive_metrics)))
    if any(not TOKEN.fullmatch(metric) for metric in positive_metrics):
        raise ValueError("--positive-metric requires a valid metric token")
    observed, counts, raw = parse_input(input_path, frozenset(positive_metrics))
    samples = sorted({key[0] for key in counts})
    metrics = sorted({key[1] for key in counts})
    if set(positive_metrics) - set(metrics):
        raise ValueError("--positive-metric names a metric absent from the input")
    if len(samples) < 2:
        raise ValueError("at least two samples are required")
    for sample in samples:
        for metric in metrics:
            if len(observed.get((sample, metric), [])) < 2:
                raise ValueError(f"{sample}/{metric} needs at least two observed values")
    summary_rows = []
    for metric in metrics:
        for sample in samples:
            values = observed[(sample, metric)]
            state = counts[(sample, metric)]
            summary_rows.append([
                sample, metric, len(values), state["left_censored"], state["right_censored"], state["missing"],
                *(f"{quantile(values, p):.6f}" for p in (0, 0.1, 0.25, 0.5, 0.75, 0.9, 1)),
                f"{sum(values) / len(values):.6f}",
            ])
    if not output_dir.name or ".." in output_dir.parts:
        raise ValueError("output must name a fresh directory without '..' components")
    parent = output_dir.parent.resolve(strict=True)
    if not parent.is_dir():
        raise ValueError("output parent must be an existing directory")
    output_dir = parent / output_dir.name
    output_dir.mkdir(mode=0o700)
    write_tsv(output_dir / "summary.tsv",
        ["sample_id", "metric", "observed_count", "left_censored_count", "right_censored_count", "missing_count", "q0_bp", "q10_bp", "q25_bp", "q50_bp", "q75_bp", "q90_bp", "q100_bp", "mean_bp"],
        summary_rows)

    qq_rows, report_pairs = [], []
    probabilities = [index / (quantile_points - 1) for index in range(quantile_points)]
    for metric in metrics:
        metric_pairs = []
        for left_sample, right_sample in combinations(samples, 2):
            left_values, right_values = observed[(left_sample, metric)], observed[(right_sample, metric)]
            points = []
            for probability in probabilities:
                left, right = quantile(left_values, probability), quantile(right_values, probability)
                points.append([probability, left, right])
                qq_rows.append([metric, left_sample, right_sample, f"{probability:.6f}", f"{left:.6f}", f"{right:.6f}", f"{right - left:.6f}"])
            pair = {"metric": metric, "left_sample": left_sample, "right_sample": right_sample, "points": points}
            metric_pairs.append(pair)
            report_pairs.append(pair)
        svg_for_metric(metric, metric_pairs, output_dir / f"qq_{metric}.svg")
    write_tsv(output_dir / "qq_pairs.tsv",
        ["metric", "left_sample", "right_sample", "probability", "left_value_bp", "right_value_bp", "right_minus_left_bp"], qq_rows)
    report = {
        "schema": SCHEMA,
        "source": {"name": input_path.name, "sha256": hashlib.sha256(raw).hexdigest()},
        "quantile_method": "type_7_linear",
        "quantile_points": quantile_points,
        "value_policy": {"default": "nonnegative_integer_bp", "positive_metrics": list(positive_metrics)},
        "samples": samples,
        "metrics": metrics,
        "counts": [{"sample_id": sample, "metric": metric, **counts[(sample, metric)]} for metric in metrics for sample in samples],
        "pair_count_per_metric": len(list(combinations(samples, 2))),
        "warnings": [
            "Descriptive comparison only; no biological replicate or significance claim is inferred.",
            "Censored and missing observations are counted but excluded, never replaced with zero.",
            "UTR metrics must be projected upstream against an explicit transcript model and strand.",
            "Library preparation, read truncation, mapping and transcript-selection effects remain possible explanations.",
        ],
        "pairs": report_pairs,
    }
    with staged_output(output_dir / "report.json") as handle:
        handle.write(json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_tsv", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--quantile-points", type=int, default=101)
    parser.add_argument("--positive-metric", action="append", default=[], metavar="METRIC",
                        help="require positive observed values for this metric (repeatable; no name-based inference)")
    args = parser.parse_args()
    try:
        run(args.input_tsv, args.output_dir, args.quantile_points, tuple(args.positive_metric))
    except (OSError, ValueError) as exc:
        parser.error(str(exc))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
