#!/usr/bin/env python3
"""Rank gene-wise read-end shifts using explicit like-library contrasts.

Exactly one contrast is primary. Differently prepared libraries are retained as
directional counterchecks and are never pooled with the primary comparison.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import html
import json
import statistics
from pathlib import Path

from compare_rna_read_end_metrics import staged_output
from rank_rna_read_end_genes import METRICS, fail, load


SCHEMA = "gentle.rna_read_end_gene_contrasts.v1"
CONTRAST_HEADER = (
    "contrast_id", "focus_sample", "comparison_sample",
    "library_preparation", "evidence_role",
)
SUMMARY_FIELDS = (
    "gene_id", "metric", "focus_label", "primary_contrast_id",
    "primary_focus_sample", "primary_comparison_sample",
    "primary_focus_observed_n", "primary_focus_right_censored_n",
    "primary_focus_missing_n", "primary_comparison_observed_n",
    "primary_comparison_right_censored_n", "primary_comparison_missing_n",
    "primary_delta_median_bp", "directional_contrasts_eligible",
    "directional_deltas_bp_json", "directional_concordant_count",
    "ranking_status",
)


def read_contrasts(path: Path, samples: list[str]) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != CONTRAST_HEADER:
            raise fail(f"unexpected contrast-sheet columns in {path}")
        rows = list(reader)
    if not rows:
        raise fail("contrast sheet is empty")
    identifiers: set[str] = set()
    primary = 0
    used_samples: set[str] = set()
    for line_number, row in enumerate(rows, 2):
        identifier = row["contrast_id"]
        if not identifier or identifier in identifiers:
            raise fail(f"empty or duplicate contrast_id at {path}:{line_number}")
        identifiers.add(identifier)
        if row["focus_sample"] == row["comparison_sample"]:
            raise fail(f"self-comparison at {path}:{line_number}")
        if row["focus_sample"] not in samples or row["comparison_sample"] not in samples:
            raise fail(f"unknown sample at {path}:{line_number}")
        if row["evidence_role"] not in {"primary", "directional_countercheck"}:
            raise fail(f"unsupported evidence role at {path}:{line_number}")
        if not row["library_preparation"]:
            raise fail(f"missing library preparation at {path}:{line_number}")
        primary += row["evidence_role"] == "primary"
        for sample in (row["focus_sample"], row["comparison_sample"]):
            if sample in used_samples:
                raise fail(f"sample {sample!r} is reused across contrasts")
            used_samples.add(sample)
    if primary != 1:
        raise fail("contrast sheet must contain exactly one primary contrast")
    return rows


def eligible(values, counts, gene, sample, metric, minimum, max_censored):
    key = (gene, sample, metric)
    observed = values[key]
    row_counts = counts[key]
    measured = row_counts["observed"] + row_counts["right_censored"]
    fraction = row_counts["right_censored"] / measured if measured else 0.0
    return observed, row_counts, len(observed) >= minimum and fraction <= max_censored


def summarize(paths: list[Path], contrast_path: Path, focus_label: str, minimum: int, max_censored: float):
    values, counts, samples = load(paths)
    contrasts = read_contrasts(contrast_path, samples)
    primary = next(row for row in contrasts if row["evidence_role"] == "primary")
    directional = [row for row in contrasts if row["evidence_role"] == "directional_countercheck"]
    genes = sorted({key[0] for key in counts})
    rows = []
    for gene in genes:
        for metric in METRICS:
            focus_values, focus_counts, focus_ok = eligible(
                values, counts, gene, primary["focus_sample"], metric, minimum, max_censored
            )
            comparison_values, comparison_counts, comparison_ok = eligible(
                values, counts, gene, primary["comparison_sample"], metric, minimum, max_censored
            )
            primary_delta = (
                statistics.median(focus_values) - statistics.median(comparison_values)
                if focus_ok and comparison_ok else None
            )
            directional_deltas = {}
            for contrast in directional:
                left, _left_counts, left_ok = eligible(
                    values, counts, gene, contrast["focus_sample"], metric, minimum, max_censored
                )
                right, _right_counts, right_ok = eligible(
                    values, counts, gene, contrast["comparison_sample"], metric, minimum, max_censored
                )
                if left_ok and right_ok:
                    directional_deltas[contrast["contrast_id"]] = statistics.median(left) - statistics.median(right)
            concordant = 0
            if primary_delta is not None and primary_delta != 0:
                concordant = sum(
                    delta != 0 and (delta > 0) == (primary_delta > 0)
                    for delta in directional_deltas.values()
                )
            rows.append({
                "gene_id": gene,
                "metric": metric,
                "focus_label": focus_label,
                "primary_contrast_id": primary["contrast_id"],
                "primary_focus_sample": primary["focus_sample"],
                "primary_comparison_sample": primary["comparison_sample"],
                "primary_focus_observed_n": len(focus_values),
                "primary_focus_right_censored_n": focus_counts["right_censored"],
                "primary_focus_missing_n": focus_counts["missing"],
                "primary_comparison_observed_n": len(comparison_values),
                "primary_comparison_right_censored_n": comparison_counts["right_censored"],
                "primary_comparison_missing_n": comparison_counts["missing"],
                "primary_delta_median_bp": "" if primary_delta is None else primary_delta,
                "directional_contrasts_eligible": ",".join(directional_deltas),
                "directional_deltas_bp_json": json.dumps(directional_deltas, sort_keys=True, separators=(",", ":")),
                "directional_concordant_count": concordant,
                "ranking_status": "eligible" if primary_delta is not None else "primary_insufficient_or_censored",
            })
    return rows, samples, contrasts


def svg(rows: list[dict], focus_label: str, top: int) -> str:
    candidates = [row for row in rows if row["metric"] == "utr3_reach_bp" and row["ranking_status"] == "eligible"]
    selected = sorted(candidates, key=lambda row: (-abs(float(row["primary_delta_median_bp"])), row["gene_id"]))[:top]
    genes = [row["gene_id"] for row in selected]
    lookup = {(row["gene_id"], row["metric"]): row for row in rows}
    width, left, panel, gap = 1180, 190, 290, 55
    height = max(235, 135 + 34 * len(genes))
    colors = {"read_length_bp": "#6c7a89", "utr5_reach_bp": "#2a9d8f", "utr3_reach_bp": "#d1495b"}
    labels = {"read_length_bp": "whole read", "utr5_reach_bp": "5′ UTR control", "utr3_reach_bp": "3′ UTR"}
    scales = {
        metric: max([1.0] + [abs(float(row["primary_delta_median_bp"])) for row in rows
                             if row["metric"] == metric and row["ranking_status"] == "eligible"])
        for metric in METRICS
    }
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
             '<rect width="100%" height="100%" fill="white"/>',
             f'<text x="24" y="31" font-family="sans-serif" font-size="20" font-weight="700">Primary regular-cDNA shifts for {html.escape(focus_label)}</text>',
             '<text x="24" y="55" font-family="sans-serif" font-size="13" fill="#444">Solid point: primary like-library pair. Small ticks: separately prepared directional counterchecks; never pooled.</text>']
    for index, metric in enumerate(METRICS):
        x0 = left + index * (panel + gap)
        center = x0 + panel / 2
        parts += [f'<text x="{center}" y="88" text-anchor="middle" font-family="sans-serif" font-size="14" fill="{colors[metric]}">{labels[metric]}</text>',
                  f'<line x1="{center}" x2="{center}" y1="101" y2="{height-36}" stroke="#777"/>',
                  f'<text x="{x0}" y="{height-12}" font-family="sans-serif" font-size="11" fill="#555">shorter ← Δ bp → longer</text>']
    for row_index, gene in enumerate(genes):
        y = 121 + row_index * 34
        parts.append(f'<text x="{left-12}" y="{y+5}" text-anchor="end" font-family="sans-serif" font-size="13">{html.escape(gene)}</text>')
        parts.append(f'<line x1="{left}" x2="{width-24}" y1="{y}" y2="{y}" stroke="#eee"/>')
        for index, metric in enumerate(METRICS):
            row = lookup.get((gene, metric))
            if not row or row["ranking_status"] != "eligible":
                continue
            center = left + index * (panel + gap) + panel / 2
            scale = panel / 2 - 12
            primary = float(row["primary_delta_median_bp"])
            x = center + primary / scales[metric] * scale
            parts.append(f'<circle cx="{x:.2f}" cy="{y}" r="5" fill="{colors[metric]}"><title>{html.escape(gene)} primary: {primary:g} bp</title></circle>')
            directional = json.loads(row["directional_deltas_bp_json"])
            for tick_index, (identifier, delta) in enumerate(sorted(directional.items())):
                tx = center + float(delta) / scales[metric] * scale
                dy = -7 if tick_index % 2 == 0 else 7
                parts.append(f'<line x1="{tx:.2f}" x2="{tx:.2f}" y1="{y+dy-3}" y2="{y+dy+3}" stroke="{colors[metric]}" stroke-width="2"><title>{html.escape(identifier)}: {float(delta):g} bp</title></line>')
    parts.append('</svg>\n')
    return "".join(parts)


def run(paths: list[Path], contrast_path: Path, focus_label: str, output: Path, minimum: int, max_censored: float, top: int) -> None:
    if not focus_label or minimum < 1 or not 0 <= max_censored <= 1 or top < 1:
        raise fail("invalid focus label, minimum, censoring threshold or top-gene count")
    rows, samples, contrasts = summarize(paths, contrast_path, focus_label, minimum, max_censored)
    output.mkdir(mode=0o700)
    with staged_output(output / "gene_contrasts.tsv") as handle:
        writer = csv.DictWriter(handle, fieldnames=SUMMARY_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    with staged_output(output / "top_genes.svg") as handle:
        handle.write(svg(rows, focus_label, top))
    with staged_output(output / "report.json") as handle:
        json.dump({
            "schema": SCHEMA,
            "focus_label": focus_label,
            "samples": samples,
            "contrasts": contrasts,
            "inputs": {
                "contrast_sheet": {"name": contrast_path.name, "sha256": hashlib.sha256(contrast_path.read_bytes()).hexdigest()},
                "projections": [
                    {"name": path.name, "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
                    for path in paths
                ],
            },
            "minimum_observed_per_sample": minimum,
            "maximum_right_censored_fraction": max_censored,
            "ranking_metric": "absolute primary regular-cDNA median difference",
            "directional_policy": "differently prepared contrasts are displayed separately and never pooled",
            "nonclaims": [
                "The ranking is descriptive and does not establish TP73, MDM2 or MDM4 causality.",
                "Enriched-library counterchecks do not replace the primary regular-cDNA contrast.",
                "Right-censored values are counts and lower-bound evidence, not exact lengths used in medians.",
            ],
        }, handle, indent=2, sort_keys=True)
        handle.write("\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--projection", action="append", required=True, type=Path)
    parser.add_argument("--contrast-sheet", required=True, type=Path)
    parser.add_argument("--focus-label", required=True)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--minimum-observed", type=int, default=10)
    parser.add_argument("--maximum-censored-fraction", type=float, default=0.25)
    parser.add_argument("--top-genes", type=int, default=20)
    args = parser.parse_args()
    try:
        run(args.projection, args.contrast_sheet, args.focus_label, args.output_dir,
            args.minimum_observed, args.maximum_censored_fraction, args.top_genes)
    except (OSError, ValueError) as exc:
        parser.error(str(exc))


if __name__ == "__main__":
    main()
