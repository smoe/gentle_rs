#!/usr/bin/env python3
"""Rank genes by sample-specific RNA read-end shifts.

The input is one or more ``projections.tsv`` files produced by
``project_rna_read_end_metrics.py``.  Rankings are descriptive: each score is
the median of focus-sample minus comparison-sample median differences.  This
keeps samples, rather than individual reads, as the comparison replicates.
"""

from __future__ import annotations

import argparse
import csv
import html
import json
import statistics
from collections import defaultdict
from pathlib import Path

from compare_rna_read_end_metrics import staged_output


SCHEMA = "gentle.rna_read_end_gene_ranking.v1"
METRICS = ("read_length_bp", "utr5_reach_bp", "utr3_reach_bp")
HEADER = (
    "sample_id", "source_sha256", "unit_id", "metric", "value_bp",
    "observation_status", "reason", "transcript_ids", "gene_ids",
)
SUMMARY_FIELDS = (
    "gene_id", "metric", "focus_sample", "focus_observed_n",
    "focus_right_censored_n", "focus_missing_n", "focus_median_bp",
    "comparison_samples_eligible", "comparison_sample_count",
    "median_sample_delta_bp", "min_sample_delta_bp", "max_sample_delta_bp",
    "ranking_status",
)


def fail(message: str) -> ValueError:
    return ValueError(f"rna_read_end_gene_ranking: {message}")


def load(paths: list[Path]):
    values: dict[tuple[str, str, str], list[int]] = defaultdict(list)
    counts: dict[tuple[str, str, str], dict[str, int]] = defaultdict(lambda: defaultdict(int))
    samples: set[str] = set()
    seen: set[tuple[str, str, str]] = set()
    for path in paths:
        with path.open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if tuple(reader.fieldnames or ()) != HEADER:
                raise fail(f"unexpected columns in {path}")
            for row_number, row in enumerate(reader, 2):
                sample, unit, metric = row["sample_id"], row["unit_id"], row["metric"]
                identity = (sample, unit, metric)
                if identity in seen:
                    raise fail(f"duplicate sample/unit/metric at {path}:{row_number}")
                seen.add(identity)
                samples.add(sample)
                if metric not in METRICS:
                    raise fail(f"unsupported metric at {path}:{row_number}")
                genes = row["gene_ids"].split(",") if row["gene_ids"] else []
                if len(genes) != 1:
                    continue
                status = row["observation_status"]
                if status not in {"observed", "right_censored", "missing"}:
                    raise fail(f"unsupported status at {path}:{row_number}")
                key = (genes[0], sample, metric)
                counts[key][status] += 1
                if status == "observed":
                    try:
                        value = int(row["value_bp"])
                    except ValueError as exc:
                        raise fail(f"non-integer observed value at {path}:{row_number}") from exc
                    if value < 0:
                        raise fail(f"negative observed value at {path}:{row_number}")
                    values[key].append(value)
                elif row["value_bp"]:
                    # Censored values are lower bounds, not exact observations.
                    try:
                        if int(row["value_bp"]) < 0:
                            raise ValueError
                    except ValueError as exc:
                        raise fail(f"invalid censored bound at {path}:{row_number}") from exc
    return values, counts, sorted(samples)


def summarize(paths: list[Path], focus: str, minimum: int, max_censored: float):
    values, counts, samples = load(paths)
    if focus not in samples:
        raise fail(f"focus sample {focus!r} is absent")
    comparisons = [sample for sample in samples if sample != focus]
    if not comparisons:
        raise fail("at least two samples are required")
    required_comparisons = (len(comparisons) + 1) // 2
    genes = sorted({key[0] for key in counts})
    rows = []
    for gene in genes:
        for metric in METRICS:
            focus_key = (gene, focus, metric)
            focus_values = values[focus_key]
            focus_counts = counts[focus_key]
            focus_total = focus_counts["observed"] + focus_counts["right_censored"]
            focus_fraction = focus_counts["right_censored"] / focus_total if focus_total else 0.0
            deltas = []
            eligible_samples = []
            for sample in comparisons:
                key = (gene, sample, metric)
                sample_values = values[key]
                sample_counts = counts[key]
                total = sample_counts["observed"] + sample_counts["right_censored"]
                censored_fraction = sample_counts["right_censored"] / total if total else 0.0
                if (
                    len(focus_values) >= minimum
                    and len(sample_values) >= minimum
                    and focus_fraction <= max_censored
                    and censored_fraction <= max_censored
                ):
                    deltas.append(statistics.median(focus_values) - statistics.median(sample_values))
                    eligible_samples.append(sample)
            rows.append({
                "gene_id": gene,
                "metric": metric,
                "focus_sample": focus,
                "focus_observed_n": len(focus_values),
                "focus_right_censored_n": focus_counts["right_censored"],
                "focus_missing_n": focus_counts["missing"],
                "focus_median_bp": "" if not focus_values else statistics.median(focus_values),
                "comparison_samples_eligible": ",".join(eligible_samples),
                "comparison_sample_count": len(eligible_samples),
                "median_sample_delta_bp": "" if not deltas else statistics.median(deltas),
                "min_sample_delta_bp": "" if not deltas else min(deltas),
                "max_sample_delta_bp": "" if not deltas else max(deltas),
                "ranking_status": "eligible" if len(deltas) >= required_comparisons else "insufficient_or_censored",
            })
    return rows, samples


def svg(rows: list[dict], focus: str, top: int) -> str:
    three = [row for row in rows if row["metric"] == "utr3_reach_bp" and row["ranking_status"] == "eligible"]
    selected = sorted(three, key=lambda row: (-abs(float(row["median_sample_delta_bp"])), row["gene_id"]))[:top]
    genes = [row["gene_id"] for row in selected]
    lookup = {(row["gene_id"], row["metric"]): row for row in rows}
    width, left, panel, gap = 1180, 190, 290, 55
    height = max(210, 115 + 34 * len(genes))
    scales = {
        metric: max(
            [1.0] + [abs(float(row["median_sample_delta_bp"])) for row in rows
                     if row["metric"] == metric and row["ranking_status"] == "eligible"]
        )
        for metric in METRICS
    }
    colors = {"read_length_bp": "#6c7a89", "utr5_reach_bp": "#2a9d8f", "utr3_reach_bp": "#d1495b"}
    labels = {"read_length_bp": "whole read", "utr5_reach_bp": "5′ UTR control", "utr3_reach_bp": "3′ UTR"}
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
             '<rect width="100%" height="100%" fill="white"/>',
             f'<text x="24" y="31" font-family="sans-serif" font-size="20" font-weight="700">Genes with the largest {html.escape(focus)} 3′-UTR shifts</text>',
             '<text x="24" y="55" font-family="sans-serif" font-size="13" fill="#444">Median of sample-level median differences; censored reads are counted but not treated as exact lengths.</text>']
    for index, metric in enumerate(METRICS):
        x0 = left + index * (panel + gap)
        center = x0 + panel / 2
        parts += [f'<text x="{center}" y="82" text-anchor="middle" font-family="sans-serif" font-size="14" fill="{colors[metric]}">{labels[metric]}</text>',
                  f'<line x1="{center}" x2="{center}" y1="94" y2="{height-32}" stroke="#777"/>',
                  f'<text x="{x0}" y="{height-10}" font-family="sans-serif" font-size="11" fill="#555">shorter ← Δ bp → longer</text>']
    for row_index, gene in enumerate(genes):
        y = 108 + row_index * 34
        parts.append(f'<text x="{left-12}" y="{y+5}" text-anchor="end" font-family="sans-serif" font-size="13">{html.escape(gene)}</text>')
        parts.append(f'<line x1="{left}" x2="{width-24}" y1="{y}" y2="{y}" stroke="#eee"/>')
        for index, metric in enumerate(METRICS):
            row = lookup.get((gene, metric))
            if not row or row["ranking_status"] != "eligible":
                continue
            delta = float(row["median_sample_delta_bp"])
            center = left + index * (panel + gap) + panel / 2
            x = center + (delta / scales[metric]) * (panel / 2 - 12)
            parts.append(f'<circle cx="{x:.2f}" cy="{y}" r="5" fill="{colors[metric]}"><title>{html.escape(gene)}: {delta:g} bp</title></circle>')
    parts.append('</svg>\n')
    return "".join(parts)


def run(paths: list[Path], focus: str, output: Path, minimum: int, max_censored: float, top: int) -> None:
    if minimum < 1 or not 0 <= max_censored <= 1 or top < 1:
        raise fail("invalid minimum, censoring threshold or top-gene count")
    rows, samples = summarize(paths, focus, minimum, max_censored)
    output.mkdir(mode=0o700)
    with staged_output(output / "gene_summary.tsv") as handle:
        writer = csv.DictWriter(handle, fieldnames=SUMMARY_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    with staged_output(output / "top_genes.svg") as handle:
        handle.write(svg(rows, focus, top))
    with staged_output(output / "report.json") as handle:
        json.dump({
            "schema": SCHEMA,
            "focus_sample": focus,
            "samples": samples,
            "minimum_observed_per_sample": minimum,
            "maximum_right_censored_fraction": max_censored,
            "minimum_eligible_comparison_samples": (len(samples) - 1 + 1) // 2,
            "method": "median across per-comparison-sample differences of read-level medians",
            "nonclaims": [
                "The ranking is descriptive and does not establish TP73, MDM2 or MDM4 causality.",
                "Whole-read and 5-prime-UTR shifts are shown as truncation/library controls, not corrections.",
                "Right-censored values are counts and lower-bound evidence, not exact lengths used in medians.",
            ],
        }, handle, indent=2, sort_keys=True)
        handle.write("\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--projection", action="append", required=True, type=Path)
    parser.add_argument("--focus-sample", required=True)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--minimum-observed", type=int, default=10)
    parser.add_argument("--maximum-censored-fraction", type=float, default=0.25)
    parser.add_argument("--top-genes", type=int, default=20)
    args = parser.parse_args()
    try:
        run(args.projection, args.focus_sample, args.output_dir, args.minimum_observed, args.maximum_censored_fraction, args.top_genes)
    except (OSError, ValueError) as exc:
        parser.error(str(exc))


if __name__ == "__main__":
    main()
