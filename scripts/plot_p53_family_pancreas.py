#!/usr/bin/env python3
"""Compatibility wrapper for the generic pancreas gene-family plotter.

Historically this script rendered the TP53/TP63/TP73 figure directly. The
canonical implementation now lives in ``plot_pancreas_gene_family.py`` so the
same plotting code can be reused for E2F1-like follow-up families without
copying the SVG renderer.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Render the p53-family pancreas grouped-bar SVG. This wrapper keeps "
            "the old TP53/TP63/TP73 command line and forwards to "
            "plot_pancreas_gene_family.py."
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
    parser.add_argument("-o", "--output", required=True, help="Output SVG path.")
    parser.add_argument(
        "--output-tsv",
        default=None,
        help="Canonical family TSV path. Defaults to OUTPUT.with_suffix('.tsv').",
    )
    parser.add_argument(
        "--metric",
        choices=["per_million", "raw"],
        default="per_million",
        help="Grouped-bar height metric.",
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


def main() -> int:
    args = parse_args()
    script = Path(__file__).with_name("plot_pancreas_gene_family.py")
    command = [sys.executable, str(script)]
    if args.tp53_batch_report:
        command.extend(["--batch-report", f"TP53={args.tp53_batch_report}"])
    else:
        command.extend(["--batch-summary", f"TP53={args.tp53_batch_summary}"])
    command.extend(
        [
            "--figure-source",
            f"TP63={args.tp63_figure_source}",
            "--figure-source",
            f"TP73={args.tp73_figure_source}",
            "--genes",
            "TP53,TP63,TP73",
            "--metric",
            args.metric,
            "--label-column",
            args.label_column,
            "--title",
            args.title,
            "-o",
            args.output,
        ]
    )
    if args.output_tsv:
        command.extend(["--output-tsv", args.output_tsv])
    return subprocess.call(command)


if __name__ == "__main__":
    raise SystemExit(main())
