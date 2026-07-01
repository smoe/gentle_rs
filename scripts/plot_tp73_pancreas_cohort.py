#!/usr/bin/env python3
"""Compatibility wrapper for the generic single-gene pancreas plotter.

The canonical renderer is ``plot_pancreas_gene_screen.py``. This script keeps the
old TP73-specific command line and default ``$COHORT_ROOT`` input path for saved
runbooks.
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path


def parse_args() -> argparse.Namespace:
    cohort_root = os.environ.get("COHORT_ROOT")
    default_input = (
        Path(cohort_root) / "figures" / "tp73_pancreas_figure_source.tsv"
        if cohort_root
        else None
    )
    parser = argparse.ArgumentParser(
        description=(
            "Render TP73 pancreas cohort read-length and TP73-support SVG plots. "
            "This wrapper forwards to plot_pancreas_gene_screen.py."
        ),
        epilog=(
            "Bar-column semantics:\n"
            "  strict_seed_passed_reads\n"
            "    Reads that pass the first seed phase before any sequence-alignment\n"
            "    confirmation. This is the default view.\n\n"
            "  accepted_tp73_reads\n"
            "    Downstream accepted TP73-support count from older cohort TSVs,\n"
            "    after target-support confirmation."
        ),
        formatter_class=argparse.RawTextHelpFormatter,
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
        "--bar-column",
        default="strict_seed_passed_reads",
        choices=["strict_seed_passed_reads", "accepted_tp73_reads"],
        help=(
            "TP73 count column used for blue bars. See the bar-column semantics below."
        ),
    )
    parser.add_argument(
        "--title",
        default="TP73 support across pancreatic cancer Nanopore cDNA samples",
        help="Figure title.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if not args.input_tsv:
        print("Input TSV required, or set COHORT_ROOT for the legacy default path.", file=sys.stderr)
        return 2
    script = Path(__file__).with_name("plot_pancreas_gene_screen.py")
    command = [
        sys.executable,
        str(script),
        args.input_tsv,
        "--gene",
        "TP73",
        "--bar-column",
        args.bar_column,
        "--label-column",
        args.label_column,
        "--title",
        args.title,
    ]
    if args.output:
        command.extend(["-o", args.output])
    return subprocess.call(command)


if __name__ == "__main__":
    raise SystemExit(main())
