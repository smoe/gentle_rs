#!/usr/bin/env python3
"""Deterministically convert every chromosome in a BigWig to bedGraph.

This is a small, auditable replacement for the UCSC ``bigWigToBedGraph``
interface used by GENtle's local-track importer.  Unlike the historical
analysis helper, it does not select or assume a chromosome.
"""

from __future__ import annotations

import math
from pathlib import Path
import sys

import pyBigWig


VERSION = "gentle-bigwig-to-bedgraph 1"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise SystemExit(message)


def convert(source: Path, destination: Path) -> int:
    require(source.is_file(), f"BigWig input does not exist: {source}")
    require(not destination.exists()
            or (destination.is_file() and destination.stat().st_size == 0),
            f"Refusing to overwrite nonempty or non-file output: {destination}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    count = 0
    try:
        bigwig = pyBigWig.open(str(source))
    except RuntimeError as error:
        raise SystemExit(f"Cannot open BigWig input: {error}") from error
    require(bigwig is not None and bigwig.isBigWig(), "Input is not a BigWig file")
    try:
        chromosomes = bigwig.chroms()
        require(bool(chromosomes), "BigWig contains no chromosomes")
        with destination.open("w", encoding="utf-8", newline="\n") as output:
            for chromosome, length in chromosomes.items():
                require(isinstance(chromosome, str) and chromosome and int(length) > 0,
                        "BigWig contains an invalid chromosome inventory")
                for start, end, score in bigwig.intervals(chromosome) or ():
                    require(0 <= start < end <= length,
                            f"Invalid interval in {chromosome}: {start}-{end}")
                    require(math.isfinite(score),
                            f"Non-finite score in {chromosome}: {start}-{end}")
                    output.write(f"{chromosome}\t{start}\t{end}\t{score:.17g}\n")
                    count += 1
    finally:
        bigwig.close()
    return count


def main() -> None:
    if len(sys.argv) == 2 and sys.argv[1] in {"--version", "-version"}:
        print(VERSION)
        return
    require(len(sys.argv) == 3,
            "usage: bigwig_to_bedgraph.py INPUT.bigWig OUTPUT.bedGraph")
    convert(Path(sys.argv[1]).resolve(), Path(sys.argv[2]).resolve())


if __name__ == "__main__":
    main()
