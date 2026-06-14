#!/usr/bin/env python3
"""Derive GENtle's tiny Clariom D hg38 gene-panel fixture from vendor ZIPs."""

from __future__ import annotations

import argparse
import csv
import io
import zipfile
from pathlib import Path


GENE_PANEL = [
    "E2F1",
    "TP73",
    "SP1",
    "PATZ1",
    "TP53",
    "TP63",
    "IL6",
    "IL10",
    "FUS",
    "TERT",
    "TARDBP",
    "MDM2",
    "CDKN1A",
    "BAX",
    "GADD45A",
    "MYC",
    "RB1",
    "ESR1",
    "GAPDH",
    "ACTB",
    "SRSF1",
]

PROBESET_MEMBER = "Clariom_D_Human.na36.hg38.probeset.csv"
TRANSCRIPT_MEMBER = "Clariom_D_Human.r1.na36.hg38.a1.transcript.csv"

PROBESET_COLUMNS = [
    "gene_symbols",
    "probeset_id",
    "seqname",
    "strand",
    "start",
    "stop",
    "probe_count",
    "transcript_cluster_id",
    "locus_type",
    "exon_id",
    "psr_id",
    "probeset_type",
    "psr_type",
    "junction_start_edge",
    "junction_stop_edge",
    "level",
    "has_cds",
]

TRANSCRIPT_COLUMNS = [
    "gene_symbols",
    "transcript_cluster_id",
    "probeset_id",
    "seqname",
    "strand",
    "start",
    "stop",
    "total_probes",
    "category",
    "locus type",
    "notes",
]


def gene_symbols(value: str | None) -> set[str]:
    symbols: set[str] = set()
    for item in (value or "").split(" /// "):
        parts = [part.strip() for part in item.split(" // ")]
        if len(parts) > 1 and parts[1] and parts[1] != "---":
            symbols.add(parts[1])
    return symbols


def iter_vendor_rows(zip_path: Path, member: str):
    with zipfile.ZipFile(zip_path) as archive:
        with archive.open(member) as handle:
            text = io.TextIOWrapper(
                handle, encoding="utf-8-sig", errors="replace", newline=""
            )
            rows = (line for line in text if not line.startswith("#"))
            yield from csv.DictReader(rows)


def write_subset(zip_path: Path, member: str, output: Path, columns: list[str]) -> int:
    wanted = set(GENE_PANEL)
    count = 0
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=columns, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for row in iter_vendor_rows(zip_path, member):
            symbols = sorted(gene_symbols(row.get("gene_assignment")) & wanted)
            if not symbols:
                continue
            out = {column: row.get(column, "") for column in columns}
            out["gene_symbols"] = ";".join(symbols)
            writer.writerow(out)
            count += 1
    return count


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--probeset-zip", required=True, type=Path)
    parser.add_argument("--transcript-zip", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()

    probeset_count = write_subset(
        args.probeset_zip,
        PROBESET_MEMBER,
        args.output_dir / "clariom_d_human_na36_hg38_gene_panel.probesets.tsv",
        PROBESET_COLUMNS,
    )
    transcript_count = write_subset(
        args.transcript_zip,
        TRANSCRIPT_MEMBER,
        args.output_dir / "clariom_d_human_na36_hg38_gene_panel.transcripts.tsv",
        TRANSCRIPT_COLUMNS,
    )
    print(f"wrote {probeset_count} probeset rows")
    print(f"wrote {transcript_count} transcript rows")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
