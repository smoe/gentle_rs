#!/usr/bin/env python3
"""Prepare a frozen CDS-aware transcript model from an explicit GTF file."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import re
from pathlib import Path

from compare_rna_read_end_metrics import staged_output
from project_rna_read_end_metrics import HASH, MODEL_SCHEMA, TOKEN, fail, merged


ATTRIBUTE = re.compile(r'(?:^|;\s*)([A-Za-z][A-Za-z0-9_]*)\s+"([^"]*)"')


def raw_and_text(path: Path) -> tuple[bytes, str]:
    raw = path.read_bytes()
    try:
        payload = gzip.decompress(raw) if path.name.endswith(".gz") else raw
        return raw, payload.decode("utf-8")
    except (OSError, UnicodeDecodeError) as exc:
        raise fail(f"GTF is not valid UTF-8 text: {path}") from exc


def parse_gtf(text: str):
    records: dict[str, dict] = {}
    for line_number, line in enumerate(text.splitlines(), 1):
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) != 9 or fields[2] not in {"exon", "CDS"}:
            continue
        contig, feature, strand = fields[0], fields[2], fields[6]
        if not TOKEN.fullmatch(contig) or strand not in {"+", "-"}:
            raise fail(f"invalid contig or strand at GTF line {line_number}")
        try:
            start, end = int(fields[3]) - 1, int(fields[4])
        except ValueError as exc:
            raise fail(f"invalid coordinates at GTF line {line_number}") from exc
        if start < 0 or start >= end:
            raise fail(f"invalid interval at GTF line {line_number}")
        attributes = dict(ATTRIBUTE.findall(fields[8]))
        transcript_id, gene_id = attributes.get("transcript_id"), attributes.get("gene_id")
        if not transcript_id or not gene_id or not TOKEN.fullmatch(transcript_id) or not TOKEN.fullmatch(gene_id):
            raise fail(f"missing or invalid gene_id/transcript_id at GTF line {line_number}")
        row = records.setdefault(transcript_id, {
            "transcript_id": transcript_id, "gene_id": gene_id,
            "contig": contig, "strand": strand, "exons": [], "cds": [],
        })
        if (row["gene_id"], row["contig"], row["strand"]) != (gene_id, contig, strand):
            raise fail(f"inconsistent transcript identity at GTF line {line_number}")
        row["exons" if feature == "exon" else "cds"].append((start, end))
    transcripts = []
    for transcript_id in sorted(records):
        row = records[transcript_id]
        exons, cds = merged(set(row["exons"])), merged(set(row["cds"]))
        if not exons or not cds:
            continue
        if any(sum(max(0, min(z, ez) - max(a, ea)) for ea, ez in exons) != z - a for a, z in cds):
            raise fail(f"CDS falls outside exons for {transcript_id}")
        transcripts.append({
            "transcript_id": transcript_id,
            "gene_id": row["gene_id"],
            "contig": row["contig"],
            "strand": row["strand"],
            "genomic_exons_0based_half_open": [list(value) for value in exons],
            "genomic_cds_0based_half_open": [list(value) for value in cds],
        })
    if not transcripts:
        raise fail("GTF contains no transcript with both exon and CDS records")
    return transcripts


def run(gtf: Path, reference: Path, assembly: str, source_id: str, annotation_sha256: str, output: Path) -> None:
    if not assembly or not TOKEN.fullmatch(source_id) or not HASH.fullmatch(annotation_sha256):
        raise fail("assembly, source_id or annotation SHA-256 is invalid")
    gtf_raw, text = raw_and_text(gtf)
    reference_raw = reference.read_bytes()
    document = {
        "schema": MODEL_SCHEMA,
        "source_id": source_id,
        "assembly": assembly,
        "reference_fasta_sha256": "sha256:" + hashlib.sha256(reference_raw).hexdigest(),
        "annotation_sha256": annotation_sha256,
        "coordinate_convention": "zero_based_half_open",
        "model_source": {
            "format": "gtf",
            "name": gtf.name,
            "sha256": "sha256:" + hashlib.sha256(gtf_raw).hexdigest(),
        },
        "transcripts": parse_gtf(text),
    }
    with staged_output(output) as handle:
        json.dump(document, handle, indent=2, sort_keys=True)
        handle.write("\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gtf", required=True, type=Path)
    parser.add_argument("--reference-fasta", required=True, type=Path)
    parser.add_argument("--assembly", required=True)
    parser.add_argument("--source-id", required=True)
    parser.add_argument("--annotation-sha256", required=True)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    try:
        run(args.gtf, args.reference_fasta, args.assembly, args.source_id, args.annotation_sha256, args.output)
    except (OSError, ValueError) as exc:
        parser.error(str(exc))


if __name__ == "__main__":
    main()
