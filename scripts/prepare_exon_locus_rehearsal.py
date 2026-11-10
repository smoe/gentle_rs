#!/usr/bin/env python3
"""General synthetic overlapping-strand inputs, never real/private reads or a GENtle build."""

from __future__ import annotations

import json
from pathlib import Path

from scripts import prepare_read_end_rehearsal as base


def replace(path: Path, value, *, ending: str = "\n") -> None:
    # Only the newly created synthetic fixture files are rewritten before they are bound.
    path.write_bytes((json.dumps(value, sort_keys=True, indent=2) + "\n").replace("\n", ending).encode("utf-8"))


def prepare(root: Path, *, ending: str = "\n") -> None:
    if ending not in ("\n", "\r\n"):
        raise ValueError("synthetic fixture requires LF or CRLF")
    base.prepare(root, oriented=True)
    annotation = json.loads((root / "annotation.json").read_bytes())
    annotation["transcripts"].append({
        "transcript_id": "SYN_OPPOSITE", "gene_ids": ["SYN_OTHER_GENE"],
        "target_name": "synthetic_contig", "target_length": 1000,
        "genomic_contig": "synthetic_contig", "genomic_strand": "-",
        "genomic_exons": [[200, 250], [400, 450]],
    })
    replace(root / "annotation.json", annotation, ending=ending)
    models = json.loads((root / "models.json").read_bytes())
    models["annotation_sha256"] = base.bound(root / "annotation.json")["sha256"]
    replace(root / "models.json", models, ending=ending)
    originals = json.loads((root / "original-sequences.json").read_bytes())
    originals["synthetic_shared"] = "A" * 12 + "C" * 20 + "T" * 18
    replace(root / "original-sequences.json", originals, ending=ending)
    with (root / "synthetic.sam").open("a", encoding="ascii", newline="\n") as target:
        sequence = originals["synthetic_shared"].translate(str.maketrans("ACGT", "TGCA"))[::-1]
        target.write(f"synthetic_shared\t16\tsynthetic_contig\t201\t60\t50M\t*\t0\t0\t{sequence}\t*\n")
    if ending == "\r\n":
        for name in ("synthetic.sam", "dictionary.json"):
            path = root / name
            path.write_bytes(path.read_bytes().replace(b"\n", b"\r\n"))
    manifest = json.loads((root / "units.json").read_bytes())
    manifest["sam_sha256"] = base.bound(root / "synthetic.sam")["sha256"]
    manifest["units"].append({"unit_id": "synthetic_shared", "unclipped_axis_length": 50,
                              "rows": [{"ordinal": 3, "qname": "synthetic_shared", "read_group": None}],
                              "masks": None, "ends": [[], []], "search": manifest["units"][0]["search"]})
    replace(root / "units.json", manifest, ending=ending)
    request = json.loads((root / "sam-request.json").read_bytes())
    request["sam"] = base.bound(root / "synthetic.sam")
    request["manifest"] = base.bound(root / "units.json")
    request["dictionary"] = base.bound(root / "dictionary.json")
    request["annotation_sources"][0]["annotation"] = base.bound(root / "annotation.json")
    replace(root / "sam-request.json", request, ending=ending)


def bind(root: Path) -> dict:
    base.bind_analysis(root)
    base.bind_originals(root)
    previous = json.loads((root / "oriented-build-command.json").read_bytes())["request"]
    old = previous["samples"][0]
    command = {"action": "locus", "command": {"action": "build", "request": {
        "schema": "gentle.rna_read_end_locus_build.v1", "source_id": previous["source_id"],
        "region": {"contig": "synthetic_contig", "contig_length": 1000, "exon": [200, 250],
                   "flank_left": 150, "flank_right": 250, "anchor_transcript_id": "SYN_TX_A"},
        "samples": [{"sample_id": old["sample_id"], "analysis": old["analysis"],
                     "annotation": old["normalized_annotation"], "model": old["model"], "original_bases": old["original_bases"]}],
        "terminal_runs": previous["terminal_runs"], "limits": previous["limits"],
        "output_dir": str(root / "locus"), "max_output_bytes": 16 * 1024 * 1024,
        "cancel_path": str(root / "locus.cancel"),
    }}}
    base.write(root / "locus-build.json", command)
    return command


def view() -> dict:
    return {"schema": "gentle.rna_read_end_locus_view.v1", "projection": {"kind": "genome"},
            "sample_ids": [], "class": None, "offset": 0, "count": 100, "max_sublanes": 100, "width": 1400}
