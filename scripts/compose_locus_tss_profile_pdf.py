#!/usr/bin/env python3
"""Compose one bound locus-context page with its selected TSS TFBS pages."""

from __future__ import annotations

import argparse
import html
import io
import zipfile
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys
from typing import Any
import xml.etree.ElementTree as ET


LOCUS_RECEIPT_SCHEMA = "gentle.tss_local_similarity_locus_report_receipt.v1"
LOCUS_SVG_SCHEMA = "gentle.gene_locus_evidence_display.v1"
TSS_REPORT_SCHEMA = "gentle.tss_tfbs_profiles.v1"
TSS_INDEX_SCHEMA = "gentle.tss_tfbs_profile_index.v1"
TSS_RECEIPT_SCHEMA = "gentle.tss_tfbs_profile_receipt.v1"
OUTPUT_SCHEMA = "gentle.integrated_locus_tss_tfbs_pdf_receipt.v2"
TARGET_FASTA_SCHEMA = "gentle.target_tss_fasta_export.v1"
LOCUS_PAGE_WIDTH = 1400.0
LOCUS_PLOT_LEFT = 255.0
LOCUS_PLOT_RIGHT = 1050.0
FONT_IDENTITY_STATUS = "glyph_used_font_sources_recorded"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_bytes())
    require(isinstance(value, dict), f"{path} is not a JSON object")
    return value


def normalized_digest(value: Any) -> str:
    require(isinstance(value, str), "expected a SHA-256 string")
    digest = value.removeprefix("sha256:")
    require(len(digest) == 64 and all(c in "0123456789abcdef" for c in digest),
            "expected a lowercase SHA-256 digest")
    return digest


def bound_output(receipt: dict[str, Any], relative_name: str, path: Path) -> None:
    outputs = receipt.get("outputs")
    require(isinstance(outputs, dict), "receipt outputs must be an object")
    require(relative_name in outputs, f"receipt does not bind output {relative_name}")
    require(normalized_digest(outputs[relative_name]) == sha256(path),
            f"output hash mismatch for {relative_name}")


def bound_input(receipt: dict[str, Any], role: str, name: str, path: Path) -> None:
    inputs = receipt.get("inputs")
    require(isinstance(inputs, list), "receipt inputs must be an array")
    matches = [item for item in inputs if item.get("role") == role and item.get("name") == name]
    require(len(matches) == 1, f"receipt must bind exactly one {role} input named {name}")
    require(normalized_digest(matches[0].get("sha256")) == sha256(path),
            f"input hash mismatch for {role} {name}")


def svg_frame(svg_path: Path) -> tuple[float, float, float]:
    root = ET.fromstring(svg_path.read_bytes())
    width = float(root.get("width", "0"))
    left = float(root.get("data-gentle-plot-left", "nan"))
    right = float(root.get("data-gentle-plot-right", "nan"))
    return width, left, right


def locus_bands(svg_path: Path) -> tuple[dict[str, str], list[tuple[str, int, int]]]:
    root = ET.fromstring(svg_path.read_bytes())
    require(root.get("data-gentle-schema") == LOCUS_SVG_SCHEMA,
            "locus SVG has the wrong scientific schema")
    bands = []
    for item in root.iter():
        name = item.get("data-gentle-tss-stretch-band")
        if name is None:
            continue
        start = int(item.get("data-gentle-genomic-start", "0"))
        end = int(item.get("data-gentle-genomic-end", "0"))
        require(name and 0 < start <= end, "locus SVG contains an invalid TSS band")
        bands.append((name, start, end))
    require(bands, "locus SVG contains no TSS background bands")
    require(float(root.get("width", "0")) == LOCUS_PAGE_WIDTH,
            "locus SVG does not use the canonical 1400-pixel page frame")
    return root.attrib, bands


def parse_fasta(path: Path) -> dict[str, dict[str, str]]:
    records: dict[str, dict[str, str]] = {}
    header: str | None = None
    sequence_lines: list[str] = []

    def finish() -> None:
        nonlocal header, sequence_lines
        if header is None:
            return
        fields: dict[str, str] = {}
        pieces = header.split("|")
        fields["gene_symbol"] = pieces[0]
        for piece in pieces[1:]:
            require("=" in piece, "FASTA header contains an unkeyed field")
            key, value = piece.split("=", 1)
            require(key and key not in fields, "FASTA header contains duplicate fields")
            require(key not in {"header", "sequence"}, "FASTA header uses a reserved field")
            fields[key] = value
        promoter_id = fields.get("promoter_id")
        require(promoter_id and promoter_id not in records,
                "FASTA promoter IDs must be present and unique")
        sequence = "".join(sequence_lines).upper()
        require(sequence and all(base in "ACGTN" for base in sequence),
                f"FASTA {promoter_id} contains unsupported sequence symbols")
        sequence_digest = hashlib.sha256(sequence.encode("ascii")).hexdigest()
        require(sequence_digest == normalized_digest(fields.get("sequence_sha256")),
                f"FASTA {promoter_id} sequence digest does not match its bases")
        records[promoter_id] = {
            **fields,
            "header": header,
            "sequence": sequence,
            "sequence_sha256": sequence_digest,
        }
        header = None
        sequence_lines = []

    for raw in path.read_text(encoding="ascii").splitlines():
        if raw.startswith(">"):
            finish()
            header = raw[1:]
        else:
            require(header is not None, "FASTA sequence precedes its header")
            sequence_lines.append("".join(raw.split()))
    finish()
    require(records, "FASTA contains no records")
    return records


def selected_fasta(
    gene: str,
    manifest_path: Path,
    tss_receipt: dict[str, Any],
    bindings: list[dict[str, Any]],
    reference: dict[str, Any],
) -> tuple[str, dict[str, Any]]:
    manifest = read_json(manifest_path)
    require(manifest.get("schema") == TARGET_FASTA_SCHEMA,
            "unexpected target TSS FASTA manifest schema")
    require(manifest.get("upstream_bp") == 500 and manifest.get("downstream_bp") == 200
            and manifest.get("sequence_orientation") == "transcript_5prime_to_3prime",
            "selected FASTA requires transcript-oriented -500/+200 source windows")
    assembly = reference.get("assembly")
    require(assembly and manifest.get("assembly_id") == assembly,
            "FASTA manifest and TSS report assembly disagree")
    bound_input(tss_receipt, "bundle_manifest", manifest_path.name, manifest_path)
    checksums = manifest_path.parent / "SHA256SUMS"
    require(checksums.is_file(), "target TSS FASTA checksum inventory is missing")
    bound_input(tss_receipt, "bundle_checksums", checksums.name, checksums)
    require(normalized_digest(manifest.get("sha256sums_sha256")) == sha256(checksums),
            "manifest/checksum-inventory binding mismatch")
    checksum_rows = {}
    for line in checksums.read_text(encoding="ascii").splitlines():
        digest, separator, name = line.partition("  ")
        require(separator and name and name not in checksum_rows,
                "checksum inventory contains a malformed or duplicate row")
        checksum_rows[name] = normalized_digest(digest)

    files = [item for item in manifest.get("files", []) if item.get("gene_symbol") == gene]
    require(len(files) == 1, f"FASTA manifest must contain exactly one {gene} file")
    file_row = files[0]
    require(Path(file_row["filename"]).name == file_row["filename"],
            "target FASTA must be a direct file in the bound bundle")
    fasta_path = manifest_path.parent / file_row["filename"]
    require(fasta_path.is_file() and fasta_path.parent == manifest_path.parent,
            "target FASTA must be a direct file in the bound bundle")
    file_digest = sha256(fasta_path)
    require(normalized_digest(file_row.get("sha256")) == file_digest,
            "FASTA hash does not match its manifest entry")
    require(checksum_rows.get(fasta_path.name) == file_digest,
            "FASTA hash does not match the checksum inventory")
    bound_input(tss_receipt, "fasta", fasta_path.name, fasta_path)
    records = parse_fasta(fasta_path)
    rows = file_row.get("records", [])
    manifest_records = {item["promoter_id"]: item for item in rows}
    require(len(manifest_records) == len(rows) and set(manifest_records) == set(records),
            "FASTA/manifest promoter IDs must be unique and equal")
    require(len(manifest_records) == file_row.get("record_count") == len(records),
            "FASTA manifest record count or promoter IDs disagree")

    output = []
    for binding in bindings:
        promoter_id = binding["promoter_id"]
        require(promoter_id in records and promoter_id in manifest_records,
                f"selected promoter {promoter_id} is absent from the bound FASTA")
        record = records[promoter_id]
        declared = manifest_records[promoter_id]
        require(record.get("gene_symbol") == gene and record.get("assembly") == assembly,
                f"selected promoter {promoter_id} gene/assembly mismatch")
        require(record.get("window") == "minus500_plus200"
                and record.get("orientation") == "transcript_5prime_to_3prime",
                f"selected promoter {promoter_id} window/orientation mismatch")
        require(record["sequence_sha256"] == binding["sequence_sha256"]
                == normalized_digest(declared.get("sequence_sha256")),
                f"selected promoter {promoter_id} sequence digest mismatch")
        require(record.get("gene_id") == binding["gene_id"] == declared.get("gene_id")
                == file_row.get("gene_id"),
                f"selected promoter {promoter_id} gene binding mismatch")
        require(record.get("chromosome") == binding["chromosome"]
                == declared.get("chromosome"),
                f"selected promoter {promoter_id} chromosome mismatch")
        require(record.get("strand") == binding["strand"] == declared.get("strand"),
                f"selected promoter {promoter_id} strand mismatch")
        require(int(record.get("tss_1based", "0")) == binding["tss_1based"]
                == declared.get("tss_1based"),
                f"selected promoter {promoter_id} TSS mismatch")
        require(len(record["sequence"]) == declared.get("sequence_length_bp") == 701,
                f"selected promoter {promoter_id} is not a 701-bp window")
        tss = binding["tss_1based"]
        expected_span = (tss - 500, tss + 200) if binding["strand"] == "+" else (tss - 200, tss + 500)
        span = re.fullmatch(r"([1-9][0-9]*)-([1-9][0-9]*)", record.get("genomic_1based", ""))
        require(span is not None and expected_span[0] > 0
                and tuple(map(int, span.groups())) == expected_span
                == (binding["start_1based"], binding["end_1based"])
                == (declared.get("genomic_start_1based"), declared.get("genomic_end_1based")),
                f"selected promoter {promoter_id} genomic span mismatch")
        transcript_lists = [record.get("transcripts", "").split(","),
                            binding["transcripts"], declared.get("transcript_ids", [])]
        require(all(isinstance(ids, list) and ids
                    and all(isinstance(item, str) and item.strip() for item in ids)
                    and len(set(ids)) == len(ids) for ids in transcript_lists)
                and set(transcript_lists[0]) == set(transcript_lists[1]) == set(transcript_lists[2]),
                f"selected promoter {promoter_id} transcript membership mismatch")
        output.append(">" + record["header"] + "\n")
        sequence = record["sequence"]
        output.extend(sequence[index:index + 80] + "\n" for index in range(0, len(sequence), 80))
    return "".join(output), {
        "source_manifest": {"path": str(manifest_path), "sha256": sha256(manifest_path)},
        "source_checksums": {"path": str(checksums), "sha256": sha256(checksums)},
        "source_fasta": {"path": str(fasta_path), "sha256": file_digest},
        "record_count": len(bindings),
        "promoter_ids": [binding["promoter_id"] for binding in bindings],
    }


def select_pages(
    gene: str,
    report: dict[str, Any],
    index: dict[str, Any],
    tss_dir: Path,
    tss_receipt: dict[str, Any],
    bands: list[tuple[str, int, int]],
) -> tuple[list[Path], list[dict[str, Any]]]:
    selected = [window for window in report.get("windows", [])
                if window.get("selected") is True
                and window.get("record", {}).get("gene_symbol") == gene]
    require(selected, f"TSS report has no selected windows for {gene}")
    order = {row["record"]["promoter_id"]: i for i, row in enumerate(selected)}
    require(len(order) == len(selected), "selected promoter IDs are not unique")
    bindings = []
    for window in selected:
        record = window["record"]
        geometry = record.get("geometry", {})
        tss = geometry.get("tss_1based")
        require(isinstance(tss, int) and tss > 0, "selected TSS coordinate is invalid")
        covering = [name for name, start, end in bands if start <= tss <= end]
        require(covering, f"selected TSS {record['promoter_id']} is outside every locus band")
        require(geometry.get("strand") in {"+", "-"}, "selected TSS strand is invalid")
        require(isinstance(geometry.get("chromosome"), str), "selected chromosome is missing")
        evidence = window.get("selection_evidence")
        require(isinstance(evidence, dict) and evidence.get("label"),
                "selected TSS lacks bound selection evidence")
        bindings.append({
            "promoter_id": record["promoter_id"],
            "gene_id": record["gene_id"],
            "chromosome": geometry["chromosome"],
            "strand": geometry["strand"],
            "tss_1based": tss,
            "start_1based": geometry["start_1based"],
            "end_1based": geometry["end_1based"],
            "sequence_sha256": record["sequence_sha256"],
            "transcripts": record["transcripts"],
            "covering_tss_bands": covering,
            "selection_label": evidence["label"],
            "selection_factor": evidence.get("factor"),
            "selection_criterion": evidence.get("criterion"),
            "historical_selection_window": evidence.get("selection_window"),
            "selection_window_join_policy": "provenance_only_not_display_geometry",
        })

    genes = [row for row in index.get("genes", []) if row.get("gene_symbol") == gene]
    require(len(genes) == 1, f"TSS index must contain exactly one {gene} entry")
    chosen: list[tuple[int, Path]] = []
    covered: set[str] = set()
    seen_pages: set[Path] = set()
    for page in genes[0].get("pages", []):
        ids = page.get("promoter_ids", [])
        require(len(ids) == len(set(ids)), "detail page repeats a promoter ID")
        overlap = [promoter_id for promoter_id in ids if promoter_id in order]
        if not overlap:
            continue
        require(len(overlap) == len(ids),
                "a selected detail page also contains unselected TSS panels")
        svgs = [name for name in page.get("files", []) if name.endswith(".svg")]
        require(len(svgs) == 1, "selected detail page must bind exactly one SVG")
        path = tss_dir / svgs[0]
        require(path.resolve() not in seen_pages, "detail index repeats the same SVG page")
        seen_pages.add(path.resolve())
        require(path.is_file(), f"selected detail SVG is missing: {path}")
        bound_output(tss_receipt, svgs[0], path)
        first = min(order[promoter_id] for promoter_id in overlap)
        chosen.append((first, path))
        # Continuation pages keep index order; FASTA bindings remain one per TSS.
        covered.update(overlap)
    require(covered == set(order), "selected detail pages do not cover the exact selected TSS set")
    chosen.sort(key=lambda row: row[0])
    bindings.sort(key=lambda row: order[row["promoter_id"]])
    return [path for _, path in chosen], bindings


def validate_locus_join(
    locus: dict[str, Any], root: dict[str, str], gene: str,
    report: dict[str, Any], bindings: list[dict[str, Any]],
    locus_sha256: str,
) -> list[dict[str, Any]]:
    require(locus.get("schema") == LOCUS_SVG_SCHEMA and locus.get("gene_symbol") == gene,
            "locus report schema/gene mismatch")
    require(locus.get("panel_id") and root.get("data-gentle-panel-id") == locus["panel_id"],
            "locus SVG/report panel mismatch")
    evidence = locus.get("isoform_evidence", {})
    reference = report.get("reference", {})
    require(evidence.get("assembly")
            and evidence["assembly"] in {reference.get("assembly"), reference.get("genome_id")},
            "locus/TSS reference assembly mismatch")
    chromosome = str(evidence.get("chromosome", "")).removeprefix("chr")
    require(chromosome and all(str(row["chromosome"]).removeprefix("chr") == chromosome
                              and row["strand"] == locus.get("gene_strand") for row in bindings),
            "locus/TSS chromosome or strand mismatch")
    selected_ids = {row["promoter_id"] for row in bindings}
    for window in report.get("windows", []):
        context = window.get("detail_context")
        if window.get("record", {}).get("promoter_id") in selected_ids:
            require(locus.get("transcript_presentation") == (context or {}).get("transcript_presentation"),
                    "overview/detail transcript presentations differ; use the same enriched locus report")
        if context is not None and window.get("record", {}).get("promoter_id") in selected_ids:
            require(context.get("schema") == "gentle.tss_detail_context.v1"
                    and normalized_digest(context.get("locus_report_sha256")) == locus_sha256,
                    "detail context and overview bind different locus reports")
    matrices = report.get("panel_resolution", {}).get("matrices", [])
    require(matrices, "TSS report has no resolved matrix identities")
    matrix_bindings = []
    for matrix in matrices:
        spec = matrix["specification"]
        source = spec["source_id"]
        matches = [track for track in locus.get("regulatory_score_tracks", [])
                   if track.get("provider_kind") == "jaspar_pwm"
                   and track.get("source_ids") == [source]
                   and track.get("provider_version") == source]
        require(len(matches) == 1,
                f"locus overview must contain exactly one resolved {source} matrix track; "
                "regenerate it with accession-preserving scoring")
        matrix_bindings.append({"source_id": source, "locus_track_id": matches[0]["track_id"]})
    return matrix_bindings


def validate_transcript_svg_pages(locus: dict[str, Any], pages: list[Path]) -> str | None:
    presentation = locus.get("transcript_presentation")
    if presentation is None:
        return None
    require(isinstance(presentation, dict) and presentation.get("schema") == "gentle.transcript_structure_presentation.v1",
            "unexpected transcript presentation schema")
    expected = normalized_digest(presentation.get("content_sha256"))
    for page in pages:
        markers = [node for node in ET.fromstring(page.read_bytes()).iter()
                   if node.get("data-role") == "source-coherent-transcripts"]
        require(len(markers) == 1 and markers[0].get("data-content-sha256") == expected,
                "SVG page does not bind the shared transcript presentation")
    return expected


def selected_genbank(index: dict[str, Any], receipt: dict[str, Any],
                     directory: Path, gene: str, bindings: list[dict[str, Any]]) -> tuple[bytes, list[dict[str, str]]]:
    return selected_annotated(index, receipt, directory, gene, bindings, "genbank")


def selected_embl(index: dict[str, Any], receipt: dict[str, Any],
                  directory: Path, gene: str, bindings: list[dict[str, Any]]) -> tuple[bytes, list[dict[str, str]]]:
    return selected_annotated(index, receipt, directory, gene, bindings, "embl")


def selected_annotated(index: dict[str, Any], receipt: dict[str, Any],
                       directory: Path, gene: str, bindings: list[dict[str, Any]],
                       file_format: str) -> tuple[bytes, list[dict[str, str]]]:
    """Concatenate engine-exported records; do not reconstruct feature biology here."""
    label, suffix = {"genbank": ("GenBank", ".gb"), "embl": ("EMBL", ".embl")}[file_format]
    entries = [entry for entry in index["genes"] if entry["gene_symbol"] == gene]
    require(len(entries) == 1, f"{label} requires one indexed gene")
    files = entries[0].get(file_format, {})
    require(isinstance(files, dict), f"{label} index must map promoter IDs to files")
    records, sources = [], []
    for binding in bindings:
        promoter = binding["promoter_id"]
        name = files.get(promoter)
        require(isinstance(name, str) and Path(name).name == name and name.endswith(suffix),
                f"selected TSS lacks indexed {label}; re-export with --formats svg,{file_format} and verified context")
        path = directory / name
        require(not path.is_symlink() and path.is_file(), f"{label} must be a direct regular file")
        bound_output(receipt, name, path)
        raw = path.read_bytes()
        source_digest = hashlib.sha256(raw).hexdigest()
        require(source_digest == normalized_digest(receipt["outputs"][name]),
                f"output hash mismatch for {name}")
        text = raw.decode("utf-8")
        if file_format == "genbank":
            require(text.startswith("LOCUS ") and text.count("\nORIGIN") == 1
                    and text.count("\n//") == 1 and text.rstrip().endswith("//"),
                    "expected one complete engine-exported GenBank record")
            origin = text.split("\nORIGIN", 1)[1].rsplit("\n//", 1)[0]
        else:
            lines = text.rstrip().splitlines()
            sq_lines = [i for i, line in enumerate(lines) if line[:2] == "SQ"]
            require(lines and lines[0].startswith("ID   ")
                    and sum(line[:2] == "ID" for line in lines) == 1
                    and len(sq_lines) == 1 and lines[sq_lines[0]].startswith("SQ   ")
                    and sum(line.strip() == "//" for line in lines) == 1 and lines[-1] == "//",
                    "expected one complete engine-exported EMBL record")
            # SQ is a summary, not sequence: only the following lines carry bases.
            origin = "\n".join(lines[sq_lines[0] + 1:-1])
        bases = re.sub(r"[0-9\s]", "", origin).upper()
        require(bases and set(bases) <= set("ACGTRYSWKMBDHVN")
                and hashlib.sha256(bases.encode("ascii")).hexdigest() == binding["sequence_sha256"],
                f"{label} bases differ from the selected TSS sequence")
        if file_format == "embl":
            records.append(raw if raw.endswith(b"\n") else raw + b"\n")
        else:
            records.append(raw.rstrip() + b"\n")
        sources.append({"promoter_id": promoter, "path": str(path), "sha256": source_digest})
    return b"".join(records), sources


def svg_sequence_bundle(receipt: dict[str, Any], page_paths: list[Path], fasta_text: str,
                        directory: Path, output_zip: Path | None) -> list[tuple[str, bytes, Path]]:
    """Preserve original SVG bytes/title hover and PDF order, with no scientific joins."""
    files: dict[str, bytes] = {}
    pages = []
    for ordinal, path in enumerate(page_paths, 1):
        raw = path.read_bytes()
        digest = hashlib.sha256(raw).hexdigest()
        require(digest == normalized_digest(receipt["inputs"]["pages"][ordinal - 1]["sha256"]),
                "SVG page changed after receipt validation")
        name = f"page_{ordinal:04d}.svg"
        files[name] = raw
        pages.append({"ordinal": ordinal, "file": name, "sha256": digest,
                      "role": receipt["page_order"][ordinal - 1]})
    manifest = {"schema": "gentle.ordered_svg_pages.v1", "pages": pages,
                "gene_symbol": receipt["gene_symbol"], "reference": receipt.get("reference"),
                "locus_report_sha256": receipt["inputs"]["locus_report"]["sha256"],
                "tss_report_sha256": receipt["inputs"]["tss_report"]["sha256"],
                "sequence_file": "selected_tss.fasta",
                "non_claims": receipt.get("non_claims")}
    if receipt["inputs"].get("transcript_presentation_sha256"):
        manifest["transcript_presentation_sha256"] = receipt["inputs"]["transcript_presentation_sha256"]
    files["pages.json"] = (json.dumps(manifest, indent=2, sort_keys=True) + "\n").encode()
    files["selected_tss.fasta"] = fasta_text.encode("ascii")
    title = html.escape(receipt["gene_symbol"] + " - locus and TSS annotations")
    navigation = "".join(f'<li><a href="#page-{p["ordinal"]}">Page {p["ordinal"]}: {html.escape(p["role"])}</a></li>' for p in pages)
    content = "".join(f'<section id="page-{p["ordinal"]}"><h2>Page {p["ordinal"]}: {html.escape(p["role"])}</h2>'
                      f'<a href="{p["file"]}">Open SVG</a><object type="image/svg+xml" data="{p["file"]}" '
                      f'aria-label="Page {p["ordinal"]}"></object></section>' for p in pages)
    files["index.html"] = (f'<!doctype html><html lang="en"><meta charset="utf-8">'
        '<meta name="viewport" content="width=device-width, initial-scale=1">'
        '<meta http-equiv="Content-Security-Policy" content="default-src \'none\'; object-src \'self\'; style-src \'unsafe-inline\'; script-src \'none\'">'
        f'<title>{title}</title><style>body{{font:16px sans-serif;margin:24px}}object{{display:block;width:100%;height:90vh}}'
        'section{break-after:page}h2{font-size:18px}@media print{nav{display:none}}</style>'
        f'<h1>{title}</h1><nav><ol>{navigation}</ol><a href="selected_tss.fasta">Selected sequences (FASTA)</a> | '
        f'<a href="pages.json">Page manifest</a></nav>{content}</html>\n').encode()
    bundle = {"schema": "gentle.ordered_svg_pages.v1", "directory": str(directory),
              "index": "index.html", "pages": pages,
              "files": {name: hashlib.sha256(raw).hexdigest() for name, raw in sorted(files.items())},
              "policy": "Original SVG bytes in exact PDF page order; title hover retained. This is a sequence of SVG files, not a multipage SVG. The outer receipt binds these files; it is not included recursively."}
    outputs = [("svg_bundle", raw, directory / name) for name, raw in sorted(files.items())]
    if output_zip:
        stream = io.BytesIO()
        with zipfile.ZipFile(stream, "w", compression=zipfile.ZIP_DEFLATED) as archive:
            for name, raw in sorted(files.items()):
                info = zipfile.ZipInfo(name, date_time=(1980, 1, 1, 0, 0, 0))
                info.compress_type = zipfile.ZIP_DEFLATED
                info.external_attr = 0o100644 << 16
                archive.writestr(info, raw)
        raw = stream.getvalue()
        bundle["zip"] = {"path": str(output_zip), "sha256": hashlib.sha256(raw).hexdigest()}
        outputs.append(("svg_bundle_zip", raw, output_zip))
    receipt["output"]["svg_bundle"] = bundle
    return outputs


def publish_composite(
    receipt: dict[str, Any], page_paths: list[Path], fasta_text: str,
    gentle_cli: Path, output_pdf: Path, output_fasta: Path, output_receipt: Path,
    genbank_bytes: bytes | None = None, output_genbank: Path | None = None,
    embl_bytes: bytes | None = None, output_embl: Path | None = None,
    output_svg_directory: Path | None = None, output_svg_zip: Path | None = None,
    pdf_representation: str = "raster",
) -> dict[str, Any]:
    require((genbank_bytes is None) == (output_genbank is None), "GenBank content/path must be supplied together")
    require((embl_bytes is None) == (output_embl is None), "EMBL content/path must be supplied together")
    annotated = [(kind, content, path) for kind, content, path in (
        ("genbank", genbank_bytes, output_genbank), ("embl", embl_bytes, output_embl),
    ) if path is not None]
    require(output_svg_zip is None or output_svg_directory is not None, "SVG ZIP requires --output-svg-directory")
    if output_svg_directory:
        require(not os.path.lexists(output_svg_directory), "SVG output directory already exists")
    bundled = (svg_sequence_bundle(receipt, page_paths, fasta_text, output_svg_directory, output_svg_zip)
               if output_svg_directory else [])
    artifacts = annotated + bundled
    outputs = (output_pdf, output_fasta, *(path for _, _, path in artifacts), output_receipt)
    require(len(set(outputs)) == len(outputs), "output paths must be distinct")
    partials = tuple(path.with_name(path.name + ".partial") for path in outputs)
    require(not any(os.path.lexists(path) for path in (*outputs, *partials)),
            "output or stale partial output/receipt exists")
    if output_svg_directory:
        output_svg_directory.mkdir(parents=True, exist_ok=False)
    for path in outputs:
        path.parent.mkdir(parents=True, exist_ok=True)
    owned_partials: list[Path] = []
    published: list[Path] = []
    committed = False
    try:
        for path in partials:
            with path.open("xb"):
                pass
            owned_partials.append(path)
        partial_pdf, partial_fasta, partial_receipt = partials[0], partials[1], partials[-1]
        partial_fasta.write_text(fasta_text, encoding="ascii")
        for (_, content, _), partial in zip(artifacts, partials[2:-1]):
            require(bool(content), "empty annotated sequence export")
            partial.write_bytes(content)
        # PDF and interactive SVG use the same staged bytes when a bundle is requested.
        staged_by_path = dict(zip(outputs, partials))
        render_pages = ([staged_by_path[output_svg_directory / p["file"]]
                         for p in receipt["output"]["svg_bundle"]["pages"]] if bundled else page_paths)
        require(pdf_representation in {"raster", "vector"},
                "PDF representation must be raster or vector")
        renderer_command = ("svg-vector-pdf-set" if pdf_representation == "vector"
                            else "svg-pdf-set")
        command = [str(gentle_cli), renderer_command, str(partial_pdf),
                   *[str(path) for path in render_pages]]
        result = subprocess.run(command, check=True, capture_output=True, text=True, timeout=900)
        summary = json.loads(result.stdout)
        require(summary.get("page_count") == len(page_paths)
                and len(summary.get("pages", [])) == len(page_paths),
                "multi-page renderer returned the wrong page count")
        for page in summary["pages"]:
            require(page.get("font_identity_status") == FONT_IDENTITY_STATUS
                    and isinstance(page.get("font_identities"), list),
                    "multi-page renderer lacks used-font audit; rebuild gentle_cli")
            for font in page["font_identities"]:
                normalized_digest(font.get("sha256"))
                require(font.get("families") and font.get("post_script_name")
                        and isinstance(font.get("face_index"), int) and font["face_index"] >= 0,
                        "multi-page renderer returned an incomplete used-font identity")
        if pdf_representation == "vector":
            require(summary.get("embedded_text") is True
                    and summary.get("svg_interactivity_preserved") is False
                    and summary.get("svg_uri_links_preserved") is False
                    and summary.get("pdf_representation")
                    == "static multipage vector PDF with embedded selectable text",
                    "vector renderer returned an incomplete representation contract")
        require(partial_pdf.is_file() and partial_pdf.stat().st_size > 0,
                "multi-page renderer produced no PDF")
        require(partial_fasta.is_file() and partial_fasta.stat().st_size > 0,
                "selected-TSS FASTA export produced no records")
        receipt["producer"]["renderer_summary"] = summary
        receipt["producer"]["pdf_representation"] = pdf_representation
        receipt["output"].update({
            "pdf_sha256": sha256(partial_pdf), "pdf_bytes": partial_pdf.stat().st_size,
            "selected_tss_fasta_sha256": sha256(partial_fasta),
            "selected_tss_fasta_bytes": partial_fasta.stat().st_size,
        })
        for (kind, _, output), partial in zip(annotated, partials[2:-1]):
            receipt["output"].update({f"selected_tss_{kind}": str(output),
                f"selected_tss_{kind}_sha256": sha256(partial),
                f"selected_tss_{kind}_bytes": partial.stat().st_size})
        partial_receipt.write_text(
            json.dumps(receipt, indent=2, sort_keys=True, ensure_ascii=False) + "\n",
            encoding="utf-8",
        )
        # Publish without overwriting racing outputs. The receipt is the last
        # commit marker; caught failures roll back only files created here.
        for partial, output in zip(partials, outputs):
            os.link(partial, output)
            published.append(output)
        committed = True
    except BaseException:
        for output in reversed(published):
            output.unlink(missing_ok=True)
        raise
    finally:
        for path in owned_partials:
            path.unlink(missing_ok=True)
        if output_svg_directory and not committed:
            try:
                output_svg_directory.rmdir()
            except OSError:
                pass  # Do not remove another writer's content.
    return receipt


def compose(args: argparse.Namespace) -> dict[str, Any]:
    locus_svg = args.locus_svg.resolve()
    locus_receipt_path = args.locus_receipt.resolve()
    locus_report_path = args.locus_report.resolve()
    tss_report_path = args.tss_report.resolve()
    tss_index_path = args.tss_index.resolve()
    tss_receipt_path = args.tss_receipt.resolve()
    tss_manifest_path = args.tss_bundle_manifest.resolve()
    tss_dir = tss_report_path.parent
    output_pdf = args.output_pdf.resolve()
    output_fasta = args.output_fasta.resolve()
    output_receipt = args.output_receipt.resolve()
    output_genbank = getattr(args, "output_genbank", None)
    output_genbank = output_genbank.resolve() if output_genbank else None
    output_embl = getattr(args, "output_embl", None)
    output_embl = output_embl.resolve() if output_embl else None
    gentle_cli = args.gentle_cli.resolve()
    require(len({output_pdf, output_fasta, output_receipt}) == 3,
            "output PDF, FASTA and receipt must be distinct")
    require(not output_pdf.exists() and not output_fasta.exists() and not output_receipt.exists(),
            "outputs must not already exist")
    require(gentle_cli.is_file(), "gentle_cli does not exist")

    locus_receipt = read_json(locus_receipt_path)
    require(locus_receipt.get("schema") == LOCUS_RECEIPT_SCHEMA,
            "unexpected locus receipt schema")
    require(locus_receipt.get("gene") == args.gene, "locus receipt belongs to another gene")
    bound_output(locus_receipt, locus_svg.name, locus_svg)
    root, bands = locus_bands(locus_svg)
    locus_report = read_json(locus_report_path)
    require(normalized_digest(locus_receipt.get("inputs", {}).get("locus_report"))
            == sha256(locus_report_path), "locus report hash mismatch")

    tss_report = read_json(tss_report_path)
    tss_index = read_json(tss_index_path)
    tss_receipt = read_json(tss_receipt_path)
    require(tss_report.get("schema") == TSS_REPORT_SCHEMA, "unexpected TSS report schema")
    require(tss_index.get("schema") == TSS_INDEX_SCHEMA, "unexpected TSS index schema")
    require(tss_receipt.get("schema") == TSS_RECEIPT_SCHEMA, "unexpected TSS receipt schema")
    bound_output(tss_receipt, tss_report_path.name, tss_report_path)
    bound_output(tss_receipt, tss_index_path.name, tss_index_path)
    require(normalized_digest(tss_index.get("report_sha256")) == sha256(tss_report_path),
            "TSS index/report binding mismatch")

    detail_pages, bindings = select_pages(
        args.gene, tss_report, tss_index, tss_dir, tss_receipt, bands,
    )
    locus_matrix_bindings = validate_locus_join(locus_report, root, args.gene, tss_report, bindings,
                                               sha256(locus_report_path))
    for detail_page in detail_pages:
        width, left, right = svg_frame(detail_page)
        require(width == LOCUS_PAGE_WIDTH
                and left == LOCUS_PLOT_LEFT and right == LOCUS_PLOT_RIGHT,
                "selected TSS page does not share the canonical locus horizontal frame")
    fasta_text, fasta_binding = selected_fasta(
        args.gene, tss_manifest_path, tss_receipt, bindings, tss_report.get("reference", {}),
    )
    genbank_bytes, genbank_sources = (selected_genbank(tss_index, tss_receipt, tss_dir, args.gene, bindings)
                                    if output_genbank else (None, []))
    embl_bytes, embl_sources = (selected_embl(tss_index, tss_receipt, tss_dir, args.gene, bindings)
                              if output_embl else (None, []))
    page_paths = [locus_svg, *detail_pages]
    transcript_presentation_sha256 = validate_transcript_svg_pages(locus_report, page_paths)
    receipt = {
        "schema": OUTPUT_SCHEMA,
        "gene_symbol": args.gene,
        "gene_id": bindings[0]["gene_id"],
        "join_policy": "selected promoter_id plus gene/chromosome/TSS/strand; TSS must fall inside a bound locus background band",
        "page_order": ["locus_context", *["selected_tss_tfbs" for _ in detail_pages]],
        "page_count": len(page_paths),
        "selected_tss_count": len(bindings),
        "selected_tss_bindings": bindings,
        "locus_matrix_bindings": locus_matrix_bindings,
        "inputs": {
            "locus_svg": {"path": str(locus_svg), "sha256": sha256(locus_svg)},
            "locus_receipt": {"path": str(locus_receipt_path), "sha256": sha256(locus_receipt_path)},
            "locus_report": {"path": str(locus_report_path), "sha256": sha256(locus_report_path)},
            "tss_report": {"path": str(tss_report_path), "sha256": sha256(tss_report_path)},
            "tss_index": {"path": str(tss_index_path), "sha256": sha256(tss_index_path)},
            "tss_receipt": {"path": str(tss_receipt_path), "sha256": sha256(tss_receipt_path)},
            "pages": [{"path": str(path), "sha256": sha256(path)} for path in page_paths],
        },
        "output": {
            "pdf": str(output_pdf),
            "selected_tss_fasta": str(output_fasta),
        },
        "selected_tss_fasta_binding": fasta_binding,
        "horizontal_alignment": {
            "page_width": LOCUS_PAGE_WIDTH,
            "plot_left": LOCUS_PLOT_LEFT,
            "plot_right": LOCUS_PLOT_RIGHT,
            "policy": "context and detail pages share the canonical locus plot frame",
        },
        "producer": {
            "script": str(Path(__file__).resolve()),
            "script_sha256": sha256(Path(__file__).resolve()),
            "revision": args.producer_revision,
            "gentle_cli": str(gentle_cli),
            "gentle_cli_sha256": sha256(gentle_cli),
        },
        "locus_panel_id": root.get("data-gentle-panel-id"),
        "score_policy": tss_report.get("score_policy"),
        "reference": tss_report.get("reference"),
        "verification": tss_report.get("verification"),
        "non_claims": tss_report.get("non_claims"),
    }
    if output_genbank:
        receipt["inputs"]["selected_tss_genbank"] = genbank_sources
    if output_embl:
        receipt["inputs"]["selected_tss_embl"] = embl_sources
    if transcript_presentation_sha256:
        receipt["inputs"]["transcript_presentation_sha256"] = transcript_presentation_sha256
    return publish_composite(receipt, page_paths, fasta_text, gentle_cli,
                             output_pdf, output_fasta, output_receipt, genbank_bytes, output_genbank,
                             embl_bytes, output_embl,
                             args.output_svg_directory.resolve() if getattr(args, "output_svg_directory", None) else None,
                             args.output_svg_zip.resolve() if getattr(args, "output_svg_zip", None) else None,
                             getattr(args, "pdf_representation", "raster"))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gene", required=True)
    parser.add_argument("--locus-svg", type=Path, required=True)
    parser.add_argument("--locus-receipt", type=Path, required=True)
    parser.add_argument("--locus-report", type=Path, required=True)
    parser.add_argument("--tss-report", type=Path, required=True)
    parser.add_argument("--tss-index", type=Path, required=True)
    parser.add_argument("--tss-receipt", type=Path, required=True)
    parser.add_argument("--tss-bundle-manifest", type=Path, required=True)
    parser.add_argument("--gentle-cli", type=Path, required=True)
    parser.add_argument("--producer-revision", required=True)
    parser.add_argument("--output-pdf", type=Path, required=True)
    parser.add_argument("--pdf-representation", choices=("raster", "vector"), default="raster",
                        help="PDF backend; vector keeps geometry/selectable text but not SVG hover")
    parser.add_argument("--output-fasta", type=Path, required=True)
    parser.add_argument("--output-genbank", type=Path, help="Optional annotated selected windows from indexed, receipt-bound GenBank exports")
    parser.add_argument("--output-embl", type=Path, help="Optional annotated selected windows from indexed, receipt-bound EMBL exports")
    parser.add_argument("--output-receipt", type=Path, required=True)
    parser.add_argument("--output-svg-directory", type=Path, help="New directory for ordered SVG pages, hover-preserving HTML index, sequences and page manifest")
    parser.add_argument("--output-svg-zip", type=Path, help="Optional deterministic ZIP of --output-svg-directory; hash bound by the composite receipt")
    return parser.parse_args()


if __name__ == "__main__":
    try:
        result = compose(parse_args())
        print(json.dumps({
            "status": "ok",
            "schema": result["schema"],
            "gene_symbol": result["gene_symbol"],
            "page_count": result["page_count"],
            "selected_tss_count": result["selected_tss_count"],
            "output": result["output"],
        }, sort_keys=True))
    except Exception as error:
        print(f"error: {error}", file=sys.stderr)
        raise SystemExit(2)
