#!/usr/bin/env python3
"""Build a compact Clariom D probe-set activity presentation.

The E-MTAB-14704 CEL-derived intermediates are intentionally not committed as
bulk data. This script derives an inspectable, reproducible view from explicit
local inputs: an APT raw PM-probe feature table, the Clariom D SQLite platform
annotation, and the vendor probeset annotation ZIP. It performs no downloads,
package installation, or external command execution.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import html
import io
import json
import math
import re
import sqlite3
import statistics
import sys
import zipfile
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
WORK = ROOT / "analysis" / "e_mtab_14704_tp73_microarray"
ANNOTATION = (
    ROOT
    / "test_files"
    / "fixtures"
    / "affymetrix_clariom_d_human_na36_hg38_subset"
    / "clariom_d_human_na36_hg38_gene_panel.probesets.tsv"
)
VENDOR_PROBESET_ZIP = (
    ROOT
    / "data"
    / "publication_resources"
    / "rostock_p73_clariomd_e_mtab_14704"
    / "library"
    / "Clariom_D_Human-na36-hg38-probeset-csv.zip"
)
VENDOR_PROBESET_MEMBER = "Clariom_D_Human.na36.hg38.probeset.csv"
SQLITE = WORK / "Rlib" / "pd.clariom.d.human" / "extdata" / "pd.clariom.d.human.sqlite"
RAW_FEATURES = WORK / "all_arrays_raw_features.tsv"
OUT = WORK / "gene_panel_probe_set_activity"

SCHEMA = "gentle.clariomd_probe_set_activity.v1"
MAX_GENES = 100
MAX_PROBE_SETS = 100_000
MAX_PM_PROBES = 1_000_000
MAX_RAW_FEATURE_ROWS = 10_000_000
MAX_ANNOTATION_LINES = 5_000_000
MAX_INPUT_LINE_CHARS = 1_000_000
MAX_VENDOR_MEMBER_BYTES = 2_000_000_000
SQLITE_QUERY_CHUNK = 500
SVG_HASH_SALT = "gentle-clariomd-probe-set-activity-v1"
FIXED_PDF_DATE = datetime(1970, 1, 1, tzinfo=timezone.utc)

GENES = ["TP73", "FUS", "PATZ1", "E2F1", "TARDBP", "PLK1", "TERT", "HDAC1", "HDAC2", "HDAC6"]
PAIRINGS = (
    {
        "pair_id": "E1",
        "GFP": "P_SKMel29_AdGFP_1.CEL",
        "DNp73beta": "P_SKMel29_AdDNp73beta_1.CEL",
        "TAp73alpha": "P_SKMel29_AdTAp73alpha_1.CEL",
    },
    {
        "pair_id": "E2",
        "GFP": "P_SKMel29_AdGFP_2.CEL",
        "DNp73beta": "P_SKMel29_AdDNp73beta_2.CEL",
        "TAp73alpha": "P_SKMel29_AdTAp73alpha_2.CEL",
    },
    {
        "pair_id": "E3",
        "GFP": "P_SKMel29_AdGFP_3.CEL",
        "DNp73beta": "P_SKMel29_AdDNp73beta_3.CEL",
        "TAp73alpha": "P_SKMel29_AdTAp73alpha_3.CEL",
    },
)
GROUPS = {
    condition: [pairing[condition] for pairing in PAIRINGS]
    for condition in ("DNp73beta", "GFP", "TAp73alpha")
}
SAMPLES = [
    pairing[condition]
    for pairing in PAIRINGS
    for condition in ("GFP", "DNp73beta", "TAp73alpha")
]
CONTRAST_COLUMNS = [
    "log2_TAp73alpha_minus_GFP",
    "log2_DNp73beta_minus_GFP",
    "log2_TAp73alpha_minus_DNp73beta",
]
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

INTERPRETATION_LIMITATIONS = [
    "Raw PM-probe intensities are descriptive and are not a formal normalized expression model.",
    "Displayed contrasts are descriptive effect views, not significance tests.",
    "Probe-set activity does not establish probe specificity or primer binding.",
    "Probe-set or junction activity is not an isoform-support verdict.",
]

CONTRAST_DEFINITIONS = {
    "log2_TAp73alpha_minus_GFP": (
        "mean across declared TAp73alpha arrays of log2(raw PM probeset mean + 1) "
        "minus the corresponding mean across declared GFP arrays"
    ),
    "log2_DNp73beta_minus_GFP": (
        "mean across declared DNp73beta arrays of log2(raw PM probeset mean + 1) "
        "minus the corresponding mean across declared GFP arrays"
    ),
    "log2_TAp73alpha_minus_DNp73beta": (
        "mean across declared TAp73alpha arrays of log2(raw PM probeset mean + 1) "
        "minus the corresponding mean across declared DNp73beta arrays"
    ),
    "paired_condition_minus_GFP": (
        "log2(condition raw PM probeset mean + 1) minus log2(matched GFP raw PM "
        "probeset mean + 1) within the explicitly declared E1/E2/E3 pair"
    ),
}

OUTPUT_FILENAMES = {
    "README.md",
    "gene_contrast_probe_set_summary.pdf",
    "gene_contrast_probe_set_summary.png",
    "gene_contrast_probe_set_summary.svg",
    "index.html",
    "manifest.json",
    "paired_gene_level_summary.tsv",
    "probe_level_activity.tsv",
    "probe_set_activity_heatmap.svg",
    "probe_set_activity_summary.tsv",
    "probe_set_individual_arrays_heatmap.pdf",
    "probe_set_individual_arrays_heatmap.png",
    "probe_set_paired_contrast_heatmap.pdf",
    "probe_set_paired_contrast_heatmap.png",
}


def clean(value: str | None) -> str:
    if value is None or value == "---":
        return ""
    return value


def split_fixture_symbols(value: str | None) -> set[str]:
    symbols: set[str] = set()
    for chunk in (value or "").replace(" /// ", ";").split(";"):
        chunk = chunk.strip()
        if chunk and chunk != "---":
            symbols.add(chunk)
    return symbols


def split_vendor_assignment(value: str | None) -> set[str]:
    symbols: set[str] = set()
    for item in (value or "").split(" /// "):
        parts = [part.strip() for part in item.split(" // ")]
        if len(parts) > 1 and parts[1] and parts[1] != "---":
            symbols.add(parts[1])
    return symbols


def parse_float(value: str) -> float:
    return float(value) if value else math.nan


def mean(values: list[float]) -> float:
    observed = [v for v in values if not math.isnan(v)]
    return sum(observed) / len(observed) if observed else math.nan


def median(values: list[float]) -> float:
    observed = [value for value in values if not math.isnan(value)]
    return float(statistics.median(observed)) if observed else math.nan


def log2(value: float) -> float:
    return math.log2(value + 1.0) if not math.isnan(value) else math.nan


def fmt(value: float) -> str:
    return "" if math.isnan(value) else f"{value:.4f}"


def display_path(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(ROOT.resolve()))
    except ValueError:
        return str(path)


def require_regular_file(path: Path, label: str) -> None:
    if not path.is_file():
        raise RuntimeError(f"{label} is not a readable regular file: {path}")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def provenance_for(path: Path, role: str) -> dict[str, object]:
    require_regular_file(path, role)
    return {
        "role": role,
        "path": display_path(path),
        "size_bytes": path.stat().st_size,
        "sha256": sha256_file(path),
    }


def require_output_separate_from_inputs(input_paths: list[Path]) -> None:
    output_paths = {(OUT / name).resolve() for name in OUTPUT_FILENAMES}
    conflicts = sorted(
        (path for path in input_paths if path.resolve() in output_paths),
        key=lambda path: str(path),
    )
    if conflicts:
        raise RuntimeError(
            "Output directory would overwrite input files: "
            + ", ".join(str(path) for path in conflicts)
        )


def bounded_lines(
    lines,
    *,
    source: Path,
    max_lines: int,
    skip_comments: bool = False,
):
    for physical_line, line in enumerate(lines, start=1):
        if physical_line > max_lines:
            raise RuntimeError(f"Input exceeds the {max_lines}-line limit: {source}")
        if len(line) > MAX_INPUT_LINE_CHARS:
            raise RuntimeError(
                f"Input line exceeds {MAX_INPUT_LINE_CHARS} characters in {source} "
                f"at physical line {physical_line}"
            )
        if skip_comments and line.startswith("#"):
            continue
        yield line


def require_columns(fieldnames: list[str] | None, required: set[str], source: Path) -> None:
    if fieldnames is None:
        raise RuntimeError(f"No header found in {source}")
    duplicates = sorted({name for name in fieldnames if fieldnames.count(name) > 1})
    if duplicates:
        raise RuntimeError(f"Duplicate columns in {source}: {', '.join(duplicates)}")
    missing = sorted(required - set(fieldnames))
    if missing:
        raise RuntimeError(f"Missing required columns in {source}: {', '.join(missing)}")


def parse_raw_intensity(value: str, *, source: Path, row_number: int, sample: str) -> float:
    if not value:
        raise RuntimeError(f"Missing raw intensity for {sample} in {source} row {row_number}")
    try:
        parsed = float(value)
    except ValueError as exc:
        raise RuntimeError(
            f"Invalid raw intensity for {sample} in {source} row {row_number}: {value!r}"
        ) from exc
    if not math.isfinite(parsed) or parsed < 0:
        raise RuntimeError(
            f"Raw intensity must be finite and non-negative for {sample} in "
            f"{source} row {row_number}: {value!r}"
        )
    return parsed


def read_probe_sets() -> tuple[list[dict[str, str]], dict[str, object], list[Path]]:
    wanted = set(GENES)
    rows: list[dict[str, str]] = []
    rows_by_probe_set: dict[str, dict[str, str]] = {}
    annotation_sources: list[Path] = []
    fixture_genes: set[str] = set()

    def append_row(row: dict[str, str | None], symbols: set[str], source: Path) -> None:
        probeset_id = (row.get("probeset_id") or "").strip()
        if not probeset_id:
            raise RuntimeError(f"Matched probe-set row without probeset_id in {source}")
        out = {column: (row.get(column) or "") for column in PROBESET_COLUMNS}
        out["probeset_id"] = probeset_id
        out["gene_symbols"] = ";".join(sorted(symbols, key=GENES.index))
        if existing := rows_by_probe_set.get(probeset_id):
            for column in PROBESET_COLUMNS:
                if column == "gene_symbols":
                    continue
                if existing[column] and out[column] and existing[column] != out[column]:
                    raise RuntimeError(
                        f"Conflicting {column} values for probeset_id {probeset_id!r} in {source}"
                    )
                if not existing[column]:
                    existing[column] = out[column]
            merged_symbols = split_fixture_symbols(existing["gene_symbols"]) | symbols
            existing["gene_symbols"] = ";".join(sorted(merged_symbols, key=GENES.index))
            return
        if len(rows) >= MAX_PROBE_SETS:
            raise RuntimeError(
                f"Selected probe sets exceed the safety limit of {MAX_PROBE_SETS}"
            )
        rows.append(out)
        rows_by_probe_set[probeset_id] = out

    if ANNOTATION.exists():
        require_regular_file(ANNOTATION, "compact probe-set annotation")
        annotation_sources.append(ANNOTATION)
        with ANNOTATION.open(encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(
                bounded_lines(
                    handle,
                    source=ANNOTATION,
                    max_lines=MAX_ANNOTATION_LINES + 1,
                ),
                delimiter="\t",
            )
            require_columns(reader.fieldnames, {"gene_symbols", "probeset_id"}, ANNOTATION)
            for row in reader:
                symbols = split_fixture_symbols(row["gene_symbols"]) & wanted
                if not symbols:
                    continue
                fixture_genes.update(symbols)
                append_row(row, symbols, ANNOTATION)
    missing = wanted - fixture_genes
    vendor_genes: set[str] = set()
    if missing:
        require_regular_file(VENDOR_PROBESET_ZIP, "vendor probe-set annotation ZIP")
        annotation_sources.append(VENDOR_PROBESET_ZIP)
        with zipfile.ZipFile(VENDOR_PROBESET_ZIP) as archive:
            if archive.namelist().count(VENDOR_PROBESET_MEMBER) != 1:
                raise RuntimeError(
                    f"Expected exactly one {VENDOR_PROBESET_MEMBER!r} member in "
                    f"{VENDOR_PROBESET_ZIP}"
                )
            member = archive.getinfo(VENDOR_PROBESET_MEMBER)
            if member.file_size > MAX_VENDOR_MEMBER_BYTES:
                raise RuntimeError(
                    f"Vendor annotation member is {member.file_size} bytes; limit is "
                    f"{MAX_VENDOR_MEMBER_BYTES}: {VENDOR_PROBESET_ZIP}"
                )
            with archive.open(VENDOR_PROBESET_MEMBER) as zipped:
                text = io.TextIOWrapper(zipped, encoding="utf-8-sig", errors="strict", newline="")
                reader = csv.DictReader(
                    bounded_lines(
                        text,
                        source=VENDOR_PROBESET_ZIP,
                        max_lines=MAX_ANNOTATION_LINES + 1,
                        skip_comments=True,
                    )
                )
                require_columns(
                    reader.fieldnames,
                    {"gene_assignment", "probeset_id"},
                    VENDOR_PROBESET_ZIP,
                )
                for row in reader:
                    symbols = split_vendor_assignment(row.get("gene_assignment")) & missing
                    if not symbols:
                        continue
                    vendor_genes.update(symbols)
                    append_row(row, symbols, VENDOR_PROBESET_ZIP)
    unresolved = sorted(wanted - fixture_genes - vendor_genes, key=GENES.index)
    if unresolved:
        raise RuntimeError(
            "No probe-set annotation found for requested genes: " + ", ".join(unresolved)
        )
    if not rows:
        raise RuntimeError("No probe sets matched the requested genes")
    order = {gene: i for i, gene in enumerate(GENES)}
    ordered = sorted(
        rows,
        key=lambda row: (
            min(order[symbol] for symbol in split_fixture_symbols(row["gene_symbols"])),
            int(clean(row.get("start"))) if clean(row.get("start")).isdigit() else sys.maxsize,
            row["probeset_id"],
        ),
    )
    usage = {
        "policy": (
            "The compact fixture is authoritative for each gene it contains; the vendor "
            "ZIP is consulted only for requested genes absent from that fixture."
        ),
        "fixture_genes": sorted(fixture_genes, key=GENES.index),
        "vendor_fallback_genes": sorted(vendor_genes, key=GENES.index),
    }
    return ordered, usage, annotation_sources


def fetch_pm_features(probe_sets: list[dict[str, str]]) -> dict[str, list[dict[str, int]]]:
    wanted = [row["probeset_id"] for row in probe_sets]
    if not wanted:
        raise RuntimeError("Cannot query PM features without selected probe sets")
    require_regular_file(SQLITE, "pd.clariom.d.human SQLite database")
    conn = sqlite3.connect(f"{SQLITE.resolve().as_uri()}?mode=ro", uri=True)
    try:
        by_probe_set: dict[str, list[dict[str, int]]] = defaultdict(list)
        feature_count = 0
        for offset in range(0, len(wanted), SQLITE_QUERY_CHUNK):
            chunk = wanted[offset : offset + SQLITE_QUERY_CHUNK]
            placeholders = ",".join("?" for _ in chunk)
            query = f"""
                SELECT fs.man_fsetid, pm.fid, pm.x, pm.y
                FROM featureSet fs
                JOIN pmfeature pm ON pm.fsetid = fs.fsetid
                WHERE fs.man_fsetid IN ({placeholders})
                ORDER BY fs.man_fsetid, pm.fid
            """
            for probeset_id, fid, x, y in conn.execute(query, chunk):
                feature_count += 1
                if feature_count > MAX_PM_PROBES:
                    raise RuntimeError(
                        f"Selected PM probes exceed the safety limit of {MAX_PM_PROBES}"
                    )
                by_probe_set[str(probeset_id)].append(
                    {"probe_id": int(fid), "x": int(x), "y": int(y)}
                )
        missing = [probeset_id for probeset_id in wanted if not by_probe_set[probeset_id]]
        if missing:
            preview = ", ".join(missing[:10])
            suffix = "..." if len(missing) > 10 else ""
            raise RuntimeError(
                f"SQLite PM-feature mapping is missing for {len(missing)} selected probe sets: "
                f"{preview}{suffix}"
            )
        return dict(by_probe_set)
    except sqlite3.Error as exc:
        raise RuntimeError(f"Could not read expected Clariom D tables from {SQLITE}: {exc}") from exc
    finally:
        conn.close()


def read_raw_intensities(probe_ids: set[int]) -> dict[int, dict[str, float]]:
    if not probe_ids:
        raise RuntimeError("Cannot read raw intensities without selected PM probe ids")
    require_regular_file(RAW_FEATURES, "APT raw PM-probe feature table")
    intensities: dict[int, dict[str, float]] = {}
    with RAW_FEATURES.open(encoding="utf-8-sig", newline="") as handle:
        header: list[str] | None = None
        for physical_line, line in enumerate(handle, start=1):
            if physical_line > MAX_RAW_FEATURE_ROWS + 1:
                raise RuntimeError(
                    f"Raw feature table exceeds the {MAX_RAW_FEATURE_ROWS}-row scan limit "
                    f"before its header: {RAW_FEATURES}"
                )
            if len(line) > MAX_INPUT_LINE_CHARS:
                raise RuntimeError(
                    f"Input line exceeds {MAX_INPUT_LINE_CHARS} characters in "
                    f"{RAW_FEATURES} at physical line {physical_line}"
                )
            if line.startswith("#"):
                continue
            header = line.rstrip("\r\n").split("\t")
            break
        if header is None:
            raise RuntimeError(f"No header found in {RAW_FEATURES}")
        require_columns(header, {"probe_id", *SAMPLES}, RAW_FEATURES)
        index = {name: i for i, name in enumerate(header)}
        for row_number, line in enumerate(handle, start=1):
            if row_number > MAX_RAW_FEATURE_ROWS:
                raise RuntimeError(
                    f"Raw feature table exceeds the {MAX_RAW_FEATURE_ROWS}-row scan limit: "
                    f"{RAW_FEATURES}"
                )
            if len(line) > MAX_INPUT_LINE_CHARS:
                raise RuntimeError(
                    f"Input line exceeds {MAX_INPUT_LINE_CHARS} characters in "
                    f"{RAW_FEATURES} at data row {row_number}"
                )
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != len(header):
                raise RuntimeError(
                    f"Expected {len(header)} tab-separated fields in {RAW_FEATURES} data row "
                    f"{row_number}, found {len(fields)}"
                )
            try:
                probe_id = int(fields[index["probe_id"]])
            except ValueError as exc:
                raise RuntimeError(
                    f"Invalid probe_id in {RAW_FEATURES} data row {row_number}: "
                    f"{fields[index['probe_id']]!r}"
                ) from exc
            if probe_id not in probe_ids:
                continue
            if probe_id in intensities:
                raise RuntimeError(
                    f"Duplicate selected probe_id {probe_id} in {RAW_FEATURES} data row {row_number}"
                )
            intensities[probe_id] = {
                sample: parse_raw_intensity(
                    fields[index[sample]],
                    source=RAW_FEATURES,
                    row_number=row_number,
                    sample=sample,
                )
                for sample in SAMPLES
            }
            if len(intensities) == len(probe_ids):
                break
    return intensities


def write_tsv(path: Path, rows: list[dict[str, str | int | float]], fieldnames: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def summarize_probe_sets(
    probe_sets: list[dict[str, str]],
    pm_features: dict[str, list[dict[str, int]]],
    intensities: dict[int, dict[str, float]],
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    probe_level_rows: list[dict[str, str]] = []
    summary_rows: list[dict[str, str]] = []
    for row in probe_sets:
        probeset_id = row["probeset_id"]
        gene_symbols = split_fixture_symbols(row["gene_symbols"])
        primary_gene = min(gene_symbols, key=GENES.index)
        features = pm_features.get(probeset_id, [])
        observed = [feature for feature in features if feature["probe_id"] in intensities]
        for feature in observed:
            values = intensities[feature["probe_id"]]
            probe_level_rows.append(
                {
                    "gene": row["gene_symbols"],
                    "primary_gene": primary_gene,
                    "probeset_id": probeset_id,
                    "probe_id": str(feature["probe_id"]),
                    "x": str(feature["x"]),
                    "y": str(feature["y"]),
                    **{sample: fmt(values[sample]) for sample in SAMPLES},
                }
            )
        sample_means = {
            sample: mean([intensities[feature["probe_id"]][sample] for feature in observed])
            for sample in SAMPLES
        }
        group_log2 = {
            group: mean([log2(sample_means[sample]) for sample in samples])
            for group, samples in GROUPS.items()
        }
        summary_rows.append(
            {
                "gene": row["gene_symbols"],
                "primary_gene": primary_gene,
                "probeset_id": probeset_id,
                "probeset_type": row["probeset_type"],
                "seqname": row["seqname"],
                "strand": row["strand"],
                "start": clean(row["start"]),
                "stop": clean(row["stop"]),
                "transcript_cluster_id": row["transcript_cluster_id"],
                "exon_id": clean(row["exon_id"]),
                "junction_start_edge": clean(row["junction_start_edge"]),
                "junction_stop_edge": clean(row["junction_stop_edge"]),
                "annotated_probe_count": row["probe_count"],
                "pm_probe_count": str(len(features)),
                "observed_pm_probe_count": str(len(observed)),
                **{sample: fmt(sample_means[sample]) for sample in SAMPLES},
                "log2_mean_DNp73beta": fmt(group_log2["DNp73beta"]),
                "log2_mean_GFP": fmt(group_log2["GFP"]),
                "log2_mean_TAp73alpha": fmt(group_log2["TAp73alpha"]),
                "log2_TAp73alpha_minus_GFP": fmt(group_log2["TAp73alpha"] - group_log2["GFP"]),
                "log2_DNp73beta_minus_GFP": fmt(group_log2["DNp73beta"] - group_log2["GFP"]),
                "log2_TAp73alpha_minus_DNp73beta": fmt(
                    group_log2["TAp73alpha"] - group_log2["DNp73beta"]
                ),
            }
        )
    return summary_rows, probe_level_rows


def colour(value: float, lo: float, hi: float) -> str:
    if math.isnan(value):
        return "#f3f4f6"
    if hi <= lo:
        t = 0.5
    else:
        t = max(0.0, min(1.0, (value - lo) / (hi - lo)))
    # Blue-white-red without a single-hue palette.
    if t < 0.5:
        f = t / 0.5
        r = int(49 + (245 - 49) * f)
        g = int(104 + (245 - 104) * f)
        b = int(176 + (245 - 176) * f)
    else:
        f = (t - 0.5) / 0.5
        r = int(245 + (190 - 245) * f)
        g = int(245 + (55 - 245) * f)
        b = int(245 + (45 - 245) * f)
    return f"#{r:02x}{g:02x}{b:02x}"


def write_svg(rows: list[dict[str, str]]) -> None:
    row_h = 14
    label_w = 292
    cell_w = 26
    top = 88
    columns = SAMPLES + CONTRAST_COLUMNS
    abundance_values = [
        log2(parse_float(row[column]))
        for row in rows
        for column in SAMPLES
        if row.get(column)
    ]
    contrast_values = [
        parse_float(row[column])
        for row in rows
        for column in CONTRAST_COLUMNS
        if row.get(column)
    ]
    if not abundance_values or not contrast_values:
        raise RuntimeError("Cannot render an activity heatmap without finite abundance and contrast values")
    abundance_lo, abundance_hi = min(abundance_values), max(abundance_values)
    contrast_limit = max(1.0, *(abs(value) for value in contrast_values))
    width = label_w + cell_w * len(columns) + 250
    height = max(214, top + row_h * len(rows) + 84)
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}" role="img">',
        "<title>TP73 overexpression probe-set raw PM activity</title>",
        "<style>"
        "text{font-family:Inter,Arial,sans-serif;font-size:10px;fill:#111827}"
        ".small{font-size:9px;fill:#374151}.tiny{font-size:8px;fill:#4b5563}"
        ".gene{font-weight:700}.axis{stroke:#d1d5db;stroke-width:1}"
        "</style>",
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        '<text x="18" y="24" style="font-size:16px;font-weight:700">Probe-set activity in TP73 isoform overexpression arrays</text>',
        '<text x="18" y="43" class="small">Descriptive raw PM-probe means and contrasts; not a formal expression model, significance test, specificity verdict, or isoform-support call.</text>',
    ]
    for i, column in enumerate(columns):
        x = label_w + i * cell_w + 18
        parts.append(
            f'<text x="{x}" y="78" class="tiny" transform="rotate(-45 {x} 78)">'
            f"{html.escape(column.replace('P_SKMel29_Ad', '').replace('.CEL', ''))}</text>"
        )
    last_gene = None
    contrast_start_x = label_w + len(SAMPLES) * cell_w
    parts.append(
        f'<line x1="{contrast_start_x - 3}" y1="50" x2="{contrast_start_x - 3}" '
        f'y2="{height - 62}" stroke="#111827" stroke-width="1.2"/>'
    )
    for r, row in enumerate(rows):
        y = top + r * row_h
        gene = row["gene"]
        if gene != last_gene:
            parts.append(f'<line x1="16" y1="{y - 5}" x2="{width - 18}" y2="{y - 5}" class="axis"/>')
            parts.append(f'<text x="18" y="{y + 9}" class="gene">{html.escape(gene)}</text>')
            last_gene = gene
        label = f"{row['probeset_id']} {row['probeset_type']} {row['start']}-{row['stop']}".strip()
        parts.append(f'<text x="66" y="{y + 9}" class="tiny">{html.escape(label[:56])}</text>')
        for c, column in enumerate(columns):
            x = label_w + c * cell_w
            value = parse_float(row[column]) if row[column] else math.nan
            if column in SAMPLES:
                displayed_value = log2(value)
                fill = colour(displayed_value, abundance_lo, abundance_hi)
                title = (
                    f"{row['probeset_id']} {column}: raw mean {fmt(value)}; "
                    f"log2(raw mean + 1) {fmt(displayed_value)}"
                )
            else:
                displayed_value = value
                fill = colour(displayed_value, -contrast_limit, contrast_limit)
                title = f"{row['probeset_id']} {column}: {fmt(displayed_value)}"
            parts.append(
                f'<rect x="{x}" y="{y}" width="{cell_w - 1}" height="{row_h - 1}" '
                f'fill="{fill}"><title>{html.escape(title)}</title></rect>'
            )
    legend_x = label_w + cell_w * len(columns) + 34
    parts.extend(
        [
            f'<text x="{legend_x}" y="92" class="small">Abundance: log2(raw PM mean + 1)</text>',
            f'<text x="{legend_x}" y="110" class="tiny">min {abundance_lo:.2f}</text>',
            f'<text x="{legend_x + 92}" y="110" class="tiny">max {abundance_hi:.2f}</text>',
        ]
    )
    for i in range(80):
        value = abundance_lo + (abundance_hi - abundance_lo) * i / 79
        parts.append(
            f'<rect x="{legend_x + i}" y="116" width="1" height="12" '
            f'fill="{colour(value, abundance_lo, abundance_hi)}"/>'
        )
    parts.extend(
        [
            f'<text x="{legend_x}" y="154" class="small">Contrasts: log2 difference</text>',
            f'<text x="{legend_x}" y="172" class="tiny">-{contrast_limit:.2f}</text>',
            f'<text x="{legend_x + 92}" y="172" class="tiny">+{contrast_limit:.2f}</text>',
        ]
    )
    for i in range(80):
        value = -contrast_limit + 2 * contrast_limit * i / 79
        parts.append(
            f'<rect x="{legend_x + i}" y="178" width="1" height="12" '
            f'fill="{colour(value, -contrast_limit, contrast_limit)}"/>'
        )
    parts.append("</svg>")
    (OUT / "probe_set_activity_heatmap.svg").write_text(
        "\n".join(parts) + "\n", encoding="utf-8"
    )


def write_html(rows: list[dict[str, str]], meta: dict[str, object]) -> None:
    cards = []
    for gene in GENES:
        gene_rows = [row for row in rows if row["primary_gene"] == gene]
        top_ta = sorted(
            gene_rows,
            key=lambda row: parse_float(row["log2_TAp73alpha_minus_GFP"]),
            reverse=True,
        )[:5]
        cards.append(f"<section><h2>{html.escape(gene)}</h2>")
        cards.append(
            f"<p>{len(gene_rows)} probesets, "
            f"{sum(int(r['observed_pm_probe_count']) for r in gene_rows)} observed PM probes.</p>"
        )
        cards.append("<table><thead><tr><th>probeset</th><th>type</th><th>TA-GFP</th><th>DN-GFP</th><th>locus</th></tr></thead><tbody>")
        for row in top_ta:
            locus = f"{row['seqname']}:{row['start']}-{row['stop']}" if row["start"] else "junction"
            cards.append(
                "<tr>"
                f"<td>{html.escape(row['probeset_id'])}</td>"
                f"<td>{html.escape(row['probeset_type'])}</td>"
                f"<td>{html.escape(row['log2_TAp73alpha_minus_GFP'])}</td>"
                f"<td>{html.escape(row['log2_DNp73beta_minus_GFP'])}</td>"
                f"<td>{html.escape(locus)}</td>"
                "</tr>"
            )
        cards.append("</tbody></table></section>")
    html_text = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>TP73 overexpression probe-set activity</title>
<style>
body {{ font-family: Inter, Arial, sans-serif; margin: 28px; color: #111827; background: #fff; }}
h1 {{ font-size: 24px; margin-bottom: 6px; }}
h2 {{ font-size: 17px; margin: 22px 0 6px; }}
p {{ max-width: 980px; color: #374151; }}
table {{ border-collapse: collapse; width: 100%; max-width: 1120px; font-size: 13px; }}
th, td {{ border-bottom: 1px solid #e5e7eb; padding: 5px 7px; text-align: left; }}
th {{ background: #f9fafb; }}
.figure {{ margin: 18px 0; border: 1px solid #e5e7eb; overflow-x: auto; }}
.meta {{ font-size: 12px; color: #4b5563; }}
</style>
</head>
<body>
<h1>Probe-set activity in TP73 isoform overexpression arrays</h1>
<p>This is a descriptive figure-preparation surface for raw PM-probe activity summarized by Clariom D probeset across three explicitly paired GFP controls, three DNp73beta arrays, and three TAp73alpha arrays. It is not a formal normalized expression model, significance test, probe-specificity or primer-binding verdict, or isoform-support claim.</p>
<p class="meta">Generated from {html.escape(display_path(RAW_FEATURES))} and {html.escape(display_path(SQLITE))}. {html.escape(json.dumps(meta, sort_keys=True))}</p>
<div class="figure"><img src="probe_set_activity_heatmap.svg" alt="Probe-set activity heatmap"></div>
{''.join(cards)}
</body>
</html>
"""
    (OUT / "index.html").write_text(html_text, encoding="utf-8")


def paired_log2(row: dict[str, str], pairing: dict[str, str], condition: str) -> float:
    control = log2(parse_float(row[pairing["GFP"]]))
    treated = log2(parse_float(row[pairing[condition]]))
    return treated - control


def write_paired_gene_summary(rows: list[dict[str, str]]) -> None:
    fieldnames = ["primary_gene", "probesets", "pm_probes"]
    for pairing in PAIRINGS:
        pair_id = pairing["pair_id"]
        fieldnames.extend(
            [
                f"{pair_id}_median_DN_minus_GFP",
                f"{pair_id}_median_TA_minus_GFP",
            ]
        )
    fieldnames.extend(
        [
            "median_DN_minus_GFP_across_pairs",
            "median_TA_minus_GFP_across_pairs",
        ]
    )
    with (OUT / "paired_gene_level_summary.tsv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for gene in GENES:
            gene_rows = [row for row in rows if row["primary_gene"] == gene]
            if not gene_rows:
                continue
            out: dict[str, str | int] = {
                "primary_gene": gene,
                "probesets": len(gene_rows),
                "pm_probes": sum(int(row["observed_pm_probe_count"]) for row in gene_rows),
            }
            dn_medians: list[float] = []
            ta_medians: list[float] = []
            for pairing in PAIRINGS:
                pair_id = pairing["pair_id"]
                dn = [paired_log2(row, pairing, "DNp73beta") for row in gene_rows]
                ta = [paired_log2(row, pairing, "TAp73alpha") for row in gene_rows]
                dn_median = median(dn)
                ta_median = median(ta)
                out[f"{pair_id}_median_DN_minus_GFP"] = fmt(dn_median)
                out[f"{pair_id}_median_TA_minus_GFP"] = fmt(ta_median)
                dn_medians.append(dn_median)
                ta_medians.append(ta_median)
            out["median_DN_minus_GFP_across_pairs"] = fmt(median(dn_medians))
            out["median_TA_minus_GFP_across_pairs"] = fmt(median(ta_medians))
            writer.writerow(out)


def save_matplotlib_figure(fig, output: Path, *, title: str, dpi: int | None = None) -> None:
    creator = "GENtle render_clariomd_probe_set_activity.py"
    if output.suffix == ".svg":
        metadata = {
            "Title": title,
            "Creator": creator,
            "Description": INTERPRETATION_LIMITATIONS[0],
            "Date": None,
        }
    elif output.suffix == ".pdf":
        metadata = {
            "Title": title,
            "Subject": "; ".join(INTERPRETATION_LIMITATIONS),
            "Creator": creator,
            "Producer": "Matplotlib via GENtle",
            "CreationDate": FIXED_PDF_DATE,
            "ModDate": FIXED_PDF_DATE,
        }
    else:
        metadata = {"Software": creator, "Description": INTERPRETATION_LIMITATIONS[0]}
    fig.savefig(output, dpi=dpi, metadata=metadata)


def load_matplotlib_stack():
    try:
        import matplotlib

        matplotlib.use("Agg")
        matplotlib.rcParams["svg.hashsalt"] = SVG_HASH_SALT
        import matplotlib.pyplot as plt
        import numpy as np
        from matplotlib.colors import TwoSlopeNorm
    except (ImportError, ModuleNotFoundError) as exc:  # pragma: no cover - host dependency.
        raise RuntimeError(
            "Matplotlib and NumPy are required to render PNG/SVG/PDF outputs; "
            "install them explicitly before running this non-installing helper"
        ) from exc
    return (plt, np, TwoSlopeNorm), {
        "matplotlib": matplotlib.__version__,
        "numpy": np.__version__,
        "python": sys.version.split()[0],
    }


def write_matplotlib_figures(rows: list[dict[str, str]], plotting_stack) -> list[str]:
    plt, np, TwoSlopeNorm = plotting_stack

    written: list[str] = []
    ordered: list[dict[str, str]] = []
    for gene in GENES:
        gene_rows = [row for row in rows if row["primary_gene"] == gene]
        gene_rows.sort(
            key=lambda row: mean(
                [paired_log2(row, pairing, "TAp73alpha") for pairing in PAIRINGS]
            ),
            reverse=True,
        )
        ordered.extend(gene_rows)

    # Gene-level contrast distribution.
    fig, axes = plt.subplots(2, 1, figsize=(12, 8), sharex=True, constrained_layout=True)
    for ax, (column, title, colour_name) in zip(
        axes,
        [
            ("log2_TAp73alpha_minus_GFP", "TAp73alpha vs GFP", "#b91c1c"),
            ("log2_DNp73beta_minus_GFP", "DNp73beta vs GFP", "#1d4ed8"),
        ],
    ):
        for idx, gene in enumerate(GENES):
            values = [parse_float(row[column]) for row in rows if row["primary_gene"] == gene]
            xs = [idx + ((n % 11) - 5) * 0.025 for n in range(len(values))]
            ax.scatter(xs, values, s=18, c=colour_name, alpha=0.55, edgecolors="none")
            if values:
                median_value = median(values)
                ax.plot(
                    [idx - 0.3, idx + 0.3],
                    [median_value, median_value],
                    color="black",
                    linewidth=2,
                )
                ax.text(
                    idx,
                    max(values) + (0.35 if max(values) > 2 else 0.22),
                    f"n={len(values)}",
                    ha="center",
                    va="bottom",
                    fontsize=8,
                    color="#374151",
                )
        ax.axhline(0, color="#6b7280", linewidth=1)
        ax.axhline(2, color="#dc2626", linewidth=1, linestyle="--")
        ax.axhline(-2, color="#2563eb", linewidth=1, linestyle="--")
        ax.set_ylabel("log2 contrast")
        ax.set_title(title, loc="left", fontsize=13, fontweight="bold")
        ax.grid(axis="y", color="#e5e7eb", linewidth=0.8)
        ax.set_ylim(-2.5, 8.1)
    axes[-1].set_xticks(range(len(GENES)))
    axes[-1].set_xticklabels(GENES, rotation=35, ha="right")
    fig.suptitle("Probe-set activity by gene: raw PM-probe means across 9 arrays", fontsize=15, fontweight="bold")
    fig.text(
        0.01,
        0.01,
        "Each point is one Clariom D probeset; black segment is median. Dashed +/-2 guides are descriptive, not significance thresholds.",
        fontsize=9,
        color="#374151",
    )
    for extension in ("png", "svg", "pdf"):
        output = OUT / f"gene_contrast_probe_set_summary.{extension}"
        save_matplotlib_figure(
            fig,
            output,
            title="Descriptive probe-set contrasts by gene",
            dpi=180 if extension == "png" else None,
        )
        written.append(output.name)
    plt.close(fig)

    height = max(10, len(ordered) * 0.045)
    sample_columns: list[tuple[str, str]] = []
    for pairing in PAIRINGS:
        sample_columns.extend(
            [
                (pairing["GFP"], f"{pairing['pair_id']} GFP"),
                (pairing["DNp73beta"], f"{pairing['pair_id']} DN"),
                (pairing["TAp73alpha"], f"{pairing['pair_id']} TA"),
            ]
        )
    matrix = np.array([[log2(parse_float(row[column])) for column, _ in sample_columns] for row in ordered])
    fig, ax = plt.subplots(figsize=(9.2, height), constrained_layout=True)
    image = ax.imshow(matrix, aspect="auto", cmap="magma")
    ax.set_xticks(range(len(sample_columns)))
    ax.set_xticklabels([label for _, label in sample_columns], rotation=45, ha="right", fontsize=8)
    ax.set_yticks([])
    ax.set_title(
        "Individual array probe-set activity: rows sorted within gene by mean paired TA-GFP",
        loc="left",
        fontsize=13,
        fontweight="bold",
    )
    add_gene_separators(ax, ordered, text_x=-0.8, line_color="white")
    for x in (2.5, 5.5):
        ax.axvline(x, color="white", linewidth=1.2)
    colourbar = fig.colorbar(image, ax=ax, shrink=0.8, pad=0.03)
    colourbar.set_label("log2 raw PM-probe mean")
    fig.text(
        0.02,
        0.005,
        "Row order: genes follow --genes; probesets are sorted within each gene by mean paired TA-GFP. Columns: E1/E2/E3 as GFP, DNp73beta, TAp73alpha.",
        fontsize=8,
    )
    for extension in ("png", "pdf"):
        output = OUT / f"probe_set_individual_arrays_heatmap.{extension}"
        save_matplotlib_figure(
            fig,
            output,
            title="Individual-array descriptive probe-set activity",
            dpi=220 if extension == "png" else None,
        )
        written.append(output.name)
    plt.close(fig)

    paired_matrix = np.array(
        [
            [
                value
                for pairing in PAIRINGS
                for value in (
                    paired_log2(row, pairing, "DNp73beta"),
                    paired_log2(row, pairing, "TAp73alpha"),
                )
            ]
            for row in ordered
        ]
    )
    paired_labels = [
        label
        for pairing in PAIRINGS
        for label in (
            f"{pairing['pair_id']} DN-GFP",
            f"{pairing['pair_id']} TA-GFP",
        )
    ]
    fig, ax = plt.subplots(figsize=(8.5, height), constrained_layout=True)
    image = ax.imshow(
        paired_matrix,
        aspect="auto",
        cmap="RdBu_r",
        norm=TwoSlopeNorm(vmin=-2.5, vcenter=0, vmax=7.5),
    )
    ax.set_xticks(range(len(paired_labels)))
    ax.set_xticklabels(paired_labels, rotation=45, ha="right", fontsize=8)
    ax.set_yticks([])
    ax.set_title(
        "Paired within-experiment probe-set contrasts: same row order as individual arrays",
        loc="left",
        fontsize=13,
        fontweight="bold",
    )
    add_gene_separators(ax, ordered, text_x=-0.7, line_color="black")
    for x in (1.5, 3.5):
        ax.axvline(x, color="black", linewidth=0.8)
    colourbar = fig.colorbar(image, ax=ax, shrink=0.8, pad=0.03)
    colourbar.set_label("paired log2 contrast")
    fig.text(
        0.02,
        0.005,
        "Row order: genes follow --genes; probesets sorted within gene by mean paired TA-GFP. Contrast: log2(condition PM mean + 1) - log2(matched GFP PM mean + 1).",
        fontsize=8,
    )
    for extension in ("png", "pdf"):
        output = OUT / f"probe_set_paired_contrast_heatmap.{extension}"
        save_matplotlib_figure(
            fig,
            output,
            title="Explicitly paired descriptive probe-set contrasts",
            dpi=220 if extension == "png" else None,
        )
        written.append(output.name)
    plt.close(fig)
    return written


def add_gene_separators(ax, ordered_rows: list[dict[str, str]], text_x: float, line_color: str) -> None:
    start = 0
    for gene in GENES:
        count = sum(1 for row in ordered_rows if row["primary_gene"] == gene)
        if not count:
            continue
        ax.axhline(start - 0.5, color=line_color, linewidth=0.7)
        ax.text(text_x, start + count / 2 - 0.5, f"{gene} ({count})", va="center", ha="right", fontsize=8)
        start += count
    ax.axhline(start - 0.5, color=line_color, linewidth=0.7)


def write_readme(meta: dict[str, object]) -> None:
    inputs = meta.get("inputs")
    if not isinstance(inputs, list):
        raise RuntimeError("Renderer metadata is missing its input provenance list")
    input_lines = []
    for source in inputs:
        if not isinstance(source, dict):
            raise RuntimeError("Renderer input provenance contains a non-object entry")
        input_lines.append(
            f"- `{source['path']}`: {source['role']}; {source['size_bytes']} bytes; "
            f"SHA-256 `{source['sha256']}`."
        )
    input_text = "\n".join(input_lines)
    text = f"""# TP73 Overexpression Probe-Set Activity

Generated local presentation for {", ".join(GENES)}.

Inputs:

{input_text}

Outputs:

- `index.html`: compact presentation page.
- `probe_set_activity_heatmap.svg`: figure-ready heatmap overview.
- `probe_set_activity_summary.tsv`: probeset-level mean raw intensity and log2 group contrasts.
- `probe_level_activity.tsv`: selected PM-probe raw intensities.
- `gene_contrast_probe_set_summary.png/.svg/.pdf`: compact per-gene contrast distribution.
- `probe_set_individual_arrays_heatmap.png/.pdf`: individual-array heatmap. Rows follow the `--genes` order and are sorted within each gene by mean paired `TAp73alpha_i - GFP_i`; columns are `E1 GFP/DN/TA`, `E2 GFP/DN/TA`, `E3 GFP/DN/TA`.
- `probe_set_paired_contrast_heatmap.png/.pdf`: within-experiment paired contrast heatmap.
- `paired_gene_level_summary.tsv`: per-gene medians for the paired contrasts.
- `manifest.json`: machine-readable provenance.

Caveat: the activity values are raw PM-probe intensities summarized by probeset. They are useful for visual inspection and figure preparation, but are not a formal normalized expression model, significance test, probe-specificity or primer-binding verdict, or isoform-support claim.

Metadata:

```json
{json.dumps(meta, indent=2, sort_keys=True)}
```
"""
    (OUT / "README.md").write_text(text, encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--raw-features",
        type=Path,
        default=RAW_FEATURES,
        help="APT raw PM-probe feature table with probe_id/x/y and one column per CEL file.",
    )
    parser.add_argument(
        "--sqlite",
        type=Path,
        default=SQLITE,
        help="pd.clariom.d.human SQLite database containing featureSet and pmfeature.",
    )
    parser.add_argument(
        "--fixture-probesets",
        type=Path,
        default=ANNOTATION,
        help="Optional compact GENtle probeset fixture to use before falling back to the vendor ZIP.",
    )
    parser.add_argument(
        "--vendor-probeset-zip",
        type=Path,
        default=VENDOR_PROBESET_ZIP,
        help="Clariom D Human na36 hg38 probeset CSV ZIP for genes missing from the compact fixture.",
    )
    parser.add_argument("--output-dir", type=Path, default=OUT, help="Directory for generated local outputs.")
    parser.add_argument(
        "--genes",
        default=",".join(GENES),
        help="Comma-separated gene symbols to render, in display order.",
    )
    return parser.parse_args()


def configure_from_args(args: argparse.Namespace) -> None:
    global RAW_FEATURES, SQLITE, ANNOTATION, VENDOR_PROBESET_ZIP, OUT, GENES
    RAW_FEATURES = args.raw_features
    SQLITE = args.sqlite
    ANNOTATION = args.fixture_probesets
    VENDOR_PROBESET_ZIP = args.vendor_probeset_zip
    OUT = args.output_dir
    genes = [gene.strip() for gene in args.genes.split(",") if gene.strip()]
    if not genes:
        raise RuntimeError("--genes must contain at least one gene symbol")
    if len(genes) > MAX_GENES:
        raise RuntimeError(f"--genes exceeds the safety limit of {MAX_GENES} symbols")
    duplicates = sorted({gene for gene in genes if genes.count(gene) > 1})
    if duplicates:
        raise RuntimeError(f"--genes contains duplicate symbols: {', '.join(duplicates)}")
    invalid = [gene for gene in genes if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]{0,63}", gene)]
    if invalid:
        raise RuntimeError(f"Invalid gene symbols in --genes: {', '.join(invalid)}")
    GENES = genes


def run(args: argparse.Namespace) -> None:
    configure_from_args(args)
    csv.field_size_limit(MAX_INPUT_LINE_CHARS)
    probe_sets, annotation_usage, annotation_sources = read_probe_sets()
    pm_features = fetch_pm_features(probe_sets)
    probe_ids = {feature["probe_id"] for features in pm_features.values() for feature in features}
    intensities = read_raw_intensities(probe_ids)
    missing = sorted(probe_ids - set(intensities))
    if missing:
        raise RuntimeError(f"Missing {len(missing)} PM probe ids in raw feature table")
    summary_rows, probe_level_rows = summarize_probe_sets(probe_sets, pm_features, intensities)
    input_paths = [RAW_FEATURES, SQLITE, *annotation_sources]
    require_output_separate_from_inputs(input_paths)
    input_provenance = [
        provenance_for(RAW_FEATURES, "APT raw PM-probe feature table"),
        provenance_for(SQLITE, "pd.clariom.d.human SQLite database"),
        *(provenance_for(path, "probe-set annotation") for path in annotation_sources),
    ]
    plotting_stack, render_runtime = load_matplotlib_stack()
    OUT.mkdir(parents=True, exist_ok=True)
    summary_fields = [
        "gene",
        "primary_gene",
        "probeset_id",
        "probeset_type",
        "seqname",
        "strand",
        "start",
        "stop",
        "transcript_cluster_id",
        "exon_id",
        "junction_start_edge",
        "junction_stop_edge",
        "annotated_probe_count",
        "pm_probe_count",
        "observed_pm_probe_count",
        *SAMPLES,
        "log2_mean_DNp73beta",
        "log2_mean_GFP",
        "log2_mean_TAp73alpha",
        *CONTRAST_COLUMNS,
    ]
    probe_fields = ["gene", "primary_gene", "probeset_id", "probe_id", "x", "y", *SAMPLES]
    write_tsv(OUT / "probe_set_activity_summary.tsv", summary_rows, summary_fields)
    write_tsv(OUT / "probe_level_activity.tsv", probe_level_rows, probe_fields)
    write_svg(summary_rows)
    write_paired_gene_summary(summary_rows)
    matplotlib_files = write_matplotlib_figures(summary_rows, plotting_stack)
    figure_files = ["probe_set_activity_heatmap.svg", *matplotlib_files]
    meta = {
        "schema": SCHEMA,
        "renderer": provenance_for(Path(__file__).resolve(), "renderer source"),
        "render_runtime": render_runtime,
        "inputs": input_provenance,
        "annotation_resolution": annotation_usage,
        "genes": GENES,
        "samples": SAMPLES,
        "probe_sets": len(summary_rows),
        "pm_probes": len(probe_level_rows),
        "groups": GROUPS,
        "sample_pairings": [dict(pairing) for pairing in PAIRINGS],
        "contrast_definitions": CONTRAST_DEFINITIONS,
        "interpretation_limitations": INTERPRETATION_LIMITATIONS,
        "input_limits": {
            "max_genes": MAX_GENES,
            "max_probe_sets": MAX_PROBE_SETS,
            "max_pm_probes": MAX_PM_PROBES,
            "max_raw_feature_rows_scanned": MAX_RAW_FEATURE_ROWS,
            "max_annotation_lines": MAX_ANNOTATION_LINES,
            "max_input_line_characters": MAX_INPUT_LINE_CHARS,
            "max_vendor_member_bytes": MAX_VENDOR_MEMBER_BYTES,
        },
        "figure_files": figure_files,
    }
    (OUT / "manifest.json").write_text(
        json.dumps(meta, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    write_html(summary_rows, meta)
    write_readme(meta)
    print(f"wrote {OUT}")
    print(json.dumps(meta, sort_keys=True))


def main() -> None:
    try:
        run(parse_args())
    except (
        csv.Error,
        OSError,
        RuntimeError,
        UnicodeError,
        ValueError,
        zipfile.BadZipFile,
    ) as exc:
        raise SystemExit(f"error: {exc}") from None


if __name__ == "__main__":
    main()
