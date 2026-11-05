#!/usr/bin/env python3
"""Project saved GENtle SAM observations to transcript-oriented UTR reach.

This is a downstream projection over a completed
``gentle.observed_sam_structures_data.v1`` artifact.  It does not reopen SAM,
invent transcript assignments, or infer CDS boundaries.  The caller supplies a
frozen, digest-bound transcript model with genomic exon and CDS intervals.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

from compare_rna_read_end_metrics import staged_output


DATA_SCHEMA = "gentle.observed_sam_structures_data.v1"
MODEL_SCHEMA = "gentle.rna_read_end_transcript_models.v1"
REPORT_SCHEMA = "gentle.rna_read_end_projection.v1"
TOKEN = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]{0,127}$")
HASH = re.compile(r"^sha256:[0-9a-f]{64}$")
ACCEPTED_RELATIONS = {"exact_full_chain", "exact_partial_chain", "within_exon"}
MAX_MODELS = 100_000
MAX_INTERVALS = 2_000_000
MAX_UNITS = 2_000_000


def fail(message: str) -> ValueError:
    return ValueError(f"rna_read_end_projection: {message}")


def load_json(path: Path) -> tuple[dict, bytes]:
    raw = path.read_bytes()
    try:
        value = json.loads(raw)
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise fail(f"invalid UTF-8 JSON: {path}") from exc
    if not isinstance(value, dict):
        raise fail(f"JSON root must be an object: {path}")
    return value, raw


def interval(value: object, field: str) -> tuple[int, int]:
    if (
        not isinstance(value, list)
        or len(value) != 2
        or any(isinstance(item, bool) or not isinstance(item, int) for item in value)
        or value[0] < 0
        or value[0] >= value[1]
    ):
        raise fail(f"{field} must be a nonempty zero-based half-open interval")
    return value[0], value[1]


def merged(values: Iterable[tuple[int, int]]) -> list[tuple[int, int]]:
    out: list[tuple[int, int]] = []
    for start, end in sorted(values):
        if out and start <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], end))
        else:
            out.append((start, end))
    return out


@dataclass(frozen=True)
class Model:
    transcript_id: str
    gene_id: str
    contig: str
    strand: str
    exons: tuple[tuple[int, int], ...]
    cds: tuple[tuple[int, int], ...]
    transcript_length: int
    cds_start: int
    cds_end: int
    offsets: tuple[tuple[int, int, int], ...]

    def project(self, ranges: Iterable[tuple[int, int]]) -> list[tuple[int, int]] | None:
        projected: list[tuple[int, int]] = []
        for start, end in ranges:
            covered = 0
            for exon_start, exon_end, offset in self.offsets:
                left, right = max(start, exon_start), min(end, exon_end)
                if left >= right:
                    continue
                covered += right - left
                if self.strand == "+":
                    projected.append((offset + left - exon_start, offset + right - exon_start))
                else:
                    projected.append((offset + exon_end - right, offset + exon_end - left))
            if covered != end - start:
                return None
        return merged(projected)


def prepare_model(row: object, index: int) -> Model:
    if not isinstance(row, dict) or set(row) != {
        "transcript_id",
        "gene_id",
        "contig",
        "strand",
        "genomic_exons_0based_half_open",
        "genomic_cds_0based_half_open",
    }:
        raise fail(f"transcripts[{index}] has unexpected fields")
    transcript_id, gene_id, contig, strand = (
        row["transcript_id"], row["gene_id"], row["contig"], row["strand"]
    )
    if (
        not all(isinstance(value, str) and TOKEN.fullmatch(value) for value in (transcript_id, gene_id, contig))
        or strand not in {"+", "-"}
    ):
        raise fail(f"transcripts[{index}] has invalid identity or strand")
    raw_exons, raw_cds = row["genomic_exons_0based_half_open"], row["genomic_cds_0based_half_open"]
    if not isinstance(raw_exons, list) or not isinstance(raw_cds, list) or not raw_exons or not raw_cds:
        raise fail(f"transcripts[{index}] needs nonempty exon and CDS intervals")
    exons = [interval(value, f"transcripts[{index}].exons") for value in raw_exons]
    cds = [interval(value, f"transcripts[{index}].cds") for value in raw_cds]
    if exons != sorted(exons) or cds != sorted(cds) or merged(exons) != exons or merged(cds) != cds:
        raise fail(f"transcripts[{index}] intervals must be sorted and nonoverlapping")
    if any(sum(max(0, min(z, ez) - max(a, ea)) for ea, ez in exons) != z - a for a, z in cds):
        raise fail(f"transcripts[{index}] CDS must be contained in exons")
    oriented = exons if strand == "+" else list(reversed(exons))
    offsets: list[tuple[int, int, int]] = []
    cursor = 0
    for start, end in oriented:
        offsets.append((start, end, cursor))
        cursor += end - start
    provisional = Model(
        transcript_id,
        gene_id,
        contig,
        strand,
        tuple(exons),
        tuple(cds),
        cursor,
        0,
        0,
        tuple(offsets),
    )
    projected_cds = provisional.project(cds)
    if not projected_cds:
        raise fail(f"transcripts[{index}] CDS projection is empty")
    return Model(
        transcript_id,
        gene_id,
        contig,
        strand,
        tuple(exons),
        tuple(cds),
        cursor,
        min(start for start, _ in projected_cds),
        max(end for _, end in projected_cds),
        tuple(offsets),
    )


def read_models(document: dict, analysis: dict) -> tuple[str, dict[str, Model]]:
    expected = {
        "schema",
        "source_id",
        "assembly",
        "reference_fasta_sha256",
        "annotation_sha256",
        "coordinate_convention",
        "model_source",
        "transcripts",
    }
    if set(document) != expected or document.get("schema") != MODEL_SCHEMA:
        raise fail("transcript-model document has an unsupported shape or schema")
    source_id = document["source_id"]
    if (
        not isinstance(source_id, str)
        or not TOKEN.fullmatch(source_id)
        or document["coordinate_convention"] != "zero_based_half_open"
        or not isinstance(document["assembly"], str)
        or not document["assembly"]
        or not HASH.fullmatch(str(document["reference_fasta_sha256"]))
        or not HASH.fullmatch(str(document["annotation_sha256"]))
        or not isinstance(document["model_source"], dict)
        or set(document["model_source"]) != {"format", "name", "sha256"}
        or document["model_source"].get("format") != "gtf"
        or not isinstance(document["model_source"].get("name"), str)
        or not document["model_source"]["name"]
        or not HASH.fullmatch(str(document["model_source"].get("sha256")))
        or not isinstance(document["transcripts"], list)
        or not 1 <= len(document["transcripts"]) <= MAX_MODELS
    ):
        raise fail("invalid transcript-model binding or limits")
    reference = analysis.get("reference")
    bindings = analysis.get("annotation_bindings")
    if not isinstance(reference, dict) or not isinstance(bindings, list):
        raise fail("saved analysis lacks reference or annotation bindings")
    if reference.get("assembly") != document["assembly"] or reference.get("fasta_sha256") != document["reference_fasta_sha256"]:
        raise fail("transcript model does not match the saved reference binding")
    matched = [row for row in bindings if isinstance(row, dict) and row.get("source_id") == source_id]
    if len(matched) != 1 or matched[0].get("annotation_sha256") != document["annotation_sha256"]:
        raise fail("transcript model does not match exactly one saved annotation binding")
    models: dict[str, Model] = {}
    total_intervals = 0
    for index, row in enumerate(document["transcripts"]):
        model = prepare_model(row, index)
        total_intervals += len(model.exons) + len(model.cds)
        if total_intervals > MAX_INTERVALS or model.transcript_id in models:
            raise fail("transcript-model interval budget or unique-ID constraint failed")
        models[model.transcript_id] = model
    return source_id, models


def comparison_models(observation: dict, source_id: str) -> set[str]:
    comparisons = observation.get("comparisons")
    if not isinstance(comparisons, list):
        return set()
    rows = [row for row in comparisons if isinstance(row, dict) and row.get("source_id") == source_id]
    if len(rows) != 1 or rows[0].get("availability") != "available":
        return set()
    comparison = rows[0].get("comparison")
    if not isinstance(comparison, dict) or comparison.get("resolution") != "resolved":
        return set()
    models = comparison.get("models")
    if not isinstance(models, list):
        return set()
    return {
        row["transcript_id"]
        for row in models
        if isinstance(row, dict)
        and isinstance(row.get("transcript_id"), str)
        and row.get("relation") in ACCEPTED_RELATIONS
    }


def reach(intervals: list[tuple[int, int]], boundary: int, side: str) -> int | None:
    for start, end in intervals:
        if start <= boundary <= end:
            return boundary - start if side == "five" else end - boundary
    return None


def one_projection(path: list[dict], model: Model) -> tuple[int | None, str, int | None, str] | None:
    ranges: list[tuple[int, int]] = []
    query_ranges: list[tuple[int, int]] = []
    axis_lengths: set[int] = set()
    orientations: set[bool] = set()
    for observation in path:
        if observation.get("contig") != model.contig:
            return None
        blocks = observation.get("blocks")
        if not isinstance(blocks, list) or not blocks:
            return None
        for block in blocks:
            if not isinstance(block, dict):
                return None
            target = interval(block.get("target"), "observation block target")
            query = interval(block.get("query"), "observation block query")
            if target[1] - target[0] != query[1] - query[0]:
                return None
            ranges.append(target)
            query_ranges.append(query)
        axis = observation.get("query_axis_length")
        strand = observation.get("query_to_reference_strand")
        if isinstance(axis, bool) or not isinstance(axis, int) or axis <= 0 or strand not in {"+", "-"}:
            return None
        axis_lengths.add(axis)
        orientations.add(strand == model.strand)
    if len(axis_lengths) != 1 or len(orientations) != 1:
        return None
    projected = model.project(ranges)
    if not projected:
        return None
    five = reach(projected, model.cds_start, "five")
    three = reach(projected, model.cds_end, "three")
    qmin = min(start for start, _ in query_ranges)
    qmax = max(end for _, end in query_ranges)
    axis = next(iter(axis_lengths))
    same = next(iter(orientations))
    five_clipped = qmin > 0 if same else qmax < axis
    three_clipped = qmax < axis if same else qmin > 0
    return (
        five,
        "right_censored" if five is not None and five_clipped else "observed" if five is not None else "missing",
        three,
        "right_censored" if three is not None and three_clipped else "observed" if three is not None else "missing",
    )


def project_unit(unit: dict, observations: dict[str, dict], source_id: str, models: dict[str, Model]):
    hypotheses = unit.get("hypotheses")
    if not isinstance(hypotheses, list) or not hypotheses:
        return None, "no_optimal_hypothesis", []
    values = []
    candidate_ids: set[str] = set()
    for hypothesis in hypotheses:
        if not isinstance(hypothesis, list) or not hypothesis:
            return None, "invalid_optimal_hypothesis", []
        try:
            path = [observations[identifier] for identifier in hypothesis]
        except (KeyError, TypeError):
            return None, "unbound_optimal_hypothesis", []
        candidates: set[str] | None = None
        for observation in path:
            current = comparison_models(observation, source_id) & models.keys()
            candidates = current if candidates is None else candidates & current
        if not candidates:
            return None, "no_compatible_cds_model", []
        for transcript_id in sorted(candidates):
            candidate_ids.add(transcript_id)
            value = one_projection(path, models[transcript_id])
            if value is None:
                return None, "model_projection_incompatible", sorted(candidate_ids)
            values.append(value)
    if not values:
        return None, "no_projected_values", sorted(candidate_ids)
    five = {(value[0], value[1]) for value in values}
    three = {(value[2], value[3]) for value in values}
    result = (
        next(iter(five)) if len(five) == 1 else (None, "missing"),
        next(iter(three)) if len(three) == 1 else (None, "missing"),
    )
    reasons = (
        "resolved" if len(five) == 1 else "ambiguous_transcript_end_geometry",
        "resolved" if len(three) == 1 else "ambiguous_transcript_end_geometry",
    )
    return (result, reasons), "resolved", sorted(candidate_ids)


def write_rows(path: Path, header: list[str], rows: list[list[object]]) -> None:
    with staged_output(path) as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def run(data_path: Path, model_path: Path, sample_id: str, output_dir: Path) -> None:
    if not TOKEN.fullmatch(sample_id):
        raise fail("sample_id must be a stable token")
    data, data_raw = load_json(data_path)
    model_document, model_raw = load_json(model_path)
    if data.get("schema") != DATA_SCHEMA or not isinstance(data.get("analysis"), dict):
        raise fail("input must be a completed observed-SAM data artifact")
    analysis = data["analysis"]
    source_id, models = read_models(model_document, analysis)
    units = analysis.get("units")
    observations_raw = analysis.get("observations")
    if not isinstance(units, list) or not isinstance(observations_raw, list) or len(units) > MAX_UNITS:
        raise fail("saved analysis unit or observation limits failed")
    observations: dict[str, dict] = {}
    for observation in observations_raw:
        identifier = observation.get("observation_id") if isinstance(observation, dict) else None
        if not isinstance(identifier, str) or identifier in observations:
            raise fail("observation IDs must be unique strings")
        observations[identifier] = observation
    metrics: list[list[object]] = []
    details: list[list[object]] = []
    counts = {"resolved": 0, "missing": 0, "right_censored": 0}
    for unit in units:
        key = unit.get("unit") if isinstance(unit, dict) else None
        if not isinstance(key, dict) or not isinstance(key.get("unit_id"), str) or not isinstance(key.get("source_sha256"), str):
            raise fail("unit identity is missing")
        unit_id, source_sha256 = key["unit_id"], key["source_sha256"]
        projection, reason, transcript_ids = project_unit(unit, observations, source_id, models)
        gene_ids = sorted({models[transcript_id].gene_id for transcript_id in transcript_ids})
        axis = unit.get("unclipped_axis_length")
        if axis is None:
            metrics.append([sample_id, "read_length_bp", "", "missing"])
            details.append([sample_id, source_sha256, unit_id, "read_length_bp", "", "missing", "axis_length_unavailable", ",".join(transcript_ids), ",".join(gene_ids)])
            counts["missing"] += 1
        elif isinstance(axis, bool) or not isinstance(axis, int) or axis <= 0:
            raise fail("unclipped_axis_length must be positive or null")
        else:
            metrics.append([sample_id, "read_length_bp", axis, "observed"])
            details.append([sample_id, source_sha256, unit_id, "read_length_bp", axis, "observed", "resolved", ",".join(transcript_ids), ",".join(gene_ids)])
            counts["resolved"] += 1
        for metric, index in (("utr5_reach_bp", 0), ("utr3_reach_bp", 1)):
            value = None if projection is None else projection[0][index][0]
            status = "missing" if projection is None else projection[0][index][1]
            metric_reason = reason if projection is None else projection[1][index]
            metrics.append([sample_id, metric, "" if value is None else value, status])
            details.append([
                sample_id,
                source_sha256,
                unit_id,
                metric,
                "" if value is None else value,
                status,
                metric_reason,
                ",".join(transcript_ids),
                ",".join(gene_ids),
            ])
            counts["missing" if status == "missing" else "right_censored" if status == "right_censored" else "resolved"] += 1
    if not output_dir.name or ".." in output_dir.parts:
        raise fail("output must name a fresh directory without '..' components")
    parent = output_dir.parent.resolve(strict=True)
    output_dir = parent / output_dir.name
    output_dir.mkdir(mode=0o700)
    write_rows(output_dir / "metrics.tsv", ["sample_id", "metric", "value_bp", "observation_status"], metrics)
    write_rows(
        output_dir / "projections.tsv",
        ["sample_id", "source_sha256", "unit_id", "metric", "value_bp", "observation_status", "reason", "transcript_ids", "gene_ids"],
        details,
    )
    report = {
        "schema": REPORT_SCHEMA,
        "sample_id": sample_id,
        "source": {"name": data_path.name, "sha256": hashlib.sha256(data_raw).hexdigest()},
        "transcript_models": {
            "name": model_path.name,
            "sha256": hashlib.sha256(model_raw).hexdigest(),
            "source_id": source_id,
            "count": len(models),
            "model_source": model_document["model_source"],
        },
        "counts": counts,
        "unit_count": len(units),
        "metric_rows": len(metrics),
        "policy": {
            "optimal_paths_only": True,
            "accepted_model_relations": sorted(ACCEPTED_RELATIONS),
            "transcript_ambiguity": "retain only when every candidate and optimal path has identical metric/status",
            "utr_anchor": "contiguous projected alignment must cross the annotated CDS boundary",
            "clipped_end": "right_censored",
            "unavailable_or_unresolved": "missing",
        },
        "nonclaims": [
            "UTR reach is annotation-relative alignment evidence, not a molecule-complete or polyadenylation-site claim.",
            "A 5-prime UTR comparison is a technical/control comparison, not a biological negative control.",
            "Library preparation, truncation, mapping, transcript choice and annotation remain possible explanations.",
        ],
    }
    with staged_output(output_dir / "report.json") as handle:
        json.dump(report, handle, indent=2, sort_keys=True)
        handle.write("\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--observed-sam-data", required=True, type=Path)
    parser.add_argument("--transcript-models", required=True, type=Path)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    try:
        run(args.observed_sam_data, args.transcript_models, args.sample_id, args.output_dir)
    except (OSError, ValueError) as exc:
        parser.error(str(exc))


if __name__ == "__main__":
    main()
