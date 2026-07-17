#!/usr/bin/env python3
"""Select reproducible CUT&RUN regions for motif-discovery pilot analyses.

All coordinates are 0-based half-open throughout: candidate intervals, BigWig
queries, merged islands, centered windows, FASTA slices, and FASTA headers.
The fixed pipeline is chromosome/score filtering, merging, exact BigWig
scoring, ranking, centering, and FASTA extraction. BigWig scoring always uses
``pyBigWig.stats(..., exact=True)``; regions without coverage score 0.0.

The default pilot chromosome policy keeps normalized chromosomes 2 through 22,
X, and Y. It deliberately excludes chromosome 1, mitochondria, and alternative
or unplaced contigs. Names are normalized only for policy matching: an optional
``chr`` prefix is removed, names are upper-cased, and M/MT/chrM all normalize to
MT. Concrete names must still agree across candidates, BigWig, and FASTA index.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
import shutil
import subprocess
import sys
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable, Iterable, Iterator, Mapping, Sequence, TextIO


DEFAULT_INCLUDE_CHROMS = tuple([str(value) for value in range(2, 23)] + ["X", "Y"])
Scorer = Callable[[str, int, int], tuple[float, float]]


class RegionSelectionError(RuntimeError):
    """A user-facing input, dependency, or consistency error."""


@dataclass(frozen=True)
class Segment:
    chrom: str
    start: int
    end: int
    score: float

    @property
    def length(self) -> int:
        return self.end - self.start


@dataclass(frozen=True)
class MergedIsland:
    chrom: str
    start: int
    end: int
    name: str
    support: float

    @property
    def length(self) -> int:
        return self.end - self.start


@dataclass(frozen=True)
class ScoredIsland:
    chrom: str
    start: int
    end: int
    name: str
    bw_max: float
    bw_mean: float
    support: float
    length: int
    rank: int


@dataclass(frozen=True)
class CenteredWindow:
    chrom: str
    start: int
    end: int
    name: str
    rank: int


@dataclass(frozen=True)
class FaiRecord:
    name: str
    length: int
    offset: int
    line_bases: int
    line_width: int


def normalize_chromosome(name: str) -> str:
    """Normalize a chromosome name for allowlist comparisons only."""

    normalized = name.strip()
    if normalized[:3].lower() == "chr":
        normalized = normalized[3:]
    normalized = normalized.upper()
    return "MT" if normalized in {"M", "MT"} else normalized


def parse_chromosome_set(value: str, *, allow_canonical: bool) -> set[str]:
    """Parse a comma-separated normalized chromosome set."""

    if allow_canonical and value.strip().lower() == "canonical":
        return set(DEFAULT_INCLUDE_CHROMS)
    if not value.strip():
        return set()
    parsed = {normalize_chromosome(item) for item in value.split(",") if item.strip()}
    if not parsed:
        raise RegionSelectionError("chromosome list must contain at least one name")
    return parsed


def resolve_chromosome_policy(include: str, exclude: str) -> tuple[set[str], set[str]]:
    included = parse_chromosome_set(include, allow_canonical=True)
    excluded = parse_chromosome_set(exclude, allow_canonical=False)
    resolved = included - excluded
    if not resolved:
        raise RegionSelectionError("chromosome policy excludes every included chromosome")
    return resolved, excluded


def _open_text(path: Path) -> TextIO:
    if path.name.lower().endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def _candidate_kind(path: Path) -> str:
    name = path.name.lower()
    if name.endswith(".gz"):
        name = name[:-3]
    return "bedgraph" if name.endswith((".bedgraph", ".bdg")) else "bed"


def _resolve_score_column(path: Path, requested: int | None, field_count: int) -> int:
    if requested is not None:
        if requested < 1:
            raise RegionSelectionError("--score-column is 1-based and must be at least 1")
        return requested
    if _candidate_kind(path) == "bedgraph":
        return 4
    return 5 if field_count >= 5 else 4


def read_candidate_segments(
    path: Path, score_column: int | None = None
) -> tuple[list[Segment], int]:
    """Read BED/bedGraph candidates and return segments plus 1-based score column."""

    if not path.is_file():
        raise RegionSelectionError(f"candidate file does not exist: {path}")

    rows: list[tuple[int, list[str]]] = []
    with _open_text(path) as handle:
        for line_number, line in enumerate(handle, start=1):
            stripped = line.strip()
            if not stripped or stripped.startswith(("#", "track ", "browser ")):
                continue
            fields = stripped.split()
            if len(fields) < 3:
                raise RegionSelectionError(
                    f"{path}:{line_number}: expected at least 3 BED columns"
                )
            rows.append((line_number, fields))

    if not rows:
        raise RegionSelectionError(f"candidate file has no data rows: {path}")
    resolved_column = _resolve_score_column(path, score_column, len(rows[0][1]))

    segments: list[Segment] = []
    for line_number, fields in rows:
        if len(fields) < resolved_column:
            raise RegionSelectionError(
                f"{path}:{line_number}: score column {resolved_column} is missing"
            )
        try:
            start = int(fields[1])
            end = int(fields[2])
            score = float(fields[resolved_column - 1])
        except ValueError as exc:
            raise RegionSelectionError(
                f"{path}:{line_number}: invalid coordinate or score: {exc}"
            ) from exc
        if start < 0 or end <= start:
            raise RegionSelectionError(
                f"{path}:{line_number}: expected a non-empty 0-based half-open interval"
            )
        if not math.isfinite(score):
            raise RegionSelectionError(f"{path}:{line_number}: score must be finite")
        segments.append(Segment(fields[0], start, end, score))
    return segments, resolved_column


def filter_segments_by_chromosome(
    segments: Iterable[Segment], included: set[str]
) -> list[Segment]:
    return [
        segment
        for segment in segments
        if normalize_chromosome(segment.chrom) in included
    ]


def nearest_rank_percentile(values: Sequence[float], percentile: float) -> float:
    """Return the nearest-rank percentile: sorted[ceil(P/100*N)-1]."""

    if not values:
        raise RegionSelectionError("cannot calculate a percentile over no segments")
    if not 0.0 <= percentile <= 100.0:
        raise RegionSelectionError("--score-percentile must be between 0 and 100")
    ordered = sorted(values)
    rank = max(1, math.ceil(percentile / 100.0 * len(ordered)))
    return ordered[rank - 1]


def filter_segments_by_score(
    segments: Sequence[Segment],
    *,
    min_score: float | None,
    score_percentile: float | None,
) -> tuple[list[Segment], float | None]:
    if min_score is not None:
        if not math.isfinite(min_score):
            raise RegionSelectionError("--min-score must be finite")
        threshold = min_score
    elif score_percentile is not None:
        threshold = nearest_rank_percentile(
            [segment.score for segment in segments], score_percentile
        )
    else:
        threshold = None
    if threshold is None:
        return list(segments), None
    return [segment for segment in segments if segment.score >= threshold], threshold


def merge_segments(
    segments: Iterable[Segment], merge_distance: int
) -> list[MergedIsland]:
    """Merge sorted-compatible segments and retain score-by-length support."""

    if merge_distance < 0:
        raise RegionSelectionError("--merge-distance must be non-negative")
    ordered = sorted(segments, key=lambda item: (item.chrom, item.start, item.end, item.score))
    merged_rows: list[tuple[str, int, int, float]] = []
    current: tuple[str, int, int, float] | None = None
    for segment in ordered:
        weighted_support = segment.score * segment.length
        if current is None:
            current = (segment.chrom, segment.start, segment.end, weighted_support)
            continue
        chrom, start, end, support = current
        if segment.chrom == chrom and segment.start - end <= merge_distance:
            current = (
                chrom,
                start,
                max(end, segment.end),
                support + weighted_support,
            )
        else:
            merged_rows.append(current)
            current = (segment.chrom, segment.start, segment.end, weighted_support)
    if current is not None:
        merged_rows.append(current)

    return [
        MergedIsland(chrom, start, end, f"island_{index:06d}", support)
        for index, (chrom, start, end, support) in enumerate(merged_rows, start=1)
    ]


def score_and_rank_islands(
    islands: Iterable[MergedIsland], scorer: Scorer
) -> list[ScoredIsland]:
    """Score and apply the documented deterministic total-order ranking."""

    scored: list[ScoredIsland] = []
    for island in islands:
        bw_max, bw_mean = scorer(island.chrom, island.start, island.end)
        if not math.isfinite(bw_max) or not math.isfinite(bw_mean):
            raise RegionSelectionError(
                f"BigWig returned a non-finite score for {island.chrom}:{island.start}-{island.end}"
            )
        scored.append(
            ScoredIsland(
                chrom=island.chrom,
                start=island.start,
                end=island.end,
                name=island.name,
                bw_max=float(bw_max),
                bw_mean=float(bw_mean),
                support=island.support,
                length=island.length,
                rank=0,
            )
        )

    scored.sort(
        key=lambda item: (
            -item.bw_max,
            -item.bw_mean,
            -item.support,
            -item.length,
            item.chrom,
            item.start,
            item.end,
        )
    )
    return [
        ScoredIsland(**{**asdict(item), "rank": rank})
        for rank, item in enumerate(scored, start=1)
    ]


def select_top_islands(
    ranked: Sequence[ScoredIsland], top_n: int | None, top_fraction: float | None
) -> list[ScoredIsland]:
    if top_n is not None:
        if top_n < 1:
            raise RegionSelectionError("--top-n must be at least 1")
        count = min(top_n, len(ranked))
    else:
        if top_fraction is None or not 0.0 < top_fraction <= 1.0:
            raise RegionSelectionError("--top-fraction must be greater than 0 and at most 1")
        count = min(len(ranked), max(1, math.ceil(len(ranked) * top_fraction))) if ranked else 0
    return list(ranked[:count])


def centered_windows(
    islands: Iterable[ScoredIsland], window_size: int, chromosome_lengths: Mapping[str, int]
) -> list[CenteredWindow]:
    if window_size < 1:
        raise RegionSelectionError("--window must be at least 1")
    half = window_size // 2
    windows: list[CenteredWindow] = []
    for island in islands:
        if island.chrom not in chromosome_lengths:
            raise RegionSelectionError(
                f"chromosome {island.chrom!r} is absent from the FASTA index"
            )
        center = (island.start + island.end) // 2
        raw_start = center - half
        raw_end = raw_start + window_size
        chrom_length = chromosome_lengths[island.chrom]
        start = max(0, min(raw_start, chrom_length))
        end = max(start, min(raw_end, chrom_length))
        windows.append(CenteredWindow(island.chrom, start, end, island.name, island.rank))
    return windows


def read_fai(path: Path) -> dict[str, FaiRecord]:
    if not path.is_file():
        raise RegionSelectionError(f"FASTA index does not exist: {path}")
    records: dict[str, FaiRecord] = {}
    with path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                raise RegionSelectionError(f"{path}:{line_number}: malformed .fai row")
            try:
                record = FaiRecord(
                    name=fields[0],
                    length=int(fields[1]),
                    offset=int(fields[2]),
                    line_bases=int(fields[3]),
                    line_width=int(fields[4]),
                )
            except ValueError as exc:
                raise RegionSelectionError(f"{path}:{line_number}: malformed .fai integer") from exc
            if (
                record.length < 0
                or record.offset < 0
                or record.line_bases < 1
                or record.line_width < record.line_bases
            ):
                raise RegionSelectionError(f"{path}:{line_number}: invalid .fai geometry")
            records[record.name] = record
    return records


def read_fasta_slice(handle, record: FaiRecord, start: int, end: int) -> str:
    """Read one 0-based half-open FASTA slice using .fai byte geometry."""

    if start < 0 or end < start or end > record.length:
        raise RegionSelectionError(
            f"FASTA slice {record.name}:{start}-{end} is outside [0,{record.length})"
        )
    chunks: list[bytes] = []
    position = start
    while position < end:
        line_index, line_offset = divmod(position, record.line_bases)
        chunk_length = min(end - position, record.line_bases - line_offset)
        byte_offset = record.offset + line_index * record.line_width + line_offset
        handle.seek(byte_offset)
        chunk = handle.read(chunk_length)
        if len(chunk) != chunk_length:
            raise RegionSelectionError(
                f"FASTA ended while reading {record.name}:{start}-{end}"
            )
        chunks.append(chunk)
        position += chunk_length
    try:
        return b"".join(chunks).decode("ascii")
    except UnicodeDecodeError as exc:
        raise RegionSelectionError(f"FASTA sequence for {record.name} is not ASCII") from exc


def write_internal_fasta(
    fasta_path: Path,
    output_path: Path,
    windows: Sequence[CenteredWindow],
    fai_records: Mapping[str, FaiRecord],
) -> None:
    with fasta_path.open("rb") as fasta, output_path.open("w", encoding="ascii") as output:
        for window in windows:
            sequence = read_fasta_slice(
                fasta, fai_records[window.chrom], window.start, window.end
            )
            output.write(f">{window.chrom}:{window.start}-{window.end}\n{sequence}\n")


def run_bedtools_getfasta(
    bedtools: str, fasta_path: Path, bed_path: Path, output_path: Path
) -> str:
    try:
        version = subprocess.run(
            [bedtools, "--version"],
            check=True,
            capture_output=True,
            text=True,
        ).stdout.strip()
        subprocess.run(
            [bedtools, "getfasta", "-fi", str(fasta_path), "-bed", str(bed_path), "-fo", str(output_path)],
            check=True,
        )
    except (OSError, subprocess.CalledProcessError) as exc:
        raise RegionSelectionError(f"bedtools getfasta failed: {exc}") from exc
    return version


class PyBigWigScorer:
    """Lazy pyBigWig adapter enforcing exact statistics."""

    def __init__(self, path: Path) -> None:
        try:
            import pyBigWig  # type: ignore[import-not-found]
        except ImportError as exc:
            raise RegionSelectionError(
                "BigWig scoring requires pyBigWig; install it with `pip install pyBigWig`"
            ) from exc
        self.version = str(getattr(pyBigWig, "__version__", "unknown"))
        try:
            self._handle = pyBigWig.open(str(path))
        except Exception as exc:
            raise RegionSelectionError(f"could not open BigWig {path}: {exc}") from exc
        if self._handle is None:
            raise RegionSelectionError(f"could not open BigWig {path}")
        self.chromosomes = dict(self._handle.chroms())

    def close(self) -> None:
        self._handle.close()

    def __enter__(self) -> "PyBigWigScorer":
        return self

    def __exit__(self, _exc_type, _exc, _traceback) -> None:
        self.close()

    def __call__(self, chrom: str, start: int, end: int) -> tuple[float, float]:
        return (
            self._stat(chrom, start, end, "max"),
            self._stat(chrom, start, end, "mean"),
        )

    def _stat(self, chrom: str, start: int, end: int, statistic: str) -> float:
        try:
            values = self._handle.stats(chrom, start, end, type=statistic, exact=True)
        except Exception as exc:
            raise RegionSelectionError(
                f"BigWig exact {statistic} failed for {chrom}:{start}-{end}: {exc}"
            ) from exc
        value = values[0] if values else None
        return 0.0 if value is None else float(value)


def guard_assembly_chromosomes(
    segments: Sequence[Segment],
    bigwig_chromosomes: Mapping[str, int],
    fai_records: Mapping[str, FaiRecord],
    *,
    allow_missing: bool,
) -> tuple[list[Segment], list[dict[str, object]]]:
    """Verify concrete candidate names exist in both sequence sources."""

    missing_rows: list[dict[str, object]] = []
    for chrom in sorted({segment.chrom for segment in segments}):
        missing_from = []
        if chrom not in bigwig_chromosomes:
            missing_from.append("BigWig header")
        if chrom not in fai_records:
            missing_from.append("FASTA .fai")
        if missing_from:
            missing_rows.append({"chromosome": chrom, "missing_from": missing_from})

    if missing_rows and not allow_missing:
        details = "; ".join(
            f"{row['chromosome']} missing from {', '.join(row['missing_from'])}"
            for row in missing_rows
        )
        raise RegionSelectionError(
            "assembly/chromosome-name mismatch: "
            f"{details}. Candidate chromosomes must exist under the same concrete name "
            "in both the BigWig header and FASTA .fai; use --allow-missing-chroms "
            "only for exploratory runs."
        )

    missing_names = {str(row["chromosome"]) for row in missing_rows}
    if missing_rows:
        print(
            "warning: --allow-missing-chroms skips unusable candidate chromosomes: "
            + "; ".join(
                f"{row['chromosome']} ({', '.join(row['missing_from'])})"
                for row in missing_rows
            ),
            file=sys.stderr,
        )
    return [segment for segment in segments if segment.chrom not in missing_names], missing_rows


def _format_number(value: float) -> str:
    return format(value, ".12g")


def _write_filtered(path: Path, segments: Sequence[Segment]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for index, segment in enumerate(
            sorted(segments, key=lambda item: (item.chrom, item.start, item.end, item.score)),
            start=1,
        ):
            handle.write(
                f"{segment.chrom}\t{segment.start}\t{segment.end}\t"
                f"segment_{index:06d}\t{_format_number(segment.score)}\n"
            )


def _write_merged(path: Path, islands: Sequence[MergedIsland]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for island in islands:
            handle.write(
                f"{island.chrom}\t{island.start}\t{island.end}\t{island.name}\t"
                f"{_format_number(island.support)}\n"
            )


def _write_ranked(path: Path, islands: Sequence[ScoredIsland]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for island in islands:
            handle.write(
                "\t".join(
                    [
                        island.chrom,
                        str(island.start),
                        str(island.end),
                        island.name,
                        _format_number(island.bw_max),
                        _format_number(island.bw_mean),
                        _format_number(island.support),
                        str(island.length),
                        str(island.rank),
                    ]
                )
                + "\n"
            )


def _write_windows(path: Path, windows: Sequence[CenteredWindow]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for window in windows:
            handle.write(
                f"{window.chrom}\t{window.start}\t{window.end}\t{window.name}\n"
            )


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def fingerprint(path: Path) -> dict[str, object]:
    resolved = path.resolve()
    return {
        "path": str(resolved),
        "size_bytes": resolved.stat().st_size,
        "sha256": sha256_file(resolved),
    }


def git_commit() -> str | None:
    repo_root = Path(__file__).resolve().parents[1]
    try:
        return subprocess.run(
            ["git", "-C", str(repo_root), "rev-parse", "HEAD"],
            check=True,
            capture_output=True,
            text=True,
        ).stdout.strip() or None
    except (OSError, subprocess.CalledProcessError):
        return None


def _flatten_provenance(prefix: str, value: object) -> Iterator[tuple[str, str]]:
    if isinstance(value, dict):
        for key in sorted(value):
            child = f"{prefix}.{key}" if prefix else str(key)
            yield from _flatten_provenance(child, value[key])
    elif isinstance(value, list):
        yield prefix, json.dumps(value, sort_keys=True, separators=(",", ":"))
    elif value is None:
        yield prefix, ""
    elif isinstance(value, bool):
        yield prefix, "true" if value else "false"
    else:
        yield prefix, str(value)


def write_provenance(json_path: Path, provenance: Mapping[str, object]) -> Path:
    json_path.parent.mkdir(parents=True, exist_ok=True)
    json_path.write_text(
        json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    tsv_path = json_path.with_suffix(".tsv")
    with tsv_path.open("w", encoding="utf-8") as handle:
        handle.write("key\tvalue\n")
        for key, value in _flatten_provenance("", provenance):
            handle.write(f"{key}\t{value}\n")
    return tsv_path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Select exact-scored CUT&RUN regions for motif-discovery pilots.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Coordinates are 0-based half-open everywhere. Chromosome names are "
            "matched to a normalized allowlist (optional chr prefix removed; upper-case; "
            "M/MT/chrM are mitochondrial). The canonical pilot policy keeps 2-22,X,Y "
            "and excludes chromosome 1, MT, and alt/unplaced contigs. BigWig max and "
            "mean use pyBigWig stats(..., exact=True); missing coverage becomes 0.0."
        ),
    )
    parser.add_argument("--bigwig", required=True, type=Path, help="BigWig signal used for exact max/mean scoring.")
    parser.add_argument(
        "--candidates",
        type=Path,
        help="Candidate BED/bedGraph. Required for this pilot; BigWig interval discovery is out of scope.",
    )
    parser.add_argument("--genome-fasta", required=True, type=Path, help="Assembly-matched indexed genome FASTA.")
    parser.add_argument("--fai", type=Path, help="FASTA .fai; defaults to <genome-fasta>.fai.")
    parser.add_argument("--output-dir", required=True, type=Path, help="Directory for all selected-region artifacts.")
    selection = parser.add_mutually_exclusive_group(required=True)
    selection.add_argument("--top-n", type=int, help="Retain at most this many ranked islands.")
    selection.add_argument("--top-fraction", type=float, help="Retain this fraction (0,1] of ranked islands.")
    parser.add_argument("--window", type=int, default=200, help="Centered window width before chromosome-end clipping (default: 200).")
    parser.add_argument("--merge-distance", type=int, default=50, help="Merge same-chromosome segments whose half-open gap is <= this value (default: 50).")
    parser.add_argument("--include-chroms", default="canonical", help="Comma-separated normalized allowlist or 'canonical' for 2-22,X,Y (default).")
    parser.add_argument("--exclude-chroms", default="", help="Comma-separated names removed after --include-chroms.")
    threshold = parser.add_mutually_exclusive_group()
    threshold.add_argument("--min-score", type=float, help="Absolute minimum candidate coverage-column score.")
    threshold.add_argument("--score-percentile", type=float, help="Nearest-rank percentile over chromosome-filtered candidate scores.")
    parser.add_argument("--score-column", type=int, help="1-based coverage column; default col4 for bedGraph and col5 for BED (col4 for BED4).")
    parser.add_argument("--fasta-tool", choices=("internal", "bedtools"), default="internal", help="FASTA extraction implementation (default: internal .fai reader).")
    parser.add_argument("--allow-missing-chroms", action="store_true", help="Warn and skip candidate chromosomes absent from BigWig or .fai instead of failing.")
    parser.add_argument("--provenance", type=Path, help="Provenance JSON path; a flat TSV sibling is also written.")
    parser.add_argument("--prefix", default="cutandrun_regions", help="Output basename stem (default: cutandrun_regions).")
    return parser


def _validate_paths(args: argparse.Namespace, fai_path: Path) -> None:
    for label, path in (
        ("BigWig", args.bigwig),
        ("genome FASTA", args.genome_fasta),
        ("FASTA index", fai_path),
    ):
        if not path.is_file():
            raise RegionSelectionError(f"{label} does not exist: {path}")
    if args.candidates is None:
        raise RegionSelectionError(
            "--candidates is required for this pilot; scanning a BigWig to derive candidates is out of scope"
        )
    if not args.candidates.is_file():
        raise RegionSelectionError(f"candidate file does not exist: {args.candidates}")
    if Path(args.prefix).name != args.prefix or args.prefix in {"", ".", ".."}:
        raise RegionSelectionError("--prefix must be a non-empty basename, not a path")


def run(args: argparse.Namespace, argv: Sequence[str]) -> dict[str, Path]:
    fai_path = args.fai or Path(f"{args.genome_fasta}.fai")
    _validate_paths(args, fai_path)
    included, excluded = resolve_chromosome_policy(args.include_chroms, args.exclude_chroms)
    all_segments, score_column = read_candidate_segments(args.candidates, args.score_column)
    chromosome_filtered = filter_segments_by_chromosome(all_segments, included)
    fai_records = read_fai(fai_path)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    stem = args.prefix
    filtered_path = args.output_dir / f"{stem}.filtered_segments.bed"
    merged_path = args.output_dir / f"{stem}.merged_islands.bed"
    ranked_path = args.output_dir / f"{stem}.scored_ranked.bed"
    columns_path = args.output_dir / f"{stem}.scored_ranked.columns.txt"

    with PyBigWigScorer(args.bigwig) as scorer:
        usable_segments, missing_chromosomes = guard_assembly_chromosomes(
            chromosome_filtered,
            scorer.chromosomes,
            fai_records,
            allow_missing=args.allow_missing_chroms,
        )
        filtered, resolved_threshold = filter_segments_by_score(
            usable_segments,
            min_score=args.min_score,
            score_percentile=args.score_percentile,
        )
        merged = merge_segments(filtered, args.merge_distance)
        ranked = score_and_rank_islands(merged, scorer)
        pybigwig_version = scorer.version

    selected = select_top_islands(ranked, args.top_n, args.top_fraction)
    windows = centered_windows(
        selected,
        args.window,
        {name: record.length for name, record in fai_records.items()},
    )
    selected_count = len(selected)
    windows_path = args.output_dir / f"{stem}.top{selected_count}.centered_w{args.window}.bed"
    fasta_path = args.output_dir / f"{stem}.top{selected_count}.fasta"

    _write_filtered(filtered_path, filtered)
    _write_merged(merged_path, merged)
    _write_ranked(ranked_path, ranked)
    columns_path.write_text(
        "chrom\tstart\tend\tname\tbw_max\tbw_mean\tsupport\tlength\trank\n",
        encoding="utf-8",
    )
    _write_windows(windows_path, windows)

    bedtools_version: str | None = None
    if args.fasta_tool == "internal":
        write_internal_fasta(args.genome_fasta, fasta_path, windows, fai_records)
    else:
        bedtools = shutil.which("bedtools")
        if bedtools is None:
            raise RegionSelectionError(
                "--fasta-tool bedtools requested, but bedtools is not on PATH"
            )
        bedtools_version = run_bedtools_getfasta(
            bedtools, args.genome_fasta, windows_path, fasta_path
        )

    provenance_path = args.provenance or args.output_dir / f"{stem}.provenance.json"
    threshold_kind = (
        "min_score"
        if args.min_score is not None
        else "nearest_rank_percentile"
        if args.score_percentile is not None
        else "none"
    )
    concrete_chromosomes = sorted({segment.chrom for segment in filtered})
    outputs = {
        "filtered_segments": filtered_path,
        "merged_islands": merged_path,
        "scored_ranked": ranked_path,
        "scored_ranked_columns": columns_path,
        "centered_windows": windows_path,
        "fasta": fasta_path,
        "provenance_json": provenance_path,
        "provenance_tsv": provenance_path.with_suffix(".tsv"),
    }
    provenance: dict[str, object] = {
        "schema": "gentle.cutandrun_region_selection_provenance.v1",
        "created_at_utc": datetime.now(timezone.utc).isoformat().replace("+00:00", "Z"),
        "argv": list(argv),
        "git_commit": git_commit(),
        "coordinate_system": "0-based_half-open",
        "pipeline": ["chromosome_and_score_filter", "merge", "exact_bigwig_score", "rank", "center_and_clip", "fasta_extract"],
        "inputs": {
            "bigwig": fingerprint(args.bigwig),
            "candidates": fingerprint(args.candidates),
            "genome_fasta": fingerprint(args.genome_fasta),
            "fai": fingerprint(fai_path),
        },
        "assembly": {"genome_fasta": str(args.genome_fasta.resolve()), "fai": str(fai_path.resolve())},
        "chromosome_policy": {
            "include_argument": args.include_chroms,
            "exclude_argument": args.exclude_chroms,
            "normalized_include_after_exclude": sorted(included),
            "normalized_exclude": sorted(excluded),
            "concrete_chromosomes_used": concrete_chromosomes,
            "missing_chromosomes_skipped": missing_chromosomes,
            "allow_missing_chromosomes": args.allow_missing_chroms,
        },
        "candidate_filter": {
            "score_column_1based": score_column,
            "threshold_kind": threshold_kind,
            "min_score_argument": args.min_score,
            "score_percentile_argument": args.score_percentile,
            "resolved_threshold": resolved_threshold,
            "nearest_rank_definition": "sorted[ceil(P/100*N)-1], with rank at least 1",
        },
        "merge_distance": args.merge_distance,
        "window_size": args.window,
        "selection": {"top_n_argument": args.top_n, "top_fraction_argument": args.top_fraction, "selected_count": selected_count},
        "ranking_key": ["bw_max_desc", "bw_mean_desc", "support_desc", "length_desc", "chrom_asc", "start_asc", "end_asc"],
        "tools": {"pyBigWig_version": pybigwig_version, "pyBigWig_exact": True, "fasta_tool": args.fasta_tool, "bedtools_version": bedtools_version},
        "counts": {"input_segments": len(all_segments), "chromosome_filtered_segments": len(chromosome_filtered), "usable_segments": len(usable_segments), "threshold_filtered_segments": len(filtered), "merged_islands": len(merged), "selected_islands": selected_count},
        "outputs": {key: str(path.resolve()) for key, path in outputs.items()},
    }
    write_provenance(provenance_path, provenance)
    return outputs


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    parsed_argv = list(argv) if argv is not None else sys.argv[1:]
    args = parser.parse_args(parsed_argv)
    exact_argv = [sys.argv[0], *parsed_argv]
    try:
        outputs = run(args, exact_argv)
    except RegionSelectionError as exc:
        parser.error(str(exc))
    for label, path in outputs.items():
        print(f"{label}: {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
