#!/usr/bin/env python3
"""Read-only source authentication for this frozen fixture, not a read analysis.

The committed manifest authenticates expected bindings and the combined GTF.
Full source files are streamed once; no extraction or output creation occurs.
"""

import argparse
from collections import Counter
import csv
import hashlib
import io
from pathlib import Path, PurePosixPath
import posixpath
import re
import sys


FIELDS = ("role", "path", "size_bytes", "sha256_or_bound_digest", "note")
ROLES = Counter(source_cram=2, cohort_receipt=2, reference_fasta=1,
                annotation_gtf=1, sample_sheet=1)
ANNOTATION = "annotation/ensembl116.tp53_tp73.gtf"
GTF_PREFIX = b"# TP53 and TP73 records copied from the Ensembl 116 per-gene extracts.\n"


def _manifest(data):
    entries = {}
    for line in (data / "SHA256SUMS").read_bytes().decode("utf-8").splitlines():
        match = re.fullmatch(r"([0-9a-f]{64})  \./(.+)", line)
        if not match:
            raise ValueError("fixture_manifest_invalid")
        digest, name = match.groups()
        if (name in entries or "\\" in name
                or any(part in ("", ".", "..") for part in name.split("/"))):
            raise ValueError("fixture_manifest_path_invalid")
        entries[name] = digest
    if not {"source_bindings.tsv", ANNOTATION} <= entries.keys():
        raise ValueError("fixture_manifest_missing_binding")
    return entries


def verify_sources(fixture, storage_root, tp53_gtf, tp73_gtf):
    """Return actual verified bindings, preserving historical row order/notes.

    Storage relocation changes only the path prefix. Inputs must remain immutable
    during extraction; this preflight does not lock the archival filesystem.
    """
    data = Path(fixture) / "data"
    expected = _manifest(data)
    payload = (data / "source_bindings.tsv").read_bytes()
    if hashlib.sha256(payload).hexdigest() != expected["source_bindings.tsv"]:
        raise ValueError("fixture_bindings_hash_mismatch")
    reader = csv.DictReader(io.StringIO(payload.decode("utf-8")), delimiter="\t")
    rows = list(reader)
    if (reader.fieldnames != list(FIELDS)
            or any(set(row) != set(FIELDS) or any(v is None for v in row.values()) for row in rows)
            or Counter(row["role"] for row in rows) != ROLES):
        raise ValueError("fixture_bindings_invalid")
    paths = [row["path"] for row in rows]
    if (len(set(paths)) != len(paths)
            or any(not path.startswith("/") or "\\" in path
                   or any(part in ("", ".", "..") for part in path[1:].split("/"))
                   for path in paths)):
        raise ValueError("fixture_source_path_invalid")
    frozen_root = PurePosixPath(posixpath.commonpath(paths))
    if frozen_root == PurePosixPath("/"):
        raise ValueError("fixture_storage_root_invalid")

    verified = []
    for row in rows:
        if (not re.fullmatch(r"0|[1-9][0-9]*", row["size_bytes"])
                or not re.fullmatch(r"[0-9a-f]{64}", row["sha256_or_bound_digest"])):
            raise ValueError("fixture_source_binding_invalid")
        relative = PurePosixPath(row["path"]).relative_to(frozen_root)
        path = Path(storage_root).joinpath(*relative.parts)
        if not path.is_file():
            raise ValueError(f"source_missing_or_nonregular: {row['role']}: {path}")
        size = path.stat().st_size
        if size != int(row["size_bytes"]):
            raise ValueError(f"source_size_mismatch: {row['role']}: {path}")
        with path.open("rb") as stream:
            digest = hashlib.file_digest(stream, "sha256").hexdigest()
        if digest != row["sha256_or_bound_digest"]:
            raise ValueError(f"source_hash_mismatch: {row['role']}: {path}")
        verified.append(dict(row, path=str(path), size_bytes=str(size),
                             sha256_or_bound_digest=digest))

    annotation = hashlib.sha256(GTF_PREFIX)
    for path in (Path(tp53_gtf), Path(tp73_gtf)):
        if not path.is_file():
            raise ValueError(f"gene_gtf_missing_or_nonregular: {path}")
        with path.open("rb") as stream:
            for block in iter(lambda: stream.read(1024 * 1024), b""):
                annotation.update(block)
    if annotation.hexdigest() != expected[ANNOTATION]:
        raise ValueError("gene_gtf_combination_hash_mismatch")

    output = io.StringIO(newline="")
    writer = csv.DictWriter(output, FIELDS, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(verified)
    return output.getvalue()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("storage_root", type=Path)
    parser.add_argument("tp53_gtf", type=Path)
    parser.add_argument("tp73_gtf", type=Path)
    args = parser.parse_args()
    try:
        verified = verify_sources(Path(__file__).resolve().parent, args.storage_root,
                                  args.tp53_gtf, args.tp73_gtf)
    except (OSError, ValueError, csv.Error) as error:
        print(f"fixture_source_verification: {error}", file=sys.stderr)
        return 2
    sys.stdout.write(verified)
    return 0


if __name__ == "__main__":
    sys.exit(main())
