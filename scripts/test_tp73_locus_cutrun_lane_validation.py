#!/usr/bin/env python3
"""Focused tests for chromosome-general BigWig conversion and locus-lane validation."""

from __future__ import annotations

from copy import deepcopy
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

import pyBigWig


ROOT = Path(__file__).resolve().parent


def load_module(name: str, filename: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / filename)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


CONVERT = load_module("bigwig_to_bedgraph", "bigwig_to_bedgraph.py")
VALIDATE = load_module("validate_tp73_locus_cutrun_lanes", "validate_tp73_locus_cutrun_lanes.py")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_bigwig(path: Path, chromosome: str, start: int, end: int) -> None:
    bigwig = pyBigWig.open(str(path), "w")
    bigwig.addHeader([("7", 1000), ("11", 1000), ("19", 1000)])
    bigwig.addEntries([chromosome], [start], ends=[end], values=[2.5])
    bigwig.close()


def gene_fixture(root: Path, gene: str, chromosome: str, strand: str) -> tuple[Path, Path]:
    tracks = []
    lanes = []
    for index in range(12):
        source = root / f"{gene}_{index}.bigWig"
        write_bigwig(source, chromosome, 100 + index, 110 + index)
        source_id = f"source:{gene}:{index}"
        track_name = f"track {gene} {index}"
        tracks.append({"source_kind": "big_wig", "path": str(source),
                       "track_name": track_name, "source_id": source_id,
                       "assembly": "GRCh38"})
        genomic_start = 101 + index
        genomic_end = 110 + index
        if strand == "+":
            local_start, local_end = genomic_start - 51 + 1, genomic_end - 51 + 1
        else:
            local_start, local_end = 250 - genomic_end + 1, 250 - genomic_start + 1
        lanes.append({
            "state": "available", "source_id": source_id,
            "source_sha256": f"sha256:{sha256(source)}", "source_assembly": "GRCh38",
            "lane": {"track_name": track_name, "interval_count": 1, "intervals": [{
                "genomic_start_1based": genomic_start, "genomic_end_1based": genomic_end,
                "local_start_1based": local_start, "local_end_1based": local_end,
                "score": 2.5,
            }]},
        })
    request = {"schema": VALIDATE.REQUEST_SCHEMA, "gene_query": gene,
               "assembly": "GRCh38", "local_tracks": tracks}
    report = {"schema": VALIDATE.LOCUS_SCHEMA, "gene_symbol": gene,
              "gene_strand": strand,
              "locus_genomic_start_1based": 51, "locus_genomic_end_1based": 250,
              "sequence_binding": {"genome_anchor": {
                  "genome_id": "GRCh38", "chromosome": chromosome,
                  "start_1based": 51, "end_1based": 250, "strand": strand}},
              "occupancy_groups": [{"lanes": lanes[index:index + 3]}
                                   for index in range(0, 12, 3)]}
    request_path = root / f"{gene}.request.json"
    report_path = root / f"{gene}.report.json"
    request_path.write_text(json.dumps(request))
    report_path.write_text(json.dumps(report))
    return request_path, report_path


class BigWigConversionTests(unittest.TestCase):
    def test_converter_emits_all_chromosomes_deterministically(self):
        with TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "three.bigWig"
            bigwig = pyBigWig.open(str(source), "w")
            bigwig.addHeader([("7", 1000), ("11", 1000), ("19", 1000)])
            bigwig.addEntries(["7", "11", "19"], [10, 20, 30],
                              ends=[15, 25, 35], values=[1.0, 2.0, 3.0])
            bigwig.close()
            first, second = root / "first.bg", root / "second.bg"
            first.touch()
            self.assertEqual(CONVERT.convert(source, first), 3)
            self.assertEqual(CONVERT.convert(source, second), 3)
            self.assertEqual(first.read_bytes(), second.read_bytes())
            self.assertEqual([line.split("\t")[0] for line in first.read_text().splitlines()],
                             ["7", "11", "19"])


class LaneValidationTests(unittest.TestCase):
    def fixtures(self, root: Path):
        pairs = {}
        for gene, chromosome, strand in (("CD44", "11", "+"),
                                         ("TGFB1", "19", "-"),
                                         ("SERPINE1", "7", "+")):
            pairs[gene] = gene_fixture(root, gene, chromosome, strand)
        return pairs

    def test_all_three_genes_and_36_lanes_validate(self):
        with TemporaryDirectory() as temporary:
            root = Path(temporary)
            pairs = self.fixtures(root)
            genes = [VALIDATE.validate_gene(gene, report, request, "GRCh38")
                     for gene, (request, report) in pairs.items()]
            self.assertEqual(sum(row["lane_count"] for row in genes), 36)
            self.assertEqual(sum(row["rendered_interval_count"] for row in genes), 36)

    def test_missing_intervals_and_wrong_hash_are_rejected(self):
        with TemporaryDirectory() as temporary:
            root = Path(temporary)
            request, report_path = self.fixtures(root)["CD44"]
            report = json.loads(report_path.read_text())
            lane = report["occupancy_groups"][0]["lanes"][0]
            for mutation in ("missing", "hash"):
                changed = deepcopy(report)
                target = changed["occupancy_groups"][0]["lanes"][0]
                if mutation == "missing":
                    target["state"] = "no_compatible_interval"
                    target["lane"]["interval_count"] = 0
                    target["lane"]["intervals"] = []
                else:
                    target["source_sha256"] = "sha256:" + "0" * 64
                changed_path = root / f"changed-{mutation}.json"
                changed_path.write_text(json.dumps(changed))
                with self.subTest(mutation=mutation), self.assertRaises(SystemExit):
                    VALIDATE.validate_gene("CD44", changed_path, request, "GRCh38")
            self.assertEqual(lane["state"], "available")

    def test_catalog_tokens_near_substrings_and_other_releases_are_rejected(self):
        with TemporaryDirectory() as temporary:
            root = Path(temporary)
            request_path, report_path = self.fixtures(root)["TGFB1"]
            for assembly in ("Human", "Ensembl", "116", "GRCh3", "GRCh37",
                             "Human GRCh38 Ensembl 115"):
                with self.subTest(assembly=assembly), self.assertRaises(SystemExit):
                    VALIDATE.validate_gene("TGFB1", report_path, request_path, assembly)

    def test_fabricated_signal_coordinates_and_duplicates_fail_native_comparison(self):
        with TemporaryDirectory() as temporary:
            root = Path(temporary)
            request, report_path = gene_fixture(root, "CD44", "11", "+")
            report = json.loads(report_path.read_text())
            for defect in ("score", "location", "duplicate"):
                changed = deepcopy(report)
                lane = changed["occupancy_groups"][0]["lanes"][0]["lane"]
                interval = lane["intervals"][0]
                if defect == "score":
                    interval["score"] = 999.0
                elif defect == "location":
                    for key in ("genomic_start_1based", "genomic_end_1based",
                                "local_start_1based", "local_end_1based"):
                        interval[key] += 50
                else:
                    lane["intervals"].append(deepcopy(interval))
                    lane["interval_count"] = 2
                report_path.write_text(json.dumps(changed))
                with self.subTest(defect=defect), self.assertRaisesRegex(SystemExit, "native BigWig"):
                    VALIDATE.validate_gene("CD44", report_path, request, "GRCh38")

    def test_projection_uses_sequence_anchor_not_negative_gene_strand(self):
        with TemporaryDirectory() as temporary:
            root = Path(temporary)
            request, report_path = gene_fixture(root, "TGFB1", "19", "+")
            report = json.loads(report_path.read_text())
            report["gene_strand"] = "-"
            report_path.write_text(json.dumps(report))
            self.assertEqual(VALIDATE.validate_gene("TGFB1", report_path, request, "GRCh38")["lane_count"], 12)

    def test_clipping_filters_and_importer_precision_preserve_valid_signal(self):
        with TemporaryDirectory() as temporary:
            root = Path(temporary)
            request_path, report_path = gene_fixture(root, "CD44", "11", "+")
            request = json.loads(request_path.read_text())
            report = json.loads(report_path.read_text())
            track = request["local_tracks"][0]
            source = Path(track["path"])
            bigwig = pyBigWig.open(str(source), "w")
            bigwig.addHeader([("11", 1000)])
            bigwig.addEntries(["11", "11"], [80, 130], ends=[120, 140], values=[2.3456789, 9.0])
            bigwig.close()
            # The first interval is clipped by the inspected locus; the second is filtered out.
            track["max_score"] = 3.0
            report["locus_genomic_start_1based"] = 101
            lane = report["occupancy_groups"][0]["lanes"][0]
            lane["source_sha256"] = f"sha256:{sha256(source)}"
            interval = lane["lane"]["intervals"][0]
            interval.update(genomic_start_1based=101, genomic_end_1based=120,
                            local_start_1based=51, local_end_1based=70, score=2.345679)
            request_path.write_text(json.dumps(request))
            report_path.write_text(json.dumps(report))
            self.assertEqual(VALIDATE.validate_gene("CD44", report_path, request_path, "GRCh38")["lane_count"], 12)
            del track["max_score"]
            request_path.write_text(json.dumps(request))
            with self.assertRaisesRegex(SystemExit, "native BigWig"):
                VALIDATE.validate_gene("CD44", report_path, request_path, "GRCh38")


if __name__ == "__main__":
    unittest.main()
