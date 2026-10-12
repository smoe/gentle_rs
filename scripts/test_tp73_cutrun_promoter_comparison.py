#!/usr/bin/env python3
"""Historical bundle integrity and hand-crafted, offline producer regressions.

Synthetic fixtures are created below in temporary directories: repeated DNA,
toy transcript rows, mock signal values and HSPs, never biological evidence.
Recreate with `python3 -m unittest scripts.test_tp73_cutrun_promoter_comparison`.
"""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from scripts import prepare_tp73_cutrun_promoter_candidates as prepare
from scripts import prepare_transcript_promoterome as promoterome
from scripts import compare_candidates_to_promoterome as compare
from scripts import render_tp73_cutrun_promoter_comparison as render


ROOT = (
    Path(__file__).resolve().parents[1]
    / "docs/examples/regulatory_region_comparison/tp73_cutrun_supported_first_comparison"
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


class Tp73PromoterComparisonExampleTests(unittest.TestCase):
    def test_candidates_retain_seven_supported_tss_windows(self) -> None:
        payload = json.loads((ROOT / "candidate_regions.json").read_text())
        self.assertEqual(
            payload["schema"], "gentle.regulatory_region_comparison_sequences.v1"
        )
        regions = payload["regions"]
        self.assertEqual(len(compare.validate_candidate_sequences(payload, ROOT / "candidate_regions.fa")), 7)
        self.assertEqual(len(regions), 7)
        self.assertEqual(
            {row["gene_query"] for row in regions}, {"CD44", "TGFB1", "SERPINE1"}
        )
        for row in regions:
            evidence = row["cutrun_support"]["evidence"]
            self.assertTrue(
                any(
                    max(item["experimental_minus_matched_gfp"].values()) > 0
                    for item in evidence
                )
            )

    def test_summary_preserves_historical_accounting_without_endorsing_completeness(self) -> None:
        with (ROOT / "summary.tsv").open(newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(len(rows), 7)
        self.assertTrue(all(row["target_cap_reached"] == "False" for row in rows))
        distal = next(
            row for row in rows
            if row["gene"] == "SERPINE1" and row["tss_1based"] == "101122158"
        )
        self.assertEqual(distal["other_genes_25pct"], "27563")
        self.assertEqual(distal["other_genes_50pct"], "0")

    def test_receipt_binds_retained_outputs(self) -> None:
        receipt = json.loads((ROOT / "receipt.json").read_text())
        self.assertFalse(receipt["target_cap_reached"])
        for name in ("candidate_regions.json", "candidate_regions.fa"):
            self.assertEqual(sha256(ROOT / name), receipt["inputs"][name], name)
        for name, expected in receipt["outputs"].items():
            self.assertEqual(sha256(ROOT / name), expected, name)


def write_json(path, value):
    path.write_text(json.dumps(value, sort_keys=True) + "\n", encoding="utf-8")


def write_tsv(path, rows, columns=None):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns or list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def reference_fixture(root):
    root.mkdir()
    windows, mappings, sequences = [], [], []
    for strand, start, end, tss, tx in [("+", 1000, 3201, 3001, "T1"), ("-", 5800, 8001, 6001, "T2")]:
        key = promoterome.window_id("1", start, end, strand)
        windows.append({"promoter_id": key, "chromosome": "1", "start_0based": start,
                        "end_0based_exclusive": end, "strand": strand, "tss_1based": tss,
                        "boundary_clipped": "false"})
        mappings.append({"promoter_id": key, "gene_id": "G1", "gene_name": "TOY", "transcript_id": tx})
        sequences.append(f'>{key}\n{"A" * 2201}\n')
    write_tsv(root / "promoter_windows.tsv", windows)
    write_tsv(root / "promoter_transcripts.tsv", mappings)
    (root / "promoter_windows.fa").write_text("".join(sequences))
    receipt = {"schema": promoterome.RECEIPT_SCHEMA, "genome_id": "Human GRCh38 Ensembl 116",
               "upstream_bp": 2000, "downstream_bp": 200, "unique_promoter_window_count": 2,
               "included_transcript_count": 2,
               "sequence_orientation": "transcript_5prime_to_3prime_via_bedtools_strand",
               "artifacts": {path.name: compare.sha256_file(path) for path in root.iterdir()}}
    write_json(root / "receipt.json", receipt)
    return receipt, windows, mappings


class CandidatePreparationTests(unittest.TestCase):
    def test_reference_and_selected_geometry_validates_both_strands(self):
        with TemporaryDirectory() as directory, patch.object(prepare, "SELECTED_TRANSCRIPTS", {"TOY": {"T1", "T2"}}):
            root = Path(directory) / "reference"
            reference_fixture(root)
            _, windows, mappings, sequences = prepare.selected_inputs(root)
            self.assertEqual({row["strand"] for row in windows.values()}, {"+", "-"})
            self.assertEqual(len(mappings), 2)
            self.assertEqual({len(sequence) for sequence in sequences.values()}, {2201})

    def test_wrong_reference_and_unbound_edits_fail_closed(self):
        with TemporaryDirectory() as directory, patch.object(prepare, "SELECTED_TRANSCRIPTS", {"TOY": {"T1", "T2"}}):
            root = Path(directory) / "reference"
            receipt, _, _ = reference_fixture(root)
            receipt["genome_id"] = "different assembly"
            write_json(root / "receipt.json", receipt)
            with self.assertRaisesRegex(RuntimeError, "requires Human"):
                prepare.selected_inputs(root)
            receipt["genome_id"] = "Human GRCh38 Ensembl 116"
            write_json(root / "receipt.json", receipt)
            with (root / "promoter_windows.fa").open("a") as handle:
                handle.write("A\n")
            with self.assertRaisesRegex(RuntimeError, "hash mismatch"):
                prepare.selected_inputs(root)

    def test_missing_mapping_or_false_negative_strand_geometry_is_rejected_even_if_rehashed(self):
        for defect in ("missing", "wrong_gene", "wrong_geometry"):
            with self.subTest(defect=defect), TemporaryDirectory() as directory, patch.object(
                    prepare, "SELECTED_TRANSCRIPTS", {"TOY": {"T1", "T2"}}):
                root = Path(directory) / "reference"
                receipt, windows, mappings = reference_fixture(root)
                if defect == "missing":
                    mappings.pop()
                    receipt["included_transcript_count"] = 1
                elif defect == "wrong_gene":
                    mappings[0]["gene_name"] = "WRONG"
                else:
                    windows[1]["tss_1based"] += 1
                write_tsv(root / "promoter_windows.tsv", windows)
                write_tsv(root / "promoter_transcripts.tsv", mappings)
                receipt["artifacts"] = {name: compare.sha256_file(root / name) for name in receipt["artifacts"]}
                write_json(root / "receipt.json", receipt)
                with self.assertRaises(RuntimeError):
                    prepare.selected_inputs(root)

    def test_preparation_writes_bound_candidates_without_bigwig_dependency(self):
        class Signal:
            def __init__(self, name):
                self.value = 0.0 if "GFP" in name else 1.0
                self.closed = False

            def values(self, chromosome, start, end):
                return [self.value] * (end - start)

            def close(self):
                self.closed = True

        with TemporaryDirectory() as directory, patch.object(prepare, "SELECTED_TRANSCRIPTS", {"TOY": {"T1", "T2"}}), patch.object(
                prepare.subprocess, "check_output", return_value="a" * 40 + "\n"):
            root = Path(directory)
            reference, _, _ = reference_fixture(root / "reference")
            indexes = root / "reference/indexes"
            indexes.mkdir()
            for suffix in ("nhr", "nin", "nsq"):
                path = indexes / f"promoterome.{suffix}"
                path.write_bytes(b"synthetic index; BLAST is mocked")
                reference["artifacts"][f"indexes/{path.name}"] = compare.sha256_file(path)
            write_json(root / "reference/receipt.json", reference)
            tracks = root / "tracks"
            tracks.mkdir()
            for conditions in prepare.TRACKS.values():
                for filename in conditions.values():
                    (tracks / filename).write_bytes(b"synthetic mocked signal; not a BigWig")
            args = SimpleNamespace(promoterome=root / "reference", track_root=tracks,
                                   output=root / "candidates", source_revision="a" * 40)
            handles = []

            def open_signal(name):
                handle = Signal(name)
                handles.append(handle)
                return handle

            prepare.prepare(args, open_bigwig=open_signal)
            payload = json.loads((args.output / "candidate_regions.json").read_text())
            self.assertEqual(len(payload["regions"]), 2)
            self.assertEqual(len(payload["sequence_equivalence_classes"]), 1)
            self.assertTrue(all(handle.closed for handle in handles))
            self.assertEqual(payload["source_bindings"]["promoterome_receipt_sha256"], sha256(root / "reference/receipt.json"))
            compare.validate_candidate_sequences(payload, args.output / "candidate_regions.fa")
            representative = payload["sequence_equivalence_classes"][0]["representative_region_id"]
            target = payload["regions"][0]["promoterome_id"]

            def fake_blast(command, **kwargs):
                self.assertIn("-max_hsps", command)
                Path(command[command.index("-out") + 1]).write_text(
                    f"{representative}\t{target}\t90\t45\t0\t0\t1\t45\t1\t45\t1e-9\t60\n")

            comparison_args = SimpleNamespace(candidates_json=args.output / "candidate_regions.json",
                                              query_fasta=args.output / "candidate_regions.fa",
                                              promoterome=args.promoterome, output=args.output / "task",
                                              task=["blastn"], blastn="mock-blastn", timeout=1,
                                              max_target_seqs=100, max_hsps=1,
                                              min_alignment_bp=40, min_identity_pct=80, max_evalue=1e-5)
            with patch.object(compare, "run_command", side_effect=fake_blast):
                compare.compare(comparison_args)
            comparison = json.loads((args.output / "task/comparison.json").read_text())
            summaries = comparison["tasks"]["blastn"]["queries"]
            self.assertEqual(len(summaries), 2)
            self.assertTrue(all(row["hsp_cap_reached"] for row in summaries))
            self.assertTrue(all(row["other_promoters"]["distinct_genes"] == 0 for row in summaries))
            rendered = render.render(args.output, root / "render", "task")
            self.assertTrue(rendered["counts_are_lower_bounds"])
            self.assertTrue((root / "render/report.md").is_file())
            with self.assertRaisesRegex(RuntimeError, "absent or empty"):
                prepare.prepare(args, open_bigwig=open_signal)

    def test_modified_or_extra_blast_index_files_are_rejected(self):
        with TemporaryDirectory() as directory:
            root = Path(directory) / "reference"
            receipt, _, _ = reference_fixture(root)
            (root / "indexes").mkdir()
            path = root / "indexes/promoterome.nin"
            path.write_bytes(b"synthetic index")
            receipt["artifacts"]["indexes/promoterome.nin"] = compare.sha256_file(path)
            write_json(root / "receipt.json", receipt)
            promoterome.validate_promoterome(root, include_blast=True)
            path.write_bytes(b"different synthetic index")
            with self.assertRaisesRegex(RuntimeError, "hash mismatch"):
                promoterome.validate_promoterome(root, include_blast=True)
            (root / "indexes/promoterome.extra").write_bytes(b"unbound")
            with self.assertRaisesRegex(RuntimeError, "files differ"):
                promoterome.validate_promoterome(root, include_blast=True)


def rendering_fixture(root, *, hits=True, lower_bounds=False, minimum_bp=40):
    root.mkdir()
    run = root / "task"
    run.mkdir()
    digest = "sha256:" + hashlib.sha256(b"A" * 100).hexdigest()
    candidate = {"schema": "gentle.regulatory_region_comparison_sequences.v1", "source_revision": "synthetic",
                 "source_bindings": {"promoterome_receipt_sha256": "f" * 64},
                 "regions": [{"region_id": "q", "gene_query": "TOY", "transcript_ids": ["T1"],
                              "sequence_length_bp": 100, "sequence_sha256": digest,
                              "genome_extraction": {"genome_id": "Synthetic reference", "gene_id": "G1",
                                                    "promoter_upstream_bp": 90, "promoter_downstream_bp": 9,
                                                    "chromosome": "1", "start_1based": 1, "end_1based": 100,
                                                    "strand": "+", "tss_1based": 91},
                              "cutrun_support": {"evidence": [{"cell_line": cell,
                                                                "experimental_minus_matched_gfp": {"TA": 1.0}}
                                                               for cell in ("SAOS-2", "SK-MEL-29-2")]}}],
                 "sequence_equivalence_classes": [{"representative_region_id": "q", "member_region_ids": ["q"],
                                                   "sequence_length_bp": 100, "sequence_sha256": digest}]}
    write_json(root / "candidate_regions.json", candidate)
    (root / "candidate_regions.fa").write_text(">q\n" + "A" * 100 + "\n")
    (run / "hits.blastn.tsv").write_text("q\tp\t90\t40\t0\t0\t1\t40\t1\t40\t1e-9\t60\n" if hits else "")
    rows = compare.parse_blast(run / "hits.blastn.tsv", minimum_bp, 85, 1e-6)
    summary = compare.summarize_task(rows, {"q": candidate["regions"][0]},
                                    {"p": {"chromosome": "2", "start_0based": "1000", "end_0based_exclusive": "1100",
                                           "strand": "+", "tss_1based": "1091"}},
                                    [{"promoter_id": "p", "gene_id": "G2", "gene_name": "PARTNER", "transcript_id": "T2"},
                                     {"promoter_id": "own", "gene_id": "G1", "gene_name": "TOY", "transcript_id": "T1"}],
                                    run / "matches.blastn.tsv")[0]
    summary.update(target_cap_reached=False, hsp_cap_reached=lower_bounds, counts_are_lower_bounds=lower_bounds)
    comparison = {"schema": compare.SCHEMA, "counting_policy_id": compare.COUNTING_POLICY,
                  "candidate_metadata_sha256": compare.sha256_file(root / "candidate_regions.json"),
                  "query_fasta_sha256": compare.sha256_file(root / "candidate_regions.fa"), "promoterome_receipt_sha256": "sha256:" + "f" * 64,
                  "background": {"genome_id": "Synthetic reference", "upstream_bp": 90, "downstream_bp": 9,
                                 "unique_promoter_window_count": 2},
                  "search_policy": {"dust": "yes", "soft_masking": True},
                  "thresholds": {"min_alignment_bp": minimum_bp, "min_identity_pct": 85, "max_evalue": 1e-6,
                                 "max_target_seqs": 10, "max_hsps_per_target": 1 if lower_bounds else 10},
                  "tasks": {"blastn": {"queries": [summary], "raw_hits_path": "hits.blastn.tsv",
                                       "raw_hits_sha256": compare.sha256_file(run / "hits.blastn.tsv"),
                                       "matches_path": "matches.blastn.tsv", "matches_sha256": compare.sha256_file(run / "matches.blastn.tsv")}}}
    write_json(run / "comparison.json", comparison)
    return comparison


class ComparisonRenderingTests(unittest.TestCase):
    def test_report_uses_inputs_and_propagates_lower_bounds_to_every_surface(self):
        with TemporaryDirectory() as directory:
            root, output = Path(directory) / "run", Path(directory) / "render"
            rendering_fixture(root, lower_bounds=True)
            receipt = render.render(root, output, "task")
            report = (output / "report.md").read_text()
            self.assertIn("2 transcript-linked", report)
            self.assertIn("PARTNER", report)
            self.assertIn("40.0%", report)
            self.assertIn("LOWER BOUNDS", report)
            self.assertNotIn("MAN2A2", report)
            self.assertNotIn("389,722", report)
            svg = (output / "first_comparison.svg").read_text()
            self.assertIn("-90/+9", svg)
            self.assertIn("85% identity", svg)
            self.assertIn("LOWER BOUNDS", svg)
            for name in ("summary.tsv", "top_matches.tsv"):
                with (output / name).open() as handle:
                    rows = list(csv.DictReader(handle, delimiter="\t"))
                self.assertTrue(all(row["counts_are_lower_bounds"] == "True" for row in rows))
            self.assertTrue(receipt["counts_are_lower_bounds"])
            self.assertEqual(receipt["inputs"]["raw_hits"], compare.sha256_file(root / "task/hits.blastn.tsv"))
            with self.assertRaisesRegex(RuntimeError, "absent or empty"):
                render.render(root, output, "task")

    def test_no_hits_and_changed_thresholds_render_empty_table_without_false_conclusions(self):
        with TemporaryDirectory() as directory:
            root, output = Path(directory) / "run", Path(directory) / "render"
            rendering_fixture(root, minimum_bp=50)
            render.render(root, output, "task")
            self.assertIn("No qualifying other-gene hit", (output / "report.md").read_text())
            self.assertNotIn("PARTNER", (output / "report.md").read_text())
            self.assertIn("≥50 bp", (output / "first_comparison.svg").read_text())
            self.assertEqual(len((output / "top_matches.tsv").read_text().splitlines()), 1)

    def test_mixed_or_historical_evidence_is_rejected_before_output(self):
        for filename in ("candidate_regions.json", "task/hits.blastn.tsv", "task/matches.blastn.tsv", "task/comparison.json"):
            with self.subTest(filename=filename), TemporaryDirectory() as directory:
                root, output = Path(directory) / "run", Path(directory) / "render"
                payload = rendering_fixture(root)
                if filename.endswith("comparison.json"):
                    payload.pop("counting_policy_id")
                    write_json(root / filename, payload)
                else:
                    with (root / filename).open("a") as handle:
                        handle.write("\n")
                with self.assertRaises(RuntimeError):
                    render.render(root, output, "task")
                self.assertFalse(output.exists())

    def test_order_breaks_respect_reverse_direction_and_orientation_changes(self):
        def hsp(q, s, reverse=False):
            return {"qstart": q, "qend": q + 19, "sstart": s + 19 if reverse else s,
                    "send": s if reverse else s + 19, "bitscore": 20}

        for blocks, expected in [([hsp(1, 1), hsp(30, 30)], [False, False]),
                                 ([hsp(81, 1, True), hsp(51, 30, True)], [False, False]),
                                 ([hsp(51, 1, True), hsp(81, 30, True)], [False, True]),
                                 ([hsp(1, 1), hsp(30, 30, True)], [False, True])]:
            self.assertEqual([broken for _, broken in render.ordered_blocks(blocks)], expected)

    def test_png_invocation_is_explicit_and_bound_in_receipt(self):
        with TemporaryDirectory() as directory:
            root, output = Path(directory) / "run", Path(directory) / "render"
            rendering_fixture(root)

            def fake_run(command, **kwargs):
                if "--version" not in command:
                    self.assertEqual(kwargs["cwd"], output)
                    self.assertIn("--format=png", command)
                    (output / "first_comparison.png").write_bytes(b"\x89PNG\r\n\x1a\nsynthetic renderer stub")
                return SimpleNamespace(stdout=b"synthetic renderer 1.0", stderr=b"")

            with patch.object(render.shutil, "which", return_value=__file__), patch.object(render.subprocess, "run", side_effect=fake_run):
                receipt = render.render(root, output, "task", png_renderer="mock-rsvg-convert")
            self.assertEqual(receipt["outputs"]["first_comparison.png"], sha256(output / "first_comparison.png"))
            self.assertEqual(receipt["png_rendering"]["version"], "synthetic renderer 1.0")
            self.assertEqual(receipt["png_rendering"]["executable_sha256"], sha256(Path(__file__)))


if __name__ == "__main__":
    unittest.main()
