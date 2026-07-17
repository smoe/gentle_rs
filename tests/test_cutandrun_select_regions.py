"""Hermetic tests for the CUT&RUN motif-pilot region selector.

All FASTA, FAI, BED, and optional BigWig inputs are synthetic and are created
inside each test's temporary directory. No external biological fixture or
ignored local data is required.
"""

from __future__ import annotations

import importlib.util
import json
import subprocess
import sys
import tempfile
import types
import unittest
from pathlib import Path
from unittest import mock


REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT_PATH = REPO_ROOT / "scripts" / "cutandrun_select_regions.py"
SPEC = importlib.util.spec_from_file_location("cutandrun_select_regions", SCRIPT_PATH)
assert SPEC is not None and SPEC.loader is not None
selector = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = selector
SPEC.loader.exec_module(selector)

try:
    import pyBigWig  # type: ignore[import-not-found]
except ImportError:
    pyBigWig = None


def write_synthetic_fasta(directory: Path) -> tuple[Path, Path, dict[str, str]]:
    sequences = {
        "1": "A" * 120,
        "2": ("acGT" * 45)[:180],
        "X": ("Tgca" * 20)[:80],
        "Y": ("GATTACA" * 10)[:70],
        "MT": ("CcAaG" * 10)[:50],
        "GL000001.1": ("ATgcCG" * 10)[:60],
    }
    fasta_path = directory / "synthetic.fa"
    fai_path = directory / "synthetic.fa.fai"
    fai_rows: list[str] = []
    with fasta_path.open("wb") as fasta:
        for name, sequence in sequences.items():
            fasta.write(f">{name}\n".encode("ascii"))
            offset = fasta.tell()
            line_bases = 7
            for position in range(0, len(sequence), line_bases):
                fasta.write(sequence[position : position + line_bases].encode("ascii"))
                fasta.write(b"\n")
            fai_rows.append(
                f"{name}\t{len(sequence)}\t{offset}\t{line_bases}\t{line_bases + 1}"
            )
    fai_path.write_text("\n".join(fai_rows) + "\n", encoding="utf-8")
    return fasta_path, fai_path, sequences


def synthetic_segments():
    return [
        selector.Segment("1", 5, 15, 12.0),
        selector.Segment("2", 5, 15, 8.0),
        selector.Segment("2", 40, 50, 9.0),
        selector.Segment("2", 110, 120, 7.0),
        selector.Segment("X", 70, 79, 6.0),
        selector.Segment("Y", 2, 8, 5.0),
        selector.Segment("MT", 3, 12, 11.0),
        selector.Segment("GL000001.1", 4, 14, 10.0),
    ]


def write_synthetic_candidates(path: Path) -> None:
    path.write_text(
        "\n".join(
            f"{row.chrom}\t{row.start}\t{row.end}\tsegment_{index}\t{row.score}"
            for index, row in enumerate(synthetic_segments(), start=1)
        )
        + "\n",
        encoding="utf-8",
    )


class CutAndRunRegionSelectionTests(unittest.TestCase):
    def test_help_documents_reproducibility_contract_without_pybigwig(self) -> None:
        result = subprocess.run(
            ["python3", str(SCRIPT_PATH), "--help"],
            cwd=REPO_ROOT,
            check=True,
            capture_output=True,
            text=True,
        )

        self.assertIn("0-based half-open", result.stdout)
        self.assertIn("2-22,X,Y", result.stdout)
        self.assertIn("excludes chromosome 1, MT", result.stdout)
        self.assertIn("exact=True", result.stdout)
        self.assertIn("--candidates", result.stdout)
        self.assertIn("out of scope", result.stdout)

    def test_default_chromosome_policy_filters_chr1_mt_and_alt_contigs(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            candidates = Path(tmp_dir) / "candidates.bed"
            write_synthetic_candidates(candidates)
            segments, score_column = selector.read_candidate_segments(candidates)
        included, excluded = selector.resolve_chromosome_policy("canonical", "")
        retained = selector.filter_segments_by_chromosome(segments, included)

        self.assertEqual(score_column, 5)
        self.assertEqual(excluded, set())
        self.assertEqual({row.chrom for row in retained}, {"2", "X", "Y"})
        self.assertNotIn("1", included)
        self.assertEqual(selector.normalize_chromosome("chrM"), "MT")
        self.assertEqual(selector.normalize_chromosome("chrX"), "X")

    def test_merge_uses_half_open_gap_and_score_length_support(self) -> None:
        included, _ = selector.resolve_chromosome_policy("canonical", "")
        retained = selector.filter_segments_by_chromosome(synthetic_segments(), included)
        merged = selector.merge_segments(retained, merge_distance=50)

        chrom2 = [row for row in merged if row.chrom == "2"]
        self.assertEqual([(row.start, row.end) for row in chrom2], [(5, 50), (110, 120)])
        self.assertEqual(chrom2[0].support, 8.0 * 10 + 9.0 * 10)
        self.assertEqual(chrom2[1].support, 7.0 * 10)

    def test_nearest_rank_threshold_is_computed_after_chromosome_filter(self) -> None:
        included, _ = selector.resolve_chromosome_policy("canonical", "")
        retained = selector.filter_segments_by_chromosome(synthetic_segments(), included)
        filtered, threshold = selector.filter_segments_by_score(
            retained, min_score=None, score_percentile=80.0
        )

        # Five retained values sort as 5,6,7,8,9; nearest rank ceil(.8*5)=4.
        self.assertEqual(threshold, 8.0)
        self.assertEqual([row.score for row in filtered], [8.0, 9.0])

    def test_centering_clips_at_chromosome_end_and_internal_fasta_preserves_case(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            root = Path(tmp_dir)
            fasta_path, fai_path, sequences = write_synthetic_fasta(root)
            records = selector.read_fai(fai_path)
            selected = [
                selector.ScoredIsland(
                    chrom="X",
                    start=70,
                    end=79,
                    name="island_x",
                    bw_max=4.0,
                    bw_mean=2.0,
                    support=54.0,
                    length=9,
                    rank=1,
                )
            ]
            windows = selector.centered_windows(
                selected, 40, {name: row.length for name, row in records.items()}
            )

            self.assertEqual((windows[0].start, windows[0].end), (54, 80))
            self.assertLess(windows[0].end - windows[0].start, 40)
            output = root / "windows.fa"
            selector.write_internal_fasta(fasta_path, output, windows, records)
            lines = output.read_text(encoding="ascii").splitlines()
            self.assertEqual(lines[0], ">X:54-80")
            self.assertEqual(lines[1], sequences["X"][54:80])
            self.assertTrue(any(character.islower() for character in lines[1]))

    def test_internal_fasta_extracts_multiple_records_with_expected_lengths(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            root = Path(tmp_dir)
            fasta_path, fai_path, sequences = write_synthetic_fasta(root)
            records = selector.read_fai(fai_path)
            windows = [
                selector.CenteredWindow("2", 3, 18, "a", 1),
                selector.CenteredWindow("Y", 60, 70, "b", 2),
            ]
            output = root / "selected.fa"
            selector.write_internal_fasta(fasta_path, output, windows, records)
            lines = output.read_text(encoding="ascii").splitlines()

            self.assertEqual(lines, [">2:3-18", sequences["2"][3:18], ">Y:60-70", sequences["Y"][60:70]])
            self.assertEqual([len(lines[1]), len(lines[3])], [15, 10])

    def test_ranking_is_deterministic_under_input_reordering(self) -> None:
        islands = [
            selector.MergedIsland("X", 20, 30, "x", 100.0),
            selector.MergedIsland("2", 20, 30, "two_late", 100.0),
            selector.MergedIsland("2", 10, 20, "two_early", 100.0),
            selector.MergedIsland("Y", 5, 25, "higher_max", 1.0),
        ]

        def scorer(chrom: str, start: int, end: int) -> tuple[float, float]:
            if (chrom, start, end) == ("Y", 5, 25):
                return 11.0, 1.0
            return 10.0, 2.0

        forward = selector.score_and_rank_islands(islands, scorer)
        reverse = selector.score_and_rank_islands(reversed(islands), scorer)

        expected = ["higher_max", "two_early", "two_late", "x"]
        self.assertEqual([row.name for row in forward], expected)
        self.assertEqual([row.name for row in reverse], expected)
        self.assertEqual([row.rank for row in forward], [1, 2, 3, 4])

    def test_pybigwig_adapter_uses_exact_stats_and_maps_none_to_zero(self) -> None:
        class FakeHandle:
            def __init__(self) -> None:
                self.calls: list[tuple[str, int, int, str, bool]] = []

            def chroms(self):
                return {"2": 100}

            def stats(self, chrom, start, end, *, type, exact):
                self.calls.append((chrom, start, end, type, exact))
                return [None] if type == "mean" else [3.5]

            def close(self) -> None:
                pass

        handle = FakeHandle()
        fake_module = types.SimpleNamespace(__version__="test", open=lambda _path: handle)
        with mock.patch.dict(sys.modules, {"pyBigWig": fake_module}):
            with selector.PyBigWigScorer(Path("synthetic.bw")) as scorer:
                observed = scorer("2", 10, 20)

        self.assertEqual(observed, (3.5, 0.0))
        self.assertEqual(
            handle.calls,
            [("2", 10, 20, "max", True), ("2", 10, 20, "mean", True)],
        )

    def test_full_pipeline_emits_parser_clean_outputs_and_provenance(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            root = Path(tmp_dir)
            fasta_path, fai_path, sequences = write_synthetic_fasta(root)
            candidates = root / "candidates.bed"
            write_synthetic_candidates(candidates)
            bigwig_path = root / "signal.bw"
            bigwig_path.write_bytes(b"synthetic BigWig placeholder for injected scorer\n")
            output_dir = root / "out"

            class FakeHandle:
                def __init__(self) -> None:
                    self.exact_arguments: list[bool] = []

                def chroms(self):
                    return {name: len(sequence) for name, sequence in sequences.items()}

                def stats(self, chrom, start, end, *, type, exact):
                    self.exact_arguments.append(exact)
                    base = float({"2": 8, "X": 6, "Y": 5}[chrom])
                    return [base + (1.0 if type == "max" else 0.0)]

                def close(self) -> None:
                    pass

            handle = FakeHandle()
            fake_module = types.SimpleNamespace(
                __version__="synthetic-test", open=lambda _path: handle
            )
            command = [
                "--bigwig",
                str(bigwig_path),
                "--candidates",
                str(candidates),
                "--genome-fasta",
                str(fasta_path),
                "--fai",
                str(fai_path),
                "--output-dir",
                str(output_dir),
                "--top-n",
                "3",
                "--prefix",
                "pilot",
            ]
            args = selector.build_parser().parse_args(command)
            with mock.patch.dict(sys.modules, {"pyBigWig": fake_module}):
                outputs = selector.run(args, [str(SCRIPT_PATH), *command])

            self.assertTrue(all(handle.exact_arguments))
            self.assertTrue(all(path.is_file() for path in outputs.values()))
            filtered_rows = (output_dir / "pilot.filtered_segments.bed").read_text(
                encoding="utf-8"
            ).splitlines()
            self.assertEqual({row.split("\t")[0] for row in filtered_rows}, {"2", "X", "Y"})
            self.assertTrue(all(not row.startswith("#") for row in filtered_rows))
            provenance = json.loads(
                (output_dir / "pilot.provenance.json").read_text(encoding="utf-8")
            )
            self.assertEqual(provenance["coordinate_system"], "0-based_half-open")
            self.assertEqual(
                provenance["chromosome_policy"]["concrete_chromosomes_used"],
                ["2", "X", "Y"],
            )
            self.assertEqual(provenance["tools"]["pyBigWig_version"], "synthetic-test")
            self.assertEqual(len(provenance["inputs"]["bigwig"]["sha256"]), 64)
            provenance_tsv = (output_dir / "pilot.provenance.tsv").read_text(
                encoding="utf-8"
            )
            self.assertIn("coordinate_system\t0-based_half-open", provenance_tsv)
            fasta_headers = [
                line
                for line in (output_dir / "pilot.top3.fasta").read_text(
                    encoding="ascii"
                ).splitlines()
                if line.startswith(">")
            ]
            self.assertEqual(len(fasta_headers), 3)

    def test_assembly_guard_names_both_sources_and_can_skip_missing(self) -> None:
        segments = [selector.Segment("chr2", 1, 5, 1.0)]
        fai = {"2": selector.FaiRecord("2", 10, 0, 10, 11)}
        with self.assertRaisesRegex(
            selector.RegionSelectionError, "BigWig header.*FASTA .fai"
        ):
            selector.guard_assembly_chromosomes(
                segments, {"2": 10}, fai, allow_missing=False
            )

        with mock.patch("sys.stderr"):
            usable, missing = selector.guard_assembly_chromosomes(
                segments, {"2": 10}, fai, allow_missing=True
            )
        self.assertEqual(usable, [])
        self.assertEqual(missing[0]["chromosome"], "chr2")

    @unittest.skipUnless(pyBigWig is not None, "pyBigWig is not installed")
    def test_end_to_end_with_synthetic_bigwig(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            root = Path(tmp_dir)
            fasta_path, fai_path, sequences = write_synthetic_fasta(root)
            candidates = root / "candidates.bed"
            candidates.write_text(
                "2\t5\t15\ta\t8\n"
                "2\t40\t50\tb\t9\n"
                "X\t70\t79\tx\t6\n"
                "1\t5\t15\tone\t20\n"
                "MT\t3\t12\tmt\t20\n",
                encoding="utf-8",
            )
            bigwig_path = root / "signal.bw"
            bigwig = pyBigWig.open(str(bigwig_path), "w")
            bigwig.addHeader([(name, len(sequence)) for name, sequence in sequences.items()])
            bigwig.addEntries(
                ["1", "2", "2", "X", "MT"],
                [5, 5, 40, 70, 3],
                ends=[15, 15, 50, 79, 12],
                values=[20.0, 8.0, 9.0, 6.0, 20.0],
            )
            bigwig.close()
            output_dir = root / "out"

            subprocess.run(
                [
                    "python3",
                    str(SCRIPT_PATH),
                    "--bigwig",
                    str(bigwig_path),
                    "--candidates",
                    str(candidates),
                    "--genome-fasta",
                    str(fasta_path),
                    "--fai",
                    str(fai_path),
                    "--output-dir",
                    str(output_dir),
                    "--top-n",
                    "2",
                    "--prefix",
                    "pilot",
                ],
                cwd=REPO_ROOT,
                check=True,
            )

            provenance = json.loads(
                (output_dir / "pilot.provenance.json").read_text(encoding="utf-8")
            )
            self.assertEqual(
                provenance["chromosome_policy"]["concrete_chromosomes_used"],
                ["2", "X"],
            )
            self.assertTrue(provenance["tools"]["pyBigWig_exact"])
            self.assertEqual(provenance["counts"]["selected_islands"], 2)
            self.assertTrue((output_dir / "pilot.top2.fasta").is_file())


if __name__ == "__main__":
    unittest.main()
