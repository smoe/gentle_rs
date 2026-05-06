import subprocess
import tempfile
import textwrap
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]


class PlotPancreasGeneScreenTests(unittest.TestCase):
    def test_generic_plotter_help_explains_seed_vs_accepted_columns(self) -> None:
        result = subprocess.run(
            [
                "python3",
                str(REPO_ROOT / "scripts" / "plot_pancreas_gene_screen.py"),
                "--help",
            ],
            check=True,
            cwd=REPO_ROOT,
            capture_output=True,
            text=True,
        )

        self.assertIn("strict_seed_passed_reads", result.stdout)
        self.assertIn("Reads that pass the first seed-matching phase only.", result.stdout)
        self.assertIn("This is the pre-alignment", result.stdout)
        self.assertIn("accepted_target_reads", result.stdout)
        self.assertIn("alignment-confirmed view", result.stdout)

    def test_tp73_wrapper_help_explains_seed_vs_accepted_columns(self) -> None:
        result = subprocess.run(
            [
                "python3",
                str(REPO_ROOT / "scripts" / "plot_tp73_pancreas_cohort.py"),
                "--help",
            ],
            check=True,
            cwd=REPO_ROOT,
            capture_output=True,
            text=True,
        )

        self.assertIn("strict_seed_passed_reads", result.stdout)
        self.assertIn("before any sequence-alignment", result.stdout)
        self.assertIn("accepted_tp73_reads", result.stdout)
        self.assertIn("after target-support confirmation", result.stdout)

    def test_generic_plotter_defaults_to_strict_seed_passed_reads(self) -> None:
        figure_tsv = textwrap.dedent(
            """\
            run_accession\tsample_id\tsample_name\ttotal_reads\tall_q0_bp\tall_q25_bp\tall_q50_bp\tall_q75_bp\tall_q90_bp\tall_q100_bp\tall_mean_bp\tstrict_seed_passed_reads\tstrict_seed_passed_max_read_bp\tstrict_seed_accepted_target_reads\taccepted_target_reads\taccepted_target_max_read_bp
            SRR1\tS1\tSample 1\t1000\t100\t200\t300\t400\t500\t600\t350\t10\t700\t50\t80\t900
            SRR2\tS2\tSample 2\t1200\t100\t220\t320\t420\t520\t620\t360\t12\t710\t60\t90\t950
            """
        )

        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_path = Path(tmp_dir)
            input_tsv = tmp_path / "e2f1_pancreas_figure_source.tsv"
            output_svg = tmp_path / "e2f1_pancreas_overview.svg"
            input_tsv.write_text(figure_tsv, encoding="utf-8")

            subprocess.run(
                [
                    "python3",
                    str(REPO_ROOT / "scripts" / "plot_pancreas_gene_screen.py"),
                    str(input_tsv),
                    "--output",
                    str(output_svg),
                ],
                check=True,
                cwd=REPO_ROOT,
            )

            svg = output_svg.read_text(encoding="utf-8")
            self.assertIn("E2F1 strict seed-passed reads", svg)
            self.assertNotIn("accepted target reads", svg)

    def test_tp73_wrapper_keeps_legacy_accepted_tp73_max_length_overlay(self) -> None:
        legacy_tsv = textwrap.dedent(
            """\
            run_accession\tsample_id\tsample_name\ttotal_reads\tall_q0_bp\tall_q25_bp\tall_q50_bp\tall_q75_bp\tall_q90_bp\tall_q100_bp\tall_mean_bp\tstrict_seed_passed_reads\tstrict_seed_passed_max_read_bp\taccepted_tp73_reads\taccepted_tp73_max_read_bp
            SRR1\tS1\tSample 1\t1000\t100\t200\t300\t400\t500\t600\t350\t10\t700\t8\t900
            SRR2\tS2\tSample 2\t1200\t100\t220\t320\t420\t520\t620\t360\t12\t710\t9\t950
            """
        )

        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_path = Path(tmp_dir)
            input_tsv = tmp_path / "legacy_tp73.tsv"
            output_svg = tmp_path / "legacy_tp73.svg"
            input_tsv.write_text(legacy_tsv, encoding="utf-8")

            subprocess.run(
                [
                    "python3",
                    str(REPO_ROOT / "scripts" / "plot_tp73_pancreas_cohort.py"),
                    str(input_tsv),
                    "--bar-column",
                    "accepted_tp73_reads",
                    "--output",
                    str(output_svg),
                ],
                check=True,
                cwd=REPO_ROOT,
            )

            svg = output_svg.read_text(encoding="utf-8")
            self.assertIn("max accepted TP73 read length (bp)", svg)
            self.assertIn('stroke="#ea580c"', svg)


if __name__ == "__main__":
    unittest.main()
