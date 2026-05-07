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

    def test_family_plotter_collects_tail_quantiles_but_hides_extremes_by_default(self) -> None:
        batch_summary = textwrap.dedent(
            """\
            sample_id\tsample_name\tstatus\treport_id\ttotal_reads\tseed_passed_reads\tall_q25_bp\tall_q50_bp\tall_q75_bp\tall_q90_bp\tall_q95_bp\tall_q99_bp\tall_q100_bp\tall_mean_bp\tstrict_seed_passed_max_read_bp
            S1\tSample 1\tok\ttp53_SRR1\t1000\t2\t500\t800\t1200\t2000\t3200\t90000\t120000\t1100\t1500
            S2\tSample 2\tok\ttp53_SRR2\t1000\t3\t600\t900\t1300\t2200\t3400\t95000\t125000\t1200\t1600
            """
        )

        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_path = Path(tmp_dir)
            input_tsv = tmp_path / "tp53_batch_summary.tsv"
            default_svg = tmp_path / "family_default.svg"
            default_out = tmp_path / "family_default.tsv"
            tail_svg = tmp_path / "family_tail.svg"
            tail_out = tmp_path / "family_tail.tsv"
            input_tsv.write_text(batch_summary, encoding="utf-8")

            subprocess.run(
                [
                    "python3",
                    str(REPO_ROOT / "scripts" / "plot_pancreas_gene_family.py"),
                    "--batch-summary",
                    f"TP53={input_tsv}",
                    "--genes",
                    "TP53",
                    "--output",
                    str(default_svg),
                    "--output-tsv",
                    str(default_out),
                ],
                check=True,
                cwd=REPO_ROOT,
            )
            subprocess.run(
                [
                    "python3",
                    str(REPO_ROOT / "scripts" / "plot_pancreas_gene_family.py"),
                    "--batch-summary",
                    f"TP53={input_tsv}",
                    "--genes",
                    "TP53",
                    "--show-all-read-max",
                    "--output",
                    str(tail_svg),
                    "--output-tsv",
                    str(tail_out),
                ],
                check=True,
                cwd=REPO_ROOT,
            )

            default_svg_text = default_svg.read_text(encoding="utf-8")
            tail_svg_text = tail_svg.read_text(encoding="utf-8")
            default_tsv_header = default_out.read_text(encoding="utf-8").splitlines()[0]
            self.assertIn(">q95<", default_svg_text)
            self.assertNotIn(">q99<", default_svg_text)
            self.assertNotIn(">q100<", default_svg_text)
            self.assertIn(">q99<", tail_svg_text)
            self.assertIn(">q100<", tail_svg_text)
            self.assertIn("all_q95_bp", default_tsv_header)
            self.assertIn("all_q99_bp", default_tsv_header)


if __name__ == "__main__":
    unittest.main()
