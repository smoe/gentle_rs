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

    def test_family_plotter_prefers_canonical_gene_screen_summary(self) -> None:
        canonical = textwrap.dedent(
            """\
            schema\tgene\trun_accession\tsample_id\tsample_name\tsource_kind\tsource_path\tanalysis_phase\treport_id\tseq_id\tseed_feature_id\tinput_path\ttotal_reads\tall_q0_bp\tall_q25_bp\tall_q50_bp\tall_q75_bp\tall_q90_bp\tall_q95_bp\tall_q99_bp\tall_q100_bp\tall_mean_bp\tseed_passed_reads\tseed_passed_per_million\tseed_passed_q90_bp\tseed_passed_q95_bp\tseed_passed_q99_bp\tseed_passed_max_bp\tseed_passed_mean_bp\taccepted_target_count\taccepted_target_per_million\taccepted_target_max_bp\taccepted_target_mean_bp\talign_selection
            gentle.rna_read_gene_screen_summary.v1\tTP53\tSRR1\tS1\tSample 1\trna_reads_batch_map\tbatch_report.json\twith_alignment\ttp53_SRR1\tseq_a\t1\treads.fa\t1000\t100\t200\t300\t400\t500\t550\t590\t600\t350\t10\t10000\t700\t750\t790\t800\t650\t8\t8000\t900\t710\tseed_passed
            gentle.rna_read_gene_screen_summary.v1\tTP63\tSRR1\tS1\tSample 1\tpancreas_gene_rna_screen\tfinal_summaries.json\tseed_only\ttp63_SRR1\tseq_b\t1\t\t1000\t100\t200\t300\t400\t500\t550\t590\t600\t350\t2\t2000\t450\t475\t490\t500\t460\t\t\t\t\tseed_passed
            """
        )

        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_path = Path(tmp_dir)
            input_tsv = tmp_path / "gene_screen_summary.tsv"
            output_svg = tmp_path / "family.svg"
            output_tsv = tmp_path / "family.tsv"
            input_tsv.write_text(canonical, encoding="utf-8")

            subprocess.run(
                [
                    "python3",
                    str(REPO_ROOT / "scripts" / "plot_pancreas_gene_family.py"),
                    "--canonical-summary",
                    str(input_tsv),
                    "--genes",
                    "TP53,TP63",
                    "--output",
                    str(output_svg),
                    "--output-tsv",
                    str(output_tsv),
                ],
                check=True,
                cwd=REPO_ROOT,
            )

            output_text = output_tsv.read_text(encoding="utf-8")
            svg = output_svg.read_text(encoding="utf-8")
            self.assertIn("schema\tgene\trun_accession", output_text.splitlines()[0])
            self.assertIn("gentle.rna_read_gene_screen_summary.v1\tTP53\tSRR1", output_text)
            self.assertIn("strict_seed_passed", output_text)
            self.assertIn("TP53 / TP63 seed-passed support", svg)


class PancreasGeneRnaScreenScriptTests(unittest.TestCase):
    def test_gene_screen_help_documents_seed_only_default(self) -> None:
        result = subprocess.run(
            [
                "bash",
                str(REPO_ROOT / "scripts" / "pancreas_gene_rna_screen.sh"),
                "help",
            ],
            check=True,
            cwd=REPO_ROOT,
            capture_output=True,
            text=True,
        )

        self.assertIn("--seed-only", result.stdout)
        self.assertIn(
            "Interpret reads and summarize seed evidence only (default)",
            result.stdout,
        )
        self.assertIn("--with-alignment", result.stdout)
        self.assertIn("--align-selection MODE", result.stdout)
        self.assertIn("stop after the seed phase", result.stdout)

    def test_gene_screen_script_defaults_to_seed_only_phase(self) -> None:
        script = (REPO_ROOT / "scripts" / "pancreas_gene_rna_screen.sh").read_text(
            encoding="utf-8"
        )

        self.assertIn('SCREEN_PHASE="${SCREEN_PHASE:-seed_only}"', script)
        self.assertIn('ALIGN_SELECTION="${ALIGN_SELECTION:-seed_passed}"', script)
        self.assertIn('--seed-only)\n        SCREEN_PHASE="seed_only"', script)
        self.assertIn('--with-alignment)\n        SCREEN_PHASE="with_alignment"', script)
        self.assertIn('if [ "$SCREEN_PHASE" = "seed_only" ]; then', script)
        self.assertIn('SCREEN_PHASE="with_alignment"', script)
        self.assertIn("--concatemer-limit is ignored in seed-only mode", script)
        self.assertIn('rna-reads align-report "$report_id"', script)

    def test_gene_screen_script_syntax_is_valid(self) -> None:
        subprocess.run(
            [
                "bash",
                "-n",
                str(REPO_ROOT / "scripts" / "pancreas_gene_rna_screen.sh"),
            ],
            check=True,
            cwd=REPO_ROOT,
        )


if __name__ == "__main__":
    unittest.main()
