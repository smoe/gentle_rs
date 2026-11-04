"""Offline wiring guards, not Rust compilation, behavior or speed acceptance."""

from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]


def check_report_reader(source):
    report = source.split("    fn primer_specificity_report_from_blast_reports(", 1)[1]
    report = report.split("    fn cdna_assay_iupac_mismatch_count(", 1)[0]
    if report.count("PreparedSequenceRegionReader::new(&catalog)") != 1:
        raise ValueError("one interpretation-local reader required")
    if report.count("&mut subject_reader,") != 2:
        raise ValueError("both hit directions must share the reader")
    guard = "subject_reader.validate_unchanged().map_err("
    if guard not in report or report.index(guard) > report.index("Ok(PrimerSpecificityReport {"):
        raise ValueError("reference guard must precede report return")


class PrimerReferenceReaderSourceTests(unittest.TestCase):
    def test_shared_reader_and_final_guard_are_wired_to_the_real_interpreter(self):
        source = (ROOT / "src/engine/ops/operation_handlers.rs").read_text()
        check_report_reader(source)
        for before, after in (
            ("PreparedSequenceRegionReader::new(&catalog)", "PreparedSequenceRegionReader::uncached(&catalog)"),
            ("&mut subject_reader,", "&catalog,"),
            ("subject_reader.validate_unchanged().map_err(", "removed_reference_guard("),
        ):
            with self.assertRaises(ValueError):
                check_report_reader(source.replace(before, after, 1))

    def test_legacy_reference_fetch_and_outer_reference_validation_are_retained(self):
        source = (ROOT / "src/engine/ops/operation_handlers.rs").read_text()
        hit = source.split("    fn primer_specificity_hit_from_blast(", 1)[1]
        hit = hit.split("    pub(crate) fn primer_specificity_aggregate_warnings(", 1)[0]
        self.assertLess(hit.index("Self::primer_specificity_three_prime_mismatches("),
                        hit.index("Self::primer_specificity_full_alignment_from_blast_strings("))
        self.assertEqual(hit.count("Self::primer_specificity_full_alignment("), 2)
        multi = (ROOT / "src/engine/analysis/primer_specificity_multi.rs").read_text()
        multi = multi.split("    pub fn import_primer_pair_multi_reference_specificity(", 1)[1]
        multi = multi.split("        let mut summary =", 1)[0]
        self.assertEqual(multi.count("self.validate_primer_specificity_handoff_database(child)"), 2)

    def test_every_read_keeps_resolution_manifest_and_file_checks(self):
        source = (ROOT / "src/genomes.rs").read_text()
        read = source.split("    fn read_prepared_sequence_region(", 1)[1]
        read = read.split("    /// List prepared chromosome/contig lengths", 1)[0]
        checks = ("self.resolve_prepared_genome_id(", "Self::load_manifest(",
                  "Self::validate_manifest_files(", "cache.bind_manifest(", "cache.index(")
        positions = [read.index(check) for check in checks]
        self.assertEqual(positions, sorted(positions))
        self.assertIn("read_fasta_slice(", read)
        self.assertIn("load_fasta_index(Path::new(&manifest.fasta_index_path))?", read)

    def test_index_reuse_is_local_and_final_validation_is_content_based(self):
        source = (ROOT / "src/genomes/prepared_sequence_reader.rs").read_text()
        for forbidden in ("static ", "OnceLock", "Mutex", "thread_local!", "unsafe"):
            self.assertNotIn(forbidden, source)
        self.assertIn("cache: Some(SequenceIndexCache::default())", source)
        self.assertIn("let before = index_fingerprint(path)?;", source)
        self.assertIn("let entries = load_fasta_index(path)?;", source)
        self.assertIn("if index_fingerprint(path)? != before", source)
        final = source.rsplit("    fn validate_unchanged(", 1)[1]
        for check in ("self.invalid", "GenomeCatalog::validate_manifest_files(manifest)?",
                      "manifest_fingerprint(path)?", "index_fingerprint(path)?"):
            self.assertIn(check, final)
        self.assertNotIn("IndexStamp", final)

    def test_new_modules_have_explicit_existing_paths_and_synthetic_provenance(self):
        for parent, relative in (
            ("src/genomes.rs", "genomes/prepared_sequence_reader.rs"),
            ("src/genomes.rs", "genomes/prepared_sequence_reader_tests.rs"),
            ("src/engine/ops/operation_handlers.rs", "../analysis/primer_specificity_reader_tests.rs"),
        ):
            source = ROOT / parent
            self.assertIn(f'#[path = "{relative}"]', source.read_text())
            target = (source.parent / relative).resolve()
            self.assertTrue(target.is_file())
            self.assertTrue(target.read_text().startswith("//!"))
            if target.name.endswith("_tests.rs"):
                self.assertIn("Synthetic", target.read_text())


if __name__ == "__main__":
    unittest.main()
