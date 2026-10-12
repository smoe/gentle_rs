"""Optional real DuckDB I/O checks for hand-crafted regulatory teaching packages.

No GENtle process, PFM scoring or biological acceptance is performed. Inputs
come from the provenance-documented fixture generator; generated files stay in
fresh temporary directories. Missing DuckDB is an explicit skip, not a pass of
this boundary.
"""

import hashlib
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "regulatory_fixture", ROOT / "test_files/fixtures/regulatory_motif_subset/make_fixture.py")
fixture = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(fixture)
DUCKDB = shutil.which(os.environ.get("GENTLE_DUCKDB_BIN", "duckdb"))


@unittest.skipUnless(DUCKDB, "real DuckDB CLI unavailable; GENtle export acceptance is separate")
class RegulatoryMotifFixtureIoTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix="gentle-motif-fixture-")
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name) / "package with spaces and 'quote'"

    def query(self, sql, database=":memory:"):
        result = subprocess.run([DUCKDB, "-no-init", "-batch", "-bail", "-json", str(database), "-c", sql],
                                capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        return json.loads(result.stdout or "[]")

    def rows(self, relative):
        return self.query(f"SELECT * FROM read_parquet({fixture.quote(self.root / relative)})")

    def check_hashes(self):
        manifest = json.loads((self.root / "manifest.json").read_bytes())
        inventory = json.loads((self.root / "file_inventory.json").read_bytes())
        annotation = json.loads((self.root / "annotation/manifest.json").read_bytes())
        for base, records in [(self.root, manifest["files"]), (self.root, inventory),
                              (self.root / "annotation", annotation["files"])]:
            for record in records:
                if record.get("path") is None:
                    self.assertIsNone(record.get("sha256"))
                    continue
                data = (base / record["path"]).read_bytes()
                self.assertEqual(len(data), record["bytes"])
                self.assertEqual(hashlib.sha256(data).hexdigest(), record["sha256"])
        return manifest, inventory

    def test_default_logical_contract_is_not_implicitly_replaced_by_tss_preset(self):
        fixture.build(self.root, DUCKDB)
        manifest, inventory = self.check_hashes()
        self.assertEqual((manifest["genome_id"], manifest["assembly"], manifest["annotation_release"]),
                         ("fixture-genome", "fixture-assembly", "fixture-gtf"))
        self.assertEqual((manifest["upstream_bp"], manifest["downstream_bp"], manifest["rows"], manifest["motif_count"]),
                         (700, 300, 7, 2))
        self.assertEqual(manifest["chromosomes"], ["1", "2", "MT"])
        self.assertEqual(len(self.rows("annotation/promoter.parquet")), 4)
        self.assertEqual(len(self.rows("annotation/transcript_tss.parquet")), 6)
        self.assertEqual([record["rows"] for record in inventory], [4, 3, 0, 0, 0, 0])
        self.assertTrue(all(record["state"] == "known_empty_intersection" and record["path"] is None
                            for record in inventory if record["chrom"] == "MT"))
        tp73 = self.rows("hits/chrom=1/MA0861.2.parquet")
        patz1 = self.rows("hits/chrom=1/MA1961.2.parquet")
        self.assertEqual(sorted((row["start"], row["strand"], row["score"]) for row in tp73),
                         [(295, "+", -4.0), (295, "-", -3.5), (1200, "+", 0.0), (1600, "-", 6.0)])
        self.assertEqual(sorted((row["start"], row["strand"], row["score"]) for row in patz1),
                         [(320, "+", 1.25), (320, "-", -0.5), (1200, "+", 0.0)])
        self.assertEqual(self.rows("hits/chrom=2/MA0861.2.parquet"), [])
        genome = self.query("SELECT genome_id, assembly_name, assembly_accession FROM genome",
                            self.root / "regulatory_tfbs.duckdb")
        self.assertEqual(genome, [{"genome_id": "fixture-genome", "assembly_name": "fixture-assembly",
                                  "assembly_accession": "fixture-accession.1"}])

    def test_tss_preset_writes_compatible_full_width_signed_hits_and_refuses_clobber(self):
        fixture.build(self.root, DUCKDB, "tss-view")
        manifest, inventory = self.check_hashes()
        self.assertEqual((manifest["genome_id"], manifest["assembly"], manifest["rows"], manifest["motif_count"]),
                         ("synthetic-genome-v1", "synthetic-assembly-v1", 5, 1))
        self.assertEqual((manifest["upstream_bp"], manifest["downstream_bp"]), (3, 7))
        self.assertEqual(manifest["chromosomes"], ["synthetic_chr"])
        self.assertEqual(inventory[0]["rows"], 5)
        promoters = self.rows("annotation/promoter.parquet")
        self.assertEqual([(row["strand"], row["promoter_start"], row["promoter_end"], row["tss_start"])
                          for row in promoters], [("+", 96, 107, 99), ("-", 292, 303, 299)])
        hits = self.rows("hits/chrom=synthetic_chr/MA0004.1.parquet")
        self.assertEqual([(row["start"], row["strand"], row["score"], row["regulation_tags"]) for row in hits],
                         fixture.tss_view_hits())
        self.assertTrue(all(row["end"] - row["start"] == 6 and row["motif_id"] == "MA0004.1"
                            and row["score_mode"] == "log2_relative_risk" for row in hits))
        metadata = self.query("SELECT motif_id, motif_length FROM motif_metadata", self.root / "regulatory_tfbs.duckdb")
        self.assertEqual(metadata, [{"motif_id": "MA0004.1", "motif_length": 6}])
        before = {path.relative_to(self.root): path.read_bytes() for path in self.root.rglob("*") if path.is_file()}
        with self.assertRaises(FileExistsError):
            fixture.build(self.root, DUCKDB, "tss-view")
        self.assertEqual({path.relative_to(self.root): path.read_bytes() for path in self.root.rglob("*") if path.is_file()}, before)


if __name__ == "__main__":
    unittest.main()
