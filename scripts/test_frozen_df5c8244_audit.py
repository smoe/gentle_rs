"""Check owner-delivered frozen review/log bytes, without Rust or GUI execution."""

import hashlib
import json
from pathlib import Path
import re
import unittest


ROOT = Path(__file__).resolve().parents[1]
BASE = Path("docs/audits/glen_followups_20261008/candidate_df5c8244")
GLEN = BASE / "glen"
NATIVE = BASE / "native-x11-supplement"
SOURCE = "df5c8244787d8c920ee171ee55167de55e8a7294"
ORIGINAL_ROOT = "/home/clawbio/.openclaw/workspace/artifacts/gentle-df5c8244-verification-20261014/"
HANDOFF_ROOT = "/tmp/glen-df5c8244-small-evidence-20261014/"
HANDOFF_HASH = "51bc128734d18f6ead063be80efffc7a24a43806a4dc83b7dab9a4b8be53a61a"
NATIVE_HASH = "867e89d7e8ce5782707d45a572bd4d116d8a06e6427b2c2c2fe05880d42a6732"
EARLY_HASH = "0b6a9d95c689a7cb0f0ce2becf118b2a9c6192ed39738879f0b867cb39b16430"
FINAL_HASH = "713a8d034e2c7964d29e6350e3fa79d026726e7d28b15686d66597f5f9754e1a"
FILTERS = {
    "dbsnp_": 11,
    "materialize_variant_allele": 2,
    "reporter_construct_handoff_": 7,
    "variant_followup_allele_pair_": 2,
    "test_support__heavy_stack_tests": 2,
    "release_version_consistency": 4,
}


def entries(payload, prefix):
    """Translate historical paths lexically, never follow them on this host."""
    result = {}
    for line in payload.decode("utf-8").splitlines():
        digest, name = line.split("  ", 1)
        if not re.fullmatch(r"[0-9a-f]{64}", digest) or not name.startswith(prefix):
            raise ValueError("Unexpected historical manifest entry")
        name = name[len(prefix):]
        parts = name.split("/")
        if any(part in ("", ".", "..") or "\\" in part or ":" in part for part in parts):
            raise ValueError("Unsafe historical manifest path")
        path = Path(*parts)
        if path.is_absolute() or path in result:
            raise ValueError("Duplicate or absolute historical manifest path")
        result[path] = digest
    return result


def retained_hashes():
    original = entries((ROOT / GLEN / "HANDOFF-SHA256SUMS").read_bytes(), HANDOFF_ROOT)
    selected = [Path(name) for name in
                ("REPORT.md", "SUMMARY.json", "identity/environment.txt", "cleanup.txt")]
    for name in (*FILTERS, "cargo_check"):
        selected.extend(Path("rust-boundary") / f"{name}.{suffix}.txt"
                        for suffix in ("command", "exit-code", "stdout", "stderr", "time"))
    result = {GLEN / path: original[path] for path in selected}
    result[GLEN / "HANDOFF-SHA256SUMS"] = HANDOFF_HASH
    native = entries((ROOT / NATIVE / "HANDOFF-SHA256SUMS").read_bytes(), "./")
    result.update({NATIVE / path: digest for path, digest in native.items()})
    result[NATIVE / "HANDOFF-SHA256SUMS"] = NATIVE_HASH
    return result


class FrozenReviewTests(unittest.TestCase):
    def test_original_bytes_membership_and_manifest_anchors(self):
        hashes = retained_hashes()
        self.assertEqual(len(hashes), 62)
        for path, digest in hashes.items():
            with self.subTest(path=path):
                self.assertFalse((ROOT / path).is_symlink())
                self.assertEqual(hashlib.sha256((ROOT / path).read_bytes()).hexdigest(), digest)
        actual = {path.relative_to(ROOT) for base in (GLEN, NATIVE)
                  for path in (ROOT / base).rglob("*") if path.is_file()}
        self.assertEqual(actual, set(hashes))
        self.assertEqual(len(entries((ROOT / GLEN / "HANDOFF-SHA256SUMS").read_bytes(),
                                     HANDOFF_ROOT)), 94)

    def test_manifest_versions_differ_only_by_retained_cleanup(self):
        early_bytes = (ROOT / NATIVE / "manifests/early-0b6a9d95/SHA256SUMS").read_bytes()
        final_bytes = (ROOT / NATIVE / "manifests/final-713a8d03/SHA256SUMS").read_bytes()
        self.assertEqual(hashlib.sha256(early_bytes).hexdigest(), EARLY_HASH)
        self.assertEqual(hashlib.sha256(final_bytes).hexdigest(), FINAL_HASH)
        early, final = entries(early_bytes, ORIGINAL_ROOT), entries(final_bytes, ORIGINAL_ROOT)
        self.assertEqual(len(early), 194)
        self.assertEqual(len(final), 195)
        self.assertEqual(set(final) - set(early), {Path("cleanup.txt")})
        self.assertEqual({path: final[path] for path in early}, early)
        self.assertEqual(hashlib.sha256((ROOT / GLEN / "cleanup.txt").read_bytes()).hexdigest(),
                         final[Path("cleanup.txt")])
        for version, digest in (("early-0b6a9d95", EARLY_HASH), ("final-713a8d03", FINAL_HASH)):
            companion = (ROOT / NATIVE / "manifests" / version /
                         "original-SHA256SUMS.sha256").read_text()
            self.assertEqual(companion.strip(), digest + "  " + ORIGINAL_ROOT + "SHA256SUMS")
        cleanup_line = next(line for line in final_bytes.splitlines()
                            if line.endswith(b"/cleanup.txt"))
        self.assertEqual(b"\n".join(line for line in final_bytes.splitlines()
                                    if line != cleanup_line) + b"\n", early_bytes)

    def test_retained_native_logs_match_both_full_versions(self):
        for version in ("early-0b6a9d95", "final-713a8d03"):
            full = entries((ROOT / NATIVE / "manifests" / version / "SHA256SUMS").read_bytes(),
                           ORIGINAL_ROOT)
            for area in ("gel", "microarray"):
                for name in ("gentle.stderr.txt", "gentle.stdout.txt", "input.sha256",
                             "openbox.log", "xvfb.log"):
                    relative = Path("native-x11") / area / name
                    with self.subTest(version=version, path=relative):
                        self.assertEqual(hashlib.sha256((ROOT / NATIVE / relative).read_bytes())
                                         .hexdigest(), full[relative])
                inputs = (ROOT / NATIVE / "native-x11" / area / "input.sha256").read_text()
                rows = inputs.splitlines()
                self.assertEqual(len(rows), 2)
                self.assertEqual(rows[0].split("  ", 1)[0], rows[1].split("  ", 1)[0])

    def test_recorded_tests_execute_without_claiming_child_stack_proof(self):
        summary = json.loads((ROOT / GLEN / "SUMMARY.json").read_bytes())
        self.assertEqual(summary["source"], {
            "sha": SOURCE, "tree": "beefc526b11d4801fc57e2189bf0bc1ec18a439d",
            "cargo_lock_sha256": "5893dd19e359d71a9a52ccaf688f952b2ef0507eda53dbdf2c3fdd60a775f84d",
        })
        self.assertEqual(summary["toolchain"]["rustc"], "1.96.1")
        self.assertIsNone(summary["toolchain"]["rust_min_stack"])
        total = 0
        for name, expected in FILTERS.items():
            stem = ROOT / GLEN / "rust-boundary" / name
            command = Path(str(stem) + ".command.txt").read_text()
            self.assertIn("env -u RUST_MIN_STACK", command)
            self.assertNotIn("default_stack_test_runner", command)
            self.assertEqual(Path(str(stem) + ".exit-code.txt").read_text().strip(), "0")
            stdout = Path(str(stem) + ".stdout.txt").read_text()
            results = re.findall(r"test result: ok\. (\d+) passed; 0 failed; 0 ignored;", stdout)
            self.assertEqual(results, [str(expected)])
            self.assertEqual(len(re.findall(r"^test .+ \.\.\. ok$", stdout, re.M)), expected)
            total += expected
        self.assertEqual(total, summary["rust_boundary"]["passed"])
        self.assertEqual(total, 28)
        dbsnp = (ROOT / GLEN / "rust-boundary/dbsnp_.stdout.txt").read_text()
        for name in ("dbsnp_spdi_genbank_round_trip_preserves_wrapped_json",
                     "dbsnp_spdi_unfolding_preserves_valid_json_invalid_data_and_other_qualifiers"):
            self.assertIn(name + " ... ok", dbsnp)
        self.assertEqual((ROOT / GLEN / "rust-boundary/cargo_check.exit-code.txt")
                         .read_text().strip(), "0")
        self.assertFalse(any(summary["mutations"].values()))

    def test_native_boundaries_and_mixed_newlines_are_preserved(self):
        summary = json.loads((ROOT / GLEN / "SUMMARY.json").read_bytes())
        screenshots = summary["screenshots"]
        self.assertEqual(screenshots["flna"], "pass_virtual_x11")
        self.assertEqual(screenshots["microarray_07_04"], "blocked_gui_unresponsive")
        self.assertEqual(screenshots["gel_10_03"], "blocked_gui_unresponsive")
        self.assertEqual(screenshots["tss_08_16"], "blocked_report_not_adopted")
        for area in ("gel", "microarray"):
            logs = ROOT / NATIVE / "native-x11" / area
            stderr = (logs / "gentle.stderr.txt").read_bytes()
            self.assertIn(b"DRI3", stderr)
            self.assertNotIn(b"panicked", stderr)
            self.assertEqual((logs / "gentle.stdout.txt").read_bytes(), b"")
            openbox = (logs / "openbox.log").read_bytes()
            self.assertIn(b"\r\n", openbox)
            self.assertIn(b"\n", openbox.replace(b"\r\n", b""))

    def test_manifest_paths_cannot_escape_or_alias(self):
        digest = "a" * 64
        for name in ("../outside", "/absolute", "x/../outside", "x//y", "x/./y",
                     "C:escape", "x\\..\\outside"):
            with self.subTest(name=name), self.assertRaises(ValueError):
                entries(f"{digest}  ./{name}\n".encode(), "./")
        with self.assertRaises(ValueError):
            entries(f"{digest}  ./same\n{digest}  ./same\n".encode(), "./")


if __name__ == "__main__":
    unittest.main()
