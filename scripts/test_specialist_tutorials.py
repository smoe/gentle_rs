"""Offline regression companions for the four existing 09/10 specialist guides.

Run: python3 scripts/test_specialist_tutorials.py
Optional real CLI: GENTLE_TUTORIAL_BIN_DIR=/absolute/bin/dir python3 \
    scripts/test_specialist_tutorials.py

Never builds Cargo, regenerates catalogs, contacts vendors, invokes agents, or
claims live GUI acceptance. Default checks inspect source units, links, shell
syntax, UI labels and existing fixture bytes, not Rust runtime behavior. The
opt-in tests use only local service/trace operations in isolated temporary
directories, including catalog-root overrides. A missing requested binary is
an error, not a silent skip. Temporary output is removed on test completion.

Fixture provenance: test_files/fixtures/sequencing_confirmation/README.md pins
the Biopython 3100.ab1/fake.ab1 sources. The tutorial input README documents the
trace-derived FASTA pair. Vendor requests are the existing synthetic
DEMO_DO_NOT_ORDER JSONs in docs/examples/external_services. No new biological
fixture is introduced; an optional renamed copy tests magic-byte detection.
"""

import copy
import csv
from datetime import date
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import struct
import subprocess
import tempfile
import unittest
from urllib.parse import unquote, urlsplit


ROOT = Path(__file__).resolve().parents[1]
TUTORIALS = ROOT / "docs/tutorial"
FIXTURES = ROOT / "test_files/fixtures/sequencing_confirmation"
BIN_DIR = os.environ.get("GENTLE_TUTORIAL_BIN_DIR")
STEMS = (
    "09-01_metabion_external_service_handoff_gui_cli",
    "09-02_geneart_external_service_handoff_gui_cli",
    "10-01_sequencing_confirmation_trace_cli",
    "10-02_sequencing_confirmation_gui",
)
REQUESTS = (
    ("metabion_oligo_single_tube_request.json", "metabion", "dna_oligo_single_tube"),
    ("metabion_mblock_request.json", "metabion", "dna_fragment"),
    ("geneart_cloned_gene_request.json", "geneart", "cloned_gene"),
    ("geneart_protein_expression_request.json", "geneart", "protein_expression"),
)


def read_json(path):
    return json.loads(path.read_text(encoding="utf-8"))


def fasta_bases(name):
    return "".join(line.strip() for line in (TUTORIALS / "inputs" / name).read_text(
        encoding="ascii").splitlines() if not line.startswith(">"))


def comparable_sequences(sequences):
    """Retain all values; ignore only unordered restriction-cache iteration."""
    result = copy.deepcopy(sequences)
    for record in result.values():
        # DNAsequence computes sites with Rayon flat_map and serializes its
        # group HashMap as pairs. Neither collection has a stable row order.
        for field in ("restriction_enzyme_sites", "restriction_enzyme_groups"):
            if field in record:
                record[field].sort(key=lambda row: json.dumps(row, sort_keys=True))
    return result


def fixture_abif_tag(tag, number=2):
    """Read a pinned fixture's byte-array tag, not a general trace importer."""
    raw = (FIXTURES / "3100.ab1").read_bytes()
    count = struct.unpack_from(">I", raw, 18)[0]
    offset = struct.unpack_from(">I", raw, 26)[0]
    for index in range(count):
        entry_offset = offset + index * 28
        entry = struct.unpack_from(">4sIHHIIII", raw, entry_offset)
        if entry[:2] == (tag, number):
            size = entry[5]
            start = entry_offset + 20 if size <= 4 else entry[6]
            data = raw[start:start + size]
            if len(data) != size:
                raise ValueError("Truncated fixture tag")
            return data
    raise ValueError(f"Missing fixture tag {tag!r}{number}")


class SpecialistSourceTests(unittest.TestCase):
    def test_sequence_comparison_preserves_bases_features_and_cut_coordinates(self):
        sequence = {"demo": {"seq": {"seq": [65, 67], "features": [{"location": [0, 2]}]},
                             "restriction_enzyme_sites": [{"pos": 1}, {"pos": 0}],
                             "restriction_enzyme_groups": [["b", ["B"]], ["a", ["A"]]]}}
        reordered = copy.deepcopy(sequence)
        for field in ("restriction_enzyme_sites", "restriction_enzyme_groups"):
            reordered["demo"][field].reverse()
        self.assertEqual(comparable_sequences(sequence), comparable_sequences(reordered))
        for field in ("base", "annotation", "cut"):
            changed = copy.deepcopy(reordered)
            if field == "base":
                changed["demo"]["seq"]["seq"][0] = 71
            elif field == "annotation":
                changed["demo"]["seq"]["features"][0]["location"] = [1, 2]
            else:
                changed["demo"]["restriction_enzyme_sites"][0]["pos"] = 2
            self.assertNotEqual(comparable_sequences(sequence), comparable_sequences(changed))
        self.assertEqual(sequence["demo"]["restriction_enzyme_sites"][0]["pos"], 1)

    def test_source_schema_identity_and_discovery_notes(self):
        rust = (ROOT / "src/workflow_examples.rs").read_text(encoding="utf-8")
        # Match current Rust fields so invented metadata cannot silently vanish
        # during the parent's later catalog generation.
        def fields(name):
            body = re.search(rf"pub struct {name} \{{(.*?)\n\}}", rust, re.S).group(1)
            return set(re.findall(r"pub (\w+):", body))

        source_fields = fields("TutorialSourceUnit")
        catalog_fields = fields("TutorialSourceCatalogSection")
        catalog_fields.remove("entry_type")
        catalog_fields.add("type")
        for stem in STEMS:
            with self.subTest(stem=stem):
                source = read_json(TUTORIALS / "sources" / f"{stem}.json")
                self.assertLessEqual(set(source), source_fields)
                self.assertEqual(source["schema"], "gentle.tutorial_source.v4")
                self.assertEqual(source["id"], stem[6:])
                self.assertIsNone(source.get("generated_chapter"))
                catalog = source["catalog"]
                self.assertLessEqual(set(catalog), catalog_fields)
                self.assertEqual(catalog["path"], f"docs/tutorial/{stem}.md")
                self.assertEqual(catalog["group"], stem[:2])
                self.assertEqual(catalog["group_position"], int(stem[3:5]))
                self.assertIs(type(catalog["order"]), int)
                self.assertEqual(catalog["status"], "manual/hybrid")
                self.assertEqual(catalog["source"], "hand_written_markdown")
                self.assertIn(catalog["type"], (
                    "gui_cli_walkthrough", "gui_walkthrough", "cli_walkthrough"))
                self.assertIn("agent_users", catalog["audiences"])
                guide = (ROOT / catalog["path"]).read_text(encoding="utf-8")
                self.assertEqual(guide.splitlines()[0], "# " + source["title"])
                notes = catalog["notes"]
                self.assertIn("Offline", notes)
                self.assertIn(f"ui open tutorial-guide {source['id']}", notes)
                reviewed = re.search(r"codex source review (\d{4}-\d{2}-\d{2})", notes)
                self.assertIsNotNone(reviewed)
                date.fromisoformat(reviewed.group(1))
                self.assertIn("32875224", notes)
                self.assertIn("not", notes)
                self.assertIn("live GUI acceptance", notes)
                self.assertIn("human scientific approval", notes)

    def test_relative_links_and_literal_command_inputs_exist(self):
        for stem in STEMS:
            guide = TUTORIALS / f"{stem}.md"
            text = guide.read_text(encoding="utf-8")
            links = re.findall(r"\]\(([^)]+)\)", text)
            self.assertTrue(links, stem)
            for link in links:
                with self.subTest(stem=stem, link=link):
                    target = urlsplit(link)
                    if target.scheme or target.netloc:
                        continue
                    path = (guide.parent / unquote(target.path)).resolve()
                    self.assertTrue(path.is_relative_to(ROOT), link)
                    self.assertTrue(path.is_file(), link)
                    if target.fragment:
                        headings = re.findall(r"^#+ (.+)$", path.read_text(
                            encoding="utf-8"), re.M)
                        anchors = {re.sub(r"[^\w\- ]", "", h.lower()).replace(" ", "-")
                                   for h in headings}
                        self.assertIn(unquote(target.fragment), anchors, link)
            for path in re.findall(r'@((?:docs|test_files|assets)/[^\s\'"`]+)', text):
                self.assertTrue((ROOT / path).is_file(), path)
            for path in re.findall(r'"path":"([^"]+)"', text):
                self.assertTrue((ROOT / path).is_file(), path)

    def test_bash_examples_are_safe_and_syntactically_valid(self):
        for stem in STEMS[:3]:
            text = (TUTORIALS / f"{stem}.md").read_text(encoding="utf-8")
            blocks = re.findall(r"```bash\n(.*?)```", text, re.S)
            self.assertTrue(blocks)
            commands = "\n".join(blocks)
            result = subprocess.run(["bash", "-n"], input=commands, text=True,
                                    capture_output=True, timeout=10)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn("mktemp -d", commands)
            self.assertIn('STATE="$RUN_DIR/', commands)
            self.assertNotRegex(commands, r"\b(?:rm|cargo|curl|wget)\b")
            for line in commands.splitlines():
                if line.startswith('"$GENTLE_CLI"'):
                    self.assertIn('--state "$STATE"', line)
                if "--output-dir" in line:
                    self.assertIn('"$RUN_DIR/', line)

    def test_vendor_requests_and_boundaries(self):
        for name, provider, kind in REQUESTS:
            request = read_json(ROOT / "docs/examples/external_services" / name)
            self.assertEqual(request["schema"], "gentle.external_service_request.v1")
            self.assertEqual((request["provider"], request["service_kind"]), (provider, kind))
            self.assertIn("DEMO_DO_NOT_ORDER", request["source_target"]["name"])
            self.assertTrue(request["request_metadata"]["human_review_only"])
            self.assertTrue(request["return_spec"]["redact_commercial_fields"])
        for stem in STEMS[:2]:
            text = (TUTORIALS / f"{stem}.md").read_text(encoding="utf-8")
            for term in ("offline", "no network", "handoff_ready", "direct_submission_available",
                         "false", "Provider Config Doctor", "does not expand shell variables"):
                self.assertIn(term, text)
        geneart = (TUTORIALS / f"{STEMS[1]}.md").read_text(encoding="utf-8")
        self.assertIn("does not guarantee", geneart)

    def test_fixture_provenance_and_low_confidence_checkpoint(self):
        for name, digest in (
            ("3100.ab1", "78588d824dd967a1e58a9667a8f7ca1687ba38266420bbc5a6562c1893888efd"),
            ("fake.ab1", "86d9ee4daee09fdcd74f8bb53b43cf5a3ed48633816f6d25892c193d99caaea6"),
        ):
            self.assertEqual(hashlib.sha256((FIXTURES / name).read_bytes()).hexdigest(), digest)
            self.assertIn(name, (FIXTURES / "README.md").read_text(encoding="utf-8"))
        expected = fasta_bases("sequencing_confirmation_trace_demo_construct.fa")
        baseline = fasta_bases("sequencing_confirmation_trace_demo_baseline.fa")
        self.assertEqual(len(expected), 48)
        self.assertEqual(len(baseline), 48)
        self.assertEqual(expected, fixture_abif_tag(b"PBAS")[:48].decode("ascii"))
        self.assertEqual([i for i, pair in enumerate(zip(expected, baseline))
                          if pair[0] != pair[1]], [23])
        self.assertEqual((expected[23], baseline[23]), ("A", "G"))
        self.assertEqual(fixture_abif_tag(b"PCON")[23], 9)
        self.assertFalse(list(FIXTURES.glob("*.scf")), "Document new SCF fixture coverage")
        for stem in STEMS[2:]:
            text = (TUTORIALS / f"{stem}.md").read_text(encoding="utf-8")
            for term in ("SCF", "magic bytes", "synthetic", "quality", "low_confidence_or_ambiguous",
                         "independent biological validation"):
                self.assertIn(term, text)

    def test_reviewed_ui_names_and_import_source_contract(self):
        vendor_ui = (ROOT / "src/app/external_services_ui.rs").read_text(encoding="utf-8")
        trace_ui = (ROOT / "src/main_area_dna/sequencing_confirmation_ui.rs").read_text(
            encoding="utf-8")
        for source, labels in (
            (vendor_ui, ("Refresh Providers", "Provider Config Doctor", "Use Selected Template",
                         "Preflight", "Prepare Quote Handoff", "Output dir",
                         "Export Handoff Bundle", "Bundle files")),
            (trace_ui, ("Baseline/reference sequence ID", "Read sequence IDs", "Imported trace IDs",
                        "Raw Trace Import", "Browse...", "Import trace", "Run confirmation",
                        "Junction breakpoints (0-based)", "Refresh reports", "Show selected",
                        "Export JSON...", "Export TSV...", "Copy summary", "Export summary...")),
        ):
            for label in labels:
                self.assertIn(json.dumps(label), source)
        labels = read_json(ROOT / "assets/i18n/en-US.json")
        gui = (TUTORIALS / f"{STEMS[3]}.md").read_text(encoding="utf-8")
        for key in ("menu.file.new_project", "menu.file.open_sequence", "menu.file.save_project",
                    "menu.patterns.sequencing_confirmation"):
            self.assertIn(labels[key], gui)
        self.assertNotIn("Save Project As...", gui)
        parser = (ROOT / "src/engine/io/sequencing_traces.rs").read_text(encoding="utf-8")
        for token in ('b"ABIF"', 'b".scf"', '"PCON"', '"PLOC"', "scf_called_base_confidence",
                      "parse_scf_trace_extracts_called_bases_and_comment_metadata"):
            self.assertIn(token, parser)
        engine = (ROOT / "src/engine/analysis/sequencing_confirmation.rs").read_text(
            encoding="utf-8")
        self.assertRegex(engine, r"fn sequencing_trace_variant_low_confidence_threshold\(\) -> u8"
                         r"\s*\{\s*20\s*\}")


@unittest.skipUnless(BIN_DIR, "Set GENTLE_TUTORIAL_BIN_DIR for real local CLI smoke; no build performed")
class SpecialistBinaryTests(unittest.TestCase):
    def setUp(self):
        suffix = ".exe" if os.name == "nt" else ""
        self.cli = Path(BIN_DIR).expanduser().resolve() / f"gentle_cli{suffix}"
        self.assertTrue(self.cli.is_file() and os.access(self.cli, os.X_OK), str(self.cli))
        temp = tempfile.TemporaryDirectory(prefix="gentle specialist ")
        self.addCleanup(temp.cleanup)
        self.run_dir = Path(temp.name).resolve()
        self.state = self.run_dir / "tutorial.gentle.json"
        self.env = dict(os.environ, GENTLE_ASSET_ROOT=str(ROOT),
                        GENTLE_SYSTEM_CONFIG_ROOT=str(self.run_dir / "system"),
                        XDG_CONFIG_HOME=str(self.run_dir / "user"),
                        GENTLE_PROJECT_ROOT=str(self.run_dir / "project"))

    def invoke(self, *args, succeeds=True):
        result = subprocess.run([str(self.cli), "--state", str(self.state), "--full-report", *args],
                                cwd=ROOT, env=self.env, text=True, capture_output=True, timeout=90)
        if succeeds:
            self.assertEqual(result.returncode, 0, result.stderr + result.stdout)
            return json.loads(result.stdout)
        self.assertNotEqual(result.returncode, 0, result.stdout)
        return result

    def shell(self, *tokens, **kwargs):
        return self.invoke("shell", shlex.join(str(token) for token in tokens), **kwargs)

    def test_vendor_preflight_and_export_all_four_requests(self):
        doctor = self.invoke("services", "providers", "doctor")
        self.assertEqual(doctor["schema"], "gentle.external_service_provider_config_doctor.v1")
        self.assertEqual(doctor["error_count"], 0)
        catalog = self.invoke("services", "providers", "list")
        self.assertEqual(catalog["schema"], "gentle.external_service_provider_catalog.v1")
        for name, provider, kind in REQUESTS:
            with self.subTest(request=name):
                request = ROOT / "docs/examples/external_services" / name
                preflight = self.invoke("services", "project-preflight", f"@{request}")
                self.assertEqual(preflight["schema"], "gentle.external_service_preflight.v1")
                self.assertEqual((preflight["provider"], preflight["service_kind"]), (provider, kind))
                self.assertTrue(preflight["eligible"])
                self.assertTrue(preflight["quote_handoff_available"])
                self.assertFalse(preflight["direct_submission_available"])
                output = self.run_dir / request.stem
                self.assertFalse(output.exists())
                quote = self.invoke("services", "project-quote", f"@{request}",
                                    "--output-dir", str(output))
                self.assertEqual(quote["schema"], "gentle.external_service_quote.v1")
                self.assertEqual(quote["quote_status"], "handoff_ready")
                self.assertFalse(quote["preflight"]["direct_submission_available"])
                self.assertEqual(read_json(output / "quote_report.json"), quote)
                bundle = quote["service_ready_bundle"]
                self.assertEqual(bundle["schema"], "gentle.external_service_artifact_bundle.v1")
                payloads = {p["payload_kind"]: p for p in bundle["inline_payloads"]}
                self.assertLessEqual({"handoff_markdown", "redacted_request_json",
                                      "normalized_line_items_json", "normalized_line_items_csv",
                                      "email_draft_markdown", "guided_wop_checklist"}, set(payloads))
                for artifact in bundle["local_files"]:
                    path = Path(artifact["path"]).resolve()
                    self.assertTrue(path.is_relative_to(output), path)
                    self.assertTrue(path.is_file(), path)
                    if artifact["artifact_kind"] in payloads:
                        self.assertEqual(path.read_text(encoding="utf-8"),
                                         payloads[artifact["artifact_kind"]]["text"])
        invalid = read_json(ROOT / "docs/examples/external_services" / REQUESTS[0][0])
        del invalid["source_target"]["line_items"]
        blocked = self.invoke("services", "project-quote", json.dumps(invalid))
        self.assertEqual(blocked["quote_status"], "blocked")
        self.assertFalse(blocked["preflight"]["eligible"])
        self.assertTrue(blocked["preflight"]["blocking_issues"])

    def load_trace_inputs(self, baseline=False):
        self.expected_id = ("sequencing_confirmation_trace_demo_construct" if baseline
                            else "trace_demo_construct")
        self.trace_id = "abi_demo_trace_gui" if baseline else "abi_demo_trace"
        for suffix in (("construct", "baseline") if baseline else ("construct",)):
            load = {"path": str(TUTORIALS / "inputs" /
                                f"sequencing_confirmation_trace_demo_{suffix}.fa")}
            if not baseline:
                load["as_id"] = self.expected_id
            self.invoke("op", json.dumps({"LoadFile": load}), "--confirm")
        sequences = read_json(self.state)["sequences"]
        self.assertEqual(set(sequences), {self.expected_id, "sequencing_confirmation_trace_demo_baseline"}
                         if baseline else {self.expected_id})
        imported = self.shell("seq-trace", "import", FIXTURES / "3100.ab1",
                              "--trace-id", self.trace_id, "--seq-id", self.expected_id)
        self.assertEqual(comparable_sequences(read_json(self.state)["sequences"]),
                         comparable_sequences(sequences))
        self.assertEqual(imported["import_report"]["format"], "abi_ab1")
        trace = self.shell("seq-trace", "show", self.trace_id)["trace"]
        self.assertEqual(trace["called_bases"][:48],
                         fasta_bases("sequencing_confirmation_trace_demo_construct.fa"))
        self.assertEqual(trace["called_base_confidence_values"][23], 9)
        self.assertEqual(len(trace["peak_locations"]), len(trace["called_bases"]))
        self.assertEqual(len(trace["channel_data"]), 4)
        self.assertTrue(all(channel["points"] for channel in trace["channel_data"]))

    def test_trace_junction_persistence_exports_and_bad_import(self):
        self.load_trace_inputs()
        report = self.shell("seq-confirm", "run", "trace_demo_construct", "--trace-id",
                            "abi_demo_trace", "--junction", "24", "--junction-flank", "12",
                            "--report-id", "trace_demo_confirm")["report"]
        self.assertEqual(report["overall_status"], "confirmed")
        self.assertEqual(report["read_seq_ids"], [])
        self.assertEqual(report["trace_ids"], ["abi_demo_trace"])
        self.assertEqual(len(report["targets"]), 1)
        target = report["targets"][0]
        self.assertEqual(target["target_id"], "junction_1")
        self.assertEqual((target["start_0based"], target["end_0based_exclusive"]), (12, 36))
        self.assertEqual(target["support_read_ids"], ["abi_demo_trace"])
        self.assertEqual(report["reads"][0]["evidence_kind"], "trace")
        self.assertEqual(self.shell("seq-confirm", "show-report", "trace_demo_confirm")["report"], report)
        before = self.state.read_bytes()
        output = self.run_dir / "trace_demo_confirm.json"
        self.shell("seq-confirm", "export-report", "trace_demo_confirm", output)
        self.assertEqual(read_json(output), report)
        tsv = self.run_dir / "trace_demo_confirm.tsv"
        self.shell("seq-confirm", "export-support-tsv", "trace_demo_confirm", tsv)
        with tsv.open(encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["target_id"], "junction_1")
        self.assertEqual(rows[0]["support_read_ids"], "abi_demo_trace")
        self.assertEqual(self.state.read_bytes(), before)
        self.shell("seq-trace", "import", FIXTURES / "fake.ab1", "--trace-id", "fake_trace",
                   succeeds=False)
        self.assertEqual(self.state.read_bytes(), before)
        self.assertEqual(self.shell("seq-trace", "list")["trace_count"], 1)
        renamed = self.run_dir / "still_abif.scf"
        shutil.copyfile(FIXTURES / "3100.ab1", renamed)
        imported = self.shell("seq-trace", "import", renamed, "--trace-id", "renamed_abi")
        self.assertEqual(imported["import_report"]["format"], "abi_ab1")

    def test_gui_equivalent_targets_keep_low_confidence_edit_unresolved(self):
        # Exercise the GUI's shared operation, not a fabricated GUI success claim.
        self.load_trace_inputs(baseline=True)
        self.invoke("op", json.dumps({"ConfirmConstructReads": {
            "expected_seq_id": self.expected_id,
            "baseline_seq_id": "sequencing_confirmation_trace_demo_baseline",
            "trace_ids": [self.trace_id], "read_seq_ids": [],
            "targets": [
                {"target_id": "full_span", "label": "Full construct span", "kind": "full_span",
                 "start_0based": 0, "end_0based_exclusive": 48, "required": True},
                {"target_id": "junction_1", "label": "Junction @ 24", "kind": "junction",
                 "start_0based": 12, "end_0based_exclusive": 36,
                 "junction_left_end_0based": 24, "required": True}],
            "alignment_mode": "local", "match_score": 2, "mismatch_score": -3,
            "gap_open": -5, "gap_extend": -1, "min_identity_fraction": 0.80,
            "min_target_coverage_fraction": 1.0, "allow_reverse_complement": True,
            "report_id": "trace_demo_gui_confirm"}}), "--confirm")
        report = self.shell("seq-confirm", "show-report", "trace_demo_gui_confirm")["report"]
        self.assertEqual(report["schema"], "gentle.sequencing_confirmation_report.v1")
        self.assertEqual(report["overall_status"], "insufficient_evidence")
        self.assertEqual(len(report["targets"]), 3)
        for target in report["targets"]:
            self.assertEqual(target["status"], "insufficient_evidence"
                             if target["kind"] == "expected_edit" else "confirmed")
        self.assertEqual(len(report["variants"]), 1)
        variant = report["variants"][0]
        self.assertEqual(variant["classification"], "low_confidence_or_ambiguous")
        self.assertEqual(variant["confidence_min"], 9)
        self.assertEqual(variant["observed_bases"], "A")
        self.assertEqual(variant["baseline_bases"], "G")
        self.assertEqual(variant["expected_bases"], "A")


if __name__ == "__main__":
    unittest.main(verbosity=2)
