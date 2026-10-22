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
import sys
import tempfile
import unittest
import xml.etree.ElementTree as ET
from urllib.parse import unquote, urlsplit


ROOT = Path(__file__).resolve().parents[1]
TUTORIALS = ROOT / "docs/tutorial"
FIXTURES = ROOT / "test_files/fixtures/sequencing_confirmation"
BIN_DIR = os.environ.get("GENTLE_TUTORIAL_BIN_DIR")
PLANNING_GUIDE = TUTORIALS / "09-04_cloning_time_cost_and_protocol_suppliers.md"
CATALOG_EXERCISE = ROOT / "docs/examples/planning/supplier_catalog_update_exercise.md"
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


def tutorial_python(text, invocation):
    """Extract one literal teaching snippet; never execute its shell commands."""
    match = re.search(re.escape(invocation) + r" <<'PY'\n(.*?)\nPY", text.replace("\r\n", "\n"), re.S)
    if match is None:
        raise ValueError(f"Missing tutorial Python block: {invocation}")
    return match.group(1)


def run_tutorial_python(text, invocation, *args):
    result = subprocess.run([sys.executable, "-", *(str(arg) for arg in args)],
                            input=tutorial_python(text, invocation), text=True,
                            capture_output=True, timeout=15, cwd=ROOT)
    if result.returncode:
        raise AssertionError(result.stdout + result.stderr)
    return result.stdout


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
    def test_parity_overrides_only_exclude_unreachable_surfaces(self):
        overrides = read_json(ROOT / "docs/parity_matrix_overrides.json")["overrides"]
        reachable = [
            f"{row['source']} {row['name']} / {row['adapter']}: {row['surfacing']}"
            for row in overrides if row["surfacing"] != "not_applicable"
        ]
        self.assertEqual([], reachable, "Reachability belongs in the canonical registry, not n/a overrides")
        for row in overrides:
            self.assertTrue(row["reason"].strip())

    def test_guarded_provider_import_guide_and_catalog_surface(self):
        guide = PLANNING_GUIDE.read_text(encoding="utf-8") + CATALOG_EXERCISE.read_text(encoding="utf-8")
        contract = (ROOT / "docs/catalog_updates.md").read_text(encoding="utf-8")
        glossary = read_json(ROOT / "docs/glossary.json")
        commands = {row["path"]: row for row in glossary["commands"]}
        overrides = read_json(ROOT / "docs/parity_matrix_overrides.json")["overrides"]
        matrix = (ROOT / "docs/gui_cli_mcp_parity.md").read_text(encoding="utf-8")
        for path in ("catalog doctor", "catalog search", "catalog updates preview", "catalog updates apply"):
            self.assertIn(path, guide)
            self.assertIn(path, contract)
            self.assertIn(path, commands)
            self.assertIn(f"| {path} | glossary-command | shell-only | prominent | shell-only | shell-only | shell-only | n/a |", matrix)
            for adapter in ("mcp", "js", "lua"):
                self.assertFalse(any(item["name"] == path and item["adapter"] == adapter
                                     for item in overrides))
        for term in ("withdrawn_providers", "previous.bin", "freshness", "not signatures", "private", "atomic", "Glen/CI"):
            self.assertIn(term, contract)
        # Parse/execute only the tutorial's Python data construction, never its CLI.
        match = re.search(r'python3 - "\$RUN/provider-update.json" <<\'PY\'\n(.*?)\nPY', guide, re.S)
        self.assertIsNotNone(match)
        with tempfile.TemporaryDirectory() as temp:
            bundle_path = Path(temp) / "bundle with spaces.json"
            subprocess.run(["python3", "-", str(bundle_path)], input=match.group(1), text=True, check=True)
            bundle = read_json(bundle_path)
            self.assertEqual(bundle["schema"], "gentle.provider_update_bundle.v1")
            self.assertIsNone(bundle["changes"][0]["base_record_sha256"])
            record = bundle["changes"][0]["record"]
            self.assertEqual(record["provider"], "synthetic_tutorial_supplier")
            self.assertEqual(record["support_status"], "catalog_only")
            self.assertFalse(record["capabilities"][0]["quote_handoff_supported"])
            self.assertEqual(record["website_url"], "https://example.invalid")
            protocol = (ROOT / "crates/gentle-protocol/src/lib.rs").read_text(encoding="utf-8")
            row_type = protocol.split("pub struct ExternalServiceProviderConfigRecord {", 1)[1].split("\n}", 1)[0]
            fields = set(re.findall(r"pub (\w+):", row_type))
            self.assertEqual(set(record), fields - {"country_code", "country_source_url", "protocols"})
            capability_type = protocol.split("pub struct ExternalServiceCapability {", 1)[1].split("\n}", 1)[0]
            self.assertEqual(set(record["capabilities"][0]), set(re.findall(r"pub (\w+):", capability_type)))
        shared = (ROOT / "src/external_service_providers/updates.rs").read_text(encoding="utf-8")
        self.assertIn("try_lock()", shared)
        self.assertIn("stale_or_tampered_preview", shared)
        self.assertIn("persist(&context.target)", shared)
        self.assertNotIn("reqwest::", shared)

    def test_protocol_overview_guide_distinguishes_planned_dependencies_from_execution(self):
        guide = (TUTORIALS / "09-04_cloning_time_cost_and_protocol_suppliers.md").read_text(encoding="utf-8")
        for term in ("Protocol Overview", "seven synthetic bench steps", "not an elapsed-time axis",
                     "planning.protocol_overview", "gentle.planning_protocol_overview.v1",
                     "Copy SVG", "Copy JSON", "source_estimate_sha256", "content_sha256",
                     "Native acceptance remains pending", "No screenshots are fabricated"):
            self.assertIn(term,guide)
        shared = (ROOT / "src/planning_protocol.rs").read_text(encoding="utf-8")
        ui = (ROOT / "src/app/planning_protocol_ui.rs").read_text(encoding="utf-8")
        self.assertIn('"gentle.planning_protocol_overview.v1"',shared)
        self.assertIn("from_explanation(output)",ui)
        self.assertIn("sha256_prefixed_bytes(&bytes)",ui)
        self.assertIn("routine.protocol.copy_svg",ui)
        self.assertNotIn("engine.apply(",ui)

    def test_protocol_overview_semantic_ids_are_static_and_subjects_are_pseudonymous(self):
        source = (ROOT / "src/app/planning_protocol_ui.rs").read_text(encoding="utf-8")
        for text in (source, source.replace("\n", "\r\n")):
            self.assertRegex(text, r"fn semantic\([^\n]*id: &'static str")
            helper = text.split("fn semantic(", 1)[1].split("pub(super) fn show(", 1)[0]
            self.assertIn("pseudonymous_subject_scope(&[id, subject])", helper)
            self.assertIn("scope.as_deref()", helper)
            self.assertRegex(text, r'&response,\s*"routine\.protocol\.step",\s*Some\(&node\.step\.step_id\)')
            self.assertNotIn('format!("routine.protocol.step.', text)
            self.assertIn("fn protocol_overview_semantic_steps_have_static_ids_and_distinct_scopes()", text)

    def test_cloning_planning_guide_and_protocol_supplier_scope(self):
        stem = "09-04_cloning_time_cost_and_protocol_suppliers"
        source = read_json(TUTORIALS / "sources" / f"{stem}.json")
        self.assertEqual(source["catalog"]["group_position"], 4)
        self.assertEqual(source["catalog"]["status"], "manual/hybrid")
        self.assertNotIn("generated_chapter", source)
        projection = next(row for row in read_json(TUTORIALS / "catalog.json")["entries"]
                          if row["id"] == source["id"])
        self.assertEqual(projection["notes"], source["catalog"]["notes"])
        review = next(row for row in read_json(TUTORIALS / "review_manifest.json")["entries"]
                      if row["tutorial_id"] == source["id"])
        self.assertEqual(projection["codex_reviewed_at"], review["codex_reviewed_at"])
        guide_path = TUTORIALS / f"{stem}.md"
        guide = guide_path.read_text(encoding="utf-8")
        for newline in ("\n", "\r\n"):
            text = guide.replace("\n", newline)
            for term in ("0.8 h", "2.6 h", "3.2 h", "87.2 h", "39 EUR",
                         "country_code", "text_pdf", "catalog_only", "explanation.time_policy",
                         "non-exhaustive", "No screenshots are fabricated",
                         "not Promega prices", "legacy_family_scalar"):
                self.assertIn(term, text)
            for target in re.findall(r"\]\(([^)]+)\)", text):
                if not urlsplit(target).scheme:
                    self.assertTrue((guide_path.parent / target.split("#")[0]).is_file(), target)
        catalog = read_json(ROOT / "assets/external_service_providers.d/protocol_suppliers.json")
        self.assertEqual(len(catalog["providers"]), 12)
        for provider in catalog["providers"]:
            self.assertRegex(provider["country_code"], r"^[A-Z]{2}$")
            self.assertTrue(provider["country_source_url"].startswith("https://"))
            self.assertEqual(provider["support_status"], "catalog_only")
            self.assertTrue(provider["protocols"])
            for protocol in provider["protocols"]:
                self.assertIn(protocol["format"], ("html", "text_pdf"))
                date.fromisoformat(protocol["checked_on"])
                self.assertTrue(protocol["url"].startswith("https://"))
            for capability in provider["capabilities"]:
                self.assertFalse(capability["quote_handoff_supported"])
                self.assertFalse(capability["direct_api_implemented"])
                self.assertEqual(capability["supported_submission_modes"], [])
        for name in ("cloning_in_stock", "cloning_waiting"):
            profile = read_json(ROOT / f"docs/examples/planning/{name}.json")
            self.assertTrue(profile["profile_id"].startswith("DEMO_DO_NOT_ORDER"))
            self.assertIn("Synthetic", profile["inventory"]["ligase"]["note"])

    def test_09_04_reuses_the_canonical_vector_workflow_without_new_operations(self):
        guide = PLANNING_GUIDE.read_text(encoding="utf-8")
        canonical = read_json(ROOT / "docs/examples/workflows/promoter_reporter_panel_planning_offline.json")
        for newline in ("\n", "\r\n"):
            with self.subTest(newline=repr(newline)), tempfile.TemporaryDirectory(prefix="09.04 O'Brien ") as temp:
                run = Path(temp).resolve()
                run_tutorial_python(guide.replace("\n", newline), 'python3 - "$ROOT" "$RUN"', ROOT, run)
                actual = read_json(run / "vector-workflow.json")
                expected = copy.deepcopy(canonical)
                expected["required_files"] = [str(ROOT / p) for p in expected["required_files"]]
                for op in expected["workflow"]["ops"]:
                    kind, data = next(iter(op.items()))
                    if kind == "LoadFile":
                        data["path"] = str(ROOT / data["path"])
                    elif kind == "RenderSequenceSvg":
                        data["path"] = str(run / Path(data["path"]).name)
                    else:
                        self.assertEqual(kind, "PlanPromoterReporterPanel")
                        request = data["request"]
                        request["helper_catalog_path"] = str(ROOT / request["helper_catalog_path"])
                        for member in request["members"]:
                            member["candidate_set_path"] = str(ROOT / member["candidate_set_path"])
                        request["output_dir"] = str(run / "planned-products")
                self.assertEqual(actual, expected)
                self.assertEqual(read_json(run / "vector-request.json"), expected["workflow"]["ops"][-1]["PlanPromoterReporterPanel"]["request"])
                self.assertFalse((run / "planned-products").exists())
                self.assertEqual(sorted(p.name for p in run.iterdir()), ["vector-request.json", "vector-workflow.json"])

    def test_09_04_supplier_table_preserves_unknown_country_and_protocols(self):
        # Hand-crafted presentation inputs, not supplier facts or eligibility evidence.
        guide = PLANNING_GUIDE.read_text(encoding="utf-8")
        with tempfile.TemporaryDirectory() as temp:
            source, output = Path(temp) / "suppliers.json", Path(temp) / "comparison.tsv"
            source.write_text(json.dumps({"providers": [
                {"provider": "unknown", "support_status": "catalog_only"},
                {"provider": "synthetic", "country_code": "DE", "support_status": "catalog_only",
                 "protocols": [{"format": "text_pdf", "title": "Synthetic PDF"}, {"format": "html", "title": "Synthetic HTML"}]},
            ]}), encoding="utf-8")
            run_tutorial_python(guide, 'python3 - "$RUN/suppliers.json" "$RUN/supplier-comparison.tsv"', source, output)
            with output.open(encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual([r["provider"] for r in rows], ["synthetic", "unknown"])
            self.assertEqual(rows[0]["formats"], "html,text_pdf")
            self.assertEqual(rows[1]["country"], "unknown")
            self.assertEqual(rows[1]["protocols"], "unassessed")

    def test_09_04_decision_sheet_copies_report_status_and_raw_bindings(self):
        # Minimal synthetic report projections: test presentation, never certify biology.
        guide = PLANNING_GUIDE.read_text(encoding="utf-8")
        estimate = dict(hands_on_hours=None, unattended_hours=None, elapsed_hours=None,
                        estimated_cost=0, routine_base_cost=0, stock_drawn_material_cost=0,
                        missing_material_purchase_proxy_cost=0, machine_run_cost=0)
        for status in ("verified", "rejected", "unavailable"):
            for newline in ("\n", "\r\n"):
                with self.subTest(status=status, newline=repr(newline)), tempfile.TemporaryDirectory() as temp:
                    run = Path(temp)
                    reports = {
                        "vector-proposal.json": dict(vector_validation=dict(helper_catalog_id="Synthetic only", status=status),
                            vector_sequence_sha256="sha256:synthetic-sequence", proposal_digest="sha256:synthetic-proposal",
                            cloning_strategy=dict(strategy="gibson"), products=[]),
                        "in-stock.json": {"planning": {"estimate": estimate}},
                        "waiting.json": {"planning": {"estimate": estimate}},
                        "promega-context.json": {"objective_summary": {"objective": {"helper_profile_id": "Synthetic commercial context"}}},
                        "promega-preflight.json": dict(capability_status="catalog_only", eligible=False,
                            quote_handoff_available=False, direct_submission_available=False),
                        "suppliers.json": {"providers": []},
                    }
                    for name, value in reports.items():
                        (run / name).write_bytes((json.dumps(value, indent=2) + "\n").replace("\n", newline).encode())
                    run_tutorial_python(guide.replace("\n", newline), 'python3 - "$RUN"', run)
                    sheet = (run / "cloning-decision.md").read_text(encoding="utf-8")
                    self.assertIn("HUMAN REVIEW REQUIRED", sheet)
                    self.assertIn("Exact-vector status (synthetic identity only): " + status, sheet)
                    self.assertIn('"hands_on_hours": null', sheet)
                    self.assertIn('"estimated_cost": 0', sheet)
                    self.assertIn('"eligible": false', sheet)
                    self.assertIn('"unpriced_material_classes": null', sheet)
                    self.assertIn("currently not_evaluated", sheet)
                    for name in reports:
                        self.assertIn(f"{name}: sha256:{hashlib.sha256((run / name).read_bytes()).hexdigest()}", sheet)

    def test_09_04_main_and_advanced_examples_have_valid_shell_and_visual_boundaries(self):
        for path in (PLANNING_GUIDE, CATALOG_EXERCISE):
            text = path.read_text(encoding="utf-8")
            blocks = re.findall(r"```bash\n(.*?)```", text, re.S)
            result = subprocess.run(["bash", "-n"], input="\n".join(blocks), text=True,
                                    capture_output=True, timeout=10)
            self.assertEqual(result.returncode, 0, result.stderr)
            for target in re.findall(r"\]\(([^)]+)\)", text):
                if not urlsplit(target).scheme:
                    self.assertTrue((path.parent / target.split("#")[0]).is_file(), target)
        guide = PLANNING_GUIDE.read_text(encoding="utf-8")
        self.assertIn("not a pGL4.10 map", guide)
        self.assertIn("not screenshot placeholders or proof", guide)
        self.assertIn("no vendor procedure has been approved", guide)
        self.assertNotIn("g promoters panel-materialize", guide)

    def test_expression_four_interface_routes_and_discovery_are_explicit(self):
        transport = (ROOT / "src/mcp_server/protein_expression.rs").read_text(encoding="utf-8")
        dispatcher = (ROOT / "src/mcp_server.rs").read_text(encoding="utf-8")
        registry = (ROOT / "crates/gentle-protocol/src/lib.rs").read_text(encoding="utf-8")
        agent = (ROOT / "src/agent_bridge.rs").read_text(encoding="utf-8")
        for token in ("run_non_mutating_shell_tool", "deny_unknown_fields",
                      '"additionalProperties": false', "value.is_object()",
                      "protein_expression_handoff_mcp_and_inner_agent_match_shared_report_without_writes",
                      "mcp_tool_call_for_capability_surface_tests", "fs::read(&state_path)"):
            self.assertIn(token, transport)
        self.assertIn("items.push(protein_expression::descriptor())", dispatcher)
        self.assertIn("protein_expression::call(default_state_path, &params.arguments)", dispatcher)
        self.assertIn('"protein_expression_handoff" => &["planning protein-expression-handoff"]', dispatcher)
        self.assertIn('"protein_expression_handoff",', registry)
        self.assertIn('"planning protein-expression-handoff",', registry)
        self.assertIn("Protein-expression review:", agent)
        self.assertIn("MCP protein_expression_handoff", agent)
        self.assertIn("no inferred product or default provider", (ROOT / "docs/protocol.md").read_text())
        inspector = (ROOT / "src/app/protein_expression_ui.rs").read_text()
        self.assertIn("fn protein_expression_inner_agent_suggestion_matches_inspector_without_mutation", inspector)
        ledger = (ROOT / "docs/agent_interface.md").read_text()
        exclusions = ledger.split("#### Intentionally MCP-excluded shell commands", 1)[1]
        exclusions = exclusions.split("### 4)", 1)[0]
        self.assertNotIn("`planning protein-expression-handoff`", exclusions)
        guide = (TUTORIALS / "09-03_protein_expression_handoff_cli.md").read_text()
        for text in (guide, guide.replace("\n", "\r\n")):
            for term in ("## Interface Parity", "structuredContent", "unsaved GUI edits",
                         "exclude", "generated_at_unix_ms"):
                self.assertIn(term, text)
            examples = [json.loads(block) for block in re.findall(r"```json\r?\n(.*?)```", text, re.S)]
            request = next(row for row in examples if row.get("name") == "protein_expression_handoff")
            self.assertEqual(request["arguments"]["format"], "text")
            self.assertNotIn("provider", request["arguments"])
            refusal = next(row for row in examples if row.get("provider") == "geneart")
            self.assertIs(refusal["requirements"]["outsourcing"]["allowed"], False)

    def test_expression_screenshot_checkpoints_require_real_scoped_capture(self):
        guide = (TUTORIALS / "09-03_protein_expression_handoff_cli.md").read_text()
        inspector = (ROOT / "src/app/protein_expression_ui.rs").read_text()
        for term in ("manual teaching captures", "raw PNG", "not a scientific oracle",
                     "no full generic-runner", "starter/project hashes",
                     "no retained semantic snapshot or typed GUI report JSON"):
            # The locator paragraph says an observation is not a scientific
            # oracle; keep whitespace-independent prose checks readable.
            self.assertIn(term, " ".join(guide.split()))
        for name in ("missing_product", "neutral_product", "provider_withheld", "provider_projected"):
            self.assertIn(f"`{name}.raw.png`", guide)
        for control in ("window.protein_expression_handoff", "protein_expression.product",
                        "protein_expression.provider", "protein_expression.requirements",
                        "protein_expression.inspect", "protein_expression.status",
                        "protein_expression.provider_binding", "protein_expression.copy_json",
                        "protein_expression.copy_text"):
            self.assertIn(f'"{control}"', inspector)
        self.assertIn("fn protein_expression_capture_targets_bind_the_actual_report_status", inspector)
        self.assertNotIn("there is no dedicated expression tool or GUI", guide)
        self.assertNotIn("natural next GUI slice", guide)

    def test_expression_manual_captures_retain_exact_images_and_honest_limits(self):
        source = read_json(TUTORIALS / "sources/09-03_protein_expression_handoff_cli.json")
        folder = ROOT / "docs/screenshots/tutorial_gui_acceptance/protein_expression_handoff_cli"
        evidence = read_json(folder / "evidence.json")
        self.assertEqual(evidence["schema"], "gentle.manual_gui_tutorial_evidence.v1")
        self.assertEqual(evidence["status"], "captured_not_accepted")
        for flag in ("automated_gui_acceptance", "four_interface_runtime_parity_verified",
                     "inner_agent_invoked", "provider_contact_performed"):
            self.assertIs(evidence[flag], False)
        self.assertEqual(evidence["producer"]["revision"], "9f01098b749c79fb5c9eccfcfb184bb085c63357")
        self.assertIs(evidence["producer"]["parity_followup_in_binary"], False)
        self.assertEqual(len(evidence["producer"]["gentle_gui_sha256"]), 64)
        self.assertEqual(len(source["graphics"]), 4)
        guide = (TUTORIALS / source["catalog"]["path"].split("/")[-1]).read_text()
        catalog = read_json(TUTORIALS / "catalog.json")
        entry = next(row for row in catalog["entries"] if row["id"] == source["id"])
        self.assertEqual(entry["graphics"], source["graphics"])
        for graphic, capture in zip(source["graphics"], evidence["captures"]):
            path = ROOT / graphic["path"]
            self.assertEqual(path, folder / capture["raw_png"])
            self.assertIn(f"/{path.name})", guide)
            self.assertEqual(graphic["capture_date"], evidence["capture_date"])
            self.assertEqual(graphic["illustrates_step"], capture["illustrates_step"])
            raw = path.read_bytes()
            self.assertEqual(hashlib.sha256(raw).hexdigest(), capture["raw_png_sha256"])
            self.assertEqual(raw[:8], b"\x89PNG\r\n\x1a\n")
            self.assertEqual(list(struct.unpack(">II", raw[16:24])), capture["png_dimensions_px"])
            self.assertGreaterEqual(min(capture["png_dimensions_px"]), 1000)
            self.assertEqual(capture["image_transforms"], [])
            self.assertIsNone(capture["semantic_snapshot"])
            self.assertIsNone(capture["retained_report_json"])
        self.assertEqual([row["observed_provider_binding"] for row in evidence["captures"]],
                         ["not_selected", "not_selected", "withheld", "projected"])
        project = evidence["project"]
        raw = (folder / project["path"]).read_bytes()
        self.assertEqual(hashlib.sha256(raw).hexdigest(), project["before_sha256"])
        self.assertEqual(project["before_sha256"], project["after_sha256"])
        record = json.loads(raw)["sequences"]["expression_review_product"]["seq"]
        self.assertEqual(bytes(record["seq"]).decode(), "ATGCACCACCACCACCACCACGCTTAA")
        for item in evidence["inputs"]:
            # Original byte hashes are historical; JSON content binding must
            # survive CRLF checkouts of the unpinned source input documents.
            text = (ROOT / item["path"]).read_text()
            for candidate in (text, text.replace("\n", "\r\n")):
                canonical = json.dumps(json.loads(candidate), sort_keys=True,
                                       separators=(",", ":")).encode()
                self.assertEqual(hashlib.sha256(canonical).hexdigest(), item["canonical_json_sha256"])

    def test_expression_inspector_tutorial_uses_shared_copy_only_report_contract(self):
        guide = (TUTORIALS / "09-03_protein_expression_handoff_cli.md").read_text(encoding="utf-8")
        inspector = (ROOT / "src/app/protein_expression_ui.rs").read_text(encoding="utf-8")
        source = read_json(TUTORIALS / "sources/09-03_protein_expression_handoff_cli.json")
        self.assertEqual(source["catalog"]["type"], "gui_cli_walkthrough")
        for text in (guide, guide.replace("\n", "\r\n")):
            for term in ("Services > Protein Expression Handoff...",
                         "ui open protein-expression-handoff",
                         "ui focus protein-expression-handoff",
                         "ui close protein-expression-handoff",
                         "Copy report JSON", "Copy report text",
                         "outsourcing_disallowed", "undo/redo", "late results",
                         "manual/hybrid", "Rust/runtime/native checks"):
                self.assertIn(term, text)
        for term in ("PlanningProteinExpressionHandoff", "execute_shell_command",
                     "clone_without_history", "binding.matches", "instance_id()",
                     "mutation_revision()", "action.shell_line.clone()"):
            self.assertIn(term, inspector)
        for forbidden in ("ServicesQuoteExport", "ServicesProjectPreflight",
                          "ServicesProvidersList", "commit_detached_execution",
                          "geneart_protein_expression_request.json"):
            self.assertNotIn(forbidden, inspector)
        for case in ("matches_shared_shell_without_live_mutation",
                     "missing_product_and_requirements_stay_absent",
                     "ui_intents_reuse_and_close_without_defaults",
                     "discards_late_results_after_request_or_project_changes",
                     "adopts_current_results_and_invalidates_displayed_results"):
            self.assertIn("fn protein_expression_inspector_" + case, inspector)

    def test_expression_walkthrough_has_a_synthetic_end_to_end_starter(self):
        workflow = read_json(ROOT / "docs/examples/workflows/protein_expression_review_starter.json")
        product = workflow["workflow"]["ops"][0]["CreateSequenceFromText"]
        self.assertEqual(product["output_id"], "expression_review_product")
        self.assertEqual(product["sequence_text"], "ATG" + "CAC" * 6 + "GCTTAA")
        requirements = read_json(ROOT / "docs/examples/plans/protein_expression_requirements_review_only.json")
        self.assertTrue(requirements["outsourcing"]["allowed"])
        self.assertEqual(requirements["outsourcing"]["approved_providers"], ["geneart"])
        self.assertEqual(len(set(requirements) & {"yield_goal", "chassis", "localization", "folding", "toxicity_induction", "tag_policy", "scale_purification", "outsourcing"}), 8)
        guide = (TUTORIALS / "09-03_protein_expression_handoff_cli.md").read_text()
        for term in ("expression_review_product", "expression-review.gentle.json", "manual/hybrid", "agent_execute_plan", "not falsely labelled", "synthetic permissions"):
            self.assertIn(term, guide)
        mcp = (ROOT / "src/mcp_server.rs").read_text()
        self.assertIn("fn protein_expression_mcp_plan_executor_matches_shared_shell_without_writes", mcp)

    def test_expression_draft_guide_matches_explicit_provider_glossary_contract(self):
        stem = "09-03_protein_expression_handoff_cli"
        source = read_json(TUTORIALS / "sources" / f"{stem}.json")
        guide = (TUTORIALS / f"{stem}.md").read_text(encoding="utf-8")
        self.assertEqual(guide.splitlines()[0], "# " + source["title"])
        glossary = read_json(ROOT / "docs/glossary.json")
        def rows(value):
            if isinstance(value, dict):
                yield value
                for child in value.values():
                    yield from rows(child)
            elif isinstance(value, list):
                for child in value:
                    yield from rows(child)
        command = next(row for row in rows(glossary)
                       if row.get("path") == "planning protein-expression-handoff")
        self.assertIn("[--provider PROVIDER]", command["usage"])
        self.assertIn("provider-neutral", command["summary"])
        for path in ("docs/cli.md", "docs/protocol.md"):
            self.assertIn(command["usage"], (ROOT / path).read_text(encoding="utf-8"))
        for text in (guide, guide.replace("\n", "\r\n")):
            self.assertIn("gentle.protein_expression_handoff.v2", text)
            self.assertIn("Legacy v1 JSON remains readable", text)
            self.assertIn("not_selected", text)
            self.assertIn("outsourcing_disallowed", text)
            self.assertIn("direct_submission_available = false", text)
            for field in ("neutral_service_draft", "provider_projection", "protein_seq_id",
                          "unresolved question IDs", "no vendor-side action",
                          "coding DNA as `unresolved`", "`not_performed`",
                          "`eligible = true`", "`quote_handoff_available = true`",
                          "`select_protein_expression_provider`"):
                self.assertIn(field, text)
            self.assertNotIn("prepare_geneart_quote_packet_after_review", text)
            blocks = re.findall(r"```bash\r?\n(.*?)```", text, re.S)
            self.assertTrue(blocks)
            result = subprocess.run(["bash", "-n"], input="\n".join(blocks).replace("\r\n", "\n"),
                                    text=True, capture_output=True, timeout=10)
            self.assertEqual(result.returncode, 0, result.stderr)
        # This is explicitly a learner-owned input, not a missing bundled file.
        self.assertIn("learner\ninput, not a bundled file", guide)
        for path in re.findall(r'@((?:docs|test_files|assets)/[^\s\'"`]+)', guide):
            self.assertTrue((ROOT / path).is_file(), path)

    def test_expression_builder_is_neutral_and_cli_no_write_cases_are_authored(self):
        builder = (ROOT / "src/engine_shell/protein_expression_service_draft.rs").read_text(encoding="utf-8")
        self.assertIn("ExternalServiceDeliveryRouteRequest", builder)
        self.assertIn("project_external_service_delivery_request", builder)
        self.assertNotIn("fs::", builder)
        self.assertNotIn("geneart_protein_expression_request.json", builder)
        self.assertNotIn('"unresolved_product"', builder)
        for field in ("cds_assessment", "coding_dna_status", "reverse_translation_status",
                      "service_intent", "nonclaims", "review_status"):
            self.assertIn(f'"{field}"', builder)
        shell = (ROOT / "src/engine_shell.rs").read_text(encoding="utf-8")
        self.assertIn("preflight.eligible && preflight.quote_handoff_available", shell)
        cli = (ROOT / "tests/protein_expression_handoff_cli.rs").read_text(encoding="utf-8")
        for case in ("neutral", "selected", "protein", "in_house", "excluded", "partial",
                     "missing_product", "missing_permission", "unsupported_capability",
                     "quote_unavailable", "blocked_preflight", "unknown_provider", "invalid_product"):
            self.assertIn(f'"{case}"', cli)
        self.assertRegex(cli, r"snapshot\(root\),\s*before")
        self.assertIn("quote_shell_arg(token)", cli)
        self.assertIn('preflight.args(["services", "project-preflight", &request])', cli)
        self.assertIn('"{case}: local preflight must not write"', cli)
        self.assertIn('"\\r\\n"', cli)
        self.assertIn("No network", (ROOT / "tests/protein_expression_handoff_README.md").read_text(encoding="utf-8"))

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

    def test_09_04_literal_walkthrough_and_catalog_update(self):
        """Opt-in exact prose examples with a prebuilt CLI; no build/network/model."""
        bash = shutil.which("bash")
        self.assertIsNotNone(bash, "Bash is required to replay this POSIX-shell chapter")
        env = dict(self.env, PATH=str(self.cli.parent) + os.pathsep + os.environ.get("PATH", ""),
                   TMPDIR=str(self.run_dir))

        def replay(path, prefix=""):
            text = path.read_text(encoding="utf-8")
            script = prefix + "\n".join(re.findall(r"```bash\n(.*?)```", text, re.S))
            result = subprocess.run([bash, "-eu"], input=script, cwd=ROOT, env=env,
                                    text=True, capture_output=True, timeout=240)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

        replay(PLANNING_GUIDE)
        run, = self.run_dir.glob("gentle-planning-09-04.*")
        self.state = run / "tutorial.gentle.json"
        self.env.update(GENTLE_SYSTEM_CONFIG_ROOT=str(run / "system"),
                        XDG_CONFIG_HOME=str(run / "user"), GENTLE_PROJECT_ROOT=str(run / "project"))
        state_before = read_json(self.state)
        self.assertEqual(set(state_before["sequences"]), {"synthetic_panel_vector", "synthetic_panel_source"})
        proposal = read_json(run / "vector-proposal.json")
        self.assertEqual(proposal["vector_validation"]["status"], "verified")
        self.assertEqual(proposal["vector_validation"]["helper_catalog_id"], "Synthetic panel vector")
        self.assertEqual(proposal["vector_validation"]["observed_length_bp"], 240)
        self.assertTrue(proposal["approval_required"])
        self.assertEqual({p["allele"] for p in proposal["products"]}, {"wild_type", "mutant"})
        self.assertFalse((run / "planned-products").exists())
        for name in ("promoter_reporter_panel_demo.synthetic_vector.svg", "promoter_reporter_panel_demo.source.svg"):
            self.assertEqual(ET.fromstring((run / name).read_bytes()).tag, "{http://www.w3.org/2000/svg}svg")
        sheet = (run / "cloning-decision.md").read_text(encoding="utf-8")
        self.assertIn("HUMAN REVIEW REQUIRED", sheet)
        self.assertIn(proposal["vector_sequence_sha256"], sheet)
        original_ids = {p["provider"] for p in read_json(run / "suppliers.json")["providers"]}

        # Restore the documented shell function, not a different test-side dispatcher.
        setup = re.search(r'^g\(\) \{.*$', PLANNING_GUIDE.read_text(encoding="utf-8"), re.M).group(0)
        env.update(RUN=str(run), STATE=str(self.state), GENTLE_CLI=str(self.cli),
                   GENTLE_PROJECT_ROOT=str(run / "project"),
                   GENTLE_SYSTEM_CONFIG_ROOT=str(run / "system"), XDG_CONFIG_HOME=str(run / "user"))
        replay(CATALOG_EXERCISE, setup + "\n")
        added = read_json(run / "provider-preview.json")
        self.assertTrue(added["applicable"])
        self.assertEqual(added["decisions"][0]["action"], "add")
        self.assertIn("stale_or_tampered_preview", (run / "stale-error.txt").read_text())
        conflict = read_json(run / "conflict-preview.json")
        self.assertFalse(conflict["applicable"])
        self.assertEqual(conflict["decisions"][0]["action"], "conflict")
        self.assertEqual(read_json(run / "withdraw-preview.json")["decisions"][0]["action"], "withdraw")
        self.assertEqual(read_json(run / "withdrawn-provider.json")["rows"], [])
        self.assertEqual({p["provider"] for p in read_json(run / "suppliers-after-update.json")["providers"]}, original_ids)
        receipt = read_json(run / "provider-receipt.json")
        history = Path(receipt["history_directory"])
        self.assertFalse((history / "previous.bin").exists())
        added_bytes = (history / "proposed.json").read_bytes()
        self.assertEqual("sha256:" + hashlib.sha256(added_bytes).hexdigest(), receipt["applied_target_sha256"])
        withdrawn = read_json(run / "withdraw-receipt.json")
        self.assertEqual((Path(withdrawn["history_directory"]) / "previous.bin").read_bytes(), added_bytes)
        target = Path(withdrawn["target_path"])
        self.assertEqual("sha256:" + hashlib.sha256(target.read_bytes()).hexdigest(), withdrawn["applied_target_sha256"])
        self.assertIn("synthetic_tutorial_supplier", read_json(target)["withdrawn_providers"])
        # Retain exact project values; only existing unordered restriction-cache rows are normalized.
        state_after = read_json(self.state)
        for state in (state_before, state_after):
            state["sequences"] = comparable_sequences(state["sequences"])
        self.assertEqual(state_after, state_before)

    def test_cloning_planning_profiles_and_protocol_supplier_boundary(self):
        estimates = []
        for name in ("cloning_in_stock", "cloning_waiting"):
            profile = ROOT / f"docs/examples/planning/{name}.json"
            self.invoke("planning", "profile", "set", f"@{profile}", "--scope", "project_override")
            report = self.invoke("routines", "explain", "restriction.digest_ligate_extract_sticky")
            estimate = report["planning"]["estimate"]
            estimates.append(estimate)
            self.assertEqual(estimate["schema"], "gentle.planning_estimate.v2")
            self.assertAlmostEqual(estimate["hands_on_hours"]["typical"], 0.8)
            self.assertAlmostEqual(estimate["unattended_hours"]["typical"], 2.6)
            self.assertAlmostEqual(estimate["estimated_cost"], 39.0)
            self.assertAlmostEqual(estimate["machine_run_cost"], 30.0)
            policy = estimate["explanation"]["time_policy"]
            self.assertEqual(policy["basis"], "legacy_family_scalar")
            self.assertEqual(policy["used_by"], ["composite_meta_score", "max_time_hours"])
            self.assertAlmostEqual(policy["hours"], 6.0 if name == "cloning_in_stock" else 90.0)
            overview = report["planning"]["protocol_overview"]
            self.assertEqual(overview["schema"], "gentle.planning_protocol_overview.v1")
            self.assertEqual(overview["status"], "available")
            self.assertEqual(len(overview["nodes"]), 7)
            self.assertTrue(overview["source_estimate_sha256"].startswith("sha256:"))
            svg = ET.fromstring(overview["svg"])
            self.assertEqual(svg.attrib["data-source-estimate-sha256"], overview["source_estimate_sha256"])
            self.assertEqual({node.attrib["data-step-id"] for node in svg.iter() if "data-step-id" in node.attrib},
                             {node["step"]["step_id"] for node in overview["nodes"]})
            self.assertEqual({(node.attrib["data-from"],node.attrib["data-to"]) for node in svg.iter() if "data-from" in node.attrib},
                             {(edge["from"],edge["to"]) for edge in overview["edges"]})
            self.assertAlmostEqual(estimate["estimated_cost"], sum(estimate[field] for field in (
                "routine_base_cost", "stock_drawn_material_cost",
                "missing_material_purchase_proxy_cost", "machine_run_cost")))
            if name == "cloning_in_stock":
                guarded = self.invoke("planning", "consult", "cloning", "--objective", json.dumps({
                    "schema": "gentle.planning_objective.v1", "max_time_hours": 5}), "--format", "json")
                restriction = next(row["estimate"] for row in guarded["strategy_candidates"]
                                   if row["family"] == "restriction")
                self.assertAlmostEqual(restriction["estimated_time_hours"], 3.2)
                self.assertFalse(restriction["passes_guardrails"])
                self.assertIn("max_time_hours_exceeded (6.0 > 5.0)", restriction["guardrail_failures"])
        self.assertAlmostEqual(estimates[0]["elapsed_hours"]["typical"], 3.2)
        self.assertAlmostEqual(estimates[1]["elapsed_hours"]["typical"], 87.2)
        self.assertAlmostEqual(estimates[0]["stock_drawn_material_cost"], 9.0)
        self.assertAlmostEqual(estimates[1]["stock_drawn_material_cost"], 2.0)
        self.assertAlmostEqual(estimates[1]["missing_material_purchase_proxy_cost"], 7.0)
        legacy_report = self.invoke("routines", "explain", "gibson.two_fragment_overlap_preview")
        legacy = legacy_report["planning"]["estimate"]
        self.assertEqual(legacy_report["planning"]["protocol_overview"]["status"], "unavailable")
        self.assertEqual(legacy_report["planning"]["protocol_overview"]["nodes"], [])
        self.assertEqual(legacy["step_model"], "absent")
        self.assertEqual(legacy["estimated_time_hours_basis"], "legacy_family_scalar")
        for field in ("hands_on_hours", "unattended_hours", "elapsed_hours"):
            self.assertIsNone(legacy[field])
        before = self.state.read_bytes()
        helper_id = "Reporter Promega Luciferase AY738222 (online)"
        consult = self.invoke("planning", "consult", "cloning", "--objective", json.dumps({
            "schema": "gentle.planning_objective.v1", "helper_profile_id": helper_id}), "--format", "json")
        self.assertEqual(consult["objective_summary"]["objective"]["helper_profile_id"], helper_id)
        self.assertTrue(any(row["helper_id"] == helper_id for row in consult["vector_candidates"]))
        catalog = self.invoke("services", "providers", "list")
        promega = next(row for row in catalog["providers"] if row["provider"] == "promega")
        self.assertEqual(promega["country_code"], "US")
        self.assertEqual(promega["protocols"][0]["format"], "text_pdf")
        german = [row for row in catalog["providers"] if row.get("country_code") == "DE" and row.get("protocols")]
        self.assertIn("iba_lifesciences", [row["provider"] for row in german])
        request = {"schema": "gentle.external_service_request.v1", "provider": "promega",
                   "service_kind": "vector_supply", "source_target": {"kind": "catalog_vector", "catalog_number": "E6651"}}
        preflight = self.invoke("services", "project-preflight", json.dumps(request))
        self.assertEqual(preflight["capability_status"], "catalog_only")
        self.assertFalse(preflight["eligible"])
        self.assertFalse(preflight["quote_handoff_available"])
        self.assertFalse(preflight["direct_submission_available"])
        self.assertIsNone(preflight["estimated_turnaround"])
        self.assertIsNone(preflight["estimated_cost_hint"])
        # Deliberately contradictory synthetic overlay: metadata must not
        # promote a catalog-only listing to an implemented adapter.
        override = copy.deepcopy(promega)
        override["default_delivery_hints"] = ["SYNTHETIC_UNVERIFIED_TURNAROUND"]
        capability = override["capabilities"][0]
        capability.update(quote_handoff_supported=True, direct_api_documented=True,
                          direct_api_implemented=True)
        overlay = self.run_dir / "project/.gentle/catalogs/external_service_providers.json"
        overlay.parent.mkdir(parents=True)
        overlay.write_text(json.dumps({"schema": "gentle.external_service_provider_config.v1",
                                       "providers": [override]}), encoding="utf-8")
        request["source_target"] = {"kind": "inline_dna", "sequence": "ATGGCTTAA"}
        preflight = self.invoke("services", "project-preflight", json.dumps(request))
        self.assertEqual(preflight["capability_status"], "catalog_only")
        self.assertFalse(preflight["eligible"])
        self.assertFalse(preflight["quote_handoff_available"])
        self.assertFalse(preflight["direct_submission_available"])
        self.assertIsNone(preflight["estimated_turnaround"])
        quote = self.invoke("services", "project-quote", json.dumps(request))
        self.assertEqual(quote["quote_status"], "blocked")
        self.assertEqual(self.state.read_bytes(), before)

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
