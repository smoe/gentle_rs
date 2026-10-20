"""Offline metadata regressions; all added evidence is synthetic temporary data.

The literal descriptor digests pin the pre-extension ac91375b contracts, not
scientific results. Evidence documents below are hand-crafted in disposable
directories and are never imported into the public evidence archive.
"""

import copy
import json
from pathlib import Path
import shutil
import tempfile
import unittest
from unittest.mock import patch

from scripts import generate_capability_maturity as generator

contract = generator.contract
ROOT = generator.ROOT
BASELINES = {
    "gentle-cloning": "sha256:f1fabf5bfba0cb5c05662c4382af5c2fa66aa238aac6401517db10990f4e89ae",
    "gentle-pcr-primer-design": "sha256:2939a343b2093db723d0744f8d0a2c336094edd59182e44bdc91781fdada32e8",
    "gentle-tss-collection": "sha256:ee83a1f6b878028cac8826252214097c3db7acacf4a7c60ba5bc47cb99c78813",
}


def source_paths():
    """Exact public generator inputs, including referenced workflow fixtures."""
    catalog = json.loads((ROOT / generator.OUTPUT_PATH).read_bytes())
    return {item["path"] for item in catalog["source_files"]} | {
        f"{generator.SKILL_DIR}/{skill}/SKILL.md" for skill in generator.SKILLS
    }


def populate(root):
    for relative in source_paths() | set(generator.generate()):
        path = root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(ROOT / relative, path)


class CapabilityMaturityTests(unittest.TestCase):
    def fixture(self):
        directory = tempfile.TemporaryDirectory(prefix="gentle-maturity-")
        self.addCleanup(directory.cleanup)
        root = Path(directory.name)
        populate(root)
        return root

    def test_generated_catalogs_are_current_and_unassessed(self):
        with patch("subprocess.run", side_effect=AssertionError("must not execute tools")):
            outputs = generator.generate()
        for path, expected in outputs.items():
            self.assertEqual((ROOT / path).read_bytes(), expected, path)
            catalog = contract.parse_json(expected)
            contract.validate_catalog(catalog)
            self.assertEqual(catalog["availability"], "not_checked")
            self.assertEqual(catalog["runtime_binding"], "not_checked")
            self.assertTrue(all(row["maturity"] == "unassessed" and row["evidence"] == []
                                for row in catalog["capabilities"]))
        self.assertEqual(len(contract.parse_json(outputs[generator.OUTPUT_PATH])["capabilities"]), 97)

    def test_discovery_is_short_and_capability_instructions_are_retained(self):
        for skill in generator.SKILLS:
            path = ROOT / generator.SKILL_DIR / skill
            text = (path / "SKILL.md").read_text(encoding="utf-8")
            name, description = generator.skill_discovery(text)
            self.assertEqual(name, skill)
            self.assertLessEqual(len(description), 1024)
            with self.assertRaisesRegex(ValueError, "1024"):
                generator.skill_discovery(text.replace("description: >-\n", "description: >-\n  " + "x" * 1025 + "\n", 1))
        scope = (ROOT / generator.SKILL_DIR / "gentle-cloning/references/capability_scope.md").read_text()
        self.assertIn("not clinical advice", scope)
        self.assertIn("No portal scraping", scope)

    def test_old_readers_keep_every_existing_route_and_execution_contract(self):
        for skill, expected in BASELINES.items():
            descriptor = contract.read_json(ROOT / generator.SKILL_DIR / skill / "INTENTS.json")[1]
            for route in descriptor["routes"]:
                self.assertEqual(route.pop("capability_id"), f"{skill}/{route['intent_id']}")
                steps = route.pop("capability_steps")
                self.assertEqual([step["step_index"] for step in steps], list(range(len(route["plan"]))))
            self.assertEqual(contract.digest(contract.canonical(descriptor)), expected)

    def test_short_discovery_does_not_replace_inner_agent_instructions(self):
        """Source-level boundary check, not a live-model or Rust acceptance test."""
        bridge = (ROOT / "src/agent_bridge.rs").read_text(encoding="utf-8")
        prompt = bridge.split('const AGENT_BRIDGE_SYSTEM_PROMPT: &str = r#"', 1)[1].split('"#;', 1)[0]
        for reference in (
            "gentle.agent_response.v1", "docs/glossary.json", "docs/cli.md",
            "docs/protocol.md", "docs/ai_prompt_contract.md", "docs/ai_cloning_primer.md",
            "docs/ai_task_playbooks.md", "docs/examples/ai_cloning_examples.md",
        ):
            self.assertIn(reference, prompt)
        context = bridge.split("pub fn build_agent_introspection_context_for_request(", 1)[1]
        context = context.split("/// One manifest-backed local reference", 1)[0]
        self.assertIn("graph: &ProjectFactGraph", context)
        self.assertIn("AGENT_INTROSPECTION_FACT_LIMIT", context)
        self.assertIn("omitted_fact_count", context)
        shell = (ROOT / "src/engine_shell.rs").read_text(encoding="utf-8")
        request = shell.split("fn execute_agents_ask_command(", 1)[1].split("let mut env_overrides", 1)[0]
        self.assertIn("build_agent_introspection_context_for_request(", request)
        self.assertIn("&engine.project_fact_graph()", request)
        for section in (prompt, context, request):
            for external_metadata in ("SKILL.md", "catalog_entry.json", "capability_maturity"):
                self.assertNotIn(external_metadata, section)

    def test_skill_info_and_runtime_projection_do_not_probe_or_change_capabilities_mode(self):
        wrapper = generator.wrapper
        script = ROOT / generator.SKILL_DIR / "gentle-cloning/gentle_cloning.py"
        with patch("subprocess.run", side_effect=AssertionError("no runtime probe")):
            info = wrapper._skill_info_payload(script)
            intents = wrapper._intents_runtime_payload(script)
        self.assertEqual(info["ui_intent_support"]["catalog_request_mode"], "capabilities")
        self.assertEqual(info["capability_maturity"]["status"], "available")
        self.assertEqual(info["capability_maturity"]["evidence_count"], 0)
        self.assertEqual(info["capability_maturity"]["runtime_binding"], "not_checked")
        self.assertEqual(wrapper._build_cli_args(wrapper.Request(mode="capabilities"), script), ["capabilities"])
        self.assertEqual(len(intents["routes"]), 59)
        self.assertTrue(all(row["capability_id"] and row["capability_steps"] for row in intents["routes"]))

    def test_missing_or_bad_sidecar_does_not_disable_discovery(self):
        root = self.fixture() / generator.SKILL_DIR / "gentle-cloning"
        sidecar = root / contract.FILENAME
        sidecar.unlink()
        self.assertEqual(contract.skill_info(root)["status"], "unavailable")
        sidecar.write_bytes(b'{"schema":"unknown"}\n')
        self.assertEqual(contract.skill_info(root)["status"], "invalid")
        info = generator.wrapper._skill_info_payload(root / "gentle_cloning.py")
        self.assertEqual(info["name"], "gentle-cloning")
        self.assertEqual(info["capability_maturity"]["status"], "invalid")

    def test_sidecar_rejects_stale_descriptor_catalog_adapter_or_request_bytes(self):
        root = self.fixture() / generator.SKILL_DIR / "gentle-cloning"
        for relative in ("INTENTS.json", "catalog_entry.json", "gentle_cloning.py",
                         "examples/request_skill_info.json"):
            path = root / relative
            before = path.read_bytes()
            try:
                path.write_bytes(before + b"\n")
                self.assertEqual(contract.skill_info(root)["status"], "invalid", relative)
            finally:
                path.write_bytes(before)
        self.assertEqual(contract.skill_info(root)["status"], "available")

    def test_generator_rejects_changed_step_targets_and_duplicate_ids(self):
        root = self.fixture()
        path = root / generator.SKILL_DIR / "gentle-cloning/INTENTS.json"
        original = path.read_bytes()
        descriptor = contract.parse_json(original)
        descriptor["routes"][0]["capability_steps"][0]["registry_refs"][0]["name"] = "primers design"
        path.write_bytes(generator.render(descriptor))
        with self.assertRaisesRegex(ValueError, "stale operation/command binding"):
            generator.generate(root)
        descriptor = contract.parse_json(original)
        descriptor["routes"].append(copy.deepcopy(descriptor["routes"][0]))
        path.write_bytes(generator.render(descriptor))
        with self.assertRaisesRegex(ValueError, "sorted and unique"):
            generator.generate(root)

    def test_workflow_operation_changes_cannot_keep_stale_refs(self):
        root = self.fixture()
        path = root / "docs/examples/workflows/simple_pcr_primer_design_offline.json"
        workflow = contract.read_json(path)[1]
        workflow["workflow"]["ops"][0] = {"ReverseComplement": {"input": "synthetic"}}
        path.write_bytes(generator.render(workflow))
        with self.assertRaisesRegex(ValueError, "stale operation/command binding"):
            generator.generate(root)

    def test_evidence_is_explicit_hashed_and_never_a_discovery_pass(self):
        root = self.fixture()
        catalog = contract.parse_json(generator.generate(root)[generator.OUTPUT_PATH])
        row = next(r for r in catalog["capabilities"] if r["intent_id"] == "primer_preflight")
        artifact = root / "synthetic-audit.json"
        artifact.write_bytes(b'{"synthetic":"not an actual audit"}\n')
        record = {
            "evidence_id": "synthetic-only", "capability_id": row["capability_id"],
            "definition_sha256": row["definition_sha256"],
            "scope": {"source_sha": "a" * 40, "tool_version": "synthetic", "binary_sha256": "sha256:" + "b" * 64,
                      "platform": "linux", "architecture": "x86_64", "build_profile": "synthetic",
                      "features": [], "input_class": "synthetic", "resources": {}},
            "kind": "execution", "verdict": "pass", "maturity": "verified",
            "artifact": {"path": "synthetic-audit.json", "sha256": contract.digest(artifact.read_bytes())},
            "assertions": ["synthetic metadata parser"], "limitations": ["not an audit or biological verdict"],
        }
        archive_path = root / generator.EVIDENCE_PATH
        archive = contract.read_json(archive_path)[1]
        archive["records"] = [record]
        archive_path.write_bytes(generator.render(archive))
        result = contract.parse_json(generator.generate(root)[generator.OUTPUT_PATH])
        updated = next(r for r in result["capabilities"] if r["capability_id"] == row["capability_id"])
        self.assertEqual(updated["maturity"], "unassessed")
        self.assertEqual(updated["evidence"], [record])
        artifact.write_bytes(b'{"synthetic":"tampered"}\n')
        with self.assertRaisesRegex(ValueError, "artifact changed"):
            generator.generate(root)
        for key, invalid in (("platform", "*"), ("platform", "unknown"), ("tool_version", "not_checked"),
                             ("source_sha", "main"), ("binary_sha256", "unknown")):
            altered = copy.deepcopy(record)
            altered["scope"][key] = invalid
            with self.assertRaises(ValueError):
                contract.validate_evidence(altered)

    def test_unknown_permission_fields_duplicate_keys_and_empty_results_are_rejected(self):
        catalog = contract.parse_json(generator.generate()[generator.OUTPUT_PATH])
        for key, value in (("auto_execute", True), ("verified", True)):
            with self.assertRaises(ValueError):
                contract.validate_catalog(dict(catalog, **{key: value}))
        for key in ("capabilities", "source_files"):
            with self.assertRaises(ValueError):
                contract.validate_catalog(dict(catalog, **{key: []}))
        for raw in (b'{"schema":1,"schema":2}', b'{"value":NaN}'):
            with self.assertRaises(ValueError):
                contract.parse_json(raw)

    def test_cross_platform_path_traversal_is_not_normalized_away(self):
        for path in ("../private", "/tmp/private", "C:/private", "C:\\private", "a/../b", "a//b"):
            self.assertFalse(contract.safe_relative_path(path), path)


if __name__ == "__main__":
    unittest.main()
