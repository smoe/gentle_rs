"""Offline source/tutorial contracts, not Rust execution or physical acceptance.

All inputs are literal synthetic data documented in docs/examples/plans/README.md.
No new byte-exact generated-file check, network, model or OpenSCAD invocation.
"""

import json
from pathlib import Path
import re
import shlex
import subprocess
import unittest

ROOT = Path(__file__).resolve().parents[1]


class StorageProjectionTutorialTests(unittest.TestCase):
    def test_literal_starter_and_explicit_request_agree(self):
        workflow = json.loads((ROOT / "docs/examples/workflows/gene_set_storage_starter.json").read_text())
        ops = workflow["workflow"]["ops"]
        self.assertEqual([next(iter(op)) for op in ops], ["CreateSequenceFromText", "CreateSequenceFromText", "ResolveGeneSet"])
        members = ops[2]["ResolveGeneSet"]["source"]["members"]
        self.assertEqual(members, ["GENEA", "GENEB"])
        text = (ROOT / "docs/tutorial/03-10_gene_set_storage_projection.md").read_text()
        for prose in (text, text.replace("\n", "\r\n")):
            blocks = re.findall(r"```bash\r?\n(.*?)```", prose, re.S)
            result = subprocess.run(["bash", "-n"], input="\n".join(blocks).replace("\r\n", "\n"), capture_output=True, text=True, timeout=10)
            self.assertEqual(result.returncode, 0, result.stderr)
            command = next(block for block in blocks if "CreateGeneSetArrangement" in block)
            request = json.loads(shlex.split(command.replace("\\\n", ""))[-1])["CreateGeneSetArrangement"]
            self.assertIs(request["confirm_physical_placement"], True)
            self.assertEqual(request["resolution_id"], "op-3")
            self.assertEqual([binding["stable_member_id"] for binding in request["member_bindings"]], ["symbol:" + member.lower() for member in members])
            self.assertEqual([binding["source_container_id"] for binding in request["member_bindings"]], ["container-1", "container-2"])
        for term in ("manual/hybrid", "not a physical inventory", "not run locally", "confirm=true", "Glen/CI", "1..24"):
            self.assertIn(term, text)
        for link in re.findall(r"\]\(([^)#]+)", text):
            self.assertTrue((ROOT / "docs/tutorial" / link).resolve().is_file(), link)

    def test_operation_registry_and_engine_share_the_arrange_contract(self):
        policy = json.loads((ROOT / "docs/collection_lift_policies.json").read_text())
        row = next(item for item in policy["policies"] if item["name"] == "CreateGeneSetArrangement")
        subject = next(item for item in row["subjects"] if item["subject_kind"] == "gene_set_resolution")
        self.assertEqual(subject["context_requirement"], "context_agnostic")
        self.assertEqual(subject["support"]["mode"], "arrange")
        source = (ROOT / "src/engine/analysis/gene_set_arrangement.rs").read_text()
        self.assertIn("self.add_serial_arrangement(", source)
        self.assertNotIn("self.apply(", source)
        self.assertIn("confirm_physical_placement", source)
        self.assertIn("collection_sequence_snapshot_sha256", source)
        self.assertIn("source_resolution_sha256", source)
        self.assertIn("gene_set_arrangement_generic_shell_and_engine_agree", source)
        self.assertIn('"CreateGeneSetArrangement"', (ROOT / "crates/gentle-protocol/src/lib.rs").read_text())


if __name__ == "__main__":
    unittest.main()
