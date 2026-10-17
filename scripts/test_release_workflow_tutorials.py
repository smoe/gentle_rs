"""Discovery contracts and opt-in runtime replay for tutorials 03.10 and 09.03.

Default: Python/source checks only. Set GENTLE_TUTORIAL_BIN_DIR to an already
built candidate directory to run the real CLI. Never builds Rust, invokes a
model/provider, runs OpenSCAD or certifies physical/scientific acceptance.
Runtime fixtures are the committed synthetic starters and reviewed requirements
documented in docs/examples/plans/README.md; outputs live only in temporary
directories. A requested but missing binary is an error, not a skip.
"""

import copy
import json
import os
from pathlib import Path
import re
import shlex
import subprocess
import sys
import tempfile
import unittest
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[1]
BIN_DIR = os.environ.get("GENTLE_TUTORIAL_BIN_DIR")
STEMS = ("03-10_gene_set_storage_projection", "09-03_protein_expression_handoff_cli")


def read_json(path):
    return json.loads(path.read_bytes())


def tutorial_cli_commands(text):
    """Read CLI argv from fenced examples without executing Cargo or a shell."""
    commands = []
    for block in re.findall(r"```bash\r?\n(.*?)```", text, re.S):
        for line in block.replace("\r\n", "\n").replace("\\\n", "").splitlines():
            args = shlex.split(line)
            if args[:2] == ["cargo", "run"] and "--" in args:
                separator = args.index("--")
                if args[separator - 2:separator] == ["--bin", "gentle_cli"]:
                    commands.append(args[separator + 1:])
    return commands


def assert_manual_registration(test, root):
    """Check the two discovery projections, not the full Rust generator."""
    catalog = read_json(root / "docs/tutorial/catalog.json")["entries"]
    review = read_json(root / "docs/tutorial/review_manifest.json")["entries"]
    manifest = read_json(root / "docs/tutorial/manifest.json")["chapters"]
    groups = {g["code"]: g for g in
              read_json(root / "docs/tutorial/sources/catalog_meta.json")["groups"]}
    index = (root / "docs/tutorial/README.md").read_text(encoding="utf-8")
    for stem in STEMS:
        with test.subTest(tutorial=stem):
            source = read_json(root / f"docs/tutorial/sources/{stem}.json")
            identifier = source["id"]
            rows = [row for row in catalog if row["id"] == identifier]
            test.assertEqual(len(rows), 1)
            row = rows[0]
            expected = dict(source["catalog"])
            expected.pop("order")
            group = groups[expected["group"]]
            expected.update(id=identifier, title=source["title"],
                            group_label=group["label"], group_order=group["order"],
                            decimal_id=stem[:5].replace("-", "."),
                            review_status="unreviewed", review_issue_template="Tutorial confusion",
                            review_issue_template_path=".github/ISSUE_TEMPLATE/tutorial-confusion.md")
            test.assertEqual(row, expected)
            test.assertEqual(row["status"], "manual/hybrid")
            test.assertEqual(row["source"], "hand_written_markdown")
            test.assertNotIn("generated_chapter", source)
            test.assertFalse(any(chapter["id"] == identifier for chapter in manifest))
            reviews = [item for item in review if item["tutorial_id"] == identifier]
            test.assertEqual(len(reviews), 1)
            test.assertEqual(reviews[0]["tutorial_kind"], "guided_walkthrough")
            test.assertEqual(reviews[0]["tutorial_status"], "active")
            test.assertIsNone(reviews[0]["human_reviewed_at"])
            guide = (root / row["path"]).read_text(encoding="utf-8")
            test.assertEqual(guide.splitlines()[0], "# " + source["title"])
            test.assertIn(f"(./{stem}.md)", index)


class ReleaseTutorialSourceTests(unittest.TestCase):
    def test_manual_discovery_matches_sources_without_claiming_generated_replay(self):
        assert_manual_registration(self, ROOT)

    def test_storage_commands_are_literal_and_newline_independent(self):
        guide = (ROOT / f"docs/tutorial/{STEMS[0]}.md").read_text(encoding="utf-8")
        commands = tutorial_cli_commands(guide)
        self.assertEqual(commands, tutorial_cli_commands(guide.replace("\n", "\r\n")))
        self.assertEqual(len(commands), 6)
        self.assertTrue(all(args[:2] == ["--state", "storage-review.gentle.json"]
                            for args in commands))
        self.assertEqual(commands[0][2:], ["workflow", "@docs/examples/workflows/gene_set_storage_starter.json"])
        request = json.loads(commands[1][-1])["CreateGeneSetArrangement"]
        self.assertIs(request["confirm_physical_placement"], True)
        self.assertEqual([args[2:4] for args in commands[2:]],
                         [["racks", verb] for verb in ("show", "labels-svg", "isometric-svg", "openscad")])
        self.assertIn("`collection_operation` in the direct CLI output", guide)

    def test_ci_replays_with_the_built_cli_on_all_three_platforms(self):
        workflow = (ROOT / ".github/workflows/ci.yml").read_text(encoding="utf-8")
        for job in ("macos", "linux", "windows"):
            # Read the job block up to the next two-space job heading.
            steps = re.split(r"\n  [a-z][a-z0-9_-]*:\n",
                             workflow.split(f"\n  {job}:\n", 1)[1], maxsplit=1)[0]
            self.assertIn("GENTLE_TUTORIAL_BIN_DIR: ${{ github.workspace }}/target/debug", steps)
            self.assertIn("scripts.test_release_workflow_tutorials.ReleaseTutorialBinaryTests -v", steps)
            self.assertLess(steps.index("Tutorial core CLI runtime smoke"),
                            steps.index("Expression and storage tutorial CLI replay"))

    def test_missing_requested_cli_fails_instead_of_skipping(self):
        with tempfile.TemporaryDirectory() as folder:
            env = dict(os.environ, GENTLE_TUTORIAL_BIN_DIR=str(Path(folder) / "absent"),
                       PYTHONDONTWRITEBYTECODE="1")
            result = subprocess.run(
                [sys.executable, "-m", "unittest",
                 "scripts.test_release_workflow_tutorials.ReleaseTutorialBinaryTests", "-v"],
                cwd=ROOT, env=env, text=True, capture_output=True, timeout=30)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("FAILED (failures=2)", result.stderr)
            self.assertNotIn("skipped", result.stderr)


@unittest.skipUnless(BIN_DIR, "Set GENTLE_TUTORIAL_BIN_DIR for real CLI replay; no build performed")
class ReleaseTutorialBinaryTests(unittest.TestCase):
    def setUp(self):
        suffix = ".exe" if os.name == "nt" else ""
        self.cli = Path(BIN_DIR).expanduser().resolve() / f"gentle_cli{suffix}"
        self.assertTrue(self.cli.is_file() and os.access(self.cli, os.X_OK), str(self.cli))
        temp = tempfile.TemporaryDirectory(prefix="gentle release tutorial's ")
        self.addCleanup(temp.cleanup)
        self.run_dir = Path(temp.name).resolve()
        self.state = self.run_dir / "review project.gentle.json"
        self.env = {key: value for key, value in os.environ.items()
                    if not key.startswith("GENTLE_") and key != "RUST_MIN_STACK"}
        for name, keys in (
            ("home", ("HOME", "USERPROFILE")), ("config", ("XDG_CONFIG_HOME",)),
            ("cache", ("XDG_CACHE_HOME",)), ("data", ("XDG_DATA_HOME",)),
            ("runtime", ("XDG_RUNTIME_DIR",)), ("tmp", ("TMPDIR", "TMP", "TEMP")),
            ("system", ("GENTLE_SYSTEM_CONFIG_ROOT",)), ("project", ("GENTLE_PROJECT_ROOT",)),
        ):
            directory = self.run_dir / name
            directory.mkdir()
            self.env.update({key: str(directory) for key in keys})
        self.env["GENTLE_ASSET_ROOT"] = str(ROOT)

    def invoke(self, *args, succeeds=True):
        args = ["@" + str(ROOT / arg[1:]) if arg.startswith("@docs/") else arg for arg in args]
        result = subprocess.run([str(self.cli), "--state", str(self.state), "--full-report", *args],
                                cwd=self.run_dir, env=self.env, text=True,
                                capture_output=True, timeout=90)
        if succeeds:
            self.assertEqual(result.returncode, 0, result.stderr + result.stdout)
            return json.loads(result.stdout)
        self.assertNotEqual(result.returncode, 0, result.stdout)
        return result

    def snapshot(self):
        return {str(path.relative_to(self.run_dir)): path.read_bytes() if path.is_file() else None
                for path in self.run_dir.rglob("*")}

    def test_expression_starter_neutral_withheld_and_selected_are_read_only(self):
        command = ["planning", "protein-expression-handoff", "--format", "json"]
        before = self.snapshot()
        missing = self.invoke(*command)
        self.assertEqual(missing["status"], "needs_product_definition")
        self.assertIsNone(missing["neutral_service_draft"])
        self.assertEqual(self.snapshot(), before)
        self.invoke("workflow", "@docs/examples/workflows/protein_expression_review_starter.json")
        before = self.snapshot()
        for fixture, provider, expected in (
            ("review_only", None, "not_selected"),
            ("in_house", None, "not_selected"),
            ("in_house", "geneart", "withheld"),
            ("review_only", "geneart", "projected"),
        ):
            with self.subTest(fixture=fixture, provider=provider):
                args = command + ["--seq-id", "expression_review_product", "--requirements",
                                  f"@docs/examples/plans/protein_expression_requirements_{fixture}.json"]
                if provider:
                    args += ["--provider", provider]
                report = self.invoke(*args)
                self.assertEqual(report["schema"], "gentle.protein_expression_handoff.v2")
                self.assertEqual(report["provider_binding"]["status"], expected)
                draft = report["neutral_service_draft"]
                self.assertEqual(draft["schema"], "gentle.external_service_delivery_route_request.v1")
                self.assertEqual(draft["source_target"]["seq_id"], "expression_review_product")
                if expected == "projected":
                    projection = report["provider_projection"]
                    self.assertEqual(projection["schema"], "gentle.external_service_request.v1")
                    self.assertEqual(projection["provider"], "geneart")
                    for field in ("source_target", "optimization_target", "vector_spec", "delivery_options"):
                        self.assertEqual(projection[field], draft[field])
                else:
                    self.assertIsNone(report["provider_projection"])
                    self.assertEqual(report["service_handoff_candidates"], [])
                    self.assertFalse(any(action["shell_line"].startswith("services project-")
                                         for action in report["suggested_next_actions"]))
                self.assertEqual(self.snapshot(), before, "Review must not write state or artifacts")

    def test_storage_starter_refusals_provenance_reopen_and_repeatable_exports(self):
        guide = (ROOT / f"docs/tutorial/{STEMS[0]}.md").read_text(encoding="utf-8")
        commands = [args[2:] for args in tutorial_cli_commands(guide)]
        results = self.invoke(*commands[0])
        self.assertEqual(results[-1]["op_id"], "op-3", "Tutorial's fresh-state identity drifted")
        starter = read_json(self.state)
        before = self.state.read_bytes()
        operation = json.loads(commands[1][1])
        for case in ("false", "omitted", "missing", "duplicate"):
            with self.subTest(refusal=case):
                invalid = copy.deepcopy(operation)
                request = invalid["CreateGeneSetArrangement"]
                if case == "false":
                    request["confirm_physical_placement"] = False
                elif case == "omitted":
                    del request["confirm_physical_placement"]
                elif case == "missing":
                    request["member_bindings"].pop()
                else:
                    request["member_bindings"][1] = copy.deepcopy(request["member_bindings"][0])
                self.invoke("op", json.dumps(invalid), succeeds=False)
                self.assertEqual(self.state.read_bytes(), before)
        result = self.invoke(*commands[1])
        report = result["collection_operation"]
        self.assertEqual(report["lifting_mode"], "arrange")
        self.assertEqual([row["member"]["stable_member_id"] for row in report["per_member_status"]],
                         ["symbol:genea", "symbol:geneb"])
        saved = read_json(self.state)
        record = saved["metadata"]["gene_set_arrangement:storage_review"]
        self.assertEqual(record["schema"], "gentle.gene_set_storage_projection.v1")
        self.assertEqual(record["report"], report)
        self.assertRegex(record["source_resolution_sha256"], r"^sha256:[0-9a-f]{64}$")
        self.assertEqual(saved["metadata"]["gene_set_artifacts"], starter["metadata"]["gene_set_artifacts"])
        self.assertEqual(saved["container_state"]["containers"], starter["container_state"]["containers"])
        self.assertEqual([row["source_seq_id"] for row in record["source_bindings"]],
                         ["storage_source_a", "storage_source_b"])
        for row in record["source_bindings"]:
            self.assertRegex(row["sequence_snapshot_sha256"], r"^sha256:[0-9a-f]{64}$")
        self.assertEqual(record["rack_id"], "rack-1")
        rack = self.invoke(*commands[2])  # A separate process reloads the saved state.
        self.assertEqual([(row["coordinate"], row["occupant"]["seq_id"]) for row in rack["placements"]],
                         [("A1", "storage_source_a"), ("A2", "storage_source_b")])
        before = self.state.read_bytes()
        self.invoke(*commands[1], succeeds=False)  # A second placement cannot reuse the ID.
        self.assertEqual(self.state.read_bytes(), before)
        for args in commands[3:]:
            with self.subTest(export=args[1]):
                path = self.run_dir / args[3]
                self.invoke(*args)
                first = path.read_bytes()
                self.assertTrue(first)
                if path.suffix == ".svg":
                    self.assertEqual(ET.fromstring(first).tag, "{http://www.w3.org/2000/svg}svg")
                else:
                    self.assertIn(b"difference()", first)
                    self.assertIn(b"storage_pcr_tube_rack", first)
                self.invoke(*args)
                self.assertEqual(path.read_bytes(), first)
        self.assertEqual(read_json(self.state)["metadata"]["gene_set_arrangement:storage_review"], record)


if __name__ == "__main__":
    unittest.main()
