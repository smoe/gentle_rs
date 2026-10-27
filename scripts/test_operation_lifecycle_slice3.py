"""Lifecycle Slice 3/4 source and opt-in actual-process synthetic file tests.

All input is invented inline or tutorial 06.07's existing literal CDS. Temporary
projects/sidecars are recreated through the supplied CLI, never private inputs.
Set GENTLE_LIFECYCLE_BIN_DIR to same-source gentle_cli/gentle_mcp binaries.
This companion never builds, downloads, calls a model or opens a native GUI.
"""

import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import tomllib
import unittest

ROOT = Path(__file__).resolve().parents[1]
BIN_DIR = os.environ.get("GENTLE_LIFECYCLE_BIN_DIR")
REQUEST = ROOT / "docs/tutorial/inputs/synthetic_sequence_design.json"


def plan(path=None, loaded=False):
    request = json.loads(REQUEST.read_bytes())
    if loaded:
        request["target"] = {"kind": "loaded_sequence", "seq_id": "synthetic_source"}
    return {"PlanDnaSequenceDesign": {"request": request, "path": path}}


def snv():
    return {"EvaluateGenomeEditIntent": {"request": {
        "schema": "gentle.genome_edit_intent_request.v1",
        "target": {"kind": "inline_sequence", "sequence": "ACGTAC"},
        "reference": {"kind": "synthetic", "namespace": "synthetic:slice3", "reference_id": "six-base-window"},
        "anchor": {"start_1based": 1001, "end_1based": 1006, "strand": "plus"},
        "expected_sequence_sha256": "sha256:e76797dcfb9f3a4639fdf013e8a72d0bd943a105bfc397387335c6540c496559",
        "design_haplotype": "reference",
        "intent": {"start_0based": 1001, "end_0based_exclusive": 1002, "deleted": "C", "inserted": "T"},
    }}}


def workflow(ops):
    return {"run_id": "synthetic-lifecycle", "ops": ops}


def quote(text):
    # GENtle's documented shared shell-quote rule, not the platform shell.
    return '"' + text.replace("\\", "\\\\").replace('"', '\\"') + '"'


class LifecycleSlice3SourceTests(unittest.TestCase):
    def test_pure_explicit_selector_and_factories(self):
        policy = (ROOT / "src/engine/operation_lifecycle.rs").read_text()
        selector = policy.split("#[cfg(test)]", 1)[0]
        for required in ("LifecycleRequest", "Operation::PlanDnaSequenceDesign", "Operation::EvaluateGenomeEditIntent",
                         "!workflow.ops.is_empty()", "preserves_project_content"):
            self.assertIn(required, selector)
        for forbidden in ("std::fs", "serde_json", "resolve_operation_payload", "ShellCommand"):
            self.assertNotIn(forbidden, selector)
        engine = (ROOT / "src/engine.rs").read_text()
        factory = engine.split("pub fn from_state_for(", 1)[1].split("fn reseed_op_counter", 1)[0]
        self.assertIn("Self::from_state(state)", factory)
        self.assertIn("engine.reseed_op_counter_from_state()", factory)
        self.assertNotIn("engine.reconcile_", factory)

    def test_no_new_payload_resolution_boundary(self):
        for file in ("src/bin/gentle_cli.rs", "src/mcp_server.rs", "src/js_interface.rs", "src/lua_interface.rs"):
            source = (ROOT / file).read_text()
            self.assertIn("GentleEngine::from_state_for", source)
            self.assertNotIn("ResolvedLifecycle", source)
        shell = (ROOT / "src/engine_shell.rs").read_text()
        self.assertNotIn("LifecycleRequest", shell, "shared shell parsing/execution stays unchanged")

    def test_synthetic_payloads_survive_lf_crlf(self):
        for payload in (plan(), snv(), workflow([plan(), snv()])):
            for ending in ("\n", "\r\n"):
                encoded = json.dumps(payload, indent=2).replace("\n", ending).encode()
                self.assertEqual(json.loads(encoded), payload)


class LifecycleSlice4SourceTests(unittest.TestCase):
    def test_plan_only_staging_and_host_binding_are_wired(self):
        writer = (ROOT / "src/engine/io/sequence_design_preview.rs").read_text()
        production = writer.split("#[cfg(test)]", 1)[0]
        for required in ("fs::symlink_metadata", "fs::canonicalize", "directory_identity",
                         "candidate_store_sidecar_dir", "eq_ignore_ascii_case", "unsupported_name_comparison",
                         "NamedTempFile::new_in", "sync_all", "persist_noclobber", "serde_json::to_string_pretty"):
            self.assertIn(required, production)
        for forbidden in ("create_dir_all", "fs::write", "temp.persist(", "resolve_operation_payload"):
            self.assertNotIn(forbidden, production)
        engine = (ROOT / "src/engine.rs").read_text()
        self.assertIn("#[serde(skip)]\n    host_project_path: Option<PathBuf>", engine)
        self.assertIn("host_project_path: self.host_project_path.clone()", engine)
        self.assertIn("self.host_project_path = live.host_project_path.clone()", engine)
        handler = (ROOT / "src/engine/ops/operation_handlers.rs").read_text()
        plan_arm = handler.split("Operation::PlanDnaSequenceDesign { request, path } =>", 1)[1].split("Operation::ValidateRnaPrimerPanel", 1)[0]
        self.assertLess(plan_arm.index("sequence_design_preview_destination"), plan_arm.index("plan_dna_sequence_design"))
        self.assertIn("output.write(&report)", plan_arm)
        self.assertNotIn("write_pretty_json_file", plan_arm)

    def test_existing_locked_identity_crate_and_documented_limits(self):
        manifest = tomllib.loads((ROOT / "Cargo.toml").read_text())
        lock = tomllib.loads((ROOT / "Cargo.lock").read_text())
        self.assertEqual(manifest["dependencies"]["same-file"], "1.0.6")
        identity, = [p for p in lock["package"] if p["name"] == "same-file"]
        self.assertEqual(identity["version"], "1.0.6")
        root, = [p for p in lock["package"] if p["name"] == "GENtle"]
        self.assertIn("same-file", root["dependencies"])
        guide = (ROOT / "docs/tutorial/06-07_synthetic_sequence_design.md").read_text()
        for required in ("fresh filename", "existing directory", "sidecar", "Save As", "Unicode", "parent-directory"):
            self.assertIn(required, guide)


@unittest.skipUnless(BIN_DIR, "Supply same-source binaries through GENTLE_LIFECYCLE_BIN_DIR; no auto-build")
class LifecycleSlice3ProcessTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="gentle-lifecycle-")
        self.addCleanup(self.temp.cleanup)
        self.work = Path(self.temp.name)
        self.state = self.work / "synthetic's project.json"
        suffix = ".exe" if os.name == "nt" else ""
        self.cli = Path(BIN_DIR).resolve() / ("gentle_cli" + suffix)
        self.mcp = Path(BIN_DIR).resolve() / ("gentle_mcp" + suffix)
        self.assertTrue(self.cli.is_file(), self.cli)
        self.assertTrue(self.mcp.is_file(), self.mcp)
        self.env = dict(os.environ, GENTLE_CANDIDATE_STORE_STRICT_LOAD="0")

    def run_cli(self, *args, ok=True):
        run = subprocess.run([str(self.cli), "--state", str(self.state), *args], cwd=self.work,
                             env=self.env, capture_output=True, timeout=90)
        self.assertEqual(run.returncode == 0, ok, run.stderr.decode(errors="replace"))
        return json.loads(run.stdout) if ok else run

    def initialize(self, ending="\n", candidates=False):
        self.state.unlink(missing_ok=True)
        self.run_cli("op", json.dumps({"CreateSequenceFromText": {
            "sequence_text": "ATGTTTAAATAA", "output_id": "synthetic_source", "name": None, "circular": False}}))
        if candidates:
            self.run_cli("op", json.dumps({"GenerateCandidateSet": {
                "set_name": "windows", "seq_id": "synthetic_source", "length_bp": 4, "step_bp": 4,
                "feature_kinds": [], "feature_label_regex": None, "max_distance_bp": None,
                "feature_geometry_mode": None, "feature_boundary_mode": None, "feature_strand_relation": None, "limit": 32}}))
        value = json.loads(self.state.read_bytes())
        value["lineage"] = {"nodes": {}, "seq_to_node": {}, "edges": [], "macro_instances": []}
        value["container_state"] = {"containers": {}, "arrangements": {}}
        self.state.write_bytes((json.dumps(value, indent=1).replace("\n", ending) + ending + " ").encode())
        return value

    def mcp_call(self, tool, arguments, ok=True, default_state=False):
        state_args = {} if default_state else {"state_path": str(self.state)}
        body = json.dumps({"jsonrpc": "2.0", "id": 1, "method": "tools/call",
                           "params": {"name": tool, "arguments": {**state_args, **arguments}}}).encode()
        framed = f"Content-Length: {len(body)}\r\n\r\n".encode() + body
        run = subprocess.run([str(self.mcp)], input=framed, capture_output=True, cwd=self.work,
                             env=self.env, timeout=90)
        self.assertEqual(run.returncode, 0, run.stderr.decode(errors="replace"))
        header, raw = run.stdout.split(b"\r\n\r\n", 1)
        length = int(header.split(b":", 1)[1].strip())
        result = json.loads(raw[:length])["result"]
        self.assertEqual(not result.get("isError", False), ok, result)
        return json.loads(result["content"][0]["text"]) if ok else result

    def test_absent_empty_existing_projects_and_shell_spellings(self):
        for ending in ("\n", "\r\n"):
            for mode in ("missing", "empty", "existing"):
                with self.subTest(ending=repr(ending), mode=mode):
                    self.state.unlink(missing_ok=True)
                    if mode == "empty":
                        self.state.write_bytes(ending.encode())
                    elif mode == "existing":
                        self.initialize(ending)
                    before = self.state.read_bytes() if self.state.exists() else None
                    op_path = self.work / "op's payload.json"
                    wf_path = self.work / "workflow's payload.json"
                    op_path.write_bytes(json.dumps(plan(), indent=2).replace("\n", ending).encode())
                    wf_path.write_bytes(json.dumps(workflow([snv(), plan()]), indent=2).replace("\n", ending).encode())
                    request_path = self.work / "request.json"
                    request_path.write_text(json.dumps(plan()["PlanDnaSequenceDesign"]["request"]))
                    commands = [("op", "@" + str(op_path)), ("workflow", "@" + str(wf_path)),
                                ("shell", "op " + quote("@" + str(op_path))),
                                ("shell", "workflow " + quote("@" + str(wf_path))),
                                ("sequence-design", "plan", "@" + str(request_path))]
                    for command in commands:
                        self.run_cli(*command)
                        self.assertEqual(self.state.read_bytes() if self.state.exists() else None, before, command)

    def test_raw_mcp_state_summary_and_pre_io_confirmation(self):
        for tool, arguments in (("op", {"operation": plan()}), ("op", {"operation": snv()}),
                                ("workflow", {"workflow": workflow([plan(), snv()])})):
            self.mcp_call(tool, dict(arguments, confirm=True))
            self.assertFalse(self.state.exists())
        value = self.initialize("\r\n")
        before = self.state.read_bytes()
        for tool, arguments in (("op", {"operation": plan(loaded=True)}),
                                ("workflow", {"workflow": workflow([plan(loaded=True), snv()])})):
            out = self.mcp_call(tool, dict(arguments, confirm=True))
            self.assertEqual(out["state_summary"]["container_count"], 0)
            self.assertEqual(out["state_summary"]["sequence_count"], len(value["sequences"]))
            self.assertEqual(self.state.read_bytes(), before)
            denied = self.mcp_call(tool, dict(arguments, confirm=False), ok=False)
            self.assertIn("confirm=true", str(denied))
        # Malformed payload plus an unreadable project (directory) must still hit approval first.
        self.state = self.work
        for tool, key in (("op", "operation"), ("workflow", "workflow")):
            for confirm in ({}, {"confirm": False}):
                denied = self.mcp_call(tool, {key: {"malformed": True}, **confirm}, ok=False)
                self.assertIn("confirm=true", str(denied))
                self.assertNotIn("Could not load", str(denied))

    def test_explicit_report_failure_prefix_and_mutating_controls(self):
        self.initialize()
        before = self.state.read_bytes()
        report = self.work / "explicit-report.json"
        self.run_cli("op", json.dumps(plan(str(report), loaded=True)))
        self.assertTrue(report.is_file())
        self.assertEqual(self.state.read_bytes(), before)
        invalid = plan(loaded=True)
        invalid["PlanDnaSequenceDesign"]["request"]["target"]["seq_id"] = "missing"
        prefix = self.work / "workflow-prefix.json"
        self.run_cli("workflow", json.dumps(workflow([plan(str(prefix), loaded=True), invalid])), ok=False)
        self.assertTrue(prefix.is_file())
        self.assertEqual(self.state.read_bytes(), before)
        mcp_prefix = self.work / "mcp-workflow-prefix.json"
        self.mcp_call("workflow", {"workflow": workflow([plan(str(mcp_prefix), loaded=True), invalid]), "confirm": True}, ok=False)
        self.assertTrue(mcp_prefix.is_file())
        self.assertEqual(self.state.read_bytes(), before)
        for ops in ([], [snv(), {"SetParameter": {"name": "max_fragments_per_container", "value": 123}}]):
            self.initialize()
            before = self.state.read_bytes()
            self.run_cli("workflow", json.dumps(workflow(ops)))
            self.assertNotEqual(self.state.read_bytes(), before, "legacy aggregate must still save")

    def test_candidate_sidecar_bytes_and_permissive_diagnostics(self):
        value = self.initialize("\r\n", candidates=True)
        files = {p: p.read_bytes() for p in self.work.rglob("*") if p.is_file()}
        self.assertIn("index_path", value["metadata"]["candidate_sets"])
        self.run_cli("op", json.dumps(plan(loaded=True)))
        self.mcp_call("workflow", {"workflow": workflow([snv(), plan(loaded=True)]), "confirm": True})
        self.assertEqual({p: p.read_bytes() for p in files}, files)
        self.assertEqual({p for p in self.work.rglob("*") if p.is_file()}, set(files))
        index = self.work / value["metadata"]["candidate_sets"]["index_path"]
        index.unlink()
        before = self.state.read_bytes()
        self.run_cli("op", json.dumps(plan(loaded=True)))
        self.assertEqual(self.state.read_bytes(), before, "permissive loader warning must not be autosaved")
        self.env["GENTLE_CANDIDATE_STORE_STRICT_LOAD"] = "1"
        self.run_cli("op", json.dumps(plan(loaded=True)), ok=False)
        self.assertEqual(self.state.read_bytes(), before)

    def test_python_real_client_op_workflow_and_shell(self):
        sys.path.insert(0, str(ROOT / "integrations/python"))
        from gentle_py.client import GentleClient
        client = GentleClient(state_path=str(self.state), cli_cmd=[str(self.cli)], cwd=str(self.work))
        client.op(plan())
        self.assertFalse(self.state.exists())
        client.workflow(workflow([snv(), plan()]))
        self.assertFalse(self.state.exists())
        self.initialize("\r\n")
        before = self.state.read_bytes()
        client.op(plan(loaded=True))
        client.workflow(workflow([snv(), plan(loaded=True)]))
        client.shell("op " + quote(json.dumps(plan(loaded=True))), expect_json=True)
        self.assertEqual(self.state.read_bytes(), before)

    def test_slice4_cli_mcp_and_defaults_refuse_project_and_sidecar_destinations(self):
        operation = plan(str(self.state))
        request = self.work / "request.json"
        request.write_text(json.dumps(operation["PlanDnaSequenceDesign"]["request"]))
        for command in [("op", json.dumps(operation)), ("workflow", json.dumps(workflow([operation]))),
                        ("shell", "op " + quote(json.dumps(operation))),
                        ("sequence-design", "plan", "@" + str(request), "--path", str(self.state))]:
            result = self.run_cli(*command, ok=False)
            self.assertIn("sequence_design_report.protected_path", result.stderr.decode(errors="replace"))
            self.assertFalse(self.state.exists())
        default = self.work / ".gentle_state.json"
        result = subprocess.run([str(self.cli), "op", json.dumps(plan(str(default)))], cwd=self.work,
                                env=self.env, capture_output=True, timeout=90)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("sequence_design_report.protected_path", result.stderr.decode(errors="replace"))
        self.mcp_call("op", {"operation": plan(str(default)), "confirm": True}, ok=False, default_state=True)
        self.assertFalse(default.exists())
        for baseline in (b"\r\n ", None):
            if baseline is None:
                self.initialize("\r\n")
                baseline = self.state.read_bytes()
            else:
                self.state.write_bytes(baseline)
            self.run_cli("op", json.dumps(operation), ok=False)
            self.mcp_call("op", {"operation": operation, "confirm": True}, ok=False)
            self.assertEqual(self.state.read_bytes(), baseline)
        self.state.unlink()  # Disposable project fixture, not a retained report.
        sidecar = self.state.parent / (self.state.name + ".candidates")
        sidecar.mkdir()
        output = sidecar / "unused.json"
        self.run_cli("op", json.dumps(plan(str(output))), ok=False)
        self.mcp_call("op", {"operation": plan(str(output)), "confirm": True}, ok=False)
        self.assertFalse(output.exists())
        self.assertFalse(self.state.exists())

    def test_slice4_fresh_only_exact_bytes_and_duplicate_workflow_preserves_prefix(self):
        report = self.work / "fresh-report.json"
        result = self.run_cli("op", json.dumps(plan(str(report))))
        before = report.read_bytes()
        self.assertEqual(json.loads(before), result["dna_sequence_design"])
        refused = self.run_cli("op", json.dumps(plan(str(report))), ok=False)
        self.assertIn("destination_exists", refused.stderr.decode(errors="replace"))
        self.assertEqual(report.read_bytes(), before)
        missing = self.work / "absent-parent" / "report.json"
        self.run_cli("op", json.dumps(plan(str(missing))), ok=False)
        self.assertFalse(missing.parent.exists())
        prefix = self.work / "duplicate-prefix.json"
        self.mcp_call("workflow", {"workflow": workflow([plan(str(prefix)), plan(str(prefix))]), "confirm": True}, ok=False)
        self.assertEqual(json.loads(prefix.read_bytes()), json.loads(before))
        self.assertFalse(self.state.exists())

    def test_slice4_python_uses_host_binding_and_native_aliases_preserve_files(self):
        sys.path.insert(0, str(ROOT / "integrations/python"))
        from gentle_py.client import GentleClient, GentleCliError
        client = GentleClient(state_path=str(self.state), cli_cmd=[str(self.cli)], cwd=str(self.work))
        with self.assertRaises(GentleCliError):
            client.op(plan(str(self.state)))
        self.assertFalse(self.state.exists())
        sentinel = self.work / "sentinel.json"
        sentinel.write_bytes(b"unchanged\r\n")
        hard = self.work / "hardlink.json"
        os.link(sentinel, hard)
        client_path = self.work / "python-report.json"
        client.op(plan(str(client_path)))
        with self.assertRaises(GentleCliError):
            client.op(plan(str(hard)))
        self.assertEqual(sentinel.read_bytes(), b"unchanged\r\n")
        self.assertTrue(os.path.samefile(sentinel, hard))


if __name__ == "__main__":
    unittest.main()
