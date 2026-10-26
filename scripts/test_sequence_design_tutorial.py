"""Executable companion for the hand-crafted synthetic tutorial 06.07.

Recreation and provenance: docs/tutorial/inputs/README.md. The independent
four-word oracle is literal teaching data, not a second production solver.
Actual CLI tests require GENTLE_TUTORIAL_BIN_DIR; they authorize only this
synthetic fixture, never a user's DNA. No builds, network, model or GUI capture.
"""

import json
import copy
import hashlib
import os
import platform
from pathlib import Path
import re
import subprocess
import tempfile
import time
import unittest

from scripts import sequence_design_audit_oracle as oracle


ROOT = Path(__file__).resolve().parents[1]
GUIDE = ROOT / "docs/tutorial/06-07_synthetic_sequence_design.md"
SOURCE = ROOT / "docs/tutorial/sources/06-07_synthetic_sequence_design.json"
REQUEST = ROOT / "docs/tutorial/inputs/synthetic_sequence_design.json"
BIN_DIR = os.environ.get("GENTLE_TUTORIAL_BIN_DIR")
OUTPUT = "ATGTTTAAGTAA"


def fenced_request(text):
    return json.loads(re.search(r"```json\r?\n(.*?)\r?\n```", text, re.S).group(1))


class SequenceDesignTutorialSourceTests(unittest.TestCase):
    def test_request_aware_artifact_approval_is_wired_without_lifecycle_claims(self):
        policy = (ROOT / "src/engine/operation_lifecycle.rs").read_text(encoding="utf-8")
        selector = policy.split("#[cfg(test)]", 1)[0]
        for text in ("Operation::PlanDnaSequenceDesign", "path: Some(_)", "path: None",
                     "ConfirmArtifactWriteAttempt", "NoAdditionalConfirmation", "_ => None"):
            self.assertIn(text, selector)
        self.assertNotIn("serde_json", selector)
        self.assertNotIn("std::fs", selector)
        shell = (ROOT / "src/engine_shell.rs").read_text(encoding="utf-8")
        approval = shell.split("pub fn is_state_mutating(", 1)[1].split("pub fn ", 1)[0]
        self.assertIn("operation_approval_policy(&operation)", approval)
        self.assertIn("policy.requires_confirmation()", approval)
        descriptor = shell.split("fn sequence_design_capability_descriptor(", 1)[1].split("fn tss_view_svg_capability_descriptor(", 1)[0]
        self.assertIn('else {"external"}', descriptor)
        self.assertIn('"requires_confirmation":true', descriptor)
        self.assertIn("write attempt", descriptor)
        bridge = (ROOT / "src/agent_bridge.rs").read_text(encoding="utf-8")
        self.assertIn('requires execution="ask" for the exact report-write attempt', bridge)
        guide = GUIDE.read_text(encoding="utf-8")
        self.assertIn("approval for the report-write attempt", guide)

    def test_discovery_and_links(self):
        source = json.loads(SOURCE.read_text(encoding="utf-8"))
        self.assertEqual(source["schema"], "gentle.tutorial_source.v4")
        self.assertNotIn("generated_chapter", source)
        self.assertEqual((source["catalog"]["group"],
                          source["catalog"]["group_position"]), ("06", 7))
        catalog = json.loads((ROOT / "docs/tutorial/catalog.json").read_text(encoding="utf-8"))
        entry, = [row for row in catalog["entries"] if row["id"] == source["id"]]
        self.assertEqual(entry["title"], source["title"])
        self.assertEqual(entry["notes"], source["catalog"]["notes"])
        self.assertEqual(entry["network"], "offline")
        text = GUIDE.read_text(encoding="utf-8")
        self.assertEqual(text.splitlines()[0], f'# {source["title"]}')
        for target in re.findall(r"\]\(([^)]+)\)", text):
            if "://" not in target:
                self.assertTrue((GUIDE.parent / target.split("#")[0]).exists(), target)

    def test_documented_request_is_exact_in_lf_and_crlf(self):
        text = GUIDE.read_bytes().decode("utf-8")
        fixture = json.loads(REQUEST.read_bytes())
        for newline in ("\n", "\r\n"):
            with self.subTest(newline=repr(newline)):
                self.assertEqual(fenced_request(text.replace("\n", newline)), fixture)
                self.assertEqual(json.loads(REQUEST.read_bytes().replace(b"\n", newline.encode())),
                                 fixture)
        self.assertEqual(fixture["purpose"], "synthetic_coding_insert")
        self.assertEqual(fixture["target"]["sequence"], "ATGTTTAAATAA")
        self.assertEqual(fixture["cds"], {"start_0based": 0, "end_0based_exclusive": 12})
        self.assertEqual(fixture["protein_sequence"], "MFK")

    def test_literal_four_word_oracle_and_teaching_table(self):
        request = json.loads(REQUEST.read_bytes())
        self.assertEqual(request["avoid_motifs"], [{"pattern": "TTTAAA", "strand": "both"}])
        self.assertEqual("TTTAAA".translate(str.maketrans("ACGT", "TGCA"))[::-1], "TTTAAA")
        bounds = request["gc_content"]
        local = request["gc_window"]
        self.assertEqual((bounds["min_basis_points"], bounds["max_basis_points"]), (1666, 1667))
        self.assertEqual(local, dict(bounds, window_bp=6))
        minimum = (bounds["min_basis_points"] * 12 + 9999) // 10000
        maximum = bounds["max_basis_points"] * 12 // 10000
        self.assertEqual((minimum, maximum), (2, 2))
        self.assertEqual(((local["min_basis_points"] * 6 + 9999) // 10000,
                          local["max_basis_points"] * 6 // 10000), (1, 1))
        # All legal words in this tiny frozen-start/stop space, independent of GENtle.
        words = ["ATGTTTAAATAA", "ATGTTCAAATAA", OUTPUT, "ATGTTCAAGTAA"]
        residues = {"ATG": "M", "TTT": "F", "TTC": "F", "AAA": "K", "AAG": "K"}
        observations = []
        for word in words:
            self.assertEqual("".join(residues[word[pos:pos + 3]] for pos in (0, 3, 6)), "MFK")
            self.assertEqual((word[:3], word[-3:]), ("ATG", "TAA"))
            gc = sum(base in "GC" for base in word)
            windows = [sum(base in "GC" for base in word[start:start + 6])
                       for start in range(7)]
            row = (gc, "TTTAAA" not in word, gc == 2, all(count == 1 for count in windows))
            observations.append(row)
            cells = [str(gc), *("yes" if value else "no" for value in row[1:])]
            self.assertIn(f'| `{word}` | ' + " | ".join(cells) + " |",
                          GUIDE.read_text(encoding="utf-8"))
        self.assertEqual(observations, [(1, False, False, False), (2, True, True, False),
                                       (2, True, True, True), (3, True, False, False)])

    def test_agent_drafts_and_review_boundaries(self):
        source = json.loads(SOURCE.read_bytes())
        text = GUIDE.read_text(encoding="utf-8")
        self.assertEqual(len(source["agent_parity"]["cases"]), 2)
        for case in source["agent_parity"]["cases"]:
            self.assertEqual(case["execution"], "ask")
            self.assertIn(case["command"], text)
        preview, = [case for case in source["agent_parity"]["cases"]
                    if case["id"] == "preview_joint_constraints"]
        self.assertEqual(preview["command"], "sequence-design plan @request.json --path preview.json")
        self.assertFalse(preview["mutating"], "Report-write approval is not project mutation")
        for required in ("ACTUAL_REVIEWED_DIGEST", "session", "not a digest to submit",
                         "search_claims_verified=false", "unverified_portable_preview",
                         "Do not pipe a digest automatically", "no applicable",
                         "not recommended", "outside", "Glen"):
            self.assertIn(required, text)
        self.assertIn("synthetic_sequence_design.json",
                      (REQUEST.parent / "README.md").read_text(encoding="utf-8"))


@unittest.skipUnless(BIN_DIR or os.environ.get("GENTLE_SEQUENCE_DESIGN_EVIDENCE_DIR"),
                     "Set GENTLE_TUTORIAL_BIN_DIR for actual sequence-design CLI replay")
class SequenceDesignTutorialCliTests(unittest.TestCase):
    def setUp(self):
        retained = os.environ.get("GENTLE_SEQUENCE_DESIGN_EVIDENCE_DIR")
        self.candidate = os.environ.get("GENTLE_SEQUENCE_DESIGN_CANDIDATE")
        self.ledger = []
        self.binary_hash = None
        if retained:
            self.assertTrue(BIN_DIR, "Retained acceptance requires GENTLE_TUTORIAL_BIN_DIR; cannot skip")
            self.assertRegex(self.candidate or "", r"^[0-9a-f]{40}$",
                             "Retained acceptance requires an exact authorized candidate")
            revision = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT,
                                      capture_output=True, check=True, timeout=10)
            self.assertEqual(revision.stdout.decode().strip(), self.candidate,
                             "Acceptance checkout differs from the authorized candidate")
            status = subprocess.run(["git", "status", "--porcelain"], cwd=ROOT,
                                    capture_output=True, check=True, timeout=10)
            self.assertEqual(status.stdout, b"", "Acceptance requires a clean isolated checkout")
            self.cwd = Path(retained).resolve() / self._testMethodName
            self.cwd.mkdir(parents=True, exist_ok=False)
            self.addCleanup(self.retain_evidence)
        else:
            self.tmp = tempfile.TemporaryDirectory(prefix="gentle design's ")
            self.addCleanup(self.tmp.cleanup)
            self.cwd = Path(self.tmp.name)
        self.state = self.cwd / "tutorial.gentle.json"
        suffix = ".exe" if os.name == "nt" else ""
        self.cli = Path(BIN_DIR).resolve() / f"gentle_cli{suffix}"
        self.assertTrue(self.cli.is_file(), f"Missing explicitly configured CLI: {self.cli}")
        self.binary_hash = hashlib.sha256(self.cli.read_bytes()).hexdigest()
        if retained:
            version = subprocess.run([str(self.cli), "--version"], cwd=self.cwd,
                                     capture_output=True, timeout=30, env=self.isolated_env())
            (self.cwd / "binary-version.stdout").write_bytes(version.stdout)
            (self.cwd / "binary-version.stderr").write_bytes(version.stderr)
            self.assertEqual(version.returncode, 0)
            self.assertRegex(version.stdout.decode(), r"(?m)^Source revision [^\r\n]*\+git\." +
                             re.escape(self.candidate) + r"\r?$")

    def isolated_env(self):
        env = {key: value for key, value in os.environ.items()
               if not key.startswith(("GENTLE_", "XDG_"))}
        for key in ("HOME", "USERPROFILE", "XDG_CONFIG_HOME", "XDG_CACHE_HOME", "XDG_DATA_HOME"):
            folder = self.cwd / "isolated-profile" / key.lower()
            folder.mkdir(parents=True, exist_ok=True)
            env[key] = str(folder)
        return env

    def retain_evidence(self):
        # Results are observations, not an auditor's scientific or performance verdict.
        current_hash = (hashlib.sha256(self.cli.read_bytes()).hexdigest()
                        if getattr(self, "cli", None) and self.cli.is_file() else None)
        record = {"schema": "gentle.sequence_design_audit_observations.v1",
                  "candidate": self.candidate, "binary": str(getattr(self, "cli", "")),
                  "binary_sha256": self.binary_hash, "platform": platform.platform(),
                  "test": self._testMethodName, "commands": self.ledger,
                  "binary_sha256_at_end": current_hash,
                  "binary_unchanged": self.binary_hash is not None and current_hash == self.binary_hash,
                  "companion_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  "oracle_sha256": hashlib.sha256(Path(oracle.__file__).read_bytes()).hexdigest(),
                  "independent_acceptance": False,
                  "artifacts": {str(p.relative_to(self.cwd)): hashlib.sha256(p.read_bytes()).hexdigest()
                                for p in sorted(self.cwd.rglob("*"))
                                if p.is_file() and p != self.cwd / "observations.json"}}
        (self.cwd / "observations.json").write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
        if self.binary_hash is not None:
            self.assertEqual(current_hash, self.binary_hash, "Audit binary changed during replay")

    def invoke(self, *args, success=True):
        command = [str(self.cli), "--state", str(self.state), *args]
        index = len(self.ledger)
        started = time.monotonic()
        row = {"argv": command, "expected_success": success, "state": str(self.state)}
        self.ledger.append(row)
        try:
            result = subprocess.run(command, cwd=self.cwd, capture_output=True, timeout=90,
                                    env=self.isolated_env())
        except subprocess.TimeoutExpired as error:
            (self.cwd / f"command-{index:03}.stdout").write_bytes(error.stdout or b"")
            (self.cwd / f"command-{index:03}.stderr").write_bytes(error.stderr or b"")
            row.update(exit_code=None, timeout=True, wall_seconds=time.monotonic() - started)
            raise
        row.update(exit_code=result.returncode, timeout=False, wall_seconds=time.monotonic() - started)
        (self.cwd / f"command-{index:03}.stdout").write_bytes(result.stdout)
        (self.cwd / f"command-{index:03}.stderr").write_bytes(result.stderr)
        if not success:
            self.assertNotEqual(result.returncode, 0, result.stdout.decode(errors="replace"))
            return result
        self.assertEqual(result.returncode, 0, result.stderr.decode(errors="replace"))
        return json.loads(result.stdout)

    def test_exhaustive_acceptance_matrix_against_literal_oracle(self):
        for name, template in oracle.audit_cases():
            observations = []
            for strategy in ("full_enumeration", "conflict_directed"):
                with self.subTest(case=name, strategy=strategy):
                    item = copy.deepcopy(template)
                    item["search_strategy"] = strategy
                    self.state = self.cwd / f"{name}-{strategy}.json"
                    _, _, report = self.plan(item, name + " " + strategy)
                    self.assertFalse(self.state.exists(), "Read-only preview created project state")
                    oracle.assert_report(item, report)
                    # Check our adversarial encoder against authentic engine bytes first.
                    if report["approval_digest"] is not None:
                        self.assertEqual(oracle.approval_digest(report), report["approval_digest"])
                    observations.append(report["output_sequence"])
            if len(observations) == 2:
                self.assertEqual(observations[0], observations[1], name)

    def test_rehashed_output_and_gc_facts_do_not_bypass_apply(self):
        for strategy in ("full_enumeration", "conflict_directed"):
            item = dict(oracle.audit_cases())["joint"]
            item["search_strategy"] = strategy
            self.state = self.cwd / f"tamper-{strategy}.json"
            self.invoke("op", json.dumps({"CreateSequenceFromText": {
                "sequence_text": "ACGT", "output_id": "keep_me", "name": None, "circular": False}}))
            baseline = self.state.read_bytes()
            _, _, report = self.plan(item, "tamper " + strategy)
            self.assertEqual(oracle.approval_digest(report), report["approval_digest"])
            for kind in ("gc_output", "window_fact", "stripped_gc", "stripped_window", "stop"):
                with self.subTest(strategy=strategy, tamper=kind):
                    changed = copy.deepcopy(report)
                    if kind in ("gc_output", "stop"):
                        changed["output_sequence"] = ("ATGTTCAAGTAA" if kind == "gc_output"
                                                       else "ATGTTTAAGTAG")
                        changed["output_sha256"] = "sha256:" + hashlib.sha256(changed["output_sequence"].encode()).hexdigest()
                        changed["edits"] = [{"position_0based": i, "before": a, "after": b}
                                            for i, (a, b) in enumerate(zip(changed["source_sequence"],
                                                                          changed["output_sequence"])) if a != b]
                    elif kind == "window_fact":
                        changed["gc_window"]["output"][0]["gc_bases"] += 1
                    else:
                        del changed["request"]["gc_content" if kind == "stripped_gc" else "gc_window"]
                    changed["approval_digest"] = oracle.approval_digest(changed)
                    path = self.cwd / f"{strategy} {kind}.json"
                    path.write_text(json.dumps(changed), encoding="utf-8")
                    self.invoke("sequence-design", "apply", "@" + str(path), "--approve",
                                changed["approval_digest"], success=False)
                    self.assertEqual(self.state.read_bytes(), baseline)

            # Valid two-edit MEF DNA may be applied, but its fabricated search proof
            # must stay explicitly unverified. This reproduces Glen's needs_fix.
            item = oracle.request("ATGGAATTCTAA", "MEF")
            item["avoid_motifs"] = [{"pattern": "GAATTC", "strand": "both"}]
            item["search_strategy"] = strategy
            _, _, submitted = self.plan(item, "nonminimum " + strategy)
            self.assertEqual(oracle.approval_digest(submitted), submitted["approval_digest"])
            submitted["output_sequence"] = "ATGGAGTTTTAA"
            submitted["output_sha256"] = "sha256:" + hashlib.sha256(b"ATGGAGTTTTAA").hexdigest()
            submitted["edits"] = [{"position_0based": 5, "before": "A", "after": "G"},
                                  {"position_0based": 8, "before": "C", "after": "T"}]
            submitted["evaluated_candidates"] = 0
            submitted["optimization_complete"] = True
            submitted["minimum_edits_proven"] = True
            submitted["approval_digest"] = oracle.approval_digest(submitted)
            path = self.cwd / f"{strategy} forged-search.json"
            path.write_text(json.dumps(submitted), encoding="utf-8")
            applied = self.invoke("sequence-design", "apply", "@" + str(path), "--approve",
                                  submitted["approval_digest"])
            receipt = applied["result"]["dna_sequence_design_receipt"]
            self.assertTrue(receipt["output_constraints_verified"])
            self.assertFalse(receipt["search_claims_verified"])
            state = json.loads(self.state.read_bytes())
            record = state["metadata"]["dna_sequence_design:" + item["output_seq_id"]]
            self.assertEqual(record["submitted_proposal"]["verification"], "unverified_portable_preview")
            self.assertEqual(record["receipt"], receipt)
            containers = [c for c in state["container_state"]["containers"].values()
                          if item["output_seq_id"] in c["members"]]
            self.assertEqual(len(containers), 1)
            self.assertEqual(containers[0]["created_by_op"], applied["result"]["op_id"])
            self.assertEqual(containers[0]["name"], "Synthetic sequence design")
            nodes = [n for n in state["lineage"]["nodes"].values() if n["seq_id"] == item["output_seq_id"]]
            self.assertEqual(len(nodes), 1)
            self.assertEqual(nodes[0]["created_by_op"], applied["result"]["op_id"])

    def plan(self, request, name):
        path = self.cwd / f"{name} request.json"
        path.write_text(json.dumps(request), encoding="utf-8")
        report_path = self.cwd / f"{name} preview.json"
        result = self.invoke("sequence-design", "plan", f"@{path}", "--path", str(report_path))
        report = json.loads(report_path.read_bytes())
        self.assertEqual(result["result"]["dna_sequence_design"], report)
        return path, report_path, report

    def test_both_strategies_direct_and_shared_shell_preview_then_reviewed_apply(self):
        for strategy in ("conflict_directed", "full_enumeration"):
            with self.subTest(strategy=strategy):
                self.state = self.cwd / f"{strategy}.gentle.json"
                # Inline sentinel DNA is synthetic, unrelated to the redesign request.
                self.invoke("op", json.dumps({"CreateSequenceFromText": {
                    "sequence_text": "ACGT", "output_id": "keep_me", "name": None,
                    "circular": False}}))
                before = self.state.read_bytes()
                initial = json.loads(before)
                request = json.loads(REQUEST.read_bytes())
                request["search_strategy"] = strategy
                path, report_path, report = self.plan(request, strategy)
                self.assertEqual(self.state.read_bytes(), before, "Preview changed project bytes")
                self.assertEqual(report["status"], "feasible")
                self.assertEqual(report["output_sequence"], OUTPUT)
                self.assertEqual(report["edits"], [{"position_0based": 8, "before": "A", "after": "G"}])
                self.assertEqual(report["initial_matches"], [{"motif_index": 0,
                    "interval": {"start_0based": 3, "end_0based_exclusive": 9}, "strand": "both"}])
                self.assertEqual(report["gc_content"]["output"]["gc_bases"], 2)
                self.assertEqual([row["gc_bases"] for row in report["gc_window"]["input"]],
                                 [1, 1, 1, 0, 0, 0, 0])
                rows = report["gc_window"]["output"]
                self.assertEqual(len(rows), 7)
                for start, row in enumerate(rows):
                    self.assertEqual(row, {"interval": {"start_0based": start,
                        "end_0based_exclusive": start + 6}, "gc_bases": 1, "satisfies_bounds": True})
                self.assertTrue(report["optimization_complete"])
                self.assertTrue(report["minimum_edits_proven"])
                self.assertTrue(report["approval_digest"].startswith("sha256:"))
                shell_result = self.invoke("shell", "sequence-design plan " +
                                           json.dumps("@" + path.as_posix()))
                self.assertEqual(shell_result["result"]["dna_sequence_design"], report)
                self.assertEqual(self.state.read_bytes(), before)
                for args in ((), ("--approve", "sha256:not-reviewed")):
                    self.invoke("sequence-design", "apply", f"@{report_path}", *args, success=False)
                    self.assertEqual(self.state.read_bytes(), before)
                applied = self.invoke("sequence-design", "apply", f"@{report_path}",
                                      "--approve", report["approval_digest"])
                receipt = applied["result"]["dna_sequence_design_receipt"]
                self.assertTrue(receipt["output_constraints_verified"])
                self.assertFalse(receipt["search_claims_verified"])
                state = json.loads(self.state.read_bytes())
                self.assertEqual(state["sequences"]["keep_me"], initial["sequences"]["keep_me"])
                product = state["sequences"][request["output_seq_id"]]
                self.assertEqual(bytes(product["seq"]["seq"]).decode("ascii"), OUTPUT)
                record = state["metadata"]["dna_sequence_design:" + request["output_seq_id"]]
                self.assertEqual(record["submitted_proposal"], {
                    "verification": "unverified_portable_preview", "report": report})
                self.assertEqual(record["receipt"], receipt)
                saved = self.state.read_bytes()
                self.invoke("sequence-design", "apply", f"@{report_path}",
                            "--approve", report["approval_digest"], success=False)
                self.assertEqual(self.state.read_bytes(), saved, "Duplicate apply changed state")
                history = self.invoke("shell", "history status")
                self.assertEqual(history["undo_count"], 0)
                self.invoke("shell", "history undo", success=False)
                self.assertEqual(self.state.read_bytes(), saved,
                                 "A fresh CLI process must not invent undo history")

    def test_exhausted_search_has_no_approval_or_created_output(self):
        request = json.loads(REQUEST.read_bytes())
        request["max_evaluations"] = 1
        for strategy in ("conflict_directed", "full_enumeration"):
            with self.subTest(strategy=strategy):
                request["search_strategy"] = strategy
                _, path, report = self.plan(request, "exhausted " + strategy)
                self.assertEqual(report["status"], "search_exhausted")
                self.assertIsNone(report["output_sequence"])
                self.assertIsNone(report["approval_digest"])
                self.assertIsNone(report["gc_window"]["output"])
                self.invoke("sequence-design", "apply", f"@{path}",
                            "--approve", "sha256:not-reviewed", success=False)
                self.assertFalse(self.state.exists())


if __name__ == "__main__":
    unittest.main()
