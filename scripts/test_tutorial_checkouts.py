"""Synthetic temporary Git fixtures for the cross-host tutorial checkout gate.

Fixtures are created below from literal bytes, committed only inside temporary
repositories, and consumed by the checkout/replay tests. The byte-preservation
regressions also copy the generated capability matrix and existing
provenance-documented TSS and probe fixtures.
No network, build, real checkout changes, or platform-specific tools are used.
"""

import hashlib
import json
import copy
import re
import os
from pathlib import Path
import subprocess
import tempfile
import tomllib
import unittest
from unittest.mock import patch
import xml.etree.ElementTree as ET
import zlib

from scripts import check_tutorial_checkouts as checker
from scripts import check_documentation_status as doc_status
from scripts import check_manual_migration as manual_migration
from scripts import manual_docs
from scripts import generate_capability_maturity as maturity_generator
from scripts.test_capability_maturity import populate as populate_maturity_fixture


def prepare_unprotected_checkout(source, destination, attributes):
    """Reproduce pre-policy CRLF drift after a caller removes a scoped rule.

    The repository-wide LF fallback would otherwise mask the missing scoped
    byte contract. This affects negative fixtures only, never the real policy.
    """
    attributes = b"\n".join(
        line for line in attributes.split(b"\n")
        if line.split() != [b"*", b"text=auto", b"eol=lf"]
    )
    return checker.prepare_checkout(source, destination, checker.MODES[1], attributes)


class DocumentationStatusTests(unittest.TestCase):
    """Copy public source inputs into synthetic disposable documentation repos."""

    def fixture(self):
        tmp = tempfile.TemporaryDirectory(prefix="gentle-doc-status-")
        self.addCleanup(tmp.cleanup)
        root = Path(tmp.name) / "source"
        paths = {
            "Cargo.toml", "README.md", "build.rs", "src/app.rs",
            "crates/gentle-shell/src/ui_intent.rs",
            "docs/reviewer_preview.md", "docs/primer_design_specialist_plan.md",
            "docs/cloning_pattern_catalog.md", "docs/cloning_routine_catalog_plan.md",
            "assets/cloning_routines.json", ".gitattributes",
        }
        for name in ["load_branch_reverse_complement_pgex_fasta",
                     "digest_ligation_extract_region_minimal"]:
            path = f"docs/examples/workflows/{name}.json"
            paths.add(path)
            paths.update(json.loads((checker.ROOT / path).read_text())["required_files"])
        paths.update(str(path.relative_to(checker.ROOT)) for path in
                     (checker.ROOT / "assets/cloning_patterns_catalog").rglob("*.json"))
        for relative in paths:
            destination = root / relative
            destination.parent.mkdir(parents=True, exist_ok=True)
            destination.write_bytes((checker.ROOT / relative).read_bytes())
        return root

    def test_current_documentation_status_and_links(self):
        self.assertEqual(doc_status.check(), [])

    def test_gui_manual_has_a_menu_first_entry_and_task_index(self):
        text = (checker.ROOT / "docs/gui.md").read_text(encoding="utf-8")
        intro = text.split("## Detailed Reference", 1)[0]
        for required in (
            "## Start Here", "## Find Your Task", "no shell commands are required",
            "File -> Open Tutorial Project...", "File -> Open Sequence...",
            "File -> Open Project...", "File -> Save Project...", "Cmd/Ctrl+K",
            "Windows", "Operation History", "session's undo stack",
            "| I Want To... | Start With |", "read-only inspector",
            "not biological validation",
        ):
            self.assertIn(required, intro)
        self.assertNotIn("cargo run", intro)
        headings = [re.sub(r"[^\w -]", "", line.lstrip("#").strip().lower()).replace(" ", "-")
                    for line in text.splitlines() if line.startswith("#")]
        for target in re.findall(r"\]\(([^)]+)\)", intro):
            path, _, fragment = target.partition("#")
            with self.subTest(link=target):
                if path:
                    self.assertTrue((checker.ROOT / "docs" / path).is_file(), target)
                elif fragment:
                    self.assertIn(fragment, headings)

    def test_gui_manual_starts_with_orientation_not_a_feature_appendix(self):
        source = (checker.ROOT / "docs/gui.md").read_bytes().decode("utf-8")
        for newline in ("\n", "\r\n"):
            text = source.replace("\r\n", "\n").replace("\n", newline)
            with self.subTest(newline=newline):
                preamble = text.split("## Start Here", 1)[0]
                blocks = re.split(r"\r?\n\s*\r?\n", preamble.strip())
                self.assertEqual(len(blocks), 3, "Only title, orientation and navigation belong before Start Here")
                self.assertEqual(blocks[0], "# GENtle GUI Manual")
                self.assertTrue(blocks[1].startswith("Open a sequence, explore its annotations"))
                self.assertTrue(blocks[2].startswith("**Quick navigation:**"))
                self.assertLess(len(preamble), 1000)
                self.assertGreater(text.index("## RNA Read Panel"), text.index("## Detailed Reference"))
                catalog = manual_docs.read_catalog()
                book = next(book for book in catalog["books"] if book["id"] == "gui")
                layout = "\n".join((checker.ROOT / topic["path"]).read_text(encoding="utf-8") for topic in book["topics"])
                for heading in ("Gene Set Inspector And Collection Operations", "Splicing Expert Workspace",
                                "RNA-Read Mapping Workspace", "Gibson Workspace", "Lineage Node Actions"):
                    self.assertIn(f"### {heading}", layout)

    def test_manual_migration_preserves_occurrences_and_indentation(self):
        # Synthetic text created here; no source-history access or Rust required.
        original = "# Manual\n\nSame paragraph.\n\nSame paragraph.\n\n```python\n    work()\n```\n"
        moved = ["# Index\r\n# Manual\r\nSame paragraph.\r\n", "Same paragraph.\n```python\n    work()\n```\n"]
        self.assertEqual(manual_migration.compare(original, moved)["missing"], [])
        duplicate_lost = manual_migration.compare(original, [moved[1]])["missing"]
        self.assertIn({"text": "Same paragraph.", "occurrences": 1}, duplicate_lost)
        changed_code = manual_migration.compare(original, [part.replace("    work()", "work()") for part in moved])["missing"]
        self.assertEqual(changed_code, [{"text": "    work()", "occurrences": 1}])

    def test_manual_migration_compares_relocated_links_without_weakening_text(self):
        original = "Read [primers](cli.md#primers), [here](#details).\n\n```md\n[x](literal.md)\n```\n"
        moved = original.replace("(cli.md#primers)", "(../cli.md#primers)").replace("(#details)", "(../gui.md#details)")
        before = manual_migration.canonical_links(original, "docs/gui.md")
        after = manual_migration.canonical_links(moved, "docs/gui/primers.md")
        self.assertEqual(before, after)
        self.assertIn("[x](literal.md)", after)
        wrong_target = manual_migration.canonical_links(moved.replace("../cli.md#primers", "../cli.md#rna"), "docs/gui/primers.md")
        self.assertTrue(manual_migration.compare(before, [wrong_target])["missing"])
        self.assertTrue(manual_migration.compare(before, [after.replace("Read", "Run")])["missing"])

    def test_manual_catalog_has_one_source_per_topic_and_current_projections(self):
        for relative, expected in manual_docs.generated().items():
            self.assertEqual((checker.ROOT / relative).read_bytes(), expected, relative)
        self.assertEqual(manual_docs.link_errors(), [])

    def test_manual_counterparts_are_reciprocal_and_missing_bindings_fail_closed(self):
        catalog = manual_docs.read_catalog()
        topics = {book["id"]: {t["id"]: t for t in book["topics"]} for book in catalog["books"]}
        self.assertEqual(set(topics), {"gui", "cli"})
        self.assertGreater(len(catalog["counterparts"]), 40)
        for gui, cli in catalog["counterparts"]:
            for book, topic, other in (("gui", gui, topics["cli"][cli]), ("cli", cli, topics["gui"][gui])):
                text = (checker.ROOT / topics[book][topic]["path"]).read_text()
                self.assertIn(f"../{'cli' if book == 'gui' else 'gui'}/{Path(other['path']).name}", text)
        for defect in ("unknown", "duplicate", "missing"):
            broken = copy.deepcopy(catalog)
            if defect == "unknown":
                broken["counterparts"][0][0] = "not-a-topic"
            elif defect == "duplicate":
                broken["counterparts"].append(broken["counterparts"][0])
            else:
                broken["counterparts"] = []
            with self.subTest(defect=defect), self.assertRaises(ValueError):
                manual_docs.validate(checker.ROOT, broken)

    def test_manual_chapters_have_balanced_fences_and_no_orphaned_nested_lists(self):
        for book in manual_docs.books(manual_docs.read_catalog()):
            for relative in manual_docs.paths(book):
                text = (checker.ROOT / relative).read_text()
                fence = None
                for line in text.splitlines():
                    mark = re.match(r"^ {0,3}(`{3,}|~{3,})", line)
                    if mark:
                        token = mark[1]
                        if fence is None:
                            fence = token
                        elif token[0] == fence[0] and len(token) >= len(fence):
                            fence = None
                self.assertIsNone(fence, relative)
                self.assertTrue(list(manual_docs.heading_rows(text)), relative)
                # Added CLI sections must give preserved nested command bullets
                # parents, not turn their four-space indentation into code.
                if book["id"] == "cli":
                    self.assertNotRegex(text, r"(?m)^## [^\n]+\n\n    - ", relative)

    def test_every_window_help_topic_has_a_catalogued_destination(self):
        catalog = manual_docs.read_catalog()
        gui = next(book for book in manual_docs.books(catalog) if book["id"] == "gui")
        topics = {topic["id"]: topic["path"] for topic in gui["topics"]}
        bindings = (checker.ROOT / "src/app/help_topics.rs").read_text(encoding="utf-8")
        destinations = re.findall(r'\w+ => \("([^"]+)", "([^"]*)"\)', bindings)
        self.assertGreaterEqual(len(destinations), 30)
        for topic, fragment in destinations:
            with self.subTest(topic=topic, fragment=fragment):
                self.assertIn(topic, topics)
                text = (checker.ROOT / topics[topic]).read_text(encoding="utf-8")
                if fragment:
                    self.assertIn(fragment, {row[3] for row in manual_docs.heading_rows(text)})

    def test_cli_manual_has_a_discovery_first_entry_and_task_index(self):
        text = (checker.ROOT / "docs/cli.md").read_text(encoding="utf-8")
        intro = text.split("## Detailed Reference", 1)[0]
        for newline in ("\n", "\r\n"):
            preamble = text.replace("\n", newline).split("## Start Here", 1)[0]
            blocks = re.split(r"\r?\n\s*\r?\n", preamble.strip())
            self.assertEqual(len(blocks), 3, "Only title, orientation and navigation belong before Start Here")
            self.assertEqual(blocks[0], "# GENtle CLI Manual")
            self.assertTrue(blocks[1].startswith("Inspect a sequence project"))
            self.assertTrue(blocks[2].startswith("**Quick navigation:**"))
            self.assertLess(len(preamble), 1000)
        self.assertIn("## Start Here", intro)
        self.assertIn("## Work With A Saved Project", intro)
        commands = re.search(r"```sh\n(.*?)\n```", intro, re.S).group(1)
        self.assertEqual(commands.splitlines(), [
            "gentle_cli --version", "gentle_cli help", "gentle_cli help primers design",
        ])
        self.assertIn('gentle_cli --state "my project.gentle.json" state-summary', intro)
        self.assertIn(".gentle_state.json", intro)
        self.assertIn("unsaved state", intro)
        index = text.split("### Commands\n", 1)[1].split("```bash", 1)[0]
        self.assertIn("| I Want To... | Start With |", index)
        self.assertIn("Reference examples, not a script", index)
        self.assertIn("disposable project", index)
        headings = [re.sub(r"[^\w -]", "", line.lstrip("#").strip().lower()).replace(" ", "-")
                    for line in text.splitlines() if line.startswith("#")]
        for target in re.findall(r"\]\(([^)]+)\)", intro + index):
            path, _, fragment = target.partition("#")
            with self.subTest(link=target):
                if path:
                    self.assertTrue((checker.ROOT / "docs" / path).is_file(), target)
                elif fragment:
                    self.assertIn(fragment, headings)

    def test_screenshot_docs_distinguish_gui_capture_from_disabled_cli_export(self):
        cli = (checker.ROOT / "docs/cli/screenshot-capture.md").read_text(encoding="utf-8")
        gui = (checker.ROOT / "docs/gui/screenshot-capture-policy.md").read_text(encoding="utf-8")
        for newline in ("\n", "\r\n"):
            with self.subTest(newline=newline):
                text = cli.replace("\n", newline)
                section = text.split("**Screenshot Capture: CLI Versus GUI**", 1)[1]
                section = section.split("Cryptic-splicing inspection:", 1)[0]
                for required in (
                    "CLI/shared-shell export remains disabled",
                    "screenshot-window is disabled by security policy",
                    "`--allow-screenshots` is rejected",
                    "GUI screenshot help is available", "Linux", "macOS",
                    "egui/eframe", "`Ask Agent`", "`Allow one screenshot`",
                    "macOS 14", "ScreenCaptureKit", "Screen Recording",
                    "`screenshot-capture`", "gui.md#screenshot-capture-policy",
                ):
                    self.assertIn(required, section)
                self.assertNotIn("Current status: disabled by security policy", section)
        self.assertIn("### Screenshot Capture Policy", gui)
        self.assertIn("Normal viewport capture does not require that permission", gui)
        self.assertNotIn("Current policy is manual documentation updates", gui)

    def test_preview_and_catalog_survive_real_lf_crlf_checkouts(self):
        root = self.fixture()
        checker.git(root, "init", "--quiet")
        checker.git(root, "config", "core.autocrlf", "false")
        checker.git(root, "config", "user.name", "Synthetic documentation test")
        checker.git(root, "config", "user.email", "fixture@example.invalid")
        checker.git(root, "add", "--all")
        checker.git(root, "-c", "commit.gpgsign=false", "commit", "--quiet", "-m", "public doc inputs")
        for mode in checker.MODES:
            with self.subTest(mode=mode[0]):
                clone = root.parent / mode[0]
                checker.prepare_checkout(root, clone, mode)
                payload = (clone / "docs/reviewer_preview.md").read_bytes()
                self.assertNotIn(b"\r", payload, "Repository policy keeps documentation LF")
                self.assertEqual(doc_status.readme_version_errors(clone), [])
                self.assertEqual(doc_status.preview_errors(clone), [])
                self.assertEqual(doc_status.catalog_errors(clone), [])

    def test_readme_release_guard_rejects_renamed_stale_and_missing_version(self):
        root = self.fixture()
        path = root / "README.md"
        original = path.read_text(encoding="utf-8")
        version = tomllib.loads((root / "Cargo.toml").read_text())["workspace"]["package"]["version"]
        marker = f"Current package version: `{version}`."
        self.assertIn(marker, original)
        for newline in ("\n", "\r\n"):
            for replacement in (f"Development version: `{version}`.",
                                "Current package version: `obsolete`.", ""):
                with self.subTest(newline=newline, replacement=replacement):
                    path.write_bytes(original.replace(marker, replacement).replace("\n", newline).encode("utf-8"))
                    self.assertEqual(doc_status.readme_version_errors(root),
                                     [f"README project status must name Cargo package version {version}"])

    def test_stale_version_and_unsafe_shared_state_are_refused(self):
        root = self.fixture()
        path = root / "docs/reviewer_preview.md"
        original = path.read_text()
        version = tomllib.loads((root / "Cargo.toml").read_text())["workspace"]["package"]["version"]
        for changed, expected in [
            (original.replace(f"GENtle `{version}`", "GENtle `obsolete`"), "workspace version"),
            (original.replace("mktemp -d", "fixed_path") + '\nrm -f "$STATE"\n', "destructive cleanup"),
            (original.replace("$WORK/digest.json", "$WORK/branch.json"), "separate fresh state"),
        ]:
            with self.subTest(expected=expected):
                path.write_text(changed)
                self.assertTrue(any(expected in error for error in doc_status.preview_errors(root)))

    def test_catalog_change_requires_documentation_reconciliation(self):
        root = self.fixture()
        path = root / "docs/cloning_pattern_catalog.md"
        count = len(list((root / "assets/cloning_patterns_catalog/gateway").rglob("*.json")))
        path.write_text(path.read_text().replace(f"| gateway | {count} |", f"| gateway | {count + 1} |"))
        self.assertTrue(any("template-family counts" in error for error in doc_status.catalog_errors(root)))

    def test_primer_plan_cannot_reintroduce_missing_gui_claim(self):
        root = self.fixture()
        path = root / "docs/primer_design_specialist_plan.md"
        path.write_text(path.read_text() + "\nThere is no dedicated, report-centric window.\n")
        self.assertTrue(any("existing PCR Designer" in error for error in doc_status.preview_errors(root)))

    def test_retired_link_is_refused_without_erasing_historical_prose(self):
        with tempfile.TemporaryDirectory(prefix="gentle-retired-doc-") as tmp:
            root = Path(tmp)
            plan = root / "docs/rna_read_origin_sparse_index_plan.md"
            plan.parent.mkdir()
            warning = "ROI capture layer is planned but not active in phase-1 scoring"
            original = (warning + "\nHistorical filename: tss_view_export_followup_codex_prompt.md\n"
                        "[archived request](https://github.com/smoe/gentle_rs/blob/"
                        + "a" * 40 + "/docs/tss_view_export_followup_codex_prompt.md)\n")
            plan.write_text(original)
            (plan.parent / "internal_12_candidate_preparation.md").write_text(
                "Historical pointer retained by newer main, not a retired active request.\n")
            source = root / "src/engine/analysis/rna_reads.rs"
            source.parent.mkdir(parents=True)
            source.write_text(warning)
            with patch.object(doc_status, "STATUS_DOCS", ()):
                self.assertEqual(doc_status.lifecycle_errors(root), [])
                for target in ("tss_view_export_followup_codex_prompt.md",
                               "tss%5Fview_export_followup_codex_prompt.md"):
                    plan.write_text(original + f"\n[old request]({target})\n")
                    self.assertTrue(any("still links to retired" in error
                                        for error in doc_status.lifecycle_errors(root)))
                plan.write_text(original)
                (plan.parent / "tss_view_export_followup_codex_prompt.md").write_text(
                    "Restored retired request.\n")
                self.assertTrue(any("completed request restored" in error
                                    for error in doc_status.lifecycle_errors(root)))

    def test_roi_capture_status_cannot_silently_become_implemented(self):
        with tempfile.TemporaryDirectory(prefix="gentle-roi-status-") as tmp:
            root = Path(tmp)
            plan = root / "docs/rna_read_origin_sparse_index_plan.md"
            plan.parent.mkdir()
            plan.write_text("ROI capture is implemented.\n")
            source = root / "src/engine/analysis/rna_reads.rs"
            source.parent.mkdir(parents=True)
            source.write_text("ROI capture layer is planned but not active in phase-1 scoring")
            with patch.object(doc_status, "STATUS_DOCS", ()):
                self.assertTrue(any("RNA capture" in error for error in doc_status.lifecycle_errors(root)))


class TutorialCheckoutTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix="gentle-checkout-test-")
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name) / "source"
        self.root.mkdir()
        checker.git(self.root, "init", "--quiet")
        checker.git(self.root, "config", "core.autocrlf", "false")
        checker.git(self.root, "config", "user.name", "Synthetic test")
        checker.git(self.root, "config", "user.email", "fixture@example.invalid")
        self.payload = b'{\n  "schema": "synthetic.fixture.v1"\n}\n'
        (self.root / "evidence.json").write_bytes(self.payload)
        (self.root / "unbound.txt").write_bytes(b"first\nsecond\n")
        (self.root / "binary.bin").write_bytes(b"\x00\xff\n\r\n")
        (self.root / ".gitattributes").write_bytes(b"evidence.json text eol=lf\n")
        checker.git(self.root, "add", "--all")
        checker.git(self.root, "-c", "commit.gpgsign=false", "commit", "--quiet", "-m", "fixture")

    def test_both_modes_preserve_bound_bytes_binary_and_history(self):
        revision = checker.git(self.root, "rev-parse", "HEAD").decode().strip()
        (self.root / "untracked-private.txt").write_text("must not be copied")
        before = checker.git(self.root, "status", "--porcelain")
        for mode in checker.MODES:
            with self.subTest(mode=mode[0]):
                target = Path(self.tmp.name) / mode[0]
                self.assertEqual(checker.prepare_checkout(self.root, target, mode), revision)
                self.assertEqual((target / "evidence.json").read_bytes(), self.payload)
                self.assertEqual((target / "binary.bin").read_bytes(), b"\x00\xff\n\r\n")
                separator = b"\r\n" if mode[0] == "crlf" else b"\n"
                self.assertEqual((target / "unbound.txt").read_bytes(),
                                 separator.join((b"first", b"second", b"")))
                self.assertFalse((target / "untracked-private.txt").exists())
                self.assertEqual(checker.git(target, "log", "-1", "--format=%H", "--", "evidence.json")
                                 .decode().strip(), revision)
        self.assertEqual(checker.git(self.root, "status", "--porcelain"), before)

    def test_manual_catalog_and_pages_survive_real_lf_crlf_checkouts(self):
        # Copy only the maintained manual inputs/projections. No image decoding,
        # GENtle binaries, Rust tools or personal runtime documentation are used.
        catalog = manual_docs.read_catalog()
        relative_paths = {manual_docs.CATALOG, manual_docs.EMBEDDED, ".gitattributes"}
        for book in manual_docs.books(catalog):
            relative_paths.update(manual_docs.paths(book))
        for relative in relative_paths:
            target = self.root / relative
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes((checker.ROOT / relative).read_bytes())
        checker.git(self.root, "add", "--all")
        checker.git(self.root, "-c", "commit.gpgsign=false", "commit", "--quiet", "-m", "manual navigation inputs")
        for mode in checker.MODES:
            target = Path(self.tmp.name) / f"manual-{mode[0]}"
            checker.prepare_checkout(self.root, target, mode)
            for relative, expected in manual_docs.generated(target).items():
                self.assertEqual((target / relative).read_bytes(), expected, relative)
                self.assertNotIn(b"\r", expected, relative)
        attributes = (checker.ROOT / ".gitattributes").read_bytes()
        for number, (rule, relative) in enumerate((
            (b"docs/gui/*.md text eol=lf", "docs/gui/rna-read-panel.md"),
            (b"docs/cli/*.md text eol=lf", "docs/cli/rna-read-panel.md"),
            (b"src/app/help_manual_sources.rs text eol=lf", manual_docs.EMBEDDED),
        )):
            unprotected = b"\n".join(line for line in attributes.split(b"\n") if line != rule)
            self.assertNotEqual(unprotected, attributes)
            target = Path(self.tmp.name) / f"manual-unprotected-{number}"
            prepare_unprotected_checkout(self.root, target, unprotected)
            generated = manual_docs.generated(target)
            self.assertNotEqual((target / relative).read_bytes(), generated[relative], relative)

    def test_repository_lf_default_preserves_source_binary_and_verbatim_bytes(self):
        # Synthetic literal fixtures exercise ordinary, previously unpinned
        # source paths and each verbatim exception from the real Git policy.
        text = b"first\nsecond\n"
        payloads = {path: text for path in (
            "src/unpinned.rs", "scripts/unpinned.sh", "docs/unpinned.md", "config.toml",
        )}
        for relative in (
            "test_files/pGEX_3X.fa",
            "test_files/fixtures/mapping/ensembl_human_tp63_all.fasta",
            "docs/figures/tert_upstream_early_coding_llr_background_tail_log10_cross_strand_spearman.tsv",
            "docs/figures/tert_upstream_early_coding_llr_background_tail_log10_cross_strand_spearman_tp73.tsv",
        ):
            payloads[relative] = (checker.ROOT / relative).read_bytes()
            self.assertNotIn(b"\r", payloads[relative], relative)
        payloads["binary.bin"] = b"\x00\xff\n\r\n"
        # Even a PDF without NULs must respect the explicit binary override.
        payloads["assets/probe.pdf"] = b"%PDF-synthetic\r\nfirst\nsecond\r\n"
        for relative in (
            "docs/audits/glen_followups_20261008/candidate_df5c8244/glen/original.log",
            "docs/audits/glen_followups_20261008/candidate_df5c8244/native-x11-supplement/original.log",
            "docs/audits/glen_followups_20261008/candidate_de9df8af/original.json",
        ):
            payloads[relative] = b"first\r\nsecond\nthird\r\n"
        (self.root / ".gitattributes").write_bytes((checker.ROOT / ".gitattributes").read_bytes())
        for relative, payload in payloads.items():
            path = self.root / relative
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(payload)
        checker.git(self.root, "add", "--all")
        checker.git(self.root, "-c", "commit.gpgsign=false", "commit", "--quiet",
                    "-m", "Synthetic repository LF and verbatim fixtures")
        for mode in checker.MODES:
            with self.subTest(mode=mode[0]):
                target = Path(self.tmp.name) / f"repository-{mode[0]}"
                checker.prepare_checkout(self.root, target, mode)
                for relative, payload in payloads.items():
                    self.assertEqual((target / relative).read_bytes(), payload, relative)
                    self.assertEqual(checker.git(target, "show", f"HEAD:{relative}"), payload, relative)
                self.assertEqual(checker.git(target, "status", "--porcelain"), b"")
                # Editors can still save CRLF: Git must normalize newly added
                # source in either ambient configuration, without touching data.
                (target / "src/added.rs").write_bytes(text.replace(b"\n", b"\r\n"))
                checker.git(target, "add", "--", "src/added.rs")
                self.assertEqual(checker.git(target, "show", ":src/added.rs"), text)
        broken = Path(self.tmp.name) / "repository-without-default"
        prepare_unprotected_checkout(self.root, broken,
                                     (checker.ROOT / ".gitattributes").read_bytes())
        self.assertEqual((broken / "src/unpinned.rs").read_bytes(),
                         text.replace(b"\n", b"\r\n"))

    def test_capability_maturity_survives_real_checkouts_and_detects_missing_rules(self):
        populate_maturity_fixture(self.root)
        attributes = (checker.ROOT / ".gitattributes").read_bytes()
        (self.root / ".gitattributes").write_bytes(attributes)
        checker.git(self.root, "add", "--all")
        checker.git(self.root, "-c", "commit.gpgsign=false", "commit", "--quiet", "-m", "public maturity inputs")
        for mode in checker.MODES:
            with self.subTest(mode=mode[0]):
                target = Path(self.tmp.name) / f"maturity-{mode[0]}"
                checker.prepare_checkout(self.root, target, mode)
                for relative, expected in maturity_generator.generate(target).items():
                    self.assertEqual((target / relative).read_bytes(), expected, relative)
                info = maturity_generator.contract.skill_info(target / maturity_generator.SKILL_DIR / "gentle-cloning")
                self.assertEqual(info["status"], "available")
                self.assertEqual(info["runtime_binding"], "not_checked")
        for index, relative in enumerate((maturity_generator.OUTPUT_PATH,
                "integrations/clawbio/skills/gentle-cloning/gentle_cloning.py")):
            unprotected = b"\n".join(line for line in attributes.split(b"\n")
                                       if not line.startswith(relative.encode() + b" "))
            self.assertNotEqual(unprotected, attributes)
            target = Path(self.tmp.name) / f"maturity-unprotected-{index}"
            prepare_unprotected_checkout(self.root, target, unprotected)
            self.assertIn(b"\r\n", (target / relative).read_bytes())
            outputs = maturity_generator.generate(target)
            self.assertTrue(any((target / path).read_bytes() != expected
                                for path, expected in outputs.items()))

    def test_missing_lf_rule_changes_digest_and_explicit_overlay_restores_it(self):
        (self.root / ".gitattributes").write_bytes(b"# No protection\n")
        checker.git(self.root, "add", "--", ".gitattributes")
        checker.git(self.root, "-c", "commit.gpgsign=false", "commit", "--quiet", "-m", "regression")
        broken = Path(self.tmp.name) / "broken"
        checker.prepare_checkout(self.root, broken, checker.MODES[1])
        self.assertNotEqual(hashlib.sha256((broken / "evidence.json").read_bytes()).digest(),
                            hashlib.sha256(self.payload).digest())
        fixed = Path(self.tmp.name) / "fixed"
        checker.prepare_checkout(self.root, fixed, checker.MODES[1],
                                 b"evidence.json text eol=lf\n")
        self.assertEqual((fixed / "evidence.json").read_bytes(), self.payload)
        self.assertEqual((self.root / ".gitattributes").read_bytes(), b"# No protection\n")

    def test_rnapkin_lock_receipt_bytes_survive_both_checkout_modes(self):
        relative = Path("docker/rnapkin/Cargo.lock")
        payload = (checker.ROOT / relative).read_bytes()
        digest = hashlib.sha256(payload).hexdigest()
        attributes = (checker.ROOT / ".gitattributes").read_bytes()
        (self.root / ".gitattributes").write_bytes(attributes)
        destination = self.root / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_bytes(payload)
        checker.git(self.root, "add", "--all")
        checker.git(self.root, "-c", "commit.gpgsign=false", "commit", "--quiet",
                    "-m", "retained RNAPKIN helper lock")
        unprotected = b"\n".join(line for line in attributes.split(b"\n")
                                   if not line.startswith(b"Cargo.lock "))
        broken = Path(self.tmp.name) / "rnapkin-lock-unprotected"
        prepare_unprotected_checkout(self.root, broken, unprotected)
        self.assertEqual((broken / relative).read_bytes(), payload.replace(b"\n", b"\r\n"))
        self.assertNotEqual(hashlib.sha256((broken / relative).read_bytes()).hexdigest(), digest)
        for mode in checker.MODES:
            with self.subTest(mode=mode[0]):
                target = Path(self.tmp.name) / f"rnapkin-lock-{mode[0]}"
                checker.prepare_checkout(self.root, target, mode)
                self.assertEqual((target / relative).read_bytes(), payload,
                                 "Cargo.lock text eol=lf must also protect the helper lock")
                self.assertEqual(hashlib.sha256((target / relative).read_bytes()).hexdigest(), digest)

    def test_generated_tutorial_json_stays_byte_exact_in_both_checkout_modes(self):
        # Use the real checkout policy with synthetic generated JSON, so removing
        # a targeted LF rule reproduces Windows' strict drift-check failure.
        (self.root / ".gitattributes").write_bytes(
            (checker.ROOT / ".gitattributes").read_bytes())
        paths = ("docs/tutorial/catalog.json", "docs/tutorial/manifest.json")
        for relative in paths:
            path = self.root / relative
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(self.payload)
        checker.git(self.root, "add", "--all")
        checker.git(self.root, "-c", "commit.gpgsign=false", "commit", "--quiet", "-m", "generated JSON")
        for mode in checker.MODES:
            target = Path(self.tmp.name) / f"generated-{mode[0]}"
            checker.prepare_checkout(self.root, target, mode)
            for relative in paths:
                with self.subTest(mode=mode[0], path=relative):
                    self.assertEqual((target / relative).read_bytes(), self.payload)

    def test_manual_tutorial_review_identity_survives_lf_and_crlf_checkouts(self):
        sources = ("07-04_microarray_evidence_gui_cli.json",
                   "10-03_measured_gel_sizing_gui_cli.json",
                   "08-18_repeat_environment_cohorts.json",
                   "09-04_cloning_time_cost_and_protocol_suppliers.json",
                   "10-04_targeted_nanopore_fastq_panel.json")
        paths = [Path("docs/tutorial/review_manifest.json"), Path("docs/tutorial/catalog.json")]
        paths.extend(Path("docs/tutorial/sources") / name for name in sources)
        payloads = {path: (checker.ROOT / path).read_bytes() for path in paths}
        (self.root / ".gitattributes").write_bytes((checker.ROOT / ".gitattributes").read_bytes())
        for path, payload in payloads.items():
            target = self.root / path
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(payload)
        checker.git(self.root, "add", "--all")
        checker.git(self.root, "-c", "commit.gpgsign=false", "commit", "--quiet",
                    "-m", "manual tutorial review identity")
        for mode in checker.MODES:
            target = Path(self.tmp.name) / f"manual-review-{mode[0]}"
            checker.prepare_checkout(self.root, target, mode)
            reviews = {entry["tutorial_id"]: entry for entry in
                       json.loads((target / paths[0]).read_bytes())["entries"]}
            catalog = {entry["id"]: entry for entry in
                       json.loads((target / paths[1]).read_bytes())["entries"]}
            known = {json.loads((target / path).read_bytes())["id"] for path in paths[2:]}
            self.assertEqual(known - set(reviews), set(), "missing entries add report warnings")
            self.assertEqual((target / paths[1]).read_bytes(), payloads[paths[1]])
            canonical = (json.dumps(json.loads(payloads[paths[1]]),
                                    ensure_ascii=False, indent=2) + "\n").encode()
            self.assertEqual((target / paths[1]).read_bytes(), canonical,
                             "catalog formatting is part of generator freshness")
            for identifier in known:
                with self.subTest(mode=mode[0], tutorial=identifier):
                    self.assertEqual(reviews[identifier]["tutorial_kind"], "guided_walkthrough")
                    self.assertEqual(reviews[identifier]["tutorial_status"], "active")
                    self.assertNotEqual(catalog[identifier]["review_status"], "missing_review_manifest_entry")
            missing = dict(reviews)
            del missing["microarray_evidence_gui_cli"]
            self.assertEqual(known - set(missing), {"microarray_evidence_gui_cli"},
                             "the original missing-entry boundary is detectable in either checkout")

    def test_release_tutorial_registration_survives_lf_and_crlf_checkouts(self):
        from scripts.test_release_workflow_tutorials import STEMS, assert_manual_registration

        # Real discovery bytes plus hand-authored sources/guides: registration
        # is not a generated runtime chapter or a fabricated review sign-off.
        paths = [Path("docs/tutorial") / name for name in
                 ("catalog.json", "review_manifest.json", "manifest.json", "README.md", "sources/catalog_meta.json")]
        for stem in STEMS:
            paths.extend((Path(f"docs/tutorial/sources/{stem}.json"), Path(f"docs/tutorial/{stem}.md")))
        payloads = {path: (checker.ROOT / path).read_bytes() for path in paths}
        (self.root / ".gitattributes").write_bytes((checker.ROOT / ".gitattributes").read_bytes())
        for path, payload in payloads.items():
            target = self.root / path
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(payload)
        checker.git(self.root, "add", "--all")
        checker.git(self.root, "-c", "commit.gpgsign=false", "commit", "--quiet",
                    "-m", "release tutorial registration")
        for mode in checker.MODES:
            target = Path(self.tmp.name) / f"release-tutorial-{mode[0]}"
            checker.prepare_checkout(self.root, target, mode)
            self.assertEqual((target / paths[0]).read_bytes(), payloads[paths[0]])
            assert_manual_registration(self, target)
            catalog = json.loads((target / paths[0]).read_bytes())
            catalog["entries"] = [row for row in catalog["entries"]
                                  if row["id"] != "gene_set_storage_projection"]
            (target / paths[0]).write_text(json.dumps(catalog), encoding="utf-8")
            with self.assertRaises(AssertionError):
                # No subTest here: failures must propagate for the negative control.
                assert_manual_registration(unittest.TestCase(), target)

    def test_measured_gel_image_requests_and_checksums_survive_lf_crlf_checkouts(self):
        relative = Path("docs/tutorial/inputs/gel_sizing")
        names = ["gel_sizing_synthetic.png", "dna-analysis.json", "protein-analysis.json", "SHA256SUMS"]
        attributes = (checker.ROOT / ".gitattributes").read_bytes()
        (self.root / ".gitattributes").write_bytes(attributes)
        payloads = {}
        for name in names:
            path = relative / name
            payloads[path] = (checker.ROOT / path).read_bytes()
            destination = self.root / path
            destination.parent.mkdir(parents=True, exist_ok=True)
            destination.write_bytes(payloads[path])
        checker.git(self.root, "add", "--all")
        checker.git(self.root, "-c", "commit.gpgsign=false", "commit", "--quiet", "-m", "drawn gel fixture")
        for mode in checker.MODES:
            target = Path(self.tmp.name) / f"gel-{mode[0]}"
            checker.prepare_checkout(self.root, target, mode)
            for path, payload in payloads.items():
                self.assertEqual((target / path).read_bytes(), payload)
            for line in (target / relative / "SHA256SUMS").read_text().splitlines():
                digest, name = line.split("  ", 1)
                self.assertEqual(hashlib.sha256((target / relative / name).read_bytes()).hexdigest(), digest)
            request = json.loads((target / relative / "dna-analysis.json").read_bytes())
            self.assertEqual(request["image_sha256"], hashlib.sha256((target / relative / names[0]).read_bytes()).hexdigest())
        unprotected = b"\n".join(line for line in attributes.split(b"\n")
                                 if not line.startswith(b"docs/tutorial/inputs/gel_sizing/*.json "))
        broken = Path(self.tmp.name) / "gel-unprotected"
        prepare_unprotected_checkout(self.root, broken, unprotected)
        self.assertNotEqual((broken / relative / "dna-analysis.json").read_bytes(), payloads[relative / "dna-analysis.json"])

    def test_public_tapestation_inputs_preserve_pdf_png_and_coordinate_hashes(self):
        relative = Path("test_files/fixtures/gel_images/tapestation_26r266")
        attributes = (checker.ROOT / ".gitattributes").read_bytes()
        (self.root / ".gitattributes").write_bytes(attributes)
        names = ("before.pdf", "after.pdf", "before.png", "after.png", "manifest.json",
                 "region_summaries.json", "SHA256SUMS")
        payloads = {}
        for name in names:
            path = relative / name
            payloads[path] = (checker.ROOT / path).read_bytes()
            destination = self.root / path
            destination.parent.mkdir(parents=True, exist_ok=True)
            destination.write_bytes(payloads[path])
        checker.git(self.root, "add", "--all")
        checker.git(self.root, "-c", "commit.gpgsign=false", "commit", "--quiet", "-m", "public real gel evidence")
        for mode in checker.MODES:
            target = Path(self.tmp.name) / f"tapestation-{mode[0]}"
            checker.prepare_checkout(self.root, target, mode)
            for path, payload in payloads.items():
                self.assertEqual((target / path).read_bytes(), payload)
            for line in (target / relative / "SHA256SUMS").read_text().splitlines():
                digest, name = line.split("  ", 1)
                self.assertEqual(hashlib.sha256((target / relative / name).read_bytes()).hexdigest(), digest)
            manifest = json.loads((target / relative / "manifest.json").read_bytes())
            for source in manifest["sources"]:
                for kind in ("image", "pdf"):
                    self.assertEqual(hashlib.sha256((target / relative / source[kind]).read_bytes()).hexdigest(),
                                     source[f"{kind}_sha256"])
        unprotected = b"\n".join(line for line in attributes.split(b"\n")
                                 if not line.startswith(b"test_files/fixtures/gel_images/tapestation_26r266/*.json "))
        broken = Path(self.tmp.name) / "tapestation-unprotected"
        prepare_unprotected_checkout(self.root, broken, unprotected)
        self.assertNotEqual((broken / relative / "manifest.json").read_bytes(), payloads[relative / "manifest.json"])

    def test_priming_return_and_plot_hashes_survive_both_checkout_modes(self):
        relative = Path("docs/figures/tapestation_priming")
        attributes = (checker.ROOT / ".gitattributes").read_bytes()
        (self.root / ".gitattributes").write_bytes(attributes)
        names = ("cell_1_custom_vs_combined.png", "cell_1_manufacturer_vs_combined.png",
                 "cell_2_custom_vs_combined.png", "cell_2_manufacturer_vs_combined.png",
                 "observations.json", "SHA256SUMS")
        for name in names:
            target = self.root / relative / name
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes((checker.ROOT / relative / name).read_bytes())
        checker.git(self.root, "add", "--all")
        checker.git(self.root, "-c", "commit.gpgsign=false", "commit", "--quiet", "-m", "priming return")
        for mode in checker.MODES:
            clone = Path(self.tmp.name) / f"priming-{mode[0]}"
            checker.prepare_checkout(self.root, clone, mode)
            for line in (clone / relative / "SHA256SUMS").read_text().splitlines():
                digest, name = line.split("  ", 1)
                self.assertEqual(hashlib.sha256((clone / relative / name).read_bytes()).hexdigest(), digest)
        unprotected = b"\n".join(line for line in attributes.split(b"\n")
                                 if not line.startswith(b"docs/figures/tapestation_priming/*.json "))
        broken = Path(self.tmp.name) / "priming-unprotected"
        prepare_unprotected_checkout(self.root, broken, unprotected)
        self.assertNotEqual((broken / relative / "observations.json").read_bytes(),
                            (checker.ROOT / relative / "observations.json").read_bytes())

    def test_repeat_lesson_raw_rmsk_hash_survives_both_checkout_modes(self):
        relative = Path("test_files/fixtures/repeat_cohort/local_project.rmsk.txt")
        payload = (checker.ROOT / relative).read_bytes()
        digest = "2fe38f29ab6ebdfb20e1ef3a0f5c28f4e01fffed5dfba32d0503c40af4722063"
        self.assertEqual(hashlib.sha256(payload).hexdigest(), digest)
        attributes = (checker.ROOT / ".gitattributes").read_bytes()
        (self.root / ".gitattributes").write_bytes(attributes)
        destination = self.root / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_bytes(payload)
        checker.git(self.root, "add", "--all")
        checker.git(self.root, "-c", "commit.gpgsign=false", "commit", "--quiet", "-m", "synthetic repeat rows")
        for mode in checker.MODES:
            target = Path(self.tmp.name) / f"repeat-{mode[0]}"
            checker.prepare_checkout(self.root, target, mode)
            self.assertEqual((target / relative).read_bytes(), payload)
            self.assertEqual(hashlib.sha256((target / relative).read_bytes()).hexdigest(), digest)
        unprotected = b"\n".join(line for line in attributes.split(b"\n")
                                 if not line.startswith(b"test_files/fixtures/repeat_cohort/*.txt "))
        broken = Path(self.tmp.name) / "repeat-unprotected"
        prepare_unprotected_checkout(self.root, broken, unprotected)
        self.assertNotEqual(hashlib.sha256((broken / relative).read_bytes()).hexdigest(), digest)

    def test_retained_vkorc1_wip_bytes_and_hashes_survive_both_checkout_modes(self):
        relative = Path("docs/tutorial/reproducibility/vkorc1_rs9923231_promoter_reporter")
        files = ["variant_promoter_context.json", "promoter_reporter_candidates.json",
                 "vkorc1_rs9923231_promoter_context.svg",
                 "vkorc1_rs9923231_reporter_reference.svg",
                 "vkorc1_rs9923231_reporter_alternate.svg"]
        provenance = (checker.ROOT / relative / "README.md").read_text(encoding="utf-8")
        self.assertIn("historical, unverified WIP", provenance)
        (self.root / ".gitattributes").write_bytes((checker.ROOT / ".gitattributes").read_bytes())
        for name in files:
            destination = self.root / relative / name
            destination.parent.mkdir(parents=True, exist_ok=True)
            payload = (checker.ROOT / relative / name).read_bytes()
            self.assertIn(hashlib.sha256(payload).hexdigest(), provenance)
            destination.write_bytes(payload)
        checker.git(self.root, "add", "--all")
        checker.git(self.root, "-c", "commit.gpgsign=false", "commit", "--quiet",
                    "-m", "retained historical VKORC1 previews")
        for mode in checker.MODES:
            target = Path(self.tmp.name) / f"vkorc1-{mode[0]}"
            checker.prepare_checkout(self.root, target, mode)
            for name in files:
                with self.subTest(mode=mode[0], file=name):
                    self.assertEqual((target / relative / name).read_bytes(),
                                     (checker.ROOT / relative / name).read_bytes())

    def test_vkorc1_guard_and_generated_chapter_survive_both_checkout_modes(self):
        fixture = Path("docs/examples/assets/allele_pair_guard/multiallelic.gb")
        chapter = Path("docs/tutorial/generated/chapters/08-04_vkorc1_warfarin_promoter_luciferase_gui.md")
        hub = Path("docs/tutorial/generated/README.md")
        payloads = {
            fixture: (checker.ROOT / fixture).read_bytes(),
            chapter: (checker.ROOT / chapter).read_bytes(),
            hub: (checker.ROOT / hub).read_bytes(),
        }
        ledger = json.loads((checker.ROOT / "docs/tutorial/generated/report.json").read_bytes())
        for path in (chapter, hub):
            name = path.relative_to("docs/tutorial/generated").as_posix()
            self.assertEqual(hashlib.sha256(payloads[path]).hexdigest(), ledger["file_checksums"][name])
        attributes = (checker.ROOT / ".gitattributes").read_bytes()
        (self.root / ".gitattributes").write_bytes(attributes)
        for relative, payload in payloads.items():
            destination = self.root / relative
            destination.parent.mkdir(parents=True, exist_ok=True)
            destination.write_bytes(payload)
        checker.git(self.root, "add", "--all")
        checker.git(self.root, "-c", "commit.gpgsign=false", "commit", "--quiet",
                    "-m", "synthetic 08.04 acceptance byte boundaries")
        unprotected = b"\n".join(
            line for line in attributes.split(b"\n")
            if not any(line.startswith(path.as_posix().encode() + b" ") for path in payloads)
        )
        broken = Path(self.tmp.name) / "vkorc1-guard-unprotected"
        prepare_unprotected_checkout(self.root, broken, unprotected)
        for relative, payload in payloads.items():
            self.assertNotEqual((broken / relative).read_bytes(), payload)
        for mode in checker.MODES:
            target = Path(self.tmp.name) / f"vkorc1-guard-{mode[0]}"
            checker.prepare_checkout(self.root, target, mode)
            for relative, payload in payloads.items():
                with self.subTest(mode=mode[0], path=relative):
                    self.assertEqual((target / relative).read_bytes(), payload)

    def test_glen_20261008_audit_report_bytes_survive_both_checkout_modes(self):
        relative = Path("docs/glen_tutorial_parity_gui_audit_20261008.md")
        expected = "186e75a987fb50cb50ed883fa2b6fa5434718fcb553ceec7aad90a8af7c1cfce"
        payload = (checker.ROOT / relative).read_bytes()
        self.assertEqual(hashlib.sha256(payload).hexdigest(), expected)
        attributes = (checker.ROOT / ".gitattributes").read_bytes()
        (self.root / ".gitattributes").write_bytes(attributes)
        destination = self.root / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_bytes(payload)
        checker.git(self.root, "add", "--all")
        checker.git(self.root, "-c", "commit.gpgsign=false", "commit", "--quiet",
                    "-m", "owner-supplied exact-candidate audit report")
        unprotected = b"\n".join(line for line in attributes.split(b"\n")
                                  if not line.startswith(relative.as_posix().encode() + b" "))
        broken = Path(self.tmp.name) / "glen-audit-unprotected"
        prepare_unprotected_checkout(self.root, broken, unprotected)
        self.assertNotEqual(hashlib.sha256((broken / relative).read_bytes()).hexdigest(), expected)
        for mode in checker.MODES:
            with self.subTest(mode=mode[0]):
                target = Path(self.tmp.name) / f"glen-audit-{mode[0]}"
                checker.prepare_checkout(self.root, target, mode)
                self.assertEqual((target / relative).read_bytes(), payload)
                self.assertEqual(hashlib.sha256((target / relative).read_bytes()).hexdigest(), expected)

    def test_f35a58e6_original_report_and_manifest_survive_both_checkout_modes(self):
        base = Path("docs/audits/glen_followups_20261008/candidate_f35a58e6")
        profile = Path("docs/audits/glen_followups_20261008/post_blast_profile_f35a58e6")
        decision = Path("docs/audits/glen_followups_20261008/flna_teaching_exon_f35a58e6")
        expected = {
            base / "REPORT.md": "63ab711e81383b407d8b4a00b90c2873553551646944696bc127f9d27ec7bc52",
            base / "SHA256SUMS": "cabb5e49d9e904d55cf0d07b9aae0879ba269c2a3ca4d9899a199e7df189d9e0",
            profile / "REPORT.md": "aef052f65bf668e8e2729b3f19ea95b05227bf3d377e6c81e03b66e7d09f795c",
            profile / "SHA256SUMS": "466e5371a53c18e856d8d28defd0a3210fd5e73265023ff01c18110afae537ff",
            decision / "decision-record.json": "5202f0605fca8cdd24bbd00e5066da855ffe41ff9aa6695ad09c8f75da4972be",
            decision / "bundle-manifest.json": "f1966fdf183c637ac571137593c6de6b926b44790c81ae23fbdbfde45516974f",
            decision / "version-label-addendum.json": "7ce94e64dbb6629dc0dbe3a6baeef6c3a567e5ca75c1a5766f145fa597f730a5",
            decision / "v2-bundle-manifest.json": "444dd108f12b9e48cf026cabfdc35846574163ac4f16908b33597e6fd598f8f0",
        }
        payloads = {path: (checker.ROOT / path).read_bytes() for path in expected}
        addendum = json.loads(payloads[decision / "version-label-addendum.json"])
        manifest = json.loads(payloads[decision / "v2-bundle-manifest.json"])
        self.assertEqual(addendum["candidate_sha"], manifest["candidate_sha"])
        self.assertEqual(addendum["original_decision_sha256"],
                         expected[decision / "decision-record.json"])
        self.assertEqual(manifest["preserves_original_decision_sha256"],
                         expected[decision / "decision-record.json"])
        self.assertEqual(addendum["loaded_annotation"]["transcript"], "ENST00000360319.10")
        self.assertEqual(addendum["retained_uniprot_cross_reference"]["transcript"],
                         "ENST00000360319.9")
        self.assertEqual(addendum["effect_on_human_choice"], "none")
        artifacts = {item["path"]: item for item in manifest["artifacts"]}
        for filename in ("decision-record.json", "version-label-addendum.json"):
            self.assertEqual(artifacts[filename]["sha256"], "sha256:" + expected[decision / filename])
            self.assertEqual(artifacts[filename]["bytes"], len(payloads[decision / filename]))
        attributes = (checker.ROOT / ".gitattributes").read_bytes()
        (self.root / ".gitattributes").write_bytes(attributes)
        for relative, payload in payloads.items():
            self.assertEqual(hashlib.sha256(payload).hexdigest(), expected[relative])
            destination = self.root / relative
            destination.parent.mkdir(parents=True, exist_ok=True)
            destination.write_bytes(payload)
        checker.git(self.root, "add", "--all")
        checker.git(self.root, "-c", "commit.gpgsign=false", "commit", "--quiet",
                    "-m", "owner-supplied frozen-main report and manifest")
        rules = tuple(path.as_posix().encode() + b" " for path in expected)
        unprotected = b"\n".join(line for line in attributes.split(b"\n")
                                   if not line.startswith(rules))
        broken = Path(self.tmp.name) / "f35a58e6-unprotected"
        prepare_unprotected_checkout(self.root, broken, unprotected)
        for relative, digest in expected.items():
            self.assertNotEqual(hashlib.sha256((broken / relative).read_bytes()).hexdigest(), digest)
        for mode in checker.MODES:
            with self.subTest(mode=mode[0]):
                target = Path(self.tmp.name) / f"f35a58e6-{mode[0]}"
                checker.prepare_checkout(self.root, target, mode)
                for relative, payload in payloads.items():
                    self.assertEqual((target / relative).read_bytes(), payload)
                    self.assertEqual(hashlib.sha256((target / relative).read_bytes()).hexdigest(),
                                     expected[relative])

    def test_original_audit_verification_receipt_survives_both_checkout_modes(self):
        relative = Path("docs/audits/glen_followups_20261008/original_archive_30f23aa4/verification.json")
        expected = "1a41c227b6e91587e6d310114d6672d875008f51170bf75428da303b3b50dbb5"
        payload = (checker.ROOT / relative).read_bytes()
        self.assertEqual(hashlib.sha256(payload).hexdigest(), expected)
        attributes = (checker.ROOT / ".gitattributes").read_bytes()
        (self.root / ".gitattributes").write_bytes(attributes)
        destination = self.root / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_bytes(payload)
        checker.git(self.root, "add", "--all")
        checker.git(self.root, "-c", "commit.gpgsign=false", "commit", "--quiet",
                    "-m", "original audit archive verification receipt")
        unprotected = b"\n".join(line for line in attributes.split(b"\n")
                                  if not line.startswith(relative.as_posix().encode() + b" "))
        broken = Path(self.tmp.name) / "original-audit-receipt-unprotected"
        prepare_unprotected_checkout(self.root, broken, unprotected)
        self.assertNotEqual(hashlib.sha256((broken / relative).read_bytes()).hexdigest(), expected)
        for mode in checker.MODES:
            with self.subTest(mode=mode[0]):
                target = Path(self.tmp.name) / f"original-audit-receipt-{mode[0]}"
                checker.prepare_checkout(self.root, target, mode)
                self.assertEqual((target / relative).read_bytes(), payload)

    def test_frozen_024ed6c0_records_survive_lf_crlf_and_detect_missing_policy(self):
        from scripts.test_frozen_024ed6_audit import BASE, retained_hashes

        expected = retained_hashes()
        expected[BASE / "SHA256SUMS"] = hashlib.sha256((checker.ROOT / BASE / "SHA256SUMS").read_bytes()).hexdigest()
        attributes = (checker.ROOT / ".gitattributes").read_bytes()
        (self.root / ".gitattributes").write_bytes(attributes)
        payloads = {path: (checker.ROOT / path).read_bytes() for path in expected}
        for path, payload in payloads.items():
            self.assertEqual(hashlib.sha256(payload).hexdigest(), expected[path])
            destination = self.root / path
            destination.parent.mkdir(parents=True, exist_ok=True)
            destination.write_bytes(payload)
        checker.git(self.root, "add", "--all")
        checker.git(self.root, "-c", "commit.gpgsign=false", "commit", "--quiet",
                    "-m", "frozen small audit records")
        unprotected = b"\n".join(line for line in attributes.split(b"\n")
                                  if not line.startswith(BASE.as_posix().encode() + b"/"))
        broken = Path(self.tmp.name) / "frozen-audit-unprotected"
        prepare_unprotected_checkout(self.root, broken, unprotected)
        for path in expected:
            self.assertNotEqual(hashlib.sha256((broken / path).read_bytes()).hexdigest(), expected[path])
        for mode in checker.MODES:
            target = Path(self.tmp.name) / f"frozen-audit-{mode[0]}"
            checker.prepare_checkout(self.root, target, mode)
            for path, payload in payloads.items():
                with self.subTest(mode=mode[0], path=path):
                    self.assertEqual((target / path).read_bytes(), payload)

    def test_frozen_df5c8244_records_preserve_mixed_bytes_in_both_checkout_modes(self):
        from scripts.test_frozen_df5c8244_audit import BASE, GLEN, NATIVE, retained_hashes

        expected = retained_hashes()
        payloads = {path: (checker.ROOT / path).read_bytes() for path in expected}
        attributes = (checker.ROOT / ".gitattributes").read_bytes()
        (self.root / ".gitattributes").write_bytes(attributes)
        for path, payload in payloads.items():
            self.assertEqual(hashlib.sha256(payload).hexdigest(), expected[path])
            destination = self.root / path
            destination.parent.mkdir(parents=True, exist_ok=True)
            destination.write_bytes(payload)
        checker.git(self.root, "add", "--all")
        checker.git(self.root, "-c", "commit.gpgsign=false", "commit", "--quiet",
                    "-m", "frozen df5c8244 review and mixed-newline logs")
        checker.git(self.root, "diff", "--check", "HEAD^", "HEAD", "--",
                    GLEN.as_posix(), NATIVE.as_posix())
        source = "scripts/test_frozen_df5c8244_audit.py"
        self.assertEqual(checker.git(self.root, "check-attr", "whitespace", "--", source).strip(),
                         f"{source}: whitespace: unspecified".encode())
        unprotected = b"\n".join(line for line in attributes.split(b"\n")
                                  if not line.startswith(BASE.as_posix().encode() + b"/"))
        broken = Path(self.tmp.name) / "df5c8244-unprotected"
        prepare_unprotected_checkout(self.root, broken, unprotected)
        for path in (GLEN / "REPORT.md", NATIVE / "INDEX.md"):
            self.assertNotEqual((broken / path).read_bytes(), payloads[path])
        for mode in checker.MODES:
            target = Path(self.tmp.name) / f"df5c8244-{mode[0]}"
            checker.prepare_checkout(self.root, target, mode)
            for path, payload in payloads.items():
                with self.subTest(mode=mode[0], path=path):
                    self.assertEqual((target / path).read_bytes(), payload)

    def test_terminal_de9df8af_ci_records_survive_lf_crlf_and_missing_policy_control(self):
        from scripts.test_frozen_de9df8af_ci import BASE, retained_hashes

        expected = retained_hashes()
        payloads = {path: (checker.ROOT / path).read_bytes() for path in expected}
        attributes = (checker.ROOT / ".gitattributes").read_bytes()
        (self.root / ".gitattributes").write_bytes(attributes)
        for path, payload in payloads.items():
            self.assertEqual(hashlib.sha256(payload).hexdigest(), expected[path])
            destination = self.root / path
            destination.parent.mkdir(parents=True, exist_ok=True)
            destination.write_bytes(payload)
        checker.git(self.root, "add", "--all")
        checker.git(self.root, "-c", "commit.gpgsign=false", "commit", "--quiet",
                    "-m", "terminal CI failure and stripped-child observations")
        unprotected = b"\n".join(line for line in attributes.split(b"\n")
                                  if not line.startswith(BASE.as_posix().encode() + b"/"))
        broken = Path(self.tmp.name) / "de9df8af-unprotected"
        prepare_unprotected_checkout(self.root, broken, unprotected)
        probe = BASE / "default_stack_from_logs.json"
        self.assertIn(b"\n", payloads[probe])
        self.assertNotEqual((broken / probe).read_bytes(), payloads[probe])
        for name in ("run.api-original.json", "jobs.api-original.json"):
            path = BASE / name
            self.assertNotIn(b"\n", payloads[path])
            self.assertEqual((broken / path).read_bytes(), payloads[path])
        for mode in checker.MODES:
            target = Path(self.tmp.name) / f"de9df8af-{mode[0]}"
            checker.prepare_checkout(self.root, target, mode)
            for path, payload in payloads.items():
                with self.subTest(mode=mode[0], path=path):
                    self.assertEqual((target / path).read_bytes(), payload)
                    attribute = checker.git(target, "check-attr", "text", "--", path.as_posix())
                    self.assertEqual(attribute.strip(),
                                     f"{path.as_posix()}: text: unset".encode())

    def test_nonpublishing_container_receipt_survives_lf_crlf_and_detects_missing_policy(self):
        relative = Path("docs/audits/glen_followups_20261008/container_validate_5e6201a6/container-build.json")
        expected = "682d7d10df1b907dd80121c7055c62b3335952147fc7a3bf41988e7f13ffc42d"
        payload = (checker.ROOT / relative).read_bytes()
        self.assertEqual(hashlib.sha256(payload).hexdigest(), expected)
        attributes = (checker.ROOT / ".gitattributes").read_bytes()
        (self.root / ".gitattributes").write_bytes(attributes)
        destination = self.root / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_bytes(payload)
        checker.git(self.root, "add", "--all")
        checker.git(self.root, "-c", "commit.gpgsign=false", "commit", "--quiet",
                    "-m", "original nonpublishing container receipt")
        unprotected = b"\n".join(line for line in attributes.split(b"\n")
                                  if not line.startswith(relative.as_posix().encode() + b" "))
        broken = Path(self.tmp.name) / "container-receipt-unprotected"
        prepare_unprotected_checkout(self.root, broken, unprotected)
        self.assertNotEqual(hashlib.sha256((broken / relative).read_bytes()).hexdigest(), expected)
        for mode in checker.MODES:
            with self.subTest(mode=mode[0]):
                target = Path(self.tmp.name) / f"container-receipt-{mode[0]}"
                checker.prepare_checkout(self.root, target, mode)
                self.assertEqual((target / relative).read_bytes(), payload)

    def test_frozen_6d8b4db7_audit_receipts_survive_both_checkout_modes(self):
        base = Path("docs/audits/glen_followups_20261008/candidate_6d8b4db7")
        expected = {
            base / "native_gui.json": "c0b0c89aaa8a2989e0c2d86de9f8327fe5e9190a4ce6634cb7866d91460efc2e",
            base / "scroll_cpu_smoke.json": "420e2c280ccb8294ecab5046672c85ffe6461104275f050fc43fdf65acaf6bc6",
            base / "agent_cli_build.json": "ff026cb592f64f8a81a472e144ad5bdd539bc09b2c0f905058fe0489a661acd1",
            base / "public.json": "d65c02e64fb429ea41662a9d5dd5336a994304780a27a0390645ec2eecedb30d",
            base / "base_comparison.json": "679e52fd00d98552c28d3dd3413c6b38e0f0fa0e93302d7fe1be754a425f434d",
            base / "base_comparison.svg": "1eee85466468d971cf88d0aad2b29841819e613755f60a81509cde42455b7bc6",
            base / "live_agent.json": "aa5a3d8e68a349c55f5764152d5e0613150c0f1a63d05ae11ca15d22d9a1ed25",
        }
        attributes = (checker.ROOT / ".gitattributes").read_bytes()
        (self.root / ".gitattributes").write_bytes(attributes)
        payloads = {path: (checker.ROOT / path).read_bytes() for path in expected}
        for relative, payload in payloads.items():
            self.assertEqual(hashlib.sha256(payload).hexdigest(), expected[relative])
            destination = self.root / relative
            destination.parent.mkdir(parents=True, exist_ok=True)
            destination.write_bytes(payload)
        checker.git(self.root, "add", "--all")
        checker.git(self.root, "-c", "commit.gpgsign=false", "commit", "--quiet",
                    "-m", "exact-candidate native, scroll, public and live evidence")
        for relative in payloads:
            unprotected = b"\n".join(line for line in attributes.split(b"\n")
                                      if not line.startswith(relative.as_posix().encode() + b" "))
            broken = Path(self.tmp.name) / f"candidate-receipt-{relative.name}-unprotected"
            prepare_unprotected_checkout(self.root, broken, unprotected)
            self.assertNotEqual(hashlib.sha256((broken / relative).read_bytes()).hexdigest(),
                                expected[relative])
        for mode in checker.MODES:
            target = Path(self.tmp.name) / f"candidate-receipts-{mode[0]}"
            checker.prepare_checkout(self.root, target, mode)
            for relative, payload in payloads.items():
                with self.subTest(mode=mode[0], path=relative):
                    self.assertEqual((target / relative).read_bytes(), payload)

    def test_sequence_design_tutorial_semantics_survive_both_checkout_modes(self):
        from scripts.test_sequence_design_tutorial import fenced_request

        request_path = Path("docs/tutorial/inputs/synthetic_sequence_design.json")
        guide_path = Path("docs/tutorial/06-07_synthetic_sequence_design.md")
        source_path = Path("docs/tutorial/sources/06-07_synthetic_sequence_design.json")
        (self.root / ".gitattributes").write_bytes((checker.ROOT / ".gitattributes").read_bytes())
        for relative in (request_path, guide_path, source_path):
            destination = self.root / relative
            destination.parent.mkdir(parents=True, exist_ok=True)
            destination.write_bytes((checker.ROOT / relative).read_bytes())
        checker.git(self.root, "add", "--all")
        checker.git(self.root, "-c", "commit.gpgsign=false", "commit", "--quiet",
                    "-m", "synthetic sequence-design teaching request")
        request = json.loads((checker.ROOT / request_path).read_bytes())
        source = json.loads((checker.ROOT / source_path).read_bytes())
        for mode in checker.MODES:
            with self.subTest(mode=mode[0]):
                target = Path(self.tmp.name) / f"sequence-design-{mode[0]}"
                checker.prepare_checkout(self.root, target, mode)
                guide_bytes = (target / guide_path).read_bytes()
                self.assertNotIn(b"\r", guide_bytes,
                                 "Repository LF policy overrides the CRLF checkout mode")
                self.assertEqual(fenced_request(guide_bytes.decode("utf-8")), request)
                self.assertIn(b"candidate_store_unavailable", guide_bytes)
                self.assertIn(b"not an empty", guide_bytes)
                self.assertIn(b"**reopen**", guide_bytes)
                self.assertEqual(json.loads((target / request_path).read_bytes()), request)
                self.assertEqual(json.loads((target / source_path).read_bytes()), source)

    def test_frozen_3ac4a02d_audit_receipts_survive_both_checkout_modes(self):
        base = Path("docs/audits/glen_followups_20261008/candidate_3ac4a02d")
        expected = {
            base / "scroll_cpu_smoke.json": "87a274adaf3c8046f4458b8d0374757ea7db3cc63fab2a5eed7e9dbc5094c00d",
            base / "scroll_verification.json": "58a16df1a020904e7b7e7ea4573e1eb314f9aec48c253094cd6bf70281fbb3ab",
            base / "agent_cli_build.json": "8f28d2f10887acd283234d40a5440a27d55fb17e88d174bd57c3b3e39986c341",
            base / "native_gui.json": "bef0707566d9f68f19f51d57cf760d47c8fbef30818b2f05dcc22bf4c75a121d",
            base / "native_verification.json": "cf6b2500740d97071f7e010f53b3fdf97d6de68f1880d9dc3b1e0b4559c4e6f0",
            base / "public.json": "4ae522db4eac214ba7877c775eb1698750424f6932c4710dbb351c5eaf9bb9b9",
            base / "base_comparison.json": "5ced393927f191ccd78bb3f81a26597688a1286eaaefb2b5ee7ce30d92dc4313",
            base / "base_comparison.svg": "9aee9045ce6e107fd43c0c3f68946907f4e8ac7aaed275367507ca2b528b2c9d",
            base / "live_agent.json": "89f4fdba74d027d813fd5daecb9a15527cb6b9a9d3fb04b23110567c441c79d7",
        }
        attributes = (checker.ROOT / ".gitattributes").read_bytes()
        (self.root / ".gitattributes").write_bytes(attributes)
        payloads = {path: (checker.ROOT / path).read_bytes() for path in expected}
        for relative, payload in payloads.items():
            self.assertEqual(hashlib.sha256(payload).hexdigest(), expected[relative])
            destination = self.root / relative
            destination.parent.mkdir(parents=True, exist_ok=True)
            destination.write_bytes(payload)
        checker.git(self.root, "add", "--all")
        checker.git(self.root, "-c", "commit.gpgsign=false", "commit", "--quiet",
                    "-m", "repaired-candidate audit evidence")
        for relative in payloads:
            unprotected = b"\n".join(line for line in attributes.split(b"\n")
                                      if not line.startswith(relative.as_posix().encode() + b" "))
            broken = Path(self.tmp.name) / f"repaired-receipt-{relative.name}-unprotected"
            prepare_unprotected_checkout(self.root, broken, unprotected)
            self.assertNotEqual(hashlib.sha256((broken / relative).read_bytes()).hexdigest(),
                                expected[relative])
        for mode in checker.MODES:
            target = Path(self.tmp.name) / f"repaired-receipts-{mode[0]}"
            checker.prepare_checkout(self.root, target, mode)
            for relative, payload in payloads.items():
                with self.subTest(mode=mode[0], path=relative):
                    self.assertEqual((target / relative).read_bytes(), payload)

    def test_public_vkorc1_audit_input_hashes_survive_both_checkout_modes(self):
        paths = (
            Path("docs/tutorial/sources/08-04_vkorc1_warfarin_promoter_luciferase_gui.json"),
            Path("docs/examples/workflows/vkorc1_rs9923231_promoter_luciferase_assay_planning.json"),
        )
        payloads = {path: (checker.ROOT / path).read_bytes() for path in paths}
        attributes = (checker.ROOT / ".gitattributes").read_bytes()
        (self.root / ".gitattributes").write_bytes(attributes)
        for relative, payload in payloads.items():
            destination = self.root / relative
            destination.parent.mkdir(parents=True, exist_ok=True)
            destination.write_bytes(payload)
        checker.git(self.root, "add", "--all")
        checker.git(self.root, "-c", "commit.gpgsign=false", "commit", "--quiet",
                    "-m", "public audit input byte contracts")
        for index, relative in enumerate(paths):
            unprotected = b"\n".join(
                line for line in attributes.split(b"\n")
                if not line.startswith(relative.as_posix().encode() + b" ")
            )
            broken = Path(self.tmp.name) / f"public-audit-unprotected-{index}"
            prepare_unprotected_checkout(self.root, broken, unprotected)
            self.assertEqual((broken / relative).read_bytes(), payloads[relative].replace(b"\n", b"\r\n"))
            self.assertNotEqual(hashlib.sha256((broken / relative).read_bytes()).digest(),
                                hashlib.sha256(payloads[relative]).digest())
        for mode in checker.MODES:
            target = Path(self.tmp.name) / f"public-audit-{mode[0]}"
            checker.prepare_checkout(self.root, target, mode)
            for relative, payload in payloads.items():
                with self.subTest(mode=mode[0], path=relative):
                    self.assertEqual((target / relative).read_bytes(), payload)
                    self.assertEqual(hashlib.sha256((target / relative).read_bytes()).digest(),
                                     hashlib.sha256(payload).digest())

    def test_scroll_followup_receipt_hash_survives_both_checkout_modes(self):
        relative = Path("docs/audits/glen_followups_20261008/scroll_cpu_smoke_d1bf3fd3.json")
        payload = (checker.ROOT / relative).read_bytes()
        expected = "60add8b059510d660c36a25c1b7115138b3fbc8693dcecf2b42fd2316ddceab5"
        self.assertEqual(hashlib.sha256(payload).hexdigest(), expected)
        attributes = (checker.ROOT / ".gitattributes").read_bytes()
        (self.root / ".gitattributes").write_bytes(attributes)
        destination = self.root / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_bytes(payload)
        checker.git(self.root, "add", "--all")
        checker.git(self.root, "-c", "commit.gpgsign=false", "commit", "--quiet", "-m", "scroll smoke receipt")
        unprotected = b"\n".join(line for line in attributes.split(b"\n")
                                  if not line.startswith(relative.as_posix().encode() + b" "))
        broken = Path(self.tmp.name) / "scroll-receipt-unprotected"
        prepare_unprotected_checkout(self.root, broken, unprotected)
        self.assertNotEqual(hashlib.sha256((broken / relative).read_bytes()).hexdigest(), expected)
        for mode in checker.MODES:
            target = Path(self.tmp.name) / f"scroll-receipt-{mode[0]}"
            checker.prepare_checkout(self.root, target, mode)
            self.assertEqual((target / relative).read_bytes(), payload)

    def test_review_dependencies_keep_commit_dates_despite_newer_graphic_mtime(self):
        source = Path("docs/tutorial/sources/synthetic.json")
        graphic = Path("docs/screenshots/synthetic.png")
        for relative, payload in ((source, b"{}\n"), (graphic, b"synthetic image")):
            path = self.root / relative
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(payload)
        checker.git(self.root, "add", "--all")
        env = dict(os.environ, GIT_AUTHOR_DATE="2024-01-02T00:00:00Z",
                   GIT_COMMITTER_DATE="2024-01-02T00:00:00Z")
        subprocess.run(["git", "-C", str(self.root), "-c", "commit.gpgsign=false",
                        "-c", "core.hooksPath=.disabled-hooks", "commit", "--quiet",
                        "-m", "Synthetic review date tie"], env=env, check=True)
        expected = checker.git(self.root, "log", "-1", "--format=%ct", "--", str(source))
        for mode in checker.MODES:
            target = Path(self.tmp.name) / f"review-dates-{mode[0]}"
            checker.prepare_checkout(self.root, target, mode)
            os.utime(target / graphic, (4102444800, 4102444800))
            for relative in (source, graphic):
                with self.subTest(mode=mode[0], path=relative):
                    self.assertEqual(checker.git(target, "log", "-1", "--format=%ct",
                                                 "--", str(relative)), expected)
            self.assertNotEqual((target / graphic).stat().st_mtime,
                                (target / source).stat().st_mtime)

    def test_uniprot_review_evidence_hashes_survive_both_checkout_modes(self):
        evidence_root = Path("docs/tutorial/reproducibility")
        records = {
            "tp53_uniprot_projection_online/tp53_grch38_ensembl116.gb":
                "a6694752cd58bd16e6e5565725cc6d2fa3b127606f8ce3c6ee7890e327fecec3",
            "tp53_uniprot_projection_online/map-result.json":
                "bbdb7dad51f305061d24770f9e5bc5554ef34a1e330b8ee9556d7086a6028be5",
            "tp53_uniprot_projection_online/feature-coding-dna.json":
                "d12fa76661ec541f77fc4c3bfda15078edbc159a39d6d278e756ad67adc45f24",
            "tp73_uniprot_projection_audit_online/tp73_uniprot_projection.svg":
                "7cb5c32e7e0c9f605e55fca781cc3cce14028bbc0fd75ba19240ef53a95d8267",
        }
        attributes = (checker.ROOT / ".gitattributes").read_bytes()
        (self.root / ".gitattributes").write_bytes(attributes)
        expected = {}
        for name, digest in records.items():
            relative = evidence_root / name
            payload = (checker.ROOT / relative).read_bytes()
            provenance = (checker.ROOT / relative.parent / "README.md").read_text(
                encoding="utf-8")
            self.assertIn(digest, provenance)
            self.assertEqual(hashlib.sha256(payload).hexdigest(), digest,
                             f"Restore exact reviewed bytes / LF policy for {relative}")
            expected[relative] = (payload, digest)
            destination = self.root / relative
            destination.parent.mkdir(parents=True, exist_ok=True)
            destination.write_bytes(payload)
        checker.git(self.root, "add", "--all")
        checker.git(self.root, "-c", "commit.gpgsign=false", "commit", "--quiet",
                    "-m", "retained UniProt review evidence")

        protected_paths = {path.as_posix().encode() for path in expected}
        unprotected = b"\n".join(
            line for line in attributes.split(b"\n")
            if not any(line.startswith(path + b" ") for path in protected_paths)
        )
        broken = Path(self.tmp.name) / "uniprot-review-unprotected"
        prepare_unprotected_checkout(self.root, broken, unprotected)
        for relative, (payload, digest) in expected.items():
            with self.subTest(mode="unprotected-crlf", path=relative):
                converted = (broken / relative).read_bytes()
                self.assertEqual(converted, payload.replace(b"\n", b"\r\n"))
                self.assertNotEqual(hashlib.sha256(converted).hexdigest(), digest)
        for mode in checker.MODES:
            target = Path(self.tmp.name) / f"uniprot-review-{mode[0]}"
            checker.prepare_checkout(self.root, target, mode)
            for relative, (payload, digest) in expected.items():
                with self.subTest(mode=mode[0], path=relative):
                    retained = (target / relative).read_bytes()
                    self.assertEqual(retained, payload,
                                     f"{relative.as_posix()} text eol=lf is required")
                    self.assertEqual(hashlib.sha256(retained).hexdigest(), digest)

    def test_historical_panel_baselines_keep_raw_hashes_in_both_checkout_modes(self):
        generated = Path("docs/tutorial/generated")
        ledger = json.loads((checker.ROOT / generated / "report.json").read_bytes())
        attributes = (checker.ROOT / ".gitattributes").read_bytes()
        (self.root / ".gitattributes").write_bytes(attributes)
        expected = {}
        for name in ("patz1_endpoint_end_matrix", "patz1_routine_common_region_screen", "patz1_sybr_juc_panel"):
            artifact = f"artifacts/patz1_transcript_assay_panels_cli/artifacts/{name}.report.json"
            relative = generated / artifact
            payload = (checker.ROOT / relative).read_bytes()
            digest = ledger["file_checksums"][artifact]
            self.assertEqual(hashlib.sha256(payload).hexdigest(), digest)
            destination = self.root / relative
            destination.parent.mkdir(parents=True, exist_ok=True)
            destination.write_bytes(payload)
            expected[relative] = (payload, digest)
        checker.git(self.root, "add", "--all")
        checker.git(self.root, "-c", "commit.gpgsign=false", "commit", "--quiet",
                    "-m", "historical synthetic panel baselines")
        unprotected = b"\n".join(
            line for line in attributes.split(b"\n")
            if not line.startswith(b"docs/tutorial/generated/artifacts/patz1_transcript_assay_panels_cli/artifacts/")
        )
        broken = Path(self.tmp.name) / "panel-baselines-unprotected"
        prepare_unprotected_checkout(self.root, broken, unprotected)
        for relative, (payload, digest) in expected.items():
            with self.subTest(mode="unprotected-crlf", path=relative):
                converted = (broken / relative).read_bytes()
                self.assertEqual(converted, payload.replace(b"\n", b"\r\n"))
                self.assertNotEqual(hashlib.sha256(converted).hexdigest(), digest)
        for mode in checker.MODES:
            target = Path(self.tmp.name) / f"panel-baselines-{mode[0]}"
            checker.prepare_checkout(self.root, target, mode)
            for relative, (payload, digest) in expected.items():
                with self.subTest(mode=mode[0], path=relative):
                    retained = (target / relative).read_bytes()
                    self.assertEqual(retained, payload,
                                     f"Restore {relative.as_posix()} text eol=lf")
                    self.assertEqual(hashlib.sha256(retained).hexdigest(), digest)

    def test_gene_assay_gui_evidence_hashes_survive_both_checkout_modes(self):
        (self.root / ".gitattributes").write_bytes(
            (checker.ROOT / ".gitattributes").read_bytes())
        evidence_dir = Path("docs/screenshots/gene_assay_study_gui")
        source_dir = checker.ROOT / evidence_dir
        evidence = json.loads((source_dir / "evidence.json").read_bytes())
        expected = {}
        for row in evidence["captures"]:
            for field in ("raw_png", "semantic_snapshot", "context_svg"):
                relative = evidence_dir / row[field]
                payload = (checker.ROOT / relative).read_bytes()
                self.assertEqual(hashlib.sha256(payload).hexdigest(), row[f"{field}_sha256"])
                expected[relative] = (payload, row[f"{field}_sha256"])
                destination = self.root / relative
                destination.parent.mkdir(parents=True, exist_ok=True)
                destination.write_bytes(payload)
        self.assertEqual(len(expected), 18)
        checker.git(self.root, "add", "--all")
        checker.git(
            self.root,
            "-c",
            "commit.gpgsign=false",
            "commit",
            "--quiet",
            "-m",
            "public GUI evidence",
        )
        attributes = (checker.ROOT / ".gitattributes").read_bytes()
        unprotected = b"\n".join(
            line for line in attributes.split(b"\n")
            if b"docs/screenshots/gene_assay_study_gui/" not in line
        )
        broken = Path(self.tmp.name) / "gene-assay-evidence-unprotected"
        prepare_unprotected_checkout(self.root, broken, unprotected)
        changed = [
            relative
            for relative, (payload, _digest) in expected.items()
            if (broken / relative).read_bytes() != payload
        ]
        self.assertEqual(
            sorted(changed),
            sorted(relative for relative in expected if relative.suffix in (".json", ".svg")),
        )
        for mode in checker.MODES:
            target = Path(self.tmp.name) / f"gene-assay-evidence-{mode[0]}"
            checker.prepare_checkout(self.root, target, mode)
            for relative, (payload, digest) in expected.items():
                with self.subTest(mode=mode[0], path=relative):
                    checked_out = (target / relative).read_bytes()
                    self.assertEqual(checked_out, payload)
                    self.assertEqual(hashlib.sha256(checked_out).hexdigest(), digest)

    def test_tss_regulatory_gui_evidence_hashes_survive_both_checkout_modes(self):
        attributes = (checker.ROOT / ".gitattributes").read_bytes()
        (self.root / ".gitattributes").write_bytes(attributes)
        evidence = json.loads((checker.ROOT /
            "docs/screenshots/tss_regulatory_view_gui/evidence.json").read_bytes())
        records = list(evidence["inputs"].values())
        records.append(evidence["agent_assistant"]["checkpoint"])
        for capture in evidence["captures"]:
            records.extend(capture[field] for field in ("raw_png", "semantic_snapshot"))
        expected = {}
        for record in records:
            relative = Path(record["path"])
            payload = (checker.ROOT / relative).read_bytes()
            self.assertEqual(hashlib.sha256(payload).hexdigest(), record["sha256"])
            expected[relative] = (payload, record["sha256"])
            destination = self.root / relative
            destination.parent.mkdir(parents=True, exist_ok=True)
            destination.write_bytes(payload)
        self.assertEqual(len(expected), 8)
        checker.git(self.root, "add", "--all")
        checker.git(self.root, "-c", "commit.gpgsign=false", "commit", "--quiet",
                    "-m", "retained synthetic TSS GUI evidence")

        unprotected = b"\n".join(
            line for line in attributes.split(b"\n")
            if b"docs/screenshots/tss_regulatory_view_gui/" not in line
        )
        broken = Path(self.tmp.name) / "tss-evidence-unprotected"
        prepare_unprotected_checkout(self.root, broken, unprotected)
        self.assertEqual(
            sorted(relative for relative, (payload, _) in expected.items()
                   if (broken / relative).read_bytes() != payload),
            sorted(Path(capture["semantic_snapshot"]["path"])
                   for capture in evidence["captures"]),
        )
        for mode in checker.MODES:
            target = Path(self.tmp.name) / f"tss-evidence-{mode[0]}"
            checker.prepare_checkout(self.root, target, mode)
            for relative, (payload, digest) in expected.items():
                with self.subTest(mode=mode[0], path=relative):
                    checked_out = (target / relative).read_bytes()
                    self.assertEqual(checked_out, payload,
                                     f"{relative} needs a scoped .gitattributes LF rule")
                    self.assertEqual(hashlib.sha256(checked_out).hexdigest(), digest)

    def test_tss_collection_gui_evidence_hashes_survive_both_checkout_modes(self):
        attributes = (checker.ROOT / ".gitattributes").read_bytes()
        (self.root / ".gitattributes").write_bytes(attributes)
        evidence = json.loads((checker.ROOT /
            "docs/screenshots/tss_collection_gui/evidence.json").read_bytes())
        records = []
        for capture in evidence["captures"]:
            records.extend(capture[field] for field in
                           ("raw_png", "semantic_snapshot", "screenshot_receipt"))

        expected = {}
        for record in records:
            relative = Path(record["path"])
            payload = (checker.ROOT / relative).read_bytes()
            self.assertEqual(hashlib.sha256(payload).hexdigest(), record["sha256"])
            expected[relative] = (payload, record["sha256"])
            destination = self.root / relative
            destination.parent.mkdir(parents=True, exist_ok=True)
            destination.write_bytes(payload)
        self.assertEqual(len(expected), 6)
        checker.git(self.root, "add", "--all")
        checker.git(self.root, "-c", "commit.gpgsign=false", "commit", "--quiet",
                    "-m", "retained TSS collection GUI evidence")

        unprotected = b"\n".join(
            line for line in attributes.split(b"\n")
            if b"docs/screenshots/tss_collection_gui/" not in line
        )
        broken = Path(self.tmp.name) / "tss-collection-evidence-unprotected"
        prepare_unprotected_checkout(self.root, broken, unprotected)
        text_paths = sorted(relative for relative in expected
                            if relative.suffix == ".json")
        self.assertEqual(len(text_paths), 4)
        self.assertEqual(
            sorted(relative for relative, (payload, _) in expected.items()
                   if (broken / relative).read_bytes() != payload),
            text_paths,
        )
        for mode in checker.MODES:
            target = Path(self.tmp.name) / f"tss-collection-evidence-{mode[0]}"
            checker.prepare_checkout(self.root, target, mode)
            for relative, (payload, digest) in expected.items():
                with self.subTest(mode=mode[0], path=relative):
                    checked_out = (target / relative).read_bytes()
                    self.assertEqual(checked_out, payload,
                                     f"{relative} needs a scoped .gitattributes LF rule")
                    self.assertEqual(hashlib.sha256(checked_out).hexdigest(), digest)

    def test_tp73_tss_and_pcr_evidence_hashes_survive_both_checkout_modes(self):
        attributes = (checker.ROOT / ".gitattributes").read_bytes()
        (self.root / ".gitattributes").write_bytes(attributes)
        evidence = json.loads((checker.ROOT /
            "docs/screenshots/tp73_dnp73_factor_curves/evidence.json").read_bytes())
        records = [evidence["input"]]
        for capture in evidence["captures"]:
            records.extend(capture[field] for field in ("raw_png", "semantic_snapshot"))

        pcr_dir = Path("docs/screenshots/tutorial_gui_acceptance/simple_pcr_selection_gui")
        pcr = json.loads((checker.ROOT / pcr_dir /
            "inspect_primer_report.screenshot.json").read_bytes())
        # The capture receipt records archived Linux paths; the retained files
        # are siblings of that receipt, with the same basenames and byte hashes.
        for record in [pcr["capture"]["raw"], *pcr["derived_views"]]:
            records.append({"path": pcr_dir / Path(record["path"]).name,
                            "sha256": record["sha256"]})
        snapshot = pcr["semantic_snapshot"]
        records.append({"path": pcr_dir / Path(snapshot["retained_path"]).name,
                        "sha256": snapshot["retained_file_sha256"]})

        expected = {}
        for record in records:
            relative = Path(record["path"])
            payload = (checker.ROOT / relative).read_bytes()
            self.assertEqual(hashlib.sha256(payload).hexdigest(), record["sha256"])
            expected[relative] = (payload, record["sha256"])
            destination = self.root / relative
            destination.parent.mkdir(parents=True, exist_ok=True)
            destination.write_bytes(payload)
        self.assertEqual(len(expected), 11)
        checker.git(self.root, "add", "--all")
        checker.git(self.root, "-c", "commit.gpgsign=false", "commit", "--quiet",
                    "-m", "retained TP73 TSS and Simple PCR evidence")

        scopes = (b"docs/screenshots/tp73_dnp73_factor_curves/",
                  b"docs/screenshots/tutorial_gui_acceptance/",
                  b"test_files/tp73.ncbi.gb")
        unprotected = b"\n".join(
            line for line in attributes.split(b"\n")
            if not any(scope in line for scope in scopes)
        )
        broken = Path(self.tmp.name) / "tp73-tss-pcr-unprotected"
        prepare_unprotected_checkout(self.root, broken, unprotected)
        text_paths = sorted(relative for relative in expected
                            if relative.suffix in (".json", ".svg", ".gb"))
        self.assertEqual(len(text_paths), 7)
        self.assertEqual(
            sorted(relative for relative, (payload, _) in expected.items()
                   if (broken / relative).read_bytes() != payload),
            text_paths,
        )
        for mode in checker.MODES:
            target = Path(self.tmp.name) / f"tp73-tss-pcr-{mode[0]}"
            checker.prepare_checkout(self.root, target, mode)
            for relative, (payload, digest) in expected.items():
                with self.subTest(mode=mode[0], path=relative):
                    checked_out = (target / relative).read_bytes()
                    self.assertEqual(checked_out, payload,
                                     f"{relative} needs a scoped .gitattributes LF rule")
                    self.assertEqual(hashlib.sha256(checked_out).hexdigest(), digest)

    def test_generated_parity_matrix_stays_byte_exact_in_both_checkout_modes(self):
        relative = "docs/gui_cli_mcp_parity.md"
        payload = (checker.ROOT / relative).read_bytes()
        self.assertIn(b"\n", payload)
        self.assertFalse(b"\r" in payload,
                         f"{relative} requires a scoped .gitattributes text eol=lf rule")
        attributes = (checker.ROOT / ".gitattributes").read_bytes()
        (self.root / ".gitattributes").write_bytes(attributes)
        destination = self.root / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_bytes(payload)
        checker.git(self.root, "add", "--all")
        checker.git(self.root, "-c", "commit.gpgsign=false", "commit", "--quiet",
                    "-m", "generated capability matrix")

        # Negative control: the missing rule must reproduce CRLF checkout drift.
        unprotected = b"\n".join(
            line for line in attributes.split(b"\n")
            if relative.encode() not in line
        )
        broken = Path(self.tmp.name) / "parity-unprotected"
        prepare_unprotected_checkout(self.root, broken, unprotected)
        self.assertTrue(
            (broken / relative).read_bytes() == payload.replace(b"\n", b"\r\n"),
            "The missing-rule control must reproduce CRLF conversion of the parity matrix",
        )
        for mode in checker.MODES:
            with self.subTest(mode=mode[0]):
                target = Path(self.tmp.name) / f"parity-{mode[0]}"
                checker.prepare_checkout(self.root, target, mode)
                self.assertTrue(
                    (target / relative).read_bytes() == payload,
                    f"{relative} changed during {mode[0]} checkout; preserve the "
                    "generator's exact bytes with .gitattributes text eol=lf",
                )

    def test_replay_uses_existing_binary_and_forces_offline(self):
        with patch.dict(os.environ, {"GENTLE_TEST_ONLINE": "1"}), \
                patch.object(checker, "check_tutorial_visuals") as visuals, \
                patch.object(checker.subprocess, "run") as run:
            checker.check_checkout(Path("existing-binary"), self.root, "crlf", 12)
        visuals.assert_called_once_with(self.root)
        self.assertEqual([call.args[0] for call in run.call_args_list],
                         [["existing-binary", "parity-matrix-check"],
                          ["existing-binary", "tutorial-catalog-check"],
                          ["existing-binary", "--check"], ["existing-binary", "tutorial-check"]])
        for call in run.call_args_list:
            self.assertNotIn("GENTLE_TEST_ONLINE", call.kwargs["env"])
            self.assertEqual(call.kwargs["cwd"], self.root)
            self.assertTrue(call.kwargs["check"])
            self.assertEqual(call.kwargs["timeout"], 12)

    def test_bound_tss_and_probe_fixtures_survive_crlf_checkout(self):
        (self.root / ".gitattributes").write_bytes(
            (checker.ROOT / ".gitattributes").read_bytes())
        fixture_dirs = (
            "test_files/fixtures/tss_profiles",
            "test_files/fixtures/probe_region_outputs/clariom_e_mtab_14704_tp73_validation",
        )
        expected = {}
        for directory in fixture_dirs:
            for source in (checker.ROOT / directory).rglob("*"):
                if not source.is_file():
                    continue
                relative = source.relative_to(checker.ROOT)
                expected[relative] = source.read_bytes()
                destination = self.root / relative
                destination.parent.mkdir(parents=True, exist_ok=True)
                destination.write_bytes(expected[relative])
        provenance_path = Path(fixture_dirs[1]) / "provenance.json"
        provenance = json.loads(expected[provenance_path])
        inputs = provenance["input_fingerprints"]
        self.assertTrue(inputs, "expected retained adapter input provenance")
        for binding in inputs:
            relative = Path(binding["path"])
            self.assertFalse(relative.is_absolute())
            self.assertNotIn("..", relative.parts)
            expected[relative] = (checker.ROOT / relative).read_bytes()
            destination = self.root / relative
            destination.parent.mkdir(parents=True, exist_ok=True)
            destination.write_bytes(expected[relative])
        self.assertTrue(expected)
        checker.git(self.root, "add", "--all")
        checker.git(self.root, "-c", "commit.gpgsign=false", "commit", "--quiet", "-m", "bound fixtures")
        for mode in checker.MODES:
            target = Path(self.tmp.name) / f"bound-{mode[0]}"
            checker.prepare_checkout(self.root, target, mode)
            for relative, payload in expected.items():
                with self.subTest(mode=mode[0], path=relative):
                    self.assertEqual((target / relative).read_bytes(), payload)
            for binding in inputs:
                with self.subTest(mode=mode[0], bound_input=binding["path"]):
                    self.assertEqual(
                        "sha256:" + hashlib.sha256(
                            (target / binding["path"]).read_bytes()).hexdigest(),
                        binding["sha256"],
                        "Preserve the exact input bytes with a scoped .gitattributes LF rule; "
                        "do not change provenance hashes or normalize input bytes")
            bundle = target / fixture_dirs[0]
            for line in (bundle / "SHA256SUMS").read_text().splitlines():
                digest, relative = line.split(maxsplit=1)
                self.assertEqual(hashlib.sha256((bundle / relative).read_bytes()).hexdigest(), digest)

    def test_failed_check_or_timeout_cannot_be_reported_as_success(self):
        for error in (subprocess.CalledProcessError(1, "validator"),
                      subprocess.TimeoutExpired("validator", 12)):
            for failed_call in (0, 1, 2, 3):
                with self.subTest(error=type(error).__name__, failed_call=failed_call), \
                        patch.object(checker, "check_tutorial_visuals"), \
                        patch.object(checker.subprocess, "run",
                                     side_effect=[None] * failed_call + [error]) as run:
                    with self.assertRaisesRegex(RuntimeError, "crlf checkout failed"):
                        checker.check_checkout(Path("existing-binary"), self.root, "crlf", 12)
                    self.assertEqual(run.call_count, failed_call + 1)

    def test_main_checks_both_modes_at_one_revision_without_building(self):
        binary = Path(self.tmp.name) / "existing-binary"
        binary.touch()
        arguments = ["--repo-root", str(self.root), "--binary", str(binary)]
        with patch.object(checker, "prepare_checkout") as prepare, \
                patch.object(checker, "check_checkout") as check:
            self.assertEqual(checker.main(arguments), 0)
        self.assertEqual([call.args[2][0] for call in prepare.call_args_list], ["lf", "crlf"])
        self.assertEqual(prepare.call_args_list[0].args[4], prepare.call_args_list[1].args[4])
        self.assertEqual([call.args[2] for call in check.call_args_list], ["lf", "crlf"])
        with patch.object(checker, "prepare_checkout"), \
                patch.object(checker, "check_checkout", side_effect=[None, RuntimeError("crlf failed")]):
            with self.assertRaisesRegex(RuntimeError, "crlf failed"):
                checker.main(arguments)


class TutorialVisualContractTests(unittest.TestCase):
    """Literal synthetic docs/images, created only in temporary directories.

    Recreation is the setup below: two unrelated generated chapters, one manual
    page and an image-free reference. No biological identifiers or GUI captures.
    Real repository assets are used only by the named archive/checkout checks.
    """

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix="gentle-visual-test-")
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name) / "source"
        self.root.mkdir()
        self.svg = b'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 20 10"><text x="1" y="5">Synthetic</text></svg>\n'
        self.write("docs/figures/synthetic.svg", self.svg)
        self.write("docs/screenshots/legacy.png", self.png())
        self.sources = []
        for identifier in ("unrelated_alpha", "unrelated_beta"):
            self.sources.append({"schema": "gentle.tutorial_source.v4", "id": identifier,
                                 "graphics": [{"kind": "generated", "path": "docs/figures/synthetic.svg",
                                               "caption": "Synthetic figure", "illustrates_step": 2}],
                                 "catalog": {"path": f"docs/tutorial/generated/chapters/{identifier}.md",
                                             "type": "executable_chapter"},
                                 "generated_chapter": {"gui_steps": ["First", "Second"]}})
            self.generated_page(identifier)
        self.sources.append({"schema": "gentle.tutorial_source.v4", "id": "manual",
                             "catalog": {"path": "docs/tutorial/manual.md", "type": "gui_walkthrough"},
                             "graphics": [{"kind": "screenshot", "path": "docs/screenshots/legacy.png",
                                           "caption": "Historical capture, not current acceptance",
                                           "illustrates_step": 1, "capture_date": "2026-01-02"}]})
        self.write("docs/tutorial/manual.md", b'# Manual\n\n![Historical capture](../screenshots/legacy.png)\n')
        self.sources.append({"schema": "gentle.tutorial_source.v4", "id": "reference",
                             "catalog": {"path": "docs/tutorial/reference.md", "type": "operational_reference"}})
        self.write("docs/tutorial/reference.md", b'# Reference\nNo decorative images needed.\n')
        self.sync()

    @staticmethod
    def png(width=1, height=1):
        def chunk(kind, data):
            return len(data).to_bytes(4, "big") + kind + data + zlib.crc32(kind + data).to_bytes(4, "big")
        header = width.to_bytes(4, "big") + height.to_bytes(4, "big") + bytes([8, 6, 0, 0, 0])
        return (b"\x89PNG\r\n\x1a\n" + chunk(b"IHDR", header)
                + chunk(b"IDAT", zlib.compress(b"\x00\xff\x00\x00\xff")) + chunk(b"IEND", b""))

    def write(self, relative, payload):
        path = self.root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(payload)
        return path

    def sync(self):
        entries, chapters = [], []
        for source in self.sources:
            self.write(f"docs/tutorial/sources/{source['id']}.json", json.dumps(source).encode() + b"\n")
            entries.append({"id": source["id"], **source["catalog"], "graphics": source.get("graphics", [])})
            if "generated_chapter" in source:
                chapters.append({"id": source["id"], **source["generated_chapter"], "graphics": source.get("graphics", [])})
        self.write("docs/tutorial/catalog.json", json.dumps({"schema": "gentle.tutorial_catalog.v2", "entries": entries}).encode() + b"\n")
        self.write("docs/tutorial/manifest.json", json.dumps({"schema": "gentle.tutorial_manifest.v2", "chapters": chapters}).encode() + b"\n")

    def generated_page(self, identifier, image_step=2, duplicate=False):
        image = '![Synthetic figure](../../../figures/synthetic.svg)\n'
        text = (f'---\nchapter_id: "{identifier}"\n---\n\n## Walkthrough: GUI, CLI and Inner Agent\n'
                + '### Step 1: First\n' + (image if image_step == 1 else '')
                + '### Step 2: Second\n' + (image if image_step == 2 or duplicate else '')
                + '## Appendix\n')
        self.write(f"docs/tutorial/generated/chapters/{identifier}.md", text.encode())

    def screenshot(self):
        raw = self.write("docs/screenshots/synthetic/capture.raw.png", self.png())
        snapshot = self.write("docs/screenshots/synthetic/capture.snapshot.json", b'{"synthetic": true}\n')
        record = {"schema": "gentle.tutorial_gui_screenshot_evidence.v1", "source_revision": "12" * 20,
                  "gentle_binary_sha256": "ab" * 32, "tutorial_manifest_sha256": "cd" * 32,
                  "acceptance_contract_sha256": "ef" * 32, "step_sha256": "34" * 32,
                  "chapter_id": "manual", "step_id": "capture", "prose_step": 1,
                  "capture": {"captured_at_unix_ms": 1767312000000,
                              "raw": {"path": "/original-private-run/capture.raw.png",
                                      "sha256": hashlib.sha256(raw.read_bytes()).hexdigest()}},
                  "derived_views": [],
                  "semantic_snapshot": {"retained_path": "/original-private-run/capture.snapshot.json",
                                        "retained_file_sha256": hashlib.sha256(snapshot.read_bytes()).hexdigest()}}
        self.write("docs/screenshots/synthetic/capture.screenshot.json", json.dumps(record).encode() + b"\n")
        self.sources[2]["graphics"][0]["path"] = "docs/screenshots/synthetic/capture.raw.png"
        self.write("docs/tutorial/manual.md", b'![Historical capture](../screenshots/synthetic/capture.raw.png)\n')
        self.sync()
        return record

    def test_synthetic_coverage_is_generic_and_reference_needs_no_image(self):
        report = checker.check_tutorial_visuals(self.root)
        self.assertEqual(report["declared_visual_count"], 3)
        rows = {row["id"]: row for row in report["tutorials"]}
        for identifier in ("unrelated_alpha", "unrelated_beta"):
            self.assertEqual(rows[identifier]["visuals"][0]["step"], 2)
        self.assertEqual(rows["reference"]["visuals"], [])
        self.assertEqual(rows["manual"]["visuals"][0]["status"], "legacy_date_only_not_acceptance")

    def test_real_catalog_and_nondefault_generated_chapters_are_checked(self):
        report = checker.check_tutorial_visuals(checker.ROOT)
        rows = {row["id"]: row for row in report["tutorials"]}
        self.assertTrue(rows["tp53_isoform_architecture_online"]["visuals"])
        self.assertTrue(rows["patz1_gene_locus_evidence_offline"]["visuals"])
        self.assertEqual(rows["measured_gel_sizing_gui_cli"]["visuals"][0]["status"],
                         "deterministic_figure_not_capture")
        self.assertIn("not current-candidate", report["non_claim"])

    def test_glen_visual_forward_port_keeps_semantic_placements_and_exclusions(self):
        # Existing assets proposed by d70bf25d, not new biological test fixtures.
        report = checker.check_tutorial_visuals(checker.ROOT)
        rows = {row["id"]: row for row in report["tutorials"]}
        expected = {
            "agent_interfaces": ["docs/figures/tp73_outer_agent_evidence_viewer.svg"],
            "gibson_specialist_testing_gui": [
                "docs/figures/gibson_single_insert_protocol_cartoon.svg",
                "docs/screenshots/tutorial_gibson_specialist_testing_08_lineage_graph.png"],
            "gibson_arrangements_gui": [
                "docs/figures/gibson_single_insert_arrangement_gel.svg",
                "docs/figures/gibson_single_insert_lineage.svg"],
            "gibson_physical_rack_gui": [
                "docs/figures/gibson_single_insert_rack_isometric.svg",
                "docs/figures/gibson_single_insert_rack_topdown.svg"],
            "qpcr_exon_junctions_gui": [
                "docs/figures/qpcr_assay_protocol_cartoon.svg",
                "docs/figures/tp73_as3_shared_cdna_qpcr.svg",
                "docs/figures/tp73_as3_nr187362_cdna_qpcr.svg"],
            "gibson_specialist_testing_baseline": [
                "docs/screenshots/screenshot_GUI_help_tutorial_testing_gibson.png"],
            "prepare_reference_genome_online": [
                "docs/screenshots/screenshot_GUI_prepare_reference_genome.png"],
        }
        for identifier, paths in expected.items():
            with self.subTest(identifier=identifier):
                visuals = {row["path"]: row for row in rows[identifier]["visuals"]}
                self.assertTrue(set(paths) <= visuals.keys())
        self.assertEqual(rows["tp73_cdna_genomic_dotplot_online"]["visuals"], [])
        self.assertNotIn("docs/screenshots/tp73_inner_agent_evidence_viewer.png",
                         [row["path"] for row in rows["agent_interfaces"]["visuals"]])
        self.assertNotIn("docs/figures/pgex_blunt_sequence_confirmation_hero.svg",
                         [row["path"] for row in rows["sequencing_confirmation_gui"]["visuals"]])
        for identifier in ("gibson_specialist_testing_baseline", "prepare_reference_genome_online"):
            self.assertEqual(rows[identifier]["visuals"][0]["status"],
                             "legacy_date_only_not_acceptance")
        preparation = rows["prepare_reference_genome_online"]["visuals"][0]
        self.assertEqual(preparation["step"], 1)
        dotplot = (checker.ROOT / "docs/tutorial/02-03_tp73_cdna_genomic_dotplot_gui.md")
        original_dotplot = ["tutorial_cdna_genomic_02_fetch_cdna_dialog.png",
                            "tutorial_cdna_genomic_03_prepare_genome_dialog.png",
                            "tutorial_cdna_genomic_04_retrieve_tp73_genomic_dialog.png",
                            "tutorial_cdna_genomic_07_dotplot_window.png",
                            "tutorial_cdna_genomic_09_dotplot_as_node.png"]
        placements = [Path(a or b).name for a, b in checker.markdown_visuals(dotplot.read_text())]
        for name in original_dotplot:
            self.assertEqual(placements.count(name), 1)
        review = (checker.ROOT / "docs/tutorial/visual_forward_port_review.md").read_text()
        reviewed = re.findall(r"^\| \d{2}\.\d{2} \| `([^`]+)` \|", review, re.M)
        expected_review = {path.removeprefix("docs/") for paths in expected.values() for path in paths}
        expected_review.update("screenshots/" + name for name in original_dotplot)
        expected_review.update(["screenshots/tp73_inner_agent_evidence_viewer.png",
                                "figures/pgex_blunt_sequence_confirmation_hero.svg"])
        self.assertEqual(len(reviewed), 19)
        self.assertEqual(set(reviewed), expected_review)

    def test_missing_asset_and_unsafe_paths_fail(self):
        for target in ("docs/figures/missing.svg", "../outside.svg", "/tmp/outside.svg",
                       "C:\\outside.png", "docs/figures/../../outside.svg", "test_files/image.png",
                       "docs/tutorial/inputs/unapproved/image.png"):
            with self.subTest(target=target), self.assertRaisesRegex(RuntimeError, "tutorial visual"):
                checker.visual_path(self.root, target)

    def test_symlink_asset_escape_is_refused(self):
        outside = Path(self.tmp.name) / "outside.svg"
        outside.write_bytes(self.svg)
        link = self.root / "docs/figures/linked.svg"
        try:
            link.symlink_to(outside)
        except OSError:
            self.skipTest("Host cannot create a symlink")
        with self.assertRaisesRegex(RuntimeError, "unsafe_path"):
            checker.visual_path(self.root, "docs/figures/linked.svg")

    def test_png_signature_dimensions_and_header_integrity(self):
        for payload in (b"not a PNG", self.png()[:24], self.png(width=0),
                        self.png()[:20] + b"\0\0\0\2" + self.png()[24:]):
            path = self.write("docs/screenshots/broken.png", payload)
            with self.subTest(payload=payload[:24]), self.assertRaisesRegex(RuntimeError, "invalid_png"):
                checker.inspect_visual_asset(self.root, path)

    def test_svg_structure_dimensions_entities_and_budget(self):
        invalid = (b"<svg>", b'<svg xmlns="http://www.w3.org/2000/svg"/>',
                   b'<svg viewBox="0 0 NaN 10"/>', b'<svg viewBox="0 0 -1 10"/>',
                   b'<svg width="100%" height="100%"/>', b'<other viewBox="0 0 20 10"/>',
                   b'<!DOCTYPE svg [<!ENTITY x "test">]><svg viewBox="0 0 20 10"/>')
        for payload in invalid:
            path = self.write("docs/figures/broken.svg", payload)
            with self.subTest(payload=payload), self.assertRaisesRegex(RuntimeError, "tutorial visual"):
                checker.inspect_visual_asset(self.root, path)
        path = self.write("docs/figures/explicit.svg", b'<svg width="20px" height="10pt"/>')
        self.assertEqual(checker.inspect_visual_asset(self.root, path)["format"], "svg")
        with patch.object(checker, "MAX_SVG_BYTES", 1), self.assertRaisesRegex(RuntimeError, "svg_size_budget"):
            checker.inspect_visual_asset(self.root, path)

    def test_remote_svg_image_is_refused_without_network(self):
        path = self.write("docs/figures/remote.svg", b'<svg viewBox="0 0 20 10"><image href="https://example.invalid/a.png"/></svg>')
        with self.assertRaisesRegex(RuntimeError, "unsafe_path"):
            checker.inspect_visual_asset(self.root, path)

    def test_missing_and_malformed_metadata_has_specific_errors(self):
        original = dict(self.sources[0]["graphics"][0])
        for field, value, code in (("caption", "", "missing_caption"), ("path", None, "unsafe_path"),
                                   ("illustrates_step", 0, "invalid_step"), ("illustrates_step", True, "invalid_step"),
                                   ("kind", "unknown", "invalid_kind")):
            self.sources[0]["graphics"][0] = {**original, field: value}
            self.sync()
            with self.subTest(field=field, value=value), self.assertRaisesRegex(RuntimeError, code):
                checker.check_tutorial_visuals(self.root)

    def test_generated_placement_is_step_bound_and_unique(self):
        self.generated_page("unrelated_beta", image_step=1)
        with self.assertRaisesRegex(RuntimeError, "wrong_step_placement"):
            checker.check_tutorial_visuals(self.root)
        self.generated_page("unrelated_beta", image_step=1, duplicate=True)
        with self.assertRaisesRegex(RuntimeError, "duplicate_step_placement"):
            checker.check_tutorial_visuals(self.root)

    def test_wrong_metadata_container_type_is_actionable(self):
        self.sources[0]["graphics"] = None
        self.sync()
        with self.assertRaisesRegex(RuntimeError, "invalid_metadata"):
            checker.check_tutorial_visuals(self.root)

    def test_manifest_cannot_drop_graphics_metadata(self):
        path = self.root / "docs/tutorial/manifest.json"
        manifest = json.loads(path.read_bytes())
        manifest["chapters"][1]["graphics"] = []
        path.write_bytes(json.dumps(manifest).encode())
        with self.assertRaisesRegex(RuntimeError, "metadata_projection_mismatch"):
            checker.check_tutorial_visuals(self.root)

    def test_figure_cannot_be_called_a_gui_screenshot(self):
        self.sources[0]["graphics"][0].update(kind="screenshot", capture_date="2026-01-02")
        self.sync()
        with self.assertRaisesRegex(RuntimeError, "figure_as_screenshot"):
            checker.check_tutorial_visuals(self.root)

    def test_generated_figure_cannot_borrow_capture_metadata(self):
        self.sources[0]["graphics"][0]["capture_date"] = "2026-01-02"
        self.sync()
        with self.assertRaisesRegex(RuntimeError, "screenshot_as_figure"):
            checker.check_tutorial_visuals(self.root)

    def test_capture_date_is_required_and_valid(self):
        for date in (None, "", "2026-02-30"):
            self.sources[2]["graphics"][0]["capture_date"] = date
            self.sync()
            with self.subTest(date=date), self.assertRaisesRegex(RuntimeError, "missing_capture_date"):
                checker.check_tutorial_visuals(self.root)

    def test_retained_capture_identity_and_hashes_are_not_current_acceptance(self):
        record = self.screenshot()
        report = checker.check_tutorial_visuals(self.root)
        visual = next(row for row in report["tutorials"] if row["id"] == "manual")["visuals"][0]
        self.assertEqual(visual["status"], "historical_provenance_bound")
        self.assertEqual(visual["revision"], record["source_revision"])
        self.write("docs/screenshots/synthetic/capture.raw.png", self.png(height=2))
        with self.assertRaisesRegex(RuntimeError, "capture_hash_mismatch"):
            checker.check_tutorial_visuals(self.root)

    def test_required_capture_provenance_and_snapshot_cannot_disappear(self):
        record = self.screenshot()
        sidecar = "docs/screenshots/synthetic/capture.screenshot.json"
        for field in ("source_revision", "gentle_binary_sha256", "step_sha256", "prose_step"):
            changed = dict(record)
            changed.pop(field)
            self.write(sidecar, json.dumps(changed).encode())
            with self.subTest(field=field), self.assertRaisesRegex(RuntimeError, "missing_capture"):
                checker.check_tutorial_visuals(self.root)
        self.write(sidecar, json.dumps(record).encode())
        (self.root / "docs/screenshots/synthetic/capture.snapshot.json").unlink()
        with self.assertRaisesRegex(RuntimeError, "missing_capture_snapshot"):
            checker.check_tutorial_visuals(self.root)
        (self.root / sidecar).unlink()
        with self.assertRaisesRegex(RuntimeError, "missing_capture_provenance"):
            checker.check_tutorial_visuals(self.root)

    def test_malformed_capture_field_types_are_actionable(self):
        record = self.screenshot()
        sidecar = "docs/screenshots/synthetic/capture.screenshot.json"
        for field, value in (("source_revision", None), ("chapter_id", 1), ("capture", None),
                             ("derived_views", None), ("derived_views", [None]),
                             ("semantic_snapshot", None)):
            self.write(sidecar, json.dumps({**record, field: value}).encode())
            with self.subTest(field=field), self.assertRaisesRegex(RuntimeError, "missing_capture"):
                checker.check_tutorial_visuals(self.root)
        capture = {**record["capture"], "captured_at_unix_ms": 10 ** 100}
        self.write(sidecar, json.dumps({**record, "capture": capture}).encode())
        with self.assertRaisesRegex(RuntimeError, "missing_capture_date"):
            checker.check_tutorial_visuals(self.root)

    def test_manual_archive_cannot_bind_an_unrelated_same_named_asset(self):
        path = self.root / "docs/screenshots/legacy.png"
        record = {"schema": "synthetic.manual_archive.v1",
                  "producer": {"revision": "12" * 20, "gentle_gui_sha256": "ab" * 32},
                  "raw_png": "docs/screenshots/unrelated/legacy.png",
                  "raw_png_sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
        self.write("docs/screenshots/evidence.json", json.dumps(record).encode())
        with self.assertRaisesRegex(RuntimeError, "capture_hash_mismatch"):
            checker.check_tutorial_visuals(self.root)
        record["raw_png"] = "legacy.png"
        self.write("docs/screenshots/evidence.json", json.dumps(record).encode())
        report = checker.check_tutorial_visuals(self.root)
        visual = next(row for row in report["tutorials"] if row["id"] == "manual")["visuals"][0]
        self.assertEqual(visual["status"], "historical_manual_bound")

    def test_capture_snapshot_symlink_cannot_escape_asset_roots(self):
        self.screenshot()
        snapshot = self.root / "docs/screenshots/synthetic/capture.snapshot.json"
        outside = Path(self.tmp.name) / "outside.snapshot.json"
        outside.write_bytes(snapshot.read_bytes())
        snapshot.unlink()
        try:
            snapshot.symlink_to(outside)
        except OSError:
            self.skipTest("Host cannot create a symlink")
        with self.assertRaisesRegex(RuntimeError, "unsafe_path"):
            checker.check_tutorial_visuals(self.root)

    def test_manual_markdown_links_and_hash_bound_crops_survive_both_checkouts(self):
        # Real current assets, source declarations and captions supplement the
        # unrelated synthetic families; no screenshot is manufactured here.
        for identifier, source_name in (("gibson_arrangements_gui", "03-07_gibson_arrangements_gui"),
                                        ("load_and_digest_pgex", "03-01_load_and_digest_pgex"),
                                        ("protein_expression_handoff_cli", "09-03_protein_expression_handoff_cli")):
            source = json.loads((checker.ROOT / f"docs/tutorial/sources/{source_name}.json").read_bytes())
            self.sources.append(source)
            guide = source["catalog"]["path"]
            self.write(guide, (checker.ROOT / guide).read_bytes())
            for graphic in source["graphics"]:
                relative = Path(graphic["path"])
                self.write(relative, (checker.ROOT / relative).read_bytes())
                if relative.parent.as_posix().startswith("docs/screenshots/"):
                    for asset in (checker.ROOT / relative.parent).glob("*"):
                        if asset.is_file():
                            self.write(relative.parent / asset.name, asset.read_bytes())
        self.sync()
        self.write(".gitattributes", (checker.ROOT / ".gitattributes").read_bytes())
        checker.git(self.root, "init", "--quiet")
        checker.git(self.root, "config", "core.autocrlf", "false")
        checker.git(self.root, "add", "--all")
        checker.git(self.root, "-c", "user.name=Synthetic test", "-c", "user.email=fixture@example.invalid",
                    "-c", "commit.gpgsign=false", "commit", "--quiet", "-m", "Synthetic visual checkout plus retained assets")
        expected = checker.check_tutorial_visuals(self.root)
        for mode in checker.MODES:
            with self.subTest(mode=mode[0]):
                target = Path(self.tmp.name) / mode[0]
                checker.prepare_checkout(self.root, target, mode)
                self.assertEqual(checker.check_tutorial_visuals(target), expected)
                expression = target / "docs/screenshots/tutorial_gui_acceptance/protein_expression_handoff_cli"
                receipt = json.loads((expression / "evidence.json").read_bytes())
                self.assertEqual(checker.visual_sha256(expression / receipt["project"]["path"]),
                                 receipt["project"]["before_sha256"])
                self.assertNotIn(b"\r", (target / "docs/tutorial/manual.md").read_bytes())
        attributes = (self.root / ".gitattributes").read_bytes()
        broken_policy = b"\n".join(line for line in attributes.split(b"\n")
                                    if not line.startswith(b"docs/screenshots/tutorial_gui_acceptance/"))
        broken = Path(self.tmp.name) / "missing-hash-policy"
        prepare_unprotected_checkout(self.root, broken, broken_policy)
        with self.assertRaisesRegex(RuntimeError, "capture_hash_mismatch"):
            checker.check_tutorial_visuals(broken)

    def test_fenced_visual_examples_are_not_resolved_as_assets(self):
        self.write("docs/tutorial/reference.md", b'# Reference\n```markdown\n![Example](../not-a-real-asset.png)\n```\n')
        checker.check_tutorial_visuals(self.root)

    def test_invalid_source_json_is_actionable(self):
        path = self.write("docs/tutorial/sources/manual.json", b"{")
        with self.assertRaisesRegex(RuntimeError, "invalid_json.*manual.json"):
            checker.check_tutorial_visuals(self.root)

    def test_visual_failure_precedes_binary_invocation(self):
        with patch.object(checker, "check_tutorial_visuals", side_effect=RuntimeError("visual refused")), \
                patch.object(checker.subprocess, "run") as run:
            with self.assertRaisesRegex(RuntimeError, "visual refused"):
                checker.check_checkout(Path("unused"), self.root, "lf", 10)
            run.assert_not_called()


class RetainedTutorialReceiptTests(unittest.TestCase):
    def test_serpine1_occupancy_fixture_matches_retained_figure_provenance(self):
        root = checker.ROOT
        fixture = root / "test_files/fixtures/genomic_regions/serpine1_offline/synthetic_serpine1_cutrun.bed"
        figure = root / "docs/tutorial/generated/artifacts/portable_genomic_regions_offline/artifacts/portable_genomic_regions.locus.svg"
        sources = [element for element in ET.parse(figure).iter()
                   if element.get("data-gentle-occupancy-source") == "synthetic_serpine1_cutrun_v1"]
        self.assertEqual(len(sources), 1, "expected retained synthetic occupancy lane")
        self.assertEqual("sha256:" + hashlib.sha256(fixture.read_bytes()).hexdigest(),
                         sources[0].get("data-gentle-occupancy-source-sha256"),
                         "Preserve the synthetic BED input's LF checkout rule")

    def test_retained_bed_bytes_match_manifest_digests(self):
        root = checker.ROOT / "docs/tutorial/generated"
        manifests = sorted(root.rglob("*.bed.manifest.json"))
        self.assertTrue(manifests, "expected retained BED receipt coverage")
        for path in manifests:
            with self.subTest(manifest=path.relative_to(root)):
                manifest = json.loads(path.read_bytes())
                self.assertEqual(manifest["schema"], "gentle.genomic_region_bed_manifest.v1")
                name = manifest["bed_file_name"]
                self.assertEqual(Path(name).name, name)
                digest = "sha256:" + hashlib.sha256((path.parent / name).read_bytes()).hexdigest()
                self.assertEqual(digest, manifest["bed_sha256"],
                                 "Preserve the retained BED export's LF checkout rule")


if __name__ == "__main__":
    unittest.main()
