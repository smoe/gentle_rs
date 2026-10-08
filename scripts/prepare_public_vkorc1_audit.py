#!/usr/bin/env python3
"""CI-only public 08.04 inputs and CLI/shared-shell parity, without publication.

GENtle prepares the complete catalogued reference and owns all biological
operations. A retained NCBI response is replayed through its supported file URL
override, not fetched again under an unrecorded identity. No agent is invoked.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
from pathlib import Path
import os
import shutil
import subprocess
import urllib.request

try:
    from tutorial_gui_acceptance import fixed_shell_command
except ModuleNotFoundError:
    from scripts.tutorial_gui_acceptance import fixed_shell_command


CONTEXT = "vkorc1_rs9923231_context"
PREFIX = "vkorc1_rs9923231_promoter_"
WORKFLOW = "docs/examples/workflows/vkorc1_rs9923231_promoter_luciferase_assay_planning.json"
SOURCE = "docs/tutorial/sources/08-04_vkorc1_warfarin_promoter_luciferase_gui.json"
REFSNP_URL = "https://api.ncbi.nlm.nih.gov/variation/v0/beta/refsnp/9923231"


def sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def report_content(report: dict) -> dict:
    """Exclude only execution identity, not scientific fields or list ordering."""
    return {key: value for key, value in report.items()
            if key not in ("generated_at_unix_ms", "op_id", "run_id")}


def pair_content(project: dict) -> dict:
    sequences = project["sequences"]
    bases = {role: bytes(sequences[PREFIX + role]["seq"]["seq"])
             for role in ("fragment", "reference", "alternate")}
    assert bases["fragment"].upper() == bases["reference"].upper(), "reference identity changed"
    assert len(bases["reference"]) == len(bases["alternate"]), "insert length changed"
    differences = [{"position_0based": index, "reference": chr(before), "alternate": chr(after)}
                   for index, (before, after) in enumerate(zip(bases["reference"], bases["alternate"]))
                   if before != after]
    assert len(differences) == 1 and differences[0]["reference"] == "C" and differences[0]["alternate"] == "T", differences
    position = differences[0]["position_0based"]
    assert bases["fragment"][:position] == bases["reference"][:position]
    assert bases["fragment"][position + 1:] == bases["reference"][position + 1:]
    return {"length_bp": len(bases["reference"]), "differences": differences,
            "sequence_sha256": {role: hashlib.sha256(value).hexdigest() for role, value in bases.items()},
            "raw_sequences": {role: value.decode("ascii") for role, value in bases.items()},
            "native_screenshot": False, "human_scientific_approval": False}


def audit(repo: Path, binary: Path, root: Path) -> int:
    repo = repo.resolve(strict=True)
    binary = binary.resolve(strict=True)
    root = root.resolve()
    root.mkdir()  # Never overwrite a previous run or historical evidence.
    revision = subprocess.check_output(["git", "-C", str(repo), "rev-parse", "HEAD"], text=True).strip()
    version = subprocess.check_output([str(binary), "--version"], text=True)
    assert any(line.startswith("Source revision ") and line.endswith("+git." + revision)
               for line in version.splitlines()), "binary source does not match clean checkout"
    receipt = {"schema": "gentle.public_vkorc1_audit.v1", "source_revision": revision,
               "profile": "dev", "binary_sha256": sha(binary), "binary_version": version,
               "cargo_lock_sha256": sha(repo / "Cargo.lock"), "workflow_sha256": sha(repo / WORKFLOW),
               "catalog_sha256": sha(repo / "assets/genomes.json"), "commands": [], "status": "running",
               "live_agent_accepted": False, "native_gui_accepted": False, "package_accepted": False,
               "human_scientific_approval": False}
    if os.environ.get("GITHUB_RUN_ID"):
        receipt["run_url"] = "https://github.com/" + os.environ["GITHUB_REPOSITORY"] + "/actions/runs/" + os.environ["GITHUB_RUN_ID"]
    env = dict(os.environ)

    def cli(state: Path, args: list[str], label: str, expected_success: bool = True) -> subprocess.CompletedProcess:
        command = [str(binary), "--state", str(state), *args]
        with (root / f"{label}.stdout.json").open("wb") as stdout, (root / f"{label}.stderr.txt").open("wb") as stderr:
            result = subprocess.run(command, cwd=repo, env=env, stdout=stdout, stderr=stderr, timeout=7200)
        receipt["commands"].append({"label": label, "argv": command, "exit_code": result.returncode})
        assert (result.returncode == 0) == expected_success, f"inspect {label} logs"
        return result

    def op(state: Path, operation: dict, label: str, expected_success: bool = True) -> None:
        cli(state, ["op", json.dumps(operation)], label, expected_success)

    try:
        seed = root / "source.project.gentle.json"
        workflow = json.loads((repo / WORKFLOW).read_bytes())["workflow"]["ops"]
        prepare = copy.deepcopy(workflow[0])
        fetch = copy.deepcopy(workflow[1])
        for operation in (prepare, fetch):
            body = next(iter(operation.values()))
            body["cache_dir"] = str(root.parent / "vkorc1-public-reference")
            body["catalog_path"] = str(repo / "assets/genomes.json")
        op(seed, prepare, "prepare-complete-reference")
        request = urllib.request.Request(REFSNP_URL, headers={"User-Agent": "GENtle-public-audit/1"})
        with urllib.request.urlopen(request, timeout=90) as response:
            raw = response.read(16 * 1024 * 1024 + 1)
            assert len(raw) <= 16 * 1024 * 1024, "unexpectedly large refSNP response"
            receipt["refsnp_source"] = {"url": REFSNP_URL, "final_url": response.url,
                                        "http_status": response.status, "date": response.headers.get("Date")}
        document = json.loads(raw)
        assert str(document["refsnp_id"]) == "9923231", "wrong public refSNP"
        refsnp = root / "refsnp-9923231.json"
        refsnp.write_bytes(raw)
        receipt["refsnp_source"]["sha256"] = sha(refsnp)
        receipt["refsnp_source"]["execution_source"] = "retained file URL replay"
        env["GENTLE_NCBI_DBSNP_REFSNP_URL"] = refsnp.as_uri()
        op(seed, fetch, "fetch-retained-public-slice")
        initial = json.loads(seed.read_bytes())
        initial_bases = initial["sequences"][CONTEXT]["seq"]["seq"]
        routes = {route: root / f"{route}.project.gentle.json" for route in ("direct", "shared")}
        for state in routes.values():
            shutil.copy2(seed, state)
        source = json.loads((repo / SOURCE).read_bytes())
        commands = {case["id"]: case["command"] for case in source["agent_parity"]["cases"]}
        direct_ops = {"annotate_promoters": copy.deepcopy(workflow[2]),
                      "inspect_context": copy.deepcopy(workflow[3]),
                      "inspect_fragments": copy.deepcopy(workflow[4])}
        report_names = {"inspect_context": "context.json", "inspect_fragments": "fragments.json"}
        for case, operation in direct_ops.items():
            direct = operation
            shared = commands[case]
            if case in report_names:
                next(iter(direct.values()))["path"] = str(root / ("direct-" + report_names[case]))
                shared += " " + fixed_shell_command(["--path", str(root / ("shared-" + report_names[case]))])
            op(routes["direct"], direct, "direct-" + case)
            cli(routes["shared"], ["shell", shared], "shared-" + case)
        for name in report_names.values():
            direct_report = json.loads((root / ("direct-" + name)).read_bytes())
            shared_report = json.loads((root / ("shared-" + name)).read_bytes())
            assert report_content(direct_report) == report_content(shared_report), f"{name} parity changed"
        context = json.loads((root / "direct-context.json").read_bytes())
        assert context["genomic_ref"] == "C" and context["genomic_alt"] == "A,G,T"
        assert context["variant_start_0based"] == 3000 and context["promoter_overlap"] is True
        assert context["genome_anchor"]["anchor_verified"] is True
        assert context["chosen_transcript_id"] == "ENST00000498155"
        extract = next(copy.deepcopy(item) for item in workflow if "ExtractRegion" in item)
        reference = next(copy.deepcopy(item) for item in workflow
                         if item.get("MaterializeVariantAllele", {}).get("allele") == "reference")
        alternate = next(copy.deepcopy(item) for item in workflow
                         if item.get("MaterializeVariantAllele", {}).get("allele") == "alternate")
        for route, state in routes.items():
            op(state, extract, route + "-extract-reviewed-fragment")
            refused = copy.deepcopy(alternate)
            refused["MaterializeVariantAllele"].pop("alternate_allele")
            refused["MaterializeVariantAllele"]["output_id"] = "ambiguous_alternate"
            before = sha(state)
            if route == "direct":
                op(state, refused, route + "-ambiguous-refusal", False)
            else:
                line = commands["reviewed_t_insert"].replace(" --alternate-base T", "")
                line = line.replace(PREFIX + "alternate", "ambiguous_alternate")
                cli(state, ["shell", line], route + "-ambiguous-refusal", False)
            assert sha(state) == before, "refusal modified persisted project"
            error = (root / (route + "-ambiguous-refusal.stderr.txt")).read_text()
            assert "multiple alternate alleles 'A,G,T'; choose one explicitly" in error
            if route == "shared":
                shutil.copy2(state, root / "public-starter.project.gentle.json")
                cli(state, ["shell", commands["reference_insert"]], route + "-reference")
                cli(state, ["shell", commands["reviewed_t_insert"]], route + "-reviewed-t")
            else:
                op(state, reference, route + "-reference")
                op(state, alternate, route + "-reviewed-t")
        projects = {route: json.loads(path.read_bytes()) for route, path in routes.items()}
        for role in ("fragment", "reference", "alternate"):
            assert projects["direct"]["sequences"][PREFIX + role]["seq"] == projects["shared"]["sequences"][PREFIX + role]["seq"], role
        for project in projects.values():
            assert project["sequences"][CONTEXT]["seq"]["seq"] == initial_bases, "source DNA changed"
        proof = pair_content(projects["shared"])
        proof["project_sha256"] = sha(routes["shared"])
        proof["synthetic"] = False
        proof["origin"] = "GENtle extraction from complete catalogued GRCh38/Ensembl 116 and retained public NCBI refSNP response"
        (root / "base-comparison.json").write_text(json.dumps(proof, indent=2) + "\n")
        cache = root.parent / "vkorc1-public-reference"
        manifests = []
        catalog_entry = json.loads((repo / "assets/genomes.json").read_bytes())["Human GRCh38 Ensembl 116"]
        for path in cache.rglob("manifest.json"):
            manifest = json.loads(path.read_bytes())
            if manifest.get("genome_id") != "Human GRCh38 Ensembl 116":
                continue
            assert manifest["sequence_source"] == catalog_entry["sequence_remote"]
            assert manifest["annotation_source"] == catalog_entry["annotations_remote"]
            destination = root / "reference-manifests" / path.relative_to(cache)
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(path, destination)
            keys = ["sequence_path", "annotation_path", "fasta_index_path"]
            keys += [key for key in ("gene_index_path", "transcript_index_path") if manifest.get(key)]
            inputs = {key: {"sha256": sha(Path(manifest[key])), "bytes": Path(manifest[key]).stat().st_size}
                      for key in keys}
            manifests.append({"manifest": str(destination.relative_to(root)), "inputs": inputs})
        assert len(manifests) == 1, "one complete matching reference manifest required"
        receipt["reference_inputs"] = manifests
        receipt["status"] = "pass"
        receipt["claim"] = "public source-bound CLI/shared-shell promoter and explicit C/T insert parity only"
    except Exception as error:
        receipt["status"] = "fail"
        receipt["message"] = str(error)
    finally:
        receipt["artifacts"] = {str(path.relative_to(root)): sha(path)
                                for path in sorted(root.rglob("*")) if path.is_file()}
        (root / "receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps(receipt, indent=2))
    return 0 if receipt["status"] == "pass" else 1


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", required=True, type=Path)
    parser.add_argument("--gentle-cli", required=True, type=Path)
    parser.add_argument("--evidence-dir", required=True, type=Path)
    args = parser.parse_args()
    raise SystemExit(audit(args.repo_root, args.gentle_cli, args.evidence_dir))
