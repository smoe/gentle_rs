#!/usr/bin/env python3
"""CI-only 08.04 generation and native replay, without publishing or local builds.

Generation happens in a disposable exact-HEAD clone. Historical panel baseline
bytes and their ledger hashes are restored before checking the new projections.
Native screenshots remain raw, and the separate base view is labelled synthetic.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess


CHAPTER = "vkorc1_warfarin_promoter_luciferase_gui"
HISTORICAL = tuple(
    f"artifacts/patz1_transcript_assay_panels_cli/artifacts/{name}.report.json"
    for name in ("patz1_endpoint_end_matrix", "patz1_routine_common_region_screen",
                 "patz1_sybr_juc_panel")
)


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run_logged(command: list[str], label: str, cwd: Path, evidence: Path,
               receipt: dict) -> None:
    path = evidence / f"{label}.log"
    with path.open("wb") as log:
        result = subprocess.run(command, cwd=cwd, stdout=log, stderr=subprocess.STDOUT, timeout=1800)
    receipt["commands"].append({"label": label, "argv": command, "exit_code": result.returncode})
    if result.returncode:
        print(f"{label} failed; retained log excerpt:\n{path.read_text(errors='replace')[-12000:]}", flush=True)
        raise RuntimeError(f"{label} failed: inspect retained log")


def pair_evidence(project: dict) -> dict:
    """Fail closed on the persisted synthetic guard, not on image similarity."""
    sequences = project["sequences"]
    prefix = "vkorc1_rs9923231_promoter_"
    bases = {role: bytes(sequences[prefix + role]["seq"]["seq"])
             for role in ("fragment", "reference", "alternate")}
    assert bases["fragment"] == b"aaaaaacaaaaaaaaaaaaa", "source changed"
    reference = bytearray(bases["fragment"])
    reference[6] = ord("C")
    assert bases["reference"] == reference, "reference changed beyond materialized base case"
    assert len(bases["alternate"]) == len(bases["reference"]), "length changed"
    differences = [{"position_0based": index, "reference": chr(before), "alternate": chr(after)}
                   for index, (before, after) in enumerate(zip(bases["reference"], bases["alternate"]))
                   if before != after]
    assert differences == [{"position_0based": 6, "reference": "C", "alternate": "T"}]
    return {"schema": "gentle.synthetic_allele_pair_evidence.v1",
            "synthetic": True, "human_scientific_approval": False,
            "online_vkorc1_accepted": False, "native_screenshot": False,
            "length_bp": 20, "differences": differences,
            "sequences": {role: value.decode("ascii").upper() for role, value in bases.items()},
            "raw_sequences": {role: value.decode("ascii") for role, value in bases.items()},
            "display_case": "uppercase; raw_sequences retains GenBank case",
            "nonclaim": "20-base allele-choice guard only, not human VKORC1 DNA or a reporter assay."}


def write_pair_view(root: Path, project_path: Path) -> None:
    evidence = pair_evidence(json.loads(project_path.read_bytes()))
    evidence["project_sha256"] = sha(project_path)
    (root / "base-comparison.json").write_text(json.dumps(evidence, indent=2) + "\n")
    rows = []
    for index, role in enumerate(("reference", "alternate")):
        rows.append(f'<text x="24" y="{100 + index * 42}">{role}: {evidence["sequences"][role]}</text>')
    (root / "base-comparison.svg").write_text(
        '<svg xmlns="http://www.w3.org/2000/svg" width="760" height="260" viewBox="0 0 760 260">'
        '<rect width="760" height="260" fill="#faf8ef"/>'
        '<g fill="#242c30" font-family="monospace" font-size="18">'
        '<text x="24" y="38">Synthetic 20-base guard, not a native screenshot</text>'
        + "".join(rows)
        + '<text x="24" y="188">Position 6 (zero-based): C to T; source unchanged.</text>'
        '<text x="24" y="225">No online-locus, assay or scientific approval claim.</text>'
        '</g></svg>\n'
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, required=True)
    parser.add_argument("--binary-dir", type=Path, required=True)
    parser.add_argument("--evidence-dir", type=Path, required=True)
    args = parser.parse_args()
    repo = args.repo_root.resolve()
    binaries = args.binary_dir.resolve()
    evidence = args.evidence_dir.resolve()
    evidence.mkdir(parents=True, exist_ok=True)
    work = evidence.parent / "vkorc1-projection-checkout"
    revision = subprocess.check_output(["git", "-C", str(repo), "rev-parse", "HEAD"], text=True).strip()
    receipt = {"schema": "gentle.vkorc1_gui_audit.v1", "source_revision": revision,
               "profile": "dev", "platform": "Linux/X11", "package_accepted": False,
               "human_scientific_approval": False, "cargo_lock_sha256": sha(repo / "Cargo.lock"),
               "binaries": {name: sha(binaries / name) for name in ("gentle", "gentle_cli", "gentle_examples_docs")},
               "commands": [], "status": "running"}

    def run(command: list[str], label: str, cwd: Path = work) -> None:
        run_logged(command, label, cwd, evidence, receipt)

    try:
        run(["git", "clone", "--quiet", "--local", "--no-hardlinks", str(repo), str(work)], "clone", repo)
        generated = work / "docs/tutorial/generated"
        old_ledger = json.loads((generated / "report.json").read_bytes())
        historical = {name: (generated / name).read_bytes() for name in HISTORICAL}
        helper = str(binaries / "gentle_examples_docs")
        for mode in ("generate", "tutorial-catalog-generate", "tutorial-manifest-generate", "tutorial-generate"):
            run([helper, mode], mode)
        ledger_path = generated / "report.json"
        ledger = json.loads(ledger_path.read_bytes())
        for name, payload in historical.items():
            (generated / name).write_bytes(payload)
            ledger["file_checksums"][name] = old_ledger["file_checksums"][name]
        ledger_path.write_text(json.dumps(ledger, indent=2))
        for mode in ("--check", "tutorial-catalog-check", "tutorial-manifest-check", "tutorial-check"):
            run([helper, mode], mode.lstrip("-"))
        changed = subprocess.check_output(["git", "diff", "--name-only", "--", "docs/examples", "docs/tutorial"], cwd=work, text=True).splitlines()
        changed += subprocess.check_output(["git", "ls-files", "--others", "--exclude-standard", "--", "docs/examples", "docs/tutorial/generated"], cwd=work, text=True).splitlines()
        projection_paths = sorted(set(changed))
        for name in projection_paths:
            path = work / name
            destination = evidence / "projections" / name
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(path, destination)
        receipt["generated_projection_paths"] = projection_paths
        native = evidence / "native"
        native.mkdir()
        script = '''set -euo pipefail
        openbox > "$1/openbox.log" 2>&1 &
        wm=$!
        trap 'kill "$wm" || true' EXIT
        sleep 1
        parent_netns=$(readlink /proc/self/ns/net)
        sudo --preserve-env=DISPLAY,XAUTHORITY unshare --net -- \
          python3 "$2/scripts/tutorial_gui_acceptance.py" --repo-root "$2" \
          --gentle "$3/gentle" --gentle-cli "$3/gentle_cli" --examples-docs "$3/gentle_examples_docs" \
          --chapter vkorc1_warfarin_promoter_luciferase_gui --evidence-dir "$1" \
          --parent-network-namespace "$parent_netns" --network-enforcement linux_network_namespace
        '''
        run(["xvfb-run", "-a", "-s", "-screen 0 1600x1000x24", "bash", "-c", script,
             "audit", str(native), str(work), str(binaries)], "native-replay")
        write_pair_view(evidence, native / CHAPTER / "starter.project.gentle.json")
        receipt["status"] = "pass"
    except Exception as error:
        receipt["status"] = "fail"
        receipt["message"] = str(error)
    finally:
        receipt["artifacts"] = {str(path.relative_to(evidence)): sha(path)
                                for path in sorted(evidence.rglob("*")) if path.is_file()}
        (evidence / "receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps(receipt, indent=2))
    return 0 if receipt["status"] == "pass" else 1


if __name__ == "__main__":
    raise SystemExit(main())
