#!/usr/bin/env python3
"""Project explicit route bindings/evidence into engine and ClawBio catalogs.

No verdict is mined from parity, tutorials or release prose. This generator
copies scoped attestations only; it never runs GENtle or changes permissions.
"""

from __future__ import annotations

import argparse
import dataclasses
import importlib.util
import json
from pathlib import Path
import re
import sys

ROOT = Path(__file__).resolve().parents[1]
SKILLS = ("gentle-cloning", "gentle-pcr-primer-design", "gentle-tss-collection")
SKILL_DIR = "integrations/clawbio/skills"
CONTRACT_PATH = ROOT / SKILL_DIR / "gentle-cloning/capability_maturity.py"
spec = importlib.util.spec_from_file_location("gentle_maturity_contract", CONTRACT_PATH)
contract = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(contract)
wrapper_spec = importlib.util.spec_from_file_location(
    "gentle_maturity_wrapper", ROOT / SKILL_DIR / "gentle-cloning/gentle_cloning.py")
wrapper = importlib.util.module_from_spec(wrapper_spec)
sys.modules[wrapper_spec.name] = wrapper
assert wrapper_spec.loader is not None
wrapper_spec.loader.exec_module(wrapper)
EVIDENCE_PATH = "assets/capability_maturity_evidence.json"
OUTPUT_PATH = "assets/capability_maturity.json"


def skill_discovery(text: str) -> tuple[str, str]:
    """Check the repository's conservative plain-name/folded-description form.

    This is deliberately not a general YAML parser. Other frontmatter belongs
    to the host; changing the discovery representation needs an explicit check.
    """
    front = text.split("---", 2)
    if len(front) != 3 or front[0].strip():
        raise ValueError("SKILL.md needs frontmatter")
    match = re.search(r"^name: ([a-z0-9-]+)\ndescription: >-\n(.*?)(?=^[^\s\n]|\Z)",
                      front[1], re.M | re.S)
    if match is None:
        raise ValueError("unsupported discovery frontmatter representation")
    name, folded = match.groups()
    lines = folded.rstrip().splitlines()
    if not lines or any(not line.startswith("  ") or not line.strip() for line in lines):
        raise ValueError("discovery description must be one folded paragraph")
    description = " ".join(line.strip() for line in folded.splitlines())
    if not re.fullmatch(r"[a-z0-9]+(?:-[a-z0-9]+)*", name) or len(name) > 64:
        raise ValueError("invalid portable skill name")
    if not 1 <= len(description) <= 1024:
        raise ValueError("skill description must contain 1..1024 characters")
    return name, description


def render(value) -> bytes:
    return (json.dumps(value, sort_keys=True, ensure_ascii=True, indent=2, allow_nan=False) + "\n").encode("utf-8")


def expected_binding(request, workflow, index, glossary):
    """Check curated references against existing pure adapter argument builders."""
    mode = request["mode"]
    refs = []
    scope = "runtime_command"
    if mode in {"skill-info", "intents"}:
        scope = "wrapper_metadata"
    elif mode == "workflow" and workflow is None:
        scope = "runtime_unregistered"
    else:
        known_fields = {field.name for field in dataclasses.fields(wrapper.Request)}
        typed = wrapper.Request(**{key: value for key, value in request.items() if key in known_fields})
        if mode == "gene-protein-2d-gel":
            workflow = wrapper._gene_protein_2d_gel_workflow(typed)
        if workflow is not None:
            scope = "workflow"
            for operation in workflow.get("workflow", workflow)["ops"]:
                if not isinstance(operation, dict) or len(operation) != 1:
                    raise ValueError("invalid workflow operation")
                ref = {"source": "engine_operation", "name": next(iter(operation))}
                if ref not in refs:
                    refs.append(ref)
        else:
            args = (["shell", wrapper.DEFAULT_DEMO_SHELL_LINE] if mode == "demo"
                    else wrapper._build_cli_args(typed, ROOT / SKILL_DIR / "gentle-cloning/gentle_cloning.py"))
            if args[:1] == ["--state"]:
                args = args[2:]
            if args[0] == "op":
                refs = [{"source": "engine_operation", "name": name} for name in json.loads(args[1])]
            elif args[0] == "shell":
                line = args[1]
                paths = sorted({row["path"] for row in glossary["commands"]}, key=lambda p: (-len(p), p))
                matches = [path for path in paths if line == path or line.startswith(path + " ")]
                if matches:
                    refs = [{"source": "glossary_command", "name": matches[0]}]
                elif line.startswith("export-lab-instructions "):
                    refs = [{"source": "engine_operation", "name": "ExportLabAssistantInstructions"}]
                else:
                    scope = "runtime_unregistered"
            elif args[0] == "capabilities":
                refs = [{"source": "glossary_command", "name": "capabilities"}]
            else:
                scope = "runtime_unregistered"
    return {"step_index": index, "request_mode": mode, "scope": scope, "registry_refs": refs}


def generate(root: Path = ROOT) -> dict[str, bytes]:
    sources = {}

    def raw_source(relative):
        if not contract.safe_relative_path(relative):
            raise ValueError(f"unsafe source path: {relative}")
        path = root / relative
        if not path.resolve().is_relative_to(root.resolve()):
            raise ValueError(f"source escapes repository: {relative}")
        raw = path.read_bytes()
        sources[relative] = contract.digest(raw)
        return raw

    def read(relative):
        return contract.parse_json(raw_source(relative))

    glossary = read("docs/glossary.json")
    adapter_path = f"{SKILL_DIR}/gentle-cloning/gentle_cloning.py"
    for relative in (adapter_path, f"{SKILL_DIR}/gentle-cloning/capability_maturity.py"):
        raw_source(relative)
    archive = read(EVIDENCE_PATH)
    contract.fields(archive, {"schema", "description", "records"}, "evidence archive")
    if archive["schema"] != contract.EVIDENCE_SCHEMA or not isinstance(archive["records"], list):
        raise ValueError("invalid evidence archive")
    evidence = {}
    for record in archive["records"]:
        contract.validate_evidence(record)
        artifact = record["artifact"]
        relative = artifact["path"]
        if not (root / relative).resolve().is_relative_to(root.resolve()):
            raise ValueError("evidence artifact escapes repository")
        raw = (root / relative).read_bytes()
        if contract.digest(raw) != artifact["sha256"]:
            raise ValueError(f"evidence artifact changed: {relative}")
        sources[relative] = artifact["sha256"]
        evidence.setdefault(record["capability_id"], []).append(record)

    capabilities = []
    for skill in SKILLS:
        prefix = f"{SKILL_DIR}/{skill}"
        entry = read(f"{prefix}/catalog_entry.json")
        descriptor = read(f"{prefix}/INTENTS.json")
        reference = {"schema": contract.SCHEMA, "path": contract.FILENAME}
        if entry.get("capability_maturity", reference) != reference:
            raise ValueError("invalid optional maturity reference")
        if descriptor.get("schema") != "clawbio.skill_intents.v1":
            raise ValueError("unsupported intent descriptor schema")
        name, description = skill_discovery((root / prefix / "SKILL.md").read_text(encoding="utf-8"))
        if name != skill or entry["name"] != skill or descriptor["skill"] != skill:
            raise ValueError("mixed skill identity")
        if not 1 <= len(entry["description"]) <= 1024:
            raise ValueError("catalog discovery description exceeds portable bound")
        if skill == "gentle-cloning" and entry["description"] != description:
            raise ValueError("generic discovery descriptions diverged")
        for route in descriptor["routes"]:
            cap_id = f"{skill}/{route['intent_id']}"
            if route.get("capability_id") != cap_id:
                raise ValueError(f"missing stable capability_id: {cap_id}")
            steps = route.get("capability_steps")
            if not isinstance(steps, list) or len(steps) != len(route["plan"]):
                raise ValueError(f"missing plan-step bindings: {cap_id}")
            inputs = []
            for index, step in enumerate(route["plan"]):
                request = step.get("input_template")
                if "input" in step:
                    request = read(f"{prefix}/{step['input']}")
                if request is None and step.get("demo") is True:
                    request = {"mode": "demo"}
                if not isinstance(request, dict):
                    raise ValueError(f"missing request for {cap_id} step {index}")
                mode = request.get("mode")
                if mode != steps[index]["request_mode"]:
                    raise ValueError(f"stale request-mode binding: {cap_id} step {index}")
                workflow = request.get("workflow")
                if request.get("workflow_path") and "{" not in request["workflow_path"]:
                    workflow = read(request["workflow_path"])
                workflow_inputs = []
                if workflow is not None:
                    for relative in workflow.get("required_files", []):
                        raw_source(relative)
                        workflow_inputs.append({"path": relative, "sha256": sources[relative]})
                if steps[index] != expected_binding(request, workflow, index, glossary):
                    raise ValueError(f"stale operation/command binding: {cap_id} step {index}")
                inputs.append({"request": request, "workflow": workflow, "workflow_inputs": workflow_inputs})
            definition = {"skill": skill, "skill_version": entry["version"],
                          "route": route, "inputs": inputs, "adapter_sha256": sources[adapter_path]}
            capabilities.append({
                "capability_id": cap_id, "skill": skill, "skill_version": entry["version"],
                "intent_id": route["intent_id"], "route_sha256": contract.digest(contract.canonical(route)),
                "definition_sha256": contract.digest(contract.canonical(definition)),
                "steps": steps, "maturity": "unassessed",
                "evidence": sorted(evidence.pop(cap_id, []), key=lambda item: item["evidence_id"]),
            })
    if evidence:
        raise ValueError("evidence refers to unknown capability IDs: " + ", ".join(sorted(evidence)))
    catalog = {"schema": contract.SCHEMA, "runtime_binding": "not_checked",
               "availability": "not_checked", "policy": "existing_execution_checks",
               "source_files": [{"path": path, "sha256": sha} for path, sha in sorted(sources.items())],
               "capabilities": sorted(capabilities, key=lambda row: row["capability_id"])}
    contract.validate_catalog(catalog)
    result = {OUTPUT_PATH: render(catalog)}
    for skill in SKILLS:
        projection = dict(catalog, capabilities=[row for row in catalog["capabilities"] if row["skill"] == skill])
        result[f"{SKILL_DIR}/{skill}/{contract.FILENAME}"] = render(projection)
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--check", action="store_true")
    mode.add_argument("--write", action="store_true")
    parser.add_argument("--root", type=Path, default=ROOT)
    args = parser.parse_args()
    try:
        outputs = generate(args.root)
        stale = []
        for relative, data in outputs.items():
            path = args.root / relative
            if args.write:
                path.write_bytes(data)
            elif not path.is_file() or path.read_bytes() != data:
                stale.append(relative)
        if stale:
            raise ValueError("stale maturity projections; run scripts/generate_capability_maturity.py --write: " + ", ".join(stale))
    except (OSError, ValueError, KeyError, TypeError) as error:
        print(f"capability maturity: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
