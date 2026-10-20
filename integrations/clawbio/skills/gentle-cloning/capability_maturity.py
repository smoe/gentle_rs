"""Read-only validation of GENtle's generated capability-maturity projection.

This module performs metadata checks, never scientific assessment, runtime
probing, permission changes or telemetry. The shared Rust contract owns scoped
evidence interpretation; discovery remains unassessed without runtime binding.
"""

from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path
from typing import Any

SCHEMA = "gentle.clawbio_capability_maturity.v1"
EVIDENCE_SCHEMA = "gentle.capability_maturity_evidence.v1"
FILENAME = "CAPABILITY_MATURITY.json"
ROOT_PREFIX = "integrations/clawbio/skills/"
KINDS = {"execution", "native_gui", "biological_review", "performance"}
LEVELS = {"unassessed", "verified", "review_required", "experimental"}
VERDICTS = {"pass", "fail", "incomplete", "withdrawn"}
UNKNOWN = {"unknown", "unavailable", "not_assessed", "not_checked", "unspecified", "any", "all"}


def digest(data: bytes) -> str:
    return "sha256:" + hashlib.sha256(data).hexdigest()


def canonical(value: Any) -> bytes:
    return json.dumps(value, ensure_ascii=True, sort_keys=True, separators=(",", ":"),
                      allow_nan=False).encode("utf-8")


def parse_json(raw: bytes) -> Any:
    """Reject duplicate keys and NaN without normalizing the retained bytes."""
    def unique(pairs):
        obj = {}
        for key, value in pairs:
            if key in obj:
                raise ValueError(f"duplicate JSON key: {key}")
            obj[key] = value
        return obj

    def invalid_constant(value):
        raise ValueError(f"invalid JSON constant: {value}")

    return json.loads(raw.decode("utf-8"), object_pairs_hook=unique,
                      parse_constant=invalid_constant)


def read_json(path: Path) -> tuple[bytes, Any]:
    """Parse and hash the same retained bytes."""
    raw = path.read_bytes()
    return raw, parse_json(raw)


def sha256(value: Any) -> bool:
    return isinstance(value, str) and re.fullmatch(r"sha256:[0-9a-f]{64}", value) is not None


def fields(value: Any, expected: set[str], context: str) -> None:
    if not isinstance(value, dict) or set(value) != expected:
        raise ValueError(f"{context}: unexpected or missing fields")


def validate_evidence(record: Any) -> None:
    fields(record, {"evidence_id", "capability_id", "definition_sha256", "scope", "kind",
                    "verdict", "maturity", "artifact", "assertions", "limitations"}, "evidence")
    if not all(isinstance(record[k], str) and record[k].strip()
               for k in ("evidence_id", "capability_id")):
        raise ValueError("evidence identity must be nonempty")
    if not sha256(record["definition_sha256"]):
        raise ValueError("invalid evidence definition digest")
    if (record["kind"] not in KINDS or record["verdict"] not in VERDICTS
            or record["maturity"] not in LEVELS):
        raise ValueError("invalid evidence kind, verdict or maturity")
    scope = record["scope"]
    fields(scope, {"source_sha", "tool_version", "binary_sha256", "platform", "architecture",
                   "build_profile", "features", "input_class", "resources"}, "scope")
    if not isinstance(scope["source_sha"], str) or not re.fullmatch(r"[0-9a-f]{40}", scope["source_sha"]):
        raise ValueError("evidence needs an exact source SHA")
    if not sha256(scope["binary_sha256"]):
        raise ValueError("evidence needs an exact binary digest")
    for key in ("tool_version", "platform", "architecture", "build_profile", "input_class"):
        if (not isinstance(scope[key], str) or not scope[key].strip() or "*" in scope[key]
                or scope[key].strip().lower() in UNKNOWN):
            raise ValueError(f"invalid evidence scope: {key}")
    features = scope["features"]
    if (not isinstance(features, list) or not all(isinstance(x, str) and x.strip() and "*" not in x for x in features)
            or features != sorted(set(features))):
        raise ValueError("features must be sorted unique strings")
    resources = scope["resources"]
    if not isinstance(resources, dict) or not all(
        isinstance(k, str) and k.strip() and "*" not in k and sha256(v) for k, v in resources.items()
    ):
        raise ValueError("resources must bind named content hashes")
    fields(record["artifact"], {"path", "sha256"}, "evidence artifact")
    if not safe_relative_path(record["artifact"]["path"]) or not sha256(record["artifact"]["sha256"]):
        raise ValueError("invalid evidence artifact binding")
    for key in ("assertions", "limitations"):
        if not isinstance(record[key], list) or not all(isinstance(x, str) and x.strip() for x in record[key]):
            raise ValueError(f"invalid evidence {key}")
    if not record["assertions"] or not record["limitations"]:
        raise ValueError("evidence needs explicit assertions and limitations")
    if record["verdict"] != "pass" and record["maturity"] != "unassessed":
        raise ValueError("non-pass evidence cannot claim verified maturity")


def safe_relative_path(value: Any) -> bool:
    return (isinstance(value, str) and bool(value) and not value.startswith("/")
            and "\\" not in value and ":" not in value
            and all(part not in ("", ".", "..") for part in value.split("/")))


def validate_catalog(catalog: Any) -> None:
    fields(catalog, {"schema", "runtime_binding", "availability", "policy", "source_files", "capabilities"}, "catalog")
    if (catalog["schema"] != SCHEMA or catalog["runtime_binding"] != "not_checked"
            or catalog["availability"] != "not_checked" or catalog["policy"] != "existing_execution_checks"):
        raise ValueError("catalog must not assert runtime binding, availability or permission")
    if not isinstance(catalog["source_files"], list) or not catalog["source_files"]:
        raise ValueError("catalog needs bound source files")
    paths = []
    for source in catalog["source_files"]:
        fields(source, {"path", "sha256"}, "source")
        if not safe_relative_path(source["path"]) or not sha256(source["sha256"]):
            raise ValueError("invalid source binding")
        paths.append(source["path"])
    if paths != sorted(set(paths)):
        raise ValueError("source files must be sorted and unique")
    if not isinstance(catalog["capabilities"], list) or not catalog["capabilities"]:
        raise ValueError("catalog must not disguise missing capabilities as an empty result")
    ids = []
    evidence_ids = set()
    for row in catalog["capabilities"]:
        fields(row, {"capability_id", "skill", "skill_version", "intent_id", "route_sha256",
                     "definition_sha256", "steps", "maturity", "evidence"}, "capability")
        for key in ("capability_id", "skill", "skill_version", "intent_id"):
            if not isinstance(row[key], str) or not row[key].strip():
                raise ValueError(f"missing capability {key}")
        if row["capability_id"] != f"{row['skill']}/{row['intent_id']}":
            raise ValueError("capability ID must bind skill and stable intent ID")
        ids.append(row["capability_id"])
        if not sha256(row["definition_sha256"]) or not sha256(row["route_sha256"]):
            raise ValueError("invalid capability definition digest")
        if row["maturity"] != "unassessed":
            raise ValueError("discovery cannot claim a runtime-scoped maturity verdict")
        steps = row["steps"]
        if not isinstance(steps, list) or not steps:
            raise ValueError("capability needs step bindings")
        for index, step in enumerate(steps):
            fields(step, {"step_index", "request_mode", "scope", "registry_refs"}, "step")
            if type(step["step_index"]) is not int or step["step_index"] != index:
                raise ValueError("step bindings must retain exact plan order")
            if step["scope"] not in {"wrapper_metadata", "runtime_command", "workflow", "runtime_unregistered"}:
                raise ValueError("invalid step scope")
            if not isinstance(step["request_mode"], str) or not step["request_mode"].strip():
                raise ValueError("invalid step request mode")
            if not isinstance(step["registry_refs"], list):
                raise ValueError("invalid registry references")
            refs = []
            for ref in step["registry_refs"]:
                fields(ref, {"source", "name"}, "registry reference")
                if ref["source"] not in {"glossary_command", "engine_operation", "mcp_tool"} or not isinstance(ref["name"], str) or not ref["name"].strip():
                    raise ValueError("invalid registry reference")
                refs.append((ref["source"], ref["name"]))
            if len(refs) != len(set(refs)) or (not refs and step["scope"] not in {"wrapper_metadata", "runtime_unregistered"}):
                raise ValueError("duplicate or missing runtime registry references")
        if not isinstance(row["evidence"], list):
            raise ValueError("evidence must be an array")
        for record in row["evidence"]:
            validate_evidence(record)
            if record["capability_id"] != row["capability_id"]:
                raise ValueError("evidence does not match capability identity")
            if record["evidence_id"] in evidence_ids:
                raise ValueError("duplicate evidence ID")
            evidence_ids.add(record["evidence_id"])
    if ids != sorted(set(ids)):
        raise ValueError("capabilities must be sorted and unique")


def skill_info(skill_root: Path) -> dict[str, Any]:
    """Optional additive discovery; a broken sidecar never disables ordinary use."""
    path = skill_root / FILENAME
    result = {"schema": SCHEMA, "path": str(path), "status": "unavailable",
              "runtime_binding": "not_checked", "availability": "not_checked",
              "source_binding": "local_skill_files_only"}
    try:
        raw, catalog = read_json(path)
        validate_catalog(catalog)
        descriptor_raw, descriptor = read_json(skill_root / "INTENTS.json")
        catalog_raw, entry = read_json(skill_root / "catalog_entry.json")
        expected_skill = entry["name"]
        reference = {"schema": SCHEMA, "path": FILENAME}
        if entry.get("capability_maturity", reference) != reference:
            raise ValueError("invalid optional maturity reference")
        if descriptor.get("schema") != "clawbio.skill_intents.v1" or descriptor.get("skill") != expected_skill:
            raise ValueError("unsupported or mixed descriptor identity")
        sources = {source["path"]: source["sha256"] for source in catalog["source_files"]}
        for name, data in (("INTENTS.json", descriptor_raw), ("catalog_entry.json", catalog_raw)):
            if sources.get(f"{ROOT_PREFIX}{expected_skill}/{name}") != digest(data):
                raise ValueError(f"stale local {name} binding")
        prefix = f"{ROOT_PREFIX}{expected_skill}/"
        for source, expected in sources.items():
            if source.startswith(prefix):
                local = skill_root / source[len(prefix):]
                if not local.resolve().is_relative_to(skill_root.resolve()) or digest(local.read_bytes()) != expected:
                    raise ValueError(f"stale local source binding: {source[len(prefix):]}")
        rows = catalog["capabilities"]
        routes = descriptor["routes"]
        if len(rows) != len(routes):
            raise ValueError("route coverage mismatch")
        by_intent = {row["intent_id"]: row for row in rows}
        for route in routes:
            row = by_intent.get(route["intent_id"])
            if (row is None or row["skill"] != expected_skill or row["skill_version"] != entry["version"]
                    or row["route_sha256"] != digest(canonical(route))
                    or row["steps"] != route["capability_steps"]):
                raise ValueError("stale or mixed route binding")
        result.update(status="available", sha256=digest(raw), capability_count=len(rows),
                      evidence_count=sum(len(row["evidence"]) for row in rows),
                      maturity="unassessed")
    except FileNotFoundError:
        result["diagnostic"] = "capability_maturity_missing"
    except (OSError, ValueError, KeyError, TypeError) as error:
        result.update(status="invalid", diagnostic="capability_maturity_invalid", detail=str(error))
    return result
