#!/usr/bin/env python3
"""Build a content-bound discussion-moderation intake packet from a GENtle run."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import re
import sys
import tempfile
from typing import Any


CONTEXT_SCHEMA = "gentle.clawbio_discussion_moderation_handoff_context.v1"
PACKET_SCHEMA = "gentle.clawbio_discussion_moderation_handoff.v1"
NATIVE_RESULT_SCHEMA = "gentle.clawbio_native_result_artifact.v1"
RESULT_SCHEMA = "gentle.clawbio_skill_result.v1"
MANIFEST_SCHEMA = "gentle.clawbio_execution_manifest.v1"
HANDOFF_ADAPTER_ID = "gentle-cloning.discussion-moderation-handoff"
HANDOFF_ADAPTER_VERSION = "1.0.0"
MAX_CONTEXT_BYTES = 256 * 1024
MAX_RESULT_BYTES = 128 * 1024 * 1024
MAX_METADATA_BYTES = 32 * 1024

_HASH_RE = re.compile(r"^sha256:[0-9a-f]{64}$")
_TYPE_URI_RE = re.compile(r"^[A-Za-z][A-Za-z0-9+.-]*:")
_TYPED_REF_RE = re.compile(r"^(?:evidence|evidence_list):[^\s:][^\s]*$")
_SENSITIVE_METADATA_KEY_RE = re.compile(
    r"(?:password|passwd|secret|token|credential|api[_-]?key|private[_-]?key|"
    r"chat[_-]?id|conversation[_-]?id|identity[_-]?(?:hash|anchor))",
    re.IGNORECASE,
)


class HandoffError(ValueError):
    """Raised when a provider bundle cannot be joined safely to caller context."""


def _reject_non_finite_json_constant(value: str) -> None:
    raise HandoffError(f"non-finite JSON constant: {value}")


def _ensure_finite(value: Any, path: str = "$") -> None:
    if isinstance(value, float) and not math.isfinite(value):
        raise HandoffError(f"non-finite JSON number at {path}")
    if isinstance(value, dict):
        for key, nested in value.items():
            _ensure_finite(nested, f"{path}.{key}")
    elif isinstance(value, list):
        for index, nested in enumerate(value):
            _ensure_finite(nested, f"{path}[{index}]")


def _canonical_json_bytes(value: Any) -> bytes:
    _ensure_finite(value)
    return json.dumps(
        value,
        allow_nan=False,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")


def _sha256_bytes(value: bytes) -> str:
    return "sha256:" + hashlib.sha256(value).hexdigest()


def _sha256_json(value: Any) -> str:
    return _sha256_bytes(_canonical_json_bytes(value))


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return "sha256:" + digest.hexdigest()


def _read_json(path: Path, *, max_bytes: int, label: str) -> dict[str, Any]:
    try:
        size = path.stat().st_size
    except FileNotFoundError as error:
        raise HandoffError(f"{label} does not exist: {path}") from error
    if not path.is_file():
        raise HandoffError(f"{label} is not a regular file: {path}")
    if size > max_bytes:
        raise HandoffError(f"{label} exceeds the {max_bytes}-byte limit")
    try:
        payload = json.loads(
            path.read_text(encoding="utf-8"),
            parse_constant=_reject_non_finite_json_constant,
        )
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise HandoffError(f"{label} is not valid UTF-8 JSON: {error}") from error
    if not isinstance(payload, dict):
        raise HandoffError(f"{label} must contain a JSON object")
    _ensure_finite(payload)
    return payload


def _atomic_write(path: Path, data: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary_name = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    temporary_path = Path(temporary_name)
    try:
        with os.fdopen(fd, "wb") as handle:
            handle.write(data)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary_path, path)
    finally:
        if temporary_path.exists():
            temporary_path.unlink()


def _reject_unknown_fields(value: dict[str, Any], allowed: set[str], label: str) -> None:
    unknown = sorted(set(value) - allowed)
    if unknown:
        raise HandoffError(f"{label} contains unsupported fields: {', '.join(unknown)}")


def _required_string(value: Any, label: str, *, max_length: int = 1024) -> str:
    if not isinstance(value, str) or not value.strip():
        raise HandoffError(f"{label} must be a non-empty string")
    normalized = value.strip()
    if len(normalized) > max_length or any(ord(character) < 32 for character in normalized):
        raise HandoffError(f"{label} must be a bounded printable string")
    return normalized


def _validate_hash(value: Any, label: str) -> str:
    normalized = _required_string(value, label, max_length=71)
    if not _HASH_RE.fullmatch(normalized):
        raise HandoffError(f"{label} must use sha256:<64 lowercase hex>")
    return normalized


def _string_list(value: Any, label: str) -> list[str]:
    if value is None:
        return []
    if not isinstance(value, list):
        raise HandoffError(f"{label} must be an array")
    result = [_required_string(item, f"{label}[{index}]") for index, item in enumerate(value)]
    if len(result) != len(set(result)):
        raise HandoffError(f"{label} must not contain duplicates")
    return result


def _safe_metadata(value: Any, label: str) -> dict[str, Any]:
    if value is None:
        return {}
    if not isinstance(value, dict):
        raise HandoffError(f"{label} must be a JSON object")
    encoded = _canonical_json_bytes(value)
    if len(encoded) > MAX_METADATA_BYTES:
        raise HandoffError(f"{label} exceeds the {MAX_METADATA_BYTES}-byte limit")

    def walk(item: Any) -> None:
        if isinstance(item, dict):
            for key, nested in item.items():
                if _SENSITIVE_METADATA_KEY_RE.search(str(key)):
                    raise HandoffError(f"{label} contains disallowed metadata key: {key}")
                walk(nested)
        elif isinstance(item, list):
            for nested in item:
                walk(nested)

    walk(value)
    return value


def _normalise_context(payload: dict[str, Any]) -> dict[str, Any]:
    _reject_unknown_fields(
        payload,
        {
            "schema",
            "topic",
            "request",
            "run_id",
            "recorded_by",
            "input_refs",
            "output_evidence",
            "parameters",
            "reference_releases",
            "environment",
        },
        "handoff context",
    )
    if payload.get("schema") != CONTEXT_SCHEMA:
        raise HandoffError(f"handoff context schema must be {CONTEXT_SCHEMA}")
    request = payload.get("request")
    if not isinstance(request, dict):
        raise HandoffError("handoff context request must be an object")
    _reject_unknown_fields(request, {"kind", "ref"}, "handoff context request")
    request_kind = _required_string(request.get("kind"), "request.kind")
    if request_kind not in {"participant", "moderator_action"}:
        raise HandoffError("request.kind must be participant or moderator_action")
    normalized_request = {
        "kind": request_kind,
        "ref": _required_string(request.get("ref"), "request.ref"),
    }

    raw_inputs = payload.get("input_refs")
    if not isinstance(raw_inputs, list) or not raw_inputs:
        raise HandoffError("handoff context input_refs must be a non-empty array")
    input_refs: list[dict[str, str]] = []
    seen_refs: set[str] = set()
    for index, raw in enumerate(raw_inputs):
        label = f"input_refs[{index}]"
        if not isinstance(raw, dict):
            raise HandoffError(f"{label} must be an object")
        _reject_unknown_fields(raw, {"ref", "content_hash", "binding_id"}, label)
        ref = _required_string(raw.get("ref"), f"{label}.ref")
        if not _TYPED_REF_RE.fullmatch(ref):
            raise HandoffError(f"{label}.ref must be evidence:<id> or evidence_list:<id>")
        if ref in seen_refs:
            raise HandoffError("handoff context input_refs must not contain duplicate refs")
        seen_refs.add(ref)
        input_refs.append(
            {
                "ref": ref,
                "content_hash": _validate_hash(raw.get("content_hash"), f"{label}.content_hash"),
                "binding_id": _required_string(
                    raw.get("binding_id"), f"{label}.binding_id", max_length=128
                ),
            }
        )

    output = payload.get("output_evidence")
    normalized_output = None
    if output is not None:
        if not isinstance(output, dict):
            raise HandoffError("output_evidence must be an object when present")
        _reject_unknown_fields(
            output,
            {
                "id",
                "version",
                "type_uri",
                "schema_version",
                "label",
                "visibility",
                "permission",
                "artifact_locator",
                "source_refs",
            },
            "output_evidence",
        )
        type_uri = _required_string(output.get("type_uri"), "output_evidence.type_uri")
        if not _TYPE_URI_RE.match(type_uri):
            raise HandoffError("output_evidence.type_uri must be an absolute URI")
        visibility = _required_string(output.get("visibility"), "output_evidence.visibility")
        if visibility not in {"shared", "private_summary", "operator_only"}:
            raise HandoffError("output_evidence.visibility has an unsupported value")
        permission = _required_string(output.get("permission"), "output_evidence.permission")
        if permission not in {"public", "sender_approved", "restricted", "unknown"}:
            raise HandoffError("output_evidence.permission has an unsupported value")
        normalized_output = {
            "id": _required_string(output.get("id"), "output_evidence.id"),
            "version": _required_string(output.get("version"), "output_evidence.version"),
            "type_uri": type_uri,
            "schema_version": _required_string(
                output.get("schema_version"), "output_evidence.schema_version"
            ),
            "label": _required_string(output.get("label"), "output_evidence.label"),
            "visibility": visibility,
            "permission": permission,
            "artifact_locator": (
                _required_string(
                    output.get("artifact_locator"), "output_evidence.artifact_locator"
                )
                if output.get("artifact_locator") is not None
                else None
            ),
            "source_refs": _string_list(output.get("source_refs"), "output_evidence.source_refs"),
        }

    return {
        "schema": CONTEXT_SCHEMA,
        "topic": _required_string(payload.get("topic"), "topic"),
        "request": normalized_request,
        "run_id": _required_string(payload.get("run_id"), "run_id"),
        "recorded_by": _required_string(payload.get("recorded_by"), "recorded_by"),
        "input_refs": input_refs,
        "output_evidence": normalized_output,
        "parameters": _safe_metadata(payload.get("parameters"), "parameters"),
        "reference_releases": _safe_metadata(
            payload.get("reference_releases"), "reference_releases"
        ),
        "environment": _safe_metadata(payload.get("environment"), "environment"),
    }


def _verify_provider_bundle(
    result_path: Path,
    manifest_path: Path,
) -> tuple[dict[str, Any], dict[str, Any]]:
    result = _read_json(result_path, max_bytes=MAX_RESULT_BYTES, label="GENtle result")
    manifest = _read_json(
        manifest_path, max_bytes=MAX_RESULT_BYTES, label="GENtle execution manifest"
    )
    if result.get("schema") != RESULT_SCHEMA:
        raise HandoffError(f"GENtle result schema must be {RESULT_SCHEMA}")
    if manifest.get("schema") != MANIFEST_SCHEMA:
        raise HandoffError(f"GENtle execution manifest schema must be {MANIFEST_SCHEMA}")
    if result.get("execution_manifest") != manifest:
        raise HandoffError("embedded and on-disk GENtle execution manifests differ")
    if result.get("status") != manifest.get("wrapper_status"):
        raise HandoffError("GENtle result status disagrees with the execution manifest")

    request_binding = manifest.get("request_binding")
    if not isinstance(request_binding, dict):
        raise HandoffError("GENtle execution manifest has no request_binding object")
    _validate_hash(
        request_binding.get("normalized_request_sha256"),
        "request_binding.normalized_request_sha256",
    )
    request_source_sha256 = request_binding.get("request_source_sha256")
    if request_source_sha256 is not None:
        _validate_hash(request_source_sha256, "request_binding.request_source_sha256")
    if not isinstance(request_binding.get("content_bound_inputs"), list):
        raise HandoffError("GENtle execution manifest has no content_bound_inputs array")

    wrapper = manifest.get("wrapper")
    if not isinstance(wrapper, dict):
        raise HandoffError("GENtle execution manifest has no wrapper identity")
    _required_string(wrapper.get("skill"), "wrapper.skill")
    _required_string(wrapper.get("skill_version"), "wrapper.skill_version")
    if wrapper.get("contract_files_complete") is not True:
        raise HandoffError("GENtle wrapper contract files were incomplete")

    native = manifest.get("native_result")
    if not isinstance(native, dict):
        raise HandoffError("GENtle execution manifest has no native_result object")
    stdout = result.get("stdout")
    stderr = result.get("stderr")
    if not isinstance(stdout, str) or not isinstance(stderr, str):
        raise HandoffError("GENtle result stdout and stderr must be strings")
    if native.get("stdout_sha256") != _sha256_bytes(stdout.encode("utf-8")):
        raise HandoffError("GENtle result stdout hash mismatch")
    if native.get("stderr_sha256") != _sha256_bytes(stderr.encode("utf-8")):
        raise HandoffError("GENtle result stderr hash mismatch")
    expected_json_hash = (
        _sha256_json(result.get("stdout_json"))
        if result.get("stdout_json") is not None
        else None
    )
    if native.get("stdout_json_sha256") != expected_json_hash:
        raise HandoffError("GENtle result parsed JSON hash mismatch")
    return result, manifest


def _bound_input_record(manifest: dict[str, Any], binding_id: str) -> dict[str, Any]:
    request_binding = manifest.get("request_binding")
    records = (
        request_binding.get("content_bound_inputs")
        if isinstance(request_binding, dict)
        else None
    )
    if not isinstance(records, list):
        raise HandoffError("GENtle manifest has no content-bound input records")
    matches = [
        record
        for record in records
        if isinstance(record, dict) and binding_id in record.get("binding_ids", [])
    ]
    if len(matches) != 1:
        raise HandoffError(
            f"caller input binding_id '{binding_id}' did not resolve exactly once"
        )
    record = matches[0]
    binding_ids = record.get("binding_ids")
    if (
        not isinstance(binding_ids, list)
        or not all(isinstance(value, str) and value for value in binding_ids)
        or len(binding_ids) != len(set(binding_ids))
    ):
        raise HandoffError(
            f"caller input binding_id '{binding_id}' has invalid binding_ids"
        )
    digest = _validate_hash(
        record.get("sha256"), f"input binding '{binding_id}' sha256"
    )
    size_bytes = record.get("size_bytes")
    if isinstance(size_bytes, bool) or not isinstance(size_bytes, int) or size_bytes < 0:
        raise HandoffError(f"caller input binding_id '{binding_id}' has invalid size_bytes")
    verification = record.get("post_execution_verification")
    if (
        not isinstance(verification, dict)
        or verification.get("status") != "unchanged"
        or verification.get("exists") is not True
        or verification.get("size_bytes") != size_bytes
        or verification.get("sha256") != digest
    ):
        raise HandoffError(f"caller input binding_id '{binding_id}' was not stable")
    return record


def _provider_parameters(
    context: dict[str, Any], manifest: dict[str, Any]
) -> dict[str, Any]:
    parameters = dict(context["parameters"])
    if "gentle_provider" in parameters:
        raise HandoffError("parameters.gentle_provider is reserved for verified provider data")
    request_binding = manifest["request_binding"]
    delegation = manifest.get("delegation")
    parameters["gentle_provider"] = {
        "request_schema": manifest.get("wrapper", {}).get("request_schema"),
        "normalized_request_sha256": request_binding.get("normalized_request_sha256"),
        "delegation": (
            {
                key: delegation.get(key)
                for key in (
                    "source_skill",
                    "source_skill_version",
                    "intent_id",
                    "plan_step_index",
                    "descriptor_sha256",
                    "catalog_sha256",
                    "route_sha256",
                    "plan_step_sha256",
                )
            }
            if isinstance(delegation, dict)
            else None
        ),
    }
    return _safe_metadata(parameters, "joined parameters")


def _provider_reference_releases(
    context: dict[str, Any], manifest: dict[str, Any]
) -> dict[str, Any]:
    releases = dict(context["reference_releases"])
    if "gentle_reference_preflight" in releases:
        raise HandoffError(
            "reference_releases.gentle_reference_preflight is reserved for provider data"
        )
    preflight = manifest.get("reference_preflight")
    if isinstance(preflight, dict) and preflight.get("payload_sha256") is not None:
        releases["gentle_reference_preflight"] = {
            "payload_sha256": preflight.get("payload_sha256")
        }
    return _safe_metadata(releases, "joined reference releases")


def _provider_environment(
    context: dict[str, Any], manifest: dict[str, Any]
) -> dict[str, Any]:
    environment = dict(context["environment"])
    reserved = {"gentle_provider", "gentle_handoff_adapter"} & set(environment)
    if reserved:
        raise HandoffError(
            "environment contains reserved verified-provider fields: "
            + ", ".join(sorted(reserved))
        )
    wrapper = manifest.get("wrapper") if isinstance(manifest.get("wrapper"), dict) else {}
    runtime = manifest.get("runtime") if isinstance(manifest.get("runtime"), dict) else {}
    environment["gentle_provider"] = {
        "execution_manifest_schema": manifest.get("schema"),
        "wrapper": {
            "skill": wrapper.get("skill"),
            "version": wrapper.get("skill_version"),
            "contract_files_complete": wrapper.get("contract_files_complete"),
            "files": [
                {
                    "contract_part": record.get("contract_part"),
                    "sha256": record.get("sha256"),
                }
                for record in wrapper.get("files", [])
                if isinstance(record, dict)
            ],
        },
        "runtime": {
            "gentle_version": runtime.get("gentle_version"),
            "files": [
                {
                    "binding_id": record.get("binding_id"),
                    "sha256": record.get("sha256"),
                }
                for record in runtime.get("content_bound_runtime_files", [])
                if isinstance(record, dict)
            ],
        },
        "state_before_sha256": manifest.get("state_binding", {})
        .get("before", {})
        .get("sha256"),
    }
    environment["gentle_handoff_adapter"] = {
        "id": HANDOFF_ADAPTER_ID,
        "version": HANDOFF_ADAPTER_VERSION,
        "sha256": _sha256_file(Path(__file__).resolve()),
    }
    return _safe_metadata(environment, "joined environment")


def _native_result_payload(
    result: dict[str, Any], manifest: dict[str, Any], manifest_sha256: str
) -> dict[str, Any]:
    return {
        "schema": NATIVE_RESULT_SCHEMA,
        "execution_manifest_sha256": manifest_sha256,
        "wrapper_status": result.get("status"),
        "command": result.get("command"),
        "exit_code": result.get("exit_code"),
        "stdout_json": result.get("stdout_json"),
        "stdout_text": result.get("stdout") if result.get("stdout_json") is None else None,
        "stderr_sha256": manifest["native_result"].get("stderr_sha256"),
        "claim_ledger": result.get("claim_ledger"),
        "artifact_summary": result.get("artifact_summary"),
        "identifier_bindings": manifest["native_result"].get("identifier_bindings", []),
        "reported_status_bindings": manifest["native_result"].get(
            "reported_status_bindings", []
        ),
        "external_primer_handoff": result.get("external_primer_handoff"),
    }


def build_handoff(
    *,
    result_path: Path,
    context_path: Path,
    output_path: Path,
    native_result_path: Path | None = None,
    manifest_path: Path | None = None,
) -> dict[str, Any]:
    result_path = result_path.resolve()
    context_path = context_path.resolve()
    output_path = output_path.resolve()
    manifest_path = (
        manifest_path.resolve()
        if manifest_path is not None
        else result_path.parent / "reproducibility" / "execution_manifest.json"
    )
    context_payload = _read_json(
        context_path, max_bytes=MAX_CONTEXT_BYTES, label="handoff context"
    )
    context = _normalise_context(context_payload)
    result, manifest = _verify_provider_bundle(result_path, manifest_path)
    outcome = manifest.get("execution_outcome")
    if outcome not in {"completed", "failed", "incomplete"}:
        raise HandoffError("GENtle execution outcome is unsupported")
    output_config = context["output_evidence"]
    if outcome == "completed" and output_config is None:
        raise HandoffError("completed GENtle runs require output_evidence context")
    if outcome == "failed" and output_config is not None:
        raise HandoffError("failed GENtle runs must not claim output_evidence")

    input_descriptors = []
    input_bindings = []
    for caller_input in context["input_refs"]:
        bound = _bound_input_record(manifest, caller_input["binding_id"])
        if bound.get("sha256") != caller_input["content_hash"]:
            raise HandoffError(
                f"caller ref '{caller_input['ref']}' content hash does not match "
                f"binding_id '{caller_input['binding_id']}'"
            )
        input_descriptors.append(
            {"ref": caller_input["ref"], "content_hash": caller_input["content_hash"]}
        )
        input_bindings.append(
            {
                "ref": caller_input["ref"],
                "binding_id": caller_input["binding_id"],
                "content_hash": caller_input["content_hash"],
                "provider_binding_ids": bound.get("binding_ids", []),
            }
        )

    manifest_sha256 = _sha256_file(manifest_path)
    result_sha256 = _sha256_file(result_path)
    context_sha256 = _sha256_file(context_path)
    routine = {
        "id": _required_string(manifest.get("wrapper", {}).get("skill"), "wrapper.skill"),
        "version": _required_string(
            manifest.get("wrapper", {}).get("skill_version"), "wrapper.skill_version"
        ),
    }
    parameters = _provider_parameters(context, manifest)
    reference_releases = _provider_reference_releases(context, manifest)
    environment = _provider_environment(context, manifest)
    request_digest = _sha256_json(
        {
            "topic": context["topic"],
            "request": context["request"],
            "routine": routine,
            "inputs": input_descriptors,
            "parameters": parameters,
            "reference_releases": reference_releases,
            "environment": environment,
        }
    )

    evidence_object = None
    output_refs: list[str] = []
    native_artifact = None
    if output_config is not None:
        native_result_path = (
            native_result_path.resolve()
            if native_result_path is not None
            else output_path.with_name(output_path.stem + ".native-result.json")
        )
        native_bytes = _canonical_json_bytes(
            _native_result_payload(result, manifest, manifest_sha256)
        )
        _atomic_write(native_result_path, native_bytes)
        native_sha256 = _sha256_bytes(native_bytes)
        artifact_locator = output_config["artifact_locator"] or str(native_result_path)
        output_ref = f"evidence:{output_config['id']}"
        output_refs = [output_ref]
        validation_status = "passed" if outcome == "completed" else "incomplete"
        recorded_at = _required_string(result.get("ended_utc"), "result.ended_utc")
        dependency_refs = sorted(
            set(context_input["ref"] for context_input in context["input_refs"])
            | set(f"source:{source_ref}" for source_ref in output_config["source_refs"])
        )
        evidence_object = {
            "id": output_config["id"],
            "topic": context["topic"],
            "version": output_config["version"],
            "content_hash": native_sha256,
            "type_uri": output_config["type_uri"],
            "schema_version": output_config["schema_version"],
            "created_at": recorded_at,
            "created_by": routine["id"],
            "creator_kind": "routine",
            "recorded_at": recorded_at,
            "recorded_by": context["recorded_by"],
            "visibility": output_config["visibility"],
            "permission": output_config["permission"],
            "retention": "external",
            "claim_level": "computed",
            "validation": {
                "status": validation_status,
                "validator": routine,
                "recorded_at": recorded_at,
                "diagnostics": {
                    "scope": "provider_integrity_and_execution_only",
                    "native_reported_status_bindings": manifest["native_result"].get(
                        "reported_status_bindings", []
                    ),
                    "biological_status_reclassified": False,
                },
            },
            "lifecycle_status": "active",
            "source_refs": output_config["source_refs"],
            "derived_from": [context_input["ref"] for context_input in context["input_refs"]],
            "dependency_refs": dependency_refs,
            "media_type": "application/json",
            "artifact": {
                "locator": artifact_locator,
                "media_type": "application/json",
                "content_hash": native_sha256,
            },
            "label": output_config["label"],
        }
        native_artifact = {
            "path": str(native_result_path),
            "locator": artifact_locator,
            "content_hash": native_sha256,
            "media_type": "application/json",
        }

    raw_warnings = result.get("warnings")
    if raw_warnings is None:
        warnings = []
    elif not isinstance(raw_warnings, list) or not all(
        isinstance(warning, str) for warning in raw_warnings
    ):
        raise HandoffError("GENtle result warnings must be an array of strings")
    else:
        warnings = raw_warnings
    diagnostics = _safe_metadata(
        {
            "provider_wrapper_status": result.get("status"),
            "provider_error": result.get("error"),
            "native_reported_status_bindings": manifest["native_result"].get(
                "reported_status_bindings", []
            ),
            "status_rule": (
                "execution_outcome_is_separate_from_native_biological_status"
            ),
        },
        "analysis-run diagnostics",
    )
    analysis_run = {
        "id": context["run_id"],
        "topic": context["topic"],
        "completed_at": _required_string(result.get("ended_utc"), "result.ended_utc"),
        "recorded_by": context["recorded_by"],
        "request": context["request"],
        "routine": routine,
        "input_refs": [item["ref"] for item in context["input_refs"]],
        "input_descriptors": input_descriptors,
        "request_digest": request_digest,
        "parameters": parameters,
        "reference_releases": reference_releases,
        "environment": environment,
        "outcome": outcome,
        "warnings": warnings,
        "output_refs": output_refs,
        "diagnostics": diagnostics,
    }
    packet = {
        "schema": PACKET_SCHEMA,
        "consumer_contract": {
            "name": "discussion-moderation record-analysis-run v3 intake",
            "request_digest_algorithm": (
                "sha256 of canonical JSON over topic, request, routine, inputs, "
                "parameters, reference_releases, and environment"
            ),
            "status": "ledger_must_recompute_and_validate",
        },
        "provider_receipt": {
            "handoff_adapter": environment["gentle_handoff_adapter"],
            "request_binding": {
                "normalized_request_sha256": manifest["request_binding"].get(
                    "normalized_request_sha256"
                ),
                "request_source_sha256": manifest["request_binding"].get(
                    "request_source_sha256"
                ),
            },
            "result": {
                "path": str(result_path),
                "schema": result.get("schema"),
                "content_hash": result_sha256,
            },
            "execution_manifest": {
                "path": str(manifest_path),
                "schema": manifest.get("schema"),
                "content_hash": manifest_sha256,
            },
            "caller_context": {
                "path": str(context_path),
                "schema": context.get("schema"),
                "content_hash": context_sha256,
            },
            "native_result_artifact": native_artifact,
        },
        "verified_input_bindings": input_bindings,
        "output_evidence_object": evidence_object,
        "analysis_run_intake": analysis_run,
        "recording_order": [
            "verify caller permissions and current input refs",
            "record output_evidence_object when present",
            "record analysis_run_intake and require the ledger to recompute request_digest",
        ],
        "ledger_owned_fields": [
            "created_at",
            "recorded_at",
            "attempt",
            "freshness_status",
            "dependency_snapshot",
            "idempotent_replay",
        ],
        "interpretation_boundary": (
            "This packet is provider evidence, not participant support, consensus, "
            "permission, or scientific interpretation."
        ),
    }
    packet_bytes = json.dumps(packet, indent=2, ensure_ascii=True).encode("utf-8") + b"\n"
    _atomic_write(output_path, packet_bytes)
    return packet


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Create a content-bound discussion-moderation intake packet from an "
            "existing gentle-cloning result bundle and explicit caller context."
        )
    )
    parser.add_argument("--result", required=True, type=Path)
    parser.add_argument("--context", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--native-result-output", type=Path)
    return parser


def main() -> int:
    args = _parser().parse_args()
    try:
        packet = build_handoff(
            result_path=args.result,
            context_path=args.context,
            output_path=args.output,
            native_result_path=args.native_result_output,
            manifest_path=args.manifest,
        )
    except HandoffError as error:
        print(f"error: {error}", file=sys.stderr)
        return 1
    print(json.dumps(packet, indent=2, ensure_ascii=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
