from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import sys

import pytest


SKILL_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = SKILL_ROOT / "discussion_moderation_handoff.py"


def _module():
    spec = importlib.util.spec_from_file_location("discussion_moderation_handoff", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _canonical_hash(value: object) -> str:
    encoded = json.dumps(
        value,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    return "sha256:" + hashlib.sha256(encoded).hexdigest()


def _bytes_hash(value: bytes) -> str:
    return "sha256:" + hashlib.sha256(value).hexdigest()


def _write_provider_bundle(
    root: Path,
    *,
    wrapper_status: str = "ok",
    execution_outcome: str = "completed",
    stdout_json: object | None = None,
    request_source_sha256: str = "sha256:" + "2" * 64,
) -> tuple[Path, Path, Path, str]:
    input_path = root / "reviewed-primer-records.json"
    input_path.write_text('{"records":["fictional-primer-pair"]}\n', encoding="utf-8")
    input_hash = _bytes_hash(input_path.read_bytes())
    if stdout_json is None and execution_outcome == "completed":
        stdout_json = {
            "schema": "gentle.synthetic_primer_review.v1",
            "report_id": "primer-report-fixture",
            "status": "specificity_fail",
            "summary_lines": ["Specificity evidence contains a conflicting product."],
        }
    stdout = json.dumps(stdout_json) + "\n" if stdout_json is not None else ""
    stderr = ""
    manifest = {
        "schema": "gentle.clawbio_execution_manifest.v1",
        "started_utc": "2026-08-02T10:00:00+00:00",
        "ended_utc": "2026-08-02T10:00:01+00:00",
        "execution_outcome": execution_outcome,
        "wrapper_status": wrapper_status,
        "request_binding": {
            "normalized_request_sha256": "sha256:" + "1" * 64,
            "request_source_sha256": request_source_sha256,
            "content_bound_inputs": [
                {
                    "resolved_path": str(input_path),
                    "size_bytes": input_path.stat().st_size,
                    "sha256": input_hash,
                    "binding_ids": ["external_primer_pairs"],
                    "post_execution_verification": {
                        "exists": True,
                        "size_bytes": input_path.stat().st_size,
                        "sha256": input_hash,
                        "status": "unchanged",
                    },
                }
            ],
        },
        "delegation": {
            "source_skill": "gentle-pcr-primer-design",
            "source_skill_version": "0.3.0",
            "intent_id": "imported_primer_review",
            "plan_step_index": 0,
            "descriptor_sha256": "sha256:" + "3" * 64,
            "catalog_sha256": "sha256:" + "4" * 64,
            "route_sha256": "sha256:" + "5" * 64,
            "plan_step_sha256": "sha256:" + "6" * 64,
        },
        "wrapper": {
            "skill": "gentle-cloning",
            "skill_version": "0.2.0",
            "request_schema": "gentle.clawbio_skill_request.v1",
            "contract_files_complete": True,
            "files": [
                {
                    "contract_part": "runner",
                    "sha256": "sha256:" + "7" * 64,
                }
            ],
        },
        "runtime": {
            "gentle_version": "gentle_cli 0.1.fixture",
            "content_bound_runtime_files": [
                {
                    "binding_id": "argv_prefix_0",
                    "sha256": "sha256:" + "8" * 64,
                }
            ],
        },
        "reference_preflight": {
            "payload": {"genome_id": "fixture-reference"},
            "payload_sha256": "sha256:" + "9" * 64,
        },
        "native_result": {
            "schema": (
                stdout_json.get("schema") if isinstance(stdout_json, dict) else None
            ),
            "stdout_sha256": _bytes_hash(stdout.encode("utf-8")),
            "stderr_sha256": _bytes_hash(stderr.encode("utf-8")),
            "stdout_json_sha256": (
                _canonical_hash(stdout_json) if stdout_json is not None else None
            ),
            "identifier_bindings": [
                {
                    "json_pointer": "/stdout_json/report_id",
                    "field": "report_id",
                    "value": "primer-report-fixture",
                }
            ] if stdout_json is not None else [],
            "reported_status_bindings": [
                {
                    "json_pointer": "/stdout_json/status",
                    "field": "status",
                    "value": stdout_json.get("status"),
                }
            ] if isinstance(stdout_json, dict) and stdout_json.get("status") else [],
        },
    }
    repro = root / "reproducibility"
    repro.mkdir()
    manifest_path = repro / "execution_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    result = {
        "schema": "gentle.clawbio_skill_result.v1",
        "status": wrapper_status,
        "ended_utc": "2026-08-02T10:00:01+00:00",
        "command": ["gentle_cli", "shell", "primers show-report fixture"],
        "exit_code": 0 if wrapper_status == "ok" else 1,
        "stdout": stdout,
        "stdout_json": stdout_json,
        "stderr": stderr,
        "warnings": [],
        "error": None if wrapper_status == "ok" else "fixture execution failed",
        "failure_summary": None,
        "external_primer_handoff": None,
        "claim_ledger": {
            "schema": "gentle.clawbio_claim_ledger.v1",
            "claims": [],
        },
        "artifact_summary": None,
        "execution_manifest": manifest,
    }
    result_path = root / "result.json"
    result_path.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    return result_path, manifest_path, input_path, input_hash


def test_canonical_analysis_request_digest_vector() -> None:
    module = _module()
    payload = {
        "topic": "fixture",
        "request": {"kind": "participant", "ref": "p"},
        "routine": {"id": "r", "version": "1"},
        "inputs": [
            {
                "ref": "evidence:e",
                "content_hash": "sha256:" + "0" * 64,
            }
        ],
        "parameters": {},
        "reference_releases": {},
        "environment": {},
    }

    assert module._sha256_json(payload) == (
        "sha256:23f8c22f6973cec4756e1cf53cadcd03c4568dede59301ea9e91967dddb3cd8e"
    )


def _write_context(
    root: Path,
    input_hash: str,
    *,
    output: bool = True,
) -> Path:
    context = {
        "schema": "gentle.clawbio_discussion_moderation_handoff_context.v1",
        "topic": "fictional-primer-discussion",
        "request": {"kind": "participant", "ref": "participant-fixture"},
        "run_id": "analysis-gentle-primer-fixture-v1",
        "recorded_by": "moderator",
        "input_refs": [
            {
                "ref": "evidence:external-primer-pairs-fixture-v1",
                "content_hash": input_hash,
                "binding_id": "external_primer_pairs",
            }
        ],
        "parameters": {"review_scope": "fictional"},
        "reference_releases": {"transcripts": "fixture-release-1"},
        "environment": {"runtime_profile": "offline-fixture"},
    }
    if output:
        context["output_evidence"] = {
            "id": "gentle-primer-result-fixture-v1",
            "version": "1",
            "type_uri": "urn:gentle:primer-review-result",
            "schema_version": "1",
            "label": "Fictional GENtle primer review",
            "visibility": "shared",
            "permission": "sender_approved",
            "source_refs": ["source-fictional-primer-record"],
        }
    context_path = root / "handoff-context.json"
    context_path.write_text(json.dumps(context, indent=2) + "\n", encoding="utf-8")
    return context_path


def test_handoff_joins_exact_inputs_and_preserves_native_biological_fail(
    tmp_path: Path,
) -> None:
    module = _module()
    result_path, manifest_path, _input_path, input_hash = _write_provider_bundle(tmp_path)
    context_path = _write_context(tmp_path, input_hash)
    output_path = tmp_path / "handoff.json"

    packet = module.build_handoff(
        result_path=result_path,
        manifest_path=manifest_path,
        context_path=context_path,
        output_path=output_path,
    )

    assert packet["schema"] == "gentle.clawbio_discussion_moderation_handoff.v1"
    assert packet["verified_input_bindings"] == [
        {
            "ref": "evidence:external-primer-pairs-fixture-v1",
            "binding_id": "external_primer_pairs",
            "content_hash": input_hash,
            "provider_binding_ids": ["external_primer_pairs"],
        }
    ]
    run = packet["analysis_run_intake"]
    assert run["outcome"] == "completed"
    assert run["output_refs"] == ["evidence:gentle-primer-result-fixture-v1"]
    assert run["input_descriptors"] == [
        {
            "ref": "evidence:external-primer-pairs-fixture-v1",
            "content_hash": input_hash,
        }
    ]
    digest_payload = {
        "topic": run["topic"],
        "request": run["request"],
        "routine": run["routine"],
        "inputs": run["input_descriptors"],
        "parameters": run["parameters"],
        "reference_releases": run["reference_releases"],
        "environment": run["environment"],
    }
    assert run["request_digest"] == _canonical_hash(digest_payload)
    adapter = packet["provider_receipt"]["handoff_adapter"]
    assert adapter["id"] == "gentle-cloning.discussion-moderation-handoff"
    assert adapter["version"] == "1.0.0"
    assert adapter == run["environment"]["gentle_handoff_adapter"]
    assert adapter["sha256"] == _bytes_hash(SCRIPT.read_bytes())
    assert run["diagnostics"]["native_reported_status_bindings"][0]["value"] == (
        "specificity_fail"
    )
    evidence = packet["output_evidence_object"]
    assert evidence["claim_level"] == "computed"
    assert evidence["validation"]["status"] == "passed"
    assert evidence["validation"]["diagnostics"]["biological_status_reclassified"] is False
    assert evidence["source_refs"] == ["source-fictional-primer-record"]
    assert evidence["derived_from"] == ["evidence:external-primer-pairs-fixture-v1"]
    native = Path(packet["provider_receipt"]["native_result_artifact"]["path"])
    assert evidence["content_hash"] == _bytes_hash(native.read_bytes())
    native_payload = json.loads(native.read_text(encoding="utf-8"))
    assert native_payload["stdout_json"]["status"] == "specificity_fail"
    assert native_payload["claim_ledger"]["schema"] == (
        "gentle.clawbio_claim_ledger.v1"
    )
    assert packet["ledger_owned_fields"] == [
        "created_at",
        "recorded_at",
        "attempt",
        "freshness_status",
        "dependency_snapshot",
        "idempotent_replay",
    ]
    assert json.loads(output_path.read_text(encoding="utf-8")) == packet


def test_request_digest_is_stable_across_run_paths_and_receipt_identity(
    tmp_path: Path,
) -> None:
    module = _module()
    digests = []
    manifest_hashes = []
    for index, run_name in enumerate(("first-run", "relocated-replay")):
        root = tmp_path / run_name
        root.mkdir()
        result_path, manifest_path, _input_path, input_hash = _write_provider_bundle(
            root,
            request_source_sha256="sha256:" + ("2" if index == 0 else "a") * 64,
        )
        context_path = _write_context(root, input_hash)
        packet = module.build_handoff(
            result_path=result_path,
            manifest_path=manifest_path,
            context_path=context_path,
            output_path=root / "handoff.json",
        )
        digests.append(packet["analysis_run_intake"]["request_digest"])
        manifest_hashes.append(
            packet["provider_receipt"]["execution_manifest"]["content_hash"]
        )

    assert digests[0] == digests[1]
    assert manifest_hashes[0] != manifest_hashes[1]


def test_handoff_rejects_internally_inconsistent_input_receipt(tmp_path: Path) -> None:
    module = _module()
    result_path, manifest_path, _input_path, input_hash = _write_provider_bundle(tmp_path)
    context_path = _write_context(tmp_path, input_hash)
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["request_binding"]["content_bound_inputs"][0][
        "post_execution_verification"
    ]["sha256"] = "sha256:" + "f" * 64
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    result = json.loads(result_path.read_text(encoding="utf-8"))
    result["execution_manifest"] = manifest
    result_path.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")

    with pytest.raises(module.HandoffError, match="was not stable"):
        module.build_handoff(
            result_path=result_path,
            manifest_path=manifest_path,
            context_path=context_path,
            output_path=tmp_path / "inconsistent.json",
        )


def test_handoff_rejects_input_hash_mismatch_and_tampered_native_output(
    tmp_path: Path,
) -> None:
    module = _module()
    result_path, manifest_path, _input_path, input_hash = _write_provider_bundle(tmp_path)
    context_path = _write_context(tmp_path, "sha256:" + "0" * 64)

    with pytest.raises(module.HandoffError, match="content hash does not match"):
        module.build_handoff(
            result_path=result_path,
            manifest_path=manifest_path,
            context_path=context_path,
            output_path=tmp_path / "mismatch.json",
        )

    context_path = _write_context(tmp_path, input_hash)
    result = json.loads(result_path.read_text(encoding="utf-8"))
    result["stdout_json"]["status"] = "specificity_pass"
    result_path.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    with pytest.raises(module.HandoffError, match="parsed JSON hash mismatch"):
        module.build_handoff(
            result_path=result_path,
            manifest_path=manifest_path,
            context_path=context_path,
            output_path=tmp_path / "tampered.json",
        )


def test_failed_run_has_no_output_evidence_or_vacuous_success(tmp_path: Path) -> None:
    module = _module()
    result_path, manifest_path, _input_path, input_hash = _write_provider_bundle(
        tmp_path,
        wrapper_status="failed",
        execution_outcome="failed",
        stdout_json=None,
    )
    context_path = _write_context(tmp_path, input_hash, output=False)

    packet = module.build_handoff(
        result_path=result_path,
        manifest_path=manifest_path,
        context_path=context_path,
        output_path=tmp_path / "failed-handoff.json",
    )

    assert packet["analysis_run_intake"]["outcome"] == "failed"
    assert packet["analysis_run_intake"]["output_refs"] == []
    assert packet["output_evidence_object"] is None
    assert packet["provider_receipt"]["native_result_artifact"] is None


def test_handoff_rejects_malformed_warning_metadata(tmp_path: Path) -> None:
    module = _module()
    result_path, manifest_path, _input_path, input_hash = _write_provider_bundle(tmp_path)
    context_path = _write_context(tmp_path, input_hash)
    result = json.loads(result_path.read_text(encoding="utf-8"))
    result["warnings"] = ["valid warning", {"unexpected": "object"}]
    result_path.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")

    with pytest.raises(module.HandoffError, match="warnings must be an array of strings"):
        module.build_handoff(
            result_path=result_path,
            manifest_path=manifest_path,
            context_path=context_path,
            output_path=tmp_path / "malformed-warning.json",
        )


def test_cli_writes_same_packet_and_rejects_sensitive_context(tmp_path: Path) -> None:
    result_path, manifest_path, _input_path, input_hash = _write_provider_bundle(tmp_path)
    context_path = _write_context(tmp_path, input_hash)
    output_path = tmp_path / "cli-handoff.json"
    run = subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--result",
            str(result_path),
            "--manifest",
            str(manifest_path),
            "--context",
            str(context_path),
            "--output",
            str(output_path),
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert run.returncode == 0, run.stderr
    assert json.loads(run.stdout) == json.loads(output_path.read_text(encoding="utf-8"))

    context = json.loads(context_path.read_text(encoding="utf-8"))
    context["environment"] = {"chat_id": "must-not-cross-provider-boundary"}
    context_path.write_text(json.dumps(context) + "\n", encoding="utf-8")
    rejected = subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--result",
            str(result_path),
            "--manifest",
            str(manifest_path),
            "--context",
            str(context_path),
            "--output",
            str(tmp_path / "rejected.json"),
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert rejected.returncode == 1
    assert "disallowed metadata key: chat_id" in rejected.stderr
