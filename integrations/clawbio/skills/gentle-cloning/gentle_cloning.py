#!/usr/bin/env python3
"""ClawBio/OpenClaw wrapper skill for deterministic GENtle CLI execution."""

from __future__ import annotations

import argparse
import base64
import dataclasses
import datetime as dt
import hashlib
import html
import json
import math
import os
from pathlib import Path
import platform
import re
import shlex
import shutil
import signal
import subprocess
import sys
import threading
from typing import Any
import xml.etree.ElementTree as ET

REQUEST_SCHEMA = "gentle.clawbio_skill_request.v1"
RESULT_SCHEMA = "gentle.clawbio_skill_result.v1"
EXECUTION_MANIFEST_SCHEMA = "gentle.clawbio_execution_manifest.v1"
EXECUTION_PROPOSAL_SCHEMA = "gentle.clawbio_execution_proposal.v1"
EXECUTION_APPROVAL_SCHEMA = "gentle.clawbio_execution_approval.v1"
APPROVED_EXECUTION_REQUEST_SCHEMA = "gentle.clawbio_approved_execution_request.v1"
DELEGATION_SCHEMA = "gentle.clawbio_skill_delegation.v1"
CLAIM_LEDGER_SCHEMA = "gentle.clawbio_claim_ledger.v1"
SKILL_CONTRACT_VERSION = "0.3.0"
APPROVAL_ENVIRONMENT_VARIABLES = (
    "PATH",
    "GENTLE_ASSET_ROOT",
    "GENTLE_PROJECT_ROOT",
    "GENTLE_SYSTEM_CONFIG_ROOT",
    "GENTLE_REFERENCE_CACHE_DIR",
    "GENTLE_HELPER_CACHE_DIR",
    "GENTLE_CUTRUN_CACHE_DIR",
    "GENTLE_BLASTN_BIN",
    "GENTLE_MAKEBLASTDB_BIN",
    "GENTLE_BLASTDBCMD_BIN",
    "GENTLE_BIGWIG_TO_BEDGRAPH_BIN",
    "GENTLE_RNAFOLD_BIN",
    "GENTLE_RNAPKIN_BIN",
    "GENTLE_SHA1_TOOL",
)
STRICT_CLAIM_ATTRIBUTION_MODE = "strict"
PCR_PRIMER_PRESENTATION_PROFILE = "pcr_primer_design"
CLAIM_SOURCE_PREFIXES = {
    "gentle_executable": "[gentle]",
    "clawbio_presentation": "[clawbio]",
    "caller_input": "[input]",
    "external_tool": "[external:<tool>]",
    "unverified": "[unverified]",
}
SKILL_INFO_SCHEMA = "gentle.clawbio_skill_info.v1"
INTENTS_RUNTIME_SCHEMA = "gentle.clawbio_skill_intents_runtime.v1"
EXTERNAL_PRIMER_HANDOFF_REQUEST_SCHEMA = "gentle.external_primer_handoff_request.v1"
EXTERNAL_PRIMER_HANDOFF_RESULT_SCHEMA = "gentle.external_primer_handoff_result.v1"
EXTERNAL_PRIMER_PAIR_BATCH_SCHEMA = "gentle.external_primer_pair_batch.v1"
EXTERNAL_PRIMER_PAIR_IMPORT_COMMAND_SCHEMA = (
    "gentle.external_primer_pair_import_command.v1"
)
EXTERNAL_PRIMER_PAIR_IMPORT_REPORT_SCHEMA = (
    "gentle.external_primer_pair_import_report.v1"
)
CDNA_ASSAY_TEST_REPORT_SCHEMA = "gentle.cdna_assay_test_report.v1"
PRIMER_SPECIFICITY_REPORT_SCHEMA = "gentle.primer_specificity_report.v4"
SKILL_NAME = "gentle-cloning"
INVOCATION_MARKER = "GENtle ClawBio skill wrapper invoked"
UI_INTENT_CATALOG_SCHEMA = "gentle.ui_intents.v1"
UI_INTENT_DISCOVERY_SHELL_LINE = "ui intents"
SUPPORTED_REQUEST_MODES = (
    "skill-info",
    "intents",
    "version",
    "capabilities",
    "state-summary",
    "shell",
    "op",
    "workflow",
    "construct-reasoning-list-inspections",
    "construct-reasoning-run-inspection",
    "primer-preflight",
    "primer-seed-from-feature",
    "primer-seed-from-splicing",
    "primer-design",
    "primer-report-list",
    "primer-report-show",
    "primer-report-export",
    "qpcr-seed-from-feature",
    "qpcr-seed-from-splicing",
    "qpcr-design",
    "qpcr-report-list",
    "qpcr-report-show",
    "qpcr-report-export",
    "cdna-pcr-test",
    "cdna-qpcr-test",
    "external-primer-handoff",
    "protein-residue-genomic-coordinates",
    "transcript-qpcr-panel",
    "restriction-cloning-pcr-handoff",
    "restriction-cloning-pcr-handoff-seed",
    "restriction-cloning-vector-suggestions",
    "restriction-cloning-handoff-list",
    "restriction-cloning-handoff-show",
    "restriction-cloning-handoff-export",
    "pcr-protocol-cartoon",
    "gene-protein-2d-gel",
    "exon-skip-plan",
    "exon-skip-materialize",
    "agent-plan",
    "agent-execute-plan",
    "raw",
)
PRIMER_SHELL_REQUEST_MODES = {
    "primer-preflight",
    "primer-seed-from-feature",
    "primer-seed-from-splicing",
    "primer-design",
    "primer-report-list",
    "primer-report-show",
    "primer-report-export",
    "qpcr-seed-from-feature",
    "qpcr-seed-from-splicing",
    "qpcr-design",
    "qpcr-report-list",
    "qpcr-report-show",
    "qpcr-report-export",
    "cdna-pcr-test",
    "cdna-qpcr-test",
    "restriction-cloning-pcr-handoff",
    "restriction-cloning-pcr-handoff-seed",
    "restriction-cloning-vector-suggestions",
    "restriction-cloning-handoff-list",
    "restriction-cloning-handoff-show",
    "restriction-cloning-handoff-export",
}
DEPRECATED_REQUEST_MODE_EQUIVALENTS = {
    "primer-preflight": "mode=shell with `primers preflight ...`",
    "primer-seed-from-feature": "mode=shell with `primers seed-from-feature ...`",
    "primer-seed-from-splicing": "mode=shell with `primers seed-from-splicing ...`",
    "primer-design": "mode=shell with `primers design ...`",
    "primer-report-list": "mode=shell with `primers list-reports ...`",
    "primer-report-show": "mode=shell with `primers show-report ...`",
    "primer-report-export": "mode=shell with `primers export-report ...`",
    "qpcr-seed-from-feature": "mode=shell with `primers seed-qpcr-from-feature ...`",
    "qpcr-seed-from-splicing": "mode=shell with `primers seed-qpcr-from-splicing ...`",
    "qpcr-design": "mode=shell with `primers design-qpcr ...`",
    "qpcr-report-list": "mode=shell with `primers list-qpcr-reports ...`",
    "qpcr-report-show": "mode=shell with `primers show-qpcr-report ...`",
    "qpcr-report-export": "mode=shell with `primers export-qpcr-report ...`",
    "cdna-pcr-test": "mode=shell with `primers test-cdna-pcr ...`",
    "cdna-qpcr-test": "mode=shell with `primers test-cdna-qpcr ...`",
    "transcript-qpcr-panel": "mode=shell with `primers transcript-qpcr-panel ...`",
    "restriction-cloning-pcr-handoff": "mode=shell with `primers prepare-restriction-cloning ...`",
    "restriction-cloning-pcr-handoff-seed": "mode=shell with `primers seed-restriction-cloning-handoff ...`",
    "restriction-cloning-vector-suggestions": "mode=shell with `primers restriction-cloning-vector-suggestions ...`",
    "restriction-cloning-handoff-list": "mode=shell with `primers list-restriction-cloning-handoffs ...`",
    "restriction-cloning-handoff-show": "mode=shell with `primers show-restriction-cloning-handoff ...`",
    "restriction-cloning-handoff-export": "mode=shell with `primers export-restriction-cloning-handoff ...`",
    "pcr-protocol-cartoon": "mode=shell with `protocol-cartoon render-svg ...`",
    "exon-skip-plan": "mode=shell with `transcripts exon-skip-plan ...`",
    "exon-skip-materialize": "mode=shell with `transcripts exon-skip-materialize ...`",
}
SVG_DIMENSION_RE = re.compile(r"^\s*([0-9]+(?:\.[0-9]+)?)")
CLAWBIO_GRAPHICS_SCALE = 2.0
CLAWBIO_SVG_PNG_TIMEOUT_SECS = 300
DEFAULT_DEMO_SHELL_LINE = (
    "protocol-cartoon render-svg gibson.two_fragment "
    "artifacts/gibson.two_fragment.protocol.svg"
)
DEFAULT_DEMO_EXPECTED_ARTIFACT = "artifacts/gibson.two_fragment.protocol.svg"
LOCAL_RUNTIME_LINEAGE = "GENtle Rust rewrite used by ClawBio"
VERSION_SCOPE = "installed_local_clawbio_runtime"
CLASSICAL_GENTLE_DISAMBIGUATION = (
    "This skill reports the locally installed ClawBio GENtle rewrite runtime, "
    "not the classical GENtle desktop release line."
)
CONTINUE_ARTIFACT_NOTICE = (
    'More figures were generated; ask "Continue" or inspect report.md/result.json.'
)
EXTERNAL_TOOL_RESOURCES = [
    {
        "resource_id": "vienna_rna",
        "display_name": "ViennaRNA RNAfold",
        "env_var": "GENTLE_RNAFOLD_BIN",
        "default_executable": "RNAfold",
        "status_shell_line": "resources status",
        "used_for": "RNA secondary-structure folding, dot-bracket output, and MFE reporting.",
    },
    {
        "resource_id": "rnapkin",
        "display_name": "rnapkin RNA structure renderer",
        "env_var": "GENTLE_RNAPKIN_BIN",
        "default_executable": "rnapkin",
        "status_shell_line": "resources status",
        "used_for": "Rendering ViennaRNA/RNAfold dot-bracket structures as SVG/PNG graphics.",
    },
]

EXTERNAL_PRIMER_HANDOFF_SUBMITTED_PURPOSES = {"qpcr", "endpoint_pcr"}
EXTERNAL_PRIMER_HANDOFF_NOT_SUBMITTED_PURPOSES = {"cloning", "sequencing"}
EXTERNAL_PRIMER_HANDOFF_SOURCE_KINDS = {
    "external",
    "commercial_catalogue",
    "literature",
    "laboratory",
}
EXTERNAL_PRIMER_HANDOFF_ANNOTATION_RECORD_ID = "gentle.handoff.record_id"
EXTERNAL_PRIMER_HANDOFF_ANNOTATION_ASSAY_PURPOSE = "gentle.handoff.assay_purpose"
EXTERNAL_PRIMER_HANDOFF_ANNOTATION_COLLECTION_ID = "gentle.handoff.collection_id"
EXTERNAL_PRIMER_HANDOFF_RESERVED_ANNOTATIONS = {
    EXTERNAL_PRIMER_HANDOFF_ANNOTATION_RECORD_ID,
    EXTERNAL_PRIMER_HANDOFF_ANNOTATION_ASSAY_PURPOSE,
    EXTERNAL_PRIMER_HANDOFF_ANNOTATION_COLLECTION_ID,
}
EXTERNAL_PRIMER_HANDOFF_IUPAC = frozenset("ACGTRYSWKMBDHVN")


class SkillError(RuntimeError):
    """Base skill error with deterministic message formatting."""


class ApprovalRequired(SkillError):
    """Stop before scientific execution after emitting a reviewable proposal."""


@dataclasses.dataclass
class Request:
    mode: str
    timeout_secs: int = 180
    state_path: str | None = None
    raw_args: list[str] | None = None
    shell_line: str | None = None
    operation: Any = None
    request_json: Any = None
    workflow: Any = None
    workflow_path: str | None = None
    system_id: str | None = None
    prompt: str | None = None
    catalog_path: str | None = None
    base_url: str | None = None
    model: str | None = None
    connect_timeout_secs: int | None = None
    read_timeout_secs: int | None = None
    max_retries: int | None = None
    max_response_bytes: int | None = None
    include_state_summary: bool | None = None
    max_candidates: int | None = None
    allow_mutating_candidates: bool | None = None
    plan: Any = None
    plan_path: str | None = None
    plan_id: str | None = None
    graph_id: str | None = None
    fact_id: str | None = None
    annotation_id: str | None = None
    candidate_id: str | None = None
    candidate_ids: list[str] | None = None
    evidence_id: str | None = None
    summary_id: str | None = None
    action_kind: str | None = None
    action_id: str | None = None
    confirm: bool | None = None
    transcript_feature_id: int | None = None
    skip_candidate_ids: list[str] | None = None
    skip_intervals_1based: list[Any] | None = None
    overlap_intervals_1based: list[Any] | None = None
    length_mod3_values: list[int] | None = None
    coding_mod3_values: list[int] | None = None
    coding_contexts: list[str] | None = None
    cds_phase_entry_kinds: list[str] | None = None
    feature_query: Any = None
    return_items: list[str] | None = None
    expected_artifacts: list[str] | None = None
    ensure_reference_prepared: Any = None
    external_primer_handoff: Any = None
    gene_symbol: str | None = None
    species: str | None = None
    source: str | None = None
    ladders: list[str] | None = None
    seq_id: str | None = None
    source_feature_id: int | None = None
    transcript_id: str | None = None
    transcript_order: str | None = None
    map_coordinate_mode: str | None = None
    qpcr_mode: str | None = None
    specificity_evidence: str | None = None
    backend: str | None = None
    primer3_executable: str | None = None
    report_id: str | None = None
    forward_primer: str | None = None
    reverse_primer: str | None = None
    probe: str | None = None
    min_amplicon_bp: int | None = None
    max_amplicon_bp: int | None = None
    max_mismatches: int | None = None
    word_size: int | None = None
    step_bp: int | None = None
    tile_bp: int | None = None
    dotplot_id: str | None = None
    render_svg_path: str | None = None
    require_3prime_exact_bases: int | None = None
    vector_seq_id: str | None = None
    pair_rank: int | None = None
    handoff_mode: str | None = None
    forward_enzyme: str | None = None
    reverse_enzyme: str | None = None
    forward_leader_5prime: str | None = None
    reverse_leader_5prime: str | None = None
    protocol_id: str | None = None
    shared_qpcr_report_id: str | None = None
    output_prefix: str | None = None
    output_path: str | None = None
    svg_path: str | None = None
    materialize_products: bool | None = None
    product_output_prefix: str | None = None
    product_gel_svg_path: str | None = None
    product_gel_ladders: list[str] | None = None
    residue_start_1based: int | None = None
    residue_end_1based: int | None = None
    claim_attribution_mode: str | None = None
    presentation_profile: str | None = None
    input_claims: list[str] | None = None
    delegation: Any = None
    input_bindings: list[Any] | None = None


@dataclasses.dataclass
class EnsureReferencePrepared:
    genome_id: str
    catalog_path: str | None = None
    cache_dir: str | None = None
    status_timeout_secs: int = 300
    prepare_timeout_secs: int = 7200


@dataclasses.dataclass
class CliResolution:
    argv_prefix: list[str]
    cwd: str | None
    label: str


@dataclasses.dataclass
class ApprovedExecution:
    proposal_path: Path
    proposal: dict[str, Any]
    approval: dict[str, Any]
    envelope_path: Path


def _format_command_text(command: list[str] | None) -> str:
    if not command:
        return "(none)"
    return " ".join(shlex.quote(v) for v in command)


def _one_line_preview(text: str, *, max_chars: int = 240) -> str | None:
    stripped = " ".join((text or "").strip().split())
    if not stripped:
        return None
    if len(stripped) <= max_chars:
        return stripped
    return stripped[: max_chars - 1].rstrip() + "..."


def _now_utc_iso() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat()


def _read_json(path: Path) -> dict[str, Any]:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except FileNotFoundError as e:
        raise SkillError(f"request JSON file '{path}' does not exist") from e
    except json.JSONDecodeError as e:
        raise SkillError(f"invalid JSON in '{path}': {e}") from e


def _canonical_json_bytes(value: Any) -> bytes:
    return json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
    ).encode("utf-8")


def _sha256_prefixed_bytes(value: bytes) -> str:
    return "sha256:" + hashlib.sha256(value).hexdigest()


def _sha256_prefixed_json(value: Any) -> str:
    return _sha256_prefixed_bytes(_canonical_json_bytes(value))


def _reject_unknown_fields(
    value: dict[str, Any], allowed: set[str], context: str
) -> None:
    unknown = sorted(set(value) - allowed)
    if unknown:
        raise SkillError(
            f"{context} contains unsupported field(s): {', '.join(unknown)}"
        )


def _required_handoff_string(value: Any, field_name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise SkillError(f"{field_name} must be a non-empty string")
    return value.strip()


def _optional_handoff_string(value: Any, field_name: str) -> str | None:
    if value is None:
        return None
    if not isinstance(value, str):
        raise SkillError(f"{field_name} must be a string when present")
    trimmed = value.strip()
    return trimmed or None


def _handoff_string_list(value: Any, field_name: str) -> list[str]:
    if value is None:
        return []
    if not isinstance(value, list):
        raise SkillError(f"{field_name} must be a string array when present")
    if not all(isinstance(item, str) for item in value):
        raise SkillError(f"{field_name} must contain only strings")
    return sorted({item.strip() for item in value if item.strip()})


def _handoff_annotations(value: Any, field_name: str) -> dict[str, str]:
    if value is None:
        return {}
    if not isinstance(value, dict):
        raise SkillError(f"{field_name} must be a string-to-string object")
    annotations: dict[str, str] = {}
    for raw_key, raw_value in value.items():
        if not isinstance(raw_key, str) or not raw_key.strip():
            raise SkillError(f"{field_name} keys must be non-empty strings")
        if not isinstance(raw_value, str):
            raise SkillError(f"{field_name}.{raw_key} must be a string")
        key = raw_key.strip()
        if key in EXTERNAL_PRIMER_HANDOFF_RESERVED_ANNOTATIONS:
            raise SkillError(
                f"{field_name}.{key} is reserved for the GENtle handoff join"
            )
        annotations[key] = raw_value
    return dict(sorted(annotations.items()))


def _normalise_handoff_sequence(value: Any, field_name: str) -> str:
    if not isinstance(value, str):
        raise SkillError(f"{field_name} must be a 5'-to-3' IUPAC sequence")
    normalized: list[str] = []
    for index, character in enumerate(value):
        if character.isspace() or (
            character.isascii() and character.isdigit()
        ):
            continue
        if not character.isascii():
            raise SkillError(
                f"{field_name} contains a non-ASCII character at position {index + 1}"
            )
        upper = character.upper()
        if upper == "U":
            upper = "T"
        if upper not in EXTERNAL_PRIMER_HANDOFF_IUPAC:
            raise SkillError(
                f"{field_name} contains invalid IUPAC character "
                f"'{character}' at position {index + 1}"
            )
        normalized.append(upper)
    if not normalized:
        raise SkillError(f"{field_name} contains no IUPAC bases after normalization")
    return "".join(normalized)


def _normalise_expected_sha256(value: Any, field_name: str) -> str | None:
    raw = _optional_handoff_string(value, field_name)
    if raw is None:
        return None
    digest = raw[7:] if raw.lower().startswith("sha256:") else raw
    if len(digest) != 64 or any(ch not in "0123456789abcdefABCDEF" for ch in digest):
        raise SkillError(f"{field_name} must be a SHA-256 digest")
    return "sha256:" + digest.lower()


def _normalise_delegation(value: Any) -> dict[str, Any] | None:
    if value is None:
        return None
    if not isinstance(value, dict):
        raise SkillError("delegation must be an object when present")
    _reject_unknown_fields(
        value,
        {
            "schema",
            "source_skill",
            "source_skill_version",
            "intent_id",
            "plan_step_index",
            "resolved_slots",
            "descriptor_sha256",
            "catalog_sha256",
        },
        "delegation",
    )
    if value.get("schema") != DELEGATION_SCHEMA:
        raise SkillError(
            f"delegation.schema must be '{DELEGATION_SCHEMA}'"
        )

    def identifier(field: str, *, max_len: int = 128) -> str:
        raw = _required_handoff_string(value.get(field), f"delegation.{field}")
        if len(raw) > max_len or not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]*", raw):
            raise SkillError(
                f"delegation.{field} must be a bounded portable identifier"
            )
        return raw

    source_skill = identifier("source_skill")
    source_skill_version = identifier("source_skill_version", max_len=64)
    intent_id = identifier("intent_id")
    raw_step_index = value.get("plan_step_index", 0)
    if isinstance(raw_step_index, bool):
        raise SkillError("delegation.plan_step_index must be a non-negative integer")
    try:
        plan_step_index = int(raw_step_index)
    except (TypeError, ValueError) as e:
        raise SkillError(
            "delegation.plan_step_index must be a non-negative integer"
        ) from e
    if plan_step_index < 0:
        raise SkillError("delegation.plan_step_index must be a non-negative integer")
    resolved_slots = value.get("resolved_slots")
    if resolved_slots is not None:
        if not isinstance(resolved_slots, dict):
            raise SkillError("delegation.resolved_slots must be an object when present")
        if len(_canonical_json_bytes(resolved_slots)) > 64 * 1024:
            raise SkillError("delegation.resolved_slots exceeds the 64 KiB limit")
    return {
        "schema": DELEGATION_SCHEMA,
        "source_skill": source_skill,
        "source_skill_version": source_skill_version,
        "intent_id": intent_id,
        "plan_step_index": plan_step_index,
        "resolved_slots": resolved_slots,
        "descriptor_sha256": _normalise_expected_sha256(
            value.get("descriptor_sha256"), "delegation.descriptor_sha256"
        ),
        "catalog_sha256": _normalise_expected_sha256(
            value.get("catalog_sha256"), "delegation.catalog_sha256"
        ),
    }


def _normalise_input_bindings(value: Any) -> list[dict[str, Any]] | None:
    if value is None:
        return None
    if not isinstance(value, list):
        raise SkillError("input_bindings must be an array when present")
    normalized: list[dict[str, Any]] = []
    seen_ids: set[str] = set()
    for index, raw in enumerate(value):
        context = f"input_bindings[{index}]"
        if not isinstance(raw, dict):
            raise SkillError(f"{context} must be an object")
        _reject_unknown_fields(
            raw,
            {"binding_id", "path", "role", "media_type", "expected_sha256"},
            context,
        )
        binding_id = _required_handoff_string(
            raw.get("binding_id"), f"{context}.binding_id"
        )
        if len(binding_id) > 128 or not re.fullmatch(
            r"[A-Za-z0-9][A-Za-z0-9._-]*", binding_id
        ):
            raise SkillError(f"{context}.binding_id must be a portable identifier")
        if binding_id in seen_ids:
            raise SkillError(f"input_bindings contains duplicate binding_id '{binding_id}'")
        seen_ids.add(binding_id)
        normalized.append(
            {
                "binding_id": binding_id,
                "path": _required_handoff_string(raw.get("path"), f"{context}.path"),
                "role": _optional_handoff_string(raw.get("role"), f"{context}.role"),
                "media_type": _optional_handoff_string(
                    raw.get("media_type"), f"{context}.media_type"
                ),
                "expected_sha256": _normalise_expected_sha256(
                    raw.get("expected_sha256"), f"{context}.expected_sha256"
                ),
            }
        )
    return normalized


def _normalise_handoff_nonnegative_int(
    value: Any, field_name: str, *, positive: bool = False
) -> int | None:
    if value is None:
        return None
    if isinstance(value, bool):
        raise SkillError(f"{field_name} must be an integer")
    try:
        parsed = int(value)
    except (TypeError, ValueError) as e:
        raise SkillError(f"{field_name} must be an integer") from e
    if positive and parsed <= 0:
        raise SkillError(f"{field_name} must be > 0")
    if not positive and parsed < 0:
        raise SkillError(f"{field_name} must be >= 0")
    return parsed


def _normalise_external_primer_handoff_record(
    value: Any, index: int
) -> dict[str, Any]:
    context = f"external_primer_handoff.records[{index}]"
    if not isinstance(value, dict):
        raise SkillError(f"{context} must be an object")
    _reject_unknown_fields(
        value,
        {
            "record_id",
            "assay_purpose",
            "source_kind",
            "provider",
            "catalogue_id",
            "source_url",
            "claimed_accession",
            "aliases",
            "claimed_target",
            "validation_claims",
            "annotations",
            "forward_sequence_5_to_3",
            "reverse_sequence_5_to_3",
            "sequence_5_to_3",
            "role",
        },
        context,
    )
    record_id = _required_handoff_string(value.get("record_id"), f"{context}.record_id")
    assay_purpose = _required_handoff_string(
        value.get("assay_purpose"), f"{context}.assay_purpose"
    ).lower()
    allowed_purposes = (
        EXTERNAL_PRIMER_HANDOFF_SUBMITTED_PURPOSES
        | EXTERNAL_PRIMER_HANDOFF_NOT_SUBMITTED_PURPOSES
    )
    if assay_purpose not in allowed_purposes:
        raise SkillError(
            f"{context}.assay_purpose must be one of: "
            + ", ".join(sorted(allowed_purposes))
        )
    source_kind = _required_handoff_string(
        value.get("source_kind", "external"), f"{context}.source_kind"
    ).lower()
    if source_kind not in EXTERNAL_PRIMER_HANDOFF_SOURCE_KINDS:
        raise SkillError(
            f"{context}.source_kind must be one of: "
            + ", ".join(sorted(EXTERNAL_PRIMER_HANDOFF_SOURCE_KINDS))
        )
    provider = _required_handoff_string(value.get("provider"), f"{context}.provider")
    catalogue_id = _optional_handoff_string(
        value.get("catalogue_id"), f"{context}.catalogue_id"
    )
    if source_kind == "commercial_catalogue" and catalogue_id is None:
        raise SkillError(
            f"{context}.catalogue_id is required for commercial_catalogue records"
        )

    standalone_raw = value.get("sequence_5_to_3")
    forward_raw = value.get("forward_sequence_5_to_3")
    reverse_raw = value.get("reverse_sequence_5_to_3")
    has_standalone = standalone_raw is not None
    has_pair = forward_raw is not None or reverse_raw is not None
    if has_standalone and has_pair:
        raise SkillError(
            f"{context} must use either sequence_5_to_3 or a forward/reverse pair, not both"
        )
    if has_pair and (forward_raw is None or reverse_raw is None):
        raise SkillError(
            f"{context} requires both forward_sequence_5_to_3 and "
            "reverse_sequence_5_to_3"
        )
    if not has_standalone and not has_pair:
        raise SkillError(
            f"{context} requires sequence_5_to_3 or a forward/reverse sequence pair"
        )
    if assay_purpose in EXTERNAL_PRIMER_HANDOFF_SUBMITTED_PURPOSES and not has_pair:
        raise SkillError(
            f"{context} uses assay_purpose={assay_purpose} and therefore requires "
            "an explicit forward/reverse primer pair"
        )

    normalized: dict[str, Any] = {
        "record_id": record_id,
        "assay_purpose": assay_purpose,
        "source_kind": source_kind,
        "provider": provider,
        "catalogue_id": catalogue_id or "",
        "source_url": _optional_handoff_string(
            value.get("source_url"), f"{context}.source_url"
        )
        or "",
        "claimed_accession": _optional_handoff_string(
            value.get("claimed_accession"), f"{context}.claimed_accession"
        )
        or "",
        "aliases": _handoff_string_list(value.get("aliases"), f"{context}.aliases"),
        "claimed_target": _optional_handoff_string(
            value.get("claimed_target"), f"{context}.claimed_target"
        )
        or "",
        "validation_claims": _handoff_string_list(
            value.get("validation_claims"), f"{context}.validation_claims"
        ),
        "annotations": _handoff_annotations(
            value.get("annotations"), f"{context}.annotations"
        ),
        "role": _optional_handoff_string(value.get("role"), f"{context}.role"),
    }
    if has_standalone:
        normalized["record_kind"] = "oligo"
        normalized["sequence_5_to_3"] = _normalise_handoff_sequence(
            standalone_raw, f"{context}.sequence_5_to_3"
        )
    else:
        normalized["record_kind"] = "primer_pair"
        normalized["forward_sequence_5_to_3"] = _normalise_handoff_sequence(
            forward_raw, f"{context}.forward_sequence_5_to_3"
        )
        normalized["reverse_sequence_5_to_3"] = _normalise_handoff_sequence(
            reverse_raw, f"{context}.reverse_sequence_5_to_3"
        )
    return normalized


def _normalise_external_primer_handoff(value: Any) -> dict[str, Any]:
    context = "external_primer_handoff"
    if not isinstance(value, dict):
        raise SkillError(f"{context} must be an object")
    _reject_unknown_fields(
        value,
        {"schema", "collection_id", "target", "evaluation", "records"},
        context,
    )
    if value.get("schema") != EXTERNAL_PRIMER_HANDOFF_REQUEST_SCHEMA:
        raise SkillError(
            f"{context}.schema must be {EXTERNAL_PRIMER_HANDOFF_REQUEST_SCHEMA}"
        )
    collection_id = _required_handoff_string(
        value.get("collection_id"), f"{context}.collection_id"
    )

    target = value.get("target")
    if not isinstance(target, dict):
        raise SkillError(f"{context}.target must be an object")
    _reject_unknown_fields(
        target,
        {
            "seq_id",
            "source_feature_id",
            "transcript_id",
            "reference_label",
            "reference_release",
            "expected_state_sha256",
        },
        f"{context}.target",
    )
    source_feature_id = _normalise_handoff_nonnegative_int(
        target.get("source_feature_id"), f"{context}.target.source_feature_id"
    )
    if source_feature_id is None:
        raise SkillError(f"{context}.target.source_feature_id is required")
    normalized_target = {
        "seq_id": _required_handoff_string(
            target.get("seq_id"), f"{context}.target.seq_id"
        ),
        "source_feature_id": source_feature_id,
        "transcript_id": _optional_handoff_string(
            target.get("transcript_id"), f"{context}.target.transcript_id"
        ),
        "reference_label": _optional_handoff_string(
            target.get("reference_label"), f"{context}.target.reference_label"
        ),
        "reference_release": _optional_handoff_string(
            target.get("reference_release"), f"{context}.target.reference_release"
        ),
        "expected_state_sha256": _normalise_expected_sha256(
            target.get("expected_state_sha256"),
            f"{context}.target.expected_state_sha256",
        ),
    }

    evaluation = value.get("evaluation", {})
    if not isinstance(evaluation, dict):
        raise SkillError(f"{context}.evaluation must be an object when present")
    _reject_unknown_fields(
        evaluation,
        {
            "min_amplicon_bp",
            "max_amplicon_bp",
            "max_mismatches",
            "require_3prime_exact_bases",
            "transcript_order",
            "map_coordinate_mode",
            "specificity_target_genome_id",
            "specificity_catalog_path",
            "specificity_cache_dir",
            "materialize_products",
            "product_gel_ladders",
        },
        f"{context}.evaluation",
    )
    normalized_evaluation = {
        "min_amplicon_bp": _normalise_handoff_nonnegative_int(
            evaluation.get("min_amplicon_bp"),
            f"{context}.evaluation.min_amplicon_bp",
            positive=True,
        ),
        "max_amplicon_bp": _normalise_handoff_nonnegative_int(
            evaluation.get("max_amplicon_bp"),
            f"{context}.evaluation.max_amplicon_bp",
            positive=True,
        ),
        "max_mismatches": _normalise_handoff_nonnegative_int(
            evaluation.get("max_mismatches"),
            f"{context}.evaluation.max_mismatches",
        ),
        "require_3prime_exact_bases": _normalise_handoff_nonnegative_int(
            evaluation.get("require_3prime_exact_bases"),
            f"{context}.evaluation.require_3prime_exact_bases",
        ),
        "transcript_order": _optional_handoff_string(
            evaluation.get("transcript_order"),
            f"{context}.evaluation.transcript_order",
        ),
        "map_coordinate_mode": _optional_handoff_string(
            evaluation.get("map_coordinate_mode"),
            f"{context}.evaluation.map_coordinate_mode",
        ),
        "specificity_target_genome_id": _optional_handoff_string(
            evaluation.get("specificity_target_genome_id"),
            f"{context}.evaluation.specificity_target_genome_id",
        ),
        "specificity_catalog_path": _optional_handoff_string(
            evaluation.get("specificity_catalog_path"),
            f"{context}.evaluation.specificity_catalog_path",
        ),
        "specificity_cache_dir": _optional_handoff_string(
            evaluation.get("specificity_cache_dir"),
            f"{context}.evaluation.specificity_cache_dir",
        ),
        "materialize_products": evaluation.get("materialize_products", False),
        "product_gel_ladders": _handoff_string_list(
            evaluation.get("product_gel_ladders"),
            f"{context}.evaluation.product_gel_ladders",
        ),
    }
    if not isinstance(normalized_evaluation["materialize_products"], bool):
        raise SkillError(f"{context}.evaluation.materialize_products must be boolean")
    if normalized_evaluation["transcript_order"] not in {
        None,
        "transcript_id",
        "genomic_first_exon",
        "genomic_last_exon",
        "antisense_first_exon",
    }:
        raise SkillError(f"{context}.evaluation.transcript_order is unsupported")
    if normalized_evaluation["map_coordinate_mode"] not in {
        None,
        "cdna",
        "genomic_aligned",
    }:
        raise SkillError(f"{context}.evaluation.map_coordinate_mode is unsupported")
    if (
        normalized_evaluation["min_amplicon_bp"] is not None
        and normalized_evaluation["max_amplicon_bp"] is not None
        and normalized_evaluation["min_amplicon_bp"]
        > normalized_evaluation["max_amplicon_bp"]
    ):
        raise SkillError(
            f"{context}.evaluation.min_amplicon_bp must be <= max_amplicon_bp"
        )
    if (
        normalized_evaluation["specificity_target_genome_id"] is None
        and (
            normalized_evaluation["specificity_catalog_path"] is not None
            or normalized_evaluation["specificity_cache_dir"] is not None
        )
    ):
        raise SkillError(
            f"{context}.evaluation specificity catalog/cache requires "
            "specificity_target_genome_id"
        )

    records = value.get("records")
    if not isinstance(records, list) or not records:
        raise SkillError(f"{context}.records must be a non-empty array")
    normalized_records = [
        _normalise_external_primer_handoff_record(record, index)
        for index, record in enumerate(records)
    ]
    record_ids = [record["record_id"] for record in normalized_records]
    if len(record_ids) != len(set(record_ids)):
        raise SkillError(f"{context}.records contains duplicate record_id values")
    normalized_records.sort(key=lambda record: record["record_id"])
    return {
        "schema": EXTERNAL_PRIMER_HANDOFF_REQUEST_SCHEMA,
        "collection_id": collection_id,
        "target": normalized_target,
        "evaluation": normalized_evaluation,
        "records": normalized_records,
    }


def _catalog_entry_path(script_path: Path) -> Path:
    return script_path.resolve().parent / "catalog_entry.json"


def _intents_descriptor_path(script_path: Path) -> Path:
    return script_path.resolve().parent / "INTENTS.json"


def _read_catalog_entry(script_path: Path) -> dict[str, Any]:
    path = _catalog_entry_path(script_path)
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (FileNotFoundError, json.JSONDecodeError):
        return {}
    return payload if isinstance(payload, dict) else {}


def _skill_info_payload(script_path: Path) -> dict[str, Any]:
    catalog_entry = _read_catalog_entry(script_path)
    catalog_path = _catalog_entry_path(script_path)
    return {
        "schema": SKILL_INFO_SCHEMA,
        "name": str(catalog_entry.get("name") or SKILL_NAME),
        "version": str(catalog_entry.get("version") or "unknown"),
        "status": str(catalog_entry.get("status") or "unknown"),
        "request_schema": REQUEST_SCHEMA,
        "result_schema": RESULT_SCHEMA,
        "execution_manifest_schema": EXECUTION_MANIFEST_SCHEMA,
        "execution_proposal_schema": EXECUTION_PROPOSAL_SCHEMA,
        "execution_approval_schema": EXECUTION_APPROVAL_SCHEMA,
        "approved_execution_request_schema": APPROVED_EXECUTION_REQUEST_SCHEMA,
        "supported_request_modes": list(SUPPORTED_REQUEST_MODES),
        "intents_runtime_schema": INTENTS_RUNTIME_SCHEMA,
        "intents_runtime_request_mode": "intents",
        "has_demo": bool(catalog_entry.get("has_demo", True)),
        "demo_command": str(
            catalog_entry.get("demo_command")
            or "python clawbio.py run gentle-cloning --demo"
        ),
        "catalog_entry_path": str(catalog_path),
        "catalog_entry_loaded": bool(catalog_entry),
        "runtime_version_command": "gentle_cli --version",
        "runtime_version_request_mode": "version",
        "runtime_lineage": LOCAL_RUNTIME_LINEAGE,
        "version_scope": VERSION_SCOPE,
        "classical_gentle_disambiguation": CLASSICAL_GENTLE_DISAMBIGUATION,
        "external_tool_resources": EXTERNAL_TOOL_RESOURCES,
        "ui_intent_support": {
            "catalog_request_mode": "capabilities",
            "catalog_shell_line": UI_INTENT_DISCOVERY_SHELL_LINE,
            "catalog_result_field": "ui_intent_catalog",
            "catalog_error_field": "ui_intent_catalog_error",
            "suggested_action_kind": "ui_intent",
            "suggested_action_ui_intent_field": "ui_intent",
        },
    }


def _skill_info_chat_summary_lines(info: dict[str, Any]) -> list[str]:
    name = str(info.get("name") or SKILL_NAME)
    version = str(info.get("version") or "unknown")
    status = str(info.get("status") or "unknown")
    return [
        f"{name} skill version {version} ({status}).",
        f"Request schema: {info.get('request_schema')}; result schema: {info.get('result_schema')}.",
        "Use request mode `version` when you need the installed local GENtle rewrite runtime version.",
        "Use request mode `intents` to compare the installed wrapper's live intent descriptor with ClawBio's on-disk snapshot.",
        "Use `resources status` or `services status` to check RNAfold/ViennaRNA and rnapkin executable-resource readiness.",
        CLASSICAL_GENTLE_DISAMBIGUATION,
    ]


def _read_intents_descriptor(script_path: Path) -> tuple[Path, str, dict[str, Any]]:
    path = _intents_descriptor_path(script_path)
    try:
        descriptor_text = path.read_text(encoding="utf-8")
    except FileNotFoundError as e:
        raise SkillError(f"INTENTS.json descriptor was not found at '{path}'") from e
    try:
        descriptor = json.loads(descriptor_text)
    except json.JSONDecodeError as e:
        raise SkillError(f"invalid JSON in INTENTS.json descriptor '{path}': {e}") from e
    if not isinstance(descriptor, dict):
        raise SkillError(f"INTENTS.json descriptor '{path}' must contain an object")
    if descriptor.get("schema") != "clawbio.skill_intents.v1":
        raise SkillError(
            "INTENTS.json descriptor has unexpected schema "
            f"`{descriptor.get('schema')}`"
        )
    return path, descriptor_text, descriptor


def _runtime_route_request_modes(
    route: dict[str, Any],
    skill_root: Path,
) -> tuple[list[str], list[str]]:
    request_modes: set[str] = set()
    warnings: list[str] = []
    for index, step in enumerate(route.get("plan") or []):
        if not isinstance(step, dict):
            warnings.append(f"route step {index} is not an object")
            continue
        if step.get("demo") is True:
            request_modes.add("demo")
        input_template = step.get("input_template")
        if isinstance(input_template, dict) and isinstance(input_template.get("mode"), str):
            request_modes.add(input_template["mode"])
        input_path = step.get("input")
        if isinstance(input_path, str) and input_path.strip():
            resolved = skill_root / input_path
            try:
                payload = json.loads(resolved.read_text(encoding="utf-8"))
            except (OSError, json.JSONDecodeError) as e:
                warnings.append(f"could not inspect {input_path}: {e}")
                continue
            mode = payload.get("mode") if isinstance(payload, dict) else None
            if isinstance(mode, str) and mode.strip():
                request_modes.add(mode)
    return sorted(request_modes), warnings


def _runtime_route_example_filenames(route: dict[str, Any]) -> list[str]:
    filenames: list[str] = []
    for step in route.get("plan") or []:
        if not isinstance(step, dict):
            continue
        input_path = step.get("input")
        if isinstance(input_path, str) and input_path.strip():
            filenames.append(Path(input_path).name)
    return filenames


def _runtime_route_requires_slots(route: dict[str, Any]) -> bool:
    for step in route.get("plan") or []:
        if isinstance(step, dict) and isinstance(step.get("slots"), dict):
            if any(
                isinstance(slot, dict) and bool(slot.get("required"))
                for slot in step["slots"].values()
            ):
                return True
    return False


def _intents_runtime_payload(script_path: Path) -> dict[str, Any]:
    descriptor_path, descriptor_text, descriptor = _read_intents_descriptor(script_path)
    skill_root = descriptor_path.parent
    routes: list[dict[str, Any]] = []
    warnings: list[str] = []
    for route in descriptor.get("routes") or []:
        if not isinstance(route, dict):
            warnings.append("skipped non-object route in INTENTS.json")
            continue
        request_modes, route_warnings = _runtime_route_request_modes(route, skill_root)
        warnings.extend(
            f"{route.get('intent_id', '(unknown)')}: {warning}"
            for warning in route_warnings
        )
        routes.append(
            {
                "intent_id": route.get("intent_id"),
                "description": route.get("description", ""),
                "trigger_terms": route.get("trigger_terms", []),
                "demo_policy": route.get("demo_policy", "never_unless_explicit"),
                "request_modes": request_modes,
                "example_filenames": _runtime_route_example_filenames(route),
                "has_input_template": any(
                    isinstance(step, dict) and isinstance(step.get("input_template"), dict)
                    for step in (route.get("plan") or [])
                ),
                "requires_slots": _runtime_route_requires_slots(route),
                "plan_step_count": len(route.get("plan") or []),
            }
        )
    return {
        "schema": INTENTS_RUNTIME_SCHEMA,
        "skill": descriptor.get("skill") or SKILL_NAME,
        "intents_schema": descriptor.get("schema"),
        "descriptor_path": str(descriptor_path),
        "descriptor_sha256": hashlib.sha256(
            descriptor_text.encode("utf-8")
        ).hexdigest(),
        "supported_request_modes": list(SUPPORTED_REQUEST_MODES),
        "route_count": len(routes),
        "routes": routes,
        "warnings": warnings,
    }


def _intents_runtime_chat_summary_lines(payload: dict[str, Any]) -> list[str]:
    return [
        f"GENtle intent descriptor runtime schema: {payload.get('schema')}.",
        f"Descriptor SHA-256: {payload.get('descriptor_sha256')}.",
        f"Routes available: {payload.get('route_count', 0)}.",
        "ClawBio can compare this hash with its local INTENTS.json snapshot before refreshing route metadata.",
    ]


def _request_mode_warnings(request: Request | None) -> list[str]:
    if request is None:
        return []
    equivalent = DEPRECATED_REQUEST_MODE_EQUIVALENTS.get(request.mode)
    if not equivalent:
        return []
    return [
        (
            f"Request mode `{request.mode}` is accepted for backward compatibility "
            f"but is deprecated; prefer {equivalent}."
        )
    ]


def _repo_root_candidates(script_path: Path) -> list[Path]:
    candidates: list[Path] = [Path.cwd()]
    candidates.extend(script_path.resolve().parents)
    unique: list[Path] = []
    seen: set[str] = set()
    for candidate in candidates:
        key = str(candidate)
        if key not in seen:
            seen.add(key)
            unique.append(candidate)
    return unique


def _find_repo_root(script_path: Path) -> Path | None:
    for candidate in _repo_root_candidates(script_path):
        if (candidate / "Cargo.toml").exists() and (
            candidate / "src" / "bin" / "gentle_cli.rs"
        ).exists():
            return candidate
    return None


def _find_repo_root_from_path(path: Path) -> Path | None:
    candidate = path.resolve()
    if candidate.is_file():
        candidate = candidate.parent
    while True:
        if (candidate / "Cargo.toml").exists() and (
            candidate / "src" / "bin" / "gentle_cli.rs"
        ).exists():
            return candidate
        if candidate.parent == candidate:
            return None
        candidate = candidate.parent


def _request_path_candidates(raw_path: str, script_path: Path) -> list[Path]:
    path = Path(raw_path)
    if path.is_absolute():
        return [path]

    candidates: list[Path] = [Path.cwd() / path]
    parts = path.parts
    if len(parts) >= 2 and parts[0] == "skills" and parts[1] == SKILL_NAME:
        candidates.append(script_path.resolve().parent / Path(*parts[2:]))
    elif parts and parts[0] == SKILL_NAME:
        candidates.append(script_path.resolve().parent / Path(*parts[1:]))
    for base in script_path.resolve().parents:
        candidates.append(base / path)
    configured_repo_root = os.environ.get("GENTLE_REPO_ROOT", "").strip()
    if configured_repo_root:
        candidates.append(Path(configured_repo_root) / path)
    repo_root = _find_repo_root(script_path)
    if repo_root is not None:
        candidates.append(repo_root / path)

    unique: list[Path] = []
    seen: set[str] = set()
    for candidate in candidates:
        key = str(candidate)
        if key not in seen:
            seen.add(key)
            unique.append(candidate)
    return unique


def _resolve_request_path(raw_path: str, script_path: Path) -> str:
    trimmed = raw_path.strip()
    if not trimmed:
        return raw_path
    for candidate in _request_path_candidates(trimmed, script_path):
        if candidate.exists():
            return str(candidate.resolve())
    return raw_path


def _resolve_existing_request_file(raw_path: str, script_path: Path) -> Path:
    trimmed = raw_path.strip()
    candidates = _request_path_candidates(trimmed, script_path)
    for candidate in candidates:
        if candidate.exists():
            return candidate.resolve()
    tried = ", ".join(f"`{candidate}`" for candidate in candidates)
    raise SkillError(
        f"request JSON file '{raw_path}' was not found. "
        f"Current cwd: `{Path.cwd()}`. Tried: {tried}"
    )


def _resolve_cli(explicit: str | None, script_path: Path) -> CliResolution:
    if explicit:
        return CliResolution(
            argv_prefix=shlex.split(explicit),
            cwd=None,
            label=f"explicit --gentle-cli: {explicit}",
        )
    env_cmd = os.environ.get("GENTLE_CLI_CMD", "").strip()
    if env_cmd:
        return CliResolution(
            argv_prefix=shlex.split(env_cmd),
            cwd=None,
            label=f"GENTLE_CLI_CMD: {env_cmd}",
        )
    path_hit = shutil.which("gentle_cli")
    if path_hit:
        return CliResolution(
            argv_prefix=[path_hit],
            cwd=None,
            label=f"PATH gentle_cli: {path_hit}",
        )
    repo_root = _find_repo_root(script_path)
    if repo_root is not None and shutil.which("cargo"):
        return CliResolution(
            argv_prefix=["cargo", "run", "--quiet", "--bin", "gentle_cli", "--"],
            cwd=str(repo_root),
            label=f"cargo run fallback in repo: {repo_root}",
        )
    raise SkillError(
        "Could not resolve gentle_cli executable. "
        "Recommended: set GENTLE_CLI_CMD to the included "
        "'./gentle_local_checkout_cli.sh' launcher (typically with "
        "GENTLE_REPO_ROOT=/absolute/path/to/GENtle), or to a Docker/OCI or "
        "Apptainer/Singularity command such as "
        "'docker run --rm -i -v \"$PWD\":/work -w /work "
        "ghcr.io/smoe/gentle_rs:cli' or "
        "'./gentle_apptainer_cli.sh /absolute/path/to/gentle.sif'. "
        "Alternatives: use --gentle-cli, install gentle_cli on PATH, or run "
        "from a local repo with cargo available."
    )


def _local_source_runtime_repo_root(
    resolution: CliResolution, script_path: Path
) -> Path | None:
    tokens = resolution.argv_prefix
    if not tokens:
        return None
    executable_name = Path(tokens[0]).name
    uses_cargo_run = executable_name in {"cargo", "cargo.exe"} and "run" in tokens[1:]
    launcher_token = next(
        (
            token
            for token in tokens
            if Path(token).name == "gentle_local_checkout_cli.sh"
        ),
        None,
    )
    uses_local_launcher = launcher_token is not None
    if not uses_cargo_run and not uses_local_launcher:
        return None
    if "--manifest-path" in tokens:
        index = tokens.index("--manifest-path")
        if index + 1 < len(tokens):
            manifest = Path(tokens[index + 1]).expanduser()
            if not manifest.is_absolute():
                manifest = Path(resolution.cwd or Path.cwd()) / manifest
            return manifest.resolve().parent
    configured = os.environ.get("GENTLE_REPO_ROOT", "").strip()
    if configured:
        return Path(configured).expanduser().resolve()
    if resolution.cwd:
        repo_root = _find_repo_root_from_path(Path(resolution.cwd))
        if repo_root is not None:
            return repo_root
    executable = Path(launcher_token or tokens[0]).expanduser()
    if executable.exists():
        repo_root = _find_repo_root_from_path(executable)
        if repo_root is not None:
            return repo_root
    return _find_repo_root(script_path)


def _pin_delegated_runtime(
    resolution: CliResolution,
    execution_cwd: Path,
    script_path: Path,
    *,
    require_local_binary: bool,
) -> CliResolution:
    tokens = resolution.argv_prefix
    if not tokens:
        raise SkillError("delegated GENtle runtime command is empty")
    executable_name = Path(tokens[0]).name
    if executable_name in {"docker", "podman", "nerdctl"} and "run" in tokens[1:]:
        immutable_image = any(
            re.search(r"@sha256:[0-9a-fA-F]{64}$", token)
            or re.fullmatch(r"sha256:[0-9a-fA-F]{64}", token)
            for token in tokens[1:]
        )
        if not immutable_image:
            raise SkillError(
                "approval-gated OCI execution requires an immutable image digest "
                "such as ghcr.io/owner/image@sha256:..."
            )
    repo_root = _local_source_runtime_repo_root(resolution, script_path)
    if repo_root is None:
        return resolution
    uses_local_launcher = any(
        Path(token).name == "gentle_local_checkout_cli.sh" for token in tokens
    )
    if uses_local_launcher:
        os.environ.setdefault(
            "GENTLE_REFERENCE_CACHE_DIR", str((repo_root / "data" / "genomes").resolve())
        )
        os.environ.setdefault(
            "GENTLE_HELPER_CACHE_DIR",
            str((repo_root / "data" / "helper_genomes").resolve()),
        )
    configured_target = os.environ.get("CARGO_TARGET_DIR", "").strip()
    if configured_target:
        target_dir = Path(configured_target).expanduser()
        if not target_dir.is_absolute():
            target_dir = execution_cwd / target_dir
    else:
        target_dir = repo_root / "target"
    executable = target_dir / "debug" / (
        "gentle_cli.exe" if os.name == "nt" else "gentle_cli"
    )
    if not executable.is_file():
        if not require_local_binary:
            return resolution
        raise SkillError(
            "the local Cargo launcher completed without producing the expected "
            f"gentle_cli executable at '{executable.resolve()}'"
        )
    return CliResolution(
        argv_prefix=[str(executable.resolve())],
        cwd=resolution.cwd,
        label=f"content-pinned executable produced by {resolution.label}",
    )


def _approval_environment_binding(execution_cwd: Path) -> dict[str, Any]:
    binding: dict[str, Any] = {}
    for name in APPROVAL_ENVIRONMENT_VARIABLES:
        value = os.environ.get(name)
        if value is None:
            continue
        record: dict[str, Any] = {"value": value}
        if name.endswith("_BIN") or name.endswith("_TOOL"):
            candidate = Path(value).expanduser()
            if not candidate.is_absolute():
                relative = execution_cwd / candidate
                if relative.is_file():
                    candidate = relative
                else:
                    located = shutil.which(value)
                    if located:
                        candidate = Path(located)
            if candidate.is_file():
                resolved = candidate.resolve()
                record.update(
                    {
                        "resolved_path": str(resolved),
                        "size_bytes": resolved.stat().st_size,
                        "sha256": _sha256_file_prefixed(resolved),
                    }
                )
        binding[name] = record
    return binding


def _wrapper_contract_binding(script_path: Path) -> dict[str, Any]:
    candidates = [
        (script_path.resolve(), "runner"),
        (_catalog_entry_path(script_path).resolve(), "catalog_entry"),
        (_intents_descriptor_path(script_path).resolve(), "intent_descriptor"),
    ]
    files = [
        _content_artifact_record(path, "wrapper_contract", script_path.parent)
        | {"contract_part": contract_part}
        for path, contract_part in candidates
        if path.is_file()
    ]
    missing_files = [
        {"contract_part": contract_part, "path": str(path)}
        for path, contract_part in candidates
        if not path.is_file()
    ]
    return {
        "skill": SKILL_NAME,
        "skill_version": SKILL_CONTRACT_VERSION,
        "contract_files_complete": not missing_files,
        "files": files,
        "missing_files": missing_files,
    }


def _safe_id_component(value: str, fallback: str) -> str:
    safe = re.sub(r"[^A-Za-z0-9]+", "_", value.strip().lower()).strip("_")
    return safe or fallback


def _gene_protein_2d_gel_svg_path(request: Request) -> str:
    if request.expected_artifacts:
        return request.expected_artifacts[0]
    gene_id = _safe_id_component(request.gene_symbol or "", "gene")
    return f"exports/{gene_id}_ensembl_protein_2d_gel.svg"


def _normalise_gene_protein_2d_gel_request(request: Request) -> None:
    if not isinstance(request.gene_symbol, str) or not request.gene_symbol.strip():
        raise SkillError(
            "mode=gene-protein-2d-gel requires non-empty string field 'gene_symbol'"
        )
    request.gene_symbol = request.gene_symbol.strip()

    if request.species is None:
        request.species = "homo_sapiens"
    if not isinstance(request.species, str) or not request.species.strip():
        raise SkillError("species must be a non-empty string when present")
    request.species = request.species.strip()

    if request.source is None:
        request.source = "ensembl"
    if not isinstance(request.source, str) or not request.source.strip():
        raise SkillError("source must be a non-empty string when present")
    request.source = request.source.strip().lower()
    if request.source != "ensembl":
        raise SkillError("mode=gene-protein-2d-gel currently supports source='ensembl'")

    if request.ladders is None:
        request.ladders = ["Protein Ladder 10-100 kDa"]
    if not isinstance(request.ladders, list) or not request.ladders:
        raise SkillError("ladders must be a non-empty string array when present")
    if not all(isinstance(v, str) and v.strip() for v in request.ladders):
        raise SkillError("ladders must contain non-empty strings")
    request.ladders = [v.strip() for v in request.ladders]

    if request.expected_artifacts is None:
        request.expected_artifacts = [_gene_protein_2d_gel_svg_path(request)]
    if len(request.expected_artifacts) != 1:
        raise SkillError(
            "mode=gene-protein-2d-gel expects exactly one SVG expected_artifacts path"
        )
    if not request.expected_artifacts[0].lower().endswith(".svg"):
        raise SkillError(
            "mode=gene-protein-2d-gel expected_artifacts path must end with .svg"
        )


def _normalise_exon_skip_interval_list(value: Any, field_name: str) -> list[dict[str, int]]:
    if value is None:
        return []
    if not isinstance(value, list):
        raise SkillError(f"{field_name} must be an array when present")
    intervals: list[dict[str, int]] = []
    for idx, item in enumerate(value):
        if not isinstance(item, dict):
            raise SkillError(f"{field_name}[{idx}] must be an object")
        try:
            start = int(item.get("start_1based"))
            end = int(item.get("end_1based"))
        except (TypeError, ValueError) as e:
            raise SkillError(
                f"{field_name}[{idx}] requires integer start_1based/end_1based"
            ) from e
        if start <= 0 or end < start:
            raise SkillError(
                f"{field_name}[{idx}] must satisfy 1 <= start_1based <= end_1based"
            )
        intervals.append({"start_1based": start, "end_1based": end})
    return intervals


def _normalise_optional_string_array(value: Any, field_name: str) -> list[str]:
    if value is None:
        return []
    if not isinstance(value, list):
        raise SkillError(f"{field_name} must be a string array when present")
    result: list[str] = []
    for idx, item in enumerate(value):
        if not isinstance(item, str) or not item.strip():
            raise SkillError(f"{field_name}[{idx}] must be a non-empty string")
        result.append(item.strip())
    return result


def _normalise_length_mod3_values(value: Any, field_name: str) -> list[int]:
    if value is None:
        return []
    if not isinstance(value, list):
        raise SkillError(f"{field_name} must be an integer array when present")
    result: list[int] = []
    for idx, item in enumerate(value):
        try:
            parsed = int(item)
        except (TypeError, ValueError) as e:
            raise SkillError(f"{field_name}[{idx}] must be an integer") from e
        if parsed not in {0, 1, 2}:
            raise SkillError(f"{field_name}[{idx}] must be one of 0, 1, or 2")
        result.append(parsed)
    return sorted(set(result))


def _normalise_exon_skip_request(request: Request) -> None:
    if request.mode == "exon-skip-plan":
        if not isinstance(request.seq_id, str) or not request.seq_id.strip():
            raise SkillError("mode=exon-skip-plan requires non-empty string field 'seq_id'")
        request.seq_id = request.seq_id.strip()
        try:
            transcript_feature_id = int(request.transcript_feature_id)
        except (TypeError, ValueError) as e:
            raise SkillError(
                "mode=exon-skip-plan requires integer field 'transcript_feature_id'"
            ) from e
        if transcript_feature_id < 0:
            raise SkillError("transcript_feature_id must be >= 0")
        request.transcript_feature_id = transcript_feature_id
        request.skip_candidate_ids = _normalise_optional_string_array(
            request.skip_candidate_ids, "skip_candidate_ids"
        )
        request.skip_intervals_1based = _normalise_exon_skip_interval_list(
            request.skip_intervals_1based, "skip_intervals_1based"
        )
        request.overlap_intervals_1based = _normalise_exon_skip_interval_list(
            request.overlap_intervals_1based, "overlap_intervals_1based"
        )
        request.length_mod3_values = _normalise_length_mod3_values(
            request.length_mod3_values, "length_mod3_values"
        )
        request.coding_mod3_values = _normalise_length_mod3_values(
            request.coding_mod3_values, "coding_mod3_values"
        )
        request.coding_contexts = _normalise_optional_string_array(
            request.coding_contexts, "coding_contexts"
        )
        request.cds_phase_entry_kinds = _normalise_optional_string_array(
            request.cds_phase_entry_kinds, "cds_phase_entry_kinds"
        )
        if request.plan_id is not None:
            if not isinstance(request.plan_id, str) or not request.plan_id.strip():
                raise SkillError("plan_id must be a non-empty string when present")
            request.plan_id = request.plan_id.strip()
        if request.feature_query is not None and not isinstance(
            request.feature_query, (dict, str)
        ):
            raise SkillError("feature_query must be an object or JSON string when present")
    elif request.mode == "exon-skip-materialize":
        if not isinstance(request.plan_id, str) or not request.plan_id.strip():
            raise SkillError(
                "mode=exon-skip-materialize requires non-empty string field 'plan_id'"
            )
        request.plan_id = request.plan_id.strip()
        if request.confirm is not True:
            raise SkillError(
                "mode=exon-skip-materialize requires confirm=true because it creates derived sequences"
            )
        request.candidate_ids = _normalise_optional_string_array(
            request.candidate_ids, "candidate_ids"
        )
        request.return_items = _normalise_optional_string_array(
            request.return_items, "return_items"
        )
        if request.output_prefix is not None:
            if not isinstance(request.output_prefix, str) or not request.output_prefix.strip():
                raise SkillError("output_prefix must be a non-empty string when present")
            request.output_prefix = request.output_prefix.strip()


def _normalise_protein_residue_genomic_coordinate_request(request: Request) -> None:
    if not isinstance(request.seq_id, str) or not request.seq_id.strip():
        raise SkillError(
            "mode=protein-residue-genomic-coordinates requires non-empty string field 'seq_id'"
        )
    request.seq_id = request.seq_id.strip()

    if request.transcript_id is not None:
        if not isinstance(request.transcript_id, str):
            raise SkillError("transcript_id must be a non-empty string when present")
        request.transcript_id = request.transcript_id.strip() or None

    def _coerce_optional_positive_int(value: Any, field_name: str) -> int | None:
        if value is None:
            return None
        if isinstance(value, str) and not value.strip():
            return None
        try:
            coerced = int(value)
        except (TypeError, ValueError) as e:
            raise SkillError(f"{field_name} must be an integer when present") from e
        if coerced <= 0:
            raise SkillError(f"{field_name} must be >= 1")
        return coerced

    residue_start = _coerce_optional_positive_int(
        request.residue_start_1based, "residue_start_1based"
    )
    if residue_start is None:
        raise SkillError(
            "mode=protein-residue-genomic-coordinates requires integer field 'residue_start_1based'"
        )
    residue_end = _coerce_optional_positive_int(
        request.residue_end_1based, "residue_end_1based"
    )
    if residue_end is None:
        residue_end = residue_start
    if residue_end < residue_start:
        raise SkillError("residue_end_1based must be >= residue_start_1based")
    request.residue_start_1based = residue_start
    request.residue_end_1based = residue_end


def _normalise_optional_str(value: Any, field_name: str) -> str | None:
    if value is None:
        return None
    if not isinstance(value, str):
        raise SkillError(f"{field_name} must be a non-empty string when present")
    trimmed = value.strip()
    if not trimmed:
        return None
    return trimmed


def _require_str(value: Any, field_name: str, mode: str) -> str:
    normalized = _normalise_optional_str(value, field_name)
    if normalized is None:
        raise SkillError(f"mode={mode} requires non-empty string field '{field_name}'")
    return normalized


def _coerce_optional_int(value: Any, field_name: str, *, minimum: int = 0) -> int | None:
    if value is None:
        return None
    if isinstance(value, str) and not value.strip():
        return None
    try:
        coerced = int(value)
    except (TypeError, ValueError) as e:
        raise SkillError(f"{field_name} must be an integer when present") from e
    if coerced < minimum:
        raise SkillError(f"{field_name} must be >= {minimum}")
    return coerced


def _normalise_construct_reasoning_inspection_request(request: Request) -> None:
    request.graph_id = _require_str(request.graph_id, "graph_id", request.mode)
    if request.mode == "construct-reasoning-list-inspections":
        request.fact_id = _normalise_optional_str(request.fact_id, "fact_id")
        request.annotation_id = _normalise_optional_str(
            request.annotation_id, "annotation_id"
        )
        request.candidate_id = _normalise_optional_str(
            request.candidate_id, "candidate_id"
        )
        request.evidence_id = _normalise_optional_str(request.evidence_id, "evidence_id")
        request.seq_id = _normalise_optional_str(request.seq_id, "seq_id")
        request.action_kind = _normalise_optional_str(request.action_kind, "action_kind")
        request.summary_id = _normalise_optional_str(request.summary_id, "summary_id")
        return

    request.action_id = _require_str(request.action_id, "action_id", request.mode)
    request.word_size = _coerce_optional_int(request.word_size, "word_size", minimum=1)
    request.step_bp = _coerce_optional_int(request.step_bp, "step_bp", minimum=1)
    request.max_mismatches = _coerce_optional_int(
        request.max_mismatches, "max_mismatches", minimum=0
    )
    request.tile_bp = _coerce_optional_int(request.tile_bp, "tile_bp", minimum=1)
    request.dotplot_id = _normalise_optional_str(request.dotplot_id, "dotplot_id")
    request.render_svg_path = _normalise_optional_str(
        request.render_svg_path, "render_svg_path"
    )
    if request.render_svg_path is None:
        request.render_svg_path = _normalise_optional_str(request.svg_path, "svg_path")


def _normalise_primer_backend_options(request: Request) -> None:
    request.backend = _normalise_optional_str(request.backend, "backend")
    if request.backend is not None and request.backend not in {"auto", "internal", "primer3"}:
        raise SkillError("backend must be one of auto|internal|primer3 when present")
    request.primer3_executable = _normalise_optional_str(
        request.primer3_executable, "primer3_executable"
    )


def _normalise_seq_feature_request(request: Request) -> None:
    request.seq_id = _require_str(request.seq_id, "seq_id", request.mode)
    source_feature_id = _coerce_optional_int(
        request.source_feature_id, "source_feature_id", minimum=0
    )
    if source_feature_id is None:
        raise SkillError(
            f"mode={request.mode} requires integer field 'source_feature_id'"
        )
    request.source_feature_id = source_feature_id


def _normalise_primer_design_request(request: Request) -> None:
    if request.request_json is None and request.operation is None:
        raise SkillError(
            f"mode={request.mode} requires 'request_json' or 'operation'"
        )
    _normalise_primer_backend_options(request)


def _normalise_primer_report_request(request: Request) -> None:
    if request.mode.endswith("-show") or request.mode.endswith("-export"):
        request.report_id = _require_str(request.report_id, "report_id", request.mode)
    if request.mode.endswith("-export"):
        request.output_path = _require_str(request.output_path, "output_path", request.mode)


def _normalise_qpcr_seed_request(request: Request) -> None:
    _normalise_seq_feature_request(request)
    request.transcript_id = _normalise_optional_str(request.transcript_id, "transcript_id")
    request.qpcr_mode = _normalise_optional_str(request.qpcr_mode, "qpcr_mode")
    if request.qpcr_mode is not None and request.qpcr_mode not in {
        "shared_gene",
        "distinguish_transcript",
    }:
        raise SkillError("qpcr_mode must be shared_gene or distinguish_transcript")
    if request.qpcr_mode == "distinguish_transcript" and request.transcript_id is None:
        raise SkillError(
            "mode=qpcr-seed-from-splicing with qpcr_mode=distinguish_transcript requires transcript_id"
        )
    request.specificity_evidence = _normalise_optional_str(
        request.specificity_evidence, "specificity_evidence"
    )
    if request.specificity_evidence is not None and request.specificity_evidence not in {
        "junction_only",
        "unique_exon_or_chain",
        "either_prefer_junction",
    }:
        raise SkillError(
            "specificity_evidence must be junction_only|unique_exon_or_chain|either_prefer_junction"
        )
    if request.specificity_evidence is not None and request.qpcr_mode != "distinguish_transcript":
        raise SkillError(
            "specificity_evidence is only valid with qpcr_mode=distinguish_transcript"
        )


def _normalise_cdna_assay_test_request(request: Request) -> None:
    _normalise_seq_feature_request(request)
    request.forward_primer = _require_str(request.forward_primer, "forward_primer", request.mode)
    request.reverse_primer = _require_str(request.reverse_primer, "reverse_primer", request.mode)
    if request.mode == "cdna-qpcr-test":
        request.probe = _require_str(request.probe, "probe", request.mode)
    request.transcript_id = _normalise_optional_str(request.transcript_id, "transcript_id")
    request.transcript_order = _normalise_optional_str(
        request.transcript_order, "transcript_order"
    )
    if request.transcript_order is not None and request.transcript_order not in {
        "transcript_id",
        "genomic_first_exon",
        "genomic_last_exon",
        "antisense_first_exon",
    }:
        raise SkillError(
            "transcript_order must be transcript_id|genomic_first_exon|genomic_last_exon|antisense_first_exon"
        )
    request.map_coordinate_mode = _normalise_optional_str(
        request.map_coordinate_mode, "map_coordinate_mode"
    )
    if request.map_coordinate_mode is not None and request.map_coordinate_mode not in {
        "cdna",
        "genomic_aligned",
    }:
        raise SkillError("map_coordinate_mode must be cdna or genomic_aligned")
    request.output_path = _normalise_optional_str(request.output_path, "output_path")
    request.svg_path = _normalise_optional_str(request.svg_path, "svg_path")
    request.product_output_prefix = _normalise_optional_str(
        request.product_output_prefix, "product_output_prefix"
    )
    request.product_gel_svg_path = _normalise_optional_str(
        request.product_gel_svg_path, "product_gel_svg_path"
    )
    if request.materialize_products is not None and not isinstance(
        request.materialize_products, bool
    ):
        raise SkillError("materialize_products must be true or false")
    if request.product_gel_ladders is not None:
        if not isinstance(request.product_gel_ladders, list) or not all(
            isinstance(item, str) and item.strip() for item in request.product_gel_ladders
        ):
            raise SkillError("product_gel_ladders must be a list of non-empty strings")
        request.product_gel_ladders = [item.strip() for item in request.product_gel_ladders]
    request.min_amplicon_bp = _coerce_optional_int(
        request.min_amplicon_bp, "min_amplicon_bp", minimum=1
    )
    request.max_amplicon_bp = _coerce_optional_int(
        request.max_amplicon_bp, "max_amplicon_bp", minimum=1
    )
    if (
        request.min_amplicon_bp is not None
        and request.max_amplicon_bp is not None
        and request.max_amplicon_bp < request.min_amplicon_bp
    ):
        raise SkillError("max_amplicon_bp must be >= min_amplicon_bp")
    request.max_mismatches = _coerce_optional_int(
        request.max_mismatches, "max_mismatches", minimum=0
    )
    request.require_3prime_exact_bases = _coerce_optional_int(
        request.require_3prime_exact_bases,
        "require_3prime_exact_bases",
        minimum=0,
    )
    if request.expected_artifacts is None:
        expected = [
            path
            for path in (request.output_path, request.svg_path, request.product_gel_svg_path)
            if isinstance(path, str) and path.strip()
        ]
        if expected:
            request.expected_artifacts = expected


def _normalise_restriction_cloning_request(request: Request) -> None:
    if request.mode == "restriction-cloning-pcr-handoff":
        if request.request_json is None and request.operation is None:
            raise SkillError(
                "mode=restriction-cloning-pcr-handoff requires 'request_json' or 'operation'"
            )
        return
    if request.mode == "restriction-cloning-pcr-handoff-seed":
        request.report_id = _require_str(request.report_id, "report_id", request.mode)
        request.vector_seq_id = _require_str(request.vector_seq_id, "vector_seq_id", request.mode)
        request.pair_rank = _coerce_optional_int(request.pair_rank, "pair_rank", minimum=1)
        request.handoff_mode = _normalise_optional_str(request.handoff_mode, "handoff_mode")
        if request.handoff_mode is not None and request.handoff_mode not in {
            "single_site",
            "directed_pair",
        }:
            raise SkillError("handoff_mode must be single_site or directed_pair")
        for field_name in (
            "forward_enzyme",
            "reverse_enzyme",
            "forward_leader_5prime",
            "reverse_leader_5prime",
        ):
            setattr(request, field_name, _normalise_optional_str(getattr(request, field_name), field_name))
        return
    if request.mode == "restriction-cloning-vector-suggestions":
        request.seq_id = _require_str(request.seq_id, "seq_id", request.mode)
        return
    if request.mode in {
        "restriction-cloning-handoff-show",
        "restriction-cloning-handoff-export",
    }:
        request.report_id = _require_str(request.report_id, "report_id", request.mode)
    if request.mode == "restriction-cloning-handoff-export":
        request.output_path = _require_str(request.output_path, "output_path", request.mode)


def _normalise_pcr_protocol_cartoon_request(request: Request) -> None:
    request.protocol_id = _require_str(request.protocol_id, "protocol_id", request.mode)
    request.output_path = _require_str(request.output_path, "output_path", request.mode)
    if request.expected_artifacts is None:
        request.expected_artifacts = [request.output_path]


def _normalise_transcript_qpcr_panel_request(request: Request) -> None:
    if not isinstance(request.seq_id, str) or not request.seq_id.strip():
        raise SkillError("mode=transcript-qpcr-panel requires non-empty string field 'seq_id'")
    request.seq_id = request.seq_id.strip()
    try:
        source_feature_id = int(request.source_feature_id)  # type: ignore[arg-type]
    except (TypeError, ValueError) as e:
        raise SkillError(
            "mode=transcript-qpcr-panel requires integer field 'source_feature_id'"
        ) from e
    if source_feature_id < 0:
        raise SkillError("source_feature_id must be a zero-based feature index >= 0")
    request.source_feature_id = source_feature_id
    if (
        not isinstance(request.shared_qpcr_report_id, str)
        or not request.shared_qpcr_report_id.strip()
    ):
        raise SkillError(
            "mode=transcript-qpcr-panel requires non-empty string field 'shared_qpcr_report_id'"
        )
    request.shared_qpcr_report_id = request.shared_qpcr_report_id.strip()
    if request.output_path is not None:
        if not isinstance(request.output_path, str) or not request.output_path.strip():
            raise SkillError("output_path must be a non-empty string when present")
        request.output_path = request.output_path.strip()


def _coerce_request(payload: dict[str, Any]) -> Request:
    schema = payload.get("schema")
    if schema != REQUEST_SCHEMA:
        raise SkillError(
            f"unsupported request schema '{schema}', expected '{REQUEST_SCHEMA}'"
        )
    mode = str(payload.get("mode", "")).strip()
    if not mode:
        raise SkillError("request missing required field 'mode'")
    timeout_secs = int(payload.get("timeout_secs", 180))
    if timeout_secs <= 0:
        raise SkillError("timeout_secs must be > 0")
    request = Request(
        mode=mode,
        timeout_secs=timeout_secs,
        state_path=payload.get("state_path"),
        raw_args=payload.get("raw_args"),
        shell_line=payload.get("shell_line"),
        operation=payload.get("operation"),
        request_json=payload.get("request_json"),
        workflow=payload.get("workflow"),
        workflow_path=payload.get("workflow_path"),
        system_id=payload.get("system_id"),
        prompt=payload.get("prompt"),
        catalog_path=payload.get("catalog_path"),
        base_url=payload.get("base_url"),
        model=payload.get("model"),
        connect_timeout_secs=payload.get("connect_timeout_secs"),
        read_timeout_secs=payload.get("read_timeout_secs"),
        max_retries=payload.get("max_retries"),
        max_response_bytes=payload.get("max_response_bytes"),
        include_state_summary=payload.get("include_state_summary"),
        max_candidates=payload.get("max_candidates"),
        allow_mutating_candidates=payload.get("allow_mutating_candidates"),
        plan=payload.get("plan"),
        plan_path=payload.get("plan_path"),
        plan_id=payload.get("plan_id"),
        graph_id=payload.get("graph_id"),
        fact_id=payload.get("fact_id"),
        annotation_id=payload.get("annotation_id"),
        candidate_id=payload.get("candidate_id"),
        candidate_ids=payload.get("candidate_ids"),
        evidence_id=payload.get("evidence_id"),
        summary_id=payload.get("summary_id"),
        action_kind=payload.get("action_kind"),
        action_id=payload.get("action_id"),
        confirm=payload.get("confirm"),
        transcript_feature_id=payload.get("transcript_feature_id"),
        skip_candidate_ids=payload.get("skip_candidate_ids"),
        skip_intervals_1based=payload.get("skip_intervals_1based"),
        overlap_intervals_1based=payload.get("overlap_intervals_1based"),
        length_mod3_values=payload.get("length_mod3_values"),
        coding_mod3_values=payload.get("coding_mod3_values"),
        coding_contexts=payload.get("coding_contexts"),
        cds_phase_entry_kinds=payload.get("cds_phase_entry_kinds"),
        feature_query=payload.get("feature_query"),
        return_items=payload.get("return_items", payload.get("return")),
        expected_artifacts=payload.get("expected_artifacts"),
        ensure_reference_prepared=payload.get("ensure_reference_prepared"),
        external_primer_handoff=payload.get("external_primer_handoff"),
        gene_symbol=payload.get("gene_symbol"),
        species=payload.get("species"),
        source=payload.get("source"),
        ladders=payload.get("ladders"),
        seq_id=payload.get("seq_id"),
        source_feature_id=payload.get("source_feature_id"),
        transcript_id=payload.get("transcript_id"),
        transcript_order=payload.get("transcript_order"),
        map_coordinate_mode=payload.get(
            "map_coordinate_mode", payload.get("transcript_map_coordinate_mode")
        ),
        qpcr_mode=payload.get("qpcr_mode"),
        specificity_evidence=payload.get("specificity_evidence"),
        backend=payload.get("backend"),
        primer3_executable=payload.get("primer3_executable"),
        report_id=payload.get("report_id"),
        forward_primer=payload.get("forward_primer"),
        reverse_primer=payload.get("reverse_primer"),
        probe=payload.get("probe"),
        min_amplicon_bp=payload.get("min_amplicon_bp"),
        max_amplicon_bp=payload.get("max_amplicon_bp"),
        max_mismatches=payload.get("max_mismatches"),
        word_size=payload.get("word_size"),
        step_bp=payload.get("step_bp"),
        tile_bp=payload.get("tile_bp"),
        dotplot_id=payload.get("dotplot_id"),
        render_svg_path=payload.get("render_svg_path"),
        require_3prime_exact_bases=payload.get("require_3prime_exact_bases"),
        vector_seq_id=payload.get("vector_seq_id"),
        pair_rank=payload.get("pair_rank"),
        handoff_mode=payload.get("handoff_mode"),
        forward_enzyme=payload.get("forward_enzyme"),
        reverse_enzyme=payload.get("reverse_enzyme"),
        forward_leader_5prime=payload.get("forward_leader_5prime"),
        reverse_leader_5prime=payload.get("reverse_leader_5prime"),
        protocol_id=payload.get("protocol_id"),
        shared_qpcr_report_id=payload.get("shared_qpcr_report_id"),
        output_prefix=payload.get("output_prefix"),
        output_path=payload.get("output_path"),
        svg_path=payload.get("svg_path"),
        materialize_products=payload.get("materialize_products"),
        product_output_prefix=payload.get("product_output_prefix"),
        product_gel_svg_path=payload.get("product_gel_svg_path"),
        product_gel_ladders=payload.get("product_gel_ladders"),
        residue_start_1based=payload.get("residue_start_1based"),
        residue_end_1based=payload.get("residue_end_1based"),
        claim_attribution_mode=payload.get("claim_attribution_mode"),
        presentation_profile=payload.get("presentation_profile"),
        input_claims=payload.get("input_claims"),
        delegation=payload.get("delegation"),
        input_bindings=payload.get("input_bindings"),
    )
    request.delegation = _normalise_delegation(request.delegation)
    request.input_bindings = _normalise_input_bindings(request.input_bindings)
    if request.claim_attribution_mode is not None:
        if request.claim_attribution_mode != STRICT_CLAIM_ATTRIBUTION_MODE:
            raise SkillError("claim_attribution_mode must be 'strict' when present")
    if request.presentation_profile is not None:
        if (
            not isinstance(request.presentation_profile, str)
            or not request.presentation_profile.strip()
        ):
            raise SkillError(
                "presentation_profile must be a non-empty string when present"
            )
        request.presentation_profile = request.presentation_profile.strip()
    if request.input_claims is not None:
        if not isinstance(request.input_claims, list):
            raise SkillError("input_claims must be a string array when present")
        if not all(
            isinstance(value, str) and value.strip()
            for value in request.input_claims
        ):
            raise SkillError("input_claims must contain non-empty strings")
        request.input_claims = [value.strip() for value in request.input_claims]
    if request.expected_artifacts is not None:
        if not isinstance(request.expected_artifacts, list):
            raise SkillError("expected_artifacts must be a string array when present")
        if not all(isinstance(v, str) and v.strip() for v in request.expected_artifacts):
            raise SkillError("expected_artifacts must contain non-empty strings")
    if request.mode not in SUPPORTED_REQUEST_MODES:
        raise SkillError(
            "unsupported mode. Use one of: " + ", ".join(SUPPORTED_REQUEST_MODES)
        )
    if request.mode == "raw":
        if not isinstance(request.raw_args, list) or not request.raw_args:
            raise SkillError("mode=raw requires non-empty string array 'raw_args'")
        if not all(isinstance(v, str) and v for v in request.raw_args):
            raise SkillError("mode=raw 'raw_args' must contain non-empty strings")
    elif request.mode == "shell":
        if not isinstance(request.shell_line, str) or not request.shell_line.strip():
            raise SkillError("mode=shell requires non-empty string field 'shell_line'")
    elif request.mode == "op":
        if request.operation is None:
            raise SkillError("mode=op requires field 'operation'")
    elif request.mode == "workflow":
        if request.workflow is None and not request.workflow_path:
            raise SkillError("mode=workflow requires 'workflow' or 'workflow_path'")
    elif request.mode == "external-primer-handoff":
        if not isinstance(request.state_path, str) or not request.state_path.strip():
            raise SkillError(
                "mode=external-primer-handoff requires explicit non-empty state_path"
            )
        request.state_path = request.state_path.strip()
        request.external_primer_handoff = _normalise_external_primer_handoff(
            request.external_primer_handoff
        )
    elif request.mode in {
        "construct-reasoning-list-inspections",
        "construct-reasoning-run-inspection",
    }:
        _normalise_construct_reasoning_inspection_request(request)
    elif request.mode in {"primer-design", "qpcr-design"}:
        _normalise_primer_design_request(request)
    elif request.mode in {
        "primer-seed-from-feature",
        "primer-seed-from-splicing",
        "qpcr-seed-from-feature",
    }:
        _normalise_seq_feature_request(request)
    elif request.mode == "qpcr-seed-from-splicing":
        _normalise_qpcr_seed_request(request)
    elif request.mode in {
        "primer-report-list",
        "primer-report-show",
        "primer-report-export",
        "qpcr-report-list",
        "qpcr-report-show",
        "qpcr-report-export",
    }:
        _normalise_primer_report_request(request)
    elif request.mode == "primer-preflight":
        _normalise_primer_backend_options(request)
    elif request.mode in {"cdna-pcr-test", "cdna-qpcr-test"}:
        _normalise_cdna_assay_test_request(request)
    elif request.mode == "protein-residue-genomic-coordinates":
        _normalise_protein_residue_genomic_coordinate_request(request)
    elif request.mode == "transcript-qpcr-panel":
        _normalise_transcript_qpcr_panel_request(request)
    elif request.mode.startswith("restriction-cloning-"):
        _normalise_restriction_cloning_request(request)
    elif request.mode == "pcr-protocol-cartoon":
        _normalise_pcr_protocol_cartoon_request(request)
    elif request.mode == "gene-protein-2d-gel":
        _normalise_gene_protein_2d_gel_request(request)
    elif request.mode in {"exon-skip-plan", "exon-skip-materialize"}:
        _normalise_exon_skip_request(request)
    elif request.mode == "agent-plan":
        if not isinstance(request.system_id, str) or not request.system_id.strip():
            raise SkillError("mode=agent-plan requires non-empty string field 'system_id'")
        if not isinstance(request.prompt, str) or not request.prompt.strip():
            raise SkillError("mode=agent-plan requires non-empty string field 'prompt'")
    elif request.mode == "agent-execute-plan":
        if request.plan is None and not request.plan_path:
            raise SkillError("mode=agent-execute-plan requires 'plan' or 'plan_path'")
        if not isinstance(request.candidate_id, str) or not request.candidate_id.strip():
            raise SkillError(
                "mode=agent-execute-plan requires non-empty string field 'candidate_id'"
            )
    for field_name in (
        "catalog_path",
        "base_url",
        "model",
        "workflow_path",
        "plan_path",
        "plan_id",
        "graph_id",
        "fact_id",
        "annotation_id",
        "system_id",
        "prompt",
        "candidate_id",
        "evidence_id",
        "summary_id",
        "action_kind",
        "action_id",
        "dotplot_id",
        "render_svg_path",
        "output_prefix",
    ):
        value = getattr(request, field_name)
        if value is not None and (not isinstance(value, str) or not value.strip()):
            raise SkillError(f"{field_name} must be a non-empty string when present")
    for field_name in (
        "connect_timeout_secs",
        "read_timeout_secs",
        "max_retries",
        "max_response_bytes",
        "max_candidates",
        "word_size",
        "step_bp",
        "tile_bp",
    ):
        value = getattr(request, field_name)
        if value is None:
            continue
        try:
            coerced = int(value)
        except (TypeError, ValueError) as e:
            raise SkillError(f"{field_name} must be an integer when present") from e
        if coerced < 0:
            raise SkillError(f"{field_name} must be >= 0 when present")
        setattr(request, field_name, coerced)
    for field_name in ("include_state_summary", "allow_mutating_candidates", "confirm"):
        value = getattr(request, field_name)
        if value is not None and not isinstance(value, bool):
            raise SkillError(f"{field_name} must be a boolean when present")
    if request.ensure_reference_prepared is not None:
        ensure = request.ensure_reference_prepared
        if not isinstance(ensure, dict):
            raise SkillError(
                "ensure_reference_prepared must be an object when present"
            )
        genome_id = str(ensure.get("genome_id", "")).strip()
        if not genome_id:
            raise SkillError(
                "ensure_reference_prepared requires non-empty string field 'genome_id'"
            )
        status_timeout_secs = int(ensure.get("status_timeout_secs", 300))
        prepare_timeout_secs = int(ensure.get("prepare_timeout_secs", 7200))
        if status_timeout_secs <= 0:
            raise SkillError(
                "ensure_reference_prepared.status_timeout_secs must be > 0"
            )
        if prepare_timeout_secs <= 0:
            raise SkillError(
                "ensure_reference_prepared.prepare_timeout_secs must be > 0"
            )
        for optional_path_field in ("catalog_path", "cache_dir"):
            value = ensure.get(optional_path_field)
            if value is not None and (not isinstance(value, str) or not value.strip()):
                raise SkillError(
                    f"ensure_reference_prepared.{optional_path_field} must be a non-empty string when present"
                )
        request.ensure_reference_prepared = EnsureReferencePrepared(
            genome_id=genome_id,
            catalog_path=ensure.get("catalog_path"),
            cache_dir=ensure.get("cache_dir"),
            status_timeout_secs=status_timeout_secs,
            prepare_timeout_secs=prepare_timeout_secs,
        )
    return request


def _json_arg(value: Any) -> str:
    if isinstance(value, str):
        return value
    return json.dumps(value, ensure_ascii=True, separators=(",", ":"))


def _gene_protein_2d_gel_workflow(request: Request) -> dict[str, Any]:
    gene_symbol = request.gene_symbol or "gene"
    species = request.species or "homo_sapiens"
    gene_id = _safe_id_component(gene_symbol, "gene")
    entry_id = f"{gene_id}_ensembl_gene"
    locus_id = f"{gene_id}_ensembl_locus"
    report_id = Path(_gene_protein_2d_gel_svg_path(request)).stem
    return {
        "run_id": f"clawbio_{report_id}",
        "ops": [
            {
                "FetchEnsemblGene": {
                    "query": gene_symbol,
                    "species": species,
                    "entry_id": entry_id,
                }
            },
            {
                "ImportEnsemblGeneSequence": {
                    "entry_id": entry_id,
                    "output_id": locus_id,
                }
            },
            {
                "DeriveProteinSequences": {
                    "seq_id": locus_id,
                    "feature_query": {
                        "seq_id": locus_id,
                        "include_source": False,
                        "include_qualifiers": False,
                        "kind_in": ["mRNA"],
                        "kind_not_in": [],
                        "range_relation": "overlap",
                        "strand": "any",
                        "qualifier_filters": [
                            {
                                "key": "biotype",
                                "value_contains": "protein_coding",
                                "value_regex": "^protein_coding$",
                                "case_sensitive": True,
                            }
                        ],
                        "sort_by": "feature_id",
                        "descending": False,
                        "offset": 0,
                    },
                    "output_prefix": f"{gene_id}_ensembl_protein",
                    "report_id": report_id,
                }
            },
            {
                "RenderProtein2dGelSvg": {
                    "report_id": report_id,
                    "path": _gene_protein_2d_gel_svg_path(request),
                    "ladders": request.ladders or ["Protein Ladder 10-100 kDa"],
                }
            },
        ],
    }


def _protein_residue_genomic_coordinate_operation(request: Request) -> dict[str, Any]:
    payload: dict[str, Any] = {
        "seq_id": request.seq_id,
        "residue_start_1based": request.residue_start_1based,
        "residue_end_1based": request.residue_end_1based,
    }
    if request.transcript_id:
        payload["transcript_id"] = request.transcript_id
    return {"QueryProteinResidueGenomicCoordinates": payload}


def _request_json_arg(request: Request) -> str:
    if request.request_json is not None:
        return _json_arg(request.request_json)
    if request.operation is not None:
        return _json_arg(request.operation)
    raise SkillError(f"mode={request.mode} requires 'request_json' or 'operation'")


def _primer_backend_tokens(request: Request) -> list[str]:
    tokens: list[str] = []
    if request.backend:
        tokens.extend(["--backend", request.backend])
    if request.primer3_executable:
        tokens.extend(["--primer3-exec", request.primer3_executable])
    return tokens


def _append_optional_int(tokens: list[str], flag: str, value: int | None) -> None:
    if value is not None:
        tokens.extend([flag, str(value)])


def _build_primer_mode_shell_line(request: Request) -> str:
    mode = request.mode
    tokens: list[str]
    if mode == "primer-preflight":
        tokens = ["primers", "preflight", *_primer_backend_tokens(request)]
    elif mode == "primer-seed-from-feature":
        tokens = ["primers", "seed-from-feature", request.seq_id or "SEQ_ID", str(request.source_feature_id)]
    elif mode == "primer-seed-from-splicing":
        tokens = ["primers", "seed-from-splicing", request.seq_id or "SEQ_ID", str(request.source_feature_id)]
    elif mode == "primer-design":
        tokens = ["primers", "design", _request_json_arg(request), *_primer_backend_tokens(request)]
    elif mode == "primer-report-list":
        tokens = ["primers", "list-reports"]
    elif mode == "primer-report-show":
        tokens = ["primers", "show-report", request.report_id or "REPORT_ID"]
    elif mode == "primer-report-export":
        tokens = [
            "primers",
            "export-report",
            request.report_id or "REPORT_ID",
            request.output_path or "primer_report.json",
        ]
    elif mode == "qpcr-seed-from-feature":
        tokens = ["primers", "seed-qpcr-from-feature", request.seq_id or "SEQ_ID", str(request.source_feature_id)]
    elif mode == "qpcr-seed-from-splicing":
        tokens = ["primers", "seed-qpcr-from-splicing", request.seq_id or "SEQ_ID", str(request.source_feature_id)]
        if request.qpcr_mode:
            tokens.extend(["--mode", request.qpcr_mode])
        if request.transcript_id:
            tokens.extend(["--transcript-id", request.transcript_id])
        if request.specificity_evidence:
            tokens.extend(["--specificity-evidence", request.specificity_evidence])
    elif mode == "qpcr-design":
        tokens = ["primers", "design-qpcr", _request_json_arg(request), *_primer_backend_tokens(request)]
    elif mode == "qpcr-report-list":
        tokens = ["primers", "list-qpcr-reports"]
    elif mode == "qpcr-report-show":
        tokens = ["primers", "show-qpcr-report", request.report_id or "QPCR_REPORT_ID"]
    elif mode == "qpcr-report-export":
        tokens = [
            "primers",
            "export-qpcr-report",
            request.report_id or "QPCR_REPORT_ID",
            request.output_path or "qpcr_report.json",
        ]
    elif mode in {"cdna-pcr-test", "cdna-qpcr-test"}:
        subcommand = "test-cdna-qpcr" if mode == "cdna-qpcr-test" else "test-cdna-pcr"
        tokens = [
            "primers",
            subcommand,
            request.seq_id or "SEQ_ID",
            str(request.source_feature_id),
            "--forward",
            request.forward_primer or "FORWARD",
            "--reverse",
            request.reverse_primer or "REVERSE",
        ]
        if mode == "cdna-qpcr-test":
            tokens.extend(["--probe", request.probe or "PROBE"])
        if request.transcript_id:
            tokens.extend(["--transcript-id", request.transcript_id])
        if request.transcript_order:
            tokens.extend(["--transcript-order", request.transcript_order])
        if request.map_coordinate_mode:
            tokens.extend(["--map-coordinate-mode", request.map_coordinate_mode])
        _append_optional_int(tokens, "--min-amplicon-bp", request.min_amplicon_bp)
        _append_optional_int(tokens, "--max-amplicon-bp", request.max_amplicon_bp)
        _append_optional_int(tokens, "--max-mismatches", request.max_mismatches)
        _append_optional_int(
            tokens,
            "--require-3prime-exact-bases",
            request.require_3prime_exact_bases,
        )
        if request.output_path:
            tokens.extend(["--path", request.output_path])
        if request.svg_path:
            tokens.extend(["--svg", request.svg_path])
        if request.materialize_products:
            tokens.append("--materialize-products")
        if request.product_output_prefix:
            tokens.extend(["--product-output-prefix", request.product_output_prefix])
        if request.product_gel_svg_path:
            tokens.extend(["--product-gel-svg", request.product_gel_svg_path])
        for ladder in request.product_gel_ladders or []:
            tokens.extend(["--product-gel-ladder", ladder])
    elif mode == "restriction-cloning-pcr-handoff":
        tokens = ["primers", "prepare-restriction-cloning", _request_json_arg(request)]
    elif mode == "restriction-cloning-pcr-handoff-seed":
        tokens = [
            "primers",
            "seed-restriction-cloning-handoff",
            request.report_id or "PRIMER_REPORT_ID",
            request.vector_seq_id or "VECTOR_SEQ_ID",
        ]
        if request.pair_rank is not None:
            tokens.extend(["--pair-rank", str(request.pair_rank)])
        if request.handoff_mode:
            tokens.extend(["--mode", request.handoff_mode])
        if request.forward_enzyme:
            tokens.extend(["--forward-enzyme", request.forward_enzyme])
        if request.reverse_enzyme:
            tokens.extend(["--reverse-enzyme", request.reverse_enzyme])
        if request.forward_leader_5prime:
            tokens.extend(["--forward-leader", request.forward_leader_5prime])
        if request.reverse_leader_5prime:
            tokens.extend(["--reverse-leader", request.reverse_leader_5prime])
    elif mode == "restriction-cloning-vector-suggestions":
        tokens = ["primers", "restriction-cloning-vector-suggestions", request.seq_id or "SEQ_ID"]
    elif mode == "restriction-cloning-handoff-list":
        tokens = ["primers", "list-restriction-cloning-handoffs"]
    elif mode == "restriction-cloning-handoff-show":
        tokens = [
            "primers",
            "show-restriction-cloning-handoff",
            request.report_id or "HANDOFF_REPORT_ID",
        ]
    elif mode == "restriction-cloning-handoff-export":
        tokens = [
            "primers",
            "export-restriction-cloning-handoff",
            request.report_id or "HANDOFF_REPORT_ID",
            request.output_path or "restriction_cloning_handoff.json",
        ]
    else:
        raise SkillError(f"Unsupported primer/PCR mode '{mode}'")
    return shlex.join(tokens)


def _build_pcr_protocol_cartoon_shell_line(request: Request) -> str:
    tokens = [
        "protocol-cartoon",
        "render-svg",
        request.protocol_id or "pcr.assay.qpcr",
        request.output_path or "artifacts/pcr.protocol.svg",
    ]
    return shlex.join(tokens)


def _build_exon_skip_shell_line(request: Request) -> str:
    if request.mode == "exon-skip-plan":
        tokens = [
            "transcripts",
            "exon-skip-plan",
            request.seq_id or "SEQ_ID",
            "--feature-id",
            str(request.transcript_feature_id if request.transcript_feature_id is not None else 0),
        ]
        for candidate_id in request.skip_candidate_ids or []:
            tokens.extend(["--skip", candidate_id])
        for interval in request.skip_intervals_1based or []:
            tokens.extend(
                ["--skip", f"{interval['start_1based']}..{interval['end_1based']}"]
            )
        for interval in request.overlap_intervals_1based or []:
            tokens.extend(
                ["--overlap", f"{interval['start_1based']}..{interval['end_1based']}"]
            )
        for value in request.length_mod3_values or []:
            tokens.extend(["--length-mod3", str(value)])
        for value in request.coding_mod3_values or []:
            tokens.extend(["--coding-mod3", str(value)])
        for context in request.coding_contexts or []:
            tokens.extend(["--coding-context", context])
        for kind in request.cds_phase_entry_kinds or []:
            tokens.extend(["--phase-entry", kind])
        if request.feature_query is not None:
            tokens.extend(
                [
                    "--feature-query-json",
                    request.feature_query
                    if isinstance(request.feature_query, str)
                    else json.dumps(request.feature_query, separators=(",", ":")),
                ]
            )
        if request.plan_id:
            tokens.extend(["--plan-id", request.plan_id])
        return shlex.join(tokens)

    tokens = ["transcripts", "exon-skip-materialize", request.plan_id or "PLAN_ID"]
    for candidate_id in request.candidate_ids or []:
        tokens.extend(["--candidate-id", candidate_id])
    if request.output_prefix:
        tokens.extend(["--output-prefix", request.output_prefix])
    for item in request.return_items or []:
        tokens.extend(["--return", item])
    return shlex.join(tokens)


def _build_construct_reasoning_inspection_shell_line(request: Request) -> str:
    if request.mode == "construct-reasoning-list-inspections":
        tokens = [
            "construct-reasoning",
            "list-inspection-actions",
            request.graph_id or "GRAPH_ID",
        ]
        for flag, value in (
            ("--fact-id", request.fact_id),
            ("--annotation-id", request.annotation_id),
            ("--candidate-id", request.candidate_id),
            ("--evidence-id", request.evidence_id),
            ("--seq-id", request.seq_id),
            ("--action-kind", request.action_kind),
            ("--summary-id", request.summary_id),
        ):
            if value:
                tokens.extend([flag, value])
        return shlex.join(tokens)

    tokens = [
        "construct-reasoning",
        "run-inspection-action",
        request.graph_id or "GRAPH_ID",
        request.action_id or "ACTION_ID",
    ]
    _append_optional_int(tokens, "--word-size", request.word_size)
    _append_optional_int(tokens, "--step", request.step_bp)
    _append_optional_int(tokens, "--max-mismatches", request.max_mismatches)
    _append_optional_int(tokens, "--tile-bp", request.tile_bp)
    if request.dotplot_id:
        tokens.extend(["--id", request.dotplot_id])
    if request.render_svg_path:
        tokens.extend(["--render-svg", request.render_svg_path])
    return shlex.join(tokens)


def _preview_text(text: str, *, max_lines: int = 6, max_chars: int = 600) -> str | None:
    trimmed = text.strip()
    if not trimmed:
        return None
    lines = trimmed.splitlines()
    preview = "\n".join(lines[:max_lines])
    if len(preview) > max_chars:
        preview = preview[: max_chars - 3].rstrip() + "..."
    elif len(lines) > max_lines or len(trimmed) > len(preview):
        preview = preview.rstrip() + "..."
    return preview


def _infer_command_compatibility_note(
    command: list[str] | None,
    stderr_text: str,
) -> str | None:
    if not command:
        return None

    stderr = stderr_text.strip()
    if not stderr:
        return None

    requested = ""
    if len(command) >= 3 and command[1] == "shell":
        requested = command[2]
    elif len(command) >= 3:
        requested = " ".join(command[1:3])
    elif len(command) >= 2:
        requested = command[1]

    requested = requested.strip()
    if not requested:
        return None

    command_unknown = (
        "Unknown shell command" in stderr
        or "unrecognized subcommand" in stderr
        or "unexpected argument" in stderr
        or "Found argument" in stderr
    )
    if not command_unknown:
        return None

    if requested.startswith(("services ", "resources ", "helpers status", "genomes status")):
        return (
            "This GENtle CLI looks older than the current ClawBio skill scaffold and "
            "may not support the requested status route yet. Update the installed "
            "`gentle_cli` on PATH, or set `GENTLE_CLI_CMD` to "
            "`gentle_local_checkout_cli.sh` for an updated GENtle checkout."
        )
    return (
        "This GENtle CLI may be older than the current ClawBio skill scaffold and "
        "may not support the requested command yet. Update the installed "
        "`gentle_cli` on PATH, or point `GENTLE_CLI_CMD` at the updated GENtle "
        "checkout launcher."
    )


def _build_failure_summary(
    *,
    stage: str,
    step: dict[str, Any] | None,
    execution_cwd: Path | None,
    note: str | None = None,
) -> dict[str, Any]:
    command = step.get("command") if isinstance(step, dict) else None
    effective_note = note
    if effective_note is None and isinstance(step, dict):
        effective_note = _infer_command_compatibility_note(
            command,
            step.get("stderr", ""),
        )
    summary = {
        "stage": stage,
        "note": effective_note,
        "command": command,
        "command_text": _format_command_text(command),
        "execution_cwd": (str(execution_cwd) if execution_cwd is not None else None),
        "exit_code": (step.get("exit_code") if isinstance(step, dict) else None),
        "stderr_preview": (
            _preview_text(step.get("stderr", "")) if isinstance(step, dict) else None
        ),
        "stdout_preview": (
            _preview_text(step.get("stdout", "")) if isinstance(step, dict) else None
        ),
    }
    return summary


def _build_failure_message(
    *,
    headline: str,
    failure_summary: dict[str, Any] | None,
) -> str:
    if not failure_summary:
        return headline

    parts = [headline]
    note = failure_summary.get("note")
    if isinstance(note, str) and note.strip():
        parts.append(note.strip())
    command_text = failure_summary.get("command_text")
    if isinstance(command_text, str) and command_text and command_text != "(none)":
        parts.append(f"Command: `{command_text}`.")
    execution_cwd = failure_summary.get("execution_cwd")
    if isinstance(execution_cwd, str) and execution_cwd:
        parts.append(f"Cwd: `{execution_cwd}`.")
    exit_code = failure_summary.get("exit_code")
    if exit_code is not None:
        parts.append(f"Exit code: `{exit_code}`.")
    stderr_preview = failure_summary.get("stderr_preview")
    stdout_preview = failure_summary.get("stdout_preview")
    if isinstance(stderr_preview, str) and stderr_preview:
        parts.append(f"Stderr preview: `{stderr_preview}`.")
    elif isinstance(stdout_preview, str) and stdout_preview:
        parts.append(f"Stdout preview: `{stdout_preview}`.")
    parts.append("Full stdout/stderr are recorded in `report.md` and `result.json`.")
    return " ".join(parts)


def _resolve_execution_cwd(
    request: Request,
    resolution: CliResolution | None,
    script_path: Path,
) -> Path:
    if request.mode == "workflow" and request.workflow_path:
        resolved_workflow_path = Path(_resolve_request_path(request.workflow_path, script_path))
        if resolved_workflow_path.exists():
            repo_root = _find_repo_root_from_path(resolved_workflow_path)
            if repo_root is not None:
                return repo_root
            return resolved_workflow_path.parent
    if resolution and resolution.cwd:
        return Path(resolution.cwd)
    return Path.cwd()


def _resolve_expected_artifact(raw_path: str, execution_cwd: Path) -> Path:
    path = Path(raw_path)
    if path.is_absolute():
        return path
    return execution_cwd / path


def _prepare_expected_artifact_parent_dirs(
    request: Request,
    execution_cwd: Path,
) -> None:
    if not request.expected_artifacts:
        return
    for raw_path in request.expected_artifacts:
        path = Path(raw_path)
        if path.is_absolute():
            continue
        parent = (execution_cwd / path).parent
        if parent == execution_cwd:
            continue
        parent.mkdir(parents=True, exist_ok=True)


def _safe_artifact_destination(raw_path: str) -> Path:
    normalized = raw_path.replace("\\", "/")
    parts: list[str] = []
    for part in Path(normalized).parts:
        if part in ("", ".", "..", "/", "\\"):
            continue
        if part.endswith(":"):
            continue
        parts.append(part)
    if not parts:
        return Path("artifact")
    return Path(*parts)


def _bundle_relative_path(output_dir: Path, destination: Path) -> str:
    try:
        return destination.resolve().relative_to(output_dir.resolve()).as_posix()
    except ValueError:
        return destination.name


def _artifact_record(
    *,
    declared_path: str,
    source_path: Path,
    destination: Path,
    output_dir: Path,
    derived_from: str | None = None,
) -> dict[str, Any]:
    record: dict[str, Any] = {
        "declared_path": declared_path,
        "bundle_path": _bundle_relative_path(output_dir, destination),
        "source_path": str(source_path.resolve()),
        "copied_path": str(destination.resolve()),
    }
    if derived_from:
        record["derived_from"] = derived_from
    return record


def _copy_collected_artifacts(
    request: Request,
    output_dir: Path,
    execution_cwd: Path,
) -> list[dict[str, Any]]:
    collected: list[dict[str, Any]] = []
    if not request.expected_artifacts:
        return collected

    generated_dir = output_dir / "generated"
    for raw_path in request.expected_artifacts:
        source_path = _resolve_expected_artifact(raw_path, execution_cwd)
        if not source_path.exists() or not source_path.is_file():
            continue

        destination = generated_dir / _safe_artifact_destination(raw_path)
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source_path, destination)
        collected.append(
            _artifact_record(
                declared_path=raw_path,
                source_path=source_path,
                destination=destination,
                output_dir=output_dir,
            )
        )
    return collected


def _png_declared_path_from_svg(raw_path: str) -> str:
    normalized = raw_path.replace("\\", "/")
    return str(Path(normalized).with_suffix(".png")).replace("\\", "/")


def _preferred_artifact_id_for_png(raw_id: str | None, bundle_path: str) -> str:
    candidate = (raw_id or "").strip()
    if candidate:
        if candidate.endswith("_svg"):
            return candidate[:-4] + "_png"
        if candidate.endswith(".svg"):
            return candidate[:-4] + ".png"
        if not candidate.endswith("_png"):
            return candidate + "_png"
        return candidate
    stem = Path(bundle_path).stem
    return stem if stem.endswith("_png") else stem + "_png"


def _artifact_reference_keys(artifact: dict[str, Any]) -> set[str]:
    keys = set()
    for key in ("declared_path", "bundle_path", "path", "derived_from"):
        value = artifact.get(key)
        if isinstance(value, str) and value.strip():
            keys.add(value.strip())
    return keys


def _collected_svg_artifacts(
    collected_artifacts: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    return [
        artifact
        for artifact in collected_artifacts
        if isinstance(artifact, dict)
        and Path(str(artifact.get("copied_path", ""))).suffix.lower() == ".svg"
        and isinstance(artifact.get("declared_path"), str)
        and str(artifact["declared_path"]).strip()
    ]


def _best_first_svg_declared_paths(
    collected_artifacts: list[dict[str, Any]],
    preferred_artifacts: list[dict[str, Any]] | None,
) -> set[str]:
    svg_artifacts = _collected_svg_artifacts(collected_artifacts)
    if not svg_artifacts:
        return set()
    if len(svg_artifacts) == 1:
        return {str(svg_artifacts[0]["declared_path"]).strip()}

    preferred = [
        dict(artifact)
        for artifact in (preferred_artifacts or [])
        if isinstance(artifact, dict)
        and isinstance(artifact.get("path"), str)
        and str(artifact["path"]).strip().lower().endswith(".svg")
    ]
    preferred.sort(
        key=lambda artifact: (
            int(artifact.get("presentation_rank", 10**9))
            if isinstance(artifact.get("presentation_rank"), int)
            else 10**9,
            str(artifact.get("artifact_id", "")),
        )
    )
    for preferred_artifact in preferred:
        preferred_path = str(preferred_artifact["path"]).strip()
        for artifact in svg_artifacts:
            if preferred_path in _artifact_reference_keys(artifact):
                return {str(artifact["declared_path"]).strip()}

    if preferred_artifacts:
        return set()
    return {str(svg_artifacts[0]["declared_path"]).strip()}


def _rasterize_collected_svg_artifacts(
    collected_artifacts: list[dict[str, Any]],
    resolution: CliResolution,
    execution_cwd: Path,
    output_dir: Path,
    rasterize_declared_paths: set[str] | None = None,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    updated_collected = list(collected_artifacts)
    rasterized_pngs: list[dict[str, Any]] = []
    for artifact in collected_artifacts:
        copied_path = Path(str(artifact.get("copied_path", "")))
        if copied_path.suffix.lower() != ".svg":
            continue
        if rasterize_declared_paths is not None:
            declared_path = str(artifact.get("declared_path", "")).strip()
            bundle_path = str(artifact.get("bundle_path", "")).strip()
            if (
                declared_path not in rasterize_declared_paths
                and bundle_path not in rasterize_declared_paths
            ):
                continue
        png_path = copied_path.with_suffix(".png")
        run_result, step = _run_cli_command(
            resolution,
            [
                "svg-png",
                str(copied_path),
                str(png_path),
                "--scale",
                str(CLAWBIO_GRAPHICS_SCALE),
            ],
            execution_cwd,
            CLAWBIO_SVG_PNG_TIMEOUT_SECS,
        )
        if run_result.returncode != 0:
            failure_summary = _build_failure_summary(
                stage="artifact_svg_png",
                step=step,
                execution_cwd=execution_cwd,
            )
            raise SkillError(
                _build_failure_message(
                    headline=(
                        f"Failed to rasterize SVG artifact '{artifact.get('declared_path', copied_path.name)}' "
                        "into the PNG-first ClawBio bundle contract."
                    ),
                    failure_summary=failure_summary,
                )
            )
        png_artifact = _artifact_record(
            declared_path=_png_declared_path_from_svg(
                str(artifact.get("declared_path", copied_path.name))
            ),
            source_path=png_path,
            destination=png_path,
            output_dir=output_dir,
            derived_from=str(artifact.get("declared_path", copied_path.name)),
        )
        updated_collected.append(png_artifact)
        rasterized_pngs.append(png_artifact)
    return updated_collected, rasterized_pngs


def _rewrite_preferred_artifacts_for_png(
    collected_artifacts: list[dict[str, Any]],
    preferred_artifacts: list[dict[str, Any]] | None,
    rasterized_pngs: list[dict[str, Any]],
) -> list[dict[str, Any]] | None:
    if not rasterized_pngs:
        return preferred_artifacts

    svg_by_declared: dict[str, dict[str, Any]] = {
        str(artifact.get("declared_path", "")): artifact
        for artifact in collected_artifacts
        if str(artifact.get("copied_path", "")).lower().endswith(".svg")
    }
    png_by_source_path: dict[str, dict[str, Any]] = {}
    for artifact in rasterized_pngs:
        derived_from = str(artifact.get("derived_from", "")).strip()
        if not derived_from:
            continue
        png_by_source_path[derived_from] = artifact
        source_svg = svg_by_declared.get(derived_from)
        if source_svg is not None:
            png_by_source_path[str(source_svg.get("bundle_path", ""))] = artifact

    rewritten: list[dict[str, Any]] = []
    for artifact in preferred_artifacts or []:
        if not isinstance(artifact, dict):
            continue
        path = artifact.get("path")
        if not isinstance(path, str) or not path.strip():
            continue
        replacement = png_by_source_path.get(path)
        if replacement is not None:
            updated = dict(artifact)
            updated["path"] = replacement["bundle_path"]
            updated["artifact_id"] = _preferred_artifact_id_for_png(
                (
                    str(updated["artifact_id"])
                    if isinstance(updated.get("artifact_id"), str)
                    else None
                ),
                replacement["bundle_path"],
            )
            updated["derived_from"] = path
            rewritten.append(updated)
            continue
        if not path.lower().endswith(".svg"):
            rewritten.append(dict(artifact))

    if rewritten:
        return _single_best_preferred_artifact(rewritten)

    synthesized: list[dict[str, Any]] = []
    for idx, artifact in enumerate(rasterized_pngs):
        synthesized.append(
            {
                "artifact_id": _preferred_artifact_id_for_png(
                    None, str(artifact["bundle_path"])
                ),
                "path": artifact["bundle_path"],
                "caption": Path(str(artifact["bundle_path"]))
                .stem.replace("_", " ")
                .replace(".", " "),
                "recommended_use": "best_first_figure" if idx == 0 else "supporting_figure",
                "presentation_rank": idx,
                "is_best_first_artifact": idx == 0,
                "derived_from": artifact.get("derived_from"),
            }
        )
    return _single_best_preferred_artifact(synthesized) or None


def _single_best_preferred_artifact(
    artifacts: list[dict[str, Any]],
) -> list[dict[str, Any]] | None:
    valid = [
        dict(artifact)
        for artifact in artifacts
        if isinstance(artifact, dict)
        and isinstance(artifact.get("path"), str)
        and artifact["path"].strip()
    ]
    if not valid:
        return None
    valid.sort(
        key=lambda artifact: (
            int(artifact.get("presentation_rank", 10**9))
            if isinstance(artifact.get("presentation_rank"), int)
            else 10**9,
            str(artifact.get("artifact_id", "")),
        )
    )
    best = valid[0]
    best["presentation_rank"] = 0
    best["is_best_first_artifact"] = True
    best.setdefault("recommended_use", "best_first_figure")
    return [best]


def _parse_svg_dimension(value: str | None) -> float | None:
    if value is None:
        return None
    match = SVG_DIMENSION_RE.match(value)
    if not match:
        return None
    try:
        return float(match.group(1))
    except ValueError:
        return None


def _storyboard_caption_for_artifact(raw_path: str, preferred: dict[str, Any] | None) -> str:
    if preferred is not None and isinstance(preferred.get("caption"), str):
        caption = preferred["caption"].strip()
        if caption:
            return caption
    stem = Path(raw_path).stem.replace("_", " ").replace(".", " ")
    return " ".join(part for part in stem.split() if part)


def _storyboard_panel_title(raw_path: str, preferred: dict[str, Any] | None) -> str:
    lowered = raw_path.lower()
    if "context" in lowered:
        return "Genomic context"
    if "reference" in lowered and ("reporter" in lowered or "construct" in lowered):
        return "Reference allele reporter"
    if "alternate" in lowered and ("reporter" in lowered or "construct" in lowered):
        return "Alternate allele reporter"
    return _storyboard_caption_for_artifact(raw_path, preferred)


def _storyboard_panel_priority(raw_path: str) -> tuple[int, str]:
    lowered = raw_path.lower()
    if "context" in lowered:
        return (0, lowered)
    if "reference" in lowered:
        return (1, lowered)
    if "alternate" in lowered:
        return (2, lowered)
    return (3, lowered)


def _storyboard_panel_record(
    copied_artifact: dict[str, Any],
    preferred: dict[str, Any] | None,
) -> dict[str, Any] | None:
    copied_path = Path(copied_artifact["copied_path"])
    if copied_path.suffix.lower() != ".svg":
        return None
    try:
        text = copied_path.read_text(encoding="utf-8")
    except OSError:
        return None
    try:
        root = ET.fromstring(text)
    except ET.ParseError:
        return None

    width = None
    height = None
    view_box = root.get("viewBox")
    if view_box:
        parts = view_box.replace(",", " ").split()
        if len(parts) == 4:
            try:
                width = float(parts[2])
                height = float(parts[3])
            except ValueError:
                width = None
                height = None
    if width is None or height is None or width <= 0.0 or height <= 0.0:
        width = _parse_svg_dimension(root.get("width")) or 800.0
        height = _parse_svg_dimension(root.get("height")) or 600.0
    if width <= 0.0 or height <= 0.0:
        return None

    return {
        "declared_path": copied_artifact["declared_path"],
        "copied_path": copied_artifact["copied_path"],
        "artifact_id": (
            preferred.get("artifact_id")
            if preferred is not None and isinstance(preferred.get("artifact_id"), str)
            else Path(copied_artifact["declared_path"]).stem
        ),
        "caption": _storyboard_caption_for_artifact(copied_artifact["declared_path"], preferred),
        "panel_title": _storyboard_panel_title(copied_artifact["declared_path"], preferred),
        "priority": _storyboard_panel_priority(copied_artifact["declared_path"]),
        "width": width,
        "height": height,
        "data_uri": (
            "data:image/svg+xml;base64,"
            + base64.b64encode(text.encode("utf-8")).decode("ascii")
        ),
    }


def _storyboard_frame_rects(
    panel_count: int,
    variant_layout: bool,
) -> list[tuple[float, float, float, float]]:
    if variant_layout and panel_count >= 3:
        return [
            (64.0, 172.0, 900.0, 860.0),
            (1000.0, 172.0, 536.0, 404.0),
            (1000.0, 628.0, 536.0, 404.0),
        ]
    columns = 2 if panel_count > 1 else 1
    rows = max(1, math.ceil(panel_count / columns))
    margin = 64.0
    gutter = 36.0
    top = 172.0
    canvas_width = 1600.0
    canvas_height = 1100.0
    frame_width = (canvas_width - margin * 2 - gutter * (columns - 1)) / columns
    frame_height = (canvas_height - top - 72.0 - gutter * (rows - 1)) / rows
    rects: list[tuple[float, float, float, float]] = []
    for idx in range(panel_count):
        row = idx // columns
        col = idx % columns
        rects.append(
            (
                margin + col * (frame_width + gutter),
                top + row * (frame_height + gutter),
                frame_width,
                frame_height,
            )
        )
    return rects


def _write_storyboard_svg(
    title: str,
    subtitle: str,
    footer: str,
    panels: list[dict[str, Any]],
    output_path: Path,
) -> None:
    variant_layout = (
        len(panels) == 3
        and "context" in panels[0]["declared_path"].lower()
        and any("reference" in panel["declared_path"].lower() for panel in panels)
        and any("alternate" in panel["declared_path"].lower() for panel in panels)
    )
    frame_rects = _storyboard_frame_rects(len(panels), variant_layout)

    lines = [
        '<svg xmlns="http://www.w3.org/2000/svg" width="1600" height="1100" viewBox="0 0 1600 1100" role="img" aria-labelledby="title subtitle">',
        "  <defs>",
        '    <linearGradient id="bg" x1="0%" y1="0%" x2="100%" y2="100%">',
        '      <stop offset="0%" stop-color="#f8fafc" />',
        '      <stop offset="100%" stop-color="#e2e8f0" />',
        "    </linearGradient>",
        '    <filter id="card-shadow" x="-10%" y="-10%" width="130%" height="140%">',
        '      <feDropShadow dx="0" dy="10" stdDeviation="14" flood-color="#0f172a" flood-opacity="0.12" />',
        "    </filter>",
        "  </defs>",
        '  <rect width="1600" height="1100" fill="url(#bg)" />',
        '  <text id="title" x="64" y="78" font-family="Inter,Segoe UI,sans-serif" font-size="34" font-weight="700" fill="#0f172a">'
        + html.escape(title)
        + "</text>",
        '  <text id="subtitle" x="64" y="114" font-family="Inter,Segoe UI,sans-serif" font-size="16" fill="#334155">'
        + html.escape(subtitle)
        + "</text>",
        '  <text x="64" y="144" font-family="Inter,Segoe UI,sans-serif" font-size="13" fill="#475569">[clawbio] Deterministic GENtle figures bundled for a ClawBio/OpenClaw experimental follow-up reply.</text>',
    ]

    for panel, (x, y, width, height) in zip(panels, frame_rects):
        header_height = 76.0
        image_pad = 20.0
        image_box_x = x + image_pad
        image_box_y = y + header_height
        image_box_width = width - image_pad * 2
        image_box_height = height - header_height - image_pad
        scale = min(
            image_box_width / panel["width"],
            image_box_height / panel["height"],
        )
        render_width = panel["width"] * scale
        render_height = panel["height"] * scale
        render_x = image_box_x + (image_box_width - render_width) / 2.0
        render_y = image_box_y + (image_box_height - render_height) / 2.0
        lines.extend(
            [
                f'  <g filter="url(#card-shadow)">',
                f'    <rect x="{x:.1f}" y="{y:.1f}" width="{width:.1f}" height="{height:.1f}" rx="28" fill="#ffffff" />',
                f'  </g>',
                f'  <text x="{x + 24.0:.1f}" y="{y + 32.0:.1f}" font-family="Inter,Segoe UI,sans-serif" font-size="20" font-weight="700" fill="#0f172a">{html.escape(panel["panel_title"])}</text>',
                f'  <text x="{x + 24.0:.1f}" y="{y + 56.0:.1f}" font-family="Inter,Segoe UI,sans-serif" font-size="13" fill="#64748b">{html.escape(panel["caption"])}</text>',
                f'  <image href="{panel["data_uri"]}" x="{render_x:.1f}" y="{render_y:.1f}" width="{render_width:.1f}" height="{render_height:.1f}" preserveAspectRatio="xMidYMid meet" />',
            ]
        )

    lines.extend(
        [
            '  <text x="64" y="1070" font-family="Inter,Segoe UI,sans-serif" font-size="12" fill="#64748b">'
            + html.escape(footer)
            + "</text>",
            "</svg>",
        ]
    )
    output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _augment_artifacts_with_storyboard(
    request: Request,
    output_dir: Path,
    collected_artifacts: list[dict[str, Any]],
    preferred_artifacts: list[dict[str, Any]] | None,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]] | None]:
    preferred_by_path = {
        artifact["path"]: artifact
        for artifact in (preferred_artifacts or [])
        if isinstance(artifact, dict) and isinstance(artifact.get("path"), str)
    }
    panels = []
    for artifact in collected_artifacts:
        preferred = preferred_by_path.get(artifact["declared_path"])
        panel = _storyboard_panel_record(artifact, preferred)
        if panel is not None:
            panels.append(panel)
    panels.sort(key=lambda panel: panel["priority"])
    if len(panels) < 2:
        return collected_artifacts, preferred_artifacts

    workflow_hint = " ".join(
        value
        for value in (
            request.workflow_path or "",
            request.shell_line or "",
        )
        if value
    ).lower()
    variant_story = (
        any("context" in panel["declared_path"].lower() for panel in panels)
        and any("reference" in panel["declared_path"].lower() for panel in panels)
        and any("alternate" in panel["declared_path"].lower() for panel in panels)
    ) or ("variant" in workflow_hint or "luciferase" in workflow_hint)
    if request.presentation_profile == PCR_PRIMER_PRESENTATION_PROFILE:
        title = "[clawbio] PCR primer-design evidence report"
        subtitle = (
            "[clawbio] GENtle-generated assay figures arranged without "
            "recomputing primer or product evidence."
        )
        footer = (
            "[clawbio] Primer pairs use comparable lanes/columns and one shared gel model; "
            "source-prefixed claims remain traceable in claim_ledger.json."
        )
        panels = [
            {
                **panel,
                "caption": _prefix_claim_line(
                    panel["caption"],
                    "gentle_executable",
                ),
            }
            for panel in panels
        ]
    elif variant_story:
        title = "Variant-to-Synthetic-Biology assay storyboard"
        subtitle = (
            "GENtle bridges genomic variant context to engineered allele-specific "
            "reporter constructs for reproducible functional follow-up."
        )
        footer = (
            "Typical answer shape: locus context, assayable promoter fragment logic, "
            "and engineered allele-specific reporter constructs for synthetic-biology follow-up."
        )
    else:
        title = "GENtle graphical storyboard"
        subtitle = "One shareable figure assembled from the deterministic graphics in this run."
        footer = (
            "The storyboard arranges existing GENtle artifacts; it does not "
            "recalculate their scientific content."
        )

    storyboard_path = output_dir / "generated" / "clawbio_storyboard.svg"
    storyboard_path.parent.mkdir(parents=True, exist_ok=True)
    _write_storyboard_svg(title, subtitle, footer, panels[:4], storyboard_path)

    storyboard_artifact = _artifact_record(
        declared_path="generated/clawbio_storyboard.svg",
        source_path=storyboard_path,
        destination=storyboard_path,
        output_dir=output_dir,
    )
    updated_collected = [*collected_artifacts, storyboard_artifact]

    storyboard_entry = {
        "artifact_id": "clawbio_storyboard_svg",
        "path": "generated/clawbio_storyboard.svg",
        "caption": title,
        "recommended_use": "best_first_figure",
        "presentation_rank": 0,
        "is_best_first_artifact": True,
    }
    merged_preferred: list[dict[str, Any]] = [storyboard_entry]
    seen_paths = {storyboard_entry["path"]}

    source_entries = preferred_artifacts or []
    if not source_entries:
        source_entries = [
            {
                "artifact_id": panel["artifact_id"],
                "path": panel["declared_path"],
                "caption": panel["caption"],
                "recommended_use": "supporting_figure",
                "presentation_rank": idx + 1,
                "is_best_first_artifact": False,
            }
            for idx, panel in enumerate(panels)
        ]

    for idx, artifact in enumerate(source_entries, start=1):
        if not isinstance(artifact, dict):
            continue
        path = artifact.get("path")
        if not isinstance(path, str) or path in seen_paths:
            continue
        updated = dict(artifact)
        updated["is_best_first_artifact"] = False
        if not isinstance(updated.get("presentation_rank"), int):
            updated["presentation_rank"] = idx
        elif updated["presentation_rank"] <= 0:
            updated["presentation_rank"] = idx
        seen_paths.add(path)
        merged_preferred.append(updated)

    return updated_collected, merged_preferred


def _external_primer_batch_record(
    record: dict[str, Any], collection_id: str
) -> dict[str, Any]:
    annotations = dict(record["annotations"])
    annotations.update(
        {
            EXTERNAL_PRIMER_HANDOFF_ANNOTATION_RECORD_ID: record["record_id"],
            EXTERNAL_PRIMER_HANDOFF_ANNOTATION_ASSAY_PURPOSE: record[
                "assay_purpose"
            ],
            EXTERNAL_PRIMER_HANDOFF_ANNOTATION_COLLECTION_ID: collection_id,
        }
    )
    return {
        "source_kind": record["source_kind"],
        "provider": record["provider"],
        "catalogue_id": record["catalogue_id"],
        "source_url": record["source_url"],
        "claimed_accession": record["claimed_accession"],
        "aliases": record["aliases"],
        "forward_sequence_5_to_3": record["forward_sequence_5_to_3"],
        "reverse_sequence_5_to_3": record["reverse_sequence_5_to_3"],
        "claimed_target": record["claimed_target"],
        "validation_claims": record["validation_claims"],
        "annotations": dict(sorted(annotations.items())),
    }


def _external_primer_not_submitted_record(record: dict[str, Any]) -> dict[str, Any]:
    return {
        "record_id": record["record_id"],
        "assay_purpose": record["assay_purpose"],
        "record_kind": record["record_kind"],
        "status": "not_submitted",
        "status_reason": (
            "assay_purpose_not_supported_by_external_primer_pair_importer"
        ),
        "record": record,
    }


def _resolve_external_primer_state_path(raw_path: str, execution_cwd: Path) -> Path:
    path = Path(raw_path).expanduser()
    if not path.is_absolute():
        path = execution_cwd / path
    resolved = path.resolve()
    if not resolved.exists():
        raise SkillError(
            f"external-primer handoff state_path '{resolved}' does not exist"
        )
    if not resolved.is_file():
        raise SkillError(
            f"external-primer handoff state_path '{resolved}' is not a file"
        )
    return resolved


def _external_primer_handoff_shell_line(context: dict[str, Any]) -> str:
    target = context["request"]["target"]
    evaluation = context["request"]["evaluation"]
    tokens = [
        "primers",
        "import-external-pairs",
        context["batch_path"],
        target["seq_id"],
        str(target["source_feature_id"]),
        "--format",
        "json",
    ]
    optional_flags = (
        ("transcript_id", "--transcript-id"),
        ("min_amplicon_bp", "--min-amplicon-bp"),
        ("max_amplicon_bp", "--max-amplicon-bp"),
        ("max_mismatches", "--max-mismatches"),
        ("require_3prime_exact_bases", "--require-3prime-exact-bases"),
        ("transcript_order", "--transcript-order"),
        ("map_coordinate_mode", "--map-coordinate-mode"),
        ("specificity_target_genome_id", "--specificity-target-genome"),
        ("specificity_catalog_path", "--specificity-catalog"),
        ("specificity_cache_dir", "--specificity-cache-dir"),
    )
    values = dict(evaluation)
    values["transcript_id"] = target["transcript_id"]
    for key, flag in optional_flags:
        value = values.get(key)
        if value is not None:
            tokens.extend([flag, str(value)])
    tokens.extend(["--artifact-output-dir", context["artifact_output_dir"]])
    if evaluation["materialize_products"]:
        tokens.append("--materialize-products")
    for ladder in evaluation["product_gel_ladders"]:
        tokens.extend(["--product-gel-ladder", ladder])
    tokens.extend(["--path", context["report_path"]])
    return shlex.join(tokens)


def _prepare_external_primer_handoff(
    request: Request, output_dir: Path, execution_cwd: Path
) -> dict[str, Any]:
    handoff = request.external_primer_handoff
    if not isinstance(handoff, dict):
        raise SkillError("external-primer handoff request was not normalized")
    state_path = _resolve_external_primer_state_path(request.state_path or "", execution_cwd)
    request.state_path = str(state_path)
    state_sha256_before = _sha256_file_prefixed(state_path)
    expected_state_sha256 = handoff["target"]["expected_state_sha256"]
    if expected_state_sha256 is not None and expected_state_sha256 != state_sha256_before:
        raise SkillError(
            "external-primer handoff target state SHA-256 mismatch: expected "
            f"{expected_state_sha256}, observed {state_sha256_before}"
        )

    handoff_dir = output_dir / "external_primer_handoff"
    artifact_output_dir = handoff_dir / "scientific_artifacts"
    handoff_dir.mkdir(parents=True, exist_ok=True)
    artifact_output_dir.mkdir(parents=True, exist_ok=True)
    request_snapshot_path = handoff_dir / "canonical_request.json"
    request_snapshot_path.write_bytes(
        json.dumps(handoff, indent=2, sort_keys=True, ensure_ascii=True).encode("utf-8")
        + b"\n"
    )

    submitted_records = [
        record
        for record in handoff["records"]
        if record["assay_purpose"] in EXTERNAL_PRIMER_HANDOFF_SUBMITTED_PURPOSES
    ]
    not_submitted_records = [
        _external_primer_not_submitted_record(record)
        for record in handoff["records"]
        if record["assay_purpose"]
        in EXTERNAL_PRIMER_HANDOFF_NOT_SUBMITTED_PURPOSES
    ]
    batch_path = handoff_dir / "submitted_external_primer_pairs.json"
    report_path = handoff_dir / "external_primer_import_report.json"
    batch_payload: dict[str, Any] | None = None
    batch_file_sha256: str | None = None
    if submitted_records:
        batch_payload = {
            "schema": EXTERNAL_PRIMER_PAIR_BATCH_SCHEMA,
            "batch_id": handoff["collection_id"],
            "pairs": [
                _external_primer_batch_record(record, handoff["collection_id"])
                for record in submitted_records
            ],
        }
        batch_path.write_bytes(
            json.dumps(
                batch_payload,
                indent=2,
                sort_keys=True,
                ensure_ascii=True,
            ).encode("utf-8")
            + b"\n"
        )
        batch_file_sha256 = _sha256_file_prefixed(batch_path)

    context = {
        "request": handoff,
        "canonical_request_sha256": _sha256_prefixed_json(handoff),
        "request_snapshot_path": str(request_snapshot_path),
        "request_snapshot_sha256": _sha256_file_prefixed(request_snapshot_path),
        "state_path": str(state_path),
        "state_sha256_before": state_sha256_before,
        "expected_state_sha256": expected_state_sha256,
        "submitted_records": submitted_records,
        "not_submitted_records": not_submitted_records,
        "batch": batch_payload,
        "batch_path": str(batch_path),
        "batch_file_sha256": batch_file_sha256,
        "report_path": str(report_path),
        "artifact_output_dir": str(artifact_output_dir),
    }
    if submitted_records:
        context["shell_line"] = _external_primer_handoff_shell_line(context)
        request.shell_line = context["shell_line"]
    return context


def _build_cli_args(request: Request, script_path: Path) -> list[str]:
    args: list[str] = []
    if request.state_path:
        args.extend(["--state", request.state_path])
    if request.mode == "capabilities":
        args.append("capabilities")
    elif request.mode == "version":
        args.append("--version")
    elif request.mode == "state-summary":
        args.append("state-summary")
    elif request.mode in {"shell", "external-primer-handoff"}:
        if not isinstance(request.shell_line, str) or not request.shell_line.strip():
            raise SkillError(f"mode={request.mode} has no prepared shell command")
        args.extend(["shell", request.shell_line.strip()])
    elif request.mode == "op":
        args.extend(["op", _json_arg(request.operation)])
    elif request.mode == "workflow":
        if request.workflow_path:
            resolved_workflow_path = _resolve_request_path(request.workflow_path, script_path)
            args.extend(["workflow", f"@{resolved_workflow_path}"])
        else:
            args.extend(["workflow", _json_arg(request.workflow)])
    elif request.mode in {
        "construct-reasoning-list-inspections",
        "construct-reasoning-run-inspection",
    }:
        args.extend(["shell", _build_construct_reasoning_inspection_shell_line(request)])
    elif request.mode in PRIMER_SHELL_REQUEST_MODES:
        args.extend(["shell", _build_primer_mode_shell_line(request)])
    elif request.mode == "protein-residue-genomic-coordinates":
        args.extend(["op", _json_arg(_protein_residue_genomic_coordinate_operation(request))])
    elif request.mode == "transcript-qpcr-panel":
        tokens = [
            "primers",
            "transcript-qpcr-panel",
            request.seq_id or "SEQ_ID",
            str(request.source_feature_id if request.source_feature_id is not None else 0),
            request.shared_qpcr_report_id or "QPCR_REPORT_ID",
        ]
        if request.output_path:
            tokens.extend(["--path", request.output_path])
        args.extend(["shell", shlex.join(tokens)])
    elif request.mode == "pcr-protocol-cartoon":
        args.extend(["shell", _build_pcr_protocol_cartoon_shell_line(request)])
    elif request.mode == "gene-protein-2d-gel":
        args.extend(["workflow", _json_arg(_gene_protein_2d_gel_workflow(request))])
    elif request.mode in {"exon-skip-plan", "exon-skip-materialize"}:
        args.extend(["shell", _build_exon_skip_shell_line(request)])
    elif request.mode == "agent-plan":
        tokens = ["agents", "plan", request.system_id.strip(), "--prompt", request.prompt.strip()]
        if request.catalog_path:
            tokens.extend(["--catalog", request.catalog_path])
        if request.base_url:
            tokens.extend(["--base-url", request.base_url])
        if request.model:
            tokens.extend(["--model", request.model])
        if request.connect_timeout_secs:
            tokens.extend(
                ["--connect-timeout-secs", str(request.connect_timeout_secs)]
            )
        if request.read_timeout_secs:
            tokens.extend(["--read-timeout-secs", str(request.read_timeout_secs)])
        if request.max_retries is not None:
            tokens.extend(["--max-retries", str(request.max_retries)])
        if request.max_response_bytes:
            tokens.extend(["--max-response-bytes", str(request.max_response_bytes)])
        if request.max_candidates:
            tokens.extend(["--max-candidates", str(request.max_candidates)])
        if request.include_state_summary is False:
            tokens.append("--no-state-summary")
        if request.allow_mutating_candidates is False:
            tokens.append("--no-mutating-candidates")
        args.extend(["shell", shlex.join(tokens)])
    elif request.mode == "agent-execute-plan":
        plan_arg = (
            f"@{_resolve_request_path(request.plan_path, script_path)}"
            if request.plan_path
            else _json_arg(request.plan)
        )
        tokens = [
            "agents",
            "execute-plan",
            plan_arg,
            "--candidate-id",
            request.candidate_id.strip(),
        ]
        if request.confirm:
            tokens.append("--confirm")
        args.extend(["shell", shlex.join(tokens)])
    elif request.mode == "raw":
        args.extend(request.raw_args or [])
    else:
        raise SkillError(f"Unsupported mode '{request.mode}'")
    return args


def _build_reference_status_args(reference: EnsureReferencePrepared) -> list[str]:
    args = ["genomes", "status", reference.genome_id]
    if reference.catalog_path:
        args.extend(["--catalog", reference.catalog_path])
    if reference.cache_dir:
        args.extend(["--cache-dir", reference.cache_dir])
    return args


def _build_reference_prepare_args(reference: EnsureReferencePrepared) -> list[str]:
    args = [
        "genomes",
        "prepare",
        reference.genome_id,
        "--timeout-secs",
        str(reference.prepare_timeout_secs),
    ]
    if reference.catalog_path:
        args.extend(["--catalog", reference.catalog_path])
    if reference.cache_dir:
        args.extend(["--cache-dir", reference.cache_dir])
    return args


def _run_cli_command(
    resolution: CliResolution,
    cli_args: list[str],
    execution_cwd: Path,
    timeout_secs: int,
) -> tuple[subprocess.CompletedProcess[str], dict[str, Any]]:
    command = resolution.argv_prefix + cli_args
    started_utc = _now_utc_iso()
    popen_kwargs: dict[str, Any] = {
        "cwd": execution_cwd,
        "stdout": subprocess.PIPE,
        "stderr": subprocess.PIPE,
        "text": False,
    }
    if os.name == "posix":
        popen_kwargs["start_new_session"] = True
    elif os.name == "nt":
        popen_kwargs["creationflags"] = getattr(
            subprocess, "CREATE_NEW_PROCESS_GROUP", 0
        )
    process = subprocess.Popen(
        command,
        **popen_kwargs,
    )
    forwarded_signals: list[dict[str, Any]] = []
    previous_handlers: dict[int, Any] = {}

    def forward_signal(signum: int, _frame: Any) -> None:
        record = {
            "signal": signum,
            "signal_name": signal.Signals(signum).name,
            "forwarded": False,
            "error": None,
        }
        try:
            if os.name == "posix":
                os.killpg(process.pid, signum)
            else:
                process.send_signal(signum)
            record["forwarded"] = True
        except (OSError, ProcessLookupError, ValueError) as error:
            record["error"] = str(error)
        forwarded_signals.append(record)

    signal_names = ("SIGINT", "SIGTERM", "SIGHUP", "SIGUSR1")
    if threading.current_thread() is threading.main_thread():
        for name in signal_names:
            signum = getattr(signal, name, None)
            if signum is None:
                continue
            previous_handlers[signum] = signal.getsignal(signum)
            signal.signal(signum, forward_signal)
    timed_out = False
    try:
        stdout_bytes, stderr_bytes = process.communicate(timeout=timeout_secs)
    except subprocess.TimeoutExpired:
        timed_out = True
        try:
            if os.name == "posix":
                os.killpg(process.pid, signal.SIGTERM)
            else:
                process.terminate()
            stdout_bytes, stderr_bytes = process.communicate(timeout=5)
        except (subprocess.TimeoutExpired, OSError):
            if os.name == "posix":
                try:
                    os.killpg(process.pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
            else:
                process.kill()
            stdout_bytes, stderr_bytes = process.communicate()
    finally:
        for signum, previous in previous_handlers.items():
            signal.signal(signum, previous)
    stdout = stdout_bytes.decode("utf-8", errors="replace")
    stderr = stderr_bytes.decode("utf-8", errors="replace")
    run_result = subprocess.CompletedProcess(
        command,
        process.returncode,
        stdout=stdout,
        stderr=stderr,
    )
    ended_utc = _now_utc_iso()
    step = {
        "command": command,
        "started_utc": started_utc,
        "ended_utc": ended_utc,
        "exit_code": run_result.returncode,
        "stdout": run_result.stdout,
        "stderr": run_result.stderr,
        "stdout_bytes": len(stdout_bytes),
        "stderr_bytes": len(stderr_bytes),
        "stdout_sha256": _sha256_prefixed_bytes(stdout_bytes),
        "stderr_sha256": _sha256_prefixed_bytes(stderr_bytes),
        "forwarded_signals": forwarded_signals,
        "timeout_secs": timeout_secs,
        "status": (
            "timed_out"
            if timed_out
            else ("ok" if run_result.returncode == 0 else "command_failed")
        ),
    }
    return run_result, step


def _reference_preflight_record(reference: EnsureReferencePrepared) -> dict[str, Any]:
    return {
        "genome_id": reference.genome_id,
        "catalog_path": reference.catalog_path,
        "cache_dir": reference.cache_dir,
        "status_timeout_secs": reference.status_timeout_secs,
        "prepare_timeout_secs": reference.prepare_timeout_secs,
        "prepared_before": None,
        "prepared_after": None,
        "prepare_attempted": False,
        "status_before": None,
        "status_after": None,
        "steps": [],
        "last_failure": None,
        "status": "pending",
    }


def _parse_json_stdout(
    stdout: str,
    context: str,
    *,
    failure_summary: dict[str, Any] | None = None,
) -> Any:
    try:
        return json.loads(stdout)
    except json.JSONDecodeError as e:
        detail = _build_failure_message(
            headline=f"{context} did not return valid JSON: {e}",
            failure_summary=failure_summary,
        )
        raise SkillError(detail) from e


def _ensure_reference_prepared(
    reference: EnsureReferencePrepared,
    resolution: CliResolution,
    execution_cwd: Path,
    record: dict[str, Any],
) -> None:
    status_args = _build_reference_status_args(reference)
    status_result, status_step = _run_cli_command(
        resolution,
        status_args,
        execution_cwd,
        reference.status_timeout_secs,
    )
    record["steps"].append(status_step)
    if status_result.returncode != 0:
        record["status"] = "failed"
        failure_summary = _build_failure_summary(
            stage="reference_status_preflight",
            step=status_step,
            execution_cwd=execution_cwd,
        )
        record["last_failure"] = failure_summary
        raise SkillError(
            _build_failure_message(
                headline=(
                    "reference-status preflight failed while checking whether the requested "
                    "reference is already prepared."
                ),
                failure_summary=failure_summary,
            )
        )
    status_payload = _parse_json_stdout(
        status_result.stdout,
        f"genomes status {reference.genome_id}",
        failure_summary=_build_failure_summary(
            stage="reference_status_preflight",
            step=status_step,
            execution_cwd=execution_cwd,
        ),
    )
    record["status_before"] = status_payload
    record["prepared_before"] = bool(status_payload.get("prepared"))
    if record["prepared_before"]:
        record["prepared_after"] = True
        record["status_after"] = status_payload
        record["status"] = "already_prepared"
        return

    prepare_args = _build_reference_prepare_args(reference)
    prepare_result, prepare_step = _run_cli_command(
        resolution,
        prepare_args,
        execution_cwd,
        reference.prepare_timeout_secs,
    )
    record["steps"].append(prepare_step)
    record["prepare_attempted"] = True
    if prepare_result.returncode != 0:
        record["status"] = "failed"
        failure_summary = _build_failure_summary(
            stage="reference_prepare_preflight",
            step=prepare_step,
            execution_cwd=execution_cwd,
        )
        record["last_failure"] = failure_summary
        raise SkillError(
            _build_failure_message(
                headline=(
                    "reference-prepare preflight failed while trying to prepare the requested "
                    "reference locally."
                ),
                failure_summary=failure_summary,
            )
        )

    status_after_result, status_after_step = _run_cli_command(
        resolution,
        status_args,
        execution_cwd,
        reference.status_timeout_secs,
    )
    record["steps"].append(status_after_step)
    if status_after_result.returncode != 0:
        record["status"] = "failed"
        failure_summary = _build_failure_summary(
            stage="post_prepare_reference_status_preflight",
            step=status_after_step,
            execution_cwd=execution_cwd,
        )
        record["last_failure"] = failure_summary
        raise SkillError(
            _build_failure_message(
                headline=(
                    "post-prepare reference-status preflight failed after the prepare step "
                    "completed."
                ),
                failure_summary=failure_summary,
            )
        )
    status_after_payload = _parse_json_stdout(
        status_after_result.stdout,
        f"post-prepare genomes status {reference.genome_id}",
        failure_summary=_build_failure_summary(
            stage="post_prepare_reference_status_preflight",
            step=status_after_step,
            execution_cwd=execution_cwd,
        ),
    )
    record["status_after"] = status_after_payload
    record["prepared_after"] = bool(status_after_payload.get("prepared"))
    if not record["prepared_after"]:
        record["status"] = "failed"
        raise SkillError(
            "reference-prepare completed but the requested reference still is not reported as prepared"
        )
    record["status"] = "prepared_during_run"


def _sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _sha256_file_prefixed(path: Path) -> str:
    return "sha256:" + _sha256_file(path)


def _resolve_bound_input_path(
    raw_path: str, execution_cwd: Path, script_path: Path
) -> Path:
    path = Path(raw_path).expanduser()
    candidates = [path] if path.is_absolute() else [execution_cwd / path]
    candidates.extend(_request_path_candidates(raw_path, script_path))
    seen: set[Path] = set()
    for candidate in candidates:
        resolved = candidate.resolve()
        if resolved in seen:
            continue
        seen.add(resolved)
        if resolved.is_file():
            return resolved
    raise SkillError(f"content-bound input file '{raw_path}' does not exist")


def _shell_input_paths(request: Request) -> list[str]:
    tokens: list[str] = []
    if request.shell_line:
        try:
            tokens.extend(shlex.split(request.shell_line))
        except ValueError as e:
            raise SkillError(f"shell_line cannot be tokenized for input binding: {e}") from e
    if request.raw_args:
        tokens.extend(request.raw_args)
    return [token[1:] for token in tokens if token.startswith("@") and len(token) > 1]


def _prepare_content_bound_inputs(
    request: Request,
    *,
    request_source_path: Path | None,
    execution_cwd: Path,
    script_path: Path,
) -> list[dict[str, Any]]:
    specs: list[dict[str, Any]] = []
    if request_source_path is not None:
        specs.append(
            {
                "binding_id": "wrapper_request",
                "path": str(request_source_path),
                "role": "wrapper_request",
                "media_type": "application/json",
                "expected_sha256": None,
            }
        )
    known_paths = (
        ("workflow_path", request.workflow_path, "workflow_definition"),
        ("plan_path", request.plan_path, "agent_plan"),
        ("catalog_path", request.catalog_path, "catalog"),
    )
    for binding_id, raw_path, role in known_paths:
        if raw_path:
            specs.append(
                {
                    "binding_id": binding_id,
                    "path": raw_path,
                    "role": role,
                    "media_type": None,
                    "expected_sha256": None,
                }
            )
    if request.primer3_executable and ("/" in request.primer3_executable or "\\" in request.primer3_executable):
        specs.append(
            {
                "binding_id": "primer3_executable",
                "path": request.primer3_executable,
                "role": "external_executable",
                "media_type": "application/octet-stream",
                "expected_sha256": None,
            }
        )
    if isinstance(request.ensure_reference_prepared, EnsureReferencePrepared):
        if request.ensure_reference_prepared.catalog_path:
            specs.append(
                {
                    "binding_id": "reference_catalog",
                    "path": request.ensure_reference_prepared.catalog_path,
                    "role": "reference_catalog",
                    "media_type": "application/json",
                    "expected_sha256": None,
                }
            )
    for index, raw_path in enumerate(_shell_input_paths(request)):
        specs.append(
            {
                "binding_id": f"shell_at_input_{index + 1}",
                "path": raw_path,
                "role": "shell_at_file",
                "media_type": None,
                "expected_sha256": None,
            }
        )
    specs.extend(request.input_bindings or [])

    by_path: dict[Path, dict[str, Any]] = {}
    for spec in specs:
        resolved = _resolve_bound_input_path(
            str(spec["path"]), execution_cwd, script_path
        )
        digest = _sha256_file_prefixed(resolved)
        expected = spec.get("expected_sha256")
        if expected is not None and expected != digest:
            raise SkillError(
                f"content-bound input '{spec['binding_id']}' SHA-256 mismatch: "
                f"expected {expected}, observed {digest}"
            )
        record = by_path.get(resolved)
        if record is None:
            record = {
                "resolved_path": str(resolved),
                "size_bytes": resolved.stat().st_size,
                "sha256": digest,
                "binding_ids": [],
                "declared_paths": [],
                "roles": [],
                "media_types": [],
                "expected_sha256_values": [],
            }
            by_path[resolved] = record
        record["binding_ids"].append(str(spec["binding_id"]))
        record["declared_paths"].append(str(spec["path"]))
        if spec.get("role"):
            record["roles"].append(str(spec["role"]))
        if spec.get("media_type"):
            record["media_types"].append(str(spec["media_type"]))
        if expected:
            record["expected_sha256_values"].append(expected)
    records = []
    for record in by_path.values():
        for field in (
            "binding_ids",
            "declared_paths",
            "roles",
            "media_types",
            "expected_sha256_values",
        ):
            record[field] = sorted(set(record[field]))
        records.append(record)
    records.sort(key=lambda row: (row["binding_ids"], row["resolved_path"]))
    return records


def _verify_content_bound_inputs_unchanged(
    records: list[dict[str, Any]],
) -> None:
    changed: list[str] = []
    for record in records:
        path = Path(str(record["resolved_path"]))
        post_execution: dict[str, Any] = {
            "exists": path.is_file(),
            "size_bytes": None,
            "sha256": None,
            "status": "missing",
        }
        if path.is_file():
            post_execution["size_bytes"] = path.stat().st_size
            post_execution["sha256"] = _sha256_file_prefixed(path)
            post_execution["status"] = (
                "unchanged"
                if post_execution["size_bytes"] == record["size_bytes"]
                and post_execution["sha256"] == record["sha256"]
                else "changed"
            )
        record["post_execution_verification"] = post_execution
        if post_execution["status"] != "unchanged":
            changed.append(
                f"{','.join(record['binding_ids'])} ({post_execution['status']})"
            )
    if changed:
        raise SkillError(
            "content-bound input changed during execution: " + "; ".join(changed)
        )


def _state_content_snapshot(request: Request | None, execution_cwd: Path) -> dict[str, Any]:
    if request is None or not request.state_path:
        return {"declared": False, "path": None, "exists": False, "sha256": None}
    raw_path = Path(request.state_path).expanduser()
    path = raw_path if raw_path.is_absolute() else execution_cwd / raw_path
    resolved = path.resolve()
    if not resolved.exists():
        return {
            "declared": True,
            "path": str(resolved),
            "exists": False,
            "size_bytes": None,
            "sha256": None,
        }
    if not resolved.is_file():
        raise SkillError(f"state_path '{resolved}' is not a regular file")
    return {
        "declared": True,
        "path": str(resolved),
        "exists": True,
        "size_bytes": resolved.stat().st_size,
        "sha256": _sha256_file_prefixed(resolved),
    }


def _verified_delegation(
    request: Request, script_path: Path
) -> dict[str, Any] | None:
    if request.delegation is None:
        return None
    delegation = request.delegation
    skill_root = script_path.resolve().parent
    skills_root = skill_root.parent
    source_root = (skills_root / delegation["source_skill"]).resolve()
    if source_root.parent != skills_root.resolve():
        raise SkillError("delegation source skill escaped the co-shipped skills root")
    descriptor_path = source_root / "INTENTS.json"
    catalog_path = source_root / "catalog_entry.json"
    if not descriptor_path.is_file() or not catalog_path.is_file():
        raise SkillError(
            "delegated skill metadata is not co-deployed with gentle-cloning: "
            f"{delegation['source_skill']}"
        )
    descriptor = _read_json(descriptor_path)
    catalog = _read_json(catalog_path)
    if descriptor.get("schema") != "clawbio.skill_intents.v1":
        raise SkillError("delegated skill returned an unsupported INTENTS schema")
    if descriptor.get("skill") != delegation["source_skill"]:
        raise SkillError("delegated skill INTENTS identity mismatch")
    if catalog.get("name") != delegation["source_skill"]:
        raise SkillError("delegated skill catalog identity mismatch")
    if catalog.get("version") != delegation["source_skill_version"]:
        raise SkillError("delegated skill catalog version mismatch")
    descriptor_sha256 = _sha256_file_prefixed(descriptor_path)
    catalog_sha256 = _sha256_file_prefixed(catalog_path)
    if delegation.get("descriptor_sha256") not in (None, descriptor_sha256):
        raise SkillError("delegated skill INTENTS SHA-256 mismatch")
    if delegation.get("catalog_sha256") not in (None, catalog_sha256):
        raise SkillError("delegated skill catalog SHA-256 mismatch")

    own_catalog = _read_catalog_entry(script_path)
    if own_catalog.get("name") != SKILL_NAME or own_catalog.get("version") != SKILL_CONTRACT_VERSION:
        raise SkillError("gentle-cloning catalog identity/version disagrees with its runtime contract")
    expected_contract = {
        "skill": SKILL_NAME,
        "skill_version": SKILL_CONTRACT_VERSION,
        "request_schema": REQUEST_SCHEMA,
        "result_schema": RESULT_SCHEMA,
        "execution_manifest_schema": EXECUTION_MANIFEST_SCHEMA,
    }
    if catalog.get("delegate_contract") != expected_contract:
        raise SkillError(
            "delegated skill is incompatible with this gentle-cloning contract"
        )
    routes = descriptor.get("routes")
    if not isinstance(routes, list):
        raise SkillError("delegated skill INTENTS descriptor has no route array")
    route = next(
        (
            candidate
            for candidate in routes
            if isinstance(candidate, dict)
            and candidate.get("intent_id") == delegation["intent_id"]
        ),
        None,
    )
    if route is None:
        raise SkillError("delegation intent_id is absent from the source descriptor")
    requires_confirmation = route.get("requires_confirmation", False)
    if not isinstance(requires_confirmation, bool):
        raise SkillError("delegation route requires_confirmation must be boolean")
    plan = route.get("plan")
    step_index = delegation["plan_step_index"]
    if not isinstance(plan, list) or step_index >= len(plan):
        raise SkillError("delegation plan_step_index is absent from the source intent")
    step = plan[step_index]
    if not isinstance(step, dict) or step.get("kind") != "skill_run" or step.get("skill") != SKILL_NAME:
        raise SkillError("delegation plan step does not target gentle-cloning")
    confirmation = step.get("confirmation")
    if requires_confirmation:
        if not isinstance(confirmation, dict) or confirmation.get("required") is not True:
            raise SkillError(
                "confirmation-gated delegation lacks a matching plan-step confirmation"
            )
        confirmation_reason = _required_handoff_string(
            confirmation.get("reason"), "delegation confirmation reason"
        )
    else:
        if isinstance(confirmation, dict) and confirmation.get("required") is True:
            raise SkillError(
                "delegation route and plan step disagree on confirmation policy"
            )
        confirmation_reason = None
    declared_request: Any = None
    if isinstance(step.get("input_template"), dict):
        declared_request = step["input_template"]
    elif isinstance(step.get("input"), str):
        declared_request_path = (source_root / step["input"]).resolve()
        if not declared_request_path.is_relative_to(source_root) or not declared_request_path.is_file():
            raise SkillError("delegation plan input is missing or escaped its skill root")
        declared_request = _read_json(declared_request_path)
    if not isinstance(declared_request, dict):
        raise SkillError("delegation plan step has no structured gentle-cloning request")
    declared_delegation = declared_request.get("delegation")
    if not isinstance(declared_delegation, dict):
        raise SkillError("delegation plan request lacks its static delegation identity")
    for field in ("source_skill", "source_skill_version", "intent_id", "plan_step_index"):
        if declared_delegation.get(field) != delegation.get(field):
            raise SkillError(f"delegation plan request disagrees on {field}")
    return {
        "schema": DELEGATION_SCHEMA,
        "source_skill": delegation["source_skill"],
        "source_skill_version": delegation["source_skill_version"],
        "intent_id": delegation["intent_id"],
        "plan_step_index": step_index,
        "resolved_slots": delegation.get("resolved_slots"),
        "resolved_slots_sha256": (
            _sha256_prefixed_json(delegation["resolved_slots"])
            if delegation.get("resolved_slots") is not None
            else None
        ),
        "descriptor_path": str(descriptor_path),
        "descriptor_sha256": descriptor_sha256,
        "catalog_path": str(catalog_path),
        "catalog_sha256": catalog_sha256,
        "route_sha256": _sha256_prefixed_json(route),
        "plan_step_sha256": _sha256_prefixed_json(step),
        "delegate_contract": expected_contract,
        "requires_confirmation": requires_confirmation,
        "confirmation_reason": confirmation_reason,
        "execution_policy": (
            "proposal_then_approved_execution"
            if requires_confirmation
            else "automatic_read_only"
        ),
        "routing_scope": "selected_route_record_only_not_natural_language_reproduction",
        "compatibility_status": "verified",
    }


def _load_approved_execution(
    payload: dict[str, Any], envelope_path: Path, script_path: Path
) -> ApprovedExecution:
    _reject_unknown_fields(
        payload,
        {"schema", "proposal_path", "approval"},
        "approved execution request",
    )
    if payload.get("schema") != APPROVED_EXECUTION_REQUEST_SCHEMA:
        raise SkillError(
            f"approved execution request schema must be '{APPROVED_EXECUTION_REQUEST_SCHEMA}'"
        )
    raw_proposal_path = _required_handoff_string(
        payload.get("proposal_path"), "approved execution proposal_path"
    )
    proposal_path = _resolve_existing_request_file(raw_proposal_path, script_path)
    proposal = _read_json(proposal_path)
    _reject_unknown_fields(
        proposal,
        {
            "schema",
            "created_utc",
            "proposal_digest",
            "proposal_digest_scope",
            "approval_basis",
            "review",
        },
        "execution proposal",
    )
    if proposal.get("schema") != EXECUTION_PROPOSAL_SCHEMA:
        raise SkillError(
            f"execution proposal schema must be '{EXECUTION_PROPOSAL_SCHEMA}'"
        )
    approval_basis = proposal.get("approval_basis")
    if not isinstance(approval_basis, dict):
        raise SkillError("execution proposal approval_basis must be an object")
    observed_digest = _sha256_prefixed_json(approval_basis)
    declared_digest = _normalise_expected_sha256(
        proposal.get("proposal_digest"), "execution proposal proposal_digest"
    )
    if declared_digest != observed_digest:
        raise SkillError(
            "execution proposal digest mismatch; the stored proposal was altered"
        )

    approval = payload.get("approval")
    if not isinstance(approval, dict):
        raise SkillError("approved execution request approval must be an object")
    _reject_unknown_fields(
        approval,
        {
            "schema",
            "approval_scope",
            "proposal_digest",
            "approval_id",
            "approved_by",
            "approved_at_utc",
        },
        "execution approval",
    )
    if approval.get("schema") != EXECUTION_APPROVAL_SCHEMA:
        raise SkillError(
            f"execution approval schema must be '{EXECUTION_APPROVAL_SCHEMA}'"
        )
    if approval.get("approval_scope") != "execute_exact_proposal":
        raise SkillError(
            "execution approval approval_scope must be 'execute_exact_proposal'"
        )
    approved_digest = _normalise_expected_sha256(
        approval.get("proposal_digest"), "execution approval proposal_digest"
    )
    if approved_digest != observed_digest:
        raise SkillError(
            "execution approval refers to a different proposal digest"
        )
    _verify_execution_proposal_review(proposal)
    for field_name in ("approval_id", "approved_by"):
        value = _required_handoff_string(
            approval.get(field_name), f"execution approval {field_name}"
        )
        if len(value) > 256:
            raise SkillError(f"execution approval {field_name} exceeds 256 characters")
        approval[field_name] = value
    approved_at_utc = _optional_handoff_string(
        approval.get("approved_at_utc"), "execution approval approved_at_utc"
    )
    approval["approved_at_utc"] = approved_at_utc
    approval["proposal_digest"] = approved_digest
    return ApprovedExecution(
        proposal_path=proposal_path.resolve(),
        proposal=proposal,
        approval=approval,
        envelope_path=envelope_path.resolve(),
    )


def _verify_execution_proposal_review(proposal: dict[str, Any]) -> None:
    basis = proposal["approval_basis"]
    review = proposal.get("review")
    if not isinstance(review, dict):
        raise SkillError("execution proposal review must be an object")
    _reject_unknown_fields(
        review,
        {
            "status",
            "source_request",
            "source_skill",
            "intent_id",
            "confirmation_reason",
            "exact_command",
            "exact_command_text",
            "input_assumptions",
            "selection_context",
            "resolved_paths",
            "backend_resolution",
            "execution_environment",
            "wrapper_contract",
            "trust_boundary",
            "execution_instruction",
        },
        "execution proposal review",
    )
    delegation = basis.get("delegation")
    normalized_request = basis.get("normalized_request")
    execution = basis.get("execution")
    if not isinstance(delegation, dict) or not isinstance(normalized_request, dict):
        raise SkillError("execution proposal review source records are incomplete")
    if not isinstance(execution, dict):
        raise SkillError("execution proposal execution binding is incomplete")
    expected = {
        "status": "approval_required",
        "source_request": basis.get("source_request"),
        "source_skill": delegation.get("source_skill"),
        "intent_id": delegation.get("intent_id"),
        "confirmation_reason": delegation.get("confirmation_reason"),
        "exact_command": execution.get("command"),
        "exact_command_text": _format_command_text(execution.get("command")),
        "input_assumptions": list(normalized_request.get("input_claims") or []),
        "selection_context": _proposal_selection_context(delegation),
        "resolved_paths": execution.get("resolved_paths"),
        "backend_resolution": basis.get("backend_resolution"),
        "execution_environment": basis.get("execution_environment"),
        "wrapper_contract": basis.get("wrapper_contract"),
        "trust_boundary": (
            "GENtle verifies the proposal and approval digest; the caller control "
            "plane is responsible for obtaining and attesting human approval."
        ),
    }
    for field_name, expected_value in expected.items():
        if review.get(field_name) != expected_value:
            raise SkillError(
                f"execution proposal review projection mismatch: {field_name}"
            )
    instruction = review.get("execution_instruction")
    if instruction != _proposal_execution_instruction(proposal["proposal_digest"]):
        raise SkillError("execution proposal review instruction is inconsistent")


def _proposal_execution_instruction(proposal_digest: str) -> dict[str, Any]:
    return {
        "schema": APPROVED_EXECUTION_REQUEST_SCHEMA,
        "proposal_path": "ABSOLUTE_PATH_TO_THIS_PROPOSAL",
        "approval": {
            "schema": EXECUTION_APPROVAL_SCHEMA,
            "approval_scope": "execute_exact_proposal",
            "proposal_digest": proposal_digest,
            "approval_id": "CALLER_ISSUED_APPROVAL_ID",
            "approved_by": "CALLER_ATTESTED_APPROVER",
            "approved_at_utc": None,
        },
    }


def _proposal_content_bindings(
    records: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    proposal_records: list[dict[str, Any]] = []
    for record in records:
        if "wrapper_request" in record.get("binding_ids", []):
            continue
        if "approval_control" in record.get("roles", []):
            continue
        proposal_records.append(
            {
                key: value
                for key, value in record.items()
                if key != "post_execution_verification"
            }
        )
    return proposal_records


def _approval_control_binding(
    path: Path, binding_id: str
) -> dict[str, Any]:
    resolved = path.resolve()
    return {
        "resolved_path": str(resolved),
        "size_bytes": resolved.stat().st_size,
        "sha256": _sha256_file_prefixed(resolved),
        "binding_ids": [binding_id],
        "declared_paths": [str(path)],
        "roles": ["approval_control"],
        "media_types": ["application/json"],
        "expected_sha256_values": [],
    }


def _resolved_execution_paths(
    request: Request, execution_cwd: Path
) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []

    def add(kind: str, declared: str) -> None:
        raw = Path(declared).expanduser()
        resolved = raw if raw.is_absolute() else execution_cwd / raw
        records.append(
            {
                "kind": kind,
                "declared_path": declared,
                "resolved_path": str(resolved.resolve()),
            }
        )

    if request.state_path:
        add("project_state", request.state_path)
    for path in request.expected_artifacts or []:
        add("expected_artifact", path)
    for field_name in (
        "workflow_path",
        "plan_path",
        "catalog_path",
        "output_path",
        "output_prefix",
        "svg_path",
        "render_svg_path",
        "product_output_prefix",
        "product_gel_svg_path",
    ):
        value = getattr(request, field_name)
        if isinstance(value, str) and value:
            add(field_name, value)

    shell_tokens: list[str] = []
    if request.shell_line:
        try:
            shell_tokens = shlex.split(request.shell_line)
        except ValueError as e:
            raise SkillError(f"shell_line cannot be tokenized for path review: {e}") from e
    path_flags = {
        "--path",
        "--output",
        "--output-dir",
        "--order-table",
        "--gel-svg",
        "--checkpoint-path",
        "--primer3-exec",
    }
    for index, token in enumerate(shell_tokens[:-1]):
        if token in path_flags:
            add(f"shell_option:{token}", shell_tokens[index + 1])
    for token in shell_tokens:
        if token.startswith("@") and len(token) > 1:
            add("shell_at_input", token[1:])

    resolved_slots = None
    if isinstance(request.delegation, dict):
        resolved_slots = request.delegation.get("resolved_slots")
    if isinstance(resolved_slots, dict):
        for key, value in sorted(resolved_slots.items()):
            if not isinstance(value, str) or not value:
                continue
            if key.endswith(("_path", "_dir", "_prefix")):
                add(f"resolved_slot:{key}", value)

    deduplicated = {
        (record["kind"], record["declared_path"], record["resolved_path"]): record
        for record in records
    }
    return [deduplicated[key] for key in sorted(deduplicated)]


def _assert_delegated_request_is_resolved(request: Request) -> None:
    payload_text = json.dumps(dataclasses.asdict(request), ensure_ascii=True)
    placeholders = sorted(set(re.findall(r"\{[A-Za-z_][A-Za-z0-9_]*\}", payload_text)))
    if placeholders:
        raise SkillError(
            "delegated request still contains unresolved slot placeholder(s): "
            + ", ".join(placeholders)
        )


def _pin_delegated_primer_backend(
    request: Request,
    resolution: CliResolution,
    execution_cwd: Path,
) -> tuple[dict[str, Any] | None, dict[str, Any] | None]:
    if request.mode != "shell" or not request.shell_line:
        return None, None
    try:
        tokens = shlex.split(request.shell_line)
    except ValueError as e:
        raise SkillError(f"shell_line cannot be tokenized for backend pinning: {e}") from e
    if "--backend" not in tokens:
        return None, None
    backend_index = tokens.index("--backend")
    if backend_index + 1 >= len(tokens):
        raise SkillError("delegated --backend option has no value")
    requested_backend = tokens[backend_index + 1].strip().lower()
    if requested_backend not in {"auto", "internal", "primer3"}:
        raise SkillError(
            "delegated primer backend must be one of auto|internal|primer3"
        )
    primer3_executable: str | None = None
    if "--primer3-exec" in tokens:
        executable_index = tokens.index("--primer3-exec")
        if executable_index + 1 >= len(tokens):
            raise SkillError("delegated --primer3-exec option has no value")
        primer3_executable = tokens[executable_index + 1]

    if requested_backend == "internal":
        return (
            {
                "requested_backend": "internal",
                "selected_backend": "internal",
                "selection_basis": "explicit_request",
                "selected_backend_version": "GENtle runtime version bound separately",
                "primer3_preflight": None,
                "primer3_executable": None,
            },
            None,
        )

    preflight_tokens = ["primers", "preflight", "--backend", requested_backend]
    if primer3_executable:
        preflight_tokens.extend(["--primer3-exec", primer3_executable])
    run_result, step = _run_cli_command(
        resolution,
        ["shell", shlex.join(preflight_tokens)],
        execution_cwd,
        min(request.timeout_secs, 180),
    )
    if run_result.returncode != 0:
        raise SkillError("delegated primer backend preflight failed")
    payload = _parse_json_stdout(
        run_result.stdout, "delegated primer backend preflight"
    )
    preflight = payload.get("preflight") if isinstance(payload, dict) else None
    if not isinstance(preflight, dict):
        raise SkillError("delegated primer backend preflight returned no preflight record")
    preflight = {
        key: value for key, value in preflight.items() if key != "probe_time_ms"
    }
    primer3_ready = bool(preflight.get("version_probe_ok")) and isinstance(
        preflight.get("resolved_path"), str
    )
    if requested_backend == "primer3" and not primer3_ready:
        raise SkillError(
            "explicit Primer3 backend is not version-probed and path-resolved; "
            "review the preflight before proposing execution"
        )
    selected_backend = "primer3" if primer3_ready else "internal"
    tokens[backend_index + 1] = selected_backend
    selected_executable: dict[str, Any] | None = None
    if selected_backend == "primer3":
        resolved_executable = Path(str(preflight["resolved_path"])).resolve()
        if not resolved_executable.is_file():
            raise SkillError(
                "Primer3 preflight resolved path is not a regular file"
            )
        if "--primer3-exec" in tokens:
            executable_index = tokens.index("--primer3-exec")
            tokens[executable_index + 1] = str(resolved_executable)
        else:
            tokens.extend(["--primer3-exec", str(resolved_executable)])
        selected_executable = {
            "path": str(resolved_executable),
            "size_bytes": resolved_executable.stat().st_size,
            "sha256": _sha256_file_prefixed(resolved_executable),
        }
    request.shell_line = shlex.join(tokens)
    return (
        {
            "requested_backend": requested_backend,
            "selected_backend": selected_backend,
            "selection_basis": (
                "primer3_version_probe" if primer3_ready else "deterministic_internal_fallback"
            ),
            "selected_backend_version": (
                preflight.get("version")
                if selected_backend == "primer3"
                else "GENtle runtime version bound separately"
            ),
            "primer3_preflight": preflight,
            "primer3_executable": selected_executable,
        },
        step,
    )


def _proposal_selection_context(delegation: dict[str, Any]) -> dict[str, Any]:
    slots = delegation.get("resolved_slots")
    if not isinstance(slots, dict):
        return {"selected_material": {}, "specificity_spaces": {}}
    selected_keys = {
        "pair_id",
        "pair_rank",
        "report_id",
        "panel_report_id",
        "seq_id",
        "seq_ids",
        "feature_id",
        "protocol_id",
    }
    selected = {
        key: slots[key]
        for key in sorted(slots)
        if key in selected_keys or key.endswith("_pair_id")
    }
    specificity = {
        key: slots[key]
        for key in sorted(slots)
        if "reference" in key or key.endswith("_genome_id")
    }
    return {
        "selected_material": selected,
        "specificity_spaces": specificity,
    }


def _proposal_request_payload(request: Request) -> dict[str, Any]:
    return {"schema": REQUEST_SCHEMA, **dataclasses.asdict(request)}


def _build_execution_proposal(
    *,
    request: Request,
    request_source_path: Path | None,
    content_bound_inputs: list[dict[str, Any]],
    delegation: dict[str, Any],
    resolution: CliResolution,
    gentle_version: str,
    execution_cwd: Path,
    state_before: dict[str, Any],
    backend_resolution: dict[str, Any] | None,
) -> dict[str, Any]:
    request_payload = _proposal_request_payload(request)
    command = resolution.argv_prefix + _build_cli_args(request, Path(__file__))
    runtime_files = _runtime_file_bindings(resolution, execution_cwd)
    source_record = None
    if request_source_path is not None:
        source_record = {
            "path": str(request_source_path.resolve()),
            "sha256": _sha256_file_prefixed(request_source_path.resolve()),
        }
    approval_basis = {
        "normalized_request": request_payload,
        "normalized_request_sha256": _sha256_prefixed_json(request_payload),
        "source_request": source_record,
        "delegation": delegation,
        "execution": {
            "command": command,
            "execution_cwd": str(execution_cwd.resolve()),
            "resolved_paths": _resolved_execution_paths(request, execution_cwd),
        },
        "runtime": {
            "gentle_version": gentle_version,
            "resolver": dataclasses.asdict(resolution),
            "content_bound_runtime_files": runtime_files,
        },
        "execution_environment": _approval_environment_binding(execution_cwd),
        "wrapper_contract": _wrapper_contract_binding(Path(__file__)),
        "state_before": state_before,
        "content_bound_inputs": _proposal_content_bindings(content_bound_inputs),
        "backend_resolution": backend_resolution,
        "reference_request": (
            dataclasses.asdict(request.ensure_reference_prepared)
            if isinstance(request.ensure_reference_prepared, EnsureReferencePrepared)
            else None
        ),
    }
    proposal_digest = _sha256_prefixed_json(approval_basis)
    return {
        "schema": EXECUTION_PROPOSAL_SCHEMA,
        "created_utc": _now_utc_iso(),
        "proposal_digest": proposal_digest,
        "proposal_digest_scope": (
            "canonical approval_basis; timestamps and proposal file location excluded"
        ),
        "approval_basis": approval_basis,
        "review": {
            "status": "approval_required",
            "source_request": source_record,
            "source_skill": delegation.get("source_skill"),
            "intent_id": delegation.get("intent_id"),
            "confirmation_reason": delegation.get("confirmation_reason"),
            "exact_command": command,
            "exact_command_text": _format_command_text(command),
            "input_assumptions": list(request.input_claims or []),
            "selection_context": _proposal_selection_context(delegation),
            "resolved_paths": approval_basis["execution"]["resolved_paths"],
            "backend_resolution": backend_resolution,
            "execution_environment": approval_basis["execution_environment"],
            "wrapper_contract": approval_basis["wrapper_contract"],
            "trust_boundary": (
                "GENtle verifies the proposal and approval digest; the caller control "
                "plane is responsible for obtaining and attesting human approval."
            ),
            "execution_instruction": _proposal_execution_instruction(
                proposal_digest
            ),
        },
    }


def _verify_backend_binding(backend_resolution: dict[str, Any] | None) -> None:
    if not isinstance(backend_resolution, dict):
        return
    executable = backend_resolution.get("primer3_executable")
    if not isinstance(executable, dict):
        return
    path = Path(str(executable.get("path", "")))
    if not path.is_file():
        raise SkillError("approved Primer3 executable is missing")
    observed = {
        "path": str(path.resolve()),
        "size_bytes": path.stat().st_size,
        "sha256": _sha256_file_prefixed(path.resolve()),
    }
    if observed != executable:
        raise SkillError("approved Primer3 executable identity changed")


def _verify_approved_execution_context(
    *,
    approved: ApprovedExecution,
    request: Request,
    content_bound_inputs: list[dict[str, Any]],
    delegation: dict[str, Any],
    resolution: CliResolution,
    gentle_version: str | None,
    execution_cwd: Path,
    state_before: dict[str, Any],
) -> None:
    basis = approved.proposal["approval_basis"]
    request_payload = _proposal_request_payload(request)
    checks = {
        "normalized request": (
            basis.get("normalized_request"),
            request_payload,
        ),
        "normalized request digest": (
            basis.get("normalized_request_sha256"),
            _sha256_prefixed_json(request_payload),
        ),
        "delegated route metadata": (basis.get("delegation"), delegation),
        "exact command": (
            (basis.get("execution") or {}).get("command"),
            resolution.argv_prefix + _build_cli_args(request, Path(__file__)),
        ),
        "execution working directory": (
            (basis.get("execution") or {}).get("execution_cwd"),
            str(execution_cwd.resolve()),
        ),
        "resolved execution paths": (
            (basis.get("execution") or {}).get("resolved_paths"),
            _resolved_execution_paths(request, execution_cwd),
        ),
        "runtime resolver": (
            (basis.get("runtime") or {}).get("resolver"),
            dataclasses.asdict(resolution),
        ),
        "runtime file identity": (
            (basis.get("runtime") or {}).get("content_bound_runtime_files"),
            _runtime_file_bindings(resolution, execution_cwd),
        ),
        "execution environment": (
            basis.get("execution_environment"),
            _approval_environment_binding(execution_cwd),
        ),
        "wrapper contract": (
            basis.get("wrapper_contract"),
            _wrapper_contract_binding(Path(__file__)),
        ),
        "project state": (basis.get("state_before"), state_before),
        "content-bound scientific inputs": (
            basis.get("content_bound_inputs"),
            _proposal_content_bindings(content_bound_inputs),
        ),
    }
    if gentle_version is not None:
        checks["GENtle runtime version"] = (
            (basis.get("runtime") or {}).get("gentle_version"),
            gentle_version,
        )
    for label, (expected, observed) in checks.items():
        if expected != observed:
            raise SkillError(
                f"approved execution context drifted after approval: {label} changed"
            )
    _verify_backend_binding(basis.get("backend_resolution"))


def _execution_step_manifest(step: dict[str, Any]) -> dict[str, Any]:
    stdout = str(step.get("stdout") or "")
    stderr = str(step.get("stderr") or "")
    return {
        "stage": step.get("manifest_stage"),
        "command": list(step.get("command") or []),
        "started_utc": step.get("started_utc"),
        "ended_utc": step.get("ended_utc"),
        "exit_code": step.get("exit_code"),
        "status": step.get("status"),
        "stdout_sha256": _sha256_prefixed_bytes(stdout.encode("utf-8")),
        "stderr_sha256": _sha256_prefixed_bytes(stderr.encode("utf-8")),
    }


def _json_pointer_token(value: str) -> str:
    return value.replace("~", "~0").replace("/", "~1")


def _native_identifier_bindings(value: Any) -> list[dict[str, Any]]:
    identifier_keys = {
        "report_id",
        "op_id",
        "run_id",
        "assay_test_id",
        "panel_id",
        "collection_id",
    }
    bindings: list[dict[str, Any]] = []

    def visit(candidate: Any, pointer: str) -> None:
        if isinstance(candidate, dict):
            for key in sorted(candidate):
                child_pointer = pointer + "/" + _json_pointer_token(str(key))
                child = candidate[key]
                if key in identifier_keys and isinstance(child, (str, int)):
                    bindings.append(
                        {"json_pointer": child_pointer, "field": key, "value": child}
                    )
                visit(child, child_pointer)
        elif isinstance(candidate, list):
            for index, child in enumerate(candidate):
                visit(child, f"{pointer}/{index}")

    visit(value, "/stdout_json")
    return bindings


def _native_status_bindings(value: Any) -> list[dict[str, Any]]:
    status_keys = {
        "status",
        "overall_status",
        "analysis_completeness",
        "confirmation_status",
    }
    bindings: list[dict[str, Any]] = []

    def visit(candidate: Any, pointer: str) -> None:
        if isinstance(candidate, dict):
            for key in sorted(candidate):
                child_pointer = pointer + "/" + _json_pointer_token(str(key))
                child = candidate[key]
                if key in status_keys and isinstance(child, (str, bool, int, float)):
                    bindings.append(
                        {"json_pointer": child_pointer, "field": key, "value": child}
                    )
                visit(child, child_pointer)
        elif isinstance(candidate, list):
            for index, child in enumerate(candidate):
                visit(child, f"{pointer}/{index}")

    visit(value, "/stdout_json")
    return bindings


def _content_artifact_record(path: Path, role: str, output_dir: Path) -> dict[str, Any]:
    resolved = path.resolve()
    try:
        bundle_path = resolved.relative_to(output_dir.resolve()).as_posix()
    except ValueError:
        bundle_path = None
    return {
        "role": role,
        "path": str(resolved),
        "bundle_path": bundle_path,
        "size_bytes": resolved.stat().st_size,
        "sha256": _sha256_file_prefixed(resolved),
    }


def _runtime_file_bindings(
    resolution: CliResolution | None, execution_cwd: Path
) -> list[dict[str, Any]]:
    if resolution is None:
        return []
    candidates: list[tuple[str, Path]] = []
    for index, token in enumerate(resolution.argv_prefix):
        if token.startswith("-"):
            continue
        raw = Path(token).expanduser()
        resolved: Path | None = None
        if raw.is_absolute() or "/" in token or "\\" in token:
            candidate = raw if raw.is_absolute() else execution_cwd / raw
            if candidate.is_file():
                resolved = candidate.resolve()
        elif index == 0:
            located = shutil.which(token)
            if located:
                resolved = Path(located).resolve()
        if resolved is not None:
            candidates.append((f"argv_prefix_{index}", resolved))
    if resolution.argv_prefix[:2] == ["cargo", "run"]:
        local_binary = execution_cwd / "target" / "debug" / "gentle_cli"
        if local_binary.is_file():
            candidates.append(("cargo_run_output", local_binary.resolve()))
    records: list[dict[str, Any]] = []
    seen: set[Path] = set()
    for binding_id, path in candidates:
        if path in seen:
            continue
        seen.add(path)
        records.append(
            {
                "binding_id": binding_id,
                "path": str(path),
                "size_bytes": path.stat().st_size,
                "sha256": _sha256_file_prefixed(path),
            }
        )
    return records


def _execution_outcome(status: str) -> str:
    if status in {"ok", "degraded_demo"}:
        return "completed"
    if status == "approval_required":
        return "not_executed"
    if status == "incomplete":
        return "incomplete"
    return "failed"


def _build_execution_manifest(
    *,
    request: Request | None,
    request_source_path: Path | None,
    content_bound_inputs: list[dict[str, Any]],
    delegation: dict[str, Any] | None,
    execution_approval: dict[str, Any] | None,
    script_path: Path,
    resolution: CliResolution | None,
    gentle_version: str | None,
    execution_cwd: Path,
    started_utc: str,
    ended_utc: str,
    status: str,
    error_message: str | None,
    state_before: dict[str, Any],
    state_after: dict[str, Any],
    reference_preflight: dict[str, Any] | None,
    execution_steps: list[dict[str, Any]],
    stdout: str,
    stderr: str,
    stdout_json: Any,
    artifacts: list[dict[str, Any]],
    claim_ledger_path: Path | None,
) -> dict[str, Any]:
    request_payload = dataclasses.asdict(request) if request is not None else None
    wrapper_identity = _wrapper_contract_binding(script_path) | {
        "request_schema": REQUEST_SCHEMA,
        "result_schema": RESULT_SCHEMA,
    }
    output_artifacts = list(artifacts)
    if claim_ledger_path is not None and claim_ledger_path.is_file():
        output_artifacts.append(
            _content_artifact_record(claim_ledger_path, "claim_ledger", claim_ledger_path.parent)
        )
    request_source_sha256 = next(
        (
            record["sha256"]
            for record in content_bound_inputs
            if "wrapper_request" in record.get("binding_ids", [])
        ),
        None,
    )
    native_schema = stdout_json.get("schema") if isinstance(stdout_json, dict) else None
    return {
        "schema": EXECUTION_MANIFEST_SCHEMA,
        "provenance_scope": "topic-neutral provider execution receipt",
        "started_utc": started_utc,
        "ended_utc": ended_utc,
        "execution_outcome": _execution_outcome(status),
        "wrapper_status": status,
        "error": error_message,
        "request_binding": {
            "normalized_request_sha256": _sha256_prefixed_json(request_payload),
            "request_source_path": str(request_source_path) if request_source_path else None,
            "request_source_sha256": request_source_sha256,
            "content_bound_inputs": content_bound_inputs,
        },
        "delegation": delegation,
        "execution_approval": execution_approval,
        "wrapper": wrapper_identity,
        "runtime": {
            "gentle_version": gentle_version,
            "resolver": dataclasses.asdict(resolution) if resolution else None,
            "execution_cwd": str(execution_cwd),
            "content_bound_runtime_files": _runtime_file_bindings(
                resolution, execution_cwd
            ),
        },
        "state_binding": {
            "before": state_before,
            "after": state_after,
        },
        "reference_preflight": {
            "payload": reference_preflight,
            "payload_sha256": (
                _sha256_prefixed_json(reference_preflight)
                if reference_preflight is not None
                else None
            ),
        },
        "execution_steps": [_execution_step_manifest(step) for step in execution_steps],
        "native_result": {
            "schema": native_schema,
            "stdout_sha256": _sha256_prefixed_bytes(stdout.encode("utf-8")),
            "stderr_sha256": _sha256_prefixed_bytes(stderr.encode("utf-8")),
            "stdout_json_sha256": (
                _sha256_prefixed_json(stdout_json) if stdout_json is not None else None
            ),
            "identifier_bindings": _native_identifier_bindings(stdout_json),
            "reported_status_bindings": _native_status_bindings(stdout_json),
            "status_interpretation": "native_fields_only_no_wrapper_scientific_reclassification",
        },
        "artifacts": sorted(
            output_artifacts,
            key=lambda artifact: (str(artifact.get("role")), str(artifact.get("path"))),
        ),
    }


def _write_repro_environment(path: Path) -> None:
    lines = [
        "name: gentle-clawbio-skill",
        "channels:",
        "  - defaults",
        "dependencies:",
        f"  - python={sys.version_info.major}.{sys.version_info.minor}",
        "variables:",
        f'  PYTHON_EXECUTABLE: "{sys.executable}"',
        f'  PLATFORM: "{platform.platform()}"',
        f'  PYTHON_VERSION: "{platform.python_version()}"',
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _parse_stdout_json(stdout: str) -> Any | None:
    trimmed = stdout.strip()
    if not trimmed:
        return None
    try:
        return json.loads(trimmed)
    except json.JSONDecodeError:
        return None


def _external_primer_handoff_base_result(
    context: dict[str, Any],
    *,
    status: str,
    gentle_version: str | None,
    state_sha256_after: str | None,
) -> dict[str, Any]:
    target = context["request"]["target"]
    return {
        "schema": EXTERNAL_PRIMER_HANDOFF_RESULT_SCHEMA,
        "status": status,
        "operation": "primers import-external-pairs",
        "collection_id": context["request"]["collection_id"],
        "canonical_request_sha256": context["canonical_request_sha256"],
        "request_snapshot_path": context["request_snapshot_path"],
        "request_snapshot_sha256": context["request_snapshot_sha256"],
        "runtime": {
            "gentle_version": gentle_version,
        },
        "target_state": {
            "state_path": context["state_path"],
            "state_sha256_before": context["state_sha256_before"],
            "state_sha256_after": state_sha256_after,
            "expected_state_sha256": context["expected_state_sha256"],
            "expected_state_match": (
                context["expected_state_sha256"] is None
                or context["expected_state_sha256"]
                == context["state_sha256_before"]
            ),
            "binding_scope": "entire_explicit_state_file",
            "seq_id": target["seq_id"],
            "source_feature_id": target["source_feature_id"],
            "transcript_id": target["transcript_id"],
            "reference_label": target["reference_label"],
            "reference_release": target["reference_release"],
            "reference_declaration_status": (
                "caller_declared_and_bound_to_state_sha256"
                if target["reference_label"] or target["reference_release"]
                else "not_declared"
            ),
        },
        "submitted_record_count": len(context["submitted_records"]),
        "not_submitted_record_count": len(context["not_submitted_records"]),
        "submitted_record_joins": [],
        "not_submitted_records": context["not_submitted_records"],
        "native_result_locator": "stdout_json",
        "report_binding": None,
        "scientific_artifacts": [],
        "warnings": [],
    }


def _external_primer_handoff_not_run_result(
    context: dict[str, Any]
) -> dict[str, Any]:
    result = _external_primer_handoff_base_result(
        context,
        status="not_run",
        gentle_version=None,
        state_sha256_after=context["state_sha256_before"],
    )
    result["execution_binding_sha256"] = _sha256_prefixed_json(
        {
            "canonical_request_sha256": context["canonical_request_sha256"],
            "state_sha256_before": context["state_sha256_before"],
            "submitted_record_count": 0,
        }
    )
    result["analysis_completeness"] = "no_pcr_compatible_records"
    result["warnings"] = [
        "No qPCR or endpoint-PCR pair was submitted to GENtle; cloning and "
        "sequencing oligos remain typed not_submitted records."
    ]
    return result


def _external_primer_handoff_failure_result(
    context: dict[str, Any], gentle_version: str | None, diagnostic: str
) -> dict[str, Any]:
    state_path = Path(context["state_path"])
    state_sha256_after = (
        _sha256_file_prefixed(state_path) if state_path.is_file() else None
    )
    result = _external_primer_handoff_base_result(
        context,
        status="verification_failed",
        gentle_version=gentle_version,
        state_sha256_after=state_sha256_after,
    )
    result["analysis_completeness"] = "unverified"
    result["warnings"] = [diagnostic]
    return result


def _external_primer_handoff_artifact_record(
    *,
    path: Path,
    artifact_kind: str,
    pair_id: str | None = None,
    media_type: str,
) -> dict[str, Any]:
    record: dict[str, Any] = {
        "artifact_kind": artifact_kind,
        "path": str(path),
        "media_type": media_type,
        "size_bytes": path.stat().st_size,
        "sha256": _sha256_file_prefixed(path),
    }
    if pair_id is not None:
        record["pair_id"] = pair_id
    return record


def _verify_external_primer_handoff(
    context: dict[str, Any], stdout_json: Any, gentle_version: str
) -> dict[str, Any]:
    if not isinstance(stdout_json, dict):
        raise SkillError("external-primer handoff returned no JSON object")
    if stdout_json.get("schema") != EXTERNAL_PRIMER_PAIR_IMPORT_COMMAND_SCHEMA:
        raise SkillError(
            "external-primer handoff returned unexpected command schema "
            f"'{stdout_json.get('schema')}'"
        )
    report = stdout_json.get("report")
    if not isinstance(report, dict):
        raise SkillError("external-primer handoff command returned no report object")
    if report.get("schema") != EXTERNAL_PRIMER_PAIR_IMPORT_REPORT_SCHEMA:
        raise SkillError(
            "external-primer handoff returned unexpected report schema "
            f"'{report.get('schema')}'"
        )
    if stdout_json.get("path") != context["report_path"]:
        raise SkillError("external-primer handoff command report path mismatch")
    report_id = report.get("report_id")
    if not isinstance(report_id, str) or not report_id.startswith(
        "external_primer_import_"
    ):
        raise SkillError(
            "external-primer handoff did not return an automatically derived report identity"
        )
    if "--report-id" in str(context.get("shell_line") or ""):
        raise SkillError("external-primer handoff command unexpectedly supplied --report-id")

    target = context["request"]["target"]
    if report.get("batch_id") != context["request"]["collection_id"]:
        raise SkillError("external-primer handoff report batch_id mismatch")
    if report.get("seq_id") != target["seq_id"]:
        raise SkillError("external-primer handoff report seq_id mismatch")
    if report.get("source_feature_id") != target["source_feature_id"]:
        raise SkillError("external-primer handoff report source_feature_id mismatch")
    if report.get("source_record_count") != len(context["submitted_records"]):
        raise SkillError("external-primer handoff report source-record count mismatch")

    provenance = report.get("input_provenance")
    if not isinstance(provenance, dict):
        raise SkillError("external-primer handoff report lacks input provenance")
    if provenance.get("input_format") != "json":
        raise SkillError("external-primer handoff report input format is not JSON")
    if provenance.get("source_path") != context["batch_path"]:
        raise SkillError("external-primer handoff report source path mismatch")
    if provenance.get("source_sha256") != context["batch_file_sha256"]:
        raise SkillError("external-primer handoff report source SHA-256 mismatch")

    pairs = report.get("pairs")
    if not isinstance(pairs, list) or not pairs:
        raise SkillError("external-primer handoff report contains no computed pair rows")
    expected_unique_pair_count = len(
        {
            (
                record["forward_sequence_5_to_3"],
                record["reverse_sequence_5_to_3"],
            )
            for record in context["submitted_records"]
        }
    )
    if report.get("unique_pair_count") != expected_unique_pair_count:
        raise SkillError("external-primer handoff report unique-pair count mismatch")
    if report.get("duplicate_source_record_count") != (
        len(context["submitted_records"]) - expected_unique_pair_count
    ):
        raise SkillError("external-primer handoff report duplicate count mismatch")

    expected_by_record_id = {
        record["record_id"]: record for record in context["submitted_records"]
    }
    joins: list[dict[str, Any]] = []
    digest_rows: list[list[str]] = []
    seen_record_ids: set[str] = set()
    seen_source_record_ids: set[str] = set()
    specificity_status_counts: dict[str, int] = {}
    artifact_root = Path(context["artifact_output_dir"]).resolve()
    scientific_artifacts: list[dict[str, Any]] = []

    for pair in pairs:
        if not isinstance(pair, dict):
            raise SkillError("external-primer handoff report contains a non-object pair row")
        pair_id = pair.get("pair_id")
        if not isinstance(pair_id, str) or not pair_id:
            raise SkillError("external-primer handoff pair row lacks pair_id")
        forward = pair.get("forward")
        reverse = pair.get("reverse")
        if not isinstance(forward, dict) or not isinstance(reverse, dict):
            raise SkillError(f"external-primer handoff pair {pair_id} lacks oligo rows")
        forward_sequence = forward.get("sequence_5_to_3")
        reverse_sequence = reverse.get("sequence_5_to_3")
        if not isinstance(forward_sequence, str) or not isinstance(
            reverse_sequence, str
        ):
            raise SkillError(f"external-primer handoff pair {pair_id} lacks sequences")
        sources = pair.get("sources")
        if not isinstance(sources, list) or not sources:
            raise SkillError(f"external-primer handoff pair {pair_id} lacks source rows")
        if pair.get("vendor_claims_used_as_biological_evidence") is not False:
            raise SkillError(
                f"external-primer handoff pair {pair_id} did not preserve claim isolation"
            )
        specificity = pair.get("specificity")
        specificity_status = (
            str(specificity.get("status") or "missing")
            if isinstance(specificity, dict)
            else "missing"
        )
        specificity_status_counts[specificity_status] = (
            specificity_status_counts.get(specificity_status, 0) + 1
        )

        assay = pair.get("cdna_assay")
        assay_test_id = assay.get("assay_test_id") if isinstance(assay, dict) else None
        if not isinstance(assay_test_id, str) or not assay_test_id:
            raise SkillError(
                f"external-primer handoff pair {pair_id} lacks cDNA assay identity"
            )
        if assay.get("schema") != CDNA_ASSAY_TEST_REPORT_SCHEMA:
            raise SkillError(
                f"external-primer handoff pair {pair_id} has unexpected cDNA assay schema"
            )
        expected_assay_fields = {
            "assay_kind": "pcr",
            "source_seq_id": target["seq_id"],
            "source_feature_id": target["source_feature_id"],
            "requested_transcript_id": target["transcript_id"],
            "forward_primer": forward_sequence,
            "reverse_primer": reverse_sequence,
        }
        evaluation = context["request"]["evaluation"]
        optional_assay_fields = {
            "min_amplicon_bp": evaluation["min_amplicon_bp"],
            "max_amplicon_bp": evaluation["max_amplicon_bp"],
            "max_mismatches": evaluation["max_mismatches"],
            "require_3prime_exact_bases": evaluation[
                "require_3prime_exact_bases"
            ],
            "transcript_order": evaluation["transcript_order"],
            "transcript_map_coordinate_mode": evaluation[
                "map_coordinate_mode"
            ],
        }
        expected_assay_fields.update(
            {
                field: value
                for field, value in optional_assay_fields.items()
                if value is not None
            }
        )
        for field, expected_value in expected_assay_fields.items():
            if assay.get(field) != expected_value:
                raise SkillError(
                    f"external-primer handoff pair {pair_id} cDNA assay "
                    f"field mismatch: {field}"
                )

        specificity_report = (
            specificity.get("report") if isinstance(specificity, dict) else None
        )
        specificity_target = evaluation["specificity_target_genome_id"]
        if specificity_target is None:
            if specificity_status != "not_run" or specificity_report is not None:
                raise SkillError(
                    f"external-primer handoff pair {pair_id} returned unrequested "
                    "specificity evidence"
                )
        elif specificity_status == "not_run":
            raise SkillError(
                f"external-primer handoff pair {pair_id} omitted requested specificity"
            )
        elif specificity_report is not None:
            if not isinstance(specificity_report, dict):
                raise SkillError(
                    f"external-primer handoff pair {pair_id} has invalid specificity report"
                )
            expected_specificity_fields = {
                "schema": PRIMER_SPECIFICITY_REPORT_SCHEMA,
                "target_genome_id": specificity_target,
                "catalog_path": evaluation["specificity_catalog_path"],
                "cache_dir": evaluation["specificity_cache_dir"],
            }
            for field, expected_value in expected_specificity_fields.items():
                if specificity_report.get(field) != expected_value:
                    raise SkillError(
                        f"external-primer handoff pair {pair_id} specificity "
                        f"field mismatch: {field}"
                    )

        materialization = pair.get("product_materialization")
        if evaluation["materialize_products"]:
            if not isinstance(materialization, dict):
                raise SkillError(
                    f"external-primer handoff pair {pair_id} omitted requested "
                    "product materialization"
                )
        elif materialization is not None:
            raise SkillError(
                f"external-primer handoff pair {pair_id} returned unrequested "
                "product materialization"
            )
        for source in sources:
            if not isinstance(source, dict):
                raise SkillError(
                    f"external-primer handoff pair {pair_id} contains invalid source row"
                )
            annotations = source.get("annotations")
            if not isinstance(annotations, dict):
                raise SkillError(
                    f"external-primer handoff pair {pair_id} source lacks annotations"
                )
            record_id = annotations.get(EXTERNAL_PRIMER_HANDOFF_ANNOTATION_RECORD_ID)
            if not isinstance(record_id, str) or record_id not in expected_by_record_id:
                raise SkillError(
                    f"external-primer handoff pair {pair_id} returned unknown record_id"
                )
            if record_id in seen_record_ids:
                raise SkillError(
                    f"external-primer handoff record_id '{record_id}' was returned more than once"
                )
            expected = expected_by_record_id[record_id]
            expected_batch = _external_primer_batch_record(
                expected, context["request"]["collection_id"]
            )
            for field in (
                "source_kind",
                "provider",
                "catalogue_id",
                "source_url",
                "claimed_accession",
                "aliases",
                "claimed_target",
                "validation_claims",
                "annotations",
            ):
                if source.get(field) != expected_batch[field]:
                    raise SkillError(
                        f"external-primer handoff source field mismatch for "
                        f"record_id '{record_id}': {field}"
                    )
            if source.get("input_format") != "json":
                raise SkillError(
                    f"external-primer handoff record_id '{record_id}' lost input format"
                )
            if source.get("source_path") != context["batch_path"]:
                raise SkillError(
                    f"external-primer handoff record_id '{record_id}' lost source path"
                )
            if source.get("source_sha256") != context["batch_file_sha256"]:
                raise SkillError(
                    f"external-primer handoff record_id '{record_id}' lost source digest"
                )
            if (
                source.get("claim_evidence_status")
                != "provenance_only_not_used_for_coverage_or_specificity"
            ):
                raise SkillError(
                    f"external-primer handoff record_id '{record_id}' lost claim isolation"
                )
            if forward_sequence != expected["forward_sequence_5_to_3"] or (
                reverse_sequence != expected["reverse_sequence_5_to_3"]
            ):
                raise SkillError(
                    f"external-primer handoff sequence mismatch for record_id '{record_id}'"
                )
            source_record_id = source.get("source_record_id")
            if not isinstance(source_record_id, str) or not source_record_id:
                raise SkillError(
                    f"external-primer handoff record_id '{record_id}' lacks source_record_id"
                )
            if source_record_id in seen_source_record_ids:
                raise SkillError(
                    f"external-primer handoff source_record_id '{source_record_id}' is duplicated"
                )
            seen_record_ids.add(record_id)
            seen_source_record_ids.add(source_record_id)
            digest_rows.append([source_record_id, forward_sequence, reverse_sequence])
            joins.append(
                {
                    "record_id": record_id,
                    "assay_purpose": expected["assay_purpose"],
                    "source_record_id": source_record_id,
                    "pair_id": pair_id,
                    "forward_sequence_5_to_3": forward_sequence,
                    "reverse_sequence_5_to_3": reverse_sequence,
                    "assay_test_id": assay_test_id,
                    "report_id": report_id,
                    "status": "computed",
                }
            )

        artifacts = pair.get("artifacts")
        if not isinstance(artifacts, dict):
            raise SkillError(f"external-primer handoff pair {pair_id} lacks artifact paths")
        artifact_specs = (
            ("cdna_report_json_path", "cdna_assay_report", "application/json"),
            ("transcript_map_svg_path", "transcript_map", "image/svg+xml"),
            ("product_gel_svg_path", "product_gel", "image/svg+xml"),
        )
        for field, artifact_kind, media_type in artifact_specs:
            raw_path = artifacts.get(field)
            if raw_path is None:
                product_gel_required = (
                    field == "product_gel_svg_path"
                    and context["request"]["evaluation"]["materialize_products"]
                )
                if field != "product_gel_svg_path" or product_gel_required:
                    raise SkillError(
                        f"external-primer handoff pair {pair_id} lacks {field}"
                    )
                continue
            if not isinstance(raw_path, str) or not raw_path:
                raise SkillError(
                    f"external-primer handoff pair {pair_id} has invalid {field}"
                )
            artifact_path = Path(raw_path).resolve()
            if not artifact_path.is_relative_to(artifact_root):
                raise SkillError(
                    f"external-primer handoff pair {pair_id} artifact escaped output root"
                )
            if not artifact_path.is_file():
                raise SkillError(
                    f"external-primer handoff pair {pair_id} artifact is missing: {artifact_path}"
                )
            if field == "cdna_report_json_path":
                try:
                    artifact_payload = json.loads(
                        artifact_path.read_text(encoding="utf-8")
                    )
                except json.JSONDecodeError as e:
                    raise SkillError(
                        f"external-primer handoff pair {pair_id} cDNA artifact is invalid JSON"
                    ) from e
                if artifact_payload != assay:
                    raise SkillError(
                        f"external-primer handoff pair {pair_id} cDNA artifact/report mismatch"
                    )
            scientific_artifacts.append(
                _external_primer_handoff_artifact_record(
                    path=artifact_path,
                    artifact_kind=artifact_kind,
                    pair_id=pair_id,
                    media_type=media_type,
                )
            )

    if seen_record_ids != set(expected_by_record_id):
        missing = sorted(set(expected_by_record_id) - seen_record_ids)
        raise SkillError(
            "external-primer handoff report omitted submitted record_id value(s): "
            + ", ".join(missing)
        )
    digest_rows.sort()
    expected_normalized_batch_sha256 = _sha256_prefixed_bytes(
        json.dumps(digest_rows, separators=(",", ":"), ensure_ascii=True).encode(
            "utf-8"
        )
    )
    if report.get("normalized_batch_sha256") != expected_normalized_batch_sha256:
        raise SkillError("external-primer handoff normalized batch SHA-256 mismatch")

    operation_result = stdout_json.get("result")
    if not isinstance(operation_result, dict):
        raise SkillError("external-primer handoff command returned no operation result")
    op_id = operation_result.get("op_id")
    if not isinstance(op_id, str) or not op_id:
        raise SkillError("external-primer handoff operation result lacks op_id")
    if operation_result.get("external_primer_pair_import_report") != report:
        raise SkillError(
            "external-primer handoff command report and operation report differ"
        )

    report_path = Path(context["report_path"])
    if not report_path.is_file():
        raise SkillError("external-primer handoff report artifact is missing")
    try:
        report_file_payload = json.loads(report_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as e:
        raise SkillError("external-primer handoff report artifact is invalid JSON") from e
    if report_file_payload != report:
        raise SkillError("external-primer handoff report artifact does not match native result")
    scientific_artifacts.insert(
        0,
        _external_primer_handoff_artifact_record(
            path=report_path,
            artifact_kind="external_primer_import_report",
            media_type="application/json",
        ),
    )

    state_path = Path(context["state_path"])
    state_sha256_after = _sha256_file_prefixed(state_path)
    report_payload_sha256 = _sha256_prefixed_json(report)
    execution_binding_sha256 = _sha256_prefixed_json(
        {
            "canonical_request_sha256": context["canonical_request_sha256"],
            "state_sha256_before": context["state_sha256_before"],
            "submitted_batch_file_sha256": context["batch_file_sha256"],
            "normalized_batch_sha256": report["normalized_batch_sha256"],
            "report_id": report_id,
            "report_payload_sha256": report_payload_sha256,
            "op_id": op_id,
            "gentle_version": gentle_version,
            "state_sha256_after": state_sha256_after,
        }
    )
    incomplete_specificity = {
        status
        for status in specificity_status_counts
        if status in {"error", "incomplete", "missing"}
    }
    result = _external_primer_handoff_base_result(
        context,
        status="incomplete" if incomplete_specificity else "complete",
        gentle_version=gentle_version,
        state_sha256_after=state_sha256_after,
    )
    result.update(
        {
            "analysis_completeness": (
                "incomplete_specificity"
                if incomplete_specificity
                else "complete_for_requested_scope"
            ),
            "execution_binding_sha256": execution_binding_sha256,
            "submitted_batch_path": context["batch_path"],
            "submitted_batch_file_sha256": context["batch_file_sha256"],
            "submitted_record_joins": sorted(
                joins, key=lambda row: row["record_id"]
            ),
            "report_binding": {
                "report_id": report_id,
                "report_schema": report["schema"],
                "normalized_batch_sha256": report["normalized_batch_sha256"],
                "report_payload_sha256": report_payload_sha256,
                "report_file_path": str(report_path),
                "report_file_sha256": _sha256_file_prefixed(report_path),
                "op_id": op_id,
            },
            "specificity_requested": bool(
                context["request"]["evaluation"]["specificity_target_genome_id"]
            ),
            "specificity_status_counts": dict(
                sorted(specificity_status_counts.items())
            ),
            "scientific_artifacts": sorted(
                scientific_artifacts,
                key=lambda artifact: (
                    artifact.get("pair_id", ""),
                    artifact["artifact_kind"],
                    artifact["path"],
                ),
            ),
        }
    )
    result["warnings"] = list(report.get("warnings") or [])
    if context["not_submitted_records"]:
        result["warnings"].append(
            "Cloning/sequencing records were retained as not_submitted and were "
            "not analyzed by the PCR importer."
        )
    return result


def _external_primer_handoff_chat_summary_lines(
    result: dict[str, Any]
) -> list[str]:
    lines = [
        "External-primer handoff status: "
        f"{result['status']} ({result.get('analysis_completeness', 'unknown')}).",
        "Submitted PCR-compatible records: "
        f"{result['submitted_record_count']}; not submitted: "
        f"{result['not_submitted_record_count']}.",
    ]
    report_binding = result.get("report_binding")
    if isinstance(report_binding, dict):
        lines.append(
            "Bound GENtle report: "
            f"{report_binding.get('report_id')} at "
            f"{report_binding.get('normalized_batch_sha256')}."
        )
    else:
        lines.append("No GENtle biological evaluation was produced.")
    return lines


def _summary_lines_from_sequence_context_view(candidate: Any) -> list[str] | None:
    if not isinstance(candidate, dict):
        return None
    if candidate.get("schema") != "gentle.sequence_context_view.v1":
        return None
    raw_lines = candidate.get("summary_lines")
    if not isinstance(raw_lines, list):
        return None
    lines = [line.strip() for line in raw_lines if isinstance(line, str) and line.strip()]
    return lines or None


def _summary_lines_from_protein_residue_genomic_coordinates(
    candidate: Any,
) -> list[str] | None:
    if not isinstance(candidate, dict):
        return None
    report = candidate
    if candidate.get("schema") != "gentle.protein_residue_genomic_coordinates.v1":
        nested = candidate.get("protein_residue_genomic_coordinates")
        if not isinstance(nested, dict):
            return None
        if nested.get("schema") != "gentle.protein_residue_genomic_coordinates.v1":
            return None
        report = nested
    seq_id = str(report.get("seq_id") or "").strip()
    match_count = int(report.get("match_count") or 0)
    residue_start = int(report.get("residue_start_1based") or 0)
    residue_end = int(report.get("residue_end_1based") or residue_start)
    transcript_ids = {
        str(match.get("transcript_id") or "").strip()
        for match in report.get("matches", [])
        if isinstance(match, dict) and str(match.get("transcript_id") or "").strip()
    }
    lines = [
        f"Mapped residue(s) {residue_start}..{residue_end} on '{seq_id}' across {match_count} transcript match(es)."
    ]
    if transcript_ids:
        lines.append(
            f"Transcript matches: {len(transcript_ids)} ({', '.join(sorted(transcript_ids)[:3])})."
        )
    matches = report.get("matches")
    if isinstance(matches, list) and matches:
        first = matches[0]
        if isinstance(first, dict):
            transcript_id = str(first.get("transcript_id") or "").strip() or "transcript"
            residue_index = int(first.get("residue_index_1based") or 0)
            amino_acid = str(first.get("amino_acid") or "X").strip() or "X"
            codon = str(first.get("codon") or "").strip()
            genomic_bases = [
                str(base.get("genomic_pos_1based"))
                for base in first.get("genomic_bases", [])
                if isinstance(base, dict) and base.get("genomic_pos_1based") is not None
            ]
            detail = (
                f"First match: {transcript_id} residue {residue_index} ({amino_acid}; codon {codon})"
            )
            if genomic_bases:
                detail += f" maps to genomic bases {', '.join(genomic_bases)}"
            detail += "."
            lines.append(detail)
            if bool(first.get("spans_exon_junction")):
                lines.append("The codon spans an exon junction.")
    return lines


def _summary_lines_from_primer_qpcr_payload(candidate: Any) -> list[str] | None:
    if not isinstance(candidate, dict):
        return None
    schema = str(candidate.get("schema") or "").strip()
    if schema == "gentle.primer_seed_request.v1":
        template = str(candidate.get("template") or "template").strip()
        roi_start = candidate.get("roi_start_0based")
        roi_end = candidate.get("roi_end_0based_exclusive")
        return [
            f"Prepared a PCR/qPCR primer seed for '{template}' over ROI {roi_start}..{roi_end}.",
            "The payload includes ready-to-run DesignPrimerPairs and DesignQpcrAssays operations.",
        ]
    if schema == "gentle.qpcr_seed_request.v1":
        template = str(candidate.get("template") or "template").strip()
        roi_start = candidate.get("roi_start_0based")
        roi_end = candidate.get("roi_end_0based_exclusive")
        rationale = candidate.get("rationale")
        summary = ""
        if isinstance(rationale, dict):
            summary = str(rationale.get("summary") or "").strip()
        lines = [
            f"Prepared a probe-based qPCR/TaqMan seed for '{template}' over ROI {roi_start}..{roi_end}.",
        ]
        if summary:
            lines.append(summary)
        lines.append("The payload carries one ready-to-run DesignQpcrAssays operation and the pcr.assay.qpcr cartoon id.")
        return lines
    if schema == "gentle.primer_design_report.v1":
        report_id = str(candidate.get("report_id") or "primer report").strip()
        pair_count = len(candidate.get("pairs", [])) if isinstance(candidate.get("pairs"), list) else 0
        return [f"Primer-design report '{report_id}' contains {pair_count} ranked PCR primer pair(s)."]
    if schema == "gentle.qpcr_design_report.v1":
        report_id = str(candidate.get("report_id") or "qPCR report").strip()
        assay_count = len(candidate.get("assays", [])) if isinstance(candidate.get("assays"), list) else 0
        lines = [
            f"qPCR/TaqMan design report '{report_id}' contains {assay_count} ranked probe-bearing assay(s)."
        ]
        best_summary = str(candidate.get("best_assay_summary") or "").strip()
        if best_summary:
            lines.append(best_summary)
        assays = candidate.get("assays")
        first_assay = assays[0] if isinstance(assays, list) and assays else None
        if isinstance(first_assay, dict):
            context = first_assay.get("transcript_context")
            if isinstance(context, dict):
                risk = str(context.get("genomic_carryover_risk") or "").strip()
                rationale = str(context.get("genomic_carryover_rationale") or "").strip()
                if risk:
                    line = f"Best assay genomic-DNA carryover risk: {risk}."
                    if rationale:
                        line += f" {rationale}"
                    lines.append(line)
        return lines
    if schema == "gentle.cdna_assay_test_command.v1":
        report = candidate.get("cdna_assay_test_report") or candidate.get("report")
        lines = _summary_lines_from_primer_qpcr_payload(report) or []
        materialization = candidate.get("materialization")
        if isinstance(materialization, dict):
            product_count = materialization.get("product_count")
            seq_ids = materialization.get("product_seq_ids")
            seq_count = len(seq_ids) if isinstance(seq_ids, list) else 0
            created_seq_ids = materialization.get("created_product_seq_ids")
            reused_seq_ids = materialization.get("reused_product_seq_ids")
            created_count = len(created_seq_ids) if isinstance(created_seq_ids, list) else 0
            reused_count = len(reused_seq_ids) if isinstance(reused_seq_ids, list) else 0
            container_id = str(materialization.get("container_id") or "").strip()
            gel_path = str(materialization.get("product_gel_svg_path") or "").strip()
            if seq_count:
                if created_count and reused_count:
                    line = (
                        f"Materialized {created_count} new and reused {reused_count} existing "
                        "cDNA assay product sequence(s)"
                    )
                elif reused_count and not created_count:
                    line = f"Reused {reused_count} existing cDNA assay product sequence(s)"
                else:
                    line = f"Materialized {seq_count} cDNA assay product sequence(s)"
                if container_id:
                    line += f" into product container '{container_id}'"
                line += "."
                lines.append(line)
            elif product_count == 0:
                lines.append("No product vial was created because the assay detected 0 products.")
            gel_summary_lines = materialization.get("gel_summary_lines")
            if isinstance(gel_summary_lines, list):
                for gel_line in gel_summary_lines[:6]:
                    text = str(gel_line or "").strip()
                    if text:
                        lines.append(text)
                if len(gel_summary_lines) > 6:
                    lines.append(
                        f"Additional gel rows omitted from chat summary: {len(gel_summary_lines) - 6}."
                    )
            if gel_path:
                lines.append(f"Product gel SVG: {gel_path}")
        return lines or None
    if schema == "gentle.cdna_assay_test_report.v1":
        assay_kind = str(candidate.get("assay_kind") or "cDNA assay").strip()
        status = str(candidate.get("overall_status") or candidate.get("status") or "unknown").strip()
        transcript_count = candidate.get("transcript_count")
        product_count = candidate.get("product_count")
        line = f"{assay_kind} test finished with status '{status}'"
        if transcript_count is not None:
            line += f" across {transcript_count} transcript template(s)"
        if product_count is not None:
            line += f" and {product_count} detected product(s)"
        line += "."
        lines = [line]
        risk = candidate.get("genomic_carryover_risk")
        if isinstance(risk, dict):
            risk_level = str(risk.get("risk_level") or "").strip()
            risk_summary = str(risk.get("summary") or "").strip()
            if risk_level:
                detail = f"Genomic-DNA carryover risk: {risk_level}."
                if risk_summary:
                    detail += f" {risk_summary}"
                lines.append(detail)
        return lines
    if schema == "gentle.transcript_qpcr_panel.v1":
        transcript_count = candidate.get("transcript_count")
        group_label = str(candidate.get("group_label") or "splicing group").strip()
        return [
            f"Built a transcript qPCR panel for {group_label} across {transcript_count} transcript row(s).",
            "Rows reuse shared qPCR components when possible and report characteristic forward-primer evidence explicitly.",
        ]
    return None


def _extract_sequence_context_bundle(stdout_json: Any) -> dict[str, Any] | None:
    if not isinstance(stdout_json, dict):
        return None
    if stdout_json.get("schema") == "gentle.sequence_context_bundle.v1":
        return stdout_json
    candidate = stdout_json.get("sequence_context_bundle")
    if isinstance(candidate, dict) and candidate.get("schema") == "gentle.sequence_context_bundle.v1":
        return candidate
    return None


def _extract_preferred_artifacts(stdout_json: Any) -> list[dict[str, Any]] | None:
    if isinstance(stdout_json, dict) and isinstance(
        stdout_json.get("preferred_artifacts"), list
    ):
        artifacts = [
            artifact
            for artifact in stdout_json["preferred_artifacts"]
            if isinstance(artifact, dict) and isinstance(artifact.get("path"), str)
        ]
        if artifacts:
            artifacts.sort(
                key=lambda artifact: (
                    int(artifact.get("presentation_rank", 10**9))
                    if isinstance(artifact.get("presentation_rank"), int)
                    else 10**9,
                    str(artifact.get("artifact_id", "")),
                )
            )
            return artifacts

    bundle = _extract_sequence_context_bundle(stdout_json)
    if not isinstance(bundle, dict):
        return None
    raw_artifacts = bundle.get("artifacts")
    if not isinstance(raw_artifacts, list):
        return None
    artifacts = [
        artifact
        for artifact in raw_artifacts
        if isinstance(artifact, dict) and isinstance(artifact.get("path"), str)
    ]
    if not artifacts:
        return None
    artifacts.sort(
        key=lambda artifact: (
            int(artifact.get("presentation_rank", 10**9))
            if isinstance(artifact.get("presentation_rank"), int)
            else 10**9,
            str(artifact.get("artifact_id", "")),
        )
    )
    return artifacts


def _summary_projection(
    lines: list[str],
    *,
    source_pointer: str,
    source_payload: Any,
    projection_id: str,
    line_pointers: list[str] | None = None,
) -> dict[str, Any]:
    return {
        "lines": lines,
        "source_document": "result.json",
        "source_pointer": source_pointer,
        "source_payload_sha256": _sha256_prefixed_json(source_payload),
        "projection_id": projection_id,
        "line_pointers": line_pointers,
    }


def _extract_chat_summary_projection(stdout_json: Any) -> dict[str, Any] | None:
    if not isinstance(stdout_json, dict):
        return None
    bundle = _extract_sequence_context_bundle(stdout_json)
    bundle_pointer = (
        "/stdout_json"
        if bundle is stdout_json
        else "/stdout_json/sequence_context_bundle"
    )
    candidates = [
        (stdout_json, "/stdout_json"),
        (stdout_json.get("sequence_context_view"), "/stdout_json/sequence_context_view"),
        (bundle, bundle_pointer),
        (
            bundle.get("sequence_context_view") if isinstance(bundle, dict) else None,
            bundle_pointer + "/sequence_context_view",
        ),
    ]
    for candidate, pointer in candidates:
        lines = _summary_lines_from_sequence_context_view(candidate)
        if lines:
            raw_lines = candidate.get("summary_lines") if isinstance(candidate, dict) else []
            line_pointers = [
                f"{pointer}/summary_lines/{index}"
                for index, line in enumerate(raw_lines)
                if isinstance(line, str) and line.strip()
            ]
            return _summary_projection(
                lines,
                source_pointer=pointer,
                source_payload=candidate,
                projection_id="gentle.clawbio.summary_lines.identity_trim.v1",
                line_pointers=line_pointers,
            )
    protein_pointer = "/stdout_json"
    protein_payload = stdout_json
    if stdout_json.get("schema") != "gentle.protein_residue_genomic_coordinates.v1":
        nested_protein = stdout_json.get("protein_residue_genomic_coordinates")
        if isinstance(nested_protein, dict):
            protein_pointer = "/stdout_json/protein_residue_genomic_coordinates"
            protein_payload = nested_protein
    protein_residue_lines = _summary_lines_from_protein_residue_genomic_coordinates(
        stdout_json
    )
    if protein_residue_lines:
        return _summary_projection(
            protein_residue_lines,
            source_pointer=protein_pointer,
            source_payload=protein_payload,
            projection_id="gentle.clawbio.protein_residue_summary.v1",
        )
    primer_lines = _summary_lines_from_primer_qpcr_payload(stdout_json)
    if primer_lines:
        return _summary_projection(
            primer_lines,
            source_pointer="/stdout_json",
            source_payload=stdout_json,
            projection_id="gentle.clawbio.primer_qpcr_summary.v1",
        )
    if isinstance(stdout_json.get("output"), dict):
        primer_lines = _summary_lines_from_primer_qpcr_payload(stdout_json["output"])
        if primer_lines:
            return _summary_projection(
                primer_lines,
                source_pointer="/stdout_json/output",
                source_payload=stdout_json["output"],
                projection_id="gentle.clawbio.primer_qpcr_summary.v1",
            )
    for key in (
        "primer_design_report",
        "qpcr_design_report",
        "cdna_assay_test_report",
        "report",
        "transcript_qpcr_panel",
    ):
        nested = stdout_json.get(key)
        primer_lines = _summary_lines_from_primer_qpcr_payload(nested)
        if primer_lines:
            return _summary_projection(
                primer_lines,
                source_pointer=f"/stdout_json/{_json_pointer_token(key)}",
                source_payload=nested,
                projection_id="gentle.clawbio.primer_qpcr_summary.v1",
            )
    raw_lines = stdout_json.get("summary_lines")
    if isinstance(raw_lines, list):
        lines = [
            line.strip()
            for line in raw_lines
            if isinstance(line, str) and line.strip()
        ]
        if lines:
            line_pointers = [
                f"/stdout_json/summary_lines/{index}"
                for index, line in enumerate(raw_lines)
                if isinstance(line, str) and line.strip()
            ]
            return _summary_projection(
                lines,
                source_pointer="/stdout_json/summary_lines",
                source_payload=raw_lines,
                projection_id="gentle.clawbio.summary_lines.identity_trim.v1",
                line_pointers=line_pointers,
            )
    return None


def _extract_chat_summary_lines(stdout_json: Any) -> list[str] | None:
    projection = _extract_chat_summary_projection(stdout_json)
    return projection["lines"] if projection is not None else None


def _runtime_version_chat_summary_lines(
    request: Request | None,
    run_result: subprocess.CompletedProcess[str] | None,
) -> list[str] | None:
    if request is None or request.mode != "version" or run_result is None:
        return None
    version_text = run_result.stdout.strip() or run_result.stderr.strip()
    if not version_text:
        return ["GENtle runtime version command completed, but did not print a version."]
    return [
        f"Installed local GENtle rewrite runtime in this ClawBio environment: {version_text}",
        CLASSICAL_GENTLE_DISAMBIGUATION,
    ]


def _append_continue_artifact_notice(
    lines: list[str] | None,
    continue_actions: list[dict[str, Any]] | None,
) -> list[str] | None:
    if not continue_actions:
        return lines
    updated = list(lines or [])
    if CONTINUE_ARTIFACT_NOTICE not in updated:
        updated.append(CONTINUE_ARTIFACT_NOTICE)
    return updated


def _strict_claim_attribution(
    *,
    request: Request | None,
    lines: list[str] | None,
    default_source_kind: str | None,
    summary_projection: dict[str, Any] | None,
    stdout_json: Any,
    command: list[str] | None,
    collected_artifacts: list[dict[str, Any]] | None = None,
    warning_lines: list[str] | None = None,
) -> tuple[list[str] | None, dict[str, Any] | None]:
    if (
        request is None
        or request.claim_attribution_mode != STRICT_CLAIM_ATTRIBUTION_MODE
    ):
        return lines, None

    prefixes = CLAIM_SOURCE_PREFIXES
    attributed_lines: list[str] = []
    claims: list[dict[str, Any]] = []
    for ordinal, line in enumerate(lines or [], start=1):
        source_kind = default_source_kind or "clawbio_presentation"
        if line == CONTINUE_ARTIFACT_NOTICE:
            source_kind = "clawbio_presentation"
        if source_kind == "gentle_executable" and summary_projection is None:
            source_kind = "unverified"
        prefix = prefixes[source_kind]
        attributed = _prefix_claim_line(line, source_kind)
        attributed_lines.append(attributed)
        evidence: dict[str, Any] | None = None
        if source_kind == "gentle_executable" and summary_projection is not None:
            line_pointers = summary_projection.get("line_pointers")
            pointer = summary_projection["source_pointer"]
            if isinstance(line_pointers, list) and ordinal <= len(line_pointers):
                pointer = line_pointers[ordinal - 1]
            pointed_value = _result_pointer_value(stdout_json, pointer)
            evidence = {
                "source_document": summary_projection["source_document"],
                "json_pointer": pointer,
                "source_value_sha256": _sha256_prefixed_json(pointed_value),
                "source_scope_pointer": summary_projection["source_pointer"],
                "source_scope_sha256": summary_projection[
                    "source_payload_sha256"
                ],
                "projection_id": summary_projection["projection_id"],
            }
        elif summary_projection is not None and line != CONTINUE_ARTIFACT_NOTICE:
            evidence = {
                "source_document": summary_projection["source_document"],
                "json_pointer": summary_projection["source_pointer"],
                "source_value_sha256": summary_projection[
                    "source_payload_sha256"
                ],
                "projection_id": summary_projection["projection_id"],
            }
        elif line == CONTINUE_ARTIFACT_NOTICE:
            evidence = {
                "source_document": "result.json",
                "json_pointer": "/suggested_actions",
                "source_value_sha256": None,
                "projection_id": "gentle.clawbio.continuation_notice.v1",
            }
        claim_identity = _canonical_json_bytes(
            {
                "ordinal": ordinal,
                "source_kind": source_kind,
                "text": line,
                "evidence": evidence,
            }
        )
        claims.append(
            {
                "claim_id": "claim_sha256_"
                + hashlib.sha256(claim_identity).hexdigest(),
                "ordinal": ordinal,
                "prefix": prefix,
                "source_kind": source_kind,
                "source_tool": (
                    "gentle_cli"
                    if source_kind == "gentle_executable"
                    else "gentle-cloning"
                ),
                "text": line,
                "display_text": attributed,
                "evidence": evidence,
            }
        )

    input_claims = []
    for index, text in enumerate(request.input_claims or []):
        evidence = {
            "source_document": "result.json",
            "json_pointer": f"/request/input_claims/{index}",
            "source_value_sha256": _sha256_prefixed_json(text),
            "projection_id": "gentle.clawbio.input_claim.identity.v1",
        }
        claim_identity = _canonical_json_bytes(
            {
                "index": index,
                "source_kind": "caller_input",
                "text": text,
                "evidence": evidence,
            },
        )
        input_claims.append(
            {
                "claim_id": "claim_sha256_"
                + hashlib.sha256(claim_identity).hexdigest(),
                "ordinal": index + 1,
                "prefix": prefixes["caller_input"],
                "source_kind": "caller_input",
                "source_tool": "clawbio_request",
                "text": text,
                "display_text": _prefix_claim_line(text, "caller_input"),
                "evidence": evidence,
            }
        )

    request_payload = dataclasses.asdict(request)
    request_sha256 = _sha256_prefixed_json(request_payload)
    processing_steps = [
        {
            "prefix": prefixes["caller_input"],
            "source_kind": "caller_input",
            "tool": "clawbio_request",
            "role": "Supplied assumptions and paths; not a GENtle finding.",
            "payload_pointer": "/request",
            "payload_sha256": request_sha256,
        },
        {
            "prefix": prefixes["clawbio_presentation"],
            "source_kind": "clawbio_presentation",
            "tool": "gentle-cloning",
            "role": "Validated delegation, artifact collection, source labeling, and graphical assembly only.",
        },
    ]
    if command:
        processing_steps.append(
            {
                "prefix": prefixes["gentle_executable"],
                "source_kind": "gentle_executable",
                "tool": "gentle_cli",
                "role": "Authoritative scientific calculation and report generation.",
                "command": list(command),
                "output_pointer": "/stdout_json",
                "output_sha256": (
                    _sha256_prefixed_json(stdout_json)
                    if stdout_json is not None
                    else None
                ),
            }
        )
    artifact_attribution = []
    for artifact in collected_artifacts or []:
        declared_path = str(artifact.get("declared_path") or "").strip()
        bundle_path = str(artifact.get("bundle_path") or "").strip()
        copied_path = str(artifact.get("copied_path") or "").strip()
        is_clawbio_storyboard = declared_path == "generated/clawbio_storyboard.svg"
        source_kind = (
            "clawbio_presentation"
            if is_clawbio_storyboard
            else "gentle_executable"
        )
        identity = json.dumps(
            {
                "declared_path": declared_path,
                "bundle_path": bundle_path,
                "source_kind": source_kind,
            },
            sort_keys=True,
            separators=(",", ":"),
        )
        artifact_attribution.append(
            {
                "artifact_id": "artifact_sha256_"
                + hashlib.sha256(identity.encode("utf-8")).hexdigest(),
                "prefix": prefixes[source_kind],
                "source_kind": source_kind,
                "source_tool": (
                    "gentle-cloning"
                    if is_clawbio_storyboard
                    else "gentle_cli"
                ),
                "declared_path": declared_path,
                "bundle_path": bundle_path,
                "derived_from": artifact.get("derived_from"),
                "content_sha256": (
                    "sha256:" + _sha256_file(Path(copied_path))
                    if copied_path and Path(copied_path).is_file()
                    else None
                ),
                "scientific_content_authority": (
                    "none_presentation_only"
                    if is_clawbio_storyboard
                    else "gentle_executable"
                ),
            }
        )
    warning_claims = []
    for ordinal, warning in enumerate(warning_lines or [], start=1):
        source_kind = "clawbio_presentation"
        warning_claims.append(
            {
                "ordinal": ordinal,
                "prefix": prefixes[source_kind],
                "source_kind": source_kind,
                "source_tool": "gentle-cloning",
                "text": warning,
                "display_text": _prefix_claim_line(warning, source_kind),
                "evidence": {
                    "source_document": "result.json",
                    "json_pointer": f"/warnings/{ordinal - 1}",
                    "source_value_sha256": _sha256_prefixed_json(warning),
                    "projection_id": "gentle.clawbio.warning.identity.v1",
                },
            }
        )
    return (
        attributed_lines or None,
        {
            "schema": CLAIM_LEDGER_SCHEMA,
            "mode": STRICT_CLAIM_ATTRIBUTION_MODE,
            "presentation_profile": request.presentation_profile,
            "normalized_request_sha256": request_sha256,
            "native_result_sha256": (
                _sha256_prefixed_json(stdout_json) if stdout_json is not None else None
            ),
            "prefixes": prefixes,
            "processing_steps": processing_steps,
            "claims": claims,
            "input_claims": input_claims,
            "warning_claims": warning_claims,
            "artifacts": artifact_attribution,
            "rules": [
                "Only [gentle] statements are copied from the GENtle executable result.",
                "Every [gentle] statement requires a JSON pointer, source-value hash, source-scope hash, and named deterministic projection.",
                "[clawbio] statements describe orchestration or presentation and must not introduce biological findings.",
                "[input] values are caller-supplied assumptions, not validated findings.",
                "Caller text that begins with a reserved source prefix is escaped before display and cannot self-assign authority.",
                "Direct external-tool claims require [external:<tool>] plus tool/version/input/output provenance; this wrapper invokes GENtle rather than external scientific tools directly.",
                "Unverified prose must use [unverified] and cannot satisfy readiness or specificity gates.",
            ],
        },
    )


def _result_pointer_value(stdout_json: Any, pointer: str) -> Any:
    prefix = "/stdout_json"
    if pointer == prefix:
        return stdout_json
    if not pointer.startswith(prefix + "/"):
        raise SkillError(f"claim evidence pointer is outside stdout_json: {pointer}")
    current = stdout_json
    for raw_token in pointer[len(prefix) + 1 :].split("/"):
        token = raw_token.replace("~1", "/").replace("~0", "~")
        if isinstance(current, dict) and token in current:
            current = current[token]
        elif isinstance(current, list) and token.isdigit() and int(token) < len(current):
            current = current[int(token)]
        else:
            raise SkillError(f"claim evidence pointer does not resolve: {pointer}")
    return current


def _claim_source_kind_from_prefix(line: str) -> str:
    if line.startswith("[gentle]"):
        return "gentle_executable"
    if line.startswith("[input]"):
        return "caller_input"
    if line.startswith("[external:"):
        return "external_tool"
    if line.startswith("[unverified]"):
        return "unverified"
    return "clawbio_presentation"


def _claim_prefix_from_line(line: str, source_kind: str) -> str:
    if source_kind == "external_tool" and line.startswith("[external:"):
        end = line.find("]")
        if end > 0:
            return line[: end + 1]
    return CLAIM_SOURCE_PREFIXES[source_kind]


def _claim_source_tool(line: str, source_kind: str) -> str:
    if source_kind == "gentle_executable":
        return "gentle_cli"
    if source_kind == "caller_input":
        return "clawbio_request"
    if source_kind == "external_tool":
        prefix = _claim_prefix_from_line(line, source_kind)
        return prefix.removeprefix("[external:").removesuffix("]")
    if source_kind == "unverified":
        return "unverified"
    return "gentle-cloning"


def _prefix_claim_line(line: str, source_kind: str) -> str:
    match = re.match(
        r"^\[(?:gentle|clawbio|input|unverified|external:[^\]]+)\]\s*",
        line,
    )
    safe_line = line
    if match is not None:
        literal_prefix = match.group(0).strip()
        safe_line = f"(literal source prefix {literal_prefix}) {line[match.end():]}"
    return f"{CLAIM_SOURCE_PREFIXES[source_kind]} {safe_line}"


def _fallback_chat_summary_lines(
    *,
    request: Request | None,
    command: list[str] | None,
    run_result: subprocess.CompletedProcess[str] | None,
    stdout_json: Any | None,
    status: str,
) -> list[str] | None:
    if run_result is None:
        return None
    stdout_preview = _one_line_preview(run_result.stdout)
    stderr_preview = _one_line_preview(run_result.stderr)
    if (
        _is_default_demo_request(request)
        and not isinstance(stdout_json, dict)
        and not stdout_preview
        and not stderr_preview
    ):
        return [
            "Generated a deterministic GENtle protocol cartoon for a two-fragment Gibson assembly.",
            (
                "The ClawBio demo now starts with a graphical export so the first reply can "
                "show an actual figure instead of only listing commands."
            ),
            "Best-first preview artifact: generated/artifacts/gibson.two_fragment.protocol.png",
        ]
    if not isinstance(stdout_json, dict) and not stdout_preview and not stderr_preview:
        return None
    command_text = _format_command_text(command)
    lines = [
        f"GENtle command completed with status `{status}` and exit code `{run_result.returncode}`.",
        f"Command: {command_text}",
    ]
    if isinstance(stdout_json, dict):
        schema = stdout_json.get("schema")
        if isinstance(schema, str) and schema.strip():
            lines.append(f"Parsed JSON output schema: {schema.strip()}")
        else:
            keys = sorted(str(key) for key in stdout_json.keys())[:8]
            if keys:
                lines.append("Parsed JSON output keys: " + ", ".join(keys))
        if (
            isinstance(request, Request)
            and request.mode == "shell"
            and (request.shell_line or "").strip() == "capabilities"
        ):
            capabilities = stdout_json.get("capabilities")
            if isinstance(capabilities, list):
                lines.append(f"Capability entries reported: {len(capabilities)}")
    else:
        if stdout_preview:
            lines.append(f"Stdout preview: {stdout_preview}")
    if stderr_preview:
        lines.append(f"Stderr preview: {stderr_preview}")
    return lines


def _action_id_from_label(label: str) -> str:
    token = re.sub(r"[^a-z0-9]+", "_", label.lower()).strip("_")
    return token or "action"


def _make_shell_request(shell_line: str, timeout_secs: int) -> dict[str, Any]:
    return {
        "schema": REQUEST_SCHEMA,
        "mode": "shell",
        "shell_line": shell_line,
        "timeout_secs": timeout_secs,
    }


def _request_rerun_shell_line(request: Request | None) -> str:
    if request is None:
        return "(unknown)"
    if request.mode in {"shell", "external-primer-handoff"} and request.shell_line:
        return request.shell_line.strip()
    if request.mode == "workflow":
        if request.workflow_path:
            return f"workflow @{request.workflow_path}"
        return "workflow <inline>"
    if request.mode == "op":
        return "op <inline>"
    if request.mode in {
        "construct-reasoning-list-inspections",
        "construct-reasoning-run-inspection",
    }:
        return _build_construct_reasoning_inspection_shell_line(request)
    if request.mode in PRIMER_SHELL_REQUEST_MODES:
        return _build_primer_mode_shell_line(request)
    if request.mode == "protein-residue-genomic-coordinates":
        seq_id = request.seq_id or "SEQ_ID"
        start = request.residue_start_1based or 1
        end = request.residue_end_1based or start
        tokens = ["transcripts", "residue-genomic-coordinates", seq_id, str(start)]
        if end != start:
            tokens.append(str(end))
        if request.transcript_id:
            tokens.extend(["--transcript", request.transcript_id])
        return shlex.join(tokens)
    if request.mode == "transcript-qpcr-panel":
        tokens = [
            "primers",
            "transcript-qpcr-panel",
            request.seq_id or "SEQ_ID",
            str(request.source_feature_id if request.source_feature_id is not None else 0),
            request.shared_qpcr_report_id or "QPCR_REPORT_ID",
        ]
        if request.output_path:
            tokens.extend(["--path", request.output_path])
        return shlex.join(tokens)
    if request.mode == "pcr-protocol-cartoon":
        return _build_pcr_protocol_cartoon_shell_line(request)
    if request.mode == "gene-protein-2d-gel":
        gene = request.gene_symbol or "GENE"
        species = request.species or "homo_sapiens"
        source = request.source or "ensembl"
        return f"gene-protein-2d-gel {gene} --species {species} --source {source}"
    if request.mode in {"exon-skip-plan", "exon-skip-materialize"}:
        return _build_exon_skip_shell_line(request)
    if request.mode == "raw" and request.raw_args:
        return shlex.join(request.raw_args)
    return request.mode


def _request_payload_for_artifact_continuation(
    request: Request | None,
    expected_artifacts: list[str],
) -> dict[str, Any]:
    timeout_secs = request.timeout_secs if request is not None else 180
    payload: dict[str, Any] = {
        "schema": REQUEST_SCHEMA,
        "mode": request.mode if request is not None else "shell",
        "timeout_secs": timeout_secs,
        "expected_artifacts": expected_artifacts,
    }
    if request is None:
        payload["shell_line"] = ""
        return payload

    if request.state_path:
        payload["state_path"] = request.state_path
    if request.mode == "shell":
        payload["shell_line"] = request.shell_line or ""
    elif request.mode == "workflow":
        if request.workflow_path:
            payload["workflow_path"] = request.workflow_path
        else:
            payload["workflow"] = request.workflow
    elif request.mode == "op":
        payload["operation"] = request.operation
    elif request.mode in {
        "construct-reasoning-list-inspections",
        "construct-reasoning-run-inspection",
    }:
        for key in (
            "graph_id",
            "fact_id",
            "annotation_id",
            "candidate_id",
            "evidence_id",
            "seq_id",
            "action_kind",
            "summary_id",
            "action_id",
            "word_size",
            "step_bp",
            "max_mismatches",
            "tile_bp",
            "dotplot_id",
            "render_svg_path",
        ):
            value = getattr(request, key)
            if value is not None:
                payload[key] = value
    elif request.mode in PRIMER_SHELL_REQUEST_MODES:
        for key in (
            "request_json",
            "operation",
            "seq_id",
            "source_feature_id",
            "transcript_id",
            "transcript_order",
            "map_coordinate_mode",
            "qpcr_mode",
            "specificity_evidence",
            "backend",
            "primer3_executable",
            "report_id",
            "forward_primer",
            "reverse_primer",
            "probe",
            "min_amplicon_bp",
            "max_amplicon_bp",
            "max_mismatches",
            "require_3prime_exact_bases",
            "vector_seq_id",
            "pair_rank",
            "handoff_mode",
            "forward_enzyme",
            "reverse_enzyme",
            "forward_leader_5prime",
            "reverse_leader_5prime",
            "output_path",
            "svg_path",
        ):
            value = getattr(request, key)
            if value is not None:
                payload[key] = value
    elif request.mode == "protein-residue-genomic-coordinates":
        payload["seq_id"] = request.seq_id
        payload["transcript_id"] = request.transcript_id
        payload["residue_start_1based"] = request.residue_start_1based
        payload["residue_end_1based"] = request.residue_end_1based
    elif request.mode == "pcr-protocol-cartoon":
        payload["protocol_id"] = request.protocol_id
        payload["output_path"] = request.output_path
    elif request.mode == "gene-protein-2d-gel":
        payload["gene_symbol"] = request.gene_symbol
        payload["species"] = request.species
        payload["source"] = request.source
        payload["ladders"] = request.ladders
    elif request.mode in {"exon-skip-plan", "exon-skip-materialize"}:
        for key in (
            "seq_id",
            "transcript_feature_id",
            "skip_candidate_ids",
            "skip_intervals_1based",
            "overlap_intervals_1based",
            "length_mod3_values",
            "coding_mod3_values",
            "coding_contexts",
            "cds_phase_entry_kinds",
            "feature_query",
            "plan_id",
            "candidate_ids",
            "output_prefix",
            "return_items",
            "confirm",
        ):
            value = getattr(request, key)
            if value is not None:
                payload[key] = value
    elif request.mode == "raw":
        payload["raw_args"] = request.raw_args
    return {key: value for key, value in payload.items() if value is not None}


def _string_list(value: Any) -> list[str]:
    if not isinstance(value, list):
        return []
    return [str(item).strip() for item in value if isinstance(item, str) and item.strip()]


def _normalize_ui_intent_metadata(value: Any) -> dict[str, Any] | None:
    if not isinstance(value, dict):
        return None
    action = str(value.get("action") or "").strip()
    target = str(value.get("target") or "").strip()
    if not action or not target:
        return None
    normalized: dict[str, Any] = {
        "action": action,
        "target": target,
    }
    for key in ("title", "detail", "keywords", "menu_path"):
        field = str(value.get(key) or "").strip()
        if field:
            normalized[key] = field
    optional_arguments = _string_list(value.get("optional_arguments"))
    if optional_arguments:
        normalized["optional_arguments"] = optional_arguments
    return normalized


def _suggested_action(
    *,
    label: str,
    kind: str,
    shell_line: str,
    timeout_secs: int,
    rationale: str,
    requires_confirmation: bool = True,
    expected_artifacts: list[str] | None = None,
    resource_key: str | None = None,
    lifecycle_status: str | None = None,
    ui_intent: dict[str, Any] | None = None,
) -> dict[str, Any]:
    action = {
        "action_id": _action_id_from_label(label),
        "label": label,
        "kind": kind,
        "shell_line": shell_line,
        "timeout_secs": timeout_secs,
        "request": _make_shell_request(shell_line, timeout_secs),
        "rationale": rationale,
        "requires_confirmation": requires_confirmation,
    }
    if expected_artifacts:
        action["expected_artifacts"] = expected_artifacts
        action["request"]["expected_artifacts"] = expected_artifacts
    if resource_key:
        action["resource_key"] = resource_key
    if lifecycle_status:
        action["lifecycle_status"] = lifecycle_status
    if ui_intent:
        action["ui_intent"] = ui_intent
    return action


def _continue_artifact_suggested_actions(
    request: Request | None,
    collected_artifacts: list[dict[str, Any]],
    preferred_artifacts: list[dict[str, Any]] | None,
) -> list[dict[str, Any]] | None:
    svg_artifacts = _collected_svg_artifacts(collected_artifacts)
    if len(svg_artifacts) < 2:
        return None

    best_references: set[str] = set()
    for artifact in preferred_artifacts or []:
        if not isinstance(artifact, dict):
            continue
        if artifact.get("is_best_first_artifact") is True:
            best_references.update(_artifact_reference_keys(artifact))
    if not best_references and preferred_artifacts:
        best_references.update(_artifact_reference_keys(preferred_artifacts[0]))

    remaining = [
        artifact
        for artifact in svg_artifacts
        if not (_artifact_reference_keys(artifact) & best_references)
    ]
    if not remaining:
        return None

    shell_line = _request_rerun_shell_line(request)
    timeout_secs = request.timeout_secs if request is not None else 180
    actions: list[dict[str, Any]] = []
    for idx, artifact in enumerate(remaining, start=1):
        source_svg = str(artifact.get("declared_path") or "").strip()
        if not source_svg:
            continue
        label = "Continue: show next figure" if idx == 1 else f"Continue: show figure {idx + 1}"
        action = {
            "action_id": f"continue_show_figure_{idx}",
            "label": label,
            "kind": "continue_artifact",
            "source_shell_line": shell_line,
            "timeout_secs": timeout_secs,
            "request": _request_payload_for_artifact_continuation(
                request,
                [source_svg],
            ),
            "rationale": (
                "This run generated more than one displayable figure. Re-run the "
                "nested request payload while collecting only this next figure."
            ),
            "requires_confirmation": False,
            "expected_artifacts": [source_svg],
            "artifact": {
                "declared_path": source_svg,
                "svg_bundle_path": artifact.get("bundle_path"),
            },
        }
        actions.append(action)
    return actions or None


def _normalize_stdout_suggested_action(action: Any) -> dict[str, Any] | None:
    if not isinstance(action, dict):
        return None
    label = str(action.get("label") or "").strip()
    shell_line = str(action.get("shell_line") or "").strip()
    kind = str(action.get("kind") or "suggested_action").strip() or "suggested_action"
    if not label or not shell_line:
        return None
    raw_timeout = action.get("timeout_secs")
    try:
        timeout_secs = int(raw_timeout)
    except (TypeError, ValueError):
        timeout_secs = 180
    expected_artifacts = [
        str(path)
        for path in action.get("expected_artifacts", [])
        if isinstance(path, str) and path.strip()
    ]
    normalized = _suggested_action(
        label=label,
        kind=kind,
        shell_line=shell_line,
        timeout_secs=timeout_secs,
        rationale=str(action.get("rationale") or "").strip(),
        requires_confirmation=bool(action.get("requires_confirmation", True)),
        expected_artifacts=expected_artifacts or None,
        resource_key=str(action.get("resource_key") or "").strip() or None,
        lifecycle_status=str(action.get("lifecycle_status") or "").strip() or None,
        ui_intent=_normalize_ui_intent_metadata(action.get("ui_intent")),
    )
    action_id = str(action.get("action_id") or "").strip()
    if action_id:
        normalized["action_id"] = action_id
    return normalized


def _extract_normalized_action_list(
    stdout_json: Any,
    field_name: str,
) -> list[dict[str, Any]] | None:
    if not isinstance(stdout_json, dict):
        return None
    raw_actions = stdout_json.get(field_name)
    if not isinstance(raw_actions, list):
        return None
    actions = [
        normalized
        for normalized in (
            _normalize_stdout_suggested_action(action) for action in raw_actions
        )
        if normalized is not None
    ]
    return actions or None


def _extract_preferred_demo_actions(stdout_json: Any) -> list[dict[str, Any]] | None:
    return _extract_normalized_action_list(stdout_json, "preferred_demo_actions")


def _normalize_stdout_blocked_action(blocked_action: Any) -> dict[str, Any] | None:
    if not isinstance(blocked_action, dict):
        return None
    action = _normalize_stdout_suggested_action(blocked_action.get("action"))
    if action is None:
        return None
    normalized: dict[str, Any] = {
        "blocked_reason": str(blocked_action.get("blocked_reason") or "").strip(),
        "unblock_hint": str(blocked_action.get("unblock_hint") or "").strip(),
        "action": action,
    }
    for key in ("download_url", "local_path_hint"):
        value = str(blocked_action.get(key) or "").strip()
        if value:
            normalized[key] = value
    return normalized


def _extract_blocked_actions(stdout_json: Any) -> list[dict[str, Any]] | None:
    if not isinstance(stdout_json, dict):
        return None
    raw_actions = stdout_json.get("blocked_actions")
    if not isinstance(raw_actions, list):
        return None
    actions = [
        normalized
        for normalized in (
            _normalize_stdout_blocked_action(action) for action in raw_actions
        )
        if normalized is not None
    ]
    return actions or None


def _extract_suggested_actions(
    stdout_json: Any,
    request: Request | None,
) -> list[dict[str, Any]] | None:
    if not isinstance(stdout_json, dict):
        return None

    actions: list[dict[str, Any]] = []
    actions.extend(_extract_normalized_action_list(stdout_json, "suggested_actions") or [])

    if stdout_json.get("schema") == "gentle.service_readiness.v1":
        has_running_prepare = False
        for reference in stdout_json.get("references", []):
            if not isinstance(reference, dict):
                continue
            genome_id = str(reference.get("genome_id", "")).strip()
            lifecycle_status = str(
                reference.get("lifecycle_status")
                or reference.get("availability_status")
                or ""
            ).strip()
            if not genome_id:
                continue
            if lifecycle_status == "running":
                has_running_prepare = True
            if lifecycle_status in {"missing", "not_prepared"}:
                actions.append(
                    _suggested_action(
                        label=f"Prepare {genome_id}",
                        kind="prepare_reference",
                        shell_line=f'genomes prepare "{genome_id}" --timeout-secs 7200',
                        timeout_secs=7500,
                        rationale=(
                            f"Reference '{genome_id}' is currently not prepared locally "
                            "and genome-backed analysis will need a prepared cache."
                        ),
                    )
                )
            elif lifecycle_status in {"failed", "cancelled", "stale"}:
                actions.append(
                    _suggested_action(
                        label=f"Retry prepare for {genome_id}",
                        kind="prepare_reference",
                        shell_line=f'genomes prepare "{genome_id}" --timeout-secs 7200',
                        timeout_secs=7500,
                        rationale=(
                            f"Reference '{genome_id}' last ended as {lifecycle_status} and is "
                            "safe to retry when you want to restore genome-backed analysis."
                        ),
                    )
                )
        for helper in stdout_json.get("helpers", []):
            if not isinstance(helper, dict):
                continue
            helper_id = str(helper.get("genome_id", "")).strip()
            lifecycle_status = str(
                helper.get("lifecycle_status")
                or helper.get("availability_status")
                or ""
            ).strip()
            if not helper_id:
                continue
            if lifecycle_status == "running":
                has_running_prepare = True
            if lifecycle_status in {"missing", "not_prepared"}:
                actions.append(
                    _suggested_action(
                        label=f"Prepare {helper_id}",
                        kind="prepare_helper",
                        shell_line=f'helpers prepare "{helper_id}" --timeout-secs 1800',
                        timeout_secs=2100,
                        rationale=(
                            f"Helper '{helper_id}' is currently not prepared locally "
                            "and helper-backed vector/plasmid workflows will need it prepared."
                        ),
                    )
                )
            elif lifecycle_status in {"failed", "cancelled", "stale"}:
                actions.append(
                    _suggested_action(
                        label=f"Retry prepare for {helper_id}",
                        kind="prepare_helper",
                        shell_line=f'helpers prepare "{helper_id}" --timeout-secs 1800',
                        timeout_secs=2100,
                        rationale=(
                            f"Helper '{helper_id}' last ended as {lifecycle_status} and is "
                            "safe to retry when helper-backed workflows need it again."
                        ),
                    )
                )

        resources = stdout_json.get("resources")
        if isinstance(resources, dict):
            attract = resources.get("attract")
            if isinstance(attract, dict):
                support_status = str(attract.get("support_status", "")).strip()
                if support_status in {"known_external_only", "runtime_invalid"}:
                    actions.append(
                        _suggested_action(
                            label="Sync ATtRACT runtime snapshot",
                            kind="sync_resource",
                            shell_line="resources sync-attract ATtRACT.zip",
                            timeout_secs=900,
                            rationale=(
                                "ATtRACT is known to GENtle, but no valid runtime snapshot "
                                "is active yet."
                            ),
                        )
                    )
        if has_running_prepare:
            actions.append(
                _suggested_action(
                    label="Re-check services status",
                    kind="refresh_status",
                    shell_line="services status",
                    timeout_secs=180,
                    rationale=(
                        "A shared prepare action is already running, so refresh the combined "
                        "readiness view instead of starting a duplicate long-running task."
                    ),
                    requires_confirmation=False,
                )
            )

    prepare_command = stdout_json.get("prepare_command")
    lifecycle_status = str(stdout_json.get("lifecycle_status") or "").strip()
    if (
        isinstance(prepare_command, str)
        and prepare_command.strip()
        and lifecycle_status not in {"running", "ready"}
    ):
        requested_key = str(
            stdout_json.get("requested_catalog_key")
            or stdout_json.get("genome_id")
            or "selected genome"
        ).strip()
        actions.append(
            _suggested_action(
                label=f"Prepare {requested_key}",
                kind="prepare_reference",
                shell_line=prepare_command.strip(),
                timeout_secs=7500,
                rationale=(
                    f"Status inspection for '{requested_key}' already provided a ready-to-run "
                    "prepare command."
                ),
            )
        )
    elif (
        lifecycle_status == "running"
        and isinstance(request, Request)
        and request.mode == "shell"
        and (request.shell_line or "").strip().startswith(("genomes status ", "helpers status "))
    ):
        requested_key = str(
            stdout_json.get("requested_catalog_key")
            or stdout_json.get("genome_id")
            or "selected genome"
        ).strip()
        refresh_shell = (request.shell_line or "").strip()
        if refresh_shell:
            actions.append(
                _suggested_action(
                    label=f"Re-check {requested_key} status",
                    kind="refresh_status",
                    shell_line=refresh_shell,
                    timeout_secs=180,
                    rationale=(
                        f"'{requested_key}' is already being prepared, so the next useful step "
                        "is to refresh its status rather than launch another prepare."
                    ),
                    requires_confirmation=False,
                )
            )

    if stdout_json.get("schema") == "gentle.cutrun_dataset_status.v1":
        dataset_id = str(stdout_json.get("dataset_id") or "").strip()
        lifecycle_status = str(stdout_json.get("lifecycle_status") or "").strip()
        requested_catalog_path = str(
            stdout_json.get("requested_catalog_path") or ""
        ).strip()
        effective_cache_dir = str(stdout_json.get("effective_cache_dir") or "").strip()
        if dataset_id:
            dataset_arg = shlex.quote(dataset_id)
            catalog_arg = (
                f" --catalog {shlex.quote(requested_catalog_path)}"
                if requested_catalog_path
                else ""
            )
            cache_arg = (
                f" --cache-dir {shlex.quote(effective_cache_dir)}"
                if effective_cache_dir
                else ""
            )
            prepare_shell = f"cutrun prepare {dataset_arg}{catalog_arg}{cache_arg}"
            refresh_shell = f"cutrun status {dataset_arg}{catalog_arg}{cache_arg}"
            if (
                isinstance(request, Request)
                and request.mode == "shell"
                and (request.shell_line or "").strip().startswith("cutrun status ")
            ):
                refresh_shell = (request.shell_line or "").strip() or refresh_shell
            if lifecycle_status in {"missing", "not_prepared"}:
                actions.append(
                    _suggested_action(
                        label=f"Prepare CUT&RUN dataset {dataset_id}",
                        kind="prepare_cutrun_dataset",
                        shell_line=prepare_shell,
                        timeout_secs=1800,
                        rationale=(
                            f"CUT&RUN dataset '{dataset_id}' is not prepared locally and "
                            "must be materialized before dataset-backed projection or read "
                            "interpretation can reuse it."
                        ),
                    )
                )
            elif lifecycle_status in {"failed", "cancelled", "stale"}:
                actions.append(
                    _suggested_action(
                        label=f"Retry prepare for CUT&RUN dataset {dataset_id}",
                        kind="prepare_cutrun_dataset",
                        shell_line=prepare_shell,
                        timeout_secs=1800,
                        rationale=(
                            f"CUT&RUN dataset '{dataset_id}' last ended as {lifecycle_status} "
                            "and is safe to retry when you want to restore dataset-backed "
                            "projection or read interpretation."
                        ),
                    )
                )
            elif lifecycle_status == "running":
                actions.append(
                    _suggested_action(
                        label=f"Re-check CUT&RUN dataset {dataset_id} status",
                        kind="refresh_status",
                        shell_line=refresh_shell,
                        timeout_secs=180,
                        rationale=(
                            f"CUT&RUN dataset '{dataset_id}' is already being prepared, so "
                            "refresh its status instead of starting another parallel prepare."
                        ),
                        requires_confirmation=False,
                    )
                )

    if isinstance(request, Request) and request.mode == "shell":
        shell_line = (request.shell_line or "").strip()
        if shell_line.startswith(
            (
                "genomes prepare ",
                "helpers prepare ",
                "cutrun prepare ",
                "resources sync-",
            )
        ):
            actions.append(
                _suggested_action(
                    label="Re-check services status",
                    kind="refresh_status",
                    shell_line="services status",
                    timeout_secs=180,
                    rationale=(
                        "After a prepare or resource-sync step, refresh the combined readiness "
                        "view to confirm the installation state."
                    ),
                    requires_confirmation=False,
                )
            )

    if not actions:
        return None

    unique: list[dict[str, Any]] = []
    seen: set[tuple[str, str]] = set()
    for action in actions:
        key = (str(action.get("kind", "")), str(action.get("shell_line", "")))
        if key in seen:
            continue
        seen.add(key)
        unique.append(action)
    return unique


def _merge_suggested_actions(
    base_actions: list[dict[str, Any]] | None,
    extra_actions: list[dict[str, Any]] | None,
) -> list[dict[str, Any]] | None:
    if not extra_actions:
        return base_actions

    def merge_key(action: dict[str, Any]) -> tuple[str, str, tuple[str, ...]]:
        expected = tuple(
            str(path)
            for path in action.get("expected_artifacts", [])
            if isinstance(path, str)
        )
        kind = str(action.get("kind") or "").strip()
        shell_line = str(action.get("shell_line") or "").strip()
        return (kind, shell_line, expected if kind == "continue_artifact" else ())

    merged = list(base_actions or [])
    seen = {merge_key(action) for action in merged if isinstance(action, dict)}
    for action in extra_actions:
        if not isinstance(action, dict):
            continue
        key = merge_key(action)
        if key in seen:
            continue
        seen.add(key)
        merged.append(action)
    return merged or None


def _stamp_action_envelope(
    actions: list[dict[str, Any]] | None,
) -> list[dict[str, Any]] | None:
    if not actions:
        return actions
    for action in actions:
        if not isinstance(action, dict):
            continue
        action.setdefault("skill_alias", SKILL_NAME)
        request_payload = action.get("request")
        if (
            isinstance(request_payload, dict)
            and request_payload.get("confirm") is True
            and not action.get("requires_confirmation", False)
        ):
            action["requires_confirmation"] = True
    return actions


def _stamp_blocked_action_envelopes(
    blocked_actions: list[dict[str, Any]] | None,
) -> list[dict[str, Any]] | None:
    if not blocked_actions:
        return blocked_actions
    for blocked_action in blocked_actions:
        if not isinstance(blocked_action, dict):
            continue
        action = blocked_action.get("action")
        if isinstance(action, dict):
            _stamp_action_envelope([action])
    return blocked_actions


def _normalize_ui_intent_catalog_row(row: Any) -> dict[str, Any] | None:
    if not isinstance(row, dict):
        return None
    target = str(row.get("target") or "").strip()
    title = str(row.get("title") or "").strip()
    if not target or not title:
        return None
    return {
        "target": target,
        "title": title,
        "detail": str(row.get("detail") or "").strip(),
        "keywords": str(row.get("keywords") or "").strip(),
        "menu_path": str(row.get("menu_path") or "").strip(),
        "actions": _string_list(row.get("actions")),
        "optional_arguments": _string_list(row.get("optional_arguments")),
    }


def _ui_intent_catalog_target_rows(
    ui_intent_catalog: Any,
) -> list[dict[str, Any]]:
    if not isinstance(ui_intent_catalog, dict):
        return []
    rows = ui_intent_catalog.get("target_details")
    if not isinstance(rows, list):
        return []
    return [
        normalized
        for normalized in (
            _normalize_ui_intent_catalog_row(row) for row in rows
        )
        if normalized is not None
    ]


def _ui_intent_catalog_suggested_actions(
    ui_intent_catalog: Any,
) -> list[dict[str, Any]] | None:
    actions: list[dict[str, Any]] = []
    for row in _ui_intent_catalog_target_rows(ui_intent_catalog):
        supported_actions = row.get("actions") or []
        if "open" not in supported_actions:
            continue
        title = str(row.get("title") or row["target"]).strip()
        detail = str(row.get("detail") or "").strip()
        menu_path = str(row.get("menu_path") or "").strip()
        rationale = detail.rstrip(".") if detail else f"Open the shared `{row['target']}` UI intent"
        if menu_path:
            rationale = f"{rationale}. Operator handoff target under the {menu_path} menu."
        else:
            rationale = f"{rationale}. Operator handoff target from the shared UI-intent catalog."
        actions.append(
            _suggested_action(
                label=f"Open {title}",
                kind="ui_intent",
                shell_line=f"ui open {row['target']}",
                timeout_secs=180,
                rationale=rationale,
                requires_confirmation=False,
                ui_intent=_normalize_ui_intent_metadata(
                    {
                        "action": "open",
                        "target": row["target"],
                        "title": title,
                        "detail": detail,
                        "keywords": row.get("keywords"),
                        "menu_path": menu_path,
                        "optional_arguments": row.get("optional_arguments"),
                    }
                ),
            )
        )
    return actions or None


def _is_default_demo_request(request: Request | None) -> bool:
    if request is None:
        return False
    return (
        request.mode == "shell"
        and (request.shell_line or "").strip() == DEFAULT_DEMO_SHELL_LINE
        and request.timeout_secs == 180
        and request.expected_artifacts == [DEFAULT_DEMO_EXPECTED_ARTIFACT]
    )


def _artifact_chat_summary_lines(
    request: Request | None,
    preferred_artifacts: list[dict[str, Any]] | None,
) -> list[str] | None:
    if not preferred_artifacts:
        return None
    best_first = next(
        (
            artifact
            for artifact in preferred_artifacts
            if artifact.get("is_best_first_artifact") is True
        ),
        preferred_artifacts[0],
    )
    best_first_path = str(best_first.get("path") or "").strip()
    if not best_first_path:
        return None
    workflow_path = str(request.workflow_path or "") if request else ""
    if "simple_pcr_primer_design_offline" in workflow_path:
        return [
            "Generated a simple PCR explanation figure and primer-design report.",
            (
                "The figure shows the selected core ROI, primer windows, chosen "
                "primers, and final amplicon; the JSON report contains the ranked "
                "primer-pair details."
            ),
            f"Best-first preview artifact: {best_first_path}",
        ]
    if not _is_default_demo_request(request):
        caption = str(best_first.get("caption") or "GENtle figure").strip()
        return [
            f"Generated GENtle display artifact: {caption}.",
            "The first preferred artifact is suitable for chat or web display.",
            f"Best-first preview artifact: {best_first_path}",
        ]
    return [
        "Generated a deterministic GENtle protocol cartoon for a two-fragment Gibson assembly.",
        (
            "The ClawBio demo now starts with a graphical export so the first reply can "
            "show an actual figure instead of only listing commands."
        ),
        f"Best-first preview artifact: {best_first_path}",
    ]


def _artifact_bundle_summary(
    collected_artifacts: list[dict[str, Any]],
    preferred_artifacts: list[dict[str, Any]] | None,
    suggested_actions: list[dict[str, Any]] | None,
) -> dict[str, Any] | None:
    if not collected_artifacts and not preferred_artifacts:
        return None
    preferred = [
        dict(artifact)
        for artifact in (preferred_artifacts or [])
        if isinstance(artifact, dict)
    ]
    best_first = next(
        (
            artifact
            for artifact in preferred
            if artifact.get("is_best_first_artifact") is True
        ),
        preferred[0] if preferred else None,
    )
    displayable_suffixes = (".png", ".svg")
    displayable_artifacts = [
        artifact
        for artifact in collected_artifacts
        if any(
            str(artifact.get(key) or "").lower().endswith(displayable_suffixes)
            for key in ("bundle_path", "copied_path", "declared_path")
        )
    ]
    continuation_actions = [
        action
        for action in (suggested_actions or [])
        if isinstance(action, dict) and action.get("kind") == "continue_artifact"
    ]
    best_first_path = (
        str(best_first.get("path") or "").strip()
        if isinstance(best_first, dict)
        else ""
    )
    summary_lines: list[str] = []
    if best_first_path:
        summary_lines.append(f"Best-first artifact: {best_first_path}")
    if displayable_artifacts:
        summary_lines.append(
            f"Displayable artifacts in bundle: {len(displayable_artifacts)}"
        )
    if continuation_actions:
        summary_lines.append(
            "Continuation actions available for "
            f"{len(continuation_actions)} additional figure(s)."
        )
    return {
        "schema": "gentle.clawbio_artifact_bundle_summary.v1",
        "best_first_artifact": best_first,
        "preferred_artifact_count": len(preferred),
        "displayable_artifact_count": len(displayable_artifacts),
        "collected_artifact_count": len(collected_artifacts),
        "continuation_action_count": len(continuation_actions),
        "summary_lines": summary_lines,
    }


def _suggested_action_command_text(action: dict[str, Any]) -> str:
    shell_line = str(action.get("shell_line") or "").strip()
    if shell_line:
        return shell_line
    request_payload = action.get("request")
    if isinstance(request_payload, dict):
        expected_artifacts = request_payload.get("expected_artifacts")
        if (
            action.get("kind") == "continue_artifact"
            and isinstance(expected_artifacts, list)
            and expected_artifacts
        ):
            return (
                "use nested request payload for "
                + ", ".join(str(path) for path in expected_artifacts)
            )
        return "use nested request payload"
    return "(unknown)"


def _ensure_default_demo_suggested_action(
    request: Request | None,
    suggested_actions: list[dict[str, Any]] | None,
) -> list[dict[str, Any]] | None:
    if not _is_default_demo_request(request):
        return suggested_actions
    actions = list(suggested_actions or [])
    if any(str(action.get("shell_line", "")).strip() == "capabilities" for action in actions):
        return actions
    actions.append(
        _suggested_action(
            label="Learn GENtle capabilities",
            kind="learn_capabilities",
            shell_line="capabilities",
            timeout_secs=180,
            rationale=(
                "The demo intentionally starts with one graphical export. Run "
                "`capabilities` next to inspect the broader deterministic GENtle "
                "CLI and engine surface."
            ),
            requires_confirmation=False,
        )
    )
    return actions


def _probe_ui_intent_catalog(
    request: Request,
    resolution: CliResolution,
    execution_cwd: Path,
    script_path: Path,
) -> tuple[dict[str, Any] | None, dict[str, Any] | None, str | None]:
    probe_request = Request(
        mode="shell",
        timeout_secs=min(request.timeout_secs, 180),
        state_path=request.state_path,
        shell_line=UI_INTENT_DISCOVERY_SHELL_LINE,
    )
    probe_result, probe_step = _run_cli_command(
        resolution,
        _build_cli_args(probe_request, script_path),
        execution_cwd,
        probe_request.timeout_secs,
    )
    failure_summary = _build_failure_summary(
        stage="ui_intent_catalog_probe",
        step=probe_step,
        execution_cwd=execution_cwd,
    )
    if probe_result.returncode != 0:
        return (
            None,
            probe_step,
            _build_failure_message(
                headline=(
                    "best-effort `ui intents` discovery probe failed; "
                    "UI-intent handoff suggestions were skipped."
                ),
                failure_summary=failure_summary,
            ),
        )
    payload = _parse_stdout_json(probe_result.stdout)
    if not isinstance(payload, dict):
        return (
            None,
            probe_step,
            _build_failure_message(
                headline=(
                    "best-effort `ui intents` discovery probe did not return valid "
                    "JSON; UI-intent handoff suggestions were skipped."
                ),
                failure_summary=failure_summary,
            ),
        )
    schema = str(payload.get("schema") or "").strip()
    if schema != UI_INTENT_CATALOG_SCHEMA:
        return (
            None,
            probe_step,
            _build_failure_message(
                headline=(
                    "best-effort `ui intents` discovery probe returned unexpected "
                    f"schema `{schema or '(missing)'}`; UI-intent handoff suggestions "
                    "were skipped."
                ),
                failure_summary=failure_summary,
            ),
        )
    return payload, probe_step, None


def _write_report(
    path: Path,
    request: Request | None,
    resolution: CliResolution | None,
    execution_cwd: Path,
    command: list[str] | None,
    run_result: subprocess.CompletedProcess[str] | None,
    stdout_json: Any | None,
    chat_summary_lines: list[str] | None,
    claim_ledger: dict[str, Any] | None,
    collected_artifacts: list[dict[str, Any]],
    reference_preflight: dict[str, Any] | None,
    started_utc: str,
    ended_utc: str,
    status: str,
    error_message: str | None,
    failure_summary: dict[str, Any] | None,
    artifact_summary: dict[str, Any] | None,
    preferred_artifacts: list[dict[str, Any]] | None,
    suggested_actions: list[dict[str, Any]] | None,
    preferred_demo_actions: list[dict[str, Any]] | None,
    blocked_actions: list[dict[str, Any]] | None,
    ui_intent_catalog: dict[str, Any] | None,
    ui_intent_catalog_error: str | None,
    external_primer_handoff: dict[str, Any] | None,
    warnings: list[str] | None,
) -> None:
    command_text = _format_command_text(command)
    stdout = run_result.stdout if run_result else ""
    stderr = run_result.stderr if run_result else ""
    exit_code = run_result.returncode if run_result else None
    stdout_preview = _one_line_preview(stdout)
    stderr_preview = _one_line_preview(stderr)
    lines = [
        "# GENtle ClawBio Skill Report",
        "",
        f"- Started (UTC): `{started_utc}`",
        f"- Ended (UTC): `{ended_utc}`",
        f"- Status: `{status}`",
        f"- Mode: `{request.mode if request is not None else 'unknown'}`",
        f"- Invocation marker: `{INVOCATION_MARKER}`",
    ]
    if request is not None and request.state_path:
        lines.append(f"- State path: `{request.state_path}`")
    if resolution is not None:
        lines.append(f"- Resolver: `{resolution.label}`")
    lines.append(f"- Execution cwd: `{execution_cwd}`")
    if exit_code is not None:
        lines.append(f"- Exit code: `{exit_code}`")
    lines.append(f"- Command text: `{command_text}`")
    lines.append(f"- Stdout preview: `{stdout_preview or '(empty)'}`")
    lines.append(f"- Stderr preview: `{stderr_preview or '(empty)'}`")
    if isinstance(stdout_json, dict) and isinstance(stdout_json.get("schema"), str):
        lines.append(f"- Parsed stdout JSON schema: `{stdout_json['schema']}`")
    if isinstance(ui_intent_catalog, dict):
        target_rows = ui_intent_catalog.get("target_details")
        target_count = len(target_rows) if isinstance(target_rows, list) else 0
        lines.append(f"- UI intent catalog targets: `{target_count}`")
    if ui_intent_catalog_error:
        lines.append(f"- UI intent catalog probe: `{ui_intent_catalog_error}`")
    if warnings:
        lines.append(f"- Warning count: `{len(warnings)}`")
    if error_message:
        lines.append(f"- Error: `{error_message}`")
    if failure_summary:
        lines.append(
            f"- Failure stage: `{failure_summary.get('stage', 'unknown')}`"
        )
    if collected_artifacts:
        lines.append("- Collected artifacts:")
        for artifact in collected_artifacts:
            line = (
                f"  - `{artifact['declared_path']}` -> `{artifact['bundle_path']}` "
                f"(`{artifact['copied_path']}`)"
            )
            if artifact.get("derived_from"):
                line += f" derived from `{artifact['derived_from']}`"
            lines.append(line)
    if preferred_artifacts:
        best_first = next(
            (
                artifact
                for artifact in preferred_artifacts
                if artifact.get("is_best_first_artifact") is True
            ),
            preferred_artifacts[0],
        )
        lines.append(
            f"- Best first artifact: `{best_first.get('path', '(unknown)')}`"
        )
    if reference_preflight:
        lines.append(
            f"- Reference preflight: `{reference_preflight.get('status', 'unknown')}`"
        )
    lines.extend(
        [
            "",
            "## Execution Summary",
            "",
            f"- Command: `{command_text}`",
            "- Content-bound execution receipt: `reproducibility/execution_manifest.json`",
        ]
    )
    if exit_code is not None:
        lines.append(f"- Exit code: `{exit_code}`")
    lines.append(f"- Stdout preview: `{stdout_preview or '(empty)'}`")
    lines.append(f"- Stderr preview: `{stderr_preview or '(empty)'}`")
    if chat_summary_lines:
        lines.extend(["", "## Chat Summary", ""])
        lines.extend(f"- {line}" for line in chat_summary_lines)
    if external_primer_handoff:
        report_binding = external_primer_handoff.get("report_binding")
        lines.extend(
            [
                "",
                "## External Primer Handoff",
                "",
                f"- Status: `{external_primer_handoff.get('status', 'unknown')}`",
                "- Canonical request SHA-256: "
                f"`{external_primer_handoff.get('canonical_request_sha256', '')}`",
                "- Submitted records: "
                f"`{external_primer_handoff.get('submitted_record_count', 0)}`",
                "- Not submitted records: "
                f"`{external_primer_handoff.get('not_submitted_record_count', 0)}`",
            ]
        )
        if isinstance(report_binding, dict):
            lines.append(
                f"- GENtle report: `{report_binding.get('report_id', '')}`"
            )
            lines.append(
                "- Normalized batch SHA-256: "
                f"`{report_binding.get('normalized_batch_sha256', '')}`"
            )
        if external_primer_handoff.get("execution_binding_sha256"):
            lines.append(
                "- Execution binding SHA-256: "
                f"`{external_primer_handoff['execution_binding_sha256']}`"
            )
        lines.append(
            "- Scientific artifact hashes: "
            f"`{len(external_primer_handoff.get('scientific_artifacts', []))}`"
        )
    if claim_ledger:
        lines.extend(
            [
                "",
                "## Claim Attribution",
                "",
                "- Strict source prefixes are active for this run.",
                "- `[gentle]`: copied or deterministically projected from the GENtle executable JSON result.",
                "- `[clawbio]`: orchestration or presentation only; not a biological finding.",
                "- `[input]`: caller-supplied assumption or path; not a validated finding.",
                "- Direct external-tool claims require `[external:<tool>]` and explicit tool/input/output provenance.",
                "- Machine-readable ledger: `claim_ledger.json`.",
            ]
        )
        if request is not None and request.input_claims:
            lines.extend(["", "### Input Assumptions", ""])
            lines.extend(
                f"- {_prefix_claim_line(text, 'caller_input')}"
                for text in request.input_claims
            )
    if warnings:
        lines.extend(["", "## Warnings", ""])
        lines.extend(f"- {line}" for line in warnings)
    if preferred_artifacts:
        lines.extend(["", "## Preferred Artifacts", ""])
        for artifact in preferred_artifacts:
            marker = " (best first)" if artifact.get("is_best_first_artifact") is True else ""
            lines.append(
                f"- `{artifact.get('artifact_id', 'artifact')}`{marker}: "
                f"`{artifact.get('path', '(unknown)')}`"
            )
            if artifact.get("caption"):
                lines.append(f"  Caption: `{artifact['caption']}`")
            if artifact.get("recommended_use"):
                lines.append(f"  Recommended use: `{artifact['recommended_use']}`")
    if artifact_summary:
        lines.extend(["", "## Artifact Bundle Summary", ""])
        for summary_line in artifact_summary.get("summary_lines", []):
            lines.append(f"- {summary_line}")
        lines.append(
            "- Displayable artifacts: "
            f"`{artifact_summary.get('displayable_artifact_count', 0)}`"
        )
        lines.append(
            "- Continuation actions: "
            f"`{artifact_summary.get('continuation_action_count', 0)}`"
        )
    if suggested_actions:
        if len(suggested_actions) == 1:
            action = suggested_actions[0]
            lines.extend(["", "## Suggested Next Step", ""])
            lines.append(f"### {action.get('label', 'Action')}")
            lines.extend(
                [
                    "",
                    "```bash",
                    _suggested_action_command_text(action),
                    "```",
                ]
            )
            if action.get("rationale"):
                lines.extend(["", f"Why: `{action['rationale']}`"])
            if action.get("requires_confirmation") is not None:
                lines.append(
                    "Confirmation: `required`"
                    if action["requires_confirmation"]
                    else "Confirmation: `not required`"
                )
        else:
            lines.extend(["", "## Suggested Actions", ""])
            for action in suggested_actions:
                lines.append(
                    f"- `{action.get('label', 'Action')}`: `{_suggested_action_command_text(action)}`"
                )
                if action.get("rationale"):
                    lines.append(f"  Why: `{action['rationale']}`")
                if action.get("requires_confirmation") is not None:
                    lines.append(
                        "  Confirmation: `required`"
                        if action["requires_confirmation"]
                        else "  Confirmation: `not required`"
                    )
    if preferred_demo_actions:
        lines.extend(["", "## Preferred Demo Actions", ""])
        for action in preferred_demo_actions:
            lines.append(
                f"- `{action.get('label', 'Action')}`: `{action.get('shell_line', '(unknown)')}`"
            )
            if action.get("rationale"):
                lines.append(f"  Why: `{action['rationale']}`")
    if blocked_actions:
        lines.extend(["", "## Blocked Actions", ""])
        for blocked in blocked_actions:
            action = blocked.get("action") if isinstance(blocked, dict) else None
            label = (
                action.get("label", "Action")
                if isinstance(action, dict)
                else "Action"
            )
            shell_line = (
                action.get("shell_line", "(unknown)")
                if isinstance(action, dict)
                else "(unknown)"
            )
            lines.append(f"- `{label}`: `{shell_line}`")
            if blocked.get("blocked_reason"):
                lines.append(f"  Blocked reason: `{blocked['blocked_reason']}`")
            if blocked.get("unblock_hint"):
                lines.append(f"  Unblock hint: `{blocked['unblock_hint']}`")
            if blocked.get("download_url"):
                lines.append(f"  Download URL: `{blocked['download_url']}`")
    if failure_summary:
        lines.extend(
            [
                "",
                "## Failure Summary",
                "",
                f"- Stage: `{failure_summary.get('stage', 'unknown')}`",
                f"- Command: `{failure_summary.get('command_text', '(none)')}`",
            ]
        )
        if failure_summary.get("execution_cwd"):
            lines.append(
                f"- Cwd: `{failure_summary['execution_cwd']}`"
            )
        if failure_summary.get("exit_code") is not None:
            lines.append(
                f"- Exit code: `{failure_summary['exit_code']}`"
            )
        if failure_summary.get("note"):
            lines.append(f"- Note: `{failure_summary['note']}`")
        if failure_summary.get("stderr_preview"):
            lines.extend(
                [
                    "",
                    "### Stderr Preview",
                    "",
                    "```text",
                    str(failure_summary["stderr_preview"]),
                    "```",
                ]
            )
        elif failure_summary.get("stdout_preview"):
            lines.extend(
                [
                    "",
                    "### Stdout Preview",
                    "",
                    "```text",
                    str(failure_summary["stdout_preview"]),
                    "```",
                ]
            )
    lines.extend(
        [
            "",
        ]
    )
    if reference_preflight:
        lines.extend(
            [
                "## Reference Preflight",
                "",
                f"- Genome: `{reference_preflight['genome_id']}`",
                f"- Prepared before run: `{reference_preflight['prepared_before']}`",
                f"- Prepare attempted: `{reference_preflight['prepare_attempted']}`",
                f"- Prepared after run: `{reference_preflight['prepared_after']}`",
                "",
            ]
        )
        status_before = reference_preflight.get("status_before")
        if isinstance(status_before, dict):
            lines.extend(
                [
                    "### Status Before",
                    "",
                    "```json",
                    json.dumps(status_before, indent=2, ensure_ascii=True),
                    "```",
                    "",
                ]
            )
        status_after = reference_preflight.get("status_after")
        if (
            isinstance(status_after, dict)
            and status_after != status_before
        ):
            lines.extend(
                [
                    "### Status After",
                    "",
                    "```json",
                    json.dumps(status_after, indent=2, ensure_ascii=True),
                    "```",
                    "",
                ]
            )
        lines.extend(["### Preflight Commands", ""])
        for idx, step in enumerate(reference_preflight.get("steps", []), start=1):
            step_command = " ".join(
                shlex.quote(v) for v in step.get("command", [])
            ) or "(none)"
            lines.extend(
                [
                    f"{idx}. `{step.get('status', 'unknown')}`",
                    f"   Command: `{step_command}`",
                    f"   Exit code: `{step.get('exit_code')}`",
                ]
            )
        lines.extend([""])
    lines.extend(
        [
            "## Command",
            "",
            "```bash",
            command_text,
            "```",
            "",
            "## Stdout",
            "",
            "```text",
            stdout.rstrip(),
            "```",
            "",
            "## Stderr",
            "",
            "```text",
            stderr.rstrip(),
            "```",
        ]
    )
    path.write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")


def _default_demo_request() -> Request:
    return Request(
        mode="shell",
        shell_line=DEFAULT_DEMO_SHELL_LINE,
        expected_artifacts=[DEFAULT_DEMO_EXPECTED_ARTIFACT],
        timeout_secs=180,
    )


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "GENtle ClawBio skill wrapper. Executes deterministic gentle_cli "
            "commands from structured request JSON, supports sequence-grounded "
            "follow-up to biological observations, and writes reproducibility artifacts."
        )
    )
    parser.add_argument("--input", help="Path to request JSON")
    parser.add_argument("--output", required=True, help="Output directory")
    parser.add_argument("--demo", action="store_true", help="Run built-in demo request")
    parser.add_argument(
        "--skill-info",
        action="store_true",
        help="Report ClawBio skill metadata without invoking gentle_cli",
    )
    parser.add_argument(
        "--gentle-cli",
        help="Explicit command used to invoke gentle_cli. Recommended runtime "
        "is Docker/OCI or Apptainer/Singularity via GENTLE_CLI_CMD; examples "
        "include 'gentle_cli', './gentle_apptainer_cli.sh /path/to/gentle.sif', "
        "or 'cargo run --bin gentle_cli --'.",
    )
    args = parser.parse_args()

    output_dir = Path(args.output).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    repro_dir = output_dir / "reproducibility"
    repro_dir.mkdir(parents=True, exist_ok=True)

    started = _now_utc_iso()
    request: Request | None = None
    request_source_path: Path | None = None
    run_result: subprocess.CompletedProcess[str] | None = None
    main_execution_step: dict[str, Any] | None = None
    command: list[str] | None = None
    resolution: CliResolution | None = None
    execution_cwd = Path.cwd()
    status = "failed"
    error_message: str | None = None
    failure_summary: dict[str, Any] | None = None
    collected_artifacts: list[dict[str, Any]] = []
    reference_preflight: dict[str, Any] | None = None
    stdout_json: Any | None = None
    chat_summary_lines: list[str] | None = None
    chat_summary_source_kind: str | None = None
    chat_summary_projection: dict[str, Any] | None = None
    claim_ledger: dict[str, Any] | None = None
    preferred_artifacts: list[dict[str, Any]] | None = None
    suggested_actions: list[dict[str, Any]] | None = None
    preferred_demo_actions: list[dict[str, Any]] | None = None
    blocked_actions: list[dict[str, Any]] | None = None
    artifact_summary: dict[str, Any] | None = None
    ui_intent_catalog: dict[str, Any] | None = None
    ui_intent_catalog_error: str | None = None
    external_primer_handoff_context: dict[str, Any] | None = None
    external_primer_handoff_result: dict[str, Any] | None = None
    external_primer_handoff_gentle_version: str | None = None
    gentle_runtime_version: str | None = None
    verified_delegation: dict[str, Any] | None = None
    approved_execution: ApprovedExecution | None = None
    execution_proposal: dict[str, Any] | None = None
    execution_proposal_path: Path | None = None
    execution_approval: dict[str, Any] | None = None
    backend_resolution: dict[str, Any] | None = None
    content_bound_inputs: list[dict[str, Any]] = []
    state_before = _state_content_snapshot(None, execution_cwd)
    state_after = _state_content_snapshot(None, execution_cwd)
    pre_execution_steps: list[dict[str, Any]] = []
    auxiliary_steps: list[dict[str, Any]] = []
    warnings: list[str] = []

    try:
        if args.skill_info:
            request = Request(mode="skill-info")
        elif args.demo:
            request = _default_demo_request()
        else:
            if not args.input:
                raise SkillError("--input is required unless --demo or --skill-info is used")
            request_source_path = _resolve_existing_request_file(
                args.input, Path(__file__)
            )
            payload = _read_json(request_source_path)
            if payload.get("schema") == APPROVED_EXECUTION_REQUEST_SCHEMA:
                approved_execution = _load_approved_execution(
                    payload, request_source_path, Path(__file__)
                )
                normalized_request = approved_execution.proposal[
                    "approval_basis"
                ].get("normalized_request")
                if not isinstance(normalized_request, dict):
                    raise SkillError(
                        "execution proposal contains no normalized request"
                    )
                request = _coerce_request(normalized_request)
                execution_approval = {
                    "status": "approved",
                    "proposal_path": str(approved_execution.proposal_path),
                    "proposal_file_sha256": _sha256_file_prefixed(
                        approved_execution.proposal_path
                    ),
                    "proposal_digest": approved_execution.proposal[
                        "proposal_digest"
                    ],
                    "approval": approved_execution.approval,
                    "approval_envelope_path": str(
                        approved_execution.envelope_path
                    ),
                    "approval_envelope_sha256": _sha256_file_prefixed(
                        approved_execution.envelope_path
                    ),
                    "trust_boundary": (
                        "caller attests approval identity; GENtle verifies the "
                        "approved proposal digest and execution context"
                    ),
                }
            else:
                request = _coerce_request(payload)
        warnings.extend(_request_mode_warnings(request))

        own_catalog = _read_catalog_entry(Path(__file__))
        if own_catalog and (
            own_catalog.get("name") != SKILL_NAME
            or own_catalog.get("version") != SKILL_CONTRACT_VERSION
        ):
            raise SkillError(
                "gentle-cloning catalog identity/version disagrees with its runtime contract"
            )

        verified_delegation = _verified_delegation(request, Path(__file__))

        if request.mode == "skill-info":
            content_bound_inputs = _prepare_content_bound_inputs(
                request,
                request_source_path=request_source_path,
                execution_cwd=execution_cwd,
                script_path=Path(__file__),
            )
            state_before = _state_content_snapshot(request, execution_cwd)
            stdout_json = _skill_info_payload(Path(__file__))
            chat_summary_lines = _skill_info_chat_summary_lines(stdout_json)
            chat_summary_source_kind = "clawbio_presentation"
            chat_summary_projection = _summary_projection(
                chat_summary_lines,
                source_pointer="/stdout_json",
                source_payload=stdout_json,
                projection_id="gentle.clawbio.skill_info_summary.v1",
            )
            status = "ok"
        elif request.mode == "intents":
            content_bound_inputs = _prepare_content_bound_inputs(
                request,
                request_source_path=request_source_path,
                execution_cwd=execution_cwd,
                script_path=Path(__file__),
            )
            state_before = _state_content_snapshot(request, execution_cwd)
            stdout_json = _intents_runtime_payload(Path(__file__))
            chat_summary_lines = _intents_runtime_chat_summary_lines(stdout_json)
            chat_summary_source_kind = "clawbio_presentation"
            chat_summary_projection = _summary_projection(
                chat_summary_lines,
                source_pointer="/stdout_json",
                source_payload=stdout_json,
                projection_id="gentle.clawbio.intents_summary.v1",
            )
            warnings.extend(_string_list(stdout_json.get("warnings")))
            status = "ok"
        else:
            try:
                resolution = _resolve_cli(args.gentle_cli, Path(__file__))
            except SkillError as e:
                if args.demo:
                    status = "degraded_demo"
                    error_message = str(e)
                    request = _default_demo_request()
                    raise
                raise

            execution_cwd = _resolve_execution_cwd(request, resolution, Path(__file__))
            confirmation_gated = bool(
                verified_delegation
                and verified_delegation.get("requires_confirmation")
            )
            if confirmation_gated:
                _assert_delegated_request_is_resolved(request)
                if approved_execution is None:
                    if request.state_path:
                        state_path = Path(request.state_path).expanduser()
                        if not state_path.is_absolute():
                            state_path = execution_cwd / state_path
                        request.state_path = str(state_path.resolve())
                    backend_resolution, backend_step = _pin_delegated_primer_backend(
                        request, resolution, execution_cwd
                    )
                    if backend_step is not None:
                        pre_execution_steps.append(backend_step)
                resolution = _pin_delegated_runtime(
                    resolution,
                    execution_cwd,
                    Path(__file__),
                    require_local_binary=approved_execution is not None,
                )
            elif approved_execution is not None:
                raise SkillError(
                    "approved execution envelope is only valid for a "
                    "confirmation-gated delegated route"
                )

            content_bound_inputs = _prepare_content_bound_inputs(
                request,
                request_source_path=request_source_path,
                execution_cwd=execution_cwd,
                script_path=Path(__file__),
            )
            if approved_execution is not None:
                content_bound_inputs.append(
                    _approval_control_binding(
                        approved_execution.proposal_path, "execution_proposal"
                    )
                )
                content_bound_inputs.sort(
                    key=lambda row: (row["binding_ids"], row["resolved_path"])
                )
            state_before = _state_content_snapshot(request, execution_cwd)
            if confirmation_gated and approved_execution is not None:
                _verify_approved_execution_context(
                    approved=approved_execution,
                    request=request,
                    content_bound_inputs=content_bound_inputs,
                    delegation=verified_delegation,
                    resolution=resolution,
                    gentle_version=None,
                    execution_cwd=execution_cwd,
                    state_before=state_before,
                )

            if verified_delegation is not None and request.mode != "external-primer-handoff":
                version_result, version_step = _run_cli_command(
                    resolution,
                    ["--version"],
                    execution_cwd,
                    min(request.timeout_secs, 180),
                )
                pre_execution_steps.append(version_step)
                if version_result.returncode != 0:
                    failure_summary = _build_failure_summary(
                        stage="delegated_runtime_version_preflight",
                        step=version_step,
                        execution_cwd=execution_cwd,
                    )
                    raise SkillError(
                        _build_failure_message(
                            headline=(
                                "GENtle runtime version preflight failed for the "
                                "delegated skill request."
                            ),
                            failure_summary=failure_summary,
                        )
                    )
                gentle_runtime_version = (
                    version_result.stdout.strip() or version_result.stderr.strip()
                )
                if not gentle_runtime_version:
                    raise SkillError(
                        "GENtle runtime version preflight returned no version text"
                    )

            if confirmation_gated:
                if approved_execution is None:
                    resolution = _pin_delegated_runtime(
                        resolution,
                        execution_cwd,
                        Path(__file__),
                        require_local_binary=True,
                    )
                    execution_proposal = _build_execution_proposal(
                        request=request,
                        request_source_path=request_source_path,
                        content_bound_inputs=content_bound_inputs,
                        delegation=verified_delegation,
                        resolution=resolution,
                        gentle_version=gentle_runtime_version or "",
                        execution_cwd=execution_cwd,
                        state_before=state_before,
                        backend_resolution=backend_resolution,
                    )
                    execution_proposal_path = (
                        repro_dir / "execution_proposal.json"
                    )
                    execution_proposal_path.write_text(
                        json.dumps(execution_proposal, indent=2, ensure_ascii=True)
                        + "\n",
                        encoding="utf-8",
                    )
                    execution_approval = {
                        "status": "approval_required",
                        "proposal_path": str(execution_proposal_path),
                        "proposal_file_sha256": _sha256_file_prefixed(
                            execution_proposal_path
                        ),
                        "proposal_digest": execution_proposal["proposal_digest"],
                        "trust_boundary": execution_proposal["review"][
                            "trust_boundary"
                        ],
                    }
                    stdout_json = execution_proposal
                    chat_summary_lines = [
                        "Prepared an exact GENtle execution proposal; no scientific command was run.",
                        (
                            "Approval is required for proposal digest "
                            + execution_proposal["proposal_digest"]
                            + "."
                        ),
                    ]
                    chat_summary_source_kind = "clawbio_presentation"
                    status = "approval_required"
                    raise ApprovalRequired("approval required before execution")
                backend_resolution = approved_execution.proposal[
                    "approval_basis"
                ].get("backend_resolution")
                _verify_approved_execution_context(
                    approved=approved_execution,
                    request=request,
                    content_bound_inputs=content_bound_inputs,
                    delegation=verified_delegation,
                    resolution=resolution,
                    gentle_version=gentle_runtime_version or "",
                    execution_cwd=execution_cwd,
                    state_before=state_before,
                )

            if request.ensure_reference_prepared is not None:
                reference_preflight = _reference_preflight_record(
                    request.ensure_reference_prepared
                )
                _ensure_reference_prepared(
                    request.ensure_reference_prepared,
                    resolution,
                    execution_cwd,
                    reference_preflight,
                )

            if request.mode == "external-primer-handoff":
                external_primer_handoff_context = _prepare_external_primer_handoff(
                    request, output_dir, execution_cwd
                )
                if not external_primer_handoff_context["submitted_records"]:
                    external_primer_handoff_result = (
                        _external_primer_handoff_not_run_result(
                            external_primer_handoff_context
                        )
                    )
                    chat_summary_lines = _external_primer_handoff_chat_summary_lines(
                        external_primer_handoff_result
                    )
                    status = "incomplete"
                    raise SkillError(
                        "No qPCR or endpoint-PCR pair was eligible for the GENtle "
                        "external-primer importer; no biological evaluation was run."
                    )
                version_result, version_step = _run_cli_command(
                    resolution,
                    ["--version"],
                    execution_cwd,
                    min(request.timeout_secs, 180),
                )
                pre_execution_steps.append(version_step)
                if version_result.returncode != 0:
                    failure_summary = _build_failure_summary(
                        stage="external_primer_handoff_version_preflight",
                        step=version_step,
                        execution_cwd=execution_cwd,
                    )
                    raise SkillError(
                        _build_failure_message(
                            headline=(
                                "GENtle runtime version preflight failed for the "
                                "external-primer handoff."
                            ),
                            failure_summary=failure_summary,
                        )
                    )
                external_primer_handoff_gentle_version = (
                    version_result.stdout.strip() or version_result.stderr.strip()
                )
                if not external_primer_handoff_gentle_version:
                    raise SkillError(
                        "GENtle runtime version preflight returned no version text"
                    )
                gentle_runtime_version = external_primer_handoff_gentle_version

            cli_args = _build_cli_args(request, Path(__file__))
            _prepare_expected_artifact_parent_dirs(request, execution_cwd)
            run_result, main_step = _run_cli_command(
                resolution,
                cli_args,
                execution_cwd,
                request.timeout_secs,
            )
            main_execution_step = main_step
            command = main_step["command"]
            status = "ok" if run_result.returncode == 0 else "command_failed"
            if run_result.returncode != 0:
                failure_summary = _build_failure_summary(
                    stage="main_command",
                    step=main_step,
                    execution_cwd=execution_cwd,
                )
                error_message = _build_failure_message(
                    headline="gentle_cli exited with a non-zero status.",
                    failure_summary=failure_summary,
                )
            stdout_json = _parse_stdout_json(run_result.stdout)
            if request.mode == "version" and run_result.returncode == 0:
                gentle_runtime_version = (
                    run_result.stdout.strip() or run_result.stderr.strip() or None
                )
            if (
                request.mode == "external-primer-handoff"
                and external_primer_handoff_context is not None
            ):
                if run_result.returncode == 0:
                    try:
                        external_primer_handoff_result = (
                            _verify_external_primer_handoff(
                                external_primer_handoff_context,
                                stdout_json,
                                external_primer_handoff_gentle_version or "",
                            )
                        )
                    except SkillError as verification_error:
                        external_primer_handoff_result = (
                            _external_primer_handoff_failure_result(
                                external_primer_handoff_context,
                                external_primer_handoff_gentle_version,
                                str(verification_error),
                            )
                        )
                        status = "verification_failed"
                        error_message = str(verification_error)
                        failure_summary = {
                            "stage": "external_primer_handoff_verification",
                            "note": str(verification_error),
                            "command": command,
                            "command_text": _format_command_text(command),
                            "execution_cwd": str(execution_cwd),
                            "exit_code": run_result.returncode,
                            "stderr_preview": _one_line_preview(run_result.stderr),
                            "stdout_preview": _one_line_preview(run_result.stdout),
                        }
                    else:
                        if external_primer_handoff_result["status"] == "incomplete":
                            status = "incomplete"
                            error_message = (
                                "GENtle returned an incomplete requested specificity "
                                "evaluation; the verified partial result was retained."
                            )
                    chat_summary_lines = _external_primer_handoff_chat_summary_lines(
                        external_primer_handoff_result
                    )
                    chat_summary_source_kind = "clawbio_presentation"
                    chat_summary_projection = _summary_projection(
                        chat_summary_lines,
                        source_pointer="/external_primer_handoff",
                        source_payload=external_primer_handoff_result,
                        projection_id="gentle.clawbio.external_primer_handoff_summary.v1",
                    )
                else:
                    external_primer_handoff_result = (
                        _external_primer_handoff_base_result(
                            external_primer_handoff_context,
                            status="command_failed",
                            gentle_version=external_primer_handoff_gentle_version,
                            state_sha256_after=_sha256_file_prefixed(
                                Path(external_primer_handoff_context["state_path"])
                            ),
                        )
                    )
                    external_primer_handoff_result["analysis_completeness"] = (
                        "not_computed"
                    )
                    external_primer_handoff_result["warnings"] = [
                        error_message or "GENtle external-primer import failed."
                    ]
                    chat_summary_lines = _external_primer_handoff_chat_summary_lines(
                        external_primer_handoff_result
                    )
                    chat_summary_source_kind = "clawbio_presentation"
                    chat_summary_projection = _summary_projection(
                        chat_summary_lines,
                        source_pointer="/external_primer_handoff",
                        source_payload=external_primer_handoff_result,
                        projection_id="gentle.clawbio.external_primer_handoff_summary.v1",
                    )
            if request.mode != "external-primer-handoff":
                chat_summary_projection = _extract_chat_summary_projection(stdout_json)
                if chat_summary_projection is not None:
                    chat_summary_lines = chat_summary_projection["lines"]
                    chat_summary_source_kind = "gentle_executable"
            if chat_summary_lines is None:
                chat_summary_lines = _runtime_version_chat_summary_lines(
                    request,
                    run_result,
                )
                if chat_summary_lines is not None:
                    chat_summary_source_kind = "clawbio_presentation"
            preferred_artifacts = _extract_preferred_artifacts(stdout_json)
            suggested_actions = _extract_suggested_actions(stdout_json, request)
            preferred_demo_actions = _extract_preferred_demo_actions(stdout_json)
            blocked_actions = _extract_blocked_actions(stdout_json)
            if request.mode == "capabilities" and run_result.returncode == 0:
                (
                    ui_intent_catalog,
                    ui_intent_step,
                    ui_intent_catalog_error,
                ) = _probe_ui_intent_catalog(
                    request,
                    resolution,
                    execution_cwd,
                    Path(__file__),
                )
                if ui_intent_step is not None:
                    auxiliary_steps.append(ui_intent_step)
                suggested_actions = _merge_suggested_actions(
                    suggested_actions,
                    _ui_intent_catalog_suggested_actions(ui_intent_catalog),
                )
            collected_artifacts = _copy_collected_artifacts(
                request, output_dir, execution_cwd
            )
            collected_artifacts, preferred_artifacts = _augment_artifacts_with_storyboard(
                request,
                output_dir,
                collected_artifacts,
                preferred_artifacts,
            )
            rasterize_declared_paths = _best_first_svg_declared_paths(
                collected_artifacts,
                preferred_artifacts,
            )
            collected_artifacts, rasterized_pngs = _rasterize_collected_svg_artifacts(
                collected_artifacts,
                resolution,
                execution_cwd,
                output_dir,
                rasterize_declared_paths,
            )
            preferred_artifacts = _rewrite_preferred_artifacts_for_png(
                collected_artifacts,
                preferred_artifacts,
                rasterized_pngs,
            )
            continue_artifact_actions = _continue_artifact_suggested_actions(
                request,
                collected_artifacts,
                preferred_artifacts,
            )
            suggested_actions = _merge_suggested_actions(
                suggested_actions,
                continue_artifact_actions,
            )
            suggested_actions = _ensure_default_demo_suggested_action(
                request,
                suggested_actions,
            )
            if chat_summary_lines is None:
                chat_summary_lines = _artifact_chat_summary_lines(
                    request,
                    preferred_artifacts,
                )
                if chat_summary_lines is not None:
                    chat_summary_source_kind = "clawbio_presentation"
                    chat_summary_projection = _summary_projection(
                        chat_summary_lines,
                        source_pointer="/preferred_artifacts",
                        source_payload=preferred_artifacts,
                        projection_id="gentle.clawbio.artifact_summary.v1",
                    )
            if chat_summary_lines is None:
                chat_summary_lines = _fallback_chat_summary_lines(
                    request=request,
                    command=command,
                    run_result=run_result,
                    stdout_json=stdout_json,
                    status=status,
                )
                if chat_summary_lines is not None:
                    chat_summary_source_kind = "clawbio_presentation"
                    chat_summary_projection = _summary_projection(
                        chat_summary_lines,
                        source_pointer=(
                            "/stdout_json"
                            if stdout_json is not None
                            else "/stdout"
                        ),
                        source_payload=(
                            stdout_json
                            if stdout_json is not None
                            else (run_result.stdout if run_result else "")
                        ),
                        projection_id="gentle.clawbio.fallback_summary.v1",
                    )
            chat_summary_lines = _append_continue_artifact_notice(
                chat_summary_lines,
                continue_artifact_actions,
            )
            artifact_summary = _artifact_bundle_summary(
                collected_artifacts,
                preferred_artifacts,
                suggested_actions,
            )
    except subprocess.TimeoutExpired as e:
        request = request if request is not None else _default_demo_request()
        error_message = f"command timed out after {e.timeout} seconds"
        status = "timeout"
    except ApprovalRequired:
        error_message = None
        status = "approval_required"
    except SkillError as e:
        if (
            failure_summary is None
            and reference_preflight is not None
            and isinstance(reference_preflight.get("last_failure"), dict)
        ):
            failure_summary = reference_preflight["last_failure"]
        error_message = str(e)
        if status not in {"degraded_demo", "incomplete", "verification_failed"}:
            status = "failed"
    except Exception as e:  # pragma: no cover - defensive boundary
        request = request if request is not None else _default_demo_request()
        error_message = f"unexpected error: {type(e).__name__}: {e}"
        status = "failed"

    try:
        state_after = _state_content_snapshot(request, execution_cwd)
    except SkillError as e:
        state_after = {
            "declared": bool(request and request.state_path),
            "path": request.state_path if request else None,
            "exists": None,
            "sha256": None,
            "error": str(e),
        }
        status = "failed"
        error_message = str(e)

    try:
        _verify_content_bound_inputs_unchanged(content_bound_inputs)
    except SkillError as e:
        provenance_error = str(e)
        error_message = (
            f"{error_message}\nProvenance verification failed: {provenance_error}"
            if error_message
            else provenance_error
        )
        status = "verification_failed"
        failure_summary = {
            "stage": "content_binding_verification",
            "note": provenance_error,
            "command": command,
            "command_text": _format_command_text(command),
            "execution_cwd": str(execution_cwd),
            "exit_code": run_result.returncode if run_result else None,
            "stderr_preview": _one_line_preview(run_result.stderr) if run_result else "",
            "stdout_preview": _one_line_preview(run_result.stdout) if run_result else "",
        }

    if claim_ledger is None:
        chat_summary_lines, claim_ledger = _strict_claim_attribution(
            request=request,
            lines=chat_summary_lines,
            default_source_kind=chat_summary_source_kind,
            summary_projection=chat_summary_projection,
            stdout_json=stdout_json,
            command=command,
            collected_artifacts=collected_artifacts,
            warning_lines=warnings,
        )

    ended = _now_utc_iso()

    report_path = output_dir / "report.md"
    result_path = output_dir / "result.json"
    claim_ledger_path = output_dir / "claim_ledger.json"
    commands_path = repro_dir / "commands.sh"
    env_path = repro_dir / "environment.yml"
    execution_manifest_path = repro_dir / "execution_manifest.json"
    checksums_path = repro_dir / "checksums.sha256"

    suggested_actions = _stamp_action_envelope(suggested_actions)
    preferred_demo_actions = _stamp_action_envelope(preferred_demo_actions)
    blocked_actions = _stamp_blocked_action_envelopes(blocked_actions)

    if claim_ledger is not None:
        claim_ledger["execution_manifest_schema"] = EXECUTION_MANIFEST_SCHEMA
        claim_ledger_path.write_text(
            json.dumps(claim_ledger, indent=2, ensure_ascii=True) + "\n",
            encoding="utf-8",
        )

    _write_report(
        path=report_path,
        request=request,
        resolution=resolution,
        execution_cwd=execution_cwd,
        command=command,
        run_result=run_result,
        stdout_json=stdout_json,
        chat_summary_lines=chat_summary_lines,
        claim_ledger=claim_ledger,
        collected_artifacts=collected_artifacts,
        reference_preflight=reference_preflight,
        started_utc=started,
        ended_utc=ended,
        status=status,
        error_message=error_message,
        failure_summary=failure_summary,
        artifact_summary=artifact_summary,
        preferred_artifacts=preferred_artifacts,
        suggested_actions=suggested_actions,
        preferred_demo_actions=preferred_demo_actions,
        blocked_actions=blocked_actions,
        ui_intent_catalog=ui_intent_catalog,
        ui_intent_catalog_error=ui_intent_catalog_error,
        external_primer_handoff=external_primer_handoff_result,
        warnings=warnings,
    )

    command_lines = [
        " ".join(shlex.quote(v) for v in step.get("command", []))
        for step in (reference_preflight or {}).get("steps", [])
        if step.get("command")
    ]
    command_lines.extend(
        " ".join(shlex.quote(v) for v in step.get("command", []))
        for step in pre_execution_steps
        if step.get("command")
    )
    if command:
        command_lines.append(" ".join(shlex.quote(v) for v in command))
    command_lines.extend(
        " ".join(shlex.quote(v) for v in step.get("command", []))
        for step in auxiliary_steps
        if step.get("command")
    )
    if not command_lines:
        command_lines = ["# no command executed"]
    commands_text = "\n".join(
        [
            "#!/usr/bin/env bash",
            "set -euo pipefail",
            f"# generated_utc: {ended}",
            f"cd {shlex.quote(str(execution_cwd))}",
            *command_lines,
            "",
        ]
    )
    commands_path.write_text(commands_text, encoding="utf-8")
    os.chmod(commands_path, 0o755)

    _write_repro_environment(env_path)

    manifest_steps: list[dict[str, Any]] = []
    for step in (reference_preflight or {}).get("steps", []):
        manifest_steps.append(dict(step, manifest_stage="reference_preflight"))
    manifest_steps.extend(
        dict(step, manifest_stage="runtime_preflight")
        for step in pre_execution_steps
    )
    if main_execution_step is not None:
        manifest_steps.append(dict(main_execution_step, manifest_stage="main_command"))
    manifest_steps.extend(
        dict(step, manifest_stage="auxiliary") for step in auxiliary_steps
    )
    manifest_artifacts = [
        _content_artifact_record(report_path, "human_report", output_dir),
        _content_artifact_record(commands_path, "replay_commands", output_dir),
        _content_artifact_record(env_path, "runtime_environment", output_dir),
    ]
    if execution_proposal_path is not None and execution_proposal_path.is_file():
        manifest_artifacts.append(
            _content_artifact_record(
                execution_proposal_path, "execution_proposal", output_dir
            )
        )
    manifest_artifacts.extend(
        _content_artifact_record(
            Path(artifact["copied_path"]), "collected_output", output_dir
        )
        for artifact in collected_artifacts
        if isinstance(artifact, dict)
        and artifact.get("copied_path")
        and Path(artifact["copied_path"]).is_file()
    )
    if external_primer_handoff_result is not None:
        manifest_artifacts.extend(
            _content_artifact_record(
                Path(artifact["path"]), "native_scientific_output", output_dir
            )
            for artifact in external_primer_handoff_result.get(
                "scientific_artifacts", []
            )
            if isinstance(artifact, dict)
            and artifact.get("path")
            and Path(artifact["path"]).is_file()
        )
    execution_manifest = _build_execution_manifest(
        request=request,
        request_source_path=request_source_path,
        content_bound_inputs=content_bound_inputs,
        delegation=verified_delegation,
        execution_approval=execution_approval,
        script_path=Path(__file__).resolve(),
        resolution=resolution,
        gentle_version=gentle_runtime_version,
        execution_cwd=execution_cwd,
        started_utc=started,
        ended_utc=ended,
        status=status,
        error_message=error_message,
        state_before=state_before,
        state_after=state_after,
        reference_preflight=reference_preflight,
        execution_steps=manifest_steps,
        stdout=(run_result.stdout if run_result else ""),
        stderr=(run_result.stderr if run_result else ""),
        stdout_json=stdout_json,
        artifacts=manifest_artifacts,
        claim_ledger_path=(claim_ledger_path if claim_ledger is not None else None),
    )
    execution_manifest_path.write_text(
        json.dumps(execution_manifest, indent=2, ensure_ascii=True) + "\n",
        encoding="utf-8",
    )

    result_payload = {
        "schema": RESULT_SCHEMA,
        "invocation_marker": INVOCATION_MARKER,
        "status": status,
        "request": (dataclasses.asdict(request) if request is not None else None),
        "started_utc": started,
        "ended_utc": ended,
        "resolver": (dataclasses.asdict(resolution) if resolution else None),
        "command": command,
        "exit_code": (run_result.returncode if run_result else None),
        "stdout": (run_result.stdout if run_result else ""),
        "stdout_json": stdout_json,
        "stderr": (run_result.stderr if run_result else ""),
        "chat_summary_lines": chat_summary_lines,
        "claim_ledger": claim_ledger,
        "execution_manifest": execution_manifest,
        "execution_proposal": execution_proposal,
        "execution_approval": execution_approval,
        "artifact_summary": artifact_summary,
        "preferred_artifacts": preferred_artifacts,
        "suggested_actions": suggested_actions,
        "preferred_demo_actions": preferred_demo_actions,
        "blocked_actions": blocked_actions,
        "ui_intent_catalog": ui_intent_catalog,
        "ui_intent_catalog_error": ui_intent_catalog_error,
        "external_primer_handoff": external_primer_handoff_result,
        "warnings": warnings or None,
        "error": error_message,
        "failure_summary": failure_summary,
        "preflight": {
            "reference_preparation": reference_preflight,
        },
        "artifacts": {
            "report_md": str(report_path),
            "result_json": str(result_path),
            "claim_ledger_json": (
                str(claim_ledger_path) if claim_ledger is not None else None
            ),
            "execution_manifest_json": str(execution_manifest_path),
            "execution_proposal_json": (
                str(execution_proposal_path)
                if execution_proposal_path is not None
                else None
            ),
            "repro_commands": str(commands_path),
            "repro_environment": str(env_path),
            "repro_checksums": str(checksums_path),
            "collected": collected_artifacts,
        },
    }
    result_path.write_text(
        json.dumps(result_payload, indent=2, ensure_ascii=True) + "\n",
        encoding="utf-8",
    )

    checksum_paths = [
        report_path,
        result_path,
        commands_path,
        env_path,
        execution_manifest_path,
    ]
    checksum_paths.extend(
        Path(artifact["copied_path"])
        for artifact in collected_artifacts
        if isinstance(artifact, dict) and artifact.get("copied_path")
    )
    if external_primer_handoff_context is not None:
        checksum_paths.append(
            Path(external_primer_handoff_context["request_snapshot_path"])
        )
        if external_primer_handoff_context.get("batch_file_sha256"):
            checksum_paths.append(Path(external_primer_handoff_context["batch_path"]))
    if external_primer_handoff_result is not None:
        checksum_paths.extend(
            Path(artifact["path"])
            for artifact in external_primer_handoff_result.get(
                "scientific_artifacts", []
            )
            if isinstance(artifact, dict) and artifact.get("path")
        )
    if claim_ledger is not None:
        checksum_paths.append(claim_ledger_path)
    if execution_proposal_path is not None:
        checksum_paths.append(execution_proposal_path)
    checksums: list[tuple[str, str]] = []
    seen_checksum_paths: set[Path] = set()
    for artifact in checksum_paths:
        resolved = artifact.resolve()
        if resolved in seen_checksum_paths or not resolved.is_file():
            continue
        seen_checksum_paths.add(resolved)
        try:
            label = resolved.relative_to(output_dir).as_posix()
        except ValueError:
            label = resolved.name
        checksums.append((label, _sha256_file(resolved)))
    checksums.sort(key=lambda row: row[0])
    checksums_lines = [f"{digest}  {name}" for name, digest in checksums]
    checksums_path.write_text("\n".join(checksums_lines) + "\n", encoding="utf-8")

    print(json.dumps(result_payload, indent=2, ensure_ascii=True))

    return 0 if status in ("ok", "degraded_demo", "approval_required") else 1


if __name__ == "__main__":
    sys.exit(main())
