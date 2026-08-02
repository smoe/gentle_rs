from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path
import sys


SKILL_ROOT = Path(__file__).resolve().parents[1]
REPO_ROOT = SKILL_ROOT.parents[3]
GENERIC_ROOT = SKILL_ROOT.parent / "gentle-cloning"


def _json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _descriptor() -> dict:
    return _json(SKILL_ROOT / "INTENTS.json")


def _routes() -> dict[str, dict]:
    return {route["intent_id"]: route for route in _descriptor()["routes"]}


def _trigger_terms(descriptor: dict) -> set[str]:
    return {
        term.strip().lower()
        for route in descriptor["routes"]
        for term in route.get("trigger_terms", [])
    }


def _load_generic_wrapper():
    path = GENERIC_ROOT / "gentle_cloning.py"
    spec = importlib.util.spec_from_file_location("gentle_cloning_delegate", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_skill_metadata_and_descriptor_are_consistent() -> None:
    descriptor = _descriptor()
    catalog = _json(SKILL_ROOT / "catalog_entry.json")
    skill_md = (SKILL_ROOT / "SKILL.md").read_text(encoding="utf-8")

    assert descriptor["schema"] == "clawbio.skill_intents.v1"
    assert descriptor["skill"] == "gentle-pcr-primer-design"
    assert catalog["name"] == descriptor["skill"]
    assert catalog["cli_alias"] == descriptor["skill"]
    assert catalog["has_script"] is False
    assert catalog["execution_delegate"] == "gentle-cloning"
    assert catalog["delegate_contract"] == {
        "skill": "gentle-cloning",
        "skill_version": "0.2.0",
        "request_schema": "gentle.clawbio_skill_request.v1",
        "result_schema": "gentle.clawbio_skill_result.v1",
        "execution_manifest_schema": "gentle.clawbio_execution_manifest.v1",
    }
    assert catalog["has_tests"] is True
    assert "name: gentle-pcr-primer-design" in skill_md
    assert f"version: {catalog['version']}" in skill_md
    assert "skill: gentle-cloning" in skill_md
    assert not list(SKILL_ROOT.glob("*.py"))


def test_every_execution_delegates_to_registered_generic_runtime() -> None:
    descriptor = _descriptor()
    assert descriptor["routes"]

    for route in descriptor["routes"]:
        assert route["demo_policy"] in {
            "never_unless_explicit",
            "only_when_explicit",
        }
        for step_index, step in enumerate(route["plan"]):
            assert step["kind"] == "skill_run"
            assert step["skill"] == "gentle-cloning"
            template = step.get("input_template")
            if template is None:
                template = _json(SKILL_ROOT / step["input"])
            assert template["schema"] == "gentle.clawbio_skill_request.v1"
            assert template["mode"] in {"shell", "workflow"}
            assert template["claim_attribution_mode"] == "strict"
            assert template["presentation_profile"] == "pcr_primer_design"
            assert template["delegation"] == {
                "schema": "gentle.clawbio_skill_delegation.v1",
                "source_skill": "gentle-pcr-primer-design",
                "source_skill_version": "0.3.0",
                "intent_id": route["intent_id"],
                "plan_step_index": step_index,
            }

    wrapper = _load_generic_wrapper()
    assert {"shell", "workflow", "op"} <= set(wrapper.SUPPORTED_REQUEST_MODES)
    for route in descriptor["routes"]:
        for step in route["plan"]:
            template = step.get("input_template")
            if template is None:
                template = _json(SKILL_ROOT / step["input"])
            parsed = wrapper._coerce_request(template)
            assert parsed.delegation["intent_id"] == route["intent_id"]


def test_experimental_assay_report_requests_gentle_gel_and_attributed_bundle() -> None:
    route = _routes()["experimental_assay_report"]
    step = route["plan"][0]
    request = step["input_template"]

    assert step["skill"] == "gentle-cloning"
    assert request["claim_attribution_mode"] == "strict"
    assert request["presentation_profile"] == "pcr_primer_design"
    assert "--gel-svg {output_prefix}.gel.svg" in request["shell_line"]
    assert request["expected_artifacts"] == [
        "{output_prefix}.json",
        "{output_prefix}.order.tsv",
        "{output_prefix}.gel.svg",
    ]


def test_transcript_panel_separates_rt_priming_from_5prime_capture_claims() -> None:
    route = _routes()["transcript_assay_panel_design"]
    step = route["plan"][0]
    request = step["input_template"]
    slots = step["slots"]

    assert request["shell_line"] == (
        "primers design-transcript-assay-panel "
        "@{operation_path} --backend {backend}"
    )
    assert request["input_claims"] == [
        "RT priming chemistry: {cdna_priming_strategy}.",
        "5-prime capture method: {five_prime_capture_method}.",
        (
            "5-prime completeness evidence or protocol: "
            "{five_prime_capture_evidence}."
        ),
    ]
    assert slots["cdna_priming_strategy"]["default"] == "unspecified"
    assert "oligo_dt" in slots["cdna_priming_strategy"]["choices"]
    assert slots["five_prime_capture_method"]["default"] == "none"
    assert {
        "cap_dependent_5prime_race",
        "rlm_race",
        "cap_trapping",
        "cage",
        "template_switching",
    } <= set(slots["five_prime_capture_method"]["choices"])
    assert slots["five_prime_capture_evidence"]["default"] == "not supplied"
    assert "five_prime_capture_method" not in request["shell_line"]
    assert "five_prime_capture_evidence" not in request["shell_line"]


def test_patz1_example_marks_oligo_dt_and_missing_cap_evidence_as_input() -> None:
    request = _json(
        SKILL_ROOT / "examples" / "request_patz1_transcript_assay_offline.json"
    )
    claims = request["input_claims"]

    assert any("oligo-dT" in claim for claim in claims)
    assert any("capture method: none" in claim for claim in claims)
    assert any("not supplied" in claim for claim in claims)
    assert request["claim_attribution_mode"] == "strict"


def test_route_inputs_and_example_workflows_resolve() -> None:
    descriptor = _descriptor()
    for route in descriptor["routes"]:
        for step in route["plan"]:
            input_path = step.get("input")
            if input_path is not None:
                assert (SKILL_ROOT / input_path).is_file(), (route["intent_id"], input_path)

    operation_source = (REPO_ROOT / "src" / "engine.rs").read_text(encoding="utf-8")
    workflow_paths = []
    for request_path in (SKILL_ROOT / "examples").glob("request_*.json"):
        request = _json(request_path)
        assert request["schema"] == "gentle.clawbio_skill_request.v1"
        workflow_path = request.get("workflow_path")
        if workflow_path:
            workflow_paths.append(REPO_ROOT / workflow_path)
    for workflow_path in workflow_paths:
        assert workflow_path.is_file(), workflow_path
        workflow = _json(workflow_path)
        assert workflow["schema"] == "gentle.workflow_example.v1"
        for operation in workflow["workflow"]["ops"]:
            assert len(operation) == 1
            operation_name = next(iter(operation))
            assert operation_name in operation_source


def test_synthetic_import_fixture_hash_and_provenance_match() -> None:
    fixture = SKILL_ROOT / "examples" / "synthetic_imported_primer_pairs.json"
    workflow = _json(
        SKILL_ROOT
        / "examples"
        / "workflows"
        / "imported_primer_review_offline.workflow.json"
    )
    digest = hashlib.sha256(fixture.read_bytes()).hexdigest()
    request = workflow["workflow"]["ops"][1]["ImportExternalPrimerPairs"]["request"]

    assert request["input_provenance"]["source_sha256"] == f"sha256:{digest}"
    assert request["input_provenance"]["source_path"] == str(
        fixture.relative_to(REPO_ROOT)
    )
    assert request["batch"] == _json(fixture)


def test_plain_primer_language_is_owned_by_specialized_descriptor() -> None:
    specialized = _descriptor()
    generic = _json(GENERIC_ROOT / "INTENTS.json")
    specialized_terms = _trigger_terms(specialized)
    generic_terms = _trigger_terms(generic)
    specialized_catalog_terms = {
        term.lower()
        for term in _json(SKILL_ROOT / "catalog_entry.json")["trigger_keywords"]
    }
    generic_catalog_terms = {
        term.lower()
        for term in _json(GENERIC_ROOT / "catalog_entry.json")["trigger_keywords"]
    }

    assert specialized_terms.isdisjoint(generic_terms)
    assert specialized_catalog_terms.isdisjoint(generic_catalog_terms)
    assert {
        "primer design",
        "endpoint rt-pcr panel",
        "cap-dependent 5-prime race",
        "rlm-race primer design",
        "full-length cdna validation",
        "sybr junction assay",
        "taqman design",
        "dual-space primer specificity",
        "collection primer specificity",
        "import primer pairs",
    } <= specialized_terms

    for term in (
        "simple pcr",
        "primer preflight",
        "transcript qpcr panel",
        "show primer report",
        "pcr protocol cartoon",
    ):
        assert f"gentle-cloning {term}" in generic_terms


def test_demo_routes_cannot_be_selected_as_ordinary_work() -> None:
    routes = _routes()
    demo_route_ids = {
        "offline_conventional_pcr_demo",
        "offline_patz1_transcript_panel_demo",
        "offline_imported_primer_demo",
        "specificity_plan_only_demo",
        "missing_primer3_preflight_demo",
    }
    assert demo_route_ids <= routes.keys()

    for route_id, route in routes.items():
        referenced = [
            step["input"]
            for step in route["plan"]
            if isinstance(step.get("input"), str)
        ]
        if route_id in demo_route_ids:
            assert route["demo_policy"] == "only_when_explicit"
            assert referenced
        else:
            assert route["demo_policy"] == "never_unless_explicit"
            assert all(
                "offline" not in input_path and "demo" not in input_path
                for input_path in referenced
            )


def test_specificity_plans_keep_genomic_and_transcriptome_spaces_separate() -> None:
    routes = _routes()
    route = routes["dual_space_specificity_plan"]
    steps = {step["id"]: step for step in route["plan"]}

    genomic = steps["genomic_dna_specificity_plan"]["input_template"]["shell_line"]
    transcriptome = steps["transcriptome_cdna_specificity_plan"]["input_template"][
        "shell_line"
    ]
    assert "{genomic_reference_id}" in genomic
    assert "{genomic_output_dir}" in genomic
    assert "{transcriptome_reference_id}" in transcriptome
    assert "{transcriptome_output_dir}" in transcriptome
    assert genomic != transcriptome
    assert (
        routes["genomic_specificity_plan"]["plan"][0]["input_template"]["shell_line"]
        == genomic
    )
    assert (
        routes["transcriptome_specificity_plan"]["plan"][0]["input_template"][
            "shell_line"
        ]
        == transcriptome
    )


def test_report_intents_use_matching_list_show_and_export_commands() -> None:
    routes = _routes()
    assert (
        routes["primer_report_list"]["plan"][0]["input_template"]["shell_line"]
        == "primers list-reports"
    )
    assert (
        routes["primer_report_show"]["plan"][0]["input_template"]["shell_line"]
        == "primers show-report {report_id}"
    )
    assert (
        routes["primer_report_export"]["plan"][0]["input_template"]["shell_line"]
        == "primers export-report {report_id} {output_path}"
    )
    assert (
        routes["qpcr_report_list"]["plan"][0]["input_template"]["shell_line"]
        == "primers list-qpcr-reports"
    )
    assert (
        routes["qpcr_report_show"]["plan"][0]["input_template"]["shell_line"]
        == "primers show-qpcr-report {report_id}"
    )
    assert (
        routes["qpcr_report_export"]["plan"][0]["input_template"]["shell_line"]
        == "primers export-qpcr-report {report_id} {output_path}"
    )


def test_committed_patz1_demo_does_not_upgrade_unrun_specificity() -> None:
    report = _json(
        REPO_ROOT
        / "docs"
        / "tutorial"
        / "generated"
        / "artifacts"
        / "patz1_transcript_assay_panels_cli"
        / "artifacts"
        / "patz1_sybr_juc_panel.report.json"
    )
    summaries = [
        assay["primer_pair_summary"]
        for assay in report["selected_assays"]
        if assay.get("primer_pair_summary")
    ]
    followups = report["specificity_followups"]

    assert summaries
    assert followups
    assert {
        summary["whole_genome_specificity_status"] for summary in summaries
    } == {"not_run"}
    assert {
        followup["genomic_confirmation_status"] for followup in followups
    } == {"not_run"}

    skill_text = "\n".join(
        [
            (SKILL_ROOT / "SKILL.md").read_text(encoding="utf-8"),
            (SKILL_ROOT / "README.md").read_text(encoding="utf-8"),
        ]
    )
    normalized_skill_text = " ".join(skill_text.split())
    assert "`not_assessed` is not a pass" in skill_text
    assert "not order-ready" in normalized_skill_text


def test_missing_primer3_example_is_diagnostic_only() -> None:
    request = _json(SKILL_ROOT / "examples" / "request_missing_primer3_preflight.json")
    assert request["mode"] == "shell"
    assert request["shell_line"].startswith("primers preflight ")
    assert "/definitely/missing/gentle-primer3-core" in request["shell_line"]
    assert "install" not in request["shell_line"].lower()


def test_specialized_routes_do_not_submit_orders_or_duplicate_science() -> None:
    descriptor_text = (SKILL_ROOT / "INTENTS.json").read_text(encoding="utf-8")
    forbidden_shell_fragments = (
        "services project-quote",
        "primers oligo-order quote",
        "submit-order",
    )
    assert all(fragment not in descriptor_text for fragment in forbidden_shell_fragments)

    production_python = [
        path
        for path in SKILL_ROOT.rglob("*.py")
        if "tests" not in path.parts
    ]
    assert production_python == []
