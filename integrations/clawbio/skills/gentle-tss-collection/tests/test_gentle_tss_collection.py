from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import re
import sys


SKILL_ROOT = Path(__file__).resolve().parents[1]
GENERIC_ROOT = SKILL_ROOT.parent / "gentle-cloning"
TUTORIAL_SOURCE = (
    SKILL_ROOT.parents[3] / "docs/tutorial/sources/08-15_tss_collection_gui.json"
)


def _json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _routes() -> dict[str, dict]:
    return {
        route["intent_id"]: route
        for route in _json(SKILL_ROOT / "INTENTS.json")["routes"]
    }


def _load_generic_wrapper():
    path = GENERIC_ROOT / "gentle_cloning.py"
    spec = importlib.util.spec_from_file_location("gentle_cloning_tss_delegate", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_metadata_and_delegate_contract_are_consistent() -> None:
    descriptor = _json(SKILL_ROOT / "INTENTS.json")
    catalog = _json(SKILL_ROOT / "catalog_entry.json")
    skill_md = (SKILL_ROOT / "SKILL.md").read_text(encoding="utf-8")

    assert descriptor["schema"] == "clawbio.skill_intents.v1"
    assert descriptor["skill"] == catalog["name"] == "gentle-tss-collection"
    assert catalog["version"] == "0.1.0"
    assert catalog["has_script"] is False
    assert catalog["execution_delegate"] == "gentle-cloning"
    assert catalog["delegate_contract"]["skill_version"] == "0.3.0"
    assert "name: gentle-tss-collection" in skill_md
    assert "version: 0.1.0" in skill_md
    assert not list(SKILL_ROOT.glob("*.py"))


def test_every_route_delegates_to_the_registered_runtime() -> None:
    wrapper = _load_generic_wrapper()
    routes = _routes()
    assert len(routes) == 6

    for route_id, route in routes.items():
        assert route["demo_policy"] == "never_unless_explicit"
        assert len(route["plan"]) == 1
        step = route["plan"][0]
        request = step["input_template"]
        assert step["kind"] == "skill_run"
        assert step["skill"] == "gentle-cloning"
        assert request["mode"] == "shell"
        assert request["claim_attribution_mode"] == "strict"
        assert request["delegation"] == {
            "schema": "gentle.clawbio_skill_delegation.v1",
            "source_skill": "gentle-tss-collection",
            "source_skill_version": "0.1.0",
            "intent_id": route_id,
            "plan_step_index": 0,
        }
        parsed = wrapper._coerce_request(request)
        assert parsed.delegation["intent_id"] == route_id


def test_only_project_mutations_require_a_second_approval() -> None:
    routes = _routes()
    gated = {
        route_id
        for route_id, route in routes.items()
        if route.get("requires_confirmation", False)
    }
    assert gated == {
        "materialize_reviewed_starts",
        "forget_collection_registry_entry",
    }
    for route_id, route in routes.items():
        confirmation = route["plan"][0].get("confirmation", {})
        assert confirmation.get("required", False) is (route_id in gated)


def test_collection_slot_extracts_the_value_not_the_slot_label() -> None:
    routes = _routes()
    for route_id in (
        "validate_stored_collection",
        "open_collection_windows",
        "forget_collection_registry_entry",
    ):
        spec = routes[route_id]["plan"][0]["slots"]["collection_id"]
        match = re.search(
            spec["pattern"],
            "inspect tss collection collection_id=tss_windows",
            flags=re.IGNORECASE,
        )
        assert match is not None
        assert match.group(1) == "tss_windows"


def test_routes_cover_the_tutorial_0815_outer_agent_contract() -> None:
    tutorial_cases = {
        case["id"]: case
        for case in _json(TUTORIAL_SOURCE)["agent_parity"]["cases"]
    }
    routes = _routes()
    assert set(routes) == set(tutorial_cases)

    expected_templates = {
        "preview_annotated_starts": "promoters tss-inventory @{request_path}",
        "materialize_reviewed_starts": "promoters tss-materialize @{request_path}",
        "list_stored_collections": "promoters tss-list",
        "validate_stored_collection": "promoters tss-collection {collection_id}",
        "open_collection_windows": "ui open tss-view --collection {collection_id}",
        "forget_collection_registry_entry": "promoters tss-forget {collection_id}",
    }
    for route_id, template in expected_templates.items():
        assert routes[route_id]["plan"][0]["input_template"]["shell_line"] == template
        assert tutorial_cases[route_id]["execution"] == "ask"
        assert tutorial_cases[route_id]["mutating"] is (
            route_id in {
                "materialize_reviewed_starts",
                "forget_collection_registry_entry",
            }
        )


def test_mutating_delegation_resolves_to_confirmation_gate() -> None:
    wrapper = _load_generic_wrapper()
    route = _routes()["materialize_reviewed_starts"]
    request = json.loads(json.dumps(route["plan"][0]["input_template"]))
    request["state_path"] = "/tmp/gentle-tss-test-state.json"
    request["shell_line"] = (
        "promoters tss-materialize "
        "@docs/examples/assets/tss_tutorial/materialize.request.json"
    )
    request["delegation"]["resolved_slots"] = {
        "request_path": "docs/examples/assets/tss_tutorial/materialize.request.json"
    }

    verified = wrapper._verified_delegation(
        wrapper._coerce_request(request), GENERIC_ROOT / "gentle_cloning.py"
    )

    assert verified is not None
    assert verified["intent_id"] == "materialize_reviewed_starts"
    assert verified["requires_confirmation"] is True
    assert "persists a named collection" in verified["confirmation_reason"]
