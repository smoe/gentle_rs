from __future__ import annotations

import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys


def _skill_script() -> Path:
    return Path(__file__).resolve().parents[1] / "gentle_cloning.py"


def _apptainer_script() -> Path:
    return Path(__file__).resolve().parents[1] / "gentle_apptainer_cli.sh"


def _local_checkout_script() -> Path:
    return Path(__file__).resolve().parents[1] / "gentle_local_checkout_cli.sh"


def test_demo_writes_expected_artifacts(tmp_path: Path) -> None:
    output_dir = tmp_path / "demo_out"
    cmd = [
        sys.executable,
        str(_skill_script()),
        "--demo",
        "--output",
        str(output_dir),
        "--gentle-cli",
        "true",
    ]
    run = subprocess.run(cmd, capture_output=True, text=True, check=False)
    assert run.returncode == 0, run.stderr

    payload = json.loads(run.stdout)
    assert payload["schema"] == "gentle.clawbio_skill_result.v1"
    assert payload["status"] in ("ok", "degraded_demo")
    assert payload["stdout_json"] is None
    assert payload["chat_summary_lines"] is None
    assert (output_dir / "report.md").exists()
    assert (output_dir / "result.json").exists()
    assert (output_dir / "reproducibility" / "commands.sh").exists()
    assert (output_dir / "reproducibility" / "environment.yml").exists()
    assert (output_dir / "reproducibility" / "checksums.sha256").exists()


def test_rejects_invalid_request_schema(tmp_path: Path) -> None:
    bad_req = tmp_path / "bad.json"
    bad_req.write_text(
        json.dumps({"schema": "bad.schema", "mode": "capabilities"}) + "\n",
        encoding="utf-8",
    )
    out_dir = tmp_path / "out"
    cmd = [
        sys.executable,
        str(_skill_script()),
        "--input",
        str(bad_req),
        "--output",
        str(out_dir),
        "--gentle-cli",
        "true",
    ]
    run = subprocess.run(cmd, capture_output=True, text=True, check=False)
    assert run.returncode != 0
    result_json = json.loads((out_dir / "result.json").read_text(encoding="utf-8"))
    assert result_json["status"] == "failed"
    assert "unsupported request schema" in result_json["error"]


def test_result_payload_promotes_sequence_context_chat_summary(tmp_path: Path) -> None:
    request_path = tmp_path / "request.json"
    request_path.write_text(
        json.dumps(
            {
                "schema": "gentle.clawbio_skill_request.v1",
                "mode": "op",
                "operation": {
                    "InspectSequenceContextView": {
                        "seq_id": "rs9923231_vkorc1",
                        "mode": "linear",
                        "viewport_start_0based": 2400,
                        "viewport_end_0based_exclusive": 3501,
                        "coordinate_mode": "genomic",
                    }
                },
            }
        )
        + "\n",
        encoding="utf-8",
    )

    fake_cli = tmp_path / "fake_cli.sh"
    fake_cli.write_text(
        "#!/usr/bin/env bash\n"
        "cat <<'JSON'\n"
        '{"schema":"gentle.sequence_context_view.v1","seq_id":"rs9923231_vkorc1",'
        '"summary_lines":["VKORC1 context around rs9923231","Visible classes: gene, mrna, variation"],'
        '"rows":[]}\n'
        "JSON\n",
        encoding="utf-8",
    )
    fake_cli.chmod(0o755)

    output_dir = tmp_path / "out"
    run = subprocess.run(
        [
            sys.executable,
            str(_skill_script()),
            "--input",
            str(request_path),
            "--output",
            str(output_dir),
            "--gentle-cli",
            str(fake_cli),
        ],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        check=False,
    )
    assert run.returncode == 0, run.stderr

    result = json.loads((output_dir / "result.json").read_text(encoding="utf-8"))
    assert result["stdout_json"]["schema"] == "gentle.sequence_context_view.v1"
    assert result["chat_summary_lines"] == [
        "VKORC1 context around rs9923231",
        "Visible classes: gene, mrna, variation",
    ]
    report = (output_dir / "report.md").read_text(encoding="utf-8")
    assert "## Chat Summary" in report
    assert "- VKORC1 context around rs9923231" in report
    assert "- Visible classes: gene, mrna, variation" in report


def test_result_payload_promotes_bundle_nested_sequence_context_chat_summary(
    tmp_path: Path,
) -> None:
    request_path = tmp_path / "request.json"
    request_path.write_text(
        json.dumps(
            {
                "schema": "gentle.clawbio_skill_request.v1",
                "mode": "op",
                "operation": {
                    "ExportSequenceContextBundle": {
                        "seq_id": "rs9923231_vkorc1",
                        "mode": "linear",
                        "viewport_start_0based": 2400,
                        "viewport_end_0based_exclusive": 3501,
                        "coordinate_mode": "genomic",
                        "include_feature_bed": True,
                        "include_text_summary": True,
                        "include_restriction_sites": False,
                        "restriction_enzymes": [],
                        "output_dir": "artifacts/rs9923231_vkorc1.context_bundle",
                    }
                },
            }
        )
        + "\n",
        encoding="utf-8",
    )

    fake_cli = tmp_path / "fake_cli.sh"
    fake_cli.write_text(
        "#!/usr/bin/env bash\n"
        "cat <<'JSON'\n"
        '{"op_id":"op-demo","messages":["bundle exported"],'
        '"sequence_context_bundle":{"schema":"gentle.sequence_context_bundle.v1",'
        '"seq_id":"rs9923231_vkorc1",'
        '"artifacts":['
        '{"artifact_id":"context_svg","path":"artifacts/context.svg","caption":"Context figure","recommended_use":"best_first_figure","presentation_rank":0,"is_best_first_artifact":true},'
        '{"artifact_id":"context_summary_text","path":"artifacts/context_summary.txt","caption":"Summary text","recommended_use":"compact_chat_summary","presentation_rank":2,"is_best_first_artifact":false}'
        '],'
        '"sequence_context_view":{"schema":"gentle.sequence_context_view.v1",'
        '"summary_lines":["VKORC1 context around rs9923231","Visible classes: gene, mrna, variation"]}}}\n'
        "JSON\n",
        encoding="utf-8",
    )
    fake_cli.chmod(0o755)

    output_dir = tmp_path / "out"
    run = subprocess.run(
        [
            sys.executable,
            str(_skill_script()),
            "--input",
            str(request_path),
            "--output",
            str(output_dir),
            "--gentle-cli",
            str(fake_cli),
        ],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        check=False,
    )
    assert run.returncode == 0, run.stderr

    result = json.loads((output_dir / "result.json").read_text(encoding="utf-8"))
    assert result["stdout_json"]["sequence_context_bundle"]["schema"] == (
        "gentle.sequence_context_bundle.v1"
    )
    assert result["chat_summary_lines"] == [
        "VKORC1 context around rs9923231",
        "Visible classes: gene, mrna, variation",
    ]
    assert result["preferred_artifacts"] == [
        {
            "artifact_id": "context_svg",
            "path": "artifacts/context.svg",
            "caption": "Context figure",
            "recommended_use": "best_first_figure",
            "presentation_rank": 0,
            "is_best_first_artifact": True,
        },
        {
            "artifact_id": "context_summary_text",
            "path": "artifacts/context_summary.txt",
            "caption": "Summary text",
            "recommended_use": "compact_chat_summary",
            "presentation_rank": 2,
            "is_best_first_artifact": False,
        },
    ]
    report = (output_dir / "report.md").read_text(encoding="utf-8")
    assert "## Chat Summary" in report
    assert "## Preferred Artifacts" in report
    assert "artifacts/context.svg" in report
    assert "- VKORC1 context around rs9923231" in report
    assert "- Visible classes: gene, mrna, variation" in report


def test_wrapper_builds_variant_storyboard_from_collected_svgs(tmp_path: Path) -> None:
    request_path = tmp_path / "request.json"
    request_path.write_text(
        json.dumps(
            {
                "schema": "gentle.clawbio_skill_request.v1",
                "mode": "workflow",
                "state_path": ".gentle_state.json",
                "workflow_path": "variant_luciferase_storyboard_demo.json",
                "timeout_secs": 1800,
                "expected_artifacts": [
                    "docs/tutorial/reproducibility/vkorc1_rs9923231_promoter_reporter/vkorc1_rs9923231_promoter_context.svg",
                    "docs/tutorial/reproducibility/vkorc1_rs9923231_promoter_reporter/vkorc1_rs9923231_reporter_reference.svg",
                    "docs/tutorial/reproducibility/vkorc1_rs9923231_promoter_reporter/vkorc1_rs9923231_reporter_alternate.svg",
                ],
            }
        )
        + "\n",
        encoding="utf-8",
    )

    figure_dir = (
        tmp_path
        / "docs"
        / "tutorial"
        / "reproducibility"
        / "vkorc1_rs9923231_promoter_reporter"
    )
    figure_dir.mkdir(parents=True, exist_ok=True)
    for name, fill in (
        ("vkorc1_rs9923231_promoter_context.svg", "#cbd5e1"),
        ("vkorc1_rs9923231_reporter_reference.svg", "#86efac"),
        ("vkorc1_rs9923231_reporter_alternate.svg", "#fca5a5"),
    ):
        (figure_dir / name).write_text(
            (
                '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 400 240">'
                f'<rect width="400" height="240" fill="{fill}" />'
                f'<text x="24" y="44" font-size="28">{name}</text>'
                "</svg>\n"
            ),
            encoding="utf-8",
        )

    fake_cli = tmp_path / "fake_cli.sh"
    fake_cli.write_text(
        "#!/usr/bin/env bash\n"
        "cat <<'JSON'\n"
        '{"schema":"gentle.workflow_run.v1","run_id":"example_vkorc1"}\n'
        "JSON\n",
        encoding="utf-8",
    )
    fake_cli.chmod(0o755)

    output_dir = tmp_path / "out"
    run = subprocess.run(
        [
            sys.executable,
            str(_skill_script()),
            "--input",
            str(request_path),
            "--output",
            str(output_dir),
            "--gentle-cli",
            str(fake_cli),
        ],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        check=False,
    )
    assert run.returncode == 0, run.stderr

    result = json.loads((output_dir / "result.json").read_text(encoding="utf-8"))
    preferred = result["preferred_artifacts"]
    assert preferred[0]["artifact_id"] == "clawbio_storyboard_svg"
    assert preferred[0]["path"] == "generated/clawbio_storyboard.svg"
    assert preferred[0]["is_best_first_artifact"] is True
    assert any(
        artifact["path"].endswith("vkorc1_rs9923231_promoter_context.svg")
        for artifact in preferred[1:]
    )
    storyboard_path = output_dir / "generated" / "clawbio_storyboard.svg"
    assert storyboard_path.exists()
    storyboard = storyboard_path.read_text(encoding="utf-8")
    assert "Variant-to-Synthetic-Biology assay storyboard" in storyboard
    assert "Genomic context" in storyboard
    assert "Reference allele reporter" in storyboard
    assert "Alternate allele reporter" in storyboard
    report = (output_dir / "report.md").read_text(encoding="utf-8")
    assert "generated/clawbio_storyboard.svg" in report


def test_apptainer_launcher_wraps_gentle_cli_with_bind_mount(tmp_path: Path) -> None:
    fake_runtime = tmp_path / "apptainer"
    capture_path = tmp_path / "apptainer_args.txt"
    fake_runtime.write_text(
        "#!/usr/bin/env bash\n"
        "printf '%s\\n' \"$@\" > \"$FAKE_APPTAINER_ARGS\"\n",
        encoding="utf-8",
    )
    fake_runtime.chmod(0o755)

    image_path = tmp_path / "gentle.sif"
    image_path.write_text("", encoding="utf-8")

    env = dict(os.environ)
    env["PATH"] = f"{tmp_path}:{env.get('PATH', '')}"
    env["FAKE_APPTAINER_ARGS"] = str(capture_path)

    run = subprocess.run(
        ["bash", str(_apptainer_script()), str(image_path), "capabilities"],
        cwd=tmp_path,
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )
    assert run.returncode == 0, run.stderr

    args = capture_path.read_text(encoding="utf-8").splitlines()
    expected_host_workdir = str(tmp_path.resolve())
    assert args == [
        "exec",
        "--bind",
        f"{expected_host_workdir}:/work",
        "--pwd",
        "/work",
        str(image_path),
        "gentle_cli",
        "capabilities",
    ]


def test_local_checkout_launcher_uses_repo_root_defaults(tmp_path: Path) -> None:
    fake_repo = tmp_path / "GENtle"
    (fake_repo / "src" / "bin").mkdir(parents=True)
    (fake_repo / "Cargo.toml").write_text("[package]\nname = 'gentle'\n", encoding="utf-8")
    (fake_repo / "src" / "bin" / "gentle_cli.rs").write_text(
        "fn main() {}\n",
        encoding="utf-8",
    )

    fake_cargo = tmp_path / "cargo"
    capture_args = tmp_path / "cargo_args.txt"
    capture_env = tmp_path / "cargo_env.txt"
    fake_cargo.write_text(
        "#!/usr/bin/env bash\n"
        "printf '%s\\n' \"$@\" > \"$FAKE_CARGO_ARGS\"\n"
        "{\n"
        "  printf 'CARGO_TARGET_DIR=%s\\n' \"$CARGO_TARGET_DIR\"\n"
        "  printf 'GENTLE_REFERENCE_CACHE_DIR=%s\\n' \"$GENTLE_REFERENCE_CACHE_DIR\"\n"
        "  printf 'GENTLE_HELPER_CACHE_DIR=%s\\n' \"$GENTLE_HELPER_CACHE_DIR\"\n"
        "} > \"$FAKE_CARGO_ENV\"\n",
        encoding="utf-8",
    )
    fake_cargo.chmod(0o755)

    env = dict(os.environ)
    env["PATH"] = f"{tmp_path}:{env.get('PATH', '')}"
    env["GENTLE_REPO_ROOT"] = str(fake_repo)
    env["FAKE_CARGO_ARGS"] = str(capture_args)
    env["FAKE_CARGO_ENV"] = str(capture_env)
    env.pop("CARGO_TARGET_DIR", None)
    env.pop("GENTLE_REFERENCE_CACHE_DIR", None)
    env.pop("GENTLE_HELPER_CACHE_DIR", None)

    run = subprocess.run(
        ["bash", str(_local_checkout_script()), "capabilities"],
        cwd=tmp_path,
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )
    assert run.returncode == 0, run.stderr

    args = capture_args.read_text(encoding="utf-8").splitlines()
    assert args == [
        "run",
        "--quiet",
        "--manifest-path",
        str(fake_repo / "Cargo.toml"),
        "--bin",
        "gentle_cli",
        "--",
        "capabilities",
    ]
    env_lines = capture_env.read_text(encoding="utf-8").splitlines()
    assert env_lines == [
        f"CARGO_TARGET_DIR={fake_repo / 'target'}",
        f"GENTLE_REFERENCE_CACHE_DIR={fake_repo / 'data' / 'genomes'}",
        f"GENTLE_HELPER_CACHE_DIR={fake_repo / 'data' / 'helper_genomes'}",
    ]


def test_workflow_request_resolves_against_gentle_repo_root(tmp_path: Path) -> None:
    fake_repo = tmp_path / "GENtle"
    (fake_repo / "src" / "bin").mkdir(parents=True)
    (fake_repo / "Cargo.toml").write_text("[package]\nname = 'gentle'\n", encoding="utf-8")
    (fake_repo / "src" / "bin" / "gentle_cli.rs").write_text(
        "fn main() {}\n",
        encoding="utf-8",
    )
    workflow_rel = Path("docs/examples/workflows/demo_workflow.json")
    workflow_abs = fake_repo / workflow_rel
    workflow_abs.parent.mkdir(parents=True)
    workflow_abs.write_text('{"schema":"gentle.workflow.v1"}\n', encoding="utf-8")

    request_path = tmp_path / "request.json"
    request_path.write_text(
        json.dumps(
            {
                "schema": "gentle.clawbio_skill_request.v1",
                "mode": "workflow",
                "workflow_path": str(workflow_rel),
                "timeout_secs": 180,
            }
        )
        + "\n",
        encoding="utf-8",
    )

    fake_cli = tmp_path / "capture_cli.sh"
    capture_path = tmp_path / "captured_args.txt"
    capture_cwd = tmp_path / "captured_cwd.txt"
    fake_cli.write_text(
        "#!/usr/bin/env bash\n"
        "printf '%s\\n' \"$@\" > \"$FAKE_CAPTURE_ARGS\"\n"
        "pwd > \"$FAKE_CAPTURE_CWD\"\n",
        encoding="utf-8",
    )
    fake_cli.chmod(0o755)

    run_cwd = tmp_path / "ClawBio" / "skills" / "gentle-cloning"
    run_cwd.mkdir(parents=True)
    output_dir = tmp_path / "out"

    env = dict(os.environ)
    env["GENTLE_REPO_ROOT"] = str(fake_repo)
    env["FAKE_CAPTURE_ARGS"] = str(capture_path)
    env["FAKE_CAPTURE_CWD"] = str(capture_cwd)

    run = subprocess.run(
        [
            sys.executable,
            str(_skill_script()),
            "--input",
            str(request_path),
            "--output",
            str(output_dir),
            "--gentle-cli",
            str(fake_cli),
        ],
        cwd=run_cwd,
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )
    assert run.returncode == 0, run.stderr
    assert capture_path.read_text(encoding="utf-8").splitlines() == [
        "workflow",
        f"@{workflow_abs.resolve()}",
    ]
    assert capture_cwd.read_text(encoding="utf-8").strip() == str(fake_repo.resolve())


def test_expected_artifacts_are_copied_into_output_bundle(tmp_path: Path) -> None:
    request_path = tmp_path / "request.json"
    request_path.write_text(
        json.dumps(
            {
                "schema": "gentle.clawbio_skill_request.v1",
                "mode": "shell",
                "shell_line": "protocol-cartoon render-svg gibson.two_fragment artifacts/demo.protocol.svg",
                "expected_artifacts": ["artifacts/demo.protocol.svg"],
                "timeout_secs": 180,
            }
        )
        + "\n",
        encoding="utf-8",
    )

    fake_cli = tmp_path / "fake_cli.sh"
    fake_cli.write_text(
        "#!/usr/bin/env bash\n"
        "set -euo pipefail\n"
        "set -- $2\n"
        "output_path=$4\n"
        "mkdir -p \"$(dirname \"$output_path\")\"\n"
        "printf '<svg>demo</svg>\\n' > \"$output_path\"\n",
        encoding="utf-8",
    )
    fake_cli.chmod(0o755)

    output_dir = tmp_path / "out"
    run = subprocess.run(
        [
            sys.executable,
            str(_skill_script()),
            "--input",
            str(request_path),
            "--output",
            str(output_dir),
            "--gentle-cli",
            str(fake_cli),
        ],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        check=False,
    )
    assert run.returncode == 0, run.stderr

    result = json.loads((output_dir / "result.json").read_text(encoding="utf-8"))
    collected = result["artifacts"]["collected"]
    assert len(collected) == 1
    assert collected[0]["declared_path"] == "artifacts/demo.protocol.svg"
    copied_path = Path(collected[0]["copied_path"])
    assert copied_path.exists()
    assert copied_path.read_text(encoding="utf-8") == "<svg>demo</svg>\n"


def test_expected_artifacts_are_sandboxed_under_generated_dir(tmp_path: Path) -> None:
    request_path = tmp_path / "request.json"
    request_path.write_text(
        json.dumps(
            {
                "schema": "gentle.clawbio_skill_request.v1",
                "mode": "shell",
                "shell_line": "protocol-cartoon render-svg gibson.two_fragment ../outside/demo.protocol.svg",
                "expected_artifacts": ["../outside/demo.protocol.svg"],
                "timeout_secs": 180,
            }
        )
        + "\n",
        encoding="utf-8",
    )

    fake_cli = tmp_path / "fake_cli.sh"
    fake_cli.write_text(
        "#!/usr/bin/env bash\n"
        "set -euo pipefail\n"
        "set -- $2\n"
        "output_path=$4\n"
        "mkdir -p \"$(dirname \"$output_path\")\"\n"
        "printf '<svg>demo</svg>\\n' > \"$output_path\"\n",
        encoding="utf-8",
    )
    fake_cli.chmod(0o755)

    run_cwd = tmp_path / "skill"
    run_cwd.mkdir(parents=True)
    output_dir = tmp_path / "out"
    run = subprocess.run(
        [
            sys.executable,
            str(_skill_script()),
            "--input",
            str(request_path),
            "--output",
            str(output_dir),
            "--gentle-cli",
            str(fake_cli),
        ],
        cwd=run_cwd,
        capture_output=True,
        text=True,
        check=False,
    )
    assert run.returncode == 0, run.stderr

    result = json.loads((output_dir / "result.json").read_text(encoding="utf-8"))
    copied_path = Path(result["artifacts"]["collected"][0]["copied_path"])
    generated_root = (output_dir / "generated").resolve()
    assert copied_path.exists()
    assert copied_path.is_relative_to(generated_root)
    assert copied_path == generated_root / "outside" / "demo.protocol.svg"


def test_reference_preflight_runs_status_prepare_and_main_command(tmp_path: Path) -> None:
    request_path = tmp_path / "request.json"
    request_path.write_text(
        json.dumps(
            {
                "schema": "gentle.clawbio_skill_request.v1",
                "mode": "capabilities",
                "timeout_secs": 180,
                "ensure_reference_prepared": {
                    "genome_id": "Human GRCh38 Ensembl 116",
                    "prepare_timeout_secs": 7200,
                },
            }
        )
        + "\n",
        encoding="utf-8",
    )

    fake_cli = tmp_path / "fake_cli.sh"
    fake_cli.write_text(
        "#!/usr/bin/env bash\n"
        "set -euo pipefail\n"
        "sentinel=\"${FAKE_PREP_SENTINEL:?}\"\n"
        "capture=\"${FAKE_CLI_CAPTURE:?}\"\n"
        "printf '%s\\n' \"$*\" >> \"$capture\"\n"
        "if [ \"$1\" = \"genomes\" ] && [ \"$2\" = \"status\" ]; then\n"
        "  if [ -f \"$sentinel\" ]; then\n"
        "    printf '{\"prepared\":true,\"genome_id\":\"%s\",\"sequence_source_type\":\"ensembl_fasta\",\"annotation_source_type\":\"ensembl_gtf\"}\\n' \"$3\"\n"
        "  else\n"
        "    printf '{\"prepared\":false,\"genome_id\":\"%s\",\"sequence_source_type\":\"ensembl_fasta\",\"annotation_source_type\":\"ensembl_gtf\"}\\n' \"$3\"\n"
        "  fi\n"
        "elif [ \"$1\" = \"genomes\" ] && [ \"$2\" = \"prepare\" ]; then\n"
        "  : > \"$sentinel\"\n"
        "  printf '{\"prepared\":true,\"genome_id\":\"%s\"}\\n' \"$3\"\n"
        "else\n"
        "  printf '{\"capabilities\":[\"demo\"]}\\n'\n"
        "fi\n",
        encoding="utf-8",
    )
    fake_cli.chmod(0o755)

    capture_path = tmp_path / "captured.txt"
    sentinel_path = tmp_path / "prepared.flag"
    output_dir = tmp_path / "out"

    env = dict(os.environ)
    env["FAKE_CLI_CAPTURE"] = str(capture_path)
    env["FAKE_PREP_SENTINEL"] = str(sentinel_path)

    run = subprocess.run(
        [
            sys.executable,
            str(_skill_script()),
            "--input",
            str(request_path),
            "--output",
            str(output_dir),
            "--gentle-cli",
            str(fake_cli),
        ],
        cwd=tmp_path,
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )
    assert run.returncode == 0, run.stderr

    result = json.loads((output_dir / "result.json").read_text(encoding="utf-8"))
    preflight = result["preflight"]["reference_preparation"]
    assert preflight["prepared_before"] is False
    assert preflight["prepare_attempted"] is True
    assert preflight["prepared_after"] is True
    assert preflight["status"] == "prepared_during_run"
    assert len(preflight["steps"]) == 3
    assert preflight["steps"][0]["command"] == [
        str(fake_cli),
        "genomes",
        "status",
        "Human GRCh38 Ensembl 116",
    ]
    assert preflight["steps"][1]["command"] == [
        str(fake_cli),
        "genomes",
        "prepare",
        "Human GRCh38 Ensembl 116",
        "--timeout-secs",
        "7200",
    ]
    assert preflight["steps"][2]["command"] == [
        str(fake_cli),
        "genomes",
        "status",
        "Human GRCh38 Ensembl 116",
    ]
    commands_text = (output_dir / "reproducibility" / "commands.sh").read_text(
        encoding="utf-8"
    )
    assert 'genomes status \'Human GRCh38 Ensembl 116\'' in commands_text
    assert (
        'genomes prepare \'Human GRCh38 Ensembl 116\' --timeout-secs 7200'
        in commands_text
    )
    assert "\n" + shlex.quote(str(fake_cli)) + " capabilities\n" in commands_text


def test_failed_command_reports_command_exit_code_and_stderr_preview(tmp_path: Path) -> None:
    request_path = tmp_path / "request.json"
    request_path.write_text(
        json.dumps(
            {
                "schema": "gentle.clawbio_skill_request.v1",
                "mode": "shell",
                "shell_line": 'genomes status "Human GRCh38 Ensembl 116"',
                "timeout_secs": 180,
            }
        )
        + "\n",
        encoding="utf-8",
    )

    fake_cli = tmp_path / "fake_cli.sh"
    fake_cli.write_text(
        "#!/usr/bin/env bash\n"
        "printf 'catalog file missing\\n' >&2\n"
        "exit 17\n",
        encoding="utf-8",
    )
    fake_cli.chmod(0o755)

    output_dir = tmp_path / "out"
    run = subprocess.run(
        [
            sys.executable,
            str(_skill_script()),
            "--input",
            str(request_path),
            "--output",
            str(output_dir),
            "--gentle-cli",
            str(fake_cli),
        ],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        check=False,
    )
    assert run.returncode != 0

    result = json.loads((output_dir / "result.json").read_text(encoding="utf-8"))
    assert result["status"] == "command_failed"
    assert result["failure_summary"] == {
        "stage": "main_command",
        "note": None,
        "command": [
            str(fake_cli),
            "shell",
            'genomes status "Human GRCh38 Ensembl 116"',
        ],
        "command_text": (
            f"{shlex.quote(str(fake_cli))} shell "
            '\'genomes status "Human GRCh38 Ensembl 116"\''
        ),
        "execution_cwd": str(tmp_path.resolve()),
        "exit_code": 17,
        "stderr_preview": "catalog file missing",
        "stdout_preview": None,
    }
    assert "gentle_cli exited with a non-zero status." in result["error"]
    assert "Exit code: `17`." in result["error"]
    assert "catalog file missing" in result["error"]
    report = (output_dir / "report.md").read_text(encoding="utf-8")
    assert "## Failure Summary" in report
    assert "catalog file missing" in report
    assert "Failure stage: `main_command`" in report


def test_unknown_services_shell_command_reports_version_mismatch_hint(tmp_path: Path) -> None:
    request_path = tmp_path / "request.json"
    request_path.write_text(
        json.dumps(
            {
                "schema": "gentle.clawbio_skill_request.v1",
                "mode": "shell",
                "shell_line": "services status",
                "timeout_secs": 180,
            }
        )
        + "\n",
        encoding="utf-8",
    )

    fake_cli = tmp_path / "fake_cli.sh"
    fake_cli.write_text(
        "#!/usr/bin/env bash\n"
        "printf \"Unknown shell command 'services'. Try: help\\n\" >&2\n"
        "exit 1\n",
        encoding="utf-8",
    )
    fake_cli.chmod(0o755)

    output_dir = tmp_path / "out"
    run = subprocess.run(
        [
            sys.executable,
            str(_skill_script()),
            "--input",
            str(request_path),
            "--output",
            str(output_dir),
            "--gentle-cli",
            str(fake_cli),
        ],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        check=False,
    )
    assert run.returncode != 0

    result = json.loads((output_dir / "result.json").read_text(encoding="utf-8"))
    assert result["status"] == "command_failed"
    note = result["failure_summary"]["note"]
    assert "older than the current ClawBio skill scaffold" in note
    assert "gentle_local_checkout_cli.sh" in note
    assert "may not support the requested status route yet" in result["error"]


def test_relative_input_path_resolves_from_copied_clawbio_skill_layout(tmp_path: Path) -> None:
    clawbio_root = tmp_path / "ClawBio"
    skill_dir = clawbio_root / "skills" / "gentle-cloning"
    examples_dir = skill_dir / "examples"
    examples_dir.mkdir(parents=True)

    copied_skill = skill_dir / "gentle_cloning.py"
    shutil.copy2(_skill_script(), copied_skill)

    request_rel = Path("skills/gentle-cloning/examples/request_services_status.json")
    request_abs = clawbio_root / request_rel
    request_abs.write_text(
        json.dumps(
            {
                "schema": "gentle.clawbio_skill_request.v1",
                "mode": "capabilities",
                "timeout_secs": 180,
            }
        )
        + "\n",
        encoding="utf-8",
    )

    output_dir = tmp_path / "out"
    unrelated_cwd = tmp_path / "somewhere-else"
    unrelated_cwd.mkdir(parents=True)

    run = subprocess.run(
        [
            sys.executable,
            str(copied_skill),
            "--input",
            str(request_rel),
            "--output",
            str(output_dir),
            "--gentle-cli",
            "true",
        ],
        cwd=unrelated_cwd,
        capture_output=True,
        text=True,
        check=False,
    )
    assert run.returncode == 0, run.stderr

    result = json.loads((output_dir / "result.json").read_text(encoding="utf-8"))
    assert result["status"] == "ok"
    assert result["request"]["mode"] == "capabilities"


def test_example_requests_cover_bootstrap_analysis_and_typical_request_routes() -> None:
    examples_dir = Path(__file__).resolve().parents[1] / "examples"
    expected = {
        "request_genomes_list_human.json": (
            "shell",
            "genomes list --filter human",
            180,
        ),
        "request_helpers_list_gst.json": (
            "shell",
            "helpers list --filter gst",
            180,
        ),
        "request_hosts_list_deor.json": (
            "shell",
            "hosts list --filter deoR",
            180,
        ),
        "request_genomes_ensembl_available_human.json": (
            "shell",
            "genomes ensembl-available --collection vertebrates --filter human",
            300,
        ),
        "request_genomes_install_ensembl_mouse.json": (
            "shell",
            "genomes install-ensembl mus_musculus --collection vertebrates",
            1800,
        ),
        "request_shell_state_summary.json": (
            "shell",
            "state-summary",
            180,
        ),
        "request_genomes_status_grch38.json": (
            "shell",
            'genomes status "Human GRCh38 Ensembl 116"',
            180,
        ),
        "request_genomes_prepare_grch38.json": (
            "shell",
            'genomes prepare "Human GRCh38 Ensembl 116" --timeout-secs 7200',
            7500,
        ),
        "request_resources_status.json": (
            "shell",
            "resources status",
            180,
        ),
        "request_services_status.json": (
            "shell",
            "services status",
            180,
        ),
        "request_helpers_status_puc19.json": (
            "shell",
            'helpers status "Plasmid pUC19 (online)"',
            180,
        ),
        "request_helpers_prepare_puc19.json": (
            "shell",
            'helpers prepare "Plasmid pUC19 (online)" --timeout-secs 1800',
            2100,
        ),
        "request_genbank_fetch_pbr322.json": (
            "shell",
            "genbank fetch J01749 --as-id pbr322_refseq",
            300,
        ),
        "request_dbsnp_fetch_rs9923231.json": (
            "shell",
            'dbsnp fetch rs9923231 "Human GRCh38 Ensembl 116" --flank-bp 3000 --output-id rs9923231_vkorc1 --annotation-scope core',
            900,
        ),
        "request_render_svg_rs9923231_vkorc1_linear.json": (
            "shell",
            "render-svg rs9923231_vkorc1 linear artifacts/rs9923231_vkorc1.context.linear.svg",
            180,
        ),
        "request_workflow_vkorc1_context_svg_auto_prepare.json": (
            "workflow",
            None,
            1200,
        ),
        "request_inspect_sequence_context_rs9923231_vkorc1.json": (
            "op",
            None,
            180,
        ),
        "request_export_sequence_context_bundle_rs9923231_vkorc1.json": (
            "op",
            None,
            300,
        ),
        "request_export_bed_rs9923231_vkorc1_context_features.json": (
            "shell",
            "features export-bed rs9923231_vkorc1 artifacts/rs9923231_vkorc1.context.features.bed --coordinate-mode genomic --kind gene --kind mRNA --kind variation --sort start --include-source --include-qualifiers",
            180,
        ),
        "request_genomes_extract_gene_tp53.json": (
            "shell",
            'genomes extract-gene "Human GRCh38 Ensembl 116" TP53 --occurrence 1 --output-id grch38_tp53',
            600,
        ),
        "request_genomes_extract_gene_tp53_auto_prepare.json": (
            "shell",
            'genomes extract-gene "Human GRCh38 Ensembl 116" TP53 --occurrence 1 --output-id grch38_tp53',
            1200,
        ),
        "request_genomes_blast_grch38_short.json": (
            "shell",
            'genomes blast "Human GRCh38 Ensembl 116" ACGTACGTACGT --task blastn-short --max-hits 10',
            300,
        ),
        "request_helpers_blast_puc19_short.json": (
            "shell",
            'helpers blast "Plasmid pUC19 (online)" ACGTACGTACGT --task blastn-short --max-hits 10',
            300,
        ),
        "request_find_restriction_sites_inline_sequence_ecori_smai.json": (
            "op",
            None,
            180,
        ),
        "request_scan_tfbs_hits_inline_sequence_sp1_tp73.json": (
            "op",
            None,
            180,
        ),
        "request_workflow_tp73_tfbs_score_tracks_summary.json": (
            "workflow",
            None,
            300,
        ),
        "request_workflow_tp73_tfbs_score_tracks_svg.json": (
            "workflow",
            None,
            300,
        ),
        "request_resources_summarize_jaspar_sp1_rest.json": (
            "shell",
            "resources summarize-jaspar --motif SP1 --motif REST --random-length 10000 --seed 123 --output artifacts/jaspar_sp1_rest.presentation.json",
            180,
        ),
        "request_render_svg_pgex_fasta_circular.json": (
            "shell",
            "render-svg pgex_fasta circular artifacts/pgex_fasta.circular.svg",
            180,
        ),
        "request_tfbs_summary_pgex_fasta.json": (
            "shell",
            "features tfbs-summary pgex_fasta --focus 1..1500 --context 1..4904 --limit 10",
            180,
        ),
        "request_inspect_feature_expert_pgex_fasta_tfbs.json": (
            "shell",
            "inspect-feature-expert pgex_fasta tfbs 0",
            180,
        ),
        "request_render_feature_expert_pgex_fasta_tfbs_svg.json": (
            "shell",
            "render-feature-expert-svg pgex_fasta tfbs 0 artifacts/pgex_fasta.tfbs.expert.svg",
            180,
        ),
        "request_inspect_feature_expert_pgex_fasta_restriction_ecori.json": (
            "shell",
            "inspect-feature-expert pgex_fasta restriction 944 --enzyme EcoRI --start 944 --end 949",
            180,
        ),
        "request_render_feature_expert_pgex_fasta_restriction_ecori_svg.json": (
            "shell",
            "render-feature-expert-svg pgex_fasta restriction 944 --enzyme EcoRI --start 944 --end 949 artifacts/pgex_fasta.restriction.ecori.expert.svg",
            180,
        ),
        "request_inspect_feature_expert_tp53_isoform.json": (
            "shell",
            "inspect-feature-expert grch38_tp53 isoform tp53_isoforms_v1",
            180,
        ),
        "request_inspect_feature_expert_tp53_splicing.json": (
            "shell",
            "inspect-feature-expert tp53_panel_source splicing 2",
            180,
        ),
        "request_export_bed_grch38_tp53_gene_models.json": (
            "shell",
            "features export-bed grch38_tp53 artifacts/grch38_tp53.gene_models.bed --coordinate-mode auto --kind gene --kind mRNA --kind exon --kind CDS --sort start --include-source --include-qualifiers",
            180,
        ),
        "request_protocol_cartoon_gibson_svg.json": (
            "shell",
            "protocol-cartoon render-svg gibson.two_fragment artifacts/gibson.two_fragment.protocol.svg",
            180,
        ),
        "request_protocol_cartoon_qpcr_svg.json": (
            "shell",
            "protocol-cartoon render-svg pcr.assay.qpcr artifacts/qpcr.assay.protocol.svg",
            180,
        ),
    }

    for name, (mode, shell_line, timeout_secs) in expected.items():
        payload = json.loads((examples_dir / name).read_text(encoding="utf-8"))
        assert payload["schema"] == "gentle.clawbio_skill_request.v1"
        assert payload["mode"] == mode
        if shell_line is not None:
            assert payload["shell_line"] == shell_line
        assert payload["timeout_secs"] == timeout_secs
        if name in {
            "request_shell_state_summary.json",
            "request_genbank_fetch_pbr322.json",
            "request_dbsnp_fetch_rs9923231.json",
            "request_inspect_sequence_context_rs9923231_vkorc1.json",
            "request_export_sequence_context_bundle_rs9923231_vkorc1.json",
            "request_render_svg_rs9923231_vkorc1_linear.json",
            "request_export_bed_rs9923231_vkorc1_context_features.json",
            "request_genomes_extract_gene_tp53.json",
            "request_genomes_extract_gene_tp53_auto_prepare.json",
            "request_tfbs_summary_pgex_fasta.json",
            "request_inspect_feature_expert_pgex_fasta_tfbs.json",
            "request_render_feature_expert_pgex_fasta_tfbs_svg.json",
            "request_inspect_feature_expert_pgex_fasta_restriction_ecori.json",
            "request_render_feature_expert_pgex_fasta_restriction_ecori_svg.json",
            "request_inspect_feature_expert_tp53_isoform.json",
            "request_inspect_feature_expert_tp53_splicing.json",
            "request_export_bed_grch38_tp53_gene_models.json",
            "request_render_svg_pgex_fasta_circular.json",
        }:
            assert payload["state_path"] == ".gentle_state.json"
        if name == "request_protocol_cartoon_gibson_svg.json":
            assert payload["expected_artifacts"] == [
                "artifacts/gibson.two_fragment.protocol.svg"
            ]
        if name == "request_render_svg_pgex_fasta_circular.json":
            assert payload["expected_artifacts"] == [
                "artifacts/pgex_fasta.circular.svg"
            ]
        if name == "request_render_svg_rs9923231_vkorc1_linear.json":
            assert payload["expected_artifacts"] == [
                "artifacts/rs9923231_vkorc1.context.linear.svg"
            ]
        if name == "request_workflow_vkorc1_context_svg_auto_prepare.json":
            assert payload["state_path"] == ".gentle_state.json"
            assert payload["ensure_reference_prepared"] == {
                "genome_id": "Human GRCh38 Ensembl 116",
                "catalog_path": "assets/genomes.json",
                "cache_dir": "data/genomes",
                "prepare_timeout_secs": 7200,
            }
            assert payload["expected_artifacts"] == [
                "artifacts/rs9923231_vkorc1.context.demo.svg"
            ]
            workflow = payload["workflow"]
            assert workflow["run_id"] == "clawbio_vkorc1_context_svg_auto_prepare"
            ops = workflow["ops"]
            assert ops[0]["FetchDbSnpRegion"]["rs_id"] == "rs9923231"
            assert ops[-1]["RenderSequenceSvg"]["path"] == (
                "artifacts/rs9923231_vkorc1.context.demo.svg"
            )
        if name == "request_genomes_extract_gene_tp53_auto_prepare.json":
            assert payload["state_path"] == ".gentle_state.json"
            assert payload["ensure_reference_prepared"] == {
                "genome_id": "Human GRCh38 Ensembl 116",
                "catalog_path": "assets/genomes.json",
                "cache_dir": "data/genomes",
                "prepare_timeout_secs": 7200,
            }
        if name == "request_workflow_vkorc1_planning.json":
            assert payload["state_path"] == ".gentle_state.json"
            assert payload["expected_artifacts"] == [
                "docs/tutorial/reproducibility/vkorc1_rs9923231_promoter_reporter/vkorc1_rs9923231_promoter_context.svg",
                "docs/tutorial/reproducibility/vkorc1_rs9923231_promoter_reporter/vkorc1_rs9923231_reporter_reference.svg",
                "docs/tutorial/reproducibility/vkorc1_rs9923231_promoter_reporter/vkorc1_rs9923231_reporter_alternate.svg",
            ]
        if name == "request_inspect_sequence_context_rs9923231_vkorc1.json":
            assert payload["operation"] == {
                "InspectSequenceContextView": {
                    "seq_id": "rs9923231_vkorc1",
                    "mode": "linear",
                    "viewport_start_0based": 2400,
                    "viewport_end_0based_exclusive": 3501,
                    "include_visible_classes": [
                        "gene",
                        "mrna",
                        "variation",
                        "tfbs",
                    ],
                    "coordinate_mode": "genomic",
                    "limit": 20,
                }
            }
        if name == "request_export_sequence_context_bundle_rs9923231_vkorc1.json":
            assert payload["operation"] == {
                "ExportSequenceContextBundle": {
                    "seq_id": "rs9923231_vkorc1",
                    "mode": "linear",
                    "viewport_start_0based": 2400,
                    "viewport_end_0based_exclusive": 3501,
                    "coordinate_mode": "genomic",
                    "include_feature_bed": True,
                    "include_text_summary": True,
                    "include_restriction_sites": False,
                    "restriction_enzymes": [],
                    "output_dir": "artifacts/rs9923231_vkorc1.context_bundle",
                }
            }
            assert payload["expected_artifacts"] == [
                "artifacts/rs9923231_vkorc1.context_bundle/bundle.json",
                "artifacts/rs9923231_vkorc1.context_bundle/context.svg",
                "artifacts/rs9923231_vkorc1.context_bundle/context_summary.json",
                "artifacts/rs9923231_vkorc1.context_bundle/context_summary.txt",
                "artifacts/rs9923231_vkorc1.context_bundle/context_features.bed",
            ]
        if name == "request_find_restriction_sites_inline_sequence_ecori_smai.json":
            assert payload["operation"] == {
                "FindRestrictionSites": {
                    "target": {
                        "kind": "inline_sequence",
                        "sequence_text": "GAATTCCCGGGATCC",
                        "topology": "linear",
                        "id_hint": "inline_ecori_smai_window",
                        "span_start_0based": 0,
                        "span_end_0based_exclusive": 15,
                    },
                    "enzymes": ["EcoRI", "SmaI"],
                    "include_cut_geometry": True,
                    "path": "artifacts/inline_ecori_smai.restriction_scan.json",
                }
            }
            assert payload["expected_artifacts"] == [
                "artifacts/inline_ecori_smai.restriction_scan.json"
            ]
        if name == "request_scan_tfbs_hits_inline_sequence_sp1_tp73.json":
            assert payload["operation"] == {
                "ScanTfbsHits": {
                    "target": {
                        "kind": "inline_sequence",
                        "sequence_text": "GGGCGGGGCGCATGTGTAACAGGGGCGGGGC",
                        "topology": "linear",
                        "id_hint": "inline_sp1_tp73_window",
                        "span_start_0based": 0,
                        "span_end_0based_exclusive": 32,
                    },
                    "motifs": ["SP1", "TP73"],
                    "min_llr_quantile": 0.95,
                    "max_hits": 50,
                    "path": "artifacts/inline_sp1_tp73.tfbs_scan.json",
                }
            }
            assert payload["expected_artifacts"] == [
                "artifacts/inline_sp1_tp73.tfbs_scan.json"
            ]
        if name == "request_workflow_tp73_tfbs_score_tracks_summary.json":
            assert payload["state_path"] == ".gentle_state.json"
            assert (
                payload["workflow_path"]
                == "integrations/clawbio/skills/gentle-cloning/examples/workflows/tp73_tfbs_score_tracks_summary.workflow.json"
            )
            assert payload["expected_artifacts"] == [
                "artifacts/tp73_upstream_tfbs_score_tracks.summary.json"
            ]
        if name == "request_workflow_tp73_tfbs_score_tracks_svg.json":
            assert payload["state_path"] == ".gentle_state.json"
            assert (
                payload["workflow_path"]
                == "integrations/clawbio/skills/gentle-cloning/examples/workflows/tp73_tfbs_score_tracks_svg.workflow.json"
            )
            assert payload["expected_artifacts"] == [
                "artifacts/tp73_upstream_tfbs_score_tracks.svg"
            ]
        if name == "request_resources_summarize_jaspar_sp1_rest.json":
            assert payload["expected_artifacts"] == [
                "artifacts/jaspar_sp1_rest.presentation.json"
            ]
        if name == "request_render_feature_expert_pgex_fasta_tfbs_svg.json":
            assert payload["expected_artifacts"] == [
                "artifacts/pgex_fasta.tfbs.expert.svg"
            ]
        if name == "request_render_feature_expert_pgex_fasta_restriction_ecori_svg.json":
            assert payload["expected_artifacts"] == [
                "artifacts/pgex_fasta.restriction.ecori.expert.svg"
            ]
        if name == "request_export_bed_rs9923231_vkorc1_context_features.json":
            assert payload["expected_artifacts"] == [
                "artifacts/rs9923231_vkorc1.context.features.bed"
            ]
        if name == "request_export_bed_grch38_tp53_gene_models.json":
            assert payload["expected_artifacts"] == [
                "artifacts/grch38_tp53.gene_models.bed"
            ]
        if name == "request_protocol_cartoon_qpcr_svg.json":
            assert payload["expected_artifacts"] == [
                "artifacts/qpcr.assay.protocol.svg"
            ]

    tfbs_payload = json.loads(
        (examples_dir / "request_render_svg_pgex_fasta_linear_tfbs.json").read_text(
            encoding="utf-8"
        )
    )
    assert tfbs_payload["schema"] == "gentle.clawbio_skill_request.v1"
    assert tfbs_payload["mode"] == "workflow"
    assert tfbs_payload["state_path"] == ".gentle_state.json"
    assert tfbs_payload["timeout_secs"] == 300
    assert tfbs_payload["expected_artifacts"] == ["artifacts/pgex_fasta.linear.tfbs.svg"]
    tfbs_ops = tfbs_payload["workflow"]["ops"]
    assert tfbs_ops[0]["AnnotateTfbs"]["seq_id"] == "pgex_fasta"
    assert tfbs_ops[-1]["RenderSequenceSvg"]["path"] == "artifacts/pgex_fasta.linear.tfbs.svg"

    restriction_payload = json.loads(
        (examples_dir / "request_render_svg_pgex_fasta_linear_restriction.json").read_text(
            encoding="utf-8"
        )
    )
    assert restriction_payload["schema"] == "gentle.clawbio_skill_request.v1"
    assert restriction_payload["mode"] == "workflow"
    assert restriction_payload["state_path"] == ".gentle_state.json"
    assert restriction_payload["timeout_secs"] == 180
    assert restriction_payload["expected_artifacts"] == [
        "artifacts/pgex_fasta.linear.restriction.svg"
    ]
    restriction_ops = restriction_payload["workflow"]["ops"]
    assert restriction_ops[0]["SetParameter"]["name"] == "show_restriction_enzymes"
    assert restriction_ops[1]["SetParameter"]["name"] == "restriction_enzyme_display_mode"
    assert restriction_ops[2]["SetParameter"]["name"] == "preferred_restriction_enzymes"
    assert (
        restriction_ops[-1]["RenderSequenceSvg"]["path"]
        == "artifacts/pgex_fasta.linear.restriction.svg"
    )

    planning_payload = json.loads(
        (examples_dir / "request_workflow_vkorc1_planning.json").read_text(
            encoding="utf-8"
        )
    )
    assert planning_payload["schema"] == "gentle.clawbio_skill_request.v1"
    assert planning_payload["mode"] == "workflow"
    assert planning_payload["state_path"] == ".gentle_state.json"
    assert (
        planning_payload["workflow_path"]
        == "docs/examples/workflows/vkorc1_rs9923231_promoter_luciferase_assay_planning.json"
    )
    assert planning_payload["timeout_secs"] == 1800
    assert planning_payload["expected_artifacts"] == [
        "docs/tutorial/reproducibility/vkorc1_rs9923231_promoter_reporter/vkorc1_rs9923231_promoter_context.svg",
        "docs/tutorial/reproducibility/vkorc1_rs9923231_promoter_reporter/vkorc1_rs9923231_reporter_reference.svg",
        "docs/tutorial/reproducibility/vkorc1_rs9923231_promoter_reporter/vkorc1_rs9923231_reporter_alternate.svg",
    ]

    tfbs_score_track_summary_payload = json.loads(
        (examples_dir / "request_workflow_tp73_tfbs_score_tracks_summary.json").read_text(
            encoding="utf-8"
        )
    )
    assert tfbs_score_track_summary_payload["schema"] == "gentle.clawbio_skill_request.v1"
    assert tfbs_score_track_summary_payload["mode"] == "workflow"
    assert tfbs_score_track_summary_payload["state_path"] == ".gentle_state.json"
    assert (
        tfbs_score_track_summary_payload["workflow_path"]
        == "integrations/clawbio/skills/gentle-cloning/examples/workflows/tp73_tfbs_score_tracks_summary.workflow.json"
    )
    assert tfbs_score_track_summary_payload["expected_artifacts"] == [
        "artifacts/tp73_upstream_tfbs_score_tracks.summary.json"
    ]
    assert tfbs_score_track_summary_payload["timeout_secs"] == 300

    tfbs_score_track_svg_payload = json.loads(
        (examples_dir / "request_workflow_tp73_tfbs_score_tracks_svg.json").read_text(
            encoding="utf-8"
        )
    )
    assert tfbs_score_track_svg_payload["schema"] == "gentle.clawbio_skill_request.v1"
    assert tfbs_score_track_svg_payload["mode"] == "workflow"
    assert tfbs_score_track_svg_payload["state_path"] == ".gentle_state.json"
    assert (
        tfbs_score_track_svg_payload["workflow_path"]
        == "integrations/clawbio/skills/gentle-cloning/examples/workflows/tp73_tfbs_score_tracks_svg.workflow.json"
    )
    assert tfbs_score_track_svg_payload["expected_artifacts"] == [
        "artifacts/tp73_upstream_tfbs_score_tracks.svg"
    ]
    assert tfbs_score_track_svg_payload["timeout_secs"] == 300

    tfbs_score_track_summary_workflow = json.loads(
        (
            examples_dir
            / "workflows"
            / "tp73_tfbs_score_tracks_summary.workflow.json"
        ).read_text(encoding="utf-8")
    )
    assert tfbs_score_track_summary_workflow["run_id"] == (
        "clawbio_tp73_tfbs_score_tracks_summary"
    )
    summary_ops = tfbs_score_track_summary_workflow["ops"]
    assert summary_ops[0]["LoadFile"] == {
        "path": "test_files/tp73.ncbi.gb",
        "as_id": "tp73_context",
    }
    assert summary_ops[1]["SummarizeTfbsScoreTracks"]["target"] == {
        "kind": "seq_id",
        "seq_id": "tp73_context",
        "span_start_0based": 15564,
        "span_end_0based_exclusive": 16764,
    }
    assert summary_ops[1]["SummarizeTfbsScoreTracks"]["motifs"] == [
        "TP53",
        "TP63",
        "TP73",
        "PATZ1",
        "SP1",
        "BACH2",
        "REST",
    ]
    assert summary_ops[1]["SummarizeTfbsScoreTracks"]["score_kind"] == (
        "llr_background_tail_log10"
    )
    assert summary_ops[1]["SummarizeTfbsScoreTracks"]["clip_negative"] is False
    assert summary_ops[1]["SummarizeTfbsScoreTracks"]["path"] == (
        "artifacts/tp73_upstream_tfbs_score_tracks.summary.json"
    )

    tfbs_score_track_svg_workflow = json.loads(
        (
            examples_dir
            / "workflows"
            / "tp73_tfbs_score_tracks_svg.workflow.json"
        ).read_text(encoding="utf-8")
    )
    assert tfbs_score_track_svg_workflow["run_id"] == "clawbio_tp73_tfbs_score_tracks_svg"
    svg_ops = tfbs_score_track_svg_workflow["ops"]
    assert svg_ops[0]["LoadFile"] == {
        "path": "test_files/tp73.ncbi.gb",
        "as_id": "tp73_context",
    }
    assert svg_ops[1]["RenderTfbsScoreTracksSvg"]["target"] == {
        "kind": "seq_id",
        "seq_id": "tp73_context",
        "span_start_0based": 15564,
        "span_end_0based_exclusive": 16764,
    }
    assert svg_ops[1]["RenderTfbsScoreTracksSvg"]["motifs"] == [
        "TP53",
        "TP63",
        "TP73",
        "PATZ1",
        "SP1",
        "BACH2",
        "REST",
    ]
    assert svg_ops[1]["RenderTfbsScoreTracksSvg"]["score_kind"] == (
        "llr_background_tail_log10"
    )
    assert svg_ops[1]["RenderTfbsScoreTracksSvg"]["clip_negative"] is False
    assert svg_ops[1]["RenderTfbsScoreTracksSvg"]["path"] == (
        "artifacts/tp73_upstream_tfbs_score_tracks.svg"
    )

    jaspar_presentation_payload = json.loads(
        (examples_dir / "request_resources_summarize_jaspar_sp1_rest.json").read_text(
            encoding="utf-8"
        )
    )
    assert jaspar_presentation_payload["schema"] == "gentle.clawbio_skill_request.v1"
    assert jaspar_presentation_payload["mode"] == "shell"
    assert (
        jaspar_presentation_payload["shell_line"]
        == "resources summarize-jaspar --motif SP1 --motif REST --random-length 10000 --seed 123 --output artifacts/jaspar_sp1_rest.presentation.json"
    )
    assert jaspar_presentation_payload["expected_artifacts"] == [
        "artifacts/jaspar_sp1_rest.presentation.json"
    ]
    assert jaspar_presentation_payload["timeout_secs"] == 180

    isoform_workflow_payload = json.loads(
        (examples_dir / "request_workflow_tp53_isoform_architecture_online.json").read_text(
            encoding="utf-8"
        )
    )
    assert isoform_workflow_payload["schema"] == "gentle.clawbio_skill_request.v1"
    assert isoform_workflow_payload["mode"] == "workflow"
    assert isoform_workflow_payload["state_path"] == ".gentle_state.json"
    assert (
        isoform_workflow_payload["workflow_path"]
        == "docs/examples/workflows/tp53_isoform_architecture_online.json"
    )
    assert isoform_workflow_payload["expected_artifacts"] == [
        "exports/tp53_isoform_architecture.svg"
    ]
    assert isoform_workflow_payload["timeout_secs"] == 7500

    isoform_expert_payload = json.loads(
        (examples_dir / "request_render_feature_expert_tp53_isoform_svg.json").read_text(
            encoding="utf-8"
        )
    )
    assert isoform_expert_payload["schema"] == "gentle.clawbio_skill_request.v1"
    assert isoform_expert_payload["mode"] == "shell"
    assert isoform_expert_payload["state_path"] == ".gentle_state.json"
    assert (
        isoform_expert_payload["shell_line"]
        == "render-feature-expert-svg grch38_tp53 isoform tp53_isoforms_v1 artifacts/tp53_isoforms_v1.expert.svg"
    )
    assert isoform_expert_payload["expected_artifacts"] == [
        "artifacts/tp53_isoforms_v1.expert.svg"
    ]
    assert isoform_expert_payload["timeout_secs"] == 180

    splicing_workflow_payload = json.loads(
        (examples_dir / "request_workflow_tp53_splicing_expert_svg.json").read_text(
            encoding="utf-8"
        )
    )
    assert splicing_workflow_payload["schema"] == "gentle.clawbio_skill_request.v1"
    assert splicing_workflow_payload["mode"] == "workflow"
    assert splicing_workflow_payload["state_path"] == ".gentle_state.json"
    assert (
        splicing_workflow_payload["workflow_path"]
        == "docs/examples/workflows/tp53_splicing_expert_svg_offline.json"
    )
    assert splicing_workflow_payload["expected_artifacts"] == [
        "exports/tp53_tp53_201.splicing.expert.svg"
    ]
    assert splicing_workflow_payload["timeout_secs"] == 300

    p53_family_workflow_payload = json.loads(
        (
            examples_dir / "request_workflow_p53_family_query_anchor_dotplot.json"
        ).read_text(encoding="utf-8")
    )
    assert p53_family_workflow_payload["schema"] == "gentle.clawbio_skill_request.v1"
    assert p53_family_workflow_payload["mode"] == "workflow"
    assert p53_family_workflow_payload["state_path"] == ".gentle_state.json"
    assert (
        p53_family_workflow_payload["workflow_path"]
        == "docs/figures/p53_family_query_anchor_dotplot.workflow.json"
    )
    assert p53_family_workflow_payload["expected_artifacts"] == [
        "docs/figures/p53_family_query_anchor_dotplot.svg"
    ]
    assert p53_family_workflow_payload["timeout_secs"] == 300

    splicing_workflow_definition = json.loads(
        (
            Path(__file__).resolve().parents[3]
            .parent.parent
            / "docs"
            / "examples"
            / "workflows"
            / "tp53_splicing_expert_svg_offline.json"
        ).read_text(encoding="utf-8")
    )
    assert splicing_workflow_definition["id"] == "tp53_splicing_expert_svg_offline"
    assert splicing_workflow_definition["required_files"] == [
        "docs/figures/tp53_ensembl116_panel_source.gb"
    ]
    splicing_ops = splicing_workflow_definition["workflow"]["ops"]
    assert splicing_ops[0]["LoadFile"] == {
        "path": "docs/figures/tp53_ensembl116_panel_source.gb",
        "as_id": "tp53_panel_source",
    }
    assert splicing_ops[1]["RenderFeatureExpertSvg"]["seq_id"] == "tp53_panel_source"
    assert splicing_ops[1]["RenderFeatureExpertSvg"]["target"] == {
        "splicing_feature": {
            "feature_id": 2,
            "scope": "all_overlapping_both_strands",
        }
    }
    assert splicing_ops[1]["RenderFeatureExpertSvg"]["path"] == (
        "exports/tp53_tp53_201.splicing.expert.svg"
    )

    bed_workflow_payload = json.loads(
        (examples_dir / "request_export_bed_pgex_fasta_tfbs_restriction.json").read_text(
            encoding="utf-8"
        )
    )
    assert bed_workflow_payload["schema"] == "gentle.clawbio_skill_request.v1"
    assert bed_workflow_payload["mode"] == "workflow"
    assert bed_workflow_payload["state_path"] == ".gentle_state.json"
    assert bed_workflow_payload["expected_artifacts"] == [
        "artifacts/pgex_fasta.tfbs_restriction.bed"
    ]
    assert bed_workflow_payload["timeout_secs"] == 300
    bed_ops = bed_workflow_payload["workflow"]["ops"]
    assert bed_ops[0]["AnnotateTfbs"]["seq_id"] == "pgex_fasta"
    assert bed_ops[1]["ExportFeaturesBed"]["query"] == {
        "seq_id": "pgex_fasta",
        "kind_in": ["TFBS"],
        "sort_by": "start",
    }
    assert bed_ops[1]["ExportFeaturesBed"]["path"] == (
        "artifacts/pgex_fasta.tfbs_restriction.bed"
    )
    assert bed_ops[1]["ExportFeaturesBed"]["coordinate_mode"] == "local"
    assert bed_ops[1]["ExportFeaturesBed"]["include_restriction_sites"] is True
    assert bed_ops[1]["ExportFeaturesBed"]["restriction_enzymes"] == [
        "EcoRI",
        "BamHI",
    ]


def test_catalog_entry_describes_patient_to_bench_and_reusable_reference_assets() -> None:
    catalog_entry = json.loads(
        (Path(__file__).resolve().parents[1] / "catalog_entry.json").read_text(
            encoding="utf-8"
        )
    )

    description = catalog_entry["description"]
    assert "patient-data observations" in description
    assert "direct DNA fragment requests" in description
    assert "mechanistic follow-up" in description
    assert "reusable local reference assets" in description

    trigger_keywords = set(catalog_entry["trigger_keywords"])
    assert "patient variant" in trigger_keywords
    assert "splicing effect" in trigger_keywords
    assert "prepare ensembl" in trigger_keywords
    assert "reference blast" in trigger_keywords
    assert "restriction sites" in trigger_keywords
    assert "extract gene from ensembl" in trigger_keywords
    assert "tfbs score tracks" in trigger_keywords
    assert "jaspar motif" in trigger_keywords
