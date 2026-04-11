from __future__ import annotations

import json
import os
from pathlib import Path
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


def test_bootstrap_example_requests_cover_status_and_prepare_routes() -> None:
    examples_dir = Path(__file__).resolve().parents[1] / "examples"
    expected = {
        "request_genomes_list_human.json": (
            "shell",
            "genomes list --filter human",
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
    }

    for name, (mode, shell_line, timeout_secs) in expected.items():
        payload = json.loads((examples_dir / name).read_text(encoding="utf-8"))
        assert payload["schema"] == "gentle.clawbio_skill_request.v1"
        assert payload["mode"] == mode
        assert payload["shell_line"] == shell_line
        assert payload["timeout_secs"] == timeout_secs
