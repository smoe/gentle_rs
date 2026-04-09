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
    assert args == [
        "exec",
        "--bind",
        f"{tmp_path}:/work",
        "--pwd",
        "/work",
        str(image_path),
        "gentle_cli",
        "capabilities",
    ]
