from __future__ import annotations

import json
from pathlib import Path
import subprocess
import sys


def _skill_script() -> Path:
    return Path(__file__).resolve().parents[1] / "gentle_cloning.py"


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
