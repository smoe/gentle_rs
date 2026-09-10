"""Offline transport demonstration for tutorial 01.03, not a scientific scorer.

Uses the adjacent hand-crafted operation fixture. GENtle computes every report;
this helper frames MCP messages, retains process output, and compares reports.
Run from the repository root with explicit same-build CLI/MCP binaries and a new
output directory. No third-party Python packages, network or model are used.
"""

import argparse
import hashlib
import json
from pathlib import Path
import subprocess


def encode_frames(messages):
    result = bytearray()
    for message in messages:
        body = json.dumps(message, separators=(",", ":")).encode("utf-8")
        result.extend(f"Content-Length: {len(body)}\r\n\r\n".encode("ascii"))
        result.extend(body)
    return bytes(result)


def decode_frames(raw):
    messages = []
    while raw:
        header, separator, body = raw.partition(b"\r\n\r\n")
        if not separator or not header.startswith(b"Content-Length: "):
            raise ValueError("Expected GENtle's Content-Length framed MCP output")
        length = int(header.removeprefix(b"Content-Length: "))
        if length <= 0 or length > len(body):
            raise ValueError("Incomplete MCP response body")
        messages.append(json.loads(body[:length]))
        raw = body[length:]
    return messages


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def tool_payload(response):
    require("error" not in response, f"MCP protocol error: {response}")
    result = response["result"]
    require(not result.get("isError", False), f"MCP tool failed: {result}")
    return json.loads(result["content"][0]["text"])


def comparable_report(report):
    # Keep raw reports unchanged; only execution metadata differs by invocation.
    return {key: value for key, value in report.items()
            if key not in {"generated_at_unix_ms", "op_id", "run_id", "report_id"}}


def run_demo(cli, mcp, output):
    cli, mcp, output = Path(cli).resolve(), Path(mcp).resolve(), Path(output).resolve()
    require(cli.is_file() and mcp.is_file(), "Build both GENtle binaries first")
    output.mkdir(parents=True, exist_ok=False)
    receipts = []

    def run(name, argv, stdin=b""):
        (output / f"{name}.stdin").write_bytes(stdin)
        completed = subprocess.run(argv, input=stdin, capture_output=True, timeout=60)
        (output / f"{name}.stdout").write_bytes(completed.stdout)
        (output / f"{name}.stderr").write_bytes(completed.stderr)
        receipts.append({"step": name, "argv": [str(arg) for arg in argv],
                         "exit_code": completed.returncode,
                         "stdout_sha256": hashlib.sha256(completed.stdout).hexdigest(),
                         "stderr_sha256": hashlib.sha256(completed.stderr).hexdigest()})
        (output / "processes.json").write_text(json.dumps(receipts, indent=2) + "\n")
        require(completed.returncode == 0, f"{name} failed; see retained stderr")
        return completed.stdout

    cli_version = run("cli-version", [cli, "--version"]).decode()
    mcp_version = run("mcp-version", [mcp, "--version"]).decode()
    require(cli_version == mcp_version, "Use CLI and MCP binaries from the same build")
    state = output / "mcp.gentle.json"
    initialize = {"jsonrpc": "2.0", "id": 1, "method": "initialize",
                  "params": {"protocolVersion": "2025-06-18"}}
    initialized = {"jsonrpc": "2.0", "method": "notifications/initialized"}

    def exchange(name, calls):
        # Each subprocess starts a fresh MCP session against the same explicit file.
        messages = [initialize, initialized, *calls]
        responses = decode_frames(run(name, [mcp, "--state", state], encode_frames(messages)))
        require([item.get("id") for item in responses] == [1, *[c["id"] for c in calls]],
                "Missing, duplicate or reordered MCP responses")
        require(responses[0]["result"]["protocolVersion"] == "2025-06-18",
                "Unexpected negotiated protocol")
        (output / f"{name}.json").write_text(json.dumps(responses, indent=2) + "\n")
        return responses[1:]

    listed = exchange("discovery", [{"jsonrpc": "2.0", "id": 2,
                                     "method": "tools/list", "params": {}}])[0]
    tools = {tool["name"]: tool for tool in listed["result"]["tools"]}
    require("op" in tools and "operation" in tools["op"]["inputSchema"]["properties"],
            "The generic typed operation route was not advertised")
    fixture = Path(__file__).with_name("mcp_restriction_operation.json")
    operation = json.loads(fixture.read_text())

    def call(confirm):
        return {"jsonrpc": "2.0", "id": 3, "method": "tools/call",
                "params": {"name": "op", "arguments": {
                    "operation": operation, "confirm": confirm}}}

    denied = exchange("without-confirmation", [call(False)])[0]
    require(denied["result"]["isError"] is True, "Expected confirmation refusal")
    require("confirm=true" in denied["result"]["content"][0]["text"],
            "Refusal did not identify the missing confirmation")
    require(not state.exists(), "Rejected operation unexpectedly created project state")
    accepted = tool_payload(exchange("confirmed-operation", [call(True)])[0])
    require(state.is_file(), "MCP op should retain its designated state file")

    # Different adapter envelopes are not scientific differences: compare the report.
    cli_result = json.loads(run("cli-operation", [cli, "--state", output / "cli.gentle.json",
                                                   "op", "@" + str(fixture)]))
    command = ("features restriction-scan --sequence-text GAATTCGGATCC --topology linear "
               "--id-hint mcp_tutorial --enzyme EcoRI --enzyme BamHI")
    shell_result = json.loads(run("cli-shell", [cli, "--state", output / "shell.gentle.json",
                                               "shell", command]))
    report = accepted["result"]["restriction_site_scan"]
    require(comparable_report(report) == comparable_report(cli_result["restriction_site_scan"])
            == comparable_report(shell_result["report"]),
            "MCP/CLI/shared-shell restriction reports differ")
    spans = [(row["enzyme_name"], row["recognition_start_0based"],
              row["recognition_end_0based_exclusive"]) for row in report["rows"]]
    require(spans == [("EcoRI", 0, 6), ("BamHI", 6, 12)],
            "Unexpected synthetic recognition spans")
    (output / "restriction-report.json").write_text(json.dumps(report, indent=2) + "\n")
    summary = {"status": "pass", "scope": "offline MCP/CLI restriction-report parity",
               "version": cli_version.strip(),
               "cli_sha256": hashlib.sha256(cli.read_bytes()).hexdigest(),
               "mcp_sha256": hashlib.sha256(mcp.read_bytes()).hexdigest(),
               "operation_sha256": hashlib.sha256(fixture.read_bytes()).hexdigest(),
               "report_sha256": hashlib.sha256((output / "restriction-report.json").read_bytes()).hexdigest(),
               "comparison_omits_execution_fields": ["generated_at_unix_ms", "op_id", "run_id", "report_id"],
               "gui_acceptance": "not_run", "biological_validation": "not_performed"}
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cli", required=True)
    parser.add_argument("--mcp", required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()
    run_demo(args.cli, args.mcp, args.output_dir)
    print(f"Shared reports agree. Evidence: {Path(args.output_dir).resolve()}")
