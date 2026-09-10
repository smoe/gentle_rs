# Discover MCP Tools and Verify a Shared GENtle Result

> Type: hand-written command-line/MCP walkthrough. Offline.
> Audience: a biologist and the person connecting their external assistant.
> Last updated: 2026-09-10.

The question is deliberately small: **where do EcoRI and BamHI recognize this
12-base synthetic DNA, `GAATTCGGATCC`?** We will ask GENtle through MCP, then
compare its scientific report with the command-line and shared-shell routes.

No AI model, subscription, API key, database or internet connection is needed.
MCP transports requests; GENtle supplies the sequence analysis. Allow about
10 minutes after building. A full rebuild can take substantially longer.

## 1. Prepare an Isolated Run

From the repository root, build the two binaries from the same revision:

```sh
cargo build --offline --locked --bin gentle_cli --bin gentle_mcp
CLI="$PWD/target/debug/gentle_cli"
MCP="$PWD/target/debug/gentle_mcp"
"$CLI" --version
"$MCP" --version
RUN=$(mktemp -d)
```

Use installed binaries instead when both report the same build and source
revision. If offline compilation reports a missing dependency, stop and arrange
an explicitly authorized dependency download; this exercise itself needs none.

Never use a working project's state file. An MCP server reads the file given
by `--state`; it does not see a running GUI's unsaved project. In the current
implementation, generic MCP **`op` saves that file even for a read-only engine
inspection**, and therefore requires `confirm=true`. The helper below creates
only disposable state and evidence under the new output directory.

## 2. Discover Before Calling

An MCP-capable client starts `gentle_mcp --state PATH`. The teaching helper
below uses GENtle's current `Content-Length` framing over standard input/output.
The following are message **bodies**, not commands to paste at a shell prompt:

```json
{"jsonrpc":"2.0","id":1,"method":"initialize","params":{"protocolVersion":"2025-06-18"}}
```

After the initialization response, send `notifications/initialized`, then:

```json
{"jsonrpc":"2.0","id":2,"method":"tools/list","params":{}}
```

Inspect `result.tools`: find `op`, its `inputSchema`, and confirmation/mutation
metadata. Do not assume every shell command has a separately named MCP tool.
In the parity documentation, **`mcp-op` is a surface label, not a tool or CLI
command**. A row with that label is reached via `tools/call`, tool name `op`,
and the exact typed engine operation. The list changes as GENtle evolves; no
fixed tool count is an acceptance criterion.

## 3. Review the Operation and the Refusal

Open the small [operation JSON](./inputs/mcp_restriction_operation.json).
It selects `FindRestrictionSites`, the inline DNA, linear topology, and exactly
EcoRI/BamHI. It has no remote source or output path. The selected span is the
whole inline sequence.

The MCP call has this shape (the helper inserts the complete operation object):

```text
tools/call
  name: op
  arguments.operation: contents of mcp_restriction_operation.json
  arguments.confirm: false, then true after review
```

The first call must return `result.isError=true`, name the missing confirmation,
and leave its state file absent. The second authorizes exactly this local
operation and its state-file write. `confirm=true` is an explicit execution
choice, not authentication and not biological approval. A process exit code of
zero alone does not mean success: check both JSON-RPC `error` and tool `isError`.

## 4. Run the Complete Example

The adjacent Python helper is ordinary standard-library transport/comparison
code. It does not scan DNA or independently calculate enzyme positions. It
performs discovery, the refused call, the confirmed call, and both CLI paths:

```sh
python3 docs/tutorial/inputs/mcp_roundtrip.py \
  --cli "$CLI" --mcp "$MCP" --output-dir "$RUN/roundtrip"
```

Inspect these retained files:

| File | What it tells you |
| --- | --- |
| `discovery.json` | Actual tools and schemas returned by this binary |
| `without-confirmation.json` | Explicit refusal, not a successful analysis |
| `confirmed-operation.json` | MCP result envelope with the engine report |
| `restriction-report.json` | The shared scientific report compared across routes |
| `processes.json` and `*.stdin/stdout/stderr` | Exact invocations, framed bytes, exits and output hashes |
| `summary.json` | Build/binary/input/report identity and narrowly scoped parity verdict |

The helper refuses an existing output directory. It stops on disagreement or a
failed tool; missing or failing steps never become a passing `summary.json`.
This is a teaching replay, not the broader
[candidate-bound Linux tutorial release gate](../testing.md#62-candidate-bound-tutorial-gate-for-external-auditors).

## 5. Check the Biology and the Adapter Boundary

Expected recognition spans in this synthetic input are:

| Enzyme | Bases | 1-based inclusive | 0-based half-open |
| --- | --- | --- | --- |
| EcoRI | GAATTC | 1..6 | [0, 6) |
| BamHI | GGATCC | 7..12 | [6, 12) |

These are recognition intervals, not the staggered cut positions. Both motifs
are palindromic, so do not expect a second duplicate reverse-strand hit.
This tiny test proves neither that all restriction-enzyme cut offsets are
correct nor that a real digest will work under particular reaction conditions.

To see the equivalent shared-shell request yourself:

```sh
"$CLI" --state "$RUN/separate-shell.gentle.json" shell \
  'features restriction-scan --sequence-text GAATTCGGATCC --topology linear --id-hint mcp_tutorial --enzyme EcoRI --enzyme BamHI'
```

Or invoke the same typed operation directly:

```sh
"$CLI" --state "$RUN/separate-cli.gentle.json" op \
  @docs/tutorial/inputs/mcp_restriction_operation.json
```

The outer envelopes differ. The helper compares MCP
`result.content[0].text` (parsed JSON), then its `result.restriction_site_scan`,
with the CLI operation's `restriction_site_scan` and the shell's `report`.
Only `generated_at_unix_ms`, `op_id`, `run_id` and `report_id` are omitted from
comparison: they describe separate invocations, and this report's ID incorporates
the operation/run IDs. Every other field, including target, coordinates and cut
geometry, must agree. Raw reports retain all metadata. Paths, envelopes and
project history are not interchangeable across interfaces; byte-identical
output from separate runs is not promised.

The shell command can also run in the GUI Shell, but doing that is a separate
manual check. MCP is not a general `shell` tool or a GUI remote control.
This example uses the CLI's shared-shell route to demonstrate reachability
parity without inventing an MCP tool that is not advertised.

## Success, Failure and Next Steps

Success requires discovery, an explicit refusal, a successful confirmed call,
the expected two recognition spans, and equality of all three reports. The
offline test runs the real binaries and asserts those coordinates as well:

```sh
GENTLE_TUTORIAL_BIN_DIR="$PWD/target/debug" \
  python3 -m unittest scripts.test_tutorial_walkthroughs.McpWalkthroughTests
```

With no binary directory configured the binary-dependent tests skip; a skip is
not acceptance. Once configured, a missing or failing binary is an error.
Retain the new run directory and exact revision with any report to Glen. No
screenshots or live GUI actions are claimed by this chapter.

Continue with the [Agent Interfaces guide](./01-01_agent_interfaces.md) for
reviewed agent proposals, or the existing
[PATZ1 endpoint/SYBR chapter](./generated/chapters/04-06_patz1_transcript_assay_panels_cli.md)
for transcript-aware primer work. Designing a panel there is not an order-ready
verdict: its specificity and evidence gates still apply.
