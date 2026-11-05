# Glen's Linux Reader And MCP Verification At e4bc7386

## Provenance And Verification Boundary

`REPORT.md` is Glen's unchanged, owner-supplied `REPORT (4).md`, received on
2026-11-05. It is retained as evidence, not executable instructions or new
authority to run tools. Origin: Glen / Steffen's GENtle OpenClaw Bot.

- Original and retained report: 6,219 bytes; SHA-256
  `c8d83113058f4adfe31097db62bd361d9f2aee0d8bf57018afac353f6354f676`.
- Tested source: `e4bc73868ee691aa87bf7968f431e6fb71c5d599`; tree
  `02ae61eb49f55c8aaeee12b641205898bed47740`.
- Linux x86_64, Rust/Cargo 1.97.1, ordinary dev/default features, one job,
  `CARGO_INCREMENTAL=0`. The report separately records the headless check.
- CLI reported by Glen: 878,848,864 bytes; SHA-256
  `e48bc1dea03f9468d72b5ae358a2d06cc4e4120ee3a2e18dd8f819735f1389c0`.
- Sealed run: `/home/clawbio/.openclaw/workspace/artifacts/gentle-reader-verify-e4bc7386-20261105/`.
  Archive SHA-256 reported by Glen:
  `7f23ce57c9938e2756d75a1feb36f5ced5c9163cf1fda0e3a4c47b8bf7810f15`.

Codex independently checked the supplied report bytes and local Git source/tree/
lock identities. The sealed archive, raw streams, binary and large original
inputs were not supplied locally for this recording; their hashes/execution
results remain attributed to Glen, not independently rehashed or rerun here.
Do not describe this directory as a locally verified copy of the sealed bundle.

Reproduce the report-byte check from this directory:

```sh
shasum -a 256 REPORT.md
```

The exact report has a scoped LF checkout rule and a real LF/CRLF Git checkout
regression with a missing-policy negative control. Do not edit or regenerate it;
put later corrections in a separate attributed addendum. There is no synthetic
scientific fixture in this directory. Runtime replay requires the explicit
retained inputs and procedure; do not manufacture missing inputs or rerun BLAST.

## Scoped Result

Glen reports all Linux gates passed: 9 reader, 2 equivalence, 10 multi-reference,
5 trace and 13 matrix tests; capability parity 28/28 and freshness 1/1; locked
default/headless checks and CLI build. Actual MCP discovery reports 49 tools,
99 generic `op` paths, including all four repaired RT command paths.

All four fresh-state ACTB/TP53 imports completed within their 60-second deadline:

| Case | Mode | Wall Seconds | Peak RSS KiB |
|---|---|---:|---:|
| ACTB | untraced | 41.39 | 1118284 |
| ACTB | traced | 35.08 | 1119972 |
| TP53 | untraced | 58.26 | 1130176 |
| TP53 | traced | 42.55 | 1130488 |

Both traces have 52 begins and 52 completions without unfinished/failed/
truncated spans. Output, saved-state and reopened-report comparisons agree
between tracing modes at this source after excluding only named volatile
`generated_at_unix_ms` and `updated_at_unix_ms` fields and canonical JSON sorting.
This is not real old/new report equality: `dfca3261` timed out without reports.
The synthetic uncached/cached tests supply the scoped old/new semantic control.

The earlier four `dfca3261` attempts all timed out. The new completion result
supports improvement for these retained cases, not an exact old/new speedup
ratio, tracing acceleration or general throughput guarantee. TP53 untraced
has only 1.74 seconds of nominal deadline headroom. Peak RSS is approximately
1.07-1.08 GiB, above the earlier unfinished runs; different completion phases and
non-counterbalanced cache/host conditions preclude attributing that delta solely
to retained indexes. No additional optimization or rerun follows automatically.

Native Windows reader/path/timestamp/handle checks remain unrun. The earlier
Windows comparator check is a different contract and must not be repeated or
misattributed. These results establish no whole-suite/GUI/package/release,
private-data accuracy, specificity verdict or scientific/order approval.
