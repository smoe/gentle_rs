# Exact-SHA reader repair and bounded import verification

Executed on 2026-11-05 for exact source
`e4bc73868ee691aa87bf7968f431e6fb71c5d599` (tree
`02ae61eb49f55c8aaeee12b641205898bed47740`, parent
`add274d78e92fe50a8a5f94c19594d06e2dd2b51`) in a clean detached worktree.
No source fix, package build, GUI run, BLAST rerun, reference preparation,
policy relaxation, publication, comparator campaign, or gel campaign was run.

## Identity and profile

- Linux x86_64, kernel `7.2.7+deb14-amd64`; Rust/Cargo 1.97.1.
- `Cargo.lock` SHA-256:
  `7ad5b2dc9485675a75d272fa945c75342dceb60ef582f2885015f3e14df41e76`.
- Explicit run-owned target:
  `/tmp/gentle-reader-verify-e4bc7386-target`.
- One job, `CARGO_INCREMENTAL=0`, ordinary dev/default-feature profile, no
  optimization override.
- CLI: 878,848,864 bytes; SHA-256
  `e48bc1dea03f9468d72b5ae358a2d06cc4e4120ee3a2e18dd8f819735f1389c0`.
  `--version` reports
  `0.1.0-internal.13+git.e4bc73868ee691aa87bf7968f431e6fb71c5d599`.

## Requested Linux gates

All commands exited zero. Complete stdout, stderr, GNU-time records, and exit
codes are retained as `gate-01` through `gate-09`.

- `cargo check -q --locked --tests -j1`: passed; 11:12.00 wall,
  8,419,880 KiB peak RSS.
- `cargo check -q --locked --no-default-features --lib --bin gentle_cli -j1`:
  passed; 4:25.56 wall, 5,920,384 KiB peak RSS.
- `prepared_sequence_reader_`: exactly 9 passed, 0 failed, 0 ignored.
- `primer_specificity_local_reader_`: exactly 2 passed, 0 failed, 0 ignored.
- `engine::primer_specificity_multi::tests::`: 10 passed, 0 failed,
  0 ignored.
- `primer_import_trace`: 5 passed, 0 failed, 0 ignored.
- `specificity_matrix_`: 13 passed, 0 failed, 0 ignored.
- Capability-registry parity: 28 passed, 0 failed, 0 ignored.
- Parity-matrix freshness: 1 passed, 0 failed, 0 ignored.
- `cargo build --locked --bin gentle_cli -j1`: passed.

The nine reader tests include LF/CRLF, Windows-path aliases, raw-input hashes,
error/handle closure, immutability, missing/corrupt input, and reader reuse.
The two equivalence tests preserve hits, warnings, products, summaries, and
reference warnings.

## MCP discovery

The parity test exercised the actual `tools/list` projection. A separate
framed stdio `tools/list` capture reported 49 tools and 99 command paths on the
generic `op` tool. All four repaired paths were present:

- `rna-reads primer-junctions`
- `primers assess-rt-candidates`
- `primers compare-rt-candidates`
- `primers redesign-rt`

## Admission

- Original retained archive: 103,172,261 bytes; SHA-256
  `35e000d086f17348ea76dbaf88712f44b808abea045a72509cabc6d6034a1dcd`.
- All 213 archive manifest hashes passed before and after the attempts.
- ACTB retained inputs: 5,999 BLAST rows / 686,525 bytes.
- TP53 retained inputs: 12,202 BLAST rows / 1,404,773 bytes.
- The moved Ensembl-116 cDNA index was restored byte-identically from its
  verified Storage Box archive: 10/10 manifest-bound files. No reference was
  built or rebound.
- Pre/post exact-source reference inspections were byte-identical:
  - cDNA fingerprint
    `sha256:3159519891bd8c0d4690e9c2f3de8eeb9621152b58c749e3d7d4920f4913263f`,
    465,769 sequences; annotation SHA-1
    `5f571d5775368aeb4e8a88c3ce4b1c8c138e3acf`.
  - GRCh38 fingerprint
    `sha256:90a64ae8ff24fdbab67286bab563c95a6b2a04be37f50ad811474106c1ba3e7e`,
    706 contigs.

## Bounded imports

Each attempt used a distinct initially absent state, a 60-second TERM deadline
and five-second kill grace. All four completed with exit 0, a saved state and
a full report. No run-owned process survived. Total wall time was 177.28
seconds, below the 300-second budget.

| Case | Trace | Exit | Wall | User | System | Peak RSS |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| ACTB | unset | 0 | 41.39 s | 31.51 s | 1.59 s | 1,118,284 KiB |
| ACTB | `1` | 0 | 35.08 s | 31.38 s | 1.45 s | 1,119,972 KiB |
| TP53 | unset | 0 | 58.26 s | 34.85 s | 2.03 s | 1,130,176 KiB |
| TP53 | `1` | 0 | 42.55 s | 34.39 s | 1.80 s | 1,130,488 KiB |

The wall-time differences are not interpreted as tracing speedups: all four
have similar user CPU time, and cache/filesystem conditions were not
counterbalanced.

Both traced attempts emitted 104 valid
`gentle.primer_import_trace.v1` records: 52 `begin` and 52 `completed`, with
no failed, abandoned, truncated, or unmatched spans. Traced import duration
was 33.576266 seconds for ACTB and 39.483314 seconds for TP53. The two
reference-specific `forward_hits` durations were 2.858517/0.301523 seconds
for ACTB and 5.134454/0.377155 seconds for TP53. Reverse-hit and amplicon-join
spans also completed.

This differs materially from `dfca3261`, where all four attempts timed out at
60 seconds and both traces remained inside their first `forward_hits` span
after entering it at about 22.6 seconds. The old residual deadline figures
are historical wall-budget observations, not CPU attribution or strict span
lower bounds.

## Scientific equivalence and reopening

All four saved states were reopened with
`primers specificity-multi-list`; all four reopen commands exited zero.
Baseline and traced output, state, and reopened-report payloads were compared
after removing exactly these volatile fields wherever present:

- `generated_at_unix_ms`
- `updated_at_unix_ms`

All six baseline/trace comparison pairs were byte-identical after canonical
JSON sorting and that normalization. No other differing field was excluded.
The normalized SHA-256 values are retained in
`normalized-content-sha256.txt`.

## Native Windows boundary

No native Windows device was connected to this gateway at verification time,
so the reader repair's native-Windows boundary was not run here. Linux reader
results do not certify Windows LF/CRLF, path, timestamp, handle, symlink, or
hard-link behavior. No remote workflow was dispatched, as requested.

## Verdict

The exact revision passes every requested Linux source gate, preserves the
full parity assertions, exposes all four repaired RT command paths in MCP
discovery, and completes both retained real imports in traced and untraced
modes with scientifically equivalent content. Native-Windows reader evidence
remains outstanding and is reported separately rather than inferred.
