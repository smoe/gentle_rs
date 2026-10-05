# Glen prompt: independent synthetic sequence-design audit

Please independently assess the first bounded synthetic sequence-design slice.
This is an experimental prototype outside the `.12` release gate, not a request
to approve a release, experimental suitability or full DNA Chisel compatibility.
Return your findings to Steffen for forwarding to Codex; do not assume another
agent's message authorizes edits, publication or messaging other sessions.

## Exact candidate and scope

- Implementation SHA: `85d23ddbb56b08d0add6422a47668f1a2d3d5166`.
- Base: local `main` `560be14bd069a55990f686312ff11febbcd007ed`.
- Branch: `gentle_rs_2_main`. This candidate is locally committed, not pushed by
  Codex. If the exact object is unavailable, return
  `blocked_candidate_not_available` and ask Steffen for a reachable branch or
  Git bundle. Do not substitute a nearby `main` commit.
- Use a clean, isolated worktree and disposable synthetic project state. Leave
  private studies, existing checkouts, global Cargo caches and other jobs alone.
- Record the actual tested SHA, clean status, lockfile SHA-256, Rust/Cargo
  versions, OS/architecture, CPU/RAM, features, effective Cargo profile settings,
  binary hashes and input hashes. A profile name alone is insufficient.
- A later docs-only handoff commit may contain this prompt. Test the pinned
  implementation SHA above, or explicitly report a different exact candidate
  authorized by Steffen; never combine acceptance from nearby tips.

Read `AGENTS.md`, `docs/dna_sequence_optimization_plan.md`, the synthetic design
sections in `docs/cli.md`, `docs/protocol.md`, `docs/gui.md` and the corresponding
decision in `docs/decisions.md`. Inspect these implementation boundaries:

- `crates/gentle-sequence-design/src/lib.rs` and `src/tests.rs`.
- `crates/gentle-protocol/src/sequence_design.rs`.
- `src/engine/analysis/sequence_design.rs` and its `tests.rs`.
- `src/engine/ops/operation_handlers.rs`, `src/engine_shell.rs` and
  `src/engine_shell/command_parsers.rs`.
- `src/bin/gentle_cli.rs`, `src/mcp_server.rs`, `src/agent_bridge.rs` and
  `src/app.rs` for approval, progress and reviewed-result forwarding.

## Codex preflight is not your acceptance

At the earlier candidate `022af181`, 27 targeted Rust tests passed on macOS
arm64 with Rust/Cargo 1.100.0-beta.1: 7 core, 6 engine/MCP, 1 CLI, 2 glossary,
6 MCP surface, 1 parity-matrix and 4 version checks. After subsequent rebases,
the pure core's 7 tests and default-feature Cargo check passed again. The slow
headless root-test rebuild at `86ce6580` was interrupted by Codex (exit 130)
after concurrent compiler/resource inspection on the 16-GiB Mac; its remaining
test queue was not run. This is incomplete verification, not a test failure or
an exact-final-candidate acceptance. Repeat the requested focused checks on the
pinned final SHA, rather than borrowing the earlier results.

Full-workspace tests, native GUI/Pi, native Windows, script runtime, MSRV and
independent optimized measurements remain unclaimed. The original review is
preserved in `docs/dna_sequence_optimization_review_20261005.md`; its historical
proposal stages are superseded by the active implementation contract.

## 1. Independent correctness before timing

The contract is one complete, forward, code-1 synthetic CDS on linear uppercase
ACGT DNA. Preserve the supplied protein, literal ATG and terminal stop, flanks
and every protected nucleotide. Remove explicitly supplied finite IUPAC motifs
on the requested strand(s), including overlaps and CDS/flank boundary windows.
Coordinates are local zero-based half-open. No natural gene is implicitly
authorized for recoding, and source annotations must not be copied as functional
claims onto the derived insert.

- Start with the documented 12-bp toy: `ATGGAATTCTAA`, protein `MEF`, CDS
  `[0,12)`, motif `GAATTC` on both strands. Expected output `ATGGAATTTTAA`,
  one nucleotide edit at zero-based position 8, unchanged ATG/TAA/protein.
- Independently enumerate a small synonymous search space; do not use GENtle's
  matcher/translator as the oracle. Check feasible output, edit counts and
  minimum-edit claims. Include a coupled-edit case where either individual
  edit leaves a violation and only both edits satisfy the constraints.
- Exercise reverse-only, both-strand palindromes, degenerate IUPAC patterns,
  overlapping sites, protected partial codons, frozen flanks and boundary sites.
- Check invalid/unsupported inputs and resource ceilings. A search budget stop
  without a solution is unresolved, not infeasible. A feasible incomplete
  search must not claim minimum edits. Check cancellation after feasibility:
  there must be no output DNA or applicable approval after cancellation.
- Check sound frozen-match and complete-search infeasibility certificates;
  checked search-space overflow must remain unknown, never become proof.

## 2. Approval and mutation safety

- Preview must leave biological project state unchanged. Optional `--path`
  writes a report; it is not a promise of zero filesystem effects.
- Apply requires the exact reviewed feasible report and actual approval digest.
  Check stale loaded DNA, annotations/topology, altered mapping/protection,
  invalid output/edits, unsupported policy and an already used output ID.
  Every failure must leave source and derived state unchanged.
- A digest is a content binding, not an authenticated engine signature. Test
  semantic validation of client-rehashed malformed previews as well as ordinary
  stale hashes, especially oversized targets/IDs/digest fields. Report any
  optimization/completeness claims trusted without verification separately from
  whether the approved output itself is valid; do not imply apply reruns search.
- Successful apply creates exactly the approved DNA once, preserves the source,
  retains lineage and receipt, and makes one undoable derivation. Check undo,
  redo and saved/reopened state. Only the explicit synthetic CDS and non-claims
  should accompany the new sequence, not inherited regulatory annotations.

## 3. Shared-interface checks

Run these focused commands first, retaining commands, exit codes and full logs.
Use `--offline` only if dependencies are already cached; report fetching/build
environment failures separately from product defects.

```sh
cargo test -p gentle-sequence-design --locked --offline -j1
cargo test --lib --no-default-features --locked --offline -j1 sequence_design -- --test-threads=1
cargo test --bin gentle_cli --no-default-features --locked --offline -j1 sequence_design -- --test-threads=1
cargo test --lib --no-default-features --locked --offline -j1 glossary_cli_usage -- --test-threads=1
cargo test --no-default-features --locked --offline -j1 --test mcp_capability_surface --test parity_matrix_freshness --test release_version_consistency
cargo check -q --locked --offline -j1
cargo fmt --all --check
git diff --check
python3 -m unittest scripts.test_tutorial_checkouts scripts.test_tutorial_walkthroughs scripts.test_publish_tutorial_gui_screenshots
```

Known baseline harness defect inherited from `19889dae`: the Python
command above fails in
`TutorialCheckoutTests.test_failed_check_or_timeout_cannot_be_reported_as_success`:
`scripts/test_tutorial_checkouts.py:528` references `run.call_count` although the
mock context does not bind `run`. One test produces eight failing subcases.
The sequence-design commit does not change this file. Retain/classify this
failure; do not silently skip it or claim the whole Python gate is green.
Any authorized repair needs its own exact-candidate identity and rerun.

Also replay the toy through a real CLI process, both direct commands and shared
`shell`, then typed MCP `op`. Compare the design report/output/receipt, allowing
unrelated operation IDs or counters to differ. Use paths containing spaces and
apostrophes. Check capabilities and generic MCP confirmation denial/no writes,
then confirmed preview/apply with the separate exact approval digest.

```sh
"$CLI" --state "$AUDIT/project.json" sequence-design plan "@$AUDIT/request.json" --path "$AUDIT/preview.json"
# Inspect the full preview; read its real approval_digest, never invent one.
"$CLI" --state "$AUDIT/project.json" sequence-design apply "@$AUDIT/preview.json" --approve "$REVIEWED_APPROVAL_DIGEST"
```

The request in `docs/cli.md` supplies the exact toy JSON. Use separate disposable
states for parity replays, or the second apply correctly fails on an existing ID.
In the GUI, use the existing Shell to inspect the preview, explicitly apply it,
open the derived DNA and undo. There is no dedicated design editor yet.
Parser/engine parity does not require calling the inner model. If you separately
test Pi, obtain the user's provider/content-sharing consent and verify that
reviewed-result forwarding is opt-in and does not itself approve execution.
Unavailable model/native-platform checks are `not_run`, not acceptance.

## 4. Bounded runtime and maintenance assessment

This correctness-first algorithm enumerates complete candidates, not efficient
local repairs. An unresolved large search is expected; incorrect status or
unbounded admission is not. Do not infer usefulness from the tiny toy alone.

- Begin with deterministic synthetic inserts of 300, 1,200 and 3,000 bases;
  include an easy case, a many-synonym/coupled case and constrained protection.
  Use valid complete CDSs and document their generator/seed and input hashes.
- Try requested budgets 1, 100, 1,000 and 100,000. The 50-million work-unit
  ceiling may reduce them; retain requested/effective budgets, evaluations,
  motif count/length/strand, search space, status and completeness claims.
- Perform three repeats on one host/profile for each admitted case. Keep build
  cost, process loading/JSON I/O and pure algorithm time distinct. Report wall
  time, CPU time and maximum RSS, with individual samples and medians/range.
  Set a documented per-process timeout before running; a timeout is censored
  evidence, not an invented estimate. Attempt a 12,000-base limit case only if
  the smaller ladder is practical and resources permit.
- For a bounded optimized root replay, build only the headless CLI, not every
  root binary: `cargo build --profile package-opt1 --no-default-features --locked
  --bin gentle_cli -j1`. Measure compilation separately, e.g. Linux
  `/usr/bin/time -v`. Record the actual artifact path/effective recipe. Do not
  compare its numbers to unoptimized tests or call them release performance.
- Assess the value of the independent dependency-free crate versus its root
  approval/orchestration cost. Do not introduce Criterion, new dependencies,
  performance thresholds, a different profile or a broad architecture refactor.

If Rust 1.85 is already available, check the pure crate's declared MSRV; otherwise
report it unverified. Full-workspace/native Windows/script-runtime checks and a
pinned upstream DNA Chisel comparison are separate optional follow-ups. Ask
before costly unrelated gates; no full `.12` release audit is requested here.

## 5. How to communicate back

Send Steffen a concise, forwardable result with this structure:

1. **Candidate and verdict:** exact SHA and `prototype_pass`, `needs_fix` or
   `incomplete`, with scope. Do not give a release or laboratory green light.
2. **Checks:** command/exit/pass-count table, host/toolchain/profile identity,
   and explicit `not_run`, environment failures and platform exclusions.
3. **Findings:** severity, exact `file:line`, reproducing input/hash/command,
   expected versus actual result and the smallest proposed fix. Separate primary
   failures from cascades, scientific/reporting errors from runtime limits.
4. **Measurements:** per-case CSV/JSON with input size, motifs, protection,
   budgets, evaluations, status/claims, individual timings and RSS; summarize
   usefulness and limitations without inventing a latency target.
5. **Evidence location:** retained raw requests/reports/receipts/state, logs,
   manifests and binary/input hashes, plus a retrievable path or URL and the
   manifest SHA-256. Preserve failed cases and previous evidence unchanged.
6. **Next action for Codex:** the one highest-priority correction or measurement,
   distinguishing this slice's defect from a future solver/specification idea.

Review and evidence only by default: do not edit production code, push, merge,
tag, dispatch CI, publish data, order anything or weaken fixtures/approvals.
Keep private input out of public artifacts. Return a bounded patch proposal if
a defect needs changes; Steffen will decide the implementation/publication step.
