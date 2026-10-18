# Sequence-Design Repair And GC Acceptance

Owner-requested gap closure, 2026-10-13. This supersedes the execution scope of
the [first-slice handoff](glen_sequence_design_handoff_20261005.md), not its
historical evidence. The [active plan](dna_sequence_optimization_plan.md) and
[durable contract](decisions.md#synthetic-sequence-design-boundaries) remain
authoritative. No `.12` gate, release, biological suitability or performance
verdict follows from this document or the added test apparatus.

## Candidate And Authority

Have the main-maintaining thread select one authorized, integrated commit that
contains these checks. Use a clean isolated checkout of that exact SHA, not a
moving `main`. Record lockfile hash, OS/CPU/RAM, Rust/Cargo versions, features,
effective profile recipe, commands, binary hashes and inputs. Do not contact
Glen from another thread, publish, dispatch CI or change private studies on the
strength of this handoff. Rust execution remains with Glen/CI.

The following criteria are evaluated by the [executed audit](#executed-audit-at-de9df8af)
below. Their scoped pass is pinned to that source, not inherited by later main.

| Gap | Executable evidence | Closure criterion |
| --- | --- | --- |
| Original `needs_fix` repair | Engine rehashed-nonminimum test; real-CLI adversarial test below | Valid non-minimal DNA has a verified-output/unverified-search receipt |
| Single-operation materialization | Engine persisted-metadata, direct/detached apply, stale/reopened-project and undo/redo tests; real-CLI container/node counts | One coherent operation/container/node; source preserved; no leaked stale writes |
| CDS-wide GC | Core integer/extrema/coupled/unrestricted-oracle tests; protocol; engine rehashed GC tests; real-CLI matrix | Exact counts and both strategies agree with an independent exhaustive oracle |
| Sliding-window GC | Core partial-codon and 36-word oracle; engine facts/tampering tests; real-CLI matrix | Every complete local window checked; global GC not replaced; no fabricated facts |
| Cancellation/budget/proofs | Core and engine cancellation-after-feasibility and exhaustion tests | No applicable cancelled DNA; incomplete search proves neither minimum nor infeasibility |
| Tutorial/adapters | Actual 06.07 CLI replay; shared-shell/direct-CLI/confirmed-MCP Rust tests | Same scientific payload and approval boundary, with explicit binary identity |
| Native/usefulness/runtime | Independent native review and bounded measurements below | Separate scoped verdicts; unavailable checks remain `not_run` |

## Executed Regression Follow-Up At `394dd5d`

On 2026-10-18 Glen executes the `8a0e2670` diagnostic-limit and `958b2238`
IUPAC follow-ups at `394dd5d0523ea91c19a8eb94d5ed8ad7bd5299b2`:
31 core, five protocol, 27 headless engine/MCP and four CLI unit tests pass.
The Python/oracle/actual-CLI suite passes 22/22 with zero skips; all four actual
CLI tests run. Retained-best/proof-limit, cancellation, lexical ties and the
180/10,800 IUPAC oracle cases pass, as do both strategies' joint GC, tampering,
materialization, save/reopen and same-session undo/redo checks. Locked check
and headless CLI/docs build pass on Debian x86_64/Rust-Cargo 1.97.1, dev/test
profiles with debug info and incremental compilation disabled.

The [verified evidence receipt](glen_audit_followups_20261008.md#focused-verification-at-394dd5d-2026-10-18)
records identity, hashes, logs and exclusions. Overall `needs_fix` is a separate
generated parity-document drift, not a failing sequence-design test. Its local
repair still needs Rust freshness/CI. This run supplies no new GUI, live-model,
runtime-performance, compiler-minimum, package or scientific verdict; the older
Linux/X11 and process measurements below remain separately attributed.

## Executed Audit At `de9df8af`

Glen executed the [issued handoff](glen_sequence_design_acceptance_handoff_20261014.md)
on 2026-10-14 at `de9df8af0a56923df01a7492f61005b1fe079d8c`, clean detached
tree `e3063e6b340e091e8d9c8079a7a5e295f82fada7`. Lockfile SHA-256:
`0b9729e60698e14d464b48230e606377a4e2e48054c2450be6f85760e999a6b9`.
His closing fetch retains that source. Verdict: sequence-design `prototype_pass`;
whole candidate/release `incomplete`; biological suitability `not_run`.

| Audited boundary | Executed result |
| --- | --- |
| Repair, single materialization, both GC contracts, cancellation/budget and adapters | 59 focused Rust passes: core 23, protocol 5, root library 27, CLI 4; locked check and parity check pass. |
| Retained Python/CLI replay | 20 passes, zero skips: sixteen source/oracle/harness checks plus four actual CLI tests. The matrix has 18 cases through both strategies, not twenty distinct CLI tests. |
| Native route | Linux/X11 under Xvfb/Openbox: Plan, reviewed Apply, created sequence window and same-session Undo/Redo. No model/provider call or Windows/macOS/Pi native result. |
| Runtime | 360 deterministic process/JSON runs: 300/1,200/3,000 bp, both strategies, five separate/joint/protected policies, four budgets, three repeats, 15-s timeout. Zero timeouts/nonzero exits; no pure-core timing or speedup claim. |

Platform: Debian forky/sid, amd64, Ryzen 5 3600, Rust/Cargo 1.96.1. The actual
GUI command builds the default dev profile; the headless CLI uses `package-opt1`
with no default features. These are distinct from installer/extracted-package
acceptance or the current `.13` `package-opt2` recipe. The parent stack was unset,
but frozen Cargo config restores 16 MiB; these are ordinary functional passes,
not stripped-child/default-stack evidence. An initial missing offline dependency
was fetched before the successful core run, not waived as a product pass.
Producer CLI SHA-256:
`bea6823292f88e9dc0bb51b76e2342a9d5b3c962f257ca8c63fad5e4bc4a8174`;
GUI SHA-256: `0ac4355e2459027bf3d0a06ff95f2a2fdb75d00fa83c6fa1a7c218118d85f2d2`.

Producer evidence root:
`/home/clawbio/.openclaw/workspace/artifacts/gentle-sequence-design-acceptance-de9df8af-20261014/`.
`REPORT.md`, `SUMMARY.json`, logs, observations and screenshots are available
locally at `/private/tmp/gentle-sequence-design-acceptance-de9df8af-20261014/`.
Original 1,600-entry `SHA256SUMS` SHA-256:
`0257840beccc62a6f9035ff066781696110bbc8b8c542fc7603bf5308e4c6df2`.
The separately delivered four-state archive has SHA-256
`a738a6f1220537852ceefe5e4137d6188c9a6a275c2e8587724fa9d9b705f245`;
its extracted root is `/private/tmp/gentle-de9df8af-four-saved-states-20261014/`
and its handoff manifest SHA-256 is
`ce53d7cbf91aeb524c6e1af11f149ade6c5760eba3b9dfa376190dc990aeccfa`.

Codex independently checks 1,576/1,600 original file hashes across those two
deliveries, including all four saved states against both original manifest and
CLI observations. Only 24 shader-cache files remain absent; this is not a claim
of complete package membership. All eight supplement hashes and its three
original-provenance copies match. Rechecking the 36 retained oracle reports and
approval digests passes offline. The four states preserve sentinel `ACGT`, exact
output/receipt binding, explicitly unverified search claims and one matching
operation/container/node per apply, including accepted valid non-minimal DNA.
Raw sequence/open/undo/redo captures are inspected, not recaptured or relabelled
as physical-platform acceptance. Original evidence stays outside the repository;
content hashes do not authenticate origin. Binary hashes are producer records,
not independently rehashed binaries in this delivery.

Glen independently re-read the pruning argument and executed the oracle matrix;
this is neither a general proof nor a claimed second-model/Claude source review.
Usefulness, biological/experimental suitability, live models and native platforms
beyond Linux/X11 remain `not_run`. The [separate terminal CI failure](audits/glen_followups_20261008/candidate_de9df8af/README.md)
blocks overall acceptance; nothing here verifies the later worker repair
`51a1634f` or newer main. The private core's [standalone compiler gate](decisions.md#workspace-compiler-compatibility)
remains open; Rust 1.96.1 and moving `stable` CI do not select a minimum.

## CI Regression Evidence (Not Independent Acceptance)

Read-only cross-check on 2026-10-13: the
[Windows job 115417977228](https://github.com/smoe/gentle_rs/actions/runs/38432492850/job/115417977228)
completed successfully at `8b6f8186c81636491dadb446e94ca68429ba0860`.
Its [decoded log source](https://api.github.com/repos/smoe/gentle_rs/actions/jobs/115417977228/logs)
records that exact checkout and binary source stamp, Rust 1.99.0
(`b940084d7`, `x86_64-pc-windows-msvc`), successful locked full-workspace tests,
and successful LF/CRLF parity-matrix checks. No local Rust command was run.

This establishes existing Windows CI regression evidence, not Glen's new audit,
an independently assessed usefulness/performance verdict, or native GUI review.
The quiet workspace output does not enumerate individual test names; do not
invent named-test counts from aggregate summaries. These audit additions were
not in that SHA. Retain the source-bound CI evidence separately and execute the
new targeted/real-CLI procedure at the integrated candidate before closing the
table above. Later SHAs do not automatically inherit this job's acceptance.

### Integrated Windows Regression Evidence

Read-only cross-check on 2026-10-13: the
[Windows job 115490639067](https://github.com/smoe/gentle_rs/actions/runs/38453643146/job/115490639067)
completed successfully at `024ed6c0520d38c4221517219d44c792084de4e5`, which
includes the repair/GC audit additions in `5e11be7c`. Its
[decoded log source](https://api.github.com/repos/smoe/gentle_rs/actions/jobs/115490639067/logs)
records the matching CLI source stamp, Rust 1.99.0 (`b940084d7`,
`x86_64-pc-windows-msvc`) and successful `cargo test -q --locked --workspace`.
The root-library summary has 4,095 passes, zero failures and 14 ignored tests;
these are aggregate counts, not a named sequence-design acceptance result.

The same run's
[safety job 115490639076](https://github.com/smoe/gentle_rs/actions/runs/38453643146/job/115490639076)
ran the companion/oracle suite: 16 checks passed and all four actual CLI tests
were explicitly skipped because `GENTLE_TUTORIAL_BIN_DIR` was not configured.
Those skipped checks cover direct/shared-shell apply, the exhaustive two-strategy
matrix, rehashed attacks and exhausted-search behavior. Ordinary workspace
regressions do not substitute for their retained replay or an independent
verdict. At inspection macOS was still running; no whole-run, native GUI,
usefulness, runtime or release acceptance was claimed at that inspection. No
local Rust was run; the table remained pending then. The separately pinned audit
above records the later scoped closure rather than changing this CI's verdict.

## Exact-Candidate Functional Checks

Run focused checks before broader gates; retain full logs and exit codes.
Earlier motif-only evidence and recorded macOS passes do not certify this SHA.
Dependency fetching failures are environmental, not proof of a product defect.

```sh
cargo test -p gentle-sequence-design --locked --offline -j 1
cargo test -p gentle-protocol --locked --offline -j 1 sequence_design
cargo test --lib --no-default-features --locked --offline -j 1 sequence_design -- --test-threads=1
cargo test --bin gentle_cli --no-default-features --locked --offline -j 1 sequence_design -- --test-threads=1
cargo check -q --locked --offline -j 1
cargo run --bin gentle_examples_docs --no-default-features --locked --offline -j 1 -- parity-matrix-check
cargo fmt --all --check
git diff --check
```

The library filter includes confirmed MCP and in-session undo/redo. Reopening
a CLI process does not persist an undo stack; the companion correctly refuses
invented cross-process history. Route parity is not model quality or script-
runtime acceptance. No inner model is needed for these checks.

Build the headless CLI once, recording the build separately from runtime:

```sh
cargo build --profile package-opt1 --no-default-features --locked --offline --bin gentle_cli -j 1
```

Review the synthetic inputs and explicitly authorize their test applies. Then
in that clean checkout supply the exact approved SHA and a fresh destination.
These tests run no compiler, model or network operation:

```sh
export GENTLE_TUTORIAL_BIN_DIR=/absolute/path/to/target/package-opt1
export GENTLE_SEQUENCE_DESIGN_CANDIDATE=THE_AUTHORIZED_FULL_40_CHARACTER_SHA
export GENTLE_SEQUENCE_DESIGN_EVIDENCE_DIR=/absolute/path/to/fresh/sequence-design-audit
python3 -m unittest scripts.test_sequence_design_audit_oracle scripts.test_sequence_design_tutorial -v
```

Retained mode fails rather than skipping when its binary path is missing,
and refuses a missing/mismatched candidate stamp or existing test
destination. It hashes the CLI and retains per-command argv, exit/timeout,
stdout/stderr, raw requests/previews, disposable states and artifact hashes.
HOME/XDG and inherited GENtle settings are isolated. Keep the unittest console
log too: command exit zero alone does not prove assertions passed.
`observations.json` records observations, never an independent verdict. Stamps
and hashes are identity evidence, not authenticated origin; Glen independently
verifies build/source identity, features, effective profile and toolchain.

There are 18 fixed 4/36-word cases, each run through both strategies: motif-only,
global-only, window-only, jointly coupled, protected partial codons, flanks,
reverse IUPAC matches, empty integer windows and proven impossibility. The oracle
uses literal standard-code synonyms and cross-multiplied integer inequalities;
it imports no engine matcher or translator. Completed results must agree on
feasibility, minimum-edit DNA and lexicographic tie-break, not visited-node count,
algorithm ID or approval digest. Validity/proof status is checked separately.

The adversarial encoder mirrors declared struct wire order only. First compare
its digest with an authentic preview; mismatch is a harness/contract failure,
not an excuse to omit honestly rehashed attacks. Invalid DNA, frozen-stop edits,
fabricated window facts and stripped bound policies must fail without writes.
Valid two-edit MEF DNA with fabricated minimum claims may apply, but the receipt
must still say `search_claims_verified=false`; refusing all non-minimal DNA would
change the intended contract rather than repair it.

Without `GENTLE_TUTORIAL_BIN_DIR`, four actual CLI tests explicitly skip. The
sixteen source/oracle/harness checks can pass offline without closing any Rust/CLI/native
acceptance. No skip may count as a green acceptance gate.

## Static Conflict-Search Cross-Read

Codex source inspection at `8b6f8186` found no additional defect. This is not a
Claude review, independent execution, or proof that the implementation is bug-free.
The unchanged inspected Git blobs are:

- `conflict_search.rs`: `3951cd838fab33d7276b860203b2c37218d91072`.
- `gc.rs`: `5563116761850d5aabc1611fc8ff2ab7ff7ef437`.
- `gc_window.rs`: `a2c4da7c9df569ae5225b317bf8b1dc1a07e9d08`.
- Engine `sequence_design.rs`: `8258f6996243aa8700c3845cdb90a5ed05f23de6`.

Load-bearing arguments to independently cross-read:

1. Unassigned codons still equal the original. A feasible target consistent with
   assigned codons must change an unassigned codon in any current motif conflict.
   For GC, its net change inside the conflicting interval requires at least one
   direction-capable codon. Branching includes its target synonym, including
   neutral/worsening alternatives; it does not require greedy repair.
2. Partial windows use only overlapping codon bases for capability. Whole-codon
   GC is not a sound substitute. All candidates receive the full predicate.
3. Assignments fix disjoint codons at nonzero nucleotide cost. Pruning an
   infeasible node already as costly as a feasible best cannot discard a better
   result or equal-cost lexicographic winner: every descendant costs more.
   Feasible nodes stop safely because additional edits cannot improve their
   primary objective. Tie-break tests remain essential.
4. `complete` requires an exhausted/pruned tree or unchanged feasible input;
   budget stops do not establish it. Cancellation clears output and proof flags,
   including after finding a best. Extrema are negative proofs only.
5. Apply revalidates output and bound facts, not search history. Its unverified-
   preview wrapper and default-false legacy receipt are deliberate. The outer
   operation owns materialization/undo; no nested fork remains.

## Independent Native And Runtime Verdicts

Replay 06.07 in native GUI Shell at the same candidate; review preview, apply
explicitly, open created DNA and verify same-session undo/redo. Retain raw
screenshots only under approved capture policy. Headless tests establish no
Windows/macOS/Pi, screenshot or live-model result.

Use the first handoff's bounded 300/1,200/3,000-bp synthetic ladder and budgets
1/100/1,000/100,000, now with both strategies and separate/joint GC policies.
Document generation, source hashes, three repeats, per-process timeout and
effective budgets. Separate build, process/JSON and core-search costs; retain
wall/CPU/RSS samples. Timeout is censored evidence, not speedup or infeasibility.
Do not add heuristics, thresholds, dependencies or runtime claims from toy counts.

The main-maintaining thread retains Glen's exact-SHA report/manifests and updates
the plan/roadmap only for demonstrated closures. Return `needs_fix` for a product
defect, `incomplete` for missing evidence/environment limits, or scoped
`prototype_pass` after the required audit. Usefulness, runtime, native and
biological suitability stay independent. Do not re-label the historical audit.
