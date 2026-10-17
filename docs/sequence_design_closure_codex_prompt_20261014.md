# Codex prompt: record sequence-design closure and the conditional MSRV defect

Prepared 2026-10-14 against `main` at
`43fb14928077eca57b27287215bb53a90fd7ece3`. Documentation and policy only. No
solver, adapter, approval, fixture or release-gate change is authorized here.

Two linked tasks in **one commit**: record Glen's completed audit, and record
the Rust-minimum defect he confirmed, without reversing the workspace
compiler decision that deliberately removed the previous floor.

## Constraints

- No local Rust build, test, benchmark or dependency update. Rust execution
  stays with Glen/CI.
- No push, tag, CI dispatch, publication or new workflow run.
- `docs/roadmap.md` must stay exactly 400 lines; the release-candidate checks
  assert it. Rewrite the existing parking-lot line in place rather than adding
  lines.
- Per `AGENTS.md`, do not add a `Done` entry to the roadmap. Completed work
  goes to `docs/CHANGELOG.md` in this same change.
- Leave `paper`, `output/`, `Cargo.lock` and historical evidence untouched.
- Do not relabel the historical
  [first-slice audit](glen_sequence_design_handoff_20261005.md). Its
  `needs_fix` at `85d23ddb` stays historical evidence.

## Evidence to record

Glen's independent audit, returned 2026-10-14, at the authorized candidate
`de9df8af0a56923df01a7492f61005b1fe079d8c` (confirmed as `origin/main`, clean
isolated worktree, verified by a closing fetch):

| Result | Value |
| --- | --- |
| Targeted Rust tests | 59 / 59 passed |
| Real revision-bound CLI acceptance | 20 / 20 passed, **no skips** |
| Two-strategy oracle agreement | Both strategies agree with the independent exhaustive oracle |
| Bounded runtime | 360 / 360 runs completed, 0 timeouts, 0 errors |
| Native Linux/X11 route | Plan → reviewed apply → sequence window → undo → redo passed |
| Sequence-design verdict | scoped `prototype_pass` |
| Overall / release verdict | `incomplete` |
| Biological suitability | not assessed, no claim |

Evidence package, held **outside this repository** on Glen's host:
`/home/clawbio/.openclaw/workspace/artifacts/gentle-sequence-design-acceptance-de9df8af-20261014/`,
principal file `REPORT.md`, manifest `SHA256SUMS`, manifest SHA-256
`0257840beccc62a6f9035ff066781696110bbc8b8c542fc7603bf5308e4c6df2`,
1,600 inventoried files, 94 MiB, including logs, oracle/CLI receipts, 360
runtime runs and native screenshots. His run-local 15 GiB worktree/cache was
removed afterwards. No push, CI dispatch, model call or publication occurred.

Record identity and location only. Do not copy the package in, restate its
counts as repository test results, or imply the manifest hash authenticates
origin rather than content.

The `incomplete` overall verdict is **not** a sequence-design fault. Glen
attributes it to the general CI at the same SHA failing in a separate
stack-runner fixture step on Linux and Windows. That failure is already
retained at
[`candidate_de9df8af`](audits/glen_followups_20261008/candidate_de9df8af/README.md),
with the narrow repair in `51a1634f`. Cross-reference it; do not reopen,
re-diagnose or merge the two records, and do not claim `51a1634f` is accepted —
that run predates it.

## Task 1: the conditional Rust-minimum defect

Glen's finding: *"Fehlende MSRV-Deklaration bleibt ein Vertragsdefekt vor
Extraktion/Veröffentlichung"* — the absent minimum is a contract defect
**before extraction or publication**, not a defect in current behavior.

This does not conflict with
[Workspace Compiler Compatibility](decisions.md#workspace-compiler-compatibility),
and you must not reverse it. That decision removed the sequence-design crate's
Rust 1.85 floor precisely because it was unverified, states that no explicit
minimum is selected now, and requires any future member minimum to inherit
`workspace.package`. `scripts/test_release_candidate.py` enforces this: a
member `rust-version` is accepted only as `{ workspace = true }` with a
workspace-level value present.

So do **not** add `rust-version` anywhere, and do not pick a number. Instead:

1. Amend the Workspace Compiler Compatibility section to state the consequence
   explicitly: with no declared minimum, a passing build certifies only the
   toolchain actually used, and extracting or publishing any workspace member
   as a standalone crate is blocked until a verified workspace minimum exists.
2. Record what would license declaring one: a CI job pinning a specific older
   toolchain that builds and tests the dependency-free
   `gentle-sequence-design` crate alone, with the passing toolchain becoming
   the declared `workspace.package` value. Declaring a floor that no job
   exercises would recreate exactly the unverified claim `2bc0d110` removed.
3. Note in the sequence-design plan and
   [Synthetic Sequence-Design Boundaries](decisions.md#synthetic-sequence-design-boundaries)
   that this gate applies to the independently-versioned-core premise: the
   crate stays private and in-workspace until then.
4. Record the only toolchain evidence that currently exists: Glen's audit ran
   on his recorded toolchain at `de9df8af`, and CI runs `dtolnay/rust-toolchain@stable`.
   Do not infer compatibility with any older release.

Adding the verifying CI job is a separate, later change. Do not add it here.

## Task 2: record the closure

**`docs/CHANGELOG.md`** — new dated entry at the top. State the audited SHA,
the table above, the external evidence identity, the scoped `prototype_pass`,
and that the `incomplete` overall verdict comes from the unrelated stack
fixture. State plainly that this is not `.12` release acceptance, not
biological validation, and not a performance claim; the 360 completed runtime
runs are bounded measurements, not a speed result.

**`docs/dna_sequence_optimization_plan.md`** — remaining-stages item 1 is now
closed by execution. Rewrite it as closed against the exact SHA and evidence
manifest, keeping items 2 and 3 scoped as they are. Item 2's GC contracts are
covered by this audit's oracle matrix, so fold what it actually closed into
item 1's record rather than leaving item 2 implying the GC contracts were never
executed; retain anything it asks for that Glen did not do, explicitly as
`not_run`. Also record that the conflict-search slice has now had independent
execution but still no second-model source review — the oracle agreement
addresses the pruning-soundness argument empirically, not by review.

**`docs/roadmap.md`** — rewrite parking-lot line 303 in place, one line. It
must stop requesting the audit Glen already completed, and instead carry only
what remains: the conditional MSRV gate before extraction/publication, the
unexecuted native checks beyond Linux/X11 (no Windows, macOS or Pi native
result), absent usefulness and biological verdicts, and that the design stays
outside `.12`. Keep the links to the plan, the acceptance procedure and the
retained first-slice audit.

**`docs/decisions.md`** — update the Synthetic Sequence-Design Boundaries
status line so it reflects independently executed motif and both GC contracts
at a pinned SHA, with the contract itself unchanged. Do not weaken any
invariant, and do not add a biological or performance claim.

## Also in this commit

`docs/glen_sequence_design_acceptance_handoff_20261014.md` is currently
untracked and carries a correction notice: it was first issued with a
fabricated 40-character tail for the `de9df8af` prefix, and the verified hash
is `de9df8af0a56923df01a7492f61005b1fe079d8c`. Glen resolved the prefix and
audited the correct commit, so no evidence is affected. Commit the file with
that notice intact as the issued-prompt record; do not silently drop it.

## Verification expected of you

Offline checks only: TOML parsing, the release-candidate and compiler-policy
tests, whitespace, release-link and roadmap line-count checks, plus
`scripts/maintenance_chore.py session-close`. Report what ran, what passed and
any warning. State explicitly that no Rust command was run and that Rust/native
verification remains with Glen and CI.
