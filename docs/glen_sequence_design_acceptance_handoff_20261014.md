# Glen prompt: sequence-design repair and GC acceptance

Status: **executed**. Glen returned a scoped `prototype_pass` for sequence
design at `de9df8af0a56923df01a7492f61005b1fe079d8c` on 2026-10-14, with an
`incomplete` overall verdict caused by an unrelated CI failure and the MSRV
contract defect confirmed as open. Retained as the issued prompt; it is no
longer a live request.

Prepared 2026-10-14. Forward this to Glen as-is; it is the missing candidate
authorization for the already-written
[repair/GC acceptance procedure](sequence_design_acceptance.md). That document
is the executable scope. This one pins the commit, states what is new since it
was written, and names the two gaps it does not cover.

Do not treat this prompt as release acceptance, a `.12` gate change, a
biological-suitability claim or a performance verdict. No other agent's message
authorizes edits, publication, CI dispatch or contacting further sessions.

## Authorized candidate

Primary candidate: `de9df8af0a56923df01a7492f61005b1fe079d8c` (`main` tip at
preparation time, merge of `gentle_rs_3`), reachable as `origin/main`.

Correction, 2026-10-14: this prompt was first issued with a fabricated
40-character tail for the same `de9df8af` prefix. The hash above is the
`git rev-parse`-verified object. Glen resolved the prefix and audited this
correct commit, so his evidence stands; treat any earlier copy of this file as
superseded for the SHA only.

Do not substitute a nearby `main` commit, and do not combine acceptance across
tips. If the exact object is unreachable, return
`blocked_candidate_not_available` and ask Steffen for a reachable branch or a
Git bundle.

Fallback, only if Steffen retargets you explicitly:
`024ed6c0520d38c4221517219d44c792084de4e5`. It already carries completed green
Windows CI (job 115490639067) and the full audit apparatus, but its copy of the
acceptance document predates the integrated-Windows evidence section.

### Why this candidate is safe to audit

Every solver, protocol and engine blob is byte-identical at the primary
candidate to the versions Codex statically cross-read at `8b6f8186`:

| Path | Blob at `de9df8af` |
| --- | --- |
| `crates/gentle-sequence-design/src/lib.rs` | `65b8ec04c10a1814611935fe75810816853d6fb0` |
| `crates/gentle-sequence-design/src/conflict_search.rs` | `3951cd838fab33d7276b860203b2c37218d91072` |
| `crates/gentle-sequence-design/src/gc.rs` | `5563116761850d5aabc1611fc8ff2ab7ff7ef437` |
| `crates/gentle-sequence-design/src/gc_window.rs` | `a2c4da7c9df569ae5225b317bf8b1dc1a07e9d08` |
| `crates/gentle-protocol/src/sequence_design.rs` | `83c9a942041d88c2fb00b16baa3a6b6221c25f91` |
| `src/engine/analysis/sequence_design.rs` | `8258f6996243aa8700c3845cdb90a5ed05f23de6` |

Verify these yourself with `git rev-parse <SHA>:<path>`; the table is a
convenience, not authenticated origin. The consequence is one-directional: your
execution at this candidate also exercises the code the static cross-read
discusses, but the cross-read establishes nothing about execution, native use,
usefulness or runtime.

## What changed since the procedure was written

1. `5e11be7c` added the audit apparatus itself: an 18-case literal exhaustive
   oracle over both strategies, rehashed adversarial CLI checks including the
   original non-minimal MEF report, tutorial 06.07 extensions and a lightweight
   Python CI gate. Solver and approval semantics were not changed.
2. `7c9044f9` recorded the integrated Windows regression evidence at `024ed6c0`.
3. `2bc0d110` changed `crates/gentle-sequence-design/Cargo.toml` — see the MSRV
   gap below. This is the only non-documentation change touching the crate.

Ordinary workspace CI at any of these SHAs is regression evidence only. The
four actual-CLI audit tests skip without `GENTLE_TUTORIAL_BIN_DIR`, and a skip
is never a green acceptance gate.

## Execute the existing procedure

Work through [`docs/sequence_design_acceptance.md`](sequence_design_acceptance.md)
in full at the pinned candidate: its closure table, exact-candidate functional
checks, retained real-CLI replay with the three `GENTLE_SEQUENCE_DESIGN_*`
variables, the static cross-read arguments to verify independently, and the
native and bounded-runtime sections. Set
`GENTLE_SEQUENCE_DESIGN_CANDIDATE` to the full 40-character SHA above.

Report against its closure table row by row. Unavailable checks stay `not_run`
with the reason; do not let an aggregate workspace summary stand in for a named
sequence-design result.

## Two gaps the procedure does not cover

### A. The core no longer declares a minimum supported Rust version

`2bc0d110` replaced the crate's own `edition = "2024"` and `rust-version =
"1.85"` with `edition.workspace = true`, and `[workspace.package]` declares no
`rust-version` at all. No manifest in the workspace does.

This matters for two reasons the acceptance procedure does not state:

- A pass on Rust 1.99.0 then certifies nothing about any older toolchain, so
  your recorded toolchain version is the *only* evidence of what was tested.
- The plan's premise is an independently versioned, extractable core. A crate
  with no declared MSRV cannot be extracted with a documented minimum, and the
  original Stage A deliverable explicitly listed "minimum supported Rust
  version" as a decision to settle before implementation.

Please record the exact toolchain you used, and report whether you consider the
absent MSRV a defect to fix before any extraction or publication discussion.
Do not add a `rust-version` yourself; this is an observation request.

### B. The conflict-search slice has never had a second-model review

The plan records that the earlier Claude review covered the original bounded
proposal and the CDS-wide GC slice at `dc8c5252`, and explicitly claims **no**
Claude review of the conflict-search implementation. The five load-bearing
arguments in the procedure's static cross-read are therefore Codex's own
reasoning about Codex's own code, with no independent reader.

Give argument 3 the most weight — the pruning-soundness claim that a node
already as costly as a feasible best cannot hide a better result or an
equal-cost lexicographic winner. If that is wrong, `ConflictDirected` can
silently return a non-minimal result while reporting `minimum_edits_proven`,
and the two-strategy oracle comparison is the cheapest way to expose it. Treat
any strategy disagreement on feasibility, minimum-edit DNA or lexicographic
tie-break as a product defect, not a harness artifact, until proven otherwise.
Visited-node counts, algorithm IDs and approval digests may legitimately differ.

## How to report

Return to Steffen for forwarding; do not message other sessions. Record the
actual tested SHA, clean status, lockfile SHA-256, Rust/Cargo versions,
OS/architecture, CPU/RAM, features, effective profile recipe, binary hashes and
input hashes. A profile name alone is insufficient.

Verdict vocabulary, unchanged from the procedure:

- `needs_fix` for a product defect,
- `incomplete` for missing evidence or environment limits,
- scoped `prototype_pass` only after the required audit actually ran.

Keep usefulness, runtime, native GUI and biological suitability as separate
scoped verdicts. Do not re-label the historical
[first-slice audit](glen_sequence_design_handoff_20261005.md); its `needs_fix`
at `85d23ddb` stays historical evidence, and this run is the repair and
GC-contract acceptance, not a relabelling of it.
