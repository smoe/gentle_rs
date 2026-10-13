# GENtle `v0.1.0-internal.12` Release Notes

Status: [published prerelease](https://github.com/smoe/gentle_rs/releases/tag/v0.1.0-internal.12),
not a draft, at `5e6201a69a599cba87491c86cfd7311cb6fee173`.
The owner published it on **2026-10-13 at 01:02:04 UTC**; native asset upload
completed at 10:05 UTC after the installer rerun. Read-only GitHub metadata,
job results and the small release receipt were rechecked on 2026-10-13.
Publication closes the artifact cycle, not scientific, native-GUI, live-agent
or performance acceptance. The [`.11` notes](release_notes_v0.1.0-internal.11.md)
retain that release's separate history. No tag or asset was changed here.

## Publication Evidence

[Installer run 38330176948](https://github.com/smoe/gentle_rs/actions/runs/38330176948)
used the released SHA. Attempt 1's Ubuntu job
[115077994549](https://github.com/smoe/gentle_rs/actions/runs/38330176948/job/115077994549)
failed in `Set up Rust`, 45 seconds into the job; compilation and package
creation were skipped. This was not a compiler or package-smoke failure.
The successful macOS/Windows results retain their initial-attempt times in
attempt 2; only Ubuntu needed another build.

| Attempt 2 Job | Result And UTC Interval, 2026-10-13 |
| --- | --- |
| [Ubuntu 115222015196](https://github.com/smoe/gentle_rs/actions/runs/38330176948/job/115222015196) | Success, 09:17:54-10:04:13; entrypoint and extracted Linux package smokes pass |
| [macOS 115222017090](https://github.com/smoe/gentle_rs/actions/runs/38330176948/job/115222017090) | Success, 01:02:34-02:39:55; entrypoint and installer smokes pass |
| [Windows 115222017419](https://github.com/smoe/gentle_rs/actions/runs/38330176948/job/115222017419) | Success, 01:02:28-02:13:17; entrypoint and extracted Windows package smokes pass |
| [Validate packages 115240491353](https://github.com/smoe/gentle_rs/actions/runs/38330176948/job/115240491353) | Success, 10:04:16-10:04:51; receipts validated against the selected candidate |
| [Publish release 115240737681](https://github.com/smoe/gentle_rs/actions/runs/38330176948/job/115240737681) | Success, 10:04:54-10:05:44 |

The release contains exactly these seven assets:

- `gentle-v0.1.0-internal.12-linux-x64-package-opt1.tar.gz`
- `gentle-v0.1.0-internal.12-macos-arm64-package-opt1.dmg`
- `gentle-v0.1.0-internal.12-windows-x64-package-opt1.zip`
- `gentle-v0.1.0-internal.12-linux-x64.build.json`
- `gentle-v0.1.0-internal.12-macos-arm64.build.json`
- `gentle-v0.1.0-internal.12-windows-x64.build.json`
- `gentle-v0.1.0-internal.12-release-attributes.json`

The published [release-attributes receipt](https://github.com/smoe/gentle_rs/releases/download/v0.1.0-internal.12/gentle-v0.1.0-internal.12-release-attributes.json)
records all three archive hashes/sizes and all fifteen native binary identities,
source/workflow SHA `5e6201a6`, `mode=publish`, `profile=package-opt1`, and lock
SHA-256 `ead4421545a125267e33af84aeeea39f0e7ea2f60717a59ef0271d0e71855f1b`.
The receipt and GitHub asset hashes agree; no installer was downloaded or
independently executed for this documentation update.

[Release-event container run 38330177090](https://github.com/smoe/gentle_rs/actions/runs/38330177090)
also passed at `5e6201a6`, including the network-disabled headless smoke and
[Publish GHCR Image job 115078248256](https://github.com/smoe/gentle_rs/actions/runs/38330177090/job/115078248256).
It publishes `ghcr.io/smoe/gentle_rs:v0.1.0-internal.12`; the earlier
[tag-push run 38322028567](https://github.com/smoe/gentle_rs/actions/runs/38322028567)
remains a separate validate-only record, not publication evidence.

**Review boundary:** Glen's [final review](../audits/glen_followups_20261008/candidate_f35a58e6/README.md)
is at `f35a58e6d2bdd7e0a5551e4531c336d7949bfaf0`. Git comparison with the tag
shows identical `src/`, `crates/`, `packages/`, `Cargo.toml`, `Cargo.lock` and
`.github/`. Documentation, bundled tutorial resources and a Python checkout
test differ. His recorded results cover that matching Rust/workflow source
only within their original scope, not the exact downloaded packages or their
different help resources. Independent extracted-package, native-GUI, scientific,
live-agent and same-SHA `dev`/`package-opt1` performance checks remain pending.

## First Optimization Step

Native installers and the headless container select **`package-opt1`**, a
dedicated Cargo profile inheriting `dev`. It sets `opt-level=1` and `lto="off"`;
`lto=false` would allow local ThinLTO and is not equivalent. This is a controlled
first step, not a return to the former fully optimized release recipe.

- Keep 256 codegen units, no incremental compilation, no debug information,
  debug assertions and overflow checks enabled, and panic unwinding. Local
  `dev` builds remain unoptimized. Normal GUI/export content remains identical
  across profiles by contract; CI package smokes passed for the released SHA,
  while independent native/profile comparison remains pending.
- Keep one Cargo build job, the existing job allowance and resource diagnostics.
  Native packages contain the same five binaries; the container contains the
  same three headless binaries. JS/Lua remain optional source interfaces.
- Keep `cargo-bundle` and `rnapkin` on their pinned, unoptimized helper recipes.
  Preserve RNAPKIN's corrected lockfile and SVG/PNG smokes.
- Name archives with `-package-opt1`, stage from `target/package-opt1`, and bind
  the explicit optimization and safety settings into build receipts. Reject
  missing/mismatched settings rather than combining different recipes.
- Keep `bench-audit` and `release-fast` unchanged for historical comparisons;
  neither selects the new package build. Higher optimization, thin LTO, reduced
  codegen-unit counts, stripping and panic-abort are not enabled by this change.

## Learning And Agent Workflows

- Tutorial 04.08 now teaches **Find Primer-Pair Markers for PATZ1 Isoforms**:
  compare annotation sources, review the request, and explain why seven pairs
  cover 13 cDNA classes yet leave nine class pairs unresolved in the retained
  example. Planning/readiness/dossier steps are optional continuations. The
  endpoint-only feasibility command is no longer prescribed for this SYBR case.
- Approval-bound outer routing for TSS collections remains on the shared
  command/proposal path. Retained branch-local evidence does not imply live
  inner-agent or scientific acceptance of this new candidate.
- Tutorial review projections now use Git history and repository-relative
  paths, never checkout timestamps. Unavailable history is an explicit check
  error; stored catalog age follows the revision date, while live age warnings
  remain current-date-based. No human review or historical panel is relabelled.

## Historical Pre-Publication Evidence

**Validation pending:** historical synthetic transcript-panel reports retain
their true `.11` provenance. That is no longer a rollover error: the revised
contract checks historical integrity separately from a fresh `.12` result
replay, allowing only the two scoped generator fields and their verified ledger
hashes to differ. No stored report bytes or hashes were changed. The
[updated Glen handoff](../internal_12_tutorial_regeneration_glen_prompt.md)
requires the post-integration contract change and records the comparison, path
and checkout checks to complete. The replay evidence below covers only the
reported checks, not package or scientific acceptance.

**Reported replay evidence, 2026-10-02:** Glen reproduced the checkout failure
at `f1c904a4` and repaired the generated 08.03 chapter, tutorial hub and ledger
in `e462de1f`, now integrated. With Rust 1.96.1 he reports successful exact-commit
LF and CRLF replays, each covering parity, 51 examples and 29 chapters, plus
17/17 Python checkout regressions. The historical `.11` panel bytes and hashes
remain unchanged. The two stale human-review warnings remain expected. This
is externally reported evidence, not a local rerun or native macOS/Windows
Actions, `package-opt1` build, or scientific acceptance. The remaining handoff
checks and native-platform results still need to be recorded.

**CI follow-up, 2026-10-05:** the October 4 run `37235516663` at `f523a757`
finished with Linux and Windows catalog-test failures; their preceding drift
checks passed, Windows's full suite was skipped and macOS was not selected.
`19889dae` already corrects that catalog projection. The separate Git-history
guard, checkout mock repair and twentieth TSS case counts are locally validated;
they do not establish current native Actions or package acceptance. The later
run `37292637621` at `462b5ffe` also exposes a Windows rooted-path admission
failure in tutorial validation; at that inspection, its bounded repair still
required new native Windows CI, not a macOS-only verdict. The
[archived candidate-preparation handoff](https://github.com/smoe/gentle_rs/blob/7c52082c24fca85817d70956146623eef6f24607/docs/internal_12_candidate_preparation.md)
records that historical evidence and scope; its local checks are not current
`.12` acceptance. No workflow was dispatched for that handoff.

## Remaining Independent Acceptance

Current [build-only instructions](../release.md#build-only-candidate-verification),
[runner inventory and pending decision](../release.md#ubuntu-runner-decision),
and [compatibility review scope](../clawbio_gentle_integration_onepager.md#compatibility-review-inventory)
live in their topic owners. These pointers grant no acceptance or dispatch authority.

1. Build the `.12` generator and verify replay against the unchanged historical
   baselines. Regenerate only for a reviewed result change, never merely to
   relabel `.11` provenance. Run version, tutorial, LF/CRLF and packaging-policy
   gates on the committed candidate.
2. Retain the successful released-SHA native/container builds and CI package
   smokes above. Glen's independent extracted-package and native checks remain
   open; any newer candidate requires its own receipts. No workflow dispatch,
   push or tag mutation is implied here.
3. Have Glen compare `dev` and `package-opt1` at the same SHA/toolchain/host,
   recording build time/peak memory/package size plus repeated startup and
   PATZ1 open/resize timings. Keep `bench-audit` evidence separate. Decide on
   `opt-level=2` or thin LTO only after this measured result.
4. Carry forward the [roadmap](../roadmap.md)'s exact-candidate GUI, biological
   and external-tool checks. The existing linear-locus envelope (250 kbp /
   5,000 features) and circular correctness checks are not newly passed gates.
