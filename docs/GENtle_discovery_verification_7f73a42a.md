# GENtle discovery/decision repair verification

Date: 2026-11-02 (Europe/Berlin)

## Verdict

The discovery/registry repair at `7f73a42ab97d0aa3c11c37182da9e7b7c4b60c6a`
passes its focused executable gates: canonical glossary descriptors, unique registry
identity, preserved per-interface JS/Lua usage, deterministic same-interface conflict
refusal, optional readiness parsing, unique parity rows/counts, and canonical generated
matrix freshness.

The broader requested `glossary_` library filter is **not green**. It has one primary,
isolated fixture failure in
`engine_shell::tests::glossary_cli_usage_flags_parse_one_by_one`: the smoke command
supplies `--binding '{}'`, while the shared parser correctly requires the typed observed
SAM binding schema. The test source blob and the affected glossary row are unchanged
from the candidate parent `9488f52d7242b910bffba7cc7c9692c37bebd7fa`, so this is
not introduced by the registry repair, but it remains a real gate failure.

GitHub Actions run 41373114753 is consequently not an overall pass. Release-safety and
Linux/headless succeeded; Windows and macOS failed at `Glossary command regressions`.
The unauthenticated API exposes the failed step but denied raw log download (HTTP 403),
so the platform failures are not independently attributed to the same assertion here.

## Identity and environment

- Tested commit: `7f73a42ab97d0aa3c11c37182da9e7b7c4b60c6a`
- Tree: `a338b931b0985a0b8cf25321f9964497914489c7`
- Parent: `9488f52d7242b910bffba7cc7c9692c37bebd7fa`
- Final fetched `origin/main`: identical to tested commit
- Lockfile SHA-256: `a393a1f3803b374bbaf5ef4f3e23813e974b1676f61a669dacb24063c059fb6a`
- Worktree: detached and clean; active dirty checkout untouched
- Target: run-owned `/home/clawbio/work/gentle-verify-7f73a42a-target`
- Platform: Debian forky, Linux x86_64, AMD Ryzen 5 3600, 12 logical CPUs, 62 GiB RAM
- Rust/Cargo: 1.97.1; ordinary dev/test profiles; repository `RUST_MIN_STACK=16777216`
- Features: first protocol test is package-local; remaining Rust gates use
  `--no-default-features`; test execution is serial (`-j1`, `--test-threads=1`)

## Executed checks

1. `cargo test --locked --offline -p gentle-protocol --lib glossary_canonicalization_tests -- --test-threads=1`
   - PASS: 4 passed, 0 failed, 0 ignored, 116 filtered.

2. `cargo test --locked --offline --no-default-features -j1 --test capability_registry_parity --test parity_matrix_freshness -- --test-threads=1`
   - PASS: 23/23 registry parity and 1/1 matrix freshness; no failures, ignores, or filtered tests.

3. `cargo test --locked --offline --no-default-features -j1 --lib glossary_ -- --test-threads=1`
   - FAIL: 12 passed, 1 failed, 0 ignored, 3396 filtered.
   - Primary failure: `engine_shell::tests::glossary_cli_usage_flags_parse_one_by_one`.
   - Isolated exact rerun: 0 passed, 1 failed, 3408 filtered; same typed-binding error.

4. `cargo run --locked --offline --no-default-features --bin gentle_examples_docs -- parity-matrix-check`
   - PASS: `status: ok`; committed `docs/gui_cli_mcp_parity.md` is canonical.

5. Focused source-metadata checks:
   `python3 -m unittest scripts.test_capability_maturity.DiscoveryCatalogTests scripts.test_tutorial_checkouts.DocumentationStatusTests.test_every_decision_section_has_one_explicit_status -v`
   - PASS: 3 passed, 0 failed, 0 skipped.

## Exact projections

- Matrix capability rows: 952 raw / 952 unique / 0 duplicates.
- Open gap keys: 164 raw / 164 unique / 0 duplicates.
- Unique gaps: MCP 42, JS 61, Lua 61.
- Summary: MCP 205 prominent / 545 shell-only / 42 gap; JS 39 / 642 / 61;
  Lua 38 / 643 / 61.
- Glossary: 614 rows / 571 unique paths / 43 duplicate raw paths. Duplicate raw paths
  are canonicalized to one registry descriptor per `(source, name/path)`.
- Decisions: 72 sections, each with exactly one classified status: 70 active,
  2 deferred, 0 superseded.

## Binary hashes

- `gentle_examples_docs`: `39e175ac660b483931bc3700a67a0107386cdcbe9eff61ba033f5fd0a9b79de1`
- `gentle_protocol` test: `8da4bec9c61dfec032fc00771ff784bedf47c04b3c449e56684ed5c20758514f`
- `capability_registry_parity`: `9b0a0c2f46b669f727399190acc0b9f0f6a426b3e271358b1266ef12469d25ed`
- `parity_matrix_freshness`: `2151f93ae445aa172d230172cd91231646d053ef3da7b80dc487649169904516`
- headless library test binary: `9315a9c24023f5734d67eb523e8c8a2b3c66530fa2da7dc6bee5dd49e96ce83b`

## Boundaries

- No source repair, commit, push, dispatch, publication, package build, model call,
  scientific acceptance, native GUI capture, or private-data analysis was performed.
- Historical `b6863907` results remain bound to that SHA and were not re-executed or
  promoted to this candidate.
- Native GUI, Windows behavior, packaging, full-suite acceptance, J1 follow-up, and
  Lifecycle Slice 6 remain separate from this scoped registry repair.

