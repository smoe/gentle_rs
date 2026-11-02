# Primer3 Interface Parity: Source Audit And Bounded Follow-Up

Audited base: `9c07188d`, 2026-11-03, in the `1ee2` worktree. The supplied
Mistral prompt names `a22698aa` and is not an accurate inventory of this tree.
This is an adapter/discovery task, not proof of Primer3/internal algorithm
equivalence or a new release gate. A read-only Claude consultation was offered;
no new Claude review is claimed.

## Findings

| Prompt claim | Current source evidence | Disposition |
| --- | --- | --- |
| No persistent configuration | `EngineParameters.primer_design_backend` and `primer3_executable`; `SetParameter` handles set/reset, serialized `ProjectState` retains them. | Keep the existing project-scoped configuration; do not create competing defaults or a global preferences store. |
| No doctor/version probe | `ShellCommand::PrimersPreflight` calls `GentleEngine::primer3_preflight_report`; CLI already forwards `primers preflight`. | Use the existing command and `gentle.primer3_preflight.v1`, not a `primers doctor` alias. |
| GUI configuration/preflight and qPCR are missing | `main_area_dna/primer_design_ui.rs` contains backend controls, Apply/Probe/Reload and `build_design_qpcr_operation`. The 2026-10-08 changelog already corrected the preflight claim. | Preserve the existing PCR Designer; no second configuration dialog or status-bar subsystem. |
| No parity mechanism | `tests/capability_registry_parity.rs`, `tests/mcp_capability_surface.rs`, `tests/parity_matrix_freshness.rs`, adapter tests and tutorial infrastructure exist. | Extend existing tests, not another shell harness with placeholder hashes. |
| No Primer3 cancellation | `test_auto_primer3_progress_cancellation_does_not_fall_back_or_commit` and bounded-runtime/process tests already exist. | Do not rewrite cancellation. Coverage of native platforms and preflight timeouts is a separate question. |
| All MCP design tools must be duplicated | Generic `op`/`workflow` routes expose `DesignPrimerPairs`, `DesignQpcrAssays` and `SetParameter`, with explicit confirmation. | Dedicated tool names are optional surfacing, not missing engine capability. |
| Inner agent needs special commands | `agents execute-plan` dispatches reviewed shell candidates through the shared executor; `primers preflight` is already a documented agent route. | Test deterministic plan execution; do not invoke a model, add `/agent` grammar or bypass confirmation. |
| Preflight needs adapter discovery | The pre-change matrix lists MCP/JS/Lua gaps for the shell-only preflight helper. Generic agent-plan execution can reach it indirectly. | Add small named readers and truthful registry metadata, not a new engine operation. |

The architectural distinction matters: reachability, prominent controls and
algorithm equivalence are different claims. Adding preflight readers improves
discovery and direct access; it does not mean every primer-related matrix gap
is closed.

## Implemented Scope

- MCP `primers_preflight` and JS/Lua `primers_preflight(project, backend?,
  primer3_executable?)` all use the existing shell parser/executor.
- Omitted arguments use project settings; overrides apply to this probe only.
  No project save, design report, configuration migration or fallback policy
  change is introduced. Probe invocation runs the selected local executable;
  read-only describes GENtle project/file persistence, not a sandbox guarantee
  for a user-selected binary.
- Existing `auto`/`internal`/`primer3` semantics and DEC-037 explicit-path
  authority remain unchanged. Missing executables remain structured diagnostic
  results, not false scientific failure/success or an implicit installation.
- Registry, parser, adapter and GUI tests exercise these routes and existing
  parameter persistence. Synthetic absent absolute paths avoid host PATH
  dependencies. JS/Lua tests include their actual VM bindings when enabled.
  Agent tests exercise the executor, not live language-model behavior.
- Compare full preflight records on the same host except elapsed probe time.
  Paths, OS errors, versions and timings make universal output hashes invalid.
  Preserve biological byte/hash contracts; no scientific fixture is normalized.
- Regenerate `docs/gui_cli_mcp_parity.md` using the existing protocol renderer.
  The glossary interface additions also require refreshing its bound SHA-256 in
  the generated maturity catalog and three ClawBio projections via
  `python3 scripts/generate_capability_maturity.py --write`. These four files
  change only that source hash, not capability maturity or scientific evidence.

## Rejected Expansion And Remaining Work

Do not adopt the prompt's new default Tm/GC/length/product values, fallback flag,
cached report in engine parameters, synchronous new GUI probe, fabricated menu
paths, placeholder SHA values or unconditional system-Primer3 CI requirement.
Those would change scientific/configuration semantics or create untested paths.
The illustrative 20-base template is not an adequate ordinary PCR design demo.

The existing GUI probe is synchronous under an engine read lock, and executable
version/help probes use blocking subprocess output. That is a real responsiveness
and timeout follow-up under DEC-048, not evidence of missing preflight. It needs
held-worker, stale-result, timeout/termination and native-platform tests before
claiming asynchronous/cancellable setup. This reader slice does not solve it.

Separate two later teaching goals: an offline synthetic setup/parity walkthrough,
and an opt-in real-Primer3 design replay with a pinned binary/input contract.
Reuse tutorial source/catalog/acceptance machinery and existing Primer3 fixtures;
do not invent tutorials 04.10/04.11 with non-existent commands or fixed hashes
of installation diagnostics. For now, the CLI setup guide teaches the real
commands and distinguishes reachability from assay validity.

Deeper Primer3 constraint mapping/thermodynamics and internal-versus-external
backend equivalence remain separate engine/report work under DEC-022.

## Changed Files

- Readers: `src/mcp_server/primer_tools.rs`, `src/mcp_server.rs`,
  `src/js_interface.rs`, `src/lua_interface.rs`.
- Test coverage: those adapter modules plus `src/main_area_dna/tests.rs`,
  `src/test_support.rs`, `tests/capability_registry_parity.rs`.
- Discovery: `crates/gentle-protocol/src/lib.rs`, `docs/glossary.json`,
  generated `docs/gui_cli_mcp_parity.md`.
- Documentation: this audit, `docs/cli/primers.md`, `docs/agent_interface.md`,
  `docs/protocol.md`, `docs/roadmap.md`, `docs/CHANGELOG.md`.
- Required source-hash refresh: `assets/capability_maturity.json`,
  `integrations/clawbio/skills/gentle-cloning/CAPABILITY_MATURITY.json`,
  `integrations/clawbio/skills/gentle-pcr-primer-design/CAPABILITY_MATURITY.json`,
  `integrations/clawbio/skills/gentle-tss-collection/CAPABILITY_MATURITY.json`.

## Validation

Run targeted checks at the final source revision, retaining platform/features:

```sh
cargo test --locked -p gentle-protocol --lib
cargo test --locked --lib primers_preflight
cargo test --locked --lib primer3_preflight
cargo test --locked --features script-interfaces --lib primers_preflight
cargo test --locked --test capability_registry_parity --test mcp_capability_surface --test parity_matrix_freshness
python3 -m unittest scripts.test_tutorial_checkouts
cargo check -q --locked
```

Native Windows, live GUI, live inner-agent and real-Primer3 results must be
reported separately from synthetic/local checks. No output equality across
different Primer3 builds or platforms is promised.

### Local Results

Checked on Darwin arm64 with `rustc 1.100.0-beta.7 (fd8e211ae 2026-10-17)`,
base `9c07188d` plus this patch:

- `cargo test -q --locked --offline -p gentle-protocol --lib`: 120 passed,
  including canonical parity-projection freshness after regeneration. A temporary
  protocol-only Rust executable called
  `gentle_protocol::render_gui_cli_mcp_parity_matrix_markdown()`; no alternative
  matrix-generation logic was introduced.
- `python3 -m unittest scripts.test_tutorial_checkouts scripts.test_capability_maturity -q`:
  98 passed, including LF/CRLF source-binding checks.
- `python3 scripts/generate_capability_maturity.py --check`: passed.
- `CARGO_INCREMENTAL=0 CARGO_PROFILE_DEV_DEBUG=0 cargo check -q --locked --offline --features script-interfaces --tests -j 2`:
  passed. The initial sandbox attempt could not write Swift caches or download
  V8; the permitted native-dependency retry completed successfully.
- Runtime tests used `CARGO_INCREMENTAL=0`, `CARGO_PROFILE_DEV_DEBUG=0`,
  `CARGO_PROFILE_TEST_DEBUG=0`, `--locked --offline --features script-interfaces`
  and `-j 2`:
  - `cargo test --lib primers_preflight -- --test-threads=1`: 7 passed,
    including actual JS/Lua VM bindings and MCP/agent execution.
  - `cargo test --lib primer3_preflight -- --test-threads=1`: 7 passed,
    including GUI helper/settings persistence. One is an opt-in external-binary
    test whose body returned without probing because
    `GENTLE_TEST_EXTERNAL_BINARIES` was unset; this is not a real-Primer3 pass.
  - `cargo test --test capability_registry_parity --test mcp_capability_surface --test parity_matrix_freshness`:
    26, 6 and 1 passed, respectively.
- The macOS linker warned about the large `__eh_frame` section exceeding compact
  unwind-table capacity. Linking and test execution succeeded; this slice does
  not change linker settings.
- `CARGO_INCREMENTAL=0 CARGO_PROFILE_DEV_DEBUG=0 cargo clippy -q --locked --offline --features script-interfaces --lib -j 2`:
  exit 0, with repository-wide warnings (including large enum variants,
  collapsible conditions and argument counts). No unrelated lint cleanup.
- `git diff --check`: passed. Scoped formatting checks passed for the new
  primer-tool tests, shared test helper, protocol registry and parity test.
  Checking the remaining touched modules still reports pre-existing observed-SAM
  and help-test formatting differences; these were left unchanged.
- Session-close check: 4 OK, 2 warnings, no failures before commit. Warnings
  cover the intended dirty/generated files and manual scope review; the audit
  above records both the regeneration and deliberate narrowing of the proposal.

Not run: the full library suite, native Windows tests, live GUI interaction,
live model calls, or a pinned real-Primer3 design replay. The scoped synthetic
tests do not establish scientific backend equivalence.
