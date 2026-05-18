# Release Notes / Changelog: `v0.1.0-internal.8` (draft)

This draft internal release is the **agent-surface formalization** cut.
Compared with `v0.1.0-internal.7`, this cut is materially narrower in new
biology (the v7 backlog absorbed most of that), but materially deeper in the
machinery that lets agents drive GENtle predictably:

- The GUI in-root decomposition started in v7 is now complete: six topical
  module extractions out of `src/main_area_dna.rs` and `src/app.rs`, with
  `.git-blame-ignore-revs` discipline preserving history.
- The workspace split has begun Phase 2: `crates/gentle-engine` is no longer a
  placeholder. `iupac_code` is its first real module, the sequence-loader
  seam now lives off the GUI top-level type, and the empty-scaffold state has
  been resolved.
- The shared engine surface gained a typed, test-enforced capability registry
  with a three-level adapter parity policy (engine capability / shell
  reachability / first-class surfacing), a portable error contract with
  `cause_chain`, and a regenerable parity matrix that distinguishes intentional
  surfacing choices from real implementation gaps.
- The "drive the working-tree binary directly" pattern is now a first-class
  agent dev loop: a wrapper script, a `doctor --agent` pre-flight, and a
  documented loop pattern for any agent (Claude, Codex, or another).
- A new Claude-facing agent interface route lands as a peer to the existing
  Codex flow.
- Splicing-expert isoform overlay and reporter-construct handoff land as
  first-version biology features that exercise the new agent surface.

219 commits in v7 was a long backlog cut. 61 commits in v8 is the
consolidation cut where the agent-surface tooling becomes load-bearing.

## Highlights

- **Engine-side capability registry** in `gentle_protocol` projects to CLI,
  MCP, JS, and Lua adapters. Per-adapter `AdapterSurfacing` ternary
  (`Prominent` / `ShellPassthrough` / `NotApplicable`) makes the three-level
  parity policy testable. `EngineError` now carries `cause_chain` and threads
  it through every adapter boundary.
- **Credible parity matrix** at `docs/gui_cli_mcp_parity.md`, regenerated
  from the registry via `scripts/regenerate_parity_matrix.sh`. Cells use the
  four-state vocabulary (`prominent` / `shell-only` / `n/a` / `gap`). The gap
  triage closed several real gaps and surfaced one residual list small enough
  to triage by hand.
- **In-tree agent dev loop** documented in `docs/agent_dev_loop.md` with
  `scripts/dev-gentle-cli` wrapper and `gentle_cli doctor --agent` pre-flight.
- **Claude agent interface** route added as a peer to the Codex flow.
- **`gentle-engine` crate re-introduced** with `iupac_code` as the first real
  module. The empty-scaffold/placeholder state from v7 is resolved.
- **Sequence-loader seam** removed `crate::app::GENtleApp` as a dependency of
  `crate::engine`. Engine no longer depends on the GUI top-level type.
- **GUI in-root decomposition complete**: six topical extractions
  (`rack_workspace_ui`, `gibson_ui`, `primer_design_ui`,
  `sequencing_confirmation_ui`, `rna_read_mapping_ui`,
  `routine_and_agent_assistant_ui`, `genome_catalog_ui`,
  `sequence_ingress_dialogs_ui`, `main_lineage_ui`). `app.rs` reduced by
  ~31%, `main_area_dna.rs` reduced by ~33% (combined: −36% across both
  monoliths).
- **Wildcard dependencies pinned** and consolidated into
  `[workspace.dependencies]` with strict version pins across the workspace.
- **Linux CI parity** added alongside the existing macOS job.
- **Splicing Expert Isoform Read Overlay V1** and **Isoform Expression
  Overlay V1** land as inspection-surface biology features.
- **Reporter Recommender Data Foundation V1** and **Synthetic Biology:
  Reporter Construct Handoff V1** land as new construct-side biology.
- **`AlignSequences`** is the first operation promoted from the inline-operand
  audit's "should accept inline" bucket to the "inline-ok" bucket.
- **TP73 genome-anchored evidence-viewer release proof** workflow lands as
  the deterministic release-acceptance path.

## Notable Changes by Area

### 1) Adapter Parity Policy and Capability Registry

- New `crates/gentle-protocol` capability registry: `CapabilityDescriptor`
  with per-adapter `AdapterSurfacing` fields (`gui`, `cli`, `mcp`, `js`,
  `lua`, `clawbio`), `mutating: CapabilityMutation`, `inline_operand_ok`,
  `engine_operations`, and a `surfacing_justifications` map. Every adapter
  projects from the same registry rather than maintaining its own list.
- New `AdapterSurfacing` ternary captures the three-level parity policy:
  - `Prominent` — first-class affordance (named tool, CLI subcommand, menu
    item, ClawBio skill op)
  - `ShellPassthrough` — reachable via shell/op route, not promoted
  - `NotApplicable` — surface-incompatible by design, with required
    justification
- Architecture doc gained a new "Adapter parity policy — three levels"
  section (`docs/architecture.md`). Phase C roadmap bullets reference the
  matrix.
- `EngineError` now carries `cause_chain: Vec<String>` plus
  `with_cause(...)` plus a `portable_payload` helper that serializes it
  for adapter boundaries. CLI, MCP, JS, Lua boundaries thread the chain.
- New parity tests: `tests/capability_registry_parity.rs` (eight asserting
  registry/listing/surfacing invariants), `tests/adapter_error_contract.rs`
  (CLI and MCP boundary error-shape and cause-chain preservation),
  `tests/mcp_capability_surface.rs` (glossary <-> MCP tools/list parity with
  an intentionally-excluded list in `docs/agent_interface.md`).

### 2) Parity Matrix Triage and Credible-Gap Behavior

- The parity matrix at `docs/gui_cli_mcp_parity.md` is now regenerated from
  the registry via `scripts/regenerate_parity_matrix.sh`. Cell vocabulary is
  the four-state schema: `prominent` / `shell-only` / `n/a` / `gap`.
- The classification rule was refined twice: first to distinguish
  `ShellPassthrough` from `NotApplicable`, then to treat engine-backed
  operations missing only a named tool as `ShellPassthrough` (reachable via
  generic `op` / `apply_operation`) rather than `gap`.
- A curated `NotApplicable` override table now carries real one-line
  justifications. Auto-generated tautologies have been replaced by
  human-authored reasons.
- The triage surfaced and closed several real parity gaps:
  - `rna-reads export-isoform-triage-tsv` was correctly classified as
    intentionally MCP-excluded.
  - The reporter catalog/recommender command family had its shared-shell
    reachability restored.
  - The glossary itself gained the entries that were inconsistent with the
    matrix.

### 3) Agent Dev Loop and Claude Agent Surface

- New documented loop at `docs/agent_dev_loop.md`: drive the working-tree
  binary directly via `scripts/dev-gentle-cli`, keep state in an explicit
  session file, discover operations via `capabilities`, validate each code
  edit before continuing. Cross-linked from `AGENTS.md`, `README.md`, and
  `docs/agent_interface.md`.
- New `scripts/dev-gentle-cli` wrapper around `cargo run -q --bin gentle_cli`:
  respects `CARGO_TARGET_DIR`, falls back to a shared workspace target so
  parallel worktrees do not all rebuild from cold, accepts `--release` for
  benchmark-style calls, tees stdout to a per-invocation file under
  `$GENTLE_DEV_LOG_DIR`.
- New `gentle_cli doctor --agent` pre-flight route: JSON report on cargo
  availability, current-binary freshness vs. source mtimes, `--state` path
  writability, and `gentle_examples_docs --check` health.
- **Claude agent interface** route added (`a2a8a6c0`, `d52b1836`) as a peer
  to the existing Codex path. Initial fix for the agent-interface window
  open path (`b7c96479`) and an endless-failure-loop fix (`8bce3359`).

### 4) Engine Extraction (Phase 2 of the workspace split)

- **Sequence-loader seam**: `GENtleApp::load_from_file` and the five
  `load_dna_from_*_file` helpers now route through the shared DNA-sequence
  module rather than the GUI top-level type. The engine no longer imports
  `crate::app::GENtleApp`. This was the seam blocking real `gentle-engine`
  extraction.
- The empty `gentle-engine` placeholder was first removed (`ee07803f`), then
  re-introduced (`6faa61ea`) as a real workspace member when ready.
- `iupac_code` is the first execution-side module to move into
  `crates/gentle-engine` (`2be52572`). A thin compatibility shim at
  `src/iupac_code.rs` preserves the existing `gentle::iupac_code::*` import
  path.
- The move commit follows the `.git-blame-ignore-revs` discipline established
  by the GUI decomposition work. Docs updated:
  `docs/workspace_split_plan.md` "Current status" notes Phase 2 has begun,
  `docs/roadmap.md` Phase C carries a bullet.

### 5) GUI In-Root Decomposition (pre-Phase-5)

In-root decomposition of the two monoliths is now complete:

- `app.rs`: 50,437 → 34,616 lines (−31%)
- `main_area_dna.rs`: 37,685 lines (decomposed in this cycle, −33% earlier)
- Combined: −36% across both monoliths

The six topical extractions in this cycle:

- `src/main_area_dna/rna_read_mapping_ui.rs` (`91dfb3b8`)
- `src/main_area_dna/primer_design_ui.rs` (`55954302`)
- `src/main_area_dna/sequencing_confirmation_ui.rs` (`2fd74ac2`)
- `src/app/rack_workspace_ui.rs` (`6090d1d3`)
- `src/app/gibson_ui.rs` (`38b3bfeb`)
- `src/app/routine_and_agent_assistant_ui.rs` (`b14536b0`) — Routine
  Assistant and Agent Assistant combined to avoid the rendering-boundary
  overlap
- `src/app/genome_catalog_ui.rs` (`2135e633`)
- `src/app/sequence_ingress_dialogs_ui.rs` (`f15d97a7`) — broader scope than
  the original UniProt prompt: GenBank, dbSNP, UniProt, Ensembl protein,
  reverse-translation, protease digest, and protein-to-DNA handoff dialogs
- `src/app/main_lineage_ui.rs` (`c8d9a667`)

Each landed as a move-only commit paired with a `docs: record … extraction`
commit. All merge SHAs appended to `.git-blame-ignore-revs`.

### 6) Dependency Hygiene

- All 12 wildcard dependencies in the root `Cargo.toml` pinned (`c3437612`).
- New `[workspace.dependencies]` block consolidates shared deps across the
  root crate and the four populated member crates
  (`gentle-{protocol,render,shell,gui}`). Member manifests switched to
  `{ workspace = true }` references.
- Codex chose strict `=X.Y.Z` pins rather than `^X.Y` caret pins, matching
  the existing `deno_core` cluster pattern. Trade-off: no version-drift
  surprise; security/patch upgrades require explicit manifest bumps.
- Routine dependency update pass (`06ccb69d`) preceded the pinning.

### 7) Inline-Operand Audit Follow-Through

- `AlignSequences` is the first operation promoted from the audit's
  "should-accept-inline-but-doesn't" bucket to the "inline-ok" bucket
  (`1a252e55`). It now accepts inline ASCII operands via `SequenceScanTarget`
  alongside the existing `seq_id` route.
- The audit doc at `docs/inline_operand_audit.md` and the roadmap follow-up
  order remain unchanged: `RenderSequenceSvg`, `RenderRnaStructureSvg`,
  `RenderTfbsScoreTrackCorrelationSvg`, `ComputeDotplot`,
  `ComputeFlexibilityTrack` are the queued follow-ups.

### 8) Splicing Expert Isoform Overlay V1

- **Splicing Expert Isoform Read Overlay V1** (`6965f2f7`): overlay of
  isoform-specific RNA reads on the splicing-expert canvas, deterministic
  across adapters.
- **Isoform Expression Overlay V1** (`b421022b`): expression-derived overlay
  rendering on isoform architecture views.
- **Conservative RNA-read isoform triage TSV export** (`6784dc2a`):
  reviewable per-read isoform-assignment exports with explicit confidence
  rules.

### 9) Reporter Recommender and Construct Handoff

- **Reporter Recommender Data Foundation V1** (`3e642129`): catalog-backed
  reporter family with deterministic shell command surface.
- **Synthetic Biology: Reporter Construct Handoff V1** (`3cb07047`):
  end-to-end handoff path for reporter construct selection through the
  shared engine/CLI route.

### 10) RNA-Read Mapping Refinements

- **`rna-reads show-alignments`** as a read-only batch route (`6b752df4`):
  enables headless agent inspection of alignments without state mutation.
- Pancreas gene-family plot helper gained a log-scale option
  (`537113cc`).

### 11) TP73 Genome-Anchored Evidence-Viewer Release Proof

- The TP73 evidence-viewer proof (`43ebd86c`) is now the deterministic
  release-acceptance workflow. Headless regeneration via
  `docs/examples/workflows/tp73_genome_evidence_viewer_release_proof.json`
  produces sequence, splicing-expert, TFBS SVG, and repeat-materialization
  JSON artifacts from `test_files/tp73.ncbi.gb` plus tiny local repeat,
  Clariom-style array, CUT&RUN-style BED, and TFBS fixtures.
- Public runbook at `docs/tp73_genome_evidence_viewer_runbook.md`. Fixture
  bundle at `test_files/fixtures/evidence_viewer/` with provenance and
  deterministic regeneration notes.
- **TP73 proof-bundle audit handoff** (`710d2ff6`) earlier in the cycle set
  up the audit path.

### 12) ClawBio Consolidation and Restriction-Site UX

- **GENtle-side ClawBio interface consolidation** (`ccf064d2`): consolidated
  ClawBio-facing routes around the shared capability registry, reducing
  duplicated tool-handler logic.
- **Sequence restriction filters to the currently painted block**
  (`9174cded`): GUI restriction-site filtering now respects the active paint
  selection, matching the shared paint-first PCR designer pattern.

### 13) Repository Tidies (post-v7 cleanups)

- **Obsolete uncompiled source tree** `src/obsolete/` removed (`98a6cfc2`).
- **Stray root `lib.rs`** removed (`c0446f6c`) — was a 146 KB orphan
  duplicate of the `gentle-protocol` crate root, unreferenced by any Cargo
  target.
- **Duplicate root generated SVGs** (`tp73_uniprot_projection.svg`,
  `tp73_uniprot_projection_multimode.svg`, `e2f1_uniprot_projection.svg`)
  removed (`0d9d0f01`).
- **Internal release notes** moved into `docs/release_notes/` (`d7f24661`).

### 14) CI, Tutorials, and Documentation

- **Linux CI** parity job added alongside the existing macOS job
  (`d3cc762a`).
- **Tutorials catalogued**: generated-but-not-human-confirmed tutorials are
  now in the catalog-backed listing with appropriate disclosure
  (`e4291a8b`).
- **GeneArt tutorial** (`9bbed40f`, `c313dbe0`).
- **`gentle_cli` source tree revisited** (`326fbd9b`) and **CLI `--help`
  output updated** (`7cafe1ec`).
- **Promoter design hosted embedded-window foreground fix** (`0bc15768`).
- **egui window management investigation** (`c36f93ba`) and **daily
  screen-found documentation drift fix** (`f62dbe49`).

## Release-Facing Known Limitations

- The capability matrix shows `GUI: 0 prominent` across all rows because the
  glossary `interfaces` vocabulary has no `gui-menu` marker and the parity
  test explicitly skips the GUI arm. The matrix correctly reports
  reachability but understates GUI prominence. Adding a `gui-menu` interface
  marker is a planned post-v8 refinement.
- The strict `=X.Y.Z` dependency pins eliminate version drift but require an
  explicit manifest bump for every patch upgrade, including security
  patches for `reqwest`, `image`, and `tempfile`. Trivial to relax if
  patch-cadence friction emerges.
- `gentle-engine` contains only `iupac_code` so far. The natural next moves
  (`amino_acids` depends on `iupac_code`; then `feature_location`,
  `exon_frame`) are queued post-tag.
- The full GUI extraction to `crates/gentle-gui` (Phase 5 in
  `docs/workspace_split_plan.md`) remains intentionally deferred. The
  in-root decomposition completed in this cycle is preparation, not a
  substitute.
- ClawBio coverage of newer biology features (splicing-expert isoform
  overlay, reporter construct handoff) is by-design narrower than the
  CLI/MCP surface. The curated skill set is intentional and tracked under
  the three-level parity policy.

## Install / Package Notes

- macOS release artifact remains `.dmg`.
- Windows release artifact remains `.zip` containing `gentle.exe`.
- Linux native installer packaging remains deferred; release metadata
  continues to default Linux distribution intent to `tarball`.
- Container/runtime story unchanged:
  - `ghcr.io/smoe/gentle_rs:cli` for headless/agent/ClawBio use
  - `ghcr.io/smoe/gentle_rs:gui` or `:latest` for interactive GUI container
    use
- The Linux CI job mirrors the macOS job's `apt-get`/feature matrix; see
  `.github/workflows/ci.yml`.

## Suggested Pre-Tag Smoke Checklist

At minimum before tagging `v0.1.0-internal.8`, run:

```bash
cargo check -q
cargo check -q --features js-interface
cargo check -q --features lua-interface
cargo check --all-features -q
cargo test -q workflow_examples -- --test-threads=1
cargo test -q rna_reads
cargo test -q "engine::tests::"
cargo test -q --test capability_registry_parity
cargo test -q --test adapter_error_contract
cargo test -q --test mcp_capability_surface
cargo run --bin gentle_examples_docs -- --check
cargo run --bin gentle_examples_docs -- tutorial-check
cargo run --release --bin gentle -- --version
cargo run --release --bin gentle_cli -- capabilities
scripts/dev-gentle-cli doctor --agent
```

Because this release touched the capability registry, parity matrix,
adapter error contract, agent dev loop, and GUI decomposition, one short
manual smoke pass should also cover:

- launching the GUI, opening the TP73 evidence-viewer release-proof state,
  exporting one evidence artifact
- opening the RNA-read Mapping workspace and inspecting one saved report
- opening the Splicing Expert and confirming the isoform-read overlay
  renders for one TP73 isoform
- running one `scripts/dev-gentle-cli capabilities` and confirming the
  payload includes a non-empty `capability_registry`
- running one MCP smoke via stdio (`cargo run -q --bin gentle_mcp --` with a
  `tools/list` framed request) and confirming tool count matches the
  registry projection
- confirming `scripts/regenerate_parity_matrix.sh` is idempotent against
  the committed matrix file

## Post-Tag Direction

This tag captures the agent-surface formalization and the start of Phase 2
engine extraction. Natural post-tag continuations:

1. Second engine move (`amino_acids` depends on `iupac_code`, then
   `feature_location` and `exon_frame` as standalone follow-ups).
2. Add a `gui-menu` interface marker to `docs/glossary.json` and an oracle
   in `prominent_projections_match_user_facing_listings` so the matrix can
   make a credible GUI prominence claim.
3. Continue inline-operand audit follow-ups in order:
   `RenderSequenceSvg`, `RenderRnaStructureSvg`,
   `RenderTfbsScoreTrackCorrelationSvg`, `ComputeDotplot`,
   `ComputeFlexibilityTrack`.
4. Splicing-expert RNA-mapping deepening, driven by what the analysis work
   in flight surfaces.
5. Cosmetic: `AGENTS.md` direct cross-link to `docs/agent_dev_loop.md`;
   optional `docs/INDEX.md` if the `docs/` tree (now 40+ files) starts to
   feel hard to navigate.

Each is intended to be one independent PR so the tag-to-tag diff remains
reviewable.

## Release-Shaped Smoke Check Results

*Pending — to be populated when the pre-tag smoke runs are executed.*

Local release-prep checks expected on the tag date:

- `git diff --check`: _pending_
- `cargo check -q`: _pending_
- `cargo test -q workflow_examples -- --test-threads=1`: _pending_
- `cargo test -q rna_reads`: _pending_
- `cargo test -q --test capability_registry_parity`: _pending_
- `cargo test -q --test adapter_error_contract`: _pending_
- `cargo test -q --test mcp_capability_surface`: _pending_
- `cargo build --release --features script-interfaces`: _pending_
  - run with `CARGO_TARGET_DIR=/tmp/gentle_release_smoke`
  - run with `CARGO_INCREMENTAL=0`
- release-artifact checks from the same built target:
  - `gentle --version`: _pending_
  - `gentle_cli capabilities`: _pending_
  - `gentle_js --version`: _pending_
  - `gentle_lua --version`: _pending_
  - `gentle_examples_docs --check`: _pending_
  - `gentle_examples_docs tutorial-check`: _pending_
  - `gentle_mcp --help`: _pending_
  - `scripts/dev-gentle-cli doctor --agent`: _pending_

Manual GUI smoke pass: _pending_.
