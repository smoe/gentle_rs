# GENtle Testing Strategy (Draft)

This document proposes automated testing for both operation correctness and
rendering correctness.

## 1. Test pyramid

1. Engine unit tests (fast, deterministic)
2. Protocol/CLI integration tests (JSON in/out, state persistence)
3. Rendering snapshot tests (linear + circular graphics)
4. Optional GUI smoke tests (window opens, toggles wired)

## 2. Engine tests (required)

Scope:

- operation semantics (`LoadFile`, `Digest`, `Ligation`, `Pcr`, `ExtractRegion`,
  `FilterByDesignConstraints`)
- structured errors/warnings
- deterministic operation log shape
- display toggle operations (`SetDisplayVisibility`)

Execution:

```bash
cargo test engine::tests::
```

## 3. CLI/protocol tests (required)

Scope:

- JSON contract compatibility
- `@file.json` input handling
- state persistence (`--state`)
- command parity (`op`, `workflow`, `state-summary`, import/export)

Recommended method:

- add integration tests that run `gentle_cli` as a subprocess
- compare normalized JSON output against golden files

### 3.1 Canonical workflow example tests

Canonical protocol examples in `docs/examples/workflows/*.json` are now part of
the test surface.

Execution:

```bash
cargo test workflow_examples -- --test-threads=1
```

Online-only examples are opt-in:

```bash
GENTLE_TEST_ONLINE=1 cargo test workflow_examples -- --test-threads=1
```

`test_mode` policy from each example file:

- `always`: parsed, validated, and executed in default runs
- `online`: executed only with `GENTLE_TEST_ONLINE=1`
- `skip`: parsed/validated only
- `optional_blast`: native tutorial checks run only with installed BLAST+;
  missing tools are skipped, while documentation is always checked (see 3.2).

Primer-design CI policy:

- Always-run workflow/tutorial examples that exercise primer design should pin
  `primer_design_backend=internal` explicitly so routine CI stays fast and does
  not depend on a locally installed `primer3_core`.
- Real external-binary primer-design checks are opt-in via:

```bash
GENTLE_TEST_EXTERNAL_BINARIES=1 cargo test real_primer3 -- --test-threads=1
```

- Those opt-in tests are intended to validate the actual `primer3_core`
  integration path separately from the default deterministic internal-backend
  CI surface.

MaxEntScan reference-parity policy:

- Default tests use synthetic unrestricted tables and do not require or bundle
  MaxEntScan resources.
- Developers with a reviewed, unpacked MaxEntScan distribution can compare
  GENtle's normalized native scorer directly with upstream `score5.pl` and
  `score3.pl`:

```bash
GENTLE_TEST_MAXENTSCAN_DIR=/path/to/MaxEntScan \
  cargo test native_maxent_scoring_matches_opt_in_reference_scripts \
  -- --test-threads=1
```

- `GENTLE_TEST_PERL_BIN` may name a non-default Perl executable. The test is
  skipped when `GENTLE_TEST_MAXENTSCAN_DIR` is unset, bounds each reference
  scorer invocation to 15 seconds, and never copies the supplied tables into
  repository fixtures.

### 3.2 Tutorial drift and runtime checks

Tutorial source + generated outputs are part of the test surface:

- source units: `docs/tutorial/sources/`
- generated runtime manifest: `docs/tutorial/manifest.json`
- generated output: `docs/tutorial/generated/`

The conservation tutorial uses `optional_blast`: its real runtime check runs
only when `makeblastdb`, `blastdbcmd`, and `blastn` are already available.
Missing tools produce an explicit skip, never an installation or a claim of
BLAST acceptance. Installed-but-broken tools and runtime failures still fail.
`tutorial-check` always checks all generated documentation, then performs this
optional runtime check separately. Documentation generation itself never runs
`optional_blast` workflows, so committed output is independent of local tools;
such chapters cannot retain runtime artifacts during generation.

A [tiny TP73 CDS database fixture](../test_files/fixtures/blast_tp73_isoforms/README.md)
ships with the source: three ENA-derived DeltaN isoforms, 4,398 bases total.
Its provenance test needs no BLAST; its native test builds a temporary index,
inspects it and checks all three exact self matches, or explicitly skips when
tools are absent. No downloads or prebuilt platform-specific indexes are used.
Run both fixture tests with `cargo test --lib tp73_blast_fixture -- --nocapture`.

Validation commands:

```bash
cargo run --bin gentle_examples_docs -- tutorial-check
cargo run --bin gentle_examples_docs -- tutorial-manifest-check
cargo run --bin gentle_examples_docs -- tutorial-catalog-check
cargo test workflow_examples -- --test-threads=1
```

`tutorial-check` regenerates tutorial outputs in a temp directory and compares
byte-for-byte with committed files under `docs/tutorial/generated`. On failure,
it prints a paste-ready tutorial feedback context with the chapter id, source
JSON, workflow JSON, generated chapter path, artifact directory, failing check,
and suggested GitHub issue-template category.

CI additionally runs a CLI smoke path for core tutorial chapters via:

- `cargo run --bin gentle_cli -- workflow @docs/examples/workflows/<core>.json`

### 3.3 Planned opt-in Mistral inner-agent conformance test

Add one live-provider integration test for the native Mistral Agent Assistant
transport. This is a provider conformance test, not a deterministic-output
snapshot and not part of default offline CI.

Activation and isolation:

- If `MISTRAL_API_KEY` is missing or blank, print one explicit skip line and
  return before model discovery or generation. Do not use `#[ignore]`, because
  that would also suppress the test when the key is deliberately supplied.
- When the key is present, use the existing `mistral_large_native` transport
  with `--model mistral-small-latest` by default. Allow
  `GENTLE_TEST_MISTRAL_MODEL` to override the model without changing the
  checked-in agent catalog.
- Use a temporary project state and report directory. Never print the API key,
  persist it in project state, or include it in failure snapshots.
- Disable retries and use a bounded request timeout so a provider outage cannot
  stall the suite indefinitely.

Test routine:

1. Run native-Mistral preflight and require successful authentication, model
   discovery, and selected-model resolution.
2. Submit the fixed, offline restriction-map prompt from the Agent Interfaces
   tutorial. Request one `execution=ask` suggestion and no network-backed
   biological operation.
3. Validate `gentle.agent_response.v1`, require exactly one non-recursive
   command, and parse it through `parse_shell_line`.
4. Resolve the suggested command's canonical capability descriptor and require
   `annotation_status=fact_annotated`; reject a command not represented by
   introspection.
5. Bind its arguments and require `introspect readiness` to report `ready`
   against the temporary state before execution.
6. Execute through the shared shell path with auto-execution disabled, then
   verify declared hard effects through `introspect verify-effects`.
7. For the restriction scan, compare normalized scientific rows with a direct
   shared-shell run: enzyme, recognition interval, recognition sequence, cut
   geometry, and matched-site count. Ignore prose, provider token usage,
   elapsed time, generated timestamps, and option ordering.
8. Run a second state-aware prompt that refers to the now-loaded `agent_toy`
   sequence. Require the suggestion to use that real sequence id rather than
   inventing another project object, then repeat parse/readiness checks.

The pass condition is deterministic containment: a variable model proposal is
accepted only when it lands inside GENtle's fixed schema, parser, capability,
readiness, execution, and effect-verification contracts. Repeating the provider
call and demanding byte-identical prose or JSON would test provider sampling,
not GENtle determinism.

## 4. Rendering tests (required)

Graphics are part of functionality. Visibility toggles and biology overlays must
be tested at the rendered output level.

### 4.1 Proposed renderer test contract

Introduce deterministic rendering export functions from a shared view model:

- `render_linear_svg(view_model) -> String`
- `render_circular_svg(view_model) -> String`

Use fixed viewport size, fonts, and deterministic layout order.

### 4.2 Snapshot workflow

For each test case:

1. Build state through engine operations
2. Build view model with selected display settings
3. Export linear and circular SVG
4. Compare against approved snapshot files

You can export from persisted engine state with:

```bash
cargo run --bin gentle_cli -- --state state.json render-svg SEQ_ID linear out.linear.svg
cargo run --bin gentle_cli -- --state state.json render-svg SEQ_ID circular out.circular.svg
```

Suggested snapshot folders:

- `tests/snapshots/linear/*.svg`
- `tests/snapshots/circular/*.svg`

### 4.3 Cases to include first

- baseline map with all overlays on
- each visibility toggle off/on individually
- dense feature labels
- linear strand direction correctness
- restriction site label crowding
- circular zero-point crossing features
- design-constraint filter pass/fail cases (GC bounds, homopolymer cap, U6
  `TTTT` rejection, forbidden motifs)

## 5. Regression gates

Recommended CI gate (minimum):

1. `cargo check -q`
2. engine unit tests
3. CLI integration tests
4. rendering snapshot tests

If snapshot output changes, require reviewer approval and explicit snapshot update.

### 5.1 External performance audit

Criterion benchmarks characterize CPU-bound, in-process GENtle algorithms;
they do not replace functional assertions or end-to-end acceptance. The
primary target exercises the actual eager DNA-window constructor, deferred
UI-thread hydration boundary, and first and steady embedded egui frames over a
feature-free control and the public TP73 locus:

```bash
cargo bench --profile bench-audit -p gentle-benchmarks \
  --bench gui_operations -- --quick --noplot
```

`bench-audit` keeps the release profile's `opt-level=3` and stripping, but uses
thin LTO, 16 codegen units, and `panic=unwind` for practical routine audits.
Cargo forces unwind for benchmark targets regardless of the release panic
setting. The dedicated, non-published `gentle-benchmarks` package depends on
GENtle only as a library with default features disabled, while explicitly
enabling the `desktop-gui` modules imported by the GUI benchmark together with
`benchmark-support`. Cargo therefore compiles the measured GUI library code but
does not prepare the root package's application binaries before sampling. Plain
`cargo bench -p gentle-benchmarks --bench gui_operations ...` remains the exact
release-like fat-LTO, one-codegen-unit mode. The two modes characterize
different binaries and their results must never be compared. A cold library
build and link can remain substantial; the custom profile is intended to bound
that pressure and support cacheable repeated audits, not to promise an
instantaneous first build. Build metadata tracks the current loose Git branch
ref; repository-wide `packed-refs` is only tracked as a fallback so unrelated
worktree maintenance does not discard the audit cache.

Implementation-side verification stops at the non-statistical smoke path:

```bash
cargo test -p gentle-benchmarks --bench gui_operations
cargo test -p gentle-benchmarks --bench specificity_finalization
```

This is a CPU-side proxy for GUI-critical preparation and painting. It cannot
measure native viewport creation, GPU/compositor behavior, event delivery, or
perceived interaction stalls. Glen is the named external auditor for the
current release assessment and therefore combines those results with an
exact-revision GUI run using the existing `gui-profiler` feature/Puffin scopes
and the DNA-viewer resize/repaint harness. The auditor, not the implementation
author, retains baselines and assigns the performance verdict.

The secondary target covers complete-query BLAST HSP parsing and
primer-specificity interpretation at 100, 1,000, and 6,600 total rows:

```bash
cargo bench --profile bench-audit -p gentle-benchmarks \
  --bench specificity_finalization -- --quick --noplot
```

The benchmark constructs its synthetic input outside timed loops, verifies the
expected hit, pairing, intended-product, and off-target counts once, reports
throughput in HSP rows, and includes a fixture SHA-256 token in every benchmark
ID. It launches no BLAST process and reads no real reference dataset. Cargo
separates compiled custom-profile artifacts under `target/bench-audit/`, but
Criterion otherwise shares `target/criterion/` across profiles. Retained audit
runs therefore set `CRITERION_HOME` to the matching profile directory and put
the profile name in every baseline label. Results should only be compared when
runner, toolchain, fixture hash, GENtle revision, and optimized profile all
match. Pull-request CI may compile or smoke the targets, but strict timing
thresholds belong to repeated external-auditor runs on a stable host after
several baselines exist. Exact commands and retained evidence are documented
in [`../benches/README.md`](../benches/README.md).

Full Ensembl FASTA retrieval, peak RSS, real BLAST/Primer3 execution, GUI frame
behavior, and Xvfb/macOS acceptance remain system-level measurements rather
than unit-test claims.

## 6. Semantic GUI acceptance identifiers

GENtle exposes a bounded, read-only semantic widget snapshot only when built
with `--features gui-test-support`. The feature is disabled by default. Set
`GENTLE_GUI_TEST_SNAPSHOT` to an explicit output path when launching `gentle`;
the app atomically replaces that JSON file after each root frame. Root-frame
publication retains the latest paint from every still-active native child
viewport and removes entries once egui no longer reports that viewport. There
is no listener, input-injection API, screenshot permission bypass, or command
route.

Semantic IDs are lowercase dotted purpose names such as
`dna.splitter.info_width` and do not derive from translated labels or screen
position. Agent-help controls receive their semantic window identity explicitly
from the rendering caller; translated or dynamic titles are used only to derive
the pseudonymous subject scope and never to guess a window kind. Repeated domain
rows share one semantic ID and carry a pseudonymous `subject_scope` derived from
stable domain fields rather than list index. V2
uses a domain-separated SHA-256 digest truncated to 128 bits; this prevents
literal-value disclosure but is not an anonymity claim for guessable inputs.
Snapshots never include visible labels, sequence strings, paths, credentials,
or arbitrary user text. Rectangles are screen-relative egui logical points,
including the native viewport origin; multiply them by the item's
`pixels_per_point` for physical pixels.

Example launch:

```bash
cargo build --locked --features gui-test-support --bin gentle
GENTLE_GUI_TEST_SNAPSHOT=/tmp/gentle-gui.json \
  target/debug/gentle --project /path/to/tp73-proof.project.gentle.json
jq '.generation, .settled, .items[] | select(.semantic_id == "dna.splitter.info_width")' \
  /tmp/gentle-gui.json
```

The snapshot schema is `gentle.gui_semantic_snapshot.v2`. `generation`
advances monotonically and `settled` is false while known deferred GUI work,
including initial feature-tree hydration, is pending. Harnesses should wait on
those fields plus the semantic controls required by the scenario rather than a
fixed sleep. Independently scheduled DNA child-viewport repaints replace that
viewport's previous rectangles in the aggregate snapshot; duplicate IDs within
one paint remain an assertion failure. Semantic IDs locate GUI affordances;
scientific correctness still comes from the typed project state and reports.
The Linux example [`../scripts/gui_semantic_xvfb_smoke.sh`](../scripts/gui_semantic_xvfb_smoke.sh)
self-skips when X11 tools or a display are unavailable, uses the read-only
`main.project.sequence.open` rectangle for one ordinary X11 double-click, then
waits for the DNA viewer and its information-width splitter. It accepts
explicit pseudonymous sequence/repeat/array row scopes plus a typed-state
verification command.

### 6.1 Tutorial GUI acceptance contracts

Tutorial source units may carry an optional
`gentle.tutorial_gui_acceptance.v1` block. The documentation generator copies
and validates this block, but the runtime GUI and engine never read it. The
contract names a deliberately incomplete starter workflow, a completed oracle,
semantic GUI targets, closed interaction kinds, and typed postconditions. Its
Rust representation has no raw command, script, shell, or executable-text
field. Every target must also exist in the closed tutorial-control catalog,
which declares the owning window, allowed interaction, persistence class,
scientific-effect authority, and any constrained replacement-text policy.
Snapshot-only or command-capable controls therefore cannot become tutorial
input targets merely by naming their semantic ids. The external runner owns
input delivery and the final verdict.

The acceptance `profile` is also closed: `smoke`, `offline-core`, or `full`.
Misspellings fail source/manifest parsing instead of silently creating another
profile.

Persistence and scientific proof are separate. A step declares
`persists_project_state: true` whenever the ordinary GUI action dirties the
project, including auxiliary GUI metadata. Only a step with
`scientific_effect: true` declares `before` and `after` fact expressions and
must include at least one non-visual typed verifier. Deterministic tests require
each non-preexisting `before` expression to be `Unsatisfied` on the starter,
not `Unknown`, and require the completion condition and every `after`
expression to be `Satisfied` on the oracle. This prevents both an already
completed starter and a view-only GUI change from producing a meaningless
green scientific result. At runtime the eventual external runner must repeat
those checks around each scientific effect and save whenever a
project-persisting step leaves `main.project.save_state` unsaved.

An explicitly declared `view_only: true` contract has the complementary
invariant: its completion condition is already `Satisfied` in both starter and
oracle, every step has view-state authority, and no step may persist project or
scientific state. This supports provenance-bound navigation and parameter-review
tutorials without dressing “the window opened” up as a scientific result. The
runner repeats the satisfied invariant after the final interaction to prove
that inspection did not replace or remove the bound project content.

Text entry uses `replace_text`, never implicit append semantics. The catalog
currently permits it only for selection formulas and bounded identifiers, with
format validation appropriate to each target. Report verifiers reject vacuous
field lists and may add typed value or field-to-field assertions. The Simple
PCR contract consequently checks a non-empty pair set, the persisted
`require_roi_flanking` constraint, and first-pair coordinates on the correct
sides of the ROI rather than treating report existence alone as completion.

GUI state and scientific state remain independent evidence channels. Visible
claims use only `gentle.gui_semantic_snapshot.v2`; facts, expected effects,
reports, and state checks use only a project saved by the ordinary GUI action.
For a project opened from a known path, `Ctrl+S` and `File -> Save Project` save
in place without invoking a native file dialog. Harnesses should wait for
`main.project.save_state` to change from `unsaved` to `saved` before inspecting
that project with `gentle_cli --project`. Semantic response rectangles use
egui's clipped interaction rectangle; controls outside a scroll viewport report
`visible=false` and zero geometry. Specialist windows report `ready` only after
their active sequence content renders; blocked windows remain visible but
disabled with outcome role `blocked`. Artifact checks inspect the artifact
directly and do not force an otherwise unnecessary project save.

The first contract is attached to `simple_pcr_selection_gui`. It opens the
starter sequence from the main project graph, uses a selection formula rather
than pixel-coordinate dragging, enters the PCR Designer through the ordinary
selection context action, and proves creation of the closed-world
`primer_design_report.exists` fact against the separate
`simple_pcr_primer_design_offline` oracle. The first main-window action is part
of the typed contract; a harness must not rely on an unrecorded setup click.

The Simple-PCR GUI starter now designs on the same 800-base TP73 extract as
the oracle, not the full locus. Its fixed formula `=201 .. 600` creates a
400-base core with 200-base flanks. The full locus remains only as provenance.
`max pairs = 5` bounds report size, not search effort. The ten-minute compute
budget is unchanged: all three smoke chapters (branch/reverse complement,
digest, Simple PCR) must complete, and a timeout remains a failure. Rerun live
Linux acceptance on the exact `.10` candidate before calling this smoke green.

`scripts/tutorial_gui_acceptance.py` is the external Linux/X11 runner. It uses
ordinary `xdotool` events and the read-only semantic rectangles; GENtle does
not inject input or certify its own GUI. Before launching the GUI, the runner
asks `gentle_examples_docs tutorial-gui-project` to materialize the incomplete
starter and separate completed oracle with path rewriting. That helper also
returns the engine-owned, pseudonymous sequence scopes, so the Python code does
not duplicate subject-identity hashing. The runner then independently repeats
starter/oracle fact checks, drives only the closed interaction vocabulary,
saves through ordinary `Ctrl+S`, and verifies facts/reports/state/artifacts
through fixed `gentle_cli` routes. Tutorial prose and fields can never supply a
shell command.

Build the GUI with semantic test support and run the smoke profile inside an
explicitly isolated Linux network namespace:

```bash
cargo build --locked --features gui-test-support \
  --bin gentle --bin gentle_cli --bin gentle_examples_docs

xvfb-run -a -s "-screen 0 1600x1000x24" \
  sh -c 'openbox >/tmp/gentle-tutorial-openbox.log 2>&1 & \
    parent_netns=$(readlink /proc/self/ns/net) && \
    exec unshare --user --map-root-user --net -- \
      python3 scripts/tutorial_gui_acceptance.py \
        --repo-root . \
        --profile smoke \
        --evidence-dir /tmp/gentle-tutorial-gui-smoke \
        --parent-network-namespace "$parent_netns" \
        --network-enforcement linux_network_namespace'
```

The runner requires `xdotool`, `xdpyinfo`, `xprop`, `xwininfo`, and a named EWMH
window manager; the example uses Openbox. It also requires `scrot` whenever a
selected contract marks a screenshot as required. The parent network namespace
identity is captured before `unshare`; the runner compares it with its own
namespace instead of assuming that `/proc/1/ns/net` is readable. Use repeated
`--chapter ID` instead of `--profile` for a bounded chapter set. `smoke`
currently contains Simple PCR, branch/reverse-complement, and BamHI/EcoRI digest
contracts.
Profiles select exactly their declared contracts, not a cumulative tier.
`offline-core` currently contains the conservation/promoter-similarity view-only
chapter; inspecting its prepared result does not prove a new BLAST computation.
`full` has no contracts yet. Online chapters remain explicit and authorized.

The two cloning contracts start from load-only workflows
`branch_gui_starter` and `digest_gui_starter`, not completed results.
`branch_gui_oracle` deliberately uses the GUI's default `_revcomp` ID;
the original scripted tutorial's explicit `_rc` name remains supported.
The digest run uses the ordinary Sequence Tools controls and the explicit
`frag` prefix. A `state` verifier with `compare_with_oracle: true` hashes
saved bases, topology, feature annotations, molecule type and end geometry
against the separately executed oracle. Labels and runtime caches are not
scientific comparison fields. Headless source/oracle tests and Python verifier
tests do not constitute a live Xvfb pass.

Each chapter receives a clean `HOME`, XDG roots, temporary directory, and
starter project. Inherited `GENTLE_*` and `*_API_KEY` variables are removed from
the GUI process and represented only by names plus value hashes in
`environment.json`. An offline contract with
`--network-enforcement not_enforced` fails as `harness_gap`; the runner never
claims isolation from Xvfb alone. The retained
`gentle.tutorial_gui_acceptance_ledger.v1` records every requested interaction,
resolved semantic item, emitted X11 event, before/after snapshot generation,
typed verifier result, project/report/artifact hash, exact workflow-source
hash, screenshot, timeout class, and failure class. `environment.json` records
the X display and window manager alongside required and optional tool versions.
`acceptance-report.json` contains the external aggregate verdict. The older
shell Xvfb script remains only a liveness example.

When a step retains a screenshot, the runner preserves one untouched X11-root
PNG and writes a `gentle.tutorial_gui_screenshot_evidence.v1` sidecar. The
sidecar binds the image to the exact source revision and GUI binary hash,
tutorial manifest and acceptance-contract hashes, chapter/prose/step identity,
semantic snapshot generation and canonical hash, pseudonymous subject scope,
logical rectangle, pixel scale, exact native X11 client id and root-screen
client geometry, egui viewport identity, coordinate transform, and physical-pixel
focus rectangle. Embedded semantic surfaces inherit the exact native client of
their owning egui viewport rather than inventing another X11 window. The
runner subtracts the semantic window origin and adds the native client origin,
so window-manager decorations cannot displace input or teaching crops. It also
records the capture backend and timestamp. Two lossless SVG teaching views are
derived from that same raw PNG: a whole-screen orientation view and a padded
interaction-context view, each with a numbered outline around the semantic
focus. Their crop transforms and hashes remain in the sidecar. The raw capture
is evidence; the SVGs are explanatory projections and never replace typed
scientific verification.

Project-metadata edits are deferred across dependent field/button interactions,
then saved immediately *before* the next scientific action (or at chapter end).
The scientific step must produce a fresh unsaved transition before the runner
saves and verifies it. An old dirty flag can therefore never masquerade as
scientific completion. The ledger records both the deferred state and the
ordinary pre-/post-scientific `Ctrl+S` events; typed verification reads only
the saved project.

Pointer actions are emitted only after the exact semantic target reports
`hovered=true` from a child/root viewport publication. Text replacement leaves
the field with `Tab` before checking its typed outcome. A tutorial may declare
a bounded `scroll` interaction when a real user must navigate a scroll area;
the direction is typed and wheel repetitions are limited to 1..=20. Navigation
is recorded in the same action ledger rather than occurring as a hidden runner
convenience. The default compute timeout is ten minutes; exceeding it is a
product-performance result, not a coordinate or focus harness gap.

### 6.2 Candidate-bound tutorial gate for external auditors

Glen can run the tutorial gate on his Linux host; it does not require a GitHub
GUI job. The gate does not fetch, build, install tools, change the source tree,
publish screenshots, repair code, or execute commands interpreted from prose.
Its verdict is **Linux/X11 offline tutorial acceptance**, not whole-release
approval, private biological acceptance, or proof of macOS/Windows behavior.

First inspect coverage, offline on any host, without launching GENtle:

```bash
python3 scripts/tutorial_acceptance.py inventory --output /tmp/tutorial-coverage.json
```

`gentle.tutorial_acceptance_coverage.v1` joins the existing generated catalog,
manifest and workflow definitions. It distinguishes `gui_scientific`,
`gui_view_only`, `workflow_only`, `manual_uncovered`, and `reference`. Workflow
test modes and required-file availability remain separate from GUI coverage.
Every inventory row starts `execution_status: not_run`; source metadata,
historical review labels and file availability are not passes. The inventory
is derived, not another authoring catalog: add contracts to tutorial source
units and regenerate through the existing helper. Chapters without GUI
contracts remain visible, including manual walkthroughs without an oracle.

For acceptance, freeze a full commit SHA, use a clean checkout, and build all
three binaries from it with the command in section 6.1. `--gentle`,
`--gentle-cli`, and `--examples-docs` can select pinned binaries outside the
checkout. All three must report the exact candidate in `--version`; the
documentation helper now supports that flag too. Bind their bytes, not just
their names. Store all evidence outside the checkout in a **new** directory.
Use a separate clean worktree when the development checkout holds private
analysis files; do not remove those files merely to satisfy the gate.

Run both populated profiles inside the isolated X11 session:

```bash
export CANDIDATE_SHA=FULL_40_CHARACTER_COMMIT_SHA
export ACCEPTANCE_DIR=/absolute/private/path/to/new-tutorial-run
xvfb-run -a -s "-screen 0 1600x1000x24" \
  sh -c 'openbox >"${ACCEPTANCE_DIR}.openbox.log" 2>&1 & \
    parent_netns=$(readlink /proc/self/ns/net) && \
    exec unshare --user --map-root-user --net -- \
      python3 scripts/tutorial_acceptance.py run \
        --repo-root . --candidate "$CANDIDATE_SHA" \
        --profile smoke --profile offline-core \
        --evidence-dir "$ACCEPTANCE_DIR" \
        --parent-network-namespace "$parent_netns"'
```

The coordinator runs the fixed `tutorial-catalog-check`,
`tutorial-manifest-check`, and `tutorial-check` commands, then delegates the
exact ordered chapter set to the existing X11 runner. Fresh HOME/XDG/temp
directories and cleared inherited `GENTLE_*`/API-key settings apply to helper
checks too. Offline namespace enforcement is mandatory before workflow checks;
environment flags alone do not prove network isolation. No online or private
study route is enabled by this coordinator.

`gentle.tutorial_acceptance_candidate.v1` binds the full revision, lockfile,
catalog/manifest hashes, three binary hashes and version output, selected
contracts, starter/oracle workflow hashes and declared input-file hashes.
Candidate identity is checked again during and after execution. The report
`gentle.tutorial_acceptance.v1` retains checks, exact command arguments, exit
codes, output bytes/hashes, required chapter results, and the underlying GUI
report/ledger references. Missing tools, missing automation support, interruption,
failed scientific checks and unexecuted chapters cannot become a green result.
The existing step timeout policies are unchanged. Termination unwinds the GUI
runner so its separate application process is stopped and receipts survive.

The final evidence check requires every selected chapter and step, the correct
starter/oracle distinction, typed verifier outcomes, saved project/oracle
hashes, and the candidate/contract/snapshot/raw-image bindings. A passing run
derives `selection.json` from the retained declared screenshot checkpoints;
no manual revision update is necessary. Recheck retained evidence offline:

```bash
python3 scripts/tutorial_acceptance.py verify \
  --candidate "$CANDIDATE_SHA" \
  --candidate-binding "$ACCEPTANCE_DIR/candidate.json" \
  --evidence-dir "$ACCEPTANCE_DIR/gui"
```

This is an integrity/binding check of externally produced receipts, not a
replacement for Glen observing the tutorial, inspecting screenshots or deciding
whether its language is understandable. Cryptographic hashes do not certify
that an untrusted producer told the truth. Retain the candidate JSON/hash with
the run report in the auditor's evidence archive.

Screenshot staging is a separate explicit command after review:

```bash
python3 scripts/publish_tutorial_gui_screenshots.py \
  --evidence-root "$ACCEPTANCE_DIR/gui" \
  --selection "$ACCEPTANCE_DIR/selection.json" \
  --candidate-binding "$ACCEPTANCE_DIR/candidate.json" \
  --expected-revision "$CANDIDATE_SHA" \
  --output-root "$ACCEPTANCE_DIR/review-images"
```

Strict staging revalidates the full run, accepts only verified checkpoints,
and refuses an existing output directory. A human may narrow the selection
for teaching, never add unverified steps. Nothing updates repository screenshots
or sends private data anywhere automatically. Review/sanitize captures before
public use; even a synthetic GUI can reveal host/display details.

The publisher's existing `--check` remains **archive integrity only**. Historical
teaching images may remain explicitly pinned to an older revision; they are
not current-candidate proof. Do not compare their revision to moving HEAD:
committing screenshots itself changes HEAD. New candidate acceptance must use
new receipts for the explicitly frozen candidate.

For a repair, retain the failed directory, fix on a development branch, freeze
and rebuild the next candidate, and rerun into a new directory. Optional
`--supersedes /path/to/prior/tutorial-acceptance-report.json` retains the prior
report hash; it never imports old passes or changes approval authority. Private
resource/biological tests remain separately authorized on copied state, with
their own data/index fingerprints and verdict. A green offline tutorial gate
must not be reported as a completed private specificity study.

Offline regression tests (synthetic receipts, no X11 capture):

```bash
python3 -m unittest scripts.test_tutorial_acceptance \
  scripts.test_tutorial_gui_acceptance scripts.test_publish_tutorial_gui_screenshots
```

## 6. Practical implementation order

Release packaging has a separate [build-only candidate path](release.md#build-only-candidate-verification).
Installer and container workflows accept the same immutable candidate SHA,
default to `publish=false`, and retain revision/lockfile-bound evidence without
publishing a release. A green ordinary CI run is not a substitute for the full
locked-workspace and external GUI/scientific gates on that exact candidate.
Offline policy regressions use synthetic temporary Git repositories and fake
installer bytes; they do not claim that Docker or native packaging has passed:

```bash
python3 -m unittest scripts.test_release_candidate -v
```

1. Keep extending engine tests alongside new operations
2. Add CLI integration tests for current protocol
3. Implement deterministic SVG exporter from view model
4. Add snapshot suite and CI wiring

## 7. Notes

- PNG snapshot testing is possible but less stable than SVG/text snapshots.
- Prefer view-model-to-SVG tests over full GUI pixel tests for determinism.
- Keep one small GUI smoke test suite for wiring confidence only.
