# Glen: Next Verification And Return

## Carrier Follow-Up, 2026-10-20

After the next ordinary owner push, verify the new carrier changes without
repeating completed sequence-design/expression/gel campaigns. Record the source
SHA, toolchain/profile, lock/binary/input hashes, commands, pass/fail/not-run
and artifact locations. No new interim freeze, publication or scientific approval.
Use a clean checkout and explicit run-owned `CARGO_TARGET_DIR`.

```sh
cargo test --locked -p gentle-protocol labware
cargo test --locked -p gentle-render labware
cargo test --locked --lib labware_ -- --test-threads=1
cargo test --locked --lib rack_openscad_all_templates_preserve_bottom_floor_and_open_top -- --test-threads=1
cargo check -q --locked
python3 -m unittest scripts.test_labware_models -v
```

Check both Safe-Lock support modes, exact 4x9 dimensions/A1/D9, refusal reasons,
generic-op/shared-shell/GUI labels, unverified disclosures and unchanged project
state. Verify the canonical parity projection: only template arguments change,
no operation or command rows are added. Retain native GUI tube/carrier export
control acceptance separately; authored GUI tests are not live acceptance.

Capture and commit the nine original carrier baselines from `1c5be317` source-
matched binaries: SCAD, simulation JSON and the embedded physical spec for each
legacy template on one fixed synthetic rack. Document exact recreation and raw
hashes; add scoped LF rules, LF/CRLF coverage in `scripts.test_tutorial_checkouts`
and Rust byte-comparison tests. Compare against new-source raw exports, without
normalizing or regenerating away differences. Glen's earlier 9/9 comparison is
historical evidence, not a committed fixture or this source's pass.

Render both new Safe-Lock carriers in OpenSCAD from real source-matched GENtle
exports. Require finite nonempty meshes, expected bounds, 2-mm floors and open
bores. Retain original SCAD/spec/mesh hashes and logs; optional tool skips are
not acceptance. Physical fit, cap clearance and thermal/material suitability
remain unverified even with successful mesh renders.

## Historical Covering Request, 2026-10-19

Updated 2026-10-19: after his push, Steffen should forward the
**[concise parity/expression follow-up](glen_parity_expression_followup_20261019.md)**.
Only that message is active for this request. The execution checklists and
earlier capture requests below are retained background, not renewed authority
or a request to repeat completed checks. No release freeze is selected.

The original covering request was prepared 2026-10-18. Its source-bound return
and historical procedures follow; continue development independently.

## Returned Results And Remaining Request

Glen executed the pre-inspector checks below at
`394dd5d0523ea91c19a8eb94d5ed8ad7bd5299b2`; see the
[verified evidence and six-item status](glen_audit_followups_20261008.md#focused-verification-at-394dd5d-2026-10-18).
README 4/4, sequence-design Rust 67/67 plus actual CLI suite 22/22, starter
suite 6/6, workflows 95/95 and rack/labware 127/127 pass. The nine historical
carrier outputs are byteidentical. Keep those checks closed at that source;
do not reissue this entire original checklist as unfinished work.

One reported defect is repaired locally: the parity matrix's two missing
tube-operation rows/counts. Codex applies Glen's canonical generated delta and
preserves newer inspector surfacing without running local Rust. On the next
ordinary owner-pushed revision, run `parity-matrix-check` with its source-matched
docs binary, `cargo test --locked -j1 --test parity_matrix_freshness`, and inspect
normal CI. Run the seven inspector regressions and `gentle-shell` UI-intent
checks in original item 3: `5893ebb0` postdates the tested source. Its updated
09.03 tutorial/catalog projection also needs source-matched freshness checking.
No new bundle or intermediate freeze is requested; final-candidate affected
checks remain on the normal release path. OpenSCAD mesh execution remains
`not_run` on Glen's host, not a failed Rust exporter. Native, agent, scroll and
package boundaries below receive no new acceptance from this return.

## Source And Isolation

Fetch published `smoe/gentle_rs main` once and record the resolved full SHA and
tree. Test that source in a clean detached checkout with an explicitly set,
run-owned absolute `CARGO_TARGET_DIR` outside the checkout. Record toolchain,
platform, features/profile and raw lockfile SHA-256. Use only binaries built
from this checkout; preserve the primary checkout, shared caches and historical
evidence. A later main commit does not invalidate this scoped report or require
another synchronization request, but it is not automatically accepted by it.

No push, workflow dispatch, tag, release publication, provider contact, live model
or package/container build is requested. The owner's follow-up additionally
requests only the four synthetic inspector tutorial screenshots described below;
no unrelated GUI/desktop/private-data capture is authorized. Existing normal CI
can be inspected without restarting it. Keep `package-opt1-fallback` separate;
the selected .13 packaging recipe remains opt2. Codex runs no local Rust.

## Original Focused Checks And Inspector Follow-Up

1. **README metadata repair:** verify the `57a2403a` repair is present. Run the
   release-version regression below and inspect the new normal push CI when
   available. Historical run `39066178748` at `1c5be317` failed the README
   label guard on Linux and Windows; it did not test the repair. Both root
   libraries and all five stripped-child families passed at that older SHA.
   Preserve those bounded results rather than calling the failed run green.

   ```sh
   cargo test --locked -j1 --test release_version_consistency
   ```

2. **Sequence-design follow-up:** verify `8a0e2670` and `958b2238` are included.
   The first retains a freshly validated best after edited-candidate match-cap
   overflow, marking search incomplete without optimality/infeasibility proofs
   instead of returning Invalid and discarding it. The second adds independent
   180 single-symbol and 10,800 paired-symbol IUPAC cases without changing the
   production matcher. Eighteen offline oracle/harness passes are not Rust
   acceptance; execution subsequently passed at `394dd5d` as recorded above.
   Execute the core/protocol tests that the failed older workspace run did not
   reach, plus the shared engine and CLI boundaries:

   ```sh
   cargo test --locked -p gentle-sequence-design -j1
   cargo test --locked -p gentle-protocol -j1 sequence_design
   cargo test --locked --lib --no-default-features -j1 sequence_design -- --test-threads=1
   cargo test --locked --bin gentle_cli --no-default-features -j1 sequence_design -- --test-threads=1
   cargo check -q --locked -j1
   cargo build --locked --no-default-features -j1 --bin gentle_cli --bin gentle_examples_docs
   ```

   Replay tutorial 06.07 and the four actual CLI acceptance tests using the
   [functional procedure](sequence_design_acceptance.md#exact-candidate-functional-checks).
   For these functional checks the dev binaries just built are sufficient;
   record that profile, without claiming package or performance acceptance.
   Set `GENTLE_TUTORIAL_BIN_DIR` to the absolute run-owned target's `debug`
   directory, `GENTLE_SEQUENCE_DESIGN_CANDIDATE` to the recorded full HEAD, and
   `GENTLE_SEQUENCE_DESIGN_EVIDENCE_DIR` to a fresh absolute evidence directory:

   ```sh
   python3 -m unittest scripts.test_sequence_design_audit_oracle scripts.test_sequence_design_tutorial -v
   ```

   All four actual CLI tests must execute; skips and the four Rust CLI unit
   tests are not substitutes. Check retained-best behavior after match-cap
   overflow, incomplete status without false proof flags, cancellation clearing
   output, lexical ties and independent IUPAC coordinates/strands/palindromes.
   Preserve both strategies' joint motif/global/window-GC, tampering,
   materialization, save/reopen and same-session undo/redo checks. Historical
   `de9df8af` prototype acceptance stays separately attributed.

3. **Integrated .13 expression/storage/tutorial paths:** follow the focused
   commands, generator freshness/parity checks and two real-CLI starter replays
   in [the preparation handoff](internal_13_preparation_handoff.md#external-completion-before-candidate-selection).
   Its owner dispatch/publication section is not part of this request. Use the
   same freshly built binaries with `scripts.test_release_workflow_tutorials`;
   no optional-binary skip is acceptance. Retain LF/CRLF and quoted-path results,
   deterministic container order, conservative expression refusals and unchanged
   review state/files. Report drift; do not regenerate historical evidence.

   The merged `gentle_rs_3` inspector at `5893ebb03b4915ba1af2f77f1d7681904266b6f9`
   supersedes the older preparation note deferring all inspector work. Include:

   ```sh
   cargo test --locked --lib -j1 protein_expression_inspector_ -- --test-threads=1
   cargo test --locked -p gentle-shell -j1 ui_intent
   cargo test --locked --lib -j1 mcp_server::protein_expression::tests:: -- --test-threads=1
   cargo test --locked --lib -j1 protein_expression_inner_agent_suggestion_ -- --test-threads=1
   cargo test --locked --lib -j1 protein_expression_agent_card_ -- --test-threads=1
   cargo test --locked --lib --features gui-test-support -j1 protein_expression_capture_targets_ -- --test-threads=1
   ```

   Require all seven inspector regressions to execute. Verify shared-report
   parity, explicit/absent inputs, typed withholding, copy-only suggestions,
   stale/closed-worker rejection and unchanged project/history/files. Check
   canonical tutorial/catalog/parity freshness at this source; static metadata
   edits are not generator acceptance. Native window review stays at release
   preparation. The later owner request additionally asks for the four
   [09.03 screenshot checkpoints](tutorial/09-03_protein_expression_handoff_cli.md#required-tutorial-screenshots)
   from a verified binary and clean synthetic profile. Four manual macOS PNGs
   at `9f01098b` are now retained with binary/input/image hashes and limits in
   the linked archive. They do not contain the uncompiled parity follow-up.
   Review them; retain typed GUI reports and semantic snapshots at the integrated
   candidate separately, without treating these images as automated acceptance.
   No provider/model contact, unrelated capture or automatic capture is allowed.
   Also execute the named MCP `protein_expression_handoff` discovery/report/
   refusal regressions, the deterministic inner-agent command comparison and
   the `gui-test-support` status/locator test. The new control card is guidance,
   not evidence of a model's behavior; compare actual reports on saved state.

4. **Merged rack safety and tube reference models:** verify the integrated
   `gentle_rs_2_main` commit `5b1da7d3981e7a2a4280cadc4b86b9345d82896c`.
   Run its protocol/renderer, headless engine/shell and default GUI regressions,
   then existing rack export and capability/parity checks:

   ```sh
   cargo test --locked -p gentle-protocol -j1 labware::tests::
   cargo test --locked -p gentle-render -j1 labware::tests::
   cargo test --locked --lib --no-default-features -j1 labware_ -- --test-threads=1
   cargo test --locked --lib -j1 labware_ -- --test-threads=1
   cargo test --locked --lib -j1 rack_ -- --test-threads=1
   python3 -m unittest scripts.test_labware_models -v
   ```

   Require positive, overflow-checked admission through 1,536 slots, refusal
   above the limit, unchanged unsupported saved racks on reopen, and refusal
   before overwriting any existing export. Exercise generic `op` inspection
   and all three tube exports with quoted paths: deterministic SCAD, legacy
   result decoding, explicit artifact-write confirmation and unchanged project
   state/history. Generator/parity checks must include the new operation IDs.
   OpenSCAD mesh checks alone do not execute the Rust exporter; missing
   OpenSCAD is a recorded skip, not mesh acceptance. Keep source URL/digest/
   revision/page and exact/approximate/unknown dimensions intact.

   Follow [the provenance and legacy-baseline boundary](../data/labware/README.md#verification-boundary).
   In separate isolated checkouts, compare raw SCAD, simulation JSON and spec
   JSON for the same fixed rack and all three existing carrier templates,
   before this feature at `1c5be3175323eec41359b980b5a41e22c7edaec1` and at
   the tested integrated source. Retain input and before/after hashes outside
   Git; no source-only or template-only result proves byte identity. Native
   Rack controls remain a release-stage check, not new capture authority here.
   Tube models are body/collar visualizations with cap/hinge/cavity omitted,
   not printable vessels or fit evidence. New carriers, support/pitch,
   96-position interpretation and storage boxes remain deferred owner choices.

Use `--locked` throughout. If dependencies are already cached, `--offline` may
also be used; otherwise record retrieval failures as environmental. Do not
change the lockfile or widen resource limits silently to obtain a pass. Missing
tools, zero selected tests, timeouts and unrun platforms remain explicit.

## Remaining Six-Gate Boundaries

Use the [six-item ledger](glen_audit_followups_20261008.md) as the status owner,
and the [release-stage request](glen_final_candidate_execution_handoff.md) only
when release preparation is authorized. Audit preservation, 08.04 contracts,
live-agent acceptance, honest visual evidence, scroll benchmarks and package
acceptance retain their individual evidence boundaries. These focused checks
do not close all six. Existing tutorial/link repairs and historical audits
must not be reissued as unfinished work. Native/physical-host, model, scroll,
package and scientific verdicts stay separate; this request authorizes none
of their additional execution.

## How To Return Results

Send **one concise response to Steffen for forwarding into the main Codex
chat**. Do not rely on direct Glen-to-Codex messaging or assume a path on one
host exists on the other. Include:

- Tested full SHA/tree, platform, Rust/Cargo versions, features/profile,
  lockfile and binary hashes, and any same-source CI URL/conclusion.
- Per-check `pass`, `fail`, `not_run` or `blocked`, actual selected test counts,
  skips, exact commands/exit codes and primary failure versus cascades.
- Remaining six-gate boundaries, distinguishing new results from historical
  passes and software execution from native/scientific/release acceptance.
- Evidence location and SHA-256 of its manifest, plus any concrete defect and
  minimal proposed repair. Do not modify or push source as part of this check.

Supply a small folder/archive containing `REPORT.md`, command/exit-code logs,
raw stdout/stderr, CI/CLI receipts and a per-file SHA-256 manifest. Steffen can
copy it to `/tmp/glen-verification-<short-sha>-<date>/` on the Codex host and
relay the path and manifest hash here. Keep bulky binaries, raw captures and
reference inputs on Glen's host with hashes/locations; transfer them only on
request. Never include credentials or private projects. Record content hashes
as identity checks, not authenticated provenance or human scientific approval.

## Preserved Prompt Records: Not Active Instructions

The following received files are committed verbatim for history:

- [Release-link repair prompt](release_handoff_link_repair_prompt.md): repair
  already closed; its preparation-time failure and old runner/version status
  are historical, not a request to restore the deleted .12 planning ledger.
- [Sequence-design closure prompt](sequence_design_closure_codex_prompt_20261014.md):
  closure was reconciled in the plan/acceptance records. The current
  [compiler decision](decisions.md#workspace-compiler-compatibility) supersedes
  its broader proposed workspace minimum/publication gate; select no MSRV.
  The [issued handoff's SHA correction](glen_sequence_design_acceptance_handoff_20261014.md)
  remains intact.
- [Original Claude VKORC1 prompt](vkorc1_acceptance_codex_prompt_20261014.original-claude.md):
  retain the original separately from the [reconciled task](vkorc1_acceptance_codex_prompt_20261014.md).
  Its in-progress CI, undecided .13 aim and other snapshot claims are not
  current status. `--offline` is not a substitute for `--locked`, nor proof
  that a lockfile was rewritten. No old embedded dispatch is authorized here.
