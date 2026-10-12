# Glen Execution Handoff: Separate Frozen Audits

The owner requests two independently reported audit targets:

- `e981469e4beb637d2def12617ff974c76661ace5`: the six follow-up gates below.
- `dfadf2a554dd8a23afae923592baa65d37162bc9`: the supplemental GenBank/variant
  regressions and experimental sequence-design recheck below.

Never substitute moving `main`, mix results between these targets or infer
merged-main/release acceptance. Read this handoff at its documentation commit;
run each audit from its own clean, isolated frozen checkout.

## Frozen Source And Bundle

Owner request, 2026-10-12: execute the remaining audit gates at
`e981469e4beb637d2def12617ff974c76661ace5`, not moving `main`.
The owner will push the committed instructions; Codex performs no push.
All required reading below is tracked in this repository. Local Codex logs,
temporary verifiers and older binaries are not prerequisites.

Expected source tree: `05246928299ebe32a6008cf38e391c1817b4ee2a`.
Expected raw `Cargo.lock` SHA-256:
`ead4421545a125267e33af84aeeea39f0e7ea2f60717a59ef0271d0e71855f1b`.

The supplied `gentle-glen-e981469e.bundle` is **46,038 bytes**, SHA-256
`8bd274c6ee4d811401100dcd584ec5bdafd5572e4e5e0e5875bbcea3771402bd`.
Its prerequisite is `f56fdf849957d03d6f232c041cc74347d5ab471b`.
Before importing it, compare both the byte count and hash. On Linux use
`sha256sum`; on macOS use `shasum -a 256`.

Use absolute paths and a new checkout, never reset Glen's primary tree:

```bash
BUNDLE=/absolute/path/gentle-glen-e981469e.bundle
CHECKOUT=/absolute/path/new/gentle-e981469e
sha256sum "$BUNDLE"
wc -c < "$BUNDLE"
git init "$CHECKOUT"
git -C "$CHECKOUT" fetch --depth=1 https://github.com/smoe/gentle_rs.git \
  f56fdf849957d03d6f232c041cc74347d5ab471b
git -C "$CHECKOUT" bundle verify "$BUNDLE"
git -C "$CHECKOUT" fetch "$BUNDLE" HEAD:refs/heads/glen-audit-source
git -C "$CHECKOUT" checkout --detach e981469e4beb637d2def12617ff974c76661ace5
git -C "$CHECKOUT" rev-parse HEAD HEAD^{tree}
git -C "$CHECKOUT" status --porcelain
sha256sum "$CHECKOUT/Cargo.lock"
```

Require the exact HEAD/tree/hash above and empty status before building. The
bundle was independently reconstructed and verified in a clean checkout.
After the owner's push, the bundle is optional: a new checkout can fetch and
detach at the exact `e981469e` SHA directly from GitHub instead. Read this
handoff at its pushed documentation commit without changing the tested checkout.
Keep new evidence and reference caches outside it; preserve historical bytes.

The later `gentle_rs_3` merge is
`5b5998086a4f98bbef17502f1191367c3f8e4f80`. It adds a production GenBank/SPDI
round-trip fix and is **not** contained in this bundle. An `e981469e` pass does
not certify that merge, current main or a release. Record failures at the
requested SHA; do not silently transplant the fix or combine different SHAs'
passes. A later release candidate requires its own complete same-source gates.

## Six Follow-Up Gates: e981469e

Run Rust only on Glen's approved external build hosts. Use the frozen
[CI workflow](../.github/workflows/ci.yml) for platform dependencies, features,
isolation and commands. Report unavailable platforms/tools/authentication as
unrun; do not hide failures, waive skips or substitute older executables.
A push-triggered CI run on newer main is not an `e981469e` run. CI's dispatched
ref must resolve to the requested source if used for these gates. Package
workflows instead have an explicit `candidate_sha` input; record their actual
workflow revision separately and verify the candidate checkout in the receipt.

1. **Source, Rust and adapter gates.** Start with the failed boundaries below,
   then default `cargo check -q --locked`, the full lib/integration suites and
   the CI matrix's headless CLI/MCP, JavaScript, Lua and all-features checks.
   Require fresh Linux/macOS/Windows results, nonzero targeted test counts,
   formatting and catalog/tutorial/parity freshness. Keep primary failures
   separate from build or mutex-poisoning cascades.

   ```bash
   cargo test -q --locked --lib -j1 engine::primer_specificity_multi::tests:: -- --test-threads=1
   cargo test -q --locked --lib -j1 engine_shell::tests::glossary_ -- --test-threads=1
   cargo test -q --locked --lib -j1 execute_introspect -- --test-threads=1
   ```

   Rerun the offline helper suite below. Codex's prior 263-test pass is helper
   evidence only, not Rust, GUI, agent or package acceptance.

   ```bash
   env -u GENTLE_TUTORIAL_BIN_DIR python3 -m unittest \
     scripts.test_release_candidate scripts.test_package_desktop \
     scripts.test_release_build_monitor scripts.test_gui_startup_summary \
     scripts.test_container scripts.test_tutorial_checkouts \
     scripts.test_primer_specificity_multi_tutorial scripts.test_ci_vkorc1_gui_audit \
     scripts.test_tutorial_gui_acceptance scripts.test_prepare_public_vkorc1_audit \
     scripts.test_tutorial_acceptance \
     scripts.test_tutorial_walkthroughs.WalkthroughSourceTests \
     scripts.test_tutorial_gap_lessons scripts.test_regulatory_motif_fixture -q
   ```

2. **08.04 native contract.** Reuse CI's `vkorc1-gui-audit` build with
   `desktop-gui,gui-test-support` and all six focused filters, followed by
   [ci_vkorc1_gui_audit.py](../scripts/ci_vkorc1_gui_audit.py) on the approved
   isolated Linux/X11 host. Retain ordinary native input, explicit T,
   multiallelic refusal, completion transition and full-record starter/oracle
   comparison. Synthetic 20-bp input is not public-locus or scientific proof.
3. **Genuine live-agent acceptance.** Prepare fresh public inputs with
   [prepare_public_vkorc1_audit.py](../scripts/prepare_public_vkorc1_audit.py)
   using a clean checkout and exact-source CLI. It invokes no agent. Follow the
   [committed genuine-agent procedure](audits/glen_followups_20261008/candidate_3ac4a02d/README.md#genuine-live-agent)
   as a procedure only: new starter, existing authorized login, clarification
   with no mutation, explicit T, two unexecuted review-first suggestions,
   ambiguous-alternate refusal, then reviewed manual execution. Compare all
   four full ordered `Seq` records with fresh direct/shared-shell oracles.
   Retain new responses and raw-state checks; no credential transfer, nested
   `agents ask` or silent automatic execution.
4. **Honest visual evidence.** Retain unedited native captures and separately
   labelled public base projections. Bind them to persisted inputs/exports and
   independently verify the sole C/T change. Do not recolour screenshots,
   reinterpret identical whole-map images or claim scientific approval.
5. **Scroll-event/next-frame smoke.** Reuse CI's `gui-scroll-audit` procedure:
   [fresh PATZ1 preparation](../scripts/prepare_real_patz1_tutorial.py), actual
   compiler-artifact executable discovery and `gui_operations --test`.
   Require 24 observations: three fixtures, four sizes, both directions,
   unchanged sequence/annotations and bounded stable next-frame viewport.
   Retain dev/profile metadata; CPU determinism is not native latency.
6. **Non-publishing .12 packages.** Follow the committed
   [build-only release instructions](release.md#build-only-candidate-verification)
   with `candidate_sha=e981469e4beb637d2def12617ff974c76661ace5`,
   `publish=false` and `package-opt1`.
   Require the exact lockfile/workflow/recipe, all three native platform packages
   and the headless container. Independently check actual extracted archive
   membership/checksums and execute all five native/three container binary
   identities/entrypoints. Retain build/resource, extraction and smoke receipts.
   A green build or policy test alone is insufficient. Workflow dispatch still
   requires explicit owner authorization; if unavailable, report it as unrun.

Do not push, tag, publish a release/image, alter original `30f23aa4` evidence,
or claim human scientific approval. No Claude consultation is requested.

## Supplemental Audits: dfadf2a5

Freeze both supplemental audits at
`dfadf2a554dd8a23afae923592baa65d37162bc9`, source tree
`9c0e7049ca1cda01c520587d85efed8bb4cb9eef`, raw `Cargo.lock` SHA-256
`ead4421545a125267e33af84aeeea39f0e7ea2f60717a59ef0271d0e71855f1b`.
The earlier `e981469e` bundle cannot supply this commit. Fetch the exact object
into a separate new checkout; if unavailable, request its reachable branch or
Git bundle rather than using nearby `main`:

```bash
CHECKOUT=/absolute/path/new/gentle-dfadf2a5
git init "$CHECKOUT"
git -C "$CHECKOUT" fetch --depth=1 https://github.com/smoe/gentle_rs.git \
  dfadf2a554dd8a23afae923592baa65d37162bc9
git -C "$CHECKOUT" checkout --detach dfadf2a554dd8a23afae923592baa65d37162bc9
git -C "$CHECKOUT" rev-parse HEAD HEAD^{tree}
git -C "$CHECKOUT" status --porcelain
sha256sum "$CHECKOUT/Cargo.lock"
```

Require the exact HEAD/tree/lockfile hash and empty status before building.
Run the following from this checkout on Glen's approved external hosts, with
evidence outside it. Codex ran **no Rust verification** of these changes.

### GenBank And Variant Safety

```bash
cargo test --locked --lib dbsnp_ -- --test-threads=1
cargo test --locked --lib materialize_variant_allele -- --test-threads=1
cargo test --locked --lib reporter_construct_handoff_ -- --test-threads=1
cargo test --locked --lib variant_followup_allele_pair_ -- --test-threads=1
cargo check -q --locked
```

Require nonzero counts for every filter and confirm both new tests execute:

- `dbsnp_spdi_genbank_round_trip_preserves_wrapped_json`: actual GenBank export
  and shared-loader import under LF and CRLF, preserving long alleles, empty
  strings, allele order, coordinates and qualifiers.
- `dbsnp_spdi_unfolding_preserves_valid_json_invalid_data_and_other_qualifiers`:
  valid JSON, escapes, invalid/non-SPDI values and ordinary notes stay unchanged.

Assembly/reference status must not change; re-import must not upgrade evidence
or manufacture VCF alleles. Successful GUI pair materialization intentionally
retains **two undo entries**, unlike a single sequence-design derivation.

Acceptance remaining separately authorized: native Windows round-trip
(Unix-only mock-fetch tests do not certify Windows); actual refSNP response and
checksum, selected placement and prepared GRCh38 reference comparison; tutorial
08.04 native C/T refusal/retry with unchanged source and correct genomic-forward
orientation. Reporter orientation/topology and package/release acceptance also
remain separate. These synthetic regressions do not certify those boundaries.

### Experimental Sequence-Design Recheck

The supplied prior read-only review agrees with the documented
[historical-prompt addendum](dna_sequence_optimization_claude_prompt.md) and
[intentional reverse-translation/design distinction](decisions.md#synthetic-sequence-design-boundaries).
It identifies no further implementation repair; Glen's independent acceptance
remains open. This is not a new Claude consultation or an acceptance verdict.
Use the [active plan](dna_sequence_optimization_plan.md) and
[tutorial 06.07](tutorial/06-07_synthetic_sequence_design.md), at this frozen SHA.

```bash
cargo test -p gentle-sequence-design --locked
cargo test -p gentle-protocol --locked sequence_design
cargo test --locked --lib sequence_design -- --test-threads=1
cargo test --locked --bin gentle_cli sequence_design -- --test-threads=1
GENTLE_TUTORIAL_BIN_DIR=/absolute/path/to/same-SHA/binaries python3 -m unittest scripts.test_sequence_design_tutorial
```

Build or obtain the ordinary `gentle_cli` executable from this exact source
before the Python replay; `cargo test --bin` produces a test harness, not that
CLI. Record its source/profile and binary hash. Missing-variable runtime skips
or stale binaries are not acceptance, even if the Python source checks pass.

Replay 06.07 and compare **both strategies** with independent tiny-space oracles,
not the production matcher/translator reused as an oracle. Audit joint motif,
CDS-wide GC and sliding-window GC constraints, client-rehashed invalid outputs
and fabricated search claims, cancellation/budget outcomes, single-operation
provenance, save/reopen and undo/redo.

The receipt verifies **output constraints**, not search history or minimum-edit
claims. Keep this experimental verdict **outside the `.12` release gate**;
native GUI, model behavior and performance each need separate evidence and
authorization. Do not transfer a design approval to reverse translation.

## Evidence To Return

Return a concise six-gate pass/fail/unrun summary for `e981469e`, plus separate
GenBank/variant and experimental sequence-design summaries for `dfadf2a5`.
Provide a retained bundle or accessible location containing:

- Full tested source SHA/tree, raw lockfile hash, clean-status observation,
  commands/exit codes, test counts and explicit skips/unavailable gates.
- OS/architecture, toolchain and external-tool versions, exact features/profile
  and package recipe; binary `--version` output and SHA-256 for each executable.
- Complete logs, native captures, public inputs/projects/exports, genuine-agent
  requests/responses/state checks, all 24 scroll observations, and native/archive
  and container receipts including actual entrypoint execution.
- A complete per-file SHA-256 manifest plus archive hash/byte count. Bind each
  receipt to its actual source/profile/input hashes; never reuse old verdicts.

Keep the two source trees' logs, binaries and receipts separately identified;
report skips and platform omissions explicitly. Distinguish product defects
from environment/data limitations. Neither supplemental summary closes the
six-gate ledger or accepts a release.

Report failures before attempting source changes; any repair is a separate
commit and candidate with rerun gates. Do not include credentials or tokens.
