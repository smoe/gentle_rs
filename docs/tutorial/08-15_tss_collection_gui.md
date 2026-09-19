# Annotated Starts To Reusable TSS Windows

An offline, synthetic tutorial for a wet-lab biologist and Glen's Linux alpha
acceptance. This teaches annotation geometry and explicit approval, not which
promoter is active. No model service, private study or human genome download is
needed. The [generated chapter](generated/chapters/08-15_tss_collection_gui.md)
contains the executable reference and typed GUI acceptance subset.

## Clean Starter

Use GUI, CLI and helper binaries from one frozen candidate revision. Record
their hashes, toolchain and [fixture](../examples/assets/tss_tutorial/README.md)
hashes. Build `gentle` with `gui-test-support` for semantic snapshots. Use fresh
HOME/XDG paths and caches, clear inherited `GENTLE_*` variables, and disable
network access. Screenshots are opt-in acceptance evidence, not a default side
effect of opening this workspace.

From the repository root:

```sh
target/debug/gentle_examples_docs example-project tss_collection_gui_starter /tmp/tss-starter.json --run-dir /tmp/tss-starter-run
target/debug/gentle_examples_docs example-project tss_collection_gui_oracle /tmp/tss-oracle.json --run-dir /tmp/tss-oracle-run
```

Open **tss-starter.json**, never the oracle. It has exactly one anchored 2,000-bp
locus, `tss_locus`, and no TSS collection. The oracle is an independent engine
replay for comparison only. Optional BLAST index warnings do not affect this
annotation-only exercise and do not constitute specificity evidence.

## Preview And Approve

1. Double-click `tss_locus`. Choose **TFBS scan > Transcript starts / TSS windows**.
2. Enter gene **TOY**, collection **tss_windows**, upstream **500**, downstream
   **200**. Choose **Inspect starts (no changes)**. Opening a view is not approval.
3. Check the three rows below. **Select available**, then explicitly
   **Approve and create selected windows**. This uses the exact preview digest
   and selected IDs; changed parameters or source content require a new preview.
4. **Refresh collections**. Expect **readable / not checked**, not validated.
   Select the row and **Inspect stored collection**. Only this action validates
   persisted members. Copy its JSON into the evidence bundle.

| Annotated start | Transcript records | Genomic window (1-based inclusive) | Strand |
| --- | --- | --- | --- |
| 601 | plus_a, plus_b | 101..801 | + |
| 901 | plus_c | 401..1101 | + |
| 1500 | minus_a | 1300..2000 | - |

Each window is **701 bp**: the start itself adds one base to 500 + 200.
The minus-strand window is reverse-complemented. All starts appear at local
base **501**; genomic coordinates decrease along the minus-strand window.
Two synthetic gene IDs share the query label TOY here; their records are not
merged. Repetitive toy DNA is unsuitable for biological TFBS or primer claims.

Inspect the independent oracle without the GUI:

```sh
target/debug/gentle_cli --project /tmp/tss-oracle.json shell 'promoters tss-list'
target/debug/gentle_cli --project /tmp/tss-oracle.json shell 'promoters tss-collection tss_windows'
```

The list is discovery metadata. The second command returns the validated
`gentle.tss_collection.v1`. Compare its approval/membership fingerprints and
member bases/annotations with the GUI's saved project, not just a screenshot.

## Automated Linux Subset

The `tss_collection_gui` contract has 18 semantic steps: open the form, enter
the gene, preview, select, approve, refresh, inspect, open the members twice,
cancel then confirm forgetting, undo and inspect again. The external runner
uses ordinary X11 input to locate controls. It never executes prose or arbitrary
shell verifiers. A fixed `GetTssCollection` verifier rejects stale members and
checks starts, strands, shared membership and sequence content against the oracle.
The window-set check expects exactly four subject-bound DNA viewers: the source
locus and three members. It detects missing or duplicate viewers, but is not a
substitute for the typed report and sequence checks. Metadata-only changes are
saved before on-disk verifiers run. Forgetting is checked as absence of the
registry entry together with retention of every sequence and viewer. Undo is
followed by fresh collection validation, not merely a restored label.

Inside a network-isolated Linux/Xvfb session with an EWMH window manager,
`xdotool`, `xdpyinfo`, `xprop`, `xwininfo` and `scrot`:

```sh
python3 scripts/tutorial_gui_acceptance.py --chapter tss_collection_gui --profile offline-core --network-enforcement external_policy --evidence-dir /tmp/tss-gui-acceptance
```

Use `external_policy` only when isolation actually is enforced externally;
record the container/network-namespace policy. Keep the checkpoint ledger,
semantic snapshots, raw screenshots, typed reports, project hashes and binary
identities. Missing/ambiguous controls are a `harness_gap`, not guessed clicks.
Authoring and unit-testing the contract does not establish a Linux GUI pass.

## Manual Lifecycle Checks

These are **not** part of the 18-step automated verdict. Use a copy of the
synthetic project and retain a checkpoint/result for each:

1. Inspect the actual sequence display for 701-bp sequences and the expected
   minus-strand orientation. The automated window count does not certify pixels
   or every coordinate label. Retain screenshots and the collection JSON.
2. Save, close and reopen. Refresh and inspect. Compare approval/membership
   fingerprints and member content with the saved receipt. Whole-file hashes
   may change with display state; record them without calling that a biological edit.
3. In one member's feature editor, move the **Annotated TSS candidate** marker
   from local base 501 to 502 and apply it. This is deliberate damage to a copy,
   not a proposed biological correction. Collection inspection must now reject the
   edited member. Listing still shows the entry as **not checked**; no scan
   may silently skip it. Retain the error rather than an invented empty result.
4. **Forget registry entry...** must initially remove nothing. Cancel once;
   request again and confirm the named ID. Only registry metadata disappears:
   all four sequences, lineage and open windows remain.
5. **Undo** restores the registry entry. Validation must still reject the edited
   member: undoing forget is not undoing the earlier biological edit. Retain
   before/after state and validation receipts.

The automated Forget/Undo cycle uses an intact collection. Steps 3-5 above test
the separate stale-member case; restoring its registry must not repair or hide
the deliberate edit. Engine regressions cover this distinction, but GUI editing
and application restart still need their own retained live evidence.

Changing collection ID cancels a pending forget confirmation. Forgetting does
not authorize sequence overwrite; a new derivation normally needs a new name.

Record separate verdicts for the automated subset, manual lifecycle and live
biology. Glen's TP73/DeltaNp73 acceptance remains a separate private-data test.
None of this fixture proves active promoters or complete biological coverage.
