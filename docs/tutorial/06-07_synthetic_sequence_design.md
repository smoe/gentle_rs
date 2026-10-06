# Redesign a synthetic coding insert: remove a motif and control GC

Can you change a coding DNA sequence without changing its protein, while removing
an unwanted recognition sequence and controlling both overall and local GC?
This ten-minute exercise makes every change small enough to check by eye.

The input is a **hand-crafted synthetic example**, not TP73, PATZ1 or another
natural gene. Its 12 bases encode MFK (methionine, phenylalanine, lysine), followed
by a stop codon. The unusual GC limits are teaching values, **not recommended
thresholds for expression or synthesis**. This experimental feature is outside
the `.12` release gate and does not claim full DNA Chisel compatibility.

## Before you start

Use a new, disposable GENtle project or a new working folder for the CLI.
No genome, JASPAR database, network access or inner agent is needed. The current
GUI route is **Shell / Agent Assistant**, not a dedicated optimization editor.

Open the [committed request](inputs/synthetic_sequence_design.json) and place a
working copy called `request.json` in your tutorial folder. Run the commands
from that folder, or replace file arguments with quoted absolute paths.
For GUI-hosted window-opening commands, use **Agent Assistant**; the DNA viewer's
sequence-local Shell and a headless CLI return an unapplied UI intent instead.
Keep one GUI Shell session open for the later undo/redo exercise.

## 1. Understand the input and the constraints

```text
Input:     ATG TTT AAA TAA
Protein:    M   F   K  stop
Expected:  ATG TTT AAG TAA
```

The synonymous choices here are `TTT/TTC` for F and `AAA/AAG` for K.
GENtle freezes the ATG start and the original stop codon, as well as any flanks
or explicitly protected bases. The CDS is `[0,12)`: zero-based, including base 0
but excluding base 12. The supplied protein `MFK` excludes the stop.

The unwanted motif `TTTAAA` crosses the F/K codon boundary at `[3,9)`.
It is palindromic, so scanning both strands reports this physical site once,
with strand `both`. These are literal DNA patterns, not enzyme-name lookups.

The request below imposes **three joint constraints**:

- Remove `TTTAAA` on both strands across the supplied sequence.
- Require exactly two G/C bases in the whole 12-base CDS.
- Require exactly one G/C base in each complete 6-base CDS window, moving one
  base at a time. There are seven windows, `[0,6)` through `[6,12)`.

GC bounds are inclusive integer basis points: 10000 means 100%.
For this deliberately tiny example, 1666..1667 basis points rounds to 2/12
globally and 1/6 locally. GENtle reports the integer counts, not an approximate
floating-point comparison. Frozen start/stop bases count; flanks and shortened
edge windows do not. All requested constraints must pass together.

Your `request.json` should contain exactly:

```json
{
  "schema": "gentle.dna_sequence_design_request.v1",
  "target": {"kind": "inline_sequence", "sequence": "ATGTTTAAATAA"},
  "purpose": "synthetic_coding_insert",
  "cds": {"start_0based": 0, "end_0based_exclusive": 12},
  "protein_sequence": "MFK",
  "genetic_code": 1,
  "protected_intervals": [],
  "avoid_motifs": [{"pattern": "TTTAAA", "strand": "both"}],
  "gc_content": {"min_basis_points": 1666, "max_basis_points": 1667},
  "gc_window": {"window_bp": 6, "min_basis_points": 1666, "max_basis_points": 1667},
  "search_strategy": "conflict_directed",
  "max_evaluations": 4096,
  "output_seq_id": "synthetic_mfk_redesign"
}
```

## 2. Preview without creating DNA

In GUI Shell or Agent Assistant's command field:

```text
sequence-design plan @request.json --path preview.json
```

CLI equivalent, using one explicit disposable project path throughout:

```sh
gentle_cli --state tutorial.gentle.json sequence-design plan @request.json --path preview.json
```

This writes `preview.json` but creates no project sequence or container.
The input is inline DNA in the request, not an existing sequence to overwrite.
Planning also leaves any already stored project DNA unchanged.

**Review before continuing.** Look for:

| Field | Expected observation |
| --- | --- |
| `status` | `feasible` |
| `source_sequence` | `ATGTTTAAATAA` |
| `output_sequence` | `ATGTTTAAGTAA`, still MFK followed by the frozen stop |
| `edits` | One A-to-G substitution at zero-based position 8 (base 9 in a 1-based viewer) |
| `initial_matches` | One `TTTAAA` site at `[3,9)`, strand `both` |
| `gc_content.input/output.gc_bases` | 1 before, 2 after; integer bounds both 2 |
| `gc_window.input` | Counts `1,1,1,0,0,0,0` over the seven windows |
| `gc_window.output` | Seven exact intervals, each with count 1 and `satisfies_bounds=true` |
| `optimization_complete` / `minimum_edits_proven` | Both true for this tiny fully resolved search |
| `approval_digest` | A real content-bound digest, needed only after review |

Keep search claims separate from DNA validation. A budget stop is unresolved,
not proof of infeasibility. Cancelled or unresolved reports have no applicable
output. Missing output measurements are absent, not zeros.

## 3. Approve and create the exact output

Only if you accept the exact DNA, edits, constraints and limitations, copy the
**actual `approval_digest` from your own preview** into the following command.
`ACTUAL_REVIEWED_DIGEST` is a placeholder, not a digest to submit.

```text
sequence-design apply @preview.json --approve ACTUAL_REVIEWED_DIGEST
```

CLI equivalent:

```sh
gentle_cli --state tutorial.gentle.json sequence-design apply @preview.json --approve ACTUAL_REVIEWED_DIGEST
```

Do not pipe a digest automatically from planning into application. The pause
is the review step; requesting a preview is not permission to create DNA.
Changing the request, strategy or output requires a new preview and approval.

The result creates `synthetic_mfk_redesign` and one singleton container named
`Synthetic sequence design`. Inline input has synthetic-import provenance;
it does not invent a parent sequence or inherit unrelated source annotations.
Its receipt must say `output_constraints_verified=true` and
`search_claims_verified=false`: apply rechecks the output, but does not rerun
search or certify the preview's minimum-edit claim. Saved metadata keeps the
submitted report tagged `unverified_portable_preview`, separate from the receipt.

In GUI-hosted Agent Assistant, open the new DNA window:

```text
ui open sequence-window synthetic_mfk_redesign
```

Verify the length is 12 bp and the bases read `ATGTTTAAGTAA`.
Headless CLI use creates the same project DNA, not a native window. You can open
`tutorial.gentle.json` in the GUI for inspection, but reopening does not restore
the previous process's undo stack.

## 4. Undo and redo in the same GUI session

If you applied through the same still-open GUI Shell/Agent Assistant session,
inspect history and explicitly undo:

```text
history status
history undo
```

The new sequence, its container and design metadata should disappear together.
The preview file stays on disk: undo is not a file-deletion operation.
Explicitly redo to restore the exact product and its provenance:

```text
history redo
```

**CLI limitation:** each `gentle_cli` invocation starts a new engine session.
Saving project state does not save undo checkpoints. A separate
`gentle_cli shell 'history undo'` cannot undo an earlier process's apply.
Use the uninterrupted GUI session for this exercise. The deterministic
same-session regression below tests these commands, without claiming live GUI
acceptance. Do not substitute deletion for undo.

## 5. Why local GC adds a real requirement

All four possible sequences retain MFK and the original ATG/TAA:

| DNA | Total G/C | Motif absent? | Whole-CDS GC passes? | Every 6-base window passes? |
| --- | --- | --- | --- | --- |
| `ATGTTTAAATAA` | 1 | no | no | no |
| `ATGTTCAAATAA` | 2 | yes | yes | no |
| `ATGTTTAAGTAA` | 2 | yes | yes | yes |
| `ATGTTCAAGTAA` | 3 | yes | no | no |

The second row demonstrates the trap: a satisfactory average does not imply
satisfactory local windows. The third row is the only sequence satisfying all
three requests. This table is independently checked in the fast tutorial tests.

As an optional second preview, omit `search_strategy` to select enumeration.
It should return the same tiny-case DNA and edit, not the same search history
or approval digest. Keep previews in separate files, review them separately,
and do not reuse an existing output ID for another apply. Neither the teaching
table nor this small example establishes runtime performance on long inserts.

## Inner-agent route

Ask GENtle's Agent Assistant:

> This is an explicitly synthetic ATGTTTAAATAA coding insert for MFK, CDS
> [0,12), genetic code 1. Preview removal of TTTAAA on both strands, GC bounds
> 1666..1667 basis points over the complete CDS and over every full 6-base
> window, using conflict_directed, budget 4096 and output ID
> synthetic_mfk_redesign. Show the actual DNA, edits, all integer counts and
> search status. Do not apply yet.

The agent should propose the same parser-valid `sequence-design plan` command.
It must not invent output DNA or hashes, infer a CDS from a label, recode a
natural assay template, or silently select another strategy. No agent is needed
to replay this tutorial. In Agent Assistant, **Use reviewed result in next
prompt** shares result content only with your consent; it is not application
approval. Execution receipts alone do not give the model the complete report.
Approve a subsequent exact apply command separately. Window opening and
undo/redo also require explicit review, not `auto` execution.

## Executable checks and independent review

Fast, standard-library-only source/teaching checks:

```sh
python3 -m unittest scripts.test_sequence_design_tutorial
```

To exercise the actual CLI, point at already built current binaries:

```sh
GENTLE_TUTORIAL_BIN_DIR=/absolute/path/to/binaries python3 -m unittest scripts.test_sequence_design_tutorial
```

Without that variable the CLI checks explicitly skip, not pass. With it set,
missing binaries or failed commands are failures. Tests build nothing, use no
network/model and write only temporary projects/reports. For the shared-parser
apply/undo/redo check in one engine session:

```sh
cargo test --lib sequence_design_tutorial --locked --offline -j 1 -- --test-threads=1
```

Glen: audit one frozen source SHA with binary/input hashes, toolchain, features
and effective profile recorded. Repeat both strategies, inspect all seven
windows and the separate receipt claims, then capture the reviewed preview and
the created DNA in the native GUI. Check same-session undo/redo. Preserve raw
captures and semantic evidence; report CLI, GUI, scientific usefulness and
performance verdicts separately. Native/platform, inner-model and human
scientific acceptance are pending until independently performed.

**Success means** that this explicit synthetic request produces the exact
approved sequence with jointly satisfied constraints and coherent provenance.
It does not establish expression, folding, splicing, synthesis, ordering or
experimental suitability. For your own insert, supply its real DNA, protein,
CDS, protected bases and justified thresholds; never borrow the toy's limits.

If a step fails, retain its command, report/status, exact revision and interface
in [tutorial feedback](../../.github/ISSUE_TEMPLATE/tutorial-confusion.md).
For missing files use quoted absolute paths; for `invalid` inspect admission
diagnostics; for `search_exhausted` do not claim infeasibility; for an existing
output ID choose a new ID and obtain a new preview and approval.
