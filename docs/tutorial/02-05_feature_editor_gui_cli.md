# Curate Annotations with Preview, Apply and Undo

> Type: hand-written GUI/CLI walkthrough. Offline; no reference downloads.
> Audience: biologists correcting annotations on an existing DNA sequence.
> Last updated: 2026-09-10. Live GUI acceptance: not yet recorded.

An annotation says what a stretch of DNA represents. Changing its boundaries
does **not** change the bases, establish expression, or repair every annotation
that refers to the same gene. This exercise makes that distinction visible.

Allow about 15 minutes. Use the tiny synthetic input, not a working project.
English control labels below are the labels in the current Feature Editor.

## 1. Open the Starter, Not a Completed Result

Open [feature_editor_demo.gb](./inputs/feature_editor_demo.gb) with
**File -> Open Sequence...**, then open its DNA viewer. Save a new tutorial
project outside the checkout. The sequence is 120 bp, with three invented
annotations in this stored order:

| Feature index | Annotation | Viewer positions | Length |
| --- | --- | --- | --- |
| 0 | gene `DEMO` | 11..100, forward | 90 bp |
| 1 | `misc_feature`, label `segment_A`, gene `DEMO` | 21..40, forward | 20 bp |
| 2 | `misc_feature`, label `reverse_control`, gene `CONTROL` | 61..80, reverse | 20 bp |

The indices are **zero-based table positions**, not permanent biological IDs.
The query result's `feature_id` identifies this stored index; its default
display order is by position. Use `--sort feature_id` when checking stored order,
and never take a row's screen position as its identity.
Coordinates entered in this editor are **1-based, inclusive**. Thus 21..40
means 20 bases and corresponds to the machine interval `[20, 40)`.
For the reverse control, base 80 is its 5-prime end and 61 its 3-prime end.
Neither a CDS nor a real organism is implied by these invented annotations.

Open **Feature Editor** from the command palette or a feature's context menu.
Select `1: misc_feature (segment_A)` in the Location tab.

## 2. Location: Review a Boundary Change

Set **Start (1-based)** to `21` and **End (1-based, inclusive)** to `45`.
Click **Preview**. The before/after interval must be 21..40 -> 21..45, still
forward. The map and saved project must not change merely because you preview.

Now change the end field to `46`. **Apply** must become disabled: the old
preview is no longer approval for this request. Return to `45`, preview again,
then click **Apply**. Only `segment_A` should extend, to 25 bp; `DEMO` remains
11..100, the reverse control remains 61..80, and the DNA remains 120 bp.

Use **Edit -> Undo**, check 21..40, then **Redo**, check 21..45. Continue from
the redone state. Undo/redo belongs to this live project session; do not expect
a new command-line process to inherit that session's undo stack.

For real annotations, review every listed boundary-sharing feature separately.
GENtle does not propagate a change to mRNA/CDS/exon records. A CDS length
warning is a reason to inspect the biological evidence, not permission to
silently adjust its reading frame.

## 3. Create: Overlap Is a Review Prompt

In **Create**, enter kind `misc_feature`, positions `25` to `35`, strand
**Forward**. Use **+** to add these ordered qualifiers, checking **Value** for
both:

| Key | Value |
| --- | --- |
| `label` | `review_patch` |
| `gene` | `DEMO` |

Click **Preview**. Under **Annotations to review**, the gene and `segment_A`
should be listed for overlap and the shared `gene=DEMO` identifier. This is
informational: neither annotation will be modified, and a shared name is not
proof of a dependency. Preview does not reserve a new feature index.

Click **Create**. In this untouched starter sequence the new record is index
3; verify its label rather than assuming that index in another project. There
are now four annotations but still 120 bp of DNA.

## 4. Split and Merge: Records, Not Molecules

In **Split**, select `segment_A` and set **Split before base (1-based)** to
`31`. Preview must show:

- genomic-left record: 21..30, 10 bp;
- genomic-right record: 31..45, 15 bp;
- both records retain the original kind, strand and ordered qualifiers.

Click **Split**. There are five records. The two pieces occupy indices 1 and
2; the reverse control and `review_patch` have shifted to 3 and 4. This did
not cut or splice the DNA, and it did not establish two biological exons.

In **Merge**, choose those two `segment_A` pieces. Preview, then **Merge**.
They return to one 21..45 record and the total returns to four. GENtle only
merges exactly touching simple records with matching kind, strand and ordered
qualifiers. It refuses gaps, overlaps and metadata conflicts rather than
inventing a reconciliation. On reverse features, genomic-left/right is not
the biological 5-prime/3-prime order.

## 5. Delete, Save and Reopen

In **Delete**, select **`review_patch`**, now index 3. Preview must identify
that complete record at 25..35. Click **Delete**. Three records remain:
`DEMO`, the extended `segment_A`, and the unchanged reverse control.
Deletion removes an annotation, not sequence bases. Later feature indices can
change after a deletion, so reselect records by their current details.

Save, close, reopen, and check the same final three records and the 120 bp DNA.
Keep the starter file unchanged. Save the edited project under a different
path, and record that path in your tutorial notes.

## Command-Line Counterpart

From the repository root, use a separate new state file. The GUI and terminal
processes do not share unsaved memory. Do not point both at an actively edited
project file.

```sh
CLI="$PWD/target/debug/gentle_cli"
RUN=$(mktemp -d)
"$CLI" --state "$RUN/editor.gentle.json" op \
  '{"LoadFile":{"path":"docs/tutorial/inputs/feature_editor_demo.gb","as_id":"editor_demo"}}'
"$CLI" --state "$RUN/editor.gentle.json" shell \
  'features query editor_demo --sort feature_id --include-qualifiers'
"$CLI" --state "$RUN/editor.gentle.json" shell \
  'features edit-location editor_demo 1 --start-1based 21 --end-1based-inclusive 45 --dry-run' \
  > "$RUN/location-preview.json"
```

Inspect `.report.before`, `.report.after`, and
`.report.before_feature_fingerprint_sha256`. With `jq` available, apply the
**same** requested interval using the returned lock:

```sh
LOCK=$(jq -r '.report.before_feature_fingerprint_sha256' "$RUN/location-preview.json")
"$CLI" --state "$RUN/editor.gentle.json" shell \
  "features edit-location editor_demo 1 --start-1based 21 --end-1based-inclusive 45 --expected-feature-fingerprint-sha256 $LOCK"
```

For the other tabs, these are the matching **preview** commands, in the same
order as the GUI exercise. Split/Merge indices below assume each preceding
apply has actually completed; preview alone does not advance the exercise.

```text
features create editor_demo --kind misc_feature --start-1based 25 --end-1based-inclusive 35 --strand forward --qualifier label=review_patch --qualifier gene=DEMO --dry-run
features split editor_demo 1 --split-before-1based 31 --dry-run
features merge editor_demo 1 2 --dry-run
features delete editor_demo 3 --dry-run
```

Every apply removes `--dry-run` and supplies the exact fingerprint(s) from its
new preview. Create locks the annotation state; Delete/Split also lock the
selected feature; Merge locks both features. See the
[curation contract](../protocol.md)
and `help features create|delete|split|merge` (choose one subcommand).
Do not invent a fingerprint or reuse one after annotations have changed.

The executable regression below replays all five operations through the real
CLI, checks preview non-mutation and stale-lock refusal, compares coordinates
and qualifiers, and reopens the saved result. It needs no `jq`:

```sh
GENTLE_TUTORIAL_BIN_DIR="$PWD/target/debug" \
  python3 -m unittest scripts.test_tutorial_walkthroughs.FeatureEditorWalkthroughTests
```

## What Counts as Success?

Record each step as pass, fail or not run. The GUI checkpoint is the visible
before/after change; the scientific oracle is the saved annotation record and
unchanged DNA. Retain the exact binary revision, starter hash, final project,
preview/apply reports and any explicitly approved screenshots.

The CLI test is not a live GUI pass. There is no typed Xvfb acceptance contract
for this chapter yet: Glen must report missing semantic controls as
`harness_gap`, not infer success from screenshots or guess click positions.
Complex/nested/fuzzy boundary edits and automatic transcript repair remain
outside this exercise; unsupported edits should give a reason, not a guess.

For feedback, include tutorial ID `feature_editor_gui_cli`, step, selected
feature label and index, expected/observed coordinates, revision and retained
reports. Never substitute a real private sequence for this fixture in a public
issue. Continue with [primer-pair determination](./generated/chapters/04-02_pcr_selection_batch_primer_pairs_offline.md)
or [portable genomic regions](./generated/chapters/08-10_portable_genomic_regions_offline.md)
when the aim is a new sequence product or a saved region rather than annotation curation.
