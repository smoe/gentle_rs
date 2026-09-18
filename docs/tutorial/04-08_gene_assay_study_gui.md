# Review a Primer-Pair Study Across Alternative Transcripts

**Question:** Which primer pairs distinguish transcript structures, what
evidence informed that choice, and what must be checked before using them?

**You will learn:** open a study from Splicing Expert, inspect both primers and
their predicted products, review one changed design requirement, and export a
canonical dossier without mistaking a planned experiment for a validated assay.

**Status:** manual/hybrid; live GUI screenshots and human review pending.
This is primer-pair study design, not the separate single-primer Nanopore
capture-pool programme. Everything below is public synthetic teaching data.
These are not orderable human PATZ1 primers.

## Prepare Once

Use GUI and CLI binaries from the **same source revision**. Run the following
from the repository root, substituting the actual binary and a new output
directory whose parent exists. The helper does not build GENtle or use a model,
network service, Primer3, BLAST or vendor.

```bash
python3 scripts/prepare_gene_assay_study_tutorial.py \
  --gentle-cli /absolute/path/to/gentle_cli \
  --output-dir /absolute/path/to/new-patz1-study
```

The helper calls GENtle, rather than calculating biology in Python. It replays
the three existing [04.06 panel operations](generated/chapters/04-06_patz1_transcript_assay_panels_cli.md),
changing only their output paths, then plans the separate
[04.07 study request](inputs/transcript_assay_followup_study.json). Its receipt
records the CLI binary hash, input hashes, exact calls and artifact hashes.
It does **not** execute that study's emitted workflow.

The [fixture provenance](../../test_files/fixtures/transcript_assay_panel/patz1/README.md)
describes an invented 240 bp minus-strand locus with three annotated mature
cDNAs of 120, 80 and 100 nt. These lengths are not expression measurements.
The GRCh38-like coordinate anchor is synthetic, not a verified human reference.

## 1. Start at the Gene, Not an Oligo List

Open `patz1-study.gentle.json` with **File > Open Project**. Open its PATZ1
sequence, select the PATZ1 gene feature (`n-2`, engine feature `1`), and use
**Open Splicing Window**. In **Structure**, choose **Gene assay study...**.
PCR Designer should open in **Gene assay study** mode. Its other modes remain
available; do not choose **RT primer pool** for this exercise.

Choose **Open study plan...** and select `study.plan.json` in your output
directory. Read the question, three transcripts/three exact-cDNA groups,
annotation release, automatic recommendation and explicit profile override.
Open **Declared evidence inventory and input digests**. In particular, find
the missing differential-threshold provenance. A geometric junction constraint
and a measured effect do not by themselves establish significance.

**Checkpoint G1:** the study is inspectable but has not authorized execution.
A saved plan is a historical record, not evidence that the currently loaded
DNA or all external inputs have been revalidated. Its source feature is `n-3`
here because the supplied ledger names that transcript-linked group; the
separate 04.06 comparison panels use `n-2`. GENtle must not silently equate
their provenance. **Inspect transcript architecture in Splicing Expert** opens
the plan's declared source feature.

## 2. Inspect Both Primers and Every Transcript Outcome

Under **Persisted panels on this source sequence**, inspect
`patz1_sybr_juc_panel`. This is explicitly comparison context, not a claim that
the panel was produced by the displayed study plan. Select different entries
in **Selected primer pair**.

Read the forward and reverse sequences, each written 5' to 3', and their
binding footprints. The small blue/orange schematic uses a mature-cDNA axis:
the minus-strand genomic locus does not reverse it again. Coordinates are
zero-based half-open intervals. Tm is melting temperature, not a recommended
PCR annealing temperature.

Below the selected pair, read all three transcript outcomes. A predicted
single product, multiple products, explicit no-product result and a missing
matrix cell have different meanings. Inspect the stored selection reasons
and JUC evidence rather than assuming that the highest score explains the
selection. Shared products do not establish transcript-specific abundance.

**Checkpoint G2:** both primer footprints and the per-transcript products are
visible for one selected pair. The JUC example names PATZ1-202's exon 1-2
junction at transcript position 40. The synthetic effect is -1.2, but the
missing threshold remains visible and is not promoted to qualified evidence.

## 3. Review the Checks, Not Just the Sequences

Choose **Open panel design, matrix and readiness checks**. Review the
**Transcript x assay product matrix** and **Requested junctions**. The old
sequence table is now labelled **Candidate primer sequences (not order approval)**.
Choose **Build experimental handoff** to run the shared readiness operation.

**Checkpoint G3:** missing genomic/transcriptome specificity is a blocker, not
a green badge. The panel's design-only result has not launched BLAST or
submitted an order. A complete transcript matrix is not whole-reference
specificity. See [04.07's specificity continuation](04-07_transcript_assay_followup_gui_cli.md#5-optional-external-specificity-plan-is-not-pass)
for real data with prepared, authentic reference resources.

## 4. Change One Requirement and Review Again

Return to **Gene assay study**, then **Open planning request...** and choose
`study.request.json`. Expand **Design and review a new iteration**.
Change `plan_id` in the JSON to a new, explicit identifier such as
`patz1_gui_review_2`. Change **Short-product maximum** from 120 to 110 bp.
This is a hard limit, not a ranking preference; the tiny fixture's limits are
not recommendations for a laboratory assay.

Choose **Normalize request** and expand **Effective normalized request**.
Inspect all defaults, input digests, coverage policy, override and missing
evidence. Only after reviewing it, tick the planning-review checkbox, supply
a **new absolute output directory**, and choose **Plan reviewed study**.

**Checkpoint G4:** GENtle writes `request.json`, `plan.json` and the exact
`workflow.json` to the new directory. The original study is untouched. Changing
the request again or changing project structure invalidates review. Automatic
recommendation and user override remain separate.

Inspect **Exact ordered operations** and both operation/workflow digests.
Planning success is not primer-design success. The planner uses stricter
default primer constraints than the permissive 04.06 comparison recipe;
neither that recipe's successful design nor its reports prove this plan feasible.

If you deliberately want to test execution, separately tick the execution-review
checkbox and choose **Execute reviewed workflow**. Observe progress and use
**Cancel study task** when needed. Record a typed design failure as a failure;
do not loosen constraints, patch the workflow or manufacture a pass for this
teaching exercise. Existing panel identities are refused rather than overwritten.
Cancellation/staleness discards engine changes, but does not erase files already
written. New approval is required after a changed request or failed iteration.

**Checkpoint G5:** no second approval means no execution. Editing workflow
bytes must fail digest verification before design. A cancelled or stale result
must not become the current project result. This checkpoint is about safe
execution, not a promised successful study with arbitrary constraints.

## 5. Export an Honest Dossier

Expand **Canonical dossier export**. Enter the absolute path to the helper's
`publication.request.json` and a new absolute output directory, for example
`.../gui-dossier`. Leave optional PDF off for the bounded offline exercise,
then choose **Export declared dossier**.

Open its `index.html` and gene page. Compare `canonical-report.json` with
the helper's `cli-dossier/canonical-report.json`; the scientific records should
agree. The request binds the **original** unexecuted study, not your new plan.
Changing which study is published requires an explicitly updated, correctly
hash-bound publication request; the GUI must not silently replace it.

**Checkpoint G6:** the page says **Pending** and explains that the study's
assays have not been executed. The separate comparison panels are not inserted
as falsely plan-bound handoffs. Real completed handoffs and reviewed order forms
can be supplied using the existing canonical contract, but this tutorial does
not fabricate them. There is no order-submission action.

## Shell and Agent Parity

The GUI Shell can open the existing designer with
`ui open pcr-design patz1_transcript_assay_demo`; choose **Gene assay study**.
The shared routes remain `primers plan-gene-isoform-study`,
`primers execute-gene-isoform-study-workflow`, and
`primers publish-gene-isoform-study`. Use the exact saved request/workflow files;
GUI Shell is not Bash. CLI processes do not inherit unsaved GUI state.

An optional inner-agent prompt is: "Help me review the declared evidence and
primer pairs in the synthetic PATZ1 gene assay study. Open PCR Designer, do not
execute design or order anything until I separately review the exact inputs."
This is a static example, not a tested provider conversation or delegated
scientific approval. Help can discover this guide with
`ui open tutorial-guide gene_assay_study_gui`.

## Troubleshooting and Review

No panels: make sure you opened the helper's project, not just its sequence or
plan JSON. Wrong-locus error: open the plan's exact source sequence rather than
renaming another sequence to match. Directory exists: choose a new run location.
Input/hash mismatch: retain the failure and review original inputs, never edit
hashes merely to make execution proceed. Missing PDF: optional printing needs a
browser; HTML export is a separate, inspectable artifact.

Glen's [screenshot and real-data acceptance request](../gene_assay_study_glen_acceptance.md)
specifies raw captures, annotations and exact-revision receipts. Automated tests
do not replace his live GUI review or the user's biological sign-off.
