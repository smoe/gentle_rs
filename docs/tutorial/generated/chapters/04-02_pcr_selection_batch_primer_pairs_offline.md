---
chapter_id: "pcr_selection_batch_primer_pairs_offline"
title: "Determine and review PCR primer pairs (offline)"
tier: "core"
example_id: "pcr_selection_batch_primer_pairs_offline"
source_example: "docs/examples/workflows/pcr_selection_batch_primer_pairs_offline.json"
example_test_mode: "always"
executed_during_generation: true
automated_status: "passing"
review_status: "codex_reviewed"
review_stale: false
codex_reviewed_at: "2026-07-28"
human_reviewed_at: null
human_reviewer: null
review_stale_reason: null
review_issue_template: null
review_issue_template_path: null
generated_artifact_dir: "docs/tutorial/generated/artifacts/pcr_selection_batch_primer_pairs_offline"
---

# Determine and review PCR primer pairs (offline)

Separate the indispensable biological ROI from primer-search flanks, generate real ranked candidates on local TP73 context, and decide what the report does and does not justify before ordering primers.

Primer-pair determination begins with geometry, not with a request for two oligos. The green ROI is the smallest biological interval that every acceptable amplicon must contain. The upstream and downstream search windows are separate intervals in which GENtle may place the forward and reverse primers. Treating a broad context interval as the ROI while also placing primers inside it can make the request impossible; this chapter deliberately avoids that common mistake.

The canonical offline workflow uses two GC-rich TP73-AS2 promoter-proximal examples and the deterministic internal backend. Each run returns five ranked pairs and a complete `gentle.primer_design_report.v1` JSON artifact. The first rank-1 pair satisfies every reported heuristic flag. The second rank-1 pair carries a homopolymer/secondary-structure advisory, which demonstrates why rank 1 means only 'best under these stated constraints and scoring rules.' It is not an order recommendation, proof of amplification, or genome-wide specificity result. Review the sequences, coordinates, Tm/GC balance, local annealing-hit counts, 3-prime clamps, self- and pair-complementarity, amplicon geometry, backend provenance, and rejection summary together. External genomic specificity and wet-lab validation remain separate follow-up steps.

**Prerequisites:** Read [Chapter 1: Load FASTA, branch, and reverse-complement](./02-01_load_branch_reverse_complement_pgex_fasta.md) first.

## Parameters That Matter

- `Core ROI: 61820..61920 or 62040..62140` (where used: green ROI and `DesignPrimerPairs.roi_start_0based` / `roi_end_0based`)
  - Why it matters: The ROI is the indispensable biology every product must include. Making the whole context the ROI can make flanking-primer geometry impossible.
  - How to derive it: Choose the smallest interval needed to answer the experiment, then verify that an acceptable product can place one primer on each side.
- `Primer-search windows` (where used: red/blue paint and forward/reverse `start_0based` / `end_0based` constraints)
  - Why it matters: They delimit candidate binding sites independently of the core ROI and keep the search small enough to evaluate exhaustively.
  - How to derive it: Place non-overlapping windows outside the ROI, with enough room for 18..24 nt primers and an allowed product no longer than 350 bp.
- `Primer limits: 18..24 nt, 55..78 C, GC 0.30..0.80, local annealing hits <= 1` (where used: forward and reverse side constraints)
  - Why it matters: They define which oligos enter pair scoring. The deliberately broad upper Tm accommodates this GC-rich teaching locus; it is not a universal laboratory default.
  - How to derive it: Choose limits from template composition, chemistry, and intended cycling conditions, then record the actual values used rather than relying on hidden defaults.
- `Pair limits: 90..350 bp product, <= 5 C Tm difference, maximum 5 returned pairs` (where used: pair constraints and candidate ranking)
  - Why it matters: The product limit keeps the assay practical; the Tm-difference bound checks pair compatibility; retaining alternatives supports review rather than blind rank-1 acceptance.
  - How to derive it: Set product length from the experimental workflow and keep enough alternatives to inspect meaningful trade-offs.
- `Heuristic rule flags and rejection summary` (where used: ranked pair rows and report-level diagnostics)
  - Why it matters: A high score can coexist with a secondary-structure or homopolymer advisory. Rejection counts explain the search but do not establish whole-genome specificity.
  - How to derive it: Review every flag and both primers before selection; preserve a concise selection/rejection rationale with the exported report.

## When This Routine Is Useful

- You want to understand why the biological ROI and the two primer-search windows must be distinct.
- You want to determine and compare several candidate pairs rather than accepting the first sequences returned.
- You want to inspect the evidence behind pair ranking: product geometry, Tm/GC balance, local uniqueness, clamps, complementary runs, and rule flags.
- You want to learn from both a clean top-ranked candidate and a top-ranked candidate that still has a heuristic advisory.
- You want deterministic report IDs and complete JSON artifacts for review or downstream handoff.
- You want to distinguish local candidate generation from later whole-genome specificity and wet-lab validation.

## What You Learn

- Separate a biologically indispensable core ROI from the two allowed primer-placement windows.
- Read a complete primer-design report and explain why a pair received its rank.
- Compare candidate pairs using amplicon geometry, Tm/GC, local annealing hits, clamps, secondary-structure proxies, dimer proxies, and rule flags.
- Recognize that an accepted or top-ranked in-silico pair can still carry advisories and still requires external specificity and experimental validation.
- Reproduce the same determination through GUI, CLI, workflow, and exported JSON.

## Applied Concepts

- **Shared Engine Contract** (`shared_engine_contract`): GUI, CLI, shell, and scripting interfaces execute the same operation semantics.
- **Deterministic Workflows** (`deterministic_workflows`): Operation chains should produce stable IDs and comparable outputs across repeated runs.
- **Artifact Exports** (`artifact_exports`): Representative outputs (CSV/protocol/SVG/text) are retained for auditability and sharing.
- **Tutorial Drift Checks** (`tutorial_drift_checks`): Tutorial content is generated from executable examples and verified in automated checks.

## At a Glance

1. Load test_files/tp73.ncbi.gb, then open Patterns -> PCR Designer... (or comma...
2. For the first example, paint only the indispensable core 61820..61920 as the ...
3. Paint the upstream search window 61720..61745 red and the downstream search w...
4. Set primer length 18..24 nt, Tm 55..78 C, GC fraction 0.30..0.80, maximum loc...
5. Run the first design and inspect report tp73_as2_promoter_batch_r01. Confirm ...
6. Open rank 1 and review both oligos plus the pair: sequences and binding coord...
7. Repeat with ROI 62040..62140, upstream 61940..61965, and downstream 62158..62...
8. Use Export to retain each full JSON report. Record why you selected or reject...

## GUI First

CLI snippets use GENtle's default `.gentle_state.json` state unless they say otherwise. Add `--state PATH` or `--project PATH` when you want an explicit sandboxed state file for copied commands.

### Step 1: Load test_files/tp73.ncbi.gb, then open Patterns -> PCR Designer... (or comma...

GUI: Load `test_files/tp73.ncbi.gb`, then open `Patterns -> PCR Designer...` (or command palette `PCR Designer`) and keep the map in linear mode.

CLI:

```bash
cargo run --bin gentle_cli -- --state /tmp/gentle-primer-pair-tutorial.json workflow @docs/examples/workflows/pcr_selection_batch_primer_pairs_offline.json
```

> Expected: The workflow creates two complete reports and exports both as retained tutorial artifacts.

### Step 2: For the first example, paint only the indispensable core 61820..61920 as the ...

GUI: For the first example, paint only the indispensable core `61820..61920` as the green ROI. The broader copied context `61720..62120` is not the ROI.

CLI:

```bash
jq '{report_id,template,required_roi:[.roi_start_0based,.roi_end_0based]}' docs/tutorial/generated/artifacts/pcr_selection_batch_primer_pairs_offline/artifacts/tp73_as2_promoter_batch_r01.report.json
```

> Expected: The first report records the indispensable core as `61820..61920`; the broader extracted context is not silently substituted as the ROI.

### Step 3: Paint the upstream search window 61720..61745 red and the downstream search w...

GUI: Paint the upstream search window `61720..61745` red and the downstream search window `62000..62040` blue. Confirm in the live cartoon that both primer windows lie outside and flank the ROI.

CLI:

```bash
jq '{required_roi:[.roi_start_0based,.roi_end_0based],forward_window:[.forward.start_0based,.forward.end_0based],reverse_window:[.reverse.start_0based,.reverse.end_0based],require_roi_flanking:.pair_constraints.require_roi_flanking}' docs/tutorial/generated/artifacts/pcr_selection_batch_primer_pairs_offline/artifacts/tp73_as2_promoter_batch_r01.report.json
```

> Expected: The upstream and downstream windows lie outside the required ROI, and `require_roi_flanking` is true.

### Step 4: Set primer length 18..24 nt, Tm 55..78 C, GC fraction 0.30..0.80, maximum loc...

GUI: Set primer length `18..24 nt`, Tm `55..78 C`, GC fraction `0.30..0.80`, maximum local annealing hits `1`, amplicon length `90..350 bp`, maximum pair Tm difference `5 C`, backend `internal`, and maximum returned pairs `5`.

CLI:

```bash
jq '{primer_length:[.forward.min_length,.forward.max_length],primer_tm_c:[.forward.min_tm_c,.forward.max_tm_c],primer_gc_fraction:[.forward.min_gc_fraction,.forward.max_gc_fraction],max_local_hits:.forward.max_anneal_hits,amplicon_bp:[.min_amplicon_bp,.max_amplicon_bp],max_tm_delta_c,max_pairs}' docs/tutorial/generated/artifacts/pcr_selection_batch_primer_pairs_offline/artifacts/tp73_as2_promoter_batch_r01.report.json
```

> Expected: The report preserves the exact primer, amplicon, Tm-difference, local-hit, and result-count limits used for determination.

### Step 5: Run the first design and inspect report tp73_as2_promoter_batch_r01. Confirm ...

GUI: Run the first design and inspect report `tp73_as2_promoter_batch_r01`. Confirm that it contains five pairs, records the internal backend, and did not skip candidate combinations because of the evaluation limit.

CLI:

```bash
jq '{pair_count,backend,rejection_summary}' docs/tutorial/generated/artifacts/pcr_selection_batch_primer_pairs_offline/artifacts/tp73_as2_promoter_batch_r01.report.json
```

> Expected: The first report has five candidates, records the internal backend, and has `pair_evaluation_limit_skipped: 0`.

### Step 6: Open rank 1 and review both oligos plus the pair: sequences and binding coord...

GUI: Open rank 1 and review both oligos plus the pair: sequences and binding coordinates, 310 bp amplicon, Tm/GC values, one local annealing hit per primer, 3-prime clamps, homopolymer/self-complementary runs, pair-complementary runs, score, and every rule flag. Compare ranks 1-3 instead of treating rank 1 as automatically order-ready.

CLI:

```bash
jq '{top_pair:(.pairs[0]|{rank,score,amplicon:[.amplicon_start_0based,.amplicon_end_0based_exclusive,.amplicon_length_bp],tm_delta_c,forward,reverse,primer_pair_complementary_run_bp,primer_pair_3prime_complementary_run_bp,rule_flags}),alternatives:[.pairs[0:3][]|{rank,score,amplicon_length_bp,tm_delta_c}]}' docs/tutorial/generated/artifacts/pcr_selection_batch_primer_pairs_offline/artifacts/tp73_as2_promoter_batch_r01.report.json
```

> Expected: Rank 1 exposes both full oligo diagnostics and pair diagnostics, while ranks 1-3 remain directly comparable.

### Step 7: Repeat with ROI 62040..62140, upstream 61940..61965, and downstream 62158..62...

GUI: Repeat with ROI `62040..62140`, upstream `61940..61965`, and downstream `62158..62188`. In `tp73_as2_promoter_batch_r02`, notice that the top forward primer has a five-base homopolymer and `forward_secondary_structure_ok: false`; the report preserves this advisory even though the pair ranks first.

CLI:

```bash
jq '{pair_count,rejection_summary,top_forward:(.pairs[0].forward|{sequence,tm_c,gc_fraction,anneal_hits,longest_homopolymer_run_bp,self_complementary_run_bp}),top_flags:.pairs[0].rule_flags}' docs/tutorial/generated/artifacts/pcr_selection_batch_primer_pairs_offline/artifacts/tp73_as2_promoter_batch_r02.report.json
```

> Expected: The second report exposes a real advisory (`forward_secondary_structure_ok: false`) rather than hiding it behind the rank.

### Step 8: Use Export to retain each full JSON report. Record why you selected or reject...

GUI: Use `Export` to retain each full JSON report. Record why you selected or rejected a candidate, then perform the separate genomic-specificity handoff and wet-lab validation before procurement.

CLI:

```bash
cargo run --bin gentle_cli -- --state /tmp/gentle-primer-pair-tutorial.json shell 'primers export-report tp73_as2_promoter_batch_r01 /tmp/tp73_as2_promoter_batch_r01.json'
```

> Expected: The shared export route writes the same engine-owned report shown by the GUI and retained by the tutorial generator.


## Follow-up Commands

```bash
cargo run --bin gentle_cli -- shell 'ui open pcr-design'
cargo run --bin gentle_cli -- shell 'ui focus pcr-design'
cargo run --bin gentle_cli -- --state /tmp/gentle-primer-pair-tutorial.json workflow @docs/examples/workflows/pcr_selection_batch_primer_pairs_offline.json
cargo run --bin gentle_cli -- --state /tmp/gentle-primer-pair-tutorial.json shell 'primers list-reports'
cargo run --bin gentle_cli -- --state /tmp/gentle-primer-pair-tutorial.json shell 'primers show-report tp73_as2_promoter_batch_r01'
cargo run --bin gentle_cli -- --state /tmp/gentle-primer-pair-tutorial.json shell 'primers export-report tp73_as2_promoter_batch_r01 /tmp/tp73_as2_promoter_batch_r01.json'
cargo run --bin gentle_cli -- --state /tmp/gentle-primer-pair-tutorial.json shell 'primers specificity-plan tp73_as2_promoter_batch_r01 --pair-rank 1 --target-genome GRCh38.p14 --output-dir /tmp/tp73_primer_specificity'
```

## Checkpoints

- The first report deterministically contains five pairs, uses the internal backend, and has `pair_evaluation_limit_skipped: 0`.
- The first rank-1 pair spans 310 bp, has one local annealing hit per primer, and has every reported rule flag set to true.
- The second report also contains five pairs without evaluation-limit skipping.
- The second rank-1 forward primer has a five-base homopolymer and `forward_secondary_structure_ok: false`; the tutorial therefore does not call it order-ready.
- Both JSON exports preserve the ROI, primer-window constraints, backend provenance, ranked alternatives, diagnostics, and rejection summary.
- The chapter explicitly distinguishes deterministic local ranking from whole-genome specificity, successful amplification, and experimental validation.

## What This Chapter Produces

- [`artifacts/pcr_selection_batch_primer_pairs_offline/artifacts/tp73_as2_promoter_batch_r01.report.json`](../artifacts/pcr_selection_batch_primer_pairs_offline/artifacts/tp73_as2_promoter_batch_r01.report.json) - schema: `gentle.primer_design_report.v1`
- [`artifacts/pcr_selection_batch_primer_pairs_offline/artifacts/tp73_as2_promoter_batch_r02.report.json`](../artifacts/pcr_selection_batch_primer_pairs_offline/artifacts/tp73_as2_promoter_batch_r02.report.json) - schema: `gentle.primer_design_report.v1`

## Tutorial Provenance

- Chapter id: `pcr_selection_batch_primer_pairs_offline`
- Tier: `core`
- Example id: `pcr_selection_batch_primer_pairs_offline`
- Tutorial source JSON: `docs/tutorial/sources/04-02_pcr_selection_batch_primer_pairs_offline.json`
- Workflow file: `docs/examples/workflows/pcr_selection_batch_primer_pairs_offline.json`
- Generated artifact dir: `docs/tutorial/generated/artifacts/pcr_selection_batch_primer_pairs_offline`
- Example test_mode: `always`
- Executed during generation: `yes`
- Automated status: `passing`
- Review status: `codex_reviewed`
- Codex reviewed at: `2026-07-28`
- Human reviewed at: `not recorded`
- Inspect the source JSON when you need full option-level detail.

## Feedback

If this tutorial is confusing, execution-stale, biologically suspect, or missing a useful figure, please open the matching tutorial issue template and include the context below.

- Tutorial title: `Determine and review PCR primer pairs (offline)`
- Tutorial/chapter id: `pcr_selection_batch_primer_pairs_offline`
- Step reached:
- Expected vs. actual:
- Interface used: GUI / CLI / Agent Assistant / ClawBio

Paste the Tutorial feedback context here:

```text

```
