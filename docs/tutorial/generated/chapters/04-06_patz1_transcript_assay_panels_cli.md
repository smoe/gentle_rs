---
chapter_id: "patz1_transcript_assay_panels_cli"
title: "Design Practical PATZ1 Endpoint and SYBR Transcript Panels from the CLI"
tier: "core"
example_id: "patz1_endpoint_sybr_transcript_assay_panel_offline"
source_example: "docs/examples/workflows/patz1_endpoint_sybr_transcript_assay_panel_offline.json"
example_test_mode: "always"
executed_during_generation: true
automated_status: "passing"
review_status: "codex_reviewed"
review_stale: false
codex_reviewed_at: "2026-07-21"
human_reviewed_at: null
human_reviewer: null
review_stale_reason: null
review_issue_template: null
review_issue_template_path: null
generated_artifact_dir: "docs/tutorial/generated/artifacts/patz1_transcript_assay_panels_cli"
---

# Design Practical PATZ1 Endpoint and SYBR Transcript Panels from the CLI

Use one synthetic minus-strand PATZ1 locus to generate an endpoint first-end x terminal-end band matrix, a Clariom-JUC-directed primer-only SYBR panel, and an annotation-confirmed routine common-region screen, then compose the persisted reports with a digest-bound isoform-evidence ledger entirely from the command line.

A multi-transcript PCR experiment has three distinct jobs. Long endpoint RT-PCR can expose transcript families as differently sized first-end x terminal-end products on a gel, although band intensity is only rough or semi-quantitative. Short SYBR assays can test selected exon-exon junctions without requiring an internal hydrolysis probe. A routine common-region assay can provide a practical screen across the intended transcript set. GENtle keeps these jobs in one `gentle.transcript_assay_panel.v2` contract while recording their experimental tier independently from the selection objective.

This chapter uses a deliberately sequence-diverse, synthetic PATZ1-like locus on the minus strand. It demonstrates exact mature-cDNA equivalence classes, a 10 kb allowed endpoint ceiling, caller-supplied preferred product ranges, oligo-dT reverse-transcription warnings, and distinct synthetic Clariom JUC and PSR evidence. Transcript annotation alone establishes whether a selected amplicon is common; PSR signal remains independent support and cannot manufacture that claim. The fixture does not produce orderable primers for the human PATZ1 locus and does not claim that a missing long product proves isoform absence. The canonical workflow is the deterministic reference; the direct CLI steps extract the same three operation objects and submit them unchanged through `primers design-transcript-assay-panel`.

The final step composes, rather than redesigns, the three persisted panels with a checked `gentle.gene_isoform_evidence.v2` ledger. The resulting `gentle.gene_transcript_assay_routine.v1` report retains panel and evidence digests, assay roles, uncovered transcript classes, candidate order rows, and a conservative experimental sequence. It cannot turn missing specificity into a pass: the synthetic panels remain visibly unaccepted until the declared external checks are completed.

Each selected assay also records a `primers specificity-plan` follow-up template. Planning creates structured BLAST jobs but does not execute them. For a real complete panel, prepare or import catalogued genomic-DNA and whole-cDNA BLAST resources, then use `primers transcript-assay-specificity-plan` and always call `primers transcript-assay-specificity-finalize` with the scheduler's execution manifest. GENtle distinguishes a completed biological `specificity_fail`, a completed but geometrically `not_assessed` result, and incomplete or failed execution evidence; only `pass` is accepted. Until a handoff is actually executed and finalized, this offline tutorial honestly retains `genomic_confirmation_status: not_run`.

**Prerequisites:** Read [Chapter 13: Determine and review PCR primer pairs (offline)](./04-02_pcr_selection_batch_primer_pairs_offline.md) first.

## Parameters That Matter

- `assay_kind=endpoint_rt_pcr / objective=isoform_end_matrix` (where used: The first `DesignTranscriptAssayPanel` operation)
  - Why it matters: This combination constructs only annotated first-end x terminal-end reactions and reports their product sizes across all transcript classes.
  - How to derive it: Use endpoint mode when differently sized gel bands are the intended readout; use a short qPCR mode for junction validation.
- `cdna_synthesis=oligo_dt / max_amplicon_bp=10000` (where used: Endpoint panel design and interpretation warnings)
  - Why it matters: PCR capacity can reach 10 kb, but oligo-dT reverse transcription may still underrepresent distant 5-prime sequence.
  - How to derive it: Set the PCR ceiling from the available long-range polymerase, then retain short end-specific assays or 5-prime RACE as confirmation for vulnerable long products.
- `coverage_policy=require_all` (where used: Both tutorial operations)
  - Why it matters: Strict coverage refuses incomplete panels instead of silently dropping a mature-cDNA class or requested junction.
  - How to derive it: Keep `require_all` for a claimed complete panel; choose `best_effort` only as an explicit exploratory decision and inspect every uncovered reason.
- `assay_tier / preferred and allowed amplicon ranges` (where used: All three transcript-panel operations)
  - Why it matters: The experimental purpose is independent of the selection objective, and product length is a structured preference after required biology rather than one hidden cutoff.
  - How to derive it: Set the allowed range from the assay chemistry and polymerase. Supply a narrower preferred range only when routine cycling and gel handling matter; inspect explicit long-range fallbacks instead of silently discarding them.
- `junction_evidence_priority=required` (where used: The primer-only SYBR operation)
  - Why it matters: Every supplied Clariom-style JUC target must be evaluated and reported rather than competing with only the automatically selected anchors.
  - How to derive it: Use `required` for junctions selected for laboratory validation and `preferred` for evidence that may guide but not block panel completion.
- `--backend internal versus --backend primer3` (where used: Direct CLI execution)
  - Why it matters: The backend override changes candidate generation while leaving the operation and report contract unchanged.
  - How to derive it: Use the internal backend for this deterministic tutorial; use `primers preflight --backend primer3` before selecting Primer3 for real designs.
- `primers specificity-plan / specificity-import` (where used: Per-assay genomic-confirmation follow-up)
  - Why it matters: The handoff lets an outer workflow own long-running BLAST processes and use their exit status as the completion signal, while GENtle remains responsible for deterministic queries, policy, provenance, and interpretation.
  - How to derive it: Generate the handoff against a prepared genome, run every returned `program` with its exact `args[]`, require exit code 0, and only then import the handoff. Do not infer completion from a non-empty output file.
- `transcript-assay-specificity-plan / transcript-assay-specificity-finalize` (where used: Complete-panel genomic or whole-cDNA confirmation after this offline tutorial)
  - Why it matters: Aggregate finalization verifies every declared command and output identity atomically, inherits the panel's allowed product ceiling for readiness, and keeps biological failure, unavailable target geometry, and execution failure separate.
  - How to derive it: Prepare a catalogued `genomic_dna` or `transcriptome_cdna` BLAST resource, generate one panel handoff, execute every declared command, and return the completed execution manifest even when a process fails. Accept only a final `pass`.

## When This Routine Is Useful

- You want a command-line-only route for endpoint RT-PCR, SYBR qPCR, or TaqMan-compatible transcript-panel design.
- You want a first-end x terminal-end reaction matrix with predicted band sizes across every annotated mature transcript.
- You want short primer-only assays directed at required or preferred exon-exon junction evidence.
- You want a routine common-region screen whose structural claim comes from transcript annotation while PSR support remains separate.
- You need to preserve every shared operation field instead of relying only on convenience flags.
- You want one digest-bound review manifest joining existing common-control, junction-validation, and endpoint reports without rerunning them.
- You want an external scheduler to own BLAST process completion while GENtle retains deterministic specificity inputs and interpretation.
- You need complete-panel acceptance to distinguish a biological specificity failure from unavailable target geometry or failed external execution.
- You want persisted reports that can be listed, shown, exported, or consumed through MCP, JavaScript, Lua, and workflows.

## What You Learn

- Distinguish endpoint isoform discovery from short primer-only junction validation.
- Run a complete externally tagged `DesignTranscriptAssayPanel` operation through the direct primer CLI route.
- Interpret exact mature-cDNA equivalence classes and transcript-by-assay product matrices conservatively.
- Read first-end x terminal-end reactions and predicted band sizes without treating a missing long product as proof of transcript absence.
- Require every supplied junction to be evaluated or returned with an explicit unresolved reason.
- Separate BLAST job planning and execution from GENtle's import-time specificity interpretation.
- Interpret `pass`, `specificity_fail`, `not_assessed`, and `incomplete` as distinct whole-panel outcomes.
- Separate transcript-annotation commonality, PSR support, product-length practicality, and the existing primer score in the selection rationale.
- Compose existing evidence and assay panels into one digest-bound review manifest without rerunning design or specificity.
- Distinguish one measured differential-junction observation from its transcript projections, and verify that a selected-pair audit is bound to the exact panel-design operation.
- Persist, list, inspect, and export transcript assay reports from one explicit project state.

## Applied Concepts

- **Shared Engine Contract** (`shared_engine_contract`): GUI, CLI, shell, and scripting interfaces execute the same operation semantics.
- **Deterministic Workflows** (`deterministic_workflows`): Operation chains should produce stable IDs and comparable outputs across repeated runs.
- **Artifact Exports** (`artifact_exports`): Representative outputs (CSV/protocol/SVG/text) are retained for auditability and sharing.
- **Tutorial Drift Checks** (`tutorial_drift_checks`): Tutorial content is generated from executable examples and verified in automated checks.

## At a Glance

1. No GUI is required. Optionally open the synthetic fixture in PCR Designer to ...
2. If using the GUI for comparison, select the PATZ1-like gene feature and confi...
3. Use the endpoint mode with oligo-dT cDNA, isoform_end_matrix, strict coverage...
4. Review Primer3 preflight only when using the external backend; the determinis...
5. Inspect the endpoint reaction and band-size matrices, including the reverse-t...
6. Use primer-only SYBR mode with required junction evidence to mirror the secon...
7. Use the routine common-region tier with a pan-transcript objective to mirror ...
8. Inspect the persisted panel rows and export the same report JSON used by the ...
9. GENtle does not yet provide a dedicated GUI composer for the compound routine...
10. Inspect the returned routine roles, evidence and panel digests, uncovered cla...
11. Inspect the per-assay specificity handoff templates; replace GENOME_ID and OU...

## GUI First

CLI snippets use GENtle's default `.gentle_state.json` state unless they say otherwise. Add `--state PATH` or `--project PATH` when you want an explicit sandboxed state file for copied commands.

### Step 1: No GUI is required. Optionally open the synthetic fixture in PCR Designer to ...

GUI: No GUI is required. Optionally open the synthetic fixture in PCR Designer to inspect the same transcript classes and assay modes graphically.

CLI:

```bash
cargo run --bin gentle_cli -- workflow @docs/examples/workflows/patz1_endpoint_sybr_transcript_assay_panel_offline.json
```

> Expected: The canonical offline workflow writes complete endpoint, discriminating SYBR, and routine common-region `gentle.transcript_assay_panel.v2` reports from committed synthetic inputs.

### Step 2: If using the GUI for comparison, select the PATZ1-like gene feature and confi...

GUI: If using the GUI for comparison, select the PATZ1-like gene feature and confirm that its transcript order follows mature 5-prime to 3-prime orientation despite the minus-strand genomic locus.

CLI:

```bash
cargo run --bin gentle_cli -- --state /tmp/gentle-patz1-transcript-panels.json op '{"LoadFile":{"path":"test_files/fixtures/transcript_assay_panel/patz1/patz1_assay_minus_strand.gb","as_id":"patz1_transcript_assay_demo"}}'
```

> Expected: The explicit project state contains the minus-strand PATZ1-like locus under the deterministic sequence id used by both operations.

### Step 3: Use the endpoint mode with oligo-dT cDNA, isoform_end_matrix, strict coverage...

GUI: Use the endpoint mode with oligo-dT cDNA, `isoform_end_matrix`, strict coverage, and a 10,000 bp maximum product to mirror the first CLI operation.

CLI:

```bash
jq '.workflow.ops[2] | .DesignTranscriptAssayPanel.path="/tmp/patz1_endpoint_end_matrix.report.json"' docs/examples/workflows/patz1_endpoint_sybr_transcript_assay_panel_offline.json > /tmp/patz1_endpoint.operation.json && jq '.workflow.ops[3] | .DesignTranscriptAssayPanel.path="/tmp/patz1_sybr_juc_panel.report.json"' docs/examples/workflows/patz1_endpoint_sybr_transcript_assay_panel_offline.json > /tmp/patz1_sybr.operation.json && jq '.workflow.ops[4] | .DesignTranscriptAssayPanel.path="/tmp/patz1_routine_common_region_screen.report.json"' docs/examples/workflows/patz1_endpoint_sybr_transcript_assay_panel_offline.json > /tmp/patz1_routine.operation.json
```

> Expected: The three `@FILE` payloads are the exact externally tagged operations from the canonical workflow, with only their output paths redirected to `/tmp`.

### Step 4: Review Primer3 preflight only when using the external backend; the determinis...

GUI: Review Primer3 preflight only when using the external backend; the deterministic tutorial itself uses the internal backend.

CLI:

```bash
cargo run --bin gentle_cli -- --state /tmp/gentle-patz1-transcript-panels.json primers preflight --backend primer3
```

> Expected: Primer3 preflight reports whether the external executable is reachable; its result does not alter the internal-backend tutorial run.

### Step 5: Inspect the endpoint reaction and band-size matrices, including the reverse-t...

GUI: Inspect the endpoint reaction and band-size matrices, including the reverse-transcription-completeness warning for long 5-prime reach.

CLI:

```bash
cargo run --bin gentle_cli -- --state /tmp/gentle-patz1-transcript-panels.json primers design-transcript-assay-panel @/tmp/patz1_endpoint.operation.json --backend internal
```

> Expected: The endpoint report is complete, primer-only, uses the 10,000 bp ceiling, and contains designed end reactions plus at least three distinct predicted band sizes.

### Step 6: Use primer-only SYBR mode with required junction evidence to mirror the secon...

GUI: Use primer-only SYBR mode with required junction evidence to mirror the second CLI operation; no probe should be created.

CLI:

```bash
cargo run --bin gentle_cli -- --state /tmp/gentle-patz1-transcript-panels.json primers design-transcript-assay-panel @/tmp/patz1_sybr.operation.json --backend internal
```

> Expected: The SYBR report is primer-only and includes a selected assay spanning the required synthetic Clariom JUC junction.

### Step 7: Use the routine common-region tier with a pan-transcript objective to mirror ...

GUI: Use the routine common-region tier with a pan-transcript objective to mirror the third CLI operation; compare its annotation evidence, preferred-range classification, PSR/JUC rows, and rejected alternatives.

CLI:

```bash
cargo run --bin gentle_cli -- --state /tmp/gentle-patz1-transcript-panels.json primers design-transcript-assay-panel @/tmp/patz1_routine.operation.json --backend internal
```

> Expected: The routine report uses `pan_transcript`, records annotation-confirmed common-region evidence, and treats any overlapping PSR row only as independent support.

### Step 8: Inspect the persisted panel rows and export the same report JSON used by the ...

GUI: Inspect the persisted panel rows and export the same report JSON used by the CLI path.

CLI:

```bash
cargo run --bin gentle_cli -- --state /tmp/gentle-patz1-transcript-panels.json primers list-transcript-assay-panels && cargo run --bin gentle_cli -- --state /tmp/gentle-patz1-transcript-panels.json primers show-transcript-assay-panel patz1_routine_common_region_screen && cargo run --bin gentle_cli -- --state /tmp/gentle-patz1-transcript-panels.json primers export-transcript-assay-panel patz1_sybr_juc_panel /tmp/patz1_sybr_juc_panel.export.json
```

> Expected: Both report ids are persisted in the selected state; list/show/export operate without reconstructing or hand-editing the reports.

### Step 9: GENtle does not yet provide a dedicated GUI composer for the compound routine...

GUI: GENtle does not yet provide a dedicated GUI composer for the compound routine. Use the GUI Shell with `primers compose-gene-assay-routine` and the same JSON request when an in-application handoff is useful; this reaches the shared operation rather than GUI-local logic.

CLI:

```bash
cargo run --bin gentle_cli -- --state /tmp/gentle-patz1-transcript-panels.json primers compose-gene-assay-routine '{"label":"Synthetic PATZ1 common-control, junction, and endpoint review","isoform_evidence_path":"test_files/fixtures/transcript_assay_panel/patz1/patz1_assay_isoform_evidence.json","expected_isoform_evidence_sha256":"sha256:267a8fc00da8fea653dd440259987b20fc68422554cfed8bbff61809d4eeaaf5","transcript_assay_panel_report_ids":["patz1_routine_common_region_screen","patz1_sybr_juc_panel","patz1_endpoint_end_matrix"]}' --path /tmp/patz1_gene_transcript_assay_routine.report.json
```

> Expected: The compound routine is a pure read over the exact evidence bytes and persisted panel reports; its three roles are common_control, junction_validation, and endpoint_structure.

### Step 10: Inspect the returned routine roles, evidence and panel digests, uncovered cla...

GUI: Inspect the returned routine roles, evidence and panel digests, uncovered classes, and specificity status. Composition is review-only and must not alter the selected sequence or silently rerun a panel.

CLI:

```bash
jq '{schema, isoform_evidence_sha256, roles:[.assay_panels[].role], specificity_accepted:[.assay_panels[].specificity_accepted], experimental_sequence:.recommended_experimental_sequence}' /tmp/patz1_gene_transcript_assay_routine.report.json
```

> Expected: The compound routine keeps every panel's specificity acceptance false in this offline run instead of promoting a planned follow-up to a completed check.

### Step 11: Inspect the per-assay specificity handoff templates; replace GENOME_ID and OU...

GUI: Inspect the per-assay specificity handoff templates; replace `GENOME_ID` and `OUTPUT_DIR` only when a prepared reference and an external BLAST runner are available.

CLI:

```bash
jq '.specificity_followups' docs/tutorial/generated/artifacts/patz1_transcript_assay_panels_cli/artifacts/patz1_sybr_juc_panel.report.json
```

> Expected: Every selected assay exposes a non-executing `primers specificity-plan` template, and the report honestly retains `genomic_confirmation_status: not_run` in this offline tutorial.


## Follow-up Commands

```bash
cargo run --bin gentle_cli -- --state /tmp/gentle-patz1-transcript-panels.json primers show-transcript-assay-panel patz1_endpoint_end_matrix
cargo run --bin gentle_cli -- --state /tmp/gentle-patz1-transcript-panels.json primers show-transcript-assay-panel patz1_sybr_juc_panel
cargo run --bin gentle_cli -- --state /tmp/gentle-patz1-transcript-panels.json primers show-transcript-assay-panel patz1_routine_common_region_screen
cargo run --bin gentle_cli -- --state /tmp/gentle-patz1-transcript-panels.json primers export-transcript-assay-panel patz1_endpoint_end_matrix /tmp/patz1_endpoint_end_matrix.export.json
cargo run --bin gentle_cli -- --state /tmp/gentle-patz1-transcript-panels.json primers export-transcript-assay-panel patz1_sybr_juc_panel /tmp/patz1_sybr_juc_panel.export.json
cargo run --bin gentle_cli -- --state /tmp/gentle-patz1-transcript-panels.json primers preflight --backend primer3
jq '.' /tmp/patz1_gene_transcript_assay_routine.report.json
jq '.specificity_followups' docs/tutorial/generated/artifacts/patz1_transcript_assay_panels_cli/artifacts/patz1_sybr_juc_panel.report.json
cargo run --bin gentle_examples_docs -- tutorial-check
```

## Checkpoints

- The endpoint report has `assay_kind: endpoint_rt_pcr`, `objective: isoform_end_matrix`, `completion_status: complete`, and no internal probes.
- Every endpoint reaction is designed and the band-size matrix contains at least three distinct predicted product sizes.
- The endpoint warnings explain that oligo-dT reverse-transcription completeness can limit long 5-prime products independently of PCR capacity.
- The endpoint warning describes gel intensity as rough/semi-quantitative rather than quantitative transcript abundance.
- The SYBR report has `assay_kind: sybr_qpcr`, no internal probes, and a non-empty short-junction assay table.
- At least one SYBR junction-evaluation row has `source_kind: clariom_juc` and `status: selected_spanning_assay`.
- The routine report has `assay_tier: routine_common_region_screen`; selected commonality is annotation-confirmed, and PSR/JUC evidence kinds remain distinct.
- The compound report has `schema: gentle.gene_transcript_assay_routine.v1`, matches the committed isoform-evidence digest, and assigns common-control, junction-validation, and endpoint-structure roles.
- A computed selected-pair audit names its method and target universe, keeps observation and projection counts separate, and carries the normalized panel-design operation digest; legacy missing audit fields are not interpreted as computed zeroes.
- Composition does not rerun assay design or specificity, and all three offline panel summaries retain `specificity_accepted: false`.
- Every selected assay has a `specificity-plan` follow-up template, while this offline run leaves genomic confirmation explicitly `not_run`.
- The fixture is identified as synthetic and is not presented as an order-ready human PATZ1 design.
- The same operation JSON remains callable through workflows, MCP `op`, JavaScript `apply_operation`, and Lua `apply_operation`.

## What This Chapter Produces

- [`artifacts/patz1_transcript_assay_panels_cli/artifacts/patz1_endpoint_end_matrix.report.json`](../artifacts/patz1_transcript_assay_panels_cli/artifacts/patz1_endpoint_end_matrix.report.json) - schema: `gentle.transcript_assay_panel.v2`
- [`artifacts/patz1_transcript_assay_panels_cli/artifacts/patz1_routine_common_region_screen.report.json`](../artifacts/patz1_transcript_assay_panels_cli/artifacts/patz1_routine_common_region_screen.report.json) - schema: `gentle.transcript_assay_panel.v2`
- [`artifacts/patz1_transcript_assay_panels_cli/artifacts/patz1_sybr_juc_panel.report.json`](../artifacts/patz1_transcript_assay_panels_cli/artifacts/patz1_sybr_juc_panel.report.json) - schema: `gentle.transcript_assay_panel.v2`

## Tutorial Provenance

- Chapter id: `patz1_transcript_assay_panels_cli`
- Tier: `core`
- Example id: `patz1_endpoint_sybr_transcript_assay_panel_offline`
- Tutorial source JSON: `docs/tutorial/sources/04-06_patz1_transcript_assay_panels_cli.json`
- Workflow file: `docs/examples/workflows/patz1_endpoint_sybr_transcript_assay_panel_offline.json`
- Generated artifact dir: `docs/tutorial/generated/artifacts/patz1_transcript_assay_panels_cli`
- Example test_mode: `always`
- Executed during generation: `yes`
- Automated status: `passing`
- Review status: `codex_reviewed`
- Codex reviewed at: `2026-07-21`
- Human reviewed at: `not recorded`
- Inspect the source JSON when you need full option-level detail.

## Feedback

If this tutorial is confusing, execution-stale, biologically suspect, or missing a useful figure, please open the matching tutorial issue template and include the context below.

- Tutorial title: `Design Practical PATZ1 Endpoint and SYBR Transcript Panels from the CLI`
- Tutorial/chapter id: `patz1_transcript_assay_panels_cli`
- Step reached:
- Expected vs. actual:
- Interface used: GUI / CLI / Agent Assistant / ClawBio

Paste the Tutorial feedback context here:

```text

```
