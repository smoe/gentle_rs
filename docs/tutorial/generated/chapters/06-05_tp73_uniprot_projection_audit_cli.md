---
chapter_id: "tp73_uniprot_projection_audit_cli"
title: "Audit a TP73 UniProt Projection Against Ensembl and Derived Coding Sequence (CLI Tutorial)"
tier: "online"
example_id: "tp73_uniprot_projection_audit_online"
source_example: "docs/examples/workflows/tp73_uniprot_projection_audit_online.json"
example_test_mode: "online"
executed_during_generation: false
automated_status: "skipped_online"
review_status: "unreviewed"
review_stale: false
codex_reviewed_at: null
human_reviewed_at: null
human_reviewer: null
review_stale_reason: null
review_issue_template: "Tutorial confusion"
review_issue_template_path: ".github/ISSUE_TEMPLATE/tutorial-confusion.md"
generated_artifact_dir: "docs/tutorial/generated/artifacts/tp73_uniprot_projection_audit_cli"
---

# Audit a TP73 UniProt Projection Against Ensembl and Derived Coding Sequence (CLI Tutorial)

Build a fully annotated TP73 locus, fetch reviewed UniProt O15350, persist the integrated audit plus parity report, and keep the matching expert SVG as an inspectable artifact.

This executable chapter turns the TP73 UniProt audit into a reproducible online workflow. It prepares GRCh38, extracts TP73 with full exon/CDS annotation, fetches reviewed UniProt O15350, projects p73 onto the locus, persists the integrated audit and direct-vs-composed parity report, and keeps one shared expert SVG. The companion guide explains the public primitives, the optional Ensembl-protein comparison, and why missing evidence must not be reported as a UniProt error.

See also: guided walkthrough [docs/tutorial/06-05_tp73_uniprot_projection_audit_cli.md](../../06-05_tp73_uniprot_projection_audit_cli.md). Use that page first when you want a human-led path; this chapter is the executable reference.

## What You Will Accomplish

- Use one executable chapter as the reproducible setup layer behind the TP73 UniProt/Ensembl audit walkthrough.
- Understand that the high-level audit and parity report are built on the same reusable primitive record families exposed to shell/CLI and future AI callers.
- Preserve both the persisted audit reports and one shared expert SVG artifact so the TP73 workflow can be inspected visually and programmatically.

## Before You Start

**Prerequisites:** Read [Chapter 17: TP53 UniProt domain mapping and feature-coding DNA query (online)](./06-04_tp53_uniprot_projection_online.md) first.

> **How to Run This Locally**
> Set `GENTLE_TEST_ONLINE=1` and run from the repository root. The workflow prepares GRCh38 Ensembl 116, extracts TP73 with full exon/CDS annotation, fetches reviewed UniProt `O15350`, and writes audit/parity reports plus SVG locally. A separate Ensembl-protein fetch is optional; during the 2026-10-01 review Ensembl lookup and sequence endpoints returned HTTP 500, so no external sequence was fabricated.

**Useful when:**

- You want one canonical online starter project that proves the TP73 UniProt/Ensembl audit path still runs end to end.
- You want a persisted TP73 audit report plus parity report before following the shell-level primitive walkthrough.
- You want a deterministic bridge from locus extraction and UniProt projection into audit, parity, and expert-view artifact export.

## At a Glance

1. Prepare Human GRCh38 Ensembl 116 and extract gene TP73 into grch38_tp73 with annotation scope full.
2. Fetch reviewed UniProt O15350 as TP73_UNIPROT from Protein Evidence..., then project it onto the TP73 locus as tp73_uniprot_o15350.
3. Run the high-level audit from Protein Evidence...; distinguish protein mismatch from missing_evidence, and keep the maintainer-email text local and unsent.
4. Run and reopen the parity report; verify that the native summary reports 0 / 10 divergent rows before trusting the outer primitive composition.
5. Open the saved projection in the Protein Expert or inspect the exported SVG artifact to verify the projected feature geometry.
6. Use the companion CLI tutorial to rebuild the same result from resolve-ensembl-links, transcript-accounting, compare-ensembl-exons, and compare-ensembl-peptide; fetch optional Ensembl protein evidence only when its REST service is healthy.

## Walkthrough: GUI, CLI and Inner Agent

Each step pairs GUI instructions with related terminal commands or guidance and a review-only inner-agent prompt. Some steps require GUI interaction; a listing command only inspects results, it does not perform the design. CLI snippets assume an installed `gentle_cli` and use GENtle's default `.gentle_state.json` unless stated otherwise. From a source checkout, replace `gentle_cli` with `cargo run --bin gentle_cli --`. Add `--state PATH` or `--project PATH` for an explicit sandbox; a separate CLI process does not inherit the open GUI project's unsaved state.

In the **GUI Shell**, enter only the shared command inside `gentle_cli shell '...'`, without the executable prefix or outer quotes. Run UI-opening commands there to open windows: a headless CLI returns the UI intent but does not open a GUI. Other terminal commands are not automatically GUI Shell commands. Inner-agent examples request a proposal for review; they are not executed during tutorial generation.

### Step 1: Prepare Human GRCh38 Ensembl 116 and extract gene TP73 into grch38_tp73 with annotation scope full

**GUI**

Prepare `Human GRCh38 Ensembl 116` and extract gene `TP73` into `grch38_tp73` with annotation scope `full`.

**CLI (terminal)**

```bash
GENTLE_TEST_ONLINE=1 gentle_cli genomes prepare "Human GRCh38 Ensembl 116" --catalog assets/genomes.json --cache-dir data/genomes --timeout-secs 7200
gentle_cli genomes extract-gene "Human GRCh38 Ensembl 116" TP73 --occurrence 1 --output-id grch38_tp73 --annotation-scope full --catalog assets/genomes.json --cache-dir data/genomes
```

**Ask the inner agent**

> Read `assets/genomes.json`, preflight exactly `Human GRCh38 Ensembl 116`, and propose preparation plus TP73 occurrence-1 extraction with `annotation_scope=full`. Expose the Ensembl URLs, cache writes, timeout and state mutation separately; core annotation is insufficient for this audit.

**Expected**

> The reference is prepared if needed and TP73 is extracted into `grch38_tp73` with non-zero transcript, exon and CDS feature counts.

### Step 2: Fetch reviewed UniProt O15350 as TP73_UNIPROT from Protein Evidence..., then project it onto the TP73 locus as tp73_uniprot_o15350

**GUI**

Fetch reviewed UniProt `O15350` as `TP73_UNIPROT` from `Protein Evidence...`, then project it onto the TP73 locus as `tp73_uniprot_o15350`.

**CLI (terminal)**

```bash
gentle_cli shell 'uniprot fetch O15350 --entry-id TP73_UNIPROT'
gentle_cli shell 'uniprot map TP73_UNIPROT grch38_tp73 --projection-id tp73_uniprot_o15350'
```

**Ask the inner agent**

> Verify gene symbol TP73 and reviewed UniProt identity `O15350` / `P73_HUMAN` before proposing fetch and projection. Reject `Q9H3D4` because it is TP63, then report the fetched source URL and ten projected transcript ids.

**Expected**

> Reviewed UniProt O15350 persists as `TP73_UNIPROT`, and its ten linked transcript projections persist as `tp73_uniprot_o15350`.

### Step 3: Run the high-level audit from Protein Evidence...; distinguish protein mismatch from missing_evidence, and keep the maintainer-email text local and unsent

**GUI**

Run the high-level audit from `Protein Evidence...`; distinguish protein mismatch from `missing_evidence`, and keep the maintainer-email text local and unsent.

**CLI (terminal)**

```bash
gentle_cli shell 'uniprot audit-projection tp73_uniprot_o15350 --report-id tp73_projection_audit'
```

**Ask the inner agent**

> Run the stored audit only after confirming every accounting row satisfies translated_nt/3=derived_aa. Separate six protein-length mismatches from four rows limited by missing Ensembl evidence, and never send or present the local draft as a validated complaint.

**Expected**

> The integrated audit is stored as `tp73_projection_audit`; all ten accounting rows use genomic CDS coordinates and separate mismatch from missing evidence.

![The native Protein Evidence audit separates six protein-length mismatches from four rows whose remaining limitation is missing Ensembl evidence. Expanded rows show corrected genomic CDS accounting, and the text box is explicitly an unsent local draft.](../../../screenshots/tp73_uniprot_projection_audit_online/01-audit-mismatch-vs-missing-evidence.png)

*Figure: The native Protein Evidence audit separates six protein-length mismatches from four rows whose remaining limitation is missing Ensembl evidence. Expanded rows show corrected genomic CDS accounting, and the text box is explicitly an unsent local draft. Screenshot captured 2026-10-01.*

### Step 4: Run and reopen the parity report; verify that the native summary reports 0 / 10 divergent rows before trusting the outer primitive composition

**GUI**

Run and reopen the parity report; verify that the native summary reports `0 / 10` divergent rows before trusting the outer primitive composition.

**CLI (terminal)**

```bash
gentle_cli shell 'uniprot audit-parity tp73_uniprot_o15350 --report-id tp73_projection_audit_parity'
```

**Ask the inner agent**

> Run direct-vs-composed parity and require zero divergent rows across status, accounting, mismatch reasons, comparison mode and draft transcript set. Preserve the structured parity payload instead of replacing it with prose.

**Expected**

> The parity report is stored as `tp73_projection_audit_parity` and reports zero divergent rows out of ten.

![The native parity summary reports zero divergent rows out of ten and confirms that the direct and composed paths selected the same transcript set for the local draft.](../../../screenshots/tp73_uniprot_projection_audit_online/02-direct-vs-composed-parity.png)

*Figure: The native parity summary reports zero divergent rows out of ten and confirms that the direct and composed paths selected the same transcript set for the local draft. Screenshot captured 2026-10-01.*

### Step 5: Open the saved projection in the Protein Expert or inspect the exported SVG artifact to verify the projected feature geometry

**GUI**

Open the saved projection in the Protein Expert or inspect the exported SVG artifact to verify the projected feature geometry.

**CLI (terminal)**

```bash
gentle_cli inspect-feature-expert grch38_tp73 uniprot-projection tp73_uniprot_o15350
gentle_cli render-feature-expert-svg grch38_tp73 uniprot-projection tp73_uniprot_o15350 exports/tp73_uniprot_projection.svg
```

**Ask the inner agent**

> Inspect or render the persisted projection without claiming that an SVG or screenshot supplies missing Ensembl evidence.

**Expected**

> Expert inspection and SVG export show the same projected TP73 feature geometry used by the GUI.

### Step 6: Use the companion CLI tutorial to rebuild the same result from resolve-ensembl-links, transcript-accounting, compare-ensembl-exons, and compare-ensembl-peptide; fetch optional Ensembl protein evidence only when its REST service is healthy

**GUI**

Use the companion CLI tutorial to rebuild the same result from `resolve-ensembl-links`, `transcript-accounting`, `compare-ensembl-exons`, and `compare-ensembl-peptide`; fetch optional Ensembl protein evidence only when its REST service is healthy.

**CLI (terminal)**

```bash
gentle_cli shell 'uniprot resolve-ensembl-links tp73_uniprot_o15350'
gentle_cli shell 'uniprot transcript-accounting tp73_uniprot_o15350'
gentle_cli shell 'uniprot compare-ensembl-exons tp73_uniprot_o15350'
gentle_cli shell 'uniprot compare-ensembl-peptide tp73_uniprot_o15350'
```

**Ask the inner agent**

> Rebuild the report from the four public primitives. Treat an optional Ensembl protein fetch as a separate network action and retain `missing_evidence` when the provider fails.

**Expected**

> The primitive CLI commands rebuild the audit evidence path behind the integrated reports.


## Ask an Outer Agent (MCP or ClawBio/OpenClaw)

An outer agent does not inherit the unsaved GUI project. Give it this chapter's canonical workflow and an explicit disposable state path; ask it to retain the structured result, artifacts and reproducibility receipt instead of replacing them with prose.

> Use GENtle's `gentle-cloning` skill to replay `docs/examples/workflows/tp73_uniprot_projection_audit_online.json` against a new disposable state. First report the exact workflow, inputs, state path, outputs and whether the selected route needs confirmation. Do not infer state from an open GUI. Return the structured result, produced artifacts and reproducibility receipt, and state any unmet prerequisite.

Equivalent direct structured request for the generic wrapper:

```json
{
  "schema": "gentle.clawbio_skill_request.v1",
  "mode": "workflow",
  "state_path": "/tmp/gentle-tp73-uniprot-projection-audit-cli.state.json",
  "workflow_path": "docs/examples/workflows/tp73_uniprot_projection_audit_online.json",
  "timeout_secs": 7200
}
```

Submitting a direct structured request is an explicit wrapper invocation. If natural language selects a narrower delegated skill and that route mutates state, selects biological material or writes artifacts, the caller must preserve that skill's proposal/approval boundary and approve only the exact bound digest. This tutorial generation step does not invoke an agent or grant approval.

## Interpretation and Reference

## Parameters That Matter

- `ExtractGenomeGene.annotation_scope` (where used: operation 2)
  - Why it matters: The audit intersects transcript exons with genomic CDS ranges; core scope omits those exon/CDS features and makes the accounting uninterpretable.
  - How to derive it: Use `full` for this audit and verify that the extraction result reports non-zero exon and CDS counts.
- `FetchUniprotSwissProt.query / entry_id` (where used: operation 3)
  - Why it matters: Pins the reviewed UniProt TP73 entry that drives the protein-coordinate system for projection and audit.
  - How to derive it: Use reviewed TP73 accession `O15350` (`P73_HUMAN`). `Q9H3D4` is TP63 and is biologically invalid here.
- `ProjectUniprotToGenome.projection_id` (where used: operation 4)
  - Why it matters: The stored projection id is the durable handle reused by the expert SVG export, the integrated audit, and the parity report.
  - How to derive it: Use a stable intent-bearing id such as `tp73_uniprot_o15350`.
- `AuditUniprotProjectionConsistency.report_id / AuditUniprotProjectionParity.report_id` (where used: operations 6 and 7)
  - Why it matters: Stable report ids make the saved audit/parity artifacts easy to reopen from GUI, CLI, or future AI orchestration.
  - How to derive it: Use descriptive ids like `tp73_projection_audit` and `tp73_projection_audit_parity`.

## Applied Concepts

- **Shared Engine Contract** (`shared_engine_contract`): GUI, CLI, shell, and scripting interfaces execute the same operation semantics.
- **UniProt Projection Mapping** (`uniprot_projection_mapping`): Reviewed UniProt protein annotations can be projected onto transcript/CDS geometry and reopened as one persisted expert-view artifact.
- **Feature Coding-DNA Attribution** (`feature_coding_dna_attribution`): A persisted UniProt projection can be queried for the exact coding DNA, exon attribution, and splice-junction exon pairs that encode one mapped protein feature.
- **Expert View Parity** (`expert_view_parity`): The same expert-view payloads should be inspectable and renderable from GUI, CLI, and other adapters without frontend-only projection logic.
- **Online Opt-in Execution** (`online_opt_in`): Network-dependent chapters remain explicit opt-in and do not break offline default CI.
- **Artifact Exports** (`artifact_exports`): Representative outputs (CSV/protocol/SVG/text) are retained for auditability and sharing.

## Follow-up Commands

```bash
gentle_cli shell 'uniprot audit-show tp73_projection_audit'
gentle_cli shell 'uniprot resolve-ensembl-links tp73_uniprot_o15350'
gentle_cli shell 'uniprot transcript-accounting tp73_uniprot_o15350'
gentle_cli shell 'uniprot compare-ensembl-exons tp73_uniprot_o15350'
gentle_cli shell 'uniprot compare-ensembl-peptide tp73_uniprot_o15350'
gentle_cli shell 'uniprot audit-parity-show tp73_projection_audit_parity'
```

## Checkpoints

- The Ensembl-116 TP73 locus carries full annotation, reviewed UniProt O15350 persists as `TP73_UNIPROT`, and ten linked transcripts persist in `tp73_uniprot_o15350`.
- Every accounting row satisfies translated_nt/3=derived_aa; the integrated audit separates six protein-length mismatches from four otherwise matching rows limited by missing Ensembl evidence.
- The local maintainer-email text remains explicitly unsent and is never treated as a validated complaint while external evidence is missing.
- The parity report stores `tp73_projection_audit_parity` with zero divergent rows out of ten and a matching draft transcript set.
- The shared expert SVG export succeeds so the projected TP73 feature geometry remains inspectable outside the live GUI.

## Tutorial Provenance

- Chapter id: `tp73_uniprot_projection_audit_cli`
- Tier: `online`
- Example id: `tp73_uniprot_projection_audit_online`
- Tutorial source JSON: `docs/tutorial/sources/06-05_tp73_uniprot_projection_audit_cli.json`
- Workflow file: `docs/examples/workflows/tp73_uniprot_projection_audit_online.json`
- Generated artifact dir: `docs/tutorial/generated/artifacts/tp73_uniprot_projection_audit_cli`
- Example test_mode: `online`
- Executed during generation: `no`
- Automated status: `skipped_online`
- Review status: `unreviewed`
- Codex reviewed at: `not recorded`
- Human reviewed at: `not recorded`
- Execution note: set `GENTLE_TEST_ONLINE=1` before `tutorial-generate` to execute this chapter.
- Inspect the source JSON when you need full option-level detail.

## Feedback

If this tutorial is confusing, execution-stale, biologically suspect, or missing a useful figure, please open the matching tutorial issue template and include the context below.

- Tutorial title: `Audit a TP73 UniProt Projection Against Ensembl and Derived Coding Sequence (CLI Tutorial)`
- Tutorial/chapter id: `tp73_uniprot_projection_audit_cli`
- Step reached:
- Expected vs. actual:
- Interface used: GUI / CLI / Agent Assistant / ClawBio

Paste the Tutorial feedback context here:

```text

```
