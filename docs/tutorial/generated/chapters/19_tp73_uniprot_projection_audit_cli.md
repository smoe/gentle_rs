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
codex_reviewed_at: null
human_reviewed_at: null
generated_artifact_dir: "docs/tutorial/generated/artifacts/tp73_uniprot_projection_audit_cli"
---

# Audit a TP73 UniProt Projection Against Ensembl and Derived Coding Sequence (CLI Tutorial)

Build a TP73 locus project, fetch UniProt and Ensembl protein evidence, persist the integrated audit plus parity report, and keep the matching expert SVG as an inspectable artifact.

This executable chapter turns the TP73 UniProt/Ensembl audit into a reproducible online workflow. It prepares GRCh38, extracts TP73, fetches the reviewed UniProt entry and matching Ensembl protein evidence, projects the protein onto the locus, persists the integrated audit and direct-vs-composed parity report, and keeps one shared expert SVG export alongside the stored report artifacts. The companion hand-written tutorial then explains how to inspect the same outcome through the public primitive commands.

See also: guided walkthrough [docs/tutorial/tp73_uniprot_projection_audit_cli.md](../../tp73_uniprot_projection_audit_cli.md). Use that page first when you want a human-led path; this chapter is the executable reference.

**Prerequisites:** Read [Chapter 17: TP53 UniProt domain mapping and feature-coding DNA query (online)](./17_tp53_uniprot_projection_online.md) first.

> **How to Run This Locally**
> Set `GENTLE_TEST_ONLINE=1` and run from the repository root. This workflow prepares/extracts GRCh38 Ensembl 116 from Ensembl FTP, fetches UniProt `Q9H3D4`, queries Ensembl protein evidence for `ENSP00000264724`, and then writes the audit/parity reports plus SVG locally.

## Parameters That Matter

- `FetchUniprotSwissProt.query / entry_id` (where used: operation 3)
  - Why it matters: Pins the reviewed UniProt TP73 entry that drives the protein-coordinate system for projection and audit.
  - How to derive it: Use the stable reviewed TP73 accession `Q9H3D4` for this canonical chapter.
- `FetchEnsemblProtein.query / entry_id` (where used: operation 4)
  - Why it matters: Keeps the Ensembl protein evidence explicit and reusable for exon/CDS and peptide comparison during audit.
  - How to derive it: Use `ENSP00000264724` with a stable project-local id such as `TP73_ENS`.
- `ProjectUniprotToGenome.projection_id` (where used: operation 5)
  - Why it matters: The stored projection id is the durable handle reused by the expert SVG export, the integrated audit, and the parity report.
  - How to derive it: Use a stable intent-bearing id such as `tp73_uniprot_q9h3d4`.
- `AuditUniprotProjectionConsistency.report_id / AuditUniprotProjectionParity.report_id` (where used: operations 7 and 8)
  - Why it matters: Stable report ids make the saved audit/parity artifacts easy to reopen from GUI, CLI, or future AI orchestration.
  - How to derive it: Use descriptive ids like `tp73_projection_audit` and `tp73_projection_audit_parity`.

## When This Routine Is Useful

- You want one canonical online starter project that proves the TP73 UniProt/Ensembl audit path still runs end to end.
- You want a persisted TP73 audit report plus parity report before following the shell-level primitive walkthrough.
- You want a deterministic bridge from locus extraction and UniProt projection into audit, parity, and expert-view artifact export.

## What You Learn

- Use one executable chapter as the reproducible setup layer behind the TP73 UniProt/Ensembl audit walkthrough.
- Understand that the high-level audit and parity report are built on the same reusable primitive record families exposed to shell/CLI and future AI callers.
- Preserve both the persisted audit reports and one shared expert SVG artifact so the TP73 workflow can be inspected visually and programmatically.

## Applied Concepts

- **Shared Engine Contract** (`shared_engine_contract`): GUI, CLI, shell, and scripting interfaces execute the same operation semantics.
- **UniProt Projection Mapping** (`uniprot_projection_mapping`): Reviewed UniProt protein annotations can be projected onto transcript/CDS geometry and reopened as one persisted expert-view artifact.
- **Feature Coding-DNA Attribution** (`feature_coding_dna_attribution`): A persisted UniProt projection can be queried for the exact coding DNA, exon attribution, and splice-junction exon pairs that encode one mapped protein feature.
- **Expert View Parity** (`expert_view_parity`): The same expert-view payloads should be inspectable and renderable from GUI, CLI, and other adapters without frontend-only projection logic.
- **Online Opt-in Execution** (`online_opt_in`): Network-dependent chapters remain explicit opt-in and do not break offline default CI.
- **Artifact Exports** (`artifact_exports`): Representative outputs (CSV/protocol/SVG/text) are retained for auditability and sharing.

## At a Glance

1. Prepare Human GRCh38 Ensembl 116 and extract gene TP73 into grch38_tp73.
2. Fetch UniProt Q9H3D4 and Ensembl protein ENSP00000264724 from Protein Evidenc...
3. Run the high-level audit and parity actions from Protein Evidence... so the s...
4. Open the saved projection in the Protein Expert or inspect the exported SVG a...
5. Use the companion CLI tutorial to rebuild the same result from resolve-ensemb...

## GUI First

CLI snippets use GENtle's default `.gentle_state.json` state unless they say otherwise. Add `--state PATH` or `--project PATH` when you want an explicit sandboxed state file for copied commands.

### Step 1: Prepare Human GRCh38 Ensembl 116 and extract gene TP73 into grch38_tp73

GUI: Prepare `Human GRCh38 Ensembl 116` and extract gene `TP73` into `grch38_tp73`.

CLI:

```bash
GENTLE_TEST_ONLINE=1 cargo run --bin gentle_cli -- genomes prepare "Human GRCh38 Ensembl 116" --catalog assets/genomes.json --cache-dir data/genomes --timeout-secs 3600
cargo run --bin gentle_cli -- genomes extract-gene "Human GRCh38 Ensembl 116" TP73 --occurrence 1 --output-id grch38_tp73 --catalog assets/genomes.json --cache-dir data/genomes
```

> Expected: The reference is prepared if needed and TP73 is extracted into the anchored sequence id `grch38_tp73`.

### Step 2: Fetch UniProt Q9H3D4 and Ensembl protein ENSP00000264724 from Protein Evidenc...

GUI: Fetch UniProt `Q9H3D4` and Ensembl protein `ENSP00000264724` from `Protein Evidence...`, then project the UniProt entry onto the TP73 locus.

CLI:

```bash
cargo run --bin gentle_cli -- shell 'uniprot fetch Q9H3D4 --entry-id Q9H3D4'
cargo run --bin gentle_cli -- shell 'ensembl-protein fetch ENSP00000264724 --entry-id TP73_ENS'
cargo run --bin gentle_cli -- shell 'uniprot map Q9H3D4 grch38_tp73 --projection-id tp73_uniprot_q9h3d4'
```

> Expected: The UniProt, Ensembl protein, and projection records persist under the ids used by the audit workflow.

### Step 3: Run the high-level audit and parity actions from Protein Evidence... so the s...

GUI: Run the high-level audit and parity actions from `Protein Evidence...` so the stored audit rows and local unsent email draft are persisted.

CLI:

```bash
cargo run --bin gentle_cli -- shell 'uniprot audit-projection tp73_uniprot_q9h3d4 --ensembl-entry TP73_ENS --report-id tp73_projection_audit'
cargo run --bin gentle_cli -- shell 'uniprot audit-parity tp73_uniprot_q9h3d4 --ensembl-entry TP73_ENS --report-id tp73_projection_audit_parity'
```

> Expected: The integrated audit and parity reports are stored as `tp73_projection_audit` and `tp73_projection_audit_parity`.

### Step 4: Open the saved projection in the Protein Expert or inspect the exported SVG a...

GUI: Open the saved projection in the Protein Expert or inspect the exported SVG artifact to verify the projected feature geometry.

CLI:

```bash
cargo run --bin gentle_cli -- inspect-feature-expert grch38_tp73 uniprot-projection tp73_uniprot_q9h3d4
cargo run --bin gentle_cli -- render-feature-expert-svg grch38_tp73 uniprot-projection tp73_uniprot_q9h3d4 exports/tp73_uniprot_projection.svg
```

> Expected: Expert inspection and SVG export show the same projected TP73 feature geometry used by the GUI.

### Step 5: Use the companion CLI tutorial to rebuild the same result from resolve-ensemb...

GUI: Use the companion CLI tutorial to rebuild the same result from `resolve-ensembl-links`, `transcript-accounting`, `compare-ensembl-exons`, and `compare-ensembl-peptide`.

CLI:

```bash
cargo run --bin gentle_cli -- shell 'uniprot resolve-ensembl-links tp73_uniprot_q9h3d4'
cargo run --bin gentle_cli -- shell 'uniprot transcript-accounting tp73_uniprot_q9h3d4'
cargo run --bin gentle_cli -- shell 'uniprot compare-ensembl-exons tp73_uniprot_q9h3d4 --ensembl-entry TP73_ENS'
cargo run --bin gentle_cli -- shell 'uniprot compare-ensembl-peptide tp73_uniprot_q9h3d4 --ensembl-entry TP73_ENS'
```

> Expected: The primitive CLI commands rebuild the audit evidence path behind the integrated reports.


## Follow-up Commands

```bash
cargo run --bin gentle_cli -- shell 'uniprot audit-show tp73_projection_audit'
cargo run --bin gentle_cli -- shell 'uniprot resolve-ensembl-links tp73_uniprot_q9h3d4'
cargo run --bin gentle_cli -- shell 'uniprot transcript-accounting tp73_uniprot_q9h3d4'
cargo run --bin gentle_cli -- shell 'uniprot compare-ensembl-exons tp73_uniprot_q9h3d4 --ensembl-entry TP73_ENS'
cargo run --bin gentle_cli -- shell 'uniprot compare-ensembl-peptide tp73_uniprot_q9h3d4 --ensembl-entry TP73_ENS'
cargo run --bin gentle_cli -- shell 'uniprot audit-parity-show tp73_projection_audit_parity'
```

## Checkpoints

- The TP73 locus, UniProt entry, and Ensembl protein evidence all persist under stable ids in one project state.
- The integrated audit stores `tp73_projection_audit` with per-transcript accounting and a local unsent maintainer-email draft.
- The parity report stores `tp73_projection_audit_parity` and confirms whether the integrated Rust audit matches the public primitive composition.
- The shared expert SVG export succeeds so the projected TP73 feature geometry remains inspectable outside the live GUI.

## Tutorial Provenance

- Chapter id: `tp73_uniprot_projection_audit_cli`
- Tier: `online`
- Example id: `tp73_uniprot_projection_audit_online`
- Tutorial source JSON: `docs/tutorial/sources/30_tp73_uniprot_projection_audit_cli.json`
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
