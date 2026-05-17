---
chapter_id: "tp73_uniprot_projection_audit_cli"
title: "Audit a TP73 UniProt Projection Against Ensembl and Derived Coding Sequence (CLI Tutorial)"
tier: "online"
example_id: "tp73_uniprot_projection_audit_online"
source_example: "docs/examples/workflows/tp73_uniprot_projection_audit_online.json"
example_test_mode: "online"
executed_during_generation: false
---

# Audit a TP73 UniProt Projection Against Ensembl and Derived Coding Sequence (CLI Tutorial)

Build a TP73 locus project, fetch UniProt and Ensembl protein evidence, persist the integrated audit plus parity report, and keep the matching expert SVG as an inspectable artifact.

This executable chapter turns the TP73 UniProt/Ensembl audit into a reproducible online workflow. It prepares GRCh38, extracts TP73, fetches the reviewed UniProt entry and matching Ensembl protein evidence, projects the protein onto the locus, persists the integrated audit and direct-vs-composed parity report, and keeps one shared expert SVG export alongside the stored report artifacts. The companion hand-written tutorial then explains how to inspect the same outcome through the public primitive commands.

See also: guided walkthrough [docs/tutorial/tp73_uniprot_projection_audit_cli.md](../../tp73_uniprot_projection_audit_cli.md). Use that page first when you want a human-led path; this chapter is the executable reference.

> **How to Run This Locally**
> Set `GENTLE_TEST_ONLINE=1` and run from the repository root. This workflow prepares/extracts GRCh38 Ensembl 116 from Ensembl FTP, fetches UniProt `Q9H3D4`, queries Ensembl protein evidence for `ENSP00000264724`, and then writes the audit/parity reports plus SVG locally.

## When This Routine Is Useful

- You want one canonical online starter project that proves the TP73 UniProt/Ensembl audit path still runs end to end.
- You want a persisted TP73 audit report plus parity report before following the shell-level primitive walkthrough.
- You want a deterministic bridge from locus extraction and UniProt projection into audit, parity, and expert-view artifact export.

## What You Learn

- Use one executable chapter as the reproducible setup layer behind the TP73 UniProt/Ensembl audit walkthrough.
- Understand that the high-level audit and parity report are built on the same reusable primitive record families exposed to shell/CLI and future AI callers.
- Preserve both the persisted audit reports and one shared expert SVG artifact so the TP73 workflow can be inspected visually and programmatically.

## Concepts

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

1. Prepare `Human GRCh38 Ensembl 116` and extract gene `TP73` into `grch38_tp73`.
2. Fetch UniProt `Q9H3D4` and Ensembl protein `ENSP00000264724` from `Protein Evidence...`, then project the UniProt entry onto the TP73 locus.
3. Run the high-level audit and parity actions from `Protein Evidence...` so the stored audit rows and local unsent email draft are persisted.
4. Open the saved projection in the Protein Expert or inspect the exported SVG artifact to verify the projected feature geometry.
5. Use the companion CLI tutorial to rebuild the same result from `resolve-ensembl-links`, `transcript-accounting`, `compare-ensembl-exons`, and `compare-ensembl-peptide`.

## Command Equivalent (After GUI)

Run the same routine non-interactively once the GUI flow is clear:

```bash
cargo run --bin gentle_cli -- workflow @docs/examples/workflows/tp73_uniprot_projection_audit_online.json
cargo run --bin gentle_cli -- shell 'workflow @docs/examples/workflows/tp73_uniprot_projection_audit_online.json'
```

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

## Canonical Source

- Chapter id: `tp73_uniprot_projection_audit_cli`
- Tier: `online`
- Example id: `tp73_uniprot_projection_audit_online`
- Workflow file: `docs/examples/workflows/tp73_uniprot_projection_audit_online.json`
- Example test_mode: `online`
- Executed during generation: `no`
- Execution note: set `GENTLE_TEST_ONLINE=1` before `tutorial-generate` to execute this chapter.
- Inspect this JSON file directly when you need full option-level detail.
