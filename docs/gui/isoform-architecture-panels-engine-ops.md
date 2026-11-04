## Isoform Architecture Panels (Engine Ops)

<!-- manual-counterparts:start -->
Related CLI reference: [Isoform Architecture and Genomic Regions](../cli/isoform-architecture.md).
These topic links do not imply identical feature coverage across interfaces.
<!-- manual-counterparts:end -->

The Engine Ops panel includes `Isoform architecture panels` for curated
transcript/protein layouts (for example TP53 panel rendering):

- `panel path`
  - path to a curated panel JSON resource
  - default: `assets/panels/tp53_isoforms_v1.json`
- `panel id`
  - optional explicit panel selector when one file contains multiple panel IDs
- `Strict`
  - when enabled, import fails on unresolved transcript mappings
  - when disabled, unresolved mappings are reported as warnings

Buttons:

- `Import Panel`
  - executes `ImportIsoformPanel` for the active sequence.
- `Open Isoform Expert`
  - opens a dedicated expert window rendered from shared engine payload
    (`FeatureExpertTarget::IsoformArchitecture`).
- `Export Isoform SVG`
  - executes `RenderIsoformArchitectureSvg` and writes a deterministic SVG.

The native expert window currently presents two scrollable tables:

- top section: transcript IDs, strand, mapping status and genomic exon coordinates
- bottom section: derived protein lengths, reference spans, translation settings,
  external opinions/status and domain coordinates
- shared isoform row ordering across both sections
- warnings from the same engine payload used by shell/CLI inspection commands.

`derived aa` reports `unavailable` when no transcript-derived product exists.
External expected lengths stay in the separate opinion column; an unattributed
legacy expectation is labelled as such, not promoted to a derivation.
Resize or scroll horizontally to inspect a wide row. The graphical exon/CDS
and protein lanes belong to the SVG export, not a native screenshot. The
[CD44 walkthrough](../tutorial/06-08_cd44_genome_transcript_protein_gui.md)
explains the distinction and the separate transcript-derived SVG route.

VCF display filtering (no recomputation needed):

- `VCF display filter` in Engine Ops offers live filtering criteria for imported
  variant features:
  - class toggles: `SNP`, `INS`, `DEL`, `SV`, `OTHER`
  - `PASS only`
  - `QUAL >= min` and `QUAL <= max` with enable checkboxes
  - required INFO keys (comma-separated)
- changes are written to shared display state and applied immediately.
- sequence SVG export uses the same VCF display filter state as the GUI display.

<!-- manual-page-nav:start -->

---

[GUI Manual](../gui.md) | [Previous: Cryptic Splicing Screen](cryptic-splicing-screen.md) | [Next: Reference Genomes (Main Menu)](reference-genomes-main-menu.md)
<!-- manual-page-nav:end -->
