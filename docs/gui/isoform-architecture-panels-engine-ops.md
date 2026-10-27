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

Expected output structure:

- top section: transcript/exon lanes on genomic coordinates
- bottom section: protein/domain lanes on amino-acid coordinates
- shared isoform row ordering across both sections
- warning/instruction footer from the same engine payload used by shell/CLI
  inspection commands.

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
