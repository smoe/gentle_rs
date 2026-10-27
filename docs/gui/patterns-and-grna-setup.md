### Patterns And gRNA Setup

<!-- manual-counterparts:start -->
Related CLI reference: [Guide Design](../cli/guides.md).
These topic links do not imply identical feature coverage across interfaces.
<!-- manual-counterparts:end -->

Patterns menu:

- `Patterns -> Import Pattern File...`
  - import one workflow-macro pattern JSON file.
- `Patterns -> Import Pattern Folder...`
  - import all pattern JSON files recursively from a folder tree.
- `Patterns -> Import Built-in Legacy Pack`
  - import `assets/cloning_patterns.json`.
- `Patterns -> Import Full Catalog`
  - import the full hierarchical catalog from
    `assets/cloning_patterns_catalog`.
- `Patterns` submenu hierarchy mirrors the directory hierarchy under
  `assets/cloning_patterns_catalog` (one JSON template file per leaf entry).
  Leaves say **Import ... template**, not run/design; they remain usable in an
  empty project. Canonical routine titles preserve names such as **gRNA**;
  directory labels use **CRISPR**, **PCR**, **DNA**, **RNA**, and **TFBS**.
  The hierarchy and routine-import browser share a background-loaded snapshot.
  **Refresh template catalog** explicitly reloads it after files change; loading
  and errors remain visible without rereading directories or evaluating project
  planning on each menu frame. Shipped catalogs resolve from the checkout or
  packaged application; user-selected import paths remain explicit.
- `Patterns -> Use gRNA routines` opens **setup**, separately from template import.
  Candidate scans require explicitly selected DNA (project graph or focused
  viewer). Anchor scan setup asks for two zero-based sequence boundaries;
  `0` and sequence length are valid endpoints, the right boundary is excluded,
  and the span must fit the candidate length (20 bp by default). Reversed
  boundary order is accepted. Run Preflight remains disabled until valid.
  These are generic candidate preselection scans, not PAM-aware guide design
  or off-target confirmation. Practical Filter and Oligos requires a stored
  guide set, selected explicitly in the form; DNA or a candidate set is not a
  substitute. No guide set is silently selected by alphabetical order.
  Existing macro preflight and transactional execution remain authoritative.
- Pilot availability is shared between egui menus and Command Palette rows:
  Undo/Redo, PCR Design, Sequencing Confirmation, Saved Genomic Regions, the
  three gRNA setup actions, and palette-only TATA, transcript-start/TSS,
  cryptic-splicing, conservation, and feature-location actions. Already-open
  PCR/confirmation windows can still be focused; empty-project conservation
  setup still offers reference retrieval. Undo/Redo explain empty session history
  or the existing background-job restriction; a busy project reports Checking
  without waiting for its lock. Disabled rows retain their reason on hover and
  in the selected-action detail area. Enter and click cannot invoke them, and
  dispatch rechecks state. The gRNA submenu is
  available when any of its children is available; imports, help and configuration
  remain reachable. This is not application-wide or native-macOS-menu coverage.
- gRNA preflight is bound to its project revision, template and bindings.
  Editing bindings, switching/replacing a project, or changing project state
  requires another preflight before Run. Catalog checks run outside menu paint
  and refresh explicitly; availability does not grant execution permission.

<!-- manual-page-nav:start -->

---

[GUI Manual](../gui.md) | [Previous: Protein Expression Handoff](protein-expression-handoff.md) | [Next: Routine Assistant](routine-assistant.md)
<!-- manual-page-nav:end -->
