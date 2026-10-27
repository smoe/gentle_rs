### Gibson Workspace

<!-- manual-counterparts:start -->
No direct CLI counterpart: This chapter describes desktop presentation or interaction; no dedicated command reference is paired with it.
<!-- manual-counterparts:end -->

- `Patterns -> Gibson...`
  - opens a dedicated destination-first Gibson specialist window for one or
    more inserts into one destination.
  - sections:
    1. `Destination`
       - choose destination sequence
       - choose `existing termini` or `defined cut/opening`
       - active DNA window pre-fills the destination by default when possible
       - optionally fill left/right cut edges from the active DNA selection
       - suggests biologically meaningful openings from MCS annotations and
         unique restriction sites when the destination already exposes them
       - default view starts with unique cutters named by the MCS annotation
         and marks them as `(MCS)`
       - `Show other unique cutters` reveals single-cutters elsewhere on the
         destination when you want to search beyond the MCS
       - when no MCS is present or no suggestions are found, `Search unique
         cleavage sites` explicitly re-scans all known restriction enzymes and
         reports when none are unique on the current destination
       - `Feature` column shows the overlapping gene when available, otherwise
         another overlapping feature label; MCS-linked rows also show the MCS
         location there
       - suggestion cells use compact wrapped columns so longer feature or
         REBASE text does not force the Gibson window excessively wide
       - shows end geometry (`blunt`, `5' overhang`, `3' overhang`) so sites
         such as `SmaI` can be chosen directly from the specialist
       - specific cutter suggestions show a compact `Cut` summary derived from
         REBASE and fill the actual cleavage window:
         - equal left/right cut edges mean a blunt cutpoint
         - different left/right cut edges mean the two primer-relevant termini
           of a sticky cut
       - `Opening sketch` renders the destination sequence at the exact cut or
         opening and then shows the destination after it has been opened into
         separate left/right arms
       - once preview is available, the same sketch also shows the two resolved
         destination 5' ends with their current Gibson overlap sizes and
         sequences
       - tooltips and quick help explain the difference between MCS-linked
         cutters and other unique cutters outside the annotated MCS
    2. `Inserts`
       - choose one or more insert sequences and orientations
       - insert rows are ordered explicitly in the GUI
       - one insert creates two terminal junctions; each additional insert
         adds one internal Gibson junction
       - current execution guardrail:
         multi-insert Gibson currently requires a defined destination opening;
         `existing_termini` remains the single-fragment handoff path
       - extra insert rows can be added, removed, and reordered directly in
         the specialist
    3. `Design Targets`
       - Gibson-specific overlap bp range
       - when overlap length is derived rather than fixed, GENtle prefers the
         shortest acceptable overlap inside that range
       - minimum overlap Tm
       - priming-segment Tm window
       - priming-segment length window
       - optional `new unique site` REBASE enzyme name:
         request one new unique restriction site on the assembled product by
         mutating one terminal overlap when possible
       - current v1 limit for that unique-site request:
         defined-site, single-insert plans with palindromic cutters
       - dedicated `Tm Model` box repeats the shared assumptions in a visually
         separate section
       - displayed Tm values use the shared Thermo Fisher-style modified
         Allawi/SantaLucia nearest-neighbor model used across GENtle:
         exact-complement duplex assumption, Allawi/SantaLucia 1997
         Watson-Crick table and terminal initiation terms, fixed 215 mM
         effective monovalent salt term, fixed 500 nM primer concentration,
         Thermo Fisher high-fidelity empirical adjustment, no
         mismatch/dangling-end/Mg correction, and fallback to the simple 2/4
         estimate for ambiguous or very short sequences
    4. `Review`
       - `Preview Gibson Plan` stays disabled until the current opening is
         actually defined enough to build a plan; for `defined opening`, that
         means both cut-edge fields must be set (or filled from a cutter
         suggestion / active selection)
       - resolves the ordered Gibson junction chain (`n + 1` junctions for
         `n` inserts)
       - derives one left and one right insert primer per fragment as
         `5' overlap + 3' priming`
       - `Opening sketch` extends into an insert-primer construction view so
         the actual overlap segment and insert-priming segment are both shown
       - `Target review` now separates overlap success from PCR-priming
         success, so it is explicit when the 5' Gibson overlaps already
         resolve cleanly and the remaining blocker is only the 3' gene-
         specific priming window
       - shows blocking errors, warnings, and the factual Gibson cartoon
       - priming-segment failures now explain the best available insert-end
         candidate explicitly, including when the insert terminus is simply too
         short / too low-`Tm` for the requested Gibson primer window
       - when overlap derivation succeeds but PCR priming still fails, the
         review now exposes structured next-step relaxations with one-click
         apply buttons:
         increase max priming length and/or lower the minimum priming `Tm`,
         then rerun the shared preview immediately
       - when a `new unique site` request is present, the review shows whether
         the product already had one unique site or which terminal overlap was
         engineered, including the enzyme name, overlap sequence, motif offset,
         and mutation count
       - visible findings use styled `Tₘ` rendering; the copyable findings box
         keeps plain `Tm` text for reliable copying
       - findings, resolved junctions, primer suggestions, and status text are
         copyable directly from the window
       - the in-window cartoon preview is rasterized from the same resolved
         deterministic SVG payload; when a platform-specific widget issue still
         prevents inline display, the exported SVG plus the textual review
         blocks remain the canonical inspection path
       - single-insert previews show both destination-insert junctions
         explicitly instead of collapsing the mechanism to one representative
         overlap
       - multi-insert previews now carry the whole ordered fragment chain
         through the same shared review/export path
       - cartoon panel text now avoids raw nucleotide strings inside the
         panels; the cartoon stays focused on mechanism while sequence details
         remain in the textual review blocks
       - the `Insert primer construction` block now shows both primers
         explicitly, including compact full-sequence plus `5' overlap` and
         `3' priming` parts
       - preview notes repeat the shared Tm-model assumptions so GUI and CLI
         users see the same explanation
    5. `Outputs`
       - shows the exact primer/product nodes the current apply step will
         create
       - `Apply Gibson Cloning` only enables after preview reports no blocking
         errors, then runs one shared engine operation and creates new sequence
         nodes for two primers per insert plus the assembled product
       - if preview is blocked, the status text now repeats the concrete
         blocking reason instead of only reporting an error count, including the
         multi-insert `defined opening` guardrail when applicable
       - in lineage `Graph` view, that apply step is projected as one explicit
         Gibson operation hub between the destination + insert inputs and all
         created outputs
       - in `Containers`, those three Gibson outputs land as separate
         singleton containers rather than one pooled multi-output container,
         matching the expected separate primer/product vials
       - in `Arrangements`, Gibson apply also records one reusable serial lane
         setup comprising the original vector, the ordered insert lane(s), and
         the assembled product, with recommended DNA ladders ready for gel
         export
       - `Cancel` closes the specialist without applying anything and keeps the
         current draft available when the window is reopened
       - the assembled product carries transferred destination/insert features
       - partially consumed destination annotations are trimmed or projected
         when a truthful product-side rewrite is available
       - surviving MCS annotations are revalidated against actual restriction
         sites on the assembled product, including new sites introduced by the
         insert
       - export plan JSON
       - export preview JSON
       - export primer summary
       - export cartoon SVG
       - hand off to Routine Assistant when current execution paths can accept it
  - uses the same shared preview command as CLI/shell:
    - `gibson preview PLAN_JSON_OR_@FILE [--output OUTPUT.json]`
  - shared mutating command parity:
    - `gibson apply PLAN_JSON_OR_@FILE`
  - GUI-first testing tutorial:
    - `docs/tutorial/03-05_gibson_specialist_testing_gui.md`
  - deliberately does not embed the generic PCR/qPCR specialist UI

<!-- manual-page-nav:start -->

---

[GUI Manual](../gui.md) | [Previous: Routine Assistant](routine-assistant.md) | [Next: Sequencing Confirmation Workspace](sequencing-confirmation-workspace.md)
<!-- manual-page-nav:end -->
