## TFBS Annotation (Engine Ops)

<!-- manual-counterparts:start -->
Related CLI reference: [TFBS and JASPAR Scans](../cli/tfbs.md).
These topic links do not imply identical feature coverage across interfaces.
<!-- manual-counterparts:end -->

The Engine Ops panel includes `TFBS annotation (log-likelihood ratio)`:

- motif selection:
  - `All known JASPAR motifs` mode (explicit all)
  - or selected motifs as comma-separated IDs/names/IUPAC
  - selection helper:
    - filter JASPAR entries
    - add entries to selected motif list with `+`
- global thresholds:
  - `min llr_bits` (absolute log-likelihood ratio score, base 2)
  - `min llr_quantile` (empirical quantile in the scanned sequence region; both strands)
- per-TF overrides:
  - `per-TF min llr_bits` as `TF=VALUE,...`
  - `per-TF min llr_quantile` as `TF=VALUE,...`
- `Clear previous TFBS`:
  - removes prior generated TFBS annotations before writing fresh results

Execution calls engine operation `AnnotateTfbs` and writes scored TFBS features
onto the active sequence. Generated TFBS qualifiers now include four scores:
`llr_bits`, `llr_quantile`, `true_log_odds_bits`, and
`true_log_odds_quantile`.

Safety behavior:

- GUI uses engine default TFBS cap (`500` accepted hits per operation) to keep
  the UI responsive on dense scans.
- CLI/JSON workflows can override this via `AnnotateTfbs.max_hits` (`0` means
  unlimited).
- The DNA-window toolbar `TFBS scan` menu reuses the same motif/threshold
  settings for non-mutating quick inspection of the current selection, visible
  span, or whole sequence, with a smaller GUI-only preview cap (`200`) to keep
  those ad-hoc scans responsive.
- The DNA-window toolbar `TFBS score tracks` menu reuses the same motif
  selection for non-mutating continuous score inspection of the current
  selection, visible span, or whole sequence.
- The same DNA window now also exposes a `Direct scan inspectors` section:
  - `RE scan` now populates a cached restriction-site table plus shared JSON
    export button
  - `TFBS scan` now populates a cached TFBS-hit table plus shared JSON export
    button
  - clicking a cached restriction-site or TFBS-hit row now selects that span in
    the active DNA window and recenters the linear view there when applicable
  - `TFBS score tracks` remains the cached score-track/SVG sibling in the same
    window
- The same `TFBS annotation` panel now also includes one cached score-track
  inspector:
  - choose the exported/displayed `value kind`
  - optionally `clip negatives` for raw bit-score families
  - inspect the shared `TfbsScoreTrackReport` immediately in the DNA window
  - each motif row now includes a `PSSM…` action that opens the shared
    `JASPAR Expert...` window focused on that motif, which is a better home
    for matrix/logo interpretation than the cramped score-track side rail
  - when the active sequence already carries imported BED-style track features
    (for example projected CUT&RUN peaks), the same inspector now shows those
    intervals as extra shared overlay rows underneath the motif score tracks
    instead of leaving them trapped in the generic feature map alone
  - export the cached view through the shared `RenderTfbsScoreTracksSvg`
    engine route instead of a GUI-only painter
- Threshold fields remain relevant for hit annotation and `TFBS scan`; the
  continuous score-track view depends on motif selection plus the chosen
  score-family/clipping mode.

The DNA-window Engine Ops panel also includes a `CUT&RUN regulatory support`
inspector for release smoke/proof runs. It is deliberately thin GUI plumbing:
dataset ids, saved ROI read-report ids, optional catalog/cache paths, an
optional promoter span, neighbor window, and species filters are sent through
the shared `InspectCutRunRegulatorySupport` engine operation. The panel
displays the returned `gentle.cutrun_regulatory_support.v1` record: evidence sources,
support windows, the `TFBS + occupancy support` table with additive
`support_status` and distance fields, motif-absent supported windows,
recurring motif context, warnings, and JSON export. It does not add GUI-only
CUT&RUN scoring or motif interpretation logic. Blank catalog/cache fields use
normal discovery; explicit values let the GUI reopen reports backed by a
project-specific prepared cache.

While TFBS annotation is running, GUI shows live progress indicators and keeps
repainting until completion:

- per-motif percentage (based on actual scan steps)
- total percentage across all selected motifs (fraction of TFs addressed)

When TFBS score tracks are being computed, the same GUI now also shows live
progress for:

- deterministic background calibration
- target-span scanning
- overall motif completion across the requested track set
- the progress panel now also exposes `Cancel ...` so long-running
  background-calibrated score-track jobs can be stopped cooperatively

TFBS display reduction (no recomputation needed):

- `TFBS display filter` offers four checkbox-enabled criteria:
  - `llr_bits`
  - `llr_quantile`
  - `true_log_odds_bits`
  - `true_log_odds_quantile`
- each enabled criterion applies its threshold live to visible TFBS features
- disabling all criteria shows all TFBS features
- sequence SVG export now uses the same TFBS filter settings as the GUI display

<!-- manual-page-nav:start -->

---

[GUI Manual](../gui.md) | [Previous: Array Tracks In DNA Viewer](array-tracks-in-dna-viewer.md) | [Next: TSS / Regulatory DNA Display](tss--regulatory-dna-display.md)
<!-- manual-page-nav:end -->
