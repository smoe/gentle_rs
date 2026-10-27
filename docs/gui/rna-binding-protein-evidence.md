### RNA-Binding Protein Evidence

<!-- manual-counterparts:start -->
No direct CLI counterpart: This chapter describes desktop presentation or interaction; no dedicated command reference is paired with it.
<!-- manual-counterparts:end -->

- The Splicing Expert now also includes an `ATtRACT / RBP evidence` section:
  - engine-owned, splice-aware motif interpretation over the selected splicing
    group; the GUI is only a viewer/filter for the shared payload
  - on large loci this section now opens collapsed by default so the Splicing
    Expert window itself stays responsive during launch
  - settings exposed in-panel:
    - requested organism
    - minimum ATtRACT quality score
    - PWM mapping policy:
      - `PWM: strict same-length` keeps PWM ranking limited to rows whose motif
        length exactly matches the linked PWM width
      - `PWM: windowed submatrix` can additionally reuse one unique
        consensus-compatible PWM subwindow when the linked PWM is longer than
        the motif; retained hits show the chosen subwindow in the hit table
    - optional `Compare strict/windowed counts` toggle to show the alternate
      policy’s aggregate counts without switching the active hit table
    - donor/acceptor flank size
    - transcript-strand-only toggle
    - explicit species-fallback toggle
  - provenance strip reports:
    - active ATtRACT source (`runtime`, `session_override`, ...)
    - active motif count
    - active snapshot fingerprint (shortened in the panel, full value retained
      in the shared payload)
    - scanned transcript/window counts
    - total hit count and unique RBP count
    - resolved organism and species-match mode
  - warnings are shown explicitly when the engine had to fall back beyond an
    exact species match or when no transcript survived the active strand filter
  - `RBP summary` table groups by factor / organism / motif model and reports:
    - hit count
    - strongest match score
    - PWM provenance counts (`E` exact-length, `W` windowed, `C` consensus)
    - region-class distribution
    - supporting transcript ids
  - `Hit rows` table lists:
    - transcript
    - factor / matrix id
    - region class
    - genomic and local coordinates
    - motif and matched sequence
    - match score plus motif quality
    - provenance column:
      - mapping policy used
      - linked PWM status
      - optional PWM subwindow range when the windowed mode was used
  - the header summary now also spells out:
    - how many active ATtRACT rows are PWM-backed vs consensus-only
    - how many retained hits were exact-length PWM, windowed PWM, or
      consensus-only
    - optional alternate-policy count summary with hit and unique-RBP deltas
  - current v1 scan is deliberately conservative:
    - normalized consensus/IUPAC motifs from `ATtRACT_db.txt` remain the
      candidate gate
    - `strict same-length` is the default because ATtRACT documents the
      `Matrix_id -> pwm.txt` linkage but does not document a gapped or
      subwindow-mapping rule
    - the optional `windowed submatrix` mode is a GENtle heuristic and is
      always labeled as such in per-hit provenance
    - rows without mapped PWM data keep exact consensus/IUPAC matching only
    - no dense lane overlay or exon-exon junction scan yet
  - practical review workflow:
    - start with `PWM: strict same-length`
    - enable `Compare strict/windowed counts`
    - if the alternate summary adds only a few hits or RBPs, the heuristic is
      likely not changing the biological picture much
    - if the alternate summary adds many hits or RBPs, switch to
      `PWM: windowed submatrix` and inspect the per-hit provenance column
      carefully before trusting the extra evidence
    - treat `windowed submatrix` rows as exploratory support unless they remain
      convincing alongside the rest of the splicing evidence

<!-- manual-page-nav:start -->

---

[GUI Manual](../gui.md) | [Previous: RNA Read Panel](rna-read-panel.md) | [Next: RNA-Read Mapping Workspace](rna-read-mapping-workspace.md)
<!-- manual-page-nav:end -->
