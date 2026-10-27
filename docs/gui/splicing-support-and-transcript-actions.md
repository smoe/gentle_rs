### Splicing Support And Transcript Actions

<!-- manual-counterparts:start -->
Related CLI reference: [Source-Coherent Transcripts](../cli/transcripts.md).
These topic links do not imply identical feature coverage across interfaces.
<!-- manual-counterparts:end -->

- The splicing expert window uses its own window-styling slot (`splicing`) so
  tint/image backdrop can be configured separately from DNA and pool windows.
- The Agent Assistant window uses its own window-styling slot (`agent assistant`)
  with `assets/backgrounds/AI_AgentInterface.png` as the default image.
- Splicing support frequencies are shown explicitly:
  - hovering exon glyphs in the lane canvas shows transcript/exon coordinates,
    support, `len%3`, and CDS flank phase details (when available)
  - hovering splice-site boundary markers now shows donor/acceptor motif,
    paired intron signature, motif class, and a short interpretation note
  - splice-site motifs are now summarized explicitly:
    - a `Splice-site motifs` table lists donor/acceptor positions and dinucleotide
      motifs, the paired intron signature (`GT-AG`, `GC-AG`, `AT-AC`, ...), and
      the best current high-level annotation GENtle can provide
    - canonical major-spliceosome `GT-AG` sites stay green/teal, known
      non-canonical classes such as `GC-AG` and `AT-AC` use distinct colors,
      and other non-canonical signatures remain red
- acceptor-proximal intron heuristics are now surfaced conservatively:
  - intron lines can highlight a branchpoint-like site and a pyrimidine-rich
    tract near the acceptor when the simple heuristics find one
  - hovering an intron line reports the current branchpoint-like motif/score
    and the best detected polypyrimidine tract span; intron lines can now be
    selected directly, Shift/Command-click toggles multi-selection, and the
    lane context menu exposes intron selection actions
    - an `Acceptor-proximal intron signals` table summarizes those heuristics
      explicitly and labels them as heuristic evidence rather than a splice
      predictor
  - rows in that intron-signal table are clickable:
    selecting one keeps the corresponding intron highlighted in the lane view
    and shows a small description block with the current branchpoint-like and
    polypyrimidine-tract interpretation
  - the same selected introns now also feed one `Intron-centered regulatory
    interpretation` table:
    - it combines those conservative intron heuristics with any already cached
      ATtRACT/RBP evidence for the same splicing target
    - each row reports donor-flank, acceptor-flank, and intron-body hit
      counts plus a short top-RBP summary
    - clicking those rows keeps the same intron focused in the lane view, so
      splice-mechanistic and regulatory context stay aligned
- exon columns include support as `n/N (%)` (plus `const` for constitutive
  exons)
  - exon lane glyphs can show CDS flank phase coloring on the left/right exon
    edges (`0/1/2`) when transcript `cds_ranges_1based` are available; the
    transcript-entry phase (`left` on forward transcripts, `right` on reverse
    transcripts) says whether the exon starts at a codon boundary/new amino
    acid (`0`) or enters a split codon (`1`/`2`)
  - exon columns also expose `len%3` color cues (`0/1/2`) based on genomic
    exon length modulo 3 (heuristic frame cue)
  - ordinary linear sequence-map exon blocks for `mRNA` / `transcript` /
    `exon` / `CDS` features keep their normal feature color but gain a thin
    top stripe using the same modulo palette: blue=`0`, amber=`1`, rose=`2`
    - hovering a striped exon block reports the specific exon block's `len%3`
      value and the same frame-neutral / frame-shifting skip hint used by the
      selected-feature details pane
  - the exon-skip planning checkbox tooltips report each candidate's `len%3`
    value, coding-only skip length modulo 3, flanking intron lengths,
    transcript position, exon support frequency, and CDS entry phase; quick-add
    buttons can select candidates by `len%3=0/1/2`, `coding%3=0/1/2`,
    coding context, codon-boundary starts, or split-codon starts
  - saved exon-skip plans carry those same frame/phase attributes in candidate
    rows, so GUI, CLI, MCP, and future adapters can inspect the same
    engine-owned facts
  - transcript-vs-exon matrix cells are color-coded by exon support frequency
    (higher support => stronger color intensity)
  - an `Exon -> exon transition matrix` shows predicted transition support
    counts with frequency-coded cells; row/column exon headers reuse `len%3`
    colors
  - a `Junction transition support` table lists donor/acceptor transitions with
    support counts and percentages
  - arc-width encoding remains active and arcs also show support-count labels
    in the lane canvas
- The splicing expert window includes a quick action:
  - `Send Group ROI -> Primer/qPCR`
  - seeds ROI start/end fields in Engine Ops primer and qPCR design forms from
    the current splicing-group genomic bounds.
- The splicing expert window also includes transcript-level quick actions:
  - `Transcript` selector (choose one lane by `n-<feature_id> <transcript_id>`)
  - `Derive group transcripts`
    - derives cDNA transcript sequences for all transcript lanes admitted by
      the current splicing scope.
  - `Derive all mRNA`
    - derives cDNA transcript sequences for every `mRNA`/`transcript` feature
      on the active sequence.
    - keeps the DNA-window status compact for large transcript sets instead of
      dumping one line per derived transcript; the derived set still becomes the
      current pool/export target in Engine Ops.
  - `Derive + Dotplot`
    - derives the selected transcript sequence, switches to transcript context,
      and opens pairwise transcript-vs-genomic dotplot.
    - strand-aware default mode:
      - transcript `+` strand -> `pair_forward`
      - transcript `-` strand -> `pair_reverse_complement`
    - step size is auto-increased when needed so pair evaluations stay within
      engine limits.
  - `Window guide [?]` hover tooltip summarizes the whole Splicing Expert:
    annotation structure, transcript quick actions, and report-driven
    RNA-read evidence for the selected locus.
- The Splicing Expert now also includes `Exon skip planning`:
  - check exon candidates from the selected transcript lane
  - add candidates that overlap the active map/text selection
  - optionally select by overlapping feature kind/label, which can cover
    imported SNP/variation, antisense RNA, repeat-masker, or other annotation
    tracks represented as sequence features
  - `Build / Refresh selection plan` stores a replayable
    `gentle.exon_skip_selection_plan.v1` in engine state
  - `Materialize skipped isoform` consumes the stored plan and creates both a
    genomic annotation product and retained-exon cDNA/mRNA sequence
- The Splicing Expert now includes a compact `RNA-read evidence` section:
  - `Report` selector filtered to the current `seq_id + seed_feature_id`
  - newest matching report is auto-selected when no explicit report is chosen
  - empty state:
    `No RNA-read report for this splicing group yet` plus
    `Open RNA-read Mapping Workspace...`
  - saved reports drive score density, thresholded cDNA support, mapped cDNA
    support, read-effects inspection, and report-driven review
  - `Isoform read support` is an IsoVis-inspired categorical overlay for saved
    RNA-read reports: each mapped isoform row shows a conservative triage chip,
    stacked read-bin bar, aligned/MSA counts, mean identity/coverage,
    secondary-mapping count, and `Audit` / `Export...` controls for the exact
    contributing reads. The bins are inspection labels, not expression values
    or novel-isoform calls. Full IsoVis-style GTF/GFF/BED stack import,
    drag-reorder, zoom-to-coordinate controls, PNG/PDF export, and live
    Ensembl/InterPro GUI querying remain outside this V1 unless covered by
    other GENtle tools.
  - the dedicated `RNA-read Mapping` workspace now stays open if you close the
    parent DNA sequence window that originally spawned it
- The Splicing Expert now includes a distinct `Array probe geometry` section:
  - consumes the cached
    `gentle.probe_region_evidence_interpretation.v2` report for the same
    splicing locus
  - lists exon-overlap PM probes, exon-exon junction-spanning probes, parent
    probeset grouping, and ambiguity/unresolved mapping tags
  - promotes junction-spanning probe rows into a compact
    `Junction-targeting array probes` preview before the full geometry table
  - this is a review-only array design/alignment constraint layer, separate
    from RNA-read evidence; it does not infer isoform support by itself

<!-- manual-page-nav:start -->

---

[GUI Manual](../gui.md) | [Previous: Splicing Expert Workspace](splicing-expert-workspace.md) | [Next: Splicing Structure And Evidence Views](splicing-structure-and-evidence-views.md)
<!-- manual-page-nav:end -->
