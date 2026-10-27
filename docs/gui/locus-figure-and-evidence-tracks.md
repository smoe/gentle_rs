### Locus Figure And Evidence Tracks

<!-- manual-counterparts:start -->
Related CLI reference: [Figure and Layout Exports](../cli/rendering.md) | [Repeats and Regulatory Annotations](../cli/repeats-regulation.md).
These topic links do not imply identical feature coverage across interfaces.
<!-- manual-counterparts:end -->

- `Locus figure` is the graphical composition surface for the same shared
  `FeatureExpertTarget::GeneLocusEvidence` operation:
  - reuses the panel, expression, RNA/cDNA, probe, occupancy, and qPCR inputs
    entered under `Evidence`, then adds strand-aware flanks, probe-effect
    tables/contrasts, an explicit coordinate system, an occupancy-layout JSON,
    legacy motif-score settings or a provider-neutral regulatory-score-track
    JSON, a hidden/automatic/fixed genomic scale bar, and an output path
  - a readiness table distinguishes genome anchors, local files, imported
    panels, and projected occupancy-track names. Missing file-backed resources
    can be browsed to their relocated copy; the engine still performs schema,
    coordinate, and provenance validation during composition
  - `Compose / refresh` calls the shared pure-read feature-expert operation and
    renders its SVG through the shared renderer into the live GUI preview;
    transcript metrics, PSR/JUC effects, occupancy/motif summaries, warnings,
    and provenance remain inspectable beside the figure
  - the interactive inspector below the publication preview keeps separate
    transcript exon/CDS models, translation boundaries, occupancy/chromatin
    lanes (including unavailable states), retained motif/site calls, Ensembl
    Regulation, predicted TF-score, saved-region, and nucleotide geometry on
    one local axis. Hover a lane for its full label, source, units, and state.
    It supports fit/zoom, per-lane visibility, nucleotide letters when the
    current viewport gives each base enough horizontal space, and opening the same span in the full DNA
    map. Visibility changes affect inspection only; changing scientific inputs
    still requires explicit recomposition
  - dragging across the inspector stages an exact local interval as a
    `reporter_candidate` genomic ROI. Saving it uses the existing
    assembly-bound region manager, where its colour is persistent and remains
    presentation metadata rather than biological evidence. A portable locus
    report JSON can be loaded to inspect the exact evidence behind an existing
    document. Newly composed reports bind the exact DNA digest, length, and
    genome anchor; a mismatch is rejected even if the sequence id matches.
    Legacy reports without this binding remain historical evidence: live DNA
    letters, selection, staging, and DNA-map navigation are disabled until
    recomposition. Staged reporter intervals are checked again when saved
  - importing a reporter-architecture comparison retains its complete envelope,
    proposed architecture rows and material legend. JSON/copy, SVG, preview,
    and PDF preserve that imported comparison, rather than silently exporting
    only its nested locus report. The interactive rows reuse the canonical
    reporter geometry; LUC attachment markers are schematic, not genomic DNA
  - regulatory-score JSON may define several JASPAR matrices/TFs and offline
    external-model outputs with independent score semantics and scales. The GUI
    does not numerically compare unlike providers and never labels an
    uncalibrated prediction as affinity. Matrix/model-specific factor bindings
    remain distinct, and only an exact typed cross-source calibration id/digest
    can authorize a shared scale; descriptive prose cannot
  - every regulatory-score lane has a visible Y axis beginning at `0`, with
    midpoint and maximum ticks plus a score-family label. In particular,
    `llr_background_tail_log10` is labeled `-log10 Ptail`: it is a modeled
    background-tail probability derived from an LLR, not the raw log-odds
    score. Negative source samples remain auditable in report JSON but are
    floored at zero in the figure, and negative site/threshold labels are not
    drawn. Compressed SVG/PNG/PDF overview traces use unconnected pixel-bucket
    min/max whiskers and mean marks, separately for each strand, rather than
    joining maxima into a misleading plateau; the aggregation is labelled
    The interactive inspector uses the same bounds and score-family labels;
    negative and missing samples interrupt traces rather than connecting
    across unevaluated positions. It does not change stored scores or thresholds
  - H3K4me3 and related marks use the explicit `chromatin_context` lane role,
    visually separate from TF occupancy. Cell line, condition, batch, and
    scaling remain caller-declared rather than inferred from filenames
  - every explicitly requested occupancy/chromatin lane remains visible even
    when no interval can be drawn. Dashed rows report `not prepared`, `no
    compatible interval`, or `assembly mismatch` and retain stable source,
    assay, mark, and factor metadata, so unavailable evidence is not mistaken
    for an observed zero signal
  - report JSON, SVG, and PDF exports reuse shared operation/render paths. The
    GUI does not implement a second locus-composition algorithm
  - when the shared report carries source-bound Ensembl Regulation evidence,
    the preview adds a distinct Ensembl-annotated regulatory-regions section.
    Feature IDs and the provider label are live links in SVG; exported PDFs
    carry corresponding PDF URI annotations. The rows are annotation evidence,
    remain separate from occupancy/chromatin/score tracks, and do not enumerate
    SCREEN-derived project features
  - existing qPCR candidates can open their persisted report or create a
    deterministic reviewed oligo-order form. Junction evidence can seed the
    existing transcript-aware qPCR designer; this is assay planning, not a
    claim that the candidate validates the displayed relationship
  - the offline PATZ1 example is generated from
    `docs/examples/workflows/patz1_gene_locus_evidence_offline.json`; its tiny
    committed tracks are synthetic and are not experimental PATZ1 evidence

<!-- manual-page-nav:start -->

---

[GUI Manual](../gui.md) | [Previous: Splicing Structure And Evidence Views](splicing-structure-and-evidence-views.md) | [Next: Saved Genomic Regions](saved-genomic-regions.md)
<!-- manual-page-nav:end -->
