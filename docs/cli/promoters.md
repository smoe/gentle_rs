# Promoter Design

<!-- manual-counterparts:start -->
Related GUI reference: [Loading and creating sequences](../gui/loading-and-creating-sequences.md) | [TSS / Regulatory DNA Display](../gui/tss--regulatory-dna-display.md).
These topic links do not imply identical feature coverage across interfaces.
<!-- manual-counterparts:end -->

- Shared shell:
  - Commands:
    - `promoters tss-list`: sorted registry metadata (`gentle.tss_collection_list.v1`), including legacy/invalid entries. Members are `not_checked`; use `promoters tss-collection ID` to validate explicitly.

- Shared shell:
  - Commands:
    - `promoters tss-view-svg SEQ_ID OUTPUT.svg [--report REPORT_JSON] [--motif ACCESSION]... [--score-kind KIND] [--keep-negative] [--span START..END] [--width PX]`:
      headless export using the native annotated-TSS renderer, as used by the
      GUI's **Export View SVG** for the same attachments. `--motif` accepts
      repeated flags or a comma list and computes local curves; `--report`
      attaches a validated profile report. `--span` is local 1-based inclusive.
      With `--motif`, only complete motif footprints inside that span are
      scanned on both strands; footprints crossing either boundary are excluded.
      For span length L and matrix length M the denominator is
      `2 * max(0, L - M + 1)` strand-windows; unavailable windows remain in that
      denominator but not the evaluated numerator. The trailing grey band marks
      starts without a complete footprint within the scored span.
      Admission still uses the **full annotated window** (50,000 bp and
      1,000,000 matrix/base combinations), not the narrow scan span. Full input
      validation and per-matrix background calibration still occur. Local scales
      and empirical quantiles (`llr_quantile`, `true_log_odds_quantile`) describe
      the scored span. This is not equivalent to cropping full-window GUI scores:
      GUI export and `--report` only clip already-scored data, retaining their
      original scoring bounds, counts and scales.
      Writes atomically and reports the lane count, exported span and SVG
      SHA-256. Requires GENtle annotated-TSS metadata on `SEQ_ID`: a plain locus
      is refused, not guessed. `--score-kind` is not inherited from the viewer —
      repeat it per export or it falls back to `llr_bits`.
      The typed result includes `tss_view_svg_export`
      (`gentle.tss_view_svg_export.v1`): exact operation/request hash, view and
      sequence hashes, geometry, span and SVG hash. This is artifact integrity,
      not scientific approval or a full environment audit.
      Every curve point and validity band is retained, but SVG hover titles are
      capped at 128 per lane/local strand, retaining strongest and sampled valid
      starts. Each lane discloses emitted/total/omitted titles. Unsampled starts
      lose their position-specific tooltip, not their scores. Top-three raw-score
      starts are summarized within each matrix lane in the displayed span, with
      TSS offsets and local strands; clipping does not change rank.
    - `promoters tss-view-svg --collection COLLECTION_ID OUTPUT_DIR` accepts the
      same optional export flags. Reuses `GetTssCollection` validation of every
      member and stored membership; no implicit selection, regeneration or
      approval. At most 32 members and 32 MiB for the complete bundle.
      The same local `--span` must fit every member. A supplied `--report` must
      match each member and is read once; no database query is performed.
      Stage ordered `001.svg`, `002.svg`, ... pages plus `index.html` and
      `receipt.json`, then publish to a **new directory** with an existing
      ordinary parent. Failures/cancellation publish no partial bundle; existing
      destinations are refused. The receipt binds the validated collection and
      its hash, exact operation, ordered pages and index bytes. Scale policy
      `per_window_per_lane` preserves each window's supplied lane scales,
      including a report's declared scaling; it does not recalibrate across
      windows or rank different matrices. Artifact-write confirmation still
      applies. ClawBio descriptor support is a separate follow-up.
    - `promoters tss-inventory REQUEST_JSON_OR_@FILE` (read-only exact-start preview)
    - `promoters tss-materialize REQUEST_JSON_OR_@FILE` (approval-bound selected windows)
    - `promoters tss-collection COLLECTION_ID` (persisted collection and member validation; [walkthrough](../tss_workspace.md))
    - `promoters tss-forget COLLECTION_ID` (metadata only; retains member sequences and lineage)

## Ortholog Promoter Cohorts and Comparisons

- `orthologs resolve-promoter-cohort --anchor-species SPECIES --anchor-genome GENOME_ID --anchor-gene QUERY --target-species SPECIES [--target-species SPECIES ...] [--target-genome SPECIES=GENOME_ID] [--transcript SPECIES=TRANSCRIPT_ID] --orthologs ORTHOLOG_RESOURCE.json [--relationship manual|co-regulated|anti-co-regulated] [--upstream-bp N] [--downstream-bp N] [--ambiguity-policy reject|first|preserve] [--catalog GENOMES.json] [--cache-dir PATH] [--path OUTPUT.json]`
  - Runs engine `ResolveOrthologPromoterCohort`.
  - Uses a local `gentle.ortholog_resource.v1` mapping table only; no live
    Ensembl or orthology API call is made.
  - Resolves the anchor promoter and target ortholog promoters with the shared
    prepared-genome promoter/TSS resolver. Species aliases in the resource are
    honored for matching.
  - Ambiguous target mappings are unresolved by default. `--ambiguity-policy
    first` chooses the stable first candidate and records a warning.
    `--ambiguity-policy preserve` keeps the target unresolved and emits
    ordered, structured `candidate_mappings[]` for review; it never promotes
    all candidates into the resolved promoter cohort.
  - Optional `--relationship` records an expected cross-species association
    without deriving evidence flags until a comparison is run.
  - Returns portable schema `gentle.ortholog_promoter_cohort.v1`.
- `orthologs promoter-comparison --cohort COHORT.json --motif TOKEN [--motif TOKEN ...|--motifs CSV] [--score-kind KIND] [--allow-negative] [--relationship manual|co-regulated|anti-co-regulated] [--expression-json JSON] [--expression-source-label LABEL] [--cutrun-dataset-id ID] [--cutrun-read-report-id ID] [--cutrun-normalization-json JSON_OR_@FILE] [--path OUTPUT.json]`
  - Runs engine `SummarizeOrthologPromoterComparison`.
  - Compares a resolved ortholog promoter cohort across separate evidence
    channels: promoter-sequence identity, TFBS score-track similarity, motif
    peak presence, optional expression rows, and CUT&RUN/occupancy
    comparability states.
  - CUT&RUN source ids can assign species/genome-matched occupancy support
    (`confirmed`, `nearby`, `motif_only`, `occupancy_only`, or `no_data`).
    Rows without matching provenance are `not_comparable`.
  - `--cutrun-normalization-json` accepts an inline JSON object, a JSON file
    path, or `@FILE`. The record supplies `normalization_method`, `unit`,
    `comparison_reference`, `provenance`, and `values[]`; every value names
    its species, normalized value, and contributing selected dataset/read
    report ids. Optional `gene_label`, `transcript_id`, and per-value
    `provenance` disambiguate assignments.
  - Quantitative output is fail-closed. Missing/incomplete metadata, unmatched
    promoter values, unselected sources, or species/genome-incompatible
    sources produce `cutrun_quantitative_comparison.status=not_comparable`
    and no numeric pairwise rows. Qualitative motif/occupancy states remain
    available, and raw peak/signal/read-count values are never compared
    directly across species.
  - Optional `--relationship` emits non-blocking expectation flags:
    co-regulated flags unexpected TFBS/CUT&RUN divergence, while
    anti-co-regulated flags unexpected concordance.
  - Returns portable schema `gentle.ortholog_promoter_comparison.v1`.
- `orthologs conservation-comparison --cohort COHORT.json [--min-conserved-bp N] [--path OUTPUT.json]`
  - Uses the resolved cohort's transcription-oriented, equal upstream/downstream
    windows; it does not repeat genome or ortholog lookup.
  - Globally aligns each target to the anchor with the shared pairwise aligner,
    so target insertions/deletions do not shift downstream comparisons.
  - Reports pairwise CIGAR/identity/coverage rows and exact intervals shared by
    every resolved target in anchor-relative and anchor-genomic coordinates.
  - Returns `gentle.ortholog_promoter_conservation.v1`. Conservation is a
    prioritization signal for regulatory inspection, not proof of function.

<!-- manual-page-nav:start -->

---

[CLI Manual](../cli.md) | [Previous: Sequencing Confirmation](sequencing-confirmation.md) | [Next: Promoter Cofactors](promoter-cofactors.md)
<!-- manual-page-nav:end -->
