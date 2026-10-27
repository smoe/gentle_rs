## Loading and creating sequences

<!-- manual-counterparts:start -->
Related CLI reference: [Project and Sequence State](../cli/projects.md) | [Protein and Structure Evidence](../cli/protein-evidence.md) | [Promoter Design](../cli/promoters.md).
These topic links do not imply identical feature coverage across interfaces.
<!-- manual-counterparts:end -->

Use the top application menu:

- `File -> New Sequence...`
  - opens a small editor for typed or pasted IUPAC DNA
  - whitespace and digits are ignored; `U` is normalized to `T`
  - optional `id`, `name`, and linear/circular topology fields are applied by
    the shared `CreateSequenceFromText` engine operation
  - creating the sequence imports it into the current project and opens a
    normal DNA sequence window
- `File -> New Sequence from Clipboard...`
  - reads the system text clipboard into the same `New Sequence` dialog
  - the clipboard text is shown for review first; the project is not mutated
    until `Create Sequence` is clicked
  - empty, non-text, or invalid clipboard contents leave the dialog open with a
    status message
- `File -> Open Sequence...`
  - the file picker supports selecting multiple sequence files at once
  - selected files are imported sequentially through the normal per-file import path
  - each successfully imported file opens its own sequence window
- `File -> Protein Evidence...`
- `File -> Open Project...`
- `File -> Open Recent Project...`
- `File -> Open Tutorial Project...`
- `File -> Save Project...`
- `File -> Close Project`
- `File -> Import REBASE Data...`
- `File -> Import JASPAR Data...`
- `File -> JASPAR Expert...`
- `File -> Export DALG SVG...`
- `File -> Quit`

Save/export filename defaults:

- `Save Project...` suggests the current project filename when the project
  already has a path.
- For unsaved projects, `Save Project...` derives the suggested filename from
  the current project name instead of using a fixed generic placeholder.
- `Export DALG SVG...` suggests the project name stem plus `.lineage.svg`.

Window/application close shortcuts:

- `Cmd+W` closes the active auxiliary/native GENtle window on macOS
- `Ctrl+W` is also accepted as a close-window shortcut
- `Cmd+Q` quits GENtle on macOS
- `Ctrl+Q` is also accepted as a quit shortcut

Genome operations are available from the dedicated `Genome` menu:

- `Genome -> Prepare Reference Genome...`
- `Genome -> Prepared References...`
- `Genome -> Retrieve Genomic Sequence...`
- `Genome -> BLAST Genome Sequence...`
- `Genome -> Import Genome Track...`
- `Genome -> Prepare Helper Genome...`
- `Genome -> Retrieve Helper Sequence...`
- `Genome -> BLAST Helper Sequence...`

`Close Project` availability:

- disabled when the app is on an untouched untitled project
  (no loaded project path and no user content yet)
- enabled once a project is loaded or any sequence/container content exists

Tutorial projects:

- `Open Tutorial Project...` builds a project from canonical tutorial workflow
  examples (`docs/examples/workflows`) in the background and opens the finished
  project for inspection when the build completes.
- The same menu also starts with a `Guided walkthroughs` submenu for
  documentation-only tutorials, such as the Agent Assistant and Agent
  Interfaces tutorial; those entries are grouped by tutorial content area and
  open the Help tutorial page directly instead of building project state.
- When a tutorial chapter declares a matching guide, the Help window opens that
  tutorial page automatically so the next GUI steps are visible right away.
- tutorial project discovery now consults `docs/tutorial/catalog.json` first and
  resolves the executable manifest/runtime chapter list from the catalog's
  `generated_runtime.manifest_path`
- Chapters are grouped by tutorial content area using the derived decimal ids
  from `docs/tutorial/catalog.json` / `docs/tutorial/manifest.json`; tier and
  online status remain visible in hover text and labels.
- Tutorial project and guided-walkthrough menus use the same window-bounded
  scrolling as `Help -> Tutorials`, so later entries remain reachable in a
  small window.
- Generated tutorial project files are written under a process- and
  invocation-specific directory below the system temp directory
  (`.../gentle_tutorial_projects/process-PID-invocation-N`) and opened without
  adding those temp files to the `Open Recent Project...` list. Separate builds
  therefore never delete or overwrite one another's workspace.
- While a tutorial project is being built, GENtle keeps the GUI responsive and
  shows live status in the `Background Jobs` panel.
- Online tutorials can therefore be cancelled explicitly:
  - from the `Open Tutorial Project...` menu while the build is running
  - or from the `Background Jobs` panel
- Progress is best-effort:
  - GENtle always reports coarse setup/build/persist phases
  - shared engine operations inside the tutorial workflow contribute their own
    progress when they support it
  - some remote/network waits still only expose the coarse phase label
### Promoter Design Window

- Additional GUI-first manual tutorial:
  - `docs/tutorial/08-04_vkorc1_warfarin_promoter_luciferase_gui.md`
  - includes a stepwise GUI workflow for `VKORC1` / `rs9923231`
    promoter-luciferase planning using the dbSNP fetch path, derived
    promoter-window / promoter-context reasoning, and one pinned local
    mammalian luciferase backbone, with per-step GUI/CLI parity mapping.
  - the sequence window now exposes a dedicated `Promoter design` window for
    `variation`, `gene`, `mRNA`, and `promoter` features via the feature
    description pane, feature-tree context menu, and map context menu; it
    drives shared `AnnotatePromoterWindows`,
    `SummarizeVariantPromoterContext`,
    `SuggestPromoterReporterFragments`, `SummarizeTfbsScoreTracks`,
    `ExtractRegion`, `MaterializeVariantAllele`, `LoadFile`, `Ligation`, and
    `Branch` operations without dropping into generic `Engine Ops`.
  - `Make reference/alternate inserts` validates both choices before creating
    either sequence. An empty multiallelic choice or invalid alternate leaves
    project sequences and output-name fields unchanged, allowing a corrected
    retry without orphan reference copies. After preflight, both shared
    operations run under one engine write lock and the existing rollback
    transaction. A later execution refusal restores the prior project,
    journal and undo/redo history before any success-driven view change;
    output names update only after both inserts succeed. Successful creation
    retains two ordinary undo entries. This is GUI composition of existing
    operations, not a new biological operation or a wet-lab validation claim.
  - analysis and export actions wrap to the visible viewport even when the
    parameter grid widens the scrollable content. Vertical scrolling reaches
    later action rows; the parameters retain their existing horizontal scroll.
    The pair action remains at the start of its separate row. Automated egui
    geometry is not proof of native label readability or screenshot acceptance.
  - the same window now also exposes positive-only TF score tracks across the
    chosen promoter span, with TERT-friendly `SP1`, Yamanaka-plus-`NANOG`, and
    `TP73/SP1/BACH2/PATZ1` presets meant for promoter interpretation before
    reporter construction.
  - the same score-track panel now also lets users switch the exported/painted
    value kind between raw bits (`llr_bits`, `true_log_odds_bits`), in-window
    quantiles (`llr_quantile`, `true_log_odds_quantile`), random-background
    percentiles (`llr_background_quantile`,
    `true_log_odds_background_quantile`), and random-background tail emphasis
    as `-log10(tail)` (`llr_background_tail_log10`,
    `true_log_odds_background_tail_log10`), so promoter inspection is not
    locked to one scoring view.
  - transcript-derived promoter windows now avoid visually stacking identical
    DNA-level promoter spans just because splice variants diverge downstream:
    exact-shared windows collapse to one promoter symbol and annotate the
    contributing transcript count in the label (`(n tx)`).
  - the map context menu also now recognizes reasoning-derived promoter spans
    and can resolve them back to the transcript feature that seeded the window,
    so `Open Promoter Design` works even before the promoter window has been
    written back as a standalone feature.
  - the same score-track panel now also shows deterministic random-background
    normalization context per motif so permissive score families can be read
    against their baseline instead of as free-floating texture:
    - raw bit views show `p99`, `Δp99`, and `bg+`
    - background-normalized views show `theory max`, `peak q`, and
      `-log10 tail`
  - each motif row now also exposes `PSSM…`, which opens the shared
    `JASPAR Expert...` window on that motif instead of squeezing the full
    matrix/logo interpretation into the score-track side rail
  - the background-normalized percentile and tail views intentionally suppress
    everything below the `0.95` deterministic random-background quantile, which
    keeps the plot focused on unusual windows rather than random low-level
    motif texture.
  - the background-normalized views are now calibrated against a quantized IID
    random-DNA window model instead of only against a finite sampled
    background, which keeps strong motifs from collapsing onto the same visual
    ceiling quite so easily.
  - the same panel now also lists a small `top peaks` summary per motif so
    users can jump from “busy trace” to the concrete strongest local windows
    without manually scanning the entire line plot first.
  - the same panel now also presents `TFBS track correlation`, which keeps
    one selectable correlation family over the displayed per-position signals:
    - `Pearson`
    - `Spearman`
    and one selectable strand-handling view:
    - `max_strands`
    - `forward_only`
    - `reverse_only`
    and shows each family as:
    - smoothed correlation is the main “do these maxima live in the same
      neighborhood?” cue
    - raw correlation stays visible as the unsmoothed sanity check
    - a per-motif strand-pairing table above the pairwise grid reports how
      each motif’s own forward and reverse curves travel together, including a
      signed reverse-minus-forward primary-peak offset
  - when score tracks are being recomputed, the window now keeps a live TFBS
    progress widget visible so deterministic background calibration and
    target-span scanning do not look like a stall on longer/randomized runs
    - a signed peak offset reports whether one motif’s strongest window leads
      or trails another
  - the same window now also exposes `TFBS similarity`, which brings the
    newer anchor-vs-candidate ranking route directly into Promoter design:
    - one anchor motif is compared against either all JASPAR motifs or one
      restricted candidate list over the exact same promoter span used for
      score tracks
    - ranking metric, optional species filters, cached remote-metadata
      enrichment, and result limit are controlled in the window
    - the ranked table is cached in-window and can be exported as JSON through
      the same shared engine report instead of becoming a GUI-only table
  - the same window now also exposes `Compare alternative promoters`, which
    groups transcript-derived promoter windows by exact DNA span:
    - rows show transcript multiplicity (`n tx`) instead of visually stacking
      identical promoter windows that only differ by downstream splice usage
    - `Use this row` retargets the Promoter design score-track / similarity
      span to that grouped promoter interval and representative transcript
  - the same window now also exposes `Compare isoform evidence`, which calls
    the shared `SummarizeIsoformPromoterComparison` route:
    - rows keep common evidence separate from promoter-group-specific evidence
      so TP73-style isoform promoter differences are not flattened into one
      undifferentiated ledger
    - `Use this row` retargets the Promoter design span to the selected
      isoform promoter group while leaving the comparison table available for
      orientation
    - `Export isoform comparison JSON...` writes the cached
      `gentle.isoform_promoter_comparison.v1` payload through the same engine
      operation used by shell/CLI workflows
  - the same window now also exposes `Build evidence matrix`, which calls the
    shared `SummarizePromoterEvidenceMatrix` route and keeps the promoter
    evidence ledger visible in-window:
    - rows group DNA-level promoter candidates with transcript support plus
      existing local TFBS, variant, repeat, enhancer/terminator,
      external-interval, and CUT&RUN-tagged feature evidence
    - `Use this row` retargets the Promoter design span to the selected
      candidate for score tracks, TFBS similarity, and reporter planning
    - `Export evidence matrix JSON...` writes the cached
      `gentle.promoter_evidence_matrix.v1` payload through the same engine
      operation instead of serializing a GUI-only table
  - the same window now also exposes pasted or loaded promoter expression rows
    through the shared `SummarizePromoterExpressionEvidence` route:
    - the table reports assigned/unassigned expression rows per promoter group,
      mean/max values, units, conditions/samples, and match evidence
    - text in the panel deliberately frames expression as association evidence,
      not promoter causation or wet-lab validation
    - warnings and unassigned expression records remain visible
    - `Export expression JSON...` reruns the same shared operation with an
      output path, keeping GUI and shell artifacts identical in shape
  - the same window now also exposes local ortholog promoter cohorts through
    the shared `ResolveOrthologPromoterCohort` and
    `SummarizeOrthologPromoterComparison` operations:
    - the ambiguity selector exposes `reject`, deterministic warning-bearing
      `first`, and `preserve`; preserved mappings remain unresolved and are
      displayed as reviewable candidate rows with identity, relationship,
      confidence, source, evidence, genome, and context
    - rows show per-species promoter windows, strand/TSS placement, pairwise
      TFBS similarity, relationship expectation flags, and CUT&RUN support
      states
    - optional normalized CUT&RUN JSON names the method, unit, shared
      comparison reference, provenance, and one source-bound value per
      resolved promoter; incomplete records remain `not_comparable` and do
      not cause raw cross-species intensity comparison
    - `Export ortholog cohort JSON...` and `Export ortholog comparison
      JSON...` rerun the same shared operations with an output path instead of
      serializing GUI-only state
    - wording frames the result as cross-species association evidence, not
      proof of conserved regulation
  - the same window now also exposes `Inspect TFBS occupancy support`, which
    calls `InspectCutRunRegulatorySupport` for the active promoter span:
    - dataset ids and saved read-report ids are passed through the engine
      operation unchanged
    - the reused support table shows motif, motif score, nearest/overlapping
      CUT&RUN support, `support_status`, and distance without inventing a GUI
      verdict
    - the engine report still preserves the legacy confirmed/unconfirmed row
      vectors while displaying additive four-state `support_status`
      (`confirmed`, `nearby`, `absent`, `motif-poor`)
  - the same window now also exposes `Reporter architecture comparison
    (read-only)` as a thin client of
    `ComparePromoterReporterArchitectures`:
    - `Prepare from current locus` seeds one editable portable request from the
      selected sequence, gene/transcript context, promoter extents, and current
      motif list; when Splicing Expert > Locus figure has a configured isoform
      panel, it also attaches that shared locus-evidence request. The status
      explains the missing configuration when only the legacy comparison can
      be prepared; invalid existing Locus figure paths/settings fail visibly
      instead of silently dropping the evidence. `Compare / refresh` runs the
      shared engine operation
    - summary, TSS classes, staged panel, warnings, non-claims, CUT&RUN
      readiness, and the integrated occupancy/score/scale summary come directly
      from the cached report
    - the collapsed `Architecture models` list exposes stable architecture id,
      transcript, modeled insert length, spliced-UTR and intron burden, explicit
      ATG policy, and a `Use transcript` navigation action without recalculating
      biology in the GUI
    - JSON and SVG paths invoke the shared exporter/renderer; leaving either
      path empty keeps that artifact in memory
    - every reporter SVG uses visible material swatches for genomic source,
      spliced exon/cDNA, and schematic `LUC` sequence; combined locus figures
      retain the same keys, and their PNG/PDF exports derive from that SVG
    - an existing laboratory construct remains unverified until exact insert,
      vector, and boundary evidence is supplied, while CUT&RUN overlap remains
      occupancy support rather than evidence of TSS use or reporter output
  - the same window now also exposes `Promoter-reporter panel` as a thin client
    of the shared `PlanPromoterReporterPanel` and
    `MaterializePromoterReporterPanel` operations:
    - request JSON selects already-reviewed promoter candidates and one loaded,
      exactly catalog-validated vector; the GUI does not resolve or rerank them
    - request JSON must name an effective mutation policy:
      `native_only_v1` creates one unchanged regulatory construct, while
      `p53_family_core_disruption_v1` creates the established
      wild-type/mutant pair and requires a motif-hit anchor; members can
      override the request default, and `unspecified` remains invalid
    - member `evidence[]` rows bind independent response or occupancy sources;
      the panel summarizes their count and the proposal JSON retains their full
      provenance without merging them semantically with motif annotations
    - study-specific interpretation limits belong in request-bound
      `scientific_caveats`; GENtle preserves them in the proposal instead of
      injecting gene names or study claims from the generic operation
    - `Plan panel` is read-only and displays member boundaries/roles, effective
      policy, anchor/evidence provenance, optional stated-rule mutation audit,
      shared restriction or Gibson strategy, primer readiness, predicted
      circular products, preserved catalog-required vector annotations, exact
      artifact paths, warnings, and biological non-claims
    - planning clears any previous approval text, and materialization remains
      disabled until the user types the complete displayed proposal digest
    - approved materialization still recomputes the proposal through the engine
      and displays the resulting sequence IDs, artifact paths, and manifest;
      drift or an existing output path is reported rather than bypassed
  - `Regulatory-fragment panel` in the same window is a thin client
    of `PlanRegulatoryFragmentPanel` and
    `RenderRegulatoryFragmentPanelSvg`:
    - selectors bind candidate A, optional partner B, required minimal promoter,
      and optional reference control from the persisted genomic-region store;
      the GUI carries each exact region-set content digest into the request
    - reference release, loaded vector/catalog identity, helper catalog, member
      and insert-length bounds remain explicit inputs
    - partner, order, orientation, and spacing comparisons are opt-in; spacing
      carries exact spacer DNA and no geometry is inferred by the GUI
    - the cached report shows ordered source instances, genomic intervals,
      pairwise contrasts, inclusion and omission reasons, and every evidence
      lane with its independent state, observations, warnings, and blockers
    - JSON export reruns the same read-only planning operation; SVG export sends
      the unchanged cached plan through the shared digest-validating renderer
    - typing the exact panel digest records review locally but does not authorize
      product creation. In `Exact design products`, enter an explicit output-ID
      prefix and choose `Prepare exact products`; this read-only action uses
      `PlanRegulatoryFragmentMaterialization`, also reachable through
      `promoters regulatory-products-plan`
    - each proposed product exposes its full DNA, topology, sequence hash,
      annotations, source instances, and omitted vector-annotation indices;
      `Copy product proposal JSON` copies the unchanged portable proposal
    - `Create approved design products` requires the separate, complete product
      digest. The shared `MaterializeRegulatoryFragmentPanel` operation rechecks
      sources and collisions and creates all products atomically, with one undo.
      The receipt lists output IDs and keeps final-product audit `not_evaluated`:
      these are designed sequences, not verified cloning reactions
    - input edits, replanning (including JSON re-export), and failed preparation
      clear product review. A creation attempt consumes GUI approval; engine
      rejection clears the stale proposal rather than silently rebuilding it
    - source-bound external locus observations retain their original typed
      payload in JSON and independent availability in the view
    - labels remain hypotheses or unresolved comparisons; the panel does not
      claim sufficiency, enhancer/silencer activity, partner dependence, or
      causal regulation
    - before fixing the ROI, GUI Shell also accepts `promoters fragment-candidates`
      for [evidence-guided boundaries](../reporter_fragment_selection.md), including
      explained extension/shortening, provenance-bound called-peak BED sources
      and optional exact-vector MCS checks. Generic coverage is not a peak seed.
      This is a read-only proposal, not a new graphical drag editor; saving a
      proposed region and approving its construct remain separate actions.
      For fixed similarity windows, GUI Shell `op @geometry_op.json` uses the
      same stateless `ComputeTssWindowGeometry` as CLI/MCP, rather than selecting
      adaptive inserts. Neither route requires or introduces a new GUI editor.
    - the executable offline tutorial
      `08-11_regulatory_fragment_panel_planning_offline` walks through the
      region selectors, exact A/B geometry, independent evidence states, and
      digest-valid JSON/SVG review without materialization
  - the same window now also exposes `Export TF score tracks SVG...`, which
    goes through the shared `RenderTfbsScoreTracksSvg` engine route instead of a
    GUI-only painter and can therefore reproduce the same stacked figure style
    from shell/CLI workflows and README assets.
  - the same window now also exposes `Export TF correlation SVG...`, which
    goes through the shared `RenderTfbsScoreTrackCorrelationSvg` route and
    writes one dual heatmap view:
    - smoothed selected-metric correlation on the left
    - raw selected-metric correlation on the right
    - selectable `max_strands` / `forward_only` / `reverse_only` signal source
    - when `max_strands` is selected, the export becomes one strand-expanded
      all-vs-all matrix:
      - every motif appears twice on each axis, ordered `F` then `R`
      - each TF-pair therefore reads as one adjacent 2x2 block
        `F-F / F-R / R-F / R-R`
      - the block values come directly from the strand-specific curves rather
        than from a prior `max(forward, reverse)` collapse
    - ranked synchronized pairs with signed peak offsets in the footer
  - stacked SVG exports now render motif labels as `TF name (JASPAR id)` and
    carry one shared dashed TSS line with a single top-level kinked arrow when
    the promoter span has either feature-derived or provenance-derived TSS
    context available; the TSS label now also gets a small outline box so the
    shared marker survives README/downsampled figure use better.
  - the same window now also exposes `Export handoff bundle`, which writes a
    portable ClawBio-facing bundle containing:
    - promoter-context JSON
    - promoter-candidate JSON
    - promoter-context SVG
    - paired reporter SVG previews
    - `report.md`
    - `result.json`
    - `commands.sh`
  - the tutorial uses a TSS-centered reverse-strand `VKORC1` fragment that
    keeps `rs9923231` inside the cloned insert instead of on its edge.
  - `Genomic alleles` exposes the marker's supplied ref/alt qualifiers.
    `Alternate base` is prefilled only for a single alternative; multiallelic
    markers require an explicit reviewed A/C/G/T choice. The field is passed
    unchanged to the shared engine; the UI never chooses an allele by gene
    strand or a pharmacogenomic label. The materialized base is retained in
    the output feature and handoff; generated Bash arguments are quoted.
  - the preferred community-facing opener for that story is the
    `VKORC1/rs9923231` hero figure in
    `docs/figures/vkorc1_rs9923231_luciferase_hero.svg`.
  - canonical workflow skeleton:
    `docs/examples/workflows/vkorc1_rs9923231_promoter_luciferase_assay_planning.json`
    (`test_mode: online`, online dbSNP lookup/reference preparation, followed
    by local backbone import and shared promoter reasoning ops). Retained
    reporter SVGs are historical WIP previews, not current-candidate or
    experimentally validated constructs; see their reproducibility README.
  - `docs/tutorial/03-05_gibson_specialist_testing_gui.md`
    - focused end-to-end test script for `Patterns -> Gibson...`
    - uses local committed inputs plus `gibson preview` parity checking
  - `docs/tutorial/03-07_gibson_arrangements_gui.md`
    - opens from a dedicated arrangement-ready starter project where the
      deterministic single-insert Gibson result and stored three-lane
      arrangement are already present, then walks through singleton output
      containers, arrangement inspection, and arrangement-level gel export
  - generated tutorial-project baseline:
    - `docs/tutorial/generated/chapters/03-04_gibson_specialist_testing_baseline.md`
    - preloads stable `gibson_destination_pgex` and `gibson_insert_demo`
      sequence IDs for the manual Gibson specialist walkthrough
    - `Gibson Specialist Starter Project (offline)` also opens the matching
      Help/Tutorial guide automatically
    - `docs/tutorial/generated/chapters/03-06_gibson_arrangements_baseline.md`
    - prebuilds the deterministic `gibson_destination_pgex` +
      `gibson_insert_demo` starter through the canonical single-insert Gibson
      apply so the assembled product and stored arrangement already exist
    - `Gibson Arrangements Starter Project (offline)` now lands on the
      arrangement walkthrough with a distinct arrangement-ready state instead
      of reusing the chapter-15 pre-apply baseline
  - `docs/tutorial/02-03_tp73_cdna_genomic_dotplot_gui.md`
    - retrieve two GenBank sequences and compare them in `Dotplot map` using
      pair modes (`pair_forward` / `pair_reverse_complement`).

Shortcut:

- `Cmd+Shift+G` opens `Retrieve Genomic Sequence...`.
- `Cmd+Shift+P` opens `Prepare Reference Genome...`.
- `Cmd+Shift+L` opens `BLAST Genome Sequence...`.

Protein-evidence behavior:

### Protein Evidence Window

- `Protein Evidence...` opens a specialist window for:
  - online fetch by accession/ID (`FetchUniprotSwissProt`)
  - offline SWISS-PROT text import (`ImportUniprotSwissProt`)
  - projection to selected sequence/transcript (`ProjectUniprotToGenome`)
  - coding-DNA lookup for one mapped protein feature against a stored
    projection (`query_uniprot_feature_coding_dna`)
  - online Ensembl protein/transcript fetch (`FetchEnsemblProtein`)
  - Ensembl protein sequence import (`ImportEnsemblProteinSequence`)
- The dialog shows a compact table of recent imported UniProt entries
  (entry/accession/source/import timestamp) and `Select` buttons to fill the
  active `entry_id`/query in the form.
- The dialog now also shows recent stored UniProt genome projections for the
  active `entry_id` / `seq_id` context and offers:
  - `Use`
    - fill the form with the persisted projection id/sequence pairing
  - `Open Protein Expert`
    - open the shared Protein Expert with transcript-native translation as the
      primary product model and the stored UniProt projection as optional
      external evidence
  - `Open Derived Protein Expert`
    - open the same Protein Expert without requiring a stored UniProt
      projection; the current transcript field acts as an optional transcript
      filter for transcript-only inspection
    - persisted transcript-native protein-derivation reports now also reopen
      through this same route from lineage, using the stored transcript filter
      automatically when the report only derived one protein
  - `Render Derived Protein SVG...`
    - export that transcript-native Protein Expert directly through the shared
      `RenderFeatureExpertSvg` route, again without requiring a stored UniProt
      projection
  - `Render Protein Mapping SVG...` / `Render SVG...`
    - export the same stored projection directly through the shared
      `RenderFeatureExpertSvg` route without first opening the expert window
- the specialist window itself now scrolls vertically, so the UniProt,
  transcript-native, and Ensembl sections stay reachable on smaller viewports
- the same specialist now also exposes external protein feature include/exclude
  fields, reusing the shared `ProteinFeatureFilter` contract so GUI-driven
  Protein Expert open/export paths can filter keys such as `DOMAIN`,
  `DNA_BIND`, or `CONFLICT` without dropping to the shell
- `Open Protein Expert` and `Open Derived Protein Expert` both reuse the
  shared isoform-architecture canvas rather than opening a separate
  protein-sequence inspector:
  - transcript-native CDS/protein derivation now defines the lower product
    geometry whenever GENtle can derive it
  - optional UniProt evidence is layered on top as one external opinion rather
    than the source of truth for which protein products exist
  - the Protein Expert details grid now shows per-row:
    - derived protein length
    - translation table / source
    - translation-speed profile / source / reference species
    - derivation mode
    - external opinion presence/source
    - comparison status
    - short mismatch summary
  - transcript-specific protein lanes now clip to each transcript's projected
    amino-acid coverage, so partial or missing domains remain tied to the
    actual transcript product instead of every lane inheriting the full
    reference annotation set
  - low-confidence disagreements between transcript-native translation and the
    external protein opinion remain visible as explicit warnings/status chips
    instead of being silently flattened away
  - the same compare surface can now also be exported in derived-only form,
    so transcript-first inspection/export no longer needs a stored UniProt
    projection
  - the shared SVG/view layout now deliberately shows both:
    - a coordinate-true exon/intron panel for locus context
    - a lower shared-genomic-exon-column transcript/product panel with
      isoform-local protein axes and shared exon-family colors
    so exon skipping can be read at the product level without losing genomic
    position context
  - those shared CDS colors are now tied to genomic exon family/position, so
    vertically aligned versions of the same exon keep one color across
    different isoforms even when local peptide segment order or exact
    boundaries change
  - when transcript features omit `cds_ranges_1based`, compatible `CDS`
    features are now preferred over whole-exon fallback so RefSeq-style locus
    records preserve transcript-to-protein linkage more faithfully
  - default protein-feature filtering now hides UniProt `CONFLICT`
    annotations from that view/export path while reporting the hidden count in
    the expert warnings block
  - membrane/topology-style UniProt annotations (`SIGNAL`, `TRANSIT`,
    `TRANSMEM`, `INTRAMEM`, `TOPO_DOM`) now render in a dedicated lower band
    beneath the protein rail so sorting signals and topology segments stay
    readable instead of colliding with ordinary domain overlays
  - the Protein Expert window itself now scrolls both vertically and
    horizontally, so wide transcript/protein comparison tables stay readable
    even when mismatch notes and domain inventories get dense
  - future external protein evidence sources should be able to plug into the
    same compare window; Ensembl protein-entry import/compare now also has a
    dedicated GUI section with:
    - `Fetch Ensembl`
    - `Open Ensembl Protein Expert`
    - `Render Ensembl Protein SVG...`
    - `Import Sequence`
    - a recent-entry table whose `Select`, `Open Protein Expert`,
      `Render SVG...`, and `Import Sequence` actions all reuse the same
      transcript-first Protein Expert path and current `seq_id`/transcript
      selection
    - the currently selected Ensembl entry now also shows an inline detail
      panel with transcript/gene/species context, feature-key summary, alias
      preview, and a one-click `Use entry transcript` affordance when the
      active transcript filter disagrees with the stored Ensembl transcript
- The direct SVG-export affordance uses the same persisted projection state and
  target syntax as CLI/shell `render-feature-expert-svg ... uniprot-projection ...`,
  so GUI and non-GUI exports stay on one deterministic rendering path.
- The same specialist now also includes a `Feature coding DNA query` section:
  - `feature`
    - case-insensitive text filter against mapped UniProt feature key/note
  - `feature transcript`
    - optional transcript filter against the stored projection
  - `mode`
    - `genomic_as_encoded`, `translation_speed_optimized`, or `both`
  - `speed profile`
    - `Auto`, `Human`, `Mouse`, `Yeast`, or `E. coli`
  - result panel:
    - requested/resolved translation-speed profile summary with source and
      reference species
    - transcript-specific amino-acid span
    - genomic coding DNA in coding-strand orientation
    - optional translation-speed optimized DNA
    - exon attribution and ordered exon-pair summary when the feature crosses a
      splice junction
  - exon numbering follows transcript order, so reverse-strand transcript exon
    `1` is the transcript 5' exon rather than the lowest genomic coordinate
- The same specialist now also includes a `UniProt projection audit` section:
  - `Audit Mapped Isoforms`
    - runs the integrated shared audit and stores a
      `gentle.uniprot_projection_audit.v1` report
  - `Audit Parity`
    - compares the direct integrated audit against the same result rebuilt
      from the reusable lower-level UniProt/Ensembl primitives and stores a
      `gentle.uniprot_projection_audit_parity.v1` report
  - `audit transcript`
    - optional transcript filter (empty = all projected transcripts)
  - `audit report_id` / `parity report_id`
    - optional persisted report-id overrides
  - recent report tables:
    - recent audits and recent parity reports are listed for the active
      `seq_id` / `entry_id` context
    - `Show` loads one stored report back into the specialist without rerunning
      the operation
  - loaded audit display:
    - per-transcript exon/CDS accounting summary
    - mismatch reasons
    - one read-only local unsent maintainer-email draft when failing isoforms
      were found
- The same specialist now also includes a `Reverse translate protein` section:
  - choose any first-class protein sequence already present in the project
  - optional controls:
    - output sequence id
    - explicit translation table override
    - requested speed profile
    - fast/slow speed mark
    - optional local annealing heuristic target/window
  - result panel:
    - created coding-sequence id
    - protein length to coding-length summary
    - resolved translation-table summary with source and
      organism/organelle context
    - requested/resolved translation-speed profile summary with source and
      reference species
    - selected speed mark and optional annealing heuristic summary
    - inline coding DNA text for the created synthetic CDS
- The same specialist now also includes a `Protease digest protein` section:
  - reuses the selected first-class protein sequence from the protein selector
  - accepts comma-separated protease names/aliases from the built-in catalog
    such as `Trypsin` or `Trypsin,Lys-C`
  - can either materialize first-class peptide sequence products or run in
    report-only prediction mode
  - result panel:
    - resolved proteases
    - cleavage-site and peptide counts
    - transcript/protein provenance when present on the source protein
    - compact cleavage-site preview with local context
    - first peptide-product rows with source amino-acid span, length, created
      sequence id, and peptide sequence preview
- The same specialist now also includes a `DNA handoff reasoning` section:
  - builds an engine-owned construct-reasoning graph for the selected target
    DNA sequence plus active protein source
  - reuses the existing reverse-translation controls instead of introducing a
    second codon-policy surface
  - lets the user choose a ranking goal:
    - `Balanced provenance`
    - `Native fidelity`
    - `Expression optimized`
  - ranks candidate-only handoff strategies rather than creating DNA
    sequences automatically:
    - transcript-native CDS reuse
    - mapped feature-coding DNA from a stored UniProt projection
    - reverse-translated synthetic fallback
  - the candidate table/detail view shows:
    - strategy
    - amino-acid coverage
    - provenance score
    - preserved and relaxed constraints
    - translation-table / speed-profile / speed-mark summary
    - warnings and next-step recommendations
- Stored UniProt projections still act as metadata/projection inputs rather
  than opening their own dedicated protein-only window, but first-class
  protein sequences imported from Ensembl or derived from transcripts can now
  drive reverse-translation workflows from the same specialist.
- Manual GUI sanity-check walkthroughs for these newer protein paths now live in:
  - [`docs/tutorial/06-02_protein_transcript_native_expert_gui.md`](../tutorial/06-02_protein_transcript_native_expert_gui.md)
  - [`docs/tutorial/06-01_protein_reverse_translation_gui.md`](../tutorial/06-01_protein_reverse_translation_gui.md)
- Use one stable `entry_id` in that window when you plan to project repeatedly.
- Persisted UniProt projections now also appear in the project lineage
  graph/table as analysis artifacts linked from the source sequence.
  - lineage actions reopen the UniProt protein expert on the selected
    projection instead of only opening the source sequence
  - lineage details expose the stored `entry_id`, transcript filter, and
    projected transcript count
- Persisted construct-reasoning graphs that contain protein-to-DNA handoff
  candidates now also reopen `Protein Evidence...` directly on that stored
  handoff graph instead of only opening the source sequence generically.

NCBI retrieval behavior:

### GenBank And dbSNP Fetch

- `Fetch GenBank / dbSNP...` opens a combined specialist window for:
  - online GenBank fetch by accession (`FetchGenBankAccession`)
  - optional project sequence ID override (`as_id`)
  - dbSNP rsID resolution plus annotated region extraction (`FetchDbSnpRegion`)
  - `+/- flank bp` control for the extracted region (default `3000`)
  - prepared-genome target selection plus optional output sequence ID override
  - automatic sequence-window open for newly imported or extracted sequence IDs
- Sequences fetched from GenBank can still participate in genome-anchored
  workflows (for example BED/BigWig/VCF projection once anchored). This is a
  complement to GenBank provenance, not a source conflict.
- dbSNP extraction in that window uses the reference-genome catalog/cache
  settings and defaults to full annotation transfer (`annotation_scope=full`,
  `max_annotation_features=0`) so tutorial/reporter workflows can start from a
  fully annotated local locus slice.
- The rsID entry is prefilled with `rs9923231` as a quick-start default for
  the VKORC1/warfarin tutorial path, so opening the dialog and pressing
  `Fetch Region` immediately exercises a real end-to-end tutorial example.
- dbSNP fetch now runs asynchronously in the dialog and streams staged status
  text while it works, including contacting the NCBI Variation refSNP service,
  waiting for the response, parsing the returned JSON, resolving the assembly
  placement, and extracting the annotated slice from the prepared genome.
- Extracted dbSNP regions now also include a one-base `variation` feature
  labeled with the rsID at the resolved position so the focal variant stays
  visible in the sequence view, including overview-scale linear maps. This
  remains a focal anchor for non-SNVs, whose full alleles/coordinates are retained
  in `dbsnp_spdi` rather than misleading single-base VCF qualifiers.
- Marker details distinguish reference match/mismatch/unavailable from
  assembly-family compatibility. Fallback or unbound assembly and missing or
  mismatching reference evidence are warned and block allele materialization;
  a matching base alone does not permit transferring coordinates between
  assemblies. Only clean A/C/G/T SNVs receive `vcf_ref`/`vcf_alt`.
- RefSeq-style dbSNP chromosome accessions such as `NC_000016.10` are matched
  against prepared contig aliases during extraction and annotation transfer, so
  human `chr16`/`16` installs still attach the expected local feature context.
- Shell parity route:
  - `genbank fetch ACCESSION [--as-id ID]`
  - `dbsnp fetch RS_ID GENOME_ID [--flank-bp N] [--output-id ID] [--annotation-scope none|core|full] [--max-annotation-features N] [--catalog PATH] [--cache-dir PATH]`

Resource import behavior:

- `Import REBASE Data...`
  - Parses a REBASE/Bairoch input file and updates
    `data/resources/rebase.enzymes.json`.
  - Loaded project sequences are refreshed so restriction-site tracks use the
    newly imported enzyme set immediately.
- `Import JASPAR Data...`
  - Parses a JASPAR PFM input file and updates
    `data/resources/jaspar.motifs.json`.
  - TF motif lookups for anchored extraction use the refreshed motif registry.
### JASPAR Expert Window

- `JASPAR Expert...`
  - Opens an app-level motif expert for the active local JASPAR registry.
  - The left side now shows a catalog table with motif ID, optional TF name,
    motif length, consensus IUPAC sequence, species names, and collection
    summaries. Every column heading is sortable; selecting the active heading
    again reverses its order.
  - Species metadata already present in the persistent JASPAR snapshot appears
    immediately. Missing metadata for the catalog rows currently rendered on
    screen is fetched automatically through the existing background operation,
    keeping the expert responsive and avoiding an all-registry request.
  - The complete expert surface has always-visible horizontal and vertical
    scroll bars, while the large catalog remains independently virtualized.
  - `Reload catalog` refreshes the local registry table from the shared engine
    while reusing any persisted JASPAR remote-metadata snapshot rows already
    available locally.
  - `Fetch visible species` refreshes and persists remote JASPAR metadata for
    only the currently visible filtered subset, so GUI species summaries come
    from the same snapshot path that shell/CLI/agents can inspect or export.
  - The detail pane shows:
    - a motif overview card (consensus, requested token, motif length, shared
      deterministic background settings),
    - a sequence logo derived from local per-column frequencies/information
      content,
    - the count matrix,
    - non-comparative score-family interpretation blocks that explain the
      best/worst sequences and how often random windows still score above zero,
    - score distributions for multiple score families (`llr_bits`,
      `true_log_odds_bits`) with histograms plus percentile/statistics tables,
    - and optional remote JASPAR metadata/species assignments when REST
      enrichment is requested.
  - The expert stays useful offline because the local matrix/logo/distribution
    view does not depend on the remote fetch succeeding.

Supported in the current flow:

- SnapGene DNA files (`.dna`)
- GenBank files (`.gb`, `.gbk` and similar)
- EMBL files (`.embl`, `.emb`)
- FASTA files (`.fa`, `.fasta`)
  - default interpretation: synthetic blunt `dsDNA`
  - optional FASTA-header metadata tokens for synthetic oligos:
    - `molecule=ssdna` for single-stranded DNA
    - `molecule=rna` for RNA (input `T` is normalized to `U`)
    - `molecule=dsdna` plus optional overhangs:
      - `forward_5=...` (alias `f5=...`)
      - `forward_3=...` (alias `f3=...`)
      - `reverse_5=...` (alias `r5=...`)
      - `reverse_3=...` (alias `r3=...`)
- NCBI sequence XML files (`.xml`, `GBSet/GBSeq`, `INSDSet/INSDSeq`)
  - other unsupported XML dialects are rejected with explicit diagnostics

Example FASTA headers:

- `>oligo_ss molecule=ssdna`
- `>oligo_rna molecule=rna`
- `>oligo_ds molecule=dsdna f5=GATC r5=CTAG`

<!-- manual-page-nav:start -->

---

[GUI Manual](../gui.md) | [Previous: Engine Ops State Persistence](engine-ops-state-persistence.md) | [Next: DNA Window Panels](dna-window-panels.md)
<!-- manual-page-nav:end -->
