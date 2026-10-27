## Array Tracks In DNA Viewer

<!-- manual-counterparts:start -->
Related CLI reference: [Tracks and Arrays](../cli/tracks-arrays.md).
These topic links do not imply identical feature coverage across interfaces.
<!-- manual-counterparts:end -->

- The command palette action `Evidence Preparation` opens a TP73 proof-path
  assistant. It loads the committed anchored TP73 locus, materializes local
  repeat fixtures, inspects/projects the Clariom D manifest, imports the
  CUT&RUN-style BED fixture, annotates TFBS, enables evidence layers, and
  renders proof SVG/report artifacts through the same shared operations used by
  CLI workflows.
- CEL/probeset preparation remains explicit: the assistant shows copyable
  `arrays probe-regions --dry-run`, `arrays run-probe-region-backend PLAN.json
  --allow-external-execution`, and `arrays inspect-probe-region-output`
  commands. Login-walled Thermo Fisher support ZIPs are only listed as expected
  local inputs; GENtle does not download them or run R package installation
  from the GUI. The same assistant also offers compact review-gated buttons for
  the explicit local-file workflow:
  `arrays import-apt-probe-region-output`,
  `arrays inspect-probe-region-output`,
  `arrays render-probe-region-output-svg`, and
  `arrays project-probe-region-output`; these buttons use the shared shell
  executor rather than GUI-local probe-region logic.
- Prepared Clariom D array output can be projected with
  `arrays project-microarray-track`. The command validates that the manifest
  coordinate system matches the open sequence's genome anchor, or that the
  manifest declares an explicit coordinate projection into that anchor build,
  before adding any features.
- Projected rows become ordinary `track` features with
  `gentle_track_source=Array` plus dataset, platform, contrast, `logFC`,
  `adj_P_Val`, probeset, transcript-cluster, exon, and optional junction
  qualifiers.
- Live feature-tree array groups include the `gentle_array_dataset` identity,
  so freshly generated `probe_region_output` evidence remains visibly distinct
  from prepared publication datasets even when platform and contrast match.
- The linear DNA viewer has an `Array` layer toggle. Array intervals draw as
  contrast heat lanes in the coordinate-aligned feature stack rather than in
  the detached regulatory band at the top of the canvas: red for positive
  `logFC`, blue for negative `logFC`, and grey near zero. Stronger colour means
  larger absolute `logFC`; in the live GUI, the outline/darkness also indicates
  `adj.P.Val <= 0.05`. The linear-map legend and SVG export legend repeat this
  rule so screenshots and exported figures remain interpretable outside GENtle.
- Junction array probes draw as two separated probe arms connected by a
  purple bridge over the skipped interval. When `junction_start_edge` and
  `junction_stop_edge` are available, the bridge is anchored to those
  exon-boundary coordinates; otherwise a junction-typed array row falls back to
  an equal split of its projected interval so it is still visually distinct
  from an ordinary exon-overlap probe.
- Selected-feature details and copied tooltip text include all projected values
  stored for a probeset via `gentle_array_value_summary`. The feature tree
  accepts scoped filters such as `source:array`, `track:Clariom`,
  `contrast:AdTAp73alpha-AdGFP`, and `gene:TP73`.
- Array details also expose the native array coordinate system, the open
  sequence anchor genome, any build-projection method/status, and both native
  and displayed genomic intervals for assembly cross-checks.
- Array feature context menus include copying the array value table and
  focusing matching probesets across all contrast lanes.
- The shared GUI Shell and command reference expose
  `arrays probe-regions` as a read-only preflight for arbitrary CEL files or
  publication-resource datasets. The JSON plan reports CEL, metadata,
  condition-contrast, annotation/library, platform, backend-candidate, local
  dependency, output, and cache-readiness checks. Users can add repeatable
  `--r-library-path PATH` flags in that shared command surface when packages
  live in an agent sandbox, user library, or system-specific tree. The plan
  reports both the requested roots and R's effective library paths; dependency
  errors repeat those paths and ask the user to verify the flag when GENtle's
  result differs from an interactive R session. The `r_oligo` candidate
  advertises the generic `scripts/probe_regions_oligo.R` helper and includes
  an advisory command for explicit RMA/CEL requests; the panel can run that
  selected backend only through the shared
  `arrays run-probe-region-backend PLAN.json --allow-external-execution`
  capability after the user enables the external R/APT confirmation checkbox.
  Missing preflight dependencies stop before R/APT is launched, and GENtle
  still does not download or install packages. For `r_oligo`, the externally
  installed platform package normally bundles its own SQLite annotation
  database. A probeset run stops rather than emitting coordinate-empty output
  if that package-local database cannot be resolved unambiguously;
  transcript-cluster-only output instead continues with a visible warning.
- For `Clariom_D_Human`, the same preflight lists the local Thermo Fisher
  na36 hg38 support ZIP paths under
  `annotation_source.vendor_support_files[]`; these login-walled files are
  manually staged under
  `data/resources/affymetrix/clariom_d_human_na36_hg38/`, never
  auto-downloaded by GENtle. Both concise canonical ZIP names and
  browser-preserved `TFS-Assets_LSG_Support-Files_...` download names are
  accepted.
- The separate sequence-tools window includes a native
  `Clariom D / probe-region evidence` inspector for completed
  `probe_regions_oligo.R` output directories and explicit APT summary +
  annotation table imports, with optional sample metadata columns to derive
  condition means/SDs and default log2FC tracks. When explicit PM probe
  coordinates are present in the annotation table, the same import can write
  and inspect a probe-coordinate table. Supplying an explicit probe-intensity
  matrix lets that table carry true PM probe-level values marked as
  `probe_level_input`; otherwise the GUI keeps the parent-probeset summary
  fallback marked as `parent_probeset_summary`. It uses the same shared
  `arrays import-apt-probe-region-output` and
  `arrays inspect-probe-region-output OUTPUT_DIR` contracts as the GUI shell,
  then previews helper outputs, sample/condition/logFC columns, chromosome/gene
  coverage, coordinate/build declarations, projection blockers, and a bounded
  chromosome-ordered row preview. The same panel can export a deterministic
  native SVG plot from `mean_log2_*` and `log2FC_*` helper-table columns, then
  project direct-compatible or explicitly mapped helper output into
  genome-anchored array features through the shared
  `arrays project-probe-region-output` route. Projection level is explicit:
  `probe_region` keeps the existing probeset/region behavior, while
  `pm_probe` projects only true probe-level rows marked as `probe_level_input`.
  After projection, the panel can run `arrays interpret-probe-region-evidence`
  and preview the resulting engine-owned evidence/transcript geometry report
  without invoking R/APT or changing the projected features. The preview
  includes compact transcript-mapping summaries for exon overlap and
  junction-spanning evidence; the full JSON report keeps the structured exon
  ordinals, ranges, junction spans, overlap base counts, geometry scores, and
  score-basis guardrails. Transcript rows show review-only labels for unique,
  shared, constraining, or absent geometry; these labels are not isoform calls.
- The inspector starts with contextual file-role help and uses the same bounded
  checkout discovery as `Configuration -> Microarrays`. When available it
  offers an explicit TP73 `arrays probe-regions` preflight command for locally
  staged E-MTAB-14704 CEL files, or direct selection of the committed synthetic
  TP73 helper output. Full-platform analyzed tables are shown with their exact
  paths but are not blindly filled into the direct-import form because they
  contain coordinate-less platform rows that require selector-aware preflight.
  Suggestions never overwrite a path without a click, and the synthetic output
  is marked as learning/test data rather than biological evidence.
- The DNA-viewer layer control distinguishes configuration from visibility.
  When no projected microarray feature exists on the active sequence it reads
  `Array setup...` and opens the separate Sequence tools window scrolled
  directly to `Clariom D / probe-region evidence` rather than toggling an empty
  layer or leaving the user at the top of a long workspace. The panel visibly offers
  the built-in `Clariom_D_Human` / `rma-sketch` / `hg38` (`GRCh38`) profile and
  a restore-defaults action plus a direct link to the general Microarrays
  configuration/help tab. Once data is projected, the control reports both
  visible and full-sequence counts when they differ, for example `Array
  (0/14)` after panning away from all array intervals.
- Right-click the Array control and choose `Focus nearest array feature` to
  navigate from the current selection midpoint, otherwise the visible-span
  midpoint, to the nearest projected array interval. This uses the shared
  engine feature query rather than a GUI-only lookup. `Open Clariom D / array
  setup` remains available from the same context menu even when array data is
  already present.

<!-- manual-page-nav:start -->

---

[GUI Manual](../gui.md) | [Previous: TATA-Box Evidence In The DNA Viewer](tata-box-evidence-in-the-dna-viewer.md) | [Next: TFBS Annotation (Engine Ops)](tfbs-annotation-engine-ops.md)
<!-- manual-page-nav:end -->
