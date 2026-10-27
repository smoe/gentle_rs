## Primer, qPCR, and Sequence-Specific RT Design

<!-- manual-counterparts:start -->
Related CLI reference: [Primer Design and Handoff](../cli/primers.md) | [Reference-Bound Primer Specificity](../cli/primer-specificity.md) | [Transcript Assay Panels](../cli/transcript-assays.md).
These topic links do not imply identical feature coverage across interfaces.
<!-- manual-counterparts:end -->

The Engine Ops panel includes `Primer and qPCR design reports` for explicit
`DesignPrimerPairs` and `DesignQpcrAssays` execution on the active sequence.
The dedicated PCR Designer additionally exposes transcript-panel and
terminal-exon RT-primer-pool modes without moving their biology into the GUI.

Layout:

- `Sequence tools…` opens sequence-specific configuration and Engine Ops in a
  separate resizable window instead of expanding the DNA-view toolbar
- `Shell…` opens the shared shell in that same tools window; horizontal and
  vertical scrollbars keep wide operation forms usable without shrinking the
  DNA map
- genome-anchor extension and repeat-index inputs also live in Sequence tools;
  the DNA view retains only compact anchor status and a `Configure…` launcher
- `Primer and qPCR design reports` is opened by default so the primer/qPCR forms
  appear first in the tools workspace
- opening a lineage sequence node reports `Opening DNA Sequence Viewer…` in the
  status bar until the queued viewer completes its first render
- first-frame fitting of lineage panels to the current window is read-only and
  does not mark a freshly opened project as modified

Primer backend/preflight controls:

- dedicated `Primer backend and Primer3 preflight` block above the design forms
- configurable fields:
  - `backend` (`auto` / `internal` / `primer3`)
    - `auto` prefers Primer3 when reachable and otherwise falls back to the
      internal backend; reports show `requested->used`
    - `internal` never calls Primer3
    - `primer3` requires a reachable executable and fails when it is missing
  - `primer3 executable` (path, default `primer3_core`)
- actions:
  - `Apply Primer Backend` (persists `primer_design_backend` and
    `primer3_executable` engine parameters)
  - `Probe Primer3` (runs availability/version probe and reports status)
  - `Reload from Engine` (refreshes controls from current engine parameters)
- preflight status line reports configured/effective executable,
  resolved path, working directory, backend, reachability, bounded-progress
  support (`supported`, `legacy/no`, or `unknown`), and probe timing.

Primer pairs form:

- simple-starter block:
  - `Use current selection as simple PCR core`
  - `max primer distance from core`
  - `Apply simple flank windows`
  - writes deterministic side-window bounds:
    - forward: `[ROI start - D .. ROI start]`
    - reverse: `[ROI end .. ROI end + D]`
  - enables `require ROI flanking`
  - keeps `min amplicon` at least the core length
  - saved-report preview now shows top accepted pairs in beginner terms:
    - amplicon coordinates/length
    - left/right distance from the core ROI
    - whether the pair cleanly flanks the core
    - Tm delta and rank score
    - the selected pair's score model/direction and compact additive
      contribution table
    - bounded rejected near misses with retained/eligible/omitted capture
      counts and completeness reason
    - the report-owned excluded-region analysis status/reason, including an
      explicit `not_run` when no region-level evidence source was consulted
    - an action that focuses the report-linked construct-reasoning overlay
- ROI and amplicon controls:
  - `ROI start`, `ROI end`
    - both accept numeric coordinates and `=` formulas
      (`=KIND.start+N`, `=KIND.end-N`, `=KIND.tss`,
      `=KIND.upstream(N)`, optional `KIND[2]`, optional
      `KIND[label=TP73]`)
    - `ROI start` also accepts range form
      (`=left .. right` or `=left to right`); range form resolves the interval
      between both coordinates, including reverse-strand upstream/TSS pairs
    - `Apply ROI formula` resolves formulas into numeric coordinates
  - `min amplicon`, `max amplicon`
    - pair-PCR now requires `min amplicon <= ROI length`
  - `max Tm delta`, `max pairs`
  - `report_id`
- side constraints (`Forward side`, `Reverse side`):
  - core fields: length bounds, optional location/start/end, Tm bounds, GC bounds,
    `max anneal hits`
  - default length bounds in GUI forms: `20..30`
  - sequence constraints:
    - `5' tail (non-annealing)` (e.g. restriction site/adaptor prefix)
    - `fixed 5'`, `fixed 3'`
    - `required motifs`, `forbidden motifs` (comma-separated IUPAC motifs)
    - `locked positions` as `offset:base,offset:base`
- pair constraints:
  - `require ROI flanking`
  - `fixed amplicon start`
  - `fixed amplicon end (exclusive)`
  - required/forbidden amplicon motifs
- PCR region queue block:
  - row source + template + `start/end/len` columns
  - `Use ROI` per row to copy queued coordinates back into active ROI fields
  - queued batch PCR also requires each queued ROI length to be at least the
    configured `min amplicon`
  - `Remove` per row and `Clear queue`
  - batch toggle: `Also create extracted region copies`
  - batch run: `Design Primer Pairs for queued regions`
    - runs one `DesignPrimerPairs` operation per queued row
    - report IDs derive from current primer `report_id` base with `_rNN` suffix
    - optional copy mode emits one `ExtractRegion` artifact per queued row
- PCR batch results block:
  - one row per queued-region run with status, region coordinates, report ID, and optional copy ID
  - quick actions per row:
    - `Show` (report summary)
    - `Export` (report JSON)
    - `Open` (opens extracted copy when present, otherwise template)
  - `Clear results` removes recorded batch rows without deleting persisted reports
- local specificity confirmation block:
  - uses the current primer `report_id` plus a 1-based pair rank
  - accepts a prepared reference genome id, max product length, and max BLAST
    hits per primer
  - `Confirm specificity` calls the shared
    `AssessPrimerPairSpecificity` engine operation, persists its report, and
    reports the stable report id, pass/fail state, intended/unintended product
    counts, primer-hit counts, and warning count in the status line
  - the confirmation block retains the latest persisted specificity-report id
    and shows a compact target/status/product/design-provenance summary
  - only annealing segments are BLASTed; saved 5' tails remain provenance in
    the shared report contract
  - the search does not use BLAST `-max_target_seqs`; the configured maximum
    hit count is a post-search review threshold, so auxiliary contigs remain in
    scope
  - intended genomic products use explicit genomic target geometry rather than
    cDNA amplicon length, and junction-spanning assays may validly report no
    contiguous intended genomic product
  - report JSON keeps genomic carryover/off-target assessment separate from
    transcriptome/cDNA cross-amplification assessment and records the validated
    database content identity

qPCR form:

- same ROI/amplicon + pair constraints as primer-pair design
- same ROI formula support as primer-pair design (`=` expressions +
  `Apply ROI formula`)
- default qPCR constraints are probe-based/TaqMan-like: 80-200 bp amplicons,
  18-24 bp primers at 63-73 C with pair delta <= 3 C, and a 20-30 bp probe at
  71-80 C
- adds a `Transcript targeting` block above ROI/probe controls:
  - `Genomic`
    - keeps the existing direct qPCR flow and emits
      `transcript_targeting = null`
  - `Shared across transcripts`
    - requires a splicing-group context from the same sequence
    - reuses the current Splicing Expert group on
      `TargetGroupTargetStrand`
    - warns that GENtle may fall back to the largest supported transcript
      subset when no fully shared exon/exon-chain span exists
  - `Specific transcript`
    - requires a splicing-group context plus the currently selected transcript
      from Splicing Expert
    - exposes `Specificity evidence`:
      - `Junction only`
      - `Unique exon / exon chain only`
      - `Either (prefer junction)`
    - the GUI does not ask the user to choose forward-vs-reverse junction side;
      GENtle selects that automatically
  - transcript-aware modes now show:
    - read-only transcript context summary
    - `Use current Splicing Expert context` / refresh action when a compatible
      splicing view is already open
    - `Open / Focus Splicing Expert`
    - inline validation when transcript-aware qPCR has not yet been seeded from
      a compatible splicing overview
- adds:
  - `Probe side` with the same side-constraint fields
  - `max probe Tm delta`
  - `max assays`
- the qPCR controls mirror the primer-pair hover help and use wider ROI
  coordinate fields so large genomic positions remain easy to enter
- in the dedicated `PCR Designer`, qPCR is now a first-class mode rather than
  an Engine-Ops-only appendix:
  - `Pair PCR | qPCR | Transcript panels | RT primer pool | Gene assay study` selector switches the right-hand
    constraints/run panel between the three workflows
  - the left-hand paint/selection workflow stays shared so the same painted or
    formula-defined ROI can feed either pair-PCR or qPCR without switching
    windows
  - lineage/report reopen on qPCR reports now lands on the qPCR mode
  - the left-hand paint/selection ROI workflow is shared, while queued
    batch-region helpers remain pair-PCR-only for now
  - Splicing Expert now also offers direct qPCR launches:
    - `Design shared-transcript qPCR`
    - `Design transcript-specific qPCR`
    - both seed the shared PCR Designer with ROI plus transcript-aware qPCR
      context instead of only sending a generic ROI
  - qPCR mode now also exposes the built-in protocol-cartoon route directly:
    - `Export qPCR Cartoon SVG...` writes the shared `pcr.assay.qpcr` strip
      without needing a separate shell/CLI round-trip
    - the same built-in protocol id remains usable from shared shell / direct
      CLI as `protocol-cartoon render-svg pcr.assay.qpcr OUTPUT.svg`
  - qPCR mode now also shows:
    - a selectable saved-report assay preview with forward / reverse / probe
      geometry for the top ranked assays
    - the selected assay's `gentle_qpcr_assay_rank_v1` contribution table,
      capture completeness, and bounded rejected assay near misses
    - the report-owned excluded-region analysis status/reason
    - an action that focuses the report-linked qPCR selection-reasoning graph
    - the currently selected saved assay now drives the live qPCR
      cartoon-geometry summary instead of always using assay rank `#1`
    - when a matching splicing expert view is active on the same template, the
      assay preview also annotates qPCR candidates as single-exon,
      junction-spanning, or junction-crossing and explains that context in
      plain language
    - when transcript-aware targeting was requested, saved reports now also
      show the requested targeting mode, requested evidence kind, realized
      evidence kind, and selected transcript/support summary without needing a
      live splicing window
    - the live qPCR cartoon-geometry summary is still derived from the selected
      saved report when available, or from the current ROI + constraint
      defaults otherwise
  - primer and qPCR report previews cache the current report revision and
    invalidate that cache when a design completes, window state is restored,
    or a report is explicitly reopened.
- `RT primer pool` is a thin GUI over the shared
  `DesignTerminalExonRtPrimerPool` operation:
  - enter the fixed 5-prime adapter, exact sequence-specific segment length,
    terminal-exon search window, retained-candidate limit, and report id
  - open a Splicing Expert group and use `Add current transcript` to append its
    selected transcript as the next priority, or enter rows manually as
    `SEQ_ID<TAB>FEATURE_ID<TAB>TRANSCRIPT_ID<TAB>LABEL`
  - `n-N` accepts the 1-based feature label shown by the GUI; a bare feature id
    is the engine's 0-based id. Row order is the explicit pool priority
  - `Design RT primer pool` runs in the existing detached background-operation
    path. The GUI does not infer terminal exons or calculate oligo sequences
  - the report view shows the selected 22-mer (or requested exact length), the
    complete adapter-bearing oligo, strand-aware terminal-exon/source geometry,
    descriptive GC/Tm, retained alternatives, and selected-pool pairwise
    complementarity
  - `Copy selected oligos TSV` copies the selected sequence text plus target,
    transcript, and source-coordinate context;
    `Open saved report` and `Export report JSON...` use the engine-owned report
    store and export path
  - Tm is visible but deliberately excluded from this workflow's ranking.
    Whole-transcriptome/genome specificity is also not implied by this report
- `Gene assay study` is a session-scoped mode in the same PCR Designer, also
  opened by Splicing Expert's `Gene assay study...` action:
  - inspect a typed saved `gentle.gene_isoform_assay_study_plan.v1` for the
    current source sequence, its declared evidence/missing assessments,
    transcript scope, recommendation/override and exact operation digests
  - inspect persisted panels on that sequence, explicitly distinguishing
    comparison panels from plan-bound results; selected-pair views show both
    primer footprints on mature cDNA and each transcript's stored products.
    The interval table retains explicitly labelled 0-based half-open bounds;
    the plots use explicitly labelled 1-based inclusive positions.
    The selected-pair axis marks the exact 1-based displayed start/end of the
    predicted amplicon with a filled green product span, while blue/orange
    arrows retain forward/reverse orientation. Primer arrows may use a minimum
    on-screen width for legibility, but the green product boundaries remain
    exact. A second, explicitly labelled expanded detail maps the primers
    within those same exact amplicon boundaries, without pretending to share
    the full-transcript scale. Both views are mature-cDNA projections, not an
    inferred genomic overlay
  - the coverage-scope box uses the experimental handoff's engine-owned
    accounting: annotated records, distinct mature cDNAs, uncovered records
    and records outside the assessed universe remain separate. Identical
    mature cDNAs cannot be distinguished by sequence-based primers; a completed
    design objective is not a specificity pass
  - selected-pair outcomes, the transcript-by-assay matrix and endpoint band
    table have independent page controls. All stored rows/assays are accessible,
    rather than truncating at 40 transcripts, 24 assays or 80 band rows. Paging
    does not change coverage, selections or exports; missing matrix cells say
    `Not assessed`, not no product
  - load a planning request, edit its short-product hard ceiling or full JSON,
    normalize through the engine, review it before planning, then separately
    review the exact workflow before execution; edits/project changes invalidate
    review, and imported plans do not inherit execution authority
  - managed background tasks expose progress and cancellation; stale results
    cannot replace project state. New output directories and panel identities
    preserve prior results; external file effects are not rolled back
  - canonical dossier export takes the same explicit hash-bound publication
    request as CLI/OpenClaw. Missing handoffs stay pending; no automatic order
    submission or alternate GUI scientific narrative is introduced
  - walkthrough: [04.08](../tutorial/04-08_gene_assay_study_gui.md). Live GUI and
    real-data acceptance are pending; the synthetic fixture is not a human assay
- `Transcript panels` is a thin GUI over the shared
  `DesignTranscriptAssayPanel` operation:
  - Splicing Expert's `Design all-transcript panel` action supplies the source
    feature/group and opens the PCR Designer in the new mode
  - endpoint RT-PCR, primer-only SYBR qPCR, and TaqMan qPCR are selectable
    without changing the engine defaults used by older callers
  - strict `Require all` coverage is the GUI default; `Best effort` remains an
    explicit choice that preserves uncovered cDNA classes in the report
  - a default-off `Pre-approve informative partial fallback` block wraps the
    exact strict one-per-class/minimal-discrimination operation in the shared
    engine contract. The user approves a whole-panel assay budget, minimum
    resolvable band-size difference, nonzero-increment rule, optional long-
    range permission, and distinct fallback report id before starting the
    strict run. Only typed coverage infeasibility can activate it; other errors
    remain fail-closed. Saved fallback panels are visibly partial and show the
    parent failure/policy identities, objective totals, uncovered classes, and
    combined-gel disposition
  - an independent experimental-tier control distinguishes routine common-
    region screening, isoform discrimination, and long-range structural
    discovery; optional preferred product limits sit inside the existing hard
    allowed range
  - Clariom JUC/probe-evidence reports can be supplied as junction-target paths,
    with required/preferred priority and Primer3 overlap controls
  - the saved-report view shows the transcript-by-assay product matrix,
    endpoint band sizes, junction evaluations, short SYBR junction assays,
    candidate primer sequences (not order approval), and the per-cell oligo-dT
    5-prime reach status
  - an empty oligo-dT threshold reports annotation-derived distances without
    making a categorical risk call; a user-supplied threshold enables the
    corresponding within/elevated matrix labels
  - the saved-report rationale shows annotation-confirmed common-region status,
    separate PSR/JUC support, routine-versus-fallback classification, the
    concise selection explanation, and bounded rejected alternatives
  - `Build experimental handoff` runs the shared
    `BuildExperimentalAssayHandoff` operation with the documented default
    readiness policy and displays one collapsible card per selected pair,
    including oligos, gate outcomes, blockers, and predicted cDNA products.
    Missing specificity remains a visible default blocker, while absent variant
    evidence is shown without being silently treated as a pass
  - `Copy advanced command` provides the equivalent
    `primers experimental-handoff PANEL_REPORT_ID ...` shell route for custom
    policy, variant evidence, order-form linkage, and JSON/TSV export; that same
    operation remains available to CLI/MCP/JavaScript/Lua

Buttons:

- `Design Primer Pairs`
- `Design Primer Pairs for queued regions`
- `Confirm specificity`
- `Design qPCR Assays`
- while the internal backend is running, PCR Designer now shows a live
  progress line for:
  - candidate enumeration
  - pair-combination evaluation
  - qPCR probe-placement screening
  - completion counts
- while a progress-capable Primer3 backend is running, the same area shows the
  current input-record number and bounded candidate work as
  `completed/bound` plus a percentage. This percentage reflects instrumented
  search work, not elapsed time or an ETA.
- GENtle detects `--progress` through Primer3's help output and can request a
  refreshed status after a confirmed progress stream becomes quiet. Older
  Primer3 installations run directly in legacy mode, but they cannot provide
  the bounded-work display.

Beginner tutorial:

- [`docs/tutorial/04-01_simple_pcr_selection_gui.md`](../tutorial/04-01_simple_pcr_selection_gui.md)
- primer-report helpers (uses current primer `report_id` field):
  - `List Primer Reports`
  - `Show report_id`
  - `Export report_id...`
  - `Confirm specificity`
  - in-panel `Primer report preview`
- qPCR-report helpers (uses current qPCR `report_id` field):
  - `List qPCR Reports`
  - `Show report_id`
  - `Export report_id...`
- transcript-panel report helpers (uses current transcript-panel `report_id`):
  - `List reports`
  - `Show report_id`
  - `Export report_id...`

All design modes persist reports into project metadata through the same stores
used by CLI/shared-shell `primers ...` commands.

Restriction-site cloning handoff:

- lives inside the PCR Designer as a post-design block, not as part of core
  primer proposal/ranking
- uses one persisted primer-pair report + selected pair rank as input and
  creates a cloning-facing handoff without modifying the original short-primer
  report artifacts
- controls:
  - `pair rank`
  - `destination vector`
  - `mode` (`single_site` / `directed_pair`)
  - `forward enzyme`, `reverse enzyme`
  - optional `forward 5' leader`, `reverse 5' leader`
  - `Vector shortcuts`
  - one-click `Recommended single-site` and `Recommended directed pairs` from
    the same MCS-first / unique-cutter ordering used by CLI/shared-shell
    reasoning helpers
  - `Suggested enzymes` fallback pickers
  - `Create restriction-tail handoff`
- behavior:
  - handoff creation now uses the same engine-owned seeded request helper as
    the CLI/shared shell, so GUI validation and auto-selection follow the
    exact same rules
  - if enzyme fields are left blank, the top recommended single-site or
    directed-pair choice is applied automatically during handoff creation
  - creates new extended forward/reverse primer artifacts plus one predicted
    tailed amplicon artifact and one per-handoff pool container
  - keeps annealing Tm semantics tied to the original annealing segment while
    recomputing full-oligo self/dimer heuristics for advisory warnings
  - blocks incompatible cloning setups when vector sites are absent/non-unique,
    insert-site counts imply internal collisions, or directed MCS order does
    not match the chosen forward/reverse site order
- saved-handoff helpers:
  - `List saved handoffs`
  - `Show saved handoff`
  - `Export saved handoff...`
- preview shows:
  - extended primer sequences
  - predicted tailed amplicon ends/length
  - vector cut positions and compatibility verdict
  - a staged pre-ligation workflow for the runnable `PcrAdvanced` + `Digest`
    steps, plus the remaining ligation placeholder hint
- vector enzyme suggestions prefer `mcs_expected_sites` from MCS-like features
  and fall back to unique cutters by cut position when no explicit MCS order is
  available

Persisted primer/qPCR design reports and primer-specificity reports now also
appear in the project lineage graph/table as analysis artifacts linked from the
template sequence.

- lineage actions reopen the PCR Designer on the selected report instead of
  only opening the source sequence
- clicking the underlying PCR-related operation glyph or `Op` cell reopens the
  same PCR Designer on the originating template sequence
- lineage details expose the stored backend plus pair/assay counts
- primer-specificity lineage details expose target database kind, result
  status, candidate-product count, and failing off-target count; reopening
  restores its compact evidence summary and cited primer-design report when
  available
- restriction-cloning handoffs also appear as lineage analysis artifacts linked
  from the template sequence and reopen the same PCR Designer with the saved
  handoff highlighted

`Design Primer Pairs` materializes lineage-visible sequence artifacts for each
accepted pair:

- one forward primer sequence (`..._fwd`)
- one reverse primer sequence (`..._rev`)
- one predicted amplicon sequence (`..._amplicon`), including configured
  non-annealing 5' tails
- one per-pair pool container containing all three artifacts

Insertion-first anchored-tail mode (MVP status):

- engine contract is available as `DesignInsertionPrimerPairs` (same persisted
  primer report schema plus `insertion_context` compensation rows)
- dedicated GUI form is not wired yet
- current GUI workaround:
  - set forward/reverse `5' tail (non-annealing)` to extension sequences
  - set forward/reverse `start/end` to anchor-adjacent primer windows
  - run `Design Primer Pairs`
  - use shell/CLI `op` payload when you need explicit `insertion_context`
    shift-compensation reporting in the stored report

Heuristic guidance surfaced in reports/warnings:

- primer ranking prefers `20..30 bp` side lengths and similar primer Tm
- 3' GC clamp (`G/C`) is preferred but not hard-required
- penalties are applied for homopolymer/self-complementary runs and primer-dimer risk

<!-- manual-page-nav:start -->

---

[GUI Manual](../gui.md) | [Previous: Design Constraints Filter (Engine Ops)](design-constraints-filter-engine-ops.md) | [Next: Sequencing Confirmation](sequencing-confirmation.md)
<!-- manual-page-nav:end -->
