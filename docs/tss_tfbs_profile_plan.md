# Accession-Pinned TSS TFBS Profile Documents

Status: core implementation and focused verification complete; Glen's committed-
producer five-gene replay and release acceptance are not yet signed off.
Prepared: 2026-09-10. Inspection baseline: local `main` at `45a9f574`;
the relevant Rust sources match development HEAD `76ac480a`.
This proposal does not change the `.10` release gate. Schedule new capability
work separately unless the release owner explicitly brings it into scope.
Claude supplied the initial W1-W10 decomposition; this revised plan has not
received a separate Claude review.

## Revised Input Programme, 2026-09-10

The later Claude/Glen work programme identifies input revision
`a106cbbd5223f8c4be55a7d846b8416bc4b3ae33`. It is available as a Git object but
is not an ancestor of this checkout. Read its named inputs without merging that
branch or changing historical figures. Implementation was rebased onto local
`main` at `a0574f96b0b9e700ec0b2918586da3ba295f5398`.

The following clarifications supersede the original missing-input assumptions:

- Read the existing `gentle.target_tss_fasta_export.v1` format. Its manifest
  binds `SHA256SUMS`, and both bind each FASTA; the outer profile receipt binds
  the manifest bytes and its `source_revision`. Require independently supplied,
  whole-string genome, assembly and dataset expectations. Do not derive a release
  by parsing words from a genome name; absent separate release metadata stays null.
- Read the real panel's `tracks[]` with per-track policies. The generic parsed
  input can express mixed score kinds; this TSS computation/rendering path must
  reject mixed score/clip grammar and any unsupported settings explicitly.
- Join the existing regulatory selection's `promoterome_id` to `promoter_id`.
  Require matching gene, chromosome, TSS and strand, and subset transcript
  membership. Selection and display windows deliberately have different lengths
  and sequence hashes; those are not join keys. FASTA/manifest transcript sets,
  in contrast, must be equal. Carry the supplied CUT&RUN factor/criterion into
  labels and conservative selected-panel legends, never a hard-coded factor.
- Steffen explicitly confirmed retaining identical DNA at distinct locations
  with a warning. This overrides the revised W2's unconditional rejection.
- Keep the new report layout in `gentle-render`, without rewriting the historical
  locus renderer. Use the existing in-process PNG/PDF path with actual used-font
  digests; do not invent external-renderer commands or require a PATH executable.

Read-only inspection confirms the input manifest lists 31/6/4/5/12 TSSs and
109/17/15/13/20 distinct transcripts for CD44/TGFB1/SERPINE1/PATZ1/TP73, the
panel has 30 tracks/28 factors, and the selection contains 13 records. These are
input facts, not evidence that the new full scientific replay has passed.

## Verification Snapshot, 2026-09-10

The implementation diff on rebased main `a0574f96` was checked on macOS arm64
with Rust `1.99.0-beta.4` (`948ec4e1e`) and Cargo `1.99.0-beta.4` (`5f94df478`).
These are development checks, not a frozen-release acceptance receipt:

| Check | Result |
| --- | --- |
| Root headless `tss_` suite | 92 passed; two explicit opt-in tests ignored |
| `gentle-protocol` complete library suite | 68 passed |
| `gentle-engine` complete library suite | 13 passed |
| `gentle-render` complete library suite | 126 passed; visual-output helper ignored |
| Raw panel parser/resolution | 8 passed, including all 30 real matrix accessions |
| Opt-in real bundle/selection reader | Passed: 58 TSSs and all 13 selections |
| PNG/PDF helpers | 10 / 6 passed, including legacy pixel/byte compatibility |
| MCP module | 39 passed, including consent and no-write failures |
| Existing motif resolver / score-track renderer | 11 / 4 passed |
| Generated GUI/CLI/MCP parity matrix | Freshness test passed |
| Default-feature `cargo check --locked --offline -q` | Passed |
| Formatting and whitespace checks | Passed |

Root tests used `CARGO_PROFILE_TEST_DEBUG=0`, `CARGO_INCREMENTAL=0` and one
compiler worker. Direct test-binary runs kept the repository's 16 MiB
`RUST_MIN_STACK` setting. Parity freshness ran the unchanged integration test
with the freshly compiled protocol library, without building every root binary.
The tiny typed-operation smoke generated all three visual formats and verified
exact report-only replay after deleting its copied source inputs. Synthetic
plus/minus/layout examples were inspected; they do not validate real promoters.
The old score-track renderer's bytes are unchanged from the rebase baseline.

Stages A-E are implemented through the shared operations and GUI Shell. Stage F
remains open: no standalone CLI binary smoke, full workspace/tutorial/release
suite, native GUI acceptance or full real five-gene scoring/render replay is
claimed by these development checks. The source formats, adapter differences,
input policy override and in-process raster choice are documented above and in
the [workflow guide](tss_tfbs_profiles.md); no historical figure was regenerated.

## Objective

Produce a generic, engine-owned document for each gene's distinct annotated
TSS windows, with an explicitly ordered, accession-pinned JASPAR panel.
Show transcript-oriented score profiles, readable scales and provenance, and
compare alternative matrices for the same factor without choosing a biological
winner. JSON, TSV, SVG and optional PNG/PDF must describe the same computation.

The five-gene, 30-matrix example is an acceptance dataset, not a runtime preset.
Historical tall locus reports, their evidence lanes and their PDFs stay intact.
New links from those reports are a separate follow-up.

## Verified Starting Point

- `src/engine/analysis/promoter_design.rs` already owns continuous scoring,
  deterministic background calibration, peaks, and Pearson/Spearman summaries.
  Its public helper currently expands motif queries and resolves matrices
  again while scanning. It cannot directly satisfy a frozen strict panel.
- `src/tf_motifs.rs` supports an active runtime registry, bundled registry and
  supplemental PFMs. It also supports aliases and consensus-only fallback.
  Hashing only `assets/jaspar.motifs.json` would not necessarily identify the
  actual matrices used by a run.
- `src/engine/protocol.rs` carries one score kind per score-track report.
  It defines each score index as the start of a window on the local sequence,
  for both motif orientations. Display clipping can already alter the arrays.
- `src/render_tfbs_score_tracks.rs` has a shared range, motif logos and TSS
  markers. Its polyline x positions currently span the width of each vector,
  rather than deriving positions from the common sequence coordinates.
- `src/engine/analysis/feature_expert_ops.rs` already checks typed calibration
  identities for shared regulatory scales. Reuse its rule, not prose matching.
- `src/svg_png.rs` uses in-process resvg; `src/svg_pdf.rs` produces single-page,
  raster-backed PDFs through that path. An external renderer is not required.
- Existing promoterome preparation retains unique genomic windows and separate
  transcript memberships. Reuse its identity/geometry semantics, not a second
  definition of a promoter.

The exact inputs were unavailable during the initial plan. The revision above
now locates them and resolves Q1-Q3. Stage F still requires a committed producer,
fresh output directory and exact-revision verification rather than inferred
acceptance from input counts.

## Contract Decisions Before Coding

### Panel Identity And Scores

- Accept exact, versioned accessions only. Reject aliases, bare factor names,
  duplicate accessions, conflicting display order, missing/invalid full PFMs,
  and case-sensitive factor-name mismatches. Never silently substitute another
  matrix, expand a family, collapse several accessions, or sort by factor name.
- Validate PFM shape, finite nonnegative counts and nonempty columns. Strict
  resolution uses a snapshot with unique exact IDs, not the alias lookup map.
  Carry expected names separately from user-facing labels and preserve available
  taxon metadata; matching capitalization alone is not species verification.
- Preserve all declared panel settings. Unknown or unsupported settings fail
  explicitly instead of being silently ignored. Freeze the JSON structure only
  against Glen's inspected `gentle.jaspar_target_panel.v1` example.
- Recommended Q3 answer: v1 requires one explicit score kind throughout the
  panel. Mixed kinds fail with the conflicting track IDs. Do not silently split
  or convert the panel. Per-track mixed units can be a later additive feature.
- Label the exact score and units. A background-tail `-log10(p)` score is not
  a log-odds score; retain background model, pseudocount, quantization, threshold
  and clipping metadata. Do not change existing score formulas for this task.
- Store computational scores and the explicit display transformation separately.
  Positive-only display must not silently redefine the signal for correlation.
  Existing legacy reports/commands retain their current behavior.
- Snapshot resolution once per run and score the resolved matrices directly.
  Bind actual registry bytes/source, any supplemental PFM source, exact matrix
  counts/digests, panel bytes, and available source URLs. Do not invent a URL for
  a local override or claim that a disk asset identifies different embedded data.

### Bundle Identity And Strand Geometry

- Recommended Q1 answer: use the existing bundle format if Glen's manifest
  already supplies the required fields; add a versioned reader/validator rather
  than creating a competing producer. Genome ID, assembly and annotation release
  are separate identities, not interchangeable strings.
- Required expectations come from a pinned manifest plus explicit caller checks.
  GRCh38, 701 bases, and -500/+200 are properties of the acceptance fixture.
  Generic length is `upstream + 1 + downstream`, with checked arithmetic and
  no silent clipping at contig ends in v1.
- Parse every required FASTA header field, rejecting duplicate keys, malformed
  coordinates, duplicate transcript IDs, duplicate promoter IDs and mismatched
  FASTA/manifest membership. Check exact file digests and normalized sequence
  digests separately, with a documented sequence-byte convention. Accept only
  A/C/G/T/N after declared case/FASTA-whitespace normalization, never inferred
  substitutions or silently discarded bases.
- For a transcript-oriented sequence index `i`, upstream extent `u`, absolute
  TSS `t`, and strand sign `s` (+1 or -1): relative position is `i - u`, and
  genomic position is `t + s * (i - u)`. Use this same transform everywhere.
  Keep genomic intervals stored ascending; separately label decreasing genomic
  axis endpoints on minus-strand figures. Never reverse-complement twice.
- One genomic TSS identity includes reference/assembly, contig, strand and TSS
  coordinate. Preserve associated gene/transcript memberships; never multiply
  it by splice variants or silently lose multi-gene membership.
- Identical sequences at different valid genomic positions are allowed and may
  share cached scores, but never share away their positional identity. Replace
  W2's unconditional rejection of equal sequence digests with consistency checks
  and an explicit repeated-sequence diagnostic.
- Bundle checks establish internal consistency, not independent reference
  authenticity. When a prepared reference is supplied, compare each oriented
  sequence with indexed extraction and record that stronger check separately.
  Hash/file links must remain inside the declared bundle; reject traversal and
  symlink escapes. Bound input sizes and record counts.
- Recommended Q2 answer: an optional, hash-bound selection file cites exact
  promoter IDs plus gene/reference context. Unknown or contradictory references
  fail. Without selection, all records are unselected; never infer selection
  from a gene name or apparent score. Selected-first order has explicit ties.

### Presentation And Comparisons

- The new TSS workflow defaults to independent matrix-specific row scales with
  readable numeric y ticks, explicit units, and visible calibration/non-claims.
  Its shared-scale option requires a compatible typed calibration binding.
  Keep the old renderer's default behavior available for existing consumers.
  Independent per-page ranges must be visible: equal trace heights across TSSs
  do not mean equal scores. A common range for the same matrix across TSS pages
  is a useful later display option, not cross-matrix calibration.
- Map every plotted point through its local window-start coordinate, not its
  vector's fractional length. Show the final `L-1` unscored positions for a
  length-L motif as unscored, not zero. Mark ambiguous-base windows as unavailable
  under the new strict profile contract, not as measured absence.
- Retain both motif-local and genomic strand in machine-readable output. A
  reverse-motif hit's genomic interval and 5-prime endpoint are not necessarily
  the same as the common-axis coordinate used to plot that window start.
  TSV must state this convention and include motif length/interval information.
- W8 compares common valid window-start coordinates. Report metric, input score
  kind, strand treatment, clipping, smoothing policy, sample count and excluded
  positions. Constant/insufficient signals have a typed undefined outcome, not
  a numeric zero correlation. Reuse shared arithmetic with these prechecks.
- Default comparison: unsmoothed Pearson and Spearman, forward-forward and
  reverse-reverse separately. Optional smoothing is explicit and bound. Maxima
  and ranked peaks remain distinct; configurable top-hit count must not be
  silently limited by the existing helper's three-peak default.
- Default digital page: one TSS, full-size readable rows, one SVG per page.
  Thirty rows can require a tall page. If fixed paper sizes are needed, agree
  on continued track pages rather than shrinking text. Multiple TSSs per page
  are optional and only allowed when minimum layout dimensions remain satisfied.

## Delivery Sequence

### A. Freeze Inputs And The Portable Contract (W1, W2, W6/W9 Definitions)

Use the pinned files and Q1-Q3 resolutions above. Document the panel, bundle,
selection, result, comparison and receipt records before building exporters.
Add tiny, clearly synthetic plus/minus fixtures with provenance. Implement exact
panel resolution and bundle validation, including actual source bindings.

Acceptance: malformed/ambiguous inputs fail before scoring or output writes;
order, case and all transcript memberships survive. Valid repeated sequences
at distinct loci are retained. The real panel counts are checked only when the
actual panel is available.

### B. Compute One Shared TSS Profile Report (W6 And Scoring Part Of W1)

Factor a narrow resolved-matrix entry into the existing scoring service, so both
legacy and strict callers share the calculation. Prepare matrix/background data
once per distinct scoring identity per run, not once for every TSS. Build the
portable per-gene reports with exact ordering, validity masks, genomic mappings,
score policies and source digests. Add progress and cancellation between bounded
stages. Stream long-form TSV rather than building it all in memory.

Acceptance: synthetic plus/minus values match the direct shared scorer at the
first valid window, TSS, last valid window and known reverse motif. Geometry
tests separately verify absolute coordinates, and a tiny hand-calculated PFM
checks score values independently of the shared scorer. Cache on/off results agree;
display options do not rerun scoring. Corrupt inputs produce no partial success.

### C. Render The Verified Reports (W3, W4, W5)

Implement per-row numeric ranges and TSS-relative axes together, then pagination.
Use a common tested coordinate/layout model for lines, ticks, logos, overlay
lanes, zero/TSS rules and unscored regions. Derive tick spacing from requested
extents. Compute or measure text bounds; long labels wrap/expand rather than
overrun a fixed gutter. Preserve readable real sequence logos, not consensus
letter substitutes. Paginated headers repeat gene, TSS, page and matrix identity.

Acceptance: -500..+200 and nondefault windows, both orientations, 100-fold
amplitude differences, all-zero/negative tracks, long labels, variable motif
lengths, many TSS markers and 30-row pages. Verify semantic SVG geometry and
rendered examples, not just label presence. Legacy SVG golden tests protect
existing output; any existing coordinate bug must not be copied into the new
mode merely to preserve old bytes.

### D. Add Within-Factor Matrix Comparison (W8)

Group only validated panel accessions sharing the declared factor ID. Consume
the reports from B; do not score again. Add maxima, top peaks, defined/undefined
correlations and a compact section/table to the same document. A single-matrix
factor has zero pairs, not an error. No factor-specific implementation branches.

Acceptance: hand-calculated two-/three-matrix signals, different motif lengths,
ties, constants, missing windows and strand handling. Changing a row's visual
scale or color cannot change comparisons. No matrix is declared the true model.

### E. Expose The Engine Workflow And Audit Exports (W7, W9)

Add one typed computation operation plus a report-consuming export operation;
the convenience shell command composes them. Names remain provisional until
contract review, for example `features tss-tfbs-profiles`. Repeated FASTA inputs,
manifest/panel/selection inputs, expected reference identity and output formats
must have shared parser/CLI/MCP/JS/Lua/Python/GUI-Shell reachability and tests.
Keep a native configuration wizard out of this first delivery, but make the
report usable from the existing GUI inspection/export surface.

No mutation of source sequences is necessary. Filesystem writes still require
the normal effect classification and MCP/agent confirmation. Preflight all
inputs and formats, write into a new staging directory, then publish a complete
output set without overwriting history. Reject nonempty destinations and unsafe
output names; cancellation/failure never leaves a success receipt.

Emit per-gene JSON/TSV/pages, a compact cross-gene index, replay instructions and
`gentle.tss_tfbs_profile_receipt.v1`. Bind all inputs, computation policies,
producer revision/binary, lockfile and every output hash. Exclude the receipt
from its own hash inventory: data files feed the index, and the receipt binds
both; the index must not hash the receipt. Separate
deterministic result identity from run timestamps and host-local execution data.

Prefer existing in-process PNG/PDF conversion. Record backend version/options
and font identities, and call PDFs raster-backed. If an external renderer is
explicitly requested, discover/validate its absolute executable, hash, version,
argv and output streams. Do not add an unconditional rsvg-convert dependency.

Acceptance: two-record/three-matrix end-to-end export, exact output inventory,
tamper detection, ordering/counts, parity and consent/no-writes tests, failed
render/cancel cleanup, and deterministic replays under a pinned font/backend.

### F. Glen's Exact-Revision Replay And Acceptance (W10)

Commit the producer first. Glen runs real replay on that frozen revision in a
fresh directory, retaining historical reports unchanged. Expected fixture-only
counts: CD44 31, TGFB1 6, SERPINE1 4, PATZ1 5 and TP73 12, total 58 TSSs; every
TSS has 30 matrices/28 distinct factors, with the three TFAP2C matrices separate,
adjacent and explicitly identified. Verify transcript-membership counts
109/17/15/13/20 against the supplied manifest and define whether these are
distinct IDs or summed memberships before asserting them.

Review one page per gene plus CD44's first/last pages, including TGFB1/PATZ1
minus-strand geometry, y scales, legends, logos and page legibility. Run focused
Rust/Python tests, locked acceptance/parity/tutorial checks and native inspection
on the same candidate. Compare failures against the exact recorded main baseline.
Record wall time, memory and responsiveness; no unmeasured performance claims.

Full TSV output can approach 2.4 million position/strand rows for this dataset.
Estimate disk use before replay. Keep large replay outputs in retained audit
artifacts; commit only approved public fixtures, small examples and provenance.
Glen can publish through his bot-owned review branch when authorized; this plan
does not authorize pushing, merging or tagging from the development checkout.

## Architecture And Parallel Work

Use `gentle-protocol` for new portable types, `gentle-engine` for small pure
identity/coordinate validators, and `gentle-render` for the new report-driven
layout. Keep filesystem orchestration and shared scoring integration in focused
root engine modules for now. Neither extracted crate may depend on the root or
GUI. Reuse existing scoring primitives; do not start a broad engine extraction.
Large dispatchers receive thin routes only, not entire implementations.

Dependency order is A -> B; C and D can proceed against B's frozen report shape;
E integrates B/C/D; F follows a committed E. Renderer tests can begin with
synthetic report fixtures after A, but W3-W5 are not three independent edits to
the same renderer. Give layout one owner and central wiring one integrator.
Receipt fields belong in A even though receipt writing lands in E.

Replace W1's repository-wide ban on strings such as SERPINE1 with a scoped
check that the new production path contains no fixed gene or matrix list.
Existing public examples, tests and unrelated gene-specific workflows remain.

Keep each delivery a reviewable commit. Update protocol/interface docs and
CHANGELOG when behavior lands, not while these commands are only proposed.
Do not rebuild after `cargo clean` just to review this planning document.

## Scientific Interpretation

Every scientific report, figure, TSV preamble and replay README must retain:

> JASPAR tracks are sequence-model predictions, not measured binding, affinity,
> occupancy, cofactor interaction, promoter activity or functional regulation.
> Cross-matrix magnitudes are not comparable without a documented calibration.

Annotated transcript starts are TSS candidates, not experimentally established
initiation sites. A factor also appearing as a target gene does not establish
autoregulation. Correlated profiles compare model outputs on the same sequence;
they do not show co-regulation or establish that one model is biologically correct.
All descriptive counts and conclusions must be derived from validated output.
