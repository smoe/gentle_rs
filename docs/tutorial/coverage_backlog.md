# Tutorial Coverage and Screenshot Backlog

Updated 2026-11-04 after reviewing Glen's microarray evidence at `37bc0cea`.
This is the current working list; `coverage_audit.md` remains a historical
snapshot. Catalog presence, executable replay, native interaction, screenshot
quality and human scientific review are separate states.

## Screenshot queue

Use only synthetic or public inputs. Retain raw captures and provenance
sidecars; label crops and annotations as derived teaching images.

1. **08.16 imported TSS evidence:** complete and partial lanes, strand
   inversion, clipped footprints, signed/zero scores, coverage warnings and
   readable normal-size presentation.
2. **07.04 microarray evidence:** the two retained Linux captures now cover a
   reopened project and its navigation menu. Opaque/high-contrast teaching
   images and cross-platform GUI/save-reopen acceptance remain pending.
3. **10.03 measured gel sizing:** native lane navigation, calibrated and null
   out-of-calibration results, and GUI save/reopen.
4. **FLNA draft:** normal-size native overview/detail and orientation. The
   2400-px render and human 24-bp teaching-exon choice are already retained;
   historical annotation-release provenance remains missing.
5. **04.09 primer specificity:** native Shell/report presentation. Defer the
   full real-reference import until the measured post-BLAST performance issue
   has a verified fix.
6. **06.07 synthetic sequence design:** preview, exact proposed output,
   created DNA and same-session undo/redo.
7. **03.05 and 10.02 handoff points:** review-before-apply and sequence-
   confirmation save/reopen, for the new Synthetic Biology learning route.

Retained 07.04 evidence (Glen / `smoe-bot`, `37bc0cea`): two checksum-bound
Linux/Xvfb observations show both synthetic-statistics contrast groups, one
selected vendor-derived feature and the nearest-feature menu. The recorded GUI
producer is `503cb391`; the CLI report producer is `3839426f`. These are not
current-main or cross-platform acceptance. Genuine CEL replay and independent
domain interpretation remain pending. CD44/TSS parent-branch captures are not
imported by this microarray-only port, so their queue entries are not closed.

Do not reopen the completed 08.18 repeat-cohort run, tutorial-CI registration
or release-link repair merely to obtain more images.

### Overnight Review Follow-Up, 2026-11-04

The forwarded review covers Glen's local `1e75d926` on top of the public
`069bd37b` boundary. The newer overview and captures are now published on
Glen's branch ending at `37bc0cea`. Only that tip's microarray tranche is
ported here; the parent overview/CD44/TSS changes need separate review and
semantic repair. Do not treat availability as an instruction to merge blindly.

- For 06.08, current main's tables and derived SVG are now explained separately.
  Verify unavailable derived lengths versus external expectations natively.
  Before integrating the published overview, reuse the shared renderer's exon-family
  identity for both exon and CDS/protein segments; test a CDS starting within
  an exon, reverse-strand ordering and the historical 742/139-aa checkpoints.
- Keep raw captures and manifests byte-identical. For new teaching captures,
  use an opaque, high-contrast profile and readable text at ordinary tutorial
  width. A crop cannot remove translucency; recapture if background detail
  obscures the labels. Crop/annotate only a separate derivative, never evidence.
- Each derivative needs its raw input hash/path, crop bounds, transformations,
  output hash and source revision. Label it as a teaching derivative and retain
  a link to the raw capture. Do not call an SVG export a native screenshot.
- Keep historical machine-local capture locators intact, labelled as
  non-replayable context when their manifests arrive. Add repository-relative
  or logical replay inputs separately; a `/home/...` or ephemeral `/tmp/...`
  path is not a public download or a reproducible input identity.

## Under-taught product areas

Address these as separate reviewable commits rather than one large tutorial
rewrite.

1. **Synthetic Biology integration — addressed first.** The new
   [Design, Build, Track and Verify route](./03-09_synthetic_biology_design_build_verify.md)
   connects existing unit tutorials while making their fixture and evidence
   boundaries explicit.
2. **Assembly-method selection.** Compare Gibson with shipped Golden Gate,
   Gateway, TOPO, TA/GC, In-Fusion and NEBuilder planning routes. Start with a
   bounded synthetic Golden Gate example; do not imply every preset has equal
   native or experimental acceptance.
3. **RNA structure.** The ssRNA structure interface needs a beginner exercise
   with a small public/synthetic sequence, parameter visibility and a clear
   distinction between prediction and measured structure.
4. **microRNA target scanning.** Add a compact walkthrough for target scan,
   coordinate interpretation and uncertainty; avoid expression or repression
   claims without external evidence.
5. **Guide design end to end.** Existing generated guide-filter and oligo-
   export chapters need a guided route from target/constraints through review
   and export. Do not present the current checks as a complete off-target or
   efficacy assessment.
6. **General restriction/ligation and physical sample tracking.** Extend the
   minimal digest/ligation replay and Gibson-specific arrangement/rack lessons
   into a method-neutral beginner route only after the higher-priority gaps.

Strong existing coverage in sequence inspection, Gibson planning, PCR/qPCR,
genome/transcript context, regulatory evidence, external handoffs and Sanger
confirmation should be linked rather than duplicated.
