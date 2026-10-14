# Tutorial Coverage and Screenshot Backlog

Updated 2026-10-14 at source `024ed6c0520d38c4221517219d44c792084de4e5`.
This is the current working list; `coverage_audit.md` remains a historical
snapshot. Catalog presence, executable replay, native interaction, screenshot
quality and human scientific review are separate states.

## Screenshot queue

Use only synthetic or public inputs. Retain raw captures and provenance
sidecars; label crops and annotations as derived teaching images.

1. **08.16 imported TSS evidence:** complete and partial lanes, strand
   inversion, clipped footprints, signed/zero scores, coverage warnings and
   readable normal-size presentation.
2. **07.04 microarray evidence:** native locus navigation, readable evidence
   configuration/results and GUI save/reopen.
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

Do not reopen the completed 08.18 repeat-cohort run, tutorial-CI registration
or release-link repair merely to obtain more images.

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

