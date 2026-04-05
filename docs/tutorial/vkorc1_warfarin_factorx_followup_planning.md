# VKORC1 / Warfarin Factor X Follow-up (Planned Companion Tutorial)

> Type: `planning note`
> Status: `planned`
> Drift note: this page is a roadmap-grounded tutorial concept, not yet a
> validated executable walkthrough.

This note records the intended second ClawBio-motivated warfarin tutorial that
should follow, not replace, the main
[`VKORC1 / rs9923231` mammalian reporter tutorial](./vkorc1_warfarin_promoter_luciferase_gui.md).

## Purpose

The first tutorial stops at a biologically sensible handoff:

- retrieve the `VKORC1` / `rs9923231` locus,
- choose one TSS-centered promoter fragment,
- place that fragment upstream of luciferase in a mammalian reporter backbone,
- export a reproducible construct-design project.

The planned second tutorial would ask a different question:

- after the promoter reporter is defined, how should we extend the construct
  family to investigate a more downstream warfarin-relevant target effect?

## Working Construct Idea

Current working concept from the ClawBio discussion:

1. keep the first `VKORC1` promoter->luciferase construct as the regulatory
   baseline,
2. extend that baseline with a second construct layer,
3. interrupt that second construct with one Factor X target site, and
4. use the resulting design to plan a later functional follow-up assay.

## Why This Stays Separate

This should remain a separate tutorial because it changes the main scientific
question:

- tutorial 1 asks whether the `VKORC1` promoter fragment containing
  `rs9923231` is worth taking into a human-cell regulatory reporter assay,
- tutorial 2 would ask how to probe a more downstream warfarin-relevant target
  mechanism once the promoter-side construct is already fixed.

Mixing both stories into one tutorial would blur the handoff point and make the
first story less reproducible.

## Open Design Questions

- what exact second cassette should be used,
- how the Factor X target site should be encoded and annotated in GENtle,
- what control constructs are needed,
- whether the readout remains luciferase-only or becomes a compound reporter
  design,
- what the minimal, bench-facing success criteria should be.

## Planned Output

When implemented, the companion tutorial should still end at a reproducible
design/handoff point:

- one construct family definition,
- one stable workflow or command sequence,
- one exported construct map,
- one short bench-facing next-actions list.

