# Synthetic Cloning Planning Defaults

`cloning_wetlab_steps.json` is hand-authored, publication-safe planning data
created for the time/cost separation slice on 2026-10-19. No values were
measured in a laboratory or derived from a published protocol. Each step's
`source.kind` and `source.note` preserve this provenance.

Only `restriction.digest_ligate_extract_sticky` has a model. These are separate
bench-planning allowances, not durations attached to `Digest`, `Ligation` or
`ExtractRegion`. The computational pattern script remains unchanged. The
modeled endpoint is recovery/review, not a verified colony or final construct;
transformation, screening, sequencing and experimental success are not modeled.

## Reproduction And Use

The model is embedded in GENtle, validated against the bundled routine IDs,
and used by the shared planning estimator for `planning consult cloning`,
`routines list`, `routines explain` and `routines compare`. It does not depend on
a network lookup or a local vendor. Recreate the model by preserving the
literal ranges and dependency list in this file; engine-shell planning tests
exercise the same data and hand-crafted malformed variants.

The dependency scheduler orders eligible step IDs lexicographically. Within
one step, hands-on and unattended durations are sequential. Independent steps
may overlap: preparation of later ligation supplies does not depend on digest
incubation. The three critical paths are calculated independently using each
step's minimum, typical or maximum duration; the ranges are scenario bounds,
not confidence intervals or probabilities. No staffing or shared-instrument
capacity is scheduled. A critical path is therefore an idealized planning
bound, not a promised completion time.

Whole-model hands-on and unattended ranges sum their respective work across
all steps; they are not both elapsed critical-path time. Existing summed
machine queue and material procurement delays are an explicitly disclosed
prefix before bench steps. The legacy business-day conversion and lack of
holiday/date anchoring remain unchanged. For the other routines, the three
new ranges are unknown and only the labelled legacy scalar is available.

Scoring and the `max_time_hours` guardrail retain the legacy family/base scalar
plus summed delays for every routine. Adding/removing this model, or changing
its synthetic durations, must not alter those policy decisions. The shared
estimate discloses `explanation.time_policy`; these uncalibrated proxies are
not cross-method laboratory-time evidence. Old reports without recorded
policy metadata stay unknown, not reinterpreted through the current model.

Material figures are one priced unit per inferred material class: existing
stock versus a missing-material purchasing proxy. They are neither a complete
bill of materials nor an orderable quote. Package sizes, quantities, expiry,
minimum orders, shipping and tax are not modeled. Unknown/unpriced classes
remain listed rather than becoming known zero-cost materials. Machine charges
retain the old family-duration basis, not the new critical-path duration.

The record supports transparent planning arithmetic only. It supplies no
enzyme-specific conditions, protocol validation, biological feasibility,
clinical validation, yield prediction, outsourcing permission or order approval.
