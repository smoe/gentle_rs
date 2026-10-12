# Poster Figure: Isoform-Aware Primer-Pair Design

Prepared 2026-09-20 against GENtle commit
`e1c7dfb221914bdfe18ba14248cf2db625503552`.
Revision 2 incorporates independently checked points from the user-supplied
Claude review. No additional Claude invocation was made.

## Files And Print Use

The committed master is `docs/figures/isoform_primer_design_inputs.svg`.
PDF/PNG companions are optional local exports, not committed analysis results.

- `isoform_primer_design_inputs.svg`: editable native vector artwork, 600 x
  450 mm (4:3 landscape). Text remains editable. Uses Avenir Next and Georgia,
  with fallbacks; another machine may substitute fonts in the SVG.
- `isoform_primer_design_inputs.pdf`: print companion with embedded fonts;
  use this version to preserve the typography when placing the figure.
- `isoform_primer_design_inputs.png`: screen preview, not the preferred print
  master. Use the SVG or PDF for scaling on a poster.

The figure can be reduced proportionally. At its nominal 600 mm width the
smallest labels are about 15.6 pt. At 400 mm width they become about 10.4 pt;
use the larger width for comfortable poster reading. Colors encode role, not
scientific confidence. Exon labels and hatching supplement color distinctions.

## Suggested Caption

**Gene-informed design of isoform-aware primer pairs in GENtle.** Explicit
reference and transcript models, assay intent, primer/product constraints and
panel-selection policy define the design problem. Optional source-bound
biological evidence informs study priorities and declared junction obligations
without replacing transcript structure. Mature-cDNA models support candidate
pair design and transcript-by-assay product assessment; identical mature cDNAs
belong to the same sequence class. No-product, single-product, multiple-product
and unassessed outcomes must remain distinct. Optional specificity gating can
run during design; an external plan/execute/finalize route also supports explicit
redesign after a bound specificity failure. Nine readiness gates cover
reference/product evidence, oligo safety/identity and experimental
interpretability. Which gates are required is explicit in the policy; gel
assessment uses named run assumptions and does not establish observed band
separation. Full-oligo thermodynamic follow-up remains separate from internal
complementary-run QC. Reports retain sequences, predicted products, selection
reasons, pair-score terms and provenance. Candidate generation, assessment,
policy-qualified ordering readiness and experimental validation are distinct.
Transcript models, matrix entries and co-migration sketches are conceptual
illustrations, not measured or computed assay results.

## Interpretation Boundaries

This is a workflow-level input map, not a claim that every input is available,
automatically retrieved, used in one numerical score, or required by every
operation. The study planner, pair designer, panel selector, locus evidence
viewer and QA routes are separate parts of the shared-engine workflow.

1. **Reference scope:** Ensembl and RefSeq can be compared, but displaying an
   additional source does not add its records to the approved design targets.
   CDS/UTR and start/end structure are annotation context; this figure does not
   claim the deferred single-primer capture/CDS-retention programme is complete.
2. **Evidence:** RNA-read, cDNA/EST, array, expression and previous qPCR records
   are optional, source-bound dimensions. Study priority, contrast thresholds,
   required/preferred junction obligations and overrides are explicit. A
   projected observation is not a new independent measurement. Occupancy and
   chromatin tracks are review context, not a universal primer-ranking term or
   proof of isoform-specific regulation. No automatic tissue-specific expression,
   protein-domain optimization, causal inference or statistical significance is
   claimed.
3. **Sequence orientation:** every mature-cDNA cartoon runs 5' to 3' from left
   to right, irrespective of the gene's genomic strand. A forward annealing
   segment is a subsequence of that cDNA; a reverse annealing segment is the
   reverse complement of its target region. Added 5' tails need not match the
   target. The arrows on class A point inward and the forward primer spans E1/E2.
4. **Illustrative matrix:** A = E1-E2-E3, B = E1-E3 and C = E0-E2-E3.
   Hypothetical assays in shared E3, across E1/E2 and across E1/E3 illustrate
   possible detection patterns. `0`, `1` and `M` represent the engine's
   `NoProduct`, `SingleProduct` and `MultipleProducts`. `?` is a display legend
   for unassessed evidence, not a fourth `TranscriptAssayDetectionStatus` variant.
   The B/E1-E3 cell deliberately shows `M`: one transcript can yield several
   products, for example if the other primer has multiple usable binding sites.
   Those sites are not depicted or computed; exon names alone do not establish
   product multiplicity. This assumes suitable binding sequences;
   no bases, product lengths, optimization result or real primers were generated.
   A common assay covers all classes but cannot discriminate them by itself.
   The cartoon omits actual sequences and excludes any assertion of measured
   abundance. In real outputs, unassessed and predicted no-product are distinct.
5. **Specificity and variants:** the search references, exact versions,
   completeness, intended-product definitions, mismatch policies and variant
   source must be explicit. Missing allele frequency remains unknown. Any
   changed physical oligo requires new assessment rather than inherited approval.
   `--specificity-check require-pass` with a prepared target can reject a
   `DesignTranscriptAssayPanel` call before panel persistence; this is the inline
   local-BLAST route, not a universal automatic check or proof that both genomic
   and transcriptome readiness gates were satisfied. The external route is
   plan -> execute -> finalize. Only `pass` is accepted. `specificity_fail`
   means complete execution evidence with biological-policy failure;
   `not_assessed` retains unavailable intended-target geometry; `incomplete`
   covers missing/failed execution, stale state or invalid provenance. Only a
   bound `specificity_fail` acceptance feeds the redesign route, with the exact
   original operation. Passing assays are retained, failed footprints excluded,
   and replacements remain candidates requiring assessment and approval. The
   return arrow depicts this conditional loop, not automatic redesign on every
   missing result. A short/top-hit list is still not exhaustive specificity.
6. **Thermodynamics:** GENtle's exact complementary-run QC is not a free-energy
   model. Hairpin/self-/cross-dimer thermodynamic follow-up is deliberately
   separated from this internal screen. It must assess full synthesized oligos,
   including adapters, and only the oligo combinations that coexist in the
   relevant reaction. Retain the actual backend/version and its temperature,
   salt/ion and oligo-concentration assumptions; do not invent conditions or
   scores. This illustration does not claim an automatic integrated
   thermodynamic gate for every transcript-panel path.
7. **Handoff:** ordering readiness depends on the named policy and bound
   evidence, not a design score, the historical `order_ready_primers` field or
   a vendor-form export. The poster's broad middle label, **Assessed**, includes
   the concrete reachable `specificity_checked` state; it does not claim that
   every gate has passed. The card constructor can assign `candidate`,
   `specificity_checked` or `order_ready`; `wet_lab_validated` is declared but
   not assigned by current production code. The arrow chain is a communication
   model, not a requirement to visit every state in order. Lab efficiency,
   specificity and biological interpretation still require experimental review.
   This figure is about primer pairs for transcript assays, not single-primer
   Nanopore enrichment.
8. **Interpretability and tiers:** common-region screening, isoform
   discrimination and long-range structure discovery are declared use tiers,
   distinct from assay chemistry and the panel-selection objective. The
   informative selector's `minimum_resolvable_size_difference_bp` is an
   explicit product-size-separation policy, not a physical gel calibration.
   Separately, `ExperimentalAssayGelAssessment` uses named `GelRunConditions`
   to project whether distinct predicted product lengths co-migrate. Without
   conditions it is not evaluated; missing usable products or failed evaluation
   is incomplete. Co-migration fails a required gel gate. The cartoon illustrates
   a possible merge of two lengths, not a computed or observed result. Passing
   a virtual gel does not prove transcript identity or quantitative abundance.
9. **Backend and explanations:** `auto`, `internal` and `primer3` are available
   design-backend choices; retain the actually resolved backend and version.
   `PrimerDesignScoreTerm` exposes additive `raw_value * weight` contributions
   and the baseline behind an existing pair score. These pair-score terms are
   not an alternative explanation of all panel selection, a new biological
   confidence score or thermodynamic validation. Evidence obligations and
   whole-panel coverage decisions retain their own reasons.

## Nine Readiness Gates

The bold lines in panel 06 name all nine gates emitted by the handoff builder.
Not all are mandatory under the default policy. A required failed, incomplete
or not-evaluated gate blocks readiness; explicitly waived/not-applicable gates
must not be presented as completed assessments.

| Gate | Poster group | Requirement at this revision |
| --- | --- | --- |
| `annotation_provenance` | Reference & product checks | Required by default; recorded annotation release, not an independent authenticity audit |
| `cdna_assay_test` | Reference & product checks | Policy-controlled, not required by default; verifies bound pair identity and product evidence |
| `genomic_carryover` | Reference & product checks | Required by default through the specificity master switch and this dimension's flag |
| `transcriptome_specificity` | Reference & product checks | Required by default through the specificity master switch and this dimension's flag |
| `critical_qc` | Oligo safety & identity | Required by default; stored critical QC, not a universal thermodynamic screen |
| `variant_evidence` | Oligo safety & identity | Policy-controlled, not required by default; provenance-bound variant report |
| `duplicate_review` | Oligo safety & identity | Required by default for linked order-form duplicates; not applicable with no linked order form |
| `experimental_practicality` | Experimental interpretability | Always evaluated as required; a long-range fallback fails if policy disallows it |
| `gel_resolution` | Experimental interpretability | Policy-controlled, not required by default; requires named run assumptions to evaluate |

The existing default allows long-range fallback, but that remains a visible
policy choice rather than permission inferred from a successful candidate.

## External-Pair Entry: Proposed, Not Implemented

The entry routes are real: `primers import-external-pairs` accepts external
sequence pairs with retained provider/source provenance, while
`primers primerbank test-cdna` explicitly tests catalogue pairs against the
selected transcript model. Both reuse shared assessment paths; catalogue or
supplier claims are not GENtle validation. They do not imply that every import
automatically becomes a finalized transcript-panel readiness card.

The poster grid remains 600 x 450 mm. These placements are proposed only:

1. **Recommended if extra height is available:** add a slim alternative-entry
   strip between the upper grid and assessment band, labelled "Existing pairs:
   PrimerBank / literature / supplier / laboratory". Connect it to cDNA testing
   and specificity, bypassing de novo generation, with source provenance
   retained. Allow approximately 25 mm extra height rather than shrinking text.
2. **Keep current dimensions:** use a separate small companion inset beside the
   main figure, showing exact F/R sequences + provider provenance -> shared
   cDNA / QC / specificity assessment. This avoids making the already-full
   reference/evidence cards imply that imported pairs are transcript annotations.

Neither entry was drawn into the figure pending the layout decision. The
backend label and additive-score wording did fit without reducing font sizes.

## Repository Source Map

All claims are grounded in the local repository at the revision above, not in
new biological measurements or an external reference search.

| Figure element | Primary repository evidence |
| --- | --- |
| Reference scope, real transcript models, source comparison | `docs/tutorial/04-08_gene_assay_study_gui.md`, sections 1-3; `docs/transcript_source_presentation.md` |
| Intent, orientation, source geometry, evidence separation | `docs/ai_task_playbooks.md`, Playbook 12; `docs/protocol.md`, cDNA PCR/qPCR test contract |
| Optional evidence types and separate specificity/abundance/responsiveness/assayability dimensions | `crates/gentle-protocol/src/isoform_evidence.rs`: `GeneIsoformEvidenceRequest`, `IsoformEvidenceSourceKind`, `GeneIsoformEvidenceComponents` |
| Study priors, contrasts, observations, threshold policy and exact-operation planning | `src/engine/protocol.rs`: `GeneIsoformAssayStudyPlanRequest`, `GeneIsoformAssayStudyPolicy`, `GeneIsoformAssayStudyEvidenceSummary` |
| Primer and pair constraints | `src/engine/protocol.rs`: `PrimerDesignSideConstraint`, `PrimerDesignPairConstraint`; `docs/architecture.md`, Primer/qPCR command contract |
| Objective, coverage, budget, exact-cDNA classes and RT reach | `src/engine/protocol.rs`: `TranscriptAssayPanelObjective`, `TranscriptAssayInformativeSelectionPolicy`, `TranscriptAssayCoveragePolicy`; `docs/architecture.md`, multi-transcript assay contract |
| Product multiplicity, use tier and minimum observable size gap | `src/engine/protocol.rs`: `TranscriptAssayDetectionStatus` at line 13723, `TranscriptAssayUseTier` at line 7392, `TranscriptAssayInformativeSelectionPolicy` at line 12295 |
| Differential junction obligations and provenance | `docs/tutorial/04-07_transcript_assay_followup_gui_cli.md`, sections 3-4 |
| Specificity and full-oligo QA separation | `docs/ai_task_playbooks.md`, Playbook 12; `docs/protocol.md`, `gentle.oligo_qc_report.v1` |
| Assembly-bound variant screening | `docs/cli.md`, `primers screen-variants`; `src/engine/protocol.rs`: `PrimerVariantScreenReport` |
| Approval, ordering and provenance | `docs/decisions.md`, gene-isoform study approval contract; `docs/gene_assay_study_gui_plan.md` |
| Nine gates and reachable readiness states | `src/engine/ops/operation_handlers.rs`: handoff gates at lines 32856-33061 and state assignment at lines 33083-33093; `src/engine/protocol.rs`: `ExperimentalAssayReadinessPolicy` at line 14586 |
| Optional specificity gate during design | `src/engine_shell/command_parsers.rs`: panel usage and parser at lines 7378 and 7617; `src/engine/ops/operation_handlers.rs`: inline assessment and rejection at lines 31642 and 31696 |
| External specificity finalization and conditional redesign | `docs/protocol.md`: transcript-assay specificity contract at lines 10000-10047; `docs/cli.md`: plan/finalize/redesign routes at lines 3345-3347 |
| Gel assumptions and co-migration | `src/engine/protocol.rs`: `ExperimentalAssayGelAssessment` at line 15123; `src/engine/ops/operation_handlers.rs`: `experimental_assay_gel_assessment` at line 32156; `crates/gentle-protocol/src/lib.rs`: `GelRunConditions` at line 2279 |
| Existing-pair route, proposed inset only | `docs/cli.md`: PrimerBank at lines 3349-3368 and `primers import-external-pairs` at lines 3369-3387 |
| Backend and additive pair-score explanations | `docs/cli.md`: panel backend at line 3421; `src/engine/protocol.rs`: `PrimerDesignScoreTerm` at line 6799 |

## Re-export

From the repository root, using librsvg (available on the authoring machine):

```bash
mkdir -p output/pdf
rsvg-convert --format pdf --output output/pdf/isoform_primer_design_inputs.pdf docs/figures/isoform_primer_design_inputs.svg
rsvg-convert --width 2400 --output output/pdf/isoform_primer_design_inputs.png docs/figures/isoform_primer_design_inputs.svg
```

No application code, scientific input files, generated tutorials, release
metadata or benchmark data were changed to create this figure.

## Revision 2 Verification

- Parsed the SVG as XML and re-exported both PDF and PNG with librsvg.
- Rendered the final PDF through Poppler and visually inspected the complete
  page for placement, legibility, primer orientation and the conditional loop.
- Checked PDF text geometry against the page and card margins; no overflow.
- Confirmed all nine gate labels in the exported PDF's text.
- Confirmed one 600 x 450 mm page, embedded fonts and zero raster images.
- Rechecked the interpretation against source types, policy defaults and
  production gate/finalization code. No biological computation, application
  rebuild, benchmark, or experimental validation was performed.
