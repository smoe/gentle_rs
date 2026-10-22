# Claude Review Prompt: Base-Editing Protocol Integration

Prepared 2026-10-22 for Steffen to forward. The owner returned Claude's
read-only software-contract review at unchanged `2ea67b1b` on the same date.
See the [review reconciliation and revised S1 contract](base_prime_editing_preparation.md#software-review-reconciliation-2026-10-22).
The prompt below remains the original request, not a current pending review or
implementation authorization. Claude did not independently access the papers
or run tests; no scientific/native acceptance or direct consultation is claimed.

---

Please review this proposed GENtle base-editing integration against the actual
source and existing contracts. This is a **read-only second opinion before
finalizing the plan**, not an implementation request. Do not edit any files,
run builds/tests or external design tools, fetch/rebase/change branches, invoke
another agent, commit, push, dispatch CI, or change tutorial/release status.
Do not inspect credentials, upload sequences or access unrelated user data.

## Repository And Exact Review Scope

Worktree: `/Users/u005069/.codex/worktrees/47dd/gentle_rs`.
Branch: `gentle_rs_2_main`.
Preparation HEAD: `2ea67b1bc35c9a3a662546fa871aad675b8cd169`.
Local `main` matched that SHA during preparation; no fetch was performed.
Do not assume it is the latest remote tip. Report any source drift that you can
verify with your permitted read-only tools; otherwise state the identity limit.

Review the **working-copy documents**, not just their versions at HEAD:

- `docs/base_prime_editing_preparation.md`, especially the proposed typed
  contracts, slice/test tables and **Base Editing Protocol Workstream** section.
- The related parking-lot entry in `docs/roadmap.md`.
- The 2026-10-22 planning entry in `docs/CHANGELOG.md`.
- This review prompt.

These changes are uncommitted planning work, not an implemented capability.
Existing untracked `output/` contains unrelated retained gel artifacts; leave
it alone. Do not switch to the primary checkout and accidentally review a
different revision or omit the working-copy plan.

Read `AGENTS.md` and the relevant sections of `docs/architecture.md`,
`docs/decisions.md`, `docs/biological_extension_guide.md`, `docs/testing.md`,
`docs/gui.md`, `docs/cli.md` and `docs/protocol.md`.

## Owner's Request And Literature

The owner asked Codex to read Porto et al., *Base editing: advances and
therapeutic opportunities* (2020), DOI `10.1038/s41573-020-0084-6`, and prepare
a plan to integrate base editing **as a protocol in GENtle**, not merely export
generic guide oligos.

Publisher: https://www.nature.com/articles/s41573-020-0084-6
Accessible full manuscript: https://pmc.ncbi.nlm.nih.gov/articles/PMC7721651/

Codex read the full manuscript and recorded the reviewed sections/access limits
as L6 in the preparation brief. L7 records the accessible methods alternative,
Osgood et al., *Genome editing with programmable base editors in human cells*
(2025), DOI `10.1016/bs.mie.2025.01.001`:
https://pmc.ncbi.nlm.nih.gov/articles/PMC11948366/

The brief also retains the older Huang/Newby/Liu protocol and its correction,
with honest subscription-access limitations. None of these sources has become
an approved editor profile, default recipe, retained/licensed source payload or
validated predictor in GENtle. The review's 2020 technology landscape must not
be presented as today's complete list or limit. If independent publication
access is unavailable to you, review the recorded requirements/source ledger
and explicitly state that you did not independently inspect the full articles.
Do not infer missing experimental parameters from an abstract.

## Original Codex Plan

Initial scope: conventional CBE/ABE research on nuclear DNA in cultured
mammalian cells, one explicit SNV and one selected reference haplotype. RNA,
mitochondrial, dual/transversion and prime editors, multiplex screens, clinical
and germline use are separate applicability contracts, not forbidden biology.
No additional current-release gate is proposed.

1. **S1: validate reference-bound intent.** Bind exact reference and desired
   alleles, coordinates/orientation, sequence hashes and annotation/reference
   versions. Start with synthetic tests; refuse stale/mismatched or unresolved
   reference/haplotype input. No implicit editor selection or source mutation.
2. **S0/S2: evaluate editor/guide combinations.** Use explicitly synthetic
   profiles until exact CBE/ABE constructs are reviewed. Validate both-strand
   spacer/PAM/window geometry, conversion compatibility, potential bystanders
   and available per-transcript consequences. Report bounded-search completeness
   and missing evidence, not efficacy/safety or invented outcome probabilities.
3. **S4a: compose a read-only protocol draft.** Bind one selected candidate to
   context/SOP evidence, exact components, ordered stages, controls, readouts,
   unresolved fields and scoped review records. Retain absent parameters;
   invent no quantities, timing, vectors or acceptance thresholds. Invalidate
   affected approval when a bound input changes. An incomplete draft is still
   inspectable and is not an executable laboratory SOP.
4. **S4b: inspect and hand off.** DNA-viewer controls display the same saved
   engine result as CLI/agents: allele, spacer/PAM, editable strand, window,
   bystanders and transcript consequences. Export JSON/Markdown/PDF without
   recalculation. Reuse existing jobs/rendering/document helpers. Later exact
   component materialization requires separate approval; no automatic ordering,
   delivery, source-DNA edit or default lenti-vector choice.
5. **S5: attach observed results.** Import local allele-count/CRISPResso2 output
   with sample/control, amplicon/reference identity, counts and denominators.
   Keep intended-only, bystander and other outcomes separate. Bulk amplicons
   do not establish clone genotype, phenotype or genome-wide specificity.
   Guide-dependent DNA, guide-independent DNA and RNA evidence stay distinct.
6. **S6: activate a tutorial only when real capabilities exist.** Extend the
   existing inactive draft with synthetic replay and then one separately
   approved public example. Keep native screenshots, computational checks,
   experimental evidence and human scientific interpretation separate.

First useful milestone: S1 plus synthetic S2 and a saved draft preview.
Recommended first implementation authorization: **S1 only**, not the entire
workbench. The existing prime-editing S3 path is not a dependency.

## Evidence Already Gathered And Source Paths

- `src/engine/protocol.rs`: `GuideCandidate`, `GuideSet`, practical-filter and
  oligo records. `src/engine.rs`: `normalize_guide_candidate` validates fields,
  but does not establish editor-specific sequence/PAM/window matching.
- `src/engine/ops/candidate_guides.rs`: guide filtering, oligo generation and
  `op_export_guide_protocol_text`. The export is a cloning checklist with
  storage/filesystem effects, not a read-only base-editing protocol evaluator.
  The generic cut offset and built-in oligo templates are not editor profiles.
- `src/engine/analysis/variant_promoter.rs`: reference-bound SNV allele checks;
  reuse appropriate validation without silently widening its mutation contract.
- `crates/gentle-protocol/src/biological_context.rs` and `genomic_regions.rs`:
  reusable context/region identities, not completed edit-intent validation.
- `src/engine/analysis/sequence_design.rs` and the corresponding protocol:
  provenance/proposal/apply precedent, not an editor-chemistry solver.
- `src/engine/state/sequence_ops.rs`, `src/engine/state/lab_assistant_export.rs`
  and `LabAssistantInstructionsExport` in `src/engine/protocol.rs`: reusable
  material/checkpoint/document presentation. `ExportLabAssistantInstructions`
  v2 takes journal/state, not an arbitrary experimental protocol draft.
- `src/engine_shell.rs`, `src/engine_shell/command_parsers.rs`, `src/ui_intent.rs`
  and current DNA-viewer/shared-availability paths: shared routing, approval and
  thin-GUI integration seams. Check actual contracts rather than assuming all
  existing approval models are interchangeable.
- `docs/rna_guides_spec.md`, `docs/sequencing_confirmation_plan.md`,
  `docs/gui_cli_mcp_parity.md` and `docs/agent_interface.md`: present capability
  limits and discovery/reachability policy.
- `docs/tutorial_drafts/base_prime_editing.md` and
  `scripts/test_tutorial_gap_lessons.py`: inactive synthetic teaching material.
  Twenty existing tests passed on the preparation working copy; these are
  teaching/source checks, not evidence that an editing engine exists.

No Rust build, new CLI replay, native GUI experiment or biological validation
was performed for this planning change. No model/tool version was selected.
Proposed type/operation names in the brief are **not callable commands**.

## Questions To Challenge

1. Is S1 the correct smallest useful primitive, or can an existing validator
   satisfy it without another model? Can the first SNV slice stay narrow while
   retaining a clean path to later indels? Identify real missing inputs.
2. Are genomic-forward alleles, local display orientation, guide direction,
   hybridized strand and chemically edited strand kept unambiguous? Check CBE
   C-to-T/G-to-A and ABE A-to-G/T-to-C projections and nonzero anchors.
3. Does the proposed window model distinguish supported activity, weaker/context
   evidence and unknown activity outside the window, without either claiming
   zero activity or pretending all positions are equivalent? Does a preference
   accidentally become a hard incompatibility?
4. Is the proposed protocol draft genuinely an editing-specific composition,
   or does it duplicate existing workflow, evidence, review/approval, document
   or process records? Which exact existing helpers should it reuse? Challenge
   any new generic registry/framework or silent v1/v2 semantic extension.
5. Is requiring context/SOP selection consistent with an inspectable incomplete
   draft and the synthetic first milestone? Define the minimum draft inputs
   versus the gates needed for a human-reviewed laboratory handoff.
6. Do retained unknowns, candidate/SOP changes, save/reopen and report retention
   have clear effect/approval semantics? Are read-only evaluation, file export,
   project persistence and component creation adequately separated?
7. Does the plan falsely imply transcript consequences, a PAM-aware search or
   an amplicon prove biological function, clone genotype, guide-independent DNA
   specificity or RNA specificity? Are the result counts/denominators adequate?
8. Is CLI/MCP/scripting/inner-agent parity concrete without inventing commands
   or a second biological implementation? Can the DNA viewer remain a thin
   report inspector using current availability and stale-job rules?
9. Which independent literal tests are indispensable before the first merge?
   Include both orientations, edge/window cases, source mismatch, synthetic
   identity, cancelled/truncated search, no-mutation and legacy compatibility.
10. Which work should be removed or deferred? Give one minimal, reviewable
    implementation request rather than turning every later slice into a gate.

## Requested Reply

Return a concise verdict, then at most six actionable findings ordered by risk,
each with a current file/symbol reference and smallest correction. Distinguish
verified defects in the plan from optional simplifications and decisions that
need Steffen or a laboratory reviewer. Say explicitly where source/literature
access prevents verification.

Then give a revised minimal implementation order and the exact first-slice
request: required inputs, shared seam, outputs/effects, fail-closed conditions,
independent tests and explicit non-goals. State what is already sound and what
can wait. Do not grant release, native GUI, efficacy or scientific acceptance.
Make no repository changes; Steffen will return your reply to Codex for an
explicit original-plan -> feedback -> revised-plan reconciliation.
