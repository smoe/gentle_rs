---
name: gentle-pcr-primer-design
description: >-
  A focused OpenClaw/ClawBio conversation and routing layer for deterministic
  PCR, RT-PCR, SYBR-qPCR, and TaqMan assay preparation. It delegates every
  scientific operation to GENtle through the registered gentle-cloning skill,
  preserving GENtle's typed reports, evidence states, provenance, and review
  gates.
version: 0.7.0
author: GENtle project
license: MIT
tags: [primer-design, pcr, rt-pcr, sybr, qpcr, taqman, isoforms, specificity, provenance]
metadata:
  openclaw:
    requires:
      bins:
        - python3
      env: []
      config: []
    always: false
    homepage: https://github.com/smoe/gentle_rs
    os: [macos, linux]
    install: []
    trigger_keywords:
      - gentle primers
      - gentle pcr
      - primer design
      - pcr primer design
      - primer preflight
      - primer3 preflight
      - conventional pcr design
      - endpoint rt-pcr panel
      - isoform pcr panel
      - transcript assay panel
      - multi-gene isoform batch approval
      - approved isoform workflow batch
      - oligo-dt cdna primer design
      - cap-dependent 5-prime race
      - rlm-race primer design
      - full-length cdna validation
      - mrna cap 5-prime assay
      - sybr junction assay
      - qpcr design
      - taqman design
      - imported primer review
      - primer specificity
      - genomic primer specificity
      - transcriptome primer specificity
      - collection primer specificity
      - compose gene assay routine
      - primer report
      - qpcr report
      - pcr protocol cartoon
      - cdna product gel
---

# GENtle PCR Primer Design

This skill is a specialized conversational surface over GENtle's existing
primer operations. It has no independent primer-design algorithm, specificity
rule, evidence interpreter, order-readiness rule, or command runner.
`INTENTS.json` delegates every executable step to the registered
`gentle-cloning` skill, which in turn invokes GENtle's parser-validated CLI or
workflow engine.

## Boundary

GENtle owns:

- sequence and transcript resolution;
- Primer3 and internal-backend invocation;
- primer, probe, product, and junction calculations;
- genomic and transcriptome/cDNA specificity evidence;
- pass, fail, not-assessed, incomplete, and review-required states;
- persisted reports, provenance, order forms, and deterministic artifacts.

This skill owns:

- recognizing a primer-design request;
- asking only for assumptions that materially change the assay;
- selecting and sequencing existing GENtle operations;
- reporting progress and missing prerequisites;
- explaining GENtle results without upgrading their evidence status;
- stopping at explicit human-review gates.

## Procedure

1. Clarify the biological target, organism/build/annotation, template kind,
   assay kind, inclusion or discrimination objective, and required specificity
   spaces. Ask concise questions only when an omitted answer changes the
   operation.
2. For transcript assays, also establish transcript or junction scope,
   cDNA-synthesis assumptions, product-size policy, and whether endpoint
   RT-PCR, primer-only SYBR, or primer-plus-probe TaqMan is intended.
   Record reverse-transcription priming separately from any cap-dependent
   5-prime capture or enrichment method.
3. Run the GENtle primer preflight before an expensive design. Report the exact
   missing executable, resource, or input. Never install it automatically.
4. Treat natural-language routing as a draft. For every route that selects
   biological material, a pair rank, or a backend, mutates state, or writes an
   artifact, ask `gentle-cloning` for an execution proposal. Show its exact
   operation, assumptions, resolved paths, selected material/reference spaces,
   pinned backend/runtime, hashes, and proposal digest. Broad phrases such as
   "primer design" never authorize immediate execution.
   For a local Cargo launcher, review the hash-bound built `gentle_cli` path;
   for OCI, require an immutable image digest rather than a mutable tag.
5. After the caller attests approval of that exact digest, execute only the
   normalized request stored in the proposal. Do not reconstruct slots or the
   command from prose a second time. List/show/preflight routes may remain
   automatic.
6. Dispatch the approved request through `gentle-cloning`:
   - `primers design` for conventional PCR;
   - `primers design-qpcr` for generic TaqMan-style assays;
   - `primers design-transcript-assay-panel` for endpoint RT-PCR and
     transcript-aware SYBR/TaqMan panels;
   - `primers import-external-pairs` for supplied primers;
   - `primers compose-gene-assay-routine` for evidence-backed composition.
   For a gene-level isoform study, first run
   `primers plan-gene-isoform-study ... --normalize-only` to materialize the
   fully defaulted, content-bound request. The command returns that request
   directly; persist stdout verbatim. Obtain approval for the proposal
   that consumes those exact normalized bytes. The resulting plan may emit an
   ordered workflow; obtain a second approval whose proposal binds that exact
   workflow, its `approved_workflow_sha256`, and its operation-batch digest.
   Execute it only through
   `primers execute-gene-isoform-study-workflow PLAN.json WORKFLOW.json`, which
   fails closed before execution on any mismatch. Never regenerate assay
   operations after the second approval.
   For multiple studies, retain each separately approved plan/workflow pair,
   compose them with
   `primers compose-gene-isoform-study-workflow-batch REQUEST.json`, and obtain
   one second-stage approval for the emitted batch. Execute only that batch via
   `primers execute-gene-isoform-study-workflow-batch BATCH.json
   --checkpoint-dir .gentle/isoform-study-checkpoints`. GENtle
   prevalidates every referenced plan/workflow and the combined ordered digest
   before the first operation, runs in a detached engine, checkpoints successful
   approved operations, and commits project state only after complete success.
   Do not loop over independently state-bound execution proposals.
   Checkpoint existence is never permission to reuse it. For an extended or
   retried batch, first call `primers inspect-gene-isoform-study-reuse ...` and
   show the returned reusable/remaining gene and operation counts plus warnings.
   Only after the user approves that exact `proposal_sha256` may the skill use
   `--reuse-proposal` with `--approve-reuse-sha256`. Never repair a mismatch,
   substitute another checkpoint, or interpret transfer approval as biological
   approval. A changed GENtle executable, baseline state, or approved prefix
   requires fresh computation.
   Set the delegated `execution_timeout_secs` before approval. Treat it as an
   operational wrapper ceiling, not a Primer3 parameter or a predicted wall
   time. Review transcript count, mature-cDNA spans, automatic and requested
   target regions, assay mode, and operation count; sequence length alone does
   not determine the number or difficulty of Primer3 calls. Single-study
   execution defaults to 7200 seconds and a multi-study execution batch to
   28800 seconds, both explicitly overridable in the proposal.
7. Inspect the persisted GENtle report. Preserve all warnings, uncovered
   transcript classes, unresolved junctions, cDNA-reach cautions, and review
   states.
8. When specificity is required, prepare genomic-DNA and transcriptome/cDNA
   handoffs separately. A `specificity-plan` result is a plan, not evidence.
   An external runner may execute BLAST; GENtle must then import/finalize the
   result before the corresponding specificity space can pass.
9. For multiple genes or sequences, use GENtle's collection specificity
   operation rather than looping through ad hoc shell commands.
10. Present order preparation only after the required GENtle review and
    specificity gates are complete. Never place an order.
    `primers oligo-order from-experimental-handoff` additionally requires the
    approved handoff SHA-256 and accepts only order-ready rows under the
    handoff's embedded named policy. Preserve the handoff/policy/readiness-row
    provenance rather than copying an `order_ready` label.
11. For a reviewed transcript panel, use `primers experimental-handoff` with
   JSON, order-table, and virtual-gel outputs. The default skill route requests
   one comparable gel lane per primer pair under one shared GENtle gel model.
12. Preserve the route's `gentle.clawbio_skill_delegation.v1` object. The
    generic runtime must verify this skill's co-deployed descriptor, catalog,
    version, intent, plan step, and delegate contract before GENtle runs.

For final presentation, ask GENtle to run
`primers publish-gene-isoform-study`, or call the equivalent MCP tool
`gene_isoform_assay_publication`. The MCP `confirm=true` is only the ordinary
file-write confirmation. The canonical report declares the only selectable
`content_blocks[]` and named `profiles[]`; ClawBio may choose those ids, host
the generated pages, or invoke the PDF projection, but it must not rewrite
scientific text or add a block. Profile/block selection is a presentation
projection of an immutable report and needs no additional scientific approval.
The projection receipt binds every emitted file.

Every run writes `gentle.clawbio_execution_manifest.v1` to
`reproducibility/execution_manifest.json`. It content-binds the resolved
request, source route, input files, state, runtime, native result, and artifacts.
It records which route was selected; it does not claim to reproduce an LLM's
natural-language routing decision.

## Source Attribution

PCR routes use strict claim attribution. Human-facing statements use these
prefixes:

- `[gentle]`: copied or deterministically projected from the GENtle executable
  result;
- `[clawbio]`: orchestration, artifact collection, or graphical arrangement
  performed by the wrapper, never a biological finding;
- `[input]`: caller-supplied assumptions and paths;
- `[external:<tool>]`: a result produced directly by a named external tool,
  only when its version, inputs, outputs, and provenance are present;
- `[unverified]`: unsupported prose, which cannot satisfy a readiness,
  specificity, or review gate.

The wrapper writes `gentle.clawbio_claim_ledger.v1` to `claim_ledger.json`.
That ledger records the exact GENtle command, request digest, processing tools,
source prefix, JSON evidence pointer/scope hash, and deterministic projection
for each displayed summary line, plus artifact-level authority
(`gentle_executable` versus `none_presentation_only`). Raw GENtle JSON is
retained unchanged; ClawBio does not rewrite it into new scientific claims.

Caller-declared preparation details use request `input_claims[]`. They appear
as individually traceable `[input]` rows in the ledger and report. A declared
oligo-dT primer, cap-capture protocol, or 5-prime completeness observation is
therefore visible but is not silently promoted to a GENtle result.

## cDNA Priming And 5-Prime-End Capture

Treat these as two independent experimental axes:

- **Reverse-transcription priming** describes how cDNA synthesis begins:
  oligo-dT, random hexamers, gene-specific primers, or a mixture. For oligo-dT
  cDNA, GENtle can report the annotation-derived distance from the mature
  transcript 3-prime end to the most 5-prime base required by an assay. That
  is geometry, not a measurement of reverse-transcription completeness.
- **5-prime capture or enrichment** describes how a capped transcript end was
  selected or tagged. Cap-dependent 5-prime RACE/RLM-RACE, cap trapping or
  oligo-capping, CAGE-like 5-prime-end capture, and template-switching
  protocols have different evidence scopes. Record the exact protocol and
  controls as `[input]`; never infer one method from another.

Cap-dependent 5-prime RACE can support a capped 5-prime end compatible with a
transcript start site. It does not by itself quantify the whole transcript or
prove the complete isoform between that end and the poly(A) tail. Cap trapping
can enrich full-length cDNA-RNA hybrids, whereas CAGE ordinarily establishes
5-prime-end/TSS evidence rather than a routine full-length PCR template.
Template switching must not be called cap-selected or full-length unless the
actual protocol and controls establish that property.

In a typical cap-dependent RNA-ligase-mediated 5-prime RACE workflow,
phosphatase treatment makes uncapped or fragmented 5-prime ends unavailable
for adapter ligation; decapping then exposes a ligatable phosphate specifically
on the previously capped RNA population. An RNA adapter, an adapter-specific
primer, and nested gene-specific primers can then enrich the intact capped
5-prime end. Exact enzymes and controls are protocol-dependent and must be
recorded rather than inferred. In cap-trapping workflows, cap selection is
combined with nuclease selection of cDNA-RNA hybrids that reached the capped
end; this is a different preparation contract from ordinary oligo-dT cDNA.

Do not fabricate a universal RACE adapter or pretend a conventional primer
pair is a complete RACE assay. A cap-dependent 5-prime RACE design needs the
protocol-specific adapter sequence plus one or more gene-specific, often
nested, primers. GENtle may design or assess the gene-specific component only
when the reviewed operation describes the available sequence and constraints;
the adapter chemistry remains an explicit `[input]` claim. Recommend short
5-prime and 3-prime diagnostic assays and sequencing of the RACE product when
they are needed to distinguish RT falloff, alternative starts, and full
isoform structure.

## Graphical Assembly

The skill asks GENtle to generate scientific graphics and only arranges the
returned artifacts. Its default experimental handoff includes a virtual gel
whose columns use one primer pair per sample lane, a shared ladder, and shared
gel conditions. Empty lanes remain visible when GENtle predicts no product.
The gel is a deterministic prediction, not a measured abundance result; default
visualization conditions do not satisfy the separate gel-resolution readiness
gate.

## Interpretation Rules

- Expression, PSR, JUC, or other array evidence may prioritize an assay region.
  It does not prove transcript usage and does not replace primer thermodynamics
  or specificity evidence.
- Keep total-expression evidence separate from feature-level residual or
  junction evidence.
- An imported primer remains an imported claim until GENtle recomputes the
  applicable assay and specificity evidence.
- `not_assessed` is not a pass. Missing evidence is not negative evidence.
- A design whose required specificity is still `not_assessed` is not
  order-ready.
- Genomic carryover/off-target evidence and transcriptome/cDNA
  cross-amplification evidence are separate target spaces and must remain
  separately identified.
- Endpoint-gel abundance is rough or semi-quantitative, not quantitative.
- A missing long 5-prime endpoint product from oligo-dT cDNA does not establish
  isoform absence.
- A cap-dependent 5-prime product supports only the evidence scope established
  by its protocol and controls; it does not automatically prove full-length
  isoform abundance.
- Never call a panel order-ready while required specificity, duplicate review,
  or human review remains incomplete.

## Operational Guardrails

- Do not install packages, download large references, delete files, publish
  results, submit orders, or store purchasing credentials.
- Do not invent organisms, genome builds, annotation releases, transcript
  identities, cDNA-synthesis methods, product limits, or specificity resources.
- Use explicit state and output paths and retain all GENtle provenance files.
- Use committed synthetic/offline examples only after an explicit request for
  a demo, example, sample, or synthetic run.
- Do not turn request prose into arbitrary shell commands. Use only the
  descriptor's fixed request templates and the parser-validated
  `gentle-cloning` modes.

## Compatibility

`gentle-cloning` remains the generic execution substrate and can still be
invoked directly for low-level GENtle work. Explicit primer-language discovery
belongs to this specialized descriptor; its plan steps deliberately retain
`skill: gentle-cloning` so no second execution engine is introduced.
