# GENtle PCR Primer Design Skill

`gentle-pcr-primer-design` is a descriptor-only OpenClaw/ClawBio skill. It
provides focused primer-design intent recognition and operation sequencing
without copying GENtle science or the `gentle-cloning` Python wrapper.

## Execution Model

```text
user request
    -> gentle-pcr-primer-design/INTENTS.json
    -> one or more skill_run steps targeting gentle-cloning
    -> digest-bound proposal and caller approval when confirmation is required
    -> gentle_cloning.py
    -> parser-validated gentle_cli command or workflow
    -> GENtle report, evidence state, provenance, and artifacts
```

ClawBio can discover this descriptor even though the directory has no Python
entrypoint. The descriptor is executable when `gentle-cloning` is registered,
because each route explicitly delegates to that skill.

Each delegated request carries `gentle.clawbio_skill_delegation.v1` metadata.
Before GENtle runs, `gentle-cloning` verifies this skill's co-deployed catalog,
descriptor, version, selected route and plan step against the exact generic
wrapper contract. A stale or partial deployment fails closed. The resulting
`gentle.clawbio_execution_manifest.v1` records those hashes and the resolved
request; it records the selected route but does not claim that natural-language
routing itself is deterministic.

Natural-language routing therefore produces a draft request, not authority to
run scientific work. Routes that select a target/pair/backend, mutate state, or
write an artifact first return `approval_required` plus a
`gentle.clawbio_execution_proposal.v1`. The caller must approve that exact
digest, after which the generic runner executes the proposal's stored request
without rerouting. Only list/show/preflight diagnostics remain automatic. See
the generic runner's
[execution-approval contract](../gentle-cloning/references/execution_approval.md).

Gene isoform-assay studies deliberately use two such approvals. The first
binds the normalized planner request, policy/version, evidence hashes, prior
plan, observations, override, and resolved defaults. The second binds the
exact ordered `DesignTranscriptAssayPanel` operation set emitted by that plan.
Approval authorizes execution of those bytes; it is not a declaration that the
biological premise is correct.

Primer backends and pair ranks have no conversational default. An explicitly
requested `auto` backend is resolved to a path/version-pinned Primer3 or the
internal backend before approval. An unknown 5-prime capture method is
`unspecified`; `none` is retained only for an explicit assertion that no such
method was used.

The specialized skill therefore cannot:

- calculate or rank primers independently;
- reinterpret a GENtle specificity verdict;
- upgrade `not_assessed` to pass;
- declare an imported primer validated;
- bypass review or order-readiness gates.

## Provenance-First Reports

Every specialized request enables strict source attribution. Display prose is
prefixed with `[gentle]`, `[clawbio]`, `[input]`, `[external:<tool>]`, or
`[unverified]`, and the wrapper writes
`gentle.clawbio_claim_ledger.v1` as `claim_ledger.json`. Only `[gentle]`
statements have an exact source pointer/scope hash and named deterministic
projection from the executable result. ClawBio may arrange figures and
describe orchestration, but it may not invent or upgrade scientific findings.
The ledger separately marks each artifact's scientific-content authority. The
original GENtle JSON remains available unchanged beside the ledger.

The claim ledger is joined to the generic execution receipt in
`reproducibility/execution_manifest.json`. The receipt keeps wrapper execution
outcome separate from GENtle's native biological verdict, so a successfully
computed `fail` remains a completed execution and a scientific fail.

Transcript-panel requests also record reverse-transcription priming,
5-prime-capture method, and supplied completeness evidence as explicit
`input_claims[]`. These appear as `[input]` in the report and ledger. GENtle's
oligo-dT reach calculation remains `[gentle]`; a caller's cap-dependent RACE,
cap-trapping, CAGE, or template-switching protocol does not become a GENtle
finding merely because the skill carried it.

## Supported Intent Families

- conventional PCR design;
- generic TaqMan-style qPCR design;
- endpoint RT-PCR and transcript/isoform panels;
- oligo-dT reach-aware and cap-dependent 5-prime-end assay planning;
- junction-aware primer-only SYBR and transcript-aware TaqMan panels;
- imported-primer review;
- separate genomic-DNA and transcriptome/cDNA specificity plans;
- collection-wide specificity;
- evidence-backed assay-routine composition;
- report review and deterministic PCR figures;
- explicit offline demos and missing-runtime diagnostics.

The `experimental_assay_report` route turns a persisted transcript panel into
one GENtle handoff JSON, an order/readiness TSV, and a virtual-gel SVG. The gel
uses one lane per primer pair under shared conditions so pair columns are
directly comparable. Predicted-empty lanes are retained, and the figure is not
presented as wet-lab evidence.

The downstream isoform-assay publication is still GENtle-owned. Its canonical
JSON embeds each content-bound plan, handoff, and order form once and declares
pointer-only presentation blocks. OpenClaw may request a declared `full`,
`review`, or `ordering` profile (or an explicit subset of declared block ids),
then host the per-gene HTML pages or browser-rendered PDF. It may not create a
new scientific block. Order sheets are accepted only from a specific handoff
digest and its named readiness policy, with card/assay/pair/oligo provenance
retained per line.

## Running The Delegated Examples

The examples are ordinary `gentle.clawbio_skill_request.v1` files consumed by
the existing generic runtime:

```bash
python clawbio.py run gentle-cloning \
  --input integrations/clawbio/skills/gentle-pcr-primer-design/examples/request_primer_preflight.json
```

```bash
python clawbio.py run gentle-cloning \
  --input integrations/clawbio/skills/gentle-pcr-primer-design/examples/request_patz1_transcript_assay_offline.json
```

The second request is synthetic and must only be offered when a user explicitly
asks for a demo or example. See [examples/README.md](examples/README.md) for
fixture provenance and recreation notes.

## Specificity Lifecycle

`primers specificity-plan` prepares a non-executing handoff. A caller or
external runner may execute the emitted BLAST command, after which GENtle
imports and finalizes the evidence. Until finalization, the applicable status
remains `not_assessed`/`not_run`; the skill must not describe the assay as
specific or order-ready.

Genomic-DNA and transcriptome/cDNA reference spaces are represented by separate
plans and output directories. Their results are not interchangeable.
