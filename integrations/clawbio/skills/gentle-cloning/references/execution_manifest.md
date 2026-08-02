# Content-Bound ClawBio Execution Receipt

Every `gentle-cloning` invocation writes
`reproducibility/execution_manifest.json` with schema
`gentle.clawbio_execution_manifest.v1`. This is a topic-neutral provider
execution receipt. It is not a discussion record, permission record, evidence
ledger, or scientific interpretation.

## What It Binds

- The normalized `gentle.clawbio_skill_request.v1` payload and its SHA-256.
- The original request-file bytes when a request file was supplied.
- Known file inputs: workflow, plan, catalog, explicit `@file` shell inputs,
  and caller-declared `input_bindings[]`.
- The state file before and after execution, including the explicit absence of
  a state file.
- The `gentle-cloning` wrapper, catalog, and intent descriptor bytes. A legacy
  script-only copied layout remains usable for ordinary direct calls but is
  marked `contract_files_complete = false`; delegated execution requires all
  three contract parts and fails closed when either metadata file is absent.
- The resolved runtime command, content hashes for local launcher/image/binary
  files that appear in the resolved command, and the GENtle version returned by
  the mandatory preflight for delegated runs.
- Every executed command step, with timing, exit state, and stdout/stderr
  hashes.
- The native GENtle JSON payload hash, native report/run/operation identifiers,
  and native status fields without reinterpreting them.
- Collected, scientific, report, replay, environment, and claim-ledger
  artifacts by path, size, and SHA-256. The outer
  `reproducibility/checksums.sha256` also binds the final `result.json` and the
  execution receipt itself.

The receipt deliberately separates `execution_outcome` from
`native_result.reported_status_bindings`. A completed computation may validly
return a biological `fail`, `unsupported`, or `not_assessed` result. The wrapper
must not convert that scientific verdict into an execution failure or a pass.

## Explicit Input Files

Use `input_bindings[]` when a request refers to an input path that is not an
`@file`, `workflow_path`, `plan_path`, or catalog path:

```json
{
  "binding_id": "external_primer_pairs",
  "path": "inputs/primer_pairs.json",
  "role": "external_primer_pair_source",
  "media_type": "application/json",
  "expected_sha256": "sha256:0123456789abcdef0123456789abcdef0123456789abcdef0123456789abcdef"
}
```

The wrapper hashes the resolved regular file before execution and verifies the
same bytes again afterward. A missing file or an `expected_sha256` mismatch
fails before GENtle is invoked; mutation or removal while GENtle runs produces
`verification_failed` and records the post-execution observation beside the
original digest. Paths remain useful locators, but only the observed content
hash identifies the input.

## Delegated Skills

A descriptor-only skill may add a `delegation` object to the normal wrapper
request:

```json
{
  "schema": "gentle.clawbio_skill_delegation.v1",
  "source_skill": "gentle-pcr-primer-design",
  "source_skill_version": "0.3.0",
  "intent_id": "primer_report_list",
  "plan_step_index": 0,
  "resolved_slots": null
}
```

Before GENtle runs, `gentle-cloning` requires the source skill's co-deployed
`INTENTS.json` and `catalog_entry.json`. It verifies their identities, versions,
hashes, selected route and plan step, and the source catalog's exact
`delegate_contract`. An incomplete or stale two-skill deployment therefore
fails closed. Optional `descriptor_sha256` and `catalog_sha256` values can pin
an expected deployment even more tightly.

`resolved_slots` is optional routing provenance. The canonical resolved wrapper
request is always bound independently. The receipt records the selected route;
it does not claim that an LLM or chat router will select the same route from the
same prose on replay.

## Exact Claim Attribution

For `claim_attribution_mode = "strict"`, `claim_ledger.json` binds every
`[gentle]` summary to:

- an exact JSON pointer or exact source scope;
- a SHA-256 of the pointed value and source scope;
- a named deterministic projection identifier.

If those bindings are unavailable, prose cannot receive `[gentle]` authority.
Caller text is always `[input]`; a caller-supplied leading `[gentle]`,
`[external:...]`, or other reserved prefix is escaped as literal text. Wrapper
layout and orchestration remain `[clawbio]` and do not become scientific
findings.

## Ownership Boundary

An external workflow may store this receipt and GENtle's native report as
content-addressed evidence, then create its own run record. It should compute
its own request digest using its owning contract and verify current permissions
and evidence refs before recording. GENtle does not emit or accept those
foreign concepts. For discussion-moderation specifically, the optional
[post-run handoff adapter](discussion_moderation_handoff.md) prepares the
expected canonical digest fields and typed output shape, while leaving the
actual validation and mutation to discussion-moderation.

Deployment caches, chat-retention policy, permissions, withdrawal,
supersession, and interpretation remain responsibilities of the caller and
operator. The receipt makes version or content drift detectable; it cannot
force an external OpenClaw process to refresh a cached skill bundle.
