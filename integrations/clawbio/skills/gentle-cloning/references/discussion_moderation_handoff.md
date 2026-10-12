# Discussion-Moderation Provenance Handoff

`discussion_moderation_handoff.py` converts one completed or attempted
`gentle-cloning` run bundle into a content-bound intake packet for
discussion-moderation's typed-evidence and `analysis-run` commands. It is a
post-run skill adapter, not a GENtle engine operation.

The helper never opens a ledger, discovers a topic, checks participant
permissions, or records a result. The caller supplies the exact permitted
references and remains responsible for validating and recording them through
discussion-moderation's owning commands.

## Why There Are Two Records

The GENtle execution manifest and a discussion-moderation analysis run describe
one event, but they have different owners:

- `gentle.clawbio_execution_manifest.v1` is the immutable provider receipt. It
  records execution paths, timestamps, request and file hashes, runtime files,
  command outcomes, native result bindings, and output artifacts.
- the ledger's `analysis-run` is the policy-aware record. It validates current
  evidence references and permissions, assigns ledger timestamps and attempt
  state, detects replay, and owns dependency freshness.

The handoff packet joins those records by hashes. It must not be copied into a
ledger as though the provider had authority over ledger-owned fields.

## Invocation

Run GENtle first, then supply a small caller context:

```bash
python3 discussion_moderation_handoff.py \
  --result /path/to/gentle-run/result.json \
  --context /path/to/handoff-context.json \
  --output /path/to/handoff-packet.json
```

By default, the helper reads
`/path/to/gentle-run/reproducibility/execution_manifest.json`. Use `--manifest`
only when that receipt is stored elsewhere. `--native-result-output` chooses
the location of the extracted, content-addressed GENtle result artifact.

The helper performs no downloads, network access, package installation, or
ledger mutation. Input JSON is size-bounded and non-finite numbers, unknown
fields, and private-chat/credential metadata keys are rejected.

## Caller Context

The context schema is
`gentle.clawbio_discussion_moderation_handoff_context.v1`:

```json
{
  "schema": "gentle.clawbio_discussion_moderation_handoff_context.v1",
  "topic": "fictional-primer-discussion",
  "request": {
    "kind": "participant",
    "ref": "participant-fixture"
  },
  "run_id": "analysis-gentle-primer-fixture-v1",
  "recorded_by": "moderator",
  "input_refs": [
    {
      "ref": "evidence:external-primer-pairs-fixture-v1",
      "content_hash": "sha256:0123456789abcdef0123456789abcdef0123456789abcdef0123456789abcdef",
      "binding_id": "external_primer_pairs"
    }
  ],
  "output_evidence": {
    "id": "gentle-primer-result-fixture-v1",
    "version": "1",
    "type_uri": "urn:gentle:primer-review-result",
    "schema_version": "1",
    "label": "Fictional GENtle primer review",
    "visibility": "shared",
    "permission": "sender_approved",
    "source_refs": [
      "source-fictional-primer-record"
    ]
  },
  "parameters": {
    "review_scope": "fictional"
  },
  "reference_releases": {
    "transcripts": "fixture-release-1"
  },
  "environment": {
    "runtime_profile": "offline"
  }
}
```

Each `input_refs[]` row must identify a current typed evidence ref and an exact
file binding already present in the provider execution manifest. Its
`content_hash` must equal the SHA-256 of the bytes GENtle bound before and after
execution. This v1 adapter intentionally rejects inferred, fuzzy, path-only,
or unrecorded inline projections. A caller that transforms evidence into a
different file must record and verify that transformation before presenting
the transformed artifact as the run input.

For a failed provider run, omit `output_evidence`. A completed run requires it.
An incomplete run may return a typed output only when a real partial artifact
exists; its validation state remains `incomplete`.

## Verification and Mapping

The output schema is `gentle.clawbio_discussion_moderation_handoff.v1`.
Verification fails closed unless:

- the result and execution-manifest schemas are recognized;
- the embedded and on-disk manifests are identical;
- wrapper status, stdout, stderr, and parsed native JSON hashes agree;
- the wrapper contract files were complete;
- every caller ref resolves to exactly one stable content-bound input whose
  before/after size and SHA-256 agree;
- each caller content hash equals the bytes consumed by GENtle.

The mapping is:

| Provider receipt | Intake packet / ledger destination |
| --- | --- |
| normalized request SHA-256 | `parameters.gentle_provider.normalized_request_sha256` |
| delegated skill/intent/plan hashes | `parameters.gentle_provider.delegation` |
| wrapper/runtime identities and hashes | `environment.gentle_provider` |
| handoff adapter ID, version, and script SHA-256 | `environment.gentle_handoff_adapter` |
| state-before SHA-256 | `environment.gentle_provider.state_before_sha256` |
| reference-preflight payload SHA-256 | `reference_releases.gentle_reference_preflight` |
| exact caller refs plus file hashes | `input_descriptors[]` and `verified_input_bindings[]` |
| extracted native result SHA-256 | `output_evidence_object.content_hash` and `artifact.content_hash` |
| native report IDs and status pointers | extracted native result artifact |

The request digest deliberately excludes filesystem paths, run timestamps, the
raw request-file hash, and the execution-manifest hash. It follows
discussion-moderation's canonical digest over `topic`, `request`, `routine`,
content-hashed inputs, parameters, reference releases, and environment. The
raw request-file and manifest hashes remain in `provider_receipt` for audit.
This keeps an exact replay idempotent even when its bundle is relocated or its
equivalent canonical request is serialized differently.

The handoff adapter's own ID, semantic version, and script SHA-256 are part of
the stable request environment and are repeated in `provider_receipt`. A change
to the transformation code therefore produces a different request identity;
relocating unchanged code does not.

## Status and Claim Semantics

Execution, evidence validation, and biological findings are separate:

- `analysis_run_intake.outcome` reports completed, failed, or incomplete
  provider execution.
- `output_evidence_object.validation` covers provider integrity and execution
  only.
- native GENtle statuses remain in `stdout_json`,
  `reported_status_bindings[]`, and the strict claim ledger.

A completed computation that reports `specificity_fail` remains a completed
run with a computed failing scientific result. It is never converted to a
scientific pass. Failed runs have no output evidence, and `not_run`,
`not_assessed`, missing, or incomplete evidence is never promoted to success.
The output evidence has `claim_level = computed`; it is not participant support
or experimental validation.

## Recording Order

The consumer must use discussion-moderation's commands rather than editing
ledger JSON:

1. Verify that every input ref is current, usable, and permitted for analysis.
2. Validate and record `output_evidence_object` when present.
3. Record `analysis_run_intake`, supplying its `request_digest` so the ledger
   recomputes and verifies it.
4. Let the ledger assign timestamps, attempts, dependency snapshots, freshness,
   and idempotent-replay state.

The packet's `input_descriptors[]` is a digest-verification aid, not an extra
stored analysis-run field. The ledger resolves the current input records and
reconstructs those descriptors itself.

## Remaining Trust Boundary

The helper proves internal hash joins; it does not cryptographically sign the
provider. Artifact authenticity still depends on trusted execution/storage or
a future signature layer. It also cannot approve permissions, infer
withdrawal/supersession, establish biological correctness, or decide how a
computed result bears on a discussion.
