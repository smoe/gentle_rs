# External Primer Evidence Handoff

Use `mode=external-primer-handoff` when a caller has already selected and is
explicitly supplying external primer or oligo records. The mode is a bounded
wrapper around GENtle's existing `primers import-external-pairs` operation. It
does not add another provenance model and does not discover records from a
caller-owned store.

## Ownership Boundary

GENtle already owns the biological identities and computed results:

```text
external source fields
  -> native source_record_id
  -> normalized native pair_id
  -> cDNA assay/result and report_id
```

The skill wrapper adds a verifiable transport join around that native chain.
It accepts explicit records, invokes the existing importer, verifies the
returned result, and exposes which caller `record_id` maps to which native
source, pair, assay, and report IDs. The caller remains responsible for data
permission, retention, withdrawal, supersession, and later interpretation.

The wrapper never:

- opens or searches a caller-owned store by a guessed path;
- treats publication or catalogue use as experimental validation;
- treats a computed result as a participant statement or vote;
- forwards private conversation text;
- chooses an experiment or writes an interpretation back to another system.

## Request Contract

The outer request uses `gentle.clawbio_skill_request.v1`, an explicit
`state_path`, and a nested `gentle.external_primer_handoff_request.v1` object.
See
[`request_external_primer_handoff_synthetic.json`](../examples/request_external_primer_handoff_synthetic.json)
for a complete, synthetic example.

The nested object contains:

- `collection_id`: stable caller identity for this supplied collection;
- `target`: required `seq_id` and `source_feature_id`, plus optional
  `transcript_id`, `reference_label`, `reference_release`, and
  `expected_state_sha256`;
- `evaluation`: bounded arguments supported by the native importer;
- `records[]`: explicit records with a caller-stable `record_id`, source
  fields, assay purpose, and one sequence or one forward/reverse pair.

Supported source kinds are `literature`, `commercial_catalogue`, `laboratory`,
and `external`. Source fields and `validation_claims[]` are retained as claims;
they are not used as evidence of coverage or specificity.

Assay-purpose behavior is deliberate:

| `assay_purpose` | Accepted shape | Wrapper behavior |
| --- | --- | --- |
| `qpcr` | forward/reverse pair | Submit to native external-pair importer |
| `endpoint_pcr` | forward/reverse pair | Submit to native external-pair importer |
| `cloning` | pair or standalone oligo | Return typed `not_submitted` record |
| `sequencing` | pair or standalone oligo | Return typed `not_submitted` record |

GENtle does not currently have a native external standalone-oligo evaluation
contract with assay-purpose semantics. The wrapper therefore fails closed
rather than silently coercing cloning or sequencing oligos into PCR evidence.

## Binding and Verification

Before execution, the wrapper normalizes IUPAC DNA sequences, sorts records
and set-like fields, and computes `canonical_request_sha256`. It writes the
canonical request and the exact native input batch into the output bundle. If
`expected_state_sha256` is supplied, a mismatch aborts before GENtle runs.

The native command omits `--report-id`, so GENtle derives its normal report
identity. After execution, the wrapper verifies:

- native command and report schemas;
- collection, target sequence, and source-feature identity;
- source fields, normalized sequences, source IDs, and claim isolation;
- input-file hash and native normalized-batch hash;
- automatically derived report identity and operation identity;
- equality of the command result, persisted report, and embedded assay JSON;
- the runtime version captured by `gentle_cli --version`;
- hashes and bounded locations of the native report and scientific artifacts;
- the entire explicit state-file hash before and after execution.

Any altered sequence, missing result, wrong target/run binding, digest
mismatch, escaped artifact path, or schema mismatch yields
`verification_failed`. A requested specificity evaluation that is missing,
incomplete, or errored yields `incomplete`; it never becomes a vacuous success.

## Returned Join

The wrapper's `result.json` retains GENtle's native JSON in `stdout_json` and
adds `external_primer_handoff` with schema
`gentle.external_primer_handoff_result.v1`. Important fields are:

- `canonical_request_sha256` and `execution_binding_sha256`;
- `target_state`, including state hashes and caller-declared reference labels;
- `submitted_record_joins[]`, binding caller `record_id` and exact normalized
  sequences to native `source_record_id`, `pair_id`, `assay_test_id`, and
  `report_id`;
- `not_submitted_records[]` for cloning/sequencing inputs;
- `report_binding`, including native batch/report/operation identities;
- `scientific_artifacts[]`, with media type, size, and SHA-256.

The result is suitable for a caller to retain as a compact ID/hash/locator
record while leaving the native GENtle report intact.

## Provenance Audit and Limits

No engine provenance type was added. Before this wrapper, GENtle already had:

- external source classes and source fields;
- deterministic source and normalized pair identities;
- normalized batch SHA-256 and input-file provenance;
- native report and operation identities;
- computed assay, specificity, warning, and artifact records;
- explicit separation between external claims and computed evidence.

The wrapper newly adds only the caller-facing canonical request digest,
target-state binding, runtime capture, artifact hash verification, and explicit
caller-record-to-native-result join.

Two general engine limitations remain visible rather than hidden:

1. Assay purpose and standalone external oligos are not yet first-class native
   external-primer engine records.
2. The wrapper binds the whole state file and preserves caller-declared
   reference labels/releases, but it does not independently prove a dedicated
   target-sequence or annotation-snapshot digest.

The checked-in example and regression tests are wholly synthetic. They do not
contain locally held laboratory primer records or imply permission to publish
such records.
