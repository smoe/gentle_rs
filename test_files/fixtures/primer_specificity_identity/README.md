# Frozen Pre-Profile Primer Specificity Identities

Origin: entirely synthetic, hand-crafted oligos, transcript, reference, HSP,
alignment, product, execution and panel records. No RefSeq sequences, external
execution, NCBI parity or biological acceptance is represented. The artificial
report ID is a fixture label, not a calculated scientific assessment ID.

`records.json` was captured from baseline
`b62d21c971c207704b96396523eec0ce61659806` with only the test harness added,
before including the unwired mismatch classifier in the engine module tree or
making any production behavior change.
Its nested values are the baseline serializers' output, not records rebuilt
with the serializers under test. `SHA256SUMS` binds the retained LF bytes.

Used by `src/engine/analysis/primer_specificity_identity_tests.rs` to reopen
typed records and check policy/report serialization hashes, actual single/multi
handoff IDs, multi-summary digests, panel design and selection hashes, and the
real selected-report lookup. One complete alignment and one absent alignment
exercise both shapes; no external reference files or tools are opened.

Recreation requires that exact baseline in a disposable checkout, adding only
the test file above and its cfg(test) module inclusion in
`src/engine/analysis/primer_specificity_multi.rs`. With an unused output path:

```sh
GENTLE_PRIMER_IDENTITY_CAPTURE_DIR=/tmp/NEW-primer-identity-baseline \
  cargo test --locked --offline --no-default-features --lib \
  engine::primer_specificity_multi::identity_tests::capture_primer_identity_baseline \
  -- --ignored --exact --test-threads=1
```

The ignored capture test refuses an existing output directory. Review and copy
the two generated files into this directory only when intentionally establishing
a new, separately identified baseline. Never regenerate these files to make a
serialization regression pass. Normal tests only read them. Checkout-byte
coverage in `scripts.test_tutorial_checkouts` checks both LF/CRLF Git modes and
a missing-policy negative control, independently of Rust builds.
