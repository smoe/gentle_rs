# Internal.12 Historical Tutorial Baseline Acceptance

Recorded result reviewed 2026-10-07. This replaces the completed forwardable
regeneration request; it is historical evidence, not current-candidate release
acceptance.

## Reported Result

Glen reproduced the checkout failure at `f1c904a4` and repaired generated
chapter 08.03, the tutorial hub and ledger in
`e462de1fe0bf1d07048610d010b3302fdae740aa`, integrated on 2026-10-02.
He reported Rust 1.96.1 exact-commit LF and CRLF replays, each covering parity,
51 examples and 29 chapters, plus 17/17 Python checkout regressions.
Two stale human-review warnings remained visible. The
[dated changelog](CHANGELOG.md#2026-10-02---separate-historical-baselines-from-current-replay)
and [release ledger](release_notes/release_notes_v0.1.0-internal.12.md)
retain this externally reported result; this documentation cleanup did not
rerun it.

The three synthetic PATZ1 panel reports keep their original `.11` generator
identities and raw hashes. Current replay may differ only in the two narrowly
scoped generator fields described by [testing](testing.md) and
[architecture](architecture.md). Comparison validates original file checksums
before the scoped ledger projection; it does not rewrite reports, strip
scientific differences or normalize hashes.

## Still Separate Or Pending

- Re-run current-candidate version, tutorial, parity and LF/CRLF gates with an
  executable built from that exact SHA. A historical pass is not a pass for
  today's source.
- Retain the remaining promoter-manifest command-path checks: all six required
  components present beside the manifest, no missing-component fallback to an
  unrelated cwd file. Do not promote a synthetic path test to real-artifact
  acceptance.
- Current-source native Windows/macOS CI, independent extracted-package and
  GUI/scientific acceptance remain separate. The later `.12` publication and
  CI package cycle are [recorded](release_notes/release_notes_v0.1.0-internal.12.md#publication-evidence),
  not reopened by this historical record. LF/CRLF simulation on Linux does not
  establish Windows file/path or packaging behavior.

The active [candidate-preparation ledger](release_notes/release_notes_v0.1.0-internal.13.md#internal13-acceptance-ledger),
[release checklist](release.md) and [roadmap](roadmap.md) own those checks and
one-SHA receipts. Keep tags/assets and retained scientific baselines unchanged.
Git history preserves the original request and exact command list.
