# GENtle `v0.1.0-internal.11` Release Candidate

Status: **unreleased; not yet green**. This document defines the next candidate,
not a claim that its acceptance or packaging has passed. `.10` is an existing
tag; its historical results do not certify `.11`.

## Intended Story

- Inspect TP73 genome-anchored tracks and PATZ1 negative-strand transcript,
  motif, occupancy and assay evidence through shared engine/GUI contracts.
- Save assembly-bound genomic regions; exchange reproducible Conservation
  requests, navigate alignment blocks and cancel active BLAST work.
- Resolve hash-bound external locus/reporter evidence before proposing
  regulatory-fragment panels. Separately approve exact, atomic multi-fragment
  products with explicit source, annotation and sequence provenance.
- Reproduce branch/reverse-complement, digest and Simple-PCR tutorials through
  ordinary Linux X11 input with native-window/coordinate binding and retained
  screenshot provenance. Starter projects must not contain completed results.
- Distribute a downloadable Linux x64 tarball alongside the Windows ZIP,
  macOS DMG and separate Linux CLI/GUI GHCR images.

Glen's `c7f3f005` runner fixes are integrated. The Simple-PCR starter now uses
the same 800-base TP73 extract as its scripted oracle, with core [200, 600)
and at most 200 bases per flank. The full source locus is retained for
provenance. This bounds the input rather than relaxing scientific gates or
raising the ten-minute compute timeout. A new live run is still required.

## Exact-Candidate Gate Ledger

Before evaluation, record the full candidate SHA, clean-tree status,
`Cargo.lock` SHA-256, Rust version, build features/profile, host and relevant
fixture/resource hashes in the evidence bundle. Each ledger row must link its
receipt and identify that same revision. Do not aggregate passes from nearby
commits. After a fix, select the new candidate and rerun its gates.

| Gate | Owner | Status |
| --- | --- | --- |
| `cargo check -q --locked` and `cargo test --locked --workspace` | Glen / CI | Pending |
| Release-shaped `script-interfaces` build and every published entrypoint smoke | Glen / platform CI | Pending |
| Examples, tutorial generation/manifest/catalog checks, interface parity | Glen / CI | Pending |
| All three Linux tutorial GUI smoke chapters, including bounded Simple PCR | Glen | Pending; new runtime measurement required |
| TP73/PATZ1 graphical inspection and Conservation navigation/cancellation | Glen | Pending |
| Copied-state IRF9/Q00978 acceptance | Glen, private copied-state evidence | Pending |
| GUI and specificity benchmark acceptance | Glen | Pending |
| TP73 CUT&RUN/evidence-viewer and PATZ1 locus-composer proof | Glen / CI | Pending |
| Both Docker targets, `runtime-cli` and `runtime-gui` | Container CI | Pending |
| Linux tarball, Windows ZIP and macOS DMG from the exact tag candidate | Release CI | Pending |
| macOS optional `screenshot-capture` compilation | macOS CI | Pending; Linux cannot validate this |
| Clean tree, version/tag/SHA consistency and generated-artifact checks | Release owner | Pending |

Use [the release checklist](../release.md) for entrypoint commands,
[GUI acceptance instructions](../testing.md#61-tutorial-gui-acceptance-contracts)
for the three-chapter smoke, and the
[benchmark contract](../testing.md#51-external-performance-audit) for
performance evidence. The TP73 proof is described in the
[evidence-viewer runbook](../tp73_genome_evidence_viewer_runbook.md) and
[CUT&RUN smoke guide](../cutrun_release_smoke.md). Retain private inputs outside
the repository; a missing private gate is pending, not silently waived.

The copied-state assay acceptance must retain the expected one-of-three linked
cDNA / two-of-four linked-record coverage, with only the actual patch record
genomically unassessed. GUI/report evidence must not rewrite that scientific
oracle to obtain a pass.

## Known Limits

- Products are **designed, not experimentally validated**. Sequence similarity
  does not establish orthology, TF binding, regulatory activity or sufficiency.
- A candidate-region BLAST database does not establish whole-genome uniqueness.
- The compact PCR tutorial is a navigation/design smoke, not a whole-locus
  performance benchmark or genome-wide specificity assessment.
- Missing whole-genome, transcriptome or private evidence remains explicit.
- The Linux tarball targets Ubuntu 24.04 x86-64 and compatible runtime libraries;
  it is not a static universal binary, AppImage or Debian package. BLAST+,
  Primer3 and other external scientific tools remain separately installed.
- Linux screenshot acceptance does not prove macOS ScreenCaptureKit behavior.
- Screenshot files explain interaction; typed sequence/report checks remain
  the scientific oracle. No background capture authority is added.

No tag or release publication is authorized by this document. Keep `.11`
unreleased until its owner has reviewed the exact-candidate evidence.
