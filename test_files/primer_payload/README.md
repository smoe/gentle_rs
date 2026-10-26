# Invented Primer Payload Evidence

Origin: hand-crafted deterministic oligos, reads, PAF and reference annotations in
`scripts/primer_payload_rehearsal.py` and inline Rust tests in
`src/engine/analysis/rna_panel/primer_payload/`. No reads from a physical assay,
external catalog, current panel or supplier order are used. The two core lengths
are 10/12 nt, not the motivating 44-mer oligos; read names deliberately repeat.
The second raw read is reverse-complemented. Literal PAF rows are teaching
evidence, not the output of a reproduced aligner/search.

Inline mask regressions additionally enumerate abstract intervals over 12
positions (12,168 ordered anchor/hit/orientation cases) and coordinate extremes.
These are geometry inputs, not biologically admissible oligos. A temporary
two-read/four-anchor fixture checks shared preparation/verification in both
orientations, including verification after the original FASTQ is removed.

`primer_payload/binding_tests.rs` uses those temporary invented sources to
deliberately alter/re-hash only new query manifests/FASTA, then checks canonical
record/policy accounting and source-byte preservation. Fixed anchor-identity
digests pin the original JSON encoding (also calculated independently in
Python); they are serialization drift guards, not source-origin approval. The
optional actual-CLI companion includes a coherently rehashed false-exclusion
refusal under LF/CRLF with spaced/apostrophe paths; no local Rust run is implied.

Recreate in a **fresh** synthetic directory using the Phase 2 addendum to
tutorial 10.04. The helper has no GENtle build, subprocess, network or implicit
input-discovery capability. It creates new fixture files with exclusive opens;
the next stage requires the actual preceding GENtle report. Capture SHA-256 for
the raw generated files and reports on the tested platform. Filesystem-dependent
absolute paths mean there is no single claimed portable raw report hash.

Consumers: shared preparation/PAF/annotation/inspection engine tests, optional
actual-CLI rehearsal and LF/CRLF/source checks in
`scripts/test_primer_payload_contract.py`. These tests do not certify native GUI,
live-model behavior, real-library accuracy, throughput, molecule counts or a
biological mechanism. Phase 1 inputs/artifacts are not edited by Phase 2.
