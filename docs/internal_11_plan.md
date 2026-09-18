# Internal .11 Integration And Transcript Capture

Planning baseline: 2026-09-18, `origin/main` at `84f34a9e`.
Development branch: `codex/internal.11`. The `.10` release remains owned by
`main`; this document neither declares it released nor waives its gates.
No version bump, tag, push or publication is part of this integration.
Claude consultation was offered; this plan is not Claude-reviewed.

## Integration Boundaries

1. Preserve Glen's tutorial history by merging
   `bot/tutorial-human-readable-20260917` (`cfd83bcc`). Review generated and
   hand-written teaching paths together. Correct fallback headings that cut
   accession versions, paths or ranges, and distinguish terminal wrappers
   from GUI Shell commands and GUI effects. Keep inner-agent examples static
   and review-only. Do not claim live GUI or biological sign-off from replay.
2. Port only `14f00836` (selectable-text vector PDF) and `00c9fa93` (multipage
   vector composition) from `glen/vector-pdf-refseq-20260915`, preserving
   attribution. Keep `pdf` as the existing raster compatibility format and
   `vector_pdf` opt-in. Retain SVG peers for hover; verify actual PDF text
   extraction, page order/dimensions, used-font provenance and failure cleanup.
3. Keep the older RefSeq commit `03304b13` outside the PDF integration.
   Current `src/transcript_presentation.rs::gff3` still requires a literal
   assembly label. The proposed accession-bound patch exception addresses a
   real gap, but hard-codes `GCF_000001405.40`/`GRCh38.p14` and reads only the
   first build/accession header. A separate change should reuse catalog-bound
   assembly identity and reject conflicting repeated headers. Do not broaden
   assembly compatibility as a side effect of enabling a PDF backend.

The starting branch also retains `ea6f9009`, the Windows ClawBio fixture fix
awaiting integration into `.10`; it is not a new `.11` feature.

## Verified Capture Baseline

See [the current contract](transcript_capture_pools.md),
`src/engine/protocol/transcript_capture.rs`,
`src/engine/analysis/transcript_capture.rs` and its adjacent tests.

- `DesignTranscriptCapturePool` and
  `primers design-transcript-capture-pool @REQUEST.json` already share a typed
  request and stored report. Shell list/show/export and introspection exist.
- Discovery searches explicit transcript windows, retains exact bindings and
  cDNA/source digests, and jointly selects a bounded pool across targets.
  Roles, tails, stage co-presence, primer budgets and sharing groups are data.
- Coverage, omitted bases, retained sequences and stage-specific full-oligo
  interactions are reported. An identical cDNA digest is not itself proof of
  primer coverage; terminal-exon binding does not imply intact CDS capture.
- There is no dedicated capture form or reference-bound capture specificity
  assessment. `specificity_status` is currently `unassessed`. An A-run warning
  is neither repeat-family annotation nor evidence of genomic specificity.
- General reuse seams already exist: primer specificity hit/expected-product
  records and local BLAST machinery, `OligoOrderForm` with duplicate review,
  and provider-neutral external-service preflight/export. Extend their joins;
  do not build a parallel aligner, procurement model or GUI-only algorithm.

## Small Implementation Slices

### A. Capture Workspace And Reproducible Request

Add a dedicated view reached from selected annotated DNA and the shared
`ui ...` intent surface. Use the existing typed capture request for target
sources, transcript universe, strand-aware windows, roles, tails, stage
co-presence, fixed oligos and budgets. Do not silently choose a longest
transcript or reinterpret terminal exon as 3'UTR.

Show source/annotation identity, missing CDS/UTR information, coverage and
omitted/retained intervals before showing the selected pool and interaction
matrix. Label specificity and order readiness as unassessed until their own
checks exist. Support request/report export, stale-input detection and the
existing managed-command path; do not run discovery in a repaint callback.
Use semantic GUI identifiers and include CLI/GUI request-equivalence tests.

Acceptance: plus/minus synthetic multi-isoform loci; unknown CDS; empty windows;
fixed degenerate oligos; joint 5'/3' stages; cancellation; save/reopen; export
parity. A GUI success is not an assay or ordering approval.

### B. Pool-Specific Reference Assessment

Bind an assessment to the selected pool/request digests and prepared genome
and transcriptome identities, assembly/release, index hashes, tool version and
effective settings. Reuse the local BLAST/preflight and full-binding machinery.
Only annealing segments are sequence-search queries; retain tails in provenance
and full-oligo interaction checks.

Define explicit allowed transcript subjects and genomic intervals per target,
including authorized multi-locus sharing. A matching gene name alone must not
permit a hit. Assess individual binding and compatible co-present oligo pairs;
do not reinterpret an RT/capture-only stage as an ordinary PCR pair. Preserve
unresolved, truncated, missing-index and stale-evidence states rather than
promoting them to pass. Keep genome and transcriptome outcomes separate.

Join assembly-bound rmsk/repeat evidence to binding spans. Distinguish intended
shared targets, incidental Alu/other repeat overlap, repeat-mediated off-targets
and missing repeat evidence. No universal repeat exclusion or silent hard
filter is implied. Resolve policy in the typed request before implementing it.

Acceptance: permitted versus unauthorized multi-locus hits, splice junctions,
opposite strands, 3' mismatches, repeated sequence at distinct loci, missing
repeat data, wrong assembly/release, stale hashes and incomplete search. Use
small synthetic references first; Glen validates real prepared references.

### C. Evidence Ranking And Reviewable Order Handoff

Attach optional microarray/expression report identities and quantitative values
to the relevant transcript/locus with explicit mapping ambiguity, tissue,
condition and normalization provenance. Initially present these as ranking
annotations; missing evidence is not zero expression. Any later weighted
optimization must be explicit in the request and must not silently change the
coverage universe or override specificity/retention constraints.

Create an explicit capture-report-to-`OligoOrderForm` handoff containing exact
full oligos, candidate/source bindings, assessment hashes, modifications and
unresolved warnings. Reuse duplicate/sequence-reuse review and provider
preflight/export. Human review binds to the exact form; changing a sequence,
tail, stage or evidence invalidates that review. Export is not submission or
experimental approval, and automatic ordering remains out of scope.

Acceptance: unavailable/ambiguous expression data leaves coverage unchanged;
duplicate ordered sequences retain all intended roles; changed pool/form
invalidates approval; missing specificity remains visible; exports match the
reviewed sequences byte-for-byte and never contact a supplier implicitly.

### D. Public Replay And Independent Acceptance

Begin input pinning alongside the GUI slice; repeat the complete replay after
specificity and handoff integration rather than waiting until then to discover
missing experimental inputs.

Prepare a public, accession/release/hash-pinned replay for TP73, E2F1, TP53,
TP63, PATZ1, FUS, MDM2 and MDM4 using the same engine request. Record download
origins/recreation, selected isoforms, CDS/UTR availability, exact adapters,
stages, biological retention requirements and constraints. Do not guess these
experimental inputs from gene names or copy a private order sheet into tests.

Optimize 5' and 3' candidates together where they are co-present. Retain all
evaluated candidate/coverage/interactions evidence, declared bounds and
uncovered targets rather than forcing a successful-looking pool. Keep the
small synthetic CI replay distinct from the full-reference system acceptance.

Glen should test one exact candidate with binary/input hashes, profile and
effective parameters: GUI responsiveness and semantic navigation; shared-shell
parity; both reference spaces; repeat distinctions; source-bound ranking and
order export. The user retains experiment design and procurement approval.

## Handoff Gates

- Tutorial generator/catalog/manifest/check tests and Python documentation
  checks, with optional/live checks explicitly reported rather than implied.
- Vector-PDF and TSS export tests, compositor rollback tests, real text
  extraction and visual inspection of public synthetic multipage output.
- `cargo fmt --check`, targeted Rust tests, `cargo check -q` and
  `git diff --check`, followed by session-close hygiene.
- No `.10` acceptance inferred from these `.11` checks. Full workspace,
  packaging, live tutorial/inner-agent and real-data gates remain separately
  revision-bound.

Visual follow-up, not silently included in the PDF port: the synthetic TSS
example's long left ruler label overlaps its first tick in both the existing
raster output and the vector output. Address label placement in the shared
renderer separately; successful PDF conversion is not publication approval.

## Integration Verification (2026-09-18)

Local macOS implementation checks used Rust/Cargo `1.99.0-beta.5`, locked
dependencies and single-job builds. The [changelog](CHANGELOG.md) records counts
and limitations. Tutorial generation/check left scientific JSON/SVG artifacts
unchanged; only chapter prose and chapter hashes changed. Inner-agent examples
were not executed. The two prior stale human-review warnings remain visible.

Repeatable entry points (no private data or live provider required):

```bash
cargo test --locked -j 1 --lib workflow_examples::tests::tutorial
cargo test --locked -j 1 --lib svg_vector_pdf::tests
cargo test --locked -j 1 --lib tss_profile_export::tests
cargo test --locked -j 1 --bin gentle_cli rendering::tests
cargo run --locked -j 1 --bin gentle_examples_docs -- tutorial-check
python3 -m unittest discover -s scripts -p 'test_*tutorial*.py'
python3 -m unittest discover -s scripts -p 'test_compose_locus_tss_profile_pdf.py'
cargo check --locked -q -j 1
cargo fmt --check
git diff --check
```

For retained PDF review, set `GENTLE_TSS_EXPORT_EXAMPLE_DIR` to a fresh absolute
directory and run the ignored
`tss_profile_export::tests::write_synthetic_export_example` test. It creates
the public hand-crafted plus/minus examples, not newly scored biological data.
Feed its ordered SVG pages to `gentle_cli svg-vector-pdf-set`, keeping the
combined PDF outside the receipt-bound export directory. Verify using
`pdfinfo`, `pdffonts`, `pdftotext`, `pdfimages -list` and rendered page images.
This exercises real conversion, unlike the compositor's mocked subprocess tests.
