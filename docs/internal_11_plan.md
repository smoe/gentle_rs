# Internal .11 Workflow Completion And Acceptance

Reconciled: 2026-09-20. The published `.10`
tag is `84f34a9e`; its [acceptance ledger](release_notes/release_notes_v0.1.0-internal.10.md#exact-candidate-gate-ledger)
remains pending. `codex/internal.11` was merged at `8597a6ee`; it no longer
needs a separate integration decision. The
[gene-informed primer-pair GUI](gene_assay_study_gui_plan.md), integrated from
`gentle_rs_2_main` at `a1f5e61e`, remains the primary `.11` scientific workflow;
package parity and exact-revision acceptance support that aim. The capture
programme below remains separate and deferred. No version bump, tag change,
upstream push, workflow dispatch or publication is part of this follow-up without
explicit approval. The owner requested the `.11` version/notes update on
2026-09-20; its final candidate SHA remains to be selected. The earlier review
started at `36ba15a0` on `smoe/gentle_rs`, not the fork's stale `main`.

## Review Reconciliation

Claude's supplied review used `85bd9fdb`. The original Codex plan below treated
capture GUI, specificity and order handoff as successive implementation slices.
Claude correctly prioritized release-status repair, TSS correctness, packaged
runtime parity and external acceptance, and identified missing lab inputs for
capture specificity and order readiness. This revision incorporates that
feedback after checking the newer tree; no fresh Claude CLI review is claimed.
His subsequent review at `b16f0453` found that T3 still lacked ID-to-label
resolution and case-independent grouping. The bounded follow-up below addresses
those omissions; the earlier reconciliation did not establish their completion.

| Review item | Current evidence and remaining action |
| --- | --- |
| R0: release status | `.10` is published without recorded exact-candidate acceptance; preserve its Pending ledger. Later fixes need their own SHA. Manual gel commits are in the tag, not an unmerged `.11` branch. |
| T1/T2: endpoints and exclusions | `085b98ec` shares `transcript_five_prime_endpoint` between extraction and inventory and adds locus-level `unassigned_transcripts`. Existing tests cover both reverse-join forms, 5'/3' clipping, fuzzy/mixed strands and legacy JSON. Do not reimplement. |
| T3: canonical gene links | `085b98ec` covered label-to-ID only. The follow-up adds unambiguous ID-to-label lookup before symbol filtering and groups by resolved ID (case-insensitive label fallback), retaining source/strand separation, original annotations and ambiguity warnings. Regressions cover mixed/case-only labels, ID-only distinct starts on both strands, ambiguous links and historical collection/approval safety. Glen must still check real TP73/DeltaNp73 symbol-versus-ID memberships and saved collections at the final SHA. |
| P1: desktop packages | Shared staging and extracted-package checks now cover all seven entrypoints and tracked resources on Windows/macOS/Linux. Offline tests cover layouts, archive round-trips, relocation and failures. Native CI must still run `gentle_cli capabilities`, MCP and tutorial-manifest validation from the actual packages; receipts and live GUI acceptance remain required. |
| A1: exact-SHA acceptance | Glen: TP73/DeltaNp73 inventory, materialization, reload, validated scans and native inspection, plus branch/reverse-complement, digest and bounded Simple-PCR tutorial chapters. Gene-identity grouping in `a3d74ca9` changed `tss_id`/`output_seq_id` for ID-bearing rows: preview afresh under a **new collection ID**; older previews, approvals and receipts are not comparable to the new grouping and cannot authorize it. Bind binary, project, input and receipt hashes to one new candidate; retain the old `.10` ledger separately. |

T3 focused verification after rebase on `main` at `d9f75837`:
`cargo test --lib --locked --offline -j 1` passed with
filters `tss_workspace::tests` (24), `engine_shell::tests::tss_` (4) and
`app::tests::tss_` (4), each with `-- --test-threads=1`.
`cargo check -q --locked --offline -j 1`, `cargo fmt --check` and
`git diff --check` also passed. These working-tree checks do not replace the
full-workspace or real-data/live-GUI acceptance at the final candidate SHA.
Main's new synthetic GUI tutorial oracle was explicitly refreshed from the
engine-emitted preview. Its independent workflow test, generated-tutorial check
and catalog/manifest consistency tests passed; no user collection was migrated.

### Verification Before Further Features

The second Claude review distinguishes implemented fixes from unverified builds.
The semantic-ID lifetime fix is committed as `d9f75837`; earlier pre-rebase
checks do not establish that the merged tree builds with `gui-test-support`.
Windows run `35386839336` at `caf84061` failed two byte-exact tutorial-index
checks and a path-separator assertion. The `8284c4ca` LF/path fix was verified
on Windows at `c2dead07` in [run 35419061304](https://github.com/smoe/gentle_rs/actions/runs/35419061304):
`workflow_examples` passed 71/71. That run's library suite exposed 210 further
Windows failures (3,667 passed, eight ignored). `36ba15a0` addresses those path,
fixture and export-root failures; all 210 affected cases passed locally across
the batch and direct rerun, not on native Windows. Its push CI
[35445426797](https://github.com/smoe/gentle_rs/actions/runs/35445426797) sampled
Linux (still running at the 2026-09-19 check); Windows and macOS were skipped.
Windows acceptance of the fix therefore remains pending. Preserve strict
byte/receipt checks; no scientific evidence was regenerated.

`fa137598` adds the missing label-to-multiple-IDs diagnostic without changing
grouping or output IDs. Local `cargo test --locked -j 1 --lib tss_` passed
167 tests with four explicit ignores; locked Cargo check, formatting and
whitespace checks passed. Both strands and nonmatching-query suppression are
covered. Warnings remain digest-bound, so affected previews require fresh
approval; this is not native-platform or real-data acceptance.

For Glen's annotated-export check, `36ba15a0` also fixes user-visible
`tss_profile_export` directory validation: a Windows drive/UNC prefix is joined
to its root before filesystem inspection. Check exports to fresh native Windows
destinations, including GenBank/EMBL/FASTA companions and receipts; local macOS
tests do not certify native Windows behavior.

The source-comparison work in `f7afe95d` was an explicit user-requested
exception to this verification freeze. It does not authorize further feature
scope before the exact-SHA gates. Windows run
[35461010366](https://github.com/smoe/gentle_rs/actions/runs/35461010366) at
`885fac49` produced ten failures: six causes, eight primary failures and two
PDF-triggered mutex-poison cascades. See the [root-cause analysis](testing.md#september-2026-failure-analysis).
The follow-up retains PDF write access and byte-exact LF protection, restores
the shared shell grammar with a quoted caller, and tests raw traversal through
export and receipt readers. Distinct catalog/assembly fixtures are restored;
unsupported cross-namespace comparisons are removed, not replaced with aliases.
Independent anchor-to-assembly authority remains a documented limitation.
Local regression checks do not certify native Windows behavior; run the new
Windows path probe and the full suite at the final candidate SHA.

Glen's `bot/gui-usability-20260919-885fac49` contribution restores the reviewed
native DNA-window sizing, bounded/resizable toolbar, 16 MiB TSS worker and root
poller wake-up. Its retained audit accepts the synthetic TSS lifecycle twice at
18/18 steps on tested code `677c796a`, while explicitly withholding general GUI
responsiveness acceptance. That evidence predates the combined Windows-fix and
merge result and therefore cannot replace the final exact-SHA GUI gate. The
package version at that point was still `0.1.0-internal.10`. The owner-requested
September 20 version/notes repair advances development metadata to `.11`;
the already-created tag at `e1c7dfb2` predates that repair. Select a corrected
SHA before candidate packaging; see the [current release notes](release_notes/release_notes_v0.1.0-internal.11.md#acceptance-status).

Before adding N1 or another scientific extension:

1. Freeze one full SHA and run `ci.yml` explicitly for `macos`, `linux` and
   `windows`. Each existing job checks all features, including
   `gui-test-support`; sampled push CI is not an all-platform verdict.
2. Run `container.yml` with that exact `candidate_sha` and `publish=false`.
   Verify the candidate receipt, not just the dispatch's workflow revision;
   a run labelled with a newer workflow can still build an older source.
3. Run build-only `release.yml` for the same source and retain extracted-package
   results and receipts for all three platforms. Use the version label matching
   Cargo.toml; do not move the published `.10` tag or publish new assets.
4. Bind Glen's TP73/DeltaNp73 and tutorial receipts to that SHA. The primer-pair
   study GUI additionally needs the original typed publication request and its
   digest-bound plans; neither synthetic tests nor a rendered paper substitute
   for those missing real-data inputs.

Use an explicitly approved pushed ref in `smoe/gentle_rs`, or a newly pushed
compatible fork ref only if separately approved; see the
[workflow-ref requirements](release.md#build-only-candidate-verification).
Git push access and Actions dispatch access are separate credentials. Report an
unavailable dispatch as blocked, not passed. No additional feature scope is
authorized by completion of the local regression fixes alone.

### Prepared Verification Gate (Not Dispatched)

After the diagnostic and documentation commits, freeze their full final SHA and
retain it in the handoff. Before any approved dispatch, ensure the named pushed
ref resolves to that SHA; do not use a moving `main` or the published `.10` tag.
Use `v0.1.0-internal.11` only with a candidate whose committed Cargo metadata
declares `.11`; build-only use does not change or certify an existing tag.
Run commands only after explicit approval, substituting that SHA and pushed ref:

```bash
CANDIDATE_SHA='REPLACE_WITH_FULL_FROZEN_SHA'
WORKFLOW_REF='REPLACE_WITH_PUSHED_REF_AT_THAT_SHA'
gh workflow run ci.yml -R smoe/gentle_rs --ref "$WORKFLOW_REF" -f platform=windows
gh workflow run ci.yml -R smoe/gentle_rs --ref "$WORKFLOW_REF" -f platform=macos
gh workflow run ci.yml -R smoe/gentle_rs --ref "$WORKFLOW_REF" -f platform=linux
gh workflow run container.yml -R smoe/gentle_rs --ref "$WORKFLOW_REF" \
  -f tag=v0.1.0-internal.11 -f candidate_sha="$CANDIDATE_SHA" -F publish=false
gh workflow run release.yml -R smoe/gentle_rs --ref "$WORKFLOW_REF" \
  -f tag=v0.1.0-internal.11 -f candidate_sha="$CANDIDATE_SHA" -F publish=false
```

| Gate | New candidate run ID | Acceptance evidence |
| --- | --- | --- |
| Windows / macOS / Linux CI | Not dispatched; blocked pending approval and pushed ref | Each run's `headSha` equals the frozen SHA and its conclusion is `success`; none of the three native jobs can be substituted by sampled CI. |
| Headless container | Not dispatched; blocked pending approval and pushed ref | Candidate `revision` in log/receipt equals the frozen SHA, build/smoke succeeds, mode is `validate_only`, image is not pushed. Workflow revision alone is insufficient. |
| Extracted desktop packages | Not dispatched; blocked pending approval and pushed ref | All three native extracted-package checks and aggregate validation succeed, receipts bind the frozen SHA; publication is skipped. |

Replace those placeholders with actual run IDs/URLs and outcomes only after
execution. Retain repository/run URLs alongside receipts (the current JSON has
no repository field). Blocked, cancelled, skipped and unavailable are not passes.

Optional follow-ups, not prerequisites silently added to this patch:

- T4: restriction-scan, digest and primer-specificity already accept
  `--tss-collection`; TFBS did before. Consider construct-reasoning separately,
  with parser exclusivity and parity tests. Do not mechanically enable physical
  pool operations on TSS windows.
- N1: migrate only TSS preview/materialize/open workers to the existing managed
  service, with held-worker cancellation/status tests and a representative
  synthetic journal-copy benchmark. General scheduling/dependencies, durable
  recovery and the other GUI workers remain deferred.
- C1: explicit CDS/isoform-retention requirements with unknown boundaries
  reported as unassessed, never pass. This does not establish specificity.
- G1: real-image memory/repaint acceptance and verified ladders for the shipped
  manual gel editor. Automatic peak detection remains deferred.
- Discuss deprecated ClawBio normalizer removal; "earliest .11" is not consent.

Capture single-primer reference assessment, Primer3 `ntthal`, order-ready state
and procurement handoff are deferred until the retained order sheet, exact
protocol oligos and co-presence stages exist. Only pair specificity currently
exists. Calibrated cofactor prioritization also remains deferred until the
experimental-anchor contract, expression binding, explicit background universe
and enrichment exist; private CUT&RUN inputs must not enter public fixtures.

## Integrated Work And Boundaries

The tutorial/PDF integration described here is merged; see verification below.

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

The starting branch retained `ea6f9009`, the Windows ClawBio fixture fix. It is
now integrated, but remains outside the published `.10` tag.

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

## Deferred Capture Design

Retain these design constraints for later review, not as committed `.11`
deliverables. The lab-input gate above precedes specificity/order implementation;
the eight-gene replay cannot substitute guessed oligos or biological targets.

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
