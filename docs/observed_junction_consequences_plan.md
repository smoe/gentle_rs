# Observed Splice-Junction Assessment And Conditional Consequences

Status: J0 recorded; J1 authorized for implementation, not yet verified.
J2/J3 are proposals only. Biological acceptance choices below remain with the
owner/Glen. Keep implemented software, synthetic regression, real RNA evidence
and unproven biology separate.

## Baseline And Scope

Glen's 2026-10-28 brief and Claude's review authorize this plan (J0), followed
by `AssessObservedSpliceJunctions` (J1). The review expected local main
`3b6eb465`; actual local main is `c768d3c1` and this branch starts at
`fd06bd26`, with the separately committed S2.1 close-out. `git rebase main`
reports up to date. No fetch/push/tag/CI dispatch. Untracked `output/` is
unrelated and remains untouched.

S2.1's earlier file reservation is resolved by its [close-out receipt](observed_transcript_structures_s2_plan.md#s21c-execution-and-close-out-2026-10-29).
Nevertheless J1 does not edit its publisher, GTF grammar, annotation normalizer
or annotation receipt. S0/S1/Phase 2 published shapes, bytes and biological
semantics remain frozen; new findings live in a separate report/operation.

| Brief request | Re-verified existing substrate |
| --- | --- |
| Registered observed-structure service | S1 summarize/inspect/cancel operations, shared shell and generic MCP; preserve-content lifecycle, write approval and no implicit save. |
| Explicit reference gene and source snapshots | Mandatory source-specific stable gene ID, separate annotation releases/hashes, zero-based half-open models. |
| Annotation preparation | S2.1 GTF preparation and close-out; empty attributes remain outside its strict grammar. |
| Exact comparison | S0 retains exact junction membership, model-chain relations and known junctions in a new combination; no tolerance/snapping. |
| Original-record support | Full `(source_sha256,record_index,read_id)` identities, legal hypotheses, definite/possible chains; no cross-record skip. |
| Ends versus chains | Saved `EndConstraint`, clipping and per-source model relations; no natural RNA end. |
| Portable saved artifacts | S1 JSON plus quoted `row_json` tables, source bindings and full EOF sidecar verification. |
| Preservation | S2.0a byte/work guards remain unchanged. |
| SAM | Private decoded/text seams exist, including inherited S2.2a; no public verified SAM ingress. J1 adds none. |
| Motifs/sample hierarchy/CDS | Not present in S1. J1 adds explicit reference/sample bindings; CDS/translation/topology wait for J2/J3. |

## J1 Contract To Implement

One engine-owned assessment consumes a saved S1 report; it does not rerun S1,
an aligner, original FASTQ/PAF or annotation preparation. Adapters only dispatch.
Operations: `AssessObservedSpliceJunctions`, `InspectObservedSpliceJunctions`,
`CancelObservedSpliceJunctions`. Shared shell:

```text
rna-reads splice-junctions-assess REQUEST_JSON_OR_@FILE
rna-reads splice-junctions-show REPORT [--table TABLE.tsv] [--rows N] [--binding JSON_OR_@FILE] [--verify-artifacts]
rna-reads splice-junctions-cancel REQUEST_JSON_OR_@FILE
```

Strict request/report/result-binding schemas use
`gentle.observed_splice_junction_assessment_{request,report,result_binding}.v1`.
The typed assessment artifact and producer-declared sample manifest have their
own v1 schema IDs. No defaults supply biology, source bindings or resource limits.

### Request And Verification

| Field | Contract |
| --- | --- |
| `source` | Tagged enum, initially `s1_report`, carrying the exact S1 result binding and raw report path/SHA-256/bytes. Verify raw bytes, semantic identity and every S1 artifact to EOF through existing inspection. A non-completed source yields unavailable with absent counts. Future tagged SAM sources require separate authorization. |
| `reference_sequence` | Optional exact raw FASTA binding, `format=fasta_plain_v1`. Full assembly must hash to S1's assembly FASTA SHA-256 at EOF; authority is verified only then. Declared excerpt additionally binds contig, zero-based offset and length; its assembly relationship remains declared. Absence means motif not_assessable. |
| `sample_manifest` | Optional raw binding of `gentle.observed_sample_manifest.v1`. Exact context `(sample_id,barcode,priming,file_id)` maps to library/prep/biological-sample IDs and optional cell line. Context mapping is unambiguous, library-to-prep and prep-to-sample are functions. Contradictions are binding errors; unused rows remain visible. |
| `near_window_bp` | Required 1..50 engineering admission window; no default. Owner confirmed 5 for the synthetic acceptance fixture only. |
| `recurrence_policy` | Optional explicit unit, min_units >=2, definite or definite_or_possible certainty, and require_qualified_flanks=true. Without a policy there is no recurrence label. |
| `budgets` | Positive explicit input/reference/manifest bytes, distinct junctions, annotated junction comparisons, near hits per junction, total work and output bytes. Overflow/exhaustion fails, never truncates. |
| `output_dir`, `cancel_path` | Fresh output; optional exact-request marker. Writes require approval even for refusal/cancellation. |

Hash the raw reference/manifest bytes, not normalized FASTA lines. Stream FASTA
to EOF, retaining only requested two-base windows; do not load an assembly into
memory. Validate unique headers, sequence grammar, declared excerpt extent,
target names/lengths and coordinate bounds. Soft masking is uppercased only for
motif classification. Missing sequence windows, N/ambiguous bases, excerpt edges
and unavailable orientation are explicit not-assessable states. Invalid hashes,
contradictory bindings, malformed data and resource violations are errors, not
empty or unannotated biological results.

### Saved-S1 Limits That Must Not Be Filled By Guessing

S1 stores model membership/relations, not the complete reference exon catalog.
For new near/event comparisons, read only the original hash-bound normalized
annotation files named by S1. There is no replacement-annotation argument.
Missing catalogs yield source-specific not_assessable; changed bytes or reference
contradictions fail binding. Saved S1 inspection itself is unchanged.

S1 read rows retain query contexts, but some excluded-only records have no query
context. Reconcile original record IDs against saved read rows and unavailable
rows. Retain unknown contexts as unassessed; never attach them to a library by
guessing from a name or a manifest-only assertion. Denominators disclose their
coverage and missing-context population; no extrapolation to absent units.

### Per Physical Junction

- Identity uses the full assembly/FASTA binding, contig and exact skipped
  `[start,end)`, with a full SHA-256. A short display prefix is never a key.
- Retain canonical intron bounds plus explicit 1-based intron first/last fields.
  Donor-exon-last and acceptor-exon-first fields follow transcription orientation
  only when assigned. Minus: donor is `end+1`, acceptor is `start`; plus: donor
  is `start`, acceptor is `end+1`. Gene-model disagreement/unresolved sources
  leave orientation unassigned; query alignment direction never selects it.
- Retain raw reference-forward first/last dinucleotides, oriented motif and
  opposite-orientation motif: gt_ag, gc_ag, at_ac, other or not_assessable.
  Reverse-complement and swap on minus; motifs never assign strand or filter.
- Each source lists all near annotated pairs within the window at both ends,
  signed `(annotated_start-observed_start, annotated_end-observed_end)` offsets
  and model IDs in deterministic order. Exact means both offsets zero. If
  only one end is shared among near hits, use shared_start_only/shared_end_only;
  if both ends are shared by different pairs without an exact pair, retain the
  full list and use near, not a fabricated exact match. Otherwise near,
  none_within_window or not_assessable. Donor/acceptor end names need orientation.
- Per-source/per-model labels are nonexclusive geometry: alternative_donor,
  alternative_acceptor, exon_skip, novel_combination, alternative_first_junction,
  alternative_terminal_acceptor, unannotated_both_ends. No winner. Exon skip
  requires same-model boundaries enclosing at least one whole intermediate exon;
  novel combination requires gene-wide annotated ends but no model with the pair.
- Count an original S1 record once per physical junction at its highest S1
  chain certainty. Do not union incompatible hypotheses or create a skip across
  observations. Preserve possible support under incomplete searches, never
  definite support. Qualification partitions records into qualified_only,
  insufficient_only, unassessed_only and mixed. Recurrence is conservatively
  based on qualified_only, not mixed evidence.
- Report support by library, cDNA prep, biological sample and cell line; these
  overlapping grouping axes are not additive. Unknown contexts stay explicit.
- Per-library denominators retain total cohort records, gene-coherent records
  (at least one non-OutsideModel observation in the reference source), excluded
  and unassessed records. Gene coherence and assessment completeness are separate
  axes, not an invented partition. Zero gene-coherent records is valid non-detection.
- Only an explicit recurrence policy can produce recurrent_candidate,
  not_recurrent_under_policy or not_assessable. Known independent qualified
  units must meet the threshold. Missing context/incomplete requisite evidence
  cannot prove a negative; no label asserts biological novelty or independence
  beyond the producer declaration.

### Chains, Output And Lifecycle

Keep saved S1 model relations alongside all end constraints/clipping and an
explicit unverified RNA-end status. No complete chain becomes a complete RNA.
Publish report.json last, with hashed artifacts `assessment.json`,
`junction_assessment.tsv`, `chain_assessment.tsv`, `sample_support.tsv` and
`unavailable.tsv`. TSVs use S1's one quoted row_json cell convention; all derive
from one typed result. No read sequences, protein fields or annotation replacement.
Use bounded staging, no-clobber publication, identity-safe nonrecursive cleanup
and explicit orphan disclosure, without modifying the old S1 publisher.

Reuse S1's preserve-content/no-save/checkpoint/approval policies and managed
artifact-only publication. Inspection/control use their own lanes. Raw apply
may retain the existing operation audit; inspection, refused writes and managed
retention do not alter project bytes, journal or undo/redo unexpectedly.
Generic MCP retains confirm:true; JS/Lua/Python reach the same shared operation.

Nonclaims: not in these bound references (never biological novel/absent), no
protein/translation/expression claim, declared sample hierarchy, explicit motif
authority, record counts not molecules and payload ends not RNA ends.

## J1 Verification And Commit Sequence

1. Commit J0 separately before code.
2. Add strict DTOs and fail-first synthetic tests; record actual failures, not
   zero-test filters. Cover exact/shifted boundaries, motifs/orientation,
   record deduplication, no cross-record skips, every event label, partial/full
   chains versus RNA ends, separate sources and unavailable/bad bindings.
3. Implement the engine and thin adapters/lifecycle. Test denominators, zero
   gene-coherent cohorts, qualification, incomplete search, hierarchy conflicts,
   missing contexts, budgets at/over each boundary, cancellation, no-overwrite,
   report-last/orphan behavior and unchanged project/history. No cryptic-splicing
   extraction is needed for J1.
4. Commit docs and canonical projections. Preserve S0/S1/Phase 2 golden bytes.

All fixtures are inline synthetic code using existing file-backed S1/Phase 2
helpers. Pin an offline-precomputed short-ID collision with its hash recipe;
never search for it in the runtime test. A synthetic CD274-like contig uses
LCG recipe recorded below; planted motifs and declared geometry are fixture-only.
Historical certainty/qualification and scientific acceptance stay gated below.

Derived assessment rows are byte-identical under reordered in-memory evidence
with bindings held fixed. Reordered raw files necessarily have different raw
hashes/request provenance; do not weaken those hashes to claim whole-report
byte equality. LF/CRLF FASTA has different raw hashes but the same motif bases.

Run locked/offline, no incremental, jobs=1: the J1 filter; existing
observed_transcript_structures, observed_transcripts, primer_payload and
lifecycle_slice3 lib filters; default capability_registry_parity and
mcp_capability_surface; headless lib/CLI/MCP check; default cargo check --tests.
Also scoped rustfmt, whitespace, canonical parity generate/check, real LF/CRLF
parity checkout, glossary/documentation guards and session-close against this
plan. Report source SHA/diff, platform/toolchain, Cargo.lock SHA-256 and each
pass/fail/ignored count without adding overlapping filters. Separate primary
failures from cascades. Native Windows/Linux, full suite, actual annotation/J1
CLI/MCP processes, real CD274 data and scientific acceptance remain separate.

## J2 Proposal: Conditional Coding Consequences

Not authorized for code. Consume a verified J1 report, one explicitly selected
versioned model per source, new hash-bound CDS intervals/frame/start/stop
preparation, an explicit genetic code and J1 reference bytes. Existing S2.1
annotation receipts/PrimerPayloadTranscript cannot carry these fields.

Follow model exons only up to the first observed donor, then supported blocks.
Rejoin a model only where the observed acceptor lies inside its exon; otherwise
continue with observed bases only. Report domain observed_supported,
model_rejoined or requires_unobserved_sequence; a C-terminus outside evidence
leaves CDS delta, frame modulo 3, first altered residue/new stop, conditional
length and divergent tail not_assessable. Reuse AMINO_ACIDS.codon2aa. Any future
cryptic-splicing helper extraction first needs a byte guard on its existing
CrypticSplicingCdsConsequence output. Test preserved frame, frameshift/early stop
and two-event frame restoration. No NMD assertion.

## J3 Proposal: Portable Protein Topology

Not authorized for code. Bind raw Swiss-Prot text by path/SHA-256, accession and
sequence version; reuse gentle-engine's parser, not project-bound operations.
Require exact equality between UniProt sequence and the bound reference CDS
translation; otherwise not_assessable/uniprot_sequence_mismatch. A reference
residue is encoded only if its same genomic codon/frame occurs before the
projected stop, including in-frame skips and restored frames.

Keep two axes with reason codes:

- membrane_anchor_state: not_assessable, membrane_anchor_preserved (all TM
  residues), membrane_anchor_altered (some), membrane_anchor_not_encoded_by_projected_orf (none).
- soluble_form_compatibility: not_assessable, compatible_with_soluble_form or
  not_compatible. Compatibility requires no encoded anchor, an ORF complete in
  evidence, single-pass topology, fully encoded SIGNAL and an explicit minimum
  encoded extracellular TOPO_DOM residue count. Multi-pass is not_assessable in v1.

Missing RNA end/CDS/phase/reference/UniProt/TM evidence is not_assessable, never
a negative. Every compatibility result says: RNA/sequence compatibility only;
no evidence of translation, secretion, stability or soluble protein.

## Owner And Glen Decisions Still Open

| Decision | Recommendation/status, not selected by Codex |
| --- | --- |
| Literal 5466827 -> 5468193 convention | Resolved by the owner on 2026-10-29: none of the original A/B/C interpretations describes this legacy token. Its producer displayed start+1 -> end+1. Store **[5466826,5468192)**, 1366 bp; label exon boundaries **5466826 -> 5468193**; inclusive intron is 5466827-5468192. Keep the old token only as a marked legacy alias, not an input convention. |
| Meaning of 3/1/3/0 | Owner confirmed order: sample 1 original cDNA=3, sample 2 original cDNA=1, sample 1 second cDNA/sequencing preparation=3, sample 2 second preparation=0. These are two biological samples, not four. Cell-line names, historical distinct-record certainty and flank qualification remain unspecified; do not relabel real counts. |
| Acceptance recurrence | Owner confirmed explicit library/min_units=2/definite/qualified flanks for the **synthetic fixture only**. This does not establish independence of real libraries. |
| Acceptance near window | Owner confirmed explicit 5 for the synthetic fixture; no service default is introduced. |
| Real CD274 ingress | Phase 2 -> S1 is available. BAM needs separately authorized verified SAM ingress plus external conversion; inherited private decoder is not a service. |
| Real reference authority | Full assembly needs explicit large byte/work budgets; excerpt has declared authority only. |
| J2 CDS format | Decide separate binding versus new S2.1 output; existing receipts stay unchanged. Decide whether to expose stop-to-last-junction distance as geometry only. |
| J3 | Decide minimum extracellular residues, SIGNAL requirement, single-pass-only v1 and portable UniProt format. |

A generic synthetic fixture uses two independently bound
CD274-like sources that lack the chosen junction, including shared-last-donor
geometry, four declared libraries with gene-coherent denominators above support,
verified gt_ag, explicit recurrence policy and incomplete-search 0 definite / 3,1,3,0
possible support. No human/GRCh38/CD274/plus/soluble default enters the service.
The fixture declares the confirmed sample/preparation order with no cell-line
names; its distinct qualified records are deliberately synthetic, not a
reclassification of the historical counts. Real scientific acceptance remains open.

Synthetic FASTA recipe: u64 LCG seed `0x4a31`, multiplier
`6364136223846793005`, increment 1, wrapping arithmetic, bases from the high two
bits in A/C/G/T order, contig length 5,468,300; plant GT at offset 5,466,826 and
AG before 5,468,192. Wrap lines at 80 bases. Assembly is
`synthetic_cd274_like_v1`, not GRCh38. All inputs are recreated in tests.
Pinned 8-hex display collision: a bounded offline search over starts 1..300000,
end=start+17, canonical JSON `(domain,reference,contig,interval)` with domain
`gentle.observed_splice_junction.v1`, assembly `synthetic_collision_v1`, FASTA
hash `sha256:` plus 64 `a`s and contig `synthetic_contig`, finds starts 57505 and
62716 sharing `69d1a3a7`. The test pins both full hashes; full identities differ.

## J1 Execution Ledger

The original J1 checkpoints below are historical. The separately authorized
[completeness compatibility follow-up](#completeness-compatibility-follow-up)
addresses their producer admission limitation without changing S1's v1 records.

J0 is committed separately as `f0c8bdb7`. The owner then clarified the legacy
coordinate convention above, then confirmed the sample/preparation order and
synthetic recurrence/window settings. Historical certainty/qualification and
J2/J3 decisions remain open. Code contains no CD274 default.

At `f0c8bdb7` plus the initial DTO/test skeleton, the fail-first command
`env CARGO_INCREMENTAL=0 CARGO_BUILD_JOBS=1 cargo test -q --locked --offline --no-default-features --lib junction_assessment -- --test-threads=1`
exits 101 with four E0599 errors for absent assess/binding/inspection methods.
This is an actual compile-red result, not executed tests or a zero-match pass.
The expanded tests additionally specify near/event/motif helpers and file-bound
reference/manifest/cancel/error behavior before their implementation.
The unchanged Cargo.lock SHA-256 is
`a393a1f3803b374bbaf5ef4f3e23813e974b1676f61a669dacb24063c059fb6a`.
At that fail-first checkpoint no successful J1 runtime verification was claimed.

The first behavioral run on `515c721d` plus engine/adapter/test edits ran 22
tests: 20 passed, 2 failed, 0 ignored. Primary failures: the complete synthetic
fixture was refused by frozen Phase 2's nonempty
`competitor_search_limitations` admission; the MCP response lacked the new
compact binding envelope. No cascade failures occurred. The envelope is added.
The frozen producer issue is **not changed**: `summarize.rs` refuses empty
limitations, whereas S1 `observed.rs` requires exactly empty limitations for
`search_complete`. Thus current end-to-end Phase 2 -> S1 sources remain
possible-only even when the producer intended an exhaustive search.
The saved-source CD274-like regression asserts that conservative behavior;
the definite 3/1/3/0 recurrence case is a separately labelled synthetic
in-memory derivation, never a rehashed/promoted source report. Completing the
requested full-producer definite acceptance requires a separately reviewed
compatibility fix; do not silently weaken frozen semantics here.

### Final Local Verification (2026-10-29)

Runtime/adapter implementation is `1c31d371`; closing publication ownership and
typed discovery-fixture repair is
`ac85d1e4d6c03feaf2ccff2ac407d2eea4645161`. J0 is `f0c8bdb7`,
DTO/fail-first specification is `515c721d`. The final runs below use `ac85d1e4`
plus the pending documented discovery/manual projections. The SHA-256 of
`git diff --binary -- src crates Cargo.lock assets` at that checkpoint is
`a83a29f11405a76a81eefbad41904748fcb41d3e829a7d947676fe38593fa3ee`.
All dirty source/crate/asset entries were ours; no unrelated source change was
included. `output/` was pre-existing, untracked and left untouched.

Every Cargo result below uses the unchanged Cargo.lock SHA-256
`a393a1f3803b374bbaf5ef4f3e23813e974b1676f61a669dacb24063c059fb6a`.
Platform: macOS 27.0.1 build 26A434, aarch64-apple-darwin;
rustc 1.100.0-beta.7 (`fd8e211ae301c84a86d13cf52a6d3ce08c589e2e`,
2026-10-17), LLVM 23.1.4. The initial local `git rebase main` reported up to
date on local main `c768d3c1`; no fetch/push/tag/CI dispatch was performed.

Command prefix for every Cargo row is
`env CARGO_INCREMENTAL=0 CARGO_BUILD_JOBS=1 cargo`;
every row uses `--locked --offline`. Filters overlap: do not sum their counts.

| Command after prefix | Passed / failed / ignored, or check outcome |
| --- | --- |
| `test -q --locked --offline --no-default-features --lib junction_assessment -- --test-threads=1` | 27 / 0 / 0 |
| `test -q --locked --offline --no-default-features --lib observed_transcript_structures -- --test-threads=1` | 70 / 0 / 0 |
| `test -q --locked --offline --no-default-features --lib observed_transcripts -- --test-threads=1` | 16 / 0 / 0 |
| `test -q --locked --offline --no-default-features --lib primer_payload -- --test-threads=1` | 101 / 0 / 0 |
| `test -q --locked --offline --no-default-features --lib lifecycle_slice3 -- --test-threads=1` | 11 / 0 / 0 |
| `test -q --locked --offline --no-default-features --lib engine_shell::tests::glossary_ -- --test-threads=1` | 12 / 0 / 0 |
| `test -q --locked --offline --test capability_registry_parity -- --test-threads=1` | 21 / 0 / 0 |
| `test -q --locked --offline --test mcp_capability_surface -- --test-threads=1` | 6 / 0 / 0 |
| `check -q --locked --offline --no-default-features --lib --bin gentle_cli --bin gentle_mcp` | exit 0; no test count |
| `check -q --locked --offline --tests` | exit 0; no test count |
| `run -q --locked --offline --bin gentle_examples_docs -- parity-matrix-generate` | canonical generation exit 0 (before the closing publication-only fix) |
| `run -q --locked --offline --bin gentle_examples_docs -- parity-matrix-check` | final freshness exit 0 |

Additional deterministic receipts, repeated after the closing code checks:

- `python3 -m unittest scripts.test_observed_transcript_contract -v`: 8/0/0.
- `python3 -m unittest scripts.test_capability_maturity -q`: 12/0/0.
- `python3 -m unittest scripts.test_tutorial_checkouts.DocumentationStatusTests -q`: 18/0/0, including LF/CRLF manual/catalog replay.
- `python3 -m unittest scripts.test_tutorial_checkouts.TutorialCheckoutTests.test_generated_parity_matrix_stays_byte_exact_in_both_checkout_modes -v`: 1/0/0, including the bad-policy control and committed LF/CRLF replay.
- `python3 scripts/check_documentation_status.py`, `python3 scripts/manual_docs.py` and `python3 scripts/generate_capability_maturity.py --check`: exit 0. The manual checker has no `--check` flag; an initial unsupported invocation was corrected to its default check mode.
- `git diff --check`: exit 0. Scoped `rustfmt --edition 2024 --config skip_children=true --check` passes for the new/changed engine, DTO, adapter and new test files. Three legacy touched files (`app.rs`, `main_area_dna/tests.rs`, `engine_shell/tests.rs`) retain unrelated existing formatter differences; formatting baseline copies from `c768d3c1` reproduces those same requested changes. No workspace formatting churn was included.
- `python3 scripts/maintenance_chore.py session-close --plan docs/observed_junction_consequences_plan.md`: 4 ok, 2 warnings, 0 failures before handoff. Warnings are intentional pending docs/generated files and manual plan-fidelity review; final worktree hygiene is checked after committing.

Primary intermediate failures were repaired, not hidden: compile-red missing
APIs/implementation type errors; incomplete synthetic producer admission;
missing MCP compact envelope; synthetic crosswalk IDs inherited from a different
toy reference (22/24 passed at that checkpoint); typed glossary smoke fixtures
that still supplied `{}` (10/12 passed); and a direct manual-index edit caught
by its canonical freshness guard (17/18 passed). The latter now derives from
`manuals.json` and canonical generation. No mutex-poisoning/cascade failures
were observed. Headless unused-code and macOS large-unwind-section linker
warnings remain warnings, not a claim of `-D warnings` acceptance.

The closing read-through also found an unknown-ownership cleanup gap in the new
publisher. `ac85d1e4` never removes an output reservation when identity capture
fails, discloses the orphan, and tests preservation of an empty replacement.
This extra small fix stays inside J1 publication safety, not S2.1 or a new feature.

Changed-file groups: `engine/protocol/observed_junctions.rs`, private
`primer_payload/junction_assessment{.rs,/}` and inline fixtures; shared operation,
shell/parser, lifecycle, managed execution, agent feedback, MCP and protocol-name
registry wiring; required optional-result constructor additions in RNA/GUI tests;
CLI/protocol/user/plan/changelog/roadmap docs and glossary. Canonical derived
changes are the parity matrix, manual index/navigation/source inventory and
four maturity projections (only the glossary source hash changes). No committed
byte fixture, dependency, biological computation outside J1 or reserved S2.1
file is changed.

**Still unverified:** native Windows/Linux, full suite, actual CLI/MCP subprocess
execution, feature-gated JS/Lua VM execution, real CD274/reference/library
evidence and scientific acceptance. JS/Lua regression tests are present but
neither no-default nor desktop-default builds enable `script-interfaces`; do not
claim those tests ran or compiled. Existing generic engine routes remain shared.
The requested full-producer definite recurrence acceptance is specifically
**not closed**, for the frozen Phase 2/S1 completeness contradiction above.
Historical 3/1/3/0 certainty/flank qualification, reference authority and J2/J3
choices still need explicit evidence/owner decisions. No J2/J3 or SAM service
implementation, protein/expression/novelty claim or downloaded reference is added.

Closing hygiene caught the newly added private module declaration out of
rustfmt's lexical order in `primer_payload.rs`. Only that new declaration is
reordered; the baseline formatting exceptions above are left alone. The runtime
counts/source overlay above predate this declaration-order-only adjustment and
the documentation commit `4fe68c42`; they are not mislabelled as new runtime
runs. The declaration and all new modules then pass scoped formatting, and a
locked/offline default Cargo check is run before final handoff. Final
session-close has 4 ok/2 warnings/0 failures: only the preserved untracked
`output/` and the manual plan-fidelity reminder remain.

### Final Rebase Receipt

Local main advanced during verification to
`e39469bdbc64b3e03ea540cd9954bc40c9744243` (documentation/tutorial-history
changes, not Rust engine changes). The closing `git rebase main` succeeds
without conflicts. Original verification SHAs above are retained as historical
receipts, not relabelled as fresh runs. Current equivalents are J0 `83ebee85`,
DTO/tests `ba4bb4b1`, engine/adapters `f5d5b7c5`, publication/discovery fix
`4044d322`, docs/projections `eb6ae3e8` and module ordering `656ce1ee`.
`git diff --quiet 803f32a5 HEAD -- src crates Cargo.lock assets` plus the three
ClawBio maturity projections exits 0: the final rebased Rust/lock/asset trees
are byte-identical to their verified pre-rebase counterparts. Main is an
ancestor of HEAD; only the original untracked `output/` remains. The final
locked/offline default `cargo check -q` exited 0 before this documentation-only
rebase, after the module-order adjustment. No new runtime/platform/scientific
acceptance is inferred from rebasing.

## Completeness Compatibility Follow-Up

Owner-authorized separately on 2026-10-29, starting at `186937bd`: allow Phase 2
to admit `declared_complete=true` with exactly empty
`competitor_search_limitations`. Retain rejection of whitespace-only text and
incomplete+empty declarations. Preserve every prior nonempty limitations string
and its conservative S1 interpretation, including `none`; no normalization or
natural-language interpretation is added. All other producer/provenance checks,
the report schemas/fields/digest recipes and S1 certainty logic stay unchanged.
Previously saved reports are not rehashed, upgraded or interpreted differently.

The existing synthetic CD274-like producer now supplies the actual declaration
and passes saved Phase 2 evidence through S1 into J1. Definite 3/1/3/0 support
and policy-bound recurrence must come from that real service path, not edited
in-memory certainty flags. Countercases retain possible-only support for limited
or incomplete producer declarations. These are still synthetic software tests;
historical library certainty/qualification, real exhaustive search and biology
are not inferred. J2/J3 and SAM ingress remain outside this fix.

Verification on macOS arm64, Rust `1.100.0-beta.7`, 2026-10-29:

- At `186937bd926449930b69a5d08d7bb33fcf5174d6` plus source overlay
  `e2f4d23b37f3d7d13cd884195c7f419e61c3ee7095c698cfc62e711bbc1fa934`
  (`git diff --binary -- src crates Cargo.lock assets | shasum -a 256`),
  locked/offline headless library filters pass: `primer_payload` 103/103 and
  `observed_transcript_structures` 70/70. Scopes overlap; they are not a full-suite
  total. Existing S1 and Phase 2 preservation oracles pass unchanged, including
  the saved producer-to-S1-to-J1 synthetic complete/limited/incomplete controls.
- These runs used `CARGO_PROFILE_TEST_DEBUG=0 CARGO_INCREMENTAL=0
  CARGO_BUILD_JOBS=1` to limit test-build overhead. The earlier normal-debug
  pre-fix build was interrupted during compilation, before tests ran; it is not
  recorded as a red or green test run. Existing headless unused-code and macOS
  large-unwind-section linker warnings remain; no `-D warnings` claim is made.
- Rebased onto local main `1fdeb1d97d15498630b88da72666fdd1eaaf8f6c`,
  then onto its newer tip `1c98980cdda4c9f841593c2c0f52e35f492a9619`.
  Retained both changelog additions when restoring the first uncommitted patch;
  the second restore merged without conflicts.
  The two changed Rust files are byte-identical to the tested pre-rebase files;
  other upstream engine changes mean the preceding runtime results must not be
  represented as post-rebase runs. S1 logic, frozen oracles and `Cargo.lock`
  remain unchanged; lockfile SHA-256 remains
  `a393a1f3803b374bbaf5ef4f3e23813e974b1676f61a669dacb24063c059fb6a`.
- At the final rebased source, 26 focused Python/source and documentation-checkout
  tests pass
  (8 contract + 18 documentation-status, with LF/CRLF replay). Documentation
  status, manual and capability-maturity freshness, scoped Rust formatting and
  whitespace checks pass. Default-feature `cargo check -q --locked --offline`
  with `CARGO_INCREMENTAL=0 CARGO_BUILD_JOBS=1` exits 0 on both rebased bases.
  Session-close has 4 ok/2 warnings/0 failures: the intentional uncommitted fix
  and preserved untracked `output/`, plus the manual plan-fidelity reminder.
  The diff remains the authorized completeness compatibility fix, its synthetic
  regressions and directly related contracts/receipts, not S2.2 or a new feature.

Still unverified: full suite, runtime rerun on the rebased engine, native
Windows/Linux, actual CLI/MCP subprocess and script-VM execution, real CD274
reference/library evidence and scientific acceptance. No new byte fixtures,
dependencies, generated artifacts, output-format changes or scientific claims.
