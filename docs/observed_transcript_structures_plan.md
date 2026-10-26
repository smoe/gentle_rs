# Observed Transcript Structures: Evidence-Layer Proposal

Status: S0 technically authorized by the owner on 2026-10-27 and implemented
as a pure Rust kernel; verification is recorded in the ledger below. S1 and
later slices remain proposals, not implemented or authorized by this S0 task.
Claude reviewed the original draft read-only at `ce4d33b5`; the reconciliation
and implementation are Codex's work, not a fresh Claude/native/scientific verdict.
The [roadmap](roadmap.md) owns priorities; this is not another release gate.

## Aim And Recommendation

Support read-structure interpretation for an explicitly selected reference gene
as a generic, downstream evidence layer:
**which splice patterns were observed, which supplied transcript models remain
compatible, and what is still unknown?** Reuse validated alignments; do not build
an aligner, replace primer tracing, or turn an RNA observation into a protein.
The reference gene is a request variable, not a TP53 default. Glen's TP53 audit
is one retained replay example, not the operation's target or naming convention.

Adopt the audit's separation of observed structures, reference aliases and
conditional consequences. Improve it by distinguishing exact splice chains
from tolerant matches, chain completeness from RNA-end completeness, and read
support from independent molecules. A chain is not a complete transcript.

Suggested operation name: `SummarizeObservedTranscriptStructures`. All new
operation names and file schemas below remain proposals, not available commands.

## Baseline And Archive Inspection

Inspected checkout: `gentle_rs_2_main`,
`ce4d33b5b990701f2ce51c85f6ce2821df9dc8ce`; local `main` and fetched
`origin/main` matched at inspection. No new fetch or remote-tip claim.
Tracked files were clean; unrelated untracked `output/` remains untouched.

Owner-supplied archive: `tp53_isoform_audit_20261026.zip`, SHA-256
`d629bc8e50a5e2429703aaef1a24fe163f7998cce015ea0ead7ef1ed7e3b1b37`.
Read-only inspection verified safe unique member paths, 25 regular files,
ZIP CRC integrity and all 24 payload entries in `MANIFEST.sha256`. The producer's
absolute manifest prefix was explicitly mapped to the archive root, not used as
an extraction destination. This checks integrity, not authorship or biology.
No supplied script, BAM parser, aligner or consequence predictor was executed.
No archive content or real primer/read sequence is copied into this repository.

The archive contains two small BAMs with indexes, read/primer-ID lists, the two
Python scripts, reports, seven result TSVs, Ensembl 116 TP53 GTF and a UCSC hg38
`ncbiRefSeq` snapshot. Inspection of the saved tables confirms 21 read rows
(9 barcode13, 12 barcode14), nine exact-primer memberships (2 and 7), ten
catalog rows totalling 21 reads, and 147 junction rows. These are table counts,
not an independent re-derivation from the BAMs.

The report distinguishes eleven canonical-compatible full chains into ten
exact and one shifted by 5 bp, five partial known chains, two reads sharing a
new junction, one retention-like case and two unspliced cases. The script's
`spliced_known_partial_chain` class combines the first sixteen; the report's
refinement is not yet a generic machine classification. The two recurrent
candidate reads share barcode14 and primer membership: independent-library or
independent-molecule replication has not been demonstrated.

`input.sha256` names an external full-cohort selection table, genome FASTA and
full GTF, but their bytes are not all bundled. The two primer-free candidate
alignments are retained as a result table, not a complete trimmed-read/alignment
replay package. The archive is therefore sufficient for contract review and
bounded BAM interpretation later, not a self-contained replay of original
selection, primer exclusion, search completeness or protein consequences.

## What To Reuse, And What Not To Port Literally

| Audit component | Useful function | Required change for GENtle |
| --- | --- | --- |
| `scripts/analyze_tp53_read_structures.py` | CIGAR-derived blocks/skips, exact grouping, per-read/catalog/junction tables | Curated-input selection becomes explicit policy; no hard-coded TP53/chr17/minus strand; bounded shared geometry rather than another private walker. |
| `compatible_transcripts` | Consecutive junction-chain comparison | Preserve all valid matches instead of greedy `choices[0]`; include supported block geometry and source coverage; exact and tolerant relations remain distinct. |
| Exonic-overlap classification | Useful descriptive coverage evidence | The hard-coded 0.8 cutoff is not an isoform or retention proof. Reference blocks containing deletions are not equivalent to read-supported bases. |
| `NOMENCLATURE_PROPOSAL.md` | Stable observed-pattern naming, aliases separate from protein | Full digests and explicit identity scope; generic core, optional P53-family display aliases. |
| `predict_novel_tp53_coding_consequence.py` | Explicitly conditional CDS interpretation | Separate later operation using versioned CDS/phase and genome bindings; do not port hard-coded intervals or the literal HGVS-like effect as a generic prediction. |

The script filters unmapped/secondary/supplementary records but depends on its
preselected BAMs for locus, strand and reported MAPQ60 selection. Its read flag
`forward/reverse` is not itself transcription direction. It primarily compares
Ensembl models; the retained RefSeq comparison is not the same shared comparer.

Coordinate admission needs particular care: the script's chain token uses
`a+1>b+1` for a skipped reference interval `[a,b)`, while its junction TSV uses
`left_exon_end_1based=a` and `right_exon_start_1based=b+1`. These are different
left-boundary conventions, not interchangeable donor fields. Preserve original
audit IDs as archive-bound aliases; reconstruct new identities from geometry.

## Existing GENtle Seams

Use symbols, not drifting line-number anchors:

| Surface | Existing evidence and proposed use |
| --- | --- |
| `src/engine/protocol/primer_payload.rs` | `PrimerPayloadQuery` binds original file/record/anchor, excluded intervals and payload hash. `PrimerPayloadSegment.alignment_operations` retains validated operations in increasing reference order. `PrimerPayloadReadEvidence` and `PrimerPayloadAlignmentReport` retain alternatives, coverage/accounting and artifact bindings. |
| `src/engine/analysis/rna_panel/primer_payload/paf.rs` | `row`, `cigar`, `difference`, `subtract_masks`: existing validated geometry and query-support boundaries. Reuse/factor the narrow kernel without changing Phase 2 bytes or interpretation. |
| `src/engine/analysis/rna_panel/primer_payload/reference.rs` | `annotate` records only exactly annotated skips in `splice_edges`. **Do not reconstruct observations from that field alone:** novel skips survive in `alignment_operations`, not necessarily in `splice_edges`. |
| `src/engine/analysis/rna_panel/primer_payload/binding.rs` | Source/query accounting and exclusions must remain verified, including original record identity and multiple primer anchors. |
| [Payload evidence contract](primer_payload_evidence.md) | Existing Phase 2 input is bound PAF, not general SAM/BAM. No general BAM-to-observed-structure route was found in the inspected paths. An explicit ingress adapter is work, not an existing capability. |
| [Source presentation](transcript_source_presentation.md) | Existing Ensembl/RefSeq source separation and physical exon identities are reusable. Current source formats do not automatically admit this archive's GTF/UCSC JSON; normalization must be explicit and versioned. |
| `src/engine/analysis/feature_expert_ops.rs`, `build_splicing_expert_view` | Reuse rendering/projection seams later; do not insert observed reads into the reference annotation or silently change its selected transcript set. |

The existing [Phase 2 plan](primer_payload_partner_phase2_plan_20261026.md)
remains unchanged. Its pending acceptance is not inherited by this proposal.

## Review Reconciliation

Claude's review was supplied by the owner and checked against the same source
baseline. Claude did not open the ZIP or run anything. Codex rechecked the
operation walker, query geometry, annotation, DAG, published rows and
project-bound splicing builder; no runtime verification was performed.

| Finding | Original Codex plan | Revised minimal contract |
| --- | --- | --- |
| F1: strand is not an observed RNA-strand fact | Put transcript strand in the exact-chain digest | Use a strand-neutral `observed_junction_chain` key over bound genomic reference identity and ascending skipped intervals. Strand evidence/alternatives and transcript-direction display are separate. |
| F2: optimal DAG membership is not published | Bind an observation to a mapping hypothesis without specifying its reconstruction | Factor the existing pure DAG calculation into one shared helper, preserving Phase 2 bytes; expose admissible edges/paths as well as segment membership. Never combine all optimal segments into one invented path. |
| F3: Phase 2 annotation/classification is not compatibility | Reuse bound evidence without explicitly excluding broad gene labels | Derive supported blocks/skips from operations; ignore `gene_ids`, `transcript_ids`, `classification` and expected-target fractions for structural admission/compatibility. Compare supplied models afresh. |
| F4: masks need projection through operations | Reuse geometry generally | Factor/reuse the operation walker and `reference::project`; never reverse-map the union `aligned_query_intervals`. Coarse detail stays unassessed. |
| F5: several anchors can cover one record | Count original records once | Reconcile anchor hypotheses at record level; exact agreement contributes once, partial containment stays a relation, and conflicting chains remain alternatives. Payload ends never become natural RNA ends. |
| F6: partial/unspliced identity | Add bounded-fragment group IDs | Keep per-record observation IDs, but no chain-family or fragment-group identity for zero-junction reads in v1. Partial observed chains remain distinct, with span/end companions and containment relations. |
| F7: known recombinations and new junctions differ | One new-junction/combination bucket | Emit distinct relations and verify terminal supported blocks as well as skips before declaring exact model compatibility. |
| F8: annotation reuse and gene selection | Allow symbol resolution and leave normalization open | S0/S1 take stable gene IDs and the existing payload annotation geometry, with source metadata in request bindings; no third annotation model or duplicate alignments for another provider. |
| F9: Splicing Expert is project/symbol-bound | Reuse its builder/projection seams | Reuse rendering/lane primitives only; construct observed lanes from the saved report, not `splicing_group_label` or a project feature. |
| F10 and ordering | General budgets and a combined SAM/BAM slice | Name structural ceilings and state/journal tests. Prefer bounded explicit SAM-text ingress for S2; binary BAM remains a later dependency decision. |

Qualifications to the review: "no identity" for unspliced reads means no
structural family, not loss of auditable observation identity. A shared DAG
needs legal path relationships, not only its optimal-node union. SAM text is a
dependency-saving ingress option, not a trivial or already available parser.
The request's reference label is an alias, not sufficient coordinate authority
for a chain hash. The owner subsequently authorized the bounded S0 technical
implementation; this does not authorize S1 or biological calibration.

## S0 Implementation Ledger (2026-10-27)

- Plan commit: `c10bbd3e`; rebase onto local `main` was a no-op at
  `ce4d33b5b990701f2ce51c85f6ce2821df9dc8ce`. No fetch/push or claim about a
  newer remote tip. S0 is a separate implementation diff on that baseline.
- `src/observed_transcript_structures.rs` exposes the pure `Input`/`Summary`
  types and `summarize`; `comparison.rs` resolves stable IDs in the existing
  annotation shape; `kernel.rs` validates observations and reconciles records.
  `geometry.rs` and `dag.rs` are shared by the live Phase 2 path. No operation,
  command, provider API, file writer or new dependency is introduced.
- The physical key binds assembly, declared reference FASTA hash, contig and
  ascending exact skips. Observation IDs additionally bind original record,
  query/alignment receipt, scope, geometry, masks and endpoint constraints.
  Source/model metadata and query MAPQ remain inspectable, not inferred truth.
  Raw-skip support is not a primer-free or biologically validated junction claim;
  flank qualification remains a separate per-skip value.
- Legal DAG paths are bounded explicitly. Identical anchors count once;
  contained partial chains remain separate and conflicting anchors remain
  alternatives. No transition is invented between distinct alignment segments.
  Empty/coarse or transcriptome-only geometry cannot become unspliced genomic
  evidence. Zero-junction reads keep observations but no family.
- `tests.rs` contains only inline synthetic loci/records/operations, including
  source-ID ambiguity, two loci/strands, overlapping antisense, multi-gene
  membership, nested intronic reads, masks, retention-like coverage versus `D`,
  partial chains, terminal coverage, counters and exhaustive DAG/basewise oracles.
  The test-only `legacy_structure_tests.rs` retains the `ce4d33b5` analyzer to
  compare serialized Phase 2 rows and work counts over ten synthetic cases.
  This is not a frozen real-data fixture or scientific replay.
- The owner authorized reducing the checkout's `target/` and focused Rust tests.
  `cargo clean --locked --offline -p GENtle` removed 26.9 GiB, reducing the
  measured local target from about 31 GiB to 8.6 GiB before rebuilding. Retain
  dependencies, source and unrelated `output/`; no other checkout was cleaned.
- Final macOS/aarch64 headless execution passes **23/23 S0 tests** and
  **35/35 primer-payload tests**, with zero failures or ignored cases. The
  latter includes ten byte/work comparisons against the retained analyzer.
  The initial 18-test checkpoint is superseded, not inherited by edited code.
  Locked/offline headless Cargo check, workspace formatting, 11 LF/CRLF-aware
  documentation-status tests, local links and whitespace pass. Session-close
  has four OK, two expected working-tree/manual-plan warnings and no failures.
- Pre-rebase executed source was the S0 working diff over `c10bbd3e`, not bare main. Its
  11 Rust-path manifest SHA-256 is
  `373aa8731630bf32170086e7dcd2aba89ae2934085b10d0a57f258184b9515cc`.
  Recreate by hashing each file, then SHA-256 the sorted UTF-8 concatenation
  `<path> NUL <lowercase file SHA-256> LF` for `src/lib.rs`,
  `src/observed_transcript_structures.rs`, all five `.rs` files directly in
  `src/observed_transcript_structures/`, and payload
  `{evidence,paf,reference,legacy_structure_tests}.rs`. No Rust file changed
  between this fingerprint and the final checks.
- Toolchain: `rustc 1.100.0-beta.5 (8ea38001b 2026-10-13)`,
  `cargo 1.100.0-beta.5 (7552a885e 2026-09-29)`, `aarch64-apple-darwin`;
  unoptimized test profile with debuginfo and no default features. `Cargo.lock`
  is unchanged, SHA-256
  `97671e924a77d61de6dad2e6c6fe3afbee440b8606d79ca0ad2b494fd32d66ad`.
  The pre-rebase test binary at
  `target/debug/build/GENtle/901e473ba9ebfd83/out/gentle-901e473ba9ebfd83`
  had SHA-256 `4ef2004ab7f4cc3a8e0d3d8a968af5714ee54895e377a045f38ab6e9ef53896a`;
  that historical receipt does not identify a later rebuilt binary at that path.
  Final measured target size is about 14 GiB after tests and Cargo check. Existing
  headless dead-code warnings and the large `__eh_frame` linker warning remain;
  no warnings were suppressed. No desktop/default-feature, full-suite,
  Linux/Windows, real-data, native, package or scientific pass is claimed.

Executed commands (the two Rust filters use the same final test binary):

```bash
CARGO_BUILD_JOBS=1 cargo test --locked --offline --no-default-features --lib observed_transcript_structures -- --test-threads=1
CARGO_BUILD_JOBS=1 cargo test --locked --offline --no-default-features --lib primer_payload -- --test-threads=1
CARGO_BUILD_JOBS=1 cargo check -q --locked --offline --no-default-features --lib
cargo fmt --all --check
python3 scripts/check_documentation_status.py
python3 -m unittest scripts.test_tutorial_checkouts.DocumentationStatusTests -q
git diff --check
python3 scripts/maintenance_chore.py session-close --plan docs/observed_transcript_structures_plan.md
```

### Commit And Rebase Integration (2026-10-27)

S0 was committed as `eaed85ff` and rebased onto local main
`f0b8371f4a99a455c60fc71a7769f3ea0380e7f3`; its new SHA is `490014c6` and
the plan commit became `b739a8a5`. Both changelog additions were retained.
The S0/Phase 2 Rust fingerprint above is unchanged, but main changed the wider
engine/adapter code, so the pre-rebase runtime verdict is not inherited.

The first post-rebase test build fails with three `E0599` errors, before any
test runs. Main's new `src/engine/lifecycle_slice3_tests.rs` and
`src/command_execution/lifecycle_tests.rs` call nonexistent `undo()`/`redo()`.
The integration repair substitutes `undo_last_operation()` once and
`redo_last_operation()` twice, with no assertion or runtime-code change.
Post-repair macOS/aarch64 verification passes 23/23 S0 tests, 35/35 primer-payload
tests and 11/11 `lifecycle_slice3` tests, with no failed or ignored cases. The
first two Rust commands above were rerun, followed by the additional filter
below. Locked/offline headless library Cargo check, formatting, documentation
status, 11 LF/CRLF-aware documentation tests and whitespace pass. Existing
headless/linker warnings remain. Lifecycle source/runtime/VM/platform acceptance
outside those filters is not inferred from these passes.

```bash
CARGO_BUILD_JOBS=1 cargo test --locked --offline --no-default-features --lib lifecycle_slice3 -- --test-threads=1
```

The executed source is `490014c6` plus the three documented test-call repairs,
not bare `f0b8371f`. The repaired source SHA-256 values are
`c391dd40e79802f3494445ac140fc0af8239428e0f8b28500ed6a516cc70da4c`
for `src/engine/lifecycle_slice3_tests.rs` and
`655b940a733b8ecf1efa3c3354a2b1b5a6140662132e746bddfef2a9560d5ac7`
for `src/command_execution/lifecycle_tests.rs`. The S0 11-path fingerprint and
lockfile hash above are unchanged. The rebuilt test binary at the same path has
SHA-256 `3cb229254da0c4fcfbd418484e2e298fc2477d3afd650fa7e1debbc834335ba7`.
No production Rust source changed after the S0 commit; these test repairs and
the documentation reconciliation are a separate integration follow-up.

## Proposed Contract

### Reference Gene Is An Explicit Request Parameter

Owner clarification, 2026-10-26: require a `reference_gene` selector in the new
request, shared unchanged by engine, CLI, MCP/agent and later GUI consumers.
S0/S1 accept a stable gene ID in a named, hash-bound annotation source, not a
symbol: the current payload annotation has no symbol field. Record the selector
and resolved source-specific IDs in the request binding and report. A later
symbol-resolution feature needs its own explicit source contract. Do not assume
HGNC, human, GRCh38 or a chromosome. Assembly/locus must agree with the bindings;
organism and transcript strand remain source-declared evidence where available,
not facts inferred from the displayed gene name or a TP53-specific rule.

Bind the selected gene separately to each available annotation source. Preserve
explicit cross-source mappings and unresolved relationships; equal symbols do
not prove Ensembl/RefSeq identity. Missing or ambiguous gene resolution makes
gene-relative interpretation unavailable with a typed reason, never a silent
fallback to another gene. Transcript/CDS selection, if requested later, must
also be explicit rather than an assumed canonical isoform.
Preserve a transcript's multiple `gene_ids`; membership in the selected gene
does not make that transcript or its reads exclusive to it. Phase 2
`expected_targets` may be cross-checked for disagreement but cannot supply a
missing selector or exclude an otherwise supported structural observation.

Keep one reference-gene context per request initially; this does not discard
off-target mappings or force ambiguous reads onto that gene. Retain competing
locus/gene evidence separately. Bind the selector to request/report identity,
not the physical splice-chain digest: switching the selected gene changes the
comparison context, not an unchanged observed chain's genomic identity.

All tables, figure headings and optional curated RNA aliases use the resolved
gene context. No `TP53`, `P53F`, p53 isoform label, fixed strand/chromosome or
gene-specific exon numbering is a core default. Such labels may occur only in
explicitly selected example data or an optional, source-bound alias profile.

### 1. Normalize Once, Preserve What Was Observed

Define a bounded engine-owned observation model consumed by one geometry and
comparison kernel. S1 input is the verified Phase 2 report, its hash-bound
`read_anchor_evidence.tsv`, source queries manifest and reference bindings.
Retain row/query/batch accounting and exclusions; a summary alone is insufficient.
An ingress-neutral internal observation is not a third annotation schema.

Reconstruct M/=/X-supported blocks and N-skips through one shared operation
walker factored from `paf::row`. Project `excluded_query_intervals` through those
blocks with `reference::project`. Do not infer genomic mask coverage from the
already unioned `aligned_query_intervals`. Empty operations or
`alignment_detail_verified=false` mean unassessed geometry, not an unspliced read.
Validate reconstructed query/reference spans and mask unions against their
stored counterparts; a matching file hash alone does not establish semantic
consistency. Contradictory geometry is an input error, not a reason to trust one
convenient serialized field and continue silently.
Do not use Phase 2 `gene_ids`, `transcript_ids`, `classification`,
`expected_dominant` or expected-target fractions as structural filters or model
matches: `reference::annotate` attaches genomic transcripts by whole-span/exon
overlap without a strand comparison, not by supported block compatibility.

`evidence::analyze` retains all deduplicated `segments` but does not serialize
`weight`, `prefix`, `suffix`, `best` or `optimal`. Factor its geometry-only DAG
calculation into a pure shared helper, used by both callers. Preserve its score,
tie handling, work accounting and old output bytes. The new result needs legal
successor/path relationships, not concatenation of the optimal-node union;
bound alternatives without enumerating an exponential number of paths. Leave
Phase 2 classification/annotation semantics unchanged. Geometry reconstruction
and this narrow factoring need byte-equivalence tests before a consumer ships.

Internally use explicit zero-based half-open genomic intervals and query
intervals. Store skipped intervals in genomic ascending order plus supported
flanking blocks; derive one-based exon labels only at presentation. SAM/BAM ingestion
must distinguish query-consuming operations, deletions and reference skips,
handle clipping and long-CIGAR metadata correctly, and never interpret a parser
placeholder as a biological intron. Prefer a maintained parser over a new BAM
decoder. See the [SAM/BAM specification](https://samtools.github.io/hts-specs/SAMv1.pdf).
Normalize inclusive GTF coordinates explicitly at ingress; see
[Ensembl's format definition](https://jun2026.archive.ensembl.org/info/website/upload/gff.html).

Keep read orientation, primer/payload orientation, genomic alignment orientation
and transcript strand as separate facts. Phase 2 query bases are normalized to
primer direction, its PAF `strand` is query-to-reference orientation, and its
annotated splice strand is supplied by the model. The current parser discards
`ts` and splice dinucleotide text; neither can be recovered from saved operations.
Strand inferences need an evidence source and alternatives. Annotation context,
declared primer design or later verified sequence context are not interchangeable
proofs. Apply transcript-direction ordering only to a supported display/relation;
keep chain identity strand-neutral. Unknown/conflicting strand or overlapping/
chimeric segments remain unresolved alternatives, not a forced linear transcript.
An alignment to a transcriptome model alone does not independently observe a
genomic splice junction. Supplementary/cross-locus evidence stays separate.

Require unmasked, query-supported flanks for qualified observed junctions;
report required flank lengths and quality filters as hashed request policy.
Never bridge primer, handle, masked or unaligned sequence into a biological
junction. Read-name membership in a primer list is not proof of primer removal.
For legacy SAM/BAM input, missing primer coordinates stay `unassessed`. Exact
internal skips, their chain and model matches may be inspected, with primer-,
flank- and end-qualification unassessed; never promote them to primer-free proof.

### 2. Reuse Annotation Geometry, Keep Comparison Sources Separate

Use `gentle.primer_payload_annotation.v1` / `PrimerPayloadTranscript` for
normalized geometry, not a parallel transcript model. A request-level
`annotation_sources` list binds each file hash, provider/resource/release and
source ID to the genomic reference FASTA/assembly. Provider metadata belongs in
that binding; add a versioned annotation extension only if a demonstrated field
cannot be represented without changing v1 meaning. Missing source metadata
stays unassessed, never inferred from an accession prefix.

Compare multiple provider annotations to the same genomic alignments; do not
require a second Phase 2 alignment batch merely to obtain RefSeq relationships.
Each source must validate its reference/assembly/contig bindings. Source gene
IDs are independently selected/mapped, including multiple gene memberships.
GTF and UCSC genePred conversion is a later explicit preparation step targeting
this same shape. Existing locus-presentation formats are not already converters.

### 3. Keep Interpretation Axes Independent

- **Topology:** spliced, unspliced, segmented/nonlinear or unresolved.
- **Reference relationship, per source:** exact full splice chain, exact partial
  consecutive chain, near-compatible with offsets,
  `all_junctions_annotated_no_single_model_contains_run`,
  `contains_unannotated_junction`, incompatible supported coverage, or unavailable.
  Exact full-chain agreement is not full transcript/end agreement or exclusivity.
- **Coverage/end status:** observed genomic span plus endpoint constraints such
  as `primer_bounded`, `masked`, `read_end_unverified`, `clipped`, `unknown`.
  Preserve contributing constraints rather than overwriting one with another.
  Payload endpoints are not natural RNA ends; 5'/3' labels require strand evidence.
- **Retention evidence:** query-supported intronic coverage, reference-relative
  candidate interval and exclusions, separately from junction compatibility.
  No `N`, a deletion span or absence of a compatible model is not retention proof.
- **Evidence adequacy:** reference/model scope, search/accounting completeness,
  mapping alternatives and primer-exclusion status; missing is not negative.

Exact model compatibility also checks supported blocks, including terminal
blocks in the corresponding first/last exons. A correct junction chain with
terminal intronic extension retains exact junction coordinates but gets a
terminal-coverage-divergent relation, not exact full/partial model compatibility.
Missing masks leave flank qualification unassessed; tolerance never excuses
unsupported block geometry. Unknown strand may support a strand-neutral match,
not identification of which sense/antisense RNA produced it.

Default compatibility tolerance is **0 bp**. Reproduce the audit's +/-5 bp only
as an explicit secondary policy, retaining exact observed coordinates, offsets
and every compatible model. Do not snap coordinates, greedily choose one match
or merge nearby chains transitively. The eleven audit chains remain ten exact
plus one near-compatible, not eleven identical observations.

Ensembl and RefSeq stay source-separated, with versioned transcript accessions,
assembly/contig authority, annotation hash, resource identifier and coverage.
An archived UCSC `ncbiRefSeq` track snapshot is not a verified direct NCBI
annotation-release/GCF binding. Do not guess aliases such as `17`/`chr17`, strip
accession versions or manufacture a consensus. With incomplete model/search
coverage, retain positive matches but withhold unsupported absence/novelty claims.
Always say **new relative to these bound references**, not newly existing biology.

### 4. Use An Identity Hierarchy, Not An Isoform Shortcut

Propose a generic `gentle.observed_transcript_structures.v1` report with:

- An **observation ID** bound to source artifact/record, alignment hypothesis,
  full block geometry, endpoints and primer/exclusion evidence.
- An **`observed_junction_chain` ID**: full SHA-256 of a versioned canonical
  encoding of the bound genomic reference identity, assembly/contig and ascending
  exact skipped intervals. The request's reference label is not coordinate
  authority by itself; retain reference FASTA/contig binding. Transcript strand,
  reference models, gene choice/symbols, read ends, counts and tolerance do not
  define this strand-neutral physical chain. Do not call it an RNA isoform.
- Each chain's required companions are observed genomic span, endpoint constraints
  and sourced strand alternatives. Partial chains are distinct observed families;
  containment/compatibility links never merge them. Common physical skips on
  opposite strands are not proof of one RNA; retain strand subgroups/ambiguity.
- Separate coverage/retention subgroups within a chain family; preserve their
  per-read evidence rather than treating them as identical complete RNAs.
- For verified zero-junction reads, retain observation IDs but no chain-family
  or fragment-group ID in v1. Report model-relative supported coverage such as
  within-exon, exon/intron-boundary-spanning or intronic, plus span/end status.
  Coarse/missing operation detail cannot establish zero junctions. Retention-like
  coverage is a source-relative subgroup, not a separate RNA family or mechanism.
- Gene assignments and literature-backed labels as versioned relations, retaining
  ambiguity. Short hash labels are display aliases with collision handling, not
  sole identifiers. Archive `P53F-RNA-v1` keys remain opaque legacy aliases.

This improves Glen's naming idea without changing the archived nomenclature.
`p53alpha-compatible RNA` may be a later curated alias, never proof of p53alpha
protein. p53 protein forms can also depend on translation initiation, not only
RNA splicing; see the [NCI TP53 isoform reference](https://tp53.cancer.gov/p53isoforms).

### 5. Reconcile Records, Then Count Support

Count unique bound original read records, not anchors, PAF/BAM rows or repeated
read names. Preserve all compatible structures for ambiguous records, with
definite/possible support and no automatic fractional or exclusive assignment.
Multiple anchors and source-reference comparisons must not multiply support.

Reconcile hypotheses by bound `(source file identity, record_index)`, validating
the recorded read ID. Exact chain agreement from several anchors contributes
once; retain their separate endpoint/mask evidence. Shorter chains caused by
different payload bounds remain containment/partial relations, not contradictory
full isoforms or a unioned inferred full chain. Genuinely conflicting chains
remain possible alternatives with shared record support. Definite junction
support requires agreement among admitted complete alternatives; incomplete
accounting cannot establish a negative. Never add counts of nonexclusive rows
to obtain a library total.

Group additionally by declared sample, barcode, library preparation, priming
condition and primer identity; missing metadata remains unknown. Two reads in
one barcode are recurrence, not two independent libraries or original molecules.
UMI/PCR deduplication, expression estimation and enrichment inference are out
of scope. Denominators describe the selected input cohort, with exclusions and
unassessed records counted, not an unprovided whole experiment.

### 6. Publish Portable Evidence, Not Hidden Session Interpretation

Prepare `report.json`, `read_structures.tsv`, `structure_catalog.tsv`,
`reference_compatibility.tsv` and `junction_catalog.tsv` with source/request/
policy/version hashes, exact coordinate labels, ambiguity, limitations and
nonclaims. Tables derive from the same typed result; read sequences are omitted.
Bind inspection/rendering to verified stored results, not a fresh alignment,
download, silent annotation substitution or execution of the supplied scripts.

Specify hard bounds including `max_input_bytes`, `max_output_bytes`,
`max_line_bytes`, `max_records`, `max_operations`, `max_skips_per_chain`,
`max_chain_families`, `max_alternatives_per_record`, `max_annotation_rows`,
`max_model_comparisons_per_family` and `max_total_comparison_work`.
Account for both per-family model work and per-record coverage/mask work; repeated
reads must not bypass limits or force reads x models x skips materialization.
Count completed/failed/unassessed units, poll request-bound cancellation, and
distinguish incomplete work from biological absence. Use fresh-directory staged
publication: no final report after a data error, and explicit incomplete status
for verified completed units on authorized cancellation. Do not silently truncate.
Reuse Phase 2 file-only publication and approval patterns, not a broad assertion
that a read-only operation cannot repair state. Filesystem writes need explicit
authorization even when project state is unchanged. Test constructor/workflow/
autosave behavior, including MCP `confirm:true` and named no-save handling:
assert project bytes and journal unchanged for both CLI/MCP and preserve the
documented constructor/state-out scope. [Lifecycle Slice 3](operation_lifecycle_slice3_plan.md)
landed independently on main at `f0b8371f`: its covered operations are
`PlanDnaSequenceDesign`, `EvaluateGenomeEditIntent` and qualifying workflows,
not a generic read-only guarantee. Future observed-transcript S1 operations
are not automatically covered and still need their own lifecycle/publication
decision and state/journal tests. S0's pure kernel adds no such operation.

## Slice Boundaries And Remaining Work

1. **S0 -- implemented pure kernel; retain its regression guards.** Synthetic observations, canonical geometry,
   stable-gene-ID resolution, strand-neutral chain identities, independent axes,
   record reconciliation and exact source comparison. Factor the operation walk
   and geometry-only DAG into shared helpers with no Phase 2 schema/byte/behavior
   change. Record an owner-reviewed decision before freezing the new schema.
   No I/O behavior change, BAM, project migration, new dependency or GUI.
2. **S1 -- consume bound Phase 2 evidence.** Verify the report, evidence TSV and
   query manifest, reconstruct geometry/hypotheses through shared helpers, and
   compare independently bound annotations. Ignore old classification/ID fields
   for structural filtering. Add shared engine/parser, bounded inspect/cancel/
   export, CLI/MCP/agent discovery and generic JS/Lua/Python reachability. Test
   file-write confirmation, named no-save behavior and unchanged project/journal.
   Keep old schemas/results byte-compatible; source-bound execution is required.
3. **S2 -- explicit SAM-text ingress and annotation preparation.** Prefer a
   bounded, hash-bound SAM input initially, with CIGAR geometry through the same
   kernel. Handle `S/H/P`, reverse flags/query coordinates, header lengths,
   unavailable MAPQ and malformed/unsupported records explicitly; do not broaden
   the existing PAF grammar. No implicit conversion, subprocess or aligner.
   Gene-scoped GTF/UCSC genePred normalization is a separate explicit prepare
   step targeting the reused annotation shape, not a third transcript schema.
   Bind owner-produced `samtools view` output to its input BAM hash, command,
   version and conversion/selection scope; conversion is not an aligner replay.
   Preserve source-record versus producer-declared QNAME grouping and duplicate
   names; never invent original FASTQ/Phase 1 identity. Missing masks stay
   unassessed. The 21-record archive replay becomes possible here with these
   additional SAM bytes; it is not currently replayable by S1. Binary BAM remains
   a later parser-dependency decision. Aggregates are cross-checks, not targets.
4. **S3 -- present saved evidence.** Distinct reference and observed-structure
   lanes in/alongside Splicing Expert, genomic/transcript orientation, clickable
   support/provenance, explicit gaps/ambiguity and a shared SVG/PDF view. Reuse
   render/lane primitives, not the project-bound `build_splicing_expert_view` or
   symbol-keyed `splicing_group_label`. Build lanes from the verified saved report
   with an explicit coordinate transform. First add a synthetic tutorial; a real
   TP53 lesson requires publication permission and source-bound native/scientific
   review. No re-analysis or project dependency introduced during painting.
5. **Later, separately reviewed:** curated RNA aliases and conditional coding
   consequences using exact CDS/phase/genetic-code/reference bindings and shared
   translation logic. The audit's assumed upstream canonical sequence and 343-aa
   prediction remain conditional, not part of initial structural inference.

The first useful delivery is S0/S1, not a TP53-specific port. S2 consumes an
explicitly prepared SAM companion, not the ZIP/BAM automatically. The missing
conversion artifact is a bounded follow-up, not an excuse to infer its contents.

## Required Deterministic Regressions

- Run the same kernel for at least two independently selected synthetic genes
  at distinct loci and on opposite strands, without any TP53 data or alias
  profile. Assert selector propagation across adapters, source-specific ID
  resolution, typed missing/ambiguous-gene results and correct figure/table labels.
  A context change alters request/report identity, not an unchanged chain ID;
  different physical loci must not collapse into one chain.
- Include overlapping sense/antisense genes, different source IDs sharing a
  label, a transcript with two `gene_ids`, a selected gene absent from Phase 2
  `expected_targets`, and supported read blocks nested in the selected gene's
  intron. No default, whole-span label or Phase 2 classification may decide the
  selected gene/strand/compatibility. Unknown sense remains ambiguous.
- Prove chain identity survives opposite selected annotation strands at the same
  physical skips; preserve strand alternatives instead of declaring one RNA.
- Byte-identical old Phase 2 results after shared helper extraction, including
  reverse geometry, masks, ties, coarse rows and budget boundaries. Multiple
  optimal paths must not create a synthetic concatenation of incompatible nodes.
- Multiple anchors from one source record: identical, contained-partial and
  conflicting chains, including reverse-normalized overlapping payloads; one
  record remains one contribution with endpoint/alternative evidence intact.
- Exact canonical chain versus a partial consecutive chain: both may match,
  neither establishes protein or exclusive isoform; full chain can lack RNA ends.
- Two reads sharing a new junction: one exact chain, two record supports; repeat
  anchors/mappings/reused names do not become independent support.
- Minus-strand and reverse-query cases, asymmetric exon lengths and explicit
  zero-/one-based boundary examples, including the archive's left-field mismatch.
- Primer/handle/internal masks, missing mask bindings and insufficient flanks;
  no artificial primer junction, no primer-free claim from a read-ID list.
- `D` versus `N`, split reference spans, insertion/clipping, retention-like
  aligned coverage versus truly unspliced fragments; no universal empty-chain ID.
- Shifted boundaries remain distinct IDs; exact/near relations are separate;
  multiple tolerance matches are retained and nontransitive cases do not merge.
- All individually known junctions in a new combination versus an unannotated
  junction; exact skips with terminal intronic extension; unsupported intervening
  blocks and ambiguous alignments; no forced first transcript or exact model.
- Missing provider, wrong assembly/contig/version/hash, coarse alignment,
  incomplete search or curated cohort: typed unassessed/limited result, no absence.
- Transcriptome-only and supplementary/cross-locus records cannot invent observed
  linear genomic transitions. Include non-TP53 and plus-strand synthetic cases.
- Stable full hashes, short-label collision, reordered input determinism, stored
  report/table verification, resource/cancel failures, no-overwrite publication,
  LF/CRLF/path portability and CLI/MCP approval/state/journal-preservation parity.
- SAM-specific tests in S2: clipping, padding, reverse orientation, unsupported
  long-CIGAR placeholders, duplicated QNAME versus distinct source records,
  conversion bindings and selected-cohort denominator; unchanged PAF admission.

Synthetic fixtures must have nearby origin/recreation/usage documentation.
Actual binary and platform executions remain future evidence, not this plan's
acceptance. Do not replace failed invariants by matching the audit's counts.

## Inputs, Decisions And Consolidated Verification

Before S2, request the explicit hash-bound SAM companion and conversion record,
then only missing provenance: original selection policy/denominator,
sample/library identities, primer-coordinate or trimmed-query bindings, and the
candidate re-alignment command/inputs if a primer-free claim is required. Genome
and full annotation hashes are present in `input.sha256`; obtain their exact
bytes only for checks that actually need them, never fetch a newer substitute.
Unrecoverable provenance remains explicitly unassessed, not a permanent demand
to reconstruct a vanished experiment. Independently verify assembly/contig and
annotation-snapshot authority before claiming cross-source novelty.

The S0 authorization covers strand-neutral physical-chain identity with sourced
strand alternatives, retained observation IDs but no unspliced family IDs,
stable-gene-ID selection, and narrow shared DAG/walker factoring under
byte-equivalence tests. S1 operation/ingress/publication remains separately
unapproved. SAM-first is the preferred S2 route,
not a requirement blocking pure S0 work; binary BAM needs a separate decision.
Nonzero tolerance, flank thresholds, real-fixture publication and conditional
CDS interpretation can be decided when their slices are reached. No universal
MAPQ or retention cutoff is adopted here.

After a committed implementation candidate, the GitHub-pushing thread should
send **one consolidated request to Glen**, not one per branch: exact source,
binary/toolchain/profile and input hashes; synthetic engine/adapter runs;
independent archive replay separating ten exact from one near canonical chain;
two-read candidate and primer-removal evidence scopes; native screenshots of
reference/observed lanes when implemented; and platform/scientific checks as
separate typed verdicts. Glen returns a manifest-bound report through Steffen.
No request is sent, no candidate freeze/dispatch is created, and no native,
release, biological or performance acceptance is claimed by this planning task.

Next step after S0 verification: obtain owner approval for S1. The updated
[Claude follow-up prompt](observed_transcript_structures_claude_review_prompt.md)
is available if a second read-only pass on the remaining qualifications is
wanted; it has not been sent or executed. S1, SAM/BAM, GUI and later consequences
are not authorized by the S0 request.
