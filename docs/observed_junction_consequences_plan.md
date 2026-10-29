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
| `near_window_bp` | Required 1..50 engineering admission window; no default. Acceptance's proposed value 5 is not selected here. |
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
LCG `state = state*1664525 + 1013904223 (mod 2^32)`, seed 274, bases selected
from the upper two bits; planted motifs and declared geometry are fixture-only.
Literal CD274 coordinate/count/recurrence acceptance stays gated below.

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
| Meaning of 3/1/3/0 | Confirm records versus reads, certainty and flank qualification; provide library/prep/sample/cell-line mapping. Do not infer independence. |
| Acceptance recurrence | Recommend explicit library/min_units=2/definite only if independent; shared preps require prep/sample units. Await owner. |
| Acceptance near window | Recommend explicit 5; no default is introduced. |
| Real CD274 ingress | Phase 2 -> S1 is available. BAM needs separately authorized verified SAM ingress plus external conversion; inherited private decoder is not a service. |
| Real reference authority | Full assembly needs explicit large byte/work budgets; excerpt has declared authority only. |
| J2 CDS format | Decide separate binding versus new S2.1 output; existing receipts stay unchanged. Decide whether to expose stop-to-last-junction distance as geometry only. |
| J3 | Decide minimum extracellular residues, SIGNAL requirement, single-pass-only v1 and portable UniProt format. |

Once confirmed, a generic synthetic fixture may use two independently bound
CD274-like sources that lack the chosen junction, including shared-last-donor
geometry, four declared libraries with gene-coherent denominators above support,
verified gt_ag, explicit recurrence policy and incomplete-search 0 definite / 3,1,3,0
possible support. No human/GRCh38/CD274/plus/soluble default enters the service.
Without confirmation this acceptance case is pending, not silently passed.

## J1 Execution Ledger

J0 is committed separately as `f0c8bdb7`. The owner then clarified the legacy
coordinate convention above; support meaning/hierarchy and synthetic recurrence
settings remain unanswered at this point. Code contains no CD274 default.

At `f0c8bdb7` plus the initial DTO/test skeleton, the fail-first command
`env CARGO_INCREMENTAL=0 CARGO_BUILD_JOBS=1 cargo test -q --locked --offline --no-default-features --lib junction_assessment -- --test-threads=1`
exits 101 with four E0599 errors for absent assess/binding/inspection methods.
This is an actual compile-red result, not executed tests or a zero-match pass.
The expanded tests additionally specify near/event/motif helpers and file-bound
reference/manifest/cancel/error behavior before their implementation.
The unchanged Cargo.lock SHA-256 is
`a393a1f3803b374bbaf5ef4f3e23813e974b1676f61a669dacb24063c059fb6a`.
No successful J1 runtime verification is claimed yet.
