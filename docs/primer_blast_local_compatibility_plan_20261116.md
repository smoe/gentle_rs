# Local Primer-BLAST-Compatible Checks Of Supplied Primer Pairs

Status: Slice 0 implemented and locally verified; later slices remain proposed. No
operational compatibility or scientific acceptance is claimed.
Prepared 2026-11-16 against `b62d21c9` (`main`, `origin/main` and this branch).
Revised 2026-11-17 after the owner supplied Claude's read-only review against
that baseline. Review disposition and the original plan's
ordering are retained below. The owner subsequently authorized only Slice 0 and
its commit; the source/test files and execution receipt are recorded below.
This is not an additional `.13` release gate.

## Goal And Non-Goals

Check supplied human primer pairs locally against a pinned human RefSeq-mRNA
snapshot, including predicted `XM_` records. Reuse GENtle's specificity engine,
reference admission, external execution receipts and inspection routes. An agent
schedules the declared work; GENtle owns product interpretation and verdicts.

This is documented compatibility with a dated public Primer-BLAST policy, not
an exact replica of NCBI's server, database or alignment implementation. It does
not add primer redesign, arbitrary-organism discovery, Entrez query support,
SNP masking, a web scraper, automatic reference updates or a new report builder.
Transcriptome specificity remains separate from genomic carryover and from
experimental PCR validation. Existing Tm and oligo-QC calculations stay unchanged.

## What The Current Source Actually Provides

| Surface | Evidence at the inspected source | Consequence |
| --- | --- | --- |
| Aligned BLAST fields | `src/genomes.rs` constants and parser; `src/engine/ops/operation_handlers.rs::primer_specificity_full_alignment_from_blast_strings` | `qlen/qseq/sseq` and complete-query HSP interpretation already exist on main. Do not assume an old branch needs cherry-picking. `4f741797` was not resolvable locally and was not reviewed. |
| Full-primer fallback | `primer_specificity_full_alignment` in the same handler file | An expanded, strand-aware subject-window semiglobal aligner exists. Extend/test it; do not create a competing aligner in an adapter. Complete HSPs currently have a separate aligned-string path, not a fresh global optimization. |
| Reference reader gap | `src/genomes/prepared_sequence_reader.rs`, `GenomeCatalog::read_prepared_sequence_region` and `list_chromosome_lengths` | The shared reader/cache/alias path assumes a genome manifest/index and compatible-assembly resolution. Sequence-only preparation writes `blast-resource-manifest.json` without a FASTA index. Add resource index preparation, resource-aware cache validation, exact resource resolution and accession lookup together; do not fabricate genome manifests. |
| Sequence-resource preparation | `prepare_blast_sequence_resource_with_options`, `adopt_blast_sequence_resource`, `normalize_genbank_blast_resource_if_needed` in `src/genomes.rs` | Reuse BLAST inspection, fingerprints and source protection. Current GenBank conversion loads all records, replaces the materialized input with FASTA and hashes the converted output; it is not the proposed provenance-preserving streaming RefSeq intake. FASTA annotation parsing specializes in Ensembl headers. |
| Search parameters and import | Direct search in `src/genomes.rs`; handoff builder/import in `operation_handlers.rs`; child validator in `primer_specificity_multi.rs` | Direct/handoff paths hard-code E-value1000; multi-import independently requires those arguments and N+1. Single-handoff identity omits command arguments. Profile-specific argument validation must cover every path, not merely the builders. |
| Project-level hit filtering | `resolve_blast_options_with_layers`, `apply_blast_thresholds_to_report` in `src/engine.rs`; handoff output parsing in `operation_handlers.rs` | Project/default thresholds can remove HSPs before interpretation, including by E-value, identity, qcovs, alignment length or bit score. The profile must explicitly resolve these, not inherit a hidden candidate veto. |
| Present mismatch policy | `PrimerSpecificityPolicy` in `crates/gentle-protocol/src/lib.rs`; hit interpretation and product finalization | Default terminal mismatch maximum is zero; product interpretation uses combined mismatches. Hits surviving the BLAST wrapper still enter pairing; terminal flags affect finalization and later compaction. The earlier claim that the terminal veto removes these hits before pairing was incorrect. |
| Product geometry and identity | `primer_specificity_collect_amplicons_for_hits_with_stats`, `_finalize_amplicons` | F-R, F-F and R-R and both spatial orientations exist. The inspected path sorts but does not deduplicate HSP-derived sites/products. Intended recognition requires exactly one matching product. New-profile site consolidation is required, not existing behavior to reuse. |
| Completeness and import | `primer_specificity_search_completeness_for_commands`; `src/engine/analysis/primer_specificity_multi.rs`; `primer_specificity_references.rs` | Existing validated subject-count, execution-receipt and reference-binding gates must remain. Multi-reference import independently validates child commands/completeness. |

Read the existing [comparison](primer_blast_comparison.md),
[CLI route guide](cli/primer-specificity.md), DEC-029 and DEC-049 before coding.
This audit was source-only; it does not prove those paths pass at this SHA.

## External Behavior Verified On 2026-11-16

The [NCBI form](https://www.ncbi.nlm.nih.gov/tools/primer-blast/) was inspected,
including raw HTML selected options/default attributes, not only rendered help.

| Setting | Observed public default |
| --- | --- |
| Database | RefSeq mRNA |
| BLAST E-value | 30000 |
| Word size | 7 |
| Ungapped BLAST | Hidden `UNGAPPED_BLAST=on` |
| Returned subjects | 50000 |
| Target amplicon ceiling | 4000 bp; distinct from design product-size controls |
| Mismatch stringency | At least one primer: total >= 2 AND last-five-base mismatches >= 2 |
| Ignore threshold | At least one primer has total mismatches >= 6 |

These are public interface settings, not proof of identical backend behavior.
The [Primer-BLAST paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC3412702/)
describes BLAST-seeded full-primer alignment. It does not establish equivalence
between GENtle's scoring/tie handling and the current server.

The [RefSeq FAQ](https://ncbi.nlm.nih.gov/books/NBK50679/) distinguishes modeled
mRNA (`XM_`) from modeled noncoding RNA (`XR_`). The
[NLM download guidance](https://support.nlm.nih.gov/kbArticle/?pn=KA-03450)
distinguishes assembly-annotated collections from cumulative species collections,
which can contain records not annotated on the reference assembly. A GRCh38 RNA
FASTA is therefore not, by itself, proof of the required sequence universe.

## Biological Contract

1. Resolve a separate named compatibility policy. Omitted policy keeps existing
   behavior. Do not merely loosen `max_3prime_mismatches`.
2. Preserve admitted BLAST candidate evidence through complete-primer alignment.
   Do not inherit generic project hit thresholds or let legacy terminal/coverage
   flags veto the new policy. The existing terminal flags are not themselves a
   pre-pairing discard; wrapper thresholds and later compaction are separate.
3. Construct inward-facing F-R, F-F and R-R candidate products through 4000 bp,
   using the existing geometry after new-profile site consolidation. Preserve
   raw HSP identities as evidence when several HSPs resolve to one site.
4. For an unintended product let `m1,m2` be per-site total mismatches and
   `t1,t2` the differences in each oligo's own 3'-terminal five bases. Apply:

   ```text
   ignored_by_total_threshold = (m1 >= 6) OR (m2 >= 6)
   excluded_by_stringency = ((m1 >= 2) AND (t1 >= 2))
                         OR ((m2 >= 2) AND (t2 >= 2))
   relevant_off_target = NOT ignored_by_total_threshold
                         AND NOT excluded_by_stringency
   ```

   The predicate belongs to the product policy after pairing, not a combined
   mismatch score. It is mathematically equivalent to the conjunction of two
   site-eligibility predicates; test that equivalence, but keep classification
   and explanation at product level as Glen requested. For the first classifier,
   use complete, unambiguous ungapped evidence and require `0 <= t <= m`.
   Read counts from the authoritative full alignment, not HSP mismatch fields
   which can remain different on the complete-HSP shortcut.
   Retain per-primer evidence and exclusion reasons for audit.
5. Bind intended transcripts explicitly by accession.version and reference;
   same-gene splice variants are not silently intended. An ignored candidate
   is not a demonstrated absence of binding or an experimental safety claim.
6. Missing/ambiguous alignment, reference drift, failed execution or unproven
   required search coverage prevents a pass. A detected off-target remains
   visible even when other evidence is incomplete. Preserve current report
   status vocabulary and separate execution failures from biological findings.

Indels, ambiguous bases, equal-scoring alignments and transcript-end windows need
explicit documented handling before claiming parity. Do not silently count gaps
as substitutions or choose the most favorable alignment to manufacture a pass.
If source evidence is insufficient, retain conservative incomplete/unresolved
evidence rather than inventing an NCBI rule. Align the declared binding segment;
preserve any explicitly declared non-annealing 5' tail and oligo provenance.
Counting gap columns in an effective difference total is an explicit possible
GENtle policy, not a verified NCBI rule. Unsupported gaps/ambiguity stay unresolved
until their treatment is specified. A clipped valid transcript window does not
justify inventing missing reference bases or calling an HSP-span substitute
complete. Intended-target recognition, compact output and summary/readiness must
consume the resolved policy rather than silently retaining the legacy veto.

Listing a candidate and classifying it as a relevant off-target are separate
outputs. Do not assume the supplied-primer web report hides every candidate
excluded by the design stringency rule. Acceptance must say which of those
behaviors is being compared with Glen's saved output.

## Completeness And The 50000-Subject Tension

Keep GENtle's validated subject-coverage requirement and existing
`validated_sequence_count + 1` limit for each primer. The original proposal's
50000 floor adds no coverage and is dropped. Record the website default and the
actual N+1 limit as a deliberate local difference. A returned-subject
cap, HSP/pairing budget, timeout or report display limit must never silently
truncate the evidence behind a pass. Existing historical evidence is not
retroactively reinterpreted.

If a caller explicitly restricts search to 50000 below the required local bound,
the result is incomplete, not exact-compatible and globally clean. No separate
web-cap emulation mode is needed in the first slice.

Subject coverage is not proof that heuristic BLAST finds every possible binding
site. Reports/documentation must distinguish corpus completeness, subject-cap
admission, candidate-search sensitivity and full-alignment availability. A small
brute-force fixture oracle can test sensitivity; it cannot prove genome-wide
absence of false negatives. Keep this limitation explicit without rewriting
unrelated historical completeness schemas.

For a contiguous exact-word search of size 7, a 20-base alignment with substitutions
at positions 7 and 14 has three exact runs of six: there is no seed on that
diagonal, although the mismatch rule may still regard its product as relevant.
Other loci/diagonals can incidentally nominate the window; that is not guaranteed.
An available exact seed is necessary for this seeding path, not a guarantee that
BLAST reports the site after extension, scoring, masking and thresholds.
The existing `cdna_assay_find_primer_hits` offers an exhaustive ungapped toy-
fixture comparator; it is not an independent global/gapped alignment oracle.
Record this sensitivity limit even when the existing subject-cap field says
`complete`. A requirement for exhaustive mismatch-site enumeration is not met
by this profile and must not acquire an unconditional specificity claim.

## Implementation Plan (Revised 2026-11-17)

### Slice 0: Legacy Identity Locks And A Pure Classifier

- First pin representative baseline stored policy, single/multi handoff, report,
  multi-summary and panel-selection JSON. Deserialize and recompute the actual
  stored identities, summary digests and selected-report hashes, not just an ID
  derived from a newly constructed default policy. Assert that selection still
  resolves, rather than silently returning `None`. Expected values come from
  the pre-change baseline, not regeneration after adding fields.
- Add an unwired pure classifier for complete ungapped evidence, its product
  conjunction, truth-table and decomposition-property tests. No request/schema
  field or new CLI parameter belongs in this first commit. No resource access,
  scientific default change or operational compatibility claim.
- New fields introduced later must be absent from legacy digest projections
  when not opted in. `Option` with `default` and `skip_serializing_if` is the
  preferred representation for new profile/evidence records, not the only valid
  Rust type. An omitted default-valued field can also work if the real legacy
  serialization/hash tests prove it. Do not change existing fields' serialization
  or use their always-serialized defaults as the pattern for new fields.

### Slice 1: Named Policy And End-To-End Interpretation

- Add an opt-in, versioned `primer_blast_compatible_v1` policy, with the requested
  `primer_blast_default` spelling only as a documented alias to that immutable
  definition. Resolve it in shared code and persist the effective parameters,
  policy version and defaults-observed date. Conflicting overrides must be
  rejected or visibly recorded as a customized profile, never silently hidden.
- Wire the tested predicate through intended recognition, off-target verdicts,
  summary counts and compact/full retention. Legacy terminal flags and combined
  mismatch thresholds must not continue to determine the new profile's verdict.
  Preserve existing omitted-policy behavior. Unknown sites cannot become negative
  evidence, and a relevant product must not disappear solely because a legacy
  `accepted_by_policy` flag is false.
- Pin both exploratory and ordinary readiness ceilings to 4000 bp in this profile;
  inherited shorter assay ceilings must not turn the proposed 2751-bp GFI1 product
  into a warning-only result. Explicit customization changes the resolved profile
  identity and its compatibility claim.
- Consolidate equivalent aligned sites deterministically before pairing, using
  exact reference/subject version, oligo role, strand, aligned coordinates and
  alignment evidence. Keep supporting HSP references. Do not collapse distinct
  sites, accessions or competing alignments by choosing fewer mismatches. Then
  ensure equivalent site pairs cannot duplicate a product or defeat unique
  intended recognition. This is opt-in work, not a legacy deduplication repair.
- Start with deterministic constructed alignment/product fixtures: one terminal
  mismatch, two terminal mismatches, 2+0 and 1+1 distributions, five/six total
  mismatches on either site, terminal-window edges, both strands and F-F/R-R.
  Add a case retained by the new policy but rejected by the legacy terminal veto.
- Gate this slice on Slice 0's stored-record identities plus unchanged legacy
  verdict fixtures. Do not expose an operational compatibility claim before
  Slices 2-3 close the search and resource boundaries.

### Slice 2: Prepared-cDNA Reader And Resolved Search Profile

- Extend the existing prepared-region reader to admit sequence-only resource
  manifests and a prepared FASTA index, with exact accession.version lookup.
  Resource preparation currently builds no such index; add/bind it without
  manufacturing a genome manifest. Extend the cache and final validation to the
  resource manifest, with exact resource-ID resolution rather than compatible
  assembly substitution. Resolve lengths/aliases from that resource's own index;
  pin actual RefSeq `sseqid` formatting with an optional real BLAST+ fixture.
  Reuse per-interpretation index caching and before/after validation. No whole-FASTA
  scan per hit, no process-global cache, invented genome artifacts or aliasing
  away accession versions. Record retained index/source identities.
- Gate new sequence-resource reads on the opt-in profile. Otherwise today's
  fallback-based legacy results can change under the same content-derived ID.
  Keep old records inspectable and old requests on their original path; any
  separate correction to legacy interpretation requires explicit approval and
  its own versioned identity/migration decision.
- Test partial HSP recovery on both strands and near transcript ends. A valid
  clipped window is different from an unavailable reference. Clamp against the
  admitted subject's indexed length, retain requested/actual bounds and reject
  invalid or insufficient evidence. The current fallback can mark HSP-span-only
  alignment `complete`; do not rely on `mode=required` to catch it. Keep this
  profile's changed interpretation separate from the legacy path.
- Reuse the aligner. For this profile, run full-query optimization over the
  admitted window; do not assume complete BLAST strings are equivalent merely
  because every query base appears. Pin scoring, tie rules, ambiguity and gap
  handling, and compare small examples with an independent oracle. Any departure
  from NCBI's documented behavior is named, not claimed away.
- Evaluate, rather than assume, Claude's alternative of extending each ungapped
  HSP diagonal across the whole primer and marking better gapped alternatives
  unresolved. A diagonal defines a placement for that HSP, not a unique site in
  the transcript. Equal-score placements and shifts still need evidence and
  deterministic treatment. Decide the alignment rule before wiring verdicts;
  neither this option nor rust-bio's traceback proves NCBI equivalence.
- One shared resolver must produce direct and handoff parameters: E-value 30000,
  word size 7, ungapped candidate search, explicit masking/scoring/effective tool
  version, 4000-bp product ceiling, required alignment and the safe subject limit.
  Retain current sensitivity limitations; do not equate `-ungapped` candidate
  search with prohibiting all gaps during later realignment.
- Include the multi-import validator's exact expected arguments in that shared
  resolver. Add profile-specific single-handoff validation for both exact roles,
  queries, required options and output bindings; its existing handoff ID does
  not cover commands. Tests must reject a compat label with E-value1000, altered
  seeding/masking/caps, a missing/duplicate command or mismatched query evidence.
- Pin generic BLAST post-filters explicitly to non-discarding settings for this
  profile, or reject conflicting settings before execution. Both direct complete
  search and imported output currently apply layered thresholds before pairing.
  Test a restrictive project/default configuration: it must not remove candidate
  HSPs silently, including via qcovs or `unique_best_hit`. Record the resolved
  settings so hidden project defaults cannot change scientific identity.
- Update import validation, not just command generation. Modified profile,
  commands, raw HSPs, references or policy must invalidate pending evidence.
  Review existing policy/report/handoff schemas before choosing additive fields
  versus a version bump; never silently change old digest semantics. Missing
  profile fields continue to deserialize with legacy meaning.

### Slice 3: Versioned Human RefSeq-mRNA Preparation

- First implement/test local intake of a tiny source with matching sequence and
  annotation records, reusing `genomes prepare-blast-resource` and inspection.
  Extend only the missing RefSeq annotation/filtering/manifest pieces. The current
  in-place, all-record GenBank converter does not meet this intake contract.
  Stream source records into distinct staged derived outputs; retain and hash
  original compressed/downloaded bytes plus the extracted/normalized products.
  Publish only validated outputs and clean up only files owned by this call.
- Define the corpus as human taxon 9606 mRNA, including `NM_` and `XM_` accession
  versions; exclude noncoding `NR_/XR_`, proteins and other taxa with recorded
  reasons/counts. Verify molecule/taxon metadata; prefixes and gene symbols alone
  are not sufficient admission evidence. Do not deduplicate distinct accessions
  merely because their sequences match.
- Prefer one documented species-wide NCBI snapshot source with matching GenBank
  annotation, rather than assuming an assembly RNA FASTA is comprehensive.
  Record exact member files, snapshot/release or capture date, source URLs,
  retained bytes/hashes, extraction rules, accession manifest and sequence counts.
  A capture date is not an invented annotation release or assembly assignment.
- Preserve accession.version, record description, NCBI GeneID, symbol/aliases,
  modeled/curated status and annotation provenance. Unknown annotations remain
  unknown, not an Ensembl-style parse guess or proof of intended identity.
- Record a taxon/molecule/accession-prefix census from admitted source records;
  catalog text alone does not establish RefSeq origin. Use typed NCBI GeneID
  linkage rather than requiring consumers to guess from mixed db_xrefs. The
  existing `genbank_record_chromosome` prefers VERSION, then ACCESSION/name;
  the new intake must require accession.version instead of using those fallbacks.
- Build/adopt the existing local BLAST resource and subject annotation sidecar;
  validate index subject count and fingerprints. Prepare new immutable resource
  IDs for updates, show added/changed/removed accessions, and retain old snapshots
  for old reports. Updating availability must not activate new panel acceptance.
- Preflight destination permissions/free space, source member sizes where known,
  projected expanded/index/staging space and required BLAST tools before asking
  to download. Storage Box is an optional configured location, not a hard-coded
  path. Report unknown estimates; get owner approval before substantial network
  transfer or replacing anything. No automatic or scheduled updater is required.

### Slice 4: Machine-Interface Integration And Acceptance

- Extend existing specificity request/handoff/import and inspection routes;
  CLI and MCP/inner agent invoke the same engine policy. No new redesign command
  or GUI workbench is needed. Preserve aggregate panel acceptance guards and
  separate genome/transcriptome verdicts; a standalone result does not silently
  become accepted panel evidence.
- Document the profile, all differences from hosted Primer-BLAST, RefSeq corpus
  contents and the exact local invocation after its syntax is implemented.
- Add CLI parser/executor, MCP pass-through and serde/digest round trips. Prove
  direct execution and imported execution yield the same scientific verdicts;
  timings/paths are not biological identity. Unannotated new capabilities must
  still fail the global introspection/parity guards.
- Finish the source-backed fail-first cases below, then optionally run a
  separately authorized prepared full-human snapshot check. Tiny CI fixtures
  cannot establish completeness of that real snapshot or live-server parity.

## Acceptance Fixtures And Required Handoff From Glen

No matching GFI1/TBCE evidence was found in the searched tracked docs, fixtures,
tests or analysis modules. Gene names plus product lengths do not determine the
primer sequences, accession versions, strand or mismatch positions.
Different database providers alone do not prove GFI1/TBCE sequence absence from
Ensembl. Compare the actual record/sequence inventories before attributing a
miss to corpus coverage rather than candidate search or product policy.

Before claiming the following real-data acceptance, obtain both exact 5'-to-3'
primer pairs with tail/binding interpretation, the PATZ1/GFI1/TBCE accession.version
records and source bytes, and Glen's saved Primer-BLAST result/settings/date.
Retain the observed database identity where available. Declare missing server
corpus/version information; do not manufacture a matching release.
Ask whether each expected product must be listed, classified as a relevant
off-target, or both, and identify the intended PATZ1 accession versions explicitly.
The existence of `allow_same_gene_splice_variants` in a request is not evidence
that an interpretation path implements it; verify actual target recognition.

| Supplied pair | Glen's expected products, pending source-backed confirmation |
| --- | --- |
| Original PATZ1 pair | PATZ1 72 bp and GFI1 2751 bp |
| Replacement PATZ1 pair | PATZ1 126 bp and TBCE 167 bp |

Use full transcript records (or explicitly delimited test windows with coordinate
offsets and no whole-transcript claim), exact mismatch positions in oligo 5'-to-3'
coordinates, product boundaries and orientations. Preserve multiple isoform
accessions. Document redistribution rights and deterministic retrieval in a
fixture README. Separate synthetic edge cases from NCBI-derived acceptance.

Required tests additionally cover 3999/4000/4001-bp products; 5/6 total mismatches;
overlapping/duplicate HSPs and multiple sites per transcript; same-sequence
different-accession records; ambiguous/gapped alignments; missing/stale FASTA or
annotation; interrupted execution, missing outputs versus valid empty output;
and insufficient subject/HSP budgets. A known intended product must be recovered
or explicitly explain why the requested assessment is incomplete, not pass.
Duplicate HSPs realigning to one site must yield the same distinct products,
counts and intended recognition as a single supporting HSP; distinct sites and
conflicting alignments must remain distinguishable. Compact/full display modes
must agree on verdicts and preserve the profile's promised audit evidence.

Byte-bound fixture inputs need scoped checkout attributes and LF/CRLF coverage
through `scripts.test_tutorial_checkouts`. Use deterministic unit fixtures and
mock/local BLAST boundaries in normal CI. Actual BLAST+ acceptance records its
version and uses the same retained corpus and request; no live web dependency.
At implementation handoff run the focused policy/reader/resource/CLI/MCP tests,
protocol/discovery guards and `cargo check --locked`; quote exact filters/counts.
Record the tested SHA, lockfile and build profile. Native Windows file/rename
behavior and the full human snapshot remain distinct acceptance items.

## Open Inputs, Not Silent Defaults

- Glen's exact source-backed two-pair fixture/settings receipt, intended accession
  versions and whether expected recovery means listing, classification or both.
- The chosen species-wide snapshot source and owner-approved storage budget.
- Verified treatment of gaps, ambiguity and alignment ties. Unsupported cases
  remain incomplete; they do not block implementation of unambiguous fixtures.
- Pin the BLAST+ version for native acceptance. NCBI backend scoring and
  supplied-primer listing semantics are not independently established by the form.
- New fields must preserve legacy projections at every stored identity site.
  Exact type/schema choices follow the baseline tests; absent profile fields
  must not reinterpret reports, pending handoffs or panel selections. New resource
  reads are opt-in here; broad legacy-reader repair is a separate decision.

## Review And Execution Record

### Original Codex Plan (2026-11-16)

The accompanying [initial Claude review request](primer_blast_local_compatibility_claude_review_prompt_20261116.md)
retains the original brief. The original implementation order was:

1. Named pair policy and fail-first synthetic tests, with legacy invariance.
2. Prepared-cDNA reader and a shared direct/handoff/import search profile.
3. Versioned human RefSeq-mRNA preparation, tiny local intake before real data.
4. Existing CLI/MCP/agent integration and Glen's source-backed acceptance.

The original plan incorrectly described the terminal veto as a pre-pairing
discard and assumed existing site/product deduplication. It also proposed an
unnecessary 50000 subject-limit floor. Those statements are superseded by the
source audit and revised implementation plan above.

### Initial CLI Review Attempt: Unavailable

The local Claude CLI was invoked with `-p --permission-mode plan`, allowed tools
`Read Grep Glob`, and the concrete original plan/evidence above. The initial
sandboxed attempt returned no review and was terminated. A network-enabled,
five-minute-bounded retry exited with `401 OAuth access token is invalid`.
No completed critique was received from that invocation; no credential
inspection, login attempt or token repair was performed.

### Owner-Supplied Claude Review And Codex Assessment (2026-11-17)

The owner subsequently supplied Claude's read-only review. Its stated method
was source reading only, with no tests, downloads or independent verification
of NCBI defaults. Claude disclosed potentially offset line numbers in its local
checkout. Codex rechecked the cited paths in the unchanged `b62d21c9` checkout;
this is plan review, not implementation or verification acceptance.

| Finding | Codex disposition |
| --- | --- |
| 1: saved identities | Accepted. Real stored-record reserialization/recomputation tests come first, covering single/multi handoffs, summary and selected-report hashes. Qualify the prescription: omitted legacy serialization is mandatory; every Rust field need not literally be `Option`. Existing serialization is not changed. |
| 2: missing deduplication | Accepted for the inspected specificity path. Sorting is not site/product deduplication. Add opt-in consolidation with retained HSP provenance and unique-intended/count tests. Do not infer that every mismatch necessarily produces duplicate HSPs. |
| 3: incomplete window labeled complete | Confirmed in `primer_specificity_full_alignment` and the reader's end check. New-profile length-clamped admission and incomplete evidence are explicit; legacy behavior does not change silently. |
| 4: argument/filter paths | Confirmed and strengthened. Multi-import has exact hard-coded arguments; single import lacks equivalent profile validation. `src/engine.rs` confirms inherited project/default thresholds really can discard HSPs before pairing in both direct and imported runs. |
| 5: four reader dependencies | Accepted: prepare/bind the missing resource FASTA index, extend typed cache validation, resolve subject lengths/accessions and require exact resource resolution. Actual BLAST RefSeq identifiers need a native local fixture. |
| 6: legacy cDNA identity | Accepted. Gate sequence-resource reads and changed alignment semantics on the explicit new profile; a broad correction needs separately approved identity/migration handling. |
| 7: predicate and legacy flow | Accepted. The predicate decomposes per site; the biological difference is not the timing of pairing. Per-product explanation remains useful. Use full-alignment counts. Gap-column counting is not accepted as a verified NCBI rule. |
| 8: flags, compaction, ceilings | Accepted. Test intended recognition/summary/retention end to end and pin both ceilings to 4000 bp. Mere presence of the splice-variant option is not proof of behavior. |
| 9: diagonal extension alternative | Consider, not adopt as settled. It offers an anchored ungapped placement per HSP, not unique transcript placement or verified NCBI equivalence. Gaps, ties and shifted placements remain explicit decisions. |
| 10: RefSeq intake gaps | Confirmed. Stream into separate derived outputs, preserve source-byte identity and typed annotation/census evidence. Additional check: `genbank_record_chromosome` does prefer VERSION, but its accession/name fallbacks are insufficient for the new contract. |
| 11: subject-limit floor | Accepted. Keep validated N+1, record its difference from the public 50000 default, and avoid the unnecessary floor. |
| 12: seed limitation/oracle | Accept the missing-exact-7-mer counterexample and bounded ungapped fixture comparator. Do not promote a seed-existence bound into a BLAST detection guarantee; no full-corpus mismatch-enumeration guarantee is claimed. |

### Revised Minimal Commit Order

1. Slice 0 only: baseline stored-identity regression tests plus the unwired pure
   ungapped classifier and truth-table/decomposition tests. No protocol change.
   This is the smallest proposed first implementation commit.
2. Slice 2's resource reader/index/length/identifier/cache admission, gated behind
   an explicit internal opt-in until the public profile is wired. Preserve the
   legacy path and test no unintended fallback or compatible-resource selection.
3. Slices 1-2's new-profile alignment, site consolidation, verdict/retention/ceiling
   plumbing and resolved direct/handoff/import settings/validation. Add omitted-
   on-legacy metadata only after the baseline digest tests exist. Split these
   concerns into smaller reviewable commits if needed, never a half-enabled claim.
4. Slice 3's streaming tiny RefSeq intake and immutable resource preparation.
   Real human-data preparation waits for explicit size/location approval.
5. Slice 4's shared-interface and source-backed acceptance gates. Glen's
   exact fixture can be acquired in parallel; until supplied, those acceptance
   rows stay pending rather than being synthesized from the expected lengths.

In the initial planning session, only the plan and review-request documents
were added. No runtime/catalog,
release status, Rust dependency or existing report was changed. No Rust builds
or tests ran and no reference dataset was downloaded. The initial documentation
hygiene check reported four OK, two warnings and no failures: intentionally
untracked documents/artifacts and the manual plan-fidelity reminder. The
pre-existing SiLA review prompt and `output/` were left untouched. The 17 November
revision changes planning documentation only; it does not implement a finding.

The owner's 17 November request authorizes Slice 0 implementation, focused
verification and commit. It does not authorize later slices, substantial data
transfer, external execution, push or scientific acceptance.

## Slice 0 Implementation (2026-11-17)

- `src/engine/analysis/primer_specificity_identity_tests.rs` reads frozen stored
  JSON, recomputes single/multi handoff identities and summary digests, checks
  report/policy serialization hashes and verifies that the actual selected-report
  lookup remains available. Rebuilding a selection from its stored acceptance
  reproduces its ID/hashes. Negative controls alter command, alignment and
  selected-report evidence. The explicit ignored capture helper refuses an
  existing output directory; normal tests never regenerate expected records.
- `test_files/fixtures/primer_specificity_identity/` retains entirely synthetic
  records captured from `b62d21c971c207704b96396523eec0ce61659806` runtime sources
  plus the test-only capture harness, before including the classifier in the
  engine module tree. Raw records SHA-256:
  `5e3695c7364418cc6b38142168d92ae47dda5a4d8a62025ffc0c592eb2f7fba6`.
  The artificial standalone report ID is a fixture label; the regression freezes
  its serialization hash rather than claiming a newly assessed scientific ID.
- `src/engine/analysis/primer_blast_mismatch.rs` is a private unwired primitive,
  with no serde types, profile/request field, engine operation or adapter route.
  Complete ungapped count inputs must be geometrically possible; evidence
  completeness and the intended/off-target decision are caller preconditions,
  not conclusions of this helper. Both exclusion reasons survive independently.
  Its four tests include every 96 x 96 valid count combination for 20-base inputs.
- `.gitattributes` and `scripts.test_tutorial_checkouts` protect both frozen
  files in LF/CRLF checkouts and detect the missing-policy negative control.
  These are required byte-contract coverage, not unrelated tutorial changes.
- The reviewed plan and historical review request are included for a durable
  scope record. The unrelated SiLA review prompt and `output/` stay untouched.

### Verification Recipe And Boundaries

Use one task-owned target directory, `--locked --offline`, one job,
`CARGO_INCREMENTAL=0`, `CARGO_PROFILE_DEV_DEBUG=0` and
`CARGO_PROFILE_TEST_DEBUG=0`. The fixture capture ran before the classifier module
inclusion; normal regression tests must run after it. The existing 26-GB checkout
target is not modified or cleaned. Remove only the task-owned target afterward.

The local host is macOS arm64, `rustc 1.100.0-beta.13 (3b0202dcd 2026-11-07)`.
Unchanged `Cargo.lock` SHA-256:
`d880f70bcb3bb8ebd8db53e0ef7a8456a060c841a891f97f8d4a2eaf159569fe`.

Results at `b62d21c9` plus this slice's final source/test changes:

| Check | Result |
| --- | --- |
| Explicit baseline capture, before classifier module inclusion | 1 capture test passed; this generates synthetic test data, not acceptance evidence |
| `cargo test --locked --offline --no-default-features --lib -- primer_specificity_multi:: primer_blast_mismatch:: specificity_reference_ --test-threads=1` | 23 passed, zero failed, 1 intentionally ignored capture helper: 5 new identity tests, 4 new classifier tests and 14 existing guards |
| `cargo check --locked --offline -q` with default features | Passed after granting normal Swift cache access; initial sandbox failure was in apple-metal's Swift cache, not GENtle source |
| `python3 -m unittest scripts.test_tutorial_checkouts -q` | 100 passed, including the new LF/CRLF/missing-policy control |
| Scoped rustfmt, staged/unstaged `git diff --check`, lockfile hash | Passed; lockfile unchanged |
| Session-close hygiene | 4 OK, 2 warnings, zero failures; intentional staged records/unrelated untracked files and manual plan-fidelity reminder |

The headless test build retains three pre-existing dead-code warnings in
`genome_track_subscription_count`, `dna_cache_diagnostics_enabled` and startup
trace `measure`, plus the macOS `__eh_frame` compact-unwind warning. None was
changed to broaden this slice. Earlier focused/standalone runs overlap the final
23-test selection and are not additional tests. The store envelope itself is
not frozen against unrelated additive collections; bound records and the real
selection/hash paths are the regression boundary.

Native Windows, Linux, full-suite, real BLAST/RefSeq, CLI/MCP compatibility-profile
and scientific acceptance are not established by these tests. Later slices remain
unimplemented. The task-owned build target is cleaned after verification; existing
checkout build artifacts and unrelated files are not removed.
