# Read-Only Claude Review: Local Primer-BLAST Compatibility

Historical initial review request, preserved for traceability. On 2026-11-17 the
owner supplied Claude's response, which Codex assessed against the source. The
plan's review-disposition section and revised minimal order now supersede the
original order/assumptions below, including the 50000 subject-limit floor and the
claim that legacy terminal flags discard sites before pairing. This file is not
a record of an additional successful CLI consultation.

Repo: `/Users/u005069/.codex/worktrees/47dd/gentle_rs`.
Baseline inspected by Codex: `b62d21c9` (local main/origin/main).

This is a plan review, not implementation. Use only Read/Grep/Glob. Do not run
shell commands, builds, tests, downloads, write tools, Git changes or contact
another agent/session. If evidence cannot be verified with those tools, say so.

## Task And Glen's Clarification

We need local checking of supplied human primer pairs, not NCBI's entire service
or a redesign workflow. Glen requests a pinned human RefSeq-mRNA snapshot
including XM records, a named Primer-BLAST-compatible profile, full-primer
alignment and per-product mismatch decisions. Specifically:

1. retain and fully align candidate primer sites;
2. construct F-R, F-F and R-R products through 4 kb;
3. only then apply the per-primer mismatch rule to each product;
4. do not merely loosen GENtle's old terminal-hit filter.

The original pair should recover PATZ1 72 bp and GFI1 2751 bp; the replacement
PATZ1 126 bp and TBCE 167 bp. Exact oligos/accession versions/source output have
not been provided in this request, so these are proposed acceptance cases, not
verified fixture facts. Glen will review the plan too.

## Original Codex Plan

Read the complete original plan in
`docs/primer_blast_local_compatibility_plan_20261116.md`, especially its current-
source audit, biological contract, completeness tension and four slices:

1. Opt-in versioned pair policy and fail-first synthetic tests. Legacy omitted
   policy/serialization/IDs remain unchanged. Apply the product predicate only
   after retaining/aligning sites and pairing; no scientific computation in an
   adapter.
2. Extend existing prepared-FASTA reader for sequence-only cDNA manifests;
   reuse the aligner and per-interpretation validation/cache. Resolve one search
   profile for direct/handoff/import, requiring complete-primer evidence and
   testing partial HSPs, transcript ends and equal-score/gap/ambiguity handling.
3. Human NM/XM taxon9606 mRNA snapshot, immutable/source-bound and annotation-
   backed, reusing sequence-resource preparation. Local tiny intake first;
   approval of size/location before any large download. Assembly RNA alone is
   not assumed to cover the species-wide collection. Preserve accession versions
   and GeneIDs; unknown metadata is not invented.
4. Existing CLI/MCP/agent/handoff inspection integration plus synthetic and
   Glen-supplied sourced acceptance fixtures, direct/import equivalence and
   unchanged aggregate readiness. No redesign/UI/report-builder scope creep.

Product rule: ignore when either site's total mismatches >=6; otherwise exclude
as sufficiently discriminating when either individual site has BOTH total>=2
AND last-five>=2. Remaining unintended products are relevant. Never substitute
sum-of-two-primers mismatch arithmetic or the legacy terminal veto.

Compatibility profile pins public defaults observed 2026-11-16: E-value30000,
wordsize7, ungapped candidate BLAST, 4000-bp target limit, the 2/2/5 and 6 rules.
Codex verified raw public NCBI form values, including website subject cap50000.
For local coverage retain at least max(50000, validated subject count+1), naming
that deviation. An explicitly lower cap cannot pass. Subject-cap coverage is
not a proof that heuristic BLAST enumerates every possible binding site.
Exact NCBI backend/global alignment scoring equivalence is not established.

## Repository Evidence To Check

- `crates/gentle-protocol/src/lib.rs`: `PrimerSpecificityPolicy`, full-alignment
  modes/scoring and legacy deserialization defaults.
- `src/engine/ops/operation_handlers.rs`:
  `primer_specificity_full_alignment`, `_from_blast_strings`,
  `primer_specificity_hit_from_blast`, `primer_specificity_collect_amplicons_*`,
  `primer_specificity_finalize_amplicons`,
  `primer_specificity_search_completeness_for_commands`, handoff args/import.
- `src/genomes.rs`: aligned HSP fields/parser;
  `prepare_blast_sequence_resource_with_options`, GenBank conversion/annotations,
  `read_prepared_sequence_region` and resource/genome manifest distinction;
  direct BLAST argument construction. Aligned HSP support already exists on main;
  old commit `4f741797` was not available locally and was not reviewed.
- `src/genomes/prepared_sequence_reader.rs` and its tests;
  `src/engine/analysis/primer_specificity_reader_tests.rs`.
- `src/engine/analysis/primer_specificity_multi.rs` and
  `primer_specificity_references.rs`: exact command/completeness validation,
  bound references, old report/selection semantics.
- `src/engine.rs`: current report.v4, handoff.v1 and policy.v2 schema constants.
- `docs/primer_blast_comparison.md`, `docs/cli/primer-specificity.md`,
  `docs/protocol.md`, DEC-029/DEC-049 in `docs/decisions.md`,
  `docs/biological_extension_guide.md`, cross-platform rules in `docs/testing.md`.

External primary sources linked in the plan: NCBI Primer-BLAST form; Ye et al.
2012 DOI10.1186/1471-2105-13-134; RefSeq FAQ; NLM KA-03450. You have no browser
or shell in this review: distinguish Codex's external checks from facts you
independently verify. Do not turn lack of external access into an endorsement.

## Requested Output

Give concise severity-ordered findings with paths/evidence, then a revised
minimal implementation order. In particular:

- Does the new product predicate preserve sites the legacy filter discards?
- Does the cDNA reader gap really reduce to the separate manifest/index path?
- Are the complete-HSP shortcut, window clipping, tied/gapped alignments and
  candidate deduplication covered without promising identical server output?
- Does max(50000,N+1) preserve existing completeness while honestly documenting
  the departure from the web default? Are heuristic limitations clear enough?
- Would profile fields/schema changes break old digests or saved handoffs?
- Is the RefSeq scope genuine human mRNA including predicted models, rather
  than assembly-only RNA, RefSeq Select, or Ensembl masquerading as RefSeq?
- What is the smallest independently testable first commit, and which facts
  must come from Glen/owner before real-case acceptance?

Do not propose unrelated refactors or invent missing primers, accessions,
expected mismatches, data sizes, native-platform receipts or biological approval.
