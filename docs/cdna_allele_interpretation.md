# Conditional Isoform-Bound cDNA Consequences

Use this after [cDNA comparison](cdna_allele_comparison.md) and
[C1 model derivation](coding_transcript_models.md), not as a new variant caller.
C2-C4 consume authenticated snapshots and preserve upstream counts, statuses,
quality policy, reference and local phase. No fixed genes, samples or barcodes.
Glen's [scoped return at `5d310958`](glen_next_verification_handoff.md#glen-return-at-5d310958)
verifies the repaired C2-C4 Linux kernel/service, actual JS/Lua/MCP and CLI/Python
routes, plus historical comparison-byte preservation. Same-source Windows and macOS
compilation/workspace tests and aggregate CI also passed; full Linux and optional
GUI/audit jobs skipped. The return distinguishes remaining host/data checks from
scientific or release acceptance.

## Read-Only Routes

```sh
gentle_cli allele-hash-screen --interpret-cdna @interpretation.request.json
gentle_cli shell 'rna-reads allele-hash-screen --interpret-cdna @interpretation.request.json'
gentle_cli allele-hash-screen --interpret-cdna @interpretation.request.json --format tsv
gentle_cli allele-hash-screen --interpret-cdna @interpretation.request.json --format summary
```

Shell/GUI Shell, generic MCP `op`, JS/Lua and Python wrappers share
`InterpretCdnaAlleles { request }`. JSON is at
`result.cdna_allele_interpretation`; text formats project that report. There is
no dedicated GUI form or publication path. No loaded project, state save, undo
checkpoint or output file is required/created. This is an expensive observation
with cooperative cancellation. MCP's existing explicit confirmation precedes I/O.

## Request And Authentication

Strict schema `gentle.cdna_allele_interpretation_request.v1`:

```json
{
  "schema": "gentle.cdna_allele_interpretation_request.v1",
  "genes": [{
    "comparison_report": {"path": "comparison.json", "sha256": "sha256:REPLACE", "bytes": 1},
    "coding_models": {"path": "coding.json", "sha256": "sha256:REPLACE", "bytes": 1}
  }],
  "linkage": "none",
  "known_variants": null,
  "policy": {},
  "budgets": {},
  "format": "json"
}
```

These bindings are placeholders, not runnable hashes. Supply actual raw bytes.
Accept 1-64 distinct `(source_id,gene_id)` pairs with identical inherited
comparison policies. Bare comparison reports and the documented CLI envelope
are supported; `inputs[].comparison_form` records which was used. Verify raw
file sizes/hashes, both comparison content hashes and typed request hashes.
Reauthenticate and rederive C1 against its original bound files; compare its
content hash. Cross-bind reference/sequence, assembly, gene, normalized
annotation and every transcript's strand/exons. No repaired/guessed reference.
Local-file and symlink rules remain strict. Declared parent-genome identity is
cross-bound, not proof obtained by rereading a whole genome.

Keep the identity namespaces separate: C1 `selector.source_id` must equal its
normalized annotation's `reference_id`, whereas the comparison's `source_id`
names its saved-analysis source. They need not have the same spelling. The
shared raw annotation, exact gene and reference/model checks above link them;
C1 rederivation and C3's saved-analysis/request-hash verification remain strict.
The synthetic rehearsal deliberately keeps `synthetic_reference` and
`invented_source` distinct. Its earlier `2f1a3062` preparation incorrectly used
the latter for C1 and was refused before consequence/VM/CLI acceptance; see
the [Glen return](glen_next_verification_handoff.md#glen-return-at-2f1a3062).

Default policy is nuclear `translation_table: 1`,
`minimum_codon_haplotype_units: null` (inherit comparison threshold) and
`maximum_out_of_model_aligned_bases: 0`. A conflicting support threshold refuses;
quality/callability is never rescored or silently relaxed.

## C2: Consequences And Same-Unit Combinations

Per gene retain the reference, inherited policy, full sample/status/phase reports
and warnings. Site rows copy counts and distinct-unit denominators, conditional
fractions and status summaries. Project every reported alternate, including
below-threshold observations, into every supplied model. Region is CDS, 5'/3'
UTR, noncoding exon, unavailable coding geometry or `not_present`; absence from
an exon never becomes an intronic/upstream coordinate.

Versioned, sequence-verified coding models supply conditional codon changes.
Missing/inconsistent/invalid accessions withhold HGVS, not usable sequence
geometry. Incomplete CDS withholds coding notation/protein effects with explicit
reasons. Transcript end incompleteness alone is not CDS incompleteness. True
noncoding models use `n.`; undefined CDS is unavailable, not noncoding by default.

Group reported candidates within the same codon or adjacent transcript positions.
Every complete definite **same-unit, same-sample** base pattern is retained,
including reference/reference; incomplete units remain separately visible even
if no complete pattern exists. Inherited threshold marks support, never selects
a winner. TCC can therefore report unchanged Ser, Tyr, synonymous Ser and Ter
in one group. No phase-block or cross-sample composition fills missing calls.
Single-site consequences assume reference elsewhere; `codon_context_unresolved`
flags missing joint evidence, while `see_codon_haplotype` points to actual patterns.

Local relations copy block IDs, fingerprints, pattern alleles and phase methods.
`same_block_candidate_patterns`/`same_block_unresolved`, `different_blocks` and
`not_in_any_block` describe local RNA evidence, **not chromosomal cis/trans**.
H1/H2 labels are never compared across samples. Disconnected blocks stay separate.

### HGVS Scope

Rules checked 2026-11-12 against official HGVS **21.1.5 (2026-09-17)**:
[DNA delins](https://hgvs-nomenclature.org/stable/recommendations/DNA/delins/),
[protein substitution](https://hgvs-nomenclature.org/stable/recommendations/protein/substitution/),
[extension](https://hgvs-nomenclature.org/stable/recommendations/protein/extension/).
Supported output includes versioned `c.N`, exonic `c.-N`/`c.*N`, `n.N`,
substitution, delins and allele brackets, with three-letter protein names and
predicted parentheses: `p.(Leu145Arg)`, `p.(Leu145=)`, `p.(Trp274Ter)`,
`p.(Met1?)`, `p.(Ter394=)` and bounded `p.(Ter394GlnextTerN|?)`.

Contrary to the original Claude proposal, transcript-consecutive changes and
two changes separated by one base in one codon can be transcript delins even
across a genomic intron. The rule applies to the named reference, not genomic
adjacency. The former stop is position 1 of a stop-loss extension; scan only the
authenticated downstream in-frame transcript. This is not a genome-wide HGVS
normalizer or a claim of alternative initiation. Multi-codon changes retain
per-codon consequences but withhold combined protein notation by default.
No SO accessions are asserted; only the documented consequence-term names.

## C3: Conditional Variant-To-Isoform Evidence

`linkage: same_unit` reopens the exact saved analyses/basecalls through the
comparison's shared lowering. The normalized kernel request hash and all input
bindings must match the original report, otherwise `linkage_input_drift` refuses.
Request identity uses the comparison kernel's shared order normalization before
hashing; reordered samples or reverse-strand block observations do not change
evidence identity. Raw input bindings remain independently strict.
No read-name join, variant recall or transitive splice assembly occurs.

A supporting unit has all query views assessed. Within **every** eligible
hypothesis, the same qualified allele must coexist with compatible full alignment
geometry. Every junction must be a consecutive model intron, and every aligned
base must lie in model exons (unless the explicit bounded policy allows extras).
Off-locus/foreign-contig blocks cannot disappear through variant-quality clipping.
Retained-intron bases remain conflicting evidence.

Report compatible units, units excluding at least one selected competitor and
unique-model units separately, with excluded-transcript unit lists. Uniqueness
is only within supplied models containing the site. Partial compatibility is not
identification of a full-length/expressed isoform. Default `linkage: none` does
not reopen the original evidence and reports no invented linkage.

## C4: Frozen Labels, Not Clinical Interpretation

Optional bound `gentle.known_variant_manifest.v1` has `source_name`,
`source_release`, `retrieved_utc` (UTC calendar timestamp), `license`,
`attribution`, `assembly`, and `entries` with
`{contig,position_1based,ref,alt,identifier}`. Require exact assembly and canonical
nonempty ACGT alleles; duplicate/invalid entries and clinical/unknown fields
refuse. Match exact genomic-forward alleles only. Composite labels require a
supported same-unit pattern; single-site labels do not upgrade support. Unrelated
indels may be present in a mixed manifest but this substitution-only caller does
not match/call them. No contig aliases, shifting, fuzzy HGVS or live lookups.
Labels retain source/release/license/attribution; they are not confirmation.

## Determinism And Resource Limits

JSON schema `gentle.cdna_allele_interpretation.v1` binds typed request, raw inputs,
HGVS rules and the existing genetic-code asset. `content_sha256` hashes compact
canonical JSON with the empty hash field omitted. TSV and Markdown summary are
pure functions of saved JSON, escaping embedded control characters. TSV includes
site, codon-group and combination records with full JSON details; JSON remains
the full contract. Identical inputs have deterministic output, without timestamps.

All genes share cumulative limits: input 128,000,000 raw bytes, work 100,000,000
units, output default 32,000,000 (hard ceiling 64,000,000) bytes. Each reread is
charged: report/model files, C1 source verification, locus reloading, optional C3
saved evidence and frozen manifest. Existing per-file/C1/projection ceilings
still apply. Output accumulation and final report/projections are bounded. Work
includes C1 rederivation, shared lowering and kernel loops, not an RSS/time SLA.
Cancellation or a limit refuses without partial publication or project changes.

Tests invent all GTF, sequences, saved alignments, labels and samples in tempdirs;
no new real/private fixture is committed. Fast LF/CRLF checkout coverage exercises
raw bindings; optional actual-process checks need
`GENTLE_CDNA_INTERPRETATION_CLI`. Skips never count as acceptance. The
Python-client replay uses that same explicit executable and the shared shell,
without a Cargo fallback. Mocked argument/JSON/text checks are not runtime proof.
The separate real TP53 receipt remains unrun: Glen's verification was synthetic,
not authorization for a new biological data run. Any later authorized receipt
must use already bound committed evidence; no private sequencing run or
other-gene analysis is implied.

No DNA genotype/zygosity, somatic/germline origin, protein expression, copy number,
NMD or pathogenicity claim follows. OD-I1-I5 remain at their conservative defaults
in the [plan](cdna_allele_interpretation_plan.md).
