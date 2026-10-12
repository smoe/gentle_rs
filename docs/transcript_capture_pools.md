# Transcript Capture Pool Discovery

This is an offline, gene-agnostic **candidate discovery** workflow. It proposes
oligos against loaded annotated transcripts; it does not approve specificity,
predict amplification yield, create an order form or order anything.

## Biological Contract

A `sense_forward` oligo has the transcript-sense sequence and binds antisense
first-strand cDNA. Its annotation-derived retained interval extends from the
binding start towards the transcript's 3-prime end. An `antisense_reverse`
oligo is the reverse complement of its transcript binding site; its retained
interval extends from the annotated 5-prime end through the binding end.
These are substrate/reach descriptions, not statements about the physical
direction of a Nanopore read or guaranteed PCR products. Oligo(dT) can provide
the other priming end; no universal adapter is required merely to discover a
sense primer. The exact protocol still governs library preparation.

Target identity, window, primer budget and sharing permission are independent:

- A target can contain multiple annotated loci/transcript sets. Each source
  uses the existing `TranscriptAssayCoverageUniverse`: all annotated cDNA
  classes or a content-bound explicit transcript list. Unresolved or ambiguous
  requested identities fail rather than disappearing. Direct UniProt inventory
  joins are not supported in this first capture slice.
- `five_prime_utr` searches wholly before the annotated CDS start in spliced
  transcript coordinates. An unknown CDS boundary is not a zero-length UTR.
  Unknown/short windows remain in the coverage denominator with an explanation.
  There is no automatic coding-sequence fallback. Explicit `transcript_range`
  or `terminal_exon_start` policies are separate user choices.
- `max_primers: 1` does not become two when coverage is incomplete. Equal,
  nonempty `sharing_group` keys permit one candidate across targets with the
  same role, tail and stage set. They do not require a shared solution.
- Exact words are enumerated from every eligible window. No longest or
  canonical representative, global alignment, gene name, or paralog allowlist
  determines the binding-site search.
- Fixed oligos retain the full sequence, provenance, stages and
  `reorder_same_sequence`. They participate in interaction diagnostics, not in
  inferred target coverage. Equal sequences belonging to separate intents are
  not silently deduplicated into an order.

## Shared Shell

The following request shape assumes an annotated `locus` is already loaded;
feature `0` must belong to the desired transcript group. It is not a claim about
any particular human reference or an actual orderable primer.

```json
{
  "schema": "gentle.transcript_capture_pool_request.v1",
  "report_id": "capture_discovery",
  "cdna_synthesis": "oligo_dt",
  "coverage_policy": "require_all",
  "targets": [{
    "target_id": "my_transcript_set",
    "sources": [{"seq_id": "locus", "source_feature_id": 0}],
    "role": "sense_forward",
    "window": {"kind": "five_prime_utr"},
    "max_primers": 1,
    "stage_ids": ["pcr"]
  }],
  "fixed_oligos": [],
  "search": {
    "min_length_bp": 20,
    "max_length_bp": 24,
    "max_candidates_per_target": 12,
    "beam_width": 128,
    "internal_a_run_min_bp": 12
  }
}
```

```text
primers design-transcript-capture-pool @capture-request.json
primers list-reports
primers show-report capture_discovery
primers export-report capture_discovery capture-report.json
```

Use these through the GUI Shell or `gentle_cli shell 'COMMAND'`. MCP `op`,
JS/Lua operation entry points and Python's operation wrapper use the same
`{"DesignTranscriptCapturePool":{"request":...}}` operation. The inner agent
discovers the shell route and its fact annotations through introspection.
There is not yet a dedicated capture-pool GUI form or reviewed order adapter.

## Read The Report

`gentle.transcript_capture_pool.v1` stores the normalized request, request
digest, operation/run IDs, complete source-record digests, available genome
anchors and the existing coverage-universe resolution. It remains in the
`gentle.primer_design_reports.v1` store in its own default-empty map; old
terminal-exon reports are not repurposed or recalculated.

Each member retains its exact annotation identity, cDNA digest, strand, UTR
boundary, permitted window and diagnostic notes. Every retained candidate
lists **all exact occurrences in the requested transcripts in its declared
binding orientation**, including out-of-window and non-permitted targets.
Those extra occurrences are not coverage successes. This is not an exhaustive
off-target search, and does not inspect omitted transcripts or the other
orientation. A later specificity screen must consider those too.

For each occurrence, coordinates are 0-based and half-open. Source ranges stay
ascending within each exon on either strand and preserve transcript exon order.
The report records upstream/downstream omitted bases, retained length and
retained-sequence digest. Canonical retained sequences form equivalence groups:
differences upstream of a forward primer are lost, but downstream differences
remain informative. Identical ambiguity symbols do not prove equivalence.
Primer-supplied bases themselves are not independent template observations.

`coverage_satisfied` means every requested member has an exact eligible site in
the proposed pool. `require_all` with false coverage is a diagnostic incomplete
proposal, not a usable all-transcript design. `best_effort` permits presenting
partial coverage but does not alter the denominator or imply specificity.
`proposed_candidate_ids` are never automatic approval.

## Ranking And Limits

Candidate retention prefers broader exact coverage, then lower worst self/fixed
3-prime and general complementary runs, lower homopolymer runs, shorter
annealing segments and stable candidate IDs. Each target's retention quota
requires actual eligible coverage of that target; permission to share does
not let another target's candidates displace all of its own candidates.
The bounded pool search ranks by
uncovered members, number of distinct proposed oligos, worst complementary
runs, homopolymer length, omitted bases and stable candidate indices. It counts
a shared oligo once physically and against every target it covers. Pairwise
interactions apply only to overlapping declared stages and use full oligos,
including tails. Canonical comparisons reuse the existing engine heuristic;
IUPAC protocol oligos use conservative possible-match runs, not probabilities
or thermodynamic energies. Fixed/fixed interactions are also reported.

The Tm estimate and its assumptions are recorded explicitly. Without
`search.tm_range`, Tm is descriptive. An explicit `{ "min_c": 55, "max_c": 70 }`
filters annealing segments under that same model; it does not establish the
reaction's annealing temperature or impose an invented oligo(dT) constraint.

The bounds are 16 targets, 16 total new primers, 256 resolved transcript
members, 1,000,000 total cDNA bases, 5,000 bp per search window, 20,000 distinct
discovery candidates, 128 retained candidates overall and beam width at most
512. The retained occurrence audit stops before exceeding 50,000 bindings or
100 million retained bases to hash, including repeated sites. Requests
exceeding hard bounds fail explicitly. Retention and beam pruning
are reported; even a complete binding result is not a proof of global optimality.

Exact internal A-runs at the declared threshold are reported as possible
oligo(dT) priming sites. They are not measured products, and the report does not
model arbitrary A-rich priming, poly(A) length, RT falloff, or competition/yield.
Retained-length spread is descriptive. Synthetic tails and poly(A) bases are
excluded from retained sequence lengths and equivalence.

## Self- And Cross-Interaction Review

This is not a new need introduced by sense primers. The existing
`primers design-terminal-exon-rt-pool` routine already ranks reverse-primer
pools using variable-segment and complete-oligo self/cross-complementarity,
including the fixed adapter. Its `selected_pool_interactions` records every
selected pair. Use that routine for its terminal-exon reverse-primer geometry,
not as a substitute for a mixed-orientation capture request.

For `primers design-transcript-capture-pool`, inspect candidate
`self_complementary_run_bp` / `self_3prime_run_bp` and the report's
`interactions`, including fixed/fixed pairs. Supply the complete oligos and
their actual shared stages; an omitted universal primer or unchosen tail
cannot be assessed. Retrieve the unchanged result with
`primers show-report REPORT_ID`. These run lengths are inexpensive sequence
heuristics, not hairpin energies, dimer melting temperatures or proof of a
compatible multiplex.

Follow shortlisted pools with a thermodynamic assessment using a tool such as
[Primer3](https://primer3.org/manual.html). Its `ntthal` utility can evaluate
hairpins, self-dimers and cross-dimers, including end-anchored interactions;
`oligotm` evaluates primer/target duplex Tm. Check full oligos as well as
annealing cores, every co-present pair, and each alternative in its proposed
pool rather than mixing mutually exclusive candidates together. Record tool
version, parameter files, monovalent/divalent salt, total dNTP, oligo
concentration and the temperature used for reported free energies. Structure
Tm is not the reaction annealing temperature or a probability of forming an
amplifiable dimer. A negative structure search is not experimental validation.

The pool routes above do **not** automatically invoke this external
thermodynamic assessment. Preserve its outputs and assumptions separately;
do not relabel complementary-run counts as Primer3 results, replace the
source report, or infer order approval from either screen alone.

## This Experiment And Next Gates

For the requested reorder: PATZ1 and FUS should each have one `sense_forward`
target with a preferred `five_prime_utr` window. MDM2/MDM4 are real additional
targets; their role/window must be supplied rather than silently forced to 5'
primers. Matching sharing keys allow comparison with separate candidates.
Do not restore the removed PATZ1 terminal primers. Supply the retained order
sheet and exact protocol oligos, with their reaction-stage membership, before
treating pool diagnostics as representative of the experiment.

Deterministic implementation tests use synthetic data. A local, non-published
Ensembl 116 candidate review has also checked the original order's orientation
and compared PATZ1/FUS/MDM2/MDM4 candidates; this is not experimental acceptance
or order approval. In particular, `terminal_exon_start` does not mean 3-prime
UTR: its candidates can truncate coding sequence. Inspect retained geometry
and CDS boundaries separately from exact binding coverage. Remaining gates are
bound transcriptome and genomic specificity (including permitted multi-locus
targets), experimental protocol review, human selection, and reviewed
oligo-order handoff. Existing pair-specificity or terminal-pool uniqueness
approval must not be copied onto this single-primer workflow.
