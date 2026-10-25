# RNA Read Origin: Current Implementation And Residual Plan

2026-10-22 addition: independent loaded loci and direct FASTQ/gzip use a separate
portable `rna-reads panel-*` contract, with all qualifying member sets and
same-pass supplier MD5/compressed SHA-256 verification. See [tutorial 10.04](tutorial/10-04_targeted_nanopore_fastq_panel.md)
and [the protocol](protocol.md#multi-locus-fastq-screening). Rust/native/real-data
acceptance is pending. This does not repair or relabel legacy annotation-order
ties or the 64-posting chain cap; the new panel does not use those winner/cap
paths. The separate `39ac4e44` repair moves matched-target bookkeeping after
strand filtering. Its rerun and the [full-output regression comparisons](glen_next_verification_handoff.md#shared-sequence-ambiguity-and-legacy-output-regressions-2026-10-25)
remain external verification, not acceptance of legacy winner semantics.
Legacy FASTA intake now
refuses FASTQ rather than treating a quality `>` as a sequence header.

Source review: 2026-10-07. Multi-gene transcript targeting, origin classes and
strand diagnostics are implemented. Annotation-independent ROI capture is
still inactive in the interpretation path; a retained option is not execution
evidence.

## Implemented Contract

`InterpretRnaReads` retains default `single_gene` behavior and accepts
`origin_mode=multi_gene_sparse` with `target_gene_ids`. The engine collects
additional transcript lanes from the loaded annotation, merges their seed
space and builds exon/transition support indexes. It reports matched and
unmatched requested targets rather than treating unrepresented genes as absent.
This is local annotated-template evidence, not whole-transcriptome coverage.

The shared `rna-reads interpret` / `batch-map` paths expose these options.
Persisted reports, summaries, sample sheets and read exports retain origin mode,
target IDs, origin counts/reasons and strand diagnostics. The RNA Mapping GUI
inspects the same reports; filtering a display does not rerun classification.

The protocol's current class names are:

- `target_coherent`
- `target_partial_local_block`
- `roi_same_strand_local_block`
- `roi_reverse_strand_local_block`
- `tp_family_ambiguous`
- `background_likely`

`origin_confidence` and `strand_confidence` are numeric heuristic diagnostics,
not calibrated probabilities. Local-block class names do not demonstrate an
independent genomic ROI capture layer, and the legacy family-like ambiguity
label does not identify a biological family by itself.

Current source: `src/engine/analysis/rna_reads.rs` and
`crates/gentle-protocol/src/lib.rs`. Canonical
[CLI](cli.md), [protocol](protocol.md) and [GUI](gui.md) document the exported
fields. Multi-gene lane expansion, missing-target persistence, origin
classification and cross-strand tie handling have deterministic Rust regression
tests; this documentation review did not execute them or establish throughput.

## Genuinely Open Work

The independent-panel slice deliberately uses strict separate FASTQ requests,
not a new format variant on the legacy single-coordinate report. Full compressed
digests are verified during parsing instead of a head/tail fingerprint or an
upfront pass. Continuation is by verified file units, not record-count resume
into possibly replaced gzip. Potential duplicate source views are refused, not
merged automatically. Its `*_candidate` classes deliberately avoid "unique
isoform" or "proven chimera" language. Demo/show/cancel routes were added to
make every tutorial stage available to the inner agent without OS shell tools.
Separate repairs of legacy same-strand winner ties and the 64-posting cap
remain open; do not count the new
panel's deterministic tests as fixing or accepting that legacy path.

1. Annotation-independent genomic ROI seed capture. The current engine warns
   `ROI capture layer is planned but not active in phase-1 scoring` when
   `roi_seed_capture_enabled=true`. The option and sample-sheet field record
   the request, not successful capture. Defining postings for genomic/exonic/
   junction context requires a separate tested contract before an antisense or
   annotation-incomplete discovery claim.
2. A general portable sparse-posting/index resource, if a measured consumer
   requires it. Current combined template seed/support indexes are not proof
   of the originally proposed comprehensive `seed_bits -> postings[]` resource
   or CSR/CSC views. Do not rebuild working multi-gene targeting merely to match
   the old design sketch.
3. Representative background/calibration and same-input performance evidence.
   Preserve explicit competing-strand diagnostics; heuristic confidence must
   not become a measured probability or uniqueness claim.
4. Cohort presentation and evidence-to-curation follow-up are owned by the
   [Nanopore region plan](rna_seq_nanopore_cloning_regions_plan.md) and
   [roadmap](roadmap.md), not new adapter-local classification.

These changes must preserve default behavior and shared-engine parity.
Any future capture fixture must be synthetic or independently provenance-bound.
No online retrieval, automatic annotation mutation or full-transcriptome
quantification is implied.
