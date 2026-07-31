# Allele Hash Screen Fixture Provenance

Status: synthetic, hand-crafted, deterministic.

This fixture is a tiny FUS-like transcript-coordinate allele screen used by
`src/allele_hash_screen.rs` tests. It is not derived from real FUS sequence,
real patient genotypes, or public sequencing reads.

Files:

- `fus_transcripts.fa`: two synthetic transcript records tagged
  `gene_symbol:FUS`.
- `fus_variants.tsv`: two phased SNV rows over `FUS_TX1` in transcript/cDNA
  coordinates.
- `fus_reads.fastq`: small reads that exercise hap1, hap2, shared ambiguous,
  off-target, and too-short uninformative buckets.
- `read_allowlist.txt`: includes all committed read ids and exercises the
  allowlist loader.
- `salmon_unmapped_names.txt`: assigns four fixture reads, including
  `fus_hap2_alt_v1`, to a synthetic Salmon-unassigned cohort.
- `salmon_mappings.sam`: assigns the remaining short read to `FUS_TX1`, so the
  two Salmon selectors together reproduce the complete fixture cohort while
  still relying on `fus_reads.fastq` for sequences.
- `generate_fixtures.py`: deterministic regeneration script.

The key regression read is `fus_hap2_alt_v1`: with `k=9`, it is exactly one
alternate-allele k-mer spanning `v1`. It has zero support in the reference
hash, but one hap2-unique k-mer after haplotype-aware materialization. This is
the anti-reference-bias guard.

Regenerate:

```bash
python3 test_files/fixtures/allele_hash_screen/generate_fixtures.py
```

Used by:

- `allele_hash_screen::tests::fus_fixture_classifies_haplotype_and_ambiguous_reads`
- `allele_hash_screen::tests::salmon_selectors_reproduce_fixture_calls_and_provenance`
- `allele_hash_screen::tests::unphased_variant_labels_report_without_inventing_phase`
