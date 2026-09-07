# Reproducible regulatory-region comparison preparation

`scripts/prepare_regulatory_region_indexes.py` prepares an explicitly declared
set of promoter/upstream regions for comparison. Sequence extraction remains a
GENtle operation against one exact prepared reference genome; the companion
only orchestrates commands, combines exported FASTA records, builds indexes,
and records hashes.

The example manifest covers the three transcript/TSS models already used for
CD44, TGFB1, and SERPINE1:

```bash
python3 scripts/prepare_regulatory_region_indexes.py \
  --manifest docs/examples/regulatory_region_comparison/cd44_tgfb1_serpine1.json \
  --output /tmp/cd44-tgfb1-serpine1-regions
```

`--catalog` and `--cache-dir` may override host-local locations without editing
the manifest's genome identity. Both resolved paths and the catalog digest are
retained in the receipt.

The output directory must be empty. It receives:

- the GENtle project and one lossless FASTA export per region;
- `candidate_regions.fa` and `candidate_regions.json`; exact-identical windows
  remain mapped to every transcript but occur only once in the comparison
  database through explicit sequence-equivalence classes;
- a BLAST v5 candidate-region database;
- complete all-vs-all `megablast` and `blastn` tables;
- deterministic canonical-k-mer signatures and pairwise Jaccard/containment
  values at several word lengths;
- exact command stdout/stderr, hashes, tool identities, a receipt, and a sorted
  checksum inventory.

The source FASTA, JSON, k-mer evidence, and comparison tables are deterministic
for identical bound inputs and tool versions. BLAST database containers may
embed build metadata; they are treated as regenerable indexed artifacts whose
exact bytes are nevertheless hashed in each run receipt.

## Updating human or using another genome

Change only the manifest:

1. select a catalog `genome_id` that binds the intended assembly and annotation
   release;
2. choose a distinct `dataset_id`;
3. list exact gene/transcript models and strand-aware upstream/downstream spans;
4. point `catalog_path` and `cache_dir` at the matching GENtle resources.

If the declared genome is not already prepared, rerun with
`--prepare-genome`. That flag is an explicit authorization for GENtle to obtain
and prepare the catalog resources; without it the companion fails closed on a
missing cache. A release update must use a new output directory and dataset ID,
so evidence from different annotation releases cannot be silently overwritten.

## Interpretation boundaries

The candidate-region BLAST database answers **which declared candidate regions
resemble one another**. It is not a whole genome and cannot establish genomic
uniqueness. Use `regions homology-screen` against explicitly validated genomic
indexes for same-genome repetition and cross-species evidence.

`blastn` provides local nucleotide alignment. The shorter canonical k-mer rows
are a deliberately more permissive discovery aid for diverged regions, but
they discard order and spacing. They therefore do not establish an alignment,
orthology, TF binding, enhancer/promoter activity, or fragment sufficiency.
Combine them with GENtle's TF-score/module evidence and explicit reporter
contrasts before forming a biological hypothesis.

Verify a completed bundle from inside its output directory:

```bash
sha256sum -c checksums.sha256
```
