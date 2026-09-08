# Reproducible regulatory-region comparison preparation

`scripts/prepare_regulatory_region_indexes.py` prepares an explicitly declared
set of promoter/upstream regions for comparison. Inputs may be transcript/TSS
windows or portable `gentle.genomic_region_set.v1` documents. Sequence
extraction remains a GENtle operation against one exact prepared reference
genome; the companion only validates/imports typed inputs, orchestrates
commands, combines exported FASTA records, builds indexes, and records hashes.

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

## Ensembl Regulation promoters versus self-defined regions

First convert the exact promoter rows from one or more bound locus reports into
a canonical region set. The helper passes each complete normalized row and its
verified Ensembl source binding through GENtle; it does not reconstruct
coordinates itself:

```bash
python3 scripts/prepare_ensembl_promoter_region_set.py \
  --report /path/to/CD44.locus.report.json \
  --report /path/to/TGFB1.locus.report.json \
  --report /path/to/SERPINE1.locus.report.json \
  --set-id ensembl_regulation_promoters_2026_08 \
  --output /tmp/ensembl-promoter-regions
```

The command fails closed if Ensembl evidence is absent, its content identity
was not verified, or the overlap result was truncated. By default it selects
only rows whose exact provider `feature_type` is `promoter`. Repeat
`--feature-type` only when another provider class is intentionally in scope.

Export self-defined promoter/ROI selections from GENtle's region manager or
with `regions export`. Keep them in a separate set: a manual span must not be
relabeled as an Ensembl promoter merely because it overlaps one.

Declare both sets in the comparison manifest:

```json
{
  "schema": "gentle.regulatory_region_index_preparation.v1",
  "dataset_id": "ensembl_vs_selected_promoters_grch38_2026_08_v1",
  "genome": {
    "genome_id": "Human GRCh38 Ensembl 116",
    "catalog_path": "assets/genomes.json",
    "cache_dir": "assets/data/genomes",
    "expected_reference": {
      "taxon_id": 9606,
      "assembly_names": ["GRCh38"],
      "assembly_accessions": ["GCA_000001405.29"]
    }
  },
  "region_sets": [
    {
      "path": "ensembl-promoter-regions/ensembl_promoters.region_set.json",
      "comparison_class": "ensembl_regulation",
      "id_prefix": "ens_",
      "purposes": ["promoter_region"],
      "selection_methods": ["ensembl_regulatory_feature"]
    },
    {
      "path": "selected-promoters.region_set.json",
      "comparison_class": "self_defined",
      "id_prefix": "selected_",
      "purposes": ["promoter_region", "reporter_candidate"]
    }
  ],
  "indexing": {
    "blast": true,
    "blast_tasks": ["megablast", "blastn", "dc-megablast"],
    "canonical_kmer_lengths": [7, 11, 15, 21]
  }
}
```

Paths to region sets are relative to the comparison manifest. Before sequence
extraction the script imports every complete set through GENtle, thereby
checking its canonical identities and content digests. It also requires a
machine-checkable expected assembly/taxon binding. Filters that select no rows,
cross-assembly inputs, duplicate IDs after prefixing, and invalid half-open
coordinates are errors.

Explicit genomic intervals are exported in assembly-reference-forward
orientation. This is suitable for candidate-pool BLAST, which searches both
strands, and for canonical k-mers, which are strand-neutral. Original strand,
purpose, selection method, evidence, persistent colour, region identities,
source-set digest, and source-file digest remain in `candidate_regions.json`.

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
4. for region-set inputs, declare the expected assembly name/accession and
   taxon, then replace the source sets with ones exported for that assembly;
5. point `catalog_path` and `cache_dir` at the matching GENtle resources.

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
