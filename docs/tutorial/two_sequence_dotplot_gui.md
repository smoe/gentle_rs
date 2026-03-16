# Retrieve a cDNA and Compare It Against the Genomic Sequence (Dotplot Tutorial)

This tutorial shows a practical and biologically meaningful comparison:

1. fetch one TP73 cDNA sequence,
2. retrieve TP73 genomic context,
3. compare cDNA (query) against genomic DNA (reference) in the Dotplot map.

## Goal

Use this pair:

- cDNA/query: `NM_001126241.3`
- genomic/reference: TP73 locus extracted from `Human GRCh38 Ensembl 116`

This gives an exon-vs-intron-aware view where cDNA blocks align to genomic exon
regions.

## GUI Workflow

### Step 1: Fetch TP73 cDNA (query)

1. Open `File -> Fetch GenBank Accession...`
2. Accession: `NM_001126241.3`
3. `as_id`: `tp73_cdna`
4. Click `Fetch and Import`

### Step 2: Prepare and retrieve TP73 genomic context (reference)

1. Open `File -> Prepare Reference Genome...`
2. Select `Human GRCh38 Ensembl 116`
3. Click `Prepare`
4. Open `File -> Retrieve Genomic Sequence...`
5. Gene query: `TP73`
6. Occurrence: `1`
7. Output ID: `tp73_genomic`
8. Click `Retrieve`

### Step 3: Open query sequence window (cDNA view) and switch to Dotplot map

1. Open sequence window for `tp73_cdna`
2. Verify the cDNA is visible in the standard sequence/DNA map view
3. Switch primary map mode to `Dotplot map`

### Step 4: Configure pairwise dotplot

Set controls:

1. `Mode`: `pair_forward`
2. `reference_seq_id`: `tp73_genomic`
3. Keep `ref_start` / `ref_end` empty for first pass
4. Use a coarse first-pass profile for large genomic spans:
   - `word`: `11`
   - `step`: `40`
   - `mismatches`: `1`
5. Click `Compute dotplot`

If you want finer structure afterward, reduce `step` gradually (for example 30,
then 20) or narrow reference span with `ref_start` / `ref_end`.

### Step 5: Graphical inspection

What to expect in cDNA-vs-genomic mode:

- cDNA exons appear as separated diagonal blocks on genomic coordinates.
- Gaps between blocks correspond to intronic distances in genomic space.
- Strong continuous diagonal across the full range is not expected for
  spliced cDNA vs unspliced genomic reference.

Interaction notes:

- Hover shows `x(query)` and `y(reference)` coordinates.
- Click locks crosshair.
- In pair mode, selection sync is from the query/x axis.
- Right-click clears crosshair.

### Step 6: Optional orientation sanity check

Set mode to `pair_reverse_complement` and recompute.

For a same-orientation mapping, this usually weakens or removes the main
forward block structure.

## Equivalent CLI Commands

Fetch cDNA:

```bash
cargo run --bin gentle_cli -- shell 'genbank fetch NM_001126241.3 --as-id tp73_cdna'
```

Prepare genome and retrieve TP73 genomic locus:

```bash
cargo run --bin gentle_cli -- genomes prepare "Human GRCh38 Ensembl 116" --catalog assets/genomes.json --cache-dir data/genomes
cargo run --bin gentle_cli -- genomes extract-gene "Human GRCh38 Ensembl 116" TP73 --occurrence 1 --output-id tp73_genomic --catalog assets/genomes.json --cache-dir data/genomes
```

Compute pairwise dotplot (coarse pass):

```bash
cargo run --bin gentle_cli -- shell 'dotplot compute tp73_cdna --reference-seq tp73_genomic --mode pair_forward --word-size 11 --step 40 --max-mismatches 1 --id tp73_cdna_vs_genomic'
```

Inspect result:

```bash
cargo run --bin gentle_cli -- shell 'dotplot show tp73_cdna_vs_genomic'
```

## Notes

- Pair modes require `reference_seq_id`.
- If `ref_start` / `ref_end` are omitted, full reference span is used.
- For long genomic references, start with coarse `step` and then refine.
