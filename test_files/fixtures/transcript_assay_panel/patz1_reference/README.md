# Authentic human PATZ1 reference snapshot

Public reference data, not a synthetic PATZ1-like sequence. Used by tutorial
04.08, `scripts/prepare_real_patz1_tutorial.py`, and the Rust
`real_patz1_sources_preserve_versions_and_verified_negative_strand_bases` test.
`manifest.json` pins exact bytes; `.gitattributes` preserves LF on Windows.

## Identity and limits

- Human PATZ1: Ensembl `ENSG00000100105.20`, NCBI Gene `23598`.
- GRCh38 chromosome 22 / `NC_000022.11`, inclusive `31325804..31346605`,
  minus strand, 20,802 genomic bp.
- Ensembl 116: 13 versioned transcript records, including noncoding and
  retained-intron models, not just the canonical coding model.
- RefSeq: `NM_014323.3`, `NM_032050.2`, `NM_032051.2`, `NM_032052.2`.
- The Ensembl gene-oriented sequence must be exactly the reverse complement
  of the independently fetched RefSeq genomic interval (checked in Rust).
- The NCBI response does not declare a release identifier. Provenance therefore
  uses accession plus retrieval date, not an inferred release. This is not a
  complete GenBank archive search. GenBank format is not a third annotation vote.
- No expression, array measurements, reference-wide specificity or experimental
  assay validation is implied by these reference annotations.

## Retrieval and recreation

Retrieved 2026-09-19. Exact URLs and transformations:

1. [Ensembl lookup](https://rest.ensembl.org/lookup/id/ENSG00000100105?expand=1;content-type=application/json)
   and [release information](https://rest.ensembl.org/info/data?content-type=application/json)
   (`{"releases":[116]}`). `ensembl_entry.json` is the existing GENtle importer's
   result, retaining raw lookup/sequence responses and their URLs. Only the
   import timestamp is normalized for storage:

   ```bash
   gentle_cli --state fetch.json op '{"FetchEnsemblGene":{"query":"ENSG00000100105","species":"homo_sapiens","assembly":"GRCh38","flank_5prime_bp":0,"flank_3prime_bp":0,"entry_id":"patz1_ensembl_116"}}'
   jq '.metadata.ensembl_gene_entries.entries.patz1_ensembl_116 | .imported_at_unix_ms = 0' fetch.json > ensembl_entry.json
   ```

2. [Ensembl GFF3](https://rest.ensembl.org/overlap/region/human/22:31325804-31346605?feature=gene;feature=transcript;feature=exon;feature=cds;content-type=text/x-gff3)
   is retained unmodified as `ensembl.gff3`. Neighbor records remain in the raw
   file; GENtle selects the explicit `gene:ENSG00000100105` parent.
3. [RefSeq chromosome GFF3](https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id=NC_000022.11&db=nuccore&report=gff3&retmode=text&from=31325804&to=31346605)
   returns the whole chromosome despite `from`/`to`. The manifest retains its
   original SHA-256. Keep headers, source-declared chromosome mapping and exact
   PATZ1 rows without changing any retained row:

   ```bash
   awk -F '\t' '/^#/ || $3 == "region" || $9 ~ /(^|;)gene=PATZ1(;|$)/' chromosome.gff3 > refseq.gff3
   ```

4. [RefSeq genomic FASTA](https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_000022.11&rettype=fasta&retmode=text&seq_start=31325804&seq_stop=31346605&strand=1)
   is retained unmodified as `refseq_genome.fasta`, on the genomic plus strand,
   including NCBI's final blank line. A file-specific `.gitattributes` exception
   excludes only that EOF whitespace warning; the manifest still checks all bytes.

Live services can change releases, annotation or serialization. Changed hashes
mean a new snapshot requiring review, not permission to silently refresh expected
values. Offline replay uses only committed data. No private inputs or invented
PSR/JUC measurements are included. Synthetic legacy fixtures elsewhere are for
regression tests, not a substitute for these PATZ1 annotations.
