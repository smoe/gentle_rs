# Rostock p73 CUT&RUN SRA Runs

This directory documents the SRA-backed CUT&RUN catalog entries in
`assets/cutrun.d/rostock_p73_sra.json`.

## Scope

The catalog entries point to the raw paired-end sequencing runs associated with
the Rostock p73 multimodal profiling study. GENtle does not bundle these raw
reads. `cutrun prepare DATASET_ID` delegates acquisition to the shared
read-acquisition/SRA Toolkit path via the entry's `reads_sra_accession`.

## Provenance

- Publication DOI: `10.3390/biom16010063`
- ArrayExpress/BioStudies accession: `E-MTAB-15709`
- SRA/ENA study: `ERP182066`
- BioProject: `PRJEB100610`
- Run accession retrieval:
  `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&term=E-MTAB-15709&retmode=json`
- RunInfo retrieval:
  `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&id=43727522,43727480,43726692,43726615,43725319,43724454,43724412,43724276,43724167,43723691,43723145,43720691&rettype=runinfo&retmode=csv`
- ENA FASTQ URL/checksum retrieval:
  `https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJEB100610&result=read_run&fields=run_accession,experiment_accession,sample_alias,library_layout,fastq_ftp,fastq_bytes,fastq_md5,submitted_ftp,submitted_bytes,submitted_md5&format=tsv`

The E-utilities query was checked on 2026-05-05 and returned twelve public
paired-end Illumina NextSeq 2000 runs. The ENA query was checked the same day
and is mirrored in `read_acquisition_manifest.tsv` with HTTPS FASTQ URLs, byte
counts, and MD5 checksums.

## Runs

| Run | Experiment | Sample label | Role |
| --- | --- | --- | --- |
| ERR15695857 | ERX15099975 | p73_TAp73a_SKMel29_p | p73 antibody, AdTAp73alpha |
| ERR15695855 | ERX15099973 | p73_DNp73b_SKMel29_p | p73 antibody, AdDNp73beta |
| ERR15695856 | ERX15099974 | p73_GFP_SKMel29_p | p73 antibody, AdGFP control |
| ERR15695854 | ERX15099972 | IgG_TAp73a_SKMel29_p | IgG control, AdTAp73alpha |
| ERR15695852 | ERX15099970 | IgG_DNp73b_SKMel29_p | IgG control, AdDNp73beta |
| ERR15695853 | ERX15099971 | IgG_GFP_SKMel29_p | IgG control, AdGFP |
| ERR15695860 | ERX15099978 | pos_TAp73a_SKMel29_p | positive-control CUT&RUN, AdTAp73alpha |
| ERR15695858 | ERX15099976 | pos_DNp73b_SKMel29_p | positive-control CUT&RUN, AdDNp73beta |
| ERR15695859 | ERX15099977 | pos_GFP_SKMel29_p | positive-control CUT&RUN, AdGFP |
| ERR15695851 | ERX15099969 | genome_TAp73a_SKMel29_p | genome/input-style control, AdTAp73alpha |
| ERR15695849 | ERX15099967 | genome_DNp73b_SKMel29_p | genome/input-style control, AdDNp73beta |
| ERR15695850 | ERX15099968 | genome_GFP_SKMel29_p | genome/input-style control, AdGFP |

## GENtle Usage

Examples:

```bash
gentle_cli shell 'cutrun list --filter E-MTAB-15709'
gentle_cli shell 'cutrun status rostock_p73_sra_err15695857_p73_tap73alpha'
gentle_cli shell 'cutrun prepare rostock_p73_sra_err15695857_p73_tap73alpha'
```

`cutrun prepare` can be large: each run has roughly 9-13 million paired spots,
and local FASTQ conversion can exceed the `.sra`/SRA-lite byte size. Keep full
downloads external to normal CI and release tests.

For a future 200 kbp TP73-focused internal fraction, derive the subset from
these SRA runs by mapping/downselecting reads against the TP73 GRCh38 envelope
and commit only a tiny provenance-documented FASTQ/FASTA fixture.
