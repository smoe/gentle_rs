# Data supporting genetic cloning with GENtle

| File               | Description                                                    |
|--------------------|----------------------------------------------------------------|
| bairoch.310        | EMBL-formatted description of proteases of the Rebase database as once downloaded from ftp.neb.com, currently not accessible. |
| JASPAR_2022.txt.gz | JASPAR 2022 release of transcription factors, downloaded from https://jaspar2022.genereg.net/download/data/2022/CORE/JASPAR2022_CORE_non-redundant_pfms_jaspar.txt |
| resources/nanopore_direct_cdna_kit14_adapters.fasta | External ONT direct-cDNA oligo signatures (SQK-LSK114, documented in ONT protocol DCS_9187_v114_revK_24Feb2026) for RNA-read concatemer/adaptor review. The shipped anchored poly(T) primer intentionally omits terminal `VN` from the matching sequence so the anchor does not overconstrain similarity search. |
| resources/nanopore_native_barcodes_nbd114_24_v14.fasta | External ONT Native Barcoding 24 v14 signature set (NB01-NB24 forward and reverse barcode sequences) for multiplex/ligation artifact review in RNA-read concatemer analysis. Added from the user-supplied barcode table referencing ONT Native Barcoding v14 documentation: https://nanoporetech.com/document/ligation-sequencing-gdna-native-barcoding-v14-sqk-nbd114-24 . |
| resources/tp73_dn_ena_transcripts.fasta | Curated TP73 DeltaN alpha/beta/gamma transcript FASTA for RNA-read transcript-catalog review. Retrieved from ENA accessions AY040827, AY040828, and AY040829 on 2026-05-01 using EBI dbfetch/ENA browser records; companion provenance and Ensembl comparison evaluation live in `assets/panels/tp73_dn_isoforms_v1.json`. |
| resources/tp73_delta_ex2_3_refseq_derived_transcripts.fasta | Derived TP73 D_{Ex2,3}Np73 alpha/beta transcript FASTA for RNA-read transcript-catalog review. Created on 2026-05-02 from the local RefSeq/NCBI GenBank TP73 locus fixture `test_files/tp73.ncbi.gb` by skipping transcript exons 2 and 3 in full-length alpha/beta source transcripts; companion provenance and derivation evaluation live in `assets/panels/tp73_delta_ex2_3_isoforms_v1.json`. |
| ../assets/publication_resources.json | Built-in catalog of publication-associated external datasets for the Rostock p73 multimodal profiling paper (DOI 10.3390/biom16010063), including ArrayExpress/BioStudies E-MTAB-14704, SRA/ENA-backed E-MTAB-15709 / PRJEB100610, and PRIDE PXD058816. `resources prepare-publication-dataset` writes per-dataset manifests/download scripts under `data/publication_resources` by default; large file download remains explicit via `--download-files`. |
| publication_resources/*/manifest.json, publication_resources/*/download_manifest.tsv, publication_resources/*/download.sh | Prepared no-byte-download manifests and executable curl scripts for the publication-associated datasets above. Generated with `cargo run --quiet --bin gentle_cli -- shell 'resources prepare-publication-dataset ACCESSION'`; regenerate after updating `../assets/publication_resources.json`. TSV manifests include declared sizes and optional archive MD5 checksums. |
| publication_resources/* downloaded payloads | Local-only external raw files fetched from the prepared manifests, intentionally ignored by git. For Clariom D CEL files, fetch with `cargo run --quiet --bin gentle_cli -- shell 'resources prepare-publication-dataset E-MTAB-14704 --categories metadata,raw_microarray --download-files'`. For the Rostock CUT&RUN FASTQs, use `--categories raw_cutrun_fastq --download-files` only when disk space is intentional. |
| resources/affymetrix/clariom_d_human_na36_hg38/README.md | Manual staging instructions for login-walled Thermo Fisher Clariom D Human na36 hg38 probeset/transcript support ZIPs. The ZIP payloads are intentionally ignored by git; `arrays probe-regions --platform Clariom_D_Human --dry-run` reports their expected local status under `annotation_source.vendor_support_files[]`. |

## Copyrights

The copyrights are with the respective authors of the data. Please reference

| Data   | Publication                                        |
|--------|----------------------------------------------------|
| Rebase | https://pmc.ncbi.nlm.nih.gov/articles/PMC102482/   |
| JASPAR | https://pmc.ncbi.nlm.nih.gov/articles/PMC10767809/ |
| Rostock p73 multimodal datasets | https://www.mdpi.com/2218-273X/16/1/63 |
