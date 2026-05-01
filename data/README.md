# Data supporting genetic cloning with GENtle

| File               | Description                                                    |
|--------------------|----------------------------------------------------------------|
| bairoch.310        | EMBL-formatted description of proteases of the Rebase database as once downloaded from ftp.neb.com, currently not accessible. |
| JASPAR_2022.txt.gz | JASPAR 2022 release of transcription factors, downloaded from https://jaspar2022.genereg.net/download/data/2022/CORE/JASPAR2022_CORE_non-redundant_pfms_jaspar.txt |
| resources/nanopore_direct_cdna_kit14_adapters.fasta | External ONT direct-cDNA oligo signatures (SQK-LSK114, documented in ONT protocol DCS_9187_v114_revK_24Feb2026) for RNA-read concatemer/adaptor review. The shipped anchored poly(T) primer intentionally omits terminal `VN` from the matching sequence so the anchor does not overconstrain similarity search. |
| resources/nanopore_native_barcodes_nbd114_24_v14.fasta | External ONT Native Barcoding 24 v14 signature set (NB01-NB24 forward and reverse barcode sequences) for multiplex/ligation artifact review in RNA-read concatemer analysis. Added from the user-supplied barcode table referencing ONT Native Barcoding v14 documentation: https://nanoporetech.com/document/ligation-sequencing-gdna-native-barcoding-v14-sqk-nbd114-24 . |
| resources/tp73_dn_ena_transcripts.fasta | Curated TP73 DeltaN alpha/beta/gamma transcript FASTA for RNA-read transcript-catalog review. Retrieved from ENA accessions AY040827, AY040828, and AY040829 on 2026-05-01 using EBI dbfetch/ENA browser records; companion provenance and Ensembl comparison evaluation live in `assets/panels/tp73_dn_isoforms_v1.json`. |
| resources/tp73_delta_ex2_3_refseq_derived_transcripts.fasta | Derived TP73 D_{Ex2,3}Np73 alpha/beta transcript FASTA for RNA-read transcript-catalog review. Created on 2026-05-02 from the local RefSeq/NCBI GenBank TP73 locus fixture `test_files/tp73.ncbi.gb` by skipping transcript exons 2 and 3 in full-length alpha/beta source transcripts; companion provenance and derivation evaluation live in `assets/panels/tp73_delta_ex2_3_isoforms_v1.json`. |

## Copyrights

The copyrights are with the respective authors of the data. Please reference

| Data   | Publication                                        |
|--------|----------------------------------------------------|
| Rebase | https://pmc.ncbi.nlm.nih.gov/articles/PMC102482/   |
| JASPAR | https://pmc.ncbi.nlm.nih.gov/articles/PMC10767809/ |
