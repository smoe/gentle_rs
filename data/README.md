# Data supporting genetic cloning with GENtle

| File               | Description                                                    |
|--------------------|----------------------------------------------------------------|
| bairoch.310        | EMBL-formatted description of proteases of the Rebase database as once downloaded from ftp.neb.com, currently not accessible. |
| JASPAR_2022.txt.gz | JASPAR 2022 release of transcription factors, downloaded from https://jaspar2022.genereg.net/download/data/2022/CORE/JASPAR2022_CORE_non-redundant_pfms_jaspar.txt |
| resources/nanopore_direct_cdna_kit14_adapters.fasta | External ONT direct-cDNA oligo signatures (SQK-LSK114, documented in ONT protocol DCS_9187_v114_revK_24Feb2026) for RNA-read concatemer/adaptor review. The shipped anchored poly(T) primer intentionally omits terminal `VN` from the matching sequence so the anchor does not overconstrain similarity search. |
| resources/nanopore_native_barcodes_nbd114_24_v14.fasta | External ONT Native Barcoding 24 v14 signature set (NB01-NB24 forward and reverse barcode sequences) for multiplex/ligation artifact review in RNA-read concatemer analysis. Added from the user-supplied barcode table referencing ONT Native Barcoding v14 documentation: https://nanoporetech.com/document/ligation-sequencing-gdna-native-barcoding-v14-sqk-nbd114-24 . |

## Copyrights

The copyrights are with the respective authors of the data. Please reference

| Data   | Publication                                        |
|--------|----------------------------------------------------|
| Rebase | https://pmc.ncbi.nlm.nih.gov/articles/PMC102482/   |
| JASPAR | https://pmc.ncbi.nlm.nih.gov/articles/PMC10767809/ |
