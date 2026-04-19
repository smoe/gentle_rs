# Ensembl transcript catalogs

Naming scheme:
- root: data/transcriptomes/ensembl/release-116/
- coding transcripts: Homo_sapiens.GRCh38.cdna.all.fa.gz
- noncoding transcripts: Homo_sapiens.GRCh38.ncrna.fa.gz

Rationale:
- keep transcript catalogs separated from genome-cache installs under data/genomes
- keep provider and release explicit in the path
- keep upstream filenames unchanged once inside a provider/release folder

Source URLs:
- https://ftp.ensembl.org/pub/release-116/vertebrates/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
- https://ftp.ensembl.org/pub/release-116/vertebrates/fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz

Primary GENtle use:
- pass either file through rna-reads inspect-concatemers --transcript-fasta
- for broad partner-gene decomposition, start with the coding catalog and add the ncrna catalog when antisense/noncoding partners matter
