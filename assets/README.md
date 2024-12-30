# GENtle Assets

This folder contains a series of accessory files that are meant to be available
at the application's start or run time. These have mostly been carried over from the C++ version of GENtle.

amino_acids.json::
blank.mysql::
blank.sqlite3::
codon_catalog.csv::
codon_tables.json::
commonvectors.db::
dna_ladders.json::
enzymes.json::
genomes.json::
helper_genomes.json::
jaspar.motifs.json:: The JASPAR database transformed into a JSON format that is meant to be mostly compatible with the JSON format offered by the JASPAR project: ```gzip -dc data/JASPAR_2022.txt.gz | perl scripts/pfm2json.pl | jq --compact-output > assets/jaspar_2022.json```
ncoils.matrix::
translations.csv::

Data notes:

- `enzymes.json` is restriction-enzyme data (REBASE-derived snapshot).
- `dna_ladders.json` is the built-in DNA-ladder catalog used for pool gel
  ladder auto-selection and rendering.
- `jaspar.motifs.json` is a built-in JASPAR CORE motif snapshot (currently
  generated from the 2026 non-redundant JASPAR-format release).
- `genomes.json` is the default reference-genome catalog.
  Entries may use explicit `sequence_remote`/`annotations_remote` URLs or
  `ncbi_assembly_accession` + `ncbi_assembly_name` for direct NCBI GenBank/RefSeq
  FTP-derived sources.
- `helper_genomes.json` is a curated starter catalog for lab helper systems
  (plasmid/lenti/adeno/AAV plus yeast/E. coli host references), intended to be
  copied and edited for local vector inventories.
  Entries can use explicit local/remote URLs or `genbank_accession` for
  NCBI EFetch-derived FASTA + GenBank annotation during prepare/index.

# "Finder" icon (seen in OS' task bar)

- icon.png: PNG with icon as used at startup by application
- icon.icns: Icon-set for MacOS

## Script to transforms the .png icon to an icon set for the MacOS
```[bash]
mkdir MyIcon.iconset
sips -z 16 16     icon.png --out MyIcon.iconset/icon_16x16.png
sips -z 32 32     icon.png --out MyIcon.iconset/icon_16x16@2x.png
sips -z 32 32     icon.png --out MyIcon.iconset/icon_32x32.png
sips -z 64 64     icon.png --out MyIcon.iconset/icon_32x32@2x.png
sips -z 128 128   icon.png --out MyIcon.iconset/icon_128x128.png
sips -z 256 256   icon.png --out MyIcon.iconset/icon_128x128@2x.png
sips -z 256 256   icon.png --out MyIcon.iconset/icon_256x256.png
sips -z 512 512   icon.png --out MyIcon.iconset/icon_256x256@2x.png
sips -z 512 512   icon.png --out MyIcon.iconset/icon_512x512.png

iconutil -c icns MyIcon.iconset
# check manually newly created MyIcon.icns file
open MyIcon.icns
# if satisfied
mv MyIcon.icns icon.icns
rm -r MyIcon.iconset
```
