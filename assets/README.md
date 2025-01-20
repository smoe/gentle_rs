# GENtle Assets

This folder contains a series of accessory files that are meant to be available
at the application's start or run time. These have mostly been carried over from the C++ version of GENtle.

- amino_acids.json: Describes properties of amino acids (one and three letter abbreviation, molecular weight, [isoelectric point](https://en.wikipedia.org/wiki/Isoelectric_point) and various scores like [Kyte-Doolittle](https://en.wikipedia.org/wiki/Hydrophobicity_scales), [Chou-Fasman](https://en.wikipedia.org/wiki/Chou%E2%80%93Fasman_method), [Hopp-Woods](https://en.wikipedia.org/wiki/Hopp%E2%80%93Woods_scale) ) 
- blank.mysql:
- blank.sqlite3:
- codon_catalog.csv: Tabular representation of nucleotide representations that should be selected by default for a given amino acid. A later development of that table should represent the codon frequencies for each species, as in https://www.creative-biostructure.com/codon-usage-frequency-table.htm .
- codon_tables.json: A codon translation table.
- commonvectors.db: Sqlite3 representation of the most common vectors.
- dna_markers.json: Collection of DNA ladders, i.e. fragments of DNA of well-defined length and relative abundance that are accompanying gels to identify the molecular weight of DNA sequences.
- enzymes.json: List of DNA restriction enzymes (endonucleases) and proteases
- ncoils.matrix: File needed for the execution of the ncoils program, mostly obsolete since alpha Fold.
- translations.csv: Translations of strings used in the user interface.

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
