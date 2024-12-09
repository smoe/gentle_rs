# GENtle Assets

This folder contains a series of accessory files that are meant to be available
at the application's start or run time.

amino_acids.json::
blank.mysql::
blank.sqlite3::
codon_catalog.csv::
codon_tables.json::
commonvectors.db::
dna_markers.json::
enzymes.json::
ncoils.matrix::
translations.csv::

# "Finder" icon (seen in OS' task bar)
icon.png:: PNG with icon as used at startup by application
icon.icns:: Icon-set for MacOS

Script to transforms the .png icon to an icon set for the MacOS
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
cp icon.png MyIcon.iconset/icon_512x512@2x.png

iconutil -c icns MyIcon.iconset
# check manually
open MyIcon.iconset
# if satisfied
mv MyIcon.icns icon.icns
rm -r MyIcon.iconset
```
