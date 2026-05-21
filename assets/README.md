# GENtle Assets

This folder contains a series of accessory files that are meant to be available
at the application's start or run time. These have mostly been carried over from the C++ version of GENtle.

amino_acids.json:: Describes properties of amino acids (one and three letter abbreviation, molecular weight, [isoelectric point](https://en.wikipedia.org/wiki/Isoelectric_point) and various scores like [Kyte-Doolittle](https://en.wikipedia.org/wiki/Hydrophobicity_scales), [Chou-Fasman](https://en.wikipedia.org/wiki/Chou%E2%80%93Fasman_method), [Hopp-Woods](https://en.wikipedia.org/wiki/Hopp%E2%80%93Woods_scale) ) 
codon_catalog.csv:: Tabular representation of nucleotide representations that should be selected by default for a given amino acid. A later development of that table should represent the codon frequencies for each species, as in https://www.creative-biostructure.com/codon-usage-frequency-table.htm .
codon_tables.json:: A codon translation table.
dna_ladders.json:: Collection of DNA ladders, i.e. fragments of DNA of well-defined length and relative abundance that are accompanying gels to identify the molecular weight of DNA sequences.
enzymes.json:: List of DNA restriction enzymes (endonucleases) and proteases
genomes.json::
helper_genomes.json::
cutrun.json:: Starter CUT&RUN dataset catalog for processed evidence (BED/BigWig) and future raw-read reuse.
cutrun.d/:: Additional built-in CUT&RUN catalog shards. The Rostock p73 `E-MTAB-15709` shard records SRA-backed paired-end runs and keeps full raw-read acquisition explicit.
host_profiles.json:: Starter host/strain catalog for construct-reasoning inspection in GUI/agent-facing workflows.
i18n/:: Embedded GUI translation catalogs. They localize visible interface
chrome only; shared shell commands, saved records, protocol fields, and
scientific identifiers remain deterministic English.
blast_defaults.json:: Default BLAST option layer (`task`, `max_hits`, optional thresholds) used when no project/request override is provided.
reporter_catalog.json:: Small V1 reporter-selection catalog used by the
offline reporter recommender. Records carry source links, license status,
benign-use scope, sequences/checksums, practical assay tags, and spectral or
color metadata when curated.
panels/tp53_isoforms_v1.json:: Curated TP53 isoform architecture panel used by isoform expert and protein-gel examples.
panels/tp73_isoforms_v1.json:: Local TP73 isoform curation seed that records lab/public hybrid transcript-class knowledge for assay design without treating public disease-transcript coverage as complete; it now also carries a local IEGT Rostock TAp73alpha coding-cDNA evidence record whose exact source cell line is still pending clarification.
panels/tp73_long_range_cdna_virtual_panel_v1.json:: Local TP73 long-range cDNA selector panel that explicitly materializes virtual 5' x 3' isoform combinations absent from the bundled public RefSeq-style annotation; the companion FASTA stores the theoretical cDNA sequences.
panels/tp73_delta_ex2_3_isoforms_v1.json:: Local TP73 DEx2,3 alpha/beta panel with RefSeq-derived exon-chain models and an IEGT plasmid-contained D_Ex2/3_p73beta coding-cDNA evidence record.
jaspar.motifs.json:: The JASPAR database transformed into a JSON format that is meant to be mostly compatible with the JSON format offered by the JASPAR project: ```gzip -dc data/JASPAR_2022.txt.gz | perl scripts/pfm2json.pl | jq --compact-output > assets/jaspar_2022.json```
ncoils.matrix::

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
- `helper_genomes.json` is a curated starter catalog for helper references
  (host genomes plus common online vector backbones), intended to be copied and
  edited for local vector inventories.
  Entries can use explicit local/remote URLs or `genbank_accession` for
  NCBI EFetch-derived FASTA + GenBank annotation during prepare/index.
  Metadata-only vector candidates may also record sequence availability,
  redistribution status, biological-safety notes, source/procurement metadata,
  marker/origin components, host compatibility, and whether the record is
  usable as an empty backbone.
  The shipped defaults avoid absolute lab-specific paths.
- `cutrun.json` is the default CUT&RUN catalog for V1 processed datasets.
  Entries can point at prepared/local or remote peaks/signal assets and are
  discoverable through the same built-in/system/project overlay chain as other
  shared catalogs.
- `cutrun.d/` contains additional built-in CUT&RUN catalog shards. The
  `rostock_p73_sra.json` shard maps `E-MTAB-15709` / `PRJEB100610` SRA runs to
  GENtle dataset IDs, while `cutrun/rostock_p73_sra/` documents ENA FASTQ URLs,
  byte counts, and checksums used to keep the raw sequences reproducible without
  committing them.
- `host_profiles.json` is the starter host-profile catalog used by construct
  reasoning and the Planning-window host browser. It is intentionally
  human-editable and source-noted rather than hidden inside compiled logic.
- `assets/panels/` contains curated isoform-panel resources. These are
  machine-readable panels rather than hidden code paths; local/lab-curated
  panels may carry structured curation metadata so public accession anchors and
  local biological insight remain distinguishable.
- `reporter_catalog.json` is a deliberately small, auditable reporter catalog
  seed rather than a public-registry mirror. Normal recommendation runs use this
  local pack offline; future source refresh or expanded annotation should be an
  explicit import/update step so provenance and license status stay visible.

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
