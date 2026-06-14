# Clariom D Human na36 hg38 support files

This directory is the canonical local staging location for Thermo Fisher
Clariom D Human hg38 support ZIPs used during probe/probeset-region
development. The files are login-walled vendor downloads, so GENtle documents
and checks for them here but does not auto-download them and does not commit
the ZIP payloads.

Expected local files:

| File | Purpose |
| --- | --- |
| `Clariom_D_Human-na36-hg38-probeset-csv.zip` | Thermo Fisher probeset coordinate/annotation CSV ZIP for Clariom D Human na36 hg38; contains `Clariom_D_Human.na36.hg38.probeset.csv` plus `XTAArray_NetAffx-CSV-Files.README.txt`. |
| `Clariom_D_Human.r1.na36.hg38.a1.transcript.csv.zip` | Thermo Fisher transcript-cluster annotation CSV ZIP for Clariom D Human r1 na36 hg38; contains `Clariom_D_Human.r1.na36.hg38.a1.transcript.csv` plus `XTAArray_NetAffx-CSV-Files.README.txt`. |

Manual download pages/files:

- `https://sec-assets.thermofisher.com/TFS-Assets/LSG/Support-Files/Clariom_D_Human-na36-hg38-probeset-csv.zip`
- `https://sec-assets.thermofisher.com/TFS-Assets/LSG/Support-Files/Clariom_D_Human.r1.na36.hg38.a1.transcript.csv.zip`

If a browser stores the files with the `TFS-Assets_LSG_Support-Files_`
prefix, normalize them into this directory with:

```bash
mkdir -p data/resources/affymetrix/clariom_d_human_na36_hg38
cp ~/Downloads/TFS-Assets_LSG_Support-Files_Clariom_D_Human-na36-hg38-probeset-csv.zip \
  data/resources/affymetrix/clariom_d_human_na36_hg38/Clariom_D_Human-na36-hg38-probeset-csv.zip
cp ~/Downloads/TFS-Assets_LSG_Support-Files_Clariom_D_Human.r1.na36.hg38.a1.transcript.csv.zip \
  data/resources/affymetrix/clariom_d_human_na36_hg38/Clariom_D_Human.r1.na36.hg38.a1.transcript.csv.zip
```

To check GENtle's view of the local support files without running R, APT, or
any downloader:

```bash
cargo run --bin gentle_cli -- arrays probe-regions \
  --dataset E-MTAB-14704 \
  --gene TP73 \
  --platform Clariom_D_Human \
  --dry-run
```

The resulting `gentle.probe_region_plan.v1` report lists these paths under
`plan.annotation_source.vendor_support_files[]`.
