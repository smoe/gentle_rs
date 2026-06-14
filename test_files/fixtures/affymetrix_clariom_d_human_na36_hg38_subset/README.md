# Clariom D Human na36 hg38 Gene-Panel Fixture

## Origin

This is a small derived fixture from Thermo Fisher/Affymetrix Clariom D Human
NetAffx CSV support files:

- `Clariom_D_Human-na36-hg38-probeset-csv.zip`
  - member: `Clariom_D_Human.na36.hg38.probeset.csv`
  - vendor header inspected by GENtle: `chip_type=Clariom_D_Human`,
    `genome-version=hg38`
- `Clariom_D_Human.r1.na36.hg38.a1.transcript.csv.zip`
  - member: `Clariom_D_Human.r1.na36.hg38.a1.transcript.csv`
  - vendor header inspected by GENtle: `chip_type=Clariom_D_Human`,
    `genome-version=hg38`

The original ZIPs are login-walled vendor files and are not committed. The
vendor README states that the CSV contents are covered by Affymetrix/Thermo
Fisher terms of use. This fixture is therefore deliberately non-substitutable:
it keeps only the minimal ID, coordinate, type, and gene-symbol fields required
for GENtle parser/projection tests, and omits broad descriptive annotations,
GO/pathway/domain text, mRNA descriptions, and the vendor README.

## Gene Panel

Rows are retained when `gene_assignment` contains one of:

`E2F1`, `TP73`, `SP1`, `PATZ1`, `TP53`, `TP63`, `IL6`, `IL10`, `FUS`, `TERT`,
`TARDBP`, `MDM2`, `CDKN1A`, `BAX`, `GADD45A`, `MYC`, `RB1`, `ESR1`, `GAPDH`,
`ACTB`, `SRSF1`.

## Deterministic Recreation

Place the vendor ZIPs locally as described in
`data/resources/affymetrix/clariom_d_human_na36_hg38/README.md`, or pass their
explicit paths:

```bash
python3 scripts/extract_clariomd_gene_panel_fixture.py \
  --probeset-zip /path/to/Clariom_D_Human-na36-hg38-probeset-csv.zip \
  --transcript-zip /path/to/Clariom_D_Human.r1.na36.hg38.a1.transcript.csv.zip \
  --output-dir test_files/fixtures/affymetrix_clariom_d_human_na36_hg38_subset
```

Expected output sizes:

- `clariom_d_human_na36_hg38_gene_panel.probesets.tsv`: 906 data rows
- `clariom_d_human_na36_hg38_gene_panel.transcripts.tsv`: 23 data rows

## Primary Usage

- Offline tests for Clariom D hg38 gene/probeset coordinate handling.
- Future probe/probeset projection tests that need realistic IDs and
  gene-panel breadth without committing the full vendor annotation set.
- Manual development around TP73-family, PATZ1, E2F, p53-pathway, cytokine,
  RNA-binding, telomerase, and housekeeping probesets.

