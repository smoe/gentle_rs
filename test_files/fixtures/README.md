# `test_files/fixtures`

Structured, committed fixtures used by parser/runtime tests and catalog smoke
checks.

## Layout

- `genomes/`
  - `AB011549.2.fa`
  - `AB011549.2.gb`
- `import_parity/`
  - `toy.small.fa`
  - `toy.small.gb`
  - `toy.small.embl`
  - `toy.multi.embl`
  - `toy.small.gbseq.xml`
  - `toy.small.insdseq.xml`
  - `toy.multi.gbseq.xml`
- `resources/`
  - `jaspar.edge.pfm`
  - `rebase.edge.withrefm`
  - `ucsc.rmsk.hg38.edge.txt`
- `primer3/`
  - `pairs.location_5_60.kv`
- `mapping/`
  - `ensembl_chimp_tp73_all.fasta`
  - `ensembl_human_tp73_all.fasta`
  - `ensembl_human_tp53_all.fasta`
  - `ensembl_human_tp63_all.fasta`
  - `ensembl_mouse_trp73_all.fasta`
- `sequencing_confirmation/`
  - `README.md`
- `microarray_tracks/`
  - `README.md`
  - `clariomd.synthetic.manifest.json`
  - `clariomd.synthetic.hg19_projected.manifest.json`
  - `clariomd.synthetic.hg19-to-hg38.tsv`
  - `clariomd.synthetic.AdTAp73alpha-AdGFP.tsv`
  - `clariomd.synthetic.AdTAp73beta-AdGFP.tsv`
  - `clariomd.tp73_vendor_subset.manifest.json`
  - `clariomd.tp73_vendor_subset.AdTAp73alpha-AdGFP.tsv`
  - `clariomd.tp73_vendor_subset.AdTAp73beta-AdGFP.tsv`
- `probe_region_outputs/`
  - `README.md`
  - `clariom_pm_probe_interpretation/`
- `affymetrix_clariom_d_human_na36_hg38_subset/`
  - `README.md`
  - `clariom_d_human_na36_hg38_gene_panel.probesets.tsv`
  - `clariom_d_human_na36_hg38_gene_panel.transcripts.tsv`

## Provenance and usage

### `genomes/AB011549.2.fa` + `genomes/AB011549.2.gb`

- Origin: E. coli reference sequence/annotation pair used as a local catalog
  fixture.
- Primary usage:
  - `assets/genomes.json` `LocalProject` entry.
  - Genome catalog validation/smoke checks in CI.
- Purpose: deterministic local test genome that avoids remote network
  dependencies.

### `import_parity/toy.small.*`

- Origin: hand-crafted synthetic 120 bp sequence represented in FASTA, GenBank,
  EMBL, NCBI GenBank XML (`GBSet/GBSeq`), and INSD XML
  (`INSDSet/INSDSeq`) forms.
- Primary usage:
  - Cross-format import/parity tests.
  - EMBL parser location/qualifier regression coverage.
  - XML import normalization tests against current GenBank baseline.
  - INSD XML adapter tests that verify it shares the same normalized feature
    semantics as GenBank XML.
- Purpose: tiny, reviewable fixture set for format-difference debugging.

### `import_parity/toy.multi.embl`

- Origin: hand-crafted two-record EMBL fixture.
- Primary usage:
  - Multi-record EMBL ingestion tests.
  - Wrapped feature-location and wrapped qualifier regression coverage.
- Purpose: ensure parser behavior stays stable for realistic EMBL line wrapping
  and record boundaries.

### `import_parity/toy.multi.gbseq.xml`

- Origin: hand-crafted two-record NCBI GenBank XML fixture paired with
  `toy.multi.embl`.
- Deterministic self-creation:
  - encode the two 24 bp synthetic records from `toy.multi.embl` as
    `GBSet/GBSeq` records
  - keep record A's `gene` location interval-only
    (`GBFeature_intervals` without `GBFeature_location`)
  - keep record A's `pseudo` qualifier value-less
  - keep record B's `misc_feature` location text-only
    (`GBFeature_location` without `GBFeature_intervals`)
- Primary usage:
  - XML multi-record import tests.
  - XML interval-derived location regression tests.
  - XML flag/value-less qualifier preservation tests.
- Purpose: keep XML edge-case coverage small and paired with the existing EMBL
  multi-record fixture.

### `resources/jaspar.edge.pfm`

- Origin: hand-crafted edge fixture from commit
  `2154eb966409b91f714c2139113c907028c42634`.
- Primary usage:
  - `src/resource_sync.rs` parser tests.
  - `src/engine_shell.rs` shared-shell execution tests for
    `resources sync-jaspar` + motif-registry reload side effects.
  - Runtime paths `resources sync-jaspar` / GUI import dialog.

### `resources/rebase.edge.withrefm`

- Origin: hand-crafted edge fixture from commit
  `2154eb966409b91f714c2139113c907028c42634`.
- Primary usage:
  - `src/resource_sync.rs` parser tests.
  - `src/engine_shell.rs` shared-shell execution tests for
    `resources sync-rebase`.
  - Runtime paths `resources sync-rebase` / GUI import dialog.

### `resources/ucsc.rmsk.hg38.edge.txt`

- Origin: small hand-curated excerpt from the UCSC hg38 `rmsk` Table Browser
  schema/sample rows for "RepeatMasker - Repeating Elements by
  RepeatMasker".
- Deterministic recreation:
  - open the UCSC schema URL:
    `https://genome.ucsc.edu/cgi-bin/hgTables?db=hg38&hgta_doSchema=describe+table+schema&hgta_group=rep&hgta_table=rmsk&hgta_track=rmsk`
  - copy the header and first four sample rows shown for hg38 `rmsk`
  - keep tabular whitespace normalized to single spaces
- Primary usage:
  - `src/ucsc_rmsk.rs` parser/resource-snapshot tests.
  - `src/engine_shell.rs` shared-shell execution tests for
    `resources sync-ucsc-rmsk`.
- Runtime/parser role:
  - exercises the 17-column UCSC RepeatMasker `.out` table contract without
    committing the full hg38 `rmsk.txt.gz` download.

### `primer3/pairs.location_5_60.kv`

- Origin: hand-crafted synthetic Primer3 key/value response fixture.
- Deterministic self-creation:
  - create a plain-text file with:
    - `PRIMER_PAIR_NUM_RETURNED=1`
    - `PRIMER_LEFT_0=5,20`
    - `PRIMER_RIGHT_0=79,20`
    - `PRIMER_PAIR_0_PRODUCT_SIZE=75`
    - terminal `=`
- Primary usage:
  - `src/engine_shell.rs` adapter-equivalence tests for
    internal-vs-Primer3 primer report normalization using a fake local
    `primer3_core` script.
- Purpose: deterministic offline coverage for Primer3 normalization/provenance
  behavior without requiring a system Primer3 installation in CI.

### `affymetrix_clariom_d_human_na36_hg38_subset/*`

- Origin: derived minimal subset from Thermo Fisher/Affymetrix Clariom D Human
  na36 hg38 NetAffx CSV support ZIPs.
- Deterministic recreation:
  - manually obtain the login-walled vendor ZIPs as described in
    `data/resources/affymetrix/clariom_d_human_na36_hg38/README.md`
  - run `scripts/extract_clariomd_gene_panel_fixture.py` with those ZIP paths
    and this fixture directory as `--output-dir`
- Primary usage:
  - offline parser/projection tests that need realistic Clariom D probeset and
    transcript-cluster IDs for a small TP73/PATZ1/p53-pathway gene panel
  - future array-evidence development without committing full vendor NetAffx
    annotation payloads
- Runtime/parser role:
  - intentionally stripped to minimal IDs, hg38 coordinates, feature types, and
    gene symbols; it is not a replacement for the full vendor annotations.

### `probe_region_outputs/clariom_pm_probe_interpretation/*`

- Origin: hand-crafted synthetic completed helper-output directory for the
  ClawBio Clariom PM-probe interpretation workflow example.
- Deterministic recreation:
  - use the committed `test_files/tp73.ncbi.gb` GRCh38.p14 chromosome-1 locus
  - choose two tiny PM-probe intervals within the first TP73 exon and one
    parent probeset-region row
  - write fixed sample, condition, and log2 fold-change values
  - mark the PM probe rows as `probe_level_input`
- Primary usage:
  - ClawBio example
    `request_workflow_clariom_pm_probe_interpretation.json`
  - offline smoke coverage for projection plus review-only probe-region
    evidence interpretation without running R/APT or committing CEL/vendor
    binary files.
- Runtime/parser role:
  - tiny stand-in for a completed
    `arrays import-apt-probe-region-output ... --probe-intensity ...` directory;
    it is not biological evidence and does not assess probe specificity,
    multi-hit status, or isoform support.

### `mapping/*.fasta` (TP73/TP53 benchmark set)

- Origin: copied from the legacy TP73 mapping corpus in
  `test_files/mapping/True_TP73/` to provide a small, deterministic benchmark pack
  under committed fixtures.
- Current provenance status:
  - biological source is Ensembl transcript FASTA exports (human/chimp TP73,
    human TP53, human TP63, mouse Trp73), as indicated by headers/filenames.
  - contributing Ensembl release for this cDNA fixture pack is `115`.
  - exact original export URLs/queries are not yet frozen in-repo; this is a
    known provenance follow-up item in `docs/roadmap.md`.
- Primary usage:
  - `src/engine.rs` seed-filter regression:
    `test_tp73_seed_filter_cross_species_and_tp53_specificity_sets`.
- Purpose:
  - keep CI/stable tests on compact mapping fixtures,
  - include TP53 as a close-family negative benchmark (as requested for better
    practical signal than TP73-only stress sets).

### `sequencing_confirmation/README.md`

- Origin: in-repo provenance shortlist for future sequencing-confirmation
  fixtures and public benchmarks.
- Primary usage:
  - plan deterministic phase-1 confirmation fixtures before raw trace intake
  - pin the first real public read benchmark and ABI parser fixture sources
  - keep future sequencing fixtures aligned with the architecture provenance
    rule before payloads are committed
- Current contents:
  - `U13852` / `pGEX-3X` synthetic-fixture plan
  - `PRJNA1066256` / `SRR27605537` real-read benchmark candidate
  - Biopython ABI fixture shortlist pinned to a repository commit

## Large exploratory XML samples

- Large one-off XML records are intentionally excluded from default committed
  fixtures to keep repository size controlled.
- If needed for local parser experiments, fetch/regenerate into a local
  scratch path and do not commit by default.
