# Repeats and Regulatory Annotations

<!-- manual-counterparts:start -->
Related GUI reference: [Locus Figure And Evidence Tracks](../gui/locus-figure-and-evidence-tracks.md).
These topic links do not imply identical feature coverage across interfaces.
<!-- manual-counterparts:end -->

- Shared shell:
  - Commands:
    - Repeat-cohort helper notes (`features repeat-query`,
      `features repeat-overlaps`, `features materialize-repeats`,
      `features repeat-cohort`, `features window-cohort-tfbs`):
      - non-mutating structured result schemas:
        `gentle.repeat_annotation_query.v1`,
        `gentle.sequence_repeat_overlap.v1`,
        `gentle.repeat_environment_cohort.v1`, and
        `gentle.window_cohort_tfbs.v1`
      - mutating materialization schema:
        `gentle.repeat_feature_materialization.v1`
      - `--rmsk PATH` accepts UCSC `rmsk.txt` or `rmsk.txt.gz` rows and keeps
        raw `repName`, `repClass`, `repFamily`, strand, genomic span, and
        optional score/divergence fields
      - filters can use raw labels or normalized aliases such as `LINE/L1`,
        `SINE/Alu`, and `LTR/ERV`
      - `repeat-overlaps` uses a prepared interval index plus the selected
        sequence's genome anchor provenance to return clipped local repeat
        spans with strand flipped correctly on reverse anchors; the prepared
        index resolves common chromosome aliases such as `1` versus `chr1`
      - repeat-overlap reports also include neutral summary metrics:
        query length, total repeat-overlap bp, union-covered query bp, coverage
        fraction, nearest-repeat distance, and class/family/alias summaries
      - `materialize-repeats` writes those overlaps back as regular
        `repeat_region` features with UCSC `rmsk_*` and
        `repName`/`repClass`/`repFamily` qualifiers plus score/divergence
        fields; omit `--append` to replace prior generated UCSC-rmsk features
        for the sequence
      - `genomes extract-region|extract-gene|extract-promoter` and matching
        `helpers ...` commands can immediately materialize repeats on the new
        anchored sequence with `--rmsk-index PATH`; `--max-repeat-features N`
        caps the generated annotations and `--append-repeat-features` keeps
        existing generated rmsk features
      - materialized repeat features can be exported as BED with
        `features export-bed SEQ_ID out.bed --kind repeat_region --coordinate-mode genomic --include-qualifiers`
      - `repeat-cohort` builds one repeat-locus row per selected annotation and
        stores all available geometry frames while marking missing transcript,
        UTR, or CDS-stop context explicitly
      - geometry modes are `repeat_midpoint`, `transcript_5utr_start`,
        `pol2_promoter_upstream`, and `cds_stop_context`; promoter-like motif
        comparisons should usually start with `pol2_promoter_upstream`
      - `window-cohort-tfbs` consumes a saved cohort JSON and summarizes motif
        score tracks over the selected windows; this is the ClawBio-friendly
        composition boundary rather than a pre-bundled GENtle dossier
      - example with a local UCSC file:
        - `gentle_cli shell 'features repeat-query "Human GRCh38 Ensembl 116" --rmsk data/rmsk.txt.gz --alias LINE/L1 --limit 100 --path /tmp/grch38_l1_repeats.json'`
        - `gentle_cli shell 'features repeat-overlaps grch38_tp53 --index data/resources/ucsc.rmsk.hg38.interval-index.json --range 0..5000 --limit 100 --path /tmp/tp53_repeat_overlaps.json'`
        - `gentle_cli shell 'features materialize-repeats grch38_tp53 --index data/resources/ucsc.rmsk.hg38.interval-index.json --max-features 1000 --path /tmp/tp53_repeat_materialization.json'`
        - `gentle_cli shell 'features repeat-cohort "Human GRCh38 Ensembl 116" --rmsk data/rmsk.txt.gz --alias LINE/L1 --geometry pol2_promoter_upstream --upstream-bp 2000 --downstream-bp 2000 --limit 50 --path /tmp/grch38_l1_promoter_context.json'`
        - `gentle_cli shell 'features window-cohort-tfbs /tmp/grch38_l1_promoter_context.json --motif SP1 --motif TP73 --score-kind llr_background_tail_log10 --path /tmp/grch38_l1_promoter_tfbs.json'`
    - Optional ENCODE SCREEN Registry V4 cCRE notes:
      - inspect the built-in species/assembly catalog without downloading:
        - `gentle_cli resources list-encode-ccre-sources`
        - `gentle_cli resources list-encode-ccre-sources --assembly mm10`
      - explicitly install the current enhancer-like-signature BED and sparse
        content-bound index:
        - human: `gentle_cli resources install-encode-ccres screen_registry_v4_grch38_els`
        - mouse: `gentle_cli resources install-encode-ccres screen_registry_v4_mm10_els`
      - an already downloaded BED can be indexed offline with
        `gentle_cli resources prepare-encode-ccre-index SOURCE_ID BED_PATH --output INDEX.json`
      - query a genome-anchored locus without changing project state:
        - `gentle_cli features encode-ccre-overlaps LOCUS_SEQ_ID --index INDEX.json --range 0..5000 --class pELS --class dELS --path /tmp/locus_ccres.json`
      - materialize verified rows for the existing regulatory-feature display:
        - `gentle_cli features materialize-encode-ccres LOCUS_SEQ_ID --index INDEX.json --max-features 1000 --path /tmp/locus_ccre_materialization.json`
      - the source id, not the filename, selects human GRCh38/hg38 versus mouse
        GRCm38/mm10. The engine cross-checks that choice against sequence anchor
        provenance and refuses incompatible species/assemblies.
      - the V4 `ELS` files contain `pELS` and `dELS`, not the complete SCREEN
        registry. Overlap is candidate regulatory evidence, not proof of
        enhancer activity or target-gene regulation.
      - the large resource files are fully optional: GENtle performs
        no SCREEN I/O or network access at startup, and missing files affect
        only an explicit cCRE command.
    - Optional Ensembl Regulation annotation notes:
      - inspect the pinned source catalog without downloading:
        - `gentle_cli resources list-ensembl-regulation-sources`
        - `gentle_cli resources list-ensembl-regulation-sources --assembly mm39`
      - explicitly fetch and index the release-bound annotation:
        - human: `gentle_cli resources install-ensembl-regulatory-features ensembl_regulation_2026_08_grch38`
        - mouse: `gentle_cli resources install-ensembl-regulatory-features ensembl_regulation_2026_08_grcm39`
      - a complete saved API response can be installed offline with `--input`;
        a partial/paginated response is rejected rather than mistaken for a
        complete annotation snapshot
      - query without changing project state:
        - `gentle_cli features ensembl-regulation-overlaps LOCUS_SEQ_ID --index INDEX.json --range 0..5000 --feature-type promoter --feature-type enhancer --path /tmp/locus_ensembl_regulation.json`
      - materialize verified rows into the existing feature display:
        - `gentle_cli features materialize-ensembl-regulation LOCUS_SEQ_ID --index INDEX.json --max-features 1000 --path /tmp/locus_ensembl_regulation_materialization.json`
      - Ensembl mouse annotations use GRCm39/mm39 and are not interchangeable
        with SCREEN Registry V4 mouse GRCm38/mm10 resources
      - source descriptors expose the separate regulatory-activity matrix and
        primary-analysis page. This annotation route does not download those
        files or claim biosample activity from overlap alone.
      - the canonical `gene-locus prepare @request.json` workflow can consume
        an installed index directly, without materializing features. Add an
        optional request member such as:

        ```json
        "ensembl_regulation": {
          "source_id": "ensembl_regulation_2026_08_grch38",
          "interval_index_path": "data/resources/ensembl-regulation.2026-08.grch38.interval-index.json",
          "feature_types": ["promoter", "enhancer", "open_chromatin_region", "ctcf", "emar"],
          "max_rows": 1000,
          "availability_policy": "retain_unavailable"
        }
        ```

        Use `required` instead when a missing or unverifiable snapshot must
        stop figure generation. Provenance-sensitive reruns may additionally
        set `expected_index_sha256` and `expected_intervals_sha256` to the exact
        `sha256:...` values from a reviewed prior report/receipt. The report and
        receipt retain source/release and observed content hashes while local
        paths follow the request's sharing policy.
      - the resulting combined SVG contains safe live feature/data-access links
        to the exact `regulation.ensembl.org` host. PDF export writes matching
        PDF link annotations rather than relying on rasterized blue text. The
        link may show the provider's current live record; plotted coordinates
        remain those of the recorded local release and digest.
      - no Ensembl Regulation file or network endpoint is touched at startup;
        absence affects only explicit Ensembl Regulation commands.

<!-- manual-page-nav:start -->

---

[CLI Manual](../cli.md) | [Previous: TFBS and JASPAR Scans](tfbs.md) | [Next: Tracks and Arrays](tracks-arrays.md)
<!-- manual-page-nav:end -->
