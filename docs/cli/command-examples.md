# Independent Command Examples

<!-- manual-counterparts:start -->
No direct GUI counterpart: This chapter documents an automation interface, example or command-only workflow; no dedicated GUI chapter is paired with it.
<!-- manual-counterparts:end -->

```bash
cargo run --bin gentle_cli -- --version
cargo run --bin gentle_cli -- capabilities
cargo run --bin gentle_cli -- state-summary
cargo run --bin gentle_cli -- containers set-exclusive container-3 false
cargo run --bin gentle_cli -- help
cargo run --bin gentle_cli -- help candidates generate
cargo run --bin gentle_cli -- help --format json
cargo run --bin gentle_cli -- op '<operation-json>'
cargo run --bin gentle_cli -- op op.json
cargo run --bin gentle_cli -- workflow '<workflow-json>'
cargo run --bin gentle_cli -- workflow docs/examples/workflows/load_branch_reverse_complement_pgex_fasta.json
cargo run --bin gentle_cli -- workflow @docs/examples/workflows/load_branch_reverse_complement_pgex_fasta.json
cargo run --bin gentle_cli -- workflow @docs/examples/workflows/rna_reads_interpret_cdna_tp73_template.json
cargo run --bin gentle_cli -- batch plan samples.tsv --template batch_templates/tp73_rna.workflow.json --out-dir runs/tp73 --state-template tp73_ensembl.seed.state.json --state-mode copy-per-row --emit local
cargo run --bin gentle_cli -- batch run samples.tsv --template batch_templates/tp73_rna.workflow.json --limit 1
cargo run --bin gentle_cli -- --progress op '<operation-json>'
cargo run --bin gentle_cli -- --progress-stdout workflow '<workflow-json>'
cargo run --bin gentle_cli -- export-state state.json
cargo run --bin gentle_cli -- import-state state.json
cargo run --bin gentle_cli -- save-project project.gentle.json
cargo run --bin gentle_cli -- load-project project.gentle.json
cargo run --bin gentle_cli -- render-svg pgex linear pgex.linear.svg
cargo run --bin gentle_cli -- render-svg pgex circular pgex.circular.svg
cargo run --bin gentle_cli -- render-dotplot-svg tp73_cdna dotplot_primary tp73.dotplot.svg
cargo run --bin gentle_cli -- render-dotplot-svg tp73_cdna dotplot_primary tp73.dotplot.svg --flex-track flex_primary --display-threshold 0.2 --intensity-gain 1.8
cargo run --bin gentle_cli -- render-dotplot-svg tp73_genomic_local tp73_multi_isoform_overlay_dotplot tp73.overlay.right.svg --overlay-x-axis right_aligned_bp --display-threshold 0.08 --intensity-gain 2.8
cargo run --bin gentle_cli -- render-dotplot-svg tp73_genomic_local tp73_multi_isoform_overlay_dotplot tp73.overlay.anchor.svg --overlay-x-axis shared_exon_anchor --overlay-anchor-exon 55034..55276 --display-threshold 0.08 --intensity-gain 2.8
cargo run --bin gentle_cli -- render-rna-svg rna_seq rna.secondary.svg
cargo run --bin gentle_cli -- rna-info rna_seq
cargo run --bin gentle_cli -- render-lineage-svg lineage.svg
cargo run --bin gentle_cli -- protocol-cartoon list
cargo run --bin gentle_cli -- protocol-cartoon render-svg gibson.two_fragment gibson.protocol.svg
cargo run --bin gentle_cli -- protocol-cartoon render-svg pcr.assay.pair pcr_pair.protocol.svg
cargo run --bin gentle_cli -- protocol-cartoon render-svg pcr.assay.pair.with_tail pcr_pair_with_tail.protocol.svg
cargo run --bin gentle_cli -- protocol-cartoon render-svg pcr.oe.substitution pcr_oe_substitution.protocol.svg
cargo run --bin gentle_cli -- protocol-cartoon render-svg pcr.assay.qpcr qpcr.protocol.svg
cargo run --bin gentle_cli -- protocol-cartoon render-template-svg docs/examples/protocol_cartoon/demo_template.json demo.protocol.svg
cargo run --bin gentle_cli -- protocol-cartoon template-validate docs/examples/protocol_cartoon/demo_template.json
cargo run --bin gentle_cli -- protocol-cartoon render-with-bindings docs/examples/protocol_cartoon/demo_template.json docs/examples/protocol_cartoon/demo_bindings.json demo.bound.protocol.svg
cargo run --bin gentle_cli -- protocol-cartoon template-export gibson.two_fragment gibson.template.json
cargo run --bin gentle_cli -- gibson preview @docs/examples/plans/gibson_destination_first_single_insert.json --output gibson.preview.json
cargo run --bin gentle_cli -- gibson apply @docs/examples/plans/gibson_destination_first_single_insert.json
cargo run --bin gentle_cli -- shell 'help'
cargo run --bin gentle_cli -- shell 'state-summary'
cargo run --bin gentle_cli -- shell 'op @op.json'
cargo run --bin gentle_cli -- render-pool-gel-svg frag_1,frag_2 digest.gel.svg
cargo run --bin gentle_cli -- render-pool-gel-svg frag_1,frag_2 digest.gel.svg --ladders "NEB 100bp DNA Ladder,NEB 1kb DNA Ladder"
cargo run --bin gentle_cli -- render-pool-gel-svg - digest.gel.svg --containers container-3,container-8
cargo run --bin gentle_cli -- render-pool-gel-svg - digest.gel.svg --arrangement arrangement-2
cargo run --bin gentle_cli -- render-gel-svg - digest.gel.svg --arrangement arrangement-2
cargo run --bin gentle_cli -- arrange-serial container-3,container-8 --id arrangement-2 --name "Digest run A" --ladders "NEB 100bp DNA Ladder"
cargo run --bin gentle_cli -- racks create-from-arrangement arrangement-2 --rack-id rack-1 --name "Bench day A"
cargo run --bin gentle_cli -- racks place-arrangement arrangement-3 --rack rack-1
cargo run --bin gentle_cli -- racks move rack-1 --from A2 --to A4
cargo run --bin gentle_cli -- racks move rack-1 --from A1 --to B1 --block
cargo run --bin gentle_cli -- racks move-samples rack-1 --from A1 --from A3 --to A2
cargo run --bin gentle_cli -- racks move-blocks rack-1 --arrangement arrangement-2 --arrangement arrangement-3 --to B2
cargo run --bin gentle_cli -- racks show rack-1
cargo run --bin gentle_cli -- racks labels-svg rack-1 arrangement-2.labels.svg --arrangement arrangement-2
cargo run --bin gentle_cli -- racks fabrication-svg rack-1 rack.fabrication.svg --template storage_pcr_tube_rack
cargo run --bin gentle_cli -- racks openscad rack-1 rack.scad --template pipetting_pcr_tube_rack
cargo run --bin gentle_cli -- racks set-profile rack-1 plate_96
cargo run --bin gentle_cli -- ladders list
cargo run --bin gentle_cli -- ladders list --filter NEB
cargo run --bin gentle_cli -- ladders list --molecule rna
cargo run --bin gentle_cli -- ladders export dna_ladders.snapshot.json
cargo run --bin gentle_cli -- ladders export dna_ladders.neb.json --filter NEB
cargo run --bin gentle_cli -- ladders export rna_ladders.snapshot.json --molecule rna
cargo run --bin gentle_cli -- export-pool frag_1,frag_2 digest.pool.gentle.json "BamHI+EcoRI digest pool"
cargo run --bin gentle_cli -- import-pool digest.pool.gentle.json imported
cargo run --bin gentle_cli -- services status
cargo run --bin gentle_cli -- services guide --channel telegram
cargo run --bin gentle_cli -- services handoff --scope clawbio --output artifacts/service_handoff.json
cargo run --bin gentle_cli -- resources status
cargo run --bin gentle_cli -- resources sync-rebase rebase.withrefm data/resources/rebase.enzymes.json --commercial-only
cargo run --bin gentle_cli -- resources sync-jaspar JASPAR2026_CORE_non-redundant_pfms_jaspar.txt data/resources/jaspar.motifs.json
cargo run --bin gentle_cli -- resources sync-ucsc-rmsk rmsk.txt.gz data/resources/ucsc.rmsk.hg38.json --assembly hg38
cargo run --bin gentle_cli -- resources install-ucsc-rmsk --assembly hg38
cargo run --bin gentle_cli -- resources prepare-ucsc-rmsk-index data/resources/ucsc.rmsk.hg38.json data/resources/ucsc.rmsk.hg38.interval-index.json
cargo run --bin gentle_cli -- resources suggest-ucsc-rmsk-index --assembly hg38 --output rmsk.indexing.json
cargo run --bin gentle_cli -- resources sync-jaspar-remote-metadata --filter TP --limit 50 data/resources/jaspar.remote_metadata.json
cargo run --bin gentle_cli -- resources summarize-jaspar --motif SP1 --motif REST --random-length 10000 --seed 123 --output jaspar.summary.json
cargo run --bin gentle_cli -- shell 'resources resolve-tf-query "Yamanaka factors" OCT4 "KLF family" --output tf_queries.json'
cargo run --bin gentle_cli -- resources benchmark-jaspar --random-length 10000 --seed 123 --output data/resources/jaspar.registry_benchmark.json
cargo run --bin gentle_cli -- resources list-jaspar --filter TP --limit 50 --output jaspar.catalog.json
cargo run --bin gentle_cli -- resources inspect-jaspar SP1 --random-length 10000 --seed 123 --fetch-remote --output jaspar.expert.json
scripts/benchmark_jaspar_interest_catalog.sh run --motif SP1 --out-dir data/resources/jaspar_interest_catalog
scripts/benchmark_jaspar_interest_catalog.sh merge --out-dir data/resources/jaspar_interest_catalog
cargo run --bin gentle_cli -- agents list
cargo run --bin gentle_cli -- agents list --catalog assets/agent_systems.json
cargo run --bin gentle_cli -- agents ask builtin_echo --prompt "summarize current project state"
cargo run --bin gentle_cli -- agents ask builtin_echo --prompt "ask: Which sequence should I use?" --execute-index 1
cargo run --bin gentle_cli -- agents ask local_llama_compat --prompt "summarize project context" --base-url http://localhost:11964 --model deepseek-r1:8b
cargo run --bin gentle_cli -- agents discover-models msty_mlx_local_compat_template
cargo run --bin gentle_cli -- agents ask msty_mlx_local_compat_template --prompt "summarize project context" --model mlx-community/granite-3.3-2b-instruct-4bit
cargo run --bin gentle_cli -- agents discover-models codex_local_stdio
cargo run --bin gentle_cli -- agents ask codex_local_stdio --prompt "summarize project context" --model gpt-5.4
cargo run --bin gentle_cli -- agents preflight pi_local_stdio
cargo run --bin gentle_cli -- agents discover-models pi_local_stdio
cargo run --bin gentle_cli -- agents ask pi_local_stdio --prompt "summarize project context" --model mistral/codestral-latest
cargo run --bin gentle_cli -- agents ask pi_local_stdio --prompt "Check the current official product information for Promega pGL4.10[luc2]" --allow-web-research
cargo run --bin gentle_cli -- op '{"PrepareGenome":{"genome_id":"ToyGenome","catalog_path":"catalog.json"}}'
cargo run --bin gentle_cli -- op '{"ExtractGenomeRegion":{"genome_id":"ToyGenome","chromosome":"chr1","start_1based":1001,"end_1based":1600,"output_id":"toy_chr1_1001_1600","annotation_scope":"core","catalog_path":"catalog.json"}}'
cargo run --bin gentle_cli -- op '{"ExtractGenomeGene":{"genome_id":"ToyGenome","gene_query":"MYGENE","occurrence":1,"output_id":"toy_mygene","catalog_path":"catalog.json"}}'
cargo run --bin gentle_cli -- op '{"DesignInsertionPrimerPairs":{"template":"toy_chr1_1001_1600","insertion":{"requested_forward_3prime_end_0based_exclusive":220,"requested_reverse_3prime_start_0based":300,"forward_extension_5prime":"GAATTC","reverse_extension_5prime":"CTCGAG","forward_window_start_0based":170,"forward_window_end_0based_exclusive":245,"reverse_window_start_0based":275,"reverse_window_end_0based_exclusive":360,"max_anchor_shift_bp":12},"forward":{"min_length":20,"max_length":30},"reverse":{"min_length":20,"max_length":30},"pair_constraints":{"require_roi_flanking":false},"min_amplicon_bp":120,"max_amplicon_bp":1200,"max_tm_delta_c":2.0,"max_pairs":50,"report_id":"insertion_demo_v1"}}'
cargo run --bin gentle_cli -- genomes list --catalog assets/genomes.json
cargo run --bin gentle_cli -- genomes list --catalog assets/helper_genomes.json
cargo run --bin gentle_cli -- genomes validate-catalog --catalog assets/genomes.json
cargo run --bin gentle_cli -- genomes update-ensembl-specs --catalog assets/genomes.json --output-catalog exports/genomes.updated.json
cargo run --bin gentle_cli -- genomes status "Human GRCh38 Ensembl 116" --catalog assets/genomes.json --cache-dir data/genomes
cargo run --bin gentle_cli -- genomes genes "Human GRCh38 Ensembl 116" --catalog assets/genomes.json --cache-dir data/genomes --filter "^TP53$" --biotype protein_coding --limit 20
cargo run --bin gentle_cli -- genomes prepare "Human GRCh38 Ensembl 116" --catalog assets/genomes.json --cache-dir data/genomes --timeout-secs 3600
cargo run --bin gentle_cli -- genomes remove-prepared "Human GRCh38 Ensembl 116" --catalog assets/genomes.json --cache-dir data/genomes
cargo run --bin gentle_cli -- genomes remove-catalog-entry "Human GRCh38 Ensembl 113" --catalog exports/genomes.custom.json
cargo run --bin gentle_cli -- genomes blast "Human GRCh38 Ensembl 116" ACGTACGTACGT --task blastn-short --max-hits 10 --options-json '{"thresholds":{"min_identity_percent":97.0,"min_query_coverage_percent":80.0}}' --catalog assets/genomes.json --cache-dir data/genomes
cargo run --bin gentle_cli -- genomes extract-region "Human GRCh38 Ensembl 116" 1 1000000 1001500 --output-id grch38_chr1_slice --annotation-scope core --catalog assets/genomes.json --cache-dir data/genomes
cargo run --bin gentle_cli -- genomes extract-gene "Human GRCh38 Ensembl 116" TP53 --occurrence 1 --output-id grch38_tp53 --catalog assets/genomes.json --cache-dir data/genomes
cargo run --bin gentle_cli -- genomes extract-gene "Human GRCh38 Ensembl 116" TP73 --occurrence 1 --output-id grch38_tp73 --catalog assets/genomes.json --cache-dir data/genomes
cargo run --bin gentle_cli -- helpers vocabulary list --filter fusion
cargo run --bin gentle_cli -- helpers vocabulary doctor --routine-catalog assets/cloning_routines.json
cargo run --bin gentle_cli -- tracks import-bed grch38_tp53 data/chipseq/peaks.bed.gz --name H3K27ac --min-score 10 --clear-existing
cargo run --bin gentle_cli -- tracks import-bigwig grch38_tp53 data/chipseq/signal.bw --name ATAC --min-score 0.2 --clear-existing
cargo run --bin gentle_cli -- tracks import-vcf grch38_tp53 data/variants/sample.vcf.gz --name Variants --min-score 20 --clear-existing
cargo run --bin gentle_cli -- arrays inspect-microarray-track data/publication_resources/rostock_p73_clariomd_e_mtab_14704/analysis/clariomd_probe_level/clariomd_microarray_track_manifest.json
cargo run --bin gentle_cli -- arrays inspect-microarray-track test_files/fixtures/microarray_tracks/clariomd.tp73_vendor_subset.manifest.json
cargo run --bin gentle_cli -- arrays project-microarray-track grch38_tp73 data/publication_resources/rostock_p73_clariomd_e_mtab_14704/analysis/clariomd_probe_level/clariomd_microarray_track_manifest.json --contrasts AdTAp73alpha-AdGFP,AdTAp73beta-AdGFP --level probeset --max-features 5000 --clear-existing
cargo run --bin gentle_cli -- arrays probe-regions --cel sample1.CEL --cel sample2.CEL --metadata samples.tsv --gene PATZ1 --gene TP73 --platform Clariom_D_Human --annotation-library path/to/NetAffx_or_APT_library --r-library-path .r-lib --condition-column condition --sample-column file --normalization rma --plot --output analysis/probe_regions --dry-run
cargo run --bin gentle_cli -- arrays probe-regions --dataset E-MTAB-14704 --gene PATZ1 --gene FUS --gene MDM2 --paired-by-replicate-suffix --platform Clariom_D_Human --plot --dry-run
cargo run --bin gentle_cli -- arrays probe-regions --dataset E-MTAB-14704 --gene TP73 --platform Clariom_D_Human --dry-run
cargo run --bin gentle_cli -- arrays inspect-probe-region-output analysis/probe_regions
cargo run --bin gentle_cli -- cutrun list --catalog assets/cutrun.json --filter TP73
cargo run --bin gentle_cli -- cutrun status tp73_release_smoke_synthetic --catalog assets/cutrun.json --cache-dir data/cutrun
cargo run --bin gentle_cli -- cutrun prepare tp73_release_smoke_synthetic --catalog assets/cutrun.json --cache-dir data/cutrun
cargo run --bin gentle_cli -- cutrun project grch38_tp73 tp73_release_smoke_synthetic --catalog assets/cutrun.json --cache-dir data/cutrun --clear-existing
cargo run --bin gentle_cli -- cutrun interpret grch38_tp53 reads_r1.fastq.gz reads_r2.fastq.gz --format fastq --layout paired_end --flank-bp 150 --report-id tp53_cutrun_reads --seed-kmer-len 9 --min-seed-matches 2 --max-mismatches 4 --min-identity 0.9 --max-fragment-span-bp 800
cargo run --bin gentle_cli -- cutrun interpret grch38_tp73 --dataset tp73_release_smoke_synthetic --catalog assets/cutrun.json --cache-dir data/cutrun --flank-bp 0 --report-id tp73_cutrun_reads_from_dataset
cargo run --bin gentle_cli -- cutrun list-read-reports grch38_tp53
cargo run --bin gentle_cli -- cutrun show-read-report tp53_cutrun_reads
cargo run --bin gentle_cli -- cutrun export-coverage tp53_cutrun_reads exports/tp53_cutrun.coverage.tsv --kind coverage
cargo run --bin gentle_cli -- cutrun export-coverage tp53_cutrun_reads exports/tp53_cutrun.cuts.tsv --kind cut_sites
cargo run --bin gentle_cli -- cutrun export-coverage tp53_cutrun_reads exports/tp53_cutrun.fragments.tsv --kind fragments
cargo run --bin gentle_cli -- cutrun inspect-regulatory-support grch38_tp73 --dataset tp73_release_smoke_synthetic --read-report tp73_cutrun_reads_from_dataset --catalog assets/cutrun.json --cache-dir data/cutrun --neighbor-window-bp 150 --species-filter human --path exports/tp73_cutrun.regulatory_support.json
cargo run --bin gentle_cli -- cutrun gene-set-regulatory-support "Human GRCh38 Ensembl 116" --group regulation_of_alternative_splicing --relationship co-regulated --dataset toy_ctcf --path exports/splicing_genes.cutrun_support.json
cargo run --bin gentle_cli -- helpers list
cargo run --bin gentle_cli -- helpers validate-catalog
cargo run --bin gentle_cli -- helpers update-ensembl-specs --catalog assets/helper_genomes.json --output-catalog exports/helper_genomes.updated.json
cargo run --bin gentle_cli -- helpers status "Plasmid pUC19 (online)"
cargo run --bin gentle_cli -- helpers prepare "Plasmid pUC19 (online)" --cache-dir data/helper_genomes --timeout-secs 600
cargo run --bin gentle_cli -- helpers remove-prepared "Plasmid pUC19 (online)" --cache-dir data/helper_genomes
cargo run --bin gentle_cli -- helpers genes "Plasmid pUC19 (online)" --filter bla --limit 20
cargo run --bin gentle_cli -- helpers blast "Plasmid pUC19 (online)" ACGTACGTACGT --task blastn-short --max-hits 10 --options-json '{"thresholds":{"min_identity_percent":95.0}}' --cache-dir data/helper_genomes
cargo run --bin gentle_cli -- cache inspect --references --cache-dir data/genomes
cargo run --bin gentle_cli -- cache clear derived-indexes-only --references --cache-dir data/genomes --prepared-id "Human GRCh38 Ensembl 116"
cargo run --bin gentle_cli -- cache clear selected-prepared --both --cache-dir data/genomes --cache-dir data/helper_genomes --prepared-path data/helper_genomes/localproject
cargo run --bin gentle_cli -- cache clear all-prepared-in-cache --both --include-orphans
cargo run --bin gentle_cli -- candidates generate sgrnas chr1_window --length 20 --step 1 --feature-kind gene --max-distance 500 --limit 5000
cargo run --bin gentle_cli -- candidates score sgrnas gc_balance "100 * (gc_fraction - at_fraction)"
cargo run --bin gentle_cli -- candidates score-weighted sgrnas priority --term gc_fraction:0.7:max --term distance_to_seq_start_bp:0.3:min --normalize
cargo run --bin gentle_cli -- candidates top-k sgrnas sgrnas_top --metric priority --k 100 --direction max
cargo run --bin gentle_cli -- candidates pareto sgrnas sgrnas_front --objective gc_fraction:max --objective distance_to_seq_start_bp:min --max-candidates 200
cargo run --bin gentle_cli -- candidates filter sgrnas sgrnas_q95 --metric gc_balance --min-quantile 0.95
cargo run --bin gentle_cli -- candidates macro @candidate_flow.gsh
cargo run --bin gentle_cli -- candidates macro --transactional --file candidate_flow.gsh
cargo run --bin gentle_cli -- candidates template-put scan_tp53 --script 'generate ${set_name} ${seq_id} --length ${len} --step 1' --param set_name --param seq_id=grch38_tp53 --param len=20
cargo run --bin gentle_cli -- candidates template-run scan_tp53 --bind set_name=tp53_candidates --transactional
cargo run --bin gentle_cli -- guides list
cargo run --bin gentle_cli -- guides put tp73_guides --json '[{"guide_id":"g1","seq_id":"tp73","start_0based":100,"end_0based_exclusive":120,"strand":"+","protospacer":"GACCTGTTGACGATGTTCCA","pam":"AGG","nuclease":"SpCas9","cut_offset_from_protospacer_start":17}]'
cargo run --bin gentle_cli -- guides filter tp73_guides --config '{"gc_min":0.3,"gc_max":0.7,"avoid_u6_terminator_tttt":true}' --output-set tp73_guides_pass
cargo run --bin gentle_cli -- guides oligos-generate tp73_guides lenti_bsmbi_u6_default --apply-5prime-g-extension --output-oligo-set tp73_lenti --passed-only
cargo run --bin gentle_cli -- guides oligos-export tp73_guides exports/tp73_guides.csv --format csv_table --oligo-set tp73_lenti
cargo run --bin gentle_cli -- guides protocol-export tp73_guides exports/tp73_guides.protocol.txt --oligo-set tp73_lenti
cargo run --bin gentle_cli -- shell 'macros template-list'
cargo run --bin gentle_cli -- shell 'macros template-import assets/cloning_patterns.json'
cargo run --bin gentle_cli -- shell 'macros template-import assets/cloning_patterns_catalog'
cargo run --bin gentle_cli -- shell 'macros template-run gibson_two_fragment_overlap_preview --bind left_seq_id=gibson_left --bind right_seq_id=gibson_right --bind overlap_bp=20 --bind assembly_prefix=gibson_demo --bind output_id=gibson_demo_forward --validate-only'
cargo run --bin gentle_cli -- shell 'macros template-run gibson_two_fragment_overlap_preview --bind left_seq_id=gibson_left --bind right_seq_id=gibson_right --bind overlap_bp=20 --bind assembly_prefix=gibson_demo --bind output_id=gibson_demo_forward --transactional'
cargo run --bin gentle_cli -- shell 'macros run --transactional --file cloning_flow.gsh'
cargo run --bin gentle_cli -- shell 'set-param vcf_display_pass_only true'
cargo run --bin gentle_cli -- shell 'set-param vcf_display_required_info_keys ["AF","DP"]'
cargo run --bin gentle_cli -- shell 'set-param tfbs_display_min_llr_quantile 0.95'
cargo run --bin gentle_cli -- shell 'set-param show_restriction_enzymes true'
cargo run --bin gentle_cli -- shell 'set-param restriction_enzyme_display_mode "preferred_only"'
cargo run --bin gentle_cli -- shell 'set-param preferred_restriction_enzymes ["EcoRI","BamHI"]'
cargo run --bin gentle_cli -- shell 'panels import-isoform grch38_tp53 assets/panels/tp53_isoforms_v1.json --panel-id tp53_isoforms_v1'
cargo run --bin gentle_cli -- shell 'panels inspect-isoform grch38_tp53 tp53_isoforms_v1'
cargo run --bin gentle_cli -- shell 'panels render-isoform-svg grch38_tp53 tp53_isoforms_v1 exports/tp53_isoform_architecture.svg'
cargo run --bin gentle_cli -- shell 'panels render-isoform-svg grch38_tp53 tp53_isoforms_v1 exports/tp53_isoform_expression.svg --expression-tsv expression.tsv'
cargo run --bin gentle_cli -- shell 'panels validate-isoform assets/panels/tp53_isoforms_v1.json --panel-id tp53_isoforms_v1'
cargo run --bin gentle_cli -- inspect-feature-expert grch38_tp53 isoform tp53_isoforms_v1
cargo run --bin gentle_cli -- render-feature-expert-svg grch38_tp53 isoform tp53_isoforms_v1 exports/tp53_isoform_architecture.svg
```

For `grna_anchor_window_scan`, the shared macro contract uses
`--bind seq_id=ID --bind anchor_a_pos=0 --bind anchor_b_pos=END`.
The Routine Assistant maps its catalog ports `anchor_a`/`anchor_b` explicitly
to those parameter names. Preflight (`--validate-only`) validates DNA kind,
zero-based boundaries within the sequence, distinct endpoints and sufficient
span for `length_bp` (default 20). Invalid bindings also block actual execution.
These templates perform generic candidate preselection, not PAM-aware guide
design. `grna_practical_filter_and_oligos` instead requires an existing
`guide_set_id`; a generic candidate set is not a guide set.

You can pass JSON from a file with `@file.json` or a bare existing file path.
Bare path-like arguments (containing `/` or `\`, or ending in `.json`) report
file-read errors when missing, rather than misleading inline-JSON syntax errors.
Use `@` to make file intent unambiguous. In GENtle shared shell lines, quote
paths with spaces or backslashes, for example `--resolution 'C:\GENtle data\resolution.json'`.
Backslash escapes the next character outside single quotes, including inside
double quotes; programmatic callers should use the shared `quote_shell_arg`
helper, not interpolate native paths. This is not the native OS shell grammar.
When loading from file path, an initial shebang line (`#!...`) is ignored so
executable script files can embed JSON payloads directly.
`workflow` accepts both raw workflow payloads (`{"run_id":"...","ops":[...]}`)
and wrapped protocol example payloads (`{"workflow":{...}}`).

<!-- manual-page-nav:start -->

---

[CLI Manual](../cli.md) | [Previous: Protocol Examples and Tutorials](examples.md) | [Next: Operation JSON Examples](operation-examples.md)
<!-- manual-page-nav:end -->
