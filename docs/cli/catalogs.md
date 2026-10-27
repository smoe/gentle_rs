# Services and Resource Catalogs

<!-- manual-counterparts:start -->
Related GUI reference: [Catalog Inspection And Provider Imports](../gui/catalog-inspection-and-provider-imports.md) | [External Services Window](../gui/external-services-window.md).
These topic links do not imply identical feature coverage across interfaces.
<!-- manual-counterparts:end -->

Catalog discovery and guarded local provider updates use the same shared shell:
`catalog doctor`, `catalog search --entity product|provider|protocol --query TEXT`,
`catalog updates preview BUNDLE_PATH --scope user|project --output PREVIEW_PATH`,
and `catalog updates apply PREVIEW_PATH --approve-digest SHA256`.
See [the contract](../catalog_updates.md) and
[tutorial 09.04](../tutorial/09-04_cloning_time_cost_and_protocol_suppliers.md#61-review-a-disposable-supplier-update-without-changing-vendor-facts).
Imports change local configuration, not project state or scientific results;
they never contact a supplier. Product search includes helper/reference rows,
not just commercially available products.



## Catalog Discovery and Paths

Catalog path note:

- `--catalog PATH` accepts either one JSON file or one directory of `*.json`
  catalog fragments.
- Directory fragments are merged in sorted filename order; duplicate entry ids
  currently fail fast instead of silently overriding each other.
- When `--catalog` is omitted, GENtle now uses a deterministic discovery chain
  instead of one hard-coded file:
  - built-in `assets/*.json` plus optional built-in `assets/*.d/`
  - system overlays under `/etc/gentle/catalogs/`
  - user overlays under `$XDG_CONFIG_HOME/gentle/catalogs/` or
    `~/.config/gentle/catalogs/`
  - project overlays under `PROJECT_ROOT/.gentle/catalogs/`
- Built-in/system/project roots can be overridden in controlled environments
  via `GENTLE_ASSET_ROOT`, `GENTLE_SYSTEM_CONFIG_ROOT`, and
  `GENTLE_PROJECT_ROOT`.
- `genomes list` and `helpers list` now also accept `--filter TEXT` and search
  ids plus catalog metadata such as aliases, tags, search terms, and the new
  helper semantic/procurement fields.
- Helper catalogs may now carry richer metadata (`summary`, `aliases`, `tags`,
  `search_terms`, `helper_kind`, `host_system`, `procurement`, `semantics`)
  while still using the same preparation/indexing pipeline as before.
- CUT&RUN V1/V2 now use the same discovery pattern:
  - built-in starter catalog: `assets/cutrun.json`
  - built-in catalog shards under `assets/cutrun.d/`, including the Rostock
    p73 `E-MTAB-15709` SRA-backed paired-end run entries
  - project/system overlays under `.gentle/catalogs/` and `/etc/gentle/catalogs/`
  - prepared-cache root defaults to `data/cutrun`
  - override cache root with `GENTLE_CUTRUN_CACHE_DIR`
  - built-in dataset `tp73_release_smoke_synthetic` is a tiny offline V1-V3
    integration fixture; it is explicitly synthetic and carries no biological
    occupancy claim
  - `cutrun prepare` and `cutrun status` now use the same shared lifecycle
    contract as reference/helper prepares:
    - `lifecycle_status = missing|running|ready|failed|cancelled|stale`
    - `current_activity` reports the active/stale dataset prepare marker when
      present
    - duplicate `cutrun prepare` callers reuse one active dataset install
      instead of launching a second parallel materialization
  - `cutrun project` only accepts genome-anchored sequences and reuses the
    shared anchored BED/BigWig track import behavior
  - `cutrun interpret` maps ad hoc `FASTA`/`FASTQ` reads only against one
    selected genome-anchored ROI plus deterministic flanks and stores a shared
    engine-owned read report for later inspection/export
  - with `--flank-bp 0`, an imported sequence whose anchor span exactly matches
    its sequence length is itself the mapping reference; this offline path does
    not require a prepared whole genome, while nonzero flanks still do
  - the same `cutrun interpret` route can also resolve prepared raw reads from
    `cutrun prepare` via `--dataset DATASET_ID [--catalog PATH] [--cache-dir PATH]`
  - catalog entries may declare raw reads with `reads_sra_accession`;
    `cutrun prepare` runs the shared read-acquisition layer first, then writes
    the existing prepared manifest shape for `cutrun interpret --dataset ...`
  - `cutrun export-coverage` writes TSV summaries for per-base coverage,
    cut-site density, or derived fragment spans from one saved CUT&RUN read
    report
  - `cutrun inspect-regulatory-support` aggregates prepared datasets and/or
    saved ROI read reports into one promoter/ROI support report with:
    - `--catalog PATH` and `--cache-dir PATH` for replaying the same custom V1
      dataset location used during prepare/project
    - merged `strong`/`moderate`/`weak` support windows
    - theoretical TFBS rows split into legacy `confirmed` vs `unconfirmed`
      vectors, with additive per-row `support_status` values:
      `confirmed`, `nearby`, `absent`, or `motif_poor`
    - motif-absent strong windows classified as
      `context_supported_by_other_motifs` or `motif_poor_supported`
    - recurring motif-context summaries across many motif-absent supported
      windows
  - `cutrun gene-set-regulatory-support` scores prepared datasets and saved
    read reports across a gene-set promoter cohort and reports evaluated versus
    unevaluated members separately; optional `--relationship` records a
    declared expectation and emits non-blocking CUT&RUN occupancy relationship
    flags over evaluated members only.


Resource update capability status:

- `gentle_cli`: supported (`resources sync-rebase`, `resources sync-jaspar`, `resources sync-ucsc-rmsk`, `resources import-gene-list-cache`, `resources import-ontology-assignment-cache`, `resources import-co-regulated-cache`, `resources install-ucsc-rmsk`, `resources prepare-ucsc-rmsk-index`, `resources suggest-ucsc-rmsk-index`, `resources list-encode-ccre-sources`, `resources install-encode-ccres`, `resources prepare-encode-ccre-index`, `resources list-ensembl-regulation-sources`, `resources install-ensembl-regulatory-features`, `resources prepare-ensembl-regulation-index`, `resources sync-jaspar-remote-metadata`, `resources summarize-jaspar`, `resources resolve-tf-query`, `gene-groups list/show/resolve/doctor/draft`, `gene-sets resolve/produce direct-list/produce ontology-assignment/produce co-regulated/promoter-cohort/regulatory-partner-screen`, `resources benchmark-jaspar`, `resources list-jaspar`, `resources list-publication-datasets`, `resources status-publication-dataset`, `resources prepare-publication-dataset`, `resources inspect-jaspar`)
- `gentle_js`: supported (`sync_rebase`, `sync_jaspar`)
- `gentle_lua`: supported (`sync_rebase`, `sync_jaspar`)


- Shared shell:
  - Commands:
    - `resources sync-rebase INPUT.withrefm_or_URL [OUTPUT.rebase.json] [--commercial-only]`
    - `resources sync-jaspar INPUT.jaspar_or_URL [OUTPUT.motifs.json]`

## Repeat and Gene-List Resources

- Shared shell:
  - Commands:
    - `resources sync-ucsc-rmsk INPUT.rmsk.txt_or_txt.gz [OUTPUT.rmsk.json] [--assembly DB] [--limit N]`
    - `resources import-gene-list-cache --input INPUT.tsv_or_csv --output CACHE.json --provider PROVIDER --version VERSION [--organism NAME|--taxon-id N|--namespace NAMESPACE]`
    - `resources import-ontology-assignment-cache --input INPUT.tsv_or_csv --output CACHE.json --provider PROVIDER --version VERSION [--ontology-namespace GO] [--organism NAME|--taxon-id N|--symbol-namespace NAMESPACE]`
    - `resources import-co-regulated-cache --input INPUT.tsv_or_csv --output CACHE.json --provider PROVIDER --version VERSION [--dataset DATASET_ID] [--contrast LABEL] [--score METHOD] [--organism NAME|--taxon-id N|--namespace NAMESPACE]`
    - `resources install-ucsc-rmsk [--assembly DB] [--input PATH_OR_URL] [--resource-output PATH] [--index-output PATH]`
    - `resources prepare-ucsc-rmsk-index RESOURCE.rmsk.json [OUTPUT.interval-index.json]`
    - `resources suggest-ucsc-rmsk-index [--assembly DB] [--output OUTPUT.json]`
    - `resources list-encode-ccre-sources [--assembly GRCh38|hg38|GRCm38|mm10] [--output OUTPUT.json]`
    - `resources install-encode-ccres SOURCE_ID [--input PATH_OR_URL] [--bed-output PATH] [--index-output PATH]`
    - `resources prepare-encode-ccre-index SOURCE_ID BED_PATH [--output INDEX.json]`
    - `resources list-ensembl-regulation-sources [--assembly GRCh38|hg38|GRCm39|mm39] [--output OUTPUT.json]`
    - `resources install-ensembl-regulatory-features SOURCE_ID [--input COMPLETE_API_RESPONSE.json] [--intervals-output PATH] [--index-output PATH]`
    - `resources prepare-ensembl-regulation-index SOURCE_ID INTERVALS_PATH [--output INDEX.json]`
    - `resources sync-jaspar-remote-metadata [--motif TOKEN ...] [--motifs CSV] [--all] [--filter TOKEN] [--limit N] [--output OUTPUT.json]`

## JASPAR Inspection and TF Queries

- Shared shell:
  - Commands:
    - `resources summarize-jaspar [--motif TOKEN ...] [--motifs CSV] [--all] [--random-length N] [--seed N] [--output OUTPUT.json]`
    - `resources resolve-tf-query QUERY [QUERY ...] [--output OUTPUT.json]`

- Shared shell:
  - Commands:
    - `resources benchmark-jaspar [--random-length N] [--seed N] [--output OUTPUT.json]`
    - `resources list-jaspar [--filter TOKEN] [--limit N] [--fetch-remote] [--output OUTPUT.json]`
    - `resources inspect-jaspar MOTIF [--random-length N] [--seed N] [--fetch-remote] [--output OUTPUT.json]`


## Service Readiness and Provider Routing

Service and resource commands:

- `services status`
  - Reports one combined service-readiness snapshot for common GENtle-backed
    analysis flows.
  - Covers canonical reference genomes (currently `Human GRCh38 Ensembl 116`),
    helper backbones (currently `Plasmid pUC19 (online)`), and active external
    resource snapshots in one record.
  - Also embeds the external-service provider catalog under
    `external_providers`; provider behavior is loaded from the overlay
    catalog rooted at `assets/external_service_providers.json` and matching
    system/user/project override directories.
  - Reference/helper rows now expose canonical lifecycle fields:
    - `resource_key`
    - `display_name`
    - `prepared`
    - `lifecycle_status = missing|running|ready|failed|cancelled|stale`
    - `current_activity` (when present)
  - `availability_status` remains as the chat/backward-compat alias
    (`not_prepared|preparing|prepared|failed|cancelled|stale`).
  - When a prepare/reindex run is active, it surfaces the live dependency
    activity plus the current phase/progress hints instead of only saying
    `prepared`/`not_prepared`.
  - Reference/helper prepare is now duplicate-safe at the shared-target level:
    one semantic target may actively prepare at a time, stale heartbeat markers
    are superseded safely, and follow-up callers see `running` instead of
    launching redundant work.
  - Includes compact `summary_lines` that chat/report layers such as
    ClawBio/OpenClaw can surface directly before deciding whether they should
    fetch, prepare, or render.
- `services providers list`
  - Returns `gentle.external_service_provider_catalog.v1`.
  - V1 catalogs `provider = geneart` plus configured overlay providers such as
    `provider = metabion`.
  - GeneArt rows include DNA fragments, cloned genes, plasmid reorder,
    mutagenesis, and protein-expression capabilities.
  - Metabion rows initially include DNA single-tube oligos and m-block DNA
    fragments/libraries as WOP/email+Excel handoff routes.
  - Direct API/order/portal submission is not enabled by this command.
  - Source-backed discovery-only suppliers additionally expose `country_code`,
    `country_source_url` and dated `protocols[]` with `html`/`text_pdf` formats.
    Filter the returned JSON, not website language or guessed delivery coverage.
    Country describes the named company/operation; unknown legacy metadata is
    absent. The global directory is non-exhaustive and does not enable quotes,
    protocol execution, downloads or automatic supplier selection. See
    [tutorial 09.04](../tutorial/09-04_cloning_time_cost_and_protocol_suppliers.md).
- `services providers doctor [--catalog PATH] [--output PATH]`
  - Returns `gentle.external_service_provider_config_doctor.v1`.
  - Validates the external-service provider config catalog discovery chain, or
    one explicit catalog file via `--catalog`.
  - Reports source scope/path, provider counts, SHA1 checksums, and
    deterministic schema/required-field errors.
  - With `--output`, also writes the JSON doctor report to disk.
  - Protocol metadata checks require uppercase two-letter country codes with
    HTTPS source URLs, readable-format enums, HTTPS protocol references and
    valid calendar review dates; this is local validation, not a live URL,
    corporate domicile or protocol-content audit.
- `services delivery-route REQUEST_JSON_OR_@FILE`
  - Parses `gentle.external_service_delivery_route_request.v1` and returns
    `gentle.external_service_delivery_route.v1`.
  - This is the generic "deliver this sequence" classifier before any vendor
    preflight or quote handoff. It uses source kind, sequence alphabet, line
    items, length hints, and vector/construct context to recommend a
    provider/service-kind candidate.
  - Current deterministic routing rules are intentionally conservative:
    short DNA oligos/primers/probes route to Metabion
    `dna_oligo_single_tube`; medium DNA fragments route to Metabion
    `dna_fragment` m-block handoff; vector-backed/cloned-gene DNA routes to
    GeneArt `cloned_gene`; long synthetic DNA without vector context routes to
    GeneArt `dna_fragment`; protein/amino-acid products route to GeneArt
    `protein_expression`.
  - RNA, missing sequence context, or unclear molecule type returns
    `needs_clarification` instead of guessing.
  - Returned candidate requests are ordinary
    `gentle.external_service_request.v1` objects suitable for human-reviewed
    `services project-preflight` or `services project-quote`.
- `services route-project-source --kind sequence|oligo-form|primer-report-rows ...`
  - Builds a delivery-route request from a real project object and immediately
    classifies it with the same logic as `services delivery-route`.
  - Sequence source:
    `services route-project-source --kind sequence --seq-id ID [--range START..END] [--as sequence|construct-output]`.
    Ranges are 0-based, end-exclusive. `--as construct-output` keeps the
    selected sequence as the source text but marks `source_target.kind =
    construct_output`, so construct-aware routing can recommend cloned-gene
    handoff.
  - Oligo form source:
    `services route-project-source --kind oligo-form --form-id FORM_ID`.
    The route request reuses the persisted oligo order form line items,
    duplicate groups, sequence-reuse groups, and duplicate-review status.
  - Primer report rows:
    `services route-project-source --kind primer-report-rows --report-id REPORT_ID --pair-rank N[,N...] [--form-id FORM_ID]`.
    GENtle first persists those pair ranks as an oligo order form, then routes
    the form. This preserves duplicate-review semantics; no transient order
    forms are used in v1.
  - Classification is allowed before review, but `services project-preflight`
    and `services project-quote` block requests whose
    `source_target.duplicate_review.status` is `review_required`.
- `services project-preflight REQUEST_JSON_OR_@FILE`
  - Parses `gentle.external_service_request.v1` and returns
    `gentle.external_service_preflight.v1`.
  - The check is deterministic and local: it validates provider/service
    support, required source-target fields, quote-handoff availability, and
    whether direct submission is currently available from GENtle.
  - Use `@request.json` for realistic requests, especially when
    `source_target`, `vector_spec`, `delivery_options`, or `return_spec`
    contain nested JSON.
- `services project-quote REQUEST_JSON_OR_@FILE [--output-dir DIR]`
  - Builds `gentle.external_service_quote.v1` with dashboard/form links and an
    inline service-ready handoff bundle.
  - No vendor order, cart submission, API upload, credential lookup, or
    purchase/shipping action is performed in this slice.
  - Configured providers add deterministic handoff payloads. For Metabion,
    quote output includes normalized line-item JSON/CSV, an email-draft
    markdown payload, and a guided WOP checklist; missing local vendor
    templates are warnings, not hard failures.
  - With `--output-dir DIR`, the same inline payloads are written as local
    bundle files and the returned quote report records them in
    `service_ready_bundle.local_files[]` alongside `quote_report.json`.
  - `return_spec.requested_payloads[]` lets ClawBio/MCP callers state whether
    they want GenBank/FASTA/protein payloads, quote metadata, vendor ids, or a
    full artifact bundle as later implementations grow.

- `services guide --channel telegram [--section overview|readiness|gene-context|tfbs|inline-dna|cloning|isoforms|follow-up] [--gene SYMBOL]`
  - Builds a Telegram-compatible bench-user guide with schema
    `gentle.telegram_guide.v1`.
  - The guide includes compact `summary_lines[]`, `readiness_summary_lines[]`,
    `menu_sections[]`, `blocked_actions[]`, and `suggested_actions[]`.
  - `summary_lines[]` also includes a text fallback menu with phrases such as
    `Continue readiness`, `Continue cloning`, and `Continue isoforms`, so chat
    adapters that do not yet render `suggested_actions[]` can still let users
    navigate the guide by sending another message.
  - Its `suggested_actions[]` entries act like navigation links:
    `kind = guide_section`, `requires_confirmation = false`, and the
    `shell_line` points back to another guide section.
  - When `--gene SYMBOL` is supplied, gene-aware sections personalize example
    actions to that symbol. Without a gene, GENtle uses section defaults such
    as `TERT`/`TP73` for promoter-TFBS and `TP73`/`TP53` for isoform demos.
  - This route is for bench-user Telegram orientation. Operator/setup status
    and confirmation-gated prepare/sync actions still come from
    `services status` and `services handoff`.
- `services handoff [--scope NAME] [--output PATH]`
  - Builds a chat-gateway handoff/doctor report with schema
    `gentle.service_handoff.v1`.
  - The report embeds the `services status` readiness payload and adds:
    `readiness[]`, `suggested_actions[]`, `running_actions[]`,
    `blocked_actions[]`, `preferred_demo_actions[]`, `preferred_artifacts[]`,
    `status_overview`, `environment_hints[]`, `warnings[]`, and
    `summary_lines[]`.
  - `status_overview` gives compact lifecycle counts, an overall
    `ready|setup_needed|setup_running|attention_needed` state, and one
    recommended next action for chat gateways that want a single primary
    button.
  - Suggested actions are deterministic GENtle shell commands for missing or
    retryable shared setup. Running shared prepares become refresh/status
    suggestions, not duplicate prepare commands.
  - Blocked actions are intentionally not offered as immediately executable
    commands. For example, ATtRACT sync is blocked until the caller supplies a
    local `ATtRACT.zip` path, then the command can be completed as
    `resources sync-attract /path/to/ATtRACT.zip`.
  - When `--output` is provided, the same JSON payload is written to disk and
    listed as the best-first `preferred_artifacts[]` entry so ClawBio can bundle
    it as an auditable installation report.

Tutorial companion:

- `docs/tutorial/09-01_metabion_external_service_handoff_gui_cli.md` walks through
  provider doctor, capability listing, preflight, quote-handoff generation, and
  GUI parity using bundled synthetic Metabion oligo and m-block request JSON.

## Resource Status

- `resources status`
  - Reports which integrated external resource snapshots are active right now.
  - Covers built-in/runtime status for `REBASE` and `JASPAR`, including item
    counts and whether GENtle is currently using a runtime override or the
    built-in snapshot.
  - JASPAR status now also reports the persisted remote-metadata snapshot used
    by the JASPAR catalog/expert cache:
    - `path`
    - `exists`
    - `valid`
    - `item_count`
    - `error`
  - Also reports normalized ATtRACT runtime status:
    - `known_external_only` = homepage/ZIP URL known, but no active normalized
      snapshot
    - `runtime_snapshot` = `data/resources/attract.motifs.json` is present and
      active
    - `session_override` = the current session has loaded an ATtRACT snapshot
      from a non-default path
  - ATtRACT status now also exposes deterministic provenance for the active and
    runtime snapshots:
    - `runtime_fingerprint`
    - `active_fingerprint`
  - Includes the current ATtRACT homepage plus published ZIP download URL:
    `https://attract.cnic.es/attract/static/ATtRACT.zip`.
  - UCSC `rmsk` status reports the default hg38 download/schema URLs, the
    default normalized snapshot path `data/resources/ucsc.rmsk.hg38.json`,
    whether that snapshot exists and validates, and the embedded index
    recommendations used by `resources suggest-ucsc-rmsk-index`.
  - Publication-associated datasets are tracked through
    `assets/publication_resources.json` and default prepared manifests under
    `data/publication_resources`:
    - Rostock p73 Clariom D transcriptome data: `E-MTAB-14704`
    - Rostock p73 CUT&RUN sequencing data: `E-MTAB-15709` / `PRJEB100610`
      with ENA paired FASTQ URLs declared as `raw_cutrun_fastq`
    - Rostock p73 Co-IP/MS proteomics data: `PXD058816`
    - `resources prepare-publication-dataset` writes `manifest.json`,
      `download_manifest.tsv`, and `download.sh`; declared sizes and optional
      archive MD5 checksums are preserved, and large file download happens only
      with `--download-files`.
    - The built-in Rostock p73 entries already ship with those no-byte-download
      manifests and scripts prepared.
  - RNA secondary-structure executable resources are reported as first-class
    readiness entries:
    - `vienna_rna`: resolves `GENTLE_RNAFOLD_BIN` or `RNAfold` and probes
      `RNAfold --version`
    - `rnapkin`: resolves `GENTLE_RNAPKIN_BIN` or `rnapkin` and probes
      `rnapkin --version` / `rnapkin -V`
    - these are executable tool resources, not normalized data snapshots;
      GENtle uses ViennaRNA/RNAfold for folding and rnapkin for rendering
  - Legacy SHA-1 verification is reported as `legacy_sha1`:
    - resolves `GENTLE_SHA1_TOOL` first, then platform defaults
      (`sha1sum`, `shasum -a 1`, and Windows `certutil`)
    - set `GENTLE_DISABLE_LEGACY_SHA1=1` to skip SHA-1 checks intentionally
    - missing tools do not block downloads by themselves; GENtle reports that
      SHA-1 verification was skipped and relies on size/gzip/parser checks
      where applicable
- `resources sync-rebase INPUT.withrefm [OUTPUT.rebase.json] [--commercial-only]`
  - Parses REBASE/Bairoch-style records (`withrefm`) into GENtle restriction-enzyme JSON.
  - `INPUT` may be a local file path or an `https://...` URL.
  - Default output: `data/resources/rebase.enzymes.json`.
  - At runtime, GENtle auto-loads this file (if present) and overrides embedded restriction enzymes.
- `resources sync-jaspar INPUT.jaspar.txt [OUTPUT.motifs.json]`
  - Parses JASPAR PFM text into motif snapshot JSON with IUPAC consensus.
  - `INPUT` may be a local file path or an `https://...` URL.
  - Default output: `data/resources/jaspar.motifs.json`.
  - `ExtractAnchoredRegion` can resolve TF motifs by ID/name from this local registry.
  - GENtle ships with built-in motifs in `assets/jaspar.motifs.json` (currently generated from JASPAR 2026 CORE non-redundant); this command provides local updates/extensions.
  - When a compact built-in entry lacks PFM rows, GENtle supplements it from
    the bundled full-PFM table where an unambiguous ID/name match exists. This
    keeps score tracks and sequence logos probabilistic instead of collapsing
    ambiguous positions into a fake consensus.

## Publication Datasets

- `resources list-publication-datasets [--filter TEXT] [--catalog PATH] [--output OUTPUT.json]`
  - Lists literature-associated external datasets from
    `assets/publication_resources.json` or an override catalog.
  - The built-in catalog currently describes the Rostock p73 multimodal paper
    datasets and records which accessions have concrete downloadable file URLs.
- `resources status-publication-dataset DATASET_ID [--catalog PATH] [--cache-dir PATH]`
  - Inspects the local prepared-manifest/cache state for one publication
    dataset id or accession.
  - Reports concrete file URLs, local paths, size checks when known, and
    warnings for accessions whose public archive file list still needs
    refreshing. The Rostock CUT&RUN entry now reports the ENA FASTQ files but
    does not download them unless asked.
- `resources prepare-publication-dataset DATASET_ID [--catalog PATH] [--cache-dir PATH] [--download-files] [--max-files N] [--category NAME|--categories CSV]`
  - Writes a deterministic `manifest.json`, `download_manifest.tsv`, and
    executable `download.sh` under `data/publication_resources/DATASET_ID`.
    `download_manifest.tsv` includes declared file sizes plus optional
    `checksum_md5` values when the source archive reports them.
  - Without `--download-files`, no raw bytes are fetched; this is the safe
    default for multi-GB PRIDE and sequencing datasets.
  - With `--download-files`, GENtle downloads the declared file URLs directly
    into the prepared dataset directory. `--max-files N` limits both manifests
    and downloads for smoke checks.
  - `--category`/`--categories` restrict the planned files before download;
    for the Rostock Clariom D data, `raw_microarray` selects just the CEL
    payloads and `metadata` selects the MAGE-TAB files. For the Rostock
    CUT&RUN data, `raw_cutrun_fastq` selects the public paired ENA FASTQs.
- `resources sync-ucsc-rmsk INPUT.rmsk.txt_or_txt.gz [OUTPUT.rmsk.json] [--assembly DB] [--limit N]`
  - Normalizes the UCSC RepeatMasker `rmsk` table into
    `gentle.ucsc_rmsk_resource.v1`.
  - `INPUT` may be a local path or an `https://...` URL and may be gzip
    compressed.
  - Default assembly: `hg38`.
  - Default output for hg38: `data/resources/ucsc.rmsk.hg38.json`.
  - The snapshot preserves UCSC 0-based half-open coordinates and the canonical
    17-column row fields (`genoName`, `genoStart`, `genoEnd`, `repName`,
    `repClass`, `repFamily`, etc.).
  - `--limit N` is intended for fixtures/smoke checks; a limited output is
    marked `truncated=true` and should not be used as a full assembly resource.
- `resources import-gene-list-cache --input INPUT.tsv_or_csv --output CACHE.json --provider PROVIDER --version VERSION [--organism NAME|--taxon-id N|--namespace NAMESPACE] [--cache-id ID] [--cache-version VERSION] [--list-id ID] [--list-label LABEL]`
  - Normalizes local TSV/CSV rows with `symbol` or `gene_id` columns into
    `gentle.gene_set_direct_list_cache.v1`.
  - Optional `list_id`/`list_label` columns group rows into named lists; CLI
    `--list-id` and `--list-label` provide defaults when the input omits them.
  - The route is local and offline-only. It does not resolve genes by itself;
    feed the emitted cache to `gene-sets produce direct-list` for a resolved
    `gentle.gene_set_resolution.v1` artifact.
- `resources import-ontology-assignment-cache --input INPUT.tsv_or_csv --output CACHE.json --provider PROVIDER --version VERSION [--ontology-namespace GO] [--organism NAME|--taxon-id N|--symbol-namespace NAMESPACE] [--cache-id ID] [--cache-version VERSION]`
  - Normalizes local TSV/CSV ontology assignment rows into
    `gentle.gene_set_ontology_assignment_cache.v1`.
  - Input rows should include a term column such as `term`, `term_id`, or
    `go_id` plus `symbol` or `gene_id`; optional evidence columns such as
    `evidence_code` and `assigned_by` are preserved for later filters.
  - The route is local and offline-only. It does not download GO annotations;
    feed the emitted cache to `gene-sets produce ontology-assignment`.
- `resources import-co-regulated-cache --input INPUT.tsv_or_csv --output CACHE.json --provider PROVIDER --version VERSION [--dataset DATASET_ID] [--contrast LABEL] [--condition LABEL] [--normalization METHOD] [--score METHOD] [--organism NAME|--taxon-id N|--namespace NAMESPACE] [--cache-id ID] [--cache-version VERSION]`
  - Normalizes local scored gene rows into
    `gentle.gene_set_co_regulated_cache.v1`.
  - Input rows should include `symbol` or `gene_id` plus the score column named
    by `--score` or a generic `score`/`value`/`logfc`/`statistic` column.
  - The route records dataset/contrast/condition metadata for later
    thresholding but does not claim the genes are regulated. Feed the emitted
    cache to `gene-sets produce co-regulated`.
- `resources install-ucsc-rmsk [--assembly DB] [--input PATH_OR_URL] [--resource-output PATH] [--index-output PATH]`
  - One-step installer that normalizes the UCSC `rmsk` table and builds the
    matching interval index sidecar.
  - With no `--input`, GENtle downloads the UCSC default
    `goldenPath/<assembly>/database/rmsk.txt.gz` source.
  - Default hg38 outputs are `data/resources/ucsc.rmsk.hg38.json` and
    `data/resources/ucsc.rmsk.hg38.interval-index.json`.
  - Limited installs are rejected for runtime use so fixture snapshots cannot be
    mistaken for full assembly resources.
- `resources prepare-ucsc-rmsk-index RESOURCE.rmsk.json [OUTPUT.interval-index.json]`
  - Builds `gentle.ucsc_rmsk_interval_index.v1` from a complete normalized
    resource snapshot.
  - Default hg38 output:
    `data/resources/ucsc.rmsk.hg38.interval-index.json`.
  - The sidecar groups records by chromosome, stores intervals sorted by
    `(start,end,row_offset)`, and carries class/family and repeat-name lookup
    dictionaries plus chromosome-alias mappings for deterministic overlap and
    display filtering.
  - Truncated snapshots are rejected so smoke-test resources cannot be used as
    full assembly indexes by accident.
- `resources suggest-ucsc-rmsk-index [--assembly DB] [--output OUTPUT.json]`
  - Emits the source URLs, column schema, and indexing recommendations without
    downloading the full table.
  - Current recommendations prioritize:
    - `interval_by_chrom_bin` for extraction-window overlap queries
    - `class_family_partition` for SINE/LINE/LTR/DNA/simple-repeat grouping
    - `repeat_name_lookup` for repeat-name expert searches
    - `quality_metric_columns` for future divergence/score shading
- `resources summarize-jaspar [--motif TOKEN ...] [--motifs CSV] [--all] [--random-length N] [--seed N] [--output OUTPUT.json]`
  - Presents local JASPAR entries through one deterministic score-oriented report.
  - For each selected entry, GENtle derives:
    - a maximizing sequence,
    - a minimizing sequence,
    - and `llr_bits` / `true_log_odds_bits` score distributions on one deterministic pseudorandom DNA background.
  - Omitting `--motif` and `--motifs` summarizes the full active local JASPAR registry.
  - `TOKEN` accepts:
    - exact motif ids or TF names,
    - common aliases such as `OCT4`,
    - catalog-backed functional groups such as `Yamanaka factors` / `stemness`,
    - and family-like queries such as `KLF family`.
  - Default background length: `10000 bp`.
  - Default seed: one built-in deterministic seed so agent/CLI runs remain reproducible.
- `resources resolve-tf-query QUERY [QUERY ...] [--output OUTPUT.json]`
  - Resolves user-facing TF queries into the exact local motif ids/names that
    GENtle will use for downstream TFBS/JASPAR routes.
  - This is the lightweight inspection route for:
    - exact ids/names,
    - aliases such as `OCT4`,
    - catalog-backed functional groups such as `Yamanaka factors` / `stemness`,
    - and family-like queries such as `KLF family`.
  - Output includes per-query resolution rows plus the expanded motif ids so
    ClawBio/CLI users can audit what a natural-language TF query mapped to.
- `resources benchmark-jaspar [--random-length N] [--seed N] [--output OUTPUT.json]`
  - Benchmarks the full active local JASPAR registry through one deterministic
    shared background and writes an export-ready drift snapshot.
  - Default output: `data/resources/jaspar.registry_benchmark.json`.
  - Includes:
    - the full per-entry presentation rows used by `summarize-jaspar`
    - one score-family aggregate summary for `llr_bits`
    - one score-family aggregate summary for `true_log_odds_bits`
    - top motifs by maximum score and by positive random-background hit rate
  - This is the intended routine artifact for regression/drift checks across
    JASPAR updates or scoring changes.
- `resources list-jaspar [--filter TOKEN] [--limit N] [--fetch-remote] [--output OUTPUT.json]`
  - Lists local JASPAR entries as a portable catalog view.
  - Includes per-row:
    - motif id
    - optional TF name
    - consensus IUPAC sequence
    - motif length
    - and, when `--fetch-remote` is requested, a compact remote summary with
      collection/class/family/tax-group fields plus species-count previews for
      the returned subset
  - `--filter` matches motif id, TF name, and consensus.
  - `--limit` caps the returned subset before optional remote lookups.
  - When `--fetch-remote` is requested, GENtle first reuses any matching rows
    from `data/resources/jaspar.remote_metadata.json` and only refreshes
    missing/selected rows live when necessary.
- `resources sync-jaspar-remote-metadata [--motif TOKEN ...] [--motifs CSV] [--all] [--filter TOKEN] [--limit N] [--output OUTPUT.json]`
  - Refreshes and persists reusable JASPAR remote metadata for one selected
    subset of the local registry.
  - Default output: `data/resources/jaspar.remote_metadata.json`.
  - Stored rows include motif id/name, consensus, motif length, full remote
    metadata payload, and a compact summary for catalog views.

## RNA-Binding Protein Resources

- `resources sync-attract INPUT.zip_or_URL [OUTPUT.attract.json]`
  - Normalizes the published ATtRACT ZIP into GENtle’s runtime motif snapshot.
  - Reads stored/deflated ZIP members directly, also on Windows; no separate
    `unzip` installation is needed. Local MaxEnt ZIP imports use the same reader.
  - `INPUT` may be a local ZIP path or an `https://...` URL.
  - Default output: `data/resources/attract.motifs.json`.
  - Current v1 normalization reads `ATtRACT_db.txt`, keeps deterministic
    consensus/IUPAC motif records, and when `pwm.txt` blocks can be mapped by
    `Matrix_id`, stores those PFM rows so splice-aware inspection can add
    PWM-backed ranking on top of the same candidate hits.
- `attract inspect-splicing SEQ_ID FEATURE_ID [--scope ...] [--organism NAME] [--flank-bp N] [--min-score X] [--min-match-quantile Q] [--pwm-mapping strict_same_length|windowed_submatrix] [--compare-policies] [--all-transcripts] [--no-fallback]`
  - Runs the shared engine-owned splice-aware ATtRACT inspection route over the
    selected Splicing Expert group.
  - Default behavior:
    - transcript-strand only
    - exact organism match first
      - optional fallback to all normalized ATtRACT entries when no exact match
      exists
      - region classification returned as `exon_body`, `donor_flank`,
      `acceptor_flank`, or `intron_body`
    - `--min-match-quantile` defaults to `0.99` and applies to PWM-backed rows;
      rows without a mapped PWM block keep exact consensus/IUPAC matching
  - Returns grouped RBP summary rows plus detailed hit rows from the same
    payload the GUI uses.
  - Payload provenance now also includes `active_resource_fingerprint` so CLI,
    GUI, and exported reports can tie the evidence back to the exact normalized
    ATtRACT snapshot content rather than only a motif count.
  - Recommended review pattern:
    - start with `--pwm-mapping strict_same_length`
    - add `--compare-policies` to see what the alternate heuristic would add
      without switching the primary result set
    - only move to `--pwm-mapping windowed_submatrix` when that alternate
      summary looks biologically interesting enough to inspect directly
  - Example:
    - `attract inspect-splicing TP73 17 --organism "Homo sapiens" --pwm-mapping strict_same_length --compare-policies`
- `resources inspect-jaspar MOTIF [--random-length N] [--seed N] [--fetch-remote] [--output OUTPUT.json]`
  - Expands one local JASPAR entry into an expert-oriented report.
  - Includes:
    - the full count matrix,
    - a sequence-logo payload,
    - `llr_bits` and `true_log_odds_bits` score distributions with compact
      histogram bins,
    - and optional remote JASPAR metadata/species assignments when
      `--fetch-remote` is requested.
  - Default background length: `10000 bp`.
  - Default seed: one built-in deterministic seed so repeated expert
    inspection stays reproducible.
  - When `--fetch-remote` is requested, GENtle reuses any matching persisted
    remote metadata row first and only then falls back to a live refresh.

<!-- manual-page-nav:start -->

---

[CLI Manual](../cli.md) | [Previous: Panel Workflows](panels.md) | [Next: Engine Settings](settings.md)
<!-- manual-page-nav:end -->
