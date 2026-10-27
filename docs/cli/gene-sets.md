# Gene Groups and Gene Sets

<!-- manual-counterparts:start -->
Related GUI reference: [Gene Set Inspector And Collection Operations](../gui/gene-set-inspector-and-collection-operations.md).
These topic links do not imply identical feature coverage across interfaces.
<!-- manual-counterparts:end -->

- Shared shell:
  - Commands:
    - `gene-groups list [--catalog PATH] [--filter TEXT] [--output OUTPUT.json]`
    - `gene-groups show GROUP_ID [--catalog PATH] [--output OUTPUT.json]`
    - `gene-groups resolve TOKEN [--catalog PATH] [--output OUTPUT.json]`
    - `gene-groups doctor [--catalog PATH] [--output OUTPUT.json]`
    - `gene-groups draft --description TEXT [--member SYMBOL] [--candidate SYMBOL=EVIDENCE] [--go GO:NNNNNNN] [--output GROUP.json]`
    - `gene-sets resolve [GROUP_ID|--group GROUP_ID|--members A,B|--go GO:NNNNNNN|--neighbors GENE --flank-genes N|--random-size N --seed N] [--genome GENOME_ID] [--output OUTPUT.json]`
    - `gene-sets produce direct-list --cache CACHE.json_or_tsv [--query LIST_ID] [--genome GENOME_ID] [--provider-id ID] [--provider-version VERSION] [--cache-version VERSION] [--organism NAME|--taxon-id N|--namespace NAMESPACE] [--output OUTPUT.json]`
    - `gene-sets produce ontology-assignment --cache CACHE.json_or_tsv --term GO:NNNNNNN [--ontology-namespace GO] [--evidence-code CODE] [--genome GENOME_ID] [--provider-id ID] [--provider-version VERSION] [--cache-version VERSION] [--organism NAME|--taxon-id N|--namespace NAMESPACE] [--output OUTPUT.json]`
    - `gene-sets produce co-regulated --cache CACHE.json_or_tsv --dataset DATASET_ID --contrast LABEL --score METHOD --threshold RULE --direction both|positive|negative [--relationship co-regulated|anti-co-regulated|manual] [--genome GENOME_ID] [--output OUTPUT.json]`
    - `gene-sets create-pool RESOLUTION_ID --member-container MEMBER_ID=CONTAINER_ID [--member-container MEMBER_ID=CONTAINER_ID ...] [--output-prefix PREFIX] [--container-name NAME] [--preview | --apply --expected-plan-fingerprint-sha256 SHA256] [--path REPORT.json]`
    - `gene-sets promoter-cohort GENOME_ID [--resolution RESOLUTION.json|--group GROUP_ID|--members A,B|--go GO:NNNNNNN] [--relationship manual|co-regulated|anti-co-regulated] [--output OUTPUT.json]`
    - `gene-sets regulatory-partner-screen GENOME_ID --resolution RESOLUTION.json --anchor-motif MOTIF --partner-motif MOTIF [--cutrun-support SUPPORT.json] [--max-distance-bp N] [--output OUTPUT.json]`

## Catalog Curation and Resolved Gene Sets

- `gene-groups list [--catalog PATH] [--filter TEXT] [--output OUTPUT.json]`
  - Lists deterministic local gene-group catalog entries.
  - Default discovery uses built-in, system, user, and project overlays.
  - The built-in catalog declares Gene Ontology as an external resource
    namespace (`GO`) and seeds groups such as `yamanaka_factors`.
  - The built-in catalog also seeds `regulation_of_alternative_splicing`,
    anchored to `GO:0000381` while remaining local reviewable GENtle knowledge.
  - `resources status` reports `gene_ontology` as a declared external mapping
    namespace; GO download/indexing is not implemented in this slice.
- `gene-groups show GROUP_ID [--catalog PATH] [--output OUTPUT.json]`
  - Shows one group with aliases, members, curation status, external mappings,
    and source provenance.
  - Release/runbook helpers such as
    `scripts/pancreas_gene_rna_screen.sh group-plan GROUP` consume this same
    shared output to write reviewable per-member RNA-screen command plans
    without hard-coding TP73, Gene Ontology, or any other single resource as the
    authoritative group definition.
- `gene-groups resolve TOKEN [--catalog PATH] [--output OUTPUT.json]`
  - Resolves a user token against group ids, labels, and aliases.
- `gene-groups doctor [--catalog PATH] [--output OUTPUT.json]`
  - Validates catalog overlays and reports duplicate ids, alias collisions,
    malformed memberships, malformed GO ids, mapping/resource issues, and
    source digests.
- `gene-groups draft --description TEXT [--id ID] [--label LABEL] [--short-description TEXT] [--organism NAME] [--taxon-id N] [--namespace NAMESPACE] [--alias TEXT] [--tag TEXT] [--usage TEXT] [--member SYMBOL] [--members A,B,C] [--candidate SYMBOL=EVIDENCE] [--unresolved-candidate TEXT] [--go GO:NNNNNNN] [--agent-provider NAME] [--agent-model NAME] [--agent-generated-at UTC] [--provenance TEXT] [--output GROUP.json]`
  - Creates a review-gated catalog fragment from an explicit long description
    and optional candidate members/mappings.
  - `--member` marks user-supplied candidates; `--candidate SYMBOL=EVIDENCE`
    marks agent/literature-suggested candidates and preserves the evidence
    note in the member row.
  - Output is `gentle.gene_group_draft.v1` with `review_required=true`, the
    description hash, agent provenance when supplied, member counts,
    candidate-member rows, unresolved candidates, the proposed group record,
    warnings, and the output path when a fragment was written.
  - The command never promotes the group into a trusted catalog by itself; add
    the generated fragment to a user/project catalog only after review.
- `gene-sets resolve [GROUP_ID|--group GROUP_ID|--members A,B|--go GO:NNNNNNN|--neighbors GENE --flank-genes N|--random-size N --seed N] [--genome GENOME_ID] [--catalog PATH] [--genome-catalog PATH] [--allow-draft] [--allow-deprecated] [--output OUTPUT.json]`
  - Resolves an analysis-ready gene-set operand from catalog groups, explicit
    members, local external mappings, same-chromosome neighbors, or
    deterministic random sampling.
  - The `Resolve new set` mode in `Genome > Gene Set Inspector...` constructs
    this same engine operation for all five source kinds, persists the report,
    and exposes resolved/unresolved rows plus provenance and JSON export.
  - `gene-groups resolve` identifies a catalog entry; `gene-sets resolve`
    expands genes, applies curation/member gating, deduplicates resolved
    identities, and records provenance/warnings.
- `gene-sets create-pool RESOLUTION_ID --member-container MEMBER_ID=CONTAINER_ID [--member-container MEMBER_ID=CONTAINER_ID ...] [--output-prefix PREFIX] [--container-name NAME] [--preview | --apply --expected-plan-fingerprint-sha256 SHA256] [--path REPORT.json]`
  - Converts explicit user intent, not logical membership alone, into one
    physical `Pool` container. The resolution must be complete and every member
    must bind to a distinct, exclusive singleton source container containing
    exactly one loaded sequence.
  - Preview is the default and does not mutate the project. It returns a
    `gentle.gene_set_pool_creation.v1` report with the membership lock, exact
    source snapshot hashes, planned derived aliquot ids, reserved container id,
    and `sha256_gene_set_pool_plan_v1` fingerprint.
  - Apply requires the exact preview fingerprint. It derives one lineage-linked
    aliquot sequence from each source, places the aliquots in the new physical
    pool, and retains the source containers and their latest mappings. Any
    unresolved member, missing/duplicate binding, stale source snapshot, or
    namespace drift rejects before mutation.
  - `--path` optionally writes the portable report JSON; it does not export a
    `gentle.pool.v1` artifact. Use `collections run export-pool` on the created
    container for that later handoff.
- `gene-sets produce direct-list --cache CACHE.json_or_tsv [--query LIST_ID] [--genome GENOME_ID] [--provider-id ID] [--provider-version VERSION] [--cache-version VERSION] [--organism NAME|--taxon-id N|--namespace NAMESPACE] [--filter FIELD=VALUE] [--output OUTPUT.json]`
  - Reads a local direct gene-list cache and emits a
    `gentle.gene_set_resolution.v1` report with
    `producer.producer_kind="direct_gene_list"`.
  - JSON caches can use schema `gentle.gene_set_direct_list_cache.v1` with
    top-level `members[]` or `lists[]` entries selected by `--query`.
    TSV/CSV caches use `symbol` or `gene_id` columns and may include
    `list_id`/`list_label` for query selection.
  - Candidate members are resolved through the existing `explicit_members`
    resolver, so genome-coordinate resolution, duplicate collapse, and
    unresolved-member reporting match `gene-sets resolve --members`.
  - The route is offline-only. It requires provider provenance plus
    provider/cache version metadata and at least one interpretation context
    field (`organism`, `taxon_id`, or `symbol_namespace`) from the cache or
    CLI flags.
- `gene-sets produce ontology-assignment --cache CACHE.json_or_tsv --term GO:NNNNNNN [--ontology-namespace GO] [--evidence-code CODE] [--genome GENOME_ID] [--provider-id ID] [--provider-version VERSION] [--cache-version VERSION] [--organism NAME|--taxon-id N|--namespace NAMESPACE] [--filter FIELD=VALUE] [--output OUTPUT.json]`
  - Reads a local ontology-assignment cache and emits a
    `gentle.gene_set_resolution.v1` report with
    `producer.producer_kind="ontology_assignment"`.
  - JSON caches can use schema
    `gentle.gene_set_ontology_assignment_cache.v1` with assignment arrays
    named `assignments`, `terms`, `ontology_terms`, `entries`, or `rows`.
    TSV/CSV caches use term columns such as `term`, `term_id`, or `go_id` plus
    `symbol` or `gene_id`.
  - `--evidence-code CODE` is recorded as a structured
    `evidence_code=CODE` filter and filters assignment rows before member
    resolution. Additional `--filter FIELD=VALUE` constraints are also stored
    and applied to assignment rows.
  - This route is distinct from `gene-sets resolve --go GO:NNNNNNN`, which
    resolves local GENtle catalog groups that cite the GO term. Ontology
    assignment retrieval resolves provider/cache membership rows through
    `explicit_members`.
  - The route is offline-only. A term with no matching assignment rows returns
    an unresolved term row plus a warning rather than silently returning an
    empty success.
- `gene-sets produce co-regulated --cache CACHE.json_or_tsv --dataset DATASET_ID --contrast LABEL --score METHOD --threshold RULE --direction both|positive|negative [--relationship co-regulated|anti-co-regulated|manual] [--genome GENOME_ID] [--provider-id ID] [--provider-version VERSION] [--cache-version VERSION] [--organism NAME|--taxon-id N|--namespace NAMESPACE] [--filter FIELD=VALUE] [--output OUTPUT.json]`
  - Reads a local expression/regulatory cohort cache and emits a
    `gentle.gene_set_resolution.v1` report with
    `producer.producer_kind="co_regulated_cohort"` and
    `co_regulated_metadata`.
  - JSON caches can use schema `gentle.gene_set_co_regulated_cache.v1` with
    row arrays named `rows`, `records`, `members`, or `entries`. TSV/CSV
    caches use `dataset_id`, `contrast`, `symbol`/`gene_id`, and a numeric
    score column such as `score`, `value`, `logfc`, or the method named by
    `--score`.
  - `--threshold` accepts `score>=N`, `score>N`, `score<=N`, `score<N`,
    `abs>=N`, or `abs>N`. `--direction` accepts `both`, `positive`, or
    `negative`.
  - The default relationship is `co-regulated` for this route, but callers can
    record `manual` or `anti-co-regulated`. The relationship is preserved as an
    expectation in metadata; it is not a computed regulatory verdict.
  - The route is offline-only and does not run expression quantification or
    regulatory inference. It materializes candidate members from an existing
    local cache.
  - Produced gene-set reports are persisted as logical lineage artifacts.
    `render-lineage-svg` renders them as `GeneSet` nodes linked to their
    producer operation and downstream promoter/CUT&RUN analysis reports, not as
    sequence or pool nodes.
- `gene-sets promoter-cohort GENOME_ID [--resolution RESOLUTION.json|--group GROUP_ID|--members A,B|--go GO:NNNNNNN|--neighbors GENE --flank-genes N|--random-size N --seed N] [--relationship manual|co-regulated|anti-co-regulated] [--upstream-bp N] [--downstream-bp N] [--catalog PATH] [--genome-catalog PATH] [--output OUTPUT.json]`
  - `--relationship` records a user-declared expectation for downstream
    evidence inspection; it does not claim the genes are co-regulated or
    anti-co-regulated.
  - Builds promoter windows for the resolved set using the shared promoter/TSS
    resolver and the default 1000 upstream / 200 downstream geometry.
  - `--relationship` records a declared cohort expectation. CUT&RUN gene-set
    regulatory-support reports can derive non-blocking occupancy consistency
    flags from that declaration; promoter cohorts simply carry the declaration
    forward.
  - V1 stays offline-first: GO resolves only from local `external_mappings`,
    and live Ensembl ortholog/paralog retrieval is not implemented.
  - Draft members are included with explicit warnings; recipes that used
    `.status // "included"` may miss those warnings.
- `gene-sets regulatory-partner-screen GENOME_ID --resolution RESOLUTION.json --anchor-motif MOTIF [--anchor-motif MOTIF ...] --partner-motif MOTIF [--partner-motif MOTIF ...] [--min-llr-bits BITS] [--motif-threshold MOTIF=BITS] [--max-distance-bp N] [--anchor-mode occupancy-preferred|motif-only] [--cutrun-support SUPPORT.json] [--upstream-bp N] [--downstream-bp N] [--genome-catalog PATH] [--cache-dir PATH] [--output OUTPUT.json]`
  - Runs `InspectRegulatoryPartnerScreen` over the exact members in an existing
    `gentle.gene_set_resolution.v1` report. It resolves strand-aware promoter
    sequences through the prepared-genome index and uses the shared uncapped
    TFBS scanner for both motif roles.
  - The default minimum LLR is `0` bits and the default maximum center distance
    is `100` bp. Repeat `--motif-threshold MOTIF=BITS` for stricter per-motif
    absolute thresholds. Anchor and partner queries must resolve to distinct
    motif ids.
  - The portable result contains exact promoter/genome motif intervals, every
    anchor-by-partner tuple, signed distance and orientation, source-report
    SHA-256 bindings, and one replayable decision-tree trace per gene.
  - `--anchor-mode occupancy-preferred` uses an optional existing
    `gentle.gene_set_cutrun_regulatory_support.v1` report. An evaluated promoter
    with zero support is a fail branch; an unevaluated promoter is an unknown
    branch. `motif-only` leaves occupancy unknown by design.
  - The Gene Set Inspector consumes the same command and opens with a compact
    candidate table; its TSV clipboard view projects the same ledger. The
    advanced tree still highlights genes by trace, and selecting a gene shows
    promoter-oriented DNA with exact anchor, partner, overlap, and TSS spans.
    The GUI's explicit `Use 150 bp` shortcut corresponds to `--max-distance-bp
    150`; no shared defaults change. This is a TSS-window motif screen, not a
    CUT&RUN-summit-centred scan, and expression is not evaluated.
  - Motif proximity and promoter occupancy are association evidence, not proof
    of direct or causal co-regulation. This first screen does not perform
    matched-control enrichment, multiple-testing correction, time-ordering, or
    distal-element analysis.

<!-- manual-page-nav:start -->

---

[CLI Manual](../cli.md) | [Previous: Host Profiles](hosts.md) | [Next: Source-Coherent Transcripts](transcripts.md)
<!-- manual-page-nav:end -->
