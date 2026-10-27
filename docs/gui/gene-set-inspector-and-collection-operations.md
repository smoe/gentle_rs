### Gene Set Inspector And Collection Operations

<!-- manual-counterparts:start -->
Related CLI reference: [Collections](../cli/collections.md) | [Gene Groups and Gene Sets](../cli/gene-sets.md).
These topic links do not imply identical feature coverage across interfaces.
<!-- manual-counterparts:end -->

- Gene sets and other multi-sequence selections should enter the GUI through
  the same collection-subject mental model rather than as a separate editor
  island:
  - a gene set is a logical collection source that resolves to auditable
    members, warnings, and provenance
  - a pool/container is a physical collection that may contain multiple
    molecules in one sample
  - an arrangement is an ordered collection for lanes, racks, plates, and
    labels
  - an alignment is an analytical collection with member order plus column
    correspondence
  - derived collections, such as promoter windows or neighboring loci, preserve
    links back to each source member
- Collection GUI affordances must stay thin over shared engine/shell contracts:
  the GUI should collect the collection operand, show per-member and aggregate
  readiness/errors/results, and call the same named operation or shell route
  that CLI/MCP/agent surfaces can invoke.
- `Genome > Gene Set Inspector...` has three shared-engine modes:
  - `Resolve new set` authors catalog-group, explicit-member, local external
    mapping, genomic-neighbor, or deterministic-random sources
  - source-specific validation and shared genome/catalog/cache gates construct
    the same `ResolveGeneSet` payload as `gene-sets resolve`
  - resolution runs in a detached engine snapshot, persists through the normal
    engine operation, and auto-selects the exact new report when its commit is
    current; stale results do not change project state
  - `Inspect & run` shows report-level review state, resolved and unresolved
    members, member-owned status/confidence metadata, warnings, provenance, and
    portable JSON copy/export
  - `Regulatory partners` accepts anchor and partner motif queries, absolute
    thresholds, promoter geometry, a distance rule, and an optional persisted
    gene-set CUT&RUN support report. It dispatches the shared
    `gene-sets regulatory-partner-screen` command on a detached engine snapshot
    and caches the portable result in the inspector
  - results open with a **Candidate list**: gene, predicted anchor and partner
    matrix IDs/sites, signed motif-centre distance, both LLR scores, proximity
    outcome, promoter-level CUT&RUN support, and expression (`Not evaluated`).
    Nearby pairs are shown first; a checkbox includes distant pairs without
    changing the report. Counts are motif pairs/member rows, not independent
    peaks or enrichment. Genes without displayed pairs remain inspectable.
  - `Use 150 bp` explicitly changes the form's proximity parameter; it does not
    change the existing shared default or select a 500 bp neighbourhood.
    The current scan window remains **TSS-centred**, not CUT&RUN-peak-centred.
    ENCODE/Ensembl regulatory annotation and promoter similarity are not gates.
  - `Copy candidate table (TSV)` copies the displayed rows with tuple/member,
    transcript, genome, operation/run identifiers and input-report digests.
    Keep `Copy report JSON` for the complete evidence and threshold policy.
  - selecting a candidate's gene opens its promoter sequence detail. The
    **Advanced evidence: decision tree and all genes** section retains the
    branching tree derived from the engine report.
    Hovering a node highlights every gene whose recorded trace traverses that
    node; selecting a gene displays its strand-aware promoter DNA with the
    exact anchor, partner, overlapping-role, and TSS bases highlighted
  - the tree and sequence renderer do not rescan DNA or derive GUI-local
    biology. Exact motif intervals, tuples, occupancy state, and node membership
    all come from `gentle.regulatory_partner_screen.v1`
  - wording remains evidence-focused: motif-family matches and promoter-level
    occupancy are cross-gene association evidence, not proof of a unique bound
    protein or causal co-regulation
  - the next peak-centred/expression-aware slice is specified separately in
    [`cofactor_candidate_interface_plan.md`](../cofactor_candidate_interface_plan.md);
    the current view does not substitute predicted motifs for measured summits
  - external mappings are resolved from configured local catalogs; opening the
    form does not perform a live ontology lookup
  - the collection-operation catalog is projected from canonical
    engine-operation lifting policies rather than maintained as a GUI-only list
  - typed adapters cover promoter-cohort derivation (`derive`),
    primer-specificity assessment (`map`), restriction-site scanning (`map`),
    non-mutating TFBS hit scanning (`map`), restriction digestion (`map`), and
    explicit physical gene-set pooling (`combine`), each with independent form
    state and operation-specific results
  - pool export, pooled-gel rendering, and serial arrangement remain visible but
    disabled with their typed physical-pool or materialization requirements
  - primer specificity binds every logical member to one exact persisted
    primer-design report, then chooses a one-based pair rank and prepared local
    BLAST target
  - restriction scanning binds every logical member to one loaded DNA
    sequence, accepts optional enzyme names/site caps/cut geometry, and shows
    aggregate site counts without guessing sequence identity from gene labels
  - TFBS hit scanning binds the same explicit sequences, requires one or more
    motif queries, and exposes shared default/per-TF likelihood thresholds.
    Its computation is unlimited by default; result rows show retained hit
    counts and whether caps, failures, or unscanned effective motifs made the
    aggregate incomplete
  - restriction digestion binds every logical member to one loaded DNA
    sequence and requires explicit known enzyme names. The first run previews
    deterministic fragment ids without changing the project; apply is enabled
    only for that exact fingerprint-locked plan. Applied fragments retain
    per-source lineage and can be opened directly from the result table. The
    operation does not create a pooled container
  - physical pool creation requires a complete resolution and binds every
    member to one distinct exclusive singleton source container. Preview locks
    the exact source snapshots, output aliquot ids, and reserved pool id;
    apply creates lineage-linked derived aliquots in one exclusive physical
    pool while retaining all source containers
  - promoter derivation chooses one prepared genome, strand-aware window sizes,
    and an optional relationship expectation
  - all adapters execute one shared shell command against one complete
    collection membership in a detached engine snapshot; `Discard result`
    always suppresses display and cooperatively prevents commit when the worker
    observes the request before detached commit, but does not remove an output
    file the shared command already wrote
  - result views keep per-member execution outcomes separate from biological
    verdicts and verify that the returned collection subject is still the one
    displayed by the inspector
- The inspector never infers an assay from a gene symbol. Missing bindings,
  unavailable pair ranks, mixed or mismatched biological context, stale
  detached results, and member-level execution failures remain explicit. The
  portable collection report can be copied as JSON, and persisted child reports
  can be opened in PCR Designer.
- Detailed implementation plan:
  [`gui_gene_set_collection_operations_plan.md`](../gui_gene_set_collection_operations_plan.md).
- Operations offered for a collection should expose the engine-declared lifting
  behavior in their controls and result view:
  - map per member
  - combine into one explicit pool/sample
  - compare members together
  - arrange/place members
  - derive new members from each source member
  - reject with a typed reason when a collection is not valid input

<!-- manual-page-nav:start -->

---

[GUI Manual](../gui.md) | [Previous: Anchored Region Extraction (Engine Ops)](anchored-region-extraction-engine-ops.md) | [Next: Splicing Expert Workspace](splicing-expert-workspace.md)
<!-- manual-page-nav:end -->
