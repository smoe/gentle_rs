# Helper Sequences

<!-- manual-counterparts:start -->
Related GUI reference: [Reference Genomes (Main Menu)](../gui/reference-genomes-main-menu.md).
These topic links do not imply identical feature coverage across interfaces.
<!-- manual-counterparts:end -->

Helper convenience commands:

- `helpers list [--catalog PATH]`
  - Same behavior as `genomes list`, but defaults to `assets/helper_genomes.json`.
  - Helper `entries[]` rows now also include an optional normalized
    `interpretation` record so ClawBio/agents/planners can consume one shared
    helper-meaning layer instead of reparsing raw catalog prose.
  - Metadata-only vector candidates may expose `sequence_availability`,
    `redistribution_status`, `biological_safety_note`, and
    `usable_as_empty_backbone` alongside procurement, markers, origins, host
    compatibility, and semantic component records.
- `helpers vocabulary list [--vocabulary PATH] [--filter TEXT]`
  - Lists the resolved helper semantics vocabulary terms that normalize helper
    catalog semantics.
  - Omitting `--vocabulary` uses the built-in/system/user/project vocabulary
    discovery chain.
  - Output includes canonical terms, aliases, descriptions, sources, and
    vocabulary-provided routine hints for debugging deterministic helper
    interpretation.
- `helpers vocabulary doctor [--vocabulary PATH] [--routine-catalog PATH]`
  - Validates helper semantics vocabulary overlays before treating them as
    stable local knowledge.
  - Reports source ordering, parsed-fragment SHA-1 digests, duplicate
    canonical terms, alias collisions, malformed routine hints, and
    missing routine-hint `source_terms` targets, plus routine-family hints that
    are unknown to the selected routine catalog.
  - Unknown routine families are reported as local-extension warnings, so labs
    can still prototype private routine vocabulary without confusing it with
    built-in GENtle support.
- `helpers ensembl-available [--collection all|vertebrates|metazoa] [--filter TEXT]`
  - Same discovery report shape as `genomes ensembl-available`, but exposed under the helper-family command tree for contract symmetry across adapters.
- `helpers install-ensembl SPECIES_DIR [--collection vertebrates|metazoa] [--catalog PATH] [--output-catalog PATH] [--genome-id ID] [--cache-dir PATH] [--timeout-secs N]`
  - Same quick-install contract as `genomes install-ensembl`, but targeting the helper-catalog/cache defaults.
- `helpers validate-catalog [--catalog PATH]`
  - Same behavior as `genomes validate-catalog`, with helper-catalog default.
- `helpers doctor-catalog [--catalog PATH]`
  - Returns `gentle.helper_vector_catalog_doctor.v1` with deterministic
    helper/vector metadata issues such as missing host system, missing sequence
    availability, missing redistribution status, missing biosafety note,
    missing empty-backbone flag, helper phage marked as an empty backbone, and
    vector records without structured semantic components.
- `helpers show-card [--catalog PATH] [--filter TEXT|--name TEXT]`
  - Returns `gentle.helper_vector_card.v1` cards as a pure projection of
    structured catalog fields: identifiers, aliases, helper kind, host system,
    sequence/redistribution/biosafety fields, empty-backbone status,
    metadata-only status, procurement, affordances, constraints, components,
    and relationships.
- `helpers update-ensembl-specs [--catalog PATH] [--output-catalog PATH]`
  - Same behavior as `genomes update-ensembl-specs`, with helper-catalog default.
- `helpers status HELPER_ID [--catalog PATH] [--cache-dir PATH]`
  - Same behavior as `genomes status`, with helper-catalog default
    (including length/mass metadata fields).

<!-- manual-page-nav:start -->

---

[CLI Manual](../cli.md) | [Previous: Reference Genomes](genomes.md) | [Next: Host Profiles](hosts.md)
<!-- manual-page-nav:end -->
