# GENtle Protocol-First Examples

This folder defines canonical workflow examples independent of any one adapter
syntax.

## Source of truth

- Canonical examples live in `docs/examples/workflows/*.json`.
- Schema: `gentle.workflow_example.v1`.
- Each file carries:
  - metadata (`id`, `title`, `summary`, `test_mode`, `required_files`)
  - canonical `workflow` JSON payload

`test_mode` values:

- `always`: validated and executed by default test runs
- `online`: executed only when `GENTLE_TEST_ONLINE=1`
- `skip`: parsed and documented, but not executed in automated tests

## Draft design-resource examples

- Draft, non-executable planning/design artifacts may live outside
  `docs/examples/workflows/`.
- Current draft resource examples:
  - `docs/examples/plans/gibson_destination_first_single_insert.json`
  - `docs/examples/plans/gibson_destination_first_multi_insert.json`
  - schema: `gentle.gibson_assembly_plan.v1`
  - `docs/examples/catalogs/helper_semantics_vocabulary_extension.json`
  - schema: `gentle.helper_semantics_vocabulary.v1`
- The single-insert Gibson plan shape is now accepted by the shared
  `gibson preview ...` path when referenced sequence ids exist in the active
  project state.
- Multi-fragment Gibson examples remain design resources and are not yet
  executed by workflow-example tests.
- Helper semantics vocabulary examples show how a project can add local
  helper-construct meaning, aliases, and routine hints under
  `.gentle/catalogs/helper_semantics_vocabulary.d/` without editing GENtle's
  built-in helper catalog.
- Folder note:
  - `docs/examples/plans/README.md`

## Example assets

- `docs/examples/assets/cdna_assay_demo.gb`
  - synthetic/hand-authored 32 bp two-exon GenBank locus for deterministic
    cDNA PCR/qPCR assay-test examples
  - deterministic recreation: use the committed file content directly; it is
    not derived from an upstream biological accession
  - used by:
    `docs/examples/workflows/cdna_pcr_qpcr_assay_test_offline.json` and the
    ClawBio request
    `integrations/clawbio/skills/gentle-cloning/examples/request_workflow_cdna_pcr_qpcr_assay_test_offline.json`
- `docs/examples/assets/cdna_assay_nonspecific_demo.gb`
  - synthetic/hand-authored 72 bp GenBank locus with two transcripts that share
    one primer pair but produce different cDNA amplicon lengths
  - deterministic recreation and usage notes live in
    `docs/examples/assets/README.md`
  - used by:
    `docs/examples/workflows/cdna_pcr_qpcr_product_gel_nonspecific_offline.json`
    and the matching ClawBio product-gel requests

## Generate adapter snippets

Generate Markdown snippets for CLI/shell/JS/Lua:

```bash
cargo run --bin gentle_examples_docs -- generate
```

Generated output is written to `docs/examples/generated`.

Validation only:

```bash
cargo run --bin gentle_examples_docs -- --check
```

## Generate tutorial output

Generated tutorial pages and retained runtime artifacts are committed under:

- `docs/tutorial/generated`

Tutorial source layer:

- `docs/tutorial/sources/catalog_meta.json`
- `docs/tutorial/sources/*.json` (schema `gentle.tutorial_source.v2`)
- generated runtime manifest: `docs/tutorial/manifest.json`
  (schema `gentle.tutorial_manifest.v1`)

Generate tutorial output:

```bash
cargo run --bin gentle_examples_docs -- tutorial-generate
```

Validate committed tutorial output against fresh generation:

```bash
cargo run --bin gentle_examples_docs -- tutorial-check
```

Relevant executable Gibson specialist setup baseline:

- `docs/examples/workflows/gibson_specialist_testing_baseline.json`
  - loads one circular destination and one linear insert under stable IDs for
    `Open Tutorial Project...` and GUI specialist testing
- `docs/examples/workflows/gibson_arrangements_baseline.json`
  - loads the same deterministic destination and insert, then applies the
    canonical single-insert Gibson plan so arrangement/gel tutorials open with
    the assembled product and stored lane arrangement already present

Relevant executable release proof baseline:

- `docs/examples/workflows/tp73_genome_evidence_viewer_release_proof.json`
  - loads the public GRCh38.p14 TP73 locus and overlays tiny local repeat,
    Clariom-style array, CUT&RUN-style BED, and TFBS fixtures for the
    genome-anchored evidence-viewer release path
  - runbook: `docs/tp73_genome_evidence_viewer_runbook.md`

## Test examples

Run default (offline-safe) example tests:

```bash
cargo test workflow_examples -- --test-threads=1
```

Run online examples (explicit opt-in):

```bash
GENTLE_TEST_ONLINE=1 cargo test workflow_examples -- --test-threads=1
```
