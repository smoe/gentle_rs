# Prepare a GeneArt Handoff from Shared External-Service Contracts

> Type: `GUI + CLI walkthrough`
> Status: `manual/hybrid`
> Drift note: this page is hand-written, but every operational step routes
> through shared `services ...` shell/CLI contracts and the provider catalog.

This tutorial rehearses GENtle's first GeneArt integration surface without
performing any vendor-side action. It is useful when you want to check whether
a construct or protein-service idea can be represented as a vendor-neutral
external-service request, preflighted locally, and exported as a reviewable
handoff bundle.

## What GeneArt Is

Invitrogen GeneArt Services are Thermo Fisher Scientific's custom DNA and
protein-expression service family. In practical GENtle terms, GeneArt is a
provider that can receive reviewed sequence designs and project context for
work such as synthetic DNA fragments, cloned gene synthesis, plasmid-related
services, mutagenesis/library projects, and protein expression or purification.

The vendor-facing workflow is normally commercial and account-bound: users
design or upload project material, review optimization and provider checks,
request quotes, and track orders through Thermo Fisher/GeneArt channels such as
the GeneArt Services Dashboard. GENtle's role in this tutorial is deliberately
smaller and safer: it prepares structured, reviewable handoff artifacts that
help a person or external orchestrator decide what should be sent to the vendor
later.

The important idea is the boundary:

- GENtle owns a provider-neutral request, preflight, quote, and artifact-bundle
  contract.
- GeneArt is one configured provider in that contract.
- ClawBio, MCP, the CLI, and the GUI can all ask for the same preflight/quote
  result.
- The user or external orchestrator must decide what should be sent back, for
  example a GenBank payload, FASTA, amino-acid sequence, quote metadata, or the
  complete bundle.

## Safety Boundary

GENtle does not do any of these in this tutorial:

- submit a GeneArt cart or order,
- upload a project through the GeneArt API,
- scrape the Thermo Fisher/GeneArt dashboard,
- check or store API credentials,
- store PO, account, shipping, billing, or procurement secrets in project state,
- claim that synthetic tutorial sequences are ready to order.

The example requests are synthetic and named `DEMO_DO_NOT_ORDER` on purpose.
Replace them with reviewed project material only after biology, vector context,
biosafety, provider representation, and local purchasing policy have all been
checked.

## Inputs

Use the bundled example requests:

- [`docs/examples/external_services/geneart_cloned_gene_request.json`](../examples/external_services/geneart_cloned_gene_request.json)
- [`docs/examples/external_services/geneart_protein_expression_request.json`](../examples/external_services/geneart_protein_expression_request.json)

The built-in provider config is:

- [`assets/external_service_providers.json`](../../assets/external_service_providers.json)

Optional project overlays can live under:

- `.gentle/catalogs/external_service_providers.json`
- `.gentle/catalogs/external_service_providers.d/*.json`

Use overlays for local policy, wording, or provider-channel updates. Do not use
them to hide the safety boundary; direct ordering remains disabled until a
future explicitly confirmed implementation adds it.

## Step 1: Check Provider Config Health

Run the doctor first. It validates the active provider config chain and reports
source provenance.

```bash
cargo run --quiet --bin gentle_cli -- services providers doctor
```

Expected outcome:

- schema is `gentle.external_service_provider_config_doctor.v1`,
- `error_count` is `0`,
- provider count includes `geneart`,
- source rows show which built-in/system/user/project config files were used.

For a pinned built-in check:

```bash
cargo run --quiet --bin gentle_cli -- services providers doctor \
  --catalog assets/external_service_providers.json
```

## Step 2: List Provider Capabilities

```bash
cargo run --quiet --bin gentle_cli -- services providers list
```

Expected GeneArt rows:

- `dna_fragment`
- `cloned_gene`
- `plasmid_reorder`
- `mutagenesis`
- `protein_expression`

Important interpretation:

- DNA fragment, cloned-gene, and plasmid-reorder rows may record documented
  direct API surfaces as planned.
- `direct_api_implemented` should remain `false` in this tutorial.
- Protein expression and mutagenesis are quote/handoff first.
- The capability catalog is the same record the GUI and automation routes
  consume.

## Step 3: Preflight the Cloned-Gene Example

```bash
cargo run --quiet --bin gentle_cli -- services project-preflight \
  @docs/examples/external_services/geneart_cloned_gene_request.json
```

Expected outcome:

- schema is `gentle.external_service_preflight.v1`,
- `provider` is `geneart`,
- `service_kind` is `cloned_gene`,
- `eligible` is `true`,
- `quote_handoff_available` is `true`,
- `direct_submission_available` is `false`,
- warnings explain that direct API support may be documented externally but is
  not implemented by GENtle yet.

If this fails, read `blocking_issues[]` first. For real cloned-gene requests,
common missing pieces are an actual DNA sequence/`seq_id`, reviewed vector
context, or a local policy overlay that requires more fields.

## Step 4: Prepare the Cloned-Gene Quote Handoff

```bash
cargo run --quiet --bin gentle_cli -- services project-quote \
  @docs/examples/external_services/geneart_cloned_gene_request.json
```

Expected quote output:

- schema is `gentle.external_service_quote.v1`,
- `quote_status` is `handoff_ready`,
- `quote_mode` is a handoff mode, not a submission mode,
- `service_ready_bundle.schema` is
  `gentle.external_service_artifact_bundle.v1`,
- inline payloads include at least:
  - `handoff_markdown`,
  - `redacted_request_json`,
  - `normalized_line_items_json`,
  - `normalized_line_items_csv`,
  - `guided_wop_checklist` or provider-channel checklist when configured.

Treat warnings as useful handoff context. They are there to keep the human
review boundary visible.

## Step 5: Export the Cloned-Gene Handoff Bundle

When the preview looks sensible, write the same handoff payloads as files. The
example path is under ignored `artifacts/` so a tutorial run does not dirty the
repository.

```bash
cargo run --quiet --bin gentle_cli -- services project-quote \
  @docs/examples/external_services/geneart_cloned_gene_request.json \
  --output-dir artifacts/external_services/geneart_cloned_gene_demo
```

Expected files:

- `quote_report.json`
- `01_handoff_markdown.md`
- `02_redacted_request_json.json`
- one or more normalized line-item/checklist payloads, depending on the active
  provider config

The returned quote report should also list those files in
`service_ready_bundle.local_files[]`.

## Step 6: Preflight and Quote the Protein-Expression Example

Protein service handoff is intentionally not treated as a direct GeneArt API
order in GENtle v1. The tutorial still uses the same request/preflight/quote
contract, with `return_spec` stating what an external caller may want back.

```bash
cargo run --quiet --bin gentle_cli -- services project-preflight \
  @docs/examples/external_services/geneart_protein_expression_request.json

cargo run --quiet --bin gentle_cli -- services project-quote \
  @docs/examples/external_services/geneart_protein_expression_request.json \
  --output-dir artifacts/external_services/geneart_protein_expression_demo
```

Expected difference from the cloned-gene example:

- `service_kind` is `protein_expression`,
- direct submission remains unavailable,
- warnings should explain that protein expression is quote/handoff first,
- `return_spec.requested_payloads[]` asks for `amino_acid_sequence`,
  `codon_targeted_dna`, `quote_metadata`, and `handoff_bundle`.

This is the place where ClawBio should be explicit. If it only needs the amino
acid sequence, it should request that. If it needs an adjusted GenBank entry or
full vendor review packet, it should request those payload families instead of
assuming one fixed return shape.

## Step 7: Repeat the Same Review in the GUI

Open GENtle and use the shared GUI inspector:

```bash
cargo run --bin gentle
```

Then:

1. Open `Services -> External Services...`.
2. Press `Refresh Providers`.
3. Confirm GeneArt appears in the provider picker.
4. Select `geneart` and `cloned_gene`.
5. Press `Use Selected Template`, or paste the cloned-gene request JSON into
   the request editor.
6. Press `Preflight`.
7. Press `Prepare Quote Handoff`.
8. Inspect payload previews, warnings, and required follow-up.
9. Set `Output dir` to
   `artifacts/external_services/geneart_cloned_gene_demo`.
10. Press `Export Handoff Bundle`.
11. Repeat with `protein_expression` and the protein request JSON.

The GUI should not contain GeneArt-specific business logic. If GUI and CLI
disagree, treat that as a bug in shared-shell plumbing or presentation, not as
a vendor-specific GUI feature to patch separately.

## Step 8: What a Real Project Should Replace

Before any real vendor interaction, replace tutorial placeholders with reviewed
project material:

- `source_target.sequence` or `seq_id` from an approved construct,
- `vector_spec` from a reviewed helper profile or vendor-onboarded vector,
- expression host, tag, purification, and QC decisions for protein projects,
- local biosafety and procurement review notes,
- a `return_spec` that says what the caller actually wants back,
- commercial/account/shipping references held outside GENtle project state.

Do not add PO numbers, account credentials, shipping addresses, or billing
details into request JSON. Use `commercial_context_ref` only as a non-secret
pointer to a separately managed process if your local workflow needs one.

## Step 9: What to Send to GeneArt Later

After the GENtle representation feels internally sane, a vendor-review packet
can be prepared from generated artifacts. Useful review material is:

- the relevant provider config row,
- one `services providers doctor` report,
- one cloned-gene `project-quote` output,
- one protein-expression `project-quote` output if relevant,
- the GenBank/FASTA/protein payloads requested by `return_spec`,
- a short explanation that GENtle v1 prepares handoff artifacts only.

Good questions for vendor or procurement review:

- Are service-kind names and product mappings understandable?
- Which fields should be first-class for cloned genes, fragments, plasmid
  reorder, mutagenesis, and protein expression?
- Which status/QAD artifacts should be imported later when direct API support
  is implemented?
- What should ClawBio request back for a given workflow: amino-acid sequence,
  adjusted GenBank, quote metadata, vendor ids, or the full bundle?
- Are any warning or follow-up statements misleading for your account route?

## Success Checklist

Mark this tutorial successful when all of these are true:

- provider doctor passes with no schema errors,
- provider list shows GeneArt rows for DNA and protein-service tracks,
- cloned-gene example preflights as eligible,
- protein-expression example preflights as eligible,
- both quote calls return `handoff_ready`,
- exported bundle directories contain `quote_report.json`,
- no vendor order, cart submission, API upload, credential lookup, or dashboard
  automation occurred,
- GUI and CLI reports agree on provider, service kind, eligibility, and direct
  submission availability.

## Follow-Up Work

Good next implementation slices are:

- derive GeneArt request JSON directly from a selected construct/helper profile,
- map `protein_to_dna_handoff` output into the protein-expression request,
- add richer GenBank/FASTA/protein payload materialization from `return_spec`,
- add manual/imported provider status records,
- add future direct GeneArt API submission only behind explicit account
  enablement, credential handling, and external-order confirmation gates.

## Feedback

If this tutorial is confusing, execution-stale, biologically suspect, or missing a useful figure, please open the matching tutorial issue template and include the context copied from GENtle Help -> Tutorial -> Copy Feedback Context.

- Tutorial title:
- Tutorial/chapter id:
- Step reached:
- Expected vs. actual:
- Interface used: GUI / CLI / Agent Assistant / ClawBio

Paste the Tutorial feedback context here:

```text

```
