# Prepare a Metabion Handoff from Shared External-Service Contracts

> Type: `GUI + CLI walkthrough`
> Status: `manual/hybrid`
> Drift note: this page is hand-written, but every operational step routes
> through shared `services ...` shell/CLI contracts and the provider catalog.

This tutorial is the first safe end-to-end rehearsal for GENtle's
external-service handoff layer. It uses Metabion as the concrete example, but
the implementation remains provider-neutral: the GUI, CLI, agents, and ClawBio
all consume the same provider catalog, request schema, preflight report, and
quote-handoff bundle.

The point is not to submit an order. The point is to verify that GENtle can
prepare clear, inspectable, vendor-reviewable artifacts:

- provider capability discovery,
- provider-config doctor output,
- local preflight,
- normalized line-item JSON/CSV,
- email-draft markdown,
- guided WOP checklist, and
- explicit follow-up warnings for human review.

## Safety Boundary

GENtle does not do any of these in this tutorial:

- scrape Metabion WOP,
- submit a cart or order,
- look up or store credentials,
- store PO, account, shipping, or billing data in project state,
- assert that a tutorial sequence is ready to order.

The example request files are synthetic and named `DEMO_DO_NOT_ORDER` on
purpose. Replace them with reviewed project material only after the biology,
provider mapping, and local purchasing policy have been checked.

## Inputs

Use the bundled example requests:

- [`docs/examples/external_services/metabion_oligo_single_tube_request.json`](../examples/external_services/metabion_oligo_single_tube_request.json)
- [`docs/examples/external_services/metabion_mblock_request.json`](../examples/external_services/metabion_mblock_request.json)

The built-in provider config is:

- [`assets/external_service_providers.json`](../../assets/external_service_providers.json)

Optional project overlays can live under:

- `.gentle/catalogs/external_service_providers.json`
- `.gentle/catalogs/external_service_providers.d/*.json`

See the example overlay:

- [`docs/examples/catalogs/external_service_providers_project_overlay.json`](../examples/catalogs/external_service_providers_project_overlay.json)

## Step 1: Check Provider Config Health

Run the doctor first. It validates the active provider config chain and reports
source provenance.

```bash
cargo run --quiet --bin gentle_cli -- services providers doctor
```

Expected outcome:

- schema is `gentle.external_service_provider_config_doctor.v1`,
- `error_count` is `0`,
- provider count includes `metabion`,
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

Expected Metabion rows:

- `dna_oligo_single_tube`
- `dna_fragment`, mapped to m-block DNA fragments/libraries
- submission modes are handoff modes such as `wop_handoff` and
  `email_excel_handoff`
- direct API/order submission is not implemented

This is the same catalog the GUI provider picker consumes.

## Step 3: Preflight the Oligo Example

```bash
cargo run --quiet --bin gentle_cli -- services project-preflight \
  @docs/examples/external_services/metabion_oligo_single_tube_request.json
```

Expected outcome:

- schema is `gentle.external_service_preflight.v1`,
- `provider` is `metabion`,
- `service_kind` is `dna_oligo_single_tube`,
- `eligible` is `true`,
- `quote_handoff_available` is `true`,
- `direct_submission_available` is `false`.

If this fails, read `blocking_issues[]` first. For real requests, the most
common failure should be a missing source field such as `source_target.sequence`
or `source_target.line_items`.

## Step 4: Prepare the Oligo Quote Handoff

```bash
cargo run --quiet --bin gentle_cli -- services project-quote \
  @docs/examples/external_services/metabion_oligo_single_tube_request.json
```

Expected quote output:

- schema is `gentle.external_service_quote.v1`,
- `quote_status` is `handoff_ready`,
- `service_ready_bundle.schema` is
  `gentle.external_service_artifact_bundle.v1`,
- inline payloads include:
  - `redacted_request_json`,
  - `normalized_line_items_json`,
  - `normalized_line_items_csv`,
  - `email_draft_markdown`,
  - `guided_wop_checklist`.

Treat warnings as useful handoff context, not noise. For example, missing local
vendor Excel templates should be warnings because the WOP/email route is still
usable with explicit human review.

## Step 5: Preflight and Quote the m-block Example

```bash
cargo run --quiet --bin gentle_cli -- services project-preflight \
  @docs/examples/external_services/metabion_mblock_request.json

cargo run --quiet --bin gentle_cli -- services project-quote \
  @docs/examples/external_services/metabion_mblock_request.json
```

Expected difference from the oligo example:

- the product name in normalized line items should refer to m-block DNA
  fragments/libraries,
- warnings or required follow-up should mention biosafety review where
  applicable,
- direct submission remains unavailable.

## Step 6: Repeat the Same Review in the GUI

Open GENtle and use the shared GUI inspector:

```bash
cargo run --bin gentle
```

Then:

1. Open `Services -> External Services...`.
2. Press `Refresh Providers`.
3. Confirm Metabion appears in the provider picker.
4. Select `metabion` and `dna_oligo_single_tube`.
5. Press `Use Selected Template`, or paste the oligo request JSON into the
   request editor.
6. Press `Preflight`.
7. Press `Prepare Quote Handoff`.
8. Inspect the payload previews and warnings.

The GUI should not contain provider-specific business logic. If the GUI and CLI
disagree, treat that as a bug in shared-shell plumbing or presentation, not as
a vendor-specific GUI feature to patch separately.

## Step 7: Optional Project Policy Overlay

Provider behavior is meant to be locally maintainable without changing the
GENtle source tree. To rehearse that path, copy the example overlay into a
project-local catalog directory:

```bash
mkdir -p .gentle/catalogs/external_service_providers.d
cp docs/examples/catalogs/external_service_providers_project_overlay.json \
  .gentle/catalogs/external_service_providers.d/metabion_local_policy.json
cargo run --quiet --bin gentle_cli -- services providers doctor
```

Expected outcome:

- the doctor reports the project-local overlay source,
- later provider ids override earlier provider ids,
- GUI and CLI provider rows update from the same catalog result.

Remove or adapt the overlay before real use if it is not your lab policy.

## Step 8: What to Send to Metabion Later

After the GENtle representation feels internally sane, a vendor-review email
can be prepared from the generated artifacts. The useful review bundle is:

- the relevant provider config rows,
- one `services providers doctor` report,
- one oligo `project-quote` output,
- one m-block `project-quote` output,
- a short explanation that GENtle v1 prepares handoff artifacts only.

Good questions for vendor review:

- Are the service-kind names and product mappings understandable?
- Are WOP/email/Excel handoff steps represented fairly?
- Are required follow-up policies missing or misleading?
- Which fields should be first-class in future normalized line items?
- Should any product family be excluded until official guidance is clearer?

## Success Checklist

Mark this tutorial successful when all of these are true:

- provider doctor passes with no schema errors,
- provider list shows Metabion DNA oligo and m-block rows,
- both example requests preflight as eligible,
- both quote calls return `handoff_ready`,
- quote output contains normalized line-item JSON/CSV and human-readable draft
  payloads,
- GUI and CLI describe the same provider/service choices,
- no step performs vendor submission or stores commercial secrets.

## Troubleshooting

If the GUI provider picker is empty:

- run `services providers doctor` first and inspect source errors,
- confirm the built-in `assets/external_service_providers.json` exists,
- check project-local overlays for duplicate provider ids with incomplete rows.

If preflight is blocked:

- inspect `blocking_issues[]`,
- make sure `source_target` has either `sequence` or `line_items`,
- make sure `provider` and `service_kind` match the provider catalog exactly.

If quote output is handoff-ready but warns about missing templates:

- this is expected unless a local vendor template fixture was configured,
- use the normalized JSON/CSV and email/WOP checklist as the deterministic
  GENtle-owned handoff layer,
- fetch official vendor templates manually under local purchasing policy.
