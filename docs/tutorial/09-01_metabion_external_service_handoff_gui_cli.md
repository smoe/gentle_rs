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

This is an offline rehearsal: no network access, email sending, portal opening,
or ordering is needed. `eligible` and `handoff_ready` mean that GENtle's local
request checks passed, not that Metabion accepted a sequence, issued a price,
approved biosafety, or guaranteed synthesis. Provider wording is the bundled
catalog's representation, not a live statement of vendor policy.

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

## Before You Start

Run Bash commands from the repository root. Links above are relative to this
Markdown file; command paths below are relative to the repository root.
Use an already-built CLI from this checkout, or set `GENTLE_CLI` to its absolute
path. Stop if the executable check fails; this walkthrough does not build it.

```bash
GENTLE_CLI="${GENTLE_CLI:-$PWD/target/debug/gentle_cli}"
test -x "$GENTLE_CLI"
RUN_DIR=$(mktemp -d "${TMPDIR:-/tmp}/gentle-metabion.XXXXXX") || exit 1
STATE="$RUN_DIR/tutorial.gentle.json"
printf 'Tutorial files: %s\n' "$RUN_DIR"
```

Keep this shell open. Each run gets new state and output paths; do not delete
an existing state to restart. The export command can overwrite files in an
existing directory, so use the fresh subdirectories below.

For reproducible CLI results, isolate catalog discovery from lab overlays:

```bash
export GENTLE_ASSET_ROOT="$PWD"
export GENTLE_SYSTEM_CONFIG_ROOT="$RUN_DIR/system-config"
export XDG_CONFIG_HOME="$RUN_DIR/user-config"
export GENTLE_PROJECT_ROOT="$RUN_DIR/project"
```

These variables affect this shell and processes launched from it only. A GUI
started elsewhere may still use lab overlays; compare its doctor source rows
before expecting identical results.

## Step 1: Check Provider Config Health

Run the doctor first. It validates the active provider config chain and reports
source provenance.

```bash
"$GENTLE_CLI" --state "$STATE" services providers doctor
```

Expected outcome:

- schema is `gentle.external_service_provider_config_doctor.v1`,
- `error_count` is `0`,
- provider count includes `metabion`,
- source rows show which built-in/system/user/project config files were used.

For a pinned built-in check:

```bash
"$GENTLE_CLI" --state "$STATE" services providers doctor \
  --catalog assets/external_service_providers.json
```

## Step 2: List Provider Capabilities

```bash
"$GENTLE_CLI" --state "$STATE" services providers list
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
"$GENTLE_CLI" --state "$STATE" services project-preflight \
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
"$GENTLE_CLI" --state "$STATE" services project-quote \
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

## Step 5: Export the Oligo Handoff Bundle

When the preview looks sensible, write the same handoff payloads into a fresh
subdirectory of `RUN_DIR`. These files stay outside the repository.

```bash
"$GENTLE_CLI" --state "$STATE" services project-quote \
  @docs/examples/external_services/metabion_oligo_single_tube_request.json \
  --output-dir "$RUN_DIR/metabion_oligo_demo"
```

Expected files:

- `quote_report.json`
- `01_handoff_markdown.md`
- `02_redacted_request_json.json`
- `03_normalized_line_items_json.json`
- `04_normalized_line_items_csv.csv`
- `05_email_draft_markdown.md`
- `06_guided_wop_checklist.md`

The returned quote report should also list those files in
`service_ready_bundle.local_files[]`.

## Step 6: Preflight and Quote the m-block Example

```bash
"$GENTLE_CLI" --state "$STATE" services project-preflight \
  @docs/examples/external_services/metabion_mblock_request.json

"$GENTLE_CLI" --state "$STATE" services project-quote \
  @docs/examples/external_services/metabion_mblock_request.json \
  --output-dir "$RUN_DIR/metabion_mblock_demo"
```

Expected difference from the oligo example:

- the product name in normalized line items should refer to m-block DNA
  fragments/libraries,
- warnings or required follow-up should mention biosafety review where
  applicable,
- direct submission remains unavailable.

## Step 7: Repeat the Same Review in the GUI

Open an already-built GENtle GUI and use the shared inspector. Do not open
vendor links or send the generated email draft during this rehearsal.

Then:

1. Open `Services -> External Services...`.
2. Press `Refresh Providers`.
   Press `Provider Config Doctor` and inspect the catalog sources as well.
3. Confirm Metabion appears in the provider picker.
4. Select `metabion` and `dna_oligo_single_tube`.
5. Paste the complete bundled oligo request JSON into the request editor after
   selecting the provider and service. `Use Selected Template` is a starter,
   not this exact fixture; changing the selection replaces the editor contents.
6. Press `Preflight`.
7. Press `Prepare Quote Handoff`.
8. Inspect the payload previews and warnings.
9. Set `Output dir` to the printed absolute `RUN_DIR` path plus a new
   `/metabion_oligo_gui` subdirectory. Do not paste the literal `$RUN_DIR`:
   this GUI field does not expand shell variables.
10. Press `Export Handoff Bundle`.
11. Confirm the generated files are listed under `Bundle files`. Repeat with
    `dna_fragment`, the m-block request, and a fresh output subdirectory.

The GUI should not contain provider-specific business logic. If the GUI and CLI
disagree, treat that as a bug in shared-shell plumbing or presentation, not as
a vendor-specific GUI feature to patch separately.

## Step 8: Optional Project Policy Overlay

Provider behavior is meant to be locally maintainable without changing the
GENtle source tree. To rehearse that path, copy the example overlay into a
temporary project-local catalog directory, not your real project's policy:

```bash
mkdir -p "$GENTLE_PROJECT_ROOT/.gentle/catalogs/external_service_providers.d"
cp docs/examples/catalogs/external_service_providers_project_overlay.json \
  "$GENTLE_PROJECT_ROOT/.gentle/catalogs/external_service_providers.d/metabion_local_policy.json"
"$GENTLE_CLI" --state "$STATE" services providers doctor
```

Expected outcome:

- the doctor reports the project-local overlay source,
- later provider ids override earlier provider ids,
- GUI and CLI provider rows update from the same catalog result.

This disposable overlay is not a recommendation to change lab policy. A GUI
must use the same project-root override to see it.

## Step 9: What to Send to Metabion Later

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
- exported quote bundles contain `quote_report.json` and generated local
  payload files,
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
