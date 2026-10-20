# 09.04 Advanced: Maintain A Disposable Supplier Catalog

Return to [tutorial 09.04](../../tutorial/09-04_cloning_time_cost_and_protocol_suppliers.md).
This optional configuration exercise is not needed to choose a vector or
review a protocol. First complete the main chapter and retain its decision
sheet. Reuse its shell session, `g` function and disposable `$RUN`/catalog roots;
do not apply this to a study or a real supplier.

Keep five lanes separate: product identity, commercial offer, protocol
document, reviewed procedure and private local observations. The record below
is hand-crafted synthetic teaching data, not Promega, an offer or a reviewed
procedure. Its invalid example domain is never contacted.

## Preview Before Approving

```bash
python3 - "$RUN/provider-update.json" <<'PY'
import json, sys
record = dict(provider="synthetic_tutorial_supplier", display_name="Synthetic tutorial only",
              support_status="catalog_only", website_url="https://example.invalid",
              dashboard_url="https://example.invalid", api_documentation_url=None)
for key in ("channels", "product_templates", "validation_rules", "default_delivery_hints",
            "default_purification_hints", "default_qc_hints", "required_followup",
            "account_enablement_notes", "warnings"):
    record[key] = []
record["capabilities"] = [dict(service_kind="synthetic_protocol", track="dna",
    display_name="Synthetic", quote_handoff_supported=False,
    direct_api_documented=False, direct_api_implemented=False,
    supported_submission_modes=[], status_tracking="none", artifact_kinds=[], notes=[])]
bundle = dict(schema="gentle.provider_update_bundle.v1", changes=[dict(
    action="replace", provider=record["provider"], base_record_sha256=None, record=record)])
with open(sys.argv[1], "x", encoding="utf-8", newline="\n") as handle:
    json.dump(bundle, handle, indent=2)
    handle.write("\n")
PY
g catalog updates preview "$RUN/provider-update.json" --scope project --output "$RUN/provider-preview.json"
```

Read `provider-preview.json` before approving: expect `add`, null previous
record hash, `protocol_document_evidence_unavailable`, `applicable=true` and
the exact destination under `$GENTLE_PROJECT_ROOT/.gentle/catalogs/`. Inspect
bindings and limitations, not just the digest. Preview writes only a fresh
artifact. Stop if a real supplier or unexpected destination appears.

After reviewing this synthetic addition only:

```bash
DIGEST="$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["preview_sha256"])' "$RUN/provider-preview.json")"
g catalog updates apply "$RUN/provider-preview.json" --approve-digest "$DIGEST" > "$RUN/provider-receipt.json"
g catalog search --entity provider --query synthetic_tutorial_supplier > "$RUN/synthetic-provider.json"
g catalog doctor
```

The receipt identifies the active raw-byte hash and retained history. Compare
the active catalog and history's `proposed.json` with that hash; `previous.bin`
must be absent because the project catalog did not previously exist. Existing
built-in suppliers must remain visible. No project undo history, stock, planning
estimate, scientific result or vendor account is changed.

## Reusing Approval Must Fail

This is a deliberate negative test, not an instruction to retry until success:

```bash
if g catalog updates apply "$RUN/provider-preview.json" --approve-digest "$DIGEST" > "$RUN/stale-output.json" 2> "$RUN/stale-error.txt"; then
    echo "ERROR: old approval was reused" >&2
    exit 1
fi
g catalog updates preview "$RUN/provider-update.json" --scope project --output "$RUN/conflict-preview.json"
```

Expect `stale_or_tampered_preview` in the retained error and a non-applicable
new preview containing `conflict`: the old null base cannot overwrite a now
existing provider. Verify that this negative test did not alter the active
catalog or project. A generic command failure is not proof of the right refusal.

## Explicitly Withdraw The Synthetic Addition

For a replacement, obtain the current **full** configuration and semantic
`record_sha256` from provider search, not the public listing projection.
For withdrawal use that current hash **without** a replacement record:

```bash
python3 - "$RUN/synthetic-provider.json" "$RUN/provider-withdraw.json" <<'PY'
import json, sys
rows = json.load(open(sys.argv[1], encoding="utf-8"))["rows"]
row, = [row for row in rows if row["id"] == "synthetic_tutorial_supplier"]
bundle = dict(schema="gentle.provider_update_bundle.v1", changes=[dict(
    action="withdraw", provider=row["id"], base_record_sha256=row["record_sha256"])])
with open(sys.argv[2], "x", encoding="utf-8", newline="\n") as handle:
    json.dump(bundle, handle, indent=2)
    handle.write("\n")
PY
g catalog updates preview "$RUN/provider-withdraw.json" --scope project --output "$RUN/withdraw-preview.json"
```

Review this new preview and its `withdraw` decision before approving:

```bash
WITHDRAW_DIGEST="$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["preview_sha256"])' "$RUN/withdraw-preview.json")"
g catalog updates apply "$RUN/withdraw-preview.json" --approve-digest "$WITHDRAW_DIGEST" > "$RUN/withdraw-receipt.json"
g catalog search --entity provider --query synthetic_tutorial_supplier > "$RUN/withdrawn-provider.json"
g services providers list > "$RUN/suppliers-after-update.json"
```

Expect no effective match for the synthetic ID, a v2 catalog with an explicit
`withdrawn_providers` tombstone, both retained transaction histories and all
original suppliers still present. Removing a row is not equivalent: it can
reveal a lower overlay. A later deliberate overlay can restore a provider.

Both writing commands require agent confirmation even if a model claims they
are read-only. Source changes, later `.d` masking or tampering require a new
preview. No-op-only previews are not applicable. See
[limits and interrupted-apply recovery](../../catalog_updates.md): never delete
receipts to make a retry pass, auto-replay history or treat a content hash as
authenticated vendor approval. Joined offers, retained/reviewed documents,
quantities/taxes/shipping, procedure extraction and inbox automation remain
later work. No Rust/native acceptance or supplier contact is implied.
