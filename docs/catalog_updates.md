# Catalog Evidence And Guarded Provider Updates

This extends the existing catalog/overlay system, not the planning engine.
All commands below use the shared shell. Direct `gentle_cli catalog ...` and
GUI Shell take the same path; MCP/JS/Lua use an explicitly selected shell
candidate through their existing agent-plan executor. There is no new named
MCP tool, catalog GUI, scheduler, inbox watcher or network fetch.

## Five Separate Lanes

| Lane | Current inspection | Boundary / next contract |
| --- | --- | --- |
| Product identity | Existing helper/vector/reference cards and semantic record hashes | Not necessarily purchasable products. Preserve accession/sequence, host, selection, license/MTA and biosafety metadata where supplied; no seller join by name or alias. |
| Commercial offer | Unavailable, not zero offers or a free product | Future explicit product-ID join to seller/SKU, region/package, currency, quote-required/public price, availability/delivery, tax/shipping, observation date and evidence. No package-aware cost model yet. |
| Protocol document | Dated official HTML/text-PDF **references**, not retained documents | Future immutable retained bytes with URL/local origin, title/revision/date, retrieval time, SHA-256, format/license and applicability claims. `checked_on` is not a document hash or freshness verdict. |
| Normalized procedure | No reviewed vendor procedure; existing synthetic bench models remain separate | Future reviewed steps cite exact source document hash/passages, quantities/equipment/times, dependencies, controls, outputs/checkpoints. Machine extraction is a draft, never self-authorizing. |
| Local profile / observations | Private; doctor/search/import do not inspect or modify it | Existing stock, negotiated unit values, delays/equipment stay local. Future observed duration/failure history must not overwrite vendor facts or silently change ranking. |

```text
catalog doctor
catalog search --entity product --query luciferas
catalog search --entity provider --query promega
catalog search --entity protocol --query luciferase --limit 20
```

`gentle.catalog_doctor.v1` includes the existing provider/helper doctors and
lane availability, not a global pass. A failed source is diagnostic, not an
empty biological/commercial result. `gentle.catalog_search.v1` uses a literal
case-insensitive substring over the selected record/reference JSON, sorted by
ID with a 1..100 limit (default 50), `total_matches` and explicit truncation.
The query is required, nonblank and at most 256 UTF-8 bytes. Product rows are
`helper_catalog_identity`, including host/reference entries. Provider rows
expose the full effective config, semantic `record_sha256`, source scope/path
and raw source hash, enabling complete replacement without guessing defaults.
Protocol rows retain format/date and `document_sha256=null`;
`procedure_review=not_evaluated`. Country is supplied company metadata, not
shipping coverage. A match does not join the lanes or establish identity.

## Provider Bundle Admission

The first writable slice admits **provider configurations only**. It does not
create products, offers, reviewed procedures or price observations. Complete
existing-record overlays remain valid; arbitrary partial JSON patches do not.

`gentle.provider_update_bundle.v1` contains 1..64 unique normalized IDs:

```json
{
  "schema": "gentle.provider_update_bundle.v1",
  "changes": [{
    "action": "replace",
    "provider": "synthetic_supplier",
    "base_record_sha256": null,
    "record": "REPLACE THIS STRING WITH A COMPLETE PROVIDER CONFIG OBJECT"
  }]
}
```

The illustrative string above is deliberately **not executable**. Use the full
object from provider search (all serialized config fields, not the public
`services providers list` projection). Unknown fields, partial/default-filled
replacements and mismatched IDs fail admission. An addition has a null base;
a replacement uses the current effective record hash. Withdrawal is explicit:

```json
{
  "schema": "gentle.provider_update_bundle.v1",
  "changes": [{
    "action": "withdraw",
    "provider": "synthetic_supplier",
    "base_record_sha256": "sha256:REPLACE_WITH_CURRENT_RECORD_HASH"
  }]
}
```

It writes a `withdrawn_providers` tombstone in the selected full config rather
than deleting the overlay and accidentally revealing an older provider. A
nonempty tombstone list requires `gentle.external_service_provider_config.v2`;
preview upgrades that file explicitly. Legacy v1 without withdrawals stays
readable and byte-compatible when serialized. Older GENtle must reject v2
instead of ignoring the withdrawal and reviving a lower record. A
later explicit full-record overlay can deliberately restore it. A file may
not both withdraw and define the same normalized ID. Provider resolution,
catalog listing and eligibility consumers use the same overlay merge.

## Preview, Review, Apply

```text
catalog updates preview provider-update.json --scope project --output preview.json
catalog updates apply preview.json --approve-digest sha256:EXACT_PREVIEW_DIGEST
```

Scope is **required** and only `user` or `project`. Destinations are the
existing discovered `catalogs/external_service_providers.json`, including
`.gentle/catalogs` under the configured project root. Built-in/system sources
are read-only. No arbitrary target path is accepted. Preview writes only a
fresh explicitly named artifact, refuses overwrite and catalog-source
destinations, and does not create its output parent directory. The source
bundle itself is retained in the preview with its raw SHA-256; parse and hash
refer to those same bytes, including LF/CRLF differences.

`gentle.provider_update_preview.v1` binds target scope/path, source candidates
(including absence and sorted `.d` membership), raw source hashes, old target
hash/absence, proposed bytes, full bundle, exact planned output paths and decisions. Per-record actions
are `add`, `replace`, `withdraw`, `no_op` or `conflict`; base and later-overlay
masking conflicts block the **whole** apply. A preview with only no-ops is not
applicable and does not rewrite a catalog. Validation errors are explicit
`catalog_update:CODE: detail` failures, not successful empty updates. The preview
lists protocol reference dates, missing evidence and unevaluated freshness;
there is no arbitrary age threshold or live document verification. Impact
statements name future provider discovery/capability/handoff use. Specific
saved planning results, biological products and protocol applicability are
**not resolved** by this tranche; saved reports are not rewritten.

Apply requires the exact preview digest and recomputes all preview fields from
the retained bundle/current sources, before writes and again under an immediate
exclusive writer lock. Tampering with proposed data and honestly recomputing a
digest is not enough: it must still reproduce the preview. A newly changed
source, directory member, scope or target requires a fresh preview. Both
artifact-writing commands are mutation-classified for agents; plan execution
requires explicit confirmation even if candidate flags falsely claim read-only.
Direct CLI/GUI Shell invocation is an explicit action. `state_changed=false`
means no project-state edit, **not** no filesystem write. Catalog changes are
external configuration, not project undo/redo history.

Limits: 512 KiB bundle, 4 MiB per source/artifact read and proposed catalog,
128 catalog files, 32 MiB combined sources. Unknown fields in existing sources
also fail update admission rather than being dropped by a rewrite. Raw parent
traversal and symlink paths are refused before normalization. Catalog roots
must therefore use non-symlink paths (on macOS, prefer canonical `/private/...`
temporary paths). This is a cooperating-writer lock, not protection against a
malicious account replacing path components concurrently; no network/untrusted
shared-filesystem transaction guarantee is made.

## Atomicity And Recovery

The `.provider-update.lock` file remains in the catalog directory; its OS lock
is released when the process/handle exits. Busy writers fail immediately.
Before replacing the active file, apply retains a new immutable directory at
`.provider-update-history/UPDATE_BINDING_DIGEST/` with `previous.bin` if previously
present, `proposed.json`, `preview.json` and `receipt.json`. The history path is
bound in the preview from source/bundle/target hashes, without a self-referential
preview digest. Absence is not an
empty file. It flushes written files; the sole commit point is an atomic
same-directory temporary-file replacement of the active catalog. No catalog
changes occur if validation or recovery staging fails. This is not a claim of
multi-file atomicity or power-loss durability on every filesystem.

`gentle.provider_update_receipt.v1` is conditional rather than a falsely final
success flag in a second file: compare the target's **raw** hash with its
`applied_target_sha256` and `previous_target_sha256`. Applied hash means the
replacement took place; previous hash/absence means it did not; anything else
requires checking later changes. The returned receipt is emitted only after
successful replacement. A crash before replacement can leave prepared history;
retry refuses an existing history directory rather than overwriting evidence.
An operator must inspect/reconcile it; never auto-replay, auto-rollback or
delete retained evidence. Old historical receipts describe their transaction,
not the current target after subsequent updates.

Content hashes are **not signatures** or authenticated vendor statements. This
slice admits hash-bound local bundles; signature trust policy, transactional
product/offer/document imports, field-level patches, resolved impact graphs,
reviewed procedure promotion and quarantined inbox discovery are later work.
No path auto-applies a bundle, downloads a sequence/document, contacts a
supplier, approves a quote, submits an order or changes scientific results,
planning times, scores or guardrails.

## Verification Handoff

Rust checks are authored, not locally executed under the owner's restriction:

```bash
cargo test --locked --lib provider_update
cargo test --locked --lib catalog_
cargo test --locked --lib provider_config
cargo test --locked --lib catalog_mcp
cargo test --locked --test capability_registry_parity
cargo test --locked --test parity_matrix_freshness
cargo test --locked --test mcp_capability_surface
cargo test --locked --test glossary_cli_usage
cargo check -q --locked
python3 -m unittest scripts.test_specialist_tutorials scripts.test_tutorial_checkouts
```

Glen/CI should run these on the integrated SHA, including native Windows lock,
atomic replacement/read-handle and raw-path boundaries. Test direct/shared
shell and agent-plan parity with the same isolated source roots. Review
tutorial 09.04's disposable synthetic import separately from commercial,
scientific and GUI acceptance. No live inner agent is needed for parser/
approval parity; no such acceptance is claimed here.
