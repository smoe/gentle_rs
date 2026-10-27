## External Services Window

<!-- manual-counterparts:start -->
Related CLI reference: [Services and Resource Catalogs](../cli/catalogs.md).
These topic links do not imply identical feature coverage across interfaces.
<!-- manual-counterparts:end -->

Access:

- Main menu: `Services -> External Services...`
- Command palette: `External Services`

The External Services window is a GUI inspector for shared-shell service
contracts, not a separate ordering implementation.

- `Refresh Providers` runs `services providers list` and populates the provider
  picker from `gentle.external_service_provider_catalog.v1`.
- `Provider Config Doctor` runs `services providers doctor` against the active
  built-in/system/user/project overlay discovery chain.
- The provider and service-kind pickers are catalog-driven. Metabion appears
  when configured, with initial DNA oligo and m-block fragment handoff rows.
- The request editor holds one editable
  `gentle.external_service_request.v1` JSON payload. `Use Selected Template`
  creates a provider-neutral starter request for the selected provider/service.
- Changing the selected provider or service kind refreshes that editable
  starter request and clears stale preflight/quote previews.
- The project-source route block calls `services route-project-source` for
  active sequence/span, persisted oligo order forms, or primer report pair
  ranks. The selected delivery-route candidate can be copied into the same
  editable request JSON for review.
- Primer report pair ranks are persisted as an oligo order form before
  classification, so exact duplicate procurement rows remain visible and
  quote handoff stays blocked until duplicate review is recorded.
- `Preflight` runs `services project-preflight`; `Prepare Quote Handoff` runs
  `services project-quote`.
- `Export Handoff Bundle` runs `services project-quote --output-dir DIR` and
  writes deterministic local handoff files, with generated paths reflected in
  the returned quote bundle's `local_files[]`.
- Quote output previews inline handoff payloads such as redacted request JSON,
  normalized line-item JSON/CSV, email draft markdown, and guided WOP checklist
  when the active provider config supplies those mappings. When a bundle has
  been exported, the same panel also lists the generated files.
- The window never scrapes portals, submits orders, looks up credentials, or
  stores PO/account/shipping/billing details in project state.

This is the intended place to review whether GENtle represents provider
products and handoff constraints fairly before deeper provider collaboration or
future direct integrations.

Tutorial companion:

- [`docs/tutorial/09-01_metabion_external_service_handoff_gui_cli.md`](../tutorial/09-01_metabion_external_service_handoff_gui_cli.md)

<!-- manual-page-nav:start -->

---

[GUI Manual](../gui.md) | [Previous: Screenshot Capture Policy](screenshot-capture-policy.md) | [Next: Configuration Window](configuration-window.md)
<!-- manual-page-nav:end -->
