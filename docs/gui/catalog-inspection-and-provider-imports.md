## Catalog Inspection And Provider Imports

<!-- manual-counterparts:start -->
Related CLI reference: [Services and Resource Catalogs](../cli/catalogs.md).
These topic links do not imply identical feature coverage across interfaces.
<!-- manual-counterparts:end -->

Catalog inspection and guarded provider imports are currently GUI **Shell**
routes, not a new catalog browser: use `catalog doctor` / `catalog search`,
then explicitly review `catalog updates preview` before approving `apply`.
See [the scope, approval and recovery contract](../catalog_updates.md).
External-services catalog views use the same overlays; refresh them explicitly
after a configuration change. Local stock/prices and project undo are separate.

<!-- manual-page-nav:start -->

---

[GUI Manual](../gui.md) | [Previous: Rack View](rack-view.md) | [Next: Screenshot Capture Policy](screenshot-capture-policy.md)
<!-- manual-page-nav:end -->
