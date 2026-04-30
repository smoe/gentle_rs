# Reviewer Quickstart (Internal Preview)

This page is for internal colleagues evaluating preview builds.

## Internal-only warning

GENtle `0.1.0` preview builds are internal engineering milestones.

- Not a public/stable release.
- Not intended for production, clinical, or regulated workflows.
- APIs, UI behavior, and file formats may still change.

## Known limitations (preview scope)

- GUI polish and publication-quality visual tuning are still in progress.
- Primer workflows are implemented, but GUI-first PCR/primer node-specific interfaces are still evolving.
- Async BLAST baseline exists for shell/MCP/agent usage; multi-job orchestration UX is still being expanded.
- Some advanced workflows remain shell-first and are not yet represented by dedicated GUI dialogs.
- Preview defaults for adaptive/helical DNA-letter rendering are optimized for clarity, but exact visual preferences may still be tuned.

## 10-minute reviewer walkthrough

Run from repository root:

```bash
STATE=/tmp/gentle_preview_state.json
rm -f "$STATE"

# 1) Load + branch + reverse-complement workflow example
cargo run --bin gentle_cli -- --state "$STATE" workflow @docs/examples/workflows/load_branch_reverse_complement_pgex_fasta.json

# 2) Digest/ligation/extract baseline workflow example
cargo run --bin gentle_cli -- --state "$STATE" workflow @docs/examples/workflows/digest_ligation_extract_region_minimal.json

# 3) Fast deterministic primer-pair report (internal backend)
cargo run --bin gentle_cli -- --state "$STATE" primers design '{"DesignPrimerPairs":{"template":"pgex_fasta","roi_start_0based":100,"roi_end_0based":300,"forward":{"min_length":18,"max_length":22,"location_0based":null,"start_0based":80,"end_0based":130,"min_tm_c":50.0,"max_tm_c":72.0,"min_gc_fraction":0.2,"max_gc_fraction":0.8,"max_anneal_hits":5},"reverse":{"min_length":18,"max_length":22,"location_0based":null,"start_0based":280,"end_0based":340,"min_tm_c":50.0,"max_tm_c":72.0,"min_gc_fraction":0.2,"max_gc_fraction":0.8,"max_anneal_hits":5},"min_amplicon_bp":120,"max_amplicon_bp":400,"max_pairs":5,"report_id":"preview_pgex_fast"}}' --backend internal
cargo run --bin gentle_cli -- --state "$STATE" primers show-report preview_pgex_fast
```

Expected outcomes:

- Workflow command returns operation records with created IDs.
- Primer design returns report ID `preview_pgex_fast` with `effective_backend.used = internal`.
- `primers show-report` returns the persisted report payload.

## GUI checks for this preview

1. Open sequence window and zoom linear map.
2. Confirm top-right routing diagnostics update (`DNA MODE`, density/columns/reason).
3. Confirm `Configuration -> Graphics` defaults:
   - mode: `Auto adaptive`
   - compressed letters enabled
4. Use `Apply` in Configuration to persist graphics defaults across restart.

## Reporting feedback

When reporting issues, include:

- command/workflow used,
- sequence ID and view span,
- screenshots (if visual),
- exact warning/error text,
- whether behavior reproduced after restart.
