# Synthetic Cloning Planning Profiles

Origin: hand-crafted on 2026-10-19 for tutorial 09.04, not a laboratory's
inventory or supplier quotation. No sequence or commercial account is used.

Recreate `cloning_in_stock.json` with restriction-enzyme and ligase allowances
of 2 and 7 EUR and a synthetic thermocycler rate of 5 EUR/hour, all available.
Copy it to `cloning_waiting.json`, mark ligase unavailable with two business
days of procurement and set thermocycler queue to 0.5 business days. IDs and
notes deliberately say `DEMO_DO_NOT_ORDER`.

Use: tutorial 09.04, its optional real-CLI regression in
`scripts/test_specialist_tutorials.py`, and shared planning-contract tests.
Expected known cost subtotal is 39 EUR in both profiles: routine base 0,
materials 9 and machine 30 (legacy six-hour occupancy basis). Waiting moves
7 from stock-draw value into the purchasing proxy; it does not buy anything.
The preserved ranking/time-guardrail proxy is 6 h in stock and 90 h with
waiting; displayed synthetic critical paths are 3.2 h and 87.2 h respectively.
Additional unknown/unpriced classes, if introduced by a catalog change, must
remain visible and prompt a tutorial review, not be silently called free.

## Vector, Decision And Catalog Exercises

Tutorial 09.04 reuses the existing
`docs/examples/workflows/promoter_reporter_panel_planning_offline.json` and
its documented 08.09 fixtures. Its Python block only rebinds input paths and
disposable output paths; it does not duplicate validation, enzyme selection,
proposal generation or scheduling. The existing generated maps are teaching
SVGs, not native screenshots or commercial vector artwork.

The supplier table and decision sheet are presentation projections of saved
reports. The worksheet preserves unknowns and raw input hashes; it is not an
engine schema, verified cross-lane join or authorization. Recreate it by running
the chapter's literal Python blocks after retaining the named shared reports.
Source tests use explicitly synthetic temporary projections, including missing
values and LF/CRLF bytes, without inventing supplier evidence.

`supplier_catalog_update_exercise.md` is the advanced, optional part of 09.04,
not another catalog entry. Its complete provider is hand-crafted synthetic
data with an invalid example URL. Recreate it with the literal Python block;
then review each preview before add/apply, stale-reuse refusal and explicit
withdrawal in the disposable catalog. The optional prebuilt-CLI companion in
`scripts/test_specialist_tutorials.py` replays the literal main/advanced commands;
it never builds Rust, downloads data, invokes a model or contacts a vendor.
