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
