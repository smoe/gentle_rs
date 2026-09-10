# TFBS Track Panel Fixture

`jaspar_30_track_panel.json` is an unmodified copy of the real input panel at
Git revision `a106cbbd5223f8c4be55a7d846b8416bc4b3ae33`, path
`docs/examples/regulatory_region_comparison/five_target_tss_integrated/jaspar_30_track_panel.json`.
It is input data, not a runtime preset or a biological validation result.

Retrieve the exact bytes with:

```sh
git show a106cbbd5223f8c4be55a7d846b8416bc4b3ae33:docs/examples/regulatory_region_comparison/five_target_tss_integrated/jaspar_30_track_panel.json
```

The tests in `src/tfbs_track_panel.rs` parse its 30 tracks, preserve all 28
factor identities and the adjacent three-matrix group, and resolve the panel
against the locally bundled JASPAR registry without replacing the runtime
registry. Mutant inputs and the small two-track examples are constructed only
in memory. The fixture is not scored and does not require the real gene FASTAs.
