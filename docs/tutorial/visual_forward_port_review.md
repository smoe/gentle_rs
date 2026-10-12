# Glen Visual Placement Review

Scope: semantic forward-port of
`d70bf25def8424aea8bf59cb851f03866c4aa7cb` from
`smoe-bot/bot/tutorial-informative-visuals-20261007`, reviewed against main
`5e6201a69a599cba87491c86cfd7311cb6fee173`. No cherry-pick, asset edits,
new capture, scientific verdict or inherited current-candidate acceptance.
Current prose and canonical source metadata take precedence over old captions.
Existing lesson-level review dates are not renewed by this placement-only review.

The nineteen proposed placements have the following individual dispositions.
Twelve assets are placed in the proposed lessons, including three already there.
The five dotplot assets are already embedded in 02.03; omitting duplicates from
02.04 is not a missing visual. Two other proposed placements are not retained.

| Proposed Lesson | Asset Under `docs/` | Disposition And Teaching Limit |
| --- | --- | --- |
| 01.01 | `screenshots/tp73_inner_agent_evidence_viewer.png` | Do not add: DNA-viewer crop shows no agent exchange/consent. Keep the newer source-bound Agent Assistant capture already in section 4. |
| 01.01 | `figures/tp73_outer_agent_evidence_viewer.svg` | Add to section 2: explicit saved-state example, local 1..1200 viewport; synthetic proof tracks are not experimental evidence. |
| 03.05 | `figures/gibson_single_insert_protocol_cartoon.svg` | Add to step 7: five mechanism stages and both junctions; template overhangs are not the plan's 30-bp overlaps. |
| 03.05 | `screenshots/tutorial_gibson_specialist_testing_08_lineage_graph.png` | Keep existing step-8 historical capture and date/acceptance limits unchanged. |
| 03.07 | `figures/gibson_single_insert_arrangement_gel.svg` | Keep existing step-5 stored lane-order example unchanged. |
| 03.07 | `figures/gibson_single_insert_lineage.svg` | Add to step 6, corrected: cloning hub and products only. The old proposed caption incorrectly implied an arrangement node appears in this SVG. |
| 03.08 | `figures/gibson_single_insert_rack_isometric.svg` | Keep existing step-5 isometric occupancy example unchanged. |
| 03.08 | `figures/gibson_single_insert_rack_topdown.svg` | Add to step 6: top-down hero occupancy, not fabrication geometry or validated printable hardware. |
| 04.03 | `figures/qpcr_assay_protocol_cartoon.svg` | Add to step 4: ROI and independent primer/probe windows; geometry is not specificity. |
| 04.03 | `figures/tp73_as3_shared_cdna_qpcr.svg` | Add to step 6: separate TP73-AS3 tester example, 100-bp product in 3/3 annotated templates, not new TP73 design outputs. |
| 04.03 | `figures/tp73_as3_nr187362_cdna_qpcr.svg` | Add to step 6: separate TP73-AS3 tester example, 147-bp product in NR_187362.1; other rows lack the forward hit, not proven genome-wide specificity. |
| 10.02 | `figures/pgex_blunt_sequence_confirmation_hero.svg` | Do not add: actual template is `readme.pcr_blunt_unknown_insert_degenerate`, illustrating degenerate PCR then vector-primer Sanger, not the claimed checkpoint/variant report. Keep the known-construct ABI lesson focused. |
| 02.04 | `screenshots/tutorial_cdna_genomic_02_fetch_cdna_dialog.png` | Keep in 02.03 step 1, figure 2; do not duplicate the accession dialog in the online replay. |
| 02.04 | `screenshots/tutorial_cdna_genomic_03_prepare_genome_dialog.png` | Keep in 02.03 step 2, figure 3; 02.04 retains its catalog/cache/network contract and link to 02.03. |
| 02.04 | `screenshots/tutorial_cdna_genomic_04_retrieve_tp73_genomic_dialog.png` | Keep in 02.03 step 2, figure 4; no duplicate retrieval screenshot. |
| 02.04 | `screenshots/tutorial_cdna_genomic_07_dotplot_window.png` | Keep in 02.03 step 4, figure 7; no duplicate workspace screenshot or new online execution claim. |
| 02.04 | `screenshots/tutorial_cdna_genomic_09_dotplot_as_node.png` | Keep in 02.03 step 5b, figure 9; no duplicate lineage screenshot or current-native verdict. |
| 03.04 | `screenshots/screenshot_GUI_help_tutorial_testing_gibson.png` | Add via source JSON at step 3: historical Help orientation only, not proof of current automatic opening or starter state. |
| 05.02 | `screenshots/screenshot_GUI_prepare_reference_genome.png` | Add via source JSON at step 1, not the old step 2. Explicitly warn that 113 is selected; the chapter requires 116 and fresh local cache status. |

## Regeneration And Verification

Change generated chapters only through `sources/*.json`, then regenerate the
catalog, manifest and chapter projections with `gentle_examples_docs`. Preserve
unrelated generated scientific artifacts and their historical generator hashes.
The existing local helper identifies as `150a339074e7a2d282567d791a8f50834b70ff9c`;
its `src/workflow_examples.rs` and `src/bin/gentle_examples_docs.rs` are byte-identical
to this main baseline (Git blobs `913482e6ae4081a5674105af6c1b4d9d517dcc5a`
and `4343f0f90f48363b3d5bb8d6190a8aa008ac62e3`). It is a projection tool here,
not a same-SHA Rust/runtime acceptance binary. Online execution stays disabled.
Full regeneration into a temporary directory also produced unrelated scientific
artifact differences with this older binary. None are imported or waived: only
the two affected chapter projections, their report records and chapter hashes
are updated. All historical artifact bytes and ledger hashes remain unchanged.

Reuse the existing visual-contract and LF/CRLF checkout tests. The added
repository regression pins the accepted placements, historical-only screenshot
status, preparation step and exclusions; it does not interpret experimental
results or certify live GUI behavior. The two legacy capture dates above retain
the old branch's date-only declarations, not new hash-bound capture evidence.
Fresh same-source tutorial/Rust/native acceptance remains Glen/CI's work.
