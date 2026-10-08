# VKORC1 Reporter Preview Provenance

## Status And Origin

These are **historical, unverified WIP previews**, retained from Glen's
`40430f98aac1830c6cd6176b253fde190423fb92`, branch
`bot/vkorc1-promoter-reporter-wip-20261008`, based on `f523a757`.
They illustrate tutorial 08.04; they are not parser fixtures, wet-lab validation,
current-candidate GUI acceptance or a pharmacogenomic recommendation.

The two reports describe a genomic-forward GRCh38 slice on chromosome 16,
`31093368..31099368`, with rs9923231 at `31096368` and a reverse-strand VKORC1
transcript. The reported reference is C, alternatives A/G/T; the tutorial
explicitly chooses T, rather than deriving a base from the gene strand.
The selected fragment is local `2412..3501` (0-based, end-exclusive).
The accompanying backbone is the explicitly synthetic
`data/tutorial_inputs/gentle_mammalian_luciferase_backbone_v1.gb`, not a
sequence-exact functional or commercial plasmid.

Glen's commit did not supply a complete input-state, binary/toolchain/lockfile
receipt or native acceptance ledger for these previews. Do not infer them from
JSON timestamps, run names, SVG titles or file presence. The exact C/T insert
sequences, matched boundaries, promoter-to-luc2 orientation, insertion site,
ligation-product choice and actual topology require an independent replay.
Drawing a product circularly does not establish circular molecular topology.

## Retained Bytes

The five files below are unchanged from Glen's commit. SHA-256 identifies the
retained bytes, not scientific correctness. Scoped LF checkout rules and
`scripts.test_tutorial_checkouts` preserve them across LF/CRLF hosts.
The tutorial, command script, `report.md` and `result.json` are maintained
instructions/status projections rather than a retained execution receipt.

| File | SHA-256 |
| --- | --- |
| `variant_promoter_context.json` | `1c074ca878e4aeb063a358e6cfeffd68705e1f9b7b052a6ffb4695ebe84e6793` |
| `promoter_reporter_candidates.json` | `2142af4f9ebe8cf2b6155721671613f7e83eba52298e2110ab6cac6a135d72a0` |
| `vkorc1_rs9923231_promoter_context.svg` | `edc5da3adec070749b689ccc0448b691793ac0660bca77d97832b5ed855cb88b` |
| `vkorc1_rs9923231_reporter_reference.svg` | `b1ae7419b44a4f2990f014d3315b82c098c3a953dad3718fce61f4e0d60b4c30` |
| `vkorc1_rs9923231_reporter_alternate.svg` | `ba011b47bee4d02b4f6f809d45b72b3bbad4b4ac28e39e82f6f8471297292de0` |

## Reproduction And Next Review

Recover exact historical bytes with `git show 40430f98:<repository-relative
file path>`; do not regenerate them and silently replace their provenance.

For a fresh replay, use the current committed SHA and the tutorial/workflow
in a disposable checkout, an explicit new state file and a separate output
directory. The supplied script writes repository-relative output paths: do
not run it over the retained artifacts in a working development checkout.
Prepare the named reference, inspect the fetched placement/ref/alt qualifiers,
then select the reviewed genomic-forward T explicitly. The initial dbSNP
lookup and reference preparation require network access; later operations are
local. Recheck the current candidate's geometry instead of assuming historical
fragment coordinates remain the recommended result.

Retain the new source SHA, lock/toolchain/profile and binary hashes, reference
preparation manifest, raw dbSNP response, input/output sequence digests, typed
reports, exact commands, SVGs and GUI screenshots separately. Compare reference
and alternate sequences: identical length/boundaries and exactly one reviewed
base change. Confirm actual construct orientation and topology before any
experimental use. A synthetic in-process test cannot replace this real-data
and native-GUI acceptance.
