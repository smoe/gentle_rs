---
name: gentle-tss-collection
description: >-
  Route annotated transcription-start previews, approved TSS-window
  materialization, collection inspection, GUI handoff, and explicit registry
  removal through GENtle's deterministic gentle-cloning runtime.
version: 0.1.0
author: GENtle project
license: MIT
tags: [tss, promoters, transcript-starts, sequence-collections, provenance]
metadata:
  openclaw:
    requires:
      bins: [python3]
      env: []
      config: []
    always: false
    homepage: https://github.com/smoe/gentle_rs
    os: [macos, linux]
    install: []
    trigger_keywords:
      - gentle tss
      - transcription start sites
      - annotated transcript starts
      - tss collection
      - tss windows
      - inspect tss collection
      - open tss collection
      - forget tss collection
---

# GENtle TSS Collections

This descriptor-only skill routes the public tutorial 08.15 collection
lifecycle through `gentle-cloning`. GENtle remains the sole owner of annotation
geometry, exact-start grouping, preview digests, sequence derivation,
collection validation, UI intents, and persisted project state.

## Procedure

1. Start from one loaded, anchored, transcript-annotated locus. Treat the
   inventory as loaded-locus annotation scope, not all possible starts and not
   experimental initiation evidence.
2. Preview with `promoters tss-inventory` using a caller-visible request JSON.
   Report the exact available rows, exclusions, warnings, selected collection
   ID, and `approval_sha256`. Never invent TSS IDs or a digest.
3. Materialize only through the confirmation-gated route. Its reviewed request
   must repeat the inventory request, exact preview digest, and explicit
   `selected_tss_ids`. Execute only the stored proposal after approval of its
   exact digest.
4. Use `promoters tss-list` only for registry discovery. Describe every row as
   not checked until `promoters tss-collection ID` validates the members.
5. Hand `ui open tss-view --collection ID` to a GUI host. A headless
   `applied=false` receipt is a handoff, not proof that windows opened.
6. Forget registry metadata only through the confirmation-gated named route.
   State explicitly that member sequences and lineage remain.

## Boundary

- Do not infer promoter activity, occupancy, preferred starts, or biological
  completeness from annotations or collection membership.
- Do not silently replace an existing collection, repair stale members, or
  reuse a changed preview under an older approval.
- Preserve GENtle warnings and typed results verbatim; the outer layer adds
  routing, approval, and reproducibility evidence only.

