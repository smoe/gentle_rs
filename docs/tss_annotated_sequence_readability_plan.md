# Readable TSS Sequence Exports And Inspection

Status: export/SnapGene work remains a proposal, not Claude-reviewed. The first
native annotated-window consumer from section 3 is implemented separately:
see [GUI usage and boundaries](gui.md#tss--regulatory-dna-display). The reviewed
report-attachment slice now adds validated full TFBS curves and imported motif
hits. Native-view SVG now exports the selected span/lanes with source binding,
gap-aware curves and hover details. Approved Slice B adds explicit local scoring
as independent, provenance-bound, gap-aware lanes using the shared engine and
validated inline DNA; see [local scoring](gui.md#compute-local-tss-scores).
New dynamic database querying and reporter rows remain follow-ups; locus TSS
selection uses the [TSS workspace](tss_workspace.md).
Source review: `6bc7c5fc37b178b972dc19a272c0fed9df2a2c73`, 2026-09-14.
Evidence: the user-supplied `CD44_selected_TSS_annotated (1).embl` and two
SnapGene screenshots. The supplied file and screenshots are not copied into
the repository. No live SnapGene import or GUI acceptance was performed.
This plan does not expand the `.10` release gate.

## Diagnosis

The unreadable names originate in GENtle, not in SnapGene. The shared
[annotated record builder](../src/tss_profile_genbank.rs#L68) names each record
`TSS_` plus a truncated hash of its promoter ID. Both GenBank and EMBL reuse
this name. A GUI-only title change would not fix another application's import
list, tabs or exported filenames.

The supplied records describe these annotation-derived TSS candidates:

| Current name | Reference and genomic TSS, 1-based | Exported genomic span, inclusive | Local TSS |
| --- | --- | --- | --- |
| `TSS_c16adb802c09` | GRCh38, chr11:35,138,721, strand + | 35,138,221..35,138,921 | 501 |
| `TSS_e1f17c03116b` | GRCh38, chr11:35,139,171, strand + | 35,138,671..35,139,371 | 501 |

Both are 701-bp transcript-oriented genomic windows, covering -500..+200 bp
relative to their selected TSS. Naming them by local base 501 would not
distinguish them. Both `ID` lines explicitly declare `linear`. The screenshot
shows a circular-display checkbox selected, but the file does not request
circular topology; whether this is remembered UI state or importer behavior
needs a separate interoperability check.

GENtle's existing single-sequence `load_from_file` path currently selects the
last record of multi-record EMBL/GenBank files (`Vec::pop`), although the parsers
return all records. A shared multi-record import/chooser is therefore a separate
required follow-up for conveniently opening both CD44 windows; native TSS mode
inspects the record that is actually loaded, not every record in the source file.
Do not change `LoadFile` silently or mistake parser coverage of both records for
GUI import coverage of both records.

There is also a data-presentation problem, not merely a font problem:

| Record | Total features | Stored motif-peak features | Signal-interval features | Exon features |
| --- | --- | --- | --- | --- |
| First | 570 | 300 | 266 | 2 |
| Second | 781 | 300 | 476 | 1 |

Each record also has source and TSS features; the second additionally has a
coding-segment annotation and a translation-start marker. The first record's
two exon features have identical displayed and original bounds but different
transcript IDs. Its comments explicitly say geometry is supplied for 2 of 78
TSS-linked transcript models: this must not become a claim that all 78 were
inspected or have identical exon structure.

Signal labels concatenate sample names with long interval/source identifiers.
Every signal step and retained motif maximum/peak becomes a `misc_feature`;
SnapGene then places them among ordinary annotations. Restriction-site labels
add another competing layer. This is not a usable default regulatory view.
All 300 motif values per record are positive and labelled
`llr_background_tail_log10`; hiding negative scores alone cannot fix this
crowding. These annotations are not all possible sites or measured binding.

## 1. Readable Identity And Short Labels

Recommended first implementation, independently useful before any new GUI:

- Generate a coordinate-led record name from verified `TssGeometry.tss_1based`,
  not from the interval start or local TSS index. Candidate short names for
  these records are `CD44_TSS35138721` and `CD44_TSS35139171`.
- Put the full identity in the definition and GUI title, for example
  `CD44 | TSS chr11:35,138,721 (+) | GRCh38 | -500..+200 bp`. Distinguish
  annotation-derived TSS candidates from experimentally mapped starts.
- Keep promoter IDs, gene IDs, input hashes and report bindings unchanged as
  machine identity/provenance. Record the display-name-to-promoter mapping in
  the export inventory. Human names must never become the scientific join key.
- Handle name collisions deterministically, including same-position records
  from different references, contigs, strands or window policies. Prefer
  informative disambiguation; use a short suffix only when necessary. Do not
  silently truncate away the coordinate or rely on record order. Confirm name
  handling with older SnapGene as well as GENtle; these are local identifiers,
  not fabricated public accessions.
- Label the point feature `TSS 35138721` and retain full genomic, local and
  TSS-relative coordinates in its details. Keep the molecular topology linear.
- Shorten ordinary feature labels using existing structured metadata, such as
  factor/mark, condition and sample. Put full source IDs, hashes and interval
  IDs in notes, not in the display label. If a field is missing, do not extract
  a supposedly authoritative value from an arbitrary filename.
- For motif labels, use the exact matrix's bound factor label plus accession
  and a concise prediction qualifier. Preserve raw score and score kind in
  details; do not relabel background-tail values as LLR bits. Imported
  jaspar-mapping hits remain a distinct evidence source, not local score peaks.

Implementation seams: `src/tss_profile_genbank.rs`,
`src/annotated_sequence_io.rs`, `src/tss_profile_export.rs`, and the existing
TSS export request/index contracts. The first increment changes presentation
only: no feature dropping, score filtering, sequence retrieval or rescoring.

### Preserve Existing Receipts

The [current verifier](../src/tss_profile_export.rs#L2288) regenerates GenBank
and EMBL bytes from the report. Changing names/labels without versioning that
projection would make earlier intact bundles fail verification.

Add an explicit annotated-projection policy/version bound to the export
request and inventory/receipt. Old recorded requests without that field retain
legacy naming/serialization; newly selected readable exports record their new
policy. Do not infer a layout version from commit hashes or silently reinterpret
an old request. Preserve strict content verification in both modes.

Keep scientific report contents unchanged. Re-export into a fresh directory
and issue new output hashes; never rename or patch the supplied receipt-bound
file in place. Cosmetic filename changes are not required for this first fix;
the indexed record name is enough to improve the import chooser.

## 2. Compact Browsing Without Losing The Audit Data

Add an explicit, independently specified browsing projection rather than
silently making the full annotated export incomplete:

- Retain the complete current feature projection as the full export, with
  clearer labels. Offer a separately named compact GenBank/EMBL companion;
  its policy and displayed/total counts are visible in its description.
- Show TSS and exon/CDS context first. Coalesce truly identical exon geometry
  for presentation only, retaining all supplied transcript memberships and
  the supplied-versus-unassessed counts. Do not merge clipped exons merely
  because their visible portions coincide.
- Summarize signal intervals per exact source/sample as disjoint feature
  segments where the format/importer supports that faithfully. Never fill
  gaps with one solid outer-span feature. Such a feature means supplied signal
  support, not a called peak, individual read, gene orientation or quantitative
  plateau. Scores and original intervals remain in the full export/report.
- Require an explicit factor/accession selection and bounded retained-peak
  display policy for compact motif annotations. Make hidden counts obvious;
  do not choose all 300 annotations as the default overview, or rank unlike
  matrices against each other. This is display selection of stored values,
  not new peak calling or a biological significance claim.
- Keep full and compact files linked to the same source report and separate
  receipt-bound projection policies. GenBank and EMBL must agree within each
  policy. The full JSON/score tables remain necessary for quantitative traces;
  a flat annotation map is not a lossless replacement for those reports.
- Colour by evidence class with a visible key. Verify any third-party colour
  qualifiers rather than inventing unsupported SnapGene directives. Readability
  must survive a monochrome importer through short labels, grouping and
  feature descriptions alone.

This needs its own small contract review after the naming fix. Do not add an
implicit compact default to ordinary `SaveFile` or change all sequence exports.

### SnapGene Compatibility: Grouping Is Not Fixed Track Placement

Follow-up documentation check: 2026-09-14. The locally installed SnapGene
Viewer app reports version 7.1.1; this does not establish that an import has
passed or that all editing features are enabled by its license.

SnapGene documents a GenBank-compatible extension for feature colour,
directionality and named/coloured segments. It requires specific header
markers and interprets the final feature note as presentation metadata.
Therefore arbitrary colour notes added to ordinary EMBL are not a verified
solution. Use an optional, explicitly marked compatibility export over the
same biological records. See the
[official format specification](https://support.snapgene.com/hc/en-us/articles/10242682237588-What-is-Genbank-SnapGene-Format).

The documented header convention includes `Exported` in the LOCUS name and a
SnapGene-format reference. Account for its map-label/alias conventions when
testing readable import names. State that GENtle generated the file using the
SnapGene GenBank format; do not fabricate a SnapGene producer version, author
or export event. Keep scientific notes separate from the final display note.

Feature types can be shown/hidden as groups; custom types are available from
SnapGene 6.0 and can be exchanged via a feature-type list. The documented
list workflow is not proof that an arbitrary custom type in EMBL will import
automatically. Test that path separately. Sources:
[feature visibility](https://support.snapgene.com/hc/en-us/articles/34537165002004-Changes-to-Features-View-in-SnapGene-8-0),
[custom type exchange](https://support.snapgene.com/hc/en-us/articles/10383958339476-Share-Custom-Feature-Types-via-a-Delimited-List-File).

Proposed adaptation:

- Use concise searchable label prefixes and consistent colours for structure,
  TSS, CUT&RUN, chromatin context and predicted motifs. Preserve exact matrix,
  assay, sample and source identity in details, with a colour key in the
  accompanying overview. Colours describe evidence classes, not confidence.
- In the compact companion, represent each signal source's support footprint
  as a non-directional multi-segment feature, preserving gaps. Retain stepwise
  values in the complete report/export. This groups intervals into one selectable
  object but does not make a quantitative signal graph.
- Where verified, offer custom presentation categories for CUT&RUN signal,
  chromatin signal and predicted TFBS, with an optional type-list companion.
  Preserve ordinary feature semantics in the portable full export. Never call
  a signal feature a CDS or reverse its biological orientation to force its
  position on the map. Custom categories are not evidence classifications.
- Do not promise fixed Y positions, automatic per-assay rows, preserved
  visibility defaults or quantitative track heights from this interchange
  format: no such import contract was found in the checked documentation.
  Type filtering and multi-segment features improve organization, while the
  importer still controls map layout.
- Verify a small synthetic file in the installed Viewer before applying the
  compatibility wrapper to real bundles. Test names, colours, segment gaps,
  directionality and category visibility. Keep plain GenBank/EMBL fallback;
  if grouping is ignored, use the explicit compact subset rather than silently
  dropping annotations from the full file.

This is a scoped extension of the browsing projection, not a request to build
a native `.dna` writer. Stable quantitative lane layout remains the shared
GENtle report/GUI responsibility. Compatibility files need their own bound
projection policy and receipt, just like the compact companion.

## 3. Native GENtle Regulatory Inspection

Connect this to the existing
[DNA-centred integration proposal](gui_regulatory_reporter_integration_plan.md),
not a competing reporter window. Its new native TSS surface remains separately
scoped; the current locus inspector and generic DNA map are not already a
`TssProfileReport` viewer.

- Use the DNA viewer as entry. Default to a linear regulatory overview for a
  positively identified TSS document, without changing ordinary cloning-map
  defaults. Put restriction sites behind an explicit layer toggle; retain an
  easy return to cloning inspection.
- Show a prominent selected TSS marker and aligned genomic, TSS-relative and
  local coordinates. On a minus-strand transcript-oriented window, genomic
  numbers decrease while local bases increase. Use the existing checked
  geometry; never reverse an already oriented sequence again.
- Separate exon models, CUT&RUN factor occupancy, H3K4me3 chromatin context,
  local motif-score traces and imported motif hits into collapsible lanes.
  Group samples/replicates by exact source metadata, not display-name guesses.
  Identical displayed exon models may share a row with a count and an expandable
  transcript list; underlying identities and TSS choices stay separate.
- Render signal as quantitative tracks when the bound report is available,
  not hundreds of equally weighted arrow features. Preserve source scale,
  unavailable states and gaps. A flat file alone permits feature inspection;
  it does not authorize reconstruction of a missing dense score trace.
- Retain detailed hover text and add concise side summaries: genomic limits,
  feature/interval counts, sample/condition, source, score units and missing
  evidence. Show a consistent evidence-class legend even on sparse pages.
- Draw only labels that fit without collisions; show selected features and
  bounded overview summaries first, with all records in an expandable table.
  Display-only visibility and label density never mutate evidence or scores.
- Reuse shared UI intents for entry and source binding. Cache presentation,
  virtualize dense rows, and keep file loading/render preparation off the
  frame loop. No implicit downloads, rescore, or construct creation.

## Verification And Approval

- Name tests: plus/minus genomic TSS, same local offset/different genomic TSS,
  repeated locus with differing reference/strand/window, long or unusual gene
  labels, deterministic collisions and record reordering.
- Export tests: GenBank/EMBL/FASTA retain identical bases; full-mode feature
  locations, clipping, strands, scores and source report remain unchanged.
  The human name uses the TSS coordinate, not either flank boundary.
- Receipt tests: historical policy verifies unchanged exports, the new policy
  verifies new exports, and tampered policy/name/annotations still fail.
- Compact tests: an explicitly counted projection never fills signal gaps,
  merges distinct clipped exons, changes scores or silently drops audit files.
  Imported sparse hits and locally computed scores cannot be conflated.
- GUI tests: label collision bounds, class visibility without data loss,
  separate coincident models, partial transcript coverage, negative strand,
  stale subject rejection and dense first/steady-frame responsiveness.
- Interoperability: public synthetic plus/minus fixtures opened in GENtle and
  the user's older SnapGene; check import names, linear topology, compound
  feature locations and readable labels. Third-party colour/layout support is
  an acceptance observation, not a promise made by the format writer.
- Glen independently re-exports and inspects the actual CD44/TGFB1 bundles on
  the exact resulting revision. Do not commit these supplied files as fixtures
  or change his previously published artifacts without a new receipt.

Recommended first approval: section 1, including legacy receipt compatibility.
Sections 2 and 3 need separate approval; they are not hidden dependencies of
renaming two records. A read-only Claude consultation is available before
finalizing the implementation scope. No implementation is authorized by this
planning document alone.
