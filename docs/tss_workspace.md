# From A Locus To TSS Windows

GENtle can enumerate exact annotated transcript starts on an anchored project
locus, derive selected windows as a named sequence collection, and open one
native TSS viewer per distinct start. This is annotation interpretation, not
evidence of experimentally established transcription initiation.

## DNA Viewer

1. Open an annotated locus with its genome anchor. Prefer an already prepared
   local reference. Include enough flanking DNA for the requested windows.
2. Choose **TFBS scan > Transcript starts / TSS windows...**, or select
   **Transcript Starts / TSS Windows** in the command palette. The palette uses
   the explicitly selected DNA sequence and reuses an open or loading viewer;
   opening the form does not inspect or materialize anything automatically.
3. Enter the gene symbol or gene ID, a new collection ID, and upstream/downstream
   sizes (defaults: 500/200 bp). Click **Inspect starts (no changes)**.
4. Review the coordinates, strands and transcript memberships. Select the
   desired rows, then **Approve and create selected windows**.
5. Click **Open TSS collection**. Up to 32 member windows open in the background;
   existing or still-loading windows are reused. Larger inventories can be
   materialized as smaller, explicitly selected collections.

For an existing collection, choose **Refresh collections**, then select its row,
or enter its ID directly. The registry browser shows the gene query, source
locus and stored window count across the project. It does not scan or hash
member sequences: **readable / not checked** is not a validation pass. Legacy
and invalid entries stay visible; an unavailable count is not zero. The list is
a snapshot: refresh it after another action changes the project.

Choose **Inspect stored collection** to validate the selected ID.
Validation runs in the background; a valid result lists member coordinates and
transcripts and offers **Copy collection JSON**. A stale/legacy result shows the
engine's diagnostic instead of claiming the members are valid. To remove that
registry entry, choose **Forget registry entry...**, review the named ID and
confirm. Changing the ID cancels this confirmation. Sequences, open windows and
lineage are retained, and the metadata removal is undoable. Normally use a new
collection ID for re-derivation: forgetting does not authorize sequence overwrite.

Missing flanks are shown as unavailable, never silently clipped. Extend the
anchored parent locus and preview again, or explicitly choose smaller flanks.
Missing transcript annotation produces an explanatory error, not an empty
biological conclusion. Clipped or uncertain 5-prime ends appear as nonselectable
diagnostics in `excluded_transcripts`. Annotations without gene linkage instead
appear in `unassigned_transcripts`: these are locus-level diagnostics, not
exclusions from every requested gene and not proven unrelated transcripts. A fuzzy 3-prime
end alone does not invalidate an exact 5-prime start. A CDS alone
does not establish a TSS. Annotation support is based on imported transcript
features, including GenBank/EMBL features, not an assumption that every genome
annotation format has a transcript index.
Both `complement(join(...))` and uniformly reverse-oriented
`join(complement(...),complement(...))` are supported; mixed-strand joins remain
unsupported. Extraction uses the same endpoint interpretation to preserve loss
of a first exon rather than inventing a TSS at a surviving boundary.

Prepared-genome transcript indexing is a separate capability: it currently
supports tabular GTF/GFF, not GenBank/XML. A whole-index request for GenBank/XML
returns `Unsupported`, naming the genome, format and annotation path; a leftover
sidecar is not accepted as a ready index. Region/gene extraction can still
succeed, but explicitly warns that transcript enrichment is unavailable.
This does not prevent TSS derivation from transcript features already imported
into a project sequence. A successfully loaded GTF/GFF index with zero records
is different from an unavailable index; neither proves biological absence.

## Shared Shell And Inner Agent

For an existing anchored `tp73_locus` sequence:

```text
promoters tss-inventory '{"seq_id":"tp73_locus","gene_query":"TP73","collection_id":"tp73_tss","upstream_bp":500,"downstream_bp":200}'
```

The result contains `tss_inventory`, including `approval_sha256`, the effective
`request`, and rows with exact `tss_id` values. The next command must use these
returned values, not guessed identifiers. In the Agent Assistant, click
**Use TSS preview in next prompt**, review the draft and send it. This explicitly
shares the gene, coordinates, transcript IDs and approval digest, not DNA bases;
the prompt is retained with the conversation. Execution receipts alone do not
disclose command results. This action neither sends a request nor approves
materialization; large previews are refused intact rather than truncated.
For example, a materialization
request file contains:

```json
{
  "inventory": {
    "seq_id": "tp73_locus",
    "gene_query": "TP73",
    "collection_id": "tp73_tss",
    "upstream_bp": 500,
    "downstream_bp": 200
  },
  "expected_approval_sha256": "COPY_THE_PREVIEW_DIGEST",
  "selected_tss_ids": ["COPY_AN_EXACT_TSS_ID_FROM_THE_PREVIEW"]
}
```

The placeholders above are explanatory, not executable approvals.

```text
promoters tss-materialize @approved_tss_windows.json
promoters tss-collection tp73_tss
ui open tss-view --collection tp73_tss
ui focus tss-view --collection tp73_tss
ui close tss-view --collection tp73_tss
```

Collection-targeted close closes member windows but retains their sequences.
The original argument-free `ui close tss-view` instead returns the active viewer
to Standard map. Headless UI commands return `applied=false`; only the GUI can
open windows. CLI, MCP and script adapters can use the same `InspectTssInventory`,
`MaterializeTssWindows`, and `GetTssCollection` typed operations through `op`.

For an optional TFBS scan, use the collection by identity rather than copying
member sequence IDs. The engine validates every member before scanning or
writing output:

```text
collections run tfbs-scan --tss-collection tp73_tss --motif AAC --max-hits 100
```

`AAC` is a small illustrative sequence motif, not a TP73 binding model; choose
the intended JASPAR model for a biological analysis. Typed map operations accept
`{"source_kind":"tss_collection","collection_id":"tp73_tss"}` as their
`collection_subject`, retaining that identity in the report.

The same identity-preserving flag also accepts the existing restriction-scan,
digest and primer-specificity operations:

```text
collections run restriction-scan --tss-collection tp73_tss --enzyme EcoRI
collections run digest --tss-collection tp73_tss --enzyme EcoRI --dry-run
```

Digest application still requires the preview's exact plan fingerprint and
explicit `--apply`. Neither example designs or validates primers. For
`collections run primer-specificity --tss-collection ID`, supply the normal
`--pair-rank` or `--pair-index` and `--target-genome`, plus `--member-report`
bindings when reports cannot be resolved uniquely. Primer-report bindings do not
replace collection membership. All four routes reject competing gene-set or
sequence subjects; TFBS/restriction/digest also reject member-sequence overrides.
The engine validates every member before mapping, including before external
specificity work. Unsupported or stale membership is never silently skipped.

To discard a stale or legacy registry entry, use
`promoters tss-forget tp73_tss` (`ForgetTssCollection` through `op`). This is a
metadata mutation requiring the usual agent/MCP confirmation, undoable in the
GUI. It does **not** delete member sequences, close windows, or remove lineage.
Retained sequence IDs still prevent overwriting; a fresh derivation therefore
normally needs a new collection ID.

`promoters tss-list` (`ListTssCollections` through `op`) returns
`gentle.tss_collection_list.v1` under `result.tss_collection_list`. Its rows are
sorted by registry ID, and `validation_status` is always `not_checked`.
`record_status` is `readable`, `legacy`, or `invalid`; it describes stored
metadata only. A corrupt registry is an error, not an empty list. A successfully
listed collection can still fail `promoters tss-collection ID` after a member
edit. Listing never rewrites, migrates or repairs stored evidence.

The inner agent is explicitly guided to use this staged route for requests such
as "alle TSS von TP73, jeweils einen in einem Fenster". When no appropriate
locus is loaded, it must first compose existing local genome extraction/extension
commands (or ask permission for retrieval). It must not calculate coordinates,
invent a loop, guess an approval hash, or claim window-opening completion while
the host still reports a queued request.

## Scientific And State Boundaries

- Starts group only when gene identity, annotation source, genomic reference,
  coordinate and strand agree. All source transcript feature memberships remain.
  A missing `gene_id` can be filled for grouping only when the same gene label
  maps to exactly one explicit ID among this locus's transcript annotations in
  the same source and strand. Conversely, an ID-only transcript participates in
  a symbol search when that ID has one unambiguous label in the same scope.
  Labels differing only in ASCII case count as one label; a deterministic
  spelling is displayed. Resolved IDs, not display labels, determine grouping;
  label-only records use a case-insensitive label key without inventing an ID.
  Conflicting labels are reported, not inherited or chosen as canonical. An
  explicit ID query can still group that ID's records, with the display label
  left unavailable and original annotations retained. Conflicting IDs and
  differing sources/strands remain separate; no gene assignment is inferred
  merely from overlap. Every inferred association is explained in the report.
  A matching label-only transcript with conflicting IDs is warned by transcript, feature, label and sorted IDs; no ID is inferred and it stays label-grouped. Warnings are digest-bound, so affected previews need fresh approval without renaming their TSS rows.
- Local and genomic strands are independent. Each output reads transcript
  5-prime to 3-prime, with TSS at local base `upstream_bp + 1`; negative-strand
  genomic labels therefore decrease along that output.
- The preview binds source DNA/annotations/anchor, effective flanks, memberships
  and output IDs. Source or parameter changes require another preview. All
  selected rows and collisions are checked before any sequences are inserted.
- `snapshot_algorithm=gentle.tss_biological_snapshot.v1` hashes the authoritative
  sequence record, overhangs and genomic anchor identity/coordinates/strand,
  excluding computed RE/ORF/GC caches, timestamps and local resource paths.
  Reloading in another process or recomputing features does not invalidate an
  approval. Legacy cache-sensitive collections remain readable as stored data
  but cannot be reused silently; inspect again under a new collection ID or
  explicitly forget their registry entry.
- Repeating the same approved selection reuses its sequences. A different
  selection requires another collection ID. No records are overwritten.
- Collections created with the previous label-sensitive grouping remain
  inspectable against their saved memberships and member snapshots; inspection
  does not regroup or rename them. New previews use gene-identity grouping and
  therefore may have different TSS IDs and approval hashes. Old approvals cannot
  authorize a changed preview. Review a fresh preview under a new collection
  ID; existing sequences and approvals are not migrated silently.
- Collections and parent-child lineage persist with the project. Edited or
  missing members cannot silently reuse stale geometry through collection
  navigation; inspect such records individually and derive a new collection.
- The `subject` is a typed `tss_collection` reference accepted by
  the supported TFBS, restriction, digest and primer-specificity collection
  operations. It uses project-sequence map policies but revalidates
  persisted membership, sequence/annotation snapshots and anchors first. Raw
  `--seq-ids` scans remain ordinary sequence scans without TSS freshness claims.
  Its read-only scan is optional. No motif, occupancy,
  reporter or cross-source consensus analysis is implied by derivation.
- Existing annotations are retained, with supporting transcript exons projected
  into the window. Metadata says **project annotation derivation**, not verified
  external report bundle. No new reference authentication is claimed.

Inventory limits: 20 Mb source locus, 100,000 features, 10,000 matching transcript
features, 256 distinct starts, 2 Mb per window and 32 Mb total candidate windows.
The inventory is exhaustive only for the transcript annotations in that loaded
source, not for all databases or every possible biological start.
Materialization also limits features to 100,000 per window and 250,000 total.
Genome-derived transcript bounds are checked against their projected 5-prime
end. Linear extraction records when the original start is removed, including
when an entire first exon disappears. Historical clipped records with neither
genomic bounds nor truncation markers cannot be reconstructed reliably; reload
them from authoritative annotation before making TSS claims.

## Acceptance

See the [offline synthetic collection tutorial](tutorial/08-15_tss_collection_gui.md)
for Glen's clean-profile walkthrough and typed repeated-opening, cancellation,
Forget/Undo checkpoints. Save/restart, visual orientation and deliberately
stale-member GUI checks remain separate manual acceptance. Neither authoring
the contract nor running its unit tests establishes a live GUI pass.

Deterministic synthetic tests exercise shared/opposite-strand starts, overlapping
genes, missing/partial annotation, flank refusal, all four local/genomic strand
combinations, stale/collision refusal, persistence/reuse, optional TFBS collection
scanning, shell parsing and deferred/reused GUI window requests.
Prepared-index regressions separately cover GenBank/XML refusal with a cached
sidecar, valid empty GTF/GFF inventories, missing/corrupt index diagnostics,
on-demand rebuilding, optional extraction warnings and promoter-background
error propagation. These are shared-engine tests, not desktop acceptance.
Prepared GenBank input is previewed, materialized and reopened in separate test
processes; regressions cover cache refresh, metadata-only removal/undo, edited
member scan refusal, fuzzy 3-prime ends and clipped first exons on both strands.
Follow-up regressions cover both reverse-join encodings and extraction, ambiguous
gene aliases, locus-level unassigned diagnostics, all four scan/map flags, and
background GUI collection inspection/confirmed metadata-only removal. Older
stored reports remain readable: an absent `unassigned_transcripts` defaults to
empty and is omitted on serialization. New previews whose grouping or diagnostics
change need a new approval; existing stored collections are not rewritten.

Glen's live TP73 acceptance remains separate: use the exact candidate binary and
loaded reference, compare inventory membership to source transcripts, repeat
opening without duplicates, inspect both strands and confirm UI responsiveness.
Include a DeltaNp73/P2-only region: truncated TAp73 internal exons must be
excluded, not offered as alternative TSSs. Repeat preview/materialization across
CLI processes and reopen after recomputing features at the same candidate SHA.
Synthetic tests are not a live GUI or real-data scientific sign-off.
