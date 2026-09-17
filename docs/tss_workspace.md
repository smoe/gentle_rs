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

Missing flanks are shown as unavailable, never silently clipped. Extend the
anchored parent locus and preview again, or explicitly choose smaller flanks.
Missing transcript annotation produces an explanatory error, not an empty
biological conclusion. Clipped or uncertain 5-prime ends and absent gene linkage
appear as nonselectable diagnostics in `excluded_transcripts`. A fuzzy 3-prime
end alone does not invalidate an exact 5-prime start. A CDS alone
does not establish a TSS. Annotation support is based on imported transcript
features, including GenBank/EMBL features, not an assumption that every genome
annotation format has a transcript index.

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

To discard a stale or legacy registry entry, use
`promoters tss-forget tp73_tss` (`ForgetTssCollection` through `op`). This is a
metadata mutation requiring the usual agent/MCP confirmation, undoable in the
GUI. It does **not** delete member sequences, close windows, or remove lineage.
Retained sequence IDs still prevent overwriting; a fresh derivation therefore
normally needs a new collection ID.

The inner agent is explicitly guided to use this staged route for requests such
as "alle TSS von TP73, jeweils einen in einem Fenster". When no appropriate
locus is loaded, it must first compose existing local genome extraction/extension
commands (or ask permission for retrieval). It must not calculate coordinates,
invent a loop, guess an approval hash, or claim window-opening completion while
the host still reports a queued request.

## Scientific And State Boundaries

- Starts group only when gene identity, annotation source, genomic reference,
  coordinate and strand agree. All source transcript feature memberships remain.
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
- Collections and parent-child lineage persist with the project. Edited or
  missing members cannot silently reuse stale geometry through collection
  navigation; inspect such records individually and derive a new collection.
- The `subject` is a typed `tss_collection` reference accepted by
  `ScanTfbsHitsCollection`. It uses project-sequence map policies but revalidates
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

Glen's live TP73 acceptance remains separate: use the exact candidate binary and
loaded reference, compare inventory membership to source transcripts, repeat
opening without duplicates, inspect both strands and confirm UI responsiveness.
Include a DeltaNp73/P2-only region: truncated TAp73 internal exons must be
excluded, not offered as alternative TSSs. Repeat preview/materialization across
CLI processes and reopen after recomputing features at the same candidate SHA.
Synthetic tests are not a live GUI or real-data scientific sign-off.
