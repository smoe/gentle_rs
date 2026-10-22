# 10.04: What Did Our Targeted cDNA Primers Recover?

Status: manual/hybrid; implementation and Rust/native acceptance pending Glen/CI.
This is a first working lesson for us to improve together, not a validated assay.

## The Question

The [TapeStation lesson](10-03_measured_gel_sizing_gui_cli.md) shows fragment
sizes, not which genes or isoforms produced those fragments. This lesson starts
with basecalled FASTQ to ask which annotated transcript sequences merit closer
inspection. Seed nomination is not proof that a transcript explains the entire
read. POD5 and the complete raw-signal download are not required here.

We will keep the custom-only libraries separate:

| Sample | Sequencing barcode | Condition | TapeStation after-cleanup lane |
| --- | --- | --- | --- |
| 26R266-1 | barcode13 | custom primers/adapter only | C1 |
| 26R266-2 | barcode14 | custom primers/adapter only | D1 |

These barcode assignments were supplied by the owner on 2026-10-22, not inferred
from filenames. TapeStation lane IDs change between reports. Folder names may
change; keep this mapping and the supplier checksum manifest. The mixed and
manufacturer-only libraries are not part of this initial screen.

The lookup is **one combined word index**, not a separate read scan for each
gene. A packed sequence word points to every matching transcript and position;
the transcript then supplies its locus and exon/junction geometry. Namespaced
member IDs prevent two loci's local feature numbers or transcript labels from
colliding. SHA-256 hashes bind the panel and sequences; they are provenance
identities, not scores or evidence of gene identity. This also keeps a panel of
about twenty enriched targets inspectable without forcing one winner per read.

## Follow This Through The Inner Agent

Every scientific/rehearsal stage below has a GENtle command. The optional Python
helper in step 3 only prepares JSON for terminal users; the inner agent drafts
the same `ScreenRnaReadPanelRequest` inline and invokes `panel-screen` directly.
No operating-system shell access or guessed file paths are needed.

Start in a disposable project and tell the Agent Assistant:

> Walk me through tutorial 10.04 one reviewed step at a time. Prepare the
> synthetic demo in NEW_ABSOLUTE_DIRECTORY. Show the actual commands and await
> approval for writes and project loading. Keep ambiguity visible. Do not use
> my sequencing downloads until I explicitly name a completed file.

| Task | Shared route the agent should propose |
| --- | --- |
| Rehearsal files | `rna-reads panel-demo OUTPUT_DIR` |
| Disposable starter | `load-project OUTPUT_DIR/artifacts/starter.json` (replaces active project) |
| Discover actual local references | `genomes list`, then `genomes status GENOME_ID` |
| Load an actual target | existing `genomes extract-gene` with the explicitly chosen reference/gene |
| Inspect transcript features | `features query SEQ_ID --kind mRNA --include-qualifiers` |
| Build or rebuild the panel | `rna-reads panel-build BUILD_REQUEST_JSON_OR_@FILE` |
| Screen explicit completed files | `rna-reads panel-screen SCREEN_REQUEST_JSON_OR_@FILE` |
| Inspect panel, report and selected nomination rows | `rna-reads panel-show PATH [--file-id ID] [--rows N]` |
| Export sample/transcript support and class-bound co-nomination sets | `rna-reads panel-export-summary SUMMARY_REQUEST_JSON_OR_@FILE` |
| Reopen the saved count summary | `rna-reads panel-show SUMMARY_REPORT_PATH` |
| Explicitly verify the three summary tables | `rna-reads panel-show SUMMARY_REPORT_PATH --verify-artifacts` |
| Request cooperative cancellation | `rna-reads panel-cancel SCREEN_REQUEST_JSON_OR_@FILE` |
| See active work | `introspect runtime` (observation, not a completion receipt) |
| Inspect nominated reads further | existing `rna-reads interpret`, `align-report`, `inspect-alignments` |

For the first real chunk, provide its exact current path, checksum-manifest
path, panel path, a fresh output/cancel path and the sample/barcode/priming
mapping. Explicitly confirm copying that file has finished. Ask the agent to
draft the full request, inspect it, and only then approve screening. Adding a
file or changing thresholds is a new request, not permission inherited from a
previous run. Model/live GUI interaction is not accepted merely by testing
these deterministic commands; Glen/CI must verify those boundaries separately.

## 1. Rehearse Without Your Sequencing Data

Use an externally verified, source-matched GENtle CLI containing the panel
commands. Do not build Rust for this tutorial on a host where builds are
delegated to Glen/CI. Replace `GENTLE` below with its absolute executable path.
Run terminal commands from the repository root, not inside the GUI Shell.

```bash
GENTLE="/absolute/path/to/gentle_cli"
RUN="$(mktemp -d)/rna-panel-demo"
STATE="$RUN/project.json"
"$GENTLE" --state "$STATE" shell "rna-reads panel-demo $RUN"
"$GENTLE" --state "$STATE" shell "load-project $RUN/artifacts/starter.json"
"$GENTLE" --state "$STATE" shell "rna-reads panel-build @$RUN/artifacts/build.json"
"$GENTLE" --state "$STATE" --progress-stderr shell "rna-reads panel-screen @$RUN/artifacts/screen.json"
"$GENTLE" --state "$STATE" shell "rna-reads panel-show $RUN/screen/report.json --rows 4"
"$GENTLE" --state "$STATE" shell "rna-reads panel-export-summary @$RUN/artifacts/summary.json"
"$GENTLE" --state "$STATE" shell "rna-reads panel-show $RUN/summary/report.json"
"$GENTLE" --state "$STATE" shell "rna-reads panel-show $RUN/summary/report.json --verify-artifacts"
```

Use a temporary path without spaces for the literal shared-shell strings above;
for paths containing spaces/apostrophes, invoke `op @FILE` directly, or use
GENtle's shared-shell quoting helper rather than concatenating unquoted paths.
The demo operation writes only artificial references, tiny compressed reads
and requests. It never accesses your downloaded data. `load-project` replaces
the active project: use only the disposable starter, never your working project.

**Checkpoint:** open `$RUN/screen/report.json`. Require `status: completed`, two
verified file receipts, four reads per file and no `files_not_completed`.
The report names the exact MD5 and SHA-256 of the compressed input bytes. The
gzip files were not unpacked onto disk or read twice for an upfront checksum.

The `demo_0` read contains only a shared first exon. Its TSV row should retain
both `locus_0:0` and `locus_0:1`, despite the two transcripts using different
terminal exons. **Also inspect `demo_1`:** it is the entire artificial A
transcript and crosses A's junction at offset 40. A has
`supported_junction_offsets: [40]`; B has `[]` but remains nominated by their
shared first exon. This screen does not reject B based on the rest of the read.
Do not call these whole-read-compatible transcripts or isoform-resolved counts.
The repeated transcript labels across different toy loci must
not merge their identities. `demo_3` has only N bases and is `unmatched`.
These are software checks on artificial DNA, not results for PATZ1.

The demo gives each artificial library a distinct read-header description so
its two compressed inputs are not identical copies. The summary must retain
two samples, four verified reads per sample and exact class-bound co-nomination sets.
It must not silently distribute the shared first-exon read across isoforms.

## 2. Prepare The Five Actual Loci

Use PATZ1, TP53, TP63, TP73 and E2F1 from one explicitly chosen, already prepared
reference release. Do not concatenate their genomic sequences or reuse feature
IDs from another loaded record. Export/save the source project for reproducibility.
Existing local reference extraction can be used, for example:

```bash
"$GENTLE" --state "$STATE" genomes extract-gene "Human GRCh38 Ensembl 116" PATZ1 \
  --output-id locus_patz1 --extract-mode gene --annotation-scope core
```

Repeat the existing command explicitly for TP53, TP63, TP73 and E2F1 with fresh
output IDs. First inspect `genomes status "Human GRCh38 Ensembl 116"` for your reference/cache; extraction
is not permission to download or prepare a missing reference. Existing loaded
loci are equally suitable if their provenance and annotation are known.

Inspect each record's actual transcript/gene feature IDs in GENtle, for example
`features query locus_patz1 --kind mRNA --include-qualifiers`. Construct
a build request like the following, adding one source row for each actual locus:

```json
{
  "schema": "gentle.rna_read_panel_build_request.v1",
  "output_path": "/absolute/run/reference-panel.json",
  "sources": [{
    "locus_id": "patz1", "seq_id": "locus_patz1", "source_feature_id": 0,
    "source": "Ensembl", "assembly": "GRCh38", "annotation_release": "116",
    "gene_id": "REPLACE_WITH_EXACT_SOURCE_GENE_ID", "gene_symbol": "PATZ1",
    "selected_transcript_ids": []
  }]
}
```

**Do not run this literal placeholder.** Replace feature 0, the gene ID, assembly
and release with the inspected source values. An empty transcript list includes
all transcripts in that explicit source feature's gene/strand group. A nonempty
list restricts it to those exact transcript IDs; it is not a canonical-isoform
selection. Selected transcripts, source loci and all annotated transcripts are
reported separately. Check these counts before screening reads.

```bash
"$GENTLE" --state "$STATE" shell 'rna-reads panel-build @/absolute/run/build.json'
```

**Checkpoint:** review the panel's source hashes, selected member IDs, exon
intervals, transcript-oriented junction offsets and independent strand per
locus. Source labels are supplied provenance, not authenticated annotation.
Potential overlapping source views are disclosed; screening refuses them until
you choose/curate non-duplicated targets. Do not load Ensembl and RefSeq versions
of the same locus as two independent genes. There is no automatic consensus.

## 3. Start With One Completed FASTQ Chunk

Set `DATA` to your current download directory and `ANALYSIS` to a new directory
outside the download. In the owner's current layout DATA is
`$HOME/GitHub/gentle_rs/data/MelanomaEnrichedSequencing`; adjust it after moving
the folders. Set `FIRST` to one exact, finished `.fastq.gz` file in barcode13.
Do not select a file merely because its size has stopped changing briefly.

```bash
DATA="$HOME/GitHub/gentle_rs/data/MelanomaEnrichedSequencing"
ANALYSIS="$(mktemp -d)"
PANEL="/absolute/run/reference-panel.json"
FIRST="/absolute/path/to/a/completed/barcode13/chunk.fastq.gz"
MANIFEST="/absolute/path/to/the/supplied/MD5SUMS"
python3 scripts/prepare_rna_panel_tutorial.py request \
  --panel "$PANEL" --md5-manifest "$MANIFEST" \
  --output-dir "$ANALYSIS/first-chunk" --output "$ANALYSIS/request.json" \
  --input 26R266-1 barcode13 custom_only "$FIRST" --declared-complete \
  --write-locus-fastas --cancel-path "$ANALYSIS/CANCEL"
"$GENTLE" --state "$STATE" --progress-stderr op \
  "{\"ScreenRnaReadPanel\":{\"request\":$(cat "$ANALYSIS/request.json")}}"
```

The input flag declares the sample, barcode and priming explicitly. The helper
reads no sequencing bytes; it only writes JSON. `declared_complete` is your copy
completion assertion. GENtle then parses gzip directly, computes SHA-256 and
MD5 during that same input pass, and checks the supplier MD5 at EOF. A corrupt
gzip, malformed FASTQ, changing file or checksum mismatch is a hard error, not
an empty biological result. No result directory is published for such a failure.
That also discards earlier unpublished units from the same failed run. Start
with one chunk per run while transfers and settings are still being checked;
do not risk an entire library on an untested multi-file request.

The supplier manifest's original absolute server paths need not exist locally.
First locate the actual supplier manifest: it may be named `MD5SUMS`,
`md5_checksums.txt` or `*.md5`; there is no required filename. Choose it explicitly,
not a wildcard selecting FASTQ files from a still-growing download. Inspect its
entries before drafting the request. Do not manufacture expected digests by
hashing an unfinished local copy, or assume the demo's manifest applies to real reads.
A unique unchanged filename is sufficient; duplicate filenames require the
exact supplier `checksum_manifest_key`. Renaming folders is fine. Renaming files
requires that explicit key or an explicit expected digest. MD5 checks transfer
integrity; it does not authenticate the source. SHA-256 binds what GENtle read.

**Progress checkpoint:** `--progress-stderr` and `introspect runtime` expose
file ordinal, stable file/sample IDs, processed reads, and compressed bytes
read/total. Compressed bytes include decoder read-ahead, not the proportion of
reads already classified; reaching the file size is not checksum success.
`streaming_unverified` remains unverified until EOF and matching checksums.
`file_checksum_verified` marks a verified unit, not completion of the whole run;
the final report remains the completion/publication oracle. Progress and marker
polls are throttled to approximately 100 ms, not emitted per sequence word.

To request cooperative cancellation from another terminal, create the exact
declared cancel marker through GENtle. Only completed, checksum-verified file units are retained;
the current partial file is discarded and the report says `cancelled`.

```bash
"$GENTLE" --state "$STATE" shell "rna-reads panel-cancel @$ANALYSIS/request.json"
```

Markers are bound to the exact original or effective screen request. Matching
pre-start cancellation is supported; a stale, malformed or different-request
marker is refused, never silently treated as cancellation. Do not create markers
with `touch` or reuse an old marker for changed inputs/settings.
Use a fresh cancel path for a new run. A command returning with a
cancelled report is not a completed analysis. New runs always need a fresh
output directory; no report or input is overwritten. Start with the terminal
CLI so this long job does not occupy an interactive GENtle window. GUI Shell,
MCP `op` (`confirm=true`), Python `op`, JS/Lua and the inner agent reach the same
engine operations; external writes require explicit approval. No dedicated
panel GUI or separate job scheduler is added here. Standard four-line FASTQ
records are supported (LF/CRLF, Phred+33, concatenated gzip members); wrapped
sequence/quality records are refused, not silently repaired.

## 4. Read The Evidence, Not Just The Largest Count

Each verified file has a nomination TSV under `artifacts/file_NNNN/` and its
own checksum receipt. The TSV retains original read IDs plus file/record identity
and **decompressed** record offsets. Those offsets are not seek positions in
the compressed file. The TSV also retains the full expected input digest,
verified before publication. Candidate JSON retains every qualifying template and
orientation, its ordered seed-chain coverage, original-read interval and any
junction words with sufficient bases on both sides.

Use `rna-reads panel-show /absolute/run/report.json --file-id file_0000 --rows 20`
to inspect the same structured rows through the Shell or an agent. Metadata
content binding is checked; this bounded prefix inspection deliberately reports
`artifact_content_verified=false`, because it does not re-hash a huge TSV.
Retained artifact receipts describe what screening wrote, not immutable storage.

| Class | What it means here |
| --- | --- |
| `single_locus_candidate` | qualifying chain support for one panel locus; not a unique isoform |
| `multi_locus_ambiguous` | more than one locus remains supported; no winner is selected |
| `multi_segment_candidate` | disjoint read segments support different loci; not proof of a fusion/chimera |
| `partial` | some indexed words match, but no chain passes the declared threshold |
| `unmatched` | no indexed canonical word matched; not proof the transcript is absent |

In a v2 screen report, `equivalence_classes` keys are canonical JSON objects
with both `class` and sorted `member_ids`. A read increments that joint key once.
The same members in overlapping cross-locus support and disjoint-segment support
remain different keys. A v1 report stored only member IDs and cannot recover
that association; an exported row is explicitly `legacy_class_unresolved`.
Do not add the same ambiguous read to every isoform's exclusive count. A read
missing the promoter-specific 5' region may still co-nominate TA and DeltaN
forms. A terminal-exon primer site alone does not recover the missing end.
There is no oligo(dT) flip, automatic primer trimming, expression estimate,
whole-transcriptome uniqueness claim or calibrated confidence score.

### 4.1 Export Supporting Read Counts, Not Molecule Copies

`panel-export-summary` turns one or more saved screen reports into tables you
can inspect in a spreadsheet. It reads their retained `artifacts/panel.json`,
not FASTQ, POD5 or the potentially huge nomination TSV. Keep each run's report
and retained panel together, even if you move the run directory.

For the rehearsal, the commands above use the generated `artifacts/summary.json`.
For actual runs, ask the inner agent to draft this request with the exact paths,
show it for approval and invoke `rna-reads panel-export-summary`:

```json
{
  "schema": "gentle.rna_read_panel_summary_request.v1",
  "report_paths": [
    "/absolute/run-one/report.json",
    "/absolute/run-two/report.json"
  ],
  "output_dir": "/absolute/new-summary"
}
```

| Saved output | How to interpret it |
| --- | --- |
| `artifacts/samples.tsv` | One row per requested sample: barcode, priming, requested/verified files and `status`. Assessed rows carry reads, `single_candidate_reads`, `co_nominated_reads`, and `no_candidate_reads`; these last three sum to reads. An unassessed sample has blank read counts, not zero. |
| `artifacts/transcripts.tsv` | One row per selected panel member and assessed sample, including observed zero support. Gene/locus/member IDs keep reused transcript labels distinct. `single_candidate_reads` nominated that member alone; `co_nominated_reads` nominated it together with alternatives through partial seed support, not whole-read fit. |
| `artifacts/co_nomination_sets.tsv` | One row per sample, `class` and exact sorted member set, with read count. A read co-nominating A and B increments `[A,B]` in that class once, never either singleton. `member_ids_json` joins to `member_id` in the transcript table. |
| `report.json` | The same structured counts, panel/policy binding, source report/input hashes, missing-file ledger, table hashes and non-claims. This is the final publication marker. Reopen it with `panel-show`. |

For example, **hypothetical, not a demo or biological result**: 12 reads support
only A, 7 only B, 30 co-nominate A and B through partial seeds and 5 nominate
neither. The sample has 54 reads, including 19 single-candidate and 30
co-nominated reads. Transcript A shows 12 single-candidate and 30 co-nominated
reads; B shows 7 and 30. If all 30 have one class, the co-nomination table
contains `[A,B]: 30` for that class; different classes get separate rows.
Do **not** turn these into 42
copies of A plus 37 copies of B, or automatically split the 30 between them.

Single-candidate means only one member in this particular panel passed its
uncalibrated seed-chain policy. It does not establish an exclusive isoform;
inspect distinguishing aligned junctions separately. These are **supporting
read counts**, not original molecule copy numbers or expression. Enrichment,
amplification, incomplete capture and repeated sequencing can distort abundance.
No PCR/UMI or read-level deduplication is performed.
The summary does not add junction-resolution counts: member-specific junction
evidence remains in each read's candidate JSON in the nomination TSV.

Combine only reports from the exact same panel and full nomination policy.
Repeated report hashes or compressed-input SHA-256 values are refused, even
under renamed paths or different sample labels: select each file once. Different
chunks may combine within one sample only with the same barcode/priming binding.
Distinct compressed files can still contain duplicate reads; the refusal is
not molecule-level deduplication. **Never include a merged delivery FASTQ and
its constituent chunks in one analysis:** their different compressed hashes
do not detect that overlap, and their reads would be counted twice.
Never pool the two actual sample identities.

A cancelled source keeps `completeness: partial`, even if another source
completed. Only its verified completed files contribute counts. Samples with
no completed files have no numerical rows in the transcript table; their
`samples.tsv` row is `unassessed` with blank read counts. JSON `sample_coverage`
lists each requested sample's barcode/priming, requested/verified file counts,
`complete|partial|unassessed` status and named `inputs_not_completed`, not just
anonymous paths. `source_runs` retains each missing input identity as well.
These are selected request-file units; including an old cancelled run keeps
its missing ledger/partial status even if a later run covers a remaining file.
Inspect and select evidence deliberately; there is no automatic reconciliation.
Do not interpret unassessed input as zero expression. `complete`
means complete for the explicit report selection, not all facility files or
all molecules. Export checks content bindings and count consistency but does
not reverify original inputs or sidecars; both revalidation flags remain false.
Legacy v1 summaries remain inspectable with their historical field names;
`ambiguous_compatible_reads` was only partial seed co-nomination despite that
name. V2 writes `co_nominated_reads` and class-bound rows. When v1 screen counts
are exported, `class_member_counts_complete=false` and the class is
`legacy_class_unresolved`; class totals cannot be joined to member sets retroactively.
An I/O failure can leave artifacts without `report.json`; those are not a
published summary. Existing output directories are never overwritten.

Limits fail explicitly rather than silently truncate: 128 reports, 256 MB of
combined report/panel metadata, 4096 completed files, 128 observed samples,
250,000 transcript rows, 65,536 co-nomination sets, two million member visits and
128 MB of summary JSON. No large-data throughput claim follows from these caps.

**Verify the exported tables explicitly:** `panel-show SUMMARY_REPORT_PATH
--verify-artifacts` reads only the three allowlisted summary TSVs, checks their
exact sizes/SHA-256 and returns `rna_read_panel_verification` bound to this
summary digest. Total table size is capped at 256 MB; traversal, symlinks and
unknown table paths are refused. The default show remains metadata-only.
Changing LF to CRLF in a TSV changes its bytes and correctly fails verification.
This receipt does not alter the saved summary, reverify FASTQ/nomination rows,
authenticate source evidence or reassess the biology. Nomination-prefix
inspection still reports `artifact_content_verified=false`.

`k=15`, minimum chain coverage 45 bases and junction flanks 5 bases are explicit
**uncalibrated starting settings**, not experimentally validated thresholds.
Index/anchor/equivalence budgets fail explicitly rather than silently discarding
later isoforms. Short reads, sequencing errors, repeats and untrimmed adapters
can change nomination rates. We will tune only after reviewing real reads and
representative background together, retaining the old policy/report.
If an anchor/chain or class-key budget is exceeded, the run still fails closed,
with file ID, zero-based record index and read ID in the diagnostic. No
`resource_limited` continuation or result directory is silently published.
Retain the diagnostic and adjust a new explicit request only after review;
previous units in that failed unpublished run are not summary evidence.

## 5. Inspect Alignments As A Separate Follow-Up

When `write_locus_fastas` is true, each file receipt lists the exact `locus_id`
for its nominated FASTA. Only nominated reads are exported, not a full FASTQ
conversion. A cross-locus read is retained in every nominated locus's FASTA.
Original reads remain untouched. Retain the nomination TSV alongside these files.
Exported bases use the screen's uppercase/U-to-T normalization; no bases are
trimmed, and quality values are validated but not used as a filtering score.
FASTA header descriptions retain `panel_class`, `ambiguous` and `panel_members`;
the original read identifier is preserved after the stable file/record token.
That token includes the sample and full verified input digest, preventing
different chunks in separate runs from acquiring the same synthetic read name.

Use the existing single-locus path on a selected nominated FASTA, with the actual
source sequence and feature IDs. Substitute the placeholders before running:

```bash
"$GENTLE" --state "$STATE" rna-reads interpret ACTUAL_SEQ_ID ACTUAL_FEATURE_ID \
  /absolute/path/to/nominated.fasta --report-id followup_one_chunk \
  --no-cdna-poly-t-flip --scope target_group_any_strand
"$GENTLE" --state "$STATE" rna-reads align-report followup_one_chunk --selection all
```

This reuses existing alignment/inspection, not a new panel aligner. Its legacy
best-transcript row **must not override** an unresolved panel class: annotation-
order ties and incomplete coverage remain limitations of that legacy path.
Review actual distinguishing junctions, aligned sequence and alternative mappings
before claiming an isoform. Full-length flags are not proof of complete capture
under this custom-priming protocol. Quality-dependent and background validation
remain separate; this screen is not an oligo-order or biological approval.

## 6. Add Chunks And Targets Deliberately

Add another `--input 26R266-2 barcode14 custom_only /exact/completed/file.fastq.gz`
to the request helper when the second library is ready. Multiple inputs may
belong to one sample; never pool the two sample identities. Select explicit
completed files, not a wildcard snapshot of a changing download directory.
The first download observed here contains `fastq_pass` chunks: conclusions from
those chunks are conditional on that delivered subset, not all molecules/run reads.

The same panel shape supports approximately twenty enriched targets without
twenty scans of each FASTQ. Distinguish genes, selected transcript targets and
all annotated isoforms; do not assume twenty genes means twenty transcripts.
Two primer targets in PATZ1 still belong to one locus; the number of enrichment
primers, genes and annotation isoforms need not be equal.
Preserve the default complete transcript group unless you deliberately restrict
it. There is no measured large-data throughput claim yet.

If cancelled, keep the verified completed-file ledger and prepare a new request
containing only the remaining files, in a new output directory. This is explicit
file-boundary continuation, not an automatic in-file/gzip checkpoint. Do not
silently count a completed file twice when combining run summaries.
The summary exporter enforces distinct compressed-input hashes; retain the
partial flag and missing-file ledger when including a cancelled source.

## Notes For Our First Joint Review

Please retain the source/binary version, panel, request, supplier manifest and
report. For the first chunk, we will inspect: integrity/completion, target counts,
the class-bound co-nomination sets, junction support, read-length/coverage examples and
nomination rates. Then record which step was unclear and what extra table or
plot would help. Neither a maximum score nor a TapeStation band proves identity.

## Verification Scope

Synthetic Rust tests are authored but not run locally; Glen/CI owns execution,
cross-platform paths/gzip, canonical parity/catalog freshness and native checks.
The demo uses literal fixed-LCG DNA and generated gzip/checksum files; no private
sequencing reads are committed. Python source/metadata checks do not establish Rust
acceptance, real-library sensitivity/specificity, throughput or scientific validity.
