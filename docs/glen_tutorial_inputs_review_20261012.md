# Glen Tutorial Inputs: Independent Archive Inspection

Inspected on 2026-10-12 on macOS/arm64, against local source
`4ad6e5b6` on `gentle_rs_2_main`. This is an offline data/receipt inspection,
not a GENtle replay, native screenshot review or scientific approval.
Glen's executed evidence is bound to the older source
`a403d2eba5c330884121dfcd4ee3956831f36d36`; it does not accept newer tutorials.

## Archive Identity And Integrity

The owner supplied
`/Users/u005069/Downloads/gentle-tutorial-inputs-a403d2eb-20261012.tar.gz`.
Its SHA-256 is
`99d333172591a39eeed0ceb29073d5aff6be44893b8e06a33946db8677a54206`.
The original archive is retained locally, not committed or publicly published.
An external reviewer must obtain these exact bytes through the owner.

Before extraction, all 114 members were checked for bounded size, duplicate or
unsafe paths, links and special files. Only regular files/directories were
accepted; extraction used Python's `data` filter into a fresh temporary parent.
The archive expands to 19,173,351 file bytes. All 89 entries in its root
`SHA256SUMS` match, including nested manifests and raw payloads. No archived
script, SQL, binary or GENtle command was executed.

To recheck integrity after independently checking/extracting the archive:

```sh
shasum -a 256 /PATH/gentle-tutorial-inputs-a403d2eb-20261012.tar.gz
```

From the extracted root:

```sh
shasum -a 256 -c SHA256SUMS
```

The archived environment record names Linux, Rust/Cargo 1.96.1, default
features, locked Cargo `dev`, DuckDB 1.5.5 and no `DISPLAY`. Executables are
not included, so their stated hashes are Glen's retained identity evidence,
not independently recomputed binary hashes. The lockfile hash is
`ead4421545a125267e33af84aeeea39f0e7ea2f60717a59ef0271d0e71855f1b`.

## Independently Checked Contents

### A1: Separate Synthetic TSS Companion

This is the original `fixture-genome` / `fixture-assembly` companion with
MA0861.2 and MA1961.2, **not** tutorial 08.16's newer `tss-view` generator
preset (`synthetic-genome-v1`, MA0004.1). It must not replace that report or
its receipts. The input and interval scores are explicitly synthetic.

- Full query: six hits, `query_complete=true`, `truncated=false`, raw scores
  `-4, -3.5, 1.25, -0.5, 0, 0` and both genomic strands. Complete means
  complete for package retention, not a complete genome scan.
- Row-limited query: two hits, `query_complete=false`, `truncated=true`.
  Its per-matrix coverage is explicitly truncated, not absence of binding.
- Query-declared package file sizes/hashes match the archived files. The
  enriched reports contain the exact parsed saved query reports. Raw query
  digests and the engine's embedded report digests are distinct roles; no
  Python JSON reserialization is substituted for the engine's canonical hash.
- The base/full/partial report `windows` are identical, including score arrays,
  geometry and sequence bindings. All 33 declared output hashes across their
  three receipts match (9 base, 12 full, 12 partial).
- Both enriched GenBank and EMBL sequences equal the input FASTA's 1,001 bases
  and normalized-sequence digest. A `[295,311)` 16-base footprint intersects
  the displayed `[300,1301)` window in 11 bases; the original footprint stays
  available in the query. This is coordinate checking, not native paint proof.
- The context retains an `available` treatment with two synthetic intervals
  and a `not_prepared` control without intervals. Imported
  `log2_relative_risk` scores remain separate from local matrix scores.

### A2: Genuine Retained FLNA Snapshot

The complete public project is not a synthetic replacement. Each embedded
raw Ensembl/UniProt payload matches its separately extracted file byte-for-byte:

| Payload | SHA-256 |
| --- | --- |
| Ensembl lookup JSON | `ca26c5e89911b2c612a8097896ece5b1c556f68876ad774e961b305edce6dfce` |
| Ensembl sequence JSON | `5b9cf22ecc1993213f4b78114676a681745dc099169ba1eda80e94bdbf247148` |
| UniProt P21333 flat file | `a0a398defb562e270eb77da541701c17d491e0161c4555bb82f1ef9db3a46c02` |
| Complete project | `88b40cf12d73d7ad1ca88a8d0fc59f9931abcf812bc0b869f81c47993d057452` |
| Replayed splicing view JSON | `cc1bc73022c074f90fce7afcb88c181cb982ef3f9d2e3253fb3ccca951e5c0d7` |

The metadata records retrieval on 2026-10-08 at 00:14:31.198 UTC for Ensembl,
gene `ENSG00000196924.21`, GRCh38 chromosome X and a genomic-minus anchor.
The 30,120-base local sequence equals the stored provider sequence and has
SHA-256 `928d71a24b4de0f0b7f0b7114eef4f6c6550417c8d7d2eefc972564a39c9415c`.
The splicing view is local plus and retains 36 transcripts, 103 unique exons,
73 junctions and six event summaries.

The stable selector `transcript_versioned_id=ENST00000369850.10` resolves to
exactly one mRNA in this project (current index 300), matching the view target.
Index 300 remains an observed result, not a reusable selector. Both reviewed
joins agree with the raw UniProt Ensembl xref lines and the view's matched xrefs:

- `ENST00000360319.9` / `ENSP00000353467.4` / `P21333-2`.
- `ENST00000369850.10` / `ENSP00000358866.3` / `P21333-1`.

Provider URLs, timestamps and UniProt's CC BY 4.0 notice are retained. This
inspection does not supply a new redistribution/legal verdict.
**No Ensembl annotation-release record or contemporaneous `/info/data`
payload is present.** Gene/transcript versions and GRCh38 are not a release
number. The snapshot is dated and hash-bound, not release-pinned. Do not fetch
today's release label and attribute it retrospectively. A teaching exon and its
observed consequence, native normal-size/2400-px readability and human review
are also missing. The [FLNA draft](tutorial/drafts/flna_splicing_interpretation.md)
therefore stays outside the active catalog.

### A3: Repeat Selection And Explicit Limits

The saved query contains four distinct UCSC excerpt loci, three genomic minus.
The incomplete-cache report preserves repeat-midpoint geometry while warning
that transcript projection is unavailable; other geometry modes are explicitly
unavailable. The prepared-reference attempt records `timeout_exit=124`, not a
completed transcript-aware cohort. No new performance reproduction is claimed.
RNA ranking remains unimplemented. This does not execute the newer LocalProject
two-locus regression or establish real repeat biology.

## Follow-Up For Glen, Through Steffen

Please use the final committed source containing the tutorial repair, record
its full SHA and rebuild/reuse only matching binaries. No push, dispatch, tag,
publication or editing private studies is requested.

1. Run current-source `gentle_examples_docs tutorial-check`, then
   `python3 scripts/check_tutorial_checkouts.py --binary /EXACT/gentle_examples_docs`.
   Keep review dates null: registration fixes missing-entry warnings, not
   acceptance. Report whether committed `generated/report.json` now matches.
2. Execute the new 08.16 `tss-view` companion and the LocalProject repeat test
   named in [the original handoff](glen_tutorial_gap_preparation_20261012.md).
   This archive's separate A1/A3 examples do not substitute for those tests.
3. Complete C/D and enriched-TSS native learner checks/screenshots; retain raw
   captures and separate derived teaching images. All are currently not run.
4. For FLNA, search retained acquisition logs for the historical release record.
   If absent, report that explicitly and propose either an honestly dated REST
   snapshot lesson for owner/Claude review or a newly pinned public bundle with
   fresh provenance. Obtain human teaching-exon selection before claiming an
   expected consequence; then run the native normal-size/2400-px checks.

Reply through Steffen with exact source/binary/input hashes, commands,
per-item pass/fail/not-run/blocked, retained artifact paths and missing platforms.
Keep your executed Linux evidence separate from this offline inspection and
from future candidate/native/scientific acceptance.
