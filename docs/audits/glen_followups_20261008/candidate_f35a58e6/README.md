# Glen Final-Candidate Evidence At `f35a58e6`

## Identity And Preservation

Glen executed the consolidated handoff on a clean detached checkout of
`f35a58e6d2bdd7e0a5551e4531c336d7949bfaf0`, tree
`f678fa505010119b0b546db116a82438e07f6405`. Raw `Cargo.lock` SHA-256:
`ead4421545a125267e33af84aeeea39f0e7ea2f60717a59ef0271d0e71855f1b`.
The exercise tag was excluded; no tag or release acceptance follows.

Steffen supplied `gentle-final-candidate-f35a58e6-20261013.tar.gz` in Downloads,
then copied it unchanged to `/private/tmp/` after macOS denied content access.
Codex independently verified the archive on 2026-10-13:

| Object | Bytes / SHA-256 |
| --- | --- |
| Archive | 103,172,261; `35e000d086f17348ea76dbaf88712f44b808abea045a72509cabc6d6034a1dcd` |
| [Original report](REPORT.md) | 8,790; `63ab711e81383b407d8b4a00b90c2873553551646944696bc127f9d27ec7bc52` |
| [Original per-file manifest](SHA256SUMS) | 27,537; `cabb5e49d9e904d55cf0d07b9aae0879ba269c2a3ca4d9899a199e7df189d9e0` |

The archive contains only relative `evidence/` regular files/directories; no
links or special files were observed. All 213 manifested file hashes pass,
and the extracted regular-file inventory exactly matches the manifest plus
the manifest itself. The original report and manifest are retained here
byte-for-byte as audit data, not instructions. Large artifacts remain in the
owner's archive; none is added to Git or executed locally. Inspection copy:
`/private/tmp/gentle-final-f35a58e6.76mWfg/evidence/`.
The two retained documents have scoped `text eol=lf` rules and an offline
LF/CRLF checkout regression, including a missing-rule negative control.
All 55 offline checkout tests pass on the documentation working tree;
`/private/tmp/gentle-f35a58e6-evidence-checkouts.log` SHA-256 is
`2e12be03af9f42a1f8dd7b387e16fdff0b818a9a3d8168ece8a7b1e3a32f3183`.
This verifies byte preservation, not the application at a newer source SHA.

To recheck an already safely extracted copy, from its `evidence/` directory:

```sh
shasum -a 256 --check --quiet SHA256SUMS
diff <(cut -c67- SHA256SUMS | sort) <(find . -type f ! -name SHA256SUMS | sort)
```

## Execution Scope

Glen's host is Debian forky/x86-64, Rust/Cargo 1.96.1, Python 3.14.8 and
DuckDB 1.5.5. Ordinary binaries use `dev`/default features; focused 08.04 and
scroll builds additionally disable incremental/debug information and use one
build job. The archive's self-identity log binds build `1791854013565` to this
SHA. Reported CLI/docs/GUI hashes are respectively
`d6e6f5de491a4f9dd8ef3a859acc2172508c700495a17b703de2d2c8d82b4e9d`,
`e6272131659ab179e85bcc575e649f6e438799444db89e1fc988bb23e15a9a81`,
`61c9c292e612e57620e57d355c201458e7e3bfbfa40622b88834797b685ee98c`.
Those executables are not present in the archive and are not independently
rehashed. The archived scroll executable is independently rehashed as
`27abc1fb5ac30112236317ab99d9291ba2292b9f202229dcc750c04993994aa1`.

| Area | Accepted Scope / Remaining Boundary |
| --- | --- |
| 08.16 TSS | One real-DuckDB regression and 3/3/2 learner queries pass; signed/zero scores, both strands and clipping retained. Fresh export has 12 outputs, report hash `104f5ad4ea7c99b74c5d97de0aef1dbc6389e96a7d61dc554293f242f222a0c9`. Exporter is exact `f35a58e6`; original producer remains `e4c94bf3`, not rewritten. Native lanes/readability not run. |
| 08.18 repeat | Two distinct genomic-forward `AGCCAG` loci, opposite repeat strands and two scored midpoint windows pass. Transcript-relative geometry remains unavailable and RNA ranking unset; no native/scientific verdict. |
| 07.04 / 10.03 | Import/projection, 26 gel tests, learner analysis/export and separate-process reopen pass. Out-of-calibration values remain null. Native navigation/readability not run. |
| FLNA | Stable selector, 36 transcripts/103 exons/73 junctions/six events and byte-identical SVG pass; retained PNG is 2400x14808, not a normal-size native screenshot. Historical Ensembl release and human exon/interpretation approval remain missing; draft stays inactive. |
| 04.09 | Focused contracts/parity and all real BLAST commands pass; real import remains performance-blocked. No completed real-reference specificity verdict, native presentation, order approval or Windows filesystem verdict follows. |
| Tutorials / provenance | 65 catalog entries, 30 executable chapters/six GUI contracts, LF/CRLF and nineteen visual-placement decisions pass. E-MTAB-15709 SDRF confirms H3K4me3; no broader DN biological claim. |

Retained logs show 60 scoped Python tests, 268 six-gate helper tests, 68
tutorial/primer/parity Rust tests and 19 focused 08.04 Rust tests passing.
These are Glen's executions, not local Codex Rust verification. The 08.04
native receipt records `status=fail`, specifically `sudo unshare` password
policy; its failure is not a product-test failure or an accepted native pass.
Live agent and `package-opt1` execution are not run.

Scroll receipt SHA-256 is
`0740931627a63106e191ae57f4552c0a2c11ad39cccff87d072f669f6013b903`.
Codex independently checks its exact SHA/profile, 24 unique fixture/size/direction
observations, three fixtures/four sizes, moved bounded viewports and unchanged
content. `native_latency_measured=false`; this closes only the deterministic
CPU event/next-frame contract at this source, not native latency.

## Independent CI Completion

After Glen assembled his report, Codex independently observes
[run 38329717079](https://github.com/smoe/gentle_rs/actions/runs/38329717079)
completed successfully at this exact SHA (attempt 1; updated
2026-10-13T03:43:15Z). Linux headless, Windows, macOS, release-safety and the
selected-platform aggregate succeed. The
[macOS job](https://github.com/smoe/gentle_rs/actions/runs/38329717079/job/115076477126)
log has 4,140 root tests passed/14 ignored and the focused primer suites
passed. Snapshot steps execute zero tests. Full Linux, scroll, authenticated
audit and 08.04 GUI jobs are skipped, not accepted. No package/container
workflow or publication was dispatched by Codex.

## Post-BLAST Import Follow-Up

The original report records CPU-bound interrupted imports for ACTB (5,999
HSP rows / 686,525 bytes) and TP53 (12,202 / 1,404,773). All three retained
execution manifests contain four completed, zero-exit, hash-bound searches;
TSVs remain untouched. No completed import state/report is archived. The
empty interrupted-command log does not provide a CPU stack profile; the
duration/CPU diagnosis is Glen's observation, not an independently profiled
root cause.

Source inspection identifies a concrete repeated-work path:
`primer_specificity_hit_from_blast` computes legacy terminal mismatches for
every hit; `get_sequence_region_with_cache` reloads the entire prepared FASTA
index on each call. The retained default handoff explicitly disables full
alignment, so merely skipping legacy fetches after complete BLAST alignment
would not address this reproducer. Do not lower hit limits, discard partial
HSPs, change terminal rules, relax fingerprint/stale checks or substitute a
nearby source's pass.

Next external diagnostic: use the same frozen source/profile and prepared
references, a fresh provenance-valid handoff and disposable state. Preserve a
bounded import-only CPU profile and file-access counts (only with approved
profiler permissions), peak RSS, elapsed/CPU time and interruption status.
Separate database hashing/annotation loading, per-hit FASTA-index/sequence
lookup and amplicon joins. Do not bypass isolation or profiler restrictions.
Only after locating the dominant cost, implement a scoped repair with exact
output/atomicity/stale-reference regressions, then rerun both real-reference
imports at the repair's new SHA. No performance repair or speedup is claimed
by this evidence-record commit.
