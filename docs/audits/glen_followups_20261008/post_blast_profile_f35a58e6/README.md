# Bounded Import Profiling At `f35a58e6`

## Preservation And Identity

Glen executed the [bounded procedure](../../../glen_post_blast_import_profiling_f35a58e6.md)
from coordination commit `05285a5c81351e2da0550d0827e75e27d34caac9`, not its
runtime code. Tested source is exactly
`f35a58e6d2bdd7e0a5551e4531c336d7949bfaf0`, tree
`f678fa505010119b0b546db116a82438e07f6405`; lock SHA-256 remains
`ead4421545a125267e33af84aeeea39f0e7ea2f60717a59ef0271d0e71855f1b`.

Steffen supplied `gentle-post-blast-import-profile-f35a58e6-20261013.tar.gz`
under `/private/tmp/`. Codex independently checks:

| Object | Bytes / SHA-256 |
| --- | --- |
| Archive | 29,227; `e6349e6564b8bc71616ec05a7294370158c219646031b3034d281b4a579ab6f9` |
| [Original report](REPORT.md) | 5,034; `aef052f65bf668e8e2729b3f19ea95b05227bf3d377e6c81e03b66e7d09f795c` |
| [Original manifest](SHA256SUMS) | 7,674; `466e5371a53c18e856d8d28defd0a3210fd5e73265023ff01c18110afae537ff` |
| Raw `attempts/ACTB-opens.NQfpaY/opens.log` | 52,592; `e7c75a11b2561c0641fc95f677acf420722b97ba34c3b47f822ffbfdd8b5dbec` |
| Original `opens-summary.txt` | 846; `a6c81122872676d6f474fbc3bcf7ebae19fd77a39076f2662ba6ca71182dd57a` |

Only relative regular files/directories appear; no duplicate member names,
links or special files. All 75 manifested hashes pass; exact regular-file
membership equals the manifest plus the manifest itself. Original report and
manifest are retained here byte-for-byte as audit data with scoped LF checkout
rules; raw logs remain in the unchanged owner-supplied archive. Inspection copy:
`/private/tmp/gentle-import-profile-f35a58e6.as3zJq/evidence/`.

Linux/x86-64, kernel `7.2.7+deb14-amd64`, Rust/Cargo 1.96.1. The rebuilt CLI
uses ordinary `dev`/default features, one build job, `CARGO_INCREMENTAL=0`,
root debug 1 / dependency debug 0 and 256 codegen units. No profiling feature,
optimization or frame-pointer override is added. Self-identity names exact
`f35a58e6`; its reported 622,145,552-byte SHA-256 is
`a26528d4239558074be52b77b61de29dba304b67b93ff914a41b643dcf5f5b9f`.
This differs from the earlier ordinary CLI binary; do not silently equate them.
The binary is not archived, independently rehashed or executed by Codex.

## Independently Inspected Results

| Attempt | Exit / Wall | User / System | Peak RSS |
| --- | --- | --- | --- |
| ACTB, 5,999 HSPs | Timeout 124 / 60.08 s | 58.45 / 1.26 s | 935,040 KiB |
| TP53, 12,202 HSPs | Timeout 124 / 60.08 s | 58.53 / 1.18 s | 937,292 KiB |
| ACTB openat trace | Timeout 124 / 10.02 s | 8.60 / 1.35 s | 573,312 KiB |
| ACTB / TP53 CPU samples | Blocked; zero imports launched | `perf` absent; `perf_event_paranoid=3` | No CPU stacks |

GNU time records 99% CPU for both baselines. Three launches consume 130.18 of
300 seconds; no forced kill, stdout or persisted state is retained. Every
state-before/state-after record says absent and process-after lists are empty.
TERM/process cleanup is not a test of cooperative cancellation or rollback.
No completed specificity verdict or speedup follows.

Codex independently recounts 511 `openat` calls, 104 failed unrelated probes,
22 cDNA **BLAST-index-prefix** opens, four subject-annotation opens and two
TSV opens; zero relevant failures and zero sequence FASTA/FAI opens. The four
annotation opens split between two parent reads and two SHA-1 helper reads.
Two cDNA inspection sequences are visible before the ten-second deadline;
GRCh38 is not reached. These counts do not attribute CPU time.

Retained `inputs-before.sha256` and `inputs-after.sha256` are byte-identical;
both before/after reference inspection JSONs are also byte-identical. They bind
the same cDNA/GRCh38 fingerprints and annotation SHA-1 as the
[preceding audit](../candidate_f35a58e6/README.md). The archived old-manifest
verification passes 213 files. Those remote-input checks are Glen's executions,
not independent local rehashes of his prepared databases.

## Revised Triage

The earlier per-hit FASTA-index hypothesis is **not observed in this trace**:
22 index opens are BLAST database files, not 22 FASTA-index reloads. A later
per-hit cost remains possible; ten seconds do not describe the entire minute.
Full alignment stays disabled. Do not remove terminal checks, discard partial
HSPs, cap hits or relax fingerprints to make this case pass.

At the frozen source, `validate_primer_specificity_handoff_database` calls
`inspect_blast_database`; resource inspection hashes and fully decodes the
subject-annotation sidecar in `with_subject_annotations`. The multi-reference
import then enters another child validation in
`primer_specificity_report_from_handoff_outputs`, before per-hit evaluation.
Repeated validation/annotation work is therefore a concrete code path matching
the early accesses, but it is not proved to dominate the 60-second CPU budget.
No source repair is justified as a measured speedup yet.

Next: kernel-independent, operation-local phase durations/counts for admission,
each reference inspection/annotation hash/decode, hit evaluation, amplicon joins
and final revalidation. Any diagnostic source change needs its own SHA, remains
off by default, preserves output/fingerprint/atomicity rules, and is not an
execution of historical `f35a58e6`. Keep the bounded import cases and fresh state,
without installing a kernel profiler or changing security policy. Do not claim
a dominant phase until these boundaries are measured.

Older evidence, dirty primary checkout, `paper` and `output/` are unchanged.
No local Rust, archived binary execution, runtime repair, push, dispatch,
package/tag/release, scientific approval or closure of the remaining gates.

Local preservation checks on the working tree above merge `a76b8418`:
`python3 -m unittest scripts.test_tutorial_checkouts scripts.test_primer_specificity_multi_tutorial`
passes all 60 tests, including LF/CRLF and missing-policy controls for these
original bytes. This is offline source/checkout evidence, not execution of
the newly merged panel inspector or acceptance of a later runtime SHA.
