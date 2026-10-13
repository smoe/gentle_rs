# Bounded Post-BLAST Import-Only Profiling

Executed 2026-10-13 following the procedure from coordination commit
`05285a5c81351e2da0550d0827e75e27d34caac9`. The tested GENtle source is the
separate frozen target `f35a58e6d2bdd7e0a5551e4531c336d7949bfaf0`, tree
`f678fa505010119b0b546db116a82438e07f6405`. No BLAST search, source repair,
permission change, CI dispatch, tag, publication, or scientific acceptance was
performed.

## Admission

- Clean detached source worktree; raw `Cargo.lock` SHA-256
  `ead4421545a125267e33af84aeeea39f0e7ea2f60717a59ef0271d0e71855f1b`.
- Linux x86_64, kernel `7.2.7+deb14-amd64`; Rust/Cargo 1.96.1.
- Ordinary `dev`/default-feature CLI build. Exact build command is recorded in
  `build-profile.txt`; `CARGO_BUILD_JOBS=1`, `CARGO_INCREMENTAL=0`, and an
  explicit run-owned `--target-dir` were used. No optimization, profiling
  feature, or frame-pointer flag was added. The repository dev profile has
  `debug=1`, `codegen-units=256`; dependencies use `debug=0`.
- CLI identifies as source `f35a58e6`; binary SHA-256
  `a26528d4239558074be52b77b61de29dba304b67b93ff914a41b643dcf5f5b9f`
  (622,145,552 bytes).
- Original evidence archive: 103,172,261 bytes, SHA-256
  `35e000d086f17348ea76dbaf88712f44b808abea045a72509cabc6d6034a1dcd`.
  All 249 members were regular files/directories with safe paths; all 213
  manifest hashes passed. The already-present original extraction was used at
  its content-bound path and was not overwritten.
- ACTB inputs: 5,999 TSV rows / 686,525 bytes. TP53 inputs: 12,202 TSV rows /
  1,404,773 bytes. Both manifests contain four completed, exit-zero commands;
  every TSV size and SHA-256 matches. Parent snapshots are null and full
  alignment is disabled exactly as retained.
- Before and after the attempts, the retained bundle file hashes were
  byte-identical. The original archive hash was also unchanged.
- Exact-source reference inspection passed before and after:
  - Ensembl 116 cDNA fingerprint `sha256:3159519891bd8c0d4690e9c2f3de8eeb9621152b58c749e3d7d4920f4913263f`,
    465,769 sequences / 1,276,790,359 bases, annotation SHA-1
    `5f571d5775368aeb4e8a88c3ce4b1c8c138e3acf`, unmasked transcriptome index.
  - Human GRCh38 Ensembl 116 fingerprint
    `sha256:90a64ae8ff24fdbab67286bab563c95a6b2a04be37f50ad811474106c1ba3e7e`,
    706 contigs / 3,291,585,349 bases, genomic index.

## Attempt results

Each launched import used a distinct initially absent state path and the exact
argument vector retained in its `command.txt`. No launched attempt wrote a
state file or left a `gentle_cli` process alive.

### ACTB baseline — TIMEOUT

- 60.08 s wall; 58.45 s user; 1.26 s system; 935,040 KiB peak RSS.
- Exit 124 after `TERM`; no forced kill, stdout, persisted state, or specificity
  verdict.

### TP53 baseline — TIMEOUT

- 60.08 s wall; 58.53 s user; 1.18 s system; 937,292 KiB peak RSS.
- Exit 124 after `TERM`; no forced kill, stdout, persisted state, or specificity
  verdict.

### ACTB CPU sample — BLOCKED

No import was launched. `perf` is not installed and the existing kernel setting
is `perf_event_paranoid=3`. No `sudo`, policy change, or retry was attempted.
Consequently no CPU stacks or phase percentages are available.

### TP53 CPU sample — BLOCKED

Same profiler boundary as ACTB; no import was launched and no CPU stacks exist.

### ACTB openat trace — TIMEOUT

- 10.02 s wall; 8.60 s user; 1.35 s system; 573,312 KiB peak RSS.
- Exit 124 after `TERM`; no forced kill, stdout, persisted state, or surviving
  process. Raw trace: 511 `openat` calls (104 failed, all unrelated missing
  locale/config/library probes; zero failed relevant input/reference opens).
- Relevant opens: 22 cDNA index-prefix files, four subject-annotation opens,
  two BLAST TSVs, one parent handoff and one execution manifest. GRCh38 was not
  reached within 10 seconds.
- The trace contains two cDNA inspection/fingerprint sequences and two
  subject-annotation reads. This is evidence of repeated access, but without
  CPU stacks it neither assigns a dominant sampled phase nor proves that file
  opening/index loading dominates runtime. The four short helper children used
  about 0.35 CPU seconds in total; parent-side parsing/evaluation is not
  separable from this trace.

Total launched import wall budget was 130.18 seconds, below the five-minute
cap. Three imports were launched; the two blocked CPU attempts were not.

## Finding and limits

Both data sizes remain CPU-bound beyond 60 seconds on this ordinary dev build,
with nearly identical one-core utilization and approximately 915 MiB peak RSS.
Neither run reached a saved result. The bounded file trace confirms repeated
cDNA resource/annotation access early in ACTB, but CPU attribution is blocked;
therefore the dominant phase among fingerprinting, annotation loading,
per-hit sequence lookup, mismatch evaluation, and amplicon joining remains
unresolved. Successful upstream BLAST is not an import or biological verdict.

No repair is proposed or attempted from this evidence.
