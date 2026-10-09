# Original Audit Archive Preservation

This is a new Codex integrity record, not a replacement for Glen's report,
screenshots, projects or benchmark evidence. The owner-approved download worked
after the hotel IP was allowlisted. Verification ran on macOS on 2026-10-08;
none of the four historical Linux executables was run.

`verification.json` is retained byte-for-byte from the local streaming verifier.
SHA-256: `1a41c227b6e91587e6d310114d6672d875008f51170bf75428da303b3b50dbb5`.
The reported source remains `30f23aa4cb084694870a139714d0f0fef1726378`.

## Origin And Retrieval

- Approved [archive](https://gentle.functional.domains/openclaw/canvas/gentle-30f23aa4-audit/tutorial-parity-30f23aa4-20261008.tar.zst): 319,263,072 bytes, SHA-256 `d257b98b95359c91d4639eb18d8a11008575d306fe7c8e890189721086437fe1`.
- Approved [per-file manifest](https://gentle.functional.domains/openclaw/canvas/gentle-30f23aa4-audit/tutorial-parity-30f23aa4-20261008.files.sha256): 139,612 bytes, SHA-256 `0f4f723115453cb1ed491654ac7671cd392116cb1f3b8c601decd31e8d128361`.
- Durable original: `/mnt/storage-box-1-gentle/Glen/handoffs/gentle/30f23aa4/tutorial-parity-30f23aa4-20261008.tar.zst`, with sibling checksum and per-file manifest. Public access is owner-controlled; a later HTTP refusal does not authorize bypassing it.
- Local retained download, metadata and verifier outputs: `/private/tmp/gentle-original-audit-30f23aa4`. The archive and binaries are intentionally not added to Git.

## Verification And Use

First check exact archive length and both pinned SHA-256 values, then run
`zstd -t`. Stream `zstd -dc` through Python's tar reader, accepting only the
declared root, safe raw relative paths, directories and regular files. Reject
links, special files, traversal, duplicates and unmanifested files. Hash each
regular file without extraction and require exact manifest membership and
length before accepting the bundle. The local command was
`python3 /private/tmp/gentle-original-audit-30f23aa4/verify_original_archive.py`;
its source hash and full 721-file receipt hash are in `verification.json`.

All 721 files and 313 directory members passed, including all four binaries,
208 `criterion-bench-audit` files and 312 `criterion-patz1` files. Decompressed
tar length is 1,406,382,080 bytes; regular-file payload is 1,404,807,210 bytes.
There were no missing, mismatched, duplicate, unmanifested or unsafe members.
The report still matches the unchanged committed report's hash.

A separate bounded inspection copied 82 already-verified data files (7,839,291
bytes) to the local `inspection/` directory, never executables. The saved
workflow has 11 sequences, nine containers and no arrangements; the 1,089-bp
inserts differ only at zero-based position 588, C/T. Original Criterion slope
estimates agree with Glen's report. Raw native and exported images were viewed
without recolouring or rewriting them; identical whole-map insert exports are
not evidence of the one-base change.

These records support item 1 of the six-item checklist. The checkout regression
preserves the new receipt's raw bytes under LF/CRLF. Integrity verification does
not rerun Glen's Linux tests, authenticate his binary source/profile, establish
a performance regression or imply package, clinical or scientific approval.
