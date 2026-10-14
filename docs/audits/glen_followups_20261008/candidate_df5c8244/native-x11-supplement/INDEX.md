# GENtle df5c8244 native-X11 supplement

This is a small, read-only handoff for inspection by Codex. It supplements the
verification of exact source revision
`df5c8244787d8c920ee171ee55167de55e8a7294`; it does not rerun, retarget, or
upgrade either blocked native-GUI verdict.

Original evidence root:
`/home/clawbio/.openclaw/workspace/artifacts/gentle-df5c8244-verification-20261014`

## Included

- Existing GENtle stderr/stdout, Xvfb and Openbox logs, and input-hash receipts
  for `native-x11/microarray/` and `native-x11/gel/`.
- The gel window listing captured after the editor action.
- The two exact launch commands recovered from the recorded tool calls.
- Recorded process/timing observations; no values were recreated by rerunning.
- Both complete-tree manifest versions and their original companion files.

## Excluded

- Screenshots, project states, binaries, credentials, raw session transcripts,
  and the large/full evidence archive.

## Manifest clarification

- `manifests/early-0b6a9d95/SHA256SUMS` is the 194-line manifest copied into
  the small handoff at 2026-10-14 14:30:45 +0200. Its SHA-256 is
  `0b6a9d95c689a7cb0f0ce2becf118b2a9c6192ed39738879f0b867cb39b16430`.
- `manifests/final-713a8d03/SHA256SUMS` is the 195-line final manifest written
  at 2026-10-14 14:31:18 +0200. Its SHA-256 is
  `713a8d034e2c7964d29e6350e3fa79d026726e7d28b15686d66597f5f9754e1a`.
- The sole line difference is the final manifest's added `cleanup.txt` entry:
  `0f7cfc444537b81e073a39068d393b80316016ad054f1500bd1a5da2c420aff8`.
- The early companion is preserved byte-for-byte. It names the original full
  manifest path, which now contains the later version, so running that old
  companion in place against today's full path fails by design. Hash the copied
  early `SHA256SUMS` itself to verify the preserved early version.

Verify this supplement from its root with:

```sh
sha256sum -c HANDOFF-SHA256SUMS
```
