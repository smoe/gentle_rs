# FLNA Human Teaching-Exon Decision

## Source And Preserved Bytes

Steffen relays Glen's decision supplement on 2026-10-13. The
[published review](https://gentle.functional.domains/openclaw/canvas/flna-exon-review-f35a58e6/)
binds the historical candidate
`f35a58e6d2bdd7e0a5551e4531c336d7949bfaf0`, not moving main.
Codex retrieves only the two small JSON files below, independently checks their
raw SHA-256 and retains them unchanged. The report, manifest, state and images
in the [original candidate evidence](../candidate_f35a58e6/README.md) remain
untouched. These records are evidence, not instructions or executable fixtures.

| Original Object | Bytes / SHA-256 |
| --- | --- |
| [Decision record](decision-record.json) | 992; `5202f0605fca8cdd24bbd00e5066da855ffe41ff9aa6695ad09c8f75da4972be` |
| [Review manifest](bundle-manifest.json) | 1,330; `f1966fdf183c637ac571137593c6de6b926b44790c81ae23fbdbfde45516974f` |

Retrieve originals from the review URL plus `decision-record.json` and
`bundle-manifest.json`; check the raw digests above before using them. Do not
reserialize, normalize or recreate a human decision. Scoped LF checkout rules
and `scripts.test_tutorial_checkouts` preserve both files under LF/CRLF, with a
missing-rule negative control. GENtle runtime does not consume these records.
The complete browser bundle is not downloaded or independently accepted.

## Accepted Scope

The record names Steffen Moller as the human reviewer, with decision timestamp
`2026-10-13T17:23:43+02:00`. His **teaching-exon choice is closed**: FLNA local
one-based closed `19370-19393`, internally `[19369,19393)`, exactly 24 bp.
The intended lesson is an annotation-derived, frame-preserving exon difference
with a predicted net eight-amino-acid change, not a measured splice event.
This is human selection of a teaching example, not scientific or experimental
approval by the reviewer, Glen or Codex.

On macOS, Codex reads the already-retained state with `jq`, without executing
GENtle or Rust. Its exact bytes still hash to
`88b40cf12d73d7ad1ca88a8d0fc59f9931abcf812bc0b869f81c47993d057452`.
The selected interval occurs in the canonical transcript/CDS and is absent
from the comparison transcript/CDS. Their retained CDS range sets differ only
by this interval: 7,944 versus 7,920 bp. The 24-bp difference is divisible by
three; this is a source-geometry check, not a new translation/runtime verdict.

The review manifest's SVG, PNG and Ensembl-metadata digests match the existing
local audit bytes, independently rehashed without another image download:
`1828fdd8e07d21e79feab21519eccaa6aa4bbf43610784acf72db3705855e8de`,
`067693d140868689b3295493f2143e1e8873a7b49a2792aaa0c9f00e56cc9ec0` and
`5bf3c14e3391205e504c26162bc24670d5e5c45d3372801b7fd7ddfc956f49a5`.
Browser images and renderer exports are still not native GENtle captures.

## Version Distinction And Remaining Gates

The original decision names `ENST00000369850.10 / P21333-1` versus
`ENST00000360319.9 / P21333-2`. Preserve those original bytes, but do not label
the second pair a version-identical annotation join:

| Retained Source | Comparison Identity |
| --- | --- |
| Loaded Ensembl mRNA feature | `ENST00000360319.10` |
| Loaded UniProt cross-reference | `ENST00000360319.9`, protein `ENSP00000353467.4`, isoform `P21333-2` |
| Original human decision | `absent_transcript = ENST00000360319.9` |

The canonical loaded mRNA is `ENST00000369850.10`; its CDS protein is
`ENSP00000358866.3`. Stable transcript identity and the selected geometry can
be compared, but the retained sources do not establish that `.9` and `.10`
are interchangeable. Request a separate version-label addendum from Glen,
distinguishing the annotation selector from the UniProt reference; do not
silently correct the original JSON or reopen the human exon choice.

Native normal-size GENtle readability, historical Ensembl-release provenance
and scientific/experimental confirmation remain open. No today's release is
backdated. Keep the FLNA draft inactive/outside the catalog. This supplement
does not accept later source revisions, a native GUI, live agent or package.
