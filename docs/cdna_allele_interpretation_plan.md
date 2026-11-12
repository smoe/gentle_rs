# Isoform-Bound cDNA Interpretation: C1-C4

Owner authorization: 2026-11-12, implement C1 through C4. The earlier C1-only
authorization is extended, not a release/scientific acceptance. No push, live
provider call, private archive regeneration or local Rust build is authorized.
Original proposal: owner-forwarded Glen/Claude prompt, 2026-11-11. Read-only
Claude consultation was offered; no new consultation has been invoked.

## Invariants And Clarifications

- Keep comparison and normalized annotation wires/bytes unchanged. New
  interpretation is a downstream read-only report over authenticated sources.
- Reuse strict GTF, shared exon geometry, original saved-analysis lowering and
  the existing genetic-code asset. No second aligner/parser/genetic code.
- CDS completeness is independent of mRNA end tags. Missing/version-conflicting
  models withhold the appropriate conditional HGVS/protein fields.
- HGVS 21.1.5 (2026-09-17): apply rules on the named reference. Consecutive
  transcript substitutions and the same-codon one-intervening-base exception
  may form a delins even across genomic introns. Never use genomic adjacency
  as a necessary transcript-HGVS condition. Keep genomic alleles separate.
- Predicted protein notation uses parentheses. Start loss is uncertain; stop
  loss scans only the authenticated downstream transcript sequence. No assumed
  translation initiation, transcript completion or biological confirmation.
- Codon combinations require definite same-unit calls. Local phase-block
  relationships are not chromosomal cis/trans facts. Do not combine samples,
  disconnected blocks or mapping hypotheses to manufacture a molecule/isoform.
- Explicit negative nonclaim text is allowed. A keyword ban is not a scientific
  validation strategy. Match frozen variant labels exactly; no clinical labels.
- Use only redistributable synthetic tests for implementation. A real TP53
  receipt needs separately available source-matched binaries and authorized
  bound artifacts; it cannot be inferred from synthetic checks.

## Delivery And Evidence Ledger

| Slice / explicit requirement | Current evidence / remaining proof |
| --- | --- |
| C1 coding models, phases, versions, protein binding | Implemented `d927ffa8` + compiler repair `36c7e954`; Glen's clean Linux run at `4966a3a5` passes compile/core/GTF/MCP/shell, real JS/Lua VMs, 4/4 CLI companion cases (zero skips), direct replay and canonical parity; separate routine ENOSPC does not invalidate that receipt. Native-host and affected newer-source checks remain separate |
| C2 authenticated 1-64-gene batch, unchanged counts/policy | Source and refusal regressions authored; raw/hash/assembly/model/cumulative limits and copied sample/status records need Rust execution |
| C2 CDS/UTR/noncoding/not-present projection, plus/minus | Source and hand-computed coordinate/codon/unavailable-model regressions authored; runtime pending |
| C2 all same-unit codon combinations, inherited support | TCC/TAC/TCA/TAA, split codon, disjoint/partial/low-depth/cross-sample tests authored; runtime pending |
| C2 local phase relations, disconnected blocks | Source copies IDs/fingerprints/method; local-pattern classes avoid chromosomal cis/trans; Rust tests pending |
| C2 substitution, delins, alleles, protein/start/stop forms | Source/tests follow official HGVS 21.1.5, including old-stop-as-1 extension; no Rust result yet |
| C2 JSON/TSV/summary projections and nonclaims | Pure saved-JSON projections and escaping/refusal/nonclaim regressions authored; execution pending |
| C3 shared authenticated lowering and comparison byte preservation | Shared lowering/hash-drift guard and same-source rerun regression authored; historical old/new binary byte proof pending externally |
| C3 compatible/discriminating/unique-model support | Same-unit/all-hypothesis/retained-intron tests authored; no QNAME/cross-hypothesis joins; runtime pending |
| C4 frozen manifest, exact SNV/composite matching | Strict synthetic-only source/UTC/license/allele/assembly/clinical-field/composite tests authored; runtime pending |
| Shared engine/shell/direct CLI/MCP/JS/Lua/Python routes | Wired to one service; typed parser/MCP confirmation/real VM/process tests authored, execution pending |
| Discovery, canonical parity, maturity and LF/CRLF | Discovery/maturity and matrix delta projected; Python source/raw-binding checkout tests pass; external canonical Rust generator/freshness pending |
| Real TP53 conditional receipt | Pending external execution; no new private data or hard-coded expected counts |
| Windows/macOS/native GUI | Pending where applicable; headless or Linux proof is not native acceptance |

The kernel belongs in `gentle-engine::allele_hash_screen`; main-crate file
authentication, C1 conversion and table injection remain outside it. C3 shares
one lowering routine rather than opening a QNAME-based adapter pipeline.
OD-I1 same-unit-only linkage, OD-I2 inherited minimum support, OD-I3 identifiers
without clinical text, OD-I4 JSON plus pure text projections and OD-I5 withheld
unverified multi-codon protein composites remain conservative defaults; richer
choices still need an explicit decision.

## Primary References

- [HGVS DNA delins](https://hgvs-nomenclature.org/stable/recommendations/DNA/delins/)
- [HGVS protein substitution](https://hgvs-nomenclature.org/stable/recommendations/protein/substitution/)
- [HGVS protein extension](https://hgvs-nomenclature.org/stable/recommendations/protein/extension/)

Official pages identified version 21.1.5 (2026-09-17) when checked 2026-11-12.
No live nomenclature dependency is used by the engine.

Work toward the entire ledger. Source/checkout checks are not Rust/VM/process
acceptance; do not close the objective while those proofs remain absent.

## Local Source Evidence, 2026-11-12

Against working changes based on `4966a3a5`, without Cargo/Rust compilation:

- `python3 -m unittest scripts.test_cdna_allele_interpretation scripts.test_tutorial_checkouts scripts.test_capability_maturity -q`:
  121 cases, 119 passed, two actual-CLI skips, including the Python-client
  JSON/text delegation regression (mocked subprocess, not runtime acceptance).
- `python3 -m unittest scripts.test_coding_transcript_models -v`: two source
  checks pass, two actual-CLI skips (separate from the command above).
- `python3 -m unittest discover -s integrations/python/tests -v`: seven wrapper
  tests pass; these use the existing synthetic process harness, not Rust binaries.
- Documentation status, Python maturity freshness and whitespace pass.
  `rustfmt` parses modified modules; this is syntax/formatting, not type checking.
- Session-close: four OK, two warnings (intentional dirty/untracked work and
  manual plan fidelity), zero failures. `paper`, `data/`, `output/`, owner prompts,
  `Cargo.toml` and `Cargo.lock` remain outside the change.

Intentional prompt corrections are recorded above and in the contract: reference-
appropriate delins, independent mRNA/CDS completeness, local-pattern rather than
chromosomal phase terminology, negative nonclaims permitted, no speculative
multi-codon protein composite. No new Claude consultation was run. Rust, VM,
same-source process, historical comparison-byte proof and native acceptance
remain pending; this ledger is not a green implementation acceptance receipt.
