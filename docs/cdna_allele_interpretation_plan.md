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
| C2 authenticated 1-64-gene batch, unchanged counts/policy | `5d310958`: kernel batch/policy and service distinct-namespace/strict-binding regressions passed; original `2f1a3062` refusal remains historical |
| C2 CDS/UTR/noncoding/not-present projection, plus/minus | `5d310958`: hand-computed mirrored coordinate/codon and unavailable-model tests passed within the 17/17 kernel suite |
| C2 all same-unit codon combinations, inherited support | `5d310958`: TCC/TAC/TCA/TAA, split codon, disjoint/partial/low-depth/cross-sample kernel cases passed |
| C2 local phase relations, disconnected blocks | `5d310958`: inherited support, disconnected blocks and sample-isolation cases passed; no chromosomal cis/trans claim |
| C2 substitution, delins, alleles, protein/start/stop forms | `5d310958`: HGVS transcript-reference, start/stop and bounded-extension cases passed; multi-codon protein composites stay withheld |
| C2 JSON/TSV/summary projections and nonclaims | `5d310958`: saved-JSON regeneration/service and both actual CLI/Python process cases passed, including deterministic formats and nonclaims |
| C3 shared authenticated lowering and comparison byte preservation | `5d310958`: service hash-drift/byte checks passed; freshly built `4966a3a5` versus `5d310958` on identical synthetic files/paths produced identical raw JSON and unchanged inputs. Earlier work-ceiling 1/1 remains at `2f1a3062` |
| C3 compatible/discriminating/unique-model support | `5d310958`: linkage, retained-intron, cross-hypothesis, incomplete-view and foreign-contig kernel cases passed; no QNAME/transitive splice join |
| C4 frozen manifest, exact SNV/composite matching | `5d310958`: exact supported-composite, timestamp/indel/assembly/tamper and fail-closed manifest cases passed; synthetic labels only |
| Shared engine/shell/direct CLI/MCP/JS/Lua/Python routes | `5d310958`: headless 9/9, script-interfaces 11/11 with actual JS/Lua VMs and MCP confirmation-before-I/O, CLI/Python 6/6 with both real processes and zero skips. Profile counts overlap |
| Discovery, canonical parity, maturity and LF/CRLF | Glen `2f1a3062`: generator zero diff, registry 32/32, freshness 1/1, glossary 16/16, checkout/maturity 116/116. Preserve these source-bound passes without repeating unaffected campaigns |
| Real TP53 conditional receipt | `not_run`: Glen's synthetic campaign did not authorize a new biological data run. Separate authorized bound-artifact receipt; no inferred counts or private-data execution |
| Windows/macOS/native GUI | Same-source Windows compile/workspace/engine job passed; raw log has zero occurrences of the eight warning names. terminal same-source macOS/overall CI passed, with full Linux and optional GUI/audit skipped; native GUI and actual Windows CLI/Python paths/handles remain separate |

The kernel belongs in `gentle-engine::allele_hash_screen`; main-crate file
authentication, C1 conversion and table injection remain outside it. C3 shares
one lowering routine rather than opening a QNAME-based adapter pipeline.
The separate kernel crate's test target failed E0282 at `2f1a3062`, independently
of the root fixture failure. Its explicit tuple type compiles in Glen's
`5d310958` run, not a local Codex build. C1 keeps its existing
`selector.source_id == annotation.reference_id`
guard; the comparison source names saved analysis and remains different. Shared
raw annotation, gene, reference and geometry bind C2, not cross-namespace string
equality. Neither schema nor fixture IDs are homogenized.
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

The C1-C4 implementation and synthetic runtime ledger now has source-bound
external receipts. Native/data/scientific/release acceptance remains distinct;
source/checkout checks alone never supplied Rust/VM/process acceptance.

## Local Source Evidence, 2026-11-12

Initial authored-source checks before Glen's rejected `2f1a3062` return,
against changes based on `4966a3a5`, without Cargo/Rust compilation:

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
multi-codon protein composite. No new Claude consultation was run. At that
original source-only checkpoint, Rust, VM, process, historical-byte and native
acceptance were pending. Use the later external receipt below for implementation
runtime, not native acceptance.

## Repair Evidence After Glen's 2f1a3062 Return

The supplied report and archive match; all 22 manifest entries were independently
verified without extraction/execution. Preserve the failed VM/CLI/crate and
blocked historical-byte verdicts in the
[return](glen_next_verification_handoff.md#glen-return-at-2f1a3062), not a generic pass.

After the ID-namespace and tuple repairs, the combined companion/checkout/maturity
command selects 122 cases: **120 passed, two real-CLI skips** (the binary variable
was explicitly unset). This includes actual LF/CRLF Git checkouts and bound
request construction keeping annotation-reference and analysis-source IDs
different. Scoped syntax/formatting, documentation status and whitespace pass.
New Rust negatives preserve wrong-selector/shared-binding refusals; they remain
unexecuted here. No Rust build, provider/private-data run or push. The repair
handoff required Glen first to
compile the separate kernel test target, then execute the failed/new service,
VM/MCP and actual-process cases before historical comparison-byte acceptance.

## Scoped Completion Audit At 5d310958, 2026-11-12

Glen's clean Linux run now supplies the missing compile-first, 17/17 kernel,
9/9 headless, 11/11 script-interfaces, 1/1 selector-refusal and 6/6 CLI/Python
receipts (both process cases ran). Actual JS/Lua VM execution and MCP confirmation
preceded route acceptance. The original 16 synthetic requirements map to the
ledger above and the named passing tests in the raw streams; the documented
HGVS/nonclaim/CDS corrections are intentional, not omitted requirements.
The distinct identity namespaces and strict shared bindings remain exercised.
Historical old/new raw JSON is byte-identical, not just semantically equal.

Codex independently checked all 64 archive manifest entries, report equality,
test streams and old/new JSON/input manifests. Same-source Windows job logs
also establish compilation/workspace tests and zero occurrences of the eight
helper/import warning names, using rustc 1.99.1 rather than Linux 1.97.1. Full
identities, counts, raw-log attribution and remaining host boundaries are in the
[Glen return](glen_next_verification_handoff.md#glen-return-at-5d310958).

No further implementation defect is evidenced by this return. The original
handoff permits a real TP53 receipt **or the reason it is missing**: this run was
synthetic-only and did not authorize new biological data execution. That receipt
remains explicitly unrun, not inferred from synthetic results. Native GUI,
unproved Windows process paths/handles and selected-candidate release gates
remain separate follow-ups. A later read-only observation confirms same-source
macOS/overall CI success; full Linux and optional GUI/audit jobs skipped. OD-I1-I5 retain their
implemented conservative defaults; this closes neither scientific nor release
acceptance and licenses no private data, push, tag or publication.
