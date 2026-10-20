# Capability Maturity And Portable Skill Discovery

## What This Answers

"Can I call this operation?" and "What evidence supports using it here?" are
different questions. A primer workflow that ran on synthetic Linux inputs has
not thereby demonstrated native Windows usability, real-reference specificity,
isoform discrimination or wet-lab performance. An absent local Primer3 binary
is an availability problem, not a retraction of a historical test.

Keep four lanes separate: interface coverage, scoped maturity evidence, current
availability/readiness and execution outcome. Existing review/approval checks
remain authoritative regardless of metadata or a model's recommendation.

## Discovery And Compatibility

The three `integrations/clawbio/skills/*/SKILL.md` files retain conservative
names and descriptions under the Agent Skills 1,024-character description
limit. The longer generic inventory lives in `references/capability_scope.md`
and the skill body. The repository checker validates its deliberately narrow
plain-name/single-folded-paragraph representation, not arbitrary YAML.

Existing OpenClaw frontmatter/triggers stay unchanged. Nested host metadata is
not a claim of strict conformance with the base specification's string-valued
metadata map; host-specific parsing/installation remains a host concern. No
Claude- or Codex-specific mandatory fields are added. Their native invocation
policies remain distinct from scientific evidence and GENtle approval.

These are external ClawBio discovery descriptions, not the inner Agent
Assistant's instructions. Its system prompt, documentation bundle, bounded
project-fact context and registered command routes are unchanged. The maturity
catalog is available through explicit introspection, not injected into every
model request. An offline source-level regression guards this separation;
it does not establish live-model behavior or native acceptance.

`INTENTS.json` remains `clawbio.skill_intents.v1`. Two additive route fields
identify the capability (`skill/intent_id`) and bind each ordered plan step to
its request mode, scope and existing registry references. Full workflows retain
their operation references, including mixed read/write steps; these references
do not grant permission. Unregistered helpers and caller-supplied workflows
are explicitly unresolved, not equated with a nearby command.

Old ClawBio readers can ignore the additions. The inspected local ClawBio
`4929f33` loader preserves unknown route fields; a local metadata-only check
accepted all 59/32/6 routes. This is not a claim that every deployed host/version
or its routing/execution was exercised. Existing descriptor and
approval digest algorithms are untouched. Descriptor bytes do change: deploy
the matching catalog/descriptor/runtime together and renew approvals where the
descriptor participates. No existing proposal is silently migrated.

## The Shared Catalog

`crates/gentle-protocol/src/capability_maturity.rs` owns the typed catalog,
validation and exact scoped-assessment helper. The engine compiles
`assets/capability_maturity.json`; shared `introspect capabilities` returns it
with missing/ambiguous registry diagnostics. `introspect all` contains the same
record. No adapter-specific scientific logic or new operation is introduced.

`scripts/generate_capability_maturity.py` projects explicit intent/step bindings,
catalog versions, referenced requests/workflows/fixture hashes and the public
`assets/capability_maturity_evidence.json` archive. Pure existing wrapper argument
builders check the curated bindings; no command is executed. Sources and outputs
are sorted and byte-bound without a generation timestamp. Hashing and parsing
use the same retained bytes, not newline-normalized text. The generator rejects
stale step bindings, duplicate IDs, altered evidence artifacts and unsafe paths.

Per-skill `CAPABILITY_MATURITY.json` uses schema
`gentle.clawbio_capability_maturity.v1`. This name intentionally does not
overload `capabilities`, which remains the existing UI-intent discovery mode.
The generic wrapper's `skill-info` and `intents` add a compact optional
`capability_maturity` reference, metadata-file hash and status, without invoking
GENtle. `available` means the catalog can be inspected, not the capability is
ready. Local descriptor, catalog, adapter and request bytes are checked;
non-local source bindings and installed-runtime applicability are not probed.
Missing/malformed/stale sidecars return explicit diagnostics while ordinary
skill discovery and execution remain available under their existing rules.

## Evidence Is Scoped, Not A Badge

The initial archive is empty. All 97 routes remain **unassessed**. Existing
parity rows, tutorial reviews and release notes are not silently converted into
attestations. Evidence must explicitly retain:

- Stable capability ID and exact definition hash.
- Tested source SHA, tool version and binary SHA-256.
- Platform, architecture, build profile and exact feature set.
- Input class and named resource hashes.
- Evidence kind: execution, native GUI, biological review or performance.
- Verdict, claimed maturity, assertions, limitations and a hashed public artifact.

Discovery never claims runtime binding. The shared `assess` helper needs an
explicit exact scope and one evidence kind. It neither probes nor authenticates
the supplied context: callers must establish actual runtime/input identities
independently. Different source, binary, platform, input, resources or kind do
not match. Missing/ambiguous registry bindings fail closed. Historical definition
mismatches remain readable but cannot contribute. Matching withdrawn, failed or
incomplete records take precedence over pass; otherwise the least mature
matching assertion wins. No timestamp selects a convenient newer pass.

A content hash is not an author's signature or scientific validation. Unknown
fields, malformed scope and empty/missing inventories are refused rather than
being interpreted as zero work or universal verification. No metadata authorizes
network access, mutations, orders or releases. Feedback receipts, automatic
sharing and telemetry are not part of this change. A future ClawBio planning/UI
consumer can display the optional records; policy enforcement beyond existing
safety checks requires its own review.

## Verification

```sh
python3 scripts/generate_capability_maturity.py --check
python3 scripts/generate_clawbio_trigger_keywords.py --check
python3 -m unittest scripts.test_capability_maturity -v
python3 -m unittest scripts.test_tutorial_checkouts -v
pytest integrations/clawbio/skills/gentle-cloning/tests/ -q -k 'not wrapper_rasterizes_svg_artifact_via_real_gentle_cli_route and not shipped_graphics_example_emits_png_first_artifacts'
```

Use `--write` to regenerate after an intentional source change. Every raw-bound
source/output has a scoped checkout rule. Real disposable LF/CRLF checkouts
verify freshness, including negative controls that remove output/input LF rules.
All added evidence fixtures are hand-crafted synthetic temporary records, not
real acceptance reports. Legacy route/plan contents are pinned against
`ac91375b`; optional additions cannot silently change their execution contract.

Glen/CI must execute the authored Rust regressions and build/check:

```sh
cargo test --locked -p gentle-protocol capability_maturity
cargo test --locked --lib capability_maturity_introspection
cargo check -q --locked
pytest integrations/clawbio/skills/gentle-cloning/tests/ -q
```

The two excluded Python graphics cases invoke Cargo indirectly; the complete
wrapper suite belongs with Glen/CI when local Rust execution is prohibited.

Python passes establish metadata behavior only. They do not establish Rust,
native GUI, live-model, biological, performance or release acceptance.
