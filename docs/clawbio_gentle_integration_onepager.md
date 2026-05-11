# ClawBio / GENtle Integration Contract

GENtle presents one ClawBio skill, `gentle-cloning`, as a thin wrapper around
deterministic `gentle_cli` command surfaces. ClawBio owns chat routing,
confirmation policy, and platform presentation. GENtle owns command execution,
artifact bundling, provenance, and machine-readable reports.

## Wire Schemas

Current schemas:

- Request: `gentle.clawbio_skill_request.v1`
- Result: `gentle.clawbio_skill_result.v1`
- Skill info: `gentle.clawbio_skill_info.v1`
- Runtime intents: `gentle.clawbio_skill_intents_runtime.v1`
- Intent descriptor: `clawbio.skill_intents.v1`

Versioning stance:

- Additive request modes, result fields, runtime-intent rows, or warnings do
  not require a v1 schema bump.
- A backward-incompatible shape change, field removal, or semantic change to an
  existing field requires a new schema.
- Deprecated request modes remain accepted for at least one release window and
  emit `warnings[]` naming the preferred thin form.

ClawBio should read `INTENTS.json` from disk at startup. It should call
`mode: "intents"` only when it wants to compare the installed wrapper's
descriptor hash with a local snapshot or during an explicit refresh/dev check.
Startup should not require a working GENtle runtime.

## Runtime Intent Refresh

Request:

```json
{
  "schema": "gentle.clawbio_skill_request.v1",
  "mode": "intents"
}
```

Result payload in `stdout_json`:

- `schema`: `gentle.clawbio_skill_intents_runtime.v1`
- `descriptor_sha256`: hash of installed `INTENTS.json`
- `supported_request_modes`: wrapper modes accepted by this installation
- `routes[]`: `intent_id`, trigger terms, request modes, examples, slot needs

This mode is wrapper-owned and does not invoke `gentle_cli`.

## Mode Audit

| Mode | Class | Stance |
| --- | --- | --- |
| `skill-info` | wrapper concern | Keep; metadata without runtime dependency. |
| `intents` | wrapper concern | Keep; runtime descriptor/hash surface without runtime dependency. |
| `version` | thin pass-through | Keep; maps to `gentle_cli --version`. |
| `capabilities` | thin pass-through plus probe | Keep; maps to `capabilities` and best-effort `ui intents`. |
| `state-summary` | thin pass-through | Keep; maps to `state-summary`. |
| `shell` | thin pass-through | Keep as preferred generic command surface. |
| `op` | thin pass-through | Keep as preferred operation surface. |
| `workflow` | wrapper concern | Keep; resolves workflow paths and bundles artifacts. |
| `raw` | thin pass-through | Keep for advanced callers that already know CLI argv. |
| `protein-residue-genomic-coordinates` | thin typed op | Keep; maps to one operation payload. |
| `gene-protein-2d-gel` | wrapper concern | Keep; parameterized workflow synthesis plus artifact naming. |
| `agent-plan` | safety wrapper | Keep for now; carries planner options and confirmation boundaries. |
| `agent-execute-plan` | safety wrapper | Keep for now; carries stored plan execution options. |
| `primer-preflight` | shell normalizer | Deprecate; prefer `mode: "shell"`. |
| `primer-seed-from-feature` | shell normalizer | Deprecate; prefer `mode: "shell"`. |
| `primer-seed-from-splicing` | shell normalizer | Deprecate; prefer `mode: "shell"`. |
| `primer-design` | shell normalizer | Deprecate; prefer `mode: "shell"`. |
| `primer-report-list` | shell normalizer | Deprecate; prefer `mode: "shell"`. |
| `primer-report-show` | shell normalizer | Deprecate; prefer `mode: "shell"`. |
| `primer-report-export` | shell normalizer | Deprecate; prefer `mode: "shell"`. |
| `qpcr-seed-from-feature` | shell normalizer | Deprecate; prefer `mode: "shell"`. |
| `qpcr-seed-from-splicing` | shell normalizer | Deprecate; prefer `mode: "shell"`. |
| `qpcr-design` | shell normalizer | Deprecate; prefer `mode: "shell"`. |
| `qpcr-report-list` | shell normalizer | Deprecate; prefer `mode: "shell"`. |
| `qpcr-report-show` | shell normalizer | Deprecate; prefer `mode: "shell"`. |
| `qpcr-report-export` | shell normalizer | Deprecate; prefer `mode: "shell"`. |
| `cdna-pcr-test` | shell normalizer | Deprecate; prefer `mode: "shell"`. |
| `cdna-qpcr-test` | shell normalizer | Deprecate; prefer `mode: "shell"`. |
| `transcript-qpcr-panel` | shell normalizer | Deprecate; prefer `mode: "shell"`. |
| `restriction-cloning-pcr-handoff` | shell normalizer | Deprecate; prefer `mode: "shell"`. |
| `restriction-cloning-pcr-handoff-seed` | shell normalizer | Deprecate; prefer `mode: "shell"`. |
| `restriction-cloning-vector-suggestions` | shell normalizer | Deprecate; prefer `mode: "shell"`. |
| `restriction-cloning-handoff-list` | shell normalizer | Deprecate; prefer `mode: "shell"`. |
| `restriction-cloning-handoff-show` | shell normalizer | Deprecate; prefer `mode: "shell"`. |
| `restriction-cloning-handoff-export` | shell normalizer | Deprecate; prefer `mode: "shell"`. |
| `pcr-protocol-cartoon` | shell normalizer | Deprecate; prefer `mode: "shell"`. |
| `exon-skip-plan` | shell normalizer | Deprecate; prefer `mode: "shell"`. |
| `exon-skip-materialize` | shell normalizer | Deprecate; prefer `mode: "shell"`. |

Class meaning:

- Thin pass-through: one request shape maps directly to one CLI shell/op/workflow
  surface.
- Shell normalizer: wrapper convenience over a shared shell family. These stay
  accepted but should not grow new biology or routing semantics.
- Wrapper concern: path resolution, artifact bundling, runtime descriptor
  exposure, or presentation metadata that belongs at the ClawBio wrapper layer.

## Deprecation Table

Deprecated in `v0.1.0-internal.8`, earliest removal `v0.1.0-internal.10`:

- Primer/qPCR/report/CDNA modes listed above as shell normalizers.
- Restriction-cloning handoff helper modes listed above as shell normalizers.
- `transcript-qpcr-panel`, `pcr-protocol-cartoon`,
  `exon-skip-plan`, and `exon-skip-materialize`.

Each deprecated mode still executes and adds a `warnings[]` entry naming the
equivalent `mode: "shell"` form.

## Drift Guards

CI checks should keep the following in sync:

- every `route.plan[].input` in `INTENTS.json` resolves under `examples/`;
- every `examples/*.json` file is either routed or explicitly allowlisted as a
  bootstrap/follow-on/dev example;
- `SKILL.md` front-matter `trigger_keywords` is generated from
  `INTENTS.json` trigger terms plus `trigger_keyword_additions.json`;
- `mode: "intents"` reports the same descriptor routes and request modes the
  installed wrapper supports.

Regenerate trigger keywords with:

```bash
python3 scripts/generate_clawbio_trigger_keywords.py --write
```

Verify without writing:

```bash
python3 scripts/generate_clawbio_trigger_keywords.py --check
```

## Service Scope Neutrality

The shared engine default for `services handoff` is provider-neutral
`scope = "default"`. ClawBio-specific calls should pass `--scope clawbio`
explicitly. Omitting `--scope` is accepted but emits a warning so callers can
update without losing backward compatibility.
