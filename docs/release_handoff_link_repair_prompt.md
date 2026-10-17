# Prompt: Repair Deleted .12 Release-Handoff Links

Prepare the smallest documentation-only repair for active links to the deleted
`.12` candidate-preparation handoff. This is not a request for broader
context-menu work, release acceptance, or runtime changes.

## Repository And Starting Evidence

Repository: `/Users/u005069/GitHub/gentle_rs` (or your equivalent checkout).
Inspection snapshot: `150a339074e7a2d282567d791a8f50834b70ff9c` on 2026-10-11.
Record your actual HEAD and working-tree status; recheck findings rather than
assuming this snapshot is still current.

Read `AGENTS.md`, `docs/architecture.md`, `docs/decisions.md`,
`docs/roadmap.md`, `docs/maintenance_chore_plan.md`, the relevant release/testing
docs, and the deletion's history before editing. Preserve unrelated work,
particularly the dirty `paper` submodule and untracked `output/`.

Cleanup commit `afbef44451ed61e5f347c3aeb586e7124201cf22` deleted
`docs/internal_12_candidate_preparation.md` (296 lines). Treat that deletion as
intentional: do not automatically restore the file or recreate its long planning
ledger under another name. Its last contents can be inspected with:

```bash
git show 7c52082c24fca85817d70956146623eef6f24607:docs/internal_12_candidate_preparation.md
rg -n 'internal_12_candidate_preparation' docs --glob '*.md'
```

At the inspected snapshot, five active Markdown links remain:

| Document | Lost destination/purpose |
| --- | --- |
| `docs/roadmap.md` | `#corrected-ci-evidence`: historical CLI/JASPAR-lock recheck |
| `docs/roadmap.md` | `#ubuntu-runner-decision`: rolling-host inventory and pending owner decision |
| `docs/roadmap.md` | `#clawbio-compatibility-decision`: retained modes and removal-review scope |
| `docs/clawbio_gentle_integration_onepager.md` | `#clawbio-compatibility-decision`: mode/call-site inventory |
| `docs/release_notes/release_notes_v0.1.0-internal.12.md` | whole handoff: historical evidence, build-only instructions and open decisions |

Other matches in the changelog, historical session notes, or this prompt are
not necessarily broken links. Do not erase valid historical mentions just to
make the search return no results.

## Minimal Repair

1. Retarget each active link to a destination that actually supplies its stated
   information. Prefer existing topic owners: `docs/testing.md` for test
   history, `docs/release.md` for build/runner guidance, and the ClawBio onepager
   for compatibility scope. Add only the small missing paragraph or inventory
   needed. A source-pinned GitHub permalink to the deleted document is also
   appropriate for explicitly historical evidence; use its full retained SHA,
   never `main`. Do not substitute an unrelated audit or generic homepage.
2. Preserve evidence distinctions. The deleted handoff records Windows job
   `109794599955`, run `36686854112`, at
   `a51bbc064ea29d4e7a03c559c5a6e0a6f8e1bd29` as closing the historical
   CLI/JASPAR-lock recheck. It is not current `.12` acceptance. The separate
   CertUtil sharing-violation note for job `107075366344` remains distinct.
   Keep dated October 5 CI notes historical rather than presenting them as
   current-source failures or passes.
3. For live inventories, verify current source, not just the deleted prose:
   - Runner consumers: `.github/workflows/ci.yml`, `container.yml`,
     `release.yml`, and `release-candidate.yml`. Distinguish rolling hosts from
     the already pinned `ubuntu-24.04` Linux installer. Preserve the pending
     pin/Ubuntu-26-qualification decision; do not make it or edit workflows.
   - ClawBio: `DEPRECATED_REQUEST_MODE_EQUIVALENTS` in
     `integrations/clawbio/skills/gentle-cloning/gentle_cloning.py`, descriptor
     routes/examples, skill docs, and tests. Earliest-removal `.11` is
     eligibility, not permission to remove modes. Do not change compatibility.
4. Keep `.12` a development draft, `.11` published history intact, and all
   scientific, GUI, performance and package acceptance limitations unchanged.
   Do not execute build-only commands merely because the docs contain them.
5. Keep the roadmap within 200–400 lines, with no `Done` entry. Record the
   completed link repair concisely in `docs/CHANGELOG.md`. Do not restructure
   the roadmap or update unrelated release status.

Allowed edits: the three affected documents, minimal supporting additions in
`docs/release.md` or `docs/testing.md` if needed, and the changelog. No code,
workflow, dependency, fixture, generated tutorial, artifact or version changes.
Do not run local Rust builds/tests/checks, Claude consultation, workflow
dispatches, pushes, tag changes or publication. Leave commits local if made.

## Verification And Handoff

Use documentation-only checks before and after the repair:

```bash
python3 scripts/maintenance_chore.py release-gate
git diff --check
python3 scripts/maintenance_chore.py session-close --plan docs/release_handoff_link_repair_prompt.md
```

Preparation-time reproduction at the snapshot above: `release-gate` reports
**3 OK / 0 warnings / 1 failure**. The failure is the three roadmap links to
this missing file. The repair itself has not been performed.

The existing release-gate checker tests local roadmap paths, not heading
fragments or the other documents: manually verify every changed relative path
and fragment, including the release notes' extra directory level. For historical
permalinks, verify the file/heading at the specified Git revision. Do not claim
that these checks certify runtime behavior or release readiness.

Return the tested HEAD, changed files, a five-link old-to-new destination map,
checks/results and any remaining ambiguity. Explain intentional historical
matches and unrelated dirty artifacts. If source cannot support a replacement
inventory or claimed CI closure, report the gap instead of guessing.
