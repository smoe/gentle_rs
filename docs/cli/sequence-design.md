# Synthetic Sequence Design

<!-- manual-counterparts:start -->
Related GUI reference: [Detailed Reference](../gui/detailed-reference.md).
These topic links do not imply identical feature coverage across interfaces.
<!-- manual-counterparts:end -->

## Synthetic Sequence Design

Start with [tutorial 06.07](../tutorial/06-07_synthetic_sequence_design.md) for an
executable, synthetic joint motif/global/window-GC example with exact-output
review, separate approval and the session-local undo boundary.

Experimental, outside the `.12` release gate. This is bounded synonymous motif
removal for explicitly synthetic coding inserts, not full DNA Chisel.

```sh
gentle_cli sequence-design plan @request.json --path preview.json
# Review status, exact DNA, nucleotide edits, completeness and non-claims first.
gentle_cli --state project.json sequence-design apply @preview.json --approve sha256:THE_REVIEWED_DIGEST
```

The same commands run in GUI Shell and through inner-agent suggestions. The
plan command with `--path` requires explicit agent confirmation for the exact
report-write attempt, including a cancelled report. Omitting `--path` adds no
operation-derived approval but does not remove candidate/tool requirements;
generic MCP retains its confirmation gate. Direct CLI invocation is explicit
user intent and gains no new universal `--confirm` flag. Successful report
writing is not search completion. Raw Plan/S1 op requests and nonempty workflows
containing only those inspections preserve post-loader project content and do
not implicitly create/rewrite project files or candidate sidecars. Named/shared-
shell inspection commands retain their existing no-change save behavior. Empty,
mixed and opaque wrappers keep legacy construction/save rules; parsing and
payload-file errors are unchanged. Plan reports require a fresh filename in an
existing directory: no overwrite, mkdir or automatic rename. The bound CLI/MCP
project filename (including absent/default files) and its derived `.candidates`
tree are excluded. A repeated preview needs a new explicit filename and required
agent approval. Generic MCP still requires `confirm=true`. Errors use existing
`EngineError` codes with `sequence_design_report.*` reasons.
Ambiguous non-ASCII same-parent protected-name comparisons are refused; use
another existing output directory. Native folder aliases are resolved for
directory-identity comparison, not a hostile parent-swap security guarantee.
Later explicit Save As and crash durability are outside the writer's guarantee.
The second command must not run until the user approves the exact feasible preview.
Use its actual `approval_digest`, never the illustrative placeholder above.
Planning never changes source/project sequences; application creates a derived
sequence without rerunning search, retains a receipt, and is undoable. Stale
input/protection/table/output bindings and existing output IDs are refused.

The approval digest binds reviewed bytes, not proof that GENtle produced their
search claims. Apply freshly verifies DNA constraints and exact edits, but does
not verify completeness or minimum edits, even for an unchanged preview.
Its receipt says `output_constraints_verified=true` and
`search_claims_verified=false`. Stored metadata places the approved report under
`submitted_proposal.report` with the mandatory tag
`submitted_proposal.verification="unverified_portable_preview"`. This is
unverified search evidence, not a new engine-certified optimization report.
Legacy receipts without verification fields and old unlabelled metadata
establish no verification claim. The single output container is named
`Synthetic sequence design`; inline input has synthetic-import provenance,
while loaded input retains a derived-source edge. A valid
two-edit output may be approved/applied even if a one-edit alternative exists;
application must not turn a submitted minimum claim into proof of optimality.

A minimal **synthetic** request (12 bp, protein MEF, one EcoRI recognition
sequence on both strands; not a natural gene or functional assay):

```json
{
  "schema": "gentle.dna_sequence_design_request.v1",
  "target": {"kind": "inline_sequence", "sequence": "ATGGAATTCTAA"},
  "purpose": "synthetic_coding_insert",
  "cds": {"start_0based": 0, "end_0based_exclusive": 12},
  "protein_sequence": "MEF",
  "genetic_code": 1,
  "protected_intervals": [],
  "avoid_motifs": [{"pattern": "GAATTC", "strand": "both"}],
  "max_evaluations": 4096,
  "output_seq_id": "synthetic_without_ecori"
}
```

For a loaded insert use `{"kind":"loaded_sequence","seq_id":"EXACT_ID"}`.
In the CLI, supply the same `--state project.json` to both planning and applying;
the loaded source must still match the preview when it is applied.
Coordinates are zero-based half-open, not shell `START..END` syntax. DNA must be
uppercase A/C/G/T, linear, with one complete forward code-1 CDS, literal ATG,
terminal stop and no internal stops. Protein excludes stop. Flanks, start/stop
triplets and protected bases are frozen. Motifs are finite IUPAC strings, not
enzyme aliases/regexes; all overlapping/boundary windows are tested.

Expected toy output is `ATGGAATTTTAA`, one nucleotide edit at position 8.
`search_exhausted` is unresolved, never infeasible; `cancelled` cannot be applied.
Incomplete feasible search does not prove minimum edits. Translation preservation
does not establish expression, folding, splicing or experimental suitability.
The original input must fit the 4,096-match diagnostic cap. An edited candidate
exceeding it is skipped with optimization incomplete, preserving any validated
best; otherwise the result is `search_exhausted`, not invalid input or proven
infeasibility. Apply still requires full output validation.
Original source annotations are omitted from the new sequence. Full limits and
algorithm: [design contract](../dna_sequence_optimization_plan.md).

Inner-agent example: "For this explicitly synthetic coding insert, preserve MEF,
the start/stop and my protected bases, and preview synonymous removal of GAATTC
on both strands. Show the exact edits and search status; do not apply yet."
The agent uses shared capabilities/readiness and the preview, not invented DNA
or hashes. Codon adaptation and broader specifications remain deferred.
Execution receipts do not expose report contents to the model. In Agent
Assistant, review the result and use **Use reviewed result in next prompt** only
if you consent to sharing its DNA and paths. This stages a draft, not approval
or execution; oversized results are refused intact, never truncated. Otherwise
review locally and supply the actual approval digest yourself.

### Explicit Conflict Search

Add `"search_strategy": "conflict_directed"` to the same request JSON to select
the experimental second solver. Commands, biological constraints and the
separate application approval are unchanged. Omission or `"full_enumeration"`
retains the original algorithm; unknown modes are refused. The inner agent
must disclose the choice, not silently change it after preview/approval.

The new algorithm tries synonymous changes in a current violating motif window,
then checks the whole sequence again. It can complete a constrained search
without enumerating unrelated codons, including when repair creates another
motif elsewhere. The report records `synonymous_conflict_search_v1`, requested
and effective candidate budgets, completeness and minimum-edit status.
Only a completed search proves minimum edits; an interrupted feasible search
is still a proposal with incomplete optimization, not a proof or speedup claim.
Changing strategy changes the approval digest: obtain and review a new preview.
The earlier toy has the same exact output under both strategies.

### CDS-Wide GC Walkthrough

For an explicitly synthetic insert, add `gc_content` to constrain the complete
declared CDS, including frozen start/stop. Bounds are inclusive basis points:
10000 means 100%, not a count of bases. Flanks do not contribute. There is no
organism default or separate strand policy; G+C is unchanged by reverse complement.
GC is enforced together with every requested motif, never as a post-search filter.

This hand-crafted 12-bp **MFK** example is not a natural gene. Save its JSON as
`request.json` before following the same preview/apply commands above:

```json
{
  "schema": "gentle.dna_sequence_design_request.v1",
  "target": {"kind": "inline_sequence", "sequence": "ATGTTTAAATAA"},
  "purpose": "synthetic_coding_insert",
  "cds": {"start_0based": 0, "end_0based_exclusive": 12},
  "protein_sequence": "MFK",
  "genetic_code": 1,
  "protected_intervals": [],
  "avoid_motifs": [],
  "gc_content": {"min_basis_points": 1666, "max_basis_points": 1667},
  "search_strategy": "conflict_directed",
  "max_evaluations": 4096,
  "output_seq_id": "synthetic_gc"
}
```

Here 16.66..16.67% becomes exactly **2 G/C bases out of 12**, using integer
rounding, not float tolerance. Input has one; expected output `ATGTTCAAATAA`
has two, still translates to MFK and changes only position 5 (zero-based).
Inspect `gc_content.input`, `gc_content.output`, count bounds, edits and search
status in `preview.json`. Only after approval apply its actual digest. GUI Shell
uses the same commands; then `ui open sequence-window synthetic_gc` opens the
derived DNA. Omit `search_strategy` to use enumeration with the same optimum.

Empty motifs are admitted only with explicit GC bounds. Invalid bounds differ
from infeasible integer count windows; a budget stop is still unresolved. Output
GC is absent without output DNA, not zero. Changing bounds requires a new
preview and approval. Satisfaction is not evidence of expression, folding or
synthesis success, and application does not verify preview search claims.

Inner-agent example: "This is a synthetic ATGTTTAAATAA coding insert for MFK.
Preview synonymous redesign with CDS [0,12), code 1 and GC bounds 1666..1667
basis points over that entire CDS. No forbidden motifs. Show integer GC counts,
edits and status; do not apply." The agent must ask for any missing request
fields, read the real preview and obtain separate exact-output approval.

### Sliding-Window GC Walkthrough

A reasonable GC fraction across the whole insert can conceal a locally extreme
stretch. Optional `gc_window` constrains **every complete overlapping window
inside the declared CDS**, stepping one base at a time. It includes frozen bases
but excludes flanks and shortened edge windows. `window_bp` must be 1..CDS length;
there is no stride or organism preset. Bounds are inclusive integer basis points,
using the same exact rounding as whole-CDS GC. Both constraints may be requested
together and are enforced jointly with motifs, translation and protection.

This is another literal synthetic MFK example, not a natural gene or a suggested
production threshold. Use its JSON as `request.json` with the preview/apply
commands above, in CLI or GUI Shell:

```json
{
  "schema": "gentle.dna_sequence_design_request.v1",
  "target": {"kind": "inline_sequence", "sequence": "ATGTTTAAATAA"},
  "purpose": "synthetic_coding_insert",
  "cds": {"start_0based": 0, "end_0based_exclusive": 12},
  "protein_sequence": "MFK",
  "genetic_code": 1,
  "protected_intervals": [],
  "avoid_motifs": [],
  "gc_content": {"min_basis_points": 1666, "max_basis_points": 1667},
  "gc_window": {"window_bp": 6, "min_basis_points": 1666, "max_basis_points": 1667},
  "search_strategy": "conflict_directed",
  "max_evaluations": 4096,
  "output_seq_id": "synthetic_window_gc"
}
```

The seven 6-bp windows are `[0,6)` through `[6,12)`. Each must contain exactly
one G/C base; the complete CDS must contain two. Expected output is
`ATGTTTAAGTAA`, changing only position 8, still encoding MFK. The previous
whole-CDS-only output, `ATGTTCAAATAA`, satisfies global GC but **fails** the local
window requirement. This demonstrates why the two constraints are not substitutes.

Review `gc_window.input` and `gc_window.output`: every window retains its exact
zero-based half-open `interval`, integer `gc_bases` and `satisfies_bounds`.
Also inspect count bounds, DNA, edits, status and non-claims before approving.
Omit `gc_content` for window-only design; omit `search_strategy` for enumeration
with the same tiny-space optimum. After separately approved application,
`ui open sequence-window synthetic_window_gc` opens the derived DNA in GENtle.

Invalid lengths/bounds differ from infeasible rounded count intervals. A budget
stop is unresolved; cancellation publishes no output/approval. Absent output
measurements are not zeros. Changing either GC constraint needs a new reviewed
preview. Passing these checks does not establish expression, folding or synthesis
success; apply checks output constraints, not submitted optimization claims.

Inner-agent example: "This is an explicitly synthetic ATGTTTAAATAA insert for
MFK, CDS [0,12), code 1. Preview GC 1666..1667 basis points globally and in every
6-bp CDS window, no forbidden motifs. Show exact window counts and edits; do not
apply." The agent must obtain missing fields and use the actual report/digest.

<!-- manual-page-nav:start -->

---

[CLI Manual](../cli.md) | [Previous: Variant and Reporter Evidence](variants-reporters.md) | [Next: Read-Only SNV Intent](snv-intent.md)
<!-- manual-page-nav:end -->
