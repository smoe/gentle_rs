## GENtle Shell (GUI)

<!-- manual-counterparts:start -->
Related CLI reference: [Interfaces and Command Discovery](../cli/interfaces.md) | [Shared Shell](../cli/shell.md) | [ClawBio and OpenClaw](../cli/clawbio.md).
These topic links do not imply identical feature coverage across interfaces.
<!-- manual-counterparts:end -->

### ClawBio Window

The separate ClawBio window reviews the current sequence/project context for
an external skill. Set the skill alias and `clawbio.py` executable, use
`Refresh Context`, then inspect `Request context` before pressing `Run`.
While a request is running, `Cancel` requests cancellation. Returned summaries,
errors and output links remain available for inspection. The window's ordinary
`Help` button only opens this manual; it neither invokes ClawBio nor captures
the screen. This is separate from the shared command shell described below.

### Shared Command Shell

The DNA-window toolbar includes a `Shell` button that opens a command panel.

Behavior:

- Uses shared command parsing/execution (`src/engine_shell.rs`).
- Same command set as CLI `gentle_cli shell`.
- Shows a command preview before execution.
- Maintains command output history in the panel.

Supported commands:

- Open `Help -> Shell Commands`.
- This reference is generated from `docs/glossary.json` and shows both:
  - usage syntax
  - task summary/description
- Use the `Interface` selector in that view to filter commands by language/access path:
  - `All`
  - `GUI shell`
  - `GUI menu`
  - `CLI shell`
  - `CLI direct`
  - `JS`
  - `Lua`

Common examples:

- `help` — show help catalog or command-specific help.
- `state-summary` — summarize loaded sequences/containers/metadata.
- `render-svg SEQ_ID linear|circular OUTPUT.svg` — export map SVG.
- `ui selection sequence-window SEQ_ID --range START..END` — set the active
  DNA-window selection in 0-based, end-exclusive coordinates.
- `display show tfbs` / `display hide restriction-enzymes` — toggle shared
  display layers used by sequence-window rendering.
- `genomes prepare GENOME_ID ...` — prepare reference genome cache/index.
- `tracks import-bed SEQ_ID PATH ...` — project BED track onto anchored sequence.
- `arrays inspect-microarray-track MANIFEST` — inspect a prepared Clariom D
  array-track manifest without changing state.
- `arrays project-microarray-track SEQ_ID MANIFEST --contrasts CSV --level probeset`
  — project verified genome-anchored array intervals into the active sequence.
- `arrays probe-regions --cel sample.CEL --gene TP73 --platform Clariom_D_Human --dry-run`
  — preflight arbitrary Affymetrix CEL probe/probeset-region inspection without
  running summarization.
- `arrays inspect-probe-region-output OUTPUT_DIR` — inspect a completed
  `probe_regions_oligo.R` output directory before plotting or projection.
- `arrays render-probe-region-output-svg OUTPUT_DIR OUTPUT.svg` — export a
  deterministic native SVG plot from completed probe-region helper output.
- `arrays project-probe-region-output SEQ_ID OUTPUT_DIR --contrasts CSV`
  — project direct-coordinate-compatible helper output into array features on
  a genome-anchored sequence.
- `candidates generate SET_NAME SEQ_ID --length N ...` — create candidate set.
- `guides oligos-export GUIDE_SET_ID OUTPUT_PATH ...` — export guide oligo set.
- `planning profile show --scope effective` — inspect merged planning profile
  (`global -> confirmed_agent_overlay -> project_override`).
- `planning consult cloning --format text` — get read-only cloning
  strategy/vector advice from the effective planning profile and catalogs.
- `planning protein-expression-handoff --format text` — get a read-only
  protein-expression handoff for high-yield requests.
- `planning sync pull @suggestion.json --source lab_manager` — register a
  pending advisory planning suggestion (apply via
  `planning suggestions accept SUGGESTION_ID`).

Status output note:

- `genomes status` / `helpers status` include optional
  `nucleotide_length_bp`, `molecular_mass_da`, and `molecular_mass_source`
  alongside source-type/source-path fields.
- `helpers status` may also include an engine-owned `interpretation` record for
  helpers with structured semantics, so GUI/ClawBio/planner flows can share
  one normalized helper-meaning layer.
- The same status payload now also includes `effective_cache_dir`,
  compatibility hints, component-level readiness, and a ready-to-run
  `prepare_command` only when the selected genome's required components are not
  valid. A stale whole-prepare activity does not override a valid BLAST
  component.
- Prepared reference/helper inspectors show a compact BLAST-component status;
  hover detail includes genomic-DNA versus transcriptome/cDNA kind,
  assembly/release, masking, prefix, BLAST database/tool versions,
  sequence/base counts, content fingerprint, validation warnings, and
  compatible operations.
- Reference/helper genome dialogs respect optional process-level shared cache
  overrides:
  - `GENTLE_REFERENCE_CACHE_DIR`
  - `GENTLE_HELPER_CACHE_DIR`
- `routines list` / `routines compare` include planning estimates when planning
  data is configured:
  `estimated_time_hours`, `estimated_cost`, `local_fit_score`,
  `composite_meta_score`, plus optional hands-on/unattended/elapsed ranges,
  time basis and stock-draw/purchasing-proxy/machine/base cost components.
  Unknown breakdowns are not displayed as zero.

Screenshot command status:

- `screenshot-window OUTPUT.png` is currently disabled by security policy.
- This reserved artifact-export command is distinct from the explicit,
  preview-before-send [GUI screenshot help](../gui.md#screenshot-capture-policy).
- Prefer deterministic SVG exports for reproducible figures; reviewed tutorial
  captures use the separate evidence-and-staging process, not this command.

<!-- manual-page-nav:start -->

---

[GUI Manual](../gui.md) | [Previous: Workspace Latency Diagnostics](workspace-latency-diagnostics.md) | [Next: Agent Assistant](agent-assistant.md)
<!-- manual-page-nav:end -->
