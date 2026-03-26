# TP73 Promoter -> Luciferase Reporter (GUI Tutorial with Matching CLI Commands)

> Type: `GUI walkthrough + CLI mapping`
> Status: `manual/hybrid`
> Drift note: this page is hand-written, but it is intentionally tied to shared
> engine/CLI operations and links to executable tutorial material where
> available.

This tutorial is a GUI-first walkthrough for planning a TP73 promoter luciferase
construct in GENtle.

It has two goals:

1. plan one concrete biological workflow (TP73 promoter fragment into a
   luciferase reporter context), and
2. learn how each GUI action maps to a functionally equivalent
   command-line command.

That GUI/CLI mapping helps in both directions:

- you can transfer the same process to other genes/targets,
- you can interpret agent-proposed shell/CLI commands with confidence.

Primary vector reference:

- Promega luciferase context sequence: [AY738222](https://www.ncbi.nlm.nih.gov/nuccore/AY738222)

## Scope

This tutorial covers computer-based planning and documentation:

1. retrieve TP73 genomic context,
2. extract promoter candidates,
3. import luciferase vector sequence,
4. preview assembly candidate(s),
5. design verification assays,
6. document and export results.

Wet-lab execution (transfection, luciferase readout, statistics) is out of
scope for this page.

## Prerequisites

1. GENtle desktop application running.
2. Genome catalog available at `assets/genomes.json`.
3. `assets/genomes.json` may be either:
   - kept at project defaults (recommended for common genomes), or
   - edited locally for custom/local mirror sources.
4. Internet connection for online fetch operations, or local fallback inputs.
5. Optional local fallback file for AY738222:
   - `data/tutorial_inputs/AY738222_promega_luciferase.gb`

## Step 1: Prepare and Retrieve TP73 Genomic Context

GUI:

1. `File -> Prepare Reference Genome...`
2. choose `Human GRCh38 Ensembl 116`
3. click `Prepare`
4. `File -> Retrieve Genomic Sequence...`
5. query gene `TP73`
6. select occurrence `1`
7. output id `grch38_tp73`

GUI/CLI mapping:

| GUI action | Engine operation | CLI equivalent |
| --- | --- | --- |
| Prepare reference | `PrepareGenome` | `gentle_cli genomes prepare "Human GRCh38 Ensembl 116" --catalog assets/genomes.json --cache-dir data/genomes` |
| Retrieve TP73 gene | `ExtractGenomeGene` | `gentle_cli genomes extract-gene "Human GRCh38 Ensembl 116" TP73 --occurrence 1 --output-id grch38_tp73 --catalog assets/genomes.json --cache-dir data/genomes` |

## Step 2: Extend Anchor and Extract Promoter Candidates

There are two GUI routes for the same underlying function:

1. sequence-window quick controls (anchor extension row directly under the
   sequence toolbar: `Extend 5'`, `Extend 3'`), and
2. full **Engine Operations panel** (shown in UI as `Engine Ops`) for complete
   anchored extraction settings.

GUI (recommended order):

1. open sequence `grch38_tp73`
2. extend `5'` by `2500 bp` (quick controls or Engine Operations panel)
3. open the anchored extraction form in the Engine Operations panel
4. set:
   - anchor = `FeatureBoundary(kind=gene,label=TP73,boundary=Start,occurrence=0)`
   - direction = `Upstream`
   - target length = `1200`
   - tolerance = `200`
   - output prefix = `tp73_promoter_fragment`
5. execute extraction

GUI/CLI mapping:

| GUI action | Engine operation | CLI equivalent |
| --- | --- | --- |
| Extend 5' anchor | `ExtendGenomeAnchor` | `gentle_cli genomes extend-anchor grch38_tp73 5p 2500 --output-id grch38_tp73_promoter_context --catalog assets/genomes.json --cache-dir data/genomes` |
| Anchored promoter extraction | `ExtractAnchoredRegion` | `gentle_cli workflow @docs/examples/workflows/tp73_promoter_luciferase_assay_planning.json` |

### Step 2b (Optional but Recommended): Dotplot Checkpoint

GUI:

1. switch sequence map to `Dotplot map`
2. start with the default local window size
3. click `Compute dotplot`
4. hover for coordinates, left-click to lock crosshair
5. use locked coordinate as reproducible reference when choosing final
   promoter candidate boundaries

Purpose:

- local repeat/structure inspection before committing a candidate for assembly
  and primer design.

## Step 3: Import Luciferase Vector AY738222 (GUI-Native)

Primary GUI route:

1. `File -> Fetch GenBank Accession...`
2. accession = `AY738222`
3. optional `as_id` = `promega_luciferase_ay738222`
4. click `Fetch and Import`

Offline fallback GUI route:

1. `File -> Open Sequence...`
2. choose local file `data/tutorial_inputs/AY738222_promega_luciferase.gb`

GUI/CLI mapping:

| GUI action | Engine operation | CLI equivalent |
| --- | --- | --- |
| Fetch AY738222 | `FetchGenBankAccession` | `gentle_cli shell 'genbank fetch AY738222 --as-id promega_luciferase_ay738222'` |
| Local file import | `LoadFile` | `gentle_cli op '{"LoadFile":{"path":"data/tutorial_inputs/AY738222_promega_luciferase.gb","as_id":"promega_luciferase_ay738222"}}' --confirm` |

## Step 4: Assembly Preview (Biology and Run Behavior)

Biological intent:

- combine a selected TP73 promoter fragment with a luciferase vector context to
  create a reporter construct candidate in the model for downstream validation
  design.

In GENtle, the routine assistant loads a stored cloning routine template, while
you provide the actual sequence inputs loaded in Steps 2 and 3.

GUI:

1. open `Patterns -> Routine Assistant...`
2. choose `Gibson Two-Fragment Overlap Preview`
3. bind:
   - `left_seq_id = tp73_promoter_fragment_1`
   - `right_seq_id = promega_luciferase_ay738222`
   - `overlap_bp = 20` (example)
   - if the vector binding is circular, use `Linearize Vector...` in this
     window before running preflight
4. run **input check** (`Validate only`)
5. run **apply run** (`Transactional`)

Biological meaning of `overlap_bp = 20`:

- this is the overlap length (20 nucleotides) of homologous sequence expected
  between promoter-fragment end and vector-fragment end,
- in overlap-based assembly (for example Gibson-style design), that homology
  determines how specifically and robustly fragments anneal at the intended
  junction,
- `20 bp` is a practical starting point: usually long enough for specific
  junction pairing while still easy to realize in designed fragment ends.

Definitions:

- **Input check** (`Validate only`): validates inputs/constraints and reports
  feasibility without changing project state.
- **Transactional run**: applies all routine steps as one all-or-nothing run;
  any failure rolls back all intermediate changes.

GUI/CLI mapping:

| GUI action | Engine route | CLI equivalent |
| --- | --- | --- |
| Routine validate-only | `macros template-run ... --validate-only` | `gentle_cli shell 'macros template-run gibson_two_fragment_overlap_preview --bind left_seq_id=tp73_promoter_fragment_1 --bind right_seq_id=promega_luciferase_ay738222 --bind overlap_bp=20 --bind assembly_prefix=tp73_luc_demo --bind output_id=tp73_luc_construct_preview --validate-only'` |
| Routine apply run (all-or-nothing) | `macros template-run ... --transactional` | same command without `--validate-only` |

## Step 5: Verification Assays (Why PCR and qPCR)

Why PCR here:

- after designing a promoter->vector construct, PCR primer pairs are used to
  verify the expected assembly (for example promoter/vector junction presence,
  insert size/orientation checks, colony/plasmid QC).

Why qPCR appears here:

- qPCR is **optional** for promoter-luciferase planning,
- it can be useful for short-amplicon quantitative checks (for example DNA
  copy-number normalization or expression-related follow-up workflows),
- standard luciferase reporter assays can be completed without qPCR if your
  experimental design does not require quantitative nucleic-acid readouts.

Background references:

- [PCR](https://en.wikipedia.org/wiki/Polymerase_chain_reaction)
- [Quantitative PCR (qPCR)](https://en.wikipedia.org/wiki/Quantitative_PCR)

GUI:

1. open construct preview sequence
2. open `Patterns -> PCR Designer...` (or command palette `PCR Designer`)
3. in linear map mode, paint one `ROI` (green) around the planned promoter/vector junction
4. paint one `Upstream` (red) and one `Downstream` (blue) primer window, then use `Add ROI to Queue` (or `Shift+drag` on ROI for immediate queue append)
5. optionally select one or more features and use DNA-window fallback `PCR ROI -> Add selected feature(s) to queue`
6. in the right pane `Design primer pairs`, set shared constraints/report base and run `Design Primer Pairs for queued regions` (or one-shot `Design Primer Pairs` for single ROI)
7. optionally run `Design qPCR Assays` from Engine Ops when quantitative assays are needed
8. review/export results (`Show` / `Export` / `Open` in PCR batch results, or report helpers by `report_id`)

For a dedicated selection-first walkthrough (multi-region queueing, optional
`ExtractRegion` copy artifacts, deterministic batch report naming), see:

- [`docs/tutorial/generated/chapters/13_pcr_selection_batch_primer_pairs_offline.md`](./generated/chapters/13_pcr_selection_batch_primer_pairs_offline.md)

GUI/CLI mapping:

| GUI action | Engine operation | CLI equivalent |
| --- | --- | --- |
| Primer design report | `DesignPrimerPairs` | `gentle_cli primers design @request.design_primer_pairs.json` |
| qPCR design report | `DesignQpcrAssays` | `gentle_cli primers design-qpcr @request.design_qpcr.json` |
| Report export | report export helpers | `gentle_cli primers export-report REPORT_ID out.json` and `gentle_cli primers export-qpcr-report REPORT_ID out.json` |

## Step 6: Documentation and Handoff

Treat this as required project documentation, not an optional extra.

GUI:

1. verify lineage includes:
   - TP73 retrieval/extension
   - promoter candidate selection
   - AY738222 import
   - assembly preview output
   - primer/qPCR report IDs
2. export sequence and SVG artifacts relevant for review
3. save project as `.gentle.json`
4. record parameter choices and IDs in your experiment notes

## Same Steps Across Interfaces (Short Summary)

The same logical workflow can be expressed via:

- GUI operations,
- direct `gentle_cli` commands,
- shared shell commands (`gentle_cli shell '...'`),
- other script integrations that call the same core engine functions.

In this tutorial, this means the same biological step is executed in the same
way regardless of GUI, CLI, or shell route.

## Recommended Success Criteria

1. The same promoter candidate can be regenerated from the same anchor/extract
   parameters.
2. The vector import ID is stable and traceable.
3. Primer/qPCR reports refer to the same construct ID.
4. GUI actions and CLI commands do the same thing at each step.

## Reference Workflow JSON (Appendix)

Reference example file:

- `docs/examples/workflows/tp73_promoter_luciferase_assay_planning.json`

Associated metadata in that file:

- format tag: `gentle.workflow_example.v1`
- `test_mode: "skip"`

Meaning:

- The JSON is a machine-readable workflow example used for reproducible
  GUI/CLI/agent consistency checks.
- It is maintained in-repository as reference documentation for this
  tutorial scenario.
- `test_mode: "skip"` indicates external dependencies are intentionally involved
  (for example online GenBank fetch and prepared-genome assumptions), so default automated tests only
  check the JSON structure for this example instead of always executing it
  end-to-end.
