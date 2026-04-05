# VKORC1 / rs9923231 PGx Alert -> Mammalian Luciferase Reporter (GUI Tutorial with Matching CLI Commands)

> Type: `GUI walkthrough + CLI mapping`
> Status: `manual/hybrid`
> Drift note: this page is hand-written, but it is intentionally tied to shared
> engine/CLI operations and links to executable workflow material where
> available.

This tutorial is a GUI-first walkthrough for taking one pharmacogenomic alert
around warfarin sensitivity and turning it into one concrete human-cell
reporter-design handoff in GENtle.

The intended division of labor is explicit:

- ClawBio interprets the alert: `rs9923231` sits upstream of reverse-strand
  `VKORC1` and suggests a regulatory follow-up rather than a coding follow-up.
- GENtle deterministically retrieves the locus, extracts one reproducible
  promoter fragment, imports one mammalian reporter backbone, previews one
  reporter construct, and exports one reviewable construct map.

It has five goals:

1. start from one concrete pharmacogenomic alert rather than an abstract gene,
2. explain why `rs9923231` motivates a regulatory follow-up around the
   reverse-strand `VKORC1` promoter,
3. extract one biologically sensible reporter fragment around that promoter and
   the SNP,
4. end on one clean mammalian reporter construct-map export rather than only a
   transient preview state, and
5. map each GUI step to a functionally equivalent CLI/engine route.

The intended assay logic is modest and explicit:

- the reporter asks whether a local `VKORC1` promoter fragment carrying
  `rs9923231` changes luciferase output in transfected human cells,
- it is a human-cell regulatory assay story, not a bacterial-expression story,
- it does **not** directly prove a drug-response phenotype by itself,
- warfarin treatment can be layered on later, but only after the promoter
  construct geometry is fixed.

Primary backbone choice used here:

- baseline assay class: one promoterless mammalian luciferase reporter plasmid
  suitable for transient transfection in human cells,
- deterministic importable record used by the current tutorial:
  [AY738222](https://www.ncbi.nlm.nih.gov/nuccore/AY738222) under the stable id
  `promega_luciferase_ay738222`,
- later escalation, not default: adenoviral delivery if transfection
  efficiency becomes the real bottleneck,
- later alternatives if the biology question changes:
  `pGL4.23[luc2/minP]`-style enhancer follow-up or NanoLuc-class mammalian
  reporters when sensitivity becomes the main constraint.

## Scope

This tutorial covers computer-based planning and documentation:

1. prepare one GRCh38 reference genome in the active GENtle instance,
2. resolve `rs9923231` through dbSNP and extract an annotated local locus slice,
3. verify the reverse-strand promoter geometry at `VKORC1`,
4. extract one default reporter fragment that keeps the SNP inside the insert,
5. import one primary mammalian reporter backbone,
6. preview one promoter->luciferase assembly candidate,
7. export one readable construct map,
8. write down one reproducibility bundle for handoff, and
9. design one verification-PCR handoff.

Wet-lab execution (cell model choice, transfection, luciferase readout,
warfarin dosing, normalization strategy, statistics) is out of scope for this
page.

Related extension idea:

- if a later assay wants to reuse the warfarin/Factor X connection more
  directly, the `factor Xa` recognition site can be embedded between other
  indicator modules so cleavage is read out as a signal ratio rather than as a
  purification step
- that extension would belong in a mammalian or adenoviral expression backbone,
  not in the bacterial `pGEX` context used elsewhere in the repository

## Default Biological Decision Used Here

To keep the tutorial reproducible and biologically cleaner, this page uses one fixed default fragment:

- reference assembly: `GRCh38`
- SNP: `rs9923231` at `chr16:31096368`
- prepared genome: `Human GRCh38 Ensembl 116`
- fetched context: `rs9923231 +/- 3000 bp`
- default reporter insert: `chr16:31095780..31096868`

Inside the fetched `6001 bp` locus slice, that default insert is:

- local GUI selection formula: `=2413 .. 3501`
- engine `ExtractRegion`: `from=2412`, `to=3501`

Why this fragment:

- it includes about `200 bp` of transcribed context past the `VKORC1` TSS,
- it keeps `rs9923231` well inside the fragment instead of at the border, and
- it adds `500 bp` of additional upstream context beyond the SNP on the
  biologically upstream side.

That is a much better starting point for a promoter reporter than "all the way
to the end of the gene body".

Assumptions used by this tutorial:

- baseline experiment type: transient transfection into one human cell line,
- baseline construct type: promoter fragment cloned upstream of luciferase in a
  mammalian reporter plasmid,
- baseline comparison: one reference promoter construct now, with the
  alternate-allele branch as the immediate next iteration,
- non-goal for this page: direct wet-lab proof of warfarin sensitivity.

Why keep `~200 bp` past the TSS at all:

- it avoids making the reporter boundary coincide exactly with the annotated
  TSS,
- it preserves a little immediate downstream / `5' UTR` context that can matter
  for promoter-proximal regulation, and
- it keeps the construct biologically modest: enough local context to make the
  promoter fragment interpretable, but not so much gene body that the insert
  stops reading like a promoter-reporter design.

## Prerequisites

1. GENtle desktop application running.
2. Genome catalog available at `assets/genomes.json`.
3. The active GENtle instance should either already have
   `Human GRCh38 Ensembl 116` prepared or be able to prepare it locally.
4. Internet connection for:
   - dbSNP rsID resolution (`FetchDbSnpRegion`)
   - online GenBank retrieval (`AY738222`)
5. Optional local fallback file for AY738222:
   - `data/tutorial_inputs/AY738222_promega_luciferase.gb`

## Step 1: Prepare the Reference and Resolve the SNP Context

GUI:

1. `File -> Prepare Reference Genome...`
2. choose `Human GRCh38 Ensembl 116`
3. click `Prepare` if needed
4. `File -> Fetch GenBank / dbSNP...`
5. rsID = `rs9923231`
6. genome = `Human GRCh38 Ensembl 116`
7. `+/- flank bp` = `3000`
8. output id = `vkorc1_rs9923231_context`
9. click `Fetch Region`

Notes:

- the dbSNP dialog already defaults to `rs9923231`, which makes this tutorial
  one-keystroke reachable in a prepared instance.
- the fetch status line should report staged progress such as contacting NCBI
  Variation, waiting, parsing, placement resolution, and annotated slice
  extraction.

GUI/CLI mapping:

| GUI action | Engine operation | CLI equivalent |
| --- | --- | --- |
| Instance preflight | `GenomeStatus` summary route | `gentle_cli genomes status "Human GRCh38 Ensembl 116" --catalog assets/genomes.json --cache-dir data/genomes` |
| Prepare reference | `PrepareGenome` | `gentle_cli genomes prepare "Human GRCh38 Ensembl 116" --catalog assets/genomes.json --cache-dir data/genomes` |
| Resolve rsID and extract locus | `FetchDbSnpRegion` | `gentle_cli op '{"FetchDbSnpRegion":{"rs_id":"rs9923231","genome_id":"Human GRCh38 Ensembl 116","flank_bp":3000,"output_id":"vkorc1_rs9923231_context","annotation_scope":"full","catalog_path":"assets/genomes.json","cache_dir":"data/genomes"}}' --confirm` |

## Step 2: Verify the Locus Biology Before Cloning

GUI:

1. open sequence `vkorc1_rs9923231_context`
2. keep at least these layers visible:
   - `Gene`
   - `mRNA`
   - `Variation`
3. note that `VKORC1` is on the reverse strand
4. note that `rs9923231` sits upstream of the `VKORC1` transcription start in
   biological terms even though the genomic coordinate is larger
5. zoom the linear map to roughly `2400..3600` with the +/- symbol, the mouse scrolling wheel or by entering these coordinates followed by pressing the apply button.

Interpretation used by this tutorial:

- `VKORC1` TSS is near local position `2613`
- `rs9923231` is at local position `3001`
- the SNP is therefore about `388 bp` upstream of the TSS in this prepared
  reference
- the local window also contains overlapping annotations on both strands;
  those are contextual, not the insert we plan to clone

This is the point where the tutorial intentionally pauses for biological sense
checking. The planned insert is not "the whole gene" and not "everything up to
the far gene boundary". It is one TSS-centered regulatory fragment.

## Step 3: Extract the Default Reporter Fragment

GUI:

1. in the DNA-window `Selection formula` field, enter:

   ```text
   =2413 .. 3501
   ```

2. click `Apply Sel`
3. confirm the highlighted span covers:
   - a small amount of transcribed `VKORC1` context near the TSS
   - the SNP marker
   - additional upstream sequence beyond the SNP
4. click `Extract Sel`
   - if the selection is not active yet, click `Select visible` or drag-select
     first; `Extract Sel` is disabled until a non-empty selection exists
   - this button lives in the wrapped action row directly below the display
     chips in the DNA window, beside `Queue PCR selection` and `PCR ROI`
   - if you do not see it immediately, widen the DNA window a little and make
     sure the sequence is still in `Mode = Region`
   - equivalent fallback: use `Engine Ops -> Extract Region` with the same
     selected bounds
5. use output id `vkorc1_rs9923231_promoter_ref`

Important coordinate note:

- the GUI selection formula uses the displayed `1-based` positions
- `ExtractRegion` uses `0-based [from,to)` coordinates
- so the matching engine extraction is `from=2412`, `to=3501`

GUI/CLI mapping:

| GUI action | Engine operation | CLI equivalent |
| --- | --- | --- |
| Extract default reporter fragment | `ExtractRegion` | `gentle_cli op '{"ExtractRegion":{"input":"vkorc1_rs9923231_context","from":2412,"to":3501,"output_id":"vkorc1_rs9923231_promoter_ref"}}' --confirm` |
| Full tutorial skeleton | workflow replay | `gentle_cli workflow @docs/examples/workflows/vkorc1_rs9923231_promoter_luciferase_assay_planning.json` |

Optional interpretation note:

- within the extracted promoter fragment itself, `rs9923231` now sits at local
  position `589`, which is much more appropriate for later reference-vs-variant
  reporter comparison than a border-adjacent SNP.
- the retained `~200 bp` of transcribed context is intentional and should be
  read as "do not cut exactly at the annotated TSS", not as "clone deep into
  the gene body"

## Step 3b (Optional): Summarize TFBS Around the SNP

If you already imported a JASPAR motif set and want one quick computational
check before cloning, summarize factor families around the SNP against the full
context window.

Shared shell example:

```bash
gentle_cli shell 'features tfbs-summary vkorc1_rs9923231_context --focus 2900..3100 --context 0..6001 --limit 25'
```

Interpretation:

- focus window `2900..3100` spans `100 bp` on either side of the SNP
- context window `0..6001` is the whole fetched locus slice
- the result helps you see which TF families cluster near the SNP versus
  elsewhere in the local context

This is advisory, not proof of mechanism.

## Step 4: Import the Primary Mammalian Reporter Backbone

Why this backbone class is appropriate here:

- a promoterless mammalian reporter keeps the regulatory question clean because
  luciferase output depends on the cloned `VKORC1` fragment rather than on a
  pre-existing vector promoter,
- luciferase is a convenient readout for human-cell transfection assays,
- this stays within a standard plasmid-transfection workflow and avoids
  escalating immediately to viral delivery.

Primary GUI route:

1. `File -> Fetch GenBank / dbSNP...`
2. GenBank accession = `AY738222`
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

Interpretation used by this tutorial:

- treat this imported sequence as the primary mammalian luciferase reporter
  backbone for the planning exercise,
- if you later substitute a different mammalian reporter vector, keep the same
  biological question and promoter-fragment geometry fixed.

## Step 5: Preview One Promoter -> Mammalian Luciferase Construct

Biological intent:

- place the selected `VKORC1` promoter fragment upstream of luciferase as the
  first human-cell reporter candidate,
- keep the fragment geometry fixed before branching into alternate alleles.

Current planning note:

- this tutorial uses the existing two-fragment preview scaffold as the shared
  parity baseline,
- once this fragment geometry is accepted, the next iteration should move into
  an explicit destination-opening Gibson plan rather than changing the promoter
  boundaries again,
- the important improvement for this page is that it should now end with a
  clean construct map you can actually show, not just with an internal
  construct id.

GUI:

1. open `Patterns -> Routine Assistant...`
2. choose `Gibson Two-Fragment Overlap Preview`
3. bind:
   - `left_seq_id = vkorc1_rs9923231_promoter_ref`
   - `right_seq_id = promega_luciferase_ay738222`
   - `overlap_bp = 20`
   - `assembly_prefix = vkorc1_luc_demo`
   - `output_id = vkorc1_rs9923231_luc_construct_preview`
4. run **input check** (`Validate only`)
5. run **apply run** (`Transactional`)

GUI/CLI mapping:

| GUI action | Engine route | CLI equivalent |
| --- | --- | --- |
| Routine validate-only | `macros template-run ... --validate-only` | `gentle_cli shell 'macros template-run gibson_two_fragment_overlap_preview --bind left_seq_id=vkorc1_rs9923231_promoter_ref --bind right_seq_id=promega_luciferase_ay738222 --bind overlap_bp=20 --bind assembly_prefix=vkorc1_luc_demo --bind output_id=vkorc1_rs9923231_luc_construct_preview --validate-only'` |
| Routine apply run | `macros template-run ... --transactional` | same command without `--validate-only` |

## Step 6: Inspect and Export a Readable Construct Map

GUI:

1. open `vkorc1_rs9923231_luc_construct_preview`
2. keep the construct in circular/standard-map view
3. for a cleaner documentation-style map, keep:
   - `Features = on`
   - `Gene = off`
   - `mRNA = off`
   - `TFBS = off`
   - `Restriction enzymes = off`
   - `GC = off`
   - `ORFs = off`
   - `Methylation = off`
4. if the sequence text panel is open, hide it so the map is the visual focus
5. use `Export SVG`
6. save the result as something stable, for example:
   - `vkorc1_rs9923231_luc_construct_preview.svg`

Why this step matters:

- it gives the tutorial a concrete endpoint artifact,
- it makes the construct readable enough for review/discussion,
- it uses the same shared `RenderSequenceSvg` route that the CLI and future
  ClawBio-facing artifact generation can reuse.

GUI/CLI mapping:

| GUI action | Engine operation | CLI equivalent |
| --- | --- | --- |
| Hide sequence panel | `SetDisplayVisibility` | `gentle_cli op '{"SetDisplayVisibility":{"target":"SequencePanel","visible":false}}' --confirm` |
| Keep map panel visible | `SetDisplayVisibility` | `gentle_cli op '{"SetDisplayVisibility":{"target":"MapPanel","visible":true}}' --confirm` |
| Hide non-essential construct tracks | repeated `SetDisplayVisibility` | repeat for `GeneFeatures`, `MrnaFeatures`, `Tfbs`, `RestrictionEnzymes`, `GcContents`, `OpenReadingFrames`, `MethylationSites` |
| Export circular construct SVG | `RenderSequenceSvg` | `gentle_cli op '{"RenderSequenceSvg":{"seq_id":"vkorc1_rs9923231_luc_construct_preview","mode":"Circular","path":"vkorc1_rs9923231_luc_construct_preview.svg"}}' --confirm` |

## Step 7: Verification PCR

GUI:

1. open `vkorc1_rs9923231_luc_construct_preview`
2. open `Patterns -> PCR Designer...`
3. paint one `ROI` around the promoter/vector junction
4. paint one `Upstream` and one `Downstream` primer window
5. run `Design Primer Pairs`
6. export the report if you want a handoff artifact

CLI parity:

| GUI action | Engine operation | CLI equivalent |
| --- | --- | --- |
| Primer design report | `DesignPrimerPairs` | `gentle_cli primers design @request.design_primer_pairs.json` |

For this tutorial, junction PCR is the default verification step. qPCR remains
optional and should only be added when your experimental plan really needs it.

## Step 8: Next Iteration After This Tutorial

Once the reference-fragment geometry is accepted, the next rational step is:

1. keep the same extracted fragment boundaries,
2. branch one alternate-allele promoter fragment from that reference insert,
3. build the matched luciferase reporter pair,
4. compare promoter activity without changing fragment geometry between alleles.

That allele-branch step is intentionally not forced into this page yet. The
first task was to stop cloning "to the end of the gene" and settle on one
biologically sensible promoter fragment.

## Reproducibility Bundle and Handoff

Treat this as the minimum handoff bundle for ClawBio -> GENtle continuity:

1. `report.md`
   - summarize the PGx alert,
   - explain why the follow-up is regulatory,
   - record the chosen promoter fragment and backbone choice,
   - state the main assumptions and non-goals.
2. `result.json`
   - record the stable ids, coordinates, backbone choice, and exported artifact
     paths in a machine-readable way.
3. commands or workflow replay
   - keep the stepwise CLI commands or the workflow path needed to replay the
     same construct preview.

This repository now carries one example bundle for that purpose under:

- `docs/tutorial/reproducibility/vkorc1_rs9923231_promoter_reporter/report.md`
- `docs/tutorial/reproducibility/vkorc1_rs9923231_promoter_reporter/result.json`
- `docs/tutorial/reproducibility/vkorc1_rs9923231_promoter_reporter/commands.sh`

Minimum fields to preserve:

1. extracted reference fragment:
   - `chr16:31095780..31096868`
   - local selection `=2413 .. 3501`
2. active prepared genome:
   - `Human GRCh38 Ensembl 116`
3. fetched locus sequence id:
   - `vkorc1_rs9923231_context`
4. extracted insert id:
   - `vkorc1_rs9923231_promoter_ref`
5. imported reporter backbone id:
   - `promega_luciferase_ay738222`
6. construct preview id:
   - `vkorc1_rs9923231_luc_construct_preview`
7. exported construct map path:
   - for example `vkorc1_rs9923231_luc_construct_preview.svg`

## Recommended Success Criteria

You should consider this tutorial successful when you can show all of the
following:

1. `rs9923231` was resolved and extracted against the prepared GRCh38 instance.
2. The resulting sequence view shows `VKORC1` as a reverse-strand promoter
   problem rather than a generic forward-strand locus.
3. The chosen reporter fragment keeps the SNP inside the insert and includes
   only a limited amount of transcribed context.
4. The mammalian luciferase reporter backbone is imported under a stable id.
5. One construct preview is exported as a readable circular map rather than
   existing only as a transient in-memory state.
6. One reproducibility bundle exists for handoff.
7. One junction-PCR planning step is recorded for handoff.

## Planned Companion Tutorial, Not This Page

One follow-on ClawBio-motivated tutorial is now explicitly reserved as a second
story rather than folded into this first one:

- first tutorial: promoter-fragment transfer into a mammalian luciferase
  reporter backbone to evaluate promoter performance,
- planned second tutorial: extend that baseline construct into one
  warfarin-on-target follow-up where a second cassette is interrupted by one
  Factor X target site.

That companion changes the biological question. The current page asks
"what promoter fragment should we build and compare?" The planned follow-up
would ask "how do we extend that construct family to probe a downstream,
warfarin-relevant functional consequence?" Keep those as separate tutorials so
the first handoff stays biologically disciplined.
