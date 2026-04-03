# VKORC1 / rs9923231 Warfarin-Response Promoter -> Luciferase Reporter (GUI Tutorial with Matching CLI Commands)

> Type: `GUI walkthrough + CLI mapping`
> Status: `manual/hybrid`
> Drift note: this page is hand-written, but it is intentionally tied to shared
> engine/CLI operations and links to executable workflow material where
> available.

This tutorial is a GUI-first walkthrough for planning a `VKORC1`
promoter-luciferase follow-up from the warfarin-sensitivity SNP `rs9923231`.

It has three goals:

1. start from one concrete pharmacogenomic alert rather than an abstract gene,
2. extract one biologically sensible reporter fragment around the reverse-strand
   `VKORC1` promoter and `rs9923231`, and
3. map each GUI step to a functionally equivalent CLI/engine route.

The intended assay logic is modest and explicit:

- the reporter tests allele-sensitive regulatory activity of the local
  `VKORC1` promoter fragment,
- it does **not** directly test warfarin response in a cell by itself,
- drug treatment can be added later as an assay condition once the promoter
  construct geometry is fixed.

Primary vector reference:

- Promega luciferase context sequence:
  [AY738222](https://www.ncbi.nlm.nih.gov/nuccore/AY738222)

## Scope

This tutorial covers computer-based planning and documentation:

1. prepare one GRCh38 reference genome in the active GENtle instance,
2. resolve `rs9923231` through dbSNP and extract an annotated local locus slice,
3. verify the reverse-strand promoter geometry at `VKORC1`,
4. extract one default reporter fragment that keeps the SNP inside the insert,
5. import one luciferase destination context,
6. preview one promoter->luciferase assembly candidate, and
7. design one verification-PCR handoff.

Wet-lab execution (cell model choice, transfection, luciferase readout,
warfarin dosing, statistics) is out of scope for this page.

## Default Biological Decision Used Here

To keep the tutorial reproducible and biologically cleaner than the earlier
TP73-style placeholder, this page uses one fixed default fragment:

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
5. zoom the linear map to roughly `2400..3600`

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

## Step 4: Import the Luciferase Destination Context

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

## Step 5: Preview One Promoter -> Luciferase Construct

Biological intent:

- place the selected `VKORC1` promoter fragment upstream of luciferase as the
  first reporter candidate,
- keep the fragment geometry fixed before branching into alternate alleles.

Current planning note:

- this tutorial uses the existing two-fragment preview scaffold as the shared
  parity baseline,
- once this fragment geometry is accepted, the next iteration should move into
  an explicit destination-opening Gibson plan rather than changing the promoter
  boundaries again.

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

## Step 6: Verification PCR

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

## Step 7: Next Iteration After This Tutorial

Once the reference-fragment geometry is accepted, the next rational step is:

1. keep the same extracted fragment boundaries,
2. branch one alternate-allele promoter fragment from that reference insert,
3. build the matched luciferase reporter pair,
4. compare promoter activity without changing fragment geometry between alleles.

That allele-branch step is intentionally not forced into this page yet. The
first task was to stop cloning "to the end of the gene" and settle on one
biologically sensible promoter fragment.

## Documentation and Handoff

Treat this as required project documentation:

1. record the extracted reference fragment:
   - `chr16:31095780..31096868`
   - local selection `=2413 .. 3501`
2. record the active prepared genome:
   - `Human GRCh38 Ensembl 116`
3. record the fetched locus sequence id:
   - `vkorc1_rs9923231_context`
4. record the extracted insert id:
   - `vkorc1_rs9923231_promoter_ref`
5. record the construct preview id:
   - `vkorc1_rs9923231_luc_construct_preview`

## Recommended Success Criteria

You should consider this tutorial successful when you can show all of the
following:

1. `rs9923231` was resolved and extracted against the prepared GRCh38 instance.
2. The resulting sequence view shows `VKORC1` as a reverse-strand promoter
   problem rather than a generic forward-strand locus.
3. The chosen reporter fragment keeps the SNP inside the insert and includes
   only a limited amount of transcribed context.
4. The luciferase destination context is imported under a stable id.
5. One construct preview and one junction-PCR planning step are recorded for
   handoff.
