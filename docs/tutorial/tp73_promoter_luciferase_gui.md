# TP73 Promoter -> Luciferase Reporter (GUI-first Tutorial)

This tutorial defines a practical planning workflow for a TP73 promoter
luciferase assay with a GUI-first path and explicit command-line parity.

Primary vector reference:

- Promega luciferase sequence (NCBI Nucleotide):
  [AY738222](https://www.ncbi.nlm.nih.gov/nuccore/AY738222)

Canonical workflow skeleton:

- `/Users/u005069/GitHub/gentle_rs/docs/examples/workflows/tp73_promoter_luciferase_assay_planning.json`
- schema: `gentle.workflow_example.v1`
- mode: `test_mode = "skip"` (contains external-input placeholders by design)

## Scope and intent

The tutorial covers in-silico planning and communication:

1. target-locus retrieval and anchored promoter extraction
2. promoter-fragment candidate generation
3. reporter-vector import and assembly preview
4. primer/qPCR design for build verification
5. lineage/report export for collaboration

Wet-lab execution (transfection/luciferase readout/statistics) is out of scope.

## Prerequisites

1. GENtle desktop app running.
2. Genome catalog available (`assets/genomes.json`).
3. Optional local fallback file for the vector:
   - `data/tutorial_inputs/AY738222_promega_luciferase.gb`
4. Optional promoter evidence tracks (`.bed`, `.bigWig`, `.vcf`) for
   prioritization.

## Step-by-step (GUI first, parity visible)

### Step 1: Prepare and retrieve TP73 genomic context

GUI:

1. `File -> Prepare Reference Genome...`
2. choose `Human GRCh38 Ensembl 116`
3. `Prepare`
4. `File -> Retrieve Genome Sequence...`
5. gene query: `TP73`
6. select occurrence `1`
7. output id: `grch38_tp73`

Parity map:

| GUI action | Engine operation | CLI parity (optional) |
| --- | --- | --- |
| Prepare reference | `PrepareGenome` | `gentle_cli genomes prepare "Human GRCh38 Ensembl 116" --catalog assets/genomes.json --cache-dir data/genomes` |
| Retrieve TP73 gene | `ExtractGenomeGene` | `gentle_cli genomes extract-gene "Human GRCh38 Ensembl 116" TP73 --occurrence 1 --output-id grch38_tp73 --catalog assets/genomes.json --cache-dir data/genomes` |

### Step 2: Extend anchor and isolate promoter candidates

GUI:

1. open sequence `grch38_tp73`
2. in Engine Ops, extend `5'` by `2500 bp`
3. open `Extract Anchored`
4. anchor = `FeatureBoundary(kind=gene,label=TP73,boundary=Start,occurrence=0)`
5. direction = `Upstream`
6. target length = `1200`, tolerance = `200`
7. output prefix = `tp73_promoter_fragment`
8. run extraction

Parity map:

| GUI action | Engine operation | CLI parity (optional) |
| --- | --- | --- |
| Extend 5' anchor | `ExtendGenomeAnchor` | `gentle_cli genomes extend-anchor grch38_tp73 5p 2500 --output-id grch38_tp73_promoter_context --catalog assets/genomes.json --cache-dir data/genomes` |
| Promoter candidate extraction | `ExtractAnchoredRegion` | Use `gentle_cli workflow` with the canonical JSON example for full payload parity |

### Step 3: Import luciferase vector sequence (AY738222)

GUI:

1. open the GUI Shell (`Shell` button)
2. run:
   - `genbank fetch AY738222 --as-id promega_luciferase_ay738222`
3. confirm the created sequence appears in project lineage/state
4. optional fallback when offline:
   - `File -> Open...` local file `data/tutorial_inputs/AY738222_promega_luciferase.gb`

Parity map:

| GUI action | Engine operation | CLI parity (optional) |
| --- | --- | --- |
| Fetch AY738222 in GUI Shell | `FetchGenBankAccession` | `gentle_cli shell 'genbank fetch AY738222 --as-id promega_luciferase_ay738222'` |
| Local fallback file import | `LoadFile` | `gentle_cli op '{"LoadFile":{"path":"data/tutorial_inputs/AY738222_promega_luciferase.gb","as_id":"promega_luciferase_ay738222"}}' --confirm` |

### Step 4: Create assembly preview (GUI routine path recommended)

Recommended GUI path:

1. `Patterns -> Routine Assistant...`
2. choose `Gibson Two-Fragment Overlap Preview`
3. bind:
   - `left_seq_id = tp73_promoter_fragment_1`
   - `right_seq_id = promega_luciferase_ay738222`
   - `overlap_bp = 20` (or your design)
4. run preflight (`validate only`)
5. run transactional execution

Parity map:

| GUI action | Engine route | CLI parity (optional) |
| --- | --- | --- |
| Gibson routine assistant | `macros template-run ... --validate-only` then transactional run | `gentle_cli shell 'macros template-run gibson_two_fragment_overlap_preview --bind left_seq_id=tp73_promoter_fragment_1 --bind right_seq_id=promega_luciferase_ay738222 --bind overlap_bp=20 --bind assembly_prefix=tp73_luc_demo --bind output_id=tp73_luc_construct_preview --validate-only'` |

Note:

- The canonical JSON skeleton currently uses a deterministic `Ligation + Branch`
  preview to keep operation-level parity explicit. For practical planning, the
  GUI routine-assistant Gibson path is preferred.

### Step 5: Design verification primers and qPCR assays

GUI:

1. open construct preview sequence
2. Engine Ops -> `Primer and qPCR design reports`
3. set ROI around promoter/vector junction
4. run `Design Primer Pairs`
5. run `Design qPCR Assays`
6. inspect and export reports

Parity map:

| GUI action | Engine operation | CLI parity (optional) |
| --- | --- | --- |
| Primer pairs | `DesignPrimerPairs` | `gentle_cli primers design @request.design_primer_pairs.json` |
| qPCR assays | `DesignQpcrAssays` | `gentle_cli primers design-qpcr @request.design_qpcr.json` |
| Report export | report store helpers | `gentle_cli primers export-report REPORT_ID out.json` / `gentle_cli primers export-qpcr-report REPORT_ID out.json` |

### Step 6: Package results for communication

GUI:

1. verify lineage graph nodes for:
   - TP73 anchor/extension
   - extracted promoter candidate(s)
   - imported AY738222 vector
   - assembly preview output
2. export sequence/SVG artifacts
3. save project `.gentle.json`

Parity map:

| GUI action | Engine/adapter route | CLI parity (optional) |
| --- | --- | --- |
| Save and export artifacts | shared engine export operations | `gentle_cli save` + `gentle_cli render-svg ...` + report-export commands |

## Recommended success criteria

1. You can regenerate one promoter candidate deterministically from the same
   anchor/extraction parameters.
2. The luciferase vector input ID is stable and referenced in lineage.
3. A single construct-preview ID is used consistently in primer and qPCR
   reports.
4. GUI and CLI routes describe the same operation semantics for each step.

## Notes on current limits

1. This tutorial is intentionally planning-focused.
2. Detailed wet-lab protocol logic (for example exact vector chemistry,
   transfection, luminometer analysis) remains user-owned.
3. Use promoter-track overlays (BED/BigWig/VCF) before final fragment selection
   when regulatory confidence is critical.
