# Inspecting TATA-Box Evidence

Start in the DNA viewer's **TFBS scan > TATA-box evidence** window. Inspect
the loaded sequence or selection, click a result's coordinates to navigate,
and optionally **Add selected to DNA map**. Adding features is explicit and
undoable. Inspection and file reading run in the background. No EPD download
or external program is required for local annotations and predictions.
The command palette also exposes **TATA-box Evidence**; it opens the current
sequence's workspace, or queues it while the sequence loads.

## Three Different Questions

| Evidence | What it tells you | What it does not establish |
|---|---|---|
| Source annotation | The loaded source labels this exact feature `TATA_box` or `TATA_signal`; qualifiers are retained. | That the annotation was experimentally validated. |
| EPDnew classification | A promoter has an evidence-backed TSS and an explicit positive/negative computational TATA classification. | The exact motif interval, occupancy, or autonomous promoter function. |
| GENtle TBP prediction | An exact versioned TBP PFM scores this oriented sequence above the chosen threshold. | Experimental binding or transcription initiation. |

EPD's [published classification rule](https://epd.expasy.org/epd/EPDnew_select.php)
uses FindM near -28 +/-3 bp from the TSS. Its BED track describes the promoter,
not the TATA motif. GENtle therefore shows that evidence at the TSS, separately
from exact source/prediction intervals. It does not multiply overlapping evidence
into an uncalibrated confidence score.

## Local Inspection

In the shared Shell, after loading a sequence with id `promoter`:

```text
promoters tata-screen '{"seq_id":"promoter"}' --output tata-report.json
```

The default uses exact transcript 5-prime boundaries as candidate TSSs, retains
each transcript, and scans motif starts from -40 to -15 bp on its strand.
Partial transcript geometry is not silently interpreted as an exact TSS.
If the sequence has no transcript annotations, provide a measured or explicitly
assumed TSS. For a first transcribed base at position 101 on the forward strand:

```text
promoters tata-screen '{"seq_id":"promoter","additional_tss":[{"id":"my_TSS","position_0based":100,"reverse":false,"source":"user-supplied TSS; specify experiment or assumption"}]}'
```

For exploratory scanning without known TSS context:

```text
promoters tata-screen '{"seq_id":"promoter","scan_without_tss":true,"minimum_llr_bits":6.0}'
```

Short AT-rich motifs are common. A hit outside a credible promoter context is
particularly weak evidence. The score threshold is not a p-value. The default
is [JASPAR TBP MA0108.3](https://jaspar.elixir.no/matrix/MA0108.3/), a 7-base
matrix; its version and actual matrix-content hash travel with every report.

## Optional EPDnew Inputs

For human GRCh38, download both files explicitly from the same release. These
commands are preparation, not part of the read-only screen:

```sh
mkdir -p data/resources/epd
curl --fail --location https://epd.expasy.org/ftp/epdnew/H_sapiens/006/Hs_EPDnew_006_hg38.bed -o data/resources/epd/Hs_EPDnew_006_hg38.bed
curl --fail --location https://epd.expasy.org/ftp/epdnew/H_sapiens/006/db/promoter_motifs.txt -o data/resources/epd/promoter_motifs.txt
```

Supply this EPD object in a screen request, or paste it in the viewer's optional
EPD source section. A button provides the same human template:

```json
{
  "bed_path": "data/resources/epd/Hs_EPDnew_006_hg38.bed",
  "motifs_path": "data/resources/epd/promoter_motifs.txt",
  "assembly": "GRCh38",
  "taxon_id": 9606,
  "release": "006",
  "source_url": "https://epd.expasy.org/ftp/epdnew/H_sapiens/006/",
  "required": true
}
```

Use source metadata appropriate to your actual download, never relabel an hg19
file as hg38. The files do not themselves encode their assembly; supplied source
metadata is an assertion whose compatibility is checked against the genome
anchor, not a sequence-level validation of the download. For repeat runs, set
`expected_bed_sha256` and `expected_motifs_sha256` from a reviewed report.
Only loaded, assembly-compatible chromosomes can be projected; no liftOver or
coordinate inference from gene names occurs. Missing motif-table IDs remain
unassessed. The adapter preserves all promoters, not just one promoter per gene.

The optional published-file parser check uses those same local downloads:

```sh
GENTLE_TEST_EPD_BED=data/resources/epd/Hs_EPDnew_006_hg38.bed GENTLE_TEST_EPD_MOTIFS=data/resources/epd/promoter_motifs.txt cargo test --locked --lib tata_epd_published_files_opt_in
```

Without these variables the test does not read external files. Synthetic tests
cover strand conversion, classification joins, missing data, and hash rejection.

## Keep And Reuse Results

`Copy report JSON` retains evidence, coordinates, model/file hashes, warnings
and the interpretation caveat. `Copy command` retains the exact input request.
To add selected rows headlessly, use:

```text
promoters tata-materialize @reviewed-tata-selection.json
```

The JSON object contains `screen` (the report's effective `request`),
`expected_report_sha256` (its `content_sha256`), and an explicit `row_ids` array.
This is a distinct mutation: changed sequence, annotations, anchor, EPD files,
model or parameters require another review. Added candidates keep their evidence
class and provenance as feature qualifiers and use the existing regulatory
feature layer in the DNA map. EPD classification markers remain labelled TSSs.

Direct CLI uses the same command following `gentle_cli --state PROJECT`.
MCP uses the advertised `op` tool with `ScreenTataBoxes` or
`MaterializeTataBoxFeatures`; normal tool confirmation still applies. JS, Lua,
Python and workflows accept these same operation JSON payloads, with no separate
scoring path in the adapters. The inner agent discovers them through capabilities
and the glossary. No existing reports or annotations are rewritten by inspection.
