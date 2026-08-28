# Inspect Cryptic-Splicing Hypotheses Without Calling Them Events

> Type: `CLI + GUI walkthrough`
> Status: `manual/hybrid`
> Data: a tiny, committed, synthetic sequence; fully offline

An intronless expression cassette can contain `GT...AG` arrangements that look
intronic enough to merit follow-up. Finding such an arrangement does **not**
show that the spliceosome uses it. This tutorial teaches GENtle's four separate
levels of information:

1. sequence structure (`structural_candidate`),
2. optional model compatibility (`model_scored`),
3. observed RNA junction support (`observed_junction` in a separate overlay),
4. projected protein-feature context (also a separate report).

The fixture was deliberately designed to produce a clear structural candidate.
It is teaching data, not a benchmark and not evidence about any real construct.

## 1. Load the synthetic cassette

From the repository root:

```bash
STATE=/tmp/gentle-cryptic-splicing-demo.gentle.json
REQUEST=/tmp/gentle-cryptic-splicing-request.json
OVERLAY=/tmp/gentle-cryptic-splicing-overlay-request.json
REPORT=/tmp/gentle-cryptic-splicing-report.json
SVG=/tmp/gentle-cryptic-splicing-report.svg

cargo run --bin gentle_cli -- --state "$STATE" \
  op '{"LoadFile":{"path":"docs/tutorial/inputs/cryptic_splicing_demo.fa","as_id":"cryptic_splicing_demo"}}' \
  --confirm
```

The input's provenance and exact sequence are recorded in
[`docs/tutorial/inputs/README.md`](./inputs/README.md).

## 2. Define one structural-only request

Create `$REQUEST` with this JSON:

```json
{
  "seq_id": "cryptic_splicing_demo",
  "start_1based": 1,
  "end_1based": 40,
  "strand": "forward",
  "model_policy": "structural_only",
  "min_pseudo_intron_bp": 20,
  "max_pseudo_intron_bp": 100,
  "max_candidate_pairs": 20
}
```

Coordinates are 1-based inclusive. There is no CDS annotation in this FASTA,
so coding consequences should remain unavailable rather than being inferred.

Run the portable inspection route:

```bash
cargo run --bin gentle_cli -- --state "$STATE" \
  shell "splicing cryptic-screen @$REQUEST"
```

Inspect these fields:

- `model.status` is `absent`, because the request explicitly chose the
  dependency-free structural policy;
- candidate `evidence_class` is `structural_candidate`, not
  `observed_junction`;
- branchpoint/polypyrimidine fields explain what was evaluable;
- `budget` distinguishes admissible, evaluated, and reported pairs;
- `source_digest`, `request_sha256`, and `effective_input_sha256` bind what was
  actually inspected.

## 3. Export the same typed result

```bash
cargo run --bin gentle_cli -- --state "$STATE" \
  shell "splicing cryptic-export @$REQUEST $REPORT"

cargo run --bin gentle_cli -- --state "$STATE" \
  shell "splicing cryptic-render @$REQUEST $SVG"
```

The JSON is the machine-readable authority. The SVG is a deterministic
projection of that screen, not an independently computed figure.

## 4. Make missing RNA evidence explicit

Create `$OVERLAY` with:

```json
{
  "screen_request": {
    "seq_id": "cryptic_splicing_demo",
    "start_1based": 1,
    "end_1based": 40,
    "strand": "forward",
    "model_policy": "structural_only",
    "min_pseudo_intron_bp": 20,
    "max_pseudo_intron_bp": 100,
    "max_candidate_pairs": 20
  },
  "rna_read_report_ids": [],
  "nearby_tolerance_bp": 2
}
```

Then run:

```bash
cargo run --bin gentle_cli -- --state "$STATE" \
  shell "splicing cryptic-overlay @$OVERLAY"
```

Every candidate is `not_assessable`. That is the correct result: no RNA report
was supplied. When real saved RNA reports are named, GENtle accepts them only
if sequence id, full source-sequence SHA-256, coordinate space, and coordinate
convention match. Junction boundaries are compared in normalized ascending
source coordinates, so reverse-strand candidates use the same exact geometry
test. A bound report with no retained matching junction yields
`no_retained_junction_support`, which is still **not evidence of biological
absence**.

## 5. Optional MaxEnt scoring

GENtle does not download or redistribute MaxEntScan probability tables. If you
have independently obtained and reviewed a local MaxEntScan directory or ZIP,
normalize it with explicit provenance:

```bash
cargo run --bin gentle_cli -- shell \
  'resources sync-maxent /path/to/local/MaxEntScan \
   data/resources/maxent_splice_model.json \
   --source-url YOUR_REVIEWED_SOURCE_URL \
   --retrieved-on YYYY-MM-DD \
   --redistribution-status user_supplied_not_redistributed'

cargo run --bin gentle_cli -- shell 'resources status'
```

Change `model_policy` to `use_active_maxent`; for a reproducible run, also copy
the reported model fingerprint into `expected_model_fingerprint_sha256`.
Separate donor and acceptor scores are the model outputs. Their sum is labelled
only as a prioritization heuristic. If the output budget truncates a scored
search, GENtle evaluates every admissible pair and returns a deterministic
global top set; `ranking_complete` records that distinction.

## 6. Inspect the same request in the GUI

1. Open the state/project and the `cryptic_splicing_demo` DNA window.
2. Open `Sequence Tools -> Cryptic screen`, or use the `Cryptic screen` section
   in Splicing Expert.
3. Enter span `1..40`, pseudo-intron bounds `20..100`, and budget `20`.
4. Leave MaxEnt off for the baseline run, then click `Run structural screen`.
5. Compare candidate coordinates, evidence labels, budget accounting, hashes,
   and warnings with the CLI report.
6. Use `Copy report JSON`, `Save JSON...`, or `Render SVG`; all dispatch the
   shared portable operations.

The GUI runs the operation on an immutable engine snapshot. It does not persist
candidate rows or mutate the project. If an active model is replaced, an old
cached report remains copyable as historical JSON, while recomputed exports and
evidence joins require a fresh screen.

## 7. Protein-feature context, when available

For a real annotated sequence, first create and review a stored UniProt genome
projection. Then provide its id through:

```text
splicing cryptic-protein REQUEST_JSON_OR_@FILE [--output OUTPUT.json]
```

The result says whether the **hypothetical removed nucleotide interval**
overlaps digest-bound projected protein-feature geometry. It does not say that
splicing occurs, that a protein product exists, or that a domain is
experimentally lost. Legacy projections without a source-sequence digest are
reported as `not_assessable`.

## Interpretation checklist

- A `GT...AG` pair is a structural candidate, not an event.
- A MaxEnt score is model compatibility, not a probability of splicing.
- `not_evaluable` is different from a negative signal.
- Nearby RNA support is not exact RNA support.
- No retained RNA support is not biological absence.
- A projected protein-feature overlap is annotation context over a hypothetical
  removal, not a biological verdict.
- Save the request, effective-input hash, model fingerprint, and evidence report
  hashes with any downstream interpretation.
