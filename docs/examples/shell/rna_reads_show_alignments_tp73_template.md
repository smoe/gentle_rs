# RNA-Read Alignment Batch Inspection Shell Template

This shell-only template demonstrates the read-only
`rna-reads show-alignments` route for reviewing exact phase-2 pairwise
alignment columns from a saved RNA-read report. It is intentionally not a
`docs/examples/workflows/*.json` file because `show-alignments` follows the
existing `show-alignment` precedent: it is a shared-shell inspector rather than
a mutating or persisted `Operation`.

The route returns `gentle.rna_read_alignment_display_batch.v1`. Each
`entries[]` item contains the same `RnaReadAlignmentDisplay` payload that a
single `rna-reads show-alignment REPORT_ID RECORD_INDEX` call would produce.

## Direct CLI

```bash
STATE=/home/clawbio/work/tp73_pancreas_benchmark/tp73_ensembl_pancreas.gentle.json
REPORT_ID=tp73_terminal_SRR32957126
GENE_ID=TP73

cargo run --quiet --bin gentle_cli -- --state "$STATE" \
  rna-reads show-alignments "$REPORT_ID" \
  --gene "$GENE_ID" --cohort all --complete-rule near --limit 5 \
  --output "$REPORT_ID.${GENE_ID}.alignment_batch.top5.json"
```

Quick shape check:

```bash
jq '{
  schema,
  report_id,
  entry_count,
  skipped_count: (.skipped_records | length),
  first_record_index: (.entries[0].record_index // null),
  first_aligned_query_bp: ((.entries[0].alignment.aligned_query // "") | length)
}' "$REPORT_ID.${GENE_ID}.alignment_batch.top5.json"
```

## Shared Shell

```bash
cargo run --quiet --bin gentle_cli -- --state "$STATE" shell \
  'rna-reads show-alignments tp73_terminal_SRR32957126 --gene TP73 --cohort all --complete-rule near --limit 5'
```

Use the same command text in GUI Shell, JS/Lua/Python shell helpers, and
MCP/agent routes that accept reviewed shared-shell commands.

## Python Wrapper

```python
import json
from pathlib import Path
import sys

STATE = "/home/clawbio/work/tp73_pancreas_benchmark/tp73_ensembl_pancreas.gentle.json"
sys.path.insert(0, str(Path("integrations/python").resolve()))
from gentle_py import GentleClient

client = GentleClient(state_path=STATE)
batch = client.shell(
    "rna-reads show-alignments "
    "tp73_terminal_SRR32957126 "
    "--gene TP73 --cohort all --complete-rule near --limit 5",
    expect_json=True,
)
print(json.dumps({
    "schema": batch["schema"],
    "entry_count": batch["entry_count"],
    "first_aligned_query_bp": len(
        batch["entries"][0]["alignment"]["aligned_query"]
    ) if batch["entries"] else 0,
}, indent=2))
```
