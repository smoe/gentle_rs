# Tiny TP73 Isoform BLAST Database Source

`cds.fasta` contains three real human TP73 DeltaN coding sequences (including
their stop codons), in transcript 5'-to-3' orientation without introns or UTRs:

| FASTA ID | ENA transcript | CDS, 1-based inclusive | Bases |
| --- | --- | --- | ---: |
| `AY040827.1_DNp73alpha_CDS` | AY040827.1 | 235..1998 | 1764 |
| `AY040828.1_DNp73beta_CDS` | AY040828.1 | 235..1587 | 1353 |
| `AY040829.1_DNp73gamma_CDS` | AY040829.1 | 235..1515 | 1281 |

Total: 4,398 bases. This is a small search/index fixture, **not** all TP73
isoforms, a whole transcriptome, or evidence of primer specificity. Shared
sequence produces legitimate cross-isoform hits; only exact full-length self
matches are required by the smoke test. These are accession-bound ENA records,
not substitutions for current Ensembl transcripts.

## Origin and Recreation

The source tree already carries the ENA transcripts, retrieved on 2026-05-01,
in [`data/resources/tp73_dn_ena_transcripts.fasta`](../../../data/resources/tp73_dn_ena_transcripts.fasta).
Their CDS coordinates and original sequence hashes come from
[`assets/panels/tp73_dn_isoforms_v1.json`](../../../assets/panels/tp73_dn_isoforms_v1.json).
`manifest.json` records the exact accession URLs, coordinates, source and CDS
sequence hashes, and FASTA file hash. No sequence was invented or downloaded
while preparing this fixture.

From the repository root, this deterministic extraction prints the committed
FASTA verbatim (70-base lines, LF endings); it does not require BLAST or a network:

```sh
python3 - <<'PY'
import json
from pathlib import Path

fixture = Path("test_files/fixtures/blast_tp73_isoforms")
manifest = json.loads((fixture / "manifest.json").read_text())
sources = {}
for block in Path(manifest["source_fasta"]).read_text().split(">")[1:]:
    lines = block.splitlines()
    sources[lines[0].split()[0]] = "".join(
        line for line in lines[1:] if not line.startswith(";")
    )
for record in manifest["records"]:
    start, end = record["cds_start_1based"], record["cds_end_1based"]
    sequence = sources[record["accession"]][start - 1:end]
    print(f'>{record["id"]} source={record["accession"]} cds={start}..{end}')
    for offset in range(0, len(sequence), 70):
        print(sequence[offset:offset + 70])
PY
```

## Use in GENtle

The `workflow_examples` tests use this fixture for tool-free provenance/CDS
validation and an optional native BLAST+ build/inspection/search smoke:

```sh
cargo test --lib tp73_blast_fixture -- --nocapture
```

If any of `makeblastdb`, `blastdbcmd`, or `blastn` is absent, the native test
prints an explicit skip without installing anything. A broken installed tool
or missing expected self match fails the test. The independent provenance test
still runs. Native indexes are created in a temporary directory and removed
afterwards; platform/tool-specific BLAST files are deliberately not committed.

For manual exploration with already-installed BLAST+, from the repository root:

```sh
db_dir="$(mktemp -d)"
makeblastdb -in test_files/fixtures/blast_tp73_isoforms/cds.fasta \
  -dbtype nucl -parse_seqids -blastdb_version 5 -out "$db_dir/tp73"
blastdbcmd -db "$db_dir/tp73" -info
blastn -task megablast -dust no \
  -query test_files/fixtures/blast_tp73_isoforms/cds.fasta \
  -db "$db_dir/tp73" -outfmt '6 qseqid sseqid pident length qstart qend sstart send'
```

The database should report three sequences. Each query must have a 100% identity
self hit spanning its complete CDS; extra cross-isoform hits are expected.
The synthetic conservation tutorial retains its separate toy-genome fixtures
because TP73 isoforms are not cross-species orthology evidence.
