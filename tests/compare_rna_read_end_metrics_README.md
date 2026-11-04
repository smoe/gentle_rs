# Read-End Metric Comparator Tests

Origin: hand-crafted synthetic scalar rows in
`test_compare_rna_read_end_metrics.py`, adapted from Glen / `smoe-bot` commit
`a1cf936d0cbccbbc6f1e2cbed8474006d740bc25`. No sequencing reads, measured UTR
lengths, transcript models or private barcode files are included.

Recreate the original fixture with `write_fixture`: six sample identifiers,
three metric labels, four observed values based on 100/200/300/400 plus literal
sample/metric offsets, and one blank missing value per cell. This creates 45
descriptive pair panels, not biological replicate evidence. The smaller
`minimal_input` fixture has two unrelated sample IDs, one arbitrary metric,
two observed values and one row of each censored/missing status. Dedicated
controls substitute observed zero, malformed quoting or ragged rows; LF and
CRLF are generated explicitly in temporary files.

These fixtures exercise only `scripts/compare_rna_read_end_metrics.py`: parser
and quantile semantics, source-byte hashes, output refusal/publication, file
closure and command-line paths. Publication faults and racing writes are
injected at the actual `os.link` boundary in temporary directories. A host
without symlink permission reports the symlink cases as skips, not acceptance.
Nothing loads a GENtle project or infers transcript geometry.

Run from the repository root:

```sh
python3 -m unittest discover -s tests -p test_compare_rna_read_end_metrics.py -v
```
