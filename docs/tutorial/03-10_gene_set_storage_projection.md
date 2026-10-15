# Place A Logical Gene Set Into A Storage Rack

> Type: `CLI / GUI Shell walkthrough`
> Status: `manual/hybrid; executable synthetic starter; Rust replay pending`

A gene set lists biological subjects, not physical tubes. An arrangement orders
explicitly chosen tubes; a rack places them at A1-style positions. This minimal
path keeps one tube per member. It does **not** mix or clone its contents, infer
gene-to-material identity, or prove that a sample exists in the laboratory.

## 1. Create A Disposable Synthetic Starter

From the repository root, choose a fresh state filename:

```bash
cargo run --bin gentle_cli -- --state storage-review.gentle.json \
  workflow @docs/examples/workflows/gene_set_storage_starter.json
```

The two literal eight-base strings represent no natural gene. A fresh starter
creates `storage_source_a`, `storage_source_b`, singleton `container-1` and
`container-2`, and unannotated logical `GENEA`/`GENEB` without a genome index.
Read the actual resolution `op_id` in the JSON (`op-3` in a fresh starter).
These are occurrence IDs, not names to reuse in an existing project. Inspect
its actual resolution and tubes first; a label never establishes the binding.

## 2. Review And Explicitly Authorize Placement

```bash
cargo run --bin gentle_cli -- --state storage-review.gentle.json op \
  '{"CreateGeneSetArrangement":{"resolution_id":"op-3","member_bindings":[{"stable_member_id":"symbol:genea","source_container_id":"container-1"},{"stable_member_id":"symbol:geneb","source_container_id":"container-2"}],"arrangement_id":"storage_review","name":"Synthetic storage review","confirm_physical_placement":true}}'
```

Omitted/false confirmation refuses. Incomplete resolutions, missing/extra or
duplicate bindings, nonexclusive/nonsingleton tubes, missing sequences and reused
arrangement IDs refuse before layout creation. This slice accepts 1..24 members
on the existing 4-by-6 small-tube rack; resolved-set order determines lane order.

Inspect `result.collection_operation`: `lifting_mode` is `arrange`, with exact
member/container/sequence-snapshot provenance. Persisted
`gene_set_arrangement:storage_review` metadata links the source-resolution digest,
report, arrangement and rack ID. The set, DNA and source tubes remain unchanged.
The existing serial creator adds its linked rack draft (normally `rack-1` in
this fresh project). This is not a physical inventory or a pool.

## 3. Inspect And Export The Existing Rack

Substitute the actual reported rack ID if it differs:

```bash
cargo run --bin gentle_cli -- --state storage-review.gentle.json racks show rack-1
cargo run --bin gentle_cli -- --state storage-review.gentle.json racks labels-svg rack-1 storage-review.labels.svg
cargo run --bin gentle_cli -- --state storage-review.gentle.json racks isometric-svg rack-1 storage-review.isometric.svg --template storage_pcr_tube_rack
cargo run --bin gentle_cli -- --state storage-review.gentle.json racks openscad rack-1 storage-review.scad --template storage_pcr_tube_rack
```

Check A1/A2 and member/tube order in all outputs. Labels, SVG and OpenSCAD consume
the same layout; SCAD does not emboss the sample labels. The existing
[24-slot reference SCAD](../figures/gibson_single_insert_storage_rack.scad)
supplies the carrier geometry. On macOS, prefer the already-installed app:

```bash
/Applications/OpenSCAD.app/Contents/MacOS/OpenSCAD \
  -o storage-review.png --imgsize=1000,750 --viewall --autocenter --render=true \
  storage-review.scad
```

The reference rendered locally after its floor repair. **The new operation and
its Rust exports have not run locally**: Glen/CI must replay and render the fresh
output. Rendering is not fit, manufacturing, material or sterility validation.

## 4. GUI, Agents And Verification

GUI Shell accepts the same `op` JSON. Then use Arrangements → Open Rack in the
existing workspace. There is no new dedicated creation button. MCP uses its
existing `op` tool with `confirm=true` plus the operation's placement flag;
no curated tool or model call is added. Workflows use the identical operation.

Within one GUI session, Undo removes arrangement/rack/projection metadata, not
the source set or tubes; Redo restores them together. Save/reopen retains the
record, not cross-process undo history. Run `gene_set_arrangement_`, rack-export,
capability/parity and workflow checks externally; source tests are not acceptance.

The starter is executable; the walkthrough stays manual/hybrid until external
runtime replay and generated tutorial registration. Freezer boxes, physical
inventory, barcode workflows and new GUI controls remain explicit follow-ups.
