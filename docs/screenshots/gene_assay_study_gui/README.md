# Public PATZ1 gene-assay GUI checkpoints

This directory retains six manual/hybrid GUI checkpoints from the exact
`.11` code candidate `3c1c32bcceacbe7d327c09f83ba57dcca2682b06`.
The session used the pinned public Ensembl 116/RefSeq PATZ1 fixture, a fresh
HOME/XDG profile, Xvfb plus Openbox, and a network namespace with no interfaces
other than loopback. No private sample, inner agent, BLAST search, order, or
experimental acceptance was involved.

Each checkpoint has:

- an untouched 1920x1080 X11-root `*.raw.png`;
- the exact `gentle.gui_semantic_snapshot.v2` snapshot written beside it; and
- a deterministic `*.context.svg` crop with a red review rectangle.

The SVGs reference the adjacent raw PNGs and are teaching views, not separate
evidence. [`evidence.json`](evidence.json) binds the candidate, binaries,
public inputs, actual design result, environment, and every retained file by
SHA-256.

The actual Primer3 2.6.1 replay produced seven primer pairs for 13 exact
mature-cDNA classes and left nine class-pair distinctions unresolved. GENtle
therefore reports `partial`. Genomic and whole-transcriptome specificity were
not run, and the study/dossier remain pending. The captures demonstrate GUI
presentation of those facts; they do not establish biological suitability or
permission to order.

The session was manual/hybrid because native file dialogs, review state and
pair inspection are not fully scripted controls. Only reached checkpoints are
published here. Exploratory and failed captures remain outside the repository.

