# Stateless sequence-inspection GUI evidence

This directory retains publication-safe manual/hybrid evidence for tutorial
02.02. The capture used only the committed synthetic 46 nt sequence and an
exact-revision debug build of the native Linux GUI.

- `01-tfbs-settings.raw.png` shows the selected `SP1,TP73` motifs,
  `min_llr_quantile=0.95`, `llr_background_tail_log10`, and disabled negative
  clipping.
- `02-direct-scan-inspectors.raw.png` shows the shared restriction-site report
  (EcoRI, SmaI and BamHI) and three retained SP1 hits. TP73 was included in the
  two-motif scan but produced no retained row at this threshold.

The whole-sequence score-track action was then attempted twice in independent
fresh Xvfb processes. Both processes aborted with exit code 134 and the same
native `thread '<unknown>' has overflowed its stack` message. The canonical
offline workflow still passed and produced its JSON/SVG artifacts. The images
therefore support only the setup and RE/TFBS-inspector claims; they do not claim
score-track GUI acceptance or biological validation.

See [`evidence.json`](evidence.json) for hashes and scope.
