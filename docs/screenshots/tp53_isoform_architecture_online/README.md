# TP53 isoform-architecture GUI evidence

These PNGs are direct Linux/X11 captures from a disposable GENtle state built
from `docs/figures/tp53_ensembl116_panel_source.gb` and
`assets/panels/tp53_isoforms_v1.json`.

- `01-isoform-panel-controls.png`: 840×602 native Sequence tools view with the
  retained source open and panel path/id populated; SHA-256
  `09125b0ffc2d6c68465110392e24f4858fd83428472c3fff0d7531144390f48f`
- `02-isoform-expert.png`: 784×490 native Isoform Expert with seven mapped
  TP53 transcript rows and seven protein-domain rows; SHA-256
  `375038e41c338800c1cb965b7de33a091c89224bb7f70c0a7726cc160c1b5d63`
- source base: `4e7adf6e9dcba577892b95439af3ec749b272c3d`
- review-branch base: `7c49ff04856b9bc728149b783880fb8d643e87a5`
- GUI binary revision shown at capture: `v0.1.1790819522996`
- environment: isolated 1920×1080 Xvfb display; child windows captured only

The disposable state imported the panel with `strict=true` through the shared
shell before the GUI reopened it. The unchecked **Strict** checkbox in the
first screenshot is the current input control for a future import, not a
read-back receipt for the existing panel. The structured inspect result is the
acceptance evidence: gene `TP53`, panel `tp53_isoforms_v1`, seven transcript
lanes, seven protein lanes and no warnings.

The current deterministic SVG is
`docs/figures/tp53_isoform_architecture.svg`, 1200×1400, SHA-256
`144e3af1f4c8a2ddd09d6691c5972f8a308d5785eb65014acf671ac4c8f9c6ca`.
It is stronger renderer evidence than a screenshot because it is the portable
shared-engine artifact. Neither the local state nor these captures execute or
accept the online Ensembl genome preparation/extraction step.
