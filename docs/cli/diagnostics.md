# Latency Diagnostics

<!-- manual-counterparts:start -->
Related GUI reference: [Workspace Latency Diagnostics](../gui/workspace-latency-diagnostics.md).
These topic links do not imply identical feature coverage across interfaces.
<!-- manual-counterparts:end -->

Optional GUI diagnostics when launching `gentle` from a terminal:

- `GENTLE_TRACE_TOPOLOGY_TRANSITIONS=1`
  - enables coarse circular/linear topology-switch trace logging from sequence
    windows
  - writes stage summaries to stderr and mirrors them into the sequence-window
    status line while enabled
  - useful when a topology switch appears to hang or repaint-loop
- `GENTLE_TRACE_SPLICING_EXPERT=1`
  - enables coarse Splicing Expert open/render trace logging from sequence
    windows
  - writes stage summaries to stderr and mirrors them into the sequence-window
    status line while enabled
  - useful when opening the Splicing Expert appears to hang or stall
- `GENTLE_TRACE_RNA_READ_MAPPING=1`
  - enables coarse RNA-read Mapping open/render trace logging from sequence
    windows
  - writes stage summaries to stderr and mirrors them into the sequence-window
    status line while enabled
  - useful when opening the RNA-read Mapping workspace appears to hang or stall
- `GENTLE_TRACE_PROMOTER_DESIGN=1`
  - enables coarse Promoter design open/render trace logging from sequence
    windows
  - writes stage summaries to stderr and mirrors them into the sequence-window
    status line while enabled
  - useful when opening the Promoter design workspace appears to hang or stall

Architecture invariant: all adapters/frontends above route cloning/business
behavior through the same shared engine.

<!-- manual-page-nav:start -->

---

[CLI Manual](../cli.md) | [Previous: Engine Settings](settings.md) | [Next: Screenshot Capture Policy](screenshot-capture.md)
<!-- manual-page-nav:end -->
