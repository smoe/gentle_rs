# Acknowledgements

The GENtle project thanks external tool authors and communities whose work helps
shape our roadmap and user experience priorities.

## Primer3 Method Reference

GENtle's primer-design backend can call the external Primer3 executable, and
GENtle's internal oligo-QC report vocabulary follows Primer3's public
distinction between broad oligo self/pair complementarity and 3'-end anchored
complementarity (`SELF_ANY`, `SELF_END`, `COMPL_ANY`, `COMPL_END`, plus the
thermodynamic `_TH` variants).

Primary references:

- Untergasser A, Cutcutache I, Koressaar T, Ye J, Faircloth BC, Remm M, Rozen
  SG. 2012. *Primer3--new capabilities and interfaces*. *Nucleic Acids
  Research* 40(15):e115.
- Koressaar T, Remm M. 2007. *Enhancements and modifications of primer design
  program Primer3*. *Bioinformatics* 23(10):1289-1291.
- Primer3 manual and source repository:
  - <https://primer3.org/manual.html>
  - <https://github.com/primer3-org/primer3>

Current status in GENtle:

- Primer3 source code is not vendored or translated into GENtle.
- The cDNA PCR/qPCR oligo-QC exact-run screen is an independent Rust
  implementation that credits Primer3's public method vocabulary and keeps
  future thermodynamic Primer3-backed checks as an optional extension.

## PCRtools Stimulus

We thank **Ruslan Kalendar** for useful stimulus in PCR-related feature planning,
especially around assay-family breadth and user workflow expectations.

Primary information source for GENtle roadmap intake:

- Peer-reviewed article (main source):
  - Kalendar R. 2025. *Comprehensive web-based platform for advanced PCR design,
    genotyping, synthetic biology, molecular diagnostics, and sequence
    analysis*. *Molecular Therapy Nucleic Acids* 36(4):102716.
  - DOI: <https://doi.org/10.1016/j.omtn.2025.102716>
  - Article page: <https://www.sciencedirect.com/science/article/pii/S2162253125002707>
  - Open-access full text (PMC): <https://pmc.ncbi.nlm.nih.gov/articles/PMC12506483/>

Secondary context source:

- Repository: <https://github.com/rkalendar/PCRtools> (consulted as supporting
  implementation context, not as primary specification source)

Additional method/figure-style reference for overlap-extension insertion/deletion mutagenesis:

- Lee J, Shin M-K, Ryu D-K, Kim S, Ryu W-S. 2010. *Insertion and Deletion
  Mutagenesis by Overlap Extension PCR*. In: Braman J (ed), *In Vitro
  Mutagenesis Protocols*. Methods in Molecular Biology, vol. 634. Humana Press.
  - DOI: <https://doi.org/10.1007/978-1-60761-652-8_10>
  - Protocol page: <https://link.springer.com/protocol/10.1007/978-1-60761-652-8_10>

Current status in GENtle:

- No PCRtools source code is currently vendored or imported directly.
- Feature ideas are tracked in `docs/roadmap.md` and implemented through
  GENtle's shared Rust engine contracts for deterministic GUI/CLI/adapter
  parity.
- This paper-first intake approach avoids treating repository code as the
  authoritative specification and helps prevent GPL-related reuse pitfalls.
