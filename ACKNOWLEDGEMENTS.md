# Acknowledgements

The GENtle project thanks external tool authors and communities whose work helps
shape our roadmap and user experience priorities.

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

Current status in GENtle:

- No PCRtools source code is currently vendored or imported directly.
- Feature ideas are tracked in `docs/roadmap.md` and implemented through
  GENtle's shared Rust engine contracts for deterministic GUI/CLI/adapter
  parity.
- This paper-first intake approach avoids treating repository code as the
  authoritative specification and helps prevent GPL-related reuse pitfalls.
