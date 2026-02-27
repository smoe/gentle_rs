# Cross-Tool Cloning Gap Matrix (GENtle vs Serial Cloner vs MacVector)

Last updated: 2026-02-27

Purpose: capture feature gaps between GENtle and commonly used cloning tools,
with emphasis on features that recur across more than one external tool.

## 1. Inputs and scope

Primary inputs used for this comparison:

- local historical feature log:
  - `Serial-Cloner-History.txt`
- Serial Cloner homepage:
  - <http://serialbasics.free.fr/Serial_Cloner.html>
- MacVector functionality/workflow pages:
  - <https://macvector.com/functionality/>
  - <https://macvector.com/functionality/gateway-and-topo-cloning/>
  - <https://macvector.com/functionality/gateway-lic-assembly/>
  - <https://macvector.com/functionality/primer-design/>
  - <https://macvector.com/functionality/auto-annotation/>
  - <https://macvector.com/functionality/database-searching/>
  - <https://macvector.com/functionality/sequence-assembly/>
  - <https://macvector.com/functionality/protein-analysis/>
  - <https://macvector.com/functionality/agarose-gel-simulation/>
  - <https://macvector.com/getting-started/useful-workflows/>
- current GENtle status:
  - `docs/roadmap.md`
  - `docs/protocol.md`
  - `docs/cli.md`

This is a product-capability gap matrix, not a wet-lab protocol document.

## 2. Gap matrix

| Capability area | Serial Cloner | MacVector | GENtle current status | Gap class |
|---|---|---|---|---|
| Protocol-guided cloning families (Gateway/TOPO/Gibson/LIC/etc.) | Gateway BP/LR workflow and ligation/build workflows | Dedicated Gateway/TOPO and Gibson/LIC workflows | Starter macro templates only; full protocol packs still incomplete | High (already on roadmap) |
| Interactive cloning workspace (clipboard/project-style fragment flow with end-compatibility confirmation) | Build-a-construct flow, ligate fragments, adaptor/shRNA-assisted insertion | Cloning Clipboard + ligation compatibility UI | Engine operations exist; no equivalent interactive fragment project/clipboard workspace yet | High |
| Primer design and validation suite | PCR assistant workflows and primer overlap handling | Primer3/QuickTest, pair testing, primer interaction checks, primer DB workflows | PCR/PcrAdvanced/PcrMutagenesis are implemented, but no first-class design/test/primer-inventory workflow | High |
| Auto-annotation / scan for missing features from curated vector collections | Feature scan and user feature collections/import/export | Auto-annotation against local annotated folders | No dedicated "scan missing features from library" workflow is documented | High |
| Sequence confirmation from reads (trace/NGS), contig editing/assembly | Local align and BLAST2Seq integration (legacy scale) | Full assembly/align-to-reference and chromatogram/contig tooling | Not present as a first-class sequencing-confirmation workflow | High |
| Database search/import UX | Built-in web access and direct parse of NCBI/EMBL pages | Entrez keyword search + integrated BLAST retrieval/download | BLAST helper/reference flows exist, but no equivalent Entrez keyword retrieval + one-click import path | Medium |
| Protein-centric analysis and reverse translation workflows | Protein windows, protein alignment, reverse translation, pI estimation | Broad protein analysis toolbox + reverse translation | DNA-centric baseline; no comparable protein analysis suite in current roadmap scope | Medium |
| CRISPR wet-lab screening flow | Not a major explicit focus in listed history | PAM scan and INDEL screening workflow surfaced in workflows | Guide design is incomplete and screening/confirmation flow is not complete | High (partly on roadmap) |
| Broad proprietary format interoperability for cloning files | Imports MacVector/VectorNTI/DNAStar/ApE/pDRAW formats | Native ecosystem + broad import/export tooling | GenBank/EMBL/FASTA/XML-first; proprietary cloning-file compatibility is limited | Medium |
| Realistic gel simulation controls | Virtual cutter and restriction-oriented simulation | Document-based agarose simulation with richer realism controls | Virtual gel rendering exists; realism and arrangement authoring remain partial | Medium (already on roadmap) |

## 3. Repeated gaps across multiple external tools

These are missing features that appear in at least two of:
Serial Cloner, MacVector, and previously-tracked SnapGene parity work.

1. Protocol-family cloning assistants beyond primitive ops.
2. Interactive primer design/test workflows (including pair interaction checks).
3. Auto-annotation from reusable feature/vector libraries.
4. Sequence confirmation workflow from Sanger/NGS reads.
5. Interactive cloning workspace abstractions (clipboard/project/fragment flow).
6. CRISPR end-to-end workflow completion (design plus screening/confirmation).

These repeated gaps should be prioritized above one-off tool-specific features.

## 4. MacVector-specific differences (not all immediate priorities)

Compared to current GENtle focus, MacVector also emphasizes:

- deep protein analysis toolbox breadth,
- broad sequence assembly stack (de novo + reference + long-read adapters),
- extensive database search UX integrated into desktop workflows.

These are important for long-term competitiveness, but they can follow the
cloning-routine and sequencing-confirmation priorities above.
