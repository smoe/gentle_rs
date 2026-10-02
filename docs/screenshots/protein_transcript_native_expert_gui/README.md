# Transcript-native Protein Expert GUI evidence

Both PNGs are direct captures of the native Linux/X11 GUI with
`test_files/tp73.project.gentle.json` open.

- `01-filter-ready.png`: the 728×476 `Protein Evidence` controls with local
  sequence `tp73.ncbi` and filter `NM_005427.4`; SHA-256
  `781643242467da725e4e41b8e6c5b33f463c5f665149bec1a26b42094d3aecbd`
- `02-filtered-derived-expert.png`: the 784×490 filtered native expert with the
  corrected TP73 title, one 636-aa derived-only row and the visible
  missing-CDS/inferred-ORF warning; SHA-256
  `075b1ad07632d43d6d4a8a66efb7936a0bb7741f178fdf21eb82fff29a4a4ff0`
- source base: `4e7adf6e9dcba577892b95439af3ec749b272c3d`
- review-branch base: `f82d59642ac57c49eae413516508cc36ce4553a4`
- GUI binary revision shown at capture: `v0.1.1790819522996`; it includes the
  then-uncommitted transcript-lane title fix
- environment: isolated 1920×1080 Xvfb display; child windows captured only

The first unpatched expert capture named an unrelated gene-like annotation
(`ATAC-STARR-seq lymphoblastoid silent region 122`) even though every displayed
transcript row belonged to TP73. The review patch makes the title prefer the
gene metadata of the transcript feature referenced by the displayed lanes. A
focused regression test places `UNRELATED_REGION` before a `TPROT` transcript
and requires the expert gene symbol to remain `TPROT`.

The GUI expert path and filter were accepted. Clicking `Render Derived Protein
SVG...` did not expose a visible save dialog or create a file in this Xvfb
session, so these images do not claim native file-dialog/export acceptance.
The shared-shell route rendered the same filtered payload to a 1200×880 SVG
with SHA-256
`70149edbbc5daf0f50f22bb34eb29ec3017f5a2ee11d7775f611262ac3653046` and
reported no changed or created sequence ids.
