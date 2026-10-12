# PATZ1 Gene-Assay Study GUI Acceptance — 2026-09-20

## Verdict

The six requested public-reference GUI checkpoints were reached manually at the
exact `.11` code candidate
`3c1c32bcceacbe7d327c09f83ba57dcca2682b06`. The captures are suitable as
tutorial illustrations of the observed state. This is **manual/hybrid GUI
acceptance**, not fully automated acceptance, biological validation, reference-
wide specificity, study execution or order approval.

The exact GUI binary SHA-256 is
`8e32e19f1cfc767d32c01c3184dd55aa045b3fc0d3c12d30ba4bd698aa4225ae`;
the exact release CLI SHA-256 is
`09439a09318184b7f7fbf9b64d8450e9876fec9e450a43f1230dd6495ed4cbd2`.
The GUI ran under Xvfb/Openbox with fresh HOME/XDG directories and a Linux user
plus network namespace. The pinned public fixture manifest SHA-256 is
`6f96dbdda26b34dbf8456559fa7dc66f7fca2802f423604d9e5318818e833c21`.

## Checkpoints

- **G1 — pass:** the loaded locus report presented 13 Ensembl transcript rows,
  the four RefSeq comparison records and the exact shared
  `ENST00000266269.10` / `NM_014323.3` chain without adding RefSeq to the design
  universe. The ordinary DNA map remained bound to the 20,802-bp minus-strand
  PATZ1 locus.
- **G2 — pass:** the saved request exposed the 13-transcript/all-annotated-cDNA
  universe, minimal-discrimination objective, best-effort coverage policy and
  constraints. RefSeq remained comparison context.
- **G3 — pass with a partial scientific result:** Primer3 2.6.1 produced seven
  selected pairs from 37 candidates. The GUI showed both 5'-to-3' primer
  sequences, the 13-by-7 product matrix and the warning that nine class-pair
  distinctions remain unresolved. All 13 exact mature-cDNA classes had at least
  one predicted product; that is coverage, not complete discrimination.
- **G4 — pass:** the separately loaded plan retained plan ID
  `patz1_real_study`, Ensembl 116 provenance and missing declared evidence. The
  persisted exploratory panel was presented as comparison context, not as an
  executed study result.
- **G5 — pass as a blocker:** the selected pair remained explicitly
  `Specificity not assessed`; all seven genomic specificity follow-ups were
  `not_run`, and no whole-transcriptome result was supplied. Readiness and
  ordering therefore remain blocked.
- **G6 — pass as pending:** canonical dossier export remained review-gated and
  required an explicit publication request and new output directory. The GUI
  did not claim the pending study was complete.

## Bound evidence

The preparation receipt SHA-256 is
`b796ddf03143c33ecb2eca6760c5aeb07bad05d7940dd1521821f3c4d7689bb4`;
the locus report is
`4f84705b05432cda7699d6fbe027f4c8fa99775b9cac904d241fd0583e0cac27`;
the study plan is
`7d4b3b2aefb18be97aed6f295dc0b2316fa687b11564a796525338f9b0401071`;
and the resulting panel is
`efe60b2fc2571aa1fe149a37acab96284981b6a9f690e223b0c31698c2d0cfa6`.

The public captures, snapshots and deterministic teaching crops are described
in [their evidence ledger](screenshots/gene_assay_study_gui/evidence.json). Raw
captures were not cropped or edited. Exploratory/failed captures were excluded
from the tutorial rather than edited into apparent success.

## Remaining boundaries

This 13-transcript fixture does not establish pagination or performance for
more than 40 transcripts, 24 assays or 80 band rows. No large-panel fixture was
silently substituted. The inner agent was not invoked. No BLAST search,
authentic genomic specificity, whole-transcriptome specificity, laboratory
validation, order form or automatic ordering occurred.

