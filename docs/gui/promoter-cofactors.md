### Promoter Cofactors

<!-- manual-counterparts:start -->
Related CLI reference: [Variant and Reporter Evidence](../cli/variants-reporters.md) | [Promoter Cofactors](../cli/promoter-cofactors.md).
These topic links do not imply identical feature coverage across interfaces.
<!-- manual-counterparts:end -->

**File > Promoter Cofactors...** and the matching command-palette entry open a
background-loaded collaborator-package browser, independent of loaded genomes.
It compares TA/DN cohort statistics and inspects strongest-per-band motif hits,
sample support/depth and promoter/gene ownership. See the
[wet-lab walkthrough and limits](../promoter_cofactor_browser.md).
Selected anchors, retained motif hits and promoters can explicitly **Save evidence
region** through the existing region ledger, or **Copy region request**. Edited
controls flag stale results; copy/save remains bound to the displayed report.
The [registered synthetic tutorial](../tutorial/08-14_promoter_cofactor_browser.md)
explains reference-free inspection and exact-reference feature attachment.
After saving a retained hit, open **Saved motif: DNA annotation**, select loaded
DNA and **Preview DNA annotation**. Review the local interval/strand and approve
before **Attach DNA annotation**. The shared engine verifies the exact prepared
reference. **Reload saved motifs** restores saved choices after reopening a
project or editing the region elsewhere, and clears any previous approval.
Stale previews are rejected; annotation is undoable. Unstranded hits,
anchors and promoters remain saved regions, not misleading directional features.
Live desktop and scientific acceptance remain separate from synthetic tests.

<!-- manual-page-nav:start -->

---

[GUI Manual](../gui.md) | [Previous: Sequencing Confirmation](sequencing-confirmation.md) | [Next: Accession-Pinned TSS Profiles](accession-pinned-tss-profiles.md)
<!-- manual-page-nav:end -->
