# RNA Secondary Structure

<!-- manual-counterparts:start -->
Related GUI reference: [RNA Structure (ssRNA)](../gui/rna-structure-ssrna.md).
These topic links do not imply identical feature coverage across interfaces.
<!-- manual-counterparts:end -->

RNA secondary-structure text command:

- `rna-info SEQ_ID`
  - Returns a JSON report from `RNAfold --noPS`.
  - Includes `stdout`/`stderr`, dot-bracket `structure`, optional `mfe_kcal_per_mol`, and command metadata.
  - Accepts only single-stranded RNA (`molecule_type` of `RNA`/`ssRNA`).
  - Uses external `RNAfold` executable (set `GENTLE_RNAFOLD_BIN` to override executable path).

<!-- manual-page-nav:start -->

---

[CLI Manual](../cli.md) | [Previous: Target Rescue Screen](rescue-screen.md) | [Next: MicroRNA Target Scans](microrna.md)
<!-- manual-page-nav:end -->
