# Affymetrix platform resources

This directory contains GENtle-owned resource specifications for Affymetrix /
Thermo Fisher expression array platforms and the local staging directories for
login-walled vendor support files. GENtle may inspect these descriptors and
check local file status, but it must not auto-download vendor payloads, raw CEL
files, CDF files, or R/Bioconductor packages.

The platform registry lives at:

- `platform_registry.json` (`gentle.affymetrix_platform_registry.v1`)

Registry entries are deliberately allowed to be `provisional`: GENtle can
recognize the chip family and report what annotation/CDF/library files are
needed without pretending that every historical platform has been locally
validated. Entries should move to `verified` only after the package names,
support-file patterns, and coordinate assumptions have been checked against real
annotation resources.

Local vendor support files belong under platform-specific subdirectories such as
`clariom_d_human_na36_hg38/`. Those payloads are not committed; nearby README
files document accepted canonical and browser-preserved download filenames.
