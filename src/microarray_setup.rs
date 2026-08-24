//! Local microarray setup discovery for GUI guidance.
//!
//! This module performs bounded checks inside the checkout that built GENtle
//! and, when different, the current working checkout. It never downloads data,
//! scans arbitrary user directories, or treats synthetic fixtures as
//! biological evidence.

use std::path::{Path, PathBuf};

/// Canonical checkout-relative staging directory for Clariom D Human support files.
pub const CLARIOM_D_HUMAN_SUPPORT_DIR: &str = "data/resources/affymetrix/clariom_d_human_na36_hg38";
/// Checkout-relative directory for the public TP73 Clariom D example.
pub const TP73_PUBLICATION_DATASET_DIR: &str =
    "data/publication_resources/rostock_p73_clariomd_e_mtab_14704";
/// Checkout-relative prepared synthetic output used for deterministic demonstrations.
pub const TP73_SYNTHETIC_OUTPUT_DIR: &str =
    "test_files/fixtures/probe_region_outputs/clariom_e_mtab_14704_tp73_validation";
/// Documentation for manually staged, login-walled vendor support files.
pub const CLARIOM_D_HUMAN_SUPPORT_README: &str =
    "data/resources/affymetrix/clariom_d_human_na36_hg38/README.md";
/// Public source page for the E-MTAB-14704 CEL dataset.
pub const E_MTAB_14704_LANDING_URL: &str =
    "https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-14704";
/// Thermo Fisher probeset annotation archive documented by the platform registry.
pub const CLARIOM_D_PROBESET_URL: &str = "https://sec-assets.thermofisher.com/TFS-Assets/LSG/Support-Files/Clariom_D_Human-na36-hg38-probeset-csv.zip";
/// Thermo Fisher transcript annotation archive documented by the platform registry.
pub const CLARIOM_D_TRANSCRIPT_URL: &str = "https://sec-assets.thermofisher.com/TFS-Assets/LSG/Support-Files/Clariom_D_Human.r1.na36.hg38.a1.transcript.csv.zip";

const CLARIOM_D_PROBESET_NAMES: &[&str] = &[
    "Clariom_D_Human-na36-hg38-probeset-csv.zip",
    "TFS-Assets_LSG_Support-Files_Clariom_D_Human-na36-hg38-probeset-csv.zip",
];
const CLARIOM_D_TRANSCRIPT_NAMES: &[&str] = &[
    "Clariom_D_Human.r1.na36.hg38.a1.transcript.csv.zip",
    "TFS-Assets_LSG_Support-Files_Clariom_D_Human.r1.na36.hg38.a1.transcript.csv.zip",
];
const SYNTHETIC_OUTPUT_REQUIRED_FILES: &[&str] = &[
    "region_intensity_chrom_order.csv",
    "normalized_feature_matrix_manifest.json",
    "provenance.json",
];

/// One bounded local resource check and the path GENtle inspected.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MicroarrayFileSuggestion {
    /// Human-readable input role.
    pub role: &'static str,
    /// Absolute path rooted in the checkout used for discovery.
    pub path: PathBuf,
    /// Whether the path currently exists as a file.
    pub exists: bool,
}

/// Local Clariom D resources that can be suggested without searching outside a checkout.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MicroarraySetupDiscovery {
    /// Checkout root used for every bounded path check.
    pub checkout_root: PathBuf,
    /// NetAffx probeset archive, accepting canonical and browser-preserved names.
    pub vendor_probeset_zip: MicroarrayFileSuggestion,
    /// NetAffx transcript archive, accepting canonical and browser-preserved names.
    pub vendor_transcript_zip: MicroarrayFileSuggestion,
    /// Public experiment metadata (IDF).
    pub publication_idf: MicroarrayFileSuggestion,
    /// Public experiment sample metadata (SDRF).
    pub publication_sdrf: MicroarrayFileSuggestion,
    /// Locally available CEL files immediately under the publication dataset directory.
    pub publication_cel_paths: Vec<PathBuf>,
    /// Existing full-platform normalized feature matrix shown as a source artifact.
    pub local_summary_tsv: MicroarrayFileSuggestion,
    /// Paired full-platform annotation, which may include coordinate-less rows.
    pub local_annotation_tsv: MicroarrayFileSuggestion,
    /// Existing sample/condition table paired with the local summary matrix.
    pub local_metadata_tsv: MicroarrayFileSuggestion,
    /// Prepared synthetic output directory for safe workflow demonstrations.
    pub synthetic_output_dir: PathBuf,
    /// Number of required helper-output files present in the synthetic directory.
    pub synthetic_output_files_present: usize,
}

impl MicroarraySetupDiscovery {
    /// Discover resources in the richer of the build and current working checkouts.
    pub fn discover_checkout() -> Self {
        let manifest_root = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        let mut roots = Vec::new();
        if let Ok(current_root) = std::env::current_dir()
            && current_root != manifest_root
        {
            roots.push(current_root);
        }
        roots.push(manifest_root);

        roots
            .into_iter()
            .map(|root| Self::discover_in(&root))
            .max_by_key(Self::readiness_score)
            .expect("at least the build checkout is inspected")
    }

    /// Discover only the documented paths below `checkout_root`.
    pub fn discover_in(checkout_root: &Path) -> Self {
        let checkout_root = checkout_root.to_path_buf();
        let support_dir = checkout_root.join(CLARIOM_D_HUMAN_SUPPORT_DIR);
        let dataset_dir = checkout_root.join(TP73_PUBLICATION_DATASET_DIR);
        let analysis_dir = dataset_dir.join("analysis/clariomd_probe_level");
        let synthetic_output_dir = checkout_root.join(TP73_SYNTHETIC_OUTPUT_DIR);

        let vendor_probeset_zip = first_named_file(
            &support_dir,
            CLARIOM_D_PROBESET_NAMES,
            "Clariom D probeset annotation ZIP",
        );
        let vendor_transcript_zip = first_named_file(
            &support_dir,
            CLARIOM_D_TRANSCRIPT_NAMES,
            "Clariom D transcript annotation ZIP",
        );
        let mut publication_cel_paths = std::fs::read_dir(&dataset_dir)
            .ok()
            .into_iter()
            .flatten()
            .filter_map(Result::ok)
            .map(|entry| entry.path())
            .filter(|path| {
                path.is_file()
                    && path
                        .extension()
                        .is_some_and(|extension| extension.eq_ignore_ascii_case("cel"))
            })
            .collect::<Vec<_>>();
        publication_cel_paths.sort();

        let synthetic_output_files_present = SYNTHETIC_OUTPUT_REQUIRED_FILES
            .iter()
            .filter(|name| synthetic_output_dir.join(name).is_file())
            .count();

        Self {
            checkout_root,
            vendor_probeset_zip,
            vendor_transcript_zip,
            publication_idf: file_suggestion(
                "E-MTAB-14704 IDF metadata",
                dataset_dir.join("E-MTAB-14704.idf.txt"),
            ),
            publication_sdrf: file_suggestion(
                "E-MTAB-14704 SDRF sample metadata",
                dataset_dir.join("E-MTAB-14704.sdrf.txt"),
            ),
            publication_cel_paths,
            local_summary_tsv: file_suggestion(
                "normalized probeset expression matrix",
                analysis_dir.join("probeset_expression_rma.tsv"),
            ),
            local_annotation_tsv: file_suggestion(
                "probeset feature annotation",
                analysis_dir.join("probeset_feature_annotation.tsv"),
            ),
            local_metadata_tsv: file_suggestion(
                "sample and condition metadata",
                analysis_dir.join("sample_table.tsv"),
            ),
            synthetic_output_dir,
            synthetic_output_files_present,
        }
    }

    /// Whether all three full-platform analysis source tables are present.
    pub fn has_local_analysis_tables(&self) -> bool {
        self.local_summary_tsv.exists
            && self.local_annotation_tsv.exists
            && self.local_metadata_tsv.exists
    }

    /// Number of the two documented vendor annotation archives currently present.
    pub fn vendor_file_count(&self) -> usize {
        usize::from(self.vendor_probeset_zip.exists)
            + usize::from(self.vendor_transcript_zip.exists)
    }

    /// Whether the committed synthetic output satisfies its minimal helper contract.
    pub fn has_synthetic_output(&self) -> bool {
        self.synthetic_output_files_present == SYNTHETIC_OUTPUT_REQUIRED_FILES.len()
    }

    /// Number of files required for the minimal prepared-output demonstration.
    pub fn synthetic_output_required_file_count(&self) -> usize {
        SYNTHETIC_OUTPUT_REQUIRED_FILES.len()
    }

    fn readiness_score(&self) -> usize {
        self.vendor_file_count()
            + usize::from(self.publication_idf.exists)
            + usize::from(self.publication_sdrf.exists)
            + self.publication_cel_paths.len()
            + usize::from(self.local_summary_tsv.exists)
            + usize::from(self.local_annotation_tsv.exists)
            + usize::from(self.local_metadata_tsv.exists)
            + self.synthetic_output_files_present
    }
}

fn file_suggestion(role: &'static str, path: PathBuf) -> MicroarrayFileSuggestion {
    MicroarrayFileSuggestion {
        role,
        exists: path.is_file(),
        path,
    }
}

fn first_named_file(
    directory: &Path,
    accepted_names: &[&str],
    role: &'static str,
) -> MicroarrayFileSuggestion {
    let path = accepted_names
        .iter()
        .map(|name| directory.join(name))
        .find(|path| path.is_file())
        .unwrap_or_else(|| directory.join(accepted_names[0]));
    file_suggestion(role, path)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn discovery_accepts_vendor_aliases_and_finds_local_tp73_inputs() {
        let temp = tempfile::tempdir().expect("temporary checkout");
        let support = temp.path().join(CLARIOM_D_HUMAN_SUPPORT_DIR);
        let dataset = temp.path().join(TP73_PUBLICATION_DATASET_DIR);
        let analysis = dataset.join("analysis/clariomd_probe_level");
        let synthetic = temp.path().join(TP73_SYNTHETIC_OUTPUT_DIR);
        std::fs::create_dir_all(&support).expect("support directory");
        std::fs::create_dir_all(&analysis).expect("analysis directory");
        std::fs::create_dir_all(&synthetic).expect("synthetic output directory");

        for path in [
            support.join(CLARIOM_D_PROBESET_NAMES[1]),
            support.join(CLARIOM_D_TRANSCRIPT_NAMES[1]),
            dataset.join("E-MTAB-14704.idf.txt"),
            dataset.join("E-MTAB-14704.sdrf.txt"),
            dataset.join("sample.CEL"),
            analysis.join("probeset_expression_rma.tsv"),
            analysis.join("probeset_feature_annotation.tsv"),
            analysis.join("sample_table.tsv"),
        ] {
            std::fs::write(path, b"fixture\n").expect("write resource fixture");
        }
        for name in SYNTHETIC_OUTPUT_REQUIRED_FILES {
            std::fs::write(synthetic.join(name), b"fixture\n").expect("write output fixture");
        }

        let found = MicroarraySetupDiscovery::discover_in(temp.path());
        assert_eq!(found.vendor_file_count(), 2);
        assert!(
            found
                .vendor_probeset_zip
                .path
                .ends_with(CLARIOM_D_PROBESET_NAMES[1])
        );
        assert_eq!(found.publication_cel_paths.len(), 1);
        assert!(found.has_local_analysis_tables());
        assert!(found.has_synthetic_output());
    }

    #[test]
    fn discovery_reports_expected_paths_without_creating_or_downloading_files() {
        let temp = tempfile::tempdir().expect("temporary checkout");
        let found = MicroarraySetupDiscovery::discover_in(temp.path());

        assert_eq!(found.vendor_file_count(), 0);
        assert!(found.publication_cel_paths.is_empty());
        assert!(!found.has_local_analysis_tables());
        assert!(!found.has_synthetic_output());
        assert!(
            found
                .vendor_probeset_zip
                .path
                .ends_with(CLARIOM_D_PROBESET_NAMES[0])
        );
    }
}
