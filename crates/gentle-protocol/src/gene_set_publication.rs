//! Portable contracts for publication-oriented multi-gene analysis reports.
//!
//! The request keeps local paths and column mappings explicit. The resolved
//! report replaces those inputs with normalized primer rows and copied asset
//! links so HTML and print renderers consume one semantic record.

use serde::{Deserialize, Serialize};

pub const GENE_SET_PUBLICATION_REQUEST_SCHEMA: &str = "gentle.gene_set_publication_request.v1";
pub const GENE_SET_PUBLICATION_REPORT_SCHEMA: &str = "gentle.gene_set_publication_report.v1";
pub const GENE_SET_PUBLICATION_GENERATION_SCHEMA: &str =
    "gentle.gene_set_publication_generation.v1";
pub const GENE_SET_PUBLICATION_BUNDLE_MANIFEST_SCHEMA: &str =
    "gentle.gene_set_publication_bundle_manifest.v1";
pub const GENE_ISOFORM_ASSAY_PUBLICATION_REQUEST_SCHEMA: &str =
    "gentle.gene_isoform_assay_publication_request.v1";
pub const GENE_ISOFORM_ASSAY_PUBLICATION_SCHEMA: &str = "gentle.gene_isoform_assay_publication.v1";
pub const GENE_ISOFORM_ASSAY_PUBLICATION_PROJECTION_SCHEMA: &str =
    "gentle.gene_isoform_assay_publication_projection.v1";

/// One explicit dashboard value shown near the report title or one gene.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct GeneSetPublicationMetric {
    pub value: String,
    pub label: String,
    pub detail: String,
}

/// One narrative section shared by the web and printable reports.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct GeneSetPublicationNarrativeSection {
    pub section_id: String,
    pub heading: String,
    pub paragraphs: Vec<String>,
    pub bullets: Vec<String>,
}

/// Mapping from a source table into GENtle's normalized primer presentation.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(default)]
pub struct GeneSetPublicationPrimerColumnMap {
    pub gene: String,
    pub pair_id: String,
    pub role: String,
    pub status: String,
    pub forward_sequence_5_to_3: String,
    pub reverse_sequence_5_to_3: String,
    pub forward_tm_c: String,
    pub reverse_tm_c: String,
    pub cdna_specificity: String,
    pub genome_assessment: String,
    pub note: String,
}

impl Default for GeneSetPublicationPrimerColumnMap {
    fn default() -> Self {
        Self {
            gene: "gene".to_string(),
            pair_id: "pair_alias".to_string(),
            role: "selection_role".to_string(),
            status: "selection_status".to_string(),
            forward_sequence_5_to_3: "forward_sequence_5_to_3".to_string(),
            reverse_sequence_5_to_3: "reverse_sequence_5_to_3".to_string(),
            forward_tm_c: "forward_tm_c_gentle_shared".to_string(),
            reverse_tm_c: "reverse_tm_c_gentle_shared".to_string(),
            cdna_specificity: "whole_cdna_specificity".to_string(),
            genome_assessment: "genome_assessment".to_string(),
            note: "genome_detail".to_string(),
        }
    }
}

/// Delimited primer-table input resolved relative to the manifest.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(default)]
pub struct GeneSetPublicationPrimerTableSource {
    pub path: String,
    pub delimiter: String,
    pub columns: GeneSetPublicationPrimerColumnMap,
}

impl Default for GeneSetPublicationPrimerTableSource {
    fn default() -> Self {
        Self {
            path: String::new(),
            delimiter: "\t".to_string(),
            columns: GeneSetPublicationPrimerColumnMap::default(),
        }
    }
}

/// One source figure attached to a gene section.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct GeneSetPublicationFigureRequest {
    pub figure_id: String,
    pub kind: String,
    pub label: String,
    pub source_path: String,
    pub pdf_source_path: Option<String>,
    pub caption: String,
    pub alt_text: String,
    pub include_in_pdf: bool,
}

/// One gene-specific report section.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct GeneSetPublicationGeneRequest {
    pub gene_symbol: String,
    pub heading: String,
    pub summary: String,
    pub metrics: Vec<GeneSetPublicationMetric>,
    pub figures: Vec<GeneSetPublicationFigureRequest>,
}

/// One downloadable supporting artifact copied into the report bundle.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct GeneSetPublicationDownloadRequest {
    pub label: String,
    pub source_path: String,
    pub description: String,
}

/// Portable request for a complete HTML/print gene-set report bundle.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(default)]
pub struct GeneSetPublicationRequest {
    pub schema: String,
    pub report_id: String,
    pub title: String,
    pub subtitle: String,
    pub eyebrow: String,
    pub generated_date: String,
    pub summary: String,
    pub conditions: Vec<String>,
    pub metrics: Vec<GeneSetPublicationMetric>,
    pub sections: Vec<GeneSetPublicationNarrativeSection>,
    pub primer_table: GeneSetPublicationPrimerTableSource,
    pub genes: Vec<GeneSetPublicationGeneRequest>,
    pub downloads: Vec<GeneSetPublicationDownloadRequest>,
    pub provenance: Vec<String>,
    pub footer: String,
    pub favicon_source_path: Option<String>,
    pub html_filename: String,
    pub markdown_filename: String,
    pub pdf_filename: String,
}

impl Default for GeneSetPublicationRequest {
    fn default() -> Self {
        Self {
            schema: GENE_SET_PUBLICATION_REQUEST_SCHEMA.to_string(),
            report_id: String::new(),
            title: String::new(),
            subtitle: String::new(),
            eyebrow: "GENtle gene-set analysis".to_string(),
            generated_date: String::new(),
            summary: String::new(),
            conditions: vec![],
            metrics: vec![],
            sections: vec![],
            primer_table: GeneSetPublicationPrimerTableSource::default(),
            genes: vec![],
            downloads: vec![],
            provenance: vec![],
            footer: String::new(),
            favicon_source_path: None,
            html_filename: "index.html".to_string(),
            markdown_filename: "report.md".to_string(),
            pdf_filename: "report.pdf".to_string(),
        }
    }
}

/// One normalized primer row consumed by both output renderers.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct GeneSetPublicationPrimerRow {
    pub gene: String,
    pub pair_id: String,
    pub role: String,
    pub status: String,
    pub forward_sequence_5_to_3: String,
    pub reverse_sequence_5_to_3: String,
    pub forward_tm_c: String,
    pub reverse_tm_c: String,
    pub cdna_specificity: String,
    pub genome_assessment: String,
    pub note: String,
}

/// Copied figure paths used by HTML and PDF assembly.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct GeneSetPublicationFigure {
    pub figure_id: String,
    pub kind: String,
    pub label: String,
    pub web_path: String,
    pub pdf_path: Option<String>,
    pub caption: String,
    pub alt_text: String,
    pub include_in_pdf: bool,
}

/// One copied download link in the resolved report.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct GeneSetPublicationDownload {
    pub label: String,
    pub web_path: String,
    pub description: String,
}

/// One gene section after paths and primer rows have been normalized.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct GeneSetPublicationGene {
    pub gene_symbol: String,
    pub heading: String,
    pub summary: String,
    pub metrics: Vec<GeneSetPublicationMetric>,
    pub figures: Vec<GeneSetPublicationFigure>,
}

/// Shared semantic record rendered to responsive HTML and printable Markdown.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct GeneSetPublicationReport {
    pub schema: String,
    pub report_id: String,
    pub title: String,
    pub subtitle: String,
    pub eyebrow: String,
    pub generated_date: String,
    pub summary: String,
    pub conditions: Vec<String>,
    pub metrics: Vec<GeneSetPublicationMetric>,
    pub sections: Vec<GeneSetPublicationNarrativeSection>,
    pub primers: Vec<GeneSetPublicationPrimerRow>,
    pub genes: Vec<GeneSetPublicationGene>,
    pub downloads: Vec<GeneSetPublicationDownload>,
    pub provenance: Vec<String>,
    pub footer: String,
    pub favicon_path: Option<String>,
    pub pdf_path: Option<String>,
}

/// Machine-readable output summary from the local bundle generator.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct GeneSetPublicationGenerationReport {
    pub schema: String,
    pub request_path: String,
    pub output_directory: String,
    pub html_path: String,
    pub markdown_path: String,
    pub resolved_report_path: String,
    pub pdf_path: Option<String>,
    pub copied_files: Vec<String>,
    pub bundle_manifest_path: String,
    pub bundle_manifest_sha256: String,
    #[serde(default)]
    pub artifacts: Vec<GeneSetPublicationBundleArtifact>,
    pub warnings: Vec<String>,
}

/// One finalized legacy publication-bundle file bound to its exact bytes.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct GeneSetPublicationBundleArtifact {
    pub path: String,
    pub sha256: String,
    pub media_type: String,
}

/// Portable content inventory for one generated legacy publication bundle.
///
/// The manifest binds the input request and every user-facing artifact. It
/// intentionally excludes itself and `generation-report.json` to avoid a
/// recursive digest contract; the generation report records the manifest hash.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct GeneSetPublicationBundleManifest {
    pub schema: String,
    pub report_id: String,
    pub request_sha256: String,
    #[serde(default)]
    pub artifacts: Vec<GeneSetPublicationBundleArtifact>,
}

/// One content-addressed GENtle report consumed by an isoform-assay dossier.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default, deny_unknown_fields)]
pub struct GeneIsoformAssayPublicationReportRef {
    pub path: String,
    pub expected_sha256: String,
}

/// One gene and its immutable GENtle planning, handoff, and procurement inputs.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default, deny_unknown_fields)]
pub struct GeneIsoformAssayPublicationGeneRequest {
    pub gene_symbol: String,
    pub study_plan: GeneIsoformAssayPublicationReportRef,
    #[serde(default)]
    pub handoffs: Vec<GeneIsoformAssayPublicationReportRef>,
    #[serde(default)]
    pub order_forms: Vec<GeneIsoformAssayPublicationReportRef>,
    #[serde(default)]
    pub figures: Vec<GeneSetPublicationFigureRequest>,
}

/// Request for a multi-page, content-addressed gene isoform-assay dossier.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(default, deny_unknown_fields)]
pub struct GeneIsoformAssayPublicationRequest {
    pub schema: String,
    pub report_id: String,
    pub title: String,
    pub subtitle: String,
    pub generated_date: String,
    #[serde(default)]
    pub genes: Vec<GeneIsoformAssayPublicationGeneRequest>,
    pub default_profile: String,
    pub footer: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub favicon_source_path: Option<String>,
    pub index_filename: String,
    pub print_filename: String,
    pub pdf_filename: String,
}

impl Default for GeneIsoformAssayPublicationRequest {
    fn default() -> Self {
        Self {
            schema: GENE_ISOFORM_ASSAY_PUBLICATION_REQUEST_SCHEMA.to_string(),
            report_id: String::new(),
            title: String::new(),
            subtitle: String::new(),
            generated_date: String::new(),
            genes: vec![],
            default_profile: "full".to_string(),
            footer: String::new(),
            favicon_source_path: None,
            index_filename: "index.html".to_string(),
            print_filename: "print.html".to_string(),
            pdf_filename: "report.pdf".to_string(),
        }
    }
}

/// One embedded source report. Its JSON value is authoritative; content blocks
/// point into it rather than copying scientific statements.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct GeneIsoformAssayPublicationBoundReport {
    pub source_path: String,
    pub sha256: String,
    pub schema: String,
    pub report_id: String,
    pub value: serde_json::Value,
}

/// One effective planner parameter and all exact source pointers that agree.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default)]
pub struct GeneIsoformAssayPublicationParameter {
    pub name: String,
    pub value: serde_json::Value,
    #[serde(default)]
    pub source_pointers: Vec<String>,
}

/// One gene-specific departure from the common planner parameters.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default)]
pub struct GeneIsoformAssayPublicationParameterOverride {
    pub gene_symbol: String,
    pub name: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub common_value: Option<serde_json::Value>,
    pub effective_value: serde_json::Value,
    pub reason: String,
    pub source_pointer: String,
}

/// One gene page in the canonical dossier data tree.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default)]
pub struct GeneIsoformAssayPublicationGene {
    pub gene_symbol: String,
    pub page_path: String,
    pub study_plan: GeneIsoformAssayPublicationBoundReport,
    #[serde(default)]
    pub handoffs: Vec<GeneIsoformAssayPublicationBoundReport>,
    #[serde(default)]
    pub order_forms: Vec<GeneIsoformAssayPublicationBoundReport>,
    #[serde(default)]
    pub figures: Vec<GeneSetPublicationFigure>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub order_sheet_path: Option<String>,
    #[serde(default)]
    pub warnings: Vec<String>,
}

/// Pointer-only description of one renderable dossier section.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct GeneIsoformAssayPublicationBlock {
    pub block_id: String,
    pub label: String,
    pub projection: String,
    pub source_pointer: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub gene_symbol: Option<String>,
}

/// Named, ordered set of blocks that a presentation client may request.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct GeneIsoformAssayPublicationProfile {
    pub profile_id: String,
    pub label: String,
    #[serde(default)]
    pub block_ids: Vec<String>,
}

/// Canonical dossier record. Renderers may select declared blocks but may not
/// introduce or rewrite scientific content.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default)]
pub struct GeneIsoformAssayPublicationReport {
    pub schema: String,
    pub report_id: String,
    pub title: String,
    pub subtitle: String,
    pub generated_date: String,
    pub default_profile: String,
    pub index_path: String,
    pub print_path: String,
    #[serde(default)]
    pub common_parameters: Vec<GeneIsoformAssayPublicationParameter>,
    #[serde(default)]
    pub parameter_overrides: Vec<GeneIsoformAssayPublicationParameterOverride>,
    #[serde(default)]
    pub genes: Vec<GeneIsoformAssayPublicationGene>,
    #[serde(default)]
    pub content_blocks: Vec<GeneIsoformAssayPublicationBlock>,
    #[serde(default)]
    pub profiles: Vec<GeneIsoformAssayPublicationProfile>,
    pub footer: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub favicon_path: Option<String>,
    #[serde(default)]
    pub warnings: Vec<String>,
}

/// One generated file bound to its exact bytes.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct GeneIsoformAssayPublicationArtifact {
    pub path: String,
    pub sha256: String,
    pub media_type: String,
}

/// Receipt for one deterministic presentation projection of a canonical report.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct GeneIsoformAssayPublicationProjectionReport {
    pub schema: String,
    pub canonical_report_path: String,
    pub canonical_report_sha256: String,
    pub selected_profile: String,
    #[serde(default)]
    pub selected_block_ids: Vec<String>,
    #[serde(default)]
    pub artifacts: Vec<GeneIsoformAssayPublicationArtifact>,
    #[serde(default)]
    pub warnings: Vec<String>,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn request_defaults_keep_output_names_and_column_contract() {
        let request: GeneSetPublicationRequest = serde_json::from_value(serde_json::json!({
            "title": "Synthetic report",
            "genes": [{"gene_symbol": "GENE1"}]
        }))
        .expect("deserialize request defaults");
        assert_eq!(request.schema, GENE_SET_PUBLICATION_REQUEST_SCHEMA);
        assert_eq!(request.html_filename, "index.html");
        assert_eq!(request.pdf_filename, "report.pdf");
        assert_eq!(request.primer_table.columns.pair_id, "pair_alias");
        assert_eq!(request.genes[0].gene_symbol, "GENE1");
    }

    #[test]
    fn legacy_generation_receipt_without_bundle_manifest_fields_still_deserializes() {
        let report: GeneSetPublicationGenerationReport =
            serde_json::from_value(serde_json::json!({
                "schema": GENE_SET_PUBLICATION_GENERATION_SCHEMA,
                "html_path": "index.html",
                "copied_files": ["data/primers.tsv"]
            }))
            .expect("deserialize pre-manifest generation receipt");
        assert!(report.bundle_manifest_path.is_empty());
        assert!(report.bundle_manifest_sha256.is_empty());
        assert!(report.artifacts.is_empty());
    }
}
