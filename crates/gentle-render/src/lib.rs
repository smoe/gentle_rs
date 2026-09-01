//! Headless GENtle export and rendering helpers.
//!
//! This crate owns the reusable non-egui render/export paths that should stay
//! adapter-neutral: feature-expert SVG, lineage SVG, protocol cartoons, and
//! the shared gel-figure renderers used by workflow-driven demos.

mod feature_expert;
mod gene_set_publication;
pub mod pool_gel;
pub mod protein_gel;
pub mod protocol_cartoon;

pub use feature_expert::{
    GeneLocusEvidenceOverlay, GeneLocusEvidenceOverlayRow, GeneLocusEvidenceOverlaySchematicTail,
    GeneLocusEvidenceOverlaySegment, SplicingExonTransition, SplicingExonTransitionMatrix,
    compute_splicing_exon_transition_matrix, compute_supported_splicing_exon_transitions,
    render_cryptic_splicing_screen, render_feature_expert_svg,
    render_gene_locus_evidence_with_overlay_svg,
};
pub use gene_set_publication::{
    render_gene_isoform_assay_publication_gene, render_gene_isoform_assay_publication_index,
    render_gene_isoform_assay_publication_print, render_gene_set_publication_html,
    render_gene_set_publication_markdown,
};
