//! Expert-view helpers for feature-centric deep-inspection UIs.
//!
//! Most portable expert-view payloads now live in `gentle-protocol`; this
//! module re-exports them for compatibility, and it also re-exports the
//! shared transition-matrix helper that now lives in `gentle-render`.

pub use gentle_protocol::{
    CDNA_EST_EVIDENCE_RESOURCE_SCHEMA, CdnaEstEvidenceKind, CdnaEstEvidenceRecord,
    CdnaEstEvidenceResource, FeatureExpertTarget, FeatureExpertView,
    GENE_ISOFORM_EVIDENCE_INSTRUCTION, GENE_ISOFORM_EVIDENCE_SCHEMA,
    GeneIsoformAssayCandidate, GeneIsoformEvidenceComponent, GeneIsoformEvidenceComponents,
    GeneIsoformEvidenceItem, GeneIsoformEvidenceProvenanceSource, GeneIsoformEvidenceReport,
    GeneIsoformEvidenceRequest, GeneIsoformExonFamilyRow, GeneIsoformFamilyRow,
    GeneIsoformJunctionRow, GeneIsoformTranscriptRow, ISOFORM_ARCHITECTURE_EXPERT_INSTRUCTION,
    IsoformEvidenceAssessmentStatus, IsoformEvidenceSourceKind,
    IsoformArchitectureCdsAaSegment, IsoformArchitectureExpertView,
    IsoformArchitectureProteinDomain, IsoformArchitectureProteinLane,
    IsoformArchitectureTranscriptLane, IsoformExpressionMatrix, IsoformExpressionRow,
    RESTRICTION_EXPERT_INSTRUCTION, RestrictionSiteExpertView, SPLICING_EXPERT_INSTRUCTION,
    SplicingBoundaryMarker, SplicingEventSummary, SplicingExonCdsPhase, SplicingExonSummary,
    SplicingExpertView, SplicingIntronSignal, SplicingJunctionArc, SplicingMatrixRow,
    SplicingRange, SplicingScopePreset, SplicingTranscriptLane, TFBS_EXPERT_INSTRUCTION,
    TfbsExpertColumn, TfbsExpertView, TranscriptProteinComparison,
    TranscriptProteinComparisonStatus, TranscriptProteinDerivation,
    TranscriptProteinExternalOpinion,
};
pub use gentle_render::{
    SplicingExonTransition, SplicingExonTransitionMatrix,
    compute_splicing_exon_transition_matrix, compute_supported_splicing_exon_transitions,
};
