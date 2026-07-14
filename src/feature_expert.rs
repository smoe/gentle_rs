//! Expert-view helpers for feature-centric deep-inspection UIs.
//!
//! Most portable expert-view payloads now live in `gentle-protocol`; this
//! module re-exports them for compatibility, and it also re-exports the
//! shared transition-matrix helper that now lives in `gentle-render`.

pub use gentle_protocol::{
    CDNA_EST_EVIDENCE_RESOURCE_SCHEMA, CdnaEstEvidenceKind, CdnaEstEvidenceRecord,
    CdnaEstEvidenceResource, FeatureExpertTarget, FeatureExpertView,
    GENE_ISOFORM_EVIDENCE_INSTRUCTION, GENE_ISOFORM_EVIDENCE_SCHEMA,
    GENE_LOCUS_EVIDENCE_DISPLAY_INSTRUCTION, GENE_LOCUS_EVIDENCE_DISPLAY_SCHEMA,
    GENE_LOCUS_OCCUPANCY_LAYOUT_SCHEMA, GeneIsoformAssayCandidate, GeneIsoformEvidenceComponent,
    GeneIsoformEvidenceComponents, GeneIsoformEvidenceItem, GeneIsoformEvidenceProvenanceSource,
    GeneIsoformEvidenceReport, GeneIsoformEvidenceRequest, GeneIsoformExonFamilyRow,
    GeneIsoformFamilyRow, GeneIsoformJunctionRow, GeneIsoformOccupancyInterval,
    GeneIsoformOccupancyLane, GeneIsoformTranscriptRow, GeneLocusAssayOverlay, GeneLocusCodonKind,
    GeneLocusCodonMarker, GeneLocusEvidenceDisplayReport, GeneLocusEvidenceDisplayRequest,
    GeneLocusMotifHit, GeneLocusMotifTrack, GeneLocusOccupancyGroup,
    GeneLocusOccupancyGroupRequest, GeneLocusOccupancyLane, GeneLocusOccupancyLaneRequest,
    GeneLocusOccupancyLaneRole, GeneLocusOccupancyLayout, GeneLocusOccupancyScaleMode,
    GeneLocusProbeClass, GeneLocusProbeEffectContrast, GeneLocusProbeEffectOverlay,
    GeneLocusProbeEffectValue, GeneLocusTranscriptMetrics,
    ISOFORM_ARCHITECTURE_EXPERT_INSTRUCTION,
    IsoformArchitectureCdsAaSegment, IsoformArchitectureExpertView,
    IsoformArchitectureProteinDomain, IsoformArchitectureProteinLane,
    IsoformArchitectureTranscriptLane, IsoformEvidenceAssessmentStatus, IsoformEvidenceSourceKind,
    IsoformExpressionMatrix, IsoformExpressionRow, RESTRICTION_EXPERT_INSTRUCTION,
    RestrictionSiteExpertView, SPLICING_EXPERT_INSTRUCTION,
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
