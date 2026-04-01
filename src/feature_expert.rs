//! Expert-view helpers for feature-centric deep-inspection UIs.
//!
//! Most portable expert-view payloads now live in `gentle-protocol`; this
//! module re-exports them for compatibility, and it also re-exports the
//! shared transition-matrix helper that now lives in `gentle-render`.

pub use gentle_protocol::{
    FeatureExpertTarget, FeatureExpertView, ISOFORM_ARCHITECTURE_EXPERT_INSTRUCTION,
    IsoformArchitectureCdsAaSegment, IsoformArchitectureExpertView,
    IsoformArchitectureProteinDomain, IsoformArchitectureProteinLane,
    IsoformArchitectureTranscriptLane, RESTRICTION_EXPERT_INSTRUCTION, RestrictionSiteExpertView,
    SPLICING_EXPERT_INSTRUCTION, SplicingBoundaryMarker, SplicingEventSummary,
    SplicingExonCdsPhase, SplicingExonSummary, SplicingExpertView, SplicingJunctionArc,
    SplicingMatrixRow, SplicingRange, SplicingScopePreset, SplicingTranscriptLane,
    TFBS_EXPERT_INSTRUCTION, TfbsExpertColumn, TfbsExpertView,
};
pub use gentle_render::{SplicingExonTransitionMatrix, compute_splicing_exon_transition_matrix};
