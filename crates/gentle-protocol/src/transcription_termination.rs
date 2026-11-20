//! Provenance-bound termination-element recognition and explicitly uncalibrated
//! bacterial/mammalian candidate heuristics. Source base intervals are zero-based,
//! half-open; related interbase-boundary sets declare their coordinate basis.
//! Mature RNA ends, processing signals and polymerase termination remain distinct.

use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

pub const TERMINATION_SCREEN_SCHEMA: &str = "gentle.transcription_termination_screen.v1";
pub const TERMINATION_CATALOG_SCHEMA: &str = "gentle.transcription_termination_catalog.v1";
pub const TERMINATION_NON_CLAIM: &str = "Reference identity and sequence heuristics do not establish expression, cleavage-site usage, termination efficiency, exact mature RNA ends or polymerase release. Authored heuristic thresholds have no measured sensitivity or specificity. No candidate is not biological absence; literal A/T runs are not polyadenylation-site proof. Bacterial profiles do not cover archaea; the mammalian Pol II profile does not cover all eukaryotes.";

/// Local-only catalog; the digest binds raw file bytes, not normalized JSON.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct TerminationCatalogSource {
    pub path: String,
    pub expected_sha256: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
/// Versioned catalog document. Entries retain sequence and rights/source assertions;
/// the engine validates entry digests before exact recognition.
pub struct TerminationCatalog {
    pub schema: String,
    pub catalog_id: String,
    pub version: String,
    pub entries: Vec<TerminationCatalogEntry>,
}

/// Full element sequence and provenance are distinct from a short signal motif.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct TerminationCatalogEntry {
    pub id: String,
    pub label: String,
    /// Canonical DNA in the element's biological 5′-to-3′ orientation.
    pub sequence: String,
    pub sequence_sha256: String,
    pub host: String,
    pub mechanism: String,
    pub element_kind: TerminationGeometryKind,
    pub source_url: String,
    pub source_version: String,
    pub recreation: String,
    /// A source declaration; hashing does not authenticate rights or authorship.
    pub redistribution_status: String,
    pub approval_provenance: String,
}

/// Authored structural heuristic policy; scores are not folding energies.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(default, deny_unknown_fields)]
pub struct BacterialIntrinsicPolicy {
    pub version: String,
    pub min_stem_bp: usize,
    pub max_stem_bp: usize,
    pub min_loop_bp: usize,
    pub max_loop_bp: usize,
    pub max_mismatches: usize,
    pub allow_gu_wobble: bool,
    pub min_gc_pairs: usize,
    pub max_tract_gap_bp: usize,
    pub tract_window_bp: usize,
    pub min_tract_t: usize,
}
impl Default for BacterialIntrinsicPolicy {
    fn default() -> Self {
        Self {
            version: "bacterial_intrinsic_heuristic.v1".into(),
            min_stem_bp: 6,
            max_stem_bp: 12,
            min_loop_bp: 3,
            max_loop_bp: 8,
            max_mismatches: 1,
            allow_gu_wobble: true,
            min_gc_pairs: 4,
            max_tract_gap_bp: 2,
            tract_window_bp: 8,
            min_tract_t: 5,
        }
    }
}

/// Narrow mammalian Pol II signal/context policy, not a site-usage predictor.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(default, deny_unknown_fields)]
pub struct MammalianPolyAPolicy {
    pub version: String,
    pub signals: Vec<String>,
    /// Cleavage-boundary candidates measured from the signal's end.
    pub min_cleavage_offset_bp: usize,
    pub max_cleavage_offset_bp: usize,
    pub downstream_window_bp: usize,
    pub min_downstream_gt: usize,
}
impl Default for MammalianPolyAPolicy {
    fn default() -> Self {
        Self {
            version: "mammalian_pol_ii_heuristic.v1".into(),
            signals: vec!["AATAAA".into(), "ATTAAA".into()],
            min_cleavage_offset_bp: 10,
            max_cleavage_offset_bp: 30,
            downstream_window_bp: 20,
            min_downstream_gt: 10,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(default, deny_unknown_fields)]
/// One source-DNA screen with explicit subject context and bounded deterministic work.
/// Predictions are opt-in and catalog reads are local-only, content-bound inputs.
pub struct TerminationScreenRequest {
    pub seq_id: String,
    pub start_0based: usize,
    pub end_0based_exclusive: Option<usize>,
    /// Explicit mRNA/CDS context. An explicit interval may include downstream flanks;
    /// absent bounds use its exact span, never an inferred mature transcript.
    pub subject_feature_id: Option<usize>,
    pub host: String,
    pub polymerase: String,
    pub catalog: Option<TerminationCatalogSource>,
    pub include_annotations: bool,
    pub bacterial_intrinsic: bool,
    pub mammalian_pol_ii: bool,
    pub requested_mechanisms: Vec<String>,
    pub bacterial_policy: BacterialIntrinsicPolicy,
    pub mammalian_policy: MammalianPolyAPolicy,
    pub max_rows: usize,
    pub max_work: usize,
}
impl Default for TerminationScreenRequest {
    fn default() -> Self {
        Self {
            seq_id: String::new(),
            start_0based: 0,
            end_0based_exclusive: None,
            subject_feature_id: None,
            host: String::new(),
            polymerase: String::new(),
            catalog: None,
            include_annotations: true,
            bacterial_intrinsic: false,
            mammalian_pol_ii: false,
            requested_mechanisms: vec![],
            bacterial_policy: BacterialIntrinsicPolicy::default(),
            mammalian_policy: MammalianPolyAPolicy::default(),
            max_rows: 5000,
            max_work: 20_000_000,
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, PartialOrd, Ord)]
#[serde(rename_all = "snake_case")]
/// Evidence origin, kept independent of geometry and experimental validation.
pub enum TerminationEvidenceKind {
    SourceAnnotation,
    ReferenceMatch,
    BacterialIntrinsicCandidate,
    MammalianPolIiCandidate,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, PartialOrd, Ord)]
#[serde(rename_all = "snake_case")]
/// Biological meaning of an interval; signal and cleavage hypotheses never imply
/// an observed mature RNA end or polymerase-release coordinate.
pub enum TerminationGeometryKind {
    ElementMatch,
    ProcessingSignal,
    CleavageWindow,
    AnnotatedRnaEnd,
    PolymeraseTerminationRegion,
    IntrinsicTerminatorCandidate,
    CandidateRegion,
}

/// Distinguish source-base spans from sets of candidate interbase cut labels.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum TerminationCoordinateBasis {
    /// Bases in `[start, end)`; reverse projection is `[length-end, length-start)`.
    BaseInterval,
    /// Integer boundaries between bases, labeled `0..=length`; `[start, end)`
    /// is a set of possible cut labels, not a span of nucleotides. Reverse maps
    /// each label `p` to `length-p`, so the set becomes
    /// `[length-(end-1), length-start+1)`.
    InterbaseBoundarySet,
}

/// A related source-DNA interval on the parent evidence row's strand. An
/// interbase set may have `end_0based_exclusive == sequence_length + 1` because
/// the boundary after the last base is a valid label.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct TerminationInterval {
    pub geometry_kind: TerminationGeometryKind,
    pub coordinate_basis: TerminationCoordinateBasis,
    pub start_0based: usize,
    pub end_0based_exclusive: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
/// One inspectable observation or candidate. Coordinates are ascending source-DNA
/// bounds; `reverse` supplies biological direction. Scores are component values,
/// not calibrated probabilities, free energies or termination efficiencies.
pub struct TerminationEvidenceRow {
    pub row_id: String,
    pub evidence_kind: TerminationEvidenceKind,
    pub geometry_kind: TerminationGeometryKind,
    pub label: String,
    pub start_0based: usize,
    pub end_0based_exclusive: usize,
    pub reverse: bool,
    pub element_id: Option<String>,
    pub source_feature_id: Option<usize>,
    pub subject_feature_id: Option<usize>,
    pub sequence_5prime_to_3prime: Option<String>,
    pub components: BTreeMap<String, f64>,
    pub related_intervals: Vec<TerminationInterval>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
/// Mechanism-specific assessment or refusal, distinct from candidate counts.
pub struct TerminationMechanismOutcome {
    pub mechanism: String,
    pub status: String,
    pub message: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
/// Portable result binding the effective request, input digests, topology and
/// catalog provenance. Empty rows alone do not establish a completed negative scan.
pub struct TerminationScreenReport {
    pub schema: String,
    pub report_id: String,
    /// Hash of this report serialized with this field empty.
    pub content_sha256: String,
    pub request: TerminationScreenRequest,
    pub sequence_sha256: String,
    pub annotation_sha256: String,
    pub catalog_sha256: Option<String>,
    /// Validated catalog retained for portable source/provenance inspection.
    pub catalog: Option<TerminationCatalog>,
    pub topology: String,
    pub status: String,
    pub mechanism_outcomes: Vec<TerminationMechanismOutcome>,
    pub rows: Vec<TerminationEvidenceRow>,
    pub scanned_windows: usize,
    pub ambiguous_windows: usize,
    pub truncated: bool,
    pub warnings: Vec<String>,
    pub non_claim: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
/// Selected rows from an exact reviewed report. Materialization recomputes inputs
/// and refuses stale, unavailable or truncated evidence before changing features.
pub struct TerminationMaterializeRequest {
    pub screen: TerminationScreenRequest,
    pub expected_report_sha256: String,
    pub row_ids: Vec<String>,
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn strict_requests_and_explicit_defaults_round_trip() {
        let request: TerminationScreenRequest = serde_json::from_str("{}").unwrap();
        assert!(!request.bacterial_intrinsic && !request.mammalian_pol_ii);
        assert_eq!(request.bacterial_policy.min_stem_bp, 6);
        assert_eq!(
            request,
            serde_json::from_str(&serde_json::to_string(&request).unwrap()).unwrap()
        );
        assert!(serde_json::from_str::<TerminationScreenRequest>(r#"{"unknown":true}"#).is_err());
        assert!(
            serde_json::from_str::<TerminationScreenRequest>(
                r#"{"catalog":{"path":"x","expected_sha256":"sha256:test","download":true}}"#
            )
            .is_err()
        );
        assert!(
            serde_json::from_str::<TerminationScreenRequest>(
                r#"{"bacterial_policy":{"unknown":1}}"#
            )
            .is_err()
        );
        assert!(
            serde_json::from_str::<TerminationMaterializeRequest>(
                r#"{"screen":{},"expected_report_sha256":"x","row_ids":[],"confirm":true}"#
            )
            .is_err()
        );
    }
}
