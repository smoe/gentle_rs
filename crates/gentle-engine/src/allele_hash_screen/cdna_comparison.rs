//! SNP observations and conservative local RNA haplotypes within the existing allele screen.
//!
//! The caller authenticates raw sources and projects supported alignment blocks. This pure
//! routine never parses a second alignment format, invents bases, joins reads by name, forces
//! diploidy, assigns parental alleles or calls DNA copy number/trans-splicing.

use super::complement_base;
use ring::digest::{SHA256, digest};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};

pub const CDNA_COMPARISON_REQUEST_SCHEMA: &str =
    "gentle.cdna_allele_comparison_evidence_request.v1";
pub const CDNA_COMPARISON_SCHEMA: &str = "gentle.cdna_allele_comparison_evidence.v1";

/// A normalized uppercase locus sequence; sha256 hashes these sequence bytes, not FASTA bytes.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CdnaReference {
    pub assembly: String,
    pub contig: String,
    pub start_0based: usize,
    pub sequence: String,
    pub sha256: String,
}

/// Explicit transcript exons in ascending genomic coordinates; strand is never guessed.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CdnaTranscript {
    pub transcript_id: String,
    pub strand: char,
    pub exons: Vec<(usize, usize)>,
}

/// One actual supported base, in both original query and reference-forward orientations.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CdnaBaseObservation {
    pub reference_position_0based: usize,
    pub query_position_0based: usize,
    pub observed_base: char,
    pub original_base: char,
    pub query_to_reference_strand: char,
    pub phred: Option<u8>,
    pub mapq: Option<u8>,
}

/// A retained mapping hypothesis. Ineligible observations remain audit-only.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CdnaHypothesis {
    pub hypothesis_id: String,
    pub eligible: bool,
    pub bases: Vec<CdnaBaseObservation>,
    /// Actual supported skips only, never inferred transitions between segments.
    pub junctions: Vec<(usize, usize)>,
}

/// One producer-bound original unit, not a QNAME-derived molecule.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CdnaUnit {
    pub unit_id: String,
    pub all_query_views_assessed: bool,
    pub search_declared_complete: bool,
    pub hypotheses: Vec<CdnaHypothesis>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CdnaSample {
    pub sample_id: String,
    pub label: String,
    pub source_sha256: String,
    pub analysis_context_sha256: String,
    pub units: Vec<CdnaUnit>,
}

/// Explicit, uncalibrated pilot policy; these are not biological genotype thresholds.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct CdnaComparisonPolicy {
    pub minimum_base_quality: u8,
    pub minimum_mapping_quality: u8,
    pub minimum_alternate_units: usize,
    pub minimum_alternate_fraction: f64,
    pub minimum_haplotype_units: usize,
    pub minimum_assignment_sites: usize,
    pub maximum_assignment_mismatches: usize,
    pub maximum_assignment_phred_penalty: u64,
    pub maximum_haplotypes: usize,
}

impl Default for CdnaComparisonPolicy {
    fn default() -> Self {
        Self {
            minimum_base_quality: 20,
            minimum_mapping_quality: 20,
            minimum_alternate_units: 3,
            minimum_alternate_fraction: 0.1,
            minimum_haplotype_units: 3,
            minimum_assignment_sites: 2,
            maximum_assignment_mismatches: 1,
            maximum_assignment_phred_penalty: 40,
            maximum_haplotypes: 8,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct CdnaComparisonBudgets {
    pub max_units: usize,
    pub max_sites: usize,
    pub max_base_observations: usize,
    pub max_hypotheses_per_unit: usize,
    pub max_phase_sites_per_block: usize,
    pub max_comparisons: u64,
    pub max_switch_explanations: usize,
    pub max_output_bytes: usize,
}

impl Default for CdnaComparisonBudgets {
    fn default() -> Self {
        Self {
            max_units: 10_000,
            max_sites: 512,
            max_base_observations: 5_000_000,
            max_hypotheses_per_unit: 128,
            max_phase_sites_per_block: 64,
            max_comparisons: 100_000_000,
            max_switch_explanations: 10_000,
            max_output_bytes: 64_000_000,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CdnaComparisonRequest {
    pub schema: String,
    pub gene: String,
    pub reference: CdnaReference,
    pub transcripts: Vec<CdnaTranscript>,
    pub samples: Vec<CdnaSample>,
    #[serde(default)]
    pub supplied_sites: Vec<usize>,
    /// Optional explicit genomic marker basis for phasing only; every marker must be supplied.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub phase_marker_positions: Option<Vec<usize>>,
    #[serde(default)]
    pub policy: CdnaComparisonPolicy,
    #[serde(default)]
    pub budgets: CdnaComparisonBudgets,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CdnaReferenceSummary {
    pub assembly: String,
    pub contig: String,
    pub start_0based: usize,
    pub end_0based_exclusive: usize,
    pub sequence_sha256: String,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum CdnaCallStatus {
    ReferenceSupported,
    AlternateSupported,
    Ambiguous,
    Unassessed,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CdnaUnitAlleleCall {
    pub unit_id: String,
    pub position_0based: usize,
    pub status: CdnaCallStatus,
    pub base: Option<char>,
    pub minimum_phred: Option<u8>,
    pub possible_bases: Vec<char>,
    pub reasons: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CdnaSiteObservation {
    pub unit_id: String,
    pub hypothesis_id: String,
    pub eligible: bool,
    pub quality_qualified: bool,
    pub observation: CdnaBaseObservation,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CdnaSiteSampleCounts {
    pub sample_id: String,
    pub total_units: usize,
    pub callable_units: usize,
    pub base_counts: BTreeMap<char, usize>,
    pub ambiguous_units: usize,
    pub unassessed_units: usize,
    pub candidate_alternate_bases: Vec<char>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CdnaTranscriptSiteAnnotation {
    pub transcript_id: String,
    pub status: String,
    pub cdna_position_1based: Option<usize>,
    pub transcript_reference_base: Option<char>,
    pub transcript_alternate_bases: Vec<char>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CdnaSiteReport {
    pub position_0based: usize,
    pub reference_base: char,
    pub supplied: bool,
    pub candidate_alternate_bases: Vec<char>,
    pub samples: Vec<CdnaSiteSampleCounts>,
    pub transcripts: Vec<CdnaTranscriptSiteAnnotation>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CdnaHaplotypeAllele {
    pub position_0based: usize,
    pub reference_base: char,
    pub observed_base: char,
    pub transcript_annotations: Vec<CdnaTranscriptSiteAnnotation>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CdnaPhaseLink {
    pub sites_0based: Vec<usize>,
    pub alleles: Vec<char>,
    /// Distinct original units witnessing this local pattern, not an assembled molecule.
    pub unit_ids: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CdnaHaplotype {
    pub haplotype_id: String,
    /// Reference/position/base fingerprint permits comparisons without equating local H1 labels.
    pub fingerprint: String,
    pub alleles: Vec<CdnaHaplotypeAllele>,
    pub exact_seed_units: Vec<String>,
    pub phase_method: String,
    pub support_links: Vec<CdnaPhaseLink>,
    pub assigned_units: usize,
    pub error_tolerant_assigned_units: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CdnaHaplotypeAssignment {
    pub unit_id: String,
    pub status: String,
    pub callable_sites: usize,
    pub compatible_haplotype_ids: Vec<String>,
    pub mismatch_count: Option<usize>,
    pub phred_penalty: Option<u64>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CdnaPhaseBlock {
    pub block_id: String,
    pub sites_0based: Vec<usize>,
    pub status: String,
    pub candidate_haplotype_count: Option<usize>,
    pub excess_expressed_haplotype_diversity: bool,
    pub copy_number_status: String,
    pub haplotypes: Vec<CdnaHaplotype>,
    pub assignments: Vec<CdnaHaplotypeAssignment>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CdnaSwitchExplanation {
    pub left_haplotype_id: String,
    pub right_haplotype_id: String,
    /// Includes the flanking markers; this is not an exact inferred molecular breakpoint.
    pub marker_bracket_0based: (usize, usize),
    pub left_informative_markers: usize,
    pub right_informative_markers: usize,
    pub mismatch_count: usize,
    pub phred_penalty: u64,
    pub same_hypothesis_junctions: BTreeMap<String, Vec<(usize, usize)>>,
    /// Genomic-left and genomic-right marker positions on each original query axis.
    pub query_marker_positions: BTreeMap<String, Vec<(usize, usize)>>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CdnaSwitchCandidate {
    pub unit_id: String,
    pub block_id: String,
    pub status: String,
    pub coordinate_order: String,
    pub whole_template_mismatch_count: Option<usize>,
    pub competing_whole_template: bool,
    pub explanations: Vec<CdnaSwitchExplanation>,
    pub mechanism_status: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CdnaHaplotypeJunctionSupport {
    pub block_id: String,
    pub haplotype_id: String,
    /// A set of observed skips, not a new physical splice-chain identity.
    pub junctions: Vec<(usize, usize)>,
    pub possible_units: usize,
    pub exclusive_units: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CdnaSampleReport {
    pub sample_id: String,
    pub label: String,
    pub source_sha256: String,
    pub analysis_context_sha256: String,
    pub total_units: usize,
    pub units_with_incomplete_views: usize,
    pub units_without_complete_search_declaration: usize,
    pub site_calls: Vec<CdnaUnitAlleleCall>,
    pub site_observations: Vec<CdnaSiteObservation>,
    pub phase_blocks: Vec<CdnaPhaseBlock>,
    pub switch_candidates: Vec<CdnaSwitchCandidate>,
    pub haplotype_junction_support: Vec<CdnaHaplotypeJunctionSupport>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CdnaSampleDifference {
    pub sample_a: String,
    pub sample_b: String,
    pub position_0based: usize,
    pub reference_base: char,
    pub sample_a_candidate_bases: Vec<char>,
    pub sample_b_candidate_bases: Vec<char>,
    pub status: String,
}

/// RNA sequence evidence only. No forced second allele, genomic genotype or copy-number estimate.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CdnaComparisonReport {
    pub schema: String,
    pub content_sha256: String,
    pub request_sha256: String,
    pub gene: String,
    pub reference: CdnaReferenceSummary,
    pub transcripts: Vec<CdnaTranscript>,
    pub policy: CdnaComparisonPolicy,
    pub policy_status: String,
    pub discovery_assessed_positions: usize,
    /// Normalized explicit phase-marker basis; omitted when every reported site is eligible.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub phase_marker_positions: Option<Vec<usize>>,
    pub sites: Vec<CdnaSiteReport>,
    pub samples: Vec<CdnaSampleReport>,
    pub sample_differences: Vec<CdnaSampleDifference>,
    pub comparisons: u64,
    pub warnings: Vec<String>,
}

fn sha(bytes: &[u8]) -> String {
    let value = digest(&SHA256, bytes);
    format!(
        "sha256:{}",
        value
            .as_ref()
            .iter()
            .map(|b| format!("{b:02x}"))
            .collect::<String>()
    )
}

fn hash<T: Serialize>(value: &T) -> Result<String, String> {
    serde_json::to_vec(value)
        .map(|v| sha(&v))
        .map_err(|e| e.to_string())
}

fn is_hash(s: &str) -> bool {
    s.strip_prefix("sha256:").is_some_and(|s| {
        s.len() == 64
            && s.bytes()
                .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
    })
}

fn identity(s: &str) -> bool {
    !s.is_empty() && s.len() <= 512 && !s.chars().any(char::is_control)
}

fn canonical(c: char) -> bool {
    matches!(c, 'A' | 'C' | 'G' | 'T')
}

struct Meter<'a> {
    work: u64,
    limit: u64,
    progress: &'a mut dyn FnMut(usize, usize) -> bool,
}

impl Meter<'_> {
    fn charge(&mut self, amount: usize) -> Result<(), String> {
        self.work = self
            .work
            .checked_add(amount as u64)
            .filter(|v| *v <= self.limit)
            .ok_or_else(|| "cdna_comparison_work_budget".to_string())?;
        Ok(())
    }
    fn poll(&mut self, completed: usize, total: usize) -> Result<(), String> {
        if (self.progress)(completed, total) {
            Ok(())
        } else {
            Err("cdna_comparison_cancelled".into())
        }
    }
}

#[derive(Debug, Clone)]
struct Call {
    base: Option<char>,
    quality: Option<u8>,
    possible: BTreeSet<char>,
    reasons: BTreeSet<String>,
    ambiguous: bool,
}

impl Call {
    fn uncovered() -> Self {
        Self {
            base: None,
            quality: None,
            possible: BTreeSet::new(),
            reasons: BTreeSet::from(["site_not_covered".into()]),
            ambiguous: false,
        }
    }
}

fn complemented(base: char) -> char {
    if canonical(base) {
        complement_base(base as u8).expect("canonical base") as char
    } else {
        crate::iupac_code::IupacCode::letter_complement(base as u8) as char
    }
}

fn validate(request: &CdnaComparisonRequest) -> Result<(), String> {
    let p = &request.policy;
    let b = &request.budgets;
    if request.schema != CDNA_COMPARISON_REQUEST_SCHEMA
        || !identity(&request.gene)
        || !identity(&request.reference.assembly)
        || !identity(&request.reference.contig)
        || request.reference.sequence.is_empty()
        || request.reference.sequence.len() > 16_000_000
        || !request
            .reference
            .sequence
            .bytes()
            .all(|c| b"ACGTRYSWKMBDHVN".contains(&c))
        || request.reference.sha256 != sha(request.reference.sequence.as_bytes())
        || !(1..=256).contains(&request.transcripts.len())
        || !(1..=16).contains(&request.samples.len())
        || p.minimum_base_quality > 93
        || p.minimum_mapping_quality > 254
        || p.minimum_alternate_units == 0
        || p.minimum_haplotype_units == 0
        || !p.minimum_alternate_fraction.is_finite()
        || !(0.0..=1.0).contains(&p.minimum_alternate_fraction)
        || !(1..=64).contains(&p.minimum_assignment_sites)
        || p.maximum_assignment_mismatches > 64
        || p.maximum_assignment_phred_penalty > 10_000
        || !(1..=16).contains(&p.maximum_haplotypes)
        || !(1..=100_000).contains(&b.max_units)
        || !(1..=1024).contains(&b.max_sites)
        || !(1..=16_000_000).contains(&b.max_base_observations)
        || !(1..=128).contains(&b.max_hypotheses_per_unit)
        || !(1..=64).contains(&b.max_phase_sites_per_block)
        || !(1..=1_000_000_000).contains(&b.max_comparisons)
        || !(1..=100_000).contains(&b.max_switch_explanations)
        || !(1..=128_000_000).contains(&b.max_output_bytes)
        || request.supplied_sites.len() > b.max_sites
    {
        return Err("cdna_comparison_schema_reference_policy_or_budget".into());
    }
    let end = request
        .reference
        .start_0based
        .checked_add(request.reference.sequence.len())
        .ok_or_else(|| "cdna_comparison_reference_overflow".to_string())?;
    let in_reference = |p| request.reference.start_0based <= p && p < end;
    if request.supplied_sites.iter().any(|&p| !in_reference(p))
        || request.supplied_sites.iter().collect::<BTreeSet<_>>().len()
            != request.supplied_sites.len()
    {
        return Err("cdna_comparison_supplied_sites".into());
    }
    if request
        .phase_marker_positions
        .as_ref()
        .is_some_and(|markers| {
            markers.is_empty()
                || markers.len() > b.max_phase_sites_per_block
                || markers.iter().collect::<BTreeSet<_>>().len() != markers.len()
                || markers
                    .iter()
                    .any(|&p| !in_reference(p) || !request.supplied_sites.contains(&p))
        })
    {
        return Err("cdna_comparison_phase_marker_positions".into());
    }
    let mut transcripts = BTreeSet::new();
    let mut exon_count = 0usize;
    for t in &request.transcripts {
        exon_count = exon_count
            .checked_add(t.exons.len())
            .ok_or("cdna_comparison_annotation_budget")?;
        if !identity(&t.transcript_id)
            || !transcripts.insert(&t.transcript_id)
            || !matches!(t.strand, '+' | '-')
            || t.exons.is_empty()
            || exon_count > 100_000
            || t.exons
                .iter()
                .any(|&(a, z)| a >= z || a < request.reference.start_0based || z > end)
            || !t.exons.windows(2).all(|w| w[0].1 <= w[1].0)
        {
            return Err("cdna_comparison_transcript_geometry".into());
        }
    }
    let mut samples = BTreeSet::new();
    let mut physical_units = BTreeSet::new();
    let mut units = 0usize;
    let mut bases = 0usize;
    for sample in &request.samples {
        units = units
            .checked_add(sample.units.len())
            .ok_or("cdna_comparison_unit_budget")?;
        if !identity(&sample.sample_id)
            || !identity(&sample.label)
            || !samples.insert(&sample.sample_id)
            || !is_hash(&sample.source_sha256)
            || !is_hash(&sample.analysis_context_sha256)
            || units > b.max_units
        {
            return Err("cdna_comparison_sample_or_unit_budget".into());
        }
        let mut unit_ids = BTreeSet::new();
        for unit in &sample.units {
            if !identity(&unit.unit_id)
                || !unit_ids.insert(&unit.unit_id)
                || !physical_units.insert((&sample.source_sha256, &unit.unit_id))
                || unit.hypotheses.len() > b.max_hypotheses_per_unit
            {
                return Err("cdna_comparison_duplicate_unit_or_hypothesis_budget".into());
            }
            let mut hypotheses = BTreeSet::new();
            for h in &unit.hypotheses {
                if !identity(&h.hypothesis_id)
                    || !hypotheses.insert(&h.hypothesis_id)
                    || h.junctions.len() > 4096
                    || h.junctions
                        .iter()
                        .any(|&(a, z)| a >= z || a < request.reference.start_0based || z > end)
                    || !h.junctions.windows(2).all(|w| w[0].1 <= w[1].0)
                {
                    return Err("cdna_comparison_hypothesis_or_junction".into());
                }
                bases = bases
                    .checked_add(h.bases.len())
                    .ok_or("cdna_comparison_base_budget")?;
                if bases > b.max_base_observations {
                    return Err("cdna_comparison_base_budget".into());
                }
                for call in &h.bases {
                    if !in_reference(call.reference_position_0based)
                        || !matches!(call.query_to_reference_strand, '+' | '-')
                        || !call.observed_base.is_ascii()
                        || !call.original_base.is_ascii()
                        || !b"ACGTRYSWKMBDHVN".contains(&(call.observed_base as u8))
                        || !b"ACGTRYSWKMBDHVN".contains(&(call.original_base as u8))
                        || call.phred.is_some_and(|q| q > 93)
                        || call.mapq == Some(255)
                        || call.observed_base
                            != if call.query_to_reference_strand == '-' {
                                complemented(call.original_base)
                            } else {
                                call.original_base
                            }
                    {
                        return Err("cdna_comparison_base_axis_or_quality".into());
                    }
                    let i = h
                        .junctions
                        .partition_point(|s| s.1 <= call.reference_position_0based);
                    if h.junctions.get(i).is_some_and(|&(a, z)| {
                        a <= call.reference_position_0based && call.reference_position_0based < z
                    }) {
                        return Err("cdna_comparison_base_inside_skip".into());
                    }
                }
            }
        }
    }
    Ok(())
}

fn quality_qualified(b: &CdnaBaseObservation, p: &CdnaComparisonPolicy) -> bool {
    canonical(b.observed_base)
        && b.phred.is_some_and(|q| q >= p.minimum_base_quality)
        && b.mapq.is_some_and(|q| q >= p.minimum_mapping_quality)
}

fn calls_for_unit(
    unit: &CdnaUnit,
    policy: &CdnaComparisonPolicy,
    meter: &mut Meter<'_>,
) -> Result<BTreeMap<usize, Call>, String> {
    let mut hypotheses = Vec::<BTreeMap<usize, Vec<&CdnaBaseObservation>>>::new();
    let mut positions = BTreeSet::new();
    for h in unit.hypotheses.iter().filter(|h| h.eligible) {
        let mut indexed = BTreeMap::<usize, Vec<&CdnaBaseObservation>>::new();
        for b in &h.bases {
            meter.charge(1)?;
            positions.insert(b.reference_position_0based);
            indexed
                .entry(b.reference_position_0based)
                .or_default()
                .push(b);
        }
        hypotheses.push(indexed);
    }
    let mut out = BTreeMap::new();
    for pos in positions {
        let mut call = Call {
            base: None,
            quality: None,
            possible: BTreeSet::new(),
            reasons: BTreeSet::new(),
            ambiguous: false,
        };
        let mut qualified = Vec::new();
        if !unit.all_query_views_assessed {
            call.reasons.insert("incomplete_query_views".into());
        }
        for h in &hypotheses {
            meter.charge(1)?;
            match h.get(&pos) {
                None => {
                    call.reasons
                        .insert("site_not_covered_in_all_hypotheses".into());
                }
                Some(values) => {
                    for b in values {
                        meter.charge(1)?;
                        if !canonical(b.observed_base) {
                            call.reasons.insert("noncanonical_base".into());
                        }
                        if b.phred.is_none() {
                            call.reasons.insert("base_quality_unavailable".into());
                        } else if b.phred.is_some_and(|q| q < policy.minimum_base_quality) {
                            call.reasons.insert("low_base_quality".into());
                        }
                        if b.mapq.is_none() {
                            call.reasons.insert("mapping_quality_unavailable".into());
                        } else if b.mapq.is_some_and(|q| q < policy.minimum_mapping_quality) {
                            call.reasons.insert("low_mapping_quality".into());
                        }
                        if quality_qualified(b, policy) {
                            call.possible.insert(b.observed_base);
                            qualified.push(b.phred.expect("qualified"));
                        }
                    }
                }
            }
        }
        if call.possible.len() > 1 {
            call.ambiguous = true;
            call.reasons.insert("conflicting_eligible_bases".into());
        } else if call.reasons.is_empty() && call.possible.len() == 1 {
            call.base = call.possible.first().copied();
            call.quality = qualified.into_iter().min();
        }
        out.insert(pos, call);
    }
    Ok(out)
}

fn annotate(
    position: usize,
    reference: char,
    alternates: &[char],
    transcripts: &[CdnaTranscript],
) -> Vec<CdnaTranscriptSiteAnnotation> {
    transcripts
        .iter()
        .map(|t| {
            let mut offset = 0usize;
            let mut local = None;
            for &(a, z) in &t.exons {
                if a <= position && position < z {
                    local = Some(offset + position - a);
                    break;
                }
                offset += z - a;
            }
            let total: usize = t.exons.iter().map(|&(a, z)| z - a).sum();
            let cdna = local.map(|i| if t.strand == '-' { total - i } else { i + 1 });
            CdnaTranscriptSiteAnnotation {
                transcript_id: t.transcript_id.clone(),
                status: if cdna.is_some() {
                    "present"
                } else {
                    "not_present"
                }
                .into(),
                cdna_position_1based: cdna,
                transcript_reference_base: cdna.map(|_| {
                    if t.strand == '-' {
                        complemented(reference)
                    } else {
                        reference
                    }
                }),
                transcript_alternate_bases: if cdna.is_some() {
                    alternates
                        .iter()
                        .map(|&b| if t.strand == '-' { complemented(b) } else { b })
                        .collect()
                } else {
                    vec![]
                },
            }
        })
        .collect()
}

fn counts(
    sample: &CdnaSample,
    calls: &[BTreeMap<usize, Call>],
    pos: usize,
    reference: char,
    policy: &CdnaComparisonPolicy,
) -> CdnaSiteSampleCounts {
    let mut base_counts = BTreeMap::new();
    let (mut ambiguous, mut unassessed) = (0, 0);
    for row in calls {
        match row.get(&pos) {
            Some(call) if call.base.is_some() => {
                *base_counts.entry(call.base.expect("present")).or_insert(0) += 1;
            }
            Some(call) if call.ambiguous => ambiguous += 1,
            _ => unassessed += 1,
        }
    }
    let callable: usize = base_counts.values().sum();
    let candidate = base_counts
        .iter()
        .filter_map(|(&base, &n)| {
            (base != reference
                && n >= policy.minimum_alternate_units
                && callable > 0
                && (n as f64) / (callable as f64) >= policy.minimum_alternate_fraction)
                .then_some(base)
        })
        .collect();
    CdnaSiteSampleCounts {
        sample_id: sample.sample_id.clone(),
        total_units: sample.units.len(),
        callable_units: callable,
        base_counts,
        ambiguous_units: ambiguous,
        unassessed_units: unassessed,
        candidate_alternate_bases: candidate,
    }
}

fn distance(
    calls: &BTreeMap<usize, Call>,
    sites: &[usize],
    template: &[char],
    meter: &mut Meter<'_>,
) -> Result<(usize, u64, usize), String> {
    let (mut mismatches, mut penalty, mut observed) = (0, 0, 0);
    for (&site, &expected) in sites.iter().zip(template) {
        meter.charge(1)?;
        if let Some(call) = calls.get(&site).filter(|c| c.base.is_some()) {
            observed += 1;
            if call.base != Some(expected) {
                mismatches += 1;
                penalty += u64::from(call.quality.unwrap_or(0).min(40).max(1));
            }
        }
    }
    Ok((mismatches, penalty, observed))
}

fn eligible_distance(mismatches: usize, penalty: u64, policy: &CdnaComparisonPolicy) -> bool {
    mismatches <= policy.maximum_assignment_mismatches
        && penalty <= policy.maximum_assignment_phred_penalty
}

/// A relation of supported local allele patterns on one exact marker scope.
#[derive(Clone)]
struct PhaseFamily {
    sites: Vec<usize>,
    patterns: BTreeMap<Vec<char>, BTreeSet<usize>>,
}

fn local_projection(pattern: &[char], from: &[usize], to: &[usize]) -> Vec<char> {
    to.iter()
        .map(|p| pattern[from.binary_search(p).expect("subset marker scope")])
        .collect()
}

/// A natural join is identifiable only when each shared-allele compatibility
/// component is a star. A many-to-many component has unobserved cross-products.
fn join_phase_families(
    left: &PhaseFamily,
    right: &PhaseFamily,
    max_haplotypes: usize,
    meter: &mut Meter<'_>,
) -> Result<Option<PhaseFamily>, String> {
    meter.charge(left.sites.len() + right.sites.len())?;
    let shared = left
        .sites
        .iter()
        .copied()
        .filter(|p| right.sites.binary_search(p).is_ok())
        .collect::<Vec<_>>();
    if shared.is_empty() {
        return Ok(None);
    }
    let left_contains_right = right
        .sites
        .iter()
        .all(|p| left.sites.binary_search(p).is_ok());
    let right_contains_left = left
        .sites
        .iter()
        .all(|p| right.sites.binary_search(p).is_ok());
    if left_contains_right && right_contains_left {
        // Equal scopes contribute positive classes; non-observation is not a
        // negative constraint. Each retained class already has its own proof.
        let mut family = left.clone();
        for (pattern, proof) in &right.patterns {
            meter.charge(pattern.len())?;
            family
                .patterns
                .entry(pattern.clone())
                .and_modify(|old| {
                    if (proof.len(), proof) < (old.len(), &*old) {
                        *old = proof.clone();
                    }
                })
                .or_insert_with(|| proof.clone());
        }
        if family.patterns.len() > max_haplotypes {
            return Err("cdna_comparison_haplotype_count_budget".into());
        }
        return Ok(Some(family));
    }
    if left_contains_right || right_contains_left {
        let (wide, narrow) = if left_contains_right {
            (left, right)
        } else {
            (right, left)
        };
        let mut projected = BTreeSet::new();
        for pattern in wide.patterns.keys() {
            meter.charge(narrow.sites.len())?;
            projected.insert(local_projection(pattern, &wide.sites, &narrow.sites));
        }
        if narrow.patterns.keys().any(|p| !projected.contains(p)) {
            return Ok(None);
        }
        return Ok(Some(wide.clone()));
    }
    let mut lgroups = BTreeMap::<Vec<char>, Vec<(&Vec<char>, &BTreeSet<usize>)>>::new();
    let mut rgroups = BTreeMap::<Vec<char>, Vec<(&Vec<char>, &BTreeSet<usize>)>>::new();
    for (pattern, proof) in &left.patterns {
        meter.charge(shared.len())?;
        lgroups
            .entry(local_projection(pattern, &left.sites, &shared))
            .or_default()
            .push((pattern, proof));
    }
    for (pattern, proof) in &right.patterns {
        meter.charge(shared.len())?;
        rgroups
            .entry(local_projection(pattern, &right.sites, &shared))
            .or_default()
            .push((pattern, proof));
    }
    if lgroups.keys().ne(rgroups.keys()) {
        return Ok(None);
    }
    for (key, lrows) in &lgroups {
        meter.charge(1)?;
        if lrows.len() > 1 && rgroups[key].len() > 1 {
            return Ok(None);
        }
    }
    let sites = left
        .sites
        .iter()
        .chain(&right.sites)
        .copied()
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect::<Vec<_>>();
    let left_contains_right = right
        .sites
        .iter()
        .all(|p| left.sites.binary_search(p).is_ok());
    let right_contains_left = left
        .sites
        .iter()
        .all(|p| right.sites.binary_search(p).is_ok());
    let mut patterns = BTreeMap::new();
    for (key, lrows) in &lgroups {
        for (lpattern, lproof) in lrows {
            for (rpattern, rproof) in &rgroups[key] {
                meter.charge(sites.len())?;
                let pattern = sites
                    .iter()
                    .map(|p| {
                        left.sites
                            .binary_search(p)
                            .map(|i| lpattern[i])
                            .unwrap_or_else(|_| {
                                rpattern[right.sites.binary_search(p).expect("union scope")]
                            })
                    })
                    .collect::<Vec<_>>();
                // A narrower relation is checked as a constraint but adds no
                // necessary witness to an already supported wider relation.
                let proof = if left_contains_right {
                    (*lproof).clone()
                } else if right_contains_left {
                    (*rproof).clone()
                } else {
                    lproof.union(rproof).copied().collect()
                };
                patterns.insert(pattern, proof);
                if patterns.len() > max_haplotypes {
                    return Err("cdna_comparison_haplotype_count_budget".into());
                }
            }
        }
    }
    Ok(Some(PhaseFamily { sites, patterns }))
}

/// Assemble supported overlaps without requiring a whole-component read.
/// All local classes must be explained; ambiguous couplings remain unresolved.
fn transitive_templates(
    sample: &CdnaSample,
    calls: &[BTreeMap<usize, Call>],
    component: &[usize],
    policy: &CdnaComparisonPolicy,
    meter: &mut Meter<'_>,
) -> Result<BTreeMap<Vec<char>, Vec<CdnaPhaseLink>>, String> {
    let mut local = BTreeMap::<Vec<usize>, BTreeMap<Vec<char>, BTreeSet<usize>>>::new();
    for (ui, row) in calls.iter().enumerate() {
        meter.charge(component.len())?;
        let observed = component
            .iter()
            .filter_map(|&p| row.get(&p).and_then(|c| c.base.map(|b| (p, b))))
            .collect::<Vec<_>>();
        if observed.len() >= 2 {
            local
                .entry(observed.iter().map(|x| x.0).collect())
                .or_default()
                .entry(observed.iter().map(|x| x.1).collect())
                .or_default()
                .insert(ui);
        }
        for (i, &(a, aa)) in observed.iter().enumerate() {
            for &(b, bb) in observed.iter().skip(i + 1) {
                meter.charge(1)?;
                local
                    .entry(vec![a, b])
                    .or_default()
                    .entry(vec![aa, bb])
                    .or_default()
                    .insert(ui);
            }
        }
    }
    let mut proofs = Vec::<CdnaPhaseLink>::new();
    let mut families = Vec::<PhaseFamily>::new();
    for (sites, rows) in local {
        let mut patterns = BTreeMap::new();
        for (alleles, units) in rows {
            meter.charge(units.len() + 1)?;
            if units.len() < policy.minimum_haplotype_units {
                continue;
            }
            let index = proofs.len();
            proofs.push(CdnaPhaseLink {
                sites_0based: sites.clone(),
                alleles: alleles.clone(),
                unit_ids: units
                    .iter()
                    .map(|&ui| sample.units[ui].unit_id.clone())
                    .collect(),
            });
            patterns.insert(alleles, BTreeSet::from([index]));
        }
        if patterns.len() > policy.maximum_haplotypes {
            return Err("cdna_comparison_haplotype_count_budget".into());
        }
        if !patterns.is_empty() {
            families.push(PhaseFamily { sites, patterns });
        }
    }
    let constraints = families.clone();
    while families.len() > 1 {
        meter.poll(0, sample.units.len())?;
        families.sort_by(|a, b| {
            a.sites
                .cmp(&b.sites)
                .then_with(|| a.patterns.keys().cmp(b.patterns.keys()))
        });
        let mut order = Vec::new();
        for (i, left) in families.iter().enumerate() {
            for (j, right) in families.iter().enumerate().skip(i + 1) {
                meter.charge(left.sites.len() + right.sites.len())?;
                let overlap = left
                    .sites
                    .iter()
                    .filter(|p| right.sites.binary_search(p).is_ok())
                    .count();
                if overlap > 0 {
                    order.push((
                        std::cmp::Reverse(overlap),
                        std::cmp::Reverse(left.sites.len() + right.sites.len() - overlap),
                        i,
                        j,
                    ));
                }
            }
        }
        order.sort_unstable();
        let mut merged = None;
        for (_, _, i, j) in order {
            if let Some(family) =
                join_phase_families(&families[i], &families[j], policy.maximum_haplotypes, meter)?
            {
                merged = Some((i, j, family));
                break;
            }
        }
        let Some((i, j, family)) = merged else {
            return Ok(BTreeMap::new());
        };
        families.remove(j);
        families.remove(i);
        families.push(family);
    }
    let Some(family) = families.pop().filter(|f| f.sites == component) else {
        return Ok(BTreeMap::new());
    };
    for constraint in &constraints {
        let mut projected = BTreeSet::new();
        for pattern in family.patterns.keys() {
            meter.charge(constraint.sites.len())?;
            projected.insert(local_projection(pattern, &family.sites, &constraint.sites));
        }
        if constraint.patterns.keys().any(|p| !projected.contains(p)) {
            return Ok(BTreeMap::new());
        }
    }
    Ok(family
        .patterns
        .into_iter()
        .map(|(pattern, links)| {
            (
                pattern,
                links.into_iter().map(|i| proofs[i].clone()).collect(),
            )
        })
        .collect())
}

fn phase_blocks(
    sample: &CdnaSample,
    calls: &[BTreeMap<usize, Call>],
    sites: &[CdnaSiteReport],
    phase_marker_positions: Option<&[usize]>,
    policy: &CdnaComparisonPolicy,
    max_component: usize,
    reference: &CdnaReference,
    transcripts: &[CdnaTranscript],
    meter: &mut Meter<'_>,
) -> Result<Vec<CdnaPhaseBlock>, String> {
    let mut graph = BTreeMap::<usize, BTreeSet<usize>>::new();
    for site in sites {
        if phase_marker_positions
            .is_some_and(|markers| markers.binary_search(&site.position_0based).is_err())
        {
            continue;
        }
        let counts = site
            .samples
            .iter()
            .find(|s| s.sample_id == sample.sample_id)
            .expect("sample counts");
        if counts.callable_units >= policy.minimum_haplotype_units {
            graph.insert(site.position_0based, BTreeSet::new());
        }
    }
    let mut pair_units = BTreeMap::<(usize, usize, char, char), BTreeSet<usize>>::new();
    for (ui, row) in calls.iter().enumerate() {
        let observed = row
            .iter()
            .filter_map(|(&p, c)| c.base.filter(|_| graph.contains_key(&p)).map(|b| (p, b)))
            .collect::<Vec<_>>();
        for (i, &(a, aa)) in observed.iter().enumerate() {
            for &(z, zz) in observed.iter().skip(i + 1) {
                meter.charge(1)?;
                pair_units.entry((a, z, aa, zz)).or_default().insert(ui);
            }
        }
    }
    for ((a, z, _, _), units) in pair_units {
        meter.charge(1)?;
        if units.len() >= policy.minimum_haplotype_units {
            graph.get_mut(&a).expect("known").insert(z);
            graph.get_mut(&z).expect("known").insert(a);
        }
    }
    let site_lookup = sites
        .iter()
        .map(|s| (s.position_0based, s))
        .collect::<BTreeMap<_, _>>();
    let transcript_strands = transcripts
        .iter()
        .map(|t| (t.transcript_id.as_str(), t.strand))
        .collect::<BTreeMap<_, _>>();
    let mut seen = BTreeSet::new();
    let mut blocks = Vec::new();
    for &root in graph.keys() {
        if seen.contains(&root) {
            continue;
        }
        let mut pending = vec![root];
        let mut component = BTreeSet::new();
        while let Some(p) = pending.pop() {
            meter.charge(1)?;
            if !seen.insert(p) {
                continue;
            }
            component.insert(p);
            if component.len() > max_component {
                return Err("cdna_comparison_phase_component_budget".into());
            }
            for &next in &graph[&p] {
                meter.charge(1)?;
                if !seen.contains(&next) {
                    pending.push(next);
                }
            }
        }
        let component = component.into_iter().collect::<Vec<_>>();
        let block_id = hash(&(&reference.sha256, &reference.contig, &component))?;
        let mut patterns = BTreeMap::<Vec<char>, Vec<usize>>::new();
        for (ui, row) in calls.iter().enumerate() {
            meter.charge(component.len())?;
            let pattern = component
                .iter()
                .map(|p| row.get(p).and_then(|c| c.base))
                .collect::<Option<Vec<_>>>();
            if let Some(pattern) = pattern {
                patterns.entry(pattern).or_default().push(ui);
            }
        }
        let templates = if component.len() >= 2 {
            transitive_templates(sample, calls, &component, policy, meter)?
        } else {
            BTreeMap::new()
        };
        if templates.len() > policy.maximum_haplotypes {
            return Err("cdna_comparison_haplotype_count_budget".into());
        }
        let mut block = CdnaPhaseBlock {
            block_id,
            sites_0based: component.clone(),
            status: if component.len() < 2 {
                "allele_only_single_site"
            } else if templates.is_empty() {
                "unresolved_non_identifiable_overlap"
            } else {
                "candidate_templates_uncalibrated"
            }
            .into(),
            candidate_haplotype_count: (!templates.is_empty()).then_some(templates.len()),
            excess_expressed_haplotype_diversity: templates.len() > 2,
            copy_number_status: "not_assessed".into(),
            haplotypes: vec![],
            assignments: vec![],
        };
        for (index, (pattern, links)) in templates.iter().enumerate() {
            meter.charge(component.len().saturating_mul(transcripts.len()))?;
            let alleles = component
                .iter()
                .zip(pattern.iter())
                .map(|(&pos, &base)| {
                    let site = site_lookup[&pos];
                    let mut annotations = site.transcripts.clone();
                    for annotation in &mut annotations {
                        annotation.transcript_alternate_bases = if annotation
                            .cdna_position_1based
                            .is_some()
                            && base != site.reference_base
                        {
                            let oriented =
                                if transcript_strands[annotation.transcript_id.as_str()] == '-' {
                                    complemented(base)
                                } else {
                                    base
                                };
                            vec![oriented]
                        } else {
                            vec![]
                        };
                    }
                    CdnaHaplotypeAllele {
                        position_0based: pos,
                        reference_base: site.reference_base,
                        observed_base: base,
                        transcript_annotations: annotations,
                    }
                })
                .collect::<Vec<_>>();
            block.haplotypes.push(CdnaHaplotype {
                haplotype_id: format!("H{}", index + 1),
                fingerprint: hash(&(&reference.sha256, &reference.contig, &component, pattern))?,
                alleles,
                exact_seed_units: patterns
                    .get(pattern)
                    .into_iter()
                    .flatten()
                    .map(|&ui| sample.units[ui].unit_id.clone())
                    .collect(),
                phase_method: if patterns
                    .get(pattern)
                    .is_some_and(|u| u.len() >= policy.minimum_haplotype_units)
                {
                    "direct_complete_vectors"
                } else {
                    "transitive_supported_local_patterns"
                }
                .into(),
                support_links: if patterns
                    .get(pattern)
                    .is_some_and(|u| u.len() >= policy.minimum_haplotype_units)
                {
                    vec![]
                } else {
                    links.clone()
                },
                assigned_units: 0,
                error_tolerant_assigned_units: 0,
            });
        }
        for (ui, row) in calls.iter().enumerate() {
            let observed = component
                .iter()
                .filter(|p| row.get(p).is_some_and(|c| c.base.is_some()))
                .count();
            let mut best = None;
            let mut matches = vec![];
            if observed >= policy.minimum_assignment_sites {
                for (index, (pattern, _)) in templates.iter().enumerate() {
                    let (m, c, _) = distance(row, &component, pattern, meter)?;
                    if !eligible_distance(m, c, policy) {
                        continue;
                    }
                    let score = (c, m);
                    if best.is_none_or(|v| score < v) {
                        best = Some(score);
                        matches.clear();
                    }
                    if best == Some(score) {
                        matches.push(index);
                    }
                }
            }
            let status = if observed < policy.minimum_assignment_sites {
                "insufficient_callable_sites"
            } else if templates.is_empty() {
                "unresolved_no_templates"
            } else if matches.is_empty() {
                "unsupported_pattern"
            } else if matches.len() > 1 {
                "ambiguous"
            } else if best.is_some_and(|(_, m)| m > 0) {
                "assigned_error_tolerant"
            } else {
                "assigned_exact"
            };
            if matches.len() == 1 {
                let h = &mut block.haplotypes[matches[0]];
                h.assigned_units += 1;
                h.error_tolerant_assigned_units += usize::from(best.is_some_and(|(_, m)| m > 0));
            }
            block.assignments.push(CdnaHaplotypeAssignment {
                unit_id: sample.units[ui].unit_id.clone(),
                status: status.into(),
                callable_sites: observed,
                compatible_haplotype_ids: matches
                    .iter()
                    .map(|&i| block.haplotypes[i].haplotype_id.clone())
                    .collect(),
                mismatch_count: best.map(|(_, m)| m),
                phred_penalty: best.map(|(c, _)| c),
            });
        }
        blocks.push(block);
    }
    Ok(blocks)
}

fn switches(
    sample: &CdnaSample,
    calls: &[BTreeMap<usize, Call>],
    blocks: &[CdnaPhaseBlock],
    policy: &CdnaComparisonPolicy,
    remaining: &mut usize,
    meter: &mut Meter<'_>,
) -> Result<Vec<CdnaSwitchCandidate>, String> {
    let mut out = vec![];
    for block in blocks {
        for (ui, unit) in sample.units.iter().enumerate() {
            meter.poll(ui, sample.units.len())?;
            let row = &calls[ui];
            let templates = block
                .haplotypes
                .iter()
                .filter(|h| {
                    if h.phase_method == "direct_complete_vectors" {
                        h.exact_seed_units
                            .len()
                            .saturating_sub(usize::from(h.exact_seed_units.contains(&unit.unit_id)))
                            >= policy.minimum_haplotype_units
                    } else {
                        !h.support_links.is_empty()
                            && h.support_links.iter().all(|link| {
                                link.unit_ids.len().saturating_sub(usize::from(
                                    link.unit_ids.contains(&unit.unit_id),
                                )) >= policy.minimum_haplotype_units
                            })
                    }
                })
                .collect::<Vec<_>>();
            if templates.len() < 2 {
                continue;
            }
            let mut whole = None;
            for h in &templates {
                let vector = h
                    .alleles
                    .iter()
                    .map(|a| a.observed_base)
                    .collect::<Vec<_>>();
                let (m, c, n) = distance(row, &block.sites_0based, &vector, meter)?;
                if n >= 2 && whole.is_none_or(|old| (c, m) < old) {
                    whole = Some((c, m));
                }
            }
            let mut best = None;
            let mut explanations = vec![];
            for (ai, a) in templates.iter().enumerate() {
                for (bi, b) in templates.iter().enumerate() {
                    if ai == bi {
                        continue;
                    }
                    let mut markers = vec![];
                    let (mut shared_mismatches, mut shared_penalty) = (0usize, 0u64);
                    for (i, &pos) in block.sites_0based.iter().enumerate() {
                        meter.charge(1)?;
                        if let Some(call) = row.get(&pos).filter(|c| c.base.is_some()) {
                            let observed = call.base.expect("observed");
                            let quality = u64::from(call.quality.unwrap_or(0).min(40).max(1));
                            if a.alleles[i].observed_base != b.alleles[i].observed_base {
                                markers.push((
                                    pos,
                                    observed,
                                    quality,
                                    a.alleles[i].observed_base,
                                    b.alleles[i].observed_base,
                                ));
                            } else if observed != a.alleles[i].observed_base {
                                // A mosaic cannot explain disagreement at a site shared by
                                // both templates. Use the same callable scope as the whole
                                // template alternative; only support counts are diagnostic-only.
                                shared_mismatches += 1;
                                shared_penalty += quality;
                            }
                        }
                    }
                    if markers.len() < 2 {
                        continue;
                    }
                    let mut prefix_a = vec![(0usize, 0u64, 0usize)];
                    let mut prefix_b = vec![(0usize, 0u64, 0usize)];
                    for &(_, observed, q, aa, bb) in &markers {
                        meter.charge(2)?;
                        let (m, c, s) = *prefix_a.last().expect("prefix");
                        prefix_a.push((
                            m + usize::from(observed != aa),
                            c + if observed != aa { q } else { 0 },
                            s + usize::from(observed == aa),
                        ));
                        let (m, c, s) = *prefix_b.last().expect("prefix");
                        prefix_b.push((
                            m + usize::from(observed != bb),
                            c + if observed != bb { q } else { 0 },
                            s + usize::from(observed == bb),
                        ));
                    }
                    let total_b = *prefix_b.last().expect("prefix");
                    for split in 1..markers.len() {
                        meter.charge(1)?;
                        let left = prefix_a[split];
                        let right = (
                            total_b.0 - prefix_b[split].0,
                            total_b.1 - prefix_b[split].1,
                            total_b.2 - prefix_b[split].2,
                        );
                        if left.2 == 0 || right.2 == 0 {
                            continue;
                        }
                        let score = (
                            left.1 + right.1 + shared_penalty,
                            left.0 + right.0 + shared_mismatches,
                        );
                        if !eligible_distance(score.1, score.0, policy)
                            || best.is_some_and(|old| score > old)
                        {
                            continue;
                        }
                        if best.is_none_or(|old| score < old) {
                            best = Some(score);
                            explanations.clear();
                        }
                        let left_pos = markers[split - 1].0;
                        let right_pos = markers[split].0;
                        let mut junctions = BTreeMap::new();
                        let mut query_markers = BTreeMap::new();
                        for h in unit.hypotheses.iter().filter(|h| h.eligible) {
                            meter.charge(h.junctions.len() + h.bases.len())?;
                            let supported = h
                                .junctions
                                .iter()
                                .copied()
                                .filter(|&(s, e)| left_pos < s && e <= right_pos)
                                .collect::<Vec<_>>();
                            let qleft = h
                                .bases
                                .iter()
                                .filter(|c| {
                                    c.reference_position_0based == left_pos
                                        && quality_qualified(c, policy)
                                })
                                .map(|c| c.query_position_0based)
                                .collect::<BTreeSet<_>>();
                            let qright = h
                                .bases
                                .iter()
                                .filter(|c| {
                                    c.reference_position_0based == right_pos
                                        && quality_qualified(c, policy)
                                })
                                .map(|c| c.query_position_0based)
                                .collect::<BTreeSet<_>>();
                            let mut pairs = BTreeSet::new();
                            for &l in &qleft {
                                for &r in &qright {
                                    meter.charge(1)?;
                                    pairs.insert((l, r));
                                }
                            }
                            if !pairs.is_empty() {
                                query_markers
                                    .insert(h.hypothesis_id.clone(), pairs.into_iter().collect());
                                if !supported.is_empty() {
                                    junctions.insert(h.hypothesis_id.clone(), supported);
                                }
                            }
                        }
                        explanations.push(CdnaSwitchExplanation {
                            left_haplotype_id: a.haplotype_id.clone(),
                            right_haplotype_id: b.haplotype_id.clone(),
                            marker_bracket_0based: (left_pos, right_pos + 1),
                            left_informative_markers: left.2,
                            right_informative_markers: right.2,
                            mismatch_count: score.1,
                            phred_penalty: score.0,
                            same_hypothesis_junctions: junctions,
                            query_marker_positions: query_markers,
                        });
                    }
                }
            }
            if explanations.is_empty() {
                continue;
            }
            if explanations.len() > *remaining {
                return Err("cdna_comparison_switch_explanation_budget".into());
            }
            *remaining -= explanations.len();
            explanations.sort_by(|a, b| {
                (
                    &a.left_haplotype_id,
                    &a.right_haplotype_id,
                    a.marker_bracket_0based,
                )
                    .cmp(&(
                        &b.left_haplotype_id,
                        &b.right_haplotype_id,
                        b.marker_bracket_0based,
                    ))
            });
            let competing = whole.is_some_and(|score| score <= best.expect("explanation score"));
            let strong = explanations
                .iter()
                .any(|e| e.left_informative_markers >= 2 && e.right_informative_markers >= 2);
            let junction = explanations
                .iter()
                .any(|e| !e.same_hypothesis_junctions.is_empty());
            out.push(CdnaSwitchCandidate {
                unit_id: unit.unit_id.clone(),
                block_id: block.block_id.clone(),
                status: if competing {
                    "whole_haplotype_alternative"
                } else if explanations.len() > 1 {
                    "ambiguous_switch_explanations"
                } else if !strong {
                    "low_information_switch_candidate"
                } else if junction {
                    "junction_associated_switch_candidate"
                } else {
                    "haplotype_switch_candidate"
                }
                .into(),
                coordinate_order: "genomic_forward".into(),
                whole_template_mismatch_count: whole.map(|(_, m)| m),
                competing_whole_template: competing,
                explanations,
                mechanism_status: "not_established".into(),
            });
        }
    }
    Ok(out)
}

fn junction_support(
    sample: &CdnaSample,
    blocks: &[CdnaPhaseBlock],
    meter: &mut Meter<'_>,
) -> Result<Vec<CdnaHaplotypeJunctionSupport>, String> {
    let units = sample
        .units
        .iter()
        .map(|u| (u.unit_id.as_str(), u))
        .collect::<BTreeMap<_, _>>();
    let mut result = vec![];
    for block in blocks {
        let mut groups = BTreeMap::<(String, Vec<(usize, usize)>), (usize, usize)>::new();
        for assignment in &block.assignments {
            if assignment.compatible_haplotype_ids.len() != 1 {
                continue;
            }
            let unit = units[assignment.unit_id.as_str()];
            let mut chains = BTreeSet::new();
            for h in unit.hypotheses.iter().filter(|h| h.eligible) {
                meter.charge(h.junctions.len() + 1)?;
                chains.insert(h.junctions.clone());
            }
            let exclusive = chains.len() == 1;
            for chain in chains {
                let entry = groups
                    .entry((assignment.compatible_haplotype_ids[0].clone(), chain))
                    .or_default();
                entry.0 += 1;
                entry.1 += usize::from(exclusive);
            }
        }
        for ((haplotype_id, junctions), (possible_units, exclusive_units)) in groups {
            result.push(CdnaHaplotypeJunctionSupport {
                block_id: block.block_id.clone(),
                haplotype_id,
                junctions,
                possible_units,
                exclusive_units,
            });
        }
    }
    Ok(result)
}

/// Compare cDNA SNP evidence for explicit samples in one exact reference locus.
///
/// Missing mapping/base-quality alternatives withhold unit consensus. Templates may phase through
/// identifiable supported overlaps; a whole-block read is not required or invented.
pub fn compare_cdna_alleles(
    mut request: CdnaComparisonRequest,
    progress: &mut dyn FnMut(usize, usize) -> bool,
) -> Result<CdnaComparisonReport, String> {
    validate(&request)?;
    // Normalize presentation order without changing any evidence or selecting an alternative.
    request
        .transcripts
        .sort_by(|a, b| a.transcript_id.cmp(&b.transcript_id));
    request
        .samples
        .sort_by(|a, b| a.sample_id.cmp(&b.sample_id));
    request.supplied_sites.sort_unstable();
    if let Some(markers) = &mut request.phase_marker_positions {
        markers.sort_unstable();
    }
    for sample in &mut request.samples {
        sample.units.sort_by(|a, b| a.unit_id.cmp(&b.unit_id));
        for unit in &mut sample.units {
            unit.hypotheses
                .sort_by(|a, b| a.hypothesis_id.cmp(&b.hypothesis_id));
            for h in &mut unit.hypotheses {
                h.bases.sort_by_key(|b| {
                    (
                        b.reference_position_0based,
                        b.query_position_0based,
                        b.query_to_reference_strand,
                        b.original_base,
                        b.observed_base,
                        b.phred,
                        b.mapq,
                    )
                });
            }
        }
    }
    let mut meter = Meter {
        work: 0,
        limit: request.budgets.max_comparisons,
        progress,
    };
    meter.poll(0, request.samples.len())?;
    let mut sample_calls = vec![];
    let mut positions = BTreeSet::new();
    for (si, sample) in request.samples.iter().enumerate() {
        let mut calls = vec![];
        for (ui, unit) in sample.units.iter().enumerate() {
            meter.poll(ui, sample.units.len())?;
            let mut row = calls_for_unit(unit, &request.policy, &mut meter)?;
            for (&p, call) in &mut row {
                positions.insert(p);
                if !canonical(
                    request.reference.sequence.as_bytes()[p - request.reference.start_0based]
                        as char,
                ) {
                    call.base = None;
                    call.quality = None;
                    call.reasons.insert("reference_base_noncanonical".into());
                }
            }
            calls.push(row);
        }
        sample_calls.push(calls);
        meter.poll(si + 1, request.samples.len())?;
    }
    let mut selected = request
        .supplied_sites
        .iter()
        .copied()
        .collect::<BTreeSet<_>>();
    let supplied = selected.clone();
    let mut assessed = 0usize;
    let total_units: usize = request.samples.iter().map(|s| s.units.len()).sum();
    for &pos in &positions {
        meter.charge(total_units)?;
        let reference =
            request.reference.sequence.as_bytes()[pos - request.reference.start_0based] as char;
        let summaries = request
            .samples
            .iter()
            .zip(&sample_calls)
            .map(|(s, c)| counts(s, c, pos, reference, &request.policy))
            .collect::<Vec<_>>();
        if summaries.iter().any(|s| s.callable_units > 0) {
            assessed += 1;
        }
        if canonical(reference)
            && summaries
                .iter()
                .any(|s| !s.candidate_alternate_bases.is_empty())
        {
            selected.insert(pos);
        }
        if selected.len() > request.budgets.max_sites {
            return Err("cdna_comparison_site_budget".into());
        }
    }
    let exon_count: usize = request.transcripts.iter().map(|t| t.exons.len()).sum();
    let mut site_reports = vec![];
    for &pos in &selected {
        meter.charge(total_units + 2 * exon_count)?;
        let reference =
            request.reference.sequence.as_bytes()[pos - request.reference.start_0based] as char;
        let summaries = request
            .samples
            .iter()
            .zip(&sample_calls)
            .map(|(s, c)| counts(s, c, pos, reference, &request.policy))
            .collect::<Vec<_>>();
        let alt = summaries
            .iter()
            .flat_map(|s| s.candidate_alternate_bases.iter().copied())
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect::<Vec<_>>();
        site_reports.push(CdnaSiteReport {
            position_0based: pos,
            reference_base: reference,
            supplied: supplied.contains(&pos),
            transcripts: annotate(pos, reference, &alt, &request.transcripts),
            candidate_alternate_bases: alt,
            samples: summaries,
        });
    }
    // Discard nonselected all-reference position matrices after complete discovery accounting.
    for calls in &mut sample_calls {
        for row in calls {
            row.retain(|p, _| selected.contains(p));
        }
    }
    let mut samples = vec![];
    let mut remaining = request.budgets.max_switch_explanations;
    for (sample, calls) in request.samples.iter().zip(&sample_calls) {
        let mut site_calls = vec![];
        let mut observations = vec![];
        for (ui, unit) in sample.units.iter().enumerate() {
            meter.poll(ui, sample.units.len())?;
            for site in &site_reports {
                meter.charge(1)?;
                let mut call = calls[ui]
                    .get(&site.position_0based)
                    .cloned()
                    .unwrap_or_else(Call::uncovered);
                if !unit.all_query_views_assessed {
                    call.reasons.insert("incomplete_query_views".into());
                }
                if !unit.hypotheses.iter().any(|h| h.eligible) {
                    call.reasons.insert("no_eligible_hypothesis".into());
                }
                site_calls.push(CdnaUnitAlleleCall {
                    unit_id: unit.unit_id.clone(),
                    position_0based: site.position_0based,
                    status: if call.ambiguous {
                        CdnaCallStatus::Ambiguous
                    } else {
                        match call.base {
                            Some(b) if b == site.reference_base => {
                                CdnaCallStatus::ReferenceSupported
                            }
                            Some(_) => CdnaCallStatus::AlternateSupported,
                            None => CdnaCallStatus::Unassessed,
                        }
                    },
                    base: call.base,
                    minimum_phred: call.quality,
                    possible_bases: call.possible.into_iter().collect(),
                    reasons: call.reasons.into_iter().collect(),
                });
            }
            for h in &unit.hypotheses {
                for base in &h.bases {
                    meter.charge(1)?;
                    if selected.contains(&base.reference_position_0based) {
                        observations.push(CdnaSiteObservation {
                            unit_id: unit.unit_id.clone(),
                            hypothesis_id: h.hypothesis_id.clone(),
                            eligible: h.eligible,
                            quality_qualified: quality_qualified(base, &request.policy),
                            observation: base.clone(),
                        });
                    }
                }
            }
        }
        let blocks = phase_blocks(
            sample,
            calls,
            &site_reports,
            request.phase_marker_positions.as_deref(),
            &request.policy,
            request.budgets.max_phase_sites_per_block,
            &request.reference,
            &request.transcripts,
            &mut meter,
        )?;
        let candidates = switches(
            sample,
            calls,
            &blocks,
            &request.policy,
            &mut remaining,
            &mut meter,
        )?;
        let junctions = junction_support(sample, &blocks, &mut meter)?;
        samples.push(CdnaSampleReport {
            sample_id: sample.sample_id.clone(),
            label: sample.label.clone(),
            source_sha256: sample.source_sha256.clone(),
            analysis_context_sha256: sample.analysis_context_sha256.clone(),
            total_units: sample.units.len(),
            units_with_incomplete_views: sample
                .units
                .iter()
                .filter(|u| !u.all_query_views_assessed)
                .count(),
            units_without_complete_search_declaration: sample
                .units
                .iter()
                .filter(|u| !u.search_declared_complete)
                .count(),
            site_calls,
            site_observations: observations,
            phase_blocks: blocks,
            switch_candidates: candidates,
            haplotype_junction_support: junctions,
        });
    }
    let mut differences = vec![];
    for (ai, a) in request.samples.iter().enumerate() {
        for (bi, b) in request.samples.iter().enumerate().skip(ai + 1) {
            for site in &site_reports {
                meter.charge(1)?;
                let aa = &site.samples[ai];
                let bb = &site.samples[bi];
                differences.push(CdnaSampleDifference {
                    sample_a: a.sample_id.clone(),
                    sample_b: b.sample_id.clone(),
                    position_0based: site.position_0based,
                    reference_base: site.reference_base,
                    sample_a_candidate_bases: aa.candidate_alternate_bases.clone(),
                    sample_b_candidate_bases: bb.candidate_alternate_bases.clone(),
                    status: if aa.callable_units == 0 || bb.callable_units == 0 {
                        "unassessed_comparison"
                    } else if aa.candidate_alternate_bases == bb.candidate_alternate_bases {
                        "same_candidate_allele_set"
                    } else {
                        "different_observed_candidate_allele_sets"
                    }
                    .into(),
                });
            }
        }
    }
    meter.poll(request.samples.len(), request.samples.len())?;
    let mut report=CdnaComparisonReport {
        schema:CDNA_COMPARISON_SCHEMA.into(),content_sha256:String::new(),request_sha256:hash(&request)?,gene:request.gene.clone(),
        reference:CdnaReferenceSummary {assembly:request.reference.assembly.clone(),contig:request.reference.contig.clone(),start_0based:request.reference.start_0based,
            end_0based_exclusive:request.reference.start_0based+request.reference.sequence.len(),sequence_sha256:request.reference.sha256.clone()},
        transcripts:request.transcripts.clone(),policy:request.policy.clone(),policy_status:"uncalibrated_pilot_thresholds".into(),
        discovery_assessed_positions:assessed,phase_marker_positions:request.phase_marker_positions.clone(),sites:site_reports,samples,sample_differences:differences,comparisons:meter.work,
        warnings:vec![
            "RNA sequence-variant candidates are not proven DNA genotypes; editing, mapping and RT/PCR errors remain alternatives. No allele is declared genomically absent from cDNA non-observation.".into(),
            "Consensus requires callable agreement of all eligible mapping/basecall alternatives. Incomplete views withhold consensus; incomplete search declarations remain a limitation without erasing local positive evidence.".into(),
            "Haplotypes are empirical supported local-pattern classes assembled transitively through uniquely identifiable overlaps, with quality-weighted assignments under uncalibrated pilot thresholds. Many-to-many overlap couplings remain unresolved; a whole-block read is not required. This is not calibrated arbitrary-ploidy statistical inference.".into(),
            "H1/H2 display names are local to a sample/block. Compare exact reference/position/base fingerprints, never local labels; do not concatenate disconnected blocks or assume an unobserved second allele.".into(),
            "Additional supported RNA template classes are excess-diversity candidates, not copy number. K caps fail explicitly; they do not prove absence of additional classes.".into(),
            "Switches use leave-one-unit-out template support and retain weak/competing explanations. A same-hypothesis junction is required for a splicing-associated label; no trans-splicing mechanism is established.".into(),
            "Junction support rows are nonadditive observed skip sets, not new physical-chain identities. Read units are not independently verified molecules or biological replicates.".into(),
            "Discovery counts only supplied normalized supported bases; they do not establish whole-genome or whole-transcript completeness. Provide sites explicitly for callable reference-negative comparisons.".into(),
            "Raw input-file and sample identities are authenticated by the calling service; this pure engine validates the declared normalized evidence and locus-sequence hash.".into(),
        ],
    };
    report.content_sha256 = hash(&report)?;
    if serde_json::to_vec(&report)
        .map_err(|e| e.to_string())?
        .len()
        > request.budgets.max_output_bytes
    {
        return Err("cdna_comparison_output_budget".into());
    }
    Ok(report)
}

#[cfg(test)]
#[path = "cdna_comparison_tests.rs"]
mod tests;
