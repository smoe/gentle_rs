//! Pure exon-locus evidence before gene assignment. Scores rank geometry, not biological probability.

use super::*;
use crate::engine::PrimerPayloadAnnotation;
use std::cmp::Reverse;

mod presentation;
pub use presentation::{Page, Projection as ViewProjection, View, page, scene, tsv};
#[cfg(test)]
mod review_tests;
#[cfg(test)]
mod tests;

pub const SCHEMA: &str = "gentle.rna_read_end_locus.v1";
pub const RULE_VERSION: &str = "gentle.locus_assignment.lexicographic.v1";

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Region {
    pub contig: String,
    pub contig_length: usize,
    pub exon: Interval,
    pub flank_left: usize,
    pub flank_right: usize,
    pub anchor_transcript_id: Option<String>,
}
impl Region {
    pub fn window(&self) -> Result<Interval, EngineError> {
        if !token(&self.contig)
            || self
                .anchor_transcript_id
                .as_ref()
                .is_some_and(|id| !token(id))
            || self.exon.0 >= self.exon.1
            || self.exon.1 > self.contig_length
            || self.contig_length > i32::MAX as usize
        {
            return Err(fail("locus_region"));
        }
        Ok((
            self.exon.0.saturating_sub(self.flank_left),
            self.exon
                .1
                .checked_add(self.flank_right)
                .ok_or_else(|| fail("locus_region_overflow"))?
                .min(self.contig_length),
        ))
    }
}

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Bindings {
    pub report_content_sha256: String,
    pub data_sha256: String,
    pub annotation_sha256: Option<String>,
    pub model_sha256: Option<String>,
}
pub struct Input<'a> {
    pub sample_id: &'a str,
    pub bindings: Bindings,
    pub analysis: &'a ObservedSamAnalysis,
    pub annotation: Option<&'a PrimerPayloadAnnotation>,
    pub model: Option<&'a ModelDocument>,
    pub original_bases: Option<oriented::Input<'a>>,
    pub diagnostics: Vec<String>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Class {
    ExactTranscriptMatch,
    CompatibleTranscripts,
    NovelIsoformCompatible,
    AntisenseCompatible,
    UninformativeOverlap,
    Chimeric,
    Unresolved,
}

#[derive(Debug, Clone, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Score {
    pub structural_conflicts: usize,
    pub matched_junctions: usize,
    pub full_chain: bool,
    pub discriminating_exon_bp: usize,
    pub endpoint_matches: usize,
    pub endpoint_distance: usize,
}
impl Score {
    fn strong(&self) -> (Reverse<usize>, usize, usize) {
        (
            Reverse(self.structural_conflicts),
            self.matched_junctions,
            self.discriminating_exon_bp,
        )
    }
    fn rank(&self) -> (Reverse<usize>, usize, bool, usize, usize, Reverse<usize>) {
        (
            Reverse(self.structural_conflicts),
            self.matched_junctions,
            self.full_chain,
            self.discriminating_exon_bp,
            self.endpoint_matches,
            Reverse(self.endpoint_distance),
        )
    }
    fn informative(&self) -> bool {
        self.structural_conflicts == 0
            && (self.matched_junctions > 0 || self.discriminating_exon_bp > 0)
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Model {
    pub key: String,
    pub sample_id: String,
    pub transcript_id: String,
    pub contig: String,
    pub gene_ids: Vec<String>,
    pub strand: char,
    pub exons: Vec<Interval>,
    pub annotation_release: String,
    pub annotation_sha256: String,
    pub model_sha256: Option<String>,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Hypothesis {
    pub index: usize,
    pub observation_ids: Vec<String>,
    pub chimeric: bool,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CandidatePath {
    pub hypothesis_index: usize,
    pub class: Class,
    pub score: Score,
    pub evidence: Vec<String>,
    pub conflicts: Vec<String>,
    pub transcript_intervals: Vec<Interval>,
    pub utr5: Metric,
    pub utr3: Metric,
    pub oriented: oriented::Evidence,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Assignment {
    pub model: Model,
    pub class: Class,
    pub score: Score,
    pub score_upper: Score,
    pub preferred: bool,
    pub evidence: Vec<String>,
    pub conflicts: Vec<String>,
    pub paths: Vec<CandidatePath>,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct MoleculeStrand {
    pub status: String,
    pub strand: Option<char>,
    pub reasons: Vec<String>,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Read {
    pub sample_id: String,
    pub read_key: String,
    pub unit: ObservedSamQueryUnitKey,
    pub membership: Membership,
    pub query_length: Option<usize>,
    pub bindings: Bindings,
    pub observations: Vec<Observation>,
    pub rows: Vec<ObservedSamRowAssessment>,
    pub alignment_directions: Vec<char>,
    pub supplementary_present: bool,
    pub search_declared_complete: bool,
    pub views_declared_complete: bool,
    pub all_query_views_assessed: bool,
    pub hypotheses: Vec<Hypothesis>,
    /// Alignment coverage in fixed genomic coordinates; transcript ends intentionally unset.
    pub genomic_basecalls: Vec<oriented::Evidence>,
    pub assignments: Vec<Assignment>,
    pub molecule_strand: MoleculeStrand,
    pub unavailable: Vec<String>,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum Membership {
    AdmittedHypothesis,
    NonAdmittedObservationOnly,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Observation {
    pub observation_id: String,
    pub on_optimal_path: bool,
    pub source_ordinal: u64,
    pub contig: String,
    pub query_to_reference_strand: char,
    pub blocks: Vec<SupportedBlock>,
    pub skips: Vec<Interval>,
    pub clipping: Option<super::super::ObservedSamClipping>,
}
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SampleSummary {
    pub sample_id: String,
    pub input_units: usize,
    pub locus_units: usize,
    pub non_admitted_observation_only_units: usize,
    pub excluded: Vec<(ObservedSamQueryUnitKey, String)>,
    pub diagnostics: Vec<String>,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Report {
    pub schema: String,
    pub identity_scope: String,
    pub source_id: String,
    pub region: Region,
    pub window: Interval,
    pub scoring_version: String,
    pub scoring_rules: Vec<String>,
    pub effective_limits: Limits,
    /// Deterministic consumed charge units, not wall time, CPU time or allocated bytes.
    pub work: usize,
    pub models: Vec<Model>,
    pub reads: Vec<Read>,
    pub samples: Vec<SampleSummary>,
    pub nonclaims: Vec<String>,
}

fn overlap(a: Interval, b: Interval) -> bool {
    a.0 < b.1 && b.0 < a.1
}
fn introns(exons: &[Interval]) -> BTreeSet<Interval> {
    exons.windows(2).map(|p| (p[0].1, p[1].0)).collect()
}
fn covered(
    ranges: &[Interval],
    exons: &[Interval],
    meter: &mut Meter<'_>,
) -> Result<usize, EngineError> {
    let mut out = 0usize;
    for &(a, z) in ranges {
        for &(l, r) in exons {
            meter.charge(1)?;
            out = out
                .checked_add(z.min(r).saturating_sub(a.max(l)))
                .ok_or_else(|| fail("locus_count_overflow"))?;
        }
    }
    Ok(out)
}
fn other_exons(
    model: &Model,
    models: &[Model],
    meter: &mut Meter<'_>,
) -> Result<Vec<Interval>, EngineError> {
    let mut other = vec![];
    for m in models
        .iter()
        .filter(|m| m.exons != model.exons || m.strand != model.strand)
    {
        meter.charge(m.exons.len())?;
        if other
            .len()
            .checked_add(m.exons.len())
            .is_none_or(|n| n > meter.limits.max_intervals)
        {
            return Err(fail("locus_interval_budget"));
        }
        other.extend_from_slice(&m.exons);
    }
    if other.is_empty() {
        return Ok(other);
    }
    meter.merge(other)
}
fn exclusive(
    ranges: &[Interval],
    model: &Model,
    other: &[Interval],
    meter: &mut Meter<'_>,
) -> Result<usize, EngineError> {
    if other.is_empty() {
        return Ok(0);
    }
    let mut count = 0;
    let mut own = vec![];
    for &(a, z) in ranges {
        for &(l, r) in &model.exons {
            meter.charge(1)?;
            if a.max(l) < z.min(r) {
                if own.len() >= meter.limits.max_intervals {
                    return Err(fail("locus_interval_budget"));
                }
                own.push((a.max(l), z.min(r)));
            }
        }
    }
    let own = meter.merge(own)?;
    for &(a, z) in &own {
        let mut at = a;
        for &(l, r) in other {
            meter.charge(1)?;
            if r <= at || l >= z {
                continue;
            }
            count += l.min(z).saturating_sub(at);
            at = at.max(r).min(z);
        }
        count += z.saturating_sub(at);
    }
    Ok(count)
}
fn chimera(path: &[&ObservedSamObservation]) -> bool {
    let contigs: BTreeSet<_> = path.iter().map(|o| &o.contig).collect();
    let strands: BTreeSet<_> = path.iter().map(|o| o.query_to_reference_strand).collect();
    if contigs.len() > 1 || strands.len() > 1 {
        return true;
    }
    let mut blocks: Vec<_> = path.iter().flat_map(|o| o.blocks.iter()).collect();
    blocks.sort_by_key(|b| b.query);
    blocks.windows(2).any(|p| {
        p[0].query.1 <= p[1].query.0
            && if strands.first() == Some(&'+') {
                p[0].target.1 > p[1].target.0
            } else {
                p[1].target.1 > p[0].target.0
            }
    })
}

fn candidate_path(
    source: &str,
    unit: &ObservedSamUnitAssessment,
    index: usize,
    path: &[&ObservedSamObservation],
    model: &Model,
    models: &[Model],
    other_exons: &[Interval],
    cds: Option<&TranscriptModel>,
    anchor: Option<char>,
    variants: &[oriented::Variant],
    originals: Option<oriented::Input<'_>>,
    policy: oriented::Policy,
    meter: &mut Meter<'_>,
) -> Result<CandidatePath, EngineError> {
    let compact: Vec<_> = path
        .iter()
        .filter(|o| o.contig == model.contig)
        .map(|o| pileup::Observation {
            observation_id: o.observation_id.clone(),
            row_ordinal: o.source_ordinal,
            blocks: o.blocks.clone(),
            connectors: vec![],
            query_to_reference_strand: o.query_to_reference_strand,
        })
        .collect();
    let oriented = oriented::project(
        source,
        &unit.unit,
        model.strand,
        &compact,
        variants,
        originals,
        policy,
        meter,
    )?;
    let mut conflicts = vec![];
    let is_chimera = chimera(path);
    if is_chimera {
        conflicts.push("chimeric_hypothesis".into());
    }
    let off_contig = path.iter().any(|o| o.contig != model.contig);
    if off_contig {
        conflicts.push("observation_on_other_contig".into());
    }
    let ranges = meter.merge(
        path.iter()
            .filter(|o| o.contig == model.contig)
            .flat_map(|o| o.blocks.iter().map(|b| b.target))
            .collect(),
    )?;
    let skips: BTreeSet<_> = path
        .iter()
        .filter(|o| o.contig == model.contig)
        .flat_map(|o| o.skips.iter().copied())
        .collect();
    let known = introns(&model.exons);
    let matched = skips.intersection(&known).count();
    let unknown: Vec<_> = skips.difference(&known).copied().collect();
    let boundary_pair = unknown.iter().all(|&(a, z)| {
        model.exons.iter().any(|e| e.1 == a) && model.exons.iter().any(|e| e.0 == z)
    });
    let bases: usize = ranges.iter().map(|i| i.1 - i.0).sum();
    let outside = bases.saturating_sub(covered(&ranges, &model.exons, meter)?);
    if outside > 0 {
        conflicts.push("aligned_bases_outside_model_exons".into());
    }
    if !unknown.is_empty() {
        conflicts.push("junction_not_in_model_chain".into());
    }
    if path.is_empty() {
        conflicts.push("no_legal_observation_path".into());
    }
    let span = (
        ranges.first().map_or(0, |i| i.0),
        ranges.last().map_or(0, |i| i.1),
    );
    let ends = (model.exons[0].0, model.exons.last().unwrap().1);
    let structural_conflicts = unknown.len()
        + usize::from(outside > 0)
        + usize::from(is_chimera)
        + usize::from(off_contig)
        + usize::from(path.is_empty());
    let full_chain = structural_conflicts == 0 && !known.is_empty() && skips == known;
    let score = Score {
        structural_conflicts,
        matched_junctions: matched,
        full_chain,
        discriminating_exon_bp: exclusive(&ranges, model, other_exons, meter)?,
        endpoint_matches: usize::from(span.0 == ends.0) + usize::from(span.1 == ends.1),
        endpoint_distance: span
            .0
            .abs_diff(ends.0)
            .checked_add(span.1.abs_diff(ends.1))
            .ok_or_else(|| fail("locus_distance_overflow"))?,
    };
    let mut evidence = vec![
        format!("matched_junctions:{matched}"),
        format!("discriminating_exon_bp:{}", score.discriminating_exon_bp),
        "bam_direction_has_zero_assignment_weight".into(),
        "uncalibrated_literal_at_has_zero_strand_weight".into(),
    ];
    if score.structural_conflicts == 0
        && score.matched_junctions == 0
        && score.discriminating_exon_bp > 0
    {
        evidence.push("exon_containment_only_no_junction_support".into());
    }
    if score.endpoint_matches > 0 {
        evidence.push("endpoint_geometry_only_not_physical_rna_end".into());
    }
    let mut opposite_overlap = false;
    for m in models.iter().filter(|m| m.strand != model.strand) {
        if covered(&ranges, &m.exons, meter)? == bases {
            opposite_overlap = true;
        }
    }
    let class = if is_chimera {
        Class::Chimeric
    } else if !unknown.is_empty() && boundary_pair && outside == 0 {
        Class::NovelIsoformCompatible
    } else if score.structural_conflicts > 0 {
        Class::Unresolved
    } else if !score.informative() && opposite_overlap {
        Class::UninformativeOverlap
    } else if anchor.is_some_and(|a| a != model.strand) && score.informative() {
        Class::AntisenseCompatible
    } else if full_chain {
        Class::ExactTranscriptMatch
    } else {
        Class::CompatibleTranscripts
    };
    let map = exon_map(&model.exons, model.strand)?;
    let mut transcript_intervals = vec![];
    for range in &ranges {
        for exon in &map {
            meter.charge(1)?;
            if let Some(i) = geometry::project(exon, *range, model.strand) {
                transcript_intervals.push(i);
            }
        }
    }
    let transcript_intervals = meter.merge(transcript_intervals)?;
    let reason = if cds.is_some() {
        Reason::ModelProjectionIncompatible
    } else {
        Reason::CandidateModelUnavailable
    };
    let mut metrics = [missing(reason), missing(reason)];
    if score.structural_conflicts == 0
        && let Some(cds) = cds
    {
        let prepared = prepare(cds, meter)?;
        if let Some(projection) = one_path(path, &prepared, meter)? {
            metrics = projection.metrics;
        }
    }
    Ok(CandidatePath {
        hypothesis_index: index,
        class,
        score,
        evidence,
        conflicts,
        transcript_intervals,
        utr5: metrics[0].clone(),
        utr3: metrics[1].clone(),
        oriented,
    })
}

pub fn build(
    source: &str,
    region: &Region,
    inputs: &[Input<'_>],
    policy: oriented::Policy,
    limits: Limits,
    control: &mut dyn FnMut() -> Result<(), EngineError>,
) -> Result<Report, EngineError> {
    control()?;
    let window = region.window()?;
    if !token(source) || inputs.is_empty() || inputs.len() > 64 {
        return Err(fail("locus_admission"));
    }
    let limits = Limits {
        max_models: limits.max_models.min(100_000),
        max_intervals: limits.max_intervals.min(1_000_000),
        max_units: limits.max_units.min(100_000),
        max_observations: limits.max_observations.min(1_000_000),
        max_work: limits.max_work.min(100_000_000),
    };
    let mut meter = Meter {
        limits,
        work: 0,
        control: Some(control),
    };
    let (mut reads, mut all_models, mut samples) = (vec![], vec![], vec![]);
    let mut run_records = 0usize;
    let mut assignment_paths = 0usize;
    let mut retained_geometry = 0usize;
    let (mut sample_ids, mut source_units) = (BTreeSet::new(), BTreeSet::new());
    for input in inputs {
        if !token(input.sample_id) || !sample_ids.insert(input.sample_id) {
            return Err(fail("locus_sample_identity"));
        }
        if !sha(&input.bindings.report_content_sha256) || !sha(&input.bindings.data_sha256) {
            return Err(fail("locus_input_binding"));
        }
        if input.analysis.units.len() > limits.max_units
            || input.analysis.observations.len() > limits.max_observations
        {
            return Err(fail("locus_input_budget"));
        }
        let variants = oriented::prepare(input.original_bases, input.analysis, policy, &mut meter)?;
        let mut models = vec![];
        if let Some(annotation) = input.annotation {
            if annotation.transcripts.len() > limits.max_models {
                return Err(fail("locus_model_budget"));
            }
            if annotation.schema != crate::engine::PRIMER_PAYLOAD_ANNOTATION_SCHEMA
                || annotation.coordinate_convention != "zero_based_half_open"
                || annotation.assembly != input.analysis.reference.assembly
                || annotation.reference_fasta_sha256 != input.analysis.reference.fasta_sha256
                || !input
                    .bindings
                    .annotation_sha256
                    .as_ref()
                    .is_some_and(|s| sha(s))
            {
                return Err(fail("locus_annotation_binding"));
            }
            let mut ids = BTreeSet::new();
            for transcript in &annotation.transcripts {
                meter.charge(1)?;
                if !ids.insert(&transcript.transcript_id) {
                    return Err(fail("locus_model_identity"));
                }
                if transcript.genomic_contig.as_deref() != Some(region.contig.as_str()) {
                    continue;
                }
                intervals(&transcript.genomic_exons, true, &mut meter)?;
                if !transcript.genomic_exons.iter().any(|i| overlap(*i, window)) {
                    continue;
                }
                let strand = transcript
                    .genomic_strand
                    .filter(|s| matches!(s, '+' | '-'))
                    .ok_or_else(|| fail("locus_model_strand"))?;
                if transcript.genomic_exons.last().unwrap().1 > region.contig_length
                    || transcript.gene_ids.is_empty()
                    || transcript.target_name != region.contig
                    || !token(&transcript.transcript_id)
                    || transcript.gene_ids.iter().any(|g| !token(g))
                {
                    return Err(fail("locus_model_geometry"));
                }
                let hash = input.bindings.annotation_sha256.as_ref().unwrap();
                models.push(Model {
                    key: pileup::hash(&(input.sample_id, hash, &transcript.transcript_id))?,
                    sample_id: input.sample_id.into(),
                    transcript_id: transcript.transcript_id.clone(),
                    contig: region.contig.clone(),
                    gene_ids: transcript.gene_ids.clone(),
                    strand,
                    exons: transcript.genomic_exons.clone(),
                    annotation_release: annotation.release.clone(),
                    annotation_sha256: hash.clone(),
                    model_sha256: input.bindings.model_sha256.clone(),
                });
                if models.len() > limits.max_models {
                    return Err(fail("locus_model_budget"));
                }
            }
        }
        models.sort_by(|a, b| a.transcript_id.cmp(&b.transcript_id));
        let others = models
            .iter()
            .map(|m| other_exons(m, &models, &mut meter))
            .collect::<Result<Vec<_>, _>>()?;
        let mut cds = BTreeMap::new();
        if let Some(document) = input.model {
            if document.schema != MODEL_SCHEMA
                || document.coordinate_convention != "zero_based_half_open"
                || document.model_source.format != "gtf"
                || document.model_source.name.is_empty()
                || !sha(&document.model_source.sha256)
                || !input.bindings.model_sha256.as_ref().is_some_and(|h| sha(h))
                || document.source_id != source
                || document.assembly != input.analysis.reference.assembly
                || document.reference_fasta_sha256 != input.analysis.reference.fasta_sha256
                || input.bindings.annotation_sha256.as_deref()
                    != Some(document.annotation_sha256.as_str())
                || document.transcripts.len() > limits.max_models
            {
                return Err(fail("locus_cds_binding"));
            }
            for m in &document.transcripts {
                meter.charge(1)?;
                if cds.insert(m.transcript_id.as_str(), m).is_some() {
                    return Err(fail("locus_model_identity"));
                }
                if let Some(normal) = models.iter().find(|n| n.transcript_id == m.transcript_id) {
                    if normal.strand != m.strand
                        || normal.exons != m.genomic_exons_0based_half_open
                        || m.contig != region.contig
                        || !normal.gene_ids.contains(&m.gene_id)
                    {
                        return Err(fail("locus_model_annotation_geometry"));
                    }
                    prepare(m, &mut meter)?;
                }
            }
        }
        let anchor = if let Some(id) = &region.anchor_transcript_id {
            if let Some(m) = models.iter().find(|m| m.transcript_id == *id) {
                if !m.exons.contains(&region.exon) {
                    return Err(fail("locus_anchor_exon"));
                }
                Some(m.strand)
            } else {
                None
            }
        } else {
            None
        };
        let observations: BTreeMap<_, _> = input
            .analysis
            .observations
            .iter()
            .map(|o| (o.observation_id.as_str(), o))
            .collect();
        if observations.len() != input.analysis.observations.len() {
            return Err(fail("locus_observation_identity"));
        }
        let mut by_unit: BTreeMap<_, Vec<_>> = BTreeMap::new();
        let known_units: BTreeSet<_> = input.analysis.units.iter().map(|u| &u.unit).collect();
        for o in &input.analysis.observations {
            meter.charge(o.blocks.len() + o.skips.len() + 1)?;
            if !known_units.contains(&o.unit)
                || !matches!(o.query_to_reference_strand, '+' | '-')
                || o.skips.iter().any(|&(a, z)| {
                    a >= z || (o.contig == region.contig && z > region.contig_length)
                })
            {
                return Err(fail("locus_observation_binding"));
            }
            for b in &o.blocks {
                if b.query.0 >= b.query.1
                    || b.target.0 >= b.target.1
                    || b.query.1 > o.query_axis_length
                    || b.query.1 - b.query.0 != b.target.1 - b.target.0
                    || (o.contig == region.contig && b.target.1 > region.contig_length)
                {
                    return Err(fail("locus_observation_geometry"));
                }
            }
            by_unit.entry(&o.unit).or_default().push(o);
        }
        let mut summary = SampleSummary {
            sample_id: input.sample_id.into(),
            input_units: input.analysis.units.len(),
            diagnostics: input.diagnostics.clone(),
            ..Default::default()
        };
        if region.anchor_transcript_id.is_some() && anchor.is_none() {
            summary
                .diagnostics
                .push("anchor_model_unavailable_no_antisense_label".into());
        }
        for unit in &input.analysis.units {
            meter.charge(1)?;
            if source_units.len() >= limits.max_units {
                return Err(fail("locus_unit_budget"));
            }
            if !source_units.insert(unit.unit.clone()) {
                return Err(fail("locus_duplicate_source_unit"));
            }
            let observed: Vec<_> = by_unit.get(&unit.unit).cloned().unwrap_or_default();
            if !observed.iter().any(|o| {
                o.contig == region.contig && o.blocks.iter().any(|b| overlap(b.target, window))
            }) {
                summary.excluded.push((
                    unit.unit.clone(),
                    if unit.unavailable.is_empty() {
                        "no_saved_block_in_locus".into()
                    } else {
                        format!(
                            "upstream_unit_unavailable:{}",
                            unavailable_names(&unit.unavailable)?.join(",")
                        )
                    },
                ));
                continue;
            }
            if reads.len() >= limits.max_units {
                return Err(fail("locus_unit_budget"));
            }
            let mut hypotheses = vec![];
            let mut paths = vec![];
            for (index, ids) in unit.hypotheses.iter().enumerate() {
                if ids.is_empty() {
                    return Err(fail("locus_empty_hypothesis"));
                }
                let mut seen = BTreeSet::new();
                let mut path = vec![];
                for id in ids {
                    meter.charge(1)?;
                    let o = *observations
                        .get(id.as_str())
                        .ok_or_else(|| fail("locus_hypothesis_binding"))?;
                    if !seen.insert(id)
                        || o.unit != unit.unit
                        || !o.on_optimal_path
                        || unit.unclipped_axis_length != Some(o.query_axis_length)
                    {
                        return Err(fail("locus_hypothesis_binding"));
                    }
                    path.push(o);
                }
                hypotheses.push(Hypothesis {
                    index,
                    observation_ids: ids.clone(),
                    chimeric: chimera(&path),
                });
                paths.push(path);
            }
            let mut assignments = vec![];
            let mut genomic_basecalls = vec![];
            let original_variants: &[oriented::Variant] =
                variants.get(&unit.unit).map_or(&[], Vec::as_slice);
            let count = original_variants
                .iter()
                .try_fold(0usize, |n, v| n.checked_add(v.runs.len()))
                .ok_or_else(|| fail("locus_run_budget"))?;
            for path in &paths {
                run_records = run_records
                    .checked_add(count)
                    .filter(|n| *n <= limits.max_observations.min(oriented::MAX_RUN_RECORDS))
                    .ok_or_else(|| fail("locus_run_budget"))?;
                let compact: Vec<_> = path
                    .iter()
                    .filter(|o| o.contig == region.contig)
                    .map(|o| pileup::Observation {
                        observation_id: o.observation_id.clone(),
                        row_ordinal: o.source_ordinal,
                        blocks: o.blocks.clone(),
                        connectors: vec![],
                        query_to_reference_strand: o.query_to_reference_strand,
                    })
                    .collect();
                let mut e = oriented::project(
                    input.sample_id,
                    &unit.unit,
                    '+',
                    &compact,
                    original_variants,
                    input.original_bases,
                    policy,
                    &mut meter,
                )?;
                e.direction = oriented::Direction::Alternatives;
                for v in &mut e.variants {
                    for r in &mut v.runs {
                        r.transcript_ends = None;
                    }
                }
                genomic_basecalls.push(e);
            }
            if paths.is_empty() {
                run_records = run_records
                    .checked_add(count)
                    .filter(|n| *n <= limits.max_observations.min(oriented::MAX_RUN_RECORDS))
                    .ok_or_else(|| fail("locus_run_budget"))?;
                let mut e = oriented::project(
                    input.sample_id,
                    &unit.unit,
                    '+',
                    &[],
                    original_variants,
                    input.original_bases,
                    policy,
                    &mut meter,
                )?;
                e.direction = oriented::Direction::Alternatives;
                e.availability = "no_legal_hypothesis_reference_coverage_unavailable".into();
                for v in &mut e.variants {
                    for r in &mut v.runs {
                        r.transcript_ends = None;
                    }
                }
                genomic_basecalls.push(e);
            }
            let membership = if paths.iter().any(|p| {
                p.iter().any(|o| {
                    o.contig == region.contig && o.blocks.iter().any(|b| overlap(b.target, window))
                })
            }) {
                Membership::AdmittedHypothesis
            } else {
                summary.non_admitted_observation_only_units += 1;
                Membership::NonAdmittedObservationOnly
            };
            for (model_index, model) in models.iter().enumerate() {
                retained_geometry = retained_geometry
                    .checked_add(model.exons.len())
                    .filter(|n| *n <= limits.max_observations)
                    .ok_or_else(|| fail("locus_retained_geometry_budget"))?;
                let mut candidate_paths = vec![];
                for (index, path) in paths.iter().enumerate() {
                    meter.charge(1)?;
                    assignment_paths = assignment_paths
                        .checked_add(1)
                        .filter(|n| *n <= limits.max_observations.min(100_000))
                        .ok_or_else(|| fail("locus_assignment_budget"))?;
                    run_records = run_records
                        .checked_add(count)
                        .filter(|n| *n <= limits.max_observations.min(oriented::MAX_RUN_RECORDS))
                        .ok_or_else(|| fail("locus_run_budget"))?;
                    let path_result = candidate_path(
                        input.sample_id,
                        unit,
                        index,
                        path,
                        model,
                        &models,
                        &others[model_index],
                        cds.get(model.transcript_id.as_str()).copied(),
                        anchor,
                        original_variants,
                        input.original_bases,
                        policy,
                        &mut meter,
                    )?;
                    candidate_paths.push(path_result);
                }
                if candidate_paths.is_empty() {
                    continue;
                }
                let low = candidate_paths
                    .iter()
                    .min_by_key(|p| p.score.rank())
                    .unwrap();
                let high = candidate_paths
                    .iter()
                    .max_by_key(|p| p.score.rank())
                    .unwrap();
                let class = if candidate_paths.iter().all(|p| p.class == low.class) {
                    low.class
                } else {
                    Class::Unresolved
                };
                assignments.push(Assignment {
                    model: model.clone(),
                    class,
                    score: low.score.clone(),
                    score_upper: high.score.clone(),
                    preferred: false,
                    evidence: candidate_paths
                        .iter()
                        .flat_map(|p| p.evidence.iter().cloned())
                        .collect::<BTreeSet<_>>()
                        .into_iter()
                        .collect(),
                    conflicts: candidate_paths
                        .iter()
                        .flat_map(|p| p.conflicts.iter().cloned())
                        .collect::<BTreeSet<_>>()
                        .into_iter()
                        .collect(),
                    paths: candidate_paths,
                });
            }
            let best = assignments
                .iter()
                .filter(|a| a.score.structural_conflicts == 0)
                .map(|a| a.score.strong())
                .max();
            for a in &mut assignments {
                a.preferred = best.is_some_and(|best| {
                    a.score.structural_conflicts == 0 && a.score.strong() == best
                });
            }
            assignments.sort_by(|a, b| {
                b.score
                    .rank()
                    .cmp(&a.score.rank())
                    .then_with(|| a.model.transcript_id.cmp(&b.model.transcript_id))
            });
            let preferred: Vec<_> = assignments.iter().filter(|a| a.preferred).collect();
            let strands: BTreeSet<_> = preferred.iter().map(|a| a.model.strand).collect();
            let supported = !preferred.is_empty()
                && preferred.iter().all(|a| a.score.informative())
                && preferred.iter().all(|a| a.score.matched_junctions > 0)
                && strands.len() == 1;
            let mut molecule_strand = MoleculeStrand {
                status: if supported {
                    "model_supported"
                } else {
                    "unresolved"
                }
                .into(),
                strand: if supported {
                    strands.first().copied()
                } else {
                    None
                },
                reasons: vec![
                    if supported {
                        "matched_junction_model_support_not_bam_direction"
                    } else {
                        "double_stranded_cdna_alignment_does_not_establish_rna_origin"
                    }
                    .into(),
                ],
            };
            if !preferred.is_empty()
                && strands.len() == 1
                && preferred
                    .iter()
                    .all(|a| a.score.structural_conflicts == 0 && a.score.matched_junctions == 0)
            {
                molecule_strand
                    .reasons
                    .push("exon_containment_only_no_junction_support".into());
            }
            if !unit.search_declared_complete
                || !unit.views_declared_complete
                || !unit.all_query_views_assessed
            {
                molecule_strand
                    .reasons
                    .push("mapping_search_or_views_incomplete_other_loci_not_excluded".into());
            }
            let alignment_directions = observed
                .iter()
                .map(|o| o.query_to_reference_strand)
                .collect::<BTreeSet<_>>()
                .into_iter()
                .collect();
            reads.push(Read {
                sample_id: input.sample_id.into(),
                read_key: oriented::read_key(input.sample_id, &unit.unit)?,
                unit: unit.unit.clone(),
                membership,
                query_length: unit.unclipped_axis_length,
                bindings: input.bindings.clone(),
                observations: observed
                    .iter()
                    .map(|o| Observation {
                        observation_id: o.observation_id.clone(),
                        on_optimal_path: o.on_optimal_path,
                        source_ordinal: o.source_ordinal,
                        contig: o.contig.clone(),
                        query_to_reference_strand: o.query_to_reference_strand,
                        blocks: o.blocks.clone(),
                        skips: o.skips.clone(),
                        clipping: o.clipping.clone(),
                    })
                    .collect(),
                rows: unit.rows.clone(),
                alignment_directions,
                supplementary_present: unit.rows.iter().any(|r| r.flag & 0x800 != 0),
                search_declared_complete: unit.search_declared_complete,
                views_declared_complete: unit.views_declared_complete,
                all_query_views_assessed: unit.all_query_views_assessed,
                hypotheses,
                genomic_basecalls,
                assignments,
                molecule_strand,
                unavailable: unavailable_names(&unit.unavailable)?,
            });
            summary.locus_units += 1;
        }
        if all_models.len() + models.len() > limits.max_models {
            return Err(fail("locus_model_budget"));
        }
        all_models.extend(models);
        samples.push(summary);
    }
    reads.sort_by(|a, b| (&a.sample_id, &a.unit).cmp(&(&b.sample_id, &b.unit)));
    let scoring_rules = [
        "lexicographic: fewer structural conflicts, matched junctions, full chain, discriminating exon bases, endpoint geometry, proximity".into(),
        "conditional model_supported strand requires matched junctions on every preferred assignment; exon containment alone leaves strand unresolved".into(),
        "full_chain means a conflict-free complete annotated junction chain with at least one intron; deletions and endpoints are separate".into(),
        "preferred flags compare the structural prefix; weak end proximity cannot break structural ties, and containment-only preference does not resolve RNA strand".into(),
        "A/T is uncalibrated and has zero strand-assignment weight; no canonical motif is inferred".into(),
    ];
    let scoring_rules = Vec::from(scoring_rules);
    let nonclaims = [
        "Scores are descriptive geometry ranks, not probabilities or proof of whole RNA, APA, physical tail length or causality.".into(),
        "Annotation strand, alignment direction and conditional molecule-strand support are distinct; pore direction is unknown.".into(),
        "Masked query bases are not removed from locus base counts.".into(),
        "A supplementary alignment overlapping the primary on the query never joins a hypothesis; classification uses the admitted path alone.".into(),
    ];
    let nonclaims = Vec::from(nonclaims);
    Ok(Report {
        schema: SCHEMA.into(),
        identity_scope: "sam_declared_query_unit".into(),
        source_id: source.into(),
        region: region.clone(),
        window,
        scoring_version: RULE_VERSION.into(),
        effective_limits: limits,
        work: meter.work,
        scoring_rules,
        models: all_models,
        reads,
        samples,
        nonclaims,
    })
}

fn unavailable_names(reasons: &[ObservedSamUnavailableReason]) -> Result<Vec<String>, EngineError> {
    reasons
        .iter()
        .map(|r| {
            serde_json::to_value(r)
                .map_err(|e| fail(&e.to_string()))?
                .as_str()
                .map(str::to_owned)
                .ok_or_else(|| fail("locus_unavailable_reason"))
        })
        .collect()
}
