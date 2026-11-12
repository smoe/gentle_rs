//! C2-C4: conditional consequences of saved cDNA observations, not a new allele caller.
//! The service authenticates files/models and injects the existing genetic code.

use super::cdna_comparison as evidence;
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};

mod hgvs;
pub use hgvs::CodonEffect;
mod linkage;
#[cfg(test)]
mod tests;

pub const SCHEMA: &str = "gentle.cdna_allele_interpretation.v1";
pub const HGVS_VERSION: &str = "21.1.5 (2026-09-17)";

/// Portable C1 projection. Sequence is mature transcript 5'-to-3', not genomic forward.
#[derive(Debug, Clone)]
pub struct CodingModel {
    pub transcript: evidence::CdnaTranscript,
    pub accession: Option<String>,
    pub protein_accession: Option<String>,
    pub sequence: String,
    pub cds: Option<(usize, usize)>,
    pub translation_status: String,
    pub translation_reasons: Vec<String>,
}

/// Full supported alignment geometry, including bases not used for variant calling.
#[derive(Debug, Clone, Default)]
pub struct AlignmentGeometry {
    pub blocks: Vec<(String, usize, usize)>,
    pub junctions: Vec<(String, usize, usize)>,
}

#[derive(Debug, Clone)]
pub struct LinkageEvidence {
    pub samples: Vec<evidence::CdnaSample>,
    pub geometry: BTreeMap<(String, String, String), AlignmentGeometry>,
}

#[derive(Debug, Clone)]
pub struct GeneInput {
    pub source_id: String,
    pub gene_id: String,
    pub report: evidence::CdnaComparisonReport,
    pub models: Vec<CodingModel>,
    pub genomic_sequence: String,
    pub linkage: Option<LinkageEvidence>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct Policy {
    pub translation_table: usize,
    /// None inherits the unchanged bound comparison threshold; conflicting overrides refuse.
    pub minimum_codon_haplotype_units: Option<usize>,
    pub maximum_out_of_model_aligned_bases: usize,
}
impl Default for Policy {
    fn default() -> Self {
        Self {
            translation_table: 1,
            minimum_codon_haplotype_units: None,
            maximum_out_of_model_aligned_bases: 0,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct Budgets {
    pub max_input_bytes: u64,
    pub max_work: u64,
    pub max_output_bytes: usize,
}
impl Default for Budgets {
    fn default() -> Self {
        Self {
            max_input_bytes: 128_000_000,
            max_work: 100_000_000,
            max_output_bytes: 32_000_000,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct KnownVariant {
    pub contig: String,
    pub position_1based: usize,
    pub r#ref: String,
    pub alt: String,
    pub identifier: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct KnownVariantManifest {
    pub schema: String,
    pub source_name: String,
    pub source_release: String,
    pub retrieved_utc: String,
    pub license: String,
    pub attribution: String,
    pub assembly: String,
    pub entries: Vec<KnownVariant>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct KnownLabel {
    pub identifier: String,
    pub source_name: String,
    pub source_release: String,
    pub license: String,
    pub attribution: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Change {
    pub genomic_position_0based: usize,
    pub cdna_position_0based: usize,
    pub reference_base: char,
    pub alternate_base: char,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Consequence {
    pub assessable: bool,
    pub region: String,
    pub cdna_position_1based: Option<usize>,
    pub nucleotide_hgvs: Option<String>,
    pub protein_hgvs: Option<String>,
    pub codon_index_1based: Option<usize>,
    pub reference_codon: Option<String>,
    pub alternate_codon: Option<String>,
    pub reference_amino_acid: Option<String>,
    pub alternate_amino_acid: Option<String>,
    pub consequence_terms: Vec<String>,
    pub reasons: Vec<String>,
    pub coding_changes: Vec<CodonEffect>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SiteRow {
    pub sample_id: String,
    pub transcript_id: String,
    pub genomic_position_0based: usize,
    pub genomic_reference_base: char,
    pub genomic_alternate_base: char,
    pub alternate_units: usize,
    pub callable_units: usize,
    pub alternate_unit_fraction: Option<f64>,
    pub inherited_counts: evidence::CdnaSiteSampleCounts,
    pub inherited_unit_statuses: BTreeMap<String, usize>,
    pub support_status: String,
    pub codon_context: String,
    pub consequence: Consequence,
    pub known_variants: Vec<KnownLabel>,
    pub linkage: Option<IsoformLinkage>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct IsoformLinkage {
    pub status: String,
    pub compatible_units: Vec<String>,
    pub discriminating_units: Vec<String>,
    pub unique_model_units: Vec<String>,
    pub excluded_transcripts: BTreeMap<String, Vec<String>>,
    pub reasons: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CodonCombination {
    pub sample_id: String,
    pub transcript_id: String,
    pub positions_0based: Vec<usize>,
    pub genomic_alleles: Vec<char>,
    pub transcript_changes: Vec<Change>,
    pub unit_ids: Vec<String>,
    pub incomplete_units: Vec<String>,
    pub minimum_support: usize,
    pub supported: bool,
    pub consequence: Consequence,
    pub known_variants: Vec<KnownLabel>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CodonGroup {
    pub sample_id: String,
    pub transcript_id: String,
    pub positions_0based: Vec<usize>,
    pub complete_units: usize,
    pub incomplete_units: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LocalPattern {
    pub fingerprint: String,
    pub phase_method: String,
    pub alleles: Vec<(usize, char)>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PhaseRelation {
    pub sample_id: String,
    pub positions_0based: (usize, usize),
    pub relation: String,
    pub block_ids: Vec<String>,
    pub local_patterns: Vec<LocalPattern>,
    pub scope: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GeneReport {
    pub source_id: String,
    pub gene_id: String,
    pub gene: String,
    pub comparison_content_sha256: String,
    pub comparison_request_sha256: String,
    pub reference: evidence::CdnaReferenceSummary,
    pub inherited_policy: evidence::CdnaComparisonPolicy,
    pub inherited_samples: Vec<evidence::CdnaSampleReport>,
    pub inherited_warnings: Vec<String>,
    pub minimum_codon_haplotype_units: usize,
    pub sites: Vec<SiteRow>,
    pub codon_combinations: Vec<CodonCombination>,
    pub codon_groups: Vec<CodonGroup>,
    pub local_phase_relations: Vec<PhaseRelation>,
}

pub(super) struct Meter<'a> {
    pub work: u64,
    pub limit: u64,
    pub progress: &'a mut dyn FnMut() -> bool,
    pub remaining_output: usize,
}
impl Meter<'_> {
    fn output<T: Serialize>(&mut self, value: &T) -> Result<(), String> {
        let bytes = serde_json::to_vec(value).map_err(|e| e.to_string())?.len();
        self.remaining_output = self
            .remaining_output
            .checked_sub(bytes)
            .ok_or_else(|| fail("output_budget"))?;
        Ok(())
    }
    pub fn charge(&mut self, n: usize) -> Result<(), String> {
        if !(self.progress)() {
            return Err("cdna_interpretation.cancelled".into());
        }
        self.work = self
            .work
            .checked_add(n as u64)
            .filter(|n| *n <= self.limit)
            .ok_or("cdna_interpretation.work_budget")?;
        Ok(())
    }
}

fn fail(reason: &str) -> String {
    format!("cdna_interpretation.{reason}")
}
fn complement(base: char, strand: char) -> char {
    if strand == '-' {
        crate::iupac_code::IupacCode::letter_complement(base as u8) as char
    } else {
        base
    }
}
fn position(model: &CodingModel, genomic: usize) -> Option<usize> {
    let mut offset = 0;
    for &(a, z) in &model.transcript.exons {
        if a <= genomic && genomic < z {
            let plus = offset + genomic - a;
            return Some(if model.transcript.strand == '-' {
                model.sequence.len() - 1 - plus
            } else {
                plus
            });
        }
        offset += z - a;
    }
    None
}

/// Consequence kernel over authenticated immutable reports; does not rerun the allele caller.
pub fn interpret(
    genes: &[GeneInput],
    policy: &Policy,
    code: &BTreeMap<String, char>,
    known: Option<&KnownVariantManifest>,
    budgets: &Budgets,
    progress: &mut dyn FnMut() -> bool,
) -> Result<(Vec<GeneReport>, u64), String> {
    if genes.is_empty()
        || genes.len() > 64
        || policy.translation_table != 1
        || code.len() != 64
        || budgets.max_work == 0
        || budgets.max_work > 100_000_000
        || budgets.max_output_bytes == 0
        || budgets.max_output_bytes > 64_000_000
        || policy.maximum_out_of_model_aligned_bases > 1_000_000
    {
        return Err(fail("schema_policy_or_budget"));
    }
    for a in b"ACGT" {
        for b in b"ACGT" {
            for c in b"ACGT" {
                let codon = String::from_utf8(vec![*a, *b, *c]).unwrap();
                if code.get(&codon).is_none_or(|a| hgvs::amino(*a).is_none()) {
                    return Err(fail("genetic_code"));
                }
            }
        }
    }
    let mut meter = Meter {
        work: 0,
        limit: budgets.max_work,
        progress,
        remaining_output: budgets.max_output_bytes,
    };
    validate_known(known, genes, &mut meter)?;
    let mut seen = BTreeSet::new();
    let mut output = Vec::new();
    let baseline = serde_json::to_value(&genes[0].report.policy).map_err(|e| e.to_string())?;
    for input in genes {
        meter.charge(1)?;
        let r = &input.report;
        if !seen.insert((&input.source_id, &input.gene_id))
            || serde_json::to_value(&r.policy).map_err(|e| e.to_string())? != baseline
        {
            return Err(fail("duplicate_gene_or_policy_mismatch"));
        }
        let threshold = r.policy.minimum_haplotype_units;
        if policy
            .minimum_codon_haplotype_units
            .is_some_and(|n| n != threshold)
        {
            return Err(fail("support_policy_mismatch"));
        }
        validate_models(input, &mut meter)?;
        meter.output(&r.samples)?;
        meter.output(&r.warnings)?;
        let mut report = GeneReport {
            source_id: input.source_id.clone(),
            gene_id: input.gene_id.clone(),
            gene: r.gene.clone(),
            comparison_content_sha256: r.content_sha256.clone(),
            comparison_request_sha256: r.request_sha256.clone(),
            reference: r.reference.clone(),
            inherited_policy: r.policy.clone(),
            inherited_samples: r.samples.clone(),
            inherited_warnings: r.warnings.clone(),
            minimum_codon_haplotype_units: threshold,
            sites: Vec::new(),
            codon_combinations: Vec::new(),
            codon_groups: Vec::new(),
            local_phase_relations: Vec::new(),
        };
        for sample in &r.samples {
            meter.charge(sample.site_calls.len())?;
            let calls = sample
                .site_calls
                .iter()
                .map(|c| ((c.unit_id.as_str(), c.position_0based), c))
                .collect::<BTreeMap<_, _>>();
            if calls.len() != sample.site_calls.len() {
                return Err(fail("duplicate_unit_site"));
            }
            for model in &input.models {
                let mut mapped = Vec::new();
                let mut complete_allele_context = BTreeSet::new();
                for site in &r.sites {
                    meter.charge(model.transcript.exons.len() + 1)?;
                    if !alternates(site).is_empty()
                        && let Some(p) = position(model, site.position_0based)
                    {
                        mapped.push((p, site.position_0based));
                    }
                }
                mapped.sort_unstable();
                let mut groups = Vec::<Vec<usize>>::new();
                for (index, &(p, genomic)) in mapped.iter().enumerate() {
                    let connected = index > 0
                        && (mapped[index - 1].0 + 1 == p
                            || model.cds.is_some_and(|(a, z)| {
                                mapped[index - 1].0 >= a
                                    && p < z
                                    && (mapped[index - 1].0 - a) / 3 == (p - a) / 3
                            }));
                    if connected {
                        groups.last_mut().unwrap().push(genomic);
                    } else {
                        groups.push(vec![genomic]);
                    }
                }
                for positions in groups.iter().filter(|v| v.len() > 1) {
                    meter.charge(sample.site_calls.len())?;
                    let mut combinations = BTreeMap::<Vec<char>, Vec<String>>::new();
                    let mut incomplete = Vec::new();
                    let units = sample
                        .site_calls
                        .iter()
                        .map(|c| c.unit_id.as_str())
                        .collect::<BTreeSet<_>>();
                    for unit in units {
                        meter.charge(positions.len())?;
                        let alleles = positions
                            .iter()
                            .map(|p| {
                                calls
                                    .get(&(unit, *p))
                                    .filter(|c| definite(c))
                                    .and_then(|c| c.base)
                            })
                            .collect::<Option<Vec<_>>>();
                        if let Some(alleles) = alleles {
                            combinations.entry(alleles).or_default().push(unit.into());
                        } else {
                            incomplete.push(unit.into());
                        }
                    }
                    let group = CodonGroup {
                        sample_id: sample.sample_id.clone(),
                        transcript_id: model.transcript.transcript_id.clone(),
                        positions_0based: positions.clone(),
                        complete_units: combinations.values().map(Vec::len).sum(),
                        incomplete_units: incomplete.clone(),
                    };
                    meter.output(&group)?;
                    report.codon_groups.push(group);
                    for (alleles, units) in combinations {
                        complete_allele_context
                            .extend(positions.iter().copied().zip(alleles.iter().copied()));
                        meter
                            .charge(positions.len().saturating_mul(model.transcript.exons.len()))?;
                        let changes = positions
                            .iter()
                            .zip(&alleles)
                            .map(|(p, alt)| change(input, model, *p, *alt))
                            .collect::<Result<Vec<_>, _>>()?;
                        let consequence = hgvs::consequence(model, &changes, code, &mut meter)?;
                        let supported = units.len() >= threshold;
                        let combination = CodonCombination {
                            sample_id: sample.sample_id.clone(),
                            transcript_id: model.transcript.transcript_id.clone(),
                            positions_0based: positions.clone(),
                            genomic_alleles: alleles.clone(),
                            transcript_changes: changes,
                            unit_ids: units,
                            incomplete_units: incomplete.clone(),
                            minimum_support: threshold,
                            supported,
                            consequence,
                            known_variants: if supported {
                                known_labels(known, input, positions, &alleles, &mut meter)?
                            } else {
                                Vec::new()
                            },
                        };
                        meter.output(&combination)?;
                        report.codon_combinations.push(combination);
                    }
                }
                for site in &r.sites {
                    let counts = site
                        .samples
                        .iter()
                        .find(|c| c.sample_id == sample.sample_id)
                        .ok_or_else(|| fail("sample_counts"))?;
                    for alt in alternates(site) {
                        meter.charge(model.transcript.exons.len() + sample.site_calls.len() + 1)?;
                        let count = counts.base_counts.get(&alt).copied().unwrap_or(0);
                        let change = position(model, site.position_0based)
                            .map(|_| change(input, model, site.position_0based, alt))
                            .transpose()?;
                        let consequence = hgvs::consequence(
                            model,
                            &change.into_iter().collect::<Vec<_>>(),
                            code,
                            &mut meter,
                        )?;
                        let mut statuses = BTreeMap::new();
                        for call in sample
                            .site_calls
                            .iter()
                            .filter(|c| c.position_0based == site.position_0based)
                        {
                            let status =
                                serde_json::to_value(call.status).map_err(|e| e.to_string())?;
                            *statuses
                                .entry(status.as_str().unwrap().to_owned())
                                .or_default() += 1;
                        }
                        let context = if groups
                            .iter()
                            .any(|g| g.len() > 1 && g.contains(&site.position_0based))
                        {
                            if complete_allele_context.contains(&(site.position_0based, alt)) {
                                "see_codon_haplotype"
                            } else {
                                "codon_context_unresolved"
                            }
                        } else {
                            "reference_assumed_no_other_candidate"
                        };
                        let row = SiteRow {
                            sample_id: sample.sample_id.clone(),
                            transcript_id: model.transcript.transcript_id.clone(),
                            genomic_position_0based: site.position_0based,
                            genomic_reference_base: site.reference_base,
                            genomic_alternate_base: alt,
                            alternate_units: count,
                            callable_units: counts.callable_units,
                            alternate_unit_fraction: (counts.callable_units > 0)
                                .then(|| count as f64 / counts.callable_units as f64),
                            inherited_counts: counts.clone(),
                            inherited_unit_statuses: statuses,
                            support_status: if counts.callable_units == 0 {
                                "unassessed"
                            } else if counts.candidate_alternate_bases.contains(&alt) {
                                "inherited_candidate"
                            } else {
                                "below_inherited_threshold"
                            }
                            .into(),
                            codon_context: context.into(),
                            consequence,
                            known_variants: known_labels(
                                known,
                                input,
                                &[site.position_0based],
                                &[alt],
                                &mut meter,
                            )?,
                            linkage: input
                                .linkage
                                .as_ref()
                                .map(|link| {
                                    linkage::link(
                                        input,
                                        link,
                                        sample,
                                        model,
                                        site.position_0based,
                                        alt,
                                        policy,
                                        &mut meter,
                                    )
                                })
                                .transpose()?,
                        };
                        meter.output(&row)?;
                        report.sites.push(row);
                    }
                }
            }
            phase_relations(
                sample,
                &r.sites,
                &mut report.local_phase_relations,
                &mut meter,
            )?;
        }
        output.push(report);
    }
    Ok((output, meter.work))
}

fn definite(call: &evidence::CdnaUnitAlleleCall) -> bool {
    matches!(
        call.status,
        evidence::CdnaCallStatus::ReferenceSupported | evidence::CdnaCallStatus::AlternateSupported
    )
}

fn alternates(site: &evidence::CdnaSiteReport) -> BTreeSet<char> {
    site.candidate_alternate_bases
        .iter()
        .copied()
        .chain(site.samples.iter().flat_map(|s| {
            s.base_counts
                .iter()
                .filter(|(_, n)| **n > 0)
                .map(|(base, _)| *base)
        }))
        .filter(|base| *base != site.reference_base)
        .collect()
}

fn change(input: &GeneInput, model: &CodingModel, p: usize, alt: char) -> Result<Change, String> {
    let local = p
        .checked_sub(input.report.reference.start_0based)
        .ok_or_else(|| fail("reference_position"))?;
    let reference = *input
        .genomic_sequence
        .as_bytes()
        .get(local)
        .ok_or_else(|| fail("reference_position"))? as char;
    Ok(Change {
        genomic_position_0based: p,
        cdna_position_0based: position(model, p).ok_or_else(|| fail("model_position"))?,
        reference_base: complement(reference, model.transcript.strand),
        alternate_base: complement(alt, model.transcript.strand),
    })
}

fn validate_models(input: &GeneInput, meter: &mut Meter<'_>) -> Result<(), String> {
    let r = &input.report;
    if input.models.len() != r.transcripts.len()
        || input.models.is_empty()
        || input.models.len() > 256
        || input.genomic_sequence.len()
            != r.reference
                .end_0based_exclusive
                .saturating_sub(r.reference.start_0based)
    {
        return Err(fail("model_geometry_mismatch"));
    }
    let mut positions = BTreeSet::new();
    let sample_ids = r
        .samples
        .iter()
        .map(|s| s.sample_id.as_str())
        .collect::<BTreeSet<_>>();
    if sample_ids.len() != r.samples.len() {
        return Err(fail("duplicate_sample"));
    }
    for site in &r.sites {
        meter.charge(site.samples.len() + 1)?;
        let local = site
            .position_0based
            .checked_sub(r.reference.start_0based)
            .ok_or_else(|| fail("site_reference_binding"))?;
        if !positions.insert(site.position_0based)
            || input
                .genomic_sequence
                .as_bytes()
                .get(local)
                .copied()
                .map(char::from)
                != Some(site.reference_base)
            || alternates(site)
                .iter()
                .any(|a| !matches!(a, 'A' | 'C' | 'G' | 'T'))
        {
            return Err(fail("site_reference_binding"));
        }
        let ids = site
            .samples
            .iter()
            .map(|s| s.sample_id.as_str())
            .collect::<BTreeSet<_>>();
        if ids != sample_ids || ids.len() != site.samples.len() {
            return Err(fail("sample_counts"));
        }
    }
    let mut ids = BTreeSet::new();
    for model in &input.models {
        meter.charge(model.sequence.len())?;
        let t = &model.transcript;
        let old = r
            .transcripts
            .iter()
            .find(|x| x.transcript_id == t.transcript_id)
            .ok_or_else(|| fail("model_geometry_mismatch"))?;
        if !ids.insert(&t.transcript_id)
            || t.strand != old.strand
            || t.exons != old.exons
            || !matches!(t.strand, '+' | '-')
            || t.exons.iter().any(|&(a, z)| a >= z)
            || t.exons.windows(2).any(|w| w[0].1 > w[1].0)
            || t.exons
                .iter()
                .try_fold(0usize, |n, &(a, z)| n.checked_add(z.saturating_sub(a)))
                != Some(model.sequence.len())
            || model
                .cds
                .is_some_and(|(a, z)| a >= z || z > model.sequence.len())
        {
            return Err(fail("model_geometry_mismatch"));
        }
        let mut offset = 0usize;
        for &(a, z) in &t.exons {
            for p in a..z {
                let local = p
                    .checked_sub(r.reference.start_0based)
                    .ok_or_else(|| fail("model_geometry_mismatch"))?;
                let base = *input
                    .genomic_sequence
                    .as_bytes()
                    .get(local)
                    .ok_or_else(|| fail("model_geometry_mismatch"))?
                    as char;
                let plus = offset + p - a;
                let q = if t.strand == '-' {
                    model.sequence.len() - 1 - plus
                } else {
                    plus
                };
                if model.sequence.as_bytes()[q] as char != complement(base, t.strand) {
                    return Err(fail("model_sequence_mismatch"));
                }
            }
            offset += z - a;
        }
        for site in &r.sites {
            meter.charge(t.exons.len() + 1)?;
            let saved = site
                .transcripts
                .iter()
                .find(|a| a.transcript_id == t.transcript_id)
                .ok_or_else(|| fail("transcript_coordinate_drift"))?;
            if saved.cdna_position_1based != position(model, site.position_0based).map(|p| p + 1) {
                return Err(fail("transcript_coordinate_drift"));
            }
        }
    }
    Ok(())
}

fn phase_relations(
    sample: &evidence::CdnaSampleReport,
    sites: &[evidence::CdnaSiteReport],
    out: &mut Vec<PhaseRelation>,
    meter: &mut Meter<'_>,
) -> Result<(), String> {
    for (i, a) in sites.iter().enumerate() {
        for b in sites.iter().skip(i + 1) {
            meter.charge(1)?;
            let left = sample
                .phase_blocks
                .iter()
                .find(|v| v.sites_0based.contains(&a.position_0based));
            let right = sample
                .phase_blocks
                .iter()
                .find(|v| v.sites_0based.contains(&b.position_0based));
            meter.charge(
                sample
                    .phase_blocks
                    .iter()
                    .map(|b| b.sites_0based.len())
                    .sum(),
            )?;
            let (relation, blocks, mut patterns) = match (left, right) {
                (Some(l), Some(r)) if l.block_id == r.block_id => (
                    if l.haplotypes.is_empty() {
                        "same_block_unresolved"
                    } else {
                        "same_block_candidate_patterns"
                    },
                    vec![l.block_id.clone()],
                    Vec::new(),
                ),
                (Some(l), Some(r)) => (
                    "different_blocks",
                    vec![l.block_id.clone(), r.block_id.clone()],
                    Vec::new(),
                ),
                _ => ("not_in_any_block", Vec::new(), Vec::new()),
            };
            if let (Some(l), Some(r)) = (left, right)
                && l.block_id == r.block_id
            {
                for h in &l.haplotypes {
                    meter.charge(h.alleles.len() + 1)?;
                    patterns.push(LocalPattern {
                        fingerprint: h.fingerprint.clone(),
                        phase_method: h.phase_method.clone(),
                        alleles: h
                            .alleles
                            .iter()
                            .filter(|allele| {
                                allele.position_0based == a.position_0based
                                    || allele.position_0based == b.position_0based
                            })
                            .map(|a| (a.position_0based, a.observed_base))
                            .collect(),
                    });
                }
            }
            let relation = PhaseRelation {
                sample_id: sample.sample_id.clone(),
                positions_0based: (a.position_0based, b.position_0based),
                relation: relation.into(),
                block_ids: blocks,
                local_patterns: patterns,
                scope: "local_RNA_candidate_patterns_not_chromosomal_phase".into(),
            };
            meter.output(&relation)?;
            out.push(relation);
        }
    }
    Ok(())
}

fn validate_known(
    known: Option<&KnownVariantManifest>,
    genes: &[GeneInput],
    meter: &mut Meter<'_>,
) -> Result<(), String> {
    let Some(k) = known else {
        return Ok(());
    };
    if k.schema != "gentle.known_variant_manifest.v1"
        || k.entries.len() > 100_000
        || [
            &k.source_name,
            &k.source_release,
            &k.license,
            &k.attribution,
        ]
        .iter()
        .any(|s| s.is_empty() || s.len() > 4096 || s.chars().any(char::is_control))
        || !valid_utc(&k.retrieved_utc)
        || genes
            .iter()
            .any(|g| g.report.reference.assembly != k.assembly)
    {
        return Err(fail("known_manifest"));
    }
    let mut seen = BTreeSet::new();
    for e in &k.entries {
        meter.charge(e.r#ref.len() + e.alt.len() + 1)?;
        if e.position_1based == 0
            || e.contig.is_empty()
            || e.contig.len() > 512
            || e.contig.chars().any(char::is_control)
            || e.position_1based.checked_add(e.r#ref.len()).is_none()
            || e.identifier.is_empty()
            || e.identifier.len() > 512
            || e.identifier.chars().any(char::is_control)
            || e.r#ref.is_empty()
            || e.alt.is_empty()
            || e.r#ref == e.alt
            || !e
                .r#ref
                .bytes()
                .chain(e.alt.bytes())
                .all(|b| b"ACGT".contains(&b))
            || !seen.insert((
                &e.contig,
                e.position_1based,
                &e.r#ref,
                &e.alt,
                &e.identifier,
            ))
        {
            return Err(fail("known_manifest_entry"));
        }
    }
    Ok(())
}

fn valid_utc(value: &str) -> bool {
    let b = value.as_bytes();
    if b.len() != 20
        || b[4] != b'-'
        || b[7] != b'-'
        || b[10] != b'T'
        || b[13] != b':'
        || b[16] != b':'
        || b[19] != b'Z'
    {
        return false;
    }
    let number = |a, z| {
        std::str::from_utf8(&b[a..z])
            .ok()
            .filter(|v| v.bytes().all(|c| c.is_ascii_digit()))
            .and_then(|v| v.parse::<usize>().ok())
    };
    let (Some(year), Some(month), Some(day), Some(hour), Some(minute), Some(second)) = (
        number(0, 4),
        number(5, 7),
        number(8, 10),
        number(11, 13),
        number(14, 16),
        number(17, 19),
    ) else {
        return false;
    };
    let leap = year % 4 == 0 && (year % 100 != 0 || year % 400 == 0);
    let days = match month {
        1 | 3 | 5 | 7 | 8 | 10 | 12 => 31,
        4 | 6 | 9 | 11 => 30,
        2 => {
            if leap {
                29
            } else {
                28
            }
        }
        _ => return false,
    };
    year > 0 && day > 0 && day <= days && hour < 24 && minute < 60 && second < 60
}

fn known_labels(
    known: Option<&KnownVariantManifest>,
    input: &GeneInput,
    positions: &[usize],
    alleles: &[char],
    meter: &mut Meter<'_>,
) -> Result<Vec<KnownLabel>, String> {
    let Some(k) = known else {
        return Ok(Vec::new());
    };
    let changed = positions
        .iter()
        .zip(alleles)
        .filter(|(p, a)| {
            input.genomic_sequence.as_bytes()[**p - input.report.reference.start_0based] as char
                != **a
        })
        .collect::<Vec<_>>();
    let Some(start) = changed.iter().map(|(p, _)| **p).min() else {
        return Ok(Vec::new());
    };
    let end = changed.iter().map(|(p, _)| **p).max().unwrap() + 1;
    meter.charge(end - start)?;
    let reference = &input.genomic_sequence
        [start - input.report.reference.start_0based..end - input.report.reference.start_0based];
    let mut alt = reference.as_bytes().to_vec();
    for (p, a) in changed {
        alt[*p - start] = *a as u8;
    }
    let mut labels = Vec::new();
    let mut label_bytes = 0usize;
    for e in &k.entries {
        meter.charge(1)?;
        if e.contig == input.report.reference.contig
            && e.position_1based == start + 1
            && e.r#ref == reference
            && e.alt.as_bytes() == alt.as_slice()
        {
            let label = KnownLabel {
                identifier: e.identifier.clone(),
                source_name: k.source_name.clone(),
                source_release: k.source_release.clone(),
                license: k.license.clone(),
                attribution: k.attribution.clone(),
            };
            label_bytes = label_bytes
                .checked_add(serde_json::to_vec(&label).map_err(|e| e.to_string())?.len())
                .ok_or_else(|| fail("output_budget"))?;
            if label_bytes > meter.remaining_output {
                return Err(fail("output_budget"));
            }
            labels.push(label);
        }
    }
    labels.sort_by(|a, b| a.identifier.cmp(&b.identifier));
    Ok(labels)
}
