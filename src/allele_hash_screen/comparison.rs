//! Explicit quality-bearing cDNA evidence ingress over saved SAM hypotheses, without another aligner.

use crate::engine::{
    EngineError, OBSERVED_SAM_DATA_SCHEMA, ObservedSamConversion, ObservedSamData,
    ObservedSamQueryUnitKey, ObservedSamResultBinding, OperationProgress,
    PRIMER_PAYLOAD_ANNOTATION_SCHEMA, PrimerPayloadAnnotation, PrimerPayloadFile,
};
use crate::observed_transcript_structures::geometry;
use gentle_engine::allele_hash_screen::cdna_comparison as core;
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};

pub const CDNA_FILE_COMPARISON_SCHEMA: &str = "gentle.cdna_allele_comparison.v1";
pub const CDNA_FILE_COMPARISON_REQUEST_SCHEMA: &str = "gentle.cdna_allele_comparison_request.v1";
pub const CDNA_BASECALL_SCHEMA: &str = "gentle.cdna_allele_basecalls.v1";

fn default_input_bytes() -> u64 {
    64_000_000
}
fn default_basecall_bases() -> usize {
    2_000_000
}
fn fail(message: impl ToString) -> EngineError {
    EngineError::invalid_input(format!("cdna_allele_comparison.{}", message.to_string()))
}
fn hash<T: Serialize>(value: &T) -> Result<String, EngineError> {
    let value = serde_json::to_value(value).map_err(fail)?;
    Ok(crate::digest_utils::sha256_prefixed_bytes(
        value.to_string().as_bytes(),
    ))
}
fn valid_hash(value: &str) -> bool {
    value
        .strip_prefix("sha256:")
        .is_some_and(|s| s.len() == 64 && s.bytes().all(|b| b.is_ascii_hexdigit()))
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CdnaLocusReference {
    pub assembly: String,
    pub contig: String,
    pub start_0based: usize,
    /// Declared whole-reference identity, distinct from the verified locus FASTA bytes.
    pub genomic_reference_sha256: String,
    pub locus_fasta: PrimerPayloadFile,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CdnaComparisonPartitionFiles {
    pub analysis: ObservedSamResultBinding,
    pub basecalls: PrimerPayloadFile,
}

/// One biological sample, optionally assembled from bounded disjoint producer partitions.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CdnaComparisonSampleFiles {
    pub sample_id: String,
    pub label: String,
    pub analysis: ObservedSamResultBinding,
    pub basecalls: PrimerPayloadFile,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub additional_partitions: Vec<CdnaComparisonPartitionFiles>,
    /// Declared original BAM identity shared by every partition's conversion receipt.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub partition_source_sha256: Option<String>,
}

impl CdnaComparisonSampleFiles {
    fn partitions(&self) -> impl Iterator<Item = (&ObservedSamResultBinding, &PrimerPayloadFile)> {
        std::iter::once((&self.analysis, &self.basecalls)).chain(
            self.additional_partitions
                .iter()
                .map(|p| (&p.analysis, &p.basecalls)),
        )
    }
}

/// File-only counterpart of the existing allele-screen module's pure cDNA comparison request.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CdnaFileComparisonRequest {
    pub schema: String,
    pub gene: String,
    pub source_id: String,
    pub gene_id: String,
    pub reference: CdnaLocusReference,
    pub normalized_annotation: PrimerPayloadFile,
    pub samples: Vec<CdnaComparisonSampleFiles>,
    #[serde(default)]
    pub supplied_sites: Vec<usize>,
    /// Optional supplied genomic marker basis for local phasing; discovery remains unchanged.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub phase_marker_positions: Option<Vec<usize>>,
    #[serde(default)]
    pub policy: core::CdnaComparisonPolicy,
    #[serde(default)]
    pub budgets: core::CdnaComparisonBudgets,
    #[serde(default = "default_input_bytes")]
    pub max_input_bytes: u64,
    #[serde(default = "default_basecall_bases")]
    pub max_basecall_bases: usize,
}

/// Original sequencing orientation and qualities are producer declarations bound to explicit units.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CdnaBasecallDocument {
    pub schema: String,
    pub analysis_context_sha256: String,
    pub orientation: String,
    pub provenance: String,
    pub reads: Vec<CdnaOriginalBasecall>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CdnaOriginalBasecall {
    pub unit: ObservedSamQueryUnitKey,
    pub read_id: String,
    pub sequence: String,
    /// Exact Phred+33 bytes in original sequencing order. None remains unavailable.
    pub qualities_phred33: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CdnaComparisonInputBinding {
    pub sample_id: String,
    pub partition_index: usize,
    pub report_content_sha256: String,
    pub analysis_data: PrimerPayloadFile,
    pub basecalls: PrimerPayloadFile,
    pub basecall_provenance: String,
    pub source_eof_verified: bool,
    /// Consensus applies to retained supplied hypotheses, not omitted alignment alternatives.
    pub inference_scope: String,
    pub units_with_incomplete_search_declaration: usize,
    pub units_with_incomplete_view_declaration: usize,
    pub units_with_unassessed_supplied_views: usize,
}

/// Portable read-only comparison plus verified raw inputs; never a DNA genotype or CNV diagnosis.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CdnaFileComparisonReport {
    pub schema: String,
    pub content_sha256: String,
    pub request_sha256: String,
    pub request: CdnaFileComparisonRequest,
    pub inputs: Vec<CdnaComparisonInputBinding>,
    pub comparison: core::CdnaComparisonReport,
    pub limitations: Vec<String>,
}

fn reserve(remaining: &mut u64, bytes: u64) -> Result<(), EngineError> {
    *remaining = remaining
        .checked_sub(bytes)
        .ok_or_else(|| fail("cumulative_input_budget"))?;
    Ok(())
}

pub(super) fn locus_sequence(
    file: &PrimerPayloadFile,
    contig: &str,
    check: &mut dyn FnMut() -> Result<(), EngineError>,
) -> Result<String, EngineError> {
    crate::observed_transcript_structures::indexed_locus::local_file(&file.path)?;
    if file.bytes == 0 || file.bytes > 16_000_000 {
        return Err(fail("locus_fasta_budget"));
    }
    let raw = crate::observed_transcript_structures::read_sam_bound_bytes(file, 16_000_000, check)?;
    let text = std::str::from_utf8(&raw).map_err(fail)?;
    let mut id = None;
    let mut sequence = String::new();
    for line in text.lines() {
        check()?;
        if let Some(header) = line.strip_prefix('>') {
            if id.is_some() {
                return Err(fail("locus_requires_one_fasta_record"));
            }
            id = header.split_whitespace().next();
        } else if !line.is_empty() {
            if id.is_none()
                || !line
                    .bytes()
                    .all(|b| b"ACGTNRYKMSWBDHVacgtnrykmswbdhv".contains(&b))
            {
                return Err(fail("locus_fasta_sequence"));
            }
            sequence.push_str(&line.to_ascii_uppercase());
        }
    }
    if id != Some(contig) || sequence.is_empty() {
        return Err(fail("locus_fasta_record"));
    }
    Ok(sequence)
}

/// Verify saved evidence, project bases through existing blocks, then call the allele-screen kernel.
/// No original SAM/BAM/FASTQ is reopened; no files or project content are written.
pub fn compare_cdna_files(
    request: CdnaFileComparisonRequest,
    progress: &mut dyn FnMut(OperationProgress) -> bool,
) -> Result<CdnaFileComparisonReport, EngineError> {
    if request.schema != CDNA_FILE_COMPARISON_REQUEST_SCHEMA
        || !(1..=16).contains(&request.samples.len())
        || !(1..=128_000_000).contains(&request.max_input_bytes)
        || !(1..=16_000_000).contains(&request.max_basecall_bases)
        || request.gene.is_empty()
        || request.gene.len() > 512
        || request.source_id.is_empty()
        || request.source_id.len() > 512
        || request.gene_id.is_empty()
        || request.gene_id.len() > 512
        || !valid_hash(&request.reference.genomic_reference_sha256)
    {
        return Err(fail("schema_or_budget"));
    }
    let mut check = || {
        if progress(OperationProgress::Workflow {
            completed: 0,
            total: request.samples.len(),
        }) {
            Ok(())
        } else {
            Err(fail("cancelled"))
        }
    };
    check()?;
    let mut remaining = request.max_input_bytes;
    reserve(&mut remaining, request.reference.locus_fasta.bytes)?;
    reserve(&mut remaining, request.normalized_annotation.bytes)?;
    for sample in &request.samples {
        if sample.additional_partitions.len() >= 64
            || (!sample.additional_partitions.is_empty()
                && sample.partition_source_sha256.is_none())
            || sample
                .partition_source_sha256
                .as_deref()
                .is_some_and(|h| !valid_hash(h) || h != h.to_ascii_lowercase())
        {
            return Err(fail("partition_schema_or_source"));
        }
        for (_, basecalls) in sample.partitions() {
            reserve(&mut remaining, basecalls.bytes)?;
        }
    }
    let reference_sequence = locus_sequence(
        &request.reference.locus_fasta,
        &request.reference.contig,
        &mut check,
    )?;
    let reference_end = request
        .reference
        .start_0based
        .checked_add(reference_sequence.len())
        .ok_or_else(|| fail("reference_interval"))?;
    let annotation: PrimerPayloadAnnotation =
        crate::rna_read_ends::bound_json(&request.normalized_annotation, &mut check)?;
    if annotation.schema != PRIMER_PAYLOAD_ANNOTATION_SCHEMA
        || annotation.assembly != request.reference.assembly
        || annotation.reference_fasta_sha256 != request.reference.genomic_reference_sha256
        || annotation.coordinate_convention != "zero_based_half_open"
        || annotation.transcripts.len() > 100_000
    {
        return Err(fail("annotation_reference"));
    }
    let mut transcript_ids = BTreeSet::new();
    let mut transcripts = Vec::new();
    for t in annotation
        .transcripts
        .iter()
        .filter(|t| t.gene_ids.contains(&request.gene_id))
    {
        if !transcript_ids.insert(&t.transcript_id)
            || t.genomic_contig.as_deref() != Some(&request.reference.contig)
            || t.target_name != request.reference.contig
            || !matches!(t.genomic_strand, Some('+' | '-'))
            || t.genomic_exons.is_empty()
            || t.genomic_exons.iter().any(|&(a, z)| {
                a >= z
                    || a < request.reference.start_0based
                    || z > reference_end
                    || z > t.target_length
            })
            || !t.genomic_exons.windows(2).all(|w| w[0].1 <= w[1].0)
        {
            return Err(fail("transcript_geometry"));
        }
        transcripts.push(core::CdnaTranscript {
            transcript_id: t.transcript_id.clone(),
            strand: t.genomic_strand.unwrap(),
            exons: t.genomic_exons.clone(),
        });
    }
    if transcripts.is_empty() || transcripts.len() > 256 {
        return Err(fail("transcript_universe"));
    }
    let mut samples = Vec::new();
    let mut inputs = Vec::new();
    let mut sample_ids = BTreeSet::new();
    let mut seen_analysis = BTreeSet::new();
    let mut seen_units = BTreeSet::new();
    let mut seen_partition_units = BTreeSet::new();
    let mut bases_left = request.max_basecall_bases;
    let mut units_left = 10_000usize;
    let mut projected_left = 2_000_000usize;
    for sample in &request.samples {
        check()?;
        if sample.sample_id.is_empty()
            || sample.sample_id.len() > 512
            || sample.label.is_empty()
            || sample.label.len() > 512
            || !sample_ids.insert(&sample.sample_id)
        {
            return Err(fail("sample_identity"));
        }
        let mut units = Vec::new();
        let mut identities = Vec::new();
        for (partition_index, (binding, basecalls)) in sample.partitions().enumerate() {
            crate::observed_transcript_structures::indexed_locus::local_file(&binding.report_path)?;
            let report_file = crate::observed_transcript_structures::indexed_locus::fingerprint(
                std::path::Path::new(&binding.report_path),
                remaining,
                &mut check,
            )?;
            reserve(&mut remaining, report_file.bytes)?;
            let preflight: crate::engine::ObservedSamReport =
                crate::rna_read_ends::bound_json(&report_file, &mut check)?;
            if preflight.content_sha256 != binding.report_content_sha256 {
                return Err(fail("report_preflight_binding"));
            }
            // Reserve all sidecar bytes before the existing inspector verifies them to EOF.
            for artifact in &preflight.artifacts {
                reserve(&mut remaining, artifact.bytes)?;
            }
            let (report, file) = crate::rna_read_ends::analysis_file(binding)?;
            if let Some(source) = &sample.partition_source_sha256 {
                if !matches!(&report.request.conversion,
                    ObservedSamConversion::BamConversion { bam_sha256, .. } if bam_sha256 == source)
                {
                    return Err(fail("partition_source_binding"));
                }
            }
            let data: ObservedSamData = crate::rna_read_ends::bound_json(&file, &mut check)?;
            let analysis = data.analysis.ok_or_else(|| fail("analysis_unavailable"))?;
            if data.schema != OBSERVED_SAM_DATA_SCHEMA
                || data.request_sha256 != report.request_sha256
                || report.semantic_context_sha256.as_deref() != Some(&analysis.context_sha256)
                || analysis.reference.assembly != request.reference.assembly
                || analysis.reference.fasta_sha256 != request.reference.genomic_reference_sha256
                || analysis.reference_gene.source_id != request.source_id
                || analysis.reference_gene.gene_id != request.gene_id
                || !analysis.annotation_bindings.iter().any(|a| {
                    a.source_id == request.source_id
                        && a.annotation_sha256 == request.normalized_annotation.sha256
                })
                || !seen_analysis.insert(file.sha256.clone())
            {
                return Err(fail("analysis_binding_or_duplicate"));
            }
            units_left = units_left
                .checked_sub(analysis.units.len())
                .ok_or_else(|| fail("unit_budget"))?;
            let calls: CdnaBasecallDocument =
                crate::rna_read_ends::bound_json(basecalls, &mut check)?;
            if calls.schema != CDNA_BASECALL_SCHEMA
                || calls.analysis_context_sha256 != analysis.context_sha256
                || calls.orientation != "original_sequencing_5_to_3"
                || calls.provenance.trim().is_empty()
                || calls.provenance.len() > 4096
                || calls.reads.len() > 16_384
            {
                return Err(fail("basecall_document"));
            }
            let unit_map = analysis
                .units
                .iter()
                .map(|u| (&u.unit, u))
                .collect::<BTreeMap<_, _>>();
            if unit_map.len() != analysis.units.len() {
                return Err(fail("duplicate_unit"));
            }
            let mut by_unit =
                BTreeMap::<ObservedSamQueryUnitKey, Vec<&CdnaOriginalBasecall>>::new();
            let mut sequences = BTreeSet::new();
            for call in &calls.reads {
                check()?;
                bases_left = bases_left
                    .checked_sub(call.sequence.len())
                    .ok_or_else(|| fail("basecall_base_budget"))?;
                let unit = unit_map
                    .get(&call.unit)
                    .ok_or_else(|| fail("foreign_basecall_unit"))?;
                if call.sequence.is_empty()
                    || unit.unclipped_axis_length != Some(call.sequence.len())
                    || unit.rows.is_empty()
                    || unit.rows.iter().any(|r| r.qname != call.read_id)
                    || call.read_id.is_empty()
                    || !call
                        .sequence
                        .bytes()
                        .all(|b| b"ACGTNRYKMSWBDHVacgtnrykmswbdhv".contains(&b))
                    || call.qualities_phred33.as_ref().is_some_and(|q| {
                        q.len() != call.sequence.len()
                            || !q.bytes().all(|b| (33..=126).contains(&b))
                    })
                    || !sequences.insert((&call.unit, call.sequence.to_ascii_uppercase()))
                {
                    return Err(fail("basecall_axis_quality_or_duplicate"));
                }
                let entry = by_unit.entry(call.unit.clone()).or_default();
                entry.push(call);
                if entry.len() > 16 {
                    return Err(fail("basecall_alternative_budget"));
                }
            }
            let obs_map = analysis
                .observations
                .iter()
                .map(|o| (&o.observation_id, o))
                .collect::<BTreeMap<_, _>>();
            if obs_map.len() != analysis.observations.len() || analysis.observations.len() > 100_000
            {
                return Err(fail("observation_population"));
            }
            for unit in &analysis.units {
                check()?;
                if !seen_units.insert(unit.unit.clone()) {
                    return Err(fail("unit_reused_across_samples"));
                }
                if let Some(source) = &sample.partition_source_sha256 {
                    if !seen_partition_units.insert((source.clone(), unit.unit.unit_id.clone())) {
                        return Err(fail("partition_unit_reused"));
                    }
                }
                let mut hypotheses = Vec::new();
                let variants = by_unit.get(&unit.unit);
                if unit
                    .hypotheses
                    .len()
                    .saturating_mul(variants.map_or(0, Vec::len))
                    > 64
                {
                    return Err(fail("hypothesis_budget"));
                }
                if let Some(variants) = variants {
                    for (hi, path) in unit.hypotheses.iter().enumerate() {
                        let mut seen_obs = BTreeSet::new();
                        for (vi, call) in variants.iter().enumerate() {
                            let mut bases = Vec::new();
                            let mut junctions = Vec::new();
                            for id in path {
                                if vi == 0 && !seen_obs.insert(id) {
                                    return Err(fail("repeated_hypothesis_observation"));
                                }
                                let obs = obs_map
                                    .get(id)
                                    .ok_or_else(|| fail("missing_hypothesis_observation"))?;
                                if obs.unit != unit.unit
                                    || !obs.on_optimal_path
                                    || obs.query_axis_length != call.sequence.len()
                                    || !matches!(obs.query_to_reference_strand, '+' | '-')
                                    || obs.aligned_query_interval.0 >= obs.aligned_query_interval.1
                                    || obs.aligned_query_interval.1 > call.sequence.len()
                                    || obs.observed_span.0 >= obs.observed_span.1
                                {
                                    return Err(fail("hypothesis_observation_binding"));
                                }
                                if obs.excluded_query_intervals.as_ref().is_some_and(|m| {
                                    m.iter().any(|&(a, z)| a >= z || z > call.sequence.len())
                                        || !m.windows(2).all(|w| w[0].1 <= w[1].0)
                                }) || obs.skips.iter().any(|&(a, z)| {
                                    a >= z || a < obs.observed_span.0 || z > obs.observed_span.1
                                }) || !obs.skips.windows(2).all(|w| w[0].1 <= w[1].0)
                                {
                                    return Err(fail("observation_masks_or_skips"));
                                }
                                if obs.blocks.iter().any(|block| {
                                    block.query.0 >= block.query.1
                                        || block.query.0 < obs.aligned_query_interval.0
                                        || block.query.1 > obs.aligned_query_interval.1
                                        || block.target.0 >= block.target.1
                                        || block.target.0 < obs.observed_span.0
                                        || block.target.1 > obs.observed_span.1
                                        || block.query.1 - block.query.0
                                            != block.target.1 - block.target.0
                                }) {
                                    return Err(fail("supported_block"));
                                }
                                if obs.contig != request.reference.contig {
                                    continue;
                                }
                                junctions.extend(obs.skips.iter().copied().filter(|&(a, z)| {
                                    a >= request.reference.start_0based && z <= reference_end
                                }));
                                for block in &obs.blocks {
                                    let start = block.target.0.max(request.reference.start_0based);
                                    let end = block.target.1.min(reference_end);
                                    for pos in start..end {
                                        projected_left = projected_left
                                            .checked_sub(1)
                                            .ok_or_else(|| fail("projection_budget"))?;
                                        if pos % 1024 == 0 {
                                            check()?;
                                        }
                                        let (q, z) = geometry::project(
                                            block,
                                            (pos, pos + 1),
                                            obs.query_to_reference_strand,
                                        )
                                        .ok_or_else(|| fail("projection"))?;
                                        if z != q + 1 {
                                            return Err(fail("projection_width"));
                                        }
                                        if obs
                                            .excluded_query_intervals
                                            .as_ref()
                                            .is_none_or(|m| m.iter().any(|&(a, z)| a <= q && q < z))
                                        {
                                            continue;
                                        }
                                        let original =
                                            call.sequence.as_bytes()[q].to_ascii_uppercase();
                                        let observed = if obs.query_to_reference_strand == '-' {
                                            crate::iupac_code::IupacCode::letter_complement(
                                                original,
                                            )
                                        } else {
                                            original
                                        };
                                        bases.push(core::CdnaBaseObservation {
                                            reference_position_0based: pos,
                                            query_position_0based: q,
                                            observed_base: observed as char,
                                            original_base: original as char,
                                            query_to_reference_strand: obs
                                                .query_to_reference_strand,
                                            phred: call
                                                .qualities_phred33
                                                .as_ref()
                                                .map(|s| s.as_bytes()[q] - 33),
                                            mapq: obs.reported_mapq,
                                        });
                                    }
                                }
                            }
                            junctions.sort();
                            junctions.dedup();
                            hypotheses.push(core::CdnaHypothesis {
                                hypothesis_id: format!("{hi}:{vi}"),
                                eligible: true,
                                bases,
                                junctions,
                            });
                        }
                    }
                }
                units.push(core::CdnaUnit {
                    unit_id: format!(
                        "{}:{}",
                        sample
                            .partition_source_sha256
                            .as_ref()
                            .unwrap_or(&unit.unit.source_sha256),
                        unit.unit.unit_id
                    ),
                    all_query_views_assessed: unit.all_query_views_assessed,
                    search_declared_complete: unit.search_declared_complete,
                    hypotheses,
                });
            }
            identities.push((file.sha256.clone(), analysis.context_sha256.clone()));
            inputs.push(CdnaComparisonInputBinding {
                sample_id: sample.sample_id.clone(),
                partition_index,
                report_content_sha256: report.content_sha256,
                analysis_data: file,
                basecalls: basecalls.clone(),
                basecall_provenance: calls.provenance,
                source_eof_verified: report.source_eof_verified,
                inference_scope: "supplied_eligible_mapping_hypotheses".into(),
                units_with_incomplete_search_declaration: analysis
                    .units
                    .iter()
                    .filter(|u| !u.search_declared_complete)
                    .count(),
                units_with_incomplete_view_declaration: analysis
                    .units
                    .iter()
                    .filter(|u| !u.views_declared_complete)
                    .count(),
                units_with_unassessed_supplied_views: analysis
                    .units
                    .iter()
                    .filter(|u| !u.all_query_views_assessed)
                    .count(),
            });
        }
        let (source_sha256, analysis_context_sha256) =
            if identities.len() == 1 && sample.partition_source_sha256.is_none() {
                identities[0].clone()
            } else {
                (
                    hash(&(
                        "cdna_partition_sources",
                        &sample.partition_source_sha256,
                        &identities,
                    ))?,
                    hash(&("cdna_partition_contexts", &identities))?,
                )
            };
        samples.push(core::CdnaSample {
            sample_id: sample.sample_id.clone(),
            label: sample.label.clone(),
            source_sha256,
            analysis_context_sha256,
            units,
        });
    }
    drop(check);
    let comparison = core::compare_cdna_alleles(
        core::CdnaComparisonRequest {
            schema: core::CDNA_COMPARISON_REQUEST_SCHEMA.into(),
            gene: request.gene.clone(),
            reference: core::CdnaReference {
                assembly: request.reference.assembly.clone(),
                contig: request.reference.contig.clone(),
                start_0based: request.reference.start_0based,
                sha256: crate::digest_utils::sha256_prefixed_bytes(reference_sequence.as_bytes()),
                sequence: reference_sequence,
            },
            transcripts,
            samples,
            supplied_sites: request.supplied_sites.clone(),
            phase_marker_positions: request.phase_marker_positions.clone(),
            policy: request.policy.clone(),
            budgets: request.budgets.clone(),
        },
        &mut |completed, total| progress(OperationProgress::Workflow { completed, total }),
    )
    .map_err(fail)?;
    let mut report=CdnaFileComparisonReport{schema:CDNA_FILE_COMPARISON_SCHEMA.into(),content_sha256:String::new(),request_sha256:hash(&request)?,request,inputs,comparison,
        limitations:vec!["Raw saved artifacts and supplied base/quality bytes are verified; original physical read identity and basecaller provenance remain producer declarations.".into(),
            "The parent genomic-reference hash is declared and cross-bound to saved annotation/analysis. Only supplied locus FASTA bytes/sequence are reverified here.".into(),
            "SNPs and haplotypes are expressed cDNA sequence candidates, not DNA genotypes. Unobserved/uncallable is not a reference allele; no second allele, molecule count, CNV or trans-splicing mechanism is invented.".into(),
            "Assessment of all supplied query views does not prove completeness of mapping alternatives. Search/view completeness declarations remain separate; inference is conditional on supplied eligible hypotheses, including primary-only selections.".into(),
            "Partitions are combined as one sample before SNP discovery and local phasing. Their original BAM hash and stable original unit IDs are producer declarations cross-bound to conversion receipts, not independently verified molecule identities.".into()]};
    report.content_sha256 = hash(&report)?;
    if serde_json::to_vec(&report).map_err(fail)?.len() > 8_000_000 {
        return Err(fail("response_budget"));
    }
    Ok(report)
}

#[cfg(test)]
mod tests;
