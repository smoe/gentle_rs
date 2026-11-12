//! Authenticated C2-C4 file service over unchanged cDNA comparisons and C1 models.
//! No aligner, live database, project mutation or filesystem publication.

use super::{coding_models, comparison};
use crate::engine::{EngineError, OperationProgress, PrimerPayloadFile};
use gentle_engine::allele_hash_screen::{cdna_comparison as evidence, cdna_interpretation as core};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};

#[cfg(test)]
mod tests;

pub const REQUEST_SCHEMA: &str = "gentle.cdna_allele_interpretation_request.v1";

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Linkage {
    #[default]
    None,
    SameUnit,
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Format {
    #[default]
    Json,
    Tsv,
    Summary,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct GeneFiles {
    pub comparison_report: PrimerPayloadFile,
    pub coding_models: PrimerPayloadFile,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Request {
    pub schema: String,
    pub genes: Vec<GeneFiles>,
    #[serde(default)]
    pub policy: core::Policy,
    #[serde(default)]
    pub linkage: Linkage,
    #[serde(default)]
    pub known_variants: Option<PrimerPayloadFile>,
    #[serde(default)]
    pub budgets: core::Budgets,
    #[serde(default)]
    pub format: Format,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct InputBinding {
    pub files: GeneFiles,
    pub comparison_form: String,
    pub comparison_content_sha256: String,
    pub coding_models_content_sha256: String,
    pub normalized_evidence_reverified: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Report {
    pub schema: String,
    #[serde(default, skip_serializing_if = "String::is_empty")]
    pub content_sha256: String,
    pub request_sha256: String,
    pub request: Request,
    pub inputs: Vec<InputBinding>,
    pub hgvs_version: String,
    pub genetic_code_sha256: String,
    pub genes: Vec<core::GeneReport>,
    pub known_variant_source: Option<BTreeMap<String, String>>,
    pub work: u64,
    pub limitations: Vec<String>,
}

fn fail(reason: impl std::fmt::Display) -> EngineError {
    let text = reason.to_string();
    EngineError::invalid_input(if text.starts_with("cdna_interpretation.") {
        text
    } else {
        format!("cdna_interpretation.{text}")
    })
}
fn digest(raw: &[u8]) -> String {
    crate::digest_utils::sha256_prefixed_bytes(raw)
}
fn canonical<T: Serialize>(value: &T) -> Result<String, EngineError> {
    Ok(digest(
        serde_json::to_value(value)
            .map_err(fail)?
            .to_string()
            .as_bytes(),
    ))
}
fn struct_hash<T: Serialize>(value: &T) -> Result<String, EngineError> {
    Ok(digest(&serde_json::to_vec(value).map_err(fail)?))
}

fn reserve(left: &mut u64, bytes: u64) -> Result<(), EngineError> {
    *left = left
        .checked_sub(bytes)
        .ok_or_else(|| fail("cumulative_input_budget"))?;
    Ok(())
}
fn read<T: serde::de::DeserializeOwned>(
    file: &PrimerPayloadFile,
    left: &mut u64,
    check: &mut dyn FnMut() -> Result<(), EngineError>,
) -> Result<T, EngineError> {
    reserve(left, file.bytes)?;
    crate::observed_transcript_structures::indexed_locus::local_file(&file.path)?;
    let bytes =
        crate::observed_transcript_structures::read_sam_bound_bytes(file, 128_000_000, check)?;
    serde_json::from_slice(&bytes).map_err(fail)
}

fn verify_comparison(report: &comparison::CdnaFileComparisonReport) -> Result<(), EngineError> {
    if report.schema != comparison::CDNA_FILE_COMPARISON_SCHEMA
        || report.comparison.schema != evidence::CDNA_COMPARISON_SCHEMA
    {
        return Err(fail("comparison_schema"));
    }
    let mut outer = report.clone();
    outer.content_sha256.clear();
    let mut inner = report.comparison.clone();
    inner.content_sha256.clear();
    if canonical(&outer)? != report.content_sha256
        || canonical(&report.request)? != report.request_sha256
        || struct_hash(&inner)? != report.comparison.content_sha256
    {
        return Err(fail("comparison_content_hash"));
    }
    if report.comparison.gene != report.request.gene
        || serde_json::to_value(&report.comparison.policy).map_err(fail)?
            != serde_json::to_value(&report.request.policy).map_err(fail)?
    {
        return Err(fail("comparison_request_binding"));
    }
    Ok(())
}

fn portable_models(
    c1: &coding_models::CodingTranscriptModelsReport,
    report: &comparison::CdnaFileComparisonReport,
    sequence: &str,
) -> Result<Vec<core::CodingModel>, EngineError> {
    if c1.schema != coding_models::CODING_MODELS_SCHEMA
        || c1.sequence_sha256 != report.comparison.reference.sequence_sha256
        || c1.sequence_sha256 != digest(sequence.as_bytes())
        || c1.request.selector.gene_id != report.request.gene_id
        || c1.request.selector.source_id != report.request.source_id
        || c1.request.reference.assembly != report.request.reference.assembly
        || c1.request.reference.contig != report.request.reference.contig
        || c1.request.reference.start_0based != report.request.reference.start_0based
        || c1.request.reference.genomic_reference_sha256
            != report.request.reference.genomic_reference_sha256
        || c1.request.normalized_annotation.sha256 != report.request.normalized_annotation.sha256
    {
        return Err(fail("reference_or_model_binding"));
    }
    let mut out = Vec::new();
    for m in &c1.models {
        let t = report
            .comparison
            .transcripts
            .iter()
            .find(|t| t.transcript_id == m.transcript.id)
            .ok_or_else(|| fail("model_geometry_mismatch"))?;
        if t.strand != m.strand || t.exons != m.genomic_exons_0based_half_open {
            return Err(fail("model_geometry_mismatch"));
        }
        let blocks =
            crate::observed_transcript_structures::read_ends::exon_map(&t.exons, t.strand)?;
        let mut transcript = Vec::new();
        for block in blocks {
            let a = block
                .target
                .0
                .checked_sub(report.request.reference.start_0based)
                .ok_or_else(|| fail("model_geometry_mismatch"))?;
            let z = block
                .target
                .1
                .checked_sub(report.request.reference.start_0based)
                .ok_or_else(|| fail("model_geometry_mismatch"))?;
            let raw = sequence
                .as_bytes()
                .get(a..z)
                .ok_or_else(|| fail("model_geometry_mismatch"))?;
            if t.strand == '+' {
                transcript.extend_from_slice(raw);
            } else {
                transcript.extend(
                    raw.iter()
                        .rev()
                        .map(|b| crate::iupac_code::IupacCode::letter_complement(*b)),
                );
            }
        }
        out.push(core::CodingModel {
            transcript: t.clone(),
            accession: m.transcript.versioned_accession.clone(),
            protein_accession: m
                .coding
                .as_ref()
                .and_then(|c| c.protein.as_ref())
                .and_then(|p| p.versioned_accession.clone()),
            sequence: String::from_utf8(transcript).map_err(fail)?,
            cds: m
                .coding
                .as_ref()
                .map(|c| c.normalized_coding_cdna_0based_half_open),
            translation_status: m.translation_status.clone(),
            translation_reasons: m.translation_reasons.clone(),
        });
    }
    Ok(out)
}

/// Authenticate bindings and delegate to the pure kernel. All results are read-only stdout data.
pub fn interpret_cdna_files(
    request: Request,
    progress: &mut dyn FnMut(OperationProgress) -> bool,
) -> Result<Report, EngineError> {
    if request.schema != REQUEST_SCHEMA
        || request.genes.is_empty()
        || request.genes.len() > 64
        || !(1..=128_000_000).contains(&request.budgets.max_input_bytes)
        || request.budgets.max_work == 0
        || request.budgets.max_work > 100_000_000
        || !(1..=64_000_000).contains(&request.budgets.max_output_bytes)
        || request.policy.translation_table != 1
        || request.policy.maximum_out_of_model_aligned_bases > 1_000_000
    {
        return Err(fail("schema_or_budget"));
    }
    let mut left = request.budgets.max_input_bytes;
    let mut work_left = request.budgets.max_work;
    let mut inputs = Vec::new();
    let mut genes = Vec::new();
    let mut seen = BTreeSet::new();
    for files in &request.genes {
        if !progress(OperationProgress::Workflow {
            completed: genes.len(),
            total: request.genes.len(),
        }) {
            return Err(fail("cancelled"));
        }
        let mut check = || {
            if progress(OperationProgress::Workflow {
                completed: genes.len(),
                total: request.genes.len(),
            }) {
                Ok(())
            } else {
                Err(fail("cancelled"))
            }
        };
        let raw: serde_json::Value = read(&files.comparison_report, &mut left, &mut check)?;
        let (raw, form) = if raw.get("schema").is_some() {
            (raw, "bare")
        } else {
            (
                raw.get("result")
                    .and_then(|r| r.get("cdna_allele_comparison"))
                    .cloned()
                    .ok_or_else(|| fail("comparison_envelope"))?,
                "cli_envelope",
            )
        };
        let comparison: comparison::CdnaFileComparisonReport =
            serde_json::from_value(raw).map_err(fail)?;
        verify_comparison(&comparison)?;
        if !seen.insert((
            comparison.request.source_id.clone(),
            comparison.request.gene_id.clone(),
        )) {
            return Err(fail("duplicate_gene"));
        }
        let c1: coding_models::CodingTranscriptModelsReport =
            read(&files.coding_models, &mut left, &mut check)?;
        let mut blank = c1.clone();
        blank.content_sha256.clear();
        if canonical(&blank)? != c1.content_sha256 || canonical(&c1.request)? != c1.request_sha256 {
            return Err(fail("coding_model_content_hash"));
        }
        for file in [
            &c1.request.gtf,
            &c1.request.normalized_annotation,
            &c1.request.reference.locus_fasta,
        ]
        .into_iter()
        .chain(c1.request.protein_fasta.iter())
        {
            reserve(&mut left, file.bytes)?;
        }
        drop(check);
        let verified = coding_models::derive_with_work_limit(
            c1.request.clone(),
            progress,
            usize::try_from(work_left).unwrap_or(usize::MAX),
        )?;
        if verified.content_sha256 != c1.content_sha256 {
            return Err(fail("coding_model_rederivation_mismatch"));
        }
        work_left = work_left
            .checked_sub(verified.work as u64)
            .ok_or_else(|| fail("work_budget"))?;
        let mut check = || {
            if progress(OperationProgress::Workflow {
                completed: genes.len(),
                total: request.genes.len(),
            }) {
                Ok(())
            } else {
                Err(fail("cancelled"))
            }
        };
        // C1 authenticated these identical locus bytes; re-reading keeps the shared FASTA grammar.
        reserve(&mut left, comparison.request.reference.locus_fasta.bytes)?;
        let sequence = comparison::locus_sequence(
            &comparison.request.reference.locus_fasta,
            &comparison.request.reference.contig,
            &mut check,
        )?;
        let models = portable_models(&c1, &comparison, &sequence)?;
        drop(check);
        if comparison.comparison.reference.assembly != comparison.request.reference.assembly
            || comparison.comparison.reference.contig != comparison.request.reference.contig
            || comparison.comparison.reference.start_0based
                != comparison.request.reference.start_0based
        {
            return Err(fail("comparison_reference_binding"));
        }
        let linkage = if request.linkage == Linkage::SameUnit {
            let lowered = comparison::lower_cdna_files(
                &comparison.request,
                progress,
                true,
                usize::try_from(work_left).unwrap_or(usize::MAX),
                left,
            )?;
            reserve(&mut left, lowered.input_bytes)?;
            work_left = work_left
                .checked_sub(lowered.projected_bases as u64)
                .ok_or_else(|| fail("work_budget"))?;
            if struct_hash(&lowered.request)? != comparison.comparison.request_sha256
                || canonical(&lowered.inputs)? != canonical(&comparison.inputs)?
            {
                return Err(fail("linkage_input_drift"));
            }
            Some(core::LinkageEvidence {
                samples: lowered.request.samples,
                geometry: lowered
                    .geometry
                    .into_iter()
                    .map(|(key, g)| {
                        (
                            key,
                            core::AlignmentGeometry {
                                blocks: g.blocks,
                                junctions: g.junctions,
                            },
                        )
                    })
                    .collect(),
            })
        } else {
            None
        };
        inputs.push(InputBinding {
            files: files.clone(),
            comparison_form: form.into(),
            comparison_content_sha256: comparison.content_sha256.clone(),
            coding_models_content_sha256: c1.content_sha256,
            normalized_evidence_reverified: linkage.is_some(),
        });
        genes.push(core::GeneInput {
            source_id: comparison.request.source_id.clone(),
            gene_id: comparison.request.gene_id.clone(),
            report: comparison.comparison,
            models,
            genomic_sequence: sequence,
            linkage,
        });
    }
    let mut check = || {
        if progress(OperationProgress::Workflow {
            completed: 0,
            total: request.genes.len(),
        }) {
            Ok(())
        } else {
            Err(fail("cancelled"))
        }
    };
    let known: Option<core::KnownVariantManifest> = request
        .known_variants
        .as_ref()
        .map(|f| read(f, &mut left, &mut check))
        .transpose()?;
    drop(check);
    let mut code = BTreeMap::new();
    for a in b"ACGT" {
        for b in b"ACGT" {
            for c in b"ACGT" {
                let aa = crate::AMINO_ACIDS.codon2aa([*a, *b, *c], Some(1));
                code.insert(
                    String::from_utf8(vec![*a, *b, *c]).unwrap(),
                    if aa == crate::amino_acids::STOP_CODON {
                        '*'
                    } else {
                        aa
                    },
                );
            }
        }
    }
    let mut budgets = request.budgets.clone();
    budgets.max_work = work_left;
    let (genes, work) = core::interpret(
        &genes,
        &request.policy,
        &code,
        known.as_ref(),
        &budgets,
        &mut || {
            progress(OperationProgress::Workflow {
                completed: 0,
                total: request.genes.len(),
            })
        },
    )
    .map_err(fail)?;
    let known_variant_source = known.as_ref().map(|k| {
        BTreeMap::from([
            ("source_name".into(), k.source_name.clone()),
            ("source_release".into(), k.source_release.clone()),
            ("retrieved_utc".into(), k.retrieved_utc.clone()),
            ("license".into(), k.license.clone()),
            ("attribution".into(), k.attribution.clone()),
        ])
    });
    let mut result=Report {schema:core::SCHEMA.into(),content_sha256:String::new(),request_sha256:canonical(&request)?,
        hgvs_version:core::HGVS_VERSION.into(),genetic_code_sha256:digest(include_bytes!("../../assets/codon_tables.json")),
        work:request.budgets.max_work-work_left+work,request,inputs,genes,known_variant_source,
        limitations:vec!["Conditional consequences of observed cDNA sequence, not DNA genotype, zygosity, expression, copy number, somatic/germline origin, pathogenicity or NMD.".into(),
            "Codon combinations use definite same-unit calls; local phase and partial transcript compatibility do not identify a complete molecule or expressed isoform.".into(),
            "Frozen known-variant labels are exact allele matches, not independent confirmation or clinical interpretation.".into()]};
    result.content_sha256 = canonical(&result)?;
    if serde_json::to_vec(&result).map_err(fail)?.len() > result.request.budgets.max_output_bytes {
        return Err(fail("output_budget"));
    }
    Ok(result)
}

fn escape(value: &str) -> String {
    value
        .replace('\\', "\\\\")
        .replace('\t', "\\t")
        .replace('\r', "\\r")
        .replace('\n', "\\n")
}

/// Pure projections from a saved JSON report: never reread source files.
pub fn render(report: &Report, format: Format) -> Result<String, EngineError> {
    if report.schema != core::SCHEMA {
        return Err(fail("projection_schema"));
    }
    let mut blank = report.clone();
    blank.content_sha256.clear();
    if canonical(&blank)? != report.content_sha256
        || canonical(&report.request)? != report.request_sha256
    {
        return Err(fail("projection_content_hash"));
    }
    let text = match format {
        Format::Json => serde_json::to_string_pretty(report).map_err(fail),
        Format::Summary => {
            let mut out = String::from("# Conditional cDNA interpretation\n\n");
            for text in &report.limitations {
                out.push_str(&format!("- {}\n", escape(text)));
            }
            for gene in &report.genes {
                out.push_str(&format!(
                    "\n## {}\n\n{} site/transcript rows; {} same-unit codon combinations.\n",
                    escape(&gene.gene),
                    gene.sites.len(),
                    gene.codon_combinations.len()
                ));
                for row in &gene.sites {
                    out.push_str(&format!(
                        "- {} / {}: {}>{}, {}/{} units, {}; {}\n",
                        escape(&row.sample_id),
                        escape(&row.transcript_id),
                        row.genomic_reference_base,
                        row.genomic_alternate_base,
                        row.alternate_units,
                        row.callable_units,
                        escape(&row.support_status),
                        escape(
                            row.consequence
                                .protein_hgvs
                                .as_deref()
                                .or(row.consequence.nucleotide_hgvs.as_deref())
                                .unwrap_or("not_assessable")
                        )
                    ));
                }
                for group in &gene.codon_groups {
                    out.push_str(&format!(
                        "- Codon group {} / {}: {} complete units, {} incomplete units.\n",
                        escape(&group.sample_id),
                        escape(&group.transcript_id),
                        group.complete_units,
                        group.incomplete_units.len()
                    ));
                }
                for c in &gene.codon_combinations {
                    out.push_str(&format!(
                        "- Combination {} / {} [{}]: {} units, {}; {}\n",
                        escape(&c.sample_id),
                        escape(&c.transcript_id),
                        c.genomic_alleles.iter().collect::<String>(),
                        c.unit_ids.len(),
                        if c.supported {
                            "supported"
                        } else {
                            "below inherited threshold"
                        },
                        escape(
                            c.consequence
                                .protein_hgvs
                                .as_deref()
                                .or(c.consequence.nucleotide_hgvs.as_deref())
                                .unwrap_or("notation unavailable")
                        )
                    ));
                }
            }
            Ok(out)
        }
        Format::Tsv => {
            let mut out = String::from(
                "# Conditional cDNA sequence evidence; no biological/clinical confirmation\n",
            );
            for limitation in &report.limitations {
                out.push_str(&format!("# {}\n", escape(limitation)));
            }
            out.push_str("record\tgene\tsample\ttranscript\tpositions_0based\talt\tunits\tcallable\tstatus\tnucleotide_hgvs\tprotein_hgvs\tconsequence\tdetails_json\n");
            for g in &report.genes {
                for group in &g.codon_groups {
                    let fields = vec![
                        "codon_group".into(),
                        g.gene.clone(),
                        group.sample_id.clone(),
                        group.transcript_id.clone(),
                        group
                            .positions_0based
                            .iter()
                            .map(usize::to_string)
                            .collect::<Vec<_>>()
                            .join(","),
                        String::new(),
                        group.complete_units.to_string(),
                        String::new(),
                        if group.incomplete_units.is_empty() {
                            "complete_patterns"
                        } else {
                            "incomplete_units_retained"
                        }
                        .into(),
                        String::new(),
                        String::new(),
                        String::new(),
                        serde_json::to_string(group).map_err(fail)?,
                    ];
                    out.push_str(
                        &fields
                            .iter()
                            .map(|v| escape(v))
                            .collect::<Vec<_>>()
                            .join("\t"),
                    );
                    out.push('\n');
                }
                for r in &g.sites {
                    let fields = vec![
                        "site".into(),
                        g.gene.clone(),
                        r.sample_id.clone(),
                        r.transcript_id.clone(),
                        r.genomic_position_0based.to_string(),
                        r.genomic_alternate_base.to_string(),
                        r.alternate_units.to_string(),
                        r.callable_units.to_string(),
                        r.support_status.clone(),
                        r.consequence.nucleotide_hgvs.clone().unwrap_or_default(),
                        r.consequence.protein_hgvs.clone().unwrap_or_default(),
                        r.consequence.consequence_terms.join(","),
                        serde_json::to_string(r).map_err(fail)?,
                    ];
                    out.push_str(
                        &fields
                            .iter()
                            .map(|v| escape(v))
                            .collect::<Vec<_>>()
                            .join("\t"),
                    );
                    out.push('\n');
                }
                for c in &g.codon_combinations {
                    let fields = vec![
                        "codon_combination".into(),
                        g.gene.clone(),
                        c.sample_id.clone(),
                        c.transcript_id.clone(),
                        c.positions_0based
                            .iter()
                            .map(usize::to_string)
                            .collect::<Vec<_>>()
                            .join(","),
                        c.genomic_alleles.iter().collect::<String>(),
                        c.unit_ids.len().to_string(),
                        String::new(),
                        if c.supported {
                            "supported"
                        } else {
                            "below_inherited_threshold"
                        }
                        .into(),
                        c.consequence.nucleotide_hgvs.clone().unwrap_or_default(),
                        c.consequence.protein_hgvs.clone().unwrap_or_default(),
                        c.consequence.consequence_terms.join(","),
                        serde_json::to_string(c).map_err(fail)?,
                    ];
                    out.push_str(
                        &fields
                            .iter()
                            .map(|v| escape(v))
                            .collect::<Vec<_>>()
                            .join("\t"),
                    );
                    out.push('\n');
                }
            }
            Ok(out)
        }
    }?;
    if text.len() > report.request.budgets.max_output_bytes {
        return Err(fail("projection_output_budget"));
    }
    Ok(text)
}
