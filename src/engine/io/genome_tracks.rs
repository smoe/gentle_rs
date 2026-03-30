//! Genome track import/parsing helpers for BED/BigWig/VCF and BLAST features.
//!
//! Parsers and generated-feature cleanup live here so every track-ingestion path
//! reuses the same normalization rules.
//!
//! Look here for:
//! - text/track reader helpers
//! - generated-feature tagging and cleanup for imported overlays
//! - BED/BigWig/VCF/BLAST projection rules shared by GUI, shell, and CLI

use super::*;

impl GentleEngine {
    pub(super) fn is_generated_genome_bed_feature(feature: &gb_io::seq::Feature) -> bool {
        feature
            .qualifier_values("gentle_generated".into())
            .any(|v| v.eq_ignore_ascii_case(GENOME_BED_TRACK_GENERATED_TAG))
    }

    pub(super) fn is_generated_genome_bigwig_feature(feature: &gb_io::seq::Feature) -> bool {
        feature
            .qualifier_values("gentle_generated".into())
            .any(|v| v.eq_ignore_ascii_case(GENOME_BIGWIG_TRACK_GENERATED_TAG))
    }

    pub(super) fn is_generated_genome_vcf_feature(feature: &gb_io::seq::Feature) -> bool {
        feature
            .qualifier_values("gentle_generated".into())
            .any(|v| v.eq_ignore_ascii_case(GENOME_VCF_TRACK_GENERATED_TAG))
    }

    pub(super) fn is_generated_blast_hit_feature(feature: &gb_io::seq::Feature) -> bool {
        feature
            .qualifier_values("gentle_generated".into())
            .any(|v| v.eq_ignore_ascii_case(BLAST_HIT_TRACK_GENERATED_TAG))
    }

    pub(super) fn is_generated_genome_signal_feature(feature: &gb_io::seq::Feature) -> bool {
        Self::is_generated_genome_bed_feature(feature)
            || Self::is_generated_genome_bigwig_feature(feature)
            || Self::is_generated_genome_vcf_feature(feature)
    }

    pub(super) fn remove_generated_genome_signal_features(features: &mut Vec<gb_io::seq::Feature>) {
        features.retain(|f| !Self::is_generated_genome_signal_feature(f));
    }

    pub(super) fn remove_generated_blast_hit_features(features: &mut Vec<gb_io::seq::Feature>) {
        features.retain(|f| !Self::is_generated_blast_hit_feature(f));
    }

    pub(super) fn open_text_reader(path: &str) -> Result<Box<dyn BufRead>, EngineError> {
        let file = std::fs::File::open(path).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!("Could not open track file '{path}': {e}"),
        })?;
        let lower = path.to_ascii_lowercase();
        if lower.ends_with(".gz") {
            let decoder = MultiGzDecoder::new(BufReader::new(file));
            Ok(Box::new(BufReader::new(decoder)))
        } else {
            Ok(Box::new(BufReader::new(file)))
        }
    }

    pub(super) fn default_track_name(path: &str) -> String {
        let file_name = Path::new(path)
            .file_name()
            .and_then(|s| s.to_str())
            .unwrap_or("track");
        let without_gz = file_name
            .strip_suffix(".gz")
            .or_else(|| file_name.strip_suffix(".GZ"))
            .unwrap_or(file_name);
        let mut base = without_gz.to_string();
        let lower = base.to_ascii_lowercase();
        for suffix in [".bed", ".bigwig", ".bw", ".vcf"] {
            if lower.ends_with(suffix) && base.len() > suffix.len() {
                base.truncate(base.len() - suffix.len());
                break;
            }
        }
        if base.trim().is_empty() {
            "track".to_string()
        } else {
            base
        }
    }

    pub(super) fn parse_bed_record(line: &str) -> Result<BedRecord, String> {
        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() < 3 {
            return Err("BED record needs at least 3 columns (chrom, start, end)".to_string());
        }
        let chromosome = fields[0].trim();
        if chromosome.is_empty() {
            return Err("BED chromosome is empty".to_string());
        }
        let start_0based = fields[1]
            .parse::<usize>()
            .map_err(|e| format!("Invalid BED start '{}': {e}", fields[1]))?;
        let end_0based = fields[2]
            .parse::<usize>()
            .map_err(|e| format!("Invalid BED end '{}': {e}", fields[2]))?;
        if end_0based <= start_0based {
            return Err(format!(
                "Invalid BED interval {}:{}-{} (end must be > start)",
                chromosome, start_0based, end_0based
            ));
        }
        let name = fields
            .get(3)
            .map(|v| v.trim())
            .filter(|v| !v.is_empty() && *v != ".")
            .map(|v| v.to_string());
        let score = fields
            .get(4)
            .map(|v| v.trim())
            .filter(|v| !v.is_empty() && *v != ".")
            .map(|v| {
                v.parse::<f64>()
                    .map_err(|e| format!("Invalid BED score '{}': {e}", v))
            })
            .transpose()?;
        let strand = fields
            .get(5)
            .and_then(|v| v.trim().chars().next())
            .filter(|c| matches!(c, '+' | '-'));

        Ok(BedRecord {
            chromosome: chromosome.to_string(),
            start_0based,
            end_0based,
            name,
            score,
            strand,
        })
    }

    pub(super) fn parse_bedgraph_record(line: &str) -> Result<BedRecord, String> {
        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() < 4 {
            return Err(
                "bedGraph record needs at least 4 columns (chrom, start, end, value)".to_string(),
            );
        }
        let chromosome = fields[0].trim();
        if chromosome.is_empty() {
            return Err("bedGraph chromosome is empty".to_string());
        }
        let start_0based = fields[1]
            .parse::<usize>()
            .map_err(|e| format!("Invalid bedGraph start '{}': {e}", fields[1]))?;
        let end_0based = fields[2]
            .parse::<usize>()
            .map_err(|e| format!("Invalid bedGraph end '{}': {e}", fields[2]))?;
        if end_0based <= start_0based {
            return Err(format!(
                "Invalid bedGraph interval {}:{}-{} (end must be > start)",
                chromosome, start_0based, end_0based
            ));
        }
        let value = fields[3].parse::<f64>().map_err(|e| {
            format!(
                "Invalid bedGraph value '{}' (expected numeric): {e}",
                fields[3]
            )
        })?;
        Ok(BedRecord {
            chromosome: chromosome.to_string(),
            start_0based,
            end_0based,
            name: None,
            score: Some(value),
            strand: None,
        })
    }

    pub(super) fn parse_vcf_record(line: &str) -> Result<VcfRecord, String> {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 8 {
            return Err(
                "VCF record needs at least 8 columns (CHROM POS ID REF ALT QUAL FILTER INFO)"
                    .to_string(),
            );
        }
        let chromosome = fields[0].trim();
        if chromosome.is_empty() {
            return Err("VCF chromosome is empty".to_string());
        }
        let pos_1based = fields[1]
            .trim()
            .parse::<usize>()
            .map_err(|e| format!("Invalid VCF POS '{}': {e}", fields[1].trim()))?;
        if pos_1based == 0 {
            return Err("VCF POS must be >= 1".to_string());
        }
        let id = fields
            .get(2)
            .map(|v| v.trim())
            .filter(|v| !v.is_empty() && *v != ".")
            .map(|v| v.to_string());
        let reference = fields
            .get(3)
            .map(|v| v.trim())
            .filter(|v| !v.is_empty() && *v != ".")
            .ok_or_else(|| "VCF REF must be non-empty".to_string())?
            .to_string();
        let alternates = fields
            .get(4)
            .map(|v| v.trim())
            .ok_or_else(|| "VCF ALT must be present".to_string())?
            .split(',')
            .map(str::trim)
            .filter(|v| !v.is_empty() && *v != ".")
            .map(|v| v.to_string())
            .collect::<Vec<_>>();
        if alternates.is_empty() {
            return Err("VCF ALT must contain at least one allele".to_string());
        }
        let qual = fields
            .get(5)
            .map(|v| v.trim())
            .filter(|v| !v.is_empty() && *v != ".")
            .map(|v| {
                v.parse::<f64>()
                    .map_err(|e| format!("Invalid VCF QUAL '{}': {e}", v))
            })
            .transpose()?;
        let filter = fields
            .get(6)
            .map(|v| v.trim())
            .filter(|v| !v.is_empty() && *v != ".")
            .map(|v| v.to_string());
        let info = fields
            .get(7)
            .map(|v| v.trim())
            .filter(|v| !v.is_empty() && *v != ".")
            .map(|v| v.to_string());
        let format = fields
            .get(8)
            .map(|v| v.trim())
            .filter(|v| !v.is_empty() && *v != ".")
            .map(|v| v.to_string());
        let sample_columns = if fields.len() > 9 {
            fields[9..]
                .iter()
                .map(|value| value.trim().to_string())
                .collect::<Vec<_>>()
        } else {
            vec![]
        };

        Ok(VcfRecord {
            chromosome: chromosome.to_string(),
            pos_1based,
            id,
            reference,
            alternates,
            qual,
            filter,
            info,
            format,
            sample_columns,
        })
    }

    pub(super) fn classify_vcf_alt(reference: &str, alt: &str) -> VcfVariantClass {
        let ref_u = reference.trim().to_ascii_uppercase();
        let alt_u = alt.trim().to_ascii_uppercase();
        if alt_u.is_empty() || alt_u == "." {
            return VcfVariantClass::Other;
        }
        if alt_u.starts_with('<')
            || alt_u.ends_with('>')
            || alt_u.contains('[')
            || alt_u.contains(']')
            || alt_u == "*"
        {
            return VcfVariantClass::Sv;
        }
        let ref_len = ref_u.len().max(1);
        let alt_len = alt_u.len().max(1);
        if ref_len == 1 && alt_len == 1 {
            return VcfVariantClass::Snp;
        }
        if alt_len > ref_len {
            return VcfVariantClass::Ins;
        }
        if alt_len < ref_len {
            return VcfVariantClass::Del;
        }
        VcfVariantClass::Other
    }

    pub(super) fn summarize_vcf_alt_genotype(
        record: &VcfRecord,
        alt_allele_index_1based: usize,
        sample_names: &[String],
    ) -> Option<VcfAltGenotypeSummary> {
        let format = record.format.as_deref()?;
        let format_fields = format.split(':').collect::<Vec<_>>();
        let gt_idx = format_fields
            .iter()
            .position(|field| field.eq_ignore_ascii_case("GT"))?;
        let mut summary = VcfAltGenotypeSummary::default();
        for (sample_idx, sample_column) in record.sample_columns.iter().enumerate() {
            let value_fields = sample_column.split(':').collect::<Vec<_>>();
            let Some(raw_gt) = value_fields.get(gt_idx).map(|v| v.trim()) else {
                continue;
            };
            if raw_gt.is_empty() || raw_gt == "." {
                continue;
            }
            let (tokens, phased) = if raw_gt.contains('|') {
                (raw_gt.split('|').collect::<Vec<_>>(), true)
            } else if raw_gt.contains('/') {
                (raw_gt.split('/').collect::<Vec<_>>(), false)
            } else {
                (vec![raw_gt], false)
            };
            let mut parsed = vec![];
            for token in tokens {
                let trimmed = token.trim();
                if trimmed.is_empty() || trimmed == "." {
                    continue;
                }
                if let Ok(index) = trimmed.parse::<usize>() {
                    parsed.push(index);
                }
            }
            if parsed.is_empty() {
                continue;
            }
            if !parsed.iter().any(|index| *index == alt_allele_index_1based) {
                continue;
            }
            summary.carriers += 1;
            if phased {
                summary.phased_carriers += 1;
            } else {
                summary.unphased_carriers += 1;
            }
            let has_ref = parsed.contains(&0);
            let has_other_alt = parsed
                .iter()
                .any(|idx| *idx > 0 && *idx != alt_allele_index_1based);
            let all_target_alt = parsed.iter().all(|idx| *idx == alt_allele_index_1based);
            if parsed.len() == 1 && all_target_alt {
                summary.haploid_alt += 1;
            } else if all_target_alt {
                summary.hom_alt += 1;
            } else if has_ref {
                summary.het += 1;
            } else if has_other_alt {
                summary.mixed_alt += 1;
            } else {
                summary.het += 1;
            }
            let sample_label = sample_names
                .get(sample_idx)
                .filter(|v| !v.trim().is_empty())
                .cloned()
                .unwrap_or_else(|| format!("sample_{}", sample_idx + 1));
            summary.carrier_samples.push(sample_label);
        }
        Some(summary)
    }

    pub(super) fn build_genome_signal_feature(
        record: &BedRecord,
        track_name: &str,
        path: &str,
        local_start_0based: usize,
        local_end_0based_exclusive: usize,
        local_strand: Option<char>,
        generated_tag: &str,
        source_label: &str,
    ) -> gb_io::seq::Feature {
        let label = record
            .name
            .as_ref()
            .filter(|v| !v.trim().is_empty())
            .cloned()
            .unwrap_or_else(|| {
                format!(
                    "{}:{}:{}-{}",
                    track_name, record.chromosome, record.start_0based, record.end_0based
                )
            });
        let mut qualifiers = vec![
            ("label".into(), Some(label.clone())),
            (
                "note".into(),
                Some(format!(
                    "{} track '{}' from '{}' [{}:{}-{}]",
                    source_label,
                    track_name,
                    path,
                    record.chromosome,
                    record.start_0based,
                    record.end_0based
                )),
            ),
            ("gentle_track_source".into(), Some(source_label.to_string())),
            ("gentle_generated".into(), Some(generated_tag.to_string())),
            ("gentle_track_name".into(), Some(track_name.to_string())),
            ("gentle_track_file".into(), Some(path.to_string())),
            ("chromosome".into(), Some(record.chromosome.clone())),
            (
                "bed_start_0based".into(),
                Some(record.start_0based.to_string()),
            ),
            ("bed_end_0based".into(), Some(record.end_0based.to_string())),
        ];
        if let Some(score) = record.score {
            qualifiers.push(("score".into(), Some(format!("{score:.6}"))));
        }
        if let Some(strand) = record.strand {
            qualifiers.push(("bed_strand".into(), Some(strand.to_string())));
        }
        if let Some(strand) = local_strand {
            qualifiers.push(("strand".into(), Some(strand.to_string())));
        }

        let base_location = gb_io::seq::Location::simple_range(
            local_start_0based as i64,
            local_end_0based_exclusive as i64,
        );
        let location = if local_strand == Some('-') {
            gb_io::seq::Location::Complement(Box::new(base_location))
        } else {
            base_location
        };

        gb_io::seq::Feature {
            kind: gb_io::FeatureKind::from("track"),
            location,
            qualifiers,
        }
    }

    pub(super) fn build_genome_bed_feature(
        record: &BedRecord,
        track_name: &str,
        path: &str,
        local_start_0based: usize,
        local_end_0based_exclusive: usize,
        local_strand: Option<char>,
    ) -> gb_io::seq::Feature {
        Self::build_genome_signal_feature(
            record,
            track_name,
            path,
            local_start_0based,
            local_end_0based_exclusive,
            local_strand,
            GENOME_BED_TRACK_GENERATED_TAG,
            "BED",
        )
    }

    pub(super) fn build_genome_vcf_feature(
        record: &VcfRecord,
        alt: &str,
        alt_allele_index_1based: usize,
        sample_names: &[String],
        track_name: &str,
        path: &str,
        local_start_0based: usize,
        local_end_0based_exclusive: usize,
    ) -> gb_io::seq::Feature {
        let variant_class = Self::classify_vcf_alt(&record.reference, alt);
        let label = record.id.clone().unwrap_or_else(|| {
            format!(
                "{}:{}:{} {}>{}",
                track_name, record.chromosome, record.pos_1based, record.reference, alt
            )
        });
        let mut qualifiers = vec![
            ("label".into(), Some(label.clone())),
            (
                "note".into(),
                Some(format!(
                    "VCF track '{}' from '{}' [{}:{} {}>{}]",
                    track_name, path, record.chromosome, record.pos_1based, record.reference, alt
                )),
            ),
            ("gentle_track_source".into(), Some("VCF".to_string())),
            (
                "gentle_generated".into(),
                Some(GENOME_VCF_TRACK_GENERATED_TAG.to_string()),
            ),
            ("gentle_track_name".into(), Some(track_name.to_string())),
            ("gentle_track_file".into(), Some(path.to_string())),
            ("chromosome".into(), Some(record.chromosome.clone())),
            ("vcf_pos_1based".into(), Some(record.pos_1based.to_string())),
            ("vcf_ref".into(), Some(record.reference.clone())),
            ("vcf_alt".into(), Some(alt.to_string())),
            (
                "vcf_alt_allele_index".into(),
                Some(alt_allele_index_1based.to_string()),
            ),
            (
                "vcf_variant_class".into(),
                Some(variant_class.as_str().to_string()),
            ),
        ];
        if let Some(id) = &record.id {
            qualifiers.push(("vcf_id".into(), Some(id.clone())));
        }
        if let Some(qual) = record.qual {
            qualifiers.push(("score".into(), Some(format!("{qual:.6}"))));
            qualifiers.push(("vcf_qual".into(), Some(format!("{qual:.6}"))));
        }
        if let Some(filter) = &record.filter {
            qualifiers.push(("vcf_filter".into(), Some(filter.clone())));
        }
        if let Some(info) = &record.info {
            qualifiers.push(("vcf_info".into(), Some(info.clone())));
        }
        if let Some(format) = &record.format {
            qualifiers.push(("vcf_format".into(), Some(format.clone())));
        }
        if !record.sample_columns.is_empty() {
            qualifiers.push((
                "vcf_sample_count".into(),
                Some(record.sample_columns.len().to_string()),
            ));
        }
        if let Some(genotype) =
            Self::summarize_vcf_alt_genotype(record, alt_allele_index_1based, sample_names)
        {
            qualifiers.push((
                "vcf_alt_carriers".into(),
                Some(genotype.carriers.to_string()),
            ));
            qualifiers.push((
                "vcf_alt_carrier_phased".into(),
                Some(genotype.phased_carriers.to_string()),
            ));
            qualifiers.push((
                "vcf_alt_carrier_unphased".into(),
                Some(genotype.unphased_carriers.to_string()),
            ));
            qualifiers.push(("vcf_gt_het".into(), Some(genotype.het.to_string())));
            qualifiers.push(("vcf_gt_hom_alt".into(), Some(genotype.hom_alt.to_string())));
            qualifiers.push((
                "vcf_gt_mixed_alt".into(),
                Some(genotype.mixed_alt.to_string()),
            ));
            qualifiers.push((
                "vcf_gt_haploid_alt".into(),
                Some(genotype.haploid_alt.to_string()),
            ));
            let zygosity = if genotype.hom_alt > 0
                && genotype.het == 0
                && genotype.mixed_alt == 0
                && genotype.haploid_alt == 0
            {
                "hom_alt"
            } else if genotype.het > 0
                && genotype.hom_alt == 0
                && genotype.mixed_alt == 0
                && genotype.haploid_alt == 0
            {
                "het"
            } else if genotype.haploid_alt > 0
                && genotype.hom_alt == 0
                && genotype.het == 0
                && genotype.mixed_alt == 0
            {
                "haploid_alt"
            } else if genotype.carriers > 0 {
                "mixed"
            } else {
                "none"
            };
            qualifiers.push(("vcf_zygosity".into(), Some(zygosity.to_string())));
            let phase = if genotype.phased_carriers > 0 && genotype.unphased_carriers == 0 {
                "phased"
            } else if genotype.unphased_carriers > 0 && genotype.phased_carriers == 0 {
                "unphased"
            } else if genotype.phased_carriers > 0 && genotype.unphased_carriers > 0 {
                "mixed"
            } else {
                "unknown"
            };
            qualifiers.push(("vcf_phase".into(), Some(phase.to_string())));
            if !genotype.carrier_samples.is_empty() {
                qualifiers.push((
                    "vcf_alt_carrier_samples".into(),
                    Some(
                        genotype
                            .carrier_samples
                            .iter()
                            .take(20)
                            .cloned()
                            .collect::<Vec<_>>()
                            .join(","),
                    ),
                ));
            }
        }

        gb_io::seq::Feature {
            kind: gb_io::FeatureKind::from("track"),
            location: gb_io::seq::Location::simple_range(
                local_start_0based as i64,
                local_end_0based_exclusive as i64,
            ),
            qualifiers,
        }
    }

    pub(super) fn build_blast_hit_feature(
        hit: &BlastHitFeatureInput,
        track_name: &str,
        local_start_0based: usize,
        local_end_0based_exclusive: usize,
        local_strand: Option<char>,
    ) -> gb_io::seq::Feature {
        let label = format!(
            "{}:{}:{}-{}",
            track_name, hit.subject_id, hit.query_start_1based, hit.query_end_1based
        );
        let mut qualifiers = vec![
            ("label".into(), Some(label)),
            (
                "note".into(),
                Some(format!(
                    "BLAST hit '{}' query={}..{} subject={}..{} identity={:.2}% bitscore={:.2} evalue={:.3e}",
                    hit.subject_id,
                    hit.query_start_1based,
                    hit.query_end_1based,
                    hit.subject_start_1based,
                    hit.subject_end_1based,
                    hit.identity_percent,
                    hit.bit_score,
                    hit.evalue
                )),
            ),
            ("gentle_track_source".into(), Some("BLAST".to_string())),
            (
                "gentle_generated".into(),
                Some(BLAST_HIT_TRACK_GENERATED_TAG.to_string()),
            ),
            ("gentle_track_name".into(), Some(track_name.to_string())),
            ("blast_subject_id".into(), Some(hit.subject_id.clone())),
            (
                "blast_query_start_1based".into(),
                Some(hit.query_start_1based.to_string()),
            ),
            (
                "blast_query_end_1based".into(),
                Some(hit.query_end_1based.to_string()),
            ),
            (
                "blast_subject_start_1based".into(),
                Some(hit.subject_start_1based.to_string()),
            ),
            (
                "blast_subject_end_1based".into(),
                Some(hit.subject_end_1based.to_string()),
            ),
            (
                "blast_identity_percent".into(),
                Some(format!("{:.6}", hit.identity_percent)),
            ),
            ("score".into(), Some(format!("{:.6}", hit.bit_score))),
            (
                "blast_bit_score".into(),
                Some(format!("{:.6}", hit.bit_score)),
            ),
            ("blast_evalue".into(), Some(format!("{:.3e}", hit.evalue))),
        ];
        if let Some(qcov) = hit.query_coverage_percent {
            qualifiers.push(("blast_qcov_percent".into(), Some(format!("{qcov:.4}"))));
        }
        if let Some(strand) = local_strand {
            qualifiers.push(("strand".into(), Some(strand.to_string())));
        }

        let base_location = gb_io::seq::Location::simple_range(
            local_start_0based as i64,
            local_end_0based_exclusive as i64,
        );
        let location = if local_strand == Some('-') {
            gb_io::seq::Location::Complement(Box::new(base_location))
        } else {
            base_location
        };

        gb_io::seq::Feature {
            kind: gb_io::FeatureKind::from("track"),
            location,
            qualifiers,
        }
    }

    pub(super) fn import_genome_bed_track(
        dna: &mut DNAsequence,
        anchor: &GenomeSequenceAnchor,
        path: &str,
        track_name: Option<&str>,
        min_score: Option<f64>,
        max_score: Option<f64>,
        clear_existing: bool,
        mut progress_cb: Option<&mut dyn FnMut(usize, usize, usize, bool) -> bool>,
    ) -> Result<GenomeBedTrackImportReport, EngineError> {
        let selected_track_name = track_name
            .map(str::trim)
            .filter(|v| !v.is_empty())
            .map(|v| v.to_string())
            .unwrap_or_else(|| Self::default_track_name(path));

        if clear_existing {
            Self::remove_generated_genome_signal_features(dna.features_mut());
        }

        let mut report = GenomeBedTrackImportReport {
            track_name: selected_track_name.clone(),
            ..Default::default()
        };
        let progress_stride = 250usize;
        let mut reader = Self::open_text_reader(path)?;
        let mut line = String::new();
        let mut line_no = 0usize;
        let mut mismatch_counts: HashMap<String, usize> = HashMap::new();

        while {
            line.clear();
            reader.read_line(&mut line).map_err(|e| EngineError {
                code: ErrorCode::Io,
                message: format!("Could not read BED file '{path}': {e}"),
            })? > 0
        } {
            line_no += 1;
            let trimmed = line.trim();
            if trimmed.is_empty() {
                continue;
            }
            if trimmed.starts_with('#')
                || trimmed.to_ascii_lowercase().starts_with("track ")
                || trimmed.to_ascii_lowercase().starts_with("browser ")
            {
                continue;
            }
            report.parsed_records += 1;
            if report.parsed_records % progress_stride == 0 {
                if let Some(cb) = progress_cb.as_mut() {
                    let should_continue = (**cb)(
                        report.parsed_records,
                        report.imported_features,
                        report.skipped_records,
                        false,
                    );
                    if !should_continue {
                        report.cancelled = true;
                        break;
                    }
                }
            }

            let record = match Self::parse_bed_record(trimmed) {
                Ok(v) => v,
                Err(e) => {
                    report.skipped_records += 1;
                    report.skipped_invalid += 1;
                    if report.warnings.len() < 20 {
                        report
                            .warnings
                            .push(format!("BED line {} skipped: {}", line_no, e));
                    }
                    continue;
                }
            };

            if !Self::chromosomes_match(&record.chromosome, &anchor.chromosome) {
                report.skipped_records += 1;
                report.skipped_wrong_chromosome += 1;
                *mismatch_counts
                    .entry(record.chromosome.clone())
                    .or_insert(0) += 1;
                continue;
            }

            let bed_start_1based = record.start_0based.saturating_add(1);
            let bed_end_1based = record.end_0based;
            if bed_end_1based < anchor.start_1based || bed_start_1based > anchor.end_1based {
                report.skipped_records += 1;
                report.skipped_non_overlap += 1;
                continue;
            }
            let overlap_start_1based = bed_start_1based.max(anchor.start_1based);
            let overlap_end_1based = bed_end_1based.min(anchor.end_1based);
            if overlap_end_1based < overlap_start_1based {
                report.skipped_records += 1;
                report.skipped_non_overlap += 1;
                continue;
            }

            if min_score.is_some() || max_score.is_some() {
                let Some(score) = record.score else {
                    report.skipped_records += 1;
                    report.skipped_missing_score += 1;
                    continue;
                };
                if min_score.map(|v| score < v).unwrap_or(false)
                    || max_score.map(|v| score > v).unwrap_or(false)
                {
                    report.skipped_records += 1;
                    report.skipped_outside_score_range += 1;
                    continue;
                }
            }

            let (local_start_0based, local_end_0based_exclusive) = if anchor.strand == Some('-') {
                (
                    anchor.end_1based.saturating_sub(overlap_end_1based),
                    anchor.end_1based.saturating_sub(overlap_start_1based) + 1,
                )
            } else {
                (
                    overlap_start_1based - anchor.start_1based,
                    overlap_end_1based - anchor.start_1based + 1,
                )
            };
            let local_strand = match (record.strand, anchor.strand) {
                (Some('+'), Some('-')) => Some('-'),
                (Some('-'), Some('-')) => Some('+'),
                (Some(strand), _) => Some(strand),
                _ => None,
            };
            if report.imported_features >= MAX_IMPORTED_SIGNAL_FEATURES {
                report.truncated_at_limit = true;
                break;
            }
            let feature = Self::build_genome_bed_feature(
                &record,
                &selected_track_name,
                path,
                local_start_0based,
                local_end_0based_exclusive,
                local_strand,
            );
            dna.features_mut().push(feature);
            report.imported_features += 1;
        }
        Self::append_chromosome_mismatch_warning(
            &mut report,
            &anchor.chromosome,
            "BED",
            &mismatch_counts,
        );
        if report.cancelled {
            report.warnings.push(format!(
                "BED import cancelled after parsed={}, imported={}, skipped={}",
                report.parsed_records, report.imported_features, report.skipped_records
            ));
        }
        if let Some(cb) = progress_cb.as_mut() {
            let _ = (**cb)(
                report.parsed_records,
                report.imported_features,
                report.skipped_records,
                true,
            );
        }

        Ok(report)
    }

    pub(super) fn resolve_bigwig_to_bedgraph_executable() -> String {
        crate::tool_overrides::resolve_tool_executable(
            BIGWIG_TO_BEDGRAPH_ENV_BIN,
            DEFAULT_BIGWIG_TO_BEDGRAPH_BIN,
        )
    }

    pub(super) fn convert_bigwig_to_bedgraph(path: &str) -> Result<NamedTempFile, EngineError> {
        let executable = Self::resolve_bigwig_to_bedgraph_executable();
        let output = NamedTempFile::new().map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!("Could not create temporary bedGraph file: {e}"),
        })?;
        let output_path = output.path().to_path_buf();
        let command_output = Command::new(&executable)
            .arg(path)
            .arg(&output_path)
            .output()
            .map_err(|e| EngineError {
                code: ErrorCode::Io,
                message: format!(
                    "Could not execute '{}' for BigWig conversion (set {} to override): {}",
                    executable, BIGWIG_TO_BEDGRAPH_ENV_BIN, e
                ),
            })?;
        if !command_output.status.success() {
            let stderr = String::from_utf8_lossy(&command_output.stderr)
                .trim()
                .to_string();
            let stdout = String::from_utf8_lossy(&command_output.stdout)
                .trim()
                .to_string();
            let detail = if !stderr.is_empty() {
                stderr
            } else if !stdout.is_empty() {
                stdout
            } else {
                format!("exit status {}", command_output.status)
            };
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "BigWig conversion failed for '{}' via '{}' (set {} to override): {}",
                    path, executable, BIGWIG_TO_BEDGRAPH_ENV_BIN, detail
                ),
            });
        }
        Ok(output)
    }

    pub(super) fn import_genome_bigwig_track(
        dna: &mut DNAsequence,
        anchor: &GenomeSequenceAnchor,
        path: &str,
        track_name: Option<&str>,
        min_score: Option<f64>,
        max_score: Option<f64>,
        clear_existing: bool,
        mut progress_cb: Option<&mut dyn FnMut(usize, usize, usize, bool) -> bool>,
    ) -> Result<GenomeBedTrackImportReport, EngineError> {
        let selected_track_name = track_name
            .map(str::trim)
            .filter(|v| !v.is_empty())
            .map(|v| v.to_string())
            .unwrap_or_else(|| Self::default_track_name(path));

        if clear_existing {
            Self::remove_generated_genome_signal_features(dna.features_mut());
        }

        let mut report = GenomeBedTrackImportReport {
            track_name: selected_track_name.clone(),
            ..Default::default()
        };
        let progress_stride = 250usize;
        let bedgraph_file = Self::convert_bigwig_to_bedgraph(path)?;
        let bedgraph_path = bedgraph_file.path().to_string_lossy().to_string();
        let mut reader = Self::open_text_reader(&bedgraph_path)?;
        let mut line = String::new();
        let mut line_no = 0usize;
        let mut mismatch_counts: HashMap<String, usize> = HashMap::new();

        while {
            line.clear();
            reader.read_line(&mut line).map_err(|e| EngineError {
                code: ErrorCode::Io,
                message: format!("Could not read converted bedGraph for '{path}': {e}"),
            })? > 0
        } {
            line_no += 1;
            let trimmed = line.trim();
            if trimmed.is_empty() {
                continue;
            }
            if trimmed.starts_with('#')
                || trimmed.to_ascii_lowercase().starts_with("track ")
                || trimmed.to_ascii_lowercase().starts_with("browser ")
            {
                continue;
            }
            report.parsed_records += 1;
            if report.parsed_records % progress_stride == 0 {
                if let Some(cb) = progress_cb.as_mut() {
                    let should_continue = (**cb)(
                        report.parsed_records,
                        report.imported_features,
                        report.skipped_records,
                        false,
                    );
                    if !should_continue {
                        report.cancelled = true;
                        break;
                    }
                }
            }

            let record = match Self::parse_bedgraph_record(trimmed) {
                Ok(v) => v,
                Err(e) => {
                    report.skipped_records += 1;
                    report.skipped_invalid += 1;
                    if report.warnings.len() < 20 {
                        report
                            .warnings
                            .push(format!("bedGraph line {} skipped: {}", line_no, e));
                    }
                    continue;
                }
            };

            if !Self::chromosomes_match(&record.chromosome, &anchor.chromosome) {
                report.skipped_records += 1;
                report.skipped_wrong_chromosome += 1;
                *mismatch_counts
                    .entry(record.chromosome.clone())
                    .or_insert(0) += 1;
                continue;
            }

            let bed_start_1based = record.start_0based.saturating_add(1);
            let bed_end_1based = record.end_0based;
            if bed_end_1based < anchor.start_1based || bed_start_1based > anchor.end_1based {
                report.skipped_records += 1;
                report.skipped_non_overlap += 1;
                continue;
            }
            let overlap_start_1based = bed_start_1based.max(anchor.start_1based);
            let overlap_end_1based = bed_end_1based.min(anchor.end_1based);
            if overlap_end_1based < overlap_start_1based {
                report.skipped_records += 1;
                report.skipped_non_overlap += 1;
                continue;
            }

            if min_score.is_some() || max_score.is_some() {
                let Some(score) = record.score else {
                    report.skipped_records += 1;
                    report.skipped_missing_score += 1;
                    continue;
                };
                if min_score.map(|v| score < v).unwrap_or(false)
                    || max_score.map(|v| score > v).unwrap_or(false)
                {
                    report.skipped_records += 1;
                    report.skipped_outside_score_range += 1;
                    continue;
                }
            }

            let (local_start_0based, local_end_0based_exclusive) = if anchor.strand == Some('-') {
                (
                    anchor.end_1based.saturating_sub(overlap_end_1based),
                    anchor.end_1based.saturating_sub(overlap_start_1based) + 1,
                )
            } else {
                (
                    overlap_start_1based - anchor.start_1based,
                    overlap_end_1based - anchor.start_1based + 1,
                )
            };
            if report.imported_features >= MAX_IMPORTED_SIGNAL_FEATURES {
                report.truncated_at_limit = true;
                break;
            }
            let feature = Self::build_genome_signal_feature(
                &record,
                &selected_track_name,
                path,
                local_start_0based,
                local_end_0based_exclusive,
                None,
                GENOME_BIGWIG_TRACK_GENERATED_TAG,
                "BigWig",
            );
            dna.features_mut().push(feature);
            report.imported_features += 1;
        }
        Self::append_chromosome_mismatch_warning(
            &mut report,
            &anchor.chromosome,
            "BigWig",
            &mismatch_counts,
        );
        if report.cancelled {
            report.warnings.push(format!(
                "BigWig import cancelled after parsed={}, imported={}, skipped={}",
                report.parsed_records, report.imported_features, report.skipped_records
            ));
        }
        if let Some(cb) = progress_cb.as_mut() {
            let _ = (**cb)(
                report.parsed_records,
                report.imported_features,
                report.skipped_records,
                true,
            );
        }

        Ok(report)
    }

    pub(super) fn import_genome_vcf_track(
        dna: &mut DNAsequence,
        anchor: &GenomeSequenceAnchor,
        path: &str,
        track_name: Option<&str>,
        min_score: Option<f64>,
        max_score: Option<f64>,
        clear_existing: bool,
        mut progress_cb: Option<&mut dyn FnMut(usize, usize, usize, bool) -> bool>,
    ) -> Result<GenomeBedTrackImportReport, EngineError> {
        let selected_track_name = track_name
            .map(str::trim)
            .filter(|v| !v.is_empty())
            .map(|v| v.to_string())
            .unwrap_or_else(|| Self::default_track_name(path));

        if clear_existing {
            Self::remove_generated_genome_signal_features(dna.features_mut());
        }

        let mut report = GenomeBedTrackImportReport {
            track_name: selected_track_name.clone(),
            ..Default::default()
        };
        let progress_stride = 250usize;
        let mut reader = Self::open_text_reader(path)?;
        let mut line = String::new();
        let mut line_no = 0usize;
        let mut sample_names: Vec<String> = vec![];
        let mut mismatch_counts: HashMap<String, usize> = HashMap::new();

        while {
            line.clear();
            reader.read_line(&mut line).map_err(|e| EngineError {
                code: ErrorCode::Io,
                message: format!("Could not read VCF file '{path}': {e}"),
            })? > 0
        } {
            line_no += 1;
            let trimmed = line.trim();
            if trimmed.is_empty() {
                continue;
            }
            if trimmed.starts_with("##") {
                continue;
            }
            if trimmed.starts_with("#CHROM") {
                let fields = trimmed.split('\t').collect::<Vec<_>>();
                sample_names = if fields.len() > 9 {
                    fields[9..]
                        .iter()
                        .map(|value| value.trim().to_string())
                        .collect::<Vec<_>>()
                } else {
                    vec![]
                };
                continue;
            }
            if trimmed.starts_with('#') {
                continue;
            }
            report.parsed_records += 1;
            if report.parsed_records % progress_stride == 0 {
                if let Some(cb) = progress_cb.as_mut() {
                    let should_continue = (**cb)(
                        report.parsed_records,
                        report.imported_features,
                        report.skipped_records,
                        false,
                    );
                    if !should_continue {
                        report.cancelled = true;
                        break;
                    }
                }
            }

            let record = match Self::parse_vcf_record(trimmed) {
                Ok(v) => v,
                Err(e) => {
                    report.skipped_records += 1;
                    report.skipped_invalid += 1;
                    if report.warnings.len() < 20 {
                        report
                            .warnings
                            .push(format!("VCF line {} skipped: {}", line_no, e));
                    }
                    continue;
                }
            };

            if !Self::chromosomes_match(&record.chromosome, &anchor.chromosome) {
                report.skipped_records += 1;
                report.skipped_wrong_chromosome += 1;
                *mismatch_counts
                    .entry(record.chromosome.clone())
                    .or_insert(0) += 1;
                continue;
            }

            let ref_len = record.reference.len().max(1);
            let variant_start_1based = record.pos_1based;
            let variant_end_1based = variant_start_1based.saturating_add(ref_len.saturating_sub(1));
            if variant_end_1based < anchor.start_1based || variant_start_1based > anchor.end_1based
            {
                report.skipped_records += 1;
                report.skipped_non_overlap += 1;
                continue;
            }

            if min_score.is_some() || max_score.is_some() {
                let Some(score) = record.qual else {
                    report.skipped_records += 1;
                    report.skipped_missing_score += 1;
                    continue;
                };
                if min_score.map(|v| score < v).unwrap_or(false)
                    || max_score.map(|v| score > v).unwrap_or(false)
                {
                    report.skipped_records += 1;
                    report.skipped_outside_score_range += 1;
                    continue;
                }
            }

            let overlap_start_1based = variant_start_1based.max(anchor.start_1based);
            let overlap_end_1based = variant_end_1based.min(anchor.end_1based);
            if overlap_end_1based < overlap_start_1based {
                report.skipped_records += 1;
                report.skipped_non_overlap += 1;
                continue;
            }

            let (local_start_0based, local_end_0based_exclusive) = if anchor.strand == Some('-') {
                (
                    anchor.end_1based.saturating_sub(overlap_end_1based),
                    anchor.end_1based.saturating_sub(overlap_start_1based) + 1,
                )
            } else {
                (
                    overlap_start_1based - anchor.start_1based,
                    overlap_end_1based - anchor.start_1based + 1,
                )
            };

            let mut stop = false;
            for (alt_idx, alt) in record.alternates.iter().enumerate() {
                if report.imported_features >= MAX_IMPORTED_SIGNAL_FEATURES {
                    report.truncated_at_limit = true;
                    stop = true;
                    break;
                }
                let feature = Self::build_genome_vcf_feature(
                    &record,
                    alt,
                    alt_idx + 1,
                    &sample_names,
                    &selected_track_name,
                    path,
                    local_start_0based,
                    local_end_0based_exclusive,
                );
                dna.features_mut().push(feature);
                report.imported_features += 1;
            }
            if stop {
                break;
            }
        }
        Self::append_chromosome_mismatch_warning(
            &mut report,
            &anchor.chromosome,
            "VCF",
            &mismatch_counts,
        );
        if report.cancelled {
            report.warnings.push(format!(
                "VCF import cancelled after parsed={}, imported={}, skipped={}",
                report.parsed_records, report.imported_features, report.skipped_records
            ));
        }
        if let Some(cb) = progress_cb.as_mut() {
            let _ = (**cb)(
                report.parsed_records,
                report.imported_features,
                report.skipped_records,
                true,
            );
        }

        Ok(report)
    }
}
