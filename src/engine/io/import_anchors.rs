//! GenBank region parsing, inferred anchor metadata, and import-origin helpers.

use super::*;

impl GentleEngine {
    pub(super) fn is_genbank_like_path(path: &str) -> bool {
        let lower = path.to_ascii_lowercase();
        lower.ends_with(".gb")
            || lower.ends_with(".gbk")
            || lower.ends_with(".genbank")
            || lower.ends_with(".gbff")
    }

    pub(super) fn parse_first_usize_tokens(raw: &str, max_items: usize) -> Vec<usize> {
        if max_items == 0 {
            return vec![];
        }
        let mut out = Vec::with_capacity(max_items);
        let mut current = String::new();
        for ch in raw.chars() {
            if ch.is_ascii_digit() || ch == ',' {
                current.push(ch);
                continue;
            }
            if !current.is_empty() {
                if let Ok(value) = current.replace(',', "").parse::<usize>() {
                    out.push(value);
                    if out.len() >= max_items {
                        return out;
                    }
                }
                current.clear();
            }
        }
        if !current.is_empty() && out.len() < max_items {
            if let Ok(value) = current.replace(',', "").parse::<usize>() {
                out.push(value);
            }
        }
        out
    }

    pub(super) fn parse_genbank_header_field_block(text: &str, key: &str) -> Option<String> {
        let mut block = String::new();
        let mut in_block = false;
        for line in text.lines() {
            if line.starts_with("ORIGIN") || line.starts_with("//") {
                break;
            }
            if let Some(rest) = line.strip_prefix(key) {
                in_block = true;
                let trimmed = rest.trim();
                if !trimmed.is_empty() {
                    block.push_str(trimmed);
                }
                continue;
            }
            if in_block {
                if line.starts_with("            ") {
                    let trimmed = line.trim();
                    if !trimmed.is_empty() {
                        if !block.is_empty() {
                            block.push(' ');
                        }
                        block.push_str(trimmed);
                    }
                    continue;
                }
                break;
            }
        }
        if block.is_empty() { None } else { Some(block) }
    }

    pub(super) fn parse_genbank_accession_region(
        path: &str,
    ) -> Option<(String, usize, usize, String, char)> {
        if !Self::is_genbank_like_path(path) {
            return None;
        }
        let text = std::fs::read_to_string(path).ok()?;
        let definition =
            Self::parse_genbank_header_field_block(&text, "DEFINITION").unwrap_or_default();
        let accession_block = Self::parse_genbank_header_field_block(&text, "ACCESSION")?;
        if accession_block.is_empty() {
            return None;
        }
        let accession = accession_block.split_whitespace().next()?.to_string();
        let lower_block = accession_block.to_ascii_lowercase();
        let region_pos = lower_block.find("region:")?;
        let region_spec = &accession_block[region_pos + "region:".len()..];
        let lower_region_spec = region_spec.to_ascii_lowercase();
        let anchor_strand = if lower_region_spec.contains("complement(") {
            '-'
        } else {
            '+'
        };
        if lower_region_spec.contains("join(") || lower_region_spec.contains("order(") {
            return None;
        }
        let numbers = Self::parse_first_usize_tokens(region_spec, 2);
        if numbers.len() < 2 {
            return None;
        }
        let mut start_1based = numbers[0];
        let mut end_1based = numbers[1];
        if start_1based == 0 || end_1based == 0 {
            return None;
        }
        if end_1based < start_1based {
            std::mem::swap(&mut start_1based, &mut end_1based);
        }
        Some((
            accession,
            start_1based,
            end_1based,
            definition,
            anchor_strand,
        ))
    }

    pub(super) fn parse_chromosome_from_definition(definition: &str) -> Option<String> {
        let lower = definition.to_ascii_lowercase();
        let marker = "chromosome ";
        let marker_pos = lower.find(marker)?;
        let raw_tail = &definition[marker_pos + marker.len()..];
        let token = raw_tail
            .chars()
            .take_while(|c| c.is_ascii_alphanumeric() || matches!(c, '_' | '-' | '.'))
            .collect::<String>();
        if token.is_empty() { None } else { Some(token) }
    }

    pub(super) fn parse_genome_id_from_definition(definition: &str) -> Option<String> {
        let lower = definition.to_ascii_lowercase();
        let marker_pos = lower.find("primary assembly")?;
        let prefix = definition[..marker_pos]
            .trim()
            .trim_end_matches(|c: char| c.is_ascii_whitespace() || c == ',' || c == '.');
        let candidate = prefix.rsplit(',').next()?.trim();
        if candidate.is_empty() {
            None
        } else {
            Some(candidate.to_string())
        }
    }

    pub(super) fn infer_imported_genbank_anchor(
        path: &str,
        dna: &DNAsequence,
    ) -> Option<GenomeSequenceAnchor> {
        let (accession, start_1based, end_1based, definition, anchor_strand) =
            Self::parse_genbank_accession_region(path)?;
        let region_len = end_1based - start_1based + 1;
        if dna.len() > 0 && dna.len() != region_len {
            return None;
        }
        let chromosome = Self::parse_chromosome_from_definition(&definition)
            .unwrap_or_else(|| accession.clone());
        let genome_id =
            Self::parse_genome_id_from_definition(&definition).unwrap_or_else(|| accession.clone());
        Some(GenomeSequenceAnchor {
            genome_id,
            chromosome,
            start_1based,
            end_1based,
            strand: Some(anchor_strand),
            anchor_verified: None,
            catalog_path: None,
            cache_dir: None,
        })
    }

    pub(super) fn verify_anchor_sequence_against_catalog(
        dna: &DNAsequence,
        anchor: &GenomeSequenceAnchor,
        catalog_path: &str,
        cache_dir: Option<&str>,
    ) -> Result<bool, String> {
        let catalog = GenomeCatalog::from_json_file(catalog_path)?;
        let mut reference = catalog.get_sequence_region_with_cache(
            &anchor.genome_id,
            &anchor.chromosome,
            anchor.start_1based,
            anchor.end_1based,
            cache_dir,
        )?;
        if anchor.strand == Some('-') {
            reference = Self::reverse_complement(&reference);
        }
        Ok(reference.eq_ignore_ascii_case(&dna.get_forward_string()))
    }

    pub(super) fn genome_anchor_fallback_policy_for_extension(
        &self,
    ) -> PreparedGenomeFallbackPolicy {
        match self.state.parameters.genome_anchor_prepared_fallback_policy {
            GenomeAnchorPreparedFallbackPolicy::Off => PreparedGenomeFallbackPolicy::Off,
            GenomeAnchorPreparedFallbackPolicy::SingleCompatible => {
                PreparedGenomeFallbackPolicy::SingleCompatible
            }
            GenomeAnchorPreparedFallbackPolicy::AlwaysExplicit => {
                PreparedGenomeFallbackPolicy::AlwaysExplicit
            }
        }
    }

    pub(super) fn classify_import_origin(path: &str, dna: &DNAsequence) -> SequenceOrigin {
        let lower = path.to_ascii_lowercase();
        if lower.ends_with(".fa")
            || lower.ends_with(".fasta")
            || lower.ends_with(".fna")
            || lower.ends_with(".ffn")
            || lower.ends_with(".faa")
        {
            // Requested policy: treat all FASTA imports as synthetic for now.
            return SequenceOrigin::ImportedSynthetic;
        }

        // For GenBank/EMBL-like records, use SOURCE/mol_type metadata if present.
        for feature in dna.features() {
            if feature.kind.to_string().to_ascii_uppercase() != "SOURCE" {
                continue;
            }
            for key in ["mol_type", "molecule_type"] {
                if let Some(value) = feature.qualifier_values(key.into()).next() {
                    let v = value.to_ascii_lowercase();
                    if v.contains("synthetic") {
                        return SequenceOrigin::ImportedSynthetic;
                    }
                    if v.contains("cdna") || v.contains("mrna") || v.contains("transcript") {
                        return SequenceOrigin::ImportedCdna;
                    }
                    if v.contains("genomic") {
                        return SequenceOrigin::ImportedGenomic;
                    }
                }
            }
        }

        SequenceOrigin::ImportedUnknown
    }
}
