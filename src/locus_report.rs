//! Portable locus documents retain their schema and verify live sequence geometry.

use crate::{
    dna_sequence::DNAsequence,
    engine::{PromoterReporterArchitectureComparisonReport, SequenceGenomeAnchorSummary},
};
use gentle_protocol::{
    GENE_LOCUS_EVIDENCE_DISPLAY_SCHEMA, GeneLocusEvidenceDisplayReport,
    GeneLocusGenomeAnchorBinding, GeneLocusSequenceBinding,
};
use std::{io::Read, path::Path};

/// Bounded file loading for potentially large score-vector reports.
pub fn read_document(path: &Path) -> Result<LocusDocument, String> {
    const MAX_BYTES: u64 = 128 * 1024 * 1024;
    let file = std::fs::File::open(path).map_err(|error| error.to_string())?;
    let mut bytes = Vec::new();
    file.take(MAX_BYTES + 1)
        .read_to_end(&mut bytes)
        .map_err(|error| error.to_string())?;
    if bytes.len() as u64 > MAX_BYTES {
        return Err("Locus report exceeds the 128 MiB import limit".into());
    }
    LocusDocument::from_json(&bytes)
}

/// Both publication envelopes supported by the locus inspector.
#[derive(Debug, Clone)]
pub struct LocusDocument(LocusEnvelope);

#[derive(Debug, Clone)]
enum LocusEnvelope {
    Locus(Box<GeneLocusEvidenceDisplayReport>),
    Reporter(Box<PromoterReporterArchitectureComparisonReport>),
}

impl LocusDocument {
    /// Dispatch by schema: both records have serde defaults and share field names.
    pub fn from_json(bytes: &[u8]) -> Result<Self, String> {
        let value: serde_json::Value = serde_json::from_slice(bytes)
            .map_err(|error| format!("Invalid locus report JSON: {error}"))?;
        let document = match value.get("schema").and_then(|value| value.as_str()) {
            Some(GENE_LOCUS_EVIDENCE_DISPLAY_SCHEMA) => serde_json::from_value(value)
                .map(LocusEnvelope::Locus)
                .map_err(|error| error.to_string())?,
            Some("gentle.promoter_reporter_architecture_comparison.v1") => {
                let report: Box<PromoterReporterArchitectureComparisonReport> =
                    serde_json::from_value(value).map_err(|error| error.to_string())?;
                let locus = report.locus_evidence.as_ref().ok_or(
                    "Reporter report has no composed locus evidence; recompose with locus evidence",
                )?;
                if locus.seq_id != report.seq_id {
                    return Err("Reporter and nested locus sequence identities disagree".into());
                }
                if let Some(binding) = &locus.sequence_binding
                    && (binding
                        .sequence_sha256
                        .strip_prefix("sha256:")
                        .unwrap_or(&binding.sequence_sha256)
                        != report
                            .source_provenance
                            .source_sequence_sha256
                            .strip_prefix("sha256:")
                            .unwrap_or(&report.source_provenance.source_sequence_sha256)
                        || binding.genome_anchor
                            != report
                                .source_provenance
                                .genome_anchor
                                .as_ref()
                                .map(anchor_binding))
                {
                    return Err(
                        "Reporter and nested locus sequence/anchor bindings disagree".into(),
                    );
                }
                LocusEnvelope::Reporter(report)
            }
            _ => {
                return Err(
                    "Expected a GENtle locus-evidence or reporter-comparison schema".into(),
                );
            }
        };
        let document = Self(document);
        let locus = document.locus();
        if locus.schema != GENE_LOCUS_EVIDENCE_DISPLAY_SCHEMA
            || locus.seq_id.is_empty()
            || locus.locus_local_start_1based == 0
            || locus.locus_local_end_1based < locus.locus_local_start_1based
        {
            return Err("Invalid locus report schema, sequence identity or local bounds".into());
        }
        Ok(document)
    }

    pub fn locus(&self) -> &GeneLocusEvidenceDisplayReport {
        match &self.0 {
            LocusEnvelope::Locus(report) => report,
            LocusEnvelope::Reporter(report) => {
                report.locus_evidence.as_ref().expect("validated envelope")
            }
        }
    }

    /// Preserve the original envelope, including reporter architectures, on export.
    pub fn to_json(&self) -> Result<Vec<u8>, String> {
        match &self.0 {
            LocusEnvelope::Locus(report) => serde_json::to_vec_pretty(report),
            LocusEnvelope::Reporter(report) => serde_json::to_vec_pretty(report),
        }
        .map_err(|error| error.to_string())
    }

    pub fn render_svg(&self) -> String {
        match &self.0 {
            LocusEnvelope::Locus(report) => gentle_render::render_gene_locus_evidence_with_overlay_svg(report, None),
            LocusEnvelope::Reporter(report) => crate::render_promoter_reporter_architecture::render_promoter_reporter_architecture_svg(report),
        }
    }

    /// Reuse the publication renderer's resolved reporter materials and legend.
    pub fn reporter_overlay(&self) -> Option<gentle_render::GeneLocusEvidenceOverlay> {
        match &self.0 {
            LocusEnvelope::Locus(_) => None,
            LocusEnvelope::Reporter(report) => {
                Some(crate::render_promoter_reporter_architecture::normalized_locus_overlay(report))
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::json;

    fn fixture() -> (GeneLocusEvidenceDisplayReport, GeneLocusSequenceBinding) {
        let dna = DNAsequence::from_sequence(&"ACGT".repeat(25)).unwrap();
        let anchor = SequenceGenomeAnchorSummary {
            seq_id: "demo".into(),
            genome_id: "GRCh38".into(),
            chromosome: "chr7".into(),
            start_1based: 1001,
            end_1based: 1100,
            strand: Some('-'),
            anchor_verified: Some(true),
        };
        let binding = sequence_binding(&dna, Some(&anchor));
        (
            GeneLocusEvidenceDisplayReport {
                schema: GENE_LOCUS_EVIDENCE_DISPLAY_SCHEMA.into(),
                seq_id: "demo".into(),
                sequence_binding: Some(binding.clone()),
                locus_local_start_1based: 1,
                locus_local_end_1based: 100,
                ..Default::default()
            },
            binding,
        )
    }

    #[test]
    fn live_binding_checks_bases_length_and_complete_anchor() {
        let (report, current) = fixture();
        assert!(verify_live_binding(&report, "demo", &current).is_ok());
        for mutation in 0..7 {
            let mut changed = current.clone();
            match mutation {
                0 => changed.sequence_sha256 = "same length, different bases".into(),
                1 => changed.sequence_length_bp += 1,
                2 => changed.genome_anchor.as_mut().unwrap().strand = Some('+'),
                3 => changed.genome_anchor.as_mut().unwrap().genome_id = "GRCh37".into(),
                4 => changed.genome_anchor.as_mut().unwrap().chromosome = "chr8".into(),
                5 => changed.genome_anchor.as_mut().unwrap().end_1based += 1,
                _ => changed.genome_anchor = None,
            }
            assert!(
                verify_live_binding(&report, "demo", &changed).is_err(),
                "mutation {mutation}"
            );
        }
        assert!(verify_live_binding(&report, "another", &current).is_err());
        let mut outside = report;
        outside.locus_local_end_1based = 101;
        assert!(verify_live_binding(&outside, "demo", &current).is_err());
    }

    #[test]
    fn legacy_locus_is_viewable_but_cannot_authorize_live_dna() {
        let (mut report, current) = fixture();
        report.sequence_binding = None;
        let bytes = serde_json::to_vec(&report).unwrap();
        assert!(!String::from_utf8_lossy(&bytes).contains("sequence_binding"));
        let loaded = LocusDocument::from_json(&bytes).unwrap();
        assert!(loaded.render_svg().contains("<svg"));
        assert!(
            verify_live_binding(loaded.locus(), "demo", &current)
                .unwrap_err()
                .contains("Historical")
        );
    }

    #[test]
    fn schema_dispatch_preserves_reporter_envelope_and_rejects_inconsistent_bindings() {
        let (locus, current) = fixture();
        let document = json!({
            "schema": "gentle.promoter_reporter_architecture_comparison.v1", "seq_id": "demo",
            "source_provenance": {
                "source_sequence_sha256": current.sequence_sha256.strip_prefix("sha256:").unwrap(),
                "genome_anchor": { "seq_id": "demo", "genome_id": "GRCh38", "chromosome": "chr7", "start_1based": 1001, "end_1based": 1100, "strand": "-" }
            },
            "locus_evidence": locus,
            "architectures": [{"architecture_id": "spliced_demo", "transcript_id": "demo_tx", "segments": [{"segment_id": "utr", "material": "spliced_cdna", "source_start_0based": 10, "source_end_0based_exclusive": 30}]}]
        });
        let loaded = LocusDocument::from_json(&serde_json::to_vec(&document).unwrap()).unwrap();
        let saved: serde_json::Value = serde_json::from_slice(&loaded.to_json().unwrap()).unwrap();
        assert_eq!(saved["schema"], document["schema"]);
        assert_eq!(saved["architectures"][0]["architecture_id"], "spliced_demo");
        let overlay = loaded.reporter_overlay().unwrap();
        assert_eq!(overlay.rows.len(), 1);
        assert_eq!(overlay.rows[0].segments[0].local_start_1based, 11);
        assert!(loaded.render_svg().contains("spliced_demo"));
        assert!(verify_live_binding(loaded.locus(), "demo", &current).is_ok());
        for (key, value) in [
            ("seq_id", json!("wrong")),
            ("source_provenance", json!({})),
            ("locus_evidence", json!(null)),
            ("schema", json!("unknown")),
        ] {
            let mut malformed = document.clone();
            malformed[key] = value;
            assert!(
                LocusDocument::from_json(&serde_json::to_vec(&malformed).unwrap()).is_err(),
                "{key}"
            );
        }
    }
}

fn anchor_binding(anchor: &SequenceGenomeAnchorSummary) -> GeneLocusGenomeAnchorBinding {
    GeneLocusGenomeAnchorBinding {
        genome_id: anchor.genome_id.clone(),
        chromosome: anchor.chromosome.clone(),
        start_1based: anchor.start_1based,
        end_1based: anchor.end_1based,
        strand: anchor.strand,
    }
}

/// Bind the loaded bases, length and coordinate authority once at composition.
pub fn sequence_binding(
    dna: &DNAsequence,
    anchor: Option<&SequenceGenomeAnchorSummary>,
) -> GeneLocusSequenceBinding {
    GeneLocusSequenceBinding {
        sequence_sha256: crate::digest_utils::sha256_prefixed_str(&dna.get_forward_string()),
        sequence_length_bp: dna.len(),
        genome_anchor: anchor.map(anchor_binding),
    }
}

/// Legacy or changed documents may be read, but cannot authorize live DNA actions.
pub fn verify_live_binding(
    report: &GeneLocusEvidenceDisplayReport,
    seq_id: &str,
    current: &GeneLocusSequenceBinding,
) -> Result<(), String> {
    if report.seq_id != seq_id {
        return Err("Report sequence identity does not match the active sequence".into());
    }
    let binding = report
        .sequence_binding
        .as_ref()
        .ok_or("Historical report has no sequence binding; recompose before selecting live DNA")?;
    if binding.sequence_length_bp != current.sequence_length_bp
        || binding.sequence_sha256 != current.sequence_sha256
    {
        return Err("Report sequence digest/length differs from the active DNA; recompose".into());
    }
    if binding.genome_anchor != current.genome_anchor {
        return Err(
            "Report genome anchor differs in assembly, contig, bounds or strand; recompose".into(),
        );
    }
    if report.locus_local_start_1based == 0
        || report.locus_local_start_1based > report.locus_local_end_1based
        || report.locus_local_end_1based > current.sequence_length_bp
    {
        return Err("Report locus bounds lie outside its source sequence".into());
    }
    Ok(())
}
