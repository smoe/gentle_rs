//! Candidate-first, read-only presentation of the engine's motif tuple ledger.
//!
//! Joins use exact report IDs, never gene symbols. Distances, ordering and
//! eligibility are copied from the ledger, not recalculated by a frontend.

use gentle_protocol::{
    GeneSetCutRunEvaluationState, RegulatoryPartnerGeneRow, RegulatoryPartnerMotifHit,
    RegulatoryPartnerMotifRole, RegulatoryPartnerScreenReport, RegulatoryPartnerTupleRow,
};
use std::collections::{BTreeMap, BTreeSet};

/// Shared GUI/clipboard columns. Every row is a motif pair, not a unique factor
/// or an independent CUT&RUN peak; expression is not assessed by this screen.
pub const REGULATORY_PARTNER_CANDIDATE_COLUMNS: [&str; 10] = [
    "Gene",
    "Predicted anchor motif",
    "Candidate motif",
    "Anchor site (genomic, 1-based)",
    "Candidate site (genomic, 1-based)",
    "Centre distance (bp)",
    "LLR bits (anchor / candidate)",
    "Distance rule",
    "Promoter CUT&RUN",
    "Expression",
];

/// References to one exact engine-owned tuple and its source rows.
pub struct RegulatoryPartnerCandidateRow<'a> {
    pub gene: &'a RegulatoryPartnerGeneRow,
    pub anchor: &'a RegulatoryPartnerMotifHit,
    pub partner: &'a RegulatoryPartnerMotifHit,
    pub tuple: &'a RegulatoryPartnerTupleRow,
}

impl RegulatoryPartnerCandidateRow<'_> {
    /// Plain table cells, suitable for the GUI and tab-separated clipboard.
    pub fn cells(&self) -> [String; 10] {
        let motif = |hit: &RegulatoryPartnerMotifHit| match hit.tf_name.as_deref() {
            Some(name) if !name.is_empty() && name != hit.tf_id => {
                format!("{name} ({})", hit.tf_id)
            }
            _ => hit.tf_id.clone(),
        };
        let site = |hit: &RegulatoryPartnerMotifHit| {
            format!(
                "{}:{}-{} ({})",
                self.gene.chromosome.as_deref().unwrap_or("unknown"),
                hit.genomic_start_1based,
                hit.genomic_end_1based,
                if hit.genomic_forward_strand { "+" } else { "-" }
            )
        };
        let occupancy = match self.tuple.occupancy_evaluation_state {
            GeneSetCutRunEvaluationState::Unevaluated => "Not evaluated",
            GeneSetCutRunEvaluationState::Evaluated if self.tuple.promoter_occupancy_supported => {
                "Promoter supported"
            }
            GeneSetCutRunEvaluationState::Evaluated => "Evaluated: no support",
        };
        [
            self.gene.gene_symbol.clone(),
            motif(self.anchor),
            motif(self.partner),
            site(self.anchor),
            site(self.partner),
            format!(
                "{:+.1}",
                self.tuple.signed_anchor_to_partner_center_distance_bp
            ),
            format!("{} / {}", self.anchor.llr_bits, self.partner.llr_bits),
            if self.tuple.within_requested_distance {
                "Within"
            } else {
                "Outside"
            }
            .into(),
            occupancy.into(),
            "Not evaluated".into(),
        ]
    }
}

/// Join ledger IDs once, preserving engine order and the engine's proximity
/// decision. Reject ambiguous or cross-member references instead of hiding
/// missing evidence, including in rows excluded by the display filter.
pub fn regulatory_partner_candidate_rows(
    report: &RegulatoryPartnerScreenReport,
    nearby_only: bool,
) -> Result<Vec<RegulatoryPartnerCandidateRow<'_>>, String> {
    let mut genes = BTreeMap::new();
    for gene in &report.ledger.genes {
        if genes.insert(gene.member_dedup_key.as_str(), gene).is_some() {
            return Err(format!("Duplicate member key '{}'", gene.member_dedup_key));
        }
    }
    let mut hits = BTreeMap::new();
    for hit in &report.ledger.motif_hits {
        if hits.insert(hit.hit_id.as_str(), hit).is_some() {
            return Err(format!("Duplicate motif hit ID '{}'", hit.hit_id));
        }
    }
    let mut tuple_ids = BTreeSet::new();
    let mut rows = Vec::new();
    for tuple in &report.ledger.tuples {
        let gene = genes.get(tuple.member_dedup_key.as_str());
        let anchor = hits.get(tuple.anchor_hit_id.as_str());
        let partner = hits.get(tuple.partner_hit_id.as_str());
        let (Some(gene), Some(anchor), Some(partner)) = (gene, anchor, partner) else {
            return Err(format!(
                "Missing gene or motif evidence for tuple '{}'",
                tuple.tuple_id
            ));
        };
        if !tuple_ids.insert(&tuple.tuple_id)
            || anchor.member_dedup_key != gene.member_dedup_key
            || partner.member_dedup_key != gene.member_dedup_key
            || anchor.role != RegulatoryPartnerMotifRole::Anchor
            || partner.role != RegulatoryPartnerMotifRole::Partner
        {
            return Err(format!(
                "Ambiguous gene or motif evidence for tuple '{}'",
                tuple.tuple_id
            ));
        }
        if !nearby_only || tuple.within_requested_distance {
            rows.push(RegulatoryPartnerCandidateRow {
                gene,
                anchor,
                partner,
                tuple,
            });
        }
    }
    Ok(rows)
}

/// Export the displayed pairs with exact join IDs and input-report digests.
/// This is a presentation of the portable report, not a new analysis or rank.
pub fn regulatory_partner_candidates_tsv(
    report: &RegulatoryPartnerScreenReport,
    nearby_only: bool,
) -> Result<String, String> {
    let rows = regulatory_partner_candidate_rows(report, nearby_only)?;
    let mut columns = REGULATORY_PARTNER_CANDIDATE_COLUMNS.to_vec();
    columns.extend([
        "tuple_id",
        "member_key",
        "transcript_id",
        "genome_id",
        "max_distance_bp",
        "op_id",
        "run_id",
        "source_report_sha256",
    ]);
    let mut output = columns.join("\t");
    output.push('\n');
    let sources = report
        .ledger
        .source_reports
        .iter()
        .map(|source| source.content_sha256.as_str())
        .collect::<Vec<_>>()
        .join(";");
    for row in rows {
        let mut cells = row.cells().to_vec();
        cells.extend([
            row.tuple.tuple_id.clone(),
            row.gene.member_dedup_key.clone(),
            row.gene.transcript_id.clone().unwrap_or_default(),
            report.ledger.request.genome_id.clone(),
            report
                .ledger
                .request
                .max_anchor_partner_distance_bp
                .to_string(),
            report.op_id.clone().unwrap_or_default(),
            report.run_id.clone().unwrap_or_default(),
            sources.clone(),
        ]);
        output.push_str(
            &cells
                .iter()
                .map(|cell| cell.replace(['\t', '\r', '\n'], " "))
                .collect::<Vec<_>>()
                .join("\t"),
        );
        output.push('\n');
    }
    Ok(output)
}

#[cfg(test)]
mod tests {
    use super::*;

    // Hand-crafted presentation fixture, not a biological scan. Deliberately
    // reverse-oriented coordinates prove that rendering preserves the ledger.
    fn report() -> RegulatoryPartnerScreenReport {
        let mut report = RegulatoryPartnerScreenReport::default();
        report.ledger.genes.push(RegulatoryPartnerGeneRow {
            member_dedup_key: "gene:one".into(),
            gene_symbol: "ONE".into(),
            chromosome: Some("chr7".into()),
            strand: Some("-".into()),
            transcript_id: Some("tx1".into()),
            ..Default::default()
        });
        report.ledger.motif_hits = vec![
            RegulatoryPartnerMotifHit {
                hit_id: "anchor".into(),
                member_dedup_key: "gene:one".into(),
                role: RegulatoryPartnerMotifRole::Anchor,
                tf_id: "anchor.1".into(),
                genomic_start_1based: 901,
                genomic_end_1based: 910,
                ..Default::default()
            },
            RegulatoryPartnerMotifHit {
                hit_id: "partner".into(),
                member_dedup_key: "gene:one".into(),
                role: RegulatoryPartnerMotifRole::Partner,
                tf_id: "partner.2".into(),
                tf_name: Some("Candidate".into()),
                genomic_start_1based: 801,
                genomic_end_1based: 810,
                llr_bits: 4.25,
                ..Default::default()
            },
        ];
        report.ledger.tuples.push(RegulatoryPartnerTupleRow {
            tuple_id: "pair1".into(),
            member_dedup_key: "gene:one".into(),
            anchor_hit_id: "anchor".into(),
            partner_hit_id: "partner".into(),
            signed_anchor_to_partner_center_distance_bp: 100.0,
            within_requested_distance: true,
            ..Default::default()
        });
        report
    }

    #[test]
    fn candidate_table_preserves_motif_ids_strand_distances_and_unknown_evidence() {
        let report = report();
        let rows = regulatory_partner_candidate_rows(&report, true).unwrap();
        let cells = rows[0].cells();
        assert_eq!(cells[2], "Candidate (partner.2)");
        assert_eq!(cells[3], "chr7:901-910 (-)");
        assert_eq!(cells[4], "chr7:801-810 (-)");
        assert_eq!(cells[5], "+100.0");
        assert_eq!(cells[6], "0 / 4.25");
        assert_eq!(&cells[8..], ["Not evaluated", "Not evaluated"]);
    }

    #[test]
    fn candidate_filter_uses_engine_outcome_without_reranking_or_merging() {
        let mut report = report();
        let mut other = report.ledger.tuples[0].clone();
        other.tuple_id = "pair2".into();
        other.within_requested_distance = false;
        report.ledger.tuples.insert(0, other);
        let all = regulatory_partner_candidate_rows(&report, false).unwrap();
        assert_eq!(
            all.iter()
                .map(|row| row.tuple.tuple_id.as_str())
                .collect::<Vec<_>>(),
            ["pair2", "pair1"]
        );
        assert_eq!(
            regulatory_partner_candidate_rows(&report, true)
                .unwrap()
                .len(),
            1
        );
        report.ledger.tuples[1].within_requested_distance = false;
        assert!(
            regulatory_partner_candidate_rows(&report, true)
                .unwrap()
                .is_empty()
        );
        assert_eq!(
            regulatory_partner_candidate_rows(&report, false)
                .unwrap()
                .len(),
            2
        );
    }

    #[test]
    fn candidate_table_distinguishes_zero_support_from_unknown() {
        let mut report = report();
        report.ledger.tuples[0].occupancy_evaluation_state =
            GeneSetCutRunEvaluationState::Evaluated;
        assert_eq!(
            regulatory_partner_candidate_rows(&report, true).unwrap()[0].cells()[8],
            "Evaluated: no support"
        );
        report.ledger.tuples[0].promoter_occupancy_supported = true;
        assert_eq!(
            regulatory_partner_candidate_rows(&report, true).unwrap()[0].cells()[8],
            "Promoter supported"
        );
    }

    #[test]
    fn candidate_table_rejects_missing_cross_member_and_duplicate_hit_references() {
        let mut missing = report();
        missing.ledger.motif_hits.pop();
        missing.ledger.tuples[0].within_requested_distance = false;
        assert!(regulatory_partner_candidate_rows(&missing, true).is_err());
        let mut crossed = report();
        crossed.ledger.motif_hits[1].member_dedup_key = "another_gene".into();
        assert!(regulatory_partner_candidate_rows(&crossed, false).is_err());
        let mut duplicate = report();
        duplicate
            .ledger
            .motif_hits
            .push(duplicate.ledger.motif_hits[0].clone());
        assert!(regulatory_partner_candidate_rows(&duplicate, false).is_err());
    }

    #[test]
    fn candidate_tsv_is_deterministic_retains_join_ids_and_cannot_inject_rows() {
        let mut report = report();
        report.op_id = Some("screen:one".into());
        report.ledger.request.genome_id = "synthetic_genome".into();
        report
            .ledger
            .source_reports
            .push(gentle_protocol::RegulatoryPartnerSourceReportRef {
                content_sha256: format!("sha256:{}", "a".repeat(64)),
                ..Default::default()
            });
        report.ledger.genes[0].gene_symbol = "Gene\twith\nwhitespace".into();
        let tsv = regulatory_partner_candidates_tsv(&report, true).unwrap();
        assert_eq!(
            tsv,
            regulatory_partner_candidates_tsv(&report, true).unwrap()
        );
        let lines = tsv.lines().collect::<Vec<_>>();
        assert_eq!(lines.len(), 2);
        assert_eq!(lines[0].split('\t').count(), lines[1].split('\t').count());
        assert!(lines[1].starts_with("Gene with whitespace\t"));
        assert!(lines[1].contains("\tpair1\tgene:one\ttx1\t"));
        assert!(lines[1].contains("\tsynthetic_genome\t0\tscreen:one\t"));
        assert!(lines[1].ends_with(&format!("sha256:{}", "a".repeat(64))));
        report.ledger.tuples.clear();
        assert_eq!(
            regulatory_partner_candidates_tsv(&report, false)
                .unwrap()
                .lines()
                .count(),
            1
        );
    }
}
