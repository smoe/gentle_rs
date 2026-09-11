-- Synthetic only. Created by the test harness, never loaded as package code.
CREATE TABLE anchors AS SELECT 1::BIGINT anchor_id, '1' chrom,
100::BIGINT anchor_start,116::BIGINT anchor_end,0.5::DOUBLE anchor_score,
true supported_tp73_saos2_TA,7.0::DOUBLE depth_tp73_saos2_TA
UNION ALL SELECT 2,'2',100,116,-0.2,false,0.0;
CREATE TABLE feature AS SELECT 1::BIGINT anchor_id,'MA9000.1' motif_id,
'gap_6_20' distance_band,133::BIGINT hit_start,148::BIGINT hit_end,
4.25::DOUBLE best_score,4.25::DOUBLE plus_score,NULL::DOUBLE minus_score,
'+' best_strand,17::BIGINT interval_distance_bp,'right' genomic_side,
1::BIGINT n_source_loci,1::BIGINT n_score_zero_loci;
CREATE TABLE anchor_promoter AS SELECT 1::BIGINT anchor_id,'promoter-a' regulatory_feature_id
UNION ALL SELECT 1,'promoter-b';
CREATE TABLE promoter AS SELECT 'promoter-a' regulatory_feature_id,'1' chrom,
50::BIGINT extended_start,200::BIGINT extended_end,'GRCh38' assembly
UNION ALL SELECT 'promoter-b','1',60,220,'GRCh38';
CREATE TABLE promoter_gene AS SELECT 'promoter-a' regulatory_feature_id,'GENE-A' gene_id,
'ensembl_native' link_source,'synthetic-v1' annotation_release
UNION ALL SELECT 'promoter-a','GENE-B','tss_derived','synthetic-v1'
UNION ALL SELECT 'promoter-b','GENE-A','ensembl_native','synthetic-v1';
CREATE TABLE cofactor_distance_isoform_comparison AS
SELECT 'MA9000.1' motif_id,'SYNTHETIC-A' motif_name,'gap_6_20' distance_band,
2::INTEGER distance_band_order,2::BIGINT anchors_total,1::BIGINT anchors_positive,
0.5::DOUBLE positive_anchor_fraction,2.5::DOUBLE ta_adjusted_odds_ratio,
0.8::DOUBLE dn_adjusted_odds_ratio,1.1::DOUBLE ta_confidence_interval_95_lower,
4.0::DOUBLE ta_confidence_interval_95_upper,0.2::DOUBLE dn_confidence_interval_95_lower,
1.8::DOUBLE dn_confidence_interval_95_upper,0.04::DOUBLE ta_q_value_bh_tax_group,
0.7::DOUBLE dn_q_value_bh_tax_group,3.125::DOUBLE ta_vs_dn_odds_ratio_ratio,
1.4::DOUBLE confidence_interval_95_lower,5.0::DOUBLE confidence_interval_95_upper,
0.03::DOUBLE q_value_bh_tax_group,1.14::DOUBLE ta_vs_dn_log_odds_difference,
'ok' ta_evaluation_status,'ok' dn_evaluation_status,'ok' evaluation_status,
'Synthetic species' source_species,'synthetic_support' class_support_flag;
INSERT INTO cofactor_distance_isoform_comparison SELECT * REPLACE (
'MA9001.1' AS motif_id,'SYNTHETIC-B' AS motif_name,1.5 AS ta_adjusted_odds_ratio)
FROM cofactor_distance_isoform_comparison;
