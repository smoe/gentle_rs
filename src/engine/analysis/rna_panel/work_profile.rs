//! Bounded deterministic observations of the existing index/nomination path; no scoring policy.

use super::*;

pub(super) const ANCHOR_UPPER_BOUNDS: [u64; 20] = [
    0, 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072,
    262144,
];

#[derive(Default)]
pub(super) struct ReadWork {
    pub anchors: u64,
    pub template_chain_evaluations: usize,
}

pub(super) fn index_profile(panel: &RnaReadPanel, index: &Index) -> RnaReadPanelIndexWorkProfile {
    let mut result = RnaReadPanelIndexWorkProfile {
        distinct_words: index.len(),
        max_postings_per_word: 0,
        max_templates_per_word: 0,
        max_loci_per_word: 0,
        words_above_posting_thresholds: [1, 8, 64, 256].into_iter().map(|n| (n, 0)).collect(),
    };
    for postings in index.values() {
        let templates = postings.iter().map(|p| p.0).collect::<BTreeSet<_>>();
        let loci = templates
            .iter()
            .map(|i| panel.templates[*i as usize].locus_id.as_str())
            .collect::<BTreeSet<_>>();
        result.max_postings_per_word = result.max_postings_per_word.max(postings.len());
        result.max_templates_per_word = result.max_templates_per_word.max(templates.len());
        result.max_loci_per_word = result.max_loci_per_word.max(loci.len());
        for (threshold, words) in &mut result.words_above_posting_thresholds {
            *words += usize::from(postings.len() > *threshold);
        }
    }
    result
}

pub(super) fn file_profile() -> RnaReadPanelFileWorkProfile {
    RnaReadPanelFileWorkProfile {
        peak_anchor_postings: None,
        anchor_histogram_upper_bounds: ANCHOR_UPPER_BOUNDS.to_vec(),
        anchor_histogram_counts: vec![0; ANCHOR_UPPER_BOUNDS.len()],
        total_anchor_postings: Some(0),
        peak_template_chain_evaluations: 0,
        equivalence_key_bytes_added: 0,
    }
}

pub(super) fn record_read(
    profile: &mut RnaReadPanelFileWorkProfile,
    work: &ReadWork,
    record_index: u64,
    read_id: &str,
) -> Result<(), EngineError> {
    let bucket = ANCHOR_UPPER_BOUNDS
        .iter()
        .position(|bound| work.anchors <= *bound)
        .ok_or_else(|| bad("panel_work_profile_anchor_range"))?;
    profile.anchor_histogram_counts[bucket] = profile.anchor_histogram_counts[bucket]
        .checked_add(1)
        .ok_or_else(|| bad("panel_work_profile_count_overflow"))?;
    profile.total_anchor_postings = profile
        .total_anchor_postings
        .and_then(|n| n.checked_add(work.anchors));
    profile.peak_template_chain_evaluations = profile
        .peak_template_chain_evaluations
        .max(work.template_chain_evaluations);
    // Equal maxima retain the first encountered record, not a template/hash iteration order.
    if profile
        .peak_anchor_postings
        .as_ref()
        .is_none_or(|peak| work.anchors > peak.anchors)
    {
        let prefix = read_id.chars().take(256).collect::<String>();
        let truncated = prefix.len() != read_id.len();
        profile.peak_anchor_postings = Some(RnaReadPanelReadWorkPeak {
            anchors: work.anchors,
            record_index,
            read_id: prefix,
            read_id_truncated: truncated,
            read_id_sha256: truncated.then(|| sha256_prefixed_str(read_id)),
        });
    }
    Ok(())
}

pub(super) fn validate(report: &RnaReadPanelScreenReport) -> Result<(), EngineError> {
    if let Some(index) = &report.index_work_profile {
        if index.distinct_words == 0
            || index.distinct_words > report.indexed_postings
            || index.max_postings_per_word == 0
            || index.max_postings_per_word > report.indexed_postings
            || index.max_templates_per_word == 0
            || index.max_templates_per_word > index.max_postings_per_word
            || index.max_loci_per_word == 0
            || index.max_loci_per_word > index.max_templates_per_word
            || index
                .words_above_posting_thresholds
                .keys()
                .copied()
                .collect::<Vec<_>>()
                != [1, 8, 64, 256]
            || index
                .words_above_posting_thresholds
                .values()
                .any(|n| *n > index.distinct_words)
        {
            return Err(bad("panel_index_work_profile_binding"));
        }
    }
    for file in &report.files_completed {
        if let Some(profile) = &file.work_profile {
            let count = profile
                .anchor_histogram_counts
                .iter()
                .try_fold(0u64, |sum, n| sum.checked_add(*n));
            if profile.anchor_histogram_upper_bounds != ANCHOR_UPPER_BOUNDS
                || profile.anchor_histogram_counts.len() != ANCHOR_UPPER_BOUNDS.len()
                || count != Some(file.reads)
                || profile.equivalence_key_bytes_added
                    != file
                        .equivalence_classes
                        .keys()
                        .map(String::len)
                        .sum::<usize>()
                || match &profile.peak_anchor_postings {
                    None => file.reads != 0 || profile.total_anchor_postings != Some(0),
                    Some(peak) => {
                        peak.record_index >= file.reads
                            || peak.anchors > report.request.policy.max_anchors_per_read as u64
                            || profile
                                .total_anchor_postings
                                .is_some_and(|total| peak.anchors > total)
                            || peak.read_id.chars().count() > 256
                            || peak.read_id_truncated != peak.read_id_sha256.is_some()
                            || peak.read_id_sha256.as_ref().is_some_and(|hash| {
                                !hash.starts_with("sha256:") || !valid_digest(&hash[7..], 64)
                            })
                            || ANCHOR_UPPER_BOUNDS.iter().position(|n| peak.anchors <= *n)
                                != profile.anchor_histogram_counts.iter().rposition(|n| *n > 0)
                    }
                }
            {
                return Err(bad("panel_file_work_profile_binding"));
            }
        }
    }
    Ok(())
}
