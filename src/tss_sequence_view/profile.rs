//! Read-only native projection of the same validated report used by TSS exports.
//! Local score arrays and sparse package hits retain separate units and provenance.

use super::*;
use gentle_protocol::{tss_motif_evidence, tss_profiles::*};
use std::{io::Read, path::Path};

/// Full motif-window-start arrays, including unavailable windows (not zero).
#[derive(Clone, Debug, Serialize)]
pub struct TssViewTrace {
    pub motif_length_bp: usize,
    pub forward: Vec<Option<f64>>,
    pub reverse: Vec<Option<f64>>,
    pub clip_negative: bool,
    pub range_is_fallback: bool,
}

/// Attachment identity is distinct from receipt verification or reference authenticity.
#[derive(Clone, Debug, Serialize)]
pub struct TssProfileAttachment {
    pub file_sha256: String,
    pub producer_revision: String,
    pub panel_id: String,
    pub warnings: Vec<String>,
}

impl TssSequenceView {
    /// Bounded file loading belongs on a worker, never in the paint callback.
    pub fn load_profile(&self, path: &Path) -> Result<Self, String> {
        const LIMIT: u64 = 256 * 1024 * 1024;
        let metadata = std::fs::metadata(path).map_err(|e| e.to_string())?;
        if !metadata.is_file() || metadata.len() > LIMIT {
            return Err("Choose a regular TSS report.json file no larger than 256 MiB".into());
        }
        let mut bytes = Vec::new();
        std::fs::File::open(path)
            .map_err(|e| e.to_string())?
            .take(LIMIT + 1)
            .read_to_end(&mut bytes)
            .map_err(|e| e.to_string())?;
        if bytes.len() as u64 > LIMIT {
            return Err("TSS report exceeds 256 MiB".into());
        }
        let report: TssProfileReport =
            serde_json::from_slice(&bytes).map_err(|e| format!("Not a TSS profile report: {e}"))?;
        let mut view = self.with_profile(&report)?;
        view.profile.as_mut().unwrap().file_sha256 = sha256_hex_bytes(&bytes);
        Ok(view)
    }

    /// Validate before attaching; equal bases alone do not establish locus identity.
    /// No scoring, database access, or changes to the annotated DNA are performed.
    pub fn with_profile(&self, report: &TssProfileReport) -> Result<Self, String> {
        crate::tss_profile_export::validate_tss_profile_report(report)
            .map_err(|e| e.to_string())?;
        if self.genome_id.as_deref() != Some(report.reference.genome_id.as_str())
            || self.assembly != report.reference.assembly
            || self.annotation_release != report.reference.annotation_release
        {
            return Err("TSS report reference/assembly/annotation release does not match the opened sequence; no aliases or liftover are inferred".into());
        }
        let matches = report
            .windows
            .iter()
            .filter(|w| w.record.promoter_id == self.promoter_id)
            .collect::<Vec<_>>();
        let [window] = matches.as_slice() else {
            return Err("TSS report must contain exactly one matching promoter ID".into());
        };
        if window.record.geometry != self.geometry
            || window.record.sequence_sha256 != self.sequence_sha256
        {
            return Err("TSS report sequence digest or genomic/TSS geometry does not match the opened sequence".into());
        }
        let mut view = self.clone();
        view.clear_profile();
        let panel = &report.panel_resolution.panel;
        let shared_range =
            (panel.scale_mode == TssScaleMode::Shared).then(|| trace_range(report, window, ""));
        for track in &window.tracks {
            let matrix = report
                .panel_resolution
                .matrices
                .iter()
                .find(|m| m.specification.source_id == track.accession)
                .ok_or("TSS track has no exact resolved matrix")?;
            let (scale_min, scale_max, range_is_fallback) =
                shared_range.unwrap_or_else(|| trace_range(report, window, &track.accession));
            let valid = track
                .forward_scores
                .iter()
                .chain(&track.reverse_scores)
                .filter(|v| v.is_some())
                .count();
            view.lanes.push(TssViewLane {
                kind: TssLaneKind::ScoreTrace,
                id: format!("profile/{}/{}/{}", panel.panel_id, track.accession, panel.score_kind),
                label: format!("{} | {}", matrix.specification.label, track.accession),
                details: format!("Saved full score arrays; matrix SHA-256 {}; panel SHA-256 {}; producer {}; normalization {}. Every x is a motif-window START on the displayed DNA, including reverse motifs. No new scoring. {}", matrix.matrix_sha256, report.panel_resolution.panel_sha256, report.producer_revision, track.normalization_reference, report.non_claims),
                units: panel.score_kind.clone(),
                state: format!("{valid}/{} valid strand-windows; {}", track.forward_scores.len() + track.reverse_scores.len(), if panel.clip_negative { "negative scores clipped for display" } else { "raw signed scores" }),
                features: vec![],
                trace: Some(TssViewTrace {
                    motif_length_bp: track.motif_length_bp,
                    forward: track.forward_scores.clone(),
                    reverse: track.reverse_scores.clone(),
                    clip_negative: panel.clip_negative,
                    range_is_fallback,
                }),
                scale_min,
                scale_max,
            });
        }
        for source in &report.imported_motif_evidence {
            let evidence = &source.report;
            if !evidence.regions.iter().any(|r| {
                r.resolved_chromosome.as_deref() == Some(self.geometry.chromosome.as_str())
                    && r.start_0based < self.geometry.end_1based
                    && r.end_0based_exclusive >= self.geometry.start_1based
            }) {
                continue;
            }
            let provider = evidence.provider.as_ref().ok_or("Missing motif provider")?;
            let subset_scope = evidence
                .regulatory_subset
                .as_ref()
                .map(|s| s.summary())
                .unwrap_or_default();
            for coverage in &evidence.motif_coverage {
                let hits =
                    tss_motif_evidence::project(evidence, &self.geometry, &coverage.motif_id);
                let (scale_min, scale_max) =
                    tss_motif_evidence::score_range(evidence, &coverage.motif_id);
                let features = hits.into_iter().map(|h| TssViewFeature {
                    feature_id: None,
                    start: h.start,
                    end: h.end,
                    reverse: h.local_strand == TssStrand::Minus,
                    clipped: h.clipped,
                    label: format!("{} | imported score {}", h.hit.motif_id, h.hit.score),
                    details: format!("Original footprint {}:{}..{}; local {} / genomic {}; raw score {} [{}]; source report {}; canonical report SHA-256 {}; source file SHA-256 {}. {}", h.hit.chromosome, h.hit.start_0based + 1, h.hit.end_0based_exclusive, h.local_strand.as_str(), h.hit.strand, h.hit.score, h.hit.score_mode, evidence.report_id, source.report_sha256, source.source.sha256, if h.clipped { "Footprint clipped to this TSS window." } else { "Full footprint." }),
                    score: Some(h.hit.score),
                }).collect();
                view.lanes.push(TssViewLane {
                    kind: TssLaneKind::ImportedMotif,
                    id: format!("imported/{}/{}", evidence.report_id, coverage.motif_id),
                    label: format!("{} | DuckDB {}", coverage.motif_id, if evidence.regulatory_subset.is_some() { "regulatory/TSS subset" } else { "retained hits" }),
                    details: format!("Provider {}; run {}; manifest SHA-256 {}; source floor {:?}; density limited {:?}; query complete {}; truncated {}. Scale fixed per source/report/matrix, NOT calibrated to local PWM scores. Missing hits are not evidence of absence. {} {}", provider.provider_kind, provider.run_id, provider.manifest_sha256, coverage.source_minimum_score, coverage.density_limited, evidence.query_complete, evidence.truncated, subset_scope, evidence.warnings.join("; ")),
                    units: provider.score_mode.clone(),
                    state: format!("{}; {}; {}", coverage.status.as_str(), if tss_motif_evidence::covers_window(evidence, &self.geometry) { "whole window queried" } else { "PARTIAL query coverage" }, if evidence.truncated { "TRUNCATED" } else if !evidence.query_complete { "query incomplete" } else { "query complete" }),
                    features,
                    trace: None,
                    scale_min,
                    scale_max,
                });
            }
        }
        let mut warnings = report.warnings.clone();
        if report.imported_motif_evidence.is_empty() {
            warnings.push("No imported DuckDB evidence attached; not evidence of absence.".into());
        } else if !view
            .lanes
            .iter()
            .any(|l| l.kind == TssLaneKind::ImportedMotif)
        {
            warnings.push(
                "Imported evidence has no query overlapping this exact chromosome/window.".into(),
            );
        }
        view.profile = Some(TssProfileAttachment {
            file_sha256: String::new(),
            producer_revision: report.producer_revision.clone(),
            panel_id: panel.panel_id.clone(),
            warnings,
        });
        Ok(view)
    }

    /// Drop only report-owned lanes; the original annotation view stays intact.
    pub fn clear_profile(&mut self) {
        self.lanes
            .retain(|l| !matches!(l.kind, TssLaneKind::ScoreTrace | TssLaneKind::ImportedMotif));
        self.profile = None;
    }
}

fn trace_range(
    report: &TssProfileReport,
    window: &TssProfileWindow,
    accession: &str,
) -> (f64, f64, bool) {
    let panel = &report.panel_resolution.panel;
    let windows = if panel.scale_mode == TssScaleMode::SharedAcrossTss {
        report.windows.as_slice()
    } else {
        std::slice::from_ref(window)
    };
    let mut min = 0.0_f64;
    let mut max = 0.0_f64;
    for score in windows
        .iter()
        .flat_map(|w| &w.tracks)
        .filter(|t| panel.scale_mode == TssScaleMode::Shared || t.accession == accession)
        .flat_map(|t| t.forward_scores.iter().chain(&t.reverse_scores))
        .flatten()
    {
        if !panel.clip_negative {
            min = min.min(*score);
        }
        max = max.max(*score);
    }
    let fallback = min == max;
    if fallback {
        max = min + 1.0;
    }
    (min, max, fallback)
}

#[cfg(test)]
pub(crate) mod tests;
