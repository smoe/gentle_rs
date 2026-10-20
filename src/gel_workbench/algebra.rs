//! Bounded, memoized replay of the frozen algebra-v1 algorithm and digest.

use super::*;
use gentle_protocol::ErrorCode;
use std::{
    collections::{BTreeMap, BTreeSet},
    sync::Arc,
};

const MAX_DEPTH: usize = 16;
const MAX_INPUTS: usize = 64;
const MAX_REPLAY_BINS: usize = 20_000 * 64;

#[cfg(test)]
mod tests;

#[derive(serde::Serialize)]
struct TermV1<'a> {
    profile_id: &'a str,
    coefficient: f64,
}

#[derive(serde::Serialize)]
struct RequestV1<'a> {
    output_id: &'a str,
    terms: Vec<TermV1<'a>>,
    range: [f64; 2],
}

fn request_v1(request: &GelAlgebraRequest) -> RequestV1<'_> {
    RequestV1 {
        output_id: &request.output_id,
        terms: request
            .terms
            .iter()
            .map(|t| TermV1 {
                profile_id: &t.profile_id,
                coefficient: t.coefficient,
            })
            .collect(),
        range: request.range,
    }
}

pub(super) fn request_json_v1(request: &GelAlgebraRequest) -> Result<String, EngineError> {
    serde_json::to_string(&request_v1(request)).map_err(|e| EngineError::internal(e.to_string()))
}

// Freeze field order/presence independently of future GelProfile extensions.
#[derive(serde::Serialize)]
struct ProfileV1<'a> {
    id: &'a str,
    label: &'a str,
    run_id: &'a str,
    lane_id: &'a str,
    source_sha256: &'a str,
    evidence_kind: GelEvidenceKind,
    size_kind: GelSizeKind,
    signal_units: &'a str,
    edges: &'a [f64],
    weights: &'a [Option<f64>],
    measure: &'a str,
    processing: &'a [String],
    warnings: &'a [String],
}

pub(super) fn input_digest_v1(
    profiles: &[&GelProfile],
    request: &GelAlgebraRequest,
) -> Result<String, EngineError> {
    let profiles: Vec<_> = profiles
        .iter()
        .map(|p| ProfileV1 {
            id: &p.id,
            label: &p.label,
            run_id: &p.run_id,
            lane_id: &p.lane_id,
            source_sha256: &p.source_sha256,
            evidence_kind: p.evidence_kind,
            size_kind: p.size_kind,
            signal_units: &p.signal_units,
            edges: &p.edges,
            weights: &p.weights,
            measure: &p.measure,
            processing: &p.processing,
            warnings: &p.warnings,
        })
        .collect();
    Ok(sha256_hex_bytes(
        &serde_json::to_vec(&(profiles, request_v1(request)))
            .map_err(|e| EngineError::internal(e.to_string()))?,
    ))
}

#[derive(Clone)]
struct Failure(GelAlgebraVerification, String);

fn unavailable(message: impl Into<String>) -> Failure {
    Failure(GelAlgebraVerification::Unverifiable, message.into())
}

fn mismatch(message: impl Into<String>) -> Failure {
    Failure(GelAlgebraVerification::Mismatch, message.into())
}

struct Replay<'a> {
    profiles: &'a BTreeMap<String, Arc<GelProfile>>,
    path: BTreeSet<String>,
    checked: BTreeMap<String, Result<(), Failure>>,
    bins: usize,
}

impl Replay<'_> {
    fn check(&mut self, id: &str, depth: usize) -> Result<(), Failure> {
        if depth > MAX_DEPTH {
            return Err(unavailable("Algebra replay exceeds depth 16"));
        }
        if self.path.contains(id) {
            return Err(mismatch("Algebra dependency cycle"));
        }
        if let Some(result) = self.checked.get(id) {
            return result.clone();
        }
        if self.checked.len() + self.path.len() >= MAX_INPUTS {
            return Err(unavailable("Algebra replay exceeds 64 profiles"));
        }
        let profile = self
            .profiles
            .get(id)
            .ok_or_else(|| unavailable(format!("Missing algebra operand '{id}'")))?
            .clone();
        if profile.id != id {
            return Err(mismatch("Profile map key and recorded identity differ"));
        }
        validate_profile(&profile).map_err(|e| mismatch(e.message))?;
        self.path.insert(id.into());
        let result = self.check_profile(&profile, depth);
        self.path.remove(id);
        self.checked.insert(id.into(), result.clone());
        result
    }

    fn check_profile(&mut self, profile: &GelProfile, depth: usize) -> Result<(), Failure> {
        self.bins = self.bins.saturating_add(profile.weights.len());
        if self.bins > MAX_REPLAY_BINS {
            return Err(unavailable("Algebra replay exceeds 1280000 bin visits"));
        }
        if profile.evidence_kind != GelEvidenceKind::SignedResidual {
            return Ok(());
        }
        let request = receipt(profile)?;
        let mut operands = Vec::new();
        for term in &request.terms {
            self.check(&term.profile_id, depth + 1)?;
            let operand = self.profiles.get(&term.profile_id).unwrap().clone();
            self.bins = self.bins.saturating_add(operand.weights.len());
            if self.bins > MAX_REPLAY_BINS {
                return Err(unavailable("Algebra replay exceeds 1280000 bin visits"));
            }
            operands.push(operand);
        }
        let refs: Vec<_> = operands.iter().map(Arc::as_ref).collect();
        let expected = super::algebra(&refs, &request).map_err(|e| mismatch(e.message))?;
        if expected != *profile {
            return Err(mismatch(format!(
                "Algebra v1 replay differs from saved residual '{}'",
                profile.id
            )));
        }
        Ok(())
    }
}

fn receipt(profile: &GelProfile) -> Result<GelAlgebraRequest, Failure> {
    if profile.processing.len() != 1 {
        return Err(unavailable("No unique algebra-v1 receipt"));
    }
    let request: GelAlgebraRequest =
        serde_json::from_str(&profile.processing[0]).map_err(|_| {
            unavailable("No readable algebra-v1 receipt; imported/legacy data remains inspectable")
        })?;
    if request.output_id != profile.id {
        return Err(mismatch("Algebra receipt output identity differs"));
    }
    if request.terms.is_empty() || request.terms.len() > 16 {
        return Err(mismatch("Algebra receipt requires 1..16 terms"));
    }
    Ok(request)
}

/// Construct a presentation from saved data only; never re-extract, normalize
/// or fit sources. Unavailable/mismatching results cannot acquire plot data.
pub fn inspect_algebra(
    store: &GelWorkbenchStore,
    id: &str,
) -> Result<GelAlgebraPresentation, EngineError> {
    let profile = store
        .profiles
        .get(id)
        .ok_or_else(|| EngineError::new(ErrorCode::NotFound, "Saved algebra profile not found"))?;
    let mut presentation = GelAlgebraPresentation {
        schema: "gentle.gel_algebra_presentation.v1".into(),
        profile_id: id.into(),
        label: profile.label.clone(),
        verification: GelAlgebraVerification::Unverifiable,
        reason: None,
        input_digest: None,
        request: None,
        series: None,
        contributions: vec![],
        inputs: vec![],
        warnings: vec![],
    };
    if profile.evidence_kind != GelEvidenceKind::SignedResidual {
        presentation.reason = Some("Selected profile is not a saved algebra result".into());
        return Ok(presentation);
    }
    let mut replay = Replay {
        profiles: &store.profiles,
        path: BTreeSet::new(),
        checked: BTreeMap::new(),
        bins: 0,
    };
    let verified = replay.check(id, 0);
    if let Err(Failure(status, reason)) = verified {
        presentation.verification = status;
        presentation.reason = Some(reason);
        return Ok(presentation);
    }
    let request = receipt(profile).map_err(|e| EngineError::internal(e.1))?;
    let density: Vec<_> = profile
        .weights
        .iter()
        .zip(profile.edges.windows(2))
        .map(|(w, edges)| w.map(|w| w / (edges[1] - edges[0])))
        .collect();
    if density.iter().flatten().any(|x| !x.is_finite()) {
        presentation.reason =
            Some("Internally consistent result exceeds finite density display range".into());
        return Ok(presentation);
    }
    let mut contributions = Vec::new();
    for (term_index, term) in request.terms.iter().enumerate() {
        if term.coefficient == 0.0 {
            continue;
        }
        let weights: Vec<_> = super::rebin(&store.profiles[&term.profile_id], &profile.edges)
            .into_iter()
            .map(|w| w.map(|w| w * term.coefficient))
            .collect();
        let density: Vec<_> = weights
            .iter()
            .zip(profile.edges.windows(2))
            .map(|(w, edges)| w.map(|w| w / (edges[1] - edges[0])))
            .collect();
        if weights
            .iter()
            .chain(&density)
            .flatten()
            .any(|x| !x.is_finite())
        {
            presentation.reason =
                Some("Internally consistent result exceeds finite operand display range".into());
            return Ok(presentation);
        }
        contributions.push(GelAlgebraContribution {
            term_index,
            profile_id: term.profile_id.clone(),
            coefficient: term.coefficient,
            weights,
            density,
        });
    }
    let (lo, hi) = density
        .iter()
        .chain(contributions.iter().flat_map(|c| c.density.iter()))
        .flatten()
        .fold((0.0_f64, 0.0_f64), |(lo, hi), v| (lo.min(*v), hi.max(*v)));
    if !(hi - lo).is_finite() {
        presentation.reason =
            Some("Internally consistent result exceeds finite density-axis range".into());
        return Ok(presentation);
    }
    let mut warnings = BTreeSet::new();
    warnings.insert("Verified means algebra-v1 internal consistency, not source authenticity or scientific acceptance. No normalization, fitting, distribution summaries or p-values.".to_string());
    warnings.insert("Every operand retains its own calibration/processing. Different ladders or registration uncertainty require review; subtraction does not align them.".into());
    warnings.insert("Operand overlays show each nonzero coefficient times its rebinned integral per size unit, on one shared axis. Each retains its own gaps; a visible operand does not fill a missing residual. No fitted or manual display scale.".into());
    if request.terms.iter().any(|t| t.coefficient == 0.0) {
        warnings.insert("Zero-coefficient operands bind provenance and grid edges but contribute neither signal nor missing-value exclusions.".into());
    }
    for key in replay.checked.keys().filter(|key| key.as_str() != id) {
        let input = store.profiles[key].as_ref().clone();
        warnings.extend(input.warnings.iter().cloned());
        presentation.inputs.push(input);
    }
    warnings.extend(profile.warnings.iter().cloned());
    if density.iter().all(Option::is_none) {
        warnings.insert(
            "No evaluated bins on the requested support; unavailable is not a measured zero."
                .into(),
        );
    }
    presentation.verification = GelAlgebraVerification::Verified;
    presentation.input_digest = Some(profile.source_sha256.clone());
    presentation.request = Some(request);
    presentation.contributions = contributions;
    presentation.series = Some(GelAlgebraSeries {
        size_kind: profile.size_kind,
        signal_units: profile.signal_units.clone(),
        measure: profile.measure.clone(),
        edges: profile.edges.clone(),
        weights: profile.weights.clone(),
        density,
    });
    presentation.warnings = warnings.into_iter().collect();
    Ok(presentation)
}
