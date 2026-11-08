//! Separate basecall-composition availability and declared-unit-bound tail imports.
//! No SEQ/FASTQ join, BAM tags, raw-signal algorithm, caller conversion or clip reinterpretation.

use super::*;

#[cfg(test)]
#[path = "evidence_tests.rs"]
mod tests;

pub const COMPOSITION_SCHEMA: &str = "gentle.rna_read_end_terminal_composition.v1";
pub const TAIL_SCHEMA: &str = "gentle.rna_read_end_tail_estimate_import.v1";

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(tag = "status", rename_all = "snake_case", deny_unknown_fields)]
pub enum Provenance {
    Undeclared,
    Declared { description: String },
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CompositionPolicy {
    pub orientation: Provenance,
    pub trimming: Provenance,
    pub adapter_sequences: Vec<String>,
    pub primer_sequences: Vec<String>,
    pub window_bp: usize,
    pub minimum_run_bp: usize,
    pub maximum_mismatches: usize,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CompositionEvidence {
    pub schema: String,
    pub policy: CompositionPolicy,
    pub availability: String,
    pub reason: String,
    pub diagnostics: Vec<String>,
    pub nonclaim: String,
}

/// OD2 default: validate/retain the requested policy, never fabricate a scan.
pub fn unavailable_composition(
    policy: CompositionPolicy,
) -> Result<CompositionEvidence, EngineError> {
    if policy.window_bp == 0
        || policy.minimum_run_bp == 0
        || policy.minimum_run_bp > policy.window_bp
        || policy.maximum_mismatches >= policy.minimum_run_bp
    {
        return Err(fail("composition_parameters"));
    }
    for provenance in [&policy.orientation, &policy.trimming] {
        if let Provenance::Declared { description } = provenance
            && (!identity(description) || description.trim().is_empty())
        {
            return Err(fail("composition_provenance"));
        }
    }
    for sequence in policy
        .adapter_sequences
        .iter()
        .chain(&policy.primer_sequences)
    {
        if sequence.is_empty()
            || sequence.len() > 4096
            || !sequence
                .bytes()
                .all(|b| b"ACGTUNRYKMSWBDHVacgtunrykmswbdhv".contains(&b))
        {
            return Err(fail("composition_sequence"));
        }
    }
    let mut diagnostics = vec![];
    if matches!(policy.orientation, Provenance::Undeclared) {
        diagnostics.push("orientation_undeclared".into());
    }
    if matches!(policy.trimming, Provenance::Undeclared) {
        diagnostics.push("trimming_undeclared".into());
    }
    Ok(CompositionEvidence {
        schema: COMPOSITION_SCHEMA.into(),
        policy,
        availability: "unavailable".into(),
        reason: "no_bound_basecalls".into(),
        diagnostics,
        nonclaim:
            "Terminal composition is uncalibrated presence evidence, never physical tail length."
                .into(),
    })
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CallerConfiguration {
    pub caller: String,
    pub version: String,
    pub model: String,
    pub kit_chemistry: String,
    pub molecule_type: String,
    pub basecaller: String,
    pub basecaller_version: String,
    pub basecaller_model: String,
    pub parameters: serde_json::Value,
    pub parameter_sha256: String,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum QcStatus {
    Passed,
    Failed,
    Unavailable,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct TailRow {
    pub unit: ObservedSamQueryUnitKey,
    pub read_id: String,
    pub tail_estimate_nt: Option<f64>,
    pub uncertainty_nt: Option<f64>,
    pub quality_score: Option<f64>,
    pub qc_status: QcStatus,
    pub qc_details: String,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct TailImport {
    pub schema: String,
    pub raw_signal_manifest_sha256: String,
    pub configuration: CallerConfiguration,
    pub rows: Vec<TailRow>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum TailAvailability {
    Available,
    RejectedQc,
    Unsupported,
    Unavailable,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct TailEvidence {
    pub row: TailRow,
    pub availability: TailAvailability,
    /// Only supported, QC-passing values enter plots/statistics.
    pub accepted_estimate_nt: Option<f64>,
    pub raw_signal_manifest_sha256: String,
    pub configuration: CallerConfiguration,
}

fn configuration(value: &CallerConfiguration) -> Result<(), EngineError> {
    for text in [
        &value.caller,
        &value.version,
        &value.model,
        &value.kit_chemistry,
        &value.molecule_type,
        &value.basecaller,
        &value.basecaller_version,
        &value.basecaller_model,
    ] {
        if !identity(text) || text.trim().is_empty() {
            return Err(fail("tail_configuration"));
        }
    }
    let parameters = crate::digest_utils::canonical_json_string(&value.parameters)
        .map_err(|e| fail(&e.to_string()))?;
    if !value.parameters.is_object()
        || parameters.len() > 65536
        || crate::digest_utils::sha256_prefixed_str(&parameters) != value.parameter_sha256
    {
        return Err(fail("tail_parameter_binding"));
    }
    Ok(())
}

/// Validate an already file-authenticated import against explicit configuration and unit witnesses.
/// Raw-source hashes bind a declared manifest; this does not read/verify POD5/FAST5 bytes.
pub fn import_tail_estimates(
    input: &TailImport,
    supported: &[CallerConfiguration],
    raw_signal_manifest_sha256: &str,
    analysis: &ObservedSamAnalysis,
    max_rows: usize,
) -> Result<Vec<TailEvidence>, EngineError> {
    if input.schema != TAIL_SCHEMA
        || !sha(raw_signal_manifest_sha256)
        || input.raw_signal_manifest_sha256 != raw_signal_manifest_sha256
        || input.rows.len() > max_rows
    {
        return Err(fail("tail_source_binding"));
    }
    configuration(&input.configuration)?;
    for value in supported {
        configuration(value)?;
    }
    let supported = supported.contains(&input.configuration);
    let units: BTreeMap<_, _> = analysis.units.iter().map(|u| (&u.unit, u)).collect();
    if units.len() != analysis.units.len() {
        return Err(fail("tail_unit_identity"));
    }
    let mut seen = BTreeSet::new();
    let mut out = Vec::new();
    for row in &input.rows {
        if !seen.insert(&row.unit) || !identity(&row.read_id) || !identity(&row.qc_details) {
            return Err(fail("tail_unit_identity"));
        }
        let unit = units
            .get(&row.unit)
            .ok_or_else(|| fail("tail_unit_binding"))?;
        if unit.rows.is_empty() || unit.rows.iter().any(|r| r.qname != row.read_id) {
            return Err(fail("tail_read_binding"));
        }
        if row
            .tail_estimate_nt
            .is_some_and(|v| !v.is_finite() || v < 0.0)
            || row
                .uncertainty_nt
                .is_some_and(|v| !v.is_finite() || v < 0.0)
            || row.quality_score.is_some_and(|v| !v.is_finite())
            || (row.qc_status == QcStatus::Passed && row.tail_estimate_nt.is_none())
        {
            return Err(fail("tail_value"));
        }
        let availability = if !supported {
            TailAvailability::Unsupported
        } else {
            match row.qc_status {
                QcStatus::Passed => TailAvailability::Available,
                QcStatus::Failed => TailAvailability::RejectedQc,
                QcStatus::Unavailable => TailAvailability::Unavailable,
            }
        };
        out.push(TailEvidence {
            row: row.clone(),
            availability,
            accepted_estimate_nt: (availability == TailAvailability::Available)
                .then_some(row.tail_estimate_nt)
                .flatten(),
            raw_signal_manifest_sha256: input.raw_signal_manifest_sha256.clone(),
            configuration: input.configuration.clone(),
        });
    }
    Ok(out)
}
