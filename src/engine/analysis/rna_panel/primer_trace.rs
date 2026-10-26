//! Independent primer trace publication using the existing FASTQ/hash/control primitives.

use super::*;
use serde_json::Value;
use std::cell::Cell;
use std::rc::Rc;

#[path = "primer_trace/matching.rs"]
mod matching;
const ALGORITHM: &str = "primer_trace_local_concordance_v1";
const BINDING_SCHEMA: &str = "gentle.rna_primer_trace_result_binding.v1";

fn digest<T: Serialize>(value: &T) -> Result<String, EngineError> {
    Ok(sha256_prefixed_bytes(&encode(value)?))
}
fn valid_hash(value: &str) -> bool {
    value
        .strip_prefix("sha256:")
        .is_some_and(|s| valid_digest(s, 64))
}

fn validate(panel: &RnaPrimerPanel) -> Result<Vec<String>, EngineError> {
    let mut basis = panel.clone();
    basis.content_sha256.clear();
    let p = &panel.policy;
    if panel.schema != RNA_PRIMER_PANEL_SCHEMA
        || digest(&basis)? != panel.content_sha256
        || panel.primers.is_empty()
        || panel.primers.len() > 128
        || p.max_full_edits > 3
        || p.max_core_edits_for_attribution > 3
        || !(64..=4096).contains(&p.end_window_nt)
        || !(4..=64).contains(&p.poly_t_min_run)
        || p.poly_t_max_mismatches > 2
        || p.poly_t_max_mismatches >= p.poly_t_min_run
        || !(10..=64).contains(&p.concordance_probe_nt)
        || p.max_probe_edits > 3
        || !(2..=16).contains(&p.concordance_seed_nt)
        || p.concordance_seed_nt
            > (p.concordance_probe_nt - p.max_probe_edits) / (p.max_probe_edits + 1)
        || !(p.concordance_probe_nt..=4096).contains(&p.min_payload_nt)
    {
        return Err(bad("primer_panel_binding_or_policy"));
    }
    let mut ids = BTreeSet::new();
    let mut sequences = BTreeSet::new();
    let mut cores = BTreeSet::new();
    let mut handles = BTreeMap::<String, usize>::new();
    for primer in &panel.primers {
        if !safe_id(&primer.primer_id)
            || !ids.insert(&primer.primer_id)
            || primer.assay_label.trim().is_empty()
            || primer.assay_label.len() > 512
            || primer.source.trim().is_empty()
            || primer.source.len() > 4096
            || primer.variant.len() > 128
            || !(10..=64).contains(&primer.sequence.len())
            || !primer.sequence.bytes().all(|b| b"ACGT".contains(&b))
            || primer
                .source_sha256
                .as_ref()
                .is_some_and(|h| !valid_hash(h))
            || primer.expected_targets.len() > 4096
            || primer.components.is_empty()
            || primer.components.len() > 16
        {
            return Err(bad("primer_panel_oligo"));
        }
        let key = |s: &[u8]| std::cmp::min(s.to_vec(), matching::rc(s));
        if !sequences.insert(key(primer.sequence.as_bytes())) {
            return Err(bad("primer_panel_duplicate_sequence"));
        }
        let mut components = primer.components.iter().collect::<Vec<_>>();
        components.sort_by_key(|c| c.start);
        let mut cursor = 0;
        let mut core_count = 0;
        let mut component_ids = BTreeSet::new();
        for c in components {
            if !safe_id(&c.component_id)
                || !component_ids.insert(&c.component_id)
                || c.start != cursor
                || c.end <= c.start
                || c.end > primer.sequence.len()
            {
                return Err(bad("primer_panel_component_tiling"));
            }
            cursor = c.end;
            match c.role.as_str() {
                "target_core" => {
                    core_count += 1;
                    if c.end != primer.sequence.len()
                        || c.end - c.start < 10
                        || !cores.insert(key(&primer.sequence.as_bytes()[c.start..c.end]))
                    {
                        return Err(bad("primer_panel_core_terminal_or_duplicate"));
                    }
                }
                "common_handle" => {
                    *handles
                        .entry(primer.sequence[c.start..c.end].into())
                        .or_default() += 1;
                }
                _ => return Err(bad("primer_panel_component_role")),
            }
        }
        if cursor != primer.sequence.len() || core_count != 1 {
            return Err(bad("primer_panel_component_tiling"));
        }
        let mut targets = BTreeSet::new();
        for t in &primer.expected_targets {
            if !safe_id(&t.locus_id) || !targets.insert((&t.locus_id, t.transcript_feature_id)) {
                return Err(bad("primer_panel_target_identity"));
            }
        }
    }
    Ok(handles
        .values()
        .any(|&n| n > 1)
        .then(|| "shared_handle_is_not_target_identity".into())
        .into_iter()
        .collect())
}

fn resolve(binding: &RnaPrimerPanelBinding) -> Result<RnaPrimerPanel, EngineError> {
    let panel = match binding {
        RnaPrimerPanelBinding::Inline { panel } => *panel.clone(),
        RnaPrimerPanelBinding::File {
            path,
            expected_sha256,
        } => {
            let panel: RnaPrimerPanel = load(path)?;
            if !valid_hash(expected_sha256) || &panel.content_sha256 != expected_sha256 {
                return Err(bad("primer_panel_expected_digest"));
            }
            panel
        }
    };
    validate(&panel)?;
    Ok(panel)
}

struct LimitedWriter {
    inner: HashedWriter,
    used: Rc<Cell<u64>>,
    limit: u64,
}
impl LimitedWriter {
    fn new(path: &Path, used: &Rc<Cell<u64>>, limit: u64) -> Result<Self, EngineError> {
        Ok(Self {
            inner: HashedWriter::create(path)?,
            used: used.clone(),
            limit,
        })
    }
    fn receipt(self, path: String) -> Result<RnaReadPanelArtifact, EngineError> {
        self.inner.receipt(path)
    }
}
impl Write for LimitedWriter {
    fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
        let total = self
            .used
            .get()
            .checked_add(bytes.len() as u64)
            .filter(|n| *n <= self.limit)
            .ok_or_else(|| std::io::Error::other("primer_trace_sidecar_budget"))?;
        let n = self.inner.write(bytes)?;
        self.used.set(total - (bytes.len() - n) as u64);
        Ok(n)
    }
    fn flush(&mut self) -> std::io::Result<()> {
        self.inner.flush()
    }
}

fn bucket(value: usize) -> String {
    if value <= 256 {
        return value.to_string();
    }
    let mut lower = 257;
    let mut upper = 512;
    while value > upper {
        lower = upper + 1;
        upper = (upper * 2).min(2_000_000);
    }
    format!("{lower}..{upper}")
}
pub(super) fn count(counts: &mut RnaPrimerTraceCounts, row: &RnaPrimerTraceRow, len: usize) {
    counts.reads += 1;
    *counts
        .anchor_classes
        .entry(row.anchor_class.clone())
        .or_default() += 1;
    *counts
        .locus_relations
        .entry(row.locus_relation.clone())
        .or_default() += 1;
    for hit in row
        .hits
        .iter()
        .filter(|h| h.anchor_class.starts_with("full_"))
    {
        let strand = if hit.strand == PrimerReadStrand::Forward {
            "forward"
        } else {
            "reverse_complement"
        };
        let key =
            serde_json::to_string(&(hit.primer_id.as_deref(), hit.anchor_class.as_str(), strand))
                .unwrap();
        *counts
            .full_hits_by_primer_class_strand
            .entry(key)
            .or_default() += 1;
        *counts
            .anchor_start_histogram
            .entry(bucket(hit.start))
            .or_default() += 1;
        *counts
            .remaining_tail_histogram
            .entry(bucket(len - hit.end))
            .or_default() += 1;
    }
}

/// Shared TSV decoding only; callers separately decide prefix inspection versus EOF verification.
pub(super) const HEADER: [&str; 21] = [
    "record_index",
    "read_id",
    "anchor_class",
    "anchor_position",
    "primer_multiplicity",
    "locus_relation",
    "evidence_tier",
    "anchor_json",
    "leader_json",
    "payload_json",
    "payload_transcript_sense",
    "payload_mean_phred",
    "payload_fraction_below_q10",
    "template_matches_json",
    "artifact_flags_json",
    "reasons_json",
    "hit_count",
    "hits_json",
    "read_length",
    "leader_length",
    "payload_length",
];

pub(super) fn parse_row(r: &csv::StringRecord) -> Result<RnaPrimerTraceRow, EngineError> {
    if r.len() != 21 {
        return Err(bad("primer_trace_inspection_columns"));
    }
    let parse = |i: usize| serde_json::from_str::<Value>(&r[i]).map_err(|e| bad(e.to_string()));
    let value = serde_json::json!({"record_index":r[0].parse::<u64>().map_err(|e|bad(e.to_string()))?,"read_id":r[1].to_string(),
        "anchor_class":r[2].to_string(),"anchor_position":if r[3].is_empty(){Value::Null}else{Value::String(r[3].to_string())},
        "primer_multiplicity":if r[4].is_empty(){Value::Null}else{Value::String(r[4].to_string())},"locus_relation":r[5].to_string(),
        "evidence_tier":r[6].to_string(),"anchor":parse(7)?,"leader":parse(8)?,"payload":parse(9)?,"payload_transcript_sense":parse(10)?,
        "payload_mean_phred":parse(11)?,"payload_fraction_below_q10":parse(12)?,"template_matches":parse(13)?,"artifact_flags":parse(14)?,"reasons":parse(15)?,"hits":parse(17)?,
        "read_length":r[18].parse::<usize>().map_err(|e|bad(e.to_string()))?,"leader_length":parse(19)?,"payload_length":parse(20)?});
    serde_json::from_value(value).map_err(|e| bad(e.to_string()))
}

fn validate_report(report: &RnaPrimerTraceReport) -> Result<(), EngineError> {
    let mut basis = report.clone();
    basis.content_sha256.clear();
    validate(&report.primer_panel)?;
    if report.schema != RNA_PRIMER_TRACE_SCHEMA
        || report.algorithm != ALGORITHM
        || digest(&basis)? != report.content_sha256
        || !matches!(
            (report.status.as_str(), report.control_error.is_some()),
            ("completed" | "cancelled", false) | ("control_error", true)
        )
        || report.request.schema != RNA_PRIMER_TRACE_REQUEST_SCHEMA
        || report.request.files.is_empty()
        || report.request.files.len() > 4096
        || !valid_hash(&report.submitted_request_sha256)
        || report.files_completed.len() + report.inputs_not_completed.len()
            != report.request.files.len()
    {
        return Err(bad("primer_trace_report_binding"));
    }
    let primer_hash = match &report.request.primer_panel {
        RnaPrimerPanelBinding::Inline { panel } => {
            validate(panel)?;
            &panel.content_sha256
        }
        RnaPrimerPanelBinding::File {
            expected_sha256, ..
        } => expected_sha256,
    };
    if primer_hash != &report.primer_panel.content_sha256
        || report.reference_panel_sha256.as_ref()
            != report
                .request
                .reference_panel
                .as_ref()
                .map(|p| &p.expected_sha256)
        || report
            .reference_panel_sha256
            .as_ref()
            .is_some_and(|h| !valid_hash(h))
        || (report.status == "completed" && !report.inputs_not_completed.is_empty())
        || encode(&report.inputs_not_completed)?
            != encode(&report.request.files[report.files_completed.len()..].to_vec())?
    {
        return Err(bad("primer_trace_report_inputs"));
    }
    let mut bytes = 0u64;
    for (i, file) in report.files_completed.iter().enumerate() {
        let id = format!("file_{i:04}");
        let sum =
            |map: &BTreeMap<String, u64>| map.values().try_fold(0u64, |a, b| a.checked_add(*b));
        if file.file_id != id
            || encode(&file.input)? != encode(&report.request.files[i])?
            || !valid_hash(&file.observed_sha256)
            || !valid_digest(&file.observed_md5, 32)
            || sum(&file.counts.anchor_classes) != Some(file.counts.reads)
            || sum(&file.counts.locus_relations) != Some(file.counts.reads)
            || sum(&file.counts.full_hits_by_primer_class_strand)
                != sum(&file.counts.anchor_start_histogram)
            || sum(&file.counts.full_hits_by_primer_class_strand).is_none()
            || sum(&file.counts.full_hits_by_primer_class_strand)
                != sum(&file.counts.remaining_tail_histogram)
            || file
                .input
                .expected_md5
                .as_ref()
                .is_some_and(|h| h != &file.observed_md5)
            || file.input.expected_sha256.as_ref().is_some_and(|h| {
                h.trim_start_matches("sha256:").to_ascii_lowercase()
                    != file.observed_sha256.trim_start_matches("sha256:")
            })
        {
            return Err(bad("primer_trace_file_binding"));
        }
        let allowed = if report.request.export_payload_fasta {
            vec![format!("{id}/reads.tsv"), format!("{id}/payload.fasta")]
        } else {
            vec![format!("{id}/reads.tsv")]
        };
        if file.artifacts.len() != allowed.len()
            || file
                .artifacts
                .iter()
                .zip(&allowed)
                .any(|(a, path)| &a.path != path || !valid_hash(&a.sha256))
        {
            return Err(bad("primer_trace_artifact_paths"));
        }
        for a in &file.artifacts {
            bytes = bytes
                .checked_add(a.bytes)
                .filter(|n| *n <= report.request.budgets.max_sidecar_bytes)
                .ok_or_else(|| bad("primer_trace_sidecar_budget"))?;
        }
    }
    let expected = summary(&report.files_completed)?;
    if encode(&expected)? != encode(&report.summary)? {
        return Err(bad("primer_trace_summary_binding"));
    }
    let allowed = if report.reference_panel_sha256.is_some() {
        vec!["primer_panel.json", "reference_panel.json", "summary.json"]
    } else {
        vec!["primer_panel.json", "summary.json"]
    };
    if report.artifacts.len() != allowed.len()
        || report
            .artifacts
            .iter()
            .zip(allowed)
            .any(|(a, p)| a.path != p || !valid_hash(&a.sha256))
    {
        return Err(bad("primer_trace_metadata_artifacts"));
    }
    if let Some(error) = &report.control_error {
        let valid_phase = match error.phase {
            RnaReadPanelControlErrorPhase::Finalization => {
                error.file_id.is_none() && report.inputs_not_completed.is_empty()
            }
            RnaReadPanelControlErrorPhase::FileBoundary
            | RnaReadPanelControlErrorPhase::Streaming => {
                error.file_id.as_deref()
                    == Some(format!("file_{:04}", report.files_completed.len()).as_str())
                    && !report.inputs_not_completed.is_empty()
            }
        };
        if error.code != "panel_control_marker_error"
            || error.message.trim().is_empty()
            || report.request.cancel_path.as_deref() != Some(error.marker_path.as_str())
            || !valid_phase
        {
            return Err(bad("primer_trace_control_error_binding"));
        }
    }
    Ok(())
}

fn summary(files: &[RnaPrimerTraceFile]) -> Result<RnaPrimerTraceSummary, EngineError> {
    let mut summary = RnaPrimerTraceSummary {
        schema: "gentle.rna_primer_trace_summary.v1".into(),
        histogram_policy:
            "all_retained_full_hits; exact_0_to_256_then_right_closed_doubling_to_2000000".into(),
        by_sample_barcode: BTreeMap::new(),
    };
    for file in files {
        let key = serde_json::to_string(&(
            &file.input.sample_id,
            &file.input.barcode,
            &file.input.priming,
        ))
        .unwrap();
        if !summary.by_sample_barcode.contains_key(&key) && summary.by_sample_barcode.len() >= 128 {
            return Err(bad("primer_trace_sample_budget"));
        }
        let total = summary.by_sample_barcode.entry(key).or_default();
        total.reads = total
            .reads
            .checked_add(file.counts.reads)
            .ok_or_else(|| bad("primer_trace_count_overflow"))?;
        for (to, from) in [
            (&mut total.anchor_classes, &file.counts.anchor_classes),
            (&mut total.locus_relations, &file.counts.locus_relations),
            (
                &mut total.full_hits_by_primer_class_strand,
                &file.counts.full_hits_by_primer_class_strand,
            ),
            (
                &mut total.anchor_start_histogram,
                &file.counts.anchor_start_histogram,
            ),
            (
                &mut total.remaining_tail_histogram,
                &file.counts.remaining_tail_histogram,
            ),
        ] {
            for (k, v) in from {
                let entry = to.entry(k.clone()).or_default();
                *entry = entry
                    .checked_add(*v)
                    .ok_or_else(|| bad("primer_trace_count_overflow"))?;
            }
        }
    }
    Ok(summary)
}

impl GentleEngine {
    /// Seal only a caller-supplied panel; never manufacture oligos or assume kit identity.
    pub fn seal_rna_primer_panel(panel: &mut RnaPrimerPanel) -> Result<(), EngineError> {
        panel.content_sha256.clear();
        panel.content_sha256 = digest(panel)?;
        validate(panel)?;
        Ok(())
    }
    /// Strict validation, or explicit read-only sealing of an inline draft with an empty digest.
    pub fn validate_rna_primer_panel(
        binding: &RnaPrimerPanelBinding,
        seal: bool,
    ) -> Result<RnaPrimerPanel, EngineError> {
        if seal {
            let RnaPrimerPanelBinding::Inline { panel } = binding else {
                return Err(bad("primer_panel_seal_requires_inline_draft"));
            };
            if !panel.content_sha256.is_empty() {
                return Err(bad("primer_panel_seal_requires_empty_digest"));
            }
            let mut panel = *panel.clone();
            Self::seal_rna_primer_panel(&mut panel)?;
            return Ok(panel);
        }
        resolve(binding)
    }

    /// Stream an independent trace into new artifacts. Errors never publish a final report;
    /// cancellation/control faults retain only verified units, without changing the project.
    pub fn trace_primer_reads(
        &self,
        mut request: TracePrimerReadsRequest,
        progress: &mut dyn FnMut(OperationProgress) -> bool,
    ) -> Result<RnaPrimerTraceReport, EngineError> {
        let submitted = digest(&request)?;
        let b = &request.budgets;
        if request.schema != RNA_PRIMER_TRACE_REQUEST_SCHEMA
            || !(1..=2_000_000).contains(&b.max_read_nt)
            || !(1..=1024).contains(&b.max_hits_per_read)
            || !(1..=8192).contains(&b.max_reference_sites)
            || !(1..=4_000_000).contains(&b.max_index_postings)
            || !(1..=500_000_000).contains(&b.max_scan_cells)
            || !(1..=4_000_000_000).contains(&b.max_sidecar_bytes)
        {
            return Err(bad("primer_trace_request_or_budget"));
        }
        let panel = resolve(&request.primer_panel)?;
        let mut warnings = validate(&panel)?;
        let reference = if let Some(binding) = &request.reference_panel {
            let panel: RnaReadPanel = load(&binding.path)?;
            validate_panel(&panel)?;
            if !valid_hash(&binding.expected_sha256)
                || panel.content_sha256 != binding.expected_sha256
                || !panel.overlapping_source_pairs.is_empty()
            {
                return Err(bad("primer_trace_reference_binding"));
            }
            Some(panel)
        } else {
            None
        };
        // Reuse the established file/checksum admission, without running its nomination policy.
        let mut admission = ScreenRnaReadPanelRequest {
            schema: RNA_PANEL_SCREEN_REQUEST_SCHEMA.into(),
            panel_path: String::new(),
            expected_panel_sha256: None,
            files: request.files.clone(),
            md5_manifest_path: request.md5_manifest_path.clone(),
            output_dir: request.output_dir.clone(),
            policy: RnaReadPanelPolicy::default(),
            cancel_path: request.cancel_path.clone(),
        };
        validate_request(&mut admission)?;
        request.files = admission.files;
        let identities = [submitted.clone(), digest(&request)?];
        let pre_cancelled = control::marker_requested(request.cancel_path.as_deref(), &identities)?;
        let mut preparation_cancelled = pre_cancelled;
        let mut preparation_error = None;
        let mut preparation_gate = control::PollGate::default();
        let mut preparation_stop = || {
            if !preparation_gate.due(std::time::Instant::now(), false) {
                return false;
            }
            {
                match control::marker_requested(request.cancel_path.as_deref(), &identities) {
                    Ok(true) => {
                        preparation_cancelled = true;
                        return true;
                    }
                    Ok(false) => {}
                    Err(error) => {
                        preparation_error = Some(control::diagnostic(
                            error,
                            request.cancel_path.as_deref(),
                            RnaReadPanelControlErrorPhase::FileBoundary,
                            Some("file_0000"),
                        ));
                        preparation_cancelled = true;
                        return true;
                    }
                }
            }
            if !progress(OperationProgress::Workflow {
                completed: 0,
                total: request.files.len(),
            }) {
                preparation_cancelled = true;
                true
            } else {
                false
            }
        };
        let sites = if let Some(reference) = reference.as_ref().filter(|_| !pre_cancelled) {
            match matching::sites(&panel, reference, &request.budgets, &mut preparation_stop) {
                Ok((sites, more)) => {
                    warnings.extend(more);
                    sites
                }
                Err(error) if error.message == "primer_trace_cancelled" => Default::default(),
                Err(error) => return Err(error),
            }
        } else {
            Default::default()
        };
        drop(preparation_stop);
        let output = PathBuf::from(&request.output_dir);
        let parent = output
            .parent()
            .filter(|p| !p.as_os_str().is_empty())
            .unwrap_or_else(|| Path::new("."));
        let staging = tempfile::Builder::new()
            .prefix(".gentle-primer-trace-")
            .tempdir_in(parent)
            .map_err(|e| bad(e.to_string()))?;
        let mut metadata_artifacts = vec![write_json(
            &staging.path().join("primer_panel.json"),
            &panel,
        )?];
        if let Some(reference) = &reference {
            metadata_artifacts.push(write_json(
                &staging.path().join("reference_panel.json"),
                reference,
            )?);
        }
        let mut report=RnaPrimerTraceReport {schema:RNA_PRIMER_TRACE_SCHEMA.into(),content_sha256:String::new(),algorithm:ALGORITHM.into(),
            submitted_request_sha256:submitted,status:if preparation_error.is_some() {"control_error"} else if preparation_cancelled {"cancelled"}else{"completed"}.into(),
            request:request.clone(),primer_panel:panel.clone(),reference_panel_sha256:reference.as_ref().map(|p|p.content_sha256.clone()),
            artifacts:metadata_artifacts,files_completed:vec![],inputs_not_completed:vec![],control_error:preparation_error,
            summary:RnaPrimerTraceSummary {schema:"gentle.rna_primer_trace_summary.v1".into(),
                histogram_policy:"all_retained_full_hits; exact_0_to_256_then_right_closed_doubling_to_2000000".into(),by_sample_barcode:BTreeMap::new()},
            warnings,nonclaims:vec![
                "Primer carry-through and short local core/payload concordance are not enrichment, exclusive transcript assignment, whole-read compatibility, expression, original molecule counts or mechanism evidence.".into(),
                "Core bases may be oligo-supplied. Locus relation describes the primary anchor. Expected-target success retains all qualifying expected sites but does not trigger the general other-template probe search; it does not exclude panel competitors.".into(),
                "Expected-target and other-template results retain ambiguity within the declared panel only. Unknown flanks, short/ambiguous probes and missing sites are not zero support or proven boundary divergence.".into(),
                "Primers and policies are caller-reviewed, uncalibrated inputs; hashes bind content, not kit identity, authenticated provenance or assay approval.".into(),
                "Payload intervals are candidates, not segmented molecules. Multiple-anchor payloads are unresolved. FASTA is opt-in, normalized/oriented, and never changes the original reads.".into(),
                "Read/class denominators count each read once; full-hit positional histograms count all retained full hits and may exceed read counts. No scientific acceptance or throughput claim.".into()]};
        let used = Rc::new(Cell::new(0));
        for (ordinal, source) in request.files.iter().enumerate() {
            if report.status != "completed" {
                break;
            }
            let file_id = format!("file_{ordinal:04}");
            let unit = tempfile::tempdir_in(staging.path()).map_err(|e| bad(e.to_string()))?;
            let mut rows =
                csv::WriterBuilder::new()
                    .delimiter(b'\t')
                    .from_writer(LimitedWriter::new(
                        &unit.path().join("reads.tsv"),
                        &used,
                        request.budgets.max_sidecar_bytes,
                    )?);
            rows.write_record(HEADER).map_err(|e| bad(e.to_string()))?;
            let mut fasta = if request.export_payload_fasta {
                Some(LimitedWriter::new(
                    &unit.path().join("payload.fasta"),
                    &used,
                    request.budgets.max_sidecar_bytes,
                )?)
            } else {
                None
            };
            let mut counts = RnaPrimerTraceCounts::default();
            let mut gate = control::PollGate::default();
            let mut marker_error = None;
            let mut stop = |reads: u64, bytes: u64, total: u64, force: bool| -> bool {
                if !gate.due(std::time::Instant::now(), force) {
                    return false;
                }
                match control::marker_requested(request.cancel_path.as_deref(), &identities) {
                    Ok(true) => return true,
                    Ok(false) => {}
                    Err(error) => {
                        marker_error = Some(control::diagnostic(
                            error,
                            request.cancel_path.as_deref(),
                            if total == 0 {
                                RnaReadPanelControlErrorPhase::FileBoundary
                            } else {
                                RnaReadPanelControlErrorPhase::Streaming
                            },
                            Some(&file_id),
                        ));
                        return true;
                    }
                }
                !progress(OperationProgress::RnaReadPanel(RnaReadPanelProgress {
                    file_id: file_id.clone(),
                    sample_id: source.sample_id.clone(),
                    file_ordinal: ordinal + 1,
                    files_total: request.files.len(),
                    processed_reads: reads,
                    compressed_bytes_read: bytes,
                    compressed_bytes_total: total,
                    phase: "primer_trace_streaming_unverified".into(),
                }))
            };
            let input = if stop(0, 0, 0, true) {
                None
            } else {
                visit_fastq(Path::new(&source.path), &mut |record| {
                    if stop(
                        counts.reads,
                        record.compressed_bytes_read,
                        record.compressed_bytes_total,
                        counts.reads == 0,
                    ) {
                        return Ok(false);
                    }
                    let row = match matching::row(
                        &panel,
                        reference.as_ref(),
                        &sites,
                        &record,
                        counts.reads,
                        &request.budgets,
                        &mut || {
                            stop(
                                counts.reads,
                                record.compressed_bytes_read,
                                record.compressed_bytes_total,
                                false,
                            )
                        },
                    ) {
                        Ok(row) => row,
                        Err(e) if e.message == "primer_trace_cancelled" => return Ok(false),
                        Err(e) => {
                            return Err(bad(format!(
                                "{}; file_id={file_id}; record_index={}; read_id={}; no report published",
                                e.message, counts.reads, record.read_id
                            )));
                        }
                    };
                    rows.write_record([
                        row.record_index.to_string(),
                        row.read_id.clone(),
                        row.anchor_class.clone(),
                        row.anchor_position.clone().unwrap_or_default(),
                        row.primer_multiplicity.clone().unwrap_or_default(),
                        row.locus_relation.clone(),
                        row.evidence_tier.clone(),
                        serde_json::to_string(&row.anchor).unwrap(),
                        serde_json::to_string(&row.leader).unwrap(),
                        serde_json::to_string(&row.payload).unwrap(),
                        serde_json::to_string(&row.payload_transcript_sense).unwrap(),
                        serde_json::to_string(&row.payload_mean_phred).unwrap(),
                        serde_json::to_string(&row.payload_fraction_below_q10).unwrap(),
                        serde_json::to_string(&row.template_matches).unwrap(),
                        serde_json::to_string(&row.artifact_flags).unwrap(),
                        serde_json::to_string(&row.reasons).unwrap(),
                        row.hits.len().to_string(),
                        serde_json::to_string(&row.hits).unwrap(),
                        row.read_length.to_string(),
                        serde_json::to_string(&row.leader_length).unwrap(),
                        serde_json::to_string(&row.payload_length).unwrap(),
                    ])
                    .map_err(|e| bad(e.to_string()))?;
                    if let (Some(fasta), Some((a, b)), Some(hit)) =
                        (&mut fasta, row.payload, &row.anchor)
                    {
                        let payload = if hit.strand == PrimerReadStrand::ReverseComplement {
                            matching::rc(&record.sequence[a..b])
                        } else {
                            record.sequence[a..b].to_vec()
                        };
                        if !payload.is_empty() {
                            writeln!(fasta,">{file_id}:{} original_read={} interval={a}..{b} candidate_payload_not_segmented\n{}",row.record_index,row.read_id,String::from_utf8_lossy(&payload)).map_err(|e|bad(e.to_string()))?;
                        }
                    }
                    count(&mut counts, &row, record.sequence.len());
                    Ok(true)
                })?
            };
            drop(stop);
            let Some(input) = input else {
                report.status = if marker_error.is_some() {
                    "control_error"
                } else {
                    "cancelled"
                }
                .into();
                report.control_error = marker_error;
                break;
            };
            if source
                .expected_md5
                .as_ref()
                .is_some_and(|h| h != &input.md5)
                || source.expected_sha256.as_ref().is_some_and(|h| {
                    h.trim_start_matches("sha256:").to_ascii_lowercase()
                        != input.sha256.trim_start_matches("sha256:")
                })
            {
                return Err(bad(
                    "primer_trace_input_digest_mismatch: no report published",
                ));
            }
            let mut artifacts = vec![
                rows.into_inner()
                    .map_err(|e| bad(e.to_string()))?
                    .receipt(format!("{file_id}/reads.tsv"))?,
            ];
            if let Some(fasta) = fasta {
                artifacts.push(fasta.receipt(format!("{file_id}/payload.fasta"))?);
            }
            let unit_path = unit.keep();
            fs::rename(&unit_path, staging.path().join(&file_id))
                .map_err(|e| bad(e.to_string()))?;
            report.files_completed.push(RnaPrimerTraceFile {
                file_id,
                input: source.clone(),
                observed_sha256: input.sha256,
                observed_md5: input.md5,
                compressed_bytes: input.bytes,
                counts,
                artifacts,
            });
        }
        report.inputs_not_completed = request
            .files
            .iter()
            .skip(report.files_completed.len())
            .cloned()
            .collect();
        if report.status == "completed" {
            match control::marker_requested(request.cancel_path.as_deref(), &identities) {
                Ok(true) => report.status = "cancelled".into(),
                Ok(false) => {}
                Err(error) => {
                    report.status = "control_error".into();
                    report.control_error = Some(control::diagnostic(
                        error,
                        request.cancel_path.as_deref(),
                        RnaReadPanelControlErrorPhase::Finalization,
                        None,
                    ));
                }
            }
        }
        if !progress(OperationProgress::Workflow {
            completed: report.files_completed.len(),
            total: request.files.len(),
        }) && report.status == "completed"
        {
            report.status = "cancelled".into();
        }
        report.summary = summary(&report.files_completed)?;
        report.artifacts.push(write_json(
            &staging.path().join("summary.json"),
            &report.summary,
        )?);
        report.content_sha256 = digest(&report)?;
        validate_report(&report)?;
        if encode(&report)?.len() > 64_000_000 {
            return Err(bad("primer_trace_report_budget"));
        }
        fs::create_dir(&output).map_err(|e| bad(format!("primer_trace_output_reserve:{e}")))?;
        if let Err(e) = fs::rename(staging.path(), output.join("artifacts")) {
            let _ = fs::remove_dir(&output);
            return Err(bad(format!("primer_trace_publish:{e}")));
        }
        let mut final_file =
            tempfile::NamedTempFile::new_in(&output).map_err(|e| bad(e.to_string()))?;
        serde_json::to_writer_pretty(final_file.as_file_mut(), &report)
            .map_err(|e| bad(e.to_string()))?;
        final_file
            .as_file()
            .sync_all()
            .map_err(|e| bad(e.to_string()))?;
        final_file
            .persist_noclobber(output.join("report.json"))
            .map_err(|e| bad(e.to_string()))?;
        Ok(report)
    }

    /// Bind the actual published report to submission, allowing only checksum normalization.
    pub fn primer_trace_result_binding(
        submitted: &TracePrimerReadsRequest,
        report: &RnaPrimerTraceReport,
    ) -> Result<RnaPrimerTraceResultBinding, EngineError> {
        validate_report(report)?;
        if digest(submitted)? != report.submitted_request_sha256 {
            return Err(bad("primer_trace_submitted_binding"));
        }
        let mut effective = submitted.clone();
        if effective.files.len() != report.request.files.len() {
            return Err(bad("primer_trace_effective_binding"));
        }
        for (file, resolved) in effective.files.iter_mut().zip(&report.request.files) {
            if file.expected_md5.as_ref().is_some_and(|h| {
                resolved.expected_md5.as_deref() != Some(h.to_ascii_lowercase().as_str())
            }) {
                return Err(bad("primer_trace_checksum_binding"));
            }
            file.expected_md5 = resolved.expected_md5.clone();
        }
        if digest(&effective)? != digest(&report.request)? {
            return Err(bad("primer_trace_effective_binding"));
        }
        Ok(RnaPrimerTraceResultBinding {
            schema: BINDING_SCHEMA.into(),
            submitted_request_sha256: report.submitted_request_sha256.clone(),
            effective_request_sha256: digest(&report.request)?,
            primer_panel_sha256: report.primer_panel.content_sha256.clone(),
            reference_panel_sha256: report.reference_panel_sha256.clone(),
            report_content_sha256: report.content_sha256.clone(),
            report_path: Path::new(&submitted.output_dir)
                .join("report.json")
                .display()
                .to_string(),
            outcome: report.status.clone(),
        })
    }

    /// Metadata and a bounded unverified TSV prefix, never an implicit FASTQ or full-sidecar replay.
    pub fn inspect_primer_trace(
        path: &str,
        binding: Option<&RnaPrimerTraceResultBinding>,
        file_id: Option<&str>,
        rows: usize,
    ) -> Result<RnaPrimerTraceInspection, EngineError> {
        if rows > 1000 || (rows > 0 && file_id.is_none()) {
            return Err(bad("primer_trace_inspect_rows"));
        }
        let report: RnaPrimerTraceReport = load(path)?;
        validate_report(&report)?;
        if let Some(b) = binding {
            if b.schema != BINDING_SCHEMA
                || b.report_path != path
                || b.report_content_sha256 != report.content_sha256
                || b.submitted_request_sha256 != report.submitted_request_sha256
                || b.effective_request_sha256 != digest(&report.request)?
                || b.primer_panel_sha256 != report.primer_panel.content_sha256
                || b.reference_panel_sha256 != report.reference_panel_sha256
                || b.outcome != report.status
            {
                return Err(bad("primer_trace_expected_result_binding"));
            }
        }
        let mut inspected = vec![];
        if let Some(id) = file_id {
            let file = report
                .files_completed
                .iter()
                .find(|f| f.file_id == id)
                .ok_or_else(|| bad("primer_trace_unknown_file"))?;
            if rows > 0 {
                let relative = format!("{id}/reads.tsv");
                if !safe_id(id) || !file.artifacts.iter().any(|a| a.path == relative) {
                    return Err(bad("primer_trace_inspection_path"));
                }
                let parent = Path::new(path).parent().unwrap_or_else(|| Path::new("."));
                let mut current = parent.to_path_buf();
                for component in ["artifacts", id, "reads.tsv"] {
                    current.push(component);
                    if fs::symlink_metadata(&current)
                        .map_err(|e| bad(e.to_string()))?
                        .file_type()
                        .is_symlink()
                    {
                        return Err(bad("primer_trace_inspection_symlink"));
                    }
                }
                let input = fs::File::open(current)
                    .map_err(|e| bad(e.to_string()))?
                    .take(4_000_000);
                let mut reader = csv::ReaderBuilder::new()
                    .delimiter(b'\t')
                    .from_reader(input);
                for record in reader.records().take(rows) {
                    let r = record.map_err(|e| bad(e.to_string()))?;
                    inspected.push(parse_row(&r)?);
                }
            }
        }
        Ok(RnaPrimerTraceInspection {
            report: Box::new(report),
            rows: inspected,
            artifact_content_verified: false,
        })
    }

    /// Publish a non-overwriting request-bound marker; requesting cancellation is not termination.
    pub fn cancel_primer_trace(request: &TracePrimerReadsRequest) -> Result<String, EngineError> {
        if request.schema != RNA_PRIMER_TRACE_REQUEST_SCHEMA {
            return Err(bad("primer_trace_cancel_schema"));
        }
        let path = request
            .cancel_path
            .as_deref()
            .ok_or_else(|| bad("primer_trace_cancel_path_required"))?;
        let identity = digest(request)?;
        if control::marker_requested(Some(path), std::slice::from_ref(&identity))? {
            return Ok("Cancellation already requested; termination not confirmed".into());
        }
        let parent = Path::new(path)
            .parent()
            .filter(|p| !p.as_os_str().is_empty())
            .unwrap_or_else(|| Path::new("."));
        let mut marker = tempfile::NamedTempFile::new_in(parent).map_err(|e| bad(e.to_string()))?;
        write!(marker, "{identity}").map_err(|e| bad(e.to_string()))?;
        marker
            .as_file()
            .sync_all()
            .map_err(|e| bad(e.to_string()))?;
        marker
            .persist_noclobber(path)
            .map_err(|e| bad(e.to_string()))?;
        Ok("Cancellation requested for this exact request; termination not confirmed".into())
    }
}

#[cfg(test)]
#[path = "primer_trace/tests.rs"]
mod tests;
