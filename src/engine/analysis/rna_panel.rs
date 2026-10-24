//! One-pass, reference-bound multi-locus FASTQ nomination, not a second aligner.
//!
//! All sequence/geometry construction uses the shared mature-transcript builder.
//! A read increments one equivalence set. No annotation-order winner, expression
//! estimate, automatic trimming or whole-transcriptome specificity is inferred.

use super::*;
use std::io::{Read, Write};
use std::path::{Path, PathBuf};

#[path = "rna_panel/control.rs"]
mod control;
#[path = "rna_panel/input.rs"]
mod input;
#[path = "rna_panel/summary.rs"]
mod summary;
#[path = "rna_panel/tools.rs"]
mod tools;
#[path = "rna_panel/verification.rs"]
mod verification;
use input::{HashedWriter, visit_fastq};

fn bad(message: impl Into<String>) -> EngineError {
    EngineError::invalid_input(message)
}
fn encode<T: Serialize>(value: &T) -> Result<Vec<u8>, EngineError> {
    let value = serde_json::to_value(value).map_err(|e| bad(format!("panel_json: {e}")))?;
    serde_json::to_vec(&value).map_err(|e| bad(format!("panel_json: {e}")))
}
fn safe_id(value: &str) -> bool {
    !value.is_empty()
        && value.len() <= 128
        && value
            .bytes()
            .all(|c| c.is_ascii_alphanumeric() || b"_-.:".contains(&c))
}
fn load<T: serde::de::DeserializeOwned>(path: &str) -> Result<T, EngineError> {
    let file = fs::File::open(path).map_err(|e| bad(format!("panel_open: {path}: {e}")))?;
    if file.metadata().map_err(|e| bad(e.to_string()))?.len() > 128_000_000 {
        return Err(bad(
            "panel_json_budget: at most 128 MB of reference metadata",
        ));
    }
    let mut bytes = vec![];
    file.take(128_000_001)
        .read_to_end(&mut bytes)
        .map_err(|e| bad(e.to_string()))?;
    if bytes.len() > 128_000_000 {
        return Err(bad("panel_json_budget"));
    }
    serde_json::from_slice(&bytes).map_err(|e| bad(format!("panel_json: {e}")))
}
fn write_json<T: Serialize>(path: &Path, value: &T) -> Result<RnaReadPanelArtifact, EngineError> {
    let mut writer = HashedWriter::create(path)?;
    serde_json::to_writer_pretty(&mut writer, value)
        .map_err(|e| bad(format!("panel_json_write: {e}")))?;
    writer.receipt(path.file_name().unwrap().to_string_lossy().into_owned())
}
fn seal_panel(panel: &mut RnaReadPanel) -> Result<(), EngineError> {
    panel.content_sha256.clear();
    panel.content_sha256 = sha256_prefixed_bytes(&encode(panel)?);
    Ok(())
}

fn validate_screen_report_binding(report: &RnaReadPanelScreenReport) -> Result<(), EngineError> {
    let mut basis = report.clone();
    basis.content_sha256.clear();
    if !matches!(
        report.schema.as_str(),
        RNA_PANEL_SCREEN_SCHEMA | RNA_PANEL_SCREEN_LEGACY_SCHEMA
    ) || report.content_sha256 != sha256_prefixed_bytes(&encode(&basis)?)
        || report.isoform_assignment_performed
        || report.alignment_confirmation_performed
    {
        return Err(bad("panel_report_binding"));
    }
    let valid_outcome = match (report.status.as_str(), &report.control_error) {
        ("completed" | "cancelled", None) => true,
        ("control_error", Some(error)) => {
            report.schema == RNA_PANEL_SCREEN_SCHEMA
                && error.code == "panel_control_marker_error"
                && !error.message.trim().is_empty()
                && report.request.cancel_path.as_deref() == Some(error.marker_path.as_str())
                && report.files_completed.len() <= report.request.files.len()
                && match error.phase {
                    RnaReadPanelControlErrorPhase::Finalization => {
                        error.file_id.is_none()
                            && report.files_completed.len() == report.request.files.len()
                    }
                    RnaReadPanelControlErrorPhase::FileBoundary
                    | RnaReadPanelControlErrorPhase::Streaming => {
                        report.files_completed.len() < report.request.files.len()
                            && error.file_id.as_deref()
                                == Some(
                                    format!("file_{:04}", report.files_completed.len()).as_str(),
                                )
                    }
                }
        }
        _ => false,
    };
    if !valid_outcome {
        return Err(bad("panel_report_control_outcome_binding"));
    }
    Ok(())
}

fn source_overlaps(loci: &[RnaReadPanelLocus]) -> Vec<(String, String)> {
    let mut pairs = vec![];
    for (i, a) in loci.iter().enumerate() {
        for b in &loci[i + 1..] {
            // Conservative admission: two seeds on the same loaded record can
            // expand to the same group even when their local feature IDs differ.
            let same_record = a.source.seq_id == b.source.seq_id;
            let same_declared_gene =
                a.source.assembly == b.source.assembly && a.source.gene_id == b.source.gene_id;
            let overlap = match (&a.genome_anchor, &b.genome_anchor) {
                (Some(x), Some(y)) => {
                    a.source.assembly == b.source.assembly
                        && x.chromosome == y.chromosome
                        && x.strand == y.strand
                        && x.start_1based <= y.end_1based
                        && y.start_1based <= x.end_1based
                }
                _ => false,
            };
            if same_record || same_declared_gene || overlap {
                pairs.push((a.source.locus_id.clone(), b.source.locus_id.clone()));
            }
        }
    }
    pairs
}

impl GentleEngine {
    /// Build a portable panel from independent loaded annotated loci, without changing project state.
    pub fn build_rna_read_panel(
        &self,
        request: &BuildRnaReadPanelRequest,
    ) -> Result<RnaReadPanel, EngineError> {
        if request.schema != RNA_PANEL_BUILD_SCHEMA
            || request.sources.is_empty()
            || request.sources.len() > 128
        {
            return Err(bad(
                "panel_sources: expected build v1 and 1..=128 explicit sources",
            ));
        }
        let mut panel = RnaReadPanel {
            schema: RNA_PANEL_SCHEMA.into(),
            content_sha256: String::new(),
            loci: vec![],
            templates: vec![],
            locus_count: 0,
            gene_count: 0,
            selected_transcript_count: 0,
            annotated_transcript_count: 0,
            overlapping_source_pairs: vec![],
        };
        let mut ids = BTreeSet::new();
        for source in &request.sources {
            if !safe_id(&source.locus_id)
                || !ids.insert(source.locus_id.clone())
                || [
                    &source.source,
                    &source.assembly,
                    &source.annotation_release,
                    &source.gene_id,
                    &source.gene_symbol,
                ]
                .iter()
                .any(|s| s.trim().is_empty())
            {
                return Err(bad(
                    "panel_source_identity: unique safe locus IDs and explicit source/assembly/release/gene identities required",
                ));
            }
            let dna = self
                .state
                .sequences
                .get(&source.seq_id)
                .ok_or_else(|| bad(format!("panel_source_missing: {}", source.seq_id)))?;
            if dna.is_protein_sequence() {
                return Err(bad("panel_source_requires_dna"));
            }
            let view = self.build_splicing_expert_view(
                &source.seq_id,
                source.source_feature_id,
                SplicingScopePreset::TargetGroupTargetStrand,
            )?;
            let templates = Self::build_qpcr_transcript_design_templates(dna, &view)?;
            if templates.is_empty() {
                return Err(bad(format!("panel_no_transcripts: {}", source.locus_id)));
            }
            let selected: BTreeSet<_> = source.selected_transcript_ids.iter().cloned().collect();
            if selected.len() != source.selected_transcript_ids.len()
                || selected
                    .iter()
                    .any(|id| !templates.iter().any(|t| &t.transcript_id == id))
            {
                return Err(bad(
                    "panel_transcript_selection: duplicate or unresolved selected transcript ID",
                ));
            }
            panel.annotated_transcript_count += templates.len();
            panel.loci.push(RnaReadPanelLocus {
                source: source.clone(),
                source_record_sha256: sha256_prefixed_bytes(&encode(&(
                    dna.get_forward_string(),
                    dna.features(),
                ))?),
                genome_anchor: self.transcript_qpcr_panel_source_anchor(&source.seq_id, dna),
                annotated_transcript_count: templates.len(),
                annotated_transcript_ids: templates
                    .iter()
                    .map(|t| t.transcript_id.clone())
                    .collect(),
            });
            for template in templates {
                if !selected.is_empty() && !selected.contains(&template.transcript_id) {
                    continue;
                }
                panel.templates.push(RnaReadPanelTemplate {
                    member_id: format!("{}:{}", source.locus_id, template.transcript_feature_id),
                    locus_id: source.locus_id.clone(),
                    transcript_feature_id: template.transcript_feature_id,
                    transcript_id: template.transcript_id,
                    strand: template.strand,
                    sequence_sha256: sha256_prefixed_str(&template.sequence),
                    sequence: template.sequence,
                    junction_offsets: template
                        .local_exon_segments
                        .iter()
                        .skip(1)
                        .map(|e| e.local_start_0based)
                        .collect(),
                    source_exons: template.exon_chain,
                });
            }
        }
        panel
            .loci
            .sort_by(|a, b| a.source.locus_id.cmp(&b.source.locus_id));
        panel
            .templates
            .sort_by(|a, b| a.member_id.cmp(&b.member_id));
        panel.locus_count = panel.loci.len();
        panel.gene_count = panel
            .loci
            .iter()
            .map(|l| (&l.source.assembly, &l.source.gene_id))
            .collect::<BTreeSet<_>>()
            .len();
        panel.selected_transcript_count = panel.templates.len();
        panel.overlapping_source_pairs = source_overlaps(&panel.loci);
        seal_panel(&mut panel)?;
        validate_panel(&panel)?;
        let destination = Path::new(&request.output_path);
        let parent = destination
            .parent()
            .filter(|p| !p.as_os_str().is_empty())
            .unwrap_or_else(|| Path::new("."));
        let mut file = tempfile::NamedTempFile::new_in(parent)
            .map_err(|e| bad(format!("panel_output: {e}")))?;
        serde_json::to_writer_pretty(file.as_file_mut(), &panel).map_err(|e| bad(e.to_string()))?;
        file.as_file().sync_all().map_err(|e| bad(e.to_string()))?;
        file.persist_noclobber(destination)
            .map_err(|e| bad(format!("panel_output_exists_or_io: {e}")))?;
        Ok(panel)
    }

    /// Nominate all supported members; false from progress requests cooperative cancellation.
    /// The legacy single-locus winner is never used here.
    pub fn screen_rna_read_panel(
        &self,
        request: ScreenRnaReadPanelRequest,
        on_progress: &mut dyn FnMut(OperationProgress) -> bool,
    ) -> Result<RnaReadPanelScreenReport, EngineError> {
        screen(request, on_progress)
    }
}

fn validate_panel(panel: &RnaReadPanel) -> Result<(), EngineError> {
    let mut copy = panel.clone();
    seal_panel(&mut copy)?;
    if panel.schema != RNA_PANEL_SCHEMA
        || panel.content_sha256 != copy.content_sha256
        || panel.templates.is_empty()
        || panel.templates.len() > 4096
        || panel.loci.is_empty()
        || panel.loci.len() > 128
        || panel.locus_count != panel.loci.len()
        || panel.selected_transcript_count != panel.templates.len()
    {
        return Err(bad(
            "panel_binding: invalid schema, digest, counts or reference limits",
        ));
    }
    let loci = panel
        .loci
        .iter()
        .map(|l| l.source.locus_id.as_str())
        .collect::<BTreeSet<_>>();
    if loci.len() != panel.loci.len()
        || loci.iter().any(|s| !safe_id(s))
        || panel.loci.iter().any(|l| {
            !l.source_record_sha256.starts_with("sha256:")
                || !valid_digest(l.source_record_sha256.trim_start_matches("sha256:"), 64)
                || l.annotated_transcript_ids.len() != l.annotated_transcript_count
                || panel
                    .templates
                    .iter()
                    .filter(|t| t.locus_id == l.source.locus_id)
                    .count()
                    > l.annotated_transcript_count
        })
        || panel.overlapping_source_pairs != source_overlaps(&panel.loci)
        || Some(panel.annotated_transcript_count)
            != panel.loci.iter().try_fold(0usize, |sum, l| {
                sum.checked_add(l.annotated_transcript_count)
            })
        || panel.gene_count
            != panel
                .loci
                .iter()
                .map(|l| (&l.source.assembly, &l.source.gene_id))
                .collect::<BTreeSet<_>>()
                .len()
    {
        return Err(bad("panel_locus_binding"));
    }
    let mut ids = BTreeSet::new();
    let mut bases = 0usize;
    for template in &panel.templates {
        bases = bases
            .checked_add(template.sequence.len())
            .ok_or_else(|| bad("panel_reference_budget"))?;
        let mut cursor = 0usize;
        let mut junctions = vec![];
        for &(start, end) in &template.source_exons {
            if end <= start {
                return Err(bad("panel_exon_geometry"));
            }
            if cursor > 0 {
                junctions.push(cursor);
            }
            cursor = cursor
                .checked_add(end - start)
                .ok_or_else(|| bad("panel_exon_geometry"))?;
        }
        if !ids.insert(&template.member_id)
            || template.member_id
                != format!("{}:{}", template.locus_id, template.transcript_feature_id)
            || !loci.contains(template.locus_id.as_str())
            || !panel.loci.iter().any(|l| {
                l.source.locus_id == template.locus_id
                    && l.annotated_transcript_ids.contains(&template.transcript_id)
            })
            || template.sequence.is_empty()
            || template.sequence.len() > 2_000_000
            || !template
                .sequence
                .bytes()
                .all(|b| b"ACGTRYSWKMBDHVN".contains(&b))
            || template.sequence_sha256 != sha256_prefixed_str(&template.sequence)
            || cursor != template.sequence.len()
            || junctions != template.junction_offsets
            || !matches!(template.strand.as_str(), "+" | "-")
        {
            return Err(bad(
                "panel_template_binding: invalid member, sequence, orientation or exon geometry",
            ));
        }
    }
    if bases > 32_000_000 {
        return Err(bad("panel_reference_budget: 32 million transcript bases"));
    }
    Ok(())
}

type Posting = (u32, u32);
type Index = HashMap<u32, Vec<Posting>>;
fn build_index(
    panel: &RnaReadPanel,
    policy: &RnaReadPanelPolicy,
) -> Result<(Index, usize), EngineError> {
    let mut index = Index::new();
    let mut count = 0;
    for (template_idx, template) in panel.templates.iter().enumerate() {
        for (pos, word) in template
            .sequence
            .as_bytes()
            .windows(policy.kmer_len)
            .enumerate()
        {
            let Some(key) = GentleEngine::encode_kmer_bits(word) else {
                continue;
            };
            count += 1;
            if count > policy.max_index_postings {
                return Err(bad(
                    "panel_posting_budget: no postings were silently discarded",
                ));
            }
            index
                .entry(key)
                .or_default()
                .push((template_idx as u32, pos as u32));
        }
    }
    if count == 0 {
        return Err(bad("panel_no_canonical_seeds"));
    }
    Ok((index, count))
}

fn nominate(
    panel: &RnaReadPanel,
    index: &Index,
    policy: &RnaReadPanelPolicy,
    read: &[u8],
    cancel: &mut dyn FnMut() -> bool,
) -> Result<(String, Vec<RnaReadPanelCandidate>), EngineError> {
    let mut candidates = vec![];
    let mut anchors_used = 0;
    let mut saw_seed = false;
    let rc = GentleEngine::reverse_complement_bytes(read);
    for (reverse, sequence) in [(false, read), (true, rc.as_slice())] {
        let mut matches = BTreeMap::<u32, Vec<(u32, u32)>>::new();
        for (read_pos, word) in sequence.windows(policy.kmer_len).enumerate() {
            if read_pos % 4096 == 0 && cancel() {
                return Err(bad("panel_cancelled"));
            }
            let Some(key) = GentleEngine::encode_kmer_bits(word) else {
                continue;
            };
            if let Some(postings) = index.get(&key) {
                saw_seed = true;
                anchors_used += postings.len();
                if anchors_used > policy.max_anchors_per_read {
                    return Err(bad(
                        "panel_anchor_budget: read incomplete, never an exclusive nomination",
                    ));
                }
                for &(template, pos) in postings {
                    matches
                        .entry(template)
                        .or_default()
                        .push((read_pos as u32, pos));
                }
            }
        }
        for (template_idx, mut anchors) in matches {
            if cancel() {
                return Err(bad("panel_cancelled"));
            }
            anchors.sort_unstable();
            anchors.dedup();
            let chain = bio::alignment::sparse::lcskpp(&anchors, policy.kmer_len);
            if chain.score < policy.min_chain_bases {
                continue;
            }
            let template = &panel.templates[template_idx as usize];
            let start = chain
                .path
                .iter()
                .map(|&i| anchors[i].0 as usize)
                .min()
                .unwrap();
            let end = chain
                .path
                .iter()
                .map(|&i| anchors[i].0 as usize + policy.kmer_len)
                .max()
                .unwrap();
            let supported_junction_offsets = template
                .junction_offsets
                .iter()
                .copied()
                .filter(|&junction| {
                    chain.path.iter().any(|&i| {
                        let pos = anchors[i].1 as usize;
                        pos + policy.min_junction_flank_bp <= junction
                            && junction + policy.min_junction_flank_bp <= pos + policy.kmer_len
                    })
                })
                .collect();
            candidates.push(RnaReadPanelCandidate {
                member_id: template.member_id.clone(),
                locus_id: template.locus_id.clone(),
                transcript_id: template.transcript_id.clone(),
                read_orientation: if reverse {
                    "reverse_complement"
                } else {
                    "as_read"
                }
                .into(),
                chain_bases: chain.score,
                read_interval_0based: if reverse {
                    (read.len() - end, read.len() - start)
                } else {
                    (start, end)
                },
                supported_junction_offsets,
            });
        }
    }
    candidates.sort_by(|a, b| {
        a.member_id
            .cmp(&b.member_id)
            .then(a.read_orientation.cmp(&b.read_orientation))
    });
    let loci = candidates
        .iter()
        .map(|c| c.locus_id.as_str())
        .collect::<BTreeSet<_>>();
    let mut class = match loci.len() {
        0 if saw_seed => "partial",
        0 => "unmatched",
        1 => "single_locus_candidate",
        _ => "multi_locus_ambiguous",
    };
    if candidates.iter().any(|a| {
        candidates.iter().any(|b| {
            a.locus_id != b.locus_id
                && (a.read_interval_0based.1 <= b.read_interval_0based.0
                    || b.read_interval_0based.1 <= a.read_interval_0based.0)
        })
    }) {
        class = "multi_segment_candidate";
    }
    Ok((class.into(), candidates))
}

struct Md5Manifest {
    exact: BTreeMap<String, Vec<String>>,
    basename: BTreeMap<String, Vec<String>>,
}

fn load_md5_manifest(path: &str) -> Result<Md5Manifest, EngineError> {
    let file = fs::File::open(path).map_err(|e| bad(format!("panel_md5_manifest: {e}")))?;
    let mut bytes = vec![];
    file.take(16_000_001)
        .read_to_end(&mut bytes)
        .map_err(|e| bad(e.to_string()))?;
    if bytes.len() > 16_000_000 {
        return Err(bad("panel_md5_manifest_budget"));
    }
    let text = std::str::from_utf8(&bytes).map_err(|_| bad("panel_md5_manifest_utf8"))?;
    let mut result = Md5Manifest {
        exact: BTreeMap::new(),
        basename: BTreeMap::new(),
    };
    for line in text.lines().filter(|line| !line.trim().is_empty()) {
        let Some((hash, raw)) = line.split_once(char::is_whitespace) else {
            return Err(bad("panel_md5_manifest_syntax"));
        };
        let key = raw
            .trim_start()
            .strip_prefix('*')
            .unwrap_or(raw.trim_start());
        if !valid_digest(hash, 32) || key.is_empty() {
            return Err(bad("panel_md5_manifest_syntax"));
        }
        result
            .exact
            .entry(key.into())
            .or_default()
            .push(hash.to_ascii_lowercase());
        result
            .basename
            .entry(key.rsplit(['/', '\\']).next().unwrap().into())
            .or_default()
            .push(hash.to_ascii_lowercase());
    }
    Ok(result)
}

fn expected_md5(
    input: &RnaReadPanelInput,
    manifest: Option<&Md5Manifest>,
) -> Result<Option<String>, EngineError> {
    if let Some(value) = &input.expected_md5 {
        return Ok(Some(value.to_ascii_lowercase()));
    }
    let Some(manifest) = manifest else {
        return Ok(None);
    };
    let filename = Path::new(&input.path)
        .file_name()
        .ok_or_else(|| bad("panel_input_filename"))?
        .to_string_lossy();
    // Folder renaming is harmless; colliding basenames require the exact supplier key.
    let found = match &input.checksum_manifest_key {
        Some(key) => manifest.exact.get(key),
        None => manifest.basename.get(filename.as_ref()),
    }
    .cloned()
    .unwrap_or_default();
    if found.len() != 1 {
        let keys = manifest
            .exact
            .keys()
            .filter(|key| {
                input.checksum_manifest_key.as_ref().map_or_else(
                    || key.rsplit(['/', '\\']).next() == Some(filename.as_ref()),
                    |expected| expected == *key,
                )
            })
            .take(8)
            .cloned()
            .collect::<Vec<_>>();
        return Err(bad(format!(
            "panel_md5_manifest_binding: expected one entry for {}, found {}; candidate keys (at most 8): {}",
            input.path,
            found.len(),
            serde_json::to_string(&keys).map_err(|e| bad(e.to_string()))?
        )));
    }
    Ok(found.into_iter().next())
}

fn valid_digest(value: &str, n: usize) -> bool {
    value.len() == n && value.bytes().all(|c| c.is_ascii_hexdigit())
}
fn validate_request(request: &mut ScreenRnaReadPanelRequest) -> Result<(), EngineError> {
    let p = &request.policy;
    if request.schema != RNA_PANEL_SCREEN_REQUEST_SCHEMA
        || request.files.is_empty()
        || request.files.len() > 4096
        || !(8..=16).contains(&p.kmer_len)
        || p.min_chain_bases < p.kmer_len as u32
        || p.min_junction_flank_bp == 0
        || p.min_junction_flank_bp > p.kmer_len / 2
        || !(1..=4_000_000).contains(&p.max_index_postings)
        || !(1..=250_000).contains(&p.max_anchors_per_read)
        || !(1..=16_384).contains(&p.max_equivalence_classes)
    {
        return Err(bad(
            "panel_admission: schema/files/policy outside declared bounds",
        ));
    }
    if request.output_dir.trim().is_empty() || Path::new(&request.output_dir).exists() {
        return Err(bad(
            "panel_output_exists_or_empty: choose a new run directory",
        ));
    }
    let mut paths = BTreeSet::new();
    let mut samples = BTreeMap::new();
    let manifest = if request.files.iter().any(|f| f.expected_md5.is_none()) {
        request
            .md5_manifest_path
            .as_deref()
            .map(load_md5_manifest)
            .transpose()?
    } else {
        None
    };
    for input in &mut request.files {
        let path = fs::canonicalize(&input.path).map_err(|e| bad(format!("panel_input: {e}")))?;
        if !paths.insert(path)
            || !safe_id(&input.sample_id)
            || input.barcode.trim().is_empty()
            || input.priming.trim().is_empty()
            || !input.declared_complete
        {
            return Err(bad(
                "panel_input_admission: distinct files, sample/barcode/priming and explicit completed-copy declaration required",
            ));
        }
        let binding = (input.barcode.clone(), input.priming.clone());
        if samples
            .insert(input.sample_id.clone(), binding.clone())
            .is_some_and(|old| old != binding)
        {
            return Err(bad("panel_sample_binding_conflict"));
        }
        input.expected_md5 = expected_md5(input, manifest.as_ref())?;
        if input
            .expected_md5
            .as_ref()
            .is_some_and(|s| !valid_digest(s, 32))
            || input
                .expected_sha256
                .as_ref()
                .is_some_and(|s| !valid_digest(s.trim_start_matches("sha256:"), 64))
            || (input.expected_md5.is_none() && input.expected_sha256.is_none())
        {
            return Err(bad(
                "panel_expected_digest: MD5 manifest/explicit MD5 or SHA256 required; verified during streaming, not upfront",
            ));
        }
    }
    Ok(())
}

fn screen(
    request: ScreenRnaReadPanelRequest,
    progress: &mut dyn FnMut(OperationProgress) -> bool,
) -> Result<RnaReadPanelScreenReport, EngineError> {
    screen_with_clock(request, progress, std::time::Instant::now)
}

// Clock injection keeps control/polling regressions deterministic without altering production cadence.
fn screen_with_clock(
    mut request: ScreenRnaReadPanelRequest,
    progress: &mut dyn FnMut(OperationProgress) -> bool,
    mut now: impl FnMut() -> std::time::Instant,
) -> Result<RnaReadPanelScreenReport, EngineError> {
    let original_cancel_id = control::cancel_identity(&request)?;
    validate_request(&mut request)?;
    let cancel_ids = [original_cancel_id, control::cancel_identity(&request)?];
    control::marker_requested(request.cancel_path.as_deref(), &cancel_ids)?;
    let panel: RnaReadPanel = load(&request.panel_path)?;
    validate_panel(&panel)?;
    if !panel.overlapping_source_pairs.is_empty() {
        return Err(bad(
            "panel_overlapping_sources: choose one annotation source per physical locus or curate disjoint targets; source views must not be counted twice",
        ));
    }
    let (index, count) = build_index(&panel, &request.policy)?;
    let output = PathBuf::from(&request.output_dir);
    let parent = output
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or_else(|| Path::new("."));
    let staging = tempfile::Builder::new()
        .prefix(".gentle-rna-panel-")
        .tempdir_in(parent)
        .map_err(|e| bad(format!("panel_staging: {e}")))?;
    write_json(&staging.path().join("panel.json"), &panel)?;
    let mut report = RnaReadPanelScreenReport { schema: RNA_PANEL_SCREEN_SCHEMA.into(), content_sha256: String::new(),
        status: "completed".into(), request: request.clone(), panel_sha256: panel.content_sha256.clone(), indexed_postings: count,
        isoform_assignment_performed: false, alignment_confirmation_performed: false,
        files_completed: vec![], files_not_completed: vec![], inputs_not_completed: Some(vec![]), control_error: None, nonclaims: vec![
            "Partial seed-chain support only co-nominates members; it does not establish compatibility with the entire read, exclusive isoform identification, full-length capture, calibrated probability, expression or transcriptome-wide specificity.".into(),
            "A read spanning one isoform's distinguishing junction can still co-nominate an alternative through a shared exon. Junction confirmation remains separate; no conflicting isoform is eliminated by this screen.".into(),
            "Every supported template/orientation is retained. No primer trimming or oligo(dT) orientation flip is applied.".into(),
            "Lookup uppercases bases and converts U to T; quality characters are validated but do not filter or weight nominations. Raw input bytes remain digest-bound.".into(),
            "Junction words are local seed evidence, not splice-aware alignment confirmation. Multi-segment candidates are not proven chimeras or fusions.".into(),
            "Empty support is not biological absence. Thresholds require real-library/background calibration; k=15 is an uncalibrated starting policy.".into(),
            "Input SHA256 binds the compressed bytes actually streamed; supplier MD5 verifies transfer integrity, not authenticated provenance.".into(),
            "Cancelled runs retain only verified completed-file units; retry remaining files in a new output directory, with no in-file resume.".into(),
            "A mid-run marker failure produces control_error, not an authorized cancellation or completed analysis; only checksum-verified completed units are retained. Parse/checksum/resource failures remain hard errors.".into(),
        ] };
    let mut equivalence_bytes = 0usize;
    for (ordinal, source) in request.files.iter().enumerate() {
        let file_id = format!("file_{ordinal:04}");
        let expected_input_digest = source
            .expected_md5
            .as_ref()
            .map(|hash| format!("md5:{hash}"))
            .unwrap_or_else(|| {
                format!(
                    "sha256:{}",
                    source
                        .expected_sha256
                        .as_ref()
                        .unwrap()
                        .trim_start_matches("sha256:")
                        .to_ascii_lowercase()
                )
            });
        let unit = tempfile::Builder::new()
            .prefix("unit-")
            .tempdir_in(staging.path())
            .map_err(|e| bad(e.to_string()))?;
        let mut rows = csv::WriterBuilder::new()
            .delimiter(b'\t')
            .from_writer(HashedWriter::create(&unit.path().join("nominations.tsv"))?);
        rows.write_record([
            "sample_id",
            "barcode",
            "file_id",
            "expected_input_digest",
            "record_index",
            "read_id",
            "decompressed_record_offset",
            "class",
            "candidates_json",
        ])
        .map_err(|e| bad(e.to_string()))?;
        let mut fastas = BTreeMap::new();
        if request.policy.write_locus_fastas {
            for (i, locus) in panel.loci.iter().enumerate() {
                let name = format!("locus_{i:03}.fasta");
                fastas.insert(
                    locus.source.locus_id.clone(),
                    (
                        name.clone(),
                        std::io::BufWriter::new(HashedWriter::create(&unit.path().join(name))?),
                    ),
                );
            }
        }
        let mut receipt = RnaReadPanelFileReceipt {
            file_id: file_id.clone(),
            input: source.clone(),
            observed_md5: String::new(),
            observed_sha256: String::new(),
            compressed_bytes: 0,
            reads: 0,
            classes: BTreeMap::new(),
            equivalence_classes: BTreeMap::new(),
            artifacts: vec![],
        };
        let mut gate = control::PollGate::default();
        let mut marker_error = None;
        let mut stop_requested = |progress: &mut dyn FnMut(OperationProgress) -> bool,
                                  processed_reads: u64,
                                  bytes_read: u64,
                                  bytes_total: u64,
                                  force: bool|
         -> bool {
            if !gate.due(now(), force) {
                return false;
            }
            match control::marker_requested(request.cancel_path.as_deref(), &cancel_ids) {
                Ok(true) => return true,
                Ok(false) => {}
                Err(error) => {
                    marker_error = Some(control::diagnostic(
                        error,
                        request.cancel_path.as_deref(),
                        if bytes_total == 0 {
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
                processed_reads,
                compressed_bytes_read: bytes_read,
                compressed_bytes_total: bytes_total,
                phase: "streaming_unverified".into(),
            }))
        };
        if stop_requested(progress, 0, 0, 0, true)
            || !progress(OperationProgress::Workflow {
                completed: ordinal,
                total: request.files.len(),
            })
        {
            drop(stop_requested);
            report.control_error = marker_error;
            report.status = if report.control_error.is_some() {
                "control_error"
            } else {
                "cancelled"
            }
            .into();
            break;
        }
        let digest = visit_fastq(Path::new(&source.path), &mut |record| {
            if stop_requested(
                progress,
                receipt.reads,
                record.compressed_bytes_read,
                record.compressed_bytes_total,
                receipt.reads == 0,
            ) {
                return Ok(false);
            }
            let nomination = nominate(
                &panel,
                &index,
                &request.policy,
                &record.sequence,
                &mut || {
                    stop_requested(
                        progress,
                        receipt.reads,
                        record.compressed_bytes_read,
                        record.compressed_bytes_total,
                        false,
                    )
                },
            );
            let (class, candidates) = match nomination {
                Err(error) if error.message == "panel_cancelled" => return Ok(false),
                Err(error) => {
                    return Err(bad(format!(
                        "{}; file_id={file_id}; record_index={}; read_id={}; input={}; no result directory published",
                        error.message, receipt.reads, record.read_id, source.path
                    )));
                }
                Ok(value) => value,
            };
            let mut members = candidates
                .iter()
                .map(|c| c.member_id.clone())
                .collect::<BTreeSet<_>>()
                .into_iter()
                .collect::<Vec<_>>();
            members.sort();
            let equivalence = serde_json::to_string(&RnaReadPanelNominationSetKey {
                class: class.clone(),
                member_ids: members.clone(),
            })
            .map_err(|e| bad(e.to_string()))?;
            if !receipt.equivalence_classes.contains_key(&equivalence)
                && receipt.equivalence_classes.len() >= request.policy.max_equivalence_classes
            {
                return Err(bad(format!(
                    "panel_equivalence_budget: no classes silently dropped; file_id={file_id}; record_index={}; read_id={}; no result directory published",
                    receipt.reads, record.read_id
                )));
            }
            if !receipt.equivalence_classes.contains_key(&equivalence) {
                equivalence_bytes += equivalence.len();
                if equivalence_bytes > 16_000_000 {
                    return Err(bad(format!(
                        "panel_equivalence_byte_budget: all file-unit keys together exceed 16 MB; file_id={file_id}; record_index={}; read_id={}; no result directory published",
                        receipt.reads, record.read_id
                    )));
                }
            }
            *receipt.equivalence_classes.entry(equivalence).or_default() += 1;
            *receipt.classes.entry(class.clone()).or_default() += 1;
            rows.write_record([
                source.sample_id.clone(),
                source.barcode.clone(),
                file_id.clone(),
                expected_input_digest.clone(),
                receipt.reads.to_string(),
                record.read_id.clone(),
                record.logical_offset.to_string(),
                class.clone(),
                serde_json::to_string(&candidates).map_err(|e| bad(e.to_string()))?,
            ])
            .map_err(|e| bad(e.to_string()))?;
            for locus in candidates
                .iter()
                .map(|c| &c.locus_id)
                .collect::<BTreeSet<_>>()
            {
                if let Some((_, writer)) = fastas.get_mut(locus) {
                    writeln!(
                        writer,
                        ">{}_{file_id}_{}_{} {} panel_class={} ambiguous={} panel_members={}\n{}",
                        source.sample_id,
                        expected_input_digest,
                        receipt.reads,
                        record.read_id,
                        class,
                        members.len() > 1,
                        serde_json::to_string(&members).map_err(|e| bad(e.to_string()))?,
                        String::from_utf8_lossy(&record.sequence)
                    )
                    .map_err(|e| bad(e.to_string()))?;
                }
            }
            receipt.reads += 1;
            Ok(true)
        })?;
        drop(stop_requested);
        let Some(digest) = digest else {
            report.control_error = marker_error;
            report.status = if report.control_error.is_some() {
                "control_error"
            } else {
                "cancelled"
            }
            .into();
            break;
        };
        if source
            .expected_md5
            .as_ref()
            .is_some_and(|expected| expected != &digest.md5)
            || source.expected_sha256.as_ref().is_some_and(|expected| {
                expected.trim_start_matches("sha256:").to_ascii_lowercase()
                    != digest.sha256.trim_start_matches("sha256:")
            })
        {
            return Err(bad(format!(
                "panel_checksum_mismatch: {}: no result directory published",
                source.path
            )));
        }
        receipt.observed_sha256 = digest.sha256;
        receipt.observed_md5 = digest.md5;
        receipt.compressed_bytes = digest.bytes;
        if !progress(OperationProgress::RnaReadPanel(RnaReadPanelProgress {
            file_id: file_id.clone(),
            sample_id: source.sample_id.clone(),
            file_ordinal: ordinal + 1,
            files_total: request.files.len(),
            processed_reads: receipt.reads,
            compressed_bytes_read: digest.bytes,
            compressed_bytes_total: digest.bytes,
            phase: "file_checksum_verified".into(),
        })) {
            report.status = "cancelled".into();
        }
        let writer = rows.into_inner().map_err(|e| bad(e.to_string()))?;
        receipt
            .artifacts
            .push(writer.receipt(format!("artifacts/{file_id}/nominations.tsv"))?);
        for (locus, (name, mut writer)) in fastas {
            writer.flush().map_err(|e| bad(e.to_string()))?;
            let writer = writer.into_inner().map_err(|e| bad(e.to_string()))?;
            let mut artifact = writer.receipt(format!("artifacts/{file_id}/{name}"))?;
            artifact.locus_id = Some(locus);
            receipt.artifacts.push(artifact);
        }
        fs::rename(unit.path(), staging.path().join(&file_id)).map_err(|e| bad(e.to_string()))?;
        report.files_completed.push(receipt);
        if report.status == "cancelled" {
            break;
        }
    }
    report.files_not_completed = request
        .files
        .iter()
        .skip(report.files_completed.len())
        .map(|i| i.path.clone())
        .collect();
    report.inputs_not_completed = Some(
        request
            .files
            .iter()
            .skip(report.files_completed.len())
            .cloned()
            .collect(),
    );
    if report.status == "completed" {
        match control::marker_requested(request.cancel_path.as_deref(), &cancel_ids) {
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
    report.content_sha256 = sha256_prefixed_bytes(&encode(&report)?);
    // Reserve a fresh destination atomically. A report is the final commit marker;
    // no existing directory, input or prior run can be overwritten.
    fs::create_dir(&output).map_err(|e| bad(format!("panel_output_reserve: {e}")))?;
    if let Err(error) = fs::rename(staging.path(), output.join("artifacts")) {
        let _ = fs::remove_dir(&output);
        return Err(bad(format!("panel_output_publish: {error}")));
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

#[cfg(test)]
#[path = "rna_panel/summary_tests.rs"]
mod summary_tests;
#[cfg(test)]
#[path = "rna_panel/tests.rs"]
mod tests;
