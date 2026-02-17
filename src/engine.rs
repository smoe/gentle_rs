use crate::{
    app::GENtleApp,
    dna_sequence::DNAsequence,
    iupac_code::IupacCode,
    methylation_sites::MethylationMode,
    ENZYMES,
};
use serde::{Deserialize, Serialize};
use std::{
    collections::{HashMap, HashSet},
    error::Error,
    fmt,
    fs::File,
    io::Write,
    path::Path,
};

pub type SeqId = String;
pub type OpId = String;
pub type RunId = String;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum DisplayTarget {
    SequencePanel,
    MapPanel,
    Features,
    RestrictionEnzymes,
    GcContents,
    OpenReadingFrames,
    MethylationSites,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default)]
pub struct DisplaySettings {
    pub show_sequence_panel: bool,
    pub show_map_panel: bool,
    pub show_features: bool,
    pub show_restriction_enzymes: bool,
    pub show_gc_contents: bool,
    pub show_open_reading_frames: bool,
    pub show_methylation_sites: bool,
}

impl Default for DisplaySettings {
    fn default() -> Self {
        Self {
            show_sequence_panel: true,
            show_map_panel: true,
            show_features: true,
            show_restriction_enzymes: true,
            show_gc_contents: true,
            show_open_reading_frames: true,
            show_methylation_sites: false,
        }
    }
}

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct ProjectState {
    pub sequences: HashMap<SeqId, DNAsequence>,
    pub metadata: HashMap<String, serde_json::Value>,
    #[serde(default)]
    pub display: DisplaySettings,
}

impl ProjectState {
    pub fn load_from_path(path: &str) -> Result<Self, EngineError> {
        let text = std::fs::read_to_string(path).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!("Could not read state file '{path}': {e}"),
        })?;
        serde_json::from_str(&text).map_err(|e| EngineError {
            code: ErrorCode::InvalidInput,
            message: format!("Could not parse state JSON '{path}': {e}"),
        })
    }

    pub fn save_to_path(&self, path: &str) -> Result<(), EngineError> {
        let text = serde_json::to_string_pretty(self).map_err(|e| EngineError {
            code: ErrorCode::Internal,
            message: format!("Could not serialize state: {e}"),
        })?;
        std::fs::write(path, text).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!("Could not write state file '{path}': {e}"),
        })
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum ExportFormat {
    GenBank,
    Fasta,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum Operation {
    LoadFile {
        path: String,
        as_id: Option<SeqId>,
    },
    SaveFile {
        seq_id: SeqId,
        path: String,
        format: ExportFormat,
    },
    Digest {
        input: SeqId,
        enzymes: Vec<String>,
        output_prefix: Option<String>,
    },
    Ligation {
        inputs: Vec<SeqId>,
        circularize_if_possible: bool,
        output_id: Option<SeqId>,
    },
    Pcr {
        template: SeqId,
        forward_primer: String,
        reverse_primer: String,
        output_id: Option<SeqId>,
    },
    ExtractRegion {
        input: SeqId,
        from: usize,
        to: usize,
        output_id: Option<SeqId>,
    },
    SetDisplayVisibility {
        target: DisplayTarget,
        visible: bool,
    },
    SetTopology {
        seq_id: SeqId,
        circular: bool,
    },
    RecomputeFeatures {
        seq_id: SeqId,
    },
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Workflow {
    pub run_id: RunId,
    pub ops: Vec<Operation>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OpResult {
    pub op_id: OpId,
    pub created_seq_ids: Vec<SeqId>,
    pub changed_seq_ids: Vec<SeqId>,
    pub warnings: Vec<String>,
    pub messages: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OperationRecord {
    pub run_id: RunId,
    pub op: Operation,
    pub result: OpResult,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub enum ErrorCode {
    InvalidInput,
    NotFound,
    Unsupported,
    Io,
    Internal,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EngineError {
    pub code: ErrorCode,
    pub message: String,
}

impl fmt::Display for EngineError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{:?}: {}", self.code, self.message)
    }
}

impl Error for EngineError {}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Capabilities {
    pub protocol_version: String,
    pub supported_operations: Vec<String>,
    pub supported_export_formats: Vec<String>,
    pub deterministic_operation_log: bool,
}

pub trait Engine {
    fn apply(&mut self, op: Operation) -> Result<OpResult, EngineError>;
    fn apply_workflow(&mut self, wf: Workflow) -> Result<Vec<OpResult>, EngineError>;
    fn snapshot(&self) -> &ProjectState;
}

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct GentleEngine {
    state: ProjectState,
    journal: Vec<OperationRecord>,
    op_counter: u64,
}

impl GentleEngine {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn from_state(state: ProjectState) -> Self {
        Self {
            state,
            ..Self::default()
        }
    }

    pub fn state(&self) -> &ProjectState {
        &self.state
    }

    pub fn state_mut(&mut self) -> &mut ProjectState {
        &mut self.state
    }

    pub fn capabilities() -> Capabilities {
        Capabilities {
            protocol_version: "v1".to_string(),
            supported_operations: vec![
                "LoadFile".to_string(),
                "SaveFile".to_string(),
                "Digest".to_string(),
                "Ligation".to_string(),
                "Pcr".to_string(),
                "ExtractRegion".to_string(),
                "SetDisplayVisibility".to_string(),
                "SetTopology".to_string(),
                "RecomputeFeatures".to_string(),
            ],
            supported_export_formats: vec!["GenBank".to_string(), "Fasta".to_string()],
            deterministic_operation_log: true,
        }
    }

    pub fn operation_log(&self) -> &[OperationRecord] {
        &self.journal
    }

    fn next_op_id(&mut self) -> OpId {
        self.op_counter += 1;
        format!("op-{}", self.op_counter)
    }

    fn derive_seq_id(path: &str) -> SeqId {
        Path::new(path)
            .file_stem()
            .map(|s| s.to_string_lossy().to_string())
            .filter(|s| !s.is_empty())
            .unwrap_or_else(|| "sequence".to_string())
    }

    fn unique_seq_id(&self, base: &str) -> SeqId {
        if !self.state.sequences.contains_key(base) {
            return base.to_string();
        }
        let mut i = 2usize;
        loop {
            let candidate = format!("{base}_{i}");
            if !self.state.sequences.contains_key(&candidate) {
                return candidate;
            }
            i += 1;
        }
    }

    fn prepare_sequence(dna: &mut DNAsequence) {
        ENZYMES
            .restriction_enzymes()
            .clone_into(dna.restriction_enzymes_mut());
        dna.set_max_restriction_enzyme_sites(Some(2));
        dna.set_methylation_mode(MethylationMode::both());
        dna.update_computed_features();
    }

    fn save_as_fasta(seq_id: &str, dna: &DNAsequence, path: &str) -> Result<(), EngineError> {
        let mut file = File::create(path).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!("Could not create FASTA file '{path}': {e}"),
        })?;

        let header = dna
            .name()
            .clone()
            .unwrap_or_else(|| seq_id.to_string())
            .replace(' ', "_");
        let seq = dna.get_forward_string();

        writeln!(file, ">{header}").map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!("Could not write FASTA header to '{path}': {e}"),
        })?;

        for chunk in seq.as_bytes().chunks(80) {
            file.write_all(chunk).map_err(|e| EngineError {
                code: ErrorCode::Io,
                message: format!("Could not write FASTA sequence to '{path}': {e}"),
            })?;
            file.write_all(b"\n").map_err(|e| EngineError {
                code: ErrorCode::Io,
                message: format!("Could not write FASTA newline to '{path}': {e}"),
            })?;
        }

        Ok(())
    }

    fn reverse_complement(seq: &str) -> String {
        seq.as_bytes()
            .iter()
            .rev()
            .map(|c| IupacCode::letter_complement(*c))
            .map(char::from)
            .collect()
    }

    fn normalize_dna_text(seq: &str) -> String {
        let cleaned = DNAsequence::validate_dna_sequence(seq.as_bytes());
        String::from_utf8_lossy(&cleaned).to_string()
    }

    fn find_subsequence(haystack: &[u8], needle: &[u8], start: usize) -> Option<usize> {
        if needle.is_empty() || haystack.len() < needle.len() || start >= haystack.len() {
            return None;
        }
        let end = haystack.len() - needle.len();
        (start..=end).find(|idx| &haystack[*idx..*idx + needle.len()] == needle)
    }

    fn apply_internal(&mut self, op: Operation) -> Result<OpResult, EngineError> {
        let op_id = self.next_op_id();
        let mut result = OpResult {
            op_id,
            created_seq_ids: vec![],
            changed_seq_ids: vec![],
            warnings: vec![],
            messages: vec![],
        };

        match op {
            Operation::LoadFile { path, as_id } => {
                let mut dna = GENtleApp::load_from_file(&path).map_err(|e| EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!("Could not load sequence file '{path}': {e}"),
                })?;
                Self::prepare_sequence(&mut dna);

                let base = as_id.unwrap_or_else(|| Self::derive_seq_id(&path));
                let seq_id = self.unique_seq_id(&base);
                self.state.sequences.insert(seq_id.clone(), dna);
                result.created_seq_ids.push(seq_id.clone());
                result.messages.push(format!("Loaded '{path}' as '{seq_id}'"));
            }
            Operation::SaveFile {
                seq_id,
                path,
                format,
            } => {
                let dna = self
                    .state
                    .sequences
                    .get(&seq_id)
                    .ok_or_else(|| EngineError {
                        code: ErrorCode::NotFound,
                        message: format!("Sequence '{seq_id}' not found"),
                    })?;

                match format {
                    ExportFormat::GenBank => {
                        dna.write_genbank_file(&path).map_err(|e| EngineError {
                            code: ErrorCode::Io,
                            message: format!("Could not write GenBank file '{path}': {e}"),
                        })?;
                    }
                    ExportFormat::Fasta => Self::save_as_fasta(&seq_id, dna, &path)?,
                }

                result.changed_seq_ids.push(seq_id.clone());
                result.messages.push(format!("Wrote '{seq_id}' to '{path}'"));
            }
            Operation::Digest {
                input,
                enzymes,
                output_prefix,
            } => {
                let dna = self
                    .state
                    .sequences
                    .get(&input)
                    .ok_or_else(|| EngineError {
                        code: ErrorCode::NotFound,
                        message: format!("Sequence '{input}' not found"),
                    })?
                    .clone();

                if enzymes.is_empty() {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: "Digest requires at least one enzyme".to_string(),
                    });
                }

                let names_ref: Vec<&str> = enzymes.iter().map(|s| s.as_str()).collect();
                let found = ENZYMES.restriction_enzymes_by_name(&names_ref);
                if found.is_empty() {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "None of the requested enzymes are known: {}",
                            enzymes.join(",")
                        ),
                    });
                }

                let found_names: HashSet<String> = found.iter().map(|e| e.name.clone()).collect();
                let missing: Vec<String> = enzymes
                    .iter()
                    .filter(|name| !found_names.contains(*name))
                    .cloned()
                    .collect();
                if !missing.is_empty() {
                    result.warnings.push(format!(
                        "Unknown enzymes ignored: {}",
                        missing.join(",")
                    ));
                }

                let fragments = dna.restriction_enzymes_full_digest(found);
                let prefix = output_prefix.unwrap_or_else(|| format!("{input}_digest"));

                for (i, mut fragment) in fragments.into_iter().enumerate() {
                    Self::prepare_sequence(&mut fragment);
                    let candidate = format!("{}_{}", prefix, i + 1);
                    let seq_id = self.unique_seq_id(&candidate);
                    self.state.sequences.insert(seq_id.clone(), fragment);
                    result.created_seq_ids.push(seq_id);
                }

                result.messages.push(format!(
                    "Digest created {} fragment(s)",
                    result.created_seq_ids.len()
                ));
            }
            Operation::Ligation {
                inputs,
                circularize_if_possible,
                output_id,
            } => {
                if inputs.is_empty() {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: "Ligation requires at least one input sequence".to_string(),
                    });
                }

                let mut merged = String::new();
                for seq_id in &inputs {
                    let dna = self
                        .state
                        .sequences
                        .get(seq_id)
                        .ok_or_else(|| EngineError {
                            code: ErrorCode::NotFound,
                            message: format!("Sequence '{seq_id}' not found"),
                        })?;
                    merged.push_str(&dna.get_forward_string());
                }

                let mut product = DNAsequence::from_sequence(&merged).map_err(|e| EngineError {
                    code: ErrorCode::Internal,
                    message: format!("Could not create ligation product: {e}"),
                })?;
                if circularize_if_possible {
                    product.set_circular(true);
                    result.warnings.push(
                        "Ligation currently concatenates inputs without sticky-end compatibility checks"
                            .to_string(),
                    );
                } else {
                    product.set_circular(false);
                }
                Self::prepare_sequence(&mut product);

                let base = output_id.unwrap_or_else(|| "ligation_product".to_string());
                let seq_id = self.unique_seq_id(&base);
                self.state.sequences.insert(seq_id.clone(), product);
                result.created_seq_ids.push(seq_id.clone());
                result.messages.push(format!(
                    "Ligation created '{}' from {} input sequence(s)",
                    seq_id,
                    inputs.len()
                ));
            }
            Operation::Pcr {
                template,
                forward_primer,
                reverse_primer,
                output_id,
            } => {
                let dna = self
                    .state
                    .sequences
                    .get(&template)
                    .ok_or_else(|| EngineError {
                        code: ErrorCode::NotFound,
                        message: format!("Sequence '{template}' not found"),
                    })?
                    .clone();

                if dna.is_circular() {
                    return Err(EngineError {
                        code: ErrorCode::Unsupported,
                        message: "PCR on circular templates is not implemented yet".to_string(),
                    });
                }

                let template_seq = dna.get_forward_string().to_ascii_uppercase();
                let template_bytes = template_seq.as_bytes();
                let fwd = Self::normalize_dna_text(&forward_primer);
                let rev = Self::normalize_dna_text(&reverse_primer);
                if fwd.is_empty() || rev.is_empty() {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: "PCR primers must not be empty".to_string(),
                    });
                }

                let rev_binding = Self::reverse_complement(&rev);
                let fwd_pos =
                    Self::find_subsequence(template_bytes, fwd.as_bytes(), 0).ok_or_else(|| {
                        EngineError {
                            code: ErrorCode::InvalidInput,
                            message: "Forward primer not found on template".to_string(),
                        }
                    })?;
                let search_from = fwd_pos + fwd.len();
                let rev_pos = Self::find_subsequence(template_bytes, rev_binding.as_bytes(), search_from)
                    .ok_or_else(|| EngineError {
                        code: ErrorCode::InvalidInput,
                        message: "Reverse primer binding site not found downstream of forward primer"
                            .to_string(),
                    })?;

                let amplicon_end = rev_pos + rev_binding.len();
                let amplicon = &template_seq[fwd_pos..amplicon_end];
                let mut pcr_product =
                    DNAsequence::from_sequence(amplicon).map_err(|e| EngineError {
                        code: ErrorCode::Internal,
                        message: format!("Could not create PCR product: {e}"),
                    })?;
                pcr_product.set_circular(false);
                Self::prepare_sequence(&mut pcr_product);

                let base = output_id.unwrap_or_else(|| format!("{template}_pcr"));
                let seq_id = self.unique_seq_id(&base);
                self.state.sequences.insert(seq_id.clone(), pcr_product);
                result.created_seq_ids.push(seq_id.clone());
                result.messages.push(format!(
                    "PCR product '{}' created: {}..{} (len {})",
                    seq_id,
                    fwd_pos,
                    amplicon_end,
                    amplicon.len()
                ));
            }
            Operation::ExtractRegion {
                input,
                from,
                to,
                output_id,
            } => {
                let dna = self
                    .state
                    .sequences
                    .get(&input)
                    .ok_or_else(|| EngineError {
                        code: ErrorCode::NotFound,
                        message: format!("Sequence '{input}' not found"),
                    })?
                    .clone();

                if from == to {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: "ExtractRegion requires from != to".to_string(),
                    });
                }
                let fragment = dna.get_range_safe(from..to).ok_or_else(|| EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!(
                        "Could not extract region {}..{} from sequence '{}'",
                        from, to, input
                    ),
                })?;
                let text = String::from_utf8_lossy(&fragment).to_string();
                let mut out = DNAsequence::from_sequence(&text).map_err(|e| EngineError {
                    code: ErrorCode::Internal,
                    message: format!("Could not create extracted sequence: {e}"),
                })?;
                out.set_circular(false);
                Self::prepare_sequence(&mut out);

                let base = output_id.unwrap_or_else(|| format!("{input}_region"));
                let seq_id = self.unique_seq_id(&base);
                self.state.sequences.insert(seq_id.clone(), out);
                result.created_seq_ids.push(seq_id.clone());
                result.messages.push(format!(
                    "Extracted region {}..{} from '{}' into '{}'",
                    from, to, input, seq_id
                ));
            }
            Operation::SetDisplayVisibility { target, visible } => {
                let (name, slot): (&str, &mut bool) = match target {
                    DisplayTarget::SequencePanel => {
                        ("sequence_panel", &mut self.state.display.show_sequence_panel)
                    }
                    DisplayTarget::MapPanel => ("map_panel", &mut self.state.display.show_map_panel),
                    DisplayTarget::Features => ("features", &mut self.state.display.show_features),
                    DisplayTarget::RestrictionEnzymes => (
                        "restriction_enzymes",
                        &mut self.state.display.show_restriction_enzymes,
                    ),
                    DisplayTarget::GcContents => {
                        ("gc_contents", &mut self.state.display.show_gc_contents)
                    }
                    DisplayTarget::OpenReadingFrames => (
                        "open_reading_frames",
                        &mut self.state.display.show_open_reading_frames,
                    ),
                    DisplayTarget::MethylationSites => (
                        "methylation_sites",
                        &mut self.state.display.show_methylation_sites,
                    ),
                };
                *slot = visible;
                result
                    .messages
                    .push(format!("Set display target '{name}' to {visible}"));
            }
            Operation::SetTopology { seq_id, circular } => {
                let dna = self
                    .state
                    .sequences
                    .get_mut(&seq_id)
                    .ok_or_else(|| EngineError {
                        code: ErrorCode::NotFound,
                        message: format!("Sequence '{seq_id}' not found"),
                    })?;
                dna.set_circular(circular);
                dna.update_computed_features();
                result.changed_seq_ids.push(seq_id.clone());
                result.messages.push(format!(
                    "Set topology of '{seq_id}' to {}",
                    if circular { "circular" } else { "linear" }
                ));
            }
            Operation::RecomputeFeatures { seq_id } => {
                let dna = self
                    .state
                    .sequences
                    .get_mut(&seq_id)
                    .ok_or_else(|| EngineError {
                        code: ErrorCode::NotFound,
                        message: format!("Sequence '{seq_id}' not found"),
                    })?;
                dna.update_computed_features();
                result.changed_seq_ids.push(seq_id.clone());
                result
                    .messages
                    .push(format!("Recomputed features for '{seq_id}'"));
            }
        }

        Ok(result)
    }
}

impl Engine for GentleEngine {
    fn apply(&mut self, op: Operation) -> Result<OpResult, EngineError> {
        let run_id = "interactive".to_string();
        let result = self.apply_internal(op.clone())?;
        self.journal.push(OperationRecord {
            run_id,
            op,
            result: result.clone(),
        });
        Ok(result)
    }

    fn apply_workflow(&mut self, wf: Workflow) -> Result<Vec<OpResult>, EngineError> {
        let mut results = Vec::new();
        for op in &wf.ops {
            let result = self.apply_internal(op.clone())?;
            self.journal.push(OperationRecord {
                run_id: wf.run_id.clone(),
                op: op.clone(),
                result: result.clone(),
            });
            results.push(result);
        }
        Ok(results)
    }

    fn snapshot(&self) -> &ProjectState {
        &self.state
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn seq(s: &str) -> DNAsequence {
        DNAsequence::from_sequence(s).unwrap()
    }

    #[test]
    fn test_set_display_visibility() {
        let mut engine = GentleEngine::new();
        let res = engine
            .apply(Operation::SetDisplayVisibility {
                target: DisplayTarget::Features,
                visible: false,
            })
            .unwrap();
        assert!(res.messages.iter().any(|m| m.contains("features")));
        assert!(!engine.state().display.show_features);
    }

    #[test]
    fn test_extract_region() {
        let mut state = ProjectState::default();
        state
            .sequences
            .insert("x".to_string(), seq(&"ATGC".repeat(100)));
        let mut engine = GentleEngine::from_state(state);
        let res = engine
            .apply(Operation::ExtractRegion {
                input: "x".to_string(),
                from: 2,
                to: 7,
                output_id: Some("part".to_string()),
            })
            .unwrap();
        assert_eq!(res.created_seq_ids, vec!["part".to_string()]);
        assert_eq!(
            engine.state().sequences.get("part").unwrap().get_forward_string(),
            "GCATG"
        );
    }

    #[test]
    fn test_ligation_simple_concatenation() {
        let mut state = ProjectState::default();
        state
            .sequences
            .insert("a".to_string(), seq(&"ATGC".repeat(40)));
        state
            .sequences
            .insert("b".to_string(), seq(&"TTAA".repeat(40)));
        let mut engine = GentleEngine::from_state(state);
        let res = engine
            .apply(Operation::Ligation {
                inputs: vec!["a".to_string(), "b".to_string()],
                circularize_if_possible: false,
                output_id: Some("ab".to_string()),
            })
            .unwrap();
        assert_eq!(res.created_seq_ids, vec!["ab".to_string()]);
        let out = engine.state().sequences.get("ab").unwrap();
        assert_eq!(
            out.get_forward_string(),
            format!("{}{}", "ATGC".repeat(40), "TTAA".repeat(40))
        );
        assert!(!out.is_circular());
    }
}
