use crate::{
    app::GENtleApp, dna_sequence::DNAsequence, iupac_code::IupacCode,
    methylation_sites::MethylationMode, restriction_enzyme::RestrictionEnzyme, ENZYMES,
};
use serde::{Deserialize, Serialize};
use std::{
    collections::hash_map::DefaultHasher,
    collections::{HashMap, HashSet},
    error::Error,
    fmt,
    fs::File,
    hash::{Hash, Hasher},
    io::Write,
    path::Path,
    time::Instant,
};

pub type SeqId = String;
pub type OpId = String;
pub type RunId = String;
pub type NodeId = String;

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

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum SequenceOrigin {
    ImportedGenomic,
    ImportedCdna,
    ImportedSynthetic,
    ImportedUnknown,
    Derived,
    InSilicoSelection,
    Branch,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LineageNode {
    pub node_id: NodeId,
    pub seq_id: SeqId,
    pub created_by_op: Option<OpId>,
    pub origin: SequenceOrigin,
    pub created_at_unix_ms: u128,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LineageEdge {
    pub from_node_id: NodeId,
    pub to_node_id: NodeId,
    pub op_id: OpId,
    pub run_id: RunId,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct LineageGraph {
    pub nodes: HashMap<NodeId, LineageNode>,
    pub seq_to_node: HashMap<SeqId, NodeId>,
    pub edges: Vec<LineageEdge>,
    pub next_node_counter: u64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default)]
pub struct EngineParameters {
    pub max_fragments_per_container: usize,
}

impl Default for EngineParameters {
    fn default() -> Self {
        Self {
            max_fragments_per_container: 80_000,
        }
    }
}

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct ProjectState {
    pub sequences: HashMap<SeqId, DNAsequence>,
    pub metadata: HashMap<String, serde_json::Value>,
    #[serde(default)]
    pub display: DisplaySettings,
    #[serde(default)]
    pub lineage: LineageGraph,
    #[serde(default)]
    pub parameters: EngineParameters,
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
pub enum PrimerLibraryMode {
    Enumerate,
    Sample,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PcrPrimerSpec {
    pub sequence: String,
    pub anneal_len: Option<usize>,
    pub max_mismatches: Option<usize>,
    pub require_3prime_exact_bases: Option<usize>,
    pub library_mode: Option<PrimerLibraryMode>,
    pub max_variants: Option<usize>,
    pub sample_seed: Option<u64>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SnpMutationSpec {
    pub zero_based_position: usize,
    pub reference: String,
    pub alternate: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum LigationProtocol {
    Sticky,
    Blunt,
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
        protocol: LigationProtocol,
        output_prefix: Option<String>,
        unique: Option<bool>,
    },
    MergeContainers {
        inputs: Vec<SeqId>,
        output_prefix: Option<String>,
    },
    Pcr {
        template: SeqId,
        forward_primer: String,
        reverse_primer: String,
        output_id: Option<SeqId>,
        unique: Option<bool>,
    },
    PcrAdvanced {
        template: SeqId,
        forward_primer: PcrPrimerSpec,
        reverse_primer: PcrPrimerSpec,
        output_id: Option<SeqId>,
        unique: Option<bool>,
    },
    PcrMutagenesis {
        template: SeqId,
        forward_primer: PcrPrimerSpec,
        reverse_primer: PcrPrimerSpec,
        mutations: Vec<SnpMutationSpec>,
        output_id: Option<SeqId>,
        unique: Option<bool>,
        require_all_mutations: Option<bool>,
    },
    ExtractRegion {
        input: SeqId,
        from: usize,
        to: usize,
        output_id: Option<SeqId>,
    },
    SelectCandidate {
        input: SeqId,
        criterion: String,
        output_id: Option<SeqId>,
    },
    FilterByMolecularWeight {
        inputs: Vec<SeqId>,
        min_bp: usize,
        max_bp: usize,
        error: f64,
        unique: bool,
        output_prefix: Option<String>,
    },
    Reverse {
        input: SeqId,
        output_id: Option<SeqId>,
    },
    Complement {
        input: SeqId,
        output_id: Option<SeqId>,
    },
    ReverseComplement {
        input: SeqId,
        output_id: Option<SeqId>,
    },
    Branch {
        input: SeqId,
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
    SetParameter {
        name: String,
        value: serde_json::Value,
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
        let mut ret = Self {
            state,
            ..Self::default()
        };
        ret.reconcile_lineage_nodes();
        ret
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
                "MergeContainers".to_string(),
                "Ligation".to_string(),
                "Pcr".to_string(),
                "PcrAdvanced".to_string(),
                "PcrMutagenesis".to_string(),
                "ExtractRegion".to_string(),
                "SelectCandidate".to_string(),
                "FilterByMolecularWeight".to_string(),
                "Reverse".to_string(),
                "Complement".to_string(),
                "ReverseComplement".to_string(),
                "Branch".to_string(),
                "SetDisplayVisibility".to_string(),
                "SetTopology".to_string(),
                "RecomputeFeatures".to_string(),
                "SetParameter".to_string(),
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

    fn now_unix_ms() -> u128 {
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map(|d| d.as_millis())
            .unwrap_or(0)
    }

    fn classify_import_origin(path: &str, dna: &DNAsequence) -> SequenceOrigin {
        let lower = path.to_ascii_lowercase();
        if lower.ends_with(".fa")
            || lower.ends_with(".fasta")
            || lower.ends_with(".fna")
            || lower.ends_with(".ffn")
            || lower.ends_with(".faa")
        {
            // Requested policy: treat all FASTA imports as synthetic for now.
            return SequenceOrigin::ImportedSynthetic;
        }

        // For GenBank/EMBL-like records, use SOURCE/mol_type metadata if present.
        for feature in dna.features() {
            if feature.kind.to_string().to_ascii_uppercase() != "SOURCE" {
                continue;
            }
            for key in ["mol_type", "molecule_type"] {
                if let Some(value) = feature.qualifier_values(key.into()).next() {
                    let v = value.to_ascii_lowercase();
                    if v.contains("synthetic") {
                        return SequenceOrigin::ImportedSynthetic;
                    }
                    if v.contains("cdna") || v.contains("mrna") || v.contains("transcript") {
                        return SequenceOrigin::ImportedCdna;
                    }
                    if v.contains("genomic") {
                        return SequenceOrigin::ImportedGenomic;
                    }
                }
            }
        }

        SequenceOrigin::ImportedUnknown
    }

    fn add_lineage_node(
        &mut self,
        seq_id: &str,
        origin: SequenceOrigin,
        created_by_op: Option<&str>,
    ) -> NodeId {
        self.state.lineage.next_node_counter += 1;
        let node_id = format!("n-{}", self.state.lineage.next_node_counter);
        let node = LineageNode {
            node_id: node_id.clone(),
            seq_id: seq_id.to_string(),
            created_by_op: created_by_op.map(|s| s.to_string()),
            origin,
            created_at_unix_ms: Self::now_unix_ms(),
        };
        self.state
            .lineage
            .seq_to_node
            .insert(seq_id.to_string(), node_id.clone());
        self.state.lineage.nodes.insert(node_id.clone(), node);
        node_id
    }

    fn ensure_lineage_node(&mut self, seq_id: &str) -> NodeId {
        if let Some(node_id) = self.state.lineage.seq_to_node.get(seq_id) {
            return node_id.clone();
        }
        self.add_lineage_node(seq_id, SequenceOrigin::ImportedUnknown, None)
    }

    fn reconcile_lineage_nodes(&mut self) {
        let seq_ids: Vec<String> = self.state.sequences.keys().cloned().collect();
        for seq_id in seq_ids {
            let _ = self.ensure_lineage_node(&seq_id);
        }
    }

    fn add_lineage_edges(
        &mut self,
        parent_seq_ids: &[SeqId],
        created_seq_ids: &[SeqId],
        op_id: &str,
        run_id: &str,
    ) {
        if parent_seq_ids.is_empty() || created_seq_ids.is_empty() {
            return;
        }
        let parent_nodes: Vec<NodeId> = parent_seq_ids
            .iter()
            .map(|seq_id| self.ensure_lineage_node(seq_id))
            .collect();
        let child_nodes: Vec<NodeId> = created_seq_ids
            .iter()
            .map(|seq_id| self.ensure_lineage_node(seq_id))
            .collect();
        for from_node_id in &parent_nodes {
            for to_node_id in &child_nodes {
                self.state.lineage.edges.push(LineageEdge {
                    from_node_id: from_node_id.clone(),
                    to_node_id: to_node_id.clone(),
                    op_id: op_id.to_string(),
                    run_id: run_id.to_string(),
                });
            }
        }
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
        Self::prepare_sequence_light(dna);
        dna.update_computed_features();
    }

    fn prepare_sequence_light(dna: &mut DNAsequence) {
        ENZYMES
            .restriction_enzymes()
            .clone_into(dna.restriction_enzymes_mut());
        dna.set_max_restriction_enzyme_sites(Some(2));
        dna.set_methylation_mode(MethylationMode::both());
    }

    fn digest_with_guard(
        dna: &DNAsequence,
        enzymes: Vec<RestrictionEnzyme>,
        max_fragments: usize,
    ) -> Result<Vec<DNAsequence>, EngineError> {
        let mut fragments = vec![dna.clone()];
        for enzyme in &enzymes {
            println!("Digesting with enzyme: {}", enzyme.name);
            let mut seen_states: HashSet<u64> = HashSet::new();
            let mut rounds: usize = 0;
            let mut last_fragment_count = fragments.len();
            let enzyme_started = Instant::now();
            // Conservative guard against non-converging digest loops.
            let max_rounds = max_fragments.min(1_024).max(64);
            loop {
                rounds += 1;
                if enzyme_started.elapsed().as_millis() > 750 {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "Digest timed out for enzyme '{}' (>{} ms)",
                            enzyme.name, 750
                        ),
                    });
                }
                if rounds > max_rounds {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "Digest did not converge for enzyme '{}' within {} rounds",
                            enzyme.name, max_rounds
                        ),
                    });
                }

                // Detect cyclic/non-converging fragment states.
                let mut state_hasher = DefaultHasher::new();
                fragments.len().hash(&mut state_hasher);
                for seq in &fragments {
                    seq.get_forward_string().hash(&mut state_hasher);
                    seq.overhang().forward_3.hash(&mut state_hasher);
                    seq.overhang().forward_5.hash(&mut state_hasher);
                    seq.overhang().reverse_3.hash(&mut state_hasher);
                    seq.overhang().reverse_5.hash(&mut state_hasher);
                }
                let state_sig = state_hasher.finish();
                if !seen_states.insert(state_sig) {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "Digest entered a repeated state for enzyme '{}'",
                            enzyme.name
                        ),
                    });
                }

                let mut found_one = false;
                let mut new_fragments: Vec<DNAsequence> = vec![];
                for seq in fragments.drain(..) {
                    if let Some(site) = enzyme.get_sites(&seq, None).first() {
                        let split = seq.split_at_restriction_enzyme_site(site);
                        found_one = true;
                        new_fragments.extend(split);
                    } else {
                        new_fragments.push(seq);
                    }

                    if new_fragments.len() > max_fragments {
                        return Err(EngineError {
                            code: ErrorCode::InvalidInput,
                            message: format!(
                                "Digest produced more than max_fragments_per_container={}",
                                max_fragments
                            ),
                        });
                    }
                }
                fragments = new_fragments;
                let current_count = fragments.len();

                if !found_one {
                    println!(
                        "Digest enzyme '{}' completed in {} round(s), fragments: {}",
                        enzyme.name, rounds, current_count
                    );
                    break;
                }
                // For linear-digest progression, total fragment count should increase.
                // If not, we are likely cutting in a pathological cycle.
                if current_count <= last_fragment_count && rounds > 1 {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "Digest stalled for enzyme '{}' (no fragment-count progress)",
                            enzyme.name
                        ),
                    });
                }
                last_fragment_count = current_count;
            }
        }
        Ok(fragments)
    }

    fn max_fragments_per_container(&self) -> usize {
        self.state.parameters.max_fragments_per_container
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

    fn reverse_complement_bytes(seq: &[u8]) -> Vec<u8> {
        seq.iter()
            .rev()
            .map(|c| IupacCode::letter_complement(*c))
            .collect()
    }

    fn complement_bytes(seq: &[u8]) -> Vec<u8> {
        seq.iter()
            .map(|c| IupacCode::letter_complement(*c))
            .collect()
    }

    fn normalize_dna_text(seq: &str) -> String {
        let cleaned = DNAsequence::validate_dna_sequence(seq.as_bytes());
        String::from_utf8_lossy(&cleaned).to_string()
    }

    fn normalize_iupac_text(seq: &str) -> Result<String, EngineError> {
        let upper = seq.trim().to_ascii_uppercase();
        if upper.is_empty() {
            return Ok(upper);
        }
        for b in upper.as_bytes() {
            if !IupacCode::is_valid_letter(*b) {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!("Invalid IUPAC nucleotide '{}'", *b as char),
                });
            }
        }
        Ok(upper)
    }

    fn primer_options(seq: &str) -> Result<Vec<Vec<u8>>, EngineError> {
        let mut ret = Vec::with_capacity(seq.len());
        for b in seq.as_bytes() {
            let opts = IupacCode::from_letter(*b).to_vec();
            if opts.is_empty() {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!("Invalid primer base '{}'", *b as char),
                });
            }
            ret.push(opts);
        }
        Ok(ret)
    }

    fn total_primer_variants(options: &[Vec<u8>]) -> usize {
        options
            .iter()
            .fold(1usize, |acc, v| acc.saturating_mul(v.len().max(1)))
    }

    fn primer_variant_by_index(options: &[Vec<u8>], mut idx: usize) -> String {
        if options.is_empty() {
            return String::new();
        }
        let mut out: Vec<u8> = vec![b'A'; options.len()];
        for pos in (0..options.len()).rev() {
            let radix = options[pos].len();
            let choice = idx % radix;
            out[pos] = options[pos][choice];
            idx /= radix;
        }
        String::from_utf8(out).unwrap_or_default()
    }

    fn expand_primer_variants(
        spec: &PcrPrimerSpec,
        cap: usize,
    ) -> Result<Vec<String>, EngineError> {
        let normalized = Self::normalize_iupac_text(&spec.sequence)?;
        if normalized.is_empty() {
            return Ok(vec![]);
        }
        let options = Self::primer_options(&normalized)?;
        let total = Self::total_primer_variants(&options);
        if total == 0 {
            return Ok(vec![]);
        }

        let max_variants = spec.max_variants.unwrap_or(total).min(cap).max(1);
        let mode = spec
            .library_mode
            .clone()
            .unwrap_or(PrimerLibraryMode::Enumerate);
        match mode {
            PrimerLibraryMode::Enumerate => {
                if total > max_variants {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "Primer variant space ({total}) exceeds max_variants ({max_variants}); use Sample mode or raise limits"
                        ),
                    });
                }
                let mut ret = Vec::with_capacity(total);
                for idx in 0..total {
                    ret.push(Self::primer_variant_by_index(&options, idx));
                }
                Ok(ret)
            }
            PrimerLibraryMode::Sample => {
                let target = max_variants.min(total);
                let mut chosen: Vec<usize> = Vec::with_capacity(target);
                let mut seen: HashSet<usize> = HashSet::with_capacity(target * 2);
                let mut state = spec.sample_seed.unwrap_or(0x9E3779B97F4A7C15);

                while chosen.len() < target {
                    state = state
                        .wrapping_mul(6364136223846793005)
                        .wrapping_add(1442695040888963407);
                    let idx = (state as usize) % total;
                    if seen.insert(idx) {
                        chosen.push(idx);
                    }
                    if seen.len() == total {
                        break;
                    }
                }
                chosen.sort_unstable();
                let mut ret = Vec::with_capacity(chosen.len());
                for idx in chosen {
                    ret.push(Self::primer_variant_by_index(&options, idx));
                }
                Ok(ret)
            }
        }
    }

    fn find_subsequence(haystack: &[u8], needle: &[u8], start: usize) -> Option<usize> {
        if needle.is_empty() || haystack.len() < needle.len() || start >= haystack.len() {
            return None;
        }
        let end = haystack.len() - needle.len();
        (start..=end).find(|idx| &haystack[*idx..*idx + needle.len()] == needle)
    }

    fn find_all_subsequences(haystack: &[u8], needle: &[u8]) -> Vec<usize> {
        let mut ret = vec![];
        if needle.is_empty() || haystack.len() < needle.len() {
            return ret;
        }
        let mut start = 0usize;
        while let Some(pos) = Self::find_subsequence(haystack, needle, start) {
            ret.push(pos);
            start = pos + 1;
            if start >= haystack.len() {
                break;
            }
        }
        ret
    }

    fn right_end_overhangs(dna: &DNAsequence) -> Vec<Vec<u8>> {
        let mut ret = Vec::new();
        if !dna.overhang().forward_3.is_empty() {
            ret.push(dna.overhang().forward_3.clone());
        }
        if !dna.overhang().reverse_5.is_empty() {
            ret.push(dna.overhang().reverse_5.clone());
        }
        ret
    }

    fn left_end_overhangs(dna: &DNAsequence) -> Vec<Vec<u8>> {
        let mut ret = Vec::new();
        if !dna.overhang().forward_5.is_empty() {
            ret.push(dna.overhang().forward_5.clone());
        }
        if !dna.overhang().reverse_3.is_empty() {
            ret.push(dna.overhang().reverse_3.clone());
        }
        ret
    }

    fn right_end_is_blunt(dna: &DNAsequence) -> bool {
        dna.overhang().forward_3.is_empty() && dna.overhang().reverse_5.is_empty()
    }

    fn left_end_is_blunt(dna: &DNAsequence) -> bool {
        dna.overhang().forward_5.is_empty() && dna.overhang().reverse_3.is_empty()
    }

    fn sticky_compatible(left: &DNAsequence, right: &DNAsequence) -> bool {
        let right_opts = Self::right_end_overhangs(left);
        let left_opts = Self::left_end_overhangs(right);
        if right_opts.is_empty() || left_opts.is_empty() {
            return false;
        }
        for r in &right_opts {
            let rc_r = Self::reverse_complement_bytes(r);
            let c_r = Self::complement_bytes(r);
            if left_opts.iter().any(|l| *l == rc_r || *l == c_r) {
                return true;
            }
        }
        false
    }

    fn find_anneal_sites(
        template: &[u8],
        anneal: &[u8],
        max_mismatches: usize,
        require_3prime_exact_bases: usize,
        three_prime_is_window_end: bool,
    ) -> Vec<usize> {
        let mut ret = vec![];
        if anneal.is_empty() || template.len() < anneal.len() {
            return ret;
        }
        let window_len = anneal.len();
        for start in 0..=(template.len() - window_len) {
            let window = &template[start..start + window_len];
            let mismatches = window
                .iter()
                .zip(anneal.iter())
                .filter(|(a, b)| a != b)
                .count();
            if mismatches > max_mismatches {
                continue;
            }
            if require_3prime_exact_bases > 0 {
                if require_3prime_exact_bases > window_len {
                    continue;
                }
                let exact_ok = if three_prime_is_window_end {
                    let from = window_len - require_3prime_exact_bases;
                    window[from..] == anneal[from..]
                } else {
                    window[..require_3prime_exact_bases] == anneal[..require_3prime_exact_bases]
                };
                if !exact_ok {
                    continue;
                }
            }
            ret.push(start);
        }
        ret
    }

    fn apply_internal(&mut self, op: Operation, run_id: &str) -> Result<OpResult, EngineError> {
        self.reconcile_lineage_nodes();
        let op_id = self.next_op_id();
        let mut parent_seq_ids: Vec<SeqId> = vec![];
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
                self.add_lineage_node(
                    &seq_id,
                    Self::classify_import_origin(
                        &path,
                        self.state
                            .sequences
                            .get(&seq_id)
                            .expect("sequence just inserted"),
                    ),
                    Some(&result.op_id),
                );
                result.created_seq_ids.push(seq_id.clone());
                result
                    .messages
                    .push(format!("Loaded '{path}' as '{seq_id}'"));
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
                result
                    .messages
                    .push(format!("Wrote '{seq_id}' to '{path}'"));
            }
            Operation::Digest {
                input,
                enzymes,
                output_prefix,
            } => {
                parent_seq_ids.push(input.clone());
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
                    result
                        .warnings
                        .push(format!("Unknown enzymes ignored: {}", missing.join(",")));
                }

                let fragments =
                    Self::digest_with_guard(&dna, found, self.max_fragments_per_container())?;
                let prefix = output_prefix.unwrap_or_else(|| format!("{input}_digest"));

                for (i, mut fragment) in fragments.into_iter().enumerate() {
                    // Keep digest interactive by deferring expensive feature recomputation.
                    Self::prepare_sequence_light(&mut fragment);
                    let candidate = format!("{}_{}", prefix, i + 1);
                    let seq_id = self.unique_seq_id(&candidate);
                    self.state.sequences.insert(seq_id.clone(), fragment);
                    self.add_lineage_node(&seq_id, SequenceOrigin::Derived, Some(&result.op_id));
                    result.created_seq_ids.push(seq_id);
                }

                result.messages.push(format!(
                    "Digest created {} fragment(s); feature recomputation deferred",
                    result.created_seq_ids.len()
                ));
            }
            Operation::MergeContainers {
                inputs,
                output_prefix,
            } => {
                if inputs.is_empty() {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: "MergeContainers requires at least one input sequence".to_string(),
                    });
                }
                if inputs.len() > self.max_fragments_per_container() {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "MergeContainers input count {} exceeds max_fragments_per_container={}",
                            inputs.len(),
                            self.max_fragments_per_container()
                        ),
                    });
                }
                parent_seq_ids.extend(inputs.clone());
                let prefix = output_prefix.unwrap_or_else(|| "merged".to_string());
                for (i, input) in inputs.iter().enumerate() {
                    let dna = self
                        .state
                        .sequences
                        .get(input)
                        .ok_or_else(|| EngineError {
                            code: ErrorCode::NotFound,
                            message: format!("Sequence '{input}' not found"),
                        })?
                        .clone();
                    let seq_id = self.unique_seq_id(&format!("{}_{}", prefix, i + 1));
                    self.state.sequences.insert(seq_id.clone(), dna);
                    self.add_lineage_node(&seq_id, SequenceOrigin::Derived, Some(&result.op_id));
                    result.created_seq_ids.push(seq_id);
                }
                result.messages.push(format!(
                    "Merged {} input sequence(s) into container prefix '{}'",
                    inputs.len(),
                    prefix
                ));
            }
            Operation::Ligation {
                inputs,
                circularize_if_possible,
                output_id,
                protocol,
                output_prefix,
                unique,
            } => {
                parent_seq_ids.extend(inputs.clone());
                if inputs.len() < 2 {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: "Ligation requires at least two input sequences".to_string(),
                    });
                }
                if 1 > self.max_fragments_per_container() {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "Ligation product count exceeds max_fragments_per_container={}",
                            self.max_fragments_per_container()
                        ),
                    });
                }
                let mut accepted: Vec<(String, String, String)> = vec![];
                for (i, left_id) in inputs.iter().enumerate() {
                    for (j, right_id) in inputs.iter().enumerate() {
                        if i == j {
                            continue;
                        }
                        let left =
                            self.state
                                .sequences
                                .get(left_id)
                                .ok_or_else(|| EngineError {
                                    code: ErrorCode::NotFound,
                                    message: format!("Sequence '{left_id}' not found"),
                                })?;
                        let right =
                            self.state
                                .sequences
                                .get(right_id)
                                .ok_or_else(|| EngineError {
                                    code: ErrorCode::NotFound,
                                    message: format!("Sequence '{right_id}' not found"),
                                })?;

                        let ok = match protocol {
                            LigationProtocol::Sticky => Self::sticky_compatible(left, right),
                            LigationProtocol::Blunt => {
                                Self::right_end_is_blunt(left) && Self::left_end_is_blunt(right)
                            }
                        };
                        if !ok {
                            continue;
                        }

                        let product = format!(
                            "{}{}",
                            left.get_forward_string(),
                            right.get_forward_string()
                        );
                        accepted.push((left_id.clone(), right_id.clone(), product));
                        if accepted.len() > self.max_fragments_per_container() {
                            return Err(EngineError {
                                code: ErrorCode::InvalidInput,
                                message: format!(
                                    "Ligation produced more than max_fragments_per_container={}",
                                    self.max_fragments_per_container()
                                ),
                            });
                        }
                    }
                }

                if accepted.is_empty() {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "No ligation products found for protocol '{:?}'",
                            protocol
                        ),
                    });
                }

                let require_unique = unique.unwrap_or(false);
                if require_unique && accepted.len() != 1 {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "Ligation unique=true requires exactly one product, found {}",
                            accepted.len()
                        ),
                    });
                }
                if output_id.is_some() && accepted.len() != 1 {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: "Ligation output_id can only be used when exactly one product is produced"
                            .to_string(),
                    });
                }

                let prefix = output_prefix.unwrap_or_else(|| "ligation".to_string());
                for (idx, (left_id, right_id, merged)) in accepted.into_iter().enumerate() {
                    let mut product =
                        DNAsequence::from_sequence(&merged).map_err(|e| EngineError {
                            code: ErrorCode::Internal,
                            message: format!("Could not create ligation product: {e}"),
                        })?;
                    product.set_circular(circularize_if_possible);
                    Self::prepare_sequence(&mut product);

                    let seq_id = if idx == 0 {
                        if let Some(id) = output_id.clone() {
                            self.unique_seq_id(&id)
                        } else {
                            self.unique_seq_id(&format!("{}_{}", prefix, idx + 1))
                        }
                    } else {
                        self.unique_seq_id(&format!("{}_{}", prefix, idx + 1))
                    };
                    self.state.sequences.insert(seq_id.clone(), product);
                    self.add_lineage_node(&seq_id, SequenceOrigin::Derived, Some(&result.op_id));
                    result.created_seq_ids.push(seq_id.clone());
                    result.messages.push(format!(
                        "Ligation product '{}' from '{}' + '{}'",
                        seq_id, left_id, right_id
                    ));
                }
            }
            Operation::Pcr {
                template,
                forward_primer,
                reverse_primer,
                output_id,
                unique,
            } => {
                parent_seq_ids.push(template.clone());
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
                let fwd_sites = Self::find_all_subsequences(template_bytes, fwd.as_bytes());
                if fwd_sites.is_empty() {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: "Forward primer not found on template".to_string(),
                    });
                }
                let rev_sites = Self::find_all_subsequences(template_bytes, rev_binding.as_bytes());
                if rev_sites.is_empty() {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: "Reverse primer binding site not found on template".to_string(),
                    });
                }

                let mut amplicon_ranges: Vec<(usize, usize)> = vec![];
                for fwd_pos in &fwd_sites {
                    for rev_pos in &rev_sites {
                        if *rev_pos < *fwd_pos {
                            continue;
                        }
                        let amplicon_end = rev_pos + rev_binding.len();
                        if amplicon_end <= template_seq.len() {
                            amplicon_ranges.push((*fwd_pos, amplicon_end));
                        }
                    }
                }
                amplicon_ranges.sort_unstable();
                amplicon_ranges.dedup();

                if amplicon_ranges.is_empty() {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: "No valid forward/reverse primer pair produced an amplicon"
                            .to_string(),
                    });
                }
                if amplicon_ranges.len() > self.max_fragments_per_container() {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "PCR produced {} amplicons, exceeding max_fragments_per_container={}",
                            amplicon_ranges.len(),
                            self.max_fragments_per_container()
                        ),
                    });
                }

                let require_unique = unique.unwrap_or(false);
                if require_unique && amplicon_ranges.len() != 1 {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "PCR unique=true requires exactly one amplicon, found {}",
                            amplicon_ranges.len()
                        ),
                    });
                }
                if output_id.is_some() && amplicon_ranges.len() != 1 {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message:
                            "PCR output_id can only be used when exactly one amplicon is produced"
                                .to_string(),
                    });
                }

                let default_base = format!("{template}_pcr");
                for (i, (start, end)) in amplicon_ranges.iter().enumerate() {
                    let amplicon = &template_seq[*start..*end];
                    let mut pcr_product =
                        DNAsequence::from_sequence(amplicon).map_err(|e| EngineError {
                            code: ErrorCode::Internal,
                            message: format!("Could not create PCR product: {e}"),
                        })?;
                    pcr_product.set_circular(false);
                    Self::prepare_sequence(&mut pcr_product);

                    let requested = if i == 0 { output_id.clone() } else { None };
                    let seq_id = if let Some(id) = requested {
                        self.unique_seq_id(&id)
                    } else if amplicon_ranges.len() == 1 {
                        self.unique_seq_id(&default_base)
                    } else {
                        self.unique_seq_id(&format!("{}_{}", default_base, i + 1))
                    };
                    self.state.sequences.insert(seq_id.clone(), pcr_product);
                    self.add_lineage_node(&seq_id, SequenceOrigin::Derived, Some(&result.op_id));
                    result.created_seq_ids.push(seq_id.clone());
                    result.messages.push(format!(
                        "PCR product '{}' created: {}..{} (len {})",
                        seq_id,
                        start,
                        end,
                        end - start
                    ));
                }
            }
            Operation::PcrAdvanced {
                template,
                forward_primer,
                reverse_primer,
                output_id,
                unique,
            } => {
                parent_seq_ids.push(template.clone());
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

                let fwd_variants = Self::expand_primer_variants(
                    &forward_primer,
                    self.max_fragments_per_container(),
                )?;
                let rev_variants = Self::expand_primer_variants(
                    &reverse_primer,
                    self.max_fragments_per_container(),
                )?;
                if fwd_variants.is_empty() || rev_variants.is_empty() {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: "PCR primers must not be empty".to_string(),
                    });
                }

                let mut candidates: Vec<(usize, usize, String)> = vec![];
                let mut seen_amplicons: HashSet<String> = HashSet::new();
                for fwd_full in &fwd_variants {
                    let fwd_anneal_len = forward_primer.anneal_len.unwrap_or(fwd_full.len());
                    if fwd_anneal_len == 0 || fwd_anneal_len > fwd_full.len() {
                        return Err(EngineError {
                            code: ErrorCode::InvalidInput,
                            message: "Invalid anneal_len in PCR primer spec".to_string(),
                        });
                    }
                    let fwd_anneal = &fwd_full[fwd_full.len() - fwd_anneal_len..];
                    let fwd_sites = Self::find_anneal_sites(
                        template_bytes,
                        fwd_anneal.as_bytes(),
                        forward_primer.max_mismatches.unwrap_or(0),
                        forward_primer.require_3prime_exact_bases.unwrap_or(0),
                        true,
                    );
                    if fwd_sites.is_empty() {
                        continue;
                    }

                    for rev_full in &rev_variants {
                        let rev_anneal_len = reverse_primer.anneal_len.unwrap_or(rev_full.len());
                        if rev_anneal_len == 0 || rev_anneal_len > rev_full.len() {
                            return Err(EngineError {
                                code: ErrorCode::InvalidInput,
                                message: "Invalid anneal_len in PCR primer spec".to_string(),
                            });
                        }
                        let rev_anneal = &rev_full[rev_full.len() - rev_anneal_len..];
                        let rev_binding = Self::reverse_complement(rev_anneal);
                        let rev_sites = Self::find_anneal_sites(
                            template_bytes,
                            rev_binding.as_bytes(),
                            reverse_primer.max_mismatches.unwrap_or(0),
                            reverse_primer.require_3prime_exact_bases.unwrap_or(0),
                            false,
                        );
                        if rev_sites.is_empty() {
                            continue;
                        }
                        let rev_full_rc = Self::reverse_complement(rev_full);
                        for fwd_pos in &fwd_sites {
                            for rev_pos in &rev_sites {
                                let interior_start = *fwd_pos + fwd_anneal_len;
                                let interior_end = *rev_pos;
                                if interior_start > interior_end {
                                    continue;
                                }
                                let interior = &template_seq[interior_start..interior_end];
                                let amplicon = format!("{fwd_full}{interior}{rev_full_rc}");
                                if seen_amplicons.insert(amplicon.clone()) {
                                    candidates.push((*fwd_pos, *rev_pos, amplicon));
                                    if candidates.len() > self.max_fragments_per_container() {
                                        return Err(EngineError {
                                            code: ErrorCode::InvalidInput,
                                            message: format!(
                                                "PCR produced more than max_fragments_per_container={}",
                                                self.max_fragments_per_container()
                                            ),
                                        });
                                    }
                                }
                            }
                        }
                    }
                }
                if candidates.is_empty() {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: "No valid advanced PCR amplicon could be formed".to_string(),
                    });
                }

                let require_unique = unique.unwrap_or(false);
                if require_unique && candidates.len() != 1 {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "PCR unique=true requires exactly one amplicon, found {}",
                            candidates.len()
                        ),
                    });
                }
                if output_id.is_some() && candidates.len() != 1 {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message:
                            "PCR output_id can only be used when exactly one amplicon is produced"
                                .to_string(),
                    });
                }

                let default_base = format!("{template}_pcr");
                for (i, (fwd_pos, rev_pos, amplicon)) in candidates.into_iter().enumerate() {
                    let mut pcr_product =
                        DNAsequence::from_sequence(&amplicon).map_err(|e| EngineError {
                            code: ErrorCode::Internal,
                            message: format!("Could not create PCR product: {e}"),
                        })?;
                    pcr_product.set_circular(false);
                    Self::prepare_sequence(&mut pcr_product);

                    let requested = if i == 0 { output_id.clone() } else { None };
                    let seq_id = if let Some(id) = requested {
                        self.unique_seq_id(&id)
                    } else if seen_amplicons.len() == 1 {
                        self.unique_seq_id(&default_base)
                    } else {
                        self.unique_seq_id(&format!("{}_{}", default_base, i + 1))
                    };
                    self.state.sequences.insert(seq_id.clone(), pcr_product);
                    self.add_lineage_node(&seq_id, SequenceOrigin::Derived, Some(&result.op_id));
                    result.created_seq_ids.push(seq_id.clone());
                    result.messages.push(format!(
                        "Advanced PCR product '{}' created from fwd@{} rev@{} (len {})",
                        seq_id,
                        fwd_pos,
                        rev_pos,
                        amplicon.len()
                    ));
                }
            }
            Operation::PcrMutagenesis {
                template,
                forward_primer,
                reverse_primer,
                mutations,
                output_id,
                unique,
                require_all_mutations,
            } => {
                parent_seq_ids.push(template.clone());
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
                if mutations.is_empty() {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: "PcrMutagenesis requires at least one mutation".to_string(),
                    });
                }

                let template_seq = dna.get_forward_string().to_ascii_uppercase();
                let template_bytes = template_seq.as_bytes();

                let mut normalized_mutations: Vec<(usize, u8, u8)> = vec![];
                for m in &mutations {
                    if m.zero_based_position >= template_bytes.len() {
                        return Err(EngineError {
                            code: ErrorCode::InvalidInput,
                            message: format!(
                                "Mutation position {} is out of bounds for template length {}",
                                m.zero_based_position,
                                template_bytes.len()
                            ),
                        });
                    }
                    let ref_nt = Self::normalize_dna_text(&m.reference).to_ascii_uppercase();
                    let alt_nt = Self::normalize_dna_text(&m.alternate).to_ascii_uppercase();
                    if ref_nt.len() != 1 || alt_nt.len() != 1 {
                        return Err(EngineError {
                            code: ErrorCode::InvalidInput,
                            message: "Mutation reference/alternate must be single nucleotides"
                                .to_string(),
                        });
                    }
                    let ref_b = ref_nt.as_bytes()[0];
                    let alt_b = alt_nt.as_bytes()[0];
                    if template_bytes[m.zero_based_position] != ref_b {
                        return Err(EngineError {
                            code: ErrorCode::InvalidInput,
                            message: format!(
                                "Mutation reference mismatch at position {}: template has '{}', expected '{}'",
                                m.zero_based_position,
                                template_bytes[m.zero_based_position] as char,
                                ref_b as char
                            ),
                        });
                    }
                    normalized_mutations.push((m.zero_based_position, ref_b, alt_b));
                }

                let fwd_variants = Self::expand_primer_variants(
                    &forward_primer,
                    self.max_fragments_per_container(),
                )?;
                let rev_variants = Self::expand_primer_variants(
                    &reverse_primer,
                    self.max_fragments_per_container(),
                )?;
                if fwd_variants.is_empty() || rev_variants.is_empty() {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: "PCR primers must not be empty".to_string(),
                    });
                }

                let require_all = require_all_mutations.unwrap_or(true);
                let mut selected: Vec<((usize, usize), String)> = vec![];
                let mut seen_amplicons: HashSet<String> = HashSet::new();
                for fwd_full in &fwd_variants {
                    let fwd_anneal_len = forward_primer.anneal_len.unwrap_or(fwd_full.len());
                    if fwd_anneal_len == 0 || fwd_anneal_len > fwd_full.len() {
                        return Err(EngineError {
                            code: ErrorCode::InvalidInput,
                            message: "Invalid anneal_len in PCR primer spec".to_string(),
                        });
                    }
                    let fwd_anneal = &fwd_full[fwd_full.len() - fwd_anneal_len..];
                    let fwd_sites = Self::find_anneal_sites(
                        template_bytes,
                        fwd_anneal.as_bytes(),
                        forward_primer.max_mismatches.unwrap_or(0),
                        forward_primer.require_3prime_exact_bases.unwrap_or(0),
                        true,
                    );
                    if fwd_sites.is_empty() {
                        continue;
                    }

                    for rev_full in &rev_variants {
                        let rev_anneal_len = reverse_primer.anneal_len.unwrap_or(rev_full.len());
                        if rev_anneal_len == 0 || rev_anneal_len > rev_full.len() {
                            return Err(EngineError {
                                code: ErrorCode::InvalidInput,
                                message: "Invalid anneal_len in PCR primer spec".to_string(),
                            });
                        }
                        let rev_anneal = &rev_full[rev_full.len() - rev_anneal_len..];
                        let rev_binding = Self::reverse_complement(rev_anneal);
                        let rev_sites = Self::find_anneal_sites(
                            template_bytes,
                            rev_binding.as_bytes(),
                            reverse_primer.max_mismatches.unwrap_or(0),
                            reverse_primer.require_3prime_exact_bases.unwrap_or(0),
                            false,
                        );
                        if rev_sites.is_empty() {
                            continue;
                        }
                        let rev_full_rc = Self::reverse_complement(rev_full);

                        for fwd_pos in &fwd_sites {
                            for rev_pos in &rev_sites {
                                let interior_start = *fwd_pos + fwd_anneal_len;
                                let interior_end = *rev_pos;
                                if interior_start > interior_end {
                                    continue;
                                }
                                let interior = &template_seq[interior_start..interior_end];
                                let amplicon = format!("{fwd_full}{interior}{rev_full_rc}");
                                let mut introduced_count = 0usize;
                                let mut valid = true;
                                for (pos, _ref_b, alt_b) in &normalized_mutations {
                                    if *pos < *fwd_pos || *pos >= (*rev_pos + rev_anneal_len) {
                                        if require_all {
                                            valid = false;
                                            break;
                                        }
                                        continue;
                                    }

                                    let observed = if *pos < (*fwd_pos + fwd_anneal_len) {
                                        let offset =
                                            fwd_full.len() - fwd_anneal_len + (*pos - *fwd_pos);
                                        fwd_full.as_bytes()[offset]
                                    } else if *pos < *rev_pos {
                                        template_bytes[*pos]
                                    } else {
                                        let offset = *pos - *rev_pos;
                                        rev_full_rc.as_bytes()[offset]
                                    };

                                    if observed == *alt_b {
                                        introduced_count += 1;
                                    } else if require_all {
                                        valid = false;
                                        break;
                                    }
                                }

                                let keep = if require_all {
                                    valid && introduced_count == normalized_mutations.len()
                                } else {
                                    introduced_count > 0
                                };
                                if keep && seen_amplicons.insert(amplicon.clone()) {
                                    selected.push(((*fwd_pos, *rev_pos), amplicon));
                                    if selected.len() > self.max_fragments_per_container() {
                                        return Err(EngineError {
                                            code: ErrorCode::InvalidInput,
                                            message: format!(
                                                "Mutagenesis PCR produced more than max_fragments_per_container={}",
                                                self.max_fragments_per_container()
                                            ),
                                        });
                                    }
                                }
                            }
                        }
                    }
                }

                if selected.is_empty() {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: if require_all {
                            "No amplicon introduced all requested mutations".to_string()
                        } else {
                            "No amplicon introduced any requested mutation".to_string()
                        },
                    });
                }
                if selected.len() > self.max_fragments_per_container() {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "Mutagenesis PCR produced {} amplicons, exceeding max_fragments_per_container={}",
                            selected.len(),
                            self.max_fragments_per_container()
                        ),
                    });
                }

                let require_unique = unique.unwrap_or(false);
                if require_unique && selected.len() != 1 {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "PCR unique=true requires exactly one amplicon, found {}",
                            selected.len()
                        ),
                    });
                }
                if output_id.is_some() && selected.len() != 1 {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message:
                            "PCR output_id can only be used when exactly one amplicon is produced"
                                .to_string(),
                    });
                }

                let default_base = format!("{template}_pcr_mut");
                for (i, ((fwd_pos, rev_pos), amplicon)) in selected.into_iter().enumerate() {
                    let mut pcr_product =
                        DNAsequence::from_sequence(&amplicon).map_err(|e| EngineError {
                            code: ErrorCode::Internal,
                            message: format!("Could not create PCR product: {e}"),
                        })?;
                    pcr_product.set_circular(false);
                    Self::prepare_sequence(&mut pcr_product);

                    let requested = if i == 0 { output_id.clone() } else { None };
                    let seq_id = if let Some(id) = requested {
                        self.unique_seq_id(&id)
                    } else if seen_amplicons.len() == 1 {
                        self.unique_seq_id(&default_base)
                    } else {
                        self.unique_seq_id(&format!("{}_{}", default_base, i + 1))
                    };
                    self.state.sequences.insert(seq_id.clone(), pcr_product);
                    self.add_lineage_node(&seq_id, SequenceOrigin::Derived, Some(&result.op_id));
                    result.created_seq_ids.push(seq_id.clone());
                    result.messages.push(format!(
                        "Mutagenesis PCR product '{}' created from fwd@{} rev@{}",
                        seq_id, fwd_pos, rev_pos
                    ));
                }
            }
            Operation::ExtractRegion {
                input,
                from,
                to,
                output_id,
            } => {
                parent_seq_ids.push(input.clone());
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
                self.add_lineage_node(&seq_id, SequenceOrigin::Derived, Some(&result.op_id));
                result.created_seq_ids.push(seq_id.clone());
                result.messages.push(format!(
                    "Extracted region {}..{} from '{}' into '{}'",
                    from, to, input, seq_id
                ));
            }
            Operation::SelectCandidate {
                input,
                criterion,
                output_id,
            } => {
                parent_seq_ids.push(input.clone());
                let dna = self
                    .state
                    .sequences
                    .get(&input)
                    .ok_or_else(|| EngineError {
                        code: ErrorCode::NotFound,
                        message: format!("Sequence '{input}' not found"),
                    })?
                    .clone();

                let base = output_id.unwrap_or_else(|| format!("{input}_selected"));
                let seq_id = self.unique_seq_id(&base);
                self.state.sequences.insert(seq_id.clone(), dna);
                self.add_lineage_node(
                    &seq_id,
                    SequenceOrigin::InSilicoSelection,
                    Some(&result.op_id),
                );
                result.created_seq_ids.push(seq_id.clone());
                result.warnings.push(
                    "Selection operation is in-silico and may not directly correspond to a unique wet-lab product"
                        .to_string(),
                );
                result.messages.push(format!(
                    "Selected candidate '{}' from '{}' using criterion '{}'",
                    seq_id, input, criterion
                ));
            }
            Operation::FilterByMolecularWeight {
                inputs,
                min_bp,
                max_bp,
                error,
                unique,
                output_prefix,
            } => {
                if inputs.is_empty() {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: "FilterByMolecularWeight requires at least one input sequence"
                            .to_string(),
                    });
                }
                if min_bp > max_bp {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!("min_bp ({min_bp}) must be <= max_bp ({max_bp})"),
                    });
                }
                if !(0.0..=1.0).contains(&error) {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: "error must be between 0.0 and 1.0".to_string(),
                    });
                }

                let min_allowed = ((min_bp as f64) * (1.0 - error)).floor() as usize;
                let max_allowed = ((max_bp as f64) * (1.0 + error)).ceil() as usize;

                let mut matches: Vec<(SeqId, DNAsequence)> = vec![];
                for input in &inputs {
                    parent_seq_ids.push(input.clone());
                    let dna = self
                        .state
                        .sequences
                        .get(input)
                        .ok_or_else(|| EngineError {
                            code: ErrorCode::NotFound,
                            message: format!("Sequence '{input}' not found"),
                        })?
                        .clone();
                    let bp = dna.len();
                    if bp >= min_allowed && bp <= max_allowed {
                        matches.push((input.clone(), dna));
                    }
                }

                if matches.len() > self.max_fragments_per_container() {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "FilterByMolecularWeight produced {} candidates, exceeding max_fragments_per_container={}",
                            matches.len(),
                            self.max_fragments_per_container()
                        ),
                    });
                }

                if unique && matches.len() != 1 {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "unique=true requires exactly one match, found {}",
                            matches.len()
                        ),
                    });
                }

                let prefix = output_prefix.unwrap_or_else(|| "mw_filter".to_string());
                for (i, (_source_id, mut dna)) in matches.into_iter().enumerate() {
                    Self::prepare_sequence(&mut dna);
                    let candidate = format!("{}_{}", prefix, i + 1);
                    let seq_id = self.unique_seq_id(&candidate);
                    self.state.sequences.insert(seq_id.clone(), dna);
                    self.add_lineage_node(
                        &seq_id,
                        SequenceOrigin::InSilicoSelection,
                        Some(&result.op_id),
                    );
                    result.created_seq_ids.push(seq_id);
                }

                result.messages.push(format!(
                    "Molecular-weight filter kept {} sequence(s) in effective range {}-{} bp (requested {}-{} bp, error={:.3})",
                    result.created_seq_ids.len(),
                    min_allowed,
                    max_allowed,
                    min_bp,
                    max_bp,
                    error
                ));
            }
            Operation::Reverse { input, output_id } => {
                parent_seq_ids.push(input.clone());
                let dna = self
                    .state
                    .sequences
                    .get(&input)
                    .ok_or_else(|| EngineError {
                        code: ErrorCode::NotFound,
                        message: format!("Sequence '{input}' not found"),
                    })?
                    .clone();
                let mut text = dna.get_forward_string();
                text = text.chars().rev().collect();
                let mut out = DNAsequence::from_sequence(&text).map_err(|e| EngineError {
                    code: ErrorCode::Internal,
                    message: format!("Could not create reverse sequence: {e}"),
                })?;
                out.set_circular(dna.is_circular());
                Self::prepare_sequence(&mut out);

                let base = output_id.unwrap_or_else(|| format!("{input}_rev"));
                let seq_id = self.unique_seq_id(&base);
                self.state.sequences.insert(seq_id.clone(), out);
                self.add_lineage_node(&seq_id, SequenceOrigin::Derived, Some(&result.op_id));
                result.created_seq_ids.push(seq_id.clone());
                result.messages.push(format!(
                    "Created reverse sequence '{}' from '{}'",
                    seq_id, input
                ));
            }
            Operation::Complement { input, output_id } => {
                parent_seq_ids.push(input.clone());
                let dna = self
                    .state
                    .sequences
                    .get(&input)
                    .ok_or_else(|| EngineError {
                        code: ErrorCode::NotFound,
                        message: format!("Sequence '{input}' not found"),
                    })?
                    .clone();
                let text: String = dna
                    .get_forward_string()
                    .as_bytes()
                    .iter()
                    .map(|c| IupacCode::letter_complement(*c))
                    .map(char::from)
                    .collect();
                let mut out = DNAsequence::from_sequence(&text).map_err(|e| EngineError {
                    code: ErrorCode::Internal,
                    message: format!("Could not create complement sequence: {e}"),
                })?;
                out.set_circular(dna.is_circular());
                Self::prepare_sequence(&mut out);

                let base = output_id.unwrap_or_else(|| format!("{input}_comp"));
                let seq_id = self.unique_seq_id(&base);
                self.state.sequences.insert(seq_id.clone(), out);
                self.add_lineage_node(&seq_id, SequenceOrigin::Derived, Some(&result.op_id));
                result.created_seq_ids.push(seq_id.clone());
                result.messages.push(format!(
                    "Created complement sequence '{}' from '{}'",
                    seq_id, input
                ));
            }
            Operation::ReverseComplement { input, output_id } => {
                parent_seq_ids.push(input.clone());
                let dna = self
                    .state
                    .sequences
                    .get(&input)
                    .ok_or_else(|| EngineError {
                        code: ErrorCode::NotFound,
                        message: format!("Sequence '{input}' not found"),
                    })?
                    .clone();
                let text = Self::reverse_complement(&dna.get_forward_string());
                let mut out = DNAsequence::from_sequence(&text).map_err(|e| EngineError {
                    code: ErrorCode::Internal,
                    message: format!("Could not create reverse-complement sequence: {e}"),
                })?;
                out.set_circular(dna.is_circular());
                Self::prepare_sequence(&mut out);

                let base = output_id.unwrap_or_else(|| format!("{input}_revcomp"));
                let seq_id = self.unique_seq_id(&base);
                self.state.sequences.insert(seq_id.clone(), out);
                self.add_lineage_node(&seq_id, SequenceOrigin::Derived, Some(&result.op_id));
                result.created_seq_ids.push(seq_id.clone());
                result.messages.push(format!(
                    "Created reverse-complement sequence '{}' from '{}'",
                    seq_id, input
                ));
            }
            Operation::Branch { input, output_id } => {
                parent_seq_ids.push(input.clone());
                let dna = self
                    .state
                    .sequences
                    .get(&input)
                    .ok_or_else(|| EngineError {
                        code: ErrorCode::NotFound,
                        message: format!("Sequence '{input}' not found"),
                    })?
                    .clone();

                let base = output_id.unwrap_or_else(|| format!("{input}_branch"));
                let seq_id = self.unique_seq_id(&base);
                self.state.sequences.insert(seq_id.clone(), dna);
                self.add_lineage_node(&seq_id, SequenceOrigin::Branch, Some(&result.op_id));
                result.created_seq_ids.push(seq_id.clone());
                result
                    .messages
                    .push(format!("Branched '{}' into '{}'", input, seq_id));
            }
            Operation::SetDisplayVisibility { target, visible } => {
                let (name, slot): (&str, &mut bool) = match target {
                    DisplayTarget::SequencePanel => (
                        "sequence_panel",
                        &mut self.state.display.show_sequence_panel,
                    ),
                    DisplayTarget::MapPanel => {
                        ("map_panel", &mut self.state.display.show_map_panel)
                    }
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
                let _ = self.ensure_lineage_node(&seq_id);
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
                let _ = self.ensure_lineage_node(&seq_id);
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
            Operation::SetParameter { name, value } => {
                match name.as_str() {
                    "max_fragments_per_container" => {
                        let raw = value.as_u64().ok_or_else(|| EngineError {
                            code: ErrorCode::InvalidInput,
                            message: "SetParameter max_fragments_per_container requires a positive integer".to_string(),
                        })?;
                        if raw == 0 {
                            return Err(EngineError {
                                code: ErrorCode::InvalidInput,
                                message: "max_fragments_per_container must be >= 1".to_string(),
                            });
                        }
                        self.state.parameters.max_fragments_per_container = raw as usize;
                        result.messages.push(format!(
                            "Set parameter '{}' to {}",
                            name, self.state.parameters.max_fragments_per_container
                        ));
                    }
                    _ => {
                        return Err(EngineError {
                            code: ErrorCode::Unsupported,
                            message: format!("Unknown parameter '{}'", name),
                        });
                    }
                }
            }
        }

        self.add_lineage_edges(
            &parent_seq_ids,
            &result.created_seq_ids,
            &result.op_id,
            run_id,
        );

        Ok(result)
    }
}

impl Engine for GentleEngine {
    fn apply(&mut self, op: Operation) -> Result<OpResult, EngineError> {
        let run_id = "interactive".to_string();
        let result = self.apply_internal(op.clone(), &run_id)?;
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
            let result = self.apply_internal(op.clone(), &wf.run_id)?;
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
            engine
                .state()
                .sequences
                .get("part")
                .unwrap()
                .get_forward_string(),
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
                output_id: None,
                protocol: LigationProtocol::Blunt,
                output_prefix: Some("ab".to_string()),
                unique: None,
            })
            .unwrap();
        assert_eq!(res.created_seq_ids.len(), 2);
        let out = engine.state().sequences.get("ab_1").unwrap();
        assert_eq!(
            out.get_forward_string(),
            format!("{}{}", "ATGC".repeat(40), "TTAA".repeat(40))
        );
        assert!(!out.is_circular());
    }

    #[test]
    fn test_merge_containers_creates_pool_copies() {
        let mut state = ProjectState::default();
        state.sequences.insert("a".to_string(), seq("ATGC"));
        state.sequences.insert("b".to_string(), seq("TTAA"));
        let mut engine = GentleEngine::from_state(state);
        let res = engine
            .apply(Operation::MergeContainers {
                inputs: vec!["a".to_string(), "b".to_string()],
                output_prefix: Some("pool".to_string()),
            })
            .unwrap();
        assert_eq!(
            res.created_seq_ids,
            vec!["pool_1".to_string(), "pool_2".to_string()]
        );
        assert_eq!(
            engine
                .state()
                .sequences
                .get("pool_1")
                .unwrap()
                .get_forward_string(),
            "ATGC"
        );
        assert_eq!(
            engine
                .state()
                .sequences
                .get("pool_2")
                .unwrap()
                .get_forward_string(),
            "TTAA"
        );
    }

    #[test]
    fn test_merge_containers_respects_max_fragments() {
        let mut state = ProjectState::default();
        state.parameters.max_fragments_per_container = 1;
        state.sequences.insert("a".to_string(), seq("ATGC"));
        state.sequences.insert("b".to_string(), seq("TTAA"));
        let mut engine = GentleEngine::from_state(state);
        let err = engine
            .apply(Operation::MergeContainers {
                inputs: vec!["a".to_string(), "b".to_string()],
                output_prefix: Some("pool".to_string()),
            })
            .unwrap_err();
        assert!(err.message.contains("max_fragments_per_container"));
    }

    #[test]
    fn test_ligation_protocol_blunt_enumerates_ordered_pairs() {
        let mut state = ProjectState::default();
        state.sequences.insert("a".to_string(), seq("ATGC"));
        state.sequences.insert("b".to_string(), seq("TTAA"));
        let mut engine = GentleEngine::from_state(state);
        let res = engine
            .apply(Operation::Ligation {
                inputs: vec!["a".to_string(), "b".to_string()],
                circularize_if_possible: false,
                output_id: None,
                protocol: LigationProtocol::Blunt,
                output_prefix: Some("lig".to_string()),
                unique: Some(false),
            })
            .unwrap();
        assert_eq!(res.created_seq_ids.len(), 2);
        assert_eq!(
            engine
                .state()
                .sequences
                .get("lig_1")
                .unwrap()
                .get_forward_string(),
            "ATGCTTAA"
        );
        assert_eq!(
            engine
                .state()
                .sequences
                .get("lig_2")
                .unwrap()
                .get_forward_string(),
            "TTAAATGC"
        );
    }

    #[test]
    fn test_ligation_protocol_sticky_uses_overhang_compatibility() {
        let mut state = ProjectState::default();
        state
            .sequences
            .insert("x".to_string(), seq("ATGGATCCGCATGGATCCGCATGGATCCGC"));
        let mut engine = GentleEngine::from_state(state);
        let digest_res = engine
            .apply(Operation::Digest {
                input: "x".to_string(),
                enzymes: vec!["BamHI".to_string()],
                output_prefix: Some("frag".to_string()),
            })
            .unwrap();
        assert!(digest_res.created_seq_ids.len() >= 2);
        let a = digest_res.created_seq_ids[0].clone();
        let b = digest_res.created_seq_ids[1].clone();

        let lig_res = engine
            .apply(Operation::Ligation {
                inputs: vec![a, b],
                circularize_if_possible: false,
                output_id: None,
                protocol: LigationProtocol::Sticky,
                output_prefix: Some("st".to_string()),
                unique: Some(false),
            })
            .unwrap();
        assert!(!lig_res.created_seq_ids.is_empty());
    }

    #[test]
    fn test_workflow_digest_merge_ligation_is_deterministic() {
        let mut base = ProjectState::default();
        base.sequences
            .insert("x".to_string(), seq("ATGGATCCGCATGGATCCGCATGGATCCGC"));

        let run_once = |state: ProjectState| {
            let mut engine = GentleEngine::from_state(state);
            let digest = engine
                .apply(Operation::Digest {
                    input: "x".to_string(),
                    enzymes: vec!["BamHI".to_string(), "EcoRI".to_string()],
                    output_prefix: Some("d".to_string()),
                })
                .unwrap();
            let merge = engine
                .apply(Operation::MergeContainers {
                    inputs: digest.created_seq_ids.clone(),
                    output_prefix: Some("m".to_string()),
                })
                .unwrap();
            let lig = engine
                .apply(Operation::Ligation {
                    inputs: merge.created_seq_ids.clone(),
                    circularize_if_possible: false,
                    output_id: None,
                    protocol: LigationProtocol::Sticky,
                    output_prefix: Some("lig".to_string()),
                    unique: Some(false),
                })
                .unwrap();
            lig.created_seq_ids
        };

        let a = run_once(base.clone());
        let b = run_once(base.clone());
        assert_eq!(a, b);
        assert!(!a.is_empty());
        assert_eq!(a.first().unwrap(), "lig_1");
    }

    #[test]
    fn test_lineage_extract_creates_parent_child_edge() {
        let mut state = ProjectState::default();
        state
            .sequences
            .insert("x".to_string(), seq(&"ATGC".repeat(100)));
        let mut engine = GentleEngine::from_state(state);
        let _ = engine
            .apply(Operation::ExtractRegion {
                input: "x".to_string(),
                from: 2,
                to: 7,
                output_id: Some("part".to_string()),
            })
            .unwrap();

        let lineage = &engine.state().lineage;
        let x_node = lineage.seq_to_node.get("x").unwrap();
        let part_node = lineage.seq_to_node.get("part").unwrap();
        assert!(lineage
            .edges
            .iter()
            .any(|e| e.from_node_id == *x_node && e.to_node_id == *part_node));
    }

    #[test]
    fn test_lineage_ligation_has_two_parents() {
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
                output_id: None,
                protocol: LigationProtocol::Blunt,
                output_prefix: Some("ab".to_string()),
                unique: None,
            })
            .unwrap();

        let lineage = &engine.state().lineage;
        let a_node = lineage.seq_to_node.get("a").unwrap();
        let b_node = lineage.seq_to_node.get("b").unwrap();
        let ab_node = lineage.seq_to_node.get(&res.created_seq_ids[0]).unwrap();
        assert!(lineage
            .edges
            .iter()
            .any(|e| e.from_node_id == *a_node && e.to_node_id == *ab_node));
        assert!(lineage
            .edges
            .iter()
            .any(|e| e.from_node_id == *b_node && e.to_node_id == *ab_node));
    }

    #[test]
    fn test_select_candidate_creates_in_silico_node() {
        let mut state = ProjectState::default();
        state
            .sequences
            .insert("frag".to_string(), seq(&"ATGC".repeat(40)));
        let mut engine = GentleEngine::from_state(state);
        let res = engine
            .apply(Operation::SelectCandidate {
                input: "frag".to_string(),
                criterion: "band_size_range:150-170bp".to_string(),
                output_id: Some("picked".to_string()),
            })
            .unwrap();

        assert_eq!(res.created_seq_ids, vec!["picked".to_string()]);
        assert!(!res.warnings.is_empty());
        let lineage = &engine.state().lineage;
        let picked_node = lineage.seq_to_node.get("picked").unwrap();
        let picked_origin = &lineage.nodes.get(picked_node).unwrap().origin;
        assert!(matches!(picked_origin, SequenceOrigin::InSilicoSelection));
    }

    #[test]
    fn test_reverse_complement_reverse_complement_and_branch() {
        let mut state = ProjectState::default();
        state.sequences.insert("s".to_string(), seq("ATGCCA"));
        let mut engine = GentleEngine::from_state(state);

        let res_rev = engine
            .apply(Operation::Reverse {
                input: "s".to_string(),
                output_id: Some("s_rev".to_string()),
            })
            .unwrap();
        assert_eq!(res_rev.created_seq_ids, vec!["s_rev".to_string()]);
        assert_eq!(
            engine
                .state()
                .sequences
                .get("s_rev")
                .unwrap()
                .get_forward_string(),
            "ACCGTA"
        );

        let res_comp = engine
            .apply(Operation::Complement {
                input: "s".to_string(),
                output_id: Some("s_comp".to_string()),
            })
            .unwrap();
        assert_eq!(res_comp.created_seq_ids, vec!["s_comp".to_string()]);
        assert_eq!(
            engine
                .state()
                .sequences
                .get("s_comp")
                .unwrap()
                .get_forward_string(),
            "TACGGT"
        );

        let res_rc = engine
            .apply(Operation::ReverseComplement {
                input: "s".to_string(),
                output_id: Some("s_rc".to_string()),
            })
            .unwrap();
        assert_eq!(res_rc.created_seq_ids, vec!["s_rc".to_string()]);
        assert_eq!(
            engine
                .state()
                .sequences
                .get("s_rc")
                .unwrap()
                .get_forward_string(),
            "TGGCAT"
        );

        let res_split = engine
            .apply(Operation::Branch {
                input: "s".to_string(),
                output_id: Some("s_branch".to_string()),
            })
            .unwrap();
        assert_eq!(res_split.created_seq_ids, vec!["s_branch".to_string()]);
        assert_eq!(
            engine
                .state()
                .sequences
                .get("s_branch")
                .unwrap()
                .get_forward_string(),
            "ATGCCA"
        );

        let lineage = &engine.state().lineage;
        let s_node = lineage.seq_to_node.get("s").unwrap();
        for derived in ["s_rev", "s_comp", "s_rc", "s_branch"] {
            let dnode = lineage.seq_to_node.get(derived).unwrap();
            assert!(lineage
                .edges
                .iter()
                .any(|e| e.from_node_id == *s_node && e.to_node_id == *dnode));
        }
    }

    #[test]
    fn test_default_max_fragments_per_container_is_80000() {
        let state = ProjectState::default();
        assert_eq!(state.parameters.max_fragments_per_container, 80_000);
    }

    #[test]
    fn test_set_parameter_max_fragments_per_container() {
        let mut engine = GentleEngine::new();
        let res = engine
            .apply(Operation::SetParameter {
                name: "max_fragments_per_container".to_string(),
                value: serde_json::json!(1234),
            })
            .unwrap();
        assert!(res
            .messages
            .iter()
            .any(|m| m.contains("max_fragments_per_container")));
        assert_eq!(engine.state().parameters.max_fragments_per_container, 1234);
    }

    #[test]
    fn test_digest_respects_max_fragments_per_container() {
        let mut state = ProjectState::default();
        state
            .sequences
            .insert("x".to_string(), seq("ATGGATCCGCATGGATCCGCATGGATCCGC"));
        state.parameters.max_fragments_per_container = 2;
        let mut engine = GentleEngine::from_state(state);
        let err = engine
            .apply(Operation::Digest {
                input: "x".to_string(),
                enzymes: vec!["BamHI".to_string()],
                output_prefix: Some("frag".to_string()),
            })
            .unwrap_err();
        assert!(err.message.contains("max_fragments_per_container"));
    }

    #[test]
    fn test_filter_by_molecular_weight_with_error_range() {
        let mut state = ProjectState::default();
        state
            .sequences
            .insert("a".to_string(), seq(&"A".repeat(100)));
        state
            .sequences
            .insert("b".to_string(), seq(&"A".repeat(180)));
        state
            .sequences
            .insert("c".to_string(), seq(&"A".repeat(260)));
        let mut engine = GentleEngine::from_state(state);

        let res = engine
            .apply(Operation::FilterByMolecularWeight {
                inputs: vec!["a".to_string(), "b".to_string(), "c".to_string()],
                min_bp: 150,
                max_bp: 200,
                error: 0.10,
                unique: false,
                output_prefix: Some("mw".to_string()),
            })
            .unwrap();

        assert_eq!(res.created_seq_ids.len(), 1);
        let out = engine.state().sequences.get("mw_1").unwrap();
        assert_eq!(out.len(), 180);
    }

    #[test]
    fn test_filter_by_molecular_weight_unique_fails_on_multiple_matches() {
        let mut state = ProjectState::default();
        state
            .sequences
            .insert("a".to_string(), seq(&"A".repeat(100)));
        state
            .sequences
            .insert("b".to_string(), seq(&"A".repeat(105)));
        state
            .sequences
            .insert("c".to_string(), seq(&"A".repeat(200)));
        let mut engine = GentleEngine::from_state(state);

        let err = engine
            .apply(Operation::FilterByMolecularWeight {
                inputs: vec!["a".to_string(), "b".to_string(), "c".to_string()],
                min_bp: 95,
                max_bp: 105,
                error: 0.10,
                unique: true,
                output_prefix: Some("mw".to_string()),
            })
            .unwrap_err();

        assert!(err.message.contains("exactly one match"));
    }

    #[test]
    fn test_pcr_single_amplicon() {
        let mut state = ProjectState::default();
        state
            .sequences
            .insert("tpl".to_string(), seq("ATGAAACCCGGGTTT"));
        let mut engine = GentleEngine::from_state(state);
        let res = engine
            .apply(Operation::Pcr {
                template: "tpl".to_string(),
                forward_primer: "ATGAAA".to_string(),
                reverse_primer: "AAACCC".to_string(),
                output_id: Some("amp".to_string()),
                unique: Some(true),
            })
            .unwrap();
        assert_eq!(res.created_seq_ids, vec!["amp".to_string()]);
        assert_eq!(
            engine
                .state()
                .sequences
                .get("amp")
                .unwrap()
                .get_forward_string(),
            "ATGAAACCCGGGTTT"
        );
    }

    #[test]
    fn test_pcr_unique_fails_on_multiple_amplicons() {
        let mut state = ProjectState::default();
        state
            .sequences
            .insert("tpl".to_string(), seq("AAAACCCCGGGGAAAACCCCGGGG"));
        let mut engine = GentleEngine::from_state(state);
        let err = engine
            .apply(Operation::Pcr {
                template: "tpl".to_string(),
                forward_primer: "AAAA".to_string(),
                reverse_primer: "CCCC".to_string(),
                output_id: None,
                unique: Some(true),
            })
            .unwrap_err();
        assert!(err.message.contains("unique=true"));
    }

    #[test]
    fn test_pcr_multiple_amplicons_without_unique() {
        let mut state = ProjectState::default();
        state
            .sequences
            .insert("tpl".to_string(), seq("AAAACCCCGGGGAAAACCCCGGGG"));
        let mut engine = GentleEngine::from_state(state);
        let res = engine
            .apply(Operation::Pcr {
                template: "tpl".to_string(),
                forward_primer: "AAAA".to_string(),
                reverse_primer: "CCCC".to_string(),
                output_id: None,
                unique: Some(false),
            })
            .unwrap();
        assert_eq!(res.created_seq_ids.len(), 3);
    }

    #[test]
    fn test_pcr_advanced_inserts_5prime_tail() {
        let mut state = ProjectState::default();
        state
            .sequences
            .insert("tpl".to_string(), seq("ATGAAACCCGGGTTT"));
        let mut engine = GentleEngine::from_state(state);
        let res = engine
            .apply(Operation::PcrAdvanced {
                template: "tpl".to_string(),
                forward_primer: PcrPrimerSpec {
                    sequence: "GGATCCATGAAA".to_string(),
                    anneal_len: Some(6),
                    max_mismatches: Some(0),
                    require_3prime_exact_bases: Some(4),
                    library_mode: None,
                    max_variants: None,
                    sample_seed: None,
                },
                reverse_primer: PcrPrimerSpec {
                    sequence: "AAACCC".to_string(),
                    anneal_len: Some(6),
                    max_mismatches: Some(0),
                    require_3prime_exact_bases: Some(4),
                    library_mode: None,
                    max_variants: None,
                    sample_seed: None,
                },
                output_id: Some("amp_adv".to_string()),
                unique: Some(true),
            })
            .unwrap();
        assert_eq!(res.created_seq_ids, vec!["amp_adv".to_string()]);
        assert_eq!(
            engine
                .state()
                .sequences
                .get("amp_adv")
                .unwrap()
                .get_forward_string(),
            "GGATCCATGAAACCCGGGTTT"
        );
    }

    #[test]
    fn test_pcr_advanced_allows_partial_match_and_introduces_mutation() {
        let mut state = ProjectState::default();
        state
            .sequences
            .insert("tpl".to_string(), seq("ATGAAACCCGGGTTT"));
        let mut engine = GentleEngine::from_state(state);
        let res = engine
            .apply(Operation::PcrAdvanced {
                template: "tpl".to_string(),
                forward_primer: PcrPrimerSpec {
                    sequence: "ATCAAA".to_string(),
                    anneal_len: Some(6),
                    max_mismatches: Some(1),
                    require_3prime_exact_bases: Some(3),
                    library_mode: None,
                    max_variants: None,
                    sample_seed: None,
                },
                reverse_primer: PcrPrimerSpec {
                    sequence: "AAACCC".to_string(),
                    anneal_len: Some(6),
                    max_mismatches: Some(0),
                    require_3prime_exact_bases: Some(4),
                    library_mode: None,
                    max_variants: None,
                    sample_seed: None,
                },
                output_id: Some("amp_mut".to_string()),
                unique: Some(true),
            })
            .unwrap();
        assert_eq!(res.created_seq_ids, vec!["amp_mut".to_string()]);
        assert_eq!(
            engine
                .state()
                .sequences
                .get("amp_mut")
                .unwrap()
                .get_forward_string(),
            "ATCAAACCCGGGTTT"
        );
    }

    #[test]
    fn test_pcr_mutagenesis_single_snp_success() {
        let mut state = ProjectState::default();
        state
            .sequences
            .insert("tpl".to_string(), seq("ATGAAACCCGGGTTT"));
        let mut engine = GentleEngine::from_state(state);
        let res = engine
            .apply(Operation::PcrMutagenesis {
                template: "tpl".to_string(),
                forward_primer: PcrPrimerSpec {
                    sequence: "ATCAAA".to_string(),
                    anneal_len: Some(6),
                    max_mismatches: Some(1),
                    require_3prime_exact_bases: Some(3),
                    library_mode: None,
                    max_variants: None,
                    sample_seed: None,
                },
                reverse_primer: PcrPrimerSpec {
                    sequence: "AAACCC".to_string(),
                    anneal_len: Some(6),
                    max_mismatches: Some(0),
                    require_3prime_exact_bases: Some(4),
                    library_mode: None,
                    max_variants: None,
                    sample_seed: None,
                },
                mutations: vec![SnpMutationSpec {
                    zero_based_position: 2,
                    reference: "G".to_string(),
                    alternate: "C".to_string(),
                }],
                output_id: Some("mut1".to_string()),
                unique: Some(true),
                require_all_mutations: Some(true),
            })
            .unwrap();
        assert_eq!(res.created_seq_ids, vec!["mut1".to_string()]);
        assert_eq!(
            engine
                .state()
                .sequences
                .get("mut1")
                .unwrap()
                .get_forward_string(),
            "ATCAAACCCGGGTTT"
        );
    }

    #[test]
    fn test_pcr_mutagenesis_fails_when_requested_snp_not_introduced() {
        let mut state = ProjectState::default();
        state
            .sequences
            .insert("tpl".to_string(), seq("ATGAAACCCGGGTTT"));
        let mut engine = GentleEngine::from_state(state);
        let err = engine
            .apply(Operation::PcrMutagenesis {
                template: "tpl".to_string(),
                forward_primer: PcrPrimerSpec {
                    sequence: "ATCAAA".to_string(),
                    anneal_len: Some(6),
                    max_mismatches: Some(1),
                    require_3prime_exact_bases: Some(3),
                    library_mode: None,
                    max_variants: None,
                    sample_seed: None,
                },
                reverse_primer: PcrPrimerSpec {
                    sequence: "AAACCC".to_string(),
                    anneal_len: Some(6),
                    max_mismatches: Some(0),
                    require_3prime_exact_bases: Some(4),
                    library_mode: None,
                    max_variants: None,
                    sample_seed: None,
                },
                mutations: vec![SnpMutationSpec {
                    zero_based_position: 2,
                    reference: "G".to_string(),
                    alternate: "T".to_string(),
                }],
                output_id: Some("mut_fail".to_string()),
                unique: Some(true),
                require_all_mutations: Some(true),
            })
            .unwrap_err();
        assert!(err
            .message
            .contains("No amplicon introduced all requested mutations"));
    }

    #[test]
    fn test_pcr_advanced_degenerate_primer_enumerate() {
        let mut state = ProjectState::default();
        state
            .sequences
            .insert("tpl".to_string(), seq("ATGAAACCCGGGTTT"));
        let mut engine = GentleEngine::from_state(state);
        let res = engine
            .apply(Operation::PcrAdvanced {
                template: "tpl".to_string(),
                forward_primer: PcrPrimerSpec {
                    sequence: "ATNAAA".to_string(),
                    anneal_len: Some(6),
                    max_mismatches: Some(1),
                    require_3prime_exact_bases: Some(3),
                    library_mode: Some(PrimerLibraryMode::Enumerate),
                    max_variants: Some(4),
                    sample_seed: None,
                },
                reverse_primer: PcrPrimerSpec {
                    sequence: "AAACCC".to_string(),
                    anneal_len: Some(6),
                    max_mismatches: Some(0),
                    require_3prime_exact_bases: Some(4),
                    library_mode: None,
                    max_variants: None,
                    sample_seed: None,
                },
                output_id: None,
                unique: Some(false),
            })
            .unwrap();
        assert_eq!(res.created_seq_ids.len(), 4);
    }

    #[test]
    fn test_pcr_advanced_degenerate_primer_sample_mode() {
        let mut state = ProjectState::default();
        state
            .sequences
            .insert("tpl".to_string(), seq("ATGAAACCCGGGTTT"));
        let mut engine = GentleEngine::from_state(state);
        let res = engine
            .apply(Operation::PcrAdvanced {
                template: "tpl".to_string(),
                forward_primer: PcrPrimerSpec {
                    sequence: "ATNAAA".to_string(),
                    anneal_len: Some(6),
                    max_mismatches: Some(1),
                    require_3prime_exact_bases: Some(3),
                    library_mode: Some(PrimerLibraryMode::Sample),
                    max_variants: Some(2),
                    sample_seed: Some(42),
                },
                reverse_primer: PcrPrimerSpec {
                    sequence: "AAACCC".to_string(),
                    anneal_len: Some(6),
                    max_mismatches: Some(0),
                    require_3prime_exact_bases: Some(4),
                    library_mode: None,
                    max_variants: None,
                    sample_seed: None,
                },
                output_id: None,
                unique: Some(false),
            })
            .unwrap();
        assert_eq!(res.created_seq_ids.len(), 2);
    }
}
