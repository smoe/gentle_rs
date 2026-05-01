//! Shared GUI intent vocabulary for shell-like adapters.
//!
//! The shell parser and external adapters use this small catalog to refer to
//! GUI destinations without depending on egui or the root application crate.

use serde::{Deserialize, Serialize};

/// GUI-facing UI-intent verb emitted by shell commands that ask an adapter to
/// open or focus a specific tool surface.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum UiIntentAction {
    Open,
    Focus,
}

impl UiIntentAction {
    /// Stable machine-readable spelling used in shell output payloads.
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Open => "open",
            Self::Focus => "focus",
        }
    }
}

/// GUI-facing UI-intent destination understood by adapter shells.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum UiIntentTarget {
    PreparedReferences,
    PrepareReferenceGenome,
    RetrieveGenomeSequence,
    BlastGenomeSequence,
    ImportGenomeTrack,
    PcrDesign,
    SequencingConfirmation,
    AgentAssistant,
    PrepareHelperGenome,
    RetrieveHelperSequence,
    BlastHelperSequence,
}

const UI_INTENT_TARGETS: [UiIntentTarget; 11] = [
    UiIntentTarget::PreparedReferences,
    UiIntentTarget::PrepareReferenceGenome,
    UiIntentTarget::RetrieveGenomeSequence,
    UiIntentTarget::BlastGenomeSequence,
    UiIntentTarget::ImportGenomeTrack,
    UiIntentTarget::PcrDesign,
    UiIntentTarget::SequencingConfirmation,
    UiIntentTarget::AgentAssistant,
    UiIntentTarget::PrepareHelperGenome,
    UiIntentTarget::RetrieveHelperSequence,
    UiIntentTarget::BlastHelperSequence,
];

const UI_INTENT_ACTION_NAMES: [&str; 2] = ["open", "focus"];
const UI_INTENT_OPTIONAL_ARGUMENTS_DEFAULT: [&str; 1] = ["genome_id"];
const UI_INTENT_OPTIONAL_ARGUMENTS_PREPARED_REFERENCES: [&str; 7] = [
    "genome_id",
    "helpers",
    "catalog_path",
    "cache_dir",
    "filter",
    "species",
    "latest",
];

/// Stable discoverability row for one GUI intent target.
#[derive(Debug, Clone, Serialize)]
pub struct UiIntentTargetCatalogRow {
    pub target: &'static str,
    pub title: &'static str,
    pub detail: &'static str,
    pub keywords: &'static str,
    pub menu_path: &'static str,
    pub actions: &'static [&'static str],
    pub optional_arguments: &'static [&'static str],
}

impl UiIntentTarget {
    /// Stable catalog order used by shell, MCP, and GUI discoverability.
    pub fn all() -> &'static [Self] {
        &UI_INTENT_TARGETS
    }

    /// Parse stable shell spellings and common aliases.
    pub fn parse(raw: &str) -> Option<Self> {
        match raw.trim().to_ascii_lowercase().as_str() {
            "prepared-references" | "prepared_references" | "prepared" => {
                Some(Self::PreparedReferences)
            }
            "prepare-reference-genome"
            | "prepare_reference_genome"
            | "prepare-genome"
            | "prepare_genome" => Some(Self::PrepareReferenceGenome),
            "retrieve-genome-sequence"
            | "retrieve_genome_sequence"
            | "retrieve-genome"
            | "retrieve_genome" => Some(Self::RetrieveGenomeSequence),
            "blast-genome-sequence" | "blast_genome_sequence" | "blast-genome" | "blast_genome" => {
                Some(Self::BlastGenomeSequence)
            }
            "import-genome-track" | "import_genome_track" | "genome-track" | "tracks-import" => {
                Some(Self::ImportGenomeTrack)
            }
            "pcr-design" | "pcr_design" | "pcr" | "pcr-designer" | "pcr_designer" => {
                Some(Self::PcrDesign)
            }
            "sequencing-confirmation"
            | "sequencing_confirmation"
            | "seq-confirm"
            | "seq_confirm"
            | "sequence-confirmation" => Some(Self::SequencingConfirmation),
            "agent-assistant" | "agent_assistant" | "agent" => Some(Self::AgentAssistant),
            "prepare-helper-genome"
            | "prepare_helper_genome"
            | "prepare-helper"
            | "prepare_helper" => Some(Self::PrepareHelperGenome),
            "retrieve-helper-sequence"
            | "retrieve_helper_sequence"
            | "retrieve-helper"
            | "retrieve_helper" => Some(Self::RetrieveHelperSequence),
            "blast-helper-sequence" | "blast_helper_sequence" | "blast-helper" | "blast_helper" => {
                Some(Self::BlastHelperSequence)
            }
            _ => None,
        }
    }

    /// Stable machine-readable spelling used in shell output payloads.
    pub fn as_str(self) -> &'static str {
        match self {
            Self::PreparedReferences => "prepared-references",
            Self::PrepareReferenceGenome => "prepare-reference-genome",
            Self::RetrieveGenomeSequence => "retrieve-genome-sequence",
            Self::BlastGenomeSequence => "blast-genome-sequence",
            Self::ImportGenomeTrack => "import-genome-track",
            Self::PcrDesign => "pcr-design",
            Self::SequencingConfirmation => "sequencing-confirmation",
            Self::AgentAssistant => "agent-assistant",
            Self::PrepareHelperGenome => "prepare-helper-genome",
            Self::RetrieveHelperSequence => "retrieve-helper-sequence",
            Self::BlastHelperSequence => "blast-helper-sequence",
        }
    }

    /// Stable human-facing title reused by UI-intent discovery surfaces.
    pub fn discoverability_title(self) -> &'static str {
        match self {
            Self::PreparedReferences => "Prepared References",
            Self::PrepareReferenceGenome => "Prepare Reference Genome",
            Self::RetrieveGenomeSequence => "Retrieve Genomic Sequence",
            Self::BlastGenomeSequence => "BLAST Genome Sequence",
            Self::ImportGenomeTrack => "Import Genome Track",
            Self::PcrDesign => "PCR Designer",
            Self::SequencingConfirmation => "Sequencing Confirmation",
            Self::AgentAssistant => "Agent Assistant",
            Self::PrepareHelperGenome => "Prepare Helper Genome",
            Self::RetrieveHelperSequence => "Retrieve Helper Sequence",
            Self::BlastHelperSequence => "BLAST Helper Sequence",
        }
    }

    /// Compatibility alias for command-catalog call sites.
    pub fn title(self) -> &'static str {
        self.discoverability_title()
    }

    /// Short discoverability copy for command palettes, agent helpers, and MCP.
    pub fn discoverability_detail(self) -> &'static str {
        match self {
            Self::PreparedReferences => "Inspect prepared reference/helper genome installations.",
            Self::PrepareReferenceGenome => "Download/index the selected reference genome.",
            Self::RetrieveGenomeSequence => {
                "Extract a region or gene from a prepared reference genome."
            }
            Self::BlastGenomeSequence => "Run BLAST against prepared reference genome indices.",
            Self::ImportGenomeTrack => "Import BED/BigWig/VCF tracks onto anchored sequences.",
            Self::PcrDesign => "Paint-first pair-PCR specialist with queue and live geometry.",
            Self::SequencingConfirmation => {
                "Confirm an expected construct from reads or imported traces."
            }
            Self::AgentAssistant => {
                "Ask configured agent systems and run suggested shared-shell commands."
            }
            Self::PrepareHelperGenome => "Download/index the selected helper genome.",
            Self::RetrieveHelperSequence => "Extract sequence from a prepared helper genome.",
            Self::BlastHelperSequence => "Run BLAST against prepared helper genome indices.",
        }
    }

    /// Compatibility alias for command-catalog call sites.
    pub fn detail(self) -> &'static str {
        self.discoverability_detail()
    }

    /// Search keywords shared across UI-intent discoverability surfaces.
    pub fn discoverability_keywords(self) -> &'static str {
        match self {
            Self::PreparedReferences => "genome prepared references helpers inspector",
            Self::PrepareReferenceGenome => "genome prepare reference",
            Self::RetrieveGenomeSequence => "genome retrieve extract anchor region gene",
            Self::BlastGenomeSequence => "genome blast reference",
            Self::ImportGenomeTrack => "genome track tracks bed bigwig vcf import anchored",
            Self::PcrDesign => "pcr primer pair roi paint queue designer",
            Self::SequencingConfirmation => {
                "sequencing confirmation sanger construct reads trace junction"
            }
            Self::AgentAssistant => "agent assistant chat automation support",
            Self::PrepareHelperGenome => "helper genome prepare",
            Self::RetrieveHelperSequence => "helper genome retrieve extract sequence",
            Self::BlastHelperSequence => "helper genome blast",
        }
    }

    /// Compatibility alias for command-catalog call sites.
    pub fn keywords(self) -> &'static str {
        self.discoverability_keywords()
    }

    /// Primary GUI menu location for the target.
    pub fn menu_path(self) -> &'static str {
        match self {
            Self::PcrDesign | Self::SequencingConfirmation => "Patterns",
            Self::AgentAssistant => "File",
            _ => "Genome",
        }
    }

    /// Stable optional arguments accepted by the target's discoverability contract.
    pub fn optional_arguments(self) -> &'static [&'static str] {
        match self {
            Self::PreparedReferences => &UI_INTENT_OPTIONAL_ARGUMENTS_PREPARED_REFERENCES,
            _ => &UI_INTENT_OPTIONAL_ARGUMENTS_DEFAULT,
        }
    }

    /// Structured discoverability row consumed by shell/MCP/GUI surfaces.
    pub fn catalog_row(self) -> UiIntentTargetCatalogRow {
        UiIntentTargetCatalogRow {
            target: self.as_str(),
            title: self.discoverability_title(),
            detail: self.discoverability_detail(),
            keywords: self.discoverability_keywords(),
            menu_path: self.menu_path(),
            actions: &UI_INTENT_ACTION_NAMES,
            optional_arguments: self.optional_arguments(),
        }
    }
}

/// Shared UI-intent catalog used by shell, MCP, and GUI discoverability layers.
pub fn ui_intent_target_catalog() -> Vec<UiIntentTargetCatalogRow> {
    UiIntentTarget::all()
        .iter()
        .copied()
        .map(UiIntentTarget::catalog_row)
        .collect()
}
