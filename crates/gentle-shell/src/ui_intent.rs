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

    /// Parse stable shell spellings and common action aliases.
    pub fn parse(raw: &str) -> Option<Self> {
        match raw.trim().to_ascii_lowercase().as_str() {
            "open" | "show" | "launch" => Some(Self::Open),
            "focus" | "activate" | "raise" => Some(Self::Focus),
            _ => None,
        }
    }
}

/// GUI-facing UI-intent destination understood by adapter shells.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum UiIntentTarget {
    OpenSequence,
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

const UI_INTENT_TARGETS: [UiIntentTarget; UiIntentTarget::COUNT] = [
    UiIntentTarget::OpenSequence,
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
const UI_INTENT_ARGUMENT_GENOME_ID: UiIntentArgument = UiIntentArgument {
    name: "genome_id",
    required: false,
    detail: "Preferred prepared genome/reference identifier to select or focus when the host supports it.",
};
const UI_INTENT_ARGUMENT_HELPERS: UiIntentArgument = UiIntentArgument {
    name: "helpers",
    required: false,
    detail:
        "Query helper-genome catalogs and caches instead of reference-genome catalogs and caches.",
};
const UI_INTENT_ARGUMENT_CATALOG_PATH: UiIntentArgument = UiIntentArgument {
    name: "catalog_path",
    required: false,
    detail: "Override path to the prepared-genome catalog used for prepared-reference selection.",
};
const UI_INTENT_ARGUMENT_CACHE_DIR: UiIntentArgument = UiIntentArgument {
    name: "cache_dir",
    required: false,
    detail: "Override prepared-genome cache directory used for prepared-reference selection.",
};
const UI_INTENT_ARGUMENT_FILTER: UiIntentArgument = UiIntentArgument {
    name: "filter",
    required: false,
    detail: "Case-insensitive text filter applied to prepared-reference catalog rows.",
};
const UI_INTENT_ARGUMENT_SPECIES: UiIntentArgument = UiIntentArgument {
    name: "species",
    required: false,
    detail: "Species token used to narrow prepared-reference selection or latest lookup.",
};
const UI_INTENT_ARGUMENT_LATEST: UiIntentArgument = UiIntentArgument {
    name: "latest",
    required: false,
    detail:
        "Prefer the latest matching prepared reference when deterministic selection is available.",
};
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
const UI_INTENT_ARGUMENTS_DEFAULT: [UiIntentArgument; 1] = [UI_INTENT_ARGUMENT_GENOME_ID];
const UI_INTENT_ARGUMENTS_PREPARED_REFERENCES: [UiIntentArgument; 7] = [
    UI_INTENT_ARGUMENT_GENOME_ID,
    UI_INTENT_ARGUMENT_HELPERS,
    UI_INTENT_ARGUMENT_CATALOG_PATH,
    UI_INTENT_ARGUMENT_CACHE_DIR,
    UI_INTENT_ARGUMENT_FILTER,
    UI_INTENT_ARGUMENT_SPECIES,
    UI_INTENT_ARGUMENT_LATEST,
];

/// Stable per-argument discoverability contract for one GUI intent target.
#[derive(Debug, Clone, Copy, Serialize)]
pub struct UiIntentArgument {
    pub name: &'static str,
    pub required: bool,
    pub detail: &'static str,
}

/// Stable discoverability row for one GUI intent target.
#[derive(Debug, Clone, Serialize)]
pub struct UiIntentTargetCatalogRow {
    pub target: &'static str,
    pub title: &'static str,
    pub detail: &'static str,
    pub keywords: &'static str,
    pub menu_path: &'static str,
    pub actions: &'static [&'static str],
    // deprecated-once-consumers-migrate: use `arguments` for structured
    // name/required/detail metadata while keeping this compatibility surface.
    pub optional_arguments: &'static [&'static str],
    pub arguments: &'static [UiIntentArgument],
}

impl UiIntentTarget {
    /// Number of stable UI-intent destinations.
    pub const COUNT: usize = 12;

    /// Stable catalog order used by shell, MCP, and GUI discoverability.
    pub fn all() -> &'static [Self] {
        &UI_INTENT_TARGETS
    }

    /// Parse stable shell spellings and common aliases.
    pub fn parse(raw: &str) -> Option<Self> {
        match raw.trim().to_ascii_lowercase().as_str() {
            "open-sequence" | "open_sequence" | "sequence-file" | "sequence_file" => {
                Some(Self::OpenSequence)
            }
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
            Self::OpenSequence => "open-sequence",
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
            Self::OpenSequence => "Open Sequence",
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
            Self::OpenSequence => "Open a FASTA, GenBank, EMBL, SnapGene, or XML sequence file.",
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
            Self::OpenSequence => "open sequence import file fasta genbank snapgene embl xml",
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
            Self::OpenSequence | Self::AgentAssistant => "File",
            Self::PcrDesign | Self::SequencingConfirmation => "Patterns",
            Self::PreparedReferences
            | Self::PrepareReferenceGenome
            | Self::RetrieveGenomeSequence
            | Self::BlastGenomeSequence
            | Self::ImportGenomeTrack
            | Self::PrepareHelperGenome
            | Self::RetrieveHelperSequence
            | Self::BlastHelperSequence => "Genome",
        }
    }

    /// Stable action names accepted for this target.
    pub fn actions(self) -> &'static [&'static str] {
        match self {
            Self::OpenSequence
            | Self::PreparedReferences
            | Self::PrepareReferenceGenome
            | Self::RetrieveGenomeSequence
            | Self::BlastGenomeSequence
            | Self::ImportGenomeTrack
            | Self::PcrDesign
            | Self::SequencingConfirmation
            | Self::AgentAssistant
            | Self::PrepareHelperGenome
            | Self::RetrieveHelperSequence
            | Self::BlastHelperSequence => &UI_INTENT_ACTION_NAMES,
        }
    }

    /// Stable optional arguments accepted by the target's discoverability contract.
    pub fn optional_arguments(self) -> &'static [&'static str] {
        match self {
            Self::PreparedReferences => &UI_INTENT_OPTIONAL_ARGUMENTS_PREPARED_REFERENCES,
            _ => &UI_INTENT_OPTIONAL_ARGUMENTS_DEFAULT,
        }
    }

    /// Structured argument contract accepted by this target.
    pub fn arguments(self) -> &'static [UiIntentArgument] {
        match self {
            Self::PreparedReferences => &UI_INTENT_ARGUMENTS_PREPARED_REFERENCES,
            Self::OpenSequence
            | Self::PrepareReferenceGenome
            | Self::RetrieveGenomeSequence
            | Self::BlastGenomeSequence
            | Self::ImportGenomeTrack
            | Self::PcrDesign
            | Self::SequencingConfirmation
            | Self::AgentAssistant
            | Self::PrepareHelperGenome
            | Self::RetrieveHelperSequence
            | Self::BlastHelperSequence => &UI_INTENT_ARGUMENTS_DEFAULT,
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
            actions: self.actions(),
            optional_arguments: self.optional_arguments(),
            arguments: self.arguments(),
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

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::{BTreeMap, BTreeSet};

    const TARGET_ALIAS_CASES: &[(&str, UiIntentTarget)] = &[
        ("open-sequence", UiIntentTarget::OpenSequence),
        ("open_sequence", UiIntentTarget::OpenSequence),
        ("sequence-file", UiIntentTarget::OpenSequence),
        ("sequence_file", UiIntentTarget::OpenSequence),
        ("prepared-references", UiIntentTarget::PreparedReferences),
        ("prepared_references", UiIntentTarget::PreparedReferences),
        ("prepared", UiIntentTarget::PreparedReferences),
        (
            "prepare-reference-genome",
            UiIntentTarget::PrepareReferenceGenome,
        ),
        (
            "prepare_reference_genome",
            UiIntentTarget::PrepareReferenceGenome,
        ),
        ("prepare-genome", UiIntentTarget::PrepareReferenceGenome),
        ("prepare_genome", UiIntentTarget::PrepareReferenceGenome),
        (
            "retrieve-genome-sequence",
            UiIntentTarget::RetrieveGenomeSequence,
        ),
        (
            "retrieve_genome_sequence",
            UiIntentTarget::RetrieveGenomeSequence,
        ),
        ("retrieve-genome", UiIntentTarget::RetrieveGenomeSequence),
        ("retrieve_genome", UiIntentTarget::RetrieveGenomeSequence),
        ("blast-genome-sequence", UiIntentTarget::BlastGenomeSequence),
        ("blast_genome_sequence", UiIntentTarget::BlastGenomeSequence),
        ("blast-genome", UiIntentTarget::BlastGenomeSequence),
        ("blast_genome", UiIntentTarget::BlastGenomeSequence),
        ("import-genome-track", UiIntentTarget::ImportGenomeTrack),
        ("import_genome_track", UiIntentTarget::ImportGenomeTrack),
        ("genome-track", UiIntentTarget::ImportGenomeTrack),
        ("tracks-import", UiIntentTarget::ImportGenomeTrack),
        ("pcr-design", UiIntentTarget::PcrDesign),
        ("pcr_design", UiIntentTarget::PcrDesign),
        ("pcr", UiIntentTarget::PcrDesign),
        ("pcr-designer", UiIntentTarget::PcrDesign),
        ("pcr_designer", UiIntentTarget::PcrDesign),
        (
            "sequencing-confirmation",
            UiIntentTarget::SequencingConfirmation,
        ),
        (
            "sequencing_confirmation",
            UiIntentTarget::SequencingConfirmation,
        ),
        ("seq-confirm", UiIntentTarget::SequencingConfirmation),
        ("seq_confirm", UiIntentTarget::SequencingConfirmation),
        (
            "sequence-confirmation",
            UiIntentTarget::SequencingConfirmation,
        ),
        ("agent-assistant", UiIntentTarget::AgentAssistant),
        ("agent_assistant", UiIntentTarget::AgentAssistant),
        ("agent", UiIntentTarget::AgentAssistant),
        ("prepare-helper-genome", UiIntentTarget::PrepareHelperGenome),
        ("prepare_helper_genome", UiIntentTarget::PrepareHelperGenome),
        ("prepare-helper", UiIntentTarget::PrepareHelperGenome),
        ("prepare_helper", UiIntentTarget::PrepareHelperGenome),
        (
            "retrieve-helper-sequence",
            UiIntentTarget::RetrieveHelperSequence,
        ),
        (
            "retrieve_helper_sequence",
            UiIntentTarget::RetrieveHelperSequence,
        ),
        ("retrieve-helper", UiIntentTarget::RetrieveHelperSequence),
        ("retrieve_helper", UiIntentTarget::RetrieveHelperSequence),
        ("blast-helper-sequence", UiIntentTarget::BlastHelperSequence),
        ("blast_helper_sequence", UiIntentTarget::BlastHelperSequence),
        ("blast-helper", UiIntentTarget::BlastHelperSequence),
        ("blast_helper", UiIntentTarget::BlastHelperSequence),
    ];

    const ACTION_ALIAS_CASES: &[(&str, UiIntentAction)] = &[
        ("open", UiIntentAction::Open),
        ("show", UiIntentAction::Open),
        ("launch", UiIntentAction::Open),
        ("focus", UiIntentAction::Focus),
        ("activate", UiIntentAction::Focus),
        ("raise", UiIntentAction::Focus),
    ];

    fn assert_exhaustive_target_match(target: UiIntentTarget) {
        match target {
            UiIntentTarget::OpenSequence
            | UiIntentTarget::PreparedReferences
            | UiIntentTarget::PrepareReferenceGenome
            | UiIntentTarget::RetrieveGenomeSequence
            | UiIntentTarget::BlastGenomeSequence
            | UiIntentTarget::ImportGenomeTrack
            | UiIntentTarget::PcrDesign
            | UiIntentTarget::SequencingConfirmation
            | UiIntentTarget::AgentAssistant
            | UiIntentTarget::PrepareHelperGenome
            | UiIntentTarget::RetrieveHelperSequence
            | UiIntentTarget::BlastHelperSequence => {}
        }
    }

    #[test]
    fn ui_intent_targets_round_trip_from_all() {
        let explicit_targets = [
            UiIntentTarget::OpenSequence,
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

        assert_eq!(UiIntentTarget::all(), explicit_targets);
        assert_eq!(UiIntentTarget::all().len(), UiIntentTarget::COUNT);
        let mut seen = BTreeSet::new();
        for target in UiIntentTarget::all() {
            assert_exhaustive_target_match(*target);
            assert!(
                seen.insert(target.as_str()),
                "duplicate target spelling {}",
                target.as_str()
            );
            assert_eq!(UiIntentTarget::parse(target.as_str()), Some(*target));
        }
    }

    #[test]
    fn ui_intent_target_aliases_are_unique() {
        let mut aliases = BTreeMap::new();
        for (alias, expected) in TARGET_ALIAS_CASES {
            let parsed = UiIntentTarget::parse(alias)
                .unwrap_or_else(|| panic!("alias {alias} did not parse"));
            assert_eq!(parsed, *expected, "alias {alias} parsed unexpectedly");
            if let Some(previous) = aliases.insert(*alias, parsed) {
                panic!("alias {alias} maps to both {previous:?} and {parsed:?}");
            }
        }
    }

    #[test]
    fn ui_intent_actions_round_trip_and_aliases_are_unique() {
        assert_eq!(
            UiIntentAction::parse(UiIntentAction::Open.as_str()),
            Some(UiIntentAction::Open)
        );
        assert_eq!(
            UiIntentAction::parse(UiIntentAction::Focus.as_str()),
            Some(UiIntentAction::Focus)
        );
        let mut aliases = BTreeMap::new();
        for (alias, expected) in ACTION_ALIAS_CASES {
            let parsed = UiIntentAction::parse(alias)
                .unwrap_or_else(|| panic!("action alias {alias} did not parse"));
            assert_eq!(
                parsed, *expected,
                "action alias {alias} parsed unexpectedly"
            );
            if let Some(previous) = aliases.insert(*alias, parsed) {
                panic!("action alias {alias} maps to both {previous:?} and {parsed:?}");
            }
        }
    }

    #[test]
    fn ui_intent_catalog_rows_are_complete() {
        let known_menus = BTreeSet::from(["File", "Genome", "Patterns"]);
        for target in UiIntentTarget::all() {
            let row = target.catalog_row();
            assert!(!row.title.trim().is_empty());
            assert!(!row.detail.trim().is_empty());
            assert!(!row.keywords.trim().is_empty());
            assert!(
                known_menus.contains(row.menu_path),
                "unexpected menu path {}",
                row.menu_path
            );
            assert!(!row.actions.is_empty());
            for action in row.actions {
                assert!(
                    UiIntentAction::parse(action).is_some(),
                    "unknown action {action}"
                );
            }
            assert_eq!(row.optional_arguments.len(), row.arguments.len());
            for (optional_name, argument) in row.optional_arguments.iter().zip(row.arguments) {
                assert_eq!(*optional_name, argument.name);
                assert!(!argument.name.trim().is_empty());
                assert!(!argument.detail.trim().is_empty());
                assert!(!argument.required);
            }
        }
    }
}
