//! Closed semantic-control vocabulary for tutorial GUI acceptance contracts.
//!
//! The semantic snapshot can expose more controls than tutorials may drive.
//! This catalog is deliberately narrower: each accepted target declares its
//! window, allowed interaction, persistence behavior, and scientific authority.
//! Tutorial contracts therefore cannot turn an arbitrary semantic text field or
//! command-capable control into an execution path.

pub const WINDOW_MAIN: &str = "window.main";
pub const WINDOW_DNA_VIEWER: &str = "window.dna_viewer";
pub const WINDOW_PCR_DESIGN: &str = "window.pcr_design";
pub const WINDOW_GENOMIC_REGIONS: &str = "window.genomic_regions";
pub const WINDOW_SPLICING_EXPERT: &str = "window.splicing_expert";
pub const MAIN_PROJECT_SEQUENCE_OPEN: &str = "main.project.sequence.open";
pub const MAIN_PROJECT_SAVE_STATE: &str = "main.project.save_state";
pub const DNA_SELECTION_FORMULA_INPUT: &str = "dna.selection_formula.input";
pub const DNA_SELECTION_FORMULA_APPLY: &str = "dna.selection_formula.apply";
pub const DNA_SELECTION_STATUS: &str = "dna.selection.status";
pub const DNA_MAP_CANVAS: &str = "dna.map.canvas";
pub const DNA_SELECTION_SIMPLE_PCR: &str = "dna.selection.simple_pcr";
pub const PCR_SIMPLE_STARTER_PANEL: &str = "pcr.simple_starter.panel";
pub const PCR_SIMPLE_STARTER_SEED_FROM_SELECTION: &str = "pcr.simple_starter.seed_from_selection";
pub const PCR_DESIGN_REPORT_ID: &str = "pcr.design.report_id";
pub const PCR_DESIGN_RUN: &str = "pcr.design.run";
pub const DNA_SELECTION_SAVE_REGION: &str = "dna.selection.save_genomic_region";
pub const GENOMIC_REGION_REFRESH: &str = "genomic_region.refresh";
pub const GENOMIC_REGION_SAVE_PENDING: &str = "genomic_region.save_pending";
pub const GENOMIC_REGION_COPY_HUMAN: &str = "genomic_region.copy_human_1based";
pub const GENOMIC_REGION_COPY_BED: &str = "genomic_region.copy_bed_0based";
pub const GENOMIC_REGION_COPY_JSON: &str = "genomic_region.copy_canonical_json";
pub const GENOMIC_REGION_EXPORT_JSON: &str = "genomic_region.export_set_json";
pub const GENOMIC_REGION_EXPORT_BED: &str = "genomic_region.export_set_bed_manifest";
pub const GENOMIC_REGION_IMPORT_JSON: &str = "genomic_region.import_set_json";
pub const GENOMIC_REGION_IMPORT_BED: &str = "genomic_region.import_set_bed_manifest";
pub const GENOMIC_REGION_CAPTURE_CUTRUN: &str = "genomic_region.capture_cutrun_window";
pub const GENOMIC_REGION_CAPTURE_ENSEMBL: &str = "genomic_region.capture_ensembl_feature";

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TutorialGuiInteractionKind {
    Click,
    RightClick,
    DoubleClick,
    ReplaceText,
    SetCheckbox,
    SelectTab,
    PressKey,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TutorialGuiControlAuthority {
    /// Snapshot-only status or window identity; never an interaction target.
    Observe,
    /// Changes only the current GUI/view state.
    ViewState,
    /// Persists project metadata but does not itself establish a scientific result.
    ProjectMetadata,
    /// Persists an engine-owned scientific result that requires typed proof.
    ScientificState,
}

impl TutorialGuiControlAuthority {
    pub fn persists_project_state(self) -> bool {
        matches!(self, Self::ProjectMetadata | Self::ScientificState)
    }

    pub fn requires_scientific_effect(self) -> bool {
        self == Self::ScientificState
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TutorialGuiTextPolicy {
    SelectionFormula,
    Identifier,
}

#[derive(Debug, Clone, Copy)]
pub struct TutorialGuiControlSpec {
    pub semantic_id: &'static str,
    pub window_id: &'static str,
    pub authority: TutorialGuiControlAuthority,
    pub allowed_interactions: &'static [TutorialGuiInteractionKind],
    pub text_policy: Option<TutorialGuiTextPolicy>,
}

const NO_INTERACTIONS: &[TutorialGuiInteractionKind] = &[];
const CLICK: &[TutorialGuiInteractionKind] = &[TutorialGuiInteractionKind::Click];
const RIGHT_CLICK: &[TutorialGuiInteractionKind] = &[TutorialGuiInteractionKind::RightClick];
const DOUBLE_CLICK: &[TutorialGuiInteractionKind] = &[TutorialGuiInteractionKind::DoubleClick];
const REPLACE_TEXT: &[TutorialGuiInteractionKind] = &[TutorialGuiInteractionKind::ReplaceText];

pub const TUTORIAL_GUI_CONTROLS: &[TutorialGuiControlSpec] = &[
    TutorialGuiControlSpec {
        semantic_id: WINDOW_MAIN,
        window_id: WINDOW_MAIN,
        authority: TutorialGuiControlAuthority::Observe,
        allowed_interactions: NO_INTERACTIONS,
        text_policy: None,
    },
    TutorialGuiControlSpec {
        semantic_id: WINDOW_DNA_VIEWER,
        window_id: WINDOW_DNA_VIEWER,
        authority: TutorialGuiControlAuthority::Observe,
        allowed_interactions: NO_INTERACTIONS,
        text_policy: None,
    },
    TutorialGuiControlSpec {
        semantic_id: MAIN_PROJECT_SEQUENCE_OPEN,
        window_id: WINDOW_MAIN,
        authority: TutorialGuiControlAuthority::ViewState,
        allowed_interactions: DOUBLE_CLICK,
        text_policy: None,
    },
    TutorialGuiControlSpec {
        semantic_id: MAIN_PROJECT_SAVE_STATE,
        window_id: WINDOW_MAIN,
        authority: TutorialGuiControlAuthority::Observe,
        allowed_interactions: NO_INTERACTIONS,
        text_policy: None,
    },
    TutorialGuiControlSpec {
        semantic_id: WINDOW_PCR_DESIGN,
        window_id: WINDOW_PCR_DESIGN,
        authority: TutorialGuiControlAuthority::Observe,
        allowed_interactions: NO_INTERACTIONS,
        text_policy: None,
    },
    TutorialGuiControlSpec {
        semantic_id: WINDOW_GENOMIC_REGIONS,
        window_id: WINDOW_GENOMIC_REGIONS,
        authority: TutorialGuiControlAuthority::Observe,
        allowed_interactions: NO_INTERACTIONS,
        text_policy: None,
    },
    TutorialGuiControlSpec {
        semantic_id: WINDOW_SPLICING_EXPERT,
        window_id: WINDOW_SPLICING_EXPERT,
        authority: TutorialGuiControlAuthority::Observe,
        allowed_interactions: NO_INTERACTIONS,
        text_policy: None,
    },
    TutorialGuiControlSpec {
        semantic_id: DNA_SELECTION_FORMULA_INPUT,
        window_id: WINDOW_DNA_VIEWER,
        authority: TutorialGuiControlAuthority::ProjectMetadata,
        allowed_interactions: REPLACE_TEXT,
        text_policy: Some(TutorialGuiTextPolicy::SelectionFormula),
    },
    TutorialGuiControlSpec {
        semantic_id: DNA_SELECTION_FORMULA_APPLY,
        window_id: WINDOW_DNA_VIEWER,
        authority: TutorialGuiControlAuthority::ViewState,
        allowed_interactions: CLICK,
        text_policy: None,
    },
    TutorialGuiControlSpec {
        semantic_id: DNA_SELECTION_STATUS,
        window_id: WINDOW_DNA_VIEWER,
        authority: TutorialGuiControlAuthority::Observe,
        allowed_interactions: NO_INTERACTIONS,
        text_policy: None,
    },
    TutorialGuiControlSpec {
        semantic_id: DNA_MAP_CANVAS,
        window_id: WINDOW_DNA_VIEWER,
        authority: TutorialGuiControlAuthority::ViewState,
        allowed_interactions: RIGHT_CLICK,
        text_policy: None,
    },
    TutorialGuiControlSpec {
        semantic_id: DNA_SELECTION_SIMPLE_PCR,
        window_id: WINDOW_DNA_VIEWER,
        authority: TutorialGuiControlAuthority::ProjectMetadata,
        allowed_interactions: CLICK,
        text_policy: None,
    },
    TutorialGuiControlSpec {
        semantic_id: PCR_SIMPLE_STARTER_PANEL,
        window_id: WINDOW_PCR_DESIGN,
        authority: TutorialGuiControlAuthority::Observe,
        allowed_interactions: NO_INTERACTIONS,
        text_policy: None,
    },
    TutorialGuiControlSpec {
        semantic_id: PCR_SIMPLE_STARTER_SEED_FROM_SELECTION,
        window_id: WINDOW_PCR_DESIGN,
        authority: TutorialGuiControlAuthority::ProjectMetadata,
        allowed_interactions: CLICK,
        text_policy: None,
    },
    TutorialGuiControlSpec {
        semantic_id: PCR_DESIGN_REPORT_ID,
        window_id: WINDOW_PCR_DESIGN,
        authority: TutorialGuiControlAuthority::ViewState,
        allowed_interactions: REPLACE_TEXT,
        text_policy: Some(TutorialGuiTextPolicy::Identifier),
    },
    TutorialGuiControlSpec {
        semantic_id: PCR_DESIGN_RUN,
        window_id: WINDOW_PCR_DESIGN,
        authority: TutorialGuiControlAuthority::ScientificState,
        allowed_interactions: CLICK,
        text_policy: None,
    },
    TutorialGuiControlSpec {
        semantic_id: DNA_SELECTION_SAVE_REGION,
        window_id: WINDOW_DNA_VIEWER,
        authority: TutorialGuiControlAuthority::ViewState,
        allowed_interactions: CLICK,
        text_policy: None,
    },
    TutorialGuiControlSpec {
        semantic_id: GENOMIC_REGION_REFRESH,
        window_id: WINDOW_GENOMIC_REGIONS,
        authority: TutorialGuiControlAuthority::ViewState,
        allowed_interactions: CLICK,
        text_policy: None,
    },
    TutorialGuiControlSpec {
        semantic_id: GENOMIC_REGION_SAVE_PENDING,
        window_id: WINDOW_GENOMIC_REGIONS,
        authority: TutorialGuiControlAuthority::ScientificState,
        allowed_interactions: CLICK,
        text_policy: None,
    },
    TutorialGuiControlSpec {
        semantic_id: GENOMIC_REGION_COPY_HUMAN,
        window_id: WINDOW_GENOMIC_REGIONS,
        authority: TutorialGuiControlAuthority::ViewState,
        allowed_interactions: CLICK,
        text_policy: None,
    },
    TutorialGuiControlSpec {
        semantic_id: GENOMIC_REGION_COPY_BED,
        window_id: WINDOW_GENOMIC_REGIONS,
        authority: TutorialGuiControlAuthority::ViewState,
        allowed_interactions: CLICK,
        text_policy: None,
    },
    TutorialGuiControlSpec {
        semantic_id: GENOMIC_REGION_COPY_JSON,
        window_id: WINDOW_GENOMIC_REGIONS,
        authority: TutorialGuiControlAuthority::ViewState,
        allowed_interactions: CLICK,
        text_policy: None,
    },
    TutorialGuiControlSpec {
        semantic_id: GENOMIC_REGION_EXPORT_JSON,
        window_id: WINDOW_GENOMIC_REGIONS,
        authority: TutorialGuiControlAuthority::ViewState,
        allowed_interactions: CLICK,
        text_policy: None,
    },
    TutorialGuiControlSpec {
        semantic_id: GENOMIC_REGION_EXPORT_BED,
        window_id: WINDOW_GENOMIC_REGIONS,
        authority: TutorialGuiControlAuthority::ViewState,
        allowed_interactions: CLICK,
        text_policy: None,
    },
    TutorialGuiControlSpec {
        semantic_id: GENOMIC_REGION_IMPORT_JSON,
        window_id: WINDOW_GENOMIC_REGIONS,
        authority: TutorialGuiControlAuthority::ScientificState,
        allowed_interactions: CLICK,
        text_policy: None,
    },
    TutorialGuiControlSpec {
        semantic_id: GENOMIC_REGION_IMPORT_BED,
        window_id: WINDOW_GENOMIC_REGIONS,
        authority: TutorialGuiControlAuthority::ScientificState,
        allowed_interactions: CLICK,
        text_policy: None,
    },
    TutorialGuiControlSpec {
        semantic_id: GENOMIC_REGION_CAPTURE_CUTRUN,
        window_id: WINDOW_DNA_VIEWER,
        authority: TutorialGuiControlAuthority::ScientificState,
        allowed_interactions: CLICK,
        text_policy: None,
    },
    TutorialGuiControlSpec {
        semantic_id: GENOMIC_REGION_CAPTURE_ENSEMBL,
        window_id: WINDOW_SPLICING_EXPERT,
        authority: TutorialGuiControlAuthority::ScientificState,
        allowed_interactions: CLICK,
        text_policy: None,
    },
];

pub fn tutorial_gui_control(semantic_id: &str) -> Option<&'static TutorialGuiControlSpec> {
    TUTORIAL_GUI_CONTROLS
        .iter()
        .find(|spec| spec.semantic_id == semantic_id)
}

/// Build the same path-free subject identity used by semantic GUI snapshots.
///
/// External acceptance runners call this through `gentle_examples_docs`
/// instead of reimplementing the domain separator and length framing.
pub fn pseudonymous_subject_scope(parts: &[&str]) -> String {
    let mut identity = b"gentle.gui.subject_scope.v1\0".to_vec();
    for part in parts {
        let length = u64::try_from(part.len()).expect("semantic scope part length fits u64");
        identity.extend_from_slice(&length.to_be_bytes());
        identity.extend_from_slice(part.as_bytes());
    }
    let digest = crate::digest_utils::sha256_hex_bytes(&identity);
    format!("subject-{}", &digest[..32])
}

pub fn validate_replacement_text(policy: TutorialGuiTextPolicy, text: &str) -> Result<(), String> {
    if text.is_empty() {
        return Err("replacement text must not be empty".to_string());
    }
    if text.len() > 512 || text.chars().any(char::is_control) {
        return Err(
            "replacement text must be at most 512 bytes and contain no controls".to_string(),
        );
    }
    match policy {
        TutorialGuiTextPolicy::SelectionFormula => {
            if !text.trim_start().starts_with('=') {
                return Err("selection-formula replacement text must start with '='".to_string());
            }
        }
        TutorialGuiTextPolicy::Identifier => {
            if text.len() > 128
                || !text
                    .bytes()
                    .all(|byte| byte.is_ascii_alphanumeric() || matches!(byte, b'_' | b'-' | b'.'))
            {
                return Err(
                    "identifier replacement text must use at most 128 ASCII letters, digits, '.', '-' or '_'"
                        .to_string(),
                );
            }
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashSet;

    #[test]
    fn tutorial_gui_control_ids_are_unique_and_windows_are_catalogued() {
        let mut ids = HashSet::new();
        for spec in TUTORIAL_GUI_CONTROLS {
            assert!(
                ids.insert(spec.semantic_id),
                "duplicate {}",
                spec.semantic_id
            );
            let window = tutorial_gui_control(spec.window_id).expect("catalogued window");
            assert_eq!(window.authority, TutorialGuiControlAuthority::Observe);
            assert!(window.allowed_interactions.is_empty());
        }
    }

    #[test]
    fn tutorial_text_policies_reject_ambiguous_or_command_like_values() {
        assert!(
            validate_replacement_text(TutorialGuiTextPolicy::SelectionFormula, "=1 .. 20").is_ok()
        );
        assert!(
            validate_replacement_text(TutorialGuiTextPolicy::SelectionFormula, "cargo run")
                .is_err()
        );
        assert!(
            validate_replacement_text(TutorialGuiTextPolicy::Identifier, "primer_report-1").is_ok()
        );
        assert!(
            validate_replacement_text(TutorialGuiTextPolicy::Identifier, "primer report; rm")
                .is_err()
        );
    }

    #[test]
    fn tutorial_subject_scope_is_length_framed_and_path_free() {
        let scope = pseudonymous_subject_scope(&["tp73_locus"]);
        assert_eq!(scope.len(), 40);
        assert!(scope.starts_with("subject-"));
        assert_ne!(
            pseudonymous_subject_scope(&["ab", "c"]),
            pseudonymous_subject_scope(&["a", "bc"])
        );
        assert!(!scope.contains("tp73"));
    }
}
