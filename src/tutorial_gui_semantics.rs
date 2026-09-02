//! Closed semantic-control vocabulary for tutorial GUI acceptance contracts.
//!
//! The semantic snapshot can expose more controls than tutorials may drive.
//! This catalog is deliberately narrower: each accepted target declares its
//! window, allowed interaction, persistence behavior, and scientific authority.
//! Tutorial contracts therefore cannot turn an arbitrary semantic text field or
//! command-capable control into an execution path.

pub const WINDOW_DNA_VIEWER: &str = "window.dna_viewer";
pub const WINDOW_PCR_DESIGN: &str = "window.pcr_design";
pub const DNA_SELECTION_FORMULA_INPUT: &str = "dna.selection_formula.input";
pub const DNA_SELECTION_FORMULA_APPLY: &str = "dna.selection_formula.apply";
pub const DNA_SELECTION_STATUS: &str = "dna.selection.status";
pub const DNA_MAP_CANVAS: &str = "dna.map.canvas";
pub const DNA_SELECTION_SIMPLE_PCR: &str = "dna.selection.simple_pcr";
pub const PCR_SIMPLE_STARTER_PANEL: &str = "pcr.simple_starter.panel";
pub const PCR_SIMPLE_STARTER_SEED_FROM_SELECTION: &str = "pcr.simple_starter.seed_from_selection";
pub const PCR_DESIGN_REPORT_ID: &str = "pcr.design.report_id";
pub const PCR_DESIGN_RUN: &str = "pcr.design.run";

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
const REPLACE_TEXT: &[TutorialGuiInteractionKind] = &[TutorialGuiInteractionKind::ReplaceText];

pub const TUTORIAL_GUI_CONTROLS: &[TutorialGuiControlSpec] = &[
    TutorialGuiControlSpec {
        semantic_id: WINDOW_DNA_VIEWER,
        window_id: WINDOW_DNA_VIEWER,
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
];

pub fn tutorial_gui_control(semantic_id: &str) -> Option<&'static TutorialGuiControlSpec> {
    TUTORIAL_GUI_CONTROLS
        .iter()
        .find(|spec| spec.semantic_id == semantic_id)
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
}
