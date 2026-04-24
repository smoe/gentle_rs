//! Feature-linked splicing/RNA-mapping action launchers for `MainAreaDna`.
//!
//! This submodule keeps the "open Splicing Expert / RNA-read Mapping /
//! transcript-derived dotplot from the current feature" behavior together so
//! the main sequence-window file does not need to interleave action-launch
//! plumbing with unrelated map/layout code.

use super::*;

impl MainAreaDna {
    pub(super) fn feature_kind_supports_splicing_expert(kind_upper: &str) -> bool {
        matches!(
            kind_upper,
            "MRNA" | "TRANSCRIPT" | "NCRNA" | "MISC_RNA" | "EXON" | "GENE" | "CDS"
        )
    }

    pub(super) fn feature_supports_splicing_expert(&self, feature_id: usize) -> bool {
        self.dna
            .read()
            .ok()
            .and_then(|dna| {
                dna.features().get(feature_id).map(|feature| {
                    let kind_upper = feature.kind.to_string().trim().to_ascii_uppercase();
                    Self::feature_kind_supports_splicing_expert(kind_upper.as_str())
                })
            })
            .unwrap_or(false)
    }

    pub(super) fn splicing_expert_view_for_feature(
        &self,
        feature_id: usize,
    ) -> Result<SplicingExpertView, String> {
        if !self.feature_supports_splicing_expert(feature_id) {
            return Err("Selected feature does not support splicing-linked actions".to_string());
        }
        match self.inspect_feature_expert_target(&FeatureExpertTarget::SplicingFeature {
            feature_id,
            scope: SplicingScopePreset::AllOverlappingAnyStrand,
        }) {
            Ok(FeatureExpertView::Splicing(view)) => Ok(view),
            Ok(_) => Err("Selected feature does not expose a splicing expert payload".to_string()),
            Err(err) => Err(err),
        }
    }

    pub(super) fn open_splicing_expert_window_for_view(&mut self, view: &SplicingExpertView) {
        self.log_splicing_expert_status(
            view,
            "open requested",
            self.splicing_expert_window_pending_initial_render,
        );
        self.splicing_expert_window_feature_id = Some(view.target_feature_id);
        self.splicing_expert_window_view = Some(Arc::new(view.clone()));
        self.splicing_expert_window_pending_initial_render = true;
        self.show_splicing_expert_window = true;
        self.log_splicing_expert_status(view, "window state stored", true);
    }

    pub(super) fn focus_splicing_expert_window_view(
        &mut self,
        ctx: &egui::Context,
        view: &SplicingExpertView,
    ) {
        let viewport_id = Self::splicing_expert_viewport_id(&view.seq_id, view.target_feature_id);
        self.log_splicing_expert_status(
            view,
            "focus requested",
            self.splicing_expert_window_pending_initial_render,
        );
        ctx.send_viewport_cmd_to(viewport_id, egui::ViewportCommand::Visible(true));
        ctx.send_viewport_cmd_to(viewport_id, egui::ViewportCommand::Focus);
    }

    pub(super) fn open_splicing_expert_for_feature(
        &mut self,
        feature_id: usize,
        source: &str,
    ) -> bool {
        match self.splicing_expert_view_for_feature(feature_id) {
            Ok(view) => {
                self.open_splicing_expert_window_for_view(&view);
                self.description_cache_expert_view = Some(FeatureExpertView::Splicing(view));
                self.description_cache_expert_error = None;
                self.op_status = format!("Opened Splicing Expert from {source}");
                true
            }
            Err(err) => {
                self.description_cache_expert_error = Some(err.clone());
                self.op_status = format!("Could not open Splicing Expert from {source}: {err}");
                false
            }
        }
    }

    pub(super) fn open_rna_read_mapping_for_feature(
        &mut self,
        feature_id: usize,
        source: &str,
    ) -> bool {
        match self.splicing_expert_view_for_feature(feature_id) {
            Ok(view) => {
                self.open_rna_read_mapping_workspace_for_view(&view);
                self.description_cache_expert_view =
                    Some(FeatureExpertView::Splicing(view.clone()));
                self.description_cache_expert_error = None;
                self.op_status = format!("Opened RNA-read Mapping from {source}");
                true
            }
            Err(err) => {
                self.description_cache_expert_error = Some(err.clone());
                self.op_status = format!("Could not open RNA-read Mapping from {source}: {err}");
                false
            }
        }
    }

    pub(super) fn open_dotplot_for_feature(&mut self, feature_id: usize, source: &str) -> bool {
        match self.splicing_expert_view_for_feature(feature_id) {
            Ok(view) => {
                self.description_cache_expert_view =
                    Some(FeatureExpertView::Splicing(view.clone()));
                self.description_cache_expert_error = None;
                self.derive_selected_transcript_and_open_dotplot(&view);
                if self.show_dotplot_window {
                    self.op_status = format!("Opened transcript dotplot from {source}");
                }
                true
            }
            Err(err) => {
                self.description_cache_expert_error = Some(err.clone());
                self.op_status = format!("Could not open transcript dotplot from {source}: {err}");
                false
            }
        }
    }
}
