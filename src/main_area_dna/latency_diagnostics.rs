//! Observation-only DNA-viewer cache diagnostics; never persisted in the project.

use super::*;

impl MainAreaDna {
    /// Read cumulative work counters without rebuilding models or waiting on a renderer.
    pub fn cache_diagnostics(&self) -> crate::gui_profiler::DnaCacheDiagnostics {
        crate::gui_profiler::DnaCacheDiagnostics {
            tree_hits: self.feature_tree_cache_hits,
            tree_builds: self.feature_tree_cache_misses,
            layer_hits: self.layer_visibility_cache_hits,
            layer_builds: self.layer_visibility_cache_misses,
            layer_feature_visits: self.layer_visibility_feature_visits,
            layer_gc_bases: self.layer_visibility_gc_bases,
            display_sync_hits: self.engine_display_sync_cache_hits,
            display_sync_builds: self.engine_display_sync_cache_misses,
            linear: match &self.map_dna {
                RenderDna::Linear(renderer) => {
                    renderer.try_read().ok().map(|r| r.cache_diagnostics())
                }
                RenderDna::Circular(_) => None,
            },
        }
    }

    pub(super) fn render_cache_diagnostics(&self, ui: &mut egui::Ui) {
        let counters = self.cache_diagnostics();
        egui::CollapsingHeader::new("DNA cache diagnostics")
            .id_salt(("dna_cache_diagnostics", self.panel_scope_key()))
            .show(ui, |ui| {
                if let Ok(json) = serde_json::to_string_pretty(&counters) {
                    ui.monospace(json);
                }
            });
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cache_diagnostics_observe_pan_rebuilds_without_mutating_state() {
        let mut area = MainAreaDna::new(
            DNAsequence::from_sequence(&"ACGT".repeat(300)).unwrap(),
            None,
            None,
        );
        area.set_linear_viewport(0, 100);
        area.compute_layer_visibility_counts();
        area.ensure_feature_tree_cache_current(area.active_linear_viewport_range());
        let first = area.cache_diagnostics();
        assert_eq!(first.layer_gc_bases, 1200);
        assert_eq!(first.tree_builds, 1);
        let before = serde_json::to_value(&*area.dna.read().unwrap()).unwrap();
        assert_eq!(first, area.cache_diagnostics());
        area.compute_layer_visibility_counts();
        area.ensure_feature_tree_cache_current(area.active_linear_viewport_range());
        assert_eq!(area.cache_diagnostics().layer_builds, first.layer_builds);
        assert_eq!(area.cache_diagnostics().tree_builds, first.tree_builds);
        area.set_linear_viewport(1, 100);
        area.compute_layer_visibility_counts();
        area.ensure_feature_tree_cache_current(area.active_linear_viewport_range());
        let panned = area.cache_diagnostics();
        assert_eq!(panned.layer_gc_bases, 2400);
        assert_eq!(panned.tree_builds, 2);
        assert_eq!(
            before,
            serde_json::to_value(&*area.dna.read().unwrap()).unwrap()
        );
        if let RenderDna::Linear(renderer) = &area.map_dna {
            let _guard = renderer.write().unwrap();
            assert!(area.cache_diagnostics().linear.is_none());
        }
    }
}
