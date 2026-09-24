//! Feature-gated CPU audit controls over the real DNA-viewer presentation path.
//!
//! These are not GUI acceptance actions or a command/agent interface.

use super::*;

/// Audit workload identity; the original nine fixture contents remain unchanged.
pub const DNA_LATENCY_WORKLOAD: &str = "density_boundary_v1";
/// Exact proposed interactive boundary, not a performance acceptance verdict.
pub const INTERACTIVE_BOUNDARY: (usize, usize) = (250_000, 5_000);

/// Independent length/count ladder followed by the exact interactive boundary.
pub fn feature_density_workloads() -> impl Iterator<Item = (usize, usize)> {
    [20_000, 250_000, 2_000_000]
        .into_iter()
        .flat_map(|length| [100, 1_000, 10_000].map(|count| (length, count)))
        .chain(std::iter::once(INTERACTIVE_BOUNDARY))
}

/// Deterministic, non-biological density ladder. No external input or RNG state.
/// Half the features cluster in the first 5 kbp; half span the complete locus.
pub fn feature_density_fixture(length: usize, count: usize) -> DNAsequence {
    use gb_io::seq::{Feature, Location};
    assert!(length >= 20_000 && count >= 1);
    let sequence = "ACGTGCAATTCG".repeat(length.div_ceil(12));
    let mut dna = DNAsequence::from_sequence(&sequence[..length]).unwrap();
    let features = dna.features_mut();
    for index in 0..count {
        let extent = if index % 2 == 0 { 5_000 } else { length };
        let start = (index * 137) % (extent - 600);
        let range = |from, to| Location::simple_range(from as i64, to as i64);
        let (kind, location) = match index % 5 {
            0 => (
                "mRNA",
                Location::Join(vec![
                    range(start, start + 90),
                    range(start + 300, start + 510),
                ]),
            ),
            1 => (
                "CDS",
                Location::Join(vec![
                    range(start + 15, start + 90),
                    range(start + 300, start + 480),
                ]),
            ),
            2 => ("exon", range(start, start + 90)),
            3 => ("regulatory", range(start, start + 24)),
            _ => ("repeat_region", range(start, start + 160)),
        };
        features.push(Feature {
            kind: kind.into(),
            location: if index % 3 == 0 {
                Location::Complement(Box::new(location))
            } else {
                location
            },
            qualifiers: vec![
                ("label".into(), Some(format!("synthetic_{index:05}"))),
                ("gene".into(), Some(format!("toy_group_{}", index / 5))),
            ],
        });
    }
    dna
}

/// Artificial boundary case with sparse cut sites, short ORFs and methylation
/// motifs. No biological source, catalog lookup or production project is used.
pub fn feature_density_boundary_fixture() -> DNAsequence {
    let (length, count) = INTERACTIVE_BOUNDARY;
    let original = feature_density_fixture(length, count);
    let mut sequence = original.forward_bytes().to_vec();
    let insert = format!("GAATTCGATCCCAGGATG{}TAA", "GCC".repeat(120));
    for start in (100..length - insert.len()).step_by(9_973) {
        sequence[start..start + insert.len()].copy_from_slice(insert.as_bytes());
    }
    let mut dna = DNAsequence::from_sequence(std::str::from_utf8(&sequence).unwrap()).unwrap();
    dna.features_mut().clone_from(original.features());
    dna.restriction_enzymes_mut().push(
        serde_json::from_value(serde_json::json!({
            "name": "synthetic_boundary_site", "sequence": "GAATTC", "cut": 1, "overlap": 4
        }))
        .unwrap(),
    );
    dna.set_methylation_mode(crate::methylation_sites::MethylationMode::both());
    dna
}

/// Counts used by the viewer's toolbar, not proof of pixel-level visibility.
#[derive(Debug, PartialEq, Eq, serde::Serialize)]
pub struct DnaLatencyLayerCount {
    pub total: usize,
    pub viewport_eligible: usize,
    pub enabled: bool,
}

/// Derived workload accompanying untimed observations. Missing/stale toolbar
/// caches produce no inventory instead of an invented zero count.
#[derive(Debug, PartialEq, Eq, serde::Serialize)]
pub struct DnaLatencyLayerInventory {
    pub annotation_features: usize,
    pub viewport: Option<(usize, usize)>,
    pub gc_bin_size_bp: usize,
    pub restriction_groups: DnaLatencyLayerCount,
    pub gc_bins: DnaLatencyLayerCount,
    pub orfs: DnaLatencyLayerCount,
    pub methylation_sites: DnaLatencyLayerCount,
}

/// Fixed view-only interactions, kept separate from fixture construction.
#[derive(Clone, Copy, Debug)]
pub enum DnaLatencyInteraction {
    Steady,
    Pan,
    Zoom,
    ToggleMrna,
    Select,
}

impl MainAreaDna {
    /// Enable derived layers explicitly for the new boundary workload only.
    pub fn enable_latency_benchmark_derived_layers(&mut self) {
        let mut display = self.dna_display.write().unwrap();
        display.set_restriction_enzyme_display_mode(RestrictionEnzymeDisplayMode::AllInView);
        display.set_show_restriction_enzyme_sites(true);
        display.set_show_gc_contents(true);
        display.set_show_open_reading_frames(true);
        display.set_show_methylation_sites(true);
    }

    /// Inspect the already-built toolbar cache without changing its counters,
    /// populating it, or modifying sequence/display state. Call outside timing.
    pub fn latency_benchmark_layer_inventory(&self) -> Option<DnaLatencyLayerInventory> {
        let viewport = self.active_linear_viewport_range();
        let dna = self.dna.read().ok()?;
        let display = self.dna_display.read().ok()?;
        let key = LayerVisibilityCacheKey {
            display_revision: display.revision(),
            feature_generation: dna.feature_generation(),
            restriction_generation: dna.restriction_enzyme_group_generation(),
            sequence_length: dna.len(),
            viewport,
        };
        let counts = &self
            .layer_visibility_cache
            .as_ref()
            .filter(|c| c.key == key)?
            .counts;
        Some(DnaLatencyLayerInventory {
            annotation_features: dna.features().len(),
            viewport,
            gc_bin_size_bp: display.gc_content_bin_size_bp(),
            restriction_groups: DnaLatencyLayerCount {
                total: dna.restriction_enzyme_groups().len(),
                viewport_eligible: counts.restriction_site_count,
                enabled: display.show_restriction_enzyme_sites(),
            },
            gc_bins: DnaLatencyLayerCount {
                total: dna.len().div_ceil(display.gc_content_bin_size_bp().max(1)),
                viewport_eligible: counts.gc_region_count,
                enabled: display.show_gc_contents(),
            },
            orfs: DnaLatencyLayerCount {
                total: dna.open_reading_frames().len(),
                viewport_eligible: counts.orf_count,
                enabled: display.show_open_reading_frames(),
            },
            methylation_sites: DnaLatencyLayerCount {
                total: dna.methylation_sites().sites().len(),
                viewport_eligible: counts.methylation_site_count,
                enabled: display.show_methylation_sites(),
            },
        })
    }

    /// Set a reproducible partial viewport and explicitly load the feature tree.
    /// No scientific operations or project writes are performed by this setup.
    pub fn prepare_latency_benchmark(&mut self, tree_loaded: bool) {
        self.feature_tree_deferred_until_interaction = !tree_loaded;
        self.dna_display
            .write()
            .unwrap()
            .set_linear_viewport(0, 5_000);
    }

    /// Probe the actual rendered map rectangle, not an assumed screen position.
    pub fn latency_benchmark_hover_position(&self) -> Option<egui::Pos2> {
        match &self.map_dna {
            RenderDna::Linear(renderer) => renderer.read().ok().map(|r| r.area().center()),
            RenderDna::Circular(_) => None,
        }
    }

    /// Exercise the same display setters used by normal viewer interactions.
    pub fn apply_latency_benchmark_interaction(&mut self, interaction: DnaLatencyInteraction) {
        let (start, span, len) = self.current_linear_viewport();
        let mut display = self.dna_display.write().unwrap();
        match interaction {
            DnaLatencyInteraction::Steady => {}
            DnaLatencyInteraction::Pan => {
                display.set_linear_viewport((start + 1).min(len.saturating_sub(span)), span);
            }
            DnaLatencyInteraction::Zoom => {
                display.set_linear_viewport(start, (span / 2).max(1));
            }
            DnaLatencyInteraction::ToggleMrna => {
                let show = !display.show_mrna_features();
                display.set_show_mrna_features(show);
            }
            DnaLatencyInteraction::Select => {
                drop(display);
                self.map_dna.select_feature(Some(0));
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn density_fixture_is_reproducible_with_exact_geometry_and_counts() {
        let first = feature_density_fixture(20_000, 100);
        let second = feature_density_fixture(20_000, 100);
        assert_eq!(first.len(), 20_000);
        assert_eq!(first.features().len(), 100);
        assert_eq!(
            serde_json::to_vec(first.features()).unwrap(),
            serde_json::to_vec(second.features()).unwrap()
        );
        assert_eq!(first.forward_bytes(), second.forward_bytes());
        for feature in first.features() {
            let (start, end) = feature.location.find_bounds().unwrap();
            assert!(0 <= start && start < end && end <= 20_000);
        }
    }

    #[test]
    fn boundary_workload_is_exact_reproducible_and_additive() {
        let workloads = feature_density_workloads().collect::<Vec<_>>();
        assert_eq!(workloads.len(), 10);
        assert_eq!(workloads.iter().copied().collect::<BTreeSet<_>>().len(), 10);
        assert_eq!(workloads.last(), Some(&INTERACTIVE_BOUNDARY));
        let dna = feature_density_boundary_fixture();
        assert_eq!((dna.len(), dna.features().len()), INTERACTIVE_BOUNDARY);
        assert_eq!(
            serde_json::to_value(&dna).unwrap(),
            serde_json::to_value(feature_density_boundary_fixture()).unwrap()
        );
        assert_eq!(
            dna.features(),
            feature_density_fixture(INTERACTIVE_BOUNDARY.0, INTERACTIVE_BOUNDARY.1).features()
        );
        assert_eq!(dna.restriction_enzymes().len(), 1);
        assert!(dna.methylation_mode().dam() && dna.methylation_mode().dcm());
        for feature in dna.features() {
            let (start, end) = feature.location.find_bounds().unwrap();
            assert!(0 <= start && start < end && end <= INTERACTIVE_BOUNDARY.0 as i64);
        }
    }

    #[test]
    fn boundary_layer_inventory_is_current_nonempty_and_observation_only() {
        let mut area = MainAreaDna::new(feature_density_boundary_fixture(), None, None);
        area.prepare_latency_benchmark(true);
        area.enable_latency_benchmark_derived_layers();
        assert!(area.latency_benchmark_layer_inventory().is_none());
        area.compute_layer_visibility_counts();
        let before = area.cache_diagnostics();
        let dna_before = serde_json::to_value(&*area.dna.read().unwrap()).unwrap();
        let inventory = area.latency_benchmark_layer_inventory().unwrap();
        assert_eq!(inventory.annotation_features, 5_000);
        assert_eq!(inventory.viewport, Some((0, 5_000)));
        for layer in [
            &inventory.restriction_groups,
            &inventory.gc_bins,
            &inventory.orfs,
            &inventory.methylation_sites,
        ] {
            assert!(layer.enabled);
            assert!(layer.viewport_eligible > 0);
            assert!(layer.total >= layer.viewport_eligible);
        }
        assert_eq!(inventory, area.latency_benchmark_layer_inventory().unwrap());
        assert_eq!(before, area.cache_diagnostics());
        assert_eq!(
            dna_before,
            serde_json::to_value(&*area.dna.read().unwrap()).unwrap()
        );
        area.apply_latency_benchmark_interaction(DnaLatencyInteraction::Pan);
        assert!(area.latency_benchmark_layer_inventory().is_none());
        area.compute_layer_visibility_counts();
        assert_eq!(
            area.latency_benchmark_layer_inventory().unwrap().viewport,
            Some((1, 5_001))
        );
    }
}
