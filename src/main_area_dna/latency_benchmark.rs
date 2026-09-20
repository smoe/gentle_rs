//! Feature-gated CPU audit controls over the real DNA-viewer presentation path.
//!
//! These are not GUI acceptance actions or a command/agent interface.

use super::*;

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
}
