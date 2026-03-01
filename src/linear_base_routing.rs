//! Adaptive routing decisions for linear DNA base-letter rendering.
//!
//! This module centralizes density-based mode selection so renderer and GUI
//! status panels report identical behavior.

use crate::engine::LinearSequenceLetterLayoutMode;

/// Approximate monospace glyph width multiplier used for density estimates.
pub const GLYPH_WIDTH_SCALE: f32 = 0.62;
/// Auto-routing threshold: dense enough to leave 1-row standard view.
///
/// We intentionally allow a small readability-overlap buffer above 1.0 because
/// the glyph-width estimator is conservative and real rendered glyph advance is
/// often narrower than nominal monospace cell width.
pub const AUTO_STANDARD_MAX_DENSITY: f32 = 1.35;
/// Auto-routing threshold: dense enough to leave 2-row helical view.
pub const AUTO_HELICAL_MAX_DENSITY: f32 = 2.0;
/// Auto-routing threshold: dense enough to leave 10-row condensed view.
pub const AUTO_CONDENSED_MAX_DENSITY: f32 = 10.0;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LinearBaseRenderMode {
    Off,
    StandardLinear,
    ContinuousHelical,
    Condensed10Row,
}

impl LinearBaseRenderMode {
    pub fn label(self) -> &'static str {
        match self {
            Self::Off => "OFF",
            Self::StandardLinear => "STANDARD",
            Self::ContinuousHelical => "HELICAL",
            Self::Condensed10Row => "CONDENSED-10",
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LinearBaseRoutePolicy {
    Auto,
    ForceStandard,
    ForceHelical,
    ForceCondensed10,
}

impl LinearBaseRoutePolicy {
    pub fn label(self) -> &'static str {
        match self {
            Self::Auto => "auto",
            Self::ForceStandard => "force-standard",
            Self::ForceHelical => "force-helical",
            Self::ForceCondensed10 => "force-condensed-10",
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub struct LinearBaseRoutingInput {
    pub span_bp: usize,
    pub viewport_width_px: f32,
    pub compression_enabled: bool,
    pub mode_setting: LinearSequenceLetterLayoutMode,
    pub font_scale: f32,
    pub min_font_size_px: f32,
    pub max_font_size_px: f32,
}

#[derive(Debug, Clone)]
pub struct LinearBaseRoutingDecision {
    pub active_mode: LinearBaseRenderMode,
    pub route_policy: LinearBaseRoutePolicy,
    pub decision_reason: String,
    pub density_ratio: f32,
    pub glyph_width_px: f32,
    pub columns_fit: f32,
    pub estimated_font_size_px: f32,
    pub tier_standard_max_density: f32,
    pub tier_helical_max_density: f32,
    pub tier_condensed_max_density: f32,
}

impl LinearBaseRoutingDecision {
    pub fn bases_visible(&self) -> bool {
        self.active_mode != LinearBaseRenderMode::Off
    }
}

pub fn mode_setting_label(mode: LinearSequenceLetterLayoutMode) -> &'static str {
    match mode {
        LinearSequenceLetterLayoutMode::AutoAdaptive => "auto-adaptive",
        LinearSequenceLetterLayoutMode::StandardLinear => "standard-linear",
        LinearSequenceLetterLayoutMode::ContinuousHelical => "continuous-helical",
        LinearSequenceLetterLayoutMode::Condensed10Row => "condensed-10-row",
    }
}

pub fn estimate_font_size_px(input: LinearBaseRoutingInput) -> f32 {
    let px_per_bp = input.viewport_width_px.max(1.0) / input.span_bp.max(1) as f32;
    (px_per_bp * input.font_scale).clamp(input.min_font_size_px, input.max_font_size_px)
}

pub fn decide_linear_base_routing(input: LinearBaseRoutingInput) -> LinearBaseRoutingDecision {
    let route_policy = match input.mode_setting {
        LinearSequenceLetterLayoutMode::AutoAdaptive => LinearBaseRoutePolicy::Auto,
        LinearSequenceLetterLayoutMode::StandardLinear => LinearBaseRoutePolicy::ForceStandard,
        LinearSequenceLetterLayoutMode::ContinuousHelical => LinearBaseRoutePolicy::ForceHelical,
        LinearSequenceLetterLayoutMode::Condensed10Row => LinearBaseRoutePolicy::ForceCondensed10,
    };
    let estimated_font_size_px = estimate_font_size_px(input);
    let glyph_width_px = (estimated_font_size_px * GLYPH_WIDTH_SCALE).max(1.0);
    let columns_fit = input.viewport_width_px.max(1.0) / glyph_width_px;
    let density_ratio = if columns_fit > 0.0 {
        input.span_bp as f32 / columns_fit
    } else {
        f32::INFINITY
    };

    let (active_mode, decision_reason) = if input.span_bp == 0 {
        (
            LinearBaseRenderMode::Off,
            "inactive: empty viewport".to_string(),
        )
    } else if !input.viewport_width_px.is_finite() || input.viewport_width_px <= 0.0 {
        (
            LinearBaseRenderMode::Off,
            "inactive: invalid viewport width".to_string(),
        )
    } else {
        match route_policy {
            LinearBaseRoutePolicy::Auto => {
                if density_ratio <= AUTO_STANDARD_MAX_DENSITY {
                    (
                        LinearBaseRenderMode::StandardLinear,
                        format!(
                            "auto: density {:.2} <= {:.2} -> standard",
                            density_ratio, AUTO_STANDARD_MAX_DENSITY
                        ),
                    )
                } else if !input.compression_enabled {
                    (
                        LinearBaseRenderMode::Off,
                        format!(
                            "inactive: density {:.2} exceeds standard {:.2} and compressed letters are disabled",
                            density_ratio, AUTO_STANDARD_MAX_DENSITY
                        ),
                    )
                } else if density_ratio <= AUTO_HELICAL_MAX_DENSITY {
                    (
                        LinearBaseRenderMode::ContinuousHelical,
                        format!(
                            "auto: density {:.2} <= {:.2} -> helical",
                            density_ratio, AUTO_HELICAL_MAX_DENSITY
                        ),
                    )
                } else if density_ratio <= AUTO_CONDENSED_MAX_DENSITY {
                    (
                        LinearBaseRenderMode::Condensed10Row,
                        format!(
                            "auto: density {:.2} <= {:.2} -> condensed-10",
                            density_ratio, AUTO_CONDENSED_MAX_DENSITY
                        ),
                    )
                } else {
                    (
                        LinearBaseRenderMode::Off,
                        format!(
                            "inactive: density {:.2} exceeds condensed capacity {:.2}",
                            density_ratio, AUTO_CONDENSED_MAX_DENSITY
                        ),
                    )
                }
            }
            LinearBaseRoutePolicy::ForceStandard => (
                LinearBaseRenderMode::StandardLinear,
                "forced: standard".to_string(),
            ),
            LinearBaseRoutePolicy::ForceHelical => (
                LinearBaseRenderMode::ContinuousHelical,
                "forced: helical".to_string(),
            ),
            LinearBaseRoutePolicy::ForceCondensed10 => (
                LinearBaseRenderMode::Condensed10Row,
                "forced: condensed-10".to_string(),
            ),
        }
    };

    LinearBaseRoutingDecision {
        active_mode,
        route_policy,
        decision_reason,
        density_ratio,
        glyph_width_px,
        columns_fit,
        estimated_font_size_px,
        tier_standard_max_density: AUTO_STANDARD_MAX_DENSITY,
        tier_helical_max_density: AUTO_HELICAL_MAX_DENSITY,
        tier_condensed_max_density: AUTO_CONDENSED_MAX_DENSITY,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn input(span_bp: usize, width_px: f32) -> LinearBaseRoutingInput {
        LinearBaseRoutingInput {
            span_bp,
            viewport_width_px: width_px,
            compression_enabled: true,
            mode_setting: LinearSequenceLetterLayoutMode::AutoAdaptive,
            font_scale: 0.85,
            min_font_size_px: 8.0,
            max_font_size_px: 14.0,
        }
    }

    #[test]
    fn auto_uses_standard_at_or_below_one_x_density() {
        let decision = decide_linear_base_routing(input(100, 1000.0));
        assert_eq!(decision.active_mode, LinearBaseRenderMode::StandardLinear);
        assert!(decision.density_ratio <= AUTO_STANDARD_MAX_DENSITY);
    }

    #[test]
    fn auto_keeps_standard_for_borderline_density_after_tolerance() {
        let decision = decide_linear_base_routing(input(320, 1224.0));
        assert_eq!(decision.active_mode, LinearBaseRenderMode::StandardLinear);
        assert!(decision.density_ratio > 1.0);
        assert!(decision.density_ratio <= AUTO_STANDARD_MAX_DENSITY);
    }

    #[test]
    fn auto_uses_helical_for_mid_density_when_enabled() {
        let decision = decide_linear_base_routing(input(300, 1000.0));
        assert_eq!(
            decision.active_mode,
            LinearBaseRenderMode::ContinuousHelical
        );
        assert!(decision.density_ratio > AUTO_STANDARD_MAX_DENSITY);
        assert!(decision.density_ratio <= AUTO_HELICAL_MAX_DENSITY);
    }

    #[test]
    fn auto_uses_condensed_for_high_density_when_enabled() {
        let decision = decide_linear_base_routing(input(650, 1000.0));
        assert_eq!(decision.active_mode, LinearBaseRenderMode::Condensed10Row);
        assert!(decision.density_ratio > AUTO_HELICAL_MAX_DENSITY);
        assert!(decision.density_ratio <= AUTO_CONDENSED_MAX_DENSITY);
    }

    #[test]
    fn auto_turns_off_when_density_exceeds_condensed_capacity() {
        let decision = decide_linear_base_routing(input(3000, 1000.0));
        assert_eq!(decision.active_mode, LinearBaseRenderMode::Off);
        assert!(decision.density_ratio > AUTO_CONDENSED_MAX_DENSITY);
    }

    #[test]
    fn force_modes_override_auto_tiers() {
        let mut params = input(3000, 1000.0);
        params.mode_setting = LinearSequenceLetterLayoutMode::StandardLinear;
        let standard = decide_linear_base_routing(params);
        assert_eq!(standard.active_mode, LinearBaseRenderMode::StandardLinear);

        params.mode_setting = LinearSequenceLetterLayoutMode::ContinuousHelical;
        let helical = decide_linear_base_routing(params);
        assert_eq!(helical.active_mode, LinearBaseRenderMode::ContinuousHelical);

        params.mode_setting = LinearSequenceLetterLayoutMode::Condensed10Row;
        let condensed = decide_linear_base_routing(params);
        assert_eq!(condensed.active_mode, LinearBaseRenderMode::Condensed10Row);
    }
}
