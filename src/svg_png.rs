//! Deterministic SVG-to-PNG rasterization helpers shared by headless tools.
//!
//! This module keeps the `resvg`-based conversion path in one reusable place so
//! wrapper/CLI/docs flows do not fork rendering behavior when they need a
//! messenger-friendly bitmap artifact from an engine-owned SVG export.

use regex::Regex;
use resvg::{self, tiny_skia, usvg};
use serde::Serialize;
use std::{
    path::{Path, PathBuf},
    sync::LazyLock,
};

static DOTPLOT_METADATA_RE: LazyLock<Regex> = LazyLock::new(|| {
    Regex::new(
        r#"<text\b[^>]*>(?:Dotplot workspace export:[^<]*|rendered_cells=[^<]*|x:\s[^<]*|y:\s[^<]*|GENtle dotplot SVG export)</text>"#,
    )
    .expect("dotplot metadata regex should compile")
});

/// Default fixed raster scale for messenger/chat-facing ClawBio figures.
pub const DEFAULT_CLAWBIO_PNG_SCALE: f32 = 2.0;

/// Options controlling deterministic SVG-to-PNG conversion.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SvgPngRenderOptions {
    /// Output scale relative to the SVG viewbox/intrinsic pixel size.
    pub scale: f32,
    /// Whether known dotplot metadata footer text should be stripped first.
    pub drop_dotplot_metadata: bool,
}

impl Default for SvgPngRenderOptions {
    fn default() -> Self {
        Self {
            scale: 1.0,
            drop_dotplot_metadata: false,
        }
    }
}

/// Machine-readable summary of one deterministic SVG-to-PNG conversion.
#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
pub struct SvgPngRenderSummary {
    /// Source SVG path that was rasterized.
    pub input_path: String,
    /// Output PNG path that was written.
    pub output_path: String,
    /// Applied scale factor.
    pub scale: String,
    /// Whether dotplot metadata stripping was enabled.
    pub drop_dotplot_metadata: bool,
    /// Written PNG width in pixels.
    pub width: u32,
    /// Written PNG height in pixels.
    pub height: u32,
}

/// Removes known dotplot metadata footer/header text while preserving axis
/// labels and figure content.
pub fn strip_dotplot_metadata_text(svg: &str) -> String {
    DOTPLOT_METADATA_RE.replace_all(svg, "").into_owned()
}

fn canonical_parent(path: &Path) -> Option<PathBuf> {
    std::fs::canonicalize(path)
        .ok()
        .and_then(|resolved| resolved.parent().map(|parent| parent.to_path_buf()))
}

/// Rasterizes one SVG file into one PNG file with deterministic `resvg`
/// behavior and an optional metadata-cleanup pass.
pub fn render_svg_file_to_png(
    input_path: &Path,
    output_path: &Path,
    options: SvgPngRenderOptions,
) -> Result<SvgPngRenderSummary, String> {
    if input_path.as_os_str().is_empty() {
        return Err("svg-png requires INPUT.svg".to_string());
    }
    if output_path.as_os_str().is_empty() {
        return Err("svg-png requires OUTPUT.png".to_string());
    }
    if !(options.scale.is_finite() && options.scale > 0.0) {
        return Err(format!(
            "svg-png requires a positive finite scale value, got {}",
            options.scale
        ));
    }

    let mut svg_text = std::fs::read_to_string(input_path)
        .map_err(|e| format!("Could not read SVG '{}': {e}", input_path.display()))?;
    if options.drop_dotplot_metadata {
        svg_text = strip_dotplot_metadata_text(&svg_text);
    }

    let mut opt = usvg::Options {
        resources_dir: canonical_parent(input_path),
        ..usvg::Options::default()
    };
    opt.fontdb_mut().load_system_fonts();

    let tree = usvg::Tree::from_str(&svg_text, &opt)
        .map_err(|e| format!("Could not parse SVG '{}': {e}", input_path.display()))?;
    let pixmap_size = tree
        .size()
        .to_int_size()
        .scale_by(options.scale)
        .ok_or_else(|| format!("Could not scale SVG size by {}", options.scale))?;
    let mut pixmap =
        tiny_skia::Pixmap::new(pixmap_size.width(), pixmap_size.height()).ok_or_else(|| {
            format!(
                "Could not allocate PNG canvas {}x{}",
                pixmap_size.width(),
                pixmap_size.height()
            )
        })?;

    let transform = if (options.scale - 1.0).abs() <= f32::EPSILON {
        tiny_skia::Transform::default()
    } else {
        tiny_skia::Transform::from_scale(options.scale, options.scale)
    };
    resvg::render(&tree, transform, &mut pixmap.as_mut());
    pixmap
        .save_png(output_path)
        .map_err(|e| format!("Could not write PNG '{}': {e}", output_path.display()))?;

    Ok(SvgPngRenderSummary {
        input_path: input_path.to_string_lossy().into_owned(),
        output_path: output_path.to_string_lossy().into_owned(),
        scale: format!("{}", options.scale),
        drop_dotplot_metadata: options.drop_dotplot_metadata,
        width: pixmap_size.width(),
        height: pixmap_size.height(),
    })
}

#[cfg(test)]
mod tests {
    use super::{
        SvgPngRenderOptions, render_svg_file_to_png, strip_dotplot_metadata_text,
    };
    use image::GenericImageView;

    #[test]
    fn strip_dotplot_metadata_preserves_axis_labels() {
        let svg = concat!(
            "<svg>",
            "<text>Dotplot workspace export: dotplot_primary</text>",
            "<text>rendered_cells=42 sampled_points=77 sample_stride=1</text>",
            "<text>x: tp73_cdna</text>",
            "<text>y: tp73_genomic</text>",
            "<text>1</text>",
            "<text>5026</text>",
            "<text>GENtle dotplot SVG export</text>",
            "</svg>"
        );
        let cleaned = strip_dotplot_metadata_text(svg);
        assert!(!cleaned.contains("Dotplot workspace export:"));
        assert!(!cleaned.contains("rendered_cells="));
        assert!(!cleaned.contains("x: tp73_cdna"));
        assert!(!cleaned.contains("y: tp73_genomic"));
        assert!(!cleaned.contains("GENtle dotplot SVG export"));
        assert!(cleaned.contains("<text>1</text>"));
        assert!(cleaned.contains("<text>5026</text>"));
    }

    #[test]
    fn render_svg_file_to_png_scales_deterministically() {
        let tmp = tempfile::tempdir().expect("tempdir");
        let input = tmp.path().join("demo.svg");
        let output = tmp.path().join("demo.png");
        std::fs::write(
            &input,
            concat!(
                "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"40\" height=\"20\" viewBox=\"0 0 40 20\">",
                "<rect x=\"0\" y=\"0\" width=\"40\" height=\"20\" fill=\"#ffffff\"/>",
                "<rect x=\"5\" y=\"5\" width=\"30\" height=\"10\" fill=\"#2563eb\"/>",
                "</svg>"
            ),
        )
        .expect("write svg");

        let summary = render_svg_file_to_png(
            &input,
            &output,
            SvgPngRenderOptions {
                scale: 2.0,
                drop_dotplot_metadata: false,
            },
        )
        .expect("render png");

        assert_eq!(summary.width, 80);
        assert_eq!(summary.height, 40);

        let image = image::open(&output).expect("open png");
        assert_eq!(image.dimensions(), (80, 40));
    }
}
