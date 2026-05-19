//! Deterministic SVG-to-PNG rasterization helpers shared by headless tools.
//!
//! This module keeps the `resvg`-based conversion path in one reusable place so
//! wrapper/CLI/docs flows do not fork rendering behavior when they need a
//! messenger-friendly bitmap artifact from an engine-owned SVG export.

use regex::Regex;
use resvg::{self, tiny_skia, usvg};
use serde::Serialize;
use std::{
    env,
    path::{Path, PathBuf},
    sync::LazyLock,
};

static DOTPLOT_METADATA_RE: LazyLock<Regex> = LazyLock::new(|| {
    Regex::new(
        r#"<text\b[^>]*>(?:Dotplot workspace export:[^<]*|rendered_cells=[^<]*|x:\s[^<]*|y:\s[^<]*|GENtle dotplot SVG export)</text>"#,
    )
    .expect("dotplot metadata regex should compile")
});
static SVG_TEXT_RE: LazyLock<Regex> =
    LazyLock::new(|| Regex::new(r"(?i)<text\b").expect("SVG text regex should compile"));

/// Optional path-list environment variable for explicit SVG rasterization font files.
pub const SVG_FONT_FILE_ENV: &str = "GENTLE_SVG_FONT_FILE";
/// Optional path-list environment variable for explicit SVG rasterization font directories.
pub const SVG_FONT_DIR_ENV: &str = "GENTLE_SVG_FONT_DIR";
/// Optional family name used for SVG generic `monospace` text.
pub const SVG_MONOSPACE_FAMILY_ENV: &str = "GENTLE_SVG_MONOSPACE_FAMILY";
/// Optional family name used for SVG generic `sans-serif` text.
pub const SVG_SANS_SERIF_FAMILY_ENV: &str = "GENTLE_SVG_SANS_SERIF_FAMILY";
/// Optional family name used for SVG generic `serif` text.
pub const SVG_SERIF_FAMILY_ENV: &str = "GENTLE_SVG_SERIF_FAMILY";

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
    /// Number of font faces visible to `usvg`/`resvg` for text rendering.
    pub font_face_count: usize,
}

/// In-memory PNG payload produced from an SVG string.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SvgPngRenderBytes {
    /// PNG bytes encoded by `tiny-skia`.
    pub bytes: Vec<u8>,
    /// Written PNG width in pixels.
    pub width: u32,
    /// Written PNG height in pixels.
    pub height: u32,
    /// Number of font faces visible to `usvg`/`resvg` for text rendering.
    pub font_face_count: usize,
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

fn configured_font_paths(var_name: &str) -> Vec<PathBuf> {
    env::var_os(var_name)
        .map(|raw| env::split_paths(&raw).collect())
        .unwrap_or_default()
}

fn load_configured_fonts(fontdb: &mut usvg::fontdb::Database) -> Result<(), String> {
    for path in configured_font_paths(SVG_FONT_FILE_ENV) {
        fontdb.load_font_file(&path).map_err(|e| {
            format!(
                "Could not load SVG rasterization font file from {SVG_FONT_FILE_ENV} '{}': {e}",
                path.display()
            )
        })?;
    }
    for path in configured_font_paths(SVG_FONT_DIR_ENV) {
        fontdb.load_fonts_dir(&path);
    }
    Ok(())
}

fn fontdb_has_family(fontdb: &usvg::fontdb::Database, family: &str) -> bool {
    fontdb.faces().any(|face| {
        face.families
            .iter()
            .any(|(name, _)| name.eq_ignore_ascii_case(family))
    })
}

fn first_font_family(fontdb: &usvg::fontdb::Database) -> Option<String> {
    fontdb
        .faces()
        .flat_map(|face| face.families.iter().map(|(name, _)| name.trim()))
        .find(|name| !name.is_empty())
        .map(ToOwned::to_owned)
}

fn choose_font_family(
    fontdb: &usvg::fontdb::Database,
    env_name: &str,
    candidates: &[&str],
) -> Option<String> {
    if let Ok(value) = env::var(env_name) {
        let requested = value.trim();
        if !requested.is_empty() && fontdb_has_family(fontdb, requested) {
            return Some(requested.to_string());
        }
    }
    candidates
        .iter()
        .find(|family| fontdb_has_family(fontdb, family))
        .map(|family| (*family).to_string())
        .or_else(|| first_font_family(fontdb))
}

fn configure_generic_font_families(fontdb: &mut usvg::fontdb::Database) {
    if let Some(family) = choose_font_family(
        fontdb,
        SVG_MONOSPACE_FAMILY_ENV,
        &[
            "DejaVu Sans Mono",
            "Liberation Mono",
            "Noto Sans Mono",
            "Courier New",
            "Menlo",
            "Monaco",
            "Consolas",
        ],
    ) {
        fontdb.set_monospace_family(family);
    }
    if let Some(family) = choose_font_family(
        fontdb,
        SVG_SANS_SERIF_FAMILY_ENV,
        &[
            "DejaVu Sans",
            "Liberation Sans",
            "Noto Sans",
            "Arial",
            "Helvetica",
        ],
    ) {
        fontdb.set_sans_serif_family(family);
    }
    if let Some(family) = choose_font_family(
        fontdb,
        SVG_SERIF_FAMILY_ENV,
        &[
            "DejaVu Serif",
            "Liberation Serif",
            "Noto Serif",
            "Times New Roman",
            "Times",
        ],
    ) {
        fontdb.set_serif_family(family);
    }
}

fn ensure_svg_text_fonts_available(svg_text: &str, font_face_count: usize) -> Result<(), String> {
    if font_face_count == 0 && SVG_TEXT_RE.is_match(svg_text) {
        return Err(format!(
            "SVG contains text, but no font faces were available to usvg/resvg. \
Install a system font package (for example fonts-dejavu-core or fonts-liberation), \
or set {SVG_FONT_FILE_ENV} / {SVG_FONT_DIR_ENV} to a readable TTF/OTF font or font directory."
        ));
    }
    Ok(())
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

    let rendered = render_svg_text_to_png_bytes(
        &svg_text,
        canonical_parent(input_path),
        options,
        &format!("'{}'", input_path.display()),
    )?;

    std::fs::write(output_path, &rendered.bytes)
        .map_err(|e| format!("Could not write PNG '{}': {e}", output_path.display()))?;

    Ok(SvgPngRenderSummary {
        input_path: input_path.to_string_lossy().into_owned(),
        output_path: output_path.to_string_lossy().into_owned(),
        scale: format!("{}", options.scale),
        drop_dotplot_metadata: options.drop_dotplot_metadata,
        width: rendered.width,
        height: rendered.height,
        font_face_count: rendered.font_face_count,
    })
}

/// Rasterizes an SVG string into deterministic PNG bytes.
pub fn render_svg_to_png_bytes(
    svg_text: &str,
    options: SvgPngRenderOptions,
) -> Result<SvgPngRenderBytes, String> {
    let mut svg_text = svg_text.to_string();
    if options.drop_dotplot_metadata {
        svg_text = strip_dotplot_metadata_text(&svg_text);
    }
    render_svg_text_to_png_bytes(&svg_text, None, options, "input SVG")
}

fn render_svg_text_to_png_bytes(
    svg_text: &str,
    resources_dir: Option<PathBuf>,
    options: SvgPngRenderOptions,
    input_label: &str,
) -> Result<SvgPngRenderBytes, String> {
    if !(options.scale.is_finite() && options.scale > 0.0) {
        return Err(format!(
            "svg-png requires a positive finite scale value, got {}",
            options.scale
        ));
    }

    let mut opt = usvg::Options {
        resources_dir,
        ..usvg::Options::default()
    };
    {
        let fontdb = opt.fontdb_mut();
        fontdb.load_system_fonts();
        load_configured_fonts(fontdb)?;
        configure_generic_font_families(fontdb);
    }
    let font_face_count = opt.fontdb.len();
    ensure_svg_text_fonts_available(svg_text, font_face_count)?;

    let tree = usvg::Tree::from_str(&svg_text, &opt)
        .map_err(|e| format!("Could not parse SVG {input_label}: {e}"))?;
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
    let bytes = pixmap
        .encode_png()
        .map_err(|e| format!("Could not encode PNG from SVG {input_label}: {e}"))?;

    Ok(SvgPngRenderBytes {
        bytes,
        width: pixmap_size.width(),
        height: pixmap_size.height(),
        font_face_count,
    })
}

#[cfg(test)]
mod tests {
    use super::{
        SvgPngRenderOptions, ensure_svg_text_fonts_available, render_svg_file_to_png,
        strip_dotplot_metadata_text,
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

    #[test]
    fn svg_text_requires_available_font_faces() {
        let err = ensure_svg_text_fonts_available(
            "<svg xmlns=\"http://www.w3.org/2000/svg\"><text>Hello</text></svg>",
            0,
        )
        .expect_err("SVG text without fonts should fail early");
        assert!(err.contains("no font faces"));
        assert!(err.contains(super::SVG_FONT_FILE_ENV));
        assert!(ensure_svg_text_fonts_available(
            "<svg xmlns=\"http://www.w3.org/2000/svg\"><rect width=\"10\" height=\"10\"/></svg>",
            0,
        )
        .is_ok());
    }
}
