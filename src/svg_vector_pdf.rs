//! Audited static SVG-to-vector-PDF conversion.
//!
//! This is deliberately separate from [`crate::svg_pdf`], whose established
//! contract is a lossless raster-backed PDF. Vector PDFs preserve paths and
//! selectable embedded text, but SVG interaction such as `<title>` hover does
//! not have a portable PDF equivalent. Callers should retain the SVG/HTML peer
//! for interaction and machine-readable data for complete evidence.

use crate::{
    digest_utils::sha256_hex_bytes,
    svg_png::{
        SVG_FONT_DIR_ENV, SVG_FONT_FILE_ENV, SVG_MONOSPACE_FAMILY_ENV, SVG_SANS_SERIF_FAMILY_ENV,
        SVG_SERIF_FAMILY_ENV, SvgPngRenderOptions, SvgUsedFontIdentity,
        strip_dotplot_metadata_text,
    },
};
use serde::Serialize;
use std::{
    collections::{BTreeMap, BTreeSet, HashSet},
    env,
    path::PathBuf,
};
use svg2pdf::usvg;

/// In-memory vector PDF plus auditable rendering facts.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct SvgVectorPdfRenderBytes {
    /// Complete PDF bytes.
    #[serde(skip)]
    pub bytes: Vec<u8>,
    /// SVG width after applying the requested scale, in CSS pixels.
    pub width: u32,
    /// SVG height after applying the requested scale, in CSS pixels.
    pub height: u32,
    /// PDF media-box width in points.
    pub page_width_pt: String,
    /// PDF media-box height in points.
    pub page_height_pt: String,
    /// All font faces available to the vector converter.
    pub font_face_count: usize,
    /// Positioned-glyph font faces actually selected by the parsed SVG tree.
    pub font_identities: Vec<SvgUsedFontIdentity>,
    /// Whether PDF text is embedded as selectable text rather than paths.
    pub embedded_text: bool,
}

fn configured_font_paths(var_name: &str) -> Vec<PathBuf> {
    env::var_os(var_name)
        .map(|raw| env::split_paths(&raw).collect())
        .unwrap_or_default()
}

fn load_configured_fonts(fontdb: &mut usvg::fontdb::Database) -> Result<(), String> {
    for path in configured_font_paths(SVG_FONT_FILE_ENV) {
        fontdb.load_font_file(&path).map_err(|error| {
            format!(
                "Could not load SVG vector-PDF font file from {SVG_FONT_FILE_ENV} '{}': {error}",
                path.display()
            )
        })?;
    }
    for path in configured_font_paths(SVG_FONT_DIR_ENV) {
        fontdb.load_fonts_dir(&path);
    }
    Ok(())
}

fn has_family(fontdb: &usvg::fontdb::Database, family: &str) -> bool {
    fontdb.faces().any(|face| {
        face.families
            .iter()
            .any(|(name, _)| name.eq_ignore_ascii_case(family))
    })
}

fn first_family(fontdb: &usvg::fontdb::Database) -> Option<String> {
    fontdb
        .faces()
        .flat_map(|face| face.families.iter().map(|(name, _)| name.trim()))
        .find(|name| !name.is_empty())
        .map(ToOwned::to_owned)
}

fn choose_family(
    fontdb: &usvg::fontdb::Database,
    env_name: &str,
    candidates: &[&str],
) -> Option<String> {
    if let Ok(value) = env::var(env_name) {
        let requested = value.trim();
        if !requested.is_empty() && has_family(fontdb, requested) {
            return Some(requested.to_string());
        }
    }
    candidates
        .iter()
        .find(|family| has_family(fontdb, family))
        .map(|family| (*family).to_string())
        .or_else(|| first_family(fontdb))
}

fn configure_generic_families(fontdb: &mut usvg::fontdb::Database) {
    if let Some(family) = choose_family(
        fontdb,
        SVG_MONOSPACE_FAMILY_ENV,
        &["DejaVu Sans Mono", "Liberation Mono", "Noto Sans Mono"],
    ) {
        fontdb.set_monospace_family(family);
    }
    if let Some(family) = choose_family(
        fontdb,
        SVG_SANS_SERIF_FAMILY_ENV,
        &["DejaVu Sans", "Liberation Sans", "Noto Sans", "Arial"],
    ) {
        fontdb.set_sans_serif_family(family);
    }
    if let Some(family) = choose_family(
        fontdb,
        SVG_SERIF_FAMILY_ENV,
        &["DejaVu Serif", "Liberation Serif", "Noto Serif", "Times"],
    ) {
        fontdb.set_serif_family(family);
    }
}

fn used_font_identity(
    fontdb: &usvg::fontdb::Database,
    id: usvg::fontdb::ID,
) -> Result<SvgUsedFontIdentity, String> {
    let face = fontdb.face(id).ok_or_else(|| {
        "Could not audit vector-PDF used font: face metadata is missing".to_string()
    })?;
    let (sha256, face_index) = fontdb
        .with_face_data(id, |bytes, index| (sha256_hex_bytes(bytes), index))
        .ok_or_else(|| {
            "Could not audit vector-PDF used font: font bytes are unreadable".to_string()
        })?;
    Ok(SvgUsedFontIdentity {
        families: face
            .families
            .iter()
            .map(|(name, _)| name.clone())
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect(),
        post_script_name: face.post_script_name.clone(),
        face_index,
        sha256,
    })
}

fn collect_used_fonts(tree: &usvg::Tree) -> Result<Vec<SvgUsedFontIdentity>, String> {
    #[derive(Default)]
    struct Audit {
        groups: HashSet<(*const usvg::Group, *const usvg::fontdb::Database)>,
        faces: HashSet<(*const usvg::fontdb::Database, usvg::fontdb::ID)>,
        identities: BTreeMap<(String, u32), SvgUsedFontIdentity>,
    }

    impl Audit {
        fn group(
            &mut self,
            group: &usvg::Group,
            fontdb: &usvg::fontdb::Database,
        ) -> Result<(), String> {
            if !self.groups.insert((group, fontdb)) {
                return Ok(());
            }
            for node in group.children() {
                match node {
                    usvg::Node::Group(group) => self.group(group, fontdb)?,
                    usvg::Node::Text(text) => {
                        for span in text.layouted() {
                            for glyph in &span.positioned_glyphs {
                                if self.faces.insert((fontdb, glyph.font)) {
                                    let identity = used_font_identity(fontdb, glyph.font)?;
                                    self.identities
                                        .entry((identity.sha256.clone(), identity.face_index))
                                        .or_insert(identity);
                                }
                            }
                        }
                    }
                    usvg::Node::Image(image) => {
                        if let usvg::ImageKind::SVG(inner) = image.kind() {
                            self.group(inner.root(), inner.fontdb())?;
                        }
                    }
                    usvg::Node::Path(_) => {}
                }
                let mut nested = Ok(());
                node.subroots(|subroot| {
                    if nested.is_ok() {
                        nested = self.group(subroot, fontdb);
                    }
                });
                nested?;
            }
            Ok(())
        }
    }

    let mut audit = Audit::default();
    audit.group(tree.root(), tree.fontdb())?;
    Ok(audit.identities.into_values().collect())
}

/// Convert an SVG string to a static vector PDF with embedded selectable text.
///
/// SVG hover/title interaction is intentionally not claimed. The SVG peer is
/// still the interactive representation; this PDF is the printable/searchable
/// projection of the same geometry.
pub fn render_svg_to_vector_pdf_bytes_audited(
    svg_text: &str,
    options: SvgPngRenderOptions,
) -> Result<SvgVectorPdfRenderBytes, String> {
    if !(options.scale.is_finite() && options.scale > 0.0) {
        return Err(format!(
            "svg-vector-pdf requires a positive finite scale value, got {}",
            options.scale
        ));
    }
    let svg = if options.drop_dotplot_metadata {
        strip_dotplot_metadata_text(svg_text)
    } else {
        svg_text.to_string()
    };
    let mut usvg_options = usvg::Options::default();
    {
        let fontdb = usvg_options.fontdb_mut();
        fontdb.load_system_fonts();
        load_configured_fonts(fontdb)?;
        configure_generic_families(fontdb);
    }
    let font_face_count = usvg_options.fontdb.len();
    if font_face_count == 0 && svg.to_ascii_lowercase().contains("<text") {
        return Err("SVG contains text, but no font faces are available for vector PDF".into());
    }
    let tree = usvg::Tree::from_str(&svg, &usvg_options)
        .map_err(|error| format!("Could not parse SVG for vector PDF: {error}"))?;
    let font_identities = collect_used_fonts(&tree)?;
    let size = tree
        .size()
        .to_int_size()
        .scale_by(options.scale)
        .ok_or_else(|| format!("Could not scale SVG size by {}", options.scale))?;
    let page_width_pt = size.width() as f32 * 72.0 / 96.0;
    let page_height_pt = size.height() as f32 * 72.0 / 96.0;
    let conversion = svg2pdf::ConversionOptions {
        compress: true,
        raster_scale: 1.5 * options.scale,
        embed_text: true,
        pdfa: false,
    };
    let page = svg2pdf::PageOptions {
        dpi: 96.0 / options.scale,
    };
    let bytes = svg2pdf::to_pdf(&tree, conversion, page)
        .map_err(|error| format!("Could not convert SVG to vector PDF: {error}"))?;
    if !bytes.starts_with(b"%PDF-") {
        return Err("vector PDF backend returned a non-PDF payload".into());
    }
    Ok(SvgVectorPdfRenderBytes {
        bytes,
        width: size.width(),
        height: size.height(),
        page_width_pt: format!("{page_width_pt:.2}"),
        page_height_pt: format!("{page_height_pt:.2}"),
        font_face_count,
        font_identities,
        embedded_text: true,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn vector_pdf_keeps_vector_geometry_and_embedded_text_policy() {
        let svg = r#"<svg xmlns="http://www.w3.org/2000/svg" width="200" height="100"><rect x="5" y="5" width="190" height="90" fill="none" stroke="black"/><text x="20" y="55" font-family="sans-serif" font-size="18">Selectable TFBS</text></svg>"#;
        let rendered =
            render_svg_to_vector_pdf_bytes_audited(svg, SvgPngRenderOptions::default()).unwrap();
        assert!(rendered.bytes.starts_with(b"%PDF-"));
        assert_eq!((rendered.width, rendered.height), (200, 100));
        assert_eq!(rendered.page_width_pt, "150.00");
        assert_eq!(rendered.page_height_pt, "75.00");
        assert!(rendered.embedded_text);
        assert!(!rendered.font_identities.is_empty());
        assert!(
            !rendered
                .bytes
                .windows(b"/Subtype /Image".len())
                .any(|window| window == b"/Subtype /Image")
        );
    }

    #[test]
    fn vector_pdf_rejects_invalid_scale() {
        let error = render_svg_to_vector_pdf_bytes_audited(
            "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"1\" height=\"1\"/>",
            SvgPngRenderOptions {
                scale: 0.0,
                drop_dotplot_metadata: false,
            },
        )
        .unwrap_err();
        assert!(error.contains("positive finite scale"));
    }
}
