//! Deterministic SVG-to-PNG rasterization helpers shared by headless tools.
//!
//! This module keeps the `resvg`-based conversion path in one reusable place so
//! wrapper/CLI/docs flows do not fork rendering behavior when they need a
//! messenger-friendly bitmap artifact from an engine-owned SVG export.

use regex::Regex;
use resvg::{self, tiny_skia, usvg};
use serde::Serialize;
use std::{
    collections::{BTreeMap, BTreeSet, HashSet},
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

/// Path-free identity of a font face selected for positioned SVG glyphs.
///
/// These are layout-used faces, not all loaded fonts or proof that every glyph
/// produced visible pixels. The inventory does not enforce font pinning or
/// establish complete glyph coverage/cross-host rendering reproducibility.
#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
pub struct SvgUsedFontIdentity {
    /// Actual family names from the selected face, sorted and deduplicated.
    pub families: Vec<String>,
    /// Actual PostScript name from the selected face.
    pub post_script_name: String,
    /// Face index within the source font file or collection.
    pub face_index: u32,
    /// Lowercase, unprefixed SHA-256 of the complete source/container bytes.
    /// The face index is recorded separately, not mixed into this digest.
    pub sha256: String,
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

    let rendered = render_svg_file_to_png_bytes(input_path, options)?;

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

/// Rasterizes an SVG file into deterministic PNG bytes while preserving the
/// source file's parent directory for relative resources.
pub fn render_svg_file_to_png_bytes(
    input_path: &Path,
    options: SvgPngRenderOptions,
) -> Result<SvgPngRenderBytes, String> {
    render_svg_file_to_png_bytes_impl(input_path, options, false).map(|(rendered, _)| rendered)
}

/// Rasterize once and bind the actual glyph-selected font faces from that tree.
///
/// Missing face metadata or unreadable font bytes fail the audited call. An SVG
/// without positioned glyphs returns an empty inventory, not all loaded fonts.
pub fn render_svg_file_to_png_bytes_audited(
    input_path: &Path,
    options: SvgPngRenderOptions,
) -> Result<(SvgPngRenderBytes, Vec<SvgUsedFontIdentity>), String> {
    render_svg_file_to_png_bytes_impl(input_path, options, true)
}

fn render_svg_file_to_png_bytes_impl(
    input_path: &Path,
    options: SvgPngRenderOptions,
    audit_fonts: bool,
) -> Result<(SvgPngRenderBytes, Vec<SvgUsedFontIdentity>), String> {
    if input_path.as_os_str().is_empty() {
        return Err("svg-png requires INPUT.svg".to_string());
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

    render_svg_text_to_png_bytes(
        &svg_text,
        canonical_parent(input_path),
        options,
        &format!("'{}'", input_path.display()),
        audit_fonts,
    )
}

/// Rasterizes an SVG string into deterministic PNG bytes.
pub fn render_svg_to_png_bytes(
    svg_text: &str,
    options: SvgPngRenderOptions,
) -> Result<SvgPngRenderBytes, String> {
    render_svg_to_png_bytes_impl(svg_text, options, false).map(|(rendered, _)| rendered)
}

/// Rasterize once, returning the unchanged PNG payload plus used font identities.
///
/// The audit covers positioned glyphs in the rendered tree and its referenced
/// subtrees. Any unreadable used face fails rather than yielding a partial audit.
pub fn render_svg_to_png_bytes_audited(
    svg_text: &str,
    options: SvgPngRenderOptions,
) -> Result<(SvgPngRenderBytes, Vec<SvgUsedFontIdentity>), String> {
    render_svg_to_png_bytes_impl(svg_text, options, true)
}

fn render_svg_to_png_bytes_impl(
    svg_text: &str,
    options: SvgPngRenderOptions,
    audit_fonts: bool,
) -> Result<(SvgPngRenderBytes, Vec<SvgUsedFontIdentity>), String> {
    let mut svg_text = svg_text.to_string();
    if options.drop_dotplot_metadata {
        svg_text = strip_dotplot_metadata_text(&svg_text);
    }
    render_svg_text_to_png_bytes(&svg_text, None, options, "input SVG", audit_fonts)
}

fn render_svg_text_to_png_bytes(
    svg_text: &str,
    resources_dir: Option<PathBuf>,
    options: SvgPngRenderOptions,
    input_label: &str,
    audit_fonts: bool,
) -> Result<(SvgPngRenderBytes, Vec<SvgUsedFontIdentity>), String> {
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
    render_svg_with_usvg_options(svg_text, options, input_label, &opt, audit_fonts)
}

fn render_svg_with_usvg_options(
    svg_text: &str,
    options: SvgPngRenderOptions,
    input_label: &str,
    opt: &usvg::Options<'_>,
    audit_fonts: bool,
) -> Result<(SvgPngRenderBytes, Vec<SvgUsedFontIdentity>), String> {
    let font_face_count = opt.fontdb.len();
    ensure_svg_text_fonts_available(svg_text, font_face_count)?;

    let tree = usvg::Tree::from_str(svg_text, opt)
        .map_err(|e| format!("Could not parse SVG {input_label}: {e}"))?;
    let fonts = if audit_fonts {
        collect_used_fonts(&tree)?
    } else {
        Vec::new()
    };
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

    Ok((
        SvgPngRenderBytes {
            bytes,
            width: pixmap_size.width(),
            height: pixmap_size.height(),
            font_face_count,
        },
        fonts,
    ))
}

fn used_font_identity(
    fontdb: &usvg::fontdb::Database,
    id: usvg::fontdb::ID,
) -> Result<SvgUsedFontIdentity, String> {
    let face = fontdb
        .face(id)
        .ok_or_else(|| "Could not audit SVG used font: face metadata is missing".to_string())?;
    let (sha256, face_index) = fontdb
        .with_face_data(id, |bytes, index| {
            (crate::digest_utils::sha256_hex_bytes(bytes), index)
        })
        .ok_or_else(|| "Could not audit SVG used font: font bytes are unreadable".to_string())?;
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
            // Subroots can themselves carry effects; follow complete clip/mask
            // chains, including effects on a tree's root group.
            let mut clip = group.clip_path();
            while let Some(current) = clip {
                self.group(current.root(), fontdb)?;
                clip = current.clip_path();
            }
            let mut mask = group.mask();
            while let Some(current) = mask {
                self.group(current.root(), fontdb)?;
                mask = current.mask();
            }
            for filter in group.filters() {
                for primitive in filter.primitives() {
                    if let usvg::filter::Kind::Image(image) = primitive.kind() {
                        self.group(image.root(), fontdb)?;
                    }
                }
            }
            for node in group.children() {
                match node {
                    usvg::Node::Group(group) => {
                        self.group(group, fontdb)?;
                        continue;
                    }
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
                            // IDs belong to a database, not globally to the outer SVG.
                            self.group(inner.root(), inner.fontdb())?;
                        }
                        continue;
                    }
                    usvg::Node::Path(_) => {}
                }
                let mut result = Ok(());
                node.subroots(|subroot| {
                    if result.is_ok() {
                        result = self.group(subroot, fontdb);
                    }
                });
                result?;
            }
            Ok(())
        }
    }

    let mut audit = Audit::default();
    audit.group(tree.root(), tree.fontdb())?;
    Ok(audit.identities.into_values().collect())
}

#[cfg(test)]
pub(crate) mod audit_test_support {
    use super::*;

    pub(crate) const ALPHA: &str = "GentleAuditAlpha";
    pub(crate) const BETA: &str = "GentleAuditBeta";

    fn u16_at(bytes: &mut [u8], offset: usize, value: u16) {
        bytes[offset..offset + 2].copy_from_slice(&value.to_be_bytes());
    }

    fn u32_at(bytes: &mut [u8], offset: usize, value: u32) {
        bytes[offset..offset + 4].copy_from_slice(&value.to_be_bytes());
    }

    fn checksum(bytes: &[u8]) -> u32 {
        bytes.chunks(4).fold(0_u32, |sum, chunk| {
            let mut word = [0; 4];
            word[..chunk.len()].copy_from_slice(chunk);
            sum.wrapping_add(u32::from_be_bytes(word))
        })
    }

    // Hand-crafted synthetic TrueType: an empty .notdef and one rectangular A.
    // No external font/artwork is copied. These literal tables regenerate the
    // exact test bytes independently of system fonts or platform font discovery.
    pub(crate) fn font(family: &str) -> Vec<u8> {
        let mut head = vec![0; 54];
        u32_at(&mut head, 0, 0x0001_0000);
        u32_at(&mut head, 4, 0x0001_0000);
        u32_at(&mut head, 12, 0x5f0f_3cf5);
        u16_at(&mut head, 18, 1_000);
        u16_at(&mut head, 40, 500);
        u16_at(&mut head, 42, 700);
        u16_at(&mut head, 46, 8);
        u16_at(&mut head, 48, 2);

        let mut hhea = vec![0; 36];
        u32_at(&mut hhea, 0, 0x0001_0000);
        for (offset, value) in [
            (4, 800),
            (6, -200),
            (10, 600),
            (14, 100),
            (16, 500),
            (18, 1),
            (34, 2),
        ] {
            u16_at(&mut hhea, offset, value as u16);
        }
        let mut hmtx = vec![0; 8];
        u16_at(&mut hmtx, 0, 600);
        u16_at(&mut hmtx, 4, 600);
        let mut maxp = vec![0; 32];
        u32_at(&mut maxp, 0, 0x0001_0000);
        for (offset, value) in [(4, 2), (6, 4), (8, 1), (14, 1)] {
            u16_at(&mut maxp, offset, value);
        }
        let mut glyf = vec![0; 36];
        for (offset, value) in [
            (0, 1),
            (6, 500),
            (8, 700),
            (10, 3),
            (20, 500),
            (24, -500),
            (30, 700),
        ] {
            u16_at(&mut glyf, offset, value as u16);
        }
        glyf[14..18].fill(1);
        let mut cmap = vec![0; 44];
        u32_at(&mut cmap, 8, 12);
        for (offset, value) in [
            (2, 1),
            (4, 3),
            (6, 1),
            (12, 4),
            (14, 32),
            (18, 4),
            (20, 4),
            (22, 1),
            (26, 65),
            (28, 65535),
            (32, 65),
            (34, 65535),
            (36, 65472),
            (38, 1),
        ] {
            u16_at(&mut cmap, offset, value);
        }
        let names = [
            (1, family.to_string()),
            (2, "Regular".into()),
            (4, format!("{family} Regular")),
            (6, format!("{family}-Regular")),
        ];
        let mut name = vec![0; 6 + names.len() * 12];
        u16_at(&mut name, 2, names.len() as u16);
        let string_offset = name.len();
        u16_at(&mut name, 4, string_offset as u16);
        for (index, (id, value)) in names.into_iter().enumerate() {
            let encoded: Vec<_> = value.encode_utf16().flat_map(u16::to_be_bytes).collect();
            let record = 6 + index * 12;
            for (offset, value) in [
                (0, 3),
                (2, 1),
                (4, 0x0409),
                (6, id),
                (8, encoded.len() as u16),
                (10, (name.len() - string_offset) as u16),
            ] {
                u16_at(&mut name, record + offset, value);
            }
            name.extend(encoded);
        }
        let mut post = vec![0; 32];
        u32_at(&mut post, 0, 0x0003_0000);
        u16_at(&mut post, 8, (-75_i16) as u16);
        u16_at(&mut post, 10, 50);
        u32_at(&mut post, 12, 1);
        let tables = BTreeMap::from([
            ("cmap", cmap),
            ("glyf", glyf),
            ("head", head),
            ("hhea", hhea),
            ("hmtx", hmtx),
            ("loca", vec![0, 0, 0, 0, 0, 18]),
            ("maxp", maxp),
            ("name", name),
            ("post", post),
        ]);
        let mut data = vec![0; 12 + tables.len() * 16];
        u32_at(&mut data, 0, 0x0001_0000);
        u16_at(&mut data, 4, tables.len() as u16);
        u16_at(&mut data, 6, 128);
        u16_at(&mut data, 8, 3);
        u16_at(&mut data, 10, 16);
        let mut head_offset = 0;
        for (index, (tag, bytes)) in tables.into_iter().enumerate() {
            let record = 12 + index * 16;
            let offset = data.len();
            data[record..record + 4].copy_from_slice(tag.as_bytes());
            u32_at(&mut data, record + 4, checksum(&bytes));
            u32_at(&mut data, record + 8, offset as u32);
            u32_at(&mut data, record + 12, bytes.len() as u32);
            if tag == "head" {
                head_offset = offset;
            }
            data.extend(bytes);
            while data.len() % 4 != 0 {
                data.push(0);
            }
        }
        let adjustment = 0xb1b0_afba_u32.wrapping_sub(checksum(&data));
        u32_at(&mut data, head_offset + 8, adjustment);
        data
    }

    pub(crate) fn collection() -> Vec<u8> {
        let mut data = vec![0; 20];
        data[..4].copy_from_slice(b"ttcf");
        u32_at(&mut data, 4, 0x0001_0000);
        u32_at(&mut data, 8, 2);
        for (index, family) in [ALPHA, BETA].into_iter().enumerate() {
            let offset = data.len() as u32;
            u32_at(&mut data, 12 + index * 4, offset);
            let mut face = font(family);
            let table_count = u16::from_be_bytes([face[4], face[5]]) as usize;
            for table in 0..table_count {
                let slot = 12 + table * 16 + 8;
                let original = u32::from_be_bytes(face[slot..slot + 4].try_into().unwrap());
                u32_at(&mut face, slot, original + offset);
            }
            data.extend(face);
        }
        data
    }

    pub(crate) fn options() -> usvg::Options<'static> {
        let mut options = usvg::Options {
            font_family: ALPHA.into(),
            ..usvg::Options::default()
        };
        let fontdb = options.fontdb_mut();
        fontdb.load_font_data(font(ALPHA));
        fontdb.load_font_data(font(BETA));
        fontdb.set_serif_family(ALPHA);
        fontdb.set_sans_serif_family(ALPHA);
        fontdb.set_monospace_family(ALPHA);
        assert_eq!(fontdb.len(), 2, "synthetic TrueType faces must load");
        options
    }

    pub(crate) fn render(
        svg: &str,
        options: SvgPngRenderOptions,
        audited: bool,
    ) -> Result<(SvgPngRenderBytes, Vec<SvgUsedFontIdentity>), String> {
        let svg = if options.drop_dotplot_metadata {
            strip_dotplot_metadata_text(svg)
        } else {
            svg.into()
        };
        render_svg_with_usvg_options(
            &svg,
            options,
            "synthetic font test SVG",
            &self::options(),
            audited,
        )
    }
}

#[cfg(test)]
mod tests {
    use super::{
        SvgPngRenderOptions, ensure_svg_text_fonts_available, render_svg_file_to_png,
        strip_dotplot_metadata_text,
    };
    use image::GenericImageView;

    #[test]
    fn audited_text_binds_only_selected_faces_and_matches_legacy_pixels() {
        use super::audit_test_support::{ALPHA, font, render};

        let svg = format!(
            r#"<svg xmlns="http://www.w3.org/2000/svg" width="80" height="40"><g font-family="{ALPHA}" font-size="18"><text y="18">AAAA</text><g><text y="38">AA</text></g></g></svg>"#
        );
        for scale in [1.0, 2.0] {
            let options = SvgPngRenderOptions {
                scale,
                drop_dotplot_metadata: false,
            };
            let (legacy, no_audit) = render(&svg, options, false).unwrap();
            let (audited, fonts) = render(&svg, options, true).unwrap();
            assert_eq!(legacy, audited);
            assert!(no_audit.is_empty());
            assert_eq!(audited.font_face_count, 2);
            assert_eq!(fonts.len(), 1);
            assert_eq!(fonts[0].families, vec![ALPHA.to_string()]);
            assert_eq!(fonts[0].post_script_name, format!("{ALPHA}-Regular"));
            assert_eq!(fonts[0].face_index, 0);
            assert_eq!(
                fonts[0].sha256,
                crate::digest_utils::sha256_hex_bytes(&font(ALPHA))
            );
            let pixels = image::load_from_memory(&audited.bytes).unwrap().to_rgba8();
            assert!(
                pixels.pixels().any(|pixel| pixel.0[3] != 0),
                "text must actually rasterize"
            );
            let metadata = serde_json::to_value(&fonts[0]).unwrap();
            assert_eq!(metadata.as_object().unwrap().len(), 4);
            assert!(metadata.get("path").is_none());
        }
    }

    #[test]
    fn audited_blank_and_stripped_text_do_not_report_loaded_fonts() {
        use super::audit_test_support::render;

        for svg in [
            r#"<svg xmlns="http://www.w3.org/2000/svg" width="20" height="20"><rect width="20" height="20"/></svg>"#,
            r#"<svg xmlns="http://www.w3.org/2000/svg" width="20" height="20"><defs><text id="unused" y="18">AAAA</text></defs></svg>"#,
            r#"<svg xmlns="http://www.w3.org/2000/svg" width="20" height="20"><text y="18">rendered_cells=AAAA</text></svg>"#,
        ] {
            let options = SvgPngRenderOptions {
                scale: 1.0,
                drop_dotplot_metadata: true,
            };
            let (legacy, _) = render(svg, options, false).unwrap();
            let (audited, fonts) = render(svg, options, true).unwrap();
            assert_eq!(legacy, audited);
            assert_eq!(audited.font_face_count, 2);
            assert!(fonts.is_empty());
        }
    }

    #[test]
    fn audited_fallback_uses_resolved_face_not_requested_family() {
        use super::audit_test_support::{ALPHA, render};

        let svg = r#"<svg xmlns="http://www.w3.org/2000/svg" width="40" height="20"><text y="18" font-family="DeliberatelyMissingFamily">AA</text></svg>"#;
        let (_, fonts) = render(svg, SvgPngRenderOptions::default(), true).unwrap();
        assert_eq!(fonts.len(), 1);
        assert_eq!(fonts[0].families, vec![ALPHA.to_string()]);
    }

    #[test]
    fn audited_pattern_and_chained_mask_subroots_include_their_text_fonts() {
        use super::audit_test_support::{ALPHA, BETA, render};

        let definitions = [
            format!(
                r#"<defs><pattern id="p" patternUnits="userSpaceOnUse" width="30" height="20"><text y="18" font-family="{ALPHA}">AA</text></pattern></defs><rect width="60" height="20" fill="url(#p)"/>"#
            ),
            format!(
                r#"<defs><mask id="m3"><text y="18" fill="white" font-family="{ALPHA}">AA</text></mask><mask id="m2" mask="url(#m3)"><rect width="60" height="20" fill="white"/></mask><mask id="m1" mask="url(#m2)"><rect width="60" height="20" fill="white"/></mask></defs><rect width="60" height="20" mask="url(#m1)"/>"#
            ),
        ];
        for definition in definitions {
            let svg = format!(
                r#"<svg xmlns="http://www.w3.org/2000/svg" width="60" height="40">{definition}<text y="38" font-family="{BETA}">AA</text></svg>"#
            );
            let (_, fonts) = render(&svg, SvgPngRenderOptions::default(), true).unwrap();
            let names: std::collections::BTreeSet<_> = fonts
                .iter()
                .flat_map(|font| font.families.iter().map(String::as_str))
                .collect();
            assert_eq!(names, std::collections::BTreeSet::from([ALPHA, BETA]));
        }
    }

    #[test]
    fn audited_nested_svg_uses_its_own_font_database() {
        use super::audit_test_support::{ALPHA, BETA, render};

        let inner = format!(
            r#"<svg xmlns="http://www.w3.org/2000/svg" width="40" height="20"><text y="18" font-family="{ALPHA}">AA</text></svg>"#
        );
        let encoded: String = inner.bytes().map(|byte| format!("%{byte:02X}")).collect();
        let svg = format!(
            r#"<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" width="40" height="40"><image width="40" height="20" xlink:href="data:image/svg+xml,{encoded}"/><text y="38" font-family="{BETA}">AA</text></svg>"#
        );
        let (_, fonts) = render(&svg, SvgPngRenderOptions::default(), true).unwrap();
        let names: std::collections::BTreeSet<_> = fonts
            .iter()
            .flat_map(|font| font.families.iter().map(String::as_str))
            .collect();
        assert_eq!(names, std::collections::BTreeSet::from([ALPHA, BETA]));
    }

    #[test]
    fn audited_collection_faces_bind_complete_bytes_and_distinct_face_indices() {
        use super::audit_test_support::{ALPHA, BETA, collection};

        let data = collection();
        let mut options = resvg::usvg::Options::default();
        options.fontdb_mut().load_font_data(data.clone());
        assert_eq!(options.fontdb.len(), 2);
        let svg = format!(
            r#"<svg xmlns="http://www.w3.org/2000/svg" width="40" height="40"><text y="18" font-family="{ALPHA}">AA</text><text y="38" font-family="{BETA}">AA</text></svg>"#
        );
        let (_, fonts) = super::render_svg_with_usvg_options(
            &svg,
            SvgPngRenderOptions::default(),
            "synthetic collection",
            &options,
            true,
        )
        .unwrap();
        assert_eq!(fonts.len(), 2);
        assert_eq!(fonts[0].face_index, 0);
        assert_eq!(fonts[1].face_index, 1);
        for font in fonts {
            assert_eq!(font.sha256, crate::digest_utils::sha256_hex_bytes(&data));
        }
    }

    #[test]
    fn audited_missing_or_unreadable_used_faces_fail_without_paths() {
        use super::audit_test_support::{ALPHA, font};

        let mut fontdb = resvg::usvg::fontdb::Database::new();
        assert!(super::used_font_identity(&fontdb, resvg::usvg::fontdb::ID::dummy()).is_err());
        let temp = tempfile::tempdir().unwrap();
        let path = temp.path().join("synthetic.ttf");
        std::fs::write(&path, font(ALPHA)).unwrap();
        fontdb.load_font_file(&path).unwrap();
        let id = fontdb.faces().next().unwrap().id;
        std::fs::remove_file(&path).unwrap();
        let error = super::used_font_identity(&fontdb, id).unwrap_err();
        assert!(error.contains("unreadable"));
        assert!(!error.contains(path.to_str().unwrap()));
    }

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
