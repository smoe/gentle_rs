//! Deterministic SVG-to-PDF conversion helpers shared by headless exports.
//!
//! The PDF path intentionally reuses GENtle's existing `resvg` rasterization
//! helper, then embeds the losslessly compressed RGB image into one PDF page. That
//! keeps PDF generation dependency-light and visually consistent with SVG/PNG
//! exports.

use flate2::{Compression, write::ZlibEncoder};
use serde::Serialize;
use std::{io::Write, path::Path};

use crate::svg_png::{
    SvgPngRenderBytes, SvgPngRenderOptions, SvgUsedFontIdentity, render_svg_file_to_png_bytes,
    render_svg_file_to_png_bytes_audited,
};

/// Machine-readable summary of one deterministic SVG-to-PDF conversion.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct SvgPdfRenderSummary {
    /// Source SVG path that was rasterized.
    pub input_path: String,
    /// Output PDF path that was written.
    pub output_path: String,
    /// Applied raster scale factor.
    pub scale: String,
    /// Whether dotplot metadata stripping was enabled.
    pub drop_dotplot_metadata: bool,
    /// Embedded image width in pixels.
    pub width: u32,
    /// Embedded image height in pixels.
    pub height: u32,
    /// Number of font faces visible to `usvg`/`resvg` for text rendering.
    pub font_face_count: usize,
    /// PDF media-box width in points.
    pub page_width_pt: String,
    /// PDF media-box height in points.
    pub page_height_pt: String,
    /// Number of validated URI link annotations written to the PDF page.
    pub uri_link_count: usize,
}

/// Machine-readable summary of one page in a deterministic multi-page PDF.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct SvgPdfSetPageSummary {
    pub input_path: String,
    pub width: u32,
    pub height: u32,
    pub font_face_count: usize,
    pub page_width_pt: String,
    pub page_height_pt: String,
}

/// Machine-readable summary of a deterministic multi-page SVG-to-PDF conversion.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct SvgPdfSetRenderSummary {
    pub output_path: String,
    pub scale: String,
    pub drop_dotplot_metadata: bool,
    pub page_count: usize,
    pub pages: Vec<SvgPdfSetPageSummary>,
}

/// One URI hotspot expressed in the source SVG's CSS-pixel coordinate system.
#[derive(Debug, Clone, PartialEq)]
pub struct SvgPdfUriLink {
    pub x_px: f32,
    pub y_px: f32,
    pub width_px: f32,
    pub height_px: f32,
    pub uri: String,
}

#[derive(Debug, Clone, PartialEq)]
struct PdfUriLink {
    left_pt: f32,
    bottom_pt: f32,
    right_pt: f32,
    top_pt: f32,
    uri: String,
}

/// Renders one SVG file into a single-page PDF.
pub fn render_svg_file_to_pdf(
    input_path: &Path,
    output_path: &Path,
    options: SvgPngRenderOptions,
) -> Result<SvgPdfRenderSummary, String> {
    render_svg_file_to_pdf_with_links(input_path, output_path, options, &[], &[])
}

/// Renders ordered SVG files into one lossless raster-backed multi-page PDF.
///
/// Every input remains a separate page at its own dimensions. This is intended
/// for composite scientific reports whose context and detailed panels have
/// different page heights; it never concatenates or rescales pages to a common
/// canvas.
pub fn render_svg_files_to_pdf(
    input_paths: &[&Path],
    output_path: &Path,
    options: SvgPngRenderOptions,
) -> Result<SvgPdfSetRenderSummary, String> {
    if input_paths.is_empty() {
        return Err("svg-pdf-set requires at least one INPUT.svg".to_string());
    }
    if output_path.as_os_str().is_empty() {
        return Err("svg-pdf-set requires OUTPUT.pdf".to_string());
    }
    if !(options.scale.is_finite() && options.scale > 0.0) {
        return Err(format!(
            "svg-pdf-set requires a positive finite scale value, got {}",
            options.scale
        ));
    }

    let mut rendered = Vec::with_capacity(input_paths.len());
    let mut pages = Vec::with_capacity(input_paths.len());
    for input_path in input_paths {
        if input_path.as_os_str().is_empty() {
            return Err("svg-pdf-set input paths must not be empty".to_string());
        }
        let page = render_svg_file_to_png_bytes(input_path, options)?;
        let page_width_pt = page.width as f32 * 72.0 / 96.0;
        let page_height_pt = page.height as f32 * 72.0 / 96.0;
        pages.push(SvgPdfSetPageSummary {
            input_path: input_path.to_string_lossy().into_owned(),
            width: page.width,
            height: page.height,
            font_face_count: page.font_face_count,
            page_width_pt: format!("{page_width_pt:.2}"),
            page_height_pt: format!("{page_height_pt:.2}"),
        });
        rendered.push(page.bytes);
    }
    let pdf = png_pages_to_pdf(&rendered)?;
    std::fs::write(output_path, pdf)
        .map_err(|e| format!("Could not write PDF '{}': {e}", output_path.display()))?;
    Ok(SvgPdfSetRenderSummary {
        output_path: output_path.to_string_lossy().into_owned(),
        scale: format!("{}", options.scale),
        drop_dotplot_metadata: options.drop_dotplot_metadata,
        page_count: pages.len(),
        pages,
    })
}

/// Render a raster-backed PDF plus path-free identities of layout-used fonts.
///
/// Font evidence comes from the same SVG-to-PNG conversion whose pixels are
/// embedded into the PDF, not a second parse. Unreadable used fonts fail before
/// writing the PDF; no positioned glyphs means an empty font inventory. The
/// existing summary and PDF bytes retain their legacy shape and defaults.
pub fn render_svg_file_to_pdf_audited(
    input_path: &Path,
    output_path: &Path,
    options: SvgPngRenderOptions,
) -> Result<(SvgPdfRenderSummary, Vec<SvgUsedFontIdentity>), String> {
    render_svg_file_to_pdf_impl(
        input_path,
        output_path,
        options,
        &[],
        &[],
        render_svg_file_to_png_bytes_audited,
    )
}

/// Renders one SVG file into a single-page raster-backed PDF with validated
/// URI annotations. `allowed_hosts` uses exact host matching.
pub fn render_svg_file_to_pdf_with_links(
    input_path: &Path,
    output_path: &Path,
    options: SvgPngRenderOptions,
    links: &[SvgPdfUriLink],
    allowed_hosts: &[&str],
) -> Result<SvgPdfRenderSummary, String> {
    render_svg_file_to_pdf_impl(
        input_path,
        output_path,
        options,
        links,
        allowed_hosts,
        |input, options| {
            render_svg_file_to_png_bytes(input, options).map(|rendered| (rendered, Vec::new()))
        },
    )
    .map(|(summary, _)| summary)
}

fn render_svg_file_to_pdf_impl(
    input_path: &Path,
    output_path: &Path,
    options: SvgPngRenderOptions,
    links: &[SvgPdfUriLink],
    allowed_hosts: &[&str],
    render_png: impl FnOnce(
        &Path,
        SvgPngRenderOptions,
    ) -> Result<(SvgPngRenderBytes, Vec<SvgUsedFontIdentity>), String>,
) -> Result<(SvgPdfRenderSummary, Vec<SvgUsedFontIdentity>), String> {
    if input_path.as_os_str().is_empty() {
        return Err("svg-pdf requires INPUT.svg".to_string());
    }
    if output_path.as_os_str().is_empty() {
        return Err("svg-pdf requires OUTPUT.pdf".to_string());
    }
    if !(options.scale.is_finite() && options.scale > 0.0) {
        return Err(format!(
            "svg-pdf requires a positive finite scale value, got {}",
            options.scale
        ));
    }

    for link in links {
        validate_https_uri(&link.uri, allowed_hosts)?;
        if !link.x_px.is_finite()
            || !link.y_px.is_finite()
            || !link.width_px.is_finite()
            || !link.height_px.is_finite()
            || link.x_px < 0.0
            || link.y_px < 0.0
            || link.width_px <= 0.0
            || link.height_px <= 0.0
        {
            return Err("svg-pdf URI link rectangles must be finite and positive".to_string());
        }
    }
    let (rendered_bytes, fonts) = render_png(input_path, options)?;
    let page_width_pt = rendered_bytes.width as f32 * 72.0 / 96.0;
    let page_height_pt = rendered_bytes.height as f32 * 72.0 / 96.0;
    let points_per_source_px = options.scale * 72.0 / 96.0;
    let mut pdf_links = Vec::with_capacity(links.len());
    for link in links {
        let left = link.x_px * points_per_source_px;
        let right = (link.x_px + link.width_px) * points_per_source_px;
        let top = page_height_pt - link.y_px * points_per_source_px;
        let bottom = page_height_pt - (link.y_px + link.height_px) * points_per_source_px;
        if left < 0.0 || bottom < 0.0 || right > page_width_pt || top > page_height_pt {
            return Err(format!(
                "svg-pdf URI link rectangle for '{}' lies outside the rendered page",
                link.uri
            ));
        }
        pdf_links.push(PdfUriLink {
            left_pt: left,
            bottom_pt: bottom,
            right_pt: right,
            top_pt: top,
            uri: link.uri.clone(),
        });
    }
    let pdf = png_bytes_to_single_page_pdf(&rendered_bytes.bytes, &pdf_links)?;
    std::fs::write(output_path, pdf)
        .map_err(|e| format!("Could not write PDF '{}': {e}", output_path.display()))?;

    Ok((
        SvgPdfRenderSummary {
            input_path: input_path.to_string_lossy().into_owned(),
            output_path: output_path.to_string_lossy().into_owned(),
            scale: format!("{}", options.scale),
            drop_dotplot_metadata: options.drop_dotplot_metadata,
            width: rendered_bytes.width,
            height: rendered_bytes.height,
            font_face_count: rendered_bytes.font_face_count,
            page_width_pt: format!("{page_width_pt:.2}"),
            page_height_pt: format!("{page_height_pt:.2}"),
            uri_link_count: pdf_links.len(),
        },
        fonts,
    ))
}

fn validate_https_uri(uri: &str, allowed_hosts: &[&str]) -> Result<(), String> {
    if !uri.is_ascii() || uri.chars().any(char::is_control) {
        return Err("svg-pdf URI links must contain printable ASCII only".to_string());
    }
    let Some(rest) = uri.strip_prefix("https://") else {
        return Err(format!("svg-pdf URI link must use HTTPS: '{uri}'"));
    };
    let authority = rest.split('/').next().unwrap_or_default();
    if authority.is_empty()
        || authority.contains('@')
        || authority.contains(':')
        || !allowed_hosts.iter().any(|host| *host == authority)
    {
        return Err(format!(
            "svg-pdf URI link host '{authority}' is not explicitly allowed"
        ));
    }
    Ok(())
}

fn escape_pdf_literal(value: &str) -> String {
    let mut escaped = String::with_capacity(value.len());
    for byte in value.bytes() {
        match byte {
            b'\\' => escaped.push_str("\\\\"),
            b'(' => escaped.push_str("\\("),
            b')' => escaped.push_str("\\)"),
            0x20..=0x7e => escaped.push(byte as char),
            _ => escaped.push_str(&format!("\\{:03o}", byte)),
        }
    }
    escaped
}

fn png_bytes_to_single_page_pdf(png_bytes: &[u8], links: &[PdfUriLink]) -> Result<Vec<u8>, String> {
    let image = image::load_from_memory(png_bytes)
        .map_err(|e| format!("Could not decode rendered PNG for PDF embedding: {e}"))?
        .to_rgb8();
    let (width, height) = image.dimensions();
    let page_width_pt = width as f32 * 72.0 / 96.0;
    let page_height_pt = height as f32 * 72.0 / 96.0;
    // PDF's FlateDecode expects a zlib stream, not raw DEFLATE or gzip. Keep
    // every pixel and the original resolution; tall scientific pages compress
    // well without introducing JPEG artifacts or changing their geometry.
    let mut encoder = ZlibEncoder::new(Vec::new(), Compression::default());
    encoder
        .write_all(image.as_raw())
        .map_err(|e| format!("Could not compress PDF image: {e}"))?;
    let compressed_rgb = encoder
        .finish()
        .map_err(|e| format!("Could not finish PDF image compression: {e}"))?;
    let content = format!("q\n{page_width_pt:.2} 0 0 {page_height_pt:.2} 0 0 cm\n/Im0 Do\nQ\n");

    let mut pdf = Vec::new();
    pdf.extend_from_slice(b"%PDF-1.4\n");
    let mut offsets = Vec::new();
    push_pdf_object(
        &mut pdf,
        &mut offsets,
        b"<< /Type /Catalog /Pages 2 0 R >>".as_slice(),
    );
    push_pdf_object(
        &mut pdf,
        &mut offsets,
        b"<< /Type /Pages /Kids [3 0 R] /Count 1 >>".as_slice(),
    );
    let annotations = if links.is_empty() {
        String::new()
    } else {
        format!(
            " /Annots [{}]",
            (0..links.len())
                .map(|index| format!("{} 0 R", 6 + index))
                .collect::<Vec<_>>()
                .join(" ")
        )
    };
    push_pdf_object(
        &mut pdf,
        &mut offsets,
        format!(
            "<< /Type /Page /Parent 2 0 R /MediaBox [0 0 {page_width_pt:.2} {page_height_pt:.2}] /Resources << /XObject << /Im0 4 0 R >> >> /Contents 5 0 R{annotations} >>"
        )
        .as_bytes(),
    );
    push_pdf_stream_object(
        &mut pdf,
        &mut offsets,
        format!(
            "<< /Type /XObject /Subtype /Image /Width {width} /Height {height} /ColorSpace /DeviceRGB /BitsPerComponent 8 /Filter /FlateDecode /Length {} >>",
            compressed_rgb.len()
        )
        .as_bytes(),
        &compressed_rgb,
    );
    push_pdf_stream_object(
        &mut pdf,
        &mut offsets,
        format!("<< /Length {} >>", content.len()).as_bytes(),
        content.as_bytes(),
    );
    for link in links {
        push_pdf_object(
            &mut pdf,
            &mut offsets,
            format!(
                "<< /Type /Annot /Subtype /Link /Rect [{:.2} {:.2} {:.2} {:.2}] /Border [0 0 0] /A << /S /URI /URI ({}) >> >>",
                link.left_pt,
                link.bottom_pt,
                link.right_pt,
                link.top_pt,
                escape_pdf_literal(&link.uri)
            )
            .as_bytes(),
        );
    }

    let xref_offset = pdf.len();
    pdf.extend_from_slice(format!("xref\n0 {}\n", offsets.len() + 1).as_bytes());
    pdf.extend_from_slice(b"0000000000 65535 f \n");
    for offset in &offsets {
        pdf.extend_from_slice(format!("{offset:010} 00000 n \n").as_bytes());
    }
    pdf.extend_from_slice(
        format!(
            "trailer\n<< /Size {} /Root 1 0 R >>\nstartxref\n{xref_offset}\n%%EOF\n",
            offsets.len() + 1
        )
        .as_bytes(),
    );
    Ok(pdf)
}

fn png_pages_to_pdf(png_pages: &[Vec<u8>]) -> Result<Vec<u8>, String> {
    if png_pages.is_empty() {
        return Err("multi-page PDF requires at least one rendered page".to_string());
    }
    struct Page {
        width: u32,
        height: u32,
        compressed_rgb: Vec<u8>,
    }
    let mut pages = Vec::with_capacity(png_pages.len());
    for png in png_pages {
        let image = image::load_from_memory(png)
            .map_err(|e| format!("Could not decode rendered PNG for PDF embedding: {e}"))?
            .to_rgb8();
        let (width, height) = image.dimensions();
        let mut encoder = ZlibEncoder::new(Vec::new(), Compression::default());
        encoder
            .write_all(image.as_raw())
            .map_err(|e| format!("Could not compress PDF page image: {e}"))?;
        pages.push(Page {
            width,
            height,
            compressed_rgb: encoder
                .finish()
                .map_err(|e| format!("Could not finish PDF page compression: {e}"))?,
        });
    }

    let mut pdf = Vec::new();
    pdf.extend_from_slice(b"%PDF-1.4\n");
    let mut offsets = Vec::new();
    push_pdf_object(
        &mut pdf,
        &mut offsets,
        b"<< /Type /Catalog /Pages 2 0 R >>".as_slice(),
    );
    let kids = (0..pages.len())
        .map(|index| format!("{} 0 R", 3 + index * 3))
        .collect::<Vec<_>>()
        .join(" ");
    push_pdf_object(
        &mut pdf,
        &mut offsets,
        format!("<< /Type /Pages /Kids [{kids}] /Count {} >>", pages.len()).as_bytes(),
    );
    for (index, page) in pages.iter().enumerate() {
        let page_object = 3 + index * 3;
        let image_object = page_object + 1;
        let content_object = page_object + 2;
        let width_pt = page.width as f32 * 72.0 / 96.0;
        let height_pt = page.height as f32 * 72.0 / 96.0;
        push_pdf_object(
            &mut pdf,
            &mut offsets,
            format!(
                "<< /Type /Page /Parent 2 0 R /MediaBox [0 0 {width_pt:.2} {height_pt:.2}] /Resources << /XObject << /Im0 {image_object} 0 R >> >> /Contents {content_object} 0 R >>"
            )
            .as_bytes(),
        );
        push_pdf_stream_object(
            &mut pdf,
            &mut offsets,
            format!(
                "<< /Type /XObject /Subtype /Image /Width {} /Height {} /ColorSpace /DeviceRGB /BitsPerComponent 8 /Filter /FlateDecode /Length {} >>",
                page.width,
                page.height,
                page.compressed_rgb.len()
            )
            .as_bytes(),
            &page.compressed_rgb,
        );
        let content = format!("q\n{width_pt:.2} 0 0 {height_pt:.2} 0 0 cm\n/Im0 Do\nQ\n");
        push_pdf_stream_object(
            &mut pdf,
            &mut offsets,
            format!("<< /Length {} >>", content.len()).as_bytes(),
            content.as_bytes(),
        );
    }
    let xref_offset = pdf.len();
    pdf.extend_from_slice(format!("xref\n0 {}\n", offsets.len() + 1).as_bytes());
    pdf.extend_from_slice(b"0000000000 65535 f \n");
    for offset in &offsets {
        pdf.extend_from_slice(format!("{offset:010} 00000 n \n").as_bytes());
    }
    pdf.extend_from_slice(
        format!(
            "trailer\n<< /Size {} /Root 1 0 R >>\nstartxref\n{xref_offset}\n%%EOF\n",
            offsets.len() + 1
        )
        .as_bytes(),
    );
    Ok(pdf)
}

fn push_pdf_object(pdf: &mut Vec<u8>, offsets: &mut Vec<usize>, body: &[u8]) {
    offsets.push(pdf.len());
    let object_number = offsets.len();
    pdf.extend_from_slice(format!("{object_number} 0 obj\n").as_bytes());
    pdf.extend_from_slice(body);
    pdf.extend_from_slice(b"\nendobj\n");
}

fn push_pdf_stream_object(
    pdf: &mut Vec<u8>,
    offsets: &mut Vec<usize>,
    dictionary: &[u8],
    stream: &[u8],
) {
    offsets.push(pdf.len());
    let object_number = offsets.len();
    pdf.extend_from_slice(format!("{object_number} 0 obj\n").as_bytes());
    pdf.extend_from_slice(dictionary);
    pdf.extend_from_slice(b"\nstream\n");
    pdf.extend_from_slice(stream);
    pdf.extend_from_slice(b"\nendstream\nendobj\n");
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use std::io::{Cursor, Read};
    use tempfile::tempdir;

    #[test]
    fn pdf_image_compression_is_lossless_deterministic_and_compact() {
        // Hand-crafted tall page: white background, coloured tracks and sparse
        // marks, independent of host fonts and external scientific inputs.
        let image = image::RgbImage::from_fn(512, 2048, |x, y| {
            if y % 100 < 2 {
                image::Rgb([20, 80, 160])
            } else if x % 97 < 3 && y % 31 < 12 {
                image::Rgb([180, 40, 60])
            } else {
                image::Rgb([255, 255, 255])
            }
        });
        let mut png = Cursor::new(Vec::new());
        image.write_to(&mut png, image::ImageFormat::Png).unwrap();
        let pdf = png_bytes_to_single_page_pdf(png.get_ref(), &[]).unwrap();
        assert_eq!(
            pdf,
            png_bytes_to_single_page_pdf(png.get_ref(), &[]).unwrap()
        );
        let stream_start = pdf.windows(7).position(|v| v == b"stream\n").unwrap() + 7;
        let header = std::str::from_utf8(&pdf[..stream_start]).unwrap();
        assert!(header.contains("/Filter /FlateDecode"));
        assert!(header.contains("/Width 512 /Height 2048"));
        assert!(header.contains("/MediaBox [0 0 384.00 1536.00]"));
        let length: usize = header
            .rsplit("/Length ")
            .next()
            .unwrap()
            .split_whitespace()
            .next()
            .unwrap()
            .parse()
            .unwrap();
        let mut restored = Vec::new();
        flate2::read::ZlibDecoder::new(&pdf[stream_start..stream_start + length])
            .read_to_end(&mut restored)
            .unwrap();
        assert_eq!(restored, *image.as_raw());
        assert!(
            pdf.len() < restored.len() / 10,
            "{} PDF bytes for {} RGB bytes",
            pdf.len(),
            restored.len()
        );
        assert!(pdf[stream_start + length..].starts_with(b"\nendstream\n"));

        let startxref = pdf.windows(10).rposition(|v| v == b"startxref\n").unwrap() + 10;
        let xref: usize = std::str::from_utf8(&pdf[startxref..])
            .unwrap()
            .lines()
            .next()
            .unwrap()
            .parse()
            .unwrap();
        let entries = std::str::from_utf8(&pdf[xref..]).unwrap();
        assert!(entries.starts_with("xref\n0 6\n"));
        for (i, entry) in entries.lines().skip(3).take(5).enumerate() {
            let offset: usize = entry.split_whitespace().next().unwrap().parse().unwrap();
            assert!(pdf[offset..].starts_with(format!("{} 0 obj\n", i + 1).as_bytes()));
        }
    }

    #[test]
    fn multi_page_pdf_preserves_order_dimensions_and_lossless_streams() {
        let mut pngs = Vec::new();
        for (width, height, colour) in [
            (20, 10, image::Rgb([12, 34, 56])),
            (13, 27, image::Rgb([210, 180, 30])),
        ] {
            let image = image::RgbImage::from_pixel(width, height, colour);
            let mut png = Cursor::new(Vec::new());
            image.write_to(&mut png, image::ImageFormat::Png).unwrap();
            pngs.push(png.into_inner());
        }
        let pdf = png_pages_to_pdf(&pngs).unwrap();
        let text = String::from_utf8_lossy(&pdf);
        assert!(text.contains("/Kids [3 0 R 6 0 R] /Count 2"));
        assert!(text.contains("/MediaBox [0 0 15.00 7.50]"));
        assert!(text.contains("/MediaBox [0 0 9.75 20.25]"));
        assert_eq!(text.matches("/Filter /FlateDecode").count(), 2);
        assert_eq!(text.matches("/Type /Page /Parent").count(), 2);
        assert_eq!(pdf, png_pages_to_pdf(&pngs).unwrap());
    }

    #[test]
    fn public_multi_page_renderer_keeps_each_svg_as_one_page() {
        let temp = tempdir().unwrap();
        let first = temp.path().join("context.svg");
        let second = temp.path().join("detail.svg");
        let output = temp.path().join("combined.pdf");
        fs::write(
            &first,
            r##"<svg xmlns="http://www.w3.org/2000/svg" width="40" height="20"><rect width="40" height="20" fill="#ffffff"/></svg>"##,
        )
        .unwrap();
        fs::write(
            &second,
            r##"<svg xmlns="http://www.w3.org/2000/svg" width="30" height="50"><rect width="30" height="50" fill="#0f766e"/></svg>"##,
        )
        .unwrap();
        let summary = render_svg_files_to_pdf(
            &[first.as_path(), second.as_path()],
            &output,
            SvgPngRenderOptions::default(),
        )
        .unwrap();
        assert_eq!(summary.page_count, 2);
        assert_eq!((summary.pages[0].width, summary.pages[0].height), (40, 20));
        assert_eq!((summary.pages[1].width, summary.pages[1].height), (30, 50));
        let pdf = fs::read(output).unwrap();
        assert!(String::from_utf8_lossy(&pdf).contains("/Count 2"));
    }

    #[test]
    fn multi_page_renderer_rejects_empty_input_without_writing() {
        let temp = tempdir().unwrap();
        let output = temp.path().join("must-not-exist.pdf");
        let error =
            render_svg_files_to_pdf(&[], &output, SvgPngRenderOptions::default()).unwrap_err();
        assert!(error.contains("at least one INPUT.svg"));
        assert!(!output.exists());
    }

    #[test]
    fn audited_pdf_keeps_legacy_bytes_summary_and_exact_embedded_font_audit() {
        use crate::svg_png::audit_test_support::{ALPHA, render};
        use std::cell::Cell;

        let temp = tempdir().unwrap();
        let input = temp.path().join("synthetic-text.svg");
        let output = temp.path().join("synthetic-text.pdf");
        let svg = format!(
            r#"<svg xmlns="http://www.w3.org/2000/svg" width="40" height="20"><text y="18" font-family="{ALPHA}">AAAA</text></svg>"#
        );
        fs::write(&input, &svg).unwrap();
        for scale in [1.0, 2.0] {
            let options = SvgPngRenderOptions {
                scale,
                drop_dotplot_metadata: false,
            };
            let (legacy, no_fonts) =
                render_svg_file_to_pdf_impl(&input, &output, options, &[], &[], |_, options| {
                    render(&svg, options, false)
                })
                .unwrap();
            assert!(no_fonts.is_empty());
            let legacy_bytes = fs::read(&output).unwrap();
            let calls = Cell::new(0);
            let (audited, fonts) =
                render_svg_file_to_pdf_impl(&input, &output, options, &[], &[], |_, options| {
                    calls.set(calls.get() + 1);
                    render(&svg, options, true)
                })
                .unwrap();
            assert_eq!(
                calls.get(),
                1,
                "PDF must use a single audited rasterization"
            );
            assert_eq!(legacy, audited);
            assert_eq!(
                serde_json::to_value(&legacy).unwrap(),
                serde_json::to_value(&audited).unwrap()
            );
            assert_eq!(legacy_bytes, fs::read(&output).unwrap());
            let (png, expected_fonts) = render(&svg, options, true).unwrap();
            assert_eq!(fonts, expected_fonts);
            assert_eq!(fonts.len(), 1);
            assert_eq!(
                legacy_bytes,
                png_bytes_to_single_page_pdf(&png.bytes, &[]).unwrap()
            );
        }
    }

    #[test]
    fn public_audited_pdf_preserves_legacy_output_for_an_svg_without_glyphs() {
        let temp = tempdir().unwrap();
        let input = temp.path().join("shapes.svg");
        let output = temp.path().join("shapes.pdf");
        fs::write(&input, r#"<svg xmlns="http://www.w3.org/2000/svg" width="20" height="20"><rect width="10" height="10"/></svg>"#).unwrap();
        let options = SvgPngRenderOptions::default();
        let legacy = render_svg_file_to_pdf(&input, &output, options).unwrap();
        let bytes = fs::read(&output).unwrap();
        let (audited, fonts) = render_svg_file_to_pdf_audited(&input, &output, options).unwrap();
        assert_eq!(legacy, audited);
        assert_eq!(bytes, fs::read(&output).unwrap());
        assert!(fonts.is_empty());
    }

    #[test]
    fn failed_font_audit_does_not_write_a_pdf() {
        let temp = tempdir().unwrap();
        let input = temp.path().join("input.svg");
        let output = temp.path().join("must-not-exist.pdf");
        let error = render_svg_file_to_pdf_impl(
            &input,
            &output,
            SvgPngRenderOptions::default(),
            &[],
            &[],
            |_, _| Err("Could not audit SVG used font: font bytes are unreadable".into()),
        )
        .unwrap_err();
        assert!(error.contains("unreadable"));
        assert!(!output.exists());
    }

    #[test]
    fn render_svg_file_to_pdf_embeds_one_image_page() {
        let temp = tempdir().expect("tempdir");
        let input = temp.path().join("demo.svg");
        let output = temp.path().join("demo.pdf");
        fs::write(
            &input,
            r##"<svg xmlns="http://www.w3.org/2000/svg" width="24" height="16"><rect width="24" height="16" fill="#ffffff"/><circle cx="12" cy="8" r="5" fill="#0f766e"/></svg>"##,
        )
        .expect("write svg");

        let summary = render_svg_file_to_pdf(&input, &output, SvgPngRenderOptions::default())
            .expect("render pdf");
        let pdf = fs::read(&output).expect("read pdf");
        assert_eq!(summary.width, 24);
        assert_eq!(summary.height, 16);
        assert_eq!(summary.uri_link_count, 0);
        assert!(pdf.starts_with(b"%PDF-1.4\n"));
        assert!(String::from_utf8_lossy(&pdf).contains("/Subtype /Image"));
    }

    #[test]
    fn render_svg_file_to_pdf_writes_valid_uri_annotations_in_page_coordinates() {
        let temp = tempdir().expect("tempdir");
        let input = temp.path().join("linked.svg");
        let output = temp.path().join("linked.pdf");
        fs::write(
            &input,
            r##"<svg xmlns="http://www.w3.org/2000/svg" width="100" height="100"><rect width="100" height="100" fill="#ffffff"/></svg>"##,
        )
        .expect("write svg");
        let summary = render_svg_file_to_pdf_with_links(
            &input,
            &output,
            SvgPngRenderOptions {
                scale: 2.0,
                drop_dotplot_metadata: false,
            },
            &[SvgPdfUriLink {
                x_px: 10.0,
                y_px: 20.0,
                width_px: 30.0,
                height_px: 10.0,
                uri: "https://regulation.ensembl.org/regulatory_features/homo_sapiens/ENSR1_958"
                    .to_string(),
            }],
            &["regulation.ensembl.org"],
        )
        .expect("render linked pdf");
        let pdf = String::from_utf8_lossy(&fs::read(&output).expect("read pdf")).into_owned();
        assert_eq!(summary.uri_link_count, 1);
        assert!(pdf.contains("/Annots [6 0 R]"));
        assert!(pdf.contains("/Subtype /Link"));
        assert!(pdf.contains("/Rect [15.00 105.00 60.00 120.00]"));
        assert!(pdf.contains(
            "/URI (https://regulation.ensembl.org/regulatory_features/homo_sapiens/ENSR1_958)"
        ));
    }

    #[test]
    fn render_svg_file_to_pdf_rejects_unapproved_or_non_https_links() {
        let temp = tempdir().expect("tempdir");
        let input = temp.path().join("unsafe.svg");
        let output = temp.path().join("unsafe.pdf");
        fs::write(
            &input,
            r##"<svg xmlns="http://www.w3.org/2000/svg" width="20" height="20"/>"##,
        )
        .expect("write svg");
        for uri in [
            "http://regulation.ensembl.org/feature/ENSR1",
            "https://evil.regulation.ensembl.org/feature/ENSR1",
        ] {
            let error = render_svg_file_to_pdf_with_links(
                &input,
                &output,
                SvgPngRenderOptions::default(),
                &[SvgPdfUriLink {
                    x_px: 1.0,
                    y_px: 1.0,
                    width_px: 5.0,
                    height_px: 5.0,
                    uri: uri.to_string(),
                }],
                &["regulation.ensembl.org"],
            )
            .expect_err("unsafe URL must fail");
            assert!(error.contains("HTTPS") || error.contains("not explicitly allowed"));
        }
    }
}
