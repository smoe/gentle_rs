//! Deterministic SVG-to-PDF conversion helpers shared by headless exports.
//!
//! The PDF path intentionally reuses GENtle's existing `resvg` rasterization
//! helper, then embeds the resulting RGB image into one simple PDF page. That
//! keeps PDF generation dependency-light and visually consistent with SVG/PNG
//! exports.

use serde::Serialize;
use std::path::Path;

use crate::svg_png::{SvgPngRenderOptions, render_svg_file_to_png_bytes};

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

/// Renders one SVG file into a single-page raster-backed PDF with validated
/// URI annotations. `allowed_hosts` uses exact host matching.
pub fn render_svg_file_to_pdf_with_links(
    input_path: &Path,
    output_path: &Path,
    options: SvgPngRenderOptions,
    links: &[SvgPdfUriLink],
    allowed_hosts: &[&str],
) -> Result<SvgPdfRenderSummary, String> {
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
    let rendered_bytes = render_svg_file_to_png_bytes(input_path, options)?;
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

    Ok(SvgPdfRenderSummary {
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
    })
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
    let rgb = image.into_raw();
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
            "<< /Type /XObject /Subtype /Image /Width {width} /Height {height} /ColorSpace /DeviceRGB /BitsPerComponent 8 /Length {} >>",
            rgb.len()
        )
        .as_bytes(),
        &rgb,
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
    use tempfile::tempdir;

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
