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
}

/// Renders one SVG file into a single-page PDF.
pub fn render_svg_file_to_pdf(
    input_path: &Path,
    output_path: &Path,
    options: SvgPngRenderOptions,
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

    let rendered_bytes = render_svg_file_to_png_bytes(input_path, options)?;
    let pdf = png_bytes_to_single_page_pdf(&rendered_bytes.bytes)?;
    std::fs::write(output_path, pdf)
        .map_err(|e| format!("Could not write PDF '{}': {e}", output_path.display()))?;

    let page_width_pt = rendered_bytes.width as f32 * 72.0 / 96.0;
    let page_height_pt = rendered_bytes.height as f32 * 72.0 / 96.0;
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
    })
}

fn png_bytes_to_single_page_pdf(png_bytes: &[u8]) -> Result<Vec<u8>, String> {
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
    push_pdf_object(
        &mut pdf,
        &mut offsets,
        format!(
            "<< /Type /Page /Parent 2 0 R /MediaBox [0 0 {page_width_pt:.2} {page_height_pt:.2}] /Resources << /XObject << /Im0 4 0 R >> >> /Contents 5 0 R >>"
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

    let xref_offset = pdf.len();
    pdf.extend_from_slice(format!("xref\n0 {}\n", offsets.len() + 1).as_bytes());
    pdf.extend_from_slice(b"0000000000 65535 f \n");
    for offset in offsets {
        pdf.extend_from_slice(format!("{offset:010} 00000 n \n").as_bytes());
    }
    pdf.extend_from_slice(
        format!("trailer\n<< /Size 6 /Root 1 0 R >>\nstartxref\n{xref_offset}\n%%EOF\n").as_bytes(),
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
        assert!(pdf.starts_with(b"%PDF-1.4\n"));
        assert!(String::from_utf8_lossy(&pdf).contains("/Subtype /Image"));
    }
}
