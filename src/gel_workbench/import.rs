//! Bounded, read-only source ingestion. PDF rasters and opt-in native thumbnails
//! remain display evidence, never reconstructed numerical electropherograms.

use std::{
    collections::BTreeSet,
    io::{self, Cursor, Read, Write},
    path::Path,
};

use base64::{Engine as _, engine::general_purpose::STANDARD};
use flate2::bufread::ZlibDecoder;
use gentle_protocol::{
    EngineError, ErrorCode,
    gel_image::GelImageRecord,
    gel_workbench::{GelSourceImportRequest, GelSourceRecord},
};
use image::{ColorType, ImageDecoder, ImageEncoder, ImageFormat, ImageReader, Limits};
use lopdf::{Dictionary, Document, LoadOptions, Object, ObjectId, Stream};

use crate::{
    digest_utils::sha256_hex_bytes,
    gel_image::{
        MAX_GEL_IMAGE_BYTES, MAX_GEL_IMAGE_PIXELS, import_gel_image_bytes, validate_label,
    },
};

const MAX_PDF_OBJECTS: usize = 4096;
const MAX_PDF_PAGES: usize = 32;
const MAX_PDF_CANDIDATES: usize = 128;
const MAX_PDF_TABLE_BYTES: usize = 1024 * 1024;
const MAX_NATIVE_MEMBERS: usize = 128;
const MAX_NATIVE_MEMBER_BYTES: u64 = 64 * 1024 * 1024;
const MAX_RASTER_ALLOC: u64 = 256 * 1024 * 1024;

fn unsupported(message: impl Into<String>) -> EngineError {
    EngineError::new(ErrorCode::Unsupported, message)
}

fn invalid_pdf() -> EngineError {
    EngineError::invalid_input("Malformed or out-of-bounds PDF source")
}

/// Retain the exact supplied source bytes and import a separately hash-bound
/// raster through the shared gel-image boundary. No files, PDF actions, external
/// resources, native XML, passwords, logs, or sample metadata are accessed.
pub fn import_source(
    request: &GelSourceImportRequest,
    bytes: &[u8],
) -> Result<(GelSourceRecord, GelImageRecord), EngineError> {
    validate_label(&request.image_id, "image_id")?;
    let source_name = Path::new(&request.path)
        .file_name()
        .and_then(|name| name.to_str())
        .ok_or_else(|| EngineError::invalid_input("Source filename must be valid UTF-8"))?;
    validate_label(source_name, "source_name")?;
    if bytes.is_empty() || bytes.len() > MAX_GEL_IMAGE_BYTES {
        return Err(EngineError::invalid_input(
            "Gel source must be nonempty and at most 16 MiB",
        ));
    }
    let native = Path::new(source_name)
        .extension()
        .and_then(|ext| ext.to_str())
        .is_some_and(|ext| ext.eq_ignore_ascii_case("d5000"));
    let (raster, image_name, extraction, warnings) = if bytes.starts_with(b"%PDF-") {
        if request.allow_native_thumbnail {
            return Err(EngineError::invalid_input(
                "allow_native_thumbnail applies only to D5000 archives",
            ));
        }
        let (raster, extraction, warnings) = pdf_raster(bytes, request.pdf_object)?;
        (raster, "embedded-raster.png", extraction, warnings)
    } else if native {
        if request.pdf_object.is_some() {
            return Err(EngineError::invalid_input(
                "pdf_object applies only to PDF sources",
            ));
        }
        if !request.allow_native_thumbnail {
            return Err(unsupported(
                "D5000 numerical data, including encrypted data, are unsupported; set allow_native_thumbnail=true for display-only Thumbnail.png import",
            ));
        }
        let (raster, warnings) = native_thumbnail(bytes)?;
        (
            raster,
            "Thumbnail.png",
            "d5000_thumbnail_v1".into(),
            warnings,
        )
    } else {
        if request.pdf_object.is_some() || request.allow_native_thumbnail {
            return Err(EngineError::invalid_input(
                "Source selection options do not apply to ordinary raster images",
            ));
        }
        let format = image::guess_format(bytes).ok();
        if !matches!(format, Some(ImageFormat::Png | ImageFormat::Jpeg)) {
            return Err(unsupported(
                "Unsupported gel source; expected PDF, PNG, JPEG or an explicitly opted-in D5000 thumbnail (TIFF requires gel-image import with tiff_page=0)",
            ));
        }
        check_raster_header(bytes, format.expect("checked above"))?;
        (
            bytes.to_vec(),
            source_name,
            "original_raster_v1".into(),
            vec![],
        )
    };
    let mut image = import_gel_image_bytes(&request.image_id, image_name, &raster, None)?;
    image.descriptor.warnings.extend(warnings.iter().cloned());
    let source = GelSourceRecord {
        request: GelSourceImportRequest {
            path: source_name.into(),
            ..request.clone()
        },
        image_id: request.image_id.clone(),
        source_name: source_name.into(),
        source_sha256: sha256_hex_bytes(bytes),
        original_base64: STANDARD.encode(bytes),
        extracted_image_sha256: image.descriptor.sha256.clone(),
        extraction,
        warnings,
    };
    Ok((source, image))
}

fn bounded_reader(bytes: &[u8], format: ImageFormat) -> ImageReader<Cursor<&[u8]>> {
    let mut reader = ImageReader::with_format(Cursor::new(bytes), format);
    let mut limits = Limits::default();
    limits.max_image_width = Some(32768);
    limits.max_image_height = Some(32768);
    limits.max_alloc = Some(MAX_RASTER_ALLOC);
    reader.limits(limits);
    reader
}

fn check_dimensions(width: u32, height: u32) -> Result<usize, EngineError> {
    let pixels = u64::from(width) * u64::from(height);
    if width < 2 || height < 2 || width > 32768 || height > 32768 || pixels > MAX_GEL_IMAGE_PIXELS {
        return Err(EngineError::invalid_input(
            "Source raster must be at least 2x2, at most 32768 per side and at most 32 megapixels",
        ));
    }
    Ok(pixels as usize)
}

fn check_raster_header(bytes: &[u8], format: ImageFormat) -> Result<(), EngineError> {
    let decoder = bounded_reader(bytes, format)
        .into_decoder()
        .map_err(|_| EngineError::invalid_input("Invalid or oversized source raster"))?;
    let (width, height) = decoder.dimensions();
    check_dimensions(width, height)?;
    if decoder.total_bytes() > MAX_RASTER_ALLOC {
        return Err(EngineError::invalid_input(
            "Source raster exceeds allocation limit",
        ));
    }
    Ok(())
}

fn whitespace(byte: u8) -> bool {
    matches!(byte, 0 | 9 | 10 | 12 | 13 | 32)
}

fn delimiter(byte: u8) -> bool {
    whitespace(byte) || b"()<>[]{}/%".contains(&byte)
}

fn word<'a>(bytes: &mut &'a [u8]) -> Result<&'a [u8], EngineError> {
    *bytes = bytes.trim_ascii_start();
    let end = bytes
        .iter()
        .position(|byte| delimiter(*byte))
        .unwrap_or(bytes.len());
    if end == 0 {
        return Err(invalid_pdf());
    }
    let value = &bytes[..end];
    *bytes = &bytes[end..];
    Ok(value)
}

fn number(bytes: &mut &[u8]) -> Result<usize, EngineError> {
    let value = word(bytes)?;
    if value.len() > 10 || !value.iter().all(u8::is_ascii_digit) {
        return Err(invalid_pdf());
    }
    std::str::from_utf8(value)
        .ok()
        .and_then(|s| s.parse().ok())
        .ok_or_else(invalid_pdf)
}

fn hex(byte: u8) -> Option<u8> {
    match byte {
        b'0'..=b'9' => Some(byte - b'0'),
        b'a'..=b'f' => Some(byte - b'a' + 10),
        b'A'..=b'F' => Some(byte - b'A' + 10),
        _ => None,
    }
}

// A deliberately conservative lexical guard, including #xx-escaped names. It
// may refuse unusual literal/binary occurrences, but never lets an encrypted or
// compressed object table reach lopdf's eager loading/decryption paths.
fn guard_pdf_names(bytes: &[u8]) -> Result<(), EngineError> {
    let mut rest = bytes;
    let mut stream_bytes = 0usize;
    while let Some(slash) = rest.iter().position(|byte| *byte == b'/') {
        rest = &rest[slash + 1..];
        let end = rest
            .iter()
            .position(|byte| delimiter(*byte))
            .unwrap_or(rest.len());
        let token = &rest[..end];
        rest = &rest[end..];
        let mut name = Vec::with_capacity(32);
        let mut index = 0;
        while index < token.len() {
            let byte = if token[index] == b'#' {
                let a = token.get(index + 1).and_then(|byte| hex(*byte));
                let b = token.get(index + 2).and_then(|byte| hex(*byte));
                match (a, b) {
                    (Some(a), Some(b)) => {
                        index += 2;
                        a * 16 + b
                    }
                    _ => break,
                }
            } else {
                token[index]
            };
            if name.len() == 128 {
                break;
            }
            name.push(byte);
            index += 1;
        }
        match name.as_slice() {
            b"Encrypt" => {
                return Err(unsupported(
                    "Encrypted PDF sources are unsupported; no password is attempted",
                ));
            }
            b"ObjStm" | b"XRef" | b"XRefStm" | b"Prev" => {
                return Err(unsupported(
                    "Only single-revision PDF sources with classic cross-reference tables are supported",
                ));
            }
            b"Length" if token != b"Length" => {
                return Err(unsupported(
                    "Escaped PDF stream-length keys are outside the bounded import subset",
                ));
            }
            b"Length" => {
                let mut value = rest;
                let length = number(&mut value)?;
                stream_bytes = stream_bytes.checked_add(length).ok_or_else(invalid_pdf)?;
                if stream_bytes > bytes.len()
                    || value
                        .trim_ascii_start()
                        .first()
                        .is_some_and(u8::is_ascii_digit)
                {
                    return Err(unsupported(
                        "PDF stream lengths must be direct and their total cannot exceed the source byte budget",
                    ));
                }
            }
            _ => {}
        }
    }
    Ok(())
}

fn pdf_preflight(bytes: &[u8]) -> Result<(), EngineError> {
    guard_pdf_names(bytes)?;
    let tail = &bytes[bytes.len().saturating_sub(512)..];
    let start = tail
        .windows(9)
        .rposition(|part| part == b"startxref")
        .ok_or_else(invalid_pdf)?;
    let mut pointer = &tail[start + 9..];
    let offset = number(&mut pointer)?;
    if pointer.trim_ascii() != b"%%EOF" {
        return Err(invalid_pdf());
    }
    let mut table = bytes.get(offset..).ok_or_else(invalid_pdf)?;
    if word(&mut table)? != b"xref" {
        return Err(unsupported("PDF cross-reference streams are unsupported"));
    }
    let mut entries = Vec::new();
    let mut ids = BTreeSet::new();
    let mut count = 0usize;
    loop {
        table = table.trim_ascii_start();
        if table.starts_with(b"trailer") {
            break;
        }
        let first = number(&mut table)?;
        let length = number(&mut table)?;
        count = count.checked_add(length).ok_or_else(invalid_pdf)?;
        if length == 0 || count > MAX_PDF_OBJECTS || first.saturating_add(length) > MAX_PDF_OBJECTS
        {
            return Err(EngineError::invalid_input(
                "PDF exceeds the 4096-object limit",
            ));
        }
        table = table.strip_prefix(b" ").unwrap_or(table);
        table = table
            .strip_prefix(b"\r\n")
            .or_else(|| table.strip_prefix(b"\n"))
            .or_else(|| table.strip_prefix(b"\r"))
            .ok_or_else(invalid_pdf)?;
        for id in first..first + length {
            // Match the library's strict, fixed-width xref grammar rather
            // than admitting malformed rows that could trigger reconstruction.
            let row = table.get(..20).ok_or_else(invalid_pdf)?;
            if row[10] != b' '
                || row[16] != b' '
                || !matches!(&row[18..], b" \n" | b" \r" | b"\r\n")
            {
                return Err(invalid_pdf());
            }
            let object_offset = number(&mut &row[..10])?;
            let generation = number(&mut &row[11..16])?;
            let state = &row[17..18];
            table = &table[20..];
            if !ids.insert(id) || generation > u16::MAX as usize {
                return Err(invalid_pdf());
            }
            match state {
                b"n" if id > 0 && object_offset < offset => {
                    let mut object = &bytes[object_offset..offset];
                    if number(&mut object)? != id
                        || number(&mut object)? != generation
                        || word(&mut object)? != b"obj"
                    {
                        return Err(invalid_pdf());
                    }
                    entries.push(object_offset);
                }
                b"f" => {}
                _ => return Err(invalid_pdf()),
            }
        }
        if bytes.len() - table.len() - offset > MAX_PDF_TABLE_BYTES {
            return Err(invalid_pdf());
        }
    }
    entries.sort_unstable();
    if entries.is_empty() || entries.windows(2).any(|pair| pair[0] == pair[1]) {
        return Err(invalid_pdf());
    }
    // Bound each declared stream allocation by its physical object's span,
    // before the library can clone any stream. Indirect lengths are outside
    // this initial subset. No stream is decompressed here.
    entries.push(offset);
    for span in entries.windows(2) {
        let object = &bytes[span[0]..span[1]];
        let mut rest = object;
        while let Some(at) = rest.windows(7).position(|part| part == b"/Length") {
            rest = &rest[at + 7..];
            if rest.first().is_some_and(|byte| !delimiter(*byte)) {
                continue;
            }
            let length = number(&mut rest)?;
            if length > object.len()
                || rest
                    .trim_ascii_start()
                    .first()
                    .is_some_and(u8::is_ascii_digit)
            {
                return Err(unsupported(
                    "PDF stream lengths must be direct and bounded by their object",
                ));
            }
        }
    }
    Ok(())
}

fn resolve<'a>(document: &'a Document, mut object: &'a Object) -> Result<&'a Object, EngineError> {
    for _ in 0..16 {
        match object {
            Object::Reference(id) => object = document.objects.get(id).ok_or_else(invalid_pdf)?,
            _ => return Ok(object),
        }
    }
    Err(invalid_pdf())
}

fn field<'a>(
    document: &'a Document,
    dict: &'a Dictionary,
    key: &[u8],
) -> Result<&'a Object, EngineError> {
    resolve(document, dict.get(key).map_err(|_| invalid_pdf())?)
}

fn pdf_page_limit(document: &Document) -> Result<(), EngineError> {
    let root = field(document, &document.trailer, b"Root")?
        .as_dict()
        .map_err(|_| invalid_pdf())?;
    let mut pending = vec![root.get(b"Pages").map_err(|_| invalid_pdf())?];
    let mut seen = BTreeSet::new();
    let mut pages = 0;
    while let Some(node) = pending.pop() {
        let id = node.as_reference().map_err(|_| invalid_pdf())?;
        if !seen.insert(id) || seen.len() > MAX_PDF_OBJECTS {
            return Err(invalid_pdf());
        }
        let dict = resolve(document, node)?
            .as_dict()
            .map_err(|_| invalid_pdf())?;
        match dict
            .get(b"Type")
            .and_then(Object::as_name)
            .map_err(|_| invalid_pdf())?
        {
            b"Page" => pages += 1,
            b"Pages" => {
                let count = field(document, dict, b"Count")?
                    .as_i64()
                    .map_err(|_| invalid_pdf())?;
                if !(0..=MAX_PDF_PAGES as i64).contains(&count) {
                    return Err(EngineError::invalid_input("PDF exceeds the 32-page limit"));
                }
                let kids = field(document, dict, b"Kids")?
                    .as_array()
                    .map_err(|_| invalid_pdf())?;
                if kids.len() > MAX_PDF_PAGES || pending.len() + kids.len() > MAX_PDF_OBJECTS {
                    return Err(invalid_pdf());
                }
                pending.extend(kids);
            }
            _ => return Err(invalid_pdf()),
        }
        if pages > MAX_PDF_PAGES {
            return Err(EngineError::invalid_input("PDF exceeds the 32-page limit"));
        }
    }
    if pages == 0 {
        return Err(invalid_pdf());
    }
    Ok(())
}

#[derive(Clone, Copy)]
struct RasterSpec {
    width: u32,
    height: u32,
    channels: usize,
    jpeg: bool,
    soft_mask: Option<ObjectId>,
}

fn raster_spec(
    document: &Document,
    stream: &Stream,
    mask: bool,
) -> Result<RasterSpec, EngineError> {
    let dict = &stream.dict;
    if dict.get(b"Subtype").and_then(Object::as_name).ok() != Some(b"Image") {
        return Err(unsupported(
            "Selected PDF object is not an embedded raster image",
        ));
    }
    for key in [
        b"Decode".as_slice(),
        b"Mask",
        b"Matte",
        b"Alternates",
        b"OPI",
        b"F",
        b"FFilter",
        b"FDecodeParms",
    ] {
        if dict.has(key) {
            return Err(unsupported(
                "PDF decode transforms, color-key masks, matte and external/alternate images are unsupported",
            ));
        }
    }
    if dict
        .get(b"ImageMask")
        .is_ok_and(|value| !matches!(value, Object::Boolean(false)))
        || dict
            .get(b"SMaskInData")
            .is_ok_and(|value| !matches!(value, Object::Integer(0)))
    {
        return Err(unsupported("PDF stencil and in-data masks are unsupported"));
    }
    if field(document, dict, b"BitsPerComponent")?.as_i64().ok() != Some(8) {
        return Err(unsupported("PDF rasters require 8 bits per component"));
    }
    let channels = match field(document, dict, b"ColorSpace")?.as_name().ok() {
        Some(b"DeviceRGB") if !mask => 3,
        Some(b"DeviceGray") => 1,
        _ => {
            return Err(unsupported(
                "PDF rasters require DeviceRGB or DeviceGray (soft masks require DeviceGray)",
            ));
        }
    };
    let width = field(document, dict, b"Width")?
        .as_i64()
        .ok()
        .and_then(|value| u32::try_from(value).ok())
        .ok_or_else(invalid_pdf)?;
    let height = field(document, dict, b"Height")?
        .as_i64()
        .ok()
        .and_then(|value| u32::try_from(value).ok())
        .ok_or_else(invalid_pdf)?;
    check_dimensions(width, height)?;
    let filter = field(document, dict, b"Filter")?;
    let filter = match filter {
        Object::Array(filters) if filters.len() == 1 => resolve(document, &filters[0])?,
        value => value,
    };
    let jpeg = match filter.as_name().ok() {
        Some(b"FlateDecode") => false,
        Some(b"DCTDecode") if !mask => true,
        _ => {
            return Err(unsupported(
                "PDF images require one FlateDecode or DCTDecode filter; soft masks require FlateDecode",
            ));
        }
    };
    if let Ok(parameters) = dict.get(b"DecodeParms") {
        match resolve(document, parameters)? {
            Object::Null => {}
            Object::Array(values)
                if values.len() == 1 && matches!(resolve(document, &values[0])?, Object::Null) => {}
            _ => {
                return Err(unsupported(
                    "PDF predictor and JPEG decode parameters are unsupported",
                ));
            }
        }
    }
    let soft_mask = match dict.get(b"SMask") {
        Err(_) | Ok(Object::Null) => None,
        Ok(Object::Name(name)) if name == b"None" => None,
        Ok(Object::Reference(id)) if !mask => Some(*id),
        _ => return Err(unsupported("Unsupported or nested PDF soft mask")),
    };
    Ok(RasterSpec {
        width,
        height,
        channels,
        jpeg,
        soft_mask,
    })
}

fn stream_at(document: &Document, id: ObjectId) -> Result<&Stream, EngineError> {
    document
        .objects
        .get(&id)
        .and_then(|object| object.as_stream().ok())
        .ok_or_else(invalid_pdf)
}

fn selected_spec(document: &Document, stream: &Stream) -> Result<RasterSpec, EngineError> {
    let spec = raster_spec(document, stream, false)?;
    if let Some(id) = spec.soft_mask {
        let mask = raster_spec(document, stream_at(document, id)?, true)?;
        if mask.width != spec.width || mask.height != spec.height {
            return Err(unsupported(
                "PDF soft-mask dimensions must match the raster",
            ));
        }
    }
    Ok(spec)
}

fn decode_pixels(stream: &Stream, spec: RasterSpec) -> Result<Vec<u8>, EngineError> {
    let expected = check_dimensions(spec.width, spec.height)? * spec.channels;
    let mut pixels = vec![0; expected];
    if spec.jpeg {
        let decoder = bounded_reader(&stream.content, ImageFormat::Jpeg)
            .into_decoder()
            .map_err(|_| EngineError::invalid_input("Invalid bounded PDF JPEG"))?;
        let color = if spec.channels == 3 {
            ColorType::Rgb8
        } else {
            ColorType::L8
        };
        if decoder.dimensions() != (spec.width, spec.height)
            || decoder.color_type() != color
            || decoder.total_bytes() != expected as u64
        {
            return Err(unsupported(
                "PDF JPEG dimensions/components disagree with its image dictionary",
            ));
        }
        decoder
            .read_image(&mut pixels)
            .map_err(|_| EngineError::invalid_input("Cannot decode bounded PDF JPEG"))?;
    } else {
        let mut decoder = ZlibDecoder::new(Cursor::new(&stream.content));
        decoder
            .read_exact(&mut pixels)
            .map_err(|_| EngineError::invalid_input("Truncated or invalid PDF Flate raster"))?;
        let mut extra = [0];
        if decoder
            .read(&mut extra)
            .map_err(|_| EngineError::invalid_input("Invalid PDF Flate checksum"))?
            != 0
            || decoder.total_in() != stream.content.len() as u64
        {
            return Err(EngineError::invalid_input(
                "PDF Flate raster length must match its dimensions exactly",
            ));
        }
    }
    Ok(pixels)
}

struct RasterOutput(Vec<u8>);

impl Write for RasterOutput {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        if bytes.len() > MAX_GEL_IMAGE_BYTES.saturating_sub(self.0.len()) {
            return Err(io::Error::other("Derived raster exceeds 16 MiB"));
        }
        self.0.extend_from_slice(bytes);
        Ok(bytes.len())
    }

    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}

fn pdf_raster(
    bytes: &[u8],
    selection: Option<u32>,
) -> Result<(Vec<u8>, String, Vec<String>), EngineError> {
    pdf_preflight(bytes)?;
    let document = Document::load_mem_with_options(
        bytes,
        LoadOptions {
            strict: true,
            max_decompressed_size: Some(MAX_PDF_TABLE_BYTES),
            ..Default::default()
        },
    )
    .map_err(|_| invalid_pdf())?;
    if document.trailer.has(b"Encrypt") || document.encryption_state.is_some() {
        return Err(unsupported("Encrypted PDF sources are unsupported"));
    }
    if document.objects.len() > MAX_PDF_OBJECTS {
        return Err(invalid_pdf());
    }
    pdf_page_limit(&document)?;
    let mut masks = BTreeSet::new();
    let mut images = Vec::new();
    for (id, object) in &document.objects {
        if let Ok(stream) = object.as_stream() {
            if stream.dict.get(b"Subtype").and_then(Object::as_name).ok() == Some(b"Image") {
                images.push(*id);
                if images.len() > MAX_PDF_CANDIDATES {
                    return Err(EngineError::invalid_input(
                        "PDF exceeds the 128-image candidate limit",
                    ));
                }
                if let Ok(Object::Reference(mask)) = stream.dict.get(b"SMask") {
                    masks.insert(*mask);
                }
            }
        }
    }
    let (id, spec) = if let Some(number) = selection {
        let id = images
            .iter()
            .find(|id| id.0 == number)
            .copied()
            .ok_or_else(|| {
                EngineError::invalid_input("Selected PDF image object does not exist")
            })?;
        if masks.contains(&id) {
            return Err(unsupported("A PDF soft mask is not an original gel raster"));
        }
        (id, selected_spec(&document, stream_at(&document, id)?)?)
    } else {
        let supported: Vec<_> = images
            .iter()
            .filter(|id| !masks.contains(*id))
            .filter_map(|id| {
                selected_spec(&document, stream_at(&document, *id).ok()?)
                    .ok()
                    .map(|spec| (*id, spec))
            })
            .collect();
        match supported.as_slice() {
            [candidate] => *candidate,
            [] => {
                return Err(unsupported(
                    "PDF contains no supported embedded original raster",
                ));
            }
            _ => {
                return Err(EngineError::invalid_input(
                    "Ambiguous PDF: select the original raster with pdf_object; no largest-image or page guess is made",
                ));
            }
        }
    };
    let mut pixels = decode_pixels(stream_at(&document, id)?, spec)?;
    let mut channels = spec.channels;
    if let Some(mask_id) = spec.soft_mask {
        let mask_stream = stream_at(&document, mask_id)?;
        let alpha = decode_pixels(mask_stream, raster_spec(&document, mask_stream, true)?)?;
        let output_len = alpha.len() * (channels + 1);
        if pixels.len() as u64 + alpha.len() as u64 + output_len as u64 > MAX_RASTER_ALLOC {
            return Err(EngineError::invalid_input(
                "PDF raster and soft mask exceed allocation limit",
            ));
        }
        let mut combined = Vec::with_capacity(output_len);
        for (color, alpha) in pixels.chunks_exact(channels).zip(alpha) {
            combined.extend_from_slice(color);
            combined.push(alpha);
        }
        pixels = combined;
        channels += 1;
    }
    let color = match channels {
        1 => ColorType::L8,
        2 => ColorType::La8,
        3 => ColorType::Rgb8,
        4 => ColorType::Rgba8,
        _ => unreachable!(),
    };
    let mut output = RasterOutput(Vec::new());
    image::codecs::png::PngEncoder::new(&mut output)
        .write_image(&pixels, spec.width, spec.height, color.into())
        .map_err(|_| {
            EngineError::invalid_input("Cannot encode derived raster within the 16 MiB limit")
        })?;
    let extraction = format!(
        "pdf_embedded_raster_v1:object={}:generation={}:filter={}:soft_mask={}",
        id.0,
        id.1,
        if spec.jpeg {
            "DCTDecode"
        } else {
            "FlateDecode"
        },
        spec.soft_mask
            .map(|id| format!("{}:{}", id.0, id.1))
            .unwrap_or_else(|| "none".into())
    );
    let mut warnings = vec![
        "Extracted embedded PDF raster, not a rendered page or screenshot. Page placement, text, annotations and other objects are not measurements.".into(),
        "Vendor-processed display pixels are not numerical detector traces, fluorescence, concentration or sample metadata; no curves or metadata were inferred.".into(),
    ];
    if spec.soft_mask.is_some() {
        warnings.push("The raster's supported 8-bit grayscale soft mask is preserved as PNG alpha; it was not discarded or composited against a guessed background.".into());
    }
    if spec.jpeg {
        warnings.push(
            "Embedded JPEG is lossy; decoding to PNG does not restore acquisition information."
                .into(),
        );
    }
    Ok((output.0, extraction, warnings))
}

fn le16(bytes: &[u8], at: usize) -> Result<usize, EngineError> {
    bytes
        .get(at..at + 2)
        .map(|part| u16::from_le_bytes([part[0], part[1]]) as usize)
        .ok_or_else(|| EngineError::invalid_input("Truncated native ZIP directory"))
}

fn le32(bytes: &[u8], at: usize) -> Result<usize, EngineError> {
    bytes
        .get(at..at + 4)
        .map(|part| u32::from_le_bytes([part[0], part[1], part[2], part[3]]) as usize)
        .ok_or_else(|| EngineError::invalid_input("Truncated native ZIP directory"))
}

// Preflight before ZipArchive allocates its directory or deduplicates names.
// Only a classic, single-disk archive is admitted. All member payloads except
// the exact root Thumbnail.png remain unopened, including unencrypted logs.
struct NativeThumbnail {
    index: usize,
    count: usize,
    local: usize,
    central: usize,
    data_start: usize,
    size: u64,
    compressed: usize,
    encrypted_members: bool,
}

fn native_directory(bytes: &[u8]) -> Result<NativeThumbnail, EngineError> {
    let bad = || EngineError::invalid_input("Malformed or out-of-bounds native ZIP directory");
    let tail_start = bytes.len().saturating_sub(65535 + 22);
    let end = bytes[tail_start..]
        .windows(4)
        .rposition(|part| part == b"PK\x05\x06")
        .map(|at| at + tail_start)
        .ok_or_else(bad)?;
    let end_record = &bytes[end..];
    let count = le16(end_record, 10)?;
    let directory_size = le32(end_record, 12)?;
    let directory_start = le32(end_record, 16)?;
    if le16(end_record, 4)? != 0
        || le16(end_record, 6)? != 0
        || le16(end_record, 8)? != count
        || count == 0
        || count > MAX_NATIVE_MEMBERS
        || directory_size > MAX_PDF_TABLE_BYTES
        || directory_start.checked_add(directory_size) != Some(end)
        || end_record.len() != 22 + le16(end_record, 20)?
    {
        return Err(bad());
    }
    let mut directory = bytes.get(directory_start..end).ok_or_else(bad)?;
    let mut names = BTreeSet::new();
    let mut thumbnail = None;
    let mut encrypted = false;
    let mut declared_total = 0u64;
    for index in 0..count {
        if !directory.starts_with(b"PK\x01\x02") {
            return Err(bad());
        }
        let flags = le16(directory, 8)?;
        let compressed = le32(directory, 20)?;
        let size = le32(directory, 24)? as u64;
        let name_len = le16(directory, 28)?;
        let extra_len = le16(directory, 30)?;
        let comment_len = le16(directory, 32)?;
        let local = le32(directory, 42)?;
        declared_total += size;
        if le16(directory, 34)? != 0
            || name_len == 0
            || name_len > 240
            || size > MAX_NATIVE_MEMBER_BYTES
            || declared_total > MAX_RASTER_ALLOC
            || compressed > bytes.len()
            || local >= directory_start
        {
            return Err(bad());
        }
        let name = directory.get(46..46 + name_len).ok_or_else(bad)?;
        if !names.insert(name.to_vec()) {
            return Err(EngineError::invalid_input(
                "Ambiguous duplicate native ZIP member names",
            ));
        }
        encrypted |= flags & (1 | 64) != 0;
        if name == b"Thumbnail.png" {
            if flags & (1 | 64) != 0 {
                return Err(unsupported(
                    "Encrypted native thumbnail is unsupported; no password is attempted",
                ));
            }
            if size == 0 || size > MAX_GEL_IMAGE_BYTES as u64 {
                return Err(bad());
            }
            // Check the local header as well; never let contradictory flags
            // route an encrypted member into the ordinary decompressor.
            let header = bytes.get(local..directory_start).ok_or_else(bad)?;
            if !header.starts_with(b"PK\x03\x04")
                || le16(header, 6)? != flags
                || le16(header, 8)? != le16(directory, 10)?
                || le16(header, 26)? != name_len
                || header.get(30..30 + name_len) != Some(name)
            {
                return Err(bad());
            }
            let data_start = local + 30 + name_len + le16(header, 28)?;
            if data_start
                .checked_add(compressed)
                .is_none_or(|end| end > directory_start)
            {
                return Err(bad());
            }
            thumbnail = Some(NativeThumbnail {
                index,
                count,
                local,
                central: end - directory.len(),
                data_start,
                size,
                compressed,
                encrypted_members: false,
            });
        }
        directory = directory
            .get(46 + name_len + extra_len + comment_len..)
            .ok_or_else(bad)?;
    }
    if !directory.is_empty() {
        return Err(bad());
    }
    let mut thumbnail = thumbnail.ok_or_else(|| unsupported("Native archive has no exact root Thumbnail.png; numerical data and logs are not used as fallbacks"))?;
    thumbnail.encrypted_members = encrypted;
    Ok(thumbnail)
}

fn native_thumbnail(bytes: &[u8]) -> Result<(Vec<u8>, Vec<String>), EngineError> {
    if !bytes.starts_with(b"PK\x03\x04") {
        return Err(unsupported(
            "Unknown native source container; only ZIP-based D5000 Thumbnail.png import is supported",
        ));
    }
    let thumbnail = native_directory(bytes)?;
    let mut archive = zip::ZipArchive::new(Cursor::new(bytes))
        .map_err(|_| EngineError::invalid_input("Invalid native ZIP archive"))?;
    if archive.len() != thumbnail.count {
        return Err(EngineError::invalid_input(
            "Native ZIP directory names were reinterpreted or deduplicated",
        ));
    }
    {
        let member = archive
            .by_index_raw(thumbnail.index)
            .map_err(|_| EngineError::invalid_input("Invalid native thumbnail member"))?;
        if member.encrypted() {
            return Err(unsupported(
                "Encrypted native thumbnail is unsupported; no password is attempted",
            ));
        }
        if member.name() != "Thumbnail.png"
            || member.name_raw() != b"Thumbnail.png"
            || member.header_start() != thumbnail.local as u64
            || member.central_header_start() != thumbnail.central as u64
            || member.data_start() != thumbnail.data_start as u64
            || member.size() != thumbnail.size
            || member.compressed_size() != thumbnail.compressed as u64
            || member.is_dir()
            || member.is_symlink()
            || member.size() > MAX_GEL_IMAGE_BYTES as u64
            || !matches!(
                member.compression(),
                zip::CompressionMethod::Stored | zip::CompressionMethod::Deflated
            )
        {
            return Err(unsupported("Unsupported native thumbnail member"));
        }
    }
    let member = archive
        .by_index(thumbnail.index)
        .map_err(|_| EngineError::invalid_input("Cannot open unencrypted native thumbnail"))?;
    let expected = member.size() as usize;
    let mut raster = Vec::with_capacity(expected.min(1024 * 1024));
    member
        .take(MAX_GEL_IMAGE_BYTES as u64 + 1)
        .read_to_end(&mut raster)
        .map_err(|_| EngineError::invalid_input("Invalid native thumbnail data or checksum"))?;
    if raster.len() != expected
        || raster.len() > MAX_GEL_IMAGE_BYTES
        || !raster.starts_with(b"\x89PNG\r\n\x1a\n")
    {
        return Err(EngineError::invalid_input(
            "Native Thumbnail.png must be a bounded PNG matching its declared size",
        ));
    }
    check_raster_header(&raster, ImageFormat::Png)?;
    let mut warnings = vec![
        "Explicitly opted-in native Thumbnail.png: display-only evidence, not raw acquisition data or a numerical electropherogram.".into(),
        "Native numerical data, including encrypted numerical data, are unsupported. No XML schema, calibration, sample metadata or curves were guessed; no password was attempted.".into(),
        "Only Thumbnail.png was read. Operational logs and all other member payloads were left unopened; the exact source archive is retained locally in the source record.".into(),
    ];
    if thumbnail.encrypted_members {
        warnings.push("ZIP directory flags report encrypted members; those members were not opened or decrypted.".into());
    }
    Ok((raster, warnings))
}

#[cfg(test)]
mod tests {
    use super::*;
    use flate2::{Compression, write::ZlibEncoder};
    use lopdf::dictionary;
    use zip::{CompressionMethod, ZipWriter, write::SimpleFileOptions};

    // Synthetic containers are assembled here, deterministically and entirely
    // in memory, to exercise import boundaries. No vendor XML schema, password,
    // metadata, operational log or biological data is fabricated as evidence.
    fn request(path: &str) -> GelSourceImportRequest {
        GelSourceImportRequest {
            image_id: "source-test".into(),
            path: path.into(),
            pdf_object: None,
            allow_native_thumbnail: false,
        }
    }

    fn deflate(bytes: &[u8]) -> Vec<u8> {
        let mut encoder = ZlibEncoder::new(Vec::new(), Compression::default());
        encoder.write_all(bytes).unwrap();
        encoder.finish().unwrap()
    }

    fn raster(bytes: &[u8], rgb: bool) -> Stream {
        Stream::new(
            dictionary! {
                "Type" => "XObject", "Subtype" => "Image",
                "Width" => 2, "Height" => 2, "BitsPerComponent" => 8,
                "ColorSpace" => if rgb { "DeviceRGB" } else { "DeviceGray" },
                "Filter" => "FlateDecode",
            },
            deflate(bytes),
        )
    }

    fn document(images: Vec<Stream>) -> (Document, Vec<ObjectId>) {
        let mut doc = Document::with_version("1.4");
        // lopdf defaults to xref streams even for a 1.4 version string. These
        // fixtures deliberately exercise the supported classic-table subset.
        doc.reference_table.cross_reference_type = lopdf::xref::XrefType::CrossReferenceTable;
        let ids: Vec<_> = images
            .into_iter()
            .map(|image| doc.add_object(image))
            .collect();
        let pages = doc.new_object_id();
        let mut xobjects = Dictionary::new();
        for (index, id) in ids.iter().enumerate() {
            xobjects.set(format!("Im{index}"), *id);
        }
        let page = doc.add_object(dictionary! {
            "Type" => "Page", "Parent" => pages,
            "MediaBox" => vec![0.into(), 0.into(), 100.into(), 100.into()],
            "Resources" => dictionary! { "XObject" => xobjects },
        });
        doc.objects.insert(
            pages,
            Object::Dictionary(dictionary! {
                "Type" => "Pages", "Kids" => vec![Object::Reference(page)], "Count" => 1,
            }),
        );
        let root = doc.add_object(dictionary! { "Type" => "Catalog", "Pages" => pages });
        doc.trailer.set("Root", root);
        (doc, ids)
    }

    fn pdf(mut doc: Document) -> Vec<u8> {
        let mut bytes = Vec::new();
        doc.save_to(&mut bytes).unwrap();
        bytes
    }

    fn png() -> Vec<u8> {
        let mut bytes = Vec::new();
        image::codecs::png::PngEncoder::new(&mut bytes)
            .write_image(&[10, 20, 30, 40], 2, 2, ColorType::L8.into())
            .unwrap();
        bytes
    }

    fn decoded(record: &GelImageRecord) -> image::DynamicImage {
        image::load_from_memory(&STANDARD.decode(&record.original_base64).unwrap()).unwrap()
    }

    fn native(entries: &[(&str, &[u8], bool)]) -> Vec<u8> {
        native_with_method(entries, CompressionMethod::Stored)
    }

    fn native_with_method(entries: &[(&str, &[u8], bool)], method: CompressionMethod) -> Vec<u8> {
        let mut writer = ZipWriter::new(Cursor::new(Vec::new()));
        for (name, bytes, _) in entries {
            writer
                .start_file(
                    name,
                    SimpleFileOptions::default().compression_method(method),
                )
                .unwrap();
            writer.write_all(bytes).unwrap();
        }
        let mut bytes = writer.finish().unwrap().into_inner();
        // Deliberately flag bogus encrypted payloads without producing real
        // encrypted data. A successful thumbnail import proves those members
        // were not opened, decrypted or used for guessed numeric traces.
        let end = bytes
            .windows(4)
            .rposition(|part| part == b"PK\x05\x06")
            .unwrap();
        let mut at = le32(&bytes[end..], 16).unwrap();
        for (_, _, encrypted) in entries {
            let name_len = le16(&bytes[at..], 28).unwrap();
            let extra_len = le16(&bytes[at..], 30).unwrap();
            let comment_len = le16(&bytes[at..], 32).unwrap();
            if *encrypted {
                let local = le32(&bytes[at..], 42).unwrap();
                bytes[at + 8] |= 1;
                bytes[local + 6] |= 1;
            }
            at += 46 + name_len + extra_len + comment_len;
        }
        bytes
    }

    fn expect_code(request: &GelSourceImportRequest, bytes: &[u8], code: ErrorCode) {
        assert_eq!(import_source(request, bytes).unwrap_err().code, code);
    }

    #[test]
    fn unique_pdf_preserves_source_and_exact_rgb_pixels() {
        let pixels = [255, 0, 2, 3, 250, 5, 6, 7, 240, 9, 10, 11];
        let (doc, _) = document(vec![raster(&pixels, true)]);
        let bytes = pdf(doc);
        let (source, image) = import_source(&request("original report.pdf"), &bytes).unwrap();
        assert_eq!(source.source_name, "original report.pdf");
        assert_eq!(STANDARD.decode(&source.original_base64).unwrap(), bytes);
        assert_eq!(source.source_sha256, sha256_hex_bytes(&bytes));
        assert_eq!(source.extracted_image_sha256, image.descriptor.sha256);
        assert_eq!((image.descriptor.width, image.descriptor.height), (2, 2));
        assert_eq!(decoded(&image).to_rgb8().as_raw(), &pixels);
        let again = import_source(&request("original report.pdf"), &bytes).unwrap();
        assert_eq!(source, again.0);
        assert_eq!(image.original_base64, again.1.original_base64);
    }

    #[test]
    fn gray_pdf_and_singleton_filter_array_are_supported() {
        let mut stream = raster(&[1, 20, 150, 255], false);
        stream
            .dict
            .set("Filter", vec![Object::Name(b"FlateDecode".to_vec())]);
        stream.dict.set("DecodeParms", vec![Object::Null]);
        let (doc, _) = document(vec![stream]);
        let (_, image) = import_source(&request("gray.pdf"), &pdf(doc)).unwrap();
        assert_eq!(decoded(&image).to_luma8().as_raw(), &[1, 20, 150, 255]);
    }

    #[test]
    fn ambiguity_requires_explicit_object_and_never_chooses_a_mask() {
        let (doc, ids) = document(vec![raster(&[1; 12], true), raster(&[2; 12], true)]);
        let bytes = pdf(doc);
        expect_code(&request("ambiguous.pdf"), &bytes, ErrorCode::InvalidInput);
        let mut selected = request("ambiguous.pdf");
        selected.pdf_object = Some(ids[1].0);
        assert_eq!(
            decoded(&import_source(&selected, &bytes).unwrap().1)
                .to_rgb8()
                .as_raw(),
            &[2; 12]
        );
        selected.pdf_object = Some(999);
        expect_code(&selected, &bytes, ErrorCode::InvalidInput);
    }

    #[test]
    fn soft_mask_is_preserved_not_dropped_or_composited() {
        let (mut doc, ids) = document(vec![
            raster(&[100; 12], true),
            raster(&[0, 17, 128, 255], false),
        ]);
        doc.objects
            .get_mut(&ids[0])
            .unwrap()
            .as_stream_mut()
            .unwrap()
            .dict
            .set("SMask", ids[1]);
        let bytes = pdf(doc);
        let (source, image) = import_source(&request("alpha.pdf"), &bytes).unwrap();
        assert_eq!(
            decoded(&image).to_rgba8().as_raw(),
            &[
                100, 100, 100, 0, 100, 100, 100, 17, 100, 100, 100, 128, 100, 100, 100, 255
            ]
        );
        assert!(source.extraction.contains("soft_mask="));
        assert!(
            source
                .warnings
                .iter()
                .any(|warning| warning.contains("preserved as PNG alpha"))
        );
        let mut selected = request("alpha.pdf");
        selected.pdf_object = Some(ids[1].0);
        expect_code(&selected, &bytes, ErrorCode::Unsupported);
    }

    #[test]
    fn unsupported_masks_transforms_and_formats_are_not_silently_dropped() {
        let changes = [
            ("Mask", Object::Array(vec![0.into(), 1.into()])),
            ("Decode", Object::Array(vec![1.into(), 0.into()])),
            (
                "DecodeParms",
                Object::Dictionary(dictionary! { "Predictor" => 12 }),
            ),
            ("BitsPerComponent", Object::Integer(16)),
            ("ColorSpace", Object::Name(b"DeviceCMYK".to_vec())),
            ("ImageMask", Object::Boolean(true)),
            ("SMaskInData", Object::Integer(1)),
            ("F", Object::string_literal("not-accessed.png")),
            (
                "Filter",
                Object::Array(vec![
                    Object::from("FlateDecode"),
                    Object::from("FlateDecode"),
                ]),
            ),
        ];
        for (key, value) in changes {
            let mut image = raster(&[1; 12], true);
            image.dict.set(key, value);
            let (doc, ids) = document(vec![image]);
            let mut selected = request("unsupported.pdf");
            selected.pdf_object = Some(ids[0].0);
            expect_code(&selected, &pdf(doc), ErrorCode::Unsupported);
        }
        for key in ["Matte", "SMask"] {
            let (mut doc, ids) = document(vec![raster(&[1; 12], true), raster(&[255; 4], false)]);
            doc.objects
                .get_mut(&ids[0])
                .unwrap()
                .as_stream_mut()
                .unwrap()
                .dict
                .set("SMask", ids[1]);
            doc.objects
                .get_mut(&ids[1])
                .unwrap()
                .as_stream_mut()
                .unwrap()
                .dict
                .set(key, Object::Array(vec![0.into()]));
            expect_code(&request("bad-mask.pdf"), &pdf(doc), ErrorCode::Unsupported);
        }
    }

    #[test]
    fn malformed_flate_has_exact_decoded_length_and_no_unbounded_read() {
        let mut truncated = deflate(&[1; 12]);
        truncated.pop();
        let mut corrupt = deflate(&[1; 12]);
        *corrupt.last_mut().unwrap() ^= 1;
        for payload in [
            deflate(&[1; 11]),
            deflate(&[1; 13]),
            deflate(&[0; 65536]),
            vec![0; 4],
            truncated,
            corrupt,
        ] {
            let mut stream = raster(&[1; 12], true);
            stream.set_content(payload);
            let (doc, _) = document(vec![stream]);
            expect_code(
                &request("bad-flate.pdf"),
                &pdf(doc),
                ErrorCode::InvalidInput,
            );
        }
        let mut payload = deflate(&[1; 12]);
        payload.extend_from_slice(b"trailing bytes");
        let mut stream = raster(&[1; 12], true);
        stream.set_content(payload);
        let (doc, _) = document(vec![stream]);
        expect_code(&request("trailing.pdf"), &pdf(doc), ErrorCode::InvalidInput);
    }

    #[test]
    fn explicit_pdf_selection_does_not_decode_other_display_objects() {
        let mut unsupported_image = raster(&[0; 4], false);
        unsupported_image.dict.set("Filter", "LZWDecode");
        unsupported_image.set_content(b"intentionally not a compressed stream".to_vec());
        let mut oversized_image = raster(&[0; 4], false);
        oversized_image.dict.set("Width", 32768);
        oversized_image.dict.set("Height", 32768);
        let (doc, ids) = document(vec![
            raster(&[23; 12], true),
            unsupported_image,
            oversized_image,
        ]);
        let mut selected = request("mixed-display.pdf");
        selected.pdf_object = Some(ids[0].0);
        let (_, image) = import_source(&selected, &pdf(doc)).unwrap();
        assert_eq!(decoded(&image).to_rgb8().as_raw(), &[23; 12]);
    }

    #[test]
    fn jpeg_is_bounded_and_dictionary_dimensions_must_match() {
        let mut jpeg = Vec::new();
        image::codecs::jpeg::JpegEncoder::new(&mut jpeg)
            .encode(&[20, 40, 60, 80], 2, 2, ColorType::L8.into())
            .unwrap();
        let expected = image::load_from_memory(&jpeg).unwrap().to_luma8();
        let mut stream = raster(&[0; 4], false);
        stream.dict.set("Filter", "DCTDecode");
        stream.set_content(jpeg);
        let (doc, _) = document(vec![stream.clone()]);
        let (source, image) = import_source(&request("jpeg.pdf"), &pdf(doc)).unwrap();
        assert_eq!(decoded(&image).to_luma8(), expected);
        assert!(
            source
                .warnings
                .iter()
                .any(|warning| warning.contains("lossy"))
        );
        stream.dict.set("Width", 3);
        let (doc, _) = document(vec![stream]);
        expect_code(
            &request("mismatched.pdf"),
            &pdf(doc),
            ErrorCode::Unsupported,
        );
    }

    #[test]
    fn source_geometry_objects_pages_and_candidates_are_bounded() {
        expect_code(
            &request("large.pdf"),
            &vec![0; MAX_GEL_IMAGE_BYTES + 1],
            ErrorCode::InvalidInput,
        );
        expect_code(&request("empty.pdf"), &[], ErrorCode::InvalidInput);
        let mut stream = raster(&[0; 4], false);
        stream.dict.set("Width", 32768);
        stream.dict.set("Height", 32768);
        let (doc, ids) = document(vec![stream]);
        let mut selected = request("oversized.pdf");
        selected.pdf_object = Some(ids[0].0);
        expect_code(&selected, &pdf(doc), ErrorCode::InvalidInput);

        let (mut doc, _) = document(vec![raster(&[0; 4], false)]);
        for object in doc.objects.values_mut() {
            if let Object::Dictionary(dict) = object {
                if dict.get(b"Type").and_then(Object::as_name).ok() == Some(b"Pages") {
                    dict.set("Count", MAX_PDF_PAGES as i64 + 1);
                }
            }
        }
        expect_code(&request("pages.pdf"), &pdf(doc), ErrorCode::InvalidInput);

        let (doc, _) = document(
            (0..=MAX_PDF_CANDIDATES)
                .map(|_| raster(&[0; 4], false))
                .collect(),
        );
        expect_code(
            &request("candidates.pdf"),
            &pdf(doc),
            ErrorCode::InvalidInput,
        );
        let (mut doc, _) = document(vec![raster(&[0; 4], false)]);
        for _ in 0..MAX_PDF_OBJECTS {
            doc.add_object(Object::Null);
        }
        expect_code(&request("objects.pdf"), &pdf(doc), ErrorCode::InvalidInput);
        let (doc, _) = document(vec![]);
        expect_code(&request("no-images.pdf"), &pdf(doc), ErrorCode::Unsupported);
    }

    #[test]
    fn encrypted_and_eagerly_decoded_pdf_tables_are_rejected_before_loading() {
        let (mut doc, _) = document(vec![raster(&[0; 4], false)]);
        // Deliberately invalid encryption dictionary: rejection must happen
        // before any automatic empty-password authentication is attempted.
        doc.trailer
            .set("Encrypt", Object::Dictionary(Dictionary::new()));
        expect_code(&request("encrypted.pdf"), &pdf(doc), ErrorCode::Unsupported);
        for name in [
            b"/Encrypt".as_slice(),
            b"/Encr#79pt",
            b"/#45ncrypt",
            b"/ObjStm",
            b"/XRef",
            b"/XRefStm",
            b"/Prev",
        ] {
            assert_eq!(
                guard_pdf_names(name).unwrap_err().code,
                ErrorCode::Unsupported
            );
        }
        let (doc, _) = document(vec![raster(&[0; 4], false)]);
        let mut malformed = pdf(doc);
        let at = malformed
            .windows(4)
            .position(|part| part == b"xref")
            .unwrap();
        malformed[at] = b'!';
        assert!(import_source(&request("bad-table.pdf"), &malformed).is_err());
    }

    #[test]
    fn real_tapestation_pdfs_match_retained_png_rgb_pixels_without_reencoding_assumptions() {
        // Real public owner-supplied reports and Poppler-extracted PNGs;
        // provenance and reproduction are in the existing fixture README.
        let cases: &[(&str, &[u8], &[u8])] = &[
            (
                "before.pdf",
                include_bytes!(
                    "../../test_files/fixtures/gel_images/tapestation_26r266/before.pdf"
                ),
                include_bytes!(
                    "../../test_files/fixtures/gel_images/tapestation_26r266/before.png"
                ),
            ),
            (
                "after.pdf",
                include_bytes!("../../test_files/fixtures/gel_images/tapestation_26r266/after.pdf"),
                include_bytes!("../../test_files/fixtures/gel_images/tapestation_26r266/after.png"),
            ),
        ];
        for (name, bytes, expected_png) in cases {
            expect_code(&request(name), bytes, ErrorCode::InvalidInput);
            let mut selected = request(name);
            selected.pdf_object = Some(5);
            let (source, image) = import_source(&selected, bytes).unwrap();
            assert_eq!(
                (image.descriptor.width, image.descriptor.height),
                (1375, 731)
            );
            assert_eq!(STANDARD.decode(&source.original_base64).unwrap(), *bytes);
            assert_eq!(source.source_sha256, sha256_hex_bytes(bytes));
            assert_eq!(source.extracted_image_sha256, image.descriptor.sha256);
            assert!(source.extraction.contains("object=5:generation=0"));
            let actual = decoded(&image);
            assert_eq!(actual.color(), ColorType::Rgba8);
            let expected = image::load_from_memory(expected_png).unwrap().to_rgb8();
            assert_eq!(
                actual.to_rgb8().as_raw(),
                expected.as_raw(),
                "{name}: compare decoded RGB samples, never PNG encoding hashes"
            );
        }
    }

    #[test]
    fn native_thumbnail_opt_in_retains_source_and_leaves_private_payloads_unread() {
        let thumbnail = png();
        let bytes = native(&[
            ("ModelData.xml", b"not encrypted data or XML", true),
            ("RunData.xml", b"not encrypted data or XML", true),
            ("Logs.zip", b"not a ZIP -- private operational data", false),
            ("Thumbnail.png", &thumbnail, false),
        ]);
        let mut selected = request("original.D5000");
        expect_code(&selected, &bytes, ErrorCode::Unsupported);
        selected.allow_native_thumbnail = true;
        let (source, image) = import_source(&selected, &bytes).unwrap();
        assert_eq!(source.source_name, "original.D5000");
        assert_eq!(STANDARD.decode(&source.original_base64).unwrap(), bytes);
        assert_eq!(STANDARD.decode(&image.original_base64).unwrap(), thumbnail);
        assert_eq!(source.source_sha256, sha256_hex_bytes(&bytes));
        assert_eq!(source.extracted_image_sha256, sha256_hex_bytes(&thumbnail));
        assert!(
            source
                .warnings
                .iter()
                .any(|warning| warning.contains("encrypted members"))
        );
        assert!(
            source
                .warnings
                .iter()
                .any(|warning| warning.contains("display-only"))
        );
        assert!(
            !source
                .warnings
                .join(" ")
                .contains("private operational data")
        );
    }

    #[test]
    fn native_missing_encrypted_nested_and_invalid_thumbnails_have_no_fallback() {
        let thumbnail = png();
        let mut selected = request("run.D5000");
        selected.allow_native_thumbnail = true;
        for entries in [
            vec![("ModelData.xml", b"unread".as_slice(), true)],
            vec![("nested/Thumbnail.png", thumbnail.as_slice(), false)],
            vec![("Thumbnail.png", thumbnail.as_slice(), true)],
        ] {
            expect_code(&selected, &native(&entries), ErrorCode::Unsupported);
        }
        let error = import_source(
            &selected,
            &native(&[("Thumbnail.png", b"private non-image contents", false)]),
        )
        .unwrap_err();
        assert_eq!(error.code, ErrorCode::InvalidInput);
        assert!(!error.message.contains("private non-image contents"));
        expect_code(
            &selected,
            b"unknown proprietary container",
            ErrorCode::Unsupported,
        );
        expect_code(
            &request("run.unknown"),
            &native(&[("Thumbnail.png", &thumbnail, false)]),
            ErrorCode::Unsupported,
        );
    }

    #[test]
    fn native_directory_is_bounded_before_zip_library_allocation_and_name_deduplication() {
        let thumbnail = png();
        let mut selected = request("run.D5000");
        selected.allow_native_thumbnail = true;
        let bytes = native(&[("Thumbnail.png", &thumbnail, false)]);
        let end = bytes
            .windows(4)
            .rposition(|part| part == b"PK\x05\x06")
            .unwrap();
        let at = le32(&bytes[end..], 16).unwrap();
        let mut oversized = bytes.clone();
        oversized[at + 24..at + 28]
            .copy_from_slice(&(MAX_GEL_IMAGE_BYTES as u32 + 1).to_le_bytes());
        expect_code(&selected, &oversized, ErrorCode::InvalidInput);
        let mut excessive = bytes.clone();
        for offset in [end + 8, end + 10] {
            excessive[offset..offset + 2]
                .copy_from_slice(&((MAX_NATIVE_MEMBERS + 1) as u16).to_le_bytes());
        }
        expect_code(&selected, &excessive, ErrorCode::InvalidInput);
        let mut inconsistent = bytes;
        inconsistent[6] |= 1;
        expect_code(&selected, &inconsistent, ErrorCode::InvalidInput);

        let mut duplicate = native(&[
            ("Thumbnail.png", &thumbnail, false),
            ("Thumbnaim.png", &thumbnail, false),
        ]);
        let positions: Vec<_> = duplicate
            .windows(13)
            .enumerate()
            .filter(|(_, bytes)| *bytes == b"Thumbnaim.png")
            .map(|(index, _)| index)
            .collect();
        assert_eq!(positions.len(), 2);
        for at in positions {
            duplicate[at..at + 13].copy_from_slice(b"Thumbnail.png");
        }
        expect_code(&selected, &duplicate, ErrorCode::InvalidInput);
    }

    #[test]
    fn ordinary_raster_is_retained_and_selection_options_are_not_ignored() {
        let bytes = png();
        let (source, image) = import_source(&request("gel.png"), &bytes).unwrap();
        assert_eq!(source.source_sha256, source.extracted_image_sha256);
        assert_eq!(source.original_base64, image.original_base64);
        let mut selected = request("gel.png");
        selected.pdf_object = Some(5);
        expect_code(&selected, &bytes, ErrorCode::InvalidInput);
        expect_code(
            &request("unknown.native"),
            b"unknown data",
            ErrorCode::Unsupported,
        );
    }

    #[test]
    fn native_deflate_is_bounded_even_when_directory_size_understates_payload() {
        let thumbnail = png();
        let valid = native_with_method(
            &[("Thumbnail.png", &thumbnail, false)],
            CompressionMethod::Deflated,
        );
        let mut selected = request("deflated.D5000");
        selected.allow_native_thumbnail = true;
        assert_eq!(
            STANDARD
                .decode(&import_source(&selected, &valid).unwrap().1.original_base64)
                .unwrap(),
            thumbnail
        );

        let bomb = vec![0; MAX_GEL_IMAGE_BYTES + 8192];
        let mut bytes = native_with_method(
            &[("Thumbnail.png", &bomb, false)],
            CompressionMethod::Deflated,
        );
        let end = bytes
            .windows(4)
            .rposition(|part| part == b"PK\x05\x06")
            .unwrap();
        let central = le32(&bytes[end..], 16).unwrap();
        bytes[central + 24..central + 28].copy_from_slice(&4u32.to_le_bytes());
        bytes[22..26].copy_from_slice(&4u32.to_le_bytes());
        expect_code(&selected, &bytes, ErrorCode::InvalidInput);
    }
}
