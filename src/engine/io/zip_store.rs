//! Minimal deterministic ZIP writer for generated OpenDocument/Office reports.
//!
//! GENtle only needs stored members for small self-contained report bundles, so
//! this module keeps ZIP generation auditable without depending on host tools.

#[derive(Debug, Clone)]
pub(crate) struct StoredZipMember {
    pub(crate) name: String,
    pub(crate) data: Vec<u8>,
}

impl StoredZipMember {
    pub(crate) fn new(name: impl Into<String>, data: impl Into<Vec<u8>>) -> Self {
        Self {
            name: name.into(),
            data: data.into(),
        }
    }
}

fn push_zip_u16(out: &mut Vec<u8>, value: u16) {
    out.extend_from_slice(&value.to_le_bytes());
}

fn push_zip_u32(out: &mut Vec<u8>, value: u32) {
    out.extend_from_slice(&value.to_le_bytes());
}

fn zip_crc32(bytes: &[u8]) -> u32 {
    let mut crc = 0xffff_ffffu32;
    for &byte in bytes {
        crc ^= byte as u32;
        for _ in 0..8 {
            let mask = (crc & 1).wrapping_neg();
            crc = (crc >> 1) ^ (0xedb8_8320 & mask);
        }
    }
    !crc
}

pub(crate) fn build_stored_zip_archive(members: &[StoredZipMember]) -> Result<Vec<u8>, String> {
    let mut out = Vec::new();
    let mut central = Vec::new();

    for member in members {
        let name_bytes = member.name.as_bytes();
        if name_bytes.len() > u16::MAX as usize {
            return Err(format!("ZIP member name '{}' is too long", member.name));
        }
        if member.data.len() > u32::MAX as usize {
            return Err(format!("ZIP member '{}' is too large", member.name));
        }
        if out.len() > u32::MAX as usize {
            return Err("ZIP local offset exceeds 32-bit ZIP limit".to_string());
        }

        let local_offset = out.len() as u32;
        let crc = zip_crc32(&member.data);
        let size = member.data.len() as u32;
        let name_len = name_bytes.len() as u16;

        push_zip_u32(&mut out, 0x0403_4b50);
        push_zip_u16(&mut out, 20);
        push_zip_u16(&mut out, 0);
        push_zip_u16(&mut out, 0);
        push_zip_u16(&mut out, 0);
        push_zip_u16(&mut out, 0);
        push_zip_u32(&mut out, crc);
        push_zip_u32(&mut out, size);
        push_zip_u32(&mut out, size);
        push_zip_u16(&mut out, name_len);
        push_zip_u16(&mut out, 0);
        out.extend_from_slice(name_bytes);
        out.extend_from_slice(&member.data);

        push_zip_u32(&mut central, 0x0201_4b50);
        push_zip_u16(&mut central, 20);
        push_zip_u16(&mut central, 20);
        push_zip_u16(&mut central, 0);
        push_zip_u16(&mut central, 0);
        push_zip_u16(&mut central, 0);
        push_zip_u16(&mut central, 0);
        push_zip_u32(&mut central, crc);
        push_zip_u32(&mut central, size);
        push_zip_u32(&mut central, size);
        push_zip_u16(&mut central, name_len);
        push_zip_u16(&mut central, 0);
        push_zip_u16(&mut central, 0);
        push_zip_u16(&mut central, 0);
        push_zip_u16(&mut central, 0);
        push_zip_u32(&mut central, 0);
        push_zip_u32(&mut central, local_offset);
        central.extend_from_slice(name_bytes);
    }

    if members.len() > u16::MAX as usize {
        return Err("ZIP archive has too many members".to_string());
    }
    if out.len() > u32::MAX as usize || central.len() > u32::MAX as usize {
        return Err("ZIP central directory exceeds 32-bit ZIP limit".to_string());
    }

    let central_offset = out.len() as u32;
    let central_size = central.len() as u32;
    out.extend_from_slice(&central);
    push_zip_u32(&mut out, 0x0605_4b50);
    push_zip_u16(&mut out, 0);
    push_zip_u16(&mut out, 0);
    push_zip_u16(&mut out, members.len() as u16);
    push_zip_u16(&mut out, members.len() as u16);
    push_zip_u32(&mut out, central_size);
    push_zip_u32(&mut out, central_offset);
    push_zip_u16(&mut out, 0);

    Ok(out)
}
