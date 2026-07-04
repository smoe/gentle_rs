//! Shared digest helpers for local fingerprints and deterministic identifiers.

use ring::digest::{Context, SHA256, digest};
use std::{
    fs::File,
    io::{self, Read},
    path::Path,
};

fn hex_encode(bytes: &[u8]) -> String {
    const HEX: &[u8; 16] = b"0123456789abcdef";
    let mut out = String::with_capacity(bytes.len() * 2);
    for byte in bytes {
        out.push(HEX[(byte >> 4) as usize] as char);
        out.push(HEX[(byte & 0x0f) as usize] as char);
    }
    out
}

pub(crate) fn sha256_hex_bytes(bytes: &[u8]) -> String {
    hex_encode(digest(&SHA256, bytes).as_ref())
}

pub(crate) fn sha256_hex_str(value: &str) -> String {
    sha256_hex_bytes(value.as_bytes())
}

pub(crate) fn sha256_prefixed_str(value: &str) -> String {
    format!("sha256:{}", sha256_hex_str(value))
}

pub(crate) fn sha256_prefixed_bytes(bytes: &[u8]) -> String {
    format!("sha256:{}", sha256_hex_bytes(bytes))
}

pub(crate) fn sha256_file_hex(path: &Path) -> io::Result<String> {
    let mut file = File::open(path)?;
    let mut context = Context::new(&SHA256);
    let mut buffer = [0u8; 8192];
    loop {
        let read = file.read(&mut buffer)?;
        if read == 0 {
            break;
        }
        context.update(&buffer[..read]);
    }
    Ok(hex_encode(context.finish().as_ref()))
}

pub(crate) fn short_sha256_id(prefix: &str, raw: &str) -> String {
    let digest = sha256_hex_str(raw);
    format!("{prefix}_{}", &digest[..16])
}
