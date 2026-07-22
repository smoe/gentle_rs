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

/// Canonical physical oligo sequence used for cross-report identities and
/// procurement comparison: remove all whitespace, uppercase, and normalize
/// RNA `U` to DNA `T`.
pub(crate) fn canonical_oligo_sequence(raw: &str) -> String {
    raw.chars()
        .filter(|ch| !ch.is_whitespace())
        .map(|ch| match ch.to_ascii_uppercase() {
            'U' => 'T',
            other => other,
        })
        .collect()
}

pub(crate) fn oligo_sequence_sha256(raw: &str) -> String {
    sha256_prefixed_str(&canonical_oligo_sequence(raw))
}

pub(crate) fn oligo_full_id(raw: &str) -> String {
    format!(
        "oligo_sha256_{}",
        sha256_hex_str(&canonical_oligo_sequence(raw))
    )
}

pub(crate) fn oligo_tube_id(raw: &str) -> String {
    let digest = sha256_hex_str(&canonical_oligo_sequence(raw));
    format!("O-{}", digest[..8].to_ascii_uppercase())
}

/// Ordered physical primer-pair identity. Forward/reverse role is part of the
/// identity; legacy assay ids can additionally include assay kind and probe.
pub(crate) fn primer_pair_full_id(forward: &str, reverse: &str) -> String {
    let identity = format!(
        "forward={}|reverse={}",
        oligo_full_id(forward),
        oligo_full_id(reverse)
    );
    format!("pair_sha256_{}", sha256_hex_str(&identity))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn canonical_oligo_identity_normalizes_case_whitespace_and_u() {
        let variants = ["acgu", " A C G U ", "ACGT\n"];
        let ids = variants.map(oligo_full_id);
        assert_eq!(ids[0], ids[1]);
        assert_eq!(ids[1], ids[2]);
        assert_eq!(canonical_oligo_sequence(variants[0]), "ACGT");
        assert_eq!(oligo_tube_id(variants[0]), oligo_tube_id(variants[2]));
    }

    #[test]
    fn primer_pair_identity_is_ordered() {
        assert_ne!(
            primer_pair_full_id("ACGT", "TGCA"),
            primer_pair_full_id("TGCA", "ACGT")
        );
    }
}
