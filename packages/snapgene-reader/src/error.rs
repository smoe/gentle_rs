use crate::PacketType;
use thiserror::Error;

/// Errors returned while parsing or adapting SnapGene `.dna` files.
#[derive(Debug, Error)]
pub enum SnapGeneError {
    #[error("Could not read SnapGene file '{path}': {source}")]
    Io {
        path: String,
        #[source]
        source: std::io::Error,
    },
    #[error("SnapGene file is empty and does not contain a cookie packet")]
    EmptyFile,
    #[error("SnapGene file does not start with a cookie packet; found {found:?}")]
    MissingCookie { found: PacketType },
    #[error("SnapGene cookie packet is too short: expected 14 bytes, got {actual}")]
    ShortCookie { actual: usize },
    #[error("SnapGene cookie mismatch: expected 'SnapGene', found '{found}'")]
    InvalidCookie { found: String },
    #[error("Unexpected end of input while reading {context}")]
    UnexpectedEof { context: String },
    #[error("SnapGene file contains more than one DNA packet")]
    DuplicateDnaPacket,
    #[error("SnapGene file does not contain a DNA packet")]
    MissingDnaPacket,
    #[error("Packet {packet_type:?} payload is invalid UTF-8: {source}")]
    InvalidUtf8 {
        packet_type: PacketType,
        #[source]
        source: std::string::FromUtf8Error,
    },
    #[error("Could not parse {context} XML: {message}")]
    Xml {
        context: &'static str,
        message: String,
    },
    #[error("Feature segment has invalid range '{range}'")]
    InvalidRange { range: String },
    #[error("Feature '{feature_name}' is missing a segment range")]
    MissingFeatureRange { feature_name: String },
    #[error("Could not build {context} location '{location}': {message}")]
    InvalidLocation {
        context: &'static str,
        location: String,
        message: String,
    },
}

