//! Strict bounded four-line FASTQ intake; raw compressed bytes are hashed as parsed.

use super::*;
use flate2::read::MultiGzDecoder;
use md5::{Digest, Md5};
use std::cell::RefCell;
use std::io::{self, BufRead, BufReader, Read, Write};
use std::rc::Rc;

pub(super) const MAX_READ_BASES: usize = 2_000_000;

pub(super) struct FastqRecord {
    pub read_id: String,
    pub sequence: Vec<u8>,
    pub quality: Vec<u8>,
    pub logical_offset: u64,
    pub compressed_bytes_read: u64,
    pub compressed_bytes_total: u64,
}

#[derive(Clone)]
struct HashState {
    sha: ring::digest::Context,
    md5: Md5,
    bytes: u64,
}

struct RawReader {
    file: fs::File,
    state: Rc<RefCell<HashState>>,
}

impl Read for RawReader {
    fn read(&mut self, buffer: &mut [u8]) -> io::Result<usize> {
        let n = self.file.read(buffer)?;
        let mut state = self.state.borrow_mut();
        state.sha.update(&buffer[..n]);
        state.md5.update(&buffer[..n]);
        state.bytes += n as u64;
        Ok(n)
    }
}

#[derive(Debug)]
pub(super) struct InputDigest {
    pub sha256: String,
    pub md5: String,
    pub bytes: u64,
}

pub(super) fn hex(bytes: &[u8]) -> String {
    bytes.iter().map(|byte| format!("{byte:02x}")).collect()
}

fn line(reader: &mut dyn BufRead, offset: &mut u64) -> Result<Option<Vec<u8>>, EngineError> {
    let mut bytes = Vec::new();
    // Bound even a malformed file with no newlines; never read an unlimited line.
    let n = reader
        .take((MAX_READ_BASES + 3) as u64)
        .read_until(b'\n', &mut bytes)
        .map_err(|e| bad(format!("fastq_io: {e}")))?;
    if n == 0 {
        return Ok(None);
    }
    if bytes.len() > MAX_READ_BASES + 2 {
        return Err(bad("fastq_record_budget: line exceeds 2000000 bases"));
    }
    *offset += n as u64;
    if bytes.last() == Some(&b'\n') {
        bytes.pop();
        if bytes.last() == Some(&b'\r') {
            bytes.pop();
        }
    }
    if bytes.len() > MAX_READ_BASES {
        return Err(bad("fastq_record_budget: line exceeds 2000000 bases"));
    }
    Ok(Some(bytes))
}

/// None means cancelled, never a verified EOF or a checksum of the whole file.
pub(super) fn visit_fastq(
    path: &Path,
    on_record: &mut dyn FnMut(FastqRecord) -> Result<bool, EngineError>,
) -> Result<Option<InputDigest>, EngineError> {
    let file =
        fs::File::open(path).map_err(|e| bad(format!("fastq_open: {}: {e}", path.display())))?;
    let before = file
        .metadata()
        .map_err(|e| bad(format!("fastq_metadata: {e}")))?;
    if !before.is_file() {
        return Err(bad("fastq_input_not_file"));
    }
    let state = Rc::new(RefCell::new(HashState {
        sha: ring::digest::Context::new(&ring::digest::SHA256),
        md5: Md5::new(),
        bytes: 0,
    }));
    let raw = RawReader {
        file,
        state: state.clone(),
    };
    let mut probe = BufReader::new(raw);
    let gzip = probe
        .fill_buf()
        .map_err(|e| bad(format!("fastq_probe: {e}")))?
        .starts_with(&[0x1f, 0x8b]);
    let mut reader: Box<dyn BufRead> = if gzip {
        Box::new(BufReader::new(MultiGzDecoder::new(probe)))
    } else {
        Box::new(probe)
    };
    let mut offset = 0;
    loop {
        let logical_offset = offset;
        let Some(header) = line(reader.as_mut(), &mut offset)? else {
            break;
        };
        if header.first() != Some(&b'@') || header.len() < 2 {
            return Err(bad(
                "fastq_header: expected @ plus a read identifier (FASTA is not FASTQ)",
            ));
        }
        let mut sequence = line(reader.as_mut(), &mut offset)?
            .ok_or_else(|| bad("fastq_truncated: missing sequence"))?;
        let plus = line(reader.as_mut(), &mut offset)?
            .ok_or_else(|| bad("fastq_truncated: missing + line"))?;
        let quality = line(reader.as_mut(), &mut offset)?
            .ok_or_else(|| bad("fastq_truncated: missing quality"))?;
        if plus.first() != Some(&b'+') || sequence.is_empty() || sequence.len() != quality.len() {
            return Err(bad(
                "fastq_record: require +, nonempty sequence and equal sequence/quality lengths",
            ));
        }
        if !quality.iter().all(|&byte| (33..=126).contains(&byte)) {
            return Err(bad("fastq_quality: expected printable Phred+33 bytes"));
        }
        let header = std::str::from_utf8(&header[1..]).map_err(|_| bad("fastq_header_utf8"))?;
        let read_id = header
            .split_ascii_whitespace()
            .next()
            .filter(|s| !s.is_empty())
            .ok_or_else(|| bad("fastq_empty_read_id"))?;
        if plus.len() > 1 {
            let repeated = std::str::from_utf8(&plus[1..]).map_err(|_| bad("fastq_plus_utf8"))?;
            if repeated.split_ascii_whitespace().next() != Some(read_id) {
                return Err(bad("fastq_plus: repeated identifier does not match header"));
            }
        }
        for byte in &mut sequence {
            *byte = byte.to_ascii_uppercase();
            if *byte == b'U' {
                *byte = b'T';
            }
            if !b"ACGTRYSWKMBDHVN".contains(byte) {
                return Err(bad("fastq_sequence: non-IUPAC base"));
            }
        }
        if !on_record(FastqRecord {
            read_id: read_id.into(),
            sequence,
            quality,
            logical_offset,
            compressed_bytes_read: state.borrow().bytes,
            compressed_bytes_total: before.len(),
        })? {
            return Ok(None);
        }
    }
    let after = fs::metadata(path).map_err(|e| bad(format!("fastq_metadata_after: {e}")))?;
    let state = state.borrow();
    if before.len() != after.len()
        || before.modified().ok() != after.modified().ok()
        || state.bytes != before.len()
    {
        return Err(bad(
            "fastq_changed_during_read: discard this file and retry after its transfer finishes",
        ));
    }
    Ok(Some(InputDigest {
        sha256: format!("sha256:{}", hex(state.sha.clone().finish().as_ref())),
        md5: hex(&state.md5.clone().finalize()),
        bytes: state.bytes,
    }))
}

/// Digest exactly the successfully written artifact bytes, without rereading a large sidecar.
pub(super) struct HashedWriter {
    pub file: fs::File,
    sha: ring::digest::Context,
    bytes: u64,
}

impl HashedWriter {
    pub fn create(path: &Path) -> Result<Self, EngineError> {
        Ok(Self {
            file: fs::OpenOptions::new()
                .write(true)
                .create_new(true)
                .open(path)
                .map_err(|e| bad(format!("panel_artifact_create: {e}")))?,
            sha: ring::digest::Context::new(&ring::digest::SHA256),
            bytes: 0,
        })
    }
    pub fn receipt(mut self, path: String) -> Result<RnaReadPanelArtifact, EngineError> {
        self.flush()
            .map_err(|e| bad(format!("panel_artifact_flush: {e}")))?;
        self.file
            .sync_all()
            .map_err(|e| bad(format!("panel_artifact_sync: {e}")))?;
        Ok(RnaReadPanelArtifact {
            path,
            sha256: format!("sha256:{}", hex(self.sha.finish().as_ref())),
            bytes: self.bytes,
            locus_id: None,
        })
    }
}

impl std::io::Write for HashedWriter {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        let n = self.file.write(bytes)?;
        self.sha.update(&bytes[..n]);
        self.bytes += n as u64;
        Ok(n)
    }
    fn flush(&mut self) -> io::Result<()> {
        self.file.flush()
    }
}
