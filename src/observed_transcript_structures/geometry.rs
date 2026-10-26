//! Shared operation geometry; no annotation, strand inference or input/output.

use crate::engine::EngineError;
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Block {
    pub query: (usize, usize),
    pub target: (usize, usize),
}

pub(crate) struct Geometry {
    pub blocks: Vec<Block>,
    pub skips: Vec<(usize, usize)>,
    pub edits: usize,
    pub matching: usize,
    pub columns: usize,
}

/// The Phase 2 walk, including its error strings and reference-order semantics.
pub(crate) fn walk(
    ops: impl IntoIterator<Item = (u8, usize)>,
    (qs, qe): (usize, usize),
    (ts, te): (usize, usize),
    strand: char,
) -> Result<Geometry, EngineError> {
    let bad = EngineError::invalid_input;
    if qs >= qe || ts >= te || !matches!(strand, '+' | '-') {
        return Err(bad("primer_payload_cigar_span"));
    }
    let (mut q, mut t, mut nm, mut matching, mut columns) = (0usize, ts, 0usize, 0usize, 0usize);
    let mut blocks = Vec::new();
    let mut skips = Vec::new();
    for (kind, len) in ops {
        if len == 0 || !b"MIDN=X".contains(&kind) {
            return Err(bad("primer_payload_cigar_operation"));
        }
        let qc = b"M=XI".contains(&kind);
        let tc = b"M=XDN".contains(&kind);
        let qnext = q
            .checked_add(if qc { len } else { 0 })
            .ok_or_else(|| bad("primer_payload_cigar_overflow"))?;
        let tnext = t
            .checked_add(if tc { len } else { 0 })
            .ok_or_else(|| bad("primer_payload_cigar_overflow"))?;
        if qnext > qe - qs || tnext > te {
            return Err(bad("primer_payload_cigar_span"));
        }
        if b"M=X".contains(&kind) {
            let query = if strand == '+' {
                (qs + q, qs + qnext)
            } else {
                (qe - qnext, qe - q)
            };
            blocks.push(Block {
                query,
                target: (t, tnext),
            });
        }
        if kind == b'N' {
            skips.push((t, tnext));
        }
        if b"XID".contains(&kind) {
            nm = nm
                .checked_add(len)
                .ok_or_else(|| bad("primer_payload_cigar_overflow"))?;
        }
        if kind == b'=' {
            matching = matching
                .checked_add(len)
                .ok_or_else(|| bad("primer_payload_cigar_overflow"))?;
        }
        if kind != b'N' {
            columns = columns
                .checked_add(len)
                .ok_or_else(|| bad("primer_payload_cigar_overflow"))?;
        }
        q = qnext;
        t = tnext;
    }
    if q != qe - qs || t != te {
        return Err(bad("primer_payload_cigar_span_or_counts"));
    }
    blocks.sort_by_key(|x| x.query);
    let mut merged: Vec<Block> = Vec::new();
    for block in blocks {
        if let Some(last) = merged.last_mut()
            && last.query.1 == block.query.0
            && (if strand == '+' {
                last.target.1 == block.target.0
            } else {
                last.target.0 == block.target.1
            })
        {
            last.query.1 = block.query.1;
            last.target = (
                last.target.0.min(block.target.0),
                last.target.1.max(block.target.1),
            );
        } else {
            merged.push(block);
        }
    }
    Ok(Geometry {
        blocks: merged,
        skips,
        edits: nm,
        matching,
        columns,
    })
}

/// Project a genomic intersection through one query-supported block.
pub(crate) fn project(
    block: &Block,
    range: (usize, usize),
    strand: char,
) -> Option<(usize, usize)> {
    let (a, z) = (block.target.0.max(range.0), block.target.1.min(range.1));
    if a >= z {
        return None;
    }
    Some(if strand == '+' {
        (
            block.query.0 + a - block.target.0,
            block.query.0 + z - block.target.0,
        )
    } else {
        (
            block.query.1 - (z - block.target.0),
            block.query.1 - (a - block.target.0),
        )
    })
}
