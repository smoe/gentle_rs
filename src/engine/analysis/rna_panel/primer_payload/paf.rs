//! Bounded, query-grouped PAF intake; tags are producer evidence, not an alignment search.

use super::*;

#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) struct Op {
    pub kind: u8,
    pub len: usize,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) struct Block {
    pub query: (usize, usize),
    pub target: (usize, usize),
}

#[derive(Debug, Clone)]
pub(super) struct Alignment {
    pub segment: PrimerPayloadSegment,
    pub blocks: Vec<Block>,
    pub skips: Vec<(usize, usize)>,
}

fn number(s: &str) -> Result<usize, EngineError> {
    if s.is_empty() || !s.bytes().all(|c| c.is_ascii_digit()) {
        return Err(bad("primer_payload_paf_integer"));
    }
    s.parse()
        .map_err(|_| bad("primer_payload_paf_integer_overflow"))
}

fn push(
    ops: &mut Vec<Op>,
    kind: u8,
    len: usize,
    b: &PrimerPayloadBudgets,
) -> Result<(), EngineError> {
    if len == 0 || ops.len() >= b.max_cigar_operations {
        return Err(bad("primer_payload_cigar_operation_budget"));
    }
    ops.push(Op { kind, len });
    Ok(())
}

pub(super) fn cigar(text: &str, b: &PrimerPayloadBudgets) -> Result<Vec<Op>, EngineError> {
    let mut ops = Vec::new();
    let mut start = 0;
    for (i, c) in text.bytes().enumerate() {
        if c.is_ascii_digit() {
            continue;
        }
        if !b"MIDN=X".contains(&c) {
            return Err(bad("primer_payload_cigar_operation"));
        }
        push(&mut ops, c, number(&text[start..i])?, b)?;
        start = i + 1;
    }
    if start != text.len() || ops.is_empty() {
        return Err(bad("primer_payload_cigar_trailing_or_empty"));
    }
    Ok(ops)
}

pub(super) fn difference(text: &str, b: &PrimerPayloadBudgets) -> Result<Vec<Op>, EngineError> {
    let bytes = text.as_bytes();
    let mut i = 0;
    let mut ops = Vec::new();
    let base = |x: u8| b"ACGTNacgtn".contains(&x);
    while i < bytes.len() {
        let token = bytes[i];
        i += 1;
        let start = i;
        let (kind, n) = match token {
            b':' => {
                while i < bytes.len() && bytes[i].is_ascii_digit() {
                    i += 1;
                }
                (b'=', number(&text[start..i])?)
            }
            b'=' | b'+' | b'-' => {
                while i < bytes.len() && base(bytes[i]) {
                    i += 1;
                }
                (
                    match token {
                        b'=' => b'=',
                        b'+' => b'I',
                        _ => b'D',
                    },
                    i - start,
                )
            }
            b'*' if i + 2 <= bytes.len() && bytes[i..i + 2].iter().all(|x| base(*x)) => {
                i += 2;
                (b'X', 1)
            }
            b'~' if i + 2 <= bytes.len() && bytes[i..i + 2].iter().all(|x| base(*x)) => {
                i += 2;
                let start = i;
                while i < bytes.len() && bytes[i].is_ascii_digit() {
                    i += 1;
                }
                let n = number(&text[start..i])?;
                if i + 2 > bytes.len() || !bytes[i..i + 2].iter().all(|x| base(*x)) {
                    return Err(bad("primer_payload_cs_skip"));
                }
                i += 2;
                (b'N', n)
            }
            _ => return Err(bad("primer_payload_cs_token")),
        };
        push(&mut ops, kind, n, b)?;
    }
    if ops.is_empty() {
        return Err(bad("primer_payload_cs_empty"));
    }
    Ok(ops)
}

fn normalized(ops: &[Op]) -> Result<Vec<Op>, EngineError> {
    let mut result: Vec<Op> = Vec::new();
    for op in ops {
        let kind = if b"=X".contains(&op.kind) {
            b'M'
        } else {
            op.kind
        };
        if let Some(last) = result.last_mut()
            && last.kind == kind
        {
            last.len = last
                .len
                .checked_add(op.len)
                .ok_or_else(|| bad("primer_payload_cigar_overflow"))?;
        } else {
            result.push(Op { kind, len: op.len });
        }
    }
    Ok(result)
}

fn detail_agrees(cg: &[Op], cs: &[Op]) -> bool {
    let (mut i, mut j, mut a, mut z) = (0usize, 0usize, 0usize, 0usize);
    while i < cg.len() && j < cs.len() {
        if cg[i].kind != cs[j].kind && !(cg[i].kind == b'M' && b"=X".contains(&cs[j].kind)) {
            return false;
        }
        let n = (cg[i].len - a).min(cs[j].len - z);
        a += n;
        z += n;
        if a == cg[i].len {
            i += 1;
            a = 0;
        }
        if z == cs[j].len {
            j += 1;
            z = 0;
        }
    }
    i == cg.len() && j == cs.len()
}

pub(super) fn row(
    text: &str,
    query: &PrimerPayloadQuery,
    strict: bool,
    b: &PrimerPayloadBudgets,
) -> Result<Option<(Alignment, usize)>, EngineError> {
    let f = text.split('\t').collect::<Vec<_>>();
    if f.len() < 12 || f[0] != query.query_id || number(f[1])? != query.query_length {
        return Err(bad("primer_payload_paf_query_binding"));
    }
    let (qs, qe) = (number(f[2])?, number(f[3])?);
    let mapq = number(f[11])?;
    if f[4] == "*" && f[5] == "*" {
        if qs != 0 || qe != 0 || f[6..11].iter().any(|s| *s != "0") || mapq > 255 || f.len() != 12 {
            return Err(bad("primer_payload_paf_no_hit"));
        }
        return Ok(None);
    }
    let (tl, ts, te, matches, block) = (
        number(f[6])?,
        number(f[7])?,
        number(f[8])?,
        number(f[9])?,
        number(f[10])?,
    );
    if !matches!(f[4], "+" | "-")
        || f[5].is_empty()
        || f[5] == "*"
        || qs >= qe
        || qe > query.query_length
        || ts >= te
        || te > tl
        || block == 0
        || matches > block
        || mapq > 255
    {
        return Err(bad("primer_payload_paf_coordinates"));
    }
    let mut tags = BTreeMap::new();
    for tag in &f[12..] {
        let mut parts = tag.splitn(3, ':');
        let (name, ty, value) = (
            parts.next().unwrap_or(""),
            parts.next().unwrap_or(""),
            parts.next().unwrap_or(""),
        );
        if name.len() != 2
            || !name.bytes().all(|c| c.is_ascii_alphanumeric())
            || ty.len() != 1
            || !b"AifZHB".contains(&ty.as_bytes()[0])
            || value.is_empty()
            || tags.insert(name, (ty, value)).is_some()
        {
            return Err(bad("primer_payload_paf_tag"));
        }
    }
    let typed = |name, ty| -> Result<Option<&str>, EngineError> {
        tags.get(name)
            .map(|(t, v)| {
                if *t == ty {
                    Ok(*v)
                } else {
                    Err(bad("primer_payload_paf_tag_type"))
                }
            })
            .transpose()
    };
    let cg = typed("cg", "Z")?.map(|s| cigar(s, b)).transpose()?;
    let cs = typed("cs", "Z")?.map(|s| difference(s, b)).transpose()?;
    if let (Some(a), Some(z)) = (&cg, &cs)
        && (normalized(a)? != normalized(z)? || !detail_agrees(a, z))
    {
        return Err(bad("primer_payload_paf_contradictory_detail"));
    }
    let edits = typed("NM", "i")?.map(number).transpose()?;
    if edits.is_some_and(|n| n > block) {
        return Err(bad("primer_payload_paf_edit_count"));
    }
    let role = typed("tp", "A")?.map(str::to_string);
    if role.as_ref().is_some_and(|s| s.len() != 1) {
        return Err(bad("primer_payload_paf_role"));
    }
    let ops = cs.as_ref().or(cg.as_ref());
    if strict && ops.is_none() {
        return Err(bad("primer_payload_alignment_detail_required"));
    }
    let mut blocks = Vec::new();
    let mut skips = Vec::new();
    if let Some(ops) = ops {
        let (mut q, mut t, mut nm, mut matching, mut columns) =
            (0usize, ts, 0usize, 0usize, 0usize);
        for op in ops {
            let qc = b"M=XI".contains(&op.kind);
            let tc = b"M=XDN".contains(&op.kind);
            let qnext = q
                .checked_add(if qc { op.len } else { 0 })
                .ok_or_else(|| bad("primer_payload_cigar_overflow"))?;
            let tnext = t
                .checked_add(if tc { op.len } else { 0 })
                .ok_or_else(|| bad("primer_payload_cigar_overflow"))?;
            if qnext > qe - qs || tnext > te {
                return Err(bad("primer_payload_cigar_span"));
            }
            if b"M=X".contains(&op.kind) {
                let query = if f[4] == "+" {
                    (qs + q, qs + qnext)
                } else {
                    (qe - qnext, qe - q)
                };
                blocks.push(Block {
                    query,
                    target: (t, tnext),
                });
            }
            if op.kind == b'N' {
                skips.push((t, tnext));
            }
            if b"XID".contains(&op.kind) {
                nm = nm
                    .checked_add(op.len)
                    .ok_or_else(|| bad("primer_payload_cigar_overflow"))?;
            }
            if op.kind == b'=' {
                matching = matching
                    .checked_add(op.len)
                    .ok_or_else(|| bad("primer_payload_cigar_overflow"))?;
            }
            if op.kind != b'N' {
                columns = columns
                    .checked_add(op.len)
                    .ok_or_else(|| bad("primer_payload_cigar_overflow"))?;
            }
            q = qnext;
            t = tnext;
        }
        if q != qe - qs
            || t != te
            || columns != block
            || edits.is_some_and(|n| n != block - matches)
            || ((cs.is_some() || ops.iter().all(|op| op.kind != b'M'))
                && (matching != matches || edits.is_some_and(|n| n != nm)))
        {
            return Err(bad("primer_payload_cigar_span_or_counts"));
        }
    }
    blocks.sort_by_key(|x| x.query);
    let mut merged: Vec<Block> = Vec::new();
    for block in blocks {
        if let Some(last) = merged.last_mut()
            && last.query.1 == block.query.0
            && (if f[4] == "+" {
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
    let blocks = merged;
    let intervals = subtract_masks(
        blocks.iter().map(|x| x.query).collect(),
        &query.excluded_query_intervals,
    );
    Ok(Some((
        Alignment {
            segment: PrimerPayloadSegment {
                query_interval: (qs, qe),
                target_name: f[5].into(),
                target_interval: (ts, te),
                strand: f[4].chars().next().unwrap(),
                mapq: (mapq != 255).then_some(mapq as u8),
                reported_role: role,
                reported_identity_per_mille: ((matches as u128 * 1000) / block as u128) as usize,
                edits,
                alignment_detail_verified: ops.is_some(),
                alignment_operations: ops
                    .into_iter()
                    .flatten()
                    .map(|op| (op.kind as char, op.len))
                    .collect(),
                aligned_query_intervals: intervals,
                gene_ids: Vec::new(),
                transcript_ids: Vec::new(),
                splice_edges: Vec::new(),
            },
            blocks,
            skips,
        },
        tl,
    )))
}

pub(super) fn subtract_masks(
    spans: Vec<(usize, usize)>,
    masks: &[(usize, usize)],
) -> Vec<(usize, usize)> {
    let mut result = Vec::new();
    for (start, end) in union(spans) {
        let mut s = start;
        for &(a, z) in masks {
            if z <= s {
                continue;
            }
            if a >= end {
                break;
            }
            if a > s {
                result.push((s, a.min(end)));
            }
            s = s.max(z);
            if s >= end {
                break;
            }
        }
        if s < end {
            result.push((s, end));
        }
    }
    result
}
