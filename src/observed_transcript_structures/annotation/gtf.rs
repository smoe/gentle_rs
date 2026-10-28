//! Strict, allocation-bounded row grammar; no permissive genome-index fallback.

use super::*;

pub(super) struct Row<'a> {
    pub feature: &'a str,
    pub contig: usize,
    pub strand: char,
    pub span: Interval,
    pub gene: String,
    pub transcript: Option<String>,
    pub gene_version: Option<String>,
    pub transcript_version: Option<String>,
}

fn decimal(value: &str) -> Result<usize, EngineError> {
    if value.is_empty() || !value.bytes().all(|c| c.is_ascii_digit()) {
        return Err(bad("observed_gtf_coordinate"));
    }
    value.parse().map_err(|_| bad("observed_gtf_coordinate"))
}

fn attributes(raw: &str) -> Result<[Option<String>; 4], EngineError> {
    let bytes = raw.as_bytes();
    let mut out: [Option<String>; 4] = [None, None, None, None];
    let mut i = 0;
    while i < bytes.len() {
        while bytes.get(i) == Some(&b' ') {
            i += 1;
        }
        if i == bytes.len() {
            break;
        }
        let start = i;
        if !bytes[i].is_ascii_alphabetic() && bytes[i] != b'_' {
            return Err(bad("observed_gtf_attribute_key"));
        }
        i += 1;
        while i < bytes.len() && (bytes[i].is_ascii_alphanumeric() || bytes[i] == b'_') {
            i += 1;
        }
        let key = &raw[start..i];
        if bytes.get(i) != Some(&b' ') {
            return Err(bad("observed_gtf_attribute_separator"));
        }
        while bytes.get(i) == Some(&b' ') {
            i += 1;
        }
        let slot = match key {
            "gene_id" => Some(0),
            "transcript_id" => Some(1),
            "gene_version" => Some(2),
            "transcript_version" => Some(3),
            _ => None,
        };
        let quoted = bytes.get(i) == Some(&b'"');
        if slot.is_some() && !quoted {
            return Err(bad("observed_gtf_attribute_quoted"));
        }
        let mut value = String::new();
        let mut value_bytes = 0usize;
        if quoted {
            i += 1;
            let mut closed = false;
            while i < bytes.len() {
                if bytes[i] == b'"' {
                    i += 1;
                    closed = true;
                    break;
                }
                let ch = if bytes[i] == b'\\' {
                    i += 1;
                    match bytes.get(i) {
                        Some(b'"') => '"',
                        Some(b'\\') => '\\',
                        _ => return Err(bad("observed_gtf_attribute_escape")),
                    }
                } else {
                    raw[i..].chars().next().unwrap()
                };
                i += ch.len_utf8();
                value_bytes += ch.len_utf8();
                if slot.is_some() {
                    if value_bytes > 512 {
                        return Err(bad("observed_gtf_attribute_length"));
                    }
                    value.push(ch);
                }
            }
            if !closed {
                return Err(bad("observed_gtf_attribute_quote"));
            }
        } else {
            let start = i;
            while i < bytes.len() && !matches!(bytes[i], b' ' | b';') {
                if matches!(bytes[i], b'"' | b'\\') {
                    return Err(bad("observed_gtf_attribute_value"));
                }
                i += 1;
            }
            value_bytes = i - start;
        }
        if value_bytes == 0 {
            return Err(bad("observed_gtf_attribute_empty"));
        }
        while bytes.get(i) == Some(&b' ') {
            i += 1;
        }
        if bytes.get(i) != Some(&b';') {
            return Err(bad("observed_gtf_attribute_terminator"));
        }
        i += 1;
        if let Some(slot) = slot {
            if out[slot].is_some() {
                return Err(bad("observed_gtf_attribute_duplicate"));
            }
            if !identifier(&value, 512) {
                return Err(bad("observed_gtf_identifier"));
            }
            out[slot] = Some(value);
        }
    }
    Ok(out)
}

pub(super) fn parse<'a>(
    line: &'a str,
    dictionary: &BTreeMap<&str, (usize, usize)>,
) -> Result<Row<'a>, EngineError> {
    if line.trim().is_empty() {
        return Err(bad("observed_gtf_blank_line"));
    }
    let mut fields = line.split('\t');
    let mut cols = [""; 9];
    for c in &mut cols {
        *c = fields.next().ok_or_else(|| bad("observed_gtf_columns"))?;
    }
    if fields.next().is_some() {
        return Err(bad("observed_gtf_columns"));
    }
    let &(contig, length) = dictionary
        .get(cols[0])
        .ok_or_else(|| bad("observed_gtf_contig"))?;
    if !identifier(cols[1], 4096) || !identifier(cols[2], 4096) {
        return Err(bad("observed_gtf_columns"));
    }
    let (start, end) = (decimal(cols[3])?, decimal(cols[4])?);
    if start == 0 || start > end || end > length {
        return Err(bad("observed_gtf_coordinate"));
    }
    if cols[5] != "."
        && (cols[5].trim() != cols[5] || cols[5].parse::<f64>().ok().is_none_or(|n| !n.is_finite()))
    {
        return Err(bad("observed_gtf_score"));
    }
    if !matches!(cols[7], "." | "0" | "1" | "2") {
        return Err(bad("observed_gtf_frame"));
    }
    let strand = match cols[6] {
        "+" => '+',
        "-" => '-',
        "." => '.',
        _ => return Err(bad("observed_gtf_strand")),
    };
    let [gene, transcript, gene_version, transcript_version] = attributes(cols[8])?;
    let gene = gene.ok_or_else(|| bad("observed_gtf_gene_id"))?;
    if (cols[2] == "gene") != transcript.is_none()
        || (cols[2] == "gene" && transcript_version.is_some())
    {
        return Err(bad("observed_gtf_transcript_id"));
    }
    Ok(Row {
        feature: cols[2],
        contig,
        strand,
        span: (start - 1, end),
        gene,
        transcript,
        gene_version,
        transcript_version,
    })
}
