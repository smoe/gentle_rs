//! Private strict SAM-text decoding over caller-supplied raw lines, not verified file ingress.
//!
//! SAM 1.6 b5341fb / optional tags a71164e. No QNAME grouping, source authority,
//! search-completeness inference, sequence/quality retention, operations or publication.
//! Memory charges bound retained representation and parser scratch, not allocator RSS.

use super::*;
use ring::digest::{Context, SHA256};

#[cfg(test)]
#[path = "sam_text_tests.rs"]
mod tests;

#[derive(Debug, Clone)]
pub(super) struct Limits {
    pub max_input_bytes: usize,
    pub max_line_bytes: usize,
    pub max_lines: usize,
    pub max_rows: usize,
    pub max_dictionary_entries: usize,
    pub max_header_ids: usize,
    pub max_fields_per_line: usize,
    /// Includes actual/SA CIGAR operators and every B-array element, even ignored tags.
    pub max_operations: usize,
    pub max_sa_entries: usize,
    pub max_retained_bytes: usize,
    pub max_work: usize,
}

/// A tag-declared companion, never an invented alignment row or completeness proof.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) struct Companion {
    pub target: String,
    pub position_1based: usize,
    pub strand: char,
    pub cigar: Vec<(char, usize)>,
    pub mapq: u8,
}

#[derive(Debug, Clone)]
pub(super) struct DecodedRow {
    pub row: sam::Row,
    pub read_group: Option<String>,
    pub companions: Vec<Companion>,
}

#[derive(Debug, Default, PartialEq, Eq)]
pub(super) struct Counts {
    pub lines: usize,
    pub header_lines: usize,
    pub alignment_rows: usize,
    pub unmapped: usize,
    pub secondary: usize,
    pub supplementary: usize,
    pub qc_fail: usize,
    pub duplicate: usize,
}

#[derive(Debug)]
pub(super) struct Document {
    pub grammar: &'static str,
    pub raw_sha256: String,
    pub raw_bytes: usize,
    pub dictionary_order: Vec<String>,
    pub read_groups: BTreeSet<String>,
    pub rows: Vec<DecodedRow>,
    /// Flag counts overlap; they are row observations, not query-unit support.
    pub counts: Counts,
    pub operations: usize,
    pub sa_entries: usize,
    pub work: usize,
    pub retained_peak_bytes: usize,
}

fn fail(field: &str) -> EngineError {
    bad(&format!("observed_sam_text_{field}"))
}
fn plus(a: usize, b: usize) -> Result<usize, EngineError> {
    a.checked_add(b).ok_or_else(|| fail("overflow"))
}
fn unsigned(s: &str, max: u64) -> Result<u64, EngineError> {
    if s.is_empty() || !s.bytes().all(|c| c.is_ascii_digit()) {
        return Err(fail("integer"));
    }
    s.parse::<u64>()
        .ok()
        .filter(|n| *n <= max)
        .ok_or_else(|| fail("integer"))
}
fn signed(s: &str, min: i64, max: i64) -> Result<i64, EngineError> {
    let digits = s.strip_prefix(['+', '-']).unwrap_or(s);
    if digits.is_empty() || !digits.bytes().all(|c| c.is_ascii_digit()) {
        return Err(fail("integer"));
    }
    s.parse::<i64>()
        .ok()
        .filter(|n| (min..=max).contains(n))
        .ok_or_else(|| fail("integer"))
}
fn float(s: &str) -> bool {
    let s = s.strip_prefix(['+', '-']).unwrap_or(s);
    let mut exponent = s.split(['e', 'E']);
    let mantissa = exponent.next().unwrap_or("");
    if let Some(e) = exponent.next() {
        let e = e.strip_prefix(['+', '-']).unwrap_or(e);
        if e.is_empty() || !e.bytes().all(|c| c.is_ascii_digit()) {
            return false;
        }
    }
    if exponent.next().is_some() {
        return false;
    }
    let valid = if let Some((a, b)) = mantissa.split_once('.') {
        a.bytes().all(|c| c.is_ascii_digit())
            && !b.is_empty()
            && b.bytes().all(|c| c.is_ascii_digit())
    } else {
        !mantissa.is_empty() && mantissa.bytes().all(|c| c.is_ascii_digit())
    };
    valid && s.parse::<f32>().is_ok_and(f32::is_finite)
}
fn reference_name(s: &str) -> bool {
    !s.is_empty()
        && s.len() <= 4096
        && !s.starts_with(['*', '='])
        && s.bytes()
            .all(|c| c.is_ascii_alphanumeric() || b"!#$%&*+./:;=?@^_|~-".contains(&c))
}
fn identifier(s: &str) -> bool {
    !s.is_empty() && s.len() <= 4096 && s.bytes().all(|c| (33..=126).contains(&c))
}
fn tag(s: &str) -> bool {
    s.len() == 2 && s.as_bytes()[0].is_ascii_alphabetic() && s.as_bytes()[1].is_ascii_alphanumeric()
}

struct Meter<'a> {
    limits: &'a Limits,
    retained: usize,
    peak: usize,
    work: usize,
    operations: usize,
    sa: usize,
    active_scratch: usize,
}
impl Meter<'_> {
    fn work(&mut self, n: usize) -> Result<(), EngineError> {
        self.work = plus(self.work, n)?;
        if self.work > self.limits.max_work {
            return Err(fail("work_budget"));
        }
        Ok(())
    }
    fn scratch(&mut self, n: usize) -> Result<(), EngineError> {
        let n = plus(self.retained, n)?;
        if n > self.limits.max_retained_bytes {
            return Err(fail("memory_budget"));
        }
        self.peak = self.peak.max(n);
        Ok(())
    }
    fn retain(&mut self, n: usize) -> Result<(), EngineError> {
        self.scratch(plus(n, self.active_scratch)?)?;
        self.retained = plus(self.retained, n)?;
        Ok(())
    }
    fn operation(&mut self) -> Result<(), EngineError> {
        self.operations = plus(self.operations, 1)?;
        self.work(1)?;
        if self.operations > self.limits.max_operations {
            return Err(fail("operation_budget"));
        }
        Ok(())
    }
}

fn cigar(
    s: &str,
    m: &mut Meter<'_>,
    check: &mut dyn FnMut() -> Result<(), EngineError>,
) -> Result<Option<Vec<(char, usize)>>, EngineError> {
    if s == "*" {
        return Ok(None);
    }
    if s.is_empty() {
        return Err(fail("cigar"));
    }
    let (mut out, mut start) = (Vec::new(), 0);
    for (i, c) in s.bytes().enumerate() {
        if c.is_ascii_digit() {
            continue;
        }
        check()?;
        if !b"MIDNSHP=X".contains(&c) {
            return Err(fail("cigar"));
        }
        let n = usize::try_from(unsigned(&s[start..i], usize::MAX as u64)?)
            .map_err(|_| fail("integer"))?;
        if n == 0 {
            return Err(fail("cigar"));
        }
        m.operation()?;
        m.retain(2 * std::mem::size_of::<(char, usize)>())?;
        out.push((c as char, n));
        start = i + 1;
    }
    if start != s.len() {
        return Err(fail("cigar"));
    }
    Ok(Some(out))
}

fn optional(
    kind: &str,
    value: &str,
    m: &mut Meter<'_>,
    check: &mut dyn FnMut() -> Result<(), EngineError>,
) -> Result<(), EngineError> {
    match kind {
        "A" if value.len() == 1 && value.bytes().all(|c| (33..=126).contains(&c)) => (),
        "Z" if value.bytes().all(|c| (32..=126).contains(&c)) => (),
        "H" if value.len() % 2 == 0
            && value
                .bytes()
                .all(|c| c.is_ascii_digit() || (b'A'..=b'F').contains(&c)) =>
        {
            ()
        }
        "i" => {
            signed(value, i32::MIN as i64, u32::MAX as i64)?;
        }
        "f" if float(value) => (),
        "B" => {
            let mut values = value.split(',');
            let subtype = values.next().unwrap_or("");
            if !matches!(subtype, "c" | "C" | "s" | "S" | "i" | "I" | "f") {
                return Err(fail("optional_type"));
            }
            for v in values {
                check()?;
                m.operation()?;
                match subtype {
                    "c" => {
                        signed(v, i8::MIN as i64, i8::MAX as i64)?;
                    }
                    "C" => {
                        signed(v, 0, u8::MAX as i64)?;
                    }
                    "s" => {
                        signed(v, i16::MIN as i64, i16::MAX as i64)?;
                    }
                    "S" => {
                        signed(v, 0, u16::MAX as i64)?;
                    }
                    "i" => {
                        signed(v, i32::MIN as i64, i32::MAX as i64)?;
                    }
                    "I" => {
                        signed(v, 0, u32::MAX as i64)?;
                    }
                    "f" if float(v) => (),
                    _ => return Err(fail("optional_value")),
                }
            }
        }
        _ => return Err(fail("optional_value")),
    }
    Ok(())
}

/// Incremental raw framing and atomic finish; failed/interrupted pushes poison the consumer.
pub(super) struct Consumer<'a> {
    expected: &'a BTreeMap<String, usize>,
    meter: Meter<'a>,
    sha: Context,
    bytes: usize,
    counts: Counts,
    seen_dictionary: BTreeSet<String>,
    dictionary_order: Vec<String>,
    groups: BTreeSet<String>,
    programs: BTreeSet<String>,
    rows: Vec<DecodedRow>,
    unterminated: bool,
    poisoned: bool,
}

impl<'a> Consumer<'a> {
    pub(super) fn new(
        expected: &'a BTreeMap<String, usize>,
        limits: &'a Limits,
        check: &mut dyn FnMut() -> Result<(), EngineError>,
    ) -> Result<Self, EngineError> {
        check()?;
        if [
            limits.max_input_bytes,
            limits.max_line_bytes,
            limits.max_lines,
            limits.max_rows,
            limits.max_dictionary_entries,
            limits.max_operations,
            limits.max_retained_bytes,
            limits.max_work,
        ]
        .contains(&0)
            || limits.max_input_bytes > 128 * 1024 * 1024
            || limits.max_line_bytes > 1024 * 1024
            || limits.max_lines > 1_000_000
            || limits.max_rows > 65_536
            || limits.max_dictionary_entries > 100_000
            || limits.max_header_ids > 100_000
            || !(11..=1035).contains(&limits.max_fields_per_line)
            || limits.max_operations > 1_000_000
            || limits.max_sa_entries > 65_536
            || limits.max_retained_bytes > 128 * 1024 * 1024
            || limits.max_work > 1_000_000_000
            || expected.is_empty()
            || expected.len() > limits.max_dictionary_entries
        {
            return Err(fail("admission"));
        }
        let mut meter = Meter {
            limits,
            retained: 0,
            peak: 0,
            work: 0,
            operations: 0,
            sa: 0,
            active_scratch: 0,
        };
        meter.retain(2048)?;
        for (name, &length) in expected {
            check()?;
            meter.work(plus(name.len(), 1)?)?;
            if !reference_name(name) || length == 0 || length > i32::MAX as usize {
                return Err(fail("dictionary"));
            }
        }
        Ok(Self {
            expected,
            meter,
            sha: Context::new(&SHA256),
            bytes: 0,
            counts: Counts::default(),
            seen_dictionary: BTreeSet::new(),
            dictionary_order: Vec::new(),
            groups: BTreeSet::new(),
            programs: BTreeSet::new(),
            rows: Vec::new(),
            unterminated: false,
            poisoned: false,
        })
    }

    pub(super) fn push_line(
        &mut self,
        raw: &[u8],
        check: &mut dyn FnMut() -> Result<(), EngineError>,
    ) -> Result<(), EngineError> {
        if self.poisoned {
            return Err(fail("poisoned"));
        }
        let result = self.push(raw, check);
        self.meter.active_scratch = 0;
        self.poisoned |= result.is_err();
        result
    }

    fn push(
        &mut self,
        raw: &[u8],
        check: &mut dyn FnMut() -> Result<(), EngineError>,
    ) -> Result<(), EngineError> {
        check()?;
        if self.unterminated {
            return Err(fail("line_after_unterminated"));
        }
        self.bytes = plus(self.bytes, raw.len())?;
        self.counts.lines = plus(self.counts.lines, 1)?;
        let b = self.meter.limits;
        if self.bytes > b.max_input_bytes
            || raw.len() > b.max_line_bytes
            || self.counts.lines > b.max_lines
        {
            return Err(fail("input_budget"));
        }
        self.meter.work(plus(raw.len(), 1)?)?;
        let scratch = raw
            .len()
            .checked_mul(2)
            .and_then(|n| n.checked_add(b.max_fields_per_line * 512 + 8192))
            .ok_or_else(|| fail("overflow"))?;
        self.meter.scratch(scratch)?;
        self.meter.active_scratch = scratch;
        self.sha.update(raw);
        self.unterminated = !raw.ends_with(b"\n");
        let raw = raw.strip_suffix(b"\n").unwrap_or(raw);
        let raw = raw.strip_suffix(b"\r").unwrap_or(raw);
        let body = std::str::from_utf8(raw).map_err(|_| fail("utf8"))?;
        if body.is_empty() || body.chars().any(|c| c.is_control() && c != '\t') {
            return Err(fail("line"));
        }
        if body.starts_with("@CO\t") {
            if !self.rows.is_empty() {
                return Err(fail("header_after_alignment"));
            }
            self.counts.header_lines += 1;
            return Ok(());
        }
        let fields = body
            .split('\t')
            .take(b.max_fields_per_line + 1)
            .collect::<Vec<_>>();
        if fields.len() > b.max_fields_per_line {
            return Err(fail("field_budget"));
        }
        if body.starts_with('@') {
            self.header(&fields)?;
        } else {
            self.alignment(&fields, check)?;
        }
        // Recheck after retention as well: a long row must not hide its parser scratch.
        self.meter.scratch(scratch)?;
        Ok(())
    }

    fn header(&mut self, f: &[&str]) -> Result<(), EngineError> {
        if !self.rows.is_empty() {
            return Err(fail("header_after_alignment"));
        }
        if !matches!(f[0], "@HD" | "@SQ" | "@RG" | "@PG") {
            return Err(fail("header_kind"));
        }
        let mut tags = BTreeMap::new();
        for field in &f[1..] {
            let (key, value) = field.split_once(':').ok_or_else(|| fail("header_field"))?;
            if !tag(key)
                || value.is_empty()
                || tags.insert(key, value).is_some()
                || (!(key == "DS" && matches!(f[0], "@SQ" | "@RG" | "@PG"))
                    && !(f[0] == "@PG" && key == "CL")
                    && !value.bytes().all(|c| (32..=126).contains(&c)))
            {
                return Err(fail("header_field"));
            }
        }
        let required = |key| {
            tags.get(key)
                .copied()
                .ok_or_else(|| fail("header_required"))
        };
        match f[0] {
            "@HD" => {
                if self.counts.lines != 1 || required("VN")? != "1.6" {
                    return Err(fail("header_version_or_order"));
                }
            }
            "@SQ" => {
                let name = required("SN")?;
                let length = unsigned(required("LN")?, i32::MAX as u64)? as usize;
                if !reference_name(name)
                    || length == 0
                    || self.expected.get(name) != Some(&length)
                    || self.seen_dictionary.contains(name)
                    || tags.get("TP").is_some_and(|v| *v != "linear")
                {
                    return Err(fail("dictionary"));
                }
                self.meter.retain(plus(1024, name.len() * 2)?)?;
                self.seen_dictionary.insert(name.into());
                self.dictionary_order.push(name.into());
            }
            "@RG" | "@PG" => {
                let id = required("ID")?;
                let total = plus(self.groups.len(), self.programs.len())?;
                if !identifier(id) || total >= self.meter.limits.max_header_ids {
                    return Err(fail("header_id_budget_or_value"));
                }
                let set = if f[0] == "@RG" {
                    &mut self.groups
                } else {
                    &mut self.programs
                };
                if set.contains(id) {
                    return Err(fail("duplicate_header_id"));
                }
                self.meter.retain(plus(512, id.len())?)?;
                set.insert(id.into());
            }
            _ => unreachable!(),
        }
        self.counts.header_lines += 1;
        Ok(())
    }

    fn alignment(
        &mut self,
        f: &[&str],
        check: &mut dyn FnMut() -> Result<(), EngineError>,
    ) -> Result<(), EngineError> {
        if f.len() < 11 {
            return Err(fail("mandatory_fields"));
        }
        if self.seen_dictionary.len() != self.expected.len() {
            return Err(fail("dictionary_incomplete"));
        }
        if self.rows.len() >= self.meter.limits.max_rows {
            return Err(fail("row_budget"));
        }
        if f[0].is_empty()
            || f[0].len() > 254
            || !f[0].bytes().all(|c| (33..=126).contains(&c) && c != b'@')
        {
            return Err(fail("qname"));
        }
        self.meter
            .retain(plus(512, plus(f[0].len(), f[2].len())?)?)?;
        let flag = unsigned(f[1], 0xfff)? as u16;
        let target = if f[2] == "*" {
            None
        } else {
            if !reference_name(f[2]) || !self.expected.contains_key(f[2]) {
                return Err(fail("target"));
            }
            Some(f[2].to_string())
        };
        let position_1based = unsigned(f[3], i32::MAX as u64)? as usize;
        let mapq = unsigned(f[4], 255)? as u8;
        let next = match f[6] {
            "*" => None,
            "=" => target.as_deref(),
            v if reference_name(v) && self.expected.contains_key(v) => Some(v),
            _ => return Err(fail("mate_reference")),
        };
        let pnext = unsigned(f[7], i32::MAX as u64)? as usize;
        if next.is_some_and(|name| pnext > self.expected[name]) {
            return Err(fail("mate_position"));
        }
        signed(f[8], -(i32::MAX as i64), i32::MAX as i64)?;
        let sequence_length = if f[9] == "*" {
            None
        } else {
            if f[9].is_empty()
                || !f[9]
                    .bytes()
                    .all(|c| c.is_ascii_alphabetic() || b"=.".contains(&c))
            {
                return Err(fail("sequence"));
            }
            Some(f[9].len())
        };
        let quality_length = if f[10] == "*" {
            None
        } else {
            if f[10].is_empty() || !f[10].bytes().all(|c| (33..=126).contains(&c)) {
                return Err(fail("quality"));
            }
            Some(f[10].len())
        };
        let mut row = sam::Row {
            ordinal: self.rows.len() as u64 + 1,
            qname: f[0].into(),
            flag,
            target,
            position_1based,
            mapq,
            cigar: cigar(f[5], &mut self.meter, check)?,
            sequence_length,
            quality_length,
            has_cg: false,
        };
        let (mut rg, mut companions, mut tags) = (None, Vec::new(), BTreeSet::new());
        for field in &f[11..] {
            check()?;
            let mut parts = field.splitn(3, ':');
            let key = parts.next().unwrap_or("");
            let kind = parts.next().ok_or_else(|| fail("optional_field"))?;
            let value = parts.next().ok_or_else(|| fail("optional_field"))?;
            if !tag(key) || !tags.insert(key) {
                return Err(fail("optional_tag"));
            }
            optional(kind, value, &mut self.meter, check)?;
            match key {
                "RG" => {
                    if kind != "Z" || !self.groups.contains(value) {
                        return Err(fail("read_group"));
                    }
                    self.meter.retain(plus(64, value.len())?)?;
                    rg = Some(value.into());
                }
                "CG" => {
                    if kind != "B" || value.split(',').next() != Some("I") {
                        return Err(fail("cg_type"));
                    }
                    row.has_cg = true;
                }
                "SA" => {
                    if kind != "Z" {
                        return Err(fail("sa_type"));
                    }
                    let entries = value
                        .strip_suffix(';')
                        .filter(|s| !s.is_empty())
                        .ok_or_else(|| fail("sa_value"))?;
                    for entry in entries.split(';') {
                        check()?;
                        self.meter.sa = plus(self.meter.sa, 1)?;
                        if self.meter.sa > self.meter.limits.max_sa_entries {
                            return Err(fail("sa_budget"));
                        }
                        let mut fields = entry.split(',');
                        let mut next = || fields.next().ok_or_else(|| fail("sa_value"));
                        let target = next()?;
                        if !reference_name(target) || !self.expected.contains_key(target) {
                            return Err(fail("sa_target"));
                        }
                        let position_1based = unsigned(next()?, i32::MAX as u64)? as usize;
                        let strand = match next()? {
                            "+" => '+',
                            "-" => '-',
                            _ => return Err(fail("sa_strand")),
                        };
                        let cigar = cigar(next()?, &mut self.meter, check)?
                            .ok_or_else(|| fail("sa_cigar"))?;
                        let mapq = unsigned(next()?, 255)? as u8;
                        unsigned(next()?, u32::MAX as u64)?; // Syntax only: no retained NM/identity score.
                        if fields.next().is_some() {
                            return Err(fail("sa_value"));
                        }
                        self.meter.retain(plus(512, target.len())?)?;
                        let mut candidate = sam::Row {
                            ordinal: row.ordinal,
                            qname: row.qname.clone(),
                            flag: if strand == '-' { 16 } else { 0 },
                            target: Some(target.into()),
                            position_1based,
                            mapq,
                            cigar: Some(cigar),
                            sequence_length: None,
                            quality_length: None,
                            has_cg: false,
                        };
                        self.validate_geometry(&candidate)?;
                        let cigar = candidate.cigar.take().unwrap();
                        companions.push(Companion {
                            target: target.into(),
                            position_1based,
                            strand,
                            cigar,
                            mapq,
                        });
                    }
                }
                _ => (), // Validate syntax, discard all unused tag values.
            }
        }
        self.validate_geometry(&row)?;
        self.counts.alignment_rows += 1;
        for (bit, count) in [
            (4, &mut self.counts.unmapped),
            (256, &mut self.counts.secondary),
            (2048, &mut self.counts.supplementary),
            (512, &mut self.counts.qc_fail),
            (1024, &mut self.counts.duplicate),
        ] {
            *count += usize::from(flag & bit != 0);
        }
        self.rows.push(DecodedRow {
            row,
            read_group: rg,
            companions,
        });
        Ok(())
    }

    pub(super) fn finish(
        self,
        check: &mut dyn FnMut() -> Result<(), EngineError>,
    ) -> Result<Document, EngineError> {
        check()?;
        if self.poisoned {
            return Err(fail("poisoned"));
        }
        if self.seen_dictionary.len() != self.expected.len() {
            return Err(fail("dictionary_incomplete"));
        }
        let hex = self
            .sha
            .finish()
            .as_ref()
            .iter()
            .map(|b| format!("{b:02x}"))
            .collect::<String>();
        Ok(Document {
            grammar: "sam_text_strict.private.v1",
            raw_sha256: format!("sha256:{hex}"),
            raw_bytes: self.bytes,
            dictionary_order: self.dictionary_order,
            read_groups: self.groups,
            rows: self.rows,
            counts: self.counts,
            operations: self.meter.operations,
            sa_entries: self.meter.sa,
            work: self.meter.work,
            retained_peak_bytes: self.meter.peak,
        })
    }

    fn validate_geometry(&mut self, row: &sam::Row) -> Result<(), EngineError> {
        let temporary = row
            .cigar
            .as_ref()
            .map_or(0, Vec::len)
            .checked_mul(std::mem::size_of::<(char, usize)>())
            .ok_or_else(|| fail("overflow"))?;
        self.meter
            .scratch(plus(self.meter.active_scratch, temporary)?)?;
        sam::validate_decoded_row(row, self.expected)
    }
}
