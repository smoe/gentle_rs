//! Pure decoded-SAM admission, clipping normalization and declared-unit identity.
//!
//! No text/file parser, reference discovery, PAF identity fabrication or public
//! schema. Reference hashes/dictionary and search coverage remain declarations.

use super::*;

#[derive(Debug, Clone, Serialize, PartialEq, Eq, PartialOrd, Ord)]
pub(super) struct UnitKey {
    pub source_sha256: String,
    pub unit_id: String,
}

/// Fields consumed after separate text/tag validation; no raw sequences or unit grouping.
#[derive(Debug, Clone, Serialize)]
pub(super) struct Row {
    pub ordinal: u64,
    pub qname: String,
    pub flag: u16,
    pub target: Option<String>,
    pub position_1based: usize,
    pub mapq: u8,
    pub cigar: Option<Vec<(char, usize)>>,
    pub sequence_length: Option<usize>,
    pub quality_length: Option<usize>,
    pub has_cg: bool,
}

#[derive(Debug, Clone, Serialize)]
pub(super) struct Unit {
    pub unit_id: String,
    pub declared_axis_length: Option<usize>,
    pub masks: Option<Vec<Interval>>,
    pub ends: [Vec<EndConstraint>; 2],
    pub search_declared_complete: bool,
    pub views_declared_complete: bool,
    pub rows: Vec<Row>,
}

#[derive(Debug, Clone, Serialize)]
pub(super) struct SamInput {
    pub source_sha256: String,
    pub reference: ReferenceBinding,
    pub reference_gene: ReferenceGene,
    pub comparison_genes: Vec<ReferenceGene>,
    pub annotation_sources: Vec<AnnotationSource>,
    pub dictionary: BTreeMap<String, usize>,
    pub unpadded_reference: bool,
    pub units: Vec<Unit>,
    pub required_flank_bases: usize,
    pub budgets: Budgets,
    pub max_rows: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Serialize)]
pub(super) enum Unavailable {
    NotProvided,
    Unmapped,
    Paired,
    PaddedReference,
    CigarNotProvided,
    NoSupportedBases,
    CgOrPlaceholder,
    QueryAxisUnsupported,
}

#[derive(Debug, PartialEq, Eq)]
pub(super) struct RowDiagnostic {
    pub ordinal: u64,
    pub unavailable: Option<Unavailable>,
    pub axis_length: Option<usize>,
    pub clipping: Option<core::Clipping>,
}

#[derive(Debug)]
pub(super) struct UnitDiagnostic {
    pub unit: UnitKey,
    pub axis_length: Option<usize>,
    pub unavailable: BTreeSet<Unavailable>,
    pub rows: Vec<RowDiagnostic>,
}

#[derive(Debug)]
pub(super) struct Analysis {
    pub coordinate_authority: &'static str,
    pub identity_scope: &'static str,
    pub summary: core::CoreSummary<UnitKey>,
    pub units: Vec<UnitDiagnostic>,
}

struct Normalized {
    axis_length: usize,
    query: Interval,
    target: Interval,
    operations: Vec<(char, usize)>,
    clipping: core::Clipping,
}

fn add(a: usize, b: usize) -> Result<usize, EngineError> {
    a.checked_add(b).ok_or_else(|| bad("observed_sam_overflow"))
}

fn name(s: &str) -> bool {
    !s.is_empty() && s.len() <= 4096 && s != "*" && !s.chars().any(char::is_whitespace)
}

fn normalize(row: &Row) -> Result<Option<Normalized>, EngineError> {
    let Some(ops) = &row.cigar else {
        if row.quality_length.is_some() && row.quality_length != row.sequence_length {
            return Err(bad("observed_sam_sequence_quality"));
        }
        return Ok(None);
    };
    if ops.is_empty()
        || ops
            .iter()
            .any(|(c, n)| *n == 0 || !"MIDNSHP=X".contains(*c))
    {
        return Err(bad("observed_sam_cigar"));
    }
    let (mut a, mut z) = (0, ops.len());
    let mut clipping = core::Clipping {
        hard: [0; 2],
        soft: [0; 2],
        padding_operations: 0,
    };
    if ops[a].0 == 'H' {
        clipping.hard[0] = ops[a].1;
        a += 1;
    }
    if a < z && ops[a].0 == 'S' {
        clipping.soft[0] = ops[a].1;
        a += 1;
    }
    if a < z && ops[z - 1].0 == 'H' {
        z -= 1;
        clipping.hard[1] = ops[z].1;
    }
    if a < z && ops[z - 1].0 == 'S' {
        z -= 1;
        clipping.soft[1] = ops[z].1;
    }
    if ops[a..z].iter().any(|(c, _)| matches!(c, 'H' | 'S')) {
        return Err(bad("observed_sam_clipping"));
    }
    let (mut q, mut t) = (0, 0);
    let mut operations = Vec::with_capacity(z - a);
    for &(c, n) in &ops[a..z] {
        if c == 'P' {
            clipping.padding_operations += 1;
            continue;
        }
        if "M=XI".contains(c) {
            q = add(q, n)?;
        }
        if "M=XDN".contains(c) {
            t = add(t, n)?;
        }
        operations.push((c, n));
    }
    let before = add(clipping.hard[0], clipping.soft[0])?;
    let after = add(clipping.soft[1], clipping.hard[1])?;
    let axis_length = add(add(before, q)?, after)?;
    let stored_length = add(add(clipping.soft[0], q)?, clipping.soft[1])?;
    if row.sequence_length.is_some_and(|n| n != stored_length)
        || (row.quality_length.is_some() && row.quality_length != row.sequence_length)
    {
        return Err(bad("observed_sam_sequence_quality"));
    }
    let query = if row.flag & 0x10 == 0 {
        (before, add(before, q)?)
    } else {
        (after, add(after, q)?)
    };
    let start = row.position_1based.saturating_sub(1);
    Ok(Some(Normalized {
        axis_length,
        query,
        target: (start, add(start, t)?),
        operations,
        clipping,
    }))
}

struct SamIdentity<'a>(&'a SamInput);

fn validate_target(
    row: &Row,
    geometry: Option<&Normalized>,
    dictionary: &BTreeMap<String, usize>,
) -> Result<(), EngineError> {
    if row.flag & 4 == 0 {
        let length = row
            .target
            .as_ref()
            .and_then(|t| dictionary.get(t))
            .ok_or_else(|| bad("observed_sam_target"))?;
        if row.position_1based == 0
            || row.position_1based > *length
            || geometry.is_some_and(|g| g.target.1 > *length)
        {
            return Err(bad("observed_sam_target"));
        }
    }
    Ok(())
}

/// Same geometry/target checks for text intake, without unit identity or interpretation.
pub(super) fn validate_decoded_row(
    row: &Row,
    dictionary: &BTreeMap<String, usize>,
) -> Result<(), EngineError> {
    let normalized = normalize(row)?;
    validate_target(row, normalized.as_ref(), dictionary)
}

impl core::Identity<UnitKey> for SamIdentity<'_> {
    fn context_id(&self) -> Result<String, EngineError> {
        digest(&("gentle.observed_sam_context.private.v1", self.0))
    }

    fn observation_id(
        &self,
        context: &core::Context<'_>,
        query: &core::Query<'_, UnitKey>,
        segment: &core::Segment<'_>,
        observation: &core::CoreObservation<UnitKey>,
    ) -> Result<String, EngineError> {
        digest(&(
            "gentle.observed_sam_observation.private.v1",
            context.reference,
            query.key,
            segment.source_ordinal,
            query.alignment_sha256,
            query.query_length,
            segment.query_interval,
            segment.target_name,
            segment.target_interval,
            segment.strand,
            segment.alignment_operations,
            segment.clipping,
            &observation.excluded_query_intervals,
            &observation.ends,
        ))
    }
}

/// Interpret only decoded, explicitly grouped units. No QNAME-derived identity.
pub(super) fn analyze(
    input: &SamInput,
    check: &mut dyn FnMut() -> Result<(), EngineError>,
) -> Result<Analysis, EngineError> {
    check()?;
    let b = &input.budgets;
    b.validate(input.required_flank_bases)?;
    if !is_hash(&input.source_sha256)
        || !is_hash(&input.reference.fasta_sha256)
        || input.reference.assembly.trim().is_empty()
        || input.units.len() > b.max_queries.min(b.max_records).min(4096)
        || !(1..=65_536).contains(&input.max_rows)
        || input.dictionary.is_empty()
        || input.dictionary.len() > 100_000
        || input.annotation_sources.len() > b.max_annotation_rows
        || input.comparison_genes.len() > b.max_annotation_rows
    {
        return Err(bad("observed_sam_admission"));
    }
    let context = core::Context {
        reference: &input.reference,
        reference_gene: &input.reference_gene,
        comparison_genes: &input.comparison_genes,
        annotation_sources: &input.annotation_sources,
        required_flank_bases: input.required_flank_bases,
        budgets: b,
    };
    let mut work = 0;
    for (target, &length) in &input.dictionary {
        check()?;
        dag::charge(&mut work, b.max_total_work)?;
        if !name(target) || length == 0 {
            return Err(bad("observed_sam_dictionary"));
        }
    }
    let mut annotation_rows = 0;
    for source in &input.annotation_sources {
        for t in &source.annotation.transcripts {
            check()?;
            annotation_rows = add(annotation_rows, 1)?;
            if annotation_rows > b.max_annotation_rows
                || input.dictionary.get(&t.target_name) != Some(&t.target_length)
            {
                return Err(bad("observed_sam_annotation_dictionary"));
            }
        }
    }
    comparison::sources(&context)?;
    // Count every supplied row/operation before copying normalized geometry,
    // including rows belonging to units that will be unavailable.
    let (mut rows, mut operations) = (0, 0);
    let mut ids = BTreeSet::new();
    let mut ordinals = BTreeSet::new();
    for unit in &input.units {
        check()?;
        dag::charge(&mut work, b.max_total_work)?;
        rows = add(rows, unit.rows.len())?;
        operations = add(operations, unit.masks.as_ref().map_or(0, Vec::len))?;
        if !name(&unit.unit_id)
            || !ids.insert(&unit.unit_id)
            || rows > input.max_rows
            || unit.rows.len() > b.max_segments_per_query
            || unit.declared_axis_length == Some(0)
            || unit.ends.iter().any(|e| e.len() > 5)
        {
            return Err(bad("observed_sam_unit_binding"));
        }
        let mut qname = None;
        for row in &unit.rows {
            check()?;
            dag::charge(&mut work, b.max_total_work)?;
            let n = row.cigar.as_ref().map_or(0, Vec::len);
            operations = add(operations, n)?;
            dag::charge_many(&mut work, b.max_total_work, n)?;
            if operations > b.max_operations {
                return Err(bad("observed_sam_operation_budget"));
            }
            if row.ordinal == 0
                || !ordinals.insert(row.ordinal)
                || !name(&row.qname)
                || qname.replace(&row.qname).is_some_and(|q| q != &row.qname)
                || row.flag & !0xfff != 0
            {
                return Err(bad("observed_sam_row_binding"));
            }
        }
        if operations > b.max_operations {
            return Err(bad("observed_sam_operation_budget"));
        }
    }
    let keys = input
        .units
        .iter()
        .map(|u| UnitKey {
            source_sha256: input.source_sha256.clone(),
            unit_id: u.unit_id.clone(),
        })
        .collect::<Vec<_>>();
    let mut diagnostics = Vec::new();
    let mut normalized = Vec::new();
    for (unit, key) in input.units.iter().zip(&keys) {
        check()?;
        let mut diagnostic = UnitDiagnostic {
            unit: key.clone(),
            axis_length: unit.declared_axis_length,
            unavailable: BTreeSet::new(),
            rows: Vec::new(),
        };
        let mut segments = Vec::new();
        if unit.rows.is_empty() {
            diagnostic.unavailable.insert(Unavailable::NotProvided);
        }
        for row in &unit.rows {
            check()?;
            let geometry = normalize(row)?;
            if let Some(g) = &geometry {
                if diagnostic.axis_length.is_some_and(|n| n != g.axis_length) {
                    return Err(bad("observed_sam_axis_binding"));
                }
                diagnostic.axis_length = Some(g.axis_length);
            }
            let axis_length = geometry.as_ref().map(|g| g.axis_length);
            let unavailable = if row.flag & 0xeb != 0 {
                Some(Unavailable::Paired)
            } else if row.flag & 4 != 0 {
                Some(Unavailable::Unmapped)
            } else if !input.unpadded_reference {
                Some(Unavailable::PaddedReference)
            } else if row.has_cg
                || geometry.as_ref().is_some_and(|g| {
                    g.operations.len() == 1
                        && g.operations[0].0 == 'N'
                        && g.clipping.soft.iter().any(|&n| n > 0)
                })
            {
                Some(Unavailable::CgOrPlaceholder)
            } else if geometry.is_none() {
                Some(Unavailable::CigarNotProvided)
            } else if axis_length.is_some_and(|n| n > 2_000_000) {
                Some(Unavailable::QueryAxisUnsupported)
            } else if !geometry
                .as_ref()
                .unwrap()
                .operations
                .iter()
                .any(|(c, _)| "M=X".contains(*c))
            {
                Some(Unavailable::NoSupportedBases)
            } else {
                None
            };
            // Even an unsupported mapped row must not conceal malformed target bounds.
            validate_target(row, geometry.as_ref(), &input.dictionary)?;
            diagnostic.rows.push(RowDiagnostic {
                ordinal: row.ordinal,
                unavailable,
                axis_length,
                clipping: geometry.as_ref().map(|g| g.clipping),
            });
            if let Some(reason) = unavailable {
                diagnostic.unavailable.insert(reason);
            }
            if let Some(g) = geometry {
                segments.push((row, g));
            }
        }
        if diagnostic.axis_length.is_some_and(|n| n > 2_000_000) {
            diagnostic
                .unavailable
                .insert(Unavailable::QueryAxisUnsupported);
        }
        if let Some(masks) = &unit.masks {
            if diagnostic.axis_length.is_none()
                || masks
                    .iter()
                    .any(|&(a, z)| a >= z || z > diagnostic.axis_length.unwrap())
                || !masks.windows(2).all(|w| w[0].1 < w[1].0)
            {
                return Err(bad("observed_sam_masks"));
            }
        }
        if !diagnostic.unavailable.is_empty() {
            segments.clear();
        }
        normalized.push(segments);
        diagnostics.push(diagnostic);
    }
    let mut queries = Vec::new();
    let mut unassessed = Vec::new();
    for (((unit, key), segments), diagnostic) in input
        .units
        .iter()
        .zip(&keys)
        .zip(&normalized)
        .zip(&diagnostics)
    {
        if !diagnostic.unavailable.is_empty() {
            unassessed.push(key.clone());
            continue;
        }
        queries.push(core::Query {
            key,
            alignment_sha256: &input.source_sha256,
            query_id: &unit.unit_id,
            query_length: diagnostic.axis_length.expect("admitted mapped unit"),
            scope: PrimerPayloadReferenceScope::Genome,
            masks: unit.masks.as_deref(),
            ends: &unit.ends,
            search_complete: unit.search_declared_complete && unit.views_declared_complete,
            segments: segments
                .iter()
                .map(|(row, g)| core::Segment {
                    query_interval: g.query,
                    target_name: row.target.as_deref().unwrap(),
                    target_interval: g.target,
                    strand: if row.flag & 0x10 == 0 { '+' } else { '-' },
                    mapq: (row.mapq != 255).then_some(row.mapq),
                    alignment_detail_verified: true,
                    alignment_operations: &g.operations,
                    stored_support: None,
                    source_ordinal: Some(row.ordinal),
                    clipping: Some(g.clipping),
                })
                .collect(),
        });
    }
    Ok(Analysis {
        coordinate_authority: "producer_declared",
        identity_scope: "sam_declared_query_unit",
        summary: core::run(
            &context,
            &queries,
            &unassessed,
            &SamIdentity(input),
            work,
            check,
        )?,
        units: diagnostics,
    })
}

#[cfg(test)]
#[path = "sam_tests.rs"]
mod tests;
