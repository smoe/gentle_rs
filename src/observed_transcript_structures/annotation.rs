//! Private incremental GTF preparation; no files, publication, operations or public schema.
//!
//! Hashes cover supplied bytes, not provider completeness or authenticated provenance.
//! Memory charges bound a conservative retained representation plus bounded parser
//! scratch, not host RSS. Limits and gap policies are engineering inputs, not assays.

use super::*;
use crate::engine::{PRIMER_PAYLOAD_ANNOTATION_SCHEMA, PrimerPayloadTranscript};
use ring::digest::{Context as ShaContext, SHA256};

mod gtf;
#[cfg(test)]
mod tests;

#[derive(Debug, Clone)]
pub(super) struct Contig {
    pub name: String,
    pub length: u64,
}

#[derive(Debug, Clone)]
pub(super) struct Dictionary {
    pub reference_id: String,
    pub reference: ReferenceBinding,
    pub contigs: Vec<Contig>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum InputScope {
    WholeSourceFile,
    ProducerExcerpt,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum TranscriptRows {
    Required,
    NotRequired,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) struct Policy {
    pub transcript_feature_rows: TranscriptRows,
    pub minimum_reference_gap_bases: usize,
}

#[derive(Debug, Clone)]
pub(super) struct Context {
    pub reference_id: String,
    pub reference: ReferenceBinding,
    pub selector: ReferenceGene,
    pub provider: Option<String>,
    pub declared_provider_release: Option<String>,
    pub scope: InputScope,
    pub policy: Policy,
}

/// No read/PAF fields and no implicit defaults. All supplied rows consume limits.
#[derive(Debug, Clone)]
pub(super) struct Limits {
    pub max_input_bytes: usize,
    pub max_line_bytes: usize,
    pub max_lines: usize,
    pub max_dictionary_entries: usize,
    pub max_indexed_genes: usize,
    pub max_indexed_transcripts: usize,
    pub max_exon_rows: usize,
    pub max_retained_bytes: usize,
    pub max_work: usize,
    pub max_diagnostics: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum Outcome {
    Completed,
    MissingGene,
    AmbiguousGeneLocus,
    DuplicateTranscriptAccession,
    TranscriptModelConflict,
    IdenticalDuplicateRow,
    ReferenceGapPolicyRefused,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub(super) enum Reason {
    MultipleTranscriptRows,
    MultipleLoci,
    CrossGeneLink,
    AmbiguousGeneLocus,
    IdenticalDuplicateRow,
    UnsupportedStrand,
    VersionConflict,
    VersionPresenceMismatch,
    NoExons,
    MissingTranscriptFeature,
    TranscriptSpanMismatch,
    ExonsOverlapOrTouch,
    GapBelowPolicy,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub(super) enum Entity {
    Gene,
    Transcript,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) struct Diagnostic {
    pub entity: Entity,
    pub id: String,
    pub selected: bool,
    pub reason: Reason,
    pub first_line: usize,
}

#[derive(Debug, Default)]
pub(super) struct Counts {
    pub lines: usize,
    pub comments: usize,
    pub data_rows: usize,
    pub selected_rows: usize,
    pub exon_rows: usize,
    pub selected_exon_rows: usize,
    pub duplicate_rows: usize,
    pub selected_duplicate_rows: usize,
    pub features: BTreeMap<String, usize>,
}

#[derive(Debug)]
pub(super) struct Preparation {
    pub grammar: &'static str,
    pub outcome: Outcome,
    pub annotation: Option<PrimerPayloadAnnotation>,
    pub selector: ReferenceGene,
    pub scope: InputScope,
    pub policy: Policy,
    pub provider: Option<String>,
    pub declared_provider_release: Option<String>,
    pub coordinate_authority: &'static str,
    pub raw_sha256: String,
    pub raw_bytes: usize,
    pub counts: Counts,
    /// Counts are entity/reason pairs, not independent biological events.
    pub selected_reasons: BTreeMap<Reason, usize>,
    pub unrelated_reasons: BTreeMap<Reason, usize>,
    pub diagnostics: Vec<Diagnostic>,
    pub diagnostic_count: usize,
    pub omitted_diagnostics: usize,
    pub indexed_genes: usize,
    pub indexed_transcripts: usize,
    pub gene_version: Option<String>,
    pub transcript_versions: BTreeMap<String, Option<String>>,
    pub included_contigs: Vec<String>,
    pub nonclaims: Vec<&'static str>,
    pub work: usize,
    pub retained_peak_bytes: usize,
}

struct Meter<'a> {
    limits: &'a Limits,
    retained: usize,
    peak: usize,
    work: usize,
}

fn plus(a: usize, b: usize) -> Result<usize, EngineError> {
    a.checked_add(b).ok_or_else(|| bad("observed_gtf_overflow"))
}

impl Meter<'_> {
    fn work(&mut self, n: usize) -> Result<(), EngineError> {
        self.work = plus(self.work, n)?;
        if self.work > self.limits.max_work {
            return Err(bad("observed_gtf_work_budget"));
        }
        Ok(())
    }
    fn scratch(&mut self, n: usize) -> Result<(), EngineError> {
        let peak = plus(self.retained, n)?;
        if peak > self.limits.max_retained_bytes {
            return Err(bad("observed_gtf_memory_budget"));
        }
        self.peak = self.peak.max(peak);
        Ok(())
    }
    fn retain(&mut self, n: usize) -> Result<(), EngineError> {
        self.scratch(n)?;
        self.retained = plus(self.retained, n)?;
        Ok(())
    }
    // Conservative node/object charge includes tree links and spare inline slots.
    fn entry(&mut self, string_bytes: usize) -> Result<(), EngineError> {
        self.retain(plus(512, string_bytes)?)
    }
}

fn text(value: &str, max: usize) -> bool {
    !value.trim().is_empty() && value.len() <= max && !value.chars().any(char::is_control)
}

fn identifier(value: &str, max: usize) -> bool {
    text(value, max) && !value.chars().any(char::is_whitespace)
}

#[derive(Default)]
struct Versions {
    missing: bool,
    first: Option<String>,
    conflict: bool,
    partial: bool,
}

impl Versions {
    fn observe(&mut self, value: Option<&str>, meter: &mut Meter<'_>) -> Result<(), EngineError> {
        if let Some(s) = value {
            if let Some(first) = &self.first {
                self.conflict |= first != s;
            } else {
                meter.retain(s.len())?;
                self.first = Some(s.into());
            }
        } else {
            self.missing = true;
        }
        self.partial = self.missing && self.first.is_some();
        Ok(())
    }
}

type Locus = (usize, char);

fn multiple_loci(loci: &BTreeSet<Locus>) -> bool {
    let contig = loci.first().map(|l| l.0);
    loci.iter().any(|l| Some(l.0) != contig)
        || (loci.iter().any(|l| l.1 == '+') && loci.iter().any(|l| l.1 == '-'))
}

#[derive(Default)]
struct Gene {
    loci: BTreeSet<Locus>,
    version: Versions,
    first_line: usize,
    exons: usize,
    duplicates: usize,
}

#[derive(Default)]
struct Transcript {
    genes: BTreeSet<String>,
    loci: BTreeSet<Locus>,
    version: Versions,
    first_line: usize,
    exons: Vec<Interval>,
    declaration: Option<Interval>,
    declaration_count: usize,
    duplicates: usize,
}

/// One caller-supplied raw line at a time, including its LF/CRLF when present.
/// A failed push poisons the consumer: no partial result can be finished later.
pub(super) struct Consumer<'a> {
    context: &'a Context,
    dictionary: &'a Dictionary,
    lookup: BTreeMap<&'a str, (usize, usize)>,
    meter: Meter<'a>,
    hash: ShaContext,
    bytes: usize,
    counts: Counts,
    genes: BTreeMap<String, Gene>,
    transcripts: BTreeMap<String, Transcript>,
    row_digests: BTreeSet<[u8; 32]>,
    unterminated: bool,
    poisoned: bool,
}

impl<'a> Consumer<'a> {
    pub(super) fn new(
        context: &'a Context,
        dictionary: &'a Dictionary,
        limits: &'a Limits,
        check: &mut dyn FnMut() -> Result<(), EngineError>,
    ) -> Result<Self, EngineError> {
        check()?;
        if [
            limits.max_input_bytes,
            limits.max_line_bytes,
            limits.max_lines,
            limits.max_dictionary_entries,
            limits.max_indexed_genes,
            limits.max_indexed_transcripts,
            limits.max_exon_rows,
            limits.max_retained_bytes,
            limits.max_work,
        ]
        .contains(&0)
            || limits.max_input_bytes > 128 * 1024 * 1024
            || limits.max_line_bytes > 1024 * 1024
            || limits.max_retained_bytes > 128 * 1024 * 1024
            || limits.max_work > 1_000_000_000
            || limits.max_lines > 1_000_000
            || limits.max_dictionary_entries > 100_000
            || limits.max_indexed_genes > 100_000
            || limits.max_indexed_transcripts > 100_000
            || limits.max_exon_rows > 1_000_000
            || limits.max_diagnostics > 1024
            || context.policy.minimum_reference_gap_bases == 0
        {
            return Err(bad("observed_gtf_admission"));
        }
        if !text(&context.reference_id, 4096)
            || !text(&context.reference.assembly, 512)
            || !is_hash(&context.reference.fasta_sha256)
            || !identifier(&context.selector.source_id, 512)
            || !identifier(&context.selector.gene_id, 512)
            || context.provider.as_ref().is_some_and(|s| !text(s, 4096))
            || context
                .declared_provider_release
                .as_ref()
                .is_some_and(|s| !text(s, 4096))
            || dictionary.reference_id != context.reference_id
            || dictionary.reference.assembly != context.reference.assembly
            || dictionary.reference.fasta_sha256 != context.reference.fasta_sha256
            || dictionary.contigs.is_empty()
            || dictionary.contigs.len() > limits.max_dictionary_entries
        {
            return Err(bad("observed_gtf_reference_binding"));
        }
        let mut meter = Meter {
            limits,
            retained: 0,
            peak: 0,
            work: 0,
        };
        meter.entry(std::mem::size_of::<Self>())?;
        let mut lookup = BTreeMap::new();
        for (i, c) in dictionary.contigs.iter().enumerate() {
            check()?;
            meter.work(plus(c.name.len(), 1)?)?;
            let n = usize::try_from(c.length).map_err(|_| bad("observed_gtf_dictionary_length"))?;
            if !identifier(&c.name, 4096)
                || c.name == "*"
                || n == 0
                || lookup.contains_key(c.name.as_str())
            {
                return Err(bad("observed_gtf_dictionary"));
            }
            meter.entry(0)?;
            lookup.insert(c.name.as_str(), (i, n));
        }
        Ok(Self {
            context,
            dictionary,
            lookup,
            meter,
            hash: ShaContext::new(&SHA256),
            bytes: 0,
            counts: Counts::default(),
            genes: BTreeMap::new(),
            transcripts: BTreeMap::new(),
            row_digests: BTreeSet::new(),
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
            return Err(bad("observed_gtf_poisoned"));
        }
        let result = self.push(raw, check);
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
            return Err(bad("observed_gtf_line_after_unterminated"));
        }
        self.bytes = plus(self.bytes, raw.len())?;
        self.counts.lines = plus(self.counts.lines, 1)?;
        if self.bytes > self.meter.limits.max_input_bytes
            || raw.len() > self.meter.limits.max_line_bytes
            || self.counts.lines > self.meter.limits.max_lines
        {
            return Err(bad("observed_gtf_input_budget"));
        }
        self.meter.work(plus(raw.len(), 1)?)?;
        self.meter.scratch(plus(
            raw.len()
                .checked_mul(2)
                .ok_or_else(|| bad("observed_gtf_overflow"))?,
            4096,
        )?)?;
        self.hash.update(raw);
        self.unterminated = !raw.ends_with(b"\n");
        let body = raw.strip_suffix(b"\n").unwrap_or(raw);
        let body = body.strip_suffix(b"\r").unwrap_or(body);
        let line = std::str::from_utf8(body).map_err(|_| bad("observed_gtf_utf8"))?;
        if line.chars().any(|c| c.is_control() && c != '\t') {
            return Err(bad("observed_gtf_control_character"));
        }
        if line.starts_with('#') {
            self.counts.comments += 1;
            return Ok(());
        }
        let row = gtf::parse(line, &self.lookup)?;
        self.counts.data_rows += 1;
        let selected = row.gene == self.context.selector.gene_id;
        self.counts.selected_rows += usize::from(selected);
        if !self.counts.features.contains_key(row.feature) {
            self.meter.entry(row.feature.len())?;
        }
        *self.counts.features.entry(row.feature.into()).or_default() += 1;
        let fingerprint: [u8; 32] = ring::digest::digest(&SHA256, body)
            .as_ref()
            .try_into()
            .expect("SHA256 length");
        let duplicate = self.row_digests.contains(&fingerprint);
        if !duplicate {
            self.meter.entry(32)?;
            self.row_digests.insert(fingerprint);
        }
        self.counts.duplicate_rows += usize::from(duplicate);
        self.counts.selected_duplicate_rows += usize::from(duplicate && selected);
        if !self.genes.contains_key(&row.gene) {
            if self.genes.len() == self.meter.limits.max_indexed_genes {
                return Err(bad("observed_gtf_gene_budget"));
            }
            self.meter
                .entry(plus(row.gene.len(), std::mem::size_of::<Gene>())?)?;
            self.genes.insert(
                row.gene.clone(),
                Gene {
                    first_line: self.counts.lines,
                    ..Gene::default()
                },
            );
        }
        let gene = self.genes.get_mut(&row.gene).unwrap();
        let locus = (row.contig, row.strand);
        if !gene.loci.contains(&locus) {
            self.meter.entry(0)?;
            gene.loci.insert(locus);
        }
        gene.version
            .observe(row.gene_version.as_deref(), &mut self.meter)?;
        gene.duplicates += usize::from(duplicate);
        if row.feature == "exon" {
            if self.counts.exon_rows == self.meter.limits.max_exon_rows {
                return Err(bad("observed_gtf_exon_budget"));
            }
            self.counts.exon_rows += 1;
            self.counts.selected_exon_rows += usize::from(selected);
            gene.exons += 1;
        }
        if let Some(id) = &row.transcript {
            if !self.transcripts.contains_key(id) {
                if self.transcripts.len() == self.meter.limits.max_indexed_transcripts {
                    return Err(bad("observed_gtf_transcript_budget"));
                }
                self.meter
                    .entry(plus(id.len(), std::mem::size_of::<Transcript>())?)?;
                self.transcripts.insert(
                    id.clone(),
                    Transcript {
                        first_line: self.counts.lines,
                        ..Transcript::default()
                    },
                );
            }
            let t = self.transcripts.get_mut(id).unwrap();
            if !t.genes.contains(&row.gene) {
                self.meter.entry(row.gene.len())?;
                t.genes.insert(row.gene.clone());
            }
            if !t.loci.contains(&locus) {
                self.meter.entry(0)?;
                t.loci.insert(locus);
            }
            t.version
                .observe(row.transcript_version.as_deref(), &mut self.meter)?;
            t.duplicates += usize::from(duplicate);
            if row.feature == "transcript" {
                t.declaration_count += 1;
                t.declaration.get_or_insert(row.span);
            }
            if row.feature == "exon" {
                self.meter.retain(4 * std::mem::size_of::<Interval>())?;
                t.exons
                    .try_reserve(1)
                    .map_err(|_| bad("observed_gtf_allocation"))?;
                t.exons.push(row.span);
            }
        }
        Ok(())
    }
}

#[derive(Default)]
struct Diagnoses {
    selected: BTreeMap<Reason, usize>,
    unrelated: BTreeMap<Reason, usize>,
    examples: Vec<Diagnostic>,
    total: usize,
}

impl Diagnoses {
    fn add(
        &mut self,
        entity: Entity,
        id: &str,
        selected: bool,
        reason: Reason,
        line: usize,
        meter: &mut Meter<'_>,
    ) -> Result<(), EngineError> {
        meter.work(1)?;
        let counts = if selected {
            &mut self.selected
        } else {
            &mut self.unrelated
        };
        if !counts.contains_key(&reason) {
            meter.entry(0)?;
        }
        *counts.entry(reason).or_default() += 1;
        self.total = plus(self.total, 1)?;
        if self.examples.len() < meter.limits.max_diagnostics {
            meter.entry(plus(id.len(), std::mem::size_of::<Diagnostic>())?)?;
            self.examples
                .try_reserve(1)
                .map_err(|_| bad("observed_gtf_allocation"))?;
            self.examples.push(Diagnostic {
                entity,
                id: id.into(),
                selected,
                reason,
                first_line: line,
            });
        }
        Ok(())
    }

    fn versions(
        &mut self,
        entity: Entity,
        id: &str,
        selected: bool,
        v: &Versions,
        line: usize,
        meter: &mut Meter<'_>,
    ) -> Result<(), EngineError> {
        if v.conflict {
            self.add(entity, id, selected, Reason::VersionConflict, line, meter)?;
        }
        if v.partial {
            self.add(
                entity,
                id,
                selected,
                Reason::VersionPresenceMismatch,
                line,
                meter,
            )?;
        }
        Ok(())
    }

    fn outcome(&self, present: bool) -> Outcome {
        if !present {
            return Outcome::MissingGene;
        }
        let has = |reason| self.selected.contains_key(&reason);
        if [
            Reason::MultipleTranscriptRows,
            Reason::MultipleLoci,
            Reason::CrossGeneLink,
        ]
        .into_iter()
        .any(has)
        {
            Outcome::DuplicateTranscriptAccession
        } else if has(Reason::AmbiguousGeneLocus) {
            Outcome::AmbiguousGeneLocus
        } else if has(Reason::IdenticalDuplicateRow) {
            Outcome::IdenticalDuplicateRow
        } else if self.selected.keys().any(|r| *r != Reason::GapBelowPolicy) {
            Outcome::TranscriptModelConflict
        } else if has(Reason::GapBelowPolicy) {
            Outcome::ReferenceGapPolicyRefused
        } else {
            Outcome::Completed
        }
    }
}

impl Consumer<'_> {
    /// Explicit end-of-supplied-input, not proof that a provider file was complete.
    pub(super) fn finish(
        self,
        check: &mut dyn FnMut() -> Result<(), EngineError>,
    ) -> Result<Preparation, EngineError> {
        check()?;
        if self.poisoned {
            return Err(bad("observed_gtf_poisoned"));
        }
        let Self {
            context,
            dictionary,
            mut meter,
            mut transcripts,
            genes,
            counts,
            bytes,
            hash,
            ..
        } = self;
        let mut d = Diagnoses::default();
        for (id, g) in &genes {
            check()?;
            meter.work(plus(id.len(), 1)?)?;
            let selected = id == &context.selector.gene_id;
            let mut reason = |r| d.add(Entity::Gene, id, selected, r, g.first_line, &mut meter);
            if multiple_loci(&g.loci) {
                reason(Reason::AmbiguousGeneLocus)?;
            }
            if g.loci.iter().any(|(_, s)| *s == '.') {
                reason(Reason::UnsupportedStrand)?;
            }
            if g.exons == 0 {
                reason(Reason::NoExons)?;
            }
            if g.duplicates > 0 {
                reason(Reason::IdenticalDuplicateRow)?;
            }
            d.versions(
                Entity::Gene,
                id,
                selected,
                &g.version,
                g.first_line,
                &mut meter,
            )?;
        }
        for (id, t) in &mut transcripts {
            check()?;
            meter.work(plus(id.len(), 1)?)?;
            let selected = t.genes.contains(&context.selector.gene_id);
            let mut reasons = BTreeSet::new();
            if t.genes.len() > 1 {
                reasons.insert(Reason::CrossGeneLink);
            }
            if multiple_loci(&t.loci) {
                reasons.insert(Reason::MultipleLoci);
            }
            if t.declaration_count > 1 {
                reasons.insert(Reason::MultipleTranscriptRows);
            }
            if t.duplicates > 0 {
                reasons.insert(Reason::IdenticalDuplicateRow);
            }
            if t.loci.iter().any(|(_, s)| *s == '.') {
                reasons.insert(Reason::UnsupportedStrand);
            }
            if t.exons.is_empty() {
                reasons.insert(Reason::NoExons);
            }
            if context.policy.transcript_feature_rows == TranscriptRows::Required
                && t.declaration_count == 0
            {
                reasons.insert(Reason::MissingTranscriptFeature);
            }
            let n = t.exons.len();
            let sort_work = n
                .checked_mul(if n == 0 { 0 } else { n.ilog2() as usize + 1 })
                .ok_or_else(|| bad("observed_gtf_overflow"))?;
            meter.work(sort_work)?;
            t.exons.sort_unstable();
            // Coordinates from different gene/locus declarations cannot form one model.
            if t.genes.len() == 1 && t.loci.len() == 1 {
                for pair in t.exons.windows(2) {
                    check()?;
                    meter.work(1)?;
                    if pair[0].1 >= pair[1].0 {
                        reasons.insert(Reason::ExonsOverlapOrTouch);
                    } else if pair[1].0 - pair[0].1 < context.policy.minimum_reference_gap_bases {
                        reasons.insert(Reason::GapBelowPolicy);
                    }
                }
                if t.declaration_count == 1
                    && !t.exons.is_empty()
                    && t.declaration != Some((t.exons[0].0, t.exons.last().unwrap().1))
                {
                    reasons.insert(Reason::TranscriptSpanMismatch);
                }
            }
            for r in reasons {
                d.add(
                    Entity::Transcript,
                    id,
                    selected,
                    r,
                    t.first_line,
                    &mut meter,
                )?;
            }
            d.versions(
                Entity::Transcript,
                id,
                selected,
                &t.version,
                t.first_line,
                &mut meter,
            )?;
        }
        let outcome = d.outcome(genes.contains_key(&context.selector.gene_id));
        let raw_sha256 = format!(
            "sha256:{}",
            hash.finish()
                .as_ref()
                .iter()
                .map(|b| format!("{b:02x}"))
                .collect::<String>()
        );
        let gene_version = genes
            .get(&context.selector.gene_id)
            .filter(|g| !g.version.conflict && !g.version.partial)
            .and_then(|g| g.version.first.as_ref());
        let mut models = Vec::new();
        let mut versions = BTreeMap::new();
        let mut included = BTreeSet::new();
        let indexed_transcripts = transcripts.len();
        for (id, t) in transcripts {
            check()?;
            if !t.genes.contains(&context.selector.gene_id) {
                continue;
            }
            let version = if t.version.conflict || t.version.partial {
                None
            } else {
                t.version.first
            };
            meter.entry(plus(id.len(), version.as_ref().map_or(0, String::len))?)?;
            versions.insert(id.clone(), version);
            if outcome != Outcome::Completed {
                continue;
            }
            let &(index, strand) = t.loci.iter().next().expect("completed locus");
            let c = &dictionary.contigs[index];
            meter.entry(plus(
                plus(id.len(), context.selector.gene_id.len())?,
                2 * c.name.len(),
            )?)?;
            if !included.contains(&index) {
                meter.entry(c.name.len())?;
                included.insert(index);
            }
            models
                .try_reserve(1)
                .map_err(|_| bad("observed_gtf_allocation"))?;
            models.push(PrimerPayloadTranscript {
                transcript_id: id,
                gene_ids: vec![context.selector.gene_id.clone()],
                target_name: c.name.clone(),
                target_length: usize::try_from(c.length)
                    .map_err(|_| bad("observed_gtf_dictionary_length"))?,
                genomic_contig: Some(c.name.clone()),
                genomic_strand: Some(strand),
                genomic_exons: t.exons,
            });
        }
        let mut nonclaims = vec![
            "The raw hash covers supplied lines, not independently verified original files or producer provenance.",
            "Scope and reference identity are declared; provider completeness and genome FASTA bytes are not verified.",
            "Missing gene means no exact gene_id match in supplied bytes, not biological absence or an alias search.",
            "Retained exon lists are complete only within supplied rows; transcript declarations cannot prove complete models.",
            "Bare v1 annotation readability does not enforce preparation policies in legacy consumers.",
        ];
        if context.policy.transcript_feature_rows == TranscriptRows::NotRequired {
            nonclaims.push("Without transcript declarations, consistent same-locus accession reuse may be undetectable.");
        }
        let context_size = context.reference_id.len()
            + context.reference.assembly.len()
            + context.reference.fasta_sha256.len()
            + context.selector.source_id.len()
            + context.selector.gene_id.len()
            + context.provider.as_ref().map_or(0, String::len)
            + context
                .declared_provider_release
                .as_ref()
                .map_or(0, String::len);
        meter.entry(plus(context_size, gene_version.map_or(0, String::len))?)?;
        let annotation = (outcome == Outcome::Completed).then(|| PrimerPayloadAnnotation {
            schema: PRIMER_PAYLOAD_ANNOTATION_SCHEMA.into(),
            reference_id: context.reference_id.clone(),
            assembly: context.reference.assembly.clone(),
            release: format!("snapshot-{}", raw_sha256),
            reference_fasta_sha256: context.reference.fasta_sha256.clone(),
            coordinate_convention: "zero_based_half_open".into(),
            transcripts: models,
        });
        check()?;
        Ok(Preparation {
            grammar: "gtf_strict.v1",
            outcome,
            annotation,
            selector: context.selector.clone(),
            scope: context.scope,
            policy: context.policy.clone(),
            provider: context.provider.clone(),
            declared_provider_release: context.declared_provider_release.clone(),
            coordinate_authority: "declared",
            raw_sha256,
            raw_bytes: bytes,
            counts,
            selected_reasons: d.selected,
            unrelated_reasons: d.unrelated,
            omitted_diagnostics: d.total - d.examples.len(),
            diagnostic_count: d.total,
            diagnostics: d.examples,
            indexed_genes: genes.len(),
            indexed_transcripts,
            gene_version: gene_version.cloned(),
            transcript_versions: versions,
            included_contigs: included
                .into_iter()
                .map(|i| dictionary.contigs[i].name.clone())
                .collect(),
            nonclaims,
            work: meter.work,
            retained_peak_bytes: meter.peak,
        })
    }
}
