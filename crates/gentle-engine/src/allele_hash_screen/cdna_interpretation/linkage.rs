//! C3: all eligible same-unit hypotheses must carry the allele and compatible geometry.

use super::*;

fn compatible(
    input: &GeneInput,
    geometry: &AlignmentGeometry,
    model: &CodingModel,
    policy: &Policy,
    meter: &mut Meter<'_>,
) -> Result<bool, String> {
    meter
        .charge(geometry.blocks.len() + geometry.junctions.len() + model.transcript.exons.len())?;
    let introns = model
        .transcript
        .exons
        .windows(2)
        .map(|w| (w[0].1, w[1].0))
        .collect::<Vec<_>>();
    let mut chain = Vec::new();
    for (contig, a, z) in &geometry.junctions {
        meter.charge(1)?;
        if contig != &input.report.reference.contig {
            return Ok(false);
        }
        let Some(index) = introns.iter().position(|i| i == &(*a, *z)) else {
            return Ok(false);
        };
        chain.push(index);
    }
    chain.sort_unstable();
    chain.dedup();
    if chain.windows(2).any(|w| w[0] + 1 != w[1]) {
        return Ok(false);
    }
    let mut ranges = geometry.blocks.clone();
    ranges.sort();
    let mut merged = Vec::<(String, usize, usize)>::new();
    for (contig, a, z) in ranges {
        if a >= z {
            return Err(fail("linkage_geometry"));
        }
        if let Some(last) = merged.last_mut()
            && last.0 == contig
            && a <= last.2
        {
            last.2 = last.2.max(z);
        } else {
            merged.push((contig, a, z));
        }
    }
    if merged.is_empty() {
        return Ok(false);
    }
    let mut outside = 0usize;
    for (contig, a, z) in merged {
        meter.charge(model.transcript.exons.len() + 1)?;
        if contig != input.report.reference.contig {
            return Ok(false);
        }
        let covered = if contig == input.report.reference.contig {
            model
                .transcript
                .exons
                .iter()
                .map(|&(x, y)| z.min(y).saturating_sub(a.max(x)))
                .sum()
        } else {
            0
        };
        outside = outside
            .checked_add((z - a).saturating_sub(covered))
            .ok_or_else(|| fail("linkage_geometry"))?;
        if outside > policy.maximum_out_of_model_aligned_bases {
            return Ok(false);
        }
    }
    Ok(true)
}

// Retained acceleration only: falling back to the same walker changes work, never evidence.
const MAX_CACHE_BYTES: usize = 8 * 1024 * 1024;
const MAX_CACHE_ENTRIES: usize = 512;
const ENTRY_OVERHEAD: usize = 512;
const INDEX_OVERHEAD: usize = 1024;
const INDEX_ENTRY_BYTES: usize = 128;

#[derive(Clone, Copy, Default)]
struct ModelMask([u64; 4]);
impl ModelMask {
    fn insert(&mut self, model: usize) {
        self.0[model / 64] |= 1 << (model % 64);
    }
    fn contains(self, model: usize) -> bool {
        self.0[model / 64] & (1 << (model % 64)) != 0
    }
    fn len(self) -> usize {
        self.0.iter().map(|word| word.count_ones() as usize).sum()
    }
}

#[derive(Clone, Copy)]
struct Witness {
    unit: usize,
    models: ModelMask,
}

struct Entry {
    competitors: Vec<usize>,
    witnesses: Vec<Witness>,
}
impl Entry {
    fn bytes(&self) -> usize {
        ENTRY_OVERHEAD
            .saturating_add(
                self.competitors
                    .capacity()
                    .saturating_mul(std::mem::size_of::<usize>()),
            )
            .saturating_add(
                self.witnesses
                    .capacity()
                    .saturating_mul(std::mem::size_of::<Witness>()),
            )
    }
    fn retain(
        &mut self,
        witness: Witness,
        available: usize,
        meter: &mut Meter<'_>,
    ) -> Result<bool, String> {
        if self.witnesses.len() == self.witnesses.capacity() {
            let capacity = self.witnesses.capacity().max(1).saturating_mul(2);
            let proposed = self.bytes().saturating_add(
                capacity
                    .saturating_sub(self.witnesses.capacity())
                    .saturating_mul(std::mem::size_of::<Witness>()),
            );
            if proposed > available {
                return Ok(false);
            }
            meter.charge(
                self.witnesses
                    .len()
                    .saturating_mul(std::mem::size_of::<Witness>())
                    .saturating_add(1),
            )?;
            self.witnesses
                .reserve_exact(capacity - self.witnesses.len());
            if self.bytes() > available {
                return Ok(false);
            }
        }
        self.witnesses.push(witness);
        Ok(true)
    }
}

fn lookup_work(entries: usize) -> usize {
    entries.checked_ilog2().map_or(1, |n| n as usize + 2)
}

/// One immutable gene/sample context; cached masks never cross a quality or geometry context.
/// Retention includes a conservative tree/index allowance and is at most eight MiB. A record
/// which cannot fit is streamed through the same walker without retaining it or dropping input.
pub(super) struct Cache<'a> {
    input: &'a GeneInput,
    evidence: &'a LinkageEvidence,
    sample: &'a evidence::CdnaSampleReport,
    source: Option<&'a evidence::CdnaSample>,
    unit_index: Option<BTreeMap<&'a str, usize>>,
    index_bytes: usize,
    entries: BTreeMap<(usize, char, usize), Entry>,
    entry_bytes: usize,
    max_bytes: usize,
    max_entries: usize,
}

impl<'a> Cache<'a> {
    pub(super) fn new(
        input: &'a GeneInput,
        evidence: &'a LinkageEvidence,
        sample: &'a evidence::CdnaSampleReport,
        meter: &mut Meter<'_>,
    ) -> Result<Self, String> {
        let max_bytes = meter.remaining_output.min(MAX_CACHE_BYTES);
        Self::initialize(input, evidence, sample, meter, max_bytes, MAX_CACHE_ENTRIES)
    }

    fn initialize(
        input: &'a GeneInput,
        evidence: &'a LinkageEvidence,
        sample: &'a evidence::CdnaSampleReport,
        meter: &mut Meter<'_>,
        max_bytes: usize,
        max_entries: usize,
    ) -> Result<Self, String> {
        meter.charge(0)?;
        if input.models.len() > 256 {
            return Err(fail("model_geometry_mismatch"));
        }
        Ok(Self {
            input,
            evidence,
            sample,
            source: None,
            unit_index: None,
            index_bytes: 0,
            entries: BTreeMap::new(),
            entry_bytes: 0,
            max_bytes: max_bytes.min(MAX_CACHE_BYTES),
            max_entries: max_entries.min(MAX_CACHE_ENTRIES),
        })
    }

    #[cfg(test)]
    pub(super) fn with_limits(
        input: &'a GeneInput,
        evidence: &'a LinkageEvidence,
        sample: &'a evidence::CdnaSampleReport,
        meter: &mut Meter<'_>,
        bytes: usize,
        entries: usize,
    ) -> Result<Self, String> {
        Self::initialize(input, evidence, sample, meter, bytes, entries)
    }

    #[cfg(test)]
    pub(super) fn retained_bytes(&self) -> usize {
        self.index_bytes + self.entry_bytes
    }

    #[cfg(test)]
    pub(super) fn retained_entries(&self) -> usize {
        self.entries.len()
    }

    fn source(&mut self, meter: &mut Meter<'_>) -> Result<&'a evidence::CdnaSample, String> {
        if let Some(source) = self.source {
            return Ok(source);
        }
        meter.charge(self.evidence.samples.len().saturating_add(1))?;
        let source = self
            .evidence
            .samples
            .iter()
            .find(|s| s.sample_id == self.sample.sample_id)
            .ok_or_else(|| fail("linkage_sample"))?;
        self.source = Some(source);
        let bytes =
            INDEX_OVERHEAD.saturating_add(source.units.len().saturating_mul(INDEX_ENTRY_BYTES));
        if bytes <= self.max_bytes {
            let mut index = BTreeMap::new();
            for (i, unit) in source.units.iter().enumerate() {
                meter.charge(lookup_work(index.len()).saturating_add(1))?;
                // Preserve the previous linear lookup's first-match behavior.
                index.entry(unit.unit_id.as_str()).or_insert(i);
            }
            self.unit_index = Some(index);
            self.index_bytes = bytes;
        }
        Ok(source)
    }

    fn evict(&mut self, meter: &mut Meter<'_>) -> Result<(), String> {
        for entry in self.entries.values() {
            meter.charge(
                entry
                    .witnesses
                    .len()
                    .saturating_add(entry.competitors.len())
                    .saturating_add(1),
            )?;
        }
        self.entries.clear();
        self.entry_bytes = 0;
        Ok(())
    }

    pub(super) fn link(
        &mut self,
        model: &CodingModel,
        p: usize,
        alt: char,
        policy: &Policy,
        meter: &mut Meter<'_>,
    ) -> Result<IsoformLinkage, String> {
        let mut result = empty();
        meter.charge(model.transcript.exons.len().saturating_add(1))?;
        if position(model, p).is_none() {
            result.status = "not_present".into();
            return Ok(result);
        }
        meter.charge(self.input.models.len().saturating_add(1))?;
        let model_index = self
            .input
            .models
            .iter()
            .position(|m| m.transcript.transcript_id == model.transcript.transcript_id)
            .ok_or_else(|| fail("model_geometry_mismatch"))?;
        let key = (p, alt, policy.maximum_out_of_model_aligned_bases);
        meter.charge(lookup_work(self.entries.len()).saturating_add(1))?;
        if let Some(entry) = self.entries.get(&key) {
            let source = self.source.ok_or_else(|| fail("linkage_sample"))?;
            for witness in &entry.witnesses {
                project(
                    &mut result,
                    self.input,
                    source,
                    &entry.competitors,
                    *witness,
                    model_index,
                    meter,
                )?;
            }
            return finish(result, meter);
        }

        let source = self.source(meter)?;
        meter.charge(
            self.input
                .models
                .iter()
                .map(|m| m.transcript.exons.len().saturating_add(1))
                .sum(),
        )?;
        let competitors = self
            .input
            .models
            .iter()
            .enumerate()
            .filter(|(_, m)| position(m, p).is_some())
            .map(|(i, _)| i)
            .collect::<Vec<_>>();
        if self.entries.len() >= self.max_entries
            || self
                .max_bytes
                .saturating_sub(self.index_bytes + self.entry_bytes)
                < ENTRY_OVERHEAD
        {
            self.evict(meter)?;
        }
        let available = self
            .max_bytes
            .saturating_sub(self.index_bytes + self.entry_bytes);
        meter.charge(competitors.len())?;
        let mut entry = Entry {
            competitors: competitors.clone(),
            witnesses: Vec::new(),
        };
        let mut retaining = self.max_entries > 0 && entry.bytes() <= available;
        walk(
            self.input,
            self.evidence,
            self.sample,
            source,
            self.unit_index.as_ref(),
            &competitors,
            p,
            alt,
            policy,
            meter,
            &mut |witness, meter| {
                project(
                    &mut result,
                    self.input,
                    source,
                    &competitors,
                    witness,
                    model_index,
                    meter,
                )?;
                if retaining {
                    // Only primitive masks/indices are retained, not repeated output strings.
                    meter.charge(1)?;
                    if !entry.retain(witness, available, meter)? {
                        retaining = false;
                        entry.witnesses = Vec::new();
                    }
                }
                Ok(())
            },
        )?;
        let result = finish(result, meter)?;
        if retaining {
            meter.charge(lookup_work(self.entries.len()).saturating_add(1))?;
            self.entry_bytes = self.entry_bytes.saturating_add(entry.bytes());
            self.entries.insert(key, entry);
        }
        Ok(result)
    }
}

fn empty() -> IsoformLinkage {
    IsoformLinkage {
        status: "no_compatible_units".into(),
        compatible_units: Vec::new(),
        discriminating_units: Vec::new(),
        unique_model_units: Vec::new(),
        excluded_transcripts: BTreeMap::new(),
        reasons: Vec::new(),
    }
}

fn walk(
    input: &GeneInput,
    evidence: &LinkageEvidence,
    sample: &evidence::CdnaSampleReport,
    source: &evidence::CdnaSample,
    index: Option<&BTreeMap<&str, usize>>,
    competitors: &[usize],
    p: usize,
    alt: char,
    policy: &Policy,
    meter: &mut Meter<'_>,
    visit: &mut dyn FnMut(Witness, &mut Meter<'_>) -> Result<(), String>,
) -> Result<(), String> {
    meter.charge(sample.site_calls.len())?;
    for call in sample
        .site_calls
        .iter()
        .filter(|c| c.position_0based == p && c.base == Some(alt) && definite(c))
    {
        let unit_index = if let Some(index) = index {
            meter.charge(lookup_work(index.len()).saturating_add(1))?;
            *index
                .get(call.unit_id.as_str())
                .ok_or_else(|| fail("linkage_unit"))?
        } else {
            meter.charge(source.units.len().saturating_add(1))?;
            source
                .units
                .iter()
                .position(|u| u.unit_id == call.unit_id)
                .ok_or_else(|| fail("linkage_unit"))?
        };
        let unit = &source.units[unit_index];
        if !unit.all_query_views_assessed {
            continue;
        }
        meter.charge(unit.hypotheses.len().saturating_add(1))?;
        let hypotheses = unit
            .hypotheses
            .iter()
            .filter(|h| h.eligible)
            .collect::<Vec<_>>();
        if hypotheses.is_empty() {
            continue;
        }
        let mut admissible = ModelMask::default();
        for &candidate in competitors {
            let mut all = true;
            for hypothesis in &hypotheses {
                meter.charge(hypothesis.bases.len().saturating_add(1))?;
                let carries = hypothesis.bases.iter().any(|b| {
                    b.reference_position_0based == p
                        && b.observed_base == alt
                        && b.phred
                            .is_some_and(|q| q >= input.report.policy.minimum_base_quality)
                        && b.mapq.is_some_and(|q| {
                            q != 255 && q >= input.report.policy.minimum_mapping_quality
                        })
                });
                meter.charge(lookup_work(evidence.geometry.len()).saturating_add(1))?;
                let geometry = evidence
                    .geometry
                    .get(&(
                        sample.sample_id.clone(),
                        unit.unit_id.clone(),
                        hypothesis.hypothesis_id.clone(),
                    ))
                    .ok_or_else(|| fail("linkage_hypothesis"))?;
                if !carries
                    || !compatible(input, geometry, &input.models[candidate], policy, meter)?
                {
                    all = false;
                    break;
                }
            }
            if all {
                admissible.insert(candidate);
            }
        }
        visit(
            Witness {
                unit: unit_index,
                models: admissible,
            },
            meter,
        )?;
    }
    Ok(())
}

fn project(
    result: &mut IsoformLinkage,
    input: &GeneInput,
    source: &evidence::CdnaSample,
    competitors: &[usize],
    witness: Witness,
    model: usize,
    meter: &mut Meter<'_>,
) -> Result<(), String> {
    meter.charge(5)?;
    if !witness.models.contains(model) {
        return Ok(());
    }
    let unit = &source.units[witness.unit].unit_id;
    meter.charge(unit.len().saturating_add(1))?;
    result.compatible_units.push(unit.clone());
    let admitted = witness.models.len();
    if admitted < competitors.len() {
        meter.charge(unit.len().saturating_add(1))?;
        result.discriminating_units.push(unit.clone());
        for &candidate in competitors {
            meter.charge(1)?;
            if !witness.models.contains(candidate) {
                let id = &input.models[candidate].transcript.transcript_id;
                meter.charge(
                    lookup_work(result.excluded_transcripts.len())
                        .saturating_add(id.len())
                        .saturating_add(unit.len())
                        .saturating_add(1),
                )?;
                result
                    .excluded_transcripts
                    .entry(id.clone())
                    .or_default()
                    .push(unit.clone());
            }
        }
    }
    if admitted == 1 {
        meter.charge(unit.len().saturating_add(1))?;
        result.unique_model_units.push(unit.clone());
    }
    Ok(())
}

fn finish(mut result: IsoformLinkage, meter: &mut Meter<'_>) -> Result<IsoformLinkage, String> {
    for units in [
        &mut result.compatible_units,
        &mut result.discriminating_units,
        &mut result.unique_model_units,
    ] {
        meter.charge(
            units
                .len()
                .saturating_mul(lookup_work(units.len()))
                .saturating_add(1),
        )?;
        units.sort();
    }
    for units in result.excluded_transcripts.values_mut() {
        meter.charge(
            units
                .len()
                .saturating_mul(lookup_work(units.len()))
                .saturating_add(1),
        )?;
        units.sort();
        units.dedup();
    }
    result.status = if !result.discriminating_units.is_empty() {
        "same_unit_discriminating_support"
    } else if !result.compatible_units.is_empty() {
        "same_unit_compatible_only"
    } else {
        "no_compatible_units"
    }
    .into();
    result.reasons.push("Compatibility is conditional on supplied eligible hypotheses, not identification of a full-length or expressed isoform.".into());
    Ok(result)
}

#[cfg(test)]
pub(super) mod legacy_oracle;
