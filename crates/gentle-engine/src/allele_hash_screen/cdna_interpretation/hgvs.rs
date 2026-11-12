//! Bounded substitution/composite notation on an explicit transcript reference (HGVS 21.1.5).
//! No intronic numbering, inferred initiation, full-length RNA or clinical classification.

use super::*;

fn accession(value: Option<&str>) -> Option<&str> {
    value.filter(|v| {
        v.rsplit_once('.').is_some_and(|(base, version)| {
            !base.is_empty()
                && !version.is_empty()
                && version.bytes().all(|b| b.is_ascii_digit())
                && base
                    .bytes()
                    .all(|b| b.is_ascii_alphanumeric() || matches!(b, b'_' | b'.'))
        })
    })
}

pub(super) fn amino(aa: char) -> Option<&'static str> {
    Some(match aa {
        'A' => "Ala",
        'R' => "Arg",
        'N' => "Asn",
        'D' => "Asp",
        'C' => "Cys",
        'Q' => "Gln",
        'E' => "Glu",
        'G' => "Gly",
        'H' => "His",
        'I' => "Ile",
        'L' => "Leu",
        'K' => "Lys",
        'M' => "Met",
        'F' => "Phe",
        'P' => "Pro",
        'S' => "Ser",
        'T' => "Thr",
        'W' => "Trp",
        'Y' => "Tyr",
        'V' => "Val",
        '*' => "Ter",
        _ => return None,
    })
}

fn coordinate(model: &CodingModel, p: usize) -> String {
    match model.cds {
        Some((start, _)) if p < start => format!("-{}", start - p),
        Some((_, end)) if p >= end => format!("*{}", p - end + 1),
        Some((start, _)) => (p - start + 1).to_string(),
        None => (p + 1).to_string(),
    }
}

fn region(model: &CodingModel, p: usize) -> &'static str {
    match model.cds {
        Some((start, _)) if p < start => "utr5",
        Some((_, end)) if p >= end => "utr3",
        Some(_) => "cds",
        None if model.translation_status == "coding_not_defined" => "coding_geometry_unavailable",
        None => "noncoding_exon",
    }
}

fn nucleotide(
    model: &CodingModel,
    changes: &[Change],
    meter: &mut Meter<'_>,
) -> Result<Option<String>, String> {
    let Some(accession) = accession(model.accession.as_deref()) else {
        return Ok(None);
    };
    if model.cds.is_some() && model.translation_status != "complete_verified" {
        return Ok(None);
    }
    if model.translation_status == "coding_not_defined" {
        return Ok(None);
    }
    let mut changes = changes
        .iter()
        .filter(|c| c.reference_base != c.alternate_base)
        .collect::<Vec<_>>();
    changes.sort_by_key(|c| c.cdna_position_0based);
    let Some(first) = changes.first() else {
        return Ok(None);
    };
    let prefix = if model.cds.is_some() { "c" } else { "n" };
    let atom = |c: &Change| {
        format!(
            "{}{}>{}",
            coordinate(model, c.cdna_position_0based),
            c.reference_base,
            c.alternate_base
        )
    };
    let body = if changes.len() == 1 {
        atom(first)
    } else {
        let last = changes.last().unwrap();
        let consecutive = changes
            .windows(2)
            .all(|w| w[0].cdna_position_0based + 1 == w[1].cdna_position_0based);
        let same_codon_one_gap = model.cds.is_some_and(|(start, end)| {
            first.cdna_position_0based >= start
                && last.cdna_position_0based < end
                && (first.cdna_position_0based - start) / 3
                    == (last.cdna_position_0based - start) / 3
                && changes.len() == 2
                && last.cdna_position_0based - first.cdna_position_0based == 2
        });
        if consecutive || same_codon_one_gap {
            let (a, z) = (first.cdna_position_0based, last.cdna_position_0based + 1);
            meter.charge(z - a)?;
            let mut replacement = model.sequence.as_bytes()[a..z].to_vec();
            for c in &changes {
                replacement[c.cdna_position_0based - a] = c.alternate_base as u8;
            }
            format!(
                "{}_{}delins{}",
                coordinate(model, a),
                coordinate(model, z - 1),
                String::from_utf8(replacement).unwrap()
            )
        } else {
            format!(
                "[{}]",
                changes
                    .iter()
                    .map(|c| atom(c))
                    .collect::<Vec<_>>()
                    .join(";")
            )
        }
    };
    Ok(Some(format!("{accession}:{prefix}.{body}")))
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CodonEffect {
    pub codon_index_1based: usize,
    pub reference_codon: String,
    pub alternate_codon: String,
    pub reference_amino_acid: String,
    pub alternate_amino_acid: String,
    pub consequence_term: String,
    pub protein_hgvs: Option<String>,
}

pub(super) fn consequence(
    model: &CodingModel,
    changes: &[Change],
    code: &BTreeMap<String, char>,
    meter: &mut Meter<'_>,
) -> Result<Consequence, String> {
    let mut result = Consequence {
        assessable: false,
        region: "not_present".into(),
        cdna_position_1based: None,
        nucleotide_hgvs: None,
        protein_hgvs: None,
        codon_index_1based: None,
        reference_codon: None,
        alternate_codon: None,
        reference_amino_acid: None,
        alternate_amino_acid: None,
        consequence_terms: Vec::new(),
        reasons: Vec::new(),
        coding_changes: Vec::new(),
    };
    if changes.is_empty() {
        result.reasons.push("not_present".into());
        return Ok(result);
    }
    let mut positions = BTreeSet::new();
    for c in changes {
        meter.charge(1)?;
        if !positions.insert(c.cdna_position_0based)
            || model
                .sequence
                .as_bytes()
                .get(c.cdna_position_0based)
                .copied()
                .map(char::from)
                != Some(c.reference_base)
            || !matches!(c.alternate_base, 'A' | 'C' | 'G' | 'T')
        {
            return Err(fail("codon_change_binding"));
        }
    }
    let regions = changes
        .iter()
        .map(|c| region(model, c.cdna_position_0based))
        .collect::<BTreeSet<_>>();
    result.region = if regions.len() == 1 {
        *regions.first().unwrap()
    } else {
        "mixed_regions"
    }
    .into();
    result.cdna_position_1based = Some(
        changes
            .iter()
            .map(|c| c.cdna_position_0based)
            .min()
            .unwrap()
            + 1,
    );
    result.nucleotide_hgvs = nucleotide(model, changes, meter)?;
    if accession(model.accession.as_deref()).is_none() {
        result
            .reasons
            .push("transcript_version_unavailable_hgvs_withheld".into());
    }
    if model.translation_status == "coding_not_defined" {
        result.reasons.push("coding_geometry_unavailable".into());
        return Ok(result);
    }
    if model.cds.is_some() && model.translation_status != "complete_verified" {
        result
            .reasons
            .extend(model.translation_reasons.iter().cloned());
        result
            .reasons
            .push(format!("translation_status:{}", model.translation_status));
        if regions.contains("cds") {
            result
                .consequence_terms
                .push("coding_sequence_variant".into());
        }
        return Ok(result);
    }
    result.assessable = true;
    for r in &regions {
        if let Some(term) = match *r {
            "utr5" => Some("5_prime_UTR_variant"),
            "utr3" => Some("3_prime_UTR_variant"),
            "noncoding_exon" => Some("non_coding_transcript_exon_variant"),
            _ => None,
        } {
            result.consequence_terms.push(term.into());
        }
    }
    if let Some((start, end)) = model.cds {
        let codons = changes
            .iter()
            .filter(|c| start <= c.cdna_position_0based && c.cdna_position_0based < end)
            .map(|c| (c.cdna_position_0based - start) / 3)
            .collect::<BTreeSet<_>>();
        for index in codons {
            meter.charge(1)?;
            let a = start + index * 3;
            let reference = model
                .sequence
                .get(a..a + 3)
                .ok_or_else(|| fail("codon_geometry"))?
                .to_owned();
            let mut alternate = reference.as_bytes().to_vec();
            for c in changes
                .iter()
                .filter(|c| a <= c.cdna_position_0based && c.cdna_position_0based < a + 3)
            {
                alternate[c.cdna_position_0based - a] = c.alternate_base as u8;
            }
            let alternate = String::from_utf8(alternate).unwrap();
            let old = *code
                .get(&reference)
                .ok_or_else(|| fail("codon_not_assessable"))?;
            let new = *code
                .get(&alternate)
                .ok_or_else(|| fail("codon_not_assessable"))?;
            let old_name = amino(old).ok_or_else(|| fail("genetic_code"))?;
            let new_name = amino(new).ok_or_else(|| fail("genetic_code"))?;
            let (term, body) = if old == new {
                (
                    if old == '*' {
                        "stop_retained_variant"
                    } else {
                        "synonymous_variant"
                    },
                    format!("{old_name}{}=", index + 1),
                )
            } else if index == 0 {
                ("start_lost", "Met1?".into())
            } else if new == '*' {
                ("stop_gained", format!("{old_name}{}Ter", index + 1))
            } else if old == '*' {
                let mut extension = None;
                for (i, codon) in model.sequence.as_bytes()[end..].chunks_exact(3).enumerate() {
                    meter.charge(changes.len() + 1)?;
                    let a = end + i * 3;
                    let mut triplet = codon.to_vec();
                    for c in changes
                        .iter()
                        .filter(|c| a <= c.cdna_position_0based && c.cdna_position_0based < a + 3)
                    {
                        triplet[c.cdna_position_0based - a] = c.alternate_base as u8;
                    }
                    let codon = std::str::from_utf8(&triplet).unwrap();
                    match code.get(codon) {
                        Some('*') => {
                            // The former stop, now an amino acid, starts the added sequence at 1.
                            extension = Some(i + 2);
                            break;
                        }
                        Some(_) => (),
                        None => break,
                    }
                }
                (
                    "stop_lost",
                    format!(
                        "Ter{}{new_name}extTer{}",
                        index + 1,
                        extension.map_or_else(|| "?".into(), |n| n.to_string())
                    ),
                )
            } else {
                (
                    "missense_variant",
                    format!("{old_name}{}{new_name}", index + 1),
                )
            };
            let protein_hgvs = accession(model.protein_accession.as_deref())
                .filter(|_| accession(model.accession.as_deref()).is_some())
                .map(|p| format!("{p}:p.({body})"));
            result.consequence_terms.push(term.into());
            result.coding_changes.push(CodonEffect {
                codon_index_1based: index + 1,
                reference_codon: reference,
                alternate_codon: alternate,
                reference_amino_acid: old_name.into(),
                alternate_amino_acid: new_name.into(),
                consequence_term: term.into(),
                protein_hgvs,
            });
        }
        if let [only] = result.coding_changes.as_slice() {
            result.codon_index_1based = Some(only.codon_index_1based);
            result.reference_codon = Some(only.reference_codon.clone());
            result.alternate_codon = Some(only.alternate_codon.clone());
            result.reference_amino_acid = Some(only.reference_amino_acid.clone());
            result.alternate_amino_acid = Some(only.alternate_amino_acid.clone());
            result.protein_hgvs = only.protein_hgvs.clone();
        } else if result.coding_changes.len() > 1 {
            result
                .reasons
                .push("multi_codon_composite_not_rendered".into());
        }
        if !result.coding_changes.is_empty()
            && accession(model.protein_accession.as_deref()).is_none()
        {
            result
                .reasons
                .push("protein_version_unavailable_hgvs_withheld".into());
        }
    }
    result.consequence_terms.sort();
    result.consequence_terms.dedup();
    Ok(result)
}
