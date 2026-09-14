//! Read-only, headless presentation of GENtle's annotated TSS sequence exports.
//!
//! This decodes the existing export grammar, not arbitrary annotation prose.
//! It does not infer TSSs, fetch evidence, score motifs, or alter source records.

use std::collections::BTreeMap;

use gb_io::seq::{Feature, Location};
use gentle_protocol::tss_profiles::{TssGeometry, TssStrand};
use serde::Serialize;

use crate::{digest_utils::sha256_hex_bytes, dna_sequence::DNAsequence};

/// Evidence classes remain separate even when their genomic intervals overlap.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Serialize)]
pub enum TssLaneKind {
    Structure,
    Signal,
    Motif,
    Other,
}

/// One original feature, in local zero-based, end-exclusive coordinates.
#[derive(Clone, Debug, Serialize)]
pub struct TssViewFeature {
    pub feature_id: usize,
    pub start: usize,
    pub end: usize,
    pub reverse: bool,
    pub clipped: bool,
    pub label: String,
    pub details: String,
    pub score: Option<f64>,
}

/// An exact source/sample, transcript, or matrix/score-kind lane.
#[derive(Clone, Debug, Serialize)]
pub struct TssViewLane {
    pub kind: TssLaneKind,
    pub id: String,
    pub label: String,
    pub details: String,
    pub units: String,
    pub state: String,
    pub features: Vec<TssViewFeature>,
    /// Retains source-locus scaling for signal, never normalizes across matrices.
    pub scale_max: f64,
}

/// Native-view input shared with tests and future adapters, not a new assay report.
#[derive(Clone, Debug, Serialize)]
pub struct TssSequenceView {
    pub title: String,
    pub promoter_id: String,
    pub assembly: String,
    pub geometry: TssGeometry,
    pub sequence_sha256: String,
    pub lanes: Vec<TssViewLane>,
    pub provenance: String,
    pub warnings: Vec<String>,
}

// Each exported field has an explicit key; never classify by arbitrary filenames.
fn field(text: &str, key: &str) -> Option<String> {
    text.split([';', '\n'])
        .map(str::trim)
        .find_map(|s| s.strip_prefix(key).and_then(|s| s.strip_prefix('=')))
        .map(|s| s.trim().to_string())
}

fn required(text: &str, key: &str) -> Result<String, String> {
    field(text, key)
        .filter(|s| !s.is_empty())
        .ok_or_else(|| format!("TSS metadata lacks {key}"))
}

fn number(text: &str, key: &str) -> Result<u64, String> {
    required(text, key)?
        .parse()
        .map_err(|_| format!("Invalid TSS {key}"))
}

fn qualifier(f: &Feature, key: &str) -> String {
    f.qualifiers
        .iter()
        .filter(|(k, _)| k.as_ref() == key)
        .filter_map(|(_, v)| v.as_deref())
        .map(|s| {
            s.trim_matches('"')
                .split_whitespace()
                .collect::<Vec<_>>()
                .join(" ")
        })
        .collect::<Vec<_>>()
        .join("; ")
}

fn bounds(loc: &Location) -> Option<(usize, usize, bool, bool)> {
    match loc {
        Location::Range((start, before), (end, after)) => Some((
            (*start).try_into().ok()?,
            (*end).try_into().ok()?,
            false,
            before.0 || after.0,
        )),
        Location::Complement(inner) => {
            bounds(inner).map(|(s, e, rev, clipped)| (s, e, !rev, clipped))
        }
        _ => None,
    }
}

fn finite_score(text: &str) -> Option<f64> {
    field(text, "raw_score")?
        .parse::<f64>()
        .ok()
        .filter(|v| v.is_finite())
}

fn validate_source_span(
    note: &str,
    g: &TssGeometry,
    start: usize,
    end: usize,
) -> Result<(), String> {
    let Some(raw) = note.strip_prefix("Source genomic interval ") else {
        return Ok(());
    };
    let (lo, hi) = raw
        .split([' ', ';'])
        .next()
        .and_then(|s| s.split_once(".."))
        .and_then(|(s, e)| Some((s.parse::<u64>().ok()?, e.parse::<u64>().ok()?)))
        .filter(|(s, e)| *s > 0 && s <= e)
        .ok_or("Invalid original genomic interval")?;
    let lo = lo.max(g.start_1based);
    let hi = hi.min(g.end_1based);
    if lo > hi {
        return Err("Original genomic interval does not overlap the TSS window".into());
    }
    let expected = if g.strand == TssStrand::Plus {
        (lo - g.start_1based, hi - g.start_1based + 1)
    } else {
        (g.end_1based - hi, g.end_1based - lo + 1)
    };
    if (start as u64, end as u64) != expected {
        return Err(
            "Feature location contradicts the original genomic interval and window strand".into(),
        );
    }
    Ok(())
}

impl TssSequenceView {
    /// Cheap recognition only. A recognized document still needs `from_dna` validation.
    pub fn recognizes(dna: &DNAsequence) -> bool {
        dna.description()
            .iter()
            .any(|s| s.starts_with("GENtle promoter_id="))
    }

    /// Decode a GENtle annotated TSS export. Missing metadata is not guessed.
    /// Hash binding protects coordinates, not authenticity or unchanged feature annotations.
    pub fn from_dna(dna: &DNAsequence) -> Result<Self, String> {
        if !Self::recognizes(dna) {
            return Err("No GENtle annotated TSS metadata. Open an annotated TSS EMBL/GenBank export; FASTA alone has no evidence lanes.".into());
        }
        if dna.is_circular() || dna.len() > 2_000_000 || dna.features().len() > 100_000 {
            return Err(
                "TSS inspection requires a linear window up to 2 Mb and 100,000 features".into(),
            );
        }
        let provenance = dna.description().join("\n");
        let metadata = dna
            .description()
            .iter()
            .map(|s| s.split_whitespace().collect::<Vec<_>>().join(" "))
            .collect::<Vec<_>>()
            .join("\n");
        let promoter_id = required(&metadata, "GENtle promoter_id")?;
        // GenBank can hard-wrap this opaque 64-digit token inside a COMMENT.
        // Whitespace has no meaning in a SHA-256 hex value; retain raw provenance separately.
        let sequence_sha256 = required(&metadata, "sequence_sha256")?
            .split_whitespace()
            .collect::<String>();
        if sha256_hex_bytes(dna.get_forward_string().to_ascii_uppercase().as_bytes())
            != sequence_sha256
        {
            return Err("TSS sequence hash mismatch: sequence was edited or reoriented; stored genomic coordinates cannot be used".into());
        }
        if required(&metadata, "local_axis")? != "transcript_5prime_to_3prime" {
            return Err("Unsupported TSS local-axis convention".into());
        }
        let (start, end) = required(&metadata, "genomic")?
            .split_once("..")
            .and_then(|(s, e)| Some((s.parse::<u64>().ok()?, e.parse::<u64>().ok()?)))
            .ok_or("Invalid TSS genomic span")?;
        if start == 0
            || end.checked_sub(start).and_then(|n| n.checked_add(1)) != Some(dna.len() as u64)
        {
            return Err("TSS genomic span does not match sequence length".into());
        }
        let local = number(&metadata, "TSS_local_1based")?;
        if local == 0 || local > dna.len() as u64 {
            return Err("TSS lies outside sequence".into());
        }
        let strand = match required(&metadata, "genomic_strand")?.as_str() {
            "+" => TssStrand::Plus,
            "-" => TssStrand::Minus,
            _ => return Err("Invalid TSS genomic strand".into()),
        };
        let geometry = TssGeometry {
            chromosome: required(&metadata, "chromosome")?,
            strand,
            tss_1based: if strand == TssStrand::Plus {
                start + (local - 1)
            } else {
                end - (local - 1)
            },
            start_1based: start,
            end_1based: end,
            upstream_bp: (local - 1) as usize,
            downstream_bp: dna.len() - local as usize,
        };
        let assembly = required(&metadata, "assembly")?;
        let gene = dna
            .features()
            .iter()
            .find(|f| f.kind.as_ref() == "source")
            .map(|f| qualifier(f, "label"))
            .filter(|s| !s.is_empty())
            .unwrap_or_else(|| "TSS".into());
        let mut lanes: BTreeMap<(TssLaneKind, String), TssViewLane> = BTreeMap::new();
        // EMBL may split a comment into several CC lines; retain whole source metadata.
        for entry in metadata.split("Signal lane ").skip(1) {
            let Some((id, rest)) = entry.split_once(": ") else {
                continue;
            };
            let label = rest
                .split(';')
                .next()
                .unwrap_or(id)
                .split_whitespace()
                .collect::<Vec<_>>()
                .join(" ");
            let identity = field(rest, "mark")
                .filter(|s| s != "not supplied")
                .or_else(|| field(rest, "factor").filter(|s| s != "not supplied"));
            let label =
                identity.map_or_else(|| label.clone(), |identity| format!("{identity} | {label}"));
            let scale_max = field(rest, "inherited_abs_max")
                .and_then(|s| s.parse::<f64>().ok())
                .filter(|v| v.is_finite() && *v > 0.0)
                .unwrap_or(1.0);
            lanes.insert(
                (TssLaneKind::Signal, id.into()),
                TssViewLane {
                    kind: TssLaneKind::Signal,
                    id: id.into(),
                    label,
                    details: rest.to_string(),
                    units: "source signal (not read endpoints)".into(),
                    features: vec![],
                    scale_max,
                    state: field(rest, "state").unwrap_or_else(|| "state not supplied".into()),
                },
            );
        }
        let mut warnings = vec!["Annotated TSS candidate, not experimentally established initiation. Motif features are retained predictions, not dense score traces or measured binding. Signal gaps are not measured zero.".into()];
        let mut tss_found = false;
        for (feature_id, f) in dna.features().iter().enumerate() {
            if f.kind.as_ref() == "source" {
                continue;
            }
            let label = qualifier(f, "label");
            let note = qualifier(f, "note");
            let Some((s, e, reverse, clipped)) = bounds(&f.location) else {
                warnings.push(format!("Feature {feature_id} ({label}) has an unsupported compound/remote location; inspect it in the standard map"));
                continue;
            };
            if s >= e || e > dna.len() {
                return Err(format!(
                    "Feature {feature_id} lies outside the TSS sequence"
                ));
            }
            validate_source_span(&note, &geometry, s, e)
                .map_err(|error| format!("Feature {feature_id}: {error}"))?;
            if label == "Annotated TSS candidate" {
                if s != geometry.upstream_bp
                    || e != s + 1
                    || reverse
                    || note.split(';').next()
                        != Some(format!("Genomic {}", geometry.tss_1based).as_str())
                {
                    return Err("TSS marker contradicts stored genomic geometry".into());
                }
                tss_found = true;
                continue;
            }
            let (kind, id, lane_label, units) = if note.starts_with("Predicted motif;") {
                let accession = required(&note, "accession")?;
                let units = required(&note, "score_kind")?;
                let expected = if reverse { "-" } else { "+" };
                let genomic = if reverse { strand.opposite() } else { strand };
                if required(&note, "motif_local_strand")? != expected
                    || required(&note, "genomic_strand")? != genomic.as_str()
                    || number(&note, "genomic_window_start")? != geometry.genomic_at(s).unwrap()
                {
                    return Err(format!(
                        "Motif feature {feature_id} contradicts its strand/coordinate metadata"
                    ));
                }
                (
                    TssLaneKind::Motif,
                    format!("{accession}/{units}"),
                    accession,
                    units,
                )
            } else if note.contains("; signal interval, not an individual read;") {
                let id = format!("{}/{}", required(&note, "group")?, required(&note, "lane")?);
                if !lanes.contains_key(&(TssLaneKind::Signal, id.clone())) {
                    warnings.push(format!("Signal source {id} lacks lane metadata; scaling uses supplied interval magnitudes"));
                }
                (
                    TssLaneKind::Signal,
                    id.clone(),
                    id,
                    "source signal (not read endpoints)".into(),
                )
            } else if f.kind.as_ref() == "exon"
                || label.contains(" CDS segment")
                || label.contains(" translation ")
            {
                let transcript = qualifier(f, "transcript_id");
                let id = if transcript.is_empty() {
                    label
                        .split_whitespace()
                        .next()
                        .unwrap_or("Structure")
                        .into()
                } else {
                    transcript
                };
                (
                    TssLaneKind::Structure,
                    id.clone(),
                    id,
                    "exon / coding segment / translation marker".into(),
                )
            } else {
                (
                    TssLaneKind::Other,
                    f.kind.to_string(),
                    f.kind.to_string(),
                    "supplied annotation".into(),
                )
            };
            let score = finite_score(&note);
            let lane = lanes
                .entry((kind, id.clone()))
                .or_insert_with(|| TssViewLane {
                    kind,
                    id,
                    label: lane_label,
                    details: String::new(),
                    units,
                    features: vec![],
                    scale_max: 1.0,
                    state: "supplied annotations".into(),
                });
            if let Some(score) = score {
                lane.scale_max = lane.scale_max.max(score.abs());
            }
            lane.features.push(TssViewFeature {
                feature_id,
                start: s,
                end: e,
                reverse,
                clipped,
                label,
                details: note,
                score,
            });
        }
        if !tss_found {
            return Err(
                "Annotated TSS point feature missing; stored metadata alone is insufficient".into(),
            );
        }
        for lane in lanes.values_mut() {
            lane.features
                .sort_by_key(|f| (f.start, f.end, f.feature_id));
        }
        Ok(Self {
            title: format!(
                "{gene} | TSS {}:{} ({}) | {assembly}",
                geometry.chromosome,
                geometry.tss_1based,
                strand.as_str()
            ),
            promoter_id,
            assembly,
            geometry,
            sequence_sha256,
            lanes: lanes.into_values().collect(),
            provenance,
            warnings,
        })
    }

    /// Three coordinate frames for the same base; never reverse an already oriented sequence.
    pub fn coordinate_label(&self, local: usize) -> String {
        match self.geometry.genomic_at(local) {
            Some(genomic) => format!(
                "local {} | TSS {:+} bp | {}:{} ({})",
                local + 1,
                local as i64 - self.geometry.upstream_bp as i64,
                self.geometry.chromosome,
                genomic,
                self.geometry.strand.as_str()
            ),
            None => "outside sequence".into(),
        }
    }
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;
    use gb_io::seq::Seq;

    // Synthetic nine-base exported-window grammar; no public/private experimental data.
    pub(crate) fn fixture(minus: bool) -> DNAsequence {
        let mut seq = Seq::empty();
        seq.seq = b"ACGTACGTA".to_vec();
        let strand = if minus { "-" } else { "+" };
        let tss = 104;
        seq.comments = vec![format!("GENtle promoter_id=toy; sequence_sha256={}", sha256_hex_bytes(&seq.seq)),
            format!("Reference=toy; assembly=toy1; chromosome=1; genomic=100..108; genomic_strand={strand}; local_axis=transcript_5prime_to_3prime; TSS_local_1based=5"),
            "Signal lane group/sample: TP73 treatment; state=available; source_id=synthetic; inherited_abs_max=8".into(),
            "Signal lane group/missing: control; state=not_prepared; source_id=missing; inherited_abs_max=1".into()];
        let feature = |kind: &str, s, e, label: &str, note: String| Feature {
            kind: kind.to_owned().into(),
            location: Location::simple_range(s, e),
            qualifiers: vec![
                ("label".into(), Some(label.into())),
                ("note".into(), Some(note)),
            ],
        };
        seq.features = vec![
            feature("source", 0, 9, "TOY", String::new()),
            feature(
                "misc_feature",
                4,
                5,
                "Annotated TSS candidate",
                format!(
                    "Genomic {tss}; annotation-derived, not experimentally established initiation"
                ),
            ),
            feature(
                "exon",
                4,
                9,
                "tx E1",
                format!(
                    "Source genomic interval {} (1-based inclusive); cropped_to_window=true",
                    if minus { "98..104" } else { "104..110" }
                ),
            ),
            feature(
                "misc_feature",
                0,
                2,
                "signal",
                format!(
                    "Source genomic interval {}; signal interval, not an individual read; group=group; lane=sample; source_id=synthetic; raw_score=2",
                    if minus { "107..108" } else { "100..101" }
                ),
            ),
            feature(
                "misc_feature",
                6,
                8,
                "signal",
                format!(
                    "Source genomic interval {}; signal interval, not an individual read; group=group; lane=sample; source_id=synthetic; raw_score=4",
                    if minus { "101..102" } else { "106..107" }
                ),
            ),
            feature(
                "misc_feature",
                0,
                3,
                "MAtoy.1 stored score peak",
                format!(
                    "Predicted motif; accession=MAtoy.1; raw_score=3; score_kind=llr_background_tail_log10; motif_local_strand=+; genomic_strand={strand}; genomic_window_start={}; no new peak calling",
                    if minus { 108 } else { 100 }
                ),
            ),
        ];
        DNAsequence::from_genbank_seq(seq)
    }

    #[test]
    fn tss_sequence_tracks_preserve_gaps_units_and_missing_sources() {
        let view = TssSequenceView::from_dna(&fixture(false)).unwrap();
        assert!(view.title.contains("TSS 1:104 (+)"));
        let signal = view.lanes.iter().find(|l| l.id == "group/sample").unwrap();
        assert_eq!(
            signal
                .features
                .iter()
                .map(|f| (f.start, f.end))
                .collect::<Vec<_>>(),
            [(0, 2), (6, 8)]
        );
        assert_eq!(signal.scale_max, 8.0);
        assert!(
            view.lanes
                .iter()
                .any(|l| l.id == "group/missing" && l.features.is_empty())
        );
        assert!(
            view.lanes
                .iter()
                .any(|l| l.units == "llr_background_tail_log10")
        );
    }

    #[test]
    fn tss_sequence_negative_axis_and_stale_bases_fail_closed() {
        let dna = fixture(true);
        let view = TssSequenceView::from_dna(&dna).unwrap();
        assert!(view.coordinate_label(0).contains("1:108 (-)"));
        assert!(view.coordinate_label(8).contains("1:100 (-)"));
        assert!(view.coordinate_label(4).contains("TSS +0 bp"));
        let mut seq = dna.clone_seq_record();
        seq.seq[0] = b'T';
        assert!(
            TssSequenceView::from_dna(&DNAsequence::from_genbank_seq(seq))
                .unwrap_err()
                .contains("hash mismatch")
        );
    }

    #[test]
    fn tss_sequence_rejects_missing_and_contradictory_markers() {
        let mut seq = fixture(false).clone_seq_record();
        seq.features[1].location = Location::single(2);
        assert!(TssSequenceView::from_dna(&DNAsequence::from_genbank_seq(seq)).is_err());
        assert!(TssSequenceView::from_dna(&DNAsequence::from_sequence("ACGT").unwrap()).is_err());
    }

    #[test]
    fn tss_sequence_groups_by_matrix_and_retains_negative_scores_without_rescoring() {
        let mut seq = fixture(true).clone_seq_record();
        let mut peak = seq.features.last().unwrap().clone();
        peak.location = Location::Complement(Box::new(Location::simple_range(0, 3)));
        peak.qualifiers = vec![("label".into(), Some("MAtoy.1 stored score peak".into())),
            ("note".into(), Some("Predicted motif; accession=MAtoy.1; raw_score=-2; score_kind=llr_background_tail_log10; motif_local_strand=-; genomic_strand=+; genomic_window_start=108; no new peak calling".into()))];
        seq.features.push(peak);
        let view = TssSequenceView::from_dna(&DNAsequence::from_genbank_seq(seq)).unwrap();
        let lane = view
            .lanes
            .iter()
            .find(|l| l.kind == TssLaneKind::Motif)
            .unwrap();
        assert_eq!(lane.features.len(), 2);
        assert!(
            lane.features
                .iter()
                .any(|f| f.reverse && f.score == Some(-2.0))
        );
        assert!(
            lane.features
                .iter()
                .any(|f| !f.reverse && f.score == Some(3.0))
        );
    }

    #[test]
    #[ignore = "explicit local-file inspection; requires GENTLE_TSS_INSPECT_PATH, no private fixture committed"]
    fn tss_sequence_inspect_local_embl() {
        let path =
            std::env::var("GENTLE_TSS_INSPECT_PATH").expect("provide local annotated TSS EMBL");
        let records = DNAsequence::from_embl_file(&path).unwrap();
        assert!(!records.is_empty());
        for dna in records {
            let view = TssSequenceView::from_dna(&dna).unwrap();
            println!(
                "{}: {} grouped lanes, {} annotated intervals",
                view.title,
                view.lanes.len(),
                view.lanes.iter().map(|l| l.features.len()).sum::<usize>()
            );
        }
    }

    #[test]
    fn tss_sequence_embl_and_genbank_roundtrip_tracks() {
        for minus in [false, true] {
            let dna = fixture(minus);
            let dir = tempfile::tempdir().unwrap();
            let embl = dir.path().join("toy.embl");
            std::fs::write(
                &embl,
                crate::annotated_sequence_io::embl_bytes(&dna.clone_seq_record()).unwrap(),
            )
            .unwrap();
            let imported = DNAsequence::from_embl_file(embl.to_str().unwrap()).unwrap();
            let gb = dir.path().join("toy.gb");
            std::fs::write(
                &gb,
                crate::annotated_sequence_io::genbank_bytes(&dna.clone_seq_record()).unwrap(),
            )
            .unwrap();
            let genbank = DNAsequence::from_genbank_file(gb.to_str().unwrap()).unwrap();
            let expected = TssSequenceView::from_dna(&dna).unwrap();
            for record in [&imported[0], &genbank[0]] {
                let view = TssSequenceView::from_dna(record).unwrap();
                assert_eq!(view.geometry, expected.geometry);
                assert_eq!(
                    view.lanes.iter().map(|l| l.features.len()).sum::<usize>(),
                    4
                );
            }
        }
    }
}
