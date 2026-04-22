use crate::{
    FeatureDirectionality, FeatureRecord, NotesRecord, PrimerRecord, PrimerStrand,
    QualifierValue, SequenceTopology, SnapGeneError, SnapGeneFile,
};
use gb_io::seq::{Feature, Location, Seq, Source, Topology};

impl TryFrom<&SnapGeneFile> for Seq {
    type Error = SnapGeneError;

    fn try_from(file: &SnapGeneFile) -> Result<Self, Self::Error> {
        let mut seq = Seq::empty();
        let notes = file.notes.as_ref();
        seq.name = notes.and_then(NotesRecord::effective_name).map(ToString::to_string);
        seq.definition = notes.and_then(|notes| notes.description.clone());
        seq.accession = notes.and_then(|notes| notes.accession_number.clone());
        seq.source = notes.and_then(|notes| {
            notes.organism.clone().map(|organism| Source {
                source: organism.clone(),
                organism: Some(organism),
            })
        });
        seq.division = notes
            .map(|notes| notes.division_code().to_string())
            .unwrap_or_default();
        seq.topology = match file.dna.topology {
            SequenceTopology::Circular => Topology::Circular,
            SequenceTopology::Linear => Topology::Linear,
        };
        seq.molecule_type = Some("DNA".to_string());
        seq.seq = file.dna.sequence.as_bytes().to_vec();
        seq.len = Some(seq.seq.len());
        if let Some(notes) = notes {
            if let Some(comments) = notes.comments.clone().filter(|value| !value.is_empty()) {
                seq.comments.push(comments);
            }
            if let Some(created_by) = notes.created_by.clone().filter(|value| !value.is_empty()) {
                seq.comments.push(format!("Created by: {created_by}"));
            }
        }

        for feature in &file.features {
            seq.features.push(feature_to_gb(feature, seq.seq.len())?);
        }
        for primer in &file.primers {
            for feature in primer_to_gb(primer, file.primer_hybridization_params.as_ref(), seq.seq.len())? {
                seq.features.push(feature);
            }
        }
        Ok(seq)
    }
}

fn feature_to_gb(feature: &FeatureRecord, sequence_len: usize) -> Result<Feature, SnapGeneError> {
    let mut qualifiers: Vec<gb_io::seq::Qualifier> = vec![];
    for qualifier in &feature.qualifiers {
        for value in &qualifier.values {
            let text = qualifier_value_text(value);
            qualifiers.push((qualifier.name.clone().into(), text));
        }
    }

    if let Some(name) = feature.name.clone().filter(|value| !value.is_empty()) {
        let has_label = qualifiers.iter().any(|(key, _)| key.as_ref() == "label");
        if !has_label {
            qualifiers.push(("label".to_string().into(), Some(name)));
        }
    }

    let location_text = feature_location_text(feature, sequence_len, false);
    let location = Location::from_gb_format(&location_text).map_err(|err| {
        SnapGeneError::InvalidLocation {
            context: "feature",
            location: location_text.clone(),
            message: err.to_string(),
        }
    })?;

    Ok(Feature {
        kind: feature.feature_type.clone().into(),
        location,
        qualifiers,
    })
}

fn primer_to_gb(
    primer: &PrimerRecord,
    params: Option<&crate::HybridizationParams>,
    sequence_len: usize,
) -> Result<Vec<Feature>, SnapGeneError> {
    let mut out = Vec::new();
    for site in primer.visible_binding_sites(params) {
        let strand_is_reverse = matches!(site.bound_strand, PrimerStrand::Reverse);
        let mut location_text = if site.location_start_1based <= site.location_end_1based {
            primer_piece_text(site.location_start_1based, site.location_end_1based, sequence_len)
        } else {
            primer_piece_text(site.location_start_1based, sequence_len, sequence_len)
                + ","
                + &primer_piece_text(1, site.location_end_1based, sequence_len)
        };
        if location_text.contains(',') {
            location_text = format!("join({location_text})");
        }
        if strand_is_reverse {
            location_text = format!("complement({location_text})");
        }
        let location = Location::from_gb_format(&location_text).map_err(|err| {
            SnapGeneError::InvalidLocation {
                context: "primer binding site",
                location: location_text.clone(),
                message: err.to_string(),
            }
        })?;

        let mut qualifiers = Vec::new();
        if let Some(name) = primer.name.clone().filter(|value| !value.is_empty()) {
            qualifiers.push(("label".to_string().into(), Some(name)));
        }
        if let Some(description) = primer.description.clone().filter(|value| !value.is_empty()) {
            qualifiers.push(("note".to_string().into(), Some(description)));
        }
        out.push(Feature {
            kind: "primer_bind".to_string().into(),
            location,
            qualifiers,
        });
    }
    Ok(out)
}

fn feature_location_text(
    feature: &FeatureRecord,
    sequence_len: usize,
    is_primer: bool,
) -> String {
    let mut pieces = Vec::new();
    for segment in &feature.segments {
        if segment.start_1based <= segment.end_1based {
            pieces.push(range_piece_text(
                segment.start_1based,
                segment.end_1based,
                sequence_len,
                is_primer,
            ));
        } else {
            pieces.push(range_piece_text(
                segment.start_1based,
                sequence_len,
                sequence_len,
                is_primer,
            ));
            pieces.push(range_piece_text(1, segment.end_1based, sequence_len, is_primer));
        }
    }

    let mut pieces: Vec<String> = match feature.directionality {
        FeatureDirectionality::Reverse => pieces
            .into_iter()
            .map(|piece| format!("complement({piece})"))
            .collect(),
        _ => pieces,
    };

    if pieces.len() == 1 {
        pieces.pop().unwrap_or_default()
    } else {
        format!("join({})", pieces.join(","))
    }
}

fn range_piece_text(
    start_1based: usize,
    end_1based: usize,
    _sequence_len: usize,
    is_primer: bool,
) -> String {
    let (start, end) = if is_primer {
        (start_1based + 1, end_1based + 1)
    } else {
        (start_1based, end_1based)
    };
    format!("{start}..{end}")
}

fn primer_piece_text(start_1based: usize, end_1based: usize, sequence_len: usize) -> String {
    range_piece_text(start_1based, end_1based, sequence_len, true)
}

fn qualifier_value_text(value: &QualifierValue) -> Option<String> {
    value
        .text
        .clone()
        .or_else(|| value.predefined.clone())
        .or_else(|| value.integer.map(|value| value.to_string()))
}
