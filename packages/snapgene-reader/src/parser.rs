use crate::{
    AdditionalSequenceProperties, DnaPacket, FeatureRecord, HybridizationParams, NotesRecord,
    PacketType, PrimerRecord, RawPacket, SequenceTopology, SnapGeneError, SnapGeneFile, xml,
};
use std::{
    fs,
    io::Read,
    path::Path,
};

pub fn parse_path<P: AsRef<Path>>(path: P) -> Result<SnapGeneFile, SnapGeneError> {
    let path_ref = path.as_ref();
    let bytes = fs::read(path_ref).map_err(|source| SnapGeneError::Io {
        path: path_ref.display().to_string(),
        source,
    })?;
    parse_bytes(&bytes)
}

pub fn parse_reader<R: Read>(mut reader: R) -> Result<SnapGeneFile, SnapGeneError> {
    let mut bytes = Vec::new();
    reader
        .read_to_end(&mut bytes)
        .map_err(|source| SnapGeneError::Io {
            path: "<reader>".to_string(),
            source,
        })?;
    parse_bytes(&bytes)
}

pub fn parse_bytes(bytes: &[u8]) -> Result<SnapGeneFile, SnapGeneError> {
    if bytes.is_empty() {
        return Err(SnapGeneError::EmptyFile);
    }

    let mut offset = 0usize;
    let mut saw_cookie = false;
    let mut sequence_kind = 0u16;
    let mut export_version = 0u16;
    let mut import_version = 0u16;
    let mut dna: Option<DnaPacket> = None;
    let mut features: Vec<FeatureRecord> = vec![];
    let mut primers: Vec<PrimerRecord> = vec![];
    let mut primer_hybridization_params: Option<HybridizationParams> = None;
    let mut notes: Option<NotesRecord> = None;
    let mut additional_sequence_properties: Option<AdditionalSequenceProperties> = None;
    let mut extra_packets: Vec<RawPacket> = vec![];

    while offset < bytes.len() {
        if bytes.len().saturating_sub(offset) < 5 {
            return Err(SnapGeneError::UnexpectedEof {
                context: format!("packet header at byte offset {offset}"),
            });
        }
        let packet_type = PacketType::from_byte(bytes[offset]);
        offset += 1;
        let length = u32::from_be_bytes(
            bytes[offset..offset + 4]
                .try_into()
                .expect("4-byte packet length"),
        ) as usize;
        offset += 4;
        if bytes.len().saturating_sub(offset) < length {
            return Err(SnapGeneError::UnexpectedEof {
                context: format!("packet payload for {packet_type:?}"),
            });
        }
        let data = &bytes[offset..offset + length];
        offset += length;

        if !saw_cookie {
            if packet_type != PacketType::Cookie {
                return Err(SnapGeneError::MissingCookie { found: packet_type });
            }
            let (kind, export, import) = parse_cookie_packet(data)?;
            sequence_kind = kind;
            export_version = export;
            import_version = import;
            saw_cookie = true;
            continue;
        }

        match packet_type {
            PacketType::Dna => {
                if dna.is_some() {
                    return Err(SnapGeneError::DuplicateDnaPacket);
                }
                dna = Some(parse_dna_packet(data)?);
            }
            PacketType::Features => {
                let xml_text = String::from_utf8(data.to_vec()).map_err(|source| {
                    SnapGeneError::InvalidUtf8 {
                        packet_type,
                        source,
                    }
                })?;
                features.extend(xml::parse_features_packet(&xml_text)?);
            }
            PacketType::Primers => {
                let xml_text = String::from_utf8(data.to_vec()).map_err(|source| {
                    SnapGeneError::InvalidUtf8 {
                        packet_type,
                        source,
                    }
                })?;
                let (params, parsed_primers) = xml::parse_primers_packet(&xml_text)?;
                if params.is_some() {
                    primer_hybridization_params = params;
                }
                primers.extend(parsed_primers);
            }
            PacketType::Notes => {
                let xml_text = String::from_utf8(data.to_vec()).map_err(|source| {
                    SnapGeneError::InvalidUtf8 {
                        packet_type,
                        source,
                    }
                })?;
                notes = Some(xml::parse_notes_packet(&xml_text)?);
            }
            PacketType::AdditionalSequenceProperties => {
                let xml_text = String::from_utf8(data.to_vec()).map_err(|source| {
                    SnapGeneError::InvalidUtf8 {
                        packet_type,
                        source,
                    }
                })?;
                additional_sequence_properties =
                    Some(xml::parse_additional_sequence_properties_packet(&xml_text)?);
            }
            PacketType::Cookie => {
                extra_packets.push(RawPacket {
                    packet_type,
                    payload: data.to_vec(),
                });
            }
            PacketType::Unknown(_) => {
                extra_packets.push(RawPacket {
                    packet_type,
                    payload: data.to_vec(),
                });
            }
        }
    }

    let dna = dna.ok_or(SnapGeneError::MissingDnaPacket)?;

    Ok(SnapGeneFile {
        sequence_kind,
        export_version,
        import_version,
        dna,
        features,
        primers,
        primer_hybridization_params,
        notes,
        additional_sequence_properties,
        extra_packets,
    })
}

fn parse_cookie_packet(data: &[u8]) -> Result<(u16, u16, u16), SnapGeneError> {
    if data.len() < 14 {
        return Err(SnapGeneError::ShortCookie { actual: data.len() });
    }
    let cookie = String::from_utf8_lossy(&data[..8]).into_owned();
    if cookie != "SnapGene" {
        return Err(SnapGeneError::InvalidCookie { found: cookie });
    }
    let sequence_kind = u16::from_be_bytes([data[8], data[9]]);
    let export_version = u16::from_be_bytes([data[10], data[11]]);
    let import_version = u16::from_be_bytes([data[12], data[13]]);
    Ok((sequence_kind, export_version, import_version))
}

fn parse_dna_packet(data: &[u8]) -> Result<DnaPacket, SnapGeneError> {
    if data.is_empty() {
        return Err(SnapGeneError::UnexpectedEof {
            context: "DNA packet flags byte".to_string(),
        });
    }
    let flags = data[0];
    let sequence = String::from_utf8(data[1..].to_vec()).map_err(|source| {
        SnapGeneError::InvalidUtf8 {
            packet_type: PacketType::Dna,
            source,
        }
    })?;
    let topology = if flags & 0x01 == 0x01 {
        SequenceTopology::Circular
    } else {
        SequenceTopology::Linear
    };
    Ok(DnaPacket {
        flags,
        sequence,
        topology,
    })
}

