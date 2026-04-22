use std::collections::HashSet;

/// Parsed SnapGene file contents.
#[derive(Clone, Debug, PartialEq)]
pub struct SnapGeneFile {
    pub sequence_kind: u16,
    pub export_version: u16,
    pub import_version: u16,
    pub dna: DnaPacket,
    pub features: Vec<FeatureRecord>,
    pub primers: Vec<PrimerRecord>,
    pub primer_hybridization_params: Option<HybridizationParams>,
    pub notes: Option<NotesRecord>,
    pub additional_sequence_properties: Option<AdditionalSequenceProperties>,
    pub extra_packets: Vec<RawPacket>,
}

/// One unparsed packet preserved verbatim from the original file.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RawPacket {
    pub packet_type: PacketType,
    pub payload: Vec<u8>,
}

/// Packet type identifier used by the SnapGene TLV container.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum PacketType {
    Cookie,
    Dna,
    Primers,
    Notes,
    AdditionalSequenceProperties,
    Features,
    Unknown(u8),
}

impl PacketType {
    pub fn from_byte(value: u8) -> Self {
        match value {
            0x09 => Self::Cookie,
            0x00 => Self::Dna,
            0x05 => Self::Primers,
            0x06 => Self::Notes,
            0x08 => Self::AdditionalSequenceProperties,
            0x0A => Self::Features,
            other => Self::Unknown(other),
        }
    }

    pub fn as_byte(self) -> u8 {
        match self {
            Self::Cookie => 0x09,
            Self::Dna => 0x00,
            Self::Primers => 0x05,
            Self::Notes => 0x06,
            Self::AdditionalSequenceProperties => 0x08,
            Self::Features => 0x0A,
            Self::Unknown(other) => other,
        }
    }
}

/// Topology encoded in the DNA packet.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SequenceTopology {
    Linear,
    Circular,
}

/// Sequence payload from the SnapGene DNA packet.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DnaPacket {
    pub flags: u8,
    pub sequence: String,
    pub topology: SequenceTopology,
}

/// SnapGene feature directionality.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FeatureDirectionality {
    None,
    Forward,
    Reverse,
}

/// One feature entry from the SnapGene Features XML packet.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FeatureRecord {
    pub name: Option<String>,
    pub feature_type: String,
    pub directionality: FeatureDirectionality,
    pub segments: Vec<FeatureSegment>,
    pub qualifiers: Vec<QualifierRecord>,
}

/// One feature segment from a SnapGene feature record.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FeatureSegment {
    pub start_1based: usize,
    pub end_1based: usize,
    pub name: Option<String>,
    pub color: Option<String>,
    pub segment_type: Option<String>,
    pub translated: bool,
}

/// One named qualifier on a feature entry.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct QualifierRecord {
    pub name: String,
    pub values: Vec<QualifierValue>,
}

/// One typed qualifier value from SnapGene XML.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct QualifierValue {
    pub text: Option<String>,
    pub predefined: Option<String>,
    pub integer: Option<i64>,
}

/// Hybridization thresholds that determine which primer binding sites SnapGene shows.
#[derive(Clone, Debug, PartialEq)]
pub struct HybridizationParams {
    pub min_continuous_match_len: usize,
    pub min_melting_temperature_c: f32,
    pub allow_mismatch: bool,
    pub show_additional_five_prime_matches: bool,
    pub minimum_five_prime_annealing: Option<usize>,
}

/// Strand a primer binds to in a SnapGene primer binding site.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum PrimerStrand {
    Forward,
    Reverse,
}

/// One primer record from the Primers XML packet.
#[derive(Clone, Debug, PartialEq)]
pub struct PrimerRecord {
    pub name: Option<String>,
    pub sequence: Option<String>,
    pub description: Option<String>,
    pub binding_sites: Vec<PrimerBindingSite>,
}

impl PrimerRecord {
    /// Return the primer binding sites that SnapGene would normally display.
    pub fn visible_binding_sites<'a>(
        &'a self,
        params: Option<&HybridizationParams>,
    ) -> Vec<&'a PrimerBindingSite> {
        let mut seen_simplified: HashSet<(usize, usize, PrimerStrand)> = HashSet::new();
        let mut visible = Vec::new();
        for site in &self.binding_sites {
            if let Some(params) = params {
                if let Some(annealed) = site.annealed_bases.as_ref() {
                    if annealed.len() < params.min_continuous_match_len {
                        continue;
                    }
                }
                if let Some(temp) = site.melting_temperature_c {
                    if temp < params.min_melting_temperature_c {
                        continue;
                    }
                }
            }
            let key = (
                site.location_start_1based,
                site.location_end_1based,
                site.bound_strand,
            );
            if site.simplified && !seen_simplified.insert(key) {
                continue;
            }
            visible.push(site);
        }
        visible
    }
}

/// One primer binding site from a SnapGene primer entry.
#[derive(Clone, Debug, PartialEq)]
pub struct PrimerBindingSite {
    pub location_start_1based: usize,
    pub location_end_1based: usize,
    pub bound_strand: PrimerStrand,
    pub simplified: bool,
    pub annealed_bases: Option<String>,
    pub melting_temperature_c: Option<f32>,
}

/// Structured contents of the SnapGene `Notes` packet.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct NotesRecord {
    pub uuid: Option<String>,
    pub sequence_type: Option<String>,
    pub confirmed_experimentally: Option<bool>,
    pub accession_number: Option<String>,
    pub custom_map_label: Option<String>,
    pub use_custom_map_label: Option<bool>,
    pub description: Option<String>,
    pub created_on: Option<String>,
    pub created_utc: Option<String>,
    pub last_modified_on: Option<String>,
    pub last_modified_utc: Option<String>,
    pub created_by: Option<String>,
    pub organism: Option<String>,
    pub sequence_class: Option<String>,
    pub transformed_into: Option<String>,
    pub comments: Option<String>,
}

impl NotesRecord {
    pub fn effective_name(&self) -> Option<&str> {
        if self.use_custom_map_label.unwrap_or(false) {
            return self.custom_map_label.as_deref();
        }
        self.accession_number.as_deref()
    }

    pub fn division_code(&self) -> &'static str {
        match self.sequence_type.as_deref() {
            Some(kind) if kind.eq_ignore_ascii_case("synthetic") => "SYN",
            _ => "UNC",
        }
    }
}

/// Structured contents of the SnapGene `AdditionalSequenceProperties` packet.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct AdditionalSequenceProperties {
    pub upstream_stickiness: Option<i32>,
    pub downstream_stickiness: Option<i32>,
    pub upstream_modification: Option<String>,
    pub downstream_modification: Option<String>,
}
