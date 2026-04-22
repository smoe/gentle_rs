use crate::{
    AdditionalSequenceProperties, FeatureDirectionality, FeatureRecord, FeatureSegment,
    HybridizationParams, NotesRecord, PrimerBindingSite, PrimerRecord, PrimerStrand,
    QualifierRecord, QualifierValue, SnapGeneError,
};
use quick_xml::{
    Reader,
    events::{BytesStart, Event},
};
use std::collections::BTreeMap;

#[derive(Clone, Debug, Default)]
struct XmlNode {
    name: String,
    attributes: BTreeMap<String, String>,
    text: String,
    children: Vec<XmlNode>,
}

impl XmlNode {
    fn child(&self, name: &str) -> Option<&XmlNode> {
        self.children.iter().find(|child| child.name == name)
    }

    fn children_named<'a>(&'a self, name: &'a str) -> impl Iterator<Item = &'a XmlNode> + 'a {
        self.children.iter().filter(move |child| child.name == name)
    }

    fn text_trimmed(&self) -> Option<String> {
        let sanitized = sanitize_text(&self.text);
        let text = sanitized.trim();
        (!text.is_empty()).then_some(text.to_string())
    }
}

pub(crate) fn parse_features_packet(xml: &str) -> Result<Vec<FeatureRecord>, SnapGeneError> {
    let root = parse_xml_document(xml, "features")?;
    let mut features = Vec::new();
    for feature_node in root.children_named("Feature") {
        let name = get_optional_attr(feature_node, "name");
        let feature_type =
            get_optional_attr(feature_node, "type").unwrap_or_else(|| "misc_feature".to_string());
        let directionality = match get_optional_attr(feature_node, "directionality")
            .as_deref()
            .map(str::trim)
        {
            Some("1") => FeatureDirectionality::Forward,
            Some("2") => FeatureDirectionality::Reverse,
            _ => FeatureDirectionality::None,
        };
        let display_name = name.clone().unwrap_or_else(|| feature_type.clone());

        let mut segments = Vec::new();
        for segment_node in feature_node.children_named("Segment") {
            let range = get_optional_attr(segment_node, "range")
                .ok_or_else(|| SnapGeneError::MissingFeatureRange {
                    feature_name: display_name.clone(),
                })?;
            let (start_1based, end_1based) = parse_range_attr(&range)?;
            let translated = get_optional_attr(segment_node, "translated")
                .map(|value| value == "1")
                .unwrap_or(false);
            segments.push(FeatureSegment {
                start_1based,
                end_1based,
                name: get_optional_attr(segment_node, "name"),
                color: get_optional_attr(segment_node, "color"),
                segment_type: get_optional_attr(segment_node, "type"),
                translated,
            });
        }

        let mut qualifiers = Vec::new();
        for qualifier_node in feature_node.children_named("Q") {
            let Some(qualifier_name) = get_optional_attr(qualifier_node, "name") else {
                continue;
            };
            let mut values = Vec::new();
            for value_node in qualifier_node.children_named("V") {
                values.push(QualifierValue {
                    text: get_optional_attr(value_node, "text"),
                    predefined: get_optional_attr(value_node, "predef"),
                    integer: get_optional_attr(value_node, "int")
                        .and_then(|value| value.parse::<i64>().ok()),
                });
            }
            qualifiers.push(QualifierRecord {
                name: qualifier_name,
                values,
            });
        }

        features.push(FeatureRecord {
            name,
            feature_type,
            directionality,
            segments,
            qualifiers,
        });
    }
    Ok(features)
}

pub(crate) fn parse_primers_packet(
    xml: &str,
) -> Result<(Option<HybridizationParams>, Vec<PrimerRecord>), SnapGeneError> {
    let root = parse_xml_document(xml, "primers")?;
    let params = root.child("HybridizationParams").map(parse_hybridization_params);
    let mut primers = Vec::new();
    for primer_node in root.children_named("Primer") {
        let mut binding_sites = Vec::new();
        for site_node in primer_node.children_named("BindingSite") {
            let range = get_optional_attr(site_node, "location").ok_or_else(|| {
                SnapGeneError::InvalidRange {
                    range: "<missing primer location>".to_string(),
                }
            })?;
            let (location_start_1based, location_end_1based) = parse_range_attr(&range)?;
            let bound_strand = match get_optional_attr(site_node, "boundStrand")
                .as_deref()
                .unwrap_or("0")
            {
                "1" => PrimerStrand::Reverse,
                _ => PrimerStrand::Forward,
            };
            let simplified = get_optional_attr(site_node, "simplified")
                .map(|value| value == "1")
                .unwrap_or(false);
            let melting_temperature_c = get_optional_attr(site_node, "meltingTemperature")
                .and_then(|value| value.parse::<f32>().ok());
            binding_sites.push(PrimerBindingSite {
                location_start_1based,
                location_end_1based,
                bound_strand,
                simplified,
                annealed_bases: get_optional_attr(site_node, "annealedBases"),
                melting_temperature_c,
            });
        }
        primers.push(PrimerRecord {
            name: get_optional_attr(primer_node, "name"),
            sequence: get_optional_attr(primer_node, "sequence"),
            description: get_optional_attr(primer_node, "description"),
            binding_sites,
        });
    }
    Ok((params, primers))
}

pub(crate) fn parse_notes_packet(xml: &str) -> Result<NotesRecord, SnapGeneError> {
    let root = parse_xml_document(xml, "notes")?;
    Ok(NotesRecord {
        uuid: child_text(&root, "UUID"),
        sequence_type: child_text(&root, "Type"),
        confirmed_experimentally: child_text(&root, "ConfirmedExperimentally")
            .map(|value| value == "1" || value.eq_ignore_ascii_case("true")),
        accession_number: child_text(&root, "AccessionNumber"),
        custom_map_label: child_text(&root, "CustomMapLabel"),
        use_custom_map_label: child_text(&root, "UseCustomMapLabel")
            .map(|value| value == "1" || value.eq_ignore_ascii_case("true")),
        description: child_text(&root, "Description"),
        created_on: root.child("Created").and_then(XmlNode::text_trimmed),
        created_utc: root
            .child("Created")
            .and_then(|node| get_optional_attr(node, "UTC")),
        last_modified_on: root.child("LastModified").and_then(XmlNode::text_trimmed),
        last_modified_utc: root
            .child("LastModified")
            .and_then(|node| get_optional_attr(node, "UTC")),
        created_by: child_text(&root, "CreatedBy"),
        organism: child_text(&root, "Organism"),
        sequence_class: child_text(&root, "SequenceClass"),
        transformed_into: child_text(&root, "TransformedInto"),
        comments: child_text(&root, "Comments"),
    })
}

pub(crate) fn parse_additional_sequence_properties_packet(
    xml: &str,
) -> Result<AdditionalSequenceProperties, SnapGeneError> {
    let root = parse_xml_document(xml, "additional sequence properties")?;
    Ok(AdditionalSequenceProperties {
        upstream_stickiness: child_text(&root, "UpstreamStickiness")
            .and_then(|value| value.parse::<i32>().ok()),
        downstream_stickiness: child_text(&root, "DownstreamStickiness")
            .and_then(|value| value.parse::<i32>().ok()),
        upstream_modification: child_text(&root, "UpstreamModification"),
        downstream_modification: child_text(&root, "DownstreamModification"),
    })
}

fn parse_hybridization_params(node: &XmlNode) -> HybridizationParams {
    HybridizationParams {
        min_continuous_match_len: get_optional_attr(node, "minContinuousMatchLen")
            .and_then(|value| value.parse::<usize>().ok())
            .unwrap_or(0),
        min_melting_temperature_c: get_optional_attr(node, "minMeltingTemperature")
            .and_then(|value| value.parse::<f32>().ok())
            .unwrap_or(0.0),
        allow_mismatch: get_optional_attr(node, "allowMismatch")
            .map(|value| value == "1" || value.eq_ignore_ascii_case("true"))
            .unwrap_or(false),
        show_additional_five_prime_matches: get_optional_attr(
            node,
            "showAdditionalFivePrimeMatches",
        )
        .map(|value| value == "1" || value.eq_ignore_ascii_case("true"))
        .unwrap_or(false),
        minimum_five_prime_annealing: get_optional_attr(node, "minimumFivePrimeAnnealing")
            .and_then(|value| value.parse::<usize>().ok()),
    }
}

fn parse_range_attr(range: &str) -> Result<(usize, usize), SnapGeneError> {
    let Some((start, end)) = range.split_once('-') else {
        return Err(SnapGeneError::InvalidRange {
            range: range.to_string(),
        });
    };
    let start_1based = start
        .trim()
        .parse::<usize>()
        .map_err(|_| SnapGeneError::InvalidRange {
            range: range.to_string(),
        })?;
    let end_1based = end
        .trim()
        .parse::<usize>()
        .map_err(|_| SnapGeneError::InvalidRange {
            range: range.to_string(),
        })?;
    if start_1based == 0 || end_1based == 0 {
        return Err(SnapGeneError::InvalidRange {
            range: range.to_string(),
        });
    }
    Ok((start_1based, end_1based))
}

fn child_text(node: &XmlNode, name: &str) -> Option<String> {
    node.child(name).and_then(XmlNode::text_trimmed)
}

fn get_optional_attr(node: &XmlNode, name: &str) -> Option<String> {
    node.attributes
        .get(name)
        .map(|value| sanitize_text(value))
        .filter(|value| !value.is_empty())
}

fn parse_xml_document(xml: &str, context: &'static str) -> Result<XmlNode, SnapGeneError> {
    let mut reader = Reader::from_str(xml);
    reader.config_mut().trim_text(false);
    let mut stack: Vec<XmlNode> = Vec::new();
    let mut root: Option<XmlNode> = None;

    loop {
        match reader.read_event() {
            Ok(Event::Start(start)) => {
                stack.push(start_node(&reader, &start, context)?);
            }
            Ok(Event::Empty(start)) => {
                let node = start_node(&reader, &start, context)?;
                if let Some(parent) = stack.last_mut() {
                    parent.children.push(node);
                } else {
                    root = Some(node);
                }
            }
            Ok(Event::End(_)) => {
                let Some(node) = stack.pop() else {
                    return Err(SnapGeneError::Xml {
                        context,
                        message: "unexpected end tag".to_string(),
                    });
                };
                if let Some(parent) = stack.last_mut() {
                    parent.children.push(node);
                } else {
                    root = Some(node);
                }
            }
            Ok(Event::Text(text)) => {
                if let Some(node) = stack.last_mut() {
                    let decoded = text.unescape().map_err(|err| SnapGeneError::Xml {
                        context,
                        message: err.to_string(),
                    })?;
                    node.text.push_str(decoded.as_ref());
                }
            }
            Ok(Event::CData(text)) => {
                if let Some(node) = stack.last_mut() {
                    let decoded = text.decode().map_err(|err| SnapGeneError::Xml {
                        context,
                        message: err.to_string(),
                    })?;
                    node.text.push_str(decoded.as_ref());
                }
            }
            Ok(Event::Decl(_))
            | Ok(Event::Comment(_))
            | Ok(Event::DocType(_))
            | Ok(Event::PI(_)) => {}
            Ok(Event::Eof) => break,
            Err(err) => {
                return Err(SnapGeneError::Xml {
                    context,
                    message: err.to_string(),
                });
            }
        }
    }

    root.ok_or(SnapGeneError::Xml {
        context,
        message: "missing XML root node".to_string(),
    })
}

fn start_node(
    reader: &Reader<&[u8]>,
    start: &BytesStart<'_>,
    context: &'static str,
) -> Result<XmlNode, SnapGeneError> {
    let name = String::from_utf8_lossy(start.name().as_ref()).into_owned();
    let mut attributes = BTreeMap::new();
    for attr in start.attributes().with_checks(false) {
        let attr = attr.map_err(|err| SnapGeneError::Xml {
            context,
            message: err.to_string(),
        })?;
        let key = String::from_utf8_lossy(attr.key.as_ref()).into_owned();
        let value = attr
            .decode_and_unescape_value(reader.decoder())
            .map_err(|err| SnapGeneError::Xml {
                context,
                message: err.to_string(),
            })?;
        attributes.insert(key, value.into_owned());
    }
    Ok(XmlNode {
        name,
        attributes,
        text: String::new(),
        children: vec![],
    })
}

fn sanitize_text(raw: &str) -> String {
    let mut out = String::with_capacity(raw.len());
    let mut inside_tag = false;
    for ch in raw.chars() {
        match ch {
            '<' => inside_tag = true,
            '>' => inside_tag = false,
            _ if inside_tag => {}
            '\r' | '\n' | '\t' => out.push(' '),
            '\u{a0}' => out.push(' '),
            other => out.push(other),
        }
    }
    out.split_whitespace().collect::<Vec<_>>().join(" ")
}
