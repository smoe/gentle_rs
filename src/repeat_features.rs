//! Repeat-feature classification helpers for display and export.
//!
//! RepeatMasker/UCSC `rmsk` data commonly arrives as ordinary sequence
//! features with `repName` / `repClass` / `repFamily` style qualifiers.  This
//! module keeps the subtype taxonomy outside GUI-only rendering code so live
//! maps, feature trees, and SVG exports can interpret those records the same
//! way.

use gb_io::seq::Feature;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum RepeatDisplayClass {
    Sine,
    Line,
    Ltr,
    DnaTransposon,
    SimpleRepeat,
    LowComplexity,
    Satellite,
    Rna,
    Retroposon,
    RollingCircle,
    Unknown,
    Other,
}

impl RepeatDisplayClass {
    pub fn label(self) -> &'static str {
        match self {
            Self::Sine => "SINE",
            Self::Line => "LINE",
            Self::Ltr => "LTR",
            Self::DnaTransposon => "DNA transposon",
            Self::SimpleRepeat => "Simple repeat",
            Self::LowComplexity => "Low complexity",
            Self::Satellite => "Satellite",
            Self::Rna => "RNA repeat",
            Self::Retroposon => "Retroposon",
            Self::RollingCircle => "Rolling-circle",
            Self::Unknown => "Unknown repeat",
            Self::Other => "Repeat",
        }
    }

    pub fn color_rgb(self) -> (u8, u8, u8) {
        match self {
            Self::Sine => (14, 116, 144),
            Self::Line => (37, 99, 235),
            Self::Ltr => (124, 58, 237),
            Self::DnaTransposon => (234, 88, 12),
            Self::SimpleRepeat => (22, 163, 74),
            Self::LowComplexity => (202, 138, 4),
            Self::Satellite => (219, 39, 119),
            Self::Rna => (8, 145, 178),
            Self::Retroposon => (147, 51, 234),
            Self::RollingCircle => (5, 150, 105),
            Self::Unknown => (100, 116, 139),
            Self::Other => (71, 85, 105),
        }
    }

    pub fn color_hex(self) -> &'static str {
        match self {
            Self::Sine => "#0e7490",
            Self::Line => "#2563eb",
            Self::Ltr => "#7c3aed",
            Self::DnaTransposon => "#ea580c",
            Self::SimpleRepeat => "#16a34a",
            Self::LowComplexity => "#ca8a04",
            Self::Satellite => "#db2777",
            Self::Rna => "#0891b2",
            Self::Retroposon => "#9333ea",
            Self::RollingCircle => "#059669",
            Self::Unknown => "#64748b",
            Self::Other => "#475569",
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RepeatFeatureDisplay {
    pub name: Option<String>,
    pub class: RepeatDisplayClass,
    pub class_label: String,
    pub family: Option<String>,
}

impl RepeatFeatureDisplay {
    pub fn class_family_label(&self) -> String {
        match self.family.as_deref() {
            Some(family) if !family.eq_ignore_ascii_case(self.class_label.trim()) => {
                format!("{} / {}", self.class_label, family)
            }
            _ => self.class_label.clone(),
        }
    }

    pub fn display_label(&self) -> String {
        let class_family = self.class_family_label();
        if let Some(name) = self.name.as_deref().filter(|name| !name.trim().is_empty()) {
            if class_family.eq_ignore_ascii_case("repeat") {
                return name.trim().to_string();
            }
            if name
                .to_ascii_lowercase()
                .contains(&class_family.to_ascii_lowercase())
            {
                return name.trim().to_string();
            }
            return format!("{} ({})", name.trim(), class_family);
        }
        if class_family.eq_ignore_ascii_case("repeat") {
            "Repeat".to_string()
        } else {
            class_family
        }
    }

    pub fn group_label(&self) -> String {
        self.class_family_label()
    }
}

const REPEAT_KIND_TOKENS: &[&str] = &[
    "REPEAT_REGION",
    "REPEAT",
    "RMSK",
    "MOBILE_ELEMENT",
    "TRANSPOSON",
];

const REPEAT_NAME_KEYS: &[&str] = &[
    "rmsk_name",
    "repName",
    "rep_name",
    "repeat_name",
    "repeatName",
    "rpt_name",
    "rptName",
    "repeat_id",
    "repeatId",
];

const REPEAT_CLASS_KEYS: &[&str] = &[
    "rmsk_class",
    "repClass",
    "rep_class",
    "repeat_class",
    "repeatClass",
    "rpt_type",
    "rpt_class",
    "repeat_type",
    "mobile_element_type",
];

const REPEAT_FAMILY_KEYS: &[&str] = &[
    "rmsk_family",
    "repFamily",
    "rep_family",
    "repeat_family",
    "repeatFamily",
    "rpt_family",
    "rptFamily",
    "repeat_subfamily",
    "subfamily",
];

const GENERIC_REPEAT_LABEL_KEYS: &[&str] =
    &["standard_name", "name", "label", "mobile_element", "note"];

pub fn repeat_feature_display(feature: &Feature) -> Option<RepeatFeatureDisplay> {
    let kind = feature.kind.to_string().trim().to_ascii_uppercase();
    let has_repeat_kind = REPEAT_KIND_TOKENS
        .iter()
        .any(|token| kind == *token || kind.contains(token));

    let specific_name = first_qualifier_text(feature, REPEAT_NAME_KEYS);
    let class_raw = first_qualifier_text(feature, REPEAT_CLASS_KEYS);
    let family_raw = first_qualifier_text(feature, REPEAT_FAMILY_KEYS);
    if !has_repeat_kind && specific_name.is_none() && class_raw.is_none() && family_raw.is_none() {
        return None;
    }

    let name = specific_name.or_else(|| {
        if has_repeat_kind {
            first_qualifier_text(feature, GENERIC_REPEAT_LABEL_KEYS)
        } else {
            None
        }
    });
    let (class_text, class_family) = class_raw
        .as_deref()
        .and_then(split_class_family)
        .unwrap_or_else(|| (class_raw.clone().unwrap_or_default(), None));
    let family = family_raw
        .or(class_family)
        .map(|value| clean_repeat_component(&value))
        .filter(|value| !value.is_empty());

    let classifier_hint = if class_text.trim().is_empty() {
        name.as_deref().unwrap_or_default()
    } else {
        class_text.as_str()
    };
    let class = classify_repeat_class(classifier_hint);
    let class_label = if class == RepeatDisplayClass::Other && !class_text.trim().is_empty() {
        clean_repeat_component(&class_text)
    } else {
        class.label().to_string()
    };

    Some(RepeatFeatureDisplay {
        name: name
            .map(|value| clean_repeat_component(&value))
            .filter(|value| !value.is_empty()),
        class,
        class_label,
        family,
    })
}

pub fn is_repeat_feature(feature: &Feature) -> bool {
    repeat_feature_display(feature).is_some()
}

fn first_qualifier_text(feature: &Feature, keys: &[&str]) -> Option<String> {
    for key in keys {
        for (qualifier_key, value) in &feature.qualifiers {
            if !qualifier_key.to_string().eq_ignore_ascii_case(key) {
                continue;
            }
            let Some(value) = value else {
                continue;
            };
            let cleaned = clean_repeat_component(value);
            if !cleaned.is_empty() {
                return Some(cleaned);
            }
        }
    }
    None
}

fn split_class_family(raw: &str) -> Option<(String, Option<String>)> {
    let cleaned = clean_repeat_component(raw);
    if cleaned.is_empty() {
        return None;
    }
    for separator in ['/', '|', ':'] {
        if let Some((class, family)) = cleaned.split_once(separator) {
            let class = clean_repeat_component(class);
            let family = clean_repeat_component(family);
            return Some((class, (!family.is_empty()).then_some(family)));
        }
    }
    Some((cleaned, None))
}

fn clean_repeat_component(raw: &str) -> String {
    raw.split_whitespace()
        .collect::<Vec<_>>()
        .join(" ")
        .trim_matches(|ch: char| matches!(ch, '"' | '\'' | ';' | ',' | '(' | ')' | '[' | ']'))
        .trim()
        .to_string()
}

fn normalize_repeat_class_token(raw: &str) -> String {
    raw.chars()
        .map(|ch| {
            if ch.is_ascii_alphanumeric() {
                ch.to_ascii_lowercase()
            } else {
                ' '
            }
        })
        .collect::<String>()
        .split_whitespace()
        .collect::<Vec<_>>()
        .join(" ")
}

fn classify_repeat_class(raw: &str) -> RepeatDisplayClass {
    let normalized = normalize_repeat_class_token(raw);
    if normalized.is_empty() {
        return RepeatDisplayClass::Other;
    }
    let tokens = normalized.split_whitespace().collect::<Vec<_>>();
    if tokens.iter().any(|token| matches!(*token, "sine" | "alu")) {
        return RepeatDisplayClass::Sine;
    }
    if tokens.iter().any(|token| matches!(*token, "line" | "l1")) {
        return RepeatDisplayClass::Line;
    }
    if tokens.contains(&"ltr") || normalized.contains("erv") {
        return RepeatDisplayClass::Ltr;
    }
    if tokens.iter().any(|token| {
        matches!(
            *token,
            "dna" | "transposon" | "hat" | "mariner" | "tc1" | "piggybac"
        )
    }) {
        return RepeatDisplayClass::DnaTransposon;
    }
    if normalized.contains("simple repeat")
        || normalized.contains("tandem repeat")
        || normalized.contains("microsatellite")
    {
        return RepeatDisplayClass::SimpleRepeat;
    }
    if normalized.contains("low complexity") {
        return RepeatDisplayClass::LowComplexity;
    }
    if normalized.contains("satellite") || tokens.contains(&"sat") {
        return RepeatDisplayClass::Satellite;
    }
    if tokens.iter().any(|token| {
        matches!(
            *token,
            "rna" | "rrna" | "trna" | "srna" | "scrna" | "snrna" | "srprna"
        )
    }) {
        return RepeatDisplayClass::Rna;
    }
    if normalized.contains("retroposon") {
        return RepeatDisplayClass::Retroposon;
    }
    if normalized.contains("rolling circle")
        || normalized.contains("helitron")
        || tokens.contains(&"rc")
    {
        return RepeatDisplayClass::RollingCircle;
    }
    if normalized.contains("unknown") || tokens.contains(&"unk") {
        return RepeatDisplayClass::Unknown;
    }
    RepeatDisplayClass::Other
}

#[cfg(test)]
mod tests {
    use super::*;
    use gb_io::seq::Location;

    fn make_feature(kind: &str, qualifiers: &[(&str, &str)]) -> Feature {
        Feature {
            kind: kind.to_string().into(),
            location: Location::simple_range(0, 10),
            qualifiers: qualifiers
                .iter()
                .map(|(key, value)| ((*key).to_string().into(), Some((*value).to_string())))
                .collect(),
        }
    }

    #[test]
    fn classifies_ucsc_rmsk_sine_family() {
        let feature = make_feature(
            "repeat_region",
            &[
                ("repName", "AluY"),
                ("repClass", "SINE"),
                ("repFamily", "Alu"),
            ],
        );
        let display = repeat_feature_display(&feature).expect("repeat display");

        assert_eq!(display.class, RepeatDisplayClass::Sine);
        assert_eq!(display.display_label(), "AluY (SINE / Alu)");
        assert_eq!(display.group_label(), "SINE / Alu");
    }

    #[test]
    fn splits_combined_repeat_class_family_values() {
        let feature = make_feature(
            "mobile_element",
            &[
                ("repeat_name", "Charlie1"),
                ("repeat_class", "DNA/hAT-Charlie"),
            ],
        );
        let display = repeat_feature_display(&feature).expect("repeat display");

        assert_eq!(display.class, RepeatDisplayClass::DnaTransposon);
        assert_eq!(display.group_label(), "DNA transposon / hAT-Charlie");
    }

    #[test]
    fn detects_genbank_style_low_complexity_repeat_type() {
        let feature = make_feature(
            "repeat_region",
            &[("rpt_type", "Low_complexity"), ("rpt_family", "AT-rich")],
        );
        let display = repeat_feature_display(&feature).expect("repeat display");

        assert_eq!(display.class, RepeatDisplayClass::LowComplexity);
        assert_eq!(display.group_label(), "Low complexity / AT-rich");
    }

    #[test]
    fn ignores_generic_non_repeat_labels() {
        let feature = make_feature("misc_feature", &[("label", "AluY")]);

        assert!(repeat_feature_display(&feature).is_none());
    }
}
