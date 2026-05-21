//! Embedded GUI translation catalogs and runtime language selection.
//!
//! This is deliberately presentation-only. Shared shell commands, protocol
//! schema fields, persisted scientific records, and adapter contracts stay in
//! deterministic English even when the GUI chrome is translated.

use std::{collections::BTreeMap, sync::OnceLock};

use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub(super) enum UiLanguage {
    #[serde(rename = "system")]
    System,
    #[serde(rename = "en-GB")]
    EnGb,
    #[serde(rename = "en-US")]
    EnUs,
    #[serde(rename = "de-DE")]
    DeDe,
    #[serde(rename = "fr-FR")]
    FrFr,
    #[serde(rename = "it-IT")]
    ItIt,
    #[serde(rename = "zh-Hans")]
    ZhHans,
    #[serde(rename = "ja-JP")]
    JaJp,
    #[serde(rename = "la")]
    La,
}

impl Default for UiLanguage {
    fn default() -> Self {
        Self::System
    }
}

impl UiLanguage {
    pub(super) const ALL: [Self; 9] = [
        Self::System,
        Self::EnGb,
        Self::EnUs,
        Self::DeDe,
        Self::FrFr,
        Self::ItIt,
        Self::ZhHans,
        Self::JaJp,
        Self::La,
    ];

    pub(super) fn id(self) -> &'static str {
        match self {
            Self::System => "system",
            Self::EnGb => "en-GB",
            Self::EnUs => "en-US",
            Self::DeDe => "de-DE",
            Self::FrFr => "fr-FR",
            Self::ItIt => "it-IT",
            Self::ZhHans => "zh-Hans",
            Self::JaJp => "ja-JP",
            Self::La => "la",
        }
    }

    pub(super) fn label(self) -> &'static str {
        match self {
            Self::System => "System default (English)",
            Self::EnGb => "English (UK)",
            Self::EnUs => "English (US)",
            Self::DeDe => "Deutsch",
            Self::FrFr => "Francais",
            Self::ItIt => "Italiano",
            Self::ZhHans => "Chinese, simplified (experimental)",
            Self::JaJp => "Japanese (experimental)",
            Self::La => "Latina (playful)",
        }
    }

    pub(super) fn effective(self) -> Self {
        match self {
            Self::System => Self::EnGb,
            other => other,
        }
    }

    fn catalog_json(self) -> &'static str {
        match self.effective() {
            Self::System | Self::EnGb => include_str!("../../assets/i18n/en-GB.json"),
            Self::EnUs => include_str!("../../assets/i18n/en-US.json"),
            Self::DeDe => include_str!("../../assets/i18n/de-DE.json"),
            Self::FrFr => include_str!("../../assets/i18n/fr-FR.json"),
            Self::ItIt => include_str!("../../assets/i18n/it-IT.json"),
            Self::ZhHans => include_str!("../../assets/i18n/zh-Hans.json"),
            Self::JaJp => include_str!("../../assets/i18n/ja-JP.json"),
            Self::La => include_str!("../../assets/i18n/la.json"),
        }
    }
}

#[derive(Debug, Clone)]
pub(super) struct I18n {
    language: UiLanguage,
}

impl Default for I18n {
    fn default() -> Self {
        Self {
            language: UiLanguage::default(),
        }
    }
}

impl I18n {
    pub(super) fn language(&self) -> UiLanguage {
        self.language
    }

    pub(super) fn set_language(&mut self, language: UiLanguage) {
        self.language = language;
    }

    pub(super) fn t(&self, key: &str) -> String {
        translate(self.language, key)
    }
}

fn parse_catalog(language: UiLanguage) -> BTreeMap<String, String> {
    serde_json::from_str(language.catalog_json()).unwrap_or_else(|err| {
        panic!(
            "embedded i18n catalog '{}' is invalid JSON: {err}",
            language.id()
        )
    })
}

fn catalog(language: UiLanguage) -> &'static BTreeMap<String, String> {
    static EN_GB: OnceLock<BTreeMap<String, String>> = OnceLock::new();
    static EN_US: OnceLock<BTreeMap<String, String>> = OnceLock::new();
    static DE_DE: OnceLock<BTreeMap<String, String>> = OnceLock::new();
    static FR_FR: OnceLock<BTreeMap<String, String>> = OnceLock::new();
    static IT_IT: OnceLock<BTreeMap<String, String>> = OnceLock::new();
    static ZH_HANS: OnceLock<BTreeMap<String, String>> = OnceLock::new();
    static JA_JP: OnceLock<BTreeMap<String, String>> = OnceLock::new();
    static LA: OnceLock<BTreeMap<String, String>> = OnceLock::new();

    match language.effective() {
        UiLanguage::System | UiLanguage::EnGb => {
            EN_GB.get_or_init(|| parse_catalog(UiLanguage::EnGb))
        }
        UiLanguage::EnUs => EN_US.get_or_init(|| parse_catalog(UiLanguage::EnUs)),
        UiLanguage::DeDe => DE_DE.get_or_init(|| parse_catalog(UiLanguage::DeDe)),
        UiLanguage::FrFr => FR_FR.get_or_init(|| parse_catalog(UiLanguage::FrFr)),
        UiLanguage::ItIt => IT_IT.get_or_init(|| parse_catalog(UiLanguage::ItIt)),
        UiLanguage::ZhHans => ZH_HANS.get_or_init(|| parse_catalog(UiLanguage::ZhHans)),
        UiLanguage::JaJp => JA_JP.get_or_init(|| parse_catalog(UiLanguage::JaJp)),
        UiLanguage::La => LA.get_or_init(|| parse_catalog(UiLanguage::La)),
    }
}

fn translate(language: UiLanguage, key: &str) -> String {
    if let Some(value) = catalog(language).get(key) {
        return value.clone();
    }
    catalog(UiLanguage::EnGb)
        .get(key)
        .cloned()
        .unwrap_or_else(|| key.to_string())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::BTreeSet;

    fn placeholders(value: &str) -> BTreeSet<String> {
        let mut found = BTreeSet::new();
        let mut cursor = value;
        while let Some(start) = cursor.find('{') {
            let after_start = &cursor[start + 1..];
            let Some(end) = after_start.find('}') else {
                break;
            };
            found.insert(after_start[..end].to_string());
            cursor = &after_start[end + 1..];
        }
        found
    }

    #[test]
    fn all_catalogs_parse_and_share_keys() {
        let base_keys = catalog(UiLanguage::EnGb)
            .keys()
            .cloned()
            .collect::<BTreeSet<_>>();
        assert!(!base_keys.is_empty());

        for language in UiLanguage::ALL {
            let keys = catalog(language).keys().cloned().collect::<BTreeSet<_>>();
            assert_eq!(keys, base_keys, "catalog key drift for {}", language.id());
        }
    }

    #[test]
    fn translated_placeholders_match_english_source() {
        let source = catalog(UiLanguage::EnGb);
        for language in UiLanguage::ALL {
            for (key, english_value) in source {
                assert_eq!(
                    placeholders(catalog(language).get(key).expect("catalog key exists")),
                    placeholders(english_value),
                    "placeholder drift for {} in {}",
                    key,
                    language.id()
                );
            }
        }
    }

    #[test]
    fn unknown_keys_fall_back_to_key_name() {
        let i18n = I18n {
            language: UiLanguage::DeDe,
        };

        assert_eq!(i18n.t("menu.file"), "Datei");
        assert_eq!(i18n.t("missing.example.key"), "missing.example.key");
    }

    #[test]
    fn language_ids_roundtrip_through_json() {
        for language in UiLanguage::ALL {
            let json = serde_json::to_string(&language).unwrap();
            assert_eq!(json, format!("\"{}\"", language.id()));
            let parsed: UiLanguage = serde_json::from_str(&json).unwrap();
            assert_eq!(parsed, language);
        }
    }
}
