//! Embedded GUI translation catalogs and runtime language selection.
//!
//! This is deliberately presentation-only. Shared shell commands, protocol
//! schema fields, persisted scientific records, and adapter contracts stay in
//! deterministic English even when the GUI chrome is translated.

use std::{
    collections::BTreeMap,
    sync::{OnceLock, RwLock},
};

use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default)]
pub(crate) enum UiLanguage {
    #[serde(rename = "system")]
    #[default]
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

impl UiLanguage {
    pub(crate) const ALL: [Self; 9] = [
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

    pub(crate) fn id(self) -> &'static str {
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

    pub(crate) fn label(self) -> &'static str {
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

    pub(crate) fn effective(self) -> Self {
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

#[derive(Debug, Clone, Default)]
pub(crate) struct I18n {
    language: UiLanguage,
}

impl I18n {
    #[cfg(test)]
    pub(crate) fn for_test_language(language: UiLanguage) -> Self {
        Self { language }
    }

    pub(crate) fn language(&self) -> UiLanguage {
        self.language
    }

    pub(crate) fn set_language(&mut self, language: UiLanguage) {
        self.language = language;
        set_current_language(language);
    }

    pub(crate) fn t(&self, key: &str) -> String {
        translate(self.language, key)
    }

    pub(crate) fn tf(&self, key: &str, values: &[(&str, &str)]) -> String {
        format_translation(&self.t(key), values)
    }

    /// Translate bundled catalog prose without replacing user-authored catalog overrides.
    pub(crate) fn catalog_text(&self, key: &str, source: &str) -> String {
        if catalog(UiLanguage::EnGb)
            .get(key)
            .is_some_and(|text| text == source)
        {
            self.t(key)
        } else {
            source.to_string()
        }
    }

    /// Match only complete, known GUI guidance; never translate a provider's free text.
    pub(crate) fn agent_hint(&self, source: &str) -> String {
        catalog(UiLanguage::EnGb)
            .iter()
            .find(|(key, value)| key.starts_with("agent.hint.") && value.as_str() == source)
            .map(|(key, _)| self.t(key))
            .unwrap_or_else(|| source.to_string())
    }
}

// Substitute only catalog placeholders, never braces contained in user/provider data.
fn format_translation(template: &str, values: &[(&str, &str)]) -> String {
    let mut result = String::with_capacity(template.len());
    let mut remaining = template;
    while let Some(start) = remaining.find('{') {
        result.push_str(&remaining[..start]);
        let Some(end) = remaining[start..].find('}') else {
            result.push_str(&remaining[start..]);
            return result;
        };
        let name = &remaining[start + 1..start + end];
        if let Some((_, value)) = values.iter().find(|(key, _)| *key == name) {
            result.push_str(value);
        } else {
            result.push_str(&remaining[start..=start + end]);
        }
        remaining = &remaining[start + end + 1..];
    }
    result.push_str(remaining);
    result
}

fn current_language_cell() -> &'static RwLock<UiLanguage> {
    static CURRENT_LANGUAGE: OnceLock<RwLock<UiLanguage>> = OnceLock::new();
    CURRENT_LANGUAGE.get_or_init(|| RwLock::new(UiLanguage::default()))
}

pub(crate) fn set_current_language(language: UiLanguage) {
    if let Ok(mut guard) = current_language_cell().write() {
        *guard = language;
    }
}

pub(crate) fn current_language() -> UiLanguage {
    #[cfg(test)]
    if let Some(language) = TEST_LANGUAGE.with(std::cell::Cell::get) {
        return language;
    }
    current_language_cell()
        .read()
        .map(|guard| *guard)
        .unwrap_or_default()
}

pub(crate) fn tr(key: &str) -> String {
    translate(current_language(), key)
}

pub(crate) fn trf(key: &str, values: &[(&str, &str)]) -> String {
    format_translation(&tr(key), values)
}

#[cfg(test)]
thread_local! {
    static TEST_LANGUAGE: std::cell::Cell<Option<UiLanguage>> = const { std::cell::Cell::new(None) };
}

/// Keep rendered labels stable while other GUI tests change the application language.
#[cfg(test)]
pub(crate) struct TestLanguageGuard {
    previous: Option<UiLanguage>,
    _thread_bound: std::marker::PhantomData<std::rc::Rc<()>>,
}

#[cfg(test)]
impl TestLanguageGuard {
    pub(crate) fn new(language: UiLanguage) -> Self {
        Self {
            previous: TEST_LANGUAGE.with(|current| current.replace(Some(language))),
            _thread_bound: std::marker::PhantomData,
        }
    }
}

#[cfg(test)]
impl Drop for TestLanguageGuard {
    fn drop(&mut self) {
        TEST_LANGUAGE.with(|current| current.set(self.previous));
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

    #[test]
    fn test_language_override_is_nested_and_thread_local() {
        let _language = TestLanguageGuard::new(UiLanguage::DeDe);
        assert_eq!(tr("menu.file"), "Datei");
        std::thread::spawn(|| {
            let _language = TestLanguageGuard::new(UiLanguage::EnGb);
            assert_eq!(tr("menu.file"), "File");
        })
        .join()
        .unwrap();
        {
            let _nested = TestLanguageGuard::new(UiLanguage::EnGb);
            assert_eq!(tr("menu.file"), "File");
        }
        assert_eq!(tr("menu.file"), "Datei");
    }

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
    fn agent_i18n_formatting_preserves_inserted_provider_text() {
        let i18n = I18n {
            language: UiLanguage::DeDe,
        };
        let raw = "model {code} / unchanged";
        assert_eq!(
            i18n.tf("agent.display.latest", &[("label", raw), ("id", "custom")]),
            "Letzte Antwort von model {code} / unchanged (custom)"
        );
        assert_eq!(
            format_translation("{a} {b}", &[("a", "{b}"), ("b", "untouched")]),
            "{b} untouched"
        );
        assert_eq!(format_translation("{missing} {", &[]), "{missing} {");
    }

    #[test]
    fn agent_i18n_covers_bundled_providers_templates_and_ui_keys() {
        let english = catalog(UiLanguage::EnGb);
        let providers: serde_json::Value =
            serde_json::from_str(include_str!("../../assets/agent_systems.json")).unwrap();
        for provider in providers["systems"].as_array().unwrap() {
            for field in ["label", "description"] {
                let key = format!(
                    "agent.provider.{}.{field}",
                    provider["id"].as_str().unwrap()
                );
                assert_eq!(
                    english.get(&key).map(String::as_str),
                    provider[field].as_str(),
                    "{key}"
                );
            }
        }
        let key_pattern = regex::Regex::new(r#"(?:\.trf?|::trf?)\(\s*"(agent\.[^"]+)""#).unwrap();
        let mut checked = BTreeSet::new();
        for source in [
            include_str!("routine_and_agent_assistant_ui.rs"),
            include_str!("../app.rs"),
        ] {
            for captures in key_pattern.captures_iter(source) {
                let key = captures[1].to_string();
                assert!(english.contains_key(&key), "missing agent GUI key: {key}");
                checked.insert(key);
            }
        }
        assert!(
            checked.len() > 180,
            "agent translation guard must cover the full surface"
        );
        for prefix in [
            "agent.ui.",
            "agent.display.",
            "agent.status.",
            "agent.provider.",
            "agent.template.",
            "agent.hint.",
        ] {
            for (key, _) in english.iter().filter(|(key, _)| key.starts_with(prefix)) {
                for language in UiLanguage::ALL {
                    assert!(
                        !catalog(language)[key].trim().is_empty(),
                        "{}: {key}",
                        language.id()
                    );
                }
            }
        }
    }

    #[test]
    fn agent_i18n_preserves_custom_catalog_text_and_unknown_diagnostics() {
        let i18n = I18n {
            language: UiLanguage::DeDe,
        };
        assert_eq!(
            i18n.catalog_text("agent.provider.builtin_echo.label", "Built-in Echo (demo)"),
            "Integriertes Echo (Demo)"
        );
        assert_eq!(
            i18n.catalog_text("agent.provider.builtin_echo.label", "My laboratory agent"),
            "My laboratory agent"
        );
        assert_eq!(
            i18n.catalog_text("agent.provider.custom.label", "Private model"),
            "Private model"
        );
        assert_eq!(
            i18n.agent_hint("Provider error {code}: raw diagnostics"),
            "Provider error {code}: raw diagnostics"
        );
        assert_ne!(
            i18n.agent_hint(crate::agent_bridge::ANTHROPIC_API_KEY_AUTH_HINT),
            crate::agent_bridge::ANTHROPIC_API_KEY_AUTH_HINT
        );
    }

    #[test]
    fn german_about_label_is_localized() {
        let i18n = I18n {
            language: UiLanguage::DeDe,
        };

        assert_eq!(i18n.t("menu.help.about"), "Über GENtle");
        assert_eq!(i18n.t("about.more_info"), "Weitere Infos …");
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
