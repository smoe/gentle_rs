use csv::ReaderBuilder;
use std::collections::HashMap;

pub struct Translations {
    values: HashMap<String, String>,
    language: String,
}

impl Translations {
    pub fn new() -> Self {
        let text = include_str!("../assets/translations.csv");

        let mut rdr = ReaderBuilder::new()
            .has_headers(true) // Specify if the CSV has headers
            .from_reader(text.as_bytes());

        let headers = rdr
            .headers()
            .expect("Could not read translations.csv headers");
        let mut languages = Self::to_vec(headers);
        let _ = languages.remove(0); // Remove the first column

        // Iterate over the records
        let mut values = HashMap::new();
        for record in rdr.records().flatten() {
            let mut record = Self::to_vec(&record);
            let key = record.remove(0);
            for (lnum, t) in record.iter().enumerate() {
                let lang_key = format!("{}:{key}", languages[lnum]);
                values.insert(lang_key, t.to_owned());
            }
        }

        Self {
            values,
            language: "en".to_owned(),
        }
    }

    pub fn set_language(&mut self, language: &str) {
        self.language = language.to_string();
    }

    pub fn get(&self, key: &str) -> Option<&str> {
        let key = format!("{}:{}", self.language, key);
        self.values.get(&key).map(|s| s.as_str())
    }

    fn to_vec(record: &csv::StringRecord) -> Vec<String> {
        record.iter().map(|s| s.to_string()).collect()
    }
}

impl Default for Translations {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default() {
        let translations = Translations::default();
        assert_eq!(translations.get("m_export_txt"), Some("Export to file"));
    }

    #[test]
    fn test_de() {
        let mut translations = Translations::default();
        translations.set_language("de");
        assert_eq!(
            translations.get("m_export_txt"),
            Some("In eine Datei exportieren")
        );
    }
}
