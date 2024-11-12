use serde_json::Value;
use std::collections::HashMap;

#[derive(Clone, Debug, Default)]
pub struct DNAmarkerPart {
    pub length: f64,
    pub strength: Option<u64>,
}

#[derive(Clone, Debug, Default)]
pub struct DNAmarker {
    name: String,
    parts: Vec<DNAmarkerPart>,
}

impl DNAmarker {
    pub fn new(name: &str, parts: &Value) -> Self {
        let parts: Vec<DNAmarkerPart> = parts
            .as_array()
            .expect("DNA marker part is not an array")
            .iter()
            .map(|p| p.as_array().expect("DNA marker subpart not an array"))
            .map(|p| DNAmarkerPart {
                length: p[0].as_f64().expect("Not an f64"),
                strength: match p.get(1) {
                    Some(s) => s.as_u64(),
                    None => None,
                },
            })
            .collect();
        Self {
            name: name.to_owned(),
            parts,
        }
    }

    pub fn name(&self) -> &str {
        &self.name
    }

    pub fn parts(&self) -> &Vec<DNAmarkerPart> {
        &self.parts
    }
}

#[derive(Clone, Debug)]
pub struct DNAMarkers(HashMap<String, DNAmarker>);

impl DNAMarkers {
    pub fn get(&self, name: &str) -> Option<&DNAmarker> {
        self.0.get(name)
    }
}

impl Default for DNAMarkers {
    fn default() -> Self {
        let mut ret = HashMap::new();
        let data = include_str!("../assets/dna_markers.json");
        let res: Value = serde_json::from_str(data).expect("Invalid JSON");
        let map = res.as_object().expect("JSON is not an object");
        for (name, parts) in map.iter() {
            let dna_marker = DNAmarker::new(name, parts);
            ret.insert(name.to_owned(), dna_marker);
        }
        Self(ret)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_dna_marker() {
        let dna_markers = DNAMarkers::default();
        let marker = dna_markers.get("GeneRuler Mix").unwrap();
        assert_eq!(marker.parts()[2].length, 8000.0);
        assert_eq!(marker.parts()[2].strength, Some(13));
    }
}
