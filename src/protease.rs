use serde::{Deserialize, Serialize};

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Protease {
    pub name: String,
    pub sequence: String,
    pub note: Option<String>,
    pub cut: isize,
}
