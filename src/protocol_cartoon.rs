//! Protocol-cartoon catalog and deterministic SVG rendering.
//!
//! This module provides a typed, engine-facing contract for protocol cartoon
//! generation. The initial scope is a conceptual two-fragment Gibson Assembly
//! strip rendered as deterministic SVG.

use serde::{Deserialize, Serialize};

/// Stable protocol-cartoon identifiers exposed through engine operations.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq, Hash)]
pub enum ProtocolCartoonKind {
    #[serde(rename = "gibson.two_fragment")]
    GibsonTwoFragment,
}

impl ProtocolCartoonKind {
    /// Canonical string identifier used in shell/CLI input and output rows.
    pub fn id(&self) -> &'static str {
        match self {
            Self::GibsonTwoFragment => "gibson.two_fragment",
        }
    }

    /// Human-readable title for UI/help surfaces.
    pub fn title(&self) -> &'static str {
        match self {
            Self::GibsonTwoFragment => "Gibson Assembly (two-fragment conceptual)",
        }
    }

    /// Short semantics summary used in list outputs.
    pub fn summary(&self) -> &'static str {
        match self {
            Self::GibsonTwoFragment => {
                "Context + chew-back -> anneal -> extend -> seal -> assembled molecule"
            }
        }
    }

    /// Parse supported shell/CLI aliases into the canonical enum.
    pub fn parse_id(raw: &str) -> Option<Self> {
        match raw.trim().to_ascii_lowercase().as_str() {
            "gibson.two_fragment" | "gibson.two-fragment" | "gibson_two_fragment" | "gibson" => {
                Some(Self::GibsonTwoFragment)
            }
            _ => None,
        }
    }

    /// Deterministic ordered catalog for list commands.
    pub fn catalog() -> Vec<Self> {
        vec![Self::GibsonTwoFragment]
    }
}

/// Machine-readable row used by list commands.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct ProtocolCartoonCatalogRow {
    pub id: String,
    pub title: String,
    pub summary: String,
}

pub fn protocol_cartoon_catalog_rows() -> Vec<ProtocolCartoonCatalogRow> {
    ProtocolCartoonKind::catalog()
        .into_iter()
        .map(|kind| ProtocolCartoonCatalogRow {
            id: kind.id().to_string(),
            title: kind.title().to_string(),
            summary: kind.summary().to_string(),
        })
        .collect()
}

/// Render one protocol cartoon as deterministic SVG.
pub fn render_protocol_cartoon_svg(kind: &ProtocolCartoonKind) -> String {
    match kind {
        ProtocolCartoonKind::GibsonTwoFragment => render_gibson_two_fragment_svg(),
    }
}

fn render_gibson_two_fragment_svg() -> String {
    r##"<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1680 760" role="img" aria-labelledby="title desc">
  <title id="title">GENtle protocol cartoon: Gibson Assembly two-fragment conceptual flow</title>
  <desc id="desc">Reaction context, canonical enzyme steps, and final assembled molecule.</desc>
  <defs>
    <linearGradient id="bg" x1="0" y1="0" x2="1" y2="1">
      <stop offset="0%" stop-color="#f5fbff"/>
      <stop offset="100%" stop-color="#edf8f5"/>
    </linearGradient>
    <linearGradient id="panel" x1="0" y1="0" x2="0" y2="1">
      <stop offset="0%" stop-color="#ffffff"/>
      <stop offset="100%" stop-color="#f6fafb"/>
    </linearGradient>
    <linearGradient id="build" x1="0" y1="0" x2="1" y2="0">
      <stop offset="0%" stop-color="#0f5d75"/>
      <stop offset="48%" stop-color="#0f5d75"/>
      <stop offset="52%" stop-color="#e3b230"/>
      <stop offset="100%" stop-color="#e3b230"/>
    </linearGradient>
    <marker id="arrow" markerWidth="12" markerHeight="12" refX="9" refY="6" orient="auto">
      <path d="M1,1 L11,6 L1,11 Z" fill="#4e6d78"/>
    </marker>
    <style>
      .title { font: 700 39px "Avenir Next", "Trebuchet MS", "Segoe UI", sans-serif; fill: #0f2a36; }
      .subtitle { font: 500 19px "Avenir Next", "Trebuchet MS", "Segoe UI", sans-serif; fill: #36535d; }
      .panel-title { font: 700 25px "Avenir Next", "Trebuchet MS", "Segoe UI", sans-serif; fill: #113242; }
      .panel-sub { font: 500 16px "Avenir Next", "Trebuchet MS", "Segoe UI", sans-serif; fill: #3d5b66; }
      .step-name { font: 700 22px "Avenir Next", "Trebuchet MS", "Segoe UI", sans-serif; fill: #113242; }
      .tiny { font: 500 14px "Avenir Next", "Trebuchet MS", "Segoe UI", sans-serif; fill: #44626d; }
      .cap { font: 600 15px "Avenir Next", "Trebuchet MS", "Segoe UI", sans-serif; fill: #2c4955; }
      .panel { fill: url(#panel); stroke: #c8dde3; stroke-width: 2; rx: 22; }
      .dna-top { fill: #0f5d75; }
      .dna-bot { fill: #e3b230; }
      .overlap { fill: #1aa3a1; }
      .hint { stroke: #88a4ae; stroke-width: 2.5; stroke-dasharray: 8 7; fill: none; }
      .flow { stroke: #4e6d78; stroke-width: 4; marker-end: url(#arrow); fill: none; }
      .nick { stroke: #8f4b2e; stroke-width: 4; stroke-linecap: round; }
    </style>
  </defs>

  <rect x="0" y="0" width="1680" height="760" fill="url(#bg)"/>
  <rect x="18" y="18" width="1644" height="724" rx="24" fill="none" stroke="#d4e4e8" stroke-width="2"/>

  <text x="44" y="68" class="title">GENtle Protocol Cartoon: Gibson Assembly</text>
  <text x="44" y="98" class="subtitle">Conceptual two-fragment flow with canonical step semantics.</text>

  <rect class="panel" x="40" y="132" width="260" height="570"/>
  <text x="64" y="175" class="panel-title">Reaction Context</text>
  <text x="64" y="205" class="panel-sub">Two linear dsDNA fragments</text>
  <text x="64" y="232" class="panel-sub">Designed overlap: 20-40 bp</text>
  <text x="64" y="259" class="panel-sub">Isothermal mix: 50 C</text>

  <rect x="66" y="300" width="165" height="14" rx="7" class="dna-top"/>
  <rect x="66" y="321" width="165" height="14" rx="7" class="dna-bot"/>
  <rect x="212" y="300" width="19" height="14" rx="7" class="overlap"/>
  <rect x="212" y="321" width="19" height="14" rx="7" class="overlap"/>
  <text x="66" y="351" class="tiny">Fragment A with overlap tail</text>

  <rect x="110" y="404" width="165" height="14" rx="7" class="dna-top"/>
  <rect x="110" y="425" width="165" height="14" rx="7" class="dna-bot"/>
  <rect x="110" y="404" width="19" height="14" rx="7" class="overlap"/>
  <rect x="110" y="425" width="19" height="14" rx="7" class="overlap"/>
  <text x="66" y="455" class="tiny">Fragment B with matching overlap</text>

  <rect class="panel" x="332" y="132" width="210" height="570"/>
  <circle cx="364" cy="177" r="18" fill="#2b8ea8"/>
  <text x="358" y="183" fill="#ffffff" font-family="Avenir Next, Trebuchet MS, Segoe UI, sans-serif" font-size="16" font-weight="700">1</text>
  <text x="393" y="185" class="step-name">Chew-back</text>
  <text x="364" y="214" class="tiny">Exonuclease exposes</text>
  <text x="364" y="234" class="tiny">single-strand overlaps</text>
  <text x="355" y="506" class="cap">Complementary tails can anneal</text>

  <rect class="panel" x="554" y="132" width="210" height="570"/>
  <circle cx="586" cy="177" r="18" fill="#2b8ea8"/>
  <text x="580" y="183" fill="#ffffff" font-family="Avenir Next, Trebuchet MS, Segoe UI, sans-serif" font-size="16" font-weight="700">2</text>
  <text x="615" y="185" class="step-name">Anneal</text>
  <text x="586" y="214" class="tiny">Complementary overhangs</text>
  <text x="586" y="234" class="tiny">base-pair in register</text>
  <path d="M648 300 C657 286, 661 286, 670 300" class="hint"/>
  <path d="M648 390 C657 404, 661 404, 670 390" class="hint"/>
  <text x="577" y="506" class="cap">Transient duplex seeds extension</text>

  <rect class="panel" x="776" y="132" width="210" height="570"/>
  <circle cx="808" cy="177" r="18" fill="#2b8ea8"/>
  <text x="802" y="183" fill="#ffffff" font-family="Avenir Next, Trebuchet MS, Segoe UI, sans-serif" font-size="16" font-weight="700">3</text>
  <text x="837" y="185" class="step-name">Extend</text>
  <text x="808" y="214" class="tiny">Polymerase fills across</text>
  <text x="808" y="234" class="tiny">annealed overlap gaps</text>
  <text x="799" y="506" class="cap">Double strands are restored</text>

  <rect class="panel" x="998" y="132" width="210" height="570"/>
  <circle cx="1030" cy="177" r="18" fill="#2b8ea8"/>
  <text x="1024" y="183" fill="#ffffff" font-family="Avenir Next, Trebuchet MS, Segoe UI, sans-serif" font-size="16" font-weight="700">4</text>
  <text x="1059" y="185" class="step-name">Seal</text>
  <text x="1030" y="214" class="tiny">Ligase closes residual</text>
  <text x="1030" y="234" class="tiny">backbone nicks</text>
  <line x1="1085" y1="329" x2="1091" y2="343" class="nick"/>
  <line x1="1117" y1="349" x2="1123" y2="363" class="nick"/>
  <text x="1021" y="506" class="cap">Covalent continuity complete</text>

  <rect class="panel" x="1220" y="132" width="420" height="570"/>
  <text x="1248" y="185" class="panel-title">Assembled Molecule</text>
  <text x="1248" y="215" class="panel-sub">Conceptual two-fragment product</text>

  <rect x="1270" y="330" width="320" height="15" rx="7" fill="url(#build)"/>
  <rect x="1270" y="352" width="320" height="15" rx="7" fill="url(#build)"/>
  <rect x="1421" y="330" width="18" height="15" rx="7" fill="#19a19f"/>
  <rect x="1421" y="352" width="18" height="15" rx="7" fill="#19a19f"/>
  <text x="1270" y="397" class="tiny">Segments are joined through designed overlap.</text>

  <path d="M305 417 H325" class="flow"/>
  <path d="M547 417 H565" class="flow"/>
  <path d="M769 417 H787" class="flow"/>
  <path d="M991 417 H1009" class="flow"/>
  <path d="M1213 417 H1231" class="flow"/>
</svg>
"##
    .to_string()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_protocol_cartoon_aliases() {
        assert_eq!(
            ProtocolCartoonKind::parse_id("gibson.two_fragment"),
            Some(ProtocolCartoonKind::GibsonTwoFragment)
        );
        assert_eq!(
            ProtocolCartoonKind::parse_id("gibson"),
            Some(ProtocolCartoonKind::GibsonTwoFragment)
        );
        assert!(ProtocolCartoonKind::parse_id("unknown.protocol").is_none());
    }

    #[test]
    fn catalog_rows_are_deterministic() {
        let rows = protocol_cartoon_catalog_rows();
        assert_eq!(rows.len(), 1);
        assert_eq!(rows[0].id, "gibson.two_fragment");
        assert!(rows[0].title.contains("Gibson"));
    }

    #[test]
    fn gibson_svg_contains_canonical_steps() {
        let svg = render_protocol_cartoon_svg(&ProtocolCartoonKind::GibsonTwoFragment);
        assert!(svg.contains("<svg"));
        assert!(svg.contains("Chew-back"));
        assert!(svg.contains("Anneal"));
        assert!(svg.contains("Extend"));
        assert!(svg.contains("Seal"));
        assert!(svg.contains("Assembled Molecule"));
    }
}
